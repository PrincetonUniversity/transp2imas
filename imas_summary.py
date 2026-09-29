#!/usr/bin/python
"""Populate the IMAS summary IDS from global TRANSP time traces and simulation metadata."""

import numpy as np
from scipy import constants

from nc2imas import nc2imas


class imas_summary:
    """Write a curated TRANSP summary IDS from high-confidence CDF quantities."""

    def __init__(self, ncimas: nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.create_ids()

    def _to_summary_time(self, data):
        """Map a 1-D profile-time quantity onto summary/TIME."""
        arr = np.asarray(data)
        if arr.size == 0:
            return arr
        if arr.ndim == 0:
            return arr
        if arr.ndim != 1:
            return np.array([])
        if len(arr) == self.im.ntimes:
            return arr.astype(np.float64, copy=False)
        if len(arr) == self.im.ntimep and self.im.is_valid(self.im.timep, self.im.timesc):
            return np.interp(self.im.timesc, self.im.timep, arr)
        return arr

    def _integrate_profile(self, variable, multiplier=1.0):
        """Volume integral for a TRANSP X-grid profile using DVOL in cm^3."""
        val = self.im.get_ncvar(variable)
        dvol = self.im.get_ncvar('DVOL')
        if not self.im.is_valid(val, dvol) or val.ndim != 2 or val.shape != dvol.shape:
            return np.array([])
        out = multiplier * np.sum(val * dvol, axis=1)
        return self._to_summary_time(out)

    def _fill(self, parent, name, data, source=None, scalar=False):
        if not hasattr(parent, name) or not self.im.is_valid(data):
            return
        node = getattr(parent, name)
        if scalar:
            self.im.fill_value(node, data, rsz=0)
        else:
            self.im.fill_value(node, self._to_summary_time(data))
        if source and hasattr(node, 'source'):
            node.source = source

    def _fill_nc(self, parent, name, cdf_name, multiplier=1.0, source=None, integer=False):
        data = self.im.get_ncvar(cdf_name)
        if not self.im.is_valid(data):
            return
        data = multiplier * data
        if integer:
            data = np.rint(data).astype(np.int32)
        self._fill(parent, name, data, source or f'TRANSP CDF:{cdf_name}')

    def create_ids(self):
        self.logger.info('\nsummary')
        sid = self.im.new_ids('summary')
        sid = self.im.fill_properties(sid)
        # Store the same TRANSP/template provenance in summary/code.parameters.
        sid = self.im.fill_code_properties(sid)

        sid.ids_properties.homogeneous_time = 1
        times = np.asarray(self.im.timesc, dtype=np.float64)
        sid.time.resize(len(times))
        sid.time[:] = times

        # DD4.1+ dataset/simulation metadata. Do not fabricate experimental
        # wall-clock timestamps; TRANSP time is relative simulation time.
        if hasattr(sid, 'type'):
            sid.type.name = 'simulation'
            sid.type.index = 2
            if hasattr(sid.type, 'description'):
                sid.type.description = 'Data coming from a simulation (may still use or refer to an experiment)'
        if hasattr(sid, 'description'):
            sid.description = f'TRANSP integrated simulation for run {self.im.device} {self.im.shot}'
        if hasattr(sid, 'machine'):
            sid.machine = str(self.im.device)
        if hasattr(sid, 'pulse'):
            try:
                sid.pulse = int(self.im.shot)
            except Exception:
                pass
        if hasattr(sid, 'simulation'):
            sim = sid.simulation
            if times.size and hasattr(sim, 'time_begin'):
                sim.time_begin = float(times[0])
            if times.size and hasattr(sim, 'time_end'):
                sim.time_end = float(times[-1])
            if times.size > 1 and hasattr(sim, 'time_step'):
                dt = np.diff(times)
                if np.all(np.isfinite(dt)) and np.allclose(dt, dt[0], rtol=1e-8, atol=1e-12):
                    sim.time_step = float(dt[0])
            if hasattr(sim, 'workflow'):
                sim.workflow = 'TRANSP integrated modeling; detailed provenance in workflow IDS and code.parameters'

        gq = sid.global_quantities
        ip_sign = self.im.ip_sign
        bt_sign = self.im.bt_sign
        q_sign = self.im.q_sign

        # Basic global quantities.
        self._fill_nc(gq, 'ip', 'PCURC', ip_sign)
        self._fill_nc(gq, 'b0', 'BZ', bt_sign)
        bz = self.im.get_ncvar('BZ')
        bzxr = self.im.get_ncvar('BZXR', self.im.cm2m)  # T*cm -> T*m
        if self.im.is_valid(bz, bzxr):
            with np.errstate(divide='ignore', invalid='ignore'):
                r0_t = np.where(np.abs(bz) > 0.0, bzxr/bz, np.nan)
            finite = r0_t[np.isfinite(r0_t)]
            if finite.size:
                self._fill(gq, 'r0', float(np.median(finite)), 'TRANSP CDF:BZXR/BZ', scalar=True)

        # li was renamed li_3 in DD4.
        li_name = 'li_3' if hasattr(gq, 'li_3') else 'li'
        self._fill_nc(gq, li_name, 'LI_3')
        self._fill_nc(gq, 'v_loop', 'VSURC', self.im.vloop_sign)
        self._fill_nc(gq, 'beta_pol', 'BETAT')
        self._fill_nc(gq, 'beta_tor', 'BTEQ')
        self._fill_nc(gq, 'volume', 'PVOL', 1.0e-6)
        self._fill_nc(gq, 'q_95', 'Q95', q_sign)
        self._fill_nc(gq, 'tau_energy', 'TAUA1')
        self._fill_nc(gq, 'h_98', 'H98Y2')
        self._fill_nc(gq, 'tau_energy_98', 'TAUE98Y2')
        self._fill_nc(gq, 'h_mode', 'LHMODE', integer=True)

        # Normalized beta: 100 beta_tor a[m] B0[T] / Ip[MA].  The physical
        # signs of B0 and Ip are retained, matching the IMAS convention.
        beta = self.im.get_ncvar('BTEQ')
        area = self.im.get_ncvar('PAREA', 1.0e-4)
        ip = ip_sign * self.im.get_ncvar('PCURC')
        b0 = bt_sign * bz
        beta_n = np.array([])
        if self.im.is_valid(beta, area, ip, b0):
            try:
                minor = np.sqrt(np.clip(area/self.im.pi, 0.0, None))
                with np.errstate(divide='ignore', invalid='ignore'):
                    beta_n = 1.0e8 * beta * minor * b0 / ip
            except Exception:
                beta_n = np.array([])
        beta_name = 'beta_tor_norm' if hasattr(gq, 'beta_tor_norm') else 'beta_normal'
        self._fill(gq, beta_name, beta_n, 'Derived from BTEQ,PAREA,BZ,PCURC')

        # Energy contents. UE/UI/UTHRM/UTOTL are already energy densities in
        # J/cm^3, so their volume integrals are energies in J -- no 3/2 factor.
        self._fill(gq, 'energy_total', self._integrate_profile('UTOTL'), 'Integral TRANSP CDF:UTOTL*DVOL')
        self._fill(gq, 'energy_thermal', self._integrate_profile('UTHRM'), 'Integral TRANSP CDF:UTHRM*DVOL')
        self._fill(gq, 'energy_electrons_thermal', self._integrate_profile('UE'), 'Integral TRANSP CDF:UE*DVOL')
        self._fill(gq, 'energy_ion_total_thermal', self._integrate_profile('UI'), 'Integral TRANSP CDF:UI*DVOL')

        # MHD pressure and diamagnetic energy use pressure definitions from the DD.
        ptowb = self.im.get_ncvar('PTOWB')
        dvol = self.im.get_ncvar('DVOL', 1.0e-6)  # m^3
        if self.im.is_valid(ptowb, dvol) and ptowb.shape == dvol.shape:
            self._fill(gq, 'energy_mhd', 1.5*np.sum(ptowb*dvol, axis=1),
                       '1.5 integral TRANSP CDF:PTOWB dV')

        pthermal = self.im.get_ncvar('PPLAS')
        ufast_perp = self.im.get_ncvar('UFASTPP', self.im.vol_inv)  # J/m3 = Pa_perp
        ufast_parallel = self.im.get_ncvar('UFASTPA', self.im.vol_inv)
        if self.im.is_valid(pthermal, ufast_perp, dvol):
            pperp = pthermal + ufast_perp
            self._fill(gq, 'energy_diamagnetic', 1.5*np.sum(pperp*dvol, axis=1),
                       '1.5 integral (PPLAS + UFASTPP) dV')
        if self.im.is_valid(ufast_perp):
            self._fill(gq, 'energy_fast_perpendicular',
                       1.5*self._integrate_profile('UFASTPP'),
                       '1.5 integral TRANSP CDF:UFASTPP*DVOL')
        if self.im.is_valid(ufast_parallel):
            # UFASTPA = p_parallel/2; DD quantity = 1.5 integral p_parallel dV.
            self._fill(gq, 'energy_fast_parallel',
                       3.0*self._integrate_profile('UFASTPA'),
                       '3 integral TRANSP CDF:UFASTPA*DVOL')

        # Poloidal magnetic-field energy.
        bpol = self.im.get_ncvar('BPOL')
        coord_x = self.im.get_ncvar('X')
        coord_xb = self.im.get_ncvar('XB')
        if self.im.is_valid(bpol, dvol, coord_x, coord_xb):
            try:
                bpol_x = self.im.interp_radial_profiles(bpol, coord_xb, coord_x, axis_value=0.0)
                ebpol = np.sum((bpol_x*bpol_x)/(2.0*constants.mu_0) * dvol, axis=1)
                self._fill(gq, 'energy_b_field_pol', ebpol, 'Integral TRANSP CDF:BPOL^2/(2 mu0) dV')
            except Exception:
                pass

        # Toroidal current components; TRANSP CUR* [A/cm^2] times DAREA [cm^2]
        # gives A directly.
        darea = self.im.get_ncvar('DAREA')
        cur = self.im.get_ncvar('CUR')
        curoh = self.im.get_ncvar('CUROH')
        curbs = self.im.get_ncvar('CURBS')
        if self.im.is_valid(darea, curoh) and darea.shape == curoh.shape:
            self._fill(gq, 'current_ohm', ip_sign*np.sum(curoh*darea, axis=1),
                       'Integral TRANSP CDF:CUROH*DAREA')
        if self.im.is_valid(darea, cur, curoh) and cur.shape == darea.shape == curoh.shape:
            self._fill(gq, 'current_non_inductive', ip_sign*np.sum((cur-curoh)*darea, axis=1),
                       'Integral TRANSP CDF:(CUR-CUROH)*DAREA')
        if self.im.is_valid(darea, curbs) and darea.shape == curbs.shape:
            self._fill(gq, 'current_bootstrap', ip_sign*np.sum(curbs*darea, axis=1),
                       'Integral TRANSP CDF:CURBS*DAREA')

        # Heating/loss powers.
        poh_total = self.im.get_ncvar('POHT')
        if self.im.is_valid(poh_total):
            self._fill(gq, 'power_ohm', poh_total, 'TRANSP CDF:POHT')
        else:
            self._fill(gq, 'power_ohm', self._integrate_profile('POH'), 'Integral TRANSP CDF:POH*DVOL')
        prad_total = self._integrate_profile('PRAD')
        self._fill(gq, 'power_radiated', prad_total, 'Integral TRANSP CDF:PRAD*DVOL')
        self._fill(gq, 'power_radiated_inside_lcfs', prad_total,
                   'Integral TRANSP CDF:PRAD*DVOL (X grid is inside LCFS)')
        self._fill(gq, 'power_bremsstrahlung', self._integrate_profile('PRAD_BR'),
                   'Integral TRANSP CDF:PRAD_BR*DVOL')
        self._fill(gq, 'power_line', self._integrate_profile('PRAD_LI'),
                   'Integral TRANSP CDF:PRAD_LI*DVOL')
        # IMAS calls this channel synchrotron; TRANSP names the same electron
        # cyclotron-radiation loss PRAD_CY.
        self._fill(gq, 'power_synchrotron', self._integrate_profile('PRAD_CY'),
                   'Integral TRANSP CDF:PRAD_CY*DVOL')

        # Heating/current-drive summary. PINJ is the launched NBI power and is
        # unambiguous.  BPCAP is deliberately not mapped to power_nbi because
        # IMAS defines coupled NBI power as also excluding fast-ion losses,
        # whereas TRANSP BPCAP is PINJ minus shine-through only.
        if hasattr(sid, 'heating_current_drive'):
            hcd = sid.heating_current_drive
            self._fill_nc(hcd, 'power_launched_nbi', 'PINJ')

        # Fusion neutron rates.  TRANSP provides both all-reaction totals and a
        # DD decomposition into thermal, beam-thermal, and beam-beam channels.
        if hasattr(sid, 'fusion') and hasattr(sid.fusion, 'neutron_rates'):
            nr = sid.fusion.neutron_rates
            self._fill_nc(nr, 'total', 'NEUTT')
            self._fill_nc(nr, 'thermal', 'NEUTX')
            if hasattr(nr, 'dd'):
                dd = nr.dd
                dd_th = self.im.get_ncvar('NEUTX_DD')
                dd_bt = self.im.get_ncvar('BTNTS_DD')
                dd_bb = self.im.get_ncvar('BBNTS_DD')
                self._fill(dd, 'thermal', dd_th, 'TRANSP CDF:NEUTX_DD')
                self._fill(dd, 'beam_thermal', dd_bt, 'TRANSP CDF:BTNTS_DD')
                self._fill(dd, 'beam_beam', dd_bb, 'TRANSP CDF:BBNTS_DD')
                if self.im.is_valid(dd_th, dd_bt, dd_bb):
                    self._fill(dd, 'total', dd_th + dd_bt + dd_bb,
                               'TRANSP DD neutron components sum')

        self.im.putdb_imas(sid)


if __name__ == '__main__':
    print('Please use: python convert.py')
