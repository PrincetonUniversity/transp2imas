#!/usr/bin/python
"""Populate the IMAS summary IDS from global TRANSP time traces and simulation metadata."""

import numpy as np
from scipy import constants

from writer_base import IdsWriter


class imas_summary(IdsWriter):
    """Write a curated TRANSP summary IDS from high-confidence CDF quantities."""

    def _to_summary_time(self, data):
        """Map a 1-D profile-time quantity onto summary/TIME."""
        arr = np.asarray(data)
        if arr.size == 0:
            return arr
        if arr.ndim == 0:
            return arr
        if arr.ndim != 1:
            return np.array([])
        if len(arr) == self.grid.ntimes:
            return arr.astype(np.float64, copy=False)
        if len(arr) == self.grid.ntimep and self.fill.is_valid(self.grid.timep, self.grid.timesc):
            return np.interp(self.grid.timesc, self.grid.timep, arr)
        return arr

    def _integrate_profile(self, variable, multiplier=1.0):
        """Volume integral for a TRANSP X-grid profile using DVOL in cm^3."""
        val = self.profile(variable)
        dvol = self.profile('DVOL')
        if not self.fill.is_valid(val, dvol) or val.ndim != 2 or val.shape != dvol.shape:
            return np.array([])
        out = multiplier * np.sum(val * dvol, axis=1)
        return self._to_summary_time(out)

    def _fill(self, parent, name, data, source=None, scalar=False):
        if not hasattr(parent, name) or not self.fill.is_valid(data):
            return
        node = getattr(parent, name)
        if scalar:
            self.fill.fill_value(node, data, rsz=0)
        else:
            self.fill.fill_value(node, self._to_summary_time(data))
        if source and hasattr(node, 'source'):
            node.source = source

    def _fill_nc(self, parent, name, cdf_name, multiplier=1.0, source=None, integer=False):
        data = self.cdf.get_ncvar(cdf_name)
        if not self.fill.is_valid(data):
            return
        data = multiplier * data
        if integer:
            data = np.rint(data).astype(np.int32)
        self._fill(parent, name, data, source or f'TRANSP CDF:{cdf_name}')

    def build(self):
        self.logger.info('\nsummary')
        sid = self.dd.new_ids('summary')
        sid = self.fill.fill_properties(sid)
        # Store the same TRANSP/template provenance in summary/code.parameters.
        sid = self.fill.fill_code_properties(sid)

        sid.ids_properties.homogeneous_time = 1
        times = np.asarray(self.grid.timesc, dtype=np.float64)
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
            sid.description = f'TRANSP integrated simulation for run {self.cdf.device} {self.cdf.shot}'
        if hasattr(sid, 'machine'):
            sid.machine = str(self.cdf.device)
        if hasattr(sid, 'pulse'):
            try:
                sid.pulse = int(self.cdf.shot)
            except (TypeError, ValueError):
                self.ctx.mapping.skipped('shot number is not an integer', source='shot', target='pulse')
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
        ip_sign = self.signs.ip_sign
        bt_sign = self.signs.bt_sign
        q_sign = self.signs.q_sign

        # Basic global quantities.
        self._fill_nc(gq, 'ip', 'PCURC', ip_sign)
        # r0 and b0 follow the convention shared with vacuum_toroidal_field in every IDS: one
        # constant r0 (middle of the limiter) and b0 = bt_sign*BZXR/r0, so r0*b0 equals the TF vacuum field times radius.
        reference = self.signs.vacuum_field_reference(self.grid.timesc)
        if reference is not None:
            r0, b0_series = reference
            self._fill(gq, 'r0', r0, 'TRANSP CDF: middle of the RLIM limiter (else mean of first/last RAXIS)', scalar=True)
            self._fill(gq, 'b0', b0_series, 'TRANSP CDF: bt_sign*BZXR/r0')

        # li was renamed li_3 in DD4.
        li_name = 'li_3' if hasattr(gq, 'li_3') else 'li'
        self._fill_nc(gq, li_name, 'LI_3')
        self._fill_nc(gq, 'v_loop', 'VSURC', self.signs.vloop_sign)
        self._fill_nc(gq, 'volume', 'PVOL', 1.0e-6)
        self._fill_nc(gq, 'q_95', 'Q95', q_sign)
        self._fill_nc(gq, 'tau_energy', 'TAUA1')
        self._fill_nc(gq, 'h_98', 'H98Y2')
        self._fill_nc(gq, 'tau_energy_98', 'TAUE98Y2')
        self._fill_nc(gq, 'h_mode', 'LHMODE', integer=True)

        # The betas by the Data Dictionary definitions (see ``plasma_betas``): against this IDS's own b0 and r0
        # and, for the normalised one, the boundary minor radius; TRANSP's BTEQ/BETAT use other normalisations.
        # The physical signs of B0 and Ip are retained in beta_tor_norm, matching the IMAS convention.
        betas = self.signs.plasma_betas()
        if not betas:
            self.ctx.mapping.skipped('betas need PTOWB, DVOL, PCURC, the vacuum field and the RMCB boundary moments',
                                     target='global_quantities/beta_tor')
        source = 'Derived from PTOWB,DVOL,PCURC,RMCB* with b0 and r0'
        self._fill(gq, 'beta_pol', betas.get('beta_pol', np.array([])), source)
        self._fill(gq, 'beta_tor', betas.get('beta_tor', np.array([])), source)
        beta_name = 'beta_tor_norm' if hasattr(gq, 'beta_tor_norm') else 'beta_normal'
        self._fill(gq, beta_name, betas.get('beta_tor_norm', np.array([])), source)

        # Energy contents. UE/UI/UTHRM/UTOTL are already energy densities in
        # J/cm^3, so their volume integrals are energies in J -- no 3/2 factor.
        self._fill(gq, 'energy_total', self._integrate_profile('UTOTL'), 'Integral TRANSP CDF:UTOTL*DVOL')
        self._fill(gq, 'energy_thermal', self._integrate_profile('UTHRM'), 'Integral TRANSP CDF:UTHRM*DVOL')
        self._fill(gq, 'energy_electrons_thermal', self._integrate_profile('UE'), 'Integral TRANSP CDF:UE*DVOL')
        self._fill(gq, 'energy_ion_total_thermal', self._integrate_profile('UI'), 'Integral TRANSP CDF:UI*DVOL')

        # MHD pressure and diamagnetic energy use pressure definitions from the DD.
        ptowb = self.profile('PTOWB')
        dvol = self.profile('DVOL', 1.0e-6)  # m^3
        if self.fill.is_valid(ptowb, dvol) and ptowb.shape == dvol.shape:
            self._fill(gq, 'energy_mhd', 1.5*np.sum(ptowb*dvol, axis=1),
                       '1.5 integral TRANSP CDF:PTOWB dV')

        pthermal = self.profile('PPLAS')
        ufast_perp = self.profile('UFASTPP', self.grid.vol_inv)  # J/m3 = Pa_perp
        ufast_parallel = self.profile('UFASTPA', self.grid.vol_inv)
        if self.fill.is_valid(pthermal, ufast_perp, dvol):
            pperp = pthermal + ufast_perp
            self._fill(gq, 'energy_diamagnetic', 1.5*np.sum(pperp*dvol, axis=1),
                       '1.5 integral (PPLAS + UFASTPP) dV')
        if self.fill.is_valid(ufast_perp):
            self._fill(gq, 'energy_fast_perpendicular',
                       1.5*self._integrate_profile('UFASTPP'),
                       '1.5 integral TRANSP CDF:UFASTPP*DVOL')
        if self.fill.is_valid(ufast_parallel):
            # UFASTPA = p_parallel/2; DD quantity = 1.5 integral p_parallel dV.
            self._fill(gq, 'energy_fast_parallel',
                       3.0*self._integrate_profile('UFASTPA'),
                       '3 integral TRANSP CDF:UFASTPA*DVOL')

        # Poloidal magnetic-field energy.
        bpol = self.profile('BPOL')
        coord_x = self.cdf.get_ncvar('X')
        coord_xb = self.cdf.get_ncvar('XB')
        if self.fill.is_valid(bpol, dvol, coord_x, coord_xb):
            try:
                bpol_x = self.grid.interp_radial_profiles(bpol, coord_xb, coord_x, axis_value=0.0)
                ebpol = np.sum((bpol_x*bpol_x)/(2.0*constants.mu_0) * dvol, axis=1)
                self._fill(gq, 'energy_b_field_pol', ebpol, 'Integral TRANSP CDF:BPOL^2/(2 mu0) dV')
            except (ValueError, TypeError, IndexError, FloatingPointError) as exc:
                self.ctx.mapping.broken(exc, source='BPOL', target='global_quantities/energy_b_field_pol')

        # Toroidal current components; TRANSP CUR* [A/cm^2] times DAREA [cm^2]
        # gives A directly.
        darea = self.profile('DAREA')
        cur = self.profile('CUR')
        curoh = self.profile('CUROH')
        curbs = self.profile('CURBS')
        if self.fill.is_valid(darea, curoh) and darea.shape == curoh.shape:
            self._fill(gq, 'current_ohm', ip_sign*np.sum(curoh*darea, axis=1),
                       'Integral TRANSP CDF:CUROH*DAREA')
        if self.fill.is_valid(darea, cur, curoh) and cur.shape == darea.shape == curoh.shape:
            self._fill(gq, 'current_non_inductive', ip_sign*np.sum((cur-curoh)*darea, axis=1),
                       'Integral TRANSP CDF:(CUR-CUROH)*DAREA')
        if self.fill.is_valid(darea, curbs) and darea.shape == curbs.shape:
            self._fill(gq, 'current_bootstrap', ip_sign*np.sum(curbs*darea, axis=1),
                       'Integral TRANSP CDF:CURBS*DAREA')

        # Heating/loss powers.
        poh_total = self.cdf.get_ncvar('POHT')
        if self.fill.is_valid(poh_total):
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
                dd_th = self.cdf.get_ncvar('NEUTX_DD')
                dd_bt = self.cdf.get_ncvar('BTNTS_DD')
                dd_bb = self.cdf.get_ncvar('BBNTS_DD')
                self._fill(dd, 'thermal', dd_th, 'TRANSP CDF:NEUTX_DD')
                self._fill(dd, 'beam_thermal', dd_bt, 'TRANSP CDF:BTNTS_DD')
                self._fill(dd, 'beam_beam', dd_bb, 'TRANSP CDF:BBNTS_DD')
                if self.fill.is_valid(dd_th, dd_bt, dd_bb):
                    self._fill(dd, 'total', dd_th + dd_bt + dd_bb,
                               'TRANSP DD neutron components sum')

        return [sid]


if __name__ == '__main__':
    print('Please use: python convert.py')
