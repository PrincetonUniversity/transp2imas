#!/usr/bin/python
"""Populate the IMAS waves IDS when spatial RF absorption/deposition data are available."""

import re
import numpy as np

from nc2imas import nc2imas


class imas_waves:
    """Translate TRANSP RF deposition/current-drive output to the IMAS waves IDS.

    The writer is deliberately conservative.  It only creates a ``waves`` IDS
    when TRANSP contains a supported *absorbed/deposited* EC or IC result.  RF
    launch-only quantities belong in the antenna/launcher IDSs and are not
    relabelled as absorbed wave power here.

    Supported mappings
    ------------------
    EC (TORAY / GENRAY / TORBEAM)
      * ``PEECHj(TIME3,X)`` -> total and thermal-electron power density.
      * ``ECCURj(TIME3,X)`` -> toroidal current inside each flux surface.
      * ``FREQECH(j)`` from TR.DAT -> wave frequency.
      * For a single EC wave only, aggregate ``PEECH`` / ``ECCUR`` are accepted
        as fallbacks when the per-gyrotron arrays are absent.

    IC (TORIC-style full-wave output)
      * ``RFPWRABj(TIME)`` -> absorbed wave power for antenna/frequency j.
      * ``FREQAj(TIME)`` -> wave frequency for antenna/frequency j.
      * ``PIC_Fj(TIME3,X)`` -> antenna/frequency-resolved absorbed wave power
        density.  Its volume integral is checked by the source-side regression
        against ``RFPWRABj``.
      * For a single IC coherent wave only, aggregate ``QICHA(TIME3,X)`` is an
        accepted fallback total absorbed-power profile when ``PIC_Fj`` is not
        available.
      * For a single IC coherent wave only, ``ICCUR(TIME3,X)`` is mapped as
        cumulative toroidal current.  Multi-wave aggregate current is not
        duplicated across coherent waves.
      * ``XGRID_NPHI`` is exported only for a single IC wave because the CDF
        spectral grid is not frequency/antenna resolved.

    ``PEICH``, ``PIICH`` and ``QRFBTOT`` are *not* wave-absorption partitions.
    They are plasma-heating/source quantities and are handled by
    ``core_sources``.  In the C-Mod validation case their sum strongly
    over-counts the absorbed RF power, while ``PIC_Fj`` and ``QICHA`` close to
    ``RFPWRABj`` at numerical precision.

    TRANSP EC current density is toroidal current density (A/cm^2), not the
    IMAS definition of ``current_parallel_density = <j.B>/B0``.  Therefore the
    cumulative toroidal current is stored in ``current_phi_inside`` (or the
    pre-3.42 ``current_tor_inside``) and ``current_parallel_density`` is left
    empty.
    """

    W_PER_CM3_TO_W_PER_M3 = 1.0e6

    def __init__(self, ncimas: nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.time = np.asarray(self.im.get_ncvar('TIME3'), dtype=np.float64)
        self.scalar_time = np.asarray(self.im.get_ncvar('TIME'), dtype=np.float64)
        self.rho = np.asarray(self.im.get_ncvar('X'), dtype=np.float64)
        self.dvol_cm3 = np.asarray(self.im.get_ncvar('DVOL'), dtype=np.float64)
        self.darea_cm2 = np.asarray(self.im.get_ncvar('DAREA'), dtype=np.float64)
        self.nt = len(self.time)
        self.nx = self.rho.shape[1] if self.rho.ndim == 2 else 0
        self.create_ids()

    # ------------------------------------------------------------------
    # CDF / namelist inventory helpers
    # ------------------------------------------------------------------
    def _has(self, name):
        return bool(self.im.check_var_name(name))

    def _profile(self, name, mult=1.0):
        """Return a TIME3 x X profile, otherwise an empty array."""
        if not self._has(name):
            return np.array([])
        arr = np.asarray(self.im.get_ncvar(name), dtype=np.float64)
        if self.rho.ndim != 2 or arr.shape != self.rho.shape:
            return np.array([])
        return arr * float(mult)

    def _series(self, name, target_time=None):
        """Return a scalar-time CDF quantity interpolated to ``target_time``."""
        if not self._has(name):
            return np.array([])
        arr = np.asarray(self.im.get_ncvar(name), dtype=np.float64)
        if arr.ndim == 0:
            arr = np.asarray([float(arr)], dtype=np.float64)
        if arr.ndim != 1 or arr.size == 0:
            return np.array([])
        if target_time is None:
            target_time = self.time
        target_time = np.asarray(target_time, dtype=np.float64)
        if arr.size == 1:
            return np.full(len(target_time), float(arr[0]), dtype=np.float64)

        dims = ()
        try:
            dims = tuple(self.im.nc.variables[name].dimensions)
        except (AttributeError, KeyError, TypeError):
            pass
        if dims and dims[0] == 'TIME3' and len(self.time) == len(arr):
            source_time = self.time
        elif dims and dims[0] == 'TIME' and len(self.scalar_time) == len(arr):
            source_time = self.scalar_time
        elif len(self.time) == len(arr):
            source_time = self.time
        elif len(self.scalar_time) == len(arr):
            source_time = self.scalar_time
        else:
            return np.array([])

        if len(source_time) == len(target_time) and np.allclose(source_time, target_time):
            return arr
        return np.interp(target_time, source_time, arr)

    @staticmethod
    def _indices_from_names(names, prefixes):
        out = set()
        for name in names:
            for prefix in prefixes:
                match = re.fullmatch(rf'{re.escape(prefix)}(\d+)', str(name))
                if match:
                    out.add(int(match.group(1)))
        return sorted(out)

    def _namelist_lookup(self, name):
        """Case-insensitive raw lookup of NAME and NAME(i) in the TR.DAT map."""
        if not getattr(self.im, 'nmlist', None):
            return None, []
        target = str(name).upper()
        base = None
        indexed = []
        for key, value in self.im.nmlist.items():
            kup = str(key).strip().upper()
            if kup == target:
                base = value
                continue
            match = re.fullmatch(rf'{re.escape(target)}\((\d+)\)', kup)
            if match:
                indexed.append((int(match.group(1)), value))
        indexed.sort(key=lambda item: item[0])
        return base, indexed

    def _namelist_numbers(self, name):
        base, indexed = self._namelist_lookup(name)

        def parse_one(value):
            if isinstance(value, (int, float, np.integer, np.floating)):
                return [float(value)]
            if value is None:
                return []
            text = str(value).strip()
            vals = []
            for token in re.split(r'[\s,;]+', text):
                token = token.strip().strip("'\"")
                if not token:
                    continue
                try:
                    vals.append(float(self.im.parse_number(token, float)))
                except (ValueError, TypeError, ZeroDivisionError, IndexError):
                    try:
                        vals.append(float(token.replace('D', 'E').replace('d', 'e')))
                    except ValueError:
                        pass
            return vals

        if base is not None:
            return np.asarray(parse_one(base), dtype=np.float64)
        if indexed:
            vals = []
            for _, value in indexed:
                parsed = parse_one(value)
                vals.append(parsed[0] if parsed else np.nan)
            return np.asarray(vals, dtype=np.float64)
        return np.array([])

    def _namelist_number(self, name, one_based_index=1):
        vals = self._namelist_numbers(name)
        if vals.size == 0:
            return None
        if vals.size == 1:
            return float(vals[0]) if np.isfinite(vals[0]) else None
        idx = int(one_based_index) - 1
        if 0 <= idx < vals.size and np.isfinite(vals[idx]):
            return float(vals[idx])
        return None

    def _namelist_text(self, name, default=''):
        base, indexed = self._namelist_lookup(name)
        value = base
        if value is None and indexed:
            value = indexed[0][1]
        if value is None:
            return default
        return str(value).strip().strip("'\"").upper()

    def _namelist_int(self, name, default=0):
        value = self._namelist_number(name, 1)
        if value is None:
            return int(default)
        try:
            return int(round(value))
        except (TypeError, ValueError, OverflowError):
            return int(default)

    def _discover_ec_indices(self):
        names = getattr(self.im, 'vars_', [])
        indices = self._indices_from_names(names, ('PEECH', 'ECCUR'))
        # PECIN is launch power only.  It can identify a gyrotron when a single
        # aggregate deposition/current profile is the only spatial result.
        launch_indices = self._indices_from_names(names, ('PECIN',))
        if indices:
            return indices

        has_total = bool(self._profile('PEECH').size or self._profile('ECCUR').size)
        if not has_total:
            return []

        nantech = self._namelist_int('NANTECH', 0)
        freq = self._namelist_numbers('FREQECH')
        candidates = launch_indices
        if not candidates and nantech == 1:
            candidates = [1]
        if not candidates and freq.size == 1:
            candidates = [1]
        # Aggregate PEECH/ECCUR cannot be split safely among several gyrotrons.
        return [1] if candidates == [1] else []

    def _discover_ic_indices(self):
        names = getattr(self.im, 'vars_', [])
        indices = self._indices_from_names(names, ('RFPWRAB', 'FREQA', 'PICHA', 'PIC_F'))
        # Require an absorbed-power trace/profile or, for a single antenna, an
        # aggregate total-heating profile. PICHA/FREQA alone are launcher data.
        absorbed = [i for i in indices
                    if self._has(f'RFPWRAB{i}') or self._profile(f'PIC_F{i}').size]
        if absorbed:
            return absorbed
        has_profile = bool(self._profile('QICHA').size)
        if has_profile and len(indices) == 1:
            return indices
        return []

    # ------------------------------------------------------------------
    # Numerical conversions
    # ------------------------------------------------------------------
    def _power_inside(self, power_density_w_m3):
        """Cumulative volume integral on the TRANSP finite-volume X grid."""
        p = np.asarray(power_density_w_m3, dtype=np.float64)
        if p.ndim != 2 or self.dvol_cm3.shape != p.shape:
            return np.array([])
        dvol_m3 = self.dvol_cm3 * 1.0e-6
        return np.cumsum(p * dvol_m3, axis=1)

    def _current_inside(self, current_density_a_cm2):
        """Cumulative toroidal current from TRANSP current density."""
        j = np.asarray(current_density_a_cm2, dtype=np.float64)
        if j.ndim != 2 or self.darea_cm2.shape != j.shape:
            return np.array([])
        sign = float(getattr(self.im, 'jphi_sign', getattr(self.im, 'ip_sign', 1.0)))
        return sign * np.cumsum(j * self.darea_cm2, axis=1)

    def _nphi_grid(self):
        if not self._has('XGRID_NPHI'):
            return np.array([])
        arr = np.asarray(self.im.get_ncvar('XGRID_NPHI'))
        if arr.ndim == 1:
            arr = np.tile(arr[None, :], (self.nt, 1))
        if arr.ndim != 2 or arr.shape[0] not in (1, self.nt):
            return np.array([])
        if arr.shape[0] == 1 and self.nt > 1:
            arr = np.tile(arr, (self.nt, 1))
        # The DD node is integer mode number.  Ignore non-finite padding.
        return arr

    def _set_nphi(self, node, values):
        values = np.asarray(values)
        if values.ndim != 1 or values.size == 0:
            return
        good = np.isfinite(values)
        vals = np.rint(values[good]).astype(np.int64)
        if vals.size == 0:
            return
        field = 'n_phi' if self.im.dd_at_least(3, 42, 0) else 'n_tor'
        target = getattr(node, field)
        target.resize(len(vals), refcheck=False)
        target[:] = vals

    def _set_current_global(self, glob, value):
        field = 'current_phi' if self.im.dd_at_least(3, 42, 0) else 'current_tor'
        setattr(glob, field, float(value))

    def _set_current_inside(self, prof, values):
        field = 'current_phi_inside' if self.im.dd_at_least(3, 42, 0) else 'current_tor_inside'
        target = getattr(prof, field)
        vals = np.asarray(values, dtype=np.float64)
        target.resize(len(vals), refcheck=False)
        target[:] = vals

    @staticmethod
    def _fill_vector(target, values):
        vals = np.asarray(values, dtype=np.float64)
        if vals.ndim != 1 or vals.size == 0:
            return
        target.resize(len(vals), refcheck=False)
        target[:] = vals

    @staticmethod
    def _set_identifier(identifier, name, index, description):
        identifier.name = str(name)
        identifier.index = int(index)
        identifier.description = str(description)

    # ------------------------------------------------------------------
    # Wave builders
    # ------------------------------------------------------------------
    def _fill_global_common(self, wave, frequency, power, current, electron_power, nphi):
        wave.global_quantities.resize(self.nt)
        for it in range(self.nt):
            glob = wave.global_quantities[it]
            glob.time = float(self.time[it])
            if frequency is not None:
                if np.isscalar(frequency):
                    if np.isfinite(frequency):
                        glob.frequency = float(frequency)
                elif len(frequency) > it and np.isfinite(frequency[it]):
                    glob.frequency = float(frequency[it])
            if isinstance(power, np.ndarray) and power.size > it and np.isfinite(power[it]):
                glob.power = float(power[it])
            if isinstance(current, np.ndarray) and current.size > it and np.isfinite(current[it]):
                self._set_current_global(glob, current[it])
            if isinstance(electron_power, np.ndarray) and electron_power.size > it and np.isfinite(electron_power[it]):
                glob.electrons.power_thermal = float(electron_power[it])
                glob.electrons.distribution_assumption = 0
            if isinstance(nphi, np.ndarray) and nphi.ndim == 2 and nphi.shape[0] > it:
                self._set_nphi(glob, nphi[it])

    def _fill_profiles(self, wave, total_pd, electron_pd=None, current_inside=None, nphi=None):
        has_total = isinstance(total_pd, np.ndarray) and total_pd.size
        has_e = isinstance(electron_pd, np.ndarray) and electron_pd.size
        has_j = isinstance(current_inside, np.ndarray) and current_inside.size
        has_nphi = isinstance(nphi, np.ndarray) and nphi.size
        if not (has_total or has_e or has_j or has_nphi):
            return

        wave.profiles_1d.resize(self.nt)
        self.im.fill_profiles1d_grid(wave)
        total_inside = self._power_inside(total_pd) if has_total else np.array([])
        e_inside = self._power_inside(electron_pd) if has_e else np.array([])

        for it in range(self.nt):
            prof = wave.profiles_1d[it]
            prof.time = float(self.time[it])
            if has_total:
                self._fill_vector(prof.power_density, total_pd[it])
                self._fill_vector(prof.power_inside, total_inside[it])
            if has_e:
                self._fill_vector(prof.electrons.power_density_thermal, electron_pd[it])
                self._fill_vector(prof.electrons.power_inside_thermal, e_inside[it])
                # distribution_assumption is a global_quantities.electrons
                # field, not a profiles_1d.electrons field (DD3 and DD4).
            if has_j:
                self._set_current_inside(prof, current_inside[it])
            if has_nphi and nphi.ndim == 2 and nphi.shape[0] > it:
                self._set_nphi(prof, nphi[it])

    def _build_ec_wave(self, wave, index, ec_indices):
        self._set_identifier(
            wave.identifier.type, 'EC', 1,
            'Wave field for electron cyclotron heating and current drive')
        wave.identifier.antenna_name = f'beam{index}'
        wave.identifier.index_in_antenna = 1

        model = self._namelist_text('EC_MODEL', 'TRANSP_EC')
        model_name = model if model and model != 'NONE' else 'TRANSP_EC'
        self._set_identifier(
            wave.wave_solver_type, model_name, 1,
            'Beam/ray-tracing EC wave deposition solver')

        # Per-gyrotron profiles are authoritative.  Aggregate fallbacks are only
        # safe when exactly one EC wave exists.
        pd = self._profile(f'PEECH{index}', self.W_PER_CM3_TO_W_PER_M3)
        j = self._profile(f'ECCUR{index}')
        if len(ec_indices) == 1:
            if not pd.size:
                pd = self._profile('PEECH', self.W_PER_CM3_TO_W_PER_M3)
            if not j.size:
                j = self._profile('ECCUR')

        current_inside = self._current_inside(j) if j.size else np.array([])
        power_inside = self._power_inside(pd) if pd.size else np.array([])
        power = power_inside[:, -1] if power_inside.size else np.array([])
        current = current_inside[:, -1] if current_inside.size else np.array([])
        frequency = self._namelist_number('FREQECH', index)

        self._fill_global_common(
            wave, frequency=frequency, power=power, current=current,
            electron_power=power, nphi=np.array([]))
        self._fill_profiles(
            wave, total_pd=pd, electron_pd=pd,
            current_inside=current_inside, nphi=np.array([]))

    def _build_ic_wave(self, wave, index, ic_indices, nphi):
        # The dedicated wave_identifier translation table is authoritative:
        # EC=1, LH=2, IC=3.  Some prose in older waves pages listed IC/LH in the
        # opposite order.
        self._set_identifier(
            wave.identifier.type, 'IC', 3,
            'Wave field for ion cyclotron frequency heating and current drive')
        wave.identifier.antenna_name = f'antenna{index}'
        wave.identifier.index_in_antenna = 1
        self._set_identifier(
            wave.wave_solver_type, 'TORIC/full_wave', 2,
            'Full-wave IC deposition (TORIC-style TRANSP output)')

        frequency = self._series(f'FREQA{index}', self.time)
        power = self._series(f'RFPWRAB{index}', self.time)
        total_pd = self._profile(f'PIC_F{index}', self.W_PER_CM3_TO_W_PER_M3)
        current_inside = np.array([])
        current = np.array([])

        # PIC_Fj is the frequency/antenna-resolved wave absorption profile and
        # is the preferred spatial source for coherent_wave j.  If only a
        # single IC wave is represented, QICHA is an acceptable aggregate
        # fallback because it is the total ICH heating profile and closes to
        # the absorbed antenna power in the validation case.
        if len(ic_indices) == 1:
            if not total_pd.size:
                total_pd = self._profile('QICHA', self.W_PER_CM3_TO_W_PER_M3)
            j = self._profile('ICCUR')
            if j.size:
                current_inside = self._current_inside(j)
                current = current_inside[:, -1]

        if total_pd.size and not power.size:
            power = self._power_inside(total_pd)[:, -1]

        nphi_this = nphi if len(ic_indices) == 1 else np.array([])
        self._fill_global_common(
            wave, frequency=frequency, power=power, current=current,
            electron_power=np.array([]), nphi=nphi_this)
        self._fill_profiles(
            wave, total_pd=total_pd, electron_pd=np.array([]),
            current_inside=current_inside, nphi=nphi_this)

    # ------------------------------------------------------------------
    # IDS construction
    # ------------------------------------------------------------------
    def create_ids(self):
        self.logger.info('\nwaves')
        if self.nt == 0 or self.rho.ndim != 2 or self.dvol_cm3.shape != self.rho.shape:
            self.logger.info('waves: invalid TIME3/X/DVOL profile grid; IDS not written')
            return

        ec_indices = self._discover_ec_indices()
        ic_indices = self._discover_ic_indices()
        if not ec_indices and not ic_indices:
            self.logger.info('waves: no supported EC/IC absorbed/deposited wave output found; IDS not written')
            return

        ids = self.im.new_ids('waves')
        ids.coherent_wave.resize(len(ec_indices) + len(ic_indices))
        iw = 0
        for index in ec_indices:
            self._build_ec_wave(ids.coherent_wave[iw], index, ec_indices)
            iw += 1

        nphi = self._nphi_grid()
        for index in ic_indices:
            self._build_ic_wave(ids.coherent_wave[iw], index, ic_indices, nphi)
            iw += 1

        self.im.fill_vac_tor_field(ids)
        ids = self.im.fill_properties(ids)
        self.im.fill_code_properties(ids)
        self.im.fill_external_model_inputs(ids, 'waves')
        ids.ids_properties.homogeneous_time = 1
        ids.time.resize(self.nt)
        ids.time[:] = self.time
        self.im.putdb_imas(ids)
        self.logger.info(f'Done: {len(ec_indices)} EC + {len(ic_indices)} IC coherent wave(s)')


if __name__ == '__main__':
    print('Please use: python convert.py ')
