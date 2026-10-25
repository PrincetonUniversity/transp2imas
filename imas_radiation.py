#!/usr/bin/python
"""Populate the IMAS radiation IDS from TRANSP radiated-power components."""

import re
import numpy as np
from nc2imas import *
from utilities import append_successful_entries
from writer_base import IdsWriter


class imas_radiation(IdsWriter):
    """Translate spatial TRANSP radiation outputs to the IMAS radiation IDS.

    Radiation quantities in this IDS are emitted power/emissivity and are
    therefore written with a positive sign.  The corresponding loss terms in
    core_sources retain the opposite (negative source) sign.
    """

    PROCESS_IDS = {
        'bremsstrahlung': (8, 'Emission from bremsstrahlung'),
        'synchrotron_radiation': (9, 'Emission from synchrotron radiation'),
        'line_radiation': (10, 'Emission from line radiation'),
        'recombination': (11, 'Emission from recombination'),
    }

    def _load_grid(self):
        """Read the time and radial grids the processes are written on (done in ``build``, not at construction)."""
        self.time = np.asarray(self.cdf.get_ncvar('TIME3'), dtype=np.float64)
        self.x = np.asarray(self.cdf.get_ncvar('X'), dtype=np.float64)
        self.xb = np.asarray(self.cdf.get_ncvar('XB'), dtype=np.float64)
        self.dvol = np.asarray(self.cdf.get_ncvar('DVOL', self.grid.vconv), dtype=np.float64)
        self.nt = len(self.time)
        self.nx = self.x.shape[1] if self.x.ndim == 2 else 0

    def _has(self, name):
        return self.cdf.check_var_name(name)

    def _get(self, name, mult=None):
        arr = self.cdf.get_ncvar(name)
        if not self.fill.is_valid(arr):
            return np.array([])
        arr = np.asarray(arr, dtype=np.float64)
        arr = np.asarray(self.check_profile(name, arr, self.nt, self.nx, grid='X'))   # wrong shape or coordinates: broken mapping
        return arr if mult is None else arr * mult

    @staticmethod
    def _sum_valid(*arrays):
        good = [np.asarray(a, dtype=np.float64) for a in arrays
                if isinstance(a, np.ndarray) and a.size]
        if not good:
            return np.array([])
        out = np.zeros_like(good[0], dtype=np.float64)
        for arr in good:
            if arr.shape == out.shape:
                out += arr
        return out

    def _inside_boundary(self, emissivity):
        """Cumulative emitted power at TRANSP XB zone boundaries."""
        if not self.fill.is_valid(emissivity, self.dvol):
            return np.array([])
        emissivity = np.asarray(emissivity, dtype=np.float64)
        if emissivity.shape != self.dvol.shape:
            return np.array([])
        return np.cumsum(emissivity * self.dvol, axis=1)

    def _inside(self, emissivity):
        """Cumulative emitted power interpolated from XB boundaries to X."""
        inside_b = self._inside_boundary(emissivity)
        if not inside_b.size:
            return np.array([])
        return self.grid.interp_radial_profiles(inside_b, self.xb, self.x, axis_value=0.0)

    def _integral(self, emissivity):
        inside_b = self._inside_boundary(emissivity)
        return inside_b[:, -1] if inside_b.size else np.array([])

    def _set_process_identifier(self, process, name):
        idx, desc = self.PROCESS_IDS[name]
        process.identifier.name = name
        process.identifier.index = idx
        process.identifier.description = desc
        if hasattr(process, 'code'):
            process.code.name = 'TRANSP'
            if hasattr(process.code, 'description'):
                process.code.description = desc

    def _prepare_process(self, process, name):
        self._set_process_identifier(process, name)
        process.global_quantities.resize(self.nt)
        process.profiles_1d.resize(self.nt)
        self.fill.fill_profiles1d_grid(process)
        for it in range(self.nt):
            process.global_quantities[it].time = float(self.time[it])
            process.profiles_1d[it].time = float(self.time[it])

    def _fill_globals(self, process, total, ion=None, neutral=None, electron=None):
        total = np.asarray(total) if isinstance(total, np.ndarray) else np.array([])
        ion = np.asarray(ion) if isinstance(ion, np.ndarray) else np.array([])
        neutral = np.asarray(neutral) if isinstance(neutral, np.ndarray) else np.array([])
        electron = np.asarray(electron) if isinstance(electron, np.ndarray) else np.array([])
        for it in range(self.nt):
            g = process.global_quantities[it].inside_lcfs
            if total.size:
                g.power = float(total[it])
            if ion.size:
                g.power_ion_total = float(ion[it])
            if neutral.size:
                g.power_neutral_total = float(neutral[it])
            if electron.size:
                g.power_electrons = float(electron[it])

    def _fill_electron_process(self, process, name, var):
        emissivity = self._get(var, self.grid.vol_inv)
        if not emissivity.size:
            return False
        self._prepare_process(process, name)
        inside = self._inside(emissivity)
        power = self._integral(emissivity)
        self._fill_globals(process, power, electron=power)
        for it in range(self.nt):
            p = process.profiles_1d[it]
            self.fill.fill_data(p.electrons.emissivity, emissivity, ind=it, chk=self.nx)
            self.fill.fill_data(p.electrons.power_inside, inside, ind=it, chk=self.nx)
        self.logger.info(f'radiation: {name} <- {var}')
        return True

    def _atom_for_token(self, token):
        token = token.upper()
        token = {'3': 'HE3', '4': 'HE4'}.get(token, token)
        try:
            return token, atom_base(token)
        except KeyError:
            self.ctx.mapping.skipped(f'unknown species token {token!r}', target='process/ion')
            return token, None

    def _fill_element_identity(self, node, token, z_ion=None):
        key, atom = self._atom_for_token(token)
        if atom is None:
            self.fill.set_species_name(node, key)
            if z_ion is not None:
                node.z_ion = float(z_ion)
            return
        node.element.resize(1)
        node.element[0].a = float(atom.a)
        self.fill.set_nuclear_charge(node.element[0], atom.z_n)
        node.element[0].atoms_n = int(atom.atoms_n)
        if z_ion is not None:
            node.z_ion = float(z_ion)
        self.fill.set_species_name(node, atom.label.strip())

    def _single_impurity_token(self):
        """Infer the single impurity element from TRANSP AIMP/XZIMP.

        Older/single-impurity TRANSP output often labels line radiation as
        PRLS_TOK / PRLX_SINGL rather than including the chemical symbol in the
        CDF variable name.  The run namelist still carries the impurity mass
        and nuclear charge.  Resolve that pair against ``atom_base.entries``
        only when it identifies a unique physical element closely enough.
        """
        if not getattr(self.namelist, 'nmlist', None):
            return None

        def scalar(name):
            key = next((k for k in self.namelist.nmlist
                        if str(k).strip().upper() == name), None)
            if key is None:
                return None
            value = self.namelist.nmlist[key]
            if isinstance(value, (list, tuple, np.ndarray)):
                if len(value) != 1:
                    return None
                value = value[0]
            if isinstance(value, str):
                value = re.split(r'[\s,;]+', value.strip())[0]
            try:
                return float(self.namelist.parse_number(value, float))
            except (ValueError, TypeError):
                try:
                    return float(str(value).replace('D', 'E').replace('d', 'e'))
                except (ValueError, TypeError):
                    return None

        aimp = scalar('AIMP')
        zimp = scalar('XZIMP')
        if aimp is None or zimp is None:
            return None

        base = atom_base()
        # Exclude TRANSP aliases/model labels and beam/source-state entries.
        excluded = {'TOK', 'Z', 'A', 'N', 'K', 'X'}
        candidates = []
        for key, (_label, aa, zz) in base.entries.items():
            if key in excluded or '_' in key:
                continue
            if abs(float(zz) - zimp) > 0.1:
                continue
            da = abs(float(aa) - aimp)
            if da <= max(0.25, 0.03 * max(abs(aimp), 1.0)):
                candidates.append((da, len(key), key))
        if not candidates:
            return None
        candidates.sort()
        best = candidates[0]
        # Nuclear charge is the primary discriminator; the mass tolerance above
        # admits natural/isotopic rounding while still rejecting unrelated ions.
        return best[2]

    def _line_inventory(self):
        """Return species-total PRLS_* and charge-state PRLX_* CDF names.

        PRLS_<element> is treated as the element-integrated impurity line
        radiation profile.  PRLX_<element>_<charge> supplies charge-state
        detail below that element.  Generic model-family totals such as TOK and
        SINGL are not treated as chemical symbols directly; for a single-
        impurity run they can be assigned to the element inferred from the
        namelist AIMP/XZIMP pair.
        """
        totals = {}
        states = {}
        ignored = {'TOK', 'SINGL'}
        for name in getattr(self.cdf, 'vars_', []):
            m = re.fullmatch(r'PRLS_([A-Z0-9]+)', name)
            if m:
                sp = m.group(1)
                if sp not in ignored:
                    _key, atom = self._atom_for_token(sp)
                    if atom is not None:
                        totals[sp] = name
                continue
            m = re.fullmatch(r'PRLX_([A-Z0-9]+)_([0-9]+)', name)
            if m:
                sp, charge = m.group(1), int(m.group(2))
                if sp not in ignored:
                    _key, atom = self._atom_for_token(sp)
                    if atom is not None:
                        states.setdefault(sp, []).append((charge, name))
        for sp in states:
            states[sp].sort(key=lambda item: item[0])

        # Single-impurity runs can use generic TOK/SINGL labels.  When AIMP and
        # XZIMP uniquely identify the element, attach exactly one of those alias
        # profiles to that species.  Prefer PRLS_TOK; PRLX_SINGL is commonly an
        # identical alias and must not be double counted.
        inferred = self._single_impurity_token()
        if inferred and inferred not in totals and inferred not in states:
            if self._has('PRLS_TOK'):
                totals[inferred] = 'PRLS_TOK'
            elif self._has('PRLX_SINGL'):
                totals[inferred] = 'PRLX_SINGL'
        return totals, states

    def _line_total_from_species(self, totals, states):
        pieces = []
        species = sorted(set(totals) | set(states))
        for sp in species:
            if sp in totals:
                arr = self._get(totals[sp], self.grid.vol_inv)
            else:
                arr = self._sum_valid(*[self._get(name, self.grid.vol_inv)
                                        for _z, name in states.get(sp, [])])
            if arr.size:
                pieces.append(arr)
        return self._sum_valid(*pieces)

    def _fill_line_radiation(self, process):
        totals, states = self._line_inventory()
        total_emissivity = self._get('PRAD_LI', self.grid.vol_inv)
        if not total_emissivity.size:
            total_emissivity = self._line_total_from_species(totals, states)
        if not total_emissivity.size:
            return False

        self._prepare_process(process, 'line_radiation')
        total_inside = self._inside(total_emissivity)
        total_power = self._integral(total_emissivity)
        # TRANSP PRAD_LI is the total line-radiation loss.  In the radiation IDS
        # it is represented as total ion/impurity emission; detailed PRLS/PRLX
        # profiles are attached to the corresponding isonuclear species below.
        self._fill_globals(process, total_power, ion=total_power)

        species = sorted(set(totals) | set(states))
        species_data = []
        for sp in species:
            state_data = [(z, name, self._get(name, self.grid.vol_inv))
                          for z, name in states.get(sp, [])]
            state_data = [(z, name, arr) for z, name, arr in state_data if arr.size]
            if sp in totals:
                sp_total = self._get(totals[sp], self.grid.vol_inv)
            else:
                sp_total = self._sum_valid(*[arr for _z, _name, arr in state_data])
            if sp_total.size or state_data:
                species_data.append((sp, sp_total, state_data))

        # Everything below that does not depend on the time slice is computed once: the enclosed power of each
        # species and charge state (a cumulative sum and a radial interpolation of the whole time-by-radius
        # array) and the species label.  They used to be recomputed for every slice.
        species_inside = []
        for sp, sp_total, state_data in species_data:
            _, atom = self._atom_for_token(sp)
            species_inside.append({
                'ion': self._inside(sp_total) if sp_total.size else None,
                'base': atom.label.strip() if atom is not None else sp,
                'states': [self._inside(state_prof) for _charge, _name, state_prof in state_data]})

        for it in range(self.nt):
            p = process.profiles_1d[it]
            self.fill.fill_data(p.emissivity_ion_total, total_emissivity, ind=it, chk=self.nx)
            self.fill.fill_data(p.power_inside_ion_total, total_inside, ind=it, chk=self.nx)
            p.ion.resize(len(species_data))
            for ji, (sp, sp_total, state_data) in enumerate(species_data):
                ion = p.ion[ji]
                cached = species_inside[ji]
                self._fill_element_identity(ion, sp)
                ion.multiple_states_flag = 1 if state_data else 0
                if sp_total.size:
                    self.fill.fill_data(ion.emissivity, sp_total, ind=it, chk=self.nx)
                    self.fill.fill_data(ion.power_inside, cached['ion'], ind=it, chk=self.nx)

                ion.state.resize(len(state_data))
                for js, (charge, _name, state_prof) in enumerate(state_data):
                    state = ion.state[js]
                    state.z_min = float(charge)
                    state.z_max = float(charge)
                    self.fill.set_species_name(state, f"{cached['base']}+{charge}")
                    self.fill.fill_data(state.emissivity, state_prof, ind=it, chk=self.nx)
                    self.fill.fill_data(state.power_inside, cached['states'][js], ind=it, chk=self.nx)

        detail = []
        if totals:
            detail.append(f'{len(totals)} PRLS species totals')
        if states:
            detail.append(f'{sum(len(v) for v in states.values())} PRLX charge states')
        suffix = ', '.join(detail) if detail else 'no species-resolved profiles'
        self.logger.info(f'radiation: line_radiation <- PRAD_LI ({suffix})')
        return True


    def build(self):
        self._load_grid()
        self.logger.info('radiation')
        ids = self.dd.new_ids('radiation')
        ids = self.fill.fill_properties(ids)
        self.fill.fill_code_properties(ids)
        ids.ids_properties.homogeneous_time = 1
        ids.time = np.asarray(self.time, dtype=np.float64)
        self.fill.fill_vac_tor_field(ids)

        builders = [
            lambda p: self._fill_electron_process(p, 'bremsstrahlung', 'PRAD_BR'),
            lambda p: self._fill_electron_process(p, 'synchrotron_radiation', 'PRAD_CY'),
            self._fill_line_radiation,
        ]
        used = append_successful_entries(lambda: self.dd.new_ids('radiation'), ids, 'process', builders)

        if used:
            self.logger.info(f'radiation: built {used} emission processes')
            return [ids]
        self.logger.info('radiation: no supported CDF radiation profiles found; IDS not written')
        return []


if __name__ == '__main__':
    print('Please use: python convert.py')
