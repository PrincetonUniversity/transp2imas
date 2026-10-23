#!/usr/bin/python
"""Populate the IMAS edge_sources IDS from TRANSP neutral and edge source terms."""

import re
import numpy as np
from nc2imas import *
from utilities import append_successful_entries


class imas_edge_sources:
    """Translate TRANSP neutral gas/recycling source terms to edge_sources.

    TRANSP's ``N0SGF_*`` and ``N0SRC_*`` arrays are *neutral densities* and
    must not be stored in edge_sources particle-source nodes.  This writer uses
    the direct source-rate arrays.  Electron-impact ionization and fast-ion
    charge exchange are kept as separate IMAS source processes using the
    corresponding ``SE*``/``PI*``/``TQI*`` and ``SFCX*``/``PCX*``/``TQCX*``
    diagnostics.  ``SIGF_*_*``/``SIRC_*_*`` remain useful consistency checks:
    in current TRANSP output they equal the ionization plus CX ion source.

    A one-dimensional GGD is described explicitly on TRANSP's ``X`` grid,
    interpreted as rho_tor_norm, and every source quantity points to the
    implicit ``nodes`` grid subset.
    """

    GRID_INDEX = 1
    GRID_SUBSET_INDEX = 1
    CM3_TO_M3 = 1.0e6

    ORIGINS = {
        'gas_flow': {
            'electron_prefix': 'SEGF_',
            'ion_prefix': 'SIGF_',
            'power_prefix': 'PIGF_',
            'torque_prefix': 'TQIGF_',
            'cx_particles_prefix': 'SFCXGF_',
            'cx_power_prefix': 'PCXGF_',
            'cx_torque_prefix': 'TQCXGF_',
            'aggregate_electron': 'SESGF',
            'aggregate_ion_prefix': 'SISGF_',
            'aggregate_power': 'PISGF',
            'aggregate_torque': 'TQISGF',
            'aggregate_cx_particles': 'SFCXSGF',
            'aggregate_cx_power': 'PCXSGF',
            'aggregate_cx_torque': 'TQCXSGF',
            'description': 'Externally supplied gas-flow neutrals',
        },
        'recycling': {
            'electron_prefix': 'SERC_',
            'ion_prefix': 'SIRC_',
            'power_prefix': 'PIRC_',
            'torque_prefix': 'TQIRC_',
            'cx_particles_prefix': 'SFCXRC_',
            'cx_power_prefix': 'PCXRC_',
            'cx_torque_prefix': 'TQCXRC_',
            'aggregate_electron': 'SESRC',
            'aggregate_ion_prefix': 'SISRC_',
            'aggregate_power': 'PISRC',
            'aggregate_torque': 'TQISRC',
            'aggregate_cx_particles': 'SFCXSRC',
            'aggregate_cx_power': 'PCXSRC',
            'aggregate_cx_torque': 'TQCXSRC',
            'description': 'Recycling neutrals',
        },
    }

    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.time = np.asarray(self.im.get_ncvar('TIME3'), dtype=np.float64)
        self.rho = np.asarray(self.im.get_ncvar('X'), dtype=np.float64)
        self.nt = len(self.time)
        self.nx = self.rho.shape[1] if self.rho.ndim == 2 else 0
        self.create_ids()

    def _has(self, name):
        return bool(self.im.check_var_name(name))

    def _get(self, name, mult=1.0):
        if not self._has(name):
            return np.array([])
        arr = self.im.get_ncvar(name)
        if not self.im.is_valid(arr):
            return np.array([])
        arr = np.asarray(arr, dtype=np.float64)
        if arr.shape != self.rho.shape:
            return np.array([])
        return arr * float(mult)

    @staticmethod
    def _sum_valid(arrays):
        arrays = [np.asarray(a, dtype=np.float64) for a in arrays
                  if isinstance(a, np.ndarray) and a.size]
        if not arrays:
            return np.array([])
        out = np.zeros_like(arrays[0], dtype=np.float64)
        for arr in arrays:
            if arr.shape == out.shape:
                out += arr
        return out

    @staticmethod
    def _canonical_token(token):
        token = str(token).strip().upper()
        return {'3': 'HE3', '4': 'HE4'}.get(token, token)

    def _atom(self, token):
        key = self._canonical_token(token)
        try:
            atom = atom_base(key)
            if not atom.label or atom.label == 'dummy':
                return key, None
            return key, atom
        except KeyError:
            self.im.mapping.skipped(f'unknown species token {key!r}', target='source/ggd/ion')
            return key, None

    def _fill_species_identity(self, node, token, ion=False):
        key, atom = self._atom(token)
        if atom is None:
            self.im.set_species_name(node, key)
            if ion and hasattr(node, 'z_ion'):
                node.z_ion = 1.0
            return False
        node.element.resize(1)
        node.element[0].a = float(atom.a)
        self.im.set_nuclear_charge(node.element[0], atom.z_n)
        node.element[0].atoms_n = int(atom.atoms_n)
        self.im.set_species_name(node, atom.label.strip())
        if ion and hasattr(node, 'z_ion'):
            # Gas/recycling ionization creates a singly charged ion at birth.
            node.z_ion = 1.0
        return True

    def _set_source_species_neutral(self, source, token):
        source.species.type.name = 'neutral'
        source.species.type.index = 4
        source.species.type.description = 'Neutral species'
        self._fill_species_identity(source.species.neutral, token, ion=False)

    def _set_source_species_ion(self, source, token):
        source.species.type.name = 'ion'
        source.species.type.index = 2
        source.species.type.description = 'Ion species in a single/average state'
        self._fill_species_identity(source.species.ion, token, ion=True)

    @staticmethod
    def _set_identifier(identifier, name, index, description):
        identifier.name = name
        identifier.index = int(index)
        identifier.description = description

    def _set_grid_coordinate_type(self, space):
        """Fill coordinates_type in both DD3 (INT_1D) and DD4 (identifier AoS)."""
        if self.im.dd_at_least(4, 0, 0):
            space.coordinates_type.resize(1)
            self._set_identifier(
                space.coordinates_type[0], 'rho_tor_norm', 12,
                'Normalized toroidal flux coordinate')
        else:
            space.coordinates_type.resize(1)
            space.coordinates_type[0] = 12

    def _fill_grid(self, grid, it):
        grid.time = float(self.time[it])
        self._set_identifier(grid.identifier, 'linear', self.GRID_INDEX,
                             'One-dimensional TRANSP rho_tor_norm grid')

        grid.space.resize(1)
        space = grid.space[0]
        self._set_identifier(space.identifier, 'primary_standard', 1,
                             'Primary space defining the standard grid')
        self._set_identifier(space.geometry_type, 'standard', 0,
                             'Standard geometry')
        self._set_grid_coordinate_type(space)

        # GGD dimension convention: objects_per_dimension(1) are nodes.
        space.objects_per_dimension.resize(1)
        objects = space.objects_per_dimension[0]
        self._set_identifier(objects.geometry_content, 'node_coordinates', 1,
                             'Geometry contains the coordinates of each node')
        objects.object.resize(self.nx)
        for ix in range(self.nx):
            objects.object[ix].geometry = np.asarray([self.rho[it, ix]], dtype=np.float64)

        # Identifier index 1 is the implicit subset containing all nodes.
        grid.grid_subset.resize(1)
        subset = grid.grid_subset[0]
        self._set_identifier(subset.identifier, 'nodes', self.GRID_SUBSET_INDEX,
                             'All nodes of the associated space')
        subset.dimension = 1

    def _set_scalar_ggd(self, aos, values, it, sign=1.0):
        if not isinstance(values, np.ndarray) or not values.size:
            return False
        aos.resize(1)
        q = aos[0]
        q.grid_index = self.GRID_INDEX
        q.grid_subset_index = self.GRID_SUBSET_INDEX
        data = sign * np.asarray(values[it], dtype=np.float64)
        q.values.resize(len(data))
        q.values[:] = data
        return True

    def _set_momentum_ggd(self, aos, values, it, sign=1.0):
        if not isinstance(values, np.ndarray) or not values.size:
            return False
        field = 'phi' if self.im.dd_at_least(3, 42, 0) else 'toroidal'
        aos.resize(1)
        q = aos[0]
        if any(str(c) for c in getattr(q, field).metadata.coordinates_same_as):
            # DD 3.39 ties the toroidal component's size to the radial one; filling radial with zeros
            # would assert a radial momentum source that TRANSP does not provide.
            aos.resize(0)
            self.im.mapping.skipped(f'{field} momentum not written: the DD ties it to the radial component',
                                    target='source/ggd/*/momentum')
            return False
        q.grid_index = self.GRID_INDEX
        q.grid_subset_index = self.GRID_SUBSET_INDEX
        data = sign * np.asarray(values[it], dtype=np.float64)
        arr = getattr(q, field)
        arr.resize(len(data))
        arr[:] = data
        return True

    def _neutral_tokens(self, spec):
        prefix = spec['electron_prefix']
        tokens = []
        for name in getattr(self.im, 'vars_', []):
            if name.startswith(prefix):
                token = name[len(prefix):]
                if token and '_' not in token and self._atom(token)[1] is not None:
                    tokens.append(token)
        if tokens:
            return sorted(set(tokens))

        # Fallback to aggregate SISGF_G / SISRC_G outputs when old CDFs do not
        # contain the neutral-resolved electron aliases.
        prefix = spec['aggregate_ion_prefix']
        for name in getattr(self.im, 'vars_', []):
            if name.startswith(prefix):
                token = name[len(prefix):]
                if token and '_' not in token and self._atom(token)[1] is not None:
                    tokens.append(token)
        return sorted(set(tokens))

    def _ion_outputs(self, spec, neutral_token):
        """Return [(ion_token, source_profile)] for one neutral-origin source."""
        prefix = spec['ion_prefix']
        pat = re.compile(rf'^{re.escape(prefix)}([A-Z0-9]+)_{re.escape(neutral_token)}$')
        outputs = []
        for name in getattr(self.im, 'vars_', []):
            match = pat.fullmatch(name)
            if not match:
                continue
            ion_token = match.group(1)
            if self._atom(ion_token)[1] is None:
                continue
            arr = self._get(name, self.CM3_TO_M3)
            if arr.size:
                outputs.append((ion_token, arr, name))

        if outputs:
            return outputs

        # Older output can have only the group-total SISGF_G / SISRC_G array.
        name = f"{spec['aggregate_ion_prefix']}{neutral_token}"
        arr = self._get(name, self.CM3_TO_M3)
        if arr.size:
            return [(neutral_token, arr, name)]
        return []

    def _origin_inventory(self, origin):
        spec = self.ORIGINS[origin]
        inventory = []
        for neutral in self._neutral_tokens(spec):
            electron_name = f"{spec['electron_prefix']}{neutral}"
            electron = self._get(electron_name, self.CM3_TO_M3)
            if not electron.size and len(self._neutral_tokens(spec)) == 1:
                electron_name = spec['aggregate_electron']
                electron = self._get(electron_name, self.CM3_TO_M3)

            ions = self._ion_outputs(spec, neutral)
            ion_sum = self._sum_valid([a for _sp, a, _name in ions])
            # For an atomic singly-ionized neutral, electron and ion source rates
            # are equal.  Prefer the species-resolved ion sum for the neutral sink.
            neutral_loss = ion_sum if ion_sum.size else electron

            power_name = f"{spec['power_prefix']}{neutral}"
            power = self._get(power_name, self.CM3_TO_M3)
            if not power.size and len(self._neutral_tokens(spec)) == 1:
                power_name = spec['aggregate_power']
                power = self._get(power_name, self.CM3_TO_M3)

            torque_name = f"{spec['torque_prefix']}{neutral}"
            torque = self._get(torque_name, self.CM3_TO_M3)
            if not torque.size and len(self._neutral_tokens(spec)) == 1:
                torque_name = spec['aggregate_torque']
                torque = self._get(torque_name, self.CM3_TO_M3)

            cx_particles_name = f"{spec['cx_particles_prefix']}{neutral}"
            cx_particles = self._get(cx_particles_name, self.CM3_TO_M3)
            if not cx_particles.size and len(self._neutral_tokens(spec)) == 1:
                cx_particles_name = spec['aggregate_cx_particles']
                cx_particles = self._get(cx_particles_name, self.CM3_TO_M3)

            cx_power_name = f"{spec['cx_power_prefix']}{neutral}"
            cx_power = self._get(cx_power_name, self.CM3_TO_M3)
            if not cx_power.size and len(self._neutral_tokens(spec)) == 1:
                cx_power_name = spec['aggregate_cx_power']
                cx_power = self._get(cx_power_name, self.CM3_TO_M3)

            cx_torque_name = f"{spec['cx_torque_prefix']}{neutral}"
            cx_torque = self._get(cx_torque_name, self.CM3_TO_M3)
            if not cx_torque.size and len(self._neutral_tokens(spec)) == 1:
                cx_torque_name = spec['aggregate_cx_torque']
                cx_torque = self._get(cx_torque_name, self.CM3_TO_M3)

            if (electron.size or ions or power.size or torque.size or
                    cx_particles.size or cx_power.size or cx_torque.size):
                inventory.append({
                    'origin': origin,
                    'neutral': neutral,
                    'electron': electron,
                    'electron_name': electron_name,
                    'ions': ions,
                    'neutral_loss': neutral_loss,
                    'power': power,
                    'power_name': power_name,
                    'torque': torque,
                    'torque_name': torque_name,
                    'cx_particles': cx_particles,
                    'cx_particles_name': cx_particles_name,
                    'cx_power': cx_power,
                    'cx_power_name': cx_power_name,
                    'cx_torque': cx_torque,
                    'cx_torque_name': cx_torque_name,
                })
        return inventory

    def _fill_ionization_source(self, source, item):
        """Write the electron-impact ionization part of a neutral-origin source.

        TRANSP ``SIGF_*_*`` / ``SIRC_*_*`` ion sources include the ion born by
        fast-ion charge exchange as well as electron-impact ionization.  The
        electron source ``SEGF_*`` / ``SERC_*`` isolates atomic ionization and
        is therefore used for the particle stoichiometry here.  The CX piece is
        written separately by :meth:`_fill_charge_exchange_source`.
        """
        origin = item['origin']
        spec = self.ORIGINS[origin]
        neutral_token = item['neutral']
        if not (item['electron'].size or item['power'].size or item['torque'].size):
            return False

        self._set_identifier(
            source.identifier, 'atomic_ionization', 707,
            f"Atomic ionization of {spec['description'].lower()} ({neutral_token})")
        self._set_source_species_neutral(source, neutral_token)
        source.ggd.resize(self.nt)

        # Atomic ionization of one neutral produces one singly-charged ion and
        # one electron.  Use SE*, not SI*, because SI*=SE*+SFCX* in TRANSP.
        particles = item['electron']
        for it in range(self.nt):
            ggd = source.ggd[it]
            ggd.time = float(self.time[it])

            self._set_scalar_ggd(ggd.electrons.particles, particles, it)
            ggd.ion.resize(1)
            ion = ggd.ion[0]
            self._fill_species_identity(ion, neutral_token, ion=True)
            self._set_scalar_ggd(ion.particles, particles, it)

            ggd.neutral.resize(1)
            neutral = ggd.neutral[0]
            self._fill_species_identity(neutral, neutral_token, ion=False)
            self._set_scalar_ggd(neutral.particles, particles, it, sign=-1.0)

            # PIGF/PIRC are the kinetic-energy transfer associated with neutral
            # ionization.  TQI* is the corresponding toroidal torque.
            if item['power'].size:
                self._set_scalar_ggd(ion.energy, item['power'], it)
                self._set_scalar_ggd(neutral.energy, item['power'], it, sign=-1.0)
            if item['torque'].size:
                self._set_momentum_ggd(ion.momentum, item['torque'], it)
                self._set_momentum_ggd(neutral.momentum, item['torque'], it, sign=-1.0)

        self.logger.info(
            f"edge_sources: {origin}/{neutral_token} atomic_ionization: "
            f"particles={item['electron_name']}, power={item['power_name']}, "
            f"torque={item['torque_name']}")
        return True

    def _fill_charge_exchange_source(self, source, item):
        """Write the fast-ion charge-exchange contribution separately.

        ``SFCX*`` is the neutral-to-ion conversion count associated with
        fast-ion CX.  ``PCX*`` and ``TQCX*`` are explicitly power/torque *to
        neutrals*, so the plasma-ion terms carry the opposite sign.  The CDF
        does not resolve the donor fast-ion species for the energy/momentum
        transfer; consequently those terms are written at total-ion level.
        """
        origin = item['origin']
        spec = self.ORIGINS[origin]
        neutral_token = item['neutral']
        particles = item['cx_particles']
        power = item['cx_power']
        torque = item['cx_torque']
        if not (particles.size or power.size or torque.size):
            return False

        self._set_identifier(
            source.identifier, 'charge_exchange', 305,
            f"Charge exchange involving {spec['description'].lower()} ({neutral_token}); "
            "TRANSP SFCX/PCX/TQCX terms")
        self._set_source_species_neutral(source, neutral_token)
        source.ggd.resize(self.nt)

        for it in range(self.nt):
            ggd = source.ggd[it]
            ggd.time = float(self.time[it])

            # SFCX is the ion born from the target neutral in a fast-ion CX
            # event.  The corresponding donor-fast-ion sink is not species-
            # resolved by these CDF arrays, so it cannot be represented as a
            # separate ion-particle sink here.
            ggd.ion.resize(1)
            ion = ggd.ion[0]
            self._fill_species_identity(ion, neutral_token, ion=True)
            self._set_scalar_ggd(ion.particles, particles, it)

            ggd.neutral.resize(1)
            neutral = ggd.neutral[0]
            self._fill_species_identity(neutral, neutral_token, ion=False)
            self._set_scalar_ggd(neutral.particles, particles, it, sign=-1.0)

            # PCX/TQCX are positive transfers to neutrals.  Use aggregate ion
            # nodes for the opposite plasma sink because the donor ion is not
            # resolved in these neutral-transport diagnostics.
            if power.size:
                self._set_scalar_ggd(ggd.total_ion_energy, power, it, sign=-1.0)
                self._set_scalar_ggd(neutral.energy, power, it)
            if torque.size:
                self._set_momentum_ggd(ggd.momentum, torque, it, sign=-1.0)
                self._set_momentum_ggd(neutral.momentum, torque, it)

        self.logger.info(
            f"edge_sources: {origin}/{neutral_token} charge_exchange: "
            f"particles={item['cx_particles_name']}, power={item['cx_power_name']}, "
            f"torque={item['cx_torque_name']}")
        return True

    def _recombination_tokens(self):
        """Return species resolved by S0V0G/SIV0G recombination diagnostics."""
        tokens = set()
        for name in getattr(self.im, 'vars_', []):
            match = re.fullmatch(r'S0V0([A-Z0-9]+)', name)
            if not match:
                match = re.fullmatch(r'SIV0([A-Z0-9]+)', name)
            if match and self._atom(match.group(1))[1] is not None:
                tokens.add(match.group(1))
        return sorted(tokens)

    def _fill_recombination_source(self, source):
        particles = self._get('S0RECO', self.CM3_TO_M3)
        power = self._get('P0RECO', self.CM3_TO_M3)
        torque = self._get('TQ0RECO', self.CM3_TO_M3)
        if not (particles.size or power.size or torque.size):
            return False

        tokens = self._recombination_tokens()
        resolved = len(tokens) == 1
        token = tokens[0] if resolved else None
        if resolved:
            species_particles = self._get(f'S0V0{token}', self.CM3_TO_M3)
            if species_particles.size:
                particles = species_particles
            self._set_source_species_ion(source, token)
            description = f'Plasma recombination producing neutral {token}'
        else:
            source.species.type.name = 'unspecified'
            source.species.type.index = 0
            source.species.type.description = 'Recombining species is not uniquely resolved'
            description = 'Neutral source produced by plasma recombination'

        self._set_identifier(source.identifier, 'recombination', 710, description)
        source.ggd.resize(self.nt)

        for it in range(self.nt):
            ggd = source.ggd[it]
            ggd.time = float(self.time[it])
            if resolved:
                # Atomic recombination consumes one ion and one electron and
                # produces one neutral of the same isotope.  Only particle
                # stoichiometry is inferred; P0RECO/TQ0RECO give the kinetic
                # content of the born neutral but do not resolve the plasma
                # energy/momentum partition of the corresponding sinks.
                self._set_scalar_ggd(ggd.electrons.particles, particles, it, sign=-1.0)
                ggd.ion.resize(1)
                ion = ggd.ion[0]
                self._fill_species_identity(ion, token, ion=True)
                self._set_scalar_ggd(ion.particles, particles, it, sign=-1.0)
                ggd.neutral.resize(1)
                neutral = ggd.neutral[0]
                self._fill_species_identity(neutral, token, ion=False)
            else:
                # Retain the total neutral production without inventing a
                # species split when TRANSP did not save it.
                ggd.neutral.resize(1)
                neutral = ggd.neutral[0]
                self.im.set_species_name(neutral, 'recombination_neutral_total')

            self._set_scalar_ggd(neutral.particles, particles, it)
            self._set_scalar_ggd(neutral.energy, power, it)
            self._set_momentum_ggd(neutral.momentum, torque, it)

        detail = f' ({token})' if resolved else ' (species unresolved)'
        self.logger.info(
            'edge_sources: recombination neutral <- S0RECO/P0RECO/TQ0RECO' + detail)
        return True

    def _fill_recombination_recapture_ionization(self, source):
        tokens = self._recombination_tokens()
        if len(tokens) != 1:
            return False
        token = tokens[0]
        particles = self._get('SERECO', self.CM3_TO_M3)
        power = self._get('PIRECO', self.CM3_TO_M3)
        torque = self._get('TQIRECO', self.CM3_TO_M3)
        if not (particles.size or power.size or torque.size):
            return False

        self._set_identifier(
            source.identifier, 'atomic_ionization', 707,
            f'Re-ionization/recapture of recombination neutrals ({token})')
        self._set_source_species_neutral(source, token)
        source.ggd.resize(self.nt)
        for it in range(self.nt):
            ggd = source.ggd[it]
            ggd.time = float(self.time[it])
            self._set_scalar_ggd(ggd.electrons.particles, particles, it)
            ggd.ion.resize(1)
            ion = ggd.ion[0]
            self._fill_species_identity(ion, token, ion=True)
            self._set_scalar_ggd(ion.particles, particles, it)
            ggd.neutral.resize(1)
            neutral = ggd.neutral[0]
            self._fill_species_identity(neutral, token, ion=False)
            self._set_scalar_ggd(neutral.particles, particles, it, sign=-1.0)
            if power.size:
                self._set_scalar_ggd(ion.energy, power, it)
                self._set_scalar_ggd(neutral.energy, power, it, sign=-1.0)
            if torque.size:
                self._set_momentum_ggd(ion.momentum, torque, it)
                self._set_momentum_ggd(neutral.momentum, torque, it, sign=-1.0)
        self.logger.info(
            f'edge_sources: recombination-neutral recapture ionization ({token}) '
            '<- SERECO/PIRECO/TQIRECO')
        return True

    def _fill_recombination_charge_exchange(self, source):
        tokens = self._recombination_tokens()
        if len(tokens) != 1:
            return False
        token = tokens[0]
        particles = self._get('SFCXRECO', self.CM3_TO_M3)
        power = self._get('PCXRECO', self.CM3_TO_M3)
        torque = self._get('TQCXRECO', self.CM3_TO_M3)
        if not (particles.size or power.size or torque.size):
            return False

        self._set_identifier(
            source.identifier, 'charge_exchange', 305,
            f'Charge exchange involving recombination neutrals ({token})')
        self._set_source_species_neutral(source, token)
        source.ggd.resize(self.nt)
        for it in range(self.nt):
            ggd = source.ggd[it]
            ggd.time = float(self.time[it])
            ggd.ion.resize(1)
            ion = ggd.ion[0]
            self._fill_species_identity(ion, token, ion=True)
            self._set_scalar_ggd(ion.particles, particles, it)
            ggd.neutral.resize(1)
            neutral = ggd.neutral[0]
            self._fill_species_identity(neutral, token, ion=False)
            self._set_scalar_ggd(neutral.particles, particles, it, sign=-1.0)
            if power.size:
                self._set_scalar_ggd(ggd.total_ion_energy, power, it, sign=-1.0)
                self._set_scalar_ggd(neutral.energy, power, it)
            if torque.size:
                self._set_momentum_ggd(ggd.momentum, torque, it, sign=-1.0)
                self._set_momentum_ggd(neutral.momentum, torque, it)
        self.logger.info(
            f'edge_sources: recombination-neutral charge exchange ({token}) '
            '<- SFCXRECO/PCXRECO/TQCXRECO')
        return True

    def create_ids(self):
        self.logger.info('\nedge_sources')
        if self.rho.ndim != 2 or self.nt == 0 or self.nx == 0:
            self.logger.info('edge_sources: missing TIME3/X radial grid; IDS not written')
            return

        inventory = []
        for origin in ('gas_flow', 'recycling'):
            inventory.extend(self._origin_inventory(origin))
        has_reco = self._has('S0RECO') or self._has('P0RECO') or self._has('TQ0RECO')
        has_reco_recapture = any(self._has(name) for name in
                                 ('SERECO', 'PIRECO', 'TQIRECO',
                                  'SFCXRECO', 'PCXRECO', 'TQCXRECO'))
        if not inventory and not has_reco and not has_reco_recapture:
            self.logger.info('edge_sources: no supported gas/recycling/recombination source profiles found')
            return

        es = self.im.new_ids('edge_sources')
        es = self.im.fill_properties(es)
        self.im.fill_code_properties(es)
        es.ids_properties.homogeneous_time = 1
        es.time = np.asarray(self.time, dtype=np.float64)

        es.midplane.name = 'magnetic_axis'
        es.midplane.index = 1
        es.midplane.description = (
            'Midplane defined by the magnetic-axis height from equilibrium/'
            'time_slice/global_quantities/magnetic_axis/z')

        es.grid_ggd.resize(self.nt)
        for it in range(self.nt):
            self._fill_grid(es.grid_ggd[it], it)

        # Each gas/recycling neutral population can generate two distinct
        # physical processes: atomic ionization and charge exchange.
        builders = []
        for item in inventory:
            builders.append(lambda source, item=item: self._fill_ionization_source(source, item))
            builders.append(lambda source, item=item: self._fill_charge_exchange_source(source, item))
        if has_reco:
            builders.append(self._fill_recombination_source)
        if has_reco_recapture:
            builders.append(self._fill_recombination_recapture_ionization)
            builders.append(self._fill_recombination_charge_exchange)
        used = append_successful_entries(lambda: self.im.new_ids('edge_sources'), es, 'source', builders)

        if used:
            self.im.putdb_imas(es)
            self.logger.info(f'edge_sources: wrote {used} spatial source entries')
        else:
            self.logger.info('edge_sources: no populated source entries; IDS not written')


if __name__ == '__main__':
    print('Please use: python convert.py')
