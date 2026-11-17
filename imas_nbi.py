#!/usr/bin/python
"""Populate the IMAS nbi IDS from TRANSP/NUBEAM injector configuration and time traces."""

from nc2imas import *
from writer_base import IdsWriter


# Relative closure, over the species total, within which the injector powers are taken to add up to it.  The bundled
# cases close to 2.4e-7; the tolerance only has to separate a sum from a different quantity.
SPECIES_CLOSURE_RTOL = 1.0e-4


class imas_nbi(IdsWriter):
    """Map the TRANSP injector time traces to the nbi IDS, one unit per ``PINJ<n>`` source.

    Inventory.  Each unit follows *original source number -> its own CDF channels -> species evidence -> output
    unit*.  The number is the ``PINJ<n>`` suffix, so a list position never stands for an injector; the channels
    ``EINJ<n>_E1`` and ``PINJ<n>`` are named from that suffix.

    Time.  The native ``TIME`` axis is written unchanged for every unit and must be finite and strictly
    increasing; anything else is reported, never sorted, de-duplicated or invented.  A source that is off is
    expressed by its power (``PINJ<n>`` is zero), never by altering the clock.  An earlier version replaced the
    times by zero wherever any of the three energy components was zero, which wrote all-zero time axes for every
    beam that has only a full-energy component (both ITER beams, all 14 DIII-D 163303U09 beams).

    Species.  TRANSP has no channel that links a unit to a species.  The one supported piece of evidence is power
    closure: when exactly one species total ``PINJ_<sp>`` exists and the injector powers add up to it
    (to ``SPECIES_CLOSURE_RTOL``), every injector injects that species.  Several species, or a total the injectors
    do not add up to, are ambiguous: the units are written without a species and the ambiguity is recorded.  A
    species total alone, with no ``PINJ<n>`` channel, is an aggregate over unknown injectors and is not written as
    a unit.
    """

    def _sources(self):
        """``[(source suffix, PINJ name)]`` of the ``PINJ<n>`` channels, in numeric order, keeping their own numbers."""
        names = self.cdf.get_var_names_re('PINJ?', nd1=1, nd2=2)
        return sorted(((name[len('PINJ'):], name) for name in names), key=lambda item: int(item[0]))

    def _native_time(self):
        """The ``TIME`` axis, or ``None`` after recording why it cannot be used (malformed times are not repaired)."""
        times = self.cdf.get_ncvar('TIME')
        if times.ndim != 1 or times.size == 0:
            self.ctx.mapping.broken('TIME is absent or not a one-dimensional axis', source='TIME', target='unit/time')
            return None
        if not np.all(np.isfinite(times)):
            self.ctx.mapping.broken('TIME has non-finite samples', source='TIME', target='unit/time')
            return None
        if times.size > 1 and not np.all(np.diff(times) > 0):
            self.ctx.mapping.broken('TIME is not strictly increasing (duplicate or reordered samples); it is not '
                                    'sorted or de-duplicated', source='TIME', target='unit/time')
            return None
        return times

    def _species(self, species_suffixes, units):
        """Species for each unit, or ``None`` where the data does not establish it (the reason is recorded).

        ``units`` is ``[(source suffix, power array)]``.
        """
        count = len(units)
        if not species_suffixes or not count:
            return [None] * count
        if len(species_suffixes) > 1:
            self.ctx.mapping.skipped(
                f'beam species {species_suffixes} cannot be matched to the {count} PINJ<n> units: no channel links '
                'a unit to a species, and list position is not evidence', target='unit/species')
            return [None] * count
        species = species_suffixes[0]
        total = self.cdf.get_ncvar(f'PINJ_{species}')
        powers = [power for _suffix, power in units]
        if self.fill.is_valid(total) and all(len(power) == len(total) for power in powers):
            reference = float(np.max(np.abs(total)))
            gap = float(np.max(np.abs(np.sum(powers, axis=0) - total)))
            if reference > 0.0 and gap <= SPECIES_CLOSURE_RTOL * reference:
                return [species] * count
            detail = (f'relative gap {gap / reference:.3g} against {SPECIES_CLOSURE_RTOL:g}' if reference > 0.0
                      else 'the species total is zero')
        else:
            detail = 'the traces are absent or differ in length'
        self.ctx.mapping.skipped(
            f'the injector powers PINJ<n> do not demonstrably add up to the species total PINJ_{species} ({detail}), '
            'so the species of the units is ambiguous', target='unit/species')
        return [None] * count

    def _set_species(self, unit, symbol):
        atoms = atom_base()
        sp_zn, sp_a, sp_label = atoms.find_beam_species_data(symbol)
        if self.fill.is_valid(sp_a, sp_zn):
            unit.species.a = float(sp_a)
            # z_n is INT_0D in DD4 (it was FLT_0D in older DDs).
            self.fill.set_nuclear_charge(unit.species, sp_zn)
            # 'label' was renamed to 'name' in DD 3.42.  Strip the TRANSP bookkeeping suffix from names such as D_BEAM.
            self.fill.set_species_name(unit.species, str(sp_label).replace('_BEAM', ''))

    def _fill_trace(self, node, name, data, times):
        """Assign ``data`` and the native ``times``; a trace of another length is a broken mapping."""
        if not self.fill.is_valid(data):
            return
        if len(data) != len(times):
            self.ctx.mapping.broken(f'{len(data)} samples on a {len(times)}-sample TIME axis', source=name,
                                    target='unit/' + name)
            return
        node.data = data
        node.time = times

    def build(self):
        self.logger.info('\nnbi')
        sources = self._sources()
        species = self.cdf.get_var_names('PINJ_?', cut='PINJ_')
        self.logger.info(f'NBeams: PINJ {[name for _s, name in sources]}; species: {species}')

        if not sources:
            if species:
                # an aggregate over injectors that are not identified is not an injector
                self.ctx.mapping.skipped(
                    f'species totals {["PINJ_" + sp for sp in species]} without per-source PINJ<n> channels are '
                    'aggregates over unknown injectors; no nbi unit is written for them', target='unit')
            return []
        times = self._native_time()
        if times is None:
            return []

        nbi = self.dd.new_ids('nbi')
        powers = [(suffix, self.cdf.get_ncvar(name)) for suffix, name in sources]
        symbols = self._species(species, powers)
        nbi.unit.resize(len(sources))
        for unit, (suffix, name), (_s, power), symbol in zip(nbi.unit, sources, powers, symbols):
            self.fill.set_nbi_unit_identity(unit, f'unit_{int(suffix)}', f'TRANSP NBI unit {int(suffix)}')
            if symbol is not None:
                self._set_species(unit, symbol)
            # energy: the full-energy component EINJ<n>_E1; power PINJ<n> (the source's own channels, by suffix)
            energy_name = f'EINJ{suffix}_E1'
            self._fill_trace(unit.energy, energy_name, self.cdf.get_ncvar(energy_name), times)
            self._fill_trace(unit.power_launched, name, power, times)

        nbi = self.fill.fill_properties(nbi)
        self.fill.fill_code_properties(nbi)
        self.fill.fill_external_model_inputs(nbi, 'nbi')
        nbi.ids_properties.homogeneous_time = 1
        nbi.time = times
        self.logger.info('Done NBI')
        return [nbi]

if __name__ == '__main__':
    print("Please use: python convert.py ")
