#!/usr/bin/python
"""Populate the IMAS nbi IDS from TRANSP/NUBEAM injector configuration and time traces."""

from nc2imas import *
from writer_base import IdsWriter


class imas_nbi(IdsWriter):
    """Map the TRANSP injector time traces to the nbi IDS, one unit per ``PINJ<n>`` source.

    The native ``TIME`` axis is written unchanged for every unit.  A source that is off is expressed by
    its power (``PINJ<n>`` is zero), never by altering the clock: an earlier version replaced the times by
    zero wherever any of the three energy components was zero, which wrote all-zero time axes for every
    beam that has only a full-energy component (both ITER beams, all 14 DIII-D 163303U09 beams).
    """

    def _sources(self):
        """``[(source suffix, PINJ name)]`` of the ``PINJ<n>`` channels, in numeric order, keeping their own numbers."""
        names = self.cdf.get_var_names_re('PINJ?', nd1=1, nd2=2)
        return sorted(((name[len('PINJ'):], name) for name in names), key=lambda item: int(item[0]))

    def _species(self, species_suffixes, unit_count):
        """Species for each unit, or ``None`` where the data does not say.

        One ``PINJ_<sp>`` species-total channel means every unit injects it.  With several species TRANSP gives
        no channel that links a unit to one of them, and the order of the names is not evidence, so none is assigned.
        """
        if len(species_suffixes) == 1:
            return [species_suffixes[0]] * unit_count
        if len(species_suffixes) > 1 and unit_count:
            self.ctx.mapping.skipped(
                f'beam species {species_suffixes} cannot be matched to the {unit_count} PINJ<n> units: '
                'no channel links a unit to a species, and list position is not evidence',
                target='unit/species')
        return [None] * unit_count

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
        node.data.resize(len(data))
        node.data = data
        node.time.resize(len(times))
        node.time = times

    def build(self):
        self.logger.info('\nnbi')
        nbi = self.dd.new_ids('nbi')
        timesc = self.cdf.get_ncvar('TIME')                 # native scalar time, written as it is
        self.fill.checkErrExit(timesc)

        sources = self._sources()
        species = self.cdf.get_var_names('PINJ_?', cut='PINJ_')
        self.logger.info(f'NBeams: PINJ {[name for _s, name in sources]}; species: {species}')

        if sources:
            units = [(f'unit_{int(suffix)}', f'TRANSP NBI unit {int(suffix)}', f'EINJ{suffix}_E1', name)
                     for suffix, name in sources]
            symbols = self._species(species, len(units))
        else:
            # no per-source channels: one unit per species total, named by its own suffix
            units = [(f'species_{sp}', f'TRANSP NBI species {sp}', 'EINJ', f'PINJ_{sp}') for sp in species]
            symbols = list(species)

        nbi.unit.resize(len(units))
        for unit, (identifier, description, energy_name, power_name), symbol in zip(nbi.unit, units, symbols):
            self.fill.set_nbi_unit_identity(unit, identifier, description)
            if symbol is not None:
                self._set_species(unit, symbol)
            # energy EINJ<n>_E1 (full-energy component) or EINJ (max initial beam energy); power PINJ<n> or PINJ_<sp>
            self._fill_trace(unit.energy, energy_name, self.cdf.get_ncvar(energy_name), timesc)
            self._fill_trace(unit.power_launched, power_name, self.cdf.get_ncvar(power_name), timesc)

        if units:
            nbi = self.fill.fill_properties(nbi)
            self.fill.fill_code_properties(nbi)
            self.fill.fill_external_model_inputs(nbi, 'nbi')
            nbi.ids_properties.homogeneous_time = 1
            nbi.time.resize(len(timesc))
            nbi.time = timesc
            self.logger.info('Done NBI')
            return [nbi]
        return []

if __name__ == '__main__':
    print("Please use: python convert.py ")
