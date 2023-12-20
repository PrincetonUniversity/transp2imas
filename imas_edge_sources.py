

from nc2imas import *

class imas_edge_sources:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger

        self.create_ids()

    def fill_source_data(self, es, i, name, times):
        if name == 'ionization':
            dstr = 'N0SGF'
        elif name == 'recombination'    :
            dstr = 'N0SRC'
        species = self.im.get_var_names(f'{dstr}_?', cut=f'{dstr}_')
        self.logger.info(f'{name}: {species}')

        #ggd neutrals
        es.source[i].ggd.resize(len(times))
        print('\n check ggd size\n\n')
        for it in range(len(times)):
            es.source[i].ggd[it].time = times[it].astype(float)
            es.source[i].ggd[it].neutral.resize(len(species))
            for ni, short in enumerate(species):
               key = 'HE' + short if short in ['3', '4'] else short
               atom = atom_base(key)
               es.source[i].ggd[it].neutral[ni].element.resize(1)
               es.source[i].ggd[it].neutral[ni].element[0].a = atom.a
               es.source[i].ggd[it].neutral[ni].element[0].z_n = atom.z_n
               es.source[i].ggd[it].neutral[ni].element[0].atoms_n = atom.atoms_n

               es.source[i].ggd[it].neutral[ni].label = key

               es.source[i].ggd[it].neutral[ni].particles.resize(1)
               self.im.fill_data_nc(es.source[i].ggd[it].neutral[ni].particles[0].values,
                        f'{dstr}_{short}', ind=it, mult=1.0e6)
        return es


    def create_ids(self):
        self.im.logger.info('\nedge_sources')
        es = self.im.ids.edge_sources

        es = self.im.fill_properties(es)
        self.im.fill_code_properties(es)
        es.ids_properties.homogeneous_time = 1
        es.time.resize(self.im.ntimep)
        es.time = self.im.timep
        ntime = self.im.ntimep

        es.midplane.name = 'magnetic_axis'
        es.midplane.index = 1
        txt = r'Midplane defined by the height of magnetic axis'
        es.midplane.description = txt + r' equilibrium/time_slice/global_quantities/magnetic_axis/z'

        #GASIZ(TIME)ATOMS/SEC; H+D+T EDGE NTRL SCE IONIZATIONS (heatlib/source)

        #E0INR(TIME) EV ; T0 (RECYCLING) @EDGE

        #DN0OUT: NEUT. DENSITY OUTSIDE THE PLASMA
        #GASIN ; ATOMS/SEC ; H+D+T GAS SOURCE AT BDY. RECYCD+GASFLD (Recycl+GASFLOW RATE DATA)

        #T0EDGE
        #T0RECYC(TIME) EV; RECYCLING GAS TEMPERATURE (AVG)
        #T0GASFL(TIME) EV; GAS FLOW TEMPERATURE (AVG)
        #OM0EDGE
        #OM0RECYC(TIME) ; RAD/SEC RECYCLING GAS Ang. Veloc. (AVG)
        #OM0GASFL(TIME) ; RAD/SEC GAS FLOW Ang. Veloc. (AVG)
        #PEDGE EXPANSION/SCRAPEOFF POWER
        #PEEDG(TIME)  WATTS; ELECTRON ENERGY VIA BDY
        #PIEDG(TIME)  WATTS; ION ENERGY VIA BDY

        #PICHEXT(TIME)
        #search EXT and EDGE

        #FUSN (P14) POWER TO CX (EXT)

        #NOTE: edge profile with density needed in IDS. The above is only integrated

        #grid_ggd(itime)

        #DN0OUT

        #source(i1)
        # list ['H', 'D', 'T','3'/'HE3', '4'/'HE4', 'L']
        sources = ['ionization', 'recombination']
        es.source.resize(len(sources))
        #ionization
        i = 0
        name = 'ionization'
        es.source[i].identifier.name = name
        es.source[i].identifierindex = 709
        es.source[i].identifier.description = r'Source from ionization (combined)'
        es.source[i].species.type.name = 'neutral'
        es.source[i].species.type.index = 4
        es.source[i].species.type.description = r'Neutral species'
        es = self.fill_source_data(es, i, name, es.time)

        i = 1
        name = 'recombination'
        es.source[i].identifier.name = name
        es.source[i].identifier.index = 710
        es.source[i].identifier.description = r'Source from recombination'
        es.source[i].species.type.name = 'unspecified'
        es.source[i].species.type.index = 0
        es.source[i].species.type.description = 'All ions'
        es = self.fill_source_data(es, i, name, es.time)


        self.im.putdb_imas(es)
        self.logger.info('Done')


if __name__ == '__main__':
    print("Please use: python convert.py ")
