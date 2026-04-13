

from nc2imas import *

class imas_edge_sources:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger

        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\nedge_sources')

        ids = self.im.ids.get('edge_sources')
        cm2m = 0.01

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        timep = self.im.get_ncvar('TIME3')
        ntimep = len(timep)
        self.im.checkErrExit(timesc, timep)


        self.logger.info(f'edge sources')
        ids.midplane.name = 'magnetic_axis'
        ids.midplane.index = 1
        txt = r'Midplane defined by the height of magnetic axis'
        esc.midplane.description = txt + r' equilibrium/time_slice/global_quantities/magnetic_axis/z'

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
        #recombination=recycling+.. ?
        sources = ['ionization', 'recombination']
        ids.source.resize(len(sources))

        #ionization
        name = 'ionization'
        es = ids.source[0]
        es.name = name
        es.index = 709
        es.description = r'Source from ionization (combined)'
        es.species.type.name = 'neutral'
        es.species.type.index = 4
        es.species.type.description = r'Neutral species of H+D+T'
        es.species.neutral.label = ''
        es.species.neutral.element.resize(3)
        for e, symb in enumerate(['H', 'D', 'T']):
            at = atom_base(symb)
            es.species.neutral.element[e].a = at.a
            es.species.neutral.element[e].z_n = at.z_n
            es.species.neutral.element[e].atoms_n = at.atoms_n

        name = 'recombination'
        es = ids.source[1]
        es.name = name
        es.index = 710
        es.description = r'Source from recombination'
        es.species.type.name = 'unspecified'
        es.species.type.index = 0
        es.species.type.description = 'All ions'







        self.im.putdb_imas(esc)
        self.logger.info('Done')


if __name__ == '__main__':
    print("Please use: python convert.py ")
