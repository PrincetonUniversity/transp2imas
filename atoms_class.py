
import numpy as np
from scipy import constants

import imas

class atom_base():
    def __init__(self, symbTransp:str=None):
        self.a = 0.0
        self.z_n = 0.0
        self.atoms_n = 1
        self.neutral_index = 0
        self.z_ion = 1.0 #ion charge; ionization state
        self.label = 'dummy'

        #Transp : IDS
        self.atoms = {'H':'H', 'D':'D', 'T':'T', 'HE3':'He3', 'HE4':'He4'}
        self.nb2 = {'N':'Ne', 'A':'Ar', 'K':'Kr', 'X':'Xe', 'Z':'HI_Z'} #'High_Z or heavy ??
        self.beams = {**self.atoms, **self.nb2} #self.atoms.copy()
        self.ions = {**self.atoms, 'LITH':'Li', 'C':'C'}
        self.neutrals = {**self.atoms, 'LITH':'Li'}
        #  NEON, ARGON, KRYPTON, XENON
        self.fusion_ions = {'P':'H', 'P14':'H', 'T':'T', '3':'He3', '4':'He4'}

        #Transp symbol : [IDS label, a, z_n]
        self.entries = {'P'   : ['H',   1.00727647,     1.0], #proton
                        'P14' : ['H',   1.00727647,     1.0], #proton
                        'H'   : ['H',   1.00784,        1.0], #hydrogen
                        'D'   : ['D',   2.013553212745, 1.0], #deuterium
                        'T'   : ['T',   3.01604928,     1.0], #tritium
                        'HE3' : ['He3', 3.016029,       2.0], #helium3
                        'HE4' : ['He4', 4.002603254,    2.0], #helium4
                        'LITH': ['Li',  6.941,          4.0], #lithium
                        'C'   : ['C',   12.011,         6.0], #carbon
                        'N'   : ['Ne',  20.1797,       10.0], #Neon
                        'A'   : ['Ar',  39.948,        18.0], #Argon
                        'K'   : ['Kr',  83.798,        36.0], #Krypton
                        'X'   : ['Xe',  131.293,       54.0], #Xeon
                        'Z'   : ['Z',   -1.0,          -1.0] #Heavy. High Z
                       }


        self.fill_constants(symbTransp)

    def fill_constants(self, symbTransp):
        lab = symbTransp
        if lab is None or lab=='':
            return
        self.label = self.entries[lab][0]
        self.a = self.entries[lab][1]
        self.z_n = self.entries[lab][2]

    #dct is one of self.atoms, beams, ions, or neutrals
    def find_species_label(self, sp_a:float, sp_zn:float, dct:dict)->str :
        eps = 0.1 #round off
        #dct: {key is Transp symbol, value is label}
        for ts, values in self.entries.items():
            lab, aa, zn = values
            if (ts in dct.keys()) and (lab in dct.values()) and \
                (abs(aa-sp_a)<eps) and (abs(zn-sp_zn)<eps):
                return lab
        return None
