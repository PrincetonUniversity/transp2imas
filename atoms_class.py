
import numpy as np
from scipy import constants
from scipy.constants import physical_constants as scipy_phyc


import imas

class atom_base():
    def __init__(self, symbTransp:str=None, ilabel:str=''):
        self.a = 0.0
        self.z_n = 0.0
        self.atoms_n = 1
        self.neutral_index = 0
        self.z_ion = 1.0 #ion charge; ionization state
        self.label = 'dummy'

        #Transp : IDS
        self.atoms = {'H':'H', 'D':'D', 'T':'T', '3':'He3', '4':'He4'}
        self.nb2 = {'N':'Ne', 'A':'Ar', 'K':'Kr', 'X':'Xe', 'Z':'HI_Z'} #'High_Z or heavy ??
        self.beams = {**self.atoms, **self.nb2} #self.atoms.copy()
        self.ions = {**self.atoms, 'LITH':'Li', 'C':'C'}
        self.neutrals = {**self.atoms, 'LITH':'Li'}
        #  NEON, ARGON, KRYPTON, XENON
        self.fusion_ions = {'P':'H', 'P14':'H', 'T':'T', '3':'He3', '4':'He4'}

        #Transp symbol : [IDS label, a, z_n] # from periodic_table_mod
        self.entries = {'P'   : ['H',   1.007276466812, 1.0], #proton
                        'P14' : ['H',   1.007276466812, 1.0], #proton
                        'H'   : ['H',   1.00782503207,  1.0], #hydrogen
                        'D'   : ['D',   2.0141017778,   1.0], #deuterium
                        'T'   : ['T',   3.0160492777,   1.0], #tritium
                        'HE3' : ['He3', 3.016029,       2.0], #helium3
                        'HE4' : ['He4', 4.002602,       2.0], #helium4
                        'LITH': ['Li',  6.015122,       3.0], #lithium #Li6
                        'LI7' : ['Li',  7.016004,       3.0], #lithium #Li7
                        'LI'  : ['Li',  6.9410,         3.0], #lithium average isotope
                        'BE'  : ['Be',  9.0122,         4.0], #beryllium average
                        'C'   : ['C',   12.011,         6.0], #carbon  #iostope average
                        'N'   : ['Ne',  20.1790,       10.0], #Neon    #iostope average
                        'NE'  : ['Ne',  20.1790,       10.0], #Neon    #iostope average
                        'A'   : ['Ar',  39.948,        18.0], #Argon   #iostope average
                        'K'   : ['Kr',  83.800,        36.0], #Krypton #iostope average
                        'X'   : ['Xe',  131.300,       54.0], #Xeon    #iostope average
                        'W'   : ['W',   183.5000,      74.0], #Tungston av
                        'Z'   : ['Z',   -1.0,          -1.0], #Heavy. High Z
                        'TOK' : ['Z',   -1.0,          -1.0] #
                       }


        self.fill_constants(symbTransp, ilabel)

    def fill_constants(self, symbTransp, ilab:str=''):
        lab = symbTransp
        if lab is None or lab=='':
            return
        self.label = self.entries[lab][0] + ilab
        self.a = self.entries[lab][1]
        self.z_n = self.entries[lab][2]

    #dct is one of self.atoms, beams, ions, or neutrals
    def find_species_label(self, sp_a:float, sp_zn:float, dct:dict)->str :
        if not (sp_a and sp_zn and dct):
            return None
        eps = 0.1 #round off
        #dct: {key is Transp symbol, value is label}
        for ts, values in self.entries.items():
            lab, aa, zn = values
            if (ts in dct.keys()) and (lab in dct.values()) and \
                (abs(aa-sp_a)<eps) and (abs(zn-sp_zn)<eps):
                return lab
        return None
