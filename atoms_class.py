
from typing import Dict
import numpy as np
from scipy import constants
from scipy.constants import physical_constants as scipy_phyc


import imas

class atom_base():
    def __init__(self, symbTransp:str=None, ilab:str=''):
        self.a = 0.0
        self.z_n = 0.0
        self.atoms_n = 1
        self.neutral_index = 0
        self.z_ion = 1.0 #ion charge; same as atomic number if only 1 element
        self.label = 'dummy'

        #Transp : key; where key is used in the dictionary below
        self.atoms = {'H':'H', 'D':'D', 'T':'T', '3':'HE3', '4':'HE4'}
        self.nb2 = {'N':'N_BEAM', 'A':'A_BEAM', 'K':'K_BEAM', 'X':'X_BEAM', 'Z':'Z_BEAM'} #'High_Z or heavy ??
        self.batoms = {'H':'H_BEAM', 'D':'D_BEAM', 'T':'T_BEAM', '3':'HE3_BEAM', '4':'HE4_BEAM'}
        self.beams = {**self.batoms, **self.nb2}
        self.ions = {**self.atoms, 'LITH':'LI', 'C':'C'} #some variables: 3=>HE3; 4=>HE4
        self.neutrals = {**self.atoms, 'LITH':'LI'}
        #  NEON, ARGON, KRYPTON, XENON

        self.fusion_ions = {'P':'P_FUSION', 'P14':'P14_FUSION', 'T':'T_FUSION', '3':'HE3_FUSION', '4':'HE4_FUSION'}
        self.impurities = {'BE':'BE', 'W':'W', 'NE':'NE', 'AR':'AR'} # not used

        #Transp key : [IDS label, a, z_n] # from transp periodic_table.f90
        self.entries = {'P_FUSION'   : ['P_FUSION',  1.007276466812, 1.0], #proton
                        'P14_FUSION' : ['P14_FUSION',1.007276466812, 1.0], #proton
                        'H'          : ['H',    1.00782503207,  1.0], #hydrogen
                        'H_BEAM'     : ['H_BEAM',    1.00782503207,  1.0], #hydrogen
                        'D'          : ['D',    2.0141017778,   1.0], #deuterium
                        'D_BEAM'     : ['D_BEAM',    2.0141017778,   1.0], #deuterium
                        'T'          : ['T',    3.0160492777,   1.0], #tritium
                        'T_BEAM'     : ['T_BEAM', 3.0160492777,   1.0], #tritium
                        'T_FUSION'   : ['T_FUSION',3.0160492777,  1.0], #tritium
                        'HE3'        : ['He3',  3.016029,       2.0], #helium3
                        'HE3_BEAM'   : ['HE3_BEAM',  3.016029,       2.0], #helium3
                        'HE3_FUSION' : ['HE3_FUSION',  3.016029,       2.0], #helium3
                        'HE4'        : ['He4',  4.002602,       2.0], #helium4
                        'HE4_BEAM'   : ['HE4_BEAM',  4.002602,       2.0], #helium4
                        'HE4_FUSION' : ['HE4_FUSION',  4.002602,       2.0], #helium4
                        'LITH': ['Li',   6.015122,       3.0], #lithium #Li6
                        'LI7' : ['Li',   7.016004,       3.0], #lithium #Li7
                        'AC'  : ['Ac',  227.0000,     89.0], #Ac iostope average
                        'AG'  : ['Ag',  107.8680,     47.0], #Ag
                        'AL'  : ['Al',   26.9815,     13.0], #Al
                        'AM'  : ['Am',  243.0000,     95.0], #Am
                        'A'   : ['Ar',   39.948,      18.0], #Argon   iostope average
                        'A_BEAM'   : ['AR_BEAM',   39.948,      18.0], #Argon   iostope average
                        'AR'  : ['Ar',   39.9480,     18.0], #Ar
                        'AS'  : ['As',   74.9216,     33.0], #As
                        'AT'  : ['At',  210.0000,     85.0], #At
                        'AU'  : ['Au',  196.9665,     79.0], #Au
                        'B'   : ['B',    10.8100,      5.0], #B
                        'BA'  : ['Ba',  137.3400,     56.0], #Ba
                        'BE'  : ['Be',    9.0122,      4.0], #Be
                        'BH'  : ['Bh',  262.0000,    107.0], #Bh
                        'BI'  : ['Bi',  208.9804,     83.0], #Bi
                        'BK'  : ['Bk',  247.0000,     97.0], #Bk
                        'BR'  : ['Br',   79.9040,     35.0], #Br
                        'C'   : ['C',    12.0110,      6.0], #C
                        'CA'  : ['Ca',   40.0800,     20.0], #Ca
                        'CD'  : ['Cd',  112.4000,     48.0], #Cd
                        'CE'  : ['Ce',  140.1200,     58.0], #Ce
                        'CF'  : ['Cf',  251.0000,     98.0], #Cf
                        'CL'  : ['Cl',   35.4530,     17.0], #Cl
                        'CM'  : ['Cm',  247.0000,     96.0], #Cm
                        'CO'  : ['Co',   58.9332,     27.0], #Co
                        'CR'  : ['Cr',   51.9960,     24.0], #Cr
                        'CS'  : ['Cs',  132.9054,     55.0], #Cs
                        'CU'  : ['Cu',   63.5460,     29.0], #Cu
                        'DB'  : ['Db',  262.0000,    105.0],
                        'BY'  : ['Dy',  162.5000,     66.0],
                        'ER'  : ['Er',  167.2600,     68.0],
                        'ES'  : ['Es',  252.0000,     99.0],
                        'EU'  : ['Eu',  151.9600,     63.0],
                        'F'   : ['F ',   18.9984,      9.0],
                        'FE'  : ['Fe',   55.8470,     26.0],
                        'FM'  : ['Fm',  257.0000,    100.0],
                        'FR'  : ['Fr',  223.0000,     87.0],
                        'GA'  : ['Ga',   69.7200,     31.0],
                        'GD'  : ['Gd',  157.2500,     64.0],
                        'GE'  : ['Ge',   72.5900,     32.0],
                        'H'   : ['H ',    1.0079,      1.0],
                        'HE'  : ['He',    4.0026,      2.0],
                        'HF'  : ['Hf',  178.4900,     72.0],
                        'HG'  : ['Hg',  200.5900,     80.0],
                        'HO'  : ['Ho',  164.9304,     67.0],
                        'HO'  : ['Hs',  265.0000,    108.0],
                        'I'   : ['I ',  126.9045,     53.0],
                        'IN'  : ['In',  114.8200,     49.0],
                        'IR'  : ['Ir',  192.2200,     77.0],
                        'K'   : ['K ',   39.0980,     19.0],
                        'KR'  : ['Kr',   83.8000,     36.0],
                        'K_BEAM'   : ['KR_BEAM',  83.800,       36.0], #Krypton iostope average
                        'LA'  : ['La',  138.9055,     57.0],
                        'LI'  : ['Li',    6.9410,      3.0],
                        'LR'  : ['Lr',  262.0000,    103.0],
                        'LU'  : ['Lu',  174.9700,     71.0],
                        'MD'  : ['Md',  258.0000,    101.0],
                        'MG'  : ['Mg' ,  24.3050,     12.0],
                        'MN'  : ['Mn',   54.9380,     25.0],
                        'MO'  : ['Mo',   95.9400,     42.0],
                        'MT'  : ['Mt',  265.0000,    109.0],
                        'N'   : ['N ',   14.0067,      7.0],
                        'N_BEAM'   : ['NE_BEAM',    20.1790, 10.0], #Neon    iostope average
                        'NA'  : ['Na',   22.9898,     11.0],
                        'NB'  : ['Nb',   92.9064,     41.0],
                        'ND'  : ['Nd',  144.2400,     60.0],
                        'NE'  : ['Ne',   20.1790,     10.0],
                        'NI'  : ['Ni',   58.7000,     28.0],
                        'NO'  : ['No',  259.0000,    102.0],
                        'NP'  : ['Np',  237.0482,     93.0],
                        'O'   : ['O ',   15.9994,      8.0],
                        'OS'  : ['Os',  190.2000,     76.0],
                        'P'   : ['P ',   30.9738,     15.0],
                        'PA'  : ['Pa',  231.0359,     91.0],
                        'PB'  : ['Pb',  207.2000,     82.0],
                        'PD'  : ['Pd',  106.4000,     46.0],
                        'PM'  : ['Pm',  145.0000,     61.0],
                        'PO'  : ['Po',  209.0000,     84.0],
                        'PR'  : ['Pr',  140.9077,     59.0],
                        'PT'  : ['Pt',  195.0900,     78.0],
                        'PU'  : ['Pu',  244.0000,     94.0],
                        'RA'  : ['Ra',  226.0254,     88.0],
                        'RB'  : ['Rb',   85.4678,     37.0],
                        'RE'  : ['Re',  186.2070,     75.0],
                        'RF'  : ['Rf',  261.0000,    104.0],
                        'RH'  : ['Rh',  102.9055,     45.0],
                        'RN'  : ['Rn',  222.0000,     86.0],
                        'RU'  : ['Ru',  101.0700,     44.0],
                        'S'   : ['S',    32.0600,     16.0], #S
                        'SB'  : ['Sb',  121.7500,     51.0], #Sb
                        'SC'  : ['Sc',   44.9559,     21.0], #Se
                        'SE'  : ['Se',   78.9600,     34.0], #Se
                        'SG'  : ['Sg',  263.00,      106.0], #Sg
                        'SI'  : ['Si',   28.0860,     14.0], #Si
                        'SM'  : ['Sm',  150.4000,     62.0], #Sm
                        'SN'  : ['Sn',  118.6900,     50.0], #Sn
                        'SR'  : ['Sr',   87.6200,     38.0], #Sr
                        'TA'  : ['Ta',  180.9479,     73.0], #Ta
                        'TB'  : ['Tb',  158.9254,     65.0], #Tb
                        'TC'  : ['Tc',   47.9000,     22.0], #Tc
                        'TE'  : ['Te',  127.6000,     52.0], #Te
                        'TH'  : ['Th',  232.0381,     90.0], #Th
                        'TI'  : ['Ti',   47.9000,     22.0], #Ti
                        'TL'  : ['Tl',  204.3700,     81.0], #Tl
                        'TM'  : ['Tm',  168.9342,     69.0], #Tm
                        'U'   : ['U',   238.0290,     92.0], #U
                        'V'   : ['V',    50.9414,     23.0], #V
                        'XE'  : ['Xe',  131.300,      54.0], #Xenon
                        'X_BEAM'  : ['XE_BEAM',  131.300, 54.0], #Xenon
                        'W'   : ['W',   183.5000,     74.0], #Tungston av
                        'Y'   : ['Y',    88.9059,     39.0], #Y
                        'YB'  : ['Yb',  173.0400,     70.0], #Yb
                        'ZN'  : ['Zn',   65.3800,     30.0], #Zn
                        'ZR'  : ['Zr',   91.2200,     40.0], #Zr
                        'Z'   : ['Z',   183.5000,     74.0], #Heavy. High Z = W ?
                        'Z_BEAM' : ['Z_BEAM', 183.5000,  74.0], #Heavy. High Z = W ?
                        'TOK' : ['Z',    -1.0,        -1.0], #
                        'X'   : ['X',    -1.0,        -1.0] #
                      }


        self.fill_constants(symbTransp, ilab)

    def fill_constants(self, symbTransp, ilab:str=''):
        key = symbTransp
        if key is None or key=='':
            return
        self.label = self.entries[key][0] + ilab
        self.a = self.entries[key][1]
        self.z_n = self.entries[key][2]


    def find_beam_species_data(self, bsymb:str):
        #bsymb is keys in self.beams
        if bsymb not in self.beams.keys():
            return 0,0,''
        key = (self.beams[bsymb]).upper()
        label = self.entries[key][0]
        a = self.entries[key][1]
        z_n = self.entries[key][2]
        return z_n, a, label


    #dct is one of self.atoms, beams, ions, or neutrals
    def find_species_label(self, sp_a:float, sp_zn:float, dct:Dict)->str :
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
