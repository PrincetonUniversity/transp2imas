
import numpy as np
from scipy import constants

import imas
from utils import *

class ion_base():
    def __init__(self, label):
        self.a = 0.0
        self.z_n = 0.0
        self.atoms_n = 0.0
        self.neutral_index = 0
        self.z_ion = 0.0
        self.label = label
        self.fill_constants(label)

    def fill_constants(self, label, a=None, z_ion=None, z_n=None):
        if label=='hydrogen':
            self.a = 1.00784 #TODO
            self.z_ion = 1.0
            self.z_n = 0.0
        elif label=='deuterium':
            self.a = 2.013553212745 #TODO
            self.z_ion = 1.0
            self.z_n = 1.0
        elif label=='carbon':
            self.a = 12.011 #TODO
            self.z_ion = 6.0
            self.z_n = 6.0