from __future__ import division
import numpy as np
from BasicSV40Body import BasicSV40Body


class MassAdjustedSV40Body(BasicSV40Body):

    def __init__(self, index=0):
        super(MassAdjustedSV40Body, self).__init__(index=0)

        self.mass = 15
        self.binding_mass = .01
