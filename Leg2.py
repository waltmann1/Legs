from __future__ import division
import numpy as np
from PolyAbs import PolyAbs


class Leg2(PolyAbs):

    def __init__(self, index=-1):

        super(Leg2, self).__init__(4, index=index)
        for i in range(0, self.length):
            self.position[i] = [1 * i, 0, 0]
            if not i == 0:
                self.bonds.append([i, i - 1])
                self.bond_names.append('longbond')
            if i > 0 and i < self.length - 2:
                self.angles.append([i, i+1, i+2])
                self.angle_names.append('stiff-leg')
            if self.length - 1 == i:
                self.type[i] = 'C'
            else:
                self.type[i] = 'H'
            if i == 0:
                self.type[i] = 'D'
