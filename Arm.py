from __future__ import division
import numpy as np
from PolyAbs import PolyAbs


class Arm(PolyAbs):

    def __init__(self, index=-1):

        super(Arm, self).__init__(11, index=index)

        for i in range(0,self.length):
            self.position[i] = [0, 0, -.5 * i]
            if not i == 0:
                self.bonds.append([i, i-1])
                self.bond_names.append('polybond')
            if self.length - i < 4:
                self.type[i] = 'qPp'
                self.charge[i] = 1
