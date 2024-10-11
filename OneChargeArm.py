from __future__ import division
import numpy as np
from PolyAbs import PolyAbs


class OneChargeArm(PolyAbs):

    def __init__(self, index=-1):

        super(OneChargeArm, self).__init__(5, index=index)

        for i in range(0,self.length):
            self.position[i] = [0, 0, -.5 * i]
            if i >= 4:
                self.type[i] = 'qPp'
                self.charge[i] = 1
            else:
                self.type[i] = 'P'
                self.charge[i] = 0
            if not i == 0:
                self.bonds.append([i, i-1])
                self.bond_names.append('polybond')