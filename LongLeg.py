from __future__ import division
import numpy as np
from PolyAbs import PolyAbs
from HelixLeg import HelixLeg


class LongLeg(HelixLeg):

    def __init__(self, index=-1):

        super(LongLeg, self).__init__(index=index)

        self.length = 11

        self.position.append(np.add(self.position[-1], [2, 0, 0]))
        self.position.append(np.add(self.position[-1], [2, 0, 0]))
        self.bond_names.append('connect_bond')
        self.bonds.append([8, 9])
        self.bond_names.append('connect_bond')
        self.bonds.append([9, 10])
        self.type.append('S')
        self.type.append('S')
        self.mass.append(self.mass[-1])
        self.mass.append(self.mass[-1])
        self.charge.append(0)
        self.charge.append(0)
