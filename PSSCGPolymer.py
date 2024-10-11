from __future__ import division
import numpy as np
import math
from ChargedPolymer import ChargedPolymer


class PSS427(ChargedPolymer):

    def __init__(self, length, index=-1, angles=None):

        super(ChargedPolymer, self).__init__(length=107, index=index)
        self.rigid_count = 0
        points = self.get_points(length)
        for ind, point in enumerate(points):
            self.position[ind] = point
            self.type[ind] = 'qPm_1'
            self.mass[ind] = 4
            if ind != 0:
                self.bonds.append([ind, ind-1])
                self.bond_names.append('pssbond')
            if angles is not None and ind > 1:
                self.angles.append([ind, ind - 1, ind - 2])
                self.angle_names.append('polyangle')
