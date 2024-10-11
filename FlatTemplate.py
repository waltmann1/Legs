from __future__ import division
import numpy as np
from PolyAbs import PolyAbs


class FlatTemplate(PolyAbs):

    def __init__(self, width, index=0):

        n = int(width)
        spacing = width / n
        x = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        y = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        xx, yy = np.meshgrid(x, y)
        pos = []
        for i in range(n):
            for j in range(n):
                    pos.append([xx[i][j], yy[i][j], 0])

        super(FlatTemplate, self).__init__(n**2, index=index)
        self.rigid_count = 1

        for ind, point in enumerate(pos):
            self.position[ind] = point
            self.mass[ind] = 100
            self.type[ind] = 'qPm'