from __future__ import division
import numpy as np
from PolyAbs import PolyAbs


class Tube(PolyAbs):

    def __init__(self, radius, length, index=0):

        width = .5
        circum = 2 * np.pi * radius
        len = .5
        num_circle = int(circum/len)
        rows = int(length/width)
        length = num_circle * rows
        circle_points = [[radius * np.cos(2 * np.pi * i /num_circle), radius * np.sin(2 * np.pi * i /num_circle ),
                          x * width]
                         for i in range(num_circle) for x in range(rows)]

        super(Tube, self).__init__(length, index=index)
        self.rigid_count = 1
        for ind, point in enumerate(circle_points):
            self.position[ind] = point
            self.mass[ind] = 1
            self.type[ind] = 'qPm'