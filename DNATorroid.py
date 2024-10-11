from __future__ import division
import numpy as np
import math
from PolyAbs import PolyAbs


class DNATorroid(PolyAbs):

    def __init__(self, length, index=0):

        #arc_length = 2 * np.pi * radius
        #length = int(arc_length / .5)

        super(DNATorroid, self).__init__(length, index=index)
        self.rigid_count = 1
        points = self.spiral_points(length)
        for i in range(self.length):
            #self.position[i] = [radius * np.cos(i * 2 * np.pi /length ), radius * np.sin(i * 2 * np.pi/ length), 0]
            self.position[i] = points[i]
            self.charge[i] = 1
            self.type[i] = 'qPm'
            self.mass[i] = 10

    def spiral_points(self,n, arc=.5, separation=1):
        """generate points on an Archimedes' spiral
        with `arc` giving the length of arc between two points
        and `separation` giving the distance between consecutive
        turnings
        - approximate arc length with circle arc at given distance
        - use a spiral equation r = b * phi
        """

        def p2c(r, phi):
            """polar to cartesian
            """
            return [r * math.cos(phi), r * math.sin(phi), 0]

        # yield a point at origin
        points=  [[0,0,0]]

        # initialize the next point in the required distance
        r = arc
        b = separation / (2 * math.pi)
        # find the first phi to satisfy distance of `arc` to the second point
        phi = float(r) / b
        count = 0
        while count < n:
            points.append(p2c(r, phi))
            # advance the variables
            # calculate phi that will give desired arc length at current radius
            # (approximating with circle)
            phi += float(arc) / r
            r = b * phi
            count += 1
        return points