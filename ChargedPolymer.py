from __future__ import division
import numpy as np
import math
from PolyAbs import PolyAbs


class ChargedPolymer(PolyAbs):

    def __init__(self, length, index=-1, angles=None):

        super(ChargedPolymer, self).__init__(length, index=index)
        self.rigid_count = 0
        points = self.get_points(length)
        for ind, point in enumerate(points):
            self.position[ind] = point
            self.type[ind] = 'qPm'
            self.mass[ind] = 1
            if ind != 0:
                self.bonds.append([ind, ind-1])
                self.bond_names.append('polybond')
            if angles is not None and ind > 1:
                self.angles.append([ind, ind - 1, ind - 2])
                self.angle_names.append('polyangle')




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
        while count < n - 1:
            points.append(p2c(r, phi))
            # advance the variables
            # calculate phi that will give desired arc length at current radius
            # (approximating with circle)
            phi += float(arc) / r
            r = b * phi
            count += 1
        return points

    def get_points(self, length, max_layer=200, seperation=1):

        per_layer = length
        layers = 1
        if length > max_layer:
            layers = int(np.ceil(length/200))
            per_layer = int(np.ceil(length/layers))
        points = self.spiral_points(per_layer)
        sep_count = 0
        print(per_layer, layers)

        for i in range(layers-1):
            sep_count += 1
            new = np.add(self.spiral_points(per_layer), [0, 0, sep_count])
            if i%2 == 0:
                new = list(reversed(new))
            for new_point in new:
                points.append(new_point)

        points = [np.add(point, [0, 0, -sep_count/2]) for point in points]
        points = points[:length]
        return points


