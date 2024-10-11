from __future__ import division
import numpy as np
from Body import Body


class PduaBody(Body):
    """ PduA protein
    points: a list containing positions of points from the body
    charges: a list of charges
    """

    def __init__(self, index=0):
        super(PduaBody, self).__init__(index=0)
        points = list(np.genfromtxt('cg_coord.txt', delimiter='\t'))
        charges = np.genfromtxt('charge.txt')
        #self.leg_sites = [[2.2749, 2.599, 1.9701], [-1.1134, 3.2696, 1.9701], [-3.3883, 0.67062, 1.9701],
         #                 [-2.2749, -2.599, 1.9701], [1.1134, -3.2696, 1.9701], [3.3883, -0.67062, 1.9701]]
        self.leg_sites=[]
        self.arm_sites = []
        z0 = 1.9701
        sq3 = 1 / 2.0 * 3 ** 0.5
        coeff = [[1, 0, 0, -1], [-0.5, sq3, sq3, 0.5], [-0.5, -sq3, -sq3, 0.5],
                 [1, 0, 0, -1], [-0.5, sq3, sq3, 0.5], [-0.5, -sq3, -sq3, 0.5]]

        for i in range(6):
            self.binding_sites.append(tuple([(coeff[i][0] * self.leg_sites[i][0] + coeff[i][1] * self.leg_sites[i][1]),
                                             (coeff[i][2] * self.leg_sites[i][0] + coeff[i][3] * self.leg_sites[i][1]),
                                             z0]))

        for i in range(len(points)):

            if list(points[i]) not in self.leg_sites:
                self.body_sites.append(list(points[i]))
                self.body_charges.append(charges[i])



        mer_body = [[0, 3, 0], [1.76, 2.43, 0]]
        mer_arm_site = [1.43, .463, 0]
        mer_leg_site = [1.76, 2.43, .5]
        mer_binding_site = [1.76, 2.43, -.5]
