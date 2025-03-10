from __future__ import division
import numpy as np
from Body import Body


class SV40Body2(Body):

    def __init__(self, index=0):
        super(SV40Body2, self).__init__(index=0)
        r = 4
        ir = 2
        self.body_sites.append([0,0,0])
        for i in range(0,5):
            self.body_sites.append([ir * np.sin(2 * np.pi/5 * i), ir * np.cos(2 * np.pi/5 * i), 0])
            self.body_sites.append([ir * np.sin(2 * np.pi / 5 * i), ir * np.cos(2 * np.pi / 5 * i), 2])
            self.body_sites.append([ir * np.sin(2 * np.pi / 5 * i + 2 * np.pi/10),
                                    ir * np.cos(2 * np.pi / 5 * i + 2 * np.pi/10), 0])
            self.body_sites.append([ir * np.sin(2 * np.pi / 5 * i + 2 * np.pi / 10),
                                    ir * np.cos(2 * np.pi / 5 * i + 2 * np.pi / 10), 2])
            self.arm_sites.append([r/2 * np.sin(2 * np.pi / 5 * i + np.pi/10), r/2 * np.cos(2 * np.pi / 5 * i + np.pi/10), -2])
            #self.leg_sites.append([r * np.sin(2 * np.pi / 5 * i + np.pi/10), r * np.cos(2 * np.pi / 5 * i + np.pi/10), .5])
            self.leg_sites.append([r * np.sin(2 * np.pi / 5 * i + np.pi / 10), r * np.cos(2 * np.pi / 5 * i + np.pi / 10), -1])