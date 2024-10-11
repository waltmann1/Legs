from __future__ import division
import numpy as np
from Body import Body


class HaganBody(Body):

    def __init__(self, index=0):
        super(HaganBody, self).__init__(index=0)
        ir = 2
        r = 4
        self.body_sites.append([0,0,0])
        for i in range(0,5):
            self.body_sites.append([ir * np.sin(2 * np.pi/5 * i), ir * np.cos(2 * np.pi/5 * i), 0])
            self.arm_sites.append([r/2 * np.sin(2 * np.pi / 5 * i + np.pi/10), r/2 * np.cos(2 * np.pi / 5 * i + np.pi/10), -2])
            self.binding_sites.append([r * np.sin(2 * np.pi / 5 * i), r * np.cos(2 * np.pi / 5 * i), 0])

        for i in range(0,4):
            self.body_sites.append(np.average(self.body_sites[i], self.body_sites[i+1]))

        self.body_sites.append(np.average(self.body_sites[0], self.body_sites[4]))