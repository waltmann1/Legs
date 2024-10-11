from __future__ import division
import numpy as np


class Body(object):

    def __init__(self, index=0):


        self.body_sites = []
        self.arm_sites = []
        self.leg_sites = []
        self.binding_sites = []
        self.type = 'X'
        self.binding_type = 'D'
        self.center = [0,0,0]
        self.mass = 1
        self.binding_mass = 1
        self.index = index
        self.orienation = [1,0,0,0]

    def shift(self, vector):

        for ind, site in enumerate(self.body_sites):
            self.body_sites[ind] = np.add(site, vector)

        for ind, site in enumerate(self.binding_sites):
            self.binding_sites[ind] = np.add(site, vector)

        for ind, site in enumerate(self.arm_sites):
            self.arm_sites[ind] = np.add(site, vector)

        for ind, site in enumerate(self.leg_sites):
            self.leg_sites[ind] = np.add(site, vector)

        self.center = np.add(self.center, vector)

    def align(self, quat):

        for ind, site in enumerate(self.body_sites):
            self.body_sites[ind] = quat.orient(site)
        for ind, site in enumerate(self.arm_sites):
            self.arm_sites[ind] = quat.orient(site)
        for ind, site in enumerate(self.leg_sites):
            self.leg_sites[ind] = quat.orient(site)
        for ind, site in enumerate(self.binding_sites):
            self.binding_sites[ind] = quat.orient(site)

        self.orientation = quat.q

        return [quat.q[3], quat.q[0], quat.q[1], quat.q[2]]






