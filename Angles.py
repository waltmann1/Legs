from __future__ import division
import numpy as np
import hoomd
from hoomd import md
from Loggable import Loggable


class Angles(Loggable):

    def __init__(self, log_list=None):

        super(Angles, self).__init__(log_list)
        self.log_values = ['angle_harmonic_energy']
        self.names = []
        self.k = []
        self.theta = []
        self.angle_ref = None

        self.names.append('stiff-leg')
        self.theta.append(np.deg2rad(180))
        self.k.append(900)

        self.names.append('polyangle')
        self.theta.append(np.deg2rad(180))
        self.k.append(-1)

    def set_all_harmonic_angles(self, system, reset=False, poly=0):

        if reset:
            self.angle_ref.disable()

        self.angle_ref = hoomd.md.angle.harmonic()
        self.add_to_logger()
        snap = system.take_snapshot(all=True)
        for a in snap.angles.types:
            name = str(a)
            if name == 'polyangle':
                self.angle_ref.angle_coeff.set(name, k=poly, t0=self.theta[self.names.index(name)])
            else:
                self.angle_ref.angle_coeff.set(name, k=self.k[self.names.index(name)],
                                               t0=self.theta[self.names.index(name)])
        return self.angle_ref



