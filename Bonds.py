from __future__ import division
import numpy as np
import hoomd
from hoomd import md
from Loggable import Loggable


class Bonds(Loggable):

    def __init__(self, log_list=None):

        super(Bonds, self).__init__(log_list)
        self.log_values = ['bond_harmonic_energy']
        self.names = []
        self.k = 300
        self.r0 = []
        self.bond_ref = None

        self.names.append('polybond')
        self.r0.append(.5)
        self.k = [300]

        self.names.append('pssbond')
        self.r0.append(1.0)
        self.k.append(300)

        self.names.append('connect_bond')
        self.r0.append(2)
        self.k.append(1)


        self.names.append('longbond')
        self.r0.append(1)
        self.k.append(300)

    def set_all_harmonic_bonds(self, system, reset=False):
        """

        :param system: the system that needs the parameters set
        :return: reference to the harmonic bond object
        """
        if reset:
            self.bond_ref.disable()

        self.bond_ref = hoomd.md.bond.harmonic()
        self.add_to_logger()
        snap = system.take_snapshot(all=True)
        for b in snap.bonds.types:
            name = str(b)
            self.bond_ref.bond_coeff.set(name, k=self.k[self.names.index(name)], r0=self.r0[self.names.index(name)])
        del snap
        return self.bond_ref


