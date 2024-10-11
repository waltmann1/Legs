from __future__ import division
import numpy as np
import hoomd
from hoomd import md
from numpy import linalg as la
from Loggable import Loggable


def qpm_size(system):
    size = 1.0
    types = system.particles.types

    for ind in range(len(system.particles)):
        part = system.particles.get(ind)
        if part.type == "qPm":
            if part.body == -1 or part.body > 10000:
                return 0.5

    return size

class LJRepulsive(Loggable):
    def __init__(self, log_list=None):

        super(LJRepulsive, self).__init__(log_list)
        self.log_values = ['pair_table_energy']

        self.epsilon = [1, 1, 1, 1, 1, 1, 1, 0]


        self.sigma = [4, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5]
        self.names = ['X', 'C', 'D', 'P', 'H', 'qPp', 'qPm', 'S']

        self.lj_repulsive_pair = None

    def set_lj_repulsive(self, neighbor_list, system, table_cd=False, helix=None, effective_charge=None, debye=None,
                         sphere=None, depletion=False, attract=None):
        cut = 2
        self.sigma[-1] = qpm_size(system)
        self.lj_repulsive_pair = hoomd.md.pair.table(width=1000, nlist=neighbor_list)
        self.add_to_logger()
        for t1 in system.particles.types:
            for t2 in system.particles.types:
                t1 = str(t1)
                t2 = str(t2)
                if t1 == 'C' and t2 == 'C':
                    self.lj_repulsive_pair.pair_coeff.set(str(t1), str(t2), rmin=10e-5, rmax=2, func=LJRepulsive_pair,
                                                          coeff=dict(sigma=2, epsilon=1))
                elif t1 == 'D' and t2 == 'D':
                    self.lj_repulsive_pair.pair_coeff.set(str(t1), str(t2), rmin=10e-5, rmax=2, func=LJRepulsive_pair,
                                                          coeff=dict(sigma=2, epsilon=1))
                elif t1 == 'X' and t2 in self.names or t2 == 'X' and t1 in self.names:
                    self.lj_repulsive_pair.pair_coeff.set(str(t1), str(t2), rmin=10e-5, rmax=2, func=LJRepulsive_pair,
                                      coeff=dict(sigma=2,
                                                 epsilon=np.sqrt(self.epsilon[self.names.index(t1)] *
                                                                  self.epsilon[self.names.index(t2)])))
                elif (t1 == 'C' and t2 == 'D' or t1 == 'D' and t2 == 'C') and table_cd:
                    self.lj_repulsive_pair.pair_coeff.set(str(t1), str(t2), rmin=10e-5, rmax=1, func=double_well2,
                                                          coeff=dict(well1=8, well2=10, barrier=8))
                    #self.lj_repulsive_pair.set_from_file(str(t1), str(t2), filename='table_potential.txt')
                elif self.is_center(t1) and t2 == 'qPp' and debye is not None or self.is_center(t2) and t1 == 'qPp' and debye is not None:
                    extra_cut = np.min([3 * debye, 5])
                    #self.lj_repulsive_pair.pair_coeff.set(str(t1), str(t2), rmin=10e-5, rmax=sphere + 1 + extra_cut,
                    #                                      func=ci_release,
                    #                                      coeff=dict(sphere_radius=sphere,
                    #                                                 effective_charge=effective_charge,
                    #                                                 debye_length=debye))
                    self.lj_repulsive_pair.pair_coeff.set(str(t1), str(t2), rmin=10e-5, rmax=sphere + 1.5,
                                                          func=short_range,
                                                          coeff=dict(sphere_radius=sphere,
                                                                     attract=attract))
                elif self.is_center(t1) and t2 == 'S' or self.is_center(
                        t2) and t1 == 'S':
                    self.lj_repulsive_pair.pair_coeff.set(str(t1), str(t2), rmin=10e-5, rmax=sphere + 2,
                                                          func=lj_sphere,
                                                          coeff=dict(sphere_radius=sphere,
                                                                     effective_charge=.5,
                                                                     debye_length=.0001))

                elif t1 == 'H' and t2 == 'H' and helix is not None:
                    self.lj_repulsive_pair.pair_coeff.set(str(t1), str(t2), rmin=10e-5, rmax=2, func=H_potential,
                                                          coeff=dict(sigma=1, epsilon=helix))
                elif t1 == 'X' and t2 == 'X' and depletion:
                    self.lj_repulsive_pair.pair_coeff.set(str(t1), str(t2), rmin=10e-5, rmax=2, func=LJRepulsive_pair,
                                                          coeff=dict(sigma=0, epsilon=0))

                elif t1 in self.names and t2 in self.names:
                    self.lj_repulsive_pair.pair_coeff.set(str(t1), str(t2), rmin=10e-5, rmax=2, func=LJRepulsive_pair,
                                      coeff=dict(sigma=(self.sigma[self.names.index(t1)] +
                                                        self.sigma[self.names.index(t2)]) / 2,
                                                 epsilon=np.sqrt(self.epsilon[self.names.index(t1)] *
                                                                  self.epsilon[self.names.index(t2)])))

                else:
                    self.lj_repulsive_pair.pair_coeff.set(str(t1), str(t2), rmin=10e-5, rmax=2, func=LJRepulsive_pair,
                                                          coeff=dict(sigma=0, epsilon=0))

        return self.lj_repulsive_pair

    def is_center(self, string):

        if len(string) > 6 and string[:6] == 'center':
            return True
        return False


class Yukawa(Loggable):

    def __init__(self, log_list=None, debye=1, total_charge=None, effective_charge=1):
        super(Yukawa, self).__init__(log_list)
        self.log_values = ['pair_yukawa_energy']
        self.lb = .7
        self.sigma = [4, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5]
        self.names = ['X', 'C', 'D', 'P', 'H', 'qPp', 'qPm', 'S']
        self.charge = [0, 0, 0, 0, 0, 1, -1, 0]
        self.kappa = 1/debye
        self.yukawa = None
        self.total_charge = total_charge
        self.effective_charge = effective_charge


    def set_yukawa(self, neighbor_list, system):

        yuk = hoomd.md.pair.yukawa(r_cut=3 / self.kappa, nlist=neighbor_list)
        self.add_to_logger()
        self.sigma[-1] = qpm_size(system)

        if self.total_charge is not None:
            count = 0
            for t in system.particles:
                t = str(t.type)
                if t == 'qPm':
                    count += 1
            if count > 0:
                self.charge[self.names.index('qPm')] = self.total_charge / count
            else:
                self.charge[self.names.index('qPm')] = 0
        for t1 in system.particles.types:
            for t2 in system.particles.types:
                t1 = str(t1)
                t2 = str(t2)
                if t1 in self.names and t2 in self.names:
                    sigma = .5 * (self.sigma[self.names.index(t1)] + self.sigma[self.names.index(t2)])
                    q1 = self.charge[self.names.index(t1)] * self.effective_charge
                    q2 = self.charge[self.names.index(t2)] * self.effective_charge
                    #eps = q1 * q2 * self.lb * np.exp(self.kappa * sigma) / (1/self.kappa + sigma)
                    eps = q1 * q2 * self.lb
                    #eps = eps * 100
                    yuk.pair_coeff.set(t1, t2, epsilon=eps, kappa=self.kappa)
                else:
                    yuk.pair_coeff.set(t1, t2, epsilon=0, kappa=self.kappa)
        self.yukawa = yuk

class LJSpecial(Loggable):
    def __init__(self, log_list=None, energy=10):

        super(LJSpecial, self).__init__(log_list)

        self.log_values = ['pair_lj_energy']

        self.names = ['X', 'C', 'D', 'P', 'H', 'qPm', 'qPp', 'S']

        self.epsilon = [1, 1, 1, 1, 1, 1, 1]

        self.sigma = [4, 0.5, 0.5, .5, .5, .5, 1, .5]

        self.lj_pair = None

        self.energy = energy

    def set_lj_special(self, neighbor_list, system, depletion=False):
        cut = 3
        self.lj_pair = hoomd.md.pair.lj(r_cut=cut, nlist=neighbor_list)
        self.add_to_logger()
        for t1 in system.particles.types:
            for t2 in system.particles.types:
                t1 = str(t1)
                t2 = str(t2)
                if t1 == 'C' and t2 == 'D' or t1 == 'D' and t2 == 'C':
                    self.lj_pair.pair_coeff.set(str(t1), str(t2), epsilon=self.energy, sigma=.5)
                elif t1 == 'S' and not self.is_center(t2) and not t2 == 'qPm' or t2 == 'S' and not self.is_center(t1) and not t1 == 'qPm':
                    self.lj_pair.pair_coeff.set(str(t1), str(t2), epsilon=1, sigma=.25 + (self.sigma[self.names.index(t2)]/2))
                elif t1 == 'X' and t2 == 'X' and depletion:
                    self.lj_pair.pair_coeff.set(str(t1), str(t2), epsilon=.5, sigma=4)
                else:
                    self.lj_pair.pair_coeff.set(str(t1), str(t2), epsilon=0, sigma=1)

    def is_center(self, string):

        if len(string) >= 6 and string[:6] == 'center':
            return True
        return False

def LJRepulsive_pair(r,rmin, rmax, sigma, epsilon):

    if r < sigma:
        V = epsilon * ((sigma / r) ** 12 - 1)
        F = epsilon * (sigma/r) ** 13
    else:
        V = 0
        F = 0
    return (V,F)


def quad(r, sigma, epsilon):
    if r < (sigma + .25) and r > (sigma - .25):
        V = 16 * epsilon * (r - 1) ** 2 - epsilon
        F = - (16 * epsilon * (2 * r - 2))
    else:
        V = 0
        F = 0
    return (V, F)


def H_potential(r, rmin, rmax, sigma, epsilon):
    if r > .5:
        return quad(r, sigma, epsilon)
    else:
        return LJRepulsive_pair(r,0,0, .5, 1)


def gaussian(std, mean, x):

    return np.exp(-(x-mean)**2 / (2 * std**2))

def gaussian_prime(std, mean, x):

    return 2 * (x-mean)/(2 * std**2) * gaussian(std, mean, x)

def double_well(r, rmin, rmax, well1, well2):

    if r > 1:
        V =0
        F =0
    elif r > .65:
        V = -well1 * gaussian(.05, .8,  r)
        F = - well1 * gaussian_prime(.05, .8, r)

    elif r > .35:
        V = -well2 * gaussian(.05, .5, r)
        F = - well2 * gaussian_prime(.05, .5, r)
    else:
        return LJRepulsive_pair(r, 0, 0, .35, 1)
    return (V,F)

def double_well2(r, rmin, rmax, well1, well2, barrier):

    if r > 1:
        V =0
        F =0
    elif r > .72:
        V = -well1 * gaussian(.04, .84,  r)
        F = - well1 * gaussian_prime(.04, .84, r)
    elif r > .48:
        V = barrier * gaussian(.04, .6, r)
        F = barrier * gaussian_prime(.04, .6, r)
    elif r > .24:
        V = -well2 * gaussian(.04, .36, r)
        F = - well2 * gaussian_prime(.04, .36, r)
    else:
        return LJRepulsive_pair(r, 0, 0, .24, 1)
    return (V,F)

def ci_release(r, rmin, rmax, sphere_radius, effective_charge, debye_length):

    V = 0
    F = 0
    if r > sphere_radius + 1:
        V = - 2 * (1 - effective_charge)  * np.exp(-(r - (sphere_radius + 1))/ debye_length)
        F = - 2 * (1 - effective_charge) / debye_length * np.exp(-(r - (sphere_radius + 1))/ debye_length)
    else:
        V = - 2 * (1 - effective_charge)
        F = 0

    print("ci_release", r,V,F, sphere_radius)
    return V, F


def short_range(r, rmin, rmax, sphere_radius, attract):

    V = 0
    F = 0
    edge1 = sphere_radius + 1.0
    edge2 = sphere_radius + 0.9
    if r > edge1:
        V = 0
        F = 0
    elif r > edge2:
        V = attract * (r - (edge1)) / (edge1 - edge2)
        F = -attract / (edge1 - edge2)
    else:
        V = - attract
        F = 0

    print("short range",r,V,F, sphere_radius)
    return V, F

def lj_sphere(r, rmin, rmax, sphere_radius, effective_charge, debye_length):

    V = 0
    F = 0
    if r > sphere_radius + 1:
        V = - 1 + np.power((r - 1 - sphere_radius), 6)
        F = -6 * np.power((r - 1 - sphere_radius), 5)
    else:
        V = - 1
        F = 0
    print("lj_sphere", r, V, F, sphere_radius)

    return V, F