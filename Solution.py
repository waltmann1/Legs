from __future__ import division
import numpy as np
from numpy import linalg as la
import hoomd
import copy as cp


class Solution(object):

    def __init__(self, mers, chains, box_length=None, flat=False):

        self.mers = [cp.deepcopy(mer) for mer in mers]
        self.chains = [cp.deepcopy(chain) for chain in chains]
        self.dump_context = None
        self.rigid_count = int(np.sum([mer.rigid_count for mer in self.mers]) + np.sum([chain.rigid_count for chain in self.chains]))
        self.num_particles = np.sum([mer.total_particles for mer in mers]) + np.sum([len(chain.mass) for chain in self.chains])
        self.num_particles = int(self.num_particles)
        self.reindex()
        self.name = "assembly"
        if box_length is not None:
            self.box_length = int(box_length)
        else:
            self.box_length = None

        self.flat = flat

        self.o_list = None

    def reindex(self):
        tag = 0
        for mer in self.mers:
            if mer.index != -1:
                mer.index = tag
                tag += 1
        for chain in self.chains:
            if chain.index != -1:
                chain.index = tag
                tag += 1

    def create_system(self):
        """

        :return: system object
        """
        if self.dump_context is None:
            self.dump_context = hoomd.context.initialize("")
        b_types = self.parse_bond_types()
        p_types = self.parse_particle_names() + ["center" + str(i) for i in range(len(self.mers), self.rigid_count)] \
                  + ['center'] + ['centerflat']
        a_types = self.parse_angle_types()
        if self.box_length is None:
            cen = self.geometric_center()
            self.shift([c * -1 for c in cen])
            l = int(self.max_radius() * 2 + 25)
        else:
            l = self.box_length

        snap = hoomd.data.make_snapshot(self.num_particles + self.rigid_count, particle_types=p_types,
                                        bond_types=b_types, angle_types=a_types,
                                        box=hoomd.data.boxdim(L=l))
        snap.bonds.resize(0)

        for x in range(self.rigid_count):
            snap.particles.position[x] = self.center_of_mass(x)
            snap.particles.mass[x] = np.sum([chain.mass[i] for chain in self.chains for i in range(chain.length)])
            if x < len(self.mers):
                snap.particles.typeid[x] = p_types.index("center")
            elif not self.flat:
                snap.particles.typeid[x] = p_types.index("center" + str(x))
            else:
                snap.particles.typeid[x] = p_types.index("centerflat")
            snap.particles.body[x] = x
            snap.particles.moment_inertia[x] = self.moment_inertia(x)
        tag = self.rigid_count

        for mer in self.mers:
            for x in range(len(mer.body.body_sites)):
                snap.particles.position[x + tag] = mer.body.body_sites[x]
                snap.particles.mass[x + tag] = mer.body.mass
                snap.particles.typeid[x + tag] = p_types.index(mer.body.type)
                snap.particles.body[x + tag] = mer.index
                snap.particles.charge[x + tag] = 0
            tag += len(mer.body.body_sites)

            for x in range(len(mer.body.binding_sites)):
                snap.particles.position[x + tag] = mer.body.binding_sites[x]
                snap.particles.mass[x + tag] = mer.body.binding_mass
                snap.particles.typeid[x + tag] = p_types.index(mer.body.binding_type)
                snap.particles.body[x + tag] = mer.index
                snap.particles.charge[x + tag] = 0

            tag += len(mer.body.binding_sites)

            for arm in mer.arms:
                for x in range(len(arm.bonds)):
                    #bodies = [chain.body[i] for i in arm.bonds[x] if chain.body[i] != -1]
                    #if len(bodies) == len(set(bodies)):
                    bond_number = snap.bonds.N + 1
                    snap.bonds.resize(bond_number)
                    snap.bonds.group[bond_number - 1] = np.add(arm.bonds[x], tag)
                    snap.bonds.typeid[bond_number - 1] = b_types.index(arm.bond_names[x])

                for x in range(arm.length):
                    snap.particles.position[x + tag] = arm.position[x]
                    snap.particles.mass[x + tag] = arm.mass[x]
                    snap.particles.typeid[x + tag] = p_types.index(arm.type[x])
                    #
                    # check that body
                    #
                    if x == 0:
                        snap.particles.body[x + tag] = mer.index
                    else:
                        snap.particles.body[x + tag] = -1
                    #snap.particles.charge[x + tag] = arm.charge[x]
                    snap.particles.charge[x + tag] = 0
                tag += arm.length

            for leg in mer.legs:
                for x in range(len(leg.bonds)):
                    #bodies = [chain.body[i] for i in arm.bonds[x] if chain.body[i] != -1]
                    #if len(bodies) == len(set(bodies)):
                    bond_number = snap.bonds.N + 1
                    snap.bonds.resize(bond_number)
                    snap.bonds.group[bond_number - 1] = np.add(leg.bonds[x], tag)
                    snap.bonds.typeid[bond_number - 1] = b_types.index(leg.bond_names[x])

                for x in range(len(leg.angles)):
                    #bodies = [chain.body[i] for i in arm.bonds[x] if chain.body[i] != -1]
                    #if len(bodies) == len(set(bodies)):
                    angle_number = snap.angles.N + 1
                    snap.angles.resize(angle_number)
                    snap.angles.group[angle_number - 1] = np.add(leg.angles[x], tag)
                    snap.angles.typeid[angle_number - 1] = a_types.index(leg.angle_names[x])


                for x in range(leg.length):
                    snap.particles.position[x + tag] = leg.position[x]
                    snap.particles.mass[x + tag] = leg.mass[x]
                    snap.particles.typeid[x + tag] = p_types.index(leg.type[x])
                    #
                    # check that body
                    #
                    if x == 0:
                        snap.particles.body[x + tag] = mer.index
                    else:
                        snap.particles.body[x + tag] = -1
                    #snap.particles.charge[x + tag] = leg.charge[x]
                    snap.particles.charge[x + tag] = 0
                tag += leg.length
        for chain in self.chains:
            for x in range(len(chain.bonds)):
                #bodies = [chain.body[i] for i in chain.bonds[x] if chain.body[i] != -1]
                #if len(bodies) == len(set(bodies)):
                bond_number = snap.bonds.N + 1
                snap.bonds.resize(bond_number)
                snap.bonds.group[bond_number - 1] = np.add(chain.bonds[x], tag)
                snap.bonds.typeid[bond_number - 1] = b_types.index(chain.bond_names[x])

            for x in range(chain.length):
                snap.particles.position[x + tag] = chain.position[x]
                snap.particles.mass[x + tag] = chain.mass[x]
                snap.particles.typeid[x + tag] = p_types.index(chain.type[x])
                snap.particles.body[x + tag] = chain.index
                #snap.particles.charge[x + tag] = chain.charge[x]
                snap.particles.charge[x + tag] = 0
            tag += chain.length

        sys = hoomd.init.read_snapshot(snap)
        if self.box_length is None:
            self.shift([c for c in cen])

        return sys


    def parse_bond_types(self):

        a_types = []
        for chain in self.chains:
            for c in chain.bond_names:
                if c not in a_types:
                    a_types.append(c)
        for mer in self.mers:
            for arm in mer.arms:
                for c in arm.bond_names:
                    if c not in a_types:
                        a_types.append(c)
            for leg in mer.legs:
                for c in leg.bond_names:
                    if not c in a_types:
                        a_types.append(c)

        return a_types

    def parse_angle_types(self):

        a_types = []
        for chain in self.chains:
            for c in chain.angle_names:
                if c not in a_types:
                    a_types.append(c)
        for mer in self.mers:
            for arm in mer.arms:
                for c in arm.angle_names:
                    if c not in a_types:
                        a_types.append(c)
            for leg in mer.legs:
                for c in leg.angle_names:
                    if not c in a_types:
                        a_types.append(c)

        return a_types

    def parse_particle_names(self):
        a_types = []
        for chain in self.chains:
            for c in chain.type:
                if c not in a_types:
                    a_types.append(c)
        for mer in self.mers:
            if mer.body.type not in a_types:
                a_types.append(mer.body.type)
            if mer.body.binding_type not in a_types:
                a_types.append(mer.body.binding_type)
            for arm in mer.arms:
                for c in arm.type:
                    if c not in a_types:
                        a_types.append(c)
            for leg in mer.legs:
                for c in leg.type:
                    if not c in a_types:
                        a_types.append(c)

        return a_types


    def geometric_center(self):

        pos = [np.mean(c.position, axis=0) for c in self.chains]
        pos_mers = [mer.body.center for mer in self.mers]
        pos = pos + pos_mers
        weights = [len(c.position) for c in self.chains]
        weights_mers = [mer.total_particles for mer in self.mers]
        weights = weights + weights_mers
        total = np.sum(weights)
        weights = [float(weight)/float(total) for weight in weights]
        center = np.array([0, 0, 0])
        for ind, p in enumerate(pos):
            #print(weights[ind])
            center = np.add(center, np.multiply(p, weights[ind]))
        return center

    def max_radius(self):

        t = [-1 * c for c in self.geometric_center()]
        #print(t)
        self.shift(t)
        rs = [la.norm(pos) for chain in self.chains for pos in chain.position]
        rs_arms = [la.norm(pos) for mer in self.mers for arm in mer.arms for pos in arm.position]
        rs_legs = [la.norm(pos) for mer in self.mers for leg in mer.legs for pos in leg.position]
        rs = rs + rs_arms + rs_legs
        self.shift([-1 * g for g in t])
        return np.max(rs)

    def shift(self, vector):

        for mer in self.mers:
            mer.shift(vector)

        for chain in self.chains:
            chain.shift(vector)

    def center_of_mass(self, center_index):
        for mer in self.mers:
            if mer.index == center_index:
                return mer.rigid_center_of_mass()
        for chain in self.chains:
            if chain.index == center_index:
                return chain.rigid_center_of_mass()

    def total_rigid_mass(self, center_index):
        for mer in self.mers:
            if mer.index == center_index:
                return mer.total_rigid_mass()
        for chain in self.chains:
            if chain.index == center_index:
                return np.sum(chain.mass)

    def moment_inertia(self, center_index):
        for mer in self.mers:
            if mer.index == center_index:
                return mer.moment_inertia()
        for chain in self.chains:
            if chain.index == center_index:
                return chain.moment_inertia()

    def dump_gsd(self, dump_file=None):
        """

        :param dump_file: name of the file to dump the xyz to
        :return: nothing
        """

        filename = dump_file
        if dump_file is None:
            filename = self.name + '.gsd'
        elif dump_file[-4:] != '.gsd':
            filename = dump_file + '.gsd'

        sys = self.create_system()
        #res_map = self.create_res_map(sys)

        hoomd.dump.gsd(filename=filename, period=None, group=hoomd.group.all(), overwrite=True)
        return filename

    def orient_quaternion(self, q):
        temp = self.geometric_center()
        self.shift(np.multiply(-1, temp))
        for chain in self.chains:
            for x in range(chain.num_particles):
                chain.position[x] = q.orient(chain.position[x])
        self.shift(temp)

    def center_at_origin(self):
        temp = self.geometric_center()
        self.shift(np.multiply(-1, temp))

    def create_maps(self, sys):

        bind_site_map = []
        invade_map = []
        types = [t for t in sys.particles.types]

        b_site = []
        invaders = []
        past_rigid = False
        for ind, part in enumerate(sys.particles):
            if part.typeid == types.index('D'):
                b_site.append(ind)
                past_rigid = True
            if part.typeid == types.index('C'):
                invaders.append(ind)
                past_rigid = True
            if part.typeid == types.index('X') and past_rigid:
                bind_site_map.append(b_site)
                invade_map.append(invaders)
                past_rigid = False
                invaders = []
                b_site = []
        bind_site_map.append(b_site)
        invade_map.append(invaders)
        return bind_site_map, invade_map

    def dump_map(self, sys, dump_file=None):

        bind_site_map, invade_map = self.create_maps(sys)
        filename = dump_file
        if dump_file is None:
            filename = self.name + '.map'
        elif dump_file[-4:] != '.map':
            filename = dump_file + '.map'
        f_out = open(filename, 'w')
        for ind, mer in enumerate(invade_map):
            f_out.write('i ')
            for site in mer:
                f_out.write(str(site) + ' ')
            f_out.write('b ')
            for site in bind_site_map[ind]:
                f_out.write(str(site) + ' ')
            f_out.write('\n')
        f_out.close()
        self.dump_map_2(sys, dump_file=dump_file)

    def dump_map_2(self, sys, dump_file=None):

        bind_site_map, invade_map, p_map, x_map, arm_map, temp_map = self.create_maps_2(sys)
        filename = dump_file
        if dump_file is None:
            filename = self.name + '.map2'
        elif dump_file[-5:] != '.map2':
            filename = dump_file + '.map2'
        f_out = open(filename, 'w')
        for ind, mer in enumerate(invade_map):
            f_out.write('i ')
            for site in mer:
                f_out.write(str(site) + ' ')
            f_out.write('b ')
            for site in bind_site_map[ind]:
                f_out.write(str(site) + ' ')
            f_out.write('p ')
            for site in p_map[ind]:
                f_out.write(str(site) + ' ')
            f_out.write('x ')
            for site in x_map[ind]:
                f_out.write(str(site) + ' ')
            f_out.write('a ')
            for site in arm_map[ind]:
                f_out.write(str(site) + ' ')
            f_out.write('\n')
        for mer in temp_map:
            f_out.write('t ')
            for site in mer:
                f_out.write(str(site) + ' ')
            f_out.write('\n')

        f_out.close()

    def create_maps_2(self, sys):

        bind_site_map = []
        invade_map = []
        p_map = []
        x_map = []
        arm_map = []
        temp_map = []

        types = [t for t in sys.particles.types]

        b_site = []
        invaders = []
        p_s = []
        x_s = []
        arm_s = []

        past_rigid = False
        for ind, part in enumerate(sys.particles):
            if part.typeid == types.index('D'):
                b_site.append(ind)
                if not past_rigid:
                    x_map.append(x_s)
                    x_s = []
                past_rigid = True
            if part.typeid == types.index('X'):
                x_s.append(ind)
            if part.typeid == types.index('C'):
                invaders.append(ind)
            if part.typeid == types.index('qPm'):
                temp_map.append(ind)
            if part.typeid == types.index('P') or part.typeid == types.index('H'):
                p_s.append(ind)
            if part.typeid == types.index('qPp'):
                arm_s.append(ind)
            if part.typeid == types.index('X') and past_rigid:
                bind_site_map.append(b_site)
                invade_map.append(invaders)
                p_map.append(p_s)
                arm_map.append(arm_s)
                past_rigid = False
                invaders = []
                b_site = []
                arm_s = []
                p_s = []
        bind_site_map.append(b_site)
        invade_map.append(invaders)
        arm_map.append(arm_s)
        p_map.append(p_s)

        per_chain = int(len(temp_map) / len(self.chains))

        temp_map = [temp_map[per_chain * i: per_chain * (i+1)] for i in range(len(self.chains))]

        return bind_site_map, invade_map, p_map, x_map, arm_map, temp_map



    def get_mers(self,positions, mer):

        num = len(positions)

        news = []

        for ind in range(0, num):
            news.append(cp.deepcopy(mer))

        mer.shift(positions[0])

        for ind, new in enumerate(news):
            new.shift(positions[ind])

        return news



class FiveCoord(Solution):

    def __init__(self, mer, chain, radius=None):

        if radius is None:
            radius = mer.max_radius()

        chain.shift([0, 0, -1 * radius])
        chains = [chain]
        mers = [cp.deepcopy(mer) for _ in range(6)]

        pos = [[radius * np.sin(2 * np.pi/5 * i), radius * np.cos(2 * np.pi/5 * i), 0] for i in range(5)]
        pos.append([0,0,0])
        for ind, m in enumerate(mers):
            m.shift(pos[ind])

        super(FiveCoord, self).__init__(mers, chains)


class Lattice(Solution):

    def __init__(self, mer, chain, n, spacing, flat=False):


        x = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        y = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        z = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        xx, yy, zz = np.meshgrid(x, y, z)
        pos = []
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    pos.append([xx[i][j][k], yy[i][j][k], zz[i][j][k]])

        mers = self.get_mers(pos, mer)


        to_remove = []

        if not flat:
            for ind, mer in enumerate(mers):
                dist = la.norm(np.subtract(mer.body.center, chain.rigid_center_of_mass()))
                if dist < (chain.max_radius() + mer.max_radius()):
                    to_remove.append(ind)
        else:
            for ind, mer in enumerate(mers):
                zloc = np.abs(mer.body.center[2])
                if zloc < 10:
                    to_remove.append(ind)


        for i in range(len(to_remove)):
            ind = len(to_remove) - 1 - i
            mers.remove(mers[to_remove[ind]])

        super(Lattice, self).__init__(mers, [chain], box_length=spacing*n, flat=flat)


class LatticeWall(Lattice):

    def __init__(self, mer, n, spacing):

        from FlatTemplate import FlatTemplate

        chain = FlatTemplate(n * spacing)

        super(LatticeWall, self).__init__(mer, chain, n, spacing, flat=True)

class EightBallLattice(Solution):


    def __init__(self, mer, chain, n, spacing):

        x = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        y = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        z = np.arange(- (n-1) * spacing / 2, (n+2) * spacing/2, spacing)
        xx, yy, zz = np.meshgrid(x, y, z)
        pos = []
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    pos.append([xx[i][j][k], yy[i][j][k], zz[i][j][k]])

        mers = self.get_mers(pos, mer)

        to_remove = []

        length = ((n-1) * spacing)/2

        spot = length/2

        chain_pos = [[i * spot, j * spot, k * spot] for i in range(-1,3,2) for j in range(-1,3,2) for k in range(-1,3,2)]

        chains = []
        for pos in chain_pos:
            newchain = cp.deepcopy(chain)
            newchain.shift(pos)
            chains.append(newchain)

        for ind, mer in enumerate(mers):
            for chane in chains:
                dist = la.norm(np.subtract(mer.body.center, chane.rigid_center_of_mass()))
                if dist < (chane.max_radius() + mer.max_radius()):
                    to_remove.append(ind)


        for i in range(len(to_remove)):
            ind = len(to_remove) - 1 - i
            mers.remove(mers[to_remove[ind]])

        super(EightBallLattice, self).__init__(mers, chains, box_length=spacing*n)


class DoubleSphere(Solution):

    def __init__(self, mer, chain, mer_number, mer_radius):

        sphere_points = np.multiply(self.unit_sphere(mer_number), mer_radius)

        mers= [cp.deepcopy(mer) for _ in range(mer_number)]
        o_list = []

        for ind, new_mer in enumerate(mers):
            align_vec = np.divide(sphere_points[ind], la.norm(sphere_points[ind]))
            q = new_mer.align(align_vec)
            o_list.append(q)
            new_mer.shift(sphere_points[ind])

        chains = [chain]

        super(DoubleSphere, self).__init__(mers, chains, box_length=mer_radius * 4)

        self.o_list = o_list


    def unit_sphere(self, n):
        points = []
        offset = 2. / n
        increment = np.pi * (3. - np.sqrt(5.));

        for i in range(n):
            y = ((i * offset) - 1) + (offset / 2);
            r = np.sqrt(1 - pow(y, 2))

            phi = i * increment

            x = np.cos(phi) * r
            z = np.sin(phi) * r

            points.append([x, y, z])
        return points