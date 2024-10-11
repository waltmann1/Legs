from __future__ import division
from Bonds import Bonds
from NonBonded import Yukawa
from NonBonded import LJSpecial
from NonBonded import LJRepulsive
from RigidBodies import Rigid
from Angles import Angles
import numpy as np
from numpy import linalg as la
import hoomd

class Simulation(object):

    def __init__(self, system, temperature=1, name="protein_sim", energy=None, total_charge=None, helix=None
                 , reinit=None, o_list=None, effective_charge=1, counter=False, debye=1, depletion=False, attract=None):
        self.system = system
        self.nlist = hoomd.md.nlist.cell(check_period=1)
        self.nlist.reset_exclusions(exclusions=['bond', 'angle', 'dihedral', 'constraint', 'body'])
        #self.log_list = ['potential_energy', 'temperature', 'kinetic_energy']
        self.log_list = ['potential_energy', 'ndof', 'kinetic_energy']
        self.log_list.append('temperature')
        self.log_period = 1000
        self.dump_period = 10000
        self.temperature = temperature
        self.name = name
        self.total_charge = total_charge
        self.effective_charge = effective_charge
        self.debye = debye
        self.counter = counter
        self.helix = helix
        self.depletion = depletion

        if o_list is not None:
            for i in range(len(o_list)):
                self.system.particles[i].orientation = o_list[i]

        self.dt = .004

        self.rigid = Rigid()
        self.rigid.set_rigid_bodies(system, reinit=reinit, o_list=o_list)

        self.bonds = Bonds(self.log_list)
        self.bonds.set_all_harmonic_bonds(system)

        self.angles = Angles(self.log_list)
        self.angles.set_all_harmonic_angles(system)


        self.dt = .004

        table_cd = True
        self.table_cd = table_cd
        if energy is not None:
            table_cd = False
            self.ljs = LJSpecial(self.log_list, energy=energy)
            self.ljs.set_lj_special(self.nlist, system, depletion=depletion)

        self.ljr = LJRepulsive(self.log_list)
        if counter:
            self.ljr.set_lj_repulsive(self.nlist, system,table_cd=table_cd, helix=helix,
                                    effective_charge=self.effective_charge, debye=self.debye,
                                    sphere=self.rigid.get_spherical_radius(system), depletion=depletion, attract=attract)
        else:
            self.ljr.set_lj_repulsive(self.nlist, system, table_cd=table_cd, helix=helix, debye=None,
                                      sphere=self.rigid.get_spherical_radius(system), depletion=depletion, attract=attract)

        self.yukawa = Yukawa(log_list=self.log_list, debye=self.debye, total_charge=total_charge, effective_charge=self.effective_charge)
        self.yukawa.set_yukawa(self.nlist, system)

        self.all = hoomd.group.all()

        self.to_integrate = hoomd.group.union(name='dof', a=hoomd.group.rigid_center(), b=hoomd.group.nonrigid())


        #if 'centerflat' in system.particles.types:
        #    group_flat = hoomd.group.type(name='group_flat', type='centerflat')
        #    self.to_integrate = hoomd.group.difference(name="", a=self.to_integrate, b=group_flat)

        if not self.check_mobile():
            self.remove_rigid_template_from_integration()

        hoomd.md.integrate.mode_standard(dt=self.dt)
        self.nve = hoomd.md.integrate.nve(group=self.to_integrate, limit=.001)
        self.nve.disable()
        self.langevin = hoomd.md.integrate.langevin(group=self.to_integrate, kT=self.temperature, seed=42)

        log_name = self.name + ".log"
        self.logger = hoomd.analyze.log(filename=log_name, quantities=self.log_list, period=self.log_period,
                          overwrite=True)

        dump_name = self.name + ".gsd"
        self.dumper = hoomd.dump.gsd(filename=dump_name, period=self.dump_period, group=self.all, overwrite=True)

    def run(self, time):

        #print(self.system.constraints)
        hoomd.run(time)

    def run_nanoseconds(self, time):

        real_time = int(time * 1e-9 / (self.time_unit * self.dt))
        self.run(real_time)

    def nve_relaxation(self, time):


         self.langevin.disable()
         self.nve.enable()

         hoomd.run(time)
         self.nve.set_params(limit=.01)
         #hoomd.run(time / 2)
         #self.nve.set_params(limit=.001)
         #self.nve.set_params(limit=.01)
         #hoomd.run(time)
         #self.nve.set_params(limit=.1)
         #hoomd.run(time)
         #self.nve.set_params(limit=1)
         self.nve.disable()
         self.langevin.enable()

    def set_dt(self, dt):
        hoomd.md.integrate.mode_standard(dt=dt)

    def run_fire(self, time):

        self.langevin.disable()
        self.nve.enable()
        fire = hoomd.md.integrate.mode_minimize_fire(dt=0.1, group=self.to_integrate, ftol=1e-2, Etol=1e-7)
        hoomd.run(time)
        del fire
        self.langevin.enable()
        self.nve.disable()
        hoomd.md.integrate.mode_standard()

    def temp_interp(self, temp1, temp2, time):

        t1 = temp1
        t2 = temp2
        self.langevin.set_params(kT=hoomd.variant.linear_interp(points=[(0, t1), (time, t2)]))
        hoomd.run(time)
        self.langevin.set_params(kT=self.temperature)

    def set_temperature(self, t):
        temp = 0.596 / 300 * t
        self.temperature = temp
        self.langevin.set_params(kT=self.temperature)

    def palace_equil(self):
        #self.run_fire(5500)
        self.temp_interp(0, 300, 100000)
        self.run(200000)

    def basic_temp_equil_no_log(self):

        self.logger.disable()
        self.dumper.disable()
        self.set_temperature(0)
        self.run(10000)
        self.temp_interp(0, 300, 100000)
        self.set_temperature(300)
        self.run(10000)
        self.logger.enable()
        self.dumper.enable()

    def set_log_period(self, period):

        self.logger.disable()
        self.log_period = period
        log_name = self.name + ".log"
        self.logger = hoomd.analyze.log(filename=log_name, quantities=self.log_list, period=self.log_period,
                                        overwrite=True)

    def set_dump_period(self, period):

        self.dumper.disable()
        self.dump_period = period
        dump_name = self.name + ".gsd"
        self.dumper = hoomd.dump.gsd(filename=dump_name, period=self.dump_period, group=self.all, overwrite=True)


    def total_kinetic_energy(self):

        ke = 0
        for part in self.system.particles:
            kin = .5 * part.mass * np.linalg.norm(part.velocity) ** 2
            print(part.type, kin)
            ke += kin


        return ke

    def ndof(self):

        return self.total_kinetic_energy() * 2 / self.temperature


    def set_total_charge(self, total_charge):

        self.yukawa.yukawa.disable()
        #del self.yukawa
        self.yukawa = Yukawa(log_list=self.log_list, debye=self.debye,  total_charge=total_charge, effective_charge= self.effective_charge)
        self.yukawa.set_yukawa(self.nlist, self.system)
        self.total_charge = total_charge

    def set_effective_charge(self, effective_charge):

        self.yukawa.yukawa.disable()
        #del self.yukawa
        self.yukawa = Yukawa(log_list=self.log_list, debye=self.debye, total_charge=self.total_charge,
                             effective_charge=effective_charge)
        self.yukawa.set_yukawa(self.nlist, self.system)
        self.effective_charge = effective_charge

    def set_debye_length(self, debye):

        self.debye = debye
        self.yukawa.yukawa.disable()
        self.yukawa = Yukawa(log_list=self.log_list, debye=self.debye,  total_charge=self.total_charge, effective_charge=self.effective_charge)
        self.yukawa.set_yukawa(self.nlist, self.system)

    def set_energy(self, energy):
        self.ljs.lj_pair.disable()
        #del self.ljs
        self.ljs = LJSpecial(self.log_list, energy=energy)
        self.ljs.set_lj_special(self.nlist, self.system)

    def rescale_sphere(self, factor=.99):

        center_pos = []
        for part in self.system.particles:
            name = part.type
            if part.type[:6] == "center":
                center_pos = part.position
            elif name == 'qPm':
                part.position = np.add(center_pos, np.multiply(np.subtract(part.position, center_pos), factor))
        self.rigid.set_rigid_bodies(self.system)
        if self.counter:
            self.ljr.lj_repulsive_pair.disable()
            self.ljr = LJRepulsive(self.log_list)
            self.ljr.set_lj_repulsive(self.nlist, self.system, table_cd=self.table_cd, helix=self.helix,
                                          effective_charge=self.effective_charge, debye=self.debye,
                                          sphere=self.rigid.get_spherical_radius(self.system), depletion=self.depletion)


    def remove_rigid_template_from_integration(self):

        for part in self.system.particles:
            name = part.type
            if name[:6] == 'center' and len(name) > 6 and part.body > -1:
                group_gone = hoomd.group.type(name=name[:6], type=name)
                self.to_integrate = hoomd.group.difference(name="", a=self.to_integrate, b=group_gone)

    def check_mobile(self):
        center_number=0
        for part in self.system.particles:
            if part.type[:6] == 'center' and len(part.type)>6:
                center_number += 1
        if center_number == 1:
             return False

        return True

    def remove_helices(self):
        self.bonds.relax_bonds()

    def rescale_sphere_remove(self, factor=.99):

        from SphericalTemplate import SphericalTemplate

        removed =0
        total = 0
        center_pos = []
        center_ind = 0
        to_remove = []
        radius = 0
        temp_ind = 0
        temps = False
        for ind, part in enumerate(self.system.particles):
            #print(ind, part.type)
            name = part.type
            if part.type[:6] == "center":
                center_pos = part.position
            elif name == 'qPm' or temps:
                temps = True
                if radius == 0:
                    radius = la.norm(np.subtract(part.position, center_pos))
                    new_radius = factor * radius
                    print(factor, new_radius, radius)
                    #quit()
                    temp = SphericalTemplate(new_radius)
                if temp_ind < len(temp.position):
                    part.position = np.add(center_pos, temp.position[temp_ind])
                    temp_ind += 1
                else:
                    #print(ind, "hey its time to change")
                    part.type = "P"
                    part.position = np.add(center_pos, np.multiply(np.subtract(part.position, center_pos), factor))
                    #to_remove.append(part.tag)
                    removed += 1
                total += 1
        #quit()
        #l = len(to_remove)
        #to_remove = [to_remove[i] for i in range(l-1,-1,-1)]
        #print(to_remove)
        #for ind in range(len(to_remove)):
        #    print(ind)
           # self.system.particles.remove(to_remove[ind])
        self.rigid.set_rigid_bodies(self.system, reset=True)
        if self.counter:
            self.ljr.lj_repulsive_pair.disable()
            self.ljr = LJRepulsive(self.log_list)
            self.ljr.set_lj_repulsive(self.nlist, self.system, table_cd=self.table_cd, helix=self.helix,
                                          effective_charge=self.effective_charge, debye=self.debye,
                                          sphere=self.rigid.get_spherical_radius(self.system),
                                      depletion=self.depletion)
        #self.set_total_charge(self.total_charge * removed/total)


    def rescale_sphere_density(self, factor=.99):

        from SphericalTemplate import SphericalTemplate

        removed = 0
        total = 0
        center_pos = []
        center_ind = 0
        to_remove = []
        radius = 0
        temp_ind = 0
        temps = False
        for ind, part in enumerate(self.system.particles):
            # print(ind, part.type)
            name = part.type
            if part.type[:6] == "center":
                center_pos = part.position
            elif name == 'qPm' or temps:
                temps = True
                if radius == 0:
                    radius = la.norm(np.subtract(part.position, center_pos))
                    new_radius = factor * radius
                    print(factor, new_radius, radius)
                    num = len(self.system.particles) - ind
                    # quit()
                    temp = SphericalTemplate(1)
                    new_points = np.multiply(temp.unit_sphere(num), new_radius)
                    #print("new points", new_points)
                if temp_ind < len(new_points):
                    part.position = np.add(center_pos, new_points[temp_ind])
                    temp_ind += 1
                else:
                    print(ind, "hey its time to change")
                    part.type = "P"
               #    part.position = np.add(center_pos, np.multiply(np.subtract(part.position, center_pos), factor))
               #    to_remove.append(part.tag)
            #    removed += 1
                total += 1
        # quit()
        # l = len(to_remove)
        # to_remove = [to_remove[i] for i in range(l-1,-1,-1)]
        # print(to_remove)
        # for ind in range(len(to_remove)):
        #    print(ind)
        # self.system.particles.remove(to_remove[ind])
        self.rigid.set_rigid_bodies(self.system, reset=True)
        self.set_total_charge(self.total_charge * (factor)**2)
        if self.counter:
            self.ljr.lj_repulsive_pair.disable()
            self.ljr = LJRepulsive(self.log_list)
            self.ljr.set_lj_repulsive(self.nlist, self.system, table_cd=self.table_cd, helix=self.helix,
                                          effective_charge=self.effective_charge, debye=self.debye,
                                          sphere=self.rigid.get_spherical_radius(self.system),
                                      depletion=self.depletion)


    def remove_sphere(self, index=12):

        to_remove = []

        for ind, part in enumerate(self.system.particles):
            name = part.type
            if name == 'qPm':
                to_remove.append(part.tag)
        l = len(to_remove)
        to_remove = [to_remove[i] for i in range(l-1, -1, -1)]
        #could use a better system to remove the center
        self.system.particles[index].type = "P"
        self.system.particles[index].body = -1
        spot = self.system.box.Lx / 2 - 1
        self.system.particles[index].position = [spot, spot, spot]
        for x in range(l):
            self.system.particles.remove(to_remove[x])

        self.remove_rigid_template_from_integration()
        self.rigid.set_rigid_bodies(self.system, reset=True)
        self.set_total_charge(0)

        self.to_integrate = hoomd.group.union(name='dof', a=hoomd.group.rigid_center(), b=hoomd.group.nonrigid())
        pgone = hoomd.group.tag_list(name='pgone', tags=[index])
        self.to_integrate = hoomd.group.difference(name='dof', a=self.to_integrate, b=pgone)
        self.langevin.disable()
        hoomd.md.integrate.mode_standard(dt=self.dt)
        self.nve = hoomd.md.integrate.nve(group=self.to_integrate, limit=.001)
        self.nve.disable()
        self.langevin = hoomd.md.integrate.langevin(group=self.to_integrate, kT=self.temperature, seed=42)


    def add_polymer(self, length, z_off=0, angles=None):

        from ChargedPolymer import ChargedPolymer

        chain = ChargedPolymer(length, angles=angles)


        snap = self.system.take_snapshot(all=True)
        n_original = snap.particles.N
        snap.particles.resize(n_original + length)
        snap.bonds.types = snap.bonds.types + [chain.bond_names[0]]
        if angles is not None:
            snap.angles.types = snap.angles.types + [chain.angle_names[0]]
        b_types = snap.bonds.types
        p_types = snap.particles.types
        a_types = snap.angles.types

        tag = n_original
        for x in range(len(chain.bonds)):
            bond_number = snap.bonds.N + 1
            snap.bonds.resize(bond_number)
            snap.bonds.group[bond_number - 1] = np.add(chain.bonds[x], tag)
            snap.bonds.typeid[bond_number - 1] = b_types.index(chain.bond_names[x])

        if angles is not None:
            for x in range(len(chain.angles)):
                angle_number = snap.angles.N + 1
                snap.angles.resize(angle_number)
                snap.angles.group[angle_number - 1] = np.add(chain.angles[x], tag)
                snap.angles.typeid[angle_number - 1] = a_types.index(chain.angle_names[x])

        for x in range(chain.length):
            snap.particles.position[x + tag] = np.add(chain.position[x], [0, 0, z_off])
            snap.particles.mass[x + tag] = chain.mass[x]
            snap.particles.typeid[x + tag] = p_types.index(chain.type[x])
            snap.particles.body[x + tag] = chain.index
            # snap.particles.charge[x + tag] = chain.charge[x]
            snap.particles.charge[x + tag] = 0
        tag += chain.length

        self.system.restore_snapshot(snap)
        self.set_total_charge(self.total_charge)
        self.bonds.set_all_harmonic_bonds(self.system, reset=True)
        self.angles.set_all_harmonic_angles(self.system, reset=True, poly=angles)

        poly = hoomd.group.tag_list(name='poly', tags=list(range(n_original, snap.particles.N)))
        self.to_integrate = hoomd.group.union(name='dof', a=self.to_integrate, b=poly)
        self.langevin.disable()
        hoomd.md.integrate.mode_standard(dt=self.dt)
        self.nve = hoomd.md.integrate.nve(group=self.to_integrate, limit=.001)
        self.nve.disable()
        self.langevin = hoomd.md.integrate.langevin(group=self.to_integrate, kT=self.temperature, seed=42)





"""
    def set_rigid_orientations(self):

        for part in self.system.particles:
            if part.body > -1:
                part.orientation = self.system.particles[part.body].orientation
"""
class InitGSD(Simulation):

    def __init__(self, name, frame, energy=None, total_charge=None, helix=False, effective_charge=1,
                 counter=False, debye=1, depletion=False):

        hoomd.context.initialize("--mode=gpu")
        system = hoomd.init.read_gsd(name, frame=frame)

        i = 0

        while not name[i].isalpha():
            i = i + 1

        name_no_loc = name[i:]

        super(InitGSD, self).__init__(system, name=name_no_loc[:-4] + '_frame' + str(frame),
                                      energy=energy,  total_charge=total_charge, helix=helix,
                                      effective_charge=effective_charge, counter=counter, debye=debye,
                                      depletion=depletion)
