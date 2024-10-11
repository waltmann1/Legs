from __future__ import division
from Mer import Mer
from Arm import Arm
from BasicSV40Body import BasicSV40Body
from Leg import Leg
from DNATorroid import DNATorroid
import copy as cp
from Solution import Solution
from Solution import FiveCoord
from Solution import Lattice
from Solution import LatticeWall
from Simulation import Simulation
from Simulation import InitGSD
from ChargedPolymer import ChargedPolymer
from SphericalTemplate import SphericalTemplate
from HaganArm import HaganArm
from SideBody import SideBody
from Leg import Leg
from HelixLeg import HelixLeg

sim = InitGSD('protein_sim.gsd', 1299, total_charge=-500, energy=14)
#sim = InitGSD('protein_sim.gsd', 5900, total_charge=-4800, energy=16)
sim.set_log_period(1000)
sim.set_dump_period(5000)
sim.set_dt(.004)
b = sim.system.box
print(b)
print(sim.nlist)
sim.remove_sphere()
sim.add_polymer(500, angles=0.0)
#sim.run(1000000)
#sim.set_energy(0)
#sim.run(100000)
#sim.remove_sphere()
sim.run(1000000)
