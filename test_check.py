from  __future__ import division
import sys
sys.path.append('/home/cwj8781/Legs')

from Mer import Mer
from Arm import Arm
from BasicSV40Body import BasicSV40Body
from Leg import Leg
from DNATorroid import DNATorroid
import copy as cp
from Solution import Solution
from Solution import FiveCoord
from Solution import Lattice
from Simulation import Simulation
from ChargedPolymer import ChargedPolymer
from SphericalTemplate import SphericalTemplate
from HaganArm import HaganArm
from HelixLeg import HelixLeg
from Solution import DoubleSphere

body = BasicSV40Body()
arm = HaganArm()
leg = HelixLeg()
mer = Mer(body, arm, leg)

N = 13

r = 14
guess = 7


tc = -25 * 2 * N * (r/guess)**2


temp = SphericalTemplate(r)

s = DoubleSphere(mer, temp,N, r+4)
sys = s.create_system()
s.dump_map(sys)

os=s.o_list

sim = Simulation(sys, energy=8, helix=0, total_charge=-427, o_list=os,effective_charge=.01,debye=.5, counter=True)
sim.set_log_period(10000)
sim.set_dump_period(5e4)
sim.set_dt(.004)
sim.run(5e6)
print(r*2, r*2-1)
for i in range(r*2, (guess-3) * 2 ,-1):
        sim.rescale_sphere_density(factor=(i-1)/i )
        sim.run(5e6)



