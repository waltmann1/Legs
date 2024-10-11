from Mer import Mer
from Arm import Arm
from MassAdjustedSV40Body import MassAdjustedSV40Body
from Leg import Leg
from DNATorroid import DNATorroid
import copy as cp
from Solution import Solution
from Solution import FiveCoord
from Solution import Lattice
from Solution import LatticeWall
from Simulation import Simulation
from ChargedPolymer import ChargedPolymer
from SphericalTemplate import SphericalTemplate
from HaganArm import HaganArm
from SideBody import SideBody
from Leg import Leg
from HelixLeg import HelixLeg
from Solution import EightBallLattice
from Solution import DoubleSphere
from PduaBody import PduaBody
from ReducedArm import ReducedArm
from BasicSV40Body import BasicSV40Body


body = BasicSV40Body()
#body = PduaBody()

#body = MassAdjustedSV40Body()

arm = HaganArm()
#arm = ReducedArm()
#arm = Arm()
leg = HelixLeg()

mer = Mer(body, arm, leg)
guess= 7
r = 14
N = 2
#tc = -N * 25 * 2 * (r/guess)**2

temp = SphericalTemplate(r)

chains= [temp]
#mer.align([0,1,0])
#mer.shift([0,15,0])
mers = [mer]

#s = Solution(mers, chains)
#s = FiveCoord(mer, c, radius=30)

#s = LatticeWall(mer, 5, 25)
#s = EightBallLattice(mer, temp, 5, 25)
s = DoubleSphere(mer, temp , N,r+4)
os = s.o_list
#s.dump_gsd('wall')
#quit()
sys = s.create_system()

#s.dump_map(sys)
sim = Simulation(sys, name=s.name, total_charge=-427, energy=8, effective_charge=1, debye=.5, counter=True)
sim.set_log_period(10000)
sim.set_dump_period(50000)
sim.set_dt(.004)
#sim.run(100000)
#sim.set_energy(0)
#sim.run(200000)
#sim.rescale_sphere(factor=.1)
#sim.remove_sphere()
sim.run(100000000)
