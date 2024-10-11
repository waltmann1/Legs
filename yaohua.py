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
from ChargedPolymer import ChargedPolymer
from SphericalTemplate import SphericalTemplate
from HaganArm import HaganArm
from SideBody import SideBody
from Leg import Leg
from HelixLeg import HelixLeg
from Solution import EightBallLattice
from Solution import DoubleSphere
from PduaBody import PduaBody

#body = BasicSV40Body()
body = PduaBody()

arm = HelixLeg()
#arm = Arm()
leg = HelixLeg()

mer = Mer(body, arm, leg)



temp = SphericalTemplate(10)

chains= [temp]
#mer.align([0,1,0])
#mer.shift([0,15,0])
mers = [mer]

#s = Solution(mers, chains)
#s = FiveCoord(mer, c, radius=30)

#s = LatticeWall(mer, 5, 25)
s = Lattice(mer, temp, 5, 17)
#s = DoubleSphere(mer, temp , 72,45)
os = s.o_list
#s.dump_gsd('wall')
#quit()
sys = s.create_system()

#s.dump_map(sys)

sim = Simulation(sys, name=s.name, total_charge=-16000, o_list=os)
sim.set_log_period(1000)
sim.set_dump_period(5000)
sim.set_dt(.004)
#sim.run(20000)
#sim.set_total_charge(-4800)
sim.run(20000000)
