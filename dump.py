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

body = BasicSV40Body()
#body = SideBody()

arm = HaganArm()
#arm = Arm()
leg = HelixLeg()

mer = Mer(body, arm, leg)

mer.align([1,1,0])
mer.dump_xyz("h.xyz")
