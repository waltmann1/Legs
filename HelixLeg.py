from __future__ import division
import numpy as np
from Leg import Leg


class HelixLeg(Leg):

    def __init__(self, index=-1):

        super(HelixLeg, self).__init__(index=index)
        for i in range(0, self.length):
            if i > 1 and i < self.length -3:
                self.type[i] = 'H'
