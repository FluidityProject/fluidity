import math
import os
import re
import sys
from math import sqrt

import vtktools
from numpy import arange
from numpy import argsort
from numpy import array
from numpy import concatenate


def values_per_node(file):

    u = vtktools.vtu(file)
    zoo = u.GetScalarField("Zooplankton")
    phyto = u.GetScalarField("Phytoplankton")
    nut = u.GetScalarField("Nutrient")
    det = u.GetScalarField("Detritus")

    return phyto, zoo, nut, det
