import assess
from assess import Y_cartesian

nu = 1.0
l = 2
m = l
rp = 1.72
solution_upper = assess.SphericalStokesSolutionDeltaFreeSlip(l, m, +1, nu=nu, rp=rp)
solution_lower = assess.SphericalStokesSolutionDeltaFreeSlip(l, m, -1, nu=nu, rp=rp)
