import assess

nu = 1.0
n = 2
solution_upper = assess.CylindricalStokesSolutionDeltaFreeSlip(n, +1, nu=nu)
solution_lower = assess.CylindricalStokesSolutionDeltaFreeSlip(n, -1, nu=nu)
