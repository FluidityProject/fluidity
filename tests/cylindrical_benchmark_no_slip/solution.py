import assess

nu = 1.0
n = 2
solution_upper = assess.CylindricalStokesSolutionDeltaZeroSlip(n, +1, nu=nu)
solution_lower = assess.CylindricalStokesSolutionDeltaZeroSlip(n, -1, nu=nu)
