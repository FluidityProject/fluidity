
import solution
from fluidity_tools import stat_parser as stat
from math import log, sqrt
from scipy.integrate import quad

def report_convergence(file1, file2):
  print(file1, "->", file2)
  
  stat1 = stat(file1)
  stat2 = stat(file2)

  print(stat1["dt"]["value"][0], "->", stat2["dt"]["value"][0])
  
  errortop_l2_1 = sqrt(sum(stat1["Fluid"]["DifferenceSquared"]["surface_integral%TopSurfaceL2Norm"][1:]*stat1["dt"]["value"][1:]))
  errortop_l2_2 = sqrt(sum(stat2["Fluid"]["DifferenceSquared"]["surface_integral%TopSurfaceL2Norm"][1:]*stat2["dt"]["value"][1:]))
  convergencetop_l2 = log((errortop_l2_1/errortop_l2_2), 2)

  print('  convergencetop_l2 = ', convergencetop_l2)
  print('    errortop_l2_1 = ', errortop_l2_1)
  print('    errortop_l2_2 = ', errortop_l2_2)
  
  errorbottom_l2_1 = sqrt(sum(stat1["Fluid"]["DifferenceSquared"]["surface_integral%BottomSurfaceL2Norm"][1:]*stat1["dt"]["value"][1:]))
  errorbottom_l2_2 = sqrt(sum(stat2["Fluid"]["DifferenceSquared"]["surface_integral%BottomSurfaceL2Norm"][1:]*stat2["dt"]["value"][1:]))
  convergencebottom_l2 = log((errorbottom_l2_1/errorbottom_l2_2), 2)

  print('  convergencebottom_l2 = ', convergencebottom_l2)
  print('    errorbottom_l2_1 = ', errorbottom_l2_1)
  print('    errorbottom_l2_2 = ', errorbottom_l2_2)
  
  error_l2_1 = sqrt(sum(stat1["Fluid"]["DifferenceSquared"]["surface_integral%SurfaceL2Norm"][1:]*stat1["dt"]["value"][1:]))
  error_l2_2 = sqrt(sum(stat2["Fluid"]["DifferenceSquared"]["surface_integral%SurfaceL2Norm"][1:]*stat2["dt"]["value"][1:]))
  convergence_l2 = log((error_l2_1/error_l2_2), 2)

  print('  convergence_l2 = ', convergence_l2)
  print('    error_l2_1 = ', error_l2_1)
  print('    error_l2_2 = ', error_l2_2)
  
  error_linf_1 = stat1["Fluid"]["FreeSurfaceDifference"]["max"].max()
  error_linf_2 = stat2["Fluid"]["FreeSurfaceDifference"]["max"].max()
  convergence_linf = log((error_linf_1/error_linf_2), 2)

  print('  convergence_linf = ', convergence_linf)
  print('    error_linf_1 = ', error_linf_1)
  print('    error_linf_2 = ', error_linf_2)

  quad1 = quad(lambda t: solution.nond_error_amp(stat1, t)**2, stat1["ElapsedTime"]["value"][0], stat1["ElapsedTime"]["value"][-1], limit=1000)
  quad2 = quad(lambda t: solution.nond_error_amp(stat2, t)**2, stat2["ElapsedTime"]["value"][0], stat2["ElapsedTime"]["value"][-1], limit=1000)
  errormaxfs_l2_1 = sqrt(quad1[0])
  errormaxfs_l2_2 = sqrt(quad2[0])
  convergencemaxfs_l2 = log((errormaxfs_l2_1/errormaxfs_l2_2), 2)

  print('  convergencemaxfs_l2 = ', convergencemaxfs_l2)
  print('    errormaxfs_l2_1 = ', errormaxfs_l2_1, '(', quad1[1], ')')
  print('    errormaxfs_l2_2 = ', errormaxfs_l2_2, '(', quad2[1], ')')

  return [convergencetop_l2, convergencebottom_l2, convergence_linf]
  
