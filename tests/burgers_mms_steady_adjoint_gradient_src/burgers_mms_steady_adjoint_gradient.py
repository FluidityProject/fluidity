from fluidity_tools import stat_parser as stat
from math import log, isnan

def get_convergence(statfileA, statfileB, field):
  dt_A = stat(statfileA)["ElapsedTime"]['value'][1] - stat(statfileA)["ElapsedTime"]['value'][0]
  dt_B = stat(statfileB)["ElapsedTime"]['value'][1] - stat(statfileB)["ElapsedTime"]['value'][0]

  a_error_l1 = sum(stat(statfileA)["Fluid"][field]["integral"])*dt_A
  b_error_l1 = sum(stat(statfileB)["Fluid"][field]["integral"])*dt_B

  a_error_l2 = sum([x**2*dt_A for x in stat(statfileA)["Fluid"][field]["l2norm"]])**0.5
  b_error_l2 = sum([x**2*dt_B for x in stat(statfileB)["Fluid"][field]["l2norm"]])**0.5

  a_error_inf = max(stat(statfileA)["Fluid"][field]["max"])
  b_error_inf = max(stat(statfileB)["Fluid"][field]["max"])

  # Velocity error calculation
  ab_ratio_l1 = a_error_l1 / b_error_l1
  ab_ratio_l2 = a_error_l2 / b_error_l2
  ab_ratio_inf = a_error_inf / b_error_inf

  ab_error = [log(ab_ratio_l1, 2), log(ab_ratio_l2, 2), log(ab_ratio_inf, 2)]
  return ab_error

def test_convergence(statfiles, fields, tol):
  for field in fields:
    for i in range(len(statfiles)):
      if i==0:
        continue
      c = min(get_convergence(statfiles[i-1], statfiles[i], field))
      if isnan(c) or c < tol:
          return False
  return True      

def test_convergence_rates(tol):
  fields = ["Error"]
  statfiles = ["mms_adjoint_B_forward.stat", "mms_adjoint_C_forward.stat", "mms_adjoint_D_forward.stat", "mms_adjoint_E_forward.stat"]
  res1 = test_convergence(statfiles, fields, tol)
  return res1

def print_convergence(statfiles, fields):
  for field in fields:
    print "Field: ", field
    for i in range(len(statfiles)):
      if i==0:
        continue
      print statfiles[i-1], " : ", statfiles[i], ' ', get_convergence(statfiles[i-1], statfiles[i], field)

def print_convergence_rates():
  print ""
  print " =============== Time and space resolution increasing ================="
  fields = ["Error"]
  statfiles = ["mms_adjoint_A_forward.stat", "mms_adjoint_B_forward.stat", "mms_adjoint_C_forward.stat", "mms_adjoint_D_forward.stat", "mms_adjoint_E_forward.stat"]
  print_convergence(statfiles, fields)
