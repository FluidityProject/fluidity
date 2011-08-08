import os
import sys
import math

### Model Constants ###
param_A_E =            -10000
param_Alpha_Chl =      7.9E-7
param_C_minS =         1.584E-8
param_C_rep =          1.76E-8
param_C_starve =       8.5E-9
param_k_AR =           1.0
param_k_S =            1.0
param_Ndis =           0.0042
param_P_ref_c =        0.14
param_Q_Nmax =         0.17
param_Q_Nmin =         0.034
param_Q_remN =         2.95
param_Q_remS =         2.27
param_Q_S_max =        0.15
param_Q_S_min =        0.04
param_R_Chl =          2E-3
param_R_maintenance =  2E-3
param_R_N =            2E-3
param_S_dis =          8.3E-4
param_S_rep =          2.1E-9
param_T_ref =          293.0
param_T_refN =         283.0
param_T_refS =         278.0
param_Theta_max_N =    4.2
param_V_ref_c =        0.01
param_V_S_ref =        0.03
param_z_sink =         0.04
param_Zeta =           2.3

def photosynthesis(vars, name, ambientTemperature, ambientVisIrrad):
  """Photsynthesis according to Sinerchia ..."""

  T_K = ambientTemperature + 273.0
  T_function = math.exp(34.12969283 - 10000/T_K)

  Q_N = ((((vars[name['AmmoniumPool']])+(vars[name['IngAmmonium']])+(vars[name['NitratePool']])+(vars[name['IngNitrate']]))) / (vars[name['CarbonPool']]))
  Q_s = ((((vars[name['SilicatePool']])+(vars[name['IngSilicate']]))) / (vars[name['CarbonPool']]))

  #P_max_c = (((Q_N) > (param_Q_Nmax))?(((param_P_ref_c)*(T_function))):(((Q_N) < (param_Q_Nmin))?(0.0):(((param_P_ref_c)*(T_function)*(((Q_N) - (param_Q_Nmin)) / ((param_Q_Nmax) - (param_Q_Nmin)))))))
  if ((Q_N) > (param_Q_Nmax)):
    P_max_c = (((param_P_ref_c)*(T_function)))
  else:
    if ((Q_N) < (param_Q_Nmin)):
      P_max_c = 0.0
    else:
      P_max_c = (((param_P_ref_c)*(T_function)*(((Q_N) - (param_Q_Nmin)) / ((param_Q_Nmax) - (param_Q_Nmin)))))

  Theta_c = ((vars[name['ChlorophyllPool']]) / (vars[name['CarbonPool']]))
  E_0 = (((4.6)*(ambientVisIrrad)))

  #P_phot_c = ((P_max_c == 0.0)||(Q_s <= param_Q_S_min))?0.0:P_max_c*(1.0 - exp(-0.284400e-2*Theta_c*E_0 / P_max_c))
  if ((P_max_c == 0.0)or(Q_s <= param_Q_S_min)):
    P_phot_c = 0.0
  else:
    P_phot_c = P_max_c*(1.0 - math.exp(-0.284400e-2*Theta_c*E_0 / P_max_c))

  return (P_phot_c, Theta_c, E_0)

def reproduction(vars, name, P_phot_c, Theta_c, E_0, stepInHours):

  Theta_N = ((vars[name['ChlorophyllPool']]) / (((vars[name['AmmoniumPool']])+(vars[name['IngAmmonium']])+(vars[name['NitratePool']])+(vars[name['IngNitrate']]))))

  #Rho_Chl = (((((E_0) > (0.0))&&((Theta_c) > (0.0))))?(((param_Theta_max_N)*((P_phot_c) / (((3600.0)*(param_Alpha_Chl)*(Theta_c)*(E_0)))))):(0.0))
  if ((((E_0) > (0.0))and((Theta_c) > (0.0)))):
    Rho_Chl = (((param_Theta_max_N)*((P_phot_c) / (((3600.0)*(param_Alpha_Chl)*(Theta_c)*(E_0))))))
  else:
    Rho_Chl = 0.0

  R_C_growth = ((((((vars[name['IngAmmonium']])+(vars[name['IngNitrate']])))*(param_Zeta))) / (((stepInHours)*(vars[name['CarbonPool']]))))
  R_C = (((param_R_maintenance)+(R_C_growth)));

  #C_d = (((((((vars[name['CarbonPool']])+(((vars[name['CarbonPool']])*((P_phot_c) - (R_C))*(stepInHours))))) >= (param_C_rep))&&((((vars[name['SilicatePool']])+(vars[name['IngSilicate']]))) >= (param_S_rep))))?(2.0):(1.0))
  if ((((((vars[name['CarbonPool']])+(((vars[name['CarbonPool']])*((P_phot_c) - (R_C))*(stepInHours))))) >= (param_C_rep))and((((vars[name['SilicatePool']])+(vars[name['IngSilicate']]))) >= (param_S_rep)))):
    C_d = 2.0
  else:
    C_d = 1.0

  return C_d


#####################################
## Utiliy functions for Hyperlight ##
#####################################

def read_chl_h42(filename):
  """This function reads the chlorophyll concentration from file in 
     Hydrolight standard format (10 header lines, followed by depth/chl pairs)
  """

  try:
    os.stat(filename)
  except:
    print "No such file: " + str(filename)
    sys.exit(1)
  f = open(filename, 'r')
  for i in range(0,10):
    # skip header lines
    line = f.readline()
  
  depths = []
  chl = []
  while line:
    line = f.readline().strip()
    if line == "":
      break
    depth = line.split()[0]
    if depth < 0.0:
      break
    depths.append(float(depth))
    chl.append(float(line.split()[1]))

  f.close()
  return (depths, chl)


def derive_PAR_irradiance(state):
  """Solar irradiance in PAR is derived from 36 individual wavebands as modelled by Hyperlight
     The spectral bands are in Wm-2nm-1, and the PAR total is in mumol phot m-2s-1
  """

  planck = 6.626E-34;
  speed = 2.998E8;
  fnmtom = 1.0E-9;
  fmoltoumol = 1.0E6;
  avagadro = 6.023E23;
  factor = fnmtom * fmoltoumol / (planck*speed*avagadro)

  par_irrad = state.scalar_fields['IrradiancePAR']
  for n in range(par_irrad.node_count):
    irrad_sum = 0.0
    for l in range(0, 35):
      wavelength = 350+l*10
      field_name = 'Irradiance_' + str(wavelength)
      spectral_irrad = state.scalar_fields[field_name]
      irrad_sum = irrad_sum + spectral_irrad.node_val(n) * factor * wavelength * 10.0;

    par_irrad.set(n, irrad_sum)
  return 0.0
