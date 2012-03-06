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

def update_living_diatom(vars, env, dt):
  """ Update kernel for living Diatom agents from LERM-ES, 
      created by M. Sinerchia nd W. Hinsley

     References:
     "Testing theories of fisheries recruitment.", Sinerchia, M., 2007
     PhD Thesis, Department of Earth Science and Engineering, Imperial College, London

     "Using an individual-based model with four trophic levels to model the effect of predation and competition on squid recruitment", 
     Sinerchia, M., Field, A. J., Woods, J. D., Vallerga, S., and Hinsley, W. R.
     ICES Journal of Marine Science, doi:10.1093/icesjms/fsr190.
  """

  # Housekeeping
  c_old = vars['Size']
  vars['AmmoniumUptake'] = vars['AmmoniumUptake'] / c_old
  vars['NitrateUptake'] = vars['NitrateUptake'] / c_old
  vars['SilicateUptake'] = vars['SilicateUptake'] / c_old

  stepInHours = dt/3600.

  # Temperature conversion
  T_K = env['Temperature'] + 273.0
  T_function = math.exp(34.12969283 - 10000/T_K)

  # Photosynthesis
  Q_N = (vars['Ammonium'] + vars['AmmoniumUptake'] + vars['Nitrate'] + vars['NitrateUptake']) / vars['Carbon']
  Q_s = (vars['Silicate'] + vars['SilicateUptake']) / vars['Carbon']

  if Q_N > param_Q_Nmax:
    P_max_c = param_P_ref_c * T_function
  else:
    if Q_N < param_Q_Nmin:
      P_max_c = 0.0
    else:
      P_max_c = param_P_ref_c * T_function * (Q_N - param_Q_Nmin) / (param_Q_Nmax - param_Q_Nmin)

  Theta_c = vars['Chlorophyll'] / vars['Carbon']
  E_0 = 4.6 * env['Irradiance']

  if P_max_c == 0.0 or Q_s <= param_Q_S_min:
    P_phot_c = 0.0
  else:
    P_phot_c = P_max_c * (1.0 - math.exp(-0.284400e-2*Theta_c*E_0 / P_max_c))

  # Chlorophyll Synthesis
  Theta_N = vars['Chlorophyll'] / (vars['Ammonium'] + vars['AmmoniumUptake'] + vars['Nitrate'] + vars['NitrateUptake'])

  if E_0 > 0.0 and Theta_c > 0.0:
    Rho_Chl = param_Theta_max_N * P_phot_c / (3600.0 * param_Alpha_Chl * Theta_c * E_0)
  else:
    Rho_Chl = 0.0

  # Respiration
  R_C_growth = (vars['AmmoniumUptake'] + vars['NitrateUptake']) * param_Zeta / (stepInHours * vars['Carbon'])
  R_C = param_R_maintenance + R_C_growth

  # Cell division
  if vars['Carbon'] + (vars['Carbon'] * (P_phot_c - R_C) * stepInHours) >= param_C_rep and (vars['Silicate'] + vars['SilicateUptake']) >= param_S_rep:
    C_d = 2.0
  else:
    C_d = 1.0

  if C_d == 2.0:
    vars['Size'] = vars['Size'] * 2

  # Nutrients Uptake
  Q_nitrate = (vars['Nitrate'] + vars['NitrateUptake']) / vars['Carbon']
  Q_ammonium = (vars['Ammonium'] + vars['AmmoniumUptake']) / vars['Carbon']
  omega = (param_k_AR / (param_k_AR + env['DissolvedAmmonium'])) * ((param_k_AR + env['DissolvedNitrate']) / (param_k_AR + env['DissolvedAmmonium'] + env['DissolvedNitrate']))

  if vars['Ammonium'] + vars['Nitrate'] < 1000.0:
    if Q_ammonium + Q_nitrate < param_Q_Nmin:
      V_max_C = param_V_ref_c * T_function
    else:
      if Q_ammonium + Q_nitrate > param_Q_Nmax:
        V_max_C = 0.0
      else:
        V_max_C = param_V_ref_c * math.pow( (param_Q_Nmax - (Q_ammonium + Q_nitrate)) / (param_Q_Nmax - param_Q_Nmin), 0.05) * T_function
  else:
    V_max_C = 0.0

  V_C_ammonium = V_max_C * env['DissolvedAmmonium'] / (param_k_AR + env['DissolvedAmmonium'])
  V_C_nitrate = V_max_C * env['DissolvedNitrate'] / (param_k_AR + env['DissolvedNitrate']) * omega 

  if vars['Carbon'] >= param_C_minS:
    if Q_s <= param_Q_S_min:
      V_S_max = param_V_S_ref * T_function
    else:
      if Q_s >= param_Q_S_max:
        V_S_max = 0.0
      else:
        V_S_max = param_V_S_ref * math.pow( (param_Q_S_max - Q_s) / (param_Q_S_max - param_Q_S_min), 0.05) * T_function
  else:
    V_S_max = 0.0

  V_S_S = V_S_max * env['DissolvedSilicate'] / (env['DissolvedSilicate'] + param_k_S)

  ammonium_uptake_rate = c_old * vars['Carbon'] * V_C_ammonium * stepInHours
  nitrate_uptake_rate = c_old * vars['Carbon'] * V_C_nitrate * stepInHours
  silicate_uptake_rate = c_old * vars['Silicate'] * V_S_S * stepInHours

  # Update Pools
  ammonium_pool_new = (vars['Ammonium'] + vars['AmmoniumUptake'] + vars['NitrateUptake'] - (vars['Ammonium'] * param_R_N * stepInHours * T_function)) / C_d
  nitrate_pool_new = 0.0
  silicate_pool_new = (vars['Silicate'] + vars['SilicateUptake']) / C_d

  if vars['Carbon'] <= param_C_starve:
    chlorophyll_pool_new = 0.0
  else:
    if Theta_N <= param_Theta_max_N:
      chlorophyll_pool_new = vars['Chlorophyll'] + Rho_Chl * (vars['AmmoniumUptake'] + vars['NitrateUptake'])
    else:
      chlorophyll_pool_new = param_Theta_max_N * (vars['Ammonium'] + vars['Nitrate'])
    chlorophyll_pool_new = chlorophyll_pool_new - (vars['Chlorophyll'] * param_R_Chl * stepInHours * T_function)
    chlorophyll_pool_new = max(chlorophyll_pool_new / C_d, 0.0)

  if vars['Carbon'] <= param_C_starve:
    carbon_pool_new = 0.0
  else:
    carbon_pool_new = ((vars['Carbon'] * (P_phot_c - (R_C * T_function)) * stepInHours) + vars['Carbon']) / C_d

  # Mortality
  if vars['Carbon'] <= param_C_starve:
    vars['Stage'] = 1.0 # is Dead

  # Remineralisation Nitrogen
  ammonium_release_rate = c_old * (vars['Ammonium'] + vars['Nitrate']) * param_R_N * stepInHours * T_function

  # Housekeeping
  vars['Carbon'] = carbon_pool_new
  vars['Chlorophyll'] = chlorophyll_pool_new
  vars['Ammonium'] = ammonium_pool_new
  vars['Nitrate'] = nitrate_pool_new
  vars['Silicate'] = silicate_pool_new

  vars['AmmoniumUptake'] = ammonium_uptake_rate
  vars['NitrateUptake'] = nitrate_uptake_rate
  vars['SilicateUptake'] = silicate_uptake_rate

  vars['AmmoniumRelease'] = ammonium_release_rate
  vars['SilicateRelease'] = 0.0

def update_dead_diatom(vars, env, dt):
  """ Update kernel for dead Diatom agents from LERM-ES, 
      created by M. Sinerchia nd W. Hinsley

     References:
     "Testing theories of fisheries recruitment.", Sinerchia, M., 2007
     PhD Thesis, Department of Earth Science and Engineering, Imperial College, London

     "Using an individual-based model with four trophic levels to model the effect of predation and competition on squid recruitment", 
     Sinerchia, M., Field, A. J., Woods, J. D., Vallerga, S., and Hinsley, W. R.
     ICES Journal of Marine Science, doi:10.1093/icesjms/fsr190.
  """

  # Housekeeping
  c_old = vars['Size']
  vars['AmmoniumUptake'] = vars['AmmoniumUptake'] / c_old
  vars['NitrateUptake'] = vars['NitrateUptake'] / c_old
  vars['SilicateUptake'] = vars['SilicateUptake'] / c_old

  stepInHours = dt/3600.

  Si_reminT = param_S_dis * math.pow(param_Q_remS, (env['Temperature'] + 273.0 - param_T_refS) / 10.0)
  N_reminT = param_Ndis * math.pow(param_Q_remN, (env['Temperature'] + 273.0 - param_T_refN) / 10.0)
  silicate_release_rate = c_old * vars['Silicate'] * Si_reminT * stepInHours
  silicate_pool_new = max(vars['Silicate'] - (vars['Silicate'] * Si_reminT * stepInHours), 0.0)
  ammonium_release_rate = c_old * (vars['Ammonium'] + vars['Nitrate']) * N_reminT * stepInHours
  ammonium_pool_new = max(vars['Ammonium'] - (vars['Ammonium'] * N_reminT * stepInHours), 0.0)
  nitrate_pool_new = max(vars['Nitrate'] - (vars['Nitrate'] * N_reminT * stepInHours), 0.0)

  # Housekeeping
  vars['Ammonium'] = ammonium_pool_new
  vars['Nitrate'] = nitrate_pool_new
  vars['Silicate'] = silicate_pool_new

  vars['AmmoniumUptake'] = 0.0
  vars['NitrateUptake'] = 0.0
  vars['SilicateUptake'] = 0.0

  vars['AmmoniumRelease'] = ammonium_release_rate
  vars['SilicateRelease'] = silicate_release_rate



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
