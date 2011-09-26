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
  
  ### Phase 1 ###
  
  # Housekeeping
  vars['AmmoniumIng'] = vars['AmmoniumIng'] / vars['Biomass']
  vars['NitrateIng'] = vars['NitrateIng'] / vars['Biomass']
  vars['SilicateIng'] = vars['SilicateIng'] / vars['Biomass']
  stepInHours = dt/3600.
  biomass_old = vars['Biomass']
  
  # Nutrient uptake 
  vars['AmmoniumIng'] = vars['AmmoniumIng'] * env['AmmoniumDepletion']
  vars['NitrateIng'] = vars['NitrateIng'] * env['NitrateDepletion']
  vars['SilicateIng'] = vars['SilicateIng'] * env['SilicateDepletion']
  
  ### Phase 2 ###

  T_K = env['Temperature'] + 273.0
  T_function = math.exp(34.12969283 - 10000/T_K)

  # Photosynthesis
  Q_N = ((((vars['AmmoniumPool'])+(vars['AmmoniumIng'])+(vars['NitratePool'])+(vars['NitrateIng']))) / (vars['CarbonPool']))
  Q_s = ((((vars['SilicatePool'])+(vars['SilicateIng']))) / (vars['CarbonPool']))

  if ((Q_N) > (param_Q_Nmax)):
    P_max_c = (((param_P_ref_c)*(T_function)))
  else:
    if ((Q_N) < (param_Q_Nmin)):
      P_max_c = 0.0
    else:
      P_max_c = (((param_P_ref_c)*(T_function)*(((Q_N) - (param_Q_Nmin)) / ((param_Q_Nmax) - (param_Q_Nmin)))))

  Theta_c = ((vars['ChlorophyllPool']) / (vars['CarbonPool']))
  E_0 = (((4.6)*(env['Irradiance'])))

  if ((P_max_c == 0.0)or(Q_s <= param_Q_S_min)):
    P_phot_c = 0.0
  else:
    P_phot_c = P_max_c*(1.0 - math.exp(-0.284400e-2*Theta_c*E_0 / P_max_c))

  # Reproduction
  Theta_N = ((vars['ChlorophyllPool']) / (((vars['AmmoniumPool'])+(vars['AmmoniumIng'])+(vars['NitratePool'])+(vars['NitrateIng']))))

  if ((((E_0) > (0.0))and((Theta_c) > (0.0)))):
    Rho_Chl = (((param_Theta_max_N)*((P_phot_c) / (((3600.0)*(param_Alpha_Chl)*(Theta_c)*(E_0))))))
  else:
    Rho_Chl = 0.0

  R_C_growth = ((((((vars['AmmoniumIng'])+(vars['NitrateIng'])))*(param_Zeta))) / (((stepInHours)*(vars['CarbonPool']))))
  R_C = (((param_R_maintenance)+(R_C_growth)));

  if ((((((vars['CarbonPool'])+(((vars['CarbonPool'])*((P_phot_c) - (R_C))*(stepInHours))))) >= (param_C_rep))and((((vars['SilicatePool'])+(vars['SilicateIng']))) >= (param_S_rep)))):
    C_d = 2.0
  else:
    C_d = 1.0

  # Update chemical pools 
  Q_nitrate = ((((vars['NitratePool'])+(vars['NitrateIng']))) / (vars['CarbonPool']))
  Q_ammonium = ((((vars['AmmoniumPool'])+(vars['AmmoniumIng']))) / (vars['CarbonPool']))
  omega = ((((param_k_AR) / (((param_k_AR)+(env['Ammonium']))))*((((param_k_AR)+(env['Nitrate']))) / (((param_k_AR)+(env['Ammonium'])+(env['Nitrate']))))))

  if ((((vars['AmmoniumPool'])+(vars['NitratePool']))) < (1000.0)):
    if ((((Q_ammonium)+(Q_nitrate))) < (param_Q_Nmin)):
      V_max_C = (((param_V_ref_c)*(T_function)))
    else:
      if ((((Q_ammonium)+(Q_nitrate))) > (param_Q_Nmax)):
        V_max_C = 0.0
      else:
        V_max_C = (((param_V_ref_c)*(math.pow(((param_Q_Nmax) - (((Q_ammonium)+(Q_nitrate)))) / ((param_Q_Nmax) - (param_Q_Nmin)), 0.05))*(T_function)))
  else:
    V_max_C = 0.0

  V_C_ammonium = (((V_max_C)*((env['Ammonium']) / (((param_k_AR)+(env['Ammonium']))))))
  V_C_nitrate = (((V_max_C)*((env['Nitrate']) / (((param_k_AR)+(env['Nitrate']))))*(omega)))

  if ((vars['CarbonPool']) >= (param_C_minS)):
    if ((Q_s) <= (param_Q_S_min)):
      V_S_max = (((param_V_S_ref)*(T_function)))
    else:
      if ((Q_s) >= (param_Q_S_max)):
        V_S_max = 0.0
      else:
         V_S_max = (((param_V_S_ref)*(math.pow(((param_Q_S_max) - (Q_s)) / ((param_Q_S_max) - (param_Q_S_min)), 0.05))*(T_function)))
  else:
    V_S_max = 0.0

  V_S_S = (((V_S_max)*((env['Silicate']) / (((env['Silicate'])+(param_k_S))))))

  # Ingestion rates
  IngAmmonium = biomass_old*((vars['CarbonPool'])*(V_C_ammonium)*(stepInHours))
  IngNitrate = biomass_old*((vars['CarbonPool'])*(V_C_nitrate)*(stepInHours))
  IngSilicate = biomass_old*((vars['SilicatePool'])*(V_S_S)*(stepInHours))

  # Ammonium excretion 
  RelAmmonium = biomass_old*((((vars['AmmoniumPool'])+(vars['NitratePool'])))*(param_R_N)*(stepInHours)*(T_function))

  AmmoniumPoolNew = (((((vars['AmmoniumPool'])+(vars['AmmoniumIng'])+(vars['NitrateIng']))) - (((vars['AmmoniumPool'])*(param_R_N)*(stepInHours)*(T_function)))) / (C_d))
  NitratePoolNew = (0.0)
  SilicatePoolNew = ((((vars['SilicatePool'])+(vars['SilicateIng']))) / (C_d))

  # Mortality, determine death...
  if ((vars['CarbonPool']) <= (param_C_starve)):
    death_flag = True
  else:
    death_flag =  False

  # ...before updating carbon pool and chlorophyll
  if (death_flag):
    ###  state change! ###
    vars['Stage'] = 1.0

    CarbonPoolNew = 0.0
    ChlorophyllPoolNew = 0.0;
  else:
    CarbonPoolNew = (((((((vars['CarbonPool'])*((P_phot_c)-(((R_C)*(T_function))))*(stepInHours)))+(vars['CarbonPool']))) - (0.0)) / (C_d))

    if (Theta_N <= param_Theta_max_N):
      ChlorophyllPoolNew = vars['ChlorophyllPool'] + Rho_Chl * (vars['AmmoniumIng'] + vars['NitrateIng'])
    else:
      ChlorophyllPoolNew = param_Theta_max_N * (vars['AmmoniumPool'] + vars['NitratePool'])
    ChlorophyllPoolNew = ChlorophyllPoolNew - (vars['ChlorophyllPool'] * param_R_Chl * stepInHours * T_function)
    if (ChlorophyllPoolNew > 0.0):
      ChlorophyllPoolNew = ChlorophyllPoolNew / C_d
    else:
      ChlorophyllPoolNew = 0.0

  ### Phase 3 (output and housekeeping) ###

  # Ensemble update (reproduction) */
  if ((C_d) == (2.0)):
    vars['Biomass'] = biomass_old*2

  # External
  if IngAmmonium < 0.0:
    IngAmmonium = 0.0
  if IngNitrate < 0.0:
    IngNitrate = 0.0
  if IngSilicate < 0.0:
    IngSilicate = 0.0
  if RelAmmonium < 0.0:
    RelAmmonium = 0.0
  vars['AmmoniumIng'] = IngAmmonium
  vars['NitrateIng'] = IngNitrate
  vars['SilicateIng'] = IngSilicate
  vars['AmmoniumRel'] = RelAmmonium
  vars['SilicateRel'] = 0.0

  # Internal
  vars['AmmoniumPool'] = AmmoniumPoolNew
  vars['NitratePool'] = NitratePoolNew
  vars['SilicatePool'] = SilicatePoolNew
  vars['CarbonPool'] = CarbonPoolNew
  vars['ChlorophyllPool'] = ChlorophyllPoolNew



def update_dead_diatom(vars, env, dt):

  ### Phase 1 ###
  
  # Housekeeping
  vars['AmmoniumIng'] = vars['AmmoniumIng'] / vars['Biomass']
  vars['NitrateIng'] = vars['NitrateIng'] / vars['Biomass']
  vars['SilicateIng'] = vars['SilicateIng'] / vars['Biomass']
  stepInHours = dt/3600.
  biomass_old = vars['Biomass']
  
  # Nutrient uptake 
  vars['AmmoniumIng'] = vars['AmmoniumIng'] * env['AmmoniumDepletion']
  vars['NitrateIng'] = vars['NitrateIng'] * env['NitrateDepletion']
  vars['SilicateIng'] = vars['SilicateIng'] * env['SilicateDepletion']

  ### Phase 2 ###

  Si_reminT = (((param_S_dis)*(math.pow(param_Q_remS, ((((env['Temperature'])+(273.0))) - (param_T_refS)) / (10.0)))))
  N_reminT = (((param_Ndis)*(math.pow(param_Q_remN, ((((env['Temperature'])+(273.0))) - (param_T_refN)) / (10.0)))))
  relAmountSi = biomass_old*((vars['SilicatePool'])*(Si_reminT)*(stepInHours));

  SilicatePoolNew = ((((vars['SilicatePool'])+(vars['SilicateIng']))) - (((vars['SilicatePool'])*(Si_reminT)*(stepInHours))))
  if (SilicatePoolNew < 0.0):
    SilicatePoolNew = 0.0
  relAmountAmm = biomass_old*((((vars['AmmoniumPool'])+(vars['NitratePool'])))*(N_reminT)*(stepInHours))

  AmmoniumPoolNew = ((((vars['AmmoniumPool'])+(vars['AmmoniumIng']))) - (((vars['AmmoniumPool'])*(N_reminT)*(stepInHours))))
  if (AmmoniumPoolNew < 0.0):
    AmmoniumPoolNew = 0.0
  NitratePoolNew = ((((vars['NitratePool'])+(vars['NitrateIng']))) - (((vars['NitratePool'])*(N_reminT)*(stepInHours))))
  if (NitratePoolNew < 0.0):
    NitratePoolNew = 0.0

  ### Phase 3 (output and housekeeping) ###

  # External
  if relAmountAmm < 0.0:
    relAmountAmm = 0.0
  if relAmountSi < 0.0:
    relAmountSi = 0.0
  vars['AmmoniumIng'] = 0.0
  vars['NitrateIng'] = 0.0
  vars['SilicateIng'] = 0.0
  vars['AmmoniumRel'] = relAmountAmm
  vars['SilicateRel'] = relAmountSi

  # Internal
  vars['AmmoniumPool'] = AmmoniumPoolNew
  vars['NitratePool'] = NitratePoolNew
  vars['SilicatePool'] = SilicatePoolNew
  vars['CarbonPool'] = 0.0
  vars['ChlorophyllPool'] = 0.0

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
