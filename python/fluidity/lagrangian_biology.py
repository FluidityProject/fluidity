import os
import sys
import math

import math
import numpy
from lebiology import stage_id, add_agent

# Parameters for FGroup Diatom
# Species: Default_Diatom_Variety
species_Default_Diatom_Variety = {
    'A_E' : -10000.0,
    'Alpha_Chl' : 7.9e-07,
    'C_minS' : 1.584e-08,
    'C_rep' : 1.76e-08,
    'C_starve' : 8.5e-09,
    'C_struct' : 8.5e-09,
    'Ndis' : 0.0042,
    'P_ref_c' : 0.14,
    'Q_Nmax' : 0.17,
    'Q_Nmin' : 0.034,
    'Q_S_max' : 0.15,
    'Q_S_min' : 0.04,
    'Q_remN' : 2.95,
    'Q_remS' : 2.27,
    'R_Chl' : 0.002,
    'R_N' : 0.002,
    'R_maintenance' : 0.002,
    'S_dis' : 0.00083,
    'S_rep' : 2.1e-09,
    'T_ref' : 293.0,
    'T_refN' : 283.0,
    'T_refS' : 278.0,
    'Theta_max_N' : 4.2,
    'V_S_ref' : 0.03,
    'V_ref_c' : 0.01,
    'Zeta' : 2.3,
    'k_AR' : 1.0,
    'k_S' : 1.0,
    'z_sink' : 0.04,
}

def update_Living_Diatom_no_sil(param, vars, env, dt):
  """ FGroup:  Diatom
      Stage:   Living
  """
  dt_in_hours = dt / 3600.0

  ### Effect of temperature ###
  T_K = (env['Temperature'] + 273.0)
  T_function = math.exp((param['A_E'] * ((1.0 / T_K) - (1.0 / param['T_ref']))))

  ### Photosynthesis ###
  Q_N = ((vars['Ammonium'] + vars['AmmoniumIngested'] + vars['Nitrate'] + vars['NitrateIngested']) / vars['Carbon'])
  #Q_s = ((vars['Silicate'] + vars['SilicateIngested']) / vars['Carbon'])
  P_max_c = (((param['P_ref_c'] * T_function)) if (Q_N > param['Q_Nmax']) else (((0.0) if (Q_N < param['Q_Nmin']) else ((param['P_ref_c'] * T_function * ((Q_N - param['Q_Nmin']) / (param['Q_Nmax'] - param['Q_Nmin'])))))))
  Theta_c = (vars['Chlorophyll'] / vars['Carbon'])
  E_0 = (4.6 * env['Irradiance'])
  #P_phot_c = ((0.0) if ((P_max_c == 0.0) or (Q_s <= param['Q_S_min'])) else ((P_max_c * (1.0 - math.exp(((-3600.0 * param['Alpha_Chl'] * Theta_c * E_0) / P_max_c))))))
  P_phot_c = ((0.0) if ((P_max_c == 0.0)) else ((P_max_c * (1.0 - math.exp(((-3600.0 * param['Alpha_Chl'] * Theta_c * E_0) / P_max_c))))))

  ### Chlorophyll Synthesis ###
  Theta_N = (vars['Chlorophyll'] / (vars['Ammonium'] + vars['AmmoniumIngested'] + vars['Nitrate'] + vars['NitrateIngested']))
  Rho_Chl = (((param['Theta_max_N'] * (P_phot_c / (3600.0 * param['Alpha_Chl'] * Theta_c * E_0)))) if ((E_0 > 0.0) and (Theta_c > 0.0)) else (0.0))

  ### Respiration ###
  R_C_growth = (((vars['AmmoniumIngested'] + vars['NitrateIngested']) * param['Zeta']) / (dt_in_hours * vars['Carbon']))
  R_C = (param['R_maintenance'] + R_C_growth)

  ### Cell Division ###
  #C_d = ((2.0) if (((vars['Carbon'] + (vars['Carbon'] * (P_phot_c - R_C) * dt_in_hours)) >= param['C_rep']) and ((vars['Silicate'] + vars['SilicateIngested']) >= param['S_rep'])) else (1.0))
  C_d = ((2.0) if (((vars['Carbon'] + (vars['Carbon'] * (P_phot_c - R_C) * dt_in_hours)) >= param['C_rep'])) else (1.0))
  if (C_d == 2.0):
    vars['Size'] = vars['Size'] * 2.0

  ### Nutrients uptake ###
  Q_nitrate = ((vars['Nitrate'] + vars['NitrateIngested']) / vars['Carbon'])
  Q_ammonium = ((vars['Ammonium'] + vars['AmmoniumIngested']) / vars['Carbon'])
  omega = ((param['k_AR'] / (param['k_AR'] + env['DissolvedAmmonium'])) * ((param['k_AR'] + env['DissolvedNitrate']) / (param['k_AR'] + env['DissolvedAmmonium'] + env['DissolvedNitrate'])))
  V_max_C = (((((param['V_ref_c'] * T_function)) if ((Q_ammonium + Q_nitrate) < param['Q_Nmin']) else (((0.0) if ((Q_ammonium + Q_nitrate) > param['Q_Nmax']) else ((param['V_ref_c'] * math.pow(((param['Q_Nmax'] - (Q_ammonium + Q_nitrate)) / (param['Q_Nmax'] - param['Q_Nmin'])), 0.05) * T_function)))))) if ((vars['Ammonium'] + vars['Nitrate']) < 1000.0) else (0.0))
  V_C_ammonium = (V_max_C * (env['DissolvedAmmonium'] / (param['k_AR'] + env['DissolvedAmmonium'])))
  V_C_nitrate = (V_max_C * (env['DissolvedNitrate'] / (param['k_AR'] + env['DissolvedNitrate'])) * omega)
  #V_S_max = (((((param['V_S_ref'] * T_function)) if (Q_s <= param['Q_S_min']) else (((0.0) if (Q_s >= param['Q_S_max']) else ((param['V_S_ref'] * math.pow(((param['Q_S_max'] - Q_s) / (param['Q_S_max'] - param['Q_S_min'])), 0.05) * T_function)))))) if (vars['Carbon'] >= param['C_minS']) else (0.0))

  # ml805 Disabling Silicate limitation
  #V_S_S = (V_S_max * (env['DissolvedSilicate'] / (env['DissolvedSilicate'] + param['k_S'])))
  vars['AmmoniumUptake'] = (vars['Carbon'] * V_C_ammonium * dt_in_hours)
  vars['NitrateUptake'] = (vars['Carbon'] * V_C_nitrate * dt_in_hours)
  #vars['SilicateUptake'] = (vars['Silicate'] * V_S_S * dt_in_hours)

  ### Update Pools ###
  Ammonium_new = ((((vars['Ammonium'] + vars['AmmoniumIngested'] + vars['NitrateIngested']) - (vars['Ammonium'] * param['R_N'] * dt_in_hours * T_function)) - (((0.0 * (Q_N - param['Q_Nmax']))) if (Q_N > param['Q_Nmax']) else (0.0))) / C_d)
  Nitrate_new = 0.0
  #Silicate_new = (((vars['Silicate'] + vars['SilicateIngested']) - 0.0) / C_d)
  C_new = max(0.0, (((vars['Carbon'] * (P_phot_c - (R_C * T_function)) * dt_in_hours) + vars['Carbon']) / C_d))
  death_flag = ((1.0) if (C_new <= param['C_starve']) else (0.0))
  if (death_flag == 1.0):
    vars['Stage'] = stage_id('Diatom', 'Dead')
  Carbon_new = C_new
  Chlorophyll_new = ((max((((((vars['Chlorophyll'] + (Rho_Chl * (vars['AmmoniumIngested'] + vars['NitrateIngested'])))) if (Theta_N <= param['Theta_max_N']) else ((vars['Chlorophyll'] - (vars['Chlorophyll'] - ((vars['Ammonium'] + vars['Nitrate']) * param['Theta_max_N']))))) - ((vars['Chlorophyll'] * param['R_Chl'] * dt_in_hours * T_function) + 0.0)) / C_d), 0.0)) if (death_flag == 0.0) else (0.0))
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])
  C_fuel_new = (vars['Carbon'] - param['C_struct'])

  ### Remineralisation Nitrogen ###
  vars['AmmoniumRelease'] = ((vars['Ammonium'] + vars['Nitrate']) * param['R_N'] * dt_in_hours * T_function)

  ### Photoadaptation ###
  ChltoC_new = (((vars['Chlorophyll'] / vars['Carbon'])) if (vars['Carbon'] > 0.0) else (0.0))
  NtoC_new = ((((vars['Ammonium'] + vars['Nitrate']) / vars['Carbon'])) if (vars['Carbon'] > 0.0) else (0.0))

  ### Setting pool variables
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  #vars['Silicate'] = Silicate_new
  vars['Carbon'] = Carbon_new
  vars['Chlorophyll'] = Chlorophyll_new
  vars['Nitrogen'] = Nitrogen_new
  vars['C_fuel'] = C_fuel_new
  vars['ChltoC'] = ChltoC_new
  vars['NtoC'] = NtoC_new

def update_Dead_Diatom_no_sil(param, vars, env, dt):
  """ FGroup:  Diatom
      Stage:   Dead
  """
  dt_in_hours = dt / 3600.0

  ### Remineralisation Dead T ###
  #Si_reminT = (param['S_dis'] * math.pow(param['Q_remS'], (((env['Temperature'] + 273.0) - param['T_refS']) / 10.0)))
  N_reminT = (param['Ndis'] * math.pow(param['Q_remN'], (((env['Temperature'] + 273.0) - param['T_refN']) / 10.0)))
  #vars['SilicateRelease'] = (vars['Silicate'] * Si_reminT * dt_in_hours)
  #Silicate_new = max((vars['Silicate'] - (vars['Silicate'] * Si_reminT * dt_in_hours)), 0.0)
  vars['AmmoniumRelease'] = ((vars['Ammonium'] + vars['Nitrate']) * N_reminT * dt_in_hours)
  Ammonium_new = max((vars['Ammonium'] - (vars['Ammonium'] * N_reminT * dt_in_hours)), 0.0)
  Nitrate_new = max((vars['Nitrate'] - (vars['Nitrate'] * N_reminT * dt_in_hours)), 0.0)

  ### Setting pool variables
  #vars['Silicate'] = Silicate_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new


###  PCZDNA  ###

def six_component_zdna(state, parameters):
    '''Calculate sources and sinks for pczdna biology model'''    

    # Based on the equations in
    # Popova, E. E.; Coward, A. C.; Nurser, G. A.; de Cuevas, B.; Fasham, M. J. R. & Anderson, T. R. 
    # Mechanisms controlling primary and new production in a global ecosystem model - Part I: 
    # Validation of the biological simulation Ocean Science, 2006, 2, 249-266. 
    # DOI: 10.5194/os-2-249-2006

    ### ml805 note:
    # Removed Chlorophyll and using P only as input for Z, 
    # so that P and C may be replaced by an agent set
    import math

    #if not check_six_component_parameters(parameters):
    #    raise TypeError("Missing Parameter")
    
    P=state.scalar_fields["Phytoplankton"]
    Z=state.scalar_fields["Zooplankton"]
    N=state.scalar_fields["Nutrient"]
    A=state.scalar_fields["Ammonium"]
    D=state.scalar_fields["Detritus"]
    I=state.scalar_fields["_PAR"] # Note: *NOT* PhotosyntheticRadation field, but the internal _PAR field, which 
                                  # has been projected onto the same mesh as phytoplankton and has been converted from
                                  # incident radiation to just the active part.
    Znew=state.scalar_fields["IteratedZooplankton"]
    Nnew=state.scalar_fields["IteratedNutrient"]
    Anew=state.scalar_fields["IteratedAmmonium"]
    Dnew=state.scalar_fields["IteratedDetritus"]
    usingDistToTop = False
    try:
        distanceToTop=state.scalar_fields["DTT"]
        usingDistToTop = True
    except KeyError:
        coords=state.vector_fields["Coordinate"]

    Z_source=state.scalar_fields["ZooplanktonSource"]
    N_source=state.scalar_fields["NutrientSource"]
    N_abs=state.scalar_fields["NutrientAbsorption"]
    A_source=state.scalar_fields["AmmoniumSource"]
    D_source=state.scalar_fields["DetritusSource"]
    try:
        PP=state.scalar_fields["PrimaryProduction"]
    except KeyError:
        PP=None
    try:
        PG=state.scalar_fields["PhytoplanktonGrazing"]
    except KeyError:
        PG=None

    alpha_c=parameters["alpha_c"]
    beta_P=parameters["beta_p"]
    beta_D=parameters["beta_d"]
    delta=parameters["delta"]
    gamma=parameters["gamma"]
    zeta=parameters["zeta"]
    epsilon=parameters["epsilon"]
    psi=parameters["psi"]
    g=parameters["g"]
    k_N=parameters["k_N"]
    k_A=parameters["k_A"]
    k_p=parameters["k_p"]
    k_z=parameters["k_z"]
    v=parameters["v"]
    mu_P=parameters["mu_P"]
    mu_Z=parameters["mu_Z"]
    mu_D=parameters["mu_D"]
    p_P=parameters["p_P"]
    theta_m=parameters["theta_m"]
    lambda_bio=parameters["lambda_bio"]
    lambda_A=parameters["lambda_A"]
    photicZoneLimit=parameters["photic_zone_limit"]
    p_D=1-p_P

    for n in range(P.node_count):
        # Values of fields on this node.
        P_n=max((P.node_val(n)), 0.0)
        Z_n=max(.5*(Z.node_val(n)+Znew.node_val(n)), 0.0)
        N_n=max(.5*(N.node_val(n)+Nnew.node_val(n)), 0.0)
        A_n=max(.5*(A.node_val(n)+Anew.node_val(n)), 0.0)
        D_n=max(.5*(D.node_val(n)+Dnew.node_val(n)), 0.0)
        I_n=max(I.node_val(n), 0.0)
        if (usingDistToTop):
            depth = distanceToTop.node_val(n)
        else:
            depth=abs(coords.node_val(n)[-1])

        if (I_n < 0.0001):
            I_n =0

        # Zooplankton grazing of phytoplankton.
	    # It looks a bit different from the original version, however
	    # it is the same function with differently normalised parameters to 
	    # simplify tuning 
        # G_P=(g * epsilon * p_P * P_n**2 * Z_n)/(g+epsilon*(p_P*P_n**2 + p_D*D_n**2))
        G_P=(g * p_P * P_n**2 * Z_n)/(epsilon + (p_P*P_n**2 + p_D*D_n**2))

        # Zooplankton grazing of detritus. (p_D - 1-p_P)
        # G_D=(g * epsilon * (1-p_P) * D_n**2 * Z_n)/(g+epsilon*(p_P*P_n**2 + p_D*D_n**2))
        G_D=(g  * (1-p_P) * D_n**2 * Z_n)/(epsilon + (p_P*P_n**2 + p_D*D_n**2))

        # Death rate of zooplankton.
	    # There is an additional linear term because we have a unified model
	    # (no below/above photoc zone distinction)
        De_Z=mu_Z*Z_n**3/(Z_n+k_z)+lambda_bio*Z_n

        # Detritus remineralisation.
        #De_D=mu_D*D_n+lambda_bio*P_n+lambda_bio*Z_n
        De_D=mu_D*D_n+lambda_bio*Z_n

        # Ammonium nitrification (only below the photic zone)
	    # This is the only above/below term
        De_A=lambda_A*A_n*(1-photic_zone(depth,100,20))

        # ml805: Disabled phytoplankton impact...
        Z_source.set(n, delta*(beta_P*G_P+beta_D*G_D) - De_Z)
        D_source.set(n, -De_D + gamma*De_Z +(1-beta_P)*G_P - beta_D*G_D)
        #D_source.set(n, -De_D + De_P + gamma*De_Z +(1-beta_P)*G_P - beta_D*G_D)
        #N_source.set(n, -J*P_n*Q_N+De_A)
        N_source.set(n, De_A)
        #A_source.set(n, -J*P_n*Q_A + De_D + (1 - delta)*(beta_P*G_P + beta_D*G_D) + (1-gamma)*De_Z-De_A)
        A_source.set(n, + De_D + (1 - delta)*(beta_P*G_P + beta_D*G_D) + (1-gamma)*De_Z-De_A)

        if PG:
            PG.set(n, G_P)



def photic_zone(z,limit,transition_length):

    depth = abs(z)
    if (depth < limit):
        return 1.
    elif (depth < limit+transition_length):
        return 1.-(depth-limit)/float(transition_length)
    else:
        return 0.0

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
