import math

# Parameters for FGroup Diatom
# Species: Default_Diatom_Variety
params_Default_Diatom_Variety = {
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

def update_Living_Diatom(param, vars, env, dt):
  """ FGroup:  Diatom
      Stage:   Living   ID: 0
  """
  dt_in_hours = dt / 3600.0

  ### Effect of temperature ###
  T_K = (env['Temperature'] + 273.0)
  T_function = math.exp((param['A_E'] * ((1.0 / T_K) - (1.0 / param['T_ref']))))

  ### Photosynthesis ###
  Q_N = ((vars['Ammonium'] + vars['AmmoniumIngested'] + vars['Nitrate'] + vars['NitrateIngested']) / vars['Carbon'])
  Q_s = ((vars['Silicate'] + vars['SilicateIngested']) / vars['Carbon'])
  P_max_c = (((param['P_ref_c'] * T_function)) if ((Q_N > param['Q_Nmax'])) else (((0.0) if ((Q_N < param['Q_Nmin'])) else ((param['P_ref_c'] * T_function * ((Q_N - param['Q_Nmin']) / (param['Q_Nmax'] - param['Q_Nmin'])))))))
  Theta_c = (vars['Chlorophyll'] / vars['Carbon'])
  E_0 = (4.6 * env['Irradiance'])
  P_phot_c = ((0.0) if (((P_max_c == 0.0)) or ((Q_s <= param['Q_S_min']))) else ((P_max_c * (1.0 - math.exp(((-3600.0 * param['Alpha_Chl'] * Theta_c * E_0) / P_max_c))))))

  ### Chlorophyll Synthesis ###
  Theta_N = (vars['Chlorophyll'] / (vars['Ammonium'] + vars['AmmoniumIngested'] + vars['Nitrate'] + vars['NitrateIngested']))
  Rho_Chl = (((param['Theta_max_N'] * (P_phot_c / (3600.0 * param['Alpha_Chl'] * Theta_c * E_0)))) if (((E_0 > 0.0)) and ((Theta_c > 0.0))) else (0.0))

  ### Respiration ###
  R_C_growth = (((vars['AmmoniumIngested'] + vars['NitrateIngested']) * param['Zeta']) / (dt_in_hours * vars['Carbon']))
  R_C = (param['R_maintenance'] + R_C_growth)

  ### Cell Division ###
  C_d = ((2.0) if ((((vars['Carbon'] + (vars['Carbon'] * (P_phot_c - R_C) * dt_in_hours)) >= param['C_rep'])) and (((vars['Silicate'] + vars['SilicateIngested']) >= param['S_rep']))) else (1.0))
  if (C_d == 2.0):
    vars['Size'] = vars['Size'] * 2.0

  ### Nutrients uptake ###
  Q_nitrate = ((vars['Nitrate'] + vars['NitrateIngested']) / vars['Carbon'])
  Q_ammonium = ((vars['Ammonium'] + vars['AmmoniumIngested']) / vars['Carbon'])
  omega = ((param['k_AR'] / (param['k_AR'] + env['DissolvedAmmonium'])) * ((param['k_AR'] + env['DissolvedNitrate']) / (param['k_AR'] + env['DissolvedAmmonium'] + env['DissolvedNitrate'])))
  V_max_C = (((((param['V_ref_c'] * T_function)) if (((Q_ammonium + Q_nitrate) < param['Q_Nmin'])) else (((0.0) if (((Q_ammonium + Q_nitrate) > param['Q_Nmax'])) else ((param['V_ref_c'] * math.pow(((param['Q_Nmax'] - (Q_ammonium + Q_nitrate)) / (param['Q_Nmax'] - param['Q_Nmin'])), 0.05) * T_function)))))) if (((vars['Ammonium'] + vars['Nitrate']) < 1000.0)) else (0.0))
  V_C_ammonium = (V_max_C * (env['DissolvedAmmonium'] / (param['k_AR'] + env['DissolvedAmmonium'])))
  V_C_nitrate = (V_max_C * (env['DissolvedNitrate'] / (param['k_AR'] + env['DissolvedNitrate'])) * omega)
  V_S_max = (((((param['V_S_ref'] * T_function)) if ((Q_s <= param['Q_S_min'])) else (((0.0) if ((Q_s >= param['Q_S_max'])) else ((param['V_S_ref'] * math.pow(((param['Q_S_max'] - Q_s) / (param['Q_S_max'] - param['Q_S_min'])), 0.05) * T_function)))))) if ((vars['Carbon'] >= param['C_minS'])) else (0.0))
  V_S_S = (V_S_max * (env['DissolvedSilicate'] / (env['DissolvedSilicate'] + param['k_S'])))
  vars['AmmoniumUptake'] = (vars['Carbon'] * V_C_ammonium * dt_in_hours)
  vars['NitrateUptake'] = (vars['Carbon'] * V_C_nitrate * dt_in_hours)
  vars['SilicateUptake'] = (vars['Silicate'] * V_S_S * dt_in_hours)

  ### Update Pools ###
  Ammonium_new = ((((vars['Ammonium'] + vars['AmmoniumIngested'] + vars['NitrateIngested']) - (vars['Ammonium'] * param['R_N'] * dt_in_hours * T_function)) - (((0.0 * (Q_N - param['Q_Nmax']))) if ((Q_N > param['Q_Nmax'])) else (0.0))) / C_d)
  Nitrate_new = 0.0
  Silicate_new = (((vars['Silicate'] + vars['SilicateIngested']) - 0.0) / C_d)
  C_new = max(0.0, (((vars['Carbon'] * (P_phot_c - (R_C * T_function)) * dt_in_hours) + vars['Carbon']) / C_d))
  death_flag = ((1.0) if ((C_new <= param['C_starve'])) else (0.0))
  if (death_flag == 1.0):
    vars['Stage'] = 1.0  # Dead
  Carbon_new = C_new
  Chlorophyll_new = ((max((((((vars['Chlorophyll'] + (Rho_Chl * (vars['AmmoniumIngested'] + vars['NitrateIngested'])))) if ((Theta_N <= param['Theta_max_N'])) else ((vars['Chlorophyll'] - (vars['Chlorophyll'] - ((vars['Ammonium'] + vars['Nitrate']) * param['Theta_max_N']))))) - ((vars['Chlorophyll'] * param['R_Chl'] * dt_in_hours * T_function) + 0.0)) / C_d), 0.0)) if ((death_flag == 0.0)) else (0.0))
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])
  C_fuel_new = (vars['Carbon'] - param['C_struct'])

  ### Remineralisation Nitrogen ###
  vars['AmmoniumRelease'] = ((vars['Ammonium'] + vars['Nitrate']) * param['R_N'] * dt_in_hours * T_function)

  ### Photoadaptation ###
  ChltoC_new = (((vars['Chlorophyll'] / vars['Carbon'])) if ((vars['Carbon'] > 0.0)) else (0.0))
  NtoC_new = ((((vars['Ammonium'] + vars['Nitrate']) / vars['Carbon'])) if ((vars['Carbon'] > 0.0)) else (0.0))

  ### Setting pool variables
  vars['NtoC'] = NtoC_new
  vars['Chlorophyll'] = Chlorophyll_new
  vars['C_fuel'] = C_fuel_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Carbon'] = Carbon_new
  vars['Nitrate'] = Nitrate_new
  vars['Silicate'] = Silicate_new
  vars['Ammonium'] = Ammonium_new
  vars['ChltoC'] = ChltoC_new

def update_Dead_Diatom(param, vars, env, dt):
  """ FGroup:  Diatom
      Stage:   Dead   ID: 1
  """
  dt_in_hours = dt / 3600.0

  ### Remineralisation Dead T ###
  Si_reminT = (param['S_dis'] * math.pow(param['Q_remS'], (((env['Temperature'] + 273.0) - param['T_refS']) / 10.0)))
  N_reminT = (param['Ndis'] * math.pow(param['Q_remN'], (((env['Temperature'] + 273.0) - param['T_refN']) / 10.0)))
  vars['SilicateRelease'] = (vars['Silicate'] * Si_reminT * dt_in_hours)
  Silicate_new = max((vars['Silicate'] - (vars['Silicate'] * Si_reminT * dt_in_hours)), 0.0)
  vars['AmmoniumRelease'] = ((vars['Ammonium'] + vars['Nitrate']) * N_reminT * dt_in_hours)
  Ammonium_new = max((vars['Ammonium'] - (vars['Ammonium'] * N_reminT * dt_in_hours)), 0.0)
  Nitrate_new = max((vars['Nitrate'] - (vars['Nitrate'] * N_reminT * dt_in_hours)), 0.0)

  ### Setting pool variables
  vars['Nitrate'] = Nitrate_new
  vars['Silicate'] = Silicate_new
  vars['Ammonium'] = Ammonium_new
