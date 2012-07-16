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

def update_Living_Diatom(param, vars, env, dt):
  """ FGroup:  Diatom
      Stage:   Living
  """
  dt_in_hours = dt / 3600.0

  ### Effect of temperature ###
  T_K = (env['Temperature'] + 273.0)
  T_function = math.exp((param['A_E'] * ((1.0 / T_K) - (1.0 / param['T_ref']))))

  ### Photosynthesis ###
  Q_N = ((vars['Ammonium'] + vars['AmmoniumIngested'] + vars['Nitrate'] + vars['NitrateIngested']) / vars['Carbon'])
  Q_s = ((vars['Silicate'] + vars['SilicateIngested']) / vars['Carbon'])
  P_max_c = (((param['P_ref_c'] * T_function)) if (Q_N > param['Q_Nmax']) else (((0.0) if (Q_N < param['Q_Nmin']) else ((param['P_ref_c'] * T_function * ((Q_N - param['Q_Nmin']) / (param['Q_Nmax'] - param['Q_Nmin'])))))))
  Theta_c = (vars['Chlorophyll'] / vars['Carbon'])
  E_0 = (4.6 * env['Irradiance'])
  P_phot_c = ((0.0) if ((P_max_c == 0.0) or (Q_s <= param['Q_S_min'])) else ((P_max_c * (1.0 - math.exp(((-3600.0 * param['Alpha_Chl'] * Theta_c * E_0) / P_max_c))))))

  ### Chlorophyll Synthesis ###
  Theta_N = (vars['Chlorophyll'] / (vars['Ammonium'] + vars['AmmoniumIngested'] + vars['Nitrate'] + vars['NitrateIngested']))
  Rho_Chl = (((param['Theta_max_N'] * (P_phot_c / (3600.0 * param['Alpha_Chl'] * Theta_c * E_0)))) if ((E_0 > 0.0) and (Theta_c > 0.0)) else (0.0))

  ### Respiration ###
  R_C_growth = (((vars['AmmoniumIngested'] + vars['NitrateIngested']) * param['Zeta']) / (dt_in_hours * vars['Carbon']))
  R_C = (param['R_maintenance'] + R_C_growth)

  ### Cell Division ###
  C_d = ((2.0) if (((vars['Carbon'] + (vars['Carbon'] * (P_phot_c - R_C) * dt_in_hours)) >= param['C_rep']) and ((vars['Silicate'] + vars['SilicateIngested']) >= param['S_rep'])) else (1.0))
  if (C_d == 2.0):
    vars['Size'] = vars['Size'] * 2.0

  ### Nutrients uptake ###
  Q_nitrate = ((vars['Nitrate'] + vars['NitrateIngested']) / vars['Carbon'])
  Q_ammonium = ((vars['Ammonium'] + vars['AmmoniumIngested']) / vars['Carbon'])
  omega = ((param['k_AR'] / (param['k_AR'] + env['DissolvedAmmonium'])) * ((param['k_AR'] + env['DissolvedNitrate']) / (param['k_AR'] + env['DissolvedAmmonium'] + env['DissolvedNitrate'])))
  V_max_C = (((((param['V_ref_c'] * T_function)) if ((Q_ammonium + Q_nitrate) < param['Q_Nmin']) else (((0.0) if ((Q_ammonium + Q_nitrate) > param['Q_Nmax']) else ((param['V_ref_c'] * math.pow(((param['Q_Nmax'] - (Q_ammonium + Q_nitrate)) / (param['Q_Nmax'] - param['Q_Nmin'])), 0.05) * T_function)))))) if ((vars['Ammonium'] + vars['Nitrate']) < 1000.0) else (0.0))
  V_C_ammonium = (V_max_C * (env['DissolvedAmmonium'] / (param['k_AR'] + env['DissolvedAmmonium'])))
  V_C_nitrate = (V_max_C * (env['DissolvedNitrate'] / (param['k_AR'] + env['DissolvedNitrate'])) * omega)
  V_S_max = (((((param['V_S_ref'] * T_function)) if (Q_s <= param['Q_S_min']) else (((0.0) if (Q_s >= param['Q_S_max']) else ((param['V_S_ref'] * math.pow(((param['Q_S_max'] - Q_s) / (param['Q_S_max'] - param['Q_S_min'])), 0.05) * T_function)))))) if (vars['Carbon'] >= param['C_minS']) else (0.0))
  V_S_S = (V_S_max * (env['DissolvedSilicate'] / (env['DissolvedSilicate'] + param['k_S'])))
  vars['AmmoniumUptake'] = (vars['Carbon'] * V_C_ammonium * dt_in_hours)
  vars['NitrateUptake'] = (vars['Carbon'] * V_C_nitrate * dt_in_hours)
  vars['SilicateUptake'] = (vars['Silicate'] * V_S_S * dt_in_hours)

  ### Update Pools ###
  Ammonium_new = ((((vars['Ammonium'] + vars['AmmoniumIngested'] + vars['NitrateIngested']) - (vars['Ammonium'] * param['R_N'] * dt_in_hours * T_function)) - (((0.0 * (Q_N - param['Q_Nmax']))) if (Q_N > param['Q_Nmax']) else (0.0))) / C_d)
  Nitrate_new = 0.0
  Silicate_new = (((vars['Silicate'] + vars['SilicateIngested']) - 0.0) / C_d)
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
  vars['Silicate'] = Silicate_new
  vars['Carbon'] = Carbon_new
  vars['Chlorophyll'] = Chlorophyll_new
  vars['Nitrogen'] = Nitrogen_new
  vars['C_fuel'] = C_fuel_new
  vars['ChltoC'] = ChltoC_new
  vars['NtoC'] = NtoC_new

def update_Dead_Diatom(param, vars, env, dt):
  """ FGroup:  Diatom
      Stage:   Dead
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
  vars['Silicate'] = Silicate_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new

# Parameters for FGroup Copepod
# Species: Default_Copepod_Variety
species_Default_Copepod_Variety = {
    'A_rep' : 480.0,
    'A_rmax' : 960.0,
    'C1_min' : 6.25e-05,
    'C2_min' : 9.2e-05,
    'C3_min' : 0.00021,
    'C4_min' : 0.00058,
    'C5_min' : 0.00125,
    'C6_min' : 0.00333,
    'C_Cal' : 20.3,
    'C_conv1' : 12000.0,
    'E_m' : 0.25,
    'E_mech' : 0.3,
    'G_max' : 0.00833,
    'G_min' : 1e-05,
    'N4_min' : 1.7e-05,
    'N5_min' : 2.5e-05,
    'N6_min' : 3.75e-05,
    'N_mp' : 0.9,
    'OW_lipid' : 0.00633,
    'PreOW4' : 0.00058,
    'PreOW5' : 0.00125,
    'QR_10' : 3.4,
    'Q_Nmax' : 0.23,
    'QnProt' : 1.0,
    'S_max' : 0.0013,
    'T_ref' : 10.0,
    'V_max' : 45.0,
    'V_mconv1' : 0.0278,
    'Vol_conv1' : 1e-09,
    'a' : 1.584,
    'b' : 1.584,
    'delta' : 0.15,
    'k' : 85.2,
    'mi' : 0.0119,
    'n' : 0.8,
    'r_bas' : 0.000417,
    'r_sda' : 0.17,
    't_max' : 1.08,
    't_min' : 0.58,
    'vPrey' : 4.2e-09,
    'vol_gut' : 1.5e-08,
    'z_startOW' : 400.0,
}
# Foodset: Default_Copepod_Variety_P
foodset_Default_Copepod_Variety_P = {
    'Living' : {
        'P_min' : 100000.0,    },
}

def update_N3_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   N3
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 1 ###
  gamma = 0.5

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Molt N3 ###
  if (vars['C_N'] >= param['N4_min']):
    vars['Stage'] = stage_id('Copepod', 'N4')

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_N4_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   N4
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 1 ###
  gamma = 0.5

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Molt N4 ###
  if (vars['C_N'] >= param['N5_min']):
    vars['Stage'] = stage_id('Copepod', 'N5')

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_N5_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   N5
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 1 ###
  gamma = 0.5

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Molt N5 ###
  if (vars['C_N'] >= param['N6_min']):
    vars['Stage'] = stage_id('Copepod', 'N6')

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_N6_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   N6
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 1 ###
  gamma = 0.5

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Molt N6 ###
  if (vars['C_N'] >= param['C1_min']):
    vars['Stage'] = stage_id('Copepod', 'C1')

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_C1_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   C1
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 1 ###
  gamma = 0.5

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Molt C1 ###
  if (vars['C_N'] >= param['C2_min']):
    vars['Stage'] = stage_id('Copepod', 'C2')

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_C2_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   C2
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 1 ###
  gamma = 0.5

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Molt C2 ###
  if (vars['C_N'] >= param['C3_min']):
    vars['Stage'] = stage_id('Copepod', 'C3')

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_C3_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   C3
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Pre_overwintering-OW4 ###
  if (vars['C_N'] >= param['PreOW4']):
    new_agent_vars = {}
    new_agent_vars.update(vars)
    new_agent_vars['Stage'] = stage_id('Copepod', 'POW4')
    new_agent_vars['Size'] = vars['Size'] * ((0.3) if (param['d_year'] <= 210.0) else (0.5))
    vars['Size'] = vars['Size'] - new_agent_vars['Size']
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
    
  if (vars['C_N'] >= param['C4_min']):
    vars['Stage'] = stage_id('Copepod', 'C4')

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 1 ###
  gamma = 0.5

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_POW4_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   POW4
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Start_overwintering POW4 ###
  if (vars['C_NN'] >= param['OW_lipid']):
    vars['Stage'] = stage_id('Copepod', 'OWD4')

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Turn POW4 into C4OW ###
  if (param['d_year'] > 351.0):
    vars['Stage'] = stage_id('Copepod', 'C4OW')

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 3 ###
  gamma = 1.0

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_POW5_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   POW5
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Start_overwintering POW5 ###
  if (vars['C_NN'] >= param['OW_lipid']):
    vars['Stage'] = stage_id('Copepod', 'OWD5')

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Turn POW5 into C5 ###
  if (param['d_year'] > 351.0):
    vars['Stage'] = stage_id('Copepod', 'C5')

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 3 ###
  gamma = 1.0

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_OWD4_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   OWD4
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Start Overwintering_descent OWD4 ###
  V_m_new = (1.0 * dt_in_hours)
  if (vars['z'] >= param['z_startOW']):
    vars['Stage'] = stage_id('Copepod', 'OW4')

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 3 ###
  gamma = 1.0

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['V_m'] = V_m_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_OWD5_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   OWD5
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Start Overwintering_descent OWD5 ###
  V_m_new = (1.0 * dt_in_hours)
  if (vars['z'] >= param['z_startOW']):
    vars['Stage'] = stage_id('Copepod', 'OW5')

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 3 ###
  gamma = 1.0

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['V_m'] = V_m_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_OW4_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   OW4
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Overwintering motion ###
  V_m_new = 0.0

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Overwintering phase OW4 ###
  R_ow = (R_bas * param['delta'])
  if (param['d_year'] == 75.0):
    vars['Stage'] = stage_id('Copepod', 'OWA4')

  ### Assimilation efficiency ###
  Gut_time = 0.0
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Energetics during OW ###
  A_C = 0.0
  growth = A_C
  respiration = R_ow
  Growth_net = (growth - respiration)

  ### CNN update non-repro ###
  gamma = 0.0
  alpha = 0.0
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  A_PelletLoss = 0.0
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['V_m'] = V_m_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new

def update_OW5_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   OW5
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Overwintering motion ###
  V_m_new = 0.0

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Overwintering phase OW5 ###
  R_ow = (R_bas * param['delta'])
  if (param['d_year'] == 75.0):
    vars['Stage'] = stage_id('Copepod', 'OWA5')

  ### Assimilation efficiency ###
  Gut_time = 0.0
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Energetics during OW ###
  A_C = 0.0
  growth = A_C
  respiration = R_ow
  Growth_net = (growth - respiration)

  ### CNN update non-repro ###
  gamma = 0.0
  alpha = 0.0
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  A_PelletLoss = 0.0
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['V_m'] = V_m_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new

def update_OWA4_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   OWA4
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Post-overwintering ascent OWA4 ###
  V_m_new = -((1.0 * dt_in_hours))
  if ((vars['z'] <= param['Max_MLD']) or (I_t <= env['Irradiance'])):
    vars['Stage'] = stage_id('Copepod', 'C4OW')

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Ontogenetic fraction of C allocated to storage 2 ###
  gamma = 0.7

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_OWA5_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   OWA5
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Post-overwintering ascent OWA5 ###
  V_m_new = -((1.0 * dt_in_hours))
  if ((vars['z'] <= param['Max_MLD']) or (I_t <= env['Irradiance'])):
    vars['Stage'] = stage_id('Copepod', 'C5')

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Ontogenetic fraction of C allocated to storage 2 ###
  gamma = 0.7

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_C4_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   C4
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Pre_overwintering-OW5 ###
  if (vars['C_N'] >= param['PreOW5']):
    new_agent_vars = {}
    new_agent_vars.update(vars)
    new_agent_vars['Stage'] = stage_id('Copepod', 'POW5')
    new_agent_vars['Size'] = vars['Size'] * ((0.3) if (param['d_year'] <= 210.0) else (0.5))
    vars['Size'] = vars['Size'] - new_agent_vars['Size']
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
    
  if (vars['C_N'] >= param['C5_min']):
    vars['Stage'] = stage_id('Copepod', 'C5')

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 2 ###
  gamma = 0.7

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_C4OW_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   C4OW
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Pre_overwintering-C4OW ###
  if (vars['C_N'] >= param['C5_min']):
    vars['Stage'] = stage_id('Copepod', 'C5')

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 2 ###
  gamma = 0.7

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_C5_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   C5
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 2 ###
  gamma = 0.7

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Molt C5 ###
  if (vars['C_N'] >= param['C6_min']):
    vars['Stage'] = stage_id('Copepod', 'C6')

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_C6_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   C6
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Become adult ###
  if (vars['C_N'] >= param['G_max']):
    vars['Stage'] = stage_id('Copepod', 'Adult')

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 2 ###
  gamma = 0.7

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_Adult_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   Adult
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Become mature ###
  A_r_new = (vars['A_r'] + dt_in_hours)
  if ((vars['A_r'] == param['A_rep']) and (vars['C_N'] >= param['G_max'])):
    vars['Stage'] = stage_id('Copepod', 'Mature')

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 2 ###
  gamma = 0.7

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['A_r'] = A_r_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new

def update_Nauplius_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   Nauplius
  """
  dt_in_hours = dt / 3600.0

  ### Naupliar mortality counter ###
  Nauplius_counter_new = (vars['Nauplius_counter'] + dt_in_hours)

  ### Naupliar mortality ###
  if (vars['Nauplius_counter'] >= 1.0):
    new_agent_vars = {}
    new_agent_vars.update(vars)
    new_agent_vars['Stage'] = stage_id('Copepod', 'Dead')
    new_agent_vars['Size'] = vars['Size'] * param['N_mp']
    vars['Size'] = vars['Size'] - new_agent_vars['Size']
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
    
  if (vars['Nauplius_counter'] >= 1.0):
    vars['Stage'] = stage_id('Copepod', 'N3')

  ### Setting pool variables
  vars['Nauplius_counter'] = Nauplius_counter_new

def update_Mature_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   Mature
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  alpha = 0.0
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Reproduction ###
  A_r_new = (vars['A_r'] + dt_in_hours)
  Reproduce = 0.0
  if ((vars['A_r'] >= param['A_rep']) and (vars['C_N'] >= param['G_max'])):
    Reproduce = 1.0
  Nauplii = (((((((vars['C_N'] - param['G_max']) / param['G_min']) * 2.0)) if (vars['C_NN'] > (vars['C_N'] - param['G_max'])) else ((2.0 * (vars['C_NN'] / param['G_min']))))) if (Reproduce == 1.0) else (0.0))
  Q_AN = (vars['Ammonium'] / (vars['Ammonium'] + vars['Nitrate']))
  if (Reproduce == 1.0):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Nauplius')
    new_agent_vars['Size'] = min(Nauplii, 800.0)
    new_agent_vars['Carbon'] = param['G_min']
    new_agent_vars['A_r'] = 0.0
    new_agent_vars['Gut_content'] = 0.0
    new_agent_vars['Ammonium'] = (param['G_min'] * Q_N)
    new_agent_vars['Nitrate'] = (param['G_min'] * Q_N * (1.0 - Q_AN))
    new_agent_vars['C_N'] = ((param['G_min'] * (1.0 - 0.05)) / 2.0)
    new_agent_vars['C_NN'] = (((1.0 - 0.05) * param['G_min']) / 2.0)
    new_agent_vars['C_shell'] = (param['G_min'] * 0.05)
    new_agent_vars['C_pmax'] = ((param['G_min'] * (1.0 - 0.05)) / 2.0)
    new_agent_vars['V_gut'] = 4.0e-6
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  if (Reproduce == 1.0):
    vars['Stage'] = stage_id('Copepod', 'Senescent')

  ### Eggs production ###
  C_NN_new = (((vars['C_NN'] - min((vars['C_N'] - param['G_max']), vars['C_NN']))) if (Reproduce == 1.0) else (vars['C_NN']))
  C_N_new = (((vars['C_N'] - min((vars['C_N'] - param['G_max']), vars['C_NN']))) if (Reproduce == 1.0) else (vars['C_N']))
  Ammonium_new = ((((vars['Ammonium'] - ((Q_N * param['G_min'] * min(Nauplii, 800.0)) + NProt_excess + Cprot + A_PelletLoss)) + vars['AmmoniumIngested'] + A_Nitrate)) if (Reproduce == 1.0) else (vars['Ammonium']))

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new
  vars['A_r'] = A_r_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['Ammonium'] = Ammonium_new

def update_Senescent_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   Senescent
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Effect of size and temperature on swimming speed ###
  W_z = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * min((S / param['S_max']), 1.0))

  ### Day-time motion ###
  I_t = ((2.0 - vars['Gut_f']) * min((param['S_max'] / S), 1.0))

  ### Night-time motion ###
  Dlocal = sum(env['P'].values())
  Kn_calc2 = (0.4 * (2.0 - vars['Gut_f']))
  k_v_night = ((((((-((vars['Direction_1'] * Kn_calc2))) if (Dlocal < vars['Dlocal_previous']) else ((vars['Direction_1'] * Kn_calc2)))) if (vars['z'] < 250.0) else (-1.0))) if (vars['z'] > param['MLDepth']) else (0.0))
  Direction_new = ((1.0) if (k_v_night > 0.0) else (-1.0))
  Dlocal_previous_new = Dlocal

  ### Swimming direction ###
  kd_calc = (0.4 * (env['Irradiance'] - I_t))
  k_v_day = ((-1.0) if (kd_calc < -1.0) else (((1.0) if (kd_calc >= 1.0) else (kd_calc))))
  k_v = ((k_v_day) if (param['surface_irradiance'] > 0.0) else (k_v_night))

  ### Motion ###
  V_m_new = (k_v * param['V_max'] * W_z * dt_in_hours)

  ### Ingestion ###
  V_gut_new = (param['vol_gut'] * L)
  I_gCells = sum(vars['PIngestedCells'].values())
  Clock_new = (((vars['Clock'] + 1.0)) if (vars['Clock'] < 48.0) else (0.0))
  Prey_vol = (param['vPrey'] * I_gCells)
  Prey_VolDaily_new = (((vars['Prey_VolDaily'] + Prey_vol)) if (vars['Clock'] < 48.0) else (0.0))

  ### Basal metabolic cost ###
  R_bas = (param['r_bas'] * math.pow(vars['C_N'], 0.8) * math.pow(param['QR_10'], ((env['Temperature'] - param['T_ref']) / 10.0)))

  ### Swimming cost ###
  P_swim = (((param['k'] / 2.0) * math.pow(env['Density'], (1.0 - param['n'])) * math.pow((L / 10000.0), -(param['n'])) * math.pow((abs( (vars['V_m'] / dt_in_hours) ) * param['V_mconv1']), (3.0 - param['n'])) * math.pow(param['mi'], param['n']) * S) / 1.0e7)
  Z_swim = (P_swim / (param['E_mech'] * param['E_m']))
  O_cons = (((Z_swim / 1000.0) * 3600.0) / (param['C_Cal'] / 1000.0))
  R_swim = (O_cons * (12.0 / 22.4) * 1000.0 * 8.33e-05)

  ### Gut content ###
  Gut_contPlusPrey = (vars['Gut_content'] + Prey_vol)
  Gut_time = ((((param['t_min'] * param['t_max']) / (((Gut_contPlusPrey / vars['V_gut']) * (param['t_max'] - param['t_min'])) + param['t_min']))) if (vars['V_gut'] > 0.0) else (param['t_max']))
  Gut_clear = (((Gut_contPlusPrey / Gut_time)) if (Gut_time > 0.0) else (0.0))
  k_C = (1.0 - math.exp(-((param['b'] * Gut_time))))
  E = ((1.0 - k_C) * Gut_clear)
  A = (k_C * Gut_clear)
  A_C = (k_C * vars['CarbonIngested'])
  E_C = ((1.0 - k_C) * vars['CarbonIngested'])
  Gut_contTemp = max(0.0, (Gut_contPlusPrey - ((A + E) * dt_in_hours)))
  Gut_ftemp = ((0.0) if ((Gut_contTemp == 0.0) and (vars['V_gut'] == 0.0)) else (math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)))
  Gut_f_new = Gut_ftemp
  I_max = (((0.67 * vars['V_gut']) - Gut_contTemp) / (param['vPrey'] * 1800.0))
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min((((math.pi * math.pow((L * 2.9e-5), 2.0) * 1.0 * env['P'][variety] * 1.0e-6 * (1.0 - math.pow((Gut_contTemp / (0.67 * vars['V_gut'])), 2.0)) * (1.0 - math.exp((-1.7e-8 * env['P'][variety]))))) if (vars['V_gut'] > 0.0) else (I_max)), I_max)
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_min']
    

  ### Assimilation efficiency ###
  k_N = (1.0 - math.exp(-((param['a'] * Gut_time))))

  ### Egestion ###
  E_N = ((1.0 - k_N) * (vars['AmmoniumIngested'] + vars['NitrateIngested']))
  E_Si = vars['SilicateIngested']

  ### Faecal pellet ###
  PV_egest = (1.4e-6 * (vars['Carbon'] / param['G_max']))
  PV_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['PV'] + (E * dt_in_hours))))
  P_amm_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['P_amm'] + E_N)))
  Pc_new = ((0.0) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else ((vars['Pc'] + E_C)))
  A_PelletLoss = (((vars['P_amm'] + E_N)) if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest) else (0.0))
  if ((vars['PV'] + (E * dt_in_hours)) >= PV_egest):
    new_agent_vars = {}
    new_agent_vars['Stage'] = stage_id('Copepod', 'Pellet')
    new_agent_vars['Size'] = 1.0
    new_agent_vars['Ammonium'] = ((vars['P_amm'] + E_N))
    new_agent_vars['PV'] = (vars['PV'] + (E * dt_in_hours))
    new_agent_vars['Carbon'] = (vars['Pc'] + E_C)
    add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  vars['SilicateRelease'] = E_Si

  ### Assimilation ###
  A_Ammonium = (k_N * vars['AmmoniumIngested'])
  A_Nitrate = (k_N * vars['NitrateIngested'])

  ### Specific Dynamic Action cost ###
  R_sda = (param['r_sda'] * A_C)

  ### Energetics ###
  growth = A_C
  respiration = (R_bas + R_sda + R_swim)
  Growth_net = (growth - respiration)

  ### Update Gut content ###
  Gut_content_new = Gut_contTemp

  ### Ontogenetic fraction of C allocated to storage 2 ###
  gamma = 0.7

  ### Fraction allocated to carapace ###
  alpha = ((0.05) if (Growth_net > 0.0) else (0.0))

  ### CNN update non-repro ###
  C_NN_new = (((vars['C_NN'] + (gamma * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else ((((vars['C_NN'] + (Growth_net * dt_in_hours))) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else (vars['C_NN']))))

  ### Lipids pool ###
  Q_N = ((vars['Ammonium'] + (A_Ammonium * dt_in_hours)) / (vars['Carbon'] + (Growth_net * dt_in_hours)))
  NProt_excess = ((((vars['Carbon'] + (Growth_net * dt_in_hours)) * (Q_N - param['Q_Nmax']))) if ((vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) and (Q_N > param['Q_Nmax'])) else (0.0))

  ### C_N Update non-repro ###
  C_N_new = (((vars['C_N'] + ((1.0 - gamma) * (1.0 - alpha) * Growth_net * dt_in_hours))) if (Growth_net >= 0.0) else (((vars['C_N']) if (vars['C_NN'] >= (abs( Growth_net ) * dt_in_hours)) else ((vars['C_N'] + (Growth_net * dt_in_hours))))))

  ### Protein Pool ###
  Cprot = (((param['QnProt'] * abs( Growth_net ) * dt_in_hours)) if ((Growth_net < 0.0) and (vars['C_NN'] < (abs( Growth_net ) * dt_in_hours))) else (0.0))

  ### Carapace C Pool ###
  C_shell_new = (vars['C_shell'] + (((Growth_net * alpha * dt_in_hours)) if (Growth_net > 0.0) else (0.0)))

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Ammonium Pool ###
  Ammonium_new = ((vars['Ammonium'] + vars['AmmoniumIngested'] + A_Nitrate) - (A_PelletLoss + NProt_excess + Cprot))

  ### Nitrate Pool ###
  Nitrate_new = 0.0

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Silicon Pool ###
  Silicate_new = 0.0

  ### Mortality due to starvation ###
  if (vars['C_N'] <= (vars['C_pmax'] / 2.0)):
    vars['Stage'] = stage_id('Copepod', 'Dead')

  ### Mortality due to senescence ###
  A_r_new = (vars['A_r'] + dt_in_hours)
  new_agent_vars = {}
  new_agent_vars.update(vars)
  new_agent_vars['Stage'] = stage_id('Copepod', 'Dead')
  new_agent_vars['Size'] = vars['Size'] * (((1.0 / (param['A_rmax'] - vars['A_r']))) if (vars['A_r'] < param['A_rmax']) else (1.0))
  vars['Size'] = vars['Size'] - new_agent_vars['Size']
  add_agent('Copepod', new_agent_vars, [vars['x'], vars['y'], -vars['z']])
  

  ### Excretion ###
  C = (NProt_excess + Cprot)
  vars['AmmoniumRelease'] = C

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['Direction'] = Direction_new
  vars['Dlocal_previous'] = Dlocal_previous_new
  vars['V_m'] = V_m_new
  vars['V_gut'] = V_gut_new
  vars['Clock'] = Clock_new
  vars['Prey_VolDaily'] = Prey_VolDaily_new
  vars['Gut_f'] = Gut_f_new
  vars['PV'] = PV_new
  vars['P_amm'] = P_amm_new
  vars['Pc'] = Pc_new
  vars['Gut_content'] = Gut_content_new
  vars['C_NN'] = C_NN_new
  vars['C_N'] = C_N_new
  vars['C_shell'] = C_shell_new
  vars['Carbon'] = Carbon_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Silicate'] = Silicate_new
  vars['A_r'] = A_r_new

def update_Pellet_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   Pellet
  """
  dt_in_hours = dt / 3600.0

  ### Pellet sinking ###
  SRpellet = (math.pow(10.0, ((0.698 * numpy.log10((vars['PV'] * 1e+12))) - 2.03)) / 48.0)
  V_m_new = SRpellet

  ### Remineralisation ###
  R_nT = (0.0042 * math.pow(2.95, (((env['Temperature'] + 273.0) - 283.0) / 10.0)))
  vars['AmmoniumRelease'] = max(((vars['Ammonium'] + vars['Nitrate']) * R_nT * dt_in_hours), 0.0)
  Ammonium_new = max(((vars['Ammonium'] - (vars['Ammonium'] * R_nT * dt_in_hours)) + vars['AmmoniumIngested']), 0.0)
  Nitrate_new = max(((vars['Nitrate'] - (vars['Nitrate'] * R_nT * dt_in_hours)) + vars['NitrateIngested']), 0.0)
  vars['SilicateRelease'] = vars['SilicateIngested']

  ### Setting pool variables
  vars['V_m'] = V_m_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new

def update_Dead_Copepod(param, vars, env, dt):
  """ FGroup:  Copepod
      Stage:   Dead
  """
  dt_in_hours = dt / 3600.0

  ### Copepod size ###
  C_pmax_new = ((vars['C_N']) if (vars['C_N'] > vars['C_pmax']) else (vars['C_pmax']))

  ### Length and Surface Area ###
  L = math.pow(10.0, ((numpy.log10((vars['C_pmax'] * param['C_conv1'])) + 8.37) / 3.07))
  S = (L * 5.4e-07)

  ### Sinking ###
  V_m_new = (100.0 * (S / param['S_max']) * dt_in_hours)

  ### Total C ###
  Carbon_new = (vars['C_N'] + vars['C_NN'] + vars['C_shell'])

  ### Total N ###
  Nitrogen_new = (vars['Ammonium'] + vars['Nitrate'])

  ### Remineralisation ###
  R_nT = (0.0042 * math.pow(2.95, (((env['Temperature'] + 273.0) - 283.0) / 10.0)))
  vars['AmmoniumRelease'] = max(((vars['Ammonium'] + vars['Nitrate']) * R_nT * dt_in_hours), 0.0)
  Ammonium_new = max(((vars['Ammonium'] - (vars['Ammonium'] * R_nT * dt_in_hours)) + vars['AmmoniumIngested']), 0.0)
  Nitrate_new = max(((vars['Nitrate'] - (vars['Nitrate'] * R_nT * dt_in_hours)) + vars['NitrateIngested']), 0.0)
  vars['SilicateRelease'] = vars['SilicateIngested']

  ### Setting pool variables
  vars['C_pmax'] = C_pmax_new
  vars['V_m'] = V_m_new
  vars['Carbon'] = Carbon_new
  vars['Nitrogen'] = Nitrogen_new
  vars['Ammonium'] = Ammonium_new
  vars['Nitrate'] = Nitrate_new

# Parameters for FGroup Predator
# Species: Default_Predator_Variety
species_Default_Predator_Variety = {
    'S_max' : 15.0,
    'T_ref' : 1.0,
    'd_0' : 90.0,
}
# Foodset: Default_Predator_Variety_P
foodset_Default_Predator_Variety_P = {
    'OWD5' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 0.0006,
        'P_speed' : 24.2,
        'S_a' : 6.9e-08,
        'k_Iv' : 1000000.0,    },
    'OWA4' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 0.00021,
        'P_speed' : 18.9,
        'S_a' : 5.4e-08,
        'k_Iv' : 1000000.0,    },
    'C6' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 0.00333,
        'P_speed' : 33.4,
        'S_a' : 9.6e-08,
        'k_Iv' : 1000000.0,    },
    'Senescent' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 0.0083,
        'P_speed' : 45.0,
        'S_a' : 1.3e-07,
        'k_Iv' : 1000000.0,    },
    'OWA5' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 0.0006,
        'P_speed' : 24.2,
        'S_a' : 6.9e-08,
        'k_Iv' : 1000000.0,    },
    'POW4' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 0.00021,
        'P_speed' : 18.9,
        'S_a' : 5.4e-08,
        'k_Iv' : 1000000.0,    },
    'POW5' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 0.0006,
        'P_speed' : 24.2,
        'S_a' : 6.9e-08,
        'k_Iv' : 1000000.0,    },
    'Adult' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 0.0075,
        'P_speed' : 43.4,
        'S_a' : 1.2e-07,
        'k_Iv' : 1000000.0,    },
    'C5' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 0.00125,
        'P_speed' : 24.2,
        'S_a' : 6.9e-08,
        'k_Iv' : 1000000.0,    },
    'C4OW' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 0.00021,
        'P_speed' : 18.9,
        'S_a' : 5.4e-08,
        'k_Iv' : 1000000.0,    },
    'OWD4' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 0.00021,
        'P_speed' : 18.9,
        'S_a' : 5.4e-08,
        'k_Iv' : 1000000.0,    },
    'Mature' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 0.0083,
        'P_speed' : 45.0,
        'S_a' : 1.3e-07,
        'k_Iv' : 1000000.0,    },
    'C3' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 0.00021,
        'P_speed' : 13.5,
        'S_a' : 3.9e-08,
        'k_Iv' : 1000000.0,    },
    'C2' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 9.2e-05,
        'P_speed' : 10.3,
        'S_a' : 3e-08,
        'k_Iv' : 1000000.0,    },
    'C1' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 6.25e-05,
        'P_speed' : 9.1,
        'S_a' : 2.6e-08,
        'k_Iv' : 1000000.0,    },
    'N3' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 1e-05,
        'P_speed' : 5.1,
        'S_a' : 1.5e-08,
        'k_Iv' : 1000000.0,    },
    'N4' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 1.7e-05,
        'P_speed' : 5.9,
        'S_a' : 1.7e-08,
        'k_Iv' : 1000000.0,    },
    'N5' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 2.5e-05,
        'P_speed' : 6.8,
        'S_a' : 1.9e-08,
        'k_Iv' : 1000000.0,    },
    'N6' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 3.75e-05,
        'P_speed' : 7.7,
        'S_a' : 2.2e-08,
        'k_Iv' : 1000000.0,    },
    'C4' : {
        'K_p' : 0.001,
        'P_minv' : 1000.0,
        'P_size' : 0.0006,
        'P_speed' : 18.9,
        'S_a' : 5.4e-08,
        'k_Iv' : 1000000.0,    },
}

def update_Existance_Predator(param, vars, env, dt):
  """ FGroup:  Predator
      Stage:   Existance
  """
  dt_in_hours = dt / 3600.0

  ### Ingestion ###
  # ml805: Prevent log10(0.0)
  W_tg = ((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * ((1.0) if (param['S_t'] >= param['S_max']) else ((param['S_t'] / param['S_max']))))
  I_gv = {}
  if (param['S_t']>0.):
    G = (((2.37 * numpy.log10(param['S_t'])) - 1.22) / 12.0)
    for variety in env['P'].keys():
      I_gv[variety] = min((W_tg * (((param[variety]['K_p'] * (param[variety]['S_a'] / 1.3e-07) * (env['Irradiance'] / 1.0) * (math.pow((env['P'][variety] - param[variety]['P_minv']), 2.0) / ((env['P'][variety] - param[variety]['P_minv']) + param[variety]['k_Iv'])) * (45.0 / param[variety]['P_speed']))) if (env['P'][variety] > param[variety]['P_minv']) else (0.0))), ((G * 0.6156 * math.exp(-((0.0321 * (param['d_year'] - param['d_0']))))) / (86400.0 * param[variety]['P_size'])))
  else:
    for variety in env['P'].keys():
      I_gv[variety] = (W_tg * (((param[variety]['K_p'] * (param[variety]['S_a'] / 1.3e-07) * (env['Irradiance'] / 1.0) * (math.pow((env['P'][variety] - param[variety]['P_minv']), 2.0) / ((env['P'][variety] - param[variety]['P_minv']) + param[variety]['k_Iv'])) * (45.0 / param[variety]['P_speed']))) if (env['P'][variety] > param[variety]['P_minv']) else (0.0)))

  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_minv']
    

  ### Egestion ###
  new_agent_vars = {}
  new_agent_vars['Stage'] = stage_id('Predator', 'Pellet')
  new_agent_vars['Size'] = 1.0
  new_agent_vars['Ammonium'] = (vars['AmmoniumIngested'] + vars['NitrateIngested'])
  add_agent('Predator', new_agent_vars, [vars['x'], vars['y'], -vars['z']])

  ### Silicon Pool ###
  Silicate_new = 0.0
  vars['SilicateRelease'] = vars['SilicateIngested']

  ### Setting pool variables
  vars['Silicate'] = Silicate_new

def update_Pellet_Predator(param, vars, env, dt):
  """ FGroup:  Predator
      Stage:   Pellet
  """
  dt_in_hours = dt / 3600.0

  ### Pellet sinking ###
  V_m_new = 5.0

  ### Remineralisation ###
  R_nTPred = (0.0042 * math.pow(2.95, (((env['Temperature'] + 273.0) - 283.0) / 10.0)))
  vars['AmmoniumRelease'] = (((vars['Ammonium'] * R_nTPred * dt_in_hours)) if (vars['Ammonium'] > 0.0) else (0.0))
  Ammonium_new = max((vars['Ammonium'] - (vars['Ammonium'] * R_nTPred * dt_in_hours)), 0.0)

  ### Silicon Pool ###
  Silicate_new = 0.0
  vars['SilicateRelease'] = vars['SilicateIngested']

  ### Setting pool variables
  vars['V_m'] = V_m_new
  vars['Ammonium'] = Ammonium_new
  vars['Silicate'] = Silicate_new

# Parameters for FGroup Basal_predator
# Species: Default_Basal_Predator_Variety
species_Default_Basal_Predator_Variety = {
    'I_max40' : 6e-05,
    'T_ref' : 1.0,
}
# Foodset: Default_Basal_Predator_Variety_P
foodset_Default_Basal_Predator_Variety_P = {
    'OWD5' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 0.0006,
        'k_Iv' : 1000000.0,    },
    'OWA4' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 0.00021,
        'k_Iv' : 1000000.0,    },
    'C6' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 0.00333,
        'k_Iv' : 1000000.0,    },
    'Senescent' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 0.0083,
        'k_Iv' : 1000000.0,    },
    'OWA5' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 0.0006,
        'k_Iv' : 1000000.0,    },
    'POW4' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 0.00021,
        'k_Iv' : 1000000.0,    },
    'POW5' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 0.0006,
        'k_Iv' : 1000000.0,    },
    'Adult' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 0.0075,
        'k_Iv' : 1000000.0,    },
    'C5' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 0.00125,
        'k_Iv' : 1000000.0,    },
    'C4OW' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 0.00021,
        'k_Iv' : 1000000.0,    },
    'OWD4' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 0.00021,
        'k_Iv' : 1000000.0,    },
    'Mature' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 0.0083,
        'k_Iv' : 1000000.0,    },
    'C3' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 0.00021,
        'k_Iv' : 1000000.0,    },
    'C2' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 9.2e-05,
        'k_Iv' : 1000000.0,    },
    'C1' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 6.25e-05,
        'k_Iv' : 1000000.0,    },
    'N3' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 5e-06,
        'k_Iv' : 1000000.0,    },
    'N4' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 1.7e-05,
        'k_Iv' : 1000000.0,    },
    'N5' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 2.5e-05,
        'k_Iv' : 1000000.0,    },
    'N6' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 3.75e-05,
        'k_Iv' : 1000000.0,    },
    'C4' : {
        'K_p' : 0.0001,
        'P_minv' : 1000.0,
        'P_size' : 0.0006,
        'k_Iv' : 1000000.0,    },
}

def update_Existance_Basal_predator(param, vars, env, dt):
  """ FGroup:  Basal_predator
      Stage:   Existance
  """
  dt_in_hours = dt / 3600.0

  ### Ingestion ###
  I_gv = {}
  for variety in env['P'].keys():
    I_gv[variety] = min(((0.3 + (0.7 * (env['Temperature'] / param['T_ref']))) * ((((env['P'][variety] - param[variety]['P_minv']) * ((env['P'][variety] - param[variety]['P_minv']) / ((env['P'][variety] - param[variety]['P_minv']) + param[variety]['k_Iv'])) * param[variety]['K_p'])) if (env['P'][variety] > param[variety]['P_minv']) else (0.0))), (param['I_max40'] / param[variety]['P_size']))
  for variety in vars['PRequest'].keys():
    vars['PRequest'][variety] = dt * I_gv[variety]
    vars['PThreshold'][variety] = param[variety]['P_minv']
    

  ### Egestion ###
  new_agent_vars = {}
  new_agent_vars['Stage'] = stage_id('Basal_predator', 'Pellet')
  new_agent_vars['Size'] = 1.0
  new_agent_vars['Ammonium'] = (vars['AmmoniumIngested'] + vars['NitrateIngested'])
  add_agent('Basal_predator', new_agent_vars, [vars['x'], vars['y'], -vars['z']])

  ### Silicon Pool ###
  Silicate_new = 0.0
  vars['SilicateRelease'] = vars['SilicateIngested']

  ### Setting pool variables
  vars['Silicate'] = Silicate_new

def update_Pellet_Basal_predator(param, vars, env, dt):
  """ FGroup:  Basal_predator
      Stage:   Pellet
  """
  dt_in_hours = dt / 3600.0

  ### Pellet sinking ###
  V_m_new = 5.0

  ### Remineralisation ###
  R_nTBP = (0.0042 * math.pow(2.95, (((env['Temperature'] + 273.0) - 283.0) / 10.0)))
  vars['AmmoniumRelease'] = (((vars['Ammonium'] * R_nTBP * dt_in_hours)) if (vars['Ammonium'] > 0.0) else (0.0))
  Ammonium_new = max((vars['Ammonium'] - (vars['Ammonium'] * R_nTBP * dt_in_hours)), 0.0)

  ### Silicon Pool ###
  Silicate_new = 0.0
  vars['SilicateRelease'] = vars['SilicateIngested']

  ### Setting pool variables
  vars['V_m'] = V_m_new
  vars['Ammonium'] = Ammonium_new
  vars['Silicate'] = Silicate_new
