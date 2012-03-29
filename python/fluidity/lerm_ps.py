import math

# Parameters for FG Diatom
A_E = -10000.0
Alpha_Chl = 7.9e-07
C_minS = 1.584e-08
C_rep = 1.76e-08
C_starve = 6.24e-09
C_struct = 8.5e-08
Ndis = 0.0042
P_ref_c = 0.14
Q_Nmax = 0.17
Q_Nmin = 0.034
Q_S_max = 0.43
Q_S_min = 0.04
Q_remN = 2.95
Q_remS = 2.27
R_Chl = 0.00072
R_N = 0.0042
R_maintenance = 0.00072
S_dis = 0.00083
S_rep = 2.1e-09
T_ref = 293.0
T_refN = 283.0
T_refS = 278.0
Theta_max_N = 4.2
V_S_ref = 0.03
V_ref_c = 0.008
Zeta = 2.3
k_AR = 0.027
k_S = 1.0
z_sink = 0.041667

## FG: Diatom;   Stage: Living;   ID: 0
def update_Living_Diatom(vars, env, dt):
  dt_in_hours = dt / 3600.0

  ### Effect of temperature ###
  # Convert to Kelvin
  T_K = (env["Temperature"] + 273.0)
  # Temperature function
  T_function = math.exp((A_E * ((1.0 / T_K) - (1.0 / T_ref))))

  ### Photosynthesis ###
  # Ratio N to C
  Q_N = ((vars["Ammonium"] + vars["AmmoniumIngested"] + vars["Nitrate"] + vars["NitrateIngested"]) / vars["Carbon"])
  # Ratio Si to C
  Q_s = ((vars["Silicate"] + vars["SilicateIngested"]) / vars["Carbon"])
  # Maximum carbon sr
  P_max_c = (((P_ref_c * T_function)) if ((Q_N > Q_Nmax)) else (((0.0) if ((Q_N < Q_Nmin)) else ((P_ref_c * T_function * ((Q_N - Q_Nmin) / (Q_Nmax - Q_Nmin)))))))
  # Calculate Theta_c
  Theta_c = (vars["Chlorophyll"] / vars["Carbon"])
  # Visible irradiance in microE/s*m_2) units
  E_0 = (4.6 * env["Irradiance"])
  # Carbon specific rate of photosynthesis
  P_phot_c = ((0.0) if (((P_max_c == 0.0)) or ((Q_s <= Q_S_min))) else ((P_max_c * (1.0 - math.exp(((-3600.0 * Alpha_Chl * Theta_c * E_0) / P_max_c))))))

  ### Chlorophyll Synthesis ###
  # Chlorophyll a to Nitrogen ratio
  Theta_N = (vars["Chlorophyll"] / (vars["Ammonium"] + vars["AmmoniumIngested"] + vars["Nitrate"] + vars["NitrateIngested"]))
  # Chlorophyll a synthesis regeneration
  Rho_Chl = (((Theta_max_N * (P_phot_c / (3600.0 * Alpha_Chl * Theta_c * E_0)))) if (((E_0 > 0.0)) and ((Theta_c > 0.0))) else (0.0))

  ### Respiration ###
  # Carbon specific rate of growth related respiration
  R_C_growth = (((vars["AmmoniumIngested"] + vars["NitrateIngested"]) * Zeta) / (dt_in_hours * vars["Carbon"]))
  # Total respiration
  R_C = (R_maintenance + R_C_growth)

  ### Cell Division ###
  # Is the diatom ready to cell divide?
  C_d = ((2.0) if ((((vars["Carbon"] + (vars["Carbon"] * (P_phot_c - R_C) * dt_in_hours)) >= C_rep)) and (((vars["Silicate"] + vars["SilicateIngested"]) >= S_rep))) else (1.0))
  # Cell division
  if (C_d == 2.0):
    vars["Size"] = vars["Size"] * 2.0

  ### Nutrients uptake ###
  # nitrate:carbon ratio
  Q_nitrate = ((vars["Nitrate"] + vars["NitrateIngested"]) / vars["Carbon"])
  # ammonium:carbon ratio
  Q_ammonium = ((vars["Ammonium"] + vars["AmmoniumIngested"]) / vars["Carbon"])
  # inhibition factor for nitrate uptake
  omega = ((k_AR / (k_AR + env["DissolvedAmmonium"])) * ((k_AR + env["DissolvedNitrate"]) / (k_AR + env["DissolvedAmmonium"] + env["DissolvedNitrate"])))
  # maximum nitrogen uptake
  V_max_C = (((((V_ref_c * T_function)) if (((Q_ammonium + Q_nitrate) < Q_Nmin)) else (((0.0) if (((Q_ammonium + Q_nitrate) > Q_Nmax)) else ((V_ref_c * math.pow(((Q_Nmax - (Q_ammonium + Q_nitrate)) / (Q_Nmax - Q_Nmin)), 0.05) * T_function)))))) if (((vars["Ammonium"] + vars["Nitrate"]) < 1000.0)) else (0.0))
  # Ammonium uptake
  V_C_ammonium = (V_max_C * (env["DissolvedAmmonium"] / (k_AR + env["DissolvedAmmonium"])))
  # Nitrate uptake
  V_C_nitrate = (V_max_C * (env["DissolvedNitrate"] / (k_AR + env["DissolvedNitrate"])) * omega)
  # Silicate uptake
  V_S_max = (((((V_S_ref * T_function)) if ((Q_s <= Q_S_min)) else (((0.0) if ((Q_s >= Q_S_max)) else ((V_S_ref * math.pow(((Q_S_max - Q_s) / (Q_S_max - Q_S_min)), 0.05) * T_function)))))) if ((vars["Carbon"] >= C_minS)) else (0.0))
  # rate of Si uptake
  V_S_S = (V_S_max * (env["DissolvedSilicate"] / (env["DissolvedSilicate"] + k_S)))
  # Ammonium absorption
  vars["AmmoniumUptake"] = (vars["Carbon"] * V_C_ammonium * dt_in_hours)
  # Nitrate absorption
  vars["NitrateUptake"] = (vars["Carbon"] * V_C_nitrate * dt_in_hours)
  # Silicate absorption
  vars["SilicateUptake"] = (vars["Silicate"] * V_S_S * dt_in_hours)

  ### Update Pools ###
  # Update Ammonium
  Ammonium_new = ((((vars["Ammonium"] + vars["AmmoniumIngested"] + vars["NitrateIngested"]) - (vars["Ammonium"] * R_N * dt_in_hours * T_function)) - (((0.0 * (Q_N - Q_Nmax))) if ((Q_N > Q_Nmax)) else (0.0))) / C_d)
  # Update Nitrate
  Nitrate_new = 0.0
  # Update Silicon
  Silicate_new = (((vars["Silicate"] + vars["SilicateIngested"]) - 0.0) / C_d)
  # Calculate Carbon
  C_new = max(0.0, (((vars["Carbon"] * (P_phot_c - (R_C * T_function)) * dt_in_hours) + vars["Carbon"]) / C_d))
  # Check Mortality
  death_flag = ((1.0) if ((C_new <= C_starve)) else (0.0))
  # Death
  if (death_flag == 1.0):
    vars["Stage"] = 1.0  # Dead
  # Update Carbon
  Carbon_new = C_new
  # Update Cholorophyll
  Chlorophyll_new = ((max((((((vars["Chlorophyll"] + (Rho_Chl * (vars["AmmoniumIngested"] + vars["NitrateIngested"])))) if ((Theta_N <= Theta_max_N)) else ((vars["Chlorophyll"] - (vars["Chlorophyll"] - ((vars["Ammonium"] + vars["Nitrate"]) * Theta_max_N))))) - ((vars["Chlorophyll"] * R_Chl * dt_in_hours * T_function) + 0.0)) / C_d), 0.0)) if ((death_flag == 0.0)) else (0.0))
  # Update total Nitrogen
  Nitrogen_new = (vars["Ammonium"] + vars["Nitrate"])
  # Update C_fuel
  C_fuel_new = (vars["Carbon"] - C_struct)

  ### Remineralisation Nitrogen ###
  # Nitrogen
  vars["AmmoniumRelease"] = ((vars["Ammonium"] + vars["Nitrate"]) * R_N * dt_in_hours * T_function)

  ### Photoadaptation ###
  # Chl to C ratio
  ChltoC_new = (((vars["Chlorophyll"] / vars["Carbon"])) if ((vars["Carbon"] > 0.0)) else (0.0))
  # N to C ratio
  NtoC_new = ((((vars["Ammonium"] + vars["Nitrate"]) / vars["Carbon"])) if ((vars["Carbon"] > 0.0)) else (0.0))

  ### Setting pool variables
  vars["NtoC"] = NtoC_new
  vars["Chlorophyll"] = Chlorophyll_new
  vars["C_fuel"] = C_fuel_new
  vars["Nitrogen"] = Nitrogen_new
  vars["Carbon"] = Carbon_new
  vars["Nitrate"] = Nitrate_new
  vars["Silicate"] = Silicate_new
  vars["Ammonium"] = Ammonium_new
  vars["ChltoC"] = ChltoC_new

## FG: Diatom;   Stage: Dead;   ID: 1
def update_Dead_Diatom(vars, env, dt):
  dt_in_hours = dt / 3600.0

  ### Remineralisation Dead T ###
  # Si dissolution rate
  Si_reminT = (S_dis * math.pow(Q_remS, (((env["Temperature"] + 273.0) - T_refS) / 10.0)))
  # N dissolution rate
  N_reminT = (Ndis * math.pow(Q_remN, (((env["Temperature"] + 273.0) - T_refN) / 10.0)))
  # Silicon
  vars["SilicateRelease"] = (vars["Silicate"] * Si_reminT * dt_in_hours)
  # Update Sipool
  Silicate_new = max((vars["Silicate"] - (vars["Silicate"] * Si_reminT * dt_in_hours)), 0.0)
  # Nitrogen
  vars["AmmoniumRelease"] = ((vars["Ammonium"] + vars["Nitrate"]) * N_reminT * dt_in_hours)
  # Update Apool
  Ammonium_new = max((vars["Ammonium"] - (vars["Ammonium"] * N_reminT * dt_in_hours)), 0.0)
  # Update Npool
  Nitrate_new = max((vars["Nitrate"] - (vars["Nitrate"] * N_reminT * dt_in_hours)), 0.0)

  ### Setting pool variables
  vars["Nitrate"] = Nitrate_new
  vars["Silicate"] = Silicate_new
  vars["Ammonium"] = Ammonium_new
