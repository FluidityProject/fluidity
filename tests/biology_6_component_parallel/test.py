import fluidity.ocean_biology as biology
import fluidity.state_types as state_types



per_day=1./(3600*24)

p={}
p["alpha_c"]=0.02*per_day
p["beta_p"]=0.75
p["beta_d"]=0.75
p["gamma"]=0.5
p["epsilon"]=3.3
p["g"]=1.3*per_day
p["k_A"]=0.5
p["k_N"]=0.5
p["mu_P"]=0.05*per_day
p["mu_Z"]=0.2*per_day
p["mu_D"]=0.05*per_day
p["p_P"]=0.75
p["v"]=1*per_day
p["theta_m"]=0.05
p["zeta"]=12.8
p["photic_zone_depth"]=116
p["psi"]=2.9

biology.pczdna(state, p)
