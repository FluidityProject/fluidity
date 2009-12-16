
def pznd(state, parameters):
    '''Calculate sources and sinks for a simple PZND model'''
    
    if not check_pznd_parameters(parameters):
        raise TypeError("Missing Parameter")
    
    P=state.scalar_fields["Phytoplankton"]
    Z=state.scalar_fields["Zooplankton"]
    N=state.scalar_fields["Nutrient"]
    D=state.scalar_fields["Detritus"]
    I=state.scalar_fields["PhotosyntheticRadiation"]
    Pnew=state.scalar_fields["IteratedPhytoplankton"]
    Znew=state.scalar_fields["IteratedZooplankton"]
    Nnew=state.scalar_fields["IteratedNutrient"]
    Dnew=state.scalar_fields["IteratedDetritus"]


    P_source=state.scalar_fields["PhytoplanktonSource"]
    Z_source=state.scalar_fields["ZooplanktonSource"]
    N_source=state.scalar_fields["NutrientSource"]
    D_source=state.scalar_fields["DetritusSource"]
    try:
        PP=state.scalar_fields["PrimaryProduction"]
    except KeyError:
        PP=None
    try:
        PG=state.scalar_fields["PhytoplanktonGrazing"]
    except KeyError:
        PG=None
    

    alpha=parameters["alpha"]
    beta=parameters["beta"]
    gamma=parameters["gamma"]
    g=parameters["g"]
    k_N=parameters["k_N"]
    k=parameters["k"]
    v=parameters["v"]
    mu_P=parameters["mu_P"]
    mu_Z=parameters["mu_Z"]
    mu_D=parameters["mu_D"]
    p_P=parameters["p_P"]
    p_D=1-p_P

    for n in range(P.node_count):
        # Values of fields on this node.
        P_n=max(.5*(P.node_val(n)+Pnew.node_val(n)), 0.0)
        Z_n=max(.5*(Z.node_val(n)+Znew.node_val(n)), 0.0)
        N_n=max(.5*(N.node_val(n)+Nnew.node_val(n)), 0.0)
        D_n=max(.5*(D.node_val(n)+Dnew.node_val(n)), 0.0)
        I_n=max(I.node_val(n), 0.0)

        # Light limited phytoplankton growth rate.
        J=(v*alpha*I_n)/(v**2+alpha**2*I_n**2)**0.5

        # Nitrate limiting factor.
        Q=N_n/(k_N+N_n)

        # Total phytoplankton growth rate.
        P_r=J*P_n*Q

        # Zooplankton grazing of phytoplankton.
        G_P=(g * p_P * P_n**2 * Z_n)/(k**2 + p_P*P_n**2 + p_D*D_n**2)

        # Zooplankton grazing of detritus.
        G_D=(g * (1-p_P) * D_n**2 * Z_n)/(k**2 + p_P*P_n**2 + p_D*D_n**2)

        # Death rate of phytoplankton.
        De_P=mu_P*P_n*P_n/(P_n+1)

        # Death rate of zooplankton.
        De_Z=mu_Z*Z_n**2

        # Detritus remineralisation.
        De_D=mu_D*D_n

        P_source.addto(n, P_r - G_P - De_P)
        
        if PP:
            PP.set(n, P_r)
	if PG:
	    PG.set(n, G_P)

        Z_source.addto(n, gamma*beta*(G_P+G_D) - De_Z)
        
        N_source.addto(n, -P_r + De_D + (1-gamma)*beta*(G_P+G_D))

        D_source.addto(n, -De_D + De_P + De_Z +(1-beta)*G_P - beta*G_D)                     

def check_pznd_parameters(parameters):
    from sys import stderr

    valid=True

    if not parameters.has_key("alpha"):
        stderr.write("PZND parameter alpha missing.\n")
        stderr.write("alpha is this initial slope of the P-I curve.\n\n")
        valid = False

    if not parameters.has_key("beta"):
        stderr.write("PZND parameter beta missing.\n")
        stderr.write("beta is the assimilation efficiency of zooplankton.\n\n")
        valid = False

    if not parameters.has_key("gamma"):
        stderr.write("PZND parameter gamma missing.\n")
        stderr.write("gamma is the zooplankton excretion parameter.\n\n")
        valid = False

    if not parameters.has_key("g"):
        stderr.write("PZND parameter g missing.\n")
        stderr.write("g is the zooplankton maximum growth rate.\n\n")
        valid = False

    if not parameters.has_key("k_N"):
        stderr.write("PZND parameter k_N missing.\n")
        stderr.write("k_N is the half-saturation constant for nutrient.\n\n")
        valid = False

    if not parameters.has_key("k"):
        stderr.write("PZND parameter k missing.\n")
        stderr.write("k is the zooplankton grazing parameter.\n\n")
        valid = False

    if not parameters.has_key("mu_P"):
        stderr.write("PZND parameter mu_P missing.\n")
        stderr.write("mu_P is the phytoplankton mortality rate.\n\n")
        valid = False

    if not parameters.has_key("mu_Z"):
        stderr.write("PZND parameter mu_Z missing.\n")
        stderr.write("mu_Z is the zooplankton mortality rate.\n\n")
        valid = False

    if not parameters.has_key("mu_D"):
        stderr.write("PZND parameter mu_D missing.\n")
        stderr.write("mu_D is the detritus remineralisation rate.\n\n")
        valid = False

    if not parameters.has_key("p_P"):
        stderr.write("PZND parameter p_P missing.\n")
        stderr.write("p_P is the relative grazing preference of zooplankton for phytoplankton.\n\n")
        valid = False

    if not parameters.has_key("v"):
        stderr.write("PZND parameter v missing.\n")
        stderr.write("v is the maximum phytoplankton growth rate.\n\n")
        valid = False

    return valid
    
def lotka_volterra(state,parameters):

    if not check_lotka_volterra_parameters(parameters):
        raise TypeError("Missing Parameter")
    
    P=state.scalar_fields["Phytoplankton"]
    Z=state.scalar_fields["Zooplankton"]
    Pnew=state.scalar_fields["IteratedPhytoplankton"]
    Znew=state.scalar_fields["IteratedZooplankton"]
    
    P_source=state.scalar_fields["PhytoplanktonSource"]
    Z_source=state.scalar_fields["ZooplanktonSource"]

    alpha=parameters["alpha"]
    beta=parameters["beta"]
    gamma=parameters["gamma"]
    delta=parameters["delta"]

    for n in range(P.node_count):
        # Values of fields on this node.
        P_n=.5*(P.node_val(n)+Pnew.node_val(n))
        Z_n=.5*(Z.node_val(n)+Znew.node_val(n))
    
        P_source.set(n, P_n*(alpha-beta*Z_n))

        Z_source.set(n, -Z_n*(gamma-delta*P_n))


def check_lotka_volterra_parameters(parameters):
    from sys import stderr

    valid=True

    if not parameters.has_key("alpha"):
        stderr.write("Lotka Voltera parameter alpha missing.\n")
        valid = False

    if not parameters.has_key("beta"):
        stderr.write("Lotka Voltera parameter beta missing.\n")
        valid = False

    if not parameters.has_key("gamma"):
        stderr.write("Lotka Voltera parameter gamma missing.\n")
        valid = False

    if not parameters.has_key("delta"):
        stderr.write("Lotka Voltera parameter delta missing.\n")
        valid = False

    if not valid:
        stderr.write(" dP/dt = P*(alpha-beta * Z)")
        stderr.write(" dZ/dt = - Z*(gamma-delta * P)")

    return valid
