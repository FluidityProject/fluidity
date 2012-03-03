#include "LERM_Diatom_c.h"

/*****************************/
/* MAIN PLANKTON UPDATE CODE */
/*****************************/

void updateLivingDiatom_c(double vars[], int n_vars, double env[], int n_env, double *dt) { 

	/* Phase 1
	 * (import and initial housekeeping)
	 */

/*
    vars[C_OLD] = vars[C_NEW];
    vars[Z_OLD] = vars[Z_NEW];
    vars[AMMONIUM_ING] /= vars[C_OLD];
    vars[NITRATE_ING] /= vars[C_OLD];
    vars[SILICATE_ING] /= vars[C_OLD];
*/

    /* External environmnet */
    //int z_old_int = (int) (vars[Z_OLD]);
    //int pLayer = getPLayer(vars[Z_NEW]);


    /* Nutrient uptake */
/*
    if (constData->chem_depletion[z_old_int][AMMONIUM]<1)
        vars[AMMONIUM_ING] *= constData->chem_depletion[z_old_int][AMMONIUM];
    if (constData->chem_depletion[z_old_int][NITRATE]<1)
        vars[NITRATE_ING] *= constData->chem_depletion[z_old_int][NITRATE];
    if (constData->chem_depletion[z_old_int][SILICATE]<1)
        vars[SILICATE_ING] *= constData->chem_depletion[z_old_int][SILICATE];
*/

    /* Phase 2
     * This is the arithmetic meat and where most of the workload lies
     * (this gets auto-generated)
     */

    double stepInHours = *dt / 3600.;

    double _ambientTemperature = env[TEMPERATURE];
    double _ambientVisIrrad = env[IRRADIANCE];
    double _ambientAmmonium = env[AMMONIUM];
    double _ambientNitrate = env[NITRATE];
    double _ambientSilicate = env[SILICATE];

    double T_K = _ambientTemperature + 273.0;
    double T_function = exp(34.12969283 - 10000/T_K);

    /* Photosynthesis */
    double Q_N = ((((vars[AMMONIUM_POOL])+(vars[AMMONIUM_ING])+(vars[NITRATE_POOL])+(vars[NITRATE_ING]))) / (vars[CARBON_POOL]));
    double Q_s = ((((vars[SILICATE_POOL])+(vars[SILICATE_ING]))) / (vars[CARBON_POOL]));
    double P_max_c = (((Q_N) > (param_Q_Nmax))?(((param_P_ref_c)*(T_function))):(((Q_N) < (param_Q_Nmin))?(0.0):(((param_P_ref_c)*(T_function)*(((Q_N) - (param_Q_Nmin)) / ((param_Q_Nmax) - (param_Q_Nmin)))))));
    double Theta_c = ((vars[CHLOROPHYLL_POOL]) / (vars[CARBON_POOL]));
    double E_0 = (((4.6)*(_ambientVisIrrad)));

    double P_phot_c = ((P_max_c == 0.0)||(Q_s <= param_Q_S_min))?0.0:P_max_c*(1.0 - exp(-0.284400e-2*Theta_c*E_0 / P_max_c));

    /* Reproduction */
    double Theta_N = ((vars[CHLOROPHYLL_POOL]) / (((vars[AMMONIUM_POOL])+(vars[AMMONIUM_ING])+(vars[NITRATE_POOL])+(vars[NITRATE_ING]))));
    double Rho_Chl = (((((E_0) > (0.0))&&((Theta_c) > (0.0))))?(((param_Theta_max_N)*((P_phot_c) / (((3600.0)*(param_Alpha_Chl)*(Theta_c)*(E_0)))))):(0.0));
    double R_C_growth = ((((((vars[AMMONIUM_ING])+(vars[NITRATE_ING])))*(param_Zeta))) / (((stepInHours)*(vars[CARBON_POOL]))));
    double R_C = (((param_R_maintenance)+(R_C_growth)));
    double C_d = (((((((vars[CARBON_POOL])+(((vars[CARBON_POOL])*((P_phot_c) - (R_C))*(stepInHours))))) >= (param_C_rep))&&((((vars[SILICATE_POOL])+(vars[SILICATE_ING]))) >= (param_S_rep))))?(2.0):(1.0));

    /* Update chemical pools */
    double Q_nitrate = ((((vars[NITRATE_POOL])+(vars[NITRATE_ING]))) / (vars[CARBON_POOL]));
    double Q_ammonium = ((((vars[AMMONIUM_POOL])+(vars[AMMONIUM_ING]))) / (vars[CARBON_POOL]));
    double omega = ((((param_k_AR) / (((param_k_AR)+(_ambientAmmonium))))*((((param_k_AR)+(_ambientNitrate))) / (((param_k_AR)+(_ambientAmmonium)+(_ambientNitrate))))));
    double V_max_C = (((((vars[AMMONIUM_POOL])+(vars[NITRATE_POOL]))) < (1000.0))?(((((Q_ammonium)+(Q_nitrate))) < (param_Q_Nmin))?(((param_V_ref_c)*(T_function))):(((((Q_ammonium)+(Q_nitrate))) > (param_Q_Nmax))?(0.0):(((param_V_ref_c)*(pow(((param_Q_Nmax) - (((Q_ammonium)+(Q_nitrate)))) / ((param_Q_Nmax) - (param_Q_Nmin)), 0.05))*(T_function))))):(0.0));
    double V_C_ammonium = (((V_max_C)*((_ambientAmmonium) / (((param_k_AR)+(_ambientAmmonium))))));
    double V_C_nitrate = (((V_max_C)*((_ambientNitrate) / (((param_k_AR)+(_ambientNitrate))))*(omega)));
    double V_S_max = (((vars[CARBON_POOL]) >= (param_C_minS))?(((Q_s) <= (param_Q_S_min))?(((param_V_S_ref)*(T_function))):(((Q_s) >= (param_Q_S_max))?(0.0):(((param_V_S_ref)*(pow(((param_Q_S_max) - (Q_s)) / ((param_Q_S_max) - (param_Q_S_min)), 0.05))*(T_function))))):(0.0));
    double V_S_S = (((V_S_max)*((_ambientSilicate) / (((_ambientSilicate)+(param_k_S))))));

    double amountAmm = ((vars[CARBON_POOL])*(V_C_ammonium)*(stepInHours));
    double amountNit = ((vars[CARBON_POOL])*(V_C_nitrate)*(stepInHours));
    double amountSi = ((vars[SILICATE_POOL])*(V_S_S)*(stepInHours));

    double _Ammonium_Pool_New = (((((vars[AMMONIUM_POOL])+(vars[AMMONIUM_ING])+
        (vars[NITRATE_ING]))) - (((vars[AMMONIUM_POOL])*(param_R_N)*(stepInHours)*(T_function)))) /
        (C_d));
    double _Nitrate_Pool_New = (0.0);
    double _Silicate_Pool_New = ((((vars[SILICATE_POOL])+(vars[SILICATE_ING]))) / (C_d));

    /* Determine death... */
    double death_flag = (((vars[CARBON_POOL]) <= (param_C_starve))?(1.0):(0.0));
    /* O if still alive, 1.0 if dead */
    vars[STAGE] = death_flag;

    /* ...before updating carbon pool */
    double _Carbon_Pool_New = (((death_flag) == (0.0))?(((((((vars[CARBON_POOL])*
               ((P_phot_c)-(((R_C)*(T_function))))*(stepInHours)))+(vars[CARBON_POOL]))) - (0.0)) / (C_d)):(0.0));

    /* Motion */
/*
    double random = rnd(physics->turbocline);
    vars[Z_NEW]=(((vars[Z_OLD]) <= (physics->turbocline))?(((random)+(((param_z_sink)*(stepInHours))))):(((vars[Z_OLD])+(((param_z_sink)*(stepInHours))))));
    int z_new_int = (int)(vars[Z_NEW] + 0.0001f);
*/

    /* External: ammonium excretion */
    double _relAmount = ((((vars[AMMONIUM_POOL])+(vars[NITRATE_POOL])))*(param_R_N)*(stepInHours)*(T_function));


    /* Phase 3
     * (output and housekeeping)
     */

    /* Ensemble update (reproduction) */
    if ((C_d) == (2.0)) {
      vars[C_NEW]=vars[C_NEW]*2;
    }

    /* Stage change (death) */
    double _Chlorophyll_Pool_New;
    if ((death_flag) == (1.0)) {
/*
		deadPhyto[z_new_int] += vars[C_NEW];
		deadPhyto_req[z_new_int][AMMONIUM] += _Ammonium_Pool_New;
		deadPhyto_req[z_new_int][NITRATE] += _Nitrate_Pool_New;
		deadPhyto_req[z_new_int][SILICATE] += _Silicate_Pool_New;
*/

		/* dead ones have no chlorophyll */
		_Chlorophyll_Pool_New = 0.0;
    } else {
    	/* if alive update chlorophyll */
        if (Theta_N <= param_Theta_max_N) {
            _Chlorophyll_Pool_New = vars[CHLOROPHYLL_POOL] + Rho_Chl * (vars[AMMONIUM_ING] + vars[NITRATE_ING]);
        }
        else {
            _Chlorophyll_Pool_New = param_Theta_max_N * (vars[AMMONIUM_POOL] + vars[NITRATE_POOL]);
        }
        _Chlorophyll_Pool_New -= vars[CHLOROPHYLL_POOL] * param_R_Chl * stepInHours * T_function;
        if (_Chlorophyll_Pool_New > 0.0) {
            _Chlorophyll_Pool_New /= C_d;
        }
        else {
            _Chlorophyll_Pool_New = 0.0;
        }
    }

    //update internal
    vars[AMMONIUM_ING] = amountAmm;
    vars[NITRATE_ING] = amountNit;
    vars[SILICATE_ING] = amountSi;
    vars[AMMONIUM_POOL] = _Ammonium_Pool_New;
    vars[NITRATE_POOL] = _Nitrate_Pool_New;
    vars[CARBON_POOL] = _Carbon_Pool_New;
    vars[SILICATE_POOL] = _Silicate_Pool_New;
    vars[CHLOROPHYLL_POOL] = _Chlorophyll_Pool_New;
    vars[AMMONIUM_REL] = _relAmount;
    vars[SILICATE_REL] = 0.0;

    //update external
/*
    chem_req[z_new_int][SILICATE] += amountSi;
   	chem_req[z_new_int][NITRATE] += amountNit;
    chem_req[z_new_int][AMMONIUM] += amountAmm;
    chem_remin[z_new_int][AMMONIUM] +=_relAmount;
*/
}

void updateDeadDiatom_c(double vars[], int n_vars, double env[], int n_env, double *dt) {
/*
    vars[C_OLD] = vars[C_NEW];
    vars[Z_OLD] = vars[Z_NEW];
    vars[AMMONIUM_ING] /= vars[C_OLD];
    vars[NITRATE_ING] /= vars[C_OLD];
    vars[SILICATE_ING] /= vars[C_OLD];
*/

/*
    if (constData->chem_depletion[z_old_int][AMMONIUM]<1)
        vars[AMMONIUM_ING]*=constData->chem_depletion[z_old_int][AMMONIUM];
    if (constData->chem_depletion[z_old_int][NITRATE]<1)
        vars[NITRATE_ING]*=constData->chem_depletion[z_old_int][NITRATE];
    if (constData->chem_depletion[z_old_int][SILICATE]<1)
        vars[SILICATE_ING]*=constData->chem_depletion[z_old_int][SILICATE];
*/

    //vars[C_NEW] += constData->phyto_dead_depletion[z_old_int];
    /* This needs to be weight averaged... */
    /*    vars[AMMONIUM] += constData->phyto_dead_chem[z_old_int][AMMONIUM];
    vars[SILICATE] += constData->phyto_dead_chem[z_old_int][SILICATE];
    vars[NITRATE] += constData->phyto_dead_chem[z_old_int][NITRATE];*/

    /*****/
    double stepInHours = *dt / 3600.;
    double _ambientTemperature = env[TEMPERATURE];

    //vars[Z_NEW]=(((vars[Z_OLD]) <= (physics->turbocline))?(((rnd(physics->turbocline))+(((param_z_sink)*(stepInHours))))):(((vars[Z_OLD])+(((param_z_sink)*(stepInHours))))));
    double Si_reminT = (((param_S_dis)*(pow(param_Q_remS, ((((_ambientTemperature)+(273.0))) - (param_T_refS)) / (10.0)))));
    double N_reminT = (((param_Ndis)*(pow(param_Q_remN, ((((_ambientTemperature)+(273.0))) - (param_T_refN)) / (10.0)))));
    double _relAmountSi = ((vars[SILICATE_POOL])*(Si_reminT)*(stepInHours));

//    chem_remin[z_old_int][SILICATE]+=_relAmountSi;

    double _Silicate_Pool_New = ((MAX(((((vars[SILICATE_POOL])+(vars[SILICATE_ING]))) - (((vars[SILICATE_POOL])*(Si_reminT)*(stepInHours)))), (0.0))));
    double _relAmountAmm = ((((vars[AMMONIUM_POOL])+(vars[NITRATE_POOL])))*(N_reminT)*(stepInHours));
//    chem_remin[z_old_int][AMMONIUM]+=_relAmountAmm;

    double _Ammonium_Pool_New = ((MAX(((((vars[AMMONIUM_POOL])+(vars[AMMONIUM_ING]))) - (((vars[AMMONIUM_POOL])*(N_reminT)*(stepInHours)))), (0.0))));
    double _Nitrate_Pool_New = ((MAX(((((vars[NITRATE_POOL])+(vars[NITRATE_ING]))) - (((vars[NITRATE_POOL])*(N_reminT)*(stepInHours)))), (0.0))));

    /*****/

    vars[CARBON_POOL]=0;
    vars[CHLOROPHYLL_POOL]=0;
    vars[AMMONIUM_ING]=0;
    vars[NITRATE_ING]=0;
    vars[SILICATE_ING]=0;
    vars[AMMONIUM_POOL] = _Ammonium_Pool_New;
    vars[NITRATE_POOL] = _Nitrate_Pool_New;
    vars[SILICATE_POOL] = _Silicate_Pool_New;
    vars[AMMONIUM_REL] = _relAmountAmm;
    vars[SILICATE_REL] = _relAmountSi;

}

