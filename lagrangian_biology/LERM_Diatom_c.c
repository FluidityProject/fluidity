#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Ambient Indexes */
#define AMMONIUM             0
#define NITRATE              1
#define SILICATE             2
#define TEMPERATURE          3
#define IRRADIANCE           4

#define STAGE                0         /* Stage (Living=0, Dead=1) */
#define C_NEW                1          /* Sub-population size next timestep */
#define C_OLD                2          /* Sub-population size last timestep */
#define CARBON_POOL          3          /* Carbon pool */
#define CHLOROPHYLL_POOL     4          /* Chlorophyll pool */
#define NITRATE_ING          5          /* nitrate uptake */
#define NITRATE_POOL         6          /* Nitrate pool */
#define SILICATE_ING         7          /* Silicate uptake */
#define SILICATE_POOL        8          /* Silicate pool */
#define AMMONIUM_ING         9          /* Ammonium Uptake */
#define AMMONIUM_POOL        10         /* Current ammonium pool */

#define _STAGE_Living        0
#define _STAGE_Dead          1

/* Macros */
#define MAX(A,B)             (A)>=(B)?(A):(B)
#define MIN(A,B)             (A)<(B)?(A):(B)

/* Model Constants */
#define param_A_E            -10000
#define param_Alpha_Chl      7.9E-7
#define param_C_minS         1.584E-8
#define param_C_rep          1.76E-8
#define param_C_starve       8.5E-9
#define param_k_AR           1.0
#define param_k_S            1.0
#define param_Ndis           0.0042
#define param_P_ref_c        0.14
#define param_Q_Nmax         0.17
#define param_Q_Nmin         0.034
#define param_Q_remN         2.95
#define param_Q_remS         2.27
#define param_Q_S_max        0.15
#define param_Q_S_min        0.04
#define param_R_Chl          2E-3
#define param_R_maintenance  2E-3
#define param_R_N            2E-3
#define param_S_dis          8.3E-4
#define param_S_rep          2.1E-9
#define param_T_ref          293.0
#define param_T_refN         283.0
#define param_T_refS         278.0
#define param_Theta_max_N    4.2
#define param_V_ref_c        0.01
#define param_V_S_ref        0.03
#define param_z_sink         0.04
#define param_Zeta           2.3




/*****************************/
/* MAIN PLANKTON UPDATE CODE */
/*****************************/

static void updateLivingDiatom_c(float *vars, float *env, float dt) { //int agent) {

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

    float stepInHours = dt / 3600.;

    float _ambientTemperature = env[TEMPERATURE];
    float _ambientVisIrrad = env[IRRADIANCE];
    float _ambientAmmonium = env[AMMONIUM];
    float _ambientNitrate = env[NITRATE];
    float _ambientSilicate = env[SILICATE];

    float T_K = _ambientTemperature + 273.0;
    float T_function = exp(34.12969283 - 10000/T_K);

    /* Photosynthesis */
    float Q_N = ((((vars[AMMONIUM_POOL])+(vars[AMMONIUM_ING])+(vars[NITRATE_POOL])+(vars[NITRATE_ING]))) / (vars[CARBON_POOL]));
    float Q_s = ((((vars[SILICATE_POOL])+(vars[SILICATE_ING]))) / (vars[CARBON_POOL]));
    float P_max_c = (((Q_N) > (param_Q_Nmax))?(((param_P_ref_c)*(T_function))):(((Q_N) < (param_Q_Nmin))?(0.0):(((param_P_ref_c)*(T_function)*(((Q_N) - (param_Q_Nmin)) / ((param_Q_Nmax) - (param_Q_Nmin)))))));
    float Theta_c = ((vars[CHLOROPHYLL_POOL]) / (vars[CARBON_POOL]));
    float E_0 = (((4.6)*(_ambientVisIrrad)));

    float P_phot_c = ((P_max_c == 0.0)||(Q_s <= param_Q_S_min))?0.0:P_max_c*(1.0 - exp(-0.284400e-2*Theta_c*E_0 / P_max_c));

    /* Reproduction */
    float Theta_N = ((vars[CHLOROPHYLL_POOL]) / (((vars[AMMONIUM_POOL])+(vars[AMMONIUM_ING])+(vars[NITRATE_POOL])+(vars[NITRATE_ING]))));
    float Rho_Chl = (((((E_0) > (0.0))&&((Theta_c) > (0.0))))?(((param_Theta_max_N)*((P_phot_c) / (((3600.0)*(param_Alpha_Chl)*(Theta_c)*(E_0)))))):(0.0));
    float R_C_growth = ((((((vars[AMMONIUM_ING])+(vars[NITRATE_ING])))*(param_Zeta))) / (((stepInHours)*(vars[CARBON_POOL]))));
    float R_C = (((param_R_maintenance)+(R_C_growth)));
    float C_d = (((((((vars[CARBON_POOL])+(((vars[CARBON_POOL])*((P_phot_c) - (R_C))*(stepInHours))))) >= (param_C_rep))&&((((vars[SILICATE_POOL])+(vars[SILICATE_ING]))) >= (param_S_rep))))?(2.0):(1.0));

    /* Update chemical pools */
    float Q_nitrate = ((((vars[NITRATE_POOL])+(vars[NITRATE_ING]))) / (vars[CARBON_POOL]));
    float Q_ammonium = ((((vars[AMMONIUM_POOL])+(vars[AMMONIUM_ING]))) / (vars[CARBON_POOL]));
    float omega = ((((param_k_AR) / (((param_k_AR)+(_ambientAmmonium))))*((((param_k_AR)+(_ambientNitrate))) / (((param_k_AR)+(_ambientAmmonium)+(_ambientNitrate))))));
    float V_max_C = (((((vars[AMMONIUM_POOL])+(vars[NITRATE_POOL]))) < (1000.0))?(((((Q_ammonium)+(Q_nitrate))) < (param_Q_Nmin))?(((param_V_ref_c)*(T_function))):(((((Q_ammonium)+(Q_nitrate))) > (param_Q_Nmax))?(0.0):(((param_V_ref_c)*(pow(((param_Q_Nmax) - (((Q_ammonium)+(Q_nitrate)))) / ((param_Q_Nmax) - (param_Q_Nmin)), 0.05))*(T_function))))):(0.0));
    float V_C_ammonium = (((V_max_C)*((_ambientAmmonium) / (((param_k_AR)+(_ambientAmmonium))))));
    float V_C_nitrate = (((V_max_C)*((_ambientNitrate) / (((param_k_AR)+(_ambientNitrate))))*(omega)));
    float V_S_max = (((vars[CARBON_POOL]) >= (param_C_minS))?(((Q_s) <= (param_Q_S_min))?(((param_V_S_ref)*(T_function))):(((Q_s) >= (param_Q_S_max))?(0.0):(((param_V_S_ref)*(pow(((param_Q_S_max) - (Q_s)) / ((param_Q_S_max) - (param_Q_S_min)), 0.05))*(T_function))))):(0.0));
    float V_S_S = (((V_S_max)*((_ambientSilicate) / (((_ambientSilicate)+(param_k_S))))));

    float amountAmm = vars[C_OLD]*((vars[CARBON_POOL])*(V_C_ammonium)*(stepInHours));
    float amountNit = vars[C_OLD]*((vars[CARBON_POOL])*(V_C_nitrate)*(stepInHours));
    float amountSi = vars[C_OLD]*((vars[SILICATE_POOL])*(V_S_S)*(stepInHours));

    float _Ammonium_Pool_New = (((((vars[AMMONIUM_POOL])+(vars[AMMONIUM_ING])+
        (vars[NITRATE_ING]))) - (((vars[AMMONIUM_POOL])*(param_R_N)*(stepInHours)*(T_function)))) /
        (C_d));
    float _Nitrate_Pool_New = (0.0);
    float _Silicate_Pool_New = ((((vars[SILICATE_POOL])+(vars[SILICATE_ING]))) / (C_d));

    /* Determine death... */
    float death_flag = (((vars[CARBON_POOL]) <= (param_C_starve))?(1.0):(0.0));
    /* O if still alive, 1.0 if dead */
    vars[STAGE] = death_flag;

    /* ...before updating carbon pool */
    float _Carbon_Pool_New = (((death_flag) == (0.0))?(((((((vars[CARBON_POOL])*
               ((P_phot_c)-(((R_C)*(T_function))))*(stepInHours)))+(vars[CARBON_POOL]))) - (0.0)) / (C_d)):(0.0));

    /* Motion */
/*
    float random = rnd(physics->turbocline);
    vars[Z_NEW]=(((vars[Z_OLD]) <= (physics->turbocline))?(((random)+(((param_z_sink)*(stepInHours))))):(((vars[Z_OLD])+(((param_z_sink)*(stepInHours))))));
    int z_new_int = (int)(vars[Z_NEW] + 0.0001f);
*/

    /* External: ammonium excretion */
    float _relAmount = vars[C_OLD]*((((vars[AMMONIUM_POOL])+(vars[NITRATE_POOL])))*(param_R_N)*(stepInHours)*(T_function));


    /* Phase 3
     * (output and housekeeping)
     */

    /* Ensemble update (reproduction) */
    if ((C_d) == (2.0)) {
      vars[C_NEW]=vars[C_OLD]*2;
    }

    /* Stage change (death) */
    float _Chlorophyll_Pool_New;
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

    //update external
/*
    chem_req[z_new_int][SILICATE] += amountSi;
   	chem_req[z_new_int][NITRATE] += amountNit;
    chem_req[z_new_int][AMMONIUM] += amountAmm;
    chem_remin[z_new_int][AMMONIUM] +=_relAmount;
*/
}

static void updateDeadDiatom(float *vars, float *env, float dt) {
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
    float stepInHours = dt / 3600.;
    float _ambientTemperature = env[TEMPERATURE];

    //vars[Z_NEW]=(((vars[Z_OLD]) <= (physics->turbocline))?(((rnd(physics->turbocline))+(((param_z_sink)*(stepInHours))))):(((vars[Z_OLD])+(((param_z_sink)*(stepInHours))))));
    float Si_reminT = (((param_S_dis)*(pow(param_Q_remS, ((((_ambientTemperature)+(273.0))) - (param_T_refS)) / (10.0)))));
    float N_reminT = (((param_Ndis)*(pow(param_Q_remN, ((((_ambientTemperature)+(273.0))) - (param_T_refN)) / (10.0)))));
    float _relAmountSi = vars[C_OLD]*((vars[SILICATE_POOL])*(Si_reminT)*(stepInHours));

//    chem_remin[z_old_int][SILICATE]+=_relAmountSi;

    float _Silicate_Pool_New = ((MAX(((((vars[SILICATE_POOL])+(vars[SILICATE_ING]))) - (((vars[SILICATE_POOL])*(Si_reminT)*(stepInHours)))), (0.0))));
    float _relAmountAmm = vars[C_OLD]*((((vars[AMMONIUM_POOL])+(vars[NITRATE_POOL])))*(N_reminT)*(stepInHours));
//    chem_remin[z_old_int][AMMONIUM]+=_relAmountAmm;

    float _Ammonium_Pool_New = ((MAX(((((vars[AMMONIUM_POOL])+(vars[AMMONIUM_ING]))) - (((vars[AMMONIUM_POOL])*(N_reminT)*(stepInHours)))), (0.0))));
    float _Nitrate_Pool_New = ((MAX(((((vars[NITRATE_POOL])+(vars[NITRATE_ING]))) - (((vars[NITRATE_POOL])*(N_reminT)*(stepInHours)))), (0.0))));

    /*****/

    vars[CARBON_POOL]=0;
    vars[CHLOROPHYLL_POOL]=0;
    vars[AMMONIUM_ING]=0;
    vars[NITRATE_ING]=0;
    vars[SILICATE_ING]=0;
    vars[AMMONIUM_POOL] = _Ammonium_Pool_New;
    vars[NITRATE_POOL] = _Nitrate_Pool_New;
    vars[SILICATE_POOL] = _Silicate_Pool_New;

}

