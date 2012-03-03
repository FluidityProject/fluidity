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
#define SILICATE_REL         8
#define SILICATE_POOL        9          /* Silicate pool */
#define AMMONIUM_ING         10          /* Ammonium Uptake */
#define AMMONIUM_REL         11
#define AMMONIUM_POOL        12         /* Current ammonium pool */

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

void updateLivingDiatom_c(double vars[], int n_vars, double env[], int n_env, double *dt);
void updateDeadDiatom_c(double vars[], int n_vars, double env[], int n_env, double *dt);
