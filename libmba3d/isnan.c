#include "confdefs.h"
//to call floating-point check from the Fortran routine
#define fpcheck_mba3d F77_FUNC(fpcheck_mba3d, FPCHECK_MBA3D)
void fpcheck_mba3d( double *a, int *flag ) 
{ 
   *flag = isnan(*a) && !isinf(*a); 
}

