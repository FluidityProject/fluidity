#include "math.h"

//to call floating-point check from the Fortran routine
void fpcheck_c_( double *a, int *flag ) 
{ 
   *flag = isnan(*a) && !isinf(*a); 
}

