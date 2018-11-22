#ifndef HAVE_VTK
#include <stdio.h>
typedef struct _Field_Info Field_Info;

struct _Field_Info
{
    int          ncomponents;
    int          identifier;
    float        interperr; // >0 means adapt+interpolate this field
                            // =0 means interpolate only
                            // <0 means ignore this field
    float        cutoff; // >0 means use rel. Hessian, with abs. cut-off
                         // =0 means use absolute Hessian
                         // <0 means use rel. Hessian, with rel. cut-off
    Field_Info   *next;
    char         *name; // a pointer to a seperately allocated string buffer
};

int readVTKFile(char const* filename, int* npoints, int* ncells, int* ncomps, int* nprop_comps, int* szenls, int* ndim, struct _Field_Info* fieldlst, double** X, double** Y, double** Z, int** ENLBAS, int** ENLIST, double** f, double** p, int onlyinfo, int onlytets) {
  
  return 0;}; 

#endif
