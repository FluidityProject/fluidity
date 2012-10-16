#include "multigridc.h"
#include "petsc.h"

PetscErrorCode MGInfo(void *vobj, const char message[]) {
  PetscErrorCode ierr=PetscInfo(vobj, message);
  return ierr;
}
