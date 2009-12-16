#ifdef HAVE_PETSC
#include "confdefs.h"
#include "petsc.h"
#include "petscis.h"

#include <string.h>

void iscopyindices_(IS *is, PetscInt *iarray,PetscErrorCode *ierr){
#if PETSC_VERSION_MAJOR==3
  const PetscInt *iptr;
#else
  PetscInt *iptr;
#endif
  PetscInt n;
  
  *ierr = ISGetIndices(*is, &iptr); if (*ierr) return;
  *ierr = ISGetSize(*is, &n); if (*ierr) return;
  
  memcpy(iarray, iptr, sizeof(PetscInt)*n);
  
  *ierr = ISRestoreIndices(*is, &iptr);
}
#endif
