#include "confdefs.h"
#include "petsc.h"

#include <string.h>

extern "C" {
  void iscopyindices_(IS *is, PetscInt *iarray,PetscErrorCode *ierr);
  void petscobjectreference_(PetscObject obj, int *__ierr );
  void pcmgsetlevels_nocomms_(PC *pc, PetscInt *levels, int *ierr);
}

// our own version of IsGetIndices, that just does a copy
// to avoid silly fortran interface issues
void iscopyindices_(IS *is, PetscInt *iarray,PetscErrorCode *ierr){
  const PetscInt *iptr;
  PetscInt n;

  *ierr = ISGetIndices(*is, &iptr); if (*ierr) return;
  *ierr = ISGetSize(*is, &n); if (*ierr) return;

  memcpy(iarray, iptr, sizeof(PetscInt)*n);

  *ierr = ISRestoreIndices(*is, &iptr);
}

// PetscObjectReference()
// missing from the fortran interface (copy of petscobjectderefence, which
// does exists, with obvious changes)

#ifdef PETSC_USE_POINTER_CONVERSION
#if defined(__cplusplus)
extern "C" {
#endif
extern void *PetscToPointer(void*);
extern int PetscFromPointer(void *);
extern void PetscRmPointer(void*);
#if defined(__cplusplus)
}
#endif

#else

#define PetscToPointer(a) (*(long *)(a))
#define PetscFromPointer(a) (long)(a)
#define PetscRmPointer(a)
#endif

void petscobjectreference_(PetscObject obj, int *__ierr ){
*__ierr = PetscObjectReference(
   (PetscObject)PetscToPointer((obj) ));
}


// This is a woraround bug in https://gitlab.com/petsc/petsc/-/issues/868
// In 3.14 not only is it impossible to "fake" a NULL argument using PETSC_NULL_KSP
// as we'd done previously, even if you do provide a comms filled with MPI_COMM_WORLD
// that because MPI_COMM_WORLD==0 it is still mistaken for a NULL argument and hits the
// bug in the petsc fortran wrapper for this routine
// So instead we just write our own wrapper without the comms argument:
void pcmgsetlevels_nocomms_(PC *pc, PetscInt *levels, int *ierr){
  *ierr = PCMGSetLevels(*pc, *levels, NULL);
}
