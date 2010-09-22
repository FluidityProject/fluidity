#include "confdefs.h"
#include "petsc.h"
#if PETSC_VERSION_MINOr==0
#include "petscfix.h"
#include "petscis.h"
#endif

#include <string.h>

extern "C" {
  void iscopyindices_(IS *is, PetscInt *iarray,PetscErrorCode *ierr);
  void PETSC_STDCALL petscobjectreference_(PetscObject obj, int *__ierr );
}

// our own version of IsGetIndices, that just does a copy
// to avoid silly fortran interface issues
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

void PETSC_STDCALL petscobjectreference_(PetscObject obj, int *__ierr ){
*__ierr = PetscObjectReference(
   (PetscObject)PetscToPointer((obj) ));
}