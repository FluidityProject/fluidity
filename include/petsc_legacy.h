! This header file should be included in fortran modules that use petsc
! directly after the implicit none. It provides legacy support for building
! Fluidity with petsc versions older than the latest released petsc. Where names
! have changed this #defines the new name as its older equivalent, so that new
! names can be used in the code everywhere. Where interfaces have changed we
! still need #ifdef PETSC_VERSION>... in the main code
#include "petscversion.h"
#include "petsc/finclude/petsc.h"

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<23)
#define PetscObjectReferenceWrapper(x, ierr) PetscObjectReference(x%v, ierr)
#else
#define PetscObjectReferenceWrapper(x, ierr) PetscObjectReference(x, ierr)
#endif

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<23)
! from 3.23 this is a parameter, so the following won't work (but is also not required)
#ifndef PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE
#define PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE PC_COMPOSITE_SYM_MULTIPLICATIVE
#endif
#endif

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<22)
#define PETSC_NULL_INTEGER_ARRAY PETSC_NULL_INTEGER
#endif

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<19)
#define PETSC_NULLPTR PETSC_NULL
#endif
