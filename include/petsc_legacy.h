! This header file should be included in fortran modules that use petsc
! directly after the implicit none. It provides legacy support for building
! Fluidity with petsc versions older than the latest released petsc. Where names
! have changed this #defines the new name as its older equivalent, so that new
! names can be used in the code everywhere. Where interfaces have changed we
! still need #ifdef PETSC_VERSION>... in the main code
#include "petscversion.h"
#include "petsc/finclude/petsc.h"

#define PetscObjectReferenceWrapper(x, ierr) PetscObjectReference(x%v, ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<9)
#define MatSolverType MatSolverPackage
#define PCFactorSetMatSolverType PCFactorSetMatSolverPackage
#endif

#ifndef PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE
#define PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE PC_COMPOSITE_SYM_MULTIPLICATIVE
#endif

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<15)
#define PCCompositeAddPCType PCCompositeAddPC
#define KSPMonitorResidual KSPMonitorDefault
#define KSPMonitorTrueResidual KSPMonitorTrueResidualNorm
#endif
