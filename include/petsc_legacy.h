! This header file should be included in fortran modules that use petsc
! directly after the implicit none. It provides legacy support for building
! Fluidity with petsc versions older than the latest released petsc. Where names
! have changed this #defines the new name as its older equivalent, so that new
! names can be used in the code everywhere. Where interfaces have changed we
! still need #ifdef PETSC_VERSION>... in the main code
#include "petscversion.h"
#ifdef HAVE_PETSC_MODULES
#if PETSC_VERSION_MINOR>=6
#include "petsc/finclude/petscdef.h"
#else
#include "finclude/petscdef.h"
#endif
#else
#if PETSC_VERSION_MINOR>=6
#include "petsc/finclude/petsc.h"
#else
#include "finclude/petsc.h"
#endif
#endif

#ifndef PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE
#define PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE PC_COMPOSITE_SYM_MULTIPLICATIVE
#endif
! Changes in petsc 3.5 PETSC_DEFAULT_DOUBLE_PRECISION -> PETSC_DEFAULT_REAL
! (can't use #ifndef cause PETSC_DEFAULT_REAL is a module variable)
#if PETSC_VERSION_MINOR<5
#define PETSC_DEFAULT_REAL PETSC_DEFAULT_DOUBLE_PRECISION
#endif
! MatStructure argument to KSP/PCSetOperators has been dropped:
! we use this macro hack which means that the call cannot be split over multiple lines
! also note the (ab)use of fortran's case insensivity to avoid recursion
! for PCSetOperators we use small caps because otherwise the preprocessor trips up on an occurence of PCSetoperators directly followed by () in a fortran comment in one of the petsc headers
#if PETSC_VERSION_MINOR<5
#define KSPSetOperators(ksp, amat, pmat, ierr) kspsetoperators(ksp, amat, pmat, DIFFERENT_NONZERO_PATTERN, ierr)
#define pcsetoperators(pc, amat, pmat, ierr) PCSetOperators(pc, amat, pmat, DIFFERENT_NONZERO_PATTERN, ierr)
! mykspgetoperators is a wrapper function defined in Petsc_Tools.F90:
#define KSPGetOperators(ksp, amat, pmat, ierr) mykspgetoperators(ksp, amat, pmat, ierr)
#endif
#if PETSC_VERSION_MINOR<6
#define MatCreateVecs MatGetVecs
#endif
! from petsc 3.7 onward all PetscOptionsXXX() calls have an additional PetscOptions first argument
#if PETSC_VERSION_MINOR<7
#define PetscOptionsGetInt(options, pre, name, ivalue, set, ierr) petscoptionsgetint(pre, name, ivalue, set, ierr)
#define PetscOptionsGetReal(options, pre, name, dvalue, set, ierr) petscoptionsgetreal(pre, name, dvalue, set, ierr)
#define PetscOptionsGetString(options, pre, name, string, set, ierr) petscoptionsgetstring(pre, name, string, set, ierr)
#define PetscOptionsHasName(options, pre, name, set, ierr) petscoptionshasname(pre, name, set, ierr)
#endif
! these PetscViewerAndFormatCreate/Destroy don't exist in petsc<3.7
! but we only use them to combine PETSC_VIEWER_STDOUT_WORLD and PETSC_VIEWER_DEFAULT
#if PETSC_VERSION_MINOR<7
#define PetscViewerAndFormatCreate NullPetscViewerAndFormatCreate
#define PetscViewerAndFormatDestroy PETSC_NULL_FUNCTION
#endif
