! This header file should be included in fortran modules that use petsc
! directly after the implicit none. It provides legacy support for building
! Fluidity with petsc versions older than the latest released petsc. Where names
! have changed this #defines the new name as its older equivalent, so that new
! names can be used in the code everywhere. Where interfaces have changed we
! still need #ifdef PETSC_VERSION>... in the main code
#include "petscversion.h"
#ifdef HAVE_PETSC_MODULES
#include "finclude/petscdef.h"
#else
#include "finclude/petsc.h"
#endif
! this is the one exception where we keep the old names for now (until
! we get rid of petsc 3.1 and 3.2 support). MatCreate{Seq|MPI}[B]AIJ()
! no longer exists in petsc>=3.3. Instead we need to call MatCreateAIJ()
! followed by MatSetUp() (now required). We provide wrapper routines
! in the petsc_tools module.
#if PETSC_VERSION_MINOR>=3
#define MatCreateSeqAIJ myMatCreateSeqAIJ
#define MatCreateMPIAIJ myMatCreateMPIAIJ
#define MatCreateSeqBAIJ myMatCreateSeqBAIJ
#define MatCreateMPIBAIJ myMatCreateMPIBAIJ
#endif

#if PETSC_VERSION_MINOR<2
#define KSP_NORM_NONE KSP_NORM_NO
#endif
#if PETSC_VERSION_MINOR<3
#define KSPCHEBYSHEV KSPCHEBYCHEV
#define KSPChebyshevSetEigenvalues KSPChebychevSetEigenvalues
#endif
#if PETSC_VERSION_MINOR<2
#define PetscBool PetscTruth
#endif
#if PETSC_VERSION_MINOR<2
#define VecSqrtAbs VecSqrt
#endif
! workaround sily bug in petsc 3.1
#ifndef PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE
#define PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE PC_COMPOSITE_SYM_MULTIPLICATIVE
#endif
! Changes in petsc-dev (master)
! should be changed to use PETSC_DEFAULT_REAL in code when petsc 3.5 is released
! (can't use #ifndef cause PETSC_DEFAULT_REAL is a module variable in petsc-dev)
#if PETSC_VERSION_MINOR>=4 && PETSC_VERSION_RELEASE==0
#define PETSC_DEFAULT_DOUBLE_PRECISION PETSC_DEFAULT_REAL
#endif
! flag argument to KSP/PCSetOperators() has been dropped:
! we use this macro hack which means that the call cannot be split over multiple lines
! also note the (ab)use of fortran's case insensivity to avoid recursion
#if PETSC_VERSION_MINOR>=4 && PETSC_VERSION_RELEASE==0
#define KSPSetOperators(ksp, amat, pmat, flag, ierr) kspsetoperators(ksp, amat, pmat, ierr)
#define PCSetOperators(pc, amat, pmat, flag, ierr) pcsetoperators(pc, amat, pmat, ierr)
#endif
