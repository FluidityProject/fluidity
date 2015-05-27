! This header file should be included in fortran modules that use petsc
! directly after the implicit none. It provides legacy support for building
! Fluidity with petsc versions older than the latest released petsc. Where names
! have changed this #defines the new name as its older equivalent, so that new
! names can be used in the code everywhere. Where interfaces have changed we
! still need #ifdef PETSC_VERSION>... in the main code
#include "petscversion.h"
#ifdef HAVE_PETSC_MODULES
#if PETSC_VERSION_RELEASE==0
#include "petsc/finclude/petscdef.h"
#else
#include "finclude/petscdef.h"
#endif
#else
#if PETSC_VERSION_RELEASE==0
#include "petsc/finclude/petsc.h"
#else
#include "finclude/petsc.h"
#endif
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
! renamed in petsc master
#if PETSC_VERSION_RELEASE==0
#define MatGetVecs MatCreateVecs
#endif
