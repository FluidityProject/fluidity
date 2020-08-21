dnl @synopsis ACX_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the BLAS
dnl linear-algebra interface (see http://www.netlib.org/blas/).
dnl On success, it sets the BLAS_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with BLAS, you should link with:
dnl
dnl 	$BLAS_LIBS $LIBS $FCLIBS
dnl
dnl in that order.  FCLIBS is the output variable of the
dnl AC_FC_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with FC libraries.
dnl Users will also need to use AC_FC_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl Many libraries are searched for, from ATLAS to CXML to ESSL.
dnl The user may also use --with-blas=<lib> in order to use some
dnl specific BLAS library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the FC env. var.) as
dnl was used to compile the BLAS library.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a BLAS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_BLAS.
dnl
dnl This macro requires autoconf 2.50 or later.
dnl
dnl @version $Id$
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
dnl Modified by Jonas Juselius <jonas@iki.fi>
dnl
AC_DEFUN([ACX_BLAS], [
AC_PREREQ(2.59)

acx_blas_ok=no
acx_blas_save_LIBS="$LIBS"
acx_blas_save_LDFLAGS="$LDFLAGS"
acx_blas_save_FFLAGS="$FFLAGS"
acx_blas_libs=""
acx_blas_dir=""

AC_ARG_WITH(blas,
	[AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])

case $with_blas in
	yes | "") ;;
	no) acx_blas_ok=disable ;;
	-l* | */* | *.a | *.so | *.so.* | *.o) acx_blas_libs="$with_blas" ;;
	*) acx_blas_libs="-l$with_blas" ;;
esac

AC_ARG_WITH(blas_dir,
	[AC_HELP_STRING([--with-blas-dir=<dir>], [look for BLAS library in <dir>])])

case $with_blas_dir in
      yes | no | "") ;;
     -L*) LDFLAGS="$LDFLAGS $with_blas_dir" 
	      acx_blas_dir="$with_blas_dir" ;;
      *) LDFLAGS="$LDFLAGS -L$with_blas_dir" 
	      acx_blas_dir="-L$with_blas_dir" ;;
esac

# Are we linking from C?
case "$ac_ext" in
  f*|F*) sgemm="sgemm" ;;
  *)
   AC_FC_FUNC([sgemm])
   LIBS="$LIBS $FCLIBS"
   ;;
esac

# If --with-blas is defined, then look for THIS AND ONLY THIS blas lib
if test $acx_blas_ok = no; then
case $with_blas in
    ""|yes) ;;
	*) save_LIBS="$LIBS"; LIBS="$acx_blas_libs $LIBS"
	AC_MSG_CHECKING([for $sgemm in $acx_blas_libs])
	AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes])
	AC_MSG_RESULT($acx_blas_ok)
	LIBS="$save_LIBS"
	acx_blas_ok=specific
	;;
esac
fi

# First, check BLAS_LIBS environment variable
if test $acx_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
	AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes; acx_blas_libs=$BLAS_LIBS])
	AC_MSG_RESULT($acx_blas_ok)
	LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
	AC_MSG_CHECKING([for builtin $sgemm])
	AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes])
	AC_MSG_RESULT($acx_blas_ok)
fi

# Intel mkl BLAS. Unfortunately some of Intel's blas routines are
# in their lapack library...
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(mkl_def, $sgemm, 
	[acx_blas_ok=yes; acx_blas_libs="-lmkl_def -lm"],
	[],[-lm])
fi
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(mkl_ipf, $sgemm, 
	[acx_blas_ok=yes; acx_blas_libs="-lmkl_ipf -lguide -lm"],
	[],[-lguide -lm])
fi
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(mkl_em64t, $sgemm,
        [acx_blas_ok=yes; acx_blas_libs="-lmkl_em64t -lguide -liomp5"],
        [],[-lguide -liomp5])
fi
# check for older mkl
if test $acx_blas_ok = no; then
	AC_MSG_NOTICE([trying Intel MKL < 7:])
	unset ac_cv_lib_mkl_def_sgemm
	AC_CHECK_LIB(mkl_lapack, lsame, [
	    acx_lapack_ok=yes;
		AC_CHECK_LIB(mkl_def, $sgemm, 
			[acx_blas_ok=yes; 
			acx_blas_libs="-lmkl_def -lmkl_lapack -lm -lpthread"],
			[],[-lm -lpthread
		])
	])
	AC_MSG_NOTICE([Intel MKL < 7... $acx_blas_ok])
fi

# BLAS in ACML (pgi)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(acml, $sgemm, [acx_blas_ok=yes; acx_blas_libs="-lacml"])
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(f77blas, $sgemm,
		[acx_blas_ok=yes; acx_blas_libs="-lf77blas -latlas"],
		[], [-latlas])
fi

# ia64-hp-hpux11.22 BLAS library?
if test $acx_blas_ok = no; then
        AC_CHECK_LIB(veclib, $sgemm, 
		[acx_blas_ok=yes; acx_blas_libs="-lveclib8"])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
    AC_MSG_NOTICE([trying PhiPACK:])
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(dgemm, dgemm,
			[AC_CHECK_LIB(sgemm, $sgemm,
			[acx_blas_ok=yes; acx_blas_libs="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
		[], [-lblas])
	])
    AC_MSG_NOTICE([PhiPACK... $acx_blas_ok])
fi

# BLAS in Alpha CXML library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(cxml, $sgemm, [acx_blas_ok=yes;acx_blas_libs="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(dxml, $sgemm, [acx_blas_ok=yes;acx_blas_libs="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $sgemm,
        			[acx_blas_libs="-xlic_lib=sunperf -lsunmath"
                    acx_blas_ok=yes],[],[-lsunmath])
		])
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(scs, $sgemm, [acx_blas_ok=yes; acx_blas_libs="-lscs"])
fi

# BLAS in SGIMATH library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $sgemm,
		     [acx_blas_ok=yes; acx_blas_libs="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
    unset ac_cv_lib_blas_sgemm
	AC_MSG_NOTICE([trying IBM ESSL:])
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(essl, $sgemm,
			[acx_blas_ok=yes; acx_blas_libs="-lessl -lblas"],
			[], [-lblas])
	])
	AC_MSG_NOTICE([IBM ESSL... $acx_blas_ok])
fi

# Generic BLAS library?
if test $acx_blas_ok = no; then
    unset ac_cv_lib_blas_sgemm
	AC_CHECK_LIB(blas, $sgemm, [acx_blas_ok=yes; acx_blas_libs="-lblas"])
fi

# blas on SGI/CRAY 
if test $acx_blas_ok = no; then
    unset ac_cv_lib_blas_sgemm
	AC_CHECK_LIB(blas, $sgemm, 
	[acx_blas_ok=yes; acx_blas_libs="-lblas -lcraylibs"],[],[-lcraylibs])
fi

# Check for vecLib framework (Darwin)
if test $acx_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="-framework vecLib $LIBS"
	AC_MSG_CHECKING([for $sgemm in vecLib])
	AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes; acx_blas_libs="-framework vecLib"])
	AC_MSG_RESULT($acx_blas_ok)
	LIBS="$save_LIBS"
fi
 
BLAS_LIBS="$acx_blas_libs"
AC_SUBST(BLAS_LIBS)

LIBS="$acx_blas_save_LIBS"
LDFLAGS="$acx_blas_save_LDFLAGS $acx_blas_dir"

test x"$acx_blas_ok" = xspecific && acx_blas_ok=yes
# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        acx_blas_ok=no
        $2
fi
])dnl ACX_BLAS

dnl @synopsis ACX_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the LAPACK
dnl linear-algebra interface (see http://www.netlib.org/lapack/).
dnl On success, it sets the LAPACK_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with LAPACK, you should link with:
dnl
dnl 	$LAPACK_LIBS $BLAS_LIBS $LIBS 
dnl
dnl in that order.  BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl The user may also use --with-lapack=<lib> in order to use some
dnl specific LAPACK library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the LAPACK and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a LAPACK
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_LAPACK.
dnl
dnl @version $Id$
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
AC_DEFUN([ACX_LAPACK], [
AC_REQUIRE([ACX_BLAS])
acx_lapack_ok=no
acx_lapack_save_LIBS="$LIBS"
acx_lapack_save_LDFLAGS="$LDFLAGS"
acx_lapack_save_FFLAGS="$FFLAGS"
acx_lapack_libs=""
acx_lapack_dir=""

AC_ARG_WITH(lapack,
	[AC_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])

case $with_lapack in
	yes | "") ;;
	no) acx_lapack_ok=disable ;;
	-l* | */* | *.a | *.so | *.so.* | *.o) acx_lapack_libs="$with_lapack" ;;
	*) acx_lapack_libs="-l$with_lapack" ;;
esac

AC_ARG_WITH(lapack_dir,
	[AC_HELP_STRING([--with-lapack-dir=<dir>], [look for LAPACK library in <dir>])])

case $with_lapack_dir in
      yes | no | "") ;;
     -L*) LDFLAGS="$LDFLAGS $with_lapack_dir" 
	      acx_lapack_dir="$with_lapack_dir" ;;
      *) LDFLAGS="$LDFLAGS -L$with_lapack_dir" 
	      acx_lapack_dir="-L$with_lapack_dir" ;;
esac

# We cannot use LAPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
	acx_lapack_ok=noblas
fi

# Are we linking from C?
case "$ac_ext" in
  f*|F*) dsyev="dsyev" ;;
  *)
   AC_FC_FUNC([dsyev])
   LIBS="$LIBS $FCLIBS"
   ;;
esac

# If --with-lapack is defined, then look for THIS AND ONLY THIS lapack lib
if test $acx_lapack_ok = no; then
case $with_lapack in
    ""|yes) ;;
	*) save_LIBS="$LIBS"; LIBS="$acx_lapack_libs $LIBS"
	AC_MSG_CHECKING([for $dsyev in $acx_lapack_libs])
	AC_TRY_LINK_FUNC($dsyev, [acx_lapack_ok=yes])
	AC_MSG_RESULT($acx_lapack_ok)
	LIBS="$save_LIBS"
	acx_lapack_ok=yes
	;;
esac
fi

# First, check LAPACK_LIBS environment variable
if test $acx_lapack_ok = no; then
if test "x$LAPACK_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $LIBS"
	AC_MSG_CHECKING([for $dsyev in $LAPACK_LIBS])
	AC_TRY_LINK_FUNC($dsyev, [acx_lapack_ok=yes; 
	     acx_lapack_libs=$LAPACK_LIBS])
	AC_MSG_RESULT($acx_lapack_ok)
	LIBS="$save_LIBS"
fi
fi

# Intel MKL LAPACK?
if test $acx_lapack_ok = no; then
	AC_CHECK_LIB(mkl_lapack, $dsyev, 
	[acx_lapack_ok=yes; acx_lapack_libs="-lmkl_lapack -lguide"],
	[],[])
fi

# Sun sunperf?
if test $acx_lapack_ok = no; then
	AC_CHECK_LIB(sunperf, $dsyev, 
	[acx_lapack_ok=yes; acx_lapack_libs="-lsunperf"],
	[],[])
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $acx_lapack_ok = no; then
	AC_MSG_CHECKING([for $dsyev in BLAS library])
	AC_TRY_LINK_FUNC($dsyev, [acx_lapack_ok=yes; acx_lapack_libs=""])
	AC_MSG_RESULT($acx_lapack_ok)
fi

# Generic LAPACK library?
if test $acx_lapack_ok = no; then
	AC_CHECK_LIB(lapack, $dsyev,
		[acx_lapack_ok=yes; acx_lapack_libs="-llapack"], [], [])
fi

LAPACK_LIBS="$LAPACK_LIBS $acx_lapack_libs"
LIBS="$acx_lapack_save_LIBS"
LDFLAGS="$acx_lapack_save_LDFLAGS $acx_lapack_dir"

AC_SUBST(LAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
        acx_lapack_ok=no
        $2
fi
])dnl ACX_LAPACK
dnl ----------------------------------------------------------------------------
dnl check for the required PETSc library
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_PETSc], [
AC_REQUIRE([ACX_BLAS])
BLAS_LIBS="$BLAS_LIBS $FLIBS"
AC_REQUIRE([ACX_LAPACK])
LAPACK_LIBS="$LAPACK_LIBS $BLAS_LIBS"
AC_PATH_XTRA

if test "x$PETSC_DIR" == "x"; then
  AC_MSG_WARN( [No PETSC_DIR set - trying to autodetect] )

  # Try to identify the obvious choice
  if test -f /usr/lib/petsc/include/petscversion.h ; then
    export PETSC_DIR=/usr/lib/petsc
  elif test -f /usr/include/petscversion.h ; then
    export PETSC_DIR=/usr
  elif test -f /usr/local/include/petscversion.h ; then
    export PETSC_DIR=/usr/local
  fi
fi
# Check again incase we failed.
if test "x$PETSC_DIR" == "x"; then 
  AC_MSG_WARN( [No PETSC_DIR set - do you need to load a petsc module?] )
  AC_MSG_ERROR( [You need to set PETSC_DIR to point at your PETSc installation... exiting] )
fi
AC_MSG_NOTICE([Using PETSC_DIR=$PETSC_DIR])

PETSC_LINK_LIBS=`make -s -f petsc_makefile_old getlinklibs 2> /dev/null || make -s -f petsc_makefile getlinklibs`
LIBS="$PETSC_LINK_LIBS $LIBS"

# need to add -I$PWD/include/ to what we get from petsc, so we can use our own petsc_legacy.h wrapper
PETSC_INCLUDE_FLAGS=`make -s -f petsc_makefile_old getincludedirs 2> /dev/null || make -s -f petsc_makefile getincludedirs`
CPPFLAGS="$CPPFLAGS $PETSC_INCLUDE_FLAGS -I$PWD/include/"
FCFLAGS="$FCFLAGS $PETSC_INCLUDE_FLAGS -I$PWD/include/"

# first check we have the right petsc version
AC_COMPUTE_INT(PETSC_VERSION_MAJOR, "PETSC_VERSION_MAJOR", [#include "petscversion.h"], 
  [AC_MSG_ERROR([Unknown petsc major version])])
AC_COMPUTE_INT(PETSC_VERSION_MINOR, "PETSC_VERSION_MINOR", [#include "petscversion.h"], 
  [AC_MSG_ERROR([Unknown petsc minor version])])
AC_MSG_NOTICE([Detected PETSc version "$PETSC_VERSION_MAJOR"."$PETSC_VERSION_MINOR"])
# if major<3 or minor<4
if test "0$PETSC_VERSION_MAJOR" -lt 3 -o "0$PETSC_VERSION_MINOR" -lt 4; then
  AC_MSG_ERROR([Fluidity needs PETSc version >=3.4])
fi

AC_LANG(Fortran)
# F90 (capital F) to invoke preprocessing
ac_ext=F90

# may need some extra flags for testing that we don't want to use in the end
SAVE_FCFLAGS="$FCFLAGS"

if test "$enable_petsc_fortran_modules" != "no" ; then
  # now try if the petsc fortran modules work:
  AC_LINK_IFELSE(
          [AC_LANG_PROGRAM([],[[
                          use petsc
                          integer :: ierr
                          print*, "hello petsc"
                          call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
                          ]])],
          [
              AC_MSG_NOTICE([PETSc modules are working.])
              AC_DEFINE(HAVE_PETSC_MODULES,1,[Define if you have petsc fortran modules.] )
              FCFLAGS="$FCFLAGS -DHAVE_PETSC_MODULES"
          ],
          [
              AC_MSG_NOTICE([PETSc modules don't work, using headers instead.])
          ])
elif test "0$PETSC_VERSION_MINOR" -ge 8; then
  AC_MSG_ERROR([PETSc v3.8 and newer requires petsc fortran modules. Do not use the --disable-petsc-fortran-modules option])
fi

# now try a more realistic program, it's a stripped down
# petsc tutorial - using the headers in the same way as we do in the code
AC_LINK_IFELSE(
[AC_LANG_SOURCE([
program test_petsc
#ifdef HAVE_PETSC_MODULES
  use petsc
#endif
implicit none
#include "petsc_legacy.h"
      double precision  norm
      PetscInt  i,j,II0,m,n,its
      PetscInt, dimension(1) ::  II, JJ
      PetscInt  Istart,Iend,ione
      PetscErrorCode ierr
      PetscBool flg
      PetscScalar, dimension(1) :: v
      PetscScalar one,neg_one
      Vec         x,b,u
      Mat         A 
      KSP         ksp

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      m = 3
      n = 3
      one  = 1.0
      neg_one = -1.0
      ione    = 1
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,'-m',m,flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,'-n',n,flg,ierr)

      call MatCreate(PETSC_COMM_WORLD,A,ierr)
      call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m*n,m*n,ierr)
      call MatSetFromOptions(A,ierr)

      call MatGetOwnershipRange(A,Istart,Iend,ierr)
      
      do 10, II0=Istart,Iend-1
        II = II0
        v = -1.0
        i = II0/n
        j = II0 - i*n  
        if (i.gt.0) then
          JJ = II - n
          call MatSetValues(A,ione,II,ione,JJ,v,INSERT_VALUES,ierr)
        endif
        if (i.lt.m-1) then
          JJ = II + n
          call MatSetValues(A,ione,II,ione,JJ,v,INSERT_VALUES,ierr)
        endif
        if (j.gt.0) then
          JJ = II - 1
          call MatSetValues(A,ione,II,ione,JJ,v,INSERT_VALUES,ierr)
        endif
        if (j.lt.n-1) then
          JJ = II + 1
          call MatSetValues(A,ione,II,ione,JJ,v,INSERT_VALUES,ierr)
        endif
        v = 4.0
        call  MatSetValues(A,ione,II,ione,II,v,INSERT_VALUES,ierr)
 10   continue

      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

      call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,m*n,u,ierr)
      call VecSetFromOptions(u,ierr)
      call VecDuplicate(u,b,ierr)
      call VecDuplicate(b,x,ierr)

      call VecSet(u,one,ierr)
      call MatMult(A,u,b,ierr)

      call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
      call KSPSetOperators(ksp,A,A,ierr)
      call KSPSetFromOptions(ksp,ierr)

      call KSPSolve(ksp,b,x,ierr)

      call VecAXPY(x,neg_one,u,ierr)
      call VecNorm(x,NORM_2,norm,ierr)
      call KSPGetIterationNumber(ksp,its,ierr)

      call KSPDestroy(ksp,ierr)
      call VecDestroy(u,ierr)
      call VecDestroy(x,ierr)
      call VecDestroy(b,ierr)
      call MatDestroy(A,ierr)

      call PetscFinalize(ierr)
end program test_petsc
])],
[
AC_MSG_NOTICE([PETSc program succesfully compiled and linked.])
],
[
cp conftest.F90 test_petsc.F90
AC_MSG_FAILURE([Failed to compile and link PETSc program.])])

FCFLAGS="$SAVE_FCFLAGS"
AC_LANG_RESTORE


AC_DEFINE(HAVE_PETSC,1,[Define if you have the PETSc library.])

])dnl ACX_PETSc

m4_include(m4/ACX_lib_automagic.m4)

dnl ----------------------------------------------------------------------------
dnl check for the optional hypre library (linked in with PETSc)
dnl ----------------------------------------------------------------------------
AC_DEFUN([ACX_hypre], [
AC_REQUIRE([ACX_PETSc])

# Ensure the comiler finds the library...
tmpLIBS=$LIBS
tmpCPPFLAGS=$CPPFLAGS
AC_LANG_SAVE
AC_LANG([Fortran])
AC_SEARCH_LIBS(
	[PCHYPRESetType],
	[HYPRE],
	[AC_DEFINE(HAVE_HYPRE,1,[Define if you have hypre library.])],)
# Save variables...
AC_LANG_RESTORE
LIBS=$tmpLIBS
CPPFLAGS=$tmpCPPFLAGS
])dnl ACX_hypre

# ===========================================================================
#     https://www.gnu.org/software/autoconf-archive/ax_python_devel.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_PYTHON_DEVEL([version])
#
# DESCRIPTION
#
#   Note: Defines as a precious variable "PYTHON_VERSION". Don't override it
#   in your configure.ac.
#
#   This macro checks for Python and tries to get the include path to
#   'Python.h'. It provides the $(PYTHON_CPPFLAGS) and $(PYTHON_LIBS) output
#   variables. It also exports $(PYTHON_EXTRA_LIBS) and
#   $(PYTHON_EXTRA_LDFLAGS) for embedding Python in your code.
#
#   You can search for some particular version of Python by passing a
#   parameter to this macro, for example ">= '2.3.1'", or "== '2.4'". Please
#   note that you *have* to pass also an operator along with the version to
#   match, and pay special attention to the single quotes surrounding the
#   version number. Don't use "PYTHON_VERSION" for this: that environment
#   variable is declared as precious and thus reserved for the end-user.
#
#   This macro should work for all versions of Python >= 2.1.0. As an end
#   user, you can disable the check for the python version by setting the
#   PYTHON_NOVERSIONCHECK environment variable to something else than the
#   empty string.
#
#   If you need to use this macro for an older Python version, please
#   contact the authors. We're always open for feedback.
#
# LICENSE
#
#   Copyright (c) 2009 Sebastian Huber <sebastian-huber@web.de>
#   Copyright (c) 2009 Alan W. Irwin
#   Copyright (c) 2009 Rafael Laboissiere <rafael@laboissiere.net>
#   Copyright (c) 2009 Andrew Collier
#   Copyright (c) 2009 Matteo Settenvini <matteo@member.fsf.org>
#   Copyright (c) 2009 Horst Knorr <hk_classes@knoda.org>
#   Copyright (c) 2013 Daniel Mullner <muellner@math.stanford.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 21

AU_ALIAS([AC_PYTHON_DEVEL], [AX_PYTHON_DEVEL])
AC_DEFUN([AX_PYTHON_DEVEL],[
	#
	# Allow the use of a (user set) custom python version
	#
	AC_ARG_VAR([PYTHON_VERSION],[The installed Python
		version to use, for example '2.3'. This string
		will be appended to the Python interpreter
		canonical name.])

	AC_PATH_PROG([PYTHON],[python[$PYTHON_VERSION]])
	if test -z "$PYTHON"; then
	   AC_MSG_ERROR([Cannot find python$PYTHON_VERSION in your system path])
	   PYTHON_VERSION=""
	fi

	#
	# Check for a version of Python >= 2.1.0
	#
	AC_MSG_CHECKING([for a version of Python >= '2.1.0'])
	ac_supports_python_ver=`$PYTHON -c "import sys; \
		ver = sys.version.split ()[[0]]; \
		print (ver >= '2.1.0')"`
	if test "$ac_supports_python_ver" != "True"; then
		if test -z "$PYTHON_NOVERSIONCHECK"; then
			AC_MSG_RESULT([no])
			AC_MSG_FAILURE([
This version of the AC@&t@_PYTHON_DEVEL macro
doesn't work properly with versions of Python before
2.1.0. You may need to re-run configure, setting the
variables PYTHON_CPPFLAGS, PYTHON_LIBS, PYTHON_SITE_PKG,
PYTHON_EXTRA_LIBS and PYTHON_EXTRA_LDFLAGS by hand.
Moreover, to disable this check, set PYTHON_NOVERSIONCHECK
to something else than an empty string.
])
		else
			AC_MSG_RESULT([skip at user request])
		fi
	else
		AC_MSG_RESULT([yes])
	fi

	#
	# if the macro parameter ``version'' is set, honour it
	#
	if test -n "$1"; then
		AC_MSG_CHECKING([for a version of Python $1])
		ac_supports_python_ver=`$PYTHON -c "import sys; \
			ver = sys.version.split ()[[0]]; \
			print (ver $1)"`
		if test "$ac_supports_python_ver" = "True"; then
		   AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
			AC_MSG_ERROR([this package requires Python $1.
If you have it installed, but it isn't the default Python
interpreter in your system path, please pass the PYTHON_VERSION
variable to configure. See ``configure --help'' for reference.
])
			PYTHON_VERSION=""
		fi
	fi

	#
	# Check if you have distutils, else fail
	#
	AC_MSG_CHECKING([for the distutils Python package])
	ac_distutils_result=`$PYTHON -c "import distutils" 2>&1`
	if test $? -eq 0; then
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
		AC_MSG_ERROR([cannot import Python module "distutils".
Please check your Python installation. The error was:
$ac_distutils_result])
		PYTHON_VERSION=""
	fi

	#
	# Check for Python include path
	#
	AC_MSG_CHECKING([for Python include path])
	if test -z "$PYTHON_CPPFLAGS"; then
		python_path=`$PYTHON -c "import distutils.sysconfig; \
			print (distutils.sysconfig.get_python_inc ());"`
		plat_python_path=`$PYTHON -c "import distutils.sysconfig; \
			print (distutils.sysconfig.get_python_inc (plat_specific=1));"`
		if test -n "${python_path}"; then
			if test "${plat_python_path}" != "${python_path}"; then
				python_path="-I$python_path -I$plat_python_path"
			else
				python_path="-I$python_path"
			fi
		fi
		PYTHON_CPPFLAGS=$python_path
	fi
	AC_MSG_RESULT([$PYTHON_CPPFLAGS])
	AC_SUBST([PYTHON_CPPFLAGS])

	#
	# Check for Python library path
	#
	AC_MSG_CHECKING([for Python library path])
	if test -z "$PYTHON_LIBS"; then
		# (makes two attempts to ensure we've got a version number
		# from the interpreter)
		ac_python_version=`cat<<EOD | $PYTHON -

# join all versioning strings, on some systems
# major/minor numbers could be in different list elements
from distutils.sysconfig import *
e = get_config_var('VERSION')
if e is not None:
	print(e)
EOD`

		if test -z "$ac_python_version"; then
			if test -n "$PYTHON_VERSION"; then
				ac_python_version=$PYTHON_VERSION
			else
				ac_python_version=`$PYTHON -c "import sys; \
					print (sys.version[[:3]])"`
			fi
		fi

		# Make the versioning information available to the compiler
		AC_DEFINE_UNQUOTED([HAVE_PYTHON], ["$ac_python_version"],
                                   [If available, contains the Python version number currently in use.])

		# First, the library directory:
		ac_python_libdir=`cat<<EOD | $PYTHON -

# There should be only one
import distutils.sysconfig
e = distutils.sysconfig.get_config_var('LIBDIR')
if e is not None:
	print (e)
EOD`

		# Now, for the library:
		ac_python_library=`cat<<EOD | $PYTHON -

import distutils.sysconfig
c = distutils.sysconfig.get_config_vars()
if 'LDVERSION' in c:
	print ('python'+c[['LDVERSION']])
else:
	print ('python'+c[['VERSION']])
EOD`

		# This small piece shamelessly adapted from PostgreSQL python macro;
		# credits goes to momjian, I think. I'd like to put the right name
		# in the credits, if someone can point me in the right direction... ?
		#
		if test -n "$ac_python_libdir" -a -n "$ac_python_library"
		then
			# use the official shared library
			ac_python_library=`echo "$ac_python_library" | sed "s/^lib//"`
			PYTHON_LIBS="-L$ac_python_libdir -l$ac_python_library"
		else
			# old way: use libpython from python_configdir
			ac_python_libdir=`$PYTHON -c \
			  "from distutils.sysconfig import get_python_lib as f; \
			  import os; \
			  print (os.path.join(f(plat_specific=1, standard_lib=1), 'config'));"`
			PYTHON_LIBS="-L$ac_python_libdir -lpython$ac_python_version"
		fi

		if test -z "PYTHON_LIBS"; then
			AC_MSG_ERROR([
  Cannot determine location of your Python DSO. Please check it was installed with
  dynamic libraries enabled, or try setting PYTHON_LIBS by hand.
			])
		fi
	fi
	AC_MSG_RESULT([$PYTHON_LIBS])
	AC_SUBST([PYTHON_LIBS])

	#
	# Check for site packages
	#
	AC_MSG_CHECKING([for Python site-packages path])
	if test -z "$PYTHON_SITE_PKG"; then
		PYTHON_SITE_PKG=`$PYTHON -c "import distutils.sysconfig; \
			print (distutils.sysconfig.get_python_lib(0,0));"`
	fi
	AC_MSG_RESULT([$PYTHON_SITE_PKG])
	AC_SUBST([PYTHON_SITE_PKG])

	#
	# libraries which must be linked in when embedding
	#
	AC_MSG_CHECKING(python extra libraries)
	if test -z "$PYTHON_EXTRA_LIBS"; then
	   PYTHON_EXTRA_LIBS=`$PYTHON -c "import distutils.sysconfig; \
                conf = distutils.sysconfig.get_config_var; \
                print (conf('LIBS') + ' ' + conf('SYSLIBS'))"`
	fi
	AC_MSG_RESULT([$PYTHON_EXTRA_LIBS])
	AC_SUBST(PYTHON_EXTRA_LIBS)

	#
	# linking flags needed when embedding
	#
	AC_MSG_CHECKING(python extra linking flags)
	if test -z "$PYTHON_EXTRA_LDFLAGS"; then
		PYTHON_EXTRA_LDFLAGS=`$PYTHON -c "import distutils.sysconfig; \
			conf = distutils.sysconfig.get_config_var; \
			print ('' if conf('PYTHONFRAMEWORK') else conf('LINKFORSHARED'))"`
	fi
	AC_MSG_RESULT([$PYTHON_EXTRA_LDFLAGS])
	AC_SUBST(PYTHON_EXTRA_LDFLAGS)

	#
	# final check to see if everything compiles alright
	#
	AC_MSG_CHECKING([consistency of all components of python development environment])
	# save current global flags
	ac_save_LIBS="$LIBS"
	ac_save_LDFLAGS="$LDFLAGS"
	ac_save_CPPFLAGS="$CPPFLAGS"
	LIBS="$ac_save_LIBS $PYTHON_LIBS $PYTHON_EXTRA_LIBS $PYTHON_EXTRA_LIBS"
	LDFLAGS="$ac_save_LDFLAGS $PYTHON_EXTRA_LDFLAGS"
	CPPFLAGS="$ac_save_CPPFLAGS $PYTHON_CPPFLAGS"
	AC_LANG_PUSH([C])
	AC_LINK_IFELSE([
		AC_LANG_PROGRAM([[#include <Python.h>]],
				[[Py_Initialize();]])
		],[pythonexists=yes],[pythonexists=no])
	AC_LANG_POP([C])
	# turn back to default flags
	CPPFLAGS="$ac_save_CPPFLAGS"
	LIBS="$ac_save_LIBS"
	LDFLAGS="$ac_save_LDFLAGS"

	AC_MSG_RESULT([$pythonexists])

        if test ! "x$pythonexists" = "xyes"; then
	   AC_MSG_FAILURE([
  Could not link test program to Python. Maybe the main Python library has been
  installed in some non-standard library path. If so, pass it to configure,
  via the LIBS environment variable.
  Example: ./configure LIBS="-L/usr/non-standard-path/python/lib"
  ============================================================================
   ERROR!
   You probably have to install the development version of the Python package
   for your distribution.  The exact name of this package varies among them.
  ============================================================================
	   ])
	  PYTHON_VERSION=""
	fi

	#
	# all done!
	#
])

# ===========================================================================
#       https://www.gnu.org/software/autoconf-archive/ax_lib_hdf5.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_HDF5([serial/parallel])
#
# DESCRIPTION
#
#   This macro provides tests of the availability of HDF5 library.
#
#   The optional macro argument should be either 'serial' or 'parallel'. The
#   former only looks for serial HDF5 installations via h5cc. The latter
#   only looks for parallel HDF5 installations via h5pcc. If the optional
#   argument is omitted, serial installations will be preferred over
#   parallel ones.
#
#   The macro adds a --with-hdf5 option accepting one of three values:
#
#     no   - do not check for the HDF5 library.
#     yes  - do check for HDF5 library in standard locations.
#     path - complete path to the HDF5 helper script h5cc or h5pcc.
#
#   If HDF5 is successfully found, this macro calls
#
#     AC_SUBST(HDF5_VERSION)
#     AC_SUBST(HDF5_CC)
#     AC_SUBST(HDF5_CFLAGS)
#     AC_SUBST(HDF5_CPPFLAGS)
#     AC_SUBST(HDF5_LDFLAGS)
#     AC_SUBST(HDF5_LIBS)
#     AC_SUBST(HDF5_FC)
#     AC_SUBST(HDF5_FFLAGS)
#     AC_SUBST(HDF5_FLIBS)
#     AC_SUBST(HDF5_TYPE)
#     AC_DEFINE(HAVE_HDF5)
#
#   and sets with_hdf5="yes".  Additionally, the macro sets
#   with_hdf5_fortran="yes" if a matching Fortran wrapper script is found.
#   Note that Autoconf's Fortran support is not used to perform this check.
#   H5CC and H5FC will contain the appropriate serial or parallel HDF5
#   wrapper script locations.
#
#   If HDF5 is disabled or not found, this macros sets with_hdf5="no" and
#   with_hdf5_fortran="no".
#
#   Your configuration script can test $with_hdf to take any further
#   actions. HDF5_{C,CPP,LD}FLAGS may be used when building with C or C++.
#   HDF5_F{FLAGS,LIBS} should be used when building Fortran applications.
#
#   To use the macro, one would code one of the following in "configure.ac"
#   before AC_OUTPUT:
#
#     1) dnl Check for HDF5 support
#        AX_LIB_HDF5()
#
#     2) dnl Check for serial HDF5 support
#        AX_LIB_HDF5([serial])
#
#     3) dnl Check for parallel HDF5 support
#        AX_LIB_HDF5([parallel])
#
#   One could test $with_hdf5 for the outcome or display it as follows
#
#     echo "HDF5 support:  $with_hdf5"
#
#   You could also for example, override the default CC in "configure.ac" to
#   enforce compilation with the compiler that HDF5 uses:
#
#     AX_LIB_HDF5([parallel])
#     if test "$with_hdf5" = "yes"; then
#             CC="$HDF5_CC"
#     else
#             AC_MSG_ERROR([Unable to find HDF5, we need parallel HDF5.])
#     fi
#
#   The HDF5_TYPE environment variable returns "parallel" or "serial",
#   depending on which type of library is found.
#
# LICENSE
#
#   Copyright (c) 2009 Timothy Brown <tbrown@freeshell.org>
#   Copyright (c) 2010 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 19

AC_DEFUN([AX_LIB_HDF5], [

AC_REQUIRE([AC_PROG_SED])
AC_REQUIRE([AC_PROG_AWK])
AC_REQUIRE([AC_PROG_GREP])

dnl Check first argument is one of the recognized values.
dnl Fail eagerly if is incorrect as this simplifies case statements below.
if   test "m4_normalize(m4_default([$1],[]))" = ""        ; then
    : # Recognized value
elif test "m4_normalize(m4_default([$1],[]))" = "serial"  ; then
    : # Recognized value
elif test "m4_normalize(m4_default([$1],[]))" = "parallel"; then
    : # Recognized value
else
    AC_MSG_ERROR([
Unrecognized value for AX[]_LIB_HDF5 within configure.ac.
If supplied, argument 1 must be either 'serial' or 'parallel'.
])
fi

dnl Add a default --with-hdf5 configuration option.
AC_ARG_WITH([hdf5],
  AS_HELP_STRING(
    [--with-hdf5=[yes/no/PATH]],
    m4_case(m4_normalize([$1]),
            [serial],   [location of h5cc for serial HDF5 configuration],
            [parallel], [location of h5pcc for parallel HDF5 configuration],
            [location of h5cc or h5pcc for HDF5 configuration])
  ),
  [if test "$withval" = "no"; then
     with_hdf5="no"
   elif test "$withval" = "yes"; then
     with_hdf5="yes"
   else
     with_hdf5="yes"
     H5CC="$withval"
   fi],
   [with_hdf5="yes"]
)

dnl Set defaults to blank
HDF5_CC=""
HDF5_VERSION=""
HDF5_CFLAGS=""
HDF5_CPPFLAGS=""
HDF5_LDFLAGS=""
HDF5_LIBS=""
HDF5_FC=""
HDF5_FFLAGS=""
HDF5_FLIBS=""
HDF5_TYPE=""

dnl Try and find hdf5 compiler tools and options.
if test "$with_hdf5" = "yes"; then
    if test -z "$H5CC"; then
        dnl Check to see if H5CC is in the path.
        AC_PATH_PROGS(
            [H5CC],
            m4_case(m4_normalize([$1]),
                [serial],   [h5cc],
                [parallel], [h5pcc],
                [h5cc h5pcc]),
            [])
    else
        AC_MSG_CHECKING([Using provided HDF5 C wrapper])
        AC_MSG_RESULT([$H5CC])
    fi
    AC_MSG_CHECKING([for HDF5 type])
    AS_CASE([$H5CC],
        [*h5pcc], [HDF5_TYPE=parallel],
        [*h5cc], [HDF5_TYPE=serial],
        [HDF5_TYPE=neither])
    AC_MSG_RESULT([$HDF5_TYPE])
    AC_MSG_CHECKING([for HDF5 libraries])
    if test ! -f "$H5CC" || test ! -x "$H5CC"; then
        AC_MSG_RESULT([no])
        AC_MSG_WARN(m4_case(m4_normalize([$1]),
            [serial],  [
Unable to locate serial HDF5 compilation helper script 'h5cc'.
Please specify --with-hdf5=<LOCATION> as the full path to h5cc.
HDF5 support is being disabled (equivalent to --with-hdf5=no).
],            [parallel],[
Unable to locate parallel HDF5 compilation helper script 'h5pcc'.
Please specify --with-hdf5=<LOCATION> as the full path to h5pcc.
HDF5 support is being disabled (equivalent to --with-hdf5=no).
],            [
Unable to locate HDF5 compilation helper scripts 'h5cc' or 'h5pcc'.
Please specify --with-hdf5=<LOCATION> as the full path to h5cc or h5pcc.
HDF5 support is being disabled (equivalent to --with-hdf5=no).
]))
        with_hdf5="no"
        with_hdf5_fortran="no"
    else
        dnl Get the h5cc output
        HDF5_SHOW=$(eval $H5CC -show)

        dnl Get the actual compiler used
        HDF5_CC=$(eval $H5CC -show | $AWK '{print $[]1}')
        if test "$HDF5_CC" = "ccache"; then
            HDF5_CC=$(eval $H5CC -show | $AWK '{print $[]2}')
        fi

        dnl h5cc provides both AM_ and non-AM_ options
        dnl depending on how it was compiled either one of
        dnl these are empty. Lets roll them both into one.

        dnl Look for "HDF5 Version: X.Y.Z"
        HDF5_VERSION=$(eval $H5CC -showconfig | $GREP 'HDF5 Version:' \
            | $AWK '{print $[]3}')

        dnl A ideal situation would be where everything we needed was
        dnl in the AM_* variables. However most systems are not like this
        dnl and seem to have the values in the non-AM variables.
        dnl
        dnl We try the following to find the flags:
        dnl (1) Look for "NAME:" tags
        dnl (2) Look for "H5_NAME:" tags
        dnl (3) Look for "AM_NAME:" tags
        dnl
        HDF5_tmp_flags=$(eval $H5CC -showconfig \
            | $GREP 'FLAGS\|Extra libraries:' \
            | $AWK -F: '{printf("%s "), $[]2}' )

        dnl Find the installation directory and append include/
        HDF5_tmp_inst=$(eval $H5CC -showconfig \
            | $GREP 'Installation point:' \
            | $AWK '{print $[]NF}' )

        dnl Add this to the CPPFLAGS
        HDF5_CPPFLAGS="-I${HDF5_tmp_inst}/include"

        dnl Now sort the flags out based upon their prefixes
        for arg in $HDF5_SHOW $HDF5_tmp_flags ; do
          case "$arg" in
            -I*) echo $HDF5_CPPFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || HDF5_CPPFLAGS="$HDF5_CPPFLAGS $arg"
              ;;
            -L*) echo $HDF5_LDFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || HDF5_LDFLAGS="$HDF5_LDFLAGS $arg"
              ;;
            -l*) echo $HDF5_LIBS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || HDF5_LIBS="$HDF5_LIBS $arg"
              ;;
          esac
        done

        HDF5_LIBS="-lhdf5 $HDF5_LIBS"
        AC_MSG_RESULT([yes (version $[HDF5_VERSION])])

        dnl See if we can compile
        AC_LANG_PUSH([C])
        ax_lib_hdf5_save_CC=$CC
        ax_lib_hdf5_save_CPPFLAGS=$CPPFLAGS
        ax_lib_hdf5_save_LIBS=$LIBS
        ax_lib_hdf5_save_LDFLAGS=$LDFLAGS
        CC=$HDF5_CC
        CPPFLAGS=$HDF5_CPPFLAGS
        LIBS=$HDF5_LIBS
        LDFLAGS=$HDF5_LDFLAGS
        AC_CHECK_HEADER([hdf5.h], [ac_cv_hadf5_h=yes], [ac_cv_hadf5_h=no])
        AC_CHECK_LIB([hdf5], [H5Fcreate], [ac_cv_libhdf5=yes],
                     [ac_cv_libhdf5=no])
        if test "$ac_cv_hadf5_h" = "no" && test "$ac_cv_libhdf5" = "no" ; then
          AC_MSG_WARN([Unable to compile HDF5 test program])
        fi
        dnl Look for HDF5's high level library
        AC_HAVE_LIBRARY([hdf5_hl], [HDF5_LIBS="-lhdf5_hl $HDF5_LIBS"], [], [])

        CC=$ax_lib_hdf5_save_CC
        CPPFLAGS=$ax_lib_hdf5_save_CPPFLAGS
        LIBS=$ax_lib_hdf5_save_LIBS
        LDFLAGS=$ax_lib_hdf5_save_LDFLAGS
        AC_LANG_POP([C])

        AC_MSG_CHECKING([for matching HDF5 Fortran wrapper])
        dnl Presume HDF5 Fortran wrapper is just a name variant from H5CC
        H5FC=$(eval echo -n $H5CC | $SED -n 's/cc$/fc/p')
        if test -x "$H5FC"; then
            AC_MSG_RESULT([$H5FC])
            with_hdf5_fortran="yes"
            AC_SUBST([H5FC])

            dnl Again, pry any remaining -Idir/-Ldir from compiler wrapper
            for arg in `$H5FC -show`
            do
              case "$arg" in #(
                -I*) echo $HDF5_FFLAGS | $GREP -e "$arg" >/dev/null \
                      || HDF5_FFLAGS="$HDF5_FFLAGS $arg"
                  ;;#(
                -L*) echo $HDF5_FFLAGS | $GREP -e "$arg" >/dev/null \
                      || HDF5_FFLAGS="$HDF5_FFLAGS $arg"
                     dnl HDF5 installs .mod files in with libraries,
                     dnl but some compilers need to find them with -I
                     echo $HDF5_FFLAGS | $GREP -e "-I${arg#-L}" >/dev/null \
                      || HDF5_FFLAGS="$HDF5_FFLAGS -I${arg#-L}"
                  ;;
              esac
            done

            dnl Make Fortran link line by inserting Fortran libraries
            for arg in $HDF5_LIBS
            do
              case "$arg" in #(
                -lhdf5_hl) HDF5_FLIBS="$HDF5_FLIBS -lhdf5hl_fortran $arg"
                  ;; #(
                -lhdf5)    HDF5_FLIBS="$HDF5_FLIBS -lhdf5_fortran $arg"
                  ;; #(
                *) HDF5_FLIBS="$HDF5_FLIBS $arg"
                  ;;
              esac
            done
        else
            AC_MSG_RESULT([no])
            with_hdf5_fortran="no"
        fi

	AC_SUBST([HDF5_VERSION])
	AC_SUBST([HDF5_CC])
	AC_SUBST([HDF5_CFLAGS])
	AC_SUBST([HDF5_CPPFLAGS])
	AC_SUBST([HDF5_LDFLAGS])
	AC_SUBST([HDF5_LIBS])
	AC_SUBST([HDF5_FC])
	AC_SUBST([HDF5_FFLAGS])
	AC_SUBST([HDF5_FLIBS])
	AC_SUBST([HDF5_TYPE])
	AC_DEFINE([HAVE_HDF5], [1], [Defined if you have HDF5 support])
    fi
fi
])

AC_DEFUN([ACX_adjoint], [
# Set variables...
AC_ARG_WITH(
	[adjoint],
	[  --with-adjoint=prefix        Prefix where libadjoint is installed],
	[adjoint="$withval"],
    [])

bakLIBS=$LIBS
tmpLIBS=$LIBS
tmpCPPFLAGS=$CPPFLAGS
if test "$adjoint" != "no"; then
  if test -d "$adjoint" ; then
    adjoint_LIBS_PATH="$adjoint/lib"
    adjoint_INCLUDES_PATH="$adjoint/include"
    # Ensure the comiler finds the library...
    tmpLIBS="$tmpLIBS -L$adjoint/lib"
    tmpCPPFLAGS="$tmpCPPFLAGS  -I$adjoint/include -I$adjoint/include/libadjoint"
  fi
  if [[ -f /usr/lib/libadjoint.a -o -f /usr/local/lib/libadjoint.a ]] ; then
    tmpLIBS="$tmpLIBS -L/usr/lib -L/usr/local/lib/ -ladjoint"
    tmpCPPFLAGS="$tmpCPPFLAGS -I/usr/include/ -I/usr/local/include/ -I/usr/include/libadjoint -I/usr/local/include/libadjoint"
  fi
fi
LIBS=$tmpLIBS
CPPFLAGS=$tmpCPPFLAGS

AC_LANG_SAVE
AC_LANG_C
AC_CHECK_LIB(
	[adjoint],
	[adj_get_adjoint_equation],
	[AC_DEFINE(HAVE_ADJOINT,1,[Define if you have libadjoint.])HAVE_ADJOINT=yes],
	[AC_MSG_WARN( [Could not link in libadjoint ... ] );HAVE_ADJOINT=no;LIBS=$bakLIBS] )
# Save variables...
AC_LANG_RESTORE

])dnl ACX_adjoint
