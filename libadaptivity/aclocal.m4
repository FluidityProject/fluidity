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
dnl @version $Id: aclocal.m4 7146 2008-07-22 22:40:14Z pfarrell $
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
dnl @version $Id: aclocal.m4 7146 2008-07-22 22:40:14Z pfarrell $
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

# add BLAS to libs
LIBS="$BLAS_LIBS $LIBS"

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
AC_REQUIRE([ACX_ParMetis])
BLAS_LIBS="$BLAS_LIBS $FLIBS"
AC_REQUIRE([ACX_LAPACK])
LAPACK_LIBS="$LAPACK_LIBS $BLAS_LIBS"
AC_PATH_XTRA

# Set variables...
AC_ARG_WITH(
	[PETSc],
	[  --with-PETSc=PFX        Prefix where PETSc is installed (use =no to disable)],
	[PETSc="$withval"],
	[if test ! -z "$PETSC_DIR"; then
		PETSc="$PETSC_DIR"
		echo "note: assuming PETSc library is in $PETSc (/lib,/include) as specified by environment variable PETSC_DIR"
	else
		PETSc="/usr/local"
		echo "note: assuming PETSc library is in /usr/local (/lib,/include)"
	fi])
AC_ARG_WITH(
	[PETSc_ARCH],
	[  --with-PETSc_ARCH=VAL   PETSc hardware architecture (optional)],
	[PETSc_ARCH="$withval"],
	[if test ! -z "$PETSC_ARCH"; then
		PETSc_ARCH="$PETSC_ARCH"
		echo "note: assuming PETSc hardware architecture to be $PETSc_ARCH as specified by environment variable PETSC_ARCH"
	else
		PETSc_ARCH=`uname -p`
		echo "note: assuming PETSc hardware architecture to be $PETSc_ARCH"
	fi])
PETSc_LIBS_PATH="$PETSc/lib/$PETSc_ARCH"
PETSc_INCLUDES_PATH="$PETSc/include"

# Check that the compiler uses the library we specified...
if test -e $PETSc_LIBS_PATH/libpetsc.a || test -e $PETSc_LIBS_PATH/libpetsc.so; then
	echo "note: using $PETSc_LIBS_PATH/libpetsc (.a/.so)"
else
	AC_MSG_NOTICE( [Could not physically find PETSc library... exiting] )
fi 
if test -e $PETSc_INCLUDES_PATH/petsc.h; then
	echo "note: using $PETSc_INCLUDES_PATH/petsc.h"
else
	AC_MSG_NOTICE( [Could not physically find PETSc header file... exiting] )
fi 

# Ensure the comiler finds the library...
tmpLIBS=$LIBS
tmpCPPFLAGS=$CPPFLAGS
AC_LANG_SAVE
AC_LANG_C
LIBS="-L$PETSc_LIBS_PATH $MPI_LIBS_PATHS $MPI_LIBS $LAPACK_LIBS $LIBS -lm"
CPPFLAGS="-I$PETSc_INCLUDES_PATH -I$PETSc/bmake/$PETSc_ARCH $CPPFLAGS"
FFLAGS="$FFLAGS -I$PETSc_INCLUDES_PATH -I$PETSc/bmake/$PETSc_ARCH"
FCFLAGS="$FCFLAGS -I$PETSc_INCLUDES_PATH -I$PETSc/bmake/$PETSc_ARCH"
AC_CHECK_LIB(
	[petsc],
	[PetscInitializePackage],
	[AC_DEFINE(HAVE_PETSC,1,[Define if you have PETSc library.])],
	[AC_CHECK_LIB(
    [craypetsc],
    [PetscInitializePackage],
    [CRAYPETSC=1; AC_DEFINE(HAVE_PETSC,1,[Define if you have PETSc library.])],
    [AC_MSG_ERROR( [Could not link in the PETSc library... exiting] )] )])
AC_CHECK_HEADER(
	[petsc.h],
	[AC_DEFINE( 
		[HAVE_PETSc_H],,
		[Define to 1 if you have the <petsc.h> header file.])],
	[AC_MSG_ERROR( [Could not compile in the PETSc headers... exiting] )] )
if test -z "$CRAYPETSC"; then
  PETSC_LIBS="-L$PETSc_LIBS_PATH -lpetscksp -lpetscmat -lpetscvec -lpetsc -lparmetis -lmetis"
else
  PETSC_LIBS="-L$PETSc_LIBS_PATH -lcraypetsc"
fi
PETSC_FLAG="-DHAVE_PETSC -I$PETSc_INCLUDES_PATH -I$PETSc/bmake/$PETSc_ARCH"
tmpCPPFLAGS="$tmpCPPFLAGS -DHAVE_PETSC -I$PETSc_INCLUDES_PATH -I$PETSc/bmake/$PETSc_ARCH"
# Save variables...
AC_LANG_RESTORE
LIBS="$PETSC_LIBS $tmpLIBS"
CPPFLAGS=$tmpCPPFLAGS
AC_SUBST(PETSC_LIBS)
AC_SUBST(PETSC_FLAG)
])dnl ACX_PETSc

AC_DEFUN([ACX_ParMetis], [
# Set variables...
AC_ARG_WITH(
	[ParMetis],
	[  --with-ParMetis=PFX        Prefix where ParMetis is installed],
	[ParMetis="$withval"],
    [])
ParMetis_LIBS_PATH="$ParMetis/lib"

# Check that the compiler uses the library we specified...
if test -e $ParMetis_LIBS_PATH/libparmetis.a; then
	echo "note: using $ParMetis_LIBS_PATH/libparmetis.a"
fi 

# Ensure the comiler finds the library...
tmpLIBS=$LIBS
tmpCPPFLAGS=$CPPFLAGS
AC_LANG_SAVE
AC_LANG_C
LIBS="$tmpLIBS -L$ParMetis_LIBS_PATH -lparmetis -lmetis -lm"
AC_CHECK_LIB(
	[parmetis],
	[ParMETIS_V3_AdaptiveRepart],
	[AC_DEFINE(HAVE_PARMETIS,1,[Define if you have ParMetis library.])],
	[AC_MSG_ERROR( [Could not link in the ParMetis library... exiting] )] )
tmpLIBS="$tmpLIBS -L$ParMetis_LIBS_PATH -lparmetis -lmetis"
# Save variables...
AC_LANG_RESTORE
LIBS=$tmpLIBS
CPPFLAGS=$tmpCPPFLAGS
])dnl ACX_ParMetis

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

dnl @synopsis AC_PYTHON_DEVEL([version])
dnl
dnl Note: Defines as a precious variable "PYTHON_VERSION". Don't
dnl override it in your configure.ac.
dnl
dnl This macro checks for Python and tries to get the include path to
dnl 'Python.h'. It provides the $(PYTHON_CPPFLAGS) and
dnl $(PYTHON_LDFLAGS) output variables. It also exports
dnl $(PYTHON_EXTRA_LIBS) and $(PYTHON_EXTRA_LDFLAGS) for embedding
dnl Python in your code.
dnl
dnl You can search for some particular version of Python by passing a
dnl parameter to this macro, for example ">= '2.3.1'", or "== '2.4'".
dnl Please note that you *have* to pass also an operator along with the
dnl version to match, and pay special attention to the single quotes
dnl surrounding the version number. Don't use "PYTHON_VERSION" for
dnl this: that environment variable is declared as precious and thus
dnl reserved for the end-user.
dnl
dnl This macro should work for all versions of Python >= 2.1.0. As an
dnl end user, you can disable the check for the python version by
dnl setting the PYTHON_NOVERSIONCHECK environment variable to something
dnl else than the empty string.
dnl
dnl If you need to use this macro for an older Python version, please
dnl contact the authors. We're always open for feedback.
dnl
dnl @category InstalledPackages
dnl @author Sebastian Huber <sebastian-huber@web.de>
dnl @author Alan W. Irwin <irwin@beluga.phys.uvic.ca>
dnl @author Rafael Laboissiere <laboissiere@psy.mpg.de>
dnl @author Andrew Collier <colliera@nu.ac.za>
dnl @author Matteo Settenvini <matteo@member.fsf.org>
dnl @author Horst Knorr <hk_classes@knoda.org>
dnl @version 2006-05-27
dnl @license GPLWithACException

AC_DEFUN([AC_PYTHON_DEVEL],[
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
	# Save LIBS and CPPFLAGS, since they are 'restored' later on
	# Adrian Umpleby, 2007/11/20
	#
	ac_save_CPPFLAGS="$CPPFLAGS"
	ac_save_LIBS="$LIBS"

	#
	# Check for a version of Python >= 2.3.0
	#
	AC_MSG_CHECKING([for a version of Python >= '2.3.0'])
	ac_supports_python_ver=`$PYTHON -c "import sys, string; \
		ver = string.split(sys.version)[[0]]; \
		print ver >= '2.3.0'"`
	if test "$ac_supports_python_ver" != "True"; then
		if test -z "$PYTHON_NOVERSIONCHECK"; then
			AC_MSG_RESULT([no])
			AC_MSG_FAILURE([
This version of the AC@&t@_PYTHON_DEVEL macro
doesn't work properly with versions of Python before
2.1.0. You may need to re-run configure, setting the
variables PYTHON_CPPFLAGS, PYTHON_LDFLAGS, PYTHON_SITE_PKG,
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
		ac_supports_python_ver=`$PYTHON -c "import sys, string; \
			ver = string.split(sys.version)[[0]]; \
			print ver $1"`
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
	if test -z "$ac_distutils_result"; then
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
           		print distutils.sysconfig.get_python_inc();"`
		if test -n "${python_path}"; then
		   	python_path="-I$python_path"
		fi
		PYTHON_CPPFLAGS=$python_path
	fi
	AC_MSG_RESULT([$PYTHON_CPPFLAGS])
	AC_SUBST([PYTHON_CPPFLAGS])

	#
	# Check for Python library path
	#
	AC_MSG_CHECKING([for Python library path])
	if test -z "$PYTHON_LDFLAGS"; then
		# (makes two attempts to ensure we've got a version number
		# from the interpreter)
		py_version=`$PYTHON -c "from distutils.sysconfig import *; \
			from string import join; \
			print join(get_config_vars('VERSION'))"`
		if test "$py_version" == "[None]"; then
			if test -n "$PYTHON_VERSION"; then
				py_version=$PYTHON_VERSION
			else
				py_version=`$PYTHON -c "import sys; \
					print sys.version[[:3]]"`
			fi
		fi

		PYTHON_LDFLAGS=`$PYTHON -c "from distutils.sysconfig import *; \
			from string import join; \
			print '-L' + get_python_lib(0,1), \
		      	'-lpython';"`$py_version
	fi
	AC_MSG_RESULT([$PYTHON_LDFLAGS])
	AC_SUBST([PYTHON_LDFLAGS])


	#
	# Check for site packages
	#
	AC_MSG_CHECKING([for Python site-packages path])
	if test -z "$PYTHON_SITE_PKG"; then
		PYTHON_SITE_PKG=`$PYTHON -c "import distutils.sysconfig; \
		        print distutils.sysconfig.get_python_lib(0,0);"`
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
                print conf('LOCALMODLIBS'), conf('LIBS')"`
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
			print conf('LINKFORSHARED')"`
	fi
	AC_MSG_RESULT([$PYTHON_EXTRA_LDFLAGS])
	AC_SUBST(PYTHON_EXTRA_LDFLAGS)

	#
	# final check to see if everything compiles alright
	#
	AC_MSG_CHECKING([consistency of all components of python development environment])
	AC_LANG_PUSH([C])
	# save current global flags
	LIBS="$ac_save_LIBS $PYTHON_LDFLAGS"
	CPPFLAGS="$ac_save_CPPFLAGS $PYTHON_CPPFLAGS"
	AC_TRY_LINK([
		#include <Python.h>
	],[
		Py_Initialize();
	],[pythonexists=yes],[pythonexists=no])

	AC_MSG_RESULT([$pythonexists])

        if test ! "$pythonexists" = "yes"; then
	   AC_MSG_WARN([
  Could not link test program to Python. Maybe the main Python library has been
  installed in some non-standard library path. If so, pass it to configure,
  via the LDFLAGS environment variable.
  Example: ./configure LDFLAGS="-L/usr/non-standard-path/python/lib"
  ============================================================================
   WARNING!
   You probably have to install the development version of the Python package
   for your distribution.  The exact name of this package varies among them.
  ============================================================================
	   ])
	  PYTHON_VERSION=""
	fi
	AC_LANG_POP
	# turn back to default flags
	CPPFLAGS="$ac_save_CPPFLAGS"
	LIBS="$ac_save_LIBS"

	#
	# all done!
	#
])
