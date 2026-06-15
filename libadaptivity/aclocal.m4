# ===========================================================================
#         https://www.gnu.org/software/autoconf-archive/ax_blas.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a library that implements the BLAS linear-algebra
#   interface (see http://www.netlib.org/blas/). On success, it sets the
#   BLAS_LIBS output variable to hold the requisite library linkages.
#
#   To link with BLAS, you should link with:
#
#     $BLAS_LIBS $LIBS $FLIBS
#
#   in that order. FLIBS is the output variable of the
#   AC_F77_LIBRARY_LDFLAGS macro (called if necessary by AX_BLAS), and is
#   sometimes necessary in order to link with F77 libraries. Users will also
#   need to use AC_F77_DUMMY_MAIN (see the autoconf manual), for the same
#   reason.
#
#   Many libraries are searched for, from ATLAS to CXML to ESSL. The user
#   may also use --with-blas=<lib> in order to use some specific BLAS
#   library <lib>. In order to link successfully, however, be aware that you
#   will probably need to use the same Fortran compiler (which can be set
#   via the F77 env. var.) as was used to compile the BLAS library.
#
#   ACTION-IF-FOUND is a list of shell commands to run if a BLAS library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_BLAS.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2019 Geoffrey M. Oxberry <goxberry@gmail.com>
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

AU_ALIAS([ACX_BLAS], [AX_BLAS])
AC_DEFUN([AX_BLAS], [
AC_PREREQ([2.55])
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
AC_REQUIRE([AC_CANONICAL_HOST])
ax_blas_ok=no

AC_ARG_WITH(blas,
	[AS_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) ax_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.dylib | *.dylib.* | *.o)
		BLAS_LIBS="$with_blas"
	;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

# Get fortran linker names of BLAS functions to check for.
AC_F77_FUNC(sgemm)
AC_F77_FUNC(dgemm)

ax_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $ax_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
	AC_LINK_IFELSE([AC_LANG_CALL([], [$sgemm])], [ax_blas_ok=yes], [BLAS_LIBS=""])
	AC_MSG_RESULT($ax_blas_ok)
	LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $ax_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_MSG_CHECKING([if $sgemm is being linked in already])
	AC_LINK_IFELSE([AC_LANG_CALL([], [$sgemm])], [ax_blas_ok=yes])
	AC_MSG_RESULT($ax_blas_ok)
	LIBS="$save_LIBS"
fi

# BLAS linked to by flexiblas? (Default on FC33+ and RHEL9+)

if test $ax_blas_ok = no; then
	AC_CHECK_LIB(flexiblas, $sgemm, [ax_blas_ok=yes
			                BLAS_LIBS="-lflexiblas"])
fi

# BLAS in OpenBLAS library? (http://xianyi.github.com/OpenBLAS/)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(openblas, $sgemm, [ax_blas_ok=yes
			                BLAS_LIBS="-lopenblas"])
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, $sgemm,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[ax_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(dgemm, $dgemm,
		[AC_CHECK_LIB(sgemm, $sgemm,
			[ax_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi

# BLAS in Intel MKL library?
if test $ax_blas_ok = no; then
	# MKL for gfortran
	if test x"$ac_cv_fc_compiler_gnu" = xyes; then
		# 64 bit
		if test $host_cpu = x86_64; then
			AC_CHECK_LIB(mkl_gf_lp64, $sgemm,
			[ax_blas_ok=yes;BLAS_LIBS="-lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread"],,
			[-lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread])
		# 32 bit
		elif test $host_cpu = i686; then
			AC_CHECK_LIB(mkl_gf, $sgemm,
				[ax_blas_ok=yes;BLAS_LIBS="-lmkl_gf -lmkl_sequential -lmkl_core -lpthread"],,
				[-lmkl_gf -lmkl_sequential -lmkl_core -lpthread])
		fi
	# MKL for other compilers (Intel, PGI, ...?)
	else
		# 64-bit
		if test $host_cpu = x86_64; then
			AC_CHECK_LIB(mkl_intel_lp64, $sgemm,
				[ax_blas_ok=yes;BLAS_LIBS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread"],,
				[-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread])
		# 32-bit
		elif test $host_cpu = i686; then
			AC_CHECK_LIB(mkl_intel, $sgemm,
				[ax_blas_ok=yes;BLAS_LIBS="-lmkl_intel -lmkl_sequential -lmkl_core -lpthread"],,
				[-lmkl_intel -lmkl_sequential -lmkl_core -lpthread])
		fi
	fi
fi
# Old versions of MKL
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(mkl, $sgemm, [ax_blas_ok=yes;BLAS_LIBS="-lmkl -lguide -lpthread"],,[-lguide -lpthread])
fi

# BLAS in Apple vecLib library?
if test $ax_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="-framework Accelerate $LIBS"
	AC_MSG_CHECKING([for $sgemm in -framework Accelerate])
	AC_LINK_IFELSE([AC_LANG_CALL([], [$sgemm])], [ax_blas_ok=yes;BLAS_LIBS="-framework Accelerate"])
	AC_MSG_RESULT($ax_blas_ok)
	LIBS="$save_LIBS"
fi

# BLAS in Alpha CXML library?
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(cxml, $sgemm, [ax_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(dxml, $sgemm, [ax_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $ax_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $sgemm,
				[BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 ax_blas_ok=yes],[],[-lsunmath])])
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(scs, $sgemm, [ax_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $sgemm,
		     [ax_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(essl, $sgemm,
			[ax_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
fi

# Generic BLAS library?
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm, [ax_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

AC_SUBST(BLAS_LIBS)

LIBS="$ax_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        ax_blas_ok=no
        $2
fi
])dnl AX_BLAS
# ===========================================================================
#        https://www.gnu.org/software/autoconf-archive/ax_lapack.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a library that implements the LAPACK linear-algebra
#   interface (see http://www.netlib.org/lapack/). On success, it sets the
#   LAPACK_LIBS output variable to hold the requisite library linkages.
#
#   To link with LAPACK, you should link with:
#
#     $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
#
#   in that order. BLAS_LIBS is the output variable of the AX_BLAS macro,
#   called automatically. FLIBS is the output variable of the
#   AC_F77_LIBRARY_LDFLAGS macro (called if necessary by AX_BLAS), and is
#   sometimes necessary in order to link with F77 libraries. Users will also
#   need to use AC_F77_DUMMY_MAIN (see the autoconf manual), for the same
#   reason.
#
#   The user may also use --with-lapack=<lib> in order to use some specific
#   LAPACK library <lib>. In order to link successfully, however, be aware
#   that you will probably need to use the same Fortran compiler (which can
#   be set via the F77 env. var.) as was used to compile the LAPACK and BLAS
#   libraries.
#
#   ACTION-IF-FOUND is a list of shell commands to run if a LAPACK library
#   is found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it
#   is not found. If ACTION-IF-FOUND is not specified, the default action
#   will define HAVE_LAPACK.
#
# LICENSE
#
#   Copyright (c) 2009 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2019 Geoffrey M. Oxberry <goxberry@gmail.com>
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

#serial 10

AU_ALIAS([ACX_LAPACK], [AX_LAPACK])
AC_DEFUN([AX_LAPACK], [
AC_REQUIRE([AX_BLAS])
ax_lapack_ok=no

AC_ARG_WITH(lapack,
        [AS_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])
case $with_lapack in
        yes | "") ;;
        no) ax_lapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.dylib | *.dylib.* | *.o)
                 LAPACK_LIBS="$with_lapack"
        ;;
        *) LAPACK_LIBS="-l$with_lapack" ;;
esac

# Get fortran linker name of LAPACK function to check for.
AC_F77_FUNC(cheev)

# We cannot use LAPACK if BLAS is not found
if test "x$ax_blas_ok" != xyes; then
        ax_lapack_ok=noblas
        LAPACK_LIBS=""
fi

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
        AC_MSG_CHECKING([for $cheev in $LAPACK_LIBS])
        AC_LINK_IFELSE([AC_LANG_CALL([], [$cheev])], [ax_lapack_ok=yes], [LAPACK_LIBS=""])
        AC_MSG_RESULT($ax_lapack_ok)
        LIBS="$save_LIBS"
        if test $ax_lapack_ok = no; then
                LAPACK_LIBS=""
        fi
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $ax_lapack_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
        AC_CHECK_FUNC($cheev, [ax_lapack_ok=yes])
        LIBS="$save_LIBS"
fi

# Generic LAPACK library?
for lapack in lapack lapack_rs6k; do
        if test $ax_lapack_ok = no; then
                save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
                AC_CHECK_LIB($lapack, $cheev,
                    [ax_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [], [$FLIBS])
                LIBS="$save_LIBS"
        fi
done

AC_SUBST(LAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
        ax_lapack_ok=no
        $2
fi
])dnl AX_LAPACK

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
