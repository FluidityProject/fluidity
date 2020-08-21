#
#  Copyright (c) 2006-2018, The Regents of the University of California,
#  through Lawrence Berkeley National Laboratory (subject to receipt of any
#  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
#  Institut (Switzerland).  All rights reserved.
#
#  License: see file COPYING in top level of source distribution.
#
AC_INIT([H5hut], [2.0.0rc6], [h5part@lists.psi.ch], H5hut)
AC_PREREQ(2.60)
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

ENABLE_DEBUG='no'
ENABLE_C='yes'
ENABLE_FORTRAN='no'
ENABLE_PARALLEL='no'
ENABLE_PYTHON='no'

###############################################################################
# --enable-xxx and --with-xxx Arguments
AC_ARG_ENABLE(
	[debug],
	[AS_HELP_STRING([--enable-debug], [Compile with debug flags [default=no]])],
	[ENABLE_DEBUG=$enableval])
AM_CONDITIONAL([ENABLE_DEBUG], [test "$ENABLE_DEBUG" = "yes"])

AC_ARG_ENABLE(
	[c],
	[AS_HELP_STRING([--enable-c],
	[Compile the C interface [default=yes]])],
	[ENABLE_C=$enableval])
AM_CONDITIONAL([ENABLE_C], [test "$ENABLE_C" = "yes"])

AC_ARG_ENABLE(
	[fortran],
	[AS_HELP_STRING([--enable-fortran],
	[Compile the Fortran interface [default=no]])],
	[ENABLE_FORTRAN=$enableval])
AM_CONDITIONAL([ENABLE_FORTRAN], [test "$ENABLE_FORTRAN" = "yes"])

AC_ARG_ENABLE(
	[python],
	[AS_HELP_STRING([--enable-python],
	[Compile the Python interface [default=no]])],
	[ENABLE_PYTHON=$enableval])
AM_CONDITIONAL([ENABLE_PYTHON], [test "$ENABLE_PYTHON" = "yes"])

AC_ARG_ENABLE(
	[parallel],
	[AS_HELP_STRING([--enable-parallel],
	[Compile the MPI/IO interface [default=no]])],
	[ENABLE_PARALLEL=$enableval])
AM_CONDITIONAL([ENABLE_PARALLEL], [test "$ENABLE_PARALLEL" = "yes"])

AC_ARG_ENABLE(
	[large-indices],
	[AC_HELP_STRING([--enable-large-indices],
	[Compile with 64bit local IDs and indices [default=no]])],
	[USE_LARGE_INDICES=$enableval])

AC_ARG_ENABLE(
        [experimental],
        [AS_HELP_STRING([--enable-experimental],
        [Compile experimental code [default=no]])],
        [ENABLE_EXPERIMENTAL=$enableval])
AM_CONDITIONAL([ENABLE_EXPERIMENTAL], [test "$ENABLE_EXPERIMENTAL" = "yes"])

AC_ARG_ENABLE(
	[vtkconverter],
	[AS_HELP_STRING([--enable-vtkconverter],
	[Compile the vtk to H5hut converter for grids [default=no]])],
	[ENABLE_VTKCONVERTER=$enableval])
AM_CONDITIONAL([ENABLE_VTKCONVERTER], [test "$ENABLE_VTKCONVERTER" = "yes"])

AC_ARG_ENABLE(
	[examples],
	[AS_HELP_STRING([--enable-examples],
	[Compile the H5hut examples [default=no]])],
	[ENABLE_EXAMPLES=$enableval])
AM_CONDITIONAL([ENABLE_EXAMPLES], [test "$ENABLE_EXAMPLES" = "yes"])

AC_ARG_WITH(
	[hdf5],
	[AS_HELP_STRING([--with-hdf5],
	[path to HDF5 installation [default=""]])],
	[HDF5_PREFIX=$withval], [])

AC_ARG_WITH(
	[mpi],
	[AS_HELP_STRING([--with-mpi],
	[path to MPI installation [default=""]])],
	[PATH=$MPI_PREFIX/bin:$PATH])

AC_ARG_WITH(
	[lustre],
	[AS_HELP_STRING([--with-lustre],
	[path to lustre user API [default=""]])],
	[LUSTREPATH=$withval], [LUSTREPATH=""])

AC_ARG_WITH(
	[parmetis],
	[AC_HELP_STRING([--with-parmetis],
	[path to ParMETIS API [default=""]])],
	[PARMETIS_PREFIX=$withval], [])

AC_ARG_WITH(
	[vtk],
	[AC_HELP_STRING([--with-vtk],
	[VTK installation path [default=""]])],
	[VTK_PREFIX=$withval], [])

# Determine a C/C++ compiler to use. 
if test "x$ENABLE_PARALLEL" = "xyes"; then
   	if test "x$MPICC" != "x"; then
	   	CC=$MPICC
	fi
   	if test "x$MPICXX" != "x"; then
	   	CXX=$MPICXX
	fi
   	if test "x$MPIFC" != "x"; then
	   	FC=$MPIFC
	else
		FC=''
		F77=''
		F90=''
	fi
   	c_compilers="mpicc cc"
	cxx_compilers="mpic++ CC"
	fc_compiler="mpif90 mpif77 mpifort"
else
	c_compilers="pgcc pathcc icc gcc cc_r cc"
	cxx_compilers="pgcc pathcc icc g++ cc_r c++"
	fc_compiler="pgf90 ifort xlf_r pathf90 g95 g90 ftn gfortran"
fi

AC_PROG_CC($c_compilers)
if test "x$ac_cv_prog_cc_c89" = "xno"; then
	AC_MSG_ERROR([C compiler not found!])
fi

AC_PROG_CXX($cxx_compilers)

if test "X$ENABLE_FORTRAN" = "Xyes"; then
   	AC_PROG_FC($fc_compiler)
	if test "x$FC" = "x"; then
	   	AC_MSG_ERROR([Fortran compiler not found!])
	fi
fi

CC=`which $CC`
CXX=`which $CXX`
FC=`which $FC`

AX_CXX_COMPILE_STDCXX_11
AC_PROG_CC_C99
if test "x$ac_cv_prog_cc_c99" = "xno"; then
	AC_MSG_ERROR([Cannot set C compiler to use C99 standard!])
fi


###############################################################################
############### PATH SERACH FUNCTION - to be used later... ####################
###############################################################################
# /*@@
#   @routine    CCTK_Search
#   @date       Wed Jul 21 11:16:35 1999
#   @author     Tom Goodale
#   @desc
#   Used to search for something in various directories
#   @enddesc
#@@*/

PATH_Search() {
	eval  $1=""
	if test $# -lt 4 ; then
		h5part_basedir=""
	else
		h5part_basedir="$4/"
	fi
	for h5part_place in $2; do
		AC_MSG_CHECKING([looking in $h5part_place ... ])
		if test -r "$h5part_basedir$h5part_place/$3" ; then
			AC_MSG_RESULT([found])
			eval $1="$h5part_place"
			break
		fi
		AC_MSG_RESULT([no])
	done
	return
}

###############################################################################
# MISC SETTINGS INCLUDING C & C++ COMPILER SETTING

AC_CANONICAL_HOST
AC_PROG_MAKE_SET
        
CFLAGS="$CFLAGS -Wall"
CXXFLAGS="$CFLAGS -Wall"
FFLAGS="$FFLAGS -Wall"

# debug enabled?
AC_MSG_CHECKING([if debug is enabled])
if test "X$ENABLE_DEBUG" = "Xyes"; then
	AC_MSG_RESULT([yes])
	CFLAGS="$CFLAGS -g3 -O0"
	CXXFLAGS="$CFLAGS -g3 -O0"
	FFLAGS="$FFLAGS -g3 -O0"
else
	AC_MSG_RESULT([no])
	CFLAGS="$CFLAGS -O3"
	CXXFLAGS="$CFLAGS -O3"
	FFLAGS="$FFLAGS -O3"
fi


AC_PROG_INSTALL
AC_PROG_AWK

# Disable shared libraries by default: can be enabled with --enable-shared
LT_INIT([disable-shared])
AC_PROG_LIBTOOL

# Default prefix for bindir, etc... (eg >> ./build/bin)
AC_PREFIX_DEFAULT(`pwd`/build)

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h limits.h stdint.h stdlib.h string.h sys/ioctl.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memmove memset pow strchr strdup strerror strstr])

# AC_DEFINE_UNQUOTED (variable, value, [description])
# Define the C preprocessor variable variable to value
# Use this macro instead of AC_DEFINE when variable or value is a shell variable. 
AC_DEFINE_UNQUOTED(MY_BUILD_VENDOR, "$host_vendor", "")
AC_DEFINE_UNQUOTED(MY_BUILD_CPU, "$host_cpu", "")
AC_DEFINE_UNQUOTED(MY_BUILD_OS, "$host_os", "")
AC_DEFINE_UNQUOTED(MY_GNUNAME, "${host_cpu}-${host_vendor}-${host_os}", "")
AC_DEFINE_UNQUOTED(MY_UNAME, "$uname", "")


###############################################################################
# PROCESS ARGUMENTS

###############################################################################
# check libraries

AC_CHECK_LIB([z], [compress2])

# hdf5
if test -z "${HDF5_PREFIX}"; then
	# alternatives to HDF5_PREFIX
	if test -n "${HDF5_DIR}"; then
		HDF5_PREFIX=${HDF5_DIR}
	fi
fi
if test -n "${HDF5_PREFIX}"; then
        CPPFLAGS="$CPPFLAGS -I$HDF5_PREFIX/include"
        LDFLAGS="$LDFLAGS -L$HDF5_PREFIX/lib"
fi
AC_CHECK_LIB([hdf5], [H5open])
AS_IF([test "x$ac_cv_lib_hdf5_H5open" != xyes], 
        [AC_MSG_ERROR([Cannot compile and link a HDF5 program])])
AC_MSG_CHECKING([checking whether hdf5 lib supports parallel IO])
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include "hdf5.h"
#ifndef H5_HAVE_PARALLEL
# error macro not defined
#endif
]])], [H5_HAVE_PARALLEL="yes"], [H5_HAVE_PARALLEL="no"])
AC_MSG_RESULT([$H5_HAVE_PARALLEL])

# lustre
AC_MSG_CHECKING([for lustre API ])
AC_MSG_RESULT([])
if test -n "${LUSTREPATH}"; then
	P=${LUSTREPATH}
elif test -n "${LUSTREROOT}" ; then
	P=${LUSTREROOT}
elif test -n "${LUSTREHOME}" ; then
	P=${LUSTREHOME}
elif test -n "${LUSTRE_DIR}" ; then
	P=${LUSTRE_DIR}
else
	P=''
	P="$P /usr"
	P="$P /usr/local"
	P="$P /usr/local/lustre"
	P="$P /opt/lustre"
fi
PATH_Search LUSTREROOT "$P" usr/include/lustre/liblustreapi.h
if test -z "$LUSTREROOT"; then
	AC_MSG_WARN([Couldn't locate the lustre API... building without support for lustre striping!])
else
        CFLAGS="$CFLAGS -DH5_USE_LUSTRE"
        CPPFLAGS="$CPPFLAGS -I$LUSTREROOT/usr/include"
        LDFLAGS="$LDFLAGS -L$LUSTREROOT/usr/lib"
        LIBS="$LIBS -llustreapi"
fi

LIBS="$LIBS -lm"



############################ enabled bindings ###############################

# C enabled?
AC_MSG_CHECKING([if C interface enabled])
if test "X$ENABLE_C" = "Xyes"; then
	AC_MSG_RESULT([yes])
else
        AC_MSG_RESULT([no])
fi

# fortran enabled?
AC_MSG_CHECKING([if fortran interface enabled])
if test "X$ENABLE_FORTRAN" = "Xyes"; then
	AC_MSG_RESULT([yes])
	AC_LANG_PUSH(Fortran)
	
	# Query Fortran symbol mangling. __AC_FC_NAME_MANGLING sets
	# ac_cv_{f77,fc}_mangling. The value contains three fields, separated
	# by commas:
	#
	# lower case / upper case:
	#    case translation of the Fortran symbols
	# underscore / no underscore:
	#    whether the compiler appends "_" to symbol names
	# extra underscore / no extra underscore:
	#    whether the compiler appends an extra "_" to symbol names already
	#    containing at least one underscore
	__AC_FC_NAME_MANGLING
	if test "X${ac_cv_fc_mangling}" == "Xlower case, no underscore, no extra underscore"; then
	    	CPPFLAGS+=" -DFC_MANGLING_LOWERCASE_NO_UNDERSCORE"
	elif test "X${ac_cv_fc_mangling}" == "Xlower case, underscore, no extra underscore"; then
	    	CPPFLAGS+=" -DFC_MANGLING_LOWERCASE_SINGLE_UNDERSCORE"
	elif test "X${ac_cv_fc_mangling}" == "Xlower case, underscore, extra underscore"; then
	    	CPPFLAGS+=" -DFC_MANGLING_LOWERCASE_DOUBLE_UNDERSCORE"
	elif test "X${ac_cv_fc_mangling}" == "Xupper case, no underscore, no extra underscore"; then
	    	CPPFLAGS+=" -DFC_MANGLING_UPPERCASE_NO_UNDERSCORE"
	elif test "X${ac_cv_fc_mangling}" == "Xupper case, underscore, no extra underscore"; then
	    	CPPFLAGS+=" -DFC_MANGLING_UPPERCASE_SINGLE_UNDERSCORE"
	elif test "X${ac_cv_fc_mangling}" == "Xupper case, underscore, extra underscore"; then
	    	CPPFLAGS+=" -DFC_MANGLING_UPPERCASE_DOUBLE_UNDERSCORE"
	fi
	AC_LANG_POP(Fortran)
else
	AC_MSG_RESULT([no])
fi

# python enabled?
if test "X$ENABLE_PYTHON" = "Xyes"; then

    	AM_PATH_PYTHON([2.5])
    	AX_PKG_SWIG([1.3], [], [ AC_MSG_ERROR([SWIG is required to build..]) ])
    	#SWIG_ENABLE_CXX
    	#SWIG_PYTHON

    	PYTHON_PREFIX=`python-config --prefix`
    	PYTHON_EXEC_PREFIX=`python-config --exec-prefix`

    	# make sure numpy is available..
    	AC_MSG_CHECKING([for numpy support])
    	NUMPY_INCLUDE=`${PYTHON} -c 'import numpy; print numpy.get_include()'`

	if test "X$?" = "X0" ; then
		AC_MSG_RESULT([$NUMPY_INCLUDE])
    	else
		AC_MSG_RESULT([no])
 		exit 1
    	fi
   	AC_SUBST([NUMPY_INCLUDE])

    	#check for parallel support
    	if test "X$ENABLE_PARALLEL" = "Xyes"; then
        	AC_MSG_CHECKING([for mpi4py support])
        	MPI4PY_INCLUDE=`${PYTHON} -c 'import mpi4py; print mpi4py.get_include()'`

        	if test "X$?" = "X0" ; then
			AC_MSG_RESULT([$MPI4PY_INCLUDE])
        	else
            		AC_MSG_RESULT([no])
            		exit 1
        	fi
		AC_SUBST([MPI4PY_INCLUDE])
    	fi

	AC_MSG_RESULT([ok])

else
    	AC_MSG_CHECKING([if python interface is enabled])
    	AC_MSG_RESULT([no])
fi

######################## large indices enabled ###########################
AC_MSG_CHECKING([if large indices are enabled])
if test "X$USE_LARGE_INDICES" = "Xyes"; then
	AC_MSG_RESULT([yes])
	CFLAGS="${CFLAGS} -DUSE_LARGE_INDICES"
else
	AC_MSG_RESULT([no])
fi

######################## parallel interface enabled ###########################
AC_MSG_CHECKING([if parallel interface enabled])
if test "X$ENABLE_PARALLEL" = "Xyes"; then
	AC_MSG_RESULT([yes])

	AS_IF([test "x$H5_HAVE_PARALLEL" = "xno"],
	        [AC_MSG_ERROR([hdf5 does not support parallel IO!])])

	if test "X$ENABLE_VTKCONVERTER" = "Xyes"; then
	   	AC_MSG_ERROR([VTK to H5hut grid converter cannot be compiled parallel!])
	fi
	CPPFLAGS="${CPPFLAGS} -DPARALLEL_IO -DMPICH_IGNORE_CXX_SEEK"

        if test "X$ENABLE_EXPERIMENTAL" = "Xyes"; then
		CPPFLAGS="$CPPFLAGS -DWITH_PARALLEL_H5GRID"
        	AC_MSG_CHECKING([for ParMETIS])
                if test "X$PARMETIS_PREFIX" != "X"; then
                        CPPFLAGS="$CPPFLAGS -I$PARMETIS_PREFIX/include"
	                LDFLAGS="$LDFLAGS -L$PARMETIS_PREFIX/lib"
	        fi

                LIBS="$LIBS -lparmetis -lmetis"
	        AC_TRY_LINK([#include <parmetis.h>], [
	                ParMETIS_V3_PartGeom (NULL, NULL, NULL, NULL, NULL);],
		        [AC_MSG_RESULT([yes]); r='yes'], [AC_MSG_RESULT([no]); r='no'] )

	        if test "X$r" = "Xno"; then
		        AC_MSG_ERROR([ParMETIS])
	        fi
        fi
else  # --enable-parallel=no
	AC_MSG_RESULT([no])
	AS_IF([test "x$H5_HAVE_PARALLEL" = "xyes"],
	        [AC_MSG_ERROR([hdf5 with support for parallel IO cannot be used for non-parallel H5hut!])])
fi

###############################################################################
# TOOLS
AC_MSG_CHECKING([if we have to compile the VTK to H5hut grid converter])
if test "X$ENABLE_VTKCONVERTER" = "Xyes"; then
	AC_MSG_RESULT([yes])
	if test -n "$VTK_PREFIX"; then
		SAVE_C_INCLUDE_PATH="$C_INCLUDE_PATH"
		SAVE_CPLUS_INCLUDE_PATH="$CPLUS_INCLUDE_PATH"
		SAVE_LIBRARY_PATH="$LIBRARY_PATH"
		C_INCLUDE_PATH="$C_INCLUDE_PATH $VTK_PREFIX/include"
		CPLUS_INCLUDE_PATH="$CPLUS_INCLUDE_PATH $VTK_PREFIX/include"
		LIBRARY_PATH="$LIBRARY_PATH $VTK_PREFIX/lib"
	fi
	AC_MSG_CHECKING([VTK inlcude files])
	VTK_INCLUDES=`cmake --find-package -DNAME=VTK -DCOMPILER_ID=GNU -DLANGUAGE=C -DMODE=COMPILE 2>/dev/null`
	if test "X$?" != "X0"; then
		AC_MSG_RESULT([no])
	        AC_MSG_ERROR([Cannot compile VTK to H5hut converter])
	else
		AC_MSG_RESULT([yes])
	fi
	AC_MSG_CHECKING([VTK libraries])
	VTK_LIBS=`cmake --find-package -DNAME=VTK -DCOMPILER_ID=GNU -DLANGUAGE=C -DMODE=LINK 2>/dev/null`
	if test "X$?" != "X0"; then
		AC_MSG_RESULT([no])
	        AC_MSG_ERROR([Cannot compile VTK to H5hut converter]); 
	else
		AC_MSG_RESULT([yes])
	fi

else
	AC_MSG_RESULT([no])
fi

###############################################################################
# EXAMPLES
AC_MSG_CHECKING([if we have to compile the H5hut examples])
if test "X$ENABLE_EXAMPLES" = "Xyes"; then
	AC_MSG_RESULT([yes])
else
	AC_MSG_RESULT([no])
fi

###############################################################################
# EXPORTING VARIABLES & CREATING OUTPUT FILES
AM_CPPFLAGS=$CPPFLAGS
AM_LDFLAGS=$LDFLAGS

AC_SUBST(HDF5_PREFIX)
AC_SUBST(CFLAGS)
AC_SUBST(FFLAGS)
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_LDFLAGS)
AC_SUBST(VTK_INCLUDES)
AC_SUBST(VTK_LIBS)
AC_CONFIG_FILES([
	Makefile
	src/Makefile
	src/h5core/Makefile
	src/include/Makefile
	src/include/h5core/Makefile
        src/C/Makefile
        src/Fortran/Makefile
        src/Python/Makefile
	test/Makefile
        doc/Makefile
	tools/Makefile
	tools/h5hutcc
        tools/vtk2h5grid/Makefile
        examples/H5/Makefile
        examples/H5Part/Makefile
        examples/H5Block/Makefile
        examples/H5Fed/Makefile
])

AC_OUTPUT

###############################################################################
# PRINTING SUMMARY
AC_MSG_RESULT([ ])
AC_MSG_RESULT([Summary:])
AC_MSG_RESULT([ ])
AC_MSG_RESULT([Host OS:             $host_os])
AC_MSG_RESULT([Host CPU:            $host_cpu])
AC_MSG_RESULT([Host vendor:         $host_vendor])
AC_MSG_RESULT([CC =                 $CC])
AC_MSG_RESULT([CXX =                $CXX])
AC_MSG_RESULT([FC =                 $FC])
AC_MSG_RESULT([CFLAGS =             $CFLAGS])
AC_MSG_RESULT([FFLAGS =             $FFLAGS])
AC_MSG_RESULT([AM_CPPFLAGS =        $AM_CPPFLAGS])
AC_MSG_RESULT([AM_LDFLAGS =         $AM_LDFLAGS])
AC_MSG_RESULT([LIBS =               $LIBS])
AC_MSG_RESULT([HDF5_PREFIX =        $HDF5_PREFIX])
AC_MSG_RESULT([LUSTREROOT =         $LUSTREROOT])

if test "x$ENABLE_PYTHON" = "xyes"; then
	AC_MSG_RESULT([PYTHON =             $PYTHON])
	AC_MSG_RESULT([SWIG =               $SWIG])
	AC_MSG_RESULT([NUMPY_INCLUDE =      $NUMPY_INCLUDE])
	if test "x$ENABLE_PARALLEL" = "xyes"; then
		AC_MSG_RESULT([MPI4PY_INCLUDE =     $MPI4PY_INCLUDE])
	fi
fi
AC_MSG_RESULT([ ])
