#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#m4_pattern_allow([AM_CPPFLAGS])
#m4_pattern_allow([AM_LDFLAGS])
#m4_pattern_allow([AM_INIT_AUTOMAKE])
AC_PREREQ(2.59)
AC_INIT([H5PartMerge],[0.1],[christof.kraus@psi.ch])
AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR([./src/H5merge.cpp])

# Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_INSTALL

# Checks for libraries.
  Z_PREFIX=""
  Z_INCL_PATH=""
  Z_LIB_PATH=""
  AC_ARG_WITH(
    [zlib],
    AC_HELP_STRING(
      [--with-zlib=PREFIX],
      [prefix, parent directory, where the z library and header files can be found]
    ),
    [Z_PREFIX=$withval
     Z_INCL_PATH=${Z_PREFIX}/include
     Z_LIB_PATH=${Z_PREFIX}/lib]
  )
  if test "${Z_PREFIX}" != "" ; then
    CPPFLAGS=${CPPFLAGS}" -I"${Z_INCL_PATH}
    LDFLAGS=${LDFLAGS}" -L"${Z_LIB_PATH}
  fi
  AC_CHECK_LIB(
    [z],
    [main],
    [],
    [AC_MSG_ERROR(z library not found !)]
  )

AC_LANG_PUSH(C)
  HDF5_PREFIX=""
  HDF5_INCL_PATH=""
  HDF5_LIB_PATH=""
  AC_ARG_WITH(
    [hdf5],
    AC_HELP_STRING([--with-hdf5=PREFIX],
      [prefix, parent directory, where the serial(!) hdf5 library and header files can be found]),
    [HDF5_PREFIX=$withval
     HDF5_INCL_PATH=${HDF5_PREFIX}"/include "
     HDF5_LIB_PATH=${HDF5_PREFIX}"/lib "]
  )
  if test "${HDF5_PREFIX}" != "" ; then
    CPPFLAGS=${CPPFLAGS}" -I"${HDF5_INCL_PATH};
    LDFLAGS=${LDFLAGS}" -L"${HDF5_LIB_PATH};
  fi
  AC_CHECK_HEADERS(
    [hdf5.h],
    [],
    [AC_MSG_ERROR(hdf5.h header file not found!)]
  )
  AC_CHECK_LIB(
    [hdf5],
    [main],
    [],
    [AC_MSG_ERROR(hdf5 library not found!)]
  )

  H5HUT_PREFIX=""
  H5HUT_INCL_PATH=""
  H5HUT_LIB_PATH=""
  AC_ARG_WITH(
    [h5hut],
    AC_HELP_STRING([--with-h5hut=PREFIX],
      [prefix, parent directory, where the serial(!) h5hut library and header files can be found]),
    [H5HUT_PREFIX=$withval
     H5HUT_INCL_PATH=${H5HUT_PREFIX}"/include "
     H5HUT_LIB_PATH=${H5HUT_PREFIX}"/lib "]
  )
  if test "${H5HUT_PREFIX}" != ""; then
    CPPFLAGS=${CPPFLAGS}" -I"${H5HUT_INCL_PATH};
    LDFLAGS=${LDFLAGS}" -L"${H5HUT_LIB_PATH};
  fi
  AC_CHECK_HEADERS(
    [H5hut.h],
    [],
    [AC_MSG_ERROR(H5hut.h not found!)]
  )
  AC_CHECK_LIB(
    [H5hut],
    [main],
    [],
    [AC_MSG_ERROR(H5hut library not found!)]
  )
  AC_CHECK_LIB(
    [H5hutC],
    [main],
    [],
    [AC_MSG_ERROR(H5hut library not found!)]
  )


AC_LANG_POP
AC_LANG(C++)

  USE_BOOST_=false
  BOOST_INCL_PATH=""
  AC_ARG_ENABLE(
    [boost],
    AC_HELP_STRING([--enable-boost],
      [[enable usage of boost::any]]),
    [USE_BOOST_=true] )

  if test "x${USE_BOOST_}" = xtrue; then
    AC_ARG_WITH(
      [boost],
      AC_HELP_STRING([--with-boost-include=DIR],
        [include path for the boost library]),
      [BOOST_INCL_PATH=$withval]
    )

    if test "${BOOST_INCL_PATH}" != ""; then
      TEMP_CPPFLAGS=${CPPFLAGS}
      CPPFLAGS=${CPPFLAGS}" -I"${BOOST_INCL_PATH}
    fi

    AC_CHECK_HEADERS(
      [boost/any.hpp],
      [AC_DEFINE_UNQUOTED(USE_BOOST)],
      [USE_BOOST_=false
       AC_MSG_ERROR(any.hpp not found!)]
    )
  fi

  if test "x${USE_BOOST_}" != xtrue ; then
    if test "${BOOST_INCL_PATH}" != "" ; then
      CPPFLAGS=${TEMP_CPPFLAGS}
    fi
  fi

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
#AC_FUNC_MALLOC

CPPFLAGS=$CPPFLAGS" -m64 -mtune=k8"
LDFLAGS=$LDFLAGS" -lm"

AC_SUBST([AM_CPPFLAGS], $CPPFLAGS)
AC_SUBST([AM_LDFLAGS], $LDFLAGS)


AC_CONFIG_FILES([makefile])
AC_OUTPUT
