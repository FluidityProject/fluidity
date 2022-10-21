#
# Copyright (c) 2006-2017, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any
# required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
# Institut (Switzerland).  All rights reserved.!
#
# License: see file COPYING in top level of source distribution.
#

AM_CPPFLAGS += -I${abs_top_srcdir}/src/include
AM_CPPFLAGS += -I${abs_top_srcdir}/examples/include
FFLAGS += -cpp $(AM_CPPFLAGS)
AM_LDFLAGS += -L${abs_top_builddir}/src/h5core/.libs

LDADD = 

if ENABLE_EXAMPLES
if ENABLE_C
LDADD +=
endif

if ENABLE_FORTRAN
LDADD += -lH5hutF
AM_LDFLAGS += -L${abs_top_builddir}/src/Fortran/.libs
endif

LDADD += -lH5hut

noinst_PROGRAMS = 

if ENABLE_C
noinst_PROGRAMS +=		\
	query			\
	read_core_vfd		\
	read_canonicalview	\
	read_setnparticles	\
	read_setview		\
	read_strided		\
	write_core_vfd		\
	write_setnparticles	\
	write_setview		\
	write_strided
endif

if ENABLE_FORTRAN
noinst_PROGRAMS +=		\
	read_core_vfdf		\
	read_canonicalviewf	\
	read_setnparticlesf	\
	read_setviewf		\
	read_stridedf		\
	write_core_vfdf		\
	write_setnparticlesf	\
	write_setviewf		\
	write_stridedf
endif

EXTRA_PROGRAMS = 

read_core_vfdf_SOURCES = read_core_vfdf.f90
read_canonicalviewf_SOURCES = read_canonicalviewf.f90
read_setnparticlesf_SOURCES = read_setnparticlesf.f90
read_setviewf_SOURCES = read_setviewf.f90
read_stridedf_SOURCES = read_stridedf.f90
write_core_vfdf_SOURCES = write_core_vfdf.f90
write_setnparticlesf_SOURCES = write_setnparticlesf.f90
write_setviewf_SOURCES = write_setviewf.f90
write_stridedf_SOURCES = write_stridedf.f90

endif

%.o : %.f90
	$(FC)  $(FFLAGS) -c $<

clean-local:
	$(RM) -f *~


