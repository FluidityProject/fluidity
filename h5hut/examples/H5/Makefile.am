#
# Copyright (c) 2006-2015, The Regents of the University of California,
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

if ENABLE_FORTRAN
LDADD += -lH5hutF
AM_LDFLAGS += -L${abs_top_builddir}/src/Fortran/.libs
endif

if ENABLE_DEBUG
AM_LDFLAGS += -static
endif

LDADD += -lH5hut

noinst_PROGRAMS = 

if ENABLE_EXAMPLES
if ENABLE_C
noinst_PROGRAMS +=		\
	attach_file		\
	openclose		\
	query			\
	read_file_attribs	\
	read_step_attribs	\
	write_file_attribs	\
	write_step_attribs
endif

if ENABLE_FORTRAN
noinst_PROGRAMS +=		\
	openclosef		\
	queryf			\
	read_file_attribsf	\
	read_step_attribsf	\
	write_file_attribsf 	\
	write_step_attribsf

openclosef_SOURCES = openclosef.f90
queryf_SOURCES = queryf.f90
read_file_attribsf_SOURCES = read_file_attribsf.f90
read_step_attribsf_SOURCES = read_step_attribsf.f90
write_file_attribsf_SOURCES = write_file_attribsf.f90
write_step_attribsf_SOURCES = write_step_attribsf.f90
endif
endif

all: ${noinst_PROGRAMS}

.f90.o:
	$(FC)  $(FFLAGS) -c $<
