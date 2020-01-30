#
# Copyright (c) 2006-2016, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any
# required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
# Institut (Switzerland).  All rights reserved.!
#
# License: see file COPYING in top level of source distribution.
#

AM_CPPFLAGS += -I${abs_top_srcdir}/src/include
FFLAGS += -cpp $(AM_CPPFLAGS)
AM_LDFLAGS += -L${abs_top_builddir}/src/h5core/.libs
if ENABLE_DEBUG
AM_LDFLAGS += -static
endif

LDADD =
LDADD += -lH5hut

if ENABLE_EXAMPLES
noinst_PROGRAMS =				\
		tetmesh_write			\
		tetmesh_write1			\
		tetmesh_write2			\
		tetmesh_read			\
		tetmesh_adjacencies		\
		tetmesh_write_tags		\
		tetmesh_read_tags		\
		trimesh_read			\
		trimesh_write			\
		trimesh_2gnuplot		\
		trimesh_adjacencies		\
		trimesh_write_dunetest

#		map_tet2globalid \
#		map_triangle2globalid

# Extra files that I wish to include in the dist tar ball.
EXTRA_DIST =					\
		tetmesh_adjacencies.c		\
		tetmesh_read.c			\
		tetmesh_read_tags.c		\
		tetmesh_write.c			\
		tetmesh_write1.c		\
		tetmesh_write2.c		\
		tetmesh_write_tags.c		\
		trimesh_adjacencies.c		\
		trimesh_read.c			\
		trimesh_write.c			\
		trimesh_2gnuplot.c		\
		trimesh_write_dunetest.c
#		map_tet2globalid.c \
#		map_triangle2globalid.c \
#		$(bin_SCRIPTS)

endif

clean-local:
	$(RM) *.h5 *~

