AM_CPPFLAGS = -I${abs_top_srcdir}/src/include -Wno-deprecated @VTK_INCLUDES@
FFLAGS += -cpp $(AM_CPPFLAGS)
AM_LDFLAGS += -L${abs_top_builddir}/src/h5core

LDADD = -lH5hut 
LDADD += @VTK_LIBS@

bin_PROGRAMS = 

bin_PROGRAMS += vtk2h5grid
vtk2h5grid_SOURCES = vtk2h5grid.cc

distclean-local:
	$(RM) *~
