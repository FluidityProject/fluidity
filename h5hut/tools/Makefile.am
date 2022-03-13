SUBDIRS =

if ENABLE_VTKCONVERTER
SUBDIRS +=	vtk2h5grid
endif

AM_CPPFLAGS += -I${abs_top_builddir}/src/include -Wno-deprecated
FFLAGS += -cpp $(AM_CPPFLAGS)
AM_LDFLAGS += -L${abs_top_builddir}/src/lib

LDADD =
LDADD += -lH5hut

if ENABLE_C
EXTRA_SCRIPTS = #h5hutcc
EXTRA_PROGRAMS = h5pAttrib
endif


#EXTRA_DIST = h5pAttrib.cc h5pToGNUplot.cc## TO BE TAILORED LATER...

h5pAttrib_SOURCES = h5pAttrib.cc
#h5pToGNUplot_SOURCES = h5pToGNUplot.cc

clean: clean-am

