##
## Programs and options...
##
TOPDIR=/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/
include $(TOPDIR)/include/make-variables
TOPDIR=/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/
HOMEPLUGINS=/d/visusers/cristina/.visit/linux-x86_64/plugins

##
## Libraries and includes...
##
VTK_INCLUDE= \
  -I$(TOPDIR)/include/vtk \
  -I$(TOPDIR)/include/vtk/Common \
  -I$(TOPDIR)/include/vtk/Filtering \
  -I$(TOPDIR)/include/vtk/Graphics \
  -I$(TOPDIR)/include/vtk/Hybrid \
  -I$(TOPDIR)/include/vtk/IO \
  -I$(TOPDIR)/include/vtk/Imaging \
  -I$(TOPDIR)/include/vtk/Rendering
MOC=$(TOPDIR)/bin/moc
CXXFLAGS=$(CXXFLAGSORIG) -I/usr/local/hdf5/include -I/d/visusers/cristina/H5Part/H5ecloud 
CPPFLAGS=$(CPPFLAGSORIG) $(VTK_INCLUDE) -I. -I/usr/local/hdf5/include -I/d/visusers/cristina/H5Part/H5ecloud -I$(TOPDIR)/include -I$(TOPDIR)/include/visit
LDFLAGS=$(LDFLAGSORIG)  -L$(HOMEPLUGINS)/databases

##
## Files...
##
PLUGINDIR=databases
PLUGINNAME=h5partDatabase

WIDGETS=
ISRC=h5partPluginInfo.C
COMMONSRC=h5partPluginInfo.C h5partCommonPluginInfo.C
MSRC=h5partMDServerPluginInfo.C avth5partFileFormat.C
ESRC=h5partEnginePluginInfo.C avth5partFileFormat.C


SRC=$(ISRC) $(COMMONSRC) $(MSRC) $(ESRC)

SHLIB_FORCED=-L/d/visusers/cristina/H5Part/H5ecloud/ -lH5Part 

ELIBS_FOR_MACOSX_PREBINDING=$(BZIP2_LIBS)
MLIBS_FOR_MACOSX_PREBINDING=$(BZIP2_LIBS)
ILIBS=
MLIBS=-lpipeline_ser -ldatabase_ser -lmir_ser -lavtmath_ser -ldbatts -lavtexceptions -lstate -lcomm -lmisc -lexpr -lparser -lplugin -lutility -lvisit_vtk -llightweight_visit_vtk  -L$(HOMEPLUGINS)/databases -lI$(PLUGINNAME) $(MLIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS)
ESERLIBS=-lpipeline_ser -ldatabase_ser -lmir_ser -lavtmath_ser -ldbatts -lavtexceptions -lstate -lcomm -lmisc -lexpr -lparser -lplugin -lutility -lvisit_vtk -llightweight_visit_vtk  -L$(HOMEPLUGINS)/databases -lI$(PLUGINNAME) $(ELIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS)
EPARLIBS=-lpipeline_par -ldatabase_par -lmir_par -lavtmath_par -ldbatts -lavtexceptions -lstate -lcomm -lmisc -lexpr -lparser -lplugin -lutility -lvisit_vtk -llightweight_visit_vtk  -L$(HOMEPLUGINS)/databases -lI$(PLUGINNAME) $(ELIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS) $(SHLIB_MPI_LIBS)

IDSO=$(HOMEPLUGINS)/databases/libIh5partDatabase.so
MDSO=$(HOMEPLUGINS)/databases/libMh5partDatabase.so
ESERDSO=$(HOMEPLUGINS)/databases/libEh5partDatabase_ser.so
EPARDSO=$(HOMEPLUGINS)/databases/libEh5partDatabase_par.so


DISTRIB=

##
## Derived objects
##
IOBJ=$(ISRC:.C=.o)
GOBJ=$(COMMONSRC:.C=.o) $(GSRC:.C=.o)
SOBJ=$(COMMONSRC:.C=.o) $(SSRC:.C=.o)
VOBJ=$(COMMONSRC:.C=.o) $(VSRC:.C=.o)
MOBJ=$(COMMONSRC:.C=.o) $(MSRC:.C=.o) $(MSPECIFICSRC:.C=_mds.o)
ESEROBJ=$(COMMONSRC:.C=.o) $(ESRC:.C=.o)
EPAROBJ=$(COMMONSRC:.C=.o) $(ESRC:.C=_par.o)

MOCSRC = $(WIDGETS:.h=_moc.C)
MOCOBJ = $(MOCSRC:.C=.o)

##
## Standard targets...
##
all: message $(IDSO) $(MDSERVERLIB) $(ENGINELIBSER) $(ENGINELIBPAR)

clean:
	$(RM) $(IOBJ) $(COMMONOBJ)
	$(RM) $(GOBJ) $(SOBJ) $(VOBJ) $(MOBJ) $(ESEROBJ) $(EPAROBJ)
	$(RM) $(MOCSRC) $(MOCOBJ)
	$(RM) $(IDSO) $(GDSO) $(SDSO) $(VDSO) $(MDSO) $(ESERDSO) $(EPARDSO)
	$(RM) $(JAVAOBJ)

##
## Other targets...
##
message:
	@echo
	@echo "****************************************************************************"
	@echo "*** Building h5part Database Plugin"
	@echo "****************************************************************************"

##
## moc
##
$(MOCSRC) or_no_widgets: $(WIDGETS)
	@rm -f $@
	$(MOC) $(@:_moc.C=.h) > $@


##
## Automatic dependency stuff
##
include $(TOPDIR)/include/make-targets
