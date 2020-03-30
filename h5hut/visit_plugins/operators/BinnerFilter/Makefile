##
## Programs and options...
##
TOPDIR=/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64
include $(TOPDIR)/include/make-variables
TOPDIR=/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64
PLUGINSHOME=$(HOME)/.visit/linux-x86_64/plugins

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
CXXFLAGS=$(CXXFLAGSORIG) $(QT_CXXFLAGS) $(PY_CXXFLAGS)
CPPFLAGS=$(CPPFLAGSORIG) $(VTK_INCLUDE) -I. -I$(TOPDIR)/include -I$(TOPDIR)/include/visit
LDFLAGS=$(LDFLAGSORIG) $(PY_LDFLAGS)  -L$(PLUGINSHOME)/operators

##
## Files...
##
PLUGINDIR=operators
PLUGINNAME=BinnerOperator

WIDGETS=QvisBinnerWindow.h
ISRC=BinnerPluginInfo.C
COMMONSRC=BinnerPluginInfo.C BinnerCommonPluginInfo.C BinnerFilter.C
GSRC=BinnerGUIPluginInfo.C QvisBinnerWindow.C
VSRC=BinnerViewerPluginInfo.C
ESRC=BinnerEnginePluginInfo.C avtBinnerFilter.C
SSRC=BinnerScriptingPluginInfo.C PyBinner.C

SRC=$(ISRC) $(COMMONSRC) $(GSRC) $(VSRC) $(ESRC) $(SSRC)


JAVASRC=Binner.java
JAVAOBJ=$(TOPDIR)/java/llnl/visit/operators/BinnerFilter.class
JAVAPLUGINFLAGS=-d $(TOPDIR)/java -classpath $(TOPDIR)/java

SHLIB_FORCED=

ELIBS_FOR_MACOSX_PREBINDING=$(BZIP2_LIBS)
VLIBS_FOR_MACOSX_PREBINDING=$(BZIP2_LIBS)
ILIBS=
GLIBS=-lgui -lmdserverproxy -lviewerproxy -lproxybase -lmdserverrpc -lviewerrpc -lwinutil -ldbatts -lavtexceptions -lstate -lcomm -lmisc -lplugin -lexpr -lparser -lutility -lI$(PLUGINNAME) $(QT_LDFLAGS) $(QT_LIBS) $(QUI_LIBS) $(X_LIBS)
SLIBS=-lstate -lmisc -lcomm -lutility $(PY_LIB) -lI$(PLUGINNAME)
VLIBS=-lpipeline_ser -lplotter_ser -lavtfilters_ser -lavtmath_ser -lavtview -ldbatts -lavtexceptions -lstate -lmisc -lcomm -lexpr -lparser -lutility -lvisit_vtk -llightweight_visit_vtk -lparallel_visit_vtk_ser $(VLIBS_FOR_MACOSX_PREBINDING) -lI$(PLUGINNAME) $(VTK_LIBS)
ESERLIBS=-lpipeline_ser -lplotter_ser -lavtfilters_ser -lavtmath_ser -lavtview -ldbatts -lavtexceptions -lstate -lmisc -lcomm -lexpr -lparser -lutility -lvisit_vtk -llightweight_visit_vtk -lparallel_visit_vtk_ser -lI$(PLUGINNAME) $(ELIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS)
EPARLIBS=-lpipeline_par -lplotter_par -lavtfilters_par -lavtmath_par -lavtview -ldbatts -lavtexceptions -lstate -lmisc -lcomm -lexpr -lparser -lutility -lvisit_vtk -llightweight_visit_vtk -lparallel_visit_vtk_par -lI$(PLUGINNAME) $(ELIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS) $(SHLIB_MPI_LIBS)

IDSO=$(PLUGINSHOME)/operators/libIBinnerOperator.so
GDSO=$(PLUGINSHOME)/operators/libGBinnerOperator.so
SDSO=$(PLUGINSHOME)/operators/libSBinnerOperator.so
VDSO=$(PLUGINSHOME)/operators/libVBinnerOperator.so
ESERDSO=$(PLUGINSHOME)/operators/libEBinnerOperator_ser.so
EPARDSO=$(PLUGINSHOME)/operators/libEBinnerOperator_par.so


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
all: message $(IDSO) $(GUILIB) $(VIEWERLIB) $(ENGINELIBSER) $(ENGINELIBPAR) $(SCRIPTINGLIB) $(JAVACLASS)

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
	@echo "*** Building Binner Operator Plugin"
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
