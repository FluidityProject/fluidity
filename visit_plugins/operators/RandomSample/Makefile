##
## Programs and options...
##
include /usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/include/make-variables
TOPDIR=/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64


##
## Libraries and includes...
##
VTK_INCLUDE= \
  -I/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/include/vtk \
  -I/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/include/vtk/Common \
  -I/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/include/vtk/Filtering \
  -I/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/include/vtk/Graphics \
  -I/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/include/vtk/Hybrid \
  -I/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/include/vtk/IO \
  -I/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/include/vtk/Imaging \
  -I/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/include/vtk/Rendering
MOC=/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/bin/moc
CXXFLAGS=$(CXXFLAGSORIG) $(QT_CXXFLAGS) $(PY_CXXFLAGS)
CPPFLAGS=$(CPPFLAGSORIG) $(VTK_INCLUDE) -I. -I/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/include -I/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/include/visit
LDFLAGS=$(LDFLAGSORIG) $(PY_LDFLAGS)  -L/d/visusers/cristina/.visit/linux-x86_64/plugins/operators

##
## Files...
##
PLUGINDIR=operators
PLUGINNAME=H5PartRandomSampleOperator

WIDGETS=QvisH5PartRandomSampleWindow.h
ISRC=H5PartRandomSamplePluginInfo.C
COMMONSRC=H5PartRandomSamplePluginInfo.C H5PartRandomSampleCommonPluginInfo.C H5PartRandomSampleAttributes.C
GSRC=H5PartRandomSampleGUIPluginInfo.C QvisH5PartRandomSampleWindow.C
VSRC=H5PartRandomSampleViewerPluginInfo.C
ESRC=H5PartRandomSampleEnginePluginInfo.C avtH5PartRandomSampleFilter.C
SSRC=H5PartRandomSampleScriptingPluginInfo.C PyH5PartRandomSampleAttributes.C

SRC=$(ISRC) $(COMMONSRC) $(GSRC) $(VSRC) $(ESRC) $(SSRC)


JAVASRC=H5PartRandomSampleAttributes.java
JAVAOBJ=/usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/java/llnl/visit/operators/H5PartRandomSampleAttributes.class
JAVAPLUGINFLAGS=-d /usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/java -classpath /usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/java

SHLIB_FORCED=

ELIBS_FOR_MACOSX_PREBINDING=$(BZIP2_LIBS)
VLIBS_FOR_MACOSX_PREBINDING=$(BZIP2_LIBS)
ILIBS=
GLIBS=-lgui -lmdserverproxy -lviewerproxy -lproxybase -lmdserverrpc -lviewerrpc -lwinutil -ldbatts -lavtexceptions -lstate -lcomm -lmisc -lplugin -lexpr -lparser -lutility -lI$(PLUGINNAME) $(QT_LDFLAGS) $(QT_LIBS) $(QUI_LIBS) $(X_LIBS)
SLIBS=-lstate -lmisc -lcomm -lutility $(PY_LIB) -lI$(PLUGINNAME)
VLIBS=-lpipeline_ser -lplotter_ser -lavtfilters_ser -lavtmath_ser -lavtview -ldbatts -lavtexceptions -lstate -lmisc -lcomm -lexpr -lparser -lutility -lvisit_vtk -llightweight_visit_vtk -lparallel_visit_vtk_ser $(VLIBS_FOR_MACOSX_PREBINDING) -lI$(PLUGINNAME) $(VTK_LIBS)
ESERLIBS=-lpipeline_ser -lplotter_ser -lavtfilters_ser -lavtmath_ser -lavtview -ldbatts -lavtexceptions -lstate -lmisc -lcomm -lexpr -lparser -lutility -lvisit_vtk -llightweight_visit_vtk -lparallel_visit_vtk_ser -lI$(PLUGINNAME) $(ELIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS)
EPARLIBS=-lpipeline_par -lplotter_par -lavtfilters_par -lavtmath_par -lavtview -ldbatts -lavtexceptions -lstate -lmisc -lcomm -lexpr -lparser -lutility -lvisit_vtk -llightweight_visit_vtk -lparallel_visit_vtk_par -lI$(PLUGINNAME) $(ELIBS_FOR_MACOSX_PREBINDING) $(VTK_LIBS) $(SHLIB_MPI_LIBS)

IDSO=/d/visusers/cristina/.visit/linux-x86_64/plugins/operators/libIH5PartRandomSampleOperator.so
GDSO=/d/visusers/cristina/.visit/linux-x86_64/plugins/operators/libGH5PartRandomSampleOperator.so
SDSO=/d/visusers/cristina/.visit/linux-x86_64/plugins/operators/libSH5PartRandomSampleOperator.so
VDSO=/d/visusers/cristina/.visit/linux-x86_64/plugins/operators/libVH5PartRandomSampleOperator.so
ESERDSO=/d/visusers/cristina/.visit/linux-x86_64/plugins/operators/libEH5PartRandomSampleOperator_ser.so
EPARDSO=/d/visusers/cristina/.visit/linux-x86_64/plugins/operators/libEH5PartRandomSampleOperator_par.so


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
	@echo "*** Building H5PartRandomSample Operator Plugin"
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
include /usr/local/visit/visit-1.5.1/1.5.1/linux-x86_64/include/make-targets
