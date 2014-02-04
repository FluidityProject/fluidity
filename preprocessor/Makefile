#    Copyright (C) 2006 Imperial College London and others.
#
#    Please see the AUTHORS file in the main source directory for a full list
#    of copyright holders.
#
#    Prof. C Pain
#    Applied Modelling and Computation Group
#    Department of Earth Science and Engineering
#    Imperial College London
#
#    amcgsoftware@imperial.ac.uk
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation,
#    version 2.1 of the License.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#    USA
SHELL = /bin/bash

PACKAGE_NAME = fluidity
FLUIDITY = $(PACKAGE_NAME)

FC      = mpif90
FCFLAGS = -I../include -J../include -I/usr/include/vtk-5.8  -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include  -I/include -I/include/libadjoint -I/usr/include/ -I/usr/local/include/ -I/usr/include/libadjoint -I/usr/local/include/libadjoint -I/usr/lib/petscdir/3.2/include -I/usr/lib/petscdir/3.2/linux-gnu-c-opt/include -I/usr/include -I/usr/include/scotch -I/usr/include/spooles -I/usr/include/suitesparse -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -DHAVE_PETSC -DHAVE_VTK=1 -ffast-math -frecord-marker=4  -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -fno-realloc-lhs -ffree-line-length-none -ffixed-line-length-none  -O3 -I/usr/lib/petscdir/3.2/include -I/usr/lib/petscdir/3.2/linux-gnu-c-opt/include -I/usr/include -I/usr/include/scotch -I/usr/include/spooles -I/usr/include/suitesparse -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -fdefault-real-8 -fdefault-double-8 -I/usr/local/include -I./

MPIF90  = mpif90

CC      = mpicc
CFLAGS  = -I/usr/include/vtk-5.8  -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include  -I/include -I/include/libadjoint -I/usr/include/ -I/usr/local/include/ -I/usr/include/libadjoint -I/usr/local/include/libadjoint -I/usr/lib/petscdir/3.2/include -I/usr/lib/petscdir/3.2/linux-gnu-c-opt/include -I/usr/include -I/usr/include/scotch -I/usr/include/spooles -I/usr/include/suitesparse -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -DHAVE_PETSC -DHAVE_VTK=1  -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I../include -I/usr/local/include 

CXX     = mpicxx
CXXFLAGS= -I/usr/include/vtk-5.8  -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include  -I/include -I/include/libadjoint -I/usr/include/ -I/usr/local/include/ -I/usr/include/libadjoint -I/usr/local/include/libadjoint -I/usr/lib/petscdir/3.2/include -I/usr/lib/petscdir/3.2/linux-gnu-c-opt/include -I/usr/include -I/usr/include/scotch -I/usr/include/spooles -I/usr/include/suitesparse -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -DHAVE_PETSC -DHAVE_VTK=1  -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I../include -I/usr/local/include

MAKE    = make

AR = ar
ARFLAGS = cr

LIB = ../lib/lib$(FLUIDITY).a

LIBS    = -L../lib ../lib/libvtkfortran.a -l$(FLUIDITY) -lmba2d -lvtkIO -lvtkGraphics -lvtkFiltering -lvtkCommon -lvtksys -ldl -L/usr/lib/petscdir/3.2/linux-gnu-c-opt/lib -L/usr/lib/petscdir/3.2/linux-gnu-c-opt/lib -L/usr/lib/petscdir/3.2/linux-gnu-c-opt/lib -lpetsc -lX11 -L/usr/lib -ldmumps -lzmumps -lsmumps -lcmumps -lmumps_common -lpord -lscalapack-openmpi -lblacsCinit-openmpi -lblacs-openmpi -lptesmumps -lptscotch -lptscotcherr -lpthread -lspooles -lHYPRE_utilities -lHYPRE_struct_mv -lHYPRE_struct_ls -lHYPRE_sstruct_mv -lHYPRE_sstruct_ls -lHYPRE_IJ_mv -lHYPRE_parcsr_ls -lcholmod -lumfpack -lamd -llapack -lblas -lhdf5_fortran -lhdf5 -lz -L/usr/lib/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.7 -L/usr/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -ldl -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -lmpi_f90 -lmpi_f77 -lgfortran -lm -lgfortran -lm -lgfortran -lm -lm -lquadmath -lm -lmpi_cxx -lstdc++ -ldl -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -ldl -lzoltan -lparmetis -lmetis -lnetcdf -llapack -ludunits2 -lpthread -lm -lstdc++ -L/usr/lib -lpython2.7 -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions -L/usr/lib -lz -lpthread -ldl -lutil -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lgfortran -lm -lquadmath -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lgfortran -lm -lquadmath -L/lib -L/usr/lib -L/usr/local/lib/ -ladjoint -L./lib  -lblas  -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lgfortran -lm -lquadmath

LINKER  = mpicxx    -I/usr/include/vtk-5.8  -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include  -I/include -I/include/libadjoint -I/usr/include/ -I/usr/local/include/ -I/usr/include/libadjoint -I/usr/local/include/libadjoint -I/usr/lib/petscdir/3.2/include -I/usr/lib/petscdir/3.2/linux-gnu-c-opt/include -I/usr/include -I/usr/include/scotch -I/usr/include/spooles -I/usr/include/suitesparse -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -DHAVE_PETSC -DHAVE_VTK=1 

TESTLIBS = $(shell echo -lmba2d -lvtkIO -lvtkGraphics -lvtkFiltering -lvtkCommon -lvtksys -ldl -L/usr/lib/petscdir/3.2/linux-gnu-c-opt/lib -L/usr/lib/petscdir/3.2/linux-gnu-c-opt/lib -L/usr/lib/petscdir/3.2/linux-gnu-c-opt/lib -lpetsc -lX11 -L/usr/lib -ldmumps -lzmumps -lsmumps -lcmumps -lmumps_common -lpord -lscalapack-openmpi -lblacsCinit-openmpi -lblacs-openmpi -lptesmumps -lptscotch -lptscotcherr -lpthread -lspooles -lHYPRE_utilities -lHYPRE_struct_mv -lHYPRE_struct_ls -lHYPRE_sstruct_mv -lHYPRE_sstruct_ls -lHYPRE_IJ_mv -lHYPRE_parcsr_ls -lcholmod -lumfpack -lamd -llapack -lblas -lhdf5_fortran -lhdf5 -lz -L/usr/lib/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.7 -L/usr/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -ldl -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -lmpi_f90 -lmpi_f77 -lgfortran -lm -lgfortran -lm -lgfortran -lm -lm -lquadmath -lm -lmpi_cxx -lstdc++ -ldl -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -ldl -lzoltan -lparmetis -lmetis -lnetcdf -llapack -ludunits2 -lpthread -lm -lstdc++ -L/usr/lib -lpython2.7 -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions -L/usr/lib -lz -lpthread -ldl -lutil -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lgfortran -lm -lquadmath -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lgfortran -lm -lquadmath -L/lib -L/usr/lib -L/usr/local/lib/ -ladjoint -L./lib  | sed 's@./lib/lib\([a-z]*\)\.a@-l\1@g')

OBJS = Populate_State.o Populate_Sub_State.o Reserve_State.o Field_Priority_Lists.o\
       Physics_From_Options.o VTK_Cache.o Initialise_Fields.o \
			 Boundary_Conditions_From_Options.o synthetic_bc.o Tidal_Modelling.o

.SUFFIXES: .f90 .F90 .cpp .c .o .a

# cancel implicit rule that tries to make .o from .mod with m2c:
%.o: %.mod

.F90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) $(GENFLAGS) -c $< 
.c.o:
	@echo "        CC $<"
	$(CC) $(CFLAGS) -c $<
.cpp.o:
	@echo "        CXX $<"
	$(CXX) $(CXXFLAGS) -c $<

$(LIB): $(OBJS)
	@mkdir -p ../lib

default: $(LIB)

# This gets run last by the fluidity build process.
check_options.o: check_options.F90
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $<
	mkdir -p ../lib


# This gets run last by the fluidity build process.
register_diagnostics.o: register_diagnostics.F90
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $<
	mkdir -p ../lib


clean:
	rm -f *.o *.d *.mod

include Makefile.dependencies
