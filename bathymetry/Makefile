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

LIB = ../lib/lib$(FLUIDITY).a

AR = ar
ARFLAGS = cr

CC  = mpicc
CXX = mpicxx
FC = mpif90

LINKER = $(CXX)

# Adding /usr/include like this is not ideal, see:
# http://gcc.gnu.org/bugzilla/show_bug.cgi?id=35707
# If /usr/include occurs too early it's going to
# interfere with other compilers - and it would completely
# stuff C++ with obsecure, ungoogleable, error messages. 
FCFLAGS   = -I../include @MODFLAG@../include -fbounds-check -Wall -Wimplicit-interface -Wno-surprising -g -O0 -frecord-marker=4 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -ffree-line-length-none -ffixed-line-length-none -finit-integer=-66666 -finit-real=nan  -I/usr/lib/petscdir/3.1/linux-gnu-c-opt/include -I/usr/lib/petscdir/3.1/include -I/usr/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/usr/include/spooles -I/usr/include/scotch -I/usr/lib -I/usr/include/suitesparse -fdefault-real-8 -fdefault-double-8 -I/usr/include/vtk-5.8  -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I/usr/include/ -I/usr/local/include/  -I/include -I/include/libadjoint -I/usr/include/ -I/usr/local/include/ -I/usr/include/libadjoint -I/usr/local/include/libadjoint -I/usr/lib/petscdir/3.1/linux-gnu-c-opt/include -I/usr/lib/petscdir/3.1/include -I/usr/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/usr/include/spooles -I/usr/include/scotch -I/usr/lib -I/usr/include/suitesparse -DHAVE_PETSC -DHAVE_VTK=1 -I/usr/include

CXXFLAGS = -I../include -fbounds-check -Wall -g -O0 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I/usr/include/vtk-5.8  -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I/usr/include/ -I/usr/local/include/  -I/include -I/include/libadjoint -I/usr/include/ -I/usr/local/include/ -I/usr/include/libadjoint -I/usr/local/include/libadjoint -I/usr/lib/petscdir/3.1/linux-gnu-c-opt/include -I/usr/lib/petscdir/3.1/include -I/usr/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/usr/include/spooles -I/usr/include/scotch -I/usr/lib -I/usr/include/suitesparse -DHAVE_PETSC -DHAVE_VTK=1

LDFLAGS   = -rdynamic   
LIBS    = -L../lib -l$(FLUIDITY) -lvtkIO -lvtkGraphics -lvtkFiltering -lvtkCommon -lvtksys -ldl -larpack -lnetcdff -lnetcdf -llapack -lblas -ludunits2 -lpthread -lm -lstdc++ -L/usr/lib -lpython2.7 -Xlinker -export-dynamic -Wl,-O0 -Wl,-Bsymbolic-functions -lssl -lcrypto -lssl -lcrypto -L/usr/lib -lz -lpthread -ldl -lutil -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lquadmath -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lquadmath -lnetcdf -lnetcdff -L/lib -lparmetis -lmetis -L/usr/lib -L/usr/local/lib/ -lzoltan -lparmetis -lmetis -L/lib -L/usr/lib -L/usr/local/lib/ -ladjoint -L/usr/lib/petscdir/3.1/linux-gnu-c-opt/lib -L/usr/lib/petscdir/3.1/linux-gnu-c-opt/lib -L/usr/lib/petscdir/3.1/linux-gnu-c-opt/lib -lpetsc -lX11 -L/usr/lib -ldmumps -lzmumps -lsmumps -lcmumps -lmumps_common -lpord -lscalapack-openmpi -lblacsCinit-openmpi -lblacs-openmpi -lspooles -lHYPRE -lHYPRE_utilities -lHYPRE_multivector -lHYPRE_krylov -lHYPRE_struct_mv -lHYPRE_struct_ls -lHYPRE_sstruct_mv -lHYPRE_sstruct_ls -lHYPRE_seq_mv -lHYPRE_parcsr_mv -lHYPRE_parcsr_block_mv -lHYPRE_DistributedMatrix -lHYPRE_MatrixMatrix -lHYPRE_IJ_mv -lHYPRE_Euclid -lHYPRE_ParaSails -lHYPRE_DistributedMatrixPilutSolver -lHYPRE_parcsr_ls -lHYPRE_FEI_fgmres -lHYPRE_mli -lHYPRE_FEI -lscotch -lhdf5_fortran -lhdf5 -lz -lumfpack -lamd -llapack -lblas -L/usr/lib/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -ldl -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -lmpi_f90 -lmpi_f77 -lgfortran -lm -lgfortran -lm -lgfortran -lm -lm -lquadmath -lm -lmpi_cxx -lstdc++ -ldl -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -ldl -L./lib  -lblas  -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lquadmath \
		  ../lib/libvtkfortran.a ../lib/libspatialindex.a ../lib/libspud.a ../lib/libjudy.a

OBJS = NetCDF_reader.o SampleNetCDF2.o import_bath_data.o read_netcdf_interface.o

.SUFFIXES: .c .o .cpp .F90 .a .so

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

DEPTH_OBJS=import_bath_data.cpp

dtest: $(DEPTH_OBJS)
	$(CXX) $(CXXFLAGS) -DDEPTH_UNITTEST -o $@ $(DEPTH_OBJS) $(LIBFLUIDITY) $(LIBS)

dtest2: $(DEPTH_OBJS)
	$(CXX) $(CXXFLAGS) -DDEPTH_UNITTEST2 -o $@ $(DEPTH_OBJS) $(LIBFLUIDITY) $(LIBS)

install:

clean:
	rm -f *.o *.mod
include Makefile.dependencies
