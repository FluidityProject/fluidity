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
FCFLAGS = -I../../include -fbounds-check -Wall -Wimplicit-interface -Wno-surprising -g -O0 -frecord-marker=4 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -ffree-line-length-none -ffixed-line-length-none -finit-integer=-66666 -finit-real=nan  -I/usr/lib/petscdir/3.1/linux-gnu-c-opt/include -I/usr/lib/petscdir/3.1/include -I/usr/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/usr/include/spooles -I/usr/include/scotch -I/usr/lib -I/usr/include/suitesparse -fdefault-real-8 -fdefault-double-8 -I../

CC  = mpicc
CFLAGS  = -fbounds-check -Wall -g -O0 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I../../include

CXX = mpicxx -I../../include
# Add -DGRAPHIC to ouput CSV files for the fluxes tests
CXXFLAGS= -I/usr/include/vtk-5.8  -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I/usr/include/ -I/usr/local/include/  -I/include -I/include/libadjoint -I/usr/include/ -I/usr/local/include/ -I/usr/include/libadjoint -I/usr/local/include/libadjoint -I/usr/lib/petscdir/3.1/linux-gnu-c-opt/include -I/usr/lib/petscdir/3.1/include -I/usr/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/usr/include/spooles -I/usr/include/scotch -I/usr/lib -I/usr/include/suitesparse -DHAVE_PETSC -DHAVE_VTK=1 -fbounds-check -Wall -g -O0 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I../../include

LIBS = ./lib/libvtkfortran.a -lvtkIO -lvtkGraphics -lvtkFiltering -lvtkCommon -lvtksys -ldl -larpack -lnetcdff -lnetcdf -llapack -lblas -ludunits2 -lpthread -lm -lstdc++ -L/usr/lib -lpython2.7 -Xlinker -export-dynamic -Wl,-O0 -Wl,-Bsymbolic-functions -lssl -lcrypto -lssl -lcrypto -L/usr/lib -lz -lpthread -ldl -lutil -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lquadmath -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lquadmath -lnetcdf -lnetcdff -L/lib -lparmetis -lmetis -L/usr/lib -L/usr/local/lib/ -lzoltan -lparmetis -lmetis -L/lib -L/usr/lib -L/usr/local/lib/ -ladjoint -L/usr/lib/petscdir/3.1/linux-gnu-c-opt/lib -L/usr/lib/petscdir/3.1/linux-gnu-c-opt/lib -L/usr/lib/petscdir/3.1/linux-gnu-c-opt/lib -lpetsc -lX11 -L/usr/lib -ldmumps -lzmumps -lsmumps -lcmumps -lmumps_common -lpord -lscalapack-openmpi -lblacsCinit-openmpi -lblacs-openmpi -lspooles -lHYPRE -lHYPRE_utilities -lHYPRE_multivector -lHYPRE_krylov -lHYPRE_struct_mv -lHYPRE_struct_ls -lHYPRE_sstruct_mv -lHYPRE_sstruct_ls -lHYPRE_seq_mv -lHYPRE_parcsr_mv -lHYPRE_parcsr_block_mv -lHYPRE_DistributedMatrix -lHYPRE_MatrixMatrix -lHYPRE_IJ_mv -lHYPRE_Euclid -lHYPRE_ParaSails -lHYPRE_DistributedMatrixPilutSolver -lHYPRE_parcsr_ls -lHYPRE_FEI_fgmres -lHYPRE_mli -lHYPRE_FEI -lscotch -lhdf5_fortran -lhdf5 -lz -lumfpack -lamd -llapack -lblas -L/usr/lib/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -ldl -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -lmpi_f90 -lmpi_f77 -lgfortran -lm -lgfortran -lm -lgfortran -lm -lm -lquadmath -lm -lmpi_cxx -lstdc++ -ldl -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -ldl -L./lib  ./lib/libspud.a ./lib/libjudy.a ./lib/libspatialindex.a  -llapack -lblas  -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lquadmath
LIBFLUIDITY=../../lib/lib$(FLUIDITY).a

FORTRANTESTS= $(basename $(wildcard *.F90))


CTESTS= $(basename $(wildcard *.c)) $(basename $(wildcard *.cpp)) 

TESTS=$(CTESTS) $(FORTRANTESTS)

# the test binaries NOT to be built
DISABLED_TESTS= test_main

# the test programs to be built:
TEST_BINARIES=$(filter-out $(DISABLED_TESTS), $(TESTS))

unittest: $(TEST_BINARIES)
	ln -sf $(addprefix ../../ocean_forcing/tests/,$(TEST_BINARIES)) ../../bin/tests/

.SUFFIXES: .f90 .F90 .c .cpp .o

%.o:	%.f90
	$(FC) $(FCFLAGS) -c $<
%.o:	%.F90
	$(FC) $(FCFLAGS) -c $<
%.o:	%.c
	$(CC) $(CFLAGS) -c $<
%.o:	%.cpp
	$(CXX) $(CXXFLAGS) -c $<
%: %.cpp

# creates a TESTNAME_main_f.o from test_main_f.cpp which calls the subroutine
# TESTNAME, that should be a subroutine in TESTNAME.F90
%_main.o: test_main.cpp
	$(CXX) $(CXXFLAGS) -D TESTNAME=$(subst _main.o,,$@)_ -o $@ -c $^ 

# link this TESTNAME_main.o with TESTNAME.o from TESTNAME.F90
test_python_state:  python_state.o test_python_state_main.o test_python_state.o python_statec.o lib/
	$(CXX) $(LDFLAGS) -o $@ $(filter %.o,$^) $(LIBS) 

# link this TESTNAME_main_f.o with TESTNAME.o from TESTNAME.F90
%: %_main.o %.o lib/
	$(CXX) $(LDFLAGS) -o $@ $(filter %.o,$^) $(LIBFLUIDITY) $(LIBS) 

# make a temp sym-link to the directory containing fluidity, adapt and sam libs
lib/:
	ln -s ../../lib

clean:
	rm -f $(TEST_BINARIES)
	rm -f *.o *.mod
	rm -f lib
