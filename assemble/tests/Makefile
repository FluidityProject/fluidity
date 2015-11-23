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
FCFLAGS = -I../../include -fbounds-check -Wall -Wimplicit-interface -Wno-surprising -g -O0 -frecord-marker=4 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -fno-realloc-lhs -ffree-line-length-none -ffixed-line-length-none -finit-integer=-66666 -finit-real=nan -I/usr/lib/petscdir/3.4.2/include -I/usr/lib/petscdir/3.4.2/linux-gnu-c-opt/include -I/usr/include -I/usr/include/suitesparse -I/usr/include/scotch -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -Iinclude/ -fdefault-real-8 -fdefault-double-8 \
  -I../

CC  = mpicc
CFLAGS  = -fbounds-check -Wall -g -O0 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include  -I../../include

CXX = mpicxx -I../../include
CXXFLAGS= -I/usr/include/vtk-5.10 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I/usr/lib/petscdir/3.4.2/include -I/usr/lib/petscdir/3.4.2/linux-gnu-c-opt/include -I/usr/include -I/usr/include/suitesparse -I/usr/include/scotch -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -Iinclude/ -DHAVE_PETSC  -I/usr/include/vtk-5.10 -DHAVE_VTK=1 -DHAVE_VTK=1 -fbounds-check -Wall -g -O0 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include  -I../include -I../../include 

LDFLAGS  = -rdynamic   

LIBS = ./lib/libadaptivity.a ./lib/libvtkfortran.a ./lib/libspatialindex.a ./lib/libspud.a ./lib/libjudy.a -lmba2d -lvtkIO -lvtkGraphics -lvtkFiltering -lvtkCommon -ldl -lzoltan -lparmetis -lmetis -L/usr/lib/petscdir/3.4.2/linux-gnu-c-opt/lib -L/usr/lib/petscdir/3.4.2/linux-gnu-c-opt/lib -L/usr/lib/petscdir/3.4.2/linux-gnu-c-opt/lib -lpetsc -L/usr/lib -lHYPRE_utilities -lHYPRE_struct_mv -lHYPRE_struct_ls -lHYPRE_sstruct_mv -lHYPRE_sstruct_ls -lHYPRE_IJ_mv -lHYPRE_parcsr_ls -ldmumps -lzmumps -lsmumps -lcmumps -lmumps_common -lpord -lscalapack-openmpi -lcholmod -lumfpack -lamd -llapack -lblas -lX11 -lpthread -lptesmumps -lptscotch -lptscotcherr -L/usr/lib/x86_64-linux-gnu -lfftw3 -lfftw3_mpi -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lm -L/usr/lib/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/5 -L/lib/x86_64-linux-gnu -lmpi_f90 -lmpi_f77 -lgfortran -lm -lgfortran -lm -lquadmath -lm -lmpi_cxx -lstdc++ -ldl -lmpi -lhwloc -lgcc_s -lpthread -ldl -lnetcdf -llapack -ludunits2 -lpthread -lm -lstdc++ -L/usr/lib -lpython2.7 -Xlinker -export-dynamic -Wl,-O0 -Wl,-Bsymbolic-functions -L/usr/lib -lz -lpthread -ldl -lutil -L/usr/lib/gcc/x86_64-linux-gnu/5 -L/usr/lib/gcc/x86_64-linux-gnu/5/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/5/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/5/../../.. -lgfortran -lm -lquadmath -L/usr/lib/gcc/x86_64-linux-gnu/5 -L/usr/lib/gcc/x86_64-linux-gnu/5/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/5/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/5/../../.. -lgfortran -lm -lquadmath -L./lib 
FLIBS =  -L/usr/lib/gcc/x86_64-linux-gnu/5 -L/usr/lib/gcc/x86_64-linux-gnu/5/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/5/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/5/../../.. -lgfortran -lm -lquadmath
LIBFLUIDITY=../../lib/lib$(FLUIDITY).a

# the test binaries NOT to be built
DISABLED_TESTS=test_mba_supermesh test_legacy_cv_faces test_lagrangian_remap test_adapt_mesh_mba3d test_pressure_solve test_matrix_free \
               test_dg_interpolation_quads test_dg_diffusion

# the test programs to be built:
TEST_BINARIES=$(filter-out $(DISABLED_TESTS), $(basename $(wildcard *.F90)))

unittest: $(TEST_BINARIES) data
	ln -sf $(addprefix ../../assemble/tests/,$(TEST_BINARIES)) ../../bin/tests/

data:
	cd ../../tests/data; \
	gmsh -3 cube_prismatic.geo; \
	../../bin/gmsh2triangle cube_prismatic.msh; \
	gmsh -3 cube_unstructured.geo; \
	../../bin/gmsh2triangle cube_unstructured.msh; \
	cd -


.SUFFIXES: .f90 .F90 .c .cpp .o .a $(.SUFFIXES)

%.o:	%.f90
	$(FC) $(FCFLAGS) -c $<
%.o:	%.F90
	$(FC) $(FCFLAGS) -c $<
%.o:	%.c
	$(CC) $(CFLAGS) -c $<
%.o:	%.cpp
	$(CXX) $(CXXFLAGS) -c $<

# creates a TESTNAME_main.o from test_main.cpp which calls the subroutine
# TESTNAME, that should be a subroutine in TESTNAME.F90
%_main.o:
	$(CXX) $(CXXFLAGS) -D TESTNAME=$(subst _main.o,,$@)_ -o $@ -c test_main.cpp 

# link this TESTNAME_main.o with TESTNAME.o from TESTNAME.F90
%: %_main.o %.o lib/
	$(CXX) $(LDFLAGS) -o $@ $(filter %.o,$^) $(LIBFLUIDITY) $(LIBS) 

# make a temp sym-link to the directory containing fluidity and adaptivity
lib/:
	ln -s ../../lib

clean:
	rm -f $(TEST_BINARIES)
	rm -f *.o *.mod
	rm -f lib
