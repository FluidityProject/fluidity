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
FCFLAGS = -I../include -J../include -I/usr/include/vtk-5.8  -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I/usr/include/ -I/usr/local/include/  -I/include -I/include/libadjoint -I/usr/include/ -I/usr/local/include/ -I/usr/include/libadjoint -I/usr/local/include/libadjoint -I/usr/lib/petscdir/3.1/linux-gnu-c-opt/include -I/usr/lib/petscdir/3.1/include -I/usr/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/usr/include/spooles -I/usr/include/scotch -I/usr/lib -I/usr/include/suitesparse -DHAVE_PETSC -DHAVE_VTK=1 -fbounds-check -Wall -Wimplicit-interface -Wno-surprising -g -O0 -frecord-marker=4 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -ffree-line-length-none -ffixed-line-length-none -finit-integer=-66666 -finit-real=nan  -I/usr/lib/petscdir/3.1/linux-gnu-c-opt/include -I/usr/lib/petscdir/3.1/include -I/usr/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/usr/include/spooles -I/usr/include/scotch -I/usr/lib -I/usr/include/suitesparse -fdefault-real-8 -fdefault-double-8 -I/usr/local/include -I./

MPIF90  = mpif90

CC      = mpicc
CFLAGS  = -I/usr/include/vtk-5.8  -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I/usr/include/ -I/usr/local/include/  -I/include -I/include/libadjoint -I/usr/include/ -I/usr/local/include/ -I/usr/include/libadjoint -I/usr/local/include/libadjoint -I/usr/lib/petscdir/3.1/linux-gnu-c-opt/include -I/usr/lib/petscdir/3.1/include -I/usr/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/usr/include/spooles -I/usr/include/scotch -I/usr/lib -I/usr/include/suitesparse -DHAVE_PETSC -DHAVE_VTK=1 -fbounds-check -Wall -g -O0 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I../include -I/usr/local/include

CXX     = mpicxx
CXXFLAGS= -I/usr/include/vtk-5.8  -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I/usr/include/ -I/usr/local/include/  -I/include -I/include/libadjoint -I/usr/include/ -I/usr/local/include/ -I/usr/include/libadjoint -I/usr/local/include/libadjoint -I/usr/lib/petscdir/3.1/linux-gnu-c-opt/include -I/usr/lib/petscdir/3.1/include -I/usr/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/usr/include/spooles -I/usr/include/scotch -I/usr/lib -I/usr/include/suitesparse -DHAVE_PETSC -DHAVE_VTK=1 -fbounds-check -Wall -g -O0 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I../include -I/usr/local/include

MAKE    = make

AR = ar
ARFLAGS = cr

LIB = ../lib/lib$(FLUIDITY).a

OBJS = C++_Debug.o Debug.o Debug_Parameters.o Diagnostic_Output.o

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


$(F90LIB): $(F90LIB_OBJS)
	@mkdir -p ../lib
	$(AR) $(ARFLAGS) $(F90LIB).a $(F90LIB_OBJS)

clean:
	rm -f *.o *.mod

include Makefile.dependencies
