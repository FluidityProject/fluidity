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
#    C.Pain@Imperial.ac.uk
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

FC      = mpif90
FCFLAGS = -I../include -I/usr/include/vtk-5.2  -I/usr/include/python2.6 -DHAVE_NUMPY -I/usr/lib/python2.6/dist-packages/numpy/core/include -DHAVE_PETSC -I/usr/lib/petscdir/3.0.0/include -I/usr/lib/petscdir/3.0.0/bmake/linux-gnu-c-opt -DHAVE_VTK=1 -ffast-math -frecord-marker=4  -DHAVE_NUMPY -I/usr/lib/python2.6/dist-packages/numpy/core/include -ffree-line-length-none -ffixed-line-length-none  -O3 -I/usr/lib/petscdir/3.0.0/linux-gnu-c-opt/include/ -I/usr/lib/petscdir/3.0.0/include -I/usr/lib/petscdir/3.0.0/bmake/linux-gnu-c-opt -fdefault-real-8 -I/usr/local/include -I./

MPIF90  = mpif90

CC      = mpicc
CFLAGS  = -I/usr/include/vtk-5.2  -I/usr/include/python2.6 -DHAVE_NUMPY -I/usr/lib/python2.6/dist-packages/numpy/core/include -DHAVE_PETSC -I/usr/lib/petscdir/3.0.0/include -I/usr/lib/petscdir/3.0.0/bmake/linux-gnu-c-opt -DHAVE_VTK=1  -I/usr/include/python2.6 -DHAVE_NUMPY -I/usr/lib/python2.6/dist-packages/numpy/core/include -I../include -I/usr/local/include -g

CXX     = mpicxx
CXXFLAGS= -I/usr/include/vtk-5.2  -I/usr/include/python2.6 -DHAVE_NUMPY -I/usr/lib/python2.6/dist-packages/numpy/core/include -DHAVE_PETSC -I/usr/lib/petscdir/3.0.0/include -I/usr/lib/petscdir/3.0.0/bmake/linux-gnu-c-opt -DHAVE_VTK=1  -I/usr/include/python2.6 -DHAVE_NUMPY -I/usr/lib/python2.6/dist-packages/numpy/core/include -I../include -I/usr/local/include

MAKE    = make
AR  = ar
ARFLAGS = cr

LIB = ../lib/libdfluidity.a

#OBJS    = ******.o \
#        ******.o

.SUFFIXES: .f90 .F90 .c .o .a

.f90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $< || (rm $(LIB) && false)
.F90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $< || (rm $(LIB) && false)
.c.o:
	@echo "        CC $<"
	$(CC) $(CFLAGS) -c $< || (rm $(LIB) && false)
.cpp.o:
	@echo "        CXX $<"
	$(CXX) $(CXXFLAGS) -c $< || (rm $(LIB) && false)

default: $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)

clean:
	rm -f *.o *.d *.mod
