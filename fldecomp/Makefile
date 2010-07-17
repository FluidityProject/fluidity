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
FCFLAGS = -I../include -ffast-math -frecord-marker=4 -O3 -I/home/angus/local/include -I/usr/include/vtk-5.2 -ffixed-line-length-none -ffree-line-length-none  -DHAVE_NUMPY -I/usr/lib/python2.6/dist-packages/numpy/core/include  -O3 -I/usr/lib64/petscdir/3.0.0//linux-gnu-c-opt/include -I/usr/lib64/petscdir/3.0.0//include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/usr/include/spooles -I/usr/include -I/usr/include/superlu -I/usr/include/suitesparse -I/usr/include/scotch -fdefault-real-8

MPIF90  = mpif90

CC      = mpicc
CFLAGS  = -m64 -O4 -I/home/angus/local/include -I/usr/include/vtk-5.2 -I/usr/include/python2.6 -DHAVE_NUMPY -I/usr/lib/python2.6/dist-packages/numpy/core/include  -I../include

CXX     = mpicxx
CXXFLAGS= -m64 -O4 -I/home/angus/local/include -I/usr/include/vtk-5.2 -I/usr/include/python2.6 -DHAVE_NUMPY -I/usr/lib/python2.6/dist-packages/numpy/core/include  -I../include

MAKE    = make
AR  = ar
ARFLAGS = cr

LINKER  = mpicxx -L/home/angus/local/lib -L/usr/lib64/petscdir/3.0.0//linux-gnu-c-opt/lib  

LIB_FLUIDITY = ../lib/libfluidity.a
LIBS    = -L../lib $(LIB_FLUIDITY)  -lf77blas -latlas  -L/home/angus/local/lib -L/usr/lib64/petscdir/3.0.0//linux-gnu-c-opt/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.3/../../../../lib -L/lib/../lib -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.3/../../.. -L/usr/lib/x86_64-linux-gnu -lgfortranbegin -lgfortran -lm  \
	../lib/libadaptivity.a ../lib/libvtkfortran.a \
	../lib/libspatialindex.a ../lib/libspud.a

OBJS = fldmain.o fldgmsh.o fldtriangle.o

FLDECOMP = ../bin/fldecomp

.SUFFIXES: .f90 .F90 .c .o .a

.f90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $< || (rm $(LIB) && false)
.F90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $< || (rm $(LIB) && false)
.cpp.o:
	@echo "        CXX $<"
	$(CXX) $(CXXFLAGS) -c $< || (rm $(LIB) && false)

fldecomp: $(LIB) $(OBJS) 
	@$(EVAL) $(LINKER)  -o $(FLDECOMP) -I../include/ \
	    $(OBJS) $(LIBS) -lmetis

clean:
	rm -f *.o *~
