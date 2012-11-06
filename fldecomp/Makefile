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

FC      = mpif90
FCFLAGS = -I../include -fbounds-check -Wall -Wimplicit-interface -Wno-surprising -g -O0 -frecord-marker=4 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -ffree-line-length-none -ffixed-line-length-none -finit-integer=-66666 -finit-real=nan  -I/usr/lib/petscdir/3.1/linux-gnu-c-opt/include -I/usr/lib/petscdir/3.1/include -I/usr/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/usr/include/spooles -I/usr/include/scotch -I/usr/lib -I/usr/include/suitesparse -fdefault-real-8 -fdefault-double-8

MPIF90  = mpif90

CC      = mpicc
CFLAGS  = -fbounds-check -Wall -g -O0 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include  -I../include

CXX     = mpicxx
CXXFLAGS= -fbounds-check -Wall -g -O0 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include  -I../include

MAKE    = make
AR  = ar
ARFLAGS = cr

LINKER  = mpicxx -rdynamic   

LIB_FLUIDITY = ../lib/libfluidity.a
LIBS    = -L../lib $(LIB_FLUIDITY)  -lblas  -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lquadmath  \
	../lib/libadaptivity.a ../lib/libvtkfortran.a \
	../lib/libspatialindex.a ../lib/libspud.a

OBJS = fldmain.o fldgmsh.o fldtriangle.o partition.o

FLDECOMP = ../bin/fldecomp

.SUFFIXES: .f90 .F90 .c .o .a

.f90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $<
.F90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $<
.cpp.o:
	@echo "        CXX $<"
	$(CXX) $(CXXFLAGS) -c $<

fldecomp: $(LIB) $(OBJS) 
	@$(EVAL) $(LINKER)  -o $(FLDECOMP) -I../include/ \
	    $(OBJS) $(LIBS) -lmetis

clean:
	rm -f *.o *~
