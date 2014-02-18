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
FCFLAGS = -I../include -ffast-math -frecord-marker=4  -fno-realloc-lhs -ffree-line-length-none -ffixed-line-length-none  -O3 -I/usr/lib/petscdir/3.2/include -I/usr/lib/petscdir/3.2/linux-gnu-c-opt/include -I/usr/include -I/usr/include/scotch -I/usr/include/spooles -I/usr/include/suitesparse -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -fdefault-real-8 -fdefault-double-8

F77      = mpif90
FFLAGS = -I../include -ffast-math -frecord-marker=4  -fno-realloc-lhs -ffree-line-length-none -ffixed-line-length-none  -O3 -fdefault-real-8 -fdefault-double-8

UTFC    = ../bin/unittest_f90

CC      = mpicc
CFLAGS  =  -I/usr/include/python2.7

CXX     = mpicxx
CXXFLAGS=  -I/usr/include/python2.7

MAKE    = make
AR  = ar
ARFLAGS = cr

LIB = ../lib/libmba2d.a

TESTLIBS = $(shell echo -lmba2d -lvtkIO -lvtkGraphics -lvtkFiltering -lvtkCommon -lvtksys -ldl -L/usr/lib/petscdir/3.2/linux-gnu-c-opt/lib -L/usr/lib/petscdir/3.2/linux-gnu-c-opt/lib -L/usr/lib/petscdir/3.2/linux-gnu-c-opt/lib -lpetsc -lX11 -L/usr/lib -ldmumps -lzmumps -lsmumps -lcmumps -lmumps_common -lpord -lscalapack-openmpi -lblacsCinit-openmpi -lblacs-openmpi -lptesmumps -lptscotch -lptscotcherr -lpthread -lspooles -lHYPRE_utilities -lHYPRE_struct_mv -lHYPRE_struct_ls -lHYPRE_sstruct_mv -lHYPRE_sstruct_ls -lHYPRE_IJ_mv -lHYPRE_parcsr_ls -lcholmod -lumfpack -lamd -llapack -lblas -lhdf5_fortran -lhdf5 -lz -L/usr/lib/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.7 -L/usr/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -ldl -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -lmpi_f90 -lmpi_f77 -lgfortran -lm -lgfortran -lm -lgfortran -lm -lm -lquadmath -lm -lmpi_cxx -lstdc++ -ldl -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -ldl -lnetcdf -llapack -lpthread -lm -lstdc++ -L/usr/lib -lpython2.7 -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions -L/usr/lib -lz -lpthread -ldl -lutil -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lgfortran -lm -lquadmath -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lgfortran -lm -lquadmath -L/lib -L/usr/lib -L/usr/local/lib/ -ladjoint -L./lib  | sed 's@./lib/lib\([a-z]*\)\.a@-l\1@g')

# These are all the objects that are supposed to go into the makefile
OBJS = ani2.o auxSE.o auxSF.o auxSP.o calcrv.o clpsF.o delaunay.o dsort.o error.o forlibmba.o \
       insrtP.o lintrp2D.o list_bsct.o loadM.o makM.o makQ_linear.o makS.o mba_analytic.o \
       mba_nodal.o minim.o moveP.o nlnfnc.o prjcrv.o saveM.o smoothing.o splitE.o \
       swapF.o time.o untangle.o update.o utils.o ZZ.o isnan.o isnan_f.o

.SUFFIXES: .f90 .F90 .c .o .a .f

.f.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $<
.f90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $<
.F90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $<
.cpp.o:
	@echo "        CXX $<"
	$(CXX) $(CXXFLAGS) -I../include -c $<

default: $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)
	cp *.mod ../include/

clean:
	rm -f *.o *.mod

