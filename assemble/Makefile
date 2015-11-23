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
FCFLAGS = -I../include -J../include -I/usr/include/vtk-5.10 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I/usr/lib/petscdir/3.4.2/include -I/usr/lib/petscdir/3.4.2/linux-gnu-c-opt/include -I/usr/include -I/usr/include/suitesparse -I/usr/include/scotch -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -Iinclude/ -DHAVE_PETSC  -I/usr/include/vtk-5.10 -DHAVE_VTK=1 -DHAVE_VTK=1 -fbounds-check -Wall -Wimplicit-interface -Wno-surprising -g -O0 -frecord-marker=4 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -fno-realloc-lhs -ffree-line-length-none -ffixed-line-length-none -finit-integer=-66666 -finit-real=nan -I/usr/lib/petscdir/3.4.2/include -I/usr/lib/petscdir/3.4.2/linux-gnu-c-opt/include -I/usr/include -I/usr/include/suitesparse -I/usr/include/scotch -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -Iinclude/ -fdefault-real-8 -fdefault-double-8 -I/usr/local/include -I./

MPIF90  = mpif90

CC      = mpicc
CFLAGS  = -I/usr/include/vtk-5.10 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I/usr/lib/petscdir/3.4.2/include -I/usr/lib/petscdir/3.4.2/linux-gnu-c-opt/include -I/usr/include -I/usr/include/suitesparse -I/usr/include/scotch -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -Iinclude/ -DHAVE_PETSC  -I/usr/include/vtk-5.10 -DHAVE_VTK=1 -DHAVE_VTK=1 -fbounds-check -Wall -g -O0 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include  -I../include -I/usr/local/include -g

CXX     = mpicxx
CXXFLAGS= -I/usr/include/vtk-5.10 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I/usr/lib/petscdir/3.4.2/include -I/usr/lib/petscdir/3.4.2/linux-gnu-c-opt/include -I/usr/include -I/usr/include/suitesparse -I/usr/include/scotch -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -Iinclude/ -DHAVE_PETSC  -I/usr/include/vtk-5.10 -DHAVE_VTK=1 -DHAVE_VTK=1 -fbounds-check -Wall -g -O0 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include  -I../include -I/usr/local/include

MAKE    = make
AR  = ar
ARFLAGS = cr

LIB = ../lib/lib$(FLUIDITY).a
OBJS    = Momentum_Diagnostic_Fields.o \
	Momentum_Equation.o Momentum_CG.o \
	qmesh.o \
	OceanSurfaceForcing.o \
	Momentum_DG.o Advection_Diffusion_CG.o Advection_Diffusion_DG.o Advection_Diffusion_FV.o \
	Solid_Configuration.o Explicit_ALE.o \
	MeshMovement.o \
	Field_Equations_CV.o Assemble_CMC.o \
	Divergence_Matrix_CG.o Gradient_Matrix_CG.o Divergence_Matrix_CV.o \
	State_Matrices.o \
	Drag.o    \
	Multimaterials.o Multiphase.o Diagnostic_fields_wrapper.o Diagnostic_Fields_Matrices.o Vorticity_Diagnostics.o \
	Compressible_Projection.o \
	Biology.o \
	Foam_Drainage.o Foam_Flow.o \
	Timeloop_utilities.o Pseudo_supermesh.o \
	Sam_integration.o Solenoidal_interpolation.o \
        Interpolation_manager.o  Discrete_Properties.o Free_Surface.o  \
	Surface_Id_Interleaving.o \
	Geostrophic_Pressure.o Hydrostatic_Pressure.o Vertical_Balance_Pressure.o \
	Adapt_Integration.o Node_Locking.o Mba3d_Integration.o \
	Adapt_State.o Adapt_State_Prescribed.o Adapt_State_Unittest.o \
	Upwind_Stabilisation.o \
        Coriolis.o Full_Projection.o \
        Mba2d_Integration.o Slope_limiters_DG.o \
        Spontaneous_Potentials.o LES.o \
        Saturation_Distribution_Search_HookeJeeves.o \
        Petsc_Solve_State.o Adaptivity_1D.o \
	Zoltan_global_variables.o Zoltan_integration.o Zoltan_callbacks.o Zoltan_detectors.o \
        Diagnostic_Children.o Reduced_Model_Runtime.o Turbine.o Implicit_Solids.o \
				Manifold_Projections.o Hybridized_Helmholtz.o Burgers_Assembly.o Pressure_Dirichlet_BCS_CV.o \
            Porous_Media.o

.SUFFIXES: .F90 .c .cpp .o .a

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

# cancel implicit rule that tries to make .o from .mod with m2c:
%.o: %.mod

clean:
	rm -f *.o *.d *.mod

# Special case mba dependencies
ifeq (yes,yes)
Mba2d_Integration.o ../include/mba2d_integration.mod: \
../include/mba2d_module.mod
endif
ifeq (,yes)
Mba3d_Integration.o ../include/mba3d_integration.mod: \
../include/mba3d_mba_nodal.mod
endif

include Makefile.dependencies

