!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

subroutine test_potential_vorticity

  use fields
  use fldebug
  use state_module
  use unittest_tools
  use vorticity_diagnostics
  use vtk_interfaces
  use spud

  implicit none

  character(len = 64) :: buffer
  integer :: i, stat
  real :: max_val
  real, dimension(3) :: pos
  type(mesh_type), pointer :: mesh
  type(scalar_field) :: perturbation_density, pv
  type(state_type) :: state
  type(vector_field) :: velocity
  type(vector_field), pointer :: positions

  call vtk_read_state("data/cube-itv5.vtu", state)

  positions => extract_vector_field(state, "Coordinate")
  assert(positions%dim == 3)
  mesh => positions%mesh
  
  call set_option("/geometry/dimension", 3, stat=stat)
  call set_option("/physical_parameters/coriolis/f_plane/f", 1.0, stat=stat)
  
  call allocate(velocity, positions%dim, mesh, "Velocity")
  call zero(velocity)
  do i = 1, node_count(velocity)
    pos = node_val(positions, i)
    call set(velocity, i, 0.5 * (/-pos(2), pos(1), 0.0/))
  end do
  call insert(state, velocity, velocity%name)
  
  call allocate(perturbation_density, mesh, "PerturbationDensity")
  call zero(perturbation_density)
  do i = 1, node_count(perturbation_density)
    pos = node_val(positions, i)
    call set(perturbation_density, i, -pos(3))
  end do
  call insert(state, perturbation_density, perturbation_density%name)

  call allocate(pv, mesh, "PotentialVorticity")
  call calculate_potential_vorticity(state, pv)

  call vtk_write_fields("data/test_potential_vorticity_out", &
    & position = positions, model = mesh, &
    & sfields = (/perturbation_density, pv/), &
    & vfields = (/velocity/))
    
  call deallocate(velocity)
  call deallocate(perturbation_density)

  max_val = 0.0
  do i = 1, node_count(pv)
    max_val = max(max_val, abs(node_val(pv, i) + 2.0))
  end do

  write(buffer, *) max_val
  call report_test("[PV == -2.0]", &
    & fnequals(max_val, 0.0, tol = 100.0 * spacing(2.0)), .false., &
    & "PV /= -2.0 - Max. abs. diff: " // buffer)
    
  call deallocate(pv)
  call deallocate(state)
  
  call report_test_no_references()
  
end subroutine test_potential_vorticity
