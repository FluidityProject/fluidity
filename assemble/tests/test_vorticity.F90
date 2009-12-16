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
!    C.Pain@Imperial.ac.uk
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

subroutine test_vorticity

  use field_derivatives
  use fields
  use fldebug
  use state_module
  use unittest_tools
  use vorticity_diagnostics
  use vtk_interfaces

  implicit none

  character(len = 64) :: buffer
  integer :: i
  real :: max_val
  real, dimension(3) :: pos
  type(mesh_type), pointer :: mesh
  type(scalar_field) :: s_field
  type(state_type) :: state
  type(vector_field) :: velocity, vorticity
  type(vector_field), pointer :: positions

  call vtk_read_state("data/cube-itv5.vtu", state)

  positions => extract_vector_field(state, "Coordinate")
  assert(positions%dim == 3)
  mesh => positions%mesh
  
  call allocate(s_field, mesh, "BaseScalar")
  call zero(s_field)
  do i = 1, node_count(positions)
    pos = node_val(positions, i)
    call set(s_field, pos(1) + pos(2) + pos(3))
  end do
  call allocate(velocity, positions%dim, mesh, "Velocity")
  call grad(s_field, positions, velocity)
  call deallocate(s_field)
  call insert(state, velocity, velocity%name)

  call allocate(vorticity, mesh_dim(mesh), mesh, "Vorticity")
  call calculate_vorticity(state, vorticity)

  call vtk_write_fields("data/test_vorticity_out", &
    & position = positions, model = mesh, &
    & vfields = (/velocity, vorticity/))
    
  call deallocate(velocity)

  max_val = 0.0
  do i = 1, node_count(vorticity)
    max_val = max(max_val, maxval(node_val(vorticity, i)))
  end do

  write(buffer, *) max_val
    call report_test("[Vorticity: Curl of gradient]", &
      & max_val .fne. 0.0, .false., &
      & "curl(grad) /= 0.0 - Max. component: " // buffer)
    
  call deallocate(vorticity)
  call deallocate(state)
  
  call report_test_no_references()
  
end subroutine test_vorticity
