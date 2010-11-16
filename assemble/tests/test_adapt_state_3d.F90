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
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
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

subroutine test_adapt_state_3d

  use adapt_state_module
  use field_options
  use fields
  use metric_assemble
  use reserve_state_module
  use spud
  use state_module
  use unittest_tools
  use vtk_interfaces
  use global_parameters
  
  implicit none

  type(mesh_type), pointer :: mesh
  type(scalar_field) :: pressure
  type(state_type) :: state, state_array(1), state_read
  type(vector_field) :: velocity
  type(vector_field), pointer :: mesh_field
  type(tensor_field) :: metric

  integer :: i, stat
  
  call vtk_read_state("data/pseudo2d.vtu", state_read)

  mesh_field => extract_vector_field(state_read, "Coordinate")

  mesh => extract_mesh(state_read, "Mesh")
  mesh%name = "CoordinateMesh"
  mesh%option_path = "/geometry/mesh"
  call add_faces(mesh)
  mesh_field%mesh = mesh

  call insert(state, mesh, "CoordinateMesh")
  call insert(state, mesh_field, "Coordinate")

  adaptivity_mesh_name = "CoordinateMesh"
  topology_mesh_name = "CoordinateMesh"

  call deallocate(state_read)

  mesh_field => extract_vector_field(state, "Coordinate")
  mesh => extract_mesh(state, "CoordinateMesh")
 
  call allocate(pressure, mesh, "Pressure")
  call allocate(velocity, mesh_dim(mesh), mesh, "Velocity")
  
  do i = 1, node_count(mesh)
    call set(pressure, i, mesh_field%val(1,i) ** 2.0)
    call set(velocity, i, node_val(mesh_field, i))
  end do
  
  call adaptivity_options(state, pressure, 1.0, .false.)

  call insert(state, pressure, "Pressure")
  call insert(state, velocity, "Velocity")
  call deallocate(pressure)
  call deallocate(velocity)

  state_array(1) = state
  call create_reserve_state(state_array)

  call set_option_attribute("/geometry/mesh/name", "CoordinateMesh", stat = stat)
  assert(stat == SPUD_NEW_KEY_WARNING)
  call add_option("/geometry/mesh/from_file", stat = stat)
  assert(stat == SPUD_NEW_KEY_WARNING)

  call set_option_attribute("/material_phase/name", "MaterialPhase", stat = stat)
  assert(stat == SPUD_NEW_KEY_WARNING)

  call set_option("/mesh_adaptivity/hr_adaptivity/maximum_number_of_nodes", 100000, stat = stat)
  assert(stat == SPUD_NEW_KEY_WARNING)

  call set_option("/geometry/dimension", 3, stat = stat)
  assert(stat == SPUD_NEW_KEY_WARNING)
  call adaptivity_bounds(state_array(1), 0.01, 1.0, name = "CoordinateMesh")

  call allocate(metric, mesh, "Metric")
  call assemble_metric(state_array, metric)
  
  call adapt_state(state_array, metric)
  call report_test("[adapt_state]", .false., .false., "adapt_state failure")
  state = state_array(1)

  mesh_field => extract_vector_field(state, "Coordinate")
  call vtk_write_fields("data/test_adapt_state_3d_out", 0, mesh_field, mesh_field%mesh) 

  call deallocate(state)
  
  call report_test_no_references()

end subroutine test_adapt_state_3d
