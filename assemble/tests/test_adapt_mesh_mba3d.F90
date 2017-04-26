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

subroutine test_adapt_mesh_mba3d

  use field_options
  use fields
  use limit_metric_module
  use mba3d_integration
  use metric_assemble
  use mesh_files
  use spud
  use state_module
  use unittest_tools
  use vtk_interfaces
  
  implicit none

#ifdef HAVE_MBA_3D
  type(mesh_type), pointer :: mesh
  type(scalar_field) :: pressure
  type(state_type) :: state, state_array(1), state_read
  type(vector_field) :: output_mesh_field, velocity
  type(vector_field), target :: input_mesh_field
  type(tensor_field) :: metric

  integer :: expected_eles, i, stat

  input_mesh_field = read_mesh_files("data/cube_unstructured", quad_degree = 1, format="gmsh")
  
  mesh => input_mesh_field%mesh
  mesh%name = "CoordinateMesh"

  call insert(state, mesh, "CoordinateMesh")
  call insert(state, mesh, "Mesh")
  call insert(state, input_mesh_field, "Coordinate")

  call deallocate(state_read)

  mesh => extract_mesh(state, "CoordinateMesh")

  call allocate(pressure, mesh, "Pressure")
  call allocate(velocity, mesh_dim(mesh), mesh, "Velocity")
  
  do i = 1, node_count(mesh)
    call set(pressure, i, input_mesh_field%val(1,i) ** 2.0)
    call set(velocity, i, node_val(input_mesh_field, i))
  end do
  
  call adaptivity_options(state, pressure, 0.1, .false.)

  call insert(state, pressure, "Pressure")
  call insert(state, velocity, "Velocity")
  call deallocate(pressure)
  call deallocate(velocity)

  call set_option("/mesh_adaptivity/hr_adaptivity/maximum_number_of_nodes", 100000, stat = stat)
  assert(stat == SPUD_NEW_KEY_WARNING)
  call adaptivity_bounds(state, 0.01, 1.0)

  call allocate(metric, mesh, "Metric")
  state_array(1) = state
  call assemble_metric(state_array, metric)
  
  call adapt_mesh_mba3d(input_mesh_field, metric, output_mesh_field)
  call report_test("[adapt_mesh_mba3d]", .false., .false., "adapt_mesh_mba3d failure")
  expected_eles = expected_elements(input_mesh_field, metric)
  call report_test("[expected_elements]", abs(float(ele_count(output_mesh_field) - expected_eles) / float(expected_eles)) > 0.25, .false., "Incorrect output mesh element count")

  call vtk_write_fields("data/test_adapt_mesh_mb3d_out", 0, output_mesh_field, output_mesh_field%mesh) 

  call deallocate(input_mesh_field)
  call deallocate(output_mesh_field)
  call deallocate(metric)
  call deallocate(state)
  
  call report_test_no_references()
#else
  call report_test("[dummy]", .false., .false., "Dummy")
#endif

end subroutine test_adapt_mesh_mba3d
