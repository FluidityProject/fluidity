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

subroutine test_cfl_number_1d

  use diagnostic_fields
  use fields
  use fldebug
  use manifold_tools
  use state_module
  use unittest_tools

  implicit none

  type(element_type) :: shape
  type(quadrature_type) :: quad
  type(mesh_type) :: coordinate_mesh, velocity_mesh
  type(scalar_field) :: cfl_no
  type(state_type) :: state, local_state
  type(vector_field) :: positions, velocity, positions_3d, local_velocity, velocity_3d
  integer :: j
  
  quad = make_quadrature(vertices = 2, dim  = 1, degree = 2)
  shape = make_element_shape(vertices = 2, dim  = 1, degree = 1, quad = quad)
  call deallocate(quad)
  
  call allocate(coordinate_mesh, nodes = 4, elements = 3, shape = shape, name = "CoordinateMesh")
  call deallocate(shape)
  
  call set_ele_nodes(coordinate_mesh, 1, (/1, 2/))
  call set_ele_nodes(coordinate_mesh, 2, (/2, 3/))
  call set_ele_nodes(coordinate_mesh, 3, (/3, 4/))

  call add_faces(coordinate_mesh)
  
  velocity_mesh = piecewise_constant_mesh(coordinate_mesh, name = "CFLNumberMesh")
  
  call allocate(positions, 1, coordinate_mesh, name = "Coordinate")
  call allocate(positions_3d, 3, coordinate_mesh, name = "Coordinate")
  call allocate(velocity, 1, velocity_mesh, name = "Velocity")
  call allocate(velocity_3d, 3, velocity_mesh, name = "CartesianVelocity")
  call allocate(local_velocity, 1, velocity_mesh, name = "NonlinearVelocity")
  call allocate(cfl_no, velocity_mesh, name = "CFLNumber")
  
  call deallocate(coordinate_mesh)
  call deallocate(velocity_mesh)
    
  call set(positions, (/1, 2, 3, 4/), spread((/0.0, 1.0, 11.0, 111.0/), 1, 1))
  call set(velocity, (/1, 2, 3/), spread((/1.0, 1.0, 5.0/), 1, 1))
  call set(positions_3d, 1, (/0.0, 0.0, 0.0/))
  call set(positions_3d, 2, (/1.0, 0.0, 0.0/))
  call set(positions_3d, 3, (/11.0, 0.0, 0.0/))
  call set(positions_3d, 4, (/111.0, 0.0, 0.0/))
  call set(velocity_3d, 1, (/1.0, 0.0, 0.0/))
  call set(velocity_3d, 2, (/1.0, 0.0, 0.0/))
  call set(velocity_3d, 3, (/ 5.0, 0.0,0.0/))
  call project_cartesian_to_local(positions_3d, velocity_3d, local_velocity)
  
  call insert(state, positions, name = positions%name)
  call insert(state, velocity, name = velocity%name)
  call insert(local_state, positions_3d, name = positions_3d%name)
  call insert(local_state, velocity_3d, name = velocity_3d%name)
  call insert(local_state, local_velocity, name = local_velocity%name)
  call deallocate(positions)
  call deallocate(velocity)
  call deallocate(positions_3d)
  call deallocate(velocity_3d)
  call deallocate(local_velocity)
  
  call calculate_cfl_number(state, cfl_no, dt = 1.0)
  call report_test("[cfl no]", node_val(cfl_no, (/1, 2, 3/)) .fne. (/1.0, 0.1, 0.05/), .false., "Incorrect CFL number")
  
  call calculate_cfl_number(state, cfl_no, dt = 10.0)
  call report_test("[cfl no]", node_val(cfl_no, (/1, 2, 3/)) .fne. (/10.0, 1.0, 0.5/), .false., "Incorrect CFL number")
  
  call calculate_courant_number_dg(state, cfl_no, dt=1.0)
  print*, node_val(cfl_no, (/1, 2, 3/))

  call calculate_courant_number_dg(local_state, cfl_no, dt=1.0)
  print*, node_val(cfl_no, (/1, 2, 3/))

  call deallocate(state)
  call deallocate(local_state)
  call deallocate(cfl_no)
  
  call report_test_no_references()

end subroutine test_cfl_number_1d
