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

subroutine test_geostrophic_pressure

  use field_options
  use elements
  use fields
  use fldebug
  use geostrophic_pressure
  use solvers
  use spud
  use state_module
  use unittest_tools
  use vtk_interfaces

  implicit none
  
  integer :: i, stat
  logical :: fail
  real :: gp_zero_point, y_zero_point, z_zero_point
  real, dimension(:), allocatable :: pos
  type(element_type) :: gp_shape
  type(mesh_type) :: gp_mesh
  type(mesh_type), pointer :: mesh
  type(scalar_field) :: gp, buoyancy
  type(state_type) :: state
  type(vector_field) :: gravity_direction, positions_remap, velocity
  type(vector_field), pointer :: positions
    
  call vtk_read_state("data/pseudo2d.vtu", state)
  
  positions => extract_vector_field(state, "Coordinate")
  mesh => positions%mesh
  ! 3D test
  assert(mesh_dim(mesh) == 3)
  call set_option("/geometry/dimension", 3, stat = stat)
  assert(stat == SPUD_NEW_KEY_WARNING)
  
  allocate(pos(positions%dim))
  
  gp_shape = make_element_shape(ele_shape(mesh, 1), degree = 2)
  gp_mesh = make_mesh(mesh, shape = gp_shape)
  call deallocate(gp_shape)
  call allocate(gp, gp_mesh, gp_name)
  call deallocate(gp_mesh)
  call zero(gp)
  call set_solver_options(gp, ksptype = "cg", pctype = "mg", atol = epsilon(0.0), rtol = 0.0, max_its = 2000, start_from_zero = .false.)
  call insert(state, gp, gp%name)
  
  ! Hydrostatic balance
  call allocate(buoyancy, mesh, "VelocityBuoyancyDensity", field_type = FIELD_TYPE_CONSTANT)
  call set(buoyancy, 1.0)
  call insert(state, buoyancy, buoyancy%name)
  call deallocate(buoyancy)
  
  call set_option("/physical_parameters/gravity/magnitude", 1.0, stat = stat)
  assert(stat == SPUD_NEW_KEY_WARNING)
  
  call allocate(gravity_direction, mesh_dim(mesh), mesh, "GravityDirection", field_type = FIELD_TYPE_CONSTANT)
  call set(gravity_direction, (/0.0, 0.0, 1.0/))
  call insert(state, gravity_direction, gravity_direction%name)
  call deallocate(gravity_direction)
  
  call calculate_geostrophic_pressure(state, gp, &
    & assemble_matrix = .true., include_buoyancy = .true., include_coriolis = .false., reference_node = 1)
      
  call allocate(positions_remap, positions%dim, gp%mesh, "Positions")
  call remap_field(positions, positions_remap)
  
  fail = .false.
  gp_zero_point = minval(gp%val)
  z_zero_point = minval(positions_remap%val(3,:))
  do i = 1, node_count(gp)
    pos = node_val(positions_remap, i)
    if(fnequals(node_val(gp, i) - gp_zero_point, pos(3) - z_zero_point, tol = 1.0e-10)) then
      ewrite(-1, *) "Error: ", abs((node_val(gp, i) - gp_zero_point) - (pos(3) - z_zero_point))
      fail = .true.
      exit
    end if
  end do
  call report_test("[Hydrostatic balance]", fail, .false., "Incorrect solution")
  
  ! Geostrophic balance
  call set_option("/physical_parameters/coriolis/f_plane/f", 1.0, stat)
  assert(stat == SPUD_NEW_KEY_WARNING)
  
  call allocate(velocity, mesh_dim(mesh), mesh, "Velocity", field_type = FIELD_TYPE_CONSTANT)
  call set(velocity, (/-1.0, 0.0, 0.0/))
  call insert(state, velocity, velocity%name)
  call insert(state, velocity%mesh, trim(velocity%name) // "Mesh")
  call deallocate(velocity)
    
  call calculate_geostrophic_pressure(state, gp, &
    & velocity_name = velocity%name, assemble_matrix = .false., include_buoyancy = .false., include_coriolis = .true.)
    
  fail = .false.
  gp_zero_point = minval(gp%val)
  y_zero_point = minval(positions_remap%val(2,:))
  do i = 1, node_count(gp)
    pos = node_val(positions_remap, i)
    if(fnequals(node_val(gp, i) - gp_zero_point, pos(2) - y_zero_point, tol = 1.0e-10)) then
      ewrite(-1, *) "Error: ", abs((node_val(gp, i) - gp_zero_point) - (pos(2) - y_zero_point))
      fail = .true.
      exit
    end if
  end do
  call report_test("[Geostrophic balance]", fail, .false., "Incorrect solution")
    
  ! Thermal wind balance
  call calculate_geostrophic_pressure(state, gp, &
    & velocity_name = velocity%name, assemble_matrix = .false., include_buoyancy = .true., include_coriolis = .true.)
  
  fail = .false.
  gp_zero_point = minval(gp%val)
  y_zero_point = minval(positions_remap%val(2,:))
  do i = 1, node_count(gp)
    pos = node_val(positions_remap, i)
    if(fnequals(node_val(gp, i) - gp_zero_point, (pos(2) - y_zero_point) + (pos(3) - z_zero_point), tol = 1.0e-10)) then
      ewrite(-1, *) "Error: ", abs((node_val(gp, i) - gp_zero_point) - ((pos(2) - y_zero_point) + (pos(3) - z_zero_point)))
      fail = .true.
      exit
    end if
  end do
  call report_test("[Thermal wind balance]", fail, .false., "Incorrect solution")
  
  call deallocate(positions_remap)
  
  call deallocate(gp)
  call deallocate(state)
  deallocate(pos)
  
  call report_test_no_references()
  
end subroutine test_geostrophic_pressure
