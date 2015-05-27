!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineeringp
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

subroutine project_vtu(input_filename_, input_filename_len, donor_basename_, donor_basename_len, &
                       &target_basename_, target_basename_len, output_filename_, output_filename_len) bind(c)

  use conservative_interpolation_module
  use fields
  use fldebug
  use global_parameters, only : current_debug_level
  use intersection_finder_module
  use linked_lists
  use solvers
  use spud
  use state_module
  use mesh_files
  use vtk_interfaces
  use iso_c_binding
  implicit none
  
  integer(kind=c_size_t), value :: input_filename_len, donor_basename_len
  integer(kind=c_size_t), value :: target_basename_len, output_filename_len
  character(kind=c_char, len=1) :: input_filename_(*), donor_basename_(*)
  character(kind=c_char, len=1) :: target_basename_(*), output_filename_(*)
  
  character(len = input_filename_len) :: input_filename
  character(len = donor_basename_len) :: donor_basename
  character(len = target_basename_len) :: target_basename
  character(len = output_filename_len) :: output_filename
  character(len = *), parameter :: fields_path = "/dummy"
  integer :: i
  integer, parameter :: quad_degree = 4
  type(element_type), pointer :: shape
  type(ilist), dimension(:), allocatable :: map_BA
  type(mesh_type) :: output_mesh, output_p0mesh
  type(mesh_type), pointer :: input_mesh
  type(state_type) :: input_state, output_state
  type(state_type), dimension(:), allocatable :: input_mesh_states, output_mesh_states
  type(scalar_field) :: output_s_field
  type(scalar_field), pointer :: input_s_field
  type(vector_field) :: donor_positions, output_v_field, target_positions
  type(vector_field), pointer :: input_v_field
  type(tensor_field) :: output_t_field
  type(tensor_field), pointer :: input_t_field
  
  ewrite(1, *) "In project_vtu"

  do i=1, input_filename_len
    input_filename(i:i)=input_filename_(i)
  end do
  do i=1, donor_basename_len
    donor_basename(i:i)=donor_basename_(i)
  end do
  do i=1, target_basename_len
    target_basename(i:i)=target_basename_(i)
  end do
  do i=1, output_filename_len
    output_filename(i:i)=output_filename_(i)
  end do

  
  call set_solver_options(fields_path // "/galerkin_projection/continuous", &
    & ksptype = "cg", pctype = "sor", atol = epsilon(0.0), rtol = 0.0, max_its = 2000, start_from_zero = .true.)
      
  call vtk_read_state(trim(input_filename), input_state, quad_degree = quad_degree)
  
  donor_positions = extract_vector_field(input_state, "Coordinate")
  input_mesh => extract_mesh(input_state, "Mesh")

  target_positions = read_mesh_files(trim(target_basename), quad_degree = quad_degree, format="gmsh")
  
  shape => ele_shape(donor_positions, 1)
  if (shape==ele_shape(target_positions,1) .and. continuity(donor_positions)==continuity(target_positions)) then
    output_mesh = target_positions%mesh
    call incref(output_mesh)
  else
    output_mesh = make_mesh(target_positions%mesh, shape, continuity = continuity(donor_positions))
  end if
  if (has_mesh(input_state, "P0Mesh")) then
    output_p0mesh = piecewise_constant_mesh(target_positions%mesh, "P0Mesh")
  end if
  
  if (donor_basename_len>0) then
    donor_positions = read_mesh_files(trim(donor_basename), quad_degree = quad_degree, format="gmsh")
  else
    ! no donor mesh specified:
    ! use the one we got from vtk_read_state
    ! this only works in serial and with a P1CG vtu
    if (continuity(donor_positions)<0 .or. shape%degree/=1) then
      FLExit("No donor mesh specified. This only works for a serial, continuous, linear input vtu")
    end if
    call incref(donor_positions)
  end if
  
  allocate(map_BA(ele_count(target_positions)))
  map_BA = rtree_intersection_finder(target_positions, donor_positions)
  
  call insert(output_state, output_mesh, "Mesh")
  if (has_mesh(input_state, "P0Mesh")) then
    call insert(output_state, output_p0mesh, "P0Mesh")
  end if
  do i = 1, scalar_field_count(input_state)
    input_s_field => extract_scalar_field(input_state, i)
    if(input_s_field%name == "vtkGhostLevels") then
      call remove_scalar_field(input_state, input_s_field%name)
      cycle
    end if
    if (input_s_field%mesh%name=="Mesh") then
      call allocate(output_s_field, output_mesh, input_s_field%name)
    else if (input_s_field%mesh%name=="P0Mesh") then
      call allocate(output_s_field, output_p0mesh, input_s_field%name)
    else
      FLAbort("State from vtk_read_state should contain Mesh and P0Mesh only")
    end if
    call zero(output_s_field)
    output_s_field%option_path = fields_path
    call insert(output_state, output_s_field, output_s_field%name)
    call deallocate(output_s_field)
  end do
  do i = 1, vector_field_count(input_state)
    input_v_field => extract_vector_field(input_state, i)
    if(input_v_field%name == "Coordinate") cycle
    if (input_v_field%mesh%name=="Mesh") then
      call allocate(output_v_field, input_v_field%dim, output_mesh, input_v_field%name)
    else if (input_v_field%mesh%name=="P0Mesh") then
      call allocate(output_v_field, input_v_field%dim, output_p0mesh, input_v_field%name)
    else
      FLAbort("State from vtk_read_state should contain Mesh and P0Mesh only")
    end if
    call zero(output_v_field)
    output_v_field%option_path = fields_path
    call insert(output_state, output_v_field, output_v_field%name)
    call deallocate(output_v_field)
  end do
  do i = 1, tensor_field_count(input_state)
    input_t_field => extract_tensor_field(input_state, i)
    if (input_t_field%mesh%name=="Mesh") then
      call allocate(output_t_field, output_mesh, input_t_field%name)
    else if (input_t_field%mesh%name=="P0Mesh") then
      call allocate(output_t_field, output_p0mesh, input_t_field%name)
    else
      FLAbort("State from vtk_read_state should contain Mesh and P0Mesh only")
    end if
    call zero(output_t_field)
    output_t_field%option_path = fields_path
    call insert(output_state, output_t_field, output_t_field%name)
    call deallocate(output_t_field)
  end do
  
  
  if(current_debug_level >= 2) then
    ewrite(2, *) "Options tree:"
    call print_options()
    ewrite(2, *) "Input state:"
    call print_state(input_state)
    ewrite(2, *) "Output state:"
    call print_state(output_state)
  end if

  call sort_states_by_mesh( (/ input_state /), input_mesh_states)
  call sort_states_by_mesh( (/ output_state /), output_mesh_states)
  do i=1, size(input_mesh_states)
    call insert(input_mesh_states(i), donor_positions, "Coordinate")
  end do
  do i=1, size(output_mesh_states)
    call insert(output_mesh_states(i), target_positions, "Coordinate")
  end do
  ! do this insert after the sort_states_by_mesh as target_positions%mesh is seen as a different mesh
  call insert(output_state, target_positions, "Coordinate")
  call deallocate(donor_positions)
  call deallocate(target_positions)
  
  call interpolation_galerkin(input_mesh_states, output_mesh_states, map_BA = map_BA)
  call deallocate(map_BA)
  deallocate(map_BA)
  call deallocate(input_mesh_states)
  deallocate(input_mesh_states)
  call deallocate(output_mesh_states)
  deallocate(output_mesh_states)
  call deallocate(input_state)
  
  call vtk_write_state(trim(output_filename), model = "Mesh", state = (/output_state/))
  call deallocate(output_mesh)
  if (has_mesh(output_state, "P0Mesh")) then
    call deallocate(output_p0mesh)
  end if
  call deallocate(output_state)
  
  call print_references(0)
  
  ewrite(1, *) "Exiting project_vtu"
       
end subroutine project_vtu
