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

subroutine project_vtu(input_filename, input_filename_len, donor_basename, donor_basename_len, target_basename, target_basename_len, output_filename, output_filename_len)

  use conservative_interpolation_module
  use fields
  use fldebug
  use global_parameters, only : current_debug_level
  use intersection_finder_module
  use linked_lists
  use solvers
  use spud
  use state_module
  use read_triangle
  use vtk_interfaces
  use write_triangle
  
  implicit none
  
  integer, intent(in) :: input_filename_len
  integer, intent(in) :: donor_basename_len
  integer, intent(in) :: target_basename_len
  integer, intent(in) :: output_filename_len
  character(len = input_filename_len), intent(in) :: input_filename
  character(len = donor_basename_len), intent(in) :: donor_basename
  character(len = target_basename_len), intent(in) :: target_basename
  character(len = output_filename_len), intent(in) :: output_filename
  
  character(len = *), parameter :: fields_path = "/dummy"
  integer :: i
  integer, parameter :: quad_degree = 4
  type(element_type), pointer :: shape
  type(ilist), dimension(:), allocatable :: map_BA
  type(mesh_type) :: output_mesh
  type(mesh_type), pointer :: input_mesh
  type(state_type) :: input_state, output_state
  type(scalar_field) :: output_s_field
  type(scalar_field), pointer :: input_s_field
  type(vector_field) :: donor_positions, output_v_field, target_positions
  type(vector_field), pointer :: input_v_field
  type(tensor_field) :: output_t_field
  type(tensor_field), pointer :: input_t_field
  
  ewrite(1, *) "In project_vtu"
  
  call set_solver_options(fields_path // "/galerkin_projection/continuous", &
    & ksptype = "cg", pctype = "sor", atol = epsilon(0.0), rtol = 0.0, max_its = 2000, start_from_zero = .true.)
      
  call vtk_read_state(trim(input_filename), input_state, quad_degree = quad_degree)
  
  donor_positions = extract_vector_field(input_state, "Coordinate")
  input_mesh => extract_mesh(input_state, "Mesh")
  target_positions = read_triangle_files(trim(target_basename), quad_degree = quad_degree)
  shape => ele_shape(donor_positions, 1)
  output_mesh = make_mesh(target_positions%mesh, shape, continuity = continuity(donor_positions))
  
  donor_positions = read_triangle_files(trim(donor_basename), quad_degree = quad_degree)
  
  allocate(map_BA(ele_count(target_positions)))
  map_BA = rtree_intersection_finder(target_positions, donor_positions)
  
  call insert(input_state, donor_positions, "Coordinate")
  
  call insert(output_state, target_positions, "Coordinate")
  call insert(output_state, output_mesh, "Mesh")
  do i = 1, scalar_field_count(input_state)
    input_s_field => extract_scalar_field(input_state, i)
    if(input_s_field%name == "vtkGhostLevels") then
      call remove_scalar_field(input_state, input_s_field%name)
      cycle
    end if
    call allocate(output_s_field, output_mesh, input_s_field%name)
    call zero(output_s_field)
    output_s_field%option_path = fields_path
    call insert(output_state, output_s_field, output_s_field%name)
    call deallocate(output_s_field)
  end do
  do i = 1, vector_field_count(input_state)
    input_v_field => extract_vector_field(input_state, i)
    if(input_v_field%name == "Coordinate") cycle
    call allocate(output_v_field, input_v_field%dim, output_mesh, input_v_field%name)
    call zero(output_v_field)
    output_v_field%option_path = fields_path
    call insert(output_state, output_v_field, output_v_field%name)
    call deallocate(output_v_field)
  end do
  do i = 1, tensor_field_count(input_state)
    input_t_field => extract_tensor_field(input_state, i)
    call allocate(output_t_field, output_mesh, input_t_field%name)
    call zero(output_t_field)
    output_t_field%option_path = fields_path
    call insert(output_state, output_t_field, output_t_field%name)
    call deallocate(output_t_field)
  end do
  
  call deallocate(donor_positions)
  call deallocate(target_positions)
  
  if(current_debug_level >= 2) then
    ewrite(2, *) "Options tree:"
    call print_options()
    ewrite(2, *) "Input state:"
    call print_state(input_state)
    ewrite(2, *) "Output state:"
    call print_state(output_state)
  end if
  
  call interpolation_galerkin(input_state, output_state, map_BA = map_BA)
  call deallocate(map_BA)
  deallocate(map_BA)
  call deallocate(input_state)
  
  call vtk_write_state(trim(output_filename), model = "Mesh", state = (/output_state/))
  call deallocate(output_state)
  call deallocate(output_mesh)
  
  call print_references(0)
  
  ewrite(1, *) "Exiting project_vtu"
       
end subroutine project_vtu
