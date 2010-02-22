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

subroutine derive_mesh(input_filename, input_filename_len, output_filename, output_filename_len, degree, cont, vtu)

  use fields
  use fldebug
  use halos
  use parallel_tools
  use read_triangle
  use shape_functions
  use vtk_interfaces
  use write_triangle
  
  implicit none
  
  integer, intent(in) :: input_filename_len
  integer, intent(in) :: output_filename_len
  character(len = input_filename_len), intent(in) :: input_filename
  character(len = output_filename_len), intent(in) :: output_filename
  integer, intent(in) :: degree
  integer, intent(in) :: cont
  integer, intent(in) :: vtu
  
  type(element_type) :: derived_shape
  type(element_type), pointer :: base_shape
  type(mesh_type) :: derived_mesh
  type(mesh_type), pointer :: base_mesh
  type(vector_field) :: derived_positions
  type(vector_field), target :: base_positions
  
  ewrite(1, *) "In derive_mesh"
  
  ewrite(2, *) "Input file: " // trim(input_filename)
  ewrite(2, *) "Output file: " // trim(output_filename)
  ewrite(2, *) "Mesh degree: ", degree
  ewrite(2, *) "Mesh is continuous? ", cont == 0
  
  base_positions = read_triangle_files(trim(input_filename), quad_degree = 1)
  base_mesh => base_positions%mesh
  if(isparallel()) call read_halos(trim(input_filename), base_mesh)
  
  base_shape => ele_shape(base_mesh, 1)
  derived_shape = make_element_shape(base_shape, degree = degree)
  derived_mesh = make_mesh(base_mesh, derived_shape, continuity = cont)
  call deallocate(derived_shape)
  
  call allocate(derived_positions, base_positions%dim, derived_mesh, name = base_positions%name)
  call remap_field(base_positions, derived_positions)
  call deallocate(derived_mesh)
  call deallocate(base_positions)
  
  if(vtu == 0) then
    call write_triangle_files(parallel_filename(output_filename), derived_positions)
  else
    if(isparallel()) then
      call vtk_write_fields(trim(output_filename) // ".pvtu", position = derived_positions, model = derived_positions%mesh)
    else
      call vtk_write_fields(trim(output_filename) // ".vtu", position = derived_positions, model = derived_positions%mesh)
    end if
  end if
  call deallocate(derived_positions)
  
  call print_references(0)
  
  ewrite(1, *) "Exiting derive_mesh"
       
end subroutine derive_mesh
