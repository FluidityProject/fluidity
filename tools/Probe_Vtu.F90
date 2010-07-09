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

subroutine probe_vtu(vtu_filename, vtu_filename_len, fieldname, fieldname_len, x, y, z, dim)

  use fields
  use fldebug
  use futils
  use pickers
  use reference_counting
  use state_module
  use vtk_interfaces
  
  implicit none
  
  integer, intent(in) :: vtu_filename_len
  integer, intent(in) :: fieldname_len
  
  character(len = vtu_filename_len), intent(in) :: vtu_filename
  character(len = fieldname_len), intent(in) :: fieldname
  real, intent(in) :: x
  real, intent(in) :: y
  real, intent(in) :: z
  integer, intent(in) :: dim
  
  character(len = 1 + 1 + real_format_len(padding = 1) + 1) :: format
  integer :: ele, stat
  logical :: allocated
  real :: s_val
  real, dimension(dim) :: coord, v_val
  real, dimension(dim + 1) :: local_coord
  real, dimension(dim, dim) :: t_val
  type(scalar_field), pointer :: s_field
  type(vector_field), pointer :: positions, v_field
  type(state_type) :: state
  type(tensor_field), pointer :: t_field
  
  ewrite(1, *) "In probe_vtu"
  
  call vtk_read_state(vtu_filename, state = state)
  
  positions => extract_vector_field(state, "Coordinate")
  if(positions%dim /= dim) then
    FLExit("Expected " // int2str(dim) // " dimensional probe coord")
  else if(ele_numbering_family(ele_shape(positions, 1)) /= FAMILY_SIMPLEX) then
    FLExit("Mesh in vtu " // vtu_filename // " is not composed of linear simplices")
  end if
  
  if(dim > 0) coord(1) = x
  if(dim > 1) coord(2) = y
  if(dim > 2) coord(3) = z
  ewrite(2, *) "Probe coord: ", coord
  call picker_inquire(positions, coord, ele, local_coord = local_coord)
  if(ele < 0) then
    FLExit("Probe point not contained in vtu " // vtu_filename)
  end if
  
  s_field => extract_scalar_field(state, fieldname, allocated = allocated, stat = stat)
  if(stat == 0) then
    format = "(" // real_format() // ")"
    s_val = eval_field(ele, s_field, local_coord)
    ewrite(2, *) fieldname // " value:"
    print format, s_val
    if(allocated) deallocate(s_field)
  else
    v_field => extract_vector_field(state, fieldname, stat = stat)
    if(stat == 0) then
      format = "(" // int2str(dim) // real_format(padding = 1) // ")"
      v_val = eval_field(ele, v_field, local_coord)
      ewrite(2, *) fieldname // " value:"
      print format, v_val
    else
      t_field => extract_tensor_field(state, fieldname, stat = stat)
      if(stat == 0) then
        format = "(" // int2str(dim ** 2) // real_format(padding = 1) // ")"
        t_val = eval_field(ele, t_field, local_coord)
        ewrite(2, *) fieldname // " value:"
        print format, t_val
      else
        FLExit("Field " // fieldname // " not found in vtu " // vtu_filename)
      end if
    end if
  end if
  
  call deallocate(state)
  
  call print_references(0)
  
  ewrite(1, *) "Exiting probe_vtu"

end subroutine probe_vtu
