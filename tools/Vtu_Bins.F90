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

subroutine vtu_bins(input_filename, input_filename_len, input_fieldname, input_fieldname_len, bounds, nbounds)

  use fields
  use fldebug
  use futils
  use mixing_statistics
  use quicksort
  use reference_counting
  use state_module
  use vtk_interfaces

  implicit none
  
  integer, intent(in) :: input_filename_len
  integer, intent(in) :: input_fieldname_len
  integer, intent(in) :: nbounds
  
  character(len = input_filename_len), intent(in) :: input_filename
  character(len = input_fieldname_len), intent(in) :: input_fieldname
  real, dimension(nbounds), intent(in) :: bounds
  
  character(len = real_format_len()) :: rformat
  character(len = real_format_len(padding = 1)) :: rformatp
  integer :: i
  integer, dimension(nbounds) :: permutation
  logical :: allocated
  real :: volume
  real, dimension(nbounds) :: lbounds
  real, dimension(nbounds + 1) :: integrals
  type(scalar_field), pointer :: field
  type(state_type) :: state
  type(vector_field), pointer :: positions
  
  ewrite(1, *) "In vtu_bins"
  
  ewrite(2, *) "Input file: ", trim(input_filename)
  ewrite(2, *) "Input field: ", trim(input_fieldname)
  ewrite(2, *) "Bounds: ", bounds
  
  call qsort(bounds, permutation)
  lbounds = bounds(permutation)
  
  call vtk_read_state(input_filename, state = state)
  positions => extract_vector_field(state, "Coordinate")
  field => extract_scalar_field(state, input_fieldname, allocated = allocated)
  
  volume = 0.0
  do i = 1, ele_count(positions)
    volume = volume + element_volume(positions, i)
  end do
  
  integrals = heaviside_integral(field, lbounds, positions)
  integrals(nbounds + 1) = heaviside_integral(field, huge(0.0), positions)
  
  rformat = real_format()
  rformatp = real_format(padding = 1)
  print "(a," // rformatp // ",a," // rformat // ")", "                      -inf - ", lbounds(1), " : ", (volume - integrals(1)) / volume
  do i = 1, nbounds - 1
    print "(" // rformatp // ",a," // rformatp // ",a," // rformat // ")", lbounds(i), " - ", lbounds(i + 1), " : ", (-integrals(i + 1) + integrals(i)) / volume
  end do
  print "(" // rformatp // "a," // rformat // ")", lbounds(nbounds), " -                        inf : ", (-integrals(nbounds + 1) + integrals(nbounds)) / volume
  
  if(allocated) deallocate(field)
  call deallocate(state)
  
  call print_references(0)
  
  ewrite(1, *) "Exiting vtu_bins"
  
end subroutine vtu_bins
