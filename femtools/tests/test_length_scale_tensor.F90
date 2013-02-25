!    Copyright (C) 2009 Imperial College London and others.
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

subroutine test_length_scale_tensor

   use unittest_tools
   use read_triangle
   use fields
   use smoothing_module
   implicit none
   
   type(vector_field) :: positions
   type(element_type) :: shape
   real, dimension(:,:,:), allocatable :: dshape
   real, dimension(:), allocatable :: detwei
   real, dimension(:,:,:), allocatable :: expected_result, computed_result
   logical :: fail
   
   positions = read_triangle_files("data/structured", quad_degree=3)

   shape = ele_shape(positions,1)
   allocate(computed_result(positions%dim, positions%dim, ele_ngi(positions, 1)))
   allocate(expected_result(positions%dim, positions%dim, ele_ngi(positions, 1)))
   allocate(detwei(ele_ngi(positions, 1)))
   allocate(dshape(ele_loc(positions, 1), ele_ngi(positions, 1), positions%dim))

   call transform_to_physical(positions, 1, shape, dshape=dshape, detwei=detwei)
   ewrite(0,*) "dshape ", dshape
   ! We'll just choose the first element here - each of them should have the same area
   computed_result = length_scale_tensor(dshape, shape)
   expected_result = 0.0
   ewrite(0,*) "tensor ", computed_result(:,:,1)
   fail = .not.fequals(computed_result(1,1,1), expected_result(1,1,1), 1.0e-9)
   call report_test("[length_scale_tensor]", fail, .false., "Result from length_scale_tensor is incorrect.")
   
   deallocate(computed_result)
   deallocate(expected_result)
   deallocate(detwei)
   deallocate(dshape)

end subroutine test_length_scale_tensor

