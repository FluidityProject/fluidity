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

subroutine test_tensor_inner_product

   use unittest_tools
   use k_epsilon
   implicit none
   
   interface
     function solution()
        real, dimension(1) :: solution
     end function
   end interface
   
   ! The two input matrices, A and B.
   ! This unit test will take the inner product of these.
   real, dimension(3, 3, 1) :: A, B
   real, dimension(1) :: expected_result, computed_result
   logical :: fail
   
   ! Note: the transpose operation gets everything into row-major order here.
   A(:,:,1) = transpose(reshape((/ 1, 2, 3, 4, 5, 6, 7, 8, 9 /), shape(A(:,:,1))))
   B(:,:,1) = transpose(reshape((/ 2, 4, 6, 8, 10, 12, 14, 16, 18 /), shape(B(:,:,1))))
   
   computed_result = tensor_inner_product(A, B)
   expected_result = solution()
   
   fail = (computed_result(1) /= expected_result(1))
   call report_test("[tensor_inner_product]", fail, .false., "Result from tensor_inner_product is incorrect.")
   
end subroutine test_tensor_inner_product

function solution()
  real, dimension(1) :: solution
  solution = (/570.0/)
end function solution
