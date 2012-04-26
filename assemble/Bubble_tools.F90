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
  
  module bubble_tools
    use fields
    use vector_tools
    use sparse_matrices_fields
    use sparsity_patterns_meshes
    use element_numbering
    use state_module
    use fldebug_parameters
    use ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use spud
    implicit none

    private
    public nodalise_bubble_basis
  contains

    subroutine nodalise_bubble_basis(shape)
      !Subroutine to transform bubble basis to a equivalent nodal one.
      type(element_type), intent(inout) :: shape
      !
      if(shape%numbering%type .ne. ELEMENT_BUBBLE) then
         FLAbort('Only applies to bubbles.')
      end if

      select case(shape%dim)
      case (2)
         select case (shape%numbering%vertices)
         case (3)
            select case( shape%loc )
            case (7)
               call nodalise_bubble_basis_P2b()
            case default
               FLAbort('Element not supported.')
            end select
         case default
            FLAbort('Family not supported')
         end select
      case default
         FLAbort('Dimension not supported.')
      end select

    contains 

      subroutine nodalise_bubble_basis_P2b()
        !         
        integer :: i,j,loc
        real, dimension(7) :: N_vals
        N_vals = eval_shape(shape, (/1.0/3.0,1.0/3.0,1.0/3.0/))

        !! n is loc x ngi, dn is loc x ngi x dim         
        shape%n(7,:) = shape%n(7,:)/N_vals(7)
        do loc = 1, 6
           shape%n(loc,:) = shape%n(loc,:) - N_vals(loc)*shape%n(7,:)
           shape%dn(loc,:,:) = shape%dn(loc,:,:) - &
                &N_vals(loc)*shape%dn(7,:,:)
        end do

        !spoly is now useless
        if(associated(shape%spoly)) then
           do i=1,size(shape%spoly,1)
              do j=1,size(shape%spoly,2)
                 shape%spoly(i,j) = (/ieee_value(0.0,ieee_quiet_nan)/)
              end do
           end do
        end if

      end subroutine nodalise_bubble_basis_P2b

    end subroutine nodalise_bubble_basis

  end module bubble_tools
