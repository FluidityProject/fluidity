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
use sparsity_patterns_meshes
use element_numbering
use state_module
use fldebug_parameters
use ieee_arithmetic, only: ieee_quiet_nan, ieee_value
use spud
implicit none

private
public nodalise_bubble_basis
!public get_p2b_mass, solve_p2b_lumped_mass

!subroutine get_p2b_mass(state,p2b_mass,p2b_lumped_mass)
!end subroutine get_p2b_mass

contains

  subroutine setup_P2b_P1dg_projection(state,P2b_mesh,P1dg_mesh)
    !Check if various matrices for P2b_P1dg projection are present,
    !if not, create them.
    type(state_type), intent(inout) :: state
    type(mesh_type), intent(inout) :: P2b_mesh, P1dg_mesh
    !
    integer :: stat, ele
    type(csr_matrix) :: p2b_p1dg_projection, p2b_mass, p1dg_mass
    type(scalar_field) :: p2b_lumped_mass, p1dg_lumped_mass
    type(csr_matrix), pointer :: matrix
    type(csr_sparsity), pointer :: sparsity
    type(scalar_field), pointer :: s_field

    !! P2b mass matrix
    matrix=>extract_csr_matrix(state,trim(P2b_mesh%name)//"Mass",stat)
    if(stat.ne.0) then
       sparsity => get_csr_sparsity_firstorder(&
            state, p2b_mesh, p2b_mesh)
       call allocate(p2b_mass,sparsity,name=trim(P2b_mesh%name)//"Mass")
       do ele = 1, element_count(P2b_mesh)
          call get_p2b_mass_ele(ele)
       end do
       call insert(state,p2b_mass,p2b_mass%name)
       call deallocate(p2b_mass)
    end if

    !! Lumped P2b mass matrix (stored as scalar field)
    s_field=>extract_scalar_field(state,trim(P2b_mesh%name)//"LumpedMass",stat)
    if(stat.ne.0) then
       call allocate(p2b_lumped_mass,P2b_mesh,&
            name=trim(P2b_mesh%name)//"LumpedMass")
       do ele = 1, element_count(P2b_mesh)
          call get_p2b_lumped_mass_ele(ele)
       end do
       call insert(state,p2b_lumped_mass,p2b_lumped_mass%name)
       call deallocate(p2b_lumped_mass)
    end if

    !! P1dg mass matrix
    matrix=>extract_csr_matrix(state,trim(P1dg_mesh%name)//"Mass",stat)
    if(stat.ne.0) then
       sparsity => get_csr_sparsity_firstorder(&
            state, p1dg_mesh, p1dg_mesh)
       call allocate(p1dg_mass,sparsity,name=trim(P1dg_mesh%name)//"Mass")
       do ele = 1, element_count(P1dg_mesh)
          call get_p1dg_mass_ele(ele)
       end do
       call insert(state,p1dg_mass,p1dg_mass%name)
       call deallocate(p1dg_mass)
    end if

    !! Lumped P1dg mass matrix (stored as scalar field)
    s_field=>extract_scalar_field(state,trim(P1dg_mesh%name)//"LumpedMass",stat)
    if(stat.ne.0) then
       call allocate(p1dg_lumped_mass,P1dg_mesh,&
            name=trim(P1dg_mesh%name)//"LumpedMass")
       do ele = 1, element_count(P1dg_mesh)
          call get_p1dg_lumped_mass_ele(ele)
       end do
       call insert(state,p1dg_lumped_mass,p1dg_lumped_mass%name)
       call deallocate(p1dg_lumped_mass)
    end if

    !!P1dg to P2b projection matrix
    matrix=>extract_csr_matrix(state,trim(P2b_mesh%name)//&
         &trim(P1dg_mesh%name)//"Projection",stat)
    if(stat.ne.0) then
       sparsity => get_csr_sparsity_firstorder(&
            state, P2b_mesh, P1dg_mesh)
       call allocate(p2b_p1dg_projection,sparsity,name=trim(P2b_mesh%name)//&
            &trim(P1dg_mesh%name)//"Projection")
       do ele = 1, element_count(P2b_mesh)
          call get_p2b_p1dg_projection_ele(ele)
       end do
       call insert(state,p2b_p1dg_projection,p2b_p1dg_projection%name)
       call deallocate(p2b_p1dg_projection)
    end if

  contains 
    subroutine get_p2b_mass_ele(ele)
      integer, intent(in) :: ele
      !
    end subroutine get_p2b_mass_ele

    subroutine get_p2b_lumped_mass_ele(ele)
      integer, intent(in) :: ele
      !
    end subroutine get_p2b_lumped_mass_ele

    subroutine get_p1dg_mass_ele(ele)
      integer, intent(in) :: ele
      !
    end subroutine get_p1dg_mass_ele

    subroutine get_p1dg_lumped_mass_ele(ele)
      integer, intent(in) :: ele
      !
    end subroutine get_p1dg_lumped_mass_ele

    subroutine get_p2b_p1dg_projection_ele(ele)
      integer, intent(in) :: ele
      !
    end subroutine get_p2b_p1dg_projection_ele

  end subroutine setup_P2b_P1dg_projection

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
