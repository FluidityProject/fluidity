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
    use solvers
    use sparsity_patterns_meshes
    use element_numbering
    use state_module
    use fldebug_parameters
    use ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use spud
    implicit none

    private
    public nodalise_bubble_basis, get_lumped_mass_p2b, project_to_p2b_lumped
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

    !! This is special cased because of the special properties of the 
    !! p2b lumped mass (it is still 3rd order, and positive).
    subroutine get_lumped_mass_p2b(state,mass,p2b_field,weight_field)
      type(csr_matrix), intent(inout) :: mass
      type(state_type), intent(inout) :: state
      type(scalar_field), intent(in) :: p2b_field
      type(scalar_field), intent(in), optional :: weight_field
      !
      integer :: ele
      type(vector_field), pointer :: X
      real, dimension(7) :: N_vals

      X=>extract_vector_field(state, "Coordinate")
      
      N_vals = eval_shape(p2b_field%mesh%shape, (/1.0/3.0,1.0/3.0,1.0/3.0/))

      if(p2b_field%mesh%shape%dim.ne.2) then
         FLAbort('Only works for 2d meshes')
      end if
      if(p2b_field%mesh%shape%loc.ne.7) then
         FLAbort('Expected p2 bubble mesh')
      end if

      call zero(mass)

      do ele = 1, ele_count(X)
         call get_lumped_mass_p2b_ele(mass,p2b_field,X,N_vals,ele,weight_field)
      end do

    end subroutine get_lumped_mass_p2b

    subroutine get_lumped_mass_p2b_ele(mass,p2b_field,X,N_vals,ele,weight_field)
      type(csr_matrix), intent(inout) :: mass
      type(scalar_field), intent(in) :: p2b_field
      type(vector_field), intent(in) :: X
      real, dimension(7), intent(in) :: N_vals
      type(scalar_field), intent(in), optional :: weight_field
      integer, intent(in) :: ele
      !
      real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
      real, dimension(ele_ngi(X,ele)) :: detwei
      real, dimension(ele_loc(p2b_field,ele),ele_loc(p2b_field,ele))&
           :: l_mass_mat
      real :: wv = 1.0/20., we = 2.0/15., wg = 9./20.
      real :: Area
      integer :: node, node2
      real, dimension(6) :: node_weights
      real, dimension(7) :: weight_vals
      call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J, detwei)
      Area = sum(detwei)

      l_mass_mat = 0.0
      node_weights = (/wv,we,wv,we,we,wv/)
      if(present(weight_field)) then
         assert(ele_loc(weight_field,ele)==7)
         weight_vals = ele_vals(weight_field,ele)
         node_weights = node_weights*weight_vals
      end if
      do node = 1, 6
         l_mass_mat(node,node) = node_weights(node)
         do node2 = 1, 6
            l_mass_mat(node,node2) = l_mass_mat(node,node2)+wg*N_vals(node)*N_vals(node2)
         end do
      end do
      do node = 1, 6
         l_mass_mat(node,7) = wg*N_vals(node)*N_vals(7)
         l_mass_mat(7,node) = wg*N_vals(node)*N_vals(7)
      end do
      l_mass_mat(7,7) = wg*N_vals(7)**2
      L_mass_mat = L_mass_mat*Area

      call addto(mass,ele_nodes(p2b_field,ele),ele_nodes(p2b_field,ele),l_mass_mat)

    end subroutine get_lumped_mass_p2b_ele

    !! This is special cased because of the special properties of the 
    !! p2b lumped mass (it is still 3rd order, and positivity preserving)
    subroutine project_to_p2b_lumped(state,field,field_projected,p2b_mass)
      type(state_type), intent(inout) :: state
      type(scalar_field), intent(in) :: field
      type(scalar_field), intent(inout) :: field_projected
      type(csr_matrix), intent(in) :: p2b_mass
      !
      type(scalar_field) :: rhs
      type(vector_field), pointer :: X
      real, dimension(field_projected%mesh%shape%loc, field%mesh%shape%loc)&
           & :: locweight
      real, dimension(field_projected%mesh%shape%loc, field_projected%mesh%shape%loc)&
           & :: p2b_vals
      integer :: toloc,fromloc,ele
      real :: wv = 1.0/20., we = 2.0/15., wg = 9./20.
      real, dimension(7) :: weights
      real, dimension(7) :: N_vals

      if(field_projected%mesh%shape%dim.ne.2) then
         FLAbort('Only works for 2d meshes')
      end if
      if(field_projected%mesh%shape%loc.ne.7) then
         FLAbort('Expected p2 bubble mesh')
      end if
      X=>extract_vector_field(state, "Coordinate")
      call allocate(rhs,field_projected%mesh,"ProjectionRHS")
      call zero(rhs)

      N_vals = eval_shape(field_projected%mesh%shape, (/1.0/3.0,1.0/3.0,1.0/3.0/))
      weights = (/wv,we,wv,we,we,wv,wg/)

      ! First construct remapping weights.
      !! Locweight_{ij} contains value of field basis function j evaluated at
      !! p2b node locations i
      do toloc=1,size(locweight,1)
         do fromloc=1,size(locweight,2)
            locweight(toloc,fromloc)=eval_shape(field%mesh%shape, fromloc, &
                 local_coords(toloc, field_projected%mesh%shape))
         end do
      end do

      !! p2b_vals_{ij} contains value of p2b basis function j at p2b node
      !!  location i
      do toloc = 1, size(locweight,1)
         do fromloc = 1, size(locweight,1)
            p2b_vals(toloc,fromloc) = eval_shape(field_projected%mesh%shape, fromloc, &
                 local_coords(toloc, field_projected%mesh%shape))
         end do
      end do

      do ele = 1, ele_count(X)
         call project_to_p2b_lumped_ele(X,locweight,p2b_vals,weights,&
              rhs,field,ele)
      end do
      call petsc_solve(field_projected,p2b_mass,rhs)
      ewrite (2,*) maxval(abs(field_projected%val))
      call deallocate(rhs)

    end subroutine project_to_p2b_lumped

    subroutine project_to_p2b_lumped_ele(X,locweight,p2b_vals,&
         &weights,rhs,field,ele)
      type(vector_field), intent(in) :: X
      type(scalar_field), intent(in) :: field
      type(scalar_field), intent(inout) :: rhs
      integer, intent(in) :: ele
      real, dimension(ele_loc(rhs,ele),ele_loc(field,ele)), intent(in) ::&
           & locweight
      real, dimension(ele_loc(rhs,ele),ele_loc(rhs,ele)), intent(in) :: p2b_vals
      real, dimension(ele_loc(rhs,ele)), intent(in) :: weights
      !
      real, dimension(ele_loc(rhs,ele)) :: l_rhs, field_vals_rhs
      real, dimension(ele_loc(field,ele)) :: field_vals
      integer :: loc
      real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
      real, dimension(ele_ngi(X,ele)) :: detwei
      real :: area

      call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J, detwei)
      Area = sum(detwei)
      !Values of field at field DOFs
      field_vals = ele_val(field,ele)
      field_vals_rhs = matmul(locweight,field_vals)
      !! Sum over p2b DOFs of p2b basis functions multiplied by
      !! weighted values of field at p2b DOF
      l_rhs = matmul(transpose(p2b_vals),field_vals_rhs*weights*Area)
      call addto(rhs,ele_nodes(rhs,ele),l_rhs)
    end subroutine project_to_p2b_lumped_ele
  end module bubble_tools
