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
    use quadrature
    use sparsity_patterns_meshes
    use element_numbering
    use state_module
    use quadrature
    use shape_functions
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

    !! This is special cased because it is only used for mass lumping
    subroutine get_p2b_lumped_mass_quadrature(quad)
      type(quadrature_type), intent(inout) :: quad
      !
      real :: wv = 1.0/20., we = 2.0/15., wg = 9./20.

      call allocate(quad, 3, 7, 3)
      quad%dim = 2
      quad%degree = 3
      quad%weight = 0.5*(/wv,we,wv,we,we,wv,wg/)
      quad%l(:,1) = (/0.,0.5,1.,0.0,0.5,1.0,1.0/3.0/)
      quad%l(:,2) = (/0.,0., 0.,0.5,0.5,0.0,1.0/3.0/)
      quad%family = 3!FAMILY_SIMPSONS

    end subroutine get_p2b_lumped_mass_quadrature

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
      type(quadrature_type) :: Simpsons_quad
      type(element_type) :: p2b_lumped_shape, X_lumped_shape,&
           & weight_lumped_shape

      X=>extract_vector_field(state, "Coordinate")

      if(p2b_field%mesh%shape%dim.ne.2) then
         FLAbort('Only works for 2d meshes')
      end if
      if(p2b_field%mesh%shape%loc.ne.7) then
         FLAbort('Expected p2 bubble mesh')
      end if

      call get_p2b_lumped_mass_quadrature(Simpsons_quad)
      p2b_lumped_shape = make_element_shape_from_element(p2b_field%mesh%shape, &
           quad=Simpsons_quad)
      X_lumped_shape = make_element_shape_from_element(X%mesh%shape, &
           quad=Simpsons_quad)

      call zero(mass)

      do ele = 1, ele_count(X)
         call get_lumped_mass_p2b_ele(mass,p2b_field,p2b_lumped_shape,&
              X,X_lumped_shape,&
              ele,weight_field)
      end do
      
      call deallocate(p2b_lumped_shape)
      call deallocate(X_lumped_shape)
      call deallocate(Simpsons_quad)

    end subroutine get_lumped_mass_p2b

    subroutine get_lumped_mass_p2b_ele(mass,p2b_field,p2b_lumped_shape,&
         X,X_lumped_shape,ele,weight_field)
      type(csr_matrix), intent(inout) :: mass
      type(scalar_field), intent(in) :: p2b_field
      type(vector_field), intent(in) :: X
      type(element_type), intent(in) :: p2b_lumped_shape,X_lumped_shape
      type(scalar_field), intent(in), optional :: weight_field
      integer, intent(in) :: ele
      !
      real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
      real, dimension(ele_ngi(X,ele)) :: detwei
      real, dimension(X_lumped_shape%ngi) :: detwei_l
      real, dimension(ele_loc(p2b_field,ele),ele_loc(p2b_field,ele))&
           :: l_mass_mat
      real :: Area
      real, dimension(p2b_lumped_shape%ngi) :: weight_vals
      call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J, detwei)
      Area = sum(detwei)
      call compute_jacobian(ele_val(X,ele), X_lumped_shape, J, detwei_l)
      assert(abs(Area-sum(detwei_l))<1.0e-10)
      assert(size(ele_val(weight_field,ele))==p2b_lumped_shape%ngi)
      if(present(weight_field)) then
         weight_vals = ele_val(weight_field,ele)
      else
         weight_vals = 1.0
      end if
      l_mass_mat = shape_shape(p2b_lumped_shape,p2b_lumped_shape,&
           weight_vals*detwei_l)
      call addto(mass,ele_nodes(p2b_field,ele),&
           ele_nodes(p2b_field,ele),l_mass_mat)

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
      integer :: ele
      type(quadrature_type) :: Simpsons_quad
      type(element_type) :: X_lumped_shape,&
           & field_projected_lumped_shape

      if(field_projected%mesh%shape%dim.ne.2) then
         FLAbort('Only works for 2d meshes')
      end if
      if(field_projected%mesh%shape%loc.ne.7) then
         FLAbort('Expected p2 bubble mesh')
      end if
      X=>extract_vector_field(state, "Coordinate")

      call get_p2b_lumped_mass_quadrature(Simpsons_quad)
      field_projected_lumped_shape = make_element_shape_from_element( &
           field_projected%mesh%shape,quad=Simpsons_quad)
      X_lumped_shape = make_element_shape_from_element(X%mesh%shape, &
           quad=Simpsons_quad)

      call allocate(rhs,field_projected%mesh,"ProjectionRHS")
      call zero(rhs)

      do ele = 1, ele_count(X)
         call project_to_p2b_lumped_ele(rhs,X,X_lumped_shape,&
              field_projected_lumped_shape,field,ele)
      end do
      call petsc_solve(field_projected,p2b_mass,rhs)
      ewrite (2,*) maxval(abs(field_projected%val))
      call deallocate(rhs)

    end subroutine project_to_p2b_lumped

    subroutine project_to_p2b_lumped_ele(rhs,X,X_lumped_shape,&
         field_projected_lumped_shape,field,ele)
      type(vector_field), intent(in) :: X
      type(scalar_field), intent(in) :: field
      type(scalar_field), intent(inout) :: rhs
      integer, intent(in) :: ele
      type(element_type), intent(in) :: X_lumped_shape,&
           field_projected_lumped_shape
      !
      real, dimension(ele_loc(rhs,ele)) :: l_rhs,field_gi
      real, dimension(ele_loc(field,ele)) :: field_vals
      integer :: loc
      real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
      real, dimension(X_lumped_shape%ngi) :: detwei

      call compute_jacobian(ele_val(X,ele), X_lumped_shape, J, detwei)

      !Values of field at field DOFs
      field_vals = ele_val(field,ele)
      assert(size(field_gi)==size(detwei))
      field_gi = matmul(transpose(field_projected_lumped_shape%n),field_vals)
      l_rhs = shape_rhs(field_projected_lumped_shape,&
           field_gi*detwei)
      
      call addto(rhs,ele_nodes(rhs,ele),l_rhs)
    end subroutine project_to_p2b_lumped_ele
  end module bubble_tools
