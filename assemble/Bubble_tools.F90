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
use sparse_matrices_fields
use sparsity_patterns_meshes
use element_numbering
use state_module
use fldebug_parameters
use ieee_arithmetic, only: ieee_quiet_nan, ieee_value
use spud
implicit none

private
public nodalise_bubble_basis, setup_Cg_dg_projection

contains

  subroutine setup_Cg_Dg_projection(state,Cg_mesh,Dg_mesh)
    !Check if various matrices for Cg_Dg projection are present,
    !if not, create them.
    type(state_type), intent(inout) :: state
    type(mesh_type), intent(inout) :: Cg_mesh, Dg_mesh
    !
    type(csr_matrix) :: cgdg_proj_mat, cg_mass, dg_mass
    type(scalar_field) :: cg_lumped_mass, dg_lumped_mass
    type(csr_matrix), pointer :: matrix
    type(csr_sparsity), pointer :: sparsity
    type(scalar_field), pointer :: s_field
    type(vector_field), pointer :: X
    integer :: ele, cg_mass_stat, cg_lumped_mass_stat, &
         & dg_mass_stat, dg_lumped_mass_stat, projection_stat

    X => extract_vector_field(state,"Coordinate")

    !! Cg mass matrix
    matrix=>extract_csr_matrix(state,trim(Cg_mesh%name)//"Mass",cg_mass_stat)
    if(cg_mass_stat.ne.0) then
       sparsity => get_csr_sparsity_firstorder(&
            state, cg_mesh, cg_mesh)
       call allocate(cg_mass,sparsity,name=trim(Cg_mesh%name)//"Mass")
       call zero(cg_mass)
    end if

    !! Lumped Cg mass matrix (stored as scalar field)
    s_field=>extract_scalar_field(state,trim(Cg_mesh%name)//"LumpedMass",&
         &cg_lumped_mass_stat)
    if(cg_lumped_mass_stat.ne.0) then
       call allocate(cg_lumped_mass,Cg_mesh,&
            name=trim(Cg_mesh%name)//"LumpedMass")
       call zero(cg_lumped_mass)
    end if

    !! Dg mass matrix
    matrix=>extract_csr_matrix(state,trim(Dg_mesh%name)//"Mass",&
         &dg_mass_stat)
    if(dg_mass_stat.ne.0) then
       sparsity => get_csr_sparsity_firstorder(&
            state, dg_mesh, dg_mesh)
       call allocate(dg_mass,sparsity,name=trim(Dg_mesh%name)//"Mass")
       call zero(dg_mass)
    end if

    !! Lumped Dg mass matrix (stored as scalar field)
    s_field=>extract_scalar_field(state,trim(Dg_mesh%name)//"LumpedMass",&
         &dg_lumped_mass_stat)
    if(dg_lumped_mass_stat.ne.0) then
       call allocate(dg_lumped_mass,Dg_mesh,&
            name=trim(Dg_mesh%name)//"LumpedMass")
       call zero(dg_lumped_mass)
    end if

    !!Dg to Cg projection matrix
    matrix=>extract_csr_matrix(state,trim(Cg_mesh%name)//&
         &trim(Dg_mesh%name)//"Projection",projection_stat)
    if(projection_stat.ne.0) then
       sparsity => get_csr_sparsity_firstorder(&
            state, Cg_mesh, Dg_mesh)
       call allocate(cgdg_proj_mat,sparsity,name=trim(Cg_mesh%name)//&
            &trim(Dg_mesh%name)//"Projection")
       call zero(cgdg_proj_mat)
    end if

    if(any((/cg_mass_stat,cg_lumped_mass_stat,dg_mass_stat&
         &,dg_lumped_mass_stat,projection_stat/).ne.0)) then
       do ele = 1, element_count(dg_mesh)
          call setup_Cg_Dg_projection_ele(ele)
       end do       
    end if

    if(cg_mass_stat.ne.0) then
       call insert(state,cg_mass,cg_mass%name)
       call deallocate(cg_mass)
    end if

    if(cg_lumped_mass_stat.ne.0) then
       call insert(state,cg_lumped_mass,cg_lumped_mass%name)
       call deallocate(cg_lumped_mass)
    end if

    if(dg_mass_stat.ne.0) then
       call insert(state,dg_mass,dg_mass%name)
       call deallocate(dg_mass)
    end if

    if(dg_lumped_mass_stat.ne.0) then
       call insert(state,dg_lumped_mass,dg_lumped_mass%name)
       call deallocate(dg_lumped_mass)
    end if

    if(projection_stat.ne.0) then
       call insert(state,cgdg_proj_mat,cgdg_proj_mat%name)
       call deallocate(cgdg_proj_mat)
    end if

  contains 
    subroutine setup_Cg_Dg_projection_ele(ele)
      integer, intent(in) :: ele
      !
      real, dimension(ele_ngi(X,ele)) :: detwei, detJ
      real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
      type(element_type) :: cg_shape, dg_shape
      integer, pointer, dimension(:) :: cg_ele, dg_ele
      real, dimension(ele_loc(cg_mesh,ele),ele_loc(cg_mesh,ele)) :: &
           & l_cg_mass
      real, dimension(ele_loc(dg_mesh,ele),ele_loc(dg_mesh,ele)) :: &
           & l_dg_mass
real, dimension(ele_loc(cg_mesh,ele),ele_loc(dg_mesh,ele)) :: &
           & l_projection
      real, dimension(ele_loc(cg_mesh,ele)) :: l_cg_lumped_mass
      real :: omega_v, omega_e, omega_b

      cg_shape = ele_shape(cg_mesh,ele)
      dg_shape = ele_shape(dg_mesh,ele)
      cg_ele => ele_nodes(cg_mesh,ele)
      dg_ele => ele_nodes(dg_mesh,ele)

      call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
           detwei=detwei)

      if(cg_mass_stat.ne.0) then
         l_cg_mass = shape_shape(cg_shape,cg_shape,detwei)
         call addto(cg_mass,cg_ele,cg_ele,l_cg_mass)
      end if
      if(cg_lumped_mass_stat.ne.0) then
         if((cg_shape%dim.ne.2) .or. (cg_shape%loc.ne.7)) then
            FLAbort('Only P2 bubble currently supported')
         end if
         omega_v = 1./20.
         omega_e = 2./15.
         omega_b = 9./20.
         l_cg_lumped_mass = (/ omega_v,omega_e,omega_v,omega_e,omega_e,&
              & omega_v,omega_b /)*detwei(1)         
         call addto(cg_lumped_mass,cg_ele,l_cg_lumped_mass)
      end if
      if(dg_mass_stat.ne.0) then
         l_dg_mass = shape_shape(dg_shape,dg_shape,detweI)
         call addto(dg_mass,dg_ele,dg_ele,l_dg_mass)
      end if
      if(dg_lumped_mass_stat.ne.0) then
         l_dg_mass = shape_shape(dg_shape,dg_shape,detweI)
         call addto(dg_lumped_mass,dg_ele,sum(l_dg_mass,2))
      end if
      if(projection_stat.ne.0) then
         l_projection = shape_shape(cg_shape,dg_shape,detweI)
         call addto(cgdg_proj_mat,cg_ele,dg_ele,l_projection)
      end if

    end subroutine setup_Cg_Dg_projection_ele
  end subroutine setup_Cg_Dg_projection

  subroutine cg_dg_projection(state,cg_field,dg_field,tol)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: cg_field
    type(scalar_field), intent(inout) :: dg_field
    real, intent(in) :: tol
    !
    type(scalar_field) :: cg_residual, cg_field_tmp, dg_field_tmp
    type(csr_matrix), pointer :: cgdg_proj_mat, cg_mass, dg_mass
    type(scalar_field), pointer :: cg_lumped_mass, dg_lumped_mass
    real :: residual

    call allocate(cg_residual,cg_field%mesh,trim(cg_field%name)//'Residual')
    call allocate(cg_field_tmp,cg_field%mesh,trim(cg_field%name)//'Tmp')
    call allocate(dg_field_tmp,dg_field%mesh,trim(dg_field%name)//'Tmp')

    cgdg_proj_mat => extract_csr_matrix(state,trim(Cg_field%mesh%name)//trim(Dg_field%mesh%name)//"Projection")
    cg_mass=>extract_csr_matrix(state,trim(Cg_field%mesh%name)//"Mass")
    dg_mass=>extract_csr_matrix(state,trim(Dg_field%mesh%name)//"Mass")
    dg_lumped_mass=>extract_scalar_field(state,trim(Dg_field%mesh%name)//"LumpedMass")
    cg_lumped_mass=>extract_scalar_field(state,trim(Cg_field%mesh%name)//"LumpedMass")
    !Initial guess
    call mult_t(dg_field,cgdg_proj_mat,cg_field)
    dg_field%val = dg_field%val/dg_lumped_mass%val

    cg_dg_solver_loop: do
       !Compute residual
       call mult(cg_residual,cg_mass,cg_field)
       call mult(cg_field_tmp,cgdg_proj_mat,dg_field)
       call addto(cg_residual, cg_field_tmp, -1.0)
       residual = maxval(abs(cg_residual%val))
       ewrite(2,*) 'residual', residual
       if(residual>tol) exit cg_dg_solver_loop
       
       !Divide by lumped mass to approximate error
       cg_residual%val = cg_residual%val/cg_lumped_mass%val

       !Compute DG approximation to error
       call mult_t(dg_field_tmp,cgdg_proj_mat,cg_residual)
       dg_field%val = dg_field%val + dg_field%val/dg_lumped_mass%val       
    end do cg_dg_solver_loop

    call deallocate(cg_residual)
    call deallocate(cg_field_tmp)
    call deallocate(dg_field_tmp)
  end subroutine cg_dg_projection
  
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
