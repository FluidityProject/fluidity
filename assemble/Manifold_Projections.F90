!    Copyright (C) 2006 Imperial College London and others.
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
module manifold_projections
  use state_module
  use fields
  use fields_base

  implicit none
  
  interface project_cartesian_to_local
    module procedure project_cartesian_to_local_generic, project_cartesian_to_local_state
  end interface
                     
  interface project_local_to_cartesian
    module procedure project_local_to_cartesian_generic, project_local_to_cartesian_state
  end interface
  private 
  public :: project_cartesian_to_local, project_local_to_cartesian

  contains 

  subroutine project_cartesian_to_local_state(state, field)
    !!< Project the cartesian velocity to local coordinates
    type(state_type), intent(inout) :: state
    type(vector_field), pointer, intent(in) :: field

    integer :: ele
    type(vector_field), pointer :: X, U_local, U_cartesian

    ewrite(1,*) "In project_cartesian_to_local"

    X=>extract_vector_field(state, "Coordinate")
    U_local=>extract_vector_field(state, "Local"//field%name)
    U_cartesian=>extract_vector_field(state, field%name)

    do ele=1, element_count(U_local)

       call project_cartesian_to_local_ele(ele, X, U_local, U_cartesian)

    end do
  end subroutine project_cartesian_to_local_state
  
  subroutine project_cartesian_to_local_generic(X, in_field_cartesian, out_field_local)
    !!< Project the cartesian velocity to local coordinates
    type(vector_field), intent(in) :: X, in_field_cartesian
    type(vector_field), intent(out) :: out_field_local

    integer :: ele

    ewrite(1,*) "In project_cartesian_to_local"

    do ele=1, element_count(out_field_local)

       call project_cartesian_to_local_ele(ele, X, out_field_local, in_field_cartesian)

    end do
  end subroutine project_cartesian_to_local_generic
  
  subroutine project_cartesian_to_local_ele(ele, X, U_local, U_cartesian)
    !!< Project the cartesian velocity to local coordinates
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: X, U_cartesian
    type(vector_field), intent(inout) :: U_local

    real, dimension(mesh_dim(U_local), mesh_dim(U_local), ele_ngi(X,ele)) :: G
    real, dimension(mesh_dim(U_local), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(ele_ngi(X,ele)) :: detwei, detJ
    real, dimension(U_cartesian%dim, ele_ngi(X,ele)) :: U_quad
    real, dimension(mesh_dim(U_local)*ele_loc(U_local,ele)) :: l_rhs
    real, dimension(mesh_dim(U_local), mesh_dim(U_local), ele_loc(U_local,ele), ele_loc(U_local,ele)) :: l_mass
    real, dimension(mesh_dim(U_local)*ele_loc(U_local,ele), mesh_dim(U_local)*ele_loc(U_local,ele)) :: l_big_mat
    type(element_type), pointer :: U_shape
    integer, dimension(:), pointer :: U_ele
    integer :: dim, dim1, dim2, gi, loc, nloc

    dim=U_local%dim

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J, detwei, detJ)

    U_shape=>ele_shape(U_local,ele)
    U_quad=ele_val_at_quad(U_cartesian,ele)
    U_ele=>ele_nodes(U_local, ele)

    nloc=ele_loc(U_local,ele)
    do dim1=1, dim
       l_rhs((dim1-1)*nloc+1:dim1*nloc)=shape_rhs(U_shape, sum(&
            J(dim1,:,:)*U_quad(:,:),dim=1)*U_shape%quadrature%weight)
    end do

    do gi=1,ele_ngi(X,ele)
       G(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)
    end do

    l_mass=shape_shape_tensor(u_shape, u_shape, &
         u_shape%quadrature%weight, G)

    do dim1 = 1, dim
       do dim2 = 1, dim
          l_big_mat(nloc*(dim1-1)+1:nloc*dim1, &
               nloc*(dim2-1)+1:nloc*dim2) = &
               l_mass(dim1,dim2,:,:)
       end do
    end do

    call solve(l_big_mat, l_rhs)

    do dim1=1, U_local%dim
       do loc=1, nloc
          call set(U_local, dim1, U_ele(loc), l_rhs((dim1-1)*nloc+loc))
       end do
    end do
    
  end subroutine project_cartesian_to_local_ele

  subroutine project_local_to_cartesian_state(state, adjoint)
    !!< Project the local velocity to cartesian coordinates
    type(state_type), intent(inout) :: state
    logical, intent(in), optional :: adjoint

    integer :: ele
    type(vector_field), pointer :: X, U_local, U_cartesian

    ewrite(1,*) "In project_local_to_cartesian"

    X=>extract_vector_field(state, "Coordinate")
    if (present_and_true(adjoint)) then
      U_local=>extract_vector_field(state, "AdjointLocalVelocity")
    else
      U_local=>extract_vector_field(state, "LocalVelocity")
    endif
    U_cartesian=>extract_vector_field(state, "Velocity")

    do ele=1, element_count(U_cartesian)

       call project_local_to_cartesian_ele(ele, X, U_cartesian, U_local)

    end do

  end subroutine project_local_to_cartesian_state

  subroutine project_local_to_cartesian_generic(X, in_field_local, out_field_cartesian)
    !!< Project the local velocity to cartesian coordinates
    type(vector_field), intent(in) :: X, in_field_local
    type(vector_field), intent(out) :: out_field_cartesian

    integer :: ele

    ewrite(1,*) "In project_local_to_cartesian_generic"

    do ele=1, element_count(out_field_cartesian)
       call project_local_to_cartesian_ele(ele, X, out_field_cartesian, in_field_local)
    end do
  end subroutine project_local_to_cartesian_generic
  
  subroutine project_local_to_cartesian_ele(ele, X, U_cartesian, U_local)
    !!< Project the local velocity to cartesian
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: X, U_local
    type(vector_field), intent(inout) :: U_cartesian

    real, dimension(ele_loc(U_cartesian,ele), ele_loc(U_cartesian,ele)) :: mass
    real, dimension(mesh_dim(U_local), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(ele_ngi(U_local,ele)) :: detwei
    real, dimension(U_local%dim, ele_ngi(X,ele)) :: U_quad
    real, dimension(X%dim, ele_ngi(X,ele)) :: U_cartesian_gi
    real, dimension(X%dim, ele_loc(U_cartesian,ele)) :: rhs
    type(element_type), pointer :: U_shape
    integer :: d, gi

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, detwei=detwei)

    U_shape=>ele_shape(U_cartesian,ele)
    U_quad=ele_val_at_quad(U_local,ele)

    mass=shape_shape(U_shape, U_shape, detwei)

    call invert(mass)

    U_cartesian_gi=0.
    do gi=1, ele_ngi(X,ele)
       U_cartesian_gi(:,gi)=matmul(transpose(J(:,:,gi)),U_quad(:,gi))
    end do

    rhs=shape_vector_rhs(U_shape, U_cartesian_gi, U_shape%quadrature%weight)

    do d=1,U_cartesian%dim
       call set(U_cartesian, d, ele_nodes(U_cartesian,ele), matmul(mass,rhs(d,:)))
    end do

  end subroutine project_local_to_cartesian_ele
end module manifold_projections
