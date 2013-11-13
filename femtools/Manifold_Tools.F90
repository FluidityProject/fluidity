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
module manifold_tools
  use spud
  use global_parameters, only:current_debug_level, OPTION_PATH_LEN
  use state_module
  use fields
  use fields_base
  use sparse_tools
  use sparsity_patterns
  use sparse_matrices_fields
  use sparsity_patterns_meshes
  use solvers

  implicit none
  
  interface project_cartesian_to_local
    module procedure project_cartesian_to_local_generic, project_cartesian_to_local_state
  end interface
                     
  interface project_local_to_cartesian
    module procedure project_local_to_cartesian_generic, project_local_to_cartesian_state
  end interface
  private 

  interface field_stats_manifold
     module procedure field_stats_scalar_manifold, &
          field_stats_vector_manifold
  end interface field_stats_manifold

  public :: project_cartesian_to_local, project_local_to_cartesian,&
       & get_local_normal, get_face_normal_manifold, get_up_vec,&
       & get_up_gi, field_stats_manifold

  contains 

  subroutine field_stats_scalar_manifold(field, X, min, max, norm2, integral)
    !!< Return scalar statistical informaion about field.
    type(scalar_field) :: field
    !! Positions field associated with field
    type(vector_field), optional :: X
    !! Minimum value in the field.
    real, intent(out), optional :: min
    !! Maximum value in the field.
    real, intent(out), optional :: max  
    !! L2 norm of the field. This requires positions to be specified as
    !! well.
    real, intent(out), optional :: norm2
    !! Integral of the field. This requires positions to be specified as
    !! well.
    real, intent(out), optional :: integral

    ewrite(1,*) 'field_Stats_scalar_manifold '//trim(field%name)
    
    if (present(min)) then
       min=minval(field%val)
       call allmin(min)
    end if

    if (present(max)) then
       max=maxval(field%val)
       call allmax(max)
    end if

    if (present(X).and.present(norm2)) then

       norm2=norm2_scalar_manifold(field, X)
       
    elseif (present(norm2)) then
       FLAbort("Cannot evaluate L2 norm without providing positions field")
    end if

    if (present(X).and.present(integral)) then

       integral=integral_scalar_manifold(field, X)
       
    elseif (present(integral)) then
       FLAbort("Cannot evaluate integral without providing positions field")
    end if

  end subroutine field_stats_scalar_manifold

  function norm2_scalar_manifold(field, X) result (norm2)
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(in) :: X
    real :: norm2
    !
    integer :: ele

    norm2 = 0.
    do ele =1, ele_count(X)
       norm2 = norm2 + norm2_scalar_manifold_ele(field,X,ele)
    end do

    !call allsum(norm2)
    
    norm2 = sqrt(norm2)

  end function norm2_scalar_manifold

  function norm2_scalar_manifold_ele(field, X, ele) result (norm2)
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(in) :: X
    integer, intent(in) :: ele
    real :: norm2
    !
    real, dimension(ele_ngi(X,ele)) :: detwei
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    type(element_type), pointer :: field_shape
    real, dimension(ele_loc(field,ele)) :: field_val

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele),&
         detwei=detwei,J=J)

    field_val=ele_val(field, ele)
    field_shape=>ele_shape(field, ele)

    norm2 = dot_product(field_val, matmul(&
         &  shape_shape(field_shape, field_shape, detwei)&
         &                                               ,field_val))
    
  end function norm2_scalar_manifold_ele

  function integral_scalar_manifold(field, X) result (integral)
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(in) :: X
    real :: integral
    !
    integer :: ele

    integral = 0.
    do ele = 1, ele_count(X)
       integral = integral + integral_scalar_manifold_ele(&
            field,X,ele)
    end do
  end function integral_scalar_manifold

  function integral_scalar_manifold_ele(field, X, ele) result (integral)
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(in) :: X
    integer, intent(in) :: ele
    real :: integral
    !
    real, dimension(ele_ngi(X,ele)) :: detwei
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele),&
         detwei=detwei,J=J)

    integral=dot_product(ele_val_at_quad(field, ele), detwei)
    
  end function integral_scalar_manifold_ele

  subroutine field_stats_vector_manifold(field, X, min, max, norm2)
    !!< Return vector statistical informaion about field.
    type(vector_field) :: field
    !! Positions field associated with field
    type(vector_field), optional :: X
    !! Minimum value in the field.
    real, intent(out), optional :: min
    !! Maximum value in the field.
    real, intent(out), optional :: max  
    !! L2 norm of the field. This requires positions to be specified as
    !! well.
    real, intent(out), optional :: norm2

    type(scalar_field) :: mag

    ewrite(1,*) 'field_Stats_vector_manifold '//trim(field%name)

    mag=magnitude(field)

    call field_stats_manifold(mag, X, min, max, norm2)

    call deallocate(mag)

  end subroutine field_stats_vector_manifold

  subroutine get_up_gi(X,ele,up_gi,orientation)
    !subroutine to replace up_gi with a normal to the surface
    !with the same orientation
    implicit none
    type(vector_field), intent(in) :: X
    integer, intent(in) :: ele
    real, dimension(X%dim,ele_ngi(X,ele)), intent(inout) :: up_gi
    integer, intent(out), optional :: orientation
    !
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    integer :: gi
    real, dimension(X%dim,ele_ngi(X,ele)) :: normal_gi
    real, dimension(ele_ngi(X,ele)) :: orientation_gi
    integer :: l_orientation
    real :: norm

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J)

    select case(mesh_dim(X)) 
    case (2)
       do gi = 1, ele_ngi(X,ele)
          normal_gi(:,gi) = cross_product(J(1,:,gi),J(2,:,gi))
          norm = sqrt(sum(normal_gi(:,gi)**2))
          normal_gi(:,gi) = normal_gi(:,gi)/norm
       end do
       do gi = 1, ele_ngi(X,ele)
          orientation_gi(gi) = dot_product(normal_gi(:,gi),up_gi(:,gi))
       end do
       do gi = 1, ele_ngi(X,ele)
          if(sign(1.0,orientation_gi(gi)).ne.sign(1.0,orientation_gi(1))) then
             ewrite(0,*) 'orientation full=',orientation_gi
             ewrite(0,*) 'gi=',gi
             ewrite(0,*) 'normal=',normal_gi(:,gi)
             ewrite(0,*) 'up=',up_gi(:,gi)
             ewrite(0,*) 'orientation=',orientation_gi(gi)
             ewrite(0,*) 'orientation(1)=',orientation_gi(1)
             FLAbort('Nasty geometry problem')
          end if
       end do


       if(orientation_gi(1)>0.0) then
          l_orientation = 1
       else
          l_orientation = -1
       end if
       if(present(orientation)) then
          orientation =l_orientation
       end if
       do gi = 1, ele_ngi(X,ele)
          up_gi(:,gi) = normal_gi(:,gi)*l_orientation
       end do
    case default
       FLAbort('not implemented')
    end select
  end subroutine get_up_gi

  subroutine project_cartesian_to_local_state(state, field, transpose)
    !!< Project the cartesian velocity to local coordinates
    type(state_type), intent(inout) :: state
    type(vector_field), pointer, intent(in) :: field
    logical, intent(in), optional :: transpose

    integer :: ele
    type(vector_field), pointer :: X, U_local, U_cartesian

    ewrite(1,*) "In project_cartesian_to_local"

    X=>extract_vector_field(state, "Coordinate")
    U_local=>extract_vector_field(state, "Local"//field%name)
    U_cartesian=>extract_vector_field(state, field%name)

    if (present_and_true(transpose)) then
      do ele=1, element_count(U_local)
         call project_cartesian_to_local_transpose_ele(ele, X, U_cartesian, U_local)
      end do
    else
      do ele=1, element_count(U_local)
         call project_cartesian_to_local_ele(ele, X, U_local, U_cartesian)
      end do
    end if
  end subroutine project_cartesian_to_local_state
 
  ! In the case tranpose=.false., the third argument is the output 
  ! In the case tranpose=.true., the second argument is the output 
  subroutine project_cartesian_to_local_generic(X, in_field_cartesian, out_field_local, transpose)
    !!< Project the cartesian velocity to local coordinates
    type(vector_field), intent(in) :: X
    type(vector_field), intent(inout) :: out_field_local, in_field_cartesian
    logical, intent(in), optional :: transpose

    integer :: ele

    ewrite(1,*) "In project_cartesian_to_local"

    if (present_and_true(transpose)) then
      do ele=1, element_count(out_field_local)
         call project_cartesian_to_local_transpose_ele(ele, X, in_field_cartesian, out_field_local)
      end do
    else
      do ele=1, element_count(out_field_local)
         call project_cartesian_to_local_ele(ele, X, out_field_local, in_field_cartesian)
      end do
    end if
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

  subroutine project_cartesian_to_local_transpose_ele(ele, X, U_cartesian, U_local)
    !!< Project the cartesian velocity to local coordinates
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: X, U_local
    type(vector_field), intent(inout) :: U_cartesian

    real, dimension(mesh_dim(U_local), mesh_dim(U_local), ele_ngi(X,ele)) :: G
    real, dimension(mesh_dim(U_local), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(ele_ngi(X,ele)) :: detwei, detJ
    real, dimension(U_cartesian%dim, ele_ngi(X,ele)) :: U_quad
    real, dimension(mesh_dim(U_local)*ele_loc(U_local,ele)) :: l_rhs
    real, dimension(mesh_dim(U_local), mesh_dim(U_local), ele_loc(U_local,ele), ele_loc(U_local,ele)) :: l_mass
    real, dimension(mesh_dim(U_local)*ele_loc(U_local,ele), mesh_dim(U_local)*ele_loc(U_local,ele)) :: l_big_mat
    real, dimension(mesh_dim(U_local)*ele_loc(U_local,ele)) :: tmp
    real, dimension(mesh_dim(U_local),ele_loc(U_local,ele)) :: U_local_ele
    real, dimension(U_local%dim, ele_ngi(X, ele)) :: tmp_at_quad
    real, dimension(X%dim, ele_ngi(X,ele)) :: U_cartesian_gi
    real, dimension(X%dim, ele_loc(U_cartesian,ele)) :: rhs
    type(element_type), pointer :: U_shape
    integer, dimension(:), pointer :: U_ele
    integer :: dim, dim1, dim2, gi, loc, nloc

    dim=U_local%dim

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J, detwei, detJ)

    U_shape=>ele_shape(U_local,ele)
    U_quad=ele_val_at_quad(U_cartesian,ele)
    U_ele=>ele_nodes(U_local, ele)

    nloc=ele_loc(U_local,ele)

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

    U_local_ele = ele_val(U_local, ele)
    do dim1=1,U_local%dim
      tmp(nloc*(dim1-1) + 1: nloc * dim1) = U_local_ele(dim1, :)
    end do

    call solve(l_big_mat, tmp)
    do dim1=1,U_local%dim
      U_local_ele(dim1,:) = tmp(nloc*(dim1-1)+1:nloc*dim1)
    end do
    tmp_at_quad = matmul(U_local_ele, U_shape%n)

    U_cartesian_gi = 0.0
    do gi=1,ele_ngi(X, ele)
      U_cartesian_gi(:, gi) = matmul(transpose(J(:, :, gi)), tmp_at_quad(:, gi))
    end do

    rhs = shape_vector_rhs(U_shape, U_cartesian_gi, U_shape%quadrature%weight)
    do dim1 = 1, U_cartesian%dim
      call set(U_cartesian, dim1, ele_nodes(U_cartesian, ele), rhs(dim1, :))
    end do

  end subroutine project_cartesian_to_local_transpose_ele

  subroutine project_local_to_cartesian_state(state, adjoint, transpose)
    !!< Project the local velocity to cartesian coordinates
    type(state_type), intent(inout) :: state
    logical, intent(in), optional :: adjoint, transpose

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

    if (present_and_true(transpose)) then
      do ele=1, element_count(U_cartesian)
         call project_local_to_cartesian_transpose_ele(ele, X, U_local, U_cartesian)
      end do
    else
      do ele=1, element_count(U_cartesian)
         call project_local_to_cartesian_ele(ele, X, U_cartesian, U_local)
      end do
    end if

  end subroutine project_local_to_cartesian_state

  subroutine project_local_to_cartesian_generic(X, in_field_local,&
       & out_field_cartesian, transpose)
    !!< Project the local velocity to cartesian coordinates
    type(vector_field), intent(in) :: X
    type(vector_field), intent(inout) :: out_field_cartesian, in_field_local
    logical, intent(in), optional :: transpose

    integer :: ele

    ewrite(1,*) "In project_local_to_cartesian_generic"

    if (present_and_true(transpose)) then
      do ele=1, element_count(out_field_cartesian)
         call project_local_to_cartesian_transpose_ele(ele, X, in_field_local, out_field_cartesian)
      end do
    else
      do ele=1, element_count(out_field_cartesian)
         call project_local_to_cartesian_ele(ele, X, out_field_cartesian,&
              & in_field_local)
      end do
    end if
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

    U_cartesian_gi=0.
    do gi=1, ele_ngi(X,ele)
       U_cartesian_gi(:,gi)=matmul(transpose(J(:,:,gi)),U_quad(:,gi))
    end do

    rhs=shape_vector_rhs(U_shape, U_cartesian_gi, U_shape%quadrature%weight)

    do d=1,U_cartesian%dim
       call solve(mass,rhs(d,:))
       call set(U_cartesian, d, ele_nodes(U_cartesian,ele),rhs(d,:))
    end do

  end subroutine project_local_to_cartesian_ele
  
  subroutine project_local_to_cartesian_transpose_ele(ele, X, U_local, U_cartesian)
    !!< Transpose of project the local velocity to cartesian
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: X, U_cartesian
    type(vector_field), intent(inout) :: U_local

    real, dimension(ele_loc(U_cartesian,ele), ele_loc(U_cartesian,ele)) :: mass
    real, dimension(mesh_dim(U_local), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(ele_ngi(U_local,ele)) :: detwei
    real, dimension(U_local%dim, ele_ngi(X,ele)) :: U_quad
    real, dimension(X%dim, ele_ngi(X,ele)) :: U_cartesian_gi
    real, dimension(X%dim, ele_loc(U_cartesian,ele)) :: rhs
    real, dimension(U_cartesian%dim, ele_loc(U_cartesian, ele)) :: tmp, U_cartesian_ele
    type(element_type), pointer :: U_shape
    integer :: d, gi

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, detwei=detwei)

    U_shape=>ele_shape(U_cartesian,ele)

    mass=shape_shape(U_shape, U_shape, detwei)
    call invert(mass)
    U_cartesian_ele = ele_val(U_cartesian, ele)

    do d=1,U_cartesian%dim
      tmp(d, :) = matmul(mass, U_cartesian_ele(d, :))
    end do

    do d=1,U_local%dim
      call set(U_local, d, ele_nodes(U_local, ele), shape_rhs(U_shape, sum(J(d,:,:)*matmul(tmp, U_shape%n),dim=1)*U_shape%quadrature%weight))
    end do
  end subroutine project_local_to_cartesian_transpose_ele

  subroutine get_local_normal(norm,weight,U,face)
    !Function returns normal to face on local 2D element
    implicit none
    type(vector_field), intent(in) :: U
    integer, intent(in) :: face
    real, dimension(U%dim, face_ngi(U,face)), intent(out) :: norm
    real, intent(out) :: weight

    integer :: i

    select case(U%mesh%shape%numbering%family)
    case (FAMILY_SIMPLEX)
       if(U%dim==1) then
          if(face==1) then
             forall(i=1:face_ngi(U,face)) norm(1,i)=1.
          else if(face==2) then
             forall(i=1:face_ngi(U,face)) norm(1,i)=-1.
          else 
             FLAbort('Funny face?')
          end if
          weight = 1.0

       else if(U%dim==2) then
          if(face==1) then
             forall(i=1:face_ngi(U,face)) norm(1:2,i)=(/-1.,0./)
          else if(face==2) then
             forall(i=1:face_ngi(U,face)) norm(1:2,i)=(/0.,-1./)
          else if(face==3) then
             forall(i=1:face_ngi(U,face)) norm(1:2,i)=(/1/sqrt(2.),1&
                  &/sqrt(2.)/)
          else 
             FLAbort('Funny face?')
          end if

          !Integral is taken on one of the edges of the local 2D element
          !This edge must be transformed to the local 1D element
          !to do numerical integration, with the following weight factors
          if(face==3) then
             weight = sqrt(2.)
          else
             weight = 1.0
          end if

       else
          FLAbort('Dimension not supported.')
       end if
    case (FAMILY_CUBE)
       if(U%dim==2) then
          if(face==1) then
             forall(i=1:face_ngi(U,face)) norm(1:2,i)=(/0.,1./)
          else if(face==2) then
             forall(i=1:face_ngi(U,face)) norm(1:2,i)=(/ 1.,0./)
          else if(face==3) then
             forall(i=1:face_ngi(U,face)) norm(1:2,i)=(/-1.,0./)
          else if(face==4) then
             forall(i=1:face_ngi(U,face)) norm(1:2,i)=(/0.,-1./)
          else
             FLAbort('Funny face?')
          end if
          weight = 1.0
       else
          FLAbort('Dimension not supported.')
       end if
    end select

  end subroutine get_local_normal

  function get_face_normal_manifold(X,ele,face) result (normal)
    implicit none
    type(vector_field), intent(in) :: X
    integer, intent(in) :: ele, face
    real, dimension(X%dim,face_ngi(X,face)) :: normal
    !
    real, dimension(ele_ngi(X,ele)) :: detwei
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(face_ngi(X,face)) :: detwei_f
    real, dimension(mesh_dim(X)-1, X%dim, face_ngi(X,face)) :: J_f
    real, dimension(ele_loc(X,ele),ele_loc(X,ele)) :: X_mass_mat
    real, dimension(mesh_dim(X), X%dim, ele_loc(X,ele)) :: J_loc
    real, dimension(ele_loc(X,ele)) :: J_loc_rhs
    real, dimension(mesh_dim(X), X%dim, face_ngi(X,face)) :: J_face_gi
    integer :: dim1, dim2, gi
    real, dimension(X%dim) :: ele_normal_gi, edge_tangent_gi, X_mid_ele,&
         &X_mid_face
    type(element_type) :: X_shape, X_face_shape
    real, dimension(X%dim,ele_loc(X,ele)) :: X_ele
    real, dimension(X%dim,face_loc(X,face)) :: X_face

    call compute_jacobian(ele_val(X,ele),ele_shape(X,ele), J=J, &
         detwei=detwei)
    call compute_jacobian(face_val(X,face),face_Shape(X,face), J=J_f, &
         detwei=detwei_f)

    !Jacobian can be expanded without error in X function space
    !so we map it to the basis function DOFs by projection
    X_shape = ele_shape(X,ele)
    X_face_shape = face_shape(X,face)
    X_mass_mat = shape_shape(X_shape,X_shape,detwei)
    do dim1 = 1, mesh_dim(X)
       do dim2 = 1, X%dim
          J_loc_rhs = shape_rhs(X_shape,J(dim1,dim2,:)*detwei)
          call solve(X_mass_mat,J_loc_rhs)
          J_loc(dim1,dim2,:) = J_loc_rhs
       end do
    end do

    do dim1 = 1, mesh_dim(X)
       do dim2 = 1, X%dim
          J_face_gi(dim1,dim2,:) = &
               & matmul(transpose(X_face_shape%n),&
               J_loc(dim1,dim2,face_local_nodes(X,face)))
       end do
    end do

    X_mid_ele = sum(ele_val(X,ele),2)/size(ele_val(X,ele),2)
    X_mid_face = sum(face_val(X,face),2)/size(face_val(X,face),2)
    select case(X%dim)
    case (3)
       select case (mesh_dim(X))
       case (2)
          do gi = 1, face_ngi(X,face)
             !Get normal to element e on face quad points
             ele_normal_gi = cross_product(J_face_gi(1,:,gi),&
                  &J_face_gi(2,:,gi))
             ele_normal_gi = ele_normal_gi/(norm2(ele_normal_gi))
             !Get tangent to face f
             edge_tangent_gi = J_f(1,:,gi)
             edge_tangent_gi = edge_tangent_gi/norm2(edge_tangent_gi)
             !Compute normal to face f in manifold
             normal(:,gi) = cross_product(ele_normal_gi,edge_tangent_gi)
             if(dot_product(normal(:,gi),X_mid_face-X_mid_ele)<0) then
                normal(:,gi)=-normal(:,gi)
             end if
          end do
       case default
          FLAbort('dimension combination not implemented')
       end select
    case default
       FLAbort('dimension combination not implemented')
    end select

  end function get_face_normal_manifold

  function get_up_vec(X_val, up) result (up_vec_out)
    implicit none
    real, dimension(:,:), intent(in) :: X_val           !(dim,loc)
    real, dimension(:), intent(in) :: up
    real, dimension(size(X_val,1)) :: up_vec_out
    !
    real, dimension(size(X_val,1)) :: t1,t2
    ! if elements are triangles:
    if(size(X_val,2)==3) then
       t1 = X_val(:,2)-X_val(:,1)
       t2 = X_val(:,3)-X_val(:,1)
       up_vec_out(1) = t1(2)*t2(3)-t1(3)*t2(2)
       up_vec_out(2) = -(t1(1)*t2(3)-t1(3)*t2(1))
       up_vec_out(3) = t1(1)*t2(2)-t1(2)*t2(1)
       up_vec_out = up_vec_out*dot_product(up_vec_out, up)
       up_vec_out = up_vec_out/sqrt(sum(up_vec_out**2))
    else
       up_vec_out = up
    end if
  end function get_up_vec

end module manifold_tools
