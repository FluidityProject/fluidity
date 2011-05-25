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

module hybridized_helmholtz
    use spud
    use fields
    use state_module
    use FLDebug
    use populate_state_module
    use write_state_module
    use populate_state_module
    use timeloop_utilities
    use sparsity_patterns_meshes
    use sparse_matrices_fields
    use solvers
    use diagnostic_variables
    use diagnostic_fields_wrapper
    use assemble_cmc
    use global_parameters, only: option_path_len
    use vector_tools, only: solve
    implicit none

contains 
  subroutine solve_hybridized_helmholtz(state,rhs)
    ! Subroutine to solve hybridized helmholtz equation
    ! If rhs (scalar pressure field) is present, then solve:
    ! <w,u> + <w,fu^\perp> - g <div w,d> + <<[w],d>> = 0
    ! <\phi,d> +  <\phi,div u> = <\phi, rhs>
    ! <<\gamma, [u]>> = 0
    ! (i.e. for unit testing)
    ! otherwise solve:
    ! <w,u> + dt*theta*<w,fu^\perp> - dt*theta*g <div w,d> + <<[w],d>> = 
    ! -dt*<w f(u^n)^\perp> + dt*g*<div w, d^n>
    ! <\phi,\eta> +  dt*theta*<\phi,div u> = <\ph
    ! <<\gamma, [u]>> = 0
    ! and then updating
    ! u = u^n + dt*u, d = d^n + dt*d
    ! (i.e. for solving wave equation)
    implicit none
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout), optional :: rhs
    !
    type(vector_field), pointer :: X, U, down
    type(scalar_field), pointer :: D,f, lambda
    type(scalar_field) :: lambda_rhs
    type(csr_sparsity) :: lambda_sparsity
    type(csr_matrix) :: lambda_mat
    real :: D0, dt, g, theta
    integer :: ele
    ewrite(1,*) '  subroutine assemble_hybridized_helmholtz('

    !Pull the fields out of state
    D=>extract_scalar_field(state, "LayerThickness")
    f=>extract_scalar_field(state, "Coriolis")
    U=>extract_vector_field(state, "LocalVelocity")
    X=>extract_vector_field(state, "Coordinate")
    down=>extract_vector_field(state, "GravityDirection")
    lambda=>extract_scalar_field(state, "LagrangeMultiplier")

    !construct/extract sparsities
    lambda_sparsity=get_csr_sparsity_firstorder(state, lambda%mesh, lambda&
         &%mesh)

    !allocate matrices
    call allocate(lambda_mat,lambda_sparsity)
    call zero(lambda_mat)

    !allocate hybridized RHS
    call allocate(lambda_rhs,lambda%mesh,"LambdaRHS")
    call zero(lambda_rhs)
    
    !get parameters
    call get_option("/physical_parameters/gravity/magnitude", g)
    !theta
    call get_option("/material_phase::Fluid/scalar_field::LayerThickness/pro&
         &gnostic/temporal_discretisation/theta",theta)
    !D0
    call get_option("/material_phase::Fluid/scalar_field::LayerThickness/p&
         &rognostic/mean_layer_thickness",D0)
    call get_option("/timestepping/timestep", dt)

    !Assemble matrices
    do ele = 1, ele_count(D)
       call assemble_hybridized_helmholtz_ele(D,f,U,X,down,ele, &
            lambda_mat,g,dt,theta,D0,lambda_rhs,rhs)
    end do

    !Solve the equations
    call petsc_solve(lambda,lambda_mat,lambda_rhs)

    !Reconstruct U and D from lambda
    do ele = 1, ele_count(D)
       call U_and_D_from_lambda_ele(D,f,U,X,down,ele, &
            g,dt,theta,D0,lambda)
    end do

    call deallocate(lambda_mat)
    call deallocate(lambda_rhs)

  end subroutine solve_hybridized_helmholtz
 
  subroutine assemble_hybridized_helmholtz_ele(D,f,U,X,down,ele, &
       lambda_mat,g,dt,theta,D0,lambda_rhs,rhs)
    implicit none
    type(scalar_field), intent(inout) :: D,f,lambda_rhs
    type(vector_field), intent(inout) :: U,X,down
    integer, intent(in) :: ele
    real, intent(in) :: g,dt,theta,D0
    type(csr_matrix), intent(inout) :: lambda_mat
    type(scalar_field), intent(inout), optional :: rhs
    !
    real, dimension(ele_loc(U,ele)*2*ele_loc(D,ele),ele_loc(U,ele)*ele_loc(D,ele))&
         &:: local_solver
    real, allocatable, dimension(:,:) :: continuity_mat, continuity_mat2
    real, allocatable, dimension(:,:) :: helmholtz_loc_mat
    real, allocatable, dimension(:,:,:) :: continuity_face_mat
    integer :: ni, lambda_ele_loc, face, ele_2, lambda_loc_count
    integer, dimension(:), pointer :: neigh
    type(element_type) :: u_shape, d_shape
    real, dimension(mesh_dim(U), X%dim, ele_ngi(U,ele)) :: J
    real, dimension(ele_ngi(x,ele)) :: f_gi
    real, dimension(X%dim, ele_ngi(X,ele)) :: up_gi
    real, dimension(X%dim) :: up_vec
    real, dimension(mesh_dim(U),ele_loc(D,ele),ele_loc(U,ele)) :: l_div_mat
    real, dimension(mesh_dim(U), mesh_dim(U), ele_ngi(U,ele)) :: Metric, &
         &Metricf
    real, dimension(X%dim, X%dim, ele_ngi(U,ele)) :: rot
    real, dimension(mesh_dim(U),mesh_dim(U),ele_loc(U,ele),ele_loc(U&
         &,ele)) :: l_u_mat
    integer :: mdim, uloc,dloc,face2,ni2,dim1,dim2,gi
    real, dimension(:), allocatable :: lambda_rhs_loc
    real, dimension(2*ele_loc(U,ele)+ele_loc(D,ele)) :: rhs_loc
    integer, dimension(:), pointer :: D_ele,U_ele
    real, dimension(ele_ngi(D,ele)) :: detwei, detJ

    mdim = mesh_dim(U)
    uloc = ele_loc(U,ele)
    dloc = ele_loc(d,ele)

    u_shape=ele_shape(u, ele)
    D_shape=ele_shape(d, ele)
    D_ele => ele_nodes(D, ele)
    U_ele => ele_nodes(U, ele)

    call get_local_solver()

    !get list of neighbours
    neigh => ele_neigh(D,ele)
    !get size of continuity_mat
    lambda_ele_loc = 0
    do ni = 1, size(neigh)
       ele_2 = neigh(ni)
       face=ele_face(U, ele, ele_2)
       lambda_ele_loc = lambda_ele_loc + face_loc(lambda_rhs,face)
    end do

    !allocate continuity_mat
    allocate(continuity_mat(ele_loc(U,ele)*2+ele_loc(D,ele),lambda_ele_loc))
    continuity_mat = 0.
    !calculate continuity_mat
    lambda_loc_count=0
    do ni = 1, size(neigh)
       face=ele_face(U, ele, neigh(ni))
       allocate(continuity_face_mat(2,face_loc(U,face),&
            &face_loc(lambda_rhs,face)))
       call get_continuity_face_mat(continuity_face_mat,X,face,&
            U,lambda_rhs)
       continuity_mat(face_local_nodes(U,face),&
            lambda_loc_count+1:lambda_loc_count+face_loc(lambda_rhs,face)) = &
            &continuity_face_mat(1,:,:)
       continuity_mat(ele_loc(U,ele)+face_local_nodes(U,face),&
            lambda_loc_count+1:lambda_loc_count+face_loc(lambda_rhs,face)) = &
            &continuity_face_mat(2,:,:)
       lambda_loc_count=lambda_loc_count+face_loc(lambda_rhs,face)
       deallocate(continuity_face_mat)
    end do

    !compute continuity_mat2 = inverse(local_solver)*continuity_mat
    allocate(continuity_mat2(ele_loc(U,ele)*2+ele_loc(D,ele),lambda_ele_loc))
    continuity_mat2 = continuity_mat
    call solve(local_solver,continuity_mat)

    !compute helmholtz_loc_mat
    allocate(helmholtz_loc_mat(lambda_ele_loc,lambda_ele_loc))
    helmholtz_loc_mat = matmul(transpose(continuity_mat),continuity_mat2)

    !construct lambda_rhs
    ! ( M    C  L )(u)   (0)
    ! ( -C^T N  0 )(h) = (j)
    ! ( L^T  0  0 )(l)   (0)
    ! 
    ! (u)   (M    C)^{-1}(0 )   (M    C)^{-1}(L)
    ! (h) = (-C^T N)     (-j) - (-C^T N)     (0)(l)
    ! so           (M    C)^{-1}(L)        (M    C)^{-1}(0 )
    ! (0) =-(L^T 0)(-C^T N)     (0) (L^T 0)(-C^T N)     (-j)
    allocate(lambda_rhs_loc(lambda_ele_loc))
    if(present(rhs)) then
       rhs_loc = ele_val(rhs,ele)
       call solve(local_solver,rhs_loc)
       lambda_rhs_loc = -matmul(transpose(continuity_mat),rhs_loc)
    else
       FLExit('Haven''t coded a timestepping version yet.')
    end if

    !insert helmholtz_loc_mat into global lambda matrix
    do ni = 1, size(neigh)
       ele_2 = neigh(ni)
       face=ele_face(U, ele, ele_2)
       call addto(lambda_rhs,ele_nodes(lambda_rhs,ele),lambda_rhs_loc)
       do ni2 = 1, size(neigh)
          ele_2 = neigh(ni)
          face2=ele_face(U, ele, ele_2)
          call addto(lambda_mat,ele_nodes(lambda_rhs,ele),&
               ele_nodes(lambda_rhs,ele),helmholtz_loc_mat)
       end do
    end do

  contains 
    subroutine get_local_solver()
      implicit none
      f_gi = ele_val_at_quad(f,ele)
      up_gi = -ele_val_at_quad(down,ele)
      call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
           detwei=detwei, detJ=detJ)
      do gi=1, ele_ngi(U,ele)
         up_vec = get_up_vec(ele_val(X,ele), up_gi(:,gi))
         up_gi(:,gi) = up_vec
      end do

      !----construct local solver
      !metrics for velocity mass and coriolis matrices
      do gi=1, ele_ngi(U,ele)
         rot(1,:,gi)=(/0.,-up_gi(3,gi),up_gi(2,gi)/)
         rot(2,:,gi)=(/up_gi(3,gi),0.,-up_gi(1,gi)/)
         rot(3,:,gi)=(/-up_gi(2,gi),up_gi(1,gi),0./)
      end do
      do gi=1,ele_ngi(U,ele)
         Metric(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)
         Metricf(:,:,gi)=matmul(J(:,:,gi), &
              matmul(rot(:,:,gi), transpose(J(:,:,gi))))/detJ(gi)
      end do

      !pressure mass matrix
      local_solver(uloc*2 + 1:uloc*2+dloc,&
           &uloc*2 + 1:uloc*2+dloc) = &
           &shape_shape(d_shape,d_shape,detwei)
      !pressure gradient matrix
      l_div_mat = dshape_shape(D_shape%dn,u_shape,D_shape%quadrature&
           &%weight)
      do dim1 = 1, mdim
         local_solver(uloc*(dim1-1)+1:uloc+dim1,uloc*2+1:uloc*2+dloc)&
              &= g*dt*theta*l_div_mat(mdim,:,:)
      end do
      !velocity mass matrix and Coriolis matrix
      l_u_mat = shape_shape_tensor(u_shape, u_shape, &
           u_shape%quadrature%weight, Metric+dt*theta*Metricf)
      do dim1 = 1, mdim
         do dim2 = 1, mdim
            local_solver(uloc*(dim1-1)+1:uloc*dim1,&
                 uloc*(dim2-1):uloc*dim2) = l_u_mat(dim1,dim2,:,:)
         end do
      end do
    end subroutine get_local_solver

  end subroutine assemble_hybridized_helmholtz_ele

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
    
    subroutine get_continuity_face_mat(continuity_face_mat,X,face,&
         U,lambda_rhs)
      ! integral is done in local coordinates to avoid computing
      ! dx/dxi on face (using properties of the Piola transform)
      ! \int_f [[w]]\lambda dS
      implicit none
      integer, intent(in) :: face
      type(scalar_field), intent(inout) :: lambda_rhs
      type(vector_field), intent(inout) :: X, U
      real, dimension(2,face_loc(U,face),face_loc(lambda_rhs,face)), intent(inout) :: continuity_face_mat
      !
      real, dimension(U%dim, face_ngi(U, face)) :: n1
      real :: weight
      type(element_type), pointer :: U_face_shape,lambda_face_shape
      real, dimension(face_ngi(U,face)) :: detwei

      U_face_shape=>face_shape(U, face)
      lambda_face_shape=>face_shape(lambda_rhs, face)

      !Get normal in local coordinates
      n1=get_normal(U,local_face_number(U%mesh,face))

      !Integral is taken on one of the edges of the local 2D element
      !This edge must be transformed to the local 1D element
      !to do numerical integration, with the following weight factors
      if(local_face_number(U%mesh, face)==3) then
         weight = sqrt(2.)
      else
         weight = 1.0
      end if
      detwei = weight*U_face_shape%quadrature%weight
      
      continuity_face_mat = shape_shape_vector(&
           U_face_shape,lambda_face_shape,detwei,n1)
     
    end subroutine get_continuity_face_mat

        function get_normal(U,face) result(norm)
          !Function returns normal to face on local 2D element
          implicit none
          type(vector_field), intent(in) :: U
          integer, intent(in) :: face
          real, dimension(U%dim, face_ngi(U,face)) :: norm
          
          integer :: i
          
          if(U%dim==1) then
             if(face==1) then
                forall(i=1:face_ngi(U,face)) norm(1,i)=1.
             else if(face==2) then
                forall(i=1:face_ngi(U,face)) norm(1,i)=-1.
             end if
          else if(U%dim==2) then
             if(face==1) then
                forall(i=1:face_ngi(U,face)) norm(1:2,i)=(/-1.,0./)
             else if(face==2) then
                forall(i=1:face_ngi(U,face)) norm(1:2,i)=(/0.,-1./)
             else if(face==3) then
                forall(i=1:face_ngi(U,face)) norm(1:2,i)=(/1/sqrt(2.),1&
                     &/sqrt(2.)/)
             else
                FLAbort('Oh dear oh dear')
             end if
          end if

      end function get_normal

end module hybridized_helmholtz
