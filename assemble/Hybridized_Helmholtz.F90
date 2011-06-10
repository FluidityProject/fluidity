
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
  subroutine solve_hybridized_helmholtz(state,D_rhs,U_Rhs,&
       &compute_cartesian,&
       &check_continuity,output_dense)
    ! Subroutine to solve hybridized helmholtz equation
    ! If D_rhs (scalar pressure field) is present, then solve:
    ! <w,u> + <w,fu^\perp> - g <div w,d> + <<[w],d>> = <w,U_rhs>
    ! <\phi,d> +  <\phi,div u> = <\phi, D_rhs>
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
    type(scalar_field), intent(inout), optional :: D_rhs
    type(vector_field), intent(inout), optional :: U_rhs
    logical, intent(in), optional :: compute_cartesian, check_continuity,&
         & output_dense
    !
    type(vector_field), pointer :: X, U, down, U_cart
    type(scalar_field), pointer :: D,f, lambda, lambda_nc
    type(scalar_field), target :: lambda_rhs
    type(csr_sparsity) :: lambda_sparsity
    type(csr_matrix) :: lambda_mat
    real :: D0, dt, g, theta
    integer :: ele,i1,i2, stat
    logical :: l_compute_cartesian,l_check_continuity, l_output_dense
    real, dimension(:,:), allocatable :: lambda_mat_dense

    ewrite(1,*) '  subroutine solve_hybridized_helmholtz('

    l_compute_cartesian = .false.
    if(present(compute_cartesian)) l_compute_cartesian = compute_cartesian
    if(present(check_continuity)) l_check_continuity = check_continuity
    if(l_check_continuity) l_compute_cartesian = .true.
    l_output_dense = .false.
    if(present(output_dense)) l_output_dense = output_dense

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
            &g,dt,theta,D0,lambda_mat=lambda_mat,&
            &lambda_rhs=lambda_rhs,D_rhs=D_rhs,U_rhs=U_rhs)
    end do

    !Solve the equations
    call petsc_solve(lambda,lambda_mat,lambda_rhs)

    !Reconstruct U and D from lambda
    do ele = 1, ele_count(D)
       call assemble_hybridized_helmholtz_ele(D,f,U,X,down,ele, &
            &g,dt,theta,D0,D_rhs=D_rhs,U_rhs=U_rhs,lambda=lambda)
    end do

    if(l_output_dense) then
       allocate(lambda_mat_dense(node_count(lambda),node_count(lambda)))
       lambda_mat_dense = dense(lambda_mat)
       do i1 = 1, node_count(lambda)
          ewrite(1,*) lambda_mat_dense(i1,:)
       end do
    end if

    call deallocate(lambda_mat)
    call deallocate(lambda_rhs)

    if(l_compute_cartesian) then
       U_cart => extract_vector_field(state, "Velocity")
       do ele = 1, ele_count(D)
          call compute_cartesian_ele(U_cart,U,X,ele)
       end do
    end if

    if(l_check_continuity) then
       ewrite(1,*) 'Checking continuity'
       do ele = 1, ele_count(D)
          call check_continuity_ele(U_cart,X,ele)
       end do
    end if

    lambda_nc=>extract_scalar_field(state, "LambdaNC",stat)
    if(stat==0) then
       call zero(lambda_nc)
       do ele = 1, ele_count(D)
          call reconstruct_lambda_nc(lambda,lambda_nc,X,ele)
       end do
    end if

    ewrite(1,*) 'END subroutine solve_hybridized_helmholtz'

  end subroutine solve_hybridized_helmholtz
 
  subroutine assemble_hybridized_helmholtz_ele(D,f,U,X,down,ele, &
       g,dt,theta,D0,lambda_mat,lambda_rhs,U_rhs,D_rhs,lambda)
    !subroutine to assemble hybridized helmholtz equation.
    !For assembly, must provide:
    !   lambda_mat,lambda_rhs
    !For assembly, may provide:
    !   D_rhs and U_rhs
    !   If neither are present, D_rhs reconstructed from D and U
    !   as part of an implicit timestepping algorithm

    !For reconstruction, must provide:
    !   lambda, D_rhs and U_rhs

    !lambda_rhs 
    !and also to reconstruct U and D from lambda 
    !(should provide 
    implicit none
    type(scalar_field), intent(inout) :: D,f
    type(vector_field), intent(inout) :: U,X,down
    integer, intent(in) :: ele
    real, intent(in) :: g,dt,theta,D0
    type(vector_field), intent(inout), optional :: U_rhs
    type(csr_matrix), intent(inout), optional :: lambda_mat
    type(scalar_field), intent(inout), target, optional :: lambda_rhs, D_rhs, lambda
    !
    real, dimension(ele_loc(U,ele)*2+ele_loc(D,ele),&
         &ele_loc(U,ele)*2+ele_loc(D,ele)) :: local_solver
    real, allocatable, dimension(:,:) :: continuity_mat, continuity_mat2
    real, allocatable, dimension(:,:) :: helmholtz_loc_mat
    real, allocatable, dimension(:,:,:) :: continuity_face_mat
    integer :: ni, lambda_ele_loc, face, ele_2, lambda_loc_count
    integer :: l_face_start, l_face_end
    integer, dimension(:), pointer :: neigh
    real, dimension(:), allocatable :: lambda_rhs_loc
    real, dimension(mesh_dim(U),ele_loc(U,ele)) :: U_rhs_loc
    real, dimension(2*ele_loc(U,ele)+ele_loc(D,ele)) :: Rhs_loc
    integer :: face2,ni2
    logical :: assembly, have_D_rhs
    type(element_type) :: U_shape
    integer :: d_start, d_end, dim1, dim2, mdim, uloc,dloc, lloc
    integer, dimension(mesh_dim(U)) :: U_start, U_end
    type(scalar_field), pointer :: l_lambda
    real, dimension(ele_ngi(D,ele)) :: detwei
    real, dimension(mesh_dim(U), X%dim, ele_ngi(U,ele)) :: J

    mdim = mesh_dim(U)
    uloc = ele_loc(U,ele)
    dloc = ele_loc(d,ele)
    U_shape = ele_shape(U,ele)

    !Check if we are in assembly mode or reconstruction mode and check if
    !the correct things have been passed in for those modes
    assembly = .false.
    have_D_rhs = .false.
    if(present(lambda_mat)) assembly = .true.
    if(present(D_rhs)) have_D_rhs = .true.
    if(assembly) then
       if(.not.present(lambda_rhs)) then
          FLAbort('Need lambda_rhs for assembly')
       end if
       if(present(lambda)) then
          FLAbort('Don''t need lambda for assembly')
       end if
       lloc = ele_loc(lambda_rhs,ele)
       l_lambda => lambda_rhs
    else
       if(present(lambda_rhs)) then
          FLAbort('Shouldn''t provide lambda_rhs for reconstruction')
       end if
       if(.not.present(D_rhs)) then
          FLAbort('Need D_rhs for reconstruction')
       end if
       if(.not.present(lambda)) then
          FLAbort('Need lambda for reconstruction')
       end if
       lloc = ele_loc(lambda,ele)
       l_lambda => lambda
    end if

    !Calculate indices in a vector containing all the U and D dofs in
    !element ele, First the u1 components, then the u2 components, then the
    !D components are stored.
    d_start = uloc*mdim + 1
    d_end   = uloc*mdim+dloc
    do dim1 = 1, mdim
       u_start(dim1) = uloc*(dim1-1)+1
       u_end(dim1) = uloc*dim1
    end do

    !Get the local_solver matrix that obtains U and D from Lambda on the
    !boundaries
    call get_local_solver(local_solver,U,X,down,D,f,ele,&
         & g,dt,theta,D0)

    !!!Construct the continuity matrix that multiplies lambda in 
    !!! the U equation
    !allocate continuity_mat
    allocate(continuity_mat(ele_loc(U,ele)*2+ele_loc(D,ele),lloc))
    continuity_mat = 0.
    !get list of neighbours
    neigh => ele_neigh(D,ele)
    !calculate continuity_mat
    do ni = 1, size(neigh)
       face=ele_face(U, ele, neigh(ni))
       allocate(continuity_face_mat(mdim,face_loc(U,face),lloc))
       continuity_face_mat = 0.
       call get_continuity_face_mat(continuity_face_mat,face,&
            U,l_lambda)
       do dim1 = 1, mdim
          continuity_mat((dim1-1)*ele_loc(U,face)+face_local_nodes(U,face),&
               &face_local_nodes(l_lambda,face))=&
               &continuity_mat((dim1-1)*ele_loc(U,face)&
               &+face_local_nodes(U,face),&
               &face_local_nodes(l_lambda,face))+&
               continuity_face_mat(dim1,:,:)
       end do
       deallocate(continuity_face_mat)
    end do

    !compute continuity_mat2 = inverse(local_solver)*continuity_mat
    allocate(continuity_mat2(ele_loc(U,ele)*2+ele_loc(D,ele),lloc))
    continuity_mat2 = continuity_mat
    call solve(local_solver,continuity_mat)

    if(assembly) then
       !compute helmholtz_loc_mat
       allocate(helmholtz_loc_mat(lloc,lloc))
       helmholtz_loc_mat = matmul(transpose(continuity_mat),continuity_mat2)

       !detwei is needed for pressure mass matrix
       call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
            detwei=detwei)
              
       !construct lambda_rhs
       ! ( M    C  -L)(u)   (0)
       ! ( -C^T N  0 )(h) = (j)
       ! ( L^T  0  0 )(l)   (0)
       ! 
       ! (u)   (M    C)^{-1}(0)   (M    C)^{-1}(L)
       ! (h) = (-C^T N)     (j) + (-C^T N)     (0)(l)
       ! so
       !        (M    C)^{-1}(L)         (M    C)^{-1}(0)
       ! (L^T 0)(-C^T N)     (0)=-(L^T 0)(-C^T N)     (j)
       allocate(lambda_rhs_loc(ele_loc(lambda_rhs,ele)))
       if(present(D_rhs)) then
          Rhs_loc = 0.
          Rhs_loc(d_start:d_end) = shape_rhs(ele_shape(D,ele),&
               &ele_val_at_quad(D_rhs,ele)*detwei)
          if(present(u_rhs)) then
             U_rhs_loc = shape_vector_rhs(ele_shape(U_rhs,ele),&
                  ele_val_at_quad(u_rhs,ele),u_shape%quadrature%weight)

             do dim1 = 1, mdim
                Rhs_loc(u_start(dim1):u_end(dim1)) = &
                     & U_rhs_loc(dim1,:)
             end do
          end if
          call solve(local_solver,Rhs_loc)
          lambda_rhs_loc = -matmul(transpose(continuity_mat),Rhs_loc)
       else
          FLExit('Haven''t coded a timestepping version yet.')
       end if
       
       !insert lambda_rhs_loc into lambda_rhs
       call addto(lambda_rhs,ele_nodes(lambda_rhs,ele),lambda_rhs_loc)
       !insert helmholtz_loc_mat into global lambda matrix
       call addto(lambda_mat,ele_nodes(lambda_rhs,ele),&
            ele_nodes(lambda_rhs,ele),helmholtz_loc_mat)

    else
       !Reconstruct U and D from lambda
       rhs_loc = matmul(continuity_mat,ele_val(lambda,ele))
       if(present(D_rhs)) then
          call assemble_rhs_reconstruction(Rhs_loc(d_start:d_end),&
               & D_rhs,X,ele)
       else
          FLExit('Haven''t coded a timestepping version yet.')             
       end if

       call solve(local_solver,Rhs_loc)
       
       do dim1 = 1, mdim
          call set(U,dim1,ele_nodes(u,ele),&
               &Rhs_loc(u_start(dim1):u_end(dim1)))
       end do
       call set(D,ele_nodes(d,ele),Rhs_loc(d_start:d_end))
    end if

  end subroutine assemble_hybridized_helmholtz_ele

  subroutine get_local_solver(local_solver,U,X,down,D,f,ele,&
       & g,dt,theta,D0)
    implicit none
    real, intent(in) :: g,dt,theta,D0
    type(vector_field), intent(inout) :: U,X,down
    type(scalar_field), intent(inout) :: D,f
    integer, intent(in) :: ele
    real, dimension(ele_loc(U,ele)*2+ele_loc(D,ele),&
         &ele_loc(U,ele)*2+ele_loc(D,ele))&
         &, intent(inout) :: local_solver
    !
    real, dimension(mesh_dim(U), X%dim, ele_ngi(U,ele)) :: J
    real, dimension(ele_ngi(x,ele)) :: f_gi
    real, dimension(X%dim, ele_ngi(X,ele)) :: up_gi
    real, dimension(X%dim) :: up_vec
    real, dimension(mesh_dim(U),ele_loc(U,ele),ele_loc(D,ele)) :: l_div_mat
    real, dimension(mesh_dim(U), mesh_dim(U), ele_ngi(U,ele)) :: Metric, &
         &Metricf
    real, dimension(X%dim, X%dim, ele_ngi(U,ele)) :: rot
    real, dimension(mesh_dim(U),mesh_dim(U),ele_loc(U,ele),ele_loc(U&
         &,ele)) :: l_u_mat
    integer :: mdim, uloc,dloc,dim1,dim2,gi
    type(element_type) :: u_shape, d_shape
    real, dimension(ele_ngi(D,ele)) :: detwei, detJ
    integer, dimension(:), pointer :: D_ele,U_ele
    integer :: d_start, d_end
    integer, dimension(mesh_dim(U)) :: U_start, U_end

    mdim = mesh_dim(U)
    uloc = ele_loc(U,ele)
    dloc = ele_loc(d,ele)
    
    u_shape=ele_shape(u, ele)
    D_shape=ele_shape(d, ele)
    D_ele => ele_nodes(D, ele)
    U_ele => ele_nodes(U, ele)

    f_gi = ele_val_at_quad(f,ele)
    up_gi = -ele_val_at_quad(down,ele)

    !do gi=1, ele_ngi(U,ele)
    !   up_vec = get_up_vec(ele_val(X,ele), up_gi(:,gi))
    !   up_gi(:,gi) = up_vec
    !end do

    !J, detJ is needed for Piola transform
    !detwei is needed for pressure mass matrix
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
         detwei=detwei, detJ=detJ)

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
            matmul(f_gi(gi)*rot(:,:,gi), transpose(J(:,:,gi))))/detJ(gi)
    end do

    !<w,u> + dt*theta*<w,fu^\perp> - g*dt*theta<div w,h> = -<w.n,l>
    !dt*theta*<\phi,div u>         + D_0<\phi,h>         = 0 

    d_start = uloc*2 + 1
    d_end   = uloc*2+dloc
    do dim1 = 1, mdim
       u_start(dim1) = uloc*(dim1-1)+1
       u_end(dim1) = uloc*dim1
    end do
    !pressure mass matrix (done in global coordinates)
    local_solver(d_start:d_end,d_start:d_end)=&
         &shape_shape(d_shape,d_shape,detwei)
    !divergence matrix (done in local coordinates)
    l_div_mat = dshape_shape(u_shape%dn,d_shape,&
         &D_shape%quadrature%weight)
    do dim1 = 1, mdim
       !pressure gradient term
       local_solver(u_start(dim1):u_end(dim1),d_start:d_end)=&
            & -g*dt*theta*l_div_mat(mdim,:,:)
       !divergence continuity term
       local_solver(d_start:d_end,u_start(dim1):u_end(dim1))=&
            & d0*dt*theta*transpose(l_div_mat(mdim,:,:))
    end do
    !velocity mass matrix and Coriolis matrix (done in local coordinates)
    l_u_mat = shape_shape_tensor(u_shape, u_shape, &
         u_shape%quadrature%weight, Metric+dt*theta*Metricf)
    do dim1 = 1, mdim
       do dim2 = 1, mdim
          local_solver(u_start(dim1):u_end(dim1),&
               u_start(dim2):u_end(dim2))=&
               & l_u_mat(dim1,dim2,:,:)
       end do
    end do
  end subroutine get_local_solver

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

  subroutine get_continuity_face_mat(continuity_face_mat,face,&
       U,lambda)
    ! integral is done in local coordinates to avoid computing
    ! dx/dxi on face (using properties of the Piola transform)
    ! \int_f [[w]]\lambda dS
    implicit none
    integer, intent(in) :: face
    type(scalar_field), intent(inout) :: lambda
    type(vector_field), intent(inout) :: U
    real, dimension(mesh_dim(U),face_loc(U,face),face_loc(lambda,face)),&
         &intent(inout) :: continuity_face_mat
    !
    real, dimension(U%dim, face_ngi(U, face)) :: n1
    real :: weight
    type(element_type), pointer :: U_face_shape,lambda_face_shape
    real, dimension(face_ngi(U,face)) :: detwei

    U_face_shape=>face_shape(U, face)
    lambda_face_shape=>face_shape(lambda, face)

    !Get normal in local coordinates
    n1=get_local_normal(U,local_face_number(U%mesh,face))

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

  function get_local_normal(U,face) result(norm)
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

  end function get_local_normal

  subroutine compute_cartesian_ele(U_cart,U,X,ele)
    implicit none
    type(vector_field), intent(inout) :: U_cart, U, X
    integer, intent(in) :: ele
    !
    real, dimension(mesh_dim(U), X%dim, ele_ngi(U,ele)) :: J
     real, dimension(ele_loc(U,ele),ele_loc(U,ele)) :: l_u_mat
    real, dimension(X%dim,ele_loc(U,ele)) :: u_rhs
    real, dimension(mesh_dim(U), ele_ngi(U,ele)) :: local_u_gi
    real, dimension(X%dim,ele_ngi(U,ele)) :: u_cart_gi
    integer :: gi, dim1, dim2
    type(element_type) :: u_shape
    integer :: d_start, d_end
    integer, dimension(mesh_dim(U)) :: U_start, U_end
    integer, dimension(:), pointer :: U_ele
    integer :: mdim, uloc
    real, dimension(ele_ngi(U,ele)) :: detwei

    mdim = mesh_dim(U)
    uloc = ele_loc(U,ele)
    do dim1 = 1, mdim
       u_start(dim1) = uloc*(dim1-1)+1
       u_end(dim1) = uloc+dim1
    end do

    U_ele => ele_nodes(U, ele)

    u_shape=ele_shape(u, ele)
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
         detwei=detwei)

    local_u_gi = ele_val_at_quad(U,ele)
    u_rhs = shape_vector_rhs(u_shape,local_u_gi,U_shape%quadrature%weight)
    l_u_mat = shape_shape(u_shape, u_shape, detwei)

    do dim1 = 1, mdim
       call solve(l_u_mat,u_rhs(dim1,:))
    end do
    
    do dim1 = 1, mdim
       call set(U_cart,dim1,u_ele,u_rhs(dim1,:))
    end do
  end subroutine compute_cartesian_ele

  subroutine check_continuity_ele(U_cart,X,ele)
    implicit none
    type(vector_field), intent(inout) :: U_cart,X
    integer, intent(in) :: ele
    !
    integer, dimension(:), pointer :: neigh
    integer :: ni,face,ele2,face2

    neigh => ele_neigh(X,ele)
    do ni = 1, size(neigh)
       ele2 = neigh(ni)
       if(ele2<ele) then
          face = ele_face(X,ele,ele2)
          if(ele2>0) then
             face2 = ele_face(X,ele2,ele)
          else
             face2 = -1
          end if
          call check_continuity_face(U_cart,X,face,face2)
       end if
    end do
  end subroutine check_continuity_ele

  subroutine check_continuity_face(U_cart,X,face,face2)
    !subroutine to check the continuity of normal component
    !of velocity at quadrature points
    implicit none
    type(vector_field), intent(inout) :: U_cart,X
    integer, intent(in) :: face,face2
    real, dimension(U_cart%dim, face_ngi(U_cart, face)) :: n1,n2
    real, dimension(U_cart%dim, face_ngi(U_cart, face)) :: u1,u2
    real, dimension(face_ngi(U_cart, face)) :: jump_at_quad
    !
    !Get normals
    call transform_facet_to_physical(X, face, &
         &                          normal=n1)
    if(face2>0) then
       call transform_facet_to_physical(X, face2, &
            &                          normal=n2)
    end if
    !
    u1 = face_val_at_quad(U_cart,face)
    if(face2>0) then
       u2 = face_val_at_quad(U_cart,face2)
    else
       u2 = 0.
    end if
    
    jump_at_quad = sum(n1*u1+n2*u2,1)

    if(maxval(abs(jump_at_quad))>1.0e-8) then
       ewrite(1,*) 'Jump at face, face2 =', jump_at_quad
    end if
    
  end subroutine check_continuity_face

  subroutine reconstruct_lambda_nc(lambda,lambda_nc,X,ele)
    type(scalar_field), intent(inout) :: lambda, lambda_nc
    type(vector_field), intent(inout) :: X
    integer, intent(in) :: ele
    !
    real, dimension(ele_loc(lambda_nc,ele)) :: nc_rhs
    type(element_type), pointer :: lambda_nc_shape
    integer, dimension(:), pointer :: neigh
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(ele_ngi(lambda,ele)) :: detwei
    integer :: ni,ele2,face
    real, dimension(ele_loc(lambda_nc,ele),ele_loc(lambda_nc,ele)) :: &
         & l_mass_mat

    nc_rhs = 0.
    lambda_nc_shape => ele_shape(lambda_nc,ele)
    neigh => ele_neigh(lambda,ele)
    do ni = 1, size(neigh)
       ele2 = neigh(ni)
       face = ele_face(X,ele,ele2)
       call get_nc_rhs_face(nc_rhs,lambda,lambda_nc,X,face)
    end do
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
         detwei=detwei)
    l_mass_mat = shape_shape(lambda_nc_shape,lambda_nc_shape,detwei)

    call solve(l_mass_mat,nc_rhs)
    call set(lambda_nc,ele_nodes(lambda_nc,ele),nc_rhs)
  end subroutine reconstruct_lambda_nc

  subroutine get_nc_rhs_face(nc_rhs,&
       lambda,lambda_nc,X,face)
    implicit none
    real, intent(inout), dimension(:) :: nc_rhs
    type(scalar_field), intent(inout) :: lambda, lambda_nc
    type(vector_field), intent(inout) :: X
    integer, intent(in) :: face
    !
    real, dimension(face_ngi(lambda,face)) :: detwei
    type(element_type), pointer :: lambda_nc_face_shape

    lambda_nc_face_shape => face_shape(lambda_nc,face)
    call transform_facet_to_physical(X,face,detwei_f=detwei)
    nc_rhs(face_local_nodes(lambda_nc,face)) = &
         &nc_rhs(face_local_nodes(lambda_nc,face)) &
         &+ shape_rhs(lambda_nc_face_shape,face_val_at_quad(lambda,face)*detwei)
  end subroutine get_nc_rhs_face

  subroutine assemble_rhs_reconstruction(Rhs_loc,&
       & D_rhs,X,ele)
    integer, intent(in) :: ele
    type(scalar_field), intent(in) :: D_rhs
    type(vector_field), intent(in) :: X
    real, dimension(ele_loc(D_rhs,ele)), intent(inout) :: Rhs_loc
    !
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(ele_ngi(X,ele)) :: detwei
    
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
         detwei=detwei)
    Rhs_loc = Rhs_loc + shape_rhs(ele_shape(D_rhs,ele),&
         ele_val_at_quad(D_rhs,ele)*detwei)
  end subroutine assemble_rhs_reconstruction

end module hybridized_helmholtz
