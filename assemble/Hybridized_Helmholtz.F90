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
    use FUtils, only : real_vector, real_matrix
    use global_parameters, only: option_path_len
    use vector_tools, only: solve
    use manifold_projections
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
    type(scalar_field), intent(in), optional :: D_rhs
    type(vector_field), intent(inout), optional :: U_rhs
    logical, intent(in), optional :: compute_cartesian, check_continuity,&
         & output_dense
    !
    type(vector_field), pointer :: X, U, down, U_cart
    type(scalar_field), pointer :: D,f, lambda, lambda_nc
    type(scalar_field) :: X1,X2,lambdaX1, lambdaX2
    type(scalar_field), target :: lambda_rhs, u_cpt
    type(csr_sparsity) :: lambda_sparsity, continuity_sparsity
    type(csr_matrix) :: lambda_mat, continuity_block_mat,continuity_block_mat1
    type(block_csr_matrix) :: continuity_mat
    real :: D0, dt, g, theta
    integer :: ele,i1, stat, dim1
    logical :: l_compute_cartesian,l_check_continuity, l_output_dense
    real, dimension(:,:), allocatable :: lambda_mat_dense
    character(len=OPTION_PATH_LEN) :: constraint_option_string

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
    continuity_sparsity=get_csr_sparsity_firstorder(state, u%mesh, lambda%mesh)

    !allocate matrices
    call allocate(lambda_mat,lambda_sparsity)
    call zero(lambda_mat)
    call allocate(continuity_mat,continuity_sparsity,(/U%dim,1/))
    call zero(continuity_mat)

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
            &lambda_rhs=lambda_rhs,D_rhs=D_rhs,U_rhs=U_rhs,&
            &continuity_mat=continuity_mat)
    end do

    ewrite(1,*)'LAMBDARHS MIN:MAX',minval(lambda_rhs%val),maxval(lambda_rhs%val)

    !Solve the equations
    call petsc_solve(lambda,lambda_mat,lambda_rhs)

    !Reconstruct U and D from lambda
    do ele = 1, ele_count(D)
       call reconstruct_u_d_ele(D,f,U,X,down,ele, &
            &g,dt,theta,D0,D_rhs=D_rhs,U_rhs=U_rhs,lambda=lambda)
    end do

    if(l_output_dense) then
       allocate(lambda_mat_dense(node_count(lambda),node_count(lambda)))
       lambda_mat_dense = dense(lambda_mat)
       do i1 = 1, node_count(lambda)
          ewrite(1,*) lambda_mat_dense(i1,:)
       end do
    end if

    if(l_compute_cartesian) then
       U_cart => extract_vector_field(state, "Velocity")
       call project_local_to_cartesian(X,U,U_cart)
    end if
    
    do ele = 1, ele_count(U)
       call check_divergence_ele(U,D,D_rhs,X,ele)
    end do

    if(l_check_continuity) then
       ewrite(1,*) 'Checking continuity'

       call zero(lambda_rhs)
       do dim1 = 1,U%dim
          u_cpt = extract_scalar_field(U,dim1)
          continuity_block_mat = block(continuity_mat,dim1,1)
          call mult_T_addto(lambda_rhs,continuity_block_mat,u_cpt)
       end do
       ewrite(1,*)'JUMPS MIN:MAX',minval(lambda_rhs%val),maxval(lambda_rhs&
            &%val)
       assert(maxval(abs(lambda_rhs%val))<1.0e-10)
       if(.not.l_compute_cartesian) then
          FLExit('Need to compute cartesian to check continuity')
       end if

       do ele = 1, ele_count(U)
          call check_continuity_ele(U_cart,X,ele)
       end do

    end if
    
    ! lambda_nc=>extract_scalar_field(state, "LambdaNC",stat)
    ! if(stat==0) then
    !    call zero(lambda_nc)
    !    do ele = 1, ele_count(D)
    !       call reconstruct_lambda_nc(lambda,lambda_nc,X,ele)
    !    end do
    ! end if
    ! ewrite(1,*) lambda_nc%val
    ! call deallocate(lambda_mat)
    ! call deallocate(lambda_rhs)

    ewrite(1,*) 'END subroutine solve_hybridized_helmholtz'

  end subroutine solve_hybridized_helmholtz
 
  subroutine assemble_hybridized_helmholtz_ele(D,f,U,X,down,ele, &
       g,dt,theta,D0,lambda_mat,lambda_rhs,U_rhs,D_rhs,continuity_mat)
    !subroutine to assemble hybridized helmholtz equation.
    !For assembly, must provide:
    !   lambda_mat,lambda_rhs
    !For assembly, may provide:
    !   D_rhs and U_rhs
    !   If neither are present, D_rhs reconstructed from D and U
    !   as part of an implicit timestepping algorithm

    implicit none
    type(scalar_field), intent(in) :: D,f
    type(scalar_field), intent(inout) :: lambda_rhs
    type(vector_field), intent(in) :: U,X,down
    type(vector_field), intent(in), optional :: U_rhs
    type(scalar_field), intent(in), optional :: D_rhs
    integer, intent(in) :: ele
    real, intent(in) :: g,dt,theta,D0
    type(csr_matrix), intent(inout) :: lambda_mat
    type(block_csr_matrix), intent(inout), optional :: continuity_mat
    !
    real, allocatable, dimension(:,:),target :: &
         &l_continuity_mat, l_continuity_mat2
    real, allocatable, dimension(:,:) :: helmholtz_loc_mat
    real, allocatable, dimension(:,:,:) :: continuity_face_mat
    real, allocatable, dimension(:,:) :: scalar_continuity_face_mat
    integer :: ni, face
    integer, dimension(:), pointer :: neigh
    real, dimension(ele_loc(lambda_rhs,ele)) :: lambda_rhs_loc,lambda_rhs_loc2
    real, dimension(:),allocatable,target :: Rhs_loc
    real, dimension(:,:), allocatable :: local_solver, local_solver_rhs
    type(element_type) :: U_shape
    integer :: stat, d_start, d_end, dim1, mdim, uloc,dloc, lloc
    integer, dimension(mesh_dim(U)) :: U_start, U_end
    type(real_vector), dimension(mesh_dim(U)) :: rhs_u_ptr
    real, dimension(:), pointer :: rhs_d_ptr
    type(real_matrix), dimension(mesh_dim(U)) :: &
         & continuity_mat_u_ptr
    logical :: have_constraint
    integer :: constraint_choice, n_constraints, constraints_start

    !Get some sizes
    lloc = ele_loc(lambda_rhs,ele)
    mdim = mesh_dim(U)
    uloc = ele_loc(U,ele)
    dloc = ele_loc(d,ele)
    U_shape = ele_shape(U,ele)

    have_constraint = &
         &have_option(trim(U%mesh%option_path)//"/from_mesh/constraint_type")
    n_constraints = 0
    if(have_constraint) then
       n_constraints = ele_n_constraints(U,ele)
    end if

    allocate(rhs_loc(2*uloc+dloc+n_constraints))
    allocate(local_solver(mdim*uloc+dloc+n_constraints,&
         mdim*uloc+dloc+n_constraints))
    allocate(local_solver_rhs(mdim*uloc+dloc,mdim*uloc+dloc))

    !Calculate indices in a vector containing all the U and D dofs in
    !element ele, First the u1 components, then the u2 components, then the
    !D components are stored.
    do dim1 = 1, mdim
       u_start(dim1) = uloc*(dim1-1)+1
       u_end(dim1) = uloc*dim1
    end do
    d_start = uloc*mdim + 1
    d_end   = uloc*mdim+dloc

    !Get pointers to different parts of rhs_loc and l_continuity_mat
    do dim1 = 1, mdim
       rhs_u_ptr(dim1)%ptr => rhs_loc(u_start(dim1):u_end(dim1))
    end do
    rhs_d_ptr => rhs_loc(d_start:d_end)
    allocate(l_continuity_mat(2*uloc+dloc+n_constraints,lloc))
    do dim1= 1,mdim
       continuity_mat_u_ptr(dim1)%ptr => &
            & l_continuity_mat(u_start(dim1):u_end(dim1),:)
    end do

    ! ( M    C  -L)(u)   (v)
    ! ( -C^T N  0 )(h) = (j)
    ! ( L^T  0  0 )(l)   (0)
    ! 
    ! (u)   (M    C)^{-1}(v)   (M    C)^{-1}(L)
    ! (h) = (-C^T N)     (j) + (-C^T N)     (0)(l)
    ! so
    !        (M    C)^{-1}(L)         (M    C)^{-1}(v)
    ! (L^T 0)(-C^T N)     (0)=-(L^T 0)(-C^T N)     (j)

    !Get the local_solver matrix that obtains U and D from Lambda on the
    !boundaries
    call get_local_solver(local_solver,U,X,down,D,f,ele,&
         & g,dt,theta,D0,have_constraint,local_solver_rhs)

    !!!Construct the continuity matrix that multiplies lambda in 
    !!! the U equation
    !allocate l_continuity_mat
    l_continuity_mat = 0.
    !get list of neighbours
    neigh => ele_neigh(D,ele)
    !calculate l_continuity_mat
    do ni = 1, size(neigh)
       face=ele_face(U, ele, neigh(ni))
       allocate(continuity_face_mat(mdim,face_loc(U,face)&
            &,face_loc(lambda_rhs,face)))
       continuity_face_mat = 0.
       call get_continuity_face_mat(continuity_face_mat,face,&
            U,lambda_rhs)
       do dim1 = 1, mdim
          continuity_mat_u_ptr(dim1)%ptr(face_local_nodes(U,face),&
               face_local_nodes(lambda_rhs,face))=&
          continuity_mat_u_ptr(dim1)%ptr(face_local_nodes(U,face),&
               face_local_nodes(lambda_rhs,face))+&
               continuity_face_mat(dim1,:,:)
       end do
       if(present(continuity_mat)) then
          do  dim1 = 1, mdim
             call addto(continuity_mat,dim1,1,face_global_nodes(U,face)&
                  &,face_global_nodes(lambda_rhs,face),&
                  &continuity_face_mat(dim1,:,:))
          end do
       end if

       deallocate(continuity_face_mat)
    end do

    !compute l_continuity_mat2 = inverse(local_solver)*l_continuity_mat
    allocate(l_continuity_mat2(uloc*2+dloc+n_constraints,lloc))
    l_continuity_mat2 = l_continuity_mat
    call solve(local_solver,l_continuity_mat2)

    !compute helmholtz_loc_mat
    allocate(helmholtz_loc_mat(lloc,lloc))
    helmholtz_loc_mat = matmul(transpose(l_continuity_mat),l_continuity_mat2)

    !construct lambda_rhs
    rhs_loc=0.
    lambda_rhs_loc = 0.
    call assemble_rhs_ele(Rhs_loc,D,U,X,ele,D_rhs,U_rhs)
    if(.not.(present(d_rhs).or.present(u_rhs))) then
       rhs_loc(1:d_end) = matmul(local_solver_rhs,rhs_loc(1:d_end))
    end if
    call solve(local_solver,Rhs_loc)
    lambda_rhs_loc = -matmul(transpose(l_continuity_mat),&
         &Rhs_loc)
    !insert lambda_rhs_loc into lambda_rhs
    call addto(lambda_rhs,ele_nodes(lambda_rhs,ele),lambda_rhs_loc)
    !insert helmholtz_loc_mat into global lambda matrix
    call addto(lambda_mat,ele_nodes(lambda_rhs,ele),&
         ele_nodes(lambda_rhs,ele),helmholtz_loc_mat)

  end subroutine assemble_hybridized_helmholtz_ele

  subroutine reconstruct_U_d_ele(D,f,U,X,down,ele, &
       g,dt,theta,D0,U_rhs,D_rhs,lambda)
    !subroutine to reconstruct U and D having solved for lambda
    implicit none
    type(scalar_field), intent(in) :: f,lambda
    type(scalar_field), intent(inout) :: D
    type(vector_field), intent(inout) :: U
    type(scalar_field), intent(in), optional :: D_rhs
    type(vector_field), intent(in) :: X,down
    type(vector_field), intent(in), optional :: U_rhs
    integer, intent(in) :: ele
    real, intent(in) :: g,dt,theta,D0
    !
    real, allocatable, dimension(:,:,:) :: continuity_face_mat
    integer :: ni, face
    integer, dimension(:), pointer :: neigh
    type(element_type) :: U_shape
    integer :: d_start, d_end, dim1, mdim, uloc,dloc,lloc
    integer, dimension(mesh_dim(U)) :: U_start, U_end
    type(real_vector), dimension(mesh_dim(U)) :: rhs_u_ptr
    real, dimension(:), pointer :: rhs_d_ptr
    type(real_matrix), dimension(mesh_dim(U)) :: &
         & continuity_mat_u_ptr
    real, dimension(ele_loc(lambda,ele)) :: lambda_val
    real, dimension(:),allocatable,target :: Rhs_loc
    real, dimension(:,:), allocatable :: local_solver, local_solver_rhs
    logical :: have_constraint
    integer :: n_constraints

    !Get some sizes
    lloc = ele_loc(lambda,ele)
    mdim = mesh_dim(U)
    uloc = ele_loc(U,ele)
    dloc = ele_loc(d,ele)
    U_shape = ele_shape(U,ele)

    have_constraint = &
         &have_option(trim(U%mesh%option_path)//"/from_mesh/constraint_type")
    n_constraints = 0
    if(have_constraint) then
       n_constraints = ele_n_constraints(U,ele)
    end if

    allocate(rhs_loc(2*uloc+dloc+n_constraints))
    allocate(local_solver(mdim*uloc+dloc+n_constraints,&
         mdim*uloc+dloc+n_constraints))
    allocate(local_solver_rhs(mdim*uloc+dloc,mdim*uloc+dloc))

    !Calculate indices in a vector containing all the U and D dofs in
    !element ele, First the u1 components, then the u2 components, then the
    !D components are stored.
    d_start = uloc*mdim + 1
    d_end   = uloc*mdim+dloc
    do dim1 = 1, mdim
       u_start(dim1) = uloc*(dim1-1)+1
       u_end(dim1) = uloc*dim1
    end do

    !Get pointers to different parts of rhs_loc and l_continuity_mat
    do dim1 = 1, mdim
       rhs_u_ptr(dim1)%ptr => rhs_loc(u_start(dim1):u_end(dim1))
    end do

    !Get the local_solver matrix that obtains U and D from Lambda on the
    !boundaries
    call get_local_solver(local_solver,U,X,down,D,f,ele,&
         & g,dt,theta,D0,have_constraint,local_solver_rhs=local_solver_rhs)

    !Construct the rhs sources for U from lambda
    call assemble_rhs_ele(Rhs_loc,D,U,X,ele,D_rhs,U_rhs)
    if(.not.(present(d_rhs).or.present(u_rhs))) then
       rhs_loc(1:d_end) = matmul(local_solver_rhs,rhs_loc(1:d_end))
    end if

    lambda_val = ele_val(lambda,ele)
    !get list of neighbours
    neigh => ele_neigh(D,ele)
    !calculate l_continuity_mat
    do ni = 1, size(neigh)
       face=ele_face(U, ele, neigh(ni))
       allocate(continuity_face_mat(mdim,face_loc(U,face),&
            face_loc(lambda,face)))
       continuity_face_mat = 0.
       call get_continuity_face_mat(continuity_face_mat,face,&
            U,lambda)
       do dim1 = 1, mdim
          rhs_u_ptr(dim1)%ptr(face_local_nodes(U,face)) = &
               & rhs_u_ptr(dim1)%ptr(face_local_nodes(U,face)) + &
               & matmul(continuity_face_mat(dim1,:,:),&
               &        face_val(lambda,face))
       end do
       deallocate(continuity_face_mat)
    end do

    ! ( M    C  -L)(u)   (0)
    ! ( -C^T N  0 )(h) = (j)
    ! ( L^T  0  0 )(l)   (0)
    ! 
    ! (u)   (M    C)^{-1}(0)   (M    C)^{-1}(L)
    ! (h) = (-C^T N)     (j) + (-C^T N)     (0)(l)

    call solve(local_solver,Rhs_loc)

    if((present(d_rhs).or.present(u_rhs))) then
       do dim1 = 1, mdim
          call set(U,dim1,ele_nodes(u,ele),&
               &Rhs_loc(u_start(dim1):u_end(dim1)))
       end do
       call set(D,ele_nodes(d,ele),Rhs_loc(d_start:d_end))
    else
       !We are in timestepping mode, so update U and D with deltaU and deltaD
       do dim1 = 1, mdim
          call addto(U,dim1,ele_nodes(u,ele),&
               &Rhs_loc(u_start(dim1):u_end(dim1)))
       end do
       call addto(D,ele_nodes(d,ele),Rhs_loc(d_start:d_end))
    end if

  end subroutine reconstruct_U_d_ele
  
  subroutine check_constraints_face(U1_loc,U2_loc,U,face)
    implicit none
    real, intent(in), dimension(:) :: U1_loc,U2_loc
    type(vector_field), intent(in) :: U
    integer, intent(in) :: face
    !
    real, dimension(face_loc(U,face)) :: U1_face, U2_face
    real, dimension(2) :: n
    real, dimension(face_loc(U,face)) :: nvals

    U1_face = U1_loc(face_local_nodes(U,face))
    U2_face = U2_loc(face_local_nodes(U,face))

    select case (local_face_number(U%mesh,face))
    case (1)
       n = (/ -1.,0. /)
    case (2)
       n = (/ 0.,-1. /)
    case (3)
       n = (/ 1./sqrt(2.),1./sqrt(2.) /)
    case default
       FLExit('!!?!??!')
    end select

    nvals = U1_face*n(1) + U2_face*n(2)
    assert(abs(nvals(2)-0.5*(nvals(1)+nvals(3)))<1.0e-10)
  end subroutine check_constraints_face

  subroutine get_local_solver(local_solver,U,X,down,D,f,ele,&
       & g,dt,theta,D0,have_constraint,local_solver_rhs)
    implicit none
    real, intent(in) :: g,dt,theta,D0
    type(vector_field), intent(in) :: U,X,down
    type(scalar_field), intent(in) :: D,f
    integer, intent(in) :: ele
    real, dimension(:,:)&
         &, intent(inout) :: local_solver
    real, dimension(:,:)&
         &, intent(inout), optional :: local_solver_rhs
    logical, intent(in) :: have_constraint
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
    type(constraints_type), pointer :: constraints
    integer :: i1, c_start, c_end

    mdim = mesh_dim(U)
    uloc = ele_loc(U,ele)
    dloc = ele_loc(d,ele)

    d_start = uloc*2 + 1
    d_end   = uloc*2+dloc
    do dim1 = 1, mdim
       u_start(dim1) = uloc*(dim1-1)+1
       u_end(dim1) = uloc*dim1
    end do

    local_solver = 0.
    
    u_shape=ele_shape(u, ele)
    D_shape=ele_shape(d, ele)
    D_ele => ele_nodes(D, ele)
    U_ele => ele_nodes(U, ele)

    f_gi = ele_val_at_quad(f,ele)
    up_gi = -ele_val_at_quad(down,ele)

    do gi=1, ele_ngi(U,ele)
       up_vec = get_up_vec(ele_val(X,ele), up_gi(:,gi))
    !   up_gi(:,gi) = up_vec
    end do

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

    !pressure mass matrix (done in global coordinates)
    local_solver(d_start:d_end,d_start:d_end)=&
         &shape_shape(d_shape,d_shape,detwei)
    !divergence matrix (done in local coordinates)
    l_div_mat = dshape_shape(u_shape%dn,d_shape,&
         &D_shape%quadrature%weight)
    do dim1 = 1, mdim
       !pressure gradient term [integrated by parts so minus sign]
       local_solver(u_start(dim1):u_end(dim1),d_start:d_end)=&
            & -g*dt*theta*l_div_mat(dim1,:,:)
       if(present(local_solver_rhs)) then
          local_solver_rhs(u_start(dim1):u_end(dim1),d_start:d_end)=&
               & dt*theta*l_div_mat(dim1,:,:)
       end if
       !divergence continuity term
       local_solver(d_start:d_end,u_start(dim1):u_end(dim1))=&
            & d0*dt*theta*transpose(l_div_mat(dim1,:,:))
       if(present(local_solver_rhs)) then
          local_solver_rhs(d_start:d_end,u_start(dim1):u_end(dim1))=&
               & -d0*dt*transpose(l_div_mat(dim1,:,:))
       end if
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

    if(present(local_solver_rhs)) then
       l_u_mat = shape_shape_tensor(u_shape, u_shape, &
            u_shape%quadrature%weight, -dt*Metricf)

       do dim1 = 1, mdim
          do dim2 = 1, mdim
             local_solver_rhs(u_start(dim1):u_end(dim1),&
                  u_start(dim2):u_end(dim2))=&
                  & l_u_mat(dim1,dim2,:,:)
          end do
       end do
    end if

    if(have_constraint) then
       constraints => U%mesh%shape%constraints
       c_start = d_end+1
       c_end = d_end + constraints%n_constraints
       do i1 = 1, constraints%n_constraints
          do dim1 = 1, mdim
             local_solver(d_end+i1,u_start(dim1):u_end(dim1))=&
                  &constraints%orthogonal(i1,:,dim1)
             local_solver(u_start(dim1):u_end(dim1),d_end+i1)=&
                  &constraints%orthogonal(i1,:,dim1)
          end do
       end do
    end if
        
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
    type(scalar_field), intent(in) :: lambda
    type(vector_field), intent(in) :: U
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
    call get_local_normal(n1,weight,U,local_face_number(U%mesh,face))
    detwei = weight*U_face_shape%quadrature%weight

    continuity_face_mat = shape_shape_vector(&
         U_face_shape,lambda_face_shape,detwei,n1)

  end subroutine get_continuity_face_mat

  subroutine get_scalar_continuity_face_mat(continuity_face_mat,face,&
       lambda)
    ! integral is done in local coordinates to avoid computing
    ! dx/dxi on face (using properties of the Piola transform)
    ! \int_f [[w]]\lambda dS
    implicit none
    integer, intent(in) :: face
    type(scalar_field), intent(in) :: lambda
    real, dimension(face_loc(lambda,face),face_loc(lambda,face)),&
         &intent(inout) :: continuity_face_mat
    !
    real :: weight
    type(element_type), pointer :: lambda_face_shape
    real, dimension(face_ngi(lambda,face)) :: detwei

    lambda_face_shape=>face_shape(lambda, face)

    !Integral is taken on one of the edges of the local 2D element
    !This edge must be transformed to the local 1D element
    !to do numerical integration, with the following weight factors
    if(face==3) then
       weight = sqrt(2.)
    else
       weight = 1.0
    end if

    !Get normal in local coordinates
    detwei = weight*lambda_face_shape%quadrature%weight

    continuity_face_mat = shape_shape(&
         lambda_face_shape,lambda_face_shape,detwei)

  end subroutine get_scalar_continuity_face_mat

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
             forall(i=1:face_ngi(U,face)) norm(1:2,i)=(/-1.,0./)
          else if(face==2) then
             forall(i=1:face_ngi(U,face)) norm(1:2,i)=(/ 1.,0./)
          else if(face==3) then
             forall(i=1:face_ngi(U,face)) norm(1:2,i)=(/0.,-1./)
          else if(face==4) then
             forall(i=1:face_ngi(U,face)) norm(1:2,i)=(/0.,1./)
          else
             FLAbort('Funny face?')
          end if
          weight = 1.0
       else
          FLAbort('Dimension not supported.')
       end if
    end select

  end subroutine get_local_normal

  subroutine compute_cartesian_ele(U_cart,U,X,ele)
    implicit none
    type(vector_field), intent(inout) :: U_cart
    type(vector_field), intent(in) :: U, X
    integer, intent(in) :: ele
    !
    real, dimension(ele_loc(U,ele),ele_loc(U,ele)) :: l_u_mat
    real, dimension(X%dim,ele_loc(U,ele)) :: u_rhs
    real, dimension(mesh_dim(U), ele_ngi(U,ele)) :: local_u_gi
    real, dimension(X%dim, ele_ngi(U,ele)) :: cart_u_gi
    integer :: dim1
    type(element_type) :: u_shape
    integer, dimension(mesh_dim(U)) :: U_start, U_end
    integer, dimension(:), pointer :: U_ele
    integer :: mdim, uloc, gi
    real, dimension(ele_ngi(U,ele)) :: detwei, detJ
    real, dimension(mesh_dim(U), X%dim, ele_ngi(U,ele)) :: J

    mdim = mesh_dim(U)
    uloc = ele_loc(U,ele)
    do dim1 = 1, mdim
       u_start(dim1) = uloc*(dim1-1)+1
       u_end(dim1) = uloc+dim1
    end do

    U_ele => ele_nodes(U, ele)

    u_shape=ele_shape(u, ele)
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
         detwei=detwei,detJ=detJ)

    local_u_gi = ele_val_at_quad(U,ele)
    do gi = 1, ele_ngi(U,ele)
       cart_u_gi(:,gi) = matmul(transpose(J(:,:,gi)),local_u_gi(:,gi))/detJ(gi)
    end do
    u_rhs = shape_vector_rhs(u_shape,cart_u_gi,detwei)
    l_u_mat = shape_shape(u_shape, u_shape, detwei)

    do dim1 = 1, mdim
       call solve(l_u_mat,u_rhs(dim1,:))
    end do
    
    do dim1 = 1, U_cart%dim
       call set(U_cart,dim1,u_ele,u_rhs(dim1,:))
    end do
  end subroutine compute_cartesian_ele

  subroutine check_continuity_ele(U_cart,X,ele)
    implicit none
    type(vector_field), intent(in) :: U_cart,X
    integer, intent(in) :: ele
    !
    integer, dimension(:), pointer :: neigh
    integer :: ni,face,ele2,face2

    neigh => ele_neigh(U_cart,ele)
    do ni = 1, size(neigh)
       ele2 = neigh(ni)
       face = ele_face(U_cart,ele,ele2)
       if(ele2>0) then
          face2 = ele_face(U_cart,ele2,ele)
       else
          face2 = -1
       end if
       call check_continuity_face(U_cart,X,ele,ele2,face,face2)
    end do
  end subroutine check_continuity_ele

  subroutine check_continuity_face(U_cart,X,ele,ele2,face,face2)
    !subroutine to check the continuity of normal component
    !of velocity at quadrature points
    implicit none
    type(vector_field), intent(in) :: U_cart,X
    integer, intent(in) :: face,face2,ele,ele2
    real, dimension(X%dim, face_ngi(U_cart, face)) :: n1,n2
    real, dimension(X%dim, face_ngi(U_cart, face)) :: u1,u2
    real, dimension(face_ngi(U_cart, face)) :: jump_at_quad
    !
    u1 = face_val_at_quad(U_cart,face)
    if(ele2>0) then
       u2 = face_val_at_quad(U_cart,face2)
    else
       u2 = 0.
    end if

    n1 = get_face_normal_manifold(X,ele,face)
    if(ele2>0) then
       n2 = get_face_normal_manifold(X,ele2,face2)
    else
       n2 = -n1
    end if
    jump_at_quad = sum(n1*u1+n2*u2,1)
    if(maxval(abs(jump_at_quad))>1.0e-8) then
       ewrite(1,*) 'Jump at face, face2 =', jump_at_quad
       FLAbort('stopping because of jumps')
    end if
    
  end subroutine check_continuity_face

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

  subroutine reconstruct_lambda_nc(lambda,lambda_nc,X,ele)
    type(scalar_field), intent(in) :: lambda
    type(scalar_field), intent(inout) :: lambda_nc
    type(vector_field), intent(in) :: X
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
    l_mass_mat = 0.
    do ni = 1, size(neigh)
       ele2 = neigh(ni)
       face = ele_face(X,ele,ele2)
       call get_nc_rhs_face(nc_rhs,l_mass_mat,lambda,lambda_nc,X,face)
    end do
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
         detwei=detwei)
    call solve(l_mass_mat,nc_rhs)
    call set(lambda_nc,ele_nodes(lambda_nc,ele),nc_rhs)
  end subroutine reconstruct_lambda_nc

  subroutine get_nc_rhs_face(nc_rhs,l_mass_mat,&
       lambda,lambda_nc,X,face)
    implicit none
    real, intent(inout), dimension(:) :: nc_rhs
    real, intent(inout), dimension(:,:) :: l_mass_mat
    type(scalar_field), intent(in) :: lambda, lambda_nc
    type(vector_field), intent(in) :: X
    integer, intent(in) :: face
    !
    real, dimension(face_ngi(lambda,face)) :: detwei
    type(element_type), pointer :: lambda_nc_face_shape
    real, dimension(X%dim,face_loc(X,face)) :: x_loc

    lambda_nc_face_shape => face_shape(lambda_nc,face)
    X_loc = face_val(X,face)

    call transform_facet_to_physical(X,face,detwei_f=detwei)
    nc_rhs(face_local_nodes(lambda_nc,face)) = &
         &nc_rhs(face_local_nodes(lambda_nc,face)) &
         &+ shape_rhs(lambda_nc_face_shape,face_val_at_quad(lambda,face)*detwei)
    l_mass_mat(face_local_nodes(lambda_nc,face), &
         &face_local_nodes(lambda_nc,face)) = &
         &l_mass_mat(face_local_nodes(lambda_nc,face), &
         &face_local_nodes(lambda_nc,face)) + &
         &shape_shape(lambda_nc_face_shape,lambda_nc_face_shape,detwei)
  end subroutine get_nc_rhs_face

  subroutine assemble_rhs_ele(Rhs_loc,D,U,X,ele,D_rhs,U_rhs)
    implicit none
    integer, intent(in) :: ele
    type(scalar_field), intent(in), optional, target :: D_rhs
    type(vector_field), intent(in), optional, target :: U_rhs
    type(vector_field), intent(in) :: X,U
    type(scalar_field), intent(in) :: D
    real, dimension(:), &
         &intent(inout) :: Rhs_loc
    !
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(mesh_dim(U),ele_loc(U,ele)) :: u_rhs_loc
    real, allocatable, dimension(:,:) :: u_cart_quad
    real, dimension(mesh_dim(U),ele_ngi(X,ele)) :: u_local_quad
    integer :: d_start, d_end, dim1, mdim, uloc,dloc, gi
    integer, dimension(mesh_dim(U)) :: U_start, U_end
    type(element_type) :: u_shape
    real, dimension(ele_ngi(D,ele)) :: detwei, detJ
    logical :: have_d_rhs,have_u_rhs
    type(scalar_field), pointer :: l_d_rhs
    type(vector_field), pointer :: l_u_rhs

    have_d_rhs = present(d_rhs)
    have_u_rhs = present(u_rhs)
    if(have_d_rhs) l_d_rhs => d_rhs
    if(have_u_rhs) l_u_rhs => u_rhs
    !Get some sizes
    mdim = mesh_dim(U)
    uloc = ele_loc(U,ele)
    dloc = ele_loc(d,ele)
    U_shape = ele_shape(U,ele)

    !Calculate indices in a vector containing all the U and D dofs in
    !element ele, First the u1 components, then the u2 components, then the
    !D components are stored.
    d_start = uloc*mdim + 1
    d_end   = uloc*mdim+dloc
    do dim1 = 1, mdim
       u_start(dim1) = uloc*(dim1-1)+1
       u_end(dim1) = uloc*dim1
    end do
    
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
         detJ=detJ,detwei=detwei)

    Rhs_loc = 0.
    if(have_d_rhs) then
       Rhs_loc(d_start:d_end) = shape_rhs(ele_shape(D,ele),&
            &ele_val_at_quad(l_D_rhs,ele)*detwei)
       if(have_u_rhs) then
          allocate(u_cart_quad(l_U_rhs%dim,ele_ngi(X,ele)))
          u_cart_quad = ele_val_at_quad(l_u_rhs,ele)
          do gi = 1, ele_ngi(D,ele)
             !Don't divide by detJ as we can use weight instead of detwei
             u_local_quad(:,gi) = matmul(J(:,:,gi),u_cart_quad(:,gi))
          end do
          U_rhs_loc = shape_vector_rhs(u_shape,&
               u_local_quad,u_shape%quadrature%weight)
          do dim1 = 1, mdim
             Rhs_loc(u_start(dim1):u_end(dim1)) = &
                  & U_rhs_loc(dim1,:)
          end do
       end if
    end if
    if(.not.(present(D_rhs).or.present(u_rhs))) then
       !We are in timestepping mode.
       rhs_loc(d_start:d_end) = ele_val(D,ele)
       do dim1 = 1, mdim
          rhs_loc(u_start(dim1):u_end(dim1)) = &
               & U_rhs_loc(dim1,:)
       end do
    end if

  end subroutine assemble_rhs_ele

  subroutine check_divergence_ele(U,D,D_rhs,X,ele)
    implicit none
    type(vector_field), intent(in) :: U, X
    type(scalar_field), intent(in) :: D, D_rhs
    integer, intent(in) :: ele
    !
    real, dimension(ele_ngi(D,ele)) :: detwei, detJ
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    integer :: dim1
    real, dimension(ele_ngi(X,ele)) :: div
    real, dimension(U%dim, ele_loc(U,ele)) :: U_loc
    type(element_type) :: u_shape, d_shape
    real, dimension(ele_loc(D,ele)) :: Div_loc
    real, dimension(mesh_dim(U),ele_loc(U,ele),ele_loc(D,ele)) :: l_div_mat
    real, dimension(ele_loc(D,ele),ele_loc(D,ele)) :: d_mass_mat
    !
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
         detJ=detJ,detwei=detwei)
    
    u_shape = ele_shape(U,ele)
    d_shape = ele_shape(D,ele)
    U_loc = ele_val(U,ele)
    
    div = 0.
    do dim1 = 1, U%dim
       div = div + &
            & matmul(U_loc(dim1,:),u_shape%dn(:,:,dim1))/detJ
    end do

    l_div_mat = dshape_shape(u_shape%dn,d_shape,&
         &D_shape%quadrature%weight)

    div_loc = 0.
    do dim1 = 1, U%dim
       div_loc = div_loc + matmul(transpose(l_div_mat(dim1,:,:))&
            &,U_loc(dim1,:))
    end do
    d_mass_mat = shape_shape(d_shape,d_shape,detwei)
    call solve(d_mass_mat,div_loc)
    !ewrite(1,*) 'div_loc', div_loc
    !ewrite(1,*) 'div', div
  end subroutine check_divergence_ele
end module hybridized_helmholtz
