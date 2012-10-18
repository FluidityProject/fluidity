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
       &D_out,U_out,&
       &compute_cartesian,&
       &check_continuity,output_dense,&
       &projection,poisson,u_rhs_local,&
       &solver_option_path)
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
    type(scalar_field), intent(inout), optional :: D_out
    type(vector_field), intent(inout), optional :: U_out
    logical, intent(in), optional :: compute_cartesian, &
         &check_continuity,output_dense, projection,poisson
    logical, intent(in), optional :: u_rhs_local !means u_rhs is in local coords
    character(len=OPTION_PATH_LEN), intent(in), optional :: solver_option_path
    !
    type(vector_field), pointer :: X, U, down, U_cart
    type(scalar_field), pointer :: D,f, lambda_nc
    type(scalar_field) :: lambda
    type(scalar_field), target :: lambda_rhs, u_cpt
    type(csr_sparsity) :: lambda_sparsity, continuity_sparsity
    type(csr_matrix) :: lambda_mat, continuity_block_mat,continuity_block_mat1
    type(block_csr_matrix) :: continuity_mat
    type(mesh_type), pointer :: lambda_mesh
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

    lambda_mesh=>extract_mesh(state, "VelocityMeshTrace")
    call allocate(lambda,lambda_mesh,name="LagrangeMultiplier")

    U_cart => extract_vector_field(state, "Velocity")

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
    call get_option("/material_phase::Fluid/scalar_field::LayerThickness/&
         &prognostic/temporal_discretisation/theta",theta)
    !D0
    call get_option("/material_phase::Fluid/scalar_field::LayerThickness/&
         &p&
         &rognostic/mean_layer_thickness",D0)
    call get_option("/timestepping/timestep", dt)
    
    !Assemble matrices
    do ele = 1, ele_count(D)
       call assemble_hybridized_helmholtz_ele(D,f,U,X,down,ele, &
            &g,dt,theta,D0,lambda_mat=lambda_mat,&
            &lambda_rhs=lambda_rhs,D_rhs=D_rhs,U_rhs=U_rhs,&
            &continuity_mat=continuity_mat,&
            &projection=projection,poisson=poisson,&
            &u_rhs_local=u_rhs_local)
    end do

    ewrite(1,*) 'LAMBDARHS', maxval(abs(lambda_rhs%val))

    !Solve the equations
    if(present(solver_option_path)) then
       call petsc_solve(lambda,lambda_mat,lambda_rhs,&
            option_path=solver_option_path)
    else
       call petsc_solve(lambda,lambda_mat,lambda_rhs,&
            option_path=trim(U_cart%mesh%option_path)//&
            &"/from_mesh/constraint_type")
    end if
    ewrite(1,*) 'LAMBDA', maxval(abs(lambda%val))

    !Reconstruct U and D from lambda
    do ele = 1, ele_count(D)
       call reconstruct_u_d_ele(D,f,U,X,down,ele, &
            &g,dt,theta,D0,D_rhs=D_rhs,U_rhs=U_rhs,lambda=lambda,&
            &D_out=D_out,U_out=U_out,&
            &projection=projection,poisson=poisson,&
            &u_rhs_local=u_rhs_local)
    end do

    if(l_output_dense) then
       allocate(lambda_mat_dense(node_count(lambda),node_count(lambda)))
       lambda_mat_dense = dense(lambda_mat)
       ewrite(1,*) '-----------'
       do i1 = 1, node_count(lambda)
          ewrite(1,*) lambda_mat_dense(i1,:)
       end do
       ewrite(1,*) '-----------'
    end if

    if(l_compute_cartesian) then
       U_cart => extract_vector_field(state, "Velocity")
       if(present(U_out)) then
          call project_local_to_cartesian(X,U_out,U_cart)
       else
          call project_local_to_cartesian(X,U,U_cart)
       end if
    end if
    if(l_check_continuity) then
       ewrite(1,*) 'Checking continuity'

       call zero(lambda_rhs)
       do dim1 = 1,U%dim
          if(present(U_out)) then
             u_cpt = extract_scalar_field(U_out,dim1)
          else
             u_cpt = extract_scalar_field(U,dim1)
          end if
          continuity_block_mat = block(continuity_mat,dim1,1)
          call mult_T_addto(lambda_rhs,continuity_block_mat,u_cpt)
          ewrite(1,*) 'U, lambda',&
               &maxval(u_cpt%val), maxval(abs(lambda_rhs%val))
       end do
       ewrite(1,*)'JUMPS MIN:MAX',minval(lambda_rhs%val),&
            &maxval(lambda_rhs%val)
       assert(maxval(abs(lambda_rhs%val))<1.0e-10)
       
       ewrite(1,*) 'D MAXABS', maxval(abs(D%val))
       do ele = 1, ele_count(U)
          call check_continuity_ele(U_cart,X,ele)
       end do
    end if
    
    call deallocate(lambda_mat)
    call deallocate(lambda_rhs)
    call deallocate(lambda)

    ewrite(1,*) 'END subroutine solve_hybridized_helmholtz'

  end subroutine solve_hybridized_helmholtz
 
  subroutine assemble_hybridized_helmholtz_ele(D,f,U,X,down,ele, &
       g,dt,theta,D0,lambda_mat,lambda_rhs,U_rhs,D_rhs,&
       continuity_mat,projection,poisson,u_rhs_local)
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
    logical, intent(in), optional :: projection, poisson,u_rhs_local
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
    real, dimension(:,:), allocatable :: local_solver_matrix, local_solver_rhs
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
    allocate(local_solver_matrix(mdim*uloc+dloc+n_constraints,&
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
    call get_local_solver(local_solver_matrix,U,X,down,D,f,ele,&
         & g,dt,theta,D0,have_constraint,local_solver_rhs,&
         projection)

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
    call solve(local_solver_matrix,l_continuity_mat2)

    !compute helmholtz_loc_mat
    allocate(helmholtz_loc_mat(lloc,lloc))
    helmholtz_loc_mat = matmul(transpose(l_continuity_mat),l_continuity_mat2)

    !construct lambda_rhs
    rhs_loc=0.
    lambda_rhs_loc = 0.
    call assemble_rhs_ele(Rhs_loc,D,U,X,ele,D_rhs,U_rhs,u_rhs_local)
    if(.not.(present(d_rhs).or.present(u_rhs)))then
       assert(.not.present_and_true(projection))
       assert(.not.present_and_true(poisson))
       rhs_loc(1:d_end) = matmul(local_solver_rhs,rhs_loc(1:d_end))
    end if
    call solve(local_solver_matrix,Rhs_loc)
    lambda_rhs_loc = -matmul(transpose(l_continuity_mat),&
         &Rhs_loc)
    !insert lambda_rhs_loc into lambda_rhs
    call addto(lambda_rhs,ele_nodes(lambda_rhs,ele),lambda_rhs_loc)
    !insert helmholtz_loc_mat into global lambda matrix
    call addto(lambda_mat,ele_nodes(lambda_rhs,ele),&
         ele_nodes(lambda_rhs,ele),helmholtz_loc_mat)

  end subroutine assemble_hybridized_helmholtz_ele

  subroutine reconstruct_U_d_ele(D,f,U,X,down,ele, &
       g,dt,theta,D0,U_rhs,D_rhs,lambda,&
       &D_out,U_out,projection,poisson,u_rhs_local)
    !subroutine to reconstruct U and D having solved for lambda
    implicit none
    type(scalar_field), intent(in) :: f,lambda
    type(scalar_field), intent(inout) :: D
    type(vector_field), intent(inout) :: U
    type(vector_field), intent(in) :: X,down
    type(scalar_field), intent(in), optional :: D_rhs
    type(vector_field), intent(in), optional :: U_rhs
    type(scalar_field), intent(inout), optional :: D_out
    type(vector_field), intent(inout), optional :: U_out
    integer, intent(in) :: ele
    real, intent(in) :: g,dt,theta,D0
    logical, intent(in), optional :: projection, poisson,u_rhs_local
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
    real, dimension(:,:), allocatable :: local_solver_matrix, local_solver_rhs
    real, dimension(mesh_dim(U),ele_loc(U,ele)) :: U_solved
    real, dimension(ele_loc(D,ele)) :: D_solved
    logical :: have_constraint
    integer :: n_constraints, i1
    type(constraints_type), pointer :: constraints
    real :: constraint_check

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
    rhs_loc = 0.
    allocate(local_solver_matrix(mdim*uloc+dloc+n_constraints,&
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
    call get_local_solver(local_solver_matrix,U,X,down,D,f,ele,&
         & g,dt,theta,D0,have_constraint,&
         & local_solver_rhs=local_solver_rhs,projection=projection,&
         & poisson=poisson)

    !Construct the rhs sources for U from lambda
    call assemble_rhs_ele(Rhs_loc,D,U,X,ele,D_rhs,U_rhs,U_rhs_local)
    if(.not.(present(d_rhs).or.present(u_rhs)))then
       assert(.not.present_and_true(projection))
       assert(.not.present_and_true(poisson))
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

    call solve(local_solver_matrix,Rhs_loc)
    rhs_loc = matmul(local_solver_matrix,rhs_loc)
    call solve(local_solver_matrix,Rhs_loc)

    do dim1 = 1, mdim
       U_solved(dim1,:) = rhs_loc(u_start(dim1):u_end(dim1))
       if(.not.(present_and_true(poisson))) then
          if(present(U_out)) then
             call set(U_out,dim1,ele_nodes(u,ele),u_solved(dim1,:))
          else
             call set(U,dim1,ele_nodes(u,ele),u_solved(dim1,:))
          end if
       end if
    end do
    
    D_solved = rhs_loc(d_start:d_end)
    if(.not.(present_and_true(projection))) then
       if(present(D_out)) then
          call set(D_out,ele_nodes(d,ele),D_solved)
       else
          call set(D,ele_nodes(d,ele),D_solved)
       end if
    end if

    !check that the constraints are satisfied
    if(have_constraint) then
       constraints => U%mesh%shape%constraints
       do i1 = 1, constraints%n_constraints
          constraint_check = 0.
          do dim1 = 1, mdim
             constraint_check = constraint_check + &
                  & sum(U_solved(dim1,:)*constraints%orthogonal(i1,:,dim1))
          end do
          if(abs(constraint_check)>1.0e-8) then
             ewrite(1,*) 'Constraint check', constraint_check
             FLAbort('Constraint not enforced')
          end if
       end do       
    end if

  end subroutine reconstruct_U_d_ele

  subroutine get_local_solver(local_solver_matrix,U,X,down,D,f,ele,&
       & g,dt,theta,D0,have_constraint, &
       & local_solver_rhs,projection,poisson)
    !Subroutine to get the matrix and rhs for obtaining U and D within
    !element ele from the lagrange multipliers on the boundaries.
    !This matrix-vector system is referred to as the "local solver" in the 
    !literature e.g.
    !Cockburn et al, Unified hybridization of discontinuous Galerkin, mixed
    ! and continuous Galerkin methods for second order elliptic problems,
    ! SIAM J. Numer. Anal., 2009
    implicit none
    !If projection is present and true, set dt to zero and just project U 
    !into div-conforming space
    real, intent(in) :: g,dt,theta,D0
    type(vector_field), intent(in) :: U,X,down
    type(scalar_field), intent(in) :: D,f
    integer, intent(in) :: ele
    real, dimension(:,:)&
         &, intent(inout) :: local_solver_matrix
    real, dimension(:,:)&
         &, intent(inout), optional :: local_solver_rhs
    logical, intent(in), optional :: projection, poisson
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
    real :: l_dt

    if(present_and_true(projection)) then
       l_dt = 0.
    else
       l_dt = dt
    end if

    mdim = mesh_dim(U)
    uloc = ele_loc(U,ele)
    dloc = ele_loc(d,ele)

    d_start = uloc*2 + 1
    d_end   = uloc*2+dloc
    do dim1 = 1, mdim
       u_start(dim1) = uloc*(dim1-1)+1
       u_end(dim1) = uloc*dim1
    end do

    local_solver_matrix = 0.
    if(present(local_solver_rhs)) then
       local_solver_rhs = 0.
    end if
    
    u_shape=ele_shape(u, ele)
    D_shape=ele_shape(d, ele)
    D_ele => ele_nodes(D, ele)
    U_ele => ele_nodes(U, ele)
    
    if(present_and_true(projection)) then
       f_gi = 0.
    else
       f_gi = ele_val_at_quad(f,ele)
    end if
    up_gi = -ele_val_at_quad(down,ele)

    call get_up_gi(X,ele,up_gi)

    !J, detJ is needed for Piola transform
    !detwei is needed for pressure mass matrix
    call compute_jacobian(X, ele, J=J, detwei=detwei, detJ=detJ)

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
    local_solver_matrix(d_start:d_end,d_start:d_end)=&
         &shape_shape(d_shape,d_shape,detwei)
       if(present(local_solver_rhs)) then
          local_solver_rhs(d_start:d_end,d_start:d_end) = &
               shape_shape(d_shape,d_shape,detwei)
       end if
    !divergence matrix (done in local coordinates)
    l_div_mat = dshape_shape(u_shape%dn,d_shape,&
         &D_shape%quadrature%weight)
    do dim1 = 1, mdim
       !pressure gradient term [integrated by parts so minus sign]
       local_solver_matrix(u_start(dim1):u_end(dim1),d_start:d_end)=&
            & -g*l_dt*theta*l_div_mat(dim1,:,:)
       if(present(local_solver_rhs)) then
          local_solver_rhs(u_start(dim1):u_end(dim1),d_start:d_end)=&
               & -g*(theta-1.0)*l_dt*l_div_mat(dim1,:,:)
       end if
       !divergence continuity term
       local_solver_matrix(d_start:d_end,u_start(dim1):u_end(dim1))=&
            & d0*l_dt*theta*transpose(l_div_mat(dim1,:,:))
       if(present(local_solver_rhs)) then
          local_solver_rhs(d_start:d_end,u_start(dim1):u_end(dim1))=&
               & d0*(theta-1.0)*l_dt*transpose(l_div_mat(dim1,:,:))
       end if
    end do
    !velocity mass matrix and Coriolis matrix (done in local coordinates)
    l_u_mat = shape_shape_tensor(u_shape, u_shape, &
         u_shape%quadrature%weight, Metric+l_dt*theta*Metricf)
    do dim1 = 1, mdim
       do dim2 = 1, mdim
          local_solver_matrix(u_start(dim1):u_end(dim1),&
               u_start(dim2):u_end(dim2))=&
               & l_u_mat(dim1,dim2,:,:)
       end do
    end do

    if(present(local_solver_rhs)) then
       l_u_mat = shape_shape_tensor(u_shape, u_shape, &
            u_shape%quadrature%weight, &
            Metric+(theta-1.0)*l_dt*Metricf)
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
             local_solver_matrix(d_end+i1,u_start(dim1):u_end(dim1))=&
                  &constraints%orthogonal(i1,:,dim1)
             local_solver_matrix(u_start(dim1):u_end(dim1),d_end+i1)=&
                  &constraints%orthogonal(i1,:,dim1)
          end do
       end do
    end if
        
  end subroutine get_local_solver

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
    real :: norm

    call compute_jacobian(X, ele, J)

    select case(mesh_dim(X)) 
    case (2)
       !Coriolis only makes sense for 2d surfaces embedded in 3d
       do gi = 1, ele_ngi(X,ele)
          normal_gi(:,gi) = cross_product(J(1,:,gi),J(2,:,gi))
          norm = sqrt(sum(normal_gi(:,gi)**2))
          normal_gi(:,gi) = normal_gi(:,gi)/norm
       end do
       do gi = 1, ele_ngi(X,ele)
          orientation_gi(gi) = dot_product(normal_gi(:,gi),up_gi(:,gi))
       end do
       if(any(abs(orientation_gi-orientation_gi(1))>1.0e-8)) then
          FLAbort('Nasty geometry problem')
       end if
       do gi = 1, ele_ngi(X,ele)
          up_gi(:,gi) = normal_gi(:,gi)*orientation_gi(gi)
       end do
       if(present(orientation)) then
          if(orientation_gi(1)>0.0) then
             orientation = 1
          else
             orientation = -1
          end if
       end if
    case default
       FLAbort('not implemented')
    end select
  end subroutine get_up_gi

  function get_orientation(X_val, up) result (orientation)
    !function compares the orientation of the element with the 
    !up direction
    implicit none
    real, dimension(:,:), intent(in) :: X_val           !(dim,loc)
    real, dimension(:,:), intent(in) :: up
    integer :: orientation
    !
    real, dimension(size(X_val,1)) :: t1,t2
    real, dimension(size(X_val,1)) :: crossprod
    integer :: gi
    ! if elements are triangles:
    if(size(X_val,2)==3) then
       t1 = X_val(:,2)-X_val(:,1)
       t2 = X_val(:,3)-X_val(:,1)
       crossprod = cross_product(t1,t2)
       if(dot_product(crossprod,up(:,1))>0.0) then
          do gi = 1, size(up,2)
             if(dot_product(crossprod,up(:,gi))<0.0) then
                FLAbort('Something nasty with down direction')
             end if             
          end do
          orientation = 1
       else
          do gi = 1, size(up,2)
             if(dot_product(crossprod,up(:,gi))>0.0) then
                FLAbort('Something nasty with down direction')
             end if             
          end do
          orientation = -1
       end if
    else
       FLAbort('Haven''t sorted out quads yet.')
    end if
  end function get_orientation

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
    call compute_jacobian(X, ele, J=J, detwei=detwei,detJ=detJ)

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

  subroutine check_continuity_local_ele(U,ele)
    implicit none
    type(vector_field), intent(in) :: U
    integer, intent(in) :: ele
    !
    integer, dimension(:), pointer :: neigh
    integer :: ni,face,ele2,face2

    neigh => ele_neigh(U,ele)
    do ni = 1, size(neigh)
       ele2 = neigh(ni)
       face = ele_face(U,ele,ele2)
       if(ele2>0) then
          face2 = ele_face(U,ele2,ele)
       else
          face2 = -1
       end if
       call check_continuity_local_face(U,ele,ele2,face,face2)
    end do
  end subroutine check_continuity_local_ele

  subroutine check_continuity_local_face(U,ele,ele2,face,face2)
    implicit none
    integer, intent(in) :: face, face2,ele,ele2
    type(vector_field), intent(in) :: U
    !
    real, dimension(U%dim, face_ngi(U, face)) :: n1,n2,u1,u2
    real :: weight, jump

    !Get normal in local coordinates
    call get_local_normal(n1,weight,U,local_face_number(U%mesh,face))
    call get_local_normal(n2,weight,U,local_face_number(U%mesh,face2))
    u1 = face_val_at_quad(U,face)
    u2 = face_val_at_quad(U,face2)
    jump = maxval(abs(sum(u1*n1+u2*n2,1)))
    ewrite(1,*) jump
    assert(jump<1.0e-8)

  end subroutine check_continuity_local_face

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
    real, dimension(X%dim, face_ngi(U_cart, face)) :: u1,u2,x1
    real, dimension(face_ngi(U_cart, face)) :: jump_at_quad
    integer :: dim1
    !
    u1 = face_val_at_quad(U_cart,face)
    x1 = face_val_at_quad(X,face)
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
       ewrite(1,*) 'Jump at quadrature face, face2 =', jump_at_quad
       ewrite(1,*) 'ELE = ',ele,ele2
       do dim1 = 1, X%dim
          ewrite(1,*) 'normal',dim1,n1(dim1,:)
          ewrite(1,*) 'normal',dim1,n2(dim1,:)
          ewrite(1,*) 'X',dim1,x1(dim1,:)
       end do
       ewrite(1,*) 'n cpt1',sum(n1*u1,1)
       ewrite(1,*) 'n cpt2',sum(n2*u2,1)
       ewrite(1,*) jump_at_quad/max(maxval(abs(u1)),maxval(abs(u2)))
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

    call compute_jacobian(X, ele, J=J, detwei=detwei)
    call compute_jacobian(X,face, J=J_f, detwei=detwei_f, facet=.true.)

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
    call compute_jacobian(X, ele, J=J, detwei=detwei)
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

  subroutine assemble_rhs_ele(Rhs_loc,D,U,X,ele,D_rhs,U_rhs,u_rhs_local)
    implicit none
    integer, intent(in) :: ele
    type(scalar_field), intent(in), optional, target :: D_rhs
    type(vector_field), intent(in), optional, target :: U_rhs
    type(vector_field), intent(in) :: X,U
    type(scalar_field), intent(in) :: D
    real, dimension(:), &
         &intent(inout) :: Rhs_loc
    logical, intent(in), optional :: u_rhs_local
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
    real, dimension(mesh_dim(U), mesh_dim(U)) :: Metric

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
    
    if(.not.(present(D_rhs).or.present(u_rhs))) then
       !We are in timestepping mode.
       rhs_loc(d_start:d_end) = ele_val(D,ele)
       u_rhs_loc = ele_val(U,ele)
       do dim1 = 1, mdim
          rhs_loc(u_start(dim1):u_end(dim1)) = &
               & U_rhs_loc(dim1,:)
       end do
    else
       have_d_rhs = present(d_rhs)
       have_u_rhs = present(u_rhs)
       if(have_d_rhs) l_d_rhs => d_rhs
       if(have_u_rhs) l_u_rhs => u_rhs
       
       call compute_jacobian(X, ele, J=J, detJ=detJ,detwei=detwei)
       
       Rhs_loc = 0.
       if(have_d_rhs) then
          Rhs_loc(d_start:d_end) = shape_rhs(ele_shape(D,ele),&
               &ele_val_at_quad(l_D_rhs,ele)*detwei)
       end if
       if(have_u_rhs) then
          if(present_and_true(u_rhs_local)) then
             u_local_quad = ele_val_at_quad(l_u_rhs,ele)
             do gi=1,ele_ngi(U,ele)
                Metric=matmul(J(:,:,gi), transpose(J(:,:,gi)))&
                     &/detJ(gi)
                u_local_quad(:,gi) = matmul(Metric,u_local_quad(:,gi))
             end do
          else
             allocate(u_cart_quad(l_U_rhs%dim,ele_ngi(X,ele)))
             u_cart_quad = ele_val_at_quad(l_u_rhs,ele)
             do gi = 1, ele_ngi(D,ele)
                !Don't divide by detJ as we can use weight instead of detwei
                u_local_quad(:,gi) = matmul(J(:,:,gi)&
                     &,u_cart_quad(:,gi))
             end do
          end if
          U_rhs_loc = shape_vector_rhs(u_shape,&
               u_local_quad,u_shape%quadrature%weight)
          do dim1 = 1, mdim
             Rhs_loc(u_start(dim1):u_end(dim1)) = &
                  & U_rhs_loc(dim1,:)
          end do
       end if
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
    call compute_jacobian(X, ele, J=J, detJ=detJ,detwei=detwei)
    
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

  subroutine compute_energy_hybridized(state,energy)
    implicit none
    type(state_type), intent(inout) :: state
    real, intent(inout) :: energy
    !
    type(scalar_field), pointer :: D
    type(vector_field), pointer :: u,X
    integer :: ele
    real :: old_energy,g,d0

    !get parameters
    call get_option("/physical_parameters/gravity/magnitude", g)
    call get_option("/material_phase::Fluid/scalar_field::LayerThickness/p&
         &rognostic/mean_layer_thickness",D0)

    U=>extract_vector_field(state, "Velocity")
    D => extract_scalar_field(state, "LayerThickness")
    X=>extract_vector_field(state, "Coordinate")

    old_energy = energy
    energy = 0.

    do ele = 1, element_count(X)
       call compute_energy_ele(energy,U,D,X,D0,g,ele)
    end do

    ewrite(1,*) 'Energy:= ', energy
    ewrite(1,*) 'Change in energy:= ', energy-old_energy

  end subroutine compute_energy_hybridized

  subroutine compute_energy_ele(energy,U,D,X,D0,g,ele)
    implicit none
    real, intent(inout) :: energy
    type(vector_field), intent(in) :: U,X
    type(scalar_field), intent(in) :: D
    integer, intent(in) :: ele
    real, intent(in) :: D0,g
    !
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(ele_ngi(D,ele)) :: detwei
    type(element_type) :: d_shape, u_shape, x_shape
    real, dimension(ele_loc(U,ele),ele_loc(U,ele)) :: u_mass
    real, dimension(ele_loc(D,ele),ele_loc(D,ele)) :: d_mass
    real, dimension(u%dim,ele_loc(U,ele)) :: U_val
    real, dimension(ele_loc(D,ele)) :: D_val
    real, dimension(X%dim,ele_loc(X,ele)) :: X_val
    integer :: dim1

    U_val = ele_val(U,ele)
    D_val = ele_val(D,ele)
    X_val = ele_val(X,ele)

    u_shape = ele_shape(u,ele)
    d_shape = ele_shape(d,ele)
    x_shape = ele_shape(X,ele)
    call compute_jacobian(X, ele, J=J, detwei=detwei)

    u_mass = shape_shape(u_shape,u_shape,detwei)
    d_mass = shape_shape(d_shape,d_shape,detwei)

    !kinetic energy
    do dim1 = 1, u%dim
       energy = energy + D0*dot_product(U_val(dim1,:),&
            &matmul(u_mass,U_val(dim1,:)))
    end do

    energy = energy + g*dot_product(D_val,&
         &matmul(D_mass,D_val))

  end subroutine compute_energy_ele

  subroutine set_velocity_from_geostrophic_balance_hybridized(&
       &state)
    implicit none
    type(state_type), intent(inout) :: state
    !
    type(scalar_field), pointer :: D,psi,f
    type(scalar_field) :: D_rhs
    type(vector_field), pointer :: U_local,down,X, U_cart
    type(vector_field) :: Coriolis_term, Balance_eqn, tmp_field
    integer :: ele,dim1
    real :: g
    logical :: elliptic_method

    D=>extract_scalar_field(state, "LayerThickness")
    psi=>extract_scalar_field(state, "Streamfunction")
    f=>extract_scalar_field(state, "Coriolis")
    U_local=>extract_vector_field(state, "LocalVelocity")
    U_cart=>extract_vector_field(state, "Velocity")
    X=>extract_vector_field(state, "Coordinate")
    down=>extract_vector_field(state, "GravityDirection")
    call get_option("/physical_parameters/gravity/magnitude", g)
    call allocate(tmp_field,mesh_dim(U_local), U_local%mesh, "tmp_field")
       call allocate(D_rhs,D%mesh,'BalancedSolverRHS')
    call allocate(Coriolis_term,mesh_dim(U_local),&
         U_local%mesh,"CoriolisTerm")
    call allocate(balance_eqn,mesh_dim(D),u_local%mesh,'BalancedEquation')

    !STAGE 1: Set velocity from streamfunction
    do ele = 1, element_count(D)
       call set_local_velocity_from_streamfunction_ele(&
            &U_local,psi,down,X,ele)
    end do

    !STAGE 1a: verify that velocity projects is div-conforming
    call project_local_to_cartesian(X,U_local,U_cart)
    do ele = 1, ele_count(U_local)
       call check_continuity_ele(U_cart,X,ele)
    end do

    !Stage 1b: verify that projection is idempotent
    ewrite(1,*) 'CHECKING CONTINUOUS', maxval(abs(u_local%val))
    call solve_hybridized_helmholtz(state,U_Rhs=U_local,&
         &U_out=tmp_field,&
         &compute_cartesian=.true.,&
         &check_continuity=.true.,projection=.true.,&
         &u_rhs_local=.true.)!verified that projection is idempotent
    assert(maxval(abs(U_local%val-tmp_field%val))<1.0e-8)

    elliptic_method = .false.

    if(elliptic_method) then

       !STAGE 2: Construct Coriolis term
       call zero(Coriolis_term)
       do ele = 1, element_count(D)
          call set_coriolis_term_ele(Coriolis_term,f,down,U_local,X,ele)
       end do!checked!signs checked

       !STAGE 3: Project Coriolis term into div-conforming space

       !debugging bits - checking if it works with cartesian instead
       call project_local_to_cartesian(X,Coriolis_Term,U_cart)
       call solve_hybridized_helmholtz(state,U_Rhs=U_cart,&
            &U_out=tmp_field,&
            &compute_cartesian=.true.,output_dense=.false.,&
            &check_continuity=.true.,projection=.true.,&
            &u_rhs_local=.false.)

       ewrite(0,*) 'REMEMBER TO REMOVE DEBUGGING TESTS'
       call solve_hybridized_helmholtz(state,U_Rhs=Coriolis_term,&
            &U_out=tmp_field,&
            &compute_cartesian=.true.,&
            &check_continuity=.true.,projection=.true.,&
            &u_rhs_local=.true.)!verified that projection is idempotent

       !STAGE 4: Construct the RHS for the balanced layer depth equation
       call zero(D_rhs)
       do ele = 1, element_count(D)
          call set_geostrophic_balance_rhs_ele(D_rhs,Coriolis_term,ele)
       end do

       !STAGE 5: Solve Poisson equation for the balanced layer depth
       ewrite(0,*) 'REMEMBER ABOUT SETTING MEAN VALUE'
       ewrite(0,*) trim(u_cart%option_path)

       call solve_hybridized_helmholtz(state,D_rhs=D_rhs,&
            &compute_cartesian=.false.,&
            &check_continuity=.false.,Poisson=.true.,&
            &solver_option_path=trim(u_cart%option_path)//'/prognostic/initial_condition::WholeMesh/balanced')

       !STAGE 6: Check if we have a balanced solution
       !Can be done by projecting balance equation into div-conforming space
       !and checking that it is equal to zero
       !STAGE 6a: Project balance equation into DG space
       do ele = 1, element_count(D)
          call set_pressure_force_ele(balance_eqn,D,X,g,ele)
       end do
       call addto(balance_eqn,coriolis_term)

       !STAGE 6b: Project balance equation into div-conforming space
       call solve_hybridized_helmholtz(state,U_Rhs=balance_eqn,&
            &U_out=balance_eqn,&
            &compute_cartesian=.true.,&
            &check_continuity=.true.,projection=.true.,&
            &u_rhs_local=.true.)

       do dim1 = 1, mesh_dim(D)
          ewrite(1,*) 'Balance equation', maxval(abs(balance_eqn%val(dim1,:)))
          assert(maxval(abs(balance_eqn%val(dim1,:)))<1.0e-8)
       end do

    else
       !Project the streamfunction into pressure space
       do ele = 1, element_count(D)
          call project_streamfunction_for_balance_ele(D,psi,X,f,g,ele)
       end do
       ewrite(1,*) maxval(abs(D%val))

       !debugging tests
       call zero(Coriolis_term)
       do ele = 1, element_count(D)
          call set_coriolis_term_ele(Coriolis_term,f,down,U_local,X,ele)
       end do
       call zero(balance_eqn)
       do ele = 1, element_count(D)
          call set_pressure_force_ele(balance_eqn,D,X,g,ele)
       end do
       call addto(balance_eqn,coriolis_term,scale=1.0)
       ewrite(1,*) 'CJC b4',maxval(abs(balance_eqn%val)),&
            & maxval(abs(coriolis_term%val))
       !Project balance equation into div-conforming space
       call solve_hybridized_helmholtz(state,U_Rhs=balance_eqn,&
            &U_out=balance_eqn,&
            &compute_cartesian=.true.,&
            &check_continuity=.true.,projection=.true.,&
            &u_rhs_local=.true.)

        do dim1 = 1, mesh_dim(D)
           ewrite(1,*) 'Balance equation', maxval(abs(balance_eqn%val(dim1,:)))
           assert(maxval(abs(balance_eqn%val(dim1,:)))<1.0e-8)
        end do
    end if
    !Clean up after yourself
    call deallocate(Coriolis_term)
    call deallocate(D_rhs)
    call deallocate(balance_eqn)
    call deallocate(tmp_field)

  end subroutine set_velocity_from_geostrophic_balance_hybridized

  subroutine project_streamfunction_for_balance_ele(D,psi,X,f,g,ele)
    implicit none
    type(scalar_field), intent(in) :: psi,f
    type(scalar_field), intent(inout) :: D
    type(vector_field), intent(in) :: X
    integer, intent(in) :: ele
    real, intent(in) :: g
    !
    real, dimension(ele_loc(d,ele),ele_loc(d,ele)) :: d_mass
    real, dimension(ele_loc(d,ele)) :: d_rhs
    type(element_type) :: psi_shape, d_shape
    real, dimension(ele_ngi(d,ele)) :: detwei, psi_quad,f_gi
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J

    f_gi = ele_val_at_quad(f,ele)
    psi_shape = ele_shape(psi,ele)
    d_shape = ele_shape(d,ele)
    psi_quad = ele_val_at_quad(psi,ele)

    call compute_jacobian(X, ele, J=J, detwei=detwei)

    d_rhs = shape_rhs(d_shape,detwei*psi_quad*f_gi/g)
    d_mass = shape_shape(d_shape,d_shape,detwei)
    call solve(d_mass,d_rhs)
    call set(D,ele_nodes(D,ele),d_rhs)

  end subroutine project_streamfunction_for_balance_ele
  
  subroutine set_pressure_force_ele(force,D,X,g,ele)
    implicit none
    type(vector_field), intent(inout) :: force
    type(scalar_field), intent(in) :: D
    type(vector_field), intent(in) :: X
    real, intent(in) :: g
    integer, intent(in) :: ele
    !
    real, dimension(ele_ngi(D,ele)) :: D_gi
    real, dimension(mesh_dim(D),ele_loc(force,ele)) :: &
         & rhs_loc
    real, dimension(ele_loc(force,ele),ele_loc(force,ele)) :: &
         & l_mass_mat    
    integer :: dim1,dim2,gi,uloc
    real, dimension(mesh_dim(force), X%dim, ele_ngi(force,ele)) :: J
    real, dimension(ele_ngi(force,ele)) :: detJ
    real, dimension(mesh_dim(force),mesh_dim(force),ele_ngi(force,ele))::&
         &Metric
    real, dimension(mesh_dim(force)*ele_loc(force,ele),&
         mesh_dim(force)*ele_loc(force,ele)) :: l_u_mat
    real, dimension(mesh_dim(force)*ele_loc(force,ele)) :: force_rhs
    type(element_type) :: force_shape

    uloc = ele_loc(force,ele)
    force_shape = ele_shape(force,ele)
    call compute_jacobian(X, ele, J=J, detJ=detJ)
    do gi=1,ele_ngi(force,ele)
       Metric(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)
    end do

    D_gi = ele_val_at_quad(D,ele)
    rhs_loc = -g*dshape_rhs(force%mesh%shape%dn,&
         D_gi*D%mesh%shape%quadrature%weight)
    l_mass_mat = shape_shape(ele_shape(force,ele),ele_shape(force,ele),&
         &force%mesh%shape%quadrature%weight)
    do dim1 = 1, mesh_dim(force)
       force_rhs((dim1-1)*uloc+1:dim1*uloc) = rhs_loc(dim1,:)
    end do
    do dim1 = 1, mesh_dim(force)
       do dim2 = 1, mesh_dim(force)
          l_u_mat((dim1-1)*uloc+1:dim1*uloc,&
               &  (dim2-1)*uloc+1:dim2*uloc ) = &
               & shape_shape(force_shape,force_shape,&
               & force_shape%quadrature%weight*Metric(dim1,dim2,:))
       end do
    end do

    call solve(l_u_mat,force_rhs)
    do dim1= 1, mesh_dim(force)
       call set(force,dim1,ele_nodes(force,ele),&
            &force_rhs((dim1-1)*uloc+1:dim1*uloc))
    end do

  end subroutine set_pressure_force_ele

  subroutine set_geostrophic_balance_rhs_ele(D_rhs,Coriolis_term,ele)
    implicit none
    type(scalar_field), intent(inout) :: D_rhs
    type(vector_field), intent(in) :: Coriolis_term
    integer, intent(in) :: ele
    !
    real, dimension(mesh_dim(Coriolis_term),ele_loc(Coriolis_term,ele)) :: &
         & Coriolis_loc
    real, dimension(ele_ngi(D_rhs,ele)) :: div_gi
    real, dimension(ele_loc(D_rhs,ele)) :: D_rhs_loc
    real, dimension(ele_loc(D_rhs,ele),ele_loc(D_rhs,ele)) :: d_mass
    integer :: dim1 
    type(element_type) :: U_shape, D_shape
    real :: g

    !Computes the divergence of projected Coriolis term
    !Can be done locally since d commutes with pullback

    call get_option("/physical_parameters/gravity/magnitude", g)
    
    U_shape = ele_shape(Coriolis_term,ele)
    D_shape = ele_shape(D_rhs,ele)
    Coriolis_loc = ele_val(Coriolis_term,ele)

    div_gi = 0.
    do dim1 = 1, mesh_dim(Coriolis_term)
       div_gi = div_gi + matmul(transpose(U_shape%dn(:,:,dim1)),&
            &Coriolis_loc(dim1,:))
    end do
    D_rhs_loc = shape_rhs(D_shape,div_gi*U_shape%quadrature%weight)
    d_mass = shape_shape(d_shape,d_shape,d_shape%quadrature%weight)
    call solve(d_mass,D_rhs_loc)
    call set(D_rhs,ele_nodes(D_rhs,ele),-D_rhs_loc/g)
  end subroutine set_geostrophic_balance_rhs_ele

  subroutine set_coriolis_term_ele(Coriolis_term,f,down,U_local,X,ele)
    implicit none
    type(vector_field), intent(inout) :: Coriolis_term
    type(vector_field), intent(in) :: U_local,X,down
    type(scalar_field), intent(in) :: f
    integer, intent(in) :: ele
    !
    real, dimension(ele_ngi(x,ele)) :: f_gi
    real, dimension(X%dim, ele_ngi(X,ele)) :: up_gi
    real, dimension(mesh_dim(U_local),mesh_dim(U_local),ele_ngi(U_local,ele))::&
         &Metric, Metricf
    real, dimension(mesh_dim(U_local), X%dim, ele_ngi(U_local,ele)) :: J
    real, dimension(ele_ngi(U_local,ele)) :: detJ
    real, dimension(X%dim, X%dim, ele_ngi(U_local,ele)) :: rot
    real, dimension(mesh_dim(U_local)*ele_loc(U_local,ele),&
         mesh_dim(U_local)*ele_loc(U_local,ele)) :: l_u_mat
    real, dimension(mesh_dim(U_local)*ele_loc(U_local,ele)) :: coriolis_rhs
    real, dimension(mesh_dim(U_local), ele_ngi(U_local,ele)) :: U_gi
    real, dimension(mesh_dim(U_local), ele_ngi(U_local,ele)) :: coriolis_gi
    real, dimension(X%dim) :: up_vec
    integer :: dim1, dim2,uloc,gi
    type(element_type) :: u_shape
    
    uloc = ele_loc(u_local,ele)
    u_shape = ele_shape(u_local,ele)

    u_gi = ele_val_at_quad(u_local,ele)
    f_gi = ele_val_at_quad(f,ele)
    up_gi = -ele_val_at_quad(down,ele)

    call get_up_gi(X,ele,up_gi)

    coriolis_rhs = 0.
    l_u_mat = 0.
    !metrics for velocity mass and coriolis matrices
    call compute_jacobian(X, ele, J=J, detJ=detJ)
    do gi=1, ele_ngi(U_local,ele)
       rot(1,:,gi)=(/0.,-up_gi(3,gi),up_gi(2,gi)/)
       rot(2,:,gi)=(/up_gi(3,gi),0.,-up_gi(1,gi)/)
       rot(3,:,gi)=(/-up_gi(2,gi),up_gi(1,gi),0./)
    end do
    do gi=1,ele_ngi(U_local,ele)
       Metric(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)
       Metricf(:,:,gi)=matmul(J(:,:,gi), &
            matmul(f_gi(gi)*rot(:,:,gi), transpose(J(:,:,gi))))/detJ(gi)
       Coriolis_gi(:,gi) = matmul(Metricf(:,:,gi),u_gi(:,gi))
    end do

    !Coriolis term is evaluated in global coordinates [hence presence of
    ! metric terms] and projected into local velocity coordinates
    do dim1 = 1, mesh_dim(U_local)
       do dim2 = 1, mesh_dim(U_local)
          l_u_mat((dim1-1)*uloc+1:dim1*uloc,&
               &  (dim2-1)*uloc+1:dim2*uloc ) = &
               & shape_shape(u_shape,u_shape,&
               & u_shape%quadrature%weight*Metric(dim1,dim2,:))
       end do
       coriolis_rhs((dim1-1)*uloc+1:dim1*uloc) = &
            & shape_rhs(u_shape,Coriolis_gi(dim1,:)*u_shape%quadrature%weight)
    end do
    call solve(l_u_mat,coriolis_rhs)
    do dim1= 1, mesh_dim(coriolis_term)
       call set(coriolis_term,dim1,ele_nodes(Coriolis_term,ele),&
            &Coriolis_rhs((dim1-1)*uloc+1:dim1*uloc))
    end do
  end subroutine set_coriolis_term_ele
  
  subroutine set_local_velocity_from_streamfunction_ele(&
       &U_local,psi,down,X,ele)
    implicit none
    type(vector_field), intent(inout) :: U_local
    type(vector_field), intent(in) :: down,X
    type(scalar_field), intent(in) :: psi
    integer, intent(in) :: ele
    !
    real, dimension(ele_loc(psi,ele)) :: psi_loc
    real, dimension(mesh_dim(psi),ele_ngi(psi,ele)) :: dpsi_gi
    real, dimension(ele_ngi(psi,ele)) :: div_gi
    real, dimension(mesh_dim(U_local),ele_loc(U_local,ele)) :: U_loc
    real, dimension(ele_loc(U_local,ele),ele_loc(U_local,ele)) :: &
         & l_mass_mat
    type(element_type) :: u_shape, psi_shape
    integer :: dim1,gi,uloc
    real, dimension(X%dim, ele_ngi(X,ele)) :: up_gi
    integer :: orientation

    uloc = ele_loc(U_local,ele)
    u_shape = ele_shape(U_local,ele)
    psi_shape = ele_shape(psi,ele)
    up_gi = -ele_val_at_quad(down,ele)
    call get_up_gi(X,ele,up_gi,orientation)

    !We can do everything in local coordinates since d commutes with pullback
    !usual tricks: dpsi lives in the U space so we can do projection

    l_mass_mat = shape_shape(u_shape,u_shape,U_shape%quadrature%weight)

    !Streamfunction at node values
    psi_loc = ele_val(psi,ele)
    !Skew gradient of streamfunction at quadrature points
    select case(mesh_dim(psi))
    case (2)
       forall(gi=1:ele_ngi(psi,ele))
          dpsi_gi(1,gi) = -sum(psi_loc*psi_shape%dn(:,gi,2))
          dpsi_gi(2,gi) =  sum(psi_loc*psi_shape%dn(:,gi,1))
       end forall
    case default
       FLAbort('Exterior derivative not implemented for given mesh dimension')
    end select
    dpsi_gi = orientation*dpsi_gi
    U_loc = shape_vector_rhs(u_shape,dpsi_gi,U_shape%quadrature%weight)

    do dim1 = 1, U_local%dim
       call solve(l_mass_mat,U_loc(dim1,:))
       call set(U_local,dim1,ele_nodes(U_local,ele),&
            u_loc(dim1,:))
    end do

    !verify divergence-free-ness
    div_gi = 0.
    do gi = 1, ele_ngi(psi,ele)
       do dim1 = 1, mesh_dim(psi)
          div_gi(gi) = div_gi(gi) + sum(u_shape%dn(:,gi,dim1)*u_loc(dim1,:))
       end do
    end do
    assert(maxval(abs(div_gi))<1.0e-8)
  end subroutine set_local_velocity_from_streamfunction_ele

  subroutine project_to_constrained_space(state,v_field)
    !wrapper for projecting vector field to constrained space
    implicit none
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field

    call solve_hybridized_helmholtz(state,&
         &U_rhs=v_field,&
         &compute_cartesian=.true.,&
         &check_continuity=.true.,projection=.true.,&
         &u_rhs_local=.true.)

  end subroutine project_to_constrained_space

end module hybridized_helmholtz
