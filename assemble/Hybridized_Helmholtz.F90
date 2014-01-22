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
    use global_parameters, only: option_path_len, PYTHON_FUNC_LEN
    use vector_tools, only: solve
    use manifold_tools
    use vtk_interfaces
    implicit none

    !Variables belonging to the Newton solver
    logical, private :: newton_initialised = .false.
    !Cached Newton solver matrices
    type(csr_matrix), private :: Newton_lambda_mat
    type(block_csr_matrix), private :: Newton_continuity_mat
    type(real_matrix), pointer, dimension(:), private &
         &:: Newton_local_solver_cache=>null()
    type(real_matrix), pointer, dimension(:), private &
         &:: Newton_local_solver_rhs_cache=>null()

contains 

  subroutine solve_hybridised_timestep_residual(state,newU,newD,&
       UResidual, DResidual)

    ! Subroutine to apply one Newton iteration to hybridised shallow
    ! water equations
    ! Written to cache all required matrices
    implicit none
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: newU !U at next timestep
    type(scalar_field), intent(inout) :: newD !D at next timestep
    type(vector_field), intent(in), optional :: UResidual
    type(scalar_field), intent(in), optional :: DResidual
    !
    type(vector_field), pointer :: X, U, down, U_cart, U_old
    type(scalar_field), pointer :: D,f, D_old
    type(scalar_field) :: lambda, D_res
    type(vector_field) :: U_res, U_res2
    type(scalar_field), target :: lambda_rhs, u_cpt
    type(csr_sparsity) :: lambda_sparsity, continuity_sparsity
    type(csr_matrix) :: continuity_block_mat
    type(mesh_type), pointer :: lambda_mesh
    real :: D0, dt, g, theta, tolerance
    integer :: ele,i1, stat, dim1, dloc, uloc,mdim,n_constraints
    logical :: have_constraint
    character(len=OPTION_PATH_LEN) :: constraint_option_string

    ewrite(1,*) 'solve_hybridised_timestep_residual(state)'

    ewrite(1,*) 'TIMESTEPPING LAMBDA'

    !get parameters
    call get_option("/physical_parameters/gravity/magnitude", g)
    call get_option("/timestepping/theta",theta)
    call get_option("/material_phase::Fluid/scalar_field::LayerThickness/&
         &prognostic/mean_layer_thickness",D0)
    call get_option("/timestepping/timestep", dt)

    !Pull the fields out of state
    D=>extract_scalar_field(state, "LayerThickness")
    f=>extract_scalar_field(state, "Coriolis")
    U=>extract_vector_field(state, "LocalVelocity")
    X=>extract_vector_field(state, "Coordinate")
    down=>extract_vector_field(state, "GravityDirection")
    U_cart => extract_vector_field(state, "Velocity")
    U_old => extract_vector_field(state, "OldLocalVelocity")
    D_old => extract_scalar_field(state, "OldLayerThickness")

    !Allocate local field variables
    lambda_mesh=>extract_mesh(state, "VelocityMeshTrace")
    call allocate(lambda,lambda_mesh,name="LagrangeMultiplier")
    call allocate(lambda_rhs,lambda%mesh,"LambdaRHS")
    call zero(lambda_rhs)
    call allocate(U_res,U%dim,U%mesh,"VelocityResidual")
    call allocate(D_res,D%mesh,"LayerThicknessResidual")
    call allocate(U_res2,U%dim,U%mesh,"VelocityResidual2")

    if(.not.newton_initialised) then
       !construct/extract sparsities
       lambda_sparsity=get_csr_sparsity_firstorder(&
            &state, lambda%mesh, lambda&
            &%mesh)
       continuity_sparsity=get_csr_sparsity_firstorder(state, u%mesh,&
            & lambda%mesh)

       !allocate matrices
       call allocate(Newton_lambda_mat,lambda_sparsity)
       call zero(Newton_lambda_mat)
       call allocate(Newton_continuity_mat,continuity_sparsity,(/U%dim,1/))
       call zero(Newton_continuity_mat)

       mdim = mesh_dim(U)
       have_constraint = &
            &have_option(trim(U%mesh%option_path)//"/from_mesh/constraint_type")

       allocate(newton_local_solver_cache(ele_count(U)))
       allocate(newton_local_solver_rhs_cache(ele_count(U)))
       do ele = 1, ele_count(D)
          uloc = ele_loc(U,ele)
          dloc = ele_loc(d,ele)
          n_constraints = 0
          if(have_constraint) then
             n_constraints = ele_n_constraints(U,ele)
             allocate(newton_local_solver_cache(ele)%ptr(&
                  mdim*uloc+dloc+n_constraints,&
                  mdim*uloc+dloc+n_constraints))
             allocate(newton_local_solver_rhs_cache(ele)%ptr(&
                  mdim*uloc+dloc,mdim*uloc+dloc))
          end if
       end do

       !Assemble matrices
       do ele = 1, ele_count(D)
          call assemble_newton_solver_ele(D,f,U,X,down,lambda_rhs,ele, &
               &g,dt,theta,D0,&
               &lambda_mat=Newton_lambda_mat,&
               &continuity_mat=Newton_continuity_mat,&
               &Local_solver_matrix=Newton_local_solver_cache(ele)%ptr,&
               &Local_solver_rhs=Newton_local_solver_rhs_cache(ele)%ptr)
       end do
       newton_initialised = .true.
    end if

    if(present(DResidual).and.present(UResidual).and..not.&
         have_option('/material_phase::Fluid/vector_field::Velocity/prognostic/wave_equation/fully_coupled/linear_debug')) then
       !Set residuals from nonlinear input
       call set(U_res,UResidual)
       call set(D_res,DResidual)
    else
       !Compute residuals for linear equation
       do ele = 1, ele_count(D)
          call get_linear_residuals_ele(U_res,D_res,&
               U_old,D_old,newU,newD,&
               newton_local_solver_cache(ele)%ptr,&
               newton_local_solver_rhs_cache(ele)%ptr,ele)
       end do
       ewrite(2,*) 'cjc U_res calc', sum(newU%val), maxval(abs(U_res%val))
       ewrite(2,*) 'cjc D_res calc', sum(newD%val), maxval(abs(D_res%val))
    end if

    do ele = 1, ele_count(D)
       call local_solve_residuals_ele(U_res,D_res,&
            newton_local_solver_cache(ele)%ptr,ele)
    end do

    call mult_t(lambda_rhs,Newton_continuity_mat,U_res)
    call scale(lambda_rhs,-1.0)
    ewrite(2,*) 'LAMBDARHS', maxval(abs(lambda_rhs%val))

    call zero(lambda)
    !Solve the equations
    call petsc_solve(lambda,newton_lambda_mat,lambda_rhs,&
         option_path=trim(U_cart%mesh%option_path)//&
         &"/from_mesh/constraint_type")
    ewrite(2,*) 'LAMBDA', maxval(abs(lambda%val))

    !Update new U and new D from lambda
    !Compute residuals 
    if(present(DResidual).and.present(UResidual).and..not.&
         have_option('/material_phase::Fluid/vector_field::Velocity/prognostic/wave_equation/fully_coupled/linear_debug')) then
       call set(U_res,UResidual)
       call set(D_res,DResidual)
    else
       do ele = 1, ele_count(D)
          call get_linear_residuals_ele(U_res,D_res,&
               U_old,D_old,newU,newD,&
               newton_local_solver_cache(ele)%ptr,&
               newton_local_solver_rhs_cache(ele)%ptr,ele)
       end do
    end if
    ewrite(2,*) 'cjc U_res recalc', sum(newU%val), maxval(abs(U_res%val))
    ewrite(2,*) 'cjc D_res recalc', sum(newD%val), maxval(abs(D_res%val))

    ! ( M    C  -L)(u)   (v)
    ! ( -C^T N  0 )(h) = (j)
    ! ( L^T  0  0 )(l)   (0)
    ! 
    ! (u)   (M    C)^{-1}(v)   (M    C)^{-1}(L)
    ! (h) = (-C^T N)     (j) + (-C^T N)     (0)(l)
    ! so
    !        (M    C)^{-1}(L)         (M    C)^{-1}(v)
    ! (L^T 0)(-C^T N)     (0)=-(L^T 0)(-C^T N)     (j)
       
    call mult(U_res2,Newton_continuity_mat,lambda)
    call addto(U_res,U_res2)
    do ele = 1, ele_count(U)
       call local_solve_residuals_ele(U_res,D_res,&
            newton_local_solver_cache(ele)%ptr,ele)
    end do
    call halo_update(U_res)
    call halo_update(D_res)

    !Negative scaling occurs here.
    call scale(U_res,-1.0)
    call scale(D_res,-1.0)

    ewrite(2,*) 'Delta U', maxval(abs(U_res%val))
    ewrite(2,*) 'Delta D', maxval(abs(D_res%val))
    
    call addto(newU,U_res)
    call addto(newD,D_res)

    !Deallocate local variables
    call deallocate(lambda_rhs)
    call deallocate(lambda)
    call deallocate(U_res)
    call deallocate(D_res)
    call deallocate(U_res2)

    ewrite(1,*) 'END solve_hybridised_timestep_residual(state)'
  
  end subroutine solve_hybridised_timestep_residual

  subroutine solve_hybridized_helmholtz(state,D_rhs,U_Rhs,&
       &D_out,U_out,&
       &dt_in,theta_in,&
       &compute_cartesian,&
       &output_dense,&
       &projection,poisson,&
       &u_rhs_local)

    ! Subroutine to solve hybridized helmholtz equation
    ! If D_rhs (scalar pressure field) is present, then solve:
    ! <w,u> + <w,fu^\perp> - g <div w,d> + <<[w],d>> = <w,U_rhs>
    ! <\phi,d> +  <\phi,div u> = <\phi, D_rhs>
    ! <<\gamma, [u]>> = 0
    ! (i.e. for unit testing)
    ! otherwise solve:
    ! <w,u> + dt*theta*<w,fu^\perp> - dt*theta*g <div w,d> + <<[w],l>> = 
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
    real, intent(in), optional :: theta_in,dt_in
    logical, intent(in), optional :: compute_cartesian, &
         &output_dense, projection,poisson
    logical, intent(in), optional :: u_rhs_local !means u_rhs is in local coords
    !
    type(vector_field), pointer :: X, U, down, U_cart
    type(scalar_field), pointer :: D,f, lambda_nc
    type(scalar_field) :: lambda
    type(scalar_field), target :: lambda_rhs, u_cpt
    type(csr_sparsity) :: lambda_sparsity, continuity_sparsity
    type(csr_matrix) :: lambda_mat, continuity_block_mat,continuity_block_mat1
    type(block_csr_matrix) :: continuity_mat
    type(mesh_type), pointer :: lambda_mesh
    real :: D0, dt, g, theta, tolerance
    integer :: ele,i1, stat, dim1
    logical :: l_compute_cartesian, l_output_dense
    real, dimension(:,:), allocatable :: lambda_mat_dense
    character(len=OPTION_PATH_LEN) :: constraint_option_string
    real :: u_max
    logical :: l_projection,l_poisson, pullback

    ewrite(1,*) '  subroutine solve_hybridized_helmholtz('

    l_compute_cartesian = .false.
    if(present(compute_cartesian)) l_compute_cartesian = compute_cartesian
    l_output_dense = .false.
    if(present(output_dense)) l_output_dense = output_dense

    if(present(projection)) then
       l_projection = projection
    else
       l_projection = .false.
    end if
    if(present(poisson)) then
       l_poisson = poisson
    end if
    ewrite(2,*) 'Projection = ', l_projection
    ewrite(2,*) 'Poisson = ', l_poisson
    if(l_poisson.and.l_projection) then
       FLAbort('Can''t do projection and poisson')
    end if

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

    if(present(theta_in)) then
       theta = theta_in
    else
       call get_option("/timestepping/theta",theta)
    end if
    !D0
    call get_option("/material_phase::Fluid/scalar_field::LayerThickness/&
         &prognostic/mean_layer_thickness",D0)
    if(present(dt_in)) then
       dt = dt_in
    else
       call get_option("/timestepping/timestep", dt)
    end if

    pullback = have_option('/material_phase::Fluid/scalar_field::LayerThickn&
         &ess/prognostic/spatial_discretisation/discontinuous_galerkin/wave_&
         &equation/pullback')

    !Assemble matrices
    do ele = 1, ele_count(D)
       call assemble_hybridized_helmholtz_ele(D,f,U,X,down,ele, &
            &g,dt,theta,D0,pullback,lambda_mat=lambda_mat,&
            &lambda_rhs=lambda_rhs,D_rhs=D_rhs,U_rhs=U_rhs,&
            &continuity_mat=continuity_mat,&
            &projection=l_projection,poisson=l_poisson,&
            &u_rhs_local=u_rhs_local)
    end do

    ewrite(2,*) 'LAMBDARHS', maxval(abs(lambda_rhs%val))
    call zero(lambda)
    !Solve the equations
    call petsc_solve(lambda,lambda_mat,lambda_rhs,&
         option_path=trim(U_cart%mesh%option_path)//&
         &"/from_mesh/constraint_type")

    !Reconstruct U and D from lambda
    do ele = 1, ele_count(D)
       call reconstruct_u_d_ele(D,f,U,X,down,ele, &
            &g,dt,theta,D0,pullback,&
            &D_rhs=D_rhs,U_rhs=U_rhs,lambda=lambda,&
            &D_out=D_out,U_out=U_out,&
            &projection=l_projection,poisson=l_poisson,&
            &u_rhs_local=u_rhs_local)
    end do

    if(l_poisson.and.present(D_out)) then
       call check_zero_level(D_out)
    end if
    ewrite(2,*) 'LAMBDA', maxval(abs(lambda%val))

    if(l_output_dense) then
       allocate(lambda_mat_dense(node_count(lambda),node_count(lambda)))
       lambda_mat_dense = dense(lambda_mat)
       ewrite(2,*) '-----------'
       do i1 = 1, node_count(lambda)
          ewrite(2,*) lambda_mat_dense(i1,:)
       end do
       ewrite(2,*) '-----------'
    end if

    if(l_compute_cartesian) then
       U_cart => extract_vector_field(state, "Velocity")
       if(present(U_out)) then
          call project_local_to_cartesian(X,U_out,U_cart)
       else
          call project_local_to_cartesian(X,U,U_cart)
       end if
    end if

    if(have_option('/geometry/mesh::VelocityMesh/check_continuity_matrix')) then
       ewrite(2,*) 'Checking continuity'

       call zero(lambda_rhs)
       u_max = 0.
       do dim1 = 1,U%dim
          if(present(U_out)) then
             u_cpt = extract_scalar_field(U_out,dim1)
          else
             u_cpt = extract_scalar_field(U,dim1)
          end if
          continuity_block_mat = block(continuity_mat,dim1,1)
          call mult_T_addto(lambda_rhs,continuity_block_mat,u_cpt)
          ewrite(2,*) 'U, lambda',&
               &maxval(u_cpt%val), maxval(abs(lambda_rhs%val))
          u_max = u_max + maxval(abs(u_cpt%val))
       end do
       ewrite(2,*)'JUMPS MIN:MAX',minval(lambda_rhs%val),&
            &maxval(lambda_rhs%val), u_max

       call get_option('/geometry/mesh::VelocityMesh/check_continuity_matrix/tolerance',tolerance)
       if(maxval(abs(lambda_rhs%val))/max(1.0,u_max/3.0)>tolerance) then
          ewrite(-1,*) 'value =', maxval(abs(lambda_rhs%val))/max(1.0,u_max/3.0)
          FLExit('Continuity matrix tolerance failure')
       end if
    end if

    if(have_option('/geometry/mesh::VelocityMesh/check_continuity')) then       
       U_cart => extract_vector_field(state, "Velocity")
       if(present(U_out)) then
          call project_local_to_cartesian(X,U_out,U_cart)
       else
          call project_local_to_cartesian(X,U,U_cart)
       end if
       call get_option('/geometry/mesh::VelocityMesh/check_continuity/tolerance', tolerance)
       do ele = 1, ele_count(U)
          call check_continuity_ele(U_cart,X,ele,tolerance)
       end do
    end if
    
    call deallocate(lambda_mat)
    call deallocate(continuity_mat)
    call deallocate(lambda_rhs)
    call deallocate(lambda)

    ewrite(1,*) 'END subroutine solve_hybridized_helmholtz'

  end subroutine solve_hybridized_helmholtz

  subroutine get_linear_residuals_ele(U_res,D_res,U,D,&
       newU,newD,local_solver,local_solver_rhs,ele)
    !Subroutine
    type(vector_field), intent(inout) :: U_res,U,newU
    type(scalar_field), intent(inout) :: D_res,D,newD
    real, dimension(:,:), intent(in) :: local_solver, local_solver_rhs
    integer, intent(in) :: ele
    !
    integer :: uloc, dloc, dim1, d_start, d_end, mdim
    integer, dimension(mesh_dim(U)) :: U_start, U_end
    real, dimension(mesh_dim(U)*ele_loc(U,ele)+ele_loc(D,ele))&
         :: rhs_loc1,rhs_loc2
    real, dimension(mesh_dim(U),ele_loc(U,ele)) :: U_val

    !Calculate indices in a vector containing all the U and D dofs in
    !element ele, First the u1 components, then the u2 components, then the
    !D components are stored.
    uloc = ele_loc(U_res,ele)
    dloc = ele_loc(D_res,ele)
    mdim = mesh_dim(U)
    do dim1 = 1, mdim
       u_start(dim1) = uloc*(dim1-1)+1
       u_end(dim1) = uloc*dim1
    end do
    d_start = uloc*mdim + 1
    d_end   = uloc*mdim+dloc

    U_val = ele_val(U,ele)
    do dim1 = 1, mesh_dim(U)
       rhs_loc1(u_start(dim1):u_end(dim1)) = u_val(dim1,:)
    end do
    rhs_loc1(d_start:d_end) = ele_val(D,ele)
    rhs_loc1 = matmul(local_solver_rhs,rhs_loc1)
    U_val = ele_val(newU,ele)
    do dim1 = 1, mesh_dim(U)
       rhs_loc2(u_start(dim1):u_end(dim1)) = u_val(dim1,:)
    end do
    rhs_loc2(d_start:d_end) = ele_val(newD,ele)
    rhs_loc2 = matmul(local_solver(1:d_end,1:d_end),rhs_loc2)
    rhs_loc2 = rhs_loc2-rhs_loc1
    do dim1 = 1, mesh_dim(U)
       call set(U_res,dim1,ele_nodes(U,ele),&
            rhs_loc2(u_start(dim1):u_end(dim1)))
    end do
    call set(D_res,ele_nodes(D,ele),&
         rhs_loc2(d_start:d_end))
  end subroutine get_linear_residuals_ele

  subroutine local_solve_residuals_ele(U_res,D_res,&
       local_solver_matrix,ele)
    type(vector_field), intent(inout) :: U_res
    type(scalar_field), intent(inout) :: D_res
    real, dimension(:,:), intent(in) :: local_solver_matrix
    integer, intent(in) :: ele
    !
    real, dimension(mesh_dim(U_res)*ele_loc(U_res,ele)+ele_loc(D_res,ele)+&
         &ele_n_constraints(U_res,ele)) :: rhs_loc
    real, dimension(mesh_dim(U_res),ele_loc(U_res,ele)) :: U_val
    integer :: uloc,dloc, d_start, d_end, dim1, mdim
    integer, dimension(mesh_dim(U_res)) :: U_start, U_end

    rhs_loc = 0.

    !Calculate indices in a vector containing all the U and D dofs in
    !element ele, First the u1 components, then the u2 components, then the
    !D components are stored.
    uloc = ele_loc(U_res,ele)
    dloc = ele_loc(D_res,ele)
    mdim = mesh_dim(U_res)
    do dim1 = 1, mdim
       u_start(dim1) = uloc*(dim1-1)+1
       u_end(dim1) = uloc*dim1
    end do
    d_start = uloc*mdim + 1
    d_end   = uloc*mdim+dloc

    U_val = ele_val(U_res,ele)
    do dim1 = 1, mesh_dim(U_res)
       rhs_loc(u_start(dim1):u_end(dim1)) = u_val(dim1,:)
    end do
    rhs_loc(d_start:d_end) = &
         & ele_val(D_res,ele)

    call solve(local_solver_matrix,Rhs_loc)
    do dim1 = 1, mesh_dim(U_res)
       call set(U_res,dim1,ele_nodes(U_res,ele),&
            rhs_loc(u_start(dim1):u_end(dim1)))
    end do
    call set(D_res,ele_nodes(D_res,ele),&
         rhs_loc(d_start:d_end))
  end subroutine local_solve_residuals_ele

  subroutine assemble_newton_solver_ele(D,f,U,X,down,lambda_rhs,ele, &
       &g,dt,theta,D0,&
       &lambda_mat,continuity_mat,&
       &Local_solver_matrix,Local_solver_rhs)
    implicit none
    type(scalar_field), intent(in) :: D,f
    type(scalar_field), intent(in) :: lambda_rhs
    type(vector_field), intent(in) :: U,X,down
    integer, intent(in) :: ele
    real, intent(in) :: g,dt,theta,D0
    type(csr_matrix), intent(inout) :: lambda_mat
    type(block_csr_matrix), intent(inout) :: continuity_mat
    real, dimension(:,:), intent(inout) ::&
         & local_solver_matrix
    real, dimension(:,:), intent(inout) ::&
         & local_solver_rhs
    !
    real, allocatable, dimension(:,:),target :: &
         &l_continuity_mat, l_continuity_mat2
    real, allocatable, dimension(:,:) :: helmholtz_loc_mat
    real, allocatable, dimension(:,:,:) :: continuity_face_mat
    real, allocatable, dimension(:,:) :: scalar_continuity_face_mat
    integer :: ni, face
    integer, dimension(:), pointer :: neigh
    type(element_type) :: U_shape
    integer :: stat, d_start, d_end, dim1, mdim, uloc,dloc, lloc, iloc
    integer, dimension(mesh_dim(U)) :: U_start, U_end
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

    !Calculate indices in a vector containing all the U and D dofs in
    !element ele, First the u1 components, then the u2 components, then the
    !D components are stored.
    do dim1 = 1, mdim
       u_start(dim1) = uloc*(dim1-1)+1
       u_end(dim1) = uloc*dim1
    end do
    d_start = uloc*mdim + 1
    d_end   = uloc*mdim+dloc

    !Get pointers to different parts of l_continuity_mat
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
         & g,dt,theta,D0,pullback=.false.,&
         & have_constraint=have_constraint,&
         & projection=.false.,poisson=.false.,&
         & local_solver_rhs=local_solver_rhs)

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
       do  dim1 = 1, mdim
          call addto(continuity_mat,dim1,1,face_global_nodes(U,face)&
               &,face_global_nodes(lambda_rhs,face),&
               &continuity_face_mat(dim1,:,:))
       end do

       deallocate(continuity_face_mat)
    end do

    !compute l_continuity_mat2 = inverse(local_solver)*l_continuity_mat
    allocate(l_continuity_mat2(uloc*2+dloc+n_constraints,lloc))
    l_continuity_mat2 = l_continuity_mat
    call solve(local_solver_matrix,l_continuity_mat2)

    !compute helmholtz_loc_mat
    allocate(helmholtz_loc_mat(lloc,lloc))
    helmholtz_loc_mat = matmul(transpose(l_continuity_mat),l_continuity_mat2)

    !insert helmholtz_loc_mat into global lambda matrix
    call addto(lambda_mat,ele_nodes(lambda_rhs,ele),&
         ele_nodes(lambda_rhs,ele),helmholtz_loc_mat)
  end subroutine assemble_newton_solver_ele
 
  subroutine assemble_hybridized_helmholtz_ele(D,f,U,X,down,ele, &
       g,dt,theta,D0,pullback,lambda_mat,lambda_rhs,U_rhs,D_rhs,&
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
    logical, intent(in) :: pullback
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
    integer :: stat, d_start, d_end, dim1, mdim, uloc,dloc, lloc, iloc
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
         & g,dt,theta,D0,pullback,have_constraint,&
         & projection=projection,poisson=poisson,&
         & local_solver_rhs=local_solver_rhs)
    if(poisson.and.(ele==1)) then
          !Fix the zero level of D
       local_solver_matrix(d_start,:) = 0.
       local_solver_matrix(:,d_start) = 0.
       local_solver_matrix(d_start,d_start) = 1.
    end if

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
    call assemble_rhs_ele(Rhs_loc,D,U,X,ele,pullback,&
         D_rhs,U_rhs,u_rhs_local)
    if(poisson.and.(ele==1)) then
       !Fix the zero level of D
       rhs_loc(d_start)=0.0
    end if
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
       g,dt,theta,D0,pullback,U_rhs,D_rhs,lambda,&
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
    logical, intent(in) :: pullback
    real, intent(in) :: g,dt,theta,D0
    logical, intent(in), optional :: projection, poisson,u_rhs_local
    !
    real, allocatable, dimension(:,:,:) :: continuity_face_mat
    integer :: ni, face
    integer, dimension(:), pointer :: neigh, d_nodes
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
    real :: constraint_check, u_max

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
         & g,dt,theta,D0,pullback,have_constraint,&
         & projection=projection,poisson=poisson,&
         & local_solver_rhs=local_solver_rhs)
    if(poisson.and.(ele==1)) then
          !Fix the zero level of D
       local_solver_matrix(d_start,:) = 0.
       local_solver_matrix(:,d_start) = 0.
       local_solver_matrix(d_start,d_start) = 1.
    end if

    !Construct the rhs sources for U from lambda
    call assemble_rhs_ele(Rhs_loc,D,U,X,ele,pullback,&
         D_rhs,U_rhs,U_rhs_local)
    if(poisson.and.(ele==1)) then
       !Fix the zero level of D
       rhs_loc(d_start)=0.0
    end if
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
          u_max = 0.
          do dim1 = 1, mdim
             u_max = max(u_max,maxval(abs(U_solved(dim1,:))))
             constraint_check = constraint_check + &
                  & sum(U_solved(dim1,:)*constraints%orthogonal(i1,:,dim1))
          end do
          if(abs(constraint_check)/(max(1.0,u_max))>1.0e-8) then
             ewrite(2,*) 'Constraint check', constraint_check
             FLAbort('Constraint not enforced')
          end if
       end do
    end if

  end subroutine reconstruct_U_d_ele

  subroutine get_local_solver(local_solver_matrix,U,X,down,D,f,ele,&
       & g,dt,theta,D0,pullback,have_constraint, &
       & projection,poisson,local_solver_rhs)
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
    logical, intent(in) :: pullback
    real, dimension(:,:)&
         &, intent(inout) :: local_solver_matrix
    real, dimension(:,:)&
         &, intent(inout), optional :: local_solver_rhs
    logical, intent(in) :: projection, poisson
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
         &,ele)) :: l_u_mat, M_lin
    real, dimension(mesh_dim(U),ele_loc(U,ele),ele_loc(U,ele)) :: A
    integer :: mdim, uloc,dloc,dim1,dim2,gi
    type(element_type) :: u_shape, d_shape
    real, dimension(ele_ngi(D,ele)) :: detwei, detJ
    integer, dimension(:), pointer :: D_ele,U_ele
    integer :: d_start, d_end
    integer, dimension(mesh_dim(U)) :: U_start, U_end
    type(constraints_type), pointer :: constraints
    integer :: i1
    real :: l_dt,l_theta,l_d0,detJ_bar,alpha

    if(projection) then
       l_dt = 0.
       l_theta = 0.
       l_d0 = 1.
    else if(poisson) then
       l_dt = 1.
       l_theta = 1.
       l_d0 = 1.
    else
       l_dt = dt
       l_theta = theta
       l_d0 = d0
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

    if(projection.or.poisson) then
       f_gi = 0.
    else
       f_gi = ele_val_at_quad(f,ele)
    end if
    up_gi = -ele_val_at_quad(down,ele)

    call get_up_gi(X,ele,up_gi)

    !J, detJ is needed for Piola transform
    !detwei is needed for pressure mass matrix
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
         detwei=detwei, detJ=detJ)
    detJ_bar = sum(detJ)/size(detJ)

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
    !not included in pressure solver
    if(.not.poisson) then
       if(pullback) then
          local_solver_matrix(d_start:d_end,d_start:d_end)=&
               &shape_shape(d_shape,d_shape,detJ_bar**2/detJ*&
               D_shape%quadrature%weight)
       else
          local_solver_matrix(d_start:d_end,d_start:d_end)=&
               &shape_shape(d_shape,d_shape,detwei)
       end if
       if(present(local_solver_rhs)) then
          if(pullback) then
             local_solver_rhs(d_start:d_end,d_start:d_end) = &
                  shape_shape(d_shape,d_shape,detJ_bar**2/detJ*&
                  D_shape%quadrature%weight)
          else
             local_solver_rhs(d_start:d_end,d_start:d_end) = &
                  shape_shape(d_shape,d_shape,detwei)
          end if
       end if
    end if
    !divergence matrix (done in local coordinates)
    if(pullback) then
       l_div_mat = dshape_shape(u_shape%dn,d_shape,&
            &detJ_bar/detJ*D_shape%quadrature%weight)
    else
       l_div_mat = dshape_shape(u_shape%dn,d_shape,D_shape%quadrature%weight)
    end if
    do dim1 = 1, mdim
       !pressure gradient term [integrated by parts so minus sign]
       local_solver_matrix(u_start(dim1):u_end(dim1),d_start:d_end)=&
            & -g*l_dt*l_theta*l_div_mat(dim1,:,:)
       if(present(local_solver_rhs)) then
          local_solver_rhs(u_start(dim1):u_end(dim1),d_start:d_end)=&
               & -g*(l_theta-1.0)*l_dt*l_div_mat(dim1,:,:)
       end if
       !divergence continuity term
       local_solver_matrix(d_start:d_end,u_start(dim1):u_end(dim1))=&
            & l_d0*l_dt*l_theta*transpose(l_div_mat(dim1,:,:))
       if(present(local_solver_rhs)) then
          local_solver_rhs(d_start:d_end,u_start(dim1):u_end(dim1))=&
               & l_d0*(l_theta-1.0)*l_dt*transpose(l_div_mat(dim1,:,:))
       end if
    end do
    !velocity mass matrix and Coriolis matrix (done in local coordinates)
    l_u_mat = shape_shape_tensor(u_shape, u_shape, &
         u_shape%quadrature%weight, Metric+l_dt*l_theta*Metricf)

    alpha=0.05
    A=0.0
    A(1,1,:)=(/ 1, 0,-1, 0, 0, 0, 0, 0, 0/)
    A(1,2,:)=(/ 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    A(1,3,:)=(/-1, 0, 1, 0, 0, 0, 0, 0, 0/)
    A(1,4,:)=(/ 0, 0, 0, 1, 0,-1, 0, 0, 0/)
    A(1,5,:)=(/ 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    A(1,6,:)=(/ 0, 0, 0,-1, 0, 1, 0, 0, 0/)
    A(1,7,:)=(/ 0, 0, 0, 0, 0, 0, 1, 0,-1/)
    A(1,8,:)=(/ 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    A(1,9,:)=(/ 0, 0, 0, 0, 0, 0,-1, 0, 1/)

    A(2,1,:)=(/ 1, 0, 0, 0, 0, 0,-1, 0, 0/)
    A(2,2,:)=(/ 0, 1, 0, 0, 0, 0, 0,-1, 0/)
    A(2,3,:)=(/ 0, 0, 1, 0, 0, 0, 0, 0,-1/)
    A(2,4,:)=(/ 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    A(2,5,:)=(/ 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    A(2,6,:)=(/ 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    A(2,7,:)=(/-1, 0, 0, 0, 0, 0, 1, 0, 0/)
    A(2,8,:)=(/ 0,-1, 0, 0, 0, 0, 0, 1, 0/)
    A(2,9,:)=(/ 0, 0,-1, 0, 0, 0, 0, 0, 1/)

    do dim1=1, mdim
       do dim2=1, mdim
          M_lin(dim1,dim2,:,:)=matmul(matmul(A(dim1,:,:),l_u_mat(dim1,dim2,:,:)),transpose(A(dim2,:,:)))
       end do
    end do
    l_u_mat=l_u_mat+alpha*M_lin


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
            Metric+(l_theta-1.0)*l_dt*Metricf)
       l_u_mat=l_u_mat+alpha*M_lin
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
    call compute_jacobian(ele_val(X,ele),ele_shape(X,ele), J=J, & 
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

  subroutine check_continuity_local_ele(U,Lambda,psi,X,ele)
    implicit none
    type(vector_field), intent(in) :: U,X
    type(scalar_field), intent(in) :: Lambda,psi
    integer, intent(in) :: ele
    !
    integer, dimension(:), pointer :: neigh
    integer :: ni,face,ele2,face2, i,dim1, loc
    type(constraints_type), pointer :: constraints
    real :: residual
    real, dimension(U%dim, ele_loc(U,ele)) :: U_loc

    !!Check internal continuity constraints    
    constraints => U%mesh%shape%constraints
    if(constraints%n_constraints>0) then
       !constraint%orthogonal(i,loc,dim1) stores the coefficient 
       !for basis function loc, dimension dim1 in equation i.
       U_loc = ele_val(U,ele)
       do i = 1, size(constraints%orthogonal,1)
          residual = 0.
          do dim1 = 1, mesh_dim(U)
             do loc = 1, ele_loc(U,ele)
                residual = residual + &
                     sum(U_loc(dim1,:)*constraints%orthogonal(i,:,dim1))
             end do
          end do
          ewrite(2,*) 'cjc residual', residual,maxval(abs(U_loc))
          assert(abs(residual/max(1.0,maxval(abs(U_loc))))<1.0e-10)
       end do
    end if

    neigh => ele_neigh(U,ele)
    do ni = 1, size(neigh)
       ele2 = neigh(ni)
       face = ele_face(U,ele,ele2)
       if(ele2>0) then
          face2 = ele_face(U,ele2,ele)
       else
          face2 = -1
       end if
       call check_continuity_local_face(U,Lambda,psi,X,ele,ele2,face,face2)
    end do
  end subroutine check_continuity_local_ele

  subroutine check_continuity_local_face(U,Lambda,psi,X,ele,ele2,face,face2)
    implicit none
    integer, intent(in) :: face,face2,ele,ele2
    type(vector_field), intent(in) :: U,X
    type(scalar_field), intent(in) :: Lambda, psi
    !
    real, dimension(U%dim, face_ngi(U, face)) :: n1,n2,u1,u2
    real, dimension(face_ngi(U, face)) :: f1,f2
    real :: weight1,weight2
    real, dimension(face_loc(Lambda,ele)) :: jump
    type(element_type), pointer :: U_face_shape
    real, pointer, dimension(:) :: detwei
    
    if(.not.all(face_global_nodes(lambda,face).eq.face_global_nodes(lambda&
         &,face2))) then
       ewrite(0,*) ele,ele2,face,face2
       ewrite(0,*) ele_nodes(lambda,ele)
       ewrite(0,*) 'ASDF'
       ewrite(0,*) ele_nodes(lambda,ele2)
       ewrite(0,*) face_global_nodes(lambda,face)
       ewrite(0,*) face_global_nodes(lambda,face2)
       ewrite(0,*) face_local_nodes(lambda,face)
       ewrite(0,*) face_local_nodes(lambda,face2)
    end if

    U_face_shape=>face_shape(U, face)
    detwei => U_face_shape%quadrature%weight

    !Get normal in local coordinates
    call get_local_normal(n1,weight1,U,local_face_number(U%mesh,face))
    call get_local_normal(n2,weight2,U,local_face_number(U%mesh,face2))
    u1 = face_val_at_quad(U,face)
    u2 = face_val_at_quad(U,face2)
    !jump = maxval(abs(sum(u1*n1+u2*n2,1)))

    jump = shape_rhs(face_shape(Lambda,ele),sum(u1*n1*weight1 &
         +u2*n2*weight2,1)*detwei)
    !ewrite(0,*) 'cjc jump', jump
    !ewrite(0,*) face_local_nodes(U,face)
    !ewrite(0,*) face_local_nodes(U,face2)
    if(maxval(abs(jump))/max(1.0,maxval(abs(u1)))>1.0e-7) then
       jump = shape_rhs(face_shape(Lambda,ele),sum(u2*n2*weight2,1)*detwei)
       ewrite(0,*) 'one side', jump
       jump = shape_rhs(face_shape(Lambda,ele),sum(u1*n1*weight1,1)*detwei)
       ewrite(0,*) 'other side', jump
       ewrite(0,*) 'Bad jump alert'
       ewrite(0,*) 'face numbers', face, face2
       ewrite(0,*) 'element numbers', ele, ele2
       ewrite(0,*) 'face X values', face_val(X,face)
       ewrite(0,*) shape_rhs(face_shape(Lambda,ele),sum( &
         u1*n1*weight1,1)*detwei)
       ewrite(0,*) shape_rhs(face_shape(Lambda,ele),sum( &
         u2*n2*weight2,1)*detwei)
       ewrite(0,*) 'psi1', face_val_at_quad(psi,face)
       ewrite(0,*) 'psi2', face_val_at_quad(psi,face2)
       FLExit('Bad jumps')
    end if

  end subroutine check_continuity_local_face

  subroutine check_continuity_ele(U_cart,X,ele,tolerance)
    implicit none
    type(vector_field), intent(in) :: U_cart,X
    integer, intent(in) :: ele
    real, intent(in) :: tolerance
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
       call check_continuity_face(U_cart,X,ele,ele2,face,face2,tolerance)
    end do
  end subroutine check_continuity_ele

  subroutine check_continuity_face(U_cart,X,ele,ele2,face,face2,tolerance)
    !subroutine to check the continuity of normal component
    !of velocity at quadrature points
    implicit none
    type(vector_field), intent(in) :: U_cart,X
    integer, intent(in) :: face,face2,ele,ele2
    real, intent(in) :: tolerance
    !
    real, dimension(X%dim, face_ngi(U_cart, face)) :: n1,n2
    real, dimension(X%dim, face_ngi(U_cart, face)) :: u1,u2,x1,x2
    real, dimension(face_ngi(U_cart, face)) :: jump_at_quad
    integer :: dim1
    !
    u1 = face_val_at_quad(U_cart,face)
    x1 = face_val_at_quad(X,face)
    if(ele2>0) then
       u2 = face_val_at_quad(U_cart,face2)
       x2 = face_val_at_quad(X,face)
       if(any(x1.ne.x2)) then
          ewrite(0,*) 'Face 1 X', x1
          ewrite(0,*) 'Face 2 X', x2
          FLExit('Something wrong with mesh?')
       end if

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
    if(maxval(abs(jump_at_quad))>tolerance) then
       ewrite(2,*) 'Jump at quadrature face, face2 =', jump_at_quad
       ewrite(2,*) 'ELE = ',ele,ele2
       do dim1 = 1, X%dim
          ewrite(2,*) 'normal',dim1,n1(dim1,:)
          ewrite(2,*) 'normal',dim1,n2(dim1,:)
          ewrite(2,*) 'X',dim1,x1(dim1,:)
       end do
       ewrite(2,*) 'n cpt1',sum(n1*u1,1)
       ewrite(2,*) 'n cpt2',sum(n2*u2,1)
       ewrite(2,*) jump_at_quad/max(maxval(abs(u1)),maxval(abs(u2)))
       FLAbort('stopping because of jumps')
    end if
    
  end subroutine check_continuity_face

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

  subroutine assemble_rhs_ele(Rhs_loc,D,U,X,ele,pullback,&
       D_rhs,U_rhs,u_rhs_local)
    implicit none
    integer, intent(in) :: ele
    logical, intent(in) :: pullback
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
    real :: detJ_bar

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
       !This will be multiplied by the local_solver_rhs matrix later.
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
       
       call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
            detJ=detJ,detwei=detwei)
       detJ_bar = sum(detJ)/size(detJ)
       
       Rhs_loc = 0.
       if(have_d_rhs) then
          if(pullback) then
             Rhs_loc(d_start:d_end) = shape_rhs(ele_shape(D,ele),&
                  &ele_val_at_quad(l_D_rhs,ele)*&
                  &detJ_bar**2/detJ*u_shape%quadrature%weight)
          else
             Rhs_loc(d_start:d_end) = shape_rhs(ele_shape(D,ele),&
                  &ele_val_at_quad(l_D_rhs,ele)*detwei)
          end if
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
                !Don't divide by detJ as we can use quadrature weight
                !instead of detwei
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

    ewrite(2,*) 'Energy:= ', energy
    ewrite(2,*) 'Percentage Change in energy:= ', (energy-old_energy)/energy

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
    call compute_jacobian(x_val,x_shape, J=J, &
         detwei=detwei)

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
    type(scalar_field) :: D_rhs,tmp_field, lambda
    type(vector_field), pointer :: U_local,down,X, U_cart
    type(vector_field) :: Coriolis_term, Balance_eqn, tmpV_field
    type(mesh_type),pointer :: Lambda_Mesh
    integer :: ele,dim1,i1, stat
    real :: g, D0
    logical :: elliptic_method, pullback
    real :: u_max, b_val, h_mean, area

    D=>extract_scalar_field(state, "LayerThickness")
    psi=>extract_scalar_field(state, "Streamfunction")
    f=>extract_scalar_field(state, "Coriolis")
    U_local=>extract_vector_field(state, "LocalVelocity")
    U_cart=>extract_vector_field(state, "Velocity")
    X=>extract_vector_field(state, "Coordinate")
    down=>extract_vector_field(state, "GravityDirection")
    call get_option("/physical_parameters/gravity/magnitude", g)
    call allocate(tmpV_field,mesh_dim(U_local), U_local%mesh, "tmpV_field")
    call allocate(D_rhs,D%mesh,'BalancedSolverRHS')
    call allocate(tmp_field,D%mesh,'TmpField')
    call allocate(Coriolis_term,mesh_dim(U_local),&
         U_local%mesh,"CoriolisTerm")
    call allocate(balance_eqn,mesh_dim(D),u_local%mesh,'BalancedEquation')

    !STAGE 1: Set velocity from streamfunction
    do ele = 1, element_count(D)
       call set_local_velocity_from_streamfunction_ele(&
            &U_local,psi,down,X,ele)
    end do

    !STAGE 1a: verify that velocity projects is div-conforming
    !call project_local_to_cartesian(X,U_local,U_cart)
    !do ele = 1, ele_count(U_local)
    !   call check_continuity_ele(U_cart,X,ele,tolerance=1.0e-8)
    !end do
    !Stage 1b: verify that projection is idempotent
    ewrite(2,*) 'CHECKING CONTINUOUS', maxval(abs(u_local%val))
    tmpV_field%val = U_local%val
    call project_to_constrained_space(state,tmpV_field)

    ewrite(2,*) maxval(abs(U_local%val-tmpV_field%val)), 'continuity'
    u_max = max(u_max,maxval(abs(U_local%val)))
    ewrite(2,*) u_max, 'u_max'

    call vtk_write_fields('ContinuityCheck', position=X, &
         model=D%mesh, &
         vfields=(/tmpV_field/))    

    lambda_mesh=>extract_mesh(state, "VelocityMeshTrace")
    
    ewrite(0,*) lambda_mesh%shape%facet2dofs(1)%dofs
    ewrite(0,*) lambda_mesh%shape%facet2dofs(2)%dofs
    ewrite(0,*) lambda_mesh%shape%facet2dofs(3)%dofs
    ewrite(0,*) lambda_mesh%shape%facet2dofs(4)%dofs
    !stop

    call allocate(lambda,lambda_mesh,"lambdacheck")
    do ele = 1, element_count(D)
       call check_continuity_local_ele(U_local,Lambda,psi,X,ele)
    end do
    do ele = 1, element_count(D)
       call check_continuity_local_ele(tmpV_field,Lambda,psi,X,ele)
    end do
    call deallocate(lambda)
    tmpV_field%val = tmpV_Field%val - u_local%val

    assert(maxval(abs(tmpV_field%val)/max(1.0,u_max))<1.0e-8)

    if(have_option("/material_phase::Fluid/vector_field::Velocity/prognostic&
         &/initial_condition::WholeMesh/balanced/elliptic_solver")) then
       ewrite(2,*) 'Using elliptic solver'
       !Construct Coriolis term
       call zero(Coriolis_term)
       do ele = 1, element_count(D)
          call set_coriolis_term_ele(Coriolis_term,f,down,U_local,X,ele)
       end do

       !Project Coriolis term into div-conforming space
       call project_to_constrained_space(state,Coriolis_term)

       !Calculate divergence of Coriolis term 
       call zero(d_rhs)

       pullback = .false.
       
       do ele = 1, element_count(D)
          call compute_divergence_ele(Coriolis_term,d_rhs,X,ele,pullback)
       end do

       !Solve for balanced layer depth

       ewrite(2,*) 'Solving elliptic problem for balanced layer depth.'
       call solve_hybridized_helmholtz(state,d_Rhs=d_rhs,&
            &U_out=tmpV_field,d_out=tmp_field,&
            &compute_cartesian=.true.,&
            &poisson=.true.,u_rhs_local=.true.)
       D%val = tmp_field%val
    else
       ewrite(2,*) 'Using streamfunction projection'
       !Project the streamfunction into pressure space
       do ele = 1, element_count(D)
          call project_streamfunction_for_balance_ele(D,psi,X,f,g,ele)
       end do
    end if

    !Subtract off the mean part
    h_mean = 0.
    area = 0.
    do ele = 1, element_count(D)
       call assemble_mean_ele(D,X,h_mean,area,ele)
    end do
    h_mean = h_mean/area
    D%val = D%val - h_mean
    !Add back on the correct mean depth
    call get_option("/material_phase::Fluid/scalar_field::LayerThickness/&
         &prognostic/mean_layer_thickness",D0)
    D%val = D%val + D0

    !debugging tests
    if(.true.) then
       call zero(Coriolis_term)
       do ele = 1, element_count(D)
          call set_coriolis_term_ele(Coriolis_term,f,down,U_local,X,ele)
       end do
       
       call zero(balance_eqn)
       do ele = 1, element_count(D)
          call set_pressure_force_ele(balance_eqn,D,X,g,ele,&
               pullback)
       end do
       call addto(balance_eqn,coriolis_term,scale=1.0)
    else
       ! call zero(balance_eqn)
       ! do ele = 1, element_count(D)
       !    call set_balance_eqn_ele(balance_eqn,f,U_local,D,X,g,ele)
       ! end do
    end if

    ewrite(2,*) 'Project balance equation into div-conforming space'
    ewrite(2,*) 'CJC b4',maxval(abs(balance_eqn%val)),&
         & maxval(abs(coriolis_term%val))
    b_val = maxval(abs(balance_eqn%val))
    !Project balance equation into div-conforming space
    call solve_hybridized_helmholtz(state,U_Rhs=balance_eqn,&
         &U_out=balance_eqn,&
         &compute_cartesian=.false.,&
         &projection=.true.,&
         &poisson=.false.,&
         &u_rhs_local=.true.)
    
    call vtk_write_fields('BalanceEqn', position=X, &
         model=D%mesh, &
         vfields=(/X,balance_eqn/))    

    do dim1 = 1, mesh_dim(D)
       ewrite(2,*) 'Balance equation', maxval(abs(balance_eqn%val(dim1,:)))/b_val
       !assert(maxval(abs(balance_eqn%val(dim1,:)))/b_val<1.0e-8)
    end do
    
 !Clean up after yourself
    call deallocate(Coriolis_term)
    call deallocate(D_rhs)
    call deallocate(balance_eqn)
    call deallocate(tmpV_field)
    call deallocate(tmp_field)

  end subroutine set_velocity_from_geostrophic_balance_hybridized
  
  ! subroutine set_balance_eqn_ele(balance_eqn,f,U_local,D,X,g,ele)
  !   type(vector_field), intent(inout) :: balance_eqn
  !   type(vector_field), intent(in) :: U_local,X
  !   type(scalar_field), intent(in) :: D,f
  !   real, intent(in) :: g
  !   integer, intent(in) :: ele
  !   !
  !   real, dimension(mesh_dim(balance_eqn), X%dim, ele_ngi(balance_eqn,ele)) :: J
  !   real, dimension(ele_ngi(balance_eqn,ele)) :: detJ
  !   real, dimension(mesh_dim(balance_eqn),mesh_dim(balance_eqn),ele_ngi(balance_eqn,ele))::&
  !        &Metric
  !   real, dimension(mesh_dim(balance_eqn)*ele_loc(balance_eqn,ele),&
  !        mesh_dim(balance_eqn)*ele_loc(balance_eqn,ele)) :: l_u_mat
  !   real, dimension(mesh_dim(balance_eqn)*ele_loc(balance_eqn,ele)) :: balance_eqn_rhs
  !   real, dimension(mesh_dim(D),ele_loc(balance_eqn,ele)) :: &
  !        & rhs_loc
  !   type(element_type) :: force_shape
  !   real, dimension(ele_ngi(D,ele)) :: D_gi
  !   integer :: uloc
  !   real, dimension(mesh_dim(U_local), ele_ngi(U_local,ele)) :: U_gi
    
  !   D_gi = ele_val_at_quad(D,ele)
  !   u_gi = ele_val_at_quad(u_local,ele)

  !   uloc = ele_loc(balance_eqn,ele)
  !   force_shape = ele_shape(balance_eqn,ele)
  !   call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
  !        detJ=detJ)
  !   do gi=1,ele_ngi(balance_eqn,ele)
  !      Metric(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)
  !   end do

  !   rhs_loc = -g*dshape_rhs(force_shape%dn,&
  !        D_gi*D%mesh%shape%quadrature%weight)
  !   rhs_loc(1,:) = rhs_loc(1,:) + -f*shape_rhs(force_shape,&
  !        u_gi(2,
  !   do dim1 = 1, mesh_dim(balance_eqn)
  !      force_rhs((dim1-1)*uloc+1:dim1*uloc) = rhs_loc(dim1,:)
  !   end do
  !   do dim1 = 1, mesh_dim(balance_eqn)
  !      do dim2 = 1, mesh_dim(balance_eqn)
  !         l_u_mat((dim1-1)*uloc+1:dim1*uloc,&
  !              &  (dim2-1)*uloc+1:dim2*uloc ) = &
  !              & shape_shape(force_shape,force_shape,&
  !              & force_shape%quadrature%weight*Metric(dim1,dim2,:))
  !      end do
  !   end do

  ! end subroutine set_balance_eqn_ele

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

    !call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
    !     detwei=detwei)

    detwei = D%mesh%shape%quadrature%weight
    d_rhs = shape_rhs(d_shape,detwei*psi_quad*f_gi/g)
    d_mass = shape_shape(d_shape,d_shape,detwei)
    call solve(d_mass,d_rhs)
    call set(D,ele_nodes(D,ele),d_rhs)

  end subroutine project_streamfunction_for_balance_ele
  
  subroutine set_pressure_force_ele(force,D,X,g,ele,pullback)
    implicit none
    type(vector_field), intent(inout) :: force
    type(scalar_field), intent(in) :: D
    type(vector_field), intent(in) :: X
    real, intent(in) :: g
    integer, intent(in) :: ele
    logical, intent(in) :: pullback
    !
    real, dimension(ele_ngi(D,ele)) :: D_gi
    real, dimension(mesh_dim(D),ele_loc(force,ele)) :: &
         & rhs_loc
    integer :: dim1,dim2,gi,uloc
    real, dimension(mesh_dim(force), X%dim, ele_ngi(force,ele)) :: J
    real, dimension(ele_ngi(force,ele)) :: detJ
    real, dimension(mesh_dim(force),mesh_dim(force),ele_ngi(force,ele))::&
         &Metric
    real, dimension(mesh_dim(force)*ele_loc(force,ele),&
         mesh_dim(force)*ele_loc(force,ele)) :: l_u_mat
    real, dimension(mesh_dim(force)*ele_loc(force,ele)) :: force_rhs
    type(element_type) :: force_shape
    real :: detJ_bar

    uloc = ele_loc(force,ele)
    force_shape = ele_shape(force,ele)
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
         detJ=detJ)
    do gi=1,ele_ngi(force,ele)
       Metric(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)
    end do

    detJ_bar = sum(detJ)/size(detJ)

    D_gi = ele_val_at_quad(D,ele)
    if(pullback) then
       D_gi = D_gi*detJ_bar/detJ
    end if
    rhs_loc = -g*dshape_rhs(force%mesh%shape%dn,&
         D_gi*D%mesh%shape%quadrature%weight)
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
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
         detJ=detJ)
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
    real, dimension(mesh_dim(psi),ele_ngi(psi,ele)) :: dpsi_gi, Uloc_gi
    real, dimension(ele_ngi(psi,ele)) :: div_gi
    real, dimension(mesh_dim(U_local),ele_loc(U_local,ele)) :: U_loc
    real, dimension(ele_loc(U_local,ele),ele_loc(U_local,ele)) :: &
         & l_mass_mat
    type(element_type) :: u_shape, psi_shape
    integer :: dim1,gi,uloc
    real, dimension(X%dim, ele_ngi(X,ele)) :: up_gi
    integer :: orientation
    real, dimension(mesh_dim(U_local), X%dim, ele_ngi(U_local,ele)) :: J
    real, dimension(ele_ngi(X,ele)) :: detJ

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
    end do

    !DEbugging check, did we get the same fields?
    do dim1= 1, U_local%dim
       Uloc_gi(dim1,:) = matmul(transpose(U_shape%n),U_loc(dim1,:))
    end do
    assert(maxval(abs(Uloc_gi-dpsi_gi))/max(1.0,maxval(abs(Uloc_gi)))<1.0e-11)

    !verify divergence-free-ness
    !This is just for debugging
    !Annoyingly it requires detJ when the rest of the subroutine doesn't
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
         detJ=detJ)
    div_gi = 0.
    do gi = 1, ele_ngi(psi,ele)
       do dim1 = 1, mesh_dim(psi)
          div_gi(gi) = div_gi(gi)+sum(u_shape%dn(:,gi,dim1)*u_loc(dim1,:))
       end do
       div_gi(gi)=div_gi(gi)/detJ(gi)
    end do
    if(maxval(abs(div_gi))>1.0e-8) then
       ewrite(0,*) 'Divergence =', maxval(abs(div_gi))
       FLExit('Divergence not small enough')
    end if

    do dim1 = 1, U_local%dim
       call set(U_local,dim1,ele_nodes(U_local,ele),&
            u_loc(dim1,:))
    end do

  end subroutine set_local_velocity_from_streamfunction_ele

  subroutine project_to_constrained_space(state,v_field)
    !wrapper for projecting vector field to constrained space
    implicit none
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
    !
    type(vector_field) :: v_field_out
    type(scalar_field) :: tmp_field
    type(scalar_field), pointer :: D

    call allocate(v_field_out,v_field%dim,v_field%mesh,"TmpVField")
    D=>extract_scalar_field(state, "LayerThickness")
    call allocate(tmp_field,D%mesh,"TmpSField")

    ewrite(1,*) 'Project to constrained space'

    call solve_hybridized_helmholtz(state,&
         &U_rhs=v_field,U_out=v_field_out,D_out=tmp_field,&
         &compute_cartesian=.true.,&
         &projection=.true.,&
         &poisson=.false.,&
         &u_rhs_local=.true.)

    call set(V_field,V_field_out)
    call deallocate(v_field_out)
    call deallocate(tmp_field)

    ewrite(1,*) 'END: Project to constrained space'

  end subroutine project_to_constrained_space

  subroutine solve_linear_timestep_hybridized(&
       &state,dt_in,theta_in)
    implicit none
    type(state_type), intent(inout) :: state
    real, intent(in) :: dt_in, theta_in
    !
    type(vector_field), pointer :: u,u_rhs
    type(scalar_field), pointer :: d,d_rhs
    
    D=>extract_scalar_field(state, "LayerThickness")
    d_rhs=>extract_scalar_field(state, "LayerThickness_RHS")
    U=>extract_vector_field(state, "LocalVelocity")
    u_rhs=>extract_vector_field(state, "LocalVelocity_RHS")
    
    call solve_hybridized_helmholtz(&
         &state,&
         &U_out=U_rhs,D_out=d_rhs,&
         &compute_cartesian=.true.,&
         &projection=.false.,poisson=.false.,&
         &output_dense=.false.)
    ewrite(1,*) 'jump in D', maxval(abs(d_rhs%val-d%val))
    ewrite(1,*) 'jump in U', maxval(abs(U_rhs%val-U%val))
    call set(d,d_rhs)
    call set(U,U_rhs)
  end subroutine solve_linear_timestep_hybridized

  subroutine assemble_mean_ele(D,X,mean,area,ele)
    type(scalar_field), intent(in) :: D
    type(vector_field), intent(in) :: X
    real, intent(inout) :: mean,area
    integer, intent(in) :: ele
    !
    real, dimension(ele_ngi(D,ele)) :: D_gi, detwei
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J

    D_gi = ele_val_at_quad(D,ele)
    call compute_jacobian(ele_val(X,ele),ele_shape(X,ele), J=J, &
         detwei=detwei)
    assert(all(detwei>0.))
    mean = mean + sum(detwei*D_gi)
    area = area + sum(detwei)
    
  end subroutine assemble_mean_ele

  subroutine compute_divergence_ele(V,Div_V,X,ele,pullback)
    !subroutine to compute the divergence of V in element ele
    !and insert the values into Div_V
    !V is represented by a div-conforming space
    !and Div_V is represented by the divergence of that space
    implicit none
    type(vector_field), intent(in) :: V, X
    type(scalar_field), intent(inout) :: Div_V
    integer, intent(in) :: ele
    logical, intent(in) :: pullback
    !
    real, dimension(ele_ngi(Div_V,ele)) :: detwei, detJ
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    integer :: dim1
    real, dimension(V%dim, ele_loc(V,ele)) :: U_loc
    type(element_type) :: u_shape, d_shape
    real, dimension(ele_loc(Div_V,ele)) :: Div_loc
    real, dimension(mesh_dim(V),ele_loc(V,ele),ele_loc(Div_V,ele)) :: l_div_mat
    real, dimension(ele_loc(Div_V,ele),ele_loc(Div_V,ele)) :: d_mass_mat
    real :: detJ_bar
    !
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
         detJ = detJ, detwei=detwei)
    
    detJ_bar = sum(detJ)/size(detJ)

    u_shape = ele_shape(V,ele)
    d_shape = ele_shape(Div_V,ele)
    U_loc = ele_val(V,ele)
    
    if(pullback) then
       l_div_mat = dshape_shape(u_shape%dn,d_shape,&
            &detJ_bar/detJ*D_shape%quadrature%weight)
    else
       l_div_mat = dshape_shape(u_shape%dn,d_shape,&
            &D_shape%quadrature%weight)
    end if

    div_loc = 0.
    do dim1 = 1, V%dim
       div_loc = div_loc + matmul(transpose(l_div_mat(dim1,:,:))&
            &,U_loc(dim1,:))
    end do
    if(pullback) then
       d_mass_mat = shape_shape(d_shape,d_shape,detJ_bar*D_shape%quadrature&
            &%weight)
    else
       d_mass_mat = shape_shape(d_shape,d_shape,detwei)
    end if
    call solve(d_mass_mat,div_loc)

    !dn is loc x ngi x dim
    !div_gi = 0.
    !do iloc = 1, ele_loc(Div_V,ele)
    !   div_gi = div_gi + matmul(u_shape%dn(iloc,:,:),U_loc(:,iloc))
    !end do
    !div_loc = shape_rhs(d_shape,div_gi) !no detwei as cancels
    !d_mass_mat = shape_shape(d_shape,d_shape,detwei)
    !call solve(d_mass_mat,div_loc)

    call set(div_V,ele_nodes(div_V,ele),div_loc)
  end subroutine compute_divergence_ele

  subroutine check_zero_level(D)
    type(scalar_field), intent(in) :: D
    !
    real, dimension(ele_loc(D,1)) :: D_val
    
    D_val = ele_val(D,1)
    assert(abs(D_val(1))<1.0e-10)
  end subroutine check_zero_level

  subroutine set_layerthickness_projection(state,name)
    type(state_type), intent(in) :: state
    character(len=*), intent(in), optional :: name
    !
    type(vector_field), pointer :: X
    type(scalar_field), pointer :: D
    character(len=PYTHON_FUNC_LEN) :: Python_Function
    integer :: ele
    real :: h_mean,area

    if(present(name)) then
       D=>extract_scalar_field(state, trim(name))
       call get_option(trim(D%option_path)//"/prescribed/value/python",&
            & Python_Function)
    else
       D=>extract_scalar_field(state, "LayerThickness")
       call get_option("/material_phase::Fluid/scalar_field&
            &::LayerThickness/prognostic/initial_condition&
            &::ProjectionFromPython/python",&
            Python_Function)
   end if
    X=>extract_vector_field(state, "Coordinate")
 
    do ele = 1, element_count(D)
       call set_layerthickness_projection_ele(D,X,Python_Function,ele)
    end do

  end subroutine set_layerthickness_projection

  subroutine set_layerthickness_projection_ele(D,X,Python_Function,ele)
    type(scalar_field), intent(inout) :: D
    type(vector_field), intent(inout) :: X
    character(len=PYTHON_FUNC_LEN), intent(in) :: Python_Function
    integer, intent(in) :: ele
    !
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(ele_ngi(D,ele)) :: detwei, D_rhs_gi
    real, dimension(X%dim,ele_ngi(D,ele)) :: X_gi
    real, dimension(ele_loc(D,ele)) :: D_rhs
    real, dimension(ele_loc(D,ele),ele_loc(D,ele)) :: mass_mat
    type(element_type), pointer :: D_shape
    integer :: stat

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
         detwei=detwei)
    D_shape => ele_shape(D,ele)
    mass_mat = shape_shape(D_shape,D_Shape,detwei)
    
    X_gi = ele_val_at_quad(X,ele)
    call set_scalar_field_from_python(python_function, len(python_function),&
         & dim=3,nodes=ele_ngi(X,ele),x=X_gi(1,:),y=X_gi(2,:)&
         &,z=x_gi(3,:),t=0.0,&
         & result=D_rhs_gi,&
         & stat=stat)
    if(stat /= 0) then
       FLAbort('Failed to set face values from Python.')
    end if
    
    D_rhs = shape_rhs(D_shape,detwei*D_rhs_gi)
    call solve(mass_mat,D_rhs)
    
    call set(D,ele_nodes(D,ele),D_rhs)

  end subroutine set_layerthickness_projection_ele

end module hybridized_helmholtz
