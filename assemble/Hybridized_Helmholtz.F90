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
    type(scalar_field), pointer :: D,f,U,X,down, lambda
    type(scalar_field) :: lambda_rhs
    type(csr_sparsity) :: lambda_sparsity
    integer :: mdim
    ewrite(1,*) '  subroutine assemble_hybridized_helmholtz('

    !Pull the fields out of state
    D=>extract_scalar_field(state, "LayerThickness")
    f=>extract_scalar_field(state, "Coriolis")
    U=>extract_vector_field(state, "LocalVelocity")
    X=>extract_vector_field(state, "Coordinate")
    down=>extract_vector_field(state, "GravityDirection")
    lambda=>extract_scalar_field(state, "LagrangeMultiplier")

    !construct/extract sparsities
    lambda_sparsity=get_csr_sparsity_firstorder(state, lambda%mesh, lambda%mesh)

    !allocate matrices
    call allocate(lambda_mat,lambda_sparsity)
    call zero(lambda_mat)

    !allocate hybridized RHS
    call allocate(lambda_rhs,lambda%mesh,"LambdaRHS")
    call zero(lambda_rhs)
    
    !get parameters
      call get_option("/physical_parameters/gravity/magnitude", g)
      !theta
      call get_option("/material_phase::Fluid/scalar_field::LayerThickness/prognostic/temporal_discretisation/theta",theta)
      !D0
      call get_option("/material_phase::Fluid/scalar_field::LayerThickness/p&
           &rognostic/mean_layer_thickness",D0)
      call get_option("/timestepping/timestep", dt)

    !Assemble matrices
    do ele = 1, ele_count(D)
       call assemble_hybridized_helmholtz_ele(D,f,U,X,down,ele, &
            lambda_mat,g,dt,theta,D0,lambda_rhs,rhs)
    end do

    !Solve the equation
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
    !implicit none
    type(scalar_field), intent(inout) :: D,f,U,X,down,lambda_rhs
    integer, intent(in) :: ele
    real, intent(in) :: g,dt,theta,D0
    type(csr_matrix), intent(inout) :: lambda_mat
    type(scalar_field), intent(inout), optional :: rhs
    !
    real, dimension(ele_loc(U)*2*ele_loc(D),ele_loc(U)*ele_loc(D))&
         &:: local_solver
    real, allocatable, dimension(:,:) :: constraint_mat, constraint_mat2
    real, allocatable, dimension(:,:) :: helmholtz_loc_mat
    integer :: ni, lambda_ele_loc, face, ele_2
    integer, dimension(:), pointer :: neigh
    
    !construct local solver
    ...

    !get list of neighbours
    neigh => ele_neigh(ele)
    !get size of constraint_mat
    lambda_ele_loc = 0
    do ni = 1, size(neigh)
       ele_2 = neigh(ni)
       face=ele_face(U, ele, ele_2)
       lambda_ele_loc = lambda_ele_loc + face_loc(U,face)
    end do
    !allocate constraint_mat
    allocate(constraint_mat(ele_loc(U)*2*ele_loc(D),lambda_ele_loc))

    !calculate constraint_mat
    ...

    !compute constraint_mat2 = inverse(local_solver)*constraint_mat
    allocate(constraint_mat2(ele_loc(U)*2*ele_loc(D),lambda_ele_loc))
    constraint_mat2 = constraint_mat
    call solve(local_solver,constraint_mat)

    !compute helmholtz_loc_mat
    allocate(helmholtz_loc_mat(lambda_ele_loc,lambda_ele_loc))
    helmholtz_loc_mat = matmul(transpose(constraint_mat),constraint_mat2)

    !construct lambda_rhs
    ...
    
    !insert helmholtz_loc_mat into global lambda matrix
    do ni = 1, size(neigh)
       ele_2 = neigh(ni)
       face=ele_face(U, ele, ele_2)
       do ni2 = 1, size(neigh)
          ele_2 = neigh(ni)
          face2=ele_face(U, ele, ele_2)
          call addto(lambda_mat,ele_nodes(lambda_rhs,ele),&
               ele_nodes(lambda_rhs,ele))
       end do
    end do
  end subroutine assemble_hybridized_helmholtz_ele

end module hybridized_helmholtz
