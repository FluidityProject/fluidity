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

module burgers_adjoint_callbacks
#ifdef HAVE_ADJOINT
#include "libadjoint/adj_fortran.h"
    use libadjoint
    use libadjoint_data_callbacks
    use iso_c_binding
    use state_module
    use fields
    use sparse_matrices_fields
    use spud
    use adjoint_global_variables, only: adj_var_lookup
    use mangle_options_tree, only: adjoint_field_path
    use mangle_dirichlet_rows_module
    implicit none

    private

    public :: register_burgers_operator_callbacks

    contains

    subroutine register_burgers_operator_callbacks(adjointer)
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int) :: ierr

      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ASSEMBLY_CB, "VelocityIdentity", c_funloc(velocity_identity_assembly_callback))
      call adj_chkierr(ierr)
      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ASSEMBLY_CB, "BurgersOperator", c_funloc(burgers_operator_assembly_callback))
      call adj_chkierr(ierr)

      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, "VelocityIdentity", c_funloc(velocity_identity_action_callback))
      call adj_chkierr(ierr)
      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, "TimesteppingOperator", c_funloc(timestepping_operator_action_callback))
      call adj_chkierr(ierr)
    end subroutine register_burgers_operator_callbacks

    subroutine velocity_identity_assembly_callback(nvar, variables, dependencies, hermitian, coefficient, context, output, rhs) bind(c)
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(c_ptr), intent(in), value :: context
      type(adj_matrix), intent(out) :: output
      type(adj_vector), intent(out) :: rhs

      type(state_type), pointer :: matrices
      type(mesh_type), pointer :: u_mesh
      type(scalar_field) :: empty_u_field

      if (coefficient /= 1.0) then
        FLAbort("velocity_identity_assembly_callback requires that the coefficient is 1.0")
      end if

      call c_f_pointer(context, matrices)
      u_mesh => extract_mesh(matrices, "VelocityMesh")

      call allocate(empty_u_field, u_mesh, trim("IdentityAdjointVelocityOutput"))
      call zero(empty_u_field)
      rhs = field_to_adj_vector(empty_u_field)
      call deallocate(empty_u_field)
      ! Note: no dirichlet boundary condition is necessary since we are assemblying the indentity matrix

      ! We're not actually going to do anything daft like allocate memory for
      ! an identity matrix. Instead, we'll fill it in with a special tag
      ! that we'll catch later on in the solution of the adjoint equations.
      output%ptr = c_null_ptr
      output%klass = ADJ_IDENTITY_MATRIX
      output%flags = 0

    end subroutine velocity_identity_assembly_callback

    subroutine burgers_operator_assembly_callback(nvar, variables, dependencies, hermitian, coefficient, context, output, rhs) bind(c)
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(c_ptr), intent(in), value :: context
      type(adj_matrix), intent(out) :: output
      type(adj_vector), intent(out) :: rhs

      type(state_type), pointer :: matrices
      type(mesh_type), pointer :: u_mesh
      type(scalar_field) :: empty_velocity_field, previous_u, iter_u
      type(scalar_field), dimension(2) :: deps
      type(csr_matrix), pointer :: mass_mat, diffusion_mat
      type(csr_matrix) :: burgers_mat, burgers_mat_T
      real :: dt, theta
      integer :: i

      if (coefficient /= 1.0) then
        FLAbort("The coefficient in burgers_operator_assembly_callback has to be 1.0")
      end if

      if (nvar==2) then
        if (variables(1)%timestep==variables(2)%timestep-1) then
          call field_from_adj_vector(dependencies(1), previous_u)
          call field_from_adj_vector(dependencies(2), iter_u)
        else if(variables(2)%timestep==variables(1)%timestep-1) then
          call field_from_adj_vector(dependencies(2), previous_u)
          call field_from_adj_vector(dependencies(1), iter_u)
        else
          FLAbort("Dependencies have no contiguous timesteps.")
        end if
      else if (nvar==1) then
        call field_from_adj_vector(dependencies(1), previous_u)
        call field_from_adj_vector(dependencies(1), iter_u)
      else
        FLAbort("Unknown number of dependencies in burgers_operator_assembly_callback")
      end if

      call get_option(trim(adjoint_field_path("/material_phase::Fluid/scalar_field::Velocity/prognostic/temporal_discretisation/theta")), theta)
      call get_option("/timestepping/timestep", dt)

      call c_f_pointer(context, matrices)
      u_mesh => extract_mesh(matrices, "VelocityMesh")
      diffusion_mat => extract_csr_matrix(matrices, "DiffusionMatrix")
      mass_mat => extract_csr_matrix(matrices, "MassMatrix")

      call allocate(empty_velocity_field, u_mesh, trim("AdjointPressureRhs"))
      call zero(empty_velocity_field)
      rhs = field_to_adj_vector(empty_velocity_field)
      call deallocate(empty_velocity_field)

      ! Assemble the output matrix
      call allocate(burgers_mat, mass_mat%sparsity, name="BurgersAdjointOutput")
      call zero(burgers_mat)

      if (.not. have_option(trim(adjoint_field_path('/material_phase::Fluid/scalar_field::Velocity/prognostic/temporal_discretisation/remove_time_term')))) then
        call addto(burgers_mat, mass_mat, 1.0/abs(dt))
      end if
      !call addto(burgers_mat, advection_mat, theta)
      call addto(burgers_mat, diffusion_mat, theta)

      ! Ensure that the input velocity has indeed dirichlet bc attached
      if (get_boundary_condition_count(previous_u)==0) then
        FLAbort("Velocity in callback must have dirichlet conditions attached")
      end if
      ! Mangle the dirichlet rows
      call mangle_dirichlet_rows(burgers_mat, previous_u, keep_diag=.true.)

      ! And transpose of desired ...
      if (hermitian == ADJ_TRUE) then
        burgers_mat_T = transpose(burgers_mat)
        output = matrix_to_adj_matrix(burgers_mat_T)
        call deallocate(burgers_mat_T)
      else
        output = matrix_to_adj_matrix(burgers_mat)
      end if
      call deallocate(burgers_mat)
      output%flags = 0
    end subroutine burgers_operator_assembly_callback

    subroutine velocity_identity_action_callback(nvar, variables, dependencies, hermitian, coefficient, input, context, output) bind(c)
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(adj_vector), intent(in), value :: input
      type(c_ptr), intent(in), value :: context
      type(adj_vector), intent(out) :: output

      type(state_type), pointer :: matrices
      type(mesh_type), pointer :: eta_mesh

      type(scalar_field) :: eta_input
      type(scalar_field) :: eta_output

      call c_f_pointer(context, matrices)
      eta_mesh => extract_mesh(matrices, "VelocityMesh")

      call field_from_adj_vector(input, eta_input)
      call allocate(eta_output, eta_mesh, "IdentityVelocityOutput")
      call set(eta_output, eta_input)
      ! Note: no dirichlet boundary condition is necessary since we are assemblying the indentity matrix
      call scale(eta_output, coefficient)
      output = field_to_adj_vector(eta_output)
      call deallocate(eta_output)
    end subroutine velocity_identity_action_callback

    subroutine timestepping_operator_action_callback(nvar, variables, dependencies, hermitian, coefficient, input, context, output) bind(c)
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(adj_vector), intent(in), value :: input
      type(c_ptr), intent(in), value :: context
      type(adj_vector), intent(out) :: output

      type(state_type), pointer :: matrices
      type(csr_matrix), pointer :: diffusion_mat, mass_mat
      type(csr_matrix) :: timestepping_mat
      type(mesh_type), pointer :: u_mesh
      ! Input/output variables
      type(scalar_field) :: u_input, previous_u, iter_u
      type(scalar_field) :: u_output
      real :: theta, dt

      call c_f_pointer(context, matrices)
      
      call get_option(trim(adjoint_field_path("/material_phase::Fluid/scalar_field::Velocity/prognostic/temporal_discretisation/theta")), &
                                        & theta, default=0.5)
      call get_option("/timestepping/timestep", dt)

      call c_f_pointer(context, matrices)
      u_mesh => extract_mesh(matrices, "VelocityMesh")
      diffusion_mat => extract_csr_matrix(matrices, "DiffusionMatrix")
      mass_mat => extract_csr_matrix(matrices, "MassMatrix")

      if (nvar==2) then
        if (variables(1)%timestep==variables(2)%timestep-1) then
          call field_from_adj_vector(dependencies(1), previous_u)
          call field_from_adj_vector(dependencies(2), iter_u)
        else if(variables(2)%timestep==variables(1)%timestep-1) then
          call field_from_adj_vector(dependencies(2), previous_u)
          call field_from_adj_vector(dependencies(1), iter_u)
        else
          FLAbort("Dependencies have no contiguous timesteps.")
        end if
      else if (nvar==1) then
        call field_from_adj_vector(dependencies(1), previous_u)
        call field_from_adj_vector(dependencies(1), iter_u)
      else
        FLAbort("Unknown number of dependencies in timestepping_operator_assembly_callback")
      end if

      ! Ensure that the input velocity has indeed dirichlet bc attached
      ! Ensure that the input velocity has indeed dirichlet bc attached
      if (get_boundary_condition_count(previous_u)==0) then
        FLAbort("Velocity in callback must have dirichlet conditions attached")
      end if

      call field_from_adj_vector(input, u_input)
      
      call allocate(u_output, u_mesh, "TimeSteppingAdjointVelocityOutput")
      call zero(u_output)

      call allocate(timestepping_mat, mass_mat%sparsity, name="TimesteppingMatrix")
      call zero(timestepping_mat)

      if (.not. have_option(trim(adjoint_field_path('/material_phase::Fluid/scalar_field::Velocity/prognostic/temporal_discretisation/remove_time_term')))) then
        call addto(timestepping_mat, mass_mat, scale=1.0/dt)
      end if

!        call mult(tmp, advection_mat, u_input)
!        call scale(tmp, theta - 1.0)
!        call addto(u_output, tmp)

      call addto(timestepping_mat, diffusion_mat, scale=theta - 1.0)

      call mangle_dirichlet_rows(timestepping_mat, previous_u, keep_diag=.false.)

      if (hermitian==ADJ_FALSE) then
        call mult(u_output, timestepping_mat, u_input)
      else
        call mult_T(u_output, timestepping_mat, u_input)
      end if
      call scale(u_output, coefficient)
      output = field_to_adj_vector(u_output)
      call deallocate(u_output)
      call deallocate(timestepping_mat)
    end subroutine timestepping_operator_action_callback


#endif
end module burgers_adjoint_callbacks
