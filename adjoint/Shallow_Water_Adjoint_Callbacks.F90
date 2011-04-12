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

module shallow_water_adjoint_callbacks
#ifdef HAVE_ADJOINT
#include "libadjoint/adj_fortran.h"
    use libadjoint
    use libadjoint_data_callbacks
    use iso_c_binding
    use state_module
    use fields
    use sparse_matrices_fields
    use spud, only: get_option
    use adjoint_global_variables, only: adj_path_lookup
    use mangle_options_tree, only: adjoint_field_path
    use manifold_projections
    implicit none

    private

    public :: register_sw_operator_callbacks

    contains

    subroutine register_sw_operator_callbacks(adjointer)
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int) :: ierr

      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ASSEMBLY_CB, "CartesianVelocityMassMatrix", c_funloc(cartesian_velocity_mass_assembly_callback))
      call adj_chkierr(ierr)
      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ASSEMBLY_CB, "LocalVelocityMassMatrix", c_funloc(local_velocity_mass_assembly_callback))
      call adj_chkierr(ierr)
      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ASSEMBLY_CB, "LayerThicknessMassMatrix", c_funloc(layerthickness_mass_assembly_callback))
      call adj_chkierr(ierr)
      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ASSEMBLY_CB, "WaveMatrix", c_funloc(wave_mat_assembly_callback))
      call adj_chkierr(ierr)
      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ASSEMBLY_CB, "Coriolis", c_funloc(coriolis_assembly_callback))
      call adj_chkierr(ierr)

      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, "CartesianVelocityMassMatrix", c_funloc(cartesian_velocity_mass_action_callback))
      call adj_chkierr(ierr)
      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, "LocalVelocityMassMatrix", c_funloc(local_velocity_mass_action_callback))
      call adj_chkierr(ierr)
      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, "LayerThicknessMassMatrix", c_funloc(layerthickness_mass_action_callback))
      call adj_chkierr(ierr)
      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, "Grad", c_funloc(grad_action_callback))
      call adj_chkierr(ierr)
      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, "GradMinusDivBigMatCoriolis", c_funloc(grad_minus_div_bigmat_coriolis_action_callback))
      call adj_chkierr(ierr)
      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, "DivBigMatGrad", c_funloc(div_bigmat_grad_action_callback))
      call adj_chkierr(ierr)
      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, "BigMatGrad", c_funloc(bigmat_grad_action_callback))
      call adj_chkierr(ierr)
      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, "BigMatCoriolis", c_funloc(bigmat_coriolis_action_callback))
      call adj_chkierr(ierr)
      ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, "MassLocalProjection", c_funloc(local_projection_action_callback))
      call adj_chkierr(ierr)
    end subroutine register_sw_operator_callbacks

    subroutine cartesian_velocity_mass_assembly_callback(nvar, variables, dependencies, hermitian, coefficient, context, output, rhs) bind(c)
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(c_ptr), intent(in), value :: context
      type(adj_matrix), intent(out) :: output
      type(adj_vector), intent(out) :: rhs

      type(state_type), pointer :: matrices
      type(vector_field), pointer :: u
      type(vector_field) :: empty_velocity_field
      type(block_csr_matrix), pointer :: u_mass_mat

      if (coefficient /= 1.0) then
        FLAbort("cartesian_velocity_mass_assembly_callback requires that the coefficient is 1.0")
      end if

      call c_f_pointer(context, matrices)
      u => extract_vector_field(matrices, "VelocityDummy")

      call allocate(empty_velocity_field, u%dim, u%mesh, trim("AdjointVelocityRhs"))
      call zero(empty_velocity_field)
      rhs = field_to_adj_vector(empty_velocity_field)
      call deallocate(empty_velocity_field)

      ! If we had dirichlet conditions, we would have to care about hermitian being true or not
      u_mass_mat => extract_block_csr_matrix(matrices, "CartesianVelocityMassMatrix")
      output = matrix_to_adj_matrix(u_mass_mat)

    end subroutine cartesian_velocity_mass_assembly_callback

    subroutine local_velocity_mass_assembly_callback(nvar, variables, dependencies, hermitian, coefficient, context, output, rhs) bind(c)
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(c_ptr), intent(in), value :: context
      type(adj_matrix), intent(out) :: output
      type(adj_vector), intent(out) :: rhs

      type(state_type), pointer :: matrices
      type(vector_field), pointer :: u
      type(vector_field) :: empty_velocity_field
      type(block_csr_matrix), pointer :: u_mass_mat

      if (coefficient /= 1.0) then
        FLAbort("local_velocity_identity_assembly_callback requires that the coefficient is 1.0")
      end if

      call c_f_pointer(context, matrices)
      u => extract_vector_field(matrices, "LocalVelocityDummy")

      call allocate(empty_velocity_field, u%dim, u%mesh, trim("AdjointVelocityRhs"))
      call zero(empty_velocity_field)
      rhs = field_to_adj_vector(empty_velocity_field)
      call deallocate(empty_velocity_field)

      ! If we had dirichlet conditions, we would have to care about hermitian being true or not
      u_mass_mat => extract_block_csr_matrix(matrices, "LocalVelocityMassMatrix")
      output = matrix_to_adj_matrix(u_mass_mat)
    end subroutine local_velocity_mass_assembly_callback

    subroutine layerthickness_mass_assembly_callback(nvar, variables, dependencies, hermitian, coefficient, context, output, rhs) bind(c)
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(c_ptr), intent(in), value :: context
      type(adj_matrix), intent(out) :: output
      type(adj_vector), intent(out) :: rhs

      type(state_type), pointer :: matrices
      type(mesh_type), pointer :: eta_mesh
      type(scalar_field) :: empty_eta_field
      type(csr_matrix), pointer :: h_mass_mat

      if (coefficient /= 1.0) then
        FLAbort("layerthickness_identity_assembly_callback requires that the coefficient is 1.0")
      end if

      call c_f_pointer(context, matrices)
      eta_mesh => extract_mesh(matrices, "LayerThicknessMesh")

      call allocate(empty_eta_field, eta_mesh, trim("AdjointLayerThicknessRhs"))
      call zero(empty_eta_field)
      rhs = field_to_adj_vector(empty_eta_field)
      call deallocate(empty_eta_field)

      ! Again, if we had dirichlet conditions we'd need to care about hermitian or not
      h_mass_mat => extract_csr_matrix(matrices, "LayerThicknessMassMatrix")
      output = matrix_to_adj_matrix(h_mass_mat)

    end subroutine layerthickness_mass_assembly_callback

    subroutine wave_mat_assembly_callback(nvar, variables, dependencies, hermitian, coefficient, context, output, rhs) bind(c)
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(c_ptr), intent(in), value :: context
      type(adj_matrix), intent(out) :: output
      type(adj_vector), intent(out) :: rhs

      type(state_type), pointer :: matrices
      type(mesh_type), pointer :: eta_mesh
      type(scalar_field) :: empty_pressure_field
      type(csr_matrix), pointer :: wave_mat
      type(csr_matrix), pointer :: wave_mat_T

      if (coefficient /= 1.0) then
        FLAbort("The coefficient in wave_mat_assembly_callback has to be 1.0")
      end if

      call c_f_pointer(context, matrices)
      eta_mesh => extract_mesh(matrices, "LayerThicknessMesh")
      wave_mat => extract_csr_matrix(matrices, "WaveMatrix")
      wave_mat_T => extract_csr_matrix(matrices, "WaveMatrixTranspose")

      call allocate(empty_pressure_field, eta_mesh, trim("AdjointPressureRhs"))
      call zero(empty_pressure_field)
      rhs = field_to_adj_vector(empty_pressure_field)
      call deallocate(empty_pressure_field)

      if (hermitian == ADJ_FALSE) then
        output = matrix_to_adj_matrix(wave_mat)
      else if (hermitian == ADJ_TRUE) then
        output = matrix_to_adj_matrix(wave_mat_T)
      end if
      output%flags = 0
    end subroutine wave_mat_assembly_callback

    subroutine coriolis_assembly_callback(nvar, variables, dependencies, hermitian, coefficient, context, output, rhs) bind(c)
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(c_ptr), intent(in), value :: context
      type(adj_matrix), intent(out) :: output
      type(adj_vector), intent(out) :: rhs

      type(state_type), pointer :: matrices
      type(vector_field), pointer :: u
      type(vector_field) :: empty_velocity_field
      type(block_csr_matrix), pointer :: inv_coriolis_mat
      type(block_csr_matrix), pointer :: inv_coriolis_mat_T

      if (coefficient /= 1.0) then
        FLAbort("The coefficient in coriolis_assembly_callback has to be 1.0")
      end if

      call c_f_pointer(context, matrices)
      u => extract_vector_field(matrices, "LocalVelocityDummy")
      inv_coriolis_mat => extract_block_csr_matrix(matrices, "InverseCoriolisMatrix")
      inv_coriolis_mat_T => extract_block_csr_matrix(matrices, "InverseCoriolisMatrixTranspose")

      call allocate(empty_velocity_field, u%dim, u%mesh, trim("AdjointVelocityRhs"))
      call zero(empty_velocity_field)
      rhs = field_to_adj_vector(empty_velocity_field)
      call deallocate(empty_velocity_field)

      if (hermitian == ADJ_FALSE) then
        output = matrix_to_adj_matrix(inv_coriolis_mat)
        output%flags = ADJ_MATRIX_INVERTED
      else if (hermitian == ADJ_TRUE) then
        output = matrix_to_adj_matrix(inv_coriolis_mat_T)
        output%flags = ADJ_MATRIX_INVERTED
      end if
    end subroutine coriolis_assembly_callback

    subroutine cartesian_velocity_mass_action_callback(nvar, variables, dependencies, hermitian, coefficient, input, context, output) bind(c)
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(adj_vector), intent(in), value :: input
      type(c_ptr), intent(in), value :: context
      type(adj_vector), intent(out) :: output

      type(state_type), pointer :: matrices
      type(vector_field), pointer :: u

      type(vector_field) :: u_input
      type(vector_field) :: u_output
      type(block_csr_matrix), pointer :: u_mass_mat

      call c_f_pointer(context, matrices)
      call field_from_adj_vector(input, u_input)
      call allocate(u_output, u_input%dim, u_input%mesh, "CartesianVelocityMassOutput")
      u_mass_mat => extract_block_csr_matrix(matrices, "CartesianVelocityMassMatrix")
      call mult(u_output, u_mass_mat, u_input)
      call scale(u_output, coefficient)
      output = field_to_adj_vector(u_output)
      call deallocate(u_output)
    end subroutine cartesian_velocity_mass_action_callback


    subroutine local_velocity_mass_action_callback(nvar, variables, dependencies, hermitian, coefficient, input, context, output) bind(c)
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(adj_vector), intent(in), value :: input
      type(c_ptr), intent(in), value :: context
      type(adj_vector), intent(out) :: output

      type(state_type), pointer :: matrices
      type(vector_field), pointer :: u

      type(vector_field) :: u_input
      type(vector_field) :: u_output
      type(block_csr_matrix), pointer :: u_mass_mat

      call c_f_pointer(context, matrices)
      call field_from_adj_vector(input, u_input)
      call allocate(u_output, u_input%dim, u_input%mesh, "LocalVelocityMassOutput")
      u_mass_mat => extract_block_csr_matrix(matrices, "LocalVelocityMassMatrix")
      call mult(u_output, u_mass_mat, u_input)
      call scale(u_output, coefficient)
      output = field_to_adj_vector(u_output)
      call deallocate(u_output)
    end subroutine local_velocity_mass_action_callback

    subroutine layerthickness_mass_action_callback(nvar, variables, dependencies, hermitian, coefficient, input, context, output) bind(c)
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
      type(csr_matrix), pointer :: h_mass_mat

      call c_f_pointer(context, matrices)
      eta_mesh => extract_mesh(matrices, "LayerThicknessMesh")

      call field_from_adj_vector(input, eta_input)
      call allocate(eta_output, eta_mesh, "LayerThicknessMassOutput")
      h_mass_mat => extract_csr_matrix(matrices, "LayerThicknessMassMatrix")
      call mult(eta_output, h_mass_mat, eta_input)
      call scale(eta_output, coefficient)
      output = field_to_adj_vector(eta_output)
      call deallocate(eta_output)
    end subroutine layerthickness_mass_action_callback

    subroutine grad_action_callback(nvar, variables, dependencies, hermitian, coefficient, input, context, output) bind(c)
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(adj_vector), intent(in), value :: input
      type(c_ptr), intent(in), value :: context
      type(adj_vector), intent(out) :: output

      type(state_type), pointer :: matrices
      type(block_csr_matrix), pointer :: div_mat, inverse_coriolis_mat
      type(vector_field), pointer :: u
      type(mesh_type), pointer :: eta_mesh
      ! Variables for hermitian == ADJ_FALSE
      type(scalar_field) :: eta_input
      type(vector_field) :: u_output
      ! Variables for hermitian == ADJ_TRUE
      type(vector_field) :: u_input
      type(scalar_field) :: eta_output

      call c_f_pointer(context, matrices)
      u => extract_vector_field(matrices, "LocalVelocityDummy")
      eta_mesh => extract_mesh(matrices, "LayerThicknessMesh")
      div_mat => extract_block_csr_matrix(matrices, "DivergenceMatrix")

      if (hermitian==ADJ_FALSE) then
        call field_from_adj_vector(input, eta_input)
        ! So, we'll mult_T with div_mat.
        call allocate(u_output, u%dim, u%mesh, "AdjointGradOutput")
        call zero(u_output)

        call mult_T(u_output, div_mat, eta_input)
        call scale(u_output, coefficient)

        output = field_to_adj_vector(u_output)
        call deallocate(u_output)
      else
        ! Do the same steps as above, but backwards and with the transposed operators
        call field_from_adj_vector(input, u_input)
        call allocate(eta_output, eta_mesh, "AdjointGradOutput")
        call zero(eta_output)
        ! We use u_output as a temporay variable
        call allocate(u_output, u%dim, u%mesh, "AdjointGradOutputTemporaryVelocityField")
        call zero(u_output)

        call mult(eta_output, div_mat, u_output)
        call scale(eta_output, coefficient)

        output = field_to_adj_vector(eta_output)
        call deallocate(eta_output)
        call deallocate(u_output)
      end if
    end subroutine grad_action_callback

    subroutine grad_minus_div_bigmat_coriolis_action_callback(nvar, variables, dependencies, hermitian, coefficient, input, context, output) bind(c)
      use global_parameters, only:  dt
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
      type(vector_field), pointer :: u
      type(block_csr_matrix), pointer :: big_mat, div_mat, coriolis_mat

      ! hermitian == ADJ_FALSE
      type(scalar_field) :: eta_output
      type(vector_field) :: u_input
      ! hermitian == ADJ_TRUE
      type(scalar_field) :: eta_input
      type(vector_field) :: u_output

      type(vector_field) :: u_tmp
      real :: theta, d0
      character(len=ADJ_DICT_LEN) :: path
      integer :: ierr

      if (coefficient /= -1.0) then
        FLAbort("The coefficient in grad_minus_div_bigmat_coriolis_action_callback has to be -1.0")
      end if

      call c_f_pointer(context, matrices)
      u => extract_vector_field(matrices, "LocalVelocityDummy")
      eta_mesh => extract_mesh(matrices, "LayerThicknessMesh")
      big_mat => extract_block_csr_matrix(matrices, "InverseBigMatrix")
      div_mat => extract_block_csr_matrix(matrices, "DivergenceMatrix")
      coriolis_mat => extract_block_csr_matrix(matrices, "CoriolisMatrix")

      ierr = adj_dict_find(adj_path_lookup, "Fluid::LayerThickness", path)
      call adj_chkierr(ierr)
      path = adjoint_field_path(path)

      call get_option(trim(path) // "/prognostic/temporal_discretisation/theta", theta)
      call get_option(trim(path) // "/prognostic/mean_layer_thickness", d0)

      if (hermitian==ADJ_FALSE) then
        call field_from_adj_vector(input, u_input)
        call allocate(u_tmp, u%dim, u%mesh, "TemporaryVelocityVariable")
        call zero(u_tmp)
        call allocate(eta_output, eta_mesh, "AdjointGradMinusDivBigmatCoriolisOutput")
        call zero(eta_output)

        call mult(u_tmp, coriolis_mat, u_input)
        call mult(u_tmp, big_mat, u_tmp)
        call scale(u_tmp, -1.0*dt*theta)
        call addto(u_tmp, u_input)
        call mult(eta_output, div_mat, u_tmp)
        call scale(eta_output, dt*d0)

        output = field_to_adj_vector(eta_output)
        call deallocate(eta_output)
        call deallocate(u_tmp)
      else
        ! Do the same steps as above, but backwards and with the transposed operators
        call field_from_adj_vector(input, eta_input)
        call allocate(u_tmp, u%dim, u%mesh, "TemporaryVelocityVariable")
        call zero(u_tmp)
        call allocate(u_output, u%dim, u%mesh, "AdjointGradMinusDivBigmatCoriolisOutput")
        call zero(u_output)

        call mult_T(u_tmp, div_mat, eta_input)
        call mult_T(u_output, big_mat, u_tmp)
        call mult_T(u_output, coriolis_mat, u_output)
        call scale(u_output, -1.0*dt*theta)
        call addto(u_output, u_tmp)
        call scale(u_output, dt*d0)

        output = field_to_adj_vector(u_output)
        call deallocate(u_output)
        call deallocate(u_tmp)
      end if
    end subroutine grad_minus_div_bigmat_coriolis_action_callback

    subroutine div_bigmat_grad_action_callback(nvar, variables, dependencies, hermitian, coefficient, input, context, output) bind(c)
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
      type(vector_field), pointer :: u
      type(block_csr_matrix), pointer :: big_mat, div_mat

      type(scalar_field) :: eta_input
      type(scalar_field) :: eta_output
      type(vector_field) :: u_tmp

      call c_f_pointer(context, matrices)
      u => extract_vector_field(matrices, "LocalVelocityDummy")
      eta_mesh => extract_mesh(matrices, "LayerThicknessMesh")
      big_mat => extract_block_csr_matrix(matrices, "InverseBigMatrix")
      div_mat => extract_block_csr_matrix(matrices, "DivergenceMatrix")

      call field_from_adj_vector(input, eta_input)
      call allocate(u_tmp, u%dim, u%mesh, "TemporaryVelocityVariable")
      call zero(u_tmp)
      call allocate(eta_output, eta_mesh, "AdjointDivBigmatGradOutput")
      call zero(eta_output)

      if (hermitian==ADJ_FALSE) then
        call mult_T(u_tmp, div_mat, eta_input)
        call mult(u_tmp, big_mat, u_tmp)
        call mult(eta_output, div_mat, u_tmp)
        call scale(eta_output, coefficient)
        output = field_to_adj_vector(eta_output)
      else
        ! Do the same steps as above, but backwards and with the transposed operators
        call mult_T(u_tmp, div_mat, eta_input)
        call mult_T(u_tmp, big_mat, u_tmp)
        call mult(eta_output, div_mat, u_tmp)
        call scale(eta_output, coefficient)
        output = field_to_adj_vector(eta_output)
      end if
        call deallocate(eta_output)
        call deallocate(u_tmp)
    end subroutine div_bigmat_grad_action_callback

    subroutine bigmat_grad_action_callback(nvar, variables, dependencies, hermitian, coefficient, input, context, output) bind(c)
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(adj_vector), intent(in), value :: input
      type(c_ptr), intent(in), value :: context
      type(adj_vector), intent(out) :: output

      ! Variables for hermitian == ADJ_FALSE
      type(scalar_field) :: eta_input
      type(vector_field) :: u_output
      ! Variables for hermitian == ADJ_TRUE
      type(vector_field) :: u_input
      type(scalar_field) :: eta_output

      type(state_type), pointer :: matrices
      type(mesh_type), pointer :: eta_mesh
      type(vector_field), pointer :: u
      type(block_csr_matrix), pointer :: big_mat, div_mat


      call c_f_pointer(context, matrices)
      u => extract_vector_field(matrices, "LocalVelocityDummy")
      eta_mesh => extract_mesh(matrices, "LayerThicknessMesh")
      big_mat => extract_block_csr_matrix(matrices, "InverseBigMatrix")
      div_mat => extract_block_csr_matrix(matrices, "DivergenceMatrix")

      if (hermitian==ADJ_FALSE) then
        call field_from_adj_vector(input, eta_input)
        ! So, we'll mult_T with div_mat, and then act big_mat on it.
        call allocate(u_output, u%dim, u%mesh, "AdjointBigMatGradOutput")
        call zero(u_output)

        call mult_T(u_output, div_mat, eta_input)
        call mult(u_output, big_mat, u_output)
        call scale(u_output, coefficient)

        output = field_to_adj_vector(u_output)
        call deallocate(u_output)
      else
        ! Do the same steps as above, but backwards and with the transposed operators
        call field_from_adj_vector(input, u_input)
        call allocate(eta_output, eta_mesh, "AdjointBigMatGradOutput")
        call allocate(u_output, u%dim, u%mesh, "AdjointBigMatGradTOutput_TemporaryVelocityField")
        call zero(u_output)
        call zero(eta_output)

        call mult_T(u_output, big_mat, u_input)
        call mult(eta_output, div_mat, u_output)
        call scale(eta_output, coefficient)

        output = field_to_adj_vector(eta_output)
        call deallocate(eta_output)
        call deallocate(u_output)
      end if
    end subroutine bigmat_grad_action_callback

    subroutine bigmat_coriolis_action_callback(nvar, variables, dependencies, hermitian, coefficient, input, context, output) bind(c)
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(adj_vector), intent(in), value :: input
      type(c_ptr), intent(in), value :: context
      type(adj_vector), intent(out) :: output

      type(vector_field) :: u_output, u_input

      type(state_type), pointer :: matrices
      type(vector_field), pointer :: u
      type(mesh_type), pointer :: eta_mesh
      type(block_csr_matrix), pointer :: coriolis_mat, big_mat

      call c_f_pointer(context, matrices)
      u => extract_vector_field(matrices, "LocalVelocityDummy")
      eta_mesh => extract_mesh(matrices, "LayerThicknessMesh")
      coriolis_mat => extract_block_csr_matrix(matrices, "CoriolisMatrix")
      big_mat => extract_block_csr_matrix(matrices, "InverseBigMatrix")

      call field_from_adj_vector(input, u_input)
      call allocate(u_output, u%dim, u%mesh, "AdjointBigMatCoriolisOutput")
      call zero(u_output)

      if (hermitian==ADJ_FALSE) then
        ! So, we'll mult with coriolis_mat and big_mat.
        call mult(u_output, coriolis_mat, u_input)
        call mult(u_output, big_mat, u_output)
      else
        call mult_T(u_output, big_mat, u_input)
        call mult_T(u_output, coriolis_mat, u_output)
      end if
      call scale(u_output, coefficient)
      output = field_to_adj_vector(u_output)
      call deallocate(u_output)
    end subroutine bigmat_coriolis_action_callback

    subroutine local_projection_action_callback(nvar, variables, dependencies, hermitian, coefficient, input, context, output) bind(c)
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(adj_vector), intent(in), value :: input
      type(c_ptr), intent(in), value :: context
      type(adj_vector), intent(out) :: output

      type(vector_field) :: u_output, u_input, tmp_u

      type(state_type), pointer :: matrices
      type(vector_field), pointer :: X, local_dummy_u, dummy_u
      type(mesh_type), pointer :: eta_mesh
      type(block_csr_matrix), pointer :: coriolis_mat, big_mat, cartesian_mass_matrix, local_mass_matrix

      call c_f_pointer(context, matrices)
      local_dummy_u => extract_vector_field(matrices, "LocalVelocityDummy")
      dummy_u => extract_vector_field(matrices, "VelocityDummy")
      X => extract_vector_field(matrices, "Coordinate")
      cartesian_mass_matrix => extract_block_csr_matrix(matrices, "CartesianVelocityMassMatrix")
      local_mass_matrix => extract_block_csr_matrix(matrices, "LocalVelocityMassMatrix")

      call field_from_adj_vector(input, u_input)

      if (hermitian==ADJ_FALSE) then
        call allocate(tmp_u, local_dummy_u%dim, local_dummy_u%mesh, "LocalProjectionTemp")
        call allocate(u_output, local_dummy_u%dim, local_dummy_u%mesh, "LocalProjectionOutput")
        call zero(u_output)
        ! So, we'll project from cartesian space to local space
        call project_cartesian_to_local(X, u_input, tmp_u)
        call mult(u_output, local_mass_matrix, tmp_u)
        call deallocate(tmp_u)
      else
        call allocate(tmp_u, dummy_u%dim, dummy_u%mesh, "CartesianProjectionTemp")
        call allocate(u_output, dummy_u%dim, dummy_u%mesh, "CartesianProjectionOutput")
        call zero(u_output)
        ! So, we'll project from local space to cartesian space
        call project_local_to_cartesian(X, u_input, tmp_u)
        call mult(u_output, cartesian_mass_matrix, tmp_u)
        call deallocate(tmp_u)
      end if
      call scale(u_output, coefficient)
      output = field_to_adj_vector(u_output)
      call deallocate(u_output)
    end subroutine local_projection_action_callback


#endif
end module shallow_water_adjoint_callbacks
