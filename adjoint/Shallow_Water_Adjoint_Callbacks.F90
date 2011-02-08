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
    implicit none

    integer, parameter :: IDENTITY_MATRIX = 30

    contains

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
      type(vector_field), pointer :: u
      type(vector_field) :: empty_velocity_field

      if (coefficient /= 1.0) then
        FLAbort("The coefficient in velocity_identity_assembly_callback has to be 1.0")
      end if

      ! context must be supplied as the address of the mesh_type we're to allocate rhs on
      call c_f_pointer(context, matrices)
      u => extract_vector_field(matrices, "VelocityDummy")

      call allocate(empty_velocity_field, u%dim, u%mesh, trim("AdjointVelocityRhs"))
      call zero(empty_velocity_field)
      rhs = field_to_adj_vector(empty_velocity_field)
      call deallocate(empty_velocity_field)

      ! We're not actually going to do anything daft like allocate memory for
      ! an identity matrix. Instead, we'll fill it in with a special tag
      ! that we'll catch later on in the solution of the adjoint equations.
      output%ptr = c_null_ptr
      output%klass = IDENTITY_MATRIX

    end subroutine velocity_identity_assembly_callback

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
      type(csr_matrix) :: wave_mat_T

      if (coefficient /= 1.0) then
        FLAbort("The coefficient in wave_mat_assembly_callback has to be 1.0")
      end if

      call c_f_pointer(context, matrices)
      eta_mesh => extract_mesh(matrices, "LayerThicknessMesh")
      wave_mat => extract_csr_matrix(matrices, "WaveMatrix")

      call allocate(empty_pressure_field, eta_mesh, trim("AdjointPressureRhs"))
      call zero(empty_pressure_field)
      rhs = field_to_adj_vector(empty_pressure_field)
      call deallocate(empty_pressure_field)

      if (hermitian == ADJ_FALSE) then
        output = matrix_to_adj_matrix(wave_mat)
      else if (hermitian == ADJ_TRUE) then
        wave_mat_T = transpose(wave_mat, symmetric_sparsity=.true.)
        output = matrix_to_adj_matrix(wave_mat_T)
        call deallocate(wave_mat_T)
      end if
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
      type(block_csr_matrix), pointer :: coriolis_mat
      type(block_csr_matrix) :: coriolis_mat_T

      if (coefficient /= 1.0) then
        FLAbort("The coefficient in coriolis_assembly_callback has to be 1.0")
      end if

      call c_f_pointer(context, matrices)
      u => extract_vector_field(matrices, "VelocityDummy")
      coriolis_mat => extract_block_csr_matrix(matrices, "CoriolisMatrix")

      call allocate(empty_velocity_field, u%dim, u%mesh, trim("AdjointVelocityRhs"))
      call zero(empty_velocity_field)
      rhs = field_to_adj_vector(empty_velocity_field)
      call deallocate(empty_velocity_field)
      
      if (hermitian == ADJ_FALSE) then
        output = matrix_to_adj_matrix(coriolis_mat)
      else if (hermitian == ADJ_TRUE) then
        coriolis_mat_T = transpose(coriolis_mat, symmetric_sparsity=.true.)
        output = matrix_to_adj_matrix(coriolis_mat_T)
        call deallocate(coriolis_mat_T)
      end if
    end subroutine coriolis_assembly_callback

    subroutine gradperp_action_callback(nvar, variables, dependencies, hermitian, coefficient, input, context, output) bind(c)
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
      u => extract_vector_field(matrices, "VelocityDummy")
      eta_mesh => extract_mesh(matrices, "LayerThicknessMesh")
      div_mat => extract_block_csr_matrix(matrices, "DivergenceMatrix")
      inverse_coriolis_mat => extract_block_csr_matrix(matrices, "InverseCoriolisMatrix")

      if (hermitian==ADJ_FALSE) then
        call field_from_adj_vector(input, eta_input)
        ! So, we'll mult_T with div_mat, then scale by g and then act inverse_coriolis_mat on it. 
        call allocate(u_output, u%dim, u%mesh, "AdjointGradPerpOutput")
        call zero(u_output)

        call mult_T(u_output, div_mat, eta_input)
        call mult(u_output, inverse_coriolis_mat, u_output)
        call scale(u_output, coefficient)

        output = field_to_adj_vector(u_output)
        call deallocate(u_output)
      else
        ! Do the same steps as above, but backwards and with the transposed operators
        call field_from_adj_vector(input, u_input)
        call allocate(eta_output, eta_mesh, "AdjointGradPerpOutput")
        call zero(eta_output)
        ! We use u_output as a temporay variable
        call allocate(u_output, u%dim, u%mesh, "AdjointGradPerpOutputTemporaryVelocityField")
        call zero(u_output)

        call mult_T(u_output, inverse_coriolis_mat, u_input)
        call mult(eta_output, div_mat, u_output)
        call scale(eta_output, coefficient)

        output = field_to_adj_vector(eta_output)
        call deallocate(eta_output)
        call deallocate(u_output)
      end if
    end subroutine gradperp_action_callback

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
      
      if (coefficient /= 1.0) then
        FLAbort("The coefficient in grad_minus_div_bigmat_coriolis_action_callback has to be 1.0")
      end if
      
      call c_f_pointer(context, matrices)
      u => extract_vector_field(matrices, "VelocityDummy")
      eta_mesh => extract_mesh(matrices, "LayerThicknessMesh")
      big_mat => extract_block_csr_matrix(matrices, "InverseBigMatrix")
      div_mat => extract_block_csr_matrix(matrices, "DivergenceMatrix")
      coriolis_mat => extract_block_csr_matrix(matrices, "CoriolisMatrix")

      call get_option("/material_phase::Fluid/scalar_field::LayerThickness/prognostic/temporal_discretisation/theta", theta) 
      call get_option("/material_phase::Fluid/scalar_field::LayerThickness/ &
                      & prognostic/mean_layer_thickness", d0)
                      
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
      u => extract_vector_field(matrices, "VelocityDummy")
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
      u => extract_vector_field(matrices, "VelocityDummy")
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
      u => extract_vector_field(matrices, "VelocityDummy")
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

#endif
end module shallow_water_adjoint_callbacks
