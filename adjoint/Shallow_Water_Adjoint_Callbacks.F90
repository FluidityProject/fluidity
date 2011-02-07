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
#if 0
    use libadjoint
    use state_module
    use fields
    implicit none

    contains

    subroutine identity_assembly_callback(nvar, variables, dependencies, hermitian, context, output, rhs) bind(c)
      use iso_c_binding
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      type(c_ptr), intent(in), value :: context
      type(adj_matrix), intent(out) :: output
      type(adj_vector), intent(out) :: rhs

      type(mesh_type), pointer :: mesh
      type(scalar_field) :: empty_pressure_field
      type(vector_field) :: empty_velocity_field

      ! context must be supplied as the address of the mesh_type we're to allocate rhs on
      call c_f_pointer(context, mesh)

      if (trim(mesh%name) == "PressureMesh") then
        call allocate(empty_pressure_field, mesh, trim("AdjointPressureRhs"))
        call zero(empty_pressure_field)
        rhs = field_to_adj_vector(empty_pressure_field)
        call deallocate(empty_pressure_field)
      else if (trim(mesh%name) == "VelocityMesh") then
        call allocate(empty_velocity_field, u%dim, mesh, trim("AdjointVelocityRhs"))
        call zero(empty_velocity_field)
        rhs = field_to_adj_vector(empty_velocity_field)
        call deallocate(empty_velocity_field)
      else 
        FLAbort ("Unknown mesh found in context of identity_assembly_callback")
      end if

      ! We're not actually going to do anything daft like allocate memory for
      ! an identity matrix. Instead, we'll fill it in with a special tag
      ! that we'll catch later on in the solution of the adjoint equations.
      output%ptr = c_null_ptr
      output%klass = IDENTITY_MATRIX

    end subroutine identity_assembly_callback



    subroutine wave_mat_assembly_callback(nvar, variables, dependencies, hermitian, context, output, rhs) bind(c)
      use iso_c_binding
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      type(c_ptr), intent(in), value :: context
      type(adj_matrix), intent(out) :: output
      type(adj_vector), intent(out) :: rhs

      type(mesh_type), pointer :: mesh
      type(scalar_field) :: empty_pressure_field
      type(csr_matrix) :: wave_mat_T

      ! context must be supplied as the address of the VelocityMesh
      call c_f_pointer(context, mesh)
      assert(trim(mesh%name) == "PressureMesh")

      call allocate(empty_pressure_field, mesh, trim("AdjointPressureRhs"))
      call zero(empty_pressure_field)
      rhs = field_to_adj_vector(empty_pressure_field)
      call deallocate(empty_pressure_field)

      if (hermitian == ADJ_FALSE) then
        output = matrix_to_adj_matrix(wave_mat)
      else if (hermitian == ADJ_TRUE .and. trim(mesh%name) == "PressureMesh") then
        wave_mat_T = transpose(wave_mat, symmetric_sparsity=.true.)
        output = matrix_to_adj_matrix(wave_mat_T)
        call deallocate(wave_mat_T)
      end if
    end subroutine wave_mat_assembly_callback

    subroutine coriolis_assembly_callback(nvar, variables, dependencies, hermitian, context, output, rhs) bind(c)
      use iso_c_binding
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      type(c_ptr), intent(in), value :: context
      type(adj_matrix), intent(out) :: output
      type(adj_vector), intent(out) :: rhs

      type(mesh_type), pointer :: mesh
      type(vector_field) :: empty_velocity_field
      type(scalar_field) :: empty_pressure_field
      type(block_csr_matrix) :: coriolis_mat_T

      ! context must be supplied as the address of the VelocityMesh
      call c_f_pointer(context, mesh)

      if (hermitian == ADJ_FALSE .and. trim(mesh%name) == "VelocityMesh") then
        call allocate(empty_velocity_field, u%dim, mesh, trim("AdjointVelocityRhs"))
        call zero(empty_velocity_field)
        rhs = field_to_adj_vector(empty_velocity_field)
        call deallocate(empty_velocity_field)
      
        output = matrix_to_adj_matrix(coriolis_mat)
      else if (hermitian == ADJ_TRUE .and. trim(mesh%name) == "PressureMesh") then
        call allocate(empty_pressure_field, mesh, trim("AdjointPressureRhs"))
        call zero(empty_pressure_field)
        rhs = field_to_adj_vector(empty_pressure_field)
        call deallocate(empty_pressure_field)

        coriolis_mat_T = transpose(coriolis_mat, symmetric_sparsity=.true.)
        output = matrix_to_adj_matrix(coriolis_mat_T)
        call deallocate(coriolis_mat_T)
      else 
        FLAbort ("Unknown mesh found in context of coriolis_assembly_callback")
      end if
    end subroutine coriolis_assembly_callback

    subroutine ggradperp_action_callback(nvar, variables, dependencies, hermitian, input, context, output) bind(c)
      use iso_c_binding
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      type(adj_vector), intent(in), value :: input
      type(c_ptr), intent(in), value :: context
      type(adj_vector), intent(out) :: output
      
      ! Variables for hermitian == ADJ_FALSE
      type(scalar_field) :: eta_input
      type(vector_field) :: u_output
      ! Variables for hermitian == ADJ_TRUE
      type(vector_field) :: u_input
      type(scalar_field) :: p_output
      
      type(vector_field) :: u_output_tmp
      type(scalar_field) :: u1, u2
      type(block_csr_matrix) :: inverse_coriolis_mat_T


      if (hermitian==ADJ_FALSE) then
        call field_from_adj_vector(input, eta_input)
        ! So, we'll mult_T with div_mat, then scale by g and then act inverse_coriolis_mat on it. 
        call allocate(u_output, u%dim, u%mesh, "AdjointGGradPerpOutput")
        call zero(u_output)

        call mult_T(u_output,div_mat,eta_input)
        call scale(u_output, -g)

        call mult(u_output, inverse_coriolis_mat, u_output)

        output = field_to_adj_vector(u_output)
        call deallocate(u_output)
      else
        ! Do the same steps as above, but backwards and with the transposed operators
        call field_from_adj_vector(input, u_input)
        call allocate(p_output, eta%mesh, "AdjointGGradPerpOutput")
        call allocate(u_output, u%dim, u%mesh, "AdjointGGradPerpOutput")
        call zero(u_output)
        call zero(p_output)

        inverse_coriolis_mat_T = transpose(inverse_coriolis_mat)
        call mult(u_output, inverse_coriolis_mat_T, u_input)
        call scale(u_output, -g)
        call mult(p_output, div_mat, u_output)

        output = field_to_adj_vector(p_output)
        call deallocate(p_output)
        call deallocate(u_output)
      end if
    end subroutine ggradperp_action_callback


    ! TODO
    subroutine div_bigmat_grad_action_callback(nvar, variables, dependencies, hermitian, input, context, output) bind(c)
      use iso_c_binding
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      type(adj_vector), intent(in), value :: input
      type(c_ptr), intent(in), value :: context
      type(adj_vector), intent(out) :: output
      
      ! Variables for hermitian == ADJ_FALSE
      type(scalar_field) :: eta_input
      type(vector_field) :: u_output
      ! Variables for hermitian == ADJ_TRUE
      type(vector_field) :: u_input
      type(scalar_field) :: p_output
      
      type(vector_field) :: u_output_tmp
      type(scalar_field) :: u1, u2
      type(block_csr_matrix) :: inverse_coriolis_mat_T


      
      if (hermitian==ADJ_FALSE) then
        call field_from_adj_vector(input, eta_input)
        ! So, we'll mult_T with div_mat, then scale by g and then act inverse_coriolis_mat on it. 
        call allocate(u_output, u%dim, u%mesh, "AdjointGGradPerpOutput")
        call zero(u_output)

        call mult_T(u_output,div_mat,eta_input)
        call scale(u_output, -g)

        call mult(u_output, inverse_coriolis_mat, u_output)

        output = field_to_adj_vector(u_output)
        call deallocate(u_output)
      else
        ! Do the same steps as above, but backwards and with the transposed operators
        call field_from_adj_vector(input, u_input)
        call allocate(p_output, eta%mesh, "AdjointGGradPerpOutput")
        call allocate(u_output, u%dim, u%mesh, "AdjointGGradPerpOutput")
        call zero(u_output)
        call zero(p_output)

        inverse_coriolis_mat_T = transpose(inverse_coriolis_mat)
        call mult(u_output, inverse_coriolis_mat_T, u_input)
        call scale(u_output, -g)
        call mult(p_output, div_mat, u_output)

        output = field_to_adj_vector(p_output)
        call deallocate(p_output)
        call deallocate(u_output)
      end if
    end subroutine div_bigmat_grad_action_callback


    subroutine bigmat_grad_action_callback(nvar, variables, dependencies, hermitian, input, context, output) bind(c)
      use iso_c_binding
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      type(adj_vector), intent(in), value :: input
      type(c_ptr), intent(in), value :: context
      type(adj_vector), intent(out) :: output
      
      ! Variables for hermitian == ADJ_FALSE
      type(scalar_field) :: eta_input
      type(vector_field) :: u_output
      ! Variables for hermitian == ADJ_TRUE
      type(vector_field) :: u_input
      type(scalar_field) :: p_output
      
      type(vector_field) :: u_output_tmp
      type(scalar_field) :: u1, u2
      type(block_csr_matrix) :: big_mat_T

      if (hermitian==ADJ_FALSE) then
        call field_from_adj_vector(input, eta_input)
        ! So, we'll mult_T with div_mat, and then act big_mat on it. 
        call allocate(u_output, u%dim, u%mesh, "AdjointBigMatGradOutput")
        call zero(u_output)

        call mult_T(u_output,div_mat,eta_input)
        call mult(u_output, big_mat, u_output)

        output = field_to_adj_vector(u_output)
        call deallocate(u_output)
      else
        ! Do the same steps as above, but backwards and with the transposed operators
        call field_from_adj_vector(input, u_input)
        call allocate(p_output, eta%mesh, "AdjointBigMatGradTOutput")
        call allocate(u_output, u%dim, u%mesh, "AdjointBigMatGradTOutput")
        call zero(u_output)
        call zero(p_output)

        call mult_T(u_output, big_mat, u_input)
        call mult(p_output, div_mat, u_output)

        output = field_to_adj_vector(p_output)
        call deallocate(p_output)
        call deallocate(u_output)
      end if
    end subroutine bigmat_grad_action_callback

    subroutine bigmat_coriolis_action_callback(nvar, variables, dependencies, hermitian, input, context, output) bind(c)
      use iso_c_binding
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      type(adj_vector), intent(in), value :: input
      type(c_ptr), intent(in), value :: context
      type(adj_vector), intent(out) :: output
      
      type(vector_field) :: u_output, u_input
      
      type(vector_field) :: u_output_tmp
      type(block_csr_matrix) :: coriolis_T

      call field_from_adj_vector(input, u_input)
      
      call allocate(u_output, u%dim, u%mesh, "AdjointBigMatCoriolisOutput")
      call zero(u_output)

      if (hermitian==ADJ_FALSE) then
        ! So, we'll mult with coriolis_mat. 
        call mult(u_output,coriolis_mat,u_input)
      else
        ! TODO
!        call mult_T(u_output,coriolis_mat,u_input)
      end if 
      output = field_to_adj_vector(u_output)
      call deallocate(u_output)
    end subroutine bigmat_coriolis_action_callback

#endif
#endif
end module shallow_water_adjoint_callbacks
