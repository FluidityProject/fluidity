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

module shallow_water_adjoint_controls
  use fields
  use global_parameters, only : PYTHON_FUNC_LEN, OPTION_PATH_LEN
  use spud
  use state_module
  use python_state
  use sparse_matrices_fields

    implicit none

    private

    public :: shallow_water_adjoint_timestep_callback 

    contains 
    
    ! Data is a void pointer used to pass variables into the callback
    subroutine shallow_water_adjoint_timestep_callback(states, timestep, functional_name, data)
      use global_parameters, only: OPTION_PATH_LEN
      use state_module
      type(state_type), dimension(:), intent(inout) :: states
      integer, intent(in) :: timestep
      character(len=*), intent(in) :: functional_name
      character(len=1), intent(inout) :: data(:)
      
      type(state_type) :: matrices
      type(state_type) :: state
      character(len=OPTION_PATH_LEN) :: field_name, control_type, control_deriv_name, field_deriv_name
      integer :: nb_controls
      integer :: i
      type(scalar_field), pointer :: sfield, adj_sfield
      type(vector_field), pointer :: vfield, adj_vfield
      type(tensor_field), pointer :: tfield, adj_tfield
      type(block_csr_matrix), pointer :: u_mass_mat
      type(csr_matrix), pointer :: h_mass_mat
     
      state = states(1)
      !print *, "In shallow_water_adjoint_timestep_callback"
      !call print_state(state)
      ! Cast the data to the matrices state 
      matrices = transfer(data, matrices)
      !print*, " ************************* Matrices: ************************** "
      !call print_state(matrices) 
      nb_controls = option_count("/adjoint/controls/control")
      do i = 0, nb_controls-1
        call get_option("/adjoint/controls/control[" // int2str(i) //"]/name", control_deriv_name)
        control_deriv_name = trim(control_deriv_name) // "_TotalDerivative"
        call get_option("/adjoint/controls/control[" // int2str(i) //"]/type/name", control_type)
        call get_option("/adjoint/controls/control[" // int2str(i) //"]/type::" // trim(control_type) // "/field_name", field_name)
        select case(trim(control_type))
          !!!!!!!!!!!!! Initial condition !!!!!!!!!!!!
          case ("initial_condition")
            if (timestep /= 0) then
             cycle
            end if 
            field_deriv_name = trim(functional_name) // "_" // control_deriv_name 
            if (has_scalar_field(state, field_deriv_name)) then
              if (trim(field_name) == "LayerThickness") then
                adj_sfield => extract_scalar_field(state, "Adjoint" // trim(field_name))
                sfield => extract_scalar_field(state, field_deriv_name)
                h_mass_mat => extract_csr_matrix(matrices, "LayerThicknessMassMatrix")
                call mult(sfield, h_mass_mat, adj_sfield)
              else
                FLAbort("Sorry, I do not know how to compute the intial condition control for " // trim(field_name) // ".")
              end if
            elseif (has_vector_field(state, field_deriv_name)) then
              if (trim(field_name) == "Velocity") then
                adj_vfield => extract_vector_field(state, "Adjoint" // trim(field_name))
                vfield => extract_vector_field(state, field_deriv_name)
                u_mass_mat => extract_block_csr_matrix(matrices, "CartesianVelocityMassMatrix")
                call mult(vfield, u_mass_mat, adj_vfield)
              else
                FLAbort("Sorry, I do not know how to compute the intial condition control for " // trim(field_name) // ".")
              end if
            elseif (has_tensor_field(state, field_deriv_name)) then
              tfield => extract_tensor_field(state, field_deriv_name)
              FLAbort("Sorry, I do not know how to compute the intial condition control for " // trim(field_name) // ".")
            else
              FLAbort("The control derivative field " // trim(field_deriv_name) // " specified for " // trim(control_deriv_name) // " is not a field in the state.")
            end if
          !!!!!!!!!!!!! Boundary condition !!!!!!!!!!!!
          case("boundary_condition")
            FLAbort("Boundary condition control not implemented yet.")
          !!!!!!!!!!!!! Source  !!!!!!!!!!!!
          case("source_term")
            FLAbort("Source control not implemented yet.")
        end select
      end do
      
      ! Cast back
      data = transfer(matrices, data)

    end subroutine shallow_water_adjoint_timestep_callback 

end module shallow_water_adjoint_controls
