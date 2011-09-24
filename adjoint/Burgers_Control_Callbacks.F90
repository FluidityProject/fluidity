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

module burgers_adjoint_controls
  use fields
  use global_parameters, only : PYTHON_FUNC_LEN, OPTION_PATH_LEN
  use spud
  use state_module
  use python_state
  use sparse_matrices_fields
  use manifold_projections
  use mangle_dirichlet_rows_module

    implicit none

    private

    public :: burgers_adjoint_timestep_callback

    contains 
    
    ! Data is a void pointer used to pass variables into the callback
    subroutine burgers_adjoint_timestep_callback(states, timestep, functional_name, data)
      use iso_c_binding, only: c_ptr
      use global_parameters, only: OPTION_PATH_LEN
      use state_module
      type(state_type), dimension(:), intent(inout) :: states
      integer, intent(in) :: timestep
      character(len=*), intent(in) :: functional_name
      type(c_ptr), intent(in) :: data
      
      type(state_type), pointer :: matrices
      character(len=OPTION_PATH_LEN) :: field_name, control_type, control_deriv_name, field_deriv_name, name, material_phase_name
      integer :: nb_controls
      integer :: i, state_id, s_idx
      type(scalar_field), pointer :: sfield, adj_sfield, adj_sfield2
      type(vector_field), pointer :: vfield, adj_vfield, positions
      type(tensor_field), pointer :: tfield, adj_tfield
      type(block_csr_matrix), pointer :: big_mat, div_mat, u_mass_mat
      type(csr_matrix), pointer :: h_mass_mat
      type(csr_matrix) :: h_mangled_mass_mat
      type(state_type), dimension(1) :: dummy_state
      
      type(scalar_field) :: dummy_u
      type(vector_field) :: local_src_tmp, local_src_tmp2, src_tmp
      real :: theta, d0, dt 
      integer :: stat
     
      ! Cast the data to the matrices state 
      call c_f_pointer(data, matrices)
      

      nb_controls = option_count("/adjoint/controls/control")
      do i = 0, nb_controls-1
        call get_option("/adjoint/controls/control[" // int2str(i) //"]/name", control_deriv_name)
        call get_option("/adjoint/controls/control[" // int2str(i) //"]/type/name", control_type)
        call get_option("/adjoint/controls/control[" // int2str(i) //"]/type::" // trim(control_type) // "/field_name", name)
        s_idx = scan(trim(name), "::")
        if (s_idx == 0) then 
          FLAbort("The control " // trim(control_deriv_name) // " uses an invalid field_name. It should be of the form Materialphase::Fieldname")
        end if
        material_phase_name = name(1:s_idx - 1)
        field_name = name(s_idx + 2:len_trim(name))
        ! Find state associated with the material phase
        do state_id = 1, size(states) 
          if (trim(states(state_id)%name) == trim(material_phase_name)) then
            exit
          end if
        end do
        ! Make sure we found the state
        if (.not. trim(states(state_id)%name) == trim(material_phase_name)) then
          FLAbort("Could not find state " // trim(material_phase_name) // " as specified in control " // trim(control_deriv_name) // ".")
        end if

        control_deriv_name = trim(control_deriv_name) // "_TotalDerivative"
        select case(trim(control_type))

          !!!!!!!!!!!!! Initial condition !!!!!!!!!!!!
          case ("initial_condition")
            if (timestep /= 0) then
             cycle
            end if 
            field_deriv_name = trim(functional_name) // "_" // control_deriv_name 
            if (has_scalar_field(states(state_id), field_deriv_name)) then
              if (trim(field_name) == "Velocity") then
                adj_sfield => extract_scalar_field(states(state_id), "Adjoint" // trim(field_name))
                sfield => extract_scalar_field(states(state_id), field_deriv_name) ! Output field
                positions => extract_vector_field(matrices, "Coordinate")
                call set(sfield, adj_sfield)

                ! Populate boundary conditions
                call allocate(dummy_u, adj_sfield%mesh, "DummyVelocity")
                dummy_u%option_path =  "/material_phase::Fluid/scalar_field::AdjointVelocity"
                call insert(dummy_state(1), positions, "Coordinate")
                call insert(dummy_state(1), dummy_u, "AdjointVelocity")
                call populate_boundary_conditions(dummy_state)

                ! Now we need to zero the derivative field everywhere we have Dirichlet conditions
                call mangle_dirichlet_rows(sfield, dummy_u)

                call deallocate(dummy_u)
                call deallocate(dummy_state(1))
              else
                FLAbort("Sorry, I do not know how to compute the intial condition control for " // trim(field_name) // ".")
              end if
            else
              FLAbort("The control derivative field " // trim(field_deriv_name) // " specified for " // trim(control_deriv_name) // " is not a field in the state.")
            end if

          !!!!!!!!!!!!! Boundary condition !!!!!!!!!!!!
          case("boundary_condition")
            FLAbort("Boundary condition control not implemented yet.")

          !!!!!!!!!!!!!       Source       !!!!!!!!!!!!
          case("source_term")
            field_deriv_name = trim(functional_name) // "_" // control_deriv_name 
            sfield => extract_scalar_field(states(state_id), field_deriv_name, stat=stat) ! Output field
            if (stat /= 0) then
              FLAbort("The control derivative field " // trim(field_deriv_name) // " specified for " // trim(control_deriv_name) // " is not a field in the state.")
            end if

            if (timestep == 0) then
             call zero(sfield)
            else
              if (trim(field_name) == "VelocitySource") then
                adj_sfield => extract_scalar_field(states(state_id), "AdjointVelocity")
                h_mass_mat => extract_csr_matrix(matrices, "MassMatrix")
                positions => extract_vector_field(matrices, "Coordinate")
                
                ! Populate boundary conditions
                call allocate(dummy_u, adj_sfield%mesh, "DummyVelocity")
                dummy_u%option_path =  "/material_phase::Fluid/scalar_field::AdjointVelocity"
                call insert(dummy_state(1), positions, "Coordinate")
                call insert(dummy_state(1), dummy_u, "AdjointVelocity")
                call populate_boundary_conditions(dummy_state)

                ! Mangle the mass matrix
                call allocate(h_mangled_mass_mat, h_mass_mat%sparsity, name="MangledMassMatrix")
                call set(h_mangled_mass_mat, h_mass_mat)
                call mangle_dirichlet_rows(h_mangled_mass_mat, dummy_u, keep_diag=.false.)
                call mult_T(sfield, h_mangled_mass_mat, adj_sfield)

                ! Deallocate temporary variables
                call deallocate(dummy_state(1))
                call deallocate(dummy_u)
                call deallocate(h_mangled_mass_mat)

              else
                FLAbort("Sorry, I do not know how to compute the source condition control for " // trim(field_name) // ".")
              end if
            end if
        end select
      end do

    end subroutine burgers_adjoint_timestep_callback 

end module burgers_adjoint_controls
