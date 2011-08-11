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

module adjoint_controls
  use fields
  use global_parameters, only : PYTHON_FUNC_LEN, OPTION_PATH_LEN
  use boundary_conditions
  use spud
  use state_module
  use python_state

    implicit none

    private

    public :: adjoint_load_controls, adjoint_write_controls, adjoint_write_control_derivatives, allocate_and_insert_control_derivative_fields

    contains 

    ! Retrieves the control details from control number control_no in the optin tree. The outputs are:
    subroutine get_control_details(states, timestep, control_no, control_name, state_id, field_name, field_type, &
                                 & active, have_lb, lb_state_id, lb_material_phase_name, lb_field_name, have_ub, &
                                 & ub_state_id, ub_material_phase_name, ub_field_name)

       ! Inputs:
      type(state_type), dimension(:), intent(in) :: states
      integer, intent(in) :: control_no, timestep ! control_no: The index of the control in the option tree
                                                  ! timestep: The timestep

      ! Outputs:
      character(len=OPTION_PATH_LEN), intent(out) :: control_name, field_type, field_name ! control_name: the name of the control 
                                                                                          ! field_name: the name of the associated field
                                                                                          ! field_type: the type of the associated field
      integer, intent(out) :: state_id ! the state number in which the associated control field exists
      logical, intent(out) :: active ! false if the control is active in this timestep (e.g. InitialConditions for timesteps >0)
      logical, intent(out), optional :: have_ub, have_lb
      integer, intent(out), optional :: lb_state_id, ub_state_id
      character(len=OPTION_PATH_LEN), intent(out), optional :: lb_material_phase_name, ub_material_phase_name
      character(len=OPTION_PATH_LEN), intent(out), optional :: ub_field_name, lb_field_name

      integer :: s_idx
      character(len=OPTION_PATH_LEN) :: control_type, material_phase_name, name 
      type(scalar_field), pointer :: sfield => null(), sfield2 => null()
      type(vector_field), pointer :: vfield => null(), vfield2 => null()
      type(tensor_field), pointer :: tfield => null(), tfield2 => null()

      call get_option("/adjoint/controls/control[" // int2str(control_no) //"]/name", control_name)
      call get_option("/adjoint/controls/control[" // int2str(control_no) //"]/type/name", control_type)
      call get_option("/adjoint/controls/control[" // int2str(control_no) //"]/type::" // trim(control_type) // "/field_name", name)

      s_idx = scan(trim(name), "::")
      if (s_idx == 0) then 
        FLExit("The control " // trim(control_name) // " uses an invalid field_name. It should be of the form Materialphase::Field")
      end if
      material_phase_name = name(1:s_idx - 1)
      field_name = name(s_idx + 2:len_trim(name))

      active = .true.

      ! Find state associated with the material phase
      do state_id = 1, size(states) 
        if (trim(states(state_id)%name) == trim(material_phase_name)) then
          exit
        end if
      end do
      if (.not. trim(states(state_id)%name) == trim(material_phase_name)) then
        FLExit("Could not find state " // trim(material_phase_name) // " as specified in control " // trim(control_name) // ".")
      end if
      
      select case (trim(control_type))
        case ("initial_condition")
          if (timestep > 0) then
            active = .false.
            field_type = ''
            field_name = ''
            return
          endif
          assert(timestep == 0)
          if (has_scalar_field(states(state_id), field_name)) then
            field_type = "scalar"
          elseif (has_vector_field(states(state_id), field_name)) then
            field_type = "vector"
          elseif (has_tensor_field(states(state_id), field_name)) then
            field_type = "tensor"
          else
            ewrite(0, *) "The control field " // trim(field_name) // " specified in control " // trim(control_name) // " is not a field in the state."
            ewrite(0, *) "The current state is: "
            call print_state(states(state_id))
            FLExit("Check your control's field settings.")
          end if

        case ("source_term")
          if (timestep == 0) then
            active = .false.
            field_type = ''
            field_name = ''
            return
          end if

          if (has_scalar_field(states(state_id), field_name)) then
            field_type = "scalar"
          elseif (has_vector_field(states(state_id), field_name)) then
            field_type = "vector"
          elseif (has_tensor_field(states(state_id), field_name)) then
            field_type = "tensor"
          else
            ewrite(0, *) "The control field " // trim(name) // " specified in control " // trim(control_name) // & 
                       & " is not a field in state " // trim(material_phase_name) // "."
            ewrite(0, *) "The current state is: "
            call print_state(states(state_id))
            FLExit("Check your control's field settings.")
          end if

        case ("boundary_condition")
          FLAbort("Boundary condition control not implemented yet.")
      end select

      !!!! Check bounds settings !!!!
      ! We only execute this when the optional files are present
      if (.not. present(have_lb) .or. .not. present(have_ub)) then
        return
      endif

      have_ub = have_option("/adjoint/controls/control[0]/bounds/upper_bound")
      have_ub = have_option("/adjoint/controls/control[" // int2str(control_no) //"]/bounds/upper_bound")
      have_lb = have_option("/adjoint/controls/control[" // int2str(control_no) //"]/bounds/lower_bound")

      ! Find material phase for the lower bound field
      if (have_lb) then 
        call get_option("/adjoint/controls/control[" // int2str(control_no) //"]/bounds/lower_bound/field_name", name)
        if (have_lb) then
          s_idx = scan(trim(name), "::")
          if (s_idx == 0) then 
            FLExit("The control lower bound for control " // trim(control_name) // " uses an invalid field_name. It should be of the form Materialphase::Field")
          end if
          lb_material_phase_name = name(1:s_idx - 1)
          lb_field_name = name(s_idx + 2:len_trim(name))
          ! Find state associated with the lower bound material phase
          do lb_state_id = 1, size(states) 
            if (trim(states(lb_state_id)%name) == trim(lb_material_phase_name)) then
              exit
            end if
          end do
          if (.not. trim(states(lb_state_id)%name) == trim(lb_material_phase_name)) then
            FLExit("Could not find state " // trim(lb_material_phase_name) // " as specified in the lower bound of control " // trim(control_name) // ".")
          end if
        end if
      end if

      ! Find material phase for the upper bound field
      if (have_ub) then
        call get_option("/adjoint/controls/control[" // int2str(control_no) //"]/bounds/upper_bound/field_name", name)
        if (have_ub) then
          s_idx = scan(trim(name), "::")
          if (s_idx == 0) then 
            FLExit("The control upper bound for control " // trim(control_name) // " uses an invalid field_name. It should be of the form Materialphase::Field")
          end if
          ub_material_phase_name = name(1:s_idx - 1)
          ub_field_name = name(s_idx + 2:len_trim(name))
          ! Find state associated with the upper bound material phase
          do ub_state_id = 1, size(states) 
            if (trim(states(ub_state_id)%name) == trim(ub_material_phase_name)) then
              exit
            end if
          end do
          if (.not. trim(states(ub_state_id)%name) == trim(ub_material_phase_name)) then
            FLExit("Could not find state " // trim(ub_material_phase_name) // " as specified in the upper bound of control " // trim(control_name) // ".")
          end if
        end if
      end if

      ! Check that the lower/upper bound fields indeed exist and have the same mesh than the control
      select case (trim(field_type))

        case ("scalar")
          ! Upper bound
          ! Check that the specified field exist in the state
          if (have_ub) then
            if (.not. has_scalar_field(states(ub_state_id), ub_field_name)) then
              ewrite(0, *) "Can not find the specified field ", trim(ub_field_name), " for the upper bound of control ", trim(control_name)
              FLExit("Check your control's field settings.")
            end if

            ! Check that the control and the upper bound use the same mesh
            sfield => extract_scalar_field(states(state_id), field_name)
            sfield2 => extract_scalar_field(states(ub_state_id), ub_field_name)
            if (.not. sfield%mesh == sfield2%mesh) then
              ewrite(0, *) "Error in control ", trim(control_name), " definition. The upper bound and the control's mesh have to be the same."
              FLExit("Check your control's field settings.")
            end if
          end if

          ! Lower bound
          if (have_lb) then
            if (.not. has_scalar_field(states(state_id), lb_field_name)) then
              ewrite(0, *) "Can not find the specified field ", trim(lb_field_name), " for the lower bound of control ", trim(control_name)
              FLExit("Check your control's field settings.")
            end if

            ! Check that the control and the lower bound use the same mesh
            sfield => extract_scalar_field(states(state_id), field_name)
            sfield2 => extract_scalar_field(states(lb_state_id), lb_field_name)
            if (.not. sfield%mesh == sfield2%mesh) then
              ewrite(0, *) "Error in control ", trim(control_name), " definition. The lower bound and the control's mesh have to be the same."
              FLExit("Check your control's field settings.")
            end if
          end if

        case ("vector")
          ! Upper bound
          ! Check that the specified field exist in the state
          if (have_ub) then
            if (.not. has_vector_field(states(ub_state_id), ub_field_name)) then
              ewrite(0, *) "Can not find the specified field ", trim(ub_field_name), " for the upper bound of control ", trim(control_name)
              FLExit("Check your control's field settings.")
            end if

            ! Check that the control and the upper bound use the same mesh
            vfield => extract_vector_field(states(state_id), field_name)
            vfield2 => extract_vector_field(states(ub_state_id), ub_field_name)
            if (.not. vfield%mesh == vfield2%mesh) then
              ewrite(0, *) "Error in control ", trim(control_name), " definition. The upper bound and the control's mesh have to be the same."
              FLExit("Check your control's field settings.")
            end if
          end if

          ! Lower bound
          if (have_lb) then
            if (.not. has_vector_field(states(state_id), lb_field_name)) then
              ewrite(0, *) "Can not find the specified field ", trim(lb_field_name), " for the lower bound of control ", trim(control_name)
              FLExit("Check your control's field settings.")
            end if

            ! Check that the control and the lower bound use the same mesh
            vfield => extract_vector_field(states(state_id), field_name)
            vfield2 => extract_vector_field(states(lb_state_id), lb_field_name)
            if (.not. vfield%mesh == vfield2%mesh) then
              ewrite(0, *) "Error in control ", trim(control_name), " definition. The lower bound and the control's mesh have to be the same."
              FLExit("Check your control's field settings.")
            end if
          end if

        case ("tensor")
          ! Upper bound
          ! Check that the specified field exist in the state
          if (have_ub) then
            if (.not. has_tensor_field(states(ub_state_id), ub_field_name)) then
              ewrite(0, *) "Can not find the specified field ", trim(ub_field_name), " for the upper bound of control ", trim(control_name)
              FLExit("Check your control's field settings.")
            end if

            ! Check that the control and the upper bound use the same mesh
            tfield => extract_tensor_field(states(state_id), field_name)
            tfield2 => extract_tensor_field(states(ub_state_id), ub_field_name)
            if (.not. tfield%mesh == tfield2%mesh) then
              ewrite(0, *) "Error in control ", trim(control_name), " definition. The lower bound and the control's mesh have to be the same."
              FLExit("Check your control's field settings.")
            end if
          end if

          ! Lower bound
          if (have_lb) then
            if (.not. has_tensor_field(states(state_id), lb_field_name)) then
              ewrite(0, *) "Can not find the specified field ", trim(lb_field_name), " for the lower bound of control ", trim(control_name)
              FLExit("Check your control's field settings.")
            end if

            ! Check that the control and the lower bound use the same mesh
            tfield => extract_tensor_field(states(state_id), field_name)
            tfield2 => extract_tensor_field(states(lb_state_id), lb_field_name)
            if (.not. sfield%mesh == sfield2%mesh) then
              ewrite(0, *) "Error in control ", trim(control_name), " definition. The lower bound and the control's mesh have to be the same."
              FLExit("Check your control's field settings.")
            end if
          end if

      end select
    end subroutine get_control_details

    ! Loop through all controls in the option file and store them in files.
    subroutine adjoint_write_controls(timestep, dt, states)
      integer, intent(in) :: timestep
      real, intent(in) :: dt
      type(state_type), dimension(:), intent(in) :: states
#ifdef HAVE_ADJOINT
      integer :: nb_controls, i, state_id
      logical :: active
      character(len=OPTION_PATH_LEN) :: simulation_name, field_name, field_type, control_name 
      logical :: have_ub, have_lb
      integer :: lb_state_id, ub_state_id
      character(len=OPTION_PATH_LEN) :: lb_material_phase_name, ub_material_phase_name
      character(len=OPTION_PATH_LEN) :: ub_field_name, lb_field_name

      if (.not. have_option("/adjoint/controls")) then
        return
      end if

      nb_controls = option_count("/adjoint/controls/control")
      call get_option("/simulation_name", simulation_name)
      ! Now loop over the controls and save them to disk
      do i = 0, nb_controls-1
        call get_control_details(states, timestep, i, control_name, state_id, field_name, field_type, active, have_lb, lb_state_id, lb_material_phase_name, lb_field_name, have_ub, ub_state_id, ub_material_phase_name, ub_field_name)
        if (active) then

          ! Save the control parameter to disk
          call python_reset()
          call python_add_state(states(state_id))
          call python_run_string("import pickle;" // &
                              &  "fname = 'control_" // trim(simulation_name) // "_" // trim(control_name) // "_" // int2str(timestep) // ".pkl';" // & 
                              &  "f = open(fname, 'wb');" // &
                              &  "field = state." // trim(field_type) // "_fields['" // trim(field_name) // "'];" // &
                              &  "pickle.dump(field.val[:], f);" // &
                              &  "f.close();")

          ! And its bounds
          call python_reset()
          if (have_lb) then
            call python_add_state(states(lb_state_id))
            call python_run_string("import pickle;" // &
                                &  "fname = 'control_" // trim(simulation_name) // "_" // trim(control_name) // "_lower_bound_" // int2str(timestep) // ".pkl';" // & 
                                &  "f = open(fname, 'wb');" // &
                                &  "field = state." // trim(field_type) // "_fields['" // trim(lb_field_name) // "'];" // &
                                &  "pickle.dump(field.val[:], f);" // &
                                &  "f.close();")
            call python_reset()
          end if

          if (have_ub) then
            call python_add_state(states(ub_state_id))
            call python_run_string("import pickle;" // &
                                &  "fname = 'control_" // trim(simulation_name) // "_" // trim(control_name) // "_upper_bound_" // int2str(timestep) // ".pkl';" // & 
                                &  "f = open(fname, 'wb');" // &
                                &  "field = state." // trim(field_type) // "_fields['" // trim(ub_field_name) // "'];" // &
                                &  "pickle.dump(field.val[:], f);" // &
                                &  "f.close();")
            call python_reset()
          endif

        end if
      end do
#endif
    end subroutine adjoint_write_controls

    ! The counterpart of adjoint_write_controls: Loop through all controls in the option file and fill the fields from the control files.
    subroutine adjoint_load_controls(timestep, dt, states)
      integer, intent(in) :: timestep
      real, intent(in) :: dt
      type(state_type), dimension(:), intent(in) :: states
#ifdef HAVE_ADJOINT
      integer :: nb_controls, i, state_id, s_idx
      logical :: active, file_exists
      character(len=OPTION_PATH_LEN) :: field_name, field_type, control_name, simulation_name, control_type

      ! Do not read anything if we are supposed to write the controls
      if (.not. have_option("/adjoint/controls/load_controls")) then
        return 
      end if

      call get_option("/simulation_name", simulation_name)
      nb_controls = option_count("/adjoint/controls/control")

      ! Now loop over the controls, read their controls files and fill the associated fields
      do i = 0, nb_controls-1
        call get_control_details(states, timestep, i, control_name, state_id, field_name, field_type, active)
        call get_option("/adjoint/controls/control[" // int2str(i) // "]/type/name", control_type)

        ! We load all OTHER controls, then apply our strong boundary conditions, then load BC controls
        if (trim(control_type) == "boundary_condition") cycle

        if (active) then
          ! Check that the control parameter file extists
          inquire(file="control_" // trim(simulation_name) // "_" // trim(control_name) // "_" // int2str(timestep) // ".pkl", exist=file_exists)
          if (.not. file_exists) then
            ewrite (0, *) "File 'control_" // trim(simulation_name) // "_" // trim(control_name) // &
                        & "_" // int2str(timestep) // ".pkl' does not exist."
            FLExit("Error while loading control files. You might want to remove '/adjoint/controls/load_controls'?")
          end if
          ! Read the control parameter from disk
          call python_reset()
          call python_add_state(states(state_id))
          call python_run_string("import pickle;" // &
                              &  "fname = 'control_" // trim(simulation_name) // "_" // trim(control_name) // "_" // int2str(timestep) // ".pkl';" // & 
                              &  "f = open(fname, 'rb');" // &
                              &  "field = state." // trim(field_type) // "_fields['" // trim(field_name) // "'];" // &
                              &  "field.val[:] = pickle.load(f);" // &
                              &  "f.close();")
          call python_reset()
        end if
      end do

      ! Now apply the BCs, because the BCs should overwrite (say) an initial_condition control
      call set_dirichlet_consistent(states)

      do i = 0, nb_controls-1
        call get_control_details(states, timestep, i, control_name, state_id, field_name, field_type, active)
        call get_option("/adjoint/controls/control[" // int2str(i) // "]/type/name", control_type)

        ! We load all OTHER controls, then apply our strong boundary conditions, then load BC controls
        if (trim(control_type) /= "boundary_condition") cycle

        if (active) then
          ! Check that the control parameter file extists
          inquire(file="control_" // trim(simulation_name) // "_" // trim(control_name) // "_" // int2str(timestep) // ".pkl", exist=file_exists)
          if (.not. file_exists) then
            ewrite (0, *) "File 'control_" // trim(simulation_name) // "_" // trim(control_name) // &
                        & "_" // int2str(timestep) // ".pkl' does not exist."
            FLExit("Error while loading control files. You might want to remove '/adjoint/controls/load_controls'?")
          end if
          ! Read the control parameter from disk
          call python_reset()
          call python_add_state(states(state_id))
          call python_run_string("import pickle;" // &
                              &  "fname = 'control_" // trim(simulation_name) // "_" // trim(control_name) // "_" // int2str(timestep) // ".pkl';" // & 
                              &  "f = open(fname, 'rb');" // &
                              &  "field = state." // trim(field_type) // "_fields['" // trim(field_name) // "'];" // &
                              &  "field.val[:] = pickle.load(f);" // &
                              &  "f.close();")
          call python_reset()
        end if
      end do
#endif
    end subroutine adjoint_load_controls

    ! Extracts the total derivatives of the functional and saves them to disk.
    subroutine adjoint_write_control_derivatives(states, timestep, functional_name)
      type(state_type), dimension(:), intent(inout) :: states
      integer, intent(in) :: timestep
      character(len=*), intent(in) :: functional_name
      character(len=OPTION_PATH_LEN) :: field_name, control_type, control_deriv_name, field_deriv_name, name, material_phase_name, simulation_name, field_type
      integer :: nb_controls
      integer :: i, state_id, s_idx
      type(scalar_field), pointer :: sfield => null()
      type(vector_field), pointer :: vfield => null()
      type(tensor_field), pointer :: tfield => null()
     
      if (.not. have_option("/adjoint/controls")) then
        return
      end if
      call get_option("/simulation_name", simulation_name)
      ! Remove the suffix _adjoint from the simulation_name

      nb_controls = option_count("/adjoint/controls/control")
      do i = 0, nb_controls-1
        call get_option("/adjoint/controls/control[" // int2str(i) //"]/name", control_deriv_name)
        call get_option("/adjoint/controls/control[" // int2str(i) //"]/type/name", control_type)
        call get_option("/adjoint/controls/control[" // int2str(i) //"]/type::" // trim(control_type) // "/field_name", name)
        s_idx = scan(trim(name), "::")
        if (s_idx == 0) then 
          FLExit("The control " // trim(control_deriv_name) // " uses an invalid field_name. It should be of the form Materialphase::Fieldname")
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
          FLExit("Could not find state " // trim(material_phase_name) // " as specified in control " // trim(control_deriv_name) // ".")
        end if

        field_name = "Adjoint" // trim(field_name)
        control_deriv_name = trim(control_deriv_name) // "_TotalDerivative"
        if (trim(control_type) == "initial_condition" .or. trim(control_type) == "source_term") then

          if (trim(control_type) == "initial_condition" .and. timestep > 0) then
            cycle
          end if

          if (trim(control_type) == "source_term" .and. timestep == 0) then
            cycle
          end if

          field_deriv_name = trim(functional_name) // "_" // control_deriv_name 
          if (has_scalar_field(states(state_id), field_name)) then
            field_type = "scalar"
          elseif (has_vector_field(states(state_id), field_name)) then
            field_type = "vector"
          elseif (has_tensor_field(states(state_id), field_name)) then
            field_type = "tensor"
          else
            ewrite(0, *) "The control field " // trim(field_deriv_name) // " specified in control " // trim(control_deriv_name) // " is not a field in the state."
            ewrite(0, *) "The current state is: "
            call print_state(states(state_id))
            FLExit("Check your control's field settings.")
          end if
        elseif (trim(control_type) == "boundary_condition") then
          FLAbort("Boundary condition control not implemented yet.")
        end if
        ! Save the control parameter to disk
        call python_reset()
        call python_add_state(states(state_id))
        call python_run_string("import pickle;" // &
                            &  "fname = 'control_" // trim(simulation_name) // "_" // trim(control_deriv_name) // "_" // int2str(timestep) // ".pkl';" // & 
                            &  "f = open(fname, 'wb');" // &
                            &  "field = state." // trim(field_type) // "_fields['" // trim(field_deriv_name) // "'];" // &
                            &  "pickle.dump(field.val[:], f);" // &
                            &  "f.close();")
        call python_reset()
      end do
    end subroutine adjoint_write_control_derivatives

    subroutine allocate_and_insert_control_derivative_fields(states)
      type(state_type), dimension(:), intent(inout) :: states
      integer :: nb_controls, nb_functionals, i, state_id, functional, s_idx
      character(len=OPTION_PATH_LEN) :: field_name, control_type, control_deriv_name, functional_name, field_deriv_name, material_phase_name, name
      type(scalar_field), pointer :: sfield => null()
      type(vector_field), pointer :: vfield => null()
      type(tensor_field), pointer :: tfield => null()
      type(scalar_field) :: sfield_deriv
      type(vector_field) :: vfield_deriv
      type(tensor_field) :: tfield_deriv

      nb_controls = option_count("/adjoint/controls/control")
      nb_functionals = option_count("/adjoint/functional")
      ! Loop over the functionals
      do functional = 0, nb_functionals-1
        call get_option("/adjoint/functional[" // int2str(functional) // "]/name", functional_name)
        ! Now loop over the controls 
        do i = 0, nb_controls-1
          call get_option("/adjoint/controls/control[" // int2str(i) //"]/name", control_deriv_name)
          call get_option("/adjoint/controls/control[" // int2str(i) //"]/type/name", control_type)
          call get_option("/adjoint/controls/control[" // int2str(i) //"]/type::" // trim(control_type) // "/field_name", name)
          s_idx = scan(trim(name), "::")
          if (s_idx == 0) then 
            FLExit("The control " // trim(control_deriv_name) // " uses an invalid field_name. It should be of the form Materialphase::Fieldname")
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
            FLExit("Could not find state " // trim(material_phase_name) // " as specified in control " // trim(control_deriv_name) // ".")
          end if

          field_name = "Adjoint" // trim(field_name)
          control_deriv_name = trim(control_deriv_name) // "_TotalDerivative"

          if (trim(control_type) == "initial_condition" .or. trim(control_type) == "source_term") then
            field_deriv_name = trim(functional_name) // "_" // control_deriv_name 

            if (has_scalar_field(states(state_id), field_name)) then
              sfield => extract_scalar_field(states(state_id), field_name)
              call allocate(sfield_deriv, sfield%mesh, field_deriv_name)
              call zero(sfield_deriv)
              call insert(states(state_id), sfield_deriv, field_deriv_name)
              call deallocate(sfield_deriv)

            elseif (has_vector_field(states(state_id), field_name)) then
              vfield => extract_vector_field(states(state_id), field_name)
              call allocate(vfield_deriv, vfield%dim, vfield%mesh, field_deriv_name)
              call zero(vfield_deriv)
              call insert(states(state_id), vfield_deriv, field_deriv_name)
              call deallocate(vfield_deriv)

            elseif (has_tensor_field(states(state_id), field_name)) then
              tfield => extract_tensor_field(states(state_id), field_name)
              call allocate(tfield_deriv, tfield%mesh, field_deriv_name, tfield%field_type, tfield%dim)
              call zero(tfield_deriv)
              call insert(states(state_id), tfield_deriv, field_deriv_name)
              call deallocate(tfield_deriv)

            else
              ewrite(0, *) "The control field " // trim(field_name) // " specified in control " // trim(control_deriv_name) // " is not a field in the state."
              ewrite(0, *) "The current state is: "
              call print_state(states(state_id))
              FLExit("Check your control's field settings.")
            end if

          elseif (trim(control_type) == "boundary_condition") then
            FLAbort("Boundary condition control not implemented yet.")
          end if

        end do
      end do
    end subroutine allocate_and_insert_control_derivative_fields

end module adjoint_controls
