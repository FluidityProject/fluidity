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
!    C.Pain@Imperial.ac.uk
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

module discrete_properties_module
  use state_module
  use solenoidal_interpolation_module
  use fields
  use spud
  use global_parameters, only : OPTION_PATH_LEN
  use field_options
  
  implicit none
  
  private
  
  public :: enforce_discrete_properties

contains

  subroutine enforce_discrete_properties(states, only_prescribed, exclude_interpolated, exclude_nonreprescribed)
    type(state_type), dimension(:), intent(inout) :: states
    ! if a field isn't prescribed then don't process it
    logical, intent(in), optional :: only_prescribed
    ! if a field has interpolation options then don't process it
    logical, intent(in), optional :: exclude_interpolated
    ! if a field hasn't been represcribed then don't process it
    logical, intent(in), optional :: exclude_nonreprescribed

    ! The fields organised by algorithm.
    type(state_type) :: alg_state
    
    integer :: state, state_cnt, mesh_i
    integer :: mesh_cnt

    character(len=255), dimension(1), parameter :: algorithms = (/&
                       & "solenoidal" /)
    integer :: alg_cnt, alg

    type(mesh_type), pointer :: mesh
    type(vector_field), pointer :: pos

    integer :: processed

    ewrite(1, *) "In enforce_discrete_properties"

    mesh_cnt = option_count("/geometry/mesh")
    alg_cnt = size(algorithms)
    state_cnt = size(states)

    do state = 1, state_cnt
      ! some fields require post processing that needs more than one mesh
      ! these are dealt with in this loop
      ewrite(2, *) "Processing fields in state " // trim(states(state)%name)

      processed = 0
      alg_loop: do alg = 1, alg_cnt
        ! only actually care about one algorithm at the moment but for futureproofing we'll do this in a loop
        ewrite(2, *) "  Considering algorithm " // trim(algorithms(alg))
        
        do mesh_i = 1, mesh_cnt
          mesh => extract_mesh(states(state), mesh_i)
          call insert(alg_state, mesh, name=trim(mesh%name))
        end do
        pos => extract_vector_field(states(state), "Coordinate")
        call insert(alg_state, pos, "Coordinate")

        select case(trim(algorithms(alg)))
          case("solenoidal")
            call collect_fields_to_process(interpolate_field_solenoidal, states(state), alg_state, &
                                            only_prescribed=only_prescribed, &
                                            exclude_interpolated=exclude_interpolated, &
                                            exclude_nonreprescribed=exclude_nonreprescribed)
          
            if(field_count(alg_state) > 1) then
              call solenoidal_interpolation(alg_state)
            end if
          case default
            FLAbort("Unknown discrete property algorithm.")
        end select

        processed = processed + field_count(alg_state) - 1
        assert(processed <= field_count(states(state)))

        call deallocate(alg_state)

        if(processed >= field_count(states(state))) then
          ewrite(2, *) "All fields in state " // trim(states(state)%name) // " processed"
          exit alg_loop
        end if
        
      end do alg_loop

    end do

    ewrite(1, *) "Exiting enforce_discrete_properties"
    
  end subroutine enforce_discrete_properties

  function interpolate_field_solenoidal(option_path, only_prescribed, exclude_interpolated, exclude_nonreprescribed) result(process)
    character(len = *), intent(in) :: option_path
    logical, intent(in), optional :: only_prescribed
    logical, intent(in), optional :: exclude_interpolated
    logical, intent(in), optional :: exclude_nonreprescribed

    logical :: process
    
    character(len = OPTION_PATH_LEN) :: base_path
    
    process = .false.
    if(len_trim(option_path) == 0) return

    base_path = trim(complete_field_path(option_path))
    
    process = have_option(trim(base_path) // "/enforce_discrete_properties/solenoidal") &
                  .or. have_option(trim(base_path) // "/enforce_discrete_properties/solenoidal_lagrange_update")
                  
    if(present_and_true(exclude_nonreprescribed)) then
      process = process .and. .not.have_option(trim(base_path)//"/do_not_recalculate")
    end if
    
    if(present_and_true(exclude_interpolated)) then
      process = process .and. .not.interpolate_field(trim(option_path))
    end if
    
    if(present_and_true(only_prescribed)) then
      process = process .and. have_option(trim(option_path)//"/prescribed")
    end if

  end function interpolate_field_solenoidal
  
  subroutine collect_fields_to_process(test, input_state, output_state, &
                                       only_prescribed, exclude_interpolated, exclude_nonreprescribed)
    !!< Collect all fields in the supplied input states that
    !!< pass the supplied test, and insert them into output_state.
    
    interface
      function test(option_path, only_prescribed, exclude_interpolated, exclude_nonreprescribed)
        implicit none
        character(len = *), intent(in) :: option_path
        logical, intent(in), optional :: only_prescribed
        logical, intent(in), optional :: exclude_interpolated
        logical, intent(in), optional :: exclude_nonreprescribed
        logical :: test
      end function test
    end interface
    type(state_type), intent(in) :: input_state
    type(state_type), intent(inout) :: output_state
    logical, intent(in), optional :: only_prescribed
    logical, intent(in), optional :: exclude_interpolated
    logical, intent(in), optional :: exclude_nonreprescribed

    integer :: i
    type(scalar_field), pointer :: s_field => null()
    type(tensor_field), pointer :: t_field => null()
    type(vector_field), pointer :: v_field => null()
    
    do i = 1, scalar_field_count(input_state)
      s_field => extract_scalar_field(input_state, i)
      if(.not.aliased(s_field).and.test(s_field%option_path, &
                                        only_prescribed=only_prescribed, &
                                        exclude_interpolated=exclude_interpolated, &
                                        exclude_nonreprescribed=exclude_nonreprescribed)) then
      
        ewrite(2, *) "    Found ", trim(s_field%name)
          
        call insert(output_state, s_field, trim(input_state%scalar_names(i)))
      end if
    end do
    
    do i = 1, vector_field_count(input_state)
      v_field => extract_vector_field(input_state, i)
      if(.not.aliased(v_field).and.test(v_field%option_path, &
                                        only_prescribed=only_prescribed, &
                                        exclude_interpolated=exclude_interpolated, &
                                        exclude_nonreprescribed=exclude_nonreprescribed)) then
        ewrite(2, *) "    Found ", trim(v_field%name)
          
        call insert(output_state, v_field, trim(input_state%vector_names(i)))
      end if
    end do
    
    do i = 1, tensor_field_count(input_state)
      t_field => extract_tensor_field(input_state, i)
      if(.not.aliased(t_field).and.test(t_field%option_path, &
                                        only_prescribed=only_prescribed, &
                                        exclude_interpolated=exclude_interpolated, &
                                        exclude_nonreprescribed=exclude_nonreprescribed)) then
        ewrite(2, *) "    Found ", trim(t_field%name)
          
        call insert(output_state, t_field, trim(input_state%tensor_names(i)))
      end if
    end do
    
  end subroutine collect_fields_to_process

end module discrete_properties_module
