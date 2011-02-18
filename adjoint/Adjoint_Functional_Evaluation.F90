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
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
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

module adjoint_functional_evaluation
#ifdef HAVE_ADJOINT
#include "libadjoint/adj_fortran.h"
  use libadjoint
  use iso_c_binding
  use fields
  use fldebug
  use global_parameters, only : PYTHON_FUNC_LEN, OPTION_PATH_LEN
  use spud
  use state_module
  use embed_python
  use python_state
  use libadjoint_data_callbacks

  public :: libadjoint_functional_derivative

  contains

  subroutine libadjoint_functional_derivative(var, ndepends, dependencies, values, name_c, start_time, end_time, output) bind(c)
    type(adj_variable), intent(in) :: var
    integer(kind=c_int), intent(in), value :: ndepends
    type(adj_variable), dimension(ndepends), intent(in) :: dependencies
    type(adj_vector), dimension(ndepends), intent(in) :: values
    character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: name_c
    adj_scalar_f, intent(in), value :: start_time
    adj_scalar_f, intent(in), value :: end_time
    type(adj_vector), intent(out) :: output

    character(len=PYTHON_FUNC_LEN) :: code
    integer :: i
    character(len=ADJ_NAME_LEN) :: name_f
    character(len = 30) :: buffer
    integer :: timestep, ierr, max_timestep, tmp_timestep, min_timestep, nmaterial_phases
    type(state_type), dimension(:,:), pointer :: states ! material_phases x timesteps
    character(len=OPTION_PATH_LEN) :: material_phase_name

    do i=1,size(name_c)
      name_f(i:i) = name_c(i)
    end do
    ierr = adj_variable_get_timestep(var, timestep)

    ! Fetch the python code the user has helpfully supplied.
    if (.not. have_option("/adjoint/functional::" // trim(name_f) // "/functional_derivative/algorithm")) then
      FLAbort("We really should have /adjoint/functional::" // trim(name_f) // "/functional_derivative/algorithm")
    endif
    call get_option("/adjoint/functional::" // trim(name_f) // "/functional_derivative/algorithm", code)

    ! Convert from the list of adj_values/adj_vectors to actual states
    if (ndepends > 0) then
      ierr = adj_variable_get_timestep(dependencies(1), max_timestep)
      min_timestep = max_timestep

      do i=2,ndepends
        ierr = adj_variable_get_timestep(dependencies(i), tmp_timestep)
        max_timestep = max(max_timestep, tmp_timestep)
        min_timestep = min(min_timestep, tmp_timestep)
      end do
      nmaterial_phases = option_count("/material_phase")
      allocate(states(nmaterial_phases, min_timestep:max_timestep))

      ! Set the names of the material phases, so that we can look them up and know
      ! which material_phase to put stuff into in convert_adj_vectors_to_states
      do i=0,nmaterial_phases-1
        call get_option("/material_phase[" // int2str(i) // "]/name", material_phase_name)
        states(i+1,:)%name = trim(material_phase_name)
      end do
      call convert_adj_vectors_to_states(dependencies, values, states)
      call python_add_states_time(states)
    endif

    ! Also set up some useful variables for the user to use
    call python_reset
    write(buffer,*) start_time
    call python_run_string("time = " // trim(buffer))

    write(buffer, *) end_time - start_time
    call python_run_string("dt = " // trim(buffer))

    write(buffer, *) timestep
    call python_run_string("n = " // trim(buffer))

    if (max_timestep >= 0) then
      call deallocate(states)
      deallocate(states)
    endif

  end subroutine libadjoint_functional_derivative

  subroutine convert_adj_vectors_to_states(dependencies, values, states)
    type(adj_variable), dimension(:), intent(in) :: dependencies
    type(adj_vector), dimension(:), intent(in) :: values
    type(state_type), dimension(:,:), intent(inout), pointer :: states ! material_phases x time

    integer :: j
    character(len=ADJ_NAME_LEN) :: name, material_phase_name, field_name
    integer :: timestep
    integer :: ierr
    type(scalar_field) :: sfield
    type(vector_field) :: vfield
    type(tensor_field) :: tfield
    type(mesh_type)    :: mesh
    integer :: s_idx
    integer :: matpas ! to loop over material_phases

    do j=1,size(dependencies)
      ierr = adj_variable_get_name(dependencies(j), name)
      ierr = adj_variable_get_timestep(dependencies(j), timestep)

      ! Do some fortran string parsing, woohoo
      s_idx = scan(trim(name), ":")
      material_phase_name = name(1:s_idx - 1)
      field_name = name(s_idx + 2:len_trim(name))

      select case(values(j)%klass)
      case (ADJ_SCALAR_FIELD)
        call field_from_adj_vector(values(j), sfield)
        do matpas=1,size(states, 1)
          if (trim(states(matpas,timestep)%name) == trim(material_phase_name)) then
            call insert(states(matpas, timestep), sfield, trim(field_name))
          endif
        end do
      case (ADJ_VECTOR_FIELD)
        call field_from_adj_vector(values(j), vfield)
        do matpas=1,size(states, 1)
          if (trim(states(matpas,timestep)%name) == trim(material_phase_name)) then
            call insert(states(matpas, timestep), vfield, trim(field_name))
          endif
        end do
      case (ADJ_TENSOR_FIELD)
        call field_from_adj_vector(values(j), tfield)
        do matpas=1,size(states, 1)
          if (trim(states(matpas,timestep)%name) == trim(material_phase_name)) then
            call insert(states(matpas, timestep), tfield, trim(field_name))
          endif
        end do
      case (ADJ_MESH_TYPE)
        call mesh_type_from_adj_vector(values(j), mesh)
        do matpas=1,size(states, 1)
          call insert(states(matpas, :), mesh, trim(field_name))
        end do
      case default
        FLAbort("Unknown adj_vector%klass")
      end select
    end do
  end subroutine convert_adj_vectors_to_states
#endif
end module adjoint_functional_evaluation
