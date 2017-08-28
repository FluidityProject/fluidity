!    Copyright (C) 2006 Imperial College London and others.   
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
#include "version.h"

module diagnostic_variables
  !!< A module to calculate and output diagnostics. This replaces the .s file.
  use iso_c_binding, only: c_long
  use fldebug 
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, &
& PYTHON_FUNC_LEN, integer_size, real_size
  use quadrature
  use futils
  use elements
  use spud
  use mpi_interfaces
  use parallel_tools
  use memory_diagnostics
  use integer_hash_table_module
  use data_structures
  use linked_lists
  use halo_data_types
  use halos_base
  use halos_debug
  use halos_allocates
  use ieee_arithmetic
  use sparse_tools
  use embed_python
  use fields_base
  use eventcounter
  use fetools
  use unittest_tools
  use halos_communications
  use halos_numbering
  use halos_ownership
  use parallel_fields, only: element_owned
  use fields
  use profiler
  use state_module
  use vtk_interfaces
  use halos_derivation
  use halos_registration
  use field_derivatives
  use field_options
  use c_interfaces
  use fefields
  use meshdiagnostics
  use sparsity_patterns
  use solvers
  use write_state_module, only: vtk_write_state_new_options
  use surface_integrals
  use detector_data_types
  use pickers
  use mixing_statistics
  use detector_tools
  use detector_parallel
  use detector_move_lagrangian
  use state_fields_module
  
  implicit none

  interface
     subroutine register_diagnostics()
     end subroutine register_diagnostics
  end interface

  private

  public :: initialise_diagnostics, initialise_convergence, &
       & initialise_steady_state, field_tag, write_diagnostics, &
       & test_and_write_convergence, initialise_detectors, write_detectors, &
       & test_and_write_steady_state, steady_state_field, convergence_field, &
       & close_diagnostic_files, run_diagnostics, &
       & diagnostic_variables_check_options, list_det_into_csr_sparsity, &
       & initialise_walltime, &
       & uninitialise_diagnostics, register_diagnostic, destroy_registered_diagnostics, set_diagnostic, &
       & get_diagnostic, initialise_constant_diagnostics, create_single_detector

  public :: default_stat
  public :: stat_type

  interface stat_field
    module procedure stat_field_scalar, stat_field_vector, stat_field_tensor
  end interface stat_field

  interface convergence_field
    module procedure convergence_field_scalar, convergence_field_vector
  end interface convergence_field

  interface steady_state_field
    module procedure steady_state_field_scalar, steady_state_field_vector
  end interface steady_state_field

  interface detector_field
     module procedure detector_field_scalar, detector_field_vector
  end interface

  ! List of registered diagnostic
  type registered_diagnostic_item
     integer :: dim
     character(len=FIELD_NAME_LEN) :: name
     character(len=FIELD_NAME_LEN) :: statistic
     character(len=FIELD_NAME_LEN) :: material_phase
     logical :: have_material_phase
     real, dimension(:), allocatable :: value
     type(registered_diagnostic_item), pointer :: next => null()
  end type registered_diagnostic_item

  type :: stat_type
    ! Idempotency variable
    logical :: initialised=.false.
    logical :: detectors_initialised = .false.
    logical :: convergence_initialised = .false.
    logical :: steady_state_initialised = .false.

    !! Output unit for diagnostics file.
    !! (assumed non-opened as long these are 0)
    integer :: diag_unit=0, conv_unit=0, &
      & detector_checkpoint_unit=0, detector_file_unit=0 

    !! Are we writing to a convergence file?
    logical :: write_convergence_file=.false.

    !! Output unit for .steady_state file (assumed non-opened as long as == 0)
    integer :: steady_state_unit = 0
    !! Are we writing to a steady state file?
    logical :: write_steady_state_file = .false.
    logical :: binary_steady_state_output = .false.

   !! Are we continuing from a detector checkpoint file?
    logical :: from_checkpoint = .false.

    !The following variable will switch to true if call to zoltan_drive for the re-load balance occur.
    logical :: zoltan_drive_call = .false.

    character(len = FIELD_NAME_LEN), dimension(:), allocatable :: mesh_list
    !! List of scalar fields to output. This is stored to ensure that
    !! additional fields inserted into state during running do not bugger up
    !! the output.
    type(stringlist), dimension(:), allocatable :: sfield_list
    type(stringlist), dimension(:), allocatable :: vfield_list
    type(stringlist), dimension(:), allocatable :: tfield_list

    ! Names of detector groups used for checkpointing; names are held in read order
    character(len = FIELD_NAME_LEN), dimension(:), allocatable :: detector_group_names
    integer, dimension(:), allocatable :: number_det_in_each_group

    type(detector_linked_list) :: detector_list

    type(registered_diagnostic_item), pointer :: registered_diagnostic_first => NULL()
    
    !! Recording wall time since the system start
    integer :: current_count, count_rate, count_max
    integer(kind = c_long) :: elapsed_count
  end type stat_type

  type(stat_type), save, target :: default_stat

contains

  function stat_mesh(mesh)
    !!< Return whether the supplied mesh should be included in the .stat file

    type(mesh_type), intent(in) :: mesh

    logical :: stat_mesh
    integer :: stat
    character(len = OPTION_PATH_LEN) :: stat_test_path

    stat_mesh = .false.
    stat_test_path=trim(complete_mesh_path(mesh%option_path,stat))
    if(stat==0) then
       stat_mesh = have_option(trim(stat_test_path) // "/stat/include_in_stat") &
            &.and..not.have_option(trim(stat_test_path) // "/stat/exclude_from_stat")
    end if

  end function stat_mesh

  function stat_field_scalar(sfield, state)
    !!< Return whether the supplied field should be included in the .stat file

    type(scalar_field), target, intent(in) :: sfield
    type(state_type), intent(in) :: state

    logical :: stat_field_scalar

    character(len = OPTION_PATH_LEN) :: stat_test_path
    logical :: include_test
    type(scalar_field), pointer :: parent_sfield

    if(sfield%name(:3) == "Old") then
      parent_sfield => extract_scalar_field(state, trim(sfield%name(4:)))
      stat_test_path = "/include_previous_time_step"
      include_test = .true.
    else if(sfield%name(:9) == "Nonlinear") then
      parent_sfield => extract_scalar_field(state, trim(sfield%name(10:)))
      stat_test_path = "/include_nonlinear_field"
      include_test = .true.
    else if(sfield%name(:8) == "Iterated") then
      parent_sfield => extract_scalar_field(state, trim(sfield%name(9:)))
      stat_test_path = "/include_nonlinear_field"
      include_test = .true.
    else
      parent_sfield => sfield
      stat_test_path = "/exclude_from_stat"
      include_test = .false.
    end if

    if((len_trim(parent_sfield%option_path) == 0).or.aliased(parent_sfield)) then
      stat_field_scalar = .false.
      return
    end if

    if(.not. have_option(trim(complete_field_path(parent_sfield%option_path)) // "/stat")) then
      stat_field_scalar = .false.
    else
      stat_field_scalar = have_option(trim(complete_field_path(parent_sfield%option_path)) // "/stat" // trim(stat_test_path))
      if(.not. include_test) then
        stat_field_scalar = (.not. stat_field_scalar)
      end if
    end if

  end function stat_field_scalar

  function stat_field_vector(vfield, state, test_for_components)
    !!< Return whether the supplied field should be included in the .stat file.

    type(vector_field), target, intent(in) :: vfield
    type(state_type), intent(in) :: state
    logical, optional, intent(in) :: test_for_components

    logical :: stat_field_vector

    character(len = OPTION_PATH_LEN) :: stat_test_path
    type(vector_field), pointer :: parent_vfield => null()

    if(vfield%name(:3) == "Old") then
      parent_vfield => extract_vector_field(state, trim(vfield%name(4:)))
      stat_test_path = "/stat/previous_time_step"
    else if(vfield%name(:9) == "Nonlinear") then
      parent_vfield => extract_vector_field(state, trim(vfield%name(10:)))
      stat_test_path = "/stat/nonlinear_field"
    else if(vfield%name(:8) == "Iterated") then
      parent_vfield => extract_vector_field(state, trim(vfield%name(9:)))
      stat_test_path = "/stat/nonlinear_field"
    else
      parent_vfield => vfield
      stat_test_path =  "/stat"
    end if

    if((len_trim(parent_vfield%option_path) == 0).or.aliased(parent_vfield)) then
      stat_field_vector = .false.
    else if(.not. have_option(trim(complete_field_path(parent_vfield%option_path)) // trim(stat_test_path))) then
      stat_field_vector = .false.
    else if(present_and_true(test_for_components)) then
      stat_field_vector = (.not. have_option(trim(complete_field_path(parent_vfield%option_path)) &
        & // trim(stat_test_path) // "/exclude_components_from_stat") &
        & .and. .not. have_option(trim(complete_field_path(parent_vfield%option_path)) &
        & // trim(stat_test_path) // "/exclude_from_stat"))
    else
      stat_field_vector = (.not. have_option(trim(complete_field_path(parent_vfield%option_path)) &
        & // trim(stat_test_path) // "/exclude_from_stat"))
    end if

  end function stat_field_vector

  function stat_field_tensor(tfield, state, test_for_components)
    !!< Return whether the supplied field should be included in the .stat file.

    type(tensor_field), target, intent(in) :: tfield
    type(state_type), intent(in) :: state
    logical, optional, intent(in) :: test_for_components

    logical :: stat_field_tensor

    character(len = OPTION_PATH_LEN) :: stat_test_path
    type(tensor_field), pointer :: parent_tfield => null()

    if(tfield%name(:3) == "Old") then
      parent_tfield => extract_tensor_field(state, trim(tfield%name(4:)))
      stat_test_path = "/stat/previous_time_step"
    else if(tfield%name(:9) == "Nonlinear") then
      parent_tfield => extract_tensor_field(state, trim(tfield%name(10:)))
      stat_test_path = "/stat/nonlinear_field"
    else if(tfield%name(:8) == "Iterated") then
      parent_tfield => extract_tensor_field(state, trim(tfield%name(9:)))
      stat_test_path = "/stat/nonlinear_field"
    else
      parent_tfield => tfield
      stat_test_path =  "/stat"
    end if

    if((len_trim(parent_tfield%option_path) == 0).or.aliased(parent_tfield)) then
      stat_field_tensor = .false.
    else if(.not. have_option(trim(complete_field_path(parent_tfield%option_path)) // trim(stat_test_path))) then
      stat_field_tensor = .false.
    else if(present_and_true(test_for_components)) then
      stat_field_tensor = (.not. have_option(trim(complete_field_path(parent_tfield%option_path)) &
        & // trim(stat_test_path) // "/exclude_components_from_stat") &
        & .and. .not. have_option(trim(complete_field_path(parent_tfield%option_path)) &
        & // trim(stat_test_path) // "/exclude_from_stat"))
    else
      stat_field_tensor = (.not. have_option(trim(complete_field_path(parent_tfield%option_path)) &
        & // trim(stat_test_path) // "/exclude_from_stat"))
    end if

  end function stat_field_tensor

  function convergence_field_scalar(sfield)
    !!< Return whether the supplied field should be included in the .convergence file

    type(scalar_field), target, intent(in) :: sfield

    logical :: convergence_field_scalar

    if(len_trim(sfield%option_path) == 0) then
      convergence_field_scalar = .false.
      return
    end if

    if (aliased(sfield)) then
       convergence_field_scalar=.false.
       return
    end if

    convergence_field_scalar=have_option(trim(complete_field_path(sfield%option_path)) // &
                            "/convergence/include_in_convergence")

  end function convergence_field_scalar

  function convergence_field_vector(vfield, test_for_components)
    !!< Return whether the supplied field should be included in the .convergence file.

    type(vector_field), target, intent(in) :: vfield
    logical, optional, intent(in) :: test_for_components

    logical :: convergence_field_vector

    if(len_trim(vfield%option_path) == 0) then
      convergence_field_vector = .false.
      return
    end if

    if (aliased(vfield)) then
       convergence_field_vector=.false.
       return
    end if

    if(present_and_true(test_for_components)) then
      convergence_field_vector = have_option(trim(complete_field_path(vfield%option_path)) // &
                            "/convergence/include_in_convergence")
    else
      convergence_field_vector=(have_option(trim(complete_field_path(vfield%option_path)) // &
                            "/convergence/include_in_convergence").or.&
                                have_option(trim(complete_field_path(vfield%option_path)) // &
                            "/convergence/exclude_components_from_convergence"))
    end if

  end function convergence_field_vector

  function steady_state_field_scalar(sfield)
    !!< Return whether the supplied field should be checked for a steady state

    type(scalar_field), target, intent(in) :: sfield

    logical :: steady_state_field_scalar

    if(len_trim(sfield%option_path) == 0) then
      steady_state_field_scalar = .false.
      return
    end if

    if (aliased(sfield)) then
       steady_state_field_scalar=.false.
       return
    end if

    steady_state_field_scalar=have_option(trim(complete_field_path(sfield%option_path)) // &
                            "/steady_state/include_in_steady_state")

  end function steady_state_field_scalar

  function steady_state_field_vector(vfield, test_for_components)
    !!< Return whether the supplied field should be checked for a steady state

    type(vector_field), target, intent(in) :: vfield
    logical, optional, intent(in) :: test_for_components

    logical :: steady_state_field_vector

    if(len_trim(vfield%option_path) == 0) then
      steady_state_field_vector = .false.
      return
    end if

    if (aliased(vfield)) then
       steady_state_field_vector=.false.
       return
    end if

    if(present_and_true(test_for_components)) then
      steady_state_field_vector = have_option(trim(complete_field_path(vfield%option_path)) // &
                            "/steady_state/include_in_steady_state")
    else
      steady_state_field_vector=(have_option(trim(complete_field_path(vfield%option_path)) // &
                            "/steady_state/include_in_steady_state").or.&
                                have_option(trim(complete_field_path(vfield%option_path)) // &
                            "/steady_state/exclude_components_from_steady_state"))
    end if

  end function steady_state_field_vector

  function detector_field_scalar(sfield)
    !!< Return whether the supplied field should be included in the .detector file
    logical :: detector_field_scalar
    type(scalar_field), target, intent(in) :: sfield
    
    if (sfield%option_path=="".or.aliased(sfield)) then
       detector_field_scalar=.false.
    else
       detector_field_scalar = have_option(&
            trim(complete_field_path(sfield%option_path)) // &
            "/detectors/include_in_detectors")
    end if

  end function detector_field_scalar

  function detector_field_vector(vfield)
    !!< Return whether the supplied field should be included in the .detector file
    logical :: detector_field_vector
    type(vector_field), target, intent(in) :: vfield

    if (vfield%option_path=="".or.aliased(vfield)) then
       detector_field_vector=.false.
    else
       detector_field_vector = have_option(&
            trim(complete_field_path(vfield%option_path)) // &
            "/detectors/include_in_detectors")
    end if

  end function detector_field_vector

  subroutine initialise_walltime
    !!< Record the initial walltime, clock_rate and maximum clock count

    call system_clock(default_stat%current_count, default_stat%count_rate, default_stat%count_max)
    default_stat%elapsed_count=0

  end subroutine initialise_walltime

  function elapsed_walltime()
    !!< Return the number of walltime seconds since the beginning of the
    !!< simulation.
    real :: elapsed_walltime
    
    integer :: new_count

    call system_clock(new_count)

    ! Deal with clock rollover. If one timestep takes more than a whole
    ! clock rollover, we have more problems than we can deal with!
    if (new_count<default_stat%current_count) then
       default_stat%elapsed_count=default_stat%elapsed_count+(new_count-default_stat%current_count)+default_stat%count_max
    else
       default_stat%elapsed_count=default_stat%elapsed_count+(new_count-default_stat%current_count)
    end if
    default_stat%current_count=new_count

    elapsed_walltime=real(default_stat%elapsed_count)/real(default_stat%count_rate)

  end function elapsed_walltime

  subroutine initialise_diagnostics(filename, state)
    !!< Set up the diagnostic file headers.

    character(len=*) :: filename
    type(state_type), dimension(:), intent(in) :: state

    integer :: column, i, j, k, s, phase, stat
    integer, dimension(2) :: shape_option
    integer :: no_mixing_bins
    real, dimension(:), pointer :: mixing_bin_bounds
    real :: current_time
    character(len = 254) :: buffer, material_phase_name, prefix
    character(len = FIELD_NAME_LEN) :: surface_integral_name, mixing_stats_name
    character(len = OPTION_PATH_LEN) :: func
    type(scalar_field) :: vfield_comp, tfield_comp
    type(mesh_type), pointer :: mesh
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield
    ! Iterator for the registered diagnostics
    type(registered_diagnostic_item), pointer :: iterator => NULL()

    ewrite(1, *) "In initialise_diagnostics"

    ! Idempotency check
    if(default_stat%initialised) then
      ewrite(2, *) "Diagnostics already initialised"
      ewrite(1, *) "Exiting initialise_diagnostics"
      return
    end if
    default_stat%initialised=.true.
    
    ! All processes must assemble the mesh and field lists
    ! Mesh field list
    allocate(default_stat%mesh_list(size(state(1)%mesh_names)))
    default_stat%mesh_list = state(1)%mesh_names
    ! Scalar field list
    allocate (default_stat%sfield_list(size(state)))
    do phase=1, size(state)
       if (associated(state(phase)%scalar_names)) then
          allocate(default_stat%sfield_list(phase)%ptr(size(state(phase)%scalar_names)))
          default_stat%sfield_list(phase)%ptr=state(phase)%scalar_names
       else
          allocate(default_stat%sfield_list(phase)%ptr(0))
       end if
    end do
    ! Vector field list
    allocate (default_stat%vfield_list(size(state)))
    do phase = 1, size(state)
       if (associated(state(phase)%vector_names)) then
          allocate(default_stat%vfield_list(phase)%ptr(size(state(phase)%vector_names)))
          default_stat%vfield_list(phase)%ptr = state(phase)%vector_names
       else
          allocate(default_stat%vfield_list(phase)%ptr(0))
       end if
    end do
    ! Tensor field list
    allocate (default_stat%tfield_list(size(state)))
    do phase = 1, size(state)
       if (associated(state(phase)%tensor_names)) then
          allocate(default_stat%tfield_list(phase)%ptr(size(state(phase)%tensor_names)))
          default_stat%tfield_list(phase)%ptr = state(phase)%tensor_names
       else
          allocate(default_stat%tfield_list(phase)%ptr(0))
       end if
    end do

    ! Only the first process should write statistics information (and hence
    ! write the headers)    
    if(getprocno() == 1) then
      default_stat%diag_unit=free_unit()
      open(unit=default_stat%diag_unit, file=trim(filename)//'.stat', action="write")

      write(default_stat%diag_unit, '(a)') "<header>"

      call initialise_constant_diagnostics(default_stat%diag_unit)

      column=0

      ! Initial columns are elapsed time and dt.
      column=column+1
      buffer=field_tag(name="ElapsedTime", column=column, statistic="value")
      write(default_stat%diag_unit, '(a)') trim(buffer)
      column=column+1
      buffer=field_tag(name="dt", column=column, statistic="value")
      write(default_stat%diag_unit, '(a)') trim(buffer)
      column=column+1
      buffer=field_tag(name="ElapsedWallTime", column=column, statistic="value")
      write(default_stat%diag_unit, '(a)') trim(buffer)

      do i = 1, size(default_stat%mesh_list)
        ! Headers for output statistics for each mesh
        mesh => extract_mesh(state(1), default_stat%mesh_list(i))

        
        if(stat_mesh(mesh)) then
          column = column + 1
          buffer = field_tag(name = mesh%name, column = column, statistic = "nodes")
          write(default_stat%diag_unit, "(a)"), trim(buffer)
          column = column + 1
          buffer = field_tag(name = mesh%name, column = column, statistic = "elements")
          write(default_stat%diag_unit, "(a)"), trim(buffer)
          column = column + 1
          buffer = field_tag(name = mesh%name, column = column, statistic = "surface_elements")
          write(default_stat%diag_unit, "(a)"), trim(buffer)
        end if
      end do

#ifdef HAVE_MEMORY_STATS
      ! Memory statistics
      do i=0, MEMORY_TYPES
          column = column + 1
          buffer = field_tag(name = memory_type_names(i), column = column,&
               & statistic = "current", material_phase_name="Memory")
          write(default_stat%diag_unit, "(a)"), trim(buffer)         
          column = column + 1
          buffer = field_tag(name = memory_type_names(i), column = column,&
               & statistic = "min", material_phase_name="Memory")
          write(default_stat%diag_unit, "(a)"), trim(buffer)         
          column = column + 1
          buffer = field_tag(name = memory_type_names(i), column = column,&
               & statistic = "max", material_phase_name="Memory")
          write(default_stat%diag_unit, "(a)"), trim(buffer)         
      end do
#endif
      
      phaseloop: do phase=1,size(state)

         material_phase_name=trim(state(phase)%name)

         do i=1, size(default_stat%sfield_list(phase)%ptr)
            ! Headers for output statistics for each scalar field
            sfield => extract_scalar_field(state(phase), default_stat%sfield_list(phase)%ptr(i))

            ! Standard scalar field stats
            if(stat_field(sfield, state(phase))) then
              column=column+1
              buffer=field_tag(name=sfield%name, column=column, statistic="min", material_phase_name=material_phase_name)
              write(default_stat%diag_unit, '(a)') trim(buffer)
              column=column+1
              buffer=field_tag(name=sfield%name, column=column, statistic="max", material_phase_name=material_phase_name)
              write(default_stat%diag_unit, '(a)') trim(buffer)
              column=column+1
              buffer=field_tag(name=sfield%name, column=column, statistic="l2norm", material_phase_name=material_phase_name)
              write(default_stat%diag_unit, '(a)') trim(buffer)
              column=column+1
              buffer=field_tag(name=sfield%name, column=column, statistic="integral", material_phase_name=material_phase_name)
              write(default_stat%diag_unit, '(a)') trim(buffer)
            end if
            
            ! Control volume stats
            if(have_option(trim(complete_field_path(sfield%option_path, stat=stat)) // "/stat/include_cv_stats")) then
              column=column+1
              buffer=field_tag(name=sfield%name, column=column, statistic="cv_l2norm", material_phase_name=material_phase_name)
              write(default_stat%diag_unit, '(a)') trim(buffer)
              column=column+1
              buffer=field_tag(name=sfield%name, column=column, statistic="cv_integral", material_phase_name=material_phase_name)
              write(default_stat%diag_unit, '(a)') trim(buffer)
            end if
            
            ! Mixing stats
            do j = 0, option_count(trim(complete_field_path(sfield%option_path, stat=stat)) // "/stat/include_mixing_stats") - 1
              call get_option(trim(complete_field_path(sfield%option_path)) &
              & // "/stat/include_mixing_stats["// int2str(j) // "]/name", mixing_stats_name)
              shape_option=option_shape(trim(complete_field_path(sfield%option_path)) // &
                   & "/stat/include_mixing_stats["// int2str(j) // "]/mixing_bin_bounds")   

              if(have_option(trim(complete_field_path(sfield%option_path)) // &
                  & "/stat/include_mixing_stats["// int2str(j) // "]/mixing_bin_bounds/constant")) then
                  shape_option=option_shape(trim(complete_field_path(sfield%option_path)) // &
                      & "/stat/include_mixing_stats["// int2str(j) // "]/mixing_bin_bounds/constant")
                  no_mixing_bins = shape_option(1)
              else if(have_option(trim(complete_field_path(sfield%option_path)) // &
                  & "/stat/include_mixing_stats["// int2str(j) // "]/mixing_bin_bounds/python")) then
                  call get_option(trim(complete_field_path(sfield%option_path)) // &
                      & "/stat/include_mixing_stats["// int2str(j) // "]/mixing_bin_bounds/python", func)
                  call get_option("/timestepping/current_time", current_time)
                  call real_vector_from_python(func, current_time, mixing_bin_bounds)
                  no_mixing_bins = size(mixing_bin_bounds)
                  deallocate(mixing_bin_bounds)
              else
                  FLExit("Unable to determine mixing bin bounds type. Check options under include_mixing_stats")                  
              end if
          
              buffer = field_tag(name=sfield%name, column=column+1, statistic="mixing_bins%" // trim(mixing_stats_name),&
                   & material_phase_name=material_phase_name, components=(no_mixing_bins))

              write(default_stat%diag_unit, '(a)') trim(buffer)
              column = column + (no_mixing_bins)

            end do
            
            ! Surface integrals
            do j = 0, option_count(trim(complete_field_path(sfield%option_path, stat = stat)) // "/stat/surface_integral") - 1
              call get_option(trim(complete_field_path(sfield%option_path)) &
              // "/stat/surface_integral[" // int2str(j) // "]/name", surface_integral_name)
              column = column + 1
              buffer = field_tag(sfield%name, column, "surface_integral%" // trim(surface_integral_name), material_phase_name)
              write(default_stat%diag_unit, "(a)") trim(buffer)
           end do
           
         end do

         do i = 1, size(default_stat%vfield_list(phase)%ptr)
           ! Headers for output statistics for each vector field
           vfield => extract_vector_field(state(phase), &
             & default_stat%vfield_list(phase)%ptr(i))

           ! Standard scalar field stats for vector field magnitude
           if(stat_field(vfield, state(phase))) then
             column = column + 1
             buffer = field_tag(name=trim(vfield%name) // "%magnitude", column=column, &
               & statistic="min", material_phase_name=material_phase_name)
             write(default_stat%diag_unit, '(a)') trim(buffer)

             column = column + 1
             buffer = field_tag(name=trim(vfield%name) // "%magnitude", column=column, &
               & statistic="max",material_phase_name= material_phase_name)
             write(default_stat%diag_unit, '(a)') trim(buffer)

             column = column + 1
             buffer = field_tag(name=trim(vfield%name) // "%magnitude", column=column, &
               & statistic="l2norm", material_phase_name=material_phase_name)
             write(default_stat%diag_unit, '(a)') trim(buffer)
           end if

           ! Standard scalar field stats for vector field components
           if(stat_field(vfield, state(phase), test_for_components = .true.)) then
             do j = 1, vfield%dim
                vfield_comp = extract_scalar_field(vfield, j)

               column = column + 1
               buffer=field_tag(name=vfield_comp%name, column=column, statistic="min", &
                 & material_phase_name=material_phase_name)
               write(default_stat%diag_unit, '(a)') trim(buffer)

               column = column + 1
               buffer=field_tag(name=vfield_comp%name, column=column, statistic="max", &
                 & material_phase_name=material_phase_name)
               write(default_stat%diag_unit, '(a)') trim(buffer)

               column = column + 1
               buffer=field_tag(name=vfield_comp%name, column=column, statistic="l2norm", &
                 & material_phase_name=material_phase_name)
               write(default_stat%diag_unit, '(a)') trim(buffer)

               column = column + 1
               buffer=field_tag(name=vfield_comp%name, column=column, statistic="integral", &
                 & material_phase_name=material_phase_name)
               write(default_stat%diag_unit, '(a)') trim(buffer)
             end do
           end if

           ! Surface integrals
           do j = 0, option_count(trim(complete_field_path(vfield%option_path, stat = stat)) // "/stat/surface_integral") - 1
             call get_option(trim(complete_field_path(vfield%option_path)) &
             // "/stat/surface_integral[" // int2str(j) // "]/name", surface_integral_name)
             column = column + 1
             buffer = field_tag(vfield%name, column, "surface_integral%" // trim(surface_integral_name), material_phase_name)
             write(default_stat%diag_unit, "(a)") trim(buffer)
           end do

           ! drag calculation
           if(have_option(trim(complete_field_path(vfield%option_path, stat=stat)) // "/stat/compute_body_forces_on_surfaces")) then
             do s = 0, option_count(trim(complete_field_path(vfield%option_path, stat = stat)) // "/stat/compute_body_forces_on_surfaces") - 1
               call get_option(trim(complete_field_path(vfield%option_path))//"/stat/compute_body_forces_on_surfaces[" // int2str(s) // "]/name", surface_integral_name)
               do j = 1, mesh_dim(vfield%mesh)
                 column = column + 1
                 buffer = field_tag(name=trim(vfield%name), column=column, statistic="force_"//trim(surface_integral_name)//"%" &
                 // int2str(j), material_phase_name=material_phase_name)
                 write(default_stat%diag_unit, '(a)') trim(buffer)
               end do
               if(have_option(trim(complete_field_path(vfield%option_path, stat=stat)) // "/stat/compute_body_forces_on_surfaces[" // int2str(s) // "]/output_terms")) then
                 do j = 1, mesh_dim(vfield%mesh)
                   column = column + 1
                   buffer = field_tag(name=trim(vfield%name), column=column, statistic="pressure_force_"//trim(surface_integral_name)//"%" &
                   // int2str(j), material_phase_name=material_phase_name)
                   write(default_stat%diag_unit, '(a)') trim(buffer)
                 end do
                 do j = 1, mesh_dim(vfield%mesh)
                   column = column + 1
                   buffer = field_tag(name=trim(vfield%name), column=column, statistic="viscous_force_"//trim(surface_integral_name)//"%" &
                   // int2str(j), material_phase_name=material_phase_name)
                   write(default_stat%diag_unit, '(a)') trim(buffer)
                 end do
               end if
             end do
           end if

           if(have_option(trim(complete_field_path(vfield%option_path, stat=stat)) // "/stat/divergence_stats")) then
              column=column+1
              buffer=field_tag(name=vfield%name, column=column, statistic="divergence%min", material_phase_name=material_phase_name)
              write(default_stat%diag_unit, '(a)') trim(buffer)
              column=column+1
              buffer=field_tag(name=vfield%name, column=column, statistic="divergence%max", material_phase_name=material_phase_name)
              write(default_stat%diag_unit, '(a)') trim(buffer)
              column=column+1
              buffer=field_tag(name=vfield%name, column=column, statistic="divergence%l2norm", material_phase_name=material_phase_name)
              write(default_stat%diag_unit, '(a)') trim(buffer)
              column=column+1
              buffer=field_tag(name=vfield%name, column=column, statistic="divergence%integral", material_phase_name=material_phase_name)
              write(default_stat%diag_unit, '(a)') trim(buffer)
           end if

           ! momentum conservation error calculation
           if(have_option(trim(complete_field_path(vfield%option_path, stat=stat)) // "/stat/calculate_momentum_conservation_error")) then
             do j = 1, mesh_dim(vfield%mesh)
               column = column + 1
               buffer = field_tag(name=trim(vfield%name), column=column, statistic="momentum_conservation%" &
               // int2str(j), material_phase_name=material_phase_name)
               write(default_stat%diag_unit, '(a)') trim(buffer)
             end do
           end if

         end do

         do i = 1, size(default_stat%tfield_list(phase)%ptr)
           ! Headers for output statistics for each tensor field
           tfield => extract_tensor_field(state(phase), &
             & default_stat%tfield_list(phase)%ptr(i))

           ! Standard scalar field stats for tensor field magnitude
           if(stat_field(tfield, state(phase))) then
             column = column + 1
             buffer = field_tag(name=trim(tfield%name) // "%magnitude", column=column, &
               & statistic="min", material_phase_name=material_phase_name)
             write(default_stat%diag_unit, '(a)') trim(buffer)

             column = column + 1
             buffer = field_tag(name=trim(tfield%name) // "%magnitude", column=column, &
               & statistic="max",material_phase_name= material_phase_name)
             write(default_stat%diag_unit, '(a)') trim(buffer)

             column = column + 1
             buffer = field_tag(name=trim(tfield%name) // "%magnitude", column=column, &
               & statistic="l2norm", material_phase_name=material_phase_name)
             write(default_stat%diag_unit, '(a)') trim(buffer)
           end if

           ! Standard scalar field stats for tensor field components
           if(stat_field(tfield, state(phase), test_for_components = .true.)) then
             do j = 1, tfield%dim(1)
               do k = 1, tfield%dim(2)
                 tfield_comp = extract_scalar_field(tfield, j, k)

                 column = column + 1
                 buffer=field_tag(name=tfield_comp%name, column=column, statistic="min", &
                   & material_phase_name=material_phase_name)
                 write(default_stat%diag_unit, '(a)') trim(buffer)

                 column = column + 1
                 buffer=field_tag(name=tfield_comp%name, column=column, statistic="max", &
                   & material_phase_name=material_phase_name)
                 write(default_stat%diag_unit, '(a)') trim(buffer)

                 column = column + 1
                 buffer=field_tag(name=tfield_comp%name, column=column, statistic="l2norm", &
                   & material_phase_name=material_phase_name)
                 write(default_stat%diag_unit, '(a)') trim(buffer)

                 column = column + 1
                 buffer=field_tag(name=tfield_comp%name, column=column, statistic="integral", &
                   & material_phase_name=material_phase_name)
                 write(default_stat%diag_unit, '(a)') trim(buffer)
               end do
             end do
           end if
         end do

      end do phaseloop
      
      ! Now add the registered diagnostics
      call register_diagnostics
      call print_registered_diagnostics

      iterator => default_stat%registered_diagnostic_first
      do while (associated(iterator)) 
        do i = 1, iterator%dim
          column = column + 1
          if (iterator%dim==1) then
            prefix=""
          else
            prefix=int2str(i)
          end if
          if (iterator%have_material_phase) then
             buffer = field_tag(name=trim(iterator%name)//trim(prefix), column=column, &
                 & statistic=iterator%statistic, material_phase_name=iterator%material_phase)
          else
             buffer = field_tag(name=trim(iterator%name)//trim(prefix), column=column, &
                 & statistic=iterator%statistic)
          end if
          write(default_stat%diag_unit, '(a)') trim(buffer)
        end do  
        iterator => iterator%next
      end do

      write(default_stat%diag_unit, '(a)') "</header>"
      flush(default_stat%diag_unit)
    end if

    call initialise_detectors(filename, state)

    ewrite(1, *) "Exiting initialise_diagnostics"

  end subroutine initialise_diagnostics


  subroutine set_diagnostic(name, statistic, material_phase, value)
    character(len=*), intent(in) :: name, statistic 
    character(len=*), intent(in), optional ::  material_phase
    real, dimension(:), intent(in) :: value
    integer :: i
    
    type(registered_diagnostic_item), pointer :: iterator => NULL()
    
    if(getprocno() == 1) then
      iterator => default_stat%registered_diagnostic_first 

      do while (.true.) 
        if (.not. associated(iterator)) then
          ewrite(0, *) "The diagnostic with name=" // trim(name) //  " statistic=" // trim(statistic)  //  &
               & "material_phase=" //  trim(material_phase) // " does not exist."
          FLAbort("Error in set_diagnostic.")
        end if
        ! Check if name and statistic match
        if (iterator%name == name .and. iterator%statistic == statistic) then
          ! Check if name of material_phase match if supplied
          if ((present(material_phase) .and. iterator%have_material_phase .and. iterator%material_phase == material_phase) &
             & .or. .not. iterator%have_material_phase) then
            ! Check that the value arrays have the same dimension
            if (size(iterator%value) /= size(value)) then
              ewrite(0, *) "The registered diagnostic with name=" // trim(name) // " statistic=" // &
                       & trim(statistic) //  "material_phase=" // trim(material_phase) //  " has dimension " // &
                       & int2str(iterator%dim) // " but a value of dimension " // int2str(size(value))  // & 
                       & " was supplied in set_diagnostic."
              FLAbort("Error in set_diagnostic.")
            end if
            ! set value
            do i = 1, iterator%dim
              iterator%value(i) = value(i)
            end do
            return
          end if
        end if
        iterator => iterator%next
      end do
    end if

  end subroutine set_diagnostic

  function get_diagnostic(name, statistic, material_phase) result(value)
    character(len=*), intent(in) :: name, statistic 
    character(len=*), intent(in), optional ::  material_phase
    real, dimension(:), pointer :: value
    integer :: i
    
    type(registered_diagnostic_item), pointer :: iterator

    iterator => null()
    value => null()
    
    if(getprocno() == 1) then
      iterator => default_stat%registered_diagnostic_first 

      do while (.true.) 
        if (.not. associated(iterator)) then
          ewrite(0, *) "The diagnostic with name=" // trim(name) //  " statistic=" // trim(statistic)  //  &
               & "material_phase=" //  trim(material_phase) // " does not exist."
          FLAbort("Error in set_diagnostic.")
        end if
        ! Check if name and statistic match
        if (iterator%name == name .and. iterator%statistic == statistic) then
          ! Check if name of material_phase match if supplied
          if ((present(material_phase) .and. iterator%have_material_phase .and. iterator%material_phase == material_phase) &
             & .or. .not. iterator%have_material_phase) then
            value => iterator%value
            return
          end if
        end if
        iterator => iterator%next
      end do
    end if
  end function get_diagnostic

  subroutine print_registered_diagnostics
  type(registered_diagnostic_item), pointer :: iterator => NULL()

  iterator => default_stat%registered_diagnostic_first

  ewrite(1, *) "Registered diagnostics:"
  do while(associated(iterator)) 
    if (iterator%have_material_phase) then
      ewrite(1, *) "Name: ", trim(iterator%name),           ", ", &
                 & "Statistic: ", trim(iterator%statistic), ", ", &             
                 & "Dimension: ", int2str(iterator%dim),    ", ", &
                 & "Material phase: ", trim(iterator%material_phase)
    else
      ewrite(1, *) "Name: ", trim(iterator%name),           ", ", &
                 & "Statistic: ", trim(iterator%statistic), ", ", &
                 & "Dimension: ", int2str(iterator%dim)
    end if
    iterator => iterator%next
  end do

  end subroutine print_registered_diagnostics

  subroutine register_diagnostic(dim, name, statistic, material_phase)
    integer, intent(in) :: dim
    character(len=*), intent(in) :: name, statistic 
    character(len=*), intent(in), optional ::  material_phase
    type(registered_diagnostic_item), pointer :: diagnostic_item, iterator => NULL()

    if(getprocno() == 1) then
      ! Allocate the new registered_diagnostic_item and fill it.
      allocate(diagnostic_item)
      diagnostic_item%dim = dim
      diagnostic_item%name = name
      diagnostic_item%statistic = statistic
      if (present(material_phase)) then
        diagnostic_item%material_phase = material_phase
        diagnostic_item%have_material_phase = .true.
      else
        diagnostic_item%have_material_phase = .false.
      end if
      allocate(diagnostic_item%value(dim))
      diagnostic_item%value = INFINITY
      nullify(diagnostic_item%next)

     ! Check if the diagnostic has not been registered yet
      if (associated(default_stat%registered_diagnostic_first)) then
        iterator => default_stat%registered_diagnostic_first
        do while (associated(iterator))
          if (iterator%dim == diagnostic_item%dim .and. iterator%name == diagnostic_item%name .and. &
            & iterator%statistic == diagnostic_item%statistic) then
            if ( (present(material_phase) .and. iterator%have_material_phase) .or. &
               & (.not. present(material_phase) .and. .not. iterator%have_material_phase) ) then
              if (present(material_phase)) then
                if (iterator%material_phase == diagnostic_item%material_phase) then
                  ewrite(0, *) "The diagnostic with name = " // trim(name) // ", statistic = " // &
                      & trim(statistic) //  ", material_phase = " // trim(material_phase) //  ", and dimension = " // &
                      & int2str(iterator%dim) // " has already been registered."
                end if
              else
                ewrite(0, *) "The diagnostic with name = " // trim(name) // ", statistic = " // trim(statistic) &
                     & //  ", and dimension = " // int2str(iterator%dim) // " has already been registered."
              end if
              FLExit("Error in register_diagnostic.")
            end if
          end if
        iterator => iterator%next
        end do
      end if

      ! Now append it to the list of registered diagnostics
      if (.not. associated(default_stat%registered_diagnostic_first)) then
        default_stat%registered_diagnostic_first => diagnostic_item
      else
        iterator => default_stat%registered_diagnostic_first
        do while(associated(iterator%next)) 
          iterator => iterator%next
        end do
        iterator%next => diagnostic_item
      end if
    end if

  end subroutine register_diagnostic

  ! Clean up the list of registered diagnostics
  subroutine destroy_registered_diagnostics
    type(registered_diagnostic_item), pointer :: next, iterator
    
    if(getprocno() == 1) then
      iterator => default_stat%registered_diagnostic_first

      do while (associated(iterator)) 
        next => iterator%next
        deallocate(iterator%value)
        deallocate(iterator)
        iterator => next
      end do
  
      ! the first registered diagnostic needs to be nullified because
      ! when adjointing, all the diagnostics are initialised/registered again
      nullify(default_stat%registered_diagnostic_first)
    end if

  end subroutine destroy_registered_diagnostics

  subroutine uninitialise_diagnostics
  ! Undo all of the initialise_diagnostics business.
  ! Necessary for adjoints, for a start ...
  ! Make sure to call close_diagnostic_files before re-initialising
    type(stat_type) :: new_default_stat

    default_stat%initialised = .false.
    deallocate(default_stat%mesh_list)
    deallocate(default_stat%sfield_list)
    deallocate(default_stat%vfield_list)

    ! The diagnostics are registered under initialise_diagnostics so they need to be destroyed here.
    call destroy_registered_diagnostics

    default_stat = new_default_stat ! new default_stat has been initialised with all the default values, na ja?
  end subroutine uninitialise_diagnostics

  subroutine initialise_constant_diagnostics(unit, binary_format)
    !!< Output constant values in the header of the stat file.
    integer, intent(in) :: unit
    !! If present and .true., indicates binary output format
    logical, optional, intent(in) :: binary_format
    
    character(len=254) :: buffer, value_buffer

#ifdef __FLUIDITY_VERSION__
    value_buffer = __FLUIDITY_VERSION__
#else
    value_buffer="Unknown"
#endif
    buffer=constant_tag(name="FluidityVersion", type="string", value=trim(value_buffer))
    write(unit, '(a)') trim(buffer)
    
    value_buffer = __DATE__ // " " // __TIME__
    buffer=constant_tag(name="CompileTime", type="string", value=trim(value_buffer))
    write(unit, '(a)') trim(buffer)

    value_buffer=date_and_time_string()
    buffer=constant_tag(name="StartTime", type="string", value=trim(value_buffer))
    write(unit, '(a)') trim(buffer)
    
    call get_environment_variable("HOSTNAME", value_buffer, default = "Unknown")
    buffer=constant_tag(name="HostName", type="string", value=trim(value_buffer))
    write(unit, '(a)') trim(buffer)
    
    ! Constant values
    if(present_and_true(binary_format)) then
      buffer = constant_tag(name = "format", type = "string", value = "binary")
      write(unit, '(a)') trim(buffer)
      buffer = constant_tag(name = "real_size", type = "integer", value = int2str(real_size))
      write(unit, '(a)') trim(buffer)
      buffer = constant_tag(name = "integer_size", type = "integer", value = int2str(integer_size))
      write(unit, '(a)') trim(buffer)
    else
      buffer = constant_tag(name = "format", type = "string", value = "plain_text")
      write(unit, '(a)') trim(buffer)
    end if
    
  contains 
    
    function date_and_time_string() 
      character(len=254) :: date_and_time_string

      character(len=8) :: date
      character(len=10) :: time
      character(len=5) :: zone

      call date_and_time(date, time, zone)
      
      date_and_time_string=date//" "//time//zone

    end function date_and_time_string

  end subroutine initialise_constant_diagnostics

  subroutine initialise_convergence(filename, state)
    !!< Set up the convergence file headers.

    character(len=*) :: filename
    type(state_type), dimension(:), intent(in) :: state

    integer :: column, i, j, phase
    character(len = 254) :: buffer, material_phase_name
    type(scalar_field) :: vfield_comp
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield

    ! Idempotency check
    if (default_stat%convergence_initialised) return
    default_stat%convergence_initialised=.true.

    if(have_option("/io/convergence/convergence_file")) then
       default_stat%write_convergence_file = .true.
    else
       default_stat%write_convergence_file = .false.
       return
    end if

    ! Only the first process should write convergence information
    if(getprocno() == 1) then
      default_stat%conv_unit=free_unit()
      open(unit=default_stat%conv_unit, file=trim(filename)//'.convergence', action="write")
    else
      return
    end if

    write(default_stat%conv_unit, '(a)') "<header>"

    column=0

    ! Initial columns are elapsed time, dt and global iteration
    column=column+1
    buffer=field_tag(name="ElapsedTime", column=column, statistic="value")
    write(default_stat%conv_unit, '(a)') trim(buffer)
    column=column+1
    buffer=field_tag(name="dt", column=column, statistic="value")
    write(default_stat%conv_unit, '(a)') trim(buffer)
    column=column+1! Vector field magnitude
    buffer=field_tag(name="Iteration", column=column, statistic="value")
    write(default_stat%conv_unit, '(a)') trim(buffer)

    phaseloop: do phase=1,size(state)

       material_phase_name=trim(state(phase)%name)

       do i=1, size(default_stat%sfield_list(phase)%ptr)
          ! Output convergence information for each scalar field.
          sfield => extract_scalar_field(state(phase), default_stat%sfield_list(phase)%ptr(i))

          if(.not. convergence_field(sfield)) then
            cycle
          end if

          column=column+1
          buffer=field_tag(name=sfield%name, column=column, statistic="error", material_phase_name=material_phase_name)
          write(default_stat%conv_unit, '(a)') trim(buffer)

       end do

         do i = 1, size(default_stat%vfield_list(phase)%ptr)
           ! Headers for output convergence information for each vector field

           vfield => extract_vector_field(state(phase), &
             & default_stat%vfield_list(phase)%ptr(i))

           if(.not. convergence_field(vfield)) then
             cycle
           end if

           column = column + 1
           buffer = field_tag(name=trim(vfield%name) // "%magnitude", column=column, &
             & statistic="error", material_phase_name=material_phase_name)
           write(default_stat%conv_unit, '(a)') trim(buffer)

           if(.not. convergence_field(vfield, test_for_components = .true.)) then
             cycle
           end if

           do j = 1, mesh_dim(vfield%mesh)
             vfield_comp = extract_scalar_field(vfield, j)

             column = column + 1
             buffer=field_tag(name=vfield_comp%name, column=column, statistic="error", &
               & material_phase_name=material_phase_name)
             write(default_stat%conv_unit, '(a)') trim(buffer)

           end do

         end do

    end do phaseloop

    write(default_stat%conv_unit, '(a)') "</header>"

  end subroutine initialise_convergence
  
  subroutine initialise_steady_state(filename, state)
    !!< Set up the steady state file headers.

    character(len=*) :: filename
    type(state_type), dimension(:), intent(in) :: state

    integer :: column, i, j, phase
    character(len = 254) :: buffer, material_phase_name
    type(scalar_field) :: vfield_comp
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield

    ! Idempotency check
    if (default_stat%steady_state_initialised) return
    default_stat%steady_state_initialised=.true.

    default_stat%write_steady_state_file = have_option("/timestepping/steady_state/steady_state_file")
    if(.not. default_stat%write_steady_state_file) return
    if(have_option("/timestepping/steady_state/steady_state_file/binary_output")) then
      default_stat%binary_steady_state_output = .true.
    else if(have_option("/timestepping/steady_state/steady_state_file/plain_text_output")) then
      default_stat%binary_steady_state_output = .false.
    else
      FLExit("Unable to determine steady state output format. Check options under /timestepping/steady_state/steady_state_file")
    end if
    
    ! Only the first process should write steady state information
    if(getprocno() /= 1) return
    
    default_stat%steady_state_unit=free_unit()
    open(unit=default_stat%steady_state_unit, file=trim(filename)//'.steady_state', action="write")

    write(default_stat%steady_state_unit, '(a)') "<header>"

    call initialise_constant_diagnostics(default_stat%steady_state_unit, binary_format = default_stat%binary_steady_state_output)

    ! Initial columns are elapsed time, dt and global iteration
    column=1
    buffer=field_tag(name="ElapsedTime", column=column, statistic="value")
    write(default_stat%steady_state_unit, '(a)') trim(buffer)
    column=column+1
    buffer=field_tag(name="dt", column=column, statistic="value")
    write(default_stat%steady_state_unit, '(a)') trim(buffer)

    phaseloop: do phase=1,size(state)
       material_phase_name = state(phase)%name

       do i = 1, scalar_field_count(state(phase))
          sfield => extract_scalar_field(state(phase), i)
          if(.not. steady_state_field(sfield)) cycle
          ! Scalar fields

          column=column+1
          buffer=field_tag(name=sfield%name, column=column, statistic="error", material_phase_name=material_phase_name)
          write(default_stat%steady_state_unit, '(a)') trim(buffer)
       end do

       do i = 1, vector_field_count(state(phase))
         vfield => extract_vector_field(state(phase), i)
         if(.not. steady_state_field(vfield)) cycle         
         ! Vector fields
         
         column = column + 1
         buffer = field_tag(name=trim(vfield%name) // "%magnitude", column=column, &
           & statistic="error", material_phase_name=material_phase_name)
         write(default_stat%steady_state_unit, '(a)') trim(buffer)

         if(.not. steady_state_field(vfield, test_for_components = .true.)) cycle         
         ! Vector field components
         
         do j = 1, mesh_dim(vfield%mesh)
           vfield_comp = extract_scalar_field(vfield, j)

           column = column + 1
           buffer=field_tag(name=vfield_comp%name, column=column, statistic="error", &
             & material_phase_name=material_phase_name)
           write(default_stat%steady_state_unit, '(a)') trim(buffer)
         end do
       end do

    end do phaseloop

    column = column + 1
    buffer = field_tag(name = "MaxChange", column=column, statistic="value")
    write(default_stat%steady_state_unit, '(a)') trim(buffer)

    write(default_stat%steady_state_unit, '(a)') "</header>"
    flush(default_stat%steady_state_unit)
    
    if(default_stat%binary_steady_state_output) then
        close(default_stat%steady_state_unit)

#ifdef STREAM_IO
      open(unit = default_stat%steady_state_unit, file = trim(filename) // '.steady_state.dat', &
        & action = "write", access = "stream", form = "unformatted", status = "replace")
#else
      FLAbort("No stream I/O support")
#endif
    end if

  end subroutine initialise_steady_state

  subroutine create_single_detector(detector_list,xfield,position,id,type,name)
    ! Allocate a single detector, populate and insert it into the given list
    ! In parallel, first check if the detector would be local and only allocate if it is
    type(detector_linked_list), intent(inout) :: detector_list
    type(vector_field), pointer :: xfield
    real, dimension(xfield%dim), intent(in) :: position
    integer, intent(in) :: id, type
    character(len=*), intent(in) :: name

    type(detector_type), pointer :: detector
    type(element_type), pointer :: shape
    real, dimension(xfield%dim+1) :: lcoords
    integer :: element

    shape=>ele_shape(xfield,1)
    assert(xfield%dim+1==local_coord_count(shape))

    ! Determine element and local_coords from position
    ! In parallel, global=.false. can often work because there will be
    ! a halo of non-owned elements in your process and so you can work out
    ! ownership without communication.  But in general it won't work.
    call picker_inquire(xfield,position,element,local_coord=lcoords,global=.true.)

    ! If we're in parallel and don't own the element, skip this detector
    if (isparallel()) then
       if (element<0) return
       if (.not.element_owned(xfield,element)) return
    else
       ! In serial make sure the detector is in the domain
       ! unless we have the write_nan_outside override
       if (element<0 .and. .not.detector_list%write_nan_outside) then
          ewrite(-1,*) "Dealing with detector ", id, " named: ", trim(name)
          FLExit("Trying to initialise detector outside of computational domain")
       end if
    end if
         
    ! Otherwise, allocate and insert detector
    allocate(detector)
    allocate(detector%position(xfield%dim))
    allocate(detector%local_coords(local_coord_count(shape)))
    call insert(detector,default_stat%detector_list)

    ! Populate detector
    detector%name=name
    detector%position=position
    detector%element=element
    detector%local_coords=lcoords
    detector%type=type
    detector%id_number=id

  end subroutine create_single_detector
  
  subroutine initialise_detectors(filename, state)
    !!< Set up the detector file headers. This has the same syntax as the
    !!< .stat file
    character(len = *), intent(in) :: filename
    type(state_type), dimension(:), intent(in) :: state

    character(len=FIELD_NAME_LEN) ::funcnam, temp_name, detector_name
    character(len=PYTHON_FUNC_LEN) :: func

    integer :: column, i, j, k, phase, m, IERROR, field_count, totaldet_global
    integer :: static_dete, python_functions_or_files, total_dete, total_dete_groups, lagrangian_dete
    integer :: python_dete, ndete, dim, str_size, type_det
    integer, dimension(2) :: shape_option
    character(len = 254) :: buffer, material_phase_name, fmt
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield, xfield
    real, allocatable, dimension(:,:) :: coords
    real, allocatable, dimension(:) :: detector_location
    real:: current_time
    character(len = OPTION_PATH_LEN) :: detectors_cp_filename, detector_file_filename

    type(detector_type), pointer :: detector
    type(element_type), pointer :: shape

    ! Idempotency check
    if (default_stat%detectors_initialised) return
    default_stat%detectors_initialised=.true.

    ewrite(2,*) "In initialise_detectors"

    ! Check whether there are actually any detectors.
    static_dete = option_count("/io/detectors/static_detector")
    lagrangian_dete = option_count("/io/detectors/lagrangian_detector")
    python_functions_or_files = option_count("/io/detectors/detector_array")
    python_dete = 0
 
    do i=1,python_functions_or_files
       write(buffer, "(a,i0,a)") "/io/detectors/detector_array[",i-1,"]"
       call get_option(trim(buffer)//"/number_of_detectors", j)
       python_dete=python_dete+j
    end do
   
    total_dete=static_dete+lagrangian_dete+python_dete
    default_stat%detector_list%total_num_det=total_dete

    total_dete_groups=static_dete+lagrangian_dete+python_functions_or_files

    allocate(default_stat%detector_group_names(total_dete_groups))
    allocate(default_stat%number_det_in_each_group(total_dete_groups))
    allocate(default_stat%detector_list%detector_names(total_dete))
    
    if (total_dete==0) return

    ! Register this I/O detector list with a global list of detector lists
    call register_detector_list(default_stat%detector_list)

    xfield=>extract_vector_field(state(1), "Coordinate")
    shape=>ele_shape(xfield,1)
    call get_option("/geometry/dimension",dim)
    call get_option("/timestepping/current_time", current_time)
    allocate(detector_location(dim))

    ! Enable detectors to drift with the mesh
    if (have_option("/io/detectors/move_with_mesh")) then
       default_stat%detector_list%move_with_mesh=.true.
    end if

    ! Set flag for NaN detector output
    if (have_option("/io/detectors/write_nan_outside_domain")) then
       default_stat%detector_list%write_nan_outside=.true.
    end if
    
    ! Retrieve the position of each detector. If the option
    ! "from_checkpoint_file" exists, it means we are continuing the simulation
    ! after checkpointing and the reading of the detector positions must be
    ! done from a file
    if (have_option("/io/detectors/static_detector/from_checkpoint_file").or. & 
         & have_option("/io/detectors/lagrangian_detector/from_checkpoint_file").or. &
         & have_option("/io/detectors/detector_array/from_checkpoint_file")) then
       default_stat%from_checkpoint=.true.
    else
       default_stat%from_checkpoint=.false.
    end if

    ! Read detectors from options
    if (.not.default_stat%from_checkpoint) then
       ewrite(2,*) "Reading detectors from options"

       ! Read all single static detector from options
       do i=1,static_dete
          write(buffer, "(a,i0,a)") "/io/detectors/static_detector[",i-1,"]"

          shape_option=option_shape(trim(buffer)//"/location")
          assert(xfield%dim==shape_option(1))
          call get_option(trim(buffer)//"/location", detector_location)

          ! The arrays below contain information about the order in which detector
          ! groups are read and how many detectors there are in each group. This is
          ! used when checkpointing detectors. In particular, when continuing a
          ! simulation from a checkpoint, with these arrays we make sure we read
          ! back the detectors from the file in the same order than at the beginning
          ! of the simulation for consistency. All the .detectors files with
          ! detector data (position, value of variables at those positions, etc.) 
          ! will have the information in the same order.
          call get_option(trim(buffer)//"/name", detector_name)
          default_stat%detector_group_names(i)=detector_name
          default_stat%number_det_in_each_group(i)=1.0
          default_stat%detector_list%detector_names(i)=detector_name

          call create_single_detector(default_stat%detector_list, xfield, &
                detector_location, i, STATIC_DETECTOR, trim(detector_name))
       end do

       ! Read all single lagrangian detector from options
       do i=1,lagrangian_dete
          write(buffer, "(a,i0,a)") "/io/detectors/lagrangian_detector[",i-1,"]"

          shape_option=option_shape(trim(buffer)//"/location")
          assert(xfield%dim==shape_option(1))
          call get_option(trim(buffer)//"/location", detector_location)

          call get_option(trim(buffer)//"/name", detector_name)
          default_stat%detector_group_names(static_dete+i)=detector_name
          default_stat%number_det_in_each_group(static_dete+i)=1.0
          default_stat%detector_list%detector_names(static_dete+i)=detector_name

          call create_single_detector(default_stat%detector_list, xfield, &
                detector_location, static_dete+i, LAGRANGIAN_DETECTOR, trim(detector_name))
       end do

       k=static_dete+lagrangian_dete+1

       do i=1,python_functions_or_files
          write(buffer, "(a,i0,a)") "/io/detectors/detector_array[",i-1,"]"

          call get_option(trim(buffer)//"/name", funcnam)
          call get_option(trim(buffer)//"/number_of_detectors", ndete)
          str_size=len_trim(int2str(ndete))
          fmt="(a,I"//int2str(str_size)//"."//int2str(str_size)//")"

          if (have_option(trim(buffer)//"/lagrangian")) then
             type_det=LAGRANGIAN_DETECTOR
          else
             type_det=STATIC_DETECTOR
          end if

          default_stat%detector_group_names(i+static_dete+lagrangian_dete)=trim(funcnam)
          default_stat%number_det_in_each_group(i+static_dete+lagrangian_dete)=ndete

          if (.not.have_option(trim(buffer)//"/from_file")) then

             ! Reading detectors from a python function
             call get_option(trim(buffer)//"/python", func)
             allocate(coords(dim,ndete))
             call set_detector_coords_from_python(coords, ndete, func, current_time)
          
             do j=1,ndete
                write(detector_name, fmt) trim(funcnam)//"_", j
                default_stat%detector_list%detector_names(k)=trim(detector_name)

                call create_single_detector(default_stat%detector_list, xfield, &
                       coords(:,j), k, type_det, trim(detector_name))
                k=k+1           
             end do
             deallocate(coords)

          else

             ! Reading from a binary file where the user has placed the detector positions
             default_stat%detector_file_unit=free_unit()
             call get_option("/io/detectors/detector_array/from_file/file_name",detector_file_filename)

#ifdef STREAM_IO
             open(unit = default_stat%detector_file_unit, file = trim(detector_file_filename), &
                  & action = "read", access = "stream", form = "unformatted")
#else
             FLAbort("No stream I/O support")
#endif
   
             do j=1,ndete
                write(detector_name, fmt) trim(funcnam)//"_", j
                default_stat%detector_list%detector_names(k)=trim(detector_name)
                read(default_stat%detector_file_unit) detector_location
                call create_single_detector(default_stat%detector_list, xfield, &
                      detector_location, k, type_det, trim(detector_name))
                k=k+1          
             end do
          end if
       end do      
    else 
       ewrite(2,*) "Reading detectors from checkpoint"

       ! If reading from checkpoint file:
       ! Detector checkpoint file names end in _det, with.groups appended for the header file
       ! and .positions.dat appended for the binary data file that holds the positions

       default_stat%detector_checkpoint_unit=free_unit()
       if (have_option("/io/detectors/static_detector")) then
          call get_option("/io/detectors/static_detector/from_checkpoint_file/file_name",detectors_cp_filename)  
       elseif (have_option("/io/detectors/lagrangian_detector")) then 
          call get_option("/io/detectors/lagrangian_detector/from_checkpoint_file/file_name",detectors_cp_filename)  
       else 
          call get_option("/io/detectors/detector_array/from_checkpoint_file/file_name",detectors_cp_filename)  
       end if 

       open(unit=default_stat%detector_checkpoint_unit, file=trim(detectors_cp_filename) // '.groups', action="read") 

       ! First we read the header of checkpoint_file to get the order in which the detectors were read initialliy
       do i=1, total_dete_groups 
          read(default_stat%detector_checkpoint_unit,'(a,i10)') default_stat%detector_group_names(i), default_stat%number_det_in_each_group(i)
       end do

       close(default_stat%detector_checkpoint_unit)

#ifdef STREAM_IO
       open(unit = default_stat%detector_checkpoint_unit, file = trim(detectors_cp_filename) // '.positions.dat', &
             & action = "read", access = "stream", form = "unformatted")
#else
       FLAbort("No stream I/O support")
#endif
 
       ! Read in order the last positions of the detectors from the binary file.

       do j=1,size(default_stat%detector_group_names)
          do i=1,static_dete
             write(buffer, "(a,i0,a)") "/io/detectors/static_detector[",i-1,"]"
             call get_option(trim(buffer)//"/name", temp_name)
       
             if (default_stat%detector_group_names(j)==temp_name) then
                read(default_stat%detector_checkpoint_unit) detector_location
                call create_single_detector(default_stat%detector_list, xfield, &
                      detector_location, i, STATIC_DETECTOR, trim(temp_name))                  
             else
                cycle
             end if
          end do
       end do

       do j=1,size(default_stat%detector_group_names)
          do i=1,lagrangian_dete
             write(buffer, "(a,i0,a)") "/io/detectors/lagrangian_detector[",i-1,"]"
             call get_option(trim(buffer)//"/name", temp_name)

             if (default_stat%detector_group_names(j)==temp_name) then
                read(default_stat%detector_checkpoint_unit) detector_location
                call create_single_detector(default_stat%detector_list, xfield, &
                      detector_location, i+static_dete, LAGRANGIAN_DETECTOR, trim(temp_name)) 
             else
                cycle
             end if
          end do
       end do

       k=static_dete+lagrangian_dete+1

       do j=1,size(default_stat%detector_group_names) 
          do i=1,python_functions_or_files
             write(buffer, "(a,i0,a)") "/io/detectors/detector_array[",i-1,"]"      
             call get_option(trim(buffer)//"/name", temp_name)

             if (default_stat%detector_group_names(j)==temp_name) then
                call get_option(trim(buffer)//"/number_of_detectors", ndete)
                str_size=len_trim(int2str(ndete))
                fmt="(a,I"//int2str(str_size)//"."//int2str(str_size)//")"

                if (have_option(trim(buffer)//"/lagrangian")) then
                   type_det=LAGRANGIAN_DETECTOR
                else
                   type_det=STATIC_DETECTOR
                end if

                do m=1,default_stat%number_det_in_each_group(j)
                   write(detector_name, fmt) trim(temp_name)//"_", m
                   read(default_stat%detector_checkpoint_unit) detector_location
                   call create_single_detector(default_stat%detector_list, xfield, &
                          detector_location, k, type_det, trim(detector_name)) 
                   default_stat%detector_list%detector_names(k)=trim(detector_name)
                   k=k+1           
                end do
             else                     
                cycle                   
             end if             
          end do
       end do

    end if  ! from_checkpoint

    default_stat%detector_list%binary_output = .true.
    if (have_option("/io/detectors/ascii_output")) then
       default_stat%detector_list%binary_output = .false.
       if(isparallel()) then
          FLAbort("No support for ascii detector output in parallel. Please use binary output.")
       end if
    end if

    ! Only the first process should write the header file
    if (getprocno() == 1) then
       default_stat%detector_list%output_unit=free_unit()
       open(unit=default_stat%detector_list%output_unit, file=trim(filename)//'.detectors', action="write")

       write(default_stat%detector_list%output_unit, '(a)') "<header>"
       call initialise_constant_diagnostics(default_stat%detector_list%output_unit, &
                binary_format = default_stat%detector_list%binary_output)

       ! Initial columns are elapsed time and dt.
       buffer=field_tag(name="ElapsedTime", column=1, statistic="value")
       write(default_stat%detector_list%output_unit, '(a)') trim(buffer)
       buffer=field_tag(name="dt", column=2, statistic="value")
       write(default_stat%detector_list%output_unit, '(a)') trim(buffer)

       ! Next columns contain the positions of all the detectors.
       column=2
       positionloop: do i=1, default_stat%detector_list%total_num_det
          buffer=field_tag(name=default_stat%detector_list%detector_names(i), column=column+1,&
              statistic="position", components=xfield%dim)
          write(default_stat%detector_list%output_unit, '(a)') trim(buffer)
          column=column+xfield%dim   ! xfield%dim == size(detector%position)
       end do positionloop
     end if

     ! Loop over all fields in state and record the ones we want to output
     allocate (default_stat%detector_list%sfield_list(size(state)))
     allocate (default_stat%detector_list%vfield_list(size(state)))
     phaseloop: do phase=1,size(state)
        material_phase_name=trim(state(phase)%name)

        ! Count the scalar fields to include in detectors
        field_count = 0
        do i = 1, size(state(phase)%scalar_names)
           sfield => extract_scalar_field(state(phase),state(phase)%scalar_names(i))   
           if (detector_field(sfield)) field_count = field_count + 1
        end do 
        allocate(default_stat%detector_list%sfield_list(phase)%ptr(field_count))
        default_stat%detector_list%num_sfields=default_stat%detector_list%num_sfields + field_count

        ! Loop over scalar fields again to store names and create header lines
        field_count = 1
        do i=1, size(state(phase)%scalar_names)

           sfield => extract_scalar_field(state(phase),state(phase)%scalar_names(i))
           if(.not. detector_field(sfield)) then
              cycle
           end if

           ! Create header for included scalar field (first proc only)
           if (getprocno() == 1) then
             do j=1, default_stat%detector_list%total_num_det
                column=column+1
                buffer=field_tag(name=sfield%name, column=column, &
                    statistic=default_stat%detector_list%detector_names(j), &
                    material_phase_name=material_phase_name)
                write(default_stat%detector_list%output_unit, '(a)') trim(buffer)
             end do
           end if

           ! Store name of included scalar field
           default_stat%detector_list%sfield_list(phase)%ptr(field_count)=state(phase)%scalar_names(i)
           field_count = field_count + 1
        end do

        ! Count the vector fields to include in detectors
        field_count = 0
        do i = 1, size(state(phase)%vector_names)
           vfield => extract_vector_field(state(phase),state(phase)%vector_names(i))   
           if (detector_field(vfield)) field_count = field_count + 1
        end do 
        allocate(default_stat%detector_list%vfield_list(phase)%ptr(field_count))
        default_stat%detector_list%num_vfields=default_stat%detector_list%num_vfields + field_count

        ! Loop over vector fields again to store names and create header lines
        field_count = 1
        do i=1, size(state(phase)%vector_names)

           vfield => extract_vector_field(state(phase),state(phase)%vector_names(i))
           if(.not. detector_field(vfield)) then
              cycle
           end if

           ! Create header for included vector field (first proc only)
           if (getprocno() == 1) then
             do j=1, default_stat%detector_list%total_num_det
                buffer=field_tag(name=vfield%name, column=column+1, &
                    statistic=default_stat%detector_list%detector_names(j), &
                    material_phase_name=material_phase_name, &
                    components=vfield%dim)
                write(default_stat%detector_list%output_unit, '(a)') trim(buffer)
                column=column+vfield%dim
             end do
           end if

           ! Store name of included vector field
           default_stat%detector_list%vfield_list(phase)%ptr(field_count)=state(phase)%vector_names(i)
           field_count = field_count + 1
        end do

     end do phaseloop

     if (getprocno() == 1) then
       write(default_stat%detector_list%output_unit, '(a)') "</header>"
       flush(default_stat%detector_list%output_unit)
       ! when using mpi_subroutines to write into the detectors file we need to close the file since 
       ! filename.detectors.dat needs to be open now with MPI_OPEN
       if ((isparallel()).or.(default_stat%detector_list%binary_output)) then
          close(default_stat%detector_list%output_unit)
       end if
    end if  

    if ((isparallel()).or.(default_stat%detector_list%binary_output)) then
       ! bit of hack to delete any existing .detectors.dat file
       ! if we don't delete the existing .detectors.dat would simply be opened for random access and 
       ! gradually overwritten, mixing detector output from the current with that of a previous run
       call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(filename) // '.detectors.dat', MPI_MODE_CREATE + MPI_MODE_RDWR + MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, default_stat%detector_list%mpi_fh, IERROR)
       call MPI_FILE_CLOSE(default_stat%detector_list%mpi_fh, IERROR)
       call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(filename) // '.detectors.dat', MPI_MODE_CREATE + MPI_MODE_RDWR, MPI_INFO_NULL, default_stat%detector_list%mpi_fh, IERROR)
       assert(ierror == MPI_SUCCESS)
    end if 

    !Get options for lagrangian detector movement
    if (check_any_lagrangian(default_stat%detector_list)) then
       call read_detector_move_options(default_stat%detector_list, "/io/detectors")
    end if

    ! And finally some sanity checks
    totaldet_global=default_stat%detector_list%length
    call allsum(totaldet_global)
    ewrite(2,*) "Found", default_stat%detector_list%length, "local and ", totaldet_global, "global detectors"

    assert(totaldet_global==default_stat%detector_list%total_num_det)

  end subroutine initialise_detectors

  function field_tag(name, column, statistic, material_phase_name, components)
    !!< Create a field tag for the given entries.
    character(len=*), intent(in) :: name
    integer, intent(in) :: column
    character(len=*), intent(in) :: statistic
    character(len=*), intent(in), optional :: material_phase_name 
    integer, intent(in), optional :: components
    character(len=254) :: field_tag

    character(len=254) :: front_buffer, material_buffer, components_buffer, end_buffer

    write(front_buffer,'(a,i0,a)') '<field column="',column,'" name="'&
            &//trim(name)//'" statistic="'//trim(statistic)//'"'

    if (present(material_phase_name)) then
        write(material_buffer,'(a)') ' material_phase="'//&
            trim(material_phase_name)//'"'
    else
        material_buffer = ''
    end if

    if (present(components)) then
        write(components_buffer,'(a,i0,a)') ' components="', components, '"' 
    else
        components_buffer = ''
    end if

    end_buffer = '/>'

    field_tag = trim(front_buffer)//trim(material_buffer)//trim(components_buffer)//trim(end_buffer)

  end function field_tag

  function constant_tag(name, type, value)
    !!< Create a field tag for the given entries.
    character(len=*), intent(in) :: name, type, value
    
    character(len=254) :: constant_tag
    
    constant_tag='<constant name="'//trim(name)&
         &//'" type="'//trim(type)//'" value="'//trim(value)//'" />'

  end function constant_tag
  
  subroutine write_diagnostics(state, time, dt, timestep, not_to_move_det_yet)
    !!< Write the diagnostics to the previously opened diagnostics file.
    type(state_type), dimension(:), intent(inout) :: state
    real, intent(in) :: time, dt
    integer, intent(in) :: timestep
    logical, intent(in), optional :: not_to_move_det_yet 

    character(len = 2 + real_format_len(padding = 1) + 1) :: format, format2, format3, format4
    character(len = OPTION_PATH_LEN) :: func
    integer :: i, j, k, phase, stat
    integer, dimension(2) :: shape_option
    integer :: nodes, elements, surface_elements
    integer :: no_mixing_bins
    real :: fmin, fmax, fnorm2, fintegral, fnorm2_cv, fintegral_cv, surface_integral
    real, dimension(:), allocatable :: f_mix_fraction
    real, dimension(:), pointer :: mixing_bin_bounds
    real :: current_time
    type(mesh_type), pointer :: mesh
    type(scalar_field) :: vfield_comp, tfield_comp
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield
    type(vector_field) :: xfield
    type(scalar_field), pointer :: cv_mass => null()
    type(registered_diagnostic_item), pointer :: iterator => NULL()
    logical :: l_move_detectors

    ewrite(1,*) 'In write_diagnostics'
    call profiler_tic("I/O")

    if(present_and_true(not_to_move_det_yet)) then
       l_move_detectors=.false.
    else
       l_move_detectors=.true.
    end if

    format="(" // real_format(padding = 1) // ")"
    format2="(2" // real_format(padding = 1) // ")"
    format3="(3" // real_format(padding = 1) // ")"
    format4="(4" // real_format(padding = 1) // ")"

    ! Only the first process should write statistics information (but all must
    ! be involved in calculating them)
    if(getprocno() == 1) then
      write(default_stat%diag_unit, trim(format), advance="no") time
      write(default_stat%diag_unit, trim(format), advance="no") dt
      write(default_stat%diag_unit, trim(format), advance="no") elapsed_walltime()
    end if

    do i = 1, size(default_stat%mesh_list)
      ! Output statistics for each mesh
      mesh => extract_mesh(state(1), default_stat%mesh_list(i))

      if(stat_mesh(mesh)) then
        call mesh_stats(mesh, nodes, elements, surface_elements)
        if(getprocno() == 1) then
          write(default_stat%diag_unit, "(a,i0,a,i0,a,i0)", advance = "no") " ", nodes, " ", elements, " ", surface_elements
        end if
      end if
    end do

#ifdef HAVE_MEMORY_STATS
    ! Memory statistics.
    call write_memory_stats(default_stat%diag_unit, format)
    call reset_memory_logs
#endif

    phaseloop: do phase=1,size(state)

       scalar_field_loop: do i=1, size(default_stat%sfield_list(phase)%ptr)
          ! Output statistics for each scalar field
          sfield=>extract_scalar_field(state(phase), default_stat%sfield_list(phase)%ptr(i))

          xfield=get_diagnostic_coordinate_field(state(phase), sfield%mesh)

          ! Standard scalar field stats
          if(stat_field(sfield, state(phase))) then
          
            call field_stats(sfield, Xfield, fmin, fmax, fnorm2, fintegral)
            if(getprocno() == 1) then
              write(default_stat%diag_unit, trim(format4), advance="no") fmin, fmax, fnorm2,&
                   & fintegral
            end if
            
          end if

          ! Control volume stats
          if(have_option(trim(complete_field_path(sfield%option_path,stat=stat)) //&
               & "/stat/include_cv_stats")) then
            
            ! Get the CV mass matrix 
            cv_mass => get_cv_mass(state(phase), sfield%mesh)

            call field_cv_stats(sfield, cv_mass, fnorm2_cv, fintegral_cv)

            ! Only the first process should write statistics information
            if(getprocno() == 1) then
              write(default_stat%diag_unit, trim(format2), advance="no") fnorm2_cv, fintegral_cv
            end if

          end if

          ! Mixing stats
          do j = 0, option_count(trim(complete_field_path(sfield%option_path, stat=stat)) // "/stat/include_mixing_stats") - 1
            if(have_option(trim(complete_field_path(sfield%option_path)) // &
                  & "/stat/include_mixing_stats["// int2str(j) // "]/mixing_bin_bounds/constant")) then
                  shape_option=option_shape(trim(complete_field_path(sfield%option_path)) // &
                      & "/stat/include_mixing_stats["// int2str(j) // "]/mixing_bin_bounds/constant")
                  ewrite(1,*) 'shape_option = ', shape_option
                  no_mixing_bins = shape_option(1)
                
            else if(have_option(trim(complete_field_path(sfield%option_path)) // &
                & "/stat/include_mixing_stats["// int2str(j) // "]/mixing_bin_bounds/python")) then
                call get_option(trim(complete_field_path(sfield%option_path)) // &
                    & "/stat/include_mixing_stats["// int2str(j) // "]/mixing_bin_bounds/python", func)
                call get_option("/timestepping/current_time", current_time)
                call real_vector_from_python(func, current_time, mixing_bin_bounds)
                no_mixing_bins = size(mixing_bin_bounds)
                deallocate(mixing_bin_bounds)
            else
                FLExit("Unable to determine mixing bin bounds type. Check options under include_mixing_stats")                  
            end if
                    
            allocate(f_mix_fraction(1:no_mixing_bins))
            f_mix_fraction = 0.0
            
            call mixing_stats(f_mix_fraction, sfield, Xfield, mixing_stats_count = j)          
         
            if(getprocno() == 1) then
               do k=1, (size(f_mix_fraction))
                  write(default_stat%diag_unit, trim(format), advance="no") f_mix_fraction(k)
               end do
            end if

            deallocate(f_mix_fraction)

          end do
         
         ! Surface integrals
         do j = 0, option_count(trim(complete_field_path(sfield%option_path, stat = stat)) // "/stat/surface_integral") - 1
           surface_integral = calculate_surface_integral(sfield, xfield, j)
           ! Only the first process should write statistics information
           if(getprocno() == 1) then
             write(default_stat%diag_unit, trim(format), advance = "no") surface_integral
           end if
         end do

         call deallocate(xfield)

       end do scalar_field_loop

       vector_field_loop: do i = 1, size(default_stat%vfield_list(phase)%ptr)
         ! Output statistics for each vector field
         vfield => extract_vector_field(state(phase), &
           & default_stat%vfield_list(phase)%ptr(i))
          
         xfield=get_diagnostic_coordinate_field(state(phase), vfield%mesh)

         ! Standard scalar field stats for vector field magnitude
         if(stat_field(vfield,state(phase))) then
           call field_stats(vfield, Xfield, fmin, fmax, fnorm2)
           ! Only the first process should write statistics information
           if(getprocno() == 1) then
             write(default_stat%diag_unit, trim(format3), advance = "no") fmin, fmax, fnorm2
           end if
         end if

         ! Standard scalar field stats for vector field components
         if(stat_field(vfield, state(phase), test_for_components = .true.)) then
           do j = 1, vfield%dim
             vfield_comp = extract_scalar_field(vfield, j)

             call field_stats(vfield_comp, Xfield, fmin, fmax, fnorm2, &
               & fintegral)
             ! Only the first process should write statistics information
             if(getprocno() == 1) then
               write(default_stat%diag_unit, trim(format4), advance = "no") fmin, fmax, fnorm2, &
                 & fintegral
             end if
           end do
         end if
         
         ! Surface integrals
         do j = 0, option_count(trim(complete_field_path(vfield%option_path, stat = stat)) // "/stat/surface_integral") - 1
           surface_integral = calculate_surface_integral(vfield, xfield, j)
           ! Only the first process should write statistics information
           if(getprocno() == 1) then
             write(default_stat%diag_unit, trim(format), advance = "no") surface_integral
           end if
         end do

         ! drag calculation
         if(have_option(trim(complete_field_path(vfield%option_path, stat=stat)) // "/stat/compute_body_forces_on_surfaces")) then
           call write_body_forces(state(phase), vfield)  
         end if

         if(have_option(trim(complete_field_path(vfield%option_path, stat=stat)) // "/stat/divergence_stats")) then
           call divergence_field_stats(vfield, Xfield, fmin, fmax, fnorm2, fintegral)
           if(getprocno() == 1) then
             write(default_stat%diag_unit, trim(format4), advance="no") fmin, fmax, fnorm2,&
                   & fintegral
           end if            
         end if

         ! momentum conservation error calculation
         if(have_option(trim(complete_field_path(vfield%option_path, stat=stat)) // "/stat/calculate_momentum_conservation_error")) then
           call write_momentum_conservation_error(state(phase), vfield)
         end if
         
         call deallocate(xfield)
         
       end do vector_field_loop

       tensor_field_loop: do i = 1, size(default_stat%tfield_list(phase)%ptr)
         ! Output statistics for each tensor field
         tfield => extract_tensor_field(state(phase), &
           & default_stat%tfield_list(phase)%ptr(i))
          
         xfield=get_diagnostic_coordinate_field(state(phase), tfield%mesh)

         ! Standard scalar field stats for tensor field magnitude
         if(stat_field(tfield,state(phase))) then
           call field_stats(tfield, Xfield, fmin, fmax, fnorm2)
           ! Only the first process should write statistics information
           if(getprocno() == 1) then
             write(default_stat%diag_unit, trim(format3), advance = "no") fmin, fmax, fnorm2
           end if
         end if

         ! Standard scalar field stats for tensor field components
         if(stat_field(tfield, state(phase), test_for_components = .true.)) then
           do j = 1, tfield%dim(1)
             do k = 1, tfield%dim(2)
               tfield_comp = extract_scalar_field(tfield, j, k)

               call field_stats(tfield_comp, Xfield, fmin, fmax, fnorm2, &
                 & fintegral)
               ! Only the first process should write statistics information
               if(getprocno() == 1) then
                 write(default_stat%diag_unit, trim(format4), advance = "no") fmin, fmax, fnorm2, &
                   & fintegral
               end if
             end do
           end do
         end if
         
         call deallocate(xfield)
         
       end do tensor_field_loop

    end do phaseloop

    ! Registered diagnostics
    call print_registered_diagnostics
    iterator => default_stat%registered_diagnostic_first
    do while (associated(iterator)) 
      ! Only the first process should write statistics information
      if(getprocno() == 1) then   
        do k=1, iterator%dim
          write(default_stat%diag_unit, trim(format), advance = "no") iterator%value(k)
        end do    
      end if
      iterator => iterator%next
    end do

    ! Output end of line
    ! Only the first process should write statistics information
    if(getprocno() == 1) then
      write(default_stat%diag_unit,'(a)') ""
      flush(default_stat%diag_unit)
    end if

    ! Move lagrangian detectors
    if ((timestep/=0).and.l_move_detectors.and.check_any_lagrangian(default_stat%detector_list)) then
       call move_lagrangian_detectors(state, default_stat%detector_list, dt, timestep)
    end if

    ! Now output any detectors.    
    call write_detectors(state, default_stat%detector_list, time, dt)

    call profiler_toc("I/O")
  
  contains
  
    subroutine write_body_forces(state, vfield)
      type(state_type), intent(in) :: state
      type(vector_field), intent(in) :: vfield
      type(tensor_field), pointer :: viscosity

      logical :: have_viscosity      
      integer :: i, s
      real :: force(vfield%dim), pressure_force(vfield%dim), viscous_force(vfield%dim)
      character(len = FIELD_NAME_LEN) :: surface_integral_name
    
      viscosity=>extract_tensor_field(state, "Viscosity", stat)
      have_viscosity = stat == 0

      do s = 0, option_count(trim(complete_field_path(vfield%option_path, stat = stat)) // "/stat/compute_body_forces_on_surfaces") - 1
        call get_option(trim(complete_field_path(vfield%option_path))//"/stat/compute_body_forces_on_surfaces[" // int2str(s) // "]/name", surface_integral_name)

        if(have_option(trim(complete_field_path(vfield%option_path, stat=stat)) // "/stat/compute_body_forces_on_surfaces[" // int2str(s) // "]/output_terms")) then
          if(have_viscosity) then
            ! calculate the forces on the surface
            call diagnostic_body_drag(state, force, surface_integral_name, pressure_force = pressure_force, viscous_force = viscous_force)
          else
            call diagnostic_body_drag(state, force, surface_integral_name, pressure_force = pressure_force)   
          end if
          if(getprocno() == 1) then
            do i=1, mesh_dim(vfield%mesh)
              write(default_stat%diag_unit, trim(format), advance="no") force(i)
            end do
            do i=1, mesh_dim(vfield%mesh)
              write(default_stat%diag_unit, trim(format), advance="no") pressure_force(i)
            end do
            if(have_viscosity) then
              do i=1, mesh_dim(vfield%mesh)
               write(default_stat%diag_unit, trim(format), advance="no") viscous_force(i)
              end do
            end if
          end if
        else
            ! calculate the forces on the surface
            call diagnostic_body_drag(state, force, surface_integral_name) 
            if(getprocno() == 1) then
             do i=1, mesh_dim(vfield%mesh)
                write(default_stat%diag_unit, trim(format), advance="no") force(i)
             end do
            end if     
        end if
      end do
      
    end subroutine write_body_forces

    subroutine write_momentum_conservation_error(state, v_field)
      type(state_type), intent(in) :: state
      type(vector_field), intent(inout) :: v_field
      
      type(vector_field), pointer :: velocity, old_velocity
      type(vector_field), pointer :: new_positions, nl_positions, old_positions
      type(scalar_field), pointer :: old_pressure, new_pressure
      type(scalar_field) :: nl_pressure, vel_comp
      real :: theta, dt
      integer :: sele, dim
      
      real, dimension(v_field%dim) :: momentum_cons, velocity_int, old_velocity_int, pressure_surface_int
      
      velocity => extract_vector_field(state, "Velocity")
      old_velocity => extract_vector_field(state, "OldVelocity")
      
      new_positions => extract_vector_field(state, "IteratedCoordinate")
      nl_positions => extract_vector_field(state, "Coordinate")
      old_positions => extract_vector_field(state, "OldCoordinate")
      
      new_pressure => extract_scalar_field(state, "Pressure")
      old_pressure => extract_scalar_field(state, "OldPressure")
      
      call get_option("/timestepping/timestep", dt)
      call get_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/theta", theta)
      
      call allocate(nl_pressure, new_pressure%mesh, "NonlinearPressure")
      call set(nl_pressure, new_pressure, old_pressure, theta)
      
      do dim = 1, velocity%dim
        vel_comp = extract_scalar_field(velocity, dim)
        call field_stats(vel_comp, new_positions, velocity_int(dim))

        vel_comp = extract_scalar_field(old_velocity, dim)
        call field_stats(vel_comp, old_positions, old_velocity_int(dim))
      end do
      
      ! pressure surface integral
      pressure_surface_int = 0.0
      do sele = 1, surface_element_count(v_field)
      
        pressure_surface_int = pressure_surface_int + pressure_surface_integral_face(sele, nl_pressure, nl_positions)
        
      end do
      
      ewrite(2,*) 'velocity_int = ', velocity_int
      ewrite(2,*) 'old_velocity_int = ', old_velocity_int
      ewrite(2,*) '(velocity_int-old_velocity_int)/dt = ', (velocity_int-old_velocity_int)/dt
      ewrite(2,*) 'pressure_surface_int = ', pressure_surface_int
      
      momentum_cons = (velocity_int-old_velocity_int)/dt - pressure_surface_int
      
      if(getprocno() == 1) then
        do dim=1, velocity%dim
          write(default_stat%diag_unit, trim(format), advance="no") momentum_cons(dim)
        end do
      end if

      call deallocate(nl_pressure)
     
    end subroutine write_momentum_conservation_error

    function pressure_surface_integral_face(sele, nl_pressure, nl_positions) result(pn)
    
      integer :: sele
      type(scalar_field) :: nl_pressure
      type(vector_field) :: nl_positions
      real, dimension(mesh_dim(nl_pressure)) :: pn
      
      real, dimension(mesh_dim(nl_pressure), face_ngi(nl_pressure, sele)) :: normal
      real, dimension(face_ngi(nl_pressure, sele)) :: detwei
      
      integer :: dim
    
      call transform_facet_to_physical( &
        nl_positions, sele, detwei_f = detwei, normal = normal)
        
      do dim = 1, size(normal, 1)
      
        pn(dim) = dot_product(face_val_at_quad(nl_pressure, sele), detwei*normal(dim, :))

      end do
      
    end function pressure_surface_integral_face

  end subroutine write_diagnostics

  subroutine test_and_write_convergence(state, time, dt, it, maxerror)
    !!< Test and write the diagnostics to the previously opened convergence file.

    type(state_type), dimension(:), intent(inout) :: state
    real, intent(in) :: time, dt
    integer, intent(in) :: it
    real, intent(out) :: maxerror

    character(len=10) :: format, iformat
    integer :: i, j, phase
    real :: error
    type(scalar_field) :: vfield_comp, nlvfield_comp
    type(scalar_field), pointer :: sfield, nlsfield
    type(vector_field), pointer :: vfield, nlvfield
    
    type(vector_field), pointer :: coordinates
    integer :: convergence_norm

    maxerror = 0.0

    format='(e15.6e3)'
    iformat='(i4)'

    if(default_stat%write_convergence_file) then
      ! Only the first process should write convergence information
      if(getprocno() == 1) then
         write(default_stat%conv_unit, format, advance="no") time
         write(default_stat%conv_unit, format, advance="no") dt
         write(default_stat%conv_unit, iformat, advance="no") it
      end if
    end if
    
    coordinates => extract_vector_field(state(1), "Coordinate")
    convergence_norm = convergence_norm_integer("/timestepping/nonlinear_iterations/tolerance")

    phaseloop: do phase=1,size(state)

       do i=1, size(default_stat%sfield_list(phase)%ptr)
          ! Output convergence information for each scalar field.
          sfield=>extract_scalar_field(state(phase), &
               &                       default_stat%sfield_list(phase)%ptr(i))

          if(.not. convergence_field(sfield)) then
            cycle
          end if

          nlsfield=>extract_scalar_field(state(phase), &
                                       "Iterated"//trim(default_stat%sfield_list(phase)%ptr(i)))

          call field_con_stats(sfield, nlsfield, error, &
                               convergence_norm, coordinates)
          maxerror = max(maxerror, error)

          if(default_stat%write_convergence_file) then
            ! Only the first process should write convergence information
            if(getprocno() == 1) then
               write(default_stat%conv_unit, format, advance="no") error
            end if
          end if

       end do

       do i = 1, size(default_stat%vfield_list(phase)%ptr)
         ! Output convergence information for each vector field

         vfield => extract_vector_field(state(phase), &
           & default_stat%vfield_list(phase)%ptr(i))

         if(.not. convergence_field(vfield)) then
           cycle
         end if

         nlvfield => extract_vector_field(state(phase), &
           & "Iterated"//default_stat%vfield_list(phase)%ptr(i))

         call field_con_stats(vfield, nlvfield, error, &
                              convergence_norm, coordinates)
         maxerror = max(maxerror, error)

         if(default_stat%write_convergence_file) then
            ! Only the first process should write convergence information
            if(getprocno() == 1) then
            write(default_stat%conv_unit, format, advance = "no") error
            end if
         end if

         if(.not. convergence_field(vfield, test_for_components = .true.)) then
           cycle
         end if

         do j = 1, vfield%dim
           vfield_comp = extract_scalar_field(vfield, j)
           nlvfield_comp = extract_scalar_field(nlvfield, j)

           call field_con_stats(vfield_comp, nlvfield_comp, error, &
                                convergence_norm, coordinates)
           maxerror = max(maxerror, error)

           if(default_stat%write_convergence_file) then
               ! Only the first process should write convergence information
               if(getprocno() == 1) then
                  write(default_stat%conv_unit, format, advance = "no") error
               end if
           end if
         end do
       end do

    end do phaseloop

    if(default_stat%write_convergence_file) then
      ! Output end of line
      ! Only the first process should write convergence information
      if(getprocno() == 1) then
         write(default_stat%conv_unit,'(a)') ""
      end if
    end if

    if(have_option("/io/convergence/convergence_vtus")) then
      call vtk_write_state_new_options(filename="convergence_test", index=it, state=state)
    end if

  end subroutine test_and_write_convergence

  subroutine test_and_write_steady_state(state, maxchange)
    !!< Test whether a steady state has been reached.

    type(state_type), dimension(:), intent(in) :: state
    real, intent(out) :: maxchange

    integer :: i, j, phase
    real :: change, dt
    type(scalar_field) :: vfield_comp, oldvfield_comp
    type(scalar_field), pointer :: sfield, oldsfield
    type(vector_field), pointer :: vfield, oldvfield

    logical :: acceleration
    
    type(vector_field), pointer :: coordinates
    integer :: convergence_norm

    character(len = *), parameter :: format = "(e15.6e3)"
    integer :: procno
    real :: elapsed_time

    ewrite(1, *) "Entering test_and_write_steady_state"

    maxchange = 0.0

    acceleration = have_option("/timestepping/steady_state/acceleration_form")
    call get_option("/timestepping/timestep", dt)

    coordinates => extract_vector_field(state(1), "Coordinate")
    convergence_norm = convergence_norm_integer("/timestepping/steady_state/tolerance")

    procno = getprocno()
    if(default_stat%write_steady_state_file .and. procno == 1) then    
      call get_option("/timestepping/current_time", elapsed_time)
      if(default_stat%binary_steady_state_output) then
        write(default_stat%steady_state_unit) elapsed_time
        write(default_stat%steady_state_unit) dt
      else
        write(default_stat%steady_state_unit, format, advance="no") elapsed_time
        write(default_stat%steady_state_unit, format, advance="no") dt
      end if
    end if

    phaseloop: do phase=1,size(state)

       do i=1, size(default_stat%sfield_list(phase)%ptr)
          ! Test steady state information for each scalar field.

          sfield=>extract_scalar_field(state(phase), i)
          if(.not. steady_state_field(sfield)) cycle
          ! Scalar fields

          oldsfield=>extract_scalar_field(state(phase), &
                                       "Old"//trim(default_stat%sfield_list(phase)%ptr(i)))

          call field_con_stats(sfield, oldsfield, change, &
                               convergence_norm, coordinates)
          if(acceleration) change = change/dt
          ewrite(2, *) trim(state(phase)%name)//"::"//trim(sfield%name), change
          maxchange = max(maxchange, change)

          if(default_stat%write_steady_state_file .and. procno == 1) then
            if(default_stat%binary_steady_state_output) then
              write(default_stat%steady_state_unit) change
            else
              write(default_stat%steady_state_unit, format, advance = "no") change
            end if
          end if

       end do

       do i = 1, vector_field_count(state(phase))
         vfield => extract_vector_field(state(phase), i)
         if(.not. steady_state_field(vfield)) cycle
         ! Vector fields


         oldvfield => extract_vector_field(state(phase), &
           & "Old"//default_stat%vfield_list(phase)%ptr(i))

         call field_con_stats(vfield, oldvfield, change, &
                              convergence_norm, coordinates)
         if(acceleration) change = change/dt
         ewrite(2, *) trim(state(phase)%name)//"::"//trim(vfield%name), change
         maxchange = max(maxchange, change)

         if(default_stat%write_steady_state_file .and. procno == 1) then
            if(default_stat%binary_steady_state_output) then
              write(default_stat%steady_state_unit) change
            else
              write(default_stat%steady_state_unit, format, advance = "no") change
            end if
         end if

         if(.not. steady_state_field(vfield, test_for_components = .true.)) cycle
         ! Vector field components

         do j = 1, vfield%dim
           vfield_comp = extract_scalar_field(vfield, j)
           oldvfield_comp = extract_scalar_field(oldvfield, j)

           call field_con_stats(vfield_comp, oldvfield_comp, change, &
                                convergence_norm, coordinates)
           if(acceleration) change = change/dt
           ewrite(2, *) trim(state(phase)%name)//"::"//trim(vfield%name), j,  change
           maxchange = max(maxchange, change)

           if(default_stat%write_steady_state_file .and. procno == 1) then
             if(default_stat%binary_steady_state_output) then
               write(default_stat%steady_state_unit) change
             else
               write(default_stat%steady_state_unit, format, advance = "no") change
             end if
           end if
         end do

       end do

    end do phaseloop

    ewrite(1, *) "maxchange = ", maxchange
    
    if(default_stat%write_steady_state_file .and. procno == 1) then
      if(default_stat%binary_steady_state_output) then
        write(default_stat%steady_state_unit) maxchange
      else
        write(default_stat%steady_state_unit, format, advance = "no") maxchange      
        ! Output end of line
        write(default_stat%steady_state_unit,'(a)') ""
      end if
      
      flush(default_stat%steady_state_unit)
    end if

    ewrite(1, *) "Exiting test_and_write_steady_state"

  end subroutine test_and_write_steady_state

  subroutine write_detectors(state, detector_list, time, dt)
    !!< Write the field values at detectors to the previously opened detectors file.
    type(state_type), dimension(:), intent(in) :: state
    type(detector_linked_list), intent(inout) :: detector_list
    real, intent(in) :: time, dt

    character(len=10) :: format_buffer
    integer :: i, j, k, phase, ele, check_no_det, totaldet_global
    real :: value
    real, dimension(:), allocatable :: vvalue
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(detector_type), pointer :: detector

    ewrite(1,*) "In write_detectors"

    !Computing the global number of detectors. This is to prevent hanging
    !when there are no detectors on any processor
    check_no_det=1
    if (detector_list%length==0) then
       check_no_det=0
    end if
    call allmax(check_no_det)
    if (check_no_det==0) then
       return
    end if

    ! If isparallel() or binary output use this:
    if ((isparallel()).or.(default_stat%detector_list%binary_output)) then    
       call write_mpi_out(state,detector_list,time,dt)
    else ! This is only for single processor with ascii output
       if(getprocno() == 1) then
          if(detector_list%binary_output) then
             write(detector_list%output_unit) time
             write(detector_list%output_unit) dt
          else
             format_buffer=reals_format(1)
             write(detector_list%output_unit, format_buffer, advance="no") time
             write(detector_list%output_unit, format_buffer, advance="no") dt
          end if
       end if

       ! Next columns contain the positions of all the detectors.
       detector => detector_list%first
       positionloop: do i=1, detector_list%length
          if(detector_list%binary_output) then
             write(detector_list%output_unit) detector%position
          else
             format_buffer=reals_format(size(detector%position))
             write(detector_list%output_unit, format_buffer, advance="no") &
                  detector%position
          end if

          detector => detector%next
       end do positionloop

       phaseloop: do phase=1,size(state)
          if (size(detector_list%sfield_list(phase)%ptr)>0) then
             do i=1, size(detector_list%sfield_list(phase)%ptr)
                ! Output statistics for each scalar field.
                sfield=>extract_scalar_field(state(phase), detector_list%sfield_list(phase)%ptr(i))

                detector => detector_list%first
                do j=1, detector_list%length
                   if (detector%element<0) then
                      if (detector_list%write_nan_outside) then
                         call cget_nan(value)
                      else
                         FLExit("Trying to write detector that is outside of domain.")
                      end if
                   else
                      value = detector_value(sfield, detector)
                   end if

                   if(detector_list%binary_output) then
                      write(detector_list%output_unit) value
                   else
                      format_buffer=reals_format(1)
                      write(detector_list%output_unit, format_buffer, advance="no") value
                   end if
                   detector => detector%next
                end do
             end do
          end if

          if (size(detector_list%vfield_list(phase)%ptr)>0) then
             do i = 1, size(detector_list%vfield_list(phase)%ptr)
                ! Output statistics for each vector field
                vfield => extract_vector_field(state(phase), &
                  & detector_list%vfield_list(phase)%ptr(i))
                allocate(vvalue(vfield%dim))

                detector => detector_list%first
                do j=1, detector_list%length
                   if (detector%element<0) then
                      if (detector_list%write_nan_outside) then
                         call cget_nan(value)
                         vvalue(:) = value
                      else
                         FLExit("Trying to write detector that is outside of domain.")
                      end if
                   else
                      vvalue = detector_value(vfield, detector)
                   end if

                   if(detector_list%binary_output) then
                      write(detector_list%output_unit) vvalue
                   else
                      format_buffer=reals_format(vfield%dim)
                      write(detector_list%output_unit, format_buffer, advance="no") vvalue
                   end if
                   detector => detector%next
                end do
                deallocate(vvalue)
             end do
          end if

       end do phaseloop

       ! Output end of line
       if (.not. detector_list%binary_output) then
          ! Output end of line
          write(detector_list%output_unit,'(a)') ""
       end if
       flush(detector_list%output_unit)
    end if

    totaldet_global=detector_list%length
    call allsum(totaldet_global)
    ewrite(2,*) "Found", detector_list%length, "local and", totaldet_global, "global detectors"

    if (totaldet_global/=detector_list%total_num_det) then
       ewrite(2,*) "We have either duplication or have lost some det"
       ewrite(2,*) "totaldet_global", totaldet_global
       ewrite(2,*) "total_num_det", detector_list%total_num_det
    end if

    ewrite(1,*) "Exiting write_detectors"

  contains

    function reals_format(reals)
      character(len=10) :: reals_format
      integer :: reals

      write(reals_format, '(a,i0,a)') '(',reals,'e15.6e3)'

    end function reals_format

  end subroutine write_detectors

  subroutine write_mpi_out(state,detector_list,time,dt)
    !!< Writes detector information (position, value of scalar and vector fields at that position, etc.) into detectors file using MPI output 
    ! commands so that when running in parallel all processors can write at the same time information into the file at the right location.       

    type(state_type), dimension(:), intent(in) :: state
    type(detector_linked_list), intent(inout) :: detector_list
    real, intent(in) :: time, dt

    integer :: i, j, phase, ierror, number_of_scalar_det_fields, realsize, dim, procno
    integer(KIND = MPI_OFFSET_KIND) :: location_to_write, offset
    integer :: number_of_vector_det_fields, number_total_columns

    real :: value
    real, dimension(:), allocatable :: vvalue
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(detector_type), pointer :: node

    ewrite(2, *) "In write_mpi_out"

    detector_list%mpi_write_count = detector_list%mpi_write_count + 1
    ewrite(2, *) "Writing detector output ", detector_list%mpi_write_count
    
    procno = getprocno()

    ewrite(2, *) "Number of detector scalar fields = ", detector_list%num_sfields
    ewrite(2, *) "Number of detector vector fields = ", detector_list%num_vfields

    call mpi_type_extent(getpreal(), realsize, ierror)
    assert(ierror == MPI_SUCCESS)

    vfield => extract_vector_field(state, "Coordinate")
    dim = vfield%dim
                           ! Time data
    number_total_columns = 2 + &
                           ! Detector coordinates
                         & detector_list%total_num_det * dim + &
                           ! Scalar detector data
                         & detector_list%total_num_det * detector_list%num_sfields + &
                           ! Vector detector data
                         & detector_list%total_num_det * detector_list%num_vfields * dim

    ! raise kind of one of the variables (each individually is a 4 byte-integer) such that the calculation is coerced to be of MPI_OFFSET_KIND (typically 8 bytes)
    ! this is necessary for files bigger than 2GB
    location_to_write = (int(detector_list%mpi_write_count, kind=MPI_OFFSET_KIND) - 1) * number_total_columns * realsize

    if(procno == 1) then
      ! Output time data
      call mpi_file_write_at(detector_list%mpi_fh, location_to_write, time, 1, getpreal(), MPI_STATUS_IGNORE, ierror)
      assert(ierror == MPI_SUCCESS)
        
      call mpi_file_write_at(detector_list%mpi_fh, location_to_write + realsize, dt, 1, getpreal(), MPI_STATUS_IGNORE, ierror)
      assert(ierror == MPI_SUCCESS)
    end if
    location_to_write = location_to_write + 2 * realsize

    node => detector_list%first
    position_loop: do i = 1, detector_list%length
      ! Output detector coordinates
      assert(size(node%position) == dim)  
    
      offset = location_to_write + (node%id_number - 1) * dim * realsize

      call mpi_file_write_at(detector_list%mpi_fh, offset, node%position, dim, getpreal(), MPI_STATUS_IGNORE, ierror)
      assert(ierror == MPI_SUCCESS)
      node => node%next
    end do position_loop
    assert(.not. associated(node))
    location_to_write = location_to_write + detector_list%total_num_det * dim * realsize

    allocate(vvalue(dim))
    state_loop: do phase = 1, size(state)

      if (allocated(detector_list%sfield_list)) then
        if (size(detector_list%sfield_list(phase)%ptr)>0) then
        scalar_loop: do i = 1, size(detector_list%sfield_list(phase)%ptr)
          ! Output statistics for each scalar field        
          sfield => extract_scalar_field(state(phase), detector_list%sfield_list(phase)%ptr(i))

          node => detector_list%first
          scalar_node_loop: do j = 1, detector_list%length
            if (node%element<0) then
               if (detector_list%write_nan_outside) then
                  call cget_nan(value)
               else
                  FLExit("Trying to write detector that is outside of domain.")
               end if
            else
               value = detector_value(sfield, node)
            end if

            offset = location_to_write + (detector_list%total_num_det * (i - 1) + (node%id_number - 1)) * realsize

            call mpi_file_write_at(detector_list%mpi_fh, offset, value, 1, getpreal(), MPI_STATUS_IGNORE, ierror)
            assert(ierror == MPI_SUCCESS)
            node => node%next
          end do scalar_node_loop
          assert(.not. associated(node))
        end do scalar_loop
        end if
      end if
      location_to_write = location_to_write + detector_list%total_num_det *detector_list%num_sfields * realsize
      
      if (allocated(detector_list%vfield_list)) then
        if (size(detector_list%vfield_list(phase)%ptr)>0) then
        vector_loop: do i = 1, size(detector_list%vfield_list(phase)%ptr)
          ! Output statistics for each vector field     
          vfield => extract_vector_field(state(phase), detector_list%vfield_list(phase)%ptr(i))

          ! We currently don't have enough information for mixed dimension
          ! vector fields (see below)
          assert(vfield%dim == dim)

          node => detector_list%first
          vector_node_loop: do j = 1, detector_list%length
            if (node%element<0) then
               if (detector_list%write_nan_outside) then
                  call cget_nan(value)
                  vvalue(:) = value
               else
                  FLExit("Trying to write detector that is outside of domain.")
               end if
            else
               vvalue = detector_value(vfield, node)
            end if

            ! Currently have to assume single dimension vector fields in
            ! order to compute the offset
            offset = location_to_write + (detector_list%total_num_det * (i - 1) + (node%id_number - 1)) * dim * realsize

            call mpi_file_write_at(detector_list%mpi_fh, offset, vvalue, dim, getpreal(), MPI_STATUS_IGNORE, ierror)
            assert(ierror == MPI_SUCCESS)
            node => node%next
          end do vector_node_loop
          assert(.not. associated(node))
        end do vector_loop
        end if
      end if
      location_to_write = location_to_write + detector_list%total_num_det * detector_list%num_vfields * dim * realsize
        
    end do state_loop
    deallocate(vvalue)

    call mpi_file_sync(detector_list%mpi_fh, ierror)
    assert(ierror == MPI_SUCCESS)

!    ! The following was used when debugging to check some of the data written
!    ! into the file
!    ! Left here in case someone would like to use the mpi_file_read_at for
!    ! debugging or checking
!   
!    number_total_columns = 2 + total_num_det * dim
!    allocate(buffer(number_total_columns))
!
!    call mpi_file_read_at(fh, 0, buffer, size(buffer), getpreal(), status, ierror)
!    call mpi_get_count(status, getpreal(), count,  ierror)
!    assert(ierror == MPI_SUCCESS)
!
!    call mpi_barrier(MPI_COMM_FEMTOOLS, ierror)
!    assert(ierror == MPI_SUCCESS)
!    
!    deallocate(buffer)
!    ewrite(2, "(a,i0,a)") "Read ", count, " reals"

    ewrite(2, *) "Exiting write_mpi_out"
   
  end subroutine write_mpi_out

  subroutine list_det_into_csr_sparsity(detector_list,ihash_sparsity,list_into_array,element_detector_list,count)
!! This subroutine creates a hash table called ihash_sparsity and a csr_sparsity matrix called element_detector_list that we use to find out 
!! how many detectors a given element has and we also obtain the location (row index) of those detectors in an array called list_into_array. 
!! This array contains the information of detector_list but in an array format, each row of the array contains the information of a detector. 
!! By accessing the array at that/those row indexes we can extract the information (position, id_number, type, etc.) of each detector present 
!! in the element (each row index corresponds to a detector)

    type(detector_linked_list), intent(inout) :: detector_list
    type(integer_hash_table), intent(inout) :: ihash_sparsity
    type(csr_sparsity), intent(inout) :: element_detector_list
    real, dimension(:,:), allocatable, intent(inout) :: list_into_array
    integer, intent(in) :: count
    
    integer, dimension(:), allocatable:: detector_count ! detectors per element
    integer, dimension(:), pointer:: detectors
    type(detector_type), pointer :: node
    type(vector_field), pointer :: vfield, xfield
    integer :: dim, i, ele, no_rows, entries, row, pos

    if (detector_list%length/=0) then
   
       node => detector_list%first

       dim=size(node%position)

       do i=1, detector_list%length

             list_into_array(i,1:dim)=node%position
             list_into_array(i,dim+1)=node%element
             list_into_array(i,dim+2)=node%id_number
             list_into_array(i,dim+3)=node%type
             list_into_array(i,dim+4)=0.0

             node => node%next
    
       end do

       ! create map between element and rows, where each row corresponds to an element 
       ! with one or more detectors
    
       ! loop over detectors: 
       ! detector is in element ele

       no_rows=count

      ! count number of detectors per row
      allocate(detector_count(1:no_rows))
      detector_count=0
      ! loop over detectors:
      node => detector_list%first

      do i=1, detector_list%length

         ele=node%element
         if  (has_key(ihash_sparsity, ele)) then
            row=fetch(ihash_sparsity, ele)
            detector_count(row)=detector_count(row)+1
         end if
         node => node%next

      end do 

      ! set up %findrm, the beginning of each row in memory
      pos=1 ! position in colm
      do row=1, no_rows
         element_detector_list%findrm(row)=pos
         pos=pos+detector_count(row)
      end do
      element_detector_list%findrm(row)=pos

      ! fill up the rows with the rom_number of the detectors in the list_into_array
      detector_count=0
      ! loop over detectors:

      do i=1, detector_list%length
      
         ele=list_into_array(i,dim+1)  
         if (has_key(ihash_sparsity, ele)) then
            row=fetch(ihash_sparsity, ele)
            detectors => row_m_ptr(element_detector_list, row)
            detector_count(row)=detector_count(row)+1
            detectors(detector_count(row))=i
         end if

      end do

      deallocate(detector_count)

    end if

  end subroutine list_det_into_csr_sparsity
    
  subroutine close_diagnostic_files()
    !! Closes .stat, .convergence and .detector file (if openened)
    !! Gives a warning for iostat/=0, no point to flabort though.

    integer:: stat, IERROR

    if (default_stat%diag_unit/=0) then
       close(default_stat%diag_unit, iostat=stat)
       if (stat/=0) then
          ewrite(0,*) "Warning: failed to close .stat file"
       end if
    end if

    if (default_stat%conv_unit/=0) then
       close(default_stat%conv_unit, iostat=stat)
       if (stat/=0) then
          ewrite(0,*) "Warning: failed to close .convergence file"
       end if
    end if

    if(default_stat%steady_state_unit /= 0) then
      close(default_stat%steady_state_unit, iostat = stat)
      if(stat /= 0) then
        ewrite(0, *) "Warning: failed to close .steady_state file"
      end if
    end if

    if (default_stat%detector_list%output_unit/=0) then
       close(default_stat%detector_list%output_unit, iostat=stat)
       if (stat/=0) then
          ewrite(0,*) "Warning: failed to close .detector file"
       end if
    end if

    if (default_stat%detector_list%mpi_fh/=0) then
       call MPI_FILE_CLOSE(default_stat%detector_list%mpi_fh, IERROR) 
       if(IERROR /= MPI_SUCCESS) then
          ewrite(0,*) "Warning: failed to close .detector file open with mpi_file_open"
       end if
    end if

  end subroutine close_diagnostic_files

  SUBROUTINE RUN_DIAGNOSTICS(state)
    !!< Initial diagnostic output.
    type(state_type), dimension(:), intent(in) :: state

    REAL ::DT,LTIME

    REAL :: VOL,MAXVOL,MINVOL
    INTEGER :: ELE, I, minvol_ele, maxvol_ele
    
    real, dimension(:), allocatable :: detwei
    type(vector_field) :: coordinate

    ! Only do this at all if there will be output.
    if (debug_level()<1) return 
    
    call get_option("/timestepping/timestep", dt)
    call get_option("/timestepping/finish_time", ltime)

    ewrite(1,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    ewrite(1,*)'% Some quantities associated with the initial set-up of this problem. %'
    ewrite(1,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    ewrite(1,*)'-'
    ewrite(1,*)'The time step (DT) is:                                 ',DT
    ewrite(1,*)'The end time (LTIME) is set to:                        ',LTIME
    ewrite(1,*)'This corresponds to this many time steps in simulation:',LTIME/DT 
    ewrite(1,*)'-'

    coordinate=extract_vector_field(state(1), "Coordinate")

    ! Edge lengths are suspended until someone generalises them
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !                    ELEMENT VOLUMES      
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAXVOL = -huge(1.0)
    MINVOL =  huge(1.0)

    allocate(detwei(ele_ngi(coordinate,1)))

    do ele = 1, element_count(coordinate)
       
       call transform_to_physical(coordinate, ele, detwei)
       vol=sum(detwei)

       if (vol>maxvol) then
          maxvol=vol
          maxvol_ele=ele
       end if
       if (vol<minvol) then
          minvol=vol
          minvol_ele=ele
       end if
       
    end do
    
    call allmax(maxvol)
    call allmin(minvol)
    ewrite(1,*)'The maximum volume of element in the mesh is:',maxvol
    ewrite(1,*)'The minimum volume of element in the mesh is:',minvol
    ewrite(1,*)'The ratio of max to min volumes is:          ',maxvol/minvol
            
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !                    Fields you've given     
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ewrite(1,*)'The following material phases and their ', &
      'contents are in existence'
    do i = 1, size(state)
       call print_state(state(i), unit=debug_unit(1))
    end do

  END SUBROUTINE RUN_DIAGNOSTICS
  
  subroutine diagnostic_variables_check_options
  
#ifndef STREAM_IO
    if(have_option("/io/detectors/binary_output")) then
      FLExit("Cannot use binary detector output format - no stream I/O support")
    end if
    
    if(have_option("/timestepping/steady_state/steady_state_file/binary_output")) then
      FLExit("Cannot use binary steady state output format - no stream I/O support")
    end if
#endif

  end subroutine diagnostic_variables_check_options

end module diagnostic_variables
