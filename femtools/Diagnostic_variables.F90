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

module diagnostic_variables
  !!< A module to calculate and output diagnostics. This replaces the .s file.
  use quadrature
  use elements
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, &
    & PYTHON_FUNC_LEN, int_16, integer_size, real_size
  use fields
  use fields_base
  use field_derivatives
  use field_options
  use state_module
  use futils
  use fetools
  use fefields
  use MeshDiagnostics
  use spud
  use parallel_tools
  use Profiler
  use sparsity_patterns
  use solvers
  use write_state_module, only: vtk_write_state_new_options
  use surface_integrals
  use vtk_interfaces
  use embed_python
  use eventcounter
  use pickers
  use sparse_tools
  use mixing_statistics
  use c_interfaces
!  use checkpoint
  use memory_diagnostics
  use data_structures
  use unittest_tools
  use integer_hash_table_module
  use halo_data_types
  use halos_allocates
  use halos_base
  use halos_debug
  use halos_numbering
  use halos_ownership
  use halos_derivation
  use halos_communications
  use halos_registration
  use mpi_interfaces
  use parallel_tools
  use fields_manipulation
  use detector_data_types
  use ieee_arithmetic, only: cget_nan
  use state_fields_module, only: get_cv_mass
  use diagnostic_tools
  
  implicit none

  private

  public :: initialise_diagnostics, initialise_convergence, &
       & initialise_steady_state, field_tag, write_diagnostics, &
       & test_and_write_convergence, &
       & test_and_write_steady_state, steady_state_field, convergence_field, &
       & close_diagnostic_files, run_diagnostics, &
       & diagnostic_variables_check_options, &
       & initialise_walltime, &
       & uninitialise_diagnostics, register_diagnostic, destroy_registered_diagnostics, set_diagnostic, &
       & get_diagnostic, initialise_constant_diagnostics

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
    integer(kind = int_16) :: elapsed_count
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
    integer :: no_mixing_bins, n_stat_fields
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
          n_stat_fields = 0
          do i=1, size(state(phase)%scalar_names)
             if (stat_field(state(phase)%scalar_fields(i)%ptr, state(phase))) then
                n_stat_fields = n_stat_fields + 1
             end if
          end do
          allocate(default_stat%sfield_list(phase)%ptr(n_stat_fields))
          n_stat_fields = 0
          do i=1, size(state(phase)%scalar_names)
             if (stat_field(state(phase)%scalar_fields(i)%ptr, state(phase))) then
                n_stat_fields = n_stat_fields + 1
                default_stat%sfield_list(phase)%ptr(n_stat_fields) = state(phase)%scalar_names(i)
             end if
          end do
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
  
  subroutine write_diagnostics(state, time, dt, timestep)
      !!< Write the diagnostics to the previously opened diagnostics file.
    type(state_type), dimension(:), intent(inout) :: state
    real, intent(in) :: time, dt
    integer, intent(in) :: timestep

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

    ewrite(1,*) 'In write_diagnostics'
    call profiler_tic("I/O")

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
          call field_stats(sfield, Xfield, fmin, fmax, fnorm2, fintegral)
          if(getprocno() == 1) then
            write(default_stat%diag_unit, trim(format4), advance="no") fmin, fmax, fnorm2,&
                 & fintegral
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
       if (default_stat%detector_list%binary_output) then
          call MPI_FILE_CLOSE(default_stat%detector_list%output_unit, IERROR) 
          if (IERROR /= MPI_SUCCESS) then
             ewrite(0,*) "Warning: failed to close .detector file open with mpi_file_open"
          end if
       else
          close(default_stat%detector_list%output_unit, iostat=stat)
          if (stat/=0) then
             ewrite(0,*) "Warning: failed to close .detector file"
          end if
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
