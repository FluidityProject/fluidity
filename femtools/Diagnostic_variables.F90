!    Copyright (C) 2006 Imperial College London and others.
   
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

module diagnostic_variables
  !!< A module to calculate and output diagnostics. This replaces the .s file.
  use quadrature
  use elements
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, &
    & PYTHON_FUNC_LEN, int_16, integer_size, real_size
  use fields
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
  use detector_data_types
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

  implicit none

  private

  public :: initialise_diagnostics, initialise_convergence, &
       & initialise_steady_state, field_tag, write_diagnostics, &
       & test_and_write_convergence, initialise_detectors, write_detectors, &
       & test_and_write_steady_state, steady_state_field, convergence_field, &
       & close_diagnostic_files, run_diagnostics, &
       & diagnostic_variables_check_options

  public ::  detector_list, name_of_detector_groups_in_read_order, number_det_in_each_group, name_of_detector_in_read_order

  interface stat_field
    module procedure stat_field_scalar, stat_field_vector
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

  interface detector_value
     module procedure detector_value_scalar, detector_value_vector
  end interface

  interface insert_det
     module procedure detector_list_insert
  end interface

  !! Output unit for diagnostics file.
  !! (assumed non-opened as long these are 0)
  integer, save :: diag_unit=0, conv_unit=0, &
    & detector_unit=0, detector_checkpoint_unit=0, detector_file_unit=0 
  logical, save :: binary_detector_output = .false.

  integer, save :: fh=0, total_num_det

  !! Are we writing to a convergence file?
  logical, save :: write_convergence_file=.false.

  !! Output unit for .steady_state file (assumed non-opened as long as == 0)
  integer, save :: steady_state_unit = 0
  !! Are we writing to a steady state file?
  logical, save :: write_steady_state_file = .false.
  logical, save :: binary_steady_state_output = .false.

  !! Are we continuing from a detector checkpoint file?
  logical, save :: detectors_checkpoint_done = .false.

  type stringlist
     !!< Container type for a list of strings.
     character(len=FIELD_NAME_LEN), dimension(:), pointer :: ptr
  end type stringlist

  character(len = FIELD_NAME_LEN), dimension(:), allocatable, save :: mesh_list
  !! List of scalar fields to output. This is stored to ensure that
  !! additional fields inserted into state during running do not bugger up
  !! the output.
  type(stringlist), dimension(:), allocatable, save :: sfield_list
  type(stringlist), dimension(:), allocatable, save :: vfield_list

  !! Similar lists for detectors
  type(stringlist), dimension(:), allocatable, save :: detector_sfield_list
  type(stringlist), dimension(:), allocatable, save :: detector_vfield_list

  character(len = FIELD_NAME_LEN), dimension(:), allocatable, save :: name_of_detector_groups_in_read_order, name_of_detector_in_read_order
  integer, dimension(:), allocatable, save :: number_det_in_each_group

! type(detector_type), dimension(:), allocatable, target, save :: detector_list

  type(detector_linked_list), target, save :: detector_list
!  type(detector_linked_list), target, save :: send_list
  real, dimension(:,:), allocatable, save :: array_det_info

  !! Recording wall time since the system start
  integer, save :: current_count, count_rate, count_max
  integer(kind = int_16), save :: elapsed_count

contains

  function stat_mesh(mesh)
    !!< Return whether the supplied mesh should be included in the .stat file

    type(mesh_type), intent(in) :: mesh

    logical :: stat_mesh

    stat_mesh = have_option(trim(complete_mesh_path(mesh%option_path)) // "/stat/include_in_stat") &
      & .and. .not. have_option(trim(complete_mesh_path(mesh%option_path)) // "/stat/exclude_from_stat")

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

    call system_clock(current_count, count_rate, count_max)
    elapsed_count=0

  end subroutine initialise_walltime

  function elapsed_walltime()
    !!< Return the number of walltime seconds since the beginning of the
    !!< simulation.
    real :: elapsed_walltime
    
    integer :: new_count

    call system_clock(new_count)

    ! Deal with clock rollover. If one timestep takes more than a whole
    ! clock rollover, we have more problems than we can deal with!
    if (new_count<current_count) then
       elapsed_count=elapsed_count+(new_count-current_count)+count_max
    else
       elapsed_count=elapsed_count+(new_count-current_count)
    end if
    current_count=new_count

    elapsed_walltime=real(elapsed_count)/real(count_rate)

  end function elapsed_walltime

  subroutine initialise_diagnostics(filename, state)
    !!< Set up the diagnostic file headers.

    character(len=*) :: filename
    type(state_type), dimension(:), intent(in) :: state

    ! Idempotency variable
    logical, save :: initialised=.false.

    integer :: column, i, j, phase, stat
    integer, dimension(2) :: shape_option
    character(len = 254) :: buffer, material_phase_name
    character(len = FIELD_NAME_LEN) :: surface_integral_name, mixing_stats_name
    type(scalar_field) :: vfield_comp
    type(mesh_type), pointer :: mesh
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield

    ewrite(1, *) "In initialise_diagnostics"

    ! Idempotency check
    if(initialised) then
      ewrite(2, *) "Diagnostics already initialised"
      ewrite(1, *) "Exiting initialise_diagnostics"
      return
    end if
    initialised=.true.
    
    ! All processes must assemble the mesh and field lists
    ! Mesh field list
    allocate(mesh_list(size(state(1)%mesh_names)))
    mesh_list = state(1)%mesh_names
    ! Scalar field list
    allocate (sfield_list(size(state)))
    do phase=1, size(state)
       if (associated(state(phase)%scalar_names)) then
          allocate(sfield_list(phase)%ptr(size(state(phase)%scalar_names)))
          sfield_list(phase)%ptr=state(phase)%scalar_names
       else
          allocate(sfield_list(phase)%ptr(0))
       end if
    end do
    ! Vector field list
    allocate (vfield_list(size(state)))
    do phase = 1, size(state)
       if (associated(state(phase)%vector_names)) then
          allocate(vfield_list(phase)%ptr(size(state(phase)%vector_names)))
          vfield_list(phase)%ptr = state(phase)%vector_names
       else
          allocate(vfield_list(phase)%ptr(0))
       end if
    end do

    ! Only the first process should write statistics information (and hence
    ! write the headers)    
    if(getprocno() == 1) then
      diag_unit=free_unit()
      open(unit=diag_unit, file=trim(filename)//'.stat', action="write")

      write(diag_unit, '(a)') "<header>"

      call initialise_constant_diagnostics(diag_unit)

      column=0

      ! Initial columns are elapsed time and dt.
      column=column+1
      buffer=field_tag(name="ElapsedTime", column=column, statistic="value")
      write(diag_unit, '(a)') trim(buffer)
      column=column+1
      buffer=field_tag(name="dt", column=column, statistic="value")
      write(diag_unit, '(a)') trim(buffer)
      column=column+1
      buffer=field_tag(name="ElapsedWallTime", column=column, statistic="value")
      write(diag_unit, '(a)') trim(buffer)
      call initialise_walltime

      do i = 1, size(mesh_list)
        ! Headers for output statistics for each mesh
        mesh => extract_mesh(state(1), mesh_list(i))

        if(stat_mesh(mesh)) then
          column = column + 1
          buffer = field_tag(name = mesh%name, column = column, statistic = "nodes")
          write(diag_unit, "(a)"), trim(buffer)
          column = column + 1
          buffer = field_tag(name = mesh%name, column = column, statistic = "elements")
          write(diag_unit, "(a)"), trim(buffer)
          column = column + 1
          buffer = field_tag(name = mesh%name, column = column, statistic = "surface_elements")
          write(diag_unit, "(a)"), trim(buffer)
        end if
      end do

#ifdef HAVE_MEMORY_STATS
      ! Memory statistics
      do i=0, MEMORY_TYPES
          column = column + 1
          buffer = field_tag(name = memory_type_names(i), column = column,&
               & statistic = "current", material_phase_name="Memory")
          write(diag_unit, "(a)"), trim(buffer)         
          column = column + 1
          buffer = field_tag(name = memory_type_names(i), column = column,&
               & statistic = "min", material_phase_name="Memory")
          write(diag_unit, "(a)"), trim(buffer)         
          column = column + 1
          buffer = field_tag(name = memory_type_names(i), column = column,&
               & statistic = "max", material_phase_name="Memory")
          write(diag_unit, "(a)"), trim(buffer)         
      end do
#endif
      
      phaseloop: do phase=1,size(state)

         material_phase_name=trim(state(phase)%name)

         do i=1, size(sfield_list(phase)%ptr)
            ! Headers for output statistics for each scalar field
            sfield => extract_scalar_field(state(phase), sfield_list(phase)%ptr(i))

            ! Standard scalar field stats
            if(stat_field(sfield, state(phase))) then
              column=column+1
              buffer=field_tag(name=sfield%name, column=column, statistic="min", material_phase_name=material_phase_name)
              write(diag_unit, '(a)') trim(buffer)
              column=column+1
              buffer=field_tag(name=sfield%name, column=column, statistic="max", material_phase_name=material_phase_name)
              write(diag_unit, '(a)') trim(buffer)
              column=column+1
              buffer=field_tag(name=sfield%name, column=column, statistic="l2norm", material_phase_name=material_phase_name)
              write(diag_unit, '(a)') trim(buffer)
              column=column+1
              buffer=field_tag(name=sfield%name, column=column, statistic="integral", material_phase_name=material_phase_name)
              write(diag_unit, '(a)') trim(buffer)
            end if
            
            ! Control volume stats
            if(have_option(trim(complete_field_path(sfield%option_path, stat=stat)) // "/stat/include_cv_stats")) then
              column=column+1
              buffer=field_tag(name=sfield%name, column=column, statistic="cv_l2norm", material_phase_name=material_phase_name)
              write(diag_unit, '(a)') trim(buffer)
              column=column+1
              buffer=field_tag(name=sfield%name, column=column, statistic="cv_integral", material_phase_name=material_phase_name)
              write(diag_unit, '(a)') trim(buffer)
            end if
            
            ! Mixing stats
            do j = 0, option_count(trim(complete_field_path(sfield%option_path, stat=stat)) // "/stat/include_mixing_stats") - 1
              call get_option(trim(complete_field_path(sfield%option_path)) &
              & // "/stat/include_mixing_stats["// int2str(j) // "]/name", mixing_stats_name)
              shape_option=option_shape(trim(complete_field_path(sfield%option_path)) // &
                   & "/stat/include_mixing_stats["// int2str(j) // "]/mixing_bin_bounds")             
              buffer = field_tag(name=sfield%name, column=column+1, statistic="mixing_bins%" // trim(mixing_stats_name),&
                   & material_phase_name=material_phase_name, components=(shape_option(1)))
              write(diag_unit, '(a)') trim(buffer)
              column = column + (shape_option(1))
            end do
            
            ! Surface integrals
            do j = 0, option_count(trim(complete_field_path(sfield%option_path, stat = stat)) // "/stat/surface_integral") - 1
              call get_option(trim(complete_field_path(sfield%option_path)) &
              // "/stat/surface_integral[" // int2str(j) // "]/name", surface_integral_name)
              column = column + 1
              buffer = field_tag(sfield%name, column, "surface_integral%" // trim(surface_integral_name), material_phase_name)
              write(diag_unit, "(a)") trim(buffer)
           end do
           
         end do

         do i = 1, size(vfield_list(phase)%ptr)
           ! Headers for output statistics for each vector field
           vfield => extract_vector_field(state(phase), &
             & vfield_list(phase)%ptr(i))

           ! Standard scalar field stats for vector field magnitude
           if(stat_field(vfield, state(phase))) then
             column = column + 1
             buffer = field_tag(name=trim(vfield%name) // "%magnitude", column=column, &
               & statistic="min", material_phase_name=material_phase_name)
             write(diag_unit, '(a)') trim(buffer)

             column = column + 1
             buffer = field_tag(name=trim(vfield%name) // "%magnitude", column=column, &
               & statistic="max",material_phase_name= material_phase_name)
             write(diag_unit, '(a)') trim(buffer)

             column = column + 1
             buffer = field_tag(name=trim(vfield%name) // "%magnitude", column=column, &
               & statistic="l2norm", material_phase_name=material_phase_name)
             write(diag_unit, '(a)') trim(buffer)
           end if

           ! Standard scalar field stats for vector field components
           if(stat_field(vfield, state(phase), test_for_components = .true.)) then
             do j = 1, vfield%dim
                vfield_comp = extract_scalar_field(vfield, j)

               column = column + 1
               buffer=field_tag(name=vfield_comp%name, column=column, statistic="min", &
                 & material_phase_name=material_phase_name)
               write(diag_unit, '(a)') trim(buffer)

               column = column + 1
               buffer=field_tag(name=vfield_comp%name, column=column, statistic="max", &
                 & material_phase_name=material_phase_name)
               write(diag_unit, '(a)') trim(buffer)

               column = column + 1
               buffer=field_tag(name=vfield_comp%name, column=column, statistic="l2norm", &
                 & material_phase_name=material_phase_name)
               write(diag_unit, '(a)') trim(buffer)

               column = column + 1
               buffer=field_tag(name=vfield_comp%name, column=column, statistic="integral", &
                 & material_phase_name=material_phase_name)
               write(diag_unit, '(a)') trim(buffer)
             end do
           end if
           
           ! Surface integrals
           do j = 0, option_count(trim(complete_field_path(vfield%option_path, stat = stat)) // "/stat/surface_integral") - 1
             call get_option(trim(complete_field_path(vfield%option_path)) &
             // "/stat/surface_integral[" // int2str(j) // "]/name", surface_integral_name)
             column = column + 1
             buffer = field_tag(vfield%name, column, "surface_integral%" // trim(surface_integral_name), material_phase_name)
             write(diag_unit, "(a)") trim(buffer)
           end do

           ! drag calculation
           if(stat_field(vfield, state(phase), test_for_components = .true.)) then
            if(have_option(trim(complete_field_path(vfield%option_path, stat)) // "/stat/compute_body_forces_on_surfaces")) then
             do j = 1, mesh_dim(vfield%mesh)
               column = column + 1
               buffer = field_tag(name=trim(vfield%name), column=column, statistic="force%" &
               // int2str(j), material_phase_name=material_phase_name)
               write(diag_unit, '(a)') trim(buffer)
             end do
             if(have_option(trim(complete_field_path(vfield%option_path, stat)) // "/stat/compute_body_forces_on_surfaces/output_terms")) then
               do j = 1, mesh_dim(vfield%mesh)
                 column = column + 1
                 buffer = field_tag(name=trim(vfield%name), column=column, statistic="pressure_force%" &
                 // int2str(j), material_phase_name=material_phase_name)
                 write(diag_unit, '(a)') trim(buffer)
               end do
               do j = 1, mesh_dim(vfield%mesh)
                 column = column + 1
                 buffer = field_tag(name=trim(vfield%name), column=column, statistic="viscous_force%" &
                 // int2str(j), material_phase_name=material_phase_name)
                 write(diag_unit, '(a)') trim(buffer)
               end do
             end if
            end if
           end if

           ! momentum conservation error calculation
           if(stat_field(vfield, state(phase))) then
            if(have_option(trim(complete_field_path(vfield%option_path, stat)) // "/stat/calculate_momentum_conservation_error")) then
             do j = 1, mesh_dim(vfield%mesh)
               column = column + 1
               buffer = field_tag(name=trim(vfield%name), column=column, statistic="momentum_conservation%" &
               // int2str(j), material_phase_name=material_phase_name)
               write(diag_unit, '(a)') trim(buffer)
             end do
            end if
           end if

         end do

      end do phaseloop

      write(diag_unit, '(a)') "</header>"
      flush(diag_unit)
    end if

    call initialise_detectors(filename, state)

    ewrite(1, *) "Exiting initialise_diagnostics"

  end subroutine initialise_diagnostics

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

    ! Idempotency variable
    logical, save :: initialised=.false.

    integer :: column, i, j, phase
    character(len = 254) :: buffer, material_phase_name
    type(scalar_field) :: vfield_comp
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield

    ! Idempotency check
    if (initialised) return
    initialised=.true.

    if(have_option("/io/convergence/convergence_file")) then
       write_convergence_file = .true.
    else
       write_convergence_file = .false.
       return
    end if

    ! Only the first process should write convergence information
    if(getprocno() == 1) then
      conv_unit=free_unit()
      open(unit=conv_unit, file=trim(filename)//'.convergence', action="write")
    else
      return
    end if

    write(conv_unit, '(a)') "<header>"

    column=0

    ! Initial columns are elapsed time, dt and global iteration
    column=column+1
    buffer=field_tag(name="ElapsedTime", column=column, statistic="value")
    write(conv_unit, '(a)') trim(buffer)
    column=column+1
    buffer=field_tag(name="dt", column=column, statistic="value")
    write(conv_unit, '(a)') trim(buffer)
    column=column+1! Vector field magnitude
    buffer=field_tag(name="Iteration", column=column, statistic="value")
    write(conv_unit, '(a)') trim(buffer)

    phaseloop: do phase=1,size(state)

       material_phase_name=trim(state(phase)%name)

       do i=1, size(sfield_list(phase)%ptr)
          ! Output convergence information for each scalar field.
          sfield => extract_scalar_field(state(phase), sfield_list(phase)%ptr(i))

          if(.not. convergence_field(sfield)) then
            cycle
          end if

          column=column+1
          buffer=field_tag(name=sfield%name, column=column, statistic="error", material_phase_name=material_phase_name)
          write(conv_unit, '(a)') trim(buffer)

       end do

         do i = 1, size(vfield_list(phase)%ptr)
           ! Headers for output convergence information for each vector field

           vfield => extract_vector_field(state(phase), &
             & vfield_list(phase)%ptr(i))

           if(.not. convergence_field(vfield)) then
             cycle
           end if

           column = column + 1
           buffer = field_tag(name=trim(vfield%name) // "%magnitude", column=column, &
             & statistic="error", material_phase_name=material_phase_name)
           write(conv_unit, '(a)') trim(buffer)

           if(.not. convergence_field(vfield, test_for_components = .true.)) then
             cycle
           end if

           do j = 1, mesh_dim(vfield%mesh)
             vfield_comp = extract_scalar_field(vfield, j)

             column = column + 1
             buffer=field_tag(name=vfield_comp%name, column=column, statistic="error", &
               & material_phase_name=material_phase_name)
             write(conv_unit, '(a)') trim(buffer)

           end do

         end do

    end do phaseloop

    write(conv_unit, '(a)') "</header>"

  end subroutine initialise_convergence
  
  subroutine initialise_steady_state(filename, state)
    !!< Set up the steady state file headers.

    character(len=*) :: filename
    type(state_type), dimension(:), intent(in) :: state

    ! Idempotency variable
    logical, save :: initialised=.false.

    integer :: column, i, j, phase
    character(len = 254) :: buffer, material_phase_name
    type(scalar_field) :: vfield_comp
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield

    ! Idempotency check
    if (initialised) return
    initialised=.true.

    write_steady_state_file = have_option("/timestepping/steady_state/steady_state_file")
    if(.not. write_steady_state_file) return
    if(have_option("/timestepping/steady_state/steady_state_file/binary_output")) then
      binary_steady_state_output = .true.
    else if(have_option("/timestepping/steady_state/steady_state_file/plain_text_output")) then
      binary_steady_state_output = .false.
    else
      FLAbort("Unable to determine steady state output format")
    end if
    
    ! Only the first process should write steady state information
    if(getprocno() /= 1) return
    
    steady_state_unit=free_unit()
    open(unit=steady_state_unit, file=trim(filename)//'.steady_state', action="write")

    write(steady_state_unit, '(a)') "<header>"

    call initialise_constant_diagnostics(steady_state_unit, binary_format = binary_steady_state_output)

    ! Initial columns are elapsed time, dt and global iteration
    column=1
    buffer=field_tag(name="ElapsedTime", column=column, statistic="value")
    write(steady_state_unit, '(a)') trim(buffer)
    column=column+1
    buffer=field_tag(name="dt", column=column, statistic="value")
    write(steady_state_unit, '(a)') trim(buffer)

    phaseloop: do phase=1,size(state)
       material_phase_name = state(phase)%name

       do i = 1, scalar_field_count(state(phase))
          sfield => extract_scalar_field(state(phase), i)
          if(.not. steady_state_field(sfield)) cycle
          ! Scalar fields

          column=column+1
          buffer=field_tag(name=sfield%name, column=column, statistic="error", material_phase_name=material_phase_name)
          write(steady_state_unit, '(a)') trim(buffer)
       end do

       do i = 1, vector_field_count(state(phase))
         vfield => extract_vector_field(state(phase), i)
         if(.not. steady_state_field(vfield)) cycle         
         ! Vector fields
         
         column = column + 1
         buffer = field_tag(name=trim(vfield%name) // "%magnitude", column=column, &
           & statistic="error", material_phase_name=material_phase_name)
         write(steady_state_unit, '(a)') trim(buffer)

         if(.not. steady_state_field(vfield, test_for_components = .true.)) cycle         
         ! Vector field components
         
         do j = 1, mesh_dim(vfield%mesh)
           vfield_comp = extract_scalar_field(vfield, j)

           column = column + 1
           buffer=field_tag(name=vfield_comp%name, column=column, statistic="error", &
             & material_phase_name=material_phase_name)
           write(steady_state_unit, '(a)') trim(buffer)
         end do
       end do

    end do phaseloop

    column = column + 1
    buffer = field_tag(name = "MaxChange", column=column, statistic="value")
    write(steady_state_unit, '(a)') trim(buffer)

    write(steady_state_unit, '(a)') "</header>"
    flush(steady_state_unit)
    
    if(binary_steady_state_output) then
        close(steady_state_unit)

#ifdef STREAM_IO
      open(unit = steady_state_unit, file = trim(filename) // '.steady_state.dat', &
        & action = "write", access = "stream", form = "unformatted", status = "replace")
#else
      FLAbort("No stream I/O support")
#endif
    end if

  end subroutine initialise_steady_state
  
  subroutine initialise_detectors(filename, state)
    !!< Set up the detector file headers. This has the same syntax as the
    !!< .stat file

    character(len = *), intent(in) :: filename
    type(state_type), dimension(:), intent(in) :: state

    character(len=FIELD_NAME_LEN) ::funcnam, temp_name
    character(len=PYTHON_FUNC_LEN) :: func

    ! Idempotency variable
    logical, save :: initialised=.false.

    integer :: column, i, j, k, phase, m, IERROR
    integer :: static_dete, python_functions_or_files, total_dete, total_dete_groups, lagrangian_dete
    integer :: python_dete, ndete, dim, str_size, type_det
    integer, dimension(2) :: shape_option
    character(len = 254) :: buffer, material_phase_name, fmt
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    real, allocatable, dimension(:,:) :: coords
    real:: current_time
    character(len = OPTION_PATH_LEN) :: detectors_cp_filename, detector_file_filename
    logical :: detectors_checkpoint_done=.false., detectors_from_file_initially=.false.

    type(detector_type), pointer :: node

    type(element_type), pointer :: shape

    ! Idempotency check
    if (initialised) return
    initialised=.true.

    ! Check whether there are actually any detectors.
    static_dete = option_count("/io/detectors/static_detector")
    lagrangian_dete = option_count("/io/detectors/lagrangian_detector")
    python_functions_or_files = option_count("/io/detectors/detector_array")
    python_dete = 0
 
    do i=1,python_functions_or_files
       write(buffer, "(a,i0,a)")  &
            "/io/detectors/detector_array[",i-1,"]"
       call get_option(trim(buffer)//"/number_of_detectors", j)
       python_dete=python_dete+j
    end do
   
    total_dete=static_dete+lagrangian_dete+python_dete
    total_num_det=total_dete

!    allocate(detector_list(total_dete))
    
    total_dete_groups=static_dete+lagrangian_dete+python_functions_or_files
 
    allocate(name_of_detector_groups_in_read_order(total_dete_groups))
    allocate(number_det_in_each_group(total_dete_groups))
    allocate(name_of_detector_in_read_order(total_dete))
    
!   if (size(detector_list)==0) return

    if (total_dete==0) return

    vfield=>extract_vector_field(state(1), "Coordinate")
    shape=>ele_shape(vfield,1)
    call get_option("/geometry/dimension",dim)
    call get_option("/timestepping/current_time", current_time)
    
    ! Retrieve the position of each detector. If the option
    ! "from_checkpoint_file" exists, it means we are continuing the simulation
    ! after checkpointing and the reading of the detector positions must be
    ! done from a file

    if (have_option("/io/detectors/static_detector/from_checkpoint_file").or. & 
& have_option("/io/detectors/lagrangian_detector/from_checkpoint_file").or. &
& have_option("/io/detectors/detector_array/from_checkpoint_file")) then
       detectors_checkpoint_done=.true.
    else
       detectors_checkpoint_done=.false.
    end if

    if (.not.detectors_checkpoint_done) then

          do i=1,static_dete
             write(buffer, "(a,i0,a)")  &
               "/io/detectors/static_detector[",i-1,"]"
         
             allocate(node)
             call insert_det(detector_list,node) 

             call get_option(trim(buffer)//"/name", node%name)
             shape_option=option_shape(trim(buffer)//"/location")
             allocate(node%position(shape_option(1)))
             call get_option(trim(buffer)//"/location", node%position)       

             node%local = .not. isparallel()
             node%type = STATIC_DETECTOR
             node%id_number = i
             
       ! The arrays below contain information about the order in which detector
       ! groups are read and how many detectors there are in each group. This is
       ! used when checkpointing detectors. In particular, when continuing a
       ! simulation from a checkpoint, with these arrays we make sure we read
       ! back the detectors from the file in the same order than at the beginning
       ! of the simulation for consistency. All the .detectors files with
       ! detector data (position, value of variables at those positions, etc.) 
       ! will have the information in the same order.

             name_of_detector_groups_in_read_order(i)=node%name
             number_det_in_each_group(i)=1.0

             name_of_detector_in_read_order(i)=node%name

             allocate(node%local_coords(local_coord_count(shape)))

          end do

          do i=1,lagrangian_dete
             write(buffer, "(a,i0,a)")  &
               "/io/detectors/lagrangian_detector[",i-1,"]"

             allocate(node)
             call insert_det(detector_list,node) 

             call get_option(trim(buffer)//"/name", node%name)
             shape_option=option_shape(trim(buffer)//"/location")
             allocate(node%position(shape_option(1)))
             call get_option(trim(buffer)//"/location", node%position)
          
             node%type = LAGRANGIAN_DETECTOR
             node%local = .true.
             node%id_number = i+static_dete
             
             name_of_detector_groups_in_read_order(i+static_dete)=node%name
             number_det_in_each_group(i+static_dete)=1.0
             name_of_detector_in_read_order(i+static_dete)=node%name

             allocate(node%local_coords(local_coord_count(shape)))

          end do

          k=static_dete+lagrangian_dete+1

          do i=1,python_functions_or_files
             write(buffer, "(a,i0,a)")  &
               "/io/detectors/detector_array[",i-1,"]"

             call get_option(trim(buffer)//"/name", funcnam)
             call get_option(trim(buffer)//"/number_of_detectors", ndete)
             str_size=len_trim(int2str(ndete))
             fmt="(a,I"//int2str(str_size)//"."//int2str(str_size)//")"

             if (have_option(trim(buffer)//"/lagrangian")) then
                  type_det=LAGRANGIAN_DETECTOR
             else
                  type_det=STATIC_DETECTOR
             end if

             name_of_detector_groups_in_read_order(i+static_dete+lagrangian_dete)=trim(funcnam)
             number_det_in_each_group(i+static_dete+lagrangian_dete)=ndete

             if (have_option(trim(buffer) // "/from_file")) then
                 detectors_from_file_initially=.true.
             else
                 detectors_from_file_initially=.false.
             end if

             if (.not.detectors_from_file_initially) then

                 ! Reading detectors from a python function

                 call get_option(trim(buffer)//"/python", func)
                 allocate(coords(dim,ndete))
                 call set_detector_coords_from_python(coords, ndete, func, current_time)
          
                 do j=1,ndete

                     allocate(node)
                     call insert_det(detector_list,node)

                     write(node%name, fmt) trim(funcnam)//"_", j

                     allocate(node%position(dim))

                     node%position=coords(:,j)
                     node%local = type_det == LAGRANGIAN_DETECTOR .or. .not. isparallel()
                     node%type=type_det

                     name_of_detector_in_read_order(k)=trim(funcnam)

                     node%id_number = k

                     allocate(node%local_coords(local_coord_count(shape)))

                     k=k+1
           
                 end do
                 deallocate(coords)

             else
                 ! Reading from a binary file where the user has placed the detector positions

                 if (getprocno() == 1) then

                 detector_file_unit=free_unit()

                 call get_option("io/detectors/detector_array/from_file/file_name",detector_file_filename)

#ifdef STREAM_IO
      open(unit = detector_file_unit, file = trim(detector_file_filename), &
        & action = "read", access = "stream", form = "unformatted")
#else
      FLAbort("No stream I/O support")
#endif
       
                 do j=1,ndete

                     allocate(node)
                     call insert_det(detector_list,node)

                     write(node%name, fmt) trim(funcnam)//"_", j

                     allocate(node%position(dim))
                     read(detector_file_unit) node%position
                     node%local = type_det == LAGRANGIAN_DETECTOR .or. .not. isparallel()
                     node%type=type_det
                     node%id_number = k

                     name_of_detector_in_read_order(k)=trim(funcnam)

                     allocate(node%local_coords(local_coord_count(shape)))
                     k=k+1
           
                 end do

                 end if
                
             end if

          end do
      
    else 
    
       !!If reading from checkpoint file
       !!!First we should read the header of checkpoint_file to get the order in which the detectors were read at the beginning of the calculation

        !IN PARALLEL DET CODE ALL PROC NEED TO READ ALL DET 

           detector_checkpoint_unit=free_unit()

           if (have_option("/io/detectors/static_detector")) then
               call get_option("io/detectors/static_detector/from_checkpoint_file/file_name",detectors_cp_filename)  !!THIS NAME ends in _det. Need to add .groups for the name of the file with the header. The binary file with the positions is called .positions.dat
           elseif (have_option("/io/detectors/lagrangian_detector")) then 
               call get_option("io/detectors/lagrangian_detector/from_checkpoint_file/file_name",detectors_cp_filename)  !!THIS NAME ends in _det. Need to add .groups for the name of the file with the header. The binary file with the positions is called .positions.dat
           else 
               call get_option("io/detectors/detector_array/from_checkpoint_file/file_name",detectors_cp_filename)  !!THIS NAME ends in _det. Need to add .groups for the name of the file with the header. The binary file with the positions is called .positions.dat
           end if 

           open(unit=detector_checkpoint_unit, file=trim(detectors_cp_filename) // '.groups', action="read") 

           do i=1, total_dete_groups  

               read(detector_checkpoint_unit,'(a,i10)') name_of_detector_groups_in_read_order(i), number_det_in_each_group(i)

           end do

           close(detector_checkpoint_unit)

#ifdef STREAM_IO
      open(unit = detector_checkpoint_unit, file = trim(detectors_cp_filename) // '.positions.dat', &
        & action = "read", access = "stream", form = "unformatted")
#else
      FLAbort("No stream I/O support")
#endif
 
       !!!Read in order the last positions of the detectors from the binary file.

          do j=1,size(name_of_detector_groups_in_read_order)

             do i=1,static_dete
                 write(buffer, "(a,i0,a)")  &
               "/io/detectors/static_detector[",i-1,"]"
                 call get_option(trim(buffer)//"/name", temp_name)
          
                 if (name_of_detector_groups_in_read_order(j)==temp_name) then

                     allocate(node)
                     call insert_det(detector_list,node)
                     node%name=temp_name
                     allocate(node%position(dim))
                     read(detector_checkpoint_unit) node%position
                     node%local = .not. isparallel()
                     node%type = STATIC_DETECTOR
                     allocate(node%local_coords(local_coord_count(shape)))   

                     node%id_number = i
                  
                 else
                     cycle
                 end if
             end do

          end do

          do j=1,size(name_of_detector_groups_in_read_order)   

             do i=1,lagrangian_dete
                 write(buffer, "(a,i0,a)")  &
               "/io/detectors/lagrangian_detector[",i-1,"]"
                 call get_option(trim(buffer)//"/name", temp_name)

                 if (name_of_detector_groups_in_read_order(j)==temp_name) then


                     allocate(node)
                     call insert_det(detector_list,node)
                     node%name=temp_name
                     allocate(node%position(dim))
                     read(detector_checkpoint_unit) node%position
                     node%local = type_det == LAGRANGIAN_DETECTOR .or. .not. isparallel()
                     node%type = LAGRANGIAN_DETECTOR

                     allocate(node%local_coords(local_coord_count(shape))) 

                     node%id_number = i+static_dete
    
                 else
                     cycle
                 end if
             end do

          end do

          k=static_dete+lagrangian_dete+1

          do j=1,size(name_of_detector_groups_in_read_order)  

             do i=1,python_functions_or_files
                 write(buffer, "(a,i0,a)")  &
               "/io/detectors/detector_array[",i-1,"]"
      
                 call get_option(trim(buffer)//"/name", temp_name)

                 if (name_of_detector_groups_in_read_order(j)==temp_name) then

                    call get_option(trim(buffer)//"/number_of_detectors", ndete)
                    str_size=len_trim(int2str(ndete))
                    fmt="(a,I"//int2str(str_size)//"."//int2str(str_size)//")"

                    if (have_option(trim(buffer)//"/lagrangian")) then
                      type_det=LAGRANGIAN_DETECTOR
                    else
                      type_det=STATIC_DETECTOR
                    end if

                    do m=1,number_det_in_each_group(j)

                      allocate(node)
                      call insert_det(detector_list,node)

                      write(node%name, fmt) trim(temp_name)//"_", m
                      allocate(node%position(dim))
                      read(detector_checkpoint_unit) node%position
                      node%local = type_det == LAGRANGIAN_DETECTOR .or. .not. isparallel()
                      node%type=type_det

                      node%id_number = k

                      allocate(node%local_coords(local_coord_count(shape)))

                      k=k+1
           
                    end do

                 else 
                    
                    cycle  
                 
                 end if
             
             end do

          end do
 
        !!IN PARALLEL DET CODE ALL PROC NEED TO READ ALL DET end if

    end if


    allocate (detector_sfield_list(size(state)))
    do phase=1, size(state)
         allocate(detector_sfield_list(phase)%ptr(&
         size(state(phase)%scalar_names)))
         detector_sfield_list(phase)%ptr=state(phase)%scalar_names
    end do

    allocate (detector_vfield_list(size(state)))
    do phase = 1, size(state)
         allocate(detector_vfield_list(phase)%ptr(&
         size(state(phase)%vector_names)))
         detector_vfield_list(phase)%ptr = state(phase)%vector_names
    end do

    binary_detector_output = have_option("/io/detectors/binary_output")

    ! Only the first process should write statistics information
    if (getprocno() == 1) then
    
    detector_unit=free_unit()
    open(unit=detector_unit, file=trim(filename)//'.detectors', action="write")

   
    write(detector_unit, '(a)') "<header>"

    call initialise_constant_diagnostics(detector_unit, binary_format = binary_detector_output)

    column=0

    ! Initial columns are elapsed time and dt.
    column=column+1
    buffer=field_tag(name="ElapsedTime", column=column, statistic="value")
    write(detector_unit, '(a)') trim(buffer)
    column=column+1
    buffer=field_tag(name="dt", column=column, statistic="value")
    write(detector_unit, '(a)') trim(buffer)

    ! Next columns contain the positions of all the detectors.

    node => detector_list%firstnode

    positionloop: do i=1, detector_list%length


         buffer=field_tag(name=node%name, column=column+1,&
            statistic="position", &
            components=size(node%position))
         write(detector_unit, '(a)') trim(buffer)
         column=column+size(node%position)

         node => node%next

    end do positionloop

    
    phaseloop: do phase=1,size(state)

            material_phase_name=trim(state(phase)%name)

        do i=1, size(detector_sfield_list(phase)%ptr)
          ! Headers for detectors for each scalar field.
            sfield => extract_scalar_field(state(phase), &
               detector_sfield_list(phase)%ptr(i))

            if(.not. detector_field(sfield)) then
                 cycle
            end if

          node => detector_list%firstnode

          do j=1, detector_list%length

            column=column+1
            buffer=field_tag(name=sfield%name, column=column, &
                  statistic=node%name, &
                  material_phase_name=material_phase_name)
            write(detector_unit, '(a)') trim(buffer)

            node => node%next

          end do

        end do

        do i = 1, size(detector_vfield_list(phase)%ptr)
          ! Headers for detectors for each vector field.
          vfield => extract_vector_field(state(phase), &
               & detector_vfield_list(phase)%ptr(i))

          if(.not. detector_field(vfield)) then
             cycle
          end if

          node => detector_list%firstnode

          do j=1, detector_list%length

             buffer=field_tag(name=vfield%name, column=column+1, &
                  statistic=node%name, &
                  material_phase_name=material_phase_name, &
                  components=vfield%dim)
             write(detector_unit, '(a)') trim(buffer)
             column=column+size(node%position)

             node => node%next

          end do
        end do

    end do phaseloop

    write(detector_unit, '(a)') "</header>"
    flush(detector_unit)

    !!!when using mpi_subroutines to write into the detectors file we need to close the file since 
    !filename.detectors.dat needs to be open now with MPI_OPEN

    if ((.not.isparallel()).and.(.not. binary_detector_output)) then
!        close(detector_unit)

!#ifdef STREAM_IO
!      open(unit = detector_unit, file = trim(filename) // '.detectors.dat', &
!        & action = "write", access = "stream", form = "unformatted")
!#else
!     FLAbort("No stream I/O support")
!#endif

    else
    
      close(detector_unit)

    end if

    end if

     if ((isparallel()).or.((.not.isparallel()).and.(binary_detector_output))) then

    call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(filename) // '.detectors.dat', MPI_MODE_CREATE + MPI_MODE_RDWR, MPI_INFO_NULL, fh, IERROR)

    end if 

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
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time, dt
    integer, intent(in) :: timestep
    logical, intent(in), optional :: not_to_move_det_yet 

    character(len = 2 + real_format_len(padding = 1) + 1) :: format, format2, format3, format4
    integer :: i, j, k, phase, stat
    integer, dimension(2) :: shape_option
    integer :: nodes, elements, surface_elements
    real :: fmin, fmax, fnorm2, fintegral, fnorm2_cv, fintegral_cv, surface_integral
    real, dimension(:), allocatable :: f_mix_fraction
    type(mesh_type), pointer :: mesh
    type(scalar_field) :: vfield_comp
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(vector_field) :: xfield

    ewrite(1,*) 'In write_diagnostics'
    call profiler_tic("I/O")

    format="(" // real_format(padding = 1) // ")"
    format2="(2" // real_format(padding = 1) // ")"
    format3="(3" // real_format(padding = 1) // ")"
    format4="(4" // real_format(padding = 1) // ")"

    ! Only the first process should write statistics information (but all must
    ! be involved in calculating them)
    if(getprocno() == 1) then
      write(diag_unit, trim(format), advance="no") time
      write(diag_unit, trim(format), advance="no") dt
      write(diag_unit, trim(format), advance="no") elapsed_walltime()
    end if

    do i = 1, size(mesh_list)
      ! Output statistics for each mesh
      mesh => extract_mesh(state(1), mesh_list(i))

      if(stat_mesh(mesh)) then
        call mesh_stats(mesh, nodes, elements, surface_elements)
        if(getprocno() == 1) then
          write(diag_unit, "(a,i0,a,i0,a,i0)", advance = "no") " ", nodes, " ", elements, " ", surface_elements
        end if
      end if
    end do

#ifdef HAVE_MEMORY_STATS
    ! Memory statistics.
    call write_memory_stats(diag_unit, format)
    call reset_memory_logs
#endif

    phaseloop: do phase=1,size(state)

       scalar_field_loop: do i=1, size(sfield_list(phase)%ptr)
          ! Output statistics for each scalar field
          sfield=>extract_scalar_field(state(phase), sfield_list(phase)%ptr(i))

          xfield=get_diagnostic_coordinate_field(state(phase), sfield%mesh)

          ! Standard scalar field stats
          if(stat_field(sfield, state(phase))) then
          
            call field_stats(sfield, Xfield, fmin, fmax, fnorm2, fintegral)
            if(getprocno() == 1) then
              write(diag_unit, trim(format4), advance="no") fmin, fmax, fnorm2,&
                   & fintegral
            end if
            
          end if

          ! Control volume stats
          if(have_option(trim(complete_field_path(sfield%option_path,stat=stat)) //&
               & "/stat/include_cv_stats")) then

            call field_cv_stats(sfield, Xfield, fnorm2_cv, fintegral_cv)

            ! Only the first process should write statistics information
            if(getprocno() == 1) then
              write(diag_unit, trim(format2), advance="no") fnorm2_cv, fintegral_cv
            end if

          end if

          ! Mixing stats
          do j = 0, option_count(trim(complete_field_path(sfield%option_path, stat=stat)) // "/stat/include_mixing_stats") - 1
            shape_option=option_shape(trim(complete_field_path(sfield%option_path)) // &
                 & "/stat/include_mixing_stats["// int2str(j) // "]/mixing_bin_bounds")           
            allocate(f_mix_fraction(1:shape_option(1)))
            f_mix_fraction = 0.0
            
            call mixing_stats(f_mix_fraction, sfield, Xfield, mixing_stats_count = j)          
         
            if(getprocno() == 1) then
               do k=1, (size(f_mix_fraction))
                  write(diag_unit, trim(format), advance="no") f_mix_fraction(k)
               end do
            end if

            deallocate(f_mix_fraction)

          end do
         
         ! Surface integrals
         do j = 0, option_count(trim(complete_field_path(sfield%option_path, stat = stat)) // "/stat/surface_integral") - 1
           surface_integral = calculate_surface_integral(sfield, xfield, j)
           ! Only the first process should write statistics information
           if(getprocno() == 1) then
             write(diag_unit, trim(format), advance = "no") surface_integral
           end if
         end do

         call deallocate(xfield)

       end do scalar_field_loop

       vector_field_loop: do i = 1, size(vfield_list(phase)%ptr)
         ! Output statistics for each vector field
         vfield => extract_vector_field(state(phase), &
           & vfield_list(phase)%ptr(i))
          
         xfield=get_diagnostic_coordinate_field(state(phase), vfield%mesh)

         ! Standard scalar field stats for vector field magnitude
         if(stat_field(vfield,state(phase))) then
           call field_stats(vfield, Xfield, fmin, fmax, fnorm2)
           ! Only the first process should write statistics information
           if(getprocno() == 1) then
             write(diag_unit, trim(format3), advance = "no") fmin, fmax, fnorm2
           end if
         end if

         ! Standard scalar field stats for vector field components
         if(stat_field(vfield, state(phase), test_for_components = .true.)) then
           do j = 1, mesh_dim(vfield%mesh)
             vfield_comp = extract_scalar_field(vfield, j)

             call field_stats(vfield_comp, Xfield, fmin, fmax, fnorm2, &
               & fintegral)
             ! Only the first process should write statistics information
             if(getprocno() == 1) then
               write(diag_unit, trim(format4), advance = "no") fmin, fmax, fnorm2, &
                 & fintegral
             end if
           end do
         end if
         
         ! Surface integrals
         do j = 0, option_count(trim(complete_field_path(vfield%option_path, stat = stat)) // "/stat/surface_integral") - 1
           surface_integral = calculate_surface_integral(vfield, xfield, j)
           ! Only the first process should write statistics information
           if(getprocno() == 1) then
             write(diag_unit, trim(format), advance = "no") surface_integral
           end if
         end do

         ! drag calculation
         if(stat_field(vfield, state(phase), test_for_components = .true.)) then
          if(have_option(trim(complete_field_path(vfield%option_path, stat)) // "/stat/compute_body_forces_on_surfaces")) then
            call write_body_forces(state(phase), vfield)  
          end if
         end if

         ! momentum conservation error calculation
         if(stat_field(vfield, state(phase))) then
          if(have_option(trim(complete_field_path(vfield%option_path, stat)) // "/stat/calculate_momentum_conservation_error")) then
            call write_momentum_conservation_error(state(phase), vfield)
          end if
         end if
         
         call deallocate(xfield)
         
       end do vector_field_loop

    end do phaseloop

    ! Output end of line
    ! Only the first process should write statistics information
    if(getprocno() == 1) then
      write(diag_unit,'(a)') ""
      flush(diag_unit)
    end if

    ! Now output any detectors.
    
    call write_detectors(state, time, dt, timestep, not_to_move_det_yet)

    call profiler_toc("I/O")
  
  contains
  
    subroutine write_body_forces(state, vfield)
      type(state_type), intent(in) :: state
      type(vector_field), intent(in) :: vfield
      
      integer :: i
      real :: force(vfield%dim), pressure_force(vfield%dim), viscous_force(vfield%dim)
    
      if(have_option(trim(complete_field_path(vfield%option_path, stat)) // "/stat/compute_body_forces_on_surfaces/output_terms")) then
        ! calculate the forces on the surface
        call diagnostic_body_drag(state, force, pressure_force = pressure_force, viscous_force = viscous_force)   
        if(getprocno() == 1) then
           do i=1, mesh_dim(vfield%mesh)
              write(diag_unit, trim(format), advance="no") force(i)
           end do
           do i=1, mesh_dim(vfield%mesh)
              write(diag_unit, trim(format), advance="no") pressure_force(i)
           end do
           do i=1, mesh_dim(vfield%mesh)
              write(diag_unit, trim(format), advance="no") viscous_force(i)
           end do
        end if
      else
        ! calculate the forces on the surface
        call diagnostic_body_drag(state, force) 
        if(getprocno() == 1) then
           do i=1, mesh_dim(vfield%mesh)
              write(diag_unit, trim(format), advance="no") force(i)
           end do
        end if     
      end if 
      
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
          write(diag_unit, trim(format), advance="no") momentum_cons(dim)
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

    if(write_convergence_file) then
      ! Only the first process should write convergence information
      if(getprocno() == 1) then
         write(conv_unit, format, advance="no") time
         write(conv_unit, format, advance="no") dt
         write(conv_unit, iformat, advance="no") it
      end if
    end if
    
    coordinates => extract_vector_field(state(1), "Coordinate")
    convergence_norm = convergence_norm_integer("/timestepping/nonlinear_iterations/tolerance")

    phaseloop: do phase=1,size(state)

       do i=1, size(sfield_list(phase)%ptr)
          ! Output convergence information for each scalar field.
          sfield=>extract_scalar_field(state(phase), &
               &                       sfield_list(phase)%ptr(i))

          if(.not. convergence_field(sfield)) then
            cycle
          end if

          nlsfield=>extract_scalar_field(state(phase), &
                                       "Iterated"//trim(sfield_list(phase)%ptr(i)))

          call field_con_stats(sfield, nlsfield, error, &
                               convergence_norm, coordinates)
          maxerror = max(maxerror, error)

          if(write_convergence_file) then
            ! Only the first process should write convergence information
            if(getprocno() == 1) then
               write(conv_unit, format, advance="no") error
            end if
          end if

       end do

       do i = 1, size(vfield_list(phase)%ptr)
         ! Output convergence information for each vector field

         vfield => extract_vector_field(state(phase), &
           & vfield_list(phase)%ptr(i))

         if(.not. convergence_field(vfield)) then
           cycle
         end if

         nlvfield => extract_vector_field(state(phase), &
           & "Iterated"//vfield_list(phase)%ptr(i))

         call field_con_stats(vfield, nlvfield, error, &
                              convergence_norm, coordinates)
         maxerror = max(maxerror, error)

         if(write_convergence_file) then
            ! Only the first process should write convergence information
            if(getprocno() == 1) then
            write(conv_unit, format, advance = "no") error
            end if
         end if

         if(.not. convergence_field(vfield, test_for_components = .true.)) then
           cycle
         end if

         do j = 1, mesh_dim(vfield%mesh)
           vfield_comp = extract_scalar_field(vfield, j)
           nlvfield_comp = extract_scalar_field(nlvfield, j)

           call field_con_stats(vfield_comp, nlvfield_comp, error, &
                                convergence_norm, coordinates)
           maxerror = max(maxerror, error)

           if(write_convergence_file) then
               ! Only the first process should write convergence information
               if(getprocno() == 1) then
                  write(conv_unit, format, advance = "no") error
               end if
           end if
         end do
       end do

    end do phaseloop

    if(write_convergence_file) then
      ! Output end of line
      ! Only the first process should write convergence information
      if(getprocno() == 1) then
         write(conv_unit,'(a)') ""
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
    if(write_steady_state_file .and. procno == 1) then    
      call get_option("/timestepping/current_time", elapsed_time)
      if(binary_steady_state_output) then
        write(steady_state_unit) elapsed_time
        write(steady_state_unit) dt
      else
        write(steady_state_unit, format, advance="no") elapsed_time
        write(steady_state_unit, format, advance="no") dt
      end if
    end if

    phaseloop: do phase=1,size(state)

       do i=1, size(sfield_list(phase)%ptr)
          ! Test steady state information for each scalar field.

          sfield=>extract_scalar_field(state(phase), i)
          if(.not. steady_state_field(sfield)) cycle
          ! Scalar fields

          oldsfield=>extract_scalar_field(state(phase), &
                                       "Old"//trim(sfield_list(phase)%ptr(i)))

          call field_con_stats(sfield, oldsfield, change, &
                               convergence_norm, coordinates)
          if(acceleration) change = change/dt
          ewrite(2, *) trim(sfield%name), change
          maxchange = max(maxchange, change)

          if(write_steady_state_file .and. procno == 1) then
            if(binary_steady_state_output) then
              write(steady_state_unit) change
            else
              write(steady_state_unit, format, advance = "no") change
            end if
          end if

       end do

       do i = 1, vector_field_count(state(phase))
         vfield => extract_vector_field(state(phase), i)
         if(.not. steady_state_field(vfield)) cycle
         ! Vector fields


         oldvfield => extract_vector_field(state(phase), &
           & "Old"//vfield_list(phase)%ptr(i))

         call field_con_stats(vfield, oldvfield, change, &
                              convergence_norm, coordinates)
         if(acceleration) change = change/dt
         ewrite(2, *) trim(vfield%name), change
         maxchange = max(maxchange, change)

         if(write_steady_state_file .and. procno == 1) then
            if(binary_steady_state_output) then
              write(steady_state_unit) change
            else
              write(steady_state_unit, format, advance = "no") change
            end if
         end if

         if(.not. steady_state_field(vfield, test_for_components = .true.)) cycle
         ! Vector field components

         do j = 1, mesh_dim(vfield%mesh)
           vfield_comp = extract_scalar_field(vfield, j)
           oldvfield_comp = extract_scalar_field(oldvfield, j)

           call field_con_stats(vfield_comp, oldvfield_comp, change, &
                                convergence_norm, coordinates)
           if(acceleration) change = change/dt
           ewrite(2, *) trim(vfield%name), j, change
           maxchange = max(maxchange, change)

           if(write_steady_state_file .and. procno == 1) then
             if(binary_steady_state_output) then
               write(steady_state_unit) change
             else
               write(steady_state_unit, format, advance = "no") change
             end if
           end if
         end do

       end do

    end do phaseloop

    ewrite(1, *) "maxchange = ", maxchange
    
    if(write_steady_state_file .and. procno == 1) then
      if(binary_steady_state_output) then
        write(steady_state_unit) maxchange
      else
        write(steady_state_unit, format, advance = "no") maxchange      
        ! Output end of line
        write(steady_state_unit,'(a)') ""
      end if
      
      flush(steady_state_unit)
    end if

    ewrite(1, *) "Exiting test_and_write_steady_state"

  end subroutine test_and_write_steady_state

  subroutine write_detectors(state, time, dt, timestep, not_to_move_det_yet)
    !!< Write the field values at detectors to the previously opened detectors file.
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time, dt
    integer, intent(in) :: timestep
    logical, intent(in), optional :: not_to_move_det_yet 

    character(len=10) :: format_buffer
    integer :: i, j, k, phase, ele, processor_number, num_proc, dimen, number_neigh_processors, all_send_lists_empty, current_proc, nprocs
    real :: value
    real, dimension(:), allocatable :: vvalue
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield, xfield
!    character(len=254), dimension(:), allocatable :: types_det
    integer, dimension(:), allocatable :: types_det
    logical :: any_lagrangian
    integer, dimension(:), allocatable :: processor_number_array

    type(detector_type), pointer :: node, temp_node
    type(integer_hash_table) :: ihash, ihash_inverse, ihash_neigh_ele

    integer, dimension(:), allocatable :: global_det_count

    type(detector_linked_list), dimension(:), allocatable :: send_list_array, receive_list_array
    integer :: len_set, target_proc_a, mapped_val_a, halo_level, communicator, nhalos, check_no_det

    type(mesh_type), pointer :: mesh_ele

    type(halo_type), pointer :: ele_halo, node_halo 

    ewrite(1,*) "Inside write_detectors subroutine"

    xfield=>extract_vector_field(state(1), "Coordinate")

    check_no_det=1
    if (detector_list%length==0) then
       check_no_det=0
    end if

    call allmax(check_no_det)

    if (check_no_det==0) then
       return
    end if

    if (detector_list%length/=0) then

       node => detector_list%firstnode

       do i = 1, detector_list%length
         
         node => node%next

       end do
       
      ! Calculate the location of the detectors in the mesh.

       call search_for_detectors(detector_list, xfield)

       node => detector_list%firstnode

       do i = 1, detector_list%length
         
         if (node%element<0) then

            node%local = .true.

         else

            node%initial_owner=getprocno()

            node%local = .true.
        
         end if 

         node => node%next

       end do

    !!! CREATE HERE THE ARRAY CALLED GLOBAL DET COUNT (GDC) THAT WE NEED BEFORE CALLING MPI_ALL_MAX() 
    !!! TO SOLVE THE CONFLICT OF OWNERSHIP OF THE DETECTORS.
    !!! AFTERWARDS UPDATE THE DETECTOR_LIST BY REMOVING THE DETECTORS (NODES IN THE LIST) WHERE GETPROCNO()/=GDC(i)
    !!! THIS IS ONLY NEEDED AT THE BEGINNING OF THE SIMULATION, at the first time step, 
    !!! ONCE WE HAVE DISTRIBUTED INITIALLY THE DETECTORS
    !!! AMONG THE DIFFERENT PROCESSORS, WE DON'T NEED TO DO THIS ANY MORE

        if (timestep==1) then

          allocate(global_det_count(detector_list%length))

          node => detector_list%firstnode

          do i = 1, detector_list%length

            global_det_count(i)=node%initial_owner

            node => node%next
           
          end do

          node => detector_list%firstnode
          
          do i = 1, size(global_det_count)

             call allmax(global_det_count(i))

             node => node%next

          end do

          node => detector_list%firstnode

          do i = 1, size(global_det_count)

             if (global_det_count(i)/=node%initial_owner) then

               temp_node => node
              
               if ((.not.associated(node%previous)).and.(detector_list%length/=1)) then
               !!this checks if the current node that we are going to remove from the list is the first 
               !!one in the list but not the only node in the list

                  node%next%previous => null()

                  node => node%next

                  temp_node%previous => null()
                  temp_node%next => null()

                  detector_list%firstnode => node
                  detector_list%firstnode%previous => null()

                  detector_list%length = detector_list%length-1   

               else 

                   if ((node%id_number==size(global_det_count)).and.(associated(node%previous))) then
                   !!this takes into account the case when the node is the last one in the list but not the only one

                        node%previous%next => null()

                        detector_list%lastnode => node%previous

                        temp_node%previous => null()
                        temp_node%next => null()

                        detector_list%lastnode%next => null()

                        detector_list%length = detector_list%length-1    

                   else    

                        if (detector_list%length==1) then
                        !!!This case takes into account if the list has only one node. 

                        temp_node%previous => null()
                        temp_node%next => null()

                        detector_list%firstnode => null()
                        detector_list%lastnode => null()

                        detector_list%length = detector_list%length-1    

                        else
                        !!case when the node is in the middle of the double linked list

                           node%previous%next => node%next

                           node%next%previous => node%previous

                           node => node%next

                           temp_node%previous => null()
                           temp_node%next => null()

                           detector_list%length = detector_list%length-1    

                        end if 
                   end if
               end if

             else 

               node => node%next

             end if
           
          end do   

          deallocate(global_det_count)

          node => detector_list%firstnode

          do i = 1, detector_list%length

            if (node%initial_owner==-1) then

               node%type = STATIC_DETECTOR
     
            end if

            node => node%next
            
          end do

       end if

    end if

    call allocate(ihash_neigh_ele)   
    processor_number=getprocno() 
    num_proc=1

    vfield => extract_vector_field(state(1),"Velocity")
    halo_level = element_halo_count(vfield%mesh)

    do ele = 1, element_count(vfield%mesh)

      processor_number=element_owner(vfield%mesh,ele)

      if ((processor_number/=getprocno()).and.(.not.has_key(ihash_neigh_ele, processor_number))) then

         call insert(ihash_neigh_ele, processor_number, num_proc)

         call fetch_pair(ihash_neigh_ele, num_proc, target_proc_a, mapped_val_a)

         num_proc=num_proc+1

      end if

    end do

    number_neigh_processors=0

    call allocate(ihash) 

    if (halo_level /= 0.) then

    ele_halo => vfield%mesh%element_halos(halo_level)

    nprocs = halo_proc_count(ele_halo)

    num_proc=1

    do i = 1, nprocs 

!!! An alternative and it seems better way to find out the neighbouring processors to a given processor is the following:   
      
      if ((halo_send_count(ele_halo, i) + halo_receive_count(ele_halo, i) > 0).and.(.not.has_key(ihash, i))) then

         call insert(ihash, i, num_proc)

         num_proc=num_proc+1

     end if

    end do

    do i=1, key_count(ihash)

       call fetch_pair(ihash, i, target_proc_a, mapped_val_a)
!       ewrite(1,*) "ihash, i is", i
!       ewrite(1,*) "ihash, pair", target_proc_a, mapped_val_a

    end do

    call allocate(ihash_inverse) 
    do i=1, key_count(ihash)

       call fetch_pair(ihash, i, target_proc_a, mapped_val_a)
       call insert(ihash_inverse, mapped_val_a, target_proc_a)

    end do

    do i=1, key_count(ihash_inverse)

       call fetch_pair(ihash_inverse, i, target_proc_a, mapped_val_a)
!       ewrite(1,*) "ihash_inverse, i is", i
!       ewrite(1,*) "ihash_inverse, pair", target_proc_a, mapped_val_a

    end do

    call deallocate(ihash_inverse) 

    do i=1, key_count(ihash_neigh_ele)

       call fetch_pair(ihash_neigh_ele, i, target_proc_a, mapped_val_a)
       ewrite(1,*) "ihash_neigh_ele, i is", i
       ewrite(1,*) "ihash_neigh_ele, pair", target_proc_a, mapped_val_a

    end do

    call deallocate(ihash_neigh_ele)

    number_neigh_processors=key_count(ihash)

    end if

    node => detector_list%firstnode

    do j=1, detector_list%length

       node%dt=dt

       node => node%next

    end do  

    !do_until_all_send_lists_between_processors_are_empty
    do  

       allocate(send_list_array(number_neigh_processors))
       allocate(receive_list_array(number_neigh_processors))

       node => detector_list%firstnode

       allocate(types_det(detector_list%length))

       any_lagrangian=.false.

       do i = 1, detector_list%length
         
         types_det(i) = node%type
   !      if (types_det(i)==LAGRANGIAN_DETECTOR)  then
         if (types_det(i)==2)  then
            any_lagrangian=.true. 
         end if

         node => node%next
         
       end do
          
       if (any_lagrangian) then

         ewrite(1,*) "SHOULD BE MOVING THE DETECTORS"

        if (.not.present(not_to_move_det_yet).and.(timestep/=0))  call move_detectors_bisection_method(state, dt, ihash, send_list_array) 

       end if

       node => detector_list%firstnode
       do i = 1, detector_list%length

         node => node%next
         
      end do

       all_send_lists_empty=number_neigh_processors
       do k=1, number_neigh_processors

          ewrite(1,*) "number_neigh_processors:", number_neigh_processors

          if (send_list_array(k)%length==0) then
            all_send_lists_empty=all_send_lists_empty-1
         end if

       end do

       ewrite(1,*) "all_send_lists_empty bef allmax. It is 0 if all list towards neigh proc. are empty:", all_send_lists_empty

       call allmax(all_send_lists_empty)

       ewrite(1,*) "all_send_lists_empty after allmax. It is 0 if all list towards neigh proc. are empty for all proc:", all_send_lists_empty

       if (all_send_lists_empty==0) exit

       node => detector_list%firstnode
       do i = 1, detector_list%length

          node => node%next
         
      end do

       if (timestep/=0) then

          call serialise_lists_exchange_receive(state,send_list_array,receive_list_array,number_neigh_processors,ihash)

       end if

       node => detector_list%firstnode
       do i = 1, detector_list%length

          node => node%next
         
       end do

       do i=1, number_neigh_processors

         if  (receive_list_array(i)%length/=0) then      
        
            call move_det_from_receive_list_to_det_list(detector_list,receive_list_array(i))

         end if

       end do

       node => detector_list%firstnode
       do i = 1, detector_list%length

          node => node%next
         
      end do
  
!!! BEFORE DEALLOCATING THE LISTS WE SHOULD MAKE SURE THEY ARE EMPTY

       do k=1, number_neigh_processors

         if (send_list_array(k)%length/=0) then  

             call flush_det(send_list_array(k))
  
         end if

       end do

       do k=1, number_neigh_processors

         if (receive_list_array(k)%length/=0) then  

             call flush_det(receive_list_array(k))
 
         end if

       end do

     deallocate(send_list_array)
     deallocate(receive_list_array)

     deallocate(types_det)
      
    end do
    !do_until_all_send_lists_between_processors_are_empty

    call deallocate(ihash) 

    if ((.not.isparallel()).and.(.not. binary_detector_output)) then

       if(getprocno() == 1) then
          if(binary_detector_output) then
            write(detector_unit) time
            write(detector_unit) dt
          else
            format_buffer=reals_format(1)
            write(detector_unit, format_buffer, advance="no") time
            write(detector_unit, format_buffer, advance="no") dt
          end if
       end if

       ! Next columns contain the positions of all the detectors.
      
       node => detector_list%firstnode

       positionloop: do i=1, detector_list%length

       if(getprocno() == 1) then
          if(binary_detector_output) then
            write(detector_unit) node%position
          else
            format_buffer=reals_format(size(node%position))
            write(detector_unit, format_buffer, advance="no") &
                  node%position
          end if
       end if

       node => node%next

       end do positionloop

       phaseloop: do phase=1,size(state)

           do i=1, size(detector_sfield_list(phase)%ptr)
              ! Output statistics for each scalar field.
              sfield=>extract_scalar_field(state(phase), &
                   &                       detector_sfield_list(phase)%ptr(i))

              if(.not. detector_field(sfield)) then
                cycle
              end if
          
              node => detector_list%firstnode

              do j=1, detector_list%length
                value =  detector_value(sfield, node)

                if(getprocno() == 1) then
                
                    if(binary_detector_output) then
                       write(detector_unit) value
                    else
                       format_buffer=reals_format(1)
                       write(detector_unit, format_buffer, advance="no") value
                    end if
                
                 end if

                 node => node%next

              end do
           end do

           allocate(vvalue(0))
 
           do i = 1, size(detector_vfield_list(phase)%ptr)
              ! Output statistics for each vector field
          
              vfield => extract_vector_field(state(phase), &
                   & detector_vfield_list(phase)%ptr(i))
          
              if(.not. detector_field(vfield)) then
                 cycle
              end if

              if (size(vvalue)/=vfield%dim) then
                 deallocate(vvalue)
                 allocate(vvalue(vfield%dim))
              end if

              node => detector_list%firstnode

              do j=1, detector_list%length
            
                 vvalue =  detector_value(vfield, node)
             
                 ! Only the first process should write statistics information
             
              if(getprocno() == 1) then
                
                if(binary_detector_output) then
                  write(detector_unit) vvalue
                else
                  format_buffer=reals_format(vfield%dim)
                  write(detector_unit, format_buffer, advance="no") vvalue
                end if

              end if

              node => node%next

              end do            
        
           end do

            deallocate(vvalue)

       end do phaseloop

       ! Output end of line
       ! Only the first process should write statistics information

       if(getprocno() == 1) then
           if(.not. binary_detector_output) then
             ! Output end of line
             write(detector_unit,'(a)') ""
           end if
           flush(detector_unit)
       end if

    else

       call write_mpi_out(state,time,dt,timestep)

    end if

   contains

    function reals_format(reals)
      character(len=10) :: reals_format
      integer :: reals

      write(reals_format, '(a,i0,a)') '(',reals,'e15.6e3)'

    end function reals_format

  end subroutine write_detectors

  subroutine flush_det(det_list)
  !Removes and deallocates all the nodes in a detector list, starting from the first node.
   
    type(detector_linked_list), intent(inout) :: det_list
    type(detector_type), pointer :: node
    integer :: i

    do i=1, det_list%length

       node => det_list%firstnode

       det_list%firstnode => node%next

       if (associated(det_list%firstnode)) then
    
          det_list%firstnode%previous => null()

       end if

       node%next => null()

       node%previous => null()

       deallocate(node)

       det_list%length = det_list%length-1  

    end do

  end subroutine flush_det

  subroutine write_mpi_out(state,time,dt,timestep)
    !!< Writes detector information (position, value of scalar and vector fields at that position, etc.) into detectors file using MPI output 
    ! commands so that when running in parallel all processors can write at the same time information into the file at the right location.       

    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time, dt
    integer, intent(in) :: timestep

    integer, ALLOCATABLE, DIMENSION(:) :: status

    integer :: i, j, phase, IERROR, nints, number_of_dt, number_of_scalar_det_fields, count, realsize, dimen
    integer(KIND=MPI_OFFSET_KIND) :: location_to_write, offset, offset_to_read
    integer :: number_of_vector_det_fields, number_total_columns

    real, dimension(:), allocatable :: buffer
    real :: value
    real, dimension(:), allocatable :: vvalue
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(detector_type), pointer :: node

    allocate( status(MPI_STATUS_SIZE) )

    node => detector_list%firstnode

    number_of_dt=timestep

    ewrite(1,*) "number of timestep when starting mpi io subroutine:", number_of_dt

    number_of_scalar_det_fields=0

    do phase=1,size(state)
   
        do i=1, size(detector_sfield_list(phase)%ptr)
           ! Output number of detector scalar fields.
           sfield=>extract_scalar_field(state(phase), &
                &                       detector_sfield_list(phase)%ptr(i))   

           if(.not. detector_field(sfield)) then
              cycle
           end if

           number_of_scalar_det_fields=number_of_scalar_det_fields+1

        end do 

    end do 

    number_of_vector_det_fields=0


    do phase=1,size(state)
   
        do i=1, size(detector_vfield_list(phase)%ptr)
           ! Output number of detector vector fields.
           vfield => extract_vector_field(state(phase), &
                & detector_vfield_list(phase)%ptr(i))
       
           if(.not. detector_field(vfield)) then
              cycle

           end if

           number_of_vector_det_fields=number_of_vector_det_fields+1

        end do 

    end do 

    call MPI_TYPE_EXTENT(getpreal(), realsize, ierror)

    vfield => extract_vector_field(state(1),"Velocity")

    dimen=vfield%dim

    number_total_columns=2+total_num_det*dimen+total_num_det*number_of_scalar_det_fields &
                       & +total_num_det*number_of_vector_det_fields*dimen

    ewrite(1,*) "total_num_det is:", total_num_det

    if(have_option("/io/stat/output_at_start")) number_of_dt=number_of_dt+1

    location_to_write=(number_of_dt-1)*number_total_columns*realsize

    if(getprocno() == 1) then

        allocate(buffer(2))
        nints=2
        buffer(1)=time
        buffer(2)=dt

        call MPI_FILE_WRITE_AT(fh,location_to_write,buffer,nints,getpreal(),status, IERROR)

        deallocate(buffer)

    end if
    
!   Offset write to where we start writing detectors
    location_to_write = location_to_write + 2*realsize

    node => detector_list%firstnode

    positionloop: do i=1, detector_list%length

        offset = location_to_write+(node%id_number-1)*size(node%position)*realsize

             if (node%initial_owner==-1) then

             if(getprocno() == 1) then

                 allocate(buffer(size(node%position)))

                 buffer=node%position
                 nints=size(node%position)

                 call MPI_FILE_WRITE_AT(fh,offset,buffer,nints,getpreal(),status,IERROR)

                 deallocate(buffer)

                 node => node%next

             else

             node => node%next

             end if

             else

                 allocate(buffer(size(node%position)))

                 buffer=node%position
                 nints=size(node%position)

                 call MPI_FILE_WRITE_AT(fh,offset,buffer,nints,getpreal(),status,IERROR)

                 deallocate(buffer)

                 node => node%next

             end if

    end do positionloop

    node => detector_list%firstnode

    location_to_write = location_to_write+total_num_det*dimen*realsize

    number_of_scalar_det_fields=0

    phaseloop: do phase=1,size(state)

        do i=1, size(detector_sfield_list(phase)%ptr)
        ! Output statistics for each scalar field.
        sfield=>extract_scalar_field(state(phase), &
        &                       detector_sfield_list(phase)%ptr(i))

        if(.not. detector_field(sfield)) then
          cycle
        end if

        number_of_scalar_det_fields=number_of_scalar_det_fields+1

        node => detector_list%firstnode

        do j=1, detector_list%length

             if (node%initial_owner==-1) then

                if(getprocno() == 1) then

                   offset = location_to_write+(total_num_det*(number_of_scalar_det_fields-1)+(node%id_number-1))*realsize

                   value =  detector_value(sfield, node)

                   allocate(buffer(1))

                   buffer=value
                   nints=1

                   call MPI_FILE_WRITE_AT(fh,offset,buffer,nints,getpreal(),status,IERROR)

                   deallocate(buffer)

                   node => node%next

                   else
                
                   node => node%next

                end if
 
             else

                offset = location_to_write+(total_num_det*(number_of_scalar_det_fields-1)+(node%id_number-1))*realsize

                value =  detector_value(sfield, node)

                allocate(buffer(1))

                buffer=value
                nints=1

                call MPI_FILE_WRITE_AT(fh,offset,buffer,nints,getpreal(),status,IERROR)

                deallocate(buffer)

                node => node%next

             end if

        end do

        end do

        location_to_write = location_to_write+(total_num_det*number_of_scalar_det_fields)*realsize

        number_of_vector_det_fields=0

        allocate(vvalue(0))

        node => detector_list%firstnode

        do i = 1, size(detector_vfield_list(phase)%ptr)
           ! Output statistics for each vector field
       
           vfield => extract_vector_field(state(phase), &
                & detector_vfield_list(phase)%ptr(i))
       
           if(.not. detector_field(vfield)) then
              cycle
           end if

           if (size(vvalue)/=vfield%dim) then
              deallocate(vvalue)
              allocate(vvalue(vfield%dim))
           end if

           number_of_vector_det_fields=number_of_vector_det_fields+1

           node => detector_list%firstnode

           do j=1, detector_list%length

              if (node%initial_owner==-1) then

                 if(getprocno() == 1) then
         
                    vvalue =  detector_value(vfield, node)

                    offset = location_to_write+(node%id_number-1)*size(node%position)*realsize

                    allocate(buffer(size(node%position)))

                    buffer=vvalue
                    nints=size(node%position)

                    call MPI_FILE_WRITE_AT(fh,offset,buffer,nints,getpreal(),status,IERROR)

                    deallocate(buffer)

                    node => node%next

                 else
                
                 node => node%next

                 end if

              else

                 vvalue =  detector_value(vfield, node)
  
                 offset = location_to_write+(node%id_number-1)*size(node%position)*realsize

                 allocate(buffer(size(node%position)))

                 buffer=vvalue
                 nints=size(node%position)

                 call MPI_FILE_WRITE_AT(fh,offset,buffer,nints,getpreal(),status,IERROR)

                 deallocate(buffer)

                 node => node%next

              end if  

           end do            
     
        end do

        deallocate(vvalue)

    end do phaseloop

    call mpi_barrier(mpi_comm_world, ierror)
   
    number_total_columns=2+total_num_det*dimen

    offset_to_read=0

    allocate(buffer(number_total_columns)) 

    call MPI_FILE_READ_AT(fh,offset_to_read,buffer,number_total_columns,getpreal(),status,IERROR)

    call MPI_GET_COUNT(status,getpreal(),count, IERROR)
  
    write(*,*) "count after reading first value is (should be 2):", count

    deallocate(buffer)

    call mpi_barrier(mpi_comm_world, ierror)
   
  end subroutine write_mpi_out

  subroutine move_detectors_bisection_method(state, dt, ihash, send_list_array)
    !!< Move Lagrangian detectors using a bisecting method, i.e., dividing the dt in smaller values that add to dt. 
    !Each smaller value is such that the detector is moved from its previous position to the boundary face between    
    !the element where it belonged to and the neighbouring element in the direction of the flow, through the minimum local 
    !coordinate. Sometimes going through the minimum coordinate is not right (does not follow the flow direction) and in 
    !that case the detector has been placed back into the previous position and previous element and care has been taken for the 
    !detector to follow the flow until the next boundary face
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: dt
    type(integer_hash_table), intent(in) :: ihash
    type(detector_linked_list), dimension(:), allocatable, intent(inout) :: send_list_array

    integer :: j, number_neigh_processors
    real, dimension(:), allocatable :: vel, old_vel, old_pos
    type(vector_field), pointer :: vfield, xfield, old_vfield
    type(detector_type), pointer :: this_det
    integer, dimension(:), pointer :: ele_number_ptr
    real :: dt_temp, dt_var_temp, dt_temp_temp, dt_var_value
    integer :: index_next_face, current_element, previous_element, & 
               & cont_repeated_situation, processor_number

    ewrite(1,*) "Inside move_detectors_bisection_method subroutine"

    vfield => extract_vector_field(state(1),"Velocity")
    old_vfield => extract_vector_field(state(1),"OldVelocity")
    xfield => extract_vector_field(state(1),"Coordinate")

    allocate(vel(vfield%dim))
    allocate(old_vel(old_vfield%dim))
    allocate(old_pos(old_vfield%dim))

    number_neigh_processors=key_count(ihash)
          
       !do_for_each_detector: 

       this_det => detector_list%firstnode

!       initial_length=detector_list%length

       do j=1, detector_list%length
             
          if (this_det%type /= LAGRANGIAN_DETECTOR) then
            
             this_det => this_det%next

             cycle

          end if

          dt_temp=0.0

          previous_element = this_det%element

          cont_repeated_situation=0.0

          !do_until_whole_dt_is_reached or detector leaves the domain towards another processor or through an outflow: 
          do  
         
             if (this_det%dt<1.0e-3*dt) exit
             
             vel =  detector_value(vfield, this_det)
             old_vel = detector_value(old_vfield, this_det)
             old_pos = this_det%position
             current_element = this_det%element

             call move_detectors_subtime_step(this_det, xfield, dt, dt_temp, old_pos, vel, old_vel, vfield, old_vfield,previous_element,index_next_face)

   !From the previous subroutine, the detector is moved until the boundary between two elements is found or all smaller time steps have reached dt. If the detector leaves the domain or we lose track of it, then it will be converted into a static one

             processor_number=getprocno()

             if (this_det%type==STATIC_DETECTOR) exit

             ele_number_ptr=>ele_neigh(xfield,this_det%element)
             previous_element=this_det%element
             this_det%element=ele_number_ptr(index_next_face) 

             processor_number=getprocno()

             call check_if_det_gone_through_domain_boundary(state, this_det, xfield, dt, dt_temp, old_pos, &
                                            vel, old_vel, vfield, old_vfield, index_next_face, current_element, cont_repeated_situation,send_list_array,ihash,processor_number)

             if (processor_number /= getprocno()) exit
     
             if (this_det%type==STATIC_DETECTOR) exit

             this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)
   
             dt_temp=dt_temp+this_det%dt
             this_det%dt=dt-dt_temp;

             if (this_det%dt<1.0e-3*dt) then

                this_det%element=previous_element
                this_det%local_coords=local_coords(xfield,this_det%element,this_det%position) 

             end if

          end do 
          !do_until_whole_dt_is_reached or detector leaves the domain towards another processor or through an outflow: 

             if (processor_number == getprocno()) then 
             !if proc_number /= getprocno() this_det has already been made to point to this_det%next inside 
             !in check_if_det_gone_through_domain_bound.

                this_det => this_det%next

             end if

       end do 
       !do_for_each_detector
   
    deallocate(vel)
    deallocate(old_vel)
    deallocate(old_pos)

  end subroutine move_detectors_bisection_method  

  subroutine serialise_lists_exchange_receive(state,send_list_array,receive_list_array,number_neigh_processors,ihash)
 
    type(state_type), dimension(:), intent(in) :: state
    type(detector_linked_list), dimension(:), intent(inout) :: send_list_array, receive_list_array
    integer, intent(inout) :: number_neigh_processors
    type(integer_hash_table), intent(in) :: ihash

    type array_ptr
       real, dimension(:,:), pointer :: ptr
    end type array_ptr

    type(array_ptr), dimension(:), allocatable :: send_list_array_serialise, receive_list_array_serialise
    type(detector_type), pointer :: node, node_rec
    type(vector_field), pointer :: vfield, xfield
    type(halo_type), pointer :: ele_halo
    type(integer_hash_table) :: gens
    integer :: global_ele, univ_ele, number_detectors_to_send, number_of_columns, total_data_received, number_detectors_received, number_data_to_send, target_proc, mapped_val, count, IERROR, dimen, i, j
    integer, PARAMETER ::TAG=12

    integer, ALLOCATABLE, DIMENSION(:) :: sendRequest
    integer, ALLOCATABLE, DIMENSION(:) :: status
    integer, ALLOCATABLE, DIMENSION(:) :: verification_arr
 
    type(integer_hash_table) :: ihash_inverse
    integer :: halo_level, gensaa, gensaaa, unn, target_proc_a, mapped_val_a
    real :: aabb, aabba, aabbaa, aabbaaa, aabbaaaa

    type(element_type), pointer :: shape

    type(mesh_type) :: pwc_mesh
    type(vector_field) :: pwc_positions

    xfield => extract_vector_field(state(1),"Coordinate")
    shape=>ele_shape(xfield,1)
    vfield => extract_vector_field(state(1),"Velocity")
 
    allocate( sendRequest(0:number_neigh_processors-1) )

    allocate(send_list_array_serialise(number_neigh_processors))

    halo_level = element_halo_count(vfield%mesh)

    if (halo_level /= 0) then

    ele_halo => vfield%mesh%element_halos(halo_level)

    end if

    call allocate(ihash_inverse) 
    do i=1, key_count(ihash)

       call fetch_pair(ihash, i, target_proc_a, mapped_val_a)
       call insert(ihash_inverse, mapped_val_a, target_proc_a)

    end do

    do i=1, key_count(ihash_inverse)

       call fetch_pair(ihash_inverse, i, target_proc_a, mapped_val_a)

    end do

    if (halo_level /= 0) then

    pwc_mesh = piecewise_constant_mesh(xfield%mesh, "PiecewiseConstantMesh")
    call allocate(pwc_positions, xfield%dim, pwc_mesh, "Coordinate")
    call deallocate(pwc_mesh)
    call remap_field(xfield, pwc_positions)
    assert(halo_verifies(ele_halo, pwc_positions))
    call deallocate(pwc_positions)

    end if

    do i=1, number_neigh_processors

       node => send_list_array(i)%firstnode

       dimen=vfield%dim

       number_detectors_to_send=send_list_array(i)%length
       number_of_columns=dimen+3

       allocate(send_list_array_serialise(i)%ptr(number_detectors_to_send,number_of_columns))

       do j=1, send_list_array(i)%length

          global_ele=node%element

          univ_ele = halo_universal_number(ele_halo, global_ele)

          send_list_array_serialise(i)%ptr(j,1:dimen)=node%position
          send_list_array_serialise(i)%ptr(j,dimen+1)=univ_ele
          send_list_array_serialise(i)%ptr(j,dimen+2)=node%dt
          send_list_array_serialise(i)%ptr(j,dimen+3)=node%id_number

          node => node%next
 
       end do

       target_proc=fetch(ihash_inverse, i)

       call MPI_ISEND(send_list_array_serialise(i)%ptr,size(send_list_array_serialise(i)%ptr), &
            & getpreal(), target_proc-1, TAG, MPI_COMM_WORLD, sendRequest(i-1), IERROR)

       !!!getprocno() returns the rank of the processor + 1, hence, for 4 proc, we have 1,2,3,4 whereas 
       !!!the ranks are 0,1,2,3. That is why I am using target_proc-1, so that for proc 4, it sends to 
       !!!proc with rank 3.

       ewrite(1,*) "IERROR is:", IERROR

    end do
  
    allocate(receive_list_array_serialise(number_neigh_processors))

    allocate( status(MPI_STATUS_SIZE) )

    call get_universal_numbering_inverse(ele_halo, gens)

    ewrite(1,*) "gens length is:", key_count(gens)    

    do i=1, key_count(gens)

      call fetch_pair(gens, i, gensaa, gensaaa)
  
    end do

    do i=1, number_neigh_processors
              
       call MPI_PROBE(MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, status(:), IERROR) 

       call MPI_GET_COUNT(status(:), getpreal(), count, IERROR) 

       number_detectors_received=count/number_of_columns

       allocate(receive_list_array_serialise(i)%ptr(number_detectors_received,number_of_columns))

       call MPI_Recv(receive_list_array_serialise(i)%ptr,count, getpreal(), status(MPI_SOURCE), TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERROR)

       do j=1, number_detectors_received

          univ_ele=receive_list_array_serialise(i)%ptr(j,dimen+1);

          global_ele=fetch(gens,univ_ele)

          allocate(node_rec)

          allocate(node_rec%position(dimen))

          node_rec%position=receive_list_array_serialise(i)%ptr(j,1:dimen)
          node_rec%element=global_ele
          node_rec%dt=receive_list_array_serialise(i)%ptr(j,dimen+2)
          node_rec%type = LAGRANGIAN_DETECTOR
          node_rec%local = .true. 
          node_rec%id_number=receive_list_array_serialise(i)%ptr(j,dimen+3)
          node_rec%initial_owner=getprocno()
 
          allocate(node_rec%local_coords(local_coord_count(shape)))          

          node_rec%local_coords=local_coords(xfield,node_rec%element,node_rec%position)

          node_rec%name=name_of_detector_in_read_order(node_rec%id_number)

          call insert_det(receive_list_array(i),node_rec) 
          
       end do

    end do    

    call MPI_WAITALL(number_neigh_processors, sendRequest, MPI_STATUSES_IGNORE, IERROR)

    call deallocate(gens)

    do i=1, number_neigh_processors

       deallocate(send_list_array_serialise(i)%ptr)
       deallocate(receive_list_array_serialise(i)%ptr)

    end do    

    deallocate(send_list_array_serialise)
    deallocate(receive_list_array_serialise)
    call deallocate(ihash_inverse) 
  end subroutine serialise_lists_exchange_receive

  subroutine move_det_from_receive_list_to_det_list(detector_list,receive_list)
   
    type(detector_linked_list), intent(inout) :: detector_list
    type(detector_linked_list), intent(inout) :: receive_list

    type(detector_type), pointer :: node
    integer :: i

    do i=1, receive_list%length

       node => receive_list%firstnode

       receive_list%firstnode => node%next

       if (associated(receive_list%firstnode)) then
    
          receive_list%firstnode%previous => null()

       end if

       call insert_det(detector_list,node) 

       receive_list%length = receive_list%length-1  

    end do

  end subroutine move_det_from_receive_list_to_det_list

  subroutine check_if_det_gone_through_domain_boundary(state, this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, vfield, old_vfield, index_next_face, current_element, cont_repeated_situation,send_list_array,ihash,processor_number)

    type(state_type), dimension(:), intent(in) :: state

    type(detector_type), pointer :: this_det, node_to_send
    type(vector_field), intent(inout) :: xfield
    real, intent(in) :: dt
    real, intent(inout) :: dt_temp
    real,  dimension(:), intent(inout) :: old_pos, vel, old_vel
    type(vector_field), pointer :: vfield, old_vfield
    integer, intent(inout) :: index_next_face, current_element, cont_repeated_situation
    type(detector_linked_list), dimension(:), intent(inout) :: send_list_array
    type(integer_hash_table), intent(in) :: ihash
    integer, intent(inout) :: processor_number

    integer :: i,k,h, univ_ele, halo_level, nhalos, ele, ele_owner, node_owner, ele_owned
    real :: dt_var_temp, dt_temp_temp, dt_var_value
    integer :: old_index_minloc, cont_aaa, cont_bbb, cont_ccc, cont_ddd, cont_eee, cont_fff, & 
               cont_ggg, cont_times_same_element, list_neigh_processor, processor_number_curr_ele
    integer, dimension(:), allocatable :: element_next_to_boundary

    integer, dimension(:), pointer :: nodes

    type(mesh_type), pointer :: mesh_ele

    type(halo_type), pointer :: ele_halo, node_halo 

    cont_aaa=0.0
    cont_bbb=0.0
    cont_ccc=0.0
    cont_ddd=0.0
    cont_eee=0.0
    cont_fff=0.0
    cont_ggg=0.0

    vfield => extract_vector_field(state(1),"Velocity")
    halo_level = element_halo_count(vfield%mesh)

    if (halo_level /= 0) then

    ele_halo => vfield%mesh%element_halos(halo_level)
  
    univ_ele = halo_universal_number(ele_halo, current_element)

    end if
    
    processor_number_curr_ele=element_owner(vfield,current_element)

    cont_times_same_element=0.0

    ele = current_element

    if (this_det%element<0.0) then

       if (element_owned(vfield,current_element)) then

          write(1,*) "CURRENT PROC OWNS the previous element:", current_element 

          !If I own the previous element where the detector was before it left the domain, it means it was not an halo element
          !of another processor, and hence, it has left through a proper boundary (outflow) and it is converted into static.
          !Or it can also be that by some issue in the geometry of the element for example, the detector is not being moved in the right 
          !direction and is leaving through the wrong face.
          !If after some checks done below, the detector is still leaving the domain through that face, it is converted into static.

          !current_element is the previous element where the detector was before moving into a negative element (outside the domain)

          cont_aaa=cont_aaa+1

          cont_repeated_situation=cont_repeated_situation+1
          if (cont_repeated_situation>99999) then
             cont_repeated_situation=1
          end if
       
          allocate(element_next_to_boundary(cont_repeated_situation))

          element_next_to_boundary(cont_repeated_situation)=current_element 

          this_det%element=current_element
          dt_var_value=this_det%dt
          this_det%dt=this_det%dt*0.6

          call update_detector_position_bisect(this_det, xfield, old_pos, vel, old_vel)
        
          do
            
            if ((all(this_det%local_coords>=0.0)).and.(this_det%local_coords(index_next_face)>1.0e-5)) exit

            this_det%dt=this_det%dt/2
          
            call update_detector_position_bisect(this_det, xfield, old_pos, vel, old_vel)

            cont_bbb=cont_bbb+1

            !If a lagrangian detector is initially in a boundary it does not work because the velocity is zero 
            !and gets stuck in this loop
            !Better not to put any detector in a boundary but in case it happens a check has been included so that 
            !the code does not get stuck.
            !The Lagrangian detector would be converted into a static one.

            if ((cont_bbb>1.0e6).or.(this_det%dt<1.0e-5*dt_var_value))  exit

          end do        

          if ((cont_bbb>1.0e6).or.(this_det%dt<1.0e-5*dt_var_value)) then

            cont_fff=cont_fff+1

            this_det%position=old_pos
            this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

            this_det%type=STATIC_DETECTOR
        
          end if 

         !Below a check is included in case the detector keeps going back to the same element next to the boundary. 
         !This means the velocity in that area is not moving the detector inside the domain away from the boundary 
         !element, and the detector continuously ends up in the boundary. In this situation the code will keep 
         !entering here to place the detector inside the same element next to the boundary and the detector 
         !does not progress any more. 

         do h=1, size(element_next_to_boundary)

            if (element_next_to_boundary(h)==current_element) then

              cont_times_same_element=cont_times_same_element+1

            end if

         end do
     
         if ((cont_repeated_situation>50.0).and.(cont_times_same_element>50.0)) then

            cont_ggg=cont_ggg+1

            this_det%position=old_pos
            this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

            this_det%type=STATIC_DETECTOR

         end if

         if (this_det%type==STATIC_DETECTOR) return

         vel =  detector_value(vfield, this_det)
         old_vel = detector_value(old_vfield, this_det)
         old_pos = this_det%position

!Next we start moving a bit the detector after placing it back into the previous element and check that this time it does not go through the same boundary and hence, ends outside the domain. If that is the case, then we make the detector static already now.             

         dt_temp_temp=dt_temp+this_det%dt
         dt_var_temp=dt-dt_temp_temp;

!As part of the RK second order algorithm, we move first the particle or detector using half of the current time step (dt_var_temp)

         this_det%position=(vel+old_vel)/2*dt_var_temp/2+old_pos

         this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

         old_index_minloc=index_next_face

!Next we check if the detector ends up having again negative local coord. through the face it came out of the domain before.
!If so, we check then if another of the local coords. is also negative and in that case it has gone out through the element via another face or boundary, but if no other local coord is also negative, then the detector is ending in the same place as before, that was outside the domain and if this is the case we make it static.

         if (this_det%local_coords(index_next_face)<0.) then

            cont_ccc=cont_ccc+1
           
            do k=1, size(this_det%local_coords)

               if (k /= old_index_minloc) then

                 cont_ddd=cont_ddd+1

                 if (this_det%local_coords(k)<0.) then

                 cont_eee=cont_eee+1

                 index_next_face=k

                 end if

               end if

            end do

            if (index_next_face == old_index_minloc) then

               cont_fff=cont_fff+1

               this_det%position=old_pos
               this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

               this_det%type=STATIC_DETECTOR
        
            end if     
   
         end if 
    
         this_det%position=old_pos
         this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

         deallocate(element_next_to_boundary)

       else

         !If I don't own the previous element where the detector was before it left the domain, it means it was an halo element
         !of another processor, and hence, the detector is sent to the send list associated with that processor.

         !Node in linked_list that contains this detector is removed and incorporated into the send_list 
         !corresponding to the neighbouring processor

         this_det%element=current_element 

         if (halo_level /= 0) then
  
         univ_ele = halo_universal_number(ele_halo, current_element)

         ewrite(1,*) "In check_if_det_gone_through_domain_boundary, else current proc do not own current ele"   
         ewrite(1,*) "detector LOCAL or GLOBAL element now is current element:", this_det%element
         ewrite(1,*) "current_element, I think this is the previous element to the one that is negative:", current_element
         ewrite(1,*) "detector UNIVERSAL element of current element:", univ_ele
         ewrite(1,*) "WE ARE IN PROC.:", getprocno()

         end if

         processor_number=element_owner(vfield,this_det%element)

         this_det%position=old_pos
         this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

         dt_temp=dt_temp+this_det%dt
         this_det%dt=dt-dt_temp

         list_neigh_processor=fetch(ihash,processor_number)

         node_to_send => this_det

         this_det => this_det%next

         call move_det_to_send_list(detector_list,node_to_send,send_list_array(list_neigh_processor))

      end if

    end if

  end subroutine check_if_det_gone_through_domain_boundary

  subroutine move_det_to_send_list(detector_list,node,send_list)

    type(detector_linked_list), intent(inout) :: detector_list
    type(detector_type), pointer :: node
    type(detector_linked_list), intent(inout) :: send_list

    if ((.not.associated(node%previous)).and.(detector_list%length/=1)) then
         !!this checks if the current node that we are going to remove from the list is the first 
         !!one in the list but not the only node in the list

            node%next%previous => null()

            detector_list%firstnode => node%next

            detector_list%firstnode%previous => null()
 
            detector_list%length = detector_list%length-1

     else 

          if ((.not.associated(node%next)).and.(associated(node%previous))) then
          !!this takes into account the case when the node is the last one in the list but not the only one

               node%previous%next => null()

               detector_list%lastnode => node%previous

               detector_list%lastnode%next => null()

               detector_list%length = detector_list%length-1

          else 

               if (detector_list%length==1) then

                  detector_list%firstnode => null()
                  detector_list%lastnode => null()

                  detector_list%length = detector_list%length-1 

              else

                  node%previous%next => node%next
                  node%next%previous => node%previous

                  detector_list%length = detector_list%length-1 

              end if

          end if

     end if

    call insert_det(send_list,node)  

  end subroutine move_det_to_send_list            

  subroutine detector_list_insert(current_list,node)

    type(detector_linked_list), intent(inout) :: current_list
    type(detector_type), pointer :: node

    if (current_list%length == 0) then

      current_list%firstnode => node 
      current_list%lastnode => node 

      current_list%firstnode%previous => null()
      current_list%lastnode%next => null()
      
      current_list%length = 1

    else
  
      node%previous => current_list%lastnode
      current_list%lastnode%next => node
      current_list%lastnode => node

      current_list%lastnode%next => null()
      
      current_list%length = current_list%length+1
    
    end if 

  end subroutine detector_list_insert

  subroutine move_detectors_subtime_step(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, vfield, old_vfield, previous_element,index_next_face)
  !!< Subroutine that makes sure the Lagrangian detector ends up on one of its boundary faces (within tolerance of
  !+/-10.0e-8) up to where it has reached using one of the smaller time steps (sutime step). 
  !Once the detector is on the boundary (within tolerance of +/-10.0e-8), in the subroutine that calls this one, the   
  !detector is asigned to the neighbouring element through that face and again calling this subroutine, the detector is 
  !moved using another smaller time step until the boundary of this new element, always following the direction of the flow. 
  !These steps are repeated until the sum of all the smaller time steps used to go from elemenet to element is equal to the 
  !total time step.
!    type(state_type), dimension(:), intent(in) :: state

    type(detector_type), pointer :: this_det
    type(vector_field), intent(inout) :: xfield
    real, intent(in) :: dt
    real, intent(inout) :: dt_temp
    real,  dimension(:), intent(inout) :: old_pos, vel, old_vel
    type(vector_field), pointer :: vfield, old_vfield
    integer, intent(in) :: previous_element
    integer, intent(inout) :: index_next_face
    real,  dimension(1:size(this_det%local_coords)) :: bbb

    real :: keep_value_this_det_dt, keep_value_this_det_dt_a
    integer :: cont, cont_a, cont_b, cont_c, cont_d, cont_p, cont_r, &
               cont_check, index_minloc_current, cont_check_a, &
               cont_static_det, bound_elem_iteration

    ewrite(1,*) "Inside move_detectors_subtime_step subroutine"

    cont=0
    cont_a=0
    cont_b=0
    cont_c=0
    cont_d=0
    cont_p=0
    cont_r=0
    cont_check=0
    cont_check_a=0

    do
   
       this_det%dt=this_det%dt/2.0
       this_det%position=((vel+old_vel)/2.0)*this_det%dt+old_pos

       this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

       keep_value_this_det_dt=this_det%dt

       cont=0

       do 
   
          if (all(this_det%local_coords>-10.0e-8)) exit

          this_det%dt=this_det%dt/2.0
       
          call update_detector_position_bisect(this_det, xfield, old_pos, vel, old_vel)

          cont=cont+1

          bbb=this_det%local_coords

          !If no matter how much this_det%dt is reduced, still the detector is not inside the new element when start moving it 
          !in the direction of the flow, it is most likely because the detector left the previous element through the wrong 
          !face. The condition of leaving the element is normally through the face with respect to which the local coordinate 
          !has the minimum value. This can lead sometimes to the detector leaving the element through a face that does not intersect 
          !the direction of the flow. In the next lines, the detector is returned to the previous element and previous position and 
          !it is checked that moving it from there in the direction of the flow, the detector stays in the element at least for a specific 
          !value of this_det%dt (smaller time step).
      
          if ((cont>1.0e6).or.(this_det%dt<1.0e-5*keep_value_this_det_dt)) then

             this_det%element=previous_element

             this_det%position=old_pos

             this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

             vel =  detector_value(vfield, this_det)
             old_vel = detector_value(old_vfield, this_det)

             index_minloc_current=minloc(this_det%local_coords, dim=1)

             cont_check=cont_check+1
                 
             this_det%dt=keep_value_this_det_dt

             this_det%position=(vel+old_vel)/2*this_det%dt+old_pos

             this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

             cont_r=0

             cont_check_a=0

             do 

                if (all(this_det%local_coords>-10.0e-8)) exit

                this_det%dt=this_det%dt/2
         
                call update_detector_position_bisect(this_det, xfield, old_pos, vel, old_vel)

                cont_r=cont_r+1

                if ((cont_r>=1.0e6).or.(this_det%dt<1.0e-5*keep_value_this_det_dt)) then

                     cont_check_a=cont_check_a+1

                end if

                if (cont_check_a/=0) exit

             end do

          end if

          if  (cont_check_a/=0) exit

       end do

       vel = detector_value(vfield, this_det)
       old_vel = detector_value(old_vfield, this_det)
       this_det%dt=this_det%dt*2.0
       this_det%position=((vel+old_vel)/2.0)*this_det%dt+old_pos

       keep_value_this_det_dt_a=this_det%dt
       this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

       if (all(this_det%local_coords>-10.0e-8)) exit

       cont_d=cont_d+1

       this_det%dt=this_det%dt/2.0
       this_det%position=old_pos
       this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)
       vel = detector_value(vfield, this_det)
       old_vel = detector_value(old_vfield, this_det)

       if  (cont_check_a/=0) exit

    end do

    !dt_var=this_det%dt

    !In the next if, since the detector seems to have gone through the wrong face and still when returning it to its   
    !previous element and moving it with the flow, the detector does not stay in the element for a specific value of dt_var 
    !(unless it is for dt_var practically zero, so the detector is in its initial position in the element), it could be that 
    !it is going through a vortex. In this case, in order to find the next element to which the detector belongs when moving 
    !in the direction of the flow, all the elements that share a node are checked. The node chosen is the one opposite to the 
    !face with respect to which the local coordinate is the maximum.

    call scenario_det_gone_through_element_vortex(this_det, xfield, dt, dt_temp, keep_value_this_det_dt, old_pos, &
                                         vel, old_vel, vfield, old_vfield, previous_element,index_next_face,cont_check,cont_check_a,bound_elem_iteration)

    !In the next if, since the detector seems to have gone through the wrong face and when returning it to its   
    !previous element and moving it with the flow, the detector stays in the element for a specific value of dt_var,
    !now the detector is moved towards the appropriate face of the element.

    call scenario_det_gone_through_wrong_face(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, index_next_face,cont_check,cont_check_a,bound_elem_iteration,index_minloc_current)
 
    !In the next if, the most straight forward scenario is dealed with, i.e., the detector previously on the boundary of an element, made to move 
    !in the direction of the flow towards the next element across that boundary, belongs to the next element for a particular value of dt_var
    !(next smaller time step).

    call scenario_det_gone_through_right_face(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, index_next_face,cont_check,bound_elem_iteration)
 
    if (bound_elem_iteration>=33554432) then 

       this_det%type=STATIC_DETECTOR
              
    end if     

    cont_static_det=0
    
    if (this_det%type==STATIC_DETECTOR) then

       cont_static_det=cont_static_det+1
          
       this_det%position=old_pos
       this_det%element=previous_element
       this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

    end if

    if (this_det%type==STATIC_DETECTOR) return

  end subroutine move_detectors_subtime_step

  subroutine scenario_det_gone_through_element_vortex(this_det, xfield, dt, dt_temp, keep_value_this_det_dt, old_pos, &
                                          vel, old_vel, vfield, old_vfield, previous_element, index_next_face,cont_check,cont_check_a,bound_elem_iteration)
  !!< Subroutine that makes sure the Lagrangian detector ends up on one of its boundary faces (within tolerance of
  !+/-10.0e-8) up to where it has reached using one of the smaller time steps (sutime step). 
  !Once the detector is on the boundary (within tolerance of +/-10.0e-8), in the subroutine that calls this one, the   
  !detector is asigned to the neighbouring element through that face and again calling this subroutine, the detector is 
  !moved using another smaller time step until the boundary of this new element, always following the direction of the flow. 
  !These steps are repeated until the sum of all the smaller time steps used to go from elemenet to element is equal to the 
  !total time step.
!    type(state_type), dimension(:), intent(in) :: state

    type(detector_type), pointer :: this_det
    type(vector_field), intent(inout) :: xfield
    real, intent(in) :: dt, keep_value_this_det_dt
    real, intent(inout) :: dt_temp
    real,  dimension(:), intent(inout) :: old_pos, vel, old_vel
    type(vector_field), pointer :: vfield, old_vfield
    integer, intent(in) :: previous_element, cont_check, cont_check_a
    integer, intent(inout) :: index_next_face, bound_elem_iteration

    real :: maxvalue
    integer :: i, k, cont_d, index_temp_maxvalue, cont_loop, node, current_element_number, number_of_elem, &
               cont_elem_neg, cont_static_det, det_inside_an_ele

    integer, dimension(:), pointer :: nodes, elements
    type(csr_sparsity), pointer :: nelist

    if ((cont_check /= 0).and.(cont_check_a /= 0)) then

       ewrite(1,*) "Inside scenario_det_gone_through_element_vortex subroutine"

       this_det%position=old_pos
       this_det%local_coords=local_coords(xfield,this_det%element,this_det%position) 
       vel = detector_value(vfield, this_det)
       old_vel = detector_value(old_vfield, this_det)
       
       nelist => extract_nelist(xfield)

       nodes => ele_nodes(xfield, this_det%element) !pointer to the nodes of a given element number this_det%element

       maxvalue=0.0

       do k=1, size(this_det%local_coords)

          if (this_det%local_coords(k)>maxvalue) then

             maxvalue=this_det%local_coords(k)

             index_temp_maxvalue=k

          end if

       end do  

       current_element_number=this_det%element

       node=nodes(index_temp_maxvalue) !node number of the node that is opposite to the face with respect to which the local 
                                       !coodinate has the highest value
       elements => row_m_ptr(nelist, node) 
 
       !pointer to the row number corresponding to the node number that returns the index of 
       !the columns that are different than zero. These indexes are the numbers of the elements 
       !that share/contain that node number
  
       number_of_elem=size(elements)

       
       !In the loop below, it is checked if the detector belongs to any of the elements that share the node
      
       do i=1, size(elements)
        
          !if (elements(i)==current_element_number) cycle 

          !If one of the elements is negative

          cont_elem_neg=0

          if ((elements(i)<0.0)) then

             cont_elem_neg=cont_elem_neg+1        

          end if

          if (cont_elem_neg/=0) cycle

          this_det%element=elements(i)

          this_det%dt=keep_value_this_det_dt

          cont_d=0

          do
  
             this_det%dt=this_det%dt/2

             this_det%position=(vel+old_vel)/2*this_det%dt+old_pos

             this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

             cont_loop=0

             do 

                if (all(this_det%local_coords>-10.0e-8)) exit

                this_det%dt=this_det%dt/2
       
                call update_detector_position_bisect(this_det, xfield, old_pos, vel, old_vel)

                cont_loop=cont_loop+1

                if ((cont_loop>1.0e6).or.(this_det%dt<1.0e-5*keep_value_this_det_dt)) exit

             end do

             if ((cont_loop>1.0e6).or.(this_det%dt<1.0e-5*keep_value_this_det_dt)) exit

             vel =  detector_value(vfield, this_det)
             old_vel = detector_value(old_vfield, this_det)
             this_det%dt=this_det%dt*2
             this_det%position=(vel+old_vel)/2*this_det%dt+old_pos

             this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

             if (all(this_det%local_coords>-10.0e-8)) exit

             cont_d=cont_d+1

             this_det%dt=this_det%dt/2
             this_det%position=old_pos
             this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)
             vel =  detector_value(vfield, this_det)
             old_vel = detector_value(old_vfield, this_det)

          end do

          if (all(this_det%local_coords>-10.0e-8)) exit

       end do

       det_inside_an_ele=0.0

       if  (all(this_det%local_coords>-10.0e-8)) then

           det_inside_an_ele=det_inside_an_ele+1

       end if

       if ((i==size(elements)).and.(det_inside_an_ele==0.0)) then

           this_det%type=STATIC_DETECTOR

       end if

       cont_static_det=0
    
       if (this_det%type==STATIC_DETECTOR) then

          cont_static_det=cont_static_det+1
          
          this_det%position=old_pos
          this_det%element=previous_element
          this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

       end if

       if (this_det%type==STATIC_DETECTOR) return

       call placing_det_in_boundary_between_elem(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, index_next_face, bound_elem_iteration)
  
    end if

  end subroutine scenario_det_gone_through_element_vortex

  subroutine placing_det_in_boundary_between_elem(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, index_next_face, bound_elem_iteration)

    type(detector_type), pointer :: this_det
    type(vector_field), intent(inout) :: xfield
    real, intent(in) :: dt
    real, intent(inout) :: dt_temp
    real,  dimension(:), intent(inout) :: old_pos, vel, old_vel
    integer, intent(inout) :: index_next_face, bound_elem_iteration
    
    real :: dt_var_check, dt_temp_check, dt_a, dt_b
    integer :: cont_p, cont_ready_out

    bound_elem_iteration=0

       do

       !dt_a is the value of the bisected dt with which the detector falls inside the element for first time

          dt_temp_check=dt_temp+this_det%dt
          dt_var_check=dt-dt_temp_check

          index_next_face=minloc(this_det%local_coords, dim=1)

          cont_ready_out=0

          cont_p=0

          !Next it is checked if the detector is on one of the element faces (within tolerance +/-10.0e-8). 

          !If not, it needs to be moved further in the direction of the flow until it hits the element face.  
  
          !Other conditions where placed before regarding if the ratio between the other local coordinates after and before moving 
          !the detector backwards was less than 1, i.e., the other local coords are generally increasing as we aproach the right 
          !face (with respect to which the local coord is minimum).

          !However there was some ambiguity when using those other conditions as well since for a few geometric cases it was not 
          !satisfied and some detectors were converted into static in the middle of the domain. Hence, they were removed. If the 
          !detector leaves through the wrong face, it is taken care of in other parts of this subroutine, making the detector
          !to go back to the previous element.

          if ((minval(this_det%local_coords)<10.0e-8).and.(minval(this_det%local_coords)>-10.0e-8)) then

                cont_ready_out=cont_ready_out+1
           
          end if 


          if (cont_ready_out /= 0)  index_next_face=minloc(this_det%local_coords, dim=1)
   
          if (cont_ready_out /= 0)  exit   

          !If dt_var_check is practically zero means that all the smaller dt_var used add up to dt so it is the final position
          !of the detector for that dt

          if (dt_var_check<1.0e-3*dt) exit 

          if (bound_elem_iteration>=33554432) exit

          dt_a=this_det%dt
          dt_b=2*dt_a

          !If not yet on one of the faces of the element (within tolerance), it needs to be moved further 
          !in the direction of the flow until it hits the element face.

          do

             !This loop terminates when the detector is on one of the element faces (within tolerance +/-10.0e-8). 

             !When that happens, the detector is moved backwards a bit, to check that the ratio of the local coordinates 
             !associated to the min local coordinate before moving backwards the detector, is not equal to 1. Same explanation 
             !as before
         
             if ((minval(this_det%local_coords)<10.0e-8).and.(minval(this_det%local_coords)>-10.0e-8)) then

                 cont_ready_out=cont_ready_out+1 

             end if 

             if (cont_ready_out /= 0)  exit  

             if (dt_var_check<1.0e-3*dt) exit 

             if (dt_b-dt_a<1.0e-8) exit

             if (bound_elem_iteration>=33554432) exit

             call iterating_for_det_in_bound_elem(this_det, xfield, old_pos, vel, old_vel, dt_a, dt_b, bound_elem_iteration)

          end do

          cont_p=cont_p+1 

       end do

  end subroutine placing_det_in_boundary_between_elem

  subroutine iterating_for_det_in_bound_elem(this_det, xfield, old_pos, vel, old_vel, dt_a, dt_b, bound_elem_iteration)

    type(detector_type), pointer :: this_det
    type(vector_field), intent(inout) :: xfield
    real, intent(inout) :: dt_a, dt_b
    real,  dimension(:), intent(inout) :: old_pos, vel, old_vel
    integer, intent(inout) :: bound_elem_iteration
    
            
    integer :: cont_a, cont_b, cont_c

       ewrite(1,*) "Inside iterating_for_det_in_bound_elem subroutine"
       

       cont_a=0
       cont_b=0
       cont_c=0

       this_det%dt=dt_a+((dt_b-dt_a)/2)

       bound_elem_iteration=2

       cont_a=cont_a+1
  
       do 
              
          call update_detector_position_bisect(this_det, xfield, old_pos, vel, old_vel)

          if (all(this_det%local_coords>-10.0e-8)) exit 

          bound_elem_iteration=bound_elem_iteration*2

          this_det%dt=dt_a+((dt_b-dt_a)/bound_elem_iteration)

          cont_b=cont_b+1

          !if i becomes too big is because no matter how much dt_var is reduced the detector is not found inside the element

          if (bound_elem_iteration>=33554432) exit

       end do

       ewrite(1,*) "cont_b:", cont_b

       ewrite(1,*) "bound_elem_iteration:", bound_elem_iteration

       ewrite(1,*) "dt_a:", dt_a

       ewrite(1,*) "dt_b:", dt_b

       dt_b=dt_a+((dt_b-dt_a)/(bound_elem_iteration/2))
       dt_a=this_det%dt

       ewrite(1,*) "dt_a:", dt_a

       ewrite(1,*) "dt_b:", dt_b

       cont_c=cont_c+1  

  end subroutine iterating_for_det_in_bound_elem

  subroutine scenario_det_gone_through_wrong_face(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, index_next_face,cont_check,cont_check_a,bound_elem_iteration, index_minloc_current)

    type(detector_type), pointer :: this_det
    type(vector_field), intent(inout) :: xfield
    real, intent(in) :: dt
    real, intent(inout) :: dt_temp
    real,  dimension(:), intent(inout) :: old_pos, vel, old_vel
    integer, intent(in) :: cont_check, cont_check_a
    integer, intent(inout) :: index_next_face, bound_elem_iteration, index_minloc_current

    real :: dt_var_check, dt_temp_check, dt_a, dt_b, minvalue
    integer :: k, cont_p, &
               cont_apq, index_temp, cont_ready

    if ((cont_check /= 0).and.(cont_check_a==0)) then

       ewrite(1,*) "Inside scenario_det_gone_through_wrong_face subroutine"

       bound_elem_iteration=0

       cont_ready=0

       do

          dt_temp_check=dt_temp+this_det%dt
          dt_var_check=dt-dt_temp_check
    
          !This loop terminates when the detector is on one of the element faces (within tolerance +/-10.0e-8) and the min local
          !coordinate occurs with respect to a different face than before when the detector left through the wrong face. 

          !if ((minval(this_det%local_coords)<10.0e-8).and.(minval(this_det%local_coords)>-10.0e-8).and.(minloc(this_det%local_coords, dim=1)/=index_minloc_current)) exit 

          !In the next if, if the min local coords still occurs in the same index as before, then it is checked if 
          !there is another local coord that is also within tolerance and the face associated with it is the new face through 
          !which the detector will leave the element.


           if (cont_ready /= 0) exit

           if (dt_var_check<1.0e-3*dt) exit 

          !The check below means we have lost the detector at some point. We dont know what has happened so we exit the loop and at the end of 
          !the subroutine we convert the detector into a static one.

          if (bound_elem_iteration>=33554432) exit

          index_next_face=minloc(this_det%local_coords, dim=1)

          cont_apq=0

          dt_a=this_det%dt
          dt_b=2*dt_a

          !If the detector is not on one of the element faces (within tolerance +/-10.0e-8), then it needs to keep moving in the
          !direction of the flow until it hits a face.

!          cont_cont=0

          do

          !As before this loop terminates when the detector is on one of the element faces (within tolerance +/-10.0e-8) 
          !and the min local coordinate occurs with respect to a different face than before when the detector left through 
          !the wrong face. Same explanations as before apply
 
             cont_ready=0

          !This loop terminates when the detector is on one of the element faces (within tolerance +/-10.0e-8) and the min local
          !coordinate occurs with respect to a different face than before when the detector left through the wrong face. 

             if ((minval(this_det%local_coords)<10.0e-8).and.(minval(this_det%local_coords)>-10.0e-8).and.(minloc(this_det%local_coords, dim=1)/=index_minloc_current)) exit 


          !In the next if, if the min local coords still occurs in the same index as before, then it is checked if 
          !there is another local coord that is also within tolerance and the face associated with it is the new face through 
          !which the detector will leave the element.

             if ((minval(this_det%local_coords)<10.0e-8).and.(minval(this_det%local_coords)>-10.0e-8).and.(minloc(this_det%local_coords, dim=1)==index_minloc_current)) then

                 minvalue=1000.0

                 do k=1, size(this_det%local_coords)

                    if (k /= index_minloc_current) then

                        if (this_det%local_coords(k)<minvalue) then

                            minvalue=this_det%local_coords(k)

                            index_temp=k

                            ewrite(1,*) "index_temp", index_temp

                        end if

                    end if

                 end do  

                 cont_apq=cont_apq+1

                 cont_ready=0

                 if ((this_det%local_coords(index_temp)<10.0e-8).and.(this_det%local_coords(index_temp)>-10.0e-8))  then

                    cont_ready=cont_ready+1

                    index_next_face=index_temp
              
                 end if

                 ewrite(1,*) "index_temp", index_temp

                 ewrite(1,*) "index_next_face", index_next_face

             end if
 
             if (cont_ready /= 0) exit
          
             dt_temp_check=dt_temp+this_det%dt
             dt_var_check=dt-dt_temp_check
               
             if (dt_var_check<1.0e-3*dt) exit 

             if (dt_b-dt_a<1.0e-8) exit

             if (bound_elem_iteration>=33554432) exit

             call iterating_for_det_in_bound_elem(this_det, xfield, old_pos, vel, old_vel, dt_a, dt_b, bound_elem_iteration)

          end do

          cont_p=cont_p+1 

       end do

    end if
  
  end subroutine scenario_det_gone_through_wrong_face

  subroutine scenario_det_gone_through_right_face(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, index_next_face, cont_check,bound_elem_iteration)

    type(detector_type), pointer :: this_det
    type(vector_field), intent(inout) :: xfield
    real, intent(in) :: dt
    real, intent(inout) :: dt_temp
    real,  dimension(:), intent(inout) :: old_pos, vel, old_vel
!    type(vector_field), pointer :: vfield, old_vfield
    integer, intent(in) :: cont_check
    integer, intent(inout) :: index_next_face, bound_elem_iteration
    
    if (cont_check == 0) then

       ewrite(1,*) "Inside scenario_det_gone_through_right_face subroutine"

       call placing_det_in_boundary_between_elem(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, index_next_face, bound_elem_iteration)
      
    end if 

  end subroutine scenario_det_gone_through_right_face

  subroutine update_detector_position_bisect(this_det, xfield, old_pos, vel, old_vel)
    !!< Moves the detector in the direction of the flow from an initial position (old_pos) during a small time step equal to dt_var 
    !(a bisection of the total time step) and updates the local coordinates 
!    type(state_type), dimension(:), intent(in) :: state

    type(detector_type), pointer :: this_det
    type(vector_field), intent(inout) :: xfield
    real,  dimension(:), intent(in) :: old_pos, vel, old_vel

       this_det%position=((vel+old_vel)/2.0)*this_det%dt+old_pos

       this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

  end subroutine update_detector_position_bisect

 function detector_value_scalar(sfield, detector) result(value)
    !!< Evaluate field at the location of the detector.
    real :: value
    type(scalar_field), intent(in) :: sfield
    type(detector_type), intent(in) :: detector

    value=0.0
    
    if(detector%element>0) then
       if(detector%element > 0) then
         value = eval_field(detector%element, sfield, detector%local_coords)
       end if
    end if

    if (.not. detector%local) call allsum(value)

  end function detector_value_scalar

  function detector_value_vector(vfield, detector) result(value)
    !!< Evaluate field at the location of the detector.
    type(vector_field), intent(in) :: vfield
    type(detector_type), intent(in) :: detector
    real, dimension(vfield%dim) :: value

    value=0.0
    
    if(detector%element>0) then
      if(detector%element > 0) then
        value = eval_field(detector%element, vfield, detector%local_coords)
      end if
    end if

    if(.not. detector%local) call allsum(value)

  end function detector_value_vector

  subroutine set_detector_coords_from_python(values, ndete, func, time)
    !!< Given a list of positions and a time, evaluate the python function
    !!< specified in the string func at those points. 
    real, dimension(:,:), target, intent(inout) :: values
    !! Func may contain any python at all but the following function must
    !! be defiled:
    !!  def val(t)
    !! where t is the time. The result must be a float. 
    character(len=*), intent(in) :: func
    real :: time
    
    real, dimension(:), pointer :: lvx,lvy,lvz
    real, dimension(0), target :: zero
    integer :: stat, dim, ndete

    call get_option("/geometry/dimension",dim)

    lvx=>values(1,:)
    lvy=>zero
    lvz=>zero
    if(dim>1) then
       lvy=>values(2,:)
       if(dim>2) then
          lvz => values(3,:)
       end if
    end if

    call set_detectors_from_python(func, len(func), dim, &
         ndete, time, dim,                               &
         lvx, lvy, lvz, stat)

    if (stat/=0) then
      ewrite(-1, *) "Python error, Python string was:"
      ewrite(-1 , *) trim(func)
      FLAbort("Dying")
    end if

  end subroutine set_detector_coords_from_python
    
  subroutine close_diagnostic_files()
    !! Closes .stat, .convergence and .detector file (if openened)
    !! Gives a warning for iostat/=0, no point to flabort though.

    integer:: stat, IERROR

    if (diag_unit/=0) then
       close(diag_unit, iostat=stat)
       if (stat/=0) then
          ewrite(0,*) "Warning: failed to close .stat file"
       end if
    end if

    if (conv_unit/=0) then
       close(conv_unit, iostat=stat)
       if (stat/=0) then
          ewrite(0,*) "Warning: failed to close .convergence file"
       end if
    end if

    if(steady_state_unit /= 0) then
      close(steady_state_unit, iostat = stat)
      if(stat /= 0) then
        ewrite(0, *) "Warning: failed to close .steady_state file"
      end if
    end if

    if (detector_unit/=0) then
       close(detector_unit, iostat=stat)
       if (stat/=0) then
          ewrite(0,*) "Warning: failed to close .detector file"
       end if
    end if

    if (fh/=0) then
       call MPI_FILE_CLOSE(fh, IERROR) 
       if (IERROR/=0) then
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
      'contents are in existance'
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
