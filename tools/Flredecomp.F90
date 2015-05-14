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
#include "confdefs.h"

subroutine flredecomp(input_basename, input_basename_len, output_basename, output_basename_len, &
  & input_nprocs, target_nprocs) bind(c)
  !!< Peform a redecomposition of an input checkpoint with input_nprocs
  !!< processes to a new checkpoint with target_nprocs processes.
  
  use checkpoint
  use fldebug
  use global_parameters, only: is_active_process, no_active_processes, topology_mesh_name
  use parallel_tools
  use populate_state_module
  use spud
  use sam_integration
  use fields
#ifdef HAVE_ZOLTAN
  use zoltan
#endif
  use zoltan_integration
  use state_module
  use initialise_ocean_forcing_module
  use iso_c_binding

  implicit none

  character(kind=c_char, len=1) :: input_basename(*)
  integer(kind=c_size_t), value :: input_basename_len
  character(kind=c_char, len=1) :: output_basename(*)
  integer(kind=c_size_t), value :: output_basename_len
  integer(kind=c_int), value :: input_nprocs
  integer(kind=c_int), value :: target_nprocs
  
  interface
    subroutine check_options()
    end subroutine check_options

#ifdef HAVE_PYTHON
    subroutine python_init()
    end subroutine python_init
#endif
  end interface
  
  character(len=input_basename_len):: input_base
  character(len=output_basename_len):: output_base
  integer :: nprocs
  type(state_type), dimension(:), pointer :: state
  type(vector_field) :: extruded_position
  logical :: skip_initial_extrusion
  integer :: i, nstates
#ifdef HAVE_ZOLTAN
  real(zoltan_float) :: ver
  integer(zoltan_int) :: ierr

  ierr = Zoltan_Initialize(ver)  
  assert(ierr == ZOLTAN_OK)
#endif
  
  ewrite(1, *) "In flredecomp"

#ifdef HAVE_PYTHON
  call python_init()
#endif
  
  nprocs = getnprocs()
  ! now turn into proper fortran strings (is there an easier way to do this?)
  do i=1, input_basename_len
    input_base(i:i)=input_basename(i)
  end do
  do i=1, output_basename_len
    output_base(i:i)=output_basename(i)
  end do
  
  ewrite(2, "(a)") "Input base name: " // trim(input_base)
  ewrite(2, "(a)") "Output base name: " // trim(output_base)
  ewrite(2, "(a,i0)") "Input number of processes: ", input_nprocs
  ewrite(2, "(a,i0)") "Target number of processes: ", target_nprocs
  ewrite(2, "(a,i0)") "Job number of processes: ", nprocs
  
  ! Input check
  if(input_nprocs < 0) then
    FLExit("Input number of processes cannot be negative!")
  else if(target_nprocs < 0) then
    FLExit("Target number of processes cannot be negative!")
  else if(input_nprocs > nprocs) then
    ewrite(-1, *) "The input number of processes must be equal or less than the number of processes currently running."
    FLExit("Running on insufficient processes.")
  else if(target_nprocs > nprocs) then
    ewrite(-1, *) "The target number of processes must be equal or less than the number of processes currently running."
    FLExit("Running on insufficient processes.")
  end if
  
  ! Load the options tree
  call load_options(trim(input_base) // ".flml")
  if(.not. have_option("/simulation_name")) then
    FLExit("Failed to find simulation name after loading options file")
  end if

  if(debug_level() >= 1) then
    ewrite(1, *) "Options tree:"
    call print_options()
  end if

#ifdef DDEBUG
  ewrite(1, *) "Performing options sanity check"
  call check_options()
  ewrite(1, *) "Options sanity check successful"
#endif

  
  ! for extruded meshes, if no checkpointed extruded mesh is present, don't bother
  ! extruding (this may be time consuming or not fit on the input_nprocs)
  skip_initial_extrusion = option_count('/geometry/mesh/from_mesh/extrude')>0 .and. & 
    option_count('/geometry/mesh/from_mesh/extrude/checkpoint_from_file')==0
        
  is_active_process = getprocno() <= input_nprocs
  no_active_processes = input_nprocs
  
  ! ! Below is a (partial) copy of the first bit of populate_state
  
  ! Find out how many states there are
  nstates=option_count("/material_phase")
  allocate(state(1:nstates))
  do i = 1, nstates
     call nullify(state(i))
  end do

  call initialise_ocean_forcing_readers
  
  call insert_external_mesh(state, save_vtk_cache = .true.)
  
  call insert_derived_meshes(state, skip_extrusion=skip_initial_extrusion)

  call compute_domain_statistics(state)

  call allocate_and_insert_fields(state)

  call initialise_prognostic_fields(state, save_vtk_cache=.true., &
    initial_mesh=.true.)

  call set_prescribed_field_values(state, initial_mesh=.true.)
  
  ! !  End populate_state calls
    
  is_active_process = .true.
  no_active_processes = target_nprocs
  
#ifdef HAVE_ZOLTAN
  call zoltan_drive(state, .true., initialise_fields=.true., ignore_extrusion=skip_initial_extrusion, &
     & flredecomping=.true., input_procs = input_nprocs, target_procs = target_nprocs)
#else
  call strip_level_2_halo(state, initialise_fields=.true.)
  call sam_drive(state, sam_options(target_nprocs))
#endif
  
  ! Output
  assert(associated(state))
  call checkpoint_simulation(state, prefix = output_base, postfix = "", protect_simulation_name = .false., &
    keep_initial_data=.true., ignore_detectors=.true., number_of_partitions=target_nprocs)

  do i = 1, size(state)
    call deallocate(state(i))
  end do
    
  ewrite(1, *) "Exiting flredecomp"
  
contains

  function sam_options(target_nparts)
    !!< Return sam options array
    
    integer, intent(in) :: target_nparts

    integer, dimension(10) :: sam_options
    
    sam_options = 0
    
    ! Target number of partitions - 0 indicates size of MPI_COMM_FEMTOOLS
    sam_options(1) = target_nparts

    ! Graph partitioning options:
    sam_options(2) = 1    ! Clean partitioning to optimise the length of the 
                          ! interface boundary.
    ! sam_options(2) = 2  ! Local diffusion
    ! sam_options(2) = 3  ! Directed diffusion
    ! sam_options(2) = 4  ! Clean partitioning to optimise the length of the 
                          ! interface boundary. This partitioning is then remapped
                          ! onto the original partitioning to maximise overlap and
                          ! therefore the volume of data migration.

    ! Heterogerious options (disabled)
    sam_options(3) = 1
    ! No node weights
    sam_options(4) = 1
    ! No edge weights
    sam_options(5) = 1
    ! Mixed formulation options
    sam_options(6) = 2 ! Enabled
                       ! Restore the level 2 halo

  end function sam_options
  
end subroutine flredecomp  
