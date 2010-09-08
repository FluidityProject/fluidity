#include "fdebug.h"

program periodise

  !! Take in an flml with an external mesh and one periodic mesh derived from it.
  !! Make it so that the periodic mesh is external, and the previously external mesh
  !! is derived from it instead.

  use populate_state_module
  use state_module
  use fields
  use reference_counting
  use spud
  use field_options
  use mesh_files
  use integer_set_module
  implicit none

  character(len=4096) :: filename, external_filename, new_external_filename, new_filename
  integer :: status
  type(state_type), dimension(:), pointer :: states
  integer :: ierr
  character(len=FIELD_NAME_LEN) :: external_name, periodic_name
  type(mesh_type), pointer :: periodic_mesh, external_mesh
  type(vector_field) :: periodic_positions, external_positions
  integer :: stat, i, nstates
  logical :: skip_initial_extrusion

  call set_debug_level(0)
  call mpi_init(ierr)
  call python_init

  call get_command_argument(1, value=filename, status=status)
  if (status > 0) then
    call usage
    stop
  else if (status < 0) then
    write(0,*) "Warning: truncating filename"
  end if

  call load_options(trim(filename))
  
  ! for extruded meshes, if no checkpointed extruded mesh is present, don't bother
  ! extruding (this may be time consuming or not fit on the input_nprocs)
  skip_initial_extrusion = option_count('/geometry/mesh/from_mesh/extrude')>0 .and. & 
    option_count('/geometry/mesh/from_mesh/extrude/checkpoint_from_file')==0
        
  ! ! Below is a (partial) copy of the first bit of populate_state
  
  ! Find out how many states there are
  nstates=option_count("/material_phase")
  allocate(states(1:nstates))
  do i = 1, nstates
     call nullify(states(i))
  end do

  call insert_external_mesh(states, save_vtk_cache = .true.)
  
  call insert_derived_meshes(states, skip_extrusion=skip_initial_extrusion)
  
  ! !  End populate_state calls
  
  call check_valid_input(states, external_filename, external_name, periodic_name)

  external_mesh => extract_mesh(states(1), trim(external_name))
  periodic_mesh => extract_mesh(states(1), trim(periodic_name))
  if (trim(external_name) /= "CoordinateMesh") then
    external_positions = extract_vector_field(states(1), trim(external_name) // 'Coordinate')
  else
    external_positions = extract_vector_field(states(1), 'Coordinate')
  end if
  call incref(external_positions)
  call allocate(periodic_positions, external_positions%dim, periodic_mesh, trim(periodic_name) // 'Coordinate')
  call remap_field(external_positions, periodic_positions, stat=stat)
  if(stat==REMAP_ERR_DISCONTINUOUS_CONTINUOUS) then
    FLAbort("Just remapped from a discontinuous to a continuous field!")
  else if(stat==REMAP_ERR_HIGHER_LOWER_CONTINUOUS) then
    FLAbort("Just remapped from a higher order to a lower order continuous field!")
  end if
  ! we've allowed it to remap from periodic to unperiodic
  
  call postprocess_periodic_mesh(external_mesh, external_positions, periodic_mesh, periodic_positions)

  ! Dump out the periodic mesh to disk:
  new_external_filename = trim(external_filename) // '_periodic'
  call write_mesh_files(new_external_filename, periodic_positions)

  ! OK! Now we need to do some setting of options.
  call manipulate_options(external_mesh, trim(external_mesh%option_path), periodic_mesh, trim(periodic_mesh%option_path), new_external_filename)

  new_filename = filename(1:len_trim(filename)-5) // '_periodised.flml'
  call write_options(new_filename)

  call deallocate(states)
  call deallocate(external_positions)
  call deallocate(periodic_positions)
  call print_references(1)
  call python_end
  call mpi_finalize(ierr)

  contains

  subroutine usage
    write(0,*) "Usage: periodise input.flml"
    write(0,*) " where input.flml has exactly one external mesh and one periodic mesh"
    write(0,*) " derived immediately from it. The tool produces input_periodised.flml"
    write(0,*) " which has the periodic mesh as the external one, and the previously"
    write(0,*) " external mesh derived from it."
  end subroutine usage

  subroutine check_valid_input(states, external_filename, external_name, periodic_name)
    type(state_type), dimension(:), intent(in) :: states
    character(len=FIELD_NAME_LEN), intent(out) :: external_name, periodic_name
    integer :: j
    integer :: seen_external_meshes, seen_periodic_meshes
    character(len=4096), intent(out) :: external_filename
    type(mesh_type), pointer :: mesh

    seen_external_meshes = 0
    seen_periodic_meshes = 0

    do j=1,mesh_count(states(1))
      mesh => extract_mesh(states(1), j)

      if (have_option(trim(mesh%option_path) // '/from_file')) then
        seen_external_meshes = seen_external_meshes + 1
        external_name = mesh%name
        call get_option(trim(mesh%option_path) // '/from_file/file_name', external_filename)
      end if

      if (have_option(trim(mesh%option_path) // '/from_mesh/periodic_boundary_conditions') .and. &
          .not. have_option(trim(mesh%option_path) // '/from_mesh/periodic_boundary_conditions/remove_periodicity')) then
        seen_periodic_meshes = seen_periodic_meshes + 1
        periodic_name = mesh%name
      end if
    end do

    if (seen_external_meshes /= 1 .or. seen_periodic_meshes /= 1) then
      call usage
      stop
    end if
  end subroutine check_valid_input

  subroutine manipulate_options(external_mesh, external_path, periodic_mesh, periodic_path, new_external_filename)
    type(mesh_type), intent(in) :: external_mesh, periodic_mesh
    character(len=*), intent(in) :: new_external_filename
    character(len=*), intent(in) :: external_path, periodic_path
    character(len=8192) :: str
    integer :: stat, periodic_bc
    integer, dimension(2) :: shape_option
    integer, dimension(:), allocatable :: boundary_ids

    ! The periodic sub-branch of the options tree has useful information that
    ! we don't want to lose just yet. So set the external mesh to be periodic
    ! first, and then overwrite

    call delete_option(external_path // '/from_file', stat=stat)
    call set_option_attribute(external_path // '/from_mesh/mesh/name', trim(periodic_mesh%name), stat=stat)
    do periodic_bc=0,option_count(periodic_path // '/from_mesh/periodic_boundary_conditions')-1

      call get_option(periodic_path // '/from_mesh/periodic_boundary_conditions[' // int2str(periodic_bc) // ']/name', str)
      call set_option_attribute(external_path // '/from_mesh/periodic_boundary_conditions[' // int2str(periodic_bc) // ']/name', trim(str), stat=stat)

      call add_option(external_path // '/from_mesh/periodic_boundary_conditions[' // int2str(periodic_bc) // ']/remove_periodicity', stat=stat)

      call get_option(periodic_path // '/from_mesh/periodic_boundary_conditions[' // int2str(periodic_bc) // ']/coordinate_map', str)
      call set_option(external_path // '/from_mesh/periodic_boundary_conditions[' // int2str(periodic_bc) // ']/coordinate_map', trim(str), stat=stat)
      if (have_option(periodic_path // '/from_mesh/periodic_boundary_conditions[' // int2str(periodic_bc) // ']/inverse_coordinate_map')) then
        call get_option(periodic_path // '/from_mesh/periodic_boundary_conditions[' // int2str(periodic_bc) // ']/inverse_coordinate_map', str)
        call set_option(external_path // '/from_mesh/periodic_boundary_conditions[' // int2str(periodic_bc) // ']/inverse_coordinate_map', trim(str), stat=stat)
      end if

      shape_option = option_shape(periodic_path // '/from_mesh/periodic_boundary_conditions['//int2str(periodic_bc)//']/aliased_boundary_ids')
      allocate(boundary_ids(shape_option(1)))
      call get_option(periodic_path // '/from_mesh/periodic_boundary_conditions['//int2str(periodic_bc)//']/aliased_boundary_ids', boundary_ids)
      call set_option(external_path // '/from_mesh/periodic_boundary_conditions['//int2str(periodic_bc)//']/aliased_boundary_ids', boundary_ids, stat=stat)
      deallocate(boundary_ids)

      shape_option = option_shape(periodic_path // '/from_mesh/periodic_boundary_conditions['//int2str(periodic_bc)//']/physical_boundary_ids')
      allocate(boundary_ids(shape_option(1)))
      call get_option(periodic_path // '/from_mesh/periodic_boundary_conditions['//int2str(periodic_bc)//']/physical_boundary_ids', boundary_ids)
      call set_option(external_path // '/from_mesh/periodic_boundary_conditions['//int2str(periodic_bc)//']/physical_boundary_ids', boundary_ids, stat=stat)
      deallocate(boundary_ids)
    end do

    call delete_option(periodic_path // '/from_mesh', stat=stat)
    call set_option_attribute(periodic_path // '/from_file/file_name', new_external_filename, stat=stat)
    call set_option_attribute(periodic_path // '/from_file/format/name', 'triangle', stat=stat)
    call add_option(periodic_path // '/from_file/stat/include_in_stat', stat=stat)
    call add_option(external_path // '/from_mesh/stat/include_in_stat', stat=stat)
  end subroutine manipulate_options

end program periodise
