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
  use write_triangle
  implicit none

  character(len=4096) :: filename, external_filename, new_filename
  integer :: status
  type(state_type), dimension(:), pointer :: states
  integer :: ierr
  character(len=FIELD_NAME_LEN) :: external_name, periodic_name
  type(mesh_type), pointer :: periodic_mesh, external_mesh
  type(vector_field) :: periodic_positions, external_positions

  call set_debug_level(3)
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
  call populate_state(states)
  call check_valid_input(states, external_filename, external_name, periodic_name)

  external_mesh => extract_mesh(states(1), trim(external_name))
  periodic_mesh => extract_mesh(states(1), trim(periodic_name))
  external_positions = get_coordinate_field(states(1), external_mesh)
  periodic_positions = get_coordinate_field(states(1), periodic_mesh)

  call postprocess_periodic_mesh(external_mesh, external_positions, periodic_mesh, periodic_positions)

  ! Dump out the periodic mesh to disk:
  call write_triangle_files(trim(external_filename) // '_periodic', states(1), periodic_mesh)

  ! OK! Now we need to do some setting of options.
  new_filename = external_filename(1:len_trim(external_filename) - 5) // '_periodised.flml'
  call manipulate_options(external_mesh, trim(external_mesh%option_path), periodic_mesh, trim(periodic_mesh%option_path), new_filename)

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

      if (have_option(trim(mesh%option_path) // '/from_mesh/periodic_boundary_conditions')) then
        seen_periodic_meshes = seen_periodic_meshes + 1
        periodic_name = mesh%name
      end if
    end do

    if (seen_external_meshes /= 1 .or. seen_periodic_meshes /= 1) then
      call usage
      stop
    end if
  end subroutine check_valid_input

  subroutine postprocess_periodic_mesh(external_mesh, external_positions, periodic_mesh, periodic_positions)
    type(mesh_type), intent(in) :: periodic_mesh, external_mesh
    type(vector_field), intent(inout) :: periodic_positions
    type(vector_field), intent(in) :: external_positions

    integer :: periodic_bc, aliased_id, j
    integer, dimension(2) :: shape_option
    integer, dimension(:), allocatable :: aliased_boundary_ids
    integer :: face

    ! Retain the original surface IDs, as we still need them --
    ! the derivation process will have set the aliased surface elements
    ! to have surface ID zero
    periodic_mesh%faces%boundary_ids = external_mesh%faces%boundary_ids

    ! We need to loop through all aliased faces and set the periodic positions to
    ! the aliased values, so that we can recover them by applying the map later

    do periodic_bc=0,option_count(trim(periodic_mesh%option_path) // '/from_mesh/periodic_boundary_conditions')-1
      shape_option = option_shape(trim(periodic_mesh%option_path) // '/from_mesh/periodic_boundary_conditions['//int2str(periodic_bc)//']/aliased_boundary_ids') 
      allocate(aliased_boundary_ids(shape_option(1)))
      call get_option(trim(periodic_mesh%option_path) // '/from_mesh/periodic_boundary_conditions[' // int2str(periodic_bc) // ']/aliased_boundary_ids', aliased_boundary_ids)
      do j=1,shape_option(1)
        aliased_id = aliased_boundary_ids(j)
        ! Now we can finally loop over faces with this id and set appropriately
        do face=1,surface_element_count(external_mesh)
          if (surface_element_id(external_mesh, face) == aliased_id) then
            call set(periodic_positions, face_global_nodes(periodic_positions, face), face_val(external_positions, face))
          end if
        end do
      end do
      deallocate(aliased_boundary_ids)
    end do

  end subroutine postprocess_periodic_mesh

  subroutine manipulate_options(external_mesh, external_path, periodic_mesh, periodic_path, new_filename)
    type(mesh_type), intent(in) :: external_mesh, periodic_mesh
    character(len=*), intent(in) :: new_filename
    character(len=*), intent(in) :: external_path, periodic_path
    character(len=8192) :: coordinate_map
    integer :: stat, periodic_bc
    integer, dimension(2) :: shape_option
    integer, dimension(:), allocatable :: boundary_ids

    ! The periodic sub-branch of the options tree has useful information that
    ! we don't want to lose just yet. So set the external mesh to be periodic
    ! first, and then overwrite

    call delete_option(external_path // '/from_file', stat=stat)
    call set_option(external_path // '/from_mesh/mesh/name', trim(periodic_mesh%name), stat=stat)
    do periodic_bc=0,option_count(periodic_path // '/from_mesh/periodic_boundary_conditions')-1
      call add_option(external_path // '/from_mesh/periodic_boundary_conditions[' // int2str(periodic_bc) // ']/remove_periodicity', stat=stat)

      call get_option(periodic_path // '/from_mesh/periodic_boundary_conditions[' // int2str(periodic_bc) // ']/coordinate_map', coordinate_map)
      call set_option(external_path // '/from_mesh/periodic_boundary_conditions[' // int2str(periodic_bc) // ']/coordinate_map', coordinate_map, stat=stat)

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
    call set_option(periodic_path // '/from_file/file_name', new_filename, stat=stat)
    call set_option(periodic_path // '/from_file/format/name', 'triangle', stat=stat)
    call add_option(periodic_path // '/from_file/stat', stat=stat)
  end subroutine manipulate_options

end program periodise
