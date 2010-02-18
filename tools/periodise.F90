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
  use write_triangle
  use integer_set_module
  implicit none

  character(len=4096) :: filename, external_filename, new_external_filename, new_filename
  integer :: status
  type(state_type), dimension(:), pointer :: states
  integer :: ierr
  character(len=FIELD_NAME_LEN) :: external_name, periodic_name
  type(mesh_type), pointer :: periodic_mesh, external_mesh
  type(vector_field) :: periodic_positions, external_positions

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
  call populate_state(states)
  call check_valid_input(states, external_filename, external_name, periodic_name)

  external_mesh => extract_mesh(states(1), trim(external_name))
  periodic_mesh => extract_mesh(states(1), trim(periodic_name))
  external_positions = get_coordinate_field(states(1), external_mesh)
  call allocate(periodic_positions, external_positions%dim, periodic_mesh, trim(periodic_name) // 'Coordinate')
  call remap_field(external_positions, periodic_positions)

  call postprocess_periodic_mesh(external_mesh, external_positions, periodic_mesh, periodic_positions)

  ! Dump out the periodic mesh to disk:
  new_external_filename = trim(external_filename) // '_periodic'
  call write_triangle_files(new_external_filename, periodic_positions)

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
    type(mesh_type), intent(inout) :: periodic_mesh
    type(mesh_type), intent(in) :: external_mesh
    type(vector_field), intent(inout) :: periodic_positions
    type(vector_field), intent(in) :: external_positions

    type(integer_set):: physical_bc_nodes
    integer :: periodic_bc, aliased_id, physical_id, j
    integer, dimension(2) :: shape_option
    integer, dimension(:), allocatable :: aliased_boundary_ids, physical_boundary_ids
    integer :: face
    integer :: ele
    real, dimension(ele_ngi(periodic_positions, 1)) :: detwei
    real :: vol

    ! Check for degenerate elements
    do ele=1,ele_count(periodic_positions)
      call transform_to_physical(periodic_positions, ele, detwei=detwei)
      vol = sum(detwei)
      assert(vol > 0)
    end do

    ! Retain the original surface IDs, as we still need them --
    ! the derivation process will have set the aliased surface elements
    ! to have surface ID zero
    periodic_mesh%faces%boundary_ids = external_mesh%faces%boundary_ids


    ! First we need to create a set of all nodes that are on physical periodic 
    ! boundaries. If these nodes (in the non-periodic mesh) also appear on 
    ! aliased boundaries (happens for double periodic), they should not be 
    ! used to determine the position in the periodic coordinate field
    
    call allocate(physical_bc_nodes)
    do periodic_bc=0,option_count(trim(periodic_mesh%option_path) // '/from_mesh/periodic_boundary_conditions')-1
      shape_option = option_shape(trim(periodic_mesh%option_path) // '/from_mesh/periodic_boundary_conditions['//int2str(periodic_bc)//']/physical_boundary_ids')
      allocate(physical_boundary_ids(shape_option(1)))
      call get_option(trim(periodic_mesh%option_path) // '/from_mesh/periodic_boundary_conditions[' // int2str(periodic_bc) // ']/physical_boundary_ids', physical_boundary_ids)
      do j=1, shape_option(1)
        physical_id = physical_boundary_ids(j)
        ! Now we can finally loop over faces with this id and set appropriately
        do face=1, surface_element_count(external_mesh)
          if (surface_element_id(external_mesh, face) == physical_id) then
            call insert(physical_bc_nodes, face_global_nodes(external_mesh, face))
          end if
        end do
      end do
      deallocate(physical_boundary_ids)
    end do

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
          call postprocess_periodic_mesh_face(periodic_positions, external_positions, face, physical_bc_nodes)
        end do
      end do
      deallocate(aliased_boundary_ids)
      
    end do
      
    call deallocate(physical_bc_nodes)

    ! Check for degenerate elements
    do ele=1,ele_count(periodic_positions)
      call transform_to_physical(periodic_positions, ele, detwei=detwei)
      vol = sum(detwei)
      assert(vol > 0)
    end do
    
  end subroutine postprocess_periodic_mesh

  subroutine postprocess_periodic_mesh_face(periodic_positions, external_positions, face, physical_bc_nodes)
    ! Write the positions of the nodes on this aliased face
    ! to the periodic_positions. The physical positions can later be
    ! recovered by applying the coordinate map. Nodes that are on both
    ! a physical and an aliased face (happens for double/triple periodic)
    ! should not be written, as they have another copy, going back
    ! applying the inverse of the coordinate map associated with the 
    ! physical face - and its that position we want to write out.
    type(vector_field), intent(inout):: periodic_positions
    type(vector_field), intent(in):: external_positions
    integer, intent(in):: face
    type(integer_set), intent(in):: physical_bc_nodes
  
    integer, dimension(face_loc(external_positions, face)):: face_non_periodic_nodes, face_periodic_nodes
    integer:: k, np_node, p_node
    
    face_non_periodic_nodes=face_global_nodes(external_positions, face)
    face_periodic_nodes=face_global_nodes(periodic_positions, face)
    
    do k=1, size(face_non_periodic_nodes)
      np_node=face_non_periodic_nodes(k)
      p_node=face_periodic_nodes(k)
      if (.not. has_value(physical_bc_nodes, np_node)) then
        call set(periodic_positions, p_node, node_val(external_positions, np_node))
      end if
    end do
  
  end subroutine postprocess_periodic_mesh_face
  
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
      call set_option_attribute(external_path // '/from_mesh/periodic_boundary_conditions[' // int2str(periodic_bc) // ']/coordinate_map/string_value/lines', '20', stat=stat)
      call set_option_attribute(external_path // '/from_mesh/periodic_boundary_conditions[' // int2str(periodic_bc) // ']/coordinate_map/string_value/type', 'python', stat=stat)

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
    call set_option(periodic_path // '/from_file/format', 'triangle', stat=stat)
    call add_option(periodic_path // '/from_file/stat/include_in_stat', stat=stat)
    call add_option(external_path // '/from_mesh/stat/include_in_stat', stat=stat)
  end subroutine manipulate_options

end program periodise
