#include "fdebug.h"

subroutine checkmesh(filename, filename_len)
  !!< Checks the validity of the supplied triangle mesh

  use fields
  use intersection_finder_module
  use linked_lists
  use meshdiagnostics
  use read_triangle
  use supermesh_construction
  use tetrahedron_intersection_module

  implicit none

  integer, intent(in) :: filename_len

  character(len = filename_len), intent(in) :: filename

  type(vector_field) :: mesh_field

  print *, "Reading in triangle mesh with base name " // trim(filename)
  mesh_field = read_triangle_files(trim(filename), quad_degree = 4)
  print *, "Read successful"

  call print_mesh_statistics(mesh_field)

  call check_elements(mesh_field)
  call check_volume_element_tangling(mesh_field)

  call deallocate(mesh_field)

  call print_references(0)

contains

  subroutine print_mesh_statistics(mesh_field)
    !!< Print some statistics for the supplied mesh field

    type(ilist) :: seeds
    type(vector_field), intent(in) :: mesh_field

    print "(a)", "Mesh statistics:"
    print "(a,i0)", "Dimension: ", mesh_field%dim
    print "(a,i0)", "Nodes: ", node_count(mesh_field)
    print "(a,i0)", "Volume elements: ", ele_count(mesh_field)
    print "(a,i0)", "Surface elements: ", surface_element_count(mesh_field)
    if(associated(mesh_field%mesh%faces)) then
      if(associated(mesh_field%mesh%faces%boundary_ids)) then
        print "(a)", "Has boundary IDs"
      else
        print "(a)", "Has no boundary IDs"
      end if
    else
      print "(a)", "Has no faces information"
    end if
    if(associated(mesh_field%mesh%region_ids)) then
      print "(a)", "Has region IDs"
    else
      print "(a)", "Has no region IDs"
    end if

    seeds = advancing_front_intersection_finder_seeds(mesh_field)
    print "(a,i0)", "Connectivity: ", seeds%length
    call deallocate(seeds)

  end subroutine print_mesh_statistics

  subroutine check_elements(mesh_field)
    !!< Check that the supplied mesh field for inverted or degenerate elements

    type(vector_field), intent(in) :: mesh_field

    integer :: i
    logical :: all_ok
    real :: min_volume, volume
    real, dimension(:), allocatable :: detwei
    type(element_type), pointer :: shape

    if(ele_count(mesh_field) > 0) then
      print "(a)", "Checking volume elements ..."
      min_volume = huge(0.0)
      all_ok = .true.
      do i = 1, ele_count(mesh_field)
        shape => ele_shape(mesh_field, i)
        if(shape%degree == 1 .and. ele_numbering_family(shape) == FAMILY_SIMPLEX .and. mesh_field%dim == 3) then
          volume = simplex_volume(mesh_field, i)
          if(volume < 0.0) then
            print "(a)", "Inverted volume element found: "
            call print_element(i, ele_val(mesh_field, i), ele_nodes(mesh_field, i))
          end if
          volume = abs(volume)
        else
          volume = element_volume(mesh_field, i)
        end if
        min_volume = min(min_volume, volume)
        if(volume < epsilon(0.0)) then
          print "(a)", "Degenerate volume element found: "
          call print_element(i, ele_val(mesh_field, i), ele_nodes(mesh_field, i))
          all_ok = .false.
        end if
      end do 
      print *, "Min volume element volume: ", min_volume
      if(all_ok) then
        print "(a)", "All volume elements are non-degenerate"
      end if
    end if

    if(surface_element_count(mesh_field) > 0) then
      print "(a)", "Checking surface elements ..."
      min_volume = huge(0.0)
      all_ok = .true.
      do i = 1, surface_element_count(mesh_field)
        allocate(detwei(face_ngi(mesh_field, i)))
        call transform_facet_to_physical(mesh_field, i, detwei_f = detwei)
        volume = abs(sum(detwei))
        deallocate(detwei)
        min_volume = min(min_volume, volume)
        if(volume < epsilon(0.0)) then
          print "(a)", "Degenerate surface element found: "
          call print_element(i, face_val(mesh_field, i), face_global_nodes(mesh_field, i))
          all_ok = .false.
        end if
      end do
      print *, "Min surface element area: ", min_volume
      if(all_ok) then
        print "(a)", "All surface elements are non-degenerate"
      end if
    end if

  end subroutine check_elements
  
  subroutine check_volume_element_tangling(mesh_field)
    !!< Check the supplied mesh field for tangling of the volume elements
    
    type(vector_field), intent(in) :: mesh_field
    
    integer :: dim, i, j, stat
    logical :: all_ok, intersection_found
    real :: ele_volume, intersection_volume
    real, parameter :: relative_tolerance = 1.0e-8
    type(inode), pointer :: llnode
    type(ilist), dimension(ele_count(mesh_field)) :: intersection_map
    type(plane_type), dimension(4) :: planes_b
    type(tet_type) :: tet_a, tet_b
    type(vector_field) :: intersection
    
    print "(a)", "Checking volume elements for tangling ..."

    dim = mesh_field%dim
    call intersector_set_dimension(dim)
    call intersector_set_exactness(.false.)

    all_ok = .true.
    intersection_map = intersection_finder(mesh_field, mesh_field)
    map_loop: do i = 1, size(intersection_map)
      ele_volume = element_volume(mesh_field, i)

      if(dim == 3 .and. (intersector_exactness .eqv. .false.)) then
        tet_b%v = ele_val(mesh_field, i)
        planes_b = get_planes(tet_b)
      end if
      
      llnode => intersection_map(i)%firstnode
      intersection_found = .false.
      do while(associated(llnode))
        if(dim == 3 .and. (intersector_exactness .eqv. .false.)) then
          tet_a%v = ele_val(mesh_field, llnode%value)
          call intersect_tets(tet_a, planes_b, ele_shape(mesh_field, llnode%value), stat = stat, output = intersection)
          if(stat /= 0) then
            llnode => llnode%next
            cycle
          end if
        else
          intersection = intersect_elements(mesh_field, llnode%value, ele_val(mesh_field, i), ele_shape(mesh_field, i))
          if(ele_count(intersection) == 0) then
            call deallocate(intersection)
            llnode => llnode%next
            cycle
          end if
        end if

        intersection_volume = 0.0
        do j = 1, ele_count(intersection)
          intersection_volume = intersection_volume + element_volume(intersection, j)
        end do
        if(abs(intersection_volume) > abs(max(relative_tolerance, relative_tolerance * ele_volume))) then
          ! intersection_found traps the equal volume overlaid elements
          ! case. The volume test traps the non-equal volume intersecting
          ! elements case.
          if(intersection_found .or. abs(ele_volume - intersection_volume) > abs(max(relative_tolerance, relative_tolerance * ele_volume))) then
            print "(a)", "Tangled volume element found: "
            call print_element(i, ele_val(mesh_field, i), ele_nodes(mesh_field, i))
            call deallocate(intersection)
            all_ok = .false.
            exit map_loop
          end if
          intersection_found = .true.
        end if

        call deallocate(intersection)
        llnode => llnode%next
      end do
    end do map_loop
      
    if(all_ok) then
      print "(a)", "No volume element tangling"
    end if

    if(dim == 3) call finalise_tet_intersector()
    
  end subroutine check_volume_element_tangling

  subroutine print_element(number, coords, numbering)
    !!< Print the supplied element information

    integer, intent(in) :: number
    real, dimension(:, :), intent(in) :: coords
    integer, dimension(size(coords, 2)), intent(in) :: numbering

    character(len = 1 + int2str_len(size(coords, 1)) + real_format_len(padding = 1) + 1) :: format_buffer
    integer :: i

    print "(a,i0)", "Element: ", number
    
    print "(a)", "Coordinates:"
    format_buffer = "(" // int2str(size(coords, 1)) // trim(real_format(padding = 1)) // ")"
    do i = 1, size(coords, 2)
      print trim(format_buffer), coords(:, i)
    end do

    print "(a)", "Numbering:"
    do i = 1, size(numbering)
      print *, numbering(i)
    end do

  end subroutine print_element

end subroutine checkmesh
