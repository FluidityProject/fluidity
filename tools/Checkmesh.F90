#include "fdebug.h"

subroutine checkmesh(filename_, filename_len) bind(c)
  !!< Checks the validity of the supplied mesh

! these 5 need to be on top and in this order, so as not to confuse silly old intel compiler 
  use quadrature
  use elements
  use sparse_tools
  use fields
  use state_module
!
  use halos
  use intersection_finder_module
  use linked_lists
  use meshdiagnostics
  use metric_tools
  use mesh_files
  use supermesh_construction
  use tetrahedron_intersection_module
  use iso_c_binding

  implicit none

  character(kind=c_char, len=1) :: filename_(*)
  integer(kind=c_size_t), value :: filename_len

  character(len = filename_len) :: filename
  character(len = real_format_len()) :: rformat
  integer :: global_ele, global_nodes, global_sele, global_facets, i
  type(vector_field) :: positions

  do i=1, filename_len
    filename(i:i)=filename_(i)
  end do


  rformat = real_format()

  print "(a)", "Reading in mesh mesh with base name " // trim(filename)
  positions = read_mesh_files(filename, quad_degree = 4, format="gmsh")
  if(isparallel()) call read_halos(filename, positions)
  print "(a)", "Read successful"

  call mesh_stats(positions%mesh, elements = global_ele, nodes = global_nodes, surface_elements = global_sele, facets=global_facets)
  
  call print_mesh_statistics(positions)

  call check_node_connectivity(positions)
  call check_elements(positions)
  call check_volume_element_tangling(positions)

  call deallocate(positions)

  call print_references(0)

contains

  subroutine print_mesh_statistics(positions)
    !!< Print some statistics for the supplied mesh
    
    type(vector_field), intent(in) :: positions

    type(ilist) :: seeds

    print "(a)", "Mesh statistics:"
    print "(a,i0)", "Dimension: ", positions%dim
    print "(a,i0)", "Nodes: ", global_nodes
    print "(a,i0)", "Volume elements: ", global_ele
    print "(a,i0)", "Surface elements: ", global_sele
    print "(a,i0)", "Facets: ", global_facets
    if(associated(positions%mesh%faces)) then
      if(associated(positions%mesh%faces%boundary_ids)) then
        if(any(positions%mesh%faces%boundary_ids /= 0)) then
          print "(a)", "Has boundary IDs"
        else
          print "(a)", "Has no boundary IDs"
        end if
      else
        print "(a)", "Has no boundary IDs"
      end if
    else
      print "(a)", "Has no faces information"
    end if
    if(associated(positions%mesh%region_ids)) then
      print "(a)", "Has region IDs"
    else
      print "(a)", "Has no region IDs"
    end if

    seeds = advancing_front_intersection_finder_seeds(positions)
    if(isparallel()) then
      print "(a,i0)", "Partition connectivity: ", seeds%length
    else
      print "(a,i0)", "Connectivity: ", seeds%length
    end if
    call deallocate(seeds)
    
    call print_mesh_edge_statistics(positions)
    call print_mesh_volume_statistics(positions)

  end subroutine print_mesh_statistics
  
  subroutine print_mesh_volume_statistics(positions)
    type(vector_field), intent(in) :: positions
    
    integer :: i
    real :: domain_volume, domain_surface_area, max_volume, min_volume, volume
    real, dimension(:), allocatable :: detwei
    
    domain_volume = 0.0
    min_volume = huge(0.0)
    max_volume = 0.0
    do i = 1, ele_count(positions)
      if(.not. element_owned(positions, i)) cycle
    
      volume = element_volume(positions, i)
      domain_volume = domain_volume + volume
      min_volume = min(min_volume, volume)
      max_volume = max(max_volume, volume)
    end do    
    domain_surface_area = 0.0
    do i = 1, surface_element_count(positions)
      if(.not. surface_element_owned(positions, i)) cycle
      
      allocate(detwei(face_ngi(positions, i)))
      call transform_facet_to_physical(positions, i, &
        detwei_f = detwei)      
      domain_surface_area = domain_surface_area + abs(sum(detwei))
      deallocate(detwei)
    end do
    
    call allsum(domain_volume)
    call allmin(min_volume)
    call allmax(max_volume)
    call allsum(domain_surface_area)
    
    print "(a," // rformat // ")", "Volume: ", domain_volume
    print "(a," // rformat // ")", "Surface area: ", domain_surface_area
    if(global_ele > 0) then
      print "(a," // rformat // ")", "Min element volume: ", min_volume
      print "(a," // rformat // ")", "Max element volume: ", max_volume
      print "(a," // rformat // ")", "Ratio of max to min element volumes: ", max_volume / min_volume
    end if
    
  end subroutine print_mesh_volume_statistics
  
  subroutine print_mesh_edge_statistics(positions)
    type(vector_field), intent(in) :: positions
    
    integer :: i, j
    integer, dimension(:), pointer :: nodes
    logical :: all_linear_simplices
    integer, dimension(2) :: edge_nodes
    integer, dimension(2, 2), parameter :: edge_nodes_2d = reshape((/3, 3, 2, 1/), (/2, 2/))
    integer, dimension(6, 2), parameter :: edge_nodes_3d = reshape((/4, 3, 2, 1, 1, 1, 2, 4, 3, 4, 3, 2/), (/6, 2/))
    real :: length, max_length, min_length
    real :: anisotropy, max_anisotropy, min_anisotropy
    real, dimension(positions%dim) :: evals
    real, dimension(positions%dim, positions%dim) :: edge_lengths, evecs
    type(element_type), pointer :: shape
    
    if(global_ele == 0) then
      return
    end if
    
    all_linear_simplices = .true.
    min_length = huge(0.0)
    max_length = 0.0
    min_anisotropy = huge(0.0)
    max_anisotropy = 0.0
    do i = 1, ele_count(positions)
      shape => ele_shape(positions, i)
      if(shape%degree /= 1 .or. ele_numbering_family(shape) /= FAMILY_SIMPLEX) then
        all_linear_simplices = .false.
        exit
      end if
      
      nodes => ele_nodes(positions, i)
      select case(positions%dim)
        case(3)
          do j = 1, size(edge_nodes_3d, 1)
            edge_nodes = (/nodes(edge_nodes_3d(j, 1)), nodes(edge_nodes_3d(j, 2))/)
            length = sqrt(sum((node_val(positions, edge_nodes(2)) - node_val(positions, edge_nodes(1))) ** 2))
            
            min_length = min(min_length, length)
            max_length = max(max_length, length)
          end do
        case(2)
          do j = 1, size(edge_nodes_2d, 1)
            edge_nodes = (/nodes(edge_nodes_2d(j, 1)), nodes(edge_nodes_2d(j, 2))/)
            length = sqrt(sum((node_val(positions, edge_nodes(2)) - node_val(positions, edge_nodes(1))) ** 2))
            
            min_length = min(min_length, length)
            max_length = max(max_length, length)
          end do
        case(1)
          edge_nodes = (/nodes(2), nodes(1)/)
          length = sqrt(sum((node_val(positions, edge_nodes(2)) - node_val(positions, edge_nodes(1))) ** 2))
          
          min_length = min(min_length, length)
          max_length = max(max_length, length)
        case default
          ewrite(-1, *) "For dimension: ", positions%dim
          FLAbort("Invalid dimension")
      end select
      
      edge_lengths = edge_lengths_from_metric(simplex_tensor(positions, i))
      call eigendecomposition_symmetric(edge_lengths, evecs, evals)
      anisotropy = maxval(evals) / minval(evals)
      min_anisotropy = min(min_anisotropy, anisotropy)
      max_anisotropy = max(max_anisotropy, anisotropy)
    end do
    
    call alland(all_linear_simplices)
    if(.not. all_linear_simplices) return    
    call allmin(min_length)
    call allmax(max_length)
    call allmin(min_anisotropy)
    call allmax(max_anisotropy)
    
    print "(a," // rformat // ")", "Min edge length: ", min_length
    print "(a," // rformat // ")", "Max edge length: ", max_length
    print "(a," // rformat // ")", "Ratio of max to min edge lengths: ", max_length / min_length
    print "(a," // rformat // ")", "Min anisotropy: ", min_anisotropy
    print "(a," // rformat // ")", "Max anisotropy: ", max_anisotropy
    
  end subroutine print_mesh_edge_statistics
  
  subroutine check_node_connectivity(positions)
    !!< Check the nodal connectivity of the supplied mesh
  
    type(vector_field), intent(in) :: positions
  
    integer :: i
    logical, dimension(node_count(positions)) :: connected_node
    
    print "(a)", "Checking nodal connectivity ..."
    
    connected_node = .false.
    do i = 1, ele_count(positions)
      connected_node(ele_nodes(positions, i)) = .true.
    end do
    if(all(connected_node)) then
      print "(a)", "All nodes are connected to volume elements"
    else
      do i = 1, size(connected_node)
        if(.not. connected_node(i)) then
          call print_node(positions, i)
          print "(a)", "Node not connected to any volume elements"
        end if
      end do
    end if
  
  end subroutine check_node_connectivity

  subroutine check_elements(positions)
    !!< Check that the supplied mesh for inverted or degenerate elements

    type(vector_field), intent(in) :: positions

    integer :: i
    logical :: all_ok
    real :: min_volume, volume
    real, dimension(:), allocatable :: detwei
    type(element_type), pointer :: shape

    if(global_ele > 0) then
      print "(a)", "Checking volume elements ..."
      min_volume = huge(0.0)
      all_ok = .true.
      do i = 1, ele_count(positions)               
        shape => ele_shape(positions, i)
        if(shape%degree == 1 .and. ele_numbering_family(shape) == FAMILY_SIMPLEX .and. positions%dim == 3) then
          volume = simplex_volume(positions, i)
          if(volume < 0.0) then
            print "(a)", "Inverted volume element found: "
            call print_element(i, ele_val(positions, i), ele_nodes(positions, i))
          end if
          volume = abs(volume)
        else
          volume = element_volume(positions, i)
        end if
        min_volume = min(min_volume, volume)
        if(volume < epsilon(0.0)) then
          print "(a)", "Degenerate volume element found: "
          call print_element(i, ele_val(positions, i), ele_nodes(positions, i))
          all_ok = .false.
        end if
      end do 
      
      call alland(all_ok)
      call allmin(min_volume)
      
      print "(a," // rformat // ")", "Min volume element volume: ", min_volume
      if(all_ok) then
        print "(a)", "All volume elements are non-degenerate"
      end if
    end if

    if(global_sele > 0) then
      print "(a)", "Checking surface elements ..."
      min_volume = huge(0.0)
      all_ok = .true.
      do i = 1, surface_element_count(positions)      
        allocate(detwei(face_ngi(positions, i)))
        call transform_facet_to_physical(positions, i, detwei_f = detwei)
        volume = abs(sum(detwei))
        deallocate(detwei)
        min_volume = min(min_volume, volume)
        if(volume < epsilon(0.0)) then
          print "(a)", "Degenerate surface element found: "
          call print_element(i, face_val(positions, i), face_global_nodes(positions, i))
          all_ok = .false.
        end if
      end do
      
      call alland(all_ok)
      call allmin(min_volume)
      
      print "(a," // rformat // ")", "Min surface element area: ", min_volume
      if(all_ok) then
        print "(a)", "All surface elements are non-degenerate"
      end if
    end if

  end subroutine check_elements
  
  subroutine check_volume_element_tangling(positions)
    !!< Check the supplied mesh for tangling of the volume elements
    
    type(vector_field), intent(in) :: positions
    
    integer :: dim, i, j, stat
    logical :: all_ok, intersection_found
    real :: ele_volume, intersection_volume
    real, parameter :: relative_tolerance = 1.0e-8
    type(inode), pointer :: llnode
    type(ilist), dimension(ele_count(positions)) :: intersection_map
    type(plane_type), dimension(4) :: planes_b
    type(tet_type) :: tet_a, tet_b
    type(vector_field) :: intersection
    
    if(global_ele == 0) then
      return
    end if
    
    print "(a)", "Checking volume elements for tangling ..."

    dim = positions%dim
    call intersector_set_dimension(dim)
    call intersector_set_exactness(.false.)

    all_ok = .true.
    intersection_map = intersection_finder(positions, positions)
    map_loop: do i = 1, size(intersection_map)
      ele_volume = element_volume(positions, i)

      if(dim == 3 .and. (intersector_exactness .eqv. .false.)) then
        tet_b%v = ele_val(positions, i)
        planes_b = get_planes(tet_b)
      end if
      
      llnode => intersection_map(i)%firstnode
      intersection_found = .false.
      do while(associated(llnode))
        if(dim == 3 .and. (intersector_exactness .eqv. .false.)) then
          tet_a%v = ele_val(positions, llnode%value)
          call intersect_tets(tet_a, planes_b, ele_shape(positions, llnode%value), stat = stat, output = intersection)
          if(stat /= 0) then
            llnode => llnode%next
            cycle
          end if
        else
          intersection = intersect_elements(positions, llnode%value, ele_val(positions, i), ele_shape(positions, i))
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
            call print_element(i, ele_val(positions, i), ele_nodes(positions, i))
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
      
    call alland(all_ok)
      
    if(all_ok) then
      print "(a)", "No volume element tangling"
    end if

    if(dim == 3) call finalise_tet_intersector()
    
  end subroutine check_volume_element_tangling

  subroutine print_node(positions, number)
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: number

    character(len = 1 + int2str_len(positions%dim) + real_format_len(padding = 1) + 1) :: format_buffer
    integer :: i
    real, dimension(positions%dim) :: coord

    print "(a,i0)", "Node: ", number
    
    print "(a)", "Coordinates:"
    coord = node_val(positions, number)
    format_buffer = "(" // trim(real_format()) // ")"
    do i = 1, size(coord)
      print trim(format_buffer), coord(i)
    end do
    
  end subroutine print_node

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
      print "(i0)", numbering(i)
    end do

  end subroutine print_element

end subroutine checkmesh
