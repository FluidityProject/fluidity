program openexofile

  use iso_c_binding, only: C_INT, C_CHAR, C_NULL_CHAR
  use fexowrapper
  implicit none

  include "netcdf.inc"

  integer :: i
  real version
  integer(kind=c_int) :: mode = 0 ! 0=read, 1=write
  integer :: exoid, ierr
  character(kind=c_char, len=100) :: filename = "my_2bricks.exo"//C_NULL_CHAR
  character(kind=c_char, len=100) :: title
  integer :: num_dim, num_nodes, num_elem, num_elem_blk
  integer :: num_node_sets, num_side_sets
  real, allocatable, dimension(:) :: x, y, z
  integer, allocatable, dimension(:) :: node_map, elem_num_map, elem_order_map
  integer, allocatable, dimension(:) :: block_ids, num_elem_in_block, num_nodes_per_elem
  integer, allocatable, dimension(:) :: elem_blk_connectivity, elem_connectivity
  integer, allocatable, dimension(:) :: node_set_ids, num_nodes_in_set
  integer, allocatable, dimension(:) :: node_set_node_list, total_node_sets_node_list



  print*, "*************************"
  print*, "test the exodusII lib"

  if (mode == 0) then
      !exoid = f_read_ex_open2(trim(filename), mode, comp_ws, io_ws, version)
      exoid = f_read_ex_open(trim(filename), version)
  end if
  print*, "exoid: ", exoid
  print*, "exodus file version: ", version
  print*, "*************************"

  if (exoid .le. 0) then
    print*, "Failed to open Mesh File ", trim(filename), "."
    print*, "ierr = ", exoid
    stop
  end if

  ! Get database parameters from exodusII file
  ierr = f_ex_get_init(exoid, title, num_dim, num_nodes, &
                       num_elem, num_elem_blk, num_node_sets, &
                       num_side_sets)
  print*, "num_dim = ", num_dim
  print*, "num_nodes = ", num_nodes
  print*, "num_elem = ", num_elem
  print*, "num_elem_blk = ", num_elem_blk
  print*, "num_node_sets = ", num_node_sets
  print*, "num_side_sets = ", num_side_sets
  if (ierr < 0) then
    print*, "Inquire database parameters failed"
    print*, "ierr = ", ierr
    stop
  end if

  ! read nodal coordinates values and names from database
  allocate(x(num_nodes))
  allocate(y(num_nodes))
  if (num_dim >= 3) then
    allocate(z(num_nodes))
  else
    z = 0.0
  end if
  x=0.0; y=0.0; z=0.0
  ! Get coordinates from the mesh:
  ierr = f_ex_get_coord(exoid, x, y, z)
  print*, "ierr = ", ierr

  ! Read node number map:
  allocate(node_map(num_nodes))
  num_nodes = 0
  ierr = f_ex_get_node_num_map(exoid, node_map)
  print *, "ierr = ", ierr

  ! read element number map
  allocate(elem_num_map(num_elem))
  elem_num_map = 0
  ierr = f_ex_get_elem_num_map(exoid, elem_num_map)
  print*, "ierr = ", ierr

  ! read element order map
  allocate(elem_order_map(num_elem))
  elem_order_map = 0
  ierr = f_ex_get_elem_order_map(exoid, elem_order_map)
  print*, "ierr = ", ierr

  ! read element block parameters (required for element connectivity)
  allocate(block_ids(num_elem_blk))
  allocate(num_elem_in_block(num_elem_blk))
  allocate(num_nodes_per_elem(num_elem_blk))
  ierr = f_ex_get_elem_block_parameters(exoid, num_elem_blk, block_ids, num_elem_in_block, num_nodes_per_elem)
  print*, "ierr = ", ierr

  ! read element connectivity:
  allocate(elem_connectivity(0))
  do i=1, num_elem_blk
     ! Get element connectivity of block 'i' and append to global element connectivity:
     allocate(elem_blk_connectivity(num_nodes_per_elem(i) * num_elem_in_block(i)))
     ierr = f_ex_get_elem_connectivity(exoid, block_ids(i), elem_blk_connectivity)
     call append_array(elem_connectivity, elem_blk_connectivity)
     deallocate(elem_blk_connectivity)
  end do
  print*, "ierr = ", ierr

  ! Get node sets (physical boundaries):
  allocate(node_set_ids(num_node_sets))
  allocate(num_nodes_in_set(num_node_sets))
  ierr = f_ex_get_node_set_param(exoid, num_node_sets, node_set_ids, num_nodes_in_set)
!  print*, "node_set_id = ", node_set_ids
!  print*, "num_nodes_in_set = ", num_nodes_in_set
  print*, "ierr = ", ierr

  ! Get node set node lists:
  ! initial length of the array holding all the nodes with an ID
  allocate(total_node_sets_node_list(0))
  do i=1, num_node_sets
     allocate(node_set_node_list(num_nodes_in_set(i)))
     ierr = f_ex_get_node_set_node_list(exoid, num_node_sets, node_set_ids(i), node_set_node_list)
     call append_array(total_node_sets_node_list, node_set_node_list)
!     print*, "total_node_sets_node_list = ", total_node_sets_node_list
     deallocate(node_set_node_list)
  end do
  print*, "ierr = ", ierr

  ! Close ExodusII File:
  ierr = f_ex_close(exoid)
  print*, "ierr = ", ierr
  if (ierr < 0) then
    print*, "Closing ExodusII Meshfile failed"
    print*, "ierr = ", ierr
    stop
  end if

  ! Copy coordinates over to Fluidity state:



  ! Deallocate arrays of read in coordinates
  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(node_map)
  deallocate(elem_num_map)
  deallocate(elem_order_map)
  deallocate(block_ids)
  deallocate(num_elem_in_block)
  deallocate(num_nodes_per_elem)
  deallocate(elem_connectivity)
  deallocate(total_node_sets_node_list)
  deallocate(node_set_ids)
  deallocate(num_nodes_in_set)















  contains
!     subroutine resize_array(array, new_size)
!        integer, allocatable, dimension(:), intent(inout) :: array
!        integer, intent(in) :: new_size
!        integer, allocatable, dimension(:) :: tmp
!        allocate(tmp(new_size))
!        tmp(1:size(array)) = array
!        deallocate(array)
!        allocate(array(size(tmp)))
!        array = tmp
!     end subroutine resize_array

     subroutine append_array(array, array2)
        integer, allocatable, dimension(:), intent(inout) :: array
        integer, allocatable, dimension(:), intent(in) :: array2
        integer, allocatable, dimension(:) :: tmp
        allocate(tmp(size(array) + size(array2)))
        tmp(1:size(array)) = array
        tmp(size(array)+1:size(array)+size(array2)) = array2
        deallocate(array)
        allocate(array(size(tmp)))
        array = tmp
     end subroutine append_array







end program openexofile
