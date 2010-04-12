#include "fdebug.h"
#include "confdefs.h"

module zoltan_integration

#ifdef HAVE_ZOLTAN
! these 5 need to be on top and in this order, so as not to confuse silly old intel compiler 
  use quadrature
  use elements
  use sparse_tools
  use fields
  use state_module

  use field_options
  use mpi_interfaces
  use halos
  use parallel_fields
  use sparsity_patterns_meshes
  use vtk_interfaces
  use zoltan
  use linked_lists
  use global_parameters, only: real_size
  use data_structures
  use populate_state_module
  use reserve_state_module
  use boundary_conditions_from_options
  use boundary_conditions
  use metric_tools
  use c_interfaces
  use surface_id_interleaving
  use memory_diagnostics
  implicit none

  ! Module-level so that the Zoltan callback functions can access it
  type(mesh_type), save :: zz_mesh
  type(vector_field), save :: zz_positions
  type(csr_sparsity), pointer :: zz_sparsity_one, zz_sparsity_two, zz_nelist
  type(halo_type), pointer :: zz_halo, zz_ele_halo
  type(vector_field), save :: new_positions

  type(integer_hash_table), save :: nodes_we_are_sending ! in old local numbers
  type(integer_set), save :: nodes_we_are_keeping ! in old local numbers
  type(integer_hash_table), save :: universal_to_new_local_numbering, universal_to_old_local_numbering
  type(integer_set), save :: new_nodes

  type(integer_set), save :: new_elements
  type(integer_set), dimension(:), allocatable :: new_nelist

  integer, parameter :: integer_size = bit_size(0_zoltan_int)/8
  integer(zoltan_int), dimension(:), pointer :: my_import_procs => null(), my_import_global_ids => null()
  integer(zoltan_int) :: my_num_import

  type(integer_set), dimension(:), allocatable :: receives
  type(state_type), dimension(:), allocatable :: source_states, target_states
  type(integer_hash_table), save :: uen_to_new_local_numbering, uen_to_old_local_numbering

  type(integer_set), dimension(:), allocatable :: old_snelist, new_snelist
  type(integer_hash_table) :: universal_surface_number_to_surface_id
  type(integer_hash_table) :: universal_surface_number_to_element_owner
  type(integer_set), save :: new_surface_elements

  type(scalar_field), save :: element_quality, node_quality
  integer, save :: max_coplanar_id

  public :: zoltan_drive
  private

  contains

  function zoltan_cb_owned_node_count(data, ierr) result(count)
    integer(zoltan_int) :: count
    integer(zoltan_int), dimension(*) :: data ! not used
    integer(zoltan_int), intent(out) :: ierr

    count = halo_nowned_nodes(zz_halo)
    ierr = ZOLTAN_OK
  end function zoltan_cb_owned_node_count

  subroutine zoltan_cb_get_owned_nodes(data, num_gid_entries, num_lid_entries, global_ids, local_ids, wgt_dim, obj_wgts, ierr)
    integer(zoltan_int), dimension(*), intent(in) :: data ! not used
    integer(zoltan_int), intent(in) :: num_gid_entries, num_lid_entries 
    integer(zoltan_int), intent(out), dimension(*) :: global_ids 
    integer(zoltan_int), intent(out), dimension(*) :: local_ids 
    integer(zoltan_int), intent(in) :: wgt_dim 
    real(zoltan_float), intent(out), dimension(*) :: obj_wgts 
    integer(zoltan_int), intent(out) :: ierr

    integer :: count, i

    assert(num_gid_entries == 1)
    assert(num_lid_entries == 1)
    assert(wgt_dim == 1)

    count = halo_nowned_nodes(zz_halo)

    call get_owned_nodes(zz_halo, local_ids(1:count))
    global_ids(1:count) = halo_universal_number(zz_halo, local_ids(1:count))

    do i=1,count
      obj_wgts(i) = 1.0
    end do

    ierr = ZOLTAN_OK
  end subroutine zoltan_cb_get_owned_nodes

  subroutine zoltan_cb_get_num_edges(data, num_gid_entries, num_lid_entries, num_obj, global_ids, local_ids, num_edges, ierr)  
    integer(zoltan_int), dimension(*), intent(in) :: data 
    integer(zoltan_int), intent(in) :: num_gid_entries, num_lid_entries, num_obj
    integer(zoltan_int), intent(in), dimension(*) :: global_ids
    integer(zoltan_int), intent(in), dimension(*) :: local_ids
    integer(zoltan_int), intent(out),dimension(*) :: num_edges
    integer(zoltan_int), intent(out) :: ierr  

    integer :: count
    integer :: node

    assert(num_gid_entries == 1)
    assert(num_lid_entries == 1)

    count = zz_halo%nowned_nodes
    assert(count == num_obj)

    do node=1,count
      num_edges(node) = row_length(zz_sparsity_two, local_ids(node))
    end do
    ierr = ZOLTAN_OK
  end subroutine zoltan_cb_get_num_edges

subroutine zoltan_cb_get_edge_list(data, num_gid_entries, num_lid_entries, num_obj, global_ids, local_ids, &
                                  &  num_edges, nbor_global_id, nbor_procs, wgt_dim, ewgts, ierr)
    integer(zoltan_int), intent(in) :: data 
    integer(zoltan_int), intent(in) :: num_gid_entries, num_lid_entries, num_obj
    integer(zoltan_int), intent(in), dimension(*) :: global_ids
    integer(zoltan_int), intent(in), dimension(*) :: local_ids
    integer(zoltan_int), intent(in), dimension(*) :: num_edges
    integer(zoltan_int), intent(out), dimension(*) :: nbor_global_id
    integer(zoltan_int), intent(out), dimension(*) :: nbor_procs
    integer(zoltan_int), intent(in) :: wgt_dim
    real(zoltan_float), intent(out), dimension(*) :: ewgts
    integer(zoltan_int), intent(out) :: ierr 

!   elements with quality greater than this value are ok
!   those with element quality below it need to be adapted
    integer, parameter :: quality_tolerance = 0.6
    integer :: count
    integer :: node, i, j, nnode
    integer :: head
    integer, dimension(:), pointer :: neighbours

!   variables for recording various element quality functional values 
    real(zoltan_float) :: quality, min_quality

!   variables for recording the local maximum edge weight and local 90th percentile edge weight
    real(zoltan_float) :: max_weight, ninety_weight

    integer, dimension(:), pointer :: my_nelist, nbor_nelist
    type(integer_set) :: nelistA, nelistB, intersection

    integer :: ele, total_num_edges

    assert(num_gid_entries == 1)
    assert(num_lid_entries == 1)
    assert(wgt_dim == 1)

    count = zz_halo%nowned_nodes
    assert(count == num_obj)

    total_num_edges = sum(num_edges(1:num_obj))

    head = 1

!   Aim is to assign high edge weights to poor quality elements
!   so that when we load balance poor quality elements are placed
!   in the centre of partitions and can be adapted

!   loop over the nodes you own
    do node=1,count
!      find nodes neighbours
       neighbours => row_m_ptr(zz_sparsity_two, local_ids(node))

!      check the number of neighbours matches the number of edges
       assert(size(neighbours) == num_edges(node))

!      find global ids for each neighbour
       nbor_global_id(head:head+size(neighbours)-1) = halo_universal_number(zz_halo, neighbours)

!      find owning proc for each neighbour
       nbor_procs(head:head+size(neighbours)-1) = halo_node_owners(zz_halo, neighbours) - 1

!      get elements associated with current node
       my_nelist => row_m_ptr(zz_nelist, local_ids(node))

!      allocate and setup a set containing the elements associated with current node
       call allocate(nelistA)
       call insert(nelistA, my_nelist)

!      loop over all neighbouring nodes
       do j=1,size(neighbours)

!         get elements associated with neighbour node
          nbor_nelist => row_m_ptr(zz_nelist, neighbours(j))

!         allocate and setup a set containing the elements associated with the neighbouring node
          call allocate(nelistB)
          call insert(nelistB, nbor_nelist)

!         get intersection of the two nelists
          call set_intersection(intersection, nelistA, nelistB)
          call deallocate(nelistB)

          min_quality = 1.0

!         loop over all the elements in the intersection
          do i=1,key_count(intersection)
             ele = fetch(intersection, i)
!            determine the quality of the element
             quality = node_val(element_quality, ele)

!            store the element quality if it's less (worse) than any previous elements
             if (quality .LT. min_quality) then
                min_quality = quality
             end if
          end do

         call deallocate(intersection)

!        check if the quality is within the tolerance         
         if (min_quality .GT. quality_tolerance) then
!           if it is
            ewgts(head + j - 1) = 1.0
         else
!           if it's not
            ewgts(head + j - 1) = (1.0 - min_quality) * (total_num_edges + 1)
         end if
      end do

      call deallocate(nelistA)
      head = head + size(neighbours)
   end do

   assert(head == sum(num_edges(1:num_obj))+1)
   
!  calculate the local maximum edge weight
   max_weight = maxval(ewgts(1:head))

!  calculate the local 90th percentile edge weight   
   ninety_weight = max_weight * 0.9

!  make the worst 10% of elements uncuttable
   do i=1,head
      if (ewgts(i) .GT. ninety_weight) then
         ewgts(i) = total_num_edges + 1
      end if
   end do
   
   ierr = ZOLTAN_OK
 end subroutine zoltan_cb_get_edge_list
  

  ! Here is how we pack nodal positions for phase one migration:
  ! ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ! | position | sz of lv-1 nnlist | lv-1 nnlist | sz of lv-2 nnlist | lv-2 nnlist | owners of level-2 nnlist | sz of nelist | nelist | sz of snelist | snelist | snelist ids | containing element of snelist |
  ! ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  subroutine zoltan_cb_pack_node_sizes(data, num_gid_entries,  num_lid_entries, num_ids, global_ids, local_ids, sizes, ierr) 
    integer(zoltan_int), dimension(*), intent(in) :: data 
    integer(zoltan_int), intent(in) :: num_gid_entries, num_lid_entries, num_ids
    integer(zoltan_int), intent(in), dimension(*) :: global_ids,  local_ids
    integer(zoltan_int), intent(out), dimension(*) :: sizes 
    integer(zoltan_int), intent(out) :: ierr  

    integer :: i, node

    ewrite(1,*) "In zoltan_cb_pack_node_sizes"
    do i=1,num_ids
      node = local_ids(i)
      sizes(i) = zz_positions%dim * real_size + &
                 1 * integer_size + row_length(zz_sparsity_one, node) * integer_size + &
                 1 * integer_size + row_length(zz_sparsity_two, node) * integer_size * 2 + &
                 1 * integer_size + row_length(zz_nelist, node) * integer_size + &
                 1 * integer_size + key_count(old_snelist(node)) * integer_size * 3
    end do
    ierr = ZOLTAN_OK
  end subroutine zoltan_cb_pack_node_sizes

  subroutine zoltan_cb_pack_nodes(data, num_gid_entries, num_lid_entries, num_ids, global_ids, local_ids, dest, sizes, idx, buf, ierr)  
    integer(zoltan_int), dimension(*), intent(in) :: data 
    integer(zoltan_int), intent(in) :: num_gid_entries, num_lid_entries, num_ids 
    integer(zoltan_int), intent(in), dimension(*) :: global_ids 
    integer(zoltan_int), intent(in), dimension(*) :: local_ids 
    integer(zoltan_int), intent(in), dimension(*) :: dest
    integer(zoltan_int), intent(in), dimension(*) :: sizes 
    integer(zoltan_int), intent(in), dimension(*) :: idx 
    integer(zoltan_int), intent(out), dimension(*) :: buf 
    integer(zoltan_int), intent(out) :: ierr

    integer :: i, j, node, ratio, head

    ewrite(1,*) "In zoltan_cb_pack_nodes"
    ratio = real_size / integer_size
    assert(int(float(real_size) / float(integer_size)) == ratio) ! that ratio really is an int

    do i=1,num_ids
      head = idx(i)
      node = local_ids(i)
      assert(halo_universal_number(zz_halo, node) == global_ids(i))
      do j=1,zz_positions%dim
        buf(head:head+ratio-1) = transfer(node_val(zz_positions, j, node), buf(head:head+ratio-1))
        head = head + ratio
      end do
      
      buf(head) = row_length(zz_sparsity_one, node)
      head = head + 1

      buf(head:head+row_length(zz_sparsity_one, node)-1) = halo_universal_number(zz_halo, row_m_ptr(zz_sparsity_one, node))
      head = head + row_length(zz_sparsity_one, node)

      buf(head) = row_length(zz_sparsity_two, node)
      head = head + 1

      buf(head:head+row_length(zz_sparsity_two, node)-1) = halo_universal_number(zz_halo, row_m_ptr(zz_sparsity_two, node))
      head = head + row_length(zz_sparsity_two, node)

      buf(head:head+row_length(zz_sparsity_two, node)-1) = halo_node_owners(zz_halo, row_m_ptr(zz_sparsity_two, node))
      head = head + row_length(zz_sparsity_two, node)

      buf(head) = row_length(zz_nelist, node)
      head = head + 1

      buf(head:head+row_length(zz_nelist,node)-1) = halo_universal_number(zz_ele_halo, row_m_ptr(zz_nelist, node))
      head = head + row_length(zz_nelist, node)

      buf(head) = key_count(old_snelist(node))
      head = head + 1
      buf(head:head + key_count(old_snelist(node)) - 1) = set2vector(old_snelist(node))
      head = head + key_count(old_snelist(node))
      buf(head:head + key_count(old_snelist(node)) - 1) = fetch(universal_surface_number_to_surface_id, set2vector(old_snelist(node)))
      head = head + key_count(old_snelist(node))
      buf(head:head + key_count(old_snelist(node)) - 1) = fetch(universal_surface_number_to_element_owner, set2vector(old_snelist(node)))
      head = head + key_count(old_snelist(node))

      !assert(head == idx(i) + (sizes(i)/integer_size) - 1)
    end do
    ierr = ZOLTAN_OK
  end subroutine zoltan_cb_pack_nodes

  subroutine zoltan_cb_unpack_nodes(data, num_gid_entries, num_ids, global_ids, sizes, idx, buf, ierr)
    integer(zoltan_int), dimension(*), intent(inout) :: data 
    integer(zoltan_int), intent(in) :: num_gid_entries 
    integer(zoltan_int), intent(in) :: num_ids 
    integer(zoltan_int), intent(in), dimension(*) :: global_ids
    integer(zoltan_int), intent(in), dimension(*) :: sizes 
    integer(zoltan_int), intent(in), dimension(*) :: idx 
    integer(zoltan_int), intent(in), dimension(*), target :: buf 
    integer(zoltan_int), intent(out) :: ierr  

    integer :: i, ratio, head, sz, j
    type(integer_set) :: new_nodes_we_have_recorded, new_nodes_we_still_need, halo_nodes_we_currently_own
    type(mesh_type) :: new_mesh
    integer, dimension(:), pointer :: neighbours
    logical :: changed
    integer :: old_local_number, new_local_number, universal_number
    real, dimension(zz_positions%dim) :: new_coord
    type(integer_hash_table) :: universal_number_to_old_owner
    integer :: old_owner, new_owner
    integer, dimension(:), pointer :: current_buf
    integer :: rank

    ewrite(1,*) "In zoltan_cb_unpack_nodes"
    ! assert new linear mesh and positions not allocated
    assert(.not. associated(new_mesh%refcount))
    assert(.not. associated(new_positions%refcount))
    rank = getrank()

    ! Figure out the nodes we are going to know about
    call allocate(new_nodes)

    ! All the nodes we currently have and still own, in universal numbering
    do i=1,key_count(nodes_we_are_keeping)
      old_local_number = fetch(nodes_we_are_keeping, i)
      call insert(new_nodes, halo_universal_number(zz_halo, old_local_number), changed=changed)
    end do

    ! All the nodes we are receiving from other people and are going to own
    do i=1,num_ids
      call insert(new_nodes, global_ids(i))
    end do

    call allocate(universal_number_to_old_owner)
    call allocate(halo_nodes_we_currently_own)

    ! All the halos of (the nodes we currently have and still own), in universal numbering
    do i=1,key_count(nodes_we_are_keeping)
      neighbours => row_m_ptr(zz_sparsity_two, fetch(nodes_we_are_keeping, i))
      do j=1,size(neighbours)
        universal_number = halo_universal_number(zz_halo, neighbours(j))
        call insert(new_nodes, universal_number, changed=changed)
        if (changed) then ! so it is a halo node
          old_owner = halo_node_owner(zz_halo, neighbours(j)) - 1
          assert(old_owner < getnprocs())
          if (old_owner == rank) then
            call insert(halo_nodes_we_currently_own, neighbours(j))
          end if
          call insert(universal_number_to_old_owner, universal_number, old_owner)
        end if
      end do
    end do

    ! We need to process these too -- nodes we own now, but will
    ! not own later, and will become our halo nodes.
    do i=1,key_count(halo_nodes_we_currently_own)
      old_local_number = fetch(halo_nodes_we_currently_own, i)
      universal_number = halo_universal_number(zz_halo, old_local_number)
      call insert(new_nodes, universal_number)
    end do

    ! All the other nodes that will form the halo of the new nodes we are receiving 
    ratio = real_size / integer_size
    do i=1,num_ids
      current_buf => buf(idx(i):idx(i) + sizes(i)/integer_size)
      head = ratio * zz_positions%dim + 1
      ! current_buf(head) is the size of the level-1 nnlist, which we want to skip over for now
      ! we will however sum up the sizes so that we can allocate the csr sparsity later
      head = head + current_buf(head) + 1
      ! now current_buf(head) is the size of the level-2 nnlist, which we want to process
      ! and we will also sum up the sizes so that we can allocate the csr sparsity later
      sz = current_buf(head)
      do j=1,sz
        universal_number = current_buf(head + j)
        old_owner = current_buf(head + j + sz) - 1
        assert(old_owner < getnprocs())
        call insert(new_nodes, universal_number) ! the sparsity for a node includes itself
        call insert(universal_number_to_old_owner, universal_number, old_owner)
      end do
    end do

    ! Now new_nodes implicitly defines a mapping
    ! between 1 .. key_count(new_nodes) [these are the new local node numbers]
    ! and the universal node numbers of the nodes.
    ! We're going to invert that to create the hash table of universal node numbers -> local node numbers
    ! to facilitate the transfer of field information later.
    call invert_set(new_nodes, universal_to_new_local_numbering)

    ! allocate the new objects
    ! We know the number of nodes, but not the number of elements .. hmm.
    ! We will allocate it with 0 elements for now, and work it out when we
    ! invert the nnlist to compute an enlist later.
    call allocate(new_mesh, key_count(new_nodes), 0, zz_mesh%shape, trim(zz_mesh%name))
    new_mesh%option_path = zz_mesh%option_path
    call allocate(new_positions, zz_positions%dim, new_mesh, trim(zz_positions%name))
    new_positions%option_path = zz_positions%option_path
    call deallocate(new_mesh)
    allocate(new_snelist(key_count(new_nodes)))
    call allocate(new_snelist)

    ! aaaand unpack, recording which universal ids we have received
    ! so that we can figure out which ones we haven't received yet
    ! so that we can ask their old owner to send on the new details
    ! a) build a set of the nodes we have recorded
    ! b) from that, build a set of the nodes we haven't yet recorded
    ! c) figure out who owns those nodes, so we can build the import list for zoltan

    allocate(new_nelist(key_count(new_nodes)))
    do i=1,key_count(new_nodes)
      call allocate(new_nelist(i))
    end do
    call allocate(new_elements)
    call allocate(new_surface_elements)

    call allocate(new_nodes_we_have_recorded)
    ! Nodes we are keeping
    do i=1,key_count(nodes_we_are_keeping)
      old_local_number = fetch(nodes_we_are_keeping, i)
      universal_number = halo_universal_number(zz_halo, old_local_number)
      new_local_number = fetch(universal_to_new_local_numbering, universal_number)
      call insert(new_nodes_we_have_recorded, universal_number)
      call set(new_positions, new_local_number, node_val(zz_positions, old_local_number))

      ! Record the nelist information
      neighbours => row_m_ptr(zz_nelist, old_local_number)
      do j=1,size(neighbours)
        call insert(new_nelist(new_local_number), halo_universal_number(zz_ele_halo, neighbours(j)))
        call insert(new_elements, halo_universal_number(zz_ele_halo, neighbours(j)))
      end do

      ! and record the snelist information
      do j=1,key_count(old_snelist(old_local_number))
        call insert(new_snelist(new_local_number), fetch(old_snelist(old_local_number), j))
        call insert(new_surface_elements, fetch(old_snelist(old_local_number), j))
        ! we don't need to add anything to the universal_surface_number_to_surface_id because we already have it
      end do
    end do

    ! Set the positions and nelist of halo_nodes_we_currently_own
    do i=1,key_count(halo_nodes_we_currently_own)
      old_local_number = fetch(halo_nodes_we_currently_own, i)
      universal_number = halo_universal_number(zz_halo, old_local_number)
      new_local_number = fetch(universal_to_new_local_numbering, universal_number)
      call set(new_positions, new_local_number, node_val(zz_positions, old_local_number))

      neighbours => row_m_ptr(zz_nelist, old_local_number)
      do j=1,size(neighbours)
        call insert(new_nelist(new_local_number), halo_universal_number(zz_ele_halo, neighbours(j)))
        call insert(new_elements, halo_universal_number(zz_ele_halo, neighbours(j)))
      end do

      call insert(new_nodes_we_have_recorded, universal_number)

      new_owner = fetch(nodes_we_are_sending, old_local_number)
      call insert(receives(new_owner+1), universal_number)

      ! and record the snelist information
      do j=1,key_count(old_snelist(old_local_number))
        call insert(new_snelist(new_local_number), fetch(old_snelist(old_local_number), j))
        call insert(new_surface_elements, fetch(old_snelist(old_local_number), j))
      end do
    end do
    call deallocate(halo_nodes_we_currently_own)

    ! Nodes we are gaining
    do i=1,num_ids
      call insert(new_nodes_we_have_recorded, global_ids(i))
      new_local_number = fetch(universal_to_new_local_numbering, global_ids(i))
      new_coord = 0
      head = idx(i)
      do j=1,zz_positions%dim
        new_coord(j) = transfer(buf(head:head+ratio-1), new_coord(j))
        head = head + ratio
      end do
      call set(new_positions, new_local_number, new_coord)

      ! Record the nelist information
      sz = buf(head) ! level-1 nnlist
      head = head + sz + 1
      sz = buf(head) ! level-2 nnlist
      head = head + 2*sz + 1
      sz = buf(head) ! nelist
      do j=1,sz
        call insert(new_nelist(new_local_number), buf(head + j))
        call insert(new_elements, buf(head + j))
      end do
      head = head + sz + 1

      ! And record the snelist information
      sz = buf(head)
      do j=1,sz
        call insert(new_snelist(new_local_number), buf(head + j))
        call insert(universal_surface_number_to_surface_id, buf(head + j), buf(head + j + sz))
        call insert(universal_surface_number_to_element_owner, buf(head + j), buf(head + j + 2*sz))
        call insert(new_surface_elements, buf(head + j))
      end do
      head = head + 3*sz + 1
    end do

    ! At this point, there might still be nodes that we have not yet recorded but
    ! we own, so we can fill them in now.
    call allocate(new_nodes_we_still_need)
    do i=1,key_count(new_nodes)
      universal_number = fetch(new_nodes, i)
      if (has_value(new_nodes_we_have_recorded, universal_number)) cycle

      old_owner = fetch(universal_number_to_old_owner, universal_number)
      if (old_owner == rank) then
        call insert(new_nodes_we_have_recorded, universal_number)
        old_local_number = fetch(universal_to_old_local_numbering, universal_number)
        new_local_number = fetch(universal_to_new_local_numbering, universal_number)
        call set(new_positions, new_local_number, node_val(zz_positions, old_local_number))

        ! Record the nelist information
        neighbours => row_m_ptr(zz_nelist, old_local_number)
        do j=1,size(neighbours)
          call insert(new_nelist(new_local_number), halo_universal_number(zz_ele_halo, neighbours(j)))
          call insert(new_elements, halo_universal_number(zz_ele_halo, neighbours(j)))
        end do

        ! and record the snelist information
        do j=1,key_count(old_snelist(old_local_number))
          call insert(new_snelist(new_local_number), fetch(old_snelist(old_local_number), j))
          call insert(new_surface_elements, fetch(old_snelist(old_local_number), j))
          ! we don't need to add anything to the universal_surface_number_to_surface_id because we already have it
        end do

        ! and record the node in the receives
        new_owner = fetch(nodes_we_are_sending, old_local_number)
        call insert(receives(new_owner+1), universal_number)
      else
        call insert(new_nodes_we_still_need, universal_number)
      end if
    end do

    ! And build the import list ...
    my_num_import = key_count(new_nodes_we_still_need)
    allocate(my_import_procs(my_num_import))
    allocate(my_import_global_ids(my_num_import))
    do i=1,my_num_import
      universal_number = fetch(new_nodes_we_still_need, i)
      my_import_global_ids(i) = universal_number
      my_import_procs(i) = fetch(universal_number_to_old_owner, universal_number)
      assert(my_import_procs(i) /= rank)
    end do

    call deallocate(new_nodes_we_have_recorded)
    call deallocate(new_nodes_we_still_need)
    call deallocate(universal_number_to_old_owner)

    ierr = ZOLTAN_OK
  end subroutine zoltan_cb_unpack_nodes

  ! Here is how we pack halo nodes for phase two migration:
  ! -----------------------------------------------------------------------------------------------------------------------------------------------------
  ! | position | new owner | size of nelist | nelist | size of snelist | snelist | surface ids | the containing volume element for each surface element |
  ! -----------------------------------------------------------------------------------------------------------------------------------------------------
  subroutine zoltan_cb_pack_halo_node_sizes(data, num_gid_entries,  num_lid_entries, num_ids, global_ids, local_ids, sizes, ierr) 
    integer(zoltan_int), dimension(*), intent(in) :: data 
    integer(zoltan_int), intent(in) :: num_gid_entries, num_lid_entries, num_ids
    integer(zoltan_int), intent(in), dimension(*) :: global_ids,  local_ids
    integer(zoltan_int), intent(out), dimension(*) :: sizes 
    integer(zoltan_int), intent(out) :: ierr  

    integer :: i, node
    ewrite(1,*) "In zoltan_cb_pack_halo_node_sizes"

    do i=1,num_ids
      node = fetch(universal_to_old_local_numbering, global_ids(i))
      sizes(i) = zz_positions%dim * real_size + &
                 2 * integer_size + row_length(zz_nelist, node) * integer_size + &
                 1 * integer_size + key_count(old_snelist(node)) * 3 * integer_size
    end do
    ierr = ZOLTAN_OK
  end subroutine zoltan_cb_pack_halo_node_sizes

  subroutine zoltan_cb_pack_halo_nodes(data, num_gid_entries, num_lid_entries, num_ids, global_ids, local_ids, dest, sizes, idx, buf, ierr)  
    integer(zoltan_int), dimension(*), intent(in) :: data 
    integer(zoltan_int), intent(in) :: num_gid_entries, num_lid_entries, num_ids 
    integer(zoltan_int), intent(in), dimension(*) :: global_ids 
    integer(zoltan_int), intent(in), dimension(*) :: local_ids 
    integer(zoltan_int), intent(in), dimension(*) :: dest
    integer(zoltan_int), intent(in), dimension(*) :: sizes 
    integer(zoltan_int), intent(in), dimension(*) :: idx 
    integer(zoltan_int), intent(out), dimension(*), target :: buf 
    integer(zoltan_int), intent(out) :: ierr

    integer :: i, j, node, ratio, head, new_owner, rank
    integer, dimension(:), pointer :: current_buf

    ewrite(1,*) "In zoltan_cb_pack_halo_nodes"
    ratio = real_size / integer_size
    rank = getrank()

    do i=1,num_ids
      current_buf => buf(idx(i):idx(i)+sizes(i)/integer_size)
      head = 1
      node = fetch(universal_to_old_local_numbering, global_ids(i))
      do j=1,zz_positions%dim
        current_buf(head:head+ratio-1) = transfer(node_val(zz_positions, j, node), current_buf(head:head+ratio-1))
        head = head + ratio
      end do
      
      ! Now compute the new owner
      if (has_key(nodes_we_are_sending, node)) then
        new_owner = fetch(nodes_we_are_sending, node)
      else
        new_owner = rank
      end if

      current_buf(head) = new_owner
      head = head + 1

      current_buf(head) = row_length(zz_nelist, node)
      head = head + 1

      current_buf(head:head+row_length(zz_nelist, node)-1) = halo_universal_number(zz_ele_halo, row_m_ptr(zz_nelist, node))
      head = head + row_length(zz_nelist, node)

      current_buf(head) = key_count(old_snelist(node))
      head = head + 1
      current_buf(head:head+key_count(old_snelist(node))-1) = set2vector(old_snelist(node))
      head = head + key_count(old_snelist(node))
      current_buf(head:head+key_count(old_snelist(node))-1) = fetch(universal_surface_number_to_surface_id, set2vector(old_snelist(node)))
      head = head + key_count(old_snelist(node))
      current_buf(head:head+key_count(old_snelist(node))-1) = fetch(universal_surface_number_to_element_owner, set2vector(old_snelist(node)))
      head = head + key_count(old_snelist(node))

      !assert(head == (sizes(i)/integer_size)+1)
    end do
    ierr = ZOLTAN_OK
  end subroutine zoltan_cb_pack_halo_nodes

  subroutine zoltan_cb_unpack_halo_nodes(data, num_gid_entries, num_ids, global_ids, sizes, idx, buf, ierr)
    integer(zoltan_int), dimension(*), intent(inout) :: data 
    integer(zoltan_int), intent(in) :: num_gid_entries 
    integer(zoltan_int), intent(in) :: num_ids 
    integer(zoltan_int), intent(in), dimension(*) :: global_ids
    integer(zoltan_int), intent(in), dimension(*) :: sizes 
    integer(zoltan_int), intent(in), dimension(*) :: idx 
    integer(zoltan_int), intent(in), dimension(*), target :: buf 
    integer(zoltan_int), intent(out) :: ierr  

    integer :: i, j
    real, dimension(zz_positions%dim) :: new_coord
    integer :: head
    integer :: ratio
    integer :: new_local_number, new_owner, sz

    ewrite(1,*) "In zoltan_cb_unpack_halo_nodes"

    ratio = real_size/integer_size

    do i=1,num_ids
      new_local_number = fetch(universal_to_new_local_numbering, global_ids(i))
      new_coord = 0
      head = idx(i)
      do j=1,zz_positions%dim
        new_coord(j) = transfer(buf(head:head+ratio-1), new_coord(j))
        head = head + ratio
      end do
      call set(new_positions, new_local_number, new_coord)

      new_owner = buf(head)
      head = head + 1

      ! record the nelist information
      sz = buf(head)
      do j=1,sz
        call insert(new_nelist(new_local_number), buf(head + j))
        call insert(new_elements, buf(head + j))
      end do
      head = head + sz + 1

      ! and record who owns this in the halo
      call insert(receives(new_owner+1), global_ids(i))

      ! and record the snelist information
      sz = buf(head)
      do j=1,sz
        call insert(new_snelist(new_local_number), buf(head + j))
        call insert(universal_surface_number_to_surface_id, buf(head + j), buf(head + j + sz))
        call insert(universal_surface_number_to_element_owner, buf(head + j), buf(head + j + 2*sz))
        call insert(new_surface_elements, buf(head + j))
      end do
      head = head + 3*sz + 1
    end do

    ierr = ZOLTAN_OK
  end subroutine zoltan_cb_unpack_halo_nodes

  subroutine zoltan_cb_pack_field_sizes(data, num_gid_entries,  num_lid_entries, num_ids, global_ids, local_ids, sizes, ierr) 
    integer(zoltan_int), dimension(*), intent(in) :: data 
    integer(zoltan_int), intent(in) :: num_gid_entries, num_lid_entries, num_ids
    integer(zoltan_int), intent(in), dimension(*) :: global_ids,  local_ids
    integer(zoltan_int), intent(out), dimension(*) :: sizes 
    integer(zoltan_int), intent(out) :: ierr  

    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield
    integer :: state_no, field_no, sz, i

    ewrite(1,*) "In zoltan_cb_pack_field_sizes"

    ! The person doing this for mixed meshes in a few years time: this is one of the things
    ! you need to change. Make it look at the loc for each element.

    sz = 0

    do state_no=1,size(source_states)
      do field_no=1,scalar_field_count(source_states(state_no))
        sfield => extract_scalar_field(source_states(state_no), field_no)
        sz = sz + ele_loc(sfield, 1)
      end do

      do field_no=1,vector_field_count(source_states(state_no))
        vfield => extract_vector_field(source_states(state_no), field_no)
        sz = sz + ele_loc(vfield, 1) * vfield%dim
      end do

      do field_no=1,tensor_field_count(source_states(state_no))
        tfield => extract_tensor_field(source_states(state_no), field_no)
        sz = sz + ele_loc(tfield, 1) * tfield%dim**2
      end do
    end do

    do i=1,num_ids
      sizes(i) = sz * real_size
    end do

    ierr = ZOLTAN_OK
  end subroutine zoltan_cb_pack_field_sizes

  subroutine zoltan_cb_pack_fields(data, num_gid_entries, num_lid_entries, num_ids, global_ids, local_ids, dest, sizes, idx, buf, ierr)  
    integer(zoltan_int), dimension(*), intent(in) :: data 
    integer(zoltan_int), intent(in) :: num_gid_entries, num_lid_entries, num_ids 
    integer(zoltan_int), intent(in), dimension(*) :: global_ids 
    integer(zoltan_int), intent(in), dimension(*) :: local_ids 
    integer(zoltan_int), intent(in), dimension(*) :: dest
    integer(zoltan_int), intent(in), dimension(*) :: sizes 
    integer(zoltan_int), intent(in), dimension(*) :: idx 
    integer(zoltan_int), intent(out), dimension(*), target :: buf 
    integer(zoltan_int), intent(out) :: ierr

    real, dimension(sizes(1) / real_size) :: rbuf ! easier to write reals to real memory
    integer :: rhead, i, state_no, field_no, loc
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield
    integer :: old_universal_element_number, old_local_element_number

    ewrite(1,*) "In zoltan_cb_pack_fields"

    do i=1,num_ids
      old_universal_element_number = global_ids(i)
      old_local_element_number = fetch(uen_to_old_local_numbering, old_universal_element_number)

      rhead = 1
      do state_no=1,size(source_states)
        do field_no=1,scalar_field_count(source_states(state_no))
          sfield => extract_scalar_field(source_states(state_no), field_no)
          loc = ele_loc(sfield, old_local_element_number)
          rbuf(rhead:rhead + loc - 1) = ele_val(sfield, old_local_element_number)
          rhead = rhead + loc
        end do

        do field_no=1,vector_field_count(source_states(state_no))
          vfield => extract_vector_field(source_states(state_no), field_no)
          if (index(vfield%name,"Coordinate")==len_trim(vfield%name)-9) cycle
          loc = ele_loc(vfield, old_local_element_number)
          rbuf(rhead:rhead + loc*vfield%dim - 1) = reshape(ele_val(vfield, old_local_element_number), (/loc*vfield%dim/))
          rhead = rhead + loc * vfield%dim
        end do

        do field_no=1,tensor_field_count(source_states(state_no))
          tfield => extract_tensor_field(source_states(state_no), field_no)
          loc = ele_loc(tfield, old_local_element_number)
          rbuf(rhead:rhead + loc*(tfield%dim**2) - 1) = reshape(ele_val(tfield, old_local_element_number), (/loc*(tfield%dim**2)/))
          rhead = rhead + loc * tfield%dim**2
        end do
      end do

      buf(idx(i):idx(i) + (sizes(i)/real_size) - 1) = transfer(rbuf, buf(idx(i):idx(i) + (sizes(i)/real_size) - 1))
    end do

    ierr = ZOLTAN_OK
  end subroutine zoltan_cb_pack_fields

  subroutine zoltan_cb_unpack_fields(data, num_gid_entries, num_ids, global_ids, sizes, idx, buf, ierr)
    integer(zoltan_int), dimension(*), intent(inout) :: data 
    integer(zoltan_int), intent(in) :: num_gid_entries 
    integer(zoltan_int), intent(in) :: num_ids 
    integer(zoltan_int), intent(in), dimension(*) :: global_ids
    integer(zoltan_int), intent(in), dimension(*) :: sizes 
    integer(zoltan_int), intent(in), dimension(*) :: idx 
    integer(zoltan_int), intent(in), dimension(*), target :: buf 
    integer(zoltan_int), intent(out) :: ierr  

    real, dimension(sizes(1) / real_size) :: rbuf ! easier to read reals 
    integer :: rhead, i, state_no, field_no, loc, ihead, ratio
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield
    integer :: old_universal_element_number, new_local_element_number

    ewrite(1,*) "In zoltan_cb_unpack_fields"

    ratio = real_size / integer_size

    do i=1,num_ids
      old_universal_element_number = global_ids(i)
      new_local_element_number = fetch(uen_to_new_local_numbering, old_universal_element_number)

      ! This is bugged in gfortran 4.4:
      !rbuf = transfer(buf(idx(i):idx(i) + (sizes(i)/real_size) - 1), rbuf, size(rbuf))
      call memcpy(rbuf, buf(idx(i):idx(i) + (sizes(i)/real_size) - 1), size(rbuf) * real_size)

      rhead = 1
      do state_no=1,size(target_states)
        do field_no=1,scalar_field_count(target_states(state_no))
          sfield => extract_scalar_field(target_states(state_no), field_no)
          loc = ele_loc(sfield, new_local_element_number)
          call set(sfield, ele_nodes(sfield, new_local_element_number), rbuf(rhead:rhead + loc - 1))
          rhead = rhead + loc
        end do

        do field_no=1,vector_field_count(target_states(state_no))
          vfield => extract_vector_field(target_states(state_no), field_no)
          if (index(vfield%name,"Coordinate")==len_trim(vfield%name)-9) cycle
          loc = ele_loc(vfield, new_local_element_number)
          call set(vfield, ele_nodes(vfield, new_local_element_number), reshape(rbuf(rhead:rhead + loc*vfield%dim - 1), (/vfield%dim, loc/)))
          rhead = rhead + loc * vfield%dim
        end do

        do field_no=1,tensor_field_count(target_states(state_no))
          tfield => extract_tensor_field(target_states(state_no), field_no)
          loc = ele_loc(tfield, new_local_element_number)
          call set(tfield, ele_nodes(tfield, new_local_element_number), reshape(rbuf(rhead:rhead + loc*(tfield%dim**2) - 1), (/tfield%dim, tfield%dim, loc/)))
          rhead = rhead + loc * tfield%dim**2
        end do
      end do

    end do
    ierr = ZOLTAN_OK
  end subroutine zoltan_cb_unpack_fields

  subroutine zoltan_drive(states, iteration, metric)
    type(state_type), dimension(:), intent(inout), target :: states
    integer, intent(in) :: iteration
    type(tensor_field), intent(inout), optional :: metric
    type(zoltan_struct), pointer :: zz

    logical :: changes
    integer(zoltan_int) :: num_gid_entries, num_lid_entries
    integer(zoltan_int), dimension(:), pointer :: p1_export_global_ids
    integer(zoltan_int), dimension(:), pointer :: p1_export_local_ids
    integer(zoltan_int), dimension(:), pointer :: p1_export_procs
    integer(zoltan_int) :: p1_num_import, p1_num_export
    integer(zoltan_int), dimension(:), pointer :: p1_import_global_ids
    integer(zoltan_int), dimension(:), pointer :: p1_import_local_ids
    integer(zoltan_int), dimension(:), pointer :: p1_import_procs
    integer, save :: dumpno = 0

    type(tensor_field) :: new_metric

    ewrite(1,*) "In zoltan_drive"

!    call vtk_write_state("input_states", index=dumpno, state=states)

    call setup_module_variables

    call set_zoltan_parameters

    call zoltan_load_balance

    if (changes .eqv. .false.) then
      ewrite(1,*) "Zoltan decided no change was necessary, exiting"
      call deallocate_zoltan_lists
      call cleanup_basic_module_variables
      dumpno = dumpno + 1
      return
    end if

    call dump_suggested_owner

    ! The general plan:
    ! Just send the nodes you own, along with a note of their dependencies
    ! The receiving process loops through all its receives and records whom
    ! it needs to receive from (the OLD owner)

    ! It builds an import list from that, then migrates again

    call are_we_keeping_or_sending_nodes

    ! Migrate here
    call zoltan_migration_phase_one ! for nodes I am going to own
    call deallocate_zoltan_lists
    call deal_with_exports
    call zoltan_migration_phase_two ! for halo nodes those nodes depend on
    call deallocate_my_lists

    call reconstruct_enlist
    call reconstruct_senlist
    call reconstruct_halo
    call dump_linear_mesh

    ! At this point, we now have the balanced linear external mesh.
    ! Get populate_state to allocate the fields and such on this new
    ! mesh.

    call initialise_transfer

    ! And now transfer the field data around.
    call transfer_fields

    call finalise_transfer

    call cleanup_basic_module_variables
    call cleanup_other_module_variables
    dumpno = dumpno + 1

    ewrite(1,*) "Exiting zoltan_drive"

    contains

    subroutine transfer_fields
      ! OK! So, here is how this is going to work. We are going to 
      ! loop through every element in which you own at least one node, and note
      ! that its information needs to be sent to the owners of its vertices.
      ! We also have to take special care of self-sends, since they don't
      ! get taken care of in the zoltan communication.

      integer :: old_ele
      integer, dimension(:), pointer :: old_local_nodes
      type(integer_set), dimension(halo_proc_count(zz_halo)) :: sends
      integer :: i, new_owner, universal_element_number
      type(integer_set) :: self_sends
      integer :: num_import, num_export
      integer, dimension(:), pointer :: import_global_ids, import_local_ids, import_procs, import_to_part
      integer, dimension(:), pointer :: export_global_ids, export_local_ids, export_procs, export_to_part
      integer :: head
      integer(zoltan_int) :: ierr

      integer :: old_universal_element_number, new_local_element_number, old_local_element_number
      integer :: state_no, field_no
      type(scalar_field), pointer :: source_sfield, target_sfield
      type(vector_field), pointer :: source_vfield, target_vfield
      type(tensor_field), pointer :: source_tfield, target_tfield

      do i=1,size(sends)
        call allocate(sends(i))
      end do
      call allocate(self_sends)

      do old_ele=1,ele_count(zz_positions)
        universal_element_number = halo_universal_number(zz_ele_halo, old_ele)
        old_local_nodes => ele_nodes(zz_positions, old_ele)
        if (.not. any(nodes_owned(zz_halo, old_local_nodes))) cycle
        do i=1,size(old_local_nodes)
          if (has_value(nodes_we_are_keeping, old_local_nodes(i))) then
            assert(node_owned(zz_halo, old_local_nodes(i)))
            call insert(self_sends, universal_element_number)
          else if (has_key(nodes_we_are_sending,  old_local_nodes(i))) then
            new_owner = fetch(nodes_we_are_sending, old_local_nodes(i))
            call insert(sends(new_owner+1), universal_element_number)
          end if
        end do
      end do

      num_export = sum(key_count(sends))
      allocate(export_global_ids(num_export))
      allocate(export_procs(num_export))

      head = 1
      do i=1,size(sends)
        export_global_ids(head:head + key_count(sends(i)) - 1) = set2vector(sends(i))
        export_procs(head:head + key_count(sends(i)) - 1) = i - 1
        head = head + key_count(sends(i))
      end do

      allocate(export_local_ids(num_export))
      export_local_ids = 666

      import_global_ids => null()
      import_local_ids => null()
      import_procs => null()
      import_to_part => null()
      export_to_part => null()

      ierr = Zoltan_Compute_Destinations(zz, &
       & num_export, export_global_ids, export_local_ids, export_procs, &
       & num_import, import_global_ids, import_local_ids, import_procs)
      assert(ierr == ZOLTAN_OK)

      ierr = Zoltan_Migrate(zz, num_import, import_global_ids, import_local_ids, import_procs, &
       & import_to_part, num_export, export_global_ids, export_local_ids, export_procs, export_to_part) 
      assert(ierr == ZOLTAN_OK)

      deallocate(export_local_ids)
      deallocate(export_procs)
      deallocate(export_global_ids)

      ierr = Zoltan_LB_Free_Part(import_global_ids, import_local_ids, import_procs, import_to_part)
      assert(ierr == ZOLTAN_OK)

      do state_no=1,size(source_states)
        assert(scalar_field_count(source_states(state_no)) == scalar_field_count(target_states(state_no)))
        assert(vector_field_count(source_states(state_no)) == vector_field_count(target_states(state_no)))
        assert(tensor_field_count(source_states(state_no)) == tensor_field_count(target_states(state_no)))

        do field_no=1,scalar_field_count(source_states(state_no))
          source_sfield => extract_scalar_field(source_states(state_no), field_no)
          target_sfield => extract_scalar_field(target_states(state_no), field_no)
          assert(trim(source_sfield%name) == trim(target_sfield%name))

          do i=1,key_count(self_sends)
            old_universal_element_number = fetch(self_sends, i)
            new_local_element_number = fetch(uen_to_new_local_numbering, old_universal_element_number)
            old_local_element_number = fetch(uen_to_old_local_numbering, old_universal_element_number)
            call set(target_sfield, ele_nodes(target_sfield, new_local_element_number), ele_val(source_sfield, old_local_element_number))
          end do
        end do

        do field_no=1,vector_field_count(source_states(state_no))
          source_vfield => extract_vector_field(source_states(state_no), field_no)
          target_vfield => extract_vector_field(target_states(state_no), field_no)
          assert(trim(source_vfield%name) == trim(target_vfield%name))
          if (source_vfield%name == new_positions%name) cycle

          do i=1,key_count(self_sends)
            old_universal_element_number = fetch(self_sends, i)
            new_local_element_number = fetch(uen_to_new_local_numbering, old_universal_element_number)
            old_local_element_number = fetch(uen_to_old_local_numbering, old_universal_element_number)
            call set(target_vfield, ele_nodes(target_vfield, new_local_element_number), ele_val(source_vfield, old_local_element_number))
          end do
        end do

        do field_no=1,tensor_field_count(source_states(state_no))
          source_tfield => extract_tensor_field(source_states(state_no), field_no)
          target_tfield => extract_tensor_field(target_states(state_no), field_no)
          assert(trim(source_tfield%name) == trim(target_tfield%name))

          do i=1,key_count(self_sends)
            old_universal_element_number = fetch(self_sends, i)
            new_local_element_number = fetch(uen_to_new_local_numbering, old_universal_element_number)
            old_local_element_number = fetch(uen_to_old_local_numbering, old_universal_element_number)
            call set(target_tfield, ele_nodes(target_tfield, new_local_element_number), ele_val(source_tfield, old_local_element_number))
          end do
        end do
      end do

      call deallocate(self_sends)
      call deallocate(sends)

      call halo_update(target_states)
      
    end subroutine transfer_fields

    subroutine finalise_transfer
      integer :: i
      call set_prescribed_field_values(states, exclude_interpolated = .true.)
      call populate_boundary_conditions(states)
      call set_boundary_conditions_values(states)
      call set_dirichlet_consistent(states)
      call alias_fields(states)

      if (present(metric)) then
        metric = new_metric
        call halo_update(metric)
      end if

      do i=1,size(source_states)
        call deallocate(source_states(i))
        call deallocate(target_states(i))
      end do
      deallocate(source_states)
      deallocate(target_states)
    end subroutine finalise_transfer

    subroutine initialise_transfer
      integer :: i
      type(state_type), dimension(size(states)) :: interpolate_states
      integer(zoltan_int) :: ierr

      ! Set up source_states
      do i=1,size(states)
        call select_fields_to_interpolate(states(i), interpolate_states(i), no_positions=.true.)
        call deallocate(states(i))
      end do

      ! Interpolate the metric, too
      if (present(metric)) then
        call insert(interpolate_states(1), metric, "ErrorMetric")
        call deallocate(metric)
      end if

      allocate(source_states(mesh_count(interpolate_states(1))))
      call halo_update(interpolate_states, level=1)
      call collect_fields_by_mesh(interpolate_states, source_states)

      do i=1,size(interpolate_states)
        call deallocate(interpolate_states(i))
      end do

      if (mesh_periodic(zz_mesh)) then
        new_positions%mesh%periodic = .true.
      end if

      call insert(states, new_positions%mesh, name = new_positions%mesh%name)
      call insert(states, new_positions, name = new_positions%name)

      if (present(metric)) then
        call allocate(new_metric, new_positions%mesh, "ErrorMetric")
        call zero(new_metric)
      end if

      call deallocate(new_positions)

      call restore_reserved_meshes(states)
      call insert_derived_meshes(states)
      call allocate_and_insert_fields(states)
      call restore_reserved_fields(states)

      ! And set up target_states
      do i=1,size(states)
        call select_fields_to_interpolate(states(i), interpolate_states(i), no_positions=.true.)
      end do

      if (present(metric)) then
        call insert(interpolate_states(1), new_metric, "ErrorMetric")
      end if

      allocate(target_states(mesh_count(interpolate_states(1))))
      call collect_fields_by_mesh(interpolate_states, target_states)

      do i=1,size(interpolate_states)
        call deallocate(interpolate_states(i))
      end do

      ierr = Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE, zoltan_cb_pack_field_sizes); assert(ierr == ZOLTAN_OK)
      ierr = Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_MULTI_FN_TYPE, zoltan_cb_pack_fields); assert(ierr == ZOLTAN_OK)
      ierr = Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE, zoltan_cb_unpack_fields); assert(ierr == ZOLTAN_OK)

    end subroutine initialise_transfer

    subroutine reconstruct_halo
      ! At this point, the receives sets have been populated with all
      ! the universal node numbers we need to receive from each process.
      ! So, we are going to use zoltan to invert this to compute
      ! the send list for each process too.
      ! Then we will allocate the l2n halo and set it.
      ! Then we will chop it down to form the l1n halo, the l1e halo, and the
      ! l2e halo. 
      ! Supply the peeps with jeeps, brick apiece, capiche?

      integer :: num_import, num_export
      integer, dimension(:), pointer :: import_global_ids, import_local_ids, import_procs
      integer, dimension(:), pointer :: export_global_ids, export_local_ids, export_procs, export_to_part
      integer :: ierr, i, head
      type(integer_set), dimension(size(receives)) :: sends
      integer, dimension(size(receives)) :: nreceives, nsends
      integer, dimension(ele_count(new_positions)) :: renumber_permutation
      integer :: universal_element_number, old_new_local_element_number, new_new_local_element_number
      integer :: stat
      
      ewrite(1,*) "In reconstruct_halo"

      num_import = 0
      do i=1,size(receives)
        nreceives(i) = key_count(receives(i))
        num_import = num_import + nreceives(i)
      end do

      allocate(import_global_ids(num_import))
      allocate(import_local_ids(num_import))
      allocate(import_procs(num_import))

      import_local_ids = 666
      head = 1
      do i=1,size(receives)
        import_global_ids(head:head + nreceives(i) - 1) = set2vector(receives(i))
        import_procs(head:head + nreceives(i) - 1) = i - 1
        head = head + nreceives(i)
      end do

      export_global_ids => null()
      export_local_ids => null()
      export_procs => null()
      export_to_part => null()

      ierr = Zoltan_Compute_Destinations(zz, &
       & num_import, import_global_ids, import_local_ids, import_procs, &
       & num_export, export_global_ids, export_local_ids, export_procs)
      assert(ierr == ZOLTAN_OK)

      ! Now we know the sends too! Thanks, Zoltan!

      deallocate(import_global_ids)
      deallocate(import_local_ids)
      deallocate(import_procs)

      ! Now create the sends sets .. easier than pulling it straight out of Zoltan's data structures
      ! as Zoltan does NOT explicitly guarantee that the sends are organised such that export_procs looks
      ! like 
      ! [ sends to process 0 | sends to process 1 | sends to process 2 ... ]
      ! If such a guarantee were available, then it would be just as easy to use that, but there
      ! is no such guarantee given in the documentation ...
      ! Is an appetite for destruction, slap a murder rap on this production

      do i=1,size(sends)
        call allocate(sends(i))
      end do

      do i=1,num_export
        call insert(sends(export_procs(i)+1), export_global_ids(i))
      end do

      do i=1,size(sends)
        nsends(i) = key_count(sends(i))
      end do

      ! Allocate the halo and such
      ! We had to grow dreads to change our description, two cops is on a milkbox, missin'

      allocate(new_positions%mesh%halos(2))
      call allocate(new_positions%mesh%halos(2), &
                    nsends = nsends, &
                    nreceives = nreceives, &
                    name = halo_name(zz_halo), &
                    communicator = halo_communicator(zz_halo), &
                    nowned_nodes = key_count(new_nodes) - num_import, &
                    data_type = halo_data_type(zz_halo))

      do i=1,size(receives)
        call set_halo_sends(new_positions%mesh%halos(2), i, fetch(universal_to_new_local_numbering, set2vector(sends(i))))
        call set_halo_receives(new_positions%mesh%halos(2), i, fetch(universal_to_new_local_numbering, set2vector(receives(i))))
      end do

      ! Now derive all the other halos ...
      ! And me teevee's always off, cause I see something truly black then

      call derive_l1_from_l2_halo(new_positions%mesh, ordering_scheme = HALO_ORDER_GENERAL, create_caches = .false.)
      call renumber_positions_trailing_receives(new_positions)
      assert(has_ownership(new_positions%mesh%halos(2)))
      assert(has_ownership(new_positions%mesh%halos(1)))
      allocate(new_positions%mesh%element_halos(2))
      call derive_element_halo_from_node_halo(new_positions%mesh, create_caches = .false.)
      call renumber_positions_elements_trailing_receives(new_positions, permutation=renumber_permutation)
      assert(has_ownership(new_positions%mesh%halos(2)))
      assert(has_ownership(new_positions%mesh%halos(1)))

      ! The previous routine has renumbered all the local elements to put the element halos in
      ! trailing receives status. However, we need the universal element number -> local element number
      ! later in the field transfer. So we need to update our records now.
      do i=1,ele_count(new_positions)
        universal_element_number = fetch(new_elements, i)
        old_new_local_element_number = i
        new_new_local_element_number = renumber_permutation(old_new_local_element_number)
        call insert(uen_to_new_local_numbering, universal_element_number, new_new_local_element_number)
      end do

      do i=1,halo_count(new_positions)
        assert(halo_verifies(new_positions%mesh%halos(i), new_positions))
      end do

      ! Now cleanup
      ! Leave your 9's at home and bring your skills to the battle

      ierr = Zoltan_LB_Free_Part(export_global_ids, export_local_ids, export_procs, export_to_part)
      assert(ierr == ZOLTAN_OK)

      do i=1,size(sends)
        call deallocate(sends(i))
      end do

      ewrite(1,*) "Exiting reconstruct_halo"

    end subroutine reconstruct_halo

    subroutine deal_with_exports
      ! The unpack routine is the ones that do most of the heavy lifting
      ! in reconstructing the new mesh, positions, etc.
      ! But they don't get called if a process is only an exporter!
      ! So we need to reconstruct the positions and nelist for the
      ! nodes we own, and mark the request for the halos in the same way.
      ! (The owner of one of our halo nodes might have changed, and we
      ! don't know that, or who owns it, so we have to ask the old owner,
      ! i.e. participate in the phase_two migration).
      integer :: i, j, old_local_number, new_local_number
      integer, dimension(:), pointer :: neighbours
      integer :: universal_number
      type(integer_set) :: halo_nodes_we_need_to_know_about
      type(integer_hash_table) :: universal_number_to_old_owner
      logical :: changed
      integer :: old_owner, new_owner
      type(mesh_type) :: new_mesh
      type(integer_set) :: halo_nodes_we_currently_own
      integer :: rank

      if (associated(new_positions%refcount)) return
      ewrite(1,*) "In deal_with_exports"

      call allocate(new_nodes)

      ! Need to allocate new_nodes, new_elements, new_positions, new_nelist, universal_to_new_local_numbering
      do i=1,key_count(nodes_we_are_keeping)
        old_local_number = fetch(nodes_we_are_keeping, i)
        call insert(new_nodes, halo_universal_number(zz_halo, old_local_number), changed=changed)
      end do

      call allocate(halo_nodes_we_need_to_know_about)
      call allocate(halo_nodes_we_currently_own)
      call allocate(universal_number_to_old_owner)

      rank = getrank()

      do i=1,key_count(nodes_we_are_keeping)
        old_local_number = fetch(nodes_we_are_keeping, i)
        neighbours => row_m_ptr(zz_sparsity_two, old_local_number)
        do j=1,size(neighbours)
          universal_number = halo_universal_number(zz_halo, neighbours(j))
          call insert(new_nodes, universal_number, changed=changed)
          if (changed) then
            old_owner = halo_node_owner(zz_halo, neighbours(j)) - 1
            if (old_owner == rank) then
              call insert(halo_nodes_we_currently_own, neighbours(j))
            else
              call insert(halo_nodes_we_need_to_know_about, universal_number)
            end if
            call insert(universal_number_to_old_owner, universal_number, old_owner)
          end if
        end do
      end do

      ! We need to process these too -- nodes we own now, but will
      ! not own later, and will become our halo nodes.
      do i=1,key_count(halo_nodes_we_currently_own)
        universal_number = halo_universal_number(zz_halo, fetch(halo_nodes_we_currently_own, i))
        call insert(new_nodes, universal_number)
      end do

      call invert_set(new_nodes, universal_to_new_local_numbering)
      call allocate(new_mesh, key_count(new_nodes), 0, zz_mesh%shape, trim(zz_mesh%name))
      new_mesh%option_path = zz_mesh%option_path
      call allocate(new_positions, zz_positions%dim, new_mesh, trim(zz_positions%name))
      new_positions%option_path = zz_positions%option_path
      call deallocate(new_mesh)
      allocate(new_snelist(key_count(new_nodes)))
      call allocate(new_snelist)
      call allocate(new_surface_elements)

      allocate(new_nelist(key_count(new_nodes)))
      do i=1,key_count(new_nodes)
        call allocate(new_nelist(i))
      end do
      call allocate(new_elements)

      do i=1,key_count(nodes_we_are_keeping)
        old_local_number = fetch(nodes_we_are_keeping, i)
        universal_number = halo_universal_number(zz_halo, old_local_number)
        new_local_number = fetch(universal_to_new_local_numbering, universal_number)
        call set(new_positions, new_local_number, node_val(zz_positions, old_local_number))

        neighbours => row_m_ptr(zz_nelist, old_local_number)
        do j=1,size(neighbours)
          call insert(new_nelist(new_local_number), halo_universal_number(zz_ele_halo, neighbours(j)))
          call insert(new_elements, halo_universal_number(zz_ele_halo, neighbours(j)))
        end do

        ! and record the snelist information
        do j=1,key_count(old_snelist(old_local_number))
          call insert(new_snelist(new_local_number), fetch(old_snelist(old_local_number), j))
          call insert(new_surface_elements, fetch(old_snelist(old_local_number), j))
          ! we don't need to add anything to the universal_surface_number_to_surface_id because we already have it
        end do
      end do

      ! Set the positions and nelist of halo_nodes_we_currently_own
      do i=1,key_count(halo_nodes_we_currently_own)
        old_local_number = fetch(halo_nodes_we_currently_own, i)
        universal_number = halo_universal_number(zz_halo, old_local_number)
        new_local_number = fetch(universal_to_new_local_numbering, universal_number)
        call set(new_positions, new_local_number, node_val(zz_positions, old_local_number))

        neighbours => row_m_ptr(zz_nelist, old_local_number)
        do j=1,size(neighbours)
          call insert(new_nelist(new_local_number), halo_universal_number(zz_ele_halo, neighbours(j)))
          call insert(new_elements, halo_universal_number(zz_ele_halo, neighbours(j)))
        end do

        new_owner = fetch(nodes_we_are_sending, old_local_number)
        call insert(receives(new_owner+1), universal_number)

        ! and record the snelist information
        do j=1,key_count(old_snelist(old_local_number))
          call insert(new_snelist(new_local_number), fetch(old_snelist(old_local_number), j))
          call insert(new_surface_elements, fetch(old_snelist(old_local_number), j))
        end do
      end do
      call deallocate(halo_nodes_we_currently_own)

      my_num_import = key_count(halo_nodes_we_need_to_know_about)
      allocate(my_import_procs(my_num_import))
      allocate(my_import_global_ids(my_num_import))
      do i=1,my_num_import
        universal_number = fetch(halo_nodes_we_need_to_know_about, i)
        my_import_global_ids(i) = universal_number
        my_import_procs(i) = fetch(universal_number_to_old_owner, universal_number)
        assert(my_import_procs(i) /= getrank())
      end do

      call deallocate(halo_nodes_we_need_to_know_about)
      call deallocate(universal_number_to_old_owner)

    end subroutine deal_with_exports

    subroutine dump_linear_mesh
      type(scalar_field) :: sends, receives, unn
      integer :: i, proc

      assert(associated(new_positions%refcount))
      assert(new_positions%refcount%count == 1)
      assert(associated(new_positions%mesh%refcount))
      assert(new_positions%mesh%refcount%count == 1)

      call allocate(sends, new_positions%mesh, "Sends")
      call zero(sends)
      call allocate(receives, new_positions%mesh, "Receives")
      call zero(receives)
      call allocate(unn, new_positions%mesh, "NewUniversalNodeNumber")
      call zero(unn)

      do proc=1,halo_proc_count(new_positions%mesh%halos(2))
        do i=1,size(new_positions%mesh%halos(2)%sends(proc)%ptr)
          call set(sends, new_positions%mesh%halos(2)%sends(proc)%ptr(i), 1.0)
        end do
        do i=1,size(new_positions%mesh%halos(2)%receives(proc)%ptr)
          call set(receives, new_positions%mesh%halos(2)%receives(proc)%ptr(i), 1.0)
        end do
      end do

      do i=1,node_count(new_positions)
        call set(unn, i, float(halo_universal_number(new_positions%mesh%halos(2), i)))
      end do
      
!      call vtk_write_fields("balanced_linear_mesh", index=dumpno, position=new_positions, model=new_positions%mesh, sfields=(/sends, receives, unn/))
!      call vtk_write_surface_mesh("output_surface_mesh", dumpno, new_positions)

      call deallocate(sends)
      call deallocate(receives)
      call deallocate(unn)
    end subroutine dump_linear_mesh

    subroutine reconstruct_senlist
      type(integer_set), dimension(key_count(new_surface_elements)) :: senlists
      integer :: i, j, expected_loc, full_elements
      integer :: universal_number, new_local_number
      type(integer_hash_table) :: universal_surface_element_to_local_numbering
      integer, dimension(:), allocatable, target :: sndgln, surface_ids, element_owners
      type(csr_sparsity), pointer :: nnlist
      integer, dimension(:), pointer :: neighbours
      integer :: k, l
      logical, dimension(key_count(new_surface_elements)) :: keep_surface_element
      integer :: universal_element_number

      ewrite(1,*) "In reconstruct_senlist"

      ! new_surface_elements currently contains the universal numbers of surface elements
      ! we don't fully have and won't be in the final mesh.
      ! So universal_surface_element_to_local_numbering here is just temporary.
      call invert_set(new_surface_elements, universal_surface_element_to_local_numbering)

      do i=1,key_count(new_surface_elements)
        call allocate(senlists(i))
      end do

      ! Invert the snelists to give the senlists

      do i=1,key_count(new_nodes)
        do j=1,key_count(new_snelist(i))
          universal_number = fetch(new_snelist(i), j)
          new_local_number = fetch(universal_surface_element_to_local_numbering, universal_number)
          call insert(senlists(new_local_number), i)
        end do
      end do

      call deallocate(universal_surface_element_to_local_numbering)

      ! Now, some of these will be degenerate, because the halo nodes will refer
      ! to elements we don't know about. We can tell these apart because they
      ! are incomplete.

      full_elements = 0
      ! For mixed meshes, this should taken from the loc of the positions mesh
      ! However, it's constant for now
      expected_loc = face_loc(zz_mesh, 1)

      ! First, count how many we have
      nnlist => extract_nnlist(new_positions%mesh)
      do i=1,key_count(new_surface_elements)
        j = key_count(senlists(i))
        assert(j <= expected_loc)
        if (j == expected_loc) then
          ! We also need to check if we have the parent volume element --
          ! it's possible to get all the information for a face, without having the
          ! corresponding element!
          universal_number = fetch(new_surface_elements, i)
          universal_element_number = fetch(universal_surface_number_to_element_owner, universal_number)
          keep_surface_element(i) = has_key(uen_to_new_local_numbering, universal_element_number)
          if (keep_surface_element(i)) full_elements = full_elements + 1
        else
          keep_surface_element(i) = .false.
!          write(0,*) "Surface element ", fetch(new_surface_elements, i), " is degenerate. Dropping .."
!          write(0,*) "Local nodes: ", set2vector(senlists(i))
        end if 
      end do

      ewrite(2,*) "Found ", key_count(new_surface_elements), " possible new surface elements."
      ewrite(2,*) "Of these, ", full_elements, " are non-degenerate."

      ! And now fill in the non-degenerate ones

      allocate(sndgln(full_elements * expected_loc))
      allocate(surface_ids(full_elements))
      allocate(element_owners(full_elements))

      j = 1
      do i=1,key_count(new_surface_elements)
        if (keep_surface_element(i)) then
          universal_number = fetch(new_surface_elements, i)
          sndgln( (j-1)*expected_loc+1 : j*expected_loc ) = set2vector(senlists(i))
          surface_ids(j) = fetch(universal_surface_number_to_surface_id, universal_number)
          universal_element_number = fetch(universal_surface_number_to_element_owner, universal_number)
          element_owners(j) = fetch(uen_to_new_local_numbering, universal_element_number)
          j = j + 1
        end if
      end do
      assert(j == full_elements + 1)

      call add_faces(new_positions%mesh, sndgln=sndgln, boundary_ids=surface_ids, element_owner=element_owners)

      do i=1,size(senlists)
        call deallocate(senlists(i))
      end do

      ! New elements is no longer valid, as we have lost the degenerate elements
      call deallocate(new_surface_elements)
      call deallocate(universal_surface_number_to_surface_id)
      call deallocate(universal_surface_number_to_element_owner)
      call deallocate(new_snelist)
      deallocate(new_snelist)

      deallocate(sndgln)
      deallocate(surface_ids)
      deallocate(element_owners)

      call deinterleave_surface_ids(new_positions%mesh, max_coplanar_id)

      ! Bingo! Our mesh has an senlist.
      ewrite(1,*) "Exiting reconstruct_senlist"
    end subroutine reconstruct_senlist

    subroutine reconstruct_enlist
      type(integer_set), dimension(key_count(new_elements)) :: enlists
      integer :: i, j, expected_loc, full_elements
      integer :: universal_number, new_local_number
      type(integer_set) :: new_elements_we_actually_have

      ewrite(1,*) "In reconstruct_enlist"

      ! new_elements currently contains the universal numbers of elements
      ! we don't fully have and won't be in the final mesh.
      ! So uen_to_new_local_numbering here is just temporary -- we will
      ! construct the proper version later.
      call invert_set(new_elements, uen_to_new_local_numbering)

      do i=1,key_count(new_elements)
        call allocate(enlists(i))
      end do

      ! Invert the nelists to give the enlists

      do i=1,key_count(new_nodes)
        do j=1,key_count(new_nelist(i))
          universal_number = fetch(new_nelist(i), j)
          new_local_number = fetch(uen_to_new_local_numbering, universal_number)
          call insert(enlists(new_local_number), i)
        end do
      end do

      call deallocate(new_nelist)
      deallocate(new_nelist)

      call deallocate(uen_to_new_local_numbering)

      ! Now, some of these will be degenerate, because the halo nodes will refer
      ! to elements we don't know about. We can tell these apart because they
      ! are incomplete.

      full_elements = 0
      ! For mixed meshes, this should be the loc of the positions mesh
      ! note we know the universal number -- fetch(new_elements, i)
      ! However, it's constant for now
      expected_loc = zz_mesh%shape%loc

      ! First, count how many we have
      do i=1,key_count(new_elements)
        if (key_count(enlists(i)) == expected_loc) then
          full_elements = full_elements + 1
!        else
!          write(0,*) "Element ", fetch(new_elements, i), " is degenerate. Dropping .."
        end if 
      end do

      ewrite(2,*) "Found ", key_count(new_elements), " possible new elements."
      ewrite(2,*) "Of these, ", full_elements, " are non-degenerate."

      ! And now fill in the non-degenerate ones

      call allocate(new_elements_we_actually_have)
      call allocate(uen_to_new_local_numbering)
      new_positions%mesh%elements = full_elements
      deallocate(new_positions%mesh%ndglno)
      allocate(new_positions%mesh%ndglno(full_elements * expected_loc))
#ifdef HAVE_MEMORY_STATS
      call register_allocation("mesh_type", "integer", full_elements * expected_loc, name=new_positions%mesh%name)
#endif
      j = 1
      do i=1,key_count(new_elements)
        if (key_count(enlists(i)) == expected_loc) then
          universal_number = fetch(new_elements, i)
          call set_ele_nodes(new_positions%mesh, j, set2vector(enlists(i)))
          call insert(new_elements_we_actually_have, universal_number)
          call insert(uen_to_new_local_numbering, universal_number, j)
          j = j + 1
        end if
      end do

      call reorder_element_numbering(new_positions, use_unns=enlists)

      do i=1,size(enlists)
        call deallocate(enlists(i))
      end do

      ! New elements is no longer valid, as we have lost the degenerate elements
      call deallocate(new_elements)
      new_elements = new_elements_we_actually_have

      ! Bingo! Our mesh has an enlist.
      ewrite(1,*) "Exiting reconstruct_enlist"
    end subroutine reconstruct_enlist

    subroutine zoltan_migration_phase_two
      integer(zoltan_int) :: ierr
      integer(zoltan_int), dimension(:), pointer :: import_to_part, export_to_part, export_procs
      integer(zoltan_int), dimension(:), pointer :: export_global_ids, export_local_ids, import_local_ids
      integer(zoltan_int) :: num_export

      ! Register the new callback functions for packing and unpacking
      ierr = Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE, zoltan_cb_pack_halo_node_sizes); assert(ierr == ZOLTAN_OK)
      ierr = Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_MULTI_FN_TYPE, zoltan_cb_pack_halo_nodes); assert(ierr == ZOLTAN_OK)
      ierr = Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE, zoltan_cb_unpack_halo_nodes); assert(ierr == ZOLTAN_OK)

      import_to_part => null()
      export_to_part => null()
      export_global_ids => null()
      export_local_ids => null()
      export_procs => null()
      num_export = -1

      if (.not. associated(my_import_procs)) then
        ! We have nothing else to receive, but we still need to take part in
        ! the communication. Set num_import to 0
        my_num_import = 0
        allocate(my_import_global_ids(0))
        allocate(my_import_procs(0))
      end if
      ewrite(1,*) "In zoltan_migration_phase_two; objects to import: ", my_num_import

      ! We should be able to do:
      ! import_local_ids => null()
      ! from my reading of the Zoltan docs. But it actually doesn't appear to be the case.
      allocate(import_local_ids(my_num_import))
      if (my_num_import > 0) then
        import_local_ids = 666
      end if

      assert(associated(my_import_global_ids))
      assert(associated(import_local_ids))
      assert(associated(my_import_procs))
      assert(all(my_import_procs >= 0))
      assert(all(my_import_procs < getnprocs()))
      ierr = Zoltan_Compute_Destinations(zz, my_num_import, my_import_global_ids, import_local_ids, my_import_procs, &
       & num_export, export_global_ids, export_local_ids, export_procs)
      assert(ierr == ZOLTAN_OK)
      ewrite(1,*) "In zoltan_migration_phase_two; objects to export: ", num_export

      ierr = Zoltan_Migrate(zz, my_num_import, my_import_global_ids, import_local_ids, my_import_procs, &
       & import_to_part, num_export, export_global_ids, export_local_ids, export_procs, export_to_part) 
      assert(ierr == ZOLTAN_OK)

      ierr = Zoltan_LB_Free_Part(export_global_ids, export_local_ids, export_procs, export_to_part)
      assert(ierr == ZOLTAN_OK)

      deallocate(import_local_ids)
    end subroutine zoltan_migration_phase_two

    subroutine deallocate_my_lists
      deallocate(my_import_global_ids)
      deallocate(my_import_procs)
    end subroutine deallocate_my_lists

    subroutine are_we_keeping_or_sending_nodes
      integer :: node, i
      integer, dimension(halo_nowned_nodes(zz_halo)) :: owned_nodes

      call allocate(nodes_we_are_sending)
      call allocate(nodes_we_are_keeping)

      do i=1,p1_num_export
        call insert(nodes_we_are_sending, p1_export_local_ids(i), p1_export_procs(i))
      end do

      call get_owned_nodes(zz_halo, owned_nodes)
      do i=1,size(owned_nodes)
        node = owned_nodes(i)
        if (.not. has_key(nodes_we_are_sending, node)) then
          call insert(nodes_we_are_keeping, node)
        end if
      end do
    end subroutine are_we_keeping_or_sending_nodes

    subroutine setup_module_variables
      integer :: nhalos
      integer, dimension(:), allocatable :: owned_nodes
      integer :: i, j, floc, eloc
      integer, dimension(:), allocatable :: sndgln
      integer :: old_element_number, universal_element_number, face_number, universal_surface_element_number
      type(mesh_type) :: pwc_mesh
      integer, dimension(:), pointer :: eles
      real :: qual
      integer, dimension(:), allocatable :: interleaved_surface_ids

      !call find_mesh_to_adapt(states(1), zz_mesh)
      zz_mesh = get_external_mesh(states)
      call incref(zz_mesh)
      if (zz_mesh%name=="CoordinateMesh") then
        zz_positions = extract_vector_field(states, "Coordinate")
      else
        zz_positions = extract_vector_field(states, trim(zz_mesh%name)//"Coordinate")
      end if
      call incref(zz_positions)
        
      zz_nelist => extract_nelist(zz_mesh)

      zz => Zoltan_Create(halo_communicator(zz_mesh))

      nhalos = halo_count(zz_mesh)
      assert(nhalos == 2)
      zz_halo => zz_mesh%halos(nhalos)

      nhalos = element_halo_count(zz_mesh)
      assert(nhalos >= 1)
      zz_ele_halo => zz_mesh%element_halos(nhalos)

      zz_sparsity_one => get_csr_sparsity_firstorder(states, zz_mesh, zz_mesh)
      zz_sparsity_two => get_csr_sparsity_secondorder(states, zz_mesh, zz_mesh)

      allocate(owned_nodes(halo_nowned_nodes(zz_halo)))
      call allocate(universal_to_old_local_numbering)
      call get_owned_nodes(zz_halo, owned_nodes)
      do i=1,size(owned_nodes)
        call insert(universal_to_old_local_numbering, halo_universal_number(zz_halo, owned_nodes(i)), owned_nodes(i))
      end do
      deallocate(owned_nodes)

      call allocate(uen_to_old_local_numbering)
      do i=1,ele_count(zz_positions)
        call insert(uen_to_old_local_numbering, halo_universal_number(zz_ele_halo, i), i)
      end do

      allocate(receives(halo_proc_count(zz_halo)))
      do i=1,size(receives)
        call allocate(receives(i))
      end do

      ! set up old_snelist
      allocate(old_snelist(node_count(zz_positions)))
      call allocate(old_snelist)
      call allocate(universal_surface_number_to_surface_id)
      call allocate(universal_surface_number_to_element_owner)
      allocate(interleaved_surface_ids(surface_element_count(zz_positions)))
      call interleave_surface_ids(zz_mesh, interleaved_surface_ids, max_coplanar_id)

      ! this is another thing that needs to be generalised for mixed meshes
      floc = face_loc(zz_positions, 1)
      eloc = ele_loc(zz_positions, 1)
      allocate(sndgln(surface_element_count(zz_positions) * floc))
      call getsndgln(zz_mesh, sndgln)

      do i=1,surface_element_count(zz_positions)
        old_element_number = face_ele(zz_positions, i)
        universal_element_number = halo_universal_number(zz_ele_halo, old_element_number)
        face_number = local_face_number(zz_positions, i)
        universal_surface_element_number = (universal_element_number-1)*eloc + face_number

        call insert(universal_surface_number_to_surface_id, universal_surface_element_number, interleaved_surface_ids(i))
        call insert(universal_surface_number_to_element_owner, universal_surface_element_number, universal_element_number)

        do j=(i-1)*floc+1,i*floc
          call insert(old_snelist(sndgln(j)), universal_surface_element_number)
        end do
      end do

      deallocate(sndgln)
      deallocate(interleaved_surface_ids)

      ! And the element quality measure
      if (present(metric)) then
        call element_quality_p0(zz_positions, metric, element_quality)
      else
        pwc_mesh = piecewise_constant_mesh(zz_mesh, "PWCMesh")
        call allocate(element_quality, pwc_mesh, "ElementQuality", field_type=FIELD_TYPE_CONSTANT)
        call set(element_quality, 1.0)
        call deallocate(pwc_mesh)
      end if

      call allocate(node_quality, zz_mesh, "NodeQuality")
      call zero(node_quality)
      do i=1,node_count(node_quality)
        eles => row_m_ptr(zz_nelist, i)  
        qual = 1.0
        do j=1,size(eles)
          qual = min(qual, node_val(element_quality, eles(j)))
        end do
        call set(node_quality, i, qual)
      end do
      call halo_update(node_quality)

    end subroutine setup_module_variables

    subroutine set_zoltan_parameters
      integer(zoltan_int) :: ierr

      if (debug_level()>1) then
         ierr = Zoltan_Set_Param(zz, "DEBUG_LEVEL", "1"); assert(ierr == ZOLTAN_OK)
      else         
         ierr = Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0"); assert(ierr == ZOLTAN_OK)
      end if
      ierr = Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH"); assert(ierr == ZOLTAN_OK)
      if (iteration == 1) then
        ierr = Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION"); assert(ierr == ZOLTAN_OK)
      else
        ierr = Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION"); assert(ierr == ZOLTAN_OK)
      end if

      ierr = Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); assert(ierr == ZOLTAN_OK)
      ierr = Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1"); assert(ierr == ZOLTAN_OK)
      ierr = Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1"); assert(ierr == ZOLTAN_OK)
      ierr = Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "1"); assert(ierr == ZOLTAN_OK)
      ierr = Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); assert(ierr == ZOLTAN_OK)

      ierr = Zoltan_Set_Fn(zz, ZOLTAN_NUM_OBJ_FN_TYPE, zoltan_cb_owned_node_count);      assert(ierr == ZOLTAN_OK)
      ierr = Zoltan_Set_Fn(zz, ZOLTAN_OBJ_LIST_FN_TYPE, zoltan_cb_get_owned_nodes);      assert(ierr == ZOLTAN_OK)
      ierr = Zoltan_Set_Fn(zz, ZOLTAN_NUM_EDGES_MULTI_FN_TYPE, zoltan_cb_get_num_edges); assert(ierr == ZOLTAN_OK)
      ierr = Zoltan_Set_Fn(zz, ZOLTAN_EDGE_LIST_MULTI_FN_TYPE, zoltan_cb_get_edge_list); assert(ierr == ZOLTAN_OK)
      ierr = Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE, zoltan_cb_pack_node_sizes); assert(ierr == ZOLTAN_OK)
      ierr = Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_MULTI_FN_TYPE, zoltan_cb_pack_nodes); assert(ierr == ZOLTAN_OK)
      ierr = Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE, zoltan_cb_unpack_nodes); assert(ierr == ZOLTAN_OK)
    end subroutine set_zoltan_parameters

    subroutine deallocate_zoltan_lists
      integer(zoltan_int) :: ierr
      ierr = Zoltan_LB_Free_Data(p1_import_global_ids, p1_import_local_ids, p1_import_procs, &
                               & p1_export_global_ids, p1_export_local_ids, p1_export_procs)
      assert(ierr == ZOLTAN_OK)
    end subroutine deallocate_zoltan_lists

    subroutine cleanup_basic_module_variables
    ! This routine deallocates everything that is guaranteed to be allocated
    ! regardless of whether Zoltan actually wants to change anything or not.
      call deallocate(element_quality)
      call deallocate(node_quality)
      call deallocate(universal_surface_number_to_surface_id)
      call deallocate(universal_surface_number_to_element_owner)
      call deallocate(old_snelist)
      deallocate(old_snelist)
      call deallocate(receives)
      deallocate(receives)

      call deallocate(universal_to_old_local_numbering)
      call deallocate(uen_to_old_local_numbering)

      new_positions%refcount => null()

      call deallocate(zz_mesh)
      call deallocate(zz_positions)
      zz_sparsity_one => null()
      zz_sparsity_two => null()
      zz_halo => null()
      zz_ele_halo => null()
      call Zoltan_Destroy(zz)
    end subroutine cleanup_basic_module_variables

    subroutine cleanup_other_module_variables
      call deallocate(nodes_we_are_sending)
      call deallocate(nodes_we_are_keeping)
      call deallocate(universal_to_new_local_numbering)
      call deallocate(new_nodes)
      call deallocate(uen_to_new_local_numbering)
    end subroutine cleanup_other_module_variables

    subroutine zoltan_load_balance
      integer(zoltan_int) :: ierr
      integer :: i, node

      ! import_* aren't used because we set RETURN_LISTS to be only EXPORT

      ierr = Zoltan_LB_Balance(zz, changes, num_gid_entries, num_lid_entries, p1_num_import, p1_import_global_ids, &
          &    p1_import_local_ids, p1_import_procs, p1_num_export, p1_export_global_ids, p1_export_local_ids, p1_export_procs)
      assert(ierr == ZOLTAN_OK)

      do i=1,p1_num_export
        node = p1_export_local_ids(i)
        assert(node_owned(zz_halo, node))
      end do
    end subroutine zoltan_load_balance

    subroutine dump_suggested_owner
      integer :: rank, i
      type(scalar_field) :: suggested_owner, unn
      type(vector_field) :: positions

      rank = getrank()
      call allocate(suggested_owner, zz_mesh, "SuggestedOwner")

      call set(suggested_owner, float(rank))
      do i=1,p1_num_export
        call set(suggested_owner, p1_export_local_ids(i), float(p1_export_procs(i)))
      end do

      call allocate(unn, zz_mesh, "OldUniversalNodeNumber")
      do i=1,node_count(unn)
        call set(unn, i, float(halo_universal_number(zz_halo, i)))
      end do

      positions = get_coordinate_field(states(1), zz_mesh)
      call halo_update(suggested_owner)
!      call vtk_write_fields("suggested_owner", index=dumpno, position=positions, model=zz_mesh, sfields=(/suggested_owner, node_quality, element_quality, unn/))
      call deallocate(positions)
      call deallocate(suggested_owner)
      call deallocate(unn)

!      call vtk_write_surface_mesh("input_surface_mesh", dumpno, zz_positions)
    end subroutine dump_suggested_owner

    subroutine zoltan_migration_phase_one
      integer(zoltan_int) :: ierr
      integer(zoltan_int), dimension(:), pointer :: import_to_part, export_to_part

      ewrite(1,*) "In zoltan_migration_phase_one; objects to import: ", p1_num_import
      ewrite(1,*) "In zoltan_migration_phase_one; objects to export: ", p1_num_export
      import_to_part => null()
      export_to_part => null()

      ierr = Zoltan_Migrate(zz, p1_num_import, p1_import_global_ids, p1_import_local_ids, p1_import_procs, &
       & import_to_part, p1_num_export, p1_export_global_ids, p1_export_local_ids, p1_export_procs, export_to_part) 
      assert(ierr == ZOLTAN_OK)
    end subroutine zoltan_migration_phase_one

  end subroutine zoltan_drive

#endif
end module zoltan_integration
