#include "fdebug.h"
#include "confdefs.h"

module zoltan_callbacks

#ifdef HAVE_ZOLTAN

  use zoltan
  use zoltan_global_variables

  ! Needed for zoltan_cb_owned_node_count
  use halos, only: halo_nowned_nodes, halo_node_owner, halo_node_owners, get_owned_nodes, halo_universal_number
  use metric_tools

  ! Needed for zoltan_cb_get_owned_nodes
  ! - added get_owned_nodes, halo_universal_number to use halos
  use sparse_tools, only: row_length

  ! Needed for zoltan_cb_get_num_edges
  use global_parameters, only: real_size, OPTION_PATH_LEN
  use parallel_tools, only: getrank, getnprocs

  ! Needed for zoltan_cb_get_edge_list
  ! - added halo_node_owners to use halos
  use mpi_interfaces

  ! Needed for zoltan_cb_pack_node_sizes
  ! - added real_size to use global_parameters
  use data_structures

  ! Needed for zoltan_cb_pack_nodes
  ! - use the whole of data structures now

  implicit none
  
  public
  
contains

  function zoltan_cb_owned_node_count(data, ierr) result(count)
    integer(zoltan_int) :: count
    integer(zoltan_int), dimension(*) :: data ! not used
    integer(zoltan_int), intent(out) :: ierr
    
    ewrite(1,*) "In zoltan_cb_owned_node_count"

    count = halo_nowned_nodes(zoltan_global_zz_halo)
    if (have_option("/mesh_adaptivity/hr_adaptivity/zoltan_options/zoltan_debug")) then
       ewrite(1,*) "zoltan_cb_owned_node_count found: ", count, " nodes"
    end if
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
    real(zoltan_float) :: max_obj_wgt
    
    ewrite(1,*) "In zoltan_cb_get_owned_nodes"
    
    assert(num_gid_entries == 1)
    assert(num_lid_entries == 1)
    assert(wgt_dim == 1)
    
    count = halo_nowned_nodes(zoltan_global_zz_halo)
    
    call get_owned_nodes(zoltan_global_zz_halo, local_ids(1:count))
    global_ids(1:count) = halo_universal_number(zoltan_global_zz_halo, local_ids(1:count))
    
    if (have_option("/mesh_adaptivity/hr_adaptivity/zoltan_options/zoltan_debug")) then
       ewrite(1,*) "zoltan_cb_get_owned nodes found local_ids: ", local_ids(1:count)
       ewrite(1,*) "zoltan_cb_get_owned nodes found global_ids: ", global_ids(1:count)
    end if
    
    if(zoltan_global_migrate_extruded_mesh) then
       ! weight the nodes according to the number of nodes in the column beneath it
       max_obj_wgt = 1.0
       do i=1,count
          obj_wgts(i) = float(row_length(zoltan_global_columns_sparsity, i))
          max_obj_wgt = max(max_obj_wgt, obj_wgts(i))
       end do
       ! normalise according to the most nodes in a column
       do i=1,count
          obj_wgts(i) = obj_wgts(i)/max_obj_wgt
       end do
    else
       do i=1,count
          obj_wgts(i) = 1.0
       end do
    end if
    
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
    character (len = OPTION_PATH_LEN) :: filename

    ewrite(1,*) "In zoltan_cb_get_num_edges"

    assert(num_gid_entries == 1)
    assert(num_lid_entries == 1)

    count = zoltan_global_zz_halo%nowned_nodes
    assert(count == num_obj)

    do node=1,count
      num_edges(node) = row_length(zoltan_global_zz_sparsity_one, local_ids(node))
    end do

    if (have_option("/mesh_adaptivity/hr_adaptivity/zoltan_options/zoltan_debug/dump_edge_counts")) then      
       write(filename, '(A,I0,A)') 'edge_counts_', getrank(),'.dat'
       open(666, file = filename)
       do node=1,count
          write(666,*) num_edges(node)
       end do
       close(666)
    end if

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
    integer :: count, err
    integer :: node, i, j
    integer :: head
    integer, dimension(:), pointer :: neighbours
    character (len = OPTION_PATH_LEN) :: filename    
    
    ! variables for recording various element quality functional values 
    real(zoltan_float) :: quality, min_quality, my_min_quality
    
    ! variables for recording the local maximum/minimum edge weights and local 90th percentile edge weight
    real(zoltan_float) :: min_weight, max_weight, ninety_weight, my_max_weight
    
    integer, dimension(:), pointer :: my_nelist, nbor_nelist
    
    integer :: total_num_edges, my_num_edges
    
    real :: value
    
    ewrite(1,*) "In zoltan_cb_get_edge_list"
    
    assert(num_gid_entries == 1)
    assert(num_lid_entries == 1)
    assert(wgt_dim == 1)
    
    count = zoltan_global_zz_halo%nowned_nodes
    assert(count == num_obj)
    
    my_num_edges = sum(num_edges(1:num_obj))
    
    if (zoltan_global_zoltan_iteration==zoltan_global_zoltan_max_adapt_iteration) then
       
       ! last iteration - hopefully the mesh is of sufficient quality by now
       ! we only want to optimize the edge cut to minimize halo communication
       ewgts(1:my_num_edges) = 1.0       
       head = 1
       do node=1,count
          ! find nodes neighbours
          neighbours => row_m_ptr(zoltan_global_zz_sparsity_one, local_ids(node))
          ! check the number of neighbours matches the number of edges
          assert(size(neighbours) == num_edges(node))
          ! find global ids for each neighbour
          nbor_global_id(head:head+size(neighbours)-1) = halo_universal_number(zoltan_global_zz_halo, neighbours)
          ! find owning proc for each neighbour
          nbor_procs(head:head+size(neighbours)-1) = halo_node_owners(zoltan_global_zz_halo, neighbours) - 1
          head = head + size(neighbours)
       end do
       ierr = ZOLTAN_OK
       return
    else
       call MPI_ALLREDUCE(my_num_edges,total_num_edges,1,MPI_INTEGER,MPI_SUM, &
            MPI_COMM_WORLD,err)
    end if
    
    
    head = 1
    
    ! Aim is to assign high edge weights to poor quality elements
    ! so that when we load balance poor quality elements are placed
    ! in the centre of partitions and can be adapted
    
    ! loop over the nodes you own
    do node=1,count
       
       ! find nodes neighbours
       neighbours => row_m_ptr(zoltan_global_zz_sparsity_one, local_ids(node))
       
       ! check the number of neighbours matches the number of edges
       assert(size(neighbours) == num_edges(node))
       
       ! find global ids for each neighbour
       nbor_global_id(head:head+size(neighbours)-1) = halo_universal_number(zoltan_global_zz_halo, neighbours)
       
       ! find owning proc for each neighbour
       nbor_procs(head:head+size(neighbours)-1) = halo_node_owners(zoltan_global_zz_halo, neighbours) - 1
       
       ! get elements associated with current node
       my_nelist => row_m_ptr(zoltan_global_zz_nelist, local_ids(node))
       
       my_min_quality = 1.0
       
       ! find quality of worst element node is associated with
       do i=1,size(my_nelist)
          quality = minval(ele_val(zoltan_global_element_quality, my_nelist(i)))
          
          if (quality .LT. my_min_quality) then
             my_min_quality = quality
          end if
       end do
       
       ! loop over all neighbouring nodes
       do j=1,size(neighbours)
          
          min_quality = my_min_quality
          
          ! get elements associated with neighbour node
          nbor_nelist => row_m_ptr(zoltan_global_zz_nelist, neighbours(j))
          
          ! loop over all the elements of the neighbour node
          do i=1, size(nbor_nelist)
             ! determine the quality of the element
             quality = minval(ele_val(zoltan_global_element_quality, nbor_nelist(i)))
             
             ! store the element quality if it's less (worse) than any previous elements
             if (quality .LT. min_quality) then
                min_quality = quality
             end if
          end do
          
          ! check if the quality is within the tolerance         
          if (min_quality .GT. zoltan_global_quality_tolerance) then
             ! if it is
             ewgts(head + j - 1) = 1.0
          else
             ! if it's not
             ewgts(head + j - 1) = ceiling((1.0 - min_quality) * 20)
          end if
       end do
       
       value = maxval(ewgts(head:head+size(neighbours)-1))
       head = head + size(neighbours)
    end do
    
    assert(head == sum(num_edges(1:num_obj))+1)
    
    ! calculate the local maximum edge weight
    my_max_weight = maxval(ewgts(1:head-1))
    
    ! calculate the local minimum edge weight
    min_weight = minval(ewgts(1:head-1))
    call MPI_ALLREDUCE(my_max_weight,max_weight,1,MPI_INTEGER,MPI_MAX, MPI_COMM_WORLD,err)
    ! calculate the local 90th percentile edge weight   
    ninety_weight = max_weight * 0.90
    
    ! don't want to adjust the weights if all the elements are of a similar quality
    if (min_weight < ninety_weight) then
       ! make the worst 10% of elements uncuttable
       do i=1,head-1
          if (ewgts(i) .GT. ninety_weight) then
             ewgts(i) = total_num_edges + 1
          end if
       end do
    end if
    
    if (zoltan_global_output_edge_weights) then
       head = 1
       do node=1,count
          neighbours => row_m_ptr(zoltan_global_zz_sparsity_one, local_ids(node))
          value = maxval(ewgts(head:head+size(neighbours)-1))
          call set(zoltan_global_max_edge_weight_on_node,local_ids(node),value)
          head = head + size(neighbours)
       end do
    end if
    
    if (have_option("/mesh_adaptivity/hr_adaptivity/zoltan_options/zoltan_debug/dump_edge_weights")) then
       write(filename, '(A,I0,A)') 'edge_weights_', getrank(),'.dat'
       open(666, file = filename)
       do i=1,head-1
          write(666,*) ewgts(i)
       end do
       close(666)
    end if
    
    
    ierr = ZOLTAN_OK
  end subroutine zoltan_cb_get_edge_list


  ! Here is how we pack nodal positions for phase one migration:
  ! --------------------------------------------------------------------------------------------------------------
  ! | position | sz of lv-1 nnlist | lv-1 nnlist | sz of lv-2 nnlist | lv-2 nnlist | owners of level-2 nnlist | 
  ! | sz of nelist | nelist | sz of snelist | snelist | snelist ids | containing element of snelist |
  ! --------------------------------------------------------------------------------------------------------------
  subroutine zoltan_cb_pack_node_sizes(data, num_gid_entries,  num_lid_entries, num_ids, global_ids, local_ids, sizes, ierr) 
    integer(zoltan_int), dimension(*), intent(in) :: data 
    integer(zoltan_int), intent(in) :: num_gid_entries, num_lid_entries, num_ids
    integer(zoltan_int), intent(in), dimension(*) :: global_ids,  local_ids
    integer(zoltan_int), intent(out), dimension(*) :: sizes 
    integer(zoltan_int), intent(out) :: ierr  

    integer :: i, node
    character (len = OPTION_PATH_LEN) :: filename

    ewrite(1,*) "In zoltan_cb_pack_node_sizes"
    do i=1,num_ids
      node = local_ids(i)
      sizes(i) = zoltan_global_zz_positions%dim * real_size + &
                1 * integer_size + row_length(zoltan_global_zz_sparsity_one, node) * integer_size + &
                1 * integer_size + row_length(zoltan_global_zz_sparsity_two, node) * integer_size * 2 + &
                1 * integer_size + row_length(zoltan_global_zz_nelist, node) * integer_size + &
                1 * integer_size + key_count(zoltan_global_old_snelist(node)) * integer_size * 3
      if(zoltan_global_preserve_mesh_regions) then
        sizes(i) = sizes(i) + row_length(zoltan_global_zz_nelist, node) * integer_size
      end if
      if(zoltan_global_preserve_columns) then
        sizes(i) = sizes(i) + integer_size
      end if
    end do

    if (have_option("/mesh_adaptivity/hr_adaptivity/zoltan_options/zoltan_debug/dump_node_sizes")) then
       write(filename, '(A,I0,A)') 'node_sizes_', getrank(),'.dat'
       open(666, file = filename)
       do i=1,num_ids
          write(666,*) sizes(i)
       end do
       close(666)
    end if

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
      assert(halo_universal_number(zoltan_global_zz_halo, node) == global_ids(i))
      do j=1,zoltan_global_zz_positions%dim
        buf(head:head+ratio-1) = transfer(node_val(zoltan_global_zz_positions, j, node), buf(head:head+ratio-1))
        head = head + ratio
      end do
      
      if(zoltan_global_preserve_columns) then
        buf(head) = zoltan_global_universal_columns(node)
        head = head + 1
      end if
      
      buf(head) = row_length(zoltan_global_zz_sparsity_one, node)
      head = head + 1

      buf(head:head+row_length(zoltan_global_zz_sparsity_one, node)-1) = halo_universal_number(zoltan_global_zz_halo, row_m_ptr(zoltan_global_zz_sparsity_one, node))
      head = head + row_length(zoltan_global_zz_sparsity_one, node)

      buf(head) = row_length(zoltan_global_zz_sparsity_two, node)
      head = head + 1

      buf(head:head+row_length(zoltan_global_zz_sparsity_two, node)-1) = halo_universal_number(zoltan_global_zz_halo, row_m_ptr(zoltan_global_zz_sparsity_two, node))
      head = head + row_length(zoltan_global_zz_sparsity_two, node)

      buf(head:head+row_length(zoltan_global_zz_sparsity_two, node)-1) = halo_node_owners(zoltan_global_zz_halo, row_m_ptr(zoltan_global_zz_sparsity_two, node))
      head = head + row_length(zoltan_global_zz_sparsity_two, node)

      buf(head) = row_length(zoltan_global_zz_nelist, node)
      head = head + 1

      buf(head:head+row_length(zoltan_global_zz_nelist,node)-1) = halo_universal_number(zoltan_global_zz_ele_halo, row_m_ptr(zoltan_global_zz_nelist, node))
      head = head + row_length(zoltan_global_zz_nelist, node)

      if(zoltan_global_preserve_mesh_regions) then
        ! put in the region_ids in the same amount of space as the nelist - this is complete overkill!
        buf(head:head+row_length(zoltan_global_zz_nelist,node)-1) = fetch(zoltan_global_universal_element_number_to_region_id, halo_universal_number(zoltan_global_zz_ele_halo, row_m_ptr(zoltan_global_zz_nelist, node)))
        head = head + row_length(zoltan_global_zz_nelist, node)
      end if

      buf(head) = key_count(zoltan_global_old_snelist(node))
      head = head + 1
      buf(head:head + key_count(zoltan_global_old_snelist(node)) - 1) = set2vector(zoltan_global_old_snelist(node))
      head = head + key_count(zoltan_global_old_snelist(node))
      buf(head:head + key_count(zoltan_global_old_snelist(node)) - 1) = fetch(zoltan_global_universal_surface_number_to_surface_id, set2vector(zoltan_global_old_snelist(node)))
      head = head + key_count(zoltan_global_old_snelist(node))
      buf(head:head + key_count(zoltan_global_old_snelist(node)) - 1) = fetch(zoltan_global_universal_surface_number_to_element_owner, set2vector(zoltan_global_old_snelist(node)))
      head = head + key_count(zoltan_global_old_snelist(node))

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
    real, dimension(zoltan_global_zz_positions%dim) :: new_coord
    type(integer_hash_table) :: universal_number_to_old_owner
    integer :: old_owner, new_owner
    integer, dimension(:), pointer :: current_buf
    integer :: rank
    
    ewrite(1,*) "In zoltan_cb_unpack_nodes"
    ! assert new linear mesh and positions not allocated
    assert(.not. associated(new_mesh%refcount))
    assert(.not. associated(zoltan_global_new_positions%refcount))
    rank = getrank()

    ! Figure out the nodes we are going to know about
    call allocate(zoltan_global_new_nodes)

    ! All the nodes we currently have and still own, in universal numbering
    do i=1,key_count(zoltan_global_nodes_we_are_keeping)
       old_local_number = fetch(zoltan_global_nodes_we_are_keeping, i)
       call insert(zoltan_global_new_nodes, halo_universal_number(zoltan_global_zz_halo, old_local_number), changed=changed)
    end do

    ! All the nodes we are receiving from other people and are going to own
    do i=1,num_ids
       call insert(zoltan_global_new_nodes, global_ids(i))
    end do

    call allocate(universal_number_to_old_owner)
    call allocate(halo_nodes_we_currently_own)

    ! All the halos of (the nodes we currently have and still own), in universal numbering
    do i=1,key_count(zoltan_global_nodes_we_are_keeping)
       neighbours => row_m_ptr(zoltan_global_zz_sparsity_two, fetch(zoltan_global_nodes_we_are_keeping, i))
       do j=1,size(neighbours)
          universal_number = halo_universal_number(zoltan_global_zz_halo, neighbours(j))
          call insert(zoltan_global_new_nodes, universal_number, changed=changed)
          if (changed) then ! so it is a halo node
             old_owner = halo_node_owner(zoltan_global_zz_halo, neighbours(j)) - 1
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
       universal_number = halo_universal_number(zoltan_global_zz_halo, old_local_number)
       call insert(zoltan_global_new_nodes, universal_number)
    end do

    ! All the other nodes that will form the halo of the new nodes we are receiving 
    ratio = real_size / integer_size
    do i=1,num_ids
       current_buf => buf(idx(i):idx(i) + sizes(i)/integer_size)
       head = ratio * zoltan_global_zz_positions%dim + 1
       if(zoltan_global_preserve_columns) then
          head = head + 1
       end if
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
          call insert(zoltan_global_new_nodes, universal_number) ! the sparsity for a node includes itself
          call insert(universal_number_to_old_owner, universal_number, old_owner)
       end do
    end do

    ! Now zoltan_global_new_nodes implicitly defines a mapping
    ! between 1 .. key_count(zoltan_global_new_nodes) [these are the new local node numbers]
    ! and the universal node numbers of the nodes.
    ! We're going to invert that to create the hash table of universal node numbers -> local node numbers
    ! to facilitate the transfer of field information later.
    call invert_set(zoltan_global_new_nodes, zoltan_global_universal_to_new_local_numbering)

    ! allocate the new objects
    ! We know the number of nodes, but not the number of elements .. hmm.
    ! We will allocate it with 0 elements for now, and work it out when we
    ! invert the nnlist to compute an enlist later.
    call allocate(new_mesh, key_count(zoltan_global_new_nodes), 0, zoltan_global_zz_mesh%shape, trim(zoltan_global_zz_mesh%name))
    new_mesh%option_path = zoltan_global_zz_mesh%option_path
    if(zoltan_global_preserve_columns) then
       allocate(new_mesh%columns(key_count(zoltan_global_new_nodes)))
    end if
    call allocate(zoltan_global_new_positions, zoltan_global_zz_positions%dim, new_mesh, trim(zoltan_global_zz_positions%name))
    zoltan_global_new_positions%option_path = zoltan_global_zz_positions%option_path
    call deallocate(new_mesh)
    allocate(zoltan_global_new_snelist(key_count(zoltan_global_new_nodes)))
    call allocate(zoltan_global_new_snelist)

    ! aaaand unpack, recording which universal ids we have received
    ! so that we can figure out which ones we haven't received yet
    ! so that we can ask their old owner to send on the new details
    ! a) build a set of the nodes we have recorded
    ! b) from that, build a set of the nodes we haven't yet recorded
    ! c) figure out who owns those nodes, so we can build the import list for zoltan

    allocate(zoltan_global_new_nelist(key_count(zoltan_global_new_nodes)))
    do i=1,key_count(zoltan_global_new_nodes)
       call allocate(zoltan_global_new_nelist(i))
    end do
    call allocate(zoltan_global_new_elements)
    call allocate(zoltan_global_new_surface_elements)

    call allocate(new_nodes_we_have_recorded)
    ! Nodes we are keeping
    do i=1,key_count(zoltan_global_nodes_we_are_keeping)
       old_local_number = fetch(zoltan_global_nodes_we_are_keeping, i)
       universal_number = halo_universal_number(zoltan_global_zz_halo, old_local_number)
       new_local_number = fetch(zoltan_global_universal_to_new_local_numbering, universal_number)
       call insert(new_nodes_we_have_recorded, universal_number)
       call set(zoltan_global_new_positions, new_local_number, node_val(zoltan_global_zz_positions, old_local_number))

       if(zoltan_global_preserve_columns) then
          zoltan_global_new_positions%mesh%columns(new_local_number) = fetch(zoltan_global_universal_to_new_local_numbering_m1d, &
               zoltan_global_universal_columns(old_local_number))
       end if

       ! Record the nelist information
       neighbours => row_m_ptr(zoltan_global_zz_nelist, old_local_number)
       do j=1,size(neighbours)
          call insert(zoltan_global_new_nelist(new_local_number), halo_universal_number(zoltan_global_zz_ele_halo, neighbours(j)))
          call insert(zoltan_global_new_elements, halo_universal_number(zoltan_global_zz_ele_halo, neighbours(j)))
          ! don't need to do anything to zoltan_global_universal_element_number_to_region_id because we already have it
       end do

       ! and record the snelist information
       do j=1,key_count(zoltan_global_old_snelist(old_local_number))
          call insert(zoltan_global_new_snelist(new_local_number), fetch(zoltan_global_old_snelist(old_local_number), j))
          call insert(zoltan_global_new_surface_elements, fetch(zoltan_global_old_snelist(old_local_number), j))
          ! we don't need to add anything to the zoltan_global_universal_surface_number_to_surface_id because we already have it
       end do
    end do

    ! Set the positions and nelist of halo_nodes_we_currently_own
    do i=1,key_count(halo_nodes_we_currently_own)
       old_local_number = fetch(halo_nodes_we_currently_own, i)
       universal_number = halo_universal_number(zoltan_global_zz_halo, old_local_number)
       new_local_number = fetch(zoltan_global_universal_to_new_local_numbering, universal_number)
       call set(zoltan_global_new_positions, new_local_number, node_val(zoltan_global_zz_positions, old_local_number))

       if(zoltan_global_preserve_columns) then
          zoltan_global_new_positions%mesh%columns(new_local_number) = fetch(zoltan_global_universal_to_new_local_numbering_m1d, &
               zoltan_global_universal_columns(old_local_number))
       end if

       neighbours => row_m_ptr(zoltan_global_zz_nelist, old_local_number)
       do j=1,size(neighbours)
          call insert(zoltan_global_new_nelist(new_local_number), halo_universal_number(zoltan_global_zz_ele_halo, neighbours(j)))
          call insert(zoltan_global_new_elements, halo_universal_number(zoltan_global_zz_ele_halo, neighbours(j)))
          ! don't need to do anything to zoltan_global_universal_element_number_to_region_id because we already have it
       end do

       call insert(new_nodes_we_have_recorded, universal_number)

       new_owner = fetch(zoltan_global_nodes_we_are_sending, old_local_number)
       call insert(zoltan_global_receives(new_owner+1), universal_number)

       ! and record the snelist information
       do j=1,key_count(zoltan_global_old_snelist(old_local_number))
          call insert(zoltan_global_new_snelist(new_local_number), fetch(zoltan_global_old_snelist(old_local_number), j))
          call insert(zoltan_global_new_surface_elements, fetch(zoltan_global_old_snelist(old_local_number), j))
       end do
    end do
    call deallocate(halo_nodes_we_currently_own)

    ! Nodes we are gaining
    do i=1,num_ids
       call insert(new_nodes_we_have_recorded, global_ids(i))
       new_local_number = fetch(zoltan_global_universal_to_new_local_numbering, global_ids(i))
       new_coord = 0
       head = idx(i)
       do j=1,zoltan_global_zz_positions%dim
          new_coord(j) = transfer(buf(head:head+ratio-1), new_coord(j))
          head = head + ratio
       end do
       call set(zoltan_global_new_positions, new_local_number, new_coord)

       if(zoltan_global_preserve_columns) then
          zoltan_global_new_positions%mesh%columns(new_local_number) = fetch(zoltan_global_universal_to_new_local_numbering_m1d, buf(head)) 
          head = head + 1
       end if

       ! Record the nelist information
       sz = buf(head) ! level-1 nnlist
       head = head + sz + 1
       sz = buf(head) ! level-2 nnlist
       head = head + 2*sz + 1
       sz = buf(head) ! nelist
       do j=1,sz
          call insert(zoltan_global_new_nelist(new_local_number), buf(head + j))
          call insert(zoltan_global_new_elements, buf(head + j))
          if(zoltan_global_preserve_mesh_regions) then
             call insert(zoltan_global_universal_element_number_to_region_id, buf(head + j), buf(head + j + sz))
          end if
       end do
       if(zoltan_global_preserve_mesh_regions) then
          head = head + 2*sz + 1
       else
          head = head + sz + 1
       end if

       ! And record the snelist information
       sz = buf(head)
       do j=1,sz
          call insert(zoltan_global_new_snelist(new_local_number), buf(head + j))
          call insert(zoltan_global_universal_surface_number_to_surface_id, buf(head + j), buf(head + j + sz))
          call insert(zoltan_global_universal_surface_number_to_element_owner, buf(head + j), buf(head + j + 2*sz))
          call insert(zoltan_global_new_surface_elements, buf(head + j))
       end do
       head = head + 3*sz + 1
    end do

    ! At this point, there might still be nodes that we have not yet recorded but
    ! we own, so we can fill them in now.
    call allocate(new_nodes_we_still_need)
    do i=1,key_count(zoltan_global_new_nodes)
       universal_number = fetch(zoltan_global_new_nodes, i)
       if (has_value(new_nodes_we_have_recorded, universal_number)) cycle

       old_owner = fetch(universal_number_to_old_owner, universal_number)
       if (old_owner == rank) then
          call insert(new_nodes_we_have_recorded, universal_number)
          old_local_number = fetch(zoltan_global_universal_to_old_local_numbering, universal_number)
          new_local_number = fetch(zoltan_global_universal_to_new_local_numbering, universal_number)
          call set(zoltan_global_new_positions, new_local_number, node_val(zoltan_global_zz_positions, old_local_number))

          if(zoltan_global_preserve_columns) then
             zoltan_global_new_positions%mesh%columns(new_local_number) = fetch(zoltan_global_universal_to_new_local_numbering_m1d, &
                  zoltan_global_universal_columns(old_local_number))
          end if

          ! Record the nelist information
          neighbours => row_m_ptr(zoltan_global_zz_nelist, old_local_number)
          do j=1,size(neighbours)
             call insert(zoltan_global_new_nelist(new_local_number), halo_universal_number(zoltan_global_zz_ele_halo, neighbours(j)))
             call insert(zoltan_global_new_elements, halo_universal_number(zoltan_global_zz_ele_halo, neighbours(j)))
             ! don't need to do anything to zoltan_global_universal_element_number_to_region_id because we already have it
          end do

          ! and record the snelist information
          do j=1,key_count(zoltan_global_old_snelist(old_local_number))
             call insert(zoltan_global_new_snelist(new_local_number), fetch(zoltan_global_old_snelist(old_local_number), j))
             call insert(zoltan_global_new_surface_elements, fetch(zoltan_global_old_snelist(old_local_number), j))
             ! we don't need to add anything to the zoltan_global_universal_surface_number_to_surface_id because we already have it
          end do

          ! and record the node in the zoltan_global_receives
          new_owner = fetch(zoltan_global_nodes_we_are_sending, old_local_number)
          call insert(zoltan_global_receives(new_owner+1), universal_number)
       else
          call insert(new_nodes_we_still_need, universal_number)
       end if
    end do

    ! And build the import list ...
    zoltan_global_my_num_import = key_count(new_nodes_we_still_need)
    allocate(zoltan_global_my_import_procs(zoltan_global_my_num_import))
    allocate(zoltan_global_my_import_global_ids(zoltan_global_my_num_import))
    do i=1,zoltan_global_my_num_import
       universal_number = fetch(new_nodes_we_still_need, i)
       zoltan_global_my_import_global_ids(i) = universal_number
       zoltan_global_my_import_procs(i) = fetch(universal_number_to_old_owner, universal_number)
       assert(zoltan_global_my_import_procs(i) /= rank)
    end do

    call deallocate(new_nodes_we_have_recorded)
    call deallocate(new_nodes_we_still_need)
    call deallocate(universal_number_to_old_owner)

    ierr = ZOLTAN_OK
  end subroutine zoltan_cb_unpack_nodes


  ! Here is how we pack halo nodes for phase two migration:
  ! -------------------------------------------------------------------------------------
  ! | position | new owner | size of nelist | nelist | size of snelist | 
  ! | snelist | surface ids | the containing volume element for each surface element |
  ! -------------------------------------------------------------------------------------
  subroutine zoltan_cb_pack_halo_node_sizes(data, num_gid_entries,  num_lid_entries, num_ids, global_ids, local_ids, sizes, ierr) 
    integer(zoltan_int), dimension(*), intent(in) :: data 
    integer(zoltan_int), intent(in) :: num_gid_entries, num_lid_entries, num_ids
    integer(zoltan_int), intent(in), dimension(*) :: global_ids,  local_ids
    integer(zoltan_int), intent(out), dimension(*) :: sizes 
    integer(zoltan_int), intent(out) :: ierr  

    integer :: i, node
    character (len = OPTION_PATH_LEN) :: filename

    ewrite(1,*) "In zoltan_cb_pack_halo_node_sizes"

    do i=1,num_ids
      node = fetch(zoltan_global_universal_to_old_local_numbering, global_ids(i))
      sizes(i) = zoltan_global_zz_positions%dim * real_size + &
                2 * integer_size + row_length(zoltan_global_zz_nelist, node) * integer_size + &
                1 * integer_size + key_count(zoltan_global_old_snelist(node)) * 3 * integer_size
      if(zoltan_global_preserve_mesh_regions) then
        sizes(i) = sizes(i) + row_length(zoltan_global_zz_nelist, node) * integer_size 
      end if
      if(zoltan_global_preserve_columns) then
        sizes(i) = sizes(i) + integer_size
      end if
    end do

    if (have_option("/mesh_adaptivity/hr_adaptivity/zoltan_options/zoltan_debug/dump_halo_node_sizes")) then
       write(filename, '(A,I0,A)') 'halo_node_sizes_', getrank(),'.dat'
       open(666, file = filename)
       do i=1,num_ids
          write(666,*) sizes(i)
       end do
       close(666)
    end if

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
       node = fetch(zoltan_global_universal_to_old_local_numbering, global_ids(i))
       do j=1,zoltan_global_zz_positions%dim
          current_buf(head:head+ratio-1) = transfer(node_val(zoltan_global_zz_positions, j, node), current_buf(head:head+ratio-1))
          head = head + ratio
       end do
       
       if(zoltan_global_preserve_columns) then
          current_buf(head) = zoltan_global_universal_columns(node)
          head = head + 1
       end if
       
       ! Now compute the new owner
       if (has_key(zoltan_global_nodes_we_are_sending, node)) then
          new_owner = fetch(zoltan_global_nodes_we_are_sending, node)
       else
          new_owner = rank
       end if
       
       current_buf(head) = new_owner
       head = head + 1
       
       current_buf(head) = row_length(zoltan_global_zz_nelist, node)
       head = head + 1
       
       current_buf(head:head+row_length(zoltan_global_zz_nelist, node)-1) = halo_universal_number(zoltan_global_zz_ele_halo, row_m_ptr(zoltan_global_zz_nelist, node))
       head = head + row_length(zoltan_global_zz_nelist, node)
       
       if(zoltan_global_preserve_mesh_regions) then
          ! put in the region_ids in the same amount of space as the nelist - this is complete overkill!
          current_buf(head:head+row_length(zoltan_global_zz_nelist, node)-1) = fetch(zoltan_global_universal_element_number_to_region_id, halo_universal_number(zoltan_global_zz_ele_halo, row_m_ptr(zoltan_global_zz_nelist, node)))
          head = head + row_length(zoltan_global_zz_nelist, node)
       end if
       
       current_buf(head) = key_count(zoltan_global_old_snelist(node))
       head = head + 1
       current_buf(head:head+key_count(zoltan_global_old_snelist(node))-1) = set2vector(zoltan_global_old_snelist(node))
       head = head + key_count(zoltan_global_old_snelist(node))
       current_buf(head:head+key_count(zoltan_global_old_snelist(node))-1) = fetch(zoltan_global_universal_surface_number_to_surface_id, set2vector(zoltan_global_old_snelist(node)))
       head = head + key_count(zoltan_global_old_snelist(node))
       current_buf(head:head+key_count(zoltan_global_old_snelist(node))-1) = fetch(zoltan_global_universal_surface_number_to_element_owner, set2vector(zoltan_global_old_snelist(node)))
       head = head + key_count(zoltan_global_old_snelist(node))
       
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
    real, dimension(zoltan_global_zz_positions%dim) :: new_coord
    integer :: head
    integer :: ratio
    integer :: new_local_number, new_owner, sz
    
    ewrite(1,*) "In zoltan_cb_unpack_halo_nodes"
    
    ratio = real_size/integer_size
    
    do i=1,num_ids
       new_local_number = fetch(zoltan_global_universal_to_new_local_numbering, global_ids(i))
       new_coord = 0
       head = idx(i)
       do j=1,zoltan_global_zz_positions%dim
          new_coord(j) = transfer(buf(head:head+ratio-1), new_coord(j))
          head = head + ratio
       end do
       call set(zoltan_global_new_positions, new_local_number, new_coord)
       
       if(zoltan_global_preserve_columns) then
          zoltan_global_new_positions%mesh%columns(new_local_number) = fetch(zoltan_global_universal_to_new_local_numbering_m1d, buf(head))
          head = head + 1
       end if
       
       new_owner = buf(head)
       head = head + 1
       
       ! record the nelist information
       sz = buf(head)
       do j=1,sz
          call insert(zoltan_global_new_nelist(new_local_number), buf(head + j))
          call insert(zoltan_global_new_elements, buf(head + j))
          if(zoltan_global_preserve_mesh_regions) then
             call insert(zoltan_global_universal_element_number_to_region_id, buf(head + j), buf(head + j + sz))
          end if
       end do
       if(zoltan_global_preserve_mesh_regions) then
          head = head + 2*sz + 1
       else
          head = head + sz + 1
       end if
       
       ! and record who owns this in the halo
       call insert(zoltan_global_receives(new_owner+1), global_ids(i))
       
       ! and record the snelist information
       sz = buf(head)
       do j=1,sz
          call insert(zoltan_global_new_snelist(new_local_number), buf(head + j))
          call insert(zoltan_global_universal_surface_number_to_surface_id, buf(head + j), buf(head + j + sz))
          call insert(zoltan_global_universal_surface_number_to_element_owner, buf(head + j), buf(head + j + 2*sz))
          call insert(zoltan_global_new_surface_elements, buf(head + j))
       end do
       head = head + 3*sz + 1
    end do
    
    ierr = ZOLTAN_OK
  end subroutine zoltan_cb_unpack_halo_nodes

#endif
  
end module zoltan_callbacks
