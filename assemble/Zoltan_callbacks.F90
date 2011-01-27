#include "fdebug.h"

module zoltan_callbacks

  use zoltan
  use zoltan_global_variables

  ! Needed for zoltan_cb_owned_node_count
  use halos, only: halo_nowned_nodes, halo_node_owners, get_owned_nodes, halo_universal_number
  use metric_tools

  ! Needed for zoltan_cb_get_owned_nodes
  ! - added get_owned_nodes, halo_universal_number to use halos
  use sparse_tools, only: row_length

  ! Needed for zoltan_cb_get_num_edges
  use global_parameters, only: real_size, OPTION_PATH_LEN
  use parallel_tools, only: getrank

  ! Needed for zoltan_cb_get_edge_list
  ! - added halo_node_owners to use halos
  use mpi_interfaces

  ! Needed for zoltan_cb_pack_node_sizes
  ! - added real_size to use global_parameters
  use data_structures, only: key_count

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

end module zoltan_callbacks
