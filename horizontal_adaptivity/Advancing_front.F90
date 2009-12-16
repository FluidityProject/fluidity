#include "fdebug.h"

module hadapt_advancing_front
  use fields
  use sparse_tools
  use quicksort
  use linked_lists
  use adjacency_lists
  use meshdiagnostics
  implicit none

  contains

  subroutine generate_layered_mesh(mesh, h_mesh)
    !! Given a columnar mesh with the positions of the vertical
    !! nodes, fill in the elements. 
    type(vector_field), intent(inout) :: mesh
    type(vector_field), intent(inout) :: h_mesh

    type(csr_sparsity) :: columns
    type(ilist), dimension(node_count(h_mesh)) :: chains

    ! node_count(mesh) - node_count(h_mesh) == number of nodes hanging in chains off horizontal mesh
    ! the heights act as priorities in a priority queue, while the indices
    ! record on which chain this height lives.
    real, dimension(node_count(mesh) - node_count(h_mesh)) :: heights
    integer, dimension(node_count(mesh) - node_count(h_mesh)) :: indices
    integer, dimension(node_count(mesh) - node_count(h_mesh)) :: sorted

 
    integer :: h_node ! horizontal node
    integer :: c_node ! chain_node
    integer :: chain_len
    integer, dimension(:), pointer :: chain_ptr, ndglno_ptr, h_elements, h_ndglno
    integer, dimension(mesh_dim(mesh) - 1) :: other_chain_heads
    integer :: i, j, k, l
    integer :: chain

    integer :: dim
    integer :: ele, h_ele

    real :: vol
    
    type(csr_sparsity), pointer :: nelist

    dim = mesh_dim(mesh)

    ! Step 0. Generate the data structures necessary:
    ! 0.1 The chain-element connectivity list
    nelist => extract_nelist(h_mesh)
    
    ! 0.1 From the node to column map, create a column to node map
    call create_columns_sparsity(columns, mesh%mesh)

    ! 0.2 The linked lists representing the chains.
    do h_node=1,node_count(h_mesh)
      chain_len = row_length(columns, h_node)
      chain_ptr => row_m_ptr(columns, h_node)
      do c_node=1,chain_len
        call insert(chains(h_node), chain_ptr(c_node))
      end do
    end do

    ! 0.3 The heights and indices of the hanging nodes.
    i = 1
    do h_node=1,node_count(h_mesh)
      chain_len = row_length(columns, h_node)
      chain_ptr => row_m_ptr(columns, h_node)
      ! Note: c_node loops from 2 to skip the node in the chain on the horizontal mesh
      do c_node=2,chain_len
        heights(i) = abs(node_val(mesh, dim, chain_ptr(c_node)))
        indices(i) = h_node
        i = i + 1
      end do
    end do
      
    ! bored with it now
    call deallocate(columns)

    ! Sort the chains by height
    call qsort(heights, sorted)

    ! The main loop.

    ele = 1
    ndglno_ptr => ele_nodes(mesh, ele)
    do i=1,size(sorted)
      ! Get the chain we're dealing with.
      chain = indices(sorted(i))
      assert(associated(chains(chain)%firstnode%next))

      ! So we're going to form an element.
      ! It's going to have nodes
      ! [chain_head, next_along_chain_from_chain_head, [others]]
      ! where others are the chain_heads of chains in elements
      ! neighbouring the current chain.

      h_elements => row_m_ptr(nelist, chain)
      do j=1,size(h_elements)
        h_ele = h_elements(j)
        ! Get the chains in sele that are NOT the current chain.
        h_ndglno => ele_nodes(h_mesh, h_ele)
        l = 1
        do k=1,size(h_ndglno)
          if (h_ndglno(k) /= chain) then
            other_chain_heads(l) = chains(h_ndglno(k))%firstnode%value
            l = l + 1
          end if
        end do

        ! So now we know the heads of the chains next to us.
        ! Let's form the element! Quick! quick!

        ndglno_ptr(1) = chains(chain)%firstnode%value
        ndglno_ptr(2) = chains(chain)%firstnode%next%value
        ndglno_ptr(3:dim+1) = other_chain_heads

        ! Now we have to orient the element. 

        vol = simplex_volume(mesh, ele)
        assert(abs(vol) /= 0.0)
        if (vol < 0.0) then
          l = ndglno_ptr(1)
          ndglno_ptr(1) = ndglno_ptr(2)
          ndglno_ptr(2) = l
        end if

        ! Get ready to process the next element

        ele = ele + 1
        if (ele <= ele_count(mesh)) then
          ndglno_ptr => ele_nodes(mesh, ele)
        end if

      end do

      ! And advance down the chain ..
      k = pop(chains(chain))
    end do

    do h_node=1,node_count(h_mesh)
      call deallocate(chains(h_node))
    end do

  end subroutine generate_layered_mesh
    
  subroutine create_columns_sparsity(columns, mesh)
    !! Auxillary routine that creates a sparsity of which the rows
    !! are the columns in a mesh, i.e. column indices of the sparsity correspond
    !! to nodes in a mesh column. This is created from the node to column map mesh%columns
    type(csr_sparsity), intent(out):: columns
    type(mesh_type), intent(in):: mesh
    
    type(csr_sparsity):: node2column_sparsity
    integer:: i, no_nodes, no_columns
    
    if (.not. associated(mesh%columns)) then
      FLAbort("Called create_columns_sparsity on a mesh without columns")
    end if
    
    no_columns=maxval(mesh%columns)
    no_nodes=node_count(mesh)
    
    ! first create trivial node to column sparsity
    call allocate(node2column_sparsity, no_nodes, no_columns, no_nodes, &
      diag=.false., name="Node2ColumnSparsity")
    
    ! each row (corresp. to a node) only has one entry
    do i=1, no_nodes+1
      node2column_sparsity%findrm(i)=i
    end do
      
    node2column_sparsity%colm=mesh%columns
    
    ! now "columns" is the transpose of that:    
    columns=transpose(node2column_sparsity)
    columns%name=trim(mesh%name)//"ColumnsSparsity"
    
    call deallocate(node2column_sparsity)
    
  end subroutine create_columns_sparsity
    
end module hadapt_advancing_front
