#include "fdebug.h"

module hadapt_advancing_front
  use fields
  use sparse_tools
  use quicksort
  use linked_lists
  use adjacency_lists
  use meshdiagnostics
  use data_structures
  use vtk_interfaces
  use spud
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
    integer :: ele, h_ele, snloc
    ! the maximum amount of faces you could possibly want to add is
    ! number of elements in the extruded mesh (ele_count(mesh))
    ! x number of faces per element (in the absence of face_count, use ele_loc)
    integer, dimension(ele_count(mesh) * ele_loc(mesh, 1)) :: element_owners, boundary_ids
    integer, dimension(ele_count(mesh) * ele_loc(mesh, 1) * mesh_dim(mesh)), target :: sndgln
    integer :: faces_seen
    integer :: bottom_count, top_count

    real :: vol
    type(integer_set) :: bottom_nodes, top_nodes
    type(csr_sparsity), pointer :: nelist
    integer, dimension(:), pointer :: nodes, faces
    integer :: top_surface_id, bottom_surface_id
    integer, dimension(1) :: other_node

    dim = mesh_dim(mesh)

    ! Step 0. Generate the data structures necessary:
    ! 0.1 The chain-element connectivity list
    nelist => extract_nelist(h_mesh)
    
    ! 0.1 From the node to column map, create a column to node map
    call create_columns_sparsity(columns, mesh%mesh)

    call allocate(bottom_nodes)
    call allocate(top_nodes)

    ! 0.2 The linked lists representing the chains.
    do h_node=1,node_count(h_mesh)
      chain_len = row_length(columns, h_node)
      chain_ptr => row_m_ptr(columns, h_node)
      do c_node=1,chain_len
        call insert(chains(h_node), chain_ptr(c_node))
      end do
      call insert(top_nodes, chain_ptr(1))
      call insert(bottom_nodes, chain_ptr(chain_len))
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

    faces_seen = 0
    if (has_faces(h_mesh%mesh)) then
      snloc = mesh_dim(mesh)
      assert( snloc==face_loc(h_mesh,1)+1 )
    end if

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

        ! if the horizontal element has any surface faces, add them to the extruded surface mesh
        if (has_faces(h_mesh%mesh)) then
          faces => ele_faces(h_mesh, h_ele)
          do k=1, size(faces)
            if (faces(k)<=surface_element_count(h_mesh)) then
              if (.not. any(chain == face_global_nodes(h_mesh, faces(k)))) cycle
              faces_seen = faces_seen + 1
              element_owners(faces_seen) = ele
              boundary_ids(faces_seen) = surface_element_id(h_mesh, faces(k))
              nodes => sndgln((faces_seen-1)*snloc+1:faces_seen*snloc)
              nodes(1) = chains(chain)%firstnode%value
              nodes(2) = chains(chain)%firstnode%next%value

              if (mesh_dim(mesh) == 3) then
                other_node = pack(face_global_nodes(h_mesh, faces(k)), (chain /= face_global_nodes(h_mesh, faces(k))))
                nodes(3) = chains(other_node(1))%firstnode%value
              end if
            end if
          end do
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

    call get_option(trim(mesh%mesh%option_path) //'/from_mesh/extrude/top_surface_id', top_surface_id, default=0)
    call get_option(trim(mesh%mesh%option_path) //'/from_mesh/extrude/bottom_surface_id', bottom_surface_id, default=0)

    if (has_faces(h_mesh%mesh)) then
      ! Add the top faces, and the bottom ones:

      do ele=1,ele_count(mesh)
        nodes => ele_nodes(mesh, ele)
        bottom_count = count(has_value(bottom_nodes, nodes))
        top_count = count(has_value(top_nodes, nodes))

        if (bottom_count == ele_loc(mesh, ele) - 1) then
          faces_seen = faces_seen + 1
          element_owners(faces_seen) = ele
          boundary_ids(faces_seen) = bottom_surface_id
          sndgln((faces_seen-1)*snloc+1:faces_seen*snloc) = pack(nodes, has_value(bottom_nodes, nodes))
        end if

        if (top_count == ele_loc(mesh, ele) - 1) then
          faces_seen = faces_seen + 1
          element_owners(faces_seen) = ele
          boundary_ids(faces_seen) = top_surface_id
          sndgln((faces_seen-1)*snloc+1:faces_seen*snloc) = pack(nodes, has_value(top_nodes, nodes))
        end if
      end do

      call add_faces(mesh%mesh, sndgln=sndgln(1:faces_seen*snloc), element_owner=element_owners(1:faces_seen), boundary_ids=boundary_ids(1:faces_seen))
    end if

    call deallocate(top_nodes)
    call deallocate(bottom_nodes)


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
