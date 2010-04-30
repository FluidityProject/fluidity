#include "fdebug.h"

module hadapt_advancing_front
  use fields
  use sparse_tools
  use quicksort
  use linked_lists
  use adjacency_lists
  use meshdiagnostics
  use data_structures
  use spud
  use halos
  use halos_derivation
  implicit none

  contains

  subroutine generate_layered_mesh(mesh, h_mesh)
    !! Given a columnar mesh with the positions of the vertical
    !! nodes, fill in the elements. 
    type(vector_field), intent(inout) :: mesh
    type(vector_field), intent(inout) :: h_mesh

    type(mesh_type) :: in_mesh
    type(csr_sparsity) :: columns
    type(ilist), dimension(node_count(h_mesh)) :: chains
    type(integer_hash_table):: old2new_ele

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

    integer, dimension(:), allocatable:: unn, unn_hanging_nodes, halo_level
    real :: vol
    type(csr_sparsity), pointer :: nelist
    integer, dimension(:), pointer :: nodes, faces, neigh
    integer, dimension(1) :: other_node
    logical:: adjacent_to_owned_element, adjacent_to_owned_column, shared_face
    
    ! stuff to do with region ids
    logical, dimension(node_count(h_mesh)) :: chain_at_top
    integer :: top_surface_id, bottom_surface_id
    integer, dimension(ele_count(h_mesh)) :: bottom_surface_ids, top_surface_ids
    logical :: top_element, bottom_element
    integer, dimension(:), allocatable :: region_ids
    integer :: n_regions, r
    logical :: apply_region_ids
    integer, dimension(2) :: shape_option
    
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
      
    ! Sort the chains by height
    call qsort(heights, sorted)
    
    if (associated(h_mesh%mesh%halos)) then
      allocate(unn(1:node_count(mesh)), unn_hanging_nodes(size(heights)))
      call get_universal_numbering(mesh%mesh%halos(2), unn)
      i = 1
      do h_node=1,node_count(h_mesh)
        chain_len = row_length(columns, h_node)
        chain_ptr => row_m_ptr(columns, h_node)
        ! Note: c_node loops from 2 to skip the node in the chain on the horizontal mesh
        do c_node=2,chain_len
          unn_hanging_nodes(i) = unn(chain_ptr(c_node))
          i = i + 1
        end do
      end do
      assert( i==size(heights)+1 )
      call enforce_consistent_parallel_ordering(heights, sorted, unn_hanging_nodes)
      deallocate(unn, unn_hanging_nodes)
      
      assert(has_faces(h_mesh%mesh)) ! needed for halo1 element recognition
      allocate(halo_level(1:element_count(mesh)))
    end if

    ! bored with it now
    call deallocate(columns)

    ! get the region id information (if any) plus the top and bottom surface ids

    n_regions = option_count(trim(mesh%mesh%option_path)//'/from_mesh/extrude/regions')
    apply_region_ids = (n_regions>1)

    do r = 0, n_regions-1
    
      if(apply_region_ids) then
        shape_option=option_shape(trim(mesh%mesh%option_path)//&
                                  "/from_mesh/extrude/regions["//int2str(r)//&
                                  "]/region_ids")
        allocate(region_ids(1:shape_option(1)))
        call get_option(trim(mesh%mesh%option_path)//&
                        "/from_mesh/extrude/regions["//int2str(r)//&
                        "]/region_ids", region_ids)
      end if
      
      call get_option(trim(mesh%mesh%option_path)// &
                      '/from_mesh/extrude/regions['//int2str(r)//&
                      ']/top_surface_id', top_surface_id, default=0)
      call get_option(trim(mesh%mesh%option_path)// &
                      '/from_mesh/extrude/regions['//int2str(r)//&
                      ']/bottom_surface_id', bottom_surface_id, default=0)
                      
      do h_ele = 1, size(top_surface_ids)
        if(apply_region_ids) then
          if(.not. any(h_mesh%mesh%region_ids(h_ele)==region_ids)) cycle
        end if
        
        top_surface_ids(h_ele) = top_surface_id
        bottom_surface_ids(h_ele) = bottom_surface_id
      end do
      
      if(apply_region_ids) deallocate(region_ids)
      
    end do

    ! The main loop.

    faces_seen = 0
    if (has_faces(h_mesh%mesh)) then
      snloc = mesh_dim(mesh)
      assert( snloc==face_loc(h_mesh,1)+1 )
    end if

    chain_at_top=.true.
    
    ele = 0
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
        top_element = all(chain_at_top(h_ndglno))
        bottom_element = chains(chain)%length==2
        do k=1,size(h_ndglno)
          if (h_ndglno(k) /= chain) then
            other_chain_heads(l) = chains(h_ndglno(k))%firstnode%value
            bottom_element = bottom_element .and. chains(h_ndglno(k))%length==1
            l = l + 1
          end if
        end do

        ! So now we know the heads of the chains next to us.
        ! Let's form the element! Quick! quick!
        ele = ele + 1
        ndglno_ptr => ele_nodes(mesh, ele)

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
          ! first the top and bottom faces
          if (top_element) then
            faces_seen = faces_seen + 1
            element_owners(faces_seen) = ele
            boundary_ids(faces_seen) = top_surface_ids(h_ele)
            nodes => sndgln((faces_seen-1)*snloc+1:faces_seen*snloc)
            nodes(1) = chains(chain)%firstnode%value
            nodes(2:dim) = other_chain_heads
          end if
          if (bottom_element) then
            faces_seen = faces_seen + 1
            element_owners(faces_seen) = ele
            boundary_ids(faces_seen) = bottom_surface_ids(h_ele)
            nodes => sndgln((faces_seen-1)*snloc+1:faces_seen*snloc)
            nodes(1) = chains(chain)%firstnode%next%value
            nodes(2:dim) = other_chain_heads
          end if

          faces => ele_faces(h_mesh, h_ele)
          neigh => ele_neigh(h_mesh, h_ele)
          adjacent_to_owned_element=.false.
          adjacent_to_owned_column=.false.
          do k=1, size(faces)
            ! whether this horizontal face is above a face of our new element
            shared_face = any(chain == face_global_nodes(h_mesh, faces(k)))
            if (shared_face .and. faces(k)<=surface_element_count(h_mesh)) then
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
            ! grab the opportunity to see whether this is a halo1 or halo2 element
            if (neigh(k)>0) then
              if (element_owned(h_mesh%mesh, neigh(k))) then
                adjacent_to_owned_column = .true.
                if (shared_face) then
                  adjacent_to_owned_element = .true.
                end if
              end if
            end if
          end do
        end if
        
        if (associated(h_mesh%mesh%halos)) then
          if(element_owned(h_mesh%mesh, h_ele)) then
            ! element ownership (based on the process with lowest rank 
            ! owning any node in the element), nicely transfers to the columns
            halo_level(ele)=0
          else if (adjacent_to_owned_element) then
            ! this is not true for halo1 - only those that directly face owned elements
            halo_level(ele)=1
          else if (adjacent_to_owned_column) then
            ! extruded halo2 elements are necessarily in a column under 
            ! a halo1 element
            halo_level(ele)=2
          else
            ! misc elements - not owned, not in any receive lists
            halo_level(ele)=3
          end if
        end if

      end do

      ! And advance down the chain ..
      k = pop(chains(chain))
      chain_at_top(chain)=.false.
    end do
      
    assert(ele==element_count(mesh))
    
    if (associated(h_mesh%mesh%halos)) then
      ! now reorder the elements according to halo level
    
      ! preserve %ndglno on in_mesh
      in_mesh = mesh%mesh
      ! get a new one for mesh
      allocate(mesh%mesh%ndglno(size(in_mesh%ndglno)))
      call allocate( old2new_ele )
      
      ! first the owned elements
      ele = 0 ! new element nr. in mesh
      do i=0, 3
        ! loop over old element nrs j:
        do j=1, element_count(in_mesh)
          if (halo_level(j)==i) then
            ele = ele + 1 ! new nr.
            call set_ele_nodes(mesh%mesh, ele, ele_nodes(in_mesh, j))
            call insert(old2new_ele, j, ele)
          end if
        end do
      end do
        
      deallocate(in_mesh%ndglno)
      deallocate(halo_level)
    
      ! renumber element ownership of faces
      do i=1, faces_seen
        element_owners(i) = fetch(old2new_ele, element_owners(i))
      end do
      call deallocate(old2new_ele)
      
      call derive_other_extruded_halos(h_mesh%mesh, mesh%mesh)
    end if

    if (has_faces(h_mesh%mesh)) then
      ! Add the top faces, and the bottom ones:
      call add_faces(mesh%mesh, sndgln=sndgln(1:faces_seen*snloc), element_owner=element_owners(1:faces_seen), boundary_ids=boundary_ids(1:faces_seen))
    end if
    
    if (associated(h_mesh%mesh%halos)) then
      ! make sure we obey zoltan's ordering convention
      call reorder_element_numbering(mesh)
    end if

    do h_node=1,node_count(h_mesh)
      call deallocate(chains(h_node))
    end do

  end subroutine generate_layered_mesh

  subroutine generate_radially_layered_mesh(mesh, shell_mesh)
    !! Given a radially columnar mesh with the positions
    !! nodes, fill in the elements. 
    type(vector_field), intent(inout) :: mesh
    type(vector_field), intent(inout) :: shell_mesh

    type(mesh_type) :: in_mesh
    type(csr_sparsity) :: columns
    type(ilist), dimension(node_count(shell_mesh)) :: chains
    type(integer_hash_table):: old2new_ele

    ! node_count(mesh) - node_count(shell_mesh) == number of nodes below the outer shell.
    ! The radial coordinate act as priorities in a priority queue, while the indices
    ! record on which chain this radius lives.
    real, dimension(node_count(mesh) - node_count(shell_mesh)) :: radius
    integer, dimension(node_count(mesh) - node_count(shell_mesh)) :: indices
    integer, dimension(node_count(mesh) - node_count(shell_mesh)) :: sorted

 
    integer :: shell_node ! node on a shell of constant radius
    integer :: c_node ! chain_node
    integer :: chain_len
    integer, dimension(:), pointer :: chain_ptr, ndglno_ptr, shell_elements, shell_ndglno
    integer, dimension(mesh_dim(mesh) - 1) :: other_chain_heads
    integer :: i, j, k, l
    integer :: chain

    integer :: dim
    integer :: ele, shell_ele, snloc
    ! the maximum amount of faces you could possibly want to add is
    ! number of elements in the extruded mesh (ele_count(mesh))
    ! x number of faces per element (in the absence of face_count, use ele_loc)
    integer, dimension(ele_count(mesh) * ele_loc(mesh, 1)) :: element_owners, boundary_ids
    integer, dimension(ele_count(mesh) * ele_loc(mesh, 1) * mesh_dim(mesh)), target :: sndgln
    integer :: faces_seen

    integer, dimension(:), allocatable:: unn, unn_hanging_nodes, halo_level
    real :: vol
    type(csr_sparsity), pointer :: nelist
    integer, dimension(:), pointer :: nodes, faces, neigh
    integer, dimension(1) :: other_node
    logical:: adjacent_to_owned_element, adjacent_to_owned_column, shared_face

    ! stuff to do with region ids
    logical, dimension(node_count(shell_mesh)) :: chain_at_top
    integer :: top_surface_id, bottom_surface_id
    integer, dimension(ele_count(shell_mesh)) :: bottom_surface_ids, top_surface_ids
    logical :: top_element, bottom_element
    integer, dimension(:), allocatable :: region_ids
    integer :: n_regions, r
    logical :: apply_region_ids
    integer, dimension(2) :: shape_option

    dim = mesh_dim(mesh)

    ! Step 0. Generate the data structures necessary:
    ! 0.1 The chain-element connectivity list
    nelist => extract_nelist(shell_mesh)
    
    ! 0.1 From the node to column map, create a column to node map
    call create_columns_sparsity(columns, mesh%mesh)

    ! 0.2 The linked lists representing the chains.
    do shell_node=1,node_count(shell_mesh)
      chain_len = row_length(columns, shell_node)
      chain_ptr => row_m_ptr(columns, shell_node)
      do c_node=1,chain_len
        call insert(chains(shell_node), chain_ptr(c_node))
      end do
    end do

    ! 0.3 The radii and indices of the hanging nodes.
    i = 1
    do shell_node=1,node_count(shell_mesh)
      chain_len = row_length(columns, shell_node)
      chain_ptr => row_m_ptr(columns, shell_node)
      ! Note: c_node loops from 2 to skip the node in the chain on the shell mesh
      do c_node=2,chain_len
        radius(i) = norm2(node_val(mesh, chain_ptr(c_node)))
        indices(i) = shell_node
        i = i + 1
      end do
    end do

    ! Sort the chains by height
    call qsort(radius, sorted)
    
    if (associated(shell_mesh%mesh%halos)) then
      allocate(unn(1:node_count(mesh)), unn_hanging_nodes(size(radius)))
      call get_universal_numbering(mesh%mesh%halos(2), unn)
      i = 1
      do shell_node=1,node_count(shell_mesh)
        chain_len = row_length(columns, shell_node)
        chain_ptr => row_m_ptr(columns, shell_node)
        ! Note: c_node loops from 2 to skip the node in the chain on the horizontal mesh
        do c_node=2,chain_len
          unn_hanging_nodes(i) = unn(chain_ptr(c_node))
          i = i + 1
        end do
      end do
      assert( i==size(radius)+1 )
      call enforce_consistent_parallel_ordering(radius, sorted, unn_hanging_nodes)
      deallocate(unn, unn_hanging_nodes)
      
      assert(has_faces(shell_mesh%mesh)) ! needed for halo1 element recognition
      allocate(halo_level(1:element_count(mesh)))
    end if
      
    ! bored with it now
    call deallocate(columns)

    ! get the region id information (if any) plus the top and bottom surface ids

    n_regions = option_count(trim(mesh%mesh%option_path)//'/from_mesh/extrude/regions')
    apply_region_ids = (n_regions>1)

    do r = 0, n_regions-1
    
      if(apply_region_ids) then
        shape_option=option_shape(trim(mesh%mesh%option_path)//&
                                  "/from_mesh/extrude/regions["//int2str(r)//&
                                  "]/region_ids")
        allocate(region_ids(1:shape_option(1)))
        call get_option(trim(mesh%mesh%option_path)//&
                        "/from_mesh/extrude/regions["//int2str(r)//&
                        "]/region_ids", region_ids)
      end if
      
      call get_option(trim(mesh%mesh%option_path)// &
                      '/from_mesh/extrude/regions['//int2str(r)//&
                      ']/top_surface_id', top_surface_id, default=0)
      call get_option(trim(mesh%mesh%option_path)// &
                      '/from_mesh/extrude/regions['//int2str(r)//&
                      ']/bottom_surface_id', bottom_surface_id, default=0)
                      
      do shell_ele = 1, size(top_surface_ids)
        if(apply_region_ids) then
          if(.not. any(shell_mesh%mesh%region_ids(shell_ele)==region_ids)) cycle
        end if
        
        top_surface_ids(shell_ele) = top_surface_id
        bottom_surface_ids(shell_ele) = bottom_surface_id
      end do
      
      if(apply_region_ids) deallocate(region_ids)
      
    end do

    ! The main loop.

    faces_seen = 0
    if (has_faces(shell_mesh%mesh)) then
      snloc = mesh_dim(mesh)
      assert( snloc==face_loc(shell_mesh,1)+1 )
    end if

    chain_at_top=.true.
    
    ele = 0
!     ndglno_ptr => ele_nodes(mesh, ele)
    do i=1,size(sorted)
      ! Get the chain we're dealing with.
      chain = indices(sorted(i))
      assert(associated(chains(chain)%firstnode%next))

      ! So we're going to form an element.
      ! It's going to have nodes
      ! [chain_head, next_along_chain_from_chain_head, [others]]
      ! where others are the chain_heads of chains in elements
      ! neighbouring the current chain.

      shell_elements => row_m_ptr(nelist, chain)
      do j=1,size(shell_elements)
        shell_ele = shell_elements(j)
        ! Get the chains in sele that are NOT the current chain.
        shell_ndglno => ele_nodes(shell_mesh, shell_ele)
        l = 1
        top_element = all(chain_at_top(shell_ndglno))
        bottom_element = chains(chain)%length==2
        do k=1,size(shell_ndglno)
          if (shell_ndglno(k) /= chain) then
            other_chain_heads(l) = chains(shell_ndglno(k))%firstnode%value
            bottom_element = bottom_element .and. chains(shell_ndglno(k))%length==1
            l = l + 1
          end if
        end do

        ! So now we know the heads of the chains next to us.
        ! Let's form the element! Quick! quick!
        ele = ele + 1
        ndglno_ptr => ele_nodes(mesh, ele)

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
        if (has_faces(shell_mesh%mesh)) then
          ! first the top and bottom faces
          if (top_element) then
            faces_seen = faces_seen + 1
            element_owners(faces_seen) = ele
            boundary_ids(faces_seen) = top_surface_ids(shell_ele)
            nodes => sndgln((faces_seen-1)*snloc+1:faces_seen*snloc)
            nodes(1) = chains(chain)%firstnode%value
            nodes(2:dim) = other_chain_heads
          end if
          if (bottom_element) then
            faces_seen = faces_seen + 1
            element_owners(faces_seen) = ele
            boundary_ids(faces_seen) = bottom_surface_ids(shell_ele)
            nodes => sndgln((faces_seen-1)*snloc+1:faces_seen*snloc)
            nodes(1) = chains(chain)%firstnode%next%value
            nodes(2:dim) = other_chain_heads
          end if

          faces => ele_faces(shell_mesh, shell_ele)
          neigh => ele_neigh(shell_mesh, shell_ele)
          adjacent_to_owned_element=.false.
          adjacent_to_owned_column=.false.
          do k=1, size(faces)
            ! whether this horizontal face is above a face of our new element
            shared_face = any(chain == face_global_nodes(shell_mesh, faces(k)))
            if (shared_face .and. faces(k)<=surface_element_count(shell_mesh)) then
              faces_seen = faces_seen + 1
              element_owners(faces_seen) = ele
              boundary_ids(faces_seen) = surface_element_id(shell_mesh, faces(k))
              nodes => sndgln((faces_seen-1)*snloc+1:faces_seen*snloc)
              nodes(1) = chains(chain)%firstnode%value
              nodes(2) = chains(chain)%firstnode%next%value

              if (mesh_dim(mesh) == 3) then
                other_node = pack(face_global_nodes(shell_mesh, faces(k)), (chain /= face_global_nodes(shell_mesh, faces(k))))
                nodes(3) = chains(other_node(1))%firstnode%value
              end if
            end if
            ! grab the opportunity to see whether this is a halo1 or halo2 element
            if (neigh(k)>0) then
              if (element_owned(shell_mesh%mesh, neigh(k))) then
                adjacent_to_owned_column = .true.
                if (shared_face) then
                  adjacent_to_owned_element = .true.
                end if
              end if
            end if
          end do
        end if
        
        if (associated(shell_mesh%mesh%halos)) then
          if(element_owned(shell_mesh%mesh, shell_ele)) then
            ! element ownership (based on the process with lowest rank 
            ! owning any node in the element), nicely transfers to the columns
            halo_level(ele)=0
          else if (adjacent_to_owned_element) then
            ! this is not true for halo1 - only those that directly face owned elements
            halo_level(ele)=1
          else if (adjacent_to_owned_column) then
            ! extruded halo2 elements are necessarily in a column under 
            ! a halo1 element
            halo_level(ele)=2
          else
            ! misc elements - not owned, not in any receive lists
            halo_level(ele)=3
          end if
        end if

      end do

      ! And advance down the chain ..
      k = pop(chains(chain))
      chain_at_top(chain)=.false.
    end do

    assert(ele==element_count(mesh))
    
    if (associated(shell_mesh%mesh%halos)) then
      ! now reorder the elements according to halo level
    
      ! preserve %ndglno on in_mesh
      in_mesh = mesh%mesh
      ! get a new one for mesh
      allocate(mesh%mesh%ndglno(size(in_mesh%ndglno)))
      call allocate( old2new_ele )
      
      ! first the owned elements
      ele = 0 ! new element nr. in mesh
      do i=0, 3
        ! loop over old element nrs j:
        do j=1, element_count(in_mesh)
          if (halo_level(j)==i) then
            ele = ele + 1 ! new nr.
            call set_ele_nodes(mesh%mesh, ele, ele_nodes(in_mesh, j))
            call insert(old2new_ele, j, ele)
          end if
        end do
      end do
        
      deallocate(in_mesh%ndglno)
      deallocate(halo_level)
    
      ! renumber element ownership of faces
      do i=1, faces_seen
        element_owners(i) = fetch(old2new_ele, element_owners(i))
      end do
      call deallocate(old2new_ele)
      
      call derive_other_extruded_halos(shell_mesh%mesh, mesh%mesh)
    end if

    if (has_faces(shell_mesh%mesh)) then
      call add_faces(mesh%mesh, sndgln=sndgln(1:faces_seen*snloc), element_owner=element_owners(1:faces_seen), boundary_ids=boundary_ids(1:faces_seen))
    end if

    if (associated(shell_mesh%mesh%halos)) then
      ! make sure we obey zoltan's ordering convention
      call reorder_element_numbering(mesh)
    end if

    do shell_node=1,node_count(shell_mesh)
      call deallocate(chains(shell_node))
    end do

  end subroutine generate_radially_layered_mesh
    
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
    
  subroutine enforce_consistent_parallel_ordering(values, ordering, unn)
    ! for a sorted array values(ordering) on the local nodes of a mesh
    ! (owned and halos nodes) make sure the ordering is consistent 
    ! between processes, i.e. if a pair of nodes is seen by more processes, 
    ! on each of those processes the pair comes in the same order in 
    ! their local values(ordering). For values that are "close" the ordering
    ! is done via the universal node numbering
    ! unsorted values
    real, dimension(:), intent(in):: values
    ! values(ordering) should be ordered locally already
    integer, dimension(:), intent(inout):: ordering
    ! unn(ordering) gives the universal node numbers
    integer, dimension(:), intent(in):: unn
    
    real:: eps, value1, value2, diff
    integer:: i, i1
    
    ! first we need to come up with a definition of what is
    ! "sufficiently close" such that we can consider 2 values equal
    eps = find_separating_eps(values, ordering)
    
    value1 = values(ordering(1))
    i1 = 1
    do i=2, size(ordering)
      value2 = values(ordering(i))
      diff = abs(value2-value1)
      if (diff>eps) then
        ! end of sequence of close values
        if (i1<i-1) then
          ! reorder entries i1 to i-1 according to unn
          call order_subsequence_on_unn(i1, i-1)
        end if
        ! start new sequence
        value1 = value2
        i1 = i
      end if
    end do
    if (i1<i-1) then
      ! reorder entries i1 to i-1 according to unn
      call order_subsequence_on_unn(i1, i-1)
    end if
      
    contains
    
    subroutine order_subsequence_on_unn(i1, i2)
      integer, intent(in):: i1, i2
      
      integer, dimension(i2-i1+1):: new_order
      
      call qsort(unn(ordering(i1:i2)), new_order)
      ! the returned new_order, starts from 1 but should start from i1
      new_order = new_order+i1-1
      ! now unn(ordering(new_order)) is sorted
      ! so replace this section of the input ordering
      ordering(i1:i2) = ordering(new_order)
      
    end subroutine order_subsequence_on_unn
    
  end subroutine enforce_consistent_parallel_ordering
    
  function find_separating_eps(values, ordering) result (eps)
    ! find a global epsilon such that for each pair of values in values
    ! either |values(i)-values(j)|<0.9*eps
    ! or |values(i)-values(j)|>1.1*eps
    ! this means we can do local comparisons that are consistent between 2 processes
    real:: eps
    real, dimension(:), intent(in):: values
    ! values(ordering) should be ordered
    integer, dimension(:), intent(inout):: ordering
  
    real:: eps0, eps_old, diff, value1, value2
    integer:: i, it1, it2
    logical:: too_close
    
    ! initial guess:
    eps0 = 100.0*epsilon(maxval(values)-minval(values))
    call allmax(eps0)
    
    ! now if any 2 subsequent values have 0.9*eps<|x_{i+1}-x_i|<1.1*eps
    ! then processes could still disagree about them being equal
    
    ! this is all a bit hideous - please feel free to think of a better
    ! way to do this
    eps = eps0
    eps_old = eps0
    do it1=1, 1000
      ! check whether this eps 
      do it2=1, 1000
        too_close = .false.
        value1 = values(ordering(1))
        do i=2, size(ordering)
          value2 = values(ordering(i))
          diff = abs(value2-value1)
          if (diff>eps) then
            ! different values with a safe distance, move on
            value1 = value2
          else if (diff>0.9*eps) then
            too_close = .true.
            exit
          end if
        end do
        if (.not. too_close) exit
        ! try a bigger eps
        eps = eps + eps0
      end do
      if (it2>1000) then
        FLAbort("Could not find a suitable epsilon for ordering values")
      end if      
      call allmax(eps)
      ! if eps has not changed we're fine - otherwise we have to
      ! check the new value of eps
      ! (note we're comparing global eps with previously agreed global eps
      !  so this branch should be safe)
      if (eps==eps_old) exit
      eps_old=eps
    end do
    if (it1>1000) then
      FLAbort("Could not find a suitable epsilon for ordering values")
    end if
    ewrite(2,*) "Processes agreed on epsilon to use for ordering:"
    ewrite(2,*) "eps, it1, it2 =", eps, it1, it2
    
  end function find_separating_eps
    
  subroutine derive_other_extruded_halos(h_mesh, out_mesh)
    ! after having derived the 2nd node halo before,
    ! now derive 1st nodal halo and the element halos
    type(mesh_type), intent(in):: h_mesh
    type(mesh_type), intent(inout):: out_mesh
    
    call derive_l1_from_l2_halo(out_mesh, &
      ordering_scheme=halo_ordering_scheme(h_mesh%halos(2)))
    assert(halo_valid_for_communication(out_mesh%halos(1)))
    
    allocate(out_mesh%element_halos(2))
    call derive_element_halo_from_node_halo(out_mesh, &
      ordering_scheme=halo_ordering_scheme(h_mesh%halos(2)))
    
    assert(halo_valid_for_communication(out_mesh%element_halos(1)))
    assert(halo_valid_for_communication(out_mesh%element_halos(2)))
    
  end subroutine derive_other_extruded_halos
  
end module hadapt_advancing_front
