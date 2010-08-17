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

    type(csr_sparsity), pointer :: nelist
    type(scalar_field) :: height_field
    type(mesh_type) :: in_mesh
    type(integer_hash_table):: old2new_ele
 
    ! the maximum amount of faces you could possibly want to add is
    ! number of elements in the extruded mesh (ele_count(mesh))
    ! x number of faces per element (in the absence of face_count, use ele_loc)
    integer, dimension(ele_count(mesh) * ele_loc(mesh, 1)) :: element_owners, boundary_ids
    integer, dimension(ele_count(mesh) * ele_loc(mesh, 1) * mesh_dim(mesh)), target :: sndgln
    integer, dimension(ele_count(h_mesh)) :: bottom_surface_ids, top_surface_ids
    integer, dimension(node_count(mesh)) :: sorted
    real, dimension(node_count(mesh)) :: heights
    integer, dimension(node_count(h_mesh)) :: hanging_node, column_size, column_count
    integer, dimension(mesh_dim(mesh) - 1) :: other_column_heads
    
    integer, dimension(:), allocatable :: region_ids
    integer, dimension(:), allocatable:: halo_level
    integer, dimension(:), pointer :: ndglno_ptr, h_elements, h_ndglno
    integer, dimension(:), pointer :: nodes, faces, neigh
    
    logical:: adjacent_to_owned_element, adjacent_to_owned_column, shared_face
    logical :: top_element, bottom_element
    logical :: apply_region_ids, propagate_region_ids
    
    integer, dimension(2) :: shape_option
    integer, dimension(1) :: other_node
    integer :: top_surface_id, bottom_surface_id
    integer :: faces_seen
    integer :: h_node, node, column, dim
    integer :: n_regions, r
    integer :: ele, h_ele, snloc
    integer :: i, j, k, l
    
    real :: vol
    
    dim = mesh_dim(mesh)
    
    nelist => extract_nelist(h_mesh)

    ! Sort the new nodes by height
    if (have_option('/geometry/spherical_earth')) then
      ! sort on radius
      height_field = magnitude(mesh)
    else
      ! sort on last coordinate
      height_field = extract_scalar_field(mesh, mesh%dim)
    end if
    ! we want it ordered top to bottom
    heights = -height_field%val
    if (associated(h_mesh%mesh%halos)) then
      call parallel_consistent_ordering(heights, mesh%mesh%halos(2), sorted)
    else
      call qsort(heights, sorted)
    end if

    if (have_option('/geometry/spherical_earth')) then
      call deallocate( height_field )
    end if
    
    if (associated(h_mesh%mesh%halos)) then
      assert(has_faces(h_mesh%mesh)) ! needed for halo1 element recognition
      allocate(halo_level(1:element_count(mesh)))
    end if

    assert(associated(mesh%mesh%columns))
    assert(associated(mesh%mesh%element_columns))

    ! if the horizontal mesh has region ids associated with it, the elements below
    ! can inherit these ids
    propagate_region_ids = associated(mesh%mesh%region_ids)
#ifdef DDEBUG
    if(propagate_region_ids) then
      assert(associated(h_mesh%mesh%region_ids))
    end if
#endif

    ! get the region id information used for extrusion (if any) plus the top and bottom surface ids
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

    faces_seen = 0
    if (has_faces(h_mesh%mesh)) then
      snloc = dim
      assert( snloc==face_loc(h_mesh,1)+1 )
    end if

    ! for each column, find out how many nodes there are
    column_size = 0
    do node = 1, node_count(mesh)
      column = mesh%mesh%columns(node)
      column_size(column) = column_size(column) + 1
    end do
      
    ! for each column keep track what the last added node is
    hanging_node = 0    
    ! we start with the top nodes
    l = 0
    do i = 1, size(sorted)
      node = sorted(i)
      column = mesh%mesh%columns(node)
      if (hanging_node(column)==0) then
        l = l + 1
        hanging_node(column) = node
      end if
      if (l==size(hanging_node)) exit
    end do
      
    ! for each column keep track how many nodes we've added already,
    ! at the moment that's just the top node in each column:
    column_count = 1
    
    ! The main loop.
    ele = 0
    do i=1, size(sorted)
      
      node = sorted(i)
      
      ! Get the column we're dealing with
      column = mesh%mesh%columns(node)
      
      if (hanging_node(column)==node) then
        if (column_count(column)==1) then
          ! this is simply the top node, we only start adding elements
          ! when we have at least 2 nodes in this column
          cycle
        else
          ! this node has already been added to the column, something
          ! must have gone wrong horribly
          FLAbort("Internal error in mesh extrustion, generate_layered_mesh")
        end if
      end if
      
      ! So we're going to form an element.
      ! It's going to have nodes
      ! [node, hanging_node(column), [others]]
      ! where others are the hanging_nodes of the neighbouring columns

      h_elements => row_m_ptr(nelist, column)
      do j=1,size(h_elements)
        h_ele = h_elements(j)
        ! Get the columns in sele that are NOT the current column.
        h_ndglno => ele_nodes(h_mesh, h_ele)
        l = 0
        ! we're adding a top element if all associated columns are still at top
        top_element = all(column_count(h_ndglno)==1)
        ! we're adding a bottom element if this column is one-but-last
        ! and the other columns are at the last node (checked in the loop below)
        bottom_element = column_count(column)==column_size(column)-1
        do k=1,size(h_ndglno)
          h_node = h_ndglno(k)
          if (h_node /= column) then
            ! add other column
            l = l +1
            other_column_heads(l) = hanging_node(h_node)
            ! as promised check that other columns are at the last node
            bottom_element = bottom_element .and. column_count(h_node)==column_size(h_node)
          end if
        end do

        ! So now we know the heads of the columns next to us.
        ! Let's form the element! Quick! quick!
        ele = ele + 1
        ndglno_ptr => ele_nodes(mesh, ele)

        ndglno_ptr(1) = hanging_node(column)
        ndglno_ptr(2) = node
        ndglno_ptr(3:dim+1) = other_column_heads

        ! Now we have to orient the element. 

        vol = simplex_volume(mesh, ele)
        assert(abs(vol) /= 0.0)
        if (vol < 0.0) then
          l = ndglno_ptr(1)
          ndglno_ptr(1) = ndglno_ptr(2)
          ndglno_ptr(2) = l
        end if
        
        ! if the horizontal mesh has region_ids these may as well be preserved here
        if(propagate_region_ids) then
          mesh%mesh%region_ids(ele) = h_mesh%mesh%region_ids(h_ele)
        end if
        
        ! we now know the relationship between the mesh element and the h_mesh surface element
        ! save this for zoltan (particularly useful if we're parallel using zoltan)...
        mesh%mesh%element_columns(ele) = h_ele

        ! if the horizontal element has any surface faces, add them to the extruded surface mesh
        if (has_faces(h_mesh%mesh)) then
          ! first the top and bottom faces
          if (top_element) then
            faces_seen = faces_seen + 1
            element_owners(faces_seen) = ele
            boundary_ids(faces_seen) = top_surface_ids(h_ele)
            nodes => sndgln((faces_seen-1)*snloc+1:faces_seen*snloc)
            nodes(1) = hanging_node(column)
            nodes(2:dim) = other_column_heads
          end if
          if (bottom_element) then
            faces_seen = faces_seen + 1
            element_owners(faces_seen) = ele
            boundary_ids(faces_seen) = bottom_surface_ids(h_ele)
            nodes => sndgln((faces_seen-1)*snloc+1:faces_seen*snloc)
            nodes(1) = node
            nodes(2:dim) = other_column_heads
          end if

          faces => ele_faces(h_mesh, h_ele)
          neigh => ele_neigh(h_mesh, h_ele)
          adjacent_to_owned_element=.false.
          adjacent_to_owned_column=.false.
          do k=1, size(faces)
            ! whether this horizontal face is above a face of our new element
            shared_face = any(column == face_global_nodes(h_mesh, faces(k)))
            if (shared_face .and. faces(k)<=surface_element_count(h_mesh)) then
              faces_seen = faces_seen + 1
              element_owners(faces_seen) = ele
              boundary_ids(faces_seen) = surface_element_id(h_mesh, faces(k))
              nodes => sndgln((faces_seen-1)*snloc+1:faces_seen*snloc)
              nodes(1) = hanging_node(column)
              nodes(2) = node

              if (dim == 3) then
                other_node = pack(face_global_nodes(h_mesh, faces(k)), mask=face_global_nodes(h_mesh, faces(k)) /= column)
                nodes(3) = hanging_node(other_node(1))
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

      ! And advance down the column
      hanging_node(column) = node
      column_count(column) = column_count(column) + 1
      assert( column_count(column)<=column_size(column) )
    end do
      
    assert(ele==element_count(mesh))
    assert(all(mesh%mesh%element_columns>0))
    
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
      call add_faces(mesh%mesh, sndgln=sndgln(1:faces_seen*snloc), &
        element_owner=element_owners(1:faces_seen), &
        boundary_ids=boundary_ids(1:faces_seen))
    end if
    
    if (associated(h_mesh%mesh%halos)) then
      ! make sure we obey zoltan's ordering convention
      call reorder_element_numbering(mesh)
    end if

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
    
  subroutine parallel_consistent_ordering(values, halo, index)

    real, dimension(:), intent(in):: values
    type(halo_type), intent(in):: halo
    integer, dimension(:), intent(out):: index
    
    real :: minv, maxv, int_base, int_range
    integer, dimension(size(values)):: ints, unn, reindex
    integer :: i, int1, int2, start
    
    minv=minval(values)
    maxv=maxval(values)
    call allmin(minv)
    call allmax(maxv)
    
    ! range of floats that can be rounded to integers
    ! we use -huge/2 to +huge/2, just to be sure
    int_base=-real(huge(1)/2)
    int_range=real(huge(1))
    ! ah what the heck, heterogenous cluster computing is all the rage
    call allmax(int_base)
    call allmin(int_range)
    
    ! round the real values to integer with as much precision as possible
    ints=floor( int_base+(values-minv)/(maxv-minv)*int_range )
    
    call halo_update(halo, ints)
    
    call qsort(ints, index)
    
    call get_universal_numbering(halo, unn)
    
    ! now go through to order equal heights on universal node number
    start=1
    int1=ints(index(1))
    do i=2, size(ints)+1
      if (i<=size(ints)) then
        int2=ints(index(i))
      else
        int2=int1+1
      end if
      if (int2>int1) then
        if (i-start>1) then
          ! sort entries start to i-1 on unn
          call qsort(unn(index(start:i-1)), reindex(1:i-start))
          index(start:i-1)=index(start+reindex(1:i-start)-1)
        end if
        ! start new series
        start=i
        int1=int2
      else if (int2<int1) then
        FLAbort("Something went wrong in sorting")
      end if
    end do
    
  end subroutine parallel_consistent_ordering
    
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
