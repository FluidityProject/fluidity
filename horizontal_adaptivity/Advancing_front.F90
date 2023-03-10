#include "fdebug.h"

module hadapt_advancing_front

   use fldebug
   use global_parameters, only: OPTION_PATH_LEN
   use futils
   use quicksort
   use data_structures
   use spud
   use parallel_tools
   use sparse_tools
   use linked_lists
   use adjacency_lists
   use parallel_fields
   use fields
   use meshdiagnostics
   use halos_derivation
   use halos

   implicit none

   private

   public :: generate_layered_mesh, create_columns_sparsity

contains

   subroutine generate_layered_mesh(mesh, h_mesh, layer_nodes)
      !! Given a columnar mesh with the positions of the vertical
      !! nodes, fill in the elements.
      type(vector_field), intent(inout) :: mesh
      type(vector_field), intent(inout) :: h_mesh
      !! Indicates first owned and recv node in each layer - owned nodes and recv nodes
      !! are (seperately) numbered consecutively within each layer - these arrays have an extra entry for convenience
      type(integer_set), dimension(:), intent(in) :: layer_nodes

      type(csr_sparsity), pointer :: nelist
      type(scalar_field) :: height_field
      type(mesh_type) :: in_mesh
      type(integer_hash_table):: old2new_ele
      character(len=OPTION_PATH_LEN) :: region_option_path, layer_path

      ! the maximum amount of faces you could possibly want to add is
      ! number of elements in the extruded mesh (ele_count(mesh))
      ! x number of faces per element (in the absence of face_count, use ele_loc)
      integer, dimension(:), allocatable :: element_owners, boundary_ids
      integer, dimension(:), target, allocatable :: sndgln
      integer, dimension(:), allocatable :: bottom_surface_ids, top_surface_ids, extruded_region_ids
      integer, dimension(:), allocatable :: sorted, unn, integer_heights
      real, dimension(:), allocatable :: heights
      integer, dimension(:), allocatable :: hanging_node, column_size, column_count
      integer, dimension(mesh_dim(mesh) - 1) :: other_column_heads

      integer, dimension(:), allocatable :: region_ids
      integer, dimension(:), allocatable:: halo_level
      integer, dimension(:), pointer :: ndglno_ptr, h_elements, h_ndglno
      integer, dimension(:), pointer :: facet_nodes, faces, neigh

      logical:: adjacent_to_owned_element, adjacent_to_owned_column, shared_face
      logical :: top_element, bottom_element
      logical :: multiple_regions
      integer, dimension(:), allocatable :: old_element_data

      integer, dimension(2) :: shape_option
      integer, dimension(1) :: other_node
      integer :: top_surface_id, bottom_surface_id, extruded_region_id, extruded_region_id_stat
      integer :: faces_seen
      integer :: h_node, node, column, dim
      integer :: n_regions, r, layer
      integer :: ele, h_ele, snloc
      integer :: i, j, k, l

      logical :: radial_layering

      real :: vol

      ! allocate our arrays
      allocate(element_owners(ele_count(mesh) * ele_loc(mesh, 1)))
      allocate(boundary_ids(ele_count(mesh) * ele_loc(mesh, 1)))
      allocate(sndgln(ele_count(mesh) * ele_loc(mesh, 1) * mesh_dim(mesh)))

      allocate(bottom_surface_ids(ele_count(h_mesh)))
      allocate(top_surface_ids(ele_count(h_mesh)))
      allocate(extruded_region_ids(ele_count(h_mesh)))

      allocate(hanging_node(node_count(h_mesh)))
      allocate(column_size(node_count(h_mesh)))
      allocate(column_count(node_count(h_mesh)))

      dim = mesh_dim(mesh)

      nelist => extract_nelist(h_mesh)

      radial_layering = have_option('/geometry/spherical_earth')
      ! Sort the new nodes by height
      if (radial_layering) then
         ! sort on radius
         height_field = magnitude(mesh)
      else
         ! sort on last coordinate
         height_field = extract_scalar_field(mesh, mesh%dim)
      end if

      if (associated(h_mesh%mesh%halos)) then
         assert(has_faces(h_mesh%mesh)) ! needed for halo1 element recognition
         allocate(halo_level(1:element_count(mesh)))
         allocate(unn(1:node_count(mesh)), integer_heights(1:node_count(mesh)))
         call get_universal_numbering(mesh%mesh%halos(2), unn)
         call create_integer_heights(height_field, mesh%mesh%halos(2), integer_heights)
      end if

      assert(associated(mesh%mesh%columns))
      assert(associated(mesh%mesh%element_columns))

      mesh%mesh%region_ids = 0

      faces_seen = 0
      if (has_faces(h_mesh%mesh)) then
         snloc = dim
         assert( snloc==face_loc(h_mesh,1)+1 )
      end if

      ele = 0

      n_regions = option_count(trim(mesh%mesh%option_path)//'/from_mesh/extrude/regions')
      if (n_regions>0) then
         ! regions directly under extrude/: single layer layer only which is specified directly under regions/
         layer_path = trim(mesh%mesh%option_path)//'/from_mesh/extrude'
         assert(size(layer_nodes)==1)
      else
         layer_path = trim(mesh%mesh%option_path)//'/from_mesh/extrude/layer[0]'
      end if

      ! The main loop.
      layers: do layer = 1, size(layer_nodes)
         ! order nodes within layer
         allocate(sorted(key_count(layer_nodes(layer))))
         if (associated(h_mesh%mesh%halos)) then
            call parallel_consistent_ordering(integer_heights, unn, sorted, layer_nodes(layer))
         else
            allocate(heights(size(sorted)))
            heights = -height_field%val(set2vector(layer_nodes(layer)))
            call qsort(heights, sorted)
            do i=1, size(sorted)
               sorted(i) = fetch(layer_nodes(layer), sorted(i))
            end do
            deallocate(heights)
         end if

         ! column_count and column_size are used to determine top and bottom elements in each column
         ! column_size: for each column, find out how many nodes there are in this layer
         column_size = 0
         do i = 1, size(sorted)
            node = sorted(i)
            column = mesh%mesh%columns(node)
            column_size(column) = column_size(column) + 1
         end do
         ! column_count: how many we've encountered so far - at the moment that's 1 (the top) node per column
         column_count = 1

         if (layer==1) then
            ! for the top layer, find what the starting node is,
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
            assert(l==size(hanging_node))
         else
            ! for layers below we continue with the nodes from the prev. layer
            ! that means we've forgotten to count it in column_size
            column_size = column_size + 1
         end if

         ! get the region id information used for extrusion (if any) plus the top and bottom surface ids
         n_regions = option_count(trim(layer_path) // '/regions')
         multiple_regions = (n_regions>1)

         if (multiple_regions .and. .not. associated(h_mesh%mesh%region_ids)) then
            FLAbort("Multiple extrude regions in options tree but no region ids in mesh")
         end if

         if (layer>1) then
            ! layers below take over the bottom_surface_id from the layer above as their top_surface_id
            top_surface_ids = bottom_surface_ids
         end if


         do r = 0, n_regions-1

            region_option_path = trim(layer_path) // "/regions[" // int2str(r) // "]"
            if(multiple_regions) then
               shape_option=option_shape(trim(region_option_path)// "/region_ids")
               allocate(region_ids(1:shape_option(1)))
               call get_option(trim(region_option_path) // "/region_ids", region_ids)
            end if

            if (layer==1) then
               ! only the top layer specifies a top_surface_id
               call get_option(trim(region_option_path) // '/top_surface_id', top_surface_id, default=0)
            end if
            call get_option(trim(region_option_path) // '/bottom_surface_id', bottom_surface_id, default=0)
            call get_option(trim(region_option_path) // '/extruded_region_id', extruded_region_id, stat=extruded_region_id_stat)

            do h_ele = 1, size(top_surface_ids)
               if(multiple_regions) then
                  if(.not. any(h_mesh%mesh%region_ids(h_ele)==region_ids)) cycle
               end if

               if (layer==1) then
                  top_surface_ids(h_ele) = top_surface_id
               end if
               bottom_surface_ids(h_ele) = bottom_surface_id
               if (extruded_region_id_stat==0) then
                  extruded_region_ids(h_ele) = extruded_region_id
               else if (associated(h_mesh%mesh%region_ids)) then
                  extruded_region_ids(h_ele) = h_mesh%mesh%region_ids(h_ele)
               end if
            end do

            if(multiple_regions) deallocate(region_ids)

         end do

         nodes: do i=1, size(sorted)

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

               mesh%mesh%region_ids(ele) = extruded_region_ids(h_ele)

               ! we now know the relationship between the mesh element and the h_mesh surface element
               ! save this for later...
               mesh%mesh%element_columns(ele) = h_ele

               ! if the horizontal element has any surface faces, add them to the extruded surface mesh
               if (has_faces(h_mesh%mesh)) then
                  ! first the top and bottom faces
                  if (top_element) then
                     faces_seen = faces_seen + 1
                     element_owners(faces_seen) = ele
                     boundary_ids(faces_seen) = top_surface_ids(h_ele)
                     facet_nodes => sndgln((faces_seen-1)*snloc+1:faces_seen*snloc)
                     facet_nodes(1) = hanging_node(column)
                     facet_nodes(2:dim) = other_column_heads
                  end if
                  if (bottom_element) then
                     faces_seen = faces_seen + 1
                     element_owners(faces_seen) = ele
                     boundary_ids(faces_seen) = bottom_surface_ids(h_ele)
                     facet_nodes => sndgln((faces_seen-1)*snloc+1:faces_seen*snloc)
                     facet_nodes(1) = node
                     facet_nodes(2:dim) = other_column_heads
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
                        facet_nodes => sndgln((faces_seen-1)*snloc+1:faces_seen*snloc)
                        facet_nodes(1) = hanging_node(column)
                        facet_nodes(2) = node

                        if (dim == 3) then
                           other_node = pack(face_global_nodes(h_mesh, faces(k)), mask=face_global_nodes(h_mesh, faces(k)) /= column)
                           facet_nodes(3) = hanging_node(other_node(1))
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
         end do nodes

         assert(all(column_count==column_size))

         deallocate(sorted)

         ! if #layers>1, set layer path for next layer
         layer_path = trim(mesh%mesh%option_path) // &
            "/from_mesh/extrude/layer[" // int2str(layer) // ']'


      end do layers

      assert(ele==element_count(mesh))
      assert(all(mesh%mesh%element_columns>0))

      if (radial_layering) then
         call deallocate( height_field )
      end if


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

         ! renumber the element columns
         allocate(old_element_data(size(mesh%mesh%element_columns)))
         old_element_data = mesh%mesh%element_columns
         do i = 1, size(mesh%mesh%element_columns)
            mesh%mesh%element_columns(fetch(old2new_ele, i)) = old_element_data(i)
         end do
         deallocate(old_element_data)

         allocate(old_element_data(size(mesh%mesh%region_ids)))
         old_element_data = mesh%mesh%region_ids
         do i = 1, size(mesh%mesh%region_ids)
            mesh%mesh%region_ids(fetch(old2new_ele, i)) = old_element_data(i)
         end do
         deallocate(old_element_data)

         call deallocate(old2new_ele)

         call derive_other_extruded_halos(h_mesh%mesh, mesh%mesh)
      end if

      if (has_faces(h_mesh%mesh)) then
         if (has_discontinuous_internal_boundaries(h_mesh%mesh)) then
            ! horizontal has element ownership information allowing internal facet pairs
            ! to have seperate surface ids (used in periodic meshes) - this means
            ! the same holds for the extruded mesh
            call add_faces(mesh%mesh, sndgln=sndgln(1:faces_seen*snloc), &
               element_owner=element_owners(1:faces_seen), &
               boundary_ids=boundary_ids(1:faces_seen))
         else
            ! no element ownership is necessary, but in this case only one of each
            ! pair of internal facets should be provided in sndgln unless we tell
            ! add_faces() to filter these out (reordering the facets)
            call add_faces(mesh%mesh, sndgln=sndgln(1:faces_seen*snloc), &
               boundary_ids=boundary_ids(1:faces_seen), &
               allow_duplicate_internal_facets=.true.)
         end if
      end if

      if (associated(h_mesh%mesh%halos)) then
         ! make sure we obey zoltan's ordering convention
         call reorder_element_numbering(mesh)
      end if

      deallocate(element_owners)
      deallocate(boundary_ids)
      deallocate(sndgln)
      deallocate(bottom_surface_ids)
      deallocate(top_surface_ids)
      deallocate(extruded_region_ids)
      deallocate(hanging_node)
      deallocate(column_size)
      deallocate(column_count)


   end subroutine generate_layered_mesh

   subroutine create_columns_sparsity(columns, mesh, positions)
      !! Auxillary routine that creates a sparsity of which the rows
      !! are the columns in a mesh, i.e. column indices of the sparsity correspond
      !! to nodes in a mesh column. This is created from the node to column map mesh%columns
      type(csr_sparsity), intent(out):: columns
      type(mesh_type), intent(in):: mesh
      ! pass in the positions if you want to guarantee that the columns are sorted in descending order
      type(vector_field), intent(in), optional :: positions

      type(csr_sparsity):: node2column_sparsity
      integer:: i, no_nodes, no_columns

      integer, dimension(:), pointer :: column_nodes
      integer, dimension(:), allocatable :: permutation

      if (.not. associated(mesh%columns)) then
         FLAbort("Called create_columns_sparsity on a mesh without columns")
      end if

      no_nodes=node_count(mesh)
      if (no_nodes==0) then
         no_columns = 0
      else
         no_columns=maxval(mesh%columns)
      end if

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

      if(present(positions)) then
         do i = 1, no_columns
            column_nodes => row_m_ptr(columns, i)
            allocate(permutation(size(column_nodes)))
            ! NOTE WELL: here we assume that we only care about the last dimension
            ! i.e. broken on the sphere like everything else!
            call qsort(node_val(positions, positions%dim, column_nodes), permutation)
            call apply_reverse_permutation(column_nodes, permutation)
            deallocate(permutation)
         end do

         columns%sorted_rows = .false.
      end if

   end subroutine create_columns_sparsity

   subroutine create_integer_heights(height_field, halo, integer_heights)
      type(scalar_field), intent(in):: height_field
      type(halo_type), intent(in):: halo
      integer, dimension(:), intent(out):: integer_heights

      real :: minv, maxv, int_base, int_range

      ! NOTE: these are the min and max of -height_field%val
      minv=-maxval(height_field)
      maxv=-minval(height_field)
      call allmin(minv)
      call allmax(maxv)

      ! range of floats that can be rounded to integers
      ! we use -huge/2 to +huge/2, just to be sure
      int_base=-real(huge(1)/2)
      int_range=real(huge(1))
      ! ah what the heck, heterogenous cluster computing is all the rage
      call allmax(int_base)
      call allmin(int_range)

      ! round the real heights to integer with as much precision as possible
      integer_heights=floor( int_base+(-height_field%val-minv)/(maxv-minv)*int_range )

      call halo_update(halo, integer_heights)

   end subroutine create_integer_heights

   subroutine parallel_consistent_ordering(ints, unn, index, layer_nodes)

      ! ints and unn are over the entire local node numbering
      integer, dimension(:), intent(in):: ints ! integer heights
      integer, dimension(:), intent(in):: unn ! for equal integer height, sort on unn secondly
      ! returns a sorted index into the local node numbering, whose length is the number of nodes in the layer
      integer, dimension(:), intent(out):: index
      type(integer_set), intent(in) :: layer_nodes ! nodes in this layer

      integer, dimension(:), allocatable :: ints_packed, unn_packed, reindex
      integer :: i, start, int1, int2

      if (key_count(layer_nodes)==0) return

      allocate(ints_packed(1:key_count(layer_nodes)), unn_packed(1:key_count(layer_nodes)), reindex(1:key_count(layer_nodes)))

      ! ints_packed: integer heights for layer nodes only
      do i=1, key_count(layer_nodes)
         ints_packed(i) = ints(fetch(layer_nodes, i))
      end do

      ! sort it (preliminary)
      call qsort(ints_packed, index)

      ! unn_packed: unn for layer nodes in preliminary order
      do i=1, size(index)
         unn_packed(i) = unn(fetch(layer_nodes, index(i)))
      end do

      ! now go through to order equal heights on universal node number
      start=1
      int1=ints_packed(index(1))
      do i=2, size(ints_packed)+1
         if (i<=size(ints_packed)) then
            int2=ints_packed(index(i))
         else
            int2=int1+1
         end if
         if (int2>int1) then
            if (i-start>1) then
               ! sort entries start to i-1 on unn
               call qsort(unn_packed(start:i-1), reindex(1:i-start))
               index(start:i-1)=index(start+reindex(1:i-start)-1)
            end if
            ! start new series
            start=i
            int1=int2
         else if (int2<int1) then
            FLAbort("Something went wrong in sorting")
         end if
      end do

      ! now change index from an index into the packed arrays to an index into the entire local node numbering
      do i=1, size(index)
         index(i) = fetch(layer_nodes, index(i))
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
