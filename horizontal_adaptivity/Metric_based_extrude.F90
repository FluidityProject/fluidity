#include "fdebug.h"

module hadapt_metric_based_extrude

  use vector_tools
  use global_parameters
  use elements
  use spud
  use quicksort
  use data_structures
  use sparse_tools
  use metric_tools
  use fields
  use meshdiagnostics
  use vtk_interfaces
  use halos
  use hadapt_advancing_front
  use hadapt_combine_meshes
  use interpolation_module

  implicit none

  public :: metric_based_extrude, recombine_metric, get_1d_mesh, get_1d_tensor

  contains

  subroutine metric_based_extrude(h_positions_new, h_positions_old, out_mesh, &
                                  full_metric, full_positions, map)
  !! Given a background mesh, and a metric on that background mesh,
  !! solve a load of 1d adaptivity problems for each column's vertical
  !! resolution.
    type(vector_field), intent(inout) :: h_positions_new
    type(vector_field), intent(inout) :: h_positions_old ! inout to allow pickers caching
    type(vector_field), intent(out) :: out_mesh
    ! the full metric, on the old full mesh
    ! (i.e. this metric is not on a mesh related to the current/new
    ! horizontal or extruded meshes)
    type(tensor_field), intent(in) :: full_metric
    type(vector_field), intent(in) :: full_positions
    ! map from new nodal positions in h_position_new to 
    integer, dimension(:), intent(in), optional, target :: map
    
    !! A bunch of 1d meshes for each column in the background mesh
    !! and for each column in the adapted mesh
    !! We could assume here that all the columns of the background mesh
    !! are the same. However, to get a better adaptive result, we might
    !! want to adapt more than once with the same metric, so I won't make
    !! that assumption. Instead, we just assume that the background mesh
    !! is columnar.
    type(vector_field), dimension(node_count(h_positions_new)) :: out_z_meshes
    type(vector_field) :: back_z_mesh
    type(mesh_type) :: mesh
    type(scalar_field) :: back_sizing

    type(element_type) :: oned_shape
    type(quadrature_type) :: oned_quad
    integer :: quadrature_degree
    integer, parameter :: loc=2
    
    integer :: i, column, h_old_ele, old_ele, prev_old_ele, old_face

    type(csr_sparsity) :: columns_sparsity
    integer, dimension(:), pointer :: column_nodes
    integer, dimension(:), pointer :: lmap
    integer, dimension(:), pointer :: old_ele_nodes, old_ele_faces, &
                                      h_old_ele_nodes, neigh_old_eles, &
                                      old_ele_neigh
    
    real, dimension(h_positions_new%dim) :: new_column_positions
    
    type(integer_set) :: top_surface_nodes, bottom_surface_nodes
    type(integer_set) :: shared_elements
    
    real :: intersection_metric_comp, intersection_pos, mesh_size
    
    integer :: old_face_count, h_old_ele_loc
    type(integer_set) :: face_columns
    type(integer_hash_table) :: column_faces

    ! all of this is broken for variable element types!
    real, dimension(full_positions%dim, face_loc(full_positions, 1)) :: old_face_positions
    real, dimension(face_loc(full_positions, 1)) :: local_coords
    integer, dimension(face_loc(full_positions, 1)) :: current_face_global_nodes, prev_face_global_nodes
    type(integer_set), dimension(face_loc(full_positions, 1)) :: all_elements
    
    ewrite(1,*) "Inside metric_based_extrude"

    call get_option("/geometry/quadrature/degree", quadrature_degree)
    oned_quad = make_quadrature(vertices=loc, dim=1, degree=quadrature_degree)
    oned_shape = make_element_shape(vertices=loc, dim=1, degree=1, quad=oned_quad)
    call deallocate(oned_quad)
    
    if (present(map)) then
      assert(node_count(h_positions_new) == size(map))
      lmap => map
    else
      allocate(lmap(node_count(h_positions_new)))
      lmap = get_element_mapping(h_positions_old, h_positions_new, only_owned=.true.)
    end if
    
    call allocate(top_surface_nodes)
    call allocate(bottom_surface_nodes)

    ! we need to make sure that the columns are descending ordered
    ! so pass in full_positions
    call create_columns_sparsity(columns_sparsity, full_positions%mesh, positions=full_positions)
    
    do column=1,node_count(h_positions_old)
      column_nodes => row_m_ptr(columns_sparsity, column)
      call insert(top_surface_nodes, column_nodes(1))
      call insert(bottom_surface_nodes, column_nodes(size(column_nodes)))
    end do
    
    ! create a 1d vertical mesh under each surface node
    do column=1,node_count(h_positions_new)
    
      if(.not.node_owned(h_positions_new, column)) cycle
    
      new_column_positions = node_val(h_positions_new, column)
    
      ! In what element of the old horizontal mesh does the new column lie?
      h_old_ele = lmap(column)
      h_old_ele_loc = ele_loc(h_positions_old, h_old_ele)
      
      ! we need to find the element at the top of the full mesh
      ! first let's find out the nodes neighbouring this elements in the horizontal mesh
      h_old_ele_nodes => ele_nodes(h_positions_old, h_old_ele)
      
      call allocate(all_elements)
      do i = 1, size(h_old_ele_nodes)
        ! then we work out the corresponding top node in the full mesh
        column_nodes => row_m_ptr(columns_sparsity, h_old_ele_nodes(i))
        ! and work out the full mesh elements connected to that node
        neigh_old_eles => node_neigh(full_positions, column_nodes(1))
        ! insert it into a set
        call insert(all_elements(i), neigh_old_eles)
      end do
      
      call set_intersection(shared_elements, all_elements)
      
      ! there should now only be one entry in shared_elements
      if(key_count(shared_elements)/=1) then
        ewrite(-1,*) 'key_count(shared_elements) = ', key_count(shared_elements)
        FLAbort("Something has gone wrong finding the shared element.")
      end if
      
      ! and that entry will be the top element... woo!
      old_ele = fetch(shared_elements, 1)
      
      call deallocate(shared_elements)
      call deallocate(all_elements)

      ! find the top facet (which won't be between any of my elements so can't
      ! be included in the main loop)
      old_ele_nodes => ele_nodes(full_positions, old_ele)
      top_ele_node_loop: do i = 1, size(old_ele_nodes)
        if(.not.has_value(top_surface_nodes, old_ele_nodes(i))) exit top_ele_node_loop
      end do top_ele_node_loop
      assert(i<=size(old_ele_nodes))
      ! i is now the local face number of the face on the top surface... woo!
      
      ! convert this to a global face number
      old_ele_faces => ele_faces(full_positions, old_ele)
      old_face = old_ele_faces(i)
      
      call allocate(column_faces)
      old_face_count = 1
      ! insert the top face into the integer hash table
      call insert(column_faces, old_face_count, old_face)

      ! record the top element as visited
      prev_old_ele = old_ele
      
      infinite_loop: do
        ! in this loop we work our way down the column of elements
        ! let's hope we get to the exit criterion eventually!
      
        old_ele_faces => ele_faces(full_positions, old_ele)
        old_ele_neigh => ele_neigh(full_positions, old_ele)
        
        face_loop: do i = 1, size(old_ele_faces)
          call allocate(face_columns)
          call insert(face_columns, full_positions%mesh%columns(face_global_nodes(full_positions, old_ele_faces(i))))
          ! the face between old_ele and the next element is the one with as many vertices as
          ! the top element and isn't facing into the abyss or the previous element
          if ((key_count(face_columns)==h_old_ele_loc).and.&
              (old_ele_neigh(i)>0).and.(old_ele_neigh(i)/=prev_old_ele)) then
            exit face_loop ! we found it
          end if
          call deallocate(face_columns)
        end do face_loop
        
        if(i==size(old_ele_faces)+1) then
          ! we didn't find a face that meets our criterion so we must be at the bottom of the column
          exit infinite_loop
        else
          ! we found a face but didn't clean up...
          call deallocate(face_columns)
        end if
        
        ! remember where we were
        prev_old_ele = old_ele
        ! where we are now
        old_ele = old_ele_neigh(i)
        
        ! the face between them - remember it!
        old_face = old_ele_faces(i)
        old_face_count = old_face_count + 1
        call insert(column_faces, old_face_count, old_face)
      
      end do infinite_loop
      
      ! the search above won't have found the bottom face so let's search for it too
      old_ele_nodes => ele_nodes(full_positions, old_ele)
      bottom_ele_node_loop: do i = 1, size(old_ele_nodes)
        if(.not.has_value(bottom_surface_nodes, old_ele_nodes(i))) exit bottom_ele_node_loop
      end do bottom_ele_node_loop
      assert(i<=size(old_ele_nodes))
      ! i is now the local face number of the face on the bottom surface... woo!
      
      ! convert this to a global face number and save it
      old_ele_faces => ele_faces(full_positions, old_ele)
      old_face = old_ele_faces(i)
      old_face_count = old_face_count + 1
      call insert(column_faces, old_face_count, old_face)
      assert(key_count(column_faces)==old_face_count)


      ! allocate the mesh and the fields we want
      call allocate(mesh, old_face_count, old_face_count-1, oned_shape, "VerticalIntersectionMesh")
      call allocate(back_z_mesh, 1, mesh, "BackVerticalMesh")
      call allocate(back_sizing, mesh, "BackSizing")
      call deallocate(mesh)


      ! work out the local coordinates on the top face
      ! only works for linear coordinates at the moment as adaptivity can
      ! only deal with this anyway
      local_coords(1:h_positions_new%dim) = new_column_positions
      local_coords(h_positions_new%dim+1) = 1.0
    
      ! and find the values of the positions
      old_face = fetch(column_faces, 1)
      old_face_positions = face_val(full_positions, old_face)
      ! wipe out the z component - i.e. assuming a horizontal face on the top surface
      old_face_positions(h_positions_new%dim+1, :) = 1.0

      ! and hey presto... we have the local_coords
      call solve(old_face_positions, local_coords)


      ! work out the position of the intersection of the column with the top face
      intersection_pos = face_eval_field(old_face, full_positions, full_positions%dim, local_coords)
      call set(back_z_mesh, 1, (/intersection_pos/))
      ! work out the value of the dim, dim component of the full_metric at the intersection point
      ! NOTE WELL: this is where we assume that up is in the vertical direction (i.e. not
      ! safe on the sphere but that's not supported yet anyway and this is assumed in
      ! lots of places)
      intersection_metric_comp = face_eval_field(old_face, full_metric, full_metric%dim(1), full_metric%dim(2), local_coords)
      mesh_size = edge_length_from_eigenvalue(intersection_metric_comp)
      call set(back_sizing, 1, mesh_size)

      ! finally find the global nodes of the first face
      prev_face_global_nodes = face_global_nodes(full_positions, old_face)

      do i = 2, old_face_count
      
        old_face = fetch(column_faces, i)
      
        current_face_global_nodes = face_global_nodes(full_positions, old_face)
        
        call permute_local_coords(local_coords, current_face_global_nodes, prev_face_global_nodes)
        
        ! work out the position of the intersection of the column with the old face
        intersection_pos = face_eval_field(old_face, full_positions, full_positions%dim, local_coords)
        call set(back_z_mesh, i, (/intersection_pos/))
        ! work out the value of the dim, dim component of the full_metric at the intersection point
        ! NOTE WELL: this is where we assume that up is in the vertical direction (i.e. not
        ! safe on the sphere but that's not supported yet anyway and this is assumed in
        ! lots of places)
        intersection_metric_comp = face_eval_field(old_face, full_metric, full_metric%dim(1), full_metric%dim(2), local_coords)
        mesh_size = edge_length_from_eigenvalue(intersection_metric_comp)
        call set(back_sizing, i, mesh_size)

        ! record the current global nodes as the previous ones for the next iteration
        prev_face_global_nodes = current_face_global_nodes
      
      end do
      
      ! and we're done... we've worked out the intersection positions in back_z_mesh and the intersection sizing
      ! function in back_sizing.
      call deallocate(column_faces)
      
      ! now let's adapt that and put the adapted mesh in out_z_meshes(column)...
      call adapt_1d(back_z_mesh, back_sizing, oned_shape, out_z_meshes(column))
      
      ! we're done with all that beautiful intersection work
      call deallocate(back_z_mesh)
      call deallocate(back_sizing)

    end do
    
    if (.not. present(map)) then
      deallocate(lmap)
    end if
    call deallocate(columns_sparsity)
    call deallocate(top_surface_nodes)
    call deallocate(bottom_surface_nodes)
      
    ! combine these into a full mesh
    call add_nelist(h_positions_new%mesh)
    call combine_z_meshes(h_positions_new, out_z_meshes, out_mesh, &
      ele_shape(full_positions, 1), full_positions%mesh%name, &
      trim(full_positions%mesh%option_path))

    call deallocate(oned_shape)
    do column=1,node_count(h_positions_new)
      if(node_owned(h_positions_new, column)) then
        call deallocate(out_z_meshes(column))
      end if
    end do
    
  contains
  
    subroutine permute_local_coords(local_coords, current_face_global_nodes, prev_face_global_nodes)
      real, dimension(:), intent(inout)  :: local_coords
      integer, dimension(:), intent(in) :: current_face_global_nodes
      integer, dimension(:), intent(in) :: prev_face_global_nodes
      
      integer :: g1, g2, missing_ind
      integer, dimension(size(current_face_global_nodes)) :: permutation, notfound
      
      assert(size(current_face_global_nodes)==size(prev_face_global_nodes))
      assert(size(current_face_global_nodes)==size(local_coords))
      
      ! because the local face numbers may have shifted between one face and the next
      ! it is necessary to work out a permutation for the local coords to be valid
      permutation = -1
      notfound = 1
      do g1 = 1, size(prev_face_global_nodes)
        do g2 = 1, size(current_face_global_nodes)
          if(prev_face_global_nodes(g1)==current_face_global_nodes(g2)) then
            permutation(g2) = g1
            notfound(g2) = 0
            exit
          end if
        end do
        if(g2==size(prev_face_global_nodes)+1) then
          missing_ind = g1
        end if
      end do
      assert(sum(notfound)==1) ! debugging check that only one position hasn't been found
      permutation(minloc(permutation)) = missing_ind
      assert(all(permutation>0))
      
      ! now permute the local coordinates so that the local node ordering is
      ! correct
      call apply_permutation(local_coords, permutation)

    end subroutine permute_local_coords
    
  end subroutine metric_based_extrude
    
  subroutine get_1d_mesh(column, back_mesh, back_columns, metric, oned_shape, z_mesh, &
                         sizing)
    integer, intent(in) :: column
    type(vector_field), intent(in) :: back_mesh
    type(csr_sparsity), intent(in) :: back_columns
    type(tensor_field), intent(in) :: metric
    type(element_type), intent(inout) :: oned_shape
    type(vector_field), intent(out) :: z_mesh
    type(scalar_field), intent(out), optional :: sizing

    type(mesh_type) :: mesh

    integer :: nodes, elements
    integer, dimension(:), pointer :: column_nodes
    integer :: i, j
    integer :: dim
    integer, parameter :: loc=2

    real, dimension(mesh_dim(back_mesh)) :: normal
    real :: mesh_size

    nodes = row_length(back_columns, column)
    elements = nodes - 1
    dim = mesh_dim(back_mesh)
    column_nodes => row_m_ptr(back_columns, column)

    call allocate(mesh, nodes, elements, oned_shape, "Mesh")
    do i=1,elements
      mesh%ndglno((i-1) * loc + 1: i*loc) = (/i, i+1/)
    end do

    call allocate(z_mesh, 1, mesh, "ZMesh")
    if(present(sizing)) then
      call allocate(sizing, mesh, "SizingFunction")
    end if
    call deallocate(mesh)

    ! normal here should be made smarter if this is on the globe.
    normal = 0.0
    normal(dim) = 1.0

    do i=1,nodes
      j = column_nodes(i)
      call set(z_mesh, i, (/node_val(back_mesh, dim, j)/))

      if(present(sizing)) then
        mesh_size = edge_length_from_eigenvalue(dot_product(matmul(normal, node_val(metric, j)), normal))
        call set(sizing, i, mesh_size)
      end if
    end do
    
  end subroutine get_1d_mesh
  
  subroutine get_1d_sizing(metric, sizing)
    ! this is a full dimensional metric defined on a 1d mesh
    type(tensor_field), intent(inout) :: metric
    ! and we want to get the sizing based on the vertical component of the metric
    type(scalar_field), intent(out) :: sizing

    real, dimension(metric%dim(1)) :: normal
    real :: mesh_size
    integer :: node

    call allocate(sizing, metric%mesh, "Back1DSizingFunction")

    ! normal here should be made smarter if this is on the globe.
    normal = 0.0
    normal(metric%dim(1)) = 1.0

    do node=1,node_count(sizing)
      mesh_size = edge_length_from_eigenvalue(dot_product(matmul(normal, node_val(metric, node)), normal))
      call set(sizing, node, mesh_size)
    end do
    
  end subroutine get_1d_sizing
  
  subroutine get_1d_tensor(column, back_tensor, oned_tensor, back_columns)
    integer, intent(in) :: column
    type(tensor_field), intent(in) :: back_tensor
    type(tensor_field), intent(inout) :: oned_tensor
    type(csr_sparsity), intent(in) :: back_columns

    integer :: nodes
    integer, dimension(:), pointer :: column_nodes
    integer :: i, j
    integer :: dim

    real, dimension(mesh_dim(back_tensor)) :: normal
    real :: oned_value

    ! normal here should be made smarter if this is on the globe.
    dim = mesh_dim(back_tensor)
    normal = 0.0
    normal(dim) = 1.0

    if(back_tensor%field_type==FIELD_TYPE_CONSTANT) then

      oned_value = dot_product(matmul(normal, node_val(back_tensor, 1)), normal)
      call set(oned_tensor, spread(spread(oned_value, 1, oned_tensor%dim(1)), 2, oned_tensor%dim(1)))

    else
    
      nodes = row_length(back_columns, column)
      column_nodes => row_m_ptr(back_columns, column)

      do i=1,nodes
        j = column_nodes(i)

        oned_value = dot_product(matmul(normal, node_val(back_tensor, j)), normal)
        call set(oned_tensor, i, spread(spread(oned_value, 1, oned_tensor%dim(1)), 2, oned_tensor%dim(1)))
      end do
      
    end if
    
  end subroutine get_1d_tensor
  
  subroutine recombine_metric(metric, column, oned_metric, back_columns)
    type(tensor_field), intent(inout) :: metric
    integer, intent(in) :: column
    type(tensor_field), intent(in) :: oned_metric
    type(csr_sparsity), intent(in) :: back_columns
    
    integer :: nodes
    integer, dimension(:), pointer :: column_nodes
    integer :: i, j
    real, dimension(1, 1) :: oned_val

    nodes = row_length(back_columns, column)
    column_nodes => row_m_ptr(back_columns, column)
    
    ! NOTE WELL: just as in get_1d_mesh, we're about to assume that the
    ! 1d metric belongs in the last entry of the full metric!
    ! We should be cleverer than this (i.e. on a sphere).
    do i = 1, nodes
      j = column_nodes(i)
      
      oned_val = node_val(oned_metric, i)
      call set(metric, metric%dim(1), metric%dim(2), j, oned_val(1,1))
    end do
  
  end subroutine recombine_metric

  function get_expected_elements(z_mesh, sizing) result(elements)
    type(vector_field), intent(inout) :: z_mesh
    type(scalar_field), intent(in) :: sizing

    type(scalar_field) :: sizing_inverse
    integer :: elements

    call allocate(sizing_inverse, z_mesh%mesh, trim(sizing%name) // "Inverse")
    call invert(sizing, sizing_inverse)
    elements = field_integral(sizing_inverse, z_mesh)
    call deallocate(sizing_inverse)
    assert(elements > 0)
  end function get_expected_elements

  subroutine adapt_1d(back_mesh, sizing, oned_shape, z_mesh, preserve_regions)
    type(vector_field), intent(in) :: back_mesh
    type(scalar_field), intent(in) :: sizing
    type(vector_field), intent(inout) :: z_mesh
    type(element_type), intent(inout) :: oned_shape
    logical, intent(in), optional :: preserve_regions

    integer :: elements
    integer :: node

    type(mesh_type) :: mesh

    real, dimension(:), allocatable :: metric_step_length
    real, dimension(element_count(back_mesh)) :: desired_ele_lengths
    real, dimension(element_count(back_mesh)) :: metric_ele_lengths
    integer :: old_node_counter, new_node_counter
    real :: old_metric_back, old_metric_front, new_metric_back, new_metric_front
    real :: new_node_position, real_step_length
    
    logical :: l_preserve_regions
    integer, dimension(:), allocatable :: nodes_per_region, tmp_region_bdy_nodes, region_bdy_nodes
    integer, dimension(:), allocatable :: tmp_region_ids, region_ids
    integer :: ele, ele_2, ni, no_region_bdys, face, region_bdy
    integer :: total_nodes, tmp_size
    integer, dimension(:), pointer :: neigh
    integer, dimension(1) :: node_array
    
    l_preserve_regions = present_and_true(preserve_regions)

    ! don't make the decision to preserve regions based on
    ! the options tree because then it would happen during
    ! mesh extrusion as well as 1d adaptivity, which would
    ! be a wasted effort (also region ids might not be available)
    if(l_preserve_regions) then
      assert(associated(back_mesh%mesh%region_ids))
      
      ! let's hope the region ids don't go too high!
      ! needs to have a minimum length of 2 (for both of the ends)
      ! but we're probably going to overestimate the size here
      ! (especially as we're catering for the possibility of negative
      !  region ids which may not even be possible!)
      tmp_size = abs(minval(back_mesh%mesh%region_ids)) + &
                 abs(maxval(back_mesh%mesh%region_ids)) + 2
      allocate(tmp_region_bdy_nodes(tmp_size))
      tmp_region_bdy_nodes = 0
      
      allocate(tmp_region_ids(tmp_size))
      tmp_region_ids = -1
      
      ! find the depths of the boundaries between region ids
      no_region_bdys = 0
      do ele = 1, ele_count(back_mesh)
        neigh => ele_neigh(back_mesh, ele)
        do ni = 1, size(neigh)
          ele_2 = neigh(ni)
          if (ele_2>0) then
            ! Internal faces only.
            if(back_mesh%mesh%region_ids(ele)/=back_mesh%mesh%region_ids(ele_2)) then
              face=ele_face(back_mesh, ele, ele_2)
              node_array=face_global_nodes(back_mesh, face)
              if(.not.(any(node_array(1)==tmp_region_bdy_nodes))) then
                ! only include this node if we haven't visited it from
                ! another element
                no_region_bdys = no_region_bdys + 1
                tmp_region_bdy_nodes(no_region_bdys) = node_array(1) ! remember if we've visited this node
                
                ! always take the region_id from the region higher up
                if((0.5*sum(ele_val(back_mesh, 1, ele)))>(0.5*sum(ele_val(back_mesh, 1, ele_2)))) then
                  tmp_region_ids(no_region_bdys) = back_mesh%mesh%region_ids(ele)
                else
                  tmp_region_ids(no_region_bdys) = back_mesh%mesh%region_ids(ele_2)
                end if
                
              end if
            end if
          else
            ! External faces get added too but they're easier - they definitely get counted as region_bdys
            face = ele_face(back_mesh, ele, ele_2)
            node_array=face_global_nodes(back_mesh, face)
            ! should be no need to check if we've visited this face already
            no_region_bdys = no_region_bdys + 1
            tmp_region_bdy_nodes(no_region_bdys) = node_array(1)
            tmp_region_ids(no_region_bdys) = back_mesh%mesh%region_ids(ele)
          end if          
        end do
      end do
      
      ! check the region bdys have been found at the ends of the domain
      ! (this uses the assumption that the back_mesh is coordinate ordered)
      assert(maxval(tmp_region_bdy_nodes)==node_count(back_mesh))
      assert(no_region_bdys>1)
      
      ! take off 1 region bdy for the first node
      no_region_bdys = no_region_bdys - 1 
            
      allocate(region_bdy_nodes(0:no_region_bdys))
      region_bdy_nodes = 0
      region_bdy_nodes(0) = 1 ! include the first node
      
      ! there's a region id for every bdy except the first node (index 0 in region_bdy_nodes)
      ! these ids correspond to the region_id from the region higher up than the boundary
      ! (i.e. we ditch a region_id from tmp_region_ids from the top node)
      allocate(region_ids(no_region_bdys))
      region_ids = -1
      
      ! sort the region bdy nodes into increasing order (corresponds to decreasing depth order)
      do region_bdy = no_region_bdys, 1, -1
        node_array = maxloc(tmp_region_bdy_nodes)
        if(tmp_region_bdy_nodes(node_array(1))>0) then
          region_bdy_nodes(region_bdy) = tmp_region_bdy_nodes(node_array(1))
          tmp_region_bdy_nodes(node_array(1)) = 0 ! blank it so we don't find it again
          ! assign the region id to this boundary from the region higher up
          region_ids(region_bdy) = tmp_region_ids(node_array(1))
        end if
      end do
      ! check the region bdys have been found at the ends of the domain
      ! (this uses the assumption that the back_mesh is coordinate ordered)
      assert(maxval(tmp_region_bdy_nodes)==1) ! should only be the first node remaining (everything else should be 0)
      assert(maxval(region_bdy_nodes)==node_count(back_mesh))
      assert(minval(region_bdy_nodes)==1)
      assert(all(region_bdy_nodes>0))
      assert(all(region_ids>=0))
            
      deallocate(tmp_region_bdy_nodes)
      deallocate(tmp_region_ids)
      
      ewrite(2,*) 'in adapt_1d'
      ewrite(2,*) 'no_region_bdys = ', no_region_bdys
      ewrite(2,*) 'region_bdy_nodes = ', region_bdy_nodes
      ewrite(2,*) 'region_ids = ', region_ids
    
    else
      no_region_bdys = 1
      allocate(region_bdy_nodes(0:no_region_bdys))
      region_bdy_nodes(0) = 1 ! include the first node
      region_bdy_nodes(no_region_bdys) = node_count(back_mesh) ! include the last node
    end if
    
    ! First we need to see how many nodes we will have.
    ! I do this by basically doing the work twice.
    ! You could be more clever and record the steps and positions,
    ! but I don't have dynamically-sized arrays :-(

    allocate(nodes_per_region(0:no_region_bdys))
    nodes_per_region = 0
    
    ! the first node
    nodes_per_region(0) = 1
    
    allocate(metric_step_length(no_region_bdys))
    metric_step_length = 1.0

    do region_bdy = 1, no_region_bdys
      do node = region_bdy_nodes(region_bdy-1), region_bdy_nodes(region_bdy)-1
        ! project the sizing function to an array over the old elements
        desired_ele_lengths(node) = 0.5*(node_val(sizing, node)+node_val(sizing, node+1))
        ! translate the current element lengths into metric space by dividing by the desired elemental length
        metric_ele_lengths(node) = (abs(node_val(back_mesh, 1, node)-node_val(back_mesh, 1, node+1))) &
                                    /desired_ele_lengths(node)
      end do
      ! work out the number of nodes in this region by summing the metric lengths (then rounding up to ensure an integer value)
      nodes_per_region(region_bdy) = ceiling(sum(metric_ele_lengths(region_bdy_nodes(region_bdy-1):(region_bdy_nodes(region_bdy)-1))))
      ! work out the step length in metric space (ideal is 1 but this will be less than that due to rounding up on previous line)
      metric_step_length(region_bdy) = (sum(metric_ele_lengths(region_bdy_nodes(region_bdy-1):(region_bdy_nodes(region_bdy)-1)))) &
                                        /nodes_per_region(region_bdy)
    end do

    total_nodes = sum(nodes_per_region)
    elements = total_nodes - 1
    call allocate(mesh, total_nodes, elements, oned_shape, "Mesh")
    if(l_preserve_regions) then
      allocate(mesh%region_ids(elements))
      mesh%region_ids = 0
    end if
    call allocate(z_mesh, 1, mesh, "AdaptedZMesh")
    call set(z_mesh, (/huge(0.0)/)) ! a bug catcher
    call deallocate(mesh)
    
    call set(z_mesh, region_bdy_nodes(0), node_val(back_mesh, 1, region_bdy_nodes(0)))

    do region_bdy = 1, no_region_bdys
      ! the node at the start of this region in the old mesh
      old_node_counter = region_bdy_nodes(region_bdy-1)
      ! the front and back positions of the next old element 
      ! (in metric space and relative to the start of this region)
      old_metric_back = 0.0
      old_metric_front = metric_ele_lengths(old_node_counter)

      ! the node at the start of this region in the new mesh
      new_node_counter = sum(nodes_per_region(0:region_bdy-1))
      
      ! the last new node position (in real space at the start of this region)
      new_node_position = node_val(back_mesh, 1, old_node_counter)
      
      ! the front and back positions of the next new element 
      ! (in metric space and relative to the start of this region)
      new_metric_back = 0.0
      new_metric_front = metric_step_length(region_bdy)
      
      do node = 1, nodes_per_region(region_bdy)-1
        new_node_counter = new_node_counter + 1
        
        real_step_length = 0.0
        
        if(l_preserve_regions) then
          ! here we assume that the elements are also ordered
          ! this subroutine doesn't take care of the node to element list
          ! so anything that does will have to reorder the region_id list as
          ! well if it doesn't have the same assumption
          z_mesh%mesh%region_ids(new_node_counter-1) = region_ids(region_bdy)
        end if
        
        do while (old_metric_front < new_metric_front)
          ! add an increment of real space to the position 
          ! (this equals the desired edge length times the metric step length)
          if((old_metric_front-old_metric_back)>(old_metric_front-new_metric_back)) then
            ! in this case the new element straddles an element boundary in the old mesh
            real_step_length = real_step_length + &
                              desired_ele_lengths(old_node_counter)*(old_metric_front-new_metric_back)
          else
            ! in this case the old element falls entirely within an element of the new mesh
            real_step_length = real_step_length + &
                              desired_ele_lengths(old_node_counter)*(old_metric_front-old_metric_back)
          end if

          ! move the back of the old element to the current position of the front
          old_metric_back = old_metric_front
          ! and then move the front on to the next element (i.e. step through the old elements)
          old_node_counter = old_node_counter + 1
          old_metric_front = old_metric_front + metric_ele_lengths(old_node_counter)
        end do
        
        ! so now the old_metric_front is ahead of the new_metric_front, we need to know by how much
        ! so we can add in that contribution to the new element edge length
        if((new_metric_front-new_metric_back)>(new_metric_front-old_metric_back)) then
          ! in this case the new element straddles an element boundary in the old mesh
          real_step_length = real_step_length + &
                             desired_ele_lengths(old_node_counter)*(new_metric_front-old_metric_back)
        else
          ! in this case the new element falls entirely within an old element
          real_step_length = real_step_length + &
                             desired_ele_lengths(old_node_counter)*(new_metric_front-new_metric_back)
        end if

        new_node_position = new_node_position - real_step_length
        
        call set(z_mesh, new_node_counter, (/new_node_position/))
        
        new_metric_back = new_metric_front
        new_metric_front = new_metric_front+metric_step_length(region_bdy)
      
      end do

      ! include the bottom node in this depth (but not the top node)
      new_node_counter = new_node_counter + 1
      call set(z_mesh, new_node_counter, (/node_val(back_mesh, 1, region_bdy_nodes(region_bdy))/))
      if(l_preserve_regions) then
        ! here we assume that the elements are also ordered
        ! this subroutine doesn't take care of the node to element list
        ! so anything that does will have to reorder the region_id list as
        ! well if it doesn't have the same assumption
        z_mesh%mesh%region_ids(new_node_counter-1) = region_ids(region_bdy)
      end if

    end do

    if(l_preserve_regions) then
      ewrite_minmax(z_mesh)
      ewrite_minmax(z_mesh%mesh%region_ids)
    end if

  end subroutine adapt_1d

end module hadapt_metric_based_extrude
