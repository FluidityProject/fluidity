#include "fdebug.h"

module hadapt_metric_based_extrude

  use elements
  use fields
  use sparse_tools
  use spud
  use metric_tools
  use vector_tools
  use meshdiagnostics
  use halos
  use vtk_interfaces
  use hadapt_extrude
  use hadapt_combine_meshes
  use global_parameters
  use interpolation_module
  implicit none

  public :: metric_based_extrude, recombine_metric, get_1d_mesh, get_1d_tensor

  contains

  subroutine metric_based_extrude(h_mesh, out_mesh, &
                                  full_metric, full_positions)
  !! Given a background mesh, and a metric on that background mesh,
  !! solve a load of 1d adaptivity problems for each column's vertical
  !! resolution.
    type(vector_field), intent(inout) :: h_mesh
    type(vector_field), intent(out) :: out_mesh
    ! the full metric, on the old full mesh
    ! (i.e. this metric is not on a mesh related to the current/new
    ! horizontal or extruded meshes)
    type(tensor_field), intent(in) :: full_metric
    type(vector_field), intent(in) :: full_positions
    
    !! A bunch of 1d meshes for each column in the background mesh
    !! and for each column in the adapted mesh
    !! We could assume here that all the columns of the background mesh
    !! are the same. However, to get a better adaptive result, we might
    !! want to adapt more than once with the same metric, so I won't make
    !! that assumption. Instead, we just assume that the background mesh
    !! is columnar.
    type(vector_field), dimension(node_count(h_mesh)) :: out_z_meshes
    type(vector_field) :: back_z_mesh, back_z_mesh_1d, out_z_mesh
    !! and the sizing field for each column:
    type(scalar_field) :: back_sizing
    type(tensor_field) :: back_z_metric

    integer :: column

    type(element_type) :: oned_shape
    type(quadrature_type) :: oned_quad
    integer :: quadrature_degree
    integer, parameter :: loc=2
    
    ! inline extrusion of the mesh
    integer :: n_regions
    logical :: apply_region_ids

    character(len=PYTHON_FUNC_LEN) :: sizing_function, depth_function
    logical:: sizing_is_constant, depth_is_constant
    real:: constant_sizing, depth, min_bottom_layer_frac
    
    integer :: d, r, stat
    integer, dimension(:), allocatable :: region_ids
    logical, dimension(node_count(h_mesh)) :: column_visited
    
    integer :: it, adapt_iterations
    
    logical :: replace_sizing_function
    character(len=PYTHON_FUNC_LEN) :: new_sizing_function
    logical:: new_sizing_is_constant
    real:: new_constant_sizing
    

    ewrite(1,*) "Inside metric_based_extrude"

    call get_option("/geometry/quadrature/degree", quadrature_degree)
    oned_quad = make_quadrature(vertices=loc, dim=1, degree=quadrature_degree)
    oned_shape = make_element_shape(vertices=loc, dim=1, degree=1, quad=oned_quad)
    call deallocate(oned_quad)

    replace_sizing_function = have_option(&
                    "/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity"//&
                    "/inhomogenous_vertical_resolution/initial_extrusion_sizing_function")
    if(replace_sizing_function) then
      ! we've selected to overrule the sizing options specified under the mesh
      ! (probably to give a finer mesh for the interpolation)
      call get_option("/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity"//&
                  "/inhomogenous_vertical_resolution/initial_extrusion_sizing_function/constant", &
                      new_constant_sizing, stat=stat)
      if (stat==0) then
        new_sizing_is_constant=.true.
      else
        new_sizing_is_constant=.false.
        call get_option("/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity"//&
                  "/inhomogenous_vertical_resolution/initial_extrusion_sizing_function/python", &
                        new_sizing_function, stat=stat)
        if (stat/=0) then
          FLAbort("Unknown way of specifying the extrusion sizing function in metric based extrusion")
        end if       
      end if
                  
    end if

    call get_option("/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity"//&
                    "/inhomogenous_vertical_resolution/adapt_iterations", &
                    adapt_iterations, default=1)

    n_regions = option_count(trim(full_positions%mesh%option_path)//&
                             '/from_mesh/extrude/regions')
    if(n_regions==0) then
      FLExit("No regions options found under extrude.")
    elseif(n_regions<0) then
      FLAbort("Negative number of regions options found under extrude.")
    end if
    apply_region_ids = (n_regions>1)
    column_visited = .false.

    do r = 0, n_regions-1
    
      call get_extrusion_options(trim(full_positions%mesh%option_path), &
                                 r, apply_region_ids, region_ids, &
                                 depth_is_constant, depth, depth_function, &
                                 sizing_is_constant, constant_sizing, sizing_function, &
                                 min_bottom_layer_frac)
                                 
      if(replace_sizing_function) then
        ! we've selected to overrule the sizing options specified under the mesh
        ! (probably to give a finer mesh for the interpolation)
        sizing_is_constant = new_sizing_is_constant
        if(sizing_is_constant) then
          constant_sizing = new_constant_sizing
        else
          sizing_function = new_sizing_function
        end if
      end if
            
      ! create a 1d vertical mesh under each surface node
      do column=1,node_count(h_mesh)
      
        if(skip_column_extrude(h_mesh%mesh, column, &
                               apply_region_ids, column_visited(column), region_ids)) cycle

        call compute_z_nodes(back_z_mesh_1d, node_val(h_mesh, column), &
                             min_bottom_layer_frac, &
                             depth_is_constant, depth, depth_function, &
                             sizing_is_constant, constant_sizing, sizing_function)

        do it = 1, adapt_iterations

          ! so now we have 1d vertical positions, we need a full version of this to do the interpolation
          call allocate(back_z_mesh, full_positions%dim, back_z_mesh_1d%mesh, "FullZMeshCoordinates")
          do d = 1, h_mesh%dim
            call set(back_z_mesh, d, node_val(h_mesh, d, column))
          end do
          call set_all(back_z_mesh, full_positions%dim, back_z_mesh_1d%val(1)%ptr)
          
          call allocate(back_z_metric, back_z_mesh%mesh, &
                        dim=full_metric%dim, name="BackgroundZMetric")
          
          ! temp. fix: old and new metric need same name for linear_interpolation -ask Patrick
          back_z_metric%name = full_metric%name
          call linear_interpolation( full_metric, full_positions, &
                                    back_z_metric, back_z_mesh )
          back_z_metric%name="BackgroundZMetric"
          
          call deallocate(back_z_mesh)
          
          call get_1d_sizing(back_z_metric, back_sizing)
          call deallocate(back_z_metric)
          
          call adapt_1d(back_z_mesh_1d, back_sizing, oned_shape, out_z_mesh)

          call deallocate(back_sizing)
          call deallocate(back_z_mesh_1d)
          
          back_z_mesh_1d = out_z_mesh
          call incref(back_z_mesh_1d)
          call deallocate(out_z_mesh)
          
        end do
        
        out_z_meshes(column) = back_z_mesh_1d
        call incref(out_z_meshes(column))
        call deallocate(back_z_mesh_1d)
        
      end do
      
      if(apply_region_ids) deallocate(region_ids)
        
    end do
      
    ! combine these into a full mesh
    call add_nelist(h_mesh%mesh)
    call combine_z_meshes(h_mesh, out_z_meshes, out_mesh, &
      ele_shape(full_positions, 1), full_positions%mesh%name, &
      trim(full_positions%mesh%option_path))

    call deallocate(oned_shape)
    do column=1,node_count(h_mesh)
      call deallocate(out_z_meshes(column))
    end do
    
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

    real, dimension(metric%dim) :: normal
    real :: mesh_size
    integer :: node

    call allocate(sizing, metric%mesh, "Back1DSizingFunction")

    ! normal here should be made smarter if this is on the globe.
    normal = 0.0
    normal(metric%dim) = 1.0

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
      call set(oned_tensor, spread(spread(oned_value, 1, oned_tensor%dim), 2, oned_tensor%dim))

    else
    
      nodes = row_length(back_columns, column)
      column_nodes => row_m_ptr(back_columns, column)

      do i=1,nodes
        j = column_nodes(i)

        oned_value = dot_product(matmul(normal, node_val(back_tensor, j)), normal)
        call set(oned_tensor, i, spread(spread(oned_value, 1, oned_tensor%dim), 2, oned_tensor%dim))
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
      call set(metric, metric%dim, metric%dim, j, oned_val(1,1))
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
      ewrite_minmax(z_mesh%val(1)%ptr)
      ewrite_minmax(z_mesh%mesh%region_ids)
    end if

  end subroutine adapt_1d

end module hadapt_metric_based_extrude
