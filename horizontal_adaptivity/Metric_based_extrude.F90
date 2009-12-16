#include "fdebug.h"

module hadapt_metric_based_extrude

  use elements
  use fields
  use sparse_tools
  use hadapt_advancing_front
  use spud
  use metric_tools
  implicit none

  public :: metric_based_extrude, combine_z_meshes

  contains

  subroutine metric_based_extrude(h_mesh, back_mesh, metric, out_mesh)
  !! Given a background mesh, and a metric on that background mesh,
  !! solve a load of 1d adaptivity problems for each column's vertical
  !! resolution.
    type(vector_field), intent(inout) :: h_mesh
    type(vector_field), intent(in) :: back_mesh
    type(tensor_field), intent(in) :: metric
    type(vector_field), intent(out) :: out_mesh

    !! A bunch of 1d meshes for each column in the background mesh
    !! and for each column in the adapted mesh
    !! We could assume here that all the columns of the background mesh
    !! are the same. However, to get a better adaptive result, we might
    !! want to adapt more than once with the same metric, so I won't make
    !! that assumption. Instead, we just assume that the background mesh
    !! is columnar.
    type(vector_field), dimension(node_count(h_mesh)) :: back_z_meshes, out_z_meshes
    !! and the sizing field for each column:
    type(scalar_field), dimension(node_count(h_mesh)) :: back_sizing

    integer :: column

    type(csr_sparsity) :: back_columns
    type(element_type) :: oned_shape
    type(quadrature_type) :: oned_quad
    integer :: quadrature_degree
    integer, parameter :: loc=2
    
    ewrite(1,*) "Inside metric_based_extrude"

    ewrite(2,*) "In metric_based_extrude"

    call get_option("/geometry/quadrature/degree", quadrature_degree)
    oned_quad = make_quadrature(vertices=loc, dim=1, degree=quadrature_degree)
    oned_shape = make_element_shape(vertices=loc, dim=1, degree=1, quad=oned_quad)
    call deallocate(oned_quad)

    call add_nelist(h_mesh%mesh)
    
    call create_columns_sparsity(back_columns, back_mesh%mesh)

    ! create a 1d vertical mesh under each surface node
    do column=1,node_count(h_mesh)
      call get_1d_mesh(column, back_mesh, back_columns, metric, oned_shape, back_z_meshes(column), back_sizing(column))
      call adapt_1d(back_z_meshes(column), back_sizing(column), oned_shape, out_z_meshes(column))
    end do
      
    call deallocate(back_columns)
      
    ! combine these into a full mesh
    call combine_z_meshes(h_mesh, out_z_meshes, out_mesh, &
      ele_shape(back_mesh, 1), back_mesh%name)

    call deallocate(oned_shape)
    do column=1,node_count(h_mesh)
      call deallocate(back_z_meshes(column))
      call deallocate(out_z_meshes(column))
      call deallocate(back_sizing(column))
    end do
    
    if (has_faces(h_mesh%mesh)) then
      ! if the horizontal mesh has a surface mesh
      ! give one to the extruded mesh as well
      call generate_extruded_surface_mesh(h_mesh, out_mesh)
    end if
    
  end subroutine metric_based_extrude
    
  subroutine combine_z_meshes(h_mesh, z_meshes, out_mesh, &
    full_shape, mesh_name)
  !! Given the h_mesh and a z_mesh under each node of it combines these
  !! into a full horiz+vertic. mesh
  type(vector_field), intent(inout):: h_mesh
  type(vector_field), dimension(:), intent(in):: z_meshes
  type(vector_field), intent(out):: out_mesh
  type(element_type), intent(in):: full_shape
  !! the name of the topol. mesh to be created, not the coordinate field
  character(len=*), intent(in):: mesh_name  
  
    type(csr_sparsity):: out_columns
    type(mesh_type):: mesh
    integer:: column, total_out_nodes, total_out_elements, z_elements, last_seen
    
    total_out_nodes = 0
    ! For each column,
    ! add (number of nodes hanging off (== number of 1d elements)) * (element connectivity of chain)
    ! to compute the number of elements the extrusion routine will produce
    total_out_elements = 0
    do column=1, size(z_meshes)
      z_elements = ele_count(z_meshes(column))
      total_out_nodes = total_out_nodes + z_elements + 1
      assert(associated(h_mesh%mesh%adj_lists))
      assert(associated(h_mesh%mesh%adj_lists%nelist))
      total_out_elements = total_out_elements + z_elements * row_length(h_mesh%mesh%adj_lists%nelist, column)
    end do

    call allocate(mesh, total_out_nodes, total_out_elements, &
      full_shape, mesh_name)
    ! allocate mapping between extruded nodes to surface node (column number)
    ! it lies under
    allocate(mesh%columns(total_out_nodes))
    if (mesh_name=="CoordinateMesh") then
      call allocate(out_mesh, mesh_dim(h_mesh)+1, mesh, "Coordinate")
    else
      call allocate(out_mesh, mesh_dim(h_mesh)+1, mesh, trim(mesh_name)//"Coordinate")
    end if
    call deallocate(mesh)

    last_seen = 0
    do column=1,node_count(h_mesh)
      call append_to_structures(column, z_meshes(column), h_mesh, out_mesh, last_seen)
    end do
      
    call create_columns_sparsity(out_columns, out_mesh%mesh)
    call generate_layered_mesh(out_mesh, h_mesh)
    call deallocate(out_columns)
    
  end subroutine combine_z_meshes

  subroutine get_1d_mesh(column, back_mesh, back_columns, metric, oned_shape, z_mesh, sizing)
    integer, intent(in) :: column
    type(vector_field), intent(in) :: back_mesh
    type(csr_sparsity), intent(in) :: back_columns
    type(tensor_field), intent(in) :: metric
    type(element_type), intent(inout) :: oned_shape
    type(vector_field), intent(out) :: z_mesh
    type(scalar_field), intent(out) :: sizing
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
    call allocate(sizing, mesh, "SizingFunction")
    call deallocate(mesh)

    ! normal here should be made smarter if this is on the globe.
    normal = 0.0
    normal(dim) = 1.0

    do i=1,nodes
      j = column_nodes(i)
      call set(z_mesh, i, (/node_val(back_mesh, dim, j)/))

      mesh_size = edge_length_from_eigenvalue(dot_product(matmul(normal, node_val(metric, j)), normal))
      call set(sizing, i, mesh_size)
    end do
    
  end subroutine get_1d_mesh

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

  subroutine adapt_1d(back_mesh, sizing, oned_shape, z_mesh)
    type(vector_field), intent(in) :: back_mesh
    type(scalar_field), intent(in) :: sizing
    type(vector_field), intent(inout) :: z_mesh
    type(element_type), intent(inout) :: oned_shape

    integer :: elements
    integer :: nodes
    integer :: node
    integer :: back_node
    real :: step
    real :: current_pos
    real :: depth

    type(mesh_type) :: mesh
    ! to prevent infinitesimally thin bottom layer if sizing function
    ! is an integer mulitple of total depth, the bottom layer needs
    ! to have at least this fraction of the layer depth above it:
    real, parameter:: MIN_BOTTOM_LAYER_FRAC=1e-3

    ! First we need to see how many nodes we will have.
    ! I do this by basically doing the work twice.
    ! You could be more clever and record the steps and positions,
    ! but I don't have dynamically-sized arrays :-(

    nodes = 1
    step = node_val(sizing, 1)
    back_node = 2
    current_pos = node_val(back_mesh, 1, 1)
    depth = node_val(back_mesh, 1, node_count(back_mesh))
    do 
      current_pos = current_pos - step

      if (current_pos <= depth + MIN_BOTTOM_LAYER_FRAC*step) then
        exit
      end if
      nodes = nodes + 1

      do while (current_pos < node_val(back_mesh, 1, back_node))
        back_node = back_node + 1
      end do
      step = compute_step(current_pos, back_mesh, back_node, sizing)
    end do
    nodes = nodes + 1

    elements = nodes - 1
    call allocate(mesh, nodes, elements, oned_shape, "Mesh")
    call allocate(z_mesh, 1, mesh, "AdaptedZMesh")
    call deallocate(mesh)

    call set(z_mesh, 1, node_val(back_mesh, 1))
    step = node_val(sizing, 1)
    back_node = 2

    do node=2,elements
      call set(z_mesh, node, (/node_val(z_mesh, 1, node-1) - step/))

      ! Now we need to compute the next step:
      do while (node_val(z_mesh, 1, node) < node_val(back_mesh, 1, back_node))
        back_node = back_node + 1
      end do
      step = compute_step(node_val(z_mesh, 1, node), back_mesh, back_node, sizing)
    end do

    ! and set the depth to be the constant, taken from the last node of back_mesh
    call set(z_mesh, nodes, node_val(back_mesh, node_count(back_mesh)))
    assert(node_val(z_mesh, 1, nodes-1) > node_val(z_mesh, 1, nodes))

    contains
      function compute_step(current_pos, back_mesh, back_node, sizing) result(step)
        real, intent(in) :: current_pos
        type(vector_field), intent(in) :: back_mesh
        integer, intent(in) :: back_node
        type(scalar_field), intent(in) :: sizing

        real :: step
        real :: local_coord
        ! Compute the local coordinate of the point we're currently at in back_ele.
        local_coord = (current_pos - node_val(back_mesh, 1, back_node - 1)) / &
                      (node_val(back_mesh, 1, back_node) - node_val(back_mesh, 1, back_node - 1))
        step = node_val(sizing, back_node - 1) + (node_val(sizing, back_node) - node_val(sizing, back_node - 1)) * local_coord
      end function compute_step

  end subroutine adapt_1d

  subroutine append_to_structures(column, z_mesh, h_mesh, out_mesh, last_seen)
    integer, intent(in) :: column
    type(vector_field), intent(in) :: z_mesh, h_mesh
    type(vector_field), intent(inout) :: out_mesh
    integer, intent(inout) :: last_seen

    integer :: j
    integer :: h_dim, v_dim

    real, dimension(mesh_dim(out_mesh)) :: pos

    h_dim = mesh_dim(h_mesh)
    v_dim = mesh_dim(out_mesh)

    do j=1,node_count(z_mesh)
      last_seen = last_seen + 1
      out_mesh%mesh%columns(last_seen)=column
      pos(1:h_dim) = node_val(h_mesh, column)
      pos(v_dim) = node_val(z_mesh, 1, j)
      call set(out_mesh, last_seen, pos)
    end do
    
  end subroutine append_to_structures

  subroutine generate_extruded_surface_mesh(h_mesh, out_mesh)
    !!< Create a surface mesh for the extruded mesh and extrapolates
    !!< the boundary ids of the surface mesh of the horizontal mesh
    !!< on it. h_mesh needs to have a surface mesh before calling this
    type(vector_field), intent(inout):: h_mesh
    type(vector_field), intent(inout):: out_mesh
    
    integer, dimension(:), allocatable:: hnode2surfnode
    integer top_boundary_id, bottom_boundary_id
    integer i, hnode
    
    if (.not. has_faces(h_mesh%mesh)) then
      FLAbort("Horizontal mesh needs to have a surface mesh for extrusion")
    end if
    
    ! that's easy isn't it:
    call add_faces(out_mesh%mesh)
    
    ! now for the complicated part, "extruding" the boundary ids
    
    ! first find out if top and bottom boundary markers have been set
    call get_option(trim(out_mesh%mesh%option_path)//"/from_mesh/extrude/bottom_surface_id", &
      bottom_boundary_id, default=0)
    call get_option(trim(out_mesh%mesh%option_path)//"/from_mesh/extrude/top_surface_id", &
      top_boundary_id, default=0)
        
    ! now we need a mapping between nodes in the horizontal mesh and its
    ! node number in the surface mesh of the horizontal mesh (i.e. the edges
    ! around the horizontal mesh), interior nodes all map to 0
    allocate( hnode2surfnode(1: node_count(h_mesh)) )
    hnode2surfnode=0
    ! loop over the surface nodes of the horizontal mesh
    do i=1, size(h_mesh%mesh%faces%surface_node_list)
      ! its node number in the full horizontal mesh
      hnode=h_mesh%mesh%faces%surface_node_list(i)
      hnode2surfnode(hnode)=i
    end do

    ! this node element list provides a mapping between horizontal
    ! surface nodes and horizontal surface elements
    call add_nelist(h_mesh%mesh%faces%surface_mesh)
    
    ! now we loop over the extruded surface mesh, and use the 2 mappings
    ! to find the horizontal surface element above each extruded surface element
    do i=1, surface_element_count(out_mesh)
      ! find the element and copy its boundary id
      call find_h_surface_element(i)
    end do
        
    contains
    
    subroutine find_h_surface_element(sele)
    ! this is put in a separate function so we can have a face_loc
    ! array
    integer, intent(in):: sele
    
      integer, dimension(1: face_loc(out_mesh, sele)):: face_nodes, hsnodes
      integer, dimension(1):: surface_elements
      integer j, n
      
      ! the nodes of this face
      face_nodes=face_global_nodes(out_mesh, sele)
      ! find the nodes in the horizontal mesh using its surface node numbering (i.e.
      ! the node numbering of the edge around the horizontal mesh), or return 0
      ! if not below an edge node
      hsnodes=hnode2surfnode(out_mesh%mesh%columns(face_nodes))
      
      do j=1, size(hsnodes)
        if (hsnodes(j)==0) then
          ! this is probably a horizontal surface element, i.e. not on one
          ! of the vertical sides, we need to find out if it's on the top
          ! or at the bottom
          if (abs(node_val(out_mesh,out_mesh%dim,face_nodes(j)))<epsilon(1.0)) then
            ! it's on top (z=0, we might need some other logic in the future)
            out_mesh%mesh%faces%boundary_ids(sele)=top_boundary_id
            return
          else
            ! bottom then (I hope)
            out_mesh%mesh%faces%boundary_ids(sele)=bottom_boundary_id
            return
          end if
        end if
      end do
      
      ! this routine returns the (surface) elements that are adjacent 
      ! to all the given nodes
      assert(associated(h_mesh%mesh%faces%surface_mesh%adj_lists))
      assert(associated(h_mesh%mesh%faces%surface_mesh%adj_lists%nelist))
      call FindCommonElements(surface_elements, n, &
        h_mesh%mesh%faces%surface_mesh%adj_lists%nelist, hsnodes)
      ! this should return just 1 element
      if (n==0) then
        ! this happens for a top or bottom surface element that has all nodes
        ! on the boundary of the horizontal mesh (i.e. a stiff triangle in the
        ! corner of the horizontal mesh)
        if (all(abs(node_val(out_mesh, out_mesh%dim, face_nodes))<epsilon(1.0))) then
          ! it's on top (z=0, we might need some other logic in the future)
          out_mesh%mesh%faces%boundary_ids(sele)=top_boundary_id
            return
        else
          ! bottom then (I hope)
          out_mesh%mesh%faces%boundary_ids(sele)=bottom_boundary_id
          return
        end if
      else if (n>1) then
        ewrite(-1,*) "Something went wrong in the extrusion of the surface mesh"
        ewrite(-1,*) "For an extruded surface element multiple horizontal surface element are found. "
        FLAbort("Check your mesh extrusion settings.")
      end if
      
      out_mesh%mesh%faces%boundary_ids(sele)= &
        h_mesh%mesh%faces%boundary_ids(surface_elements(1))
      
    end subroutine find_h_surface_element
    
  end subroutine generate_extruded_surface_mesh

end module hadapt_metric_based_extrude
