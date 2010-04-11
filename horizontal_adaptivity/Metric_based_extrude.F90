#include "fdebug.h"

module hadapt_metric_based_extrude

  use elements
  use fields
  use sparse_tools
  use hadapt_advancing_front
  use spud
  use metric_tools
  use vector_tools
  use meshdiagnostics
  use halos
  implicit none

  public :: metric_based_extrude, combine_z_meshes, combine_r_meshes

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
      ele_shape(back_mesh, 1), back_mesh%name, trim(back_mesh%mesh%option_path))

    call deallocate(oned_shape)
    do column=1,node_count(h_mesh)
      call deallocate(back_z_meshes(column))
      call deallocate(out_z_meshes(column))
      call deallocate(back_sizing(column))
    end do
    
  end subroutine metric_based_extrude
    
  subroutine combine_z_meshes(h_mesh, z_meshes, out_mesh, &
    full_shape, mesh_name, option_path)
  !! Given the h_mesh and a z_mesh under each node of it combines these
  !! into a full horiz+vertic. mesh
  type(vector_field), intent(inout):: h_mesh
  type(vector_field), dimension(:), intent(in):: z_meshes
  type(vector_field), intent(out):: out_mesh
  type(element_type), intent(in):: full_shape
  !! the name of the topol. mesh to be created, not the coordinate field
  character(len=*), intent(in):: mesh_name, option_path
  
    type(csr_sparsity):: out_columns
    type(mesh_type):: mesh
    integer, dimension(:), allocatable:: no_hanging_nodes
    integer:: column, total_out_nodes, total_out_elements, z_elements, last_seen
    
    allocate(no_hanging_nodes(1:node_count(h_mesh)))
    no_hanging_nodes=0
    do column=1, size(z_meshes)
      if (node_owned(h_mesh, column)) then
        no_hanging_nodes(column)=ele_count(z_meshes(column))
      end if
    end do
    if (associated(h_mesh%mesh%halos)) then
      call halo_update(h_mesh%mesh%halos(2), no_hanging_nodes)
    end if
    
    total_out_nodes = 0
    ! For each column,
    ! add (number of nodes hanging off (== number of 1d elements)) * (element connectivity of chain)
    ! to compute the number of elements the extrusion routine will produce
    total_out_elements = 0
    do column=1, size(z_meshes)
      z_elements = no_hanging_nodes(column)
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

    out_mesh%mesh%option_path=option_path
    out_mesh%option_path=""
    out_mesh%mesh%periodic = mesh_periodic(h_mesh)

    last_seen = 0
    do column=1,node_count(h_mesh)
      if (node_owned(h_mesh, column)) then
        call append_to_structures(column, z_meshes(column), h_mesh, out_mesh, last_seen)
      else
        ! for non-owned columns we reserve node numbers, 
        ! but don't fill in out_mesh positions yet
        ! note that in this way out_mesh will obtain the same halo ordering
        ! convention as h_mesh
        out_mesh%mesh%columns(last_seen+1:last_seen+no_hanging_nodes(column)+1) = column
        last_seen = last_seen + no_hanging_nodes(column)+1
      end if
    end do
      
    call create_columns_sparsity(out_columns, out_mesh%mesh)
    
    if (associated(h_mesh%mesh%halos)) then
      ! derive l2 node halo for the out_mesh
      call derive_extruded_l2_node_halo(h_mesh%mesh, out_mesh%mesh, out_columns)
      ! positions in the non-owned columns can now simply be halo-updated
      call halo_update(out_mesh)
    end if
      
    call generate_layered_mesh(out_mesh, h_mesh)
    call deallocate(out_columns)
    
  end subroutine combine_z_meshes
    
  subroutine derive_extruded_l2_node_halo(h_mesh, out_mesh, columns)
    ! derive the l2 node halo for the extruded mesh
    type(mesh_type), intent(in):: h_mesh
    type(mesh_type), intent(inout):: out_mesh
    type(csr_sparsity), intent(in):: columns
    
    integer, dimension(:), allocatable:: nsends, nreceives, sends, receives
    integer:: proc_count, nowned_nodes
    integer:: i, j, l, node, proc
    
    assert(halo_count(h_mesh) == 2)
    assert(halo_valid_for_communication(h_mesh%halos(1)))
    assert(halo_valid_for_communication(h_mesh%halos(2)))
    
    assert(halo_count(out_mesh) == 0)
    assert(size(columns,1)==node_count(h_mesh))
    assert(size(columns,2)==node_count(out_mesh))
    
    ! this is easy, ownership of a node is determined by ownership
    ! of the column it is in, where columns correspond to nodes in h_mesh
    allocate(out_mesh%halos(2))    
    
    nowned_nodes = 0
    do i=1, node_count(h_mesh)
      if (node_owned(h_mesh%halos(2), i)) then
        nowned_nodes = nowned_nodes + row_length(columns, i)
      end if
    end do
    
    proc_count = halo_proc_count(h_mesh%halos(2))
    allocate(nsends(proc_count), nreceives(proc_count))
    do proc=1, proc_count

      nsends(proc) = 0
      do i=1, halo_send_count(h_mesh%halos(2), proc)
        node = halo_send(h_mesh%halos(2), proc, i)
        nsends(proc) = nsends(proc) + row_length(columns, node)
      end do
      
      nreceives(proc) = 0
      do i=1, halo_receive_count(h_mesh%halos(2), proc)
        node = halo_receive(h_mesh%halos(2), proc, i)
        nreceives(proc) = nreceives(proc) + row_length(columns, node)
      end do
    end do
    
    call allocate(out_mesh%halos(2), &
                  nsends = nsends, &
                  nreceives = nreceives, &
                  name = trim(out_mesh%name) // "Level2Halo", &
                  communicator = halo_communicator(h_mesh%halos(2)), &
                  nowned_nodes = nowned_nodes, &
                  data_type = halo_data_type(h_mesh%halos(2)), &
                  ordering_scheme = halo_ordering_scheme(h_mesh%halos(2)))
                  
    do proc=1, proc_count

      allocate( sends(1:nsends(proc)) )
      j=1
      do i=1, halo_send_count(h_mesh%halos(2), proc)
        node = halo_send(h_mesh%halos(2), proc, i)
        l = row_length(columns, node)
        sends(j:j+l-1) = row_m(columns, node)
        j=j+l
      end do
      assert( j==nsends(proc)+1 )
      call set_halo_sends(out_mesh%halos(2), proc, sends)
      deallocate(sends)
        
      allocate( receives(1:nreceives(proc)) )
      j=1
      do i=1, halo_receive_count(h_mesh%halos(2), proc)
        node = halo_receive(h_mesh%halos(2), proc, i)
        l = row_length(columns, node)
        receives(j:j+l-1) = row_m(columns, node)
        j=j+l
      end do
      assert( j==nreceives(proc)+1 )
      call set_halo_receives(out_mesh%halos(2), proc, receives)
      deallocate(receives)
      
    end do
    
    assert(halo_valid_for_communication(out_mesh%halos(2)))
    
    deallocate( nsends, nreceives )
    
    call create_global_to_universal_numbering(out_mesh%halos(2))
    call create_ownership(out_mesh%halos(2))    
    
  end subroutine derive_extruded_l2_node_halo
  
  subroutine combine_r_meshes(shell_mesh, r_meshes, out_mesh, &
    full_shape, mesh_name, option_path)
  !! Given the shell_mesh and a r_mesh under each node of it combines these
  !! into a full layered mesh
  type(vector_field), intent(inout):: shell_mesh
  type(vector_field), dimension(:), intent(in):: r_meshes
  type(vector_field), intent(out):: out_mesh
  type(element_type), intent(in):: full_shape
  !! the name of the topol. mesh to be created, not the coordinate field
  character(len=*), intent(in):: mesh_name, option_path
  
    type(csr_sparsity):: out_columns
    type(mesh_type):: mesh
    integer:: column, total_out_nodes, total_out_elements, r_elements, last_seen
    
    total_out_nodes = 0
    ! For each column,
    ! add (number of nodes hanging off (== number of 1d elements)) * (element connectivity of chain)
    ! to compute the number of elements the extrusion routine will produce
    total_out_elements = 0
    do column=1, size(r_meshes)
      r_elements = ele_count(r_meshes(column))
      total_out_nodes = total_out_nodes + r_elements + 1
      assert(associated(shell_mesh%mesh%adj_lists))
      assert(associated(shell_mesh%mesh%adj_lists%nelist))
      total_out_elements = total_out_elements + r_elements * row_length(shell_mesh%mesh%adj_lists%nelist, column)
    end do

    call allocate(mesh, total_out_nodes, total_out_elements, &
      full_shape, mesh_name)
    ! allocate mapping between extruded nodes to surface node (column number)
    ! it lies under
    allocate(mesh%columns(total_out_nodes))
    if (mesh_name=="CoordinateMesh") then
      call allocate(out_mesh, mesh_dim(shell_mesh)+1, mesh, "Coordinate")
    else
      call allocate(out_mesh, mesh_dim(shell_mesh)+1, mesh, trim(mesh_name)//"Coordinate")
    end if
    call deallocate(mesh)

    out_mesh%mesh%option_path=option_path
    out_mesh%option_path=""
    out_mesh%mesh%periodic = mesh_periodic(shell_mesh) ! Leave this here for now
    
    last_seen = 0
    do column=1,node_count(shell_mesh)
      call append_to_structures_radial(column, r_meshes(column), shell_mesh, out_mesh, last_seen)
    end do
      
    call create_columns_sparsity(out_columns, out_mesh%mesh)
    call generate_radially_layered_mesh(out_mesh, shell_mesh)
    call deallocate(out_columns)
    
  end subroutine combine_r_meshes

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

  subroutine append_to_structures_radial(column, r_mesh, shell_mesh, out_mesh, last_seen)
    integer, intent(in) :: column
    type(vector_field), intent(in) :: r_mesh, shell_mesh
    type(vector_field), intent(inout) :: out_mesh
    integer, intent(inout) :: last_seen

    integer :: j
    integer :: shell_dim, r_dim

    real, dimension(mesh_dim(out_mesh)) :: pos

    shell_dim = mesh_dim(shell_mesh)+1 ! This is the number of coordinates and not the topological dimesnion
    r_dim = mesh_dim(out_mesh)

    do j=1,node_count(r_mesh)
      last_seen = last_seen + 1
      out_mesh%mesh%columns(last_seen)=column
      pos(1:r_dim) = node_val(r_mesh, j)             ! note that in z-dir extrusion case only the z value was loaded into the out mesh from the z mesh. Here however,
      call set(out_mesh, last_seen, pos)
    end do
    
  end subroutine append_to_structures_radial

end module hadapt_metric_based_extrude
