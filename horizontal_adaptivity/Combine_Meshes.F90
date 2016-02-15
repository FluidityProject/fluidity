#include "fdebug.h"

module hadapt_combine_meshes
  !!< Extrude a given 2D mesh to a full 3D mesh.
  !!< The layer depths are specified by a sizing function
  !!< which can be arbitrary python.
  use elements
  use fields
  use spud
  use quadrature
  use global_parameters
  use sparse_tools
  use hadapt_advancing_front
  use linked_lists
  use halos
  implicit none

  private
  
  public :: combine_z_meshes

  contains

  subroutine combine_z_meshes(h_mesh, z_meshes, out_mesh, &
    full_shape, mesh_name, option_path, sl)
  !! Given the h_mesh and a z_mesh under each node of it combines these
  !! into a full horiz+vertic. mesh
  type(vector_field), intent(inout):: h_mesh
  type(vector_field), dimension(:), intent(in):: z_meshes
  type(vector_field), intent(out):: out_mesh
  type(element_type), intent(in):: full_shape
  !! the name of the topol. mesh to be created, not the coordinate field
  character(len=*), intent(in):: mesh_name, option_path
  logical, intent(in), optional :: sl
  logical :: sigma_layers
  
    type(csr_sparsity):: out_columns
    type(mesh_type):: mesh
    integer, dimension(:), allocatable:: no_hanging_nodes
    integer:: column, total_out_nodes, total_out_elements, z_elements, last_seen

    sigma_layers=present_and_true(sl)

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
    mesh%columns = 0
    ! allocate mapping between extruded elements to surface elements in horizontal mesh
    ! it lies under
    allocate(mesh%element_columns(total_out_elements))
    mesh%element_columns = 0
    ! if the horizontal mesh has region ids these can be propagated down into the full
    ! mesh.  allocate space for this...
    if(associated(h_mesh%mesh%region_ids)) then
      allocate(mesh%region_ids(total_out_elements))
    end if
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
        if (sigma_layers) then
          ! If we have sigma layers we will first use one chain to create a dummy
          ! 'mesh' that has the correct topological properties. We will later over write the
          !  vector field with the correct one. We always have chain '1', so we'll use that now.
          call append_to_structures(column, z_meshes(1), h_mesh, out_mesh, last_seen)
        else
          call append_to_structures(column, z_meshes(column), h_mesh, out_mesh, last_seen)
        end if
      else
        ! for non-owned columns we reserve node numbers, 
        ! but don't fill in out_mesh positions yet
        ! note that in this way out_mesh will obtain the same halo ordering
        ! convention as h_mesh
        out_mesh%mesh%columns(last_seen+1:last_seen+no_hanging_nodes(column)+1) = column
        last_seen = last_seen + no_hanging_nodes(column)+1
      end if
    end do
    assert(all(out_mesh%mesh%columns>0))
      
    call create_columns_sparsity(out_columns, out_mesh%mesh)
    
    if (associated(h_mesh%mesh%halos)) then
      ! derive l2 node halo for the out_mesh
      call derive_extruded_l2_node_halo(h_mesh%mesh, out_mesh%mesh, out_columns)
      ! positions in the non-owned columns can now simply be halo-updated
      call halo_update(out_mesh)
    end if
      
    call generate_layered_mesh(out_mesh, h_mesh)

    ! If we have sigma layers we now populate the vector field with the actual
    ! node positions.
    if (sigma_layers) then
      last_seen = 0
      do column=1,node_count(h_mesh)
        if (node_owned(h_mesh, column)) then
          call append_to_structures(column, z_meshes(column), h_mesh, out_mesh, last_seen)
        end if
      end do
      call halo_update(out_mesh)
    end if

    call deallocate(out_columns)
    
  end subroutine combine_z_meshes
  
  subroutine append_to_structures(column, z_mesh, h_mesh, out_mesh, last_seen)
    integer, intent(in) :: column
    type(vector_field), intent(in) :: z_mesh, h_mesh
    type(vector_field), intent(inout) :: out_mesh
    integer, intent(inout) :: last_seen

    integer :: j
    integer :: v_dim
    logical :: radial_extrusion

    real, dimension(mesh_dim(out_mesh)) :: pos, origin, direction

    radial_extrusion = have_option("/geometry/spherical_earth")
    
    v_dim = mesh_dim(out_mesh)

    origin = 0.0
    origin(1:h_mesh%dim) = node_val(h_mesh, column)
    direction = 0.0
    if (radial_extrusion) direction = origin/norm2(origin)

    do j=1,node_count(z_mesh)
      last_seen = last_seen + 1
      out_mesh%mesh%columns(last_seen)=column
      if (radial_extrusion) then
        pos = origin + direction*node_val(z_mesh, 1, j)
      else
        pos = origin
        pos(v_dim) = node_val(z_mesh, 1, j)
      end if
      call set(out_mesh, last_seen, pos)
    end do
    
  end subroutine append_to_structures

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


end module hadapt_combine_meshes
