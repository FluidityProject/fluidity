#include "fdebug.h"

module project_metric_to_surface_module

  use fields
  use sparse_tools
  use vector_tools
  use merge_tensors
  use quicksort
  use state_module
  use hadapt_advancing_front
  use boundary_conditions
  use field_options
  use edge_length_module
  use vtk_interfaces
  implicit none

  public :: project_metric_to_surface, vertically_align_metric, incorporate_bathymetric_metric
  contains

  subroutine project_metric_to_surface(volume_metric, h_mesh, surface_metric)
    !! Given a volume metric, project the metric up the columns to the surface mesh
    type(tensor_field), intent(in) :: volume_metric
    type(vector_field), intent(in) :: h_mesh
    type(tensor_field), intent(out) :: surface_metric

    type(csr_sparsity):: columns
    integer :: column
    real, dimension(mesh_dim(volume_metric), mesh_dim(volume_metric)) :: tmp_metric
    real, dimension(mesh_dim(volume_metric)-1, mesh_dim(volume_metric)-1) :: tmp_smetric

    call allocate(surface_metric, h_mesh%mesh, "Surface"//trim(volume_metric%name), field_type=volume_metric%field_type)
    
    if(volume_metric%field_type==FIELD_TYPE_CONSTANT) then
      tmp_metric = node_val(volume_metric, 1)
      tmp_smetric = reduce_metric_dimension(tmp_metric)
      call set(surface_metric, tmp_smetric)
    else  
      call create_columns_sparsity(columns, volume_metric%mesh)

      do column=1,node_count(h_mesh)
        call merge_up_columns(volume_metric, columns, column, tmp_metric)
        tmp_smetric = reduce_metric_dimension(tmp_metric)
        call set(surface_metric, column, tmp_smetric)
      end do
        
      call deallocate(columns)
    end if
    
    contains

      subroutine merge_up_columns(volume_metric, columns, column, tmp_metric)
      !! Descend down the column, merging as you go
        type(tensor_field), intent(in) :: volume_metric
        type(csr_sparsity), intent(in) :: columns
        integer, intent(in) :: column
        real, dimension(:, :), intent(out) :: tmp_metric

        integer, dimension(:), pointer :: column_nodes
        integer :: column_len
        integer :: node
        real, dimension(mesh_dim(volume_metric), mesh_dim(volume_metric)) :: tmp_val

        column_nodes => row_m_ptr(columns, column)
        column_len = row_length(columns, column)

        tmp_metric = node_val(volume_metric, column_nodes(1))
        do node=2,column_len
          tmp_val = node_val(volume_metric, column_nodes(node))
          call merge_tensor(tmp_metric, tmp_val)
        end do

      end subroutine merge_up_columns

  end subroutine project_metric_to_surface

  function reduce_metric_dimension(volume_metric, normal) result(surface_metric)
    !! Given a (e.g.) 3D metric, squash it down to a 2D metric
    real, dimension(:, :), intent(in) :: volume_metric
    real, dimension(:), intent(in), optional :: normal
    
    real, dimension(size(volume_metric, 1)-1, size(volume_metric, 2)-1) :: surface_metric

    integer :: vdim, sdim
    
    real, dimension(size(normal)) :: tangent1, tangent2, tangent1_dir
    real :: project
    
    real, dimension(size(volume_metric, 1), size(volume_metric, 2)) :: rotated_metric, &
                                                                       rotation_matrix
    assert(size(volume_metric, 1)==size(volume_metric, 2))
     
    vdim = size(volume_metric, 1)
    sdim = vdim - 1

    tangent1 = 0.0
    tangent2 = 0.0
    rotation_matrix = 0.0
    
    if(present(normal)) then
      if(vdim>1) then
      
        ! we want the 1st tangent to align with the x axis
        ! which we assume is the first index
        ! NOTE: this clearly won't work for normals that are in the x-direction
        tangent1_dir = 0.0
        tangent1_dir(1) = 1.0
        
        project = dot_product(normal, tangent1_dir)
        tangent1 = tangent1_dir - project*normal
        
        tangent1 = tangent1/sqrt(sum(tangent1**2))
        
        rotation_matrix(:,1) = tangent1
        
        if(vdim>2) then
        
          tangent2 = cross_product(normal, tangent1)
          ! we want the 2nd tangent to align with the y axis
          ! which we assume is the second index
          if(tangent2(2)<0.0) tangent2 = -tangent2
          
          rotation_matrix(:,2) = tangent2
        end if
      end if
      
      rotation_matrix(:,vdim) = normal
      
      rotated_metric = matmul(transpose(rotation_matrix), matmul(volume_metric, rotation_matrix))
      
    else
      rotated_metric = volume_metric
    end if

    surface_metric = rotated_metric(1:sdim, 1:sdim)

  end function reduce_metric_dimension
    
  subroutine vertically_align_metric(state, error_metric)
  !!< This separates out the metric projected to the horizontal plane (1 or 2D), 
  !!< and the metric projected in the vertical (gravity) direction, and combines 
  !!< this 1 or 2D hor. metric and 1D vert. metric back into a full resp. 2 or 3D metric. 
  !!< This ensures that for large aspect ratio problems the horizontal and vertical
  !!< metric are completely independent. Typically the metric for large aspect ratio
  !!< problems already decomposes in an (almost) vertical eigenvector and 2 horizontal
  !!< ones, however even the slightest tilt causes vertical error bounds to be "leaked"
  !!< into the horizontal leading to unexpected strict horizontal bounds.
  type(state_type), intent(in):: state
  type(tensor_field), intent(inout):: error_metric
    
    type(vector_field):: down_here
    type(vector_field), pointer:: down
    integer:: i, stat
      
    down => extract_vector_field(state, "GravityDirection", stat=stat)

    if (stat /= 0) return
    
    call allocate(down_here, down%dim, error_metric%mesh, "GravityDirectionOnCoordinateMesh")
    call remap_field(down, down_here)
    
    do i=1, node_count(error_metric)
      
      call set( error_metric, i, vertically_align_metric_node( &
                                    node_val(error_metric, i), &
                                    node_val(down_here, i) ) )
      
    end do
      
    call deallocate(down_here)
    
  end subroutine vertically_align_metric
    
  function vertically_align_metric_node(metric, down) result (projected_metric)
    real, dimension(:,:):: metric
    real, dimension(:):: down
    real, dimension(1:size(down), 1:size(down)):: projected_metric
    
    real, dimension(1:size(down), 1:size(down)):: pv, ph
    integer:: i
    
    ! vertical projection
    call outer_product( down, down, pv)
    
    ! ph = identity -pv
    ph=-pv
    forall(i=1:size(down)) ph(i,i)=ph(i,i)+1.0
    
    projected_metric=matmul( ph, matmul(metric, ph)) + matmul( pv, matmul(metric, pv))
    
  end function vertically_align_metric_node
  
  subroutine incorporate_bathymetric_metric(state, volume_metric, surface_positions, surface_metric)
  type(state_type), intent(in) :: state
  type(tensor_field), intent(in) :: volume_metric
  type(vector_field), intent(in) :: surface_positions ! only passed in for debugging output
  type(tensor_field), intent(inout) :: surface_metric
    
    type(scalar_field), pointer :: bottomdis
    integer, dimension(:), pointer :: surface_element_list
    type(mesh_type), pointer :: surface_mesh
    type(mesh_type) :: bottom_metric_mesh
    type(vector_field) :: normal

    type(vector_field) :: coordinate
    real, dimension(mesh_dim(volume_metric), face_ngi(volume_metric, 1)) :: normal_bdy
    real, dimension(face_ngi(volume_metric, 1)) :: detwei_bdy
    real, dimension(mesh_dim(volume_metric)) :: n
    
    integer, dimension(:), allocatable :: base_volume_index
    type(tensor_field) :: tmp_metric
    
    integer :: i, sele
    
    type(scalar_field) :: edge_lengths
    integer, save :: adaptcnt=0
    
    ewrite(1,*) 'Entering incorporate_bathymetric_metric'
    
    bottomdis => extract_scalar_field(state, "DistanceToBottom")
    call get_boundary_condition(bottomdis, 1, surface_mesh = surface_mesh, &
                                surface_element_list = surface_element_list)
    
    bottom_metric_mesh = make_mesh(surface_mesh, shape = face_shape(volume_metric, 1), &
                                   continuity = continuity(volume_metric))
    
    call allocate(normal, mesh_dim(volume_metric), bottom_metric_mesh, "BottomNormals")
    call zero(normal)
    
    call deallocate(bottom_metric_mesh)

    allocate(base_volume_index(node_count(normal)))
    base_volume_index = 0
    
    coordinate = get_coordinate_field(state, volume_metric%mesh)
    
    do i = 1, size(surface_element_list)
      sele = surface_element_list(i)
      
      call transform_facet_to_physical(coordinate, sele, &
                                       detwei_f=detwei_bdy, normal=normal_bdy)
                                       
      call addto(normal, ele_nodes(normal, i), &
                 shape_vector_rhs(ele_shape(normal, i), normal_bdy, detwei_bdy))
                 
      ! record the volume node numbers of each node in the bottom mesh
      ! (this will get overwritten by the same value if continuous, as we expect)
      base_volume_index(ele_nodes(normal, i)) = face_global_nodes(volume_metric, sele)
    end do
    
    assert(all(base_volume_index > 0))
  
    assert(associated(volume_metric%mesh%columns))

    call allocate(tmp_metric, surface_metric%mesh, "TemporarySurfaceMetric")
    call zero(tmp_metric)

    assert(node_count(tmp_metric)==node_count(normal))

    do i = 1, node_count(normal)
      ! get node normal
      n=node_val(normal,i)
      ! normalise it
      n=n/sqrt(sum(n**2))
      
      call set(tmp_metric, volume_metric%mesh%columns(base_volume_index(i)), &
               reduce_metric_dimension(node_val(volume_metric, base_volume_index(i)), n))
      
    end do

    if (have_option('/mesh_adaptivity/hr_adaptivity/debug/write_metric_stages')) then
      call allocate(edge_lengths, tmp_metric%mesh, "EdgeLengths")
      
      call get_edge_lengths(tmp_metric, edge_lengths)
      call vtk_write_fields('bathymetric_metric', adaptcnt, &
        surface_positions, surface_positions%mesh, &
        sfields=(/ edge_lengths /), tfields=(/ tmp_metric /) )

      call get_edge_lengths(surface_metric, edge_lengths)
      call vtk_write_fields('initial_horizontal_metric', adaptcnt, &
        surface_positions, surface_positions%mesh, &
        sfields=(/ edge_lengths /), tfields=(/ surface_metric /) )
        
      adaptcnt=adaptcnt+1
      
      call deallocate(edge_lengths)
    end if

    call merge_tensor_fields(surface_metric, tmp_metric)

    call deallocate(tmp_metric)
    call deallocate(normal)
    call deallocate(coordinate)

    ewrite(2,*) 'Leaving incorporate_bathymetric_metric'
  
  end subroutine incorporate_bathymetric_metric
  
end module project_metric_to_surface_module
