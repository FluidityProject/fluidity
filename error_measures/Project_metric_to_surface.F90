#include "fdebug.h"

module project_metric_to_surface_module

  use fields
  use sparse_tools
  use vector_tools
  use merge_tensors
  use quicksort
  use state_module
  use hadapt_advancing_front
  implicit none

  public :: project_metric_to_surface, vertically_align_metric
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

    call allocate(surface_metric, h_mesh%mesh, "SurfaceMetric")
    call create_columns_sparsity(columns, volume_metric%mesh)

    do column=1,node_count(h_mesh)
      call merge_up_columns(volume_metric, columns, column, tmp_metric)
      tmp_smetric = reduce_metric_dimension(tmp_metric)
      call set(surface_metric, column, tmp_smetric)
    end do
      
    call deallocate(columns)

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

  function reduce_metric_dimension(volume_metric) result(surface_metric)
    !! Given a (e.g.) 3D metric, squash it down to a 2D metric
    !! by taking out the two most horizontal eigen components
    real, dimension(:, :), intent(in) :: volume_metric
    real, dimension(size(volume_metric, 1), size(volume_metric, 1)) :: volume_evecs
    real, dimension(size(volume_metric, 1)-1, size(volume_metric, 1)-1) :: surface_metric, surface_evecs
    real, dimension(size(volume_metric, 1)) :: volume_evals, norms
    real, dimension(size(volume_metric, 1)-1) :: surface_evals
    integer, dimension(size(volume_metric, 1)) :: sorted_idx

    integer :: vdim, sdim
    integer :: i

    vdim = size(volume_metric, 1)
    sdim = vdim - 1

    call eigendecomposition_symmetric(volume_metric, volume_evecs, volume_evals)

    ! We loop over the eigenvectors of the metric and project them to the surface.
    ! We then discard the one with the smallest norm -- this is the projection of the vertical
    ! component.
    ! We record the infinity norm of the first sdim components of the vector.
    do i=1,vdim
      norms(i) = maxval(abs(volume_evecs(1:sdim, i)))
    end do

    call qsort(norms, sorted_idx)

    do i=1,sdim
      surface_evecs(:, i) = volume_evecs(1:sdim, sorted_idx(vdim - i + 1))
      assert(norm2(surface_evecs(:, i)) > 0.0)
      surface_evecs(:, i) = surface_evecs(:, i) / norm2(surface_evecs(:, i)) ! and normalise
      surface_evals(i) = volume_evals(sorted_idx(vdim - i + 1))
    end do

    call eigenrecomposition(surface_metric, surface_evecs, surface_evals)
    
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
  
end module project_metric_to_surface_module
