#include "fdebug.h"

module edge_length_module

  use fldebug
  use vector_tools
  use unittest_tools
  use metric_tools
  use fields

  implicit none

  private
  
  public :: get_edge_lengths, get_directional_edge_lengths

  interface get_edge_lengths
    module procedure get_edge_lengths_field, get_edge_lengths_tensor, get_edge_lengths_ele
  end interface

  contains

  subroutine get_edge_lengths_field(metric, field)
    !!< Given a metric, calculate the desired edge length
    !!< at each node.
    type(tensor_field), intent(in) :: metric
    type(scalar_field), intent(inout) :: field
    integer :: i, j, dim
    real, dimension(mesh_dim(metric%mesh), mesh_dim(metric%mesh)) :: evectors
    real, dimension(mesh_dim(metric%mesh)) :: evalues, desired_lengths
    
    dim = mesh_dim(metric%mesh)

    do i=1,metric%mesh%nodes
      call eigendecomposition_symmetric(node_val(metric, i), evectors, evalues)
      desired_lengths = -1.0
      do j=1,dim
        if (evalues(j) /= 0.0) then
          desired_lengths(j) = 1.0/sqrt(abs(evalues(j))) ! compute the desired edge length in each direction
        end if
      end do
      field%val(i) = sum(desired_lengths, mask=(desired_lengths > 0.0)) / dim ! now take the average
    end do

    ewrite(2,*) "maxval(edgelengths) == ", maxval(field%val)
    ewrite(2,*) "maxloc(edgelengths) == ", maxloc(field%val)
    ewrite(2,*) "minval(edgelengths) == ", minval(field%val)
    ewrite(2,*) "minloc(edgelengths) == ", minloc(field%val)
  end subroutine get_edge_lengths_field
  
  subroutine get_edge_lengths_tensor(metric, tensor)
    !!< Replace the eigenvalues with the edge length corresponding to those
    !!< eigenvalues.
    type(tensor_field), intent(in) :: metric
    type(tensor_field), intent(inout) :: tensor
    real, dimension(mesh_dim(metric%mesh), mesh_dim(metric%mesh)) :: evectors
    real, dimension(mesh_dim(metric%mesh)) :: evalues

    integer :: i

    do i=1,metric%mesh%nodes
      call eigendecomposition_symmetric(node_val(metric, i), evectors, evalues)
      evalues = edge_length_from_eigenvalue(evalues)
      call eigenrecomposition(tensor%val(:, :, i), evectors, evalues)
    end do
  end subroutine get_edge_lengths_tensor

  subroutine get_edge_lengths_ele(mesh, coordinates, ele, edge_lengths)
    !!< Given an element, calculate the desired edge length at each node.
    type(mesh_type), intent(in) :: mesh
    type(vector_field), intent(in) :: coordinates
    integer, intent(in) :: ele
    real, dimension(mesh_dim(mesh), mesh_dim(mesh), ele_ngi(mesh, 1)),  intent(inout) :: edge_lengths
    real, dimension(mesh_dim(mesh), mesh_dim(mesh)) :: ele_tensor

    ele_tensor = simplex_tensor(coordinates, ele)
    edge_lengths = spread(edge_length_from_eigenvalue(ele_tensor), 3, size(edge_lengths, 3))

  end subroutine get_edge_lengths_ele

  subroutine get_directional_edge_lengths(metric, field, vec)
    !!< Get the edge lengths of metric
    !!< in the direction of the vector vec.

    type(tensor_field), intent(in) :: metric
    type(scalar_field), intent(inout) :: field
    real, dimension(metric%dim(1)), intent(in) :: vec
    real, dimension(metric%dim(1)) :: evals
    real, dimension(metric%dim(1), metric%dim(2)) :: evecs
    real :: len

    integer :: i, j, dim

    dim = metric%dim(1)

    do i=1,metric%mesh%nodes
      call eigendecomposition_symmetric(metric%val(:, :, i), evecs, evals)
      len = 0.0
      do j=1,dim
        len = len + abs(dot_product(vec, evecs(:, j))) * evals(j)
      end do
      field%val(i) = edge_length_from_eigenvalue(len)
    end do
  end subroutine get_directional_edge_lengths

end module edge_length_module
