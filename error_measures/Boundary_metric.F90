#include "fdebug.h"

module boundary_metric

  use global_parameters, only: domain_bbox
  use unittest_tools, only: get_mat_diag
  use fields
  use node_boundary, only: initialise_boundcount, node_lies_on_boundary

  implicit none

  logical, save :: use_boundary_metric = .false.

  contains

  subroutine initialise_boundary_metric
    use_boundary_metric = .false.
  end subroutine

  subroutine form_boundary_metric(error_metric, positions)
    type(tensor_field), intent(inout) :: error_metric
    type(vector_field), intent(in) :: positions
    integer :: i, node
    real, dimension(positions%dim) :: domain_width

    call initialise_boundcount(error_metric%mesh, positions)

    do i=1,positions%dim
      domain_width(i) = abs(domain_bbox(i,2)-domain_bbox(i,1))
    end do

    do node=1,node_count(error_metric)
      if (node_lies_on_boundary(node)) then
        error_metric%val(:, :, node) = get_mat_diag(domain_width)
      end if
    end do
  end subroutine form_boundary_metric

end module boundary_metric
