#include "fdebug.h"

module aspect_ratios_module

  use metric_tools
  use vector_tools
  use fields
  implicit none

  contains

  subroutine get_aspect_ratios(metric, field)
    !!< Given a metric, calculate the desired edge length
    !!< at each node.
    type(tensor_field), intent(in) :: metric
    type(scalar_field), intent(inout) :: field
    integer :: i, dim
    real, dimension(mesh_dim(metric%mesh), mesh_dim(metric%mesh)) :: evectors
    real, dimension(mesh_dim(metric%mesh)) :: evalues
    
    dim = mesh_dim(metric%mesh)

    do i=1,metric%mesh%nodes
      call eigendecomposition_symmetric(node_val(metric, i), evectors, evalues)
      field%val(i) = aspect_ratio(abs(evalues))
    end do
  end subroutine get_aspect_ratios

end module aspect_ratios_module
