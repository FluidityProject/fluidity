#include "fdebug.h"

module aspect_ratios_module

  use spud
  use vector_tools
  use metric_tools
  use fields

  implicit none

  private

  public :: bound_metric_aspect_ratio, get_aspect_ratios
  
  interface bound_metric_aspect_ratio
    module procedure bound_metric_aspect_ratio_options, &
      & bound_metric_aspect_ratio_ratio
  end interface bound_metric_aspect_ratio

contains

  subroutine bound_metric_aspect_ratio_options(metric)
    !!< Apply a metric aspect ratio bound

    type(tensor_field), intent(inout) :: metric

    integer :: stat
    real :: aspect_ratio_bound

    call get_option("/mesh_adaptivity/hr_adaptivity/aspect_ratio_bound", aspect_ratio_bound, stat = stat)
    if(stat /= SPUD_NO_ERROR) then
      ewrite(1, *) "No aspect ratio bound"
      return
    end if

    if(aspect_ratio_bound <= 0.0) then
      FLExit("Aspect ratio bound must be positive")
    end if
    
    call bound_metric_aspect_ratio(metric, aspect_ratio_bound)

  end subroutine bound_metric_aspect_ratio_options
  
  subroutine bound_metric_aspect_ratio_ratio(metric, aspect_ratio_bound)
    !!< Apply a metric aspect ratio bound

    type(tensor_field), intent(inout) :: metric
    real, intent(in) :: aspect_ratio_bound

    integer :: i
    real :: evals_ratio_bound
    real, dimension(metric%dim(1)) :: evals
    real, dimension(metric%dim(1), metric%dim(2)) :: evecs

    ewrite(1, *) "In bound_metric_aspect_ratio_ratio"

    ewrite(2, *) "Aspect ratio bound: ", aspect_ratio_bound
    assert(aspect_ratio_bound > 0.0)
    
    if(aspect_ratio_bound < 1.0) then
      evals_ratio_bound = aspect_ratio_bound ** 2
    else
      evals_ratio_bound = (1.0 / aspect_ratio_bound) ** 2
    end if
    ewrite(2, *) "Eigenvalues ratio bound: ", evals_ratio_bound
    
    do i = 1, node_count(metric)
      call eigendecomposition_symmetric(metric%val(:, :, i), evecs, evals)
      evals = max(evals, evals_ratio_bound * maxval(evals))
      call eigenrecomposition(metric%val(:, :, i), evecs, evals)
    end do
    
    ewrite(1, *) "Exiting bound_metric_aspect_ratio_ratio"
    
  end subroutine bound_metric_aspect_ratio_ratio

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
