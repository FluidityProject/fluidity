#include "fdebug.h"

module actuator_line_metric

  use spud 
  use fldebug
  use unittest_tools
  use metric_tools
  use fields
  use vtk_interfaces
  use edge_length_module
  use merge_tensors
  use actuator_line_model

  implicit none
  
  private

  public :: initialise_actuator_line_metric, form_actuator_line_metric, & 
            use_actuator_line_metric, actuator_line_initialised
            
  logical, save :: use_actuator_line_metric = .true.
  logical, save :: actuator_line_initialised = .false.
  real :: alm_factor = 5.0
  
  contains
 
  subroutine initialise_actuator_line_metric

      if(.not. actuator_line_initialised) then
          use_actuator_line_metric= .true.
        if (have_option("/mesh_adaptivity/hr_adaptivity/actuator_line_factor")) then
            call get_option("mesh_adapticity/hr_adaptivity/actuator_line_factor",alm_factor)
        else
            alm_factor=5.0
        endif
        actuator_line_initialised =.true.  
    endif
  end subroutine initialise_actuator_line_metric
  
  subroutine form_actuator_line_metric(positions, error_metric, max_metric)
    type(tensor_field), intent(inout) :: error_metric !!< The metric formed so far
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(in) :: max_metric

    real, dimension(positions%dim, positions%dim) :: domain_metric
    integer :: i
    integer, save :: adaptcnt = 0
    real, dimension(positions%dim) :: domain_width

    type(scalar_field) :: edgelen
    logical :: debug_metric
    real, dimension(positions%dim, positions%dim) :: max_metric_nodes

    ewrite(2,*) "++: Constraining metric to actuator line model"

    debug_metric = have_option("/mesh_adaptivity/hr_adaptivity/debug/write_metric_stages")

  !  do i = 1, positions%dim
  !    domain_width(i) = abs(domain_bbox(i,2) - domain_bbox(i,1))
  !  end do
  !  domain_width = width_factor * domain_width
  !  domain_width = eigenvalue_from_edge_length(domain_width)

  !  ! Now make the diagonal matrix out of it and merge.
  !  do i=1,error_metric%mesh%nodes
  !    domain_metric = get_mat_diag(domain_width) ! domain_metric might change in merge_tensor
  !    call merge_tensor(error_metric%val(:, :, i), domain_metric)
  !    max_metric_nodes = node_val(max_metric, i)
  !    call merge_tensor(error_metric%val(:, :, i), max_metric_nodes)
  !  end do

    if (debug_metric) then
      call allocate(edgelen, error_metric%mesh, "Desired edge lengths")
      call get_edge_lengths(error_metric, edgelen)
      call vtk_write_fields("actuator_line", adaptcnt, positions, positions%mesh, &
                            sfields=(/edgelen/), tfields=(/error_metric, max_metric/))
      call deallocate(edgelen)
    endif

    adaptcnt = adaptcnt + 1
  end subroutine form_actuator_line_metric

end module actuator_line_metric
