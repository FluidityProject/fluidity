#include "fdebug.h"

module limit_metric_module

  use fldebug
  use vector_tools, only : determinant => det
  use elements
  use spud
  use fields
  use meshdiagnostics

  implicit none

  private

  public :: limit_metric, limit_metric_elements, expected_elements, &
    & expected_nodes, determinant

  interface expected_nodes
    module procedure expected_nodes_expected_elements, expected_nodes_metric
  end interface expected_nodes
  
  interface limit_metric
    module procedure limit_metric_nodes_options, limit_metric_nodes_minmax, &
      & limit_metric_nodes_target
  end interface limit_metric
  
  interface limit_metric_elements
    module procedure limit_metric_elements_minmax, limit_metric_elements_target
  end interface limit_metric_elements

contains

  subroutine limit_metric_nodes_options(positions, metric)
    type(tensor_field), intent(inout) :: metric
    type(vector_field), intent(in) :: positions

    character(len = *), parameter :: base_path = "/mesh_adaptivity/hr_adaptivity"
    integer :: max_nodes, min_nodes, nodes, stat
    real :: increase_tolerance
    
    call mesh_stats(positions, nodes = nodes)

    call get_option(base_path // "/minimum_number_of_nodes", min_nodes, default = 1)
    if(have_option(base_path // "/minimum_number_of_nodes/per_process")) then
      min_nodes = min_nodes * getnprocs()
    end if
    call get_option(base_path // "/maximum_number_of_nodes", max_nodes, default = 100000)
    if(have_option(base_path // "/maximum_number_of_nodes/per_process")) then
      max_nodes = max_nodes * getnprocs()
    end if
    call get_option(base_path // "/max_node_increase", increase_tolerance, stat = stat)
    if(stat == SPUD_NO_ERROR) then
      max_nodes = min(max_nodes, int(nodes * increase_tolerance))
    end if

    call limit_metric(positions, metric, min_nodes, max_nodes)

  end subroutine limit_metric_nodes_options
  
  subroutine limit_metric_nodes_minmax(positions, metric, min_nodes, max_nodes)
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(inout) :: metric
    integer, intent(in) :: min_nodes
    integer, intent(in) :: max_nodes
    
    integer :: elements, max_eles, min_eles, nodes
    ! The ratio of elements to nodes
    real :: eles_per_node
    
    assert(min_nodes > 0)
    assert(max_nodes >= min_nodes)
    
    call mesh_stats(positions, nodes = nodes, elements = elements)
    
    ! FIXME: maybe a better way to do this?
    eles_per_node = float(elements) / float(nodes)
    
    min_eles = eles_per_node * min_nodes
    max_eles = eles_per_node * max_nodes

    if(min_eles < 0) then
      ewrite(-1, *) "Minimum elements: ", min_eles
      FLAbort("Invalid minimum number of elements")
    end if
    if(max_eles < 0) then
      ewrite(-1, *) "Maximum elements: ", max_eles
      FLAbort("Invalid maximum number of elements")
    end if

    call limit_metric_elements(positions, metric, min_eles, max_eles)
    
  end subroutine limit_metric_nodes_minmax
  
  subroutine limit_metric_elements_minmax(positions, metric, min_eles, max_eles)
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(inout) :: metric
    integer, intent(in) :: min_eles
    integer, intent(in) :: max_eles
    
    integer :: expected_eles
    ! the scaling factor to divide the metric by
    real :: beta
    
    expected_eles = expected_elements(positions, metric, global = .true.)
    
    if(expected_eles > max_eles) then
      beta = ((1.0 / expected_eles) * max_eles) ** (2.0 / positions%dim)
      ewrite(2,*) "Scaling factor to conform to maximum node limit: ", beta
      call scale(metric, beta)
    else if(expected_eles < min_eles) then
      beta = ((1.0 / expected_eles) * min_eles) ** (2.0 / positions%dim)
      ewrite(2,*) "Scaling factor to conform to minimum node limit: ", beta
      call scale(metric, beta)
    end if
    
  end subroutine limit_metric_elements_minmax
  
  subroutine limit_metric_nodes_target(positions, metric, target_nodes)
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(inout) :: metric
    integer, intent(in) :: target_nodes
    
    call limit_metric(positions, metric, min_nodes = target_nodes, max_nodes = target_nodes)
    
  end subroutine limit_metric_nodes_target
  
  subroutine limit_metric_elements_target(positions, metric, target_eles)
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(inout) :: metric
    integer, intent(in) :: target_eles
    
    call limit_metric_elements(positions, metric, min_eles = target_eles, max_eles = target_eles)
    
  end subroutine limit_metric_elements_target

  function expected_elements(old_positions, metric, global) result(xpct)
    type(vector_field), intent(in) :: old_positions
    type(tensor_field), intent(in) :: metric
    !! If present and .true., calculate the global number of expected elements
    logical, optional, intent(in) :: global
    integer :: ele
    real, dimension(mesh_dim(metric), mesh_dim(metric)) :: avg_metric
    real :: sumvol, det
    integer :: xpct
    real :: gamma
    
    logical :: lglobal
    
    lglobal = present_and_true(global)

    sumvol = 0.0

    ! Gamma is the volume of an optimal element
    ! (in metric space)
    select case(mesh_dim(metric))
      case(3)
        gamma = 1.0 / sqrt(72.0)
      case(2)
        gamma = sqrt(3.0) / 4.0
      case(1)
        gamma = 1.0
      case default
        FLAbort("Invalid dimension")
    end select

    do ele=1,ele_count(old_positions)
      if(lglobal) then
        if(.not. element_owned(old_positions, ele)) cycle
      end if
      avg_metric = sum(ele_val(metric, ele), 3) / ele_loc(metric, ele)
      det = determinant(avg_metric)
      sumvol = sumvol + abs(sqrt(det) * simplex_volume(old_positions, ele))
    end do
    if(lglobal) call allsum(sumvol)

    if(sumvol < 0.0) then
      ewrite(-1, *) "Total volume in metric space: ", sumvol
      FLAbort("Negative volume")
    end if
    if((sumvol/gamma)>=huge(xpct)) then
      ewrite(-1, *) "ERROR: The error metric &
        & indicates that number of elements required is ", sumvol/gamma
      ewrite(-1, *) "If this is what you want then, great, congratulations, &
        &this is a record. Please ask the developers to get rid of this &
        &32bit integer that's causing trouble. Otherwise, please review your &
        &error targets"
        FLExit("integer overflow")
    end if
    xpct = int(sumvol / gamma)
    if (xpct == 0) xpct = 1

    if(xpct < 0) then
      ewrite(-1, *) "Expected elements: ", xpct
      FLAbort("Invalid number of expected elements")
    end if
    if (lglobal) then
      ewrite(2, *) "Expected global n/o elements: ", xpct
    else
      ewrite(2, *) "Expected n/o elements: ", xpct
    end if

  end function expected_elements

  function expected_nodes_expected_elements(old_positions, expected_eles, global) result(expected_nods)
    !!< Return the expected number of nodes based upon the supplied expected
    !!< number of elements

    type(vector_field), intent(in) :: old_positions
    integer, intent(in) :: expected_eles
    !! If present and .true., calculate the global number of expected elements
    logical, optional, intent(in) :: global

    integer :: elements, expected_nods, nodes

    if(present_and_true(global)) then
      call mesh_stats(old_positions, nodes = nodes, elements = elements)
    else
      nodes = node_count(old_positions)
      elements = ele_count(old_positions)
    end if
      
    ! FIXME: maybe a better way to do this?
    expected_nods = (float(nodes) / float(elements)) * expected_eles

  end function expected_nodes_expected_elements

  function expected_nodes_metric(old_positions, metric, global) result(expected_nods)
    !!< Return the expected number of nodes based upon the supplied metric

    type(vector_field), intent(in) :: old_positions
    type(tensor_field), intent(in) :: metric
    !! If present and .true., calculate the global number of expected elements
    logical, optional, intent(in) :: global

    integer :: expected_nods

    integer :: expected_eles

    expected_eles = expected_elements(old_positions, metric, global)
    expected_nods = expected_nodes(old_positions, expected_eles, global)

  end function expected_nodes_metric
  
  subroutine limit_metric_module_check_options()
    !!< Check metric limiting specific options
    
    character(len = *), parameter :: base_path = "/mesh_adaptivity/hr_adaptivity"
    integer :: max_nodes, min_nodes, stat
   
    if(.not. have_option(base_path)) then
      ! Nothing to check
      return
    end if
   
    call get_option(base_path // "/minimum_number_of_nodes", min_nodes, default = 0)
    call get_option(base_path // "/maximum_number_of_nodes", max_nodes, stat = stat)
    if(stat /= SPUD_NO_ERROR) then
      FLExit("Maximum number of nodes must be specified when using hr mesh adaptivity")
    else if(min_nodes > max_nodes) then
      FLExit("The minimum number of nodes cannot be greater than the maximum number of nodes")
    end if
    
  end subroutine limit_metric_module_check_options

end module limit_metric_module
