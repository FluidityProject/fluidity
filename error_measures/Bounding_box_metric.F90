#include "fdebug.h"

module bounding_box_metric
!!< This module implements code to restrict
!!< the desired edge lengths of the metric
!!< to the dimensions of a bounding box
!!< for the domain. Adaptivity really, really
!!< doesn't like it, otherwise.

  use vtk_interfaces
  use edge_length_module
  use merge_tensors
  use unittest_tools
  use metric_tools
  use fields
  use spud
  implicit none

  logical, save :: use_bounding_box_metric = .true.
  logical, save :: bounding_box_initialised = .false.
  real :: width_factor = 2.0
  real, dimension(:), allocatable, save :: domain_width

  contains

  subroutine initialise_bounding_box_metric(positions)
    type(vector_field), intent(in) :: positions
    integer :: i, dim

    if (.not. bounding_box_initialised) then
      use_bounding_box_metric = .true.
      if (have_option("/mesh_adaptivity/hr_adaptivity/bounding_box_factor")) then
        call get_option("/mesh_adaptivity/hr_adaptivity/bounding_box_factor", width_factor)
      else
        width_factor = 2.0
      end if
      bounding_box_initialised = .true.
      dim = positions%dim
      allocate(domain_width(dim))
      do i=1,dim
        ! Get the width of the bounding box.
        domain_width(i) = abs(maxval(positions%val(i)%ptr) - minval(positions%val(i)%ptr))
      end do
    end if

  end subroutine initialise_bounding_box_metric

  subroutine form_bounding_box_metric(positions, error_metric, max_metric)
    type(tensor_field), intent(inout) :: error_metric !!< The metric formed so far
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(in) :: max_metric

    real, dimension(positions%dim, positions%dim) :: domain_metric
    integer :: i
    integer, save :: adaptcnt = 0
    real, dimension(positions%dim) :: tmp_domain_width

    type(scalar_field) :: edgelen
    logical :: debug_metric
    real, dimension(positions%dim, positions%dim) :: max_metric_nodes

    ewrite(2,*) "++: Constraining metric to bounding box"

    debug_metric = have_option("/mesh_adaptivity/hr_adaptivity/debug/write_metric_stages")

    tmp_domain_width = width_factor * domain_width
    tmp_domain_width = eigenvalue_from_edge_length(tmp_domain_width)

    ! Now make the diagonal matrix out of it and merge.
    do i=1,error_metric%mesh%nodes
      domain_metric = get_mat_diag(tmp_domain_width) ! domain_metric might change in merge_tensor
      call merge_tensor(error_metric%val(:, :, i), domain_metric)
      max_metric_nodes = node_val(max_metric, i)
      call merge_tensor(error_metric%val(:, :, i), max_metric_nodes)
    end do

    if (debug_metric) then
      call allocate(edgelen, error_metric%mesh, "Desired edge lengths")
      call get_edge_lengths(error_metric, edgelen)
      call vtk_write_fields("bounding_box", adaptcnt, positions, positions%mesh, &
                            sfields=(/edgelen/), tfields=(/error_metric/))
      call deallocate(edgelen)
    endif

    adaptcnt = adaptcnt + 1
  end subroutine form_bounding_box_metric

end module bounding_box_metric
