#include "fdebug.h"

module richardson_metric_module

  use fldebug
  use global_parameters, only : OPTION_PATH_LEN
  use spud
  use vector_tools
  use unittest_tools
  use metric_tools
  use fields
  use edge_length_module
  use state_module
  use field_options
  use limit_metric_module
  use merge_tensors
  
  implicit none

  private
  
  public :: initialise_richardson_number_metric, form_richardson_number_metric
  
  logical, public :: use_richardson_number_metric = .false.

contains

  subroutine initialise_richardson_number_metric
  
    use_richardson_number_metric = option_count("/material_phase/scalar_field::RichardsonNumber" // &
      & "/diagnostic/adaptivity_options/richardson_number_metric") > 0
      
  end subroutine initialise_richardson_number_metric

  subroutine form_richardson_number_metric(states, metric)
    type(state_type), dimension(:), intent(in) :: states
    type(tensor_field), intent(inout) :: metric

    character(len = OPTION_PATH_LEN) :: base_path
    integer :: i, stat
    real :: ri_min, ri_max, h_min, h_max
    type(scalar_field), pointer :: ri
    type(scalar_field) :: ri_metric
    type(tensor_field) :: ri_metric_tensor
    
    do i = 1, size(states)
      ri => extract_scalar_field(states(i), "RichardsonNumber", stat)
      if(stat /= SPUD_NO_ERROR) cycle
      base_path = trim(complete_field_path(ri%option_path)) // "/adaptivity_options/richardson_number_metric"
      if(.not. have_option(base_path)) cycle
      
      call get_option(trim(base_path) // "/min_ri", ri_min, default = 0.0)
      call get_option(trim(base_path) // "/max_ri", ri_max)
      call get_option(trim(base_path) // "/min_edge_length", h_min)
      call get_option(trim(base_path) // "/max_edge_length", h_max)
            
      call form_richardson_number_metric_internal(states(i), metric, ri_metric, ri_min, ri_max, h_min, h_max)
      
#ifdef DDEBUG
      call form_anisotropic_metric_from_isotropic_metric(ri_metric, ri_metric_tensor)
      call check_metric(ri_metric_tensor)
      call deallocate(ri_metric_tensor)
#endif
      
      if(have_option(trim(base_path) // "/anisotropy_preserving_merge")) then
        ewrite(2, *) "Using anisotropy preserving metric merge"
        call merge_isotropic_anisotropic_metrics(ri_metric, metric)
      else
        ewrite(2, *) "Using direct metric merge"
        call form_anisotropic_metric_from_isotropic_metric(ri_metric, ri_metric_tensor)
        call merge_tensor_fields(metric, ri_metric_tensor)
        call deallocate(ri_metric_tensor)
      end if
      call deallocate(ri_metric)
    end do
    
#ifdef DDEBUG
    call check_metric(metric)
#endif

  end subroutine form_richardson_number_metric

  subroutine form_richardson_number_metric_internal(state, metric, ri_metric, ri_min, ri_max, h_min, h_max)
  !! Form a length scale by a linear scaling with RichardsonNumber.
  !! If ri <= ri_min, then let h = h_min
  !! If ri >= ri_max, then let h = h_max
  !! Else, do a linear fit.
    type(state_type), intent(in) :: state
    type(tensor_field), intent(inout) :: metric
    type(scalar_field), intent(out) :: ri_metric
    real, intent(in) :: ri_min, ri_max, h_min, h_max

    type(scalar_field), pointer :: richardson_number
    integer :: node
    real :: eigenval, ri, h, m

    ewrite(1, *) "In form_richardson_number_metric_internal"

    call allocate(ri_metric, metric%mesh, "RichardsonNumberMetric")
    call zero(ri_metric)

    richardson_number => extract_scalar_field(state, "RichardsonNumber")

    do node=1,node_count(metric)
      ri = node_val(richardson_number, node)
      if (ri <= ri_min) then
        h = h_min
      else if (ri >= ri_max) then
        h = h_max
      else
        m = (h_max - h_min) / (ri_max - ri_min)
        h = m * (ri - ri_max) + h_max
      end if
      if(is_nan(h)) then
        ! Oops, seems we have a NaN in our Richardson number - let's assume we
        ! need a maximum edge length here
        h = h_max
      end if
#ifdef DDEBUG
      if(h <= 0.0) then
        ewrite(-1, *) "Edge length: ", h
        FLAbort("Negative edge length")
      end if
#endif

      eigenval = eigenvalue_from_edge_length(h)
      call set(ri_metric, node, eigenval)
    end do
        
    ewrite(1, *) "Exiting form_richardson_number_metric_internal"
    
  end subroutine form_richardson_number_metric_internal

  subroutine merge_isotropic_anisotropic_metrics(isotropic, anisotropic)
    type(scalar_field), intent(in) :: isotropic
    type(tensor_field), intent(inout) :: anisotropic

    integer :: i
    real :: isotropic_edge_length
    real, dimension(mesh_dim(anisotropic)) :: anisotropic_edge_lengths, eigenvals
    real, dimension(mesh_dim(anisotropic), mesh_dim(anisotropic)) :: eigenvecs
    type(scalar_field) :: edge_lengths

    assert(isotropic%mesh == anisotropic%mesh)
    
    call allocate(edge_lengths, anisotropic%mesh, "EdgeLengths")
    call get_edge_lengths(anisotropic, edge_lengths)

    do i = 1, node_count(isotropic)
      call eigendecomposition_symmetric(node_val(anisotropic, i), eigenvecs, eigenvals)
      anisotropic_edge_lengths = edge_length_from_eigenvalue(eigenvals)
      isotropic_edge_length = edge_length_from_eigenvalue(node_val(isotropic, i))
      if(node_val(edge_lengths, i) > isotropic_edge_length) then
        anisotropic_edge_lengths = anisotropic_edge_lengths * (isotropic_edge_length / node_val(edge_lengths, i))
        eigenvals = eigenvalue_from_edge_length(anisotropic_edge_lengths)
        call eigenrecomposition(anisotropic%val(:, :, i), eigenvecs, eigenvals)
      end if
    end do
    
    call deallocate(edge_lengths)
    
  end subroutine merge_isotropic_anisotropic_metrics
  
end module richardson_metric_module
