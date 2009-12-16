#include "fdebug.h"

module reference_meshes
  
  use conformity_measurement
  use fields
  use fldebug
  use global_parameters, only : OPTION_PATH_LEN
  use interpolation_module
  use merge_tensors
  use spud
  use state_module
  
  implicit none
  
  private
  
  public :: enforce_reference_meshes
    
contains

  subroutine enforce_reference_meshes(states, metric_positions, metric)
    type(state_type), dimension(:), intent(in) :: states
    type(vector_field), intent(in) :: metric_positions
    type(tensor_field), intent(inout) :: metric
    
    character(len = *), parameter :: base_name = "/mesh_adaptivity/hr_adaptivity/reference_mesh"
    character(len = OPTION_PATH_LEN) :: mesh_name
    integer :: i, nreference_meshes
    logical :: minimum
    type(mesh_type), pointer :: mesh => null()
    type(vector_field) :: reference_positions

    ewrite(1, *) "In enforce_reference_meshes"

    nreference_meshes = option_count(base_name)
    ewrite(2, *) "Number of reference meshes: ", nreference_meshes
    if(nreference_meshes == 0) return
    
    do i = 0, nreference_meshes - 1         
      call get_option(base_name // "[" // int2str(i) // "]" // "/mesh_name", mesh_name)
      ewrite(2, *) "Enforcing reference mesh: " // trim(mesh_name) 
      assert(size(states) > 0)
      mesh => extract_mesh(states(1), mesh_name)
      reference_positions = get_coordinate_field(states(1), mesh)
      
      minimum = have_option(base_name // "[" // int2str(i) // "]" // "/minimum")
#ifdef DDEBUG
      if(.not. minimum) then
        assert(have_option(base_name // "[" // int2str(i) // "]" // "/maximum"))
      end if
#endif
      ewrite(2, *) "Reference mesh is a minimum?", minimum
      
      call enforce_reference_mesh(metric_positions, reference_positions, metric, minimum)
      
      call deallocate(reference_positions)
    end do

    ewrite(1, *) "Exiting enforce_reference_meshes"
  
  end subroutine enforce_reference_meshes
  
  subroutine enforce_reference_mesh(metric_positions, reference_positions, metric, minimum)
    type(vector_field), intent(in) :: metric_positions
    type(vector_field), intent(in) :: reference_positions
    type(tensor_field), intent(inout) :: metric
    logical, intent(in) :: minimum
    
    type(tensor_field) :: mesh_metric, interpolated_mesh_metric
    
    call compute_mesh_metric(reference_positions, mesh_metric)
    call allocate(interpolated_mesh_metric, metric_positions%mesh, "MetricCoordinate")
    call linear_interpolation(mesh_metric, reference_positions, interpolated_mesh_metric, metric_positions)
    call merge_tensor_fields(metric, interpolated_mesh_metric, .not. minimum)
    
    call deallocate(mesh_metric)
    call deallocate(interpolated_mesh_metric)
    
  end subroutine enforce_reference_mesh

end module reference_meshes
