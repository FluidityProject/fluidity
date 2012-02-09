#include "confdefs.h"
#include "fdebug.h"

subroutine test_solenoidal_interpolation

  use fields
  use populate_state_module
  use spud
  use state_module
  use form_metric_field
  use metric_assemble
  use adapt_state_module
  use field_derivatives
  use vtk_interfaces
  use solenoidal_interpolation_module
  use global_parameters
  use interpolation_module
  use unittest_tools
  use reference_counting
  use diagnostic_fields_wrapper

  implicit none

  type(state_type), dimension(:), pointer :: states_old => null()
  type(state_type), dimension(:), pointer :: states_new => null()
  type(mesh_type), pointer :: u_mesh
  type(vector_field), pointer :: u_old, u_new, x_old, x_new
  type(scalar_field), pointer :: p_new, p_old
  type(scalar_field), pointer :: div_new, div_old
  type(state_type) :: interpolation_state_old, interpolation_state_new
  logical :: fail
  
  call test_solenoidal_interpolation_file_A_to_B(file_A   = "solenoidal_interpolation_A", &
                                                 file_B   = "solenoidal_interpolation_B", &
                                                 div_name = "FiniteElementDivergence")
  
  call test_solenoidal_interpolation_file_A_to_B(file_A   = "solenoidal_interpolation_press_cg_test_div_cv_A", &
                                                 file_B   = "solenoidal_interpolation_press_cg_test_div_cv_B", &
                                                 div_name = "ControlVolumeDivergence")
  
 contains
 
  subroutine test_solenoidal_interpolation_file_A_to_B(file_A, file_B, div_name)
  
     character(len=*), intent(in) :: file_A, file_B, div_name
     
     ewrite(1,*) "Testing solenoidal interpolation between flml files: "
     ewrite(1,*) "   1/ ",trim(file_A)//".flml"
     ewrite(1,*) "   2/ ",trim(file_B)//".flml"
     ewrite(1,*) "Using the divergence field: ",trim(div_name)     
     
     call load_options("data/"//trim(file_A)//".flml")
     call populate_state(states_old)
     call clear_options
     call load_options("data/"//trim(file_B)//".flml")
     call populate_state(states_new)

     u_mesh => extract_mesh(states_old(1), "VelocityMesh")
     u_old => extract_vector_field(states_old(1), "Velocity")
     x_old => extract_vector_field(states_old(1), "Coordinate")
     p_old => extract_scalar_field(states_old(1), "Pressure")
     div_old => extract_scalar_field(states_old(1), trim(div_name))

     u_new => extract_vector_field(states_new(1), "Velocity")
     x_new => extract_vector_field(states_new(1), "Coordinate")
     p_new => extract_scalar_field(states_new(1), "Pressure")
     div_new => extract_scalar_field(states_new(1), trim(div_name))

     call insert(interpolation_state_old, u_old, "Velocity")
     call insert(interpolation_state_old, u_old%mesh, "Mesh")
     call insert(interpolation_state_old, x_old, "Coordinate")
     call insert(interpolation_state_new, u_new, "Velocity")
     call insert(interpolation_state_new, u_new%mesh, "Mesh")
     call insert(interpolation_state_new, x_new, "Coordinate")
     call linear_interpolation(interpolation_state_old, interpolation_state_new)
     call deallocate(interpolation_state_old)
     call deallocate(interpolation_state_new)

     call calculate_diagnostic_variables(states_old)
     
     ewrite(1,*) "Initial max value of velocity divergence",maxval(div_old)
     
     call calculate_diagnostic_variables(states_new)

     ewrite(1,*) "After linear interpolation max value of velocity divergence",maxval(div_new)

     call vtk_write_state("data/"//trim(file_B), 0, state=states_old)
     call vtk_write_state("data/"//trim(file_B), 1, state=states_new)
    
     call solenoidal_interpolation(u_new, x_new, p_new%mesh, p_new)
  
     call calculate_diagnostic_variables(states_new)

     ewrite(1,*) "After solenoidal interpolation max value of velocity divergence",maxval(div_new)

     call vtk_write_state("data/"//trim(file_B), 2, state=states_new)

     fail = maxval(div_new) > epsilon(0.0_4)
     call report_test("[solenoidal interpolation: divergence free]", fail, .false., "Should be divergence-free!")

     call deallocate(states_old(1))
     call deallocate(states_new(1))
     call print_references(-1)
     call clear_options
   
   end subroutine test_solenoidal_interpolation_file_A_to_B
   
end subroutine test_solenoidal_interpolation
