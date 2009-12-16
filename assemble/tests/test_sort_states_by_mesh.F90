#include "confdefs.h"

subroutine test_sort_states_by_mesh

  use fields
  use populate_state_module
  use spud
  use state_module
  use form_metric_field
  use metric_assemble
  use adapt_state_module
  use field_derivatives
  use vtk_interfaces
  use conservative_interpolation_module
  use global_parameters
  use interpolation_module
  use unittest_tools
  implicit none

  type(state_type), dimension(:), pointer :: states => null()
  type(state_type), dimension(:), allocatable :: sorted_states
  integer :: state, sfield
  type(scalar_field), pointer :: sfield_A, sfield_B
  logical :: fail

  call load_options("data/cg_interpolation_A.flml")
  call populate_state(states)

  call sort_states_by_mesh(states, sorted_states)

  fail = .false.
  do state=1,size(sorted_states)
    do sfield=2,scalar_field_count(sorted_states(state))
      sfield_A => extract_scalar_field(sorted_states(state), sfield-1)
      sfield_B => extract_scalar_field(sorted_states(state), sfield)
      fail = fail .or. (sfield_A%mesh%refcount%id /= sfield_B%mesh%refcount%id)
    end do
  end do

  call report_test("[sort_states_by_mesh]", fail, .false., "Ah, whatever")
  
  do state = 1, size(states)
    call deallocate(states(state))
  end do
  deallocate(states)
  
  do state = 1, size(sorted_states)
    call deallocate(sorted_states(state))
  end do
  deallocate(sorted_states)
  
  call report_test_no_references()

end subroutine test_sort_states_by_mesh
