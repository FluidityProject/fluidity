! Offline vtu diagnostics tools
!
! James Maddison

#include "fdebug.h"

module fldiagnostics_module

  use fields_data_types
  use fields
  use fldebug
  use state_module

  implicit none

  private

  public :: find_existing_field_rank, find_mesh_field

contains

  function find_mesh_field(state, meshfield_name) result(mesh)
    !!< Return the mesh for the field with name meshfield_name in state
  
    type(state_type), intent(in) :: state
    character(len = *), intent(in) :: meshfield_name
  
    type(mesh_type), pointer :: mesh
    
    integer :: s_stat, v_stat, t_stat
    type(scalar_field), pointer :: s_field
    type(tensor_field), pointer :: t_field
    type(vector_field), pointer :: v_field
    
    s_field => extract_scalar_field(state, meshfield_name, s_stat)
    v_field => extract_vector_field(state, meshfield_name, v_stat)
    t_field => extract_tensor_field(state, meshfield_name, t_stat)
  
    if((s_stat == 0 .and. v_stat == 0) .or. (s_stat == 0  .and. t_stat == 0) &
      & .or. (v_stat == 0 .and. t_stat == 0)) then
      FLAbort("Multiple field types found for mesh field")
    else if(s_stat == 0) then
      mesh => s_field%mesh
    else if(v_stat == 0) then
      mesh => v_field%mesh
    else if(t_stat == 0) then
      mesh => t_field%mesh
    else
      FLAbort("Mesh field not found")
    end if
  
  end function find_mesh_field  

end module fldiagnostics_module
