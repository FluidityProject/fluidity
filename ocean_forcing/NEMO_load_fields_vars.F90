! This module contains the temporary fields used in the NEMO data loading process
module NEMO_load_fields_vars

  use fldebug
  use global_parameters, only: OPTION_PATH_LEN, pi, current_debug_level
  use spud
  use fields
  use state_module
 
  type(scalar_field), save :: salinity_t, temperature_t, pressure_t
  type(vector_field), save :: velocity_t
  type(vector_field), save :: position

  contains

subroutine deallocate_temp_fields

  call deallocate(temperature_t)
  call deallocate(salinity_t)
  call deallocate(pressure_t)
  call deallocate(velocity_t)
  call deallocate(position)

end subroutine

end module NEMO_load_fields_vars
