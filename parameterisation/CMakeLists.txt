target_sources(fluidity PRIVATE
  Equation_of_State.F90
  Gls_vertical_turbulence_model.F90
  iceshelf_meltrate_surf_normal.F90
  k_epsilon.F90
)

add_subdirectory(tests)
