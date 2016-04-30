add_library_sources(Equation_of_State.F90 Gls_vertical_turbulence_model.F90 k_epsilon.F90 iceshelf_meltrate_surf_normal.F90)
add_subdirectory(tests EXCLUDE_FROM_ALL)

set(unittests ${parameterisation_unittests}${unittests} PARENT_SCOPE)
