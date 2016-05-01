
subroutine check_options

  use hadapt_extrude, only: hadapt_extrude_check_options
  use fluids_module, only: fluids_module_check_options
  use bulk_parameterisations, only: bulk_parameterisations_check_options
  use sediment, only: sediment_check_options
  use diagnostic_variables, only: diagnostic_variables_check_options
  use adaptive_timestepping, only: adaptive_timestepping_check_options
  use checkpoint, only: checkpoint_check_options
  use field_options, only: field_options_check_options
  use Coordinates, only: Coordinates_check_options
  use write_state_module, only: write_state_module_check_options
  use shape_functions_test, only: shape_functions_test_check_options
  use vertical_extrapolation_module, only: vertical_extrapolation_module_check_options
  use diagnostic_fields_new, only: diagnostic_fields_new_check_options
  use populate_state_module, only: populate_state_module_check_options
  use populate_sub_state_module, only: populate_sub_state_module_check_options
  use k_epsilon, only: k_epsilon_check_options
  use gls, only: gls_check_options
  use momentum_equation, only: momentum_equation_check_options
  use geostrophic_pressure, only: geostrophic_pressure_check_options
  use adaptivity_1d, only: adaptivity_1d_check_options
  use biology, only: biology_check_options
  use compressible_projection, only: compressible_projection_check_options
  use implicit_solids, only: implicit_solids_check_options
  use momentum_DG, only: momentum_DG_check_options
  use adapt_state_prescribed_module, only: adapt_state_prescribed_module_check_options
  use advection_diffusion_fv, only: advection_diffusion_fv_check_options
  use coriolis_module, only: coriolis_module_check_options
  use field_equations_cv, only: field_equations_cv_check_options
  use turbine, only: turbine_check_options
  use adapt_state_module, only: adapt_state_module_check_options
  use advection_diffusion_cg, only: advection_diffusion_cg_check_options
  use qmesh_module, only: qmesh_module_check_options
  use shallow_water_equations, only: shallow_water_equations_check_options
  use mba3d_integration, only: mba3d_integration_check_options
  use sam_integration, only: sam_integration_check_options
  use advection_diffusion_DG, only: advection_diffusion_DG_check_options
  use mba2d_integration, only: mba2d_integration_check_options
  use free_surface_module, only: free_surface_module_check_options
  use interpolation_manager, only: interpolation_manager_check_options
  use adapt_integration, only: adapt_integration_check_options
  use limit_metric_module, only: limit_metric_module_check_options

   continue
     call hadapt_extrude_check_options
  call fluids_module_check_options
  call bulk_parameterisations_check_options
  call sediment_check_options
  call diagnostic_variables_check_options
  call adaptive_timestepping_check_options
  call checkpoint_check_options
  call field_options_check_options
  call Coordinates_check_options
  call write_state_module_check_options
  call shape_functions_test_check_options
  call vertical_extrapolation_module_check_options
  call diagnostic_fields_new_check_options
  call populate_state_module_check_options
  call populate_sub_state_module_check_options
  call k_epsilon_check_options
  call gls_check_options
  call momentum_equation_check_options
  call geostrophic_pressure_check_options
  call adaptivity_1d_check_options
  call biology_check_options
  call compressible_projection_check_options
  call implicit_solids_check_options
  call momentum_DG_check_options
  call adapt_state_prescribed_module_check_options
  call advection_diffusion_fv_check_options
  call coriolis_module_check_options
  call field_equations_cv_check_options
  call turbine_check_options
  call adapt_state_module_check_options
  call advection_diffusion_cg_check_options
  call qmesh_module_check_options
  call shallow_water_equations_check_options
  call mba3d_integration_check_options
  call sam_integration_check_options
  call advection_diffusion_DG_check_options
  call mba2d_integration_check_options
  call free_surface_module_check_options
  call interpolation_manager_check_options
  call adapt_integration_check_options
  call limit_metric_module_check_options

end subroutine check_options
