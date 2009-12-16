subroutine test_dg_diffusion
  use spud
  use populate_state_module
  use global_parameters
  use state_module
  use unittest_tools
  use vtk_interfaces
  use advection_diffusion_dg
  use sparsity_patterns
  implicit none

  type(state_type), dimension(:), pointer :: states => null()

  type(csr_matrix) :: big_m, mass
  type(csr_sparsity) :: sparsity
  type(scalar_field) :: rhs, tracer
  
  call load_options("test_dg_diffusion_1d.flml")
  call populate_state(states)
  call allocate_and_insert_auxilliary_fields(states)

  tracer = extract_scalar_field(states(1), "Tracer")
  
  call allocate(rhs, tracer%mesh, "RHS")

  sparsity = make_sparsity_transpose(tracer%mesh, tracer%mesh, "Sparsity")

  call allocate(big_m, sparsity, name="Big_m")
  call allocate(mass, sparsity, name="Mass")

  call construct_advection_diffusion_dg(big_m, rhs, "Tracer",&
       & states(1), mass)

  call mmwrite("diffusion.mm", big_m)
  call mmwrite("mass.mm", mass)

end subroutine test_dg_diffusion
