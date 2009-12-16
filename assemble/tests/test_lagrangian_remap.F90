subroutine test_lagrangian_remap
  use lagrangian_remap
  use fields
  use spud
  use populate_state_module
  use vtk_interfaces
  use state_module
  use solvers

  implicit none

  type(state_type), dimension(:), pointer :: states
  type(vector_field), pointer :: velocity, coordinate
  type(scalar_field), dimension(1) :: new_fields, old_fields
  real :: dt
  integer :: i

  call load_options("data/explicit-hyperc-shear-input.flml")
  call populate_state(states)

  call get_option("/timestepping/timestep", dt)

  velocity => extract_vector_field(states(1), "Velocity")
  coordinate => extract_vector_field(states(1), "Coordinate")

  old_fields(1) = extract_scalar_field(states(1), "MaterialVolumeFraction")

  call allocate(new_fields(1), old_fields(1)%mesh, name="NewMaterialVolumeFraction")
  new_fields(1)%option_path = "/field/prognostic/galerkin_projection/continuous"
  call set_solver_options(new_fields(1), ksptype='cg', pctype='eisenstat', rtol=1.0e-7, max_its=10000)
  new_fields(1)%option_path = "/field/prognostic/lagrangian_remap"
  call set_solver_options(new_fields(1), ksptype='cg', pctype='eisenstat', rtol=1.0e-7, max_its=10000)
  new_fields(1)%option_path = "/field"
  call set(new_fields(1), old_fields(1))

  do i = 0, 3000
!     if(mod(i,10)==0) then
      call vtk_write_fields(filename="data/lagrangian_remap", index=i, position=coordinate, model=coordinate%mesh, &
                            sfields=(/old_fields(1), new_fields(1)/), vfields=(/velocity/))
!     end if

    call lagrangian_advection(coordinate, velocity, dt, &
                                    old_fields, new_fields)

    call set(old_fields(1), new_fields(1))
  end do

  call deallocate(new_fields(1))

end subroutine test_lagrangian_remap
