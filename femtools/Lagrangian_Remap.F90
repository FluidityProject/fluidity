

module lagrangian_remap
  use conservative_interpolation_module
  use fields
  use vtk_interfaces
  use interpolation_module
  use solvers
  use sparsity_patterns
  use sparse_tools
  use sparse_matrices_fields

  implicit none

contains

  subroutine lagrangian_advection(old_position, velocity, dt, &
                                  old_fields, new_fields)

    type(vector_field), intent(in) :: old_position
    type(vector_field), intent(in) :: velocity
    real, intent(in) :: dt
    type(scalar_field), dimension(:), intent(in), target :: old_fields
    type(scalar_field), dimension(:), intent(inout) :: new_fields

    type(state_type) :: old_state, new_state

    type(vector_field) :: new_position, unadvected_velocity
    integer :: i
    integer :: subcycle_factor
    type(element_type), pointer :: shape
    type(mesh_type), pointer :: mesh
    type(csr_matrix) :: mass_matrix_new, mass_matrix_old
    type(csr_sparsity) :: sparsity
    type(scalar_field), dimension(size(old_fields)) :: rhs, advected_fields
    integer :: field_count, field, ele
    real, dimension(ele_ngi(old_fields(1), 1)) :: detwei
    real, dimension(ele_loc(old_fields(1), 1), ele_loc(old_fields(1), 1)) :: little_mass_matrix
    
    type(state_type), dimension(1) :: old_interpolation_state, new_interpolation_state

    call insert(old_state, old_position, "Coordinate")
    call insert(old_state, velocity, "Velocity")
    call insert(old_state, old_position%mesh, "Mesh")

    call allocate(new_position, old_position%dim, old_position%mesh, name="Coordinate")
    call zero(new_position)

    call allocate(unadvected_velocity, velocity%dim, velocity%mesh, name="Velocity")
    call zero(unadvected_velocity)
    call insert(new_state, new_position, "Coordinate")
    call insert(new_state, unadvected_velocity, "Velocity")
    call insert(new_state, new_position%mesh, "Mesh")

    call addto(new_position, old_position)
    subcycle_factor = 100
    do i=1,subcycle_factor
      call linear_interpolation(old_state, new_state)
      call addto(new_position, unadvected_velocity, dt/subcycle_factor)
    end do
!     call addto(new_position, velocity, dt)

    ! Now we solve the mass matrix equation to make it conservative, na ja?
    mesh => old_fields(1)%mesh
    field_count = size(old_fields)
    shape => ele_shape(mesh, 1)

    sparsity = make_sparsity(mesh, mesh, name="MassMatrixSparsity")
    call allocate(mass_matrix_old, sparsity, name="OldMassMatrix")
    call zero(mass_matrix_old)
    call allocate(mass_matrix_new, sparsity, name="NewMassMatrix")
    call zero(mass_matrix_new)
    call deallocate(sparsity)

    do ele=1,ele_count(mesh)
      call transform_to_physical(old_position, ele, detwei=detwei)
      little_mass_matrix = shape_shape(shape, shape, detwei)
      call addto(mass_matrix_old, ele_nodes(old_position, ele), ele_nodes(old_position, ele), little_mass_matrix)

      call transform_to_physical(new_position, ele, detwei=detwei)
      little_mass_matrix = shape_shape(shape, shape, detwei)
      call addto(mass_matrix_new, ele_nodes(new_position, ele), ele_nodes(new_position, ele), little_mass_matrix)
    end do

    do field=1,field_count
      call allocate(rhs(field), mesh, "Rhs" // int2str(field))
      call mult(rhs(field), mass_matrix_old, old_fields(field))

      call allocate(advected_fields(field), mesh, "AdvectedField" // int2str(field))
      call set(advected_fields(field), old_fields(field))
    end do

    call deallocate(mass_matrix_old)

    call petsc_solve(advected_fields, mass_matrix_new, rhs, option_path=trim(new_fields(1)%option_path) &
                                              & // "/prognostic/lagrangian_remap")
    call deallocate(mass_matrix_new)

    do field=1,field_count
      call deallocate(rhs(field))
    end do

    do field = 1, field_count
      call insert(old_interpolation_state(1), advected_fields(field), name=trim(advected_fields(field)%name))
      call insert(new_interpolation_state(1), new_fields(field), name=trim(new_fields(field)%name))
    end do

    call interpolation_galerkin(old_interpolation_state, new_position, &
                                new_interpolation_state, old_position, &
                                force_bounded=.true.)
    
    call deallocate(old_interpolation_state(1))
    call deallocate(new_interpolation_state(1))

    do field=1,field_count
      call deallocate(advected_fields(field))
    end do

    call deallocate(new_position)
    call deallocate(unadvected_velocity)
    call deallocate(new_state)
    call deallocate(old_state)

  end subroutine lagrangian_advection

end module lagrangian_remap
