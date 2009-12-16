subroutine test_pseudo2d_gradation

  use vtk_interfaces
  use field_derivatives
  use unittest_tools
  use state_module
  use global_parameters, only: pseudo2d_coord
  use gradation_metric
  use form_metric_field
  use metric_tools
  implicit none
  
  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  type(scalar_field) :: pressure_field, node_field
  type(tensor_field) :: hessian
  integer :: i, j
  logical :: fail
  real :: x, y, z, angle
  real, dimension(3, 3) :: evecs
  real, dimension(3) :: evals, domvec

  pseudo2d_coord = 3

  call vtk_read_state("data/squat_cube_front.vtu", state)
  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  call allocate(pressure_field, mesh, "Pressure")
  call allocate(node_field, mesh, "Node numbering")
  call allocate(hessian, mesh, "Hessian")

  do i=1,mesh%nodes
    x = positions%val(1)%ptr(i)
    y = positions%val(2)%ptr(i)
    z = positions%val(3)%ptr(i)
    pressure_field%val(i) = x * x
  end do

  call get_node_field(mesh, node_field)

  call compute_hessian(pressure_field, positions, hessian)
  call vtk_write_fields("data/pseudo2d_gradation", 0, positions, mesh, sfields=(/pressure_field, node_field/), tfields=(/hessian/))

  ! At this point the tensors are slightly misaligned due to numerical errors
  ! in computing the hessian. We want out gradation algorithm to fix that.

  do i=1,mesh%nodes
    call eigendecomposition_symmetric(hessian%val(:, :, i), evecs, evals)
    do j=1,3
      if (evals(j) == 0.0) evals(j) = minval(evals, mask=(evals /= 0.0))
    end do
    call eigenrecomposition(hessian%val(:, :, i), evecs, evals)
  end do

  call vtk_write_fields("data/pseudo2d_gradation", 1, positions, mesh, sfields=(/pressure_field, node_field/), tfields=(/hessian/))

  gamma0 = huge(0.0) ! don't change eigenvalues, only eigenvectors
  call form_gradation_metric(positions, hessian)

  call vtk_write_fields("data/pseudo2d_gradation", 2, positions, mesh, sfields=(/pressure_field, node_field/), tfields=(/hessian/))

  fail = .false.
  do i=1,mesh%nodes
    x = positions%val(1)%ptr(i)
    y = positions%val(2)%ptr(i)
    z = positions%val(3)%ptr(i)

    if (x == minval(positions%val(1)%ptr) .or. x == maxval(positions%val(1)%ptr)) cycle
    if (y == minval(positions%val(2)%ptr) .or. y == maxval(positions%val(2)%ptr)) cycle
    
    call eigendecomposition_symmetric(hessian%val(:, :, i), evecs, evals)
    domvec = dominant_eigenvector(evecs, evals)
    angle = get_angle(domvec, (/1.0, 0.0, 0.0/))
    if (angle > 0.17453292519943295) then ! 20 degrees
      write(0,*) "i == ", i
      call write_vector(domvec, "Dominant eigenvector")
      call write_vector(node_val(positions, i), "Position")
      write(0,*) "Angle with x axis:", angle
      fail = .true.
    end if
  end do

  call report_test("[gradation of a hessian]", fail, .false., "Pass!")

  call deallocate(hessian)
  call deallocate(pressure_field)
  call deallocate(state)

end subroutine test_pseudo2d_gradation
