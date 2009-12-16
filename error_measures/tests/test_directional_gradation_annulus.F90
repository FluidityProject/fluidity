subroutine test_directional_gradation_annulus

  use gradation_metric
  use unittest_tools
  use state_module
  use vtk_interfaces
  use vector_tools
  use metric_tools
  use edge_length_module
  implicit none

  integer :: i
  type(tensor_field) :: metric
  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  type(scalar_field) :: edgelen
  real :: x, y, z, pi, rand, rand2, sign
  real, dimension(3) :: evals
  real, dimension(3, 3) :: evecs, rotmat
  logical :: fail

  call vtk_read_state("data/squashed_annulus.vtu", state)
  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")

  call allocate(metric, mesh, "Error metric")
  call allocate(edgelen, mesh, "Desired edge lengths")

  pi = 4.0 * atan(1.0)

  do i=1,mesh%nodes
    x = positions%val(1)%ptr(i)
    y = positions%val(2)%ptr(i)
    z = positions%val(3)%ptr(i)

    call random_number(rand)
    call random_number(rand2)
    if (rand2 > 0.5) then
      sign = 1.0
    else
      sign = -1.0
    end if
    evals = (/0.10, 0.25, 0.25/)
    evals = eigenvalue_from_edge_length(evals)

    ! Make the metric point radially
    evecs(:, 1) = (/x, y, 0.0/)
    evecs(:, 2) = (/-y, x, 0.0/)
    evecs(:, 3) = (/0.0, 0.0, 1.0/)

    rotmat = get_rotation_matrix_cross(-1 * evecs(:, 3), sign * rand * pi / 24.0)
    evecs = matmul(rotmat, evecs)

    call eigenrecomposition(metric%val(:, :, i), evecs, evals)
  end do

  call get_edge_lengths(metric, edgelen)
  call vtk_write_fields("data/directional_gradation_annulus", 0, positions, mesh, &
                        sfields=(/edgelen/), &
                        tfields=(/metric/))

  call form_gradation_metric(positions, metric)

  call get_edge_lengths(metric, edgelen)
  call vtk_write_fields("data/directional_gradation_annulus", 1, positions, mesh, &
                        sfields=(/edgelen/), &
                        tfields=(/metric/))

  fail = .false.
  do i=1,mesh%nodes
    x = positions%val(1)%ptr(i)
    y = positions%val(2)%ptr(i)
    z = positions%val(3)%ptr(i)

    call eigendecomposition_symmetric(metric%val(:, :, i), evecs, evals)
    if (get_angle(dominant_eigenvector(evecs, evals), (/x, y, 0.0/)) > 0.4) then
      write(0,*) "i == ", i
      call write_matrix(metric%val(:, :, i), "metric")
      call write_matrix(evecs, "eigenvectors")
      call write_vector(evals, "eigenvalues")
      write(0,*) "dominant_eigenvector(evecs, evals) == ", dominant_eigenvector(evecs, evals)
      write(0,*) "angle == ", get_angle(dominant_eigenvector(evecs, evals), (/x, y, 0.0/))
      fail = .true.
    end if
  end do

  call report_test("[directional gradation on the annulus]", fail, .false., &
                   "The metric field should have constant unchanged edge length, &
                   & and its principal eigenvector should be radial.")

  call deallocate(metric)
  call deallocate(edgelen)

end subroutine test_directional_gradation_annulus
