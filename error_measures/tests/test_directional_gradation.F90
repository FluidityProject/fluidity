subroutine test_directional_gradation

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
  real :: pi, edgelen_val
  real, dimension(3) :: evals
  real, dimension(3, 3) :: evecs, rotmat
  logical :: fail
  real :: ang_sum, abs_ang_sum
  real, dimension(50) :: rand_angle

  call vtk_read_state("data/anisotropic.vtu", state)
  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")

  call allocate(metric, mesh, "Error metric")
  call allocate(edgelen, mesh, "Desired edge lengths")

  pi = 4.0 * atan(1.0)
  ang_sum = 0.0
  abs_ang_sum = 0.0

  rand_angle = (/1.34709966759, 0.311742316029, 10.2424266028, 10.7693084058, -0.936104796565, -0.506890236888, 5.3967205881, &
  -8.29583306098, 6.7317160712, 11.3374319432, 14.8021552737, -6.30494787702, -6.37278249459, -13.2016459733, -12.1030066266, &
  5.22772243646, -9.59888033136, -1.46060726774, -2.31689708218, 10.94115068, 6.83609523304, 10.6183307185, -4.90686553471, &
  -6.24381839029, -11.4155691205, 5.89798704379, 13.2602266611, 11.1242519935, -12.1051829544, 6.99237904103, 12.5013449907, &
  -13.7481918613, 6.35625870005, 12.9701054974, -5.86275722242, -5.35945709316, -5.66120593237, -5.44199708288, -1.96228646539, &
  -7.14006554468, -14.8433497271, -2.24701922775, -11.7559198616, -14.499166861, -10.908333849, 4.11292104312, 7.56519137309, &
  7.15819852375, -10.3970222832, 1.53848011904/)

  do i=1,mesh%nodes
    evals = (/0.10, 0.25, 0.25/)
    evals = eigenvalue_from_edge_length(evals)

    evecs = get_matrix_identity(3)
    rotmat = get_rotation_matrix_cross(-1 * evecs(:, 3), rand_angle(i) * pi / 180.0)
    evecs = matmul(rotmat, evecs)
    call eigenrecomposition(metric%val(:, :, i), evecs, evals)

    ang_sum = ang_sum + rand_angle(i)
    abs_ang_sum = abs_ang_sum + abs(rand_angle(i))
  end do

  !write(0,*) "avg ang_sum: ", (ang_sum / mesh%nodes)
  !write(0,*) "avg abs_ang_sum: ", (abs_ang_sum / mesh%nodes)
  !write(0,*) "mesh%nodes == ", mesh%nodes

  call get_edge_lengths(metric, edgelen)
  edgelen_val = edgelen%val(1) ! constant
  call vtk_write_fields("data/directional_gradation", 0, positions, mesh, &
                        sfields=(/edgelen/), &
                        vfields=(/positions/), &
                        tfields=(/metric/))

  call form_gradation_metric(positions, metric)

  call get_edge_lengths(metric, edgelen)
  call vtk_write_fields("data/directional_gradation", 1, positions, mesh, &
                        sfields=(/edgelen/), &
                        vfields=(/positions/), &
                        tfields=(/metric/))

  fail = .false.
  do i=1,mesh%nodes
    if (edgelen%val(i) .fne. edgelen_val) fail = .true.
    call eigendecomposition_symmetric(metric%val(:, :, i), evecs, evals)
    if (get_angle(dominant_eigenvector(evecs, evals), (/1.0, 0.0, 0.0/)) > 0.10) then
      write(0,*) "i == ", i
      call write_matrix(metric%val(:, :, i), "metric")
      call write_matrix(evecs, "eigenvectors")
      call write_vector(evals, "eigenvalues")
      write(0,*) "dominant_eigenvector(evecs, evals) == ", dominant_eigenvector(evecs, evals)
      write(0,*) "angle == ", get_angle(dominant_eigenvector(evecs, evals), (/1.0, 0.0, 0.0/))
      fail = .true.
    end if
  end do

  call report_test("[directional gradation]", fail, .false., &
                   "The metric field should have constant unchanged edge length, &
                   & and its principal eigenvector should be (1, 0, 0).")

  call deallocate(metric)
  call deallocate(edgelen)

end subroutine test_directional_gradation
