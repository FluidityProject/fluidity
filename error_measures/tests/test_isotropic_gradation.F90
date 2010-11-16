subroutine test_isotropic_gradation

  use vtk_interfaces
  use metric_assemble
  use unittest_tools
  use adapt_state_module
  use state_module
  use mpi
  use edge_length_module
  use gradation_metric
  use global_parameters, only: pseudo2d_coord
  implicit none

  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  type(tensor_field) :: metric
  type(scalar_field) :: edgelen
  integer :: i, ierr
  logical :: fail
  real :: x, y, z
  real, dimension(3, 3) :: x_metric, other_metric

  pseudo2d_coord = 3

  gamma0 = 1.0

  x_metric(1, :) = (/100.0, 0.0, 0.0/)
  x_metric(2, :) = (/0.0, 100.0, 0.0/)
  x_metric(3, :) = (/0.0, 0.0, 100.0/)

  other_metric = x_metric / 25.0 ! 5 times as big

  call vtk_read_state("data/pseudo2d.vtu", state)
  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  
  call allocate(metric, mesh, "Metric")
  call allocate(edgelen, mesh, "Edge lengths")

  do i=1,mesh%nodes
    x = positions%val(1,i)
    y = positions%val(2,i)
    z = positions%val(3,i)
    if (x == 0.0) then
      metric%val(:, :, i) = x_metric
    else
      metric%val(:, :, i) = other_metric
    end if
  end do


  call get_edge_lengths(metric, edgelen)
  call vtk_write_fields("data/gradation_isotropic", 0, positions, mesh, sfields=(/edgelen/), tfields=(/metric/))
  call form_gradation_metric(positions, metric)
  call get_edge_lengths(metric, edgelen)
  call vtk_write_fields("data/gradation_isotropic", 1, positions, mesh, sfields=(/edgelen/), tfields=(/metric/))

  fail = .false.
  do i=1,mesh%nodes
    if (edgelen%val(i) .fne. 0.1) fail = .true.
  end do
  call report_test("[isotropic gradation 1]", fail, .false., "Gradation applied with a smoothening factor of 1.0 &
                  & should give a constant edge length field.")

  do i=1,mesh%nodes
    x = positions%val(1,i)
    y = positions%val(2,i)
    z = positions%val(3,i)
    if (x == 0.0) then
      metric%val(:, :, i) = x_metric
    else
      metric%val(:, :, i) = other_metric
    end if
  end do

  gamma0 = 1.1
  call form_gradation_metric(positions, metric)
  call get_edge_lengths(metric, edgelen)
  call vtk_write_fields("data/gradation_isotropic", 2, positions, mesh, sfields=(/edgelen/), tfields=(/metric/))

  fail = .false.
  do i=1,mesh%nodes
    x = positions%val(1,i)
    y = positions%val(2,i)
    z = positions%val(3,i)
    if (x == 0.0 .and. (edgelen%val(i) .fne. 0.1)) fail = .true.
    if (x >  0.0 .and. x < 1.25 .and. edgelen%val(i) > 0.3) fail = .true.
    if (x >= 5.0 .and. (edgelen%val(i) .fne. 0.5)) fail = .true.
  end do
  call report_test("[isotropic gradation 2]", fail, .false., "Gradation applied with a smoothening factor of 1.1 &
                   & should give a smoothly varying edge length field.")

end subroutine test_isotropic_gradation
