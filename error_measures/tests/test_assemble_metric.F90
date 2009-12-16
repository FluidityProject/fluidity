subroutine test_assemble_metric

  use metric_assemble
  use form_metric_field
  use state_module
  use vtk_interfaces
  use vector_tools
  use unittest_tools
  use mpi

  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: position_field
  type(scalar_field), pointer :: pressure_field
  type(tensor_field) :: metric
  type(metric_options) :: opts
  logical :: fail = .false., warn = .false.
  integer :: i, j
  real :: x, y, z 
  real :: max_eigenbound, min_eigenbound
  real, dimension(3) :: evalues
  real, dimension(3,3) :: evectors
  integer :: ierr

  call vtk_read_state("data/test_spr.vtu", state)

  mesh => extract_mesh(state, "Mesh")
  position_field => extract_vector_field(state, "Coordinate")
  pressure_field => extract_scalar_field(state, "Pressure")
  pressure_field%options%relative = .false.
  pressure_field%options%error = 1.0
  pressure_field%options%min_psi = 1e-5

  do i=1,mesh%nodes
    x = position_field%val(1)%ptr(i)
    y = position_field%val(2)%ptr(i)
    z = position_field%val(3)%ptr(i)
    pressure_field%val(i) = 0.5 * x * x + 0.5 * y * y
  end do

  call allocate(metric, mesh, "Metric")

  opts%min_edge_length = 0.1
  opts%max_edge_length = 1.0
  opts%use_anisotropic_edge_length = .false.
  ! the corresponding eigenbounds are:
  max_eigenbound = 1.0/(opts%min_edge_length * opts%min_edge_length)
  min_eigenbound = 1.0/(opts%max_edge_length * opts%max_edge_length)

  call assemble_metric((/state/), metric, opts)

  do i=1,mesh%nodes
    call eigendecomposition_symmetric(node_val(metric, i), evectors, evalues)
    do j=1,3
      if (evalues(j) .flt. min_eigenbound) then
        print *, "min: i == ", i, "; j == ", j, "; evalue == ", evalues(j), "; min_eigenbound == ", min_eigenbound
        fail = .true.
      end if
      if (evalues(j) .fgt. max_eigenbound) then
        print *, "max: i == ", i, "; j == ", j, "; evalue == ", evalues(j)
        fail = .true.
      end if
    end do
  end do

  call report_test("[eigenbounds]", fail, warn, "Eigenvalues should lie between the bounds set by the options.")
  call vtk_write_fields("data/metric", 0, position_field, mesh, sfields=(/pressure_field/), tfields=(/metric/))
end subroutine test_assemble_metric
