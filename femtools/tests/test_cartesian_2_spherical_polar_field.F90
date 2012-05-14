subroutine test_cartesian_2_spherical_polar_field
  !Test for routine test_cartesian_2_spherical_polar_field, ensuring conversion of a vector field
  ! containing the position vector in Cartesian coordiantes into a vector field containing
  ! the position vector in spherical-polar coordinates is correct.

  use fields
  use vtk_interfaces
  use state_module
  use Coordinates
  use unittest_tools
  implicit none

  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: CartesianCoordinate
  type(vector_field), pointer :: PolarCoordinate
  type(vector_field) :: difference
  logical :: fail

  call vtk_read_state("data/on_sphere_rotations/spherical_shell_withFields.vtu", state)
  mesh => extract_mesh(state, "Mesh")
  CartesianCoordinate => extract_vector_field(state, "CartesianCoordinate")
  PolarCoordinate => extract_vector_field(state, "PolarCoordinate")

  !Apply transformation to Cartesian field obtained from vtu file and 
  ! compare with spherical-polar position-vector.
  call allocate(difference, 3 , mesh, 'difference')
  call cartesian_2_spherical_polar(CartesianCoordinate, difference)
  call addto(difference, PolarCoordinate, -1.0)
  call vtk_write_fields("data/cartesian_2_spherical_polar_field_out", 0, &
                        CartesianCoordinate, mesh, vfields=(/difference/))
  fail = any(difference%val > 1e-8)
  call report_test("[Coordinate change of whole field: Cartesian to spherical-polar.]", &
                   fail, .false., "Position vector components not transformed correctly.")
  
  call deallocate(difference)

end subroutine
