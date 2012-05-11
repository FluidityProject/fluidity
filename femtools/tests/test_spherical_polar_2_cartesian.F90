subroutine test_spherical_polar_2_cartesian
  !Test for ensuring the position vector at a point is correctly converted from a spherical-polar
  ! basis to a cartesian basis.
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
  integer :: node
  real, dimension(3) :: X, RTP
  logical :: fail

  call vtk_read_state("data/on_sphere_rotations/spherical_shell_withFields.vtu", state)
  mesh => extract_mesh(state, "Mesh")
  CartesianCoordinate => extract_vector_field(state, "CartesianCoordinate")
  PolarCoordinate => extract_vector_field(state, "PolarCoordinate")

  !Extract the components of points in vtu file in polar coordinates, apply transformation and 
  ! compare with cartesian position-vector.
  call allocate(difference, 3 , mesh, 'difference')
  do node=1,node_count(PolarCoordinate)
    RTP = node_val(PolarCoordinate, node)
    call spherical_polar_2_cartesian(RTP(1), RTP(2), RTP(3), X(1), X(2), X(3))
    call set(difference, node, X)
  enddo
  call addto(difference, CartesianCoordinate, -1.0)
  call vtk_write_fields("data/test_spherical_polar_2_cartesian_difference", 0, &
                        CartesianCoordinate, mesh, vfields=(/difference/))
  fail = any(difference%val > 1e-8)
  call report_test("[Coordinate change: Spherical-polar to Cartesian.]", &
                   fail, .false., "Position vector components not transformed correctly.")
  
  call deallocate(difference)

end subroutine
