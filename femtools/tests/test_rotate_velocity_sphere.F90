subroutine test_rotate_velocity_sphere

  use fields
  use vtk_interfaces
  use state_module
  use Coordinates
  use unittest_tools
  implicit none

  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: CartesianCoordinate
  type(vector_field), pointer :: UnitRadialVector_inCartesian
  type(vector_field), pointer :: UnitPolarVector_inCartesian
  type(vector_field), pointer :: UnitAzimuthalVector_inCartesian
  type(vector_field), pointer :: UnitRadialVector_inPolar
  type(vector_field), pointer :: UnitPolarVector_inPolar
  type(vector_field), pointer :: UnitAzimuthalVector_inPolar
  type(vector_field) :: diff_RadialVector_inPolar, diff_PolarVector_inPolar, diff_AzimuthalVector_inPolar
  logical :: fail

  call vtk_read_state("data/on_sphere_rotations/spherical_shell_withFields.vtu", state)
  mesh => extract_mesh(state, "Mesh")
  CartesianCoordinate => extract_vector_field(state, "CartesianCoordinate")
  UnitRadialVector_inCartesian => extract_vector_field(state, "UnitRadialVector_inCartesian")
  UnitPolarVector_inCartesian => extract_vector_field(state, "UnitPolarVector_inCartesian")
  UnitAzimuthalVector_inCartesian => extract_vector_field(state, "UnitAzimuthalVector_inCartesian")
  UnitRadialVector_inPolar => extract_vector_field(state, "UnitRadialVector_inPolar")
  UnitPolarVector_inPolar => extract_vector_field(state, "UnitPolarVector_inPolar")
  UnitAzimuthalVector_inPolar => extract_vector_field(state, "UnitAzimuthalVector_inPolar")

  !Test the change of basis from Cartesian to spherical-polar.

  call allocate(diff_RadialVector_inPolar, 3 , mesh, 'difference')
  call allocate(diff_PolarVector_inPolar, 3 , mesh, 'difference')
  call allocate(diff_AzimuthalVector_inPolar, 3 , mesh, 'difference')

  !Set the components difference-vector equal to the unit radial vector, and then apply
  ! transformation to sperical-polar basis. Then compare with vector already in sperical-polar
  ! basis, obtained from vtu.
  call set(diff_RadialVector_inPolar, UnitRadialVector_inCartesian)
  call rotate_velocity_sphere(diff_RadialVector_inPolar, state)
  call addto(diff_RadialVector_inPolar, UnitRadialVector_inPolar, -1.0)
  fail = any(diff_RadialVector_inPolar%val > 1e-12)
  call report_test("[vector basis change: cartesian to spherical-polar of unit-radial vector.]", &
                   fail, .false., "Radial unit vector components not transformed correctly.")
  
  !Set the components difference-vector equal to the unit-polar vector, and then apply
  ! transformation to sperical-polar basis. Then compare with vector already in sperical-polar
  ! basis, obtained from vtu.
  call set(diff_PolarVector_inPolar, UnitPolarVector_inCartesian)
  call rotate_velocity_sphere(diff_PolarVector_inPolar, state)
  call addto(diff_PolarVector_inPolar, UnitPolarVector_inPolar, -1.0)
  fail = any(diff_PolarVector_inPolar%val > 1e-12)
  call report_test("[vector basis change: cartesian to spherical-polar of unit-polar vector.]", &
                   fail, .false., "Polar unit vector components not transformed correctly.")

  !Set the components difference-vector equal to the unit-azimuthal vector, and then apply
  ! transformation to sperical-polar basis. Then compare with vector already in sperical-polar
  ! basis, obtained from vtu.
  call set(diff_AzimuthalVector_inPolar, UnitAzimuthalVector_inCartesian)
  call rotate_velocity_sphere(diff_AzimuthalVector_inPolar, state)
  call addto(diff_AzimuthalVector_inPolar, UnitAzimuthalVector_inPolar, -1.0)
  fail = any(diff_AzimuthalVector_inPolar%val > 1e-12)
  call report_test("[vector basis change: cartesian to spherical-polar of unit-azimuthal vector.]", &
                   fail, .false., "Azimuthal unit vector components not transformed correctly.")

  call deallocate(diff_RadialVector_inPolar)
  call deallocate(diff_PolarVector_inPolar)
  call deallocate(diff_AzimuthalVector_inPolar)

end subroutine
