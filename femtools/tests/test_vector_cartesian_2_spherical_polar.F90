!    Copyright (C) 2012 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
subroutine test_vector_cartesian_2_spherical_polar
  !Subroutine/unit-test of correct transformation of vector components from a 
  ! spherical-polar basis to a Cartesian basis. This subroitine obtains the components
  ! of a unit vector in the polar direction, a unit vector in the polar direction and 
  ! a unit vector in the azimuthal direction, in both of the aforementioned bases.
  ! The transformation (tested routine) is applied to the componets in spherical-polar
  ! basis and compared to the cartesian-basis components, one vector at the time.

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
  type(vector_field), pointer :: UnitRadialVector_inCartesian
  type(vector_field), pointer :: UnitPolarVector_inCartesian
  type(vector_field), pointer :: UnitAzimuthalVector_inCartesian
  type(vector_field), pointer :: UnitRadialVector_inPolar
  type(vector_field), pointer :: UnitPolarVector_inPolar
  type(vector_field), pointer :: UnitAzimuthalVector_inPolar
  type(vector_field) :: URVdifference, UPVdifference, UAVdifference
  real, dimension(3) :: XYZ, RTP !Arrays containing a signel node's position vector
                                 ! components in cartesian & spherical-polar coordinates.
  real, dimension(3) :: vector_components_polar     !Array containing components of a
                                                    ! vector in polar basis
  real, dimension(3) :: vector_components_cartesian !Array containing components of a
                                                    ! vector in cartesian basis
  logical :: fail
  integer :: node

  !Extract vector fields from file.
  call vtk_read_state("data/on_sphere_rotations/spherical_shell_withFields.vtu", state)
  mesh => extract_mesh(state, "Mesh")
  CartesianCoordinate => extract_vector_field(state, "CartesianCoordinate")
  PolarCoordinate => extract_vector_field(state, "PolarCoordinate")
  UnitRadialVector_inCartesian => extract_vector_field(state, "UnitRadialVector_inCartesian")
  UnitPolarVector_inCartesian => extract_vector_field(state, "UnitPolarVector_inCartesian")
  UnitAzimuthalVector_inCartesian => extract_vector_field(state, "UnitAzimuthalVector_inCartesian")
  UnitRadialVector_inPolar => extract_vector_field(state, "UnitRadialVector_inPolar")
  UnitPolarVector_inPolar => extract_vector_field(state, "UnitPolarVector_inPolar")
  UnitAzimuthalVector_inPolar => extract_vector_field(state, "UnitAzimuthalVector_inPolar")

  call allocate(URVdifference, 3 , mesh, 'UnitRadialVectorDifference')
  call allocate(UPVdifference, 3 , mesh, 'UnitPolarVectorDifference')
  call allocate(UAVdifference, 3 , mesh, 'UnitAzimuthalVectorDifference')

  !Test correct transforamtion of unit-radial vector.
  do node=1,node_count(URVdifference)
    XYZ = node_val(CartesianCoordinate, node)
    vector_components_cartesian = node_val(UnitRadialVector_inCartesian, node)
    call vector_cartesian_2_spherical_polar(vector_components_cartesian(1), &
                                            vector_components_cartesian(2), &
                                            vector_components_cartesian(3), &
                                            XYZ(1), XYZ(2), XYZ(3), &
                                            vector_components_polar(1), &
                                            vector_components_polar(2), &
                                            vector_components_polar(3), &
                                            RTP(1), RTP(2), RTP(3))
    call set(URVdifference, node, vector_components_polar)
  enddo
  call addto(URVdifference, UnitRadialVector_inPolar, -1.0)
  fail = any(URVdifference%val > 1e-12)
  call report_test( &
          "[vector basis change: Cartesian to spherical-polar of unit-radial vector.]", &
          fail, .false., "Radial unit vector components not transformed correctly.")

  !Test correct transforamtion of unit-polar vector.
  do node=1,node_count(UPVdifference)
    XYZ = node_val(CartesianCoordinate, node)
    vector_components_cartesian = node_val(UnitPolarVector_inCartesian, node)
    call vector_cartesian_2_spherical_polar(vector_components_cartesian(1), &
                                            vector_components_cartesian(2), &
                                            vector_components_cartesian(3), &
                                            XYZ(1), XYZ(2), XYZ(3), &
                                            vector_components_polar(1), &
                                            vector_components_polar(2), &
                                            vector_components_polar(3), &
                                            RTP(1), RTP(2), RTP(3))
    call set(UPVdifference, node, vector_components_polar)
  enddo
  call addto(UPVdifference, UnitPolarVector_inPolar, -1.0)
  fail = any(UPVdifference%val > 1e-12)
  call report_test( &
          "[vector basis change: Cartesian to spherical-polar of unit-polar vector.]", &
          fail, .false., "Polar unit vector components not transformed correctly.")

  !Test correct transforamtion of unit-azimuthal vector.
  do node=1,node_count(UAVdifference)
    XYZ = node_val(CartesianCoordinate, node)
    vector_components_cartesian = node_val(UnitAzimuthalVector_inCartesian, node)
    call vector_cartesian_2_spherical_polar(vector_components_cartesian(1), &
                                            vector_components_cartesian(2), &
                                            vector_components_cartesian(3), &
                                            XYZ(1), XYZ(2), XYZ(3), &
                                            vector_components_polar(1), &
                                            vector_components_polar(2), &
                                            vector_components_polar(3), &
                                            RTP(1), RTP(2), RTP(3))
    call set(UAVdifference, node, vector_components_polar)
  enddo
  call addto(UAVdifference, UnitAzimuthalVector_inPolar, -1.0)
  fail = any(UAVdifference%val > 1e-12)
  call report_test( &
          "[vector basis change: Cartesian to spherical-polar of unit-azimuthal vector.]", &
          fail, .false., "Azimuthal unit vector components not transformed correctly.")

  !Output difference vectors for visualisation.
  call vtk_write_fields("data/test_vector_spherical_polar_2_cartesian_out", 0, &
                        CartesianCoordinate, mesh, &
                        vfields=(/URVdifference, UPVdifference, UAVdifference/))

  call deallocate(URVdifference)
  call deallocate(UPVdifference)
  call deallocate(UAVdifference)

end
