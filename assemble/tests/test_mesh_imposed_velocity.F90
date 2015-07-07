!    Copyright (C) 2006 Imperial College London and others.
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
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
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

#include "fdebug.h" 

subroutine test_mesh_imposed_velocity
  
  use state_module
  use fields
  use vtk_interfaces
  use meshmovement
  use spud
  use unittest_tools
  use futils
  
  ! constants
  character(*), parameter :: FIXTURE_NAME = 'mesh_imposed_velocity'
  logical, parameter :: DEBUG = .true.
  real, parameter :: DT = 1./12.

  ! prescribed variables
  type(state_type), dimension(1) :: states
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  type(vector_field) :: grid_velocity
  type(vector_field) :: old_coordinate

  ! variables under test
  type(vector_field) :: new_coordinate
  
  ! work variables
  type(vector_field) :: expected_coordinate
  integer :: stat

  !==========================================
  ! call tests here
  !==========================================
  call test_linear_velocity()
  call test_rotational_velocity_3d()

  
contains
  
  subroutine setup
    ! mesh_0.vtu is a unit cube
    call vtk_read_state("data/mesh_0.vtu", states(1))

    positions => extract_vector_field(states(1), "Coordinate")

    mesh => positions%mesh
    mesh%name = "CoordinateMesh"

    call insert(states(1), mesh, "CoordinateMesh")
    call insert(states(1), positions, "Coordinate")

    call allocate(old_coordinate, mesh_dim(mesh), mesh, "OldCoordinate")
    call allocate(grid_velocity, mesh_dim(mesh), mesh, "GridVelocity")
    call allocate(new_coordinate, mesh_dim(mesh), mesh, "IteratedCoordinate")
    call allocate(expected_coordinate, mesh_dim(mesh), mesh, "ExpectedCoordinate")

    call set(old_coordinate, positions)
    
    call insert(states(1), old_coordinate, "OldCoordinate")
    call insert(states(1), grid_velocity, "GridVelocity")
    call insert(states(1), new_coordinate, "IteratedCoordinate")
    call insert(states(1), expected_coordinate, "ExpectedCoordinate")
    
    call deallocate(old_coordinate)
    call deallocate(grid_velocity)
    call deallocate(new_coordinate)
    call deallocate(expected_coordinate)

    call add_option("/timestepping/timestep", stat = stat)
    assert(stat == SPUD_NEW_KEY_WARNING)
    call set_option("/timestepping/timestep", DT, stat = stat)
    assert(stat == SPUD_NEW_KEY_WARNING)

    call add_option("/mesh_adaptivity/mesh_movement/imposed_grid_velocity",&
         stat = stat)
    assert(stat == SPUD_NEW_KEY_WARNING)
  end subroutine setup


  subroutine teardown
    call deallocate(states)
    call report_test_no_references()
  end subroutine teardown


  subroutine test_linear_velocity
    call setup()

    ! invent a uniform grid velocity
    do i = 1, node_count(mesh)
       call set(grid_velocity, i, (/ 1., 2., 3./))
    end do
    call assert_nonzero('linear_velocity', 'grid_velocity', grid_velocity)

    ! what are the expected coordinates?  Make sure we are not
    ! inadvertently doing nothing
    call set(expected_coordinate, old_coordinate)
    call addto(expected_coordinate, grid_velocity, scale=DT)
    call assert_not_equal('linear_velocity', &
         'expected_coordinate', expected_coordinate, &
         'old_coordinate', old_coordinate)

    ! compute new_coordinate field from old_coordinate and
    ! grid_velocity, and check it matches
    call move_mesh_imposed_velocity(states)
    call assert_equal('linear_velocity', &
         'new_coordinate', new_coordinate, &
         'expected_coordinate', expected_coordinate)

    call teardown()
  end subroutine test_linear_velocity

  
  subroutine test_rotational_velocity_3d
    real, parameter :: OMEGA = 2*3.1415926535897931
    real :: dtheta
    real, dimension(2, 2) :: rotation
    real, dimension(3) :: r, u
    call setup()

    dtheta = OMEGA*DT
    rotation(1, 1) = cos(dtheta)
    rotation(1, 2) = -sin(dtheta)
    rotation(2, 1) = sin(dtheta)
    rotation(2, 2) = cos(dtheta)

    ! add option to transform to cylindrical coordinates
    call add_option("/mesh_adaptivity/mesh_movement/transform_coordinates",&
         stat = stat)

    ! loop over nodes
    do i = 1, node_count(mesh)
       r = positions%val(:, i)
       
       ! set velocity
       u = (/ -OMEGA*r(2), OMEGA*r(1), 0. /)
       call set(grid_velocity, i, u)

       ! what are the expected coordinates?
       r(1:2) = matmul(rotation, r(1:2))
       call set(expected_coordinate, i, r)
    end do
    
    call assert_nonzero('rotational_velocity', &
         'expected_coordinate', expected_coordinate)

    ! compute new_coordinate field from old_coordinate and
    ! grid_velocity, and check it matches
    call move_mesh_imposed_velocity(states)
    call assert_equal('rotational_velocity_3d', &
         'new_coordinate', new_coordinate, &
         'expected_coordinate', expected_coordinate)

    call teardown()
  end subroutine test_rotational_velocity_3d


  subroutine assert_nonzero(test_name, field_name, field)
    character(*), intent(in) :: test_name
    character(*), intent(in) :: field_name
    type(vector_field), intent(inout) :: field
    real :: maxabsval
    character(20) :: buffer
    logical :: failure

    maxabsval = maxval(abs(field%val))
    write (buffer, '(e10.4)'), maxabsval
    
    failure = (maxabsval .feq. 0.)
    call report_test("["//FIXTURE_NAME//"::"//test_name//"]", &
         failure, .false., trim(field_name)//&
         " value should be nonzero, but max abs value = "//trim(buffer))
  end subroutine assert_nonzero

  
  subroutine assert_equal(test_name, field1_name, field1, field2_name, field2) 
    character(*), intent(in) :: test_name
    character(*), intent(in) :: field1_name, field2_name
    type(vector_field), intent(inout) :: field1, field2
    call assert_with_vector_fields(&
       test_name, field1_name, field1, field2_name, field2, .true.) 
  end subroutine assert_equal


  subroutine assert_not_equal(test_name, field1_name, field1, field2_name, field2) 
    character(*), intent(in) :: test_name
    character(*), intent(in) :: field1_name, field2_name
    type(vector_field), intent(inout) :: field1, field2
    call assert_with_vector_fields(&
       test_name, field1_name, field1, field2_name, field2, .false.) 
  end subroutine assert_not_equal


  subroutine assert_with_vector_fields(&
       test_name, field1_name, field1, field2_name, field2, equals) 
    character(*), intent(in) :: test_name
    character(*), intent(in) :: field1_name, field2_name
    type(vector_field), intent(inout) :: field1, field2
    logical, intent(in) :: equals
    type(vector_field) :: vector_work
    real :: maxabsval
    character(40) :: stipulation
    character(40) :: buffer
    logical :: failure
    call allocate(vector_work, mesh_dim(mesh), mesh, "VectorWork")
    
    call set(vector_work, field1)
    call addto(vector_work, field2, scale=-1.)
    
    maxabsval = maxval(abs(vector_work%val))
    write (buffer, '(e10.4)'), maxabsval

    if (equals) then
       failure = .not. (maxabsval .feq. 0.)
       stipulation = "should be equal to"
    else
       failure = (maxabsval .feq. 0.)
       stipulation = "should not be equal to"
    end if
    
    call report_test("["//FIXTURE_NAME//"::"//test_name//"]", &
         failure, .false., &
         trim(field1_name)//" "//trim(stipulation)//" "//&
         trim(field2_name)//", but max abs value of difference = "//&
         trim(buffer))

    if (failure) then
       call print_vectors(3, field1_name, field1, field2_name, field2)
    end if

    call deallocate(vector_work)
  end subroutine assert_with_vector_fields

  
  subroutine print_vectors( dim, field1_name, field1, field2_name, field2)
    integer, intent(in) :: dim
    character(*), intent(in) :: field1_name
    character(*), intent(in), optional :: field2_name
    type(vector_field), intent(inout) :: field1
    type(vector_field), intent(inout), optional :: field2
    integer, parameter :: NSTART = 7
    integer, parameter :: NEND = 3
    integer, parameter :: COLWIDTH = 12
    character(78) :: buffer
    character(39) :: semi
    character(COLWIDTH) :: segment
    integer :: nnodes, nfields, i, j, k
    real :: val
    
    nnodes = size(field1%val, 2) 
    if (present(field2_name) .and. present(field2)) then
       nfields = 2
    else
       nfields = 1
    end if

    write (0, '(a)') ''
    nodes: do k = 0, nnodes
       buffer(:) = ' '
       
       ! skip the possibly large bunch of middle nodes
       if ((k .gt. NSTART+1) .and. (nnodes-k+1 .gt. NEND)) cycle

       fields: do j = 1, nfields
          semi(:) = ' '

          ! field separator
          if (j .eq. 2) then
             buffer = buffer(1:3*COLWIDTH)//'  |'
          end if
          
          ! header
          if (k .eq. 0) then
             if (j .eq. 1) then
                semi = field1_name
             else
                semi = field2_name
             end if

          else
             dims: do i = 1, 3
                segment(:) = ' '

                if (i .le. dim) then
                   if ((k .eq. NSTART+1) .and. (nnodes-k .gt. NEND)) then
                      ! before skipping the middle, write symbol to say 'etc.'
                      segment(6:6) = ':'
                   else
                      ! within domain dimensions, write value
                      if (j .eq. 1) then
                         val = field1%val(i, k)
                      else
                         val = field2%val(i, k)
                      end if
                      write (segment, '(e'//int2str(COLWIDTH)//'.4)'), val
                   end if
                end if

                semi = semi(1:(i-1)*COLWIDTH)//segment
             end do dims
          end if
          
          buffer = buffer(1:(j-1)*40)//semi
       end do fields
       
       write (0, '(a)') '  '//buffer
    end do nodes
    
    write (0, '(a)') ''
  end subroutine print_vectors

  
end subroutine test_mesh_imposed_velocity
