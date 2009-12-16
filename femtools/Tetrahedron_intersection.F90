#define BUF_SIZE 100
#include "fdebug.h"

module tetrahedron_intersection_module

  use elements
  use vector_tools
  use fields_data_types
  use fields_base
  use fields_allocates
  use fields_manipulation
  implicit none

  type tet_type
    real, dimension(3, 4) :: V ! vertices of the tet
  end type tet_type

  type plane_type
    real, dimension(3) :: normal
    real :: c
  end type plane_type

  type(tet_type), dimension(BUF_SIZE) :: tet_array, tet_array_tmp
  integer :: tet_cnt = 0, tet_cnt_tmp = 0
  type(mesh_type), save :: intersection_mesh
  logical, save :: mesh_allocated = .false.

  public :: tet_type, plane_type, intersect_tets, get_planes, finalise_tet_intersector

  interface intersect_tets
    module procedure intersect_tets_dt
  end interface

  interface get_planes
    module procedure get_planes_tet, get_planes_hex
  end interface

  contains

  subroutine finalise_tet_intersector
    if (mesh_allocated) then
      call deallocate(intersection_mesh)
      mesh_allocated = .false.
    end if
  end subroutine finalise_tet_intersector

  subroutine intersect_tets_dt(tetA, planesB, shape, stat, output)
    type(tet_type), intent(in) :: tetA
    type(plane_type), dimension(:), intent(in) :: planesB
    type(element_type), intent(in) :: shape
    type(vector_field), intent(inout) :: output
    integer :: ele
    integer, intent(out) :: stat

    integer :: i, j
    real :: vol
    real, dimension(3) :: vec_tmp

    tet_cnt = 1
    tet_array(1) = tetA

    if (.not. mesh_allocated) then
      call allocate(intersection_mesh, BUF_SIZE * 4, BUF_SIZE, shape, name="IntersectionMesh")
      intersection_mesh%ndglno = (/ (i, i=1,BUF_SIZE*4) /)
      mesh_allocated = .true.
    end if

    do i=1,size(planesB)
      ! Clip the tet_array against the i'th plane
      tet_cnt_tmp = 0

      do j=1,tet_cnt
        call clip(planesB(i), tet_array(j))
      end do

      if (i /= size(planesB)) then
        tet_cnt = tet_cnt_tmp
        tet_array(1:tet_cnt) = tet_array_tmp(1:tet_cnt)
      else
        ! Copy the result if the volume is > epsilon
        tet_cnt = 0
        do j=1,tet_cnt_tmp
          vol = tet_volume(tet_array_tmp(j))
          if (vol < 0.0) then
            vec_tmp = tet_array_tmp(j)%V(:, 1)
            tet_array_tmp(j)%V(:, 1) = tet_array_tmp(j)%V(:, 2)
            tet_array_tmp(j)%V(:, 2) = vec_tmp
            vol = -vol
          end if

          if (vol > epsilon(0.0)) then
            tet_cnt = tet_cnt + 1
            tet_array(tet_cnt) = tet_array_tmp(j)
          end if
        end do
      end if
    end do

    if (tet_cnt == 0) then
      stat=1
      return
    end if

    stat = 0
    intersection_mesh%nodes = tet_cnt*4
    intersection_mesh%elements = tet_cnt
    call allocate(output, 3, intersection_mesh, "IntersectionCoordinates")

    do ele=1,tet_cnt
      call set(output, ele_nodes(output, ele), tet_array(ele)%V)
    end do
  end subroutine intersect_tets_dt

  subroutine clip(plane, tet)
  ! Clip tet against the plane
  ! and append any output to tet_array_tmp.
    type(plane_type), intent(in) :: plane
    type(tet_type), intent(in) :: tet

    real, dimension(4) :: dists
    integer :: neg_cnt, pos_cnt, zer_cnt
    integer, dimension(4) :: neg_idx, pos_idx, zer_idx
    integer :: i

    real :: invdiff, w0, w1
    type(tet_type) :: tet_tmp

    ! Negative == inside
    ! Positive == outside

    neg_cnt = 0
    pos_cnt = 0
    zer_cnt = 0

    dists = distances_to_plane(plane, tet)
    do i=1,4
      if (abs(dists(i)) < epsilon(0.0)) then
        zer_cnt = zer_cnt + 1
        zer_idx(zer_cnt) = i
      else if (dists(i) < 0.0) then
        neg_cnt = neg_cnt + 1
        neg_idx(neg_cnt) = i
      else if (dists(i) > 0.0) then
        pos_cnt = pos_cnt + 1
        pos_idx(pos_cnt) = i
      end if
    end do

    if (neg_cnt == 0) then
      ! tet is completely on positive side of plane, full clip
      return
    end if

    if (pos_cnt == 0) then
      ! tet is completely on negative side of plane, no clip
      tet_cnt_tmp = tet_cnt_tmp + 1
      tet_array_tmp(tet_cnt_tmp) = tet
      return
    end if

    ! The tet is split by the plane, so we have more work to do.

    select case(pos_cnt)
    case(3)
      ! +++-
      tet_cnt_tmp = tet_cnt_tmp + 1
      tet_array_tmp(tet_cnt_tmp) = tet
      do i=1,pos_cnt
        invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
        w0 = -dists(neg_idx(1)) * invdiff
        w1 =  dists(pos_idx(i)) * invdiff
        tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(i)) = &
           w0 * tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(i)) + &
           w1 * tet_array_tmp(tet_cnt_tmp)%V(:, neg_idx(1))
      end do
    case(2)
      select case(neg_cnt)
      case(2)
        ! ++--
        do i=1,pos_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          w0 = -dists(neg_idx(1)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tet_tmp%V(:, i) = w0 * tet%V(:, pos_idx(i)) + w1 * tet%V(:, neg_idx(1))
        end do
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(2)) )
          w0 = -dists(neg_idx(2)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tet_tmp%V(:, i+2) = w0 * tet%V(:, pos_idx(i)) + w1 * tet%V(:, neg_idx(2))
        end do

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp) = tet
        tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(1)) = tet_tmp%V(:, 3)
        tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(2)) = tet_tmp%V(:, 2)

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp)%V(:, 1) = tet%V(:, neg_idx(2))
        tet_array_tmp(tet_cnt_tmp)%V(:, 2) = tet_tmp%V(:, 4)
        tet_array_tmp(tet_cnt_tmp)%V(:, 3) = tet_tmp%V(:, 3)
        tet_array_tmp(tet_cnt_tmp)%V(:, 4) = tet_tmp%V(:, 2)

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp)%V(:, 1) = tet%V(:, neg_idx(1))
        tet_array_tmp(tet_cnt_tmp)%V(:, 2) = tet_tmp%V(:, 1)
        tet_array_tmp(tet_cnt_tmp)%V(:, 3) = tet_tmp%V(:, 2)
        tet_array_tmp(tet_cnt_tmp)%V(:, 4) = tet_tmp%V(:, 3)
      case(1)
        ! ++-0
        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp) = tet
        do i=1,pos_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          w0 = -dists(neg_idx(1)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(i)) = &
             w0 * tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(i)) + &
             w1 * tet_array_tmp(tet_cnt_tmp)%V(:, neg_idx(1))
        end do
      end select
    case(1)
      select case(neg_cnt)
      case(3)
        ! +---
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          w0 = -dists(neg_idx(i)) * invdiff
          w1 =  dists(pos_idx(1)) * invdiff
          tet_tmp%V(:, i) = w0 * tet%V(:, pos_idx(1)) + w1 * tet%V(:, neg_idx(i))
        end do

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp) = tet
        tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(1)) = tet_tmp%V(:, 1)

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp)%V(:, 1) = tet_tmp%V(:, 1)
        tet_array_tmp(tet_cnt_tmp)%V(:, 2) = tet%V(:, neg_idx(2))
        tet_array_tmp(tet_cnt_tmp)%V(:, 3) = tet%V(:, neg_idx(3))
        tet_array_tmp(tet_cnt_tmp)%V(:, 4) = tet_tmp%V(:, 2)

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp)%V(:, 1) = tet%V(:, neg_idx(3))
        tet_array_tmp(tet_cnt_tmp)%V(:, 2) = tet_tmp%V(:, 2)
        tet_array_tmp(tet_cnt_tmp)%V(:, 3) = tet_tmp%V(:, 3)
        tet_array_tmp(tet_cnt_tmp)%V(:, 4) = tet_tmp%V(:, 1)
      case(2)
        ! +--0
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          w0 = -dists(neg_idx(i)) * invdiff
          w1 =  dists(pos_idx(1)) * invdiff
          tet_tmp%V(:, i) = w0 * tet%V(:, pos_idx(1)) + w1 * tet%V(:, neg_idx(i))
        end do

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp) = tet
        tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(1)) = tet_tmp%V(:, 1)

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp)%V(:, 1) = tet_tmp%V(:, 2)
        tet_array_tmp(tet_cnt_tmp)%V(:, 2) = tet%V(:, zer_idx(1))
        tet_array_tmp(tet_cnt_tmp)%V(:, 3) = tet%V(:, neg_idx(2))
        tet_array_tmp(tet_cnt_tmp)%V(:, 4) = tet_tmp%V(:, 1)
      case(1)
        ! +-00
        invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(1)) )
        w0 = -dists(neg_idx(1)) * invdiff
        w1 =  dists(pos_idx(1)) * invdiff

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp) = tet
        tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(1)) = w0 * tet%V(:, pos_idx(1)) + w1 * tet%V(:, neg_idx(1))
      end select
    end select

  end subroutine clip

  pure function get_planes_tet(tet) result(plane)
    type(tet_type), intent(in) :: tet
    type(plane_type), dimension(4) :: plane

    real, dimension(3) :: edge10, edge20, edge30, edge21, edge31
    real :: det
    integer :: i

    edge10 = tet%V(:, 2) - tet%V(:, 1);
    edge20 = tet%V(:, 3) - tet%V(:, 1);
    edge30 = tet%V(:, 4) - tet%V(:, 1);
    edge21 = tet%V(:, 3) - tet%V(:, 2);
    edge31 = tet%V(:, 4) - tet%V(:, 2);
 
    plane(1)%normal = unit_cross(edge20, edge10)
    plane(2)%normal = unit_cross(edge10, edge30)
    plane(3)%normal = unit_cross(edge30, edge20)
    plane(4)%normal = unit_cross(edge21, edge31)

    det = dot_product(edge10, plane(4)%normal)
    if (det < 0) then
      do i=1,4
        plane(i)%normal = -plane(i)%normal
      end do
    end if

    ! And calibrate what is the zero of this plane by dotting with
    ! a point we know to be on it
    do i=1,4
      plane(i)%c = dot_product(tet%V(:, i), plane(i)%normal)
    end do

  end function get_planes_tet

  function get_planes_hex(positions, ele) result(plane)
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: ele
    type(plane_type), dimension(8) :: plane
    integer, dimension(:), pointer :: faces
    integer :: i, face, j
    ! dim x f_loc
    real, dimension(3, 4) :: pos
    real, dimension(3) :: edgeA, edgeB
    integer, dimension(4) :: fnodes
    integer, dimension(8) :: enodes
    integer :: node_not_on_face
    real :: det

    ! This could be done much more efficiently by exploiting
    ! more information about how we number faces and such on a hex

    assert(positions%mesh%shape%numbering%family == FAMILY_CUBE)
    assert(positions%mesh%shape%degree == 1)

    faces => ele_faces(positions, ele)
    assert(size(faces) == 8)
    enodes = ele_nodes(positions, ele)
    do i=1,8
      face = faces(i)
      pos = face_val(positions, face)
      fnodes = face_global_nodes(positions, face)

      ! Only need 3 points to constrain the plane, so pos(:, 4) is never used
      edgeA = pos(:, 2) - pos(:, 1)
      edgeB = pos(:, 3) - pos(:, 1)
      plane(i)%normal = unit_cross(edgeA, edgeB)

      ! Find a node not on the face
      do j=1,8
        if (.not. any(fnodes == enodes(j)) ) then
          node_not_on_face = enodes(j)
          exit
        end if
      end do

      ! Now we need to make sure that we have the correct sense of the normal
      det = dot_product(plane(i)%normal, node_val(positions, node_not_on_face))
      ! The node of the element not on the face should be on the positive side
      ! of the face, so if it is on the negative side we need to invert the normal
      if (det < 0) then
        plane(i)%normal = -plane(i)%normal
      end if

      ! Now we calibrate the constant (setting the 'zero level' of the plane, as it were)
      ! with a node we know is on the face
      plane(i)%c = dot_product(plane(i)%normal, node_val(positions, fnodes(1)))

    end do
  end function get_planes_hex

  pure function unit_cross(vecA, vecB) result(cross)
    real, dimension(3), intent(in) :: vecA, vecB
    real, dimension(3) :: cross
    cross(1) = vecA(2) * vecB(3) - vecA(3) * vecB(2)
    cross(2) = vecA(3) * vecB(1) - vecA(1) * vecB(3)
    cross(3) = vecA(1) * vecB(2) - vecA(2) * vecB(1)

    cross = cross / norm2(cross)
  end function unit_cross

  pure function distances_to_plane(plane, tet) result(dists)
    type(plane_type), intent(in) :: plane
    type(tet_type), intent(in) :: tet
    real, dimension(4) :: dists
    integer :: i

    forall(i=1:4)
      dists(i) = dot_product(plane%normal, tet%V(:, i)) - plane%c
    end forall
  end function distances_to_plane

  pure function tet_volume(tet) result(vol)
    type(tet_type), intent(in) :: tet
    real :: vol
    real, dimension(3) :: cross, vecA, vecB, vecC

    vecA = tet%V(:, 1) - tet%V(:, 4)
    vecB = tet%V(:, 2) - tet%V(:, 4)
    vecC = tet%V(:, 3) - tet%V(:, 4)

    cross(1) = vecB(2) * vecC(3) - vecB(3) * vecC(2)
    cross(2) = vecB(3) * vecC(1) - vecB(1) * vecC(3)
    cross(3) = vecB(1) * vecC(2) - vecB(2) * vecC(1)

    vol = dot_product(vecA, cross) / 6.0
  end function tet_volume

end module tetrahedron_intersection_module
