#include "fdebug.h"

module supermesh_construction
  use iso_c_binding, only: c_float, c_double
  use fldebug
  use futils
  use sparse_tools
  use elements
  use fields_data_types
  use fields_base
  use linked_lists
  use fields_allocates
  use fields_manipulation
  use metric_tools
  use unify_meshes_module
  use transform_elements
  use tetrahedron_intersection_module
#ifdef HAVE_SUPERMESH
  use libsupermesh_tri_intersection, only : tri_type, tri_buf_size, &
    & intersect_tris
  use libsupermesh_tet_intersection
  use libsupermesh_interval_intersection, only : intersect_intervals
  use libsupermesh_supermesh, only : intersect_quads, &
    & lib_intersect_elements => intersect_elements
  use libsupermesh_intersection_finder, &
    & cintersector_set_input => rtree_intersection_finder_set_input
#endif
  implicit none

#ifndef HAVE_SUPERMESH
  interface cintersector_set_input
    module procedure intersector_set_input_sp
  
    subroutine cintersector_set_input(nodes_A, nodes_B, ndim, loc)
      use iso_c_binding, only: c_double
      implicit none
      real(kind = c_double), dimension(ndim, loc), intent(in) :: nodes_A, nodes_B
      integer, intent(in) :: ndim, loc
    end subroutine cintersector_set_input
  end interface cintersector_set_input

  interface 
    subroutine cintersector_drive
    end subroutine cintersector_drive
  end interface

  interface
    subroutine cintersector_query(nonods, totele)
      implicit none
      integer, intent(out) :: nonods, totele
    end subroutine cintersector_query
  end interface

  interface cintersector_get_output
    module procedure intersector_get_output_sp
  
    subroutine cintersector_get_output(nonods, totele, ndim, loc, nodes, enlist)
      use iso_c_binding, only: c_double
      implicit none
      integer, intent(in) :: nonods, totele, ndim, loc
      real(kind = c_double), dimension(nonods * ndim), intent(out) :: nodes
      integer, dimension(totele * loc), intent(out) :: enlist
    end subroutine cintersector_get_output
  end interface cintersector_get_output

  interface intersector_set_dimension
    subroutine cintersector_set_dimension(ndim)
      implicit none
      integer, intent(in) :: ndim
    end subroutine cintersector_set_dimension
  end interface intersector_set_dimension

  interface 
    subroutine cintersector_set_exactness(exact)
      implicit none
      integer, intent(in) :: exact
    end subroutine cintersector_set_exactness
  end interface
#endif

#ifndef HAVE_SUPERMESH
  ! I hope this is big enough ...
  real, dimension(1024), save :: nodes_tmp
#endif
  logical, save :: intersector_exactness = .false.
  integer, save, public :: returned_intervals = 0, returned_tris = 0, returned_quads = 0, returned_tets = 0

  private

  public :: intersect_elements, intersector_set_dimension, intersector_set_exactness
  public :: construct_supermesh, compute_projection_error, intersector_exactness

  contains

#ifdef HAVE_SUPERMESH
  subroutine intersector_set_dimension(ndim)
    integer, intent(in) :: ndim
  end subroutine intersector_set_dimension

  function intersect_elements(positions_A, ele_A, posB, shape) result(intersection)
    type(vector_field), intent(in) :: positions_A
    integer, intent(in) :: ele_A
    type(vector_field) :: intersection
    type(element_type), intent(in) :: shape
    real, dimension(:, :), intent(in) :: posB
    
    type(mesh_type) :: intersection_mesh
!    type(tri_type), dimension(tri_buf_size) :: trisC !using this we get hit with a 12% slowdown
!    type(tet_type), dimension(tet_buf_size) :: tetsC !using this we get hit with a 12% slowdown
    real, dimension(1, 2, 1), save            ::  intsC_real
    real, dimension(2, 3, tri_buf_size), save ::  trisC_real
    real, dimension(2, 3, 4 * (2 ** 4) - 2), save ::  quadsC_real
    real, dimension(3, 4, tet_buf_size), save ::  tetsC_real
    type(tet_type) :: tet_A, tet_B
    type(tri_type) :: tri_A, tri_B
    real, dimension(2, 4) :: quad_A, quad_B, quad_temp
    
    integer :: n_trisC, n_intsC, n_quadsC, n_tetsC, i, dim, loc, nonods, totele
#ifdef HAVE_SUPERMESH
    integer, dimension(:, :), allocatable :: ndglno
#endif

    dim = positions_A%dim
    assert(shape%loc == dim + 1)

     if ( dim == 1 ) then
       ! 1D
       call intersect_intervals(ele_val(positions_A, ele_A), posB, intsC_real, n_intsC)
       call allocate(intersection_mesh, n_intsC * 2, n_intsC, shape, "IntersectionMeshInt")
       intersection_mesh%continuity = -1

       if ( n_intsC > 0 ) then
         intersection_mesh%ndglno = (/ (i, i=1,2 * n_intsC) /)
       end if

       call allocate(intersection, positions_A%dim, intersection_mesh, "IntersectionCoordinatesInt")

       if ( n_intsC > 0 ) then
         do i = 1, n_intsC
           call set(intersection, ele_nodes(intersection, i), intsC_real(:,:,i))
         end do
         returned_intervals = returned_intervals + n_intsC
      end if

      call deallocate(intersection_mesh)
     else if ( dim == 2 ) then
       if ( ele_loc(positions_A, ele_A) == 3 ) then
         ! Triangles (2D)
         tri_A%v = ele_val(positions_A, ele_A)
         tri_B%v = posB
         call intersect_tris(tri_A%v, tri_B%v, trisC_real, n_trisC)
         call allocate(intersection_mesh, n_trisC * 3, n_trisC, shape, "IntersectionMeshTri")
         intersection_mesh%continuity = -1

        if ( n_trisC > 0 ) then
          intersection_mesh%ndglno = (/ (i, i=1,3 * n_trisC) /)
        end if

        call allocate(intersection, positions_A%dim, intersection_mesh, "IntersectionCoordinatesTri")

        if ( n_trisC > 0 ) then
          do i = 1, n_trisC
            call set(intersection, ele_nodes(intersection, i), trisC_real(:,:,i)) !trisC(i)%v)
          end do
          returned_tris = returned_tris + n_trisC
        end if

        call deallocate(intersection_mesh)
      else if ( ele_loc(positions_A, ele_A) == 4 ) then
        ! Quads (2D) - Rearrange from Z order to a GMSH consistent order
        quad_temp = ele_val(positions_A, ele_A)
        quad_A(:,1:2) = quad_temp(:,1:2) ; quad_A(:,3) = quad_temp(:,4) ; quad_A(:,4) = quad_temp(:,3)
        quad_temp = posB
        quad_B(:,1:2) = quad_temp(:,1:2) ; quad_B(:,3) = quad_temp(:,4) ; quad_B(:,4) = quad_temp(:,3)
        call lib_intersect_elements(quad_A, quad_B, quadsC_real, n_quadsC)
        call allocate(intersection_mesh, n_quadsC * 3, n_quadsC, shape, "IntersectionMeshQuad")
        intersection_mesh%continuity = -1

        if ( n_quadsC > 0 ) then
          intersection_mesh%ndglno = (/ (i, i=1,3 * n_quadsC) /)
        end if

        call allocate(intersection, positions_A%dim, intersection_mesh, "IntersectionCoordinatesQuad")

        if ( n_quadsC > 0 ) then
          do i = 1, n_quadsC
            call set(intersection, ele_nodes(intersection, i), quadsC_real(:,:,i))
          end do
          returned_quads = returned_quads + n_quadsC
        end if

        call deallocate(intersection_mesh)
       end if
     else if ( dim == 3 ) then
       ! Tets (3D)
       tet_A%v = ele_val(positions_A, ele_A)
       tet_B%v = posB
       call intersect_tets(tet_A%v, tet_B%v, tetsC_real, n_tetsC)
       call allocate(intersection_mesh, n_tetsC * 4, n_tetsC, shape)
       intersection_mesh%continuity = -1

       if ( n_tetsC > 0 ) then
         intersection_mesh%ndglno = (/ (i, i=1,4 * n_tetsC) /)
       end if
 
       call allocate(intersection, positions_A%dim, intersection_mesh, "IntersectionCoordinatesTet")
       if ( n_tetsC > 0 ) then
         do i = 1, n_tetsC
           call set(intersection, ele_nodes(intersection, i), tetsC_real(:,:,i)) !tetsC(i)%v)
         end do
       end if
       returned_tets = returned_tets + n_tetsC
       call deallocate(intersection_mesh)
     else
       ! Unkown D
       FLAbort("Unknown input dimensions.")
     end if

  end function intersect_elements
#else
  subroutine intersector_set_input_sp(nodes_A, nodes_B, ndim, loc)
    real(kind = c_float), dimension(ndim, loc), intent(in) :: nodes_A
    real(kind = c_float), dimension(ndim, loc), intent(in) :: nodes_B
    integer, intent(in) :: ndim
    integer, intent(in) :: loc
    
    call cintersector_set_input(real(nodes_A, kind = c_double), real(nodes_B, kind = c_double), ndim, loc)
  
  end subroutine intersector_set_input_sp
  
  subroutine intersector_get_output_sp(nonods, totele, ndim, loc, nodes, enlist)
    integer, intent(in) :: nonods
    integer, intent(in) :: totele
    integer, intent(in) :: ndim
    integer, intent(in) :: loc
    real(kind = c_float), dimension(nonods * ndim), intent(out) :: nodes
    integer, dimension(totele * loc), intent(out) :: enlist
    
    real(kind = c_double), dimension(size(nodes)) :: lnodes

    call cintersector_get_output(nonods, totele, ndim, loc, lnodes, enlist)
    nodes = lnodes

  end subroutine intersector_get_output_sp

  function intersect_elements(positions_A, ele_A, posB, shape) result(intersection)
    type(vector_field), intent(in) :: positions_A
    integer, intent(in) :: ele_A
    type(vector_field) :: intersection
    type(mesh_type) :: intersection_mesh
    type(element_type), intent(in) :: shape
    real, dimension(:, :), intent(in) :: posB

    integer :: dim, loc
    integer :: nonods, totele
    integer :: i

    dim = positions_A%dim
#ifdef DDEBUG
    select case(dim)
      case(2)
        assert(shape%loc == 3)
      case(3)
        assert(shape%loc == 4)
    end select
#endif

    loc = ele_loc(positions_A, ele_A)

    call cintersector_set_input(ele_val(positions_A, ele_A), posB, dim, loc)
    call cintersector_drive
    call cintersector_query(nonods, totele)
    call allocate(intersection_mesh, nonods, totele, shape, "IntersectionMesh")
    intersection_mesh%continuity = -1
    call allocate(intersection, dim, intersection_mesh, "IntersectionCoordinates")
    if (nonods > 0) then
#ifdef DDEBUG
      intersection_mesh%ndglno = -1
#endif
      call cintersector_get_output(nonods, totele, dim, dim + 1, nodes_tmp, intersection_mesh%ndglno)

      do i = 1, dim
        intersection%val(i,:) = nodes_tmp((i - 1) * nonods + 1:i * nonods)
      end do
      select case(dim)
        case(1)
          returned_intervals = returned_intervals + ele_count(intersection)
        case(2)
          returned_tris = returned_tris + ele_count(intersection)
        case(3)
          returned_tets = returned_tets + ele_count(intersection)
        case default
          FLAbort("Invalid dimension")
      end select
    end if

    call deallocate(intersection_mesh)

  end function intersect_elements
#endif

#ifdef HAVE_SUPERMESH
  subroutine intersector_set_exactness(exactness)
    logical, intent(in) :: exactness

    if(exactness) then
      FLAbort("Arbitrary precision arithmetic not supported by libsupermesh")
    end if

  end subroutine intersector_set_exactness
#else
  subroutine intersector_set_exactness(exactness)
    logical, intent(in) :: exactness
    integer :: exact

    if (exactness) then
      exact = 1
    else
      exact = 0
    end if
    intersector_exactness = exactness

    call cintersector_set_exactness(exact)
  end subroutine intersector_set_exactness
#endif

  ! A higher-level interface to supermesh construction.
  subroutine construct_supermesh(new_positions, ele_B, old_positions, map_BA, supermesh_shape, supermesh)
    type(vector_field), intent(in) :: new_positions, old_positions
    integer, intent(in) :: ele_B
    type(ilist) :: map_BA
    type(element_type), intent(in) :: supermesh_shape
    type(vector_field), intent(out) :: supermesh
    integer :: ele_A
    type(inode), pointer :: llnode
    type(vector_field), dimension(map_BA%length) :: intersection
    real, dimension(new_positions%dim, ele_loc(new_positions, ele_B)) :: pos_B
    type(plane_type), dimension(4) :: planes_B
    type(tet_type) :: tet_A, tet_B
    integer :: lstat, dim, j, i

    dim = new_positions%dim

    if (dim == 3) then
      tet_B%V = ele_val(new_positions, ele_B)
      planes_B = get_planes(tet_B)
    else
      pos_B = ele_val(new_positions, ele_B)
    end if

    j = 1

    llnode => map_BA%firstnode
    do while(associated(llnode))
      ele_A = llnode%value
      if (dim == 3) then
        tet_A%V = ele_val(old_positions, ele_A)
        call intersect_tets(tet_A, planes_B, supermesh_shape, stat=lstat, output=intersection(j))
        if (lstat == 1) then
          llnode => llnode%next
          cycle
        end if
        assert(continuity(intersection(j)) < 0)
      else
        intersection(j) = intersect_elements(old_positions, ele_A, pos_B, supermesh_shape)
        assert(continuity(intersection(j)) < 0)
      end if


      if (ele_count(intersection(j)) > 0) then
        allocate(intersection(j)%mesh%region_ids(ele_count(intersection(j))))
        intersection(j)%mesh%region_ids = ele_A
        j = j + 1
      else
        call deallocate(intersection(j))
      end if

      llnode => llnode%next
    end do

    supermesh = unify_meshes(intersection(1:j-1))
    supermesh%name = "Coordinate"

    do i=1,j-1
      deallocate(intersection(i)%mesh%region_ids)
      call deallocate(intersection(i))
    end do
    call finalise_tet_intersector

  end subroutine construct_supermesh

  subroutine compute_projection_error(old_field, old_positions, supermesh_field_shape, element_value, new_positions, ele_B, supermesh,&
                                    & inversion_matrices_A, inversion_matrix_B, error)
    type(scalar_field), intent(in) :: old_field
    type(vector_field), intent(in) :: old_positions, new_positions, supermesh
    type(element_type), intent(in), target :: supermesh_field_shape
    real, dimension(:), intent(in) :: element_value
    integer, intent(in) :: ele_B
    real, dimension(:, :, :), intent(in) :: inversion_matrices_A
    real, dimension(:, :), intent(in) :: inversion_matrix_B
    real, intent(out) :: error

    type(mesh_type) :: supermesh_field_mesh
    type(scalar_field) :: new_field_on_supermesh, old_field_on_supermesh, projection_error
    real, dimension(ele_loc(old_positions, 1), ele_loc(old_positions, 1)) :: inversion_matrix_A

    integer :: ele_C
    real, dimension(ele_ngi(supermesh, 1)) :: detwei_C
    integer, dimension(:), pointer :: nodes
    integer :: node_C, i, j

    integer :: ele_A
    real, dimension(ele_loc(new_positions, ele_B)) :: local_coords
    integer :: dim
    real :: val
    type(vector_field) :: supermesh_positions_remapped
    real, dimension(ele_loc(old_field, 1)) :: old_values

    dim = old_positions%dim

    supermesh_field_mesh = make_mesh(supermesh%mesh, supermesh_field_shape, -1, "SupermeshFieldMesh")
    call allocate(supermesh_positions_remapped, dim, supermesh_field_mesh, "SupermeshPositionsRemapped")
    call remap_field(supermesh, supermesh_positions_remapped)

    call allocate(old_field_on_supermesh, supermesh_field_mesh, "SubmeshField")
    call zero(old_field_on_supermesh)
    call allocate(new_field_on_supermesh, supermesh_field_mesh, "SubmeshField")
    call zero(new_field_on_supermesh)

    do ele_C=1,ele_count(supermesh)
      nodes => ele_nodes(old_field_on_supermesh, ele_C)
      ele_A = ele_region_id(supermesh, ele_C)
      inversion_matrix_A = inversion_matrices_A(:, :, ele_A)
      do i=1,size(nodes)
        node_C = nodes(i)
        
        val = 0.0
        local_coords(1:dim) = node_val(supermesh_positions_remapped, node_C); local_coords(dim+1) = 1.0
        local_coords = matmul(inversion_matrix_A, local_coords)
        old_values = ele_val(old_field, ele_A)
        do j=1,ele_loc(old_field, ele_A)
          val = val + eval_shape(ele_shape(old_field, ele_A), j, local_coords) * old_values(j)
        end do
        call set(old_field_on_supermesh, node_C, val)

        val = 0.0
        local_coords(1:dim) = node_val(supermesh_positions_remapped, node_C); local_coords(dim+1) = 1.0
        local_coords = matmul(inversion_matrix_B, local_coords)
        do j=1,supermesh_field_shape%loc
          val = val + eval_shape(supermesh_field_shape, j, local_coords) * element_value(j)
        end do
        call set(new_field_on_supermesh, node_C, val)
      end do
    end do

    call deallocate(supermesh_positions_remapped)

    call allocate(projection_error, supermesh_field_mesh, "ProjectionError")
    call set(projection_error, new_field_on_supermesh)
    call addto(projection_error, old_field_on_supermesh, -1.0)

    error = 0.0
    do ele_C=1,ele_count(supermesh)
      call transform_to_physical(supermesh, ele_C, detwei_C)
      assert(ele_ngi(supermesh, ele_C) == ele_ngi(projection_error, ele_C))
      error = error + dot_product(ele_val_at_quad(projection_error, ele_C)**2, detwei_C)
    end do

    call deallocate(supermesh_field_mesh)
    call deallocate(old_field_on_supermesh)
    call deallocate(new_field_on_supermesh)
    call deallocate(projection_error)
  end subroutine compute_projection_error

end module supermesh_construction
