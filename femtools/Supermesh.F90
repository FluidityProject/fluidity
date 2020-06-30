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
#ifdef HAVE_LIBSUPERMESH
  use libsupermesh, only : libsupermesh_intersect_elements => intersect_elements
#endif
  use tetrahedron_intersection_module
  implicit none

#ifdef HAVE_LIBSUPERMESH
  real, dimension(:, :, :), allocatable, save :: elements_c
#else
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
  
  ! I hope this is big enough ...
  real, dimension(1024), save :: nodes_tmp
#endif
  logical, save :: intersector_exactness = .false.

  private

  public :: intersect_elements, intersector_set_dimension, intersector_set_exactness
  public :: construct_supermesh, compute_projection_error, intersector_exactness

  contains

#ifdef HAVE_LIBSUPERMESH
  subroutine intersector_set_dimension(ndim)
    integer, intent(in) :: ndim
    
    if(allocated(elements_c)) then
      if(size(elements_c, 1) == ndim) return
      deallocate(elements_c)
    end if
    
    select case(ndim)
      case(1)
        allocate(elements_c(1, 2, 2))
      case(2)
        allocate(elements_c(2, 3, 62))
      case(3)
        allocate(elements_c(3, 4, 3645))
      case default
        FLAbort("Invalid dimension")
    end select
    
  end subroutine intersector_set_dimension

  function intersect_elements(positions_A, ele_A, posB, shape, empty_intersection) result(intersection)
    type(vector_field), intent(in) :: positions_A
    integer, intent(in) :: ele_A
    real, dimension(:, :), intent(in) :: posB
    type(element_type), intent(in) :: shape
    ! if present, returns whether the intersection is empty or not
    ! if present and the intersection is empty, no intersection mesh is returned (should be discarded and not deallocated)
    ! if not present and the intersection is empty, a valid 0-element mesh is returned
    ! this is an optimisation that avoids allocating lots of empty meshes inside supermesh loops
    logical, optional, intent(out) :: empty_intersection
    
    type(vector_field) :: intersection
    
    integer :: i, n_elements_c
    type(mesh_type) :: intersection_mesh

    call libsupermesh_intersect_elements(reordered(ele_val(positions_A, ele_A)), reordered(posB), elements_c, n_elements_c)

    if (present(empty_intersection)) then
      if (n_elements_c==0) then
        empty_intersection = .true.
        return
      else
        empty_intersection = .false.
      end if
    end if

    call allocate(intersection_mesh, size(elements_c, 2) * n_elements_c, n_elements_c, shape, "IntersectionMesh")
    intersection_mesh%continuity = -1
    forall(i = 1:size(elements_c, 2) * n_elements_c)
      intersection_mesh%ndglno(i) = i
    end forall
    
    call allocate(intersection, size(elements_c, 1), intersection_mesh, "IntersectionCoordinates")
    do i = 1, n_elements_c
      call set(intersection, ele_nodes(intersection, i), elements_c(:, :, i))
    end do
    
    call deallocate(intersection_mesh)
    
  contains
  
    function reordered(element)
      ! dim x loc
      real, dimension(:, :), intent(in) :: element
      
      real, dimension(size(element, 1), size(element, 2)) :: reordered
      
      ! See toFluidityElementNodeOrdering in femtools/GMSH_Common.F90
      if(size(element, 1) == 2 .and. size(element, 2) == 4) then
        reordered = element(:, (/1, 2, 4, 3/))
      else if(size(element, 1) == 3 .and. size(element, 2) == 8) then
        reordered = element(:, (/1, 2, 4, 3, 5, 6, 8, 7/))
      else
        reordered = element
      end if
    
    end function reordered

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

  function intersect_elements(positions_A, ele_A, posB, shape, empty_intersection) result(intersection)
    type(vector_field), intent(in) :: positions_A
    integer, intent(in) :: ele_A
    type(vector_field) :: intersection
    type(mesh_type) :: intersection_mesh
    type(element_type), intent(in) :: shape
    real, dimension(:, :), intent(in) :: posB
    ! if present, returns whether the intersection is empty or not
    ! if present and the intersection is empty, no intersection mesh is returned (should be discarded and not deallocated)
    ! if not present and the intersection is empty, a valid 0-element mesh is returned
    ! this is an optimisation that avoids allocating lots of empty meshes inside supermesh loops
    logical, optional, intent(out) :: empty_intersection

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

    if (present(empty_intersection)) then
      if (totele==0) then
         empty_intersection = .true.
         return
      else
         empty_intersection = .false.
      end if
    end if
    
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
    end if

    call deallocate(intersection_mesh)

  end function intersect_elements
#endif

  subroutine intersector_set_exactness(exactness)
    logical, intent(in) :: exactness

#ifdef HAVE_LIBSUPERMESH
    if(exactness) then
      FLAbort("Arbitrary precision arithmetic not supported by libsupermesh")
    end if
#else
    integer :: exact

    if (exactness) then
      exact = 1
    else
      exact = 0
    end if
    intersector_exactness = exactness

    call cintersector_set_exactness(exact)
#endif

  end subroutine intersector_set_exactness

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
    logical :: empty_intersection

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
        intersection(j) = intersect_elements(old_positions, ele_A, pos_B, supermesh_shape, empty_intersection)
        if (empty_intersection) then
          llnode => llnode%next
          cycle
        end if
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
        
        local_coords(1:dim) = node_val(supermesh_positions_remapped, node_C); local_coords(dim+1) = 1.0
        local_coords = matmul(inversion_matrix_A, local_coords)
        old_values = ele_val(old_field, ele_A)
        val = dot_product(eval_shape(ele_shape(old_field, ele_A), local_coords), old_values)
        call set(old_field_on_supermesh, node_C, val)

        local_coords(1:dim) = node_val(supermesh_positions_remapped, node_C); local_coords(dim+1) = 1.0
        local_coords = matmul(inversion_matrix_B, local_coords)
        val = dot_product(eval_shape(supermesh_field_shape, local_coords), element_value)
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
