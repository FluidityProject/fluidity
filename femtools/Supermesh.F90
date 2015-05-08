#include "fdebug.h"

module supermesh_construction
  use fields_data_types
#ifdef HAVE_SUPERMESH
!  use elements, only: allocate, local_coord_count, local_vertices, &
!                     boundary_numbering, operator(==)
!  use libsupermesh_elements, only: element_type
!  use elements
!  use libsupermesh_elements, only: element_type_lib => element_type, &
!                  real_format_lib => real_format

  use libsupermesh_elements, only: element_type_lib => element_type
  use libsupermesh_fields_data_types, only: mesh_faces_lib => mesh_faces, &
                    mesh_subdomain_mesh_lib => mesh_subdomain_mesh, &
                    scalar_field_lib => scalar_field, &
                    vector_field_lib => vector_field, &
                    tensor_field_lib => tensor_field, &
                    mesh_type_lib => mesh_type, &
  scalar_boundary_condition_lib => scalar_boundary_condition, &
  vector_boundary_condition_lib => vector_boundary_condition, &
  scalar_boundary_conditions_ptr_lib => scalar_boundary_conditions_ptr, &
  vector_boundary_conditions_ptr => vector_boundary_conditions_ptr
  use libsupermesh_fields_allocates, only : allocate, deallocate
  use libsupermesh_fields_base, only : ele_count_lib => ele_count, &
                    node_count_lib => node_count, &
                    ele_val_lib => ele_val, &
                    ele_loc_lib => ele_loc
  use libsupermesh_construction
  use libsupermesh_quadrature, only : &
                    quadrature_type_lib => quadrature_type, &
                    make_quadrature_lib => make_quadrature, &
                    deallocate
  use libsupermesh_shape_functions, only : &
                    make_element_shape_lib => make_element_shape, &
                    deallocate
#endif
  use fields_allocates
  use fields_base
#ifdef HAVE_SUPERMESH
!  use fields_allocates, only: deallocate, make_mesh, allocate, zero
!  use libsupermesh_fields_allocates, make_mesh_lib => make_mesh
  
!  use fields_base, only: face_ngi, ele_faces, face_global_nodes, node_val, &
!                         ele_nodes, ele_val, ele_ngi, ele_loc, &
!                         ele_region_id, ele_shape, ele_count, &
!                         ele_val_at_quad
!  use libsupermesh_fields_base
  
!  use libsupermesh_futils, only: real_format_lib => real_format
#endif
  use sparse_tools
  use futils
  use metric_tools
  use unify_meshes_module
  use linked_lists
  use transform_elements
  use global_parameters, only : real_4, real_8
  use tetrahedron_intersection_module
  implicit none

  interface cintersector_set_input
    module procedure intersector_set_input_sp
  
    subroutine cintersector_set_input(nodes_A, nodes_B, ndim, loc)
      use global_parameters, only : real_8
      implicit none
      real(kind = real_8), dimension(ndim, loc), intent(in) :: nodes_A, nodes_B
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
      use global_parameters, only : real_8
      implicit none
      integer, intent(in) :: nonods, totele, ndim, loc
      real(kind = real_8), dimension(nonods * ndim), intent(out) :: nodes
      integer, dimension(totele * loc), intent(out) :: enlist
    end subroutine cintersector_get_output
  end interface cintersector_get_output

#ifdef HAVE_SUPERMESH
  interface intersector_set_dimension
    subroutine libsupermesh_cintersector_set_dimension(ndim)
      implicit none
      integer, intent(in) :: ndim
    end subroutine libsupermesh_cintersector_set_dimension
  end interface intersector_set_dimension

  interface 
    subroutine libsupermesh_cintersector_set_exactness(exact)
      implicit none
      integer, intent(in) :: exact
    end subroutine libsupermesh_cintersector_set_exactness
  end interface
#else
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
  
  ! I hope this is big enough ...
  real, dimension(1024) :: nodes_tmp
  logical :: intersector_exactness = .false.

  public :: intersect_elements, intersector_set_dimension, intersector_set_exactness
  public :: construct_supermesh, compute_projection_error, intersector_exactness
  private

  contains
  
  subroutine intersector_set_input_sp(nodes_A, nodes_B, ndim, loc)
    real(kind = real_4), dimension(ndim, loc), intent(in) :: nodes_A
    real(kind = real_4), dimension(ndim, loc), intent(in) :: nodes_B
    integer, intent(in) :: ndim
    integer, intent(in) :: loc
    
    call cintersector_set_input(real(nodes_A, kind = real_8), real(nodes_B, kind = real_8), ndim, loc)
  
  end subroutine intersector_set_input_sp
  
  subroutine intersector_get_output_sp(nonods, totele, ndim, loc, nodes, enlist)
    integer, intent(in) :: nonods
    integer, intent(in) :: totele
    integer, intent(in) :: ndim
    integer, intent(in) :: loc
    real(kind = real_4), dimension(nonods * ndim), intent(out) :: nodes
    integer, dimension(totele * loc), intent(out) :: enlist
    
    real(kind = real_8), dimension(size(nodes)) :: lnodes
    
    call cintersector_get_output(nonods, totele, ndim, loc, lnodes, enlist)
    nodes = lnodes
    
  end subroutine intersector_get_output_sp

#ifdef HAVE_SUPERMESH
  function intersect_elements(positions_A, ele_A, posB, shape) result(intersection)
    type(vector_field), intent(in) :: positions_A
    integer, intent(in) :: ele_A
    type(vector_field) :: intersection
    type(mesh_type) :: intersection_mesh
    type(element_type), intent(in) :: shape
    real, dimension(:, :), intent(in) :: posB
    
    type(vector_field_lib) :: positions_A_lib
    real, dimension(:,:), allocatable :: positions_A_lib_val
    type(vector_field_lib) :: intersection_lib
    type(element_type_lib) :: shape_lib
!    type(element_type_lib) :: shape_lib_temp
    type(quadrature_type_lib) :: quad_lib
!    type(mesh_type_lib) :: mesh_lib
    type(mesh_type) :: new_mesh
    
!    integer :: dim, loc, dim_lib, loc_lib, dimA, n_count
    integer :: dimA, n_count
    
!    integer :: dim_libbb, vertices_libbb, ele_count_libbb, node_count_libbb, ele_count_libbb_mesh
!    dim_libbb = mesh_dim(positions_A)
!    vertices_libbb = ele_loc(positions_A, 1)
!    ele_count_libbb = ele_count(positions_A)
!    node_count_libbb = node_count(positions_A)
    
!!!    quad_lib = make_quadrature_lib(vertices = positions_A%dim + 2, dim = positions_A%dim + 1, degree = positions_A%mesh%shape%degree * 2) ! JAMES
!!!    shape_lib = make_element_shape_lib(vertices = positions_A%dim + 1, dim = positions_A%dim, degree = positions_A%mesh%shape%degree, quad = quad_lib) ! JAMES

!!!    quad_lib = make_quadrature_lib(vertices = positions_A%mesh%shape%quadrature%vertices, dim = positions_A%mesh%shape%quadrature%dim, ngi = positions_A%mesh%shape%quadrature%ngi, degree = positions_A%mesh%shape%quadrature%degree * 2)
!    quad_lib = make_quadrature_lib(vertices = positions_A%mesh%shape%quadrature%vertices, dim = positions_A%mesh%shape%quadrature%dim, ngi = positions_A%mesh%shape%quadrature%ngi, degree = positions_A%mesh%shape%degree * 2)
!    shape_lib_temp = make_element_shape_lib(vertices = positions_A%mesh%shape%loc, dim = positions_A%mesh%shape%dim, degree = positions_A%mesh%shape%degree, quad = quad_lib)
!    call deallocate(quad_lib)
    
!    call allocate(mesh_lib, node_count(positions_A), ele_count(positions_A), shape_lib_temp)
!    mesh_lib%ndglno = positions_A%mesh%ndglno
!!!    call allocate(positions_A_lib, positions_A%dim, mesh_lib)
 
!!!    allocate(positions_A_lib%val(size(positions_A%val, 1), size(positions_A%val, 2)))
!!!    positions_A_lib%val = positions_A%val 
!!!    positions_A_lib%dim = positions_A%dim
    
!!!    quad_lib = make_quadrature_lib(vertices = shape%quadrature%vertices, dim = shape%quadrature%dim, ngi = shape%quadrature%ngi, degree = shape%quadrature%degree)
    quad_lib = make_quadrature_lib(vertices = shape%quadrature%vertices, dim = shape%quadrature%dim, ngi = shape%quadrature%ngi, degree = shape%degree * 2)
    shape_lib = make_element_shape_lib(vertices = shape%loc, dim = shape%dim, degree = shape%degree, quad = quad_lib)
       
    dimA = positions_A%dim
!    dimB = positions_b%dim
    n_count = 0
    
    select case(positions_A%field_type)
    case(FIELD_TYPE_NORMAL)
      n_count = node_count(positions_A%mesh)
      allocate(positions_a_lib_val(dimA,n_count))
    case(FIELD_TYPE_CONSTANT)
      allocate(positions_a_lib_val(dimA,1))
    case(FIELD_TYPE_DEFERRED)
      allocate(positions_a_lib_val(0,0))
    end select
    
!    positions_A_lib_val = positions_A%val
    positions_A_lib_val = ele_val(positions_A, ele_A)
    
    intersection_lib = libsupermesh_intersect_elements(positions_A_lib_val, ele_count(positions_A), &
        positions_A%mesh%shape%quadrature%vertices, positions_A%mesh%shape%quadrature%dim, &
        ele_A, posB, shape_lib, ele_loc(positions_A, ele_A), dimA, node_count(positions_A), &
        positions_A%mesh%shape%loc, positions_A%field_type, positions_A%mesh%ndglno)
    
    call deallocate(shape_lib)
!    call deallocate(mesh_lib)
!    call deallocate(positions_A_lib)

    call allocate(new_mesh, node_count_lib(intersection_lib), ele_count_lib(intersection_lib), shape)
    new_mesh%ndglno = intersection_lib%mesh%ndglno
    new_mesh%continuity = intersection_lib%mesh%continuity
    
    call allocate(intersection, intersection_lib%dim, new_mesh, name = intersection_lib%name)
    intersection%val = intersection_lib%val
    call deallocate(new_mesh)
    
  end function intersect_elements
#else
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
      call cintersector_get_output(nonods, totele, dim, loc, nodes_tmp, intersection_mesh%ndglno)

      do i = 1, dim
        intersection%val(i,:) = nodes_tmp((i - 1) * nonods + 1:i * nonods)
      end do
    end if

    call deallocate(intersection_mesh)

  end function intersect_elements
#endif

#ifdef HAVE_SUPERMESH
  subroutine intersector_set_exactness(exactness)
    logical, intent(in) :: exactness
    
    call libsupermesh_intersector_set_exactness(exactness)
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
