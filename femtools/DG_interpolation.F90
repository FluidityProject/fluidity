#include "fdebug.h"

module dg_interpolation_module

  use vector_tools
  use quadrature
  use futils
  use spud
  use sparse_tools
  use tensors
  use sparse_tools
  use transform_elements
  use adjacency_lists
  use unittest_tools
  use linked_lists
  use supermesh_construction
  use fetools
  use intersection_finder_module
  use fields
  use meshdiagnostics
  use sparsity_patterns
  use vtk_interfaces
  use state_module
  use interpolation_module
  use sparse_matrices_fields
  use solvers
  implicit none

  private
  public :: dg_interpolation_galerkin_supermesh_free

  interface dg_interpolation_galerkin_supermesh_free
    module procedure dg_interpolation_galerkin_scalars_supermesh_free, dg_interpolation_galerkin_state_supermesh_free
  end interface

  contains

  subroutine dg_interpolation_galerkin_scalars_supermesh_free(old_fields, old_position, new_fields, new_position, map_BA)
    type(scalar_field), dimension(:), intent(in) :: old_fields
    type(vector_field), intent(in) :: old_position

    type(scalar_field), dimension(:), intent(inout), target :: new_fields
    type(vector_field), intent(in) :: new_position
    type(ilist), dimension(:), intent(in), optional, target :: map_BA

    integer :: ele_B
    integer :: ele_A
    integer :: j, k, l, field, field_count
    real, dimension(:), allocatable :: detwei_B

    ! We want to compute the mixed mass matrix M^{BA}.
    ! But that's huge. So, we compute a part of a row (not even the whole row)
    ! and multiply it by a part of the solution on the old mesh A
    ! to get a component of the RHS of the matrix we want to solve.
    real, dimension(ele_loc(new_fields(1), 1), ele_loc(new_fields(1), 1)) :: mat_int
    real, dimension(ele_loc(new_fields(1), 1)) :: val_A
    real, dimension(ele_loc(new_fields(1), 1), size(new_fields)) :: rhs
    ! For each element in B, we will need to identify the local coordinates in B
    ! of the positions of the gauss points of all its children C elements.
    ! So we'll need to assemble and invert that matrix (the global-to-local inversion matrix):
    real, dimension(ele_loc(new_position, 1), ele_loc(new_position, 1)) :: inversion_matrix_B, inversion_matrix_A
    real, dimension(ele_loc(new_position, 1), ele_loc(new_position, 1), ele_count(old_position)) :: inversion_matrices_A
    real, dimension(ele_loc(new_fields(1), 1), ele_loc(new_fields(1), 1)) :: mat, little_mass_matrix

    real, dimension(:, :), allocatable :: pos_at_quad_A
    real, dimension(:, :), allocatable :: basis_at_quad_B, basis_at_quad_A
    integer :: nloc
    logical :: P1
    integer :: dim
    type(ilist), dimension(:), pointer :: lmap_BA
    type(inode), pointer :: llnode
    type(quadrature_type) :: new_quad, old_quad
    type(mesh_type), pointer :: B_mesh
    integer :: ngi
    type(element_type) :: new_shape, old_shape, new_position_inc_quad_shape
    type(vector_field) :: new_position_inc_quad

    ! Linear positions -- definitely linear positions.
    assert(old_position%mesh%shape%degree == 1)
    assert(continuity(old_position) >= 0)
    assert(continuity(new_position) >= 0)

    field_count = size(old_fields)
    nloc = ele_loc(new_fields(1), 1)
    dim = new_position%dim

    B_mesh => new_fields(1)%mesh
    old_shape = new_fields(1)%mesh%shape

    new_quad = make_quadrature(vertices=ele_loc(new_position, 1), dim=dim, degree=8)
    new_shape = make_element_shape(vertices=ele_loc(new_position, 1), dim=dim, degree=old_shape%degree, quad=new_quad)
    ngi = new_shape%ngi

    ! I can't believe how ugly this is - someone please write a "swap-out ngi"
    ! routine
    new_position_inc_quad_shape = make_element_shape(vertices=ele_loc(new_position, 1), dim=dim, degree=new_position%mesh%shape%degree, quad=new_quad)
    new_position_inc_quad = new_position
    new_position_inc_quad%mesh%shape = new_position_inc_quad_shape

    assert(continuity(new_fields(1)) < 0)
    P1 = (new_shape%degree == 1)
    call intersector_set_dimension(dim)
    if (present(map_BA)) then
      lmap_BA => map_BA
    else
      allocate(lmap_BA(ele_count(new_position_inc_quad)))
      lmap_BA = intersection_finder(new_position_inc_quad, old_position)
    end if

    old_quad = B_mesh%shape%quadrature
    B_mesh%shape = new_shape

    do field=1,field_count
      call zero(new_fields(field))
    end do

    allocate(pos_at_quad_A(new_position_inc_quad%dim+1, ngi))
    allocate(basis_at_quad_A(nloc, ngi))
    allocate(basis_at_quad_B(nloc, ngi))
    allocate(detwei_B(ngi))

    do ele_A=1,ele_count(old_position)
      call local_coords_matrix(old_position, ele_A, inversion_matrices_A(:, :, ele_A))
    end do

    do ele_B=1,ele_count(new_position_inc_quad)

      ! First thing: assemble and invert the inversion matrix.
      call local_coords_matrix(new_position_inc_quad, ele_B, inversion_matrix_B)

      ! Second thing: assemble the mass matrix of B on the left.
      call transform_to_physical(new_position_inc_quad, ele_B, detwei=detwei_B)
      little_mass_matrix = shape_shape(new_shape, new_shape, detwei_B)

      llnode => lmap_BA(ele_B)%firstnode
      rhs = 0

      do while (associated(llnode))
        ele_A = llnode%value
        inversion_matrix_A = inversion_matrices_A(:, :, ele_A)

        basis_at_quad_B = new_shape%n

        pos_at_quad_A(1:dim, :) = ele_val_at_quad(new_position_inc_quad, ele_B)
        pos_at_quad_A(dim+1, :) = 1.0
        pos_at_quad_A = matmul(inversion_matrix_A, pos_at_quad_A)

        ! Evaluate the basis functions at the local coordinates
        if (P1) then
          do j=1,ngi
            ! Check if it's inside the element or not
            if (min(minval(pos_at_quad_A(:, j)), minval(1 - pos_at_quad_A(:, j))) >= 0.0) then
              basis_at_quad_A(:, j) = pos_at_quad_A(:, j)
            else
              basis_at_quad_A(:, j) = 0.0
            end if
          end do
        else
          do j=1,ngi
            if (min(minval(pos_at_quad_A(:, j)), minval(1 - pos_at_quad_A(:, j))) >= 0.0) then
              basis_at_quad_A(:, j) = eval_shape(new_shape, pos_at_quad_A(:, j))
            else
              basis_at_quad_A(:, j) = 0.0
            end if
          end do
        end if
          
        ! Combined outer_product and tensormul_3_1 to see if it is faster.
        mat = 0.0
        do j=1,ngi
          forall (k=1:nloc,l=1:nloc)
            mat_int(k, l) = basis_at_quad_B(k, j) * basis_at_quad_A(l, j)
          end forall
          mat = mat + mat_int * detwei_B(j)
        end do

        do field=1,field_count
          val_A = ele_val(old_fields(field), ele_A)
          rhs(:, field) = rhs(:, field) + matmul(mat, val_A)
        end do

        llnode => llnode%next
      end do

      call solve(little_mass_matrix, rhs)
      do field=1,field_count
        call set(new_fields(field), ele_nodes(new_fields(field), ele_B), rhs(:, field))
      end do
    end do

    if (.not. present(map_BA)) then
      do ele_B=1,ele_count(new_position_inc_quad)
        call deallocate(lmap_BA(ele_B))
      end do
      deallocate(lmap_BA)
    end if

    deallocate(pos_at_quad_A)
    deallocate(detwei_B)

    B_mesh%shape = old_shape
    call deallocate(new_quad)
    call deallocate(new_shape)

  end subroutine dg_interpolation_galerkin_scalars_supermesh_free

  subroutine dg_interpolation_galerkin_state_supermesh_free(old_state, new_state, map_BA)
    type(state_type), intent(inout) :: old_state, new_state

    type(scalar_field), dimension(:), pointer :: old_fields, new_fields
    type(vector_field), pointer :: old_position, new_position
    type(ilist), dimension(:), intent(in), optional :: map_BA

    call collapse_state(old_state, old_fields)
    call collapse_state(new_state, new_fields)

    old_position => extract_vector_field(old_state, "Coordinate")
    new_position => extract_vector_field(new_state, "Coordinate")

    if(present(map_BA)) then
      call dg_interpolation_galerkin_scalars_supermesh_free(old_fields, old_position, new_fields, new_position, map_BA=map_BA)
    else
      call dg_interpolation_galerkin_scalars_supermesh_free(old_fields, old_position, new_fields, new_position)
    end if

    deallocate(old_fields)
    deallocate(new_fields)

  end subroutine dg_interpolation_galerkin_state_supermesh_free
end module dg_interpolation_module
