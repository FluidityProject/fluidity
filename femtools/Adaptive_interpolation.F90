#include "fdebug.h"

module adaptive_interpolation_module

  use fldebug
  use vector_tools
  use quadrature
  use futils
  use elements
  use spud
  use quicksort
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
  use sparse_matrices_fields
  use solvers
  implicit none

  integer :: max_ai_degree=14

  private
  public :: adaptive_interpolation, max_ai_degree

  contains
    ! Interpolate a field onto new_positions in such a way that the L2 error squared
    ! of the projection is less than error_tolerance.
    subroutine adaptive_interpolation(old_field, old_positions, new_field, new_positions, error_tolerance, achieved_error, no_refinements)
      type(scalar_field), intent(in) :: old_field
      type(vector_field), intent(in) :: old_positions
      type(scalar_field), intent(inout) :: new_field
      type(vector_field), intent(inout) :: new_positions
      real, intent(in) :: error_tolerance
      real, intent(out) :: achieved_error
      integer, intent(out) :: no_refinements

      integer :: ele_A, ele_B
      real, dimension(ele_loc(new_positions, 1), ele_loc(new_positions, 1), ele_count(old_positions)) :: inversion_matrices_A

      integer :: dim
      type(ilist), dimension(ele_count(new_positions)) :: map_BA
      type(quadrature_type) :: supermesh_quad
      type(element_type) :: supermesh_shape

      type(vector_field) :: supermesh
      real :: ele_error, domain_volume, ele_volume, ele_tol
      type(element_type), dimension(max_ai_degree) :: shape_fns
      real, dimension(:), allocatable :: element_value
      integer :: p
      real, dimension(ele_loc(new_positions, 1), ele_loc(new_positions, 1)) :: inversion_matrix_B

      assert(old_positions%mesh%shape%degree == 1)
      assert(continuity(old_positions) >= 0)
      assert(continuity(new_positions) >= 0)
      assert(continuity(new_field) < 0)
      assert(error_tolerance > 0.0)

      dim = mesh_dim(new_positions)
      call intersector_set_dimension(dim)
      call zero(new_field)

      achieved_error = 0.0
      no_refinements = 0

      supermesh_quad = make_quadrature(vertices=ele_loc(new_positions, 1), dim=dim, degree=2*max_ai_degree, family=FAMILY_GM)
      supermesh_shape = make_element_shape(vertices=ele_loc(new_positions, 1), dim=dim, degree=1, quad=supermesh_quad)
      map_BA = intersection_finder(new_positions, old_positions)

      if (associated(new_field%mesh%region_ids)) then
        deallocate(new_field%mesh%region_ids)
      end if
      allocate(new_field%mesh%region_ids(ele_count(new_field)))
      new_field%mesh%region_ids = -10000

      do p=1,max_ai_degree
        shape_fns(p) = make_element_shape(vertices=ele_loc(new_positions, 1), dim=dim, degree=p, quad=supermesh_quad)
      end do
      allocate(element_value(shape_fns(max_ai_degree)%loc))

      domain_volume = 0.0
      do ele_A=1,ele_count(old_positions)
        call local_coords_matrix(old_positions, ele_A, inversion_matrices_A(:, :, ele_A))
        domain_volume = domain_volume + simplex_volume(old_positions, ele_A)
      end do

      do ele_B=1,ele_count(new_positions)
        call local_coords_matrix(new_positions, ele_B, inversion_matrix_B)
        p = element_degree(new_field, ele_B)

        call construct_supermesh(new_positions, ele_B, old_positions, map_BA(ele_B), supermesh_shape, supermesh)
        call galerkin_projection_ele_exact(old_field, inversion_matrices_A, shape_fns(p), new_positions, ele_B, &
                                         & inversion_matrix_B, supermesh, element_value)
        ele_volume = simplex_volume(new_positions, ele_B)
        ele_tol = error_tolerance * (ele_volume / domain_volume)
        ewrite(3,*) "ele_B: ", ele_B, "; ele_tol: ", ele_tol
        call compute_projection_error(old_field, old_positions, shape_fns(p), element_value, new_positions, ele_B, &
                                    & supermesh, inversion_matrices_A, inversion_matrix_B, ele_error)
        ewrite(3,*) "  p: ", p, "; ele_error: ", ele_error
        do while (ele_error > ele_tol)
          p = p + 1
          if (p > max_ai_degree) then
            ewrite(0,*) "Warning: cannot refine element ", ele_B, " to degree ", p
            exit
          end if
          no_refinements = no_refinements + 1
          call galerkin_projection_ele_exact(old_field, inversion_matrices_A, shape_fns(p), new_positions, ele_B, &
                                           & inversion_matrix_B, supermesh, element_value)
          call compute_projection_error(old_field, old_positions, shape_fns(p), element_value, new_positions, ele_B, &
                                      & supermesh, inversion_matrices_A, inversion_matrix_B, ele_error)
          ewrite(3,*) "  p: ", p, "; ele_error: ", ele_error
        end do

        achieved_error = achieved_error + ele_error
        new_field%mesh%region_ids(ele_B) = min(p, max_ai_degree)
        !write(0,*) "ele_B: ", ele_B, "; p: ", min(p, max_ai_degree - 1)

        call deallocate(supermesh)
      end do

      call deallocate(supermesh_shape)
      call deallocate(supermesh_quad)
      do ele_B=1,ele_count(new_positions)
        call deallocate(map_BA(ele_B))
      end do

    end subroutine adaptive_interpolation

    subroutine galerkin_projection_ele(old_field, inversion_matrices_A, shape_B, new_positions, ele_B, &
                                     & inversion_matrix_B, supermesh, element_value, stat)
      type(scalar_field), intent(in) :: old_field
      type(vector_field), intent(in) :: new_positions
      real, dimension(:, :, :), intent(in) :: inversion_matrices_A
      type(element_type), intent(in) :: shape_B
      type(vector_field), intent(in) :: supermesh
      integer, intent(in) :: ele_B
      real, dimension(:, :), intent(in) :: inversion_matrix_B
      integer, intent(out) :: stat
      real, dimension(:), intent(out) :: element_value

      integer :: ele_C, ele_A
      real :: vol_B, vols_C
      real, dimension(ele_loc(new_positions, ele_B), ele_loc(new_positions, ele_B)) :: inversion_matrix_A

      real, dimension(shape_B%loc, shape_B%loc) :: little_mass_matrix
      real, dimension(shape_B%loc, ele_loc(old_field, 1)) :: little_mixed_mass_matrix, little_mixed_mass_matrix_int
      real, dimension(shape_B%loc) :: little_rhs

      real, dimension(shape_B%ngi) :: detwei_B
      real, dimension(ele_ngi(supermesh, 1)) :: detwei_C

      type(element_type), pointer :: shape_A

      real, dimension(new_positions%dim, ele_ngi(supermesh, 1)) :: intersection_val_at_quad
      real, dimension(new_positions%dim+1, ele_ngi(supermesh, 1)) :: pos_at_quad_A, pos_at_quad_B
      real, dimension(shape_B%loc, ele_ngi(supermesh, 1)) :: basis_at_quad_B
      real, dimension(ele_loc(old_field, 1), ele_ngi(supermesh, 1)) :: basis_at_quad_A

      integer :: dim
      integer :: j, k, l

      dim = new_positions%dim

      assert(new_positions%mesh%shape%quadrature%degree == shape_B%quadrature%degree)
      call transform_to_physical(new_positions, ele_B, detwei=detwei_B)
      little_mass_matrix = shape_shape(shape_B, shape_B, detwei_B)

      vol_B = sum(detwei_B)
      vols_C = 0.0

      little_rhs = 0.0

      do ele_C=1,ele_count(supermesh)
        ele_A = ele_region_id(supermesh, ele_C)
        shape_A => ele_shape(old_field, ele_A)
        inversion_matrix_A = inversion_matrices_A(:, :, ele_A)

        intersection_val_at_quad = ele_val_at_quad(supermesh, ele_C)
        pos_at_quad_B(1:dim, :) = intersection_val_at_quad
        pos_at_quad_B(dim+1, :) = 1.0
        pos_at_quad_B = matmul(inversion_matrix_B, pos_at_quad_B)

        pos_at_quad_A(1:dim, :) = intersection_val_at_quad
        pos_at_quad_A(dim+1, :) = 1.0
        pos_at_quad_A = matmul(inversion_matrix_A, pos_at_quad_A)

        call transform_to_physical(supermesh, ele_C, detwei_C)
        vols_C = vols_C + sum(detwei_C)

        if (shape_B%degree==0) then
          basis_at_quad_B = 1.0
        elseif (shape_B%degree==1) then
          basis_at_quad_B = pos_at_quad_B 
        else
          do j=1,ele_ngi(supermesh, ele_C)
            basis_at_quad_B(:, j) = eval_shape(shape_B, pos_at_quad_B(:, j))
          end do
        end if

        if (shape_A%degree==0) then
          basis_at_quad_A = 1.0
        elseif (shape_A%degree==1) then
          basis_at_quad_A = pos_at_quad_A 
        else
          do j=1,ele_ngi(supermesh, ele_C)
            basis_at_quad_A(:, j) = eval_shape(shape_A, pos_at_quad_A(:, j))
          end do
        end if

        little_mixed_mass_matrix = 0.0
        little_mixed_mass_matrix_int = 0.0
        do j=1,ele_ngi(supermesh, ele_C)
          forall (k=1:shape_B%loc,l=1:ele_loc(old_field, ele_A))
            little_mixed_mass_matrix_int(k, l) = basis_at_quad_B(k, j) * basis_at_quad_A(l, j)
          end forall
          little_mixed_mass_matrix = little_mixed_mass_matrix + little_mixed_mass_matrix_int * detwei_C(j)
        end do

        little_rhs = little_rhs + matmul(little_mixed_mass_matrix, ele_val(old_field, ele_A))
    end do

    if (abs(vol_B - vols_C)/vol_B > 0.01) then
      stat = 1
      return
    else
      stat = 0
      call solve(little_mass_matrix, little_rhs)
      element_value(1:shape_B%loc) = little_rhs
    end if

    end subroutine galerkin_projection_ele

    subroutine galerkin_projection_ele_exact(old_field, inversion_matrices_A, shape_fn, new_positions, ele_B, &
                                           & inversion_matrix_B, supermesh, element_value)
      type(scalar_field), intent(in) :: old_field
      type(vector_field), intent(in) :: new_positions
      real, dimension(:, :, :), intent(in) :: inversion_matrices_A
      type(element_type), intent(in) :: shape_fn
      type(vector_field), intent(in) :: supermesh
      integer, intent(in) :: ele_B
      real, dimension(:, :), intent(in) :: inversion_matrix_B
      real, dimension(:), intent(out) :: element_value

      integer :: stat

      call galerkin_projection_ele(old_field, inversion_matrices_A, shape_fn, new_positions, ele_B, inversion_matrix_B, supermesh, element_value, stat)
      if (stat /= 0) then
        call intersector_set_exactness(.true.)
        call galerkin_projection_ele(old_field, inversion_matrices_A, shape_fn, new_positions, ele_B, inversion_matrix_B, supermesh, element_value, stat)
        call intersector_set_exactness(.false.)
      end if

    end subroutine galerkin_projection_ele_exact
end module adaptive_interpolation_module
