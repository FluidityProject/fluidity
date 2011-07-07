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

#include "fdebug.h"
module mangle_dirichlet_rows_module
  use sparse_matrices_fields
  use fields
  use field_options
  use boundary_conditions_from_options
  use boundary_conditions
  implicit none

  interface mangle_dirichlet_rows
    module procedure mangle_dirichlet_rows_csr_scalar, mangle_dirichlet_rows_scalar
  end interface

  interface set_inactive_rows
    module procedure set_inactive_rows_csr_scalar
  end interface

  interface compute_inactive_rows
    module procedure compute_inactive_rows_csr_scalar, compute_inactive_rows_csr_vector, compute_inactive_rows_block_csr_vector
  end interface


  contains

    subroutine mangle_dirichlet_rows_csr_scalar(matrix, field, keep_diag, rhs)
        type(csr_matrix), intent(inout) :: matrix
        type(scalar_field), intent(in) :: field
        logical, intent(in) :: keep_diag
        type(scalar_field), intent(inout), optional :: rhs

        integer :: i, j, node
        character(len=FIELD_NAME_LEN) :: bctype
        integer, dimension(:), pointer :: node_list
        type(scalar_field), pointer :: bc_field

        real,  dimension(:), pointer :: row
        real, pointer :: diag_ptr
        real :: diag

        do i=1,get_boundary_condition_count(field)
          call get_boundary_condition(field, i, type=bctype, surface_node_list=node_list)

          if (bctype /= "dirichlet") cycle

          if (present(rhs)) then
            bc_field => extract_surface_field(field, i, "value")
          end if

          do j=1,size(node_list)
            node = node_list(j)
            diag_ptr => diag_val_ptr(matrix, node)
            diag = diag_ptr
            row  => row_val_ptr(matrix, node)
            row = 0
            if (keep_diag) then
              diag_ptr = 1.0
            end if

            if (present(rhs)) then
              call set(rhs, node, node_val(bc_field, j))
            end if
          end do
        end do
    end subroutine mangle_dirichlet_rows_csr_scalar

    subroutine mangle_dirichlet_rows_scalar(scalar, bc_field)
        type(scalar_field), intent(inout) :: scalar
        type(scalar_field), intent(in) :: bc_field

        integer :: i, j, node
        character(len=FIELD_NAME_LEN) :: bctype
        integer, dimension(:), pointer :: node_list

        assert(node_count(scalar) == node_count(bc_field))

        do i=1,get_boundary_condition_count(bc_field)
          call get_boundary_condition(bc_field, i, type=bctype, surface_node_list=node_list)

          if (bctype /= "dirichlet") cycle

          do j=1,size(node_list)
            node = node_list(j)
            call set(scalar, node, 0.0)
          end do
        end do
    end subroutine mangle_dirichlet_rows_scalar

    subroutine set_inactive_rows_csr_scalar(matrix, field)
      ! Any nodes associated with Dirichlet BCs, we make them inactive
      type(csr_matrix), intent(inout) :: matrix
      type(scalar_field), intent(in) :: field

      integer :: i
      character(len=FIELD_NAME_LEN) :: bctype
      integer, dimension(:), pointer :: node_list

      do i=1,get_boundary_condition_count(field)
        call get_boundary_condition(field, i, type=bctype, surface_node_list=node_list)
        if (bctype /= "dirichlet") cycle
        call set_inactive(matrix, node_list)
      end do

    end subroutine set_inactive_rows_csr_scalar

    subroutine compute_inactive_rows_csr_scalar(x, A, rhs)
      type(scalar_field), intent(inout):: x
      type(csr_matrix), intent(in):: A
      type(scalar_field), intent(inout):: rhs
      logical, dimension(:), pointer:: inactive_mask
      integer:: i

      inactive_mask => get_inactive_mask(A)

      if (.not. associated(inactive_mask)) return

      do i=1, size(A,1)
        if (inactive_mask(i)) then

          ! we want [rhs_i - ((L+U)*x)_i]/A_ii, where L+U is off-diag part of A
          ! by setting x_i to zero before this is the same as [rhs_i- (A*x)_i]/A_ii
          call set(x, i, 0.0)
          call set(x, i, (node_val(rhs,i)- &
          dot_product(row_val_ptr(A,i),node_val(x, row_m_ptr(A,i))) &
          &  )/val(A,i,i) )

        end if
      end do
    end subroutine compute_inactive_rows_csr_scalar

    subroutine compute_inactive_rows_csr_vector(x, A, rhs)
      type(vector_field), intent(inout):: x
      type(csr_matrix), intent(in):: A
      type(vector_field), intent(inout):: rhs

      FLAbort("Not implemented yet")
    end subroutine compute_inactive_rows_csr_vector

    subroutine compute_inactive_rows_block_csr_vector(x, A, rhs)
      type(vector_field), intent(inout):: x
      type(block_csr_matrix), intent(in):: A
      type(vector_field), intent(inout):: rhs

      FLAbort("Not implemented yet")
    end subroutine compute_inactive_rows_block_csr_vector

end module mangle_dirichlet_rows_module
