!    Copyright (C) 2006 Imperial College London and others.                                                                                                                                                    
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineeringp
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
#include "confdefs.h"
#ifdef HAVE_ADJOINT

#include "fdebug.h"
#include "libadjoint/adj_fortran.h"
#endif
module libadjoint_data_callbacks
#ifdef HAVE_ADJOINT
use libadjoint
use boundary_conditions
use fields
use fields_manipulation
use sparse_tools

implicit none

  interface field_to_adj_vector
    module procedure scalar_field_to_adj_vector, vector_field_to_adj_vector, tensor_field_to_adj_vector
  end interface

  interface field_from_adj_vector
    module procedure scalar_field_from_adj_vector, vector_field_from_adj_vector, tensor_field_from_adj_vector
  end interface

  interface matrix_to_adj_matrix
    module procedure block_csr_matrix_to_adj_matrix, csr_matrix_to_adj_matrix
  end interface
  
  interface matrix_from_adj_matrix
    module procedure block_csr_matrix_from_adj_matrix, csr_matrix_from_adj_matrix
  end interface

  integer, parameter :: ADJ_SCALAR_FIELD=1, ADJ_VECTOR_FIELD=2, ADJ_TENSOR_FIELD=3, ADJ_CSR_MATRIX=11, ADJ_BLOCK_CSR_MATRIX=12, ADJ_MESH_TYPE=21
  integer, parameter :: ADJ_IDENTITY_MATRIX = 30
  ! things that can live in %flags
  integer, parameter :: ADJ_MATRIX_INVERTED = 1

  private 
  public :: field_to_adj_vector, field_from_adj_vector, adj_register_femtools_data_callbacks
  public :: matrix_to_adj_matrix, matrix_from_adj_matrix, mesh_type_to_adj_vector, mesh_type_from_adj_vector
  public :: femtools_vec_destroy_proc
  public :: femtools_mat_destroy_proc
  public :: ADJ_SCALAR_FIELD, ADJ_VECTOR_FIELD, ADJ_TENSOR_FIELD, ADJ_CSR_MATRIX, ADJ_BLOCK_CSR_MATRIX, ADJ_MESH_TYPE
  public :: ADJ_IDENTITY_MATRIX, ADJ_MATRIX_INVERTED
  public :: femtools_vec_axpy_proc

  contains

  subroutine adj_register_femtools_data_callbacks(adjointer)
    type(adj_adjointer), intent(inout) :: adjointer
    integer(kind=c_int) :: ierr

    ierr = adj_register_data_callback(adjointer, ADJ_VEC_DUPLICATE_CB, c_funloc(femtools_vec_duplicate_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_VEC_AXPY_CB, c_funloc(femtools_vec_axpy_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_VEC_DESTROY_CB, c_funloc(femtools_vec_destroy_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_VEC_SET_VALUES_CB, c_funloc(femtools_vec_setvalues_proc))  
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_VEC_DIVIDE_CB, c_funloc(femtools_vec_divide_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_VEC_GET_NORM_CB, c_funloc(femtools_vec_getnorm_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_VEC_GET_SIZE_CB, c_funloc(femtools_vec_get_size_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_VEC_DOT_PRODUCT_CB, c_funloc(femtools_vec_dot_product_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_VEC_SET_RANDOM_CB, c_funloc(femtools_vec_set_random_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_MAT_AXPY_CB, c_funloc(femtools_mat_axpy_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_MAT_DESTROY_CB, c_funloc(femtools_mat_destroy_proc))
    call adj_chkierr(ierr)
!    ierr = adj_register_data_callback(adjointer, ADJ_MAT_DUPLICATE_CB, c_funloc(femtools_mat_duplicate_proc))
!    call adj_chkierr(ierr)
end subroutine

  subroutine femtools_vec_duplicate_proc(x, newx) bind(c)
    ! Creates a new vector of the same type as an existing vector and set its entries to zero.
    use iso_c_binding
    use libadjoint_data_structures
    type(adj_vector), intent(in), value :: x
    type(adj_vector), intent(out) :: newx
    type(scalar_field), pointer :: x_scalar
    type(vector_field), pointer :: x_vector
    type(tensor_field), pointer :: x_tensor
    type(mesh_type), pointer :: x_mesh
    type(scalar_field) :: newx_scalar
    type(vector_field) :: newx_vector
    type(tensor_field) :: newx_tensor
    type(mesh_type) :: newx_mesh
    integer :: i
    character(len=30) :: bc_type
    character(len=FIELD_NAME_LEN) :: bc_name
    integer, dimension(:), pointer:: surface_element_list

    select case(x%klass)
      case(ADJ_SCALAR_FIELD)
        call c_f_pointer(x%ptr, x_scalar)
        call allocate(newx_scalar, x_scalar%mesh, trim(x_scalar%name))
        call zero(newx_scalar)
        do i=1,get_boundary_condition_count(x_scalar)
          call get_boundary_condition(x_scalar, i, type=bc_type, surface_element_list=surface_element_list, name=bc_name)
          if (trim(bc_type) == "dirichlet") then
            call add_boundary_condition_surface_elements(newx_scalar, bc_name, bc_type, surface_element_list)
          end if
        end do
        newx = field_to_adj_vector(newx_scalar)
        call deallocate(newx_scalar)
      case(ADJ_VECTOR_FIELD)
        call c_f_pointer(x%ptr, x_vector)
        call allocate(newx_vector, x_vector%dim, x_vector%mesh, trim(x_vector%name))
        call zero(newx_vector)
        do i=1,get_boundary_condition_count(x_vector)
          call get_boundary_condition(x_vector, i, type=bc_type, surface_element_list=surface_element_list, name=bc_name)
          if (trim(bc_type) == "dirichlet") then
            call add_boundary_condition_surface_elements(newx_vector, bc_name, bc_type, surface_element_list)
          end if
        end do
        newx = field_to_adj_vector(newx_vector)
        call deallocate(newx_vector)
      case(ADJ_TENSOR_FIELD)
        call c_f_pointer(x%ptr, x_tensor)
        call allocate(newx_tensor, x_tensor%mesh, trim(x_tensor%name), dim=x_tensor%dim)
        call zero(newx_tensor)
        newx = field_to_adj_vector(newx_tensor)
        call deallocate(newx_tensor)
      case default
        FLAbort("Unknown x%klass")
    end select

  end subroutine femtools_vec_duplicate_proc

  subroutine femtools_vec_axpy_proc(y, alpha, x) bind(c)
    use iso_c_binding
    use libadjoint_data_structures
    type(adj_vector), intent(inout) :: y
    adj_scalar_f, intent(in), value :: alpha
    type(adj_vector), intent(in), value :: x

    type(scalar_field), pointer :: x_scalar, y_scalar
    type(vector_field), pointer :: x_vector, y_vector
    type(tensor_field), pointer :: x_tensor, y_tensor

    assert(x%klass==y%klass)
    if (y%klass==ADJ_SCALAR_FIELD) then
      call c_f_pointer(y%ptr, y_scalar)
      call c_f_pointer(x%ptr, x_scalar)
      call addto(y_scalar, x_scalar, alpha)
    else if (y%klass==ADJ_VECTOR_FIELD) then
      call c_f_pointer(y%ptr, y_vector)
      call c_f_pointer(x%ptr, x_vector)
      call addto(y_vector, x_vector, alpha)
    else if (y%klass==ADJ_TENSOR_FIELD) then
      call c_f_pointer(y%ptr, y_tensor)
      call c_f_pointer(x%ptr, x_tensor)
      call addto(y_tensor, x_tensor, alpha)
    else
      FLAbort("adj_vector class not supported.")
    end if
  end subroutine femtools_vec_axpy_proc

  subroutine femtools_vec_destroy_proc(x) bind(c)
    ! Destroys a vector.
    use iso_c_binding
    type(adj_vector), intent(inout) :: x
    type(scalar_field), pointer :: scalar
    type(vector_field), pointer :: vector
    type(tensor_field), pointer :: tensor
    type(mesh_type), pointer :: mesh

    if (x%klass == ADJ_SCALAR_FIELD) then
      call c_f_pointer(x%ptr, scalar)
      call deallocate(scalar)
      deallocate(scalar)
    else if (x%klass == ADJ_VECTOR_FIELD) then
      call c_f_pointer(x%ptr, vector)
      call deallocate(vector)
      deallocate(vector)
    else if (x%klass == ADJ_TENSOR_FIELD) then
      call c_f_pointer(x%ptr, tensor)
      call deallocate(tensor)
      deallocate(tensor) 
    else if (x%klass == ADJ_MESH_TYPE) then
      call c_f_pointer(x%ptr, mesh)
      call deallocate(mesh)
      deallocate(mesh)
    else
      FLAbort("adj_vector class not supported.")
    end if
  end subroutine femtools_vec_destroy_proc

  subroutine femtools_vec_setvalues_proc(vec, scalars) bind(c)
    type(adj_vector), intent(inout) :: vec
    adj_scalar_f, dimension(*), intent(in) :: scalars
    type(scalar_field) :: scalar

    if (vec%klass == ADJ_SCALAR_FIELD) then
      call field_from_adj_vector(vec, scalar)
      scalar%val=scalars(1:size(scalar%val))
    else
      FLAbort("adj_vector class not supported (yet).")
    end if

  end subroutine femtools_vec_setvalues_proc

  subroutine femtools_vec_divide_proc(numerator, denominator) bind(c)
    type(adj_vector), intent(in), value :: numerator
    type(adj_vector), intent(inout) :: denominator

    type(scalar_field) :: scalar_numerator, scalar_denominator

    assert(numerator%klass == denominator%klass)
    if (numerator%klass == ADJ_SCALAR_FIELD) then
      call field_from_adj_vector(denominator, scalar_denominator)
      call field_from_adj_vector(numerator, scalar_numerator)
      if (scalar_numerator%field_type==FIELD_TYPE_CONSTANT) then
        assert(scalar_denominator%field_type==FIELD_TYPE_CONSTANT)
      end if
      assert(scalar_numerator%mesh==scalar_denominator%mesh)
      scalar_numerator%val=scalar_numerator%val/scalar_denominator%val
    else
      FLAbort("adj_vector class not supported.")
    end if
  end subroutine femtools_vec_divide_proc

  subroutine femtools_vec_getnorm_proc(vec, norm) bind(c)
    type(adj_vector), intent(in), value :: vec
    adj_scalar_f, intent(out) :: norm
    type(scalar_field), pointer :: scalar
    type(vector_field), pointer :: vector
    type(tensor_field), pointer :: tensor

    if (vec%klass == ADJ_SCALAR_FIELD) then
      call c_f_pointer(vec%ptr, scalar)
      norm = norm2(scalar%val)
    else if (vec%klass == ADJ_VECTOR_FIELD) then
      call c_f_pointer(vec%ptr, vector)
      norm = norm2(reshape(vector%val, (/size(vector%val)/)))
    else if (vec%klass == ADJ_TENSOR_FIELD) then
      call c_f_pointer(vec%ptr, tensor)
      norm = norm2(reshape(tensor%val, (/size(tensor%val)/)))
    else
      FLAbort("adj_vector class not supported.")
    end if
  end subroutine femtools_vec_getnorm_proc

  subroutine femtools_vec_dot_product_proc(x, y, val) bind(c)
    use libadjoint_data_structures
    type(adj_vector), intent(in), value :: x, y
    adj_scalar_f, intent(out) :: val
    type(scalar_field), pointer :: x_s, y_s
    type(vector_field), pointer :: x_v, y_v
    type(tensor_field), pointer :: x_t, y_t

    assert(x%klass == y%klass)
    if (x%klass == ADJ_SCALAR_FIELD) then
      call c_f_pointer(x%ptr, x_s)
      call c_f_pointer(y%ptr, y_s)
      val = dot_product(x_s%val, y_s%val)
    else if (x%klass == ADJ_VECTOR_FIELD) then
      call c_f_pointer(x%ptr, x_v)
      call c_f_pointer(y%ptr, y_v)
      val = sum(x_v%val * y_v%val)
    else if (x%klass == ADJ_TENSOR_FIELD) then
      call c_f_pointer(x%ptr, x_t)
      call c_f_pointer(y%ptr, y_t)
      val = sum(x_t%val * y_t%val)
    else
      FLAbort("adj_vector class not supported.")
    end if
  end subroutine femtools_vec_dot_product_proc

  subroutine femtools_vec_set_random_proc(vec) bind(c)
    use libadjoint_data_structures
    type(adj_vector), intent(inout) :: vec
    type(scalar_field), pointer :: scalar
    type(vector_field), pointer :: vector
    type(tensor_field), pointer :: tensor

    if (vec%klass == ADJ_SCALAR_FIELD) then
      call c_f_pointer(vec%ptr, scalar)
      call random_number(scalar%val)
    else if (vec%klass == ADJ_VECTOR_FIELD) then
      call c_f_pointer(vec%ptr, vector)
      call random_number(vector%val)
    else if (vec%klass == ADJ_TENSOR_FIELD) then
      call c_f_pointer(vec%ptr, tensor)
      call random_number(tensor%val)
    else
      FLAbort("adj_vector class not supported.")
    end if
  end subroutine femtools_vec_set_random_proc

  subroutine femtools_vec_get_size_proc(vec, sz) bind(c)
    type(adj_vector), intent(in), value :: vec
    integer(kind=c_int), intent(out) :: sz
    type(scalar_field), pointer :: scalar
    type(vector_field), pointer :: vector
    type(tensor_field), pointer :: tensor

    if (vec%klass == ADJ_SCALAR_FIELD) then
      call c_f_pointer(vec%ptr, scalar)
      sz = node_count(scalar)
    else if (vec%klass == ADJ_VECTOR_FIELD) then
      call c_f_pointer(vec%ptr, vector)
      sz = node_count(vector) * vector%dim
    else if (vec%klass == ADJ_TENSOR_FIELD) then
      call c_f_pointer(vec%ptr, tensor)
      sz = node_count(tensor) * tensor%dim(1) * tensor%dim(2)
    else
      FLAbort("adj_vector class not supported.")
    end if
  end subroutine femtools_vec_get_size_proc

  subroutine femtools_mat_axpy_proc(Y, alpha, X) bind(c)  
    ! Computes Y = alpha*X + Y.
    use iso_c_binding
    use libadjoint_data_structures
    type(adj_matrix), intent(inout) :: Y
    adj_scalar_f, intent(in), value :: alpha
    type(adj_matrix), intent(in), value :: X

    type(csr_matrix), pointer :: csr_y, csr_x
    type(block_csr_matrix), pointer :: block_csr_x, block_csr_y
    integer :: i 

    assert(x%klass==y%klass)
    if (x%klass==ADJ_CSR_MATRIX) then
      call c_f_pointer(x%ptr, csr_x)
      call c_f_pointer(y%ptr, csr_y)
      call addto(csr_y, csr_x, alpha)
    else if (x%klass==ADJ_BLOCK_CSR_MATRIX) then
      call c_f_pointer(x%ptr, block_csr_x)
      call c_f_pointer(y%ptr, block_csr_y)
      assert(all(block_csr_x%blocks==block_csr_y%blocks))
      call addto(block_csr_y,  (/ ( i, i=1,block_csr_x%blocks(1) ) /),  (/ ( i, i=1,block_csr_x%blocks(2) ) /), block_csr_x, alpha)
    else
      FLAbort("adj_matrix class not supported.")
    end if

  end subroutine femtools_mat_axpy_proc


  subroutine femtools_mat_destroy_proc(mat) bind(c)
    ! Frees space taken by a matrix.
    use iso_c_binding
    type(adj_matrix), intent(inout) :: mat
    type(csr_matrix), pointer :: csr_mat
    type(block_csr_matrix), pointer :: block_csr_mat

    if (mat%klass==ADJ_CSR_MATRIX) then
      call c_f_pointer(mat%ptr, csr_mat)
      call deallocate(csr_mat)
      deallocate(csr_mat)
    else if (mat%klass==ADJ_BLOCK_CSR_MATRIX) then
      call c_f_pointer(mat%ptr, block_csr_mat)
      call deallocate(block_csr_mat)
      deallocate(block_csr_mat)
    else
      FLAbort("adj_matrix class not supported.")
    end if
  end subroutine femtools_mat_destroy_proc

 subroutine scalar_field_from_adj_vector(input, output) 
    type(adj_vector), intent(in) :: input
    type(scalar_field), intent(out) :: output
    type(scalar_field), pointer :: tmp

    assert(input%klass==ADJ_SCALAR_FIELD)
    call c_f_pointer(input%ptr, tmp)
    output = tmp
  end subroutine scalar_field_from_adj_vector

  function scalar_field_to_adj_vector(input) result(output)
    type(scalar_field), intent(in), target :: input
    type(scalar_field), pointer :: input_ptr
    type(adj_vector) :: output

    call incref(input)
    output%klass = ADJ_SCALAR_FIELD
    output%flags = 0 
    allocate(input_ptr)
    input_ptr = input
    output%ptr = c_loc(input_ptr)
  end function scalar_field_to_adj_vector

  subroutine vector_field_from_adj_vector(input, output) 
    type(adj_vector), intent(in) :: input
    type(vector_field), intent(out) :: output
    type(vector_field), pointer :: tmp

    assert(input%klass==ADJ_VECTOR_FIELD)
    call c_f_pointer(input%ptr, tmp)
    output = tmp
  end subroutine vector_field_from_adj_vector

  function vector_field_to_adj_vector(input) result(output)
    type(vector_field), intent(in), target :: input
    type(vector_field), pointer :: input_ptr
    type(adj_vector) :: output

    call incref(input)
    output%klass = ADJ_VECTOR_FIELD
    output%flags = 0
    allocate(input_ptr)
    input_ptr = input
    output%ptr = c_loc(input_ptr)
  end function vector_field_to_adj_vector


  subroutine tensor_field_from_adj_vector(input, output) 
    type(adj_vector), intent(in) :: input
    type(tensor_field), intent(out) :: output
    type(tensor_field), pointer :: tmp

    assert(input%klass==ADJ_TENSOR_FIELD)
    call c_f_pointer(input%ptr, tmp)
    output = tmp
  end subroutine tensor_field_from_adj_vector

  function tensor_field_to_adj_vector(input) result(output)
    type(tensor_field), intent(in), target :: input
    type(tensor_field), pointer :: input_ptr
    type(adj_vector) :: output

    call incref(input)
    output%klass = ADJ_TENSOR_FIELD
    output%flags = 0
    allocate(input_ptr)
    input_ptr = input
    output%ptr = c_loc(input_ptr)
  end function tensor_field_to_adj_vector

  function csr_matrix_to_adj_matrix(input) result(output)
    type(csr_matrix), intent(in), target :: input
    type(csr_matrix), pointer :: input_ptr
    type(adj_matrix) :: output

    call incref(input)
    output%klass = ADJ_CSR_MATRIX 
    output%flags = 0
    allocate(input_ptr)
    input_ptr = input
    output%ptr = c_loc(input_ptr)
  end function csr_matrix_to_adj_matrix

  subroutine csr_matrix_from_adj_matrix(input, output) 
    type(adj_matrix), intent(in) :: input
    type(csr_matrix), intent(out) :: output
    type(csr_matrix), pointer :: tmp

    assert(input%klass==ADJ_CSR_MATRIX)
    call c_f_pointer(input%ptr, tmp)
    output = tmp
  end subroutine csr_matrix_from_adj_matrix

  function block_csr_matrix_to_adj_matrix(input) result(output)
    type(block_csr_matrix), intent(in), target :: input
    type(block_csr_matrix), pointer :: input_ptr
    type(adj_matrix) :: output

    call incref(input)
    output%klass = ADJ_BLOCK_CSR_MATRIX
    output%flags = 0
    allocate(input_ptr)
    input_ptr = input
    output%ptr = c_loc(input_ptr)
  end function block_csr_matrix_to_adj_matrix

  subroutine block_csr_matrix_from_adj_matrix(input, output) 
    type(adj_matrix), intent(in) :: input
    type(block_csr_matrix), intent(out) :: output
    type(block_csr_matrix), pointer :: tmp

    assert(input%klass==ADJ_BLOCK_CSR_MATRIX)
    call c_f_pointer(input%ptr, tmp)
    output = tmp
  end subroutine block_csr_matrix_from_adj_matrix

  function mesh_type_to_adj_vector(input) result(output)
    type(mesh_type), intent(in), target :: input
    type(mesh_type), pointer :: input_ptr
    type(adj_vector) :: output

    call incref(input)
    output%klass = ADJ_MESH_TYPE
    output%flags = 0
    allocate(input_ptr)
    input_ptr = input
    output%ptr = c_loc(input_ptr)
  end function mesh_type_to_adj_vector

  subroutine mesh_type_from_adj_vector(input, output) 
    type(adj_vector), intent(in) :: input
    type(mesh_type), intent(out) :: output
    type(mesh_type), pointer :: tmp

    assert(input%klass==ADJ_MESH_TYPE)
    call c_f_pointer(input%ptr, tmp)
    output = tmp
  end subroutine mesh_type_from_adj_vector
#endif

end module libadjoint_data_callbacks

