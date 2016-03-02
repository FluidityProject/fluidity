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
#include "fdebug.h"
module sparse_matrices_fields
use iso_c_binding
use fldebug
use sparse_tools
use fields
use sparse_tools_petsc
use c_interfaces
implicit none

  interface mult
     module procedure csr_mult_scalar, csr_mult_scalar_vector,&
          csr_mult_vector_scalar, csr_mult_vector_vector, csr_mult_vector
  end interface
  
  interface mult_addto
     module procedure block_csr_mult_addto_vector, csr_mult_addto_scalar
  end interface

  interface mult_T
     module procedure csr_mult_T_scalar, csr_mult_T_vector_scalar, csr_mult_T_vector_vector
  end interface

  interface mult_T_addto
     module procedure csr_mult_T_addto_scalar
  end interface

  interface mult_diag
     module procedure csr_diag_mult_scalar, csr_diag_mult_scalar_v
  end interface

  interface addto_diag
     module procedure csr_diag_addto_scalar, block_csr_diag_addto_scalar, &
       petsc_csr_diag_addto_scalar, petsc_csr_diag_addto_vector
  end interface

  interface extract_diagonal
     module procedure block_csr_extract_diagonal
  end interface

private

public :: mult, mult_addto, mult_T, mult_T_addto, mult_diag, addto_diag,&
          extract_diagonal, mult_div_tensorinvscalar_div_t,&
	  mult_div_tensorinvscalar_vector, mult_div_invscalar_div_t,&
	  mult_div_vector_div_t, mult_div_invvector_div_t
  
contains

  subroutine csr_mult_scalar(x, A, b)
    !!< Calculate x=A*b
    type(scalar_field), intent(inout), target :: x
    type(csr_matrix), intent(in) :: A
    type(scalar_field), intent(in), target :: b
    real, dimension(:), allocatable :: tmp

    if (compare_pointers(c_loc(x), c_loc(b))) then
      FLAbort("You can't pass the same field in for x and b.")
    end if

    select case(b%field_type)
    case(FIELD_TYPE_NORMAL)
      call mult(x%val, A, b%val)
    case(FIELD_TYPE_CONSTANT)
      allocate(tmp(size(x%val)))
      tmp=b%val(1)
      call mult(x%val, A, tmp)
      deallocate(tmp)
    end select
    
  end subroutine csr_mult_scalar

  subroutine csr_mult_addto_scalar(x, A, b)
    !!< Replace x with x+A*b
    type(scalar_field), intent(inout) :: x
    type(csr_matrix), intent(in) :: A
    type(scalar_field), intent(in) :: b
    real, dimension(:), allocatable :: tmp

    select case(b%field_type)
    case(FIELD_TYPE_NORMAL)
      call mult_addto(x%val, A, b%val)
    case(FIELD_TYPE_CONSTANT)
      allocate(tmp(size(x%val)))
      tmp=b%val(1)
      call mult_addto(x%val, A, tmp)
      deallocate(tmp)
    end select
    
  end subroutine csr_mult_addto_scalar

  subroutine csr_mult_vector(x, A, b)
    !!< Calculate x=A*b
    type(vector_field), intent(inout) :: x
    type(csr_matrix), intent(in) :: A
    type(vector_field), intent(in) :: b
    
    integer :: i
    type(scalar_field) :: x_comp, b_comp

    assert(x%dim==b%dim)
    do i = 1, b%dim
      x_comp = extract_scalar_field(x, i)
      b_comp = extract_scalar_field(b, i)
      call mult(x_comp, A, b_comp)
    end do

  end subroutine csr_mult_vector

  subroutine block_csr_mult_addto_vector(x, A, b)
    !!< Calculate x=A*b
    type(vector_field), intent(inout) :: x
    type(block_csr_matrix), intent(in) :: A
    type(vector_field), intent(in) :: b
    
    integer :: i, j
    type(scalar_field) :: x_comp, b_comp
    type(csr_matrix) :: A_block

    assert(x%dim==b%dim)
    call zero(x)
    do i = 1, b%dim
       b_comp = extract_scalar_field(b, i)
       if(A%diagonal) then
          x_comp = extract_scalar_field(x, i)
          A_block = block(A,i,i)
          call mult_addto(x_comp, A_block, b_comp)
       else
          do j = 1, b%dim
             x_comp = extract_scalar_field(x, j)
             A_block = block(A,i,j)
             call mult_addto(x_comp, A_block, b_comp)
          end do
       end if
    end do
  end subroutine block_csr_mult_addto_vector

  subroutine csr_diag_mult_scalar(A, b)
    !!< Calculate A_{i,i}=A_{i,i}*b_{i}
    !!<           A_{i,j}=A_{i,j}, i/=j
    type(csr_matrix), intent(inout) :: A
    type(scalar_field), intent(in) :: b
    real, dimension(:), allocatable :: tmp
    integer :: i

    assert(size(A,1)==node_count(b))
    assert(size(A,2)==node_count(b))

    select case(b%field_type)
    case(FIELD_TYPE_NORMAL)
      do i = 1, size(A,1)
        call set(A, i, i, val(A, i, i)*b%val(i))
      end do
    case(FIELD_TYPE_CONSTANT)
      allocate(tmp(node_count(b)))
      tmp=b%val(1)
      do i = 1, size(A,1)
        call set(A, i, i, val(A, i, i)*tmp(i))
      end do
      deallocate(tmp)
    end select
    
  end subroutine csr_diag_mult_scalar

  subroutine csr_diag_mult_scalar_v(A, b)
    !!< Calculate A_{i,i}=A_{i,i}*b_{i}
    !!<           A_{i,j}=A_{i,j}, i/=j
    type(csr_matrix), intent(inout) :: A
    real, dimension(:), intent(in) :: b
    integer :: i

    assert(size(A,1)==size(b))
    assert(size(A,2)==size(b))

    do i = 1, size(A,1)
       call set(A, i, i, val(A, i, i)*b(i))
    end do

  end subroutine csr_diag_mult_scalar_v

  subroutine csr_diag_addto_scalar(A, b, scale)
    !!< Calculate X_{i,i}=A_{i,i}*b_{i}
    !!<           X_{i,j}=A_{i,j}, i/=j
    type(csr_matrix), intent(inout) :: A
    type(scalar_field), intent(in) :: b
    real, optional, intent(in) :: scale
    real, dimension(:), allocatable :: tmp
    integer :: i

    assert(size(A,1)==node_count(b))
    assert(size(A,2)==node_count(b))

    select case(b%field_type)
    case(FIELD_TYPE_NORMAL)
      if(present(scale)) then
        do i = 1, size(A,1)
          call addto(A, i, i, b%val(i)*scale)
        end do
      else
        do i = 1, size(A,1)
          call addto(A, i, i, b%val(i))
        end do
      end if
    case(FIELD_TYPE_CONSTANT)
      allocate(tmp(node_count(b)))
      tmp=b%val(1)
      if(present(scale)) then
        do i = 1, size(A,1)
          call addto(A, i, i, tmp(i)*scale)
        end do
      else
        do i = 1, size(A,1)
          call addto(A, i, i, tmp(i))
        end do
      end if
      deallocate(tmp)
    end select
    
  end subroutine csr_diag_addto_scalar

  subroutine block_csr_diag_addto_scalar(A, b, scale)
    !!< Calculate X_{i,i}=A_{i,i}*b_{i}
    !!<           X_{i,j}=A_{i,j}, i/=j
    type(block_csr_matrix), intent(inout) :: A
    type(scalar_field), intent(in) :: b
    real, optional, intent(in) :: scale
    real, dimension(:), allocatable :: tmp
    integer :: i, j
    real :: l_scale
    type(csr_matrix) :: block_A

    if(present(scale)) then
      l_scale = scale
    else
      l_scale = 1.0
    end if

    assert(A%blocks(1)==A%blocks(2))

    do j = 1, A%blocks(1)

      block_A = block(A, j, j)

      assert(size(block_A,1)==node_count(b))
      assert(size(block_A,2)==node_count(b))

      select case(b%field_type)
      case(FIELD_TYPE_NORMAL)
        do i = 1, size(block_A,1)
          call addto(block_A, i, i, b%val(i)*l_scale)
        end do
      case(FIELD_TYPE_CONSTANT)
        allocate(tmp(node_count(b)))
        tmp=b%val(1)
        do i = 1, size(block_A,1)
          call addto(block_A, i, i, tmp(i)*l_scale)
        end do
        deallocate(tmp)
      end select

    end do

  end subroutine block_csr_diag_addto_scalar

  subroutine petsc_csr_diag_addto_scalar(A, b, scale)
    !!< Calculate X_{i,i}=A_{i,i}*b_{i}
    !!<           X_{i,j}=A_{i,j}, i/=j
    type(petsc_csr_matrix), intent(inout) :: A
    type(scalar_field), intent(in) :: b
    real, optional, intent(in) :: scale
    real, dimension(:), allocatable :: tmp
    integer :: i, j
    real :: l_scale

    if(present(scale)) then
      l_scale = scale
    else
      l_scale = 1.0
    end if

    assert( blocks(A,1)==blocks(A,2) )
    assert( block_size(A,1)==node_count(b) )
    assert( block_size(A,2)==node_count(b) )

    do j = 1, blocks(A,1)

      select case(b%field_type)
      case(FIELD_TYPE_NORMAL)
        do i = 1, block_size(A,1)
          call addto(A, j, j, i, i, b%val(i)*l_scale)
        end do
      case(FIELD_TYPE_CONSTANT)
        allocate(tmp(node_count(b)))
        tmp=b%val(1)
        do i = 1, block_size(A,1)
          call addto(A, j, j, i, i, tmp(i)*l_scale)
        end do
        deallocate(tmp)
      end select

    end do

  end subroutine petsc_csr_diag_addto_scalar

  subroutine petsc_csr_diag_addto_vector(A, b, scale)
    !!< Calculate X_{i,i}=A_{i,i}*b_{i}
    !!<           X_{i,j}=A_{i,j}, i/=j
    type(petsc_csr_matrix), intent(inout) :: A
    type(vector_field), intent(in) :: b
    real, optional, intent(in) :: scale
    real, dimension(:), allocatable :: tmp
    integer :: i, j
    real :: l_scale

    if(present(scale)) then
      l_scale = scale
    else
      l_scale = 1.0
    end if

    assert( blocks(A,1)==blocks(A,2) )
    assert( block_size(A,1)==node_count(b) )
    assert( block_size(A,2)==node_count(b) )

    do j = 1, blocks(A,1)

      select case(b%field_type)
      case(FIELD_TYPE_NORMAL)
        do i = 1, block_size(A,1)
          call addto(A, j, j, i, i, b%val(j,i)*l_scale)
        end do
      case(FIELD_TYPE_CONSTANT)
        allocate(tmp(node_count(b)))
        tmp=b%val(j,:)
        do i = 1, block_size(A,1)
          call addto(A, j, j, i, i, tmp(i)*l_scale)
        end do
        deallocate(tmp)
      end select

    end do

  end subroutine petsc_csr_diag_addto_vector
  
  subroutine block_csr_extract_diagonal(A,diagonal)
    !!< Extracts diagonal components of a block_csr matrix.
    !!< The vector field diagonal needs to be allocated before the call.

    type(block_csr_matrix), intent(in) :: A
    type(vector_field), intent(inout) :: diagonal
    integer i, j
    
    assert( diagonal%dim==blocks(A,1) )
    assert( node_count(diagonal)==block_size(A,1))
    assert( block_size(A,1)==block_size(A,2))
    assert( blocks(A,1)==blocks(A,2))
    
    if (associated(A%sparsity%centrm)) then
       do i=1, blocks(A,1)
          call set_all(diagonal, i, A%val(i,i)%ptr(A%sparsity%centrm))
       end do
    else
       do i=1, blocks(A,1)
          do j=1, block_size(A,1)
             call set(diagonal, i, j, val(A, i, i, j, j))
          end do
       end do       
    end if

  end subroutine block_csr_extract_diagonal

  subroutine csr_mult_scalar_vector(x, A, b)
    !!< Calculate x=A*b Where b is a vector field, A is a 1*dim block
    !!< block_csr_matrix and x is a scalar field.
    type(scalar_field), intent(inout) :: x
    type(block_csr_matrix), intent(in) :: A
    type(vector_field), intent(in) :: b

    real, dimension(:), allocatable :: tmpb, tmpx
    integer :: dim

    assert(all(A%blocks==(/1,b%dim/)))

    allocate(tmpx(size(x%val)))
    call zero(x)
    
    do dim=1,b%dim
       select case(b%field_type)
       case(FIELD_TYPE_NORMAL)
          call mult(tmpx, block(A,1,dim), b%val(dim,:))
       case(FIELD_TYPE_CONSTANT)
          allocate(tmpb(size(x%val)))
          tmpb=b%val(dim,1)
          call mult(tmpx, block(A,1,dim), tmpb)
          deallocate(tmpb)
       end select
       x%val=x%val+tmpx
    end do

  end subroutine csr_mult_scalar_vector

  subroutine csr_mult_vector_scalar(x, A, b)
    !!< Calculate x=A*b Where b is a scalar field, A is a dim*1 block
    !!< block_csr_matrix and x is a vector field.
    type(vector_field), intent(inout) :: x
    type(block_csr_matrix), intent(in) :: A
    type(scalar_field), intent(in) :: b

    real, dimension(:), allocatable :: tmpb
    integer :: dim

    assert(all(A%blocks==(/x%dim,1/)))

    call zero(x)
    
    do dim=1,x%dim
       select case(b%field_type)
       case(FIELD_TYPE_NORMAL)
          call mult(x%val(dim,:), block(A,dim,1), b%val)
       case(FIELD_TYPE_CONSTANT)
          allocate(tmpb(size(x%val(dim,:))))
          tmpb=b%val
          call mult(x%val(dim,:), block(A,dim,1), tmpb)
          deallocate(tmpb)
       end select
    end do

  end subroutine csr_mult_vector_scalar

  subroutine csr_mult_vector_vector(x, A, b)
    !!< Calculate x=A*b Where b is a vector field, A is a dim_x*dim_b block
    !!< block_csr_matrix and x is a vector field.
    type(vector_field), intent(inout) :: x
    type(block_csr_matrix), intent(in) :: A
    type(vector_field), intent(in) :: b

    real, dimension(:), allocatable :: tmpb, tmpx
    integer :: dim_x, dim_b

    assert(all(A%blocks==(/x%dim,b%dim/)))

    allocate(tmpx(size(x%val(1,:))))
    call zero(x)
    
    do dim_x=1,x%dim
       do dim_b=1,b%dim
          if (A%diagonal .and. dim_b/=dim_x) cycle
          
          select case(b%field_type)
          case(FIELD_TYPE_NORMAL)
             call mult(tmpx, block(A,dim_x,dim_b), b%val(dim_b,:))
          case(FIELD_TYPE_CONSTANT)
             allocate(tmpb(size(x%val)))
             tmpb=b%val(dim_b,dim_x)
             call mult(tmpx, block(A,dim_x,dim_b), tmpb)
             deallocate(tmpb)
          end select
          x%val(dim_x,:)=x%val(dim_x,:)+tmpx
       end do
    end do

  end subroutine csr_mult_vector_vector

  subroutine csr_mult_T_scalar(x, A, b)
    !!< Calculate x=A^T*b
    type(scalar_field), intent(inout) :: x
    type(csr_matrix), intent(in) :: A
    type(scalar_field), intent(in) :: b
    real, dimension(:), allocatable :: tmp

    select case(b%field_type)
    case(FIELD_TYPE_NORMAL)
      call mult_T(x%val, A, b%val)
    case(FIELD_TYPE_CONSTANT)
      allocate(tmp(size(x%val)))
      tmp=b%val(1)
      call mult_T(x%val, A, tmp)
      deallocate(tmp)
    end select
    
  end subroutine csr_mult_T_scalar

  subroutine csr_mult_T_addto_scalar(x, A, b)
    !!< Calculate x=A^T*b
    type(scalar_field), intent(inout) :: x
    type(csr_matrix), intent(in) :: A
    type(scalar_field), intent(in) :: b
    real, dimension(:), allocatable :: tmp

    select case(b%field_type)
    case(FIELD_TYPE_NORMAL)
      call mult_T_addto(x%val, A, b%val)
    case(FIELD_TYPE_CONSTANT)
      allocate(tmp(size(x%val)))
      tmp=b%val(1)
      call mult_T_addto(x%val, A, tmp)
      deallocate(tmp)
    end select

  end subroutine csr_mult_T_addto_scalar

  subroutine csr_mult_T_vector_scalar(x, A, b)
    !!< Calculate x=A^T*b Where b is a scalar field, A is a 1*dim block
    !!< block_csr_matrix and x is a vector field.
    type(vector_field), intent(inout) :: x
    type(block_csr_matrix), intent(in) :: A
    type(scalar_field), intent(in) :: b

    real, dimension(:), allocatable :: tmpb
    integer :: dim

    assert(all(A%blocks==(/1,x%dim/)))

    call zero(x)
    
    do dim=1,x%dim
       select case(b%field_type)
       case(FIELD_TYPE_NORMAL)
          call mult_T(x%val(dim,:), block(A,1,dim), b%val)
       case(FIELD_TYPE_CONSTANT)
          allocate(tmpb(node_count(b)))
          tmpb=b%val(1)
          call mult_T(x%val(dim,:), block(A,1,dim), tmpb)
          deallocate(tmpb)
       end select
    end do

  end subroutine csr_mult_T_vector_scalar

  subroutine csr_mult_T_vector_vector(x, A, b)
    !!< Calculate x=A^T*b, where b is a vector fields, A is a dim_b*dim_x block
    !!< block_csr_matrix and x is a vector field.
    type(vector_field), intent(inout) :: x
    type(block_csr_matrix), intent(in) :: A
    type(vector_field), intent(in) :: b

    real, dimension(:), allocatable :: tmpb, tmpx
    integer :: dim_x, dim_b

    assert(all(A%blocks==(/b%dim,x%dim/)))

    allocate(tmpx(size(x%val(1,:))))
    call zero(x)
    
    do dim_x=1,x%dim
       do dim_b=1,b%dim
          if (A%diagonal .and. dim_b/=dim_x) cycle
          
          select case(b%field_type)
          case(FIELD_TYPE_NORMAL)
             call mult_T(tmpx, block(A,dim_b,dim_x), b%val(dim_b,:))
          case(FIELD_TYPE_CONSTANT)
             allocate(tmpb(size(x%val)))
             tmpb=b%val(dim_b,dim_x)
             call mult_T(tmpx, block(A,dim_b,dim_x), tmpb)
             deallocate(tmpb)
          end select
          x%val(dim_x,:)=x%val(dim_x,:)+tmpx
       end do
    end do

  end subroutine csr_mult_T_vector_vector
  
  subroutine mult_div_vector_div_T(product, matrix1, vfield, matrix2)
    !!< Perform the matrix multiplication:
    !!< 
    !!<   product = matrix1*diag(vector)*matrix2^T
    !!< 
    !!< Only works on block csr matrices with monotonic row entries in colm
    type(csr_matrix) :: product
    type(block_csr_matrix), intent(in) :: matrix1, matrix2
    type(vector_field), intent(in) :: vfield

    integer, dimension(:), pointer :: row, col, row_product
    type(real_vector), dimension(blocks(matrix1,2)) :: row_val, col_val
    integer :: i,j,k1,k2,jcol,dim,ndim
    real :: entry0
    integer :: nentry0

    ewrite(1,*) 'Entering mult_div_invvector_div_T'

    call zero(product)

    nentry0 = 0
    assert(size(matrix1,2)==size(matrix2,2))
    assert(size(matrix1,1)==size(product,1))
    assert(size(matrix2,1)==size(product,2))
    assert(blocks(matrix1,1)==blocks(matrix2,1))
    assert(blocks(matrix1,2)==blocks(matrix2,2))
    assert(blocks(matrix1,2)==vfield%dim)
    assert(blocks(matrix1,1)==1)

    ndim = blocks(matrix1,2)

    if(.not.matrix1%sparsity%sorted_rows.or..not.matrix2%sparsity%sorted_rows) then
      FLAbort("mult_div_invvector_div_T assumes sorted rows")
    end if

    ! multiplication M_ij=A_ik * D_k * B_jk

    do i=1, size(matrix1%sparsity,1)
       row=>row_m_ptr(matrix1, i)

       do dim = 1, ndim
         row_val(dim)%ptr=>row_val_ptr(matrix1, 1, dim, i)
       end do

       row_product=>row_m_ptr(product, i)
       do jcol = 1, size(row_product)
         j = row_product(jcol)

         col=>row_m_ptr(matrix2, j)
         do dim = 1, ndim
           col_val(dim)%ptr=>row_val_ptr(matrix2, 1, dim, j)
         end do

         ! to compute entry M_ij find common column indices k in rows
         ! A_ik and B_jk - this is done by walking through them together
         ! from left to right, where we're using that both are stored
         ! in sorted order
         
         entry0=0.0

         k1 = 1
         k2 = 1
         do while (k1<=size(row) .and. k2<=size(col))
           if(row(k1)<col(k2)) then
             k1 = k1 + 1
           else if(row(k1)==col(k2)) then
             do dim = 1, ndim
               entry0=entry0+row_val(dim)%ptr(k1)* &
                 col_val(dim)%ptr(k2)*node_val(vfield, dim, row(k1))
             end do
             k1 = k1 + 1
             k2 = k2 + 1
           else
             k2 = k2 + 1
           end if
         end do
         nentry0 = nentry0 + 1
         product%val(nentry0) = entry0
       end do
    end do

  end subroutine mult_div_vector_div_T
  
  subroutine mult_div_invscalar_div_T(product, matrix1, sfield, matrix2)
    !!< Perform the matrix multiplication:
    !!< 
    !!<   product = matrix1*diag(1./scalar)*matrix2^T
    !!< 
    !!< Only works on csr matrices with monotonic row entries in colm
    type(csr_matrix), intent(inout) :: product
    type(csr_matrix), intent(in) :: matrix1
    type(scalar_field), intent(in) :: sfield
    type(csr_matrix), intent(in) :: matrix2
    
    integer, dimension(:), pointer :: row, col, row_product
    real, dimension(:), pointer :: row_val, col_val
    integer :: i,j,k1,k2,jcol
    real :: entry0
    integer :: nentry0
    logical :: addflag

    ewrite(1,*) 'Entering mult_div_invscalar_div_T'

    call zero(product)

    nentry0 = 0
    assert(size(matrix1,2)==size(matrix2,2))

    if(.not.matrix1%sparsity%sorted_rows.or..not.matrix2%sparsity%sorted_rows) then
      FLAbort("mult_div_invscalar_div_T assumes sorted rows")
    end if

    do i=1, size(matrix1%sparsity,1)
       row=>row_m_ptr(matrix1, i)

       if(size(row)>0) then
          row_val=>row_val_ptr(matrix1, i)

          row_product=>row_m_ptr(product, i)
          do jcol = 1, size(row_product)
             j = row_product(jcol)

             col=>row_m_ptr(matrix2, j)
             col_val=>row_val_ptr(matrix2, j)

             if(size(col)>0) then
                entry0=0.0

                addflag = .false.

                k1 = 1
                k2 = 1
                do
                   if((k1.gt.size(row)).or.(k2.gt.size(col))) exit
                   if(row(k1)<col(k2)) then
                      k1 = k1 + 1
                   else
                      if(row(k1)==col(k2)) then
                         ! Note the transpose in the second val call.
                         entry0=entry0+row_val(k1)* &
                              col_val(k2)/node_val(sfield, row(k1))
                         addflag = .true.
                         k1 = k1 + 1
                         k2 = k2 + 1
                      else
                         k2 = k2 + 1
                      end if
                   end if
                end do
                if(addflag) then
                   nentry0 = nentry0 + 1
                   product%val(nentry0) = entry0
                end if
             end if
          end do
       end if
    end do

  end subroutine mult_div_invscalar_div_T

  subroutine mult_div_invvector_div_T(product, matrix1, vfield, matrix2)
    !!< Perform the matrix multiplication:
    !!< 
    !!<   product = matrix1*diag(1./vector)*matrix2^T
    !!< 
    !!< Only works on block csr matrices with monotonic row entries in colm
    type(csr_matrix) :: product
    type(block_csr_matrix), intent(in) :: matrix1, matrix2
    type(vector_field), intent(in) :: vfield

    integer, dimension(:), pointer :: row, col, row_product
    type(real_vector), dimension(blocks(matrix1,2)) :: row_val, col_val
    integer :: i,j,k1,k2,jcol,dim,ndim
    real :: entry0
    integer :: nentry0
    logical :: addflag

    ewrite(1,*) 'Entering mult_div_invvector_div_T'

    call zero(product)

    nentry0 = 0
    assert(size(matrix1,2)==size(matrix2,2))
    assert(blocks(matrix1,1)==blocks(matrix2,1))
    assert(blocks(matrix1,2)==blocks(matrix2,2))
    assert(blocks(matrix1,2)==vfield%dim)
    assert(blocks(matrix1,1)==1)

    ndim = blocks(matrix1,2)

    if(.not.matrix1%sparsity%sorted_rows.or..not.matrix2%sparsity%sorted_rows) then
      FLAbort("mult_div_invvector_div_T assumes sorted rows")
    end if

    do i=1, size(matrix1%sparsity,1)
       row=>row_m_ptr(matrix1, i)

       if(size(row)>0) then
          do dim = 1, ndim
            row_val(dim)%ptr=>row_val_ptr(matrix1, 1, dim, i)
          end do

          row_product=>row_m_ptr(product, i)
          do jcol = 1, size(row_product)
             j = row_product(jcol)

             col=>row_m_ptr(matrix2, j)
             do dim = 1, ndim
               col_val(dim)%ptr=>row_val_ptr(matrix2, 1, dim, j)
             end do

             if(size(col)>0) then
                entry0=0.0

                addflag = .false.

                k1 = 1
                k2 = 1
                do
                   if((k1.gt.size(row)).or.(k2.gt.size(col))) exit
                   if(row(k1)<col(k2)) then
                      k1 = k1 + 1
                   else
                      if(row(k1)==col(k2)) then
                         ! Note the transpose in the second val call.
                         do dim = 1, ndim
                           entry0=entry0+row_val(dim)%ptr(k1)* &
                                col_val(dim)%ptr(k2)/node_val(vfield, dim, row(k1))
                         end do
                         addflag = .true.
                         k1 = k1 + 1
                         k2 = k2 + 1
                      else
                         k2 = k2 + 1
                      end if
                   end if
                end do
                if(addflag) then
                   nentry0 = nentry0 + 1
                   product%val(nentry0) = entry0
                end if
             end if
          end do
       end if
    end do

  end subroutine mult_div_invvector_div_T

  subroutine mult_div_tensorinvscalar_div_T(product, matrix1, tfield, sfield, matrix2, isotropic)
    !!< Perform the matrix multiplication:
    !!< 
    !!<     product = matrix1*tensor*diag(1./vector)*matrix2^T
    !!< 
    !!< Only works on block csr matrices with monotonic row entries in colm
    type(csr_matrix), intent(inout) :: product
    type(block_csr_matrix), intent(in) :: matrix1, matrix2
    type(tensor_field), intent(in) :: tfield
    type(scalar_field), intent(in) :: sfield
    logical, optional :: isotropic

    integer, dimension(:), pointer :: row, col, row_product
    type(real_vector), dimension(blocks(matrix1,2)) :: row_val, col_val
    integer :: i,j,k1,k2,jcol,dim1,dim2,ndim
    real :: entry0
    integer :: nentry0
    logical :: addflag
    logical :: l_isotropic

    ewrite(1,*) 'Entering mult_div_tensorinvscalar_div_T'

    if(present(isotropic)) then
      l_isotropic = isotropic
    else
      l_isotropic = .false.
    end if

    call zero(product)

    nentry0 = 0
    assert(size(matrix1,2)==size(matrix2,2))
    assert(blocks(matrix1,1)==blocks(matrix2,1))
    assert(blocks(matrix1,2)==blocks(matrix2,2))
    assert(blocks(matrix1,2)==tfield%dim(1))
    assert(blocks(matrix2,2)==tfield%dim(2))
    assert(blocks(matrix1,1)==1)

    ndim = blocks(matrix1,2)

    if(.not.matrix1%sparsity%sorted_rows.or..not.matrix2%sparsity%sorted_rows) then
      FLAbort("mult_div_tensorinvscalar_div_T assumes sorted rows")
    end if

    do i=1, size(matrix1%sparsity,1)
       row=>row_m_ptr(matrix1, i)

       if(size(row)>0) then
          do dim1 = 1, ndim
            row_val(dim1)%ptr=>row_val_ptr(matrix1, 1, dim1, i)
          end do

          row_product=>row_m_ptr(product, i)
          do jcol = 1, size(row_product)
             j = row_product(jcol)

             col=>row_m_ptr(matrix2, j)
             do dim1 = 1, ndim
               col_val(dim1)%ptr=>row_val_ptr(matrix2, 1, dim1, j)
             end do

             if(size(col)>0) then
                entry0=0.0

                addflag = .false.

                k1 = 1
                k2 = 1
                do
                   if((k1.gt.size(row)).or.(k2.gt.size(col))) exit
                   if(row(k1)<col(k2)) then
                      k1 = k1 + 1
                   else
                      if(row(k1)==col(k2)) then
                         ! Note the transpose in the second val call.
                         if(l_isotropic) then
                            do dim1 = 1, ndim
                                  entry0=entry0+row_val(dim1)%ptr(k1)* &
                                        col_val(dim1)%ptr(k2)*node_val(tfield, dim1, dim1, row(k1))&
                                        /node_val(sfield, row(k1))
                            end do
                         else
                            do dim1 = 1, ndim
                              do dim2 = 1, ndim
                                  entry0=entry0+row_val(dim1)%ptr(k1)* &
                                        col_val(dim2)%ptr(k2)*node_val(tfield, dim1, dim2, row(k1))&
                                        /node_val(sfield, row(k1))
                              end do
                            end do
                         end if
                         addflag = .true.
                         k1 = k1 + 1
                         k2 = k2 + 1
                      else
                         k2 = k2 + 1
                      end if
                   end if
                end do
                if(addflag) then
                   nentry0 = nentry0 + 1
                   product%val(nentry0) = entry0
                end if
             end if
          end do
       end if
    end do

  end subroutine mult_div_tensorinvscalar_div_T

  subroutine mult_div_tensorinvscalar_vector(product, matrix1, tfield, sfield, vfield, isotropic)
    !!< Perform the matrix multiplication:
    !!< 
    !!<     product = matrix1*tensor*diag(1./vector)*vfield
    !!< 
    !!< Note that product is not zeroed by this subroutine!
    type(scalar_field), intent(inout) :: product
    type(block_csr_matrix), intent(in) :: matrix1
    type(tensor_field), intent(in) :: tfield
    type(scalar_field), intent(in) :: sfield
    type(vector_field), intent(in) :: vfield
    logical, optional :: isotropic

    integer :: dim1, dim2, ndim, row
    type(real_vector), dimension(blocks(matrix1, 2)) :: row_val
    integer, dimension(:), pointer :: row_indices
    real :: entry0

    logical :: l_isotropic

    ewrite(1,*) 'Entering mult_div_tensorinvscalar_vector'

    if(present(isotropic)) then
      l_isotropic = isotropic
    else
      l_isotropic = .false.
    end if

    ndim = blocks(matrix1, 2)

    do row = 1, size(matrix1, 1)
      row_indices=>row_m_ptr(matrix1, row)

      if(size(row_indices)>0) then
        do dim1 = 1, ndim
          row_val(dim1)%ptr=>row_val_ptr(matrix1, 1, dim1, row)
        end do

        entry0 = 0.0

        if(l_isotropic) then
          do dim1 = 1, tfield%dim(1)
              entry0 = entry0 + &
                        sum((row_val(dim1)%ptr*node_val(tfield, dim1, dim1, row_indices)&
                            /node_val(sfield, row_indices))&
                        *node_val(vfield, dim1, row_indices))
          end do
        else
          do dim1 = 1, tfield%dim(1)
            do dim2 = 1, tfield%dim(2)
              entry0 = entry0 + &
                        sum((row_val(dim1)%ptr*node_val(tfield, dim1, dim2, row_indices)&
                            /node_val(sfield, row_indices))&
                        *node_val(vfield, dim2, row_indices))
            end do
          end do
        end if

        call addto(product, row, entry0)
      end if

    end do

  end subroutine mult_div_tensorinvscalar_vector

end module sparse_matrices_fields
