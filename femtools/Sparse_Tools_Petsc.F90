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
!    C.Pain@Imperial.ac.uk
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
module sparse_tools_petsc
  !!< This module is an extension to the sparse_tools module that 
  !!< implements a csr matrix type 'petsc_csr_matrix' that directly
  !!< stores the matrix in petsc format.
#include "petscversion.h"
  use FLDebug
  use Sparse_Tools
  use Reference_Counting
  use parallel_tools
  use halo_data_types
  use halos_allocates
  use fields_base
  use fields_manipulation
  use petsc_tools
#ifdef HAVE_PETSC_MODULES
  use petsc
#if PETSC_VERSION_MINOR==0
  use petscvec
  use petscmat
  use petscksp
  use petscpc
  use petscis
  use petscmg
#endif
#endif
  implicit none
#ifdef HAVE_PETSC_MODULES
#if PETSC_VERSION_MINOR==0
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscpcdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscisdef.h"
#else
#include "finclude/petscdef.h"
#endif
#else
#include "finclude/petsc.h"
#if PETSC_VERSION_MINOR==0
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscviewer.h"
#include "finclude/petscis.h"
#endif
#endif

  private
  
  type petsc_csr_matrix
     !! the matrix in PETSc format
     Mat :: M
     !! petsc numbering for rows and columns
     type(petsc_numbering_type) :: row_numbering, column_numbering
     
     !! The halos associated with the rows and columns of the matrix.
     type(halo_type), pointer :: row_halo => null(), column_halo => null()
     !! Reference counting
     type(refcount_type), pointer :: refcount => null()
     !! Name
     character(len=FIELD_NAME_LEN) :: name=""
     !! if .false. we need to call assemble before extracting info
     !! and or go into a solve:
     logical:: is_assembled=.false.
  end type petsc_csr_matrix

  type petsc_csr_matrix_pointer
    type(petsc_csr_matrix), pointer :: ptr => null()
  end type petsc_csr_matrix_pointer
  
    
  interface allocate
     module procedure allocate_petsc_csr_matrix_from_sparsity, &
       allocate_petsc_csr_matrix_from_nnz, &
       allocate_petsc_csr_matrix_from_petsc_matrix
  end interface
  
  interface deallocate
     module procedure deallocate_petsc_csr_matrix
  end interface

  interface size
     module procedure petsc_csr_size
  end interface

  interface block_size
     module procedure petsc_csr_block_size
  end interface

  interface blocks
     module procedure petsc_csr_blocks_withdim !, petsc_csr_blocks_nodim
  end interface
  
  interface entries
     module procedure petsc_csr_entries
  end interface

  interface zero
     module procedure petsc_csr_zero
  end interface

  interface addto
     module procedure petsc_csr_addto, petsc_csr_vaddto, &
       petsc_csr_blocks_addto_withmask, petsc_csr_block_addto, &
       petsc_csr_blocks_addto
  end interface

  interface addto_diag
     module procedure petsc_csr_addto_diag, petsc_csr_vaddto_diag
  end interface

  interface scale
     module procedure petsc_csr_scale
  end interface
  
  interface extract_diagonal
     ! this one would be more logical in Sparse_Matrices_Fields
     ! but it depends on petsc headers
     module procedure petsc_csr_extract_diagonal
  end interface

  interface mult_T
     module procedure petsc_csr_mult_T_vector, petsc_csr_mult_T_scalar_to_vector
  end interface
    
  interface mult
     module procedure petsc_csr_mult_vector, petsc_csr_mult_vector_to_scalar
  end interface
    
  interface assemble
     module procedure petsc_csr_assemble
  end interface
  
#include "Reference_count_interface_petsc_csr_matrix.F90"

  public :: petsc_csr_matrix, petsc_csr_matrix_pointer, &
     allocate, deallocate, &
     size, block_size, blocks, entries, &
     zero, addto, addto_diag, scale, &
     extract_diagonal, assemble, incref_petsc_csr_matrix, &
     ptap, mult, mult_T, dump_matrix

contains

  subroutine allocate_petsc_csr_matrix_from_sparsity(matrix, sparsity, blocks, &
      name, diagonal)
    !!< Allocates a petsc_csr_matrix, i.e. a csr_matrix variant
    !!< that directly stores in petsc format. The provided sparsity
    !!< is only used to workout the number of nonzeros per row and may be
    !!< thrown away after this call.
    type(petsc_csr_matrix), intent(out) :: matrix
    type(csr_sparsity), optional, intent(in):: sparsity
    integer, dimension(2), intent(in):: blocks
    character(len=*) :: name
    !! only take diagonal blocks into account when estimating matrix sparsity
    !! does not change the matrix structure otherwise
    logical, optional, intent(in):: diagonal

    PetscErrorCode:: ierr
    logical:: ldiagonal
    integer:: nprows
    
    ldiagonal=present_and_true(diagonal)

    matrix%name = name
    
    call allocate( matrix%row_numbering, &
      nnodes=size(sparsity,1), &
      nfields=blocks(1), &
      halo=sparsity%row_halo )
      
    if (size(sparsity,1)==size(sparsity,2) .and. blocks(1)==blocks(2) .and. &
      associated(sparsity%row_halo, sparsity%column_halo)) then
        
      ! row and column numbering are the same
      matrix%column_numbering=matrix%row_numbering
      call incref(matrix%column_numbering)
      
    else
    
      ! create seperate column numbering
      call allocate( matrix%column_numbering, &
        nnodes=size(sparsity,2), &
        nfields=blocks(2), &
        halo=sparsity%column_halo )
    end if
    
    if (.not. IsParallel()) then

      ! Create serial matrix:
      matrix%M=csr2petsc_CreateSeqAIJ(sparsity, matrix%row_numbering, &
        matrix%column_numbering, ldiagonal)
      
    else
    
      if (associated(sparsity%row_halo)) then
        if (sparsity%row_halo%data_type==HALO_TYPE_CG_NODE) then
          ! Mask out non-local rows, as the assembled bits in those
          ! will be incomplete and need to be thrown out. In the case
          ! of DG assembly however the local bits are proper contributions
          ! and need to be added in the global matrix
          nprows=matrix%row_numbering%nprivatenodes
          matrix%row_numbering%gnn2unn(nprows+1:,:)=-1
        end if
      end if
        
      ! Create parallel matrix:
      matrix%M=csr2petsc_CreateMPIAIJ(sparsity, matrix%row_numbering, &
        matrix%column_numbering, ldiagonal)
      
    endif
    
    ! this is very important for assembly routines (e.g. DG IP viscosity)
    ! that try to add zeros outside the provided sparsity; if we go outside
    ! the provided n/o nonzeros the assembly will become very slow!!!
    call MatSetOption(matrix%M, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr)

    ! to make sure we're not underestimating the number of nonzeros ever, make
    ! petsc fail if new allocations are necessary. If uncommenting the setting of this
    ! option fixes your problem the number of no
    call MatSetOption(matrix%M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE, ierr)

    ! saves us from doing a transpose for block inserts
    call MatSetOption(matrix%M, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

    if (associated(sparsity%row_halo)) then
      ! these are also pointed to in the row_numbering
      ! but only refcounted here
      matrix%row_halo => sparsity%row_halo
      call incref(matrix%row_halo)
    end if
    
    if (associated(sparsity%column_halo)) then
      ! these are also pointed to in the column_numbering
      ! but only refcounted here
      matrix%column_halo => sparsity%column_halo
      call incref(matrix%column_halo)
    end if
    
    nullify(matrix%refcount) ! Hack for gfortran component initialisation
    !                         bug.
    call addref(matrix)

  end subroutine allocate_petsc_csr_matrix_from_sparsity

  subroutine allocate_petsc_csr_matrix_from_nnz(matrix, rows, columns, &
      dnnz, onnz, blocks, name, halo, row_halo, column_halo, &
      element_size)
    !!< Allocates a petsc_csr_matrix, i.e. a csr_matrix variant
    !!< that directly stores in petsc format. For this version the number
    !!< of nonzeros in each row needs to be provided explicitly. This allows
    !!< for a more fine-grained nnz estimate in case of different sparsities
    !!< on off-diagonal blocks (i.e. for component to component coupling)
    !!< dnnz is the number of local entries in the row (i.e. owned by this
    !!< process) onnz is the number of non-local/not-owned entries. This has
    !!< to be specified for all rows of the different vertical blocks, i.e.
    !!< size(dnnz)==size(onnz)==nprows*blocks(1), where nprows is the number
    !!< of private rows (contiguisly numbered). In serial only dnnz is used
    !!< and size(dnnz)==rows*blocks(1).
    type(petsc_csr_matrix), intent(out) :: matrix
    integer, intent(in):: rows, columns
    integer, dimension(:), intent(in):: dnnz, onnz
    integer, dimension(2), intent(in):: blocks
    character(len=*), intent(in) :: name
    type(halo_type), pointer, optional:: halo, row_halo, column_halo
    !!< If provided, the actual PETSc matrix will employ a block structure
    !!< where each local block consists of the degrees of freedom 
    !!< of 'element_size' subsequent indices. Note that these blocks are
    !!< something entirely different than the blocks of the previous 'blocks'
    !!< argument, which only refer to dim argument in the set/addto interface
    !!< i.e. call addto(M, dim1, dim2, i1, i2, val) where
    !!< 1<=dim1<=blocks(1) and 1<=dim2<=blocks(2).
    !!< To distinguish the element_size blocks will be referred to as
    !!< element blocks. They consist of entries that have the same
    !!< value for (i1-1)/element_size and (i2-1)/element_size (integer division)
    !!< At the moment entries which have a different dim1 or dim2 are not
    !!< included in the element block (this might be a future option).
    !!< If provided the arguments rows and columns change their meaning
    !!< to be the number of element block rows and columns. The size of dnnz and onnz
    !!< is still nprows*blocks(1), but they now contain the number of 
    !!< nonzero element blocks in an element block row.
    integer, intent(in), optional:: element_size

    type(halo_type), pointer:: lrow_halo, lcolumn_halo
    PetscErrorCode:: ierr
    integer:: nprows, npcols, urows, ucols
    integer:: index_rows, index_columns
    logical:: use_element_blocks
    
    matrix%name = name
    
    nullify( lrow_halo )
    nullify( lcolumn_halo )
    if (present(halo)) then
      lrow_halo => halo
      lcolumn_halo => halo
    end if
    if(present(row_halo)) then
      lrow_halo => row_halo
    end if
    if(present(column_halo)) then
      lcolumn_halo => column_halo
    end if
    
    if (present(element_size)) then
      index_rows=rows*element_size
      index_columns=columns*element_size
      use_element_blocks= element_size>1
    else
      index_rows=rows
      index_columns=columns
      use_element_blocks=.false.
    end if
    
    call allocate( matrix%row_numbering, &
      nnodes=index_rows, &
      nfields=blocks(1), &
      halo=lrow_halo )
      
    if (rows==columns .and. blocks(1)==blocks(2) .and. &
      associated(lrow_halo, lcolumn_halo)) then
        
      ! row and column numbering are the same
      matrix%column_numbering=matrix%row_numbering
      call incref(matrix%column_numbering)
      
    else
    
      ! create seperate column numbering
      call allocate( matrix%column_numbering, &
        nnodes=index_columns, &
        nfields=blocks(2), &
        halo=lcolumn_halo )
    end if
      
    urows=matrix%row_numbering%universal_length
    ucols=matrix%column_numbering%universal_length
    
    if (IsParallel()) then
      nprows=matrix%row_numbering%nprivatenodes
      npcols=matrix%column_numbering%nprivatenodes
      if (associated(lrow_halo)) then
        if (lrow_halo%data_type==HALO_TYPE_CG_NODE) then
          ! Mask out non-local rows, as the assembled bits in those
          ! will be incomplete and need to be thrown out. In the case
          ! of DG assembly however the local bits are proper contributions
          ! and need to be added in the global matrix
          matrix%row_numbering%gnn2unn(nprows+1:,:)=-1
        end if
      end if
    end if
    
    if (use_element_blocks .and. .not. IsParallel()) then
      
      assert( size(dnnz)==urows/element_size )
      
      ! Create serial block matrix:
      call MatCreateSeqBAIJ(MPI_COMM_SELF, element_size, &
         urows, ucols, PETSC_NULL_INTEGER, &
         dnnz, matrix%M, ierr)
         
    elseif (use_element_blocks) then
      
      assert( size(dnnz)==nprows*blocks(1)/element_size )
      assert( size(onnz)==nprows*blocks(1)/element_size )
      
      call MatCreateMPIBAIJ(MPI_COMM_WORLD, element_size, &
         nprows*blocks(1), npcols*blocks(2), &
         urows, ucols, &
         PETSC_NULL_INTEGER, dnnz, PETSC_NULL_INTEGER, onnz, matrix%M, ierr)
    
    else if (.not. IsParallel()) then

      assert( size(dnnz)==urows )
      
      ! Create serial matrix:
      call MatCreateSeqAIJ(MPI_COMM_SELF, urows, ucols, PETSC_NULL_INTEGER, &
         dnnz, matrix%M, ierr)
      
    else
      
      assert( size(dnnz)==nprows*blocks(1) )
      assert( size(onnz)==nprows*blocks(1) )
      
      call MatCreateMPIAIJ(MPI_COMM_WORLD, nprows*blocks(1), npcols*blocks(2), &
         urows, ucols, &
         PETSC_NULL_INTEGER, dnnz, PETSC_NULL_INTEGER, onnz, matrix%M, ierr)
      
    endif
    
    if (.not. use_element_blocks) then
      ! this is very important for assembly routines (e.g. DG IP viscosity)
      ! that try to add zeros outside the provided sparsity; if we go outside
      ! the provided n/o nonzeros the assembly will become very slow!!!
      call MatSetOption(matrix%M, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr)
    end if

    ! to make sure we're not underestimating the number of nonzeros ever, make
    ! petsc fail if new allocations are necessary. If uncommenting the setting of this
    ! option fixes your problem the number of no
    call MatSetOption(matrix%M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE, ierr)

    ! saves us from doing a transpose for block inserts
    call MatSetOption(matrix%M, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

    if (associated(lrow_halo)) then
      ! these are also pointed to in the row_numbering
      ! but only refcounted here
      matrix%row_halo => lrow_halo
      call incref(matrix%row_halo)
    end if
    
    if (associated(lcolumn_halo)) then
      ! these are also pointed to in the column_numbering
      ! but only refcounted here
      matrix%column_halo => lcolumn_halo
      call incref(matrix%column_halo)
    end if
    
    nullify(matrix%refcount) ! Hack for gfortran component initialisation
    !                         bug.
    call addref(matrix)

  end subroutine allocate_petsc_csr_matrix_from_nnz
    
  subroutine allocate_petsc_csr_matrix_from_petsc_matrix(matrix, &
    M, row_numbering, column_numbering, name)
    !!< Allocates a petsc_csr_matrix using an already created real petsc matrix
    !!< row_numbering and column_numbering have to be supplied to specify the
    !!< relation between the numbering used inside the petsc matrix and the
    !!< numbering to be used for the interface. References to those numberings
    !!< will be added. The supplied petsc matrix must be in assembled state.
    !!< After this it should be possible to add new entries and zero the matrix 
    !!< through the petsc_csr_matrix interface, but only if all nonzero
    !!< entries have been preallocated.
    type(petsc_csr_matrix), intent(out):: matrix
    Mat, intent(in):: M
    type(petsc_numbering_type), intent(in):: row_numbering, column_numbering
    character(len=*), intent(in):: name
    
    PetscErrorCode:: ierr
      
    matrix%M=M
    matrix%name=name
    matrix%is_assembled=.true.
    
    matrix%row_numbering=row_numbering
    call incref(matrix%row_numbering)
    
    matrix%column_numbering=column_numbering
    call incref(matrix%column_numbering)
    
    if (associated(row_numbering%halo)) then
      matrix%row_halo => row_numbering%halo
      call incref(row_numbering%halo)
    end if
    
    if (associated(column_numbering%halo)) then
      matrix%column_halo => column_numbering%halo
      call incref(column_numbering%halo)
    end if
    
    ! make sure the matrix options are consistent with petsc_csr_matrix interface
    
    ! this is very important for assembly routines (e.g. DG IP viscosity)
    ! that try to add zeros outside the provided sparsity; if we go outside
    ! the provided n/o nonzeros the assembly will become very slow!!!
    call MatSetOption(matrix%M, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr)

    ! to make sure we're not underestimating the number of nonzeros ever, make
    ! petsc fail if new allocations are necessary. If uncommenting the setting of this
    ! option fixes your problem the number of no
    call MatSetOption(matrix%M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE, ierr)

    ! saves us from doing a transpose for block inserts
    call MatSetOption(matrix%M, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
    
    nullify(matrix%refcount) ! Hack for gfortran component initialisation
    !                         bug.
    call addref(matrix)
    
  end subroutine allocate_petsc_csr_matrix_from_petsc_matrix
  
  subroutine deallocate_petsc_csr_matrix(matrix, stat)
    type(petsc_csr_matrix), intent(inout) :: matrix
    integer, intent(out), optional :: stat

    integer :: lstat

    lstat=0

    call decref(matrix)
    if (has_references(matrix)) then
       goto 42
    end if
    
    call MatDestroy(matrix%M, lstat)
    if (lstat/=0) goto 42
    
    call deallocate(matrix%row_numbering)
    
    call deallocate(matrix%column_numbering)
    
    if (associated(matrix%row_halo)) then
       call deallocate(matrix%row_halo)
    end if
    
    if (associated(matrix%column_halo)) then
       call deallocate(matrix%column_halo)
    end if
    
42  if (present(stat)) then
       stat=lstat
    else
       if (lstat/=0) then
          FLAbort("Failed to deallocate matrix")
       end if
    end if

  end subroutine deallocate_petsc_csr_matrix

  pure function petsc_csr_size(matrix, dim)
    !!< Clone of size function.
    integer :: petsc_csr_size
    type(petsc_csr_matrix), intent(in) :: matrix
    integer, optional, intent(in) :: dim

    integer :: rows, cols
    
    rows = size(matrix%row_numbering%gnn2unn)
    cols = size(matrix%column_numbering%gnn2unn)
    
    if (.not.present(dim)) then
       petsc_csr_size = rows * cols
    else if (dim==1) then
       petsc_csr_size = rows
    else if (dim==2) then
       petsc_csr_size = cols
    else
       ! not allowed to flabort in pure function
       petsc_csr_size = 0
    end if
    
  end function petsc_csr_size

  pure function petsc_csr_block_size(matrix, dim)
    !!< size of each block
    integer :: petsc_csr_block_size
    type(petsc_csr_matrix), intent(in) :: matrix
    integer, optional, intent(in) :: dim

    integer :: rows, cols
    
    rows = size(matrix%row_numbering%gnn2unn,1)
    cols = size(matrix%column_numbering%gnn2unn,1)
    
    if (.not.present(dim)) then
       petsc_csr_block_size = rows * cols
    else if (dim==1) then
       petsc_csr_block_size = rows
    else if (dim==2) then
       petsc_csr_block_size = cols
    else
       ! not allowed to flabort in pure function
       petsc_csr_block_size = 0
    end if
    
  end function petsc_csr_block_size
  
  pure function petsc_csr_blocks_withdim(matrix, dim) result (blocks)
    !!< Number of blocks
    integer :: blocks
    type(petsc_csr_matrix), intent(in) :: matrix
    integer, optional, intent(in) :: dim

    integer :: rows, cols
    
    rows = size(matrix%row_numbering%gnn2unn,2)
    cols = size(matrix%column_numbering%gnn2unn,2)
    
    if (.not.present(dim)) then
       blocks = rows * cols
    else if (dim==1) then
       blocks = rows
    else if (dim==2) then
       blocks = cols
    else
       ! not allowed to flabort in pure function
       blocks = 0
    end if
    
  end function petsc_csr_blocks_withdim
  
! causes gfortran to complain about ambiguous generic interface:
! 
!    pure function petsc_csr_blocks_nodim(matrix) result (blocks)
!      !!< Number of blocks
!      integer, dimension(2) :: blocks
!      type(petsc_csr_matrix), intent(in) :: matrix
!  
!      integer :: rows, cols
!      
!      rows = size(matrix%row_numbering%gnn2unn,2)
!      cols = size(matrix%column_numbering%gnn2unn,2)
!      
!      blocks = (/ rows, cols /)
!      
!    end function petsc_csr_blocks_nodim
! ============================================================================
  
  function petsc_csr_entries(matrix) result (entries)
    !!< Return the number of (potentially) non-zero entries in matrix.
    integer :: entries
    type(petsc_csr_matrix), intent(in) :: matrix
      
    double precision, dimension(MAT_INFO_SIZE):: matrixinfo
    PetscErrorCode:: ierr

    ! get the necessary info about the matrix:
    call myMatGetInfo(matrix%M, MAT_LOCAL, matrixinfo, ierr)
    entries=matrixinfo(MAT_INFO_NZ_USED)

  end function petsc_csr_entries
  
  subroutine petsc_csr_zero(matrix)
    !!< Zero the entries of a csr matrix.
    type(petsc_csr_matrix), intent(inout) :: matrix

    PetscErrorCode:: ierr
    
    call MatZeroEntries(matrix%M, ierr)
    matrix%is_assembled=.true.
    
  end subroutine petsc_csr_zero
  
  subroutine petsc_csr_addto(matrix, blocki, blockj, i, j, val)
    !!< Add value to matrix(blocki, blockj, i,j)
    type(petsc_csr_matrix), intent(inout) :: matrix
    integer, intent(in) :: blocki,blockj,i,j
    real, intent(in) :: val

    PetscErrorCode:: ierr
    integer:: row, col
    
    row=matrix%row_numbering%gnn2unn(i,blocki)
    col=matrix%column_numbering%gnn2unn(j,blockj)
    
    call MatSetValue(matrix%M, row, col, val, ADD_VALUES, ierr)
    matrix%is_assembled=.false.

  end subroutine petsc_csr_addto

  subroutine petsc_csr_vaddto(matrix, blocki, blockj, i, j, val)
    !!< Add multiple values to matrix(blocki, blockj, i,j)
    type(petsc_csr_matrix), intent(inout) :: matrix
    integer, intent(in) :: blocki,blockj
    integer, dimension(:), intent(in) :: i,j
    real, dimension(size(i),size(j)), intent(in) :: val
    
    PetscInt, dimension(size(i)):: idxm
    PetscInt, dimension(size(j)):: idxn
    PetscErrorCode:: ierr
    
    integer:: k
    
    idxm=matrix%row_numbering%gnn2unn(i,blocki)
    idxn=matrix%column_numbering%gnn2unn(j,blockj)
    
    if (IsParallel() .and. maxval(i)>matrix%row_numbering%nprivatenodes) then
      ! due to a bug in petsc, the option MAT_IGNORE_ZERO_ENTRIES causes
      ! an entire non-owned row to be dropped if only its first entry is zero
      ! therefore we have to assemble column by column
      do k=1, size(j)
        ! luckily columns are contiguous in memory in fortran, so no copy, required
        call MatSetValues(matrix%M, size(i), idxm, 1, idxn(k:k), val(:,k), &
          ADD_VALUES, ierr)
      end do
    else
    
      call MatSetValues(matrix%M, size(i), idxm, size(j), idxn, val, &
        ADD_VALUES, ierr)
    end if
    matrix%is_assembled=.false.

  end subroutine petsc_csr_vaddto
  
  subroutine petsc_csr_block_addto(matrix, i, j, val)
    !!< Adds a local matrix for all components of entry (i,j) in the matrix
    
    type(petsc_csr_matrix), intent(inout) :: matrix
    integer, intent(in) :: i
    integer, intent(in) :: j
    real, dimension(:,:), intent(in) :: val
    
    PetscInt, dimension(size(matrix%row_numbering%gnn2unn,2)):: idxm
    PetscInt, dimension(size(matrix%column_numbering%gnn2unn,2)):: idxn
    PetscErrorCode:: ierr
    integer:: blockj
    logical:: insert_per_column
    
    ! due to a bug in petsc, the option MAT_IGNORE_ZERO_ENTRIES causes
    ! an entire non-owned row to be dropped if only its first entry is zero
    ! therefore we have to assemble column by column, if:
    insert_per_column=IsParallel() .and. &
        i>matrix%row_numbering%nprivatenodes
        
    idxm=matrix%row_numbering%gnn2unn(i,:)
    idxn=matrix%column_numbering%gnn2unn(j,:)
    
    if (insert_per_column) then
      do blockj=1, size(matrix%column_numbering%gnn2unn,2)
        call MatSetValues(matrix%M, size(idxm), idxm, 1, idxn(blockj:blockj), &
                  val(:,blockj), ADD_VALUES, ierr)
      end do
    else
      call MatSetValues(matrix%M, size(idxm), idxm, size(idxn), idxn, &
                  val, ADD_VALUES, ierr)
    end if
    matrix%is_assembled=.false.

  end subroutine petsc_csr_block_addto
  
  subroutine petsc_csr_blocks_addto(matrix, i, j, val)
    !!< Add the (blocki, blockj, :, :) th matrix of val onto the (blocki, blockj) th
    !!< block of the block csr matrix, for all blocks of the block csr matrix.
    
    type(petsc_csr_matrix), intent(inout) :: matrix
    integer, dimension(:), intent(in) :: i
    integer, dimension(:), intent(in) :: j
    real, dimension(:,:,:,:), intent(in) :: val
    
    PetscScalar, dimension(size(i), size(j)):: value
    PetscInt, dimension(size(i)):: idxm
    PetscInt, dimension(size(j)):: idxn
    PetscErrorCode:: ierr
    integer:: blocki, blockj, k
    logical:: insert_per_column
    
    ! due to a bug in petsc, the option MAT_IGNORE_ZERO_ENTRIES causes
    ! an entire non-owned row to be dropped if only its first entry is zero
    ! therefore we have to assemble column by column, if:
    insert_per_column=IsParallel() .and. &
        maxval(i)>matrix%row_numbering%nprivatenodes
    
    do blocki=1, size(matrix%row_numbering%gnn2unn,2)
      idxm=matrix%row_numbering%gnn2unn(i,blocki)
      do blockj=1, size(matrix%column_numbering%gnn2unn,2)
        idxn=matrix%column_numbering%gnn2unn(j,blockj)
        ! unfortunately we need a copy here to pass contiguous memory
        value=val(blocki, blockj, :, :)
        if (insert_per_column) then
          do k=1, size(j)
            ! luckily columns are contiguous in memory in fortran, so no copy, required
            call MatSetValues(matrix%M, size(i), idxm, 1, idxn(k:k), &
                value(:,k), ADD_VALUES, ierr)
          end do
        else
          call MatSetValues(matrix%M, size(i), idxm, size(j), idxn, &
              value, ADD_VALUES, ierr)
        end if
      end do
    end do
    matrix%is_assembled=.false.

  end subroutine petsc_csr_blocks_addto
  
  subroutine petsc_csr_blocks_addto_withmask(matrix, i, j, val, block_mask)
    !!< Add the (blocki, blockj, :, :) th matrix of val onto the (blocki, blockj) th
    !!< block of the block csr matrix, for all blocks of the block csr matrix.
    
    type(petsc_csr_matrix), intent(inout) :: matrix
    integer, dimension(:), intent(in) :: i
    integer, dimension(:), intent(in) :: j
    real, dimension(:,:,:,:), intent(in) :: val
    logical, dimension(:,:), intent(in) :: block_mask
    
    PetscScalar, dimension(size(i), size(j)):: value
    PetscInt, dimension(size(i)):: idxm
    PetscInt, dimension(size(j)):: idxn
    PetscErrorCode:: ierr
    integer:: blocki, blockj, k
    logical:: insert_per_column
    
    ! due to a bug in petsc, the option MAT_IGNORE_ZERO_ENTRIES causes
    ! an entire non-owned row to be dropped if only its first entry is zero
    ! therefore we have to assemble column by column, if:
    insert_per_column=IsParallel() .and. &
        maxval(i)>matrix%row_numbering%nprivatenodes
    
    do blocki=1, size(matrix%row_numbering%gnn2unn,2)
      idxm=matrix%row_numbering%gnn2unn(i,blocki)
      do blockj=1, size(matrix%column_numbering%gnn2unn,2)
        if (block_mask(blocki,blockj)) then
          idxn=matrix%column_numbering%gnn2unn(j,blockj)
          ! unfortunately we need a copy here to pass contiguous memory
          value=val(blocki, blockj, :, :)
          if (insert_per_column) then
            do k=1, size(j)
              ! luckily columns are contiguous in memory in fortran, so no copy, required
              call MatSetValues(matrix%M, size(i), idxm, 1, idxn(k:k), &
                  value(:,k), ADD_VALUES, ierr)
            end do
          else
            call MatSetValues(matrix%M, size(i), idxm, size(j), idxn, &
                value, ADD_VALUES, ierr)
          end if
        end if
      end do
    end do
    matrix%is_assembled=.false.

  end subroutine petsc_csr_blocks_addto_withmask
  
  subroutine petsc_csr_scale(matrix, scale)
    !!< Scale matrix by scale.
    type(petsc_csr_matrix), intent(inout) :: matrix
    real, intent(in) :: scale
    
    PetscErrorCode:: ierr
    
    call MatScale(matrix%M, scale, ierr)
    matrix%is_assembled=.false. ! I think?
    
  end subroutine petsc_csr_scale

  subroutine petsc_csr_addto_diag(matrix, blocki, blockj, i, val)
    !!< Add val to matrix(i,i)
    !!< Adding to the diagonal of a non-diagonal block is supported.
    type(petsc_csr_matrix), intent(inout) :: matrix
    integer, intent(in) :: blocki,blockj, i
    real, intent(in) :: val

    call addto(matrix, blocki, blockj, i, i, val)
    matrix%is_assembled=.false.

  end subroutine petsc_csr_addto_diag

  subroutine petsc_csr_vaddto_diag(matrix, blocki, blockj, i, val)
    !!< Add val to matrix(i,i)
    type(petsc_csr_matrix), intent(inout) :: matrix
    integer, intent(in) :: blocki, blockj
    integer, dimension(:), intent(in) :: i
    real, dimension(size(i)), intent(in) :: val
    
    integer:: k
    
    ! can't think of a more efficient way
    do k=1, size(i)
      call addto(matrix, blocki, blockj, i(k), i(k), val(k))
    end do
    matrix%is_assembled=.false.

  end subroutine petsc_csr_vaddto_diag

  subroutine petsc_csr_extract_diagonal(matrix,diagonal)
    !!< Extracts diagonal components of a block_csr matrix.
    !!< The vector field diagonal needs to be allocated before the call.

    type(petsc_csr_matrix), intent(inout) :: matrix
    type(vector_field), intent(inout) :: diagonal
      
    PetscErrorCode:: ierr
    Vec:: diagonal_vec
    
    assert( diagonal%dim==blocks(matrix,1) )
    assert( node_count(diagonal)==block_size(matrix,1))
    assert( block_size(matrix,1)==block_size(matrix,2))
    assert( blocks(matrix,1)==blocks(matrix,2))
    
    call petsc_csr_assemble(matrix)
    
    diagonal_vec=PetscNumberingCreateVec(matrix%row_numbering)
    call MatGetDiagonal(matrix%M, diagonal_vec, ierr)
    call petsc2field(diagonal_vec, matrix%row_numbering, diagonal)

  end subroutine petsc_csr_extract_diagonal
    
  subroutine petsc_csr_assemble(matrix)
    !!< if necessary assemble the matrix
    type(petsc_csr_matrix), intent(inout) :: matrix
    
    PetscErrorCode:: ierr
    
    call alland(matrix%is_assembled)
    if (.not. matrix%is_assembled) then
      call MatAssemblyBegin(matrix%M, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(matrix%M, MAT_FINAL_ASSEMBLY, ierr)
    end if
    matrix%is_assembled=.true.
    
  end subroutine petsc_csr_assemble
    
  subroutine ptap(c, a, p)
    !!< Perform the matrix multiplication A=P^T A P
    type(petsc_csr_matrix), intent(out):: c
    type(petsc_csr_matrix), intent(inout):: a
    type(petsc_csr_matrix), intent(inout):: p
    
    PetscErrorCode:: ierr
    
    call assemble(a)
    call assemble(p)
    
    ! this creates the petsc ptap matrix and computes it
    call MatPTAP(a%M, p%M, MAT_INITIAL_MATRIX, 1.5, c%M, ierr)
    
    ! rest of internals for c is copied from A
    c%row_numbering=a%row_numbering
    call incref(c%row_numbering)
    c%column_numbering=a%column_numbering
    call incref(c%column_numbering)
    
    if (associated(a%row_halo)) then
      c%row_halo => a%row_halo
      call incref(c%row_halo)
    else
      nullify(c%row_halo)
    end if
    if (associated(a%column_halo)) then
      c%column_halo => a%column_halo
      call incref(c%column_halo)
    else
      nullify(c%column_halo)
    end if
    
    ! I think it is assembled now?
    c%is_assembled=.true.
    
    ! make up a name
    c%name=trim(a%name)//"_"//trim(p%name)//"_ptap"
    
    ! the new c get its own reference:
    nullify(c%refcount) ! Hack for gfortran component initialisation
    !                         bug.
    call addref(c)
    
  end subroutine ptap
    
  subroutine petsc_csr_mult_vector(x, A, b)
    !!< Performs the matrix-vector multiplication x=Ab
    type(vector_field), intent(inout):: x
    type(petsc_csr_matrix), intent(inout):: A
    type(vector_field), intent(in):: b
    
    PetscErrorCode:: ierr
    Vec:: bvec, xvec
 
    assert( node_count(x)==block_size(A, 1) )
    assert( node_count(b)==block_size(A, 2) )
    assert( x%dim==blocks(A,1) )
    assert( b%dim==blocks(A,2) )

    call petsc_csr_assemble(A)
    
    ! copy b to petsc vector
    bvec=PetscNumberingCreateVec(A%column_numbering)
    call field2petsc(b, A%column_numbering, bvec)
    ! creates PETSc solution vec of the right size:
    xvec=PetscNumberingCreateVec(A%row_numbering)

    ! perform the multiply
    call MatMult(A%M, bvec, xvec, ierr)
    ! copy answer back to vector_field
    call petsc2field(xvec, A%row_numbering, x)
    ! destroy the PETSc vecs
    call VecDestroy(bvec, ierr)
    call VecDestroy(xvec, ierr)
    
  end subroutine petsc_csr_mult_vector
  
  subroutine petsc_csr_mult_vector_to_scalar(x, A, b)
    !!< Performs the matrix-scalar multiplication x=Ab
    type(scalar_field), intent(inout):: x
    type(petsc_csr_matrix), intent(inout):: A
    type(vector_field), intent(in):: b
    
    PetscErrorCode:: ierr
    Vec:: bvec, xvec
    
    assert( node_count(x)==block_size(A, 1) )
    assert( node_count(b)==block_size(A, 2) )
    assert( 1==blocks(A,1) )
    assert( b%dim==blocks(A,2) )
    
    call petsc_csr_assemble(A)

    ! copy b to petsc vector
    bvec=PetscNumberingCreateVec(A%column_numbering)
    call field2petsc(b, A%column_numbering, bvec)
    ! creates PETSc solution vec of the right size:
    xvec=PetscNumberingCreateVec(A%row_numbering)
    ! perform the multiply
    call MatMult(A%M, bvec, xvec, ierr)
    ! copy answer back to vector_field
    call petsc2field(xvec, A%row_numbering, x)
    ! destroy the PETSc vecs
    call VecDestroy(bvec, ierr)
    call VecDestroy(xvec, ierr)
    
  end subroutine petsc_csr_mult_vector_to_scalar
  
  subroutine petsc_csr_mult_T_vector(x, A, b)
    !!< Performs the matrix-vector multiplication x=Ab
    type(vector_field), intent(inout):: x
    type(petsc_csr_matrix), intent(inout):: A
    type(vector_field), intent(in):: b
    
    PetscErrorCode:: ierr
    Vec:: bvec, xvec
    
    assert( node_count(x)==block_size(A, 2) )
    assert( node_count(b)==block_size(A, 1) )
    assert( x%dim==blocks(A,2) )
    assert( b%dim==blocks(A,1) )

    call petsc_csr_assemble(A)    

    ! copy b to petsc vector
    bvec=PetscNumberingCreateVec(A%row_numbering)
    call field2petsc(b, A%row_numbering, bvec)
    ! creates PETSc solution vec of the right size:
    xvec=PetscNumberingCreateVec(A%column_numbering)
    ! perform the multiply
    call MatMultTranspose(A%M, bvec, xvec, ierr)
    ! copy answer back to vector_field
    call petsc2field(xvec, A%column_numbering, x)
    ! destroy the PETSc vecs
    call VecDestroy(bvec, ierr)
    call VecDestroy(xvec, ierr)
    
  end subroutine petsc_csr_mult_T_vector
  
  subroutine petsc_csr_mult_T_scalar_to_vector(x, A, b)
    !!< Performs the matrix-scalar multiplication x=Ab
    type(vector_field), intent(inout):: x
    type(petsc_csr_matrix), intent(inout):: A
    type(scalar_field), intent(in):: b
    
    PetscErrorCode:: ierr
    Vec:: bvec, xvec
    
    assert( node_count(x)==block_size(A, 2) )
    assert( node_count(b)==block_size(A, 1) )
    assert( x%dim==blocks(A,2) )
    assert( 1==blocks(A,1) )
    
    call petsc_csr_assemble(A)

    ! copy b to petsc vector
    bvec=PetscNumberingCreateVec(A%row_numbering)
    call field2petsc(b, A%row_numbering, bvec)
    ! creates PETSc solution vec of the right size:
    xvec=PetscNumberingCreateVec(A%column_numbering)
    ! perform the multiply
    call MatMultTranspose(A%M, bvec, xvec, ierr)
    ! copy answer back to vector_field
    call petsc2field(xvec, A%column_numbering, x)
    ! destroy the PETSc vecs
    call VecDestroy(bvec, ierr)
    call VecDestroy(xvec, ierr)
    
  end subroutine petsc_csr_mult_T_scalar_to_vector
  
  subroutine dump_matrix(name,A)
    character(len=*), intent(in):: name
    type(petsc_csr_matrix):: A
    Vec:: x0, b
    
    x0=PetscNumberingCreateVec(A%column_numbering)
    b=PetscNumberingCreateVec(A%row_numbering)
    call DumpMatrixEquation(name, x0, A%M, b)
    
    end subroutine dump_matrix
  
#include "Reference_count_petsc_csr_matrix.F90"
end module sparse_tools_petsc
