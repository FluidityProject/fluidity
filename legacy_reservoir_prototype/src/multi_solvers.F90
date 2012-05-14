
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

module solvers_module

  use fldebug
  use sparse_tools_petsc
  use solvers
  use fields

  implicit none
    
  private
  
  public  :: solver
  
  interface solver
     module procedure solve_via_copy_to_petsc_csr_matrix
  end interface solver
  
contains

! -----------------------------------------------------------------------------

  subroutine solve_via_copy_to_petsc_csr_matrix(A, &
                                                x, &
                                                b, &
                                                findfe, &
                                                colfe, &
                                                option_path)
  
     !!< Solve a matrix Ax = b system via copying over to the
     !!< petsc_csr_matrix type and calling the femtools solver
     !!< using the spud options given by the field options path
     
    integer, dimension(:), intent(in) :: findfe 
    integer, dimension(:), intent(in) :: colfe
    real, dimension(:), intent(in) :: a,b
    real, dimension(:), intent(inout) :: x
    character(len=*), intent(in) :: option_path
    
    ! local variables
    integer :: i,j,k
    integer :: rows
    integer, dimension(:), allocatable :: dnnz
    type(petsc_csr_matrix) :: matrix
    type(scalar_field) :: rhs
    type(scalar_field) :: solution
    
    rows = size(x)
    
    print *, size(x)+1, size(findfe)

    assert(size(x) == size(b))
    assert(size(a) == size(colfe))
    assert((size(x)+1) == size(findfe))
    
    ! find the number of non zeros per row
    allocate(dnnz(rows))
    
    dnnz = 0
    
    do i = 1,rows
       
       dnnz(i) = findfe(i+1) - findfe(i)
       
    end do
    
    ! allocate the petsc_csr_matrix using nnz (pass dnnz also for onnz) 
    call allocate(matrix, &
                  rows, &
                  rows, &
                  dnnz, &
                  dnnz, &
                  (/1,1/), &
                  name = 'dummy')
    
    call zero(matrix)
    
    ! add in the entries to petsc matrix
    do i = 1, rows
    
      do j = findfe(i), findfe(i+1) - 1

          k = colfe(j)
          
          call addto(matrix, &
                     blocki = 1, &
                     blockj = 1, &
                     i = i, &
                     j = k, &
                     val = a(j))
       
       end do 
    
    end do 

    call assemble(matrix)    
    
    ! Set up rhs and initial guess which are scalar field types.        
    allocate(rhs%val(rows))
    allocate(solution%val(rows))
        
    do i = 1,rows
    
       call set(rhs, &
                i, &
                b(i))
                
       call set(solution, &
                i, &
                x(i))
    
    end do
        
    ! solve matrix
    call petsc_solve(solution, &
                     matrix, &
                     rhs, &
                     option_path = trim(option_path))
    
    ! copy solution back    
    do i = 1,rows
                
       x(i) =  node_val(solution, &
                        i)
    
    end do
    
    ! deallocate as needed
    
    deallocate(dnnz)    
    deallocate(rhs%val)
    deallocate(solution%val)
    
    call deallocate(matrix)
      
  end subroutine solve_via_copy_to_petsc_csr_matrix

! -----------------------------------------------------------------------------

end module solvers_module


