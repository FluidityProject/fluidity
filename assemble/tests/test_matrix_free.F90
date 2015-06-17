!! test matrix-free PETSc solve 
#include "fdebug.h"
subroutine test_matrix_free

  use unittest_tools
  use solvers
  use fields
  use state_module
  use elements
  use sparse_tools
  use mesh_files
  use vtk_interfaces
  implicit none
  
  logical :: fail=.false., warn=.false.

  type(state_type) :: state
  type(vector_field), target:: positions
  type(scalar_field) :: psi
  type(mesh_type) :: psi_mesh
  type(mesh_type), pointer :: x_mesh
  type(element_type) :: psi_shape
  type(quadrature_type) :: quad
  integer :: dim
  integer :: quad_degree=4

  interface 
     function rhs_func(X)       
       ! A function which evaluates the right hand side at a number of
       ! points. Each column of X describes a set of points at which the
       ! right hand side is to be evaluated.
       real, dimension(:,:), intent(in) :: X
       real, dimension(size(X,2)) :: rhs_func
     end function rhs_func
  end interface

  call set_global_debug_level(3)

  positions=read_mesh_files("data/cube_unstructured", &
       quad_degree=QUAD_DEGREE, format="gmsh")
  x_mesh => positions%mesh
  
  call insert(state, positions, name="Coordinate")
  call insert(state, positions%mesh, "Coordinate_mesh")

  ! Shape functions for psi
  dim=mesh_dim(positions)
  quad=make_quadrature(vertices = dim+1, dim =dim, degree=quad_degree)

  psi_shape=make_element_shape(vertices = dim+1, dim =dim, degree=2, quad=quad)
  psi_mesh=make_mesh(model=positions%mesh, shape=psi_shape)

  call insert(state, psi_mesh, "Psi_Mesh")
  call allocate(psi, psi_mesh, "Psi")
  call insert(state, psi, "Psi")

  call run_model(state, rhs_func)

  call vtk_write_fields('test', 0, positions, positions%mesh, &
       sfields=(/ psi /) )
  
  call report_test("[id matrix solve]", fail, warn, "Solving Ix = b should yield x == b.")

end subroutine test_matrix_free

subroutine run_model(state, rhs_func)
  use unittest_tools
  use solvers
  use fields
  use state_module
  use elements
  use sparse_tools
  use mesh_files
  use sparsity_patterns
  implicit none
  type(state_type), intent(inout) :: state
  interface 
     function rhs_func(X)       
         ! A function which evaluates the right hand side at a number of
       ! points. Each column of X describes a set of points at which the
       ! right hand side is to be evaluated.
       real, dimension(:,:), intent(in) :: X
       real, dimension(size(X,2)) :: rhs_func
     end function rhs_func
  end interface
  
  type(vector_field), pointer :: positions
  type(scalar_field), pointer :: psi
  
  ! We form and solve the equation A*psi=rhs
  type(csr_matrix) :: A
  type(csr_sparsity) :: A_sparsity
  type(scalar_field) :: RHS
  integer :: ele
  
  ! Extract the required fields from state.
  positions=>extract_vector_field(state, "Coordinate")
  psi=>extract_scalar_field(state, "Psi")
  
  ! Calculate the sparsity of A based on the connectivity of psi.
  A_sparsity=make_sparsity(psi%mesh, psi%mesh, name='LaplacianSparsity')
  call allocate(A, A_sparsity)
  
  call zero(A)
  
  call allocate(rhs, psi%mesh, "RHS")
  call zero(rhs)
  
  ! Assemble A element by element.
  do ele=1, element_count(psi)
     call assemble_element_contribution(A, rhs, positions, psi, rhs_func,&
          & ele)  
  end do

  ewrite(1,*) 'sum A', sum(A%val), maxval(A%val)
  
  ! It is necessary to fix the value of one node in the solution.
  ! We choose node 1.
  call set(A, 1, 1, INFINITY)
  
  ewrite(1,*) 'sum A', sum(A%val), maxval(A%val)

  psi%options%abs_error = 1.0e-8
  
  psi%options%max_its = 10000

  !call set_solver_options(psi, &
  !  ksptype="cg", pctype="sor", atol=1.0e-8, rtol=1.0e-8, max_its=1000, &
  !  start_from_zero=.true.)
  
  call zero(psi)

  ewrite(1,*) maxval(rhs%val)
  
  !call petsc_solve(psi,A,rhs)
  call petsc_solve_matrix_free(psi, A, rhs)
  
  call deallocate(A)
  call deallocate(A_sparsity)
  call deallocate(rhs)
  
end subroutine run_model

subroutine assemble_element_contribution(A, rhs, positions, psi, rhs_func&
       &, ele) 
  use unittest_tools
  use solvers
  use fields
  use state_module
  use elements
  use sparse_tools
  use mesh_files
  implicit none
  type(csr_matrix), intent(inout) :: A
  type(scalar_field), intent(inout) :: rhs
  type(vector_field), intent(in) :: positions
  type(scalar_field), intent(in) :: psi
  interface 
     function rhs_func(X)       
       ! A function which evaluates the right hand side at a number of
       ! points. Each column of X describes a set of points at which the
       ! right hand side is to be evaluated.
       real, dimension(:,:), intent(in) :: X
       real, dimension(size(X,2)) :: rhs_func
     end function rhs_func
  end interface
  integer, intent(in) :: ele
  
  ! Locations of quadrature points
  real, dimension(positions%dim,ele_ngi(positions,ele)) :: X_quad
  ! Derivatives of shape function:
  real, dimension(ele_loc(psi,ele), &
       ele_ngi(psi,ele), positions%dim) :: dshape_psi
  ! Coordinate transform * quadrature weights.
  real, dimension(ele_ngi(positions,ele)) :: detwei    
  ! Node numbers of psi element.
  integer, dimension(:), pointer :: ele_psi
  ! Shape functions.
  type(element_type), pointer :: shape_psi
  ! Local Laplacian matrix 
  real, dimension(ele_loc(psi, ele), ele_loc(psi, ele)) :: psi_mat
  ! Local right hand side.
  real, dimension(ele_loc(psi, ele)) :: lrhs
  
  ele_psi=>ele_nodes(psi, ele)
  shape_psi=>ele_shape(psi, ele)
  
  ! Locations of quadrature points.
  X_quad=ele_val_at_quad(positions, ele)
  
  ! Transform derivatives and weights into physical space.
  call transform_to_physical(positions, ele, shape_psi, dshape=dshape_psi,&
       & detwei=detwei)
  
  ! Local assembly:
  psi_mat=dshape_dot_dshape(dshape_psi, dshape_psi, detwei)
  
  lrhs=shape_rhs(shape_psi, rhs_func(X_quad)*detwei)
  
  ! Global assembly:
  call addto(A, ele_psi, ele_psi, psi_mat)
  
  call addto(rhs, ele_psi, lrhs)
  
end subroutine assemble_element_contribution

function rhs_func(X)
  ! Right hand side function for laplacian operator.
  !
  ! Each column of X is interpretted as a position at which RHS should be
  ! evaluated. 
  use fetools
  implicit none
  real, dimension(:,:), intent(in) :: X
  real, dimension(size(X,2)) :: rhs_func
  real, parameter :: PI=4.0*atan(1.0)
  integer :: i,dim

  dim=size(X,1)

  rhs_func=2*(dim-1+0.25)*PI*cos(0.5*PI*X(1,:))

  do i=2,dim
     rhs_func=rhs_func*cos(PI*X(i,:))
  end do
    
  rhs_func=rhs_func+0.5*PI*sin(0.5*PI*X(1,:))

end function rhs_func

subroutine petsc_solve_matrix_free(x, matrix, rhs)
  !!< Solve a linear system
  use fields
  use sparse_tools
  use solvers
  use petsc_tools
  use matrix_free_solvers
  implicit none
  type(scalar_field), intent(inout) :: x
  type(scalar_field), intent(in) :: rhs
  type(csr_matrix), intent(in) :: matrix
  integer, dimension(:), allocatable :: petsc_numbering
#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#endif
  KSP :: ksp 
  Mat :: Amat, Pmat
  Vec :: y, b
  PCType :: pctype = PCNONE
  KSPType :: ksptype = KSPCG
  PC :: pc
  PetscReal rtol, atol, dtol
  PetscInt max_its
  PetscErrorCode :: ierr
  KSPConvergedReason :: reason

  integer literations,i
  logical :: lstartfromzero=.true.

  assert(size(x%val)==size(rhs%val))
  assert(size(x%val)==size(matrix,2))
  assert(size(rhs%val)==size(matrix,1))
  
  !  call allocate(petsc_numbering, &
  !     size(rhs%val), 1)

  allocate(petsc_numbering(size(rhs%val)) )
  petsc_numbering = -1
  do i = 1, size(rhs%val)
     petsc_numbering(i) = i
  end do

  !set up preconditioner matrix
  Pmat = csr2petsc(matrix)!,petsc_numbering)
  !Amat = csr2petsc(matrix,petsc_numbering)

  !set up matrix-free matrix
  call petsc_mult_setup(Pmat)
  call MatCreateShell(MPI_COMM_SELF,size(rhs%val),size(rhs%val), &
       size(rhs%val),size(rhs%val),PETSC_NULL_INTEGER,Amat,ierr)
  call MatShellSetOperation(AMat,MATOP_MULT,petsc_mult,ierr)

  !set up KSP
  call KSPCreate(MPI_COMM_SELF, ksp, ierr)
  call KSPSetOperators(ksp, amat, pmat, DIFFERENT_NONZERO_PATTERN, ierr)
  call KSPGetPC(ksp, pc, ierr)
  pctype=PCNONE
  ksptype=KSPCG
  call PCSetType(pc, pctype, ierr)
  call KSPSetType(ksp, ksptype, ierr)
  dtol=PETSC_DEFAULT_DOUBLE_PRECISION
  max_its = 1000
  rtol =1.0e-7
  atol =1.0e-7

  call KSPSetInitialGuessNonzero(ksp, PETSC_FALSE, ierr)
  call KSPSetUp(ksp, ierr)

  call VecCreateSeq(MPI_COMM_SELF,size(rhs%val), b, ierr)

  call VecDuplicate(b, y, ierr)

  ! copy array into PETSc vecs
  call VecSetValues(b, size(rhs%val), &
       petsc_numbering(1:size(rhs%val)), &
       rhs%val( 1:size(rhs%val) ), INSERT_VALUES, ierr)

  call VecSetValues(y, size(rhs%val), &
       petsc_numbering(1:size(rhs%val)), &
       rhs%val( 1:size(rhs%val) )*0.0, INSERT_VALUES, ierr)
     
  ! the solve and convergence check
  call KSPSolve(ksp, b, y, ierr)
  call KSPGetConvergedReason(ksp, reason, ierr)
  call KSPGetIterationNumber(ksp, literations, ierr)

  ewrite(1,*) reason, literations
        
  ! Copy back the result using the petsc numbering:
  !call petsc2field(y, petsc_numbering, x)
  
  ! destroy all PETSc objects and the petsc_numbering
  !call petsc_solve_destroy(y, A, b, ksp, petsc_numbering)
  
end subroutine petsc_solve_matrix_free
