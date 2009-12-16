#include "confdefs.h"
#define DIMENSION 2

program test_laplacian
  ! A small program to solve laplacian psi = f
  ! 
  ! This tests and illustrates the use of fields and shape functions to
  ! solve finite element problems.
  use read_triangle
  use fields
  use FEtools
  use elements
  use sparse_tools
  use vtk_interfaces
  use transform_elements
  use sparsity_patterns
  use solvers
  use state_module
  use adapt_state_module 

  implicit none
#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif
  
  type(vector_field), target :: positions
  type(scalar_field), target :: psi, rhs, adj
  type(mesh_type) :: psi_mesh
  integer :: degree, quad_degree
  type(quadrature_type), target :: quad
  type(element_type), target :: psi_shape
  type(state_type) :: state
  interface 
     function rhs_func(X)
       ! A function which evaluates the right hand side at a number of
       ! points. Each column of X describes a set of points at which the
       ! right hand side is to be evaluated.
       real, dimension(:,:), intent(in) :: X
       real, dimension(size(X,2)) :: rhs_func
     end function rhs_func
  end interface
  interface
     subroutine set_adjoint(adjoint, positions)
       ! Compute the adjoint RHS.
       use fields_data_types
       type(scalar_field), intent(inout) :: adjoint
       type(vector_field), intent(in) :: positions
     end subroutine set_adjoint
  end interface
  interface 
     function solution(X)
       ! A function which evaluates the analytic solution at a number of
       ! points. Each column of X describes a set of points at which the
       ! right hand side is to be evaluated.
       real, dimension(:,:), intent(in) :: X
       real, dimension(size(X,2)) :: solution
     end function solution
  end interface
  interface 
     function adjoint(X)
       real, dimension(:,:), intent(in) :: X
       real, dimension(size(X,2)) :: adjoint
     end function adjoint
  end interface
  ! Arguments for handling the command line
  character(len=256) :: filename, degree_buffer
  integer :: status

#ifdef HAVE_PETSC
  integer :: ierr
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

  call get_command_argument(1, value=filename, status=status)
  
  select case(status)
  case(1:)
     call usage
     stop
  case(:-1)
     write(0,*) "Warning: truncating filename"
  end select

  call get_command_argument(2, value=degree_buffer, status=status)
  
  select case(status)
  case(1:)
     ! No degree specified.
     degree=1

  case default
     read(degree_buffer, *, iostat=status) degree

     if (status/=0) then
        write (0,*) trim(degree_buffer)//" is not an integer."
        call usage
        stop
     end if

  end select

  call get_command_argument(3, value=degree_buffer, status=status)
  
  select case(status)
  case(1:)
     ! No quadrature degree specified.
     quad_degree=2*degree

  case default
     read(degree_buffer, *, iostat=status) quad_degree

     if (status/=0) then
        write (0,*) trim(degree_buffer)//" is not an integer."
        call usage
        stop
     end if

  end select

  positions=read_triangle_files(filename, quad_degree)
  call insert(state, positions, "Coordinate")

  ! Shape functions for psi
  psi_shape=make_element_shape(loc=mesh_dim(positions)+1, dimension=mesh_dim(positions), degree=degree, &
   & quad=positions%mesh%shape%quadrature)
  psi_mesh=make_mesh(model=positions%mesh, shape=psi_shape)
  call insert(state, psi_mesh, "Mesh")
  call allocate(psi, psi_mesh, "Forward")
  call allocate(rhs, psi_mesh, "Right-hand side")
  call allocate(adj, psi_mesh, "Adjoint")
  
  ! Do the actual finite element calculation.
  call set_rhs(rhs, positions, rhs_func)
  call run_model(positions, psi, rhs)
  call set_adjoint(rhs, positions)
  call run_model(positions, adj, rhs)

  call insert(state, psi, "Forward")
  call insert(state, adj, "Adjoint")

  call set_constants(state, solution, adjoint)

  call analyse_output(positions, psi, solution)

  if (degree<=2.0) then
     ! Output to a vtk file.
     call vtk_write_state(trim(filename), index=1, state=(/state/))
  end if

contains

  subroutine set_constants(state, solution, adjoint)
    type(state_type), intent(inout) :: state
    interface 
       function solution(X)
         real, dimension(:,:), intent(in) :: X
         real, dimension(size(X,2)) :: solution
       end function solution
    end interface
    interface 
       function adjoint(X)
         real, dimension(:,:), intent(in) :: X
         real, dimension(size(X,2)) :: adjoint
       end function adjoint
    end interface

    type(scalar_field), pointer :: psi, adj
    type(vector_field), pointer :: positions
    real, dimension(1) :: tmp

    psi => extract_scalar_field(state, "Forward")
    adj =>  extract_scalar_field(state, "Adjoint")
    positions => extract_vector_field(state, "Coordinate")

    tmp = (solution(reshape(node_val(positions, 1), (/positions%dim, 1/))) - node_val(psi, 1))
    psi%val = psi%val + tmp(1)
    tmp = (adjoint(reshape(node_val(positions, 1), (/positions%dim, 1/))) - node_val(adj, 1))
    adj%val = adj%val + tmp(1)

  end subroutine set_constants
  
  subroutine run_model(positions, psi, rhs)
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: psi
    type(scalar_field), intent(in) :: rhs

    ! We form and solve the equation A*psi=b
    type(csr_matrix) :: A
    type(csr_sparsity) :: A_sparsity
    integer :: ele
    
    ! Calculate the sparsity of A based on the connectivity of psi.
    A_sparsity=make_sparsity(psi%mesh, psi%mesh, name='LaplacianSparsity')
    call allocate(A, A_sparsity)
    
    ! Assemble A element by element.
    do ele=1, element_count(psi)
       call assemble_element_contribution(A, positions, psi, ele) 
    end do

    ! It is necessary to fix the value of one node in the solution.
    ! We choose the first node of the first element.
    call set(A, 1, 1, INFINITY)

    psi%options%abs_error = 1.0e-8
    psi%options%max_its = 10000
    call zero(psi)
    call petsc_solve(psi, A, rhs)
  end subroutine run_model

  subroutine assemble_element_contribution(A, positions, psi, ele) 
    type(csr_matrix), intent(inout) :: A
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: psi
    integer, intent(in) :: ele

    ! Locations of nodes.
    real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
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
    type(element_type), pointer :: shape_psi, shape_X

    ele_psi=>ele_nodes(psi, ele)
    shape_psi=>ele_shape(psi, ele)
    shape_X=>ele_shape(positions, ele)

    ! Locations of local vertices.
    X_ele=ele_val(positions, ele)

    ! Locations of quadrature points.
    X_quad=ele_val_at_quad(positions, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(X_ele, shape_X, shape_psi, dm_t=dshape_psi,&
         & detwei=detwei)

    ! Matrix entry:
    call addto(A, ele_psi, ele_psi, &
         -dshape_dot_dshape(dshape_psi, dshape_psi, detwei))

  end subroutine assemble_element_contribution

  subroutine set_rhs(rhs, positions, rhs_func)
    type(scalar_field), intent(inout) :: rhs
    type(vector_field), intent(in) :: positions
    interface 
      function rhs_func(X)
        ! A function which evaluates the right hand side at a number of
        ! points. Each column of X describes a set of points at which the
        ! right hand side is to be evaluated.
        real, dimension(:,:), intent(in) :: X
        real, dimension(size(X,2)) :: rhs_func
      end function rhs_func
    end interface
    ! Locations of nodes.
    real, dimension(positions%dim,ele_loc(positions,1)) :: X_ele
    ! Locations of quadrature points
    real, dimension(positions%dim,ele_ngi(positions,1)) :: X_quad
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,1)) :: detwei    
    ! Node numbers of rhs element.
    integer, dimension(:), pointer :: ele_rhs
    ! Shape functions.
    type(element_type), pointer :: shape_psi, shape_X
    integer :: ele

    call zero(rhs)

    do ele=1,ele_count(rhs)
      ele_rhs=>ele_nodes(rhs, ele)
      shape_psi=>ele_shape(rhs, ele)
      shape_X=>ele_shape(positions, ele)
      X_ele=ele_val(positions, ele)
      X_quad=ele_val_at_quad(positions, ele)
      call transform_to_physical(X_ele, shape_X, shape_psi, detwei=detwei)
      call addto(rhs, ele_rhs, shape_rhs(shape_psi, detwei * rhs_func(X_quad)))
    end do
  end subroutine set_rhs
  
  subroutine analyse_output(positions, psi, solution)
    type(vector_field), target :: positions
    type(scalar_field), target :: psi
    interface 
       function solution(X)
         ! A function which evaluates the analytic solution at a number of
         ! points. Each column of X describes a set of points at which the
         ! right hand side is to be evaluated.
         real, dimension(:,:), intent(in) :: X
         real, dimension(size(X,2)) :: solution
       end function solution
    end interface
    
    type(vector_field) :: newpos
    real, dimension(:, :), allocatable :: val
    ! Total error in the solution.
    real :: error
    
    integer :: ele

    call allocate(newpos, positions%dim, psi%mesh, "Coordinate")
    call remap_field(from_field=positions, to_field=newpos)
    
    error=0
    do ele=1,element_count(positions)
       error=error+element_error(newpos, psi, solution, ele)
    end do
    
    write(0,*) "Degree = ", psi%mesh%shape%degree
    write(0,*) "Quad Degree = ", psi%mesh%shape%quadrature%degree
    write(0,'(a, i0, e22.8)') "Nodes, Error = ", node_count(psi), error
    
  end subroutine analyse_output

  function element_error(positions, psi, solution, ele)
    real :: element_error
    type(vector_field), target :: positions
    type(scalar_field), target :: psi
    interface 
       function solution(X)
         ! A function which evaluates the analytic solution at a number of
         ! points. Each column of X describes a set of points at which the
         ! right hand side is to be evaluated.
         real, dimension(:,:), intent(in) :: X
         real, dimension(size(X,2)) :: solution
       end function solution
    end interface
    integer, intent(in) :: ele

    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei    
    ! Locations of nodes.
    real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
    ! Shape functions.
    type(element_type), pointer :: shape_X
    real, dimension(1) :: offset
    
    shape_X=>ele_shape(positions, ele)

    ! Locations of local vertices.
    X_ele=ele_val(positions, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(X_ele, shape_X, detwei=detwei)

    ! Offset for zero pressure node
    offset=-solution(spread(node_val(positions,1),2,1))
    
    element_error=dot_product(detwei, &
         abs(ele_val_at_quad(psi, ele) &
         &   - solution(ele_val_at_quad(positions, ele)) &
             - offset(1) ))

  end function element_error

  subroutine usage
    
    write (0,*) "usage: test_laplacian <triangle_file_name> [<degree> [<quad_degree>]]"
    
  end subroutine usage

end program test_laplacian

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
  
  !rhs_func=-8.0*PI**2*cos(X(X_,:)*(2.0*PI))*cos(X(Y_,:)*(2.0*PI))
  rhs_func = -121 * PI * PI * cos(11 * PI * X(X_,:))
  rhs_func = rhs_func + -1*PI*sin(PI*(2*X(Y_,:) + 1)/2)
end function rhs_func

function solution(X)
  ! Analytic solution of scheme at X.
  use fetools
  implicit none
  real, dimension(:,:), intent(in) :: X
  real, dimension(size(X,2)) :: solution
  real, parameter :: PI=4.0*atan(1.0)
  
  !solution=cos(X(X_,:)*(2.0*PI))*cos(X(Y_,:)*(2.0*PI))
  solution = cos(11*PI * X(X_, :))
  solution = solution + sin(PI*(2*X(Y_, :) + 1)/2)/PI
end function solution

function adjoint(X)
  ! Analytic adjoint of scheme at X.
  use fetools
  implicit none
  real, dimension(:,:), intent(in) :: X
  real, dimension(size(X,2)) :: adjoint
  real, parameter :: PI=4.0*atan(1.0)
  
  adjoint = 1 - X(Y_, :)
end function adjoint

subroutine set_adjoint(adjoint, positions)
  use fields
  use fetools
  use transform_elements
  type(scalar_field), intent(inout) :: adjoint
  type(vector_field), intent(in) :: positions
  ! Locations of nodes.
  real, dimension(positions%dim,ele_loc(positions,1)) :: X_ele
  ! Locations of quadrature points
  real, dimension(positions%dim,ele_ngi(positions,1)) :: X_quad
  ! Derivatives of shape function:
  real, dimension(ele_loc(adjoint,1), &
       ele_ngi(adjoint,1), positions%dim) :: dshape_adjoint
  ! Coordinate transform * quadrature weights.
  real, dimension(ele_ngi(positions,1)) :: detwei    
  ! Node numbers of adjoint element.
  integer, dimension(:), pointer :: ele_adjoint
  ! Shape functions.
  type(element_type), pointer :: shape_psi, shape_X
  real, dimension(DIMENSION, ele_ngi(positions, 1)) :: v
  integer :: ele

  call zero(adjoint)

  forall(i=1:ele_ngi(positions,1))
    v(:, i) = (/0, 1/)
  end forall

  do ele=1,ele_count(adjoint)
    ele_adjoint=>ele_nodes(adjoint, ele)
    shape_psi=>ele_shape(adjoint, ele)
    shape_X=>ele_shape(positions, ele)
    X_ele=ele_val(positions, ele)
    X_quad=ele_val_at_quad(positions, ele)
    call transform_to_physical(X_ele, shape_X, shape_psi, dn_t=dshape_adjoint, detwei=detwei)
    call addto(adjoint, ele_adjoint, dshape_dot_vector_rhs(dshape_adjoint, v, detwei))
  end do
end subroutine set_adjoint
