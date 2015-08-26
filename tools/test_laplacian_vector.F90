#include "confdefs.h"
program test_laplacian_vector
  ! A small program to solve laplacian psi = f
  ! 
  ! This tests and illustrates the use of fields and shape functions to
  ! solve finite element problems.
  !
  ! The analytical solution provided is valid on a 1x1 square or 
  ! 1x1x1 cube with the boundary conditions 1 and 2 applied on the
  ! sides of the first coordinates direction (e.g. test_laplacian.poly
  ! - create a mesh with 'triangle -a0.01 -e test_laplacian.poly' )
  use mesh_files
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
  use boundary_conditions
  use fldebug

#ifdef HAVE_PETSC_MODULES
  use petsc
#endif
  implicit none
#include "petsc_legacy.h"
  
  character, parameter:: NEWLINE_CHAR=achar(10)
  character(len=*), parameter:: BC_PYTHON_FUNCTION= &
     "def val(X, t):"//NEWLINE_CHAR// &
     "  import math"//NEWLINE_CHAR// &
     "  return math.cos(math.pi*2*X[1])"
  real, parameter:: BC_CONST_VALUE=1.0
  
  type(vector_field), target :: positions
  type(vector_field), target :: psi
  type(mesh_type) :: psi_mesh
  integer :: degree, quad_degree, msh_dim, vertices, psi_dim
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
     function solution(X)
       ! A function which evaluates the analytic solution at a number of
       ! points. Each column of X describes a set of points at which the
       ! right hand side is to be evaluated.
       real, dimension(:,:), intent(in) :: X
       real, dimension(size(X,2)) :: solution
     end function solution
  end interface
  character(len=256) :: filename

#ifdef HAVE_PETSC
  integer :: ierr
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

  call set_debug_level(3)

  call python_init

  call read_command_line(filename, psi_dim, degree, quad_degree)

  positions=read_mesh_files(filename, quad_degree=quad_degree)

  call insert(state, positions, "Coordinate")

  call insert(state, positions%mesh, "Coordinate_mesh")

  ! Shape functions for psi
  msh_dim=mesh_dim(positions)
  vertices=msh_dim+1
  quad=make_quadrature(vertices, msh_dim, degree=quad_degree)

  psi_shape=make_element_shape(vertices, msh_dim, degree=degree, quad=quad)
  psi_mesh=make_mesh(model=positions%mesh, shape=psi_shape)

  call insert(state, psi_mesh, "Psi_Mesh")

  call allocate(psi, psi_dim, psi_mesh, "Psi")
  
  ! Add two boundary conditions, specify their name, something to identify
  ! the type of b.c. with, and to which part of the surface to apply
  ! this b.c., i.e to all surface elements with a boundary_ids in the list
  call add_boundary_condition(psi, name="my_first_boundary", &
    type="dirichlet", boundary_ids=(/ 1 /) )
  call add_boundary_condition(psi, name="my_second_boundary", &
    type="dirichlet", boundary_ids=(/ 2 /) )
  
  call set_boundary_condition_values(psi, positions)
  
  call insert(state, psi, "Psi")
  
  ! Do the actual finite element calculation.
  call run_model(state, rhs_func)

  !call analyse_output(positions, psi, solution)

  if (degree<=2) then
     ! Output to a vtk file.
     call vtk_write_state(trim(filename), index=1, state=(/ state /), &
        model="Psi_Mesh")
  end if

  call python_end

contains

  subroutine set_boundary_condition_values(psi, positions)
    type(vector_field), intent(inout):: psi
    type(vector_field), intent(in):: positions
      
    type(mesh_type), pointer:: bc_surface_mesh
    type(vector_field) bc_field
    type(vector_field) bc_positions
    integer, dimension(:), pointer:: bc_surface_elements
    
    integer :: i
    type(scalar_field) :: bc_comp

    ! First boundary condition, let's keep it simple and constant.
    ! ask for a mesh of this part of the surface only
    call get_boundary_condition(psi, name="my_first_boundary", &
      surface_mesh=bc_surface_mesh)
      
    ! allocate a field on it to set b.c. values on
    call allocate(bc_field, psi%dim, bc_surface_mesh, name='value')
    
    ! set it to some constant
    do i = 1, psi%dim
      bc_comp = extract_scalar_field(bc_field, i)
      call set(bc_comp, BC_CONST_VALUE)
    end do
    
    ! insert it to the boundary condition:
    call insert_surface_field(psi, "my_first_boundary", bc_field)
    
    ! deallocate our reference
    call deallocate(bc_field)
    
    ! For the second boundary condition, let's do something more audacious:
    ! we're gonna initialise it with a python function!
    ! we start the same, but also ask for list of surface_elements
    call get_boundary_condition(psi, name="my_second_boundary", &
      surface_mesh=bc_surface_mesh, surface_element_list=bc_surface_elements)
    call allocate(bc_field, psi%dim, bc_surface_mesh, name='value')
    
    ! need a positions field on the surface mesh only
    ! note: this is a dim vector field on a dim-1 mesh
    call allocate(bc_positions, positions%dim, bc_surface_mesh)
    call remap_field_to_surface(positions, bc_positions, bc_surface_elements)
    
    ! initialise with python function
    do i = 1, psi%dim
      bc_comp = extract_scalar_field(bc_field, i)
      call set_from_python_function(bc_comp, BC_PYTHON_FUNCTION, bc_positions, 0.0)
    end do
    call deallocate(bc_positions)
    
    ! same as above insert it to the boundary condition and deallocate our reference:
    call insert_surface_field(psi, "my_second_boundary", bc_field)
    call deallocate(bc_field)
    
  end subroutine set_boundary_condition_values
  
  subroutine run_model(state, rhs_func)
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
    type(vector_field), pointer :: psi
    
    ! We form and solve the equation A*psi=rhs
    type(block_csr_matrix) :: A
    type(csr_sparsity) :: A_sparsity
    type(vector_field) :: RHS
    integer :: ele
    
    ! Extract the required fields from state.
    positions=>extract_vector_field(state, "Coordinate")
    psi=>extract_vector_field(state, "Psi")

    ! Calculate the sparsity of A based on the connectivity of psi.
    A_sparsity=make_sparsity(psi%mesh, psi%mesh, name='LaplacianSparsity')
    call allocate(A, A_sparsity, (/psi%dim,psi%dim/))
    
    call zero(A)
    
    call allocate(rhs, psi%dim, psi%mesh, "RHS")
    call zero(rhs)

    ! Assemble A element by element.
    do ele=1, element_count(psi)
       call assemble_element_contribution(A, rhs, positions, psi, rhs_func,&
            & ele)  
    end do
      
    !call assemble_boundary_contribution(rhs, positions, psi, "my_first_boundary")

    !call assemble_boundary_contribution(rhs, positions, psi, "my_second_boundary")

    ! It is necessary to fix the value of one node in the solution.
    ! We choose node 1.
    !call set(A, 1, 1, INFINITY)

    call apply_dirichlet_conditions(A, rhs, psi)

    call zero(psi)

    !call matrix2file("A", A)
    
    call set_solver_options(psi, ksptype='cg', pctype='mg', rtol=1e-7)

    call petsc_solve(psi, A, rhs)
    
    call vtk_write_fields('test', 0, positions, positions%mesh, &
      vfields=(/ psi /) )
        
    call deallocate(A)
    call deallocate(A_sparsity)
    call deallocate(rhs)

  end subroutine run_model

  subroutine assemble_element_contribution(A, rhs, positions, psi, rhs_func&
       &, ele) 
    type(block_csr_matrix), intent(inout) :: A
    type(vector_field), intent(inout) :: rhs
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: psi
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
    ! Local Laplacian matrix 
    real, dimension(ele_loc(psi, ele), ele_loc(psi, ele)) :: psi_mat
    ! Local right hand side.
    real, dimension(ele_loc(psi, ele)) :: lrhs

    integer :: i

    ele_psi=>ele_nodes(psi, ele)
    shape_psi=>ele_shape(psi, ele)

    ! Locations of quadrature points.
    X_quad=ele_val_at_quad(positions, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_psi, dshape=dshape_psi,&
         & detwei=detwei)

    ! Local assembly:
    psi_mat=dshape_dot_dshape(dshape_psi, dshape_psi, detwei)

    lrhs=-shape_rhs(shape_psi, rhs_func(X_quad)*detwei)

    ! Global assembly:
    do i = 1, psi%dim
      call addto(A, i, i, ele_psi, ele_psi, psi_mat)

      call addto(rhs, i, ele_psi, lrhs)
    end do

  end subroutine assemble_element_contribution
    
  subroutine assemble_boundary_contribution(rhs, positions, psi, bc_name)
    type(vector_field), intent(inout) :: rhs
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: psi
    character(len=*), intent(in):: bc_name
    
    type(vector_field), pointer:: bc_surface_field
    integer, dimension(:), pointer:: surface_element_list
    type(element_type), pointer:: psi_face_shape
    ! note that we assume all shapes to be the same in each element
    real, dimension(face_ngi(positions,1)) :: detwei_face
    real, dimension(psi%dim, face_loc(psi,1)) :: bc_val_face
    real, dimension(face_loc(psi,1), face_loc(psi,1)) :: face_mat
    integer, dimension(face_loc(psi,1)) :: psi_face_nodes
      
    integer ele, face, i
    type(scalar_field) :: rhs_comp
      
    ! pull out the b.c value field again:
    bc_surface_field => extract_surface_field(psi, bc_name, 'value')
    ! retrieve the list of surface elements/faces where this bc is applied
    call get_boundary_condition(psi, bc_name, &
      surface_element_list=surface_element_list)
    
    ! now do the surface integral:
    
    ! loop over the elements of the surface field
    ! NOTE: for a 3D (resp. 2D) problem each element of the 2D (resp. 1D)
    ! surface field is a face in the 3D (resp. 2D)  mesh(es) of psi and 
    ! positions
    do ele=1, ele_count(bc_surface_field)
      
      ! get the values in the 2D element in the usual way:
      bc_val_face=ele_val(bc_surface_field, ele)
      
      ! to acces value in the 3D fields, we need to know the face number
      face=surface_element_list(ele)
      
      ! given those calculate the quadrature weights
      call transform_facet_to_physical(positions, face, detwei_face)

      ! integral over the face of the form \int N_i N_j
      ! where N_i and N_j are shape functions of psi
      psi_face_shape => face_shape(psi, face)
      face_mat=shape_shape(psi_face_shape, psi_face_shape, detwei_face)

      ! global node numbers of nodes of this face in psi%mesh
      psi_face_nodes=face_global_nodes(psi, face)
      
      do i = 1, psi%dim
        rhs_comp = extract_scalar_field(rhs, i)
        call addto(rhs_comp, psi_face_nodes, matmul(face_mat, bc_val_face(i,:)))
      end do
    end do
    
  end subroutine assemble_boundary_contribution

  subroutine read_command_line(filename, solution_dimension, degree, quad_degree)
    ! Read the input filename, degree and quadrature degree on the command
    ! line.
    character(len=*), intent(out) :: filename
    integer, intent(out) :: solution_dimension, degree, quad_degree
    character(len=256) :: buffer
    integer :: status
    
    call get_command_argument(1, value=filename, status=status)
  
    select case(status)
    case(1:)
       call usage
       stop
    case(:-1)
       write(0,*) "Warning: truncating filename"
    end select
    
    call get_command_argument(2, value=buffer, status=status)
    
    select case(status)
    case(1:)
       ! No solution_dimension specified.
       solution_dimension=1
       
    case default
       read(buffer, *, iostat=status) solution_dimension
       
       if (status/=0) then
          write (0,*) trim(buffer)//" is not an integer."
          call usage
          stop
       end if

       if(solution_dimension > 3) then
          write (0,*) "solution_dimension "//trim(buffer)//" requested."
          write (0,*) "Vector dimensions greater than 3 not supported."
          call usage
          stop
       end if
       
       if(solution_dimension < 1) then
          write (0,*) "solution_dimension "//trim(buffer)//" requested."
          write (0,*) "Vector dimensions less than 1 are not supported."
          call usage
          stop
       end if
       
    end select
    
    call get_command_argument(3, value=buffer, status=status)
    
    select case(status)
    case(1:)
       ! No degree specified.
       degree=1
       
    case default
       read(buffer, *, iostat=status) degree
       
       if (status/=0) then
          write (0,*) trim(buffer)//" is not an integer."
          call usage
          stop
       end if
       
    end select

    call get_command_argument(4, value=buffer, status=status)
    
    select case(status)
    case(1:)
       ! No quadrature degree specified.
       quad_degree=2*degree
       
    case default
       read(buffer, *, iostat=status) quad_degree
       
       if (status/=0) then
          write (0,*) trim(buffer)//" is not an integer."
          call usage
          stop
       end if
       
    end select

  end subroutine read_command_line
  
  subroutine analyse_output(positions, psi, solution)
    type(vector_field), target :: positions
    type(vector_field), target :: psi
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
      
    call field2file("Coordinate", newpos)
    call field2file("Psi", psi)

    write(0,*) "Degree = ", psi%mesh%shape%degree
    write(0,*) "Quad Degree = ", psi%mesh%shape%quadrature%degree
    write(0,'(a, i0, e22.8)') "Nodes, Error = ", node_count(psi), error
    
  end subroutine analyse_output

  function element_error(positions, psi, solution, ele)
    real :: element_error
    type(vector_field), target :: positions
    type(vector_field), target :: psi
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
    ! Offset for zero pressure node
    real, dimension(1) :: offset

    ! Transform weights into physical space.
    call transform_to_physical(positions, ele, detwei=detwei)

    offset=-solution(spread(node_val(positions,1),2,1))
    
    element_error=dot_product(detwei, &
         abs(ele_val_at_quad(psi, ele, 1) &
         &   - solution(ele_val_at_quad(positions, ele)) &
             - offset(1) ))

  end function element_error

  subroutine usage
    
    write (0,*) "usage: test_laplacian <mesh_file_name> [<solution_dimension> <degree> [<quad_degree>]]"
    
  end subroutine usage

end program test_laplacian_vector

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

function solution(X)
  ! Analytic solution of scheme at a set of points, X.
  use fetools
  implicit none
  real, dimension(:,:), intent(in) :: X
  real, dimension(size(X,2)) :: solution
  real, parameter :: PI=4.0*atan(1.0)
  integer :: i,dim
  
  dim=size(X,1)

  solution=-2*cos(0.5*PI*X(1,:))/PI

  do i=2,dim
     solution=solution*cos(PI*X(i,:))
  end do
    
  solution=solution-2*sin(0.5*PI*X(1,:))/PI

end function solution
