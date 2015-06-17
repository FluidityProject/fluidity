!! test matrix-free PETSc solve 
!! WE SOLVE MINUS LAPLACE EQUATION!
#include "fdebug.h"

  subroutine test_pressure_solve

    use unittest_tools
    use solvers
    use fields
    use state_module
    use elements
    use sparse_tools
    use mesh_files
    use vtk_interfaces
    use boundary_conditions
    use global_parameters, only: OPTION_PATH_LEN
    use free_surface_module
    implicit none

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"

    logical :: fail=.false., warn=.false.

    type(state_type) :: state,state_in
    type(vector_field), target:: positions
    type(scalar_field) :: psi, DistanceToTop
    type(mesh_type) :: psi_mesh
    type(mesh_type), pointer :: x_mesh
    type(element_type) :: psi_shape
    type(quadrature_type) :: quad
    integer :: quad_degree=4
    integer, parameter:: DIM=3
    real :: eps0
    real, dimension(DIM), parameter:: LENGTH=(/ 1.0, 1.0, 1.0 /)
    integer, dimension(DIM), parameter:: WAVENUMBER=(/ 2, 2 , 0 /)
    character(len=OPTION_PATH_LEN) solver_option_path
    integer :: unit,io1

    namelist/epsilon_data/eps0

    call set_global_debug_level(3)

    unit = free_unit()
    open(unit=unit, file="epsilon.dat", status='old', &
         iostat=io1)

    if(io1.ne.0) then
       ewrite(-1,*) 'Looked for ', "epsilon.dat"
       FLExit('Could not read from .dat file')
    end if

    read(unit, epsilon_data)
    close(unit) 

    positions=read_mesh_files("cube_unstructured", &
         quad_degree=QUAD_DEGREE, format="gmsh")

    !call vtk_read_state("cube-1_1.vtu", state_in, quad_degree)
    !positions=extract_vector_field(state_in,name="Coordinate")

    !positions=read_mesh_files("test_laplacian.1", &
    !     quad_degree=QUAD_DEGREE, format="gmsh")
    x_mesh => positions%mesh

    call insert(state, positions, name="Coordinate")
    call insert(state, positions%mesh, "Coordinate_mesh")

    ! Shape functions for psi
    assert(dim==mesh_dim(positions))
    quad=make_quadrature(vertices = dim+1, dim =dim, degree=quad_degree)

    psi_shape=make_element_shape(vertices = dim+1, dim =dim, degree=1, quad=quad)
    psi_mesh=make_mesh(model=positions%mesh, shape=psi_shape)

    call insert(state, psi_mesh, "Psi_Mesh")
    call allocate(psi, psi_mesh, "Psi")
    call allocate(DistanceToTop, psi_mesh, "DistanceToTop")
    call add_boundary_condition(DistanceToTop,"top","surface",(/1/) )

    call insert(state, psi, "Psi")
    call insert(state, DistanceToTop,"DistanceToTop")

    call run_model(state)

    call report_test("[id matrix solve]", fail, warn, "Solving Ix = b should yield x == b.")

  contains

    subroutine run_model(state)
      use global_parameters, only: PYTHON_FUNC_LEN
      use unittest_tools
      use solvers
      use boundary_conditions
      use fields
      use state_module
      use elements
      use sparse_tools
      use mesh_files
      use sparsity_patterns
      use boundary_conditions
      implicit none
      type(state_type), intent(inout) :: state

      type(vector_field), pointer :: positions
      type(scalar_field), pointer :: psi
      type(scalar_field) :: exact, error
      type(scalar_field), pointer :: topdis
      type(mesh_type) :: top_surface_mesh
      type(csr_matrix) :: vprolongator
      integer, dimension(:), pointer :: top_surface_node_list => null()
      integer, dimension(:), pointer :: top_surface_element_list => null()

      ! We form and solve the equation A*psi=rhs
      type(csr_sparsity) :: A_sparsity
      type(csr_matrix) :: A
      type(scalar_field) :: RHS
      integer :: ele, unit, stat

      character(len=PYTHON_FUNC_LEN) :: func, buffer
      logical :: file_exists

      ! Extract the required fields from state.
      positions=>extract_vector_field(state, "Coordinate")
      positions%val(1,:) = positions%val(1,:)/eps0
      positions%val(2,:) = positions%val(2,:)/eps0

      psi=>extract_scalar_field(state, "Psi")


      !=======================
      !Get Boundary conditions
      !=======================
      topdis => extract_scalar_field(state,"DistanceToTop",stat=stat)
      if(stat/=0) then
         FLExit('DistanceToTop is not present')
      end if      
      call get_boundary_condition(topdis,1,&
           surface_element_list=top_surface_element_list)
      call create_surface_mesh(top_surface_mesh, &
           top_surface_node_list, psi%mesh, &
           top_surface_element_list, name="PsiTopSurfaceMesh")
    !=======================

      !Calculate the sparsity of A based on the connectivity of psi.
      A_sparsity=make_sparsity(psi%mesh, psi%mesh, name='LaplacianSparsity')
      call allocate(A, A_sparsity)
      call zero(A)

      call get_laplacian(A,positions,psi)
      call set(A, top_surface_node_list(1),top_surface_node_list(1), INFINITY)

      call allocate(rhs, psi%mesh, "RHS")
      call zero(rhs)

      call allocate(error, psi%mesh, "Error")
      call zero(error)

      inquire(file="exact.py",exist=file_exists)  
      if (.not.file_exists) then
          ewrite(-1,*) "Looking for file", "exact.py"
          FLExit('Couldnt find file')
      end if
      unit=free_unit() 
      open(unit, file="exact.py", action="read",&
           & status="old")
      read(unit, '(a)', end=43) func
      ! Read all the lines of the file and put in newlines between them.
      do
         read(unit, '(a)', end=43) buffer
         func=trim(func)//achar(10)//trim(buffer)
      end do
43    func=trim(func)//achar(10)
      close(unit)

      call allocate(exact,psi_mesh,name='Exact')
      call set_from_python_function(exact, trim(func), positions, 0.0)

      solver_option_path=""
      call set_solver_options(solver_option_path, "my_solver", &
           ksptype=KSPPREONLY, PCTYPE="mg", &
           rtol=1.0e-100, atol=1.0e-30, max_its=1000, &
           start_from_zero=.true.)

      exact%val = exact%val - exact%val(top_surface_node_list(1)) 
      exact%val(top_surface_node_list(1)) = 0.0

      call mult(rhs%val,A,exact%val)

      call zero(psi)

      ! supplying the prolongator to petsc_solve makes 'mg' 
      ! use the vertical_lumping option
      vprolongator = &
           vertical_prolongator_from_free_surface(state, psi%mesh)

      call set_solver_options(psi, &
           ksptype="cg", pctype="mg", &
           atol=1.0e-100, rtol=1.0e-20, max_its=1000, &
           start_from_zero=.true.)

      ewrite(1,*) 'with vertical lumping and internal smoother'

      call petsc_solve(psi, A, rhs, prolongator=vprolongator, &
           & exact=exact, surface_node_list=top_surface_node_list)

      ewrite(1,*) 'with vertical lumping, no additive smoother'

      call petsc_solve(psi, A, rhs, prolongator=vprolongator, &
           & exact=exact)

      ewrite(1,*) 'without vertical lumping'

      call petsc_solve(psi, A, rhs, &
           & exact=exact)

      !error%val = psi%val-exact%val
      !call vtk_write_fields('test', 0, positions, positions%mesh, &
      !     sfields=(/ psi, rhs, exact, error /) )

    end subroutine run_model

    subroutine get_laplacian(A,positions,psi,H0,vertical)
      implicit none
      type(csr_matrix), intent(inout) :: A
      type(vector_field), intent(in) :: positions
      type(scalar_field), intent(in) :: psi
      type(scalar_field), intent(in), optional :: H0
      logical, intent(in), optional :: vertical
      !
      integer :: ele
      logical :: lvertical

      if(present(vertical)) then
         lvertical = vertical
      else
         lvertical = .false.
      end if

      do ele = 1, element_count(psi)
         if(present(H0)) then
            call assemble_laplacian_element_contribution(&
                 &A, positions, psi, ele, lvertical, H0)
         else
            call assemble_laplacian_element_contribution(&
                 &A, positions, psi, ele, lvertical)
         end if
      end do

    end subroutine get_laplacian

subroutine assemble_laplacian_element_contribution(A, positions, psi, ele, &
     vertical,H0) 
  use unittest_tools
  use solvers
  use fields
  use state_module
  use elements
  use sparse_tools
  use mesh_files
  implicit none
  type(csr_matrix), intent(inout) :: A
  type(vector_field), intent(in) :: positions
  type(scalar_field), intent(in) :: psi
  logical, intent(in) :: vertical
  type(scalar_field), optional, intent(in) :: H0
  integer, intent(in) :: ele
  
  ! Locations of quadrature points
  real, dimension(positions%dim,ele_ngi(positions,ele)) :: X_quad
  ! H0 at quadrature points
  real, dimension(ele_ngi(positions,ele)) :: H0_quad
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
  ! tensor
  real, dimension(positions%dim,positions%dim,ele_ngi(positions,ele)) :: &
       tau_quad
  integer :: i
  
  ele_psi=>ele_nodes(psi, ele)
  shape_psi=>ele_shape(psi, ele)
  
  !tau_quad = 0.0
  !tau_quad(3,3,:) = 1.0
  !if(.not.vertical .or. positions%dim==2) then
  !   tau_quad(1,1,:) = 1.0
  !   tau_quad(2,2,:) = 1.0
  !end if

  ! Locations of quadrature points.
  X_quad=ele_val_at_quad(positions, ele)
  
  ! value of H0 at quadrature points
  if(present(H0)) then
     H0_quad=ele_val_at_quad(H0, ele)
  else
     H0_quad = 1.0
  end if

  ! Transform derivatives and weights into physical space.
  call transform_to_physical(positions, ele, shape_psi, dshape=dshape_psi,&
       & detwei=detwei)
  
  ! Local assembly:
  !psi_mat=dshape_tensor_dshape(dshape_psi, tau_quad, dshape_psi, detwei)
  psi_mat=dshape_dot_dshape(dshape_psi, dshape_psi, detwei*H0_quad)
  
  ! Global assembly:
  call addto(A, ele_psi, ele_psi, psi_mat)
  
end subroutine assemble_laplacian_element_contribution

  end subroutine test_pressure_solve

