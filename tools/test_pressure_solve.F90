  !! test pressure solvers on large aspect ratio domains
  !! WE SOLVE MINUS LAPLACE EQUATION -- i.e. geometric Laplacian
  !! This is to make a positive definite matrix (instead of negative)
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
    use global_parameters, only: OPTION_PATH_LEN, PYTHON_FUNC_LEN
    use free_surface_module
    use FLDebug
#ifdef HAVE_PETSC_MODULES
  use petsc
#endif
  implicit none
#include "petsc_legacy.h"

    type(state_type) :: state
    type(vector_field), target:: positions, vertical_normal
    type(scalar_field) :: psi, DistanceToTop, exact
    type(mesh_type) :: psi_mesh
    type(mesh_type), pointer :: x_mesh
    type(element_type) :: psi_shape
    type(quadrature_type) :: quad
    integer :: quad_degree=4, unit
    integer, parameter:: DIM=3
    real :: eps0
    character(len=999) filename, exact_sol_filename
    character(len=PYTHON_FUNC_LEN) :: func, buffer
    logical :: vl_as, vl, no_vl, sor, vl_as_wsor
    logical :: file_exists

    call set_global_debug_level(3)

    ewrite(1,*) 'test_pressure_solve'

    call pressure_solve_options(filename, eps0, &
         & exact_sol_filename, vl_as, vl_as_wsor, vl, no_vl, sor)

    ewrite(2,*) 'Using mesh files:',trim(filename)
    ewrite(2,*) 'epsilon =', eps0
    ewrite(2,*) 'using exact solution file', trim(exact_sol_filename)
    ewrite(2,*) vl_as, vl_as_wsor, vl, no_vl, sor

    if(vl_as.or.(vl.or.(no_vl.or.(vl_as_wsor.or.sor)))) then

       positions=read_mesh_files(trim(filename), &
            quad_degree=QUAD_DEGREE, format="gmsh")

       x_mesh => positions%mesh

       call insert(state, positions, name="Coordinate")
       call insert(state, positions%mesh, "Coordinate_mesh")
       
       call allocate(vertical_normal, mesh_dim(x_mesh), x_mesh, &
         field_type=FIELD_TYPE_CONSTANT, name="GravityDirection")
       call zero(vertical_normal)
       call set(vertical_normal, mesh_dim(x_mesh), -1.0)
       call insert(state, vertical_normal, name="GravityDirection")
       call deallocate(vertical_normal)

       ! Shape functions for psi
       assert(dim==mesh_dim(positions))
       quad=make_quadrature(vertices=dim+1, dim=dim, degree=quad_degree)

       psi_shape=make_element_shape(vertices=dim+1, dim=dim, degree=1, quad=quad)
       psi_mesh=make_mesh(model=positions%mesh, shape=psi_shape)

       call insert(state, psi_mesh, "Psi_Mesh")
       call allocate(psi, psi_mesh, "Psi")
       call allocate(DistanceToTop, psi_mesh, "DistanceToTop")
       call add_boundary_condition(DistanceToTop,"top","surface",(/1/) )

       call insert(state, psi, "Psi")
       call insert(state, DistanceToTop,"DistanceToTop")

       inquire(file=trim(exact_sol_filename),exist=file_exists)  
       if (.not.file_exists) FLAbort('Couldnt find exact_sol_filename file')
       unit=free_unit() 
       open(unit, file=trim(exact_sol_filename), action="read",&
            & status="old")
       read(unit, '(a)', end=43) func
       ! Read all the lines of the file and put in newlines between them.
       do
          read(unit, '(a)', end=43) buffer
          func=trim(func)//achar(10)//trim(buffer)
       end do
43     func=trim(func)//achar(10)
       close(unit)
       
       call allocate(exact,psi%mesh,name='Exact')
       call set_from_python_function(exact, trim(func), positions, 0.0)
       
       positions%val(1,:) = positions%val(1,:)/eps0
       positions%val(2,:) = positions%val(2,:)/eps0
       
       call insert(state,exact,'Exact')
       
       call run_model(state,vl_as,vl_as_wsor,vl,no_vl,sor)
    end if

  end subroutine test_pressure_solve

  subroutine run_model(state,vl_as,vl_as_wsor,vl,no_vl,sor)
    use global_parameters, only: PYTHON_FUNC_LEN
    use unittest_tools
    use solvers
    use boundary_conditions
    use fields
    use state_module
    use elements
    use sparse_tools_petsc
    use sparsity_patterns
    use boundary_conditions
    use free_surface_module
    use FLDebug
    use multigrid
    use spud
    implicit none
    type(state_type), intent(inout) :: state
    logical, intent(in) :: vl_as, vl_as_wsor, vl, no_vl, sor

    type(vector_field), pointer :: positions
    type(scalar_field), pointer :: psi, exact
    type(scalar_field) :: error
    type(scalar_field), pointer :: topdis
    type(mesh_type) :: top_surface_mesh
    type(petsc_csr_matrix) :: vprolongator
    integer, dimension(:), pointer :: top_surface_node_list => null()
    integer, dimension(:), pointer :: top_surface_element_list => null()
    !character(len=OPTION_PATH_LEN) solver_option_path
    integer :: stat

    ! We form and solve the equation A*psi=rhs
    type(csr_sparsity) :: A_sparsity
    type(csr_matrix) :: A
    type(scalar_field) :: RHS

    ! Extract the required fields from state.
    psi=>extract_scalar_field(state, "Psi")
    positions=>extract_vector_field(state, "Coordinate")
    exact=>extract_scalar_field(state, "Exact")

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

    call allocate(rhs, psi%mesh, "RHS")
    call zero(rhs)
    
    call get_laplacian(A,positions,psi)
    !call set_reference_node(A, top_surface_node_list(1), rhs, 0.0)
    call set(A, top_surface_node_list(1),top_surface_node_list(1), INFINITY)

    call allocate(error, psi%mesh, "Error")
    call zero(error)

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
         atol=1.0e-100, rtol=1.0e-20, max_its=10000, &
         start_from_zero=.true.)
         
    call add_option(trim(psi%option_path)//'/solver/diagnostics/monitors/true_error', stat=stat)

    if(vl_as) then
       ewrite(1,*) 'with vertical lumping and internal smoother'

       call petsc_solve_monitor_exact(exact, error_filename='with_vl_and_is.dat')
       call petsc_solve(psi, A, rhs, prolongators=(/ vprolongator /), &
            & surface_node_list=top_surface_node_list, &
            & internal_smoothing_option=INTERNAL_SMOOTHING_SEPARATE_SOR)
    end if

    if(vl_as_wsor) then       
       ewrite(1,*) 'with vertical lumping and internal smoother and wrapped &
            &sor'
       call petsc_solve_monitor_exact(exact, error_filename='with_vl_and_is_wrap_sor.dat')
       call petsc_solve(psi, A, rhs, prolongators=(/ vprolongator /), &
            & surface_node_list=top_surface_node_list, &
            & internal_smoothing_option=INTERNAL_SMOOTHING_WRAP_SOR)
    end if

    if(vl) then
       ewrite(1,*) 'with vertical lumping, no additive smoother'

       call petsc_solve_monitor_exact(exact, error_filename='with_vl_without_is.dat')
       call petsc_solve(psi, A, rhs, prolongators=(/ vprolongator /))
    end if

    if(no_vl) then
       ewrite(1,*) 'without vertical lumping'

       call petsc_solve_monitor_exact(exact, error_filename='without_vl.dat')
       call petsc_solve(psi, A, rhs)
    end if

    if(sor) then
       ewrite(1,*) 'Using SOR'

       call petsc_solve_monitor_exact(exact, error_filename='sor.dat')
       call petsc_solve(psi, A, rhs)
    end if

  end subroutine run_model

  subroutine get_laplacian(A,positions,psi)
    use sparse_tools
    use fields
    implicit none
    type(csr_matrix), intent(inout) :: A
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: psi
    !
    integer :: ele

    do ele = 1, element_count(psi)
       call assemble_laplacian_element_contribution(&
            &A, positions, psi, ele)
    end do

  end subroutine get_laplacian

  subroutine assemble_laplacian_element_contribution(A, positions, psi, ele)
    use unittest_tools
    use solvers
    use fields
    use state_module
    use elements
    use sparse_tools
    implicit none
    type(csr_matrix), intent(inout) :: A
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: psi
    integer, intent(in) :: ele

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

    ele_psi=>ele_nodes(psi, ele)
    shape_psi=>ele_shape(psi, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_psi, dshape=dshape_psi,&
         & detwei=detwei)

    ! Local assembly:
    psi_mat=dshape_dot_dshape(dshape_psi, dshape_psi, detwei)

    ! Global assembly:
    call addto(A, ele_psi, ele_psi, psi_mat)

  end subroutine assemble_laplacian_element_contribution

  subroutine pressure_solve_options(filename, eps0, &
       & exact_sol_filename, vl_as, vl_as_wsor, vl, no_vl, sor)
    use Fldebug
    use petsc_tools
#ifdef HAVE_PETSC_MODULES
  use petsc
#endif
  implicit none
#include "petsc_legacy.h"

    character(len=*), intent(out):: filename, exact_sol_filename
    logical, intent(out) :: vl_as, vl, no_vl, sor, vl_as_wsor
    real, intent(out) :: eps0

#if PETSC_VERSION_MINOR>=2
    PetscBool:: flag
#else
    PetscTruth:: flag
#endif
    PetscErrorCode :: ierr
    PetscReal :: number_in=0.0

    call PetscOptionsGetString(&
         &PETSC_NULL_CHARACTER, '-filename', filename, flag, ierr)
    if (.not. flag) then
       call usage()
    end if

    call PetscOptionsGetReal(PETSC_NULL_CHARACTER, '-epsilon', number_in, flag, ierr)
    if(.not. flag) then
       call usage()
    end if
    eps0 = number_in

    call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-exact_solution', &
         & exact_sol_filename, flag, ierr) 
    if (.not. flag) then 
       call usage()
    end if

    call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-vl_as', vl_as, ierr)

    call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-vl_as_wsor', vl_as_wsor, ierr)

    call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-vl', vl, ierr)

    call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-no_vl', no_vl, ierr)

    call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-sor', sor, ierr)

  end subroutine pressure_solve_options

  subroutine usage()
    use FLDebug    
    ewrite(0,*) 'Usage: test_pressure_solve -filename <filename> &
         &-exact_solution <exact_solution_python_filename> -epsilon &
         &<epsilon> [options ...]'
    ewrite(0,*) 'Options:'
    ewrite(0,*) '-vl_as' 
    ewrite(0,*) '       Performs a solve using vertical lumping with additive &
         &smoother'
    ewrite(0,*) '-vl_as_wsor' 
    ewrite(0,*) '       Performs a solve using vertical lumping with additive &
         &smoother and wrapped sor'
    ewrite(0,*) '-vl'
    ewrite(0,*) '       Performs a solve using vertical lumping without additi&
         &ve smoother'
    ewrite(0,*) '-no_vl'
    ewrite(0,*) '       Performs a solve using regular mg'
    stop
  end subroutine usage


