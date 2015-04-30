subroutine test_laplacian
  use quadrature
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
  use unittest_tools
  use field_derivatives, only: compute_hessian
  use interpolation_module, only: linear_interpolation
  use node_boundary, only: deallocate_boundcount
  use vector_tools

  implicit none
  
  type(vector_field), target :: positions
  type(scalar_field), target :: psi, rhs, soln, load
  type(mesh_type) :: psi_mesh
  integer :: degree, quad_degree
  integer :: meshes
  integer :: solution_mesh, mesh
  type(element_type), target :: psi_shape
  type(state_type), dimension(:), pointer :: state

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
       real, dimension(:), intent(in) :: X
       real :: solution
     end function solution
  end interface
  interface 
     function loaddata(X)
       real, dimension(:), intent(in) :: X
       real :: loaddata
     end function loaddata
  end interface
  ! Arguments for handling the command line
  character(len=256) :: filename, buf

  filename = "data/laplacian_grid"
  degree = 1
  quad_degree=2*degree
  ! meshes is the number of meshes to compare computed residuals on.
  ! by convention, the higher the number, the finer it is.
  ! the problem is solved on the chosen grid and the residuals are computed
  ! both ways on the finer meshes.
  meshes = 2
  solution_mesh = 2

  allocate(state(meshes))

  do mesh=1,meshes
    write(buf, '(i0)') mesh
    positions=read_mesh_files(trim(trim(filename) // "." // trim(buf)), quad_degree=quad_degree, format="gmsh")
    call insert(state(mesh), positions, "Coordinate")
    call deallocate(positions)

  ! Shape functions for psi
    psi_shape=make_element_shape(vertices=mesh_dim(positions)+1, dim=mesh_dim(positions), degree=degree, &
    quad=positions%mesh%shape%quadrature)
    psi_mesh=make_mesh(model=positions%mesh, shape=psi_shape)
    call deallocate(psi_shape)
    call insert(state(mesh), psi_mesh, "Mesh")
    call deallocate(psi_mesh)
    call allocate(psi, psi_mesh, "ForwardSolution")
    call zero(psi)
    call insert(state(mesh), psi, "ForwardSolution")
    call deallocate(psi)
    call allocate(rhs, psi_mesh, "RightHandSide")
    call set_rhs(rhs, positions, rhs_func)
    call insert(state(mesh), rhs, "RightHandSide")
    call deallocate(rhs)
    call allocate(load, psi_mesh, "LoadData")
    call set_from_function(load, loaddata, positions)
    call insert(state(mesh), load, "LoadData")
    call deallocate(load)
    call allocate(soln, psi_mesh, "AnalyticalSolution")
    call set_from_function(soln, solution, positions)
    call insert(state(mesh), soln, "AnalyticalSolution")
    call deallocate(soln)
  end do
  
  ! Do the actual finite element calculation.
  call run_model(state(solution_mesh))
  call analyse_output(state(solution_mesh), solution)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Below: Patrick's residual estimation code
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  soln_A => extract_scalar_field(state(solution_mesh), "Forward solution")
!  positions_A => extract_vector_field(state(solution_mesh), "Coordinate")
!  do mesh=solution_mesh,meshes
!    if (mesh /= solution_mesh) then
!      soln_B => extract_scalar_field(state(mesh), "Forward solution")
!      positions_B => extract_vector_field(state(mesh), "Coordinate")
!      call interpolate_fields(positions_A, soln_A, positions_B, soln_B)  
!    end if
!
!    call discrete_residual(state(mesh))
!    call zz_residual(state(mesh))
!  end do

  if (degree<=2.0) then
     ! Output to a vtk file.
     do mesh=solution_mesh,meshes
       call vtk_write_state(trim(filename), index=mesh, state=(/state(mesh)/))
     end do
  end if

contains
  
  subroutine run_model(state)
    type(state_type), intent(in) :: state

    type(csr_matrix) :: A
    type(scalar_field), pointer :: psi, rhs

    call assemble_equations(state, A)

    psi => extract_scalar_field(state, "ForwardSolution")
    rhs => extract_scalar_field(state, "RightHandSide")

    call set_debug_level(3)
    call zero(psi)
    call set_solver_options(psi, ksptype='cg', pctype="sor", rtol=1e-7)
    call petsc_solve(psi, A, rhs)
    
    ! since A_ij=\int \nabla N_i \nabla N_j we actually wanted to solve
    ! -A psi=rhs
    call scale(psi, -1.0)
    
    call deallocate(A)
    
  end subroutine run_model

  subroutine assemble_equations(state, A)
    type(state_type), intent(in) :: state
    type(csr_matrix), intent(out) :: A

    ! We form and solve the equation A*psi=b
    type(vector_field), pointer :: positions
    type(scalar_field), pointer :: psi, rhs
    type(csr_sparsity) :: A_sparsity
    integer :: ele

    positions => extract_vector_field(state, "Coordinate")
    psi => extract_scalar_field(state, "ForwardSolution")
    rhs => extract_scalar_field(state, "RightHandSide")
    
    ! Calculate the sparsity of A based on the connectivity of psi.
    A_sparsity=make_sparsity(psi%mesh, psi%mesh, name='LaplacianSparsity')
    call allocate(A, A_sparsity)
    call zero(A)
    
    ! Assemble A element by element.
    do ele=1, element_count(psi)
       call assemble_element_contribution(A, positions, psi, ele)
    end do

    call set(A, find_zero_zero(positions, psi%mesh), find_zero_zero(positions, psi%mesh), INFINITY)
    
  end subroutine assemble_equations

  subroutine assemble_element_contribution(A, positions, psi, ele) 
    type(csr_matrix), intent(inout) :: A
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: psi
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

    ele_psi=>ele_nodes(psi, ele)
    shape_psi=>ele_shape(psi, ele)

    ! Locations of quadrature points.
    X_quad=ele_val_at_quad(positions, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_psi, dshape=dshape_psi,&
         & detwei=detwei)

    ! Matrix entry:
    call addto(A, ele_psi, ele_psi, &
         dshape_dot_dshape(dshape_psi, dshape_psi, detwei))

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
    ! Locations of quadrature points
    real, dimension(positions%dim,ele_ngi(positions,1)) :: X_quad
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,1)) :: detwei    
    ! Node numbers of rhs element.
    integer, dimension(:), pointer :: ele_rhs
    ! Shape functions.
    type(element_type), pointer :: shape_psi
    integer :: ele

    call zero(rhs)

    do ele=1,ele_count(rhs)
      ele_rhs=>ele_nodes(rhs, ele)
      shape_psi=>ele_shape(rhs, ele)
      X_quad=ele_val_at_quad(positions, ele)
      call transform_to_physical(positions, ele, detwei=detwei)
      call addto(rhs, ele_rhs, shape_rhs(shape_psi, detwei * rhs_func(X_quad)))
    end do
  end subroutine set_rhs
  
  subroutine analyse_output(state, solution)
    type(state_type), intent(inout) :: state
    interface 
       function solution(X)
         real, dimension(:), intent(in) :: X
         real :: solution
       end function solution
    end interface

    type(scalar_field), pointer :: psi, soln
    type(scalar_field) :: err
    type(vector_field), pointer :: positions
    ! Coordinate transform * quadrature weights.
    real, dimension(:), allocatable :: detwei    
    ! Shape functions.
    type(element_type), pointer :: shape_psi
    integer :: ele
    real :: integral_error
    logical :: fail

    positions => extract_vector_field(state, "Coordinate")

    psi => extract_scalar_field(state, "ForwardSolution")
    soln => extract_scalar_field(state, "AnalyticalSolution")

    ! Compute the error in the solution
    call allocate(err, psi%mesh, "Error")
    err%val = psi%val - soln%val
    call insert(state, err, "Error")
    call deallocate(err)

    allocate(detwei(ele_ngi(positions, 1)))

    integral_error = 0.0
    do ele=1,ele_count(err)
      shape_psi=>ele_shape(psi, ele)
      call transform_to_physical(positions, ele, detwei=detwei)
      integral_error = integral_error + dot_product(detwei, ele_val_at_quad(err, ele))
    end do

    deallocate(detwei)

    fail = (integral_error > 0.1)
    call report_test("[laplacian]", fail, .false., "The error in the Laplacian simulation has increased .. ")
  end subroutine analyse_output

  subroutine discrete_residual(state)
    type(state_type), intent(inout) :: state

    type(csr_matrix) :: A
    type(scalar_field) :: discrete_residual_field, soln_residual
    type(scalar_field), pointer :: psi, rhs, soln
    type(vector_field), pointer :: positions

    psi => extract_scalar_field(state, "Forward solution")
    rhs => extract_scalar_field(state, "Right-hand side")
    soln => extract_scalar_field(state, "Analytical solution")

    positions => extract_vector_field(state, "Coordinate")

    call assemble_equations(state, A)

    call allocate(discrete_residual_field, psi%mesh, "DiscreteResidual")
    call zero(discrete_residual_field)
    call allocate(soln_residual, psi%mesh, "DiscreteAnalyticalResidual")
    call zero(soln_residual)

    call mult(discrete_residual_field%val, A, psi%val)
    call mult(soln_residual%val, A, soln%val)

    discrete_residual_field%val = discrete_residual_field%val - rhs%val
    soln_residual%val = soln_residual%val - rhs%val
    call insert(state, discrete_residual_field, "DiscreteResidual")
    call deallocate(discrete_residual_field)
    call insert(state, soln_residual, "DiscreteAnalyticalResidual")
    call deallocate(soln_residual)
    call deallocate(A)
  end subroutine discrete_residual

  subroutine zz_residual(state)
    type(state_type), intent(inout) :: state

    type(tensor_field) :: hessian
    type(scalar_field) :: zz_residual_field, soln_residual

    type(scalar_field), pointer :: psi, load, soln
    type(vector_field), pointer :: positions

    psi => extract_scalar_field(state, "ForwardSolution")
    soln => extract_scalar_field(state, "AnalyticalSolution")
    load => extract_scalar_field(state, "LoadData")
    positions => extract_vector_field(state, "Coordinate")

    call allocate(hessian, psi%mesh, "Hessian")
    call allocate(zz_residual_field, psi%mesh, "ZZResidual")
    call allocate(soln_residual, psi%mesh, "ZZAnalyticalResidual")

    call compute_hessian(psi, positions, hessian)
    call trace(hessian, zz_residual_field)

    call compute_hessian(soln, positions, hessian)
    call trace(hessian, soln_residual)

    zz_residual_field%val = - 1 * zz_residual_field%val - load%val
    soln_residual%val = -1 * soln_residual%val - load%val

    call insert(state, zz_residual_field, "ZZResidual")
    call insert(state, soln_residual, "ZZAnalyticalResidual")
    call deallocate(zz_residual_field)
    call deallocate(soln_residual)

    call deallocate_boundcount
    call deallocate(hessian)
  end subroutine zz_residual

  subroutine interpolate_fields(inpositions, infield, outpositions, outfield)
    type(vector_field), intent(in) :: inpositions, outpositions
    type(scalar_field), intent(in) :: infield
    type(scalar_field), intent(inout) :: outfield

    type(state_type) :: state_in, state_out
    type(vector_field) :: inpositions_mapped, outpositions_mapped

    call nullify(state_in)
    call nullify(state_out)

    call allocate(inpositions_mapped, inpositions%dim, infield%mesh, "Coordinate")
    call remap_field(inpositions, inpositions_mapped)
    call allocate(outpositions_mapped, outpositions%dim, outfield%mesh, "Coordinate")
    call remap_field(outpositions, outpositions_mapped)

    call insert(state_in, inpositions_mapped, "Coordinate")
    call insert(state_in, infield, trim(infield%name))
    call insert(state_in, infield%mesh, "Mesh")

    call insert(state_out, outpositions_mapped, "Coordinate")
    call insert(state_out, outfield, trim(outfield%name))
    call insert(state_out, outfield%mesh, "Mesh")

    call linear_interpolation(state_in, state_out)

    call deallocate(state_in)
    call deallocate(state_out)
    call deallocate(inpositions_mapped)
    call deallocate(outpositions_mapped)
  end subroutine interpolate_fields

  function find_zero_zero(positions, mesh) result(node)
    ! find the node *closest* to (0,0)
    integer :: node
    type(vector_field), intent(in) :: positions
    type(mesh_type), intent(in) :: mesh

    type(vector_field) :: p_model
    real min_distance
    integer :: i

    call allocate(p_model, positions%dim, mesh, "Coordinate")
    call remap_field(positions, p_model)

    node = 0
    min_distance=huge(1.0)

    do i=1,node_count(p_model)
      if (norm2(node_val(p_model, i))<min_distance) then
        node = i
        min_distance=norm2(node_val(p_model, i))
      end if
    end do

    call deallocate(p_model)

  end function find_zero_zero

end subroutine test_laplacian

function rhs_func(X)
  ! Right hand side function for laplacian operator.
  !
  ! Each column of X is interpreted as a position at which RHS should be
  ! evaluated. 
  use fetools
  implicit none
  real, dimension(:,:), intent(in) :: X
  real, dimension(size(X,2)) :: rhs_func
  real, parameter :: PI=3.1415926535897931
  
  !rhs_func=-8.0*PI**2*cos(X(X_,:)*(2.0*PI))*cos(X(Y_,:)*(2.0*PI))
  rhs_func = -121 * PI * PI * cos(11 * PI * X(X_,:))
  rhs_func = rhs_func  -1*PI*sin(PI*(2*X(Y_,:) + 1)/2)
end function rhs_func

function solution(X)
  ! Analytic solution of scheme at X.
  use fetools
  implicit none
  real, dimension(:), intent(in) :: X
  real :: solution
  real, parameter :: PI=3.1415926535897931
  
  !solution=cos(X(X_,:)*(2.0*PI))*cos(X(Y_,:)*(2.0*PI))
  solution = cos(11*PI * X(X_))
  solution = solution + sin(PI*(2*X(Y_) + 1)/2)/PI
  solution = solution - 1.31830988618379
end function solution

function loaddata(X)
  ! Analytic loaddata of scheme at X.
  use fetools
  implicit none
  real, dimension(:), intent(in) :: X
  real :: loaddata
  real, parameter :: PI=3.1415926535897931
  
  !loaddata=cos(X(X_,:)*(2.0*PI))*cos(X(Y_,:)*(2.0*PI))
  loaddata = -121 * PI * PI * cos(11 * PI * X(X_))
  loaddata = loaddata  -1*PI*sin(PI*(2*X(Y_) + 1)/2)
end function loaddata
