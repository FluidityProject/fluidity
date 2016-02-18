subroutine test_vertical_integration()
use quadrature
use fields
use state_module
use mesh_files
use boundary_conditions
use vertical_extrapolation_module
use fields_calculations
use vtk_interfaces
use unittest_tools
implicit none
  
  character, parameter:: NEWLINE_CHAR=achar(10)
  ! vertically constant:
  character(len=*), parameter:: PYTHON_FUNCTION1= &
     "def val(X, t):"//NEWLINE_CHAR// &
     "  import math"//NEWLINE_CHAR// &
     "  return math.cos(math.pi*X[0])*math.sin(math.pi*X[1])"
  ! function with vertical derivative:
  character(len=*), parameter:: PYTHON_FUNCTION2= &
     "def val(X, t):"//NEWLINE_CHAR// &
     "  import math"//NEWLINE_CHAR// &
     "  return math.cos(math.pi*X[0])*math.sin(math.pi*X[1])*math.cos(math.pi*X[2])"
  ! its vertical derivative (note gravity is in -X[2] direction):
  character(len=*), parameter:: PYTHON_FUNCTION_DERIVATIVE= &
     "def val(X, t):"//NEWLINE_CHAR// &
     "  import math"//NEWLINE_CHAR// &
     "  return math.pi*math.cos(math.pi*X[0])*math.sin(math.pi*X[1])*math.sin(math.pi*X[2])"
  ! the mesh files should have the top surface marked with boundary id 1
  integer, dimension(1), parameter:: TOP_BOUNDARY_IDS=(/ 1 /)
  integer, parameter:: QUAD_DEGREE=4, NVERTICES=4, DIM=3
  integer, parameter:: POLY_DEGREE=2 ! degree of fields we're integrating
  real, dimension(1:DIM), parameter:: DOWN=(/ 0.0, 0.0, -1.0 /)
  
  logical fail
  real l2error
  
  call test_vertical_integration_from_file("data/cube_prismatic", &
    PYTHON_FUNCTION1, l2error=l2error)
  fail= l2error>1e-10
  call report_test("[test_vertical_integration_prismatic]", fail, .false., &
    "Too large error in vertical integration on prismatic mesh.")
  
  call test_vertical_integration_from_file("data/cube_unstructured", &
    PYTHON_FUNCTION1, l2error=l2error)
  fail= l2error>2e-3
  call report_test("[test_vertical_integration_unstructured]", fail, .false., &
    "Too large error in vertical integration on unstructured mesh.")

  call test_vertical_integration_from_file("data/cube_prismatic", &
    PYTHON_FUNCTION2, python_function_derivative=PYTHON_FUNCTION_DERIVATIVE, &
    l2error=l2error)
  fail= l2error>1e-3
  call report_test("[test_vertical_integration_prismatic_gradient]", fail, .false., &
    "Too large error in vertical integration on prismatic mesh.")
    
  call test_vertical_integration_from_file("data/cube_unstructured", &
    PYTHON_FUNCTION2, python_function_derivative=PYTHON_FUNCTION_DERIVATIVE, &
    l2error=l2error)
  fail= l2error>5e-3
  call report_test("[test_vertical_integration_unstructured_gradient]", fail, .false., &
    "Too large error in vertical integration on unstructured mesh.")
  
contains

subroutine test_vertical_integration_from_file(mesh_file, &
  python_function, python_function_derivative, &
  l2error)
character(len=*), intent(in):: mesh_file
character(len=*), intent(in):: python_function
character(len=*), optional, intent(in):: python_function_derivative
real, intent(out):: l2error

  type(state_type) state
  type(vector_field), target:: positions, vertical_normal, bc_positions
  type(scalar_field) to_field, from_field, error_field, rhs
  type(mesh_type), pointer:: x_mesh, surface_mesh
  type(mesh_type) dg_quad_mesh
  type(element_type) quad_shape
  type(quadrature_type) quad
  integer, dimension(:), pointer:: surface_element_list
    
  ! vertical integration() needs a "Coordinate" and  "GravityDirection" field
  positions = read_mesh_files(mesh_file, quad_degree=QUAD_DEGREE, format="gmsh")
  x_mesh => positions%mesh
  
  call allocate(vertical_normal, mesh_dim(x_mesh), x_mesh, &
    name="GravityDirection", &
    field_type=FIELD_TYPE_CONSTANT)
  call set(vertical_normal, DOWN)

  call insert(state, positions, name="Coordinate")
  call insert(state, vertical_normal, name="GravityDirection")
  
  ! we're going to integrate downwards over an quadratic DG field:
  quad=x_mesh%shape%quadrature
  quad_shape=make_element_shape(NVERTICES, DIM, POLY_DEGREE, quad)
  dg_quad_mesh=make_mesh(x_mesh, shape=quad_shape, continuity=-1, &
    name="DGMesh")
  call deallocate(quad_shape)
  call allocate(to_field, dg_quad_mesh, name="ToField")
  
  ! set the function on the top of the mesh:
  call add_boundary_condition(to_field, "Top", type="Dummy", &
    boundary_ids=TOP_BOUNDARY_IDS )
  call get_boundary_condition(to_field, "Top", &
    surface_mesh=surface_mesh, &      
    surface_element_list=surface_element_list)
  call allocate(from_field, surface_mesh, name="FromField")
  call allocate(bc_positions, DIM, surface_mesh, name="BCPositions")
  call remap_field_to_surface(positions, bc_positions, surface_element_list)
  call set_from_python_function(from_field, python_function, bc_positions, &
      time=0.0)
      
  if (present(python_function_derivative)) then
    call allocate(rhs, dg_quad_mesh)
    call set_from_python_function(rhs, python_function_derivative, &
      positions, time=0.0)
    
    call vertical_integration(from_field, to_field, positions, &
      vertical_normal, surface_element_list, rhs=rhs)
  else
    call vertical_integration(from_field, to_field, positions, &
      vertical_normal, surface_element_list)
  end if
    
  call allocate(error_field, dg_quad_mesh, name="ErrorField")
  ! first set it to the reference solution using the same python function
  call set_from_python_function(error_field, python_function, positions, &
      time=0.0)
  ! then subtract the found solution:
  call addto(error_field, to_field, scale=-1.0)
  
  !call vtk_write_fields(mesh_file, 0, &
  !  positions, dg_quad_mesh, sfields=(/ to_field, error_field /))

  l2error=norm2(error_field, positions)
  print *, l2error
  
  call deallocate(state)
  call deallocate(vertical_normal)
  call deallocate(to_field)
  call deallocate(from_field)
  call deallocate(error_field)
  call deallocate(bc_positions)
  if (present(python_function_derivative)) then
    call deallocate(rhs)
  end if
  
  call deallocate(positions)
  call deallocate(dg_quad_mesh)
      
end subroutine test_vertical_integration_from_file
  
end subroutine test_vertical_integration
