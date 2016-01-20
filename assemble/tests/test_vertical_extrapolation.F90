subroutine test_vertical_extrapolation()
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
  ! the mesh files should have the top surface marked with boundary id 1
  integer, dimension(1), parameter:: TOP_BOUNDARY_IDS=(/ 2 /)
  integer, parameter:: QUAD_DEGREE=4, NVERTICES=4, DIM=3
  integer, parameter:: POLY_DEGREE=2 ! degree of fields we're integrating
  
  logical fail
  real l2error
  
  call test_vertical_extrapolation_from_file("data/cube_prismatic", &
    PYTHON_FUNCTION1, l2error)
  fail= l2error>1e-10
  call report_test("[test_vertical_extrapolation_prismatic]", fail, .false., &
    "Too large error in vertical extrapolation on prismatic mesh.")
  
  call test_vertical_extrapolation_from_file("data/cube_unstructured", &
    PYTHON_FUNCTION1, l2error)
  fail= l2error>1e-3  
  call report_test("[test_vertical_extrapolation_unstructured]", fail, .false., &
    "Too large error in vertical extrapolation on unstructured mesh.")

  
contains

subroutine test_vertical_extrapolation_from_file(mesh_file, &
  python_function, l2error)
character(len=*), intent(in):: mesh_file
character(len=*), intent(in):: python_function
real, intent(out):: l2error

  type(vector_field), target:: positions, bc_positions, vertical_normal
  type(scalar_field) to_field, from_field, error_field, from_surface_field
  type(mesh_type), pointer:: x_mesh, surface_mesh
  type(mesh_type) dg_quad_mesh
  type(element_type) quad_shape
  type(quadrature_type) quad
  integer, dimension(:), pointer:: surface_element_list
  integer:: i, sele
    
  positions=read_mesh_files(mesh_file, quad_degree=QUAD_DEGREE, format="gmsh")
  x_mesh => positions%mesh
  
  ! we're going to extrapolate downwards onto a quadratic DG field:
  quad=x_mesh%shape%quadrature
  quad_shape=make_element_shape(NVERTICES, DIM, POLY_DEGREE, quad)
  dg_quad_mesh=make_mesh(x_mesh, shape=quad_shape, continuity=-1, &
    name="DGMesh")
  call allocate(to_field, dg_quad_mesh, name="ToField")
  
  ! set the function on the top of the mesh:
  call add_boundary_condition(to_field, "Top", type="Dummy", &
    boundary_ids=TOP_BOUNDARY_IDS )
  call get_boundary_condition(to_field, "Top", &
    surface_mesh=surface_mesh, &      
    surface_element_list=surface_element_list)
  call allocate(from_surface_field, surface_mesh, name="FromSurfaceField")
  call allocate(bc_positions, DIM, surface_mesh, name="BCPositions")
  call remap_field_to_surface(positions, bc_positions, surface_element_list)
  call set_from_python_function(from_surface_field, python_function, bc_positions, &
      time=0.0)
  ! map values of from_surface_field on faces of from_field
  call allocate(from_field, dg_quad_mesh, name="FromField")
  ! make sure other values are not used:
  call set(from_field, huge(0.0))
  do i=1, size(surface_element_list)
    sele=surface_element_list(i)
    call set(from_field, face_global_nodes(dg_quad_mesh, i), &
      ele_val(from_surface_field, i))
  end do
      
  call allocate(vertical_normal, positions%dim, x_mesh, &
    field_type=FIELD_TYPE_CONSTANT, name="VerticalNormal")
  call set(vertical_normal, (/ 0., 0., -1. /) )
      
  call VerticalExtrapolation(from_field, to_field, positions, &
    vertical_normal, surface_element_list)
    
  call allocate(error_field, dg_quad_mesh, name="ErrorField")
  ! first set it to the reference solution using the same python function
  call set_from_python_function(error_field, python_function, positions, &
      time=0.0)
  ! then subtract the found solution:
  call addto(error_field, to_field, scale=-1.0)
  
!  call vtk_write_fields(mesh_file, 0, &
!    positions, dg_quad_mesh, sfields=(/ from_field, to_field, error_field /))

  l2error=norm2(error_field, positions)
  print *, l2error
  
  call deallocate(to_field)
  call deallocate(from_field)
  call deallocate(from_surface_field)
  call deallocate(error_field)
  call deallocate(bc_positions)
  call deallocate(vertical_normal)
      
end subroutine test_vertical_extrapolation_from_file
  
end subroutine test_vertical_extrapolation
