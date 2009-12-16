subroutine test_vertical_integration
  use vertical_integration
  use read_triangle
  use Fields_manipulation
  use boundary_conditions
  use vtk_interfaces
  implicit none
  
  !field to integrate, integrated field
  type(scalar_field) :: field, vfield
  !positions field
  type(vector_field) :: positions, s_positions
  type(mesh_type), pointer :: s_mesh
  integer, dimension(:), pointer :: s_e_list

  positions=read_triangle_files("vert_test",quad_degree=4)

  call set_from_python_function(field, func="def val(X,t):\n   return&
       & X[0]*X[1]", position=positions, time=0.0)

  call add_boundary_condition(field, name="TopSurface", &
       type="neumann", boundary_ids=(/1/))
  call get_boundary_condition(field,"TopSurface",surface_mesh=s_mesh)

  call allocate(vfield,s_mesh,name='Integral')

  call get_vertical_integral(field, positions, vfield)

  call allocate(s_positions, 3, s_mesh)

  call remap_vector_field_to_surface(&
       positions, s_positions, s_e_list)

  call vtk_write_fields("vertical_out", &
       & 0, S_Positions, S_Positions%mesh, (/vfield/))

end subroutine test_vertical_integration
