#include "fdebug.h"

subroutine vertical_integration(target_basename, target_basename_len, &
  & integrated_filename, integrated_filename_len, &
  & integrated_fieldname, integrated_fieldname_len, &
  & output_basename, output_basename_len, &
  & top, bottom, sizing, field_b_degree)
  
  use fields
  use fldebug
  use global_parameters, only : current_debug_level, OPTION_PATH_LEN
  use hadapt_extrude
  use intersection_finder_module
  use linked_lists
  use read_triangle
  use reference_counting
  use solvers
  use sparse_tools
  use sparsity_patterns
  use spud
  use state_module
  use supermesh_assembly
  use supermesh_construction
  use tetrahedron_intersection_module
  use vtk_interfaces
  use write_triangle
  
  implicit none
  
  integer, intent(in) :: target_basename_len
  integer, intent(in) :: integrated_filename_len
  integer, intent(in) :: integrated_fieldname_len
  integer, intent(in) :: output_basename_len
  
  character(len = target_basename_len), intent(in) :: target_basename
  character(len = integrated_filename_len), intent(in) :: integrated_filename
  character(len = integrated_fieldname_len), intent(in) :: integrated_fieldname
  character(len = output_basename_len), intent(in) :: output_basename
  
  real, intent(in) :: top
  real, intent(in) :: bottom
  real, intent(in) :: sizing
  integer, intent(in) :: field_b_degree
  
#ifdef DUMP_INTERSECTIONS
  integer :: ndumped_intersections = 0
#endif
  character(len = OPTION_PATH_LEN) :: base_path, mesh_path
  integer :: dim, i, intersecting_element, j, k, l, nintersecting_elements, quad_degree, stat, ntests
  logical :: allocated_field_a
  type(csr_matrix) :: matrix
  type(csr_sparsity) :: sparsity
  type(element_type) :: field_b_shape, field_b_shape_ext
  type(element_type), pointer :: field_c_shape, positions_b_surf_shape
  type(mesh_type) :: field_b_mesh, positions_b_mesh_ele
  type(plane_type), dimension(4) :: planes_b
  type(scalar_field) :: field_c, field_b, rhs
  type(scalar_field), pointer :: field_a
  type(state_type) :: state_a
  type(tet_type) :: tet_a, tet_b
  type(vector_field) :: positions_b_surf_ele, positions_b_ext, positions_b_surf, &
    & positions_c
  type(vector_field), pointer :: positions_a
  
  ewrite(1, *) "In vertical_integration"
  
  quad_degree = max(field_b_degree * 2, 1)
  
  ! Step 1: Read in the data
  
  positions_b_surf = read_triangle_files(target_basename, quad_degree = quad_degree)
  dim = positions_b_surf%dim + 1

  assert(ele_count(positions_b_surf) > 0)
  positions_b_surf_shape => ele_shape(positions_b_surf, 1)
  field_b_shape = make_element_shape(positions_b_surf_shape%loc, positions_b_surf_shape%dim, field_b_degree, quad = positions_b_surf_shape%quadrature)
  if(field_b_degree == 0) then
    field_b_mesh = make_mesh(positions_b_surf%mesh,field_b_shape, name = "VerticalIntegralMesh", continuity = -1)
  else
    field_b_mesh = make_mesh(positions_b_surf%mesh,field_b_shape, name = "VerticalIntegralMesh")
  end if
  call deallocate(field_b_shape)
  call allocate(field_b, field_b_mesh, trim(integrated_fieldname) // "VerticalIntegral")
  call deallocate(field_b_mesh)
  
  call vtk_read_state(integrated_filename, state_a, quad_degree = quad_degree)
  positions_a => extract_vector_field(state_a, "Coordinate")
  if(positions_a%dim /= positions_b_surf%dim + 1) then
    FLAbort("Integrated mesh must have dimension one higher than the target mesh")
  end if
  field_a => extract_scalar_field(state_a, integrated_fieldname, allocated = allocated_field_a)
  
  ! Step 2: Set up the options
  
  mesh_path = "/geometry/mesh::ExtrudedCoordinateMesh"
  base_path = trim(mesh_path) // "/from_mesh/extrude/regions[0]"
  call set_option("/geometry/quadrature/degree", quad_degree, stat)
  assert(stat == SPUD_NEW_KEY_WARNING)
  call set_option(trim(base_path) // "/bottom_depth/constant", top - bottom, stat)
  assert(stat == SPUD_NEW_KEY_WARNING)
  call set_option(trim(base_path) // "/sizing_function/constant", sizing, stat)
  assert(stat == SPUD_NEW_KEY_WARNING)
  call set_solver_options(field_b, ksptype = "cg", pctype = "sor", rtol = 1.0e-10, max_its = 3000)
  if(current_debug_level >= 2) then
    ewrite(2, *) "Options tree:"
    call print_options()
  end if

  ! Step 3: Extrude each element in the target mesh, supermesh and assemble
 
#ifdef DUMP_EXTRUSION
  ! Using positions_b_ext as working memory here
  call extrude(positions_b_surf, mesh_path, positions_b_ext)
  ! and apply the offset
  positions_b_ext%val(dim)%ptr = positions_b_ext%val(dim)%ptr + top
  call write_triangle_files("extruded_vertical_integration_mesh", positions_b_ext)
  call deallocate(positions_b_ext)
#endif
  
  sparsity = make_sparsity(field_b%mesh, field_b%mesh, name = "Sparsity")
  call allocate(matrix, sparsity, name = "Matrix")
  call deallocate(sparsity)
  call allocate(rhs, field_b%mesh, "Rhs")
  
  call rtree_intersection_finder_set_input(positions_a)
  call intersector_set_dimension(dim)
  call intersector_set_exactness(.false.)
    
  call zero(matrix)
  call zero(rhs)
  do i = 1, ele_count(positions_b_surf)
    positions_b_surf_shape => ele_shape(positions_b_surf, i)
    field_b_shape = ele_shape(field_b, i)
    
    ! Assemble the LHS matrix for this element
    call assemble_mass_ele(i, positions_b_surf, field_b, matrix)
  
    ! Generate a small surface mesh for just this element
    call allocate(positions_b_mesh_ele, ele_loc(positions_b_surf, i), 1, positions_b_surf_shape, "ElementIntegralMesh")
    call set_ele_nodes(positions_b_mesh_ele, 1, (/(j, j = 1, ele_loc(positions_b_surf, i))/))
    call allocate(positions_b_surf_ele, dim - 1, positions_b_mesh_ele, "ElementIntegralCoordinate")
    call set(positions_b_surf_ele, ele_nodes(positions_b_surf_ele, 1), ele_val(positions_b_surf, i))
    
    ! Extrude the elemental surface mesh
    call extrude(positions_b_surf_ele, mesh_path, positions_b_ext)
    ! and apply the offset
    positions_b_ext%val(dim)%ptr = positions_b_ext%val(dim)%ptr + top
    
    ! Loop over the extruded elemental mesh
    do j = 1, ele_count(positions_b_ext)
      ! Find intersections with the mesh we are vertically integrating
      call rtree_intersection_finder_find(positions_b_ext, j)
      call rtree_intersection_finder_query_output(nintersecting_elements)
      
      if(nintersecting_elements > 0 .and. dim == 3) then
        tet_b%v = ele_val(positions_b_ext, j)
        planes_b = get_planes(tet_b)
      end if
      
      do k = 1, nintersecting_elements
        call rtree_intersection_finder_get_output(intersecting_element, k)
        
        ! Supermesh each intersection
        if(dim == 3 .and. (intersector_exactness .eqv. .false.)) then
          tet_A%v = ele_val(positions_a, intersecting_element)
          call intersect_tets(tet_a, planes_b, ele_shape(positions_a, intersecting_element), stat = stat, output = positions_c)
          if(stat == 1) then
            ! No intersection to integrate
            cycle
          end if
        else
          positions_c = intersect_elements(positions_b_ext, j, ele_val(positions_a, intersecting_element), ele_shape(positions_a, intersecting_element))
        end if
        
        if(ele_count(positions_c) == 0) then
          ! No intersection to integrate
          call deallocate(positions_c)
          cycle
        end if
        if(volume(positions_c) < epsilon(0.0)) then
          ! Negligable intersection to integrate
          call deallocate(positions_c)
          cycle
        end if
        
        
#ifdef DUMP_INTERSECTIONS
        call vtk_write_fields("vertical_integration_intersections", index = ndumped_intersections, position = positions_c, model = positions_c%mesh)
        ndumped_intersections = ndumped_intersections + 1
#endif
        
        ! Project the donor field onto the supermesh
        call allocate(field_c, positions_c%mesh, "IntersectionIntegratedField")
        do l = 1, node_count(field_c)
          call set(field_c, l, &
              & dot_product( &
                & ele_val(field_a, intersecting_element), &
                & eval_shape(ele_shape(field_a, intersecting_element), local_coords(positions_a, intersecting_element, node_val(positions_c, l))) &
              & ) &
            & )
        end do
        
        do l = 1, ele_count(field_c)
          ! Project the extruded target test function onto the supermesh
          field_c_shape => ele_shape(field_c, l)
          field_b_shape_ext = extruded_shape_function(i, l, positions_b_surf, positions_c, field_b_shape, field_c_shape, &
            & form_dn = .false.)
                     
          ! Assemble the RHS for this intersection for this element
          call assemble_rhs_ele(l, i, field_b_shape_ext, positions_c, field_c, rhs)
          
          call deallocate(field_b_shape_ext)
        end do
        
        call deallocate(positions_c)
        call deallocate(field_c)
      end do
    end do

    call deallocate(positions_b_mesh_ele)
    call deallocate(positions_b_surf_ele)
    call deallocate(positions_b_ext) 
  end do

  call rtree_intersection_finder_reset(ntests)
  if(dim == 3) call finalise_tet_intersector()

  ewrite_minmax(rhs%val)
  
  ! Step 4: Solve
  
  call zero(field_b)
  call petsc_solve(field_b, matrix, rhs)
  ewrite_minmax(field_b%val)
  call deallocate(matrix)
  call deallocate(rhs)
  
  ! Step 5: Write output
  
  call vtk_write_fields(output_basename, position = positions_b_surf, model = field_b%mesh, sfields = (/field_b/))
  
  ! Step 6: Cleanup
  
  call deallocate(positions_b_surf)
  call deallocate(field_b)
  call deallocate(state_a)
  if(allocated_field_a) deallocate(field_a)
  
  call print_references(0)
  
  ewrite(1, *) "Exiting vertical_integration"
  
contains
  
  function volume(positions)
    type(vector_field), intent(in) :: positions
    
    integer :: i
    real :: volume
    
    volume = 0.0
    do i = 1, ele_count(positions)
      volume = volume + element_volume(positions, i)
    end do
  
  end function volume
  
  subroutine assemble_mass_ele(ele, positions, field, matrix)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: field
    type(csr_matrix), intent(inout) :: matrix
    
    integer, dimension(:), pointer :: element_nodes
    real, dimension(ele_ngi(positions, ele)) :: detwei
    type(element_type), pointer :: shape
    
    shape => ele_shape(field, ele)
    
    call transform_to_physical(positions, ele, detwei = detwei)
      
    element_nodes => ele_nodes(rhs, ele)
    call addto(matrix, element_nodes, element_nodes, shape_shape(shape, shape, detwei))
  
  end subroutine assemble_mass_ele
  
  subroutine assemble_rhs_ele(ele, ele_out, test_function, positions, field, rhs)
    integer, intent(in) :: ele
    integer, intent(in) :: ele_out
    type(element_type), intent(in) :: test_function
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: field
    type(scalar_field), intent(inout) :: rhs
    
    real, dimension(ele_ngi(positions, ele)) :: detwei
    
    call transform_to_physical(positions, ele, detwei = detwei)
    
    call addto(rhs, ele_nodes(rhs, ele_out), shape_rhs(test_function, detwei * ele_val_at_quad(field, ele)))
  
  end subroutine assemble_rhs_ele

end subroutine vertical_integration
