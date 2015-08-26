#include "fdebug.h"

subroutine vertical_integration(target_basename_, target_basename_len, &
  & integrated_filename_, integrated_filename_len, &
  & output_basename_, output_basename_len, &
  & top, bottom, sizing, field_b_continuity, field_b_degree) bind(c)

  use fields
  use fldebug
  use global_parameters, only : current_debug_level, OPTION_PATH_LEN
  use hadapt_extrude
  use halos
  use intersection_finder_module
  use linked_lists
  use read_gmsh
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
  use mesh_files
  use iso_c_binding

  implicit none

  integer(kind=c_size_t), value :: target_basename_len
  integer(kind=c_size_t), value :: integrated_filename_len
  integer(kind=c_size_t), value :: output_basename_len
  real(kind=c_double), value :: top, bottom, sizing
  integer(kind=c_int), value :: field_b_continuity
  integer(kind=c_int), value :: field_b_degree
  character(kind=c_char, len=1) :: target_basename_(*)
  character(kind=c_char, len=1) :: integrated_filename_(*), output_basename_(*)


  character(len = target_basename_len) :: target_basename
  character(len = integrated_filename_len) :: integrated_filename
  character(len = output_basename_len) :: output_basename
  character(len = *), parameter :: solver_path = "/temporary/solver/path"
  character(len = OPTION_PATH_LEN) :: base_path, mesh_path
  integer :: dim, ele_a, ele_b, ele_b_surf, i, index, j, k, nele_as, stat, ntests
  integer, parameter :: quad_degree = 4
  real, dimension(:, :), allocatable :: lshape
  type(csr_matrix) :: matrix
  type(csr_sparsity) :: sparsity
  type(element_type) :: field_b_shape, field_b_shape_ext
  type(element_type), pointer :: field_c_shape, positions_b_surf_shape
  type(mesh_type), pointer :: field_a_mesh, field_c_mesh
  type(mesh_type) :: field_b_mesh
  type(plane_type), dimension(4) :: planes_b
  type(scalar_field), dimension(:), allocatable :: field_a, field_b, field_c, rhs
  type(state_type) :: state_a, state_b
  type(tet_type) :: tet_a, tet_b
  type(vector_field) :: positions_b_ext, positions_b_surf, vfield_b
  type(vector_field), pointer :: positions_a, vfield_a
  type(vector_field), target :: positions_c

  ewrite(-1, *) "In vertical_integration"

  do i=1, target_basename_len
    target_basename(i:i)=target_basename_(i)
  end do
  do i=1, integrated_filename_len
    integrated_filename(i:i)=integrated_filename_(i)
  end do
  do i=1, output_basename_len
    output_basename(i:i)=output_basename_(i)
  end do
  ! Step 1: Read in the data
  write(*,*) target_basename, output_basename, integrated_filename

  positions_b_surf = read_gmsh_file(target_basename, quad_degree = quad_degree)
  dim = positions_b_surf%dim + 1

  call vtk_read_state(integrated_filename, state_a, quad_degree = quad_degree)
  positions_a => extract_vector_field(state_a, "Coordinate")
  field_a_mesh => extract_mesh(state_a, "Mesh")
  if(positions_a%dim /= positions_b_surf%dim + 1) then
    ewrite(-1, *) "Integrated mesh dimension: ", positions_a%dim
    ewrite(-1, *) "Target mesh dimension: ", positions_b_surf%dim
    FLExit("Integrated mesh must have dimension one higher than the target mesh")
  end if

  assert(ele_count(positions_b_surf) > 0)
  positions_b_surf_shape => ele_shape(positions_b_surf, 1)
  field_b_shape = make_element_shape(positions_b_surf_shape%numbering%vertices, positions_b_surf_shape%dim, field_b_degree, quad = positions_b_surf_shape%quadrature)
  if(field_b_degree == 0) then
    field_b_mesh = make_mesh(positions_b_surf%mesh, field_b_shape, name = "VerticalIntegralMesh", continuity = field_b_continuity)
  else
    field_b_mesh = make_mesh(positions_b_surf%mesh, field_b_shape, name = "VerticalIntegralMesh", continuity = field_b_continuity)
  end if
  call deallocate(field_b_shape)

  call insert(state_b, field_b_mesh, field_b_mesh%name)
  call insert(state_b, positions_b_surf, "Coordinate")
  assert(has_vector_field(state_a, "Coordinate"))
  allocate(field_a(scalar_field_count(state_a) + dim * (vector_field_count(state_a) - 1)))
  allocate(field_b(size(field_a)))
  do i = 1, scalar_field_count(state_a)
    field_a(i) = extract_scalar_field(state_a, i)
    if(.not. field_a(i)%mesh == field_a_mesh) then
      if(field_a(i)%mesh%shape%degree == 0) then
        FLExit("vertical_integration does not support degree zero donor fields")
      else
        ewrite(-1, *) "Donor field mesh: " // trim(field_a(i)%mesh%name)
        ewrite(-1, *) "With degree: ", field_a(i)%mesh%shape%degree
        FLAbort("Unexpected donor field mesh")
      end if
    end if
    call allocate(field_b(i), field_b_mesh, field_a(i)%name)
    call zero(field_b(i))
    call insert(state_b, field_b(i), field_b(i)%name)
    call deallocate(field_b(i))
    field_b(i) = extract_scalar_field(state_b, i)
  end do
  index = scalar_field_count(state_a)
  do i = 1, vector_field_count(state_a)
    vfield_a => extract_vector_field(state_a, i)
    if(.not. vfield_a%mesh == field_a_mesh) then
      if(vfield_a%mesh%shape%degree == 0) then
        FLExit("vertical_integration does not support degree zero donor fields")
      else
        ewrite(-1, *) "Donor field mesh: " // trim(vfield_a%mesh%name)
        ewrite(-1, *) "With degree: ", vfield_a%mesh%shape%degree
        FLAbort("Unexpected donor field mesh")
      end if
    end if
    if(vfield_a%name == "Coordinate") cycle
    call allocate(vfield_b, dim, field_b_mesh, vfield_a%name)
    call zero(vfield_b)
    call insert(state_b, vfield_b, vfield_b%name)
    call deallocate(vfield_b)
    vfield_b = extract_vector_field(state_b, vfield_a%name)
    do j = 1, dim
      field_a(index + j) = extract_scalar_field(vfield_a, j)
      field_b(index + j) = extract_scalar_field(vfield_b, j)
    end do
    index = index + dim
  end do
  ! Currently no way to handle dim dimensional tensors on a (dim - 1)
  ! dimensional mesh
  assert(index == size(field_a))

  ! Step 2: Set up the options

  mesh_path = "/geometry/mesh::ExtrudedCoordinateMesh"
  base_path = trim(mesh_path) // "/from_mesh/extrude/regions[0]"
  call set_option("/geometry/quadrature/degree", quad_degree, stat)
  assert(stat == SPUD_NEW_KEY_WARNING)
  call set_option(trim(base_path) // "/bottom_depth/constant", top - bottom, stat)
  assert(stat == SPUD_NEW_KEY_WARNING)
  call set_option(trim(base_path) // "/sizing_function/constant", sizing, stat)
  assert(stat == SPUD_NEW_KEY_WARNING)
  call set_solver_options(solver_path, ksptype = "cg", pctype = "sor", rtol = 1.0e-10, max_its = 3000)
  if(current_debug_level >= 2) then
    ewrite(2, *) "Options tree:"
    call print_options()
  end if

  ! Step 3: Extrude the target mesh, supermesh and assemble

  ! Extrude the surface mesh
  call extrude(positions_b_surf, mesh_path, positions_b_ext)
  ! and apply the offset
  positions_b_ext%val(dim,:) = positions_b_ext%val(dim,:) + top
#ifdef DUMP_EXTRUSION
  call write_mesh_files("extruded_vertical_integration_mesh", format="gmsh", positions=positions_b_ext)
#endif

  sparsity = make_sparsity(field_b_mesh, field_b_mesh, name = "Sparsity")
  call allocate(matrix, sparsity, name = "Matrix")
  call deallocate(sparsity)
  call zero(matrix)
  allocate(rhs(size(field_b)))
  do i = 1, size(field_b)
    call allocate(rhs(i), field_b_mesh, "RHS")
    call zero(rhs(i))
  end do

  call rtree_intersection_finder_set_input(positions_a)
  call intersector_set_dimension(dim)
  call intersector_set_exactness(.false.)

  ! Assemble the mass matrix
  do ele_b_surf = 1, ele_count(positions_b_surf)
    call assemble_mass_ele(ele_b_surf, positions_b_surf, field_b_mesh, matrix)
  end do

  ! Assemble the RHS
  allocate(field_c(size(rhs)))
  do ele_b = 1, ele_count(positions_b_ext)
    ewrite(2, "(a,i0,a,i0)") "Processing element ", ele_b, " of ", ele_count(positions_b_ext)

    assert(associated(positions_b_ext%mesh%element_columns))
    ele_b_surf = positions_b_ext%mesh%element_columns(ele_b)

    ! Find intersections with the mesh we are vertically integrating
    call rtree_intersection_finder_find(positions_b_ext, ele_b)
    call rtree_intersection_finder_query_output(nele_as)

    if(nele_as > 0 .and. dim == 3 .and. (intersector_exactness .eqv. .false.)) then
      tet_b%v = ele_val(positions_b_ext, ele_b)
      planes_b = get_planes(tet_b)
    end if

    do i = 1, nele_as
      call rtree_intersection_finder_get_output(ele_a, i)

      ! Supermesh each intersection
      if(dim == 3 .and. (intersector_exactness .eqv. .false.)) then
        tet_A%v = ele_val(positions_a, ele_a)
        call intersect_tets(tet_a, planes_b, ele_shape(positions_a, ele_a), stat = stat, output = positions_c)
        if(stat == 1) then
          ! No intersection to integrate
          cycle
        end if
      else
        positions_c = intersect_elements(positions_b_ext, ele_b, ele_val(positions_a, ele_a), ele_shape(positions_a, ele_a))
      end if

      if(ele_count(positions_c) == 0) then
        ! No intersection to integrate
        call deallocate(positions_c)
        cycle
      end if
      if(dim /= 3 .or. (intersector_exactness .eqv. .true.)) then  ! The stat argument to intersect_tets checks this
        if(volume(positions_c) < epsilon(0.0)) then
          ! Negligable intersection to integrate
          call deallocate(positions_c)
          cycle
        end if
      end if

      ! Assume here that the donor Coordinate and fields are on the same mesh
      assert(positions_a%mesh == field_a_mesh)
      field_c_mesh => positions_c%mesh
      allocate(lshape(ele_loc(field_c_mesh, 1), node_count(field_c_mesh)))
      do j = 1, node_count(field_c_mesh)
        lshape(:, j) = eval_shape(ele_shape(field_a_mesh, 1), local_coords(positions_a, ele_a, node_val(positions_c, j)))
      end do

      do j = 1, size(rhs)
        call allocate(field_c(j), field_c_mesh, "IntersectionIntegratedField")
        ! Project the donor field onto the supermesh
        do k = 1, node_count(field_c(j))
          call set(field_c(j), k, &
              & dot_product( &
                & ele_val(field_a(j), ele_a), &
                  & lshape(:, k) &
              & ) &
            & )
        end do
      end do
      deallocate(lshape)

      do j = 1, ele_count(positions_c)
        ! Project the extruded target test function onto the supermesh
        field_c_shape => ele_shape(field_c(1), j)
        field_b_shape_ext = extruded_shape_function(ele_b_surf, j, positions_b_surf, positions_c, field_b_mesh%shape, field_c_shape, &
          & form_dn = .false.)

        do k = 1, size(rhs)
          ! Assemble the RHS for this intersection for this element
          call assemble_rhs_ele(j, ele_b_surf, field_b_shape_ext, positions_c, field_c(k), rhs(k))
        end do

        call deallocate(field_b_shape_ext)
      end do

      do j = 1, size(rhs)
        call deallocate(field_c(j))
      end do
      call deallocate(positions_c)
    end do
  end do
  deallocate(field_c)

  call deallocate(positions_b_ext)
  call rtree_intersection_finder_reset(ntests)
  if(dim == 3) call finalise_tet_intersector()

  ! Step 4: Solve

  call petsc_solve(field_b, matrix, rhs, option_path = solver_path)
  call deallocate(matrix)
  do i = 1, size(field_a)
    call deallocate(rhs(i))
  end do

  ! Step 5: Write output

  call vtk_write_state(output_basename, model = field_b_mesh%name, state = (/state_b/))

  ! Step 6: Cleanup

  call deallocate(positions_b_surf)
  call deallocate(field_b_mesh)
  call deallocate(state_a)
  call deallocate(state_b)
  deallocate(field_a)
  deallocate(field_b)

  call print_references(0)

  ewrite(1, *) "Exiting vertical_integration"

contains

  function volume(positions)
    type(vector_field), intent(in) :: positions
    
    real :: volume

    integer :: i

    volume = 0.0
    do i = 1, ele_count(positions)
      volume = volume + element_volume(positions, i)
    end do

  end function volume

  subroutine assemble_mass_ele(ele, positions, mesh, matrix)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(mesh_type), intent(in) :: mesh
    type(csr_matrix), intent(inout) :: matrix

    integer, dimension(:), pointer :: element_nodes
    real, dimension(ele_ngi(positions, ele)) :: detwei
    type(element_type), pointer :: shape

    shape => ele_shape(mesh, ele)

    call transform_to_physical(positions, ele, detwei = detwei)

    element_nodes => ele_nodes(mesh, ele)
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
