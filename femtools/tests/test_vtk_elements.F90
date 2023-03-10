#include "fdebug.h"

subroutine test_vtk_elements

   !> @brief
   !! Unit test testing that we write and read VTK elements correctly.
   !! This is mostly an issue because VTK uses a different node ordering
   !! sometimes depending on the element shape

   use fields
   use fldebug
   use state_module
   use unittest_tools
   use vtk_interfaces
   use quadrature
   use elements


   implicit none

   type(quadrature_type) :: quad
   type(element_type) :: ele
   type(mesh_type) :: mesh
   type(mesh_type), pointer :: mesh_in
   integer :: i, n_nodes, n_elements = 4
   type(vector_field) :: pos
   type(scalar_field), dimension(1) :: sfields
   type(state_type) :: state

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! VTK_TRIANGLE (=5)
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   quad = make_quadrature(3, 2, 1); ele  = make_element_shape(3, 2, 1, quad)
   n_nodes = ele%loc * n_elements
   call allocate(mesh, n_nodes, n_elements, ele, "Mesh")
   mesh%ndglno = vtk_mesh2fluidity_numbering([(i, i=1, n_nodes)], ele)

   call allocate(pos, 3, mesh, "Pos")
   pos%val(1,:) = [0.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 0.0, 0.5, 0.0, 0.0, 0.5]
   pos%val(2,:) = [0.0, 0.0, 0.5, 0.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 0.0, 0.5]
   pos%val(3,:) = [0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2]

   call allocate(sfields(1), mesh, "TestField")
   sfields(1)%val = [(i, i=1, n_nodes)]
   call vtk_write_fields("vtk_triangle", 1, pos, mesh, sfields)

   call vtk_read_state("data/ref_vtk_triangle.vtu", state)
   mesh_in => extract_mesh(state, "Mesh")

   call report_test("[vtk_triangle]", .not. all(mesh%ndglno == mesh_in%ndglno),&
      .false., "DOFs are not in Fluidity ordering")

   call deallocate(quad); call deallocate(ele); call deallocate(mesh)
   call deallocate(pos); call deallocate(sfields(1)); call deallocate(state)

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! VTK_QUAD (=9)
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   quad = make_quadrature(4, 2, 1); ele  = make_element_shape(4, 2, 1, quad)
   n_nodes = ele%loc * n_elements
   call allocate(mesh, n_nodes, n_elements, ele, "Mesh")
   mesh%ndglno = vtk_mesh2fluidity_numbering([(i, i=1, n_nodes)], ele)

   call allocate(pos, 3, mesh, "Pos")
   pos%val(1,:) = [0.0, 1.0, 1.0, 0.0, 1.0, 2.0, 2.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 2.0, 2.0, 1.0]
   pos%val(2,:) = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 2.0, 2.0]
   pos%val(3,:) = [0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2]

   call allocate(sfields(1), mesh, "TestField")
   sfields(1)%val = [(i, i=1, n_nodes)]
   call vtk_write_fields("vtk_quad", 1, pos, mesh, sfields)

   call vtk_read_state("data/ref_vtk_quad.vtu", state)
   mesh_in => extract_mesh(state, "Mesh")

   call report_test("[vtk_quad]", .not. all(mesh%ndglno == mesh_in%ndglno), &
      .false., "DOFs are not in Fluidity ordering")

   call deallocate(quad); call deallocate(ele); call deallocate(mesh)
   call deallocate(pos); call deallocate(sfields(1)); call deallocate(state)


   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! VTK_TETRA (=10)
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   quad = make_quadrature(4, 3, 1); ele  = make_element_shape(4, 3, 1, quad)
   n_nodes = ele%loc * n_elements
   call allocate(mesh, n_nodes, n_elements, ele, "Mesh")
   mesh%ndglno = vtk_mesh2fluidity_numbering([(i, i=1, n_nodes)], ele)

   call allocate(pos, 3, mesh, "Pos")
   pos%val(1,:) = [0.0, 1.0, 0.5, 0.5, 1.0, 1.0, 0.5, 0.5, 1.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5]
   pos%val(2,:) = [0.0, 0.0, 0.5, 0.5, 0.0, 1.0, 0.5, 0.5, 1.0, 1.0, 0.5, 0.5, 1.0, 0.0, 0.5, 0.5]
   pos%val(3,:) = [0.0, 0.0, 0.0, 1.0, 0.1, 0.1, 0.1, 1.0, 0.3, 0.3, 0.3, 1.0, 0.2, 0.2, 0.2, 1.0]

   call allocate(sfields(1), mesh, "TestField")
   sfields(1)%val = [(i, i=1, n_nodes)]
   call vtk_write_fields("vtk_tetra", 1, pos, mesh, sfields)

   call vtk_read_state("data/ref_vtk_tetra.vtu", state)
   mesh_in => extract_mesh(state, "Mesh")

   call report_test("[vtk_tetra]", .not. all(mesh%ndglno == mesh_in%ndglno), &
      .false., "DOFs are not in Fluidity ordering")

   call deallocate(quad); call deallocate(ele); call deallocate(mesh)
   call deallocate(pos); call deallocate(sfields(1)); call deallocate(state)

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! VTK_HEXAHEDRON (=12)
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   quad = make_quadrature(8, 3, 1); ele  = make_element_shape(8, 3, 1, quad)
   n_nodes = ele%loc * n_elements
   call allocate(mesh, n_nodes, n_elements, ele, "Mesh")
   mesh%ndglno = vtk_mesh2fluidity_numbering([(i, i=1, n_nodes)], ele)

   call allocate(pos, 3, mesh, "Pos")
   pos%val(1,:) = [0.5, 0.6, 0.6, 0.5, 0.0, 1.2, 1.2, 0.0, &
      0.7, 1.9, 1.9, 0.7, 1.2, 1.3, 1.3, 1.2, &
      1.9, 2.0, 2.0, 1.9, 1.4, 2.6, 2.6, 1.4, &
      2.1, 3.3, 3.3, 2.1, 2.6, 2.7, 2.7, 2.6 ]
   pos%val(2,:) = [0.0, 0.0, 1.2, 1.2, 0.5, 0.5, 0.6, 0.6, &
      0.5, 0.5, 0.6, 0.6, 0.0, 0.0, 1.2, 1.2, &
      0.0, 0.0, 1.2, 1.2, 0.5, 0.5, 0.6, 0.6, &
      0.5, 0.5, 0.6, 0.6, 0.0, 0.0, 1.2, 1.2 ]
   pos%val(3,:) = [0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, &
      0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, &
      0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, &
      0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0 ]

   call allocate(sfields(1), mesh, "TestField")
   sfields(1)%val = [(i, i=1, n_nodes)]
   call vtk_write_fields("vtk_hexahedra", 1, pos, mesh, sfields)

   call vtk_read_state("data/ref_vtk_hexahedra.vtu", state)
   mesh_in => extract_mesh(state, "Mesh")

   call report_test("[vtk_hexahedra]", .not. all(mesh%ndglno == mesh_in%ndglno), &
      .false., "DOFs are not in Fluidity ordering")

   call deallocate(quad); call deallocate(ele); call deallocate(mesh)
   call deallocate(pos); call deallocate(sfields(1)); call deallocate(state)


end subroutine test_vtk_elements
