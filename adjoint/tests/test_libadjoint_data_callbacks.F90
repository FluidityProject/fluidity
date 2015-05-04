subroutine test_libadjoint_data_callbacks

  use unittest_tools

#include "confdefs.h"
#ifndef HAVE_ADJOINT
    call report_test("[libadjoint_datacallbacks]", .false., .false., "")
#else

  use mesh_files
  use fields
  use libadjoint
  use libadjoint_data_callbacks

  implicit none

  type(vector_field) :: positions
  type(mesh_type) :: mesh
  type(element_type) :: shape
  type(scalar_field) :: fieldA, fieldB, fieldC
  type(adj_vector) :: adj_fieldA, adj_fieldB
  integer :: i
  logical :: fail = .false.

  positions = read_mesh_files("data/pslgA", quad_degree=4, format="gmsh")
  shape = make_element_shape(vertices = ele_loc(positions, 1), dim  = positions%dim, degree = 3, quad = positions%mesh%shape%quadrature)
  mesh = make_mesh(positions%mesh, shape, name="Mesh")
  call allocate(fieldA, mesh, "FieldA")
  call allocate(fieldB, mesh, "FieldB")

  do i=1,size(fieldB%val)
    call random_number(fieldB%val(i))
  end do
  call set(fieldA, 0.0)  

  adj_fieldA = field_to_adj_vector(fieldA)
  adj_fieldB = field_to_adj_vector(fieldB)
  ! A = 2*B 
  call femtools_vec_axpy_proc(adj_fieldA, 2.0, adj_fieldB)
  ! Undo this operation manually:
  call field_from_adj_vector(adj_fieldA, fieldC)
  call addto(fieldC, fieldB, -2.0)

  if (any(fieldC%val/=0.0)) then
    fail=.true.
  end if
  call deallocate(fieldA)
  call deallocate(fieldB)

  call femtools_vec_destroy_proc(adj_fieldA)  
  call femtools_vec_destroy_proc(adj_fieldB)

  call deallocate(mesh)
  call deallocate(shape)

  call report_test("[libadjoint_datacallbacks]", fail, .false., "|C| = 0")

#endif

end subroutine test_libadjoint_data_callbacks


