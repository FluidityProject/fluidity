#include "confdefs.h"
#include "fdebug.h"

subroutine unifiedmesh(filename1, filename1_len, &
                        filename2, filename2_len, &
                        output, output_len)

  use mpi_interfaces
  use fldebug
  use read_triangle
  use write_triangle
  use fields
  use linked_lists
  use intersection_finder_module
  use transform_elements
  use elements
  use supermesh_construction
  use vtk_interfaces
  use unify_meshes_module

  implicit none
  
  integer, intent(in) :: filename1_len, filename2_len, output_len
  character(len=filename1_len), intent(in) :: filename1
  character(len=filename2_len), intent(in) :: filename2
  character(len=output_len), intent(in) :: output

  type(vector_field) :: positionsA, positionsB
  type(ilist), dimension(:), allocatable :: map_BA
  real, dimension(:), allocatable :: tri_detwei
  integer :: ele_A, ele_B
  type(inode), pointer :: llnode
  type(vector_field) :: intersection
  type(element_type) :: supermesh_shape
  type(quadrature_type) :: supermesh_quad
  integer :: dim

  type(mesh_type) :: accum_mesh
  type(vector_field) :: accum_positions, accum_positions_tmp
  
  call set_global_debug_level(0)

  positionsA = read_triangle_files(trim(filename1), quad_degree=1)
  positionsB = read_triangle_files(trim(filename2), quad_degree=1)

  dim = positionsA%dim

  allocate(map_BA(ele_count(positionsB)))

  supermesh_quad = make_quadrature(vertices=dim+1, dim=dim, degree=5)
  supermesh_shape = make_element_shape(vertices=dim+1, dim=dim, degree=1, quad=supermesh_quad)
  allocate(tri_detwei(supermesh_shape%ngi))

  call allocate(accum_mesh, 0, 0, supermesh_shape, "AccumulatedMesh")
  call allocate(accum_positions, dim, accum_mesh, "AccumulatedPositions")

  map_BA = intersection_finder(positionsB, positionsA)
  call intersector_set_dimension(dim)

  do ele_B=1,ele_count(positionsB)
    llnode => map_BA(ele_B)%firstnode
    do while(associated(llnode))
      ele_A = llnode%value
      intersection = intersect_elements(positionsA, ele_A, ele_val(positionsB, ele_B), supermesh_shape)
      call unify_meshes(accum_positions, intersection, accum_positions_tmp)
      call deallocate(accum_positions)
      accum_positions = accum_positions_tmp

      llnode => llnode%next
      
      call deallocate(intersection)
    end do

  end do
  
!   call allocate(accum_positions_tmp, ! need a way to make this continuous!
!   call remap_field(accum_positions, accum_positions_tmp)
  call write_triangle_files(trim(output), accum_positions_tmp)
!   call deallocate(accum_positions_tmp)

  call vtk_write_fields(trim(output), 0, accum_positions, accum_positions%mesh)

end subroutine unifiedmesh
