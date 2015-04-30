#include "confdefs.h"
#include "fdebug.h"

subroutine unifiedmesh(filename1, filename1_len, &
                        filename2, filename2_len, &
                        output, output_len)

  use mpi_interfaces
  use fldebug
  use mesh_files
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

  positionsA = read_mesh_files(trim(filename1), quad_degree=1, format="gmsh")
  positionsB = read_mesh_files(trim(filename2), quad_degree=1, format="gmsh")

  dim = positionsA%dim

  allocate(map_BA(ele_count(positionsB)))

  supermesh_quad = make_quadrature(vertices=dim+1, dim=dim, degree=5)
  supermesh_shape = make_element_shape(vertices=dim+1, dim=dim, degree=1, quad=supermesh_quad)
  allocate(tri_detwei(supermesh_shape%ngi))

  call allocate(accum_mesh, 0, 0, supermesh_shape, "AccumulatedMesh")
  call allocate(accum_positions, dim, accum_mesh, "AccumulatedPositions")

  map_BA = intersection_finder(positionsB, positionsA)
  call intersector_set_dimension(dim)

  ! inputs: positionsB, map_BA, positionsA, supermesh_shape
  ! output: the supermesh!
  call recursive_supermesh(positionsA, positionsB, map_BA, supermesh_shape, 1, ele_count(positionsB), accum_positions)
!  do ele_B=1,ele_count(positionsB)
!    llnode => map_BA(ele_B)%firstnode
!    do while(associated(llnode))
!      ele_A = llnode%value
!      intersection = intersect_elements(positionsA, ele_A, ele_val(positionsB, ele_B), supermesh_shape)
!      call unify_meshes(accum_positions, intersection, accum_positions_tmp)
!      call deallocate(accum_positions)
!      accum_positions = accum_positions_tmp

!      llnode => llnode%next
!      call deallocate(intersection)
!    end do
!  end do

!   call allocate(accum_positions_tmp, ! need a way to make this continuous!
!   call remap_field(accum_positions, accum_positions_tmp)
!   call write_triangle_files(trim(output), accum_positions_tmp)
!   call deallocate(accum_positions_tmp)

  call write_mesh_files(trim(output), format="gmsh", positions=accum_positions)
  call vtk_write_fields(trim(output), 0, accum_positions, accum_positions%mesh)

  contains
  recursive subroutine recursive_supermesh(positionsA, positionsB, map_BA, supermesh_shape, start_ele, end_ele, supermesh)
    type(vector_field), intent(in) :: positionsA, positionsB
    type(ilist), dimension(:), intent(in) :: map_BA
    type(element_type), intent(in) :: supermesh_shape
    integer, intent(in) :: start_ele, end_ele
    type(vector_field), intent(out) :: supermesh

    integer, parameter :: blocksize = 4
    integer :: i

    type(vector_field) :: supermesh_tmp
    type(mesh_type) :: supermesh_mesh
    type(vector_field) :: supermesh_accum
    type(inode), pointer :: llnode

    integer :: ele_A, ele_B
    integer :: new_start, new_end, step

    call allocate(supermesh_mesh, 0, 0, supermesh_shape, "AccumulatedMesh")
    call allocate(supermesh, positionsA%dim, supermesh_mesh, "AccumulatedPositions")
    call deallocate(supermesh_mesh)

    if ((end_ele - start_ele) <= blocksize) then
      do ele_B=start_ele,end_ele
        llnode => map_BA(ele_B)%firstnode
        do while(associated(llnode))
          ele_A = llnode%value
          supermesh_tmp = intersect_elements(positionsA, ele_A, ele_val(positionsB, ele_B), supermesh_shape)
          call unify_meshes_quadratic(supermesh, supermesh_tmp, supermesh_accum)
          call deallocate(supermesh)
          call deallocate(supermesh_tmp)
          supermesh = supermesh_accum
          llnode => llnode%next
        end do
      end do
    else
      do i=1,blocksize
        step = (end_ele - start_ele)/blocksize
        new_start = start_ele + (i-1)*step
        if (i /= blocksize) then
          new_end   = start_ele + (i)*step - 1
        else
          new_end = end_ele ! step might not divide exactly
        end if
        write(0,*) "Calling recursive_supermesh with element range ", new_start, new_end
        call recursive_supermesh(positionsA, positionsB, map_BA, supermesh_shape, new_start, new_end, supermesh_tmp)
        call unify_meshes_quadratic(supermesh, supermesh_tmp, supermesh_accum)
        call deallocate(supermesh)
        call deallocate(supermesh_tmp)
        supermesh = supermesh_accum
      end do
    endif
  end subroutine recursive_supermesh

end subroutine unifiedmesh
