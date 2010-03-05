#include "fdebug.h"

subroutine test_extrude_azimuthal
  use fldebug
  use hadapt_extrude
  use fields
  use spud
  use unittest_tools
  use vtk_interfaces
  use read_triangle
  use write_triangle
  use sparse_tools
  use global_parameters
  implicit none

  integer, parameter:: quad_degree = 1
  type(ilist) :: seeds
  type(vector_field) :: v_mesh, out_mesh
  type(scalar_field) :: paint
  
  v_mesh = read_triangle_files("data/triangle.2", quad_degree = quad_degree)
  v_mesh%val(1)%ptr = v_mesh%val(1)%ptr + 1.0
  
  call extrude_azimuthal(v_mesh, out_mesh, ndivisions = 3)
  
  call write_triangle_files("data/test_extrude_azimuthal_out", out_mesh)
  
  call report_test("[nodes]", node_count(out_mesh) /= 4 * 3, .false., "Incorrect number of nodes")
  call report_test("[elements]", ele_count(out_mesh) /= 3 * 3 * 3, .false., "Incorrect number of elements")
  seeds = advancing_front_intersection_finder_seeds(out_mesh)
  call paint_connectivity(out_mesh, seeds, paint)
  call vtk_write_fields("data/test_extrude_azimuthal_out.vtu", position = out_mesh, model = paint%mesh, sfields = (/paint/))
  call deallocate(paint)
  call deallocate(seeds)
  
  call deallocate(v_mesh)
  call deallocate(out_mesh)
  
  call report_test_no_references()
  
contains
  
  subroutine paint_connectivity(positions, seeds, paint) 
    type(vector_field), intent(in) :: positions
    type(ilist), intent(in) :: seeds
    type(scalar_field), intent(out) :: paint
    
    integer :: ele, i
    type(csr_sparsity), pointer :: eelist
    logical, dimension(:), allocatable :: tested
    integer, dimension(:), pointer :: neigh
    type(ilist) :: next
    type(inode), pointer :: seed
    type(mesh_type) :: pwc_mesh
    
    pwc_mesh = piecewise_constant_mesh(out_mesh%mesh, name = "PiecewiseConstantMesh")
    call allocate(paint, pwc_mesh, "Connectivity")
    call deallocate(pwc_mesh)
  
    eelist => extract_eelist(positions)
    
    allocate(tested(ele_count(positions)))
    tested = .false.
    
    seed => seeds%firstnode
    do while(associated(seed))
      ele = seed%value
      assert(ele > 0)
      assert(ele <= ele_count(positions))
      assert(.not. tested(ele))
      
      call set(paint, ele, float(seed%value))
      tested(ele) = .true.
      neigh => row_m_ptr(eelist, ele)
      do i = 1, size(neigh)
        if(neigh(i) <= 0) cycle
        if(tested(neigh(i))) cycle
      
        call insert(next, neigh(i))
      end do
      
      do while(next%length > 0)
        ele = pop(next)
        if(tested(ele)) cycle
        
        call set(paint, ele, float(seed%value))
        tested(ele) = .true.
        neigh => row_m_ptr(eelist, ele)
        do i = 1, size(neigh)
          if(neigh(i) <= 0) cycle
          if(tested(neigh(i))) cycle
          ! Should check if neigh(i) is already in the list
          
          call insert(next, neigh(i))
        end do
      end do
      seed => seed%next
    end do
    assert(all(tested))

    deallocate(tested)
    
  end subroutine paint_connectivity
  
end subroutine test_extrude_azimuthal
