!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineeringp
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

subroutine supermesh_difference(vtu1_filename_, vtu1_filename_len, vtu2_filename_, &
            &vtu2_filename_len, output_filename_, output_filename_len) bind(c)
  
  use fields
  use fldebug
  use interpolation_module
  use intersection_finder_module
  use linked_lists
  use reference_counting
  use state_module
  use supermesh_construction
  use unify_meshes_module
  use vtk_interfaces
  use iso_c_binding
  
  implicit none
  
  integer(kind=c_size_t), value  :: vtu1_filename_len, vtu2_filename_len, output_filename_len
  character(kind=c_char, len=1) :: vtu1_filename_(*), vtu2_filename_(*), output_filename_(*)
  
  character(len = vtu1_filename_len) :: vtu1_filename
  character(len = vtu2_filename_len) :: vtu2_filename
  character(len = output_filename_len) :: output_filename
  integer :: dim, ele_1, ele_2, ele_3, i, j, nintersections, nintersections_max
  integer, dimension(:), allocatable :: ele_map_13, ele_map_23, node_map_13, &
    & node_map_23
  integer, dimension(:, :), allocatable :: intersection_parents
  logical :: p0
  type(element_type), pointer :: shape
  type(ilist), dimension(:), allocatable :: ele_map_21
  type(inode), pointer :: node
  type(mesh_type) :: pwc_mesh_3
  type(mesh_type), pointer :: mesh_3
  type(scalar_field) :: s_field_3
  type(scalar_field), pointer :: s_field, s_field_13, s_field_23
  type(state_type) :: state_1, state_2, state_3
  type(state_type), dimension(2) :: state_13, state_1_split, state_23, &
    & state_2_split, state_3_split
  type(tensor_field) :: t_field_3
  type(tensor_field), pointer :: t_field, t_field_13, t_field_23
  type(vector_field), pointer :: positions_1, positions_2, v_field, &
    & v_field_13, v_field_23
  type(vector_field) :: intersection, v_field_3
  type(vector_field), target :: positions_3
  type(vector_field), dimension(:), allocatable :: intersections
  
  ewrite(1, *) "In supermesh_difference"

  do i=1, vtu1_filename_len
    vtu1_filename(i:i)=vtu1_filename_(i)
  end do
  do i=1, vtu2_filename_len
    vtu2_filename(i:i)=vtu2_filename_(i)
  end do
  do i=1, output_filename_len
    output_filename(i:i)=output_filename_(i)
  end do

  
  call vtk_read_state(vtu1_filename, state_1)
  call vtk_read_state(vtu2_filename, state_2)
  
  p0 = has_mesh(state_1, "P0Mesh")
  
  positions_1 => extract_vector_field(state_1, "Coordinate")
  positions_2 => extract_vector_field(state_2, "Coordinate")
  
  dim = positions_1%dim
  if(positions_2%dim /= dim) then
    FLExit("Input vtu dimensions do not match")
  end if
  call intersector_set_dimension(dim)
  call intersector_set_exactness(.false.)
  
  allocate(ele_map_21(ele_count(positions_1)))
  ! Use the rtree to avoid continuity assumptions
  ele_map_21 = rtree_intersection_finder(positions_1, positions_2)
  
  nintersections_max = sum(ele_map_21%length)
  ewrite(2, "(a,i0)") "Maximum number of intersections: ", nintersections_max
  allocate(intersections(nintersections_max))
  allocate(intersection_parents(nintersections_max, 2))
  nintersections = 0
  
  do ele_1 = 1, ele_count(positions_1)
    node => ele_map_21(ele_1)%firstnode
    do while(associated(node))
      ele_2 = node%value
      
      ! TODO: Integrate the tet intersector
      intersection = intersect_elements(positions_1, ele_1, ele_val(positions_2, ele_2), ele_shape(positions_1, ele_1))
      if(ele_count(intersection) > 0) then
        nintersections = nintersections + 1
        intersections(nintersections) = intersection
        intersection_parents(nintersections, :) = (/ele_1, ele_2/)
      else
        call deallocate(intersection)
      end if
      
      node => node%next
    end do
  end do  
  call deallocate(ele_map_21)
  deallocate(ele_map_21)
  
  ewrite(2, "(a,i0)") "Number of intersections: ", nintersections
  
  positions_3 = unify_meshes(intersections(:nintersections))
  positions_3%name = "Coordinate"
  allocate(node_map_13(node_count(positions_3)))
  allocate(node_map_23(node_count(positions_3)))
  if(p0) then
    allocate(ele_map_13(ele_count(positions_3)))
    allocate(ele_map_23(ele_count(positions_3)))
  end if
  ele_3 = 0
  do i = 1, nintersections
    do j = 1, ele_count(intersections(i))
      ele_3 = ele_3 + 1
      node_map_13(ele_nodes(positions_3, ele_3)) = intersection_parents(i, 1)
      node_map_23(ele_nodes(positions_3, ele_3)) = intersection_parents(i, 2)
      if(p0) then
        ele_map_13(ele_3) = intersection_parents(i, 1)
        ele_map_23(ele_3) = intersection_parents(i, 2)
      end if
    end do
    call deallocate(intersections(i))
  end do
  deallocate(intersections)
  deallocate(intersection_parents)
  
  mesh_3 => positions_3%mesh
  pwc_mesh_3 = piecewise_constant_mesh(mesh_3, "PiecewiseConstantMesh")
    
  call insert(state_3, positions_3, "Coordinate")
  call insert(state_1_split(1), positions_1, "Coordinate")
  call insert(state_2_split(1), positions_2, "Coordinate")
  call insert(state_13(1), positions_3, "Coordinate")
  call insert(state_23(1), positions_3, "Coordinate")
  call insert(state_3_split(1), positions_3, "Coordinate")
  call insert(state_1_split(1), positions_1%mesh, "Mesh")
  call insert(state_13(1), mesh_3, "Mesh")
  call insert(state_2_split(1), positions_2%mesh, "Mesh")
  call insert(state_23(1), mesh_3, "Mesh")
  call insert(state_3, mesh_3, "Mesh")
  call insert(state_3_split(1), mesh_3, "Mesh")
  if(p0) then
    call insert(state_1_split(2), positions_1, "Coordinate")
    call insert(state_2_split(2), positions_2, "Coordinate")
    call insert(state_13(2), positions_3, "Coordinate")
    call insert(state_23(2), positions_3, "Coordinate")
    call insert(state_3_split(2), positions_3, "Coordinate")
    call insert(state_1_split(2), extract_mesh(state_1, "P0Mesh"), "Mesh")
    call insert(state_2_split(2), extract_mesh(state_2, "P0Mesh"), "Mesh")
    call insert(state_13(2), pwc_mesh_3, "Mesh")
    call insert(state_23(2), pwc_mesh_3, "Mesh")
    call insert(state_3_split(2), pwc_mesh_3, "Mesh")
  end if
  call deallocate(positions_3)
    
  do i = 1, scalar_field_count(state_1)
    s_field => extract_scalar_field(state_1, i)
    shape => ele_shape(s_field, 1)
    
    if(shape%degree == 0) then    
      assert(p0)      
      call insert(state_1_split(2), s_field, s_field%name)
      call insert(state_2_split(2), extract_scalar_field(state_2, s_field%name), s_field%name)
      
      call allocate(s_field_3, pwc_mesh_3, s_field%name)
      call insert(state_13(2), s_field_3, s_field_3%name)
    else
      call insert(state_1_split(1), s_field, s_field%name)
      call insert(state_2_split(1), extract_scalar_field(state_2, s_field%name), s_field%name)
      
      call allocate(s_field_3, mesh_3, s_field%name)
      call insert(state_13(1), s_field_3, s_field_3%name)
    end if
    call deallocate(s_field_3)
    
    if(shape%degree == 0) then
      call allocate(s_field_3, pwc_mesh_3, s_field%name)
      call insert(state_23(2), s_field_3, s_field_3%name)
    else
      call allocate(s_field_3, mesh_3, s_field%name)
      call insert(state_23(1), s_field_3, s_field_3%name)
    end if
    call deallocate(s_field_3)
    
    if(shape%degree == 0) then
      call allocate(s_field_3, pwc_mesh_3, s_field%name)
      call insert(state_3_split(2), s_field_3, s_field_3%name)
    else
      call allocate(s_field_3, mesh_3, s_field%name)
      call insert(state_3_split(1), s_field_3, s_field_3%name)
    end if
    call insert(state_3, s_field_3, s_field_3%name)
    call deallocate(s_field_3)    
  end do    
  do i = 1, vector_field_count(state_1)
    v_field => extract_vector_field(state_1, i)
    if(v_field%name == "Coordinate") cycle
    shape => ele_shape(v_field, 1)
    
    if(shape%degree == 0) then    
      assert(p0)      
      call insert(state_1_split(2), v_field, v_field%name)
      call insert(state_2_split(2), extract_scalar_field(state_2, v_field%name), v_field%name)
      
      call allocate(v_field_3, dim, pwc_mesh_3, v_field%name)
      call insert(state_13(2), v_field_3, v_field_3%name)
    else
      call insert(state_1_split(1), v_field, v_field%name)
      call insert(state_2_split(1), extract_scalar_field(state_2, v_field%name), v_field%name)
      
      call allocate(v_field_3, dim, mesh_3, v_field%name)
      call insert(state_13(1), v_field_3, v_field_3%name)
    end if
    call deallocate(v_field_3)
    
    if(shape%degree == 0) then
      call allocate(v_field_3, dim, pwc_mesh_3, v_field%name)
      call insert(state_23(2), v_field_3, v_field_3%name)
    else
      call allocate(v_field_3, dim, mesh_3, v_field%name)
      call insert(state_23(1), v_field_3, v_field_3%name)
    end if
    call deallocate(v_field_3)
    
    if(shape%degree == 0) then
      call allocate(v_field_3, dim, pwc_mesh_3, v_field%name)
      call insert(state_3_split(2), v_field_3, v_field_3%name)
    else
      call allocate(v_field_3, dim, mesh_3, v_field%name)
      call insert(state_3_split(1), v_field_3, v_field_3%name)
    end if
    call insert(state_3, v_field_3, v_field_3%name)
    call deallocate(v_field_3)    
  end do
  do i = 1, tensor_field_count(state_1)
    t_field => extract_tensor_field(state_1, i)
    shape => ele_shape(t_field, 1)
    
    if(shape%degree == 0) then    
      assert(p0)      
      call insert(state_1_split(2), t_field, t_field%name)
      call insert(state_2_split(2), extract_scalar_field(state_2, t_field%name), t_field%name)
      
      call allocate(t_field_3, pwc_mesh_3, t_field%name)
      call insert(state_13(2), t_field_3, t_field_3%name)
    else
      call insert(state_1_split(1), t_field, t_field%name)
      call insert(state_2_split(1), extract_scalar_field(state_2, t_field%name), t_field%name)
      
      call allocate(t_field_3, mesh_3, t_field%name)
      call insert(state_13(1), t_field_3, t_field_3%name)
    end if
    call deallocate(t_field_3)
    
    if(shape%degree == 0) then
      call allocate(t_field_3, pwc_mesh_3, t_field%name)
      call insert(state_23(2), t_field_3, t_field_3%name)
    else
      call allocate(t_field_3, mesh_3, t_field%name)
      call insert(state_23(1), t_field_3, t_field_3%name)
    end if
    call deallocate(t_field_3)
    
    if(shape%degree == 0) then
      call allocate(t_field_3, pwc_mesh_3, t_field%name)
      call insert(state_3_split(2), t_field_3, t_field_3%name)
    else
      call allocate(t_field_3, mesh_3, t_field%name)
      call insert(state_3_split(1), t_field_3, t_field_3%name)
    end if
    call insert(state_3, t_field_3, t_field_3%name)
    call deallocate(t_field_3)    
  end do  
  call deallocate(pwc_mesh_3)
  
  call linear_interpolation(state_1_split(1), state_13(1), map = node_map_13)
  deallocate(node_map_13)  
  call linear_interpolation(state_2_split(1), state_23(1), map = node_map_23)
  deallocate(node_map_23)   
  if(p0) then 
    call linear_interpolation(state_1_split(2), state_13(2), map = ele_map_13)
    deallocate(ele_map_13)
    call linear_interpolation(state_2_split(2), state_23(2), map = ele_map_23)
    deallocate(ele_map_23)
  end if
  call deallocate(state_1_split)
  call deallocate(state_1)
  call deallocate(state_2_split)
  call deallocate(state_2)
  
  do i = 1, scalar_field_count(state_13(1))
    s_field_13 => extract_scalar_field(state_13(1), i)
    s_field_23 => extract_scalar_field(state_23(1), s_field_13%name)
    s_field_3 = extract_scalar_field(state_3_split(1), s_field_13%name)
    s_field_3%val = s_field_13%val - s_field_23%val
  end do
  do i = 1, vector_field_count(state_13(1))
    v_field_13 => extract_vector_field(state_13(1), i)
    if(v_field_13%name == "Coordinate") cycle
    v_field_23 => extract_vector_field(state_23(1), v_field_13%name)
    v_field_3 = extract_vector_field(state_3_split(1), v_field_13%name)
    do j = 1, dim
      v_field_3%val(i,:) = v_field_13%val(i,:) - v_field_23%val(i,:)
    end do
  end do
  do i = 1, tensor_field_count(state_13(1))
    t_field_13 => extract_tensor_field(state_13(1), i)
    t_field_23 => extract_tensor_field(state_23(1), t_field_13%name)
    t_field_3 = extract_tensor_field(state_3_split(1), t_field_13%name)    
    t_field_3%val = t_field_13%val - t_field_23%val
  end do
  if(p0) then
    do i = 1, scalar_field_count(state_13(2))
      s_field_13 => extract_scalar_field(state_13(2), i)
      s_field_23 => extract_scalar_field(state_23(2), s_field_13%name)
      s_field_3 = extract_scalar_field(state_3_split(2), s_field_13%name)
      s_field_3%val = s_field_13%val - s_field_23%val
    end do
    do i = 1, vector_field_count(state_13(2))
      v_field_13 => extract_vector_field(state_13(2), i)
      if(v_field_13%name == "Coordinate") cycle
      v_field_23 => extract_vector_field(state_23(2), v_field_13%name)
      v_field_3 = extract_vector_field(state_3_split(2), v_field_13%name)
      do j = 1, dim
        v_field_3%val(i,:) = v_field_13%val(i,:) - v_field_23%val(i,:)
      end do
    end do
    do i = 1, tensor_field_count(state_13(2))
      t_field_13 => extract_tensor_field(state_13(2), i)
      t_field_23 => extract_tensor_field(state_23(2), t_field_13%name)
      t_field_3 = extract_tensor_field(state_3_split(2), t_field_13%name)    
      t_field_3%val = t_field_13%val - t_field_23%val
    end do
  end if
  call deallocate(state_13)
  call deallocate(state_23)
  call deallocate(state_3_split)
  
  call vtk_write_state(output_filename, state = (/state_3/))
  call deallocate(state_3)
  
  call print_references(0)
  
  ewrite(1, *) "Exiting supermesh_difference"
  
end subroutine supermesh_difference
