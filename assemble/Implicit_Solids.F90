!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    C.Pain@Imperial.ac.uk
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

module implicit_solids

  use state_module
  use vtk_interfaces
  use fields
  use linked_lists
  use intersection_finder_module
  use tetrahedron_intersection_module
  use unify_meshes_module
  use unittest_tools
  use global_parameters, only: FIELD_NAME_LEN
  use spud
  use read_triangle

  implicit none

  private
  public:: solids

contains

  subroutine solids(state)
    implicit none

    type(state_type),intent(inout) :: state
    character(len=field_name_len) :: external_mesh_name
    type(vector_field), pointer :: positions
    type(vector_field) :: external_positions
    type(scalar_field), pointer :: solid
    integer :: ele_A, ele_B, ele_C
    type(tet_type) :: tet_A, tet_B
    type(plane_type), dimension(4) :: planes_A
    integer :: stat, nintersections, i, j, k
    integer :: ntests
    integer, dimension(:), pointer :: ele_A_nodes
    type(vector_field) :: intersection
    real, dimension(:), allocatable :: detwei_v, detwei
    real, dimension(:,:), allocatable :: pos_A
    real :: vol, max_solid

    ewrite(3, *) "inside femdem"

    solid => extract_scalar_field(state, "SolidConcentration")
    call zero(solid)

    positions => extract_vector_field(state, "Coordinate")

    call get_option("/implicit_solids/mesh_name", external_mesh_name)
    external_positions = &
         read_triangle_files(trim(external_mesh_name), quad_degree=1)

    assert(positions%dim >= 2)
    assert(positions%dim == external_positions%dim)

    allocate(detwei_v(face_ngi(external_positions, 1)))
    allocate(detwei(ele_ngi(external_positions, 1)))

    call rtree_intersection_finder_set_input(positions)

    do ele_B = 1, ele_count(external_positions)

       call rtree_intersection_finder_find(external_positions, ele_B)
       call rtree_intersection_finder_query_output(nintersections)

       if (positions%dim == 3) then
          tet_B%v = ele_val(external_positions, ele_B)
       end if

       do j = 1, nintersections

          call rtree_intersection_finder_get_output(ele_A, j)

          if (positions%dim == 3) then

             tet_A%v = ele_val(positions, ele_A)
             planes_A = get_planes(tet_A)

             call intersect_tets(tet_B, planes_A, &
                  ele_shape(external_positions, ele_B), &
                  stat=stat, output=intersection)

          else

             call intersector_set_dimension(positions%dim)
             allocate(pos_A(positions%dim, ele_loc(positions, ele_A)))
             pos_A = ele_val(positions, ele_A)
             intersection = intersect_elements(external_positions, &
                  ele_B, pos_A, ele_shape(positions, ele_A) )
             deallocate(pos_A)

          end if

          if (stat == 1) cycle

          vol = 0.
          do ele_C = 1, ele_count(intersection)
             vol = vol + abs(simplex_volume(intersection, ele_C))
          end do

          ele_A_nodes => ele_nodes(positions, ele_A)
          do k = 1, size(ele_A_nodes)
             call addto(solid, ele_A_nodes(k), vol)
          end do

          call deallocate(intersection)

       end do

    end do

    max_solid = max( tiny(0.), maxval(solid) )
    do i = 1, node_count(solid)
       call set(solid, i, node_val(solid, i)/max_solid)
    end do

    ewrite_minmax(solid)

    call deallocate(external_positions)
    call finalise_tet_intersector
    call rtree_intersection_finder_reset(ntests)
    deallocate(detwei, detwei_v)

    ewrite(3, *) "leaving femdem"

  end subroutine solids

end module implicit_solids
