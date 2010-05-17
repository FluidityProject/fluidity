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

  implicit none
  
  private

contains
  
  subroutine solids(state)
    implicit none

    type(state_type):: state, external_state
    character(len=field_name_len):: external_mesh_name

    type(vector_field), pointer :: positions, external_positions
    type(scalar_field) :: volume, surface
    type(mesh_type) :: pwc_mesh

    integer :: ele_B, ele_A, ele_C
    real :: vol, area
    type(tet_type) :: tet_B
    type(plane_type), dimension(6) :: planes_A
    integer :: stat, nintersections, k
    integer :: ntests
    integer :: l
    integer, dimension(:), pointer :: neighs, faces
    type(vector_field) :: sm_surface, intersection
    type(scalar_field) :: sm_colours
    integer :: face_C
    real, dimension(:), allocatable :: detwei_v, detwei
    integer :: dumpno
    real :: input_volume, output_volume
    real :: input_area, output_area
    real, dimension(3, 3) :: face_pos
    real :: total_pipe_volume, total_pipe_area
    real :: total_distance_integral

    positions=> extract_vector_field(state, "Coordinate")
    call add_faces(positions%mesh)

    call get_option("/Implicit_Solids/mesh_name", external_mesh_name)

    call vtk_read_state(trim(external_mesh_name), external_state)
    external_positions => extract_vector_field(external_state, "Coordinate")
    call add_faces(external_positions%mesh)

    pwc_mesh = piecewise_constant_mesh(positions%mesh, "PWCMesh")
    call allocate(volume, pwc_mesh, "Volume")
    call zero(volume)
    call allocate(surface, pwc_mesh, "Surface")
    call zero(surface)
    call deallocate(pwc_mesh)
    allocate(detwei_v(face_ngi(external_positions, 1)))
    allocate(detwei(ele_ngi(external_positions, 1)))

    call rtree_intersection_finder_set_input(positions)

    dumpno = 0
    total_pipe_volume = 0.
    total_pipe_area = 0.
    total_distance_integral = 0.

    do ele_B = 1, ele_count(external_positions)

       call rtree_intersection_finder_find(external_positions, ele_B)
       call rtree_intersection_finder_query_output(nintersections)
       call transform_to_physical(external_positions, ele_B, detwei=detwei)
       input_volume = sum(detwei)
       total_pipe_volume = total_pipe_volume + input_volume

       output_volume = 0.0
       tet_B%v = ele_val(external_positions, ele_B)
       neighs => ele_neigh(external_positions, ele_B)
       faces => ele_faces(external_positions, ele_B)
       input_area = 0.0
       output_area = 0.0
       do l = 1, size(neighs)
          if (neighs(l) < 0) then
             tet_B%colours(l) = ele_B
             face_pos = face_val(external_positions, faces(l))
             input_area = input_area + 0.5 * norm2(cross_product(face_pos(:, 2) - &
                  face_pos(:, 1), face_pos(:, 1) - face_pos(:, 3)))
          else
             tet_B%colours(l) = 0
          end if
       end do
       total_pipe_area = total_pipe_area + input_area

       do k = 1, nintersections
          call rtree_intersection_finder_get_output(ele_A, k)
          planes_A = get_planes(positions, ele_A)
          call intersect_tets(tet_B, planes_A, ele_shape(external_positions, ele_B), &
               stat=stat, output=intersection, surface_shape=face_shape(external_positions, 1), &
               surface_positions=sm_surface, surface_colours=sm_colours)
          if (stat == 1) then
             cycle
          end if

          !call vtk_write_fields("colours", dumpno, position=sm_surface, model=sm_surface%mesh, sfields=(/sm_colours/))
          !call vtk_write_fields("supermesh", dumpno, position=intersection, model=intersection%mesh)
          !dumpno = dumpno + 1

          ! Area and volume counting here
          vol = 0.0
          do ele_C = 1, ele_count(intersection)
             vol = vol + abs(simplex_volume(intersection, ele_C))
          end do
          output_volume = output_volume + vol
          call addto(volume, ele_A, vol)


          area = 0.0
          do face_C=1,ele_count(sm_surface)
             if (node_val(sm_colours, face_C) > 0) then
                area = area + abs(simplex_volume(sm_surface, face_C))
             end if
          end do
          output_area = output_area + area
          call addto(surface, ele_A, area)

          call deallocate(sm_surface)
          call deallocate(sm_colours)
          call deallocate(intersection)
       end do

    end do


    write(0, "(a, " // real_format(padding=1) // ")") "Total pipe volume:       ", total_pipe_volume
    write(0, "(a, " // real_format(padding=1) // ")") "Total found volume:      ", sum(volume%val)
    write(0, "(a)") ""
    write(0, "(a, " // real_format(padding=1) // ")") "Total pipe area:         ", total_pipe_area
    write(0, "(a, " // real_format(padding=1) // ")") "Total found area:        ", sum(surface%val)

    call vtk_write_fields("pipe_dreams", position=positions, model=positions%mesh, sfields=(/volume, surface/))

    call deallocate_faces(positions%mesh)
    call deallocate_faces(external_positions%mesh)
    call deallocate(external_state)
    call deallocate(volume)
    call deallocate(surface)
    call finalise_tet_intersector
    call rtree_intersection_finder_reset(ntests)
    deallocate(detwei_v)

  end subroutine solids

end module implicit_solids
