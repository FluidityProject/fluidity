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
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
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

subroutine test_trace_space

  use fields
  use fldebug
  use read_triangle
  use unittest_tools

  implicit none
  
  integer :: degree, ele, node
  integer, parameter :: min_degree = 1, max_degree = 9
  logical :: fail
  real, dimension(:, :), allocatable :: l_coords, quad_l_coords
  type(quadrature_type) :: quad
  type(element_type) :: trace_shape
  type(element_type), pointer :: base_shape
  type(mesh_type) :: trace_mesh
  type(mesh_type), pointer :: base_mesh
  type(vector_field), target :: positions
  type(vector_field) :: trace_positions
  type(scalar_field) :: D, L
  
  positions = read_triangle_files("quad_grid", quad_degree = 1)
  base_mesh => positions%mesh
  base_shape => ele_shape(base_mesh, 1)
  call report_test("[Linear quad input mesh]", &
    & ele_numbering_family(base_shape) /= FAMILY_CUBE .or. base_shape%degree /= 1 .or. base_shape%dim /= 2, .false., &
    & "Input mesh not composed of linear quads")

  do degree = min_degree, max_degree

     print "(a,i0)", "Degree = ", degree

     ! make trace mesh 
     quad=make_quadrature(base_shape%ndof, base_shape%dim, degree=4, family=FAMILY_COOLS)
     trace_shape=make_element_shape(base_shape%ndof, base_shape%dim, degree, quad, type=ELEMENT_TRACE)
     call report_test("[Derived loc]", &
      & trace_shape%ndof /= 4*(degree+1), .false., &
      & "Incorrect local node count")

     trace_mesh=make_mesh(base_mesh, trace_shape, continuity=-1, name="TraceMesh")
     call report_test("[Derived ele_count]", &
          & ele_count(trace_mesh) /= ele_count(base_mesh), .false., &
          & "Incorrect element count")

     ! test local coords
     allocate(quad_l_coords(base_shape%dim, trace_shape%ndof))
     quad_l_coords = quad_local_coords(degree)
     allocate(l_coords(base_shape%dim, trace_shape%ndof))
     fail = .false.
     ele_loop: do ele = 1, ele_count(trace_mesh)      
        fail = ele_loc(trace_mesh, ele) /= trace_shape%ndof
        if(fail) exit ele_loop
      
        do node = 1, size(l_coords, 2)
           l_coords(:, node) = local_coords(node, trace_shape%numbering)
        end do
        fail = fnequals(l_coords, quad_l_coords, tol = 1.0e3 * epsilon(0.0))
        if(fail) then
           do node = 1, size(l_coords, 2)
              print *, node, l_coords(:, node)
              print *, node, quad_l_coords(:, node)
           end do
           exit ele_loop
        end if
     end do ele_loop
     deallocate(l_coords)
     deallocate(quad_l_coords)

     ! setup field on base mesh and on trace mesh
     call allocate(D, positions%mesh, 'D')
     call set_from_python_function(D, "def val(X,t): return X[0]+X[1]", positions, 0.0)
     call allocate(L, trace_mesh, 'L')
     call allocate(trace_positions, positions%dim, trace_mesh, "TraceCoordinate")
     call remap_field(positions, trace_positions)
     call set_from_python_function(L, "def val(X,t): return X[0]+X[1]", trace_positions, 0.0)

     ! test face local nodes
     call test_face_local_nodes(L)

     ! test face local numbering



     ! test trace values

     ! test projection

    
    call report_test("[Derived mesh numbering]", fail, .false., "Invalid derived mesh numbering, failed on element " // int2str(ele))
    
    call deallocate(trace_shape)
    call deallocate(trace_mesh)
  end do
  
  call deallocate(positions)
  call report_test_no_references()
  
contains

  function quad_local_coords(degree) result(l_coords)
    !!< Return the node local coords
  
    integer, intent(in) :: degree
    
    integer :: i, index, j
    real, dimension(2, 4*(degree + 1)) :: l_coords
        
    index = 1
    ! face 1:
    do i=0, degree
       assert(index <= size(l_coords, 2))
       l_coords(1, index)=float(i) / float(degree)
       l_coords(2, index)=0.
       index = index + 1
    end do
    ! face 2:
    do i=0, degree
       assert(index <= size(l_coords, 2))
       l_coords(1, index)=1.
       l_coords(2, index)=float(i) / float(degree)
       index = index + 1
    end do
    ! face 3:
    do i=0, degree
       assert(index <= size(l_coords, 2))
       l_coords(1, index)=0.
       l_coords(2, index)=float(i) / float(degree)
       index = index + 1
    end do
    ! face 4:
    do i=0, degree
       assert(index <= size(l_coords, 2))
       l_coords(1, index)=float(i) / float(degree)
       l_coords(2, index)=1.
       index = index + 1
    end do

    assert(index == size(l_coords, 2) + 1)
    
  end function quad_local_coords

    subroutine test_face_local_nodes(L)
      implicit none
      type(scalar_field), intent(in) :: L

      integer :: ele

      do ele = 1, ele_count(L)
         call test_face_local_nodes_ele(L,ele)
      end do
      
    end subroutine test_face_local_nodes

    subroutine test_face_local_nodes_ele(L,ele)
      implicit none
      type(scalar_field), intent(in) :: L
      integer, intent(in) :: ele
      !
      integer, dimension(:), pointer :: neigh
      integer :: ni,ele2,face,face2
      
      neigh => ele_neigh(L,ele)
      do ni = 1, size(neigh)
         ele2 = neigh(ni)
         face = ele_face(L,ele,ele2)
         if(ele2>0) then
            face2 = ele_face(L,ele2,ele)
         else
            face2 = -1
         end if
         call test_face_local_nodes_face(L,ele,face,face2)
      end do
    end subroutine test_face_local_nodes_ele

    subroutine test_face_local_nodes_face(L,ele,face,face2)
      implicit none
      type(scalar_field), intent(in) :: L
      integer, intent(in) :: ele,face,face2
      !
      integer, dimension(face_loc(L,face)) :: nods1, nods2, nods3, nods_loc
      integer, dimension(:), pointer :: L_ele

      L_ele => ele_nodes(L,ele)
      nods1 = face_global_nodes(L,face)
      if(face2>0) then
         nods2 = face_global_nodes(L,face2)
      end if
      nods_loc = face_local_nodes(L,face)
      nods3 = L_ele(nods_loc)

      if(face2>0) then
         call report_test('[Global node numbers on face]', any(nods1/=nods2), .false., 'Global node numbers don''t agree on face')
      end if
      call report_test('[Face global node numbers]', any(nods1/=nods3), .false., 'Face global nodes doesn''t match global numbers')
    end subroutine test_face_local_nodes_face

end subroutine test_trace_space
