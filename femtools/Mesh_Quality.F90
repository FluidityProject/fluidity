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
module mesh_quality

  use iso_c_binding
  use FLdebug
  use vector_tools, only: trace_mat, det
  use element_numbering, only: FAMILY_SIMPLEX
  use fields

  implicit none

  interface
     subroutine mesh_quality_c(dim, n_nodes, n_elements, connectivity_len,&
       measure, points, connectivity, quality) bind(c)
       use iso_c_binding
       integer (c_int) :: dim, n_nodes, n_elements, connectivity_len, measure
       real (c_double) :: points(dim, n_nodes)
       integer (c_int) :: connectivity(connectivity_len)
       real (c_double) :: quality(n_elements)
     end subroutine mesh_quality_c
  end interface

  private

  public :: get_mesh_quality


  integer, public, parameter :: VTK_QUALITY_EDGE_RATIO = 0, &
       VTK_QUALITY_ASPECT_RATIO = 1, &
       VTK_QUALITY_RADIUS_RATIO = 2, &
       VTK_QUALITY_ASPECT_FROBENIUS = 3, &
       VTK_QUALITY_MED_ASPECT_FROBENIUS = 4, &
       VTK_QUALITY_MAX_ASPECT_FROBENIUS = 5, &
       VTK_QUALITY_MIN_ANGLE = 6, &
       VTK_QUALITY_COLLAPSE_RATIO = 1, &
       VTK_QUALITY_MAX_ANGLE = 8, &
       VTK_QUALITY_CONDITION = 9, &
       VTK_QUALITY_SCALED_JACOBIAN = 10, &
       VTK_QUALITY_SHEAR = 11, &
       VTK_QUALITY_RELATIVE_SIZE_SQUARED = 12, &
       VTK_QUALITY_SHAPE = 13, &
       VTK_QUALITY_SHAPE_AND_SIZE = 14, &
       VTK_QUALITY_DISTORTION = 15, &
       VTK_QUALITY_MAX_EDGE_RATIO = 16, &
       VTK_QUALITY_SKEW = 17, &
       VTK_QUALITY_TAPER = 18, &
       VTK_QUALITY_ASPECT_VOLUME = 19, &
       VTK_QUALITY_ASPECT_STRETCH = 20, &
       VTK_QUALITY_ASPECT_DIAGONAL = 21, &
       VTK_QUALITY_ASPECT_DIMENSION = 22, &
       VTK_QUALITY_ASPECT_ODDY = 23, &
       VTK_QUALITY_ASPECT_SHEAR_AND_SIZE = 24, &
       VTK_QUALITY_ASPECT_JACOBIAN = 25, &
       VTK_QUALITY_ASPECT_WARPAGE = 26, &
       VTK_QUALITY_ASPECT_GAMMA = 27, &
       VTK_QUALITY_AREA = 28, &
       VTK_QUALITY_ASPECT_BETA = 29,&
       FLUIDITY_QUALITY_SIMPLE_ANISOTROPIC = -1, &
       FLUIDITY_QUALITY_ALGEBRAIC_SHAPE_ONE = -2, &
       FLUIDITY_QUALITY_ALGEBRAIC_SHAPE_TWO = -3

contains

  subroutine get_mesh_quality(positions, s_field, quality_measure)
    integer, intent(inout) :: quality_measure
    type(vector_field), intent(in) :: positions 
    type(scalar_field), intent(inout) :: s_field

    assert(element_count(positions) == element_count(s_field))
    assert(node_count(s_field) == element_count(s_field))

    if (positions%mesh%shape%numbering%family /= FAMILY_SIMPLEX&
         .or. positions%mesh%shape%loc /= positions%mesh%shape%dim+1) then
       FLAbort("Trying to get mesh quality for a mesh which isn't linear simplicial. This isn't currently supported.")
    endif

    if (quality_measure>=0) then
       call mesh_quality_c(positions%dim, node_count(positions), ele_count(positions),&
            size(positions%mesh%ndglno), quality_measure,&
            positions%val, positions%mesh%ndglno,&
            s_field%val)
    else

       select case(quality_measure)
       case(FLUIDITY_QUALITY_SIMPLE_ANISOTROPIC)
          call simple_anisotropic_measure(positions,s_field)
       case(FLUIDITY_QUALITY_ALGEBRAIC_SHAPE_ONE)
          call algebraic_shape_one_measure(positions,s_field)
       case(FLUIDITY_QUALITY_ALGEBRAIC_SHAPE_TWO)
          call algebraic_shape_two_measure(positions,s_field)
       case default
          FLAbort("Unknown mesh quality function.")
       end select
    end if

  end subroutine get_mesh_quality

  subroutine simple_anisotropic_measure(positions, s_field, metric)

    !!! This function implements the quality measure used in papers by
    !!! Frederic Alauzet's group. 
    !!!
    !!! Q_E = C_N (\sum_i |e^E_i|^2)^(N/2)/V_E
    !!! 
    !!! where E is an element index, N the dimension, e_i are edge vectors,
    !!! V is the volume/area of the simplex and the C_N is chosen such that
    !!! Q_E = 1 for a regular simplex.
    !!!
    !!! See Mesh Generation: Application to Finite Elements by
    !!! Pascal Jean Frey & Paul-Louis George for references.

    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: s_field
    type(tensor_field), intent(in), optional :: metric

    integer :: ele, i, j, k
    real    :: n, coeff, quality
    real    :: X(mesh_dim(positions), ele_loc(positions,1))
    real    :: M(mesh_dim(positions),mesh_dim(positions))
    real    :: l(mesh_dim(positions)*(mesh_dim(positions)+1)/2)
    

    !!! This measure is trivial in 1d

    select case(mesh_dim(positions))
    case(1)
       call set(s_field,1.0)
       return
    case(2)
       n = 2.0
       coeff = 1.0/6.0
    case(3)
       n = 3.0
       coeff = sqrt(3.0)/216.0
    end select

    if (.not. present(metric)) then
       !!! if no metric is supplied, we use the identity
       M=0
       do i=1, size(M,1)
          M(i,i) = 1.0
       end do
    end if

    do ele=1,element_count(positions)

       if (present(metric)) then
          !!! average the metric over the element in case it's not P0
          M = sum(ele_val(metric, ele))/ele_loc(metric, ele)
       end if

       quality=0.0
       k=1
       X = ele_val(positions,ele)
       do i=1, ele_loc(positions,ele)
          do j=i+1,ele_loc(positions,ele)
             l(k)=length2(X(:,i),X(:,j),M)
             quality = quality + l(k)
             l(k)=sqrt(l(k))
             k=k+1
          end do
       end do

       quality = coeff*quality**(n/2.0)/cell_measure(L)

       call set(s_field, ele, quality)
    end do

    contains

      real function length2(p1,p2,M)
        real, intent(in) :: p1(:),p2(:),M(:,:)

        length2 = dot_product(matmul(p2-p1,M),p2-p1)
      end function length2

      real function cell_measure(lengths)
        real, dimension(:), intent(in) :: lengths
        real :: s
        
        select case(size(lengths))
        case(3)
           s= sum(lengths)/2.0
           !! In 2d use Heron's formula for the area
           cell_measure = sqrt(abs(s*(s-l(1))*(s-l(2))*(s-l(3))))
        case(6)
           !! In 3D use the (hideous) Cayley-Menger determinant.
           !!
           !! | 0    1     1     1     1   |
           !! | 1    0   d_1^2 d_2^2 d_3^2 |
           !! | 1  d_1^2   0   d_4^2 d_5^2 |  = 288V^2
           !! | 1  d_2^2 d_4^2   0   d_6^2 |
           !! | 1  d_3^2 d_5^2 d_6^2   0   |
           cell_measure = sqrt(abs(-2.0*(l(1)**2*l(6))**2&
                -2.0*(l(1)*l(2)*l(4))**2&
                +2.0*(l(1)*l(2)*l(5))**2&
                +2.0*(l(1)*l(3)*l(4))**2&
                -2.0*(l(1)*l(3)*l(5))**2&
                +2.0*(l(1)*l(3)*l(6))**2&
                +2.0*(l(1)*l(4)*l(6))**2&
                +2.0*(l(1)*l(5)*l(6))**2&
                -2.0*(l(1)*l(6)**2)**2&
                -2.0*(l(2)**2*l(5))**2&
                +2.0*(l(2)*l(3)*l(5))**2&
                -2.0*(l(2)*l(3)*l(6))**2&
                +2.0*(l(2)*l(3)*l(5))**2&
                -2.0*(l(2)*l(5)**2)**2&
                +2.0*(l(2)*l(5)*l(6))**2&
                -2.0*(l(3)**2*l(4))**2&
                -2.0*(l(3)*l(4)**2)**2&
                +2.0*(l(3)*l(4)*l(5))**2&
                +2.0*(l(3)*l(4)*l(6))**2&
                -2.0*(l(4)*l(5)*l(6))**2)/288.0)
        case default
        end select
           
      end function cell_measure

    end subroutine simple_anisotropic_measure

    subroutine algebraic_shape_one_measure(positions, s_field)

    !!! This function implements the quality measure

    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: s_field

    integer :: ele, i, j, dim

    integer, pointer, dimension(:) :: nodes

    real, dimension(mesh_dim(positions)) :: edge
    real, dimension(mesh_dim(positions), mesh_dim(positions)) :: m

    dim = mesh_dim(positions)

    do ele = 1, element_count(s_field)
       nodes => ele_nodes(positions, ele)
       m = 0.0
       do i = 1, ele_loc(positions, ele)-1
          do j = i, ele_loc(positions, ele)
             edge = node_val(positions, nodes(i)) - node_val(positions, nodes(j))

             m = m + outer_product(edge, edge)
          end do
       end do

       m = 2.0*m/(dim+1.0)

       call set(s_field, ele, dim*det(m)**(1.0/(1.0*dim))/trace_mat(m))

    end do

  end subroutine algebraic_shape_one_measure

    subroutine algebraic_shape_two_measure(positions, s_field)

    !!! This function implements the quality measure

    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: s_field

    integer :: ele, i, j, dim

    integer, pointer, dimension(:) :: nodes

    real :: q_k
    real, dimension(mesh_dim(positions)) :: edge
    real, dimension(mesh_dim(positions), mesh_dim(positions)) :: m

    dim = mesh_dim(positions)

    do ele = 1, element_count(s_field)
       nodes => ele_nodes(positions, ele)
       m = 0.0
       do i = 1, ele_loc(positions, ele)-1
          do j = i, ele_loc(positions, ele)
             edge = node_val(positions, nodes(i)) - node_val(positions, nodes(j))

             m = m + outer_product(edge, edge)
          end do
       end do

       q_k = trace_mat(m)
       call invert(m)
       q_k = dim**2/(q_k*trace_mat(m))

       call set(s_field, ele, q_k)

    end do

  end subroutine algebraic_shape_two_measure

end module mesh_quality
