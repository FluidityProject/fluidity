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
  use element_numbering, only: FAMILY_SIMPLEX
  use fields


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


  integer, public :: VTK_QUALITY_EDGE_RATIO = 0, &
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
       VTK_QUALITY_ASPECT_BETA = 29

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

    call mesh_quality_c(positions%dim, node_count(positions), ele_count(positions),&
         size(positions%mesh%ndglno), quality_measure,&
         positions%val, positions%mesh%ndglno,&
         s_field%val)

  end subroutine get_mesh_quality

end module mesh_quality
