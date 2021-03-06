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

subroutine test_1d_pickers

  use fields
  use fldebug
  use pickers
  use mesh_files
  use unittest_tools
  
  implicit none
  
  integer :: ele
  type(vector_field) :: positions

  positions = read_mesh_files("data/interval", quad_degree = 1, format="gmsh")

  call report_test("[Picker pointer allocated]", .not. associated(positions%picker), .false., "Picker pointer not allocated")
  call report_test("[No picker attached]", associated(positions%picker%ptr), .false., "Picker already attached")
  
  call picker_inquire(positions, (/-1.0 /), ele)
  call report_test("[Point not contained]", ele > 0, .false., "Incorrectly reported point contained in mesh")

  call report_test("[Picker cached]", .not. associated(positions%picker%ptr), .false., "Picker not cached")

  call picker_inquire(positions, (/ 0.25 /), ele)
  call report_test("[Point contained]", ele /= 3, .false., "Reported incorrect containing element")
  
  call deallocate(positions)
  
  call report_test_no_references()
  
end subroutine test_1d_pickers
