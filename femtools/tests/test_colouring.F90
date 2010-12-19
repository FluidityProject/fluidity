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
  
  subroutine test_colouring(sparsity, colour_sets)
  use fields_manipulation
  use state_module
  use vtk_interfaces
  use colouring
  use sparsity_patterns_meshes
  use unittest_tools
  implicit none

  type(state_type) :: state
  type(mesh_type), pointer  :: mesh, p0_mesh
  type(csr_sparsity), pointer :: sparsity
  type(integer_set), dimension(:), pointer, intent(out) :: colour_sets
  integer :: maxdgr, i
  logical :: fail


!  call vtk_read_state("data/anisotropic.vtu", state)
!  mesh => extract_mesh(state, "Mesh")
!  p0_mesh = piecewise_constant_mesh(mesh, "P0Mesh")
!  sparsity => get_csr_sparsity_compactdgdouble(state, p0_mesh)

!  maxdgr=0
!  do i=1, size(sparsity, 1)
!     maxdgr=max(maxdgr, row_length(sparsity, i))
!  enddo

!  call allocate(colour_sets)
!  call colour_sparsity(sparsity, colour_sets)

!  if (size(colour_sets) > maxdgr+1 ) fail = .true.
!  call report_test("colour sets", fail, .false., "colour sets should not be greater than max degree")
  
  end subroutine test_colouring
