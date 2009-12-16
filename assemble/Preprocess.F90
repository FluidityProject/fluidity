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
module preprocess_module
  !! Calculate shape functions and sparsity patterns.
  use sparse_tools
  use sparsity_patterns
  use shape_functions
  use fields
  use state_module
  use fldebug
  use flcomms_module
  use global_parameters, only : halo_tag, halo_tag_p
  use sparsity_patterns_meshes
  implicit none

  private
  public form_matrix_sparsity

contains

  subroutine form_matrix_sparsity(state)
    !!< Generate the sparsity of the matrices used in the model
    !!< and insert in all states
    type(state_type), dimension(:), intent(inout):: state
    
    type(csr_sparsity) CT_sparsity, U_sparsity, CMC_sparsity
    type(mesh_type), pointer :: U_mesh, P_mesh
    integer ph
      
    U_mesh => extract_mesh(state(1), "VelocityMesh")
    P_mesh => extract_mesh(state(1), "PressureMesh")
    
    ! Generate sparsity patterns.
    ! First order for CG but second order for DG due to viscosity.
    if (U_mesh%continuity>=0) then
       U_sparsity=make_sparsity(U_mesh, U_mesh, "VelocitySparsity")
!!$       if(isparallel()) then
!!$         U_sparsity%halo_tag = halo_tag
!!$         U_sparsity%private_columns = get_nowned_nodes(halo_tag)
!!$         U_sparsity%private_rows = U_sparsity%private_columns
!!$       end if
    else
       U_sparsity=make_sparsity_transpose(U_mesh, U_mesh, "VelocitySparsity")
!!$       if(isparallel()) then
!!$         U_sparsity%halo_tag = halo_tag_p
!!$         U_sparsity%private_columns = get_nowned_nodes(halo_tag_p)
!!$         U_sparsity%private_rows = U_sparsity%private_columns
!!$       end if
    end if

    ! Sparsity of CMC in fluidity notation. This takes into account the
    ! second order pressure operator.
    CMC_sparsity=make_sparsity_transpose(P_mesh, U_mesh, "PressurePoissonSparsity")
!!$    if(isparallel()) then
!!$      CMC_sparsity%halo_tag = halo_tag_p
!!$      CMC_sparsity%private_columns = get_nowned_nodes(halo_tag_p)
!!$      CMC_sparsity%private_rows = CMC_sparsity%private_columns
!!$    end if

    ! Sparsity of C^T - the transpose of the pressure gradient operator.
    CT_sparsity=make_sparsity(P_mesh, U_mesh, "PressureGradientSparsity")

    do ph=1, size(state)    
       call insert(state(ph), U_sparsity, U_sparsity%name)
       call insert(state(ph), CMC_sparsity, CMC_sparsity%name)
       call insert(state(ph), CT_sparsity, CT_sparsity%name)
    end do
    
    ! drop our references to the sparsities
    call deallocate(U_sparsity)
    call deallocate(CMC_sparsity)
    call deallocate(CT_sparsity)
    
  end subroutine form_matrix_sparsity

end module preprocess_module
