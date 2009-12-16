#include "fdebug.h"
module data_structures
  use sparse_tools, only : csr_matrix,dynamic_csr_matrix,&
       block_dynamic_csr_matrix
  use elements, only : element_type
  use transform_elements
  use fields
  use dgtools
  use global_parameters_gallopede

  implicit none

  integer, parameter :: INTERIOR=1, SLIP=2, NOSLIP=0

  type dg_mesh
     integer :: n_element, N_verts, N_vels, N_moms, n_dens
     real, dimension(:), pointer :: X=>null(),Y=>null()
     integer, dimension(:), pointer :: EVlist_m=>null(), EVlist_h=>null()
     integer, dimension(:), pointer :: EVList_X=>null(), EVList_u=>null()
     integer, dimension(:), pointer :: EVList_cX=>null()
     type(csr_matrix) :: mass_u,mass_h,MLCM_mat_block, BARO_CMC,u_mat,D_mat,&
          bdy_list,sparse_m, mass_cx
     type(dynamic_csr_matrix) :: CMT,CMC,mass_u_inv
     type(block_dynamic_csr_matrix) :: M_inv,C
     type(vector_field), pointer :: positions=>null(),connectivity=>null()
!     type(csr_matrix) :: mass_field
     integer, dimension(:), pointer :: bdy_nu_lno=>null(), bdy_nh_lno=>null()
     integer, dimension(:), pointer :: bdy_nm_lno=>null(), bdy_nx_lno=>null()
     type(element_type), pointer :: nu, nh, nu_f, nh_f, nm, nm_f, nx
     
  end type dg_mesh

  type bc_info
     integer :: N_interior, N_tangents
     integer, dimension(:), pointer :: interior_list=>null()
     integer, dimension(:), pointer :: tangent_list=>null()
     real, dimension(:,:), pointer :: tangents=>null()
     integer, dimension(:), pointer :: bc_marker=>null()
     integer, dimension(:), pointer :: lifted_ordering=>null()
  end type bc_info

  contains

    subroutine get_bc_marker(bcs,n)
      type(bc_info), intent(inout) :: bcs
      integer, intent(in) :: n

      !locals
      integer :: i

      if(associated(bcs%bc_marker)) then
         deallocate(bcs%bc_marker)
         bcs%bc_marker => null()
      end if
      if(associated(bcs%lifted_ordering)) then
         deallocate(bcs%lifted_ordering)
         bcs%lifted_ordering=> null()
      end if
      allocate( bcs%bc_marker(n) )
      allocate( bcs%lifted_ordering(n) )
      bcs%bc_marker = 0

      ewrite(3,*) n
      ewrite(3,*) bcs%N_interior
      ewrite(3,*) bcs%N_tangents

      do i = 1,bcs%N_interior
         bcs%bc_marker(bcs%interior_list(i)) = 1
         bcs%lifted_ordering(bcs%interior_list(i)) = i
      end do
      do i = 1,bcs%N_tangents
         bcs%bc_marker(bcs%tangent_list(i)) = 2
         bcs%lifted_ordering(bcs%tangent_list(i)) = i
      end do

    end subroutine get_bc_marker

end module data_structures

