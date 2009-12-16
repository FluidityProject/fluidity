#include "fdebug.h"

module GN_operator
  
  use elements
  use sparse_tools
!  use quadrature
!  use global_numbering
!  use shape_functions
  use global_parameters_gallopede
!  use adjacency_lists
  use transform_elements
!  use dgtools
!  use text_io
  use data_structures
  use fetools
!  use vector_tools
  use gallopede_solvers
!  use vtk_io
!  use mesh_tools
  use fldebug
  
  implicit none

  public :: assemble_GN_alpha_operator, get_vels_alpha
  public :: assemble_GN_operator, get_vels_GN

  private
  
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
      
contains
  
  subroutine assemble_GN_alpha_operator(mom11,mom12,mom21,mom22, &
       rhs1,rhs2,m1,m2,mesh,u1,u2)
    implicit none
    type(csr_matrix), intent(inout) :: mom11, mom12, mom21, mom22
    real, dimension(n_vels), intent(inout) :: rhs1, rhs2
    real, dimension(n_moms), intent(in) :: m1, m2
    type(dg_mesh) :: mesh
    real, dimension(n_vels), intent(in), optional :: u1,u2

    !locals
    integer :: ele !index for elements
    real, dimension(2,3) :: ele_X !coordinates for triangle vertices
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t !gradient of shape fn
    real, dimension(mesh%nu%ngi) :: detwei !Jacobian * weights
    integer, dimension(:), pointer :: u_ele,X_ele,m_ele
    real, dimension(mesh%nu%loc,mesh%nu%loc) :: Q
    real, dimension(2,2,mesh%nu%loc,mesh%nu%loc) :: QT
    real, dimension(mesh%nu%loc,mesh%nm%loc) :: Qrhs

    ewrite(1,*)("subroutine assemble_GN_alpha_operator")

    call zero(mom11)
    call zero(mom12)
    call zero(mom21)
    call zero(mom22)

    if(present(u1)) then 
       rhs1 = 0.
       rhs2 = 0.
    end if

    ele_loop: do ele = 1, n_elements

       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)       
       m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)

       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)

       call transform_to_physical(ele_X, mesh%nm, &
            dn_t = dnu_t, detwei = detwei)

       !Right-hand side
       if(present(u1)) then
          Qrhs  = shape_shape(mesh%nu,mesh%nm,detwei)
          rhs1(u_ele) = rhs1(u_ele) + matmul(Qrhs,m1(m_ele))
          rhs2(u_ele) = rhs2(u_ele) + matmul(Qrhs,m2(m_ele))
       end if

       !Matrix
       Q = shape_shape(mesh%nu,mesh%nu,detwei)
       QT = dshape_outer_dshape(dnu_t,dnu_t,detwei)
       call addto(mom11,u_ele,u_ele,Q+alpha*alpha*QT(1,1,:,:))
       call addto(mom12,u_ele,u_ele,alpha*alpha*QT(1,2,:,:))
       call addto(mom21,u_ele,u_ele,alpha*alpha*QT(2,1,:,:))
       call addto(mom22,u_ele,u_ele,Q+alpha*alpha*QT(2,2,:,:))

    end do ele_loop

    ewrite(1,*)("end subroutine assemble_GN_alpha_operator")

  end subroutine assemble_GN_alpha_operator

  subroutine get_vels_alpha(Mesh, u, m, bcs)
    type(dg_mesh), intent(in) :: Mesh
    real, dimension(N_vels*2), target, intent(inout) :: u
    real, dimension(N_moms*2), target, intent(in) :: m
    type(bc_info), intent(in) :: bcs

    !locals
    real, dimension(:), pointer :: u1,u2,m1,m2,rhs1,rhs2
    real, allocatable, dimension(:), target :: rhs
    type(csr_matrix) :: mom11, mom12, mom21, mom22


    KSPType :: ksp_type
    PCType :: pc_type

    ksp_type = KSPCG
    PC_type = PCICC

    ewrite(1,*) 'subroutine get_vels_alpha'

    allocate( rhs(N_vels*2) )

    mom11 = clone(mesh%Mass_u)
    mom12 = clone(mesh%Mass_u)
    mom21 = clone(mesh%Mass_u)
    mom22 = clone(mesh%Mass_u)

    u1 => u(1:N_vels)
    u2 => u(N_vels+1:2*N_vels)
    rhs1 => rhs(1:N_vels)
    rhs2 => rhs(N_vels+1:2*N_vels)
    m1 => m(1:N_moms)
    m2 => m(N_moms+1:2*N_moms)

    call assemble_GN_alpha_operator(mom11,mom12,mom21,mom22, &
         rhs1,rhs2,m1,m2,mesh,u1,u2)  
    
    call gallopede_solve_lift(u1,u2, &
         mom11, mom12, mom21, mom22, &
         rhs1, rhs2, &
         bcs, &
         ksp_type, &
         pc_type, 1.0e-10, 1000)

    deallocate ( rhs )
    call deallocate( mom11 )    
    call deallocate( mom12 )
    call deallocate( mom21 )    
    call deallocate( mom22 )

    ewrite(1,*) 'END subroutine get_vels_alpha'

  end subroutine get_vels_alpha

  subroutine assemble_GN_operator(mom11,mom12,mom21,mom22, &
       rhs1,rhs2,m1,m2,D,mesh,u1,u2)
    implicit none
    type(csr_matrix), intent(inout) :: mom11, mom12, mom21, mom22
    real, dimension(n_vels), intent(inout) :: rhs1, rhs2
    real, dimension(n_moms), intent(in) :: m1, m2
    real, dimension(n_dens), intent(in) :: D
    type(dg_mesh) :: mesh
    real, dimension(n_vels), intent(in), optional :: u1,u2

    !locals
    integer :: ele !index for elements
    real, dimension(2,3) :: ele_X !coordinates for triangle vertices
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t !gradient of shape fn
    real, dimension(mesh%nu%ngi) :: detwei !Jacobian * weights
    integer, dimension(:), pointer :: u_ele,X_ele,m_ele,h_ele
    real, dimension(mesh%nu%loc,mesh%nu%loc) :: Q
    real, dimension(2,2,mesh%nu%loc,mesh%nu%loc) :: QT
    real, dimension(mesh%nu%loc,mesh%nm%loc) :: Qrhs
    real, dimension(mesh%nh%ngi) :: dlocgi

    ewrite(1,*)("subroutine assemble_GN_operator")

    call zero(mom11)
    call zero(mom12)
    call zero(mom21)
    call zero(mom22)

    if(present(u1)) then 
       rhs1 = 0.
       rhs2 = 0.
    end if

    ele_loop: do ele = 1, n_elements

       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)       
       m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)

       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)
       call transform_to_physical(ele_X, mesh%nm, m= mesh%nu, &
            dm_t = dnu_t, detwei = detwei)

       !Right-hand side
       if(present(u1)) then
          Qrhs  = shape_shape(mesh%nu,mesh%nm,detwei)
          rhs1(u_ele) = rhs1(u_ele) + matmul(Qrhs,m1(m_ele))
          rhs2(u_ele) = rhs2(u_ele) + matmul(Qrhs,m2(m_ele))
       end if

       !Matrix
       dlocgi = matmul(transpose(mesh%nh%n),D(h_ele))
       Q = shape_shape(mesh%nu,mesh%nu,detwei*Dlocgi)
       QT = dshape_outer_dshape(dnu_t,dnu_t,detwei*Dlocgi*Dlocgi*Dlocgi/3.0)
       call addto(mom11,u_ele,u_ele,Q+QT(1,1,:,:))
       call addto(mom12,u_ele,u_ele,QT(1,2,:,:))
       call addto(mom21,u_ele,u_ele,QT(2,1,:,:))
       call addto(mom22,u_ele,u_ele,Q+QT(2,2,:,:))

    end do ele_loop

    ewrite(1,*)("end subroutine assemble_GN_operator")

  end subroutine assemble_GN_operator

  subroutine get_vels_GN(Mesh, u, m, D, bcs)
    type(dg_mesh), intent(in) :: Mesh
    real, dimension(N_vels*2), target, intent(inout) :: u
    real, dimension(N_moms*2), target, intent(in) :: m
    real, dimension(N_dens), target, intent(in) :: D
    type(bc_info), intent(in) :: bcs

    !locals
    real, dimension(:), pointer :: u1,u2,m1,m2,rhs1,rhs2
    real, allocatable, dimension(:), target :: rhs
    type(csr_matrix) :: mom11, mom12, mom21, mom22


    KSPType :: ksp_type
    PCType :: pc_type

    ksp_type = KSPCG
    PC_type = PCICC

    ewrite(1,*) 'subroutine get_vels_GN'

    allocate( rhs(N_vels*2) )

    mom11 = clone(mesh%Mass_u)
    mom12 = clone(mesh%Mass_u)
    mom21 = clone(mesh%Mass_u)
    mom22 = clone(mesh%Mass_u)

    u1 => u(1:N_vels)
    u2 => u(N_vels+1:2*N_vels)
    rhs1 => rhs(1:N_vels)
    rhs2 => rhs(N_vels+1:2*N_vels)
    m1 => m(1:N_moms)
    m2 => m(N_moms+1:2*N_moms)

    call assemble_GN_operator(mom11,mom12,mom21,mom22, &
         rhs1,rhs2,m1,m2,D,mesh,u1,u2)  
    
    call gallopede_solve_lift(u1,u2, &
         mom11, mom12, mom21, mom22, &
         rhs1, rhs2, &
         bcs, &
         ksp_type, &
         pc_type, 1.0e-10, 1000)

    deallocate ( rhs )
    call deallocate( mom11 )    
    call deallocate( mom12 )
    call deallocate( mom21 )    
    call deallocate( mom22 )

    ewrite(1,*) 'END subroutine get_vels_GN'

  end subroutine get_vels_GN

end module GN_operator
