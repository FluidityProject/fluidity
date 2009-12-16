#include "fdebug.h"
module gn_momentum_equation

  use elements
  use sparse_tools
  !use quadrature
  !use global_numbering
  !use shape_functions
  use global_parameters_gallopede
  use adjacency_lists
  use transform_elements
  use dgtools
  use solvers
  use data_structures
  use text_io
  use FETools
  use gallopede_solvers
  use vector_tools
  use gn_operator
  use vtk_io
  use fldebug
  
  implicit none

  !=============================================================
  !TODO
  !=============================================================
  !take internal boundary integrals out of get_mom
  !fix Dirichlet condition
  !=============================================================

  public :: solve_gn_momentum_equation

  private 

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"


contains

  subroutine solve_gn_momentum_equation(m,u,D,mesh,u_bcs)
    implicit none

    real, dimension(2*n_moms), target, intent(inout) :: m
    real, dimension(2*n_vels), target, intent(inout) :: u
    real, dimension(n_Dens), intent(in) :: D
    type(dg_mesh), intent(in) :: mesh
    type(bc_info), intent(in) :: u_bcs

    !local variables
    real, allocatable, dimension(:), target :: nonu, rhs
    real, dimension(:), pointer :: rhs1, rhs2
    real, dimension(:), pointer :: u1,u2,nonu1,nonu2, m1, m2

    integer :: nits, globi, globj
    type(block_csr_matrix) :: mommat
    type(csr_matrix) :: mat11, mat12, mat21, mat22
    real :: val1

    KSPType :: ksp_type
    PCType :: pc_type

    ewrite(1,*)("Subroutine solve_gn_momentum_equation")

    ewrite(2,*) shape(mesh%nu%dn)

    ewrite(2,*)("allocating rhs")
    allocate( rhs(2*N_moms) )
    rhs1 => rhs(1:N_moms)
    rhs2 => rhs(N_moms+1:2*N_moms)

    u1 => u(1:N_vels)
    u2 => u(N_vels+1:2*N_vels)
    m1 => m(1:N_moms)
    m2 => m(N_moms+1:2*N_moms)

    ewrite(2,*)("allocating nonu")
    allocate( nonu(2*N_vels) )
    nonu1 => nonu(1:N_vels)
    nonu2 => nonu(N_vels+1:2*N_vels)

    ewrite(2,*)("Constructing mommat")
    mommat = block_clone(mesh%Mass_m, (/2, 2/), type=CSR_REAL)
    mat11 = block(mommat,1,1)
    mat12 = block(mommat,1,2)
    mat21 = block(mommat,2,1)
    mat22 = block(mommat,2,2)

    ewrite(2,*)(size(mat11))
    ewrite(2,*)(size(mat12))
    ewrite(2,*)(size(mat21))
    ewrite(2,*)(size(mat22))

    nonu = u

    nits = 0

    ewrite(2,*)("starting loop")

    nlinear_loop: do
       if(nits == mom_maxnits) exit

       nits = nits + 1
       ewrite(2,*)(nits)
       ewrite(2,*)(mom_maxnits)

       call assemble_gn_momentum_equation(mat11,mat12,mat21,mat22, &
           rhs1,rhs2, &
           m1, m2, &
           0.5*(nonu1+u1),0.5*(nonu2+u2), &
           D,mesh)

       !solve GN equation
       ksp_type = KSPGMRES
       pc_type = PCSOR

       !need to check that we really have coupling
       !call gallopede_block_solve(m, mommat, rhs, &
       !     ksp_type, pc_type, 1.0e-10, 5000)
       call gallopede_solve(m1,mat11,rhs1, ksp_type, &
            pc_type,1.0e-10,5000)
       call gallopede_solve(m2,mat22,rhs2, ksp_type, &
            pc_type,1.0e-10,5000)
    

       ewrite(2,*)(sum(m))
       ewrite(2,*)(sum(D))

       call get_vels_GN(Mesh, nonu, m, D, u_bcs)

       ewrite(2,*) 'max(u)=',maxval(abs(u))

    end do nlinear_loop

    u = nonu

    deallocate( rhs )
    deallocate( nonu )
    ewrite(2,*)("deallocating mommat")
    call deallocate( mommat )
    
    ewrite(1,*)("END Subroutine solve_gn_momentum_equation")

  end subroutine solve_gn_momentum_equation

  subroutine assemble_gn_momentum_equation(mat11,mat12,mat21,mat22, &
       rhs1,rhs2,m1,m2,u1,u2, &
       D,mesh)
    implicit none
    real, intent(out), dimension(n_moms) :: rhs1,rhs2
    real, intent(in), dimension(n_vels) :: u1,u2
    real, intent(in), dimension(n_moms) :: m1,m2
    real, intent(in), dimension(N_dens) :: D
    type(csr_matrix) :: mat11, mat12, mat21, mat22
    type(dg_mesh), intent(in) :: mesh

    !locals
    integer :: ele, iloc, jloc, globi, globj
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele,m_ele_2
    integer, dimension(:), pointer :: u_ele_2,h_ele_2,X_ele_2
    real, dimension(mesh%nu%loc,2) :: mloc
    real, dimension(mesh%nu%loc,2) :: mloc_2
    real, dimension(mesh%nu%loc,2,2) :: graduloc_2
    real, dimension(mesh%nu%ngi,2) :: mlocgi
    real, dimension(mesh%nu%ngi,2,2) :: gradulocgi,gradulocgi_2
    real, dimension(mesh%nu%ngi,2) :: ulocgi
    real, dimension(mesh%nh%ngi) :: dlocgi
    real, dimension(mesh%nh_f%ngi) :: dlocgi_f
    real, dimension(mesh%nu_f%ngi,2) :: mlocgi_f, ulocgi_f
    real, dimension(2,2) :: gradulocgi_f,gradulocgi_f_2

    ! List of neighbours of current element.
    integer, dimension(:), allocatable :: neigh
    integer, dimension(:), pointer :: bdyu, bdyu_2, bdyh, bdyh_2,bdym,bdym_2
    integer, dimension(:), pointer :: bdyx
    integer :: bdy_i, ele_2, gi, nod
    real, dimension(2,mesh%nm%loc) :: ele_X, ele_X_2
    real, dimension(2,mesh%nm_f%loc) :: ele_Xf, ele_Xf_2
    real :: kmat, kmat11, kmat12, kmat21, kmat22
    integer :: ni
    real :: detwei(mesh%nm%ngi)
    real :: detwei_2(mesh%nm%ngi), detwei_f(mesh%nm_f%ngi)
    real :: normal(2,mesh%nm_f%ngi)
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t,dnu_t_2
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t,dnm_t_2
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t,dnh_t_2
    real, dimension(mesh%nu_f%ngi) :: udotnE, udotnF
    integer :: i,j

    real, dimension(mesh%nm%loc,mesh%nm%loc) :: Q
    real, dimension(2,mesh%nm%loc,mesh%nh%loc) :: dQmh
    real, dimension(mesh%nm_f%loc,mesh%nm_f%loc) :: Qf

    real, dimension(2,2,mesh%nm%ngi) :: Stress !stress tensor
    real, dimension(2,2,mesh%nm_f%ngi) :: Flux_1, Flux_2 !stress tensor

    integer, dimension(mesh%nh_f%ngi), target :: surface_h_lno
    integer, dimension(mesh%nm_f%ngi), target :: surface_m_lno
    integer, dimension(mesh%nu_f%ngi), target :: surface_u_lno

    logical :: PRESSURE_PARTS=.true.

    real :: avg_cst

    ewrite(1,*)("subroutine assemble_gn_momentum_equation")

    call zero(mat11)
    call zero(mat12)
    call zero(mat21)
    call zero(mat22)

    rhs1 = 0.
    rhs2 = 0.

    ! Initial (and for meshes with only one shape of element, final) size
    ! of neigh.
    ewrite(2,*)("Allocating memory for neigh");
    allocate(neigh(row_length(mesh%bdy_list,1)))

    ele_loop: do ele = 1, n_elements

       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)

       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)

       !call set_global_debug_level(4);

       call transform_to_physical(ele_X, mesh%nm, m = mesh%nu, dn_t = dnm_t, &
            dm_t = dnu_t, detwei = detwei)
       call transform_to_physical(ele_X, mesh%nx, m = mesh%nh, dm_t = dnh_t)

       !construct height on gauss points
       Dlocgi = matmul(transpose(mesh%nh%n),D(h_ele))
       !construct vels on gauss points
       ulocgi(:,1) = matmul(transpose(mesh%nu%n),u1(u_ele))
       ulocgi(:,2) = matmul(transpose(mesh%nu%n),u2(u_ele))
       !construct vels gradient on gauss points
       gradulocgi(:,1,1) = matmul(transpose(dnu_t(:,:,1)),u1(u_ele))
       gradulocgi(:,1,2) = matmul(transpose(dnu_t(:,:,2)),u1(u_ele))
       gradulocgi(:,2,1) = matmul(transpose(dnu_t(:,:,1)),u2(u_ele))
       gradulocgi(:,2,2) = matmul(transpose(dnu_t(:,:,2)),u2(u_ele))

       !VOLUME INTEGRALS

       !==========================================================
       !Mass matrix

       Q  = shape_shape(mesh%nm,mesh%nm,detwei) 

       rhs1(m_ele) = rhs1(m_ele) + matmul(Q,m1(m_ele))
       rhs2(m_ele) = rhs2(m_ele) + matmul(Q,m2(m_ele))

       call addto(mat11,m_ele,m_ele,Q)
       call addto(mat22,m_ele,m_ele,Q)

       !==========================================================
       !Advection bulk integral

       if(ADVECTION_FLAG) then

          ! div(um)
          !SIGN CORRECT

          Q = dshape_dot_vector_shape(dnm_t,transpose(ulocgi),mesh%nm,detwei)

          rhs1(m_ele) = rhs1(m_ele) - (1-theta)*dt*matmul(Q,m1(m_ele))
          rhs2(m_ele) = rhs2(m_ele) - (1-theta)*dt*matmul(Q,m2(m_ele))

          call addto(mat11,m_ele,m_ele,dt*theta*Q)
          call addto(mat22,m_ele,m_ele,dt*theta*Q)
       end if

       !============================================================
       !============================================================
       !nonlinear rhs contributions - bulk

       stress = 0.

       if(NONLINEAR_FLAG) then

          !-d_j dl/du^j_k u_{i,k}
          != - d_j ( D^3/3 div u delta_{j,k} u_{k,i} )
          != - d_j ( D^3/3 div u u_{j,i} )
          !-div dldgradu part -- SIGN CORRECT
          ! d_i ( - D^3/3 (div u)^2)
          !minus signs because we integrated by parts and 
          !moved to the other side of the equation
          !D-stress term -- SIGN CORRECT
          !minus signs because integrate by parts and move to other side
          !factor of 1/6-1/2 = -1/3

          !construct stress tensor

          stress(1,1,:) = stress(1,1,:)-dlocgi*dlocgi*dlocgi* &
               (gradulocgi(:,1,1) + gradulocgi(:,2,2)) * &
               gradulocgi(:,1,1)/3.0
          stress(1,2,:) = stress(1,2,:)-dlocgi*dlocgi*dlocgi* &
               (gradulocgi(:,1,1) + gradulocgi(:,2,2)) * &
               gradulocgi(:,2,1)/3.0
          stress(2,1,:) = stress(2,1,:)-dlocgi*dlocgi*dlocgi* &
               (gradulocgi(:,1,1) + gradulocgi(:,2,2)) * &
               gradulocgi(:,1,2)/3.0
          stress(2,2,:) = stress(2,2,:)-dlocgi*dlocgi*dlocgi* &
               (gradulocgi(:,1,1) + gradulocgi(:,2,2)) * &
               gradulocgi(:,2,2)/3.0
          stress(1,1,:) = stress(1,1,:)-dlocgi*dlocgi*dlocgi* &
               (gradulocgi(:,1,1) + gradulocgi(:,2,2)) * &
               (gradulocgi(:,1,1) + gradulocgi(:,2,2))/3.0
          stress(2,2,:) = stress(2,2,:)-dlocgi*dlocgi*dlocgi* &
               (gradulocgi(:,1,1) + gradulocgi(:,2,2)) * &
               (gradulocgi(:,1,1) + gradulocgi(:,2,2))/3.0
       end if

       !==============================================================
       !==============================================================
       !Pressure gradient bulk terms

       if(PRESSURE_PARTS.and.PRESSURE_FLAG) then
          stress(1,1,:) = stress(1,1,:) + 0.5*g0*dlocgi*dlocgi
          stress(2,2,:) = stress(2,2,:) + 0.5*g0*dlocgi*dlocgi
       end if

       !combine stresses on right-hand side

       !ewrite(2,*) stress

       rhs1(m_ele) = rhs1(m_ele) + &
            dt*(matmul(dnm_t(:,:,1),detwei * Stress(1,1,:)) + &
            matmul(dnm_t(:,:,2),detwei * Stress(2,1,:)))
       rhs2(m_ele) = rhs2(m_ele) + &
            dt*(matmul(dnm_t(:,:,1),detwei * Stress(1,2,:)) + &
            matmul(dnm_t(:,:,2),detwei * Stress(2,2,:)))

       if((.not.PRESSURE_PARTS).and.PRESSURE_FLAG) then
 
          dQmh = shape_dshape(mesh%nm,dnh_t,detwei)  
          
          !pressure w/o integrating by parts
          rhs1(m_ele) = rhs1(m_ele) - &
               dt*0.5*g0*matmul(dQmh(1,:,:),D(h_ele)*D(h_ele))
          rhs2(m_ele) = rhs2(m_ele) - &
               dt*0.5*g0*matmul(dQmh(2,:,:),D(h_ele)*D(h_ele))
       end if

       !==============================================================

       !SURFACE INTEGRALS

       !==============================================================

       if (size(neigh)/=row_length(mesh%bdy_list,ele)) then
          ewrite(2,*)("reallocating neigh")
          deallocate(neigh)
          allocate(neigh(row_length(mesh%bdy_list,ele)))
       end if

       neigh=row_m(mesh%bdy_list,ele)

       !surface integrals
       neighbourloop2: do ni=1,size(neigh)
          ele_2=neigh(ni)

          ! check for external bdy
          if (ele_2==0) then
             if(.false.) then
                cycle neighbourloop2
             else
                avg_cst = 1.0
                surface_h_lno = &
                     boundary_local_num(offnods(ni,:),mesh%nh%numbering)
                bdyh => surface_h_lno
                surface_u_lno = &
                     boundary_local_num(offnods(ni,:),mesh%nu%numbering)
                bdyu => surface_u_lno
                !surface_m_lno = &
                !     boundary_local_num(offnods(ni,:),mesh%nm%numbering)
                !bdym => surface_m_lno
                bdym => offnods(ni,:)
                bdyx => offnods(ni,:)
             end if
          else

             avg_cst = 0.5
             u_ele_2=>mesh%EVList_u((ELE_2-1)*mesh%nu%LOC+1:ELE_2*mesh%nu%LOC)
             m_ele_2=>mesh%EVList_m((ELE_2-1)*mesh%nm%LOC+1:ELE_2*mesh%nm%LOC)
             h_ele_2=>mesh%EVList_h((ELE_2-1)*mesh%Nh%LOC+1:ELE_2*mesh%Nh%LOC)
             X_ele_2=>mesh%EVList_X((ELE_2-1)*3+1:ELE_2*3)

             bdy_i=ival(mesh%bdy_list, ele, ele_2)
             bdyu=> mesh%bdy_nu_lno( &
                  (bdy_i-1)*mesh%nu_f%loc+1:bdy_i*mesh%nu_f%loc)
             bdym=> mesh%bdy_nm_lno( &
                  (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)
             bdyh=>mesh%bdy_nh_lno( &
                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)

             bdy_i=ival(mesh%bdy_list, ele_2, ele)
             bdyu_2=> mesh%bdy_nu_lno( &
                  (bdy_i-1)*mesh%nu_f%loc+1:bdy_i*mesh%nu_f%loc)
             bdym_2=> mesh%bdy_nm_lno( &
                  (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)
             bdyh_2=>mesh%bdy_nh_lno( &
                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)

             ! Locations of local vertices.
             ele_X_2(1,:)=mesh%X(X_ele_2)
             ele_X_2(2,:)=mesh%Y(X_ele_2)

          end if

          ! Locations of bdy vertices.
          ele_Xf=ele_X(:,bdym)

          if(ele_2.ne.0) then
             call transform_to_physical(ele_X, mesh%nm, &
                  m = mesh%nu, dn_t = dnm_t_2, &
                  dm_t = dnu_t_2, detwei = detwei_2)
             gradulocgi_2(:,1,1) = &
                  matmul(transpose(dnu_t_2(:,:,1)),u1(u_ele_2))
             gradulocgi_2(:,1,2) = &
                  matmul(transpose(dnu_t_2(:,:,2)),u1(u_ele_2))
             gradulocgi_2(:,2,1) = &
                  matmul(transpose(dnu_t_2(:,:,1)),u2(u_ele_2))
             gradulocgi_2(:,2,2) = &
                  matmul(transpose(dnu_t_2(:,:,2)),u2(u_ele_2))
          end if

          ! Change of coordinates on face
          call transform_bdy_to_physical(ele_X, ele_Xf, mesh%nm, mesh%nm_f, &
               detw_f = detwei_f,normal = normal)

          !==========================================================
          !Advection surface integral

          if(ADVECTION_FLAG) then
             if(ele_2.ne.0) then
                !inside face

                ulocgi_f(:,1) = matmul(transpose(mesh%nu_f%n),u1(u_ele(bdyu)))
                ulocgi_f(:,2) = matmul(transpose(mesh%nu_f%n),u2(u_ele(bdyu)))

                Qf = shape_shape(mesh%nm_f,mesh%nm_f,detwei_f*( &
                     ulocgi_f(:,1)*normal(1,:) + &
                     ulocgi_f(:,2)*normal(2,:)))

                call addto(mat11,m_ele(bdym),m_ele(bdym),0.5*dt*theta*Qf)
                call addto(mat22,m_ele(bdym),m_ele(bdym),0.5*dt*theta*Qf)

                rhs1(m_ele(bdym)) = rhs1(m_ele(bdym)) - &
                     0.5*dt*(1-theta)*matmul(Qf,m1(m_ele(bdym)))
                rhs2(m_ele(bdym)) = rhs2(m_ele(bdym)) - &
                     0.5*dt*(1-theta)*matmul(Qf,m2(m_ele(bdym)))

                !outside face

                ulocgi_f(:,1) = &
                     matmul(transpose(mesh%nu_f%n),u1(u_ele_2(bdyu_2)))
                ulocgi_f(:,2) = &
                     matmul(transpose(mesh%nu_f%n),u2(u_ele_2(bdyu_2)))

                Qf = shape_shape(mesh%nm_f,mesh%nm_f,detwei_f*( &
                     ulocgi_f(:,1)*normal(1,:) + &
                     ulocgi_f(:,2)*normal(2,:)))

                call addto(mat11,m_ele(bdym),m_ele_2(bdym_2),0.5*dt*theta*Qf)
                call addto(mat22,m_ele(bdym),m_ele_2(bdym_2),0.5*dt*theta*Qf)

                rhs1(m_ele(bdym)) = rhs1(m_ele(bdym)) - &
                     0.5*dt*(1-theta)*matmul(Qf,m1(m_ele_2(bdym_2)))
                rhs2(m_ele(bdym)) = rhs2(m_ele(bdym)) - &
                     0.5*dt*(1-theta)*matmul(Qf,m2(m_ele_2(bdym_2)))
             end if
          end if

          flux_1 = 0.
          flux_2 = 0.
          
          !============================================================
          !nonlinear rhs contributions - surface

          if(NONLINEAR_FLAG) then

             !construct height on gauss points on surface
             Dlocgi_f = matmul(transpose(mesh%nh_f%n),D(h_ele(bdyh)))
             
             !inside face

             !construct grad u on gauss points on surface
             !we just take an average here
             forall(i=1:2,j=1:2)
                gradulocgi_f(i,j) = sum(detwei*gradulocgi(:,i,j))/ &
                     sum(detwei)
             end forall
             
             flux_1(1,1,:) = flux_1(1,1,:) - &
                  dlocgi_f**3*(gradulocgi_f(1,1) + gradulocgi_f(2,2)) &
                  * (gradulocgi_f(1,1) + gradulocgi_f(2,2) &
                  + gradulocgi_f(1,1))
             flux_1(1,2,:) = flux_1(1,2,:) - &
                  dlocgi_f**3*(gradulocgi_f(1,1) + gradulocgi_f(2,2)) &
                  * gradulocgi_f(2,1)
             flux_1(2,1,:) = flux_1(2,1,:) - &
                  dlocgi_f**3*(gradulocgi_f(1,1) + gradulocgi_f(2,2)) &
                  * gradulocgi_f(1,2)
             flux_1(2,2,:) = flux_1(2,2,:) - &
                  dlocgi_f**3*(gradulocgi_f(1,1) + gradulocgi_f(2,2)) &
                  * (gradulocgi_f(1,1) + gradulocgi_f(2,2) &
                  + gradulocgi_f(2,2))
             
             !outside face
             if(ele_2.ne.0) then
                Dlocgi_f = matmul(transpose(mesh%nh_f%n),D(h_ele_2(bdyh_2)))
                
                !construct grad u on gauss points on surface
                !we just take an average here
                forall(i=1:2,j=1:2)
                   gradulocgi_f(i,j) = sum(detwei*gradulocgi_2(:,i,j))/ &
                        sum(detwei_2)
                end forall
                
                flux_2(1,1,:) = flux_2(1,1,:) - &
                     dlocgi_f**3*(gradulocgi_f(1,1) + gradulocgi_f(2,2)) &
                     * (gradulocgi_f(1,1) + gradulocgi_f(2,2) &
                     + gradulocgi_f(1,1))
                flux_2(1,2,:) = flux_2(1,2,:) - &
                     dlocgi_f**3*(gradulocgi_f(1,1) + gradulocgi_f(2,2)) &
                     * gradulocgi_f(2,1)
                flux_2(2,1,:) = flux_2(2,1,:) - &
                     dlocgi_f**3*(gradulocgi_f(1,1) + gradulocgi_f(2,2)) &
                     * gradulocgi_f(1,2)
                flux_2(2,2,:) = flux_2(2,2,:) - &
                     dlocgi_f**3*(gradulocgi_f(1,1) + gradulocgi_f(2,2)) &
                     * (gradulocgi_f(1,1) + gradulocgi_f(2,2) &
                     + gradulocgi_f(2,2))
                
             end if
          end if

          !==============================================================
          !==============================================================

          ewrite(4,*) 'Pressure gradient surface term'

          if(PRESSURE_FLAG.and.PRESSURE_PARTS) then
             Dlocgi_f = matmul(transpose(mesh%nh_f%n),D(h_ele(bdyh)))
             flux_1(1,1,:) = flux_1(1,1,:) + 0.5*g0*dlocgi_f*dlocgi_f
             flux_1(2,2,:) = flux_1(2,2,:) + 0.5*g0*dlocgi_f*dlocgi_f
             flux_2(1,1,:) = flux_2(1,1,:) + 0.5*g0*dlocgi_f*dlocgi_f
             flux_2(2,2,:) = flux_2(2,2,:) + 0.5*g0*dlocgi_f*dlocgi_f
          end if

          ewrite(4,*) 'Flux integral'
          
          !inside face
          rhs1(m_ele(bdym)) = rhs1(m_ele(bdym)) - &
               dt*avg_cst*shape_rhs(mesh%nm_f, &
               (normal(1,:)*flux_1(1,1,:)+ &
               normal(2,:)*flux_1(2,1,:) )*detwei_f)
          rhs2(m_ele(bdym)) = rhs2(m_ele(bdym)) - &
               dt*avg_cst*shape_rhs(mesh%nm_f, &
               (normal(1,:)*flux_1(1,2,:)+ &
               normal(2,:)*flux_1(2,2,:) )*detwei_f)
          !outside face
          rhs1(m_ele(bdym)) = rhs1(m_ele(bdym)) - &
               dt*(1-avg_cst)*shape_rhs(mesh%nm_f, &
               (normal(1,:)*flux_2(1,1,:)+ &
               normal(2,:)*flux_2(2,1,:) )*detwei_f)
          rhs2(m_ele(bdym)) = rhs2(m_ele(bdym)) - &
               dt*(1-avg_cst)*shape_rhs(mesh%nm_f, &
               (normal(1,:)*flux_2(1,2,:)+ &
               normal(2,:)*flux_2(2,2,:) )*detwei_f)

          !ewrite(2,*)("next neigh")

       end do neighbourloop2

    end do ele_loop
    deallocate( neigh )

    ewrite(1,*)("END subroutine assemble_gn_momentum_equation")

  end subroutine assemble_gn_momentum_equation

end module gn_momentum_equation
