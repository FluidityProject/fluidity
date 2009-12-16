#include "fdebug.h"
module gn_momentum_equation

  use elements
  use sparse_tools
  use quadrature
  use global_numbering
  use shape_functions
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

  integer, parameter :: NO_BC = 0

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

  subroutine solve_gn_momentum_equation(m,u,D,mesh)
    implicit none

    real, dimension(2*n_vels), target, intent(inout) :: m
    real, dimension(2*n_vels), target, intent(inout) :: u
    real, dimension(n_Dens), intent(in) :: D
    type(dg_mesh), intent(in) :: mesh

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

    ewrite(2,*)("Subroutine solve_gn_momentum_equation")

    ewrite(2,*)("allocating rhs")
    allocate( rhs(2*N_vels) )
    rhs1 => rhs(1:N_vels)
    rhs2 => rhs(N_vels+1:2*N_vels)

    u1 => u(1:N_vels)
    u2 => u(N_vels+1:2*N_vels)
    m1 => m(1:N_vels)
    m2 => m(N_vels+1:2*N_vels)

    ewrite(2,*)("allocating nonu")
    allocate( nonu(2*N_vels) )
    nonu1 => nonu(1:N_vels)
    nonu2 => nonu(N_vels+1:2*N_vels)

    ewrite(2,*)("Constructing mommat")
    mommat = block_clone(mesh%Mass_u, (/2, 2/), type=CSR_REAL)
    mat11 = block(mommat,1,1)
    mat12 = block(mommat,1,2)
    mat21 = block(mommat,2,1)
    mat22 = block(mommat,2,2)

    ewrite(3,*)(size(mat11))
    ewrite(3,*)(size(mat12))
    ewrite(3,*)(size(mat21))
    ewrite(3,*)(size(mat22))

    nonu = u

    nits = 0

    ewrite(2,*)("starting loop")

    nlinear_loop: do
       if(nits == mom_maxnits) exit

       nits = nits + 1
       ewrite(3,*)(nits)
       ewrite(3,*)(mom_maxnits)

       call assemble_gn_momentum_equation(mat11,mat12,mat21,mat22, &
           rhs1,rhs2, &
           m1, m2, &
           0.5*(nonu1+u1),0.5*(nonu2+u2), &
           D,mesh)

       !solve GN equation
       ksp_type = KSPGMRES
       pc_type = PCSOR

       call gallopede_block_solve(m, mommat, rhs, &
            ksp_type, pc_type, 1.0e-8, 5000)

       ewrite(3,*)(sum(m))
       ewrite(3,*)(sum(D))

       call get_vels(Mesh, D, nonu, m)

    end do nlinear_loop

    u = nonu

    deallocate( rhs )
    deallocate( nonu )
    ewrite(2,*)("deallocating mommat")
    call deallocate( mommat )
    
    ewrite(2,*)("END Subroutine solve_gn_momentum_equation")

  end subroutine solve_gn_momentum_equation

  subroutine assemble_gn_momentum_equation(mat11,mat12,mat21,mat22, &
       rhs1,rhs2,m1,m2,u1,u2, &
       D,mesh)
    implicit none
    real, intent(out), dimension(n_vels) :: rhs1,rhs2
    real, intent(in), dimension(n_vels) :: u1,u2
    real, intent(in), dimension(n_vels) :: m1,m2
    real, intent(in), dimension(N_dens) :: D
    type(csr_matrix) :: mat11, mat12, mat21, mat22
    type(dg_mesh), intent(in) :: mesh

    !locals
    integer :: ele, iloc, jloc, globi, globj
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele 
    integer, dimension(:), pointer :: u_ele_2,h_ele_2,X_ele_2
    real, dimension(mesh%nu%loc,2) :: mloc
    real, dimension(mesh%nu%loc,2,2) :: graduloc
    real, dimension(mesh%nu%loc,2) :: mloc_2
    real, dimension(mesh%nu%loc,2,2) :: graduloc_2
    real, dimension(mesh%nu%ngi,2) :: mlocgi
    real, dimension(mesh%nu%ngi,2,2) :: gradulocgi
    real, dimension(mesh%nu%ngi) :: divulocgi, u1locgi, u2locgi
    real, dimension(mesh%nh%ngi) :: dlocgi
    real, dimension(mesh%nh_f%ngi) :: dlocgi_f, u1locgi_f, u2locgi_f
    real, dimension(mesh%nu_f%ngi,2) :: mlocgi_f
    real, dimension(mesh%nu_f%ngi,2,2) :: gradulocgi_f
    real, dimension(mesh%nu_f%ngi) :: divulocgi_f

    ! List of neighbours of current element.
    integer, dimension(:), allocatable :: neigh
    integer, dimension(:), pointer :: bdyu, bdyu_2, bdyh, bdyh_2
    integer :: bdy_i, ele_2, gi, nod
    real, dimension(2,mesh%nu%loc) :: ele_X, ele_X_2
    real, dimension(2,mesh%nu_f%loc) :: ele_Xf, ele_Xf_2
    real :: kmat, kmat11, kmat12, kmat21, kmat22
    integer :: ni
    real :: detwei(mesh%nu%ngi)
    real :: detwei_2(mesh%nu%ngi), detwei_f(mesh%nu_f%ngi)
    real :: normal(2,mesh%nu_f%ngi)
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t,dnu_t_2
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t,dnh_t_2
    real, dimension(mesh%nu_f%ngi) :: udotnE, udotnF
    integer :: bcnt,i,j

    integer, dimension(mesh%nh_f%ngi), target :: surface_h_lno

    ! Local node number map for big NC element.
    integer, dimension(9) :: local_glno

    ! Local integration matrices for big NC element.
    real, dimension(2,mesh%nh%loc,mesh%nu%loc+3*mesh%nu_f%loc) :: B
    real, dimension(2,2,mesh%nu%loc+3*mesh%nu_f%loc, &
         mesh%nu%loc+3*mesh%nu_f%loc) :: BQB
    real, dimension(mesh%nu%loc,mesh%nu%loc) :: Quinv
    real, dimension(mesh%nh%loc,mesh%nh%loc) :: Qhinv
    real, dimension(2,mesh%nh_f%loc,mesh%nu_f%loc) :: Q_surfB
    real, dimension(2,mesh%nu_f%loc,mesh%nu_f%loc) :: Q_surf
    integer :: option
    real :: avg_cst

    ewrite(2,*)("subroutine assemble_gn_momentum_equation")

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
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)

       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)

       call transform_to_physical(ele_X, mesh%nu, m = mesh%nh, dn_t = dnu_t, &
            dm_t = dnh_t, detwei = detwei)

       Quinv = shape_shape(mesh%nu,mesh%nu,detwei)
       call invert(Quinv)

       !construct height on gauss points
       !ewrite(2,*)("Making dlocgi")
       Dlocgi = matmul(transpose(mesh%nh%n),D(h_ele))
       u1locgi = matmul(transpose(mesh%nu%n),u1(u_ele))
       u2locgi = matmul(transpose(mesh%nu%n),u2(u_ele))

       !construct momentum on cell nodes
       Mloc(:,1) = m1(u_ele)
       Mloc(:,2) = m2(u_ele)

       !construct momentum on gauss points

       forall(i = 1:2)
          mlocgi(:,i) = matmul(transpose(mesh%nu%n),mloc(:,i))
       end forall

       !construct gradu on cell nodes
       graduloc = get_graduloc(u1(u_ele),u2(u_ele),ele_X,mesh,Quinv)
       !construct gradu on gauss points
       forall(i = 1:2, j = 1:2)
          gradulocgi(:,i,j) = matmul(transpose(mesh%nu%n),graduloc(:,i,j))
       end forall
       divulocgi = gradulocgi(:,1,1) + gradulocgi(:,2,2)

       !ewrite(2,*)("volume integrals")

       if(.true.) then

          !volume integrals
          iloc_loop: do iloc = 1, mesh%nu%loc
             globi = u_ele(iloc)
             jloc_loop: do jloc = 1, mesh%nu%loc
                globj = u_ele(jloc)

                !==========================================================
                !Mass matrix

                kmat = sum( detwei * mesh%nu%n(iloc,:) * mesh%nu%n(jloc,:) )

                rhs1(globi) = rhs1(globi) + kmat * m1(globj)
                rhs2(globi) = rhs2(globi) + kmat * m2(globj)

                call addto(mat11,globi,globj,kmat)
                call addto(mat22,globi,globj,kmat)

                !==========================================================
                !Advection bulk integral

                if(ADVECTION_FLAG) then

                   ! div(um)
                   !SIGN CORRECT
                   kmat = -sum(detwei * mesh%nu%n(jloc,:) * &
                        ( dnu_t(iloc,:,1) * u1(globj) + &
                        dnu_t(iloc,:,2) * u2(globj) &
                        ) )

                   rhs1(globi) = rhs1(globi) - (1-theta)*dt* &
                        kmat*m1(globj)
                   rhs2(globi) = rhs2(globi) - (1-theta)*dt* &
                        kmat*m2(globj)

                   call addto(mat11,globi,globj,dt*theta*kmat)
                   call addto(mat22,globi,globj,dt*theta*kmat) !\/

                end if
             end do jloc_loop

             !============================================================
             !============================================================
             !nonlinear rhs contributions - bulk

             if(NONLINEAR_FLAG) then

                !-d_j dl/du^j_k u_{i,k}
                != - d_j ( D^3/3 div u delta_{j,k} u_{k,i} )
                != - d_j ( D^3/3 div u u_{j,i} )
                !-div dldgradu part -- SIGN CORRECT
                !minus signs because we integrated by parts and 
                !moved to the other side of the equation
                rhs1(globi)= rhs1(globi) &
                     -sum(dt*(1.0/3.0)*(dlocgi(:)**3)* &
                     detwei(:)*divulocgi(:) * &
                     (dnu_t(iloc,:,1) * gradulocgi(:,1,1) &
                     +dnu_t(iloc,:,2) * gradulocgi(:,2,1)) )
                rhs2(globi)= rhs2(globi) &
                     -sum(dt*(1.0/3.0)*(dlocgi(:)**3)* &
                     detwei(:)*divulocgi(:) * &
                     (dnu_t(iloc,:,1) * gradulocgi(:,2,1) &
                     +dnu_t(iloc,:,2) * gradulocgi(:,2,2)) )

                ! d_i ( - D^3/3 (div u)^2)
                !D-stress term -- SIGN CORRECT
                !minus signs because integrate by parts and move to other side
                !factor of 1/6-1/2 = -1/3
                rhs1(globi)=rhs1(globi)- &
                     sum(dt*(1.0/3.0)*(dlocgi(:)**3)*detwei(:)* &
                     dnu_t(iloc,:,1) * (divulocgi(:)**2) )
                rhs2(globi)=rhs2(globi)- &
                     sum(dt*(1.0/3.0)*(dlocgi(:)**3)*detwei(:)* &
                     dnu_t(iloc,:,2) * (divulocgi(:)**2) )
             end if

             !==============================================================
             !==============================================================
             !Pressure gradient bulk terms

             if(PRESSURE_FLAG) then
                ! d_i ( g_0 D^2/2)
                !pressure gradient term -- SIGN correct
                !minus sign because no integration by parts
                if(.true.) then
                   do jloc = 1, mesh%nh%loc
                      rhs1(globi) = rhs1(globi) - dt*g0 * &
                           sum( mesh%nu%n(iloc,:) * detwei & 
                           * dnh_t(jloc,:,1) ) &
                           * 0.5* D(h_ele(jloc))**2
                      rhs2(globi) = rhs2(globi) - dt*g0 * &
                           sum( mesh%nu%n(iloc,:) * detwei & 
                           * dnh_t(jloc,:,2) ) &
                           * 0.5* D(h_ele(jloc))**2
                   end do
                else
                   !plus sign because integrate by parts and move to other side
                   do jloc = 1, mesh%nh%loc
                      rhs1(globi)=rhs1(globi)+sum(dt*g0* &
                           dnu_t(iloc,:,1)*detwei(:)* 0.5 *&
                           mesh%nh%n(jloc,:)) * D(h_ele(jloc))**2
                      rhs2(globi)=rhs2(globi)+sum(dt*g0* &
                           dnu_t(iloc,:,2)*detwei(:)* 0.5 *&
                           mesh%nh%n(jloc,:)) * D(h_ele(jloc))**2
                   end do
                end if
             end if

             !==============================================================

          end do iloc_loop

       end if

       !ewrite(2,*)("surface integrals")

       !===================================================================

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
                surface_h_lno = boundary_local_num(offnods(ni,:), &
                     mesh%nh%numbering)
                bdyh => surface_h_lno
                bdyu => offnods(ni,:)
             end if
          else

             avg_cst = 0.5
             u_ele_2=>mesh%EVList_u((ELE_2-1)*mesh%nu%LOC+1:ELE_2*mesh%nu%LOC)
             h_ele_2=>mesh%EVList_h((ELE_2-1)*mesh%Nh%LOC+1:ELE_2*mesh%Nh%LOC)
             X_ele_2=>mesh%EVList_X((ELE_2-1)*mesh%Nu%LOC+1:ELE_2*mesh%Nu%LOC)

             bdy_i=ival(mesh%bdy_list, ele, ele_2)
             bdyu=> mesh%bdy_nu_lno( &
                  (bdy_i-1)*mesh%nu_f%loc+1:bdy_i*mesh%nu_f%loc)
             bdyh=>mesh%bdy_nh_lno( &
                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)

             bdy_i=ival(mesh%bdy_list, ele_2, ele)
             bdyu_2=> mesh%bdy_nu_lno( &
                  (bdy_i-1)*mesh%nu_f%loc+1:bdy_i*mesh%nu_f%loc)
             bdyh_2=>mesh%bdy_nh_lno( &
                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)

             ! Locations of local vertices.
             ele_X_2(1,:)=mesh%X(X_ele_2)
             ele_X_2(2,:)=mesh%Y(X_ele_2)

             !construct momentum on cell nodes
             Mloc_2 = get_Mloc(u1(u_ele_2),u2(u_ele_2),ele_X_2, &
                  D(h_ele_2),mesh,Quinv)

             !construct gradu on cell nodes
             graduloc_2 = get_graduloc(u1(u_ele_2),u2(u_ele_2), &
                  ele_X_2,mesh,Quinv)

          end if

          ! Locations of bdy vertices.
          ele_Xf=ele_X(:,bdyu)

          ! Change of coordinates on face
          call transform_bdy_to_physical(ele_X, ele_Xf, mesh%nu, mesh%nu_f, &
               detw_f = detwei_f,normal = normal)

          !==========================================================
          !Advection surface integral

          !ewrite(2,*)("advection integral")

          if(ADVECTION_FLAG) then

             if(ele_2.ne.0) then

                !inside face

                u1locgi_f = matmul(transpose(mesh%nu_f%n), &
                     u1(u_ele(bdyu)))
                u2locgi_f = matmul(transpose(mesh%nu_f%n), &
                     u2(u_ele(bdyu)))

                do iloc = 1, mesh%nu_f%loc
                   globi = u_ele(bdyu(iloc))
                   do jloc = 1, mesh%nu_f%loc
                      globj = u_ele(bdyu(jloc))

                      kmat = sum( mesh%nu_f%n(iloc,:) * &
                           mesh%nu_f%n(jloc,:) * detwei * &
                           (u1locgi_f * normal(1, :) + &
                           u2locgi_f * normal(2, :)) )

                      call addto(mat11,globi,globj, &
                           dt*theta*0.5*kmat)
                      call addto(mat22,globi,globj, &
                           dt*theta*0.5*kmat)
                      rhs1(globi) = rhs1(globi) - &
                           dt*(1-theta)*0.5*kmat*m1(globj)
                      rhs2(globi) = rhs2(globi) - &
                           dt*(1-theta)*0.5*kmat*m2(globj)
                   end do
                end do

                !outside face

                u1locgi_f = matmul(transpose(mesh%nu_f%n), &
                     u1(u_ele_2(bdyu_2)))
                u2locgi_f = matmul(transpose(mesh%nu_f%n), &
                     u2(u_ele_2(bdyu_2)))

                do iloc = 1, mesh%nu_f%loc
                   globi = u_ele(bdyu(iloc))
                   do jloc = 1, mesh%nu_f%loc
                      globj = u_ele_2(bdyu_2(jloc))

                      kmat = sum( mesh%nu_f%n(iloc,:) * &
                           mesh%nu_f%n(jloc,:) * detwei * &
                           (u1locgi_f * normal(1, :) + &
                           u2locgi_f * normal(2, :)) )

                      call addto(mat11,globi,globj, &
                           dt*theta*0.5*kmat)
                      call addto(mat22,globi,globj, &
                           dt*theta*0.5*kmat)
                      rhs1(globi) = rhs1(globi) - &
                           dt*(1-theta)*0.5*kmat*m1(globj)
                      rhs2(globi) = rhs2(globi) - &
                           dt*(1-theta)*0.5*kmat*m2(globj)
                   end do
                end do

             end if

          end if

          !============================================================
          !nonlinear rhs contributions - surface

          !ewrite(2,*)("nonlinear surface integrals")

          if(NONLINEAR_FLAG) then

             if(ele_2.ne.0) then

                !-d_j dl/du^j_k u_{i,k}
                != - d_j ( D^3/3 div u delta_{j,k} u_{k,i} )
                != - d_j ( D^3/3 div u u_{j,i} )
                !-div dldgradu part
                !plus signs because we integrated by parts and 
                !moved to the other side of the equation

                ! d_i ( - D^3/3 (div u)^2)
                !D-stress term -- SIGN CORRECT
                !minus signs because integrate by parts and move to other side
                !factor of 1/6-1/2 = -1/3

                !construct height on gauss points on surface
                Dlocgi_f = matmul(transpose(mesh%nh_f%n),D(h_ele(bdyh)))

                !inside face

                !construct grad u on gauss points on surface
                forall(i = 1:2, j = 1:2)
                   gradulocgi_f(:,i,j) = &
                        matmul(transpose(mesh%nu_f%n),graduloc(bdyu,i,j))
                end forall
                divulocgi_f = gradulocgi_f(:,1,1) + gradulocgi_f(:,2,2)
                !ewrite(2,*)("after divulocgi")

                do iloc = 1, mesh%nu_f%loc
                   globi = u_ele(bdyu(iloc))
                   do jloc = 1, mesh%nu_f%loc
                      globj = u_ele(bdyu(jloc))
                      rhs1(globi)= rhs1(globi) &
                           + avg_cst*sum(dt*(1.0/3.0)*(dlocgi_f(:)**3)* &
                           detwei_f(:)*divulocgi_f(:) * mesh%nu_f%n(iloc,:) * &
                           (normal(1,:) * gradulocgi_f(:,1,1) &
                           +normal(2,:) * gradulocgi_f(:,2,1)) )
                      rhs2(globi)= rhs2(globi) &
                           + avg_cst*sum(dt*(1.0/3.0)*(dlocgi_f(:)**3)* &
                           detwei_f(:)*divulocgi_f(:) * mesh%nu_f%n(iloc,:) * &
                           (normal(1,:) * gradulocgi_f(:,2,1) &
                           +normal(2,:) * gradulocgi_f(:,2,2)) )
                      !D stress
                      rhs1(globi)=rhs1(globi) + &
                           avg_cst*sum(dt*(1.0/3.0)* &
                           (dlocgi_f(:)**3)*detwei_f(:)* &
                           mesh%nu_f%n(iloc,:) *normal(1,:) *&
                           (divulocgi_f(:)**2) )
                      rhs2(globi)=rhs2(globi) + &
                           avg_cst*sum(dt*(1.0/3.0)* &
                           (dlocgi_f(:)**3)*detwei_f(:)* &
                           mesh%nu_f%n(iloc,:) *normal(2,:) * &
                           (divulocgi_f(:)**2) )
                   end do
                end do

                !outside face
                if(ele_2.ne.0) then
                   !construct grad u on gauss points on surface
                   forall(i = 1:2, j = 1:2)
                      gradulocgi_f(:,i,j) = &
                           matmul(transpose(mesh%nu_f%n),graduloc_2(bdyu_2,i,j))
                   end forall
                   divulocgi_f = gradulocgi_f(:,1,1) + gradulocgi_f(:,2,2)

                   do iloc = 1, mesh%nu_f%loc
                      globi = u_ele(bdyu(iloc))
                      do jloc = 1, mesh%nu_f%loc
                         globj = u_ele_2(bdyu_2(jloc))
                         rhs1(globi)= rhs1(globi) &
                              + avg_cst*sum(dt*(1.0/3.0)*(dlocgi_f(:)**3)* &
                              detwei_f(:)*divulocgi_f(:) * mesh%nu_f%n(iloc,:) * &
                              (normal(1,:) * gradulocgi_f(:,1,1) &
                              +normal(2,:) * gradulocgi_f(:,2,1)) )
                         rhs2(globi)= rhs2(globi) &
                              + avg_cst*sum(dt*(1.0/3.0)*(dlocgi_f(:)**3)* &
                              detwei_f(:)*divulocgi_f(:) * mesh%nu_f%n(iloc,:) * &
                              (normal(1,:) * gradulocgi_f(:,2,1) &
                              +normal(2,:) * gradulocgi_f(:,2,2)) )
                         !D stress
                         rhs1(globi)=rhs1(globi) + &
                              avg_cst*sum(dt*(1.0/3.0)* &
                              (dlocgi_f(:)**3)*detwei_f(:)* &
                              mesh%nu_f%n(iloc,:) *normal(1,:) * &
                              (divulocgi_f(:)**2) )
                         rhs2(globi)=rhs2(globi) + &
                              avg_cst*sum(dt*(1.0/3.0)* &
                              (dlocgi_f(:)**3)*detwei_f(:)* &
                              mesh%nu_f%n(iloc,:) *normal(2,:) * &
                              (divulocgi_f(:)**2) )
                      end do
                   end do
                end if
             end if
          end if

          !==============================================================
          !==============================================================

          !Pressure gradient surface term

          !ewrite(2,*)("pressure")

          if(PRESSURE_FLAG) then
             if(.false.) then
                do iloc = 1, mesh%nu_f%loc
                   globi = u_ele(bdyu(iloc))
                   do jloc = 1, mesh%nh_f%loc
                      globj = h_ele(bdyh(jloc))
                      rhs1(globi)=rhs1(globi)-sum(dt*0.5*g0* &
                           mesh%nu_f%n(iloc,:)*normal(1,:)* &
                           detwei_f(:) * mesh%nh_f%n(jloc,:))* &
                           D(globj) ** 2
                      rhs2(globi)=rhs2(globi)-sum(dt*0.5*g0* &
                           mesh%nu_f%n(iloc,:)*normal(2,:)* &
                           detwei_f(:) * mesh%nh_f%n(jloc,:))* &
                           D(globj) ** 2
                   end do
                end do
             end if
          end if

          !ewrite(2,*)("next neigh")

       end do neighbourloop2

    end do ele_loop
    deallocate( neigh )

    ewrite(2,*)("END subroutine assemble_gn_momentum_equation")

  end subroutine assemble_gn_momentum_equation

  function get_Mloc(u1,u2,ele_X,D,mesh,Qinv) result (Mloc) 
    implicit none

    type(dg_mesh), intent(in) :: mesh
    real, intent(in), dimension(2,mesh%nu%loc) :: ele_X
    real, intent(in), dimension(mesh%nu%loc) :: u1
    real, intent(in), dimension(mesh%nh%loc) :: D
    real, intent(in), dimension(mesh%nu%loc) :: u2
    real, dimension(mesh%nu%loc,2) :: Mloc
    real, intent(in), dimension(mesh%nu%loc,mesh%nu%loc) :: Qinv

    !locals
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nu%ngi) :: Dgi
    real, dimension(mesh%nu_f%ngi) :: Dgi_f
    real, dimension(mesh%nu%loc) :: rhs1,rhs2
    integer :: iloc, surf, globi, globj,jloc
    integer, dimension(:), pointer :: bdyu,bdyh
    real, dimension(mesh%nu%ngi) :: detwei
    real, dimension(mesh%nu_f%ngi) :: detwei_f
    real, dimension(2,mesh%nu_f%ngi) :: normal
    real, dimension(2,2) :: ele_xf
    integer, dimension(mesh%nh_f%ngi), target :: surface_h_lno

    rhs1 = 0.
    rhs2 = 0.

    !ewrite(2,*)("Function get_mloc")

    !volume integration

    call transform_to_physical(ele_X, mesh%nu, m = mesh%nu, dn_t = dnu_t, &
         detwei = detwei)

    !construct height at Gauss pts
    Dgi = matmul(transpose(mesh%nh%n),D)

    do iloc = 1, mesh%nu%loc
       do jloc = 1 , mesh%nu%loc
          rhs1(iloc) = rhs1(iloc) + &
               sum(mesh%nu%n(iloc,:)*mesh%nu%n(iloc,:)*dgi) * u1(jloc)
          rhs2(iloc) = rhs2(iloc) + &
               sum(mesh%nu%n(iloc,:)*mesh%nu%n(iloc,:)*dgi) * u2(jloc)
          rhs1(iloc) = rhs1(iloc) + &
               sum(dnu_t(iloc,:,1)*detwei*(dgi**3)*( &
               dnu_t(jloc,:,1) * u1(jloc) + &
               dnu_t(jloc,:,2) * u2(jloc))) /3.0
          rhs2(iloc) = rhs2(iloc) + &
               sum(dnu_t(iloc,:,2)*detwei*(dgi**3)*( &
               dnu_t(jloc,:,1) * u1(jloc) + &
               dnu_t(jloc,:,2) * u2(jloc))) /3.0
       end do
    end do

    !surface integration -- cheating because gradients constant
    do surf = 1, 3
       surface_h_lno = boundary_local_num(offnods(surf,:), &
            mesh%nh%numbering)
       bdyh => surface_h_lno
       bdyu => offnods(surf,:)
       ele_Xf = ele_X(:,bdyu)
       
       call transform_bdy_to_physical(ele_X, ele_Xf, mesh%nu, mesh%nu_f, &
            detw_f = detwei_f,normal = normal)

       !construct height at Gauss pts on surface
       Dgi_f = matmul(transpose(mesh%nh_f%n),D(bdyh))

       do iloc = 1, mesh%nu_f%loc
          do jloc = 1, mesh%nu_f%loc
             rhs1(bdyu(iloc)) = rhs1(bdyu(iloc)) - &
                  sum(normal(1,:) * detwei_f * (dgi_f**3) * &
                  mesh%nu_f%n(iloc,:) * ( &
                  dnu_t(jloc,1,1) * u1(bdyu(jloc)) + &
                  dnu_t(jloc,1,2) * u2(bdyu(jloc))) )/3.0
             rhs2(bdyu(iloc)) = rhs2(bdyu(iloc)) - &
                  sum(normal(2,:) * detwei_f * (dgi_f**3) * &
                  mesh%nu_f%n(iloc,:) * ( &
                  dnu_t(jloc,1,1) * u1(bdyu(jloc)) + &
                  dnu_t(jloc,1,2) * u2(bdyu(jloc))) )/3.0
          end do
       end do
    end do

    mloc(:,1) = matmul(Qinv,rhs1)
    mloc(:,2) = matmul(Qinv,rhs2)

    !ewrite(2,*)("END Function get_mloc")

  end function get_Mloc

  function get_graduloc(u1,u2,ele_X,mesh,Qinv) result (graduloc) 
    implicit none

    type(dg_mesh), intent(in) :: mesh
    real, intent(in), dimension(2,3) :: ele_X
    real, intent(in), dimension(mesh%nu%loc) :: u1
    real, intent(in), dimension(mesh%nu%loc) :: u2
    real, dimension(mesh%nu%loc,2,2) :: graduloc
    real, intent(in), dimension(mesh%nu%loc,mesh%nu%loc) :: Qinv

    !locals
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nu%loc) :: rhs11,rhs12,rhs21,rhs22
    integer :: iloc, surf, globi, globj,jloc
    integer, dimension(:), pointer :: bdyu,bdyh
    real, dimension(mesh%nu%ngi) :: detwei
    real, dimension(mesh%nu_f%ngi) :: detwei_f
    real, dimension(2,mesh%nu_f%ngi) :: normal
    real, dimension(2,2) :: ele_xf
    integer, dimension(mesh%nh_f%ngi), target :: surface_h_lno

    rhs11 = 0.
    rhs12 = 0.
    rhs21 = 0.
    rhs22 = 0.

    !volume integration

    call transform_to_physical(ele_X, mesh%nu, m = mesh%nu, dn_t = dnu_t, &
         detwei = detwei)

    do iloc = 1, mesh%nu%loc
       do jloc = 1 , mesh%nu%loc
          rhs11(iloc) = rhs11(iloc) + &
               sum(dnu_t(jloc,:,1)*detwei* &
               mesh%nu%n(iloc,:)) * u1(jloc)
          rhs12(iloc) = rhs12(iloc) + &
               sum(dnu_t(jloc,:,2)*detwei* &
               mesh%nu%n(iloc,:)) * u1(jloc)
          rhs21(iloc) = rhs21(iloc) + &
               sum(dnu_t(jloc,:,1)*detwei* &
               mesh%nu%n(iloc,:)) * u2(jloc)
          rhs22(iloc) = rhs22(iloc) + &
               sum(dnu_t(jloc,:,2)*detwei* &
               mesh%nu%n(iloc,:)) * u2(jloc)
       end do
    end do

    graduloc(:,1,1) = matmul(Qinv,rhs11)
    graduloc(:,1,2) = matmul(Qinv,rhs12)
    graduloc(:,2,1) = matmul(Qinv,rhs21)
    graduloc(:,2,2) = matmul(Qinv,rhs22)

  end function get_graduloc

end module gn_momentum_equation
