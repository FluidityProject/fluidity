#include "fdebug.h"
module MLCM_momentum_equation

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
!  use gallopede_solvers
  use vector_tools
  use MLCM_operator
  use vtk_io
  use fldebug
  use density_equation  
!  use gallopede_barotropic_equation
  use mesh_tools, only : adapt_timestep
  implicit none

  public :: solve_MLCM_momentum_equation, get_nonlinear_non_layer2

  private 

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"

contains

  subroutine solve_MLCM_momentum_equation(m,u,D,bottom,rho,mesh,bcs)
    implicit none

    real, dimension(n_dens), intent(in) :: bottom
    real, dimension(n_layers), intent(in) :: rho
    real, dimension(2*n_layers*n_moms), intent(inout) :: m
    real, dimension(:), intent(inout) :: u
    real, dimension(n_layers*n_dens), intent(inout) :: D
    type(dg_mesh) :: mesh
    type(bc_info), intent(in) :: bcs
!    real, dimension(n_layers), intent(inout), optional :: h0,mbar

    !local variables
    real, dimension(2*n_layers*n_vels) :: u_init
    real, dimension(2*n_layers*n_moms) :: newm, oldm, m_rk
    real, dimension(n_layers*n_dens) :: newD, D_rk
    real, dimension(n_layers*n_dens) :: oldD
    real, dimension(2*n_moms) ::rhs
    real, dimension(n_dens) :: h0, h0old
    real, dimension(2*n_moms):: mbar,mbarold, bpr
    real, dimension(n_dens*n_layers):: h_tot
    real, dimension(n_dens*n_layers) :: p
    real, dimension(n_pres*n_layers) :: np
    real, dimension(4) :: rk_mult

    integer :: nits, globi, globj, smooth, maxnits, n_min
    type(block_csr_matrix) :: bmommat
    real :: val1, dt_out, theta_pass
    integer :: layer_i,ix,jx

    KSPType :: ksp_type
    PCType :: pc_type


    assert(size(m) == 2*n_layers*n_moms)
    assert(size(D) == n_layers*n_dens)
    assert(size(bottom) == n_dens)
    assert(size(rho) == n_layers)

    ksp_type = KSPGMRES
    pc_type = PCSOR

    rk_mult= (/1.0, 0.5,  0.5, 1.0/)

    ewrite(2,*)("Subroutine solve_MLCM_momentum_equation")

    ewrite(2,*)("allocating mem")

    call allocate(bmommat,mesh%mass_h%sparsity,(/2,2/))


    if (mom_maxnits>0) then
       maxnits=mom_maxnits
    else
       maxnits=4
    end if

    ewrite(2,*)("Constructing mommat")

    if (RIGID_LID) then
       n_min=2
    else
       n_min=1
    end if


    if (BAROTROPIC_SPLIT) then

       call  get_pressure_remainder(D,u,bottom,rho,mesh,bpr,bcs)
       call  get_h0_mbar(D,u,mesh,h0,mbar)

       h0old=h0
       mbarold=mbar

       do nits=1,10
          call solve_barotropic_eqn(h0,h0old,mbar,bottom,mesh,bpr,dt/10.0)
       end do

    end if

!   initialize additional variables used in RK timestepping

    oldm = m
    newm = m
    m_rk=m
    nits = 0
    oldD=D
    newD=D
    D_rk=D
    
    
    ! Theta to zero for the first guess explicit timestep

    theta_pass=0.0

    ewrite(2,*)("starting non-linear loop")

    nits=0
    nlinear_loop: do
       if(nits == maxnits) exit

       nits = nits + 1
       ewrite(3,*) 'Nonlinear iteration', nits, 'of ', mom_maxnits

       if (DENSITY_FLAG) then

          if (mom_maxnits>0) then
             newD=oldD
          else
             newD=oldD+rk_mult(nits)*(newD-oldD)
          end if

          if (BAROTROPIC_SPLIT) then 

             do jx=1,3

             call u_bar_fix(theta*newD+(1-theta)*oldD,&
                  u,mesh,&
                  theta*h0+(1-theta)*h0old,&
                  theta*mbar+(1-theta)*mbarold)

          do layer_i = 2, N_Layers
             
             call solve_density_equation(&
                  newD((layer_i-1)*n_dens+1:layer_i*n_dens),&
                  u(2*(layer_i-1)*n_vels+1:2*(layer_i-1)*n_vels+n_vels),&
                  u(2*(layer_i-1)*n_vels+n_vels+1:&
                  2*(layer_i)*n_vels),&
                  mesh)
          end do
             
             do layer_i=1,n_dens
                newD(layer_i)=&
                     h0(layer_i)&
                     -sum(newD((/ (layer_i+(jx-1)*n_dens,jx=2,n_layers)/)))
             end do

          end do

            do layer_i = n_min, N_Layers
             
             call solve_density_equation(&
                  newD((layer_i-1)*n_dens+1:layer_i*n_dens),&
                  u(2*(layer_i-1)*n_vels+1:2*(layer_i-1)*n_vels+n_vels),&
                  u(2*(layer_i-1)*n_vels+n_vels+1:&
                  2*(layer_i)*n_vels),&
                  mesh)
          end do

          end if

          if (RIGID_LID) then

             do layer_i=1,n_dens
                newD(layer_i)=&
                     sum(oldD((/ (layer_i+(jx-1)*n_dens,jx=1,n_layers)/)))&
                     -sum(newD((/ (layer_i+(jx-1)*n_dens,jx=2,n_layers)/)))
             end do

          end if

!        call transport equation for interface heights

!          call solve_interface_equation(&
!                h_tot,&
!                u((/ ((2*(jx)*n_vels+ix, ix=1,n_vels), jx=0,n_layers-1)/)),&
!                u((/ (((2*(jx)+1)*n_vels+ix, ix=1,n_vels),&
!                jx=0,n_layers-1)/)),&
!                mesh,theta_pass,bottom)



!          do layer_i=1,n_layers-1
!             newD((layer_i-1)*N_dens+1:(layer_i)*N_dens)=&
!                  h_tot((layer_i-1)*N_dens+1:(layer_i)*N_dens)&
!                  -h_tot((layer_i)*N_dens+1:(layer_i+1)*N_dens)
!          end do
!          newD((n_layers-1)*n_dens+1:n_dens*n_layers)=&
!               h_tot((n_layers-1)*n_dens+1:n_dens*n_layers)+bottom


          ! If RK4 need D_new= D_old + f(D)

          if (mom_maxnits==0) newD=newD+(oldD-D)


         
       end if


       if (MOMENTUM_FLAG) then


          ! get the full non-hydrostatic pressure


          if (mom_maxnits==0) theta_pass=rk_mult(nits)
          if (BAROTROPIC_SPLIT) then
             call get_nonlinear_non_layer2((1-theta_pass)*oldD+(theta_pass)*D&
                  ,u,bottom,rho,mesh,p,np,h0=(theta*h0+(1-theta)*h0old))
          else
             call get_nonlinear_non_layer2((1-theta_pass)*oldD+(theta_pass)*D&
                  ,u,bottom,rho,mesh,p,np)
          end if
          if (mom_maxnits==0) theta_pass=0.0


          if (mom_maxnits>0) then
             newm=oldm
          else
             newm=oldm+rk_mult(nits)*(newm-oldm)
          end if

       do layer_i = n_min, N_Layers

          ! get matrix for momentum problem

             call assemble_MLCM_momentum_equation(bmommat, &
                  rhs(1:n_moms), rhs(n_moms+1:2*n_moms),&
                  oldm((layer_i-1)*2*n_moms+1:(layer_i-1)*2*n_moms+n_moms),&
                  oldm((layer_i-1)*2*n_moms+1+n_moms:(layer_i)*2*n_moms),&
                  p((layer_i-1)*N_dens+1:layer_i*N_dens),&
                  np((layer_i-1)*N_pres+1:layer_i*N_pres),&
               u((/ ((2*(jx)*n_vels+ix, ix=1,n_vels), jx=0,n_layers-1)/)),&
                  u((/ (((2*(jx)+1)*n_vels+ix, ix=1,n_vels),&
                           jx=0,n_layers-1)/)),&
                  (1.0-theta_pass)*oldD+(theta_pass)*D,&
                  bottom,rho,mesh,layer_i,theta_pass)

             ewrite(1,*) 'Solving Momementum equation'

             
            call gallopede_block_solve(&
                 newm(2*(layer_i-1)*n_moms+1:2*layer_i*n_moms),&
                 bmommat, rhs, &
                  ksp_type, pc_type, 1.0e-10, 20000)           




            ! If RK need m_new = m_old +f(m)
             if (mom_maxnits==0) newm=newm+(oldm-m)

             
             ewrite(2,*) 'max m', 0.5*maxval(m+newm)
             
             ewrite(2,*) 'max u', maxval(u)

          end do



       end if

       if (nits<maxnits) then


          ! check that theta gets initialized for next run through
          if (mom_maxnits>0) then
             theta_pass=theta
          else
             theta_pass=rk_mult(nits+1)
          end if

          ewrite(2,*) 'Getting new vels'
          call get_vels_layer(Mesh, bcs,&
               (1.0-theta_pass)*oldD+(theta_pass)*newD,&
               bottom, rho, u,  &
               (1-theta_pass)*oldm + theta_pass*newm,u)

       end if

       ewrite(3,*)(sum(m))
       ewrite(3,*)(sum(D))

       D=newD
       m=newm

       if (mom_maxnits==0) then
          theta_pass=0.0
          m_rk=m_rk+(1.0/(6.0*rk_mult(nits)))*(newm-oldm)
          D_rk=D_rk+(1.0/(6.0*rk_mult(nits)))*(newD-oldD)
       end if
    end do nlinear_loop

    if (mom_maxnits==0) then
       m=m_rk
       D=D_rk
    end if

    call get_vels_layer(Mesh, bcs,D,&
         bottom, rho, u,&
         m,u)

    ewrite(2,*)("deallocating mommat")


    call deallocate( bmommat )
    ewrite(2,*)("END Subroutine solve_MLCM_momentum_equation")

    call energy_find(D,u,bottom,mesh,rho)

  end subroutine solve_MLCM_momentum_equation

  subroutine assemble_MLCM_momentum_equation(mommat, &
       rhs1, rhs2, m1, m2,p,np,  &
       u1, u2, D,bottom, rho, mesh, level,theta_in)
    implicit none
    real, intent(out), dimension(n_moms) :: rhs1, rhs2
    real, intent(in), dimension(n_layers*n_vels) :: u1, u2
    real, intent(in), dimension(n_layers) :: rho
    real, intent(in), dimension(n_moms) :: m1, m2
    real, intent(in), dimension(n_dens) :: p
    real, intent(in), dimension(n_pres) :: np
    real, intent(in), dimension(n_dens*n_layers) :: D
    real, intent(in), dimension(n_dens) :: bottom
    integer, intent(in) :: level
    real, intent(in) :: theta_in
    type(block_csr_matrix), intent(inout) :: mommat
    type(dg_mesh), intent(in) :: mesh

    !locals
    integer :: ele
    integer, dimension(:), pointer :: X_ele=>null(),u_ele=>null()
    integer, dimension(:), pointer :: h_ele=>null(),m_ele=>null()
    integer, dimension(:), pointer :: n_ele=>null(), cX_ele=>null()
    integer, dimension(:), pointer :: u_ele_2=>null(),h_ele_2=>null()
    integer, dimension(:), pointer :: x_ele_2=>null(),m_ele_2=>null()
    real, dimension(2,2,mesh%nu%ngi) :: gradulocgi
    real, dimension(mesh%nh%ngi) :: pnhlocgi,dlocgi, elllocgi
    real, dimension(2,mesh%nu%ngi) :: ulocgi, du_tenslocgi
    real, dimension(2,mesh%nu%ngi) :: graddlocgi
    real, dimension(mesh%nu%ngi) :: divulocgi, xlocgi, ylocgi,djlocgi
    real, dimension(mesh%nu_f%ngi) :: u1locgi_f, u2locgi_f
    real, dimension(mesh%nu_f%ngi) ::  xlocgi_f, ylocgi_f,djlocgi_f
    real, dimension(mesh%nh_f%ngi) :: pnhlocgi_f, elllocgi_f
    real, dimension(2,mesh%nu_f%ngi) :: ulocgi_f,du_tenslocgi_f

    ! List of neighbours of current element.
    integer, dimension(:), allocatable :: neigh
    integer, dimension(:), pointer :: bdyu=>null(), bdyu_2=>null()
    integer, dimension(:), pointer :: bdyh=>null(), bdyh_2=>null()
    integer, dimension(:), pointer :: bdym=>null(), bdym_2=>null()
    integer :: bdy_i, ele_2, gi, nod
    real, dimension(2,mesh%nm%loc) :: ele_X, ele_X_2
    real, dimension(2,mesh%nm_f%loc) :: ele_Xf, ele_Xf_2
    real :: kmat, kmat11, kmat12, kmat21, kmat22
    integer :: ni , iloc ,jloc, nj
    real :: detwei(mesh%nu%ngi), dlocgi_f(mesh%nh_f%ngi)
    real :: detwei_2(mesh%nu%ngi), detwei_f(mesh%nu_f%ngi)
    real :: normal(2,mesh%nu_f%ngi)
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t,dnu_t_2
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t,dnh_t_2
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t,dnm_t_2

    integer, dimension(mesh%nh_f%loc), target :: surface_h_lno
    integer, dimension(mesh%nu_f%loc), target :: surface_u_lno
    integer, dimension(mesh%nm_f%loc), target :: surface_m_lno

    real :: avg_cst
    real, dimension(mesh%nm_f%ngi) :: upwinder

    !pointer stuff
    type(layer), dimension(N_Layers) :: Dl
    real, dimension(2,mesh%nm%loc) :: dRHS
    real, dimension(2,mesh%nm_f%loc) :: dRHS_f

    real, dimension(2,2,mesh%nm%loc,mesh%nm%loc) :: TQmm
    real, dimension(mesh%nm%loc,mesh%nm%loc) :: Qmm    
    real, dimension(2,mesh%nm%loc,mesh%nh%loc) :: Qmh    
    real, dimension(mesh%nm_f%loc,mesh%nm_f%loc) :: Qmm_f
    real, dimension(2,2,mesh%nm_f%loc,mesh%nm_f%loc) :: TQmm_f

    ewrite(1,*)("subroutine assemble_MLCM_momentum_equation")

    ewrite(2,*)("zeroing stuff")
    call zero(mommat)

    rhs1 = 0.0
    rhs2 = 0.0

    assert(size(D)  == n_layers*n_dens)
    assert(size(u1) == n_layers*n_vels)
    assert(size(u2) == n_layers*n_vels)
    assert(size(RHS1)== n_moms)
    assert(size(RHS2)== n_moms)
    assert(size(m1)== n_moms)
    assert(size(m2)== n_moms)
    assert(size(rho)== n_layers)

    print*, 'max |u|=', maxval(abs(u1((level-1)*n_vels+1:level*n_vels)))


    ! Initial (and for meshes with only one shape of element, final) size
    ! of neigh.
    ewrite(2,*)("Allocating memory for neigh");
    allocate(neigh(row_length(mesh%bdy_list,1)))

    ele_loop: do ele = 1, n_elements

       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       cX_ele=>mesh%EVList_cX((ELE-1)*3+1:ELE*3)
       n_ele=>h_ele


       assert(maxval(u_ele).le.n_vels)
       assert(maxval(m_ele).le.n_moms)
       assert(maxval(h_ele).le.n_dens)
       assert(maxval(X_ele).le.n_verts)

       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)

       call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nu,&
            dn_t = dnm_t, dm_t = dnu_t, detwei = detwei)
       call transform_to_physical(ele_X, n=mesh%nx, m = mesh%nh,&
            dm_t = dnh_t)

       !volume integrals

       Qmm=0.0
       Qmm = shape_shape(mesh%nm,mesh%nm,detwei)
       rhs1(m_ele) = rhs1(m_ele) + matmul(Qmm,m1(m_ele))
       rhs2(m_ele) = rhs2(m_ele) + matmul(Qmm,m2(m_ele))

       call addto(mommat,1,1,m_ele,m_ele,Qmm)
       call addto(mommat,2,2,m_ele,m_ele,Qmm)

       ulocgi(1,:) = matmul(u1(n_vels*(level-1)+u_ele),mesh%nu%n)
       ulocgi(2,:) = matmul(u2(n_vels*(level-1)+u_ele),mesh%nu%n)
       divulocgi=matmul(u1(n_vels*(level-1)+u_ele),dnu_t(:,:,1))+&
            matmul(u2(n_vels*(level-1)+u_ele),dnu_t(:,:,2))

       if(ADVECTION_FLAG) then

          ! div(um)

          Qmm=-dshape_dot_vector_shape(dnm_t,ulocgi,&
               mesh%nm,detwei)
!          forall(iloc=1:mesh%nm%loc,jloc=1:mesh%nm%loc)
!          Qmm(iloc,jloc)=sum(mesh%nm%n(iloc,:)*ulocgi(1,:)&
!               *dnm_t(jloc,:,1)*detwei)
!          end forall
          rhs1(m_ele)=rhs1(m_ele) - (1-theta_in)*dt*matmul(Qmm,m1(m_ele))

!          forall(iloc=1:mesh%nm%loc,jloc=1:mesh%nm%loc)
!          Qmm(iloc,jloc)=sum(mesh%nm%n(iloc,:)*ulocgi(2,:)&
!               *dnm_t(jloc,:,2)*detwei)
!          end forall

          rhs2(m_ele)=rhs2(m_ele) - (1-theta_in)*dt*matmul(Qmm,m2(m_ele))
          !  grad(u.m)

          TQmm=-dshape_outer_vector_shape(dnm_t,ulocgi,mesh%nm,&
               detwei*(-IMPLICIT_FACTOR))

          rhs1(m_ele)=rhs1(m_ele) - (1-theta_in)*dt*&
               (matmul(TQmm(1,1,:,:),m1(m_ele))+matmul(TQmm(1,2,:,:),m2(m_ele)))
          rhs2(m_ele)=rhs2(m_ele) - (1-theta_in)*dt*&
               (matmul(TQmm(2,1,:,:),m1(m_ele))+matmul(TQmm(2,2,:,:),m2(m_ele)))

          call addto(mommat,1,1,m_ele,m_ele,dt*theta_in*TQmm(1,1,:,:))
          call addto(mommat,1,2,m_ele,m_ele,dt*theta_in*TQmm(1,2,:,:))
          call addto(mommat,2,1,m_ele,m_ele,dt*theta_in*TQmm(2,1,:,:))
          call addto(mommat,2,2,m_ele,m_ele,dt*theta_in*TQmm(2,2,:,:))



          call addto(mommat,1,1,m_ele,m_ele,dt*theta_in*Qmm)
          call addto(mommat,2,2,m_ele,m_ele,dt*theta_in*Qmm)

       end if

       !============================================================

       !real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t,dnu_t_2

       if(NONLINEAR_FLAG) then

          !nonlinear term
          ! first find grad u

          gradulocgi(1,1,:) = matmul(u1((level-1)*n_vels+u_ele),dnu_t(:,:,1))
          gradulocgi(1,2,:) = matmul(u1((level-1)*n_vels+u_ele),dnu_t(:,:,2))
          gradulocgi(2,1,:) = matmul(u2((level-1)*n_vels+u_ele),dnu_t(:,:,1))
          gradulocgi(2,2,:) = matmul(u2((level-1)*n_vels+u_ele),dnu_t(:,:,2))
          dlocgi = matmul(d(h_ele),mesh%nh%n)

          TQmm(1,2,:,:) = shape_shape(mesh%nm,mesh%nm,detwei*gradulocgi(1,2,:))
          TQmm(1,1,:,:) = shape_shape(mesh%nm,mesh%nm,detwei*gradulocgi(1,1,:))
          TQmm(2,2,:,:) = shape_shape(mesh%nm,mesh%nm,detwei*gradulocgi(2,2,:))
          TQmm(2,1,:,:) = shape_shape(mesh%nm,mesh%nm,detwei*gradulocgi(2,1,:))
     

          ! need M_j d_i u_j - m_i d_j u_j

          rhs1(m_ele) = rhs1(m_ele) - (1-theta_in)*dt*&
               (matmul(TQmm(2,1,:,:),m2(m_ele)) - &
               matmul(TQmm(2,2,:,:),m1(m_ele)))
          rhs2(m_ele) = rhs2(m_ele) - (1-theta_in)*dt*&
               (matmul(TQmm(1,2,:,:),m1(m_ele)) - &
               matmul(TQmm(1,1,:,:),m2(m_ele)))
     
          call addto(mommat,1,1,m_ele,m_ele,-dt*theta_in*TQmm(2,2,:,:))
          call addto(mommat,1,2,m_ele,m_ele,dt*theta_in*TQmm(2,1,:,:))
          call addto(mommat,2,1,m_ele,m_ele,dt*theta_in*TQmm(1,2,:,:))
          call addto(mommat,2,2,m_ele,m_ele,-dt*theta_in*TQmm(1,1,:,:))

       end if

       !pressure term      

          

       if(.true.) then
          dlocgi = matmul(p(h_ele),mesh%nh%n)

          drhs = -dshape_rhs(dnm_t,dlocgi*detwei)
          Qmh=shape_dshape(mesh%nm,dnh_t,detwei)

!       rhs1(m_ele) = rhs1(m_ele) + dt*matmul(Qmh(1,:,:),p(n_ele))
!       rhs2(m_ele) = rhs2(m_ele) + dt*matmul(Qmh(2,:,:),p(n_ele))

          rhs1(m_ele) = rhs1(m_ele) + dt*drhs(1,:)
         rhs2(m_ele) = rhs2(m_ele) + dt*drhs(2,:)

         dlocgi = matmul(np(cX_ele),mesh%nm%n)

          drhs = -dshape_rhs(dnm_t,dlocgi*detwei)
          rhs1(m_ele) = rhs1(m_ele) + dt*drhs(1,:)
         rhs2(m_ele) = rhs2(m_ele) + dt*drhs(2,:)
          
       end if

       if (TIDAL_FORCING_FLAG) then

          rhs1(m_ele)=rhs1(m_ele)&
               +Tidal_mag*sin(2*3.1415927*t/Tidal_T)

       end if


       if (ROTATION_FLAG) then
          xlocgi=matmul(ele_X(1,:),mesh%nx%n)
          ylocgi=matmul(ele_X(2,:),mesh%nx%n)
           
!          drhs = dshape_rhs(dnm_t,0.5*dlocgi(:)*dlocgi(:)*rho&
!               *(f_0+beta*ylocgi)*xlocgi*ulocgi(2,:))

          drhs(1,:)=-shape_rhs(mesh%nm,rho(level)*detwei*&
               dlocgi*(f_0+beta*ylocgi)*ulocgi(2,:))
          drhs(2,:)=shape_rhs(mesh%nm,rho(level)*detwei*&
                dlocgi*(f_0+beta*ylocgi)*ulocgi(1,:))

          rhs1(m_ele) = rhs1(m_ele) - dt*drhs(1,:)
          rhs2(m_ele) = rhs2(m_ele) - dt*drhs(2,:)

          drhs=-dshape_rhs(dnm_t,0.5*rho(level)*detwei&
               *dlocgi*dlocgi*ulocgi(1,:)*b_0)

          rhs1(m_ele) = rhs1(m_ele) - dt*drhs(1,:)
          rhs2(m_ele) = rhs2(m_ele) - dt*drhs(2,:)

          drhs(1,:)=shape_rhs(mesh%nm,0.5*rho(level)&
               *detwei*b_0*dlocgi*divulocgi)
          rhs1(m_ele) = rhs1(m_ele) - dt*drhs(1,:)

          graddlocgi(1,:)=matmul(bottom(h_ele),dnh_t(:,:,1))
          graddlocgi(2,:)=matmul(bottom(h_ele),dnh_t(:,:,2))

          drhs(1,:)=-shape_rhs(mesh%nm,0.5*rho(level)*detwei*dlocgi&
               *sum(ulocgi*graddlocgi,1))
          rhs1(m_ele) = rhs1(m_ele) - dt*drhs(1,:)

          drhs=shape_vector_rhs(mesh%nm,graddlocgi,rho(level)*detwei&
               *dlocgi*ulocgi(1,:)*b_0)
          rhs1(m_ele) = rhs1(m_ele) - dt*drhs(1,:)
          rhs2(m_ele) = rhs2(m_ele) - dt*drhs(2,:)


          do nj=level+1,n_layers
             
             graddlocgi(1,:)=matmul(d(n_dens*(nj-1)+h_ele),dnh_t(:,:,1))
             graddlocgi(2,:)=matmul(d(n_dens*(nj-1)+h_ele),dnh_t(:,:,2))

             drhs(1,:)=shape_rhs(mesh%nm,0.5*rho(level)*detwei*dlocgi&
                  *sum(ulocgi*graddlocgi,1))
             rhs1(m_ele) = rhs1(m_ele) - dt*drhs(1,:)

             drhs=-shape_vector_rhs(mesh%nm,graddlocgi,rho(level)*detwei&
                  *dlocgi*ulocgi(1,:)*b_0)
             rhs1(m_ele) = rhs1(m_ele) - dt*drhs(1,:)
             rhs2(m_ele) = rhs2(m_ele) - dt*drhs(2,:)

          end do

          do nj=level+1,n_layers
             
             graddlocgi(1,:)=matmul(d(n_dens*(nj-1)+h_ele),dnh_t(:,:,1))
             graddlocgi(2,:)=matmul(d(n_dens*(nj-1)+h_ele),dnh_t(:,:,2))
             ulocgi(1,:) = matmul(u1(n_vels*(nj-1)+u_ele),mesh%nu%n)
             ulocgi(2,:) = matmul(u2(n_vels*(nj-1)+u_ele),mesh%nu%n)
             djlocgi = matmul(d(n_dens*(nj-1)+h_ele),mesh%nh%n)
             divulocgi=matmul(u1(n_vels*(nj-1)+u_ele),dnu_t(:,:,1))+&
                  matmul(u2(n_vels*(nj-1)+u_ele),dnu_t(:,:,2))

             drhs(1,:)=shape_rhs(mesh%nm,rho(level)*detwei*dlocgi&
                  *b_0*sum(ulocgi*graddlocgi,1))
             rhs1(m_ele) = rhs1(m_ele) - dt*drhs(1,:)
             drhs(1,:)=shape_rhs(mesh%nm,rho(level)*detwei*dlocgi&
                  *b_0*djlocgi*divulocgi)
             rhs1(m_ele) = rhs1(m_ele) - dt*drhs(1,:)

          end do
             

          do nj=1,level-1
             
             graddlocgi(1,:)=matmul(d(n_dens*(level-1)+h_ele),dnh_t(:,:,1))
             graddlocgi(2,:)=matmul(d(n_dens*(level-1)+h_ele),dnh_t(:,:,2))
             ulocgi(1,:) = matmul(u1(n_vels*(nj-1)+u_ele),mesh%nu%n)
             ulocgi(2,:) = matmul(u2(n_vels*(nj-1)+u_ele),mesh%nu%n)
             djlocgi = matmul(d(n_dens*(nj-1)+h_ele),mesh%nh%n)

             
             drhs=dshape_rhs(dnm_t,rho(nj)*&
                  detwei*dlocgi*djlocgi*ulocgi(1,:)*b_0)
             rhs1(m_ele) = rhs1(m_ele) - dt*drhs(1,:)
             rhs2(m_ele) = rhs2(m_ele) - dt*drhs(2,:)

             drhs=shape_vector_rhs(mesh%nm,graddlocgi,rho(nj)*&
                  detwei*djlocgi*ulocgi(1,:)*b_0)
             rhs1(m_ele) = rhs1(m_ele) - dt*drhs(1,:)
             rhs2(m_ele) = rhs2(m_ele) - dt*drhs(2,:)
          end do
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
                surface_u_lno = boundary_local_num(offnods(ni,:), &
                     mesh%nu%numbering)
                bdyu => surface_u_lno
                surface_m_lno = boundary_local_num(offnods(ni,:), &
                     mesh%nm%numbering)
                bdym => surface_m_lno
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
             bdyh=>mesh%bdy_nh_lno( &
                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)
             bdym=>mesh%bdy_nm_lno( &
                  (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)

             bdy_i=ival(mesh%bdy_list, ele_2, ele)
             bdyu_2=> mesh%bdy_nu_lno( &
                  (bdy_i-1)*mesh%nu_f%loc+1:bdy_i*mesh%nu_f%loc)
             bdyh_2=>mesh%bdy_nh_lno( &
                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)
             bdym_2=>mesh%bdy_nm_lno( &
                  (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)


             ! Locations of local vertices.
             ele_X_2(1,:)=mesh%X(X_ele_2)
             ele_X_2(2,:)=mesh%Y(X_ele_2)

           end if

          ! Locations of bdy vertices.

          ele_Xf=ele_X(:,bdym)

          ! Change of coordinates on face
          call transform_bdy_to_physical(ele_X, ele_Xf, mesh%nm, mesh%nm_f, &
               detwei_f = detwei_f,normal = normal)

          !==========================================================
          !Advection surface integral

          !ewrite(2,*)("advection integral")

          if(ADVECTION_FLAG) then

             if(ele_2.ne.0) then

                ulocgi_f(1,:) = matmul(u1(n_vels*(level-1)+u_ele(bdyu))&
                     ,mesh%nu_f%n)
                ulocgi_f(2,:) = matmul(u2(n_vels*(level-1)+u_ele(bdyu))&
                     ,mesh%nu_f%n) 

                upwinder=0.0;
                where(sum(normal*ulocgi_f,1) > 0.0) upwinder=1.0
                

                u1locgi_f = matmul(u1(n_vels*(level-1)+u_ele(bdyu)),&
                     mesh%nu_f%n)
                u2locgi_f = matmul(u2(n_vels*(level-1)+u_ele(bdyu)),&
                     mesh%nu_f%n)

                Qmm_f =shape_shape(mesh%nm_f,mesh%nm_f, &
                     upwinder*detwei_f *sum(normal*ulocgi_f,1))
                TQmm_f=shape_shape_vector_outer_vector(mesh%nm_f,&
                     mesh%nm_f, upwinder*detwei_f*(-IMPLICIT_FACTOR)&
                     ,normal,ulocgi_f)


                call addto(mommat,1,1,m_ele(bdym),m_ele(bdym),&
                     dt*theta_in*TQmm_f(1,1,:,:))
                call addto(mommat,2,2,m_ele(bdym),m_ele(bdym),&
                     dt*theta_in*TQmm_f(2,2,:,:))
                call addto(mommat,2,1,m_ele(bdym),m_ele(bdym),&
                     dt*theta_in*TQmm_f(2,1,:,:))
                call addto(mommat,1,2,m_ele(bdym),m_ele(bdym),&
                     dt*theta_in*TQmm_f(1,2,:,:))

                call addto(mommat,1,1,m_ele(bdym),m_ele(bdym),&
                     dt*theta_in*Qmm_f(:,:))
                call addto(mommat,2,2,m_ele(bdym),m_ele(bdym),&
                     dt*theta_in*Qmm_f(:,:))

                rhs1(m_ele(bdym)) = rhs1(m_ele(bdym)) - dt*(1-theta_in)*&
                     (matmul(TQmm_f(1,1,:,:),m1(m_ele(bdym)))&
                     +matmul(TQmm_f(1,2,:,:),m2(m_ele(bdym))))
                rhs2(m_ele(bdym)) = rhs2(m_ele(bdym)) - dt*(1-theta_in)*&
                     (matmul(TQmm_f(2,1,:,:),m1(m_ele(bdym)))&
                     +matmul(TQmm_f(2,2,:,:),m2(m_ele(bdym))))

                rhs1(m_ele(bdym)) = rhs1(m_ele(bdym)) - dt*(1-theta_in)*&
                     matmul(Qmm_f(:,:),m1(m_ele(bdym)))
                rhs2(m_ele(bdym)) = rhs2(m_ele(bdym)) - dt*(1-theta_in)*&
                     matmul(Qmm_f(:,:),m2(m_ele(bdym)))

                !outside face

                u1locgi_f = matmul(u1(n_vels*(level-1)+u_ele_2(bdyu_2)),&
                     mesh%nu_f%n)
                u2locgi_f = matmul(u2(n_vels*(level-1)+u_ele_2(bdyu_2)),&
                     mesh%nu_f%n)

                Qmm_f = shape_shape(mesh%nm_f,mesh%nm_f, &
                     (1.0-upwinder)*detwei_f *sum(normal*ulocgi_f,1) )
                
               TQmm_f= shape_shape_vector_outer_vector(mesh%nm_f,&
                     mesh%nm_f,  (1.0-upwinder)*detwei_f*(-IMPLICIT_FACTOR)&
                     ,normal,ulocgi_f)


                call addto(mommat,1,1,m_ele(bdym),m_ele(bdym_2),&
                     dt*theta_in*TQmm_f(1,1,:,:))
                call addto(mommat,2,2,m_ele(bdym),m_ele(bdym_2),&
                     dt*theta_in*TQmm_f(2,2,:,:))
                call addto(mommat,2,1,m_ele(bdym),m_ele_2(bdym_2),&
                     dt*theta_in*TQmm_f(2,1,:,:))
                call addto(mommat,1,2,m_ele(bdym),m_ele_2(bdym_2),&
                     dt*theta_in*TQmm_f(1,2,:,:))

               call addto(mommat,1,1,m_ele(bdym),m_ele_2(bdym_2),&
                     dt*theta_in*Qmm_f(:,:))
                call addto(mommat,2,2,m_ele(bdym),m_ele_2(bdym_2),&
                     dt*theta_in*Qmm_f(:,:))

                rhs1(m_ele(bdym)) = rhs1(m_ele(bdym)) - dt*(1-theta_in)*&
                     (matmul(TQmm_f(1,1,:,:),m1(m_ele_2(bdym_2)))&
                     +matmul(TQmm_f(1,2,:,:),m2(m_ele_2(bdym_2))))
                rhs2(m_ele(bdym)) = rhs2(m_ele(bdym)) - dt*(1-theta_in)*&
                     (matmul(TQmm_f(2,1,:,:),m1(m_ele_2(bdym_2)))&
                     +matmul(TQmm_f(2,2,:,:),m2(m_ele_2(bdym_2))))

                rhs1(m_ele(bdym)) = rhs1(m_ele(bdym)) - dt*(1-theta_in)*&
                     matmul(Qmm_f(:,:),m1(m_ele_2(bdym_2)))
                rhs2(m_ele(bdym)) = rhs2(m_ele(bdym)) - dt*(1-theta_in)*&
                     matmul(Qmm_f(:,:),m2(m_ele_2(bdym_2)))

             end if
          end if

            if (ROTATION_FLAG) then
          xlocgi_f=matmul(ele_X(1,bdym),mesh%nm_f%n)
          ylocgi_f=matmul(ele_X(2,bdym),mesh%nm_f%n)
          dlocgi_f= matmul(d(n_dens*(level-1)+h_ele(bdyh)),mesh%nh_f%n)

          
          drhs_f=shape_vector_rhs(mesh%nm_f,normal,&
               0.5*rho(level)*detwei_f&
               *dlocgi_f*dlocgi_f*ulocgi_f(1,:)*b_0)

          rhs1(m_ele(bdym)) = rhs1(m_ele(bdym)) - &
                     dt*drhs_f(1,:)
          rhs2(m_ele(bdym)) = rhs2(m_ele(bdym)) - &
               dt*drhs_f(2,:)
          

          do nj=1,level-1
             
             ulocgi_f(1,:) = matmul(u1(n_vels*(nj-1)+u_ele(bdyu))&
                  ,mesh%nu_f%n)
             ulocgi_f(2,:) = matmul(u2(n_vels*(nj-1)+u_ele(bdyu))&
                  ,mesh%nu_f%n) 
             djlocgi_f = matmul(d(n_dens*(nj-1)+h_ele(bdyh)),mesh%nh_f%n)

             
             drhs=-shape_vector_rhs(mesh%nm_f,normal,&
                  rho(nj)*detwei_f*dlocgi_f*djlocgi_f*ulocgi_f(1,:)*b_0)
                   
             rhs1(m_ele(bdym)) = rhs1(m_ele(bdym)) - &
                  dt*drhs_f(1,:)
             rhs2(m_ele(bdym)) = rhs2(m_ele(bdym)) - &
                  dt*drhs_f(2,:)



          end do

       end if

!            pressure term

             dlocgi_f=matmul(p(n_ele(bdyh)),&
                  mesh%nh_f%n)

             drhs_f=shape_vector_rhs(mesh%nm_f,normal,&
                  detwei_f*dlocgi_f)                   
             rhs1(m_ele(bdym)) = rhs1(m_ele(bdym)) + &
                  dt*drhs_f(1,:)
             rhs2(m_ele(bdym)) = rhs2(m_ele(bdym)) + &
                  dt*drhs_f(2,:)
          
             dlocgi_f=matmul(np(cX_ele(bdym)),&
                  mesh%nm_f%n)
          

             drhs_f=shape_vector_rhs(mesh%nm_f,normal,&
                  detwei_f*dlocgi_f)                   
             rhs1(m_ele(bdym)) = rhs1(m_ele(bdym)) + &
                  dt*drhs_f(1,:)
             rhs2(m_ele(bdym)) = rhs2(m_ele(bdym)) + &
                  dt*drhs_f(2,:)
          

       end do neighbourloop2


    end do ele_loop
    deallocate( neigh )

    ewrite(2,*)("END subroutine assemble_MLCM_momentum_equation")
 
  end subroutine assemble_MLCM_momentum_equation

  subroutine get_pressure_MLCM(u,D,bottom,rho,mesh,P)
    real, dimension(N_vels*N_Layers*2), intent(in) :: u
    real, dimension(N_dens*N_layers), intent(in) :: D
    real, dimension(N_verts), intent(in) :: bottom
    real, dimension(N_Layers), intent(in) :: rho
    type(dg_mesh), intent(inout) :: mesh
    real, dimension(N_dens*N_Layers), intent(out) :: P

    !locals
    real, dimension(:), allocatable :: A, B
    real, dimension(:), allocatable :: divCu, udotgradc
    integer :: i,j

    ewrite(1,*) 'subroutine get_pressure_mlcm'

    call get_hydrostatic(rho,u,bottom,D,mesh,P)


    if(NONHYDROSTATIC_PRESSURE_FLAG) then
      

       allocate( A(N_Layers*N_dens) )
       allocate( B(N_Layers*N_dens) )
       allocate( divCu(N_Layers*N_dens), udotgradc(N_Layers*N_dens) )
       
       call getA(u,mesh,A)
       call getB(u,bottom,D,mesh,B)
       call getdivCu2(D,bottom,u,mesh,divCu)
       call getudotgradc2(rho,D,bottom,u,mesh,udotgradc)


       do i = 1, N_Layers
          P((i-1)*N_dens+1:i*N_dens) = P((i-1)*N_dens+1:i*N_dens) + &
               0.5*rho(i)*( &
               D((i-1)*N_dens+1:i*N_dens) * A((i-1)*N_dens+1:i*N_dens) + &
               B((i-1)*N_dens+1:i*N_dens) &
               )**2.0 &
            - udotgradC((i-1)*N_dens+1:i*N_dens)
          if(i>1) then
             do j = 1, i-1
                P((i-1)*N_dens+1:i*N_dens) = P((i-1)*N_dens+1:i*N_dens) + &
                     rho(j)*( &
                     divCu((i-1)*N_dens+1:i*N_dens) &
                     )
             end do
          end if
       end do
    end if

    ewrite(1,*) 'end subroutine get_pressure_mlcm'

  end subroutine get_pressure_MLCM

 subroutine get_pressure_remainder(D,u,bottom,rho,mesh,pr_out,bcs)
    implicit none
    real, dimension(:),intent(in):: u
    real, dimension(:), intent(in) :: D
    real, dimension(:), intent(in) :: rho
    real, dimension(:), intent(in) :: bottom
    type(dg_mesh), intent(inout) :: mesh
    real, dimension(:) :: pr_out
     type(bc_info), intent(in) :: bcs

    !locals

    integer :: ele, i,j,k,jj
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele,n_ele,cX_ele
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(n_dens,2) :: nrhs
    real, dimension(n_dens*n_layers) :: R,S
    real, dimension(n_dens) :: h
    real, dimension(2*n_moms) :: Kf
    real, dimension(2*n_vels*n_layers) :: po,pr
    real, dimension(2*n_moms*n_layers) :: prt
    real, dimension(2,mesh%nh%loc) :: drhs, ele_X2,ele_X_2

    real, dimension(mesh%nu%loc) :: u1l,u2l
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(n_layers,mesh%nu%ngi) :: Djlocgi,Rjlocgi,Sjlocgi
    real, dimension(mesh%nu%ngi) :: detwei, lamlocgi,hlocgi
    real, dimension(2,2,mesh%nu%ngi) :: gradulocgi
    real, dimension(2,n_layers,mesh%nu%ngi) :: gradDlocgi,ulocgi
    real, dimension(2,n_layers,mesh%nu%ngi) :: gradRlocgi, gradSlocgi
    real, dimension(mesh%nu%ngi) :: udotgradClocgi, divCulocgi
    real, dimension(2,0:n_layers,mesh%nh%ngi) :: gradhlocgi
    real, dimension(mesh%nh%loc) :: Xloc
    real, dimension(mesh%nh%loc,0:n_layers) ::  Hloc
    integer, dimension(:), pointer :: u_ele_2=>null(),h_ele_2=>null()
    integer, dimension(:), pointer :: x_ele_2=>null(),m_ele_2=>null()

    ! List of neighbours of current element.
    integer, dimension(:), allocatable :: neigh
    integer, dimension(:), pointer :: bdyu=>null(), bdyu_2=>null()
    integer, dimension(:), pointer :: bdyh=>null(), bdyh_2=>null()
    integer, dimension(:), pointer :: bdym=>null(), bdym_2=>null()
    integer :: bdy_i, ele_2, ni, nod
    real, dimension(2,mesh%nm_f%loc) :: ele_Xf, ele_Xf_2
    real, dimension(mesh%nh_f%ngi) :: detwei_f, hlocgi_f
    real :: normal(2,mesh%nu_f%ngi)
    
    integer, dimension(mesh%nh_f%loc), target :: surface_h_lno
    integer, dimension(mesh%nu_f%loc), target :: surface_u_lno
    integer, dimension(mesh%nm_f%loc), target :: surface_m_lno


    real, dimension(2,mesh%nm%loc,mesh%nh%loc):: dQmh
    real, dimension(mesh%nm%loc,mesh%nm%loc):: Qmm
    real, dimension(mesh%nh%loc,mesh%nh%loc):: Qhh
    real, dimension(2,mesh%nm%loc,mesh%nm%loc):: dQmm
    KSPType :: ksp_type
    PCType  :: pc_type
    type(csr_matrix) :: mom_eqn_mat


    call allocate(mom_eqn_mat,mesh%mass_h%sparsity)
    allocate(neigh(row_length(mesh%bdy_list,1)))

    assert(size(D)==n_dens*n_layers)
    assert(size(rho)==n_layers)
    assert(size(bottom)==n_dens)

    ewrite(1,*) 'subroutine get_pressure_remainder'

    call get_R(D,u,bottom,rho,mesh,R)
    call get_S(D,u,bottom,rho,mesh,S)
    call get_Kf(D,u,bottom,rho,mesh,Kf)
    print*, 'MAXVALS=', maxval(R),maxval(S), maxval(kf)

    do i=1,n_dens
       h=sum(d((/(i+(j-1)*n_dens,j=1,n_layers)/)))-bottom(i)
    end do

    if (all(mesh%mass_h%val==0)) print*, 'Mass matrix zeroed'

    pc_type = PCSOR  ! Preconditioner context
    ksp_type = KSPGMRES ! Krylov subspace context

    do i=1,n_layers

    !zero RHS
    nrhs=0.0
    call zero(mom_eqn_mat)

    do ele = 1, N_elements

       m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       cX_ele=>mesh%EVList_cX((ELE-1)*3+1:ELE*3)
       n_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
         
       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)

       call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nh, &
            dn_t=dnm_t, dm_t = dnh_t, detwei = detwei)
       call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nu, &
            dm_t = dnu_t)

       do j=1,n_layers
          gradhlocgi(1,j,:)=-matmul(bottom(h_ele),dnh_t(:,:,1))
          gradhlocgi(2,j,:)=-matmul(bottom(h_ele),dnh_t(:,:,2))
       end do
     
       gradhlocgi(1,0,:)=matmul(h(h_ele),dnh_t(:,:,1))
       gradhlocgi(2,0,:)=matmul(h(h_ele),dnh_t(:,:,2))

       hlocgi(:)=matmul(h(h_ele),mesh%nh%n)

       do j=1,n_layers
          ulocgi(j,1,:)=matmul(u(2*(j-1)*n_vels+u_ele),mesh%nu%n)
          ulocgi(j,2,:)=matmul(u((2*(j-1)+1)*n_vels+u_ele),mesh%nu%n)
          djlocgi(j,:)=matmul(d((j-1)*n_dens+h_ele),mesh%nh%n)
          rjlocgi(j,:)=matmul(R((j-1)*n_dens+h_ele),mesh%nh%n)
          sjlocgi(j,:)=matmul(S((j-1)*n_dens+h_ele),mesh%nh%n)
          gradrlocgi(1,j,:)=matmul(r((j-1)*n_dens+h_ele),dnh_t(:,:,1))
          gradrlocgi(2,j,:)=matmul(r((j-1)*n_dens+h_ele),dnh_t(:,:,2))
          gradslocgi(1,j,:)=matmul(s((j-1)*n_dens+h_ele),dnh_t(:,:,1))
          gradslocgi(2,j,:)=matmul(s((j-1)*n_dens+h_ele),dnh_t(:,:,2))
          graddlocgi(1,j,:)=matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,1))
          graddlocgi(2,j,:)=matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,2))
          do k=j+1,n_layers
             gradhlocgi(1,k,:)=gradhlocgi(1,k,:)+&
                matmul(d((k-1)*n_dens+h_ele),dnh_t(:,:,1))
             gradhlocgi(2,k,:)=gradhlocgi(2,k,:)+&
               matmul(d((k-1)*n_dens+h_ele),dnh_t(:,:,2))
          end do
       end do

         nrhs(h_ele,1)=nrhs(h_ele,1)-(shape_rhs(mesh%nh,&
            g0*rho(i)*gradhlocgi(1,0,:)*detwei))
          nrhs(h_ele,2)=nrhs(h_ele,2)-(shape_rhs(mesh%nh,&
            g0*rho(i)*gradhlocgi(2,0,:)*detwei))

!       nrhs(m_ele,:)=nrhs(m_ele,:)+transpose(dshape_rhs(dnm_t,&
!            g0*rho(i)*hlocgi*detwei))

!      nrhs(m_ele,:)=nrhs(m_ele,:)-transpose(shape_vector_rhs(mesh%nm,&
!           graddlocgi(:,i,:),&
!           2.0*rho(i)*rjlocgi(i,:)*detwei))


!       nrhs(m_ele,:)=nrhs(m_ele,:)-transpose(shape_vector_rhs(mesh%nm,&
!            gradrlocgi(:,j,:),&
!            rho(i)*djlocgi(i,:)*detwei))

!       nrhs(m_ele,:)=nrhs(m_ele,:)&
!            -transpose(shape_vector_rhs(mesh%nm,gradhlocgi(:,i,:),&
!            rho(i)*sjlocgi(i,:)*detwei))

       do j=1,i-1

          nrhs(h_ele,1)=nrhs(h_ele,1)-(shape_rhs(mesh%nh,&
            g0*(rho(j)-rho(i))*graddlocgi(1,j,:)*detwei))
          nrhs(h_ele,2)=nrhs(h_ele,2)+(shape_rhs(mesh%nh,&
            g0*(rho(j)-rho(i))*graddlocgi(2,j,:)*detwei))
!           nrhs(m_ele,:)=nrhs(m_ele,:)+transpose(dshape_rhs(dnm_t,&
!            djlocgi(j,:)*(rho(j)-rho(i))*g0*detwei))

!          nrhs(m_ele,:)=nrhs(m_ele,:)-transpose(shape_vector_rhs(mesh%nm,&
!            graddlocgi(:,j,:),&
!            sjlocgi(j,:)*detwei*rho(j)))

!           nrhs(m_ele,:)=nrhs(m_ele,:)-transpose(shape_vector_rhs(mesh%nm,&
!            gradslocgi(:,j,:),&
!            djlocgi(j,:)*detwei*rho(j)))
         
           
        end do

        Qhh = shape_shape(mesh%nh,mesh%nh,detwei)
        call addto(mom_eqn_mat,h_ele,h_ele,Qhh)

   !ewrite(2,*)("surface integrals")

       !===================================================================

       if (size(neigh)/=row_length(mesh%bdy_list,ele)) then
          ewrite(2,*)("reallocating neigh")
          deallocate(neigh)
          allocate(neigh(row_length(mesh%bdy_list,ele)))
       end if

       neigh=row_m(mesh%bdy_list,ele)

       !surface integrals
       neighbourloop4: do ni=1,size(neigh)
          ele_2=neigh(ni)

          ! check for external bdy
          if (ele_2==0) then
             if(.false.) then
                cycle neighbourloop4
             else
                surface_h_lno = boundary_local_num(offnods(ni,:), &
                     mesh%nh%numbering)
                bdyh => surface_h_lno
                surface_u_lno = boundary_local_num(offnods(ni,:), &
                     mesh%nu%numbering)
                bdyu => surface_u_lno
                surface_m_lno = boundary_local_num(offnods(ni,:), &
                     mesh%nm%numbering)
                bdym => surface_m_lno
             end if
          else

             u_ele_2=>mesh%EVList_u((ELE_2-1)*mesh%nu%LOC+1:ELE_2*mesh%nu%LOC)
             m_ele_2=>mesh%EVList_m((ELE_2-1)*mesh%nm%LOC+1:ELE_2*mesh%nm%LOC)
             h_ele_2=>mesh%EVList_h((ELE_2-1)*mesh%Nh%LOC+1:ELE_2*mesh%Nh%LOC)
             X_ele_2=>mesh%EVList_X((ELE_2-1)*3+1:ELE_2*3)



             bdy_i=ival(mesh%bdy_list, ele, ele_2)
             bdyu=> mesh%bdy_nu_lno( &
                  (bdy_i-1)*mesh%nu_f%loc+1:bdy_i*mesh%nu_f%loc)
             bdyh=>mesh%bdy_nh_lno( &
                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)
             bdym=>mesh%bdy_nm_lno( &
                  (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)

             bdy_i=ival(mesh%bdy_list, ele_2, ele)
             bdyu_2=> mesh%bdy_nu_lno( &
                  (bdy_i-1)*mesh%nu_f%loc+1:bdy_i*mesh%nu_f%loc)
             bdyh_2=>mesh%bdy_nh_lno( &
                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)
             bdym_2=>mesh%bdy_nm_lno( &
                  (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)


             ! Locations of local vertices.
             ele_X_2(1,:)=mesh%X(X_ele_2)
             ele_X_2(2,:)=mesh%Y(X_ele_2)

           end if

          ! Locations of bdy vertices.

          ele_Xf=ele_X(:,bdym)

          ! Change of coordinates on face
          call transform_bdy_to_physical(ele_X, ele_Xf, mesh%nm, mesh%nm_f, &
               detwei_f = detwei_f,normal = normal)

          hlocgi_f=matmul(h(h_ele(bdyh)),&
                  mesh%nh_f%n)

!          nrhs(m_ele(bdym),:)=nrhs(m_ele(bdym),:)&
!               -transpose(shape_vector_rhs(mesh%nm_f,normal,&
!               detwei_f*g0*rho(i)*hlocgi_f))

          do j=1,i-1
!             hlocgi_f=matmul(d(h_ele(bdyh+(j-1)*n_dens)),&
!                  mesh%nh_f%n)

!             nrhs(m_ele(bdym),:)=nrhs(m_ele(bdym),:)&
!               -transpose(shape_vector_rhs(mesh%nm_f,normal,&
!               detwei_f*g0*(rho(j)-rho(i))*hlocgi_f))
          end do

       end do neighbourloop4

     end do

   pr(2*(i-1)*N_dens+1:2*i*N_dens)=0.0

   call gallopede_solve(pr(2*(i-1)*N_dens+1:(2*(i-1)+1)*N_dens), &
                  mom_eqn_mat,NRHS(:,1),ksp_type,pc_type,1e-16,20000)
   call gallopede_solve(pr((2*(i-1)+1)*N_dens+1:2*i*N_dens), &
                  mom_eqn_mat,NRHS(:,2),ksp_type,pc_type,1e-16,20000)
end do

   po=0.0

!   call pressure_operator_solve(Mesh, D, bottom, rho, po,&
!      pr,bcs, pr)
         call get_vels_layer(Mesh, bcs, D, bottom, rho, po,pr)

    print*, 'P maxs',  maxval(pr), maxval(po)

   call py_out(mesh,po,2*n_layers,6,.True.,"po")
   call py_out(mesh,pr,2*n_layers,6,.False.,"pr")

    call get_h0_mbar(D,po,mesh,mbar=pr_out)

    print*, maxval(pr_out)

      pr_out=pr_out-kf

       print*, maxval(pr_out)

       call  deallocate(mom_eqn_mat)

  end subroutine get_pressure_remainder

  subroutine get_R(D,u,bottom,rho,mesh,R)
    implicit none
    real, dimension(:),intent(in):: u
    real, dimension(:), intent(in) :: D
    real, dimension(:), intent(in) :: rho
    real, dimension(:), intent(in) :: bottom
    type(dg_mesh), intent(inout) :: mesh
    real, dimension(:), intent(out):: R

    !locals

    integer :: ele, i,j,k,jj
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele,n_ele,cX_ele,&
                                       r_ele
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(n_dens) :: rhs
    real, dimension(2,mesh%nh%loc) :: drhs, ele_X2

    real, dimension(mesh%nu%loc) :: u1l,u2l
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(n_layers,mesh%nh%ngi) :: Djlocgi
    real, dimension(2,mesh%nh%ngi) :: gradhlocgi
     real, dimension(2,mesh%nu%ngi) :: vlocgi
    real, dimension(mesh%nu%ngi) :: detwei
    real, dimension(n_layers,2,2,mesh%nu%ngi) :: gradulocgi
    real, dimension(n_layers,2,mesh%nu%ngi) :: gradDlocgi,ulocgi
    real, dimension(mesh%nh%loc) :: Xloc
    real, dimension(mesh%nh%loc,0:n_layers) ::  Hloc
    type(element_type), pointer::nr

    real, dimension(2,mesh%nm%loc,mesh%nh%loc):: dQmh
    real, dimension(:,:), allocatable:: Qhh
    real, dimension(2,mesh%nm%loc,mesh%nm%loc):: dQmm
    logical :: Mass_make
    type(csr_matrix) :: loc_mat
 

    KSPType :: ksp_type
    PCType  :: pc_type

    assert(size(D)==n_dens*n_layers)
    assert(size(rho)==n_layers)
    assert(size(bottom)==n_dens)

    ewrite(1,*) 'subroutine get_R'

    if (size(R)==size(D)) then
       call allocate(loc_mat,mesh%mass_h%sparsity)
       nr=>mesh%nh
       allocate(Qhh(mesh%nh%loc,mesh%nh%loc))
    else
        call allocate(loc_mat,mesh%mass_u%sparsity)
        nr=> mesh%nm
        allocate(Qhh(mesh%nm%loc,mesh%nm%loc))
    end if


    call zero(loc_mat)
    pc_type = PCSOR  ! Preconditioner context
    ksp_type = KSPGMRES ! Krylov subspace context

    !zero RHS

    do i=1,n_layers

    rhs=0.0

    do ele = 1, N_elements

       m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       cX_ele=>mesh%EVList_cX((ELE-1)*3+1:ELE*3)
       n_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)

       if (size(R)==size(D)) then
          r_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
       else
          r_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       end if
         
       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)

       call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nh, &
            dn_t=dnm_t, dm_t = dnh_t, detwei = detwei)
       call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nu, &
            dm_t = dnu_t)

       gradhlocgi(1,:)=-matmul(bottom(h_ele),dnh_t(:,:,1))
       gradhlocgi(2,:)=-matmul(bottom(h_ele),dnh_t(:,:,2))

       do j=i+1,n_layers
          gradhlocgi(1,:)=gradhlocgi(1,:)&
               +matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,1))
          gradhlocgi(2,:)=gradhlocgi(2,:)+&
               matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,2))
       end do

       do j=1,n_layers
          gradulocgi(j,1,1,:)=matmul(u(2*(j-1)*n_vels+u_ele),dnu_t(:,:,1))
          gradulocgi(j,2,1,:)=matmul(u(2*(j-1)*n_vels+u_ele),dnu_t(:,:,2))
          gradulocgi(j,1,2,:)=matmul(u((2*(j-1)+1)*n_vels+u_ele),dnu_t(:,:,1))
          gradulocgi(j,2,2,:)=matmul(u((2*(j-1)+1)*n_vels+u_ele),dnu_t(:,:,2))
          ulocgi(j,1,:)=matmul(u(2*(j-1)*n_vels+u_ele),mesh%nu%n)
          ulocgi(j,2,:)=matmul(u((2*(j-1)+1)*n_vels+u_ele),mesh%nu%n)
          djlocgi(j,:)=matmul(d((j-1)*n_dens+h_ele),mesh%nh%n)
          graddlocgi(j,1,:)=matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,1))
          graddlocgi(j,2,:)=matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,2))
       end do
     


           rhs(r_ele)=rhs(r_ele)+shape_rhs(nr,&
               1.0/3.0*djlocgi(i,:)*(gradulocgi(i,1,1,:)*gradulocgi(i,1,1,:)&
               +gradulocgi(i,1,2,:)*gradulocgi(i,2,1,:)&
               +gradulocgi(i,2,1,:)*gradulocgi(i,1,2,:)&
               +gradulocgi(i,2,2,:)*gradulocgi(i,2,2,:)&
               +(gradulocgi(i,1,1,:)+gradulocgi(i,2,2,:))&
               *(gradulocgi(i,1,1,:)+gradulocgi(i,2,2,:)))&
               *detwei)

       !1/2 u .(u . grad grad h)

       rhs(r_ele)=rhs(r_ele)-0.5*shape_rhs(nr,&
            (gradulocgi(i,1,1,:)+gradulocgi(i,2,2,:))&
            *(ulocgi(i,1,:)*gradhlocgi(1,:)&
            +ulocgi(i,2,:)*gradhlocgi(2,:))*detwei)

       rhs(r_ele)=rhs(r_ele)-0.5*dshape_dot_vector_rhs(nr%dn,&
       ulocgi(i,:,:),(ulocgi(i,1,:)*gradhlocgi(1,:)&
       +ulocgi(i,2,:)*gradhlocgi(2,:))*detwei)
       
       rhs(r_ele)=rhs(r_ele)-0.5*shape_rhs(nr,&
            ((ulocgi(i,1,:)*gradulocgi(i,1,1,:)+&
            ulocgi(i,2,:)*gradulocgi(i,2,1,:))*gradhlocgi(1,:)&
            +(ulocgi(i,1,:)*gradulocgi(i,1,2,:)+&
            ulocgi(i,2,:)*gradulocgi(i,2,2,:))*gradhlocgi(2,:))&
            *detwei)

       do j=1,n_layers


          vlocgi=ulocgi(j,:,:)

            rhs(h_ele)=rhs(h_ele)-0.5*dshape_dot_vector_rhs(nr%dn,&
                 vlocgi,&
                  djlocgi(j,:)*(gradulocgi(j,1,1,:)+gradulocgi(j,2,2,:))&
               +(ulocgi(j,1,:)*graddlocgi(j,1,:)&
               +(ulocgi(j,1,:)*graddlocgi(j,1,:)))*detwei)

            vlocgi(1,:)=ulocgi(j,1,:)*gradulocgi(j,1,1,:)&
                 +ulocgi(j,2,:)*gradulocgi(j,2,1,:)
             vlocgi(2,:)=ulocgi(j,1,:)*gradulocgi(j,1,2,:)&
                 +ulocgi(j,2,:)*gradulocgi(j,2,2,:)

             rhs(h_ele)=rhs(h_ele)-0.5*dshape_dot_vector_rhs(nr%dn,&
                  vlocgi,djlocgi(j,:)*detwei)

          end do

          if(i==1) then
             Qhh=shape_shape(nr,nr,detwei)
                 call addto(loc_mat,r_ele,r_ele,Qhh)
          end if

       end do

       call gallopede_solve(R((i-1)*N_dens+1:i*N_dens), &
                  loc_mat,RHS,ksp_type,pc_type,1e-30,20000)

    end do

       call deallocate(loc_mat)

  end subroutine get_R



  subroutine get_S(D,u,bottom,rho,mesh,S)
    implicit none
    real, dimension(:),intent(in):: u
    real, dimension(:), intent(in) :: D
    real, dimension(:), intent(in) :: rho
    real, dimension(:), intent(in) :: bottom
    type(dg_mesh), intent(inout) :: mesh
    type(element_type), pointer :: ns
    real, dimension(:) :: s

    !locals

    integer :: ele, i,j,k,jj
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele,n_ele,&
         s_ele
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(n_dens) :: rhs
    real, dimension(2,mesh%nh%loc) :: drhs, ele_X2

    real, dimension(mesh%nu%loc) :: u1l,u2l
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(n_layers,mesh%nu%ngi) :: Djlocgi
    real, dimension(2,mesh%nh%ngi) :: gradhlocgi
    real, dimension(2,mesh%nu%ngi) :: vlocgi
    real, dimension(mesh%nu%ngi) :: detwei
    real, dimension(n_layers,2,2,mesh%nu%ngi) :: gradulocgi
    real, dimension(n_layers,2,mesh%nu%ngi) :: gradDlocgi,ulocgi
    real, dimension(mesh%nh%loc) :: Xloc
    real, dimension(mesh%nh%loc,0:n_layers) ::  Hloc

    real, dimension(2,mesh%nm%loc,mesh%nh%loc):: dQmh
    real, allocatable, dimension(:,:):: Qmm
    real, dimension(2,mesh%nm%loc,mesh%nm%loc):: dQmm
    KSPType :: ksp_type
    PCType  :: pc_type
    type(csr_matrix) :: loc_mat

    assert(size(D)==n_dens*n_layers)
    assert(size(rho)==n_layers)
    assert(size(bottom)==n_dens)

    ewrite(1,*) 'subroutine get_S'

!    call zero (mesh%sparse_m)
!    call zero (mesh%mass_cX)

    if (all(mesh%mass_h%val==0)) print*, 'Mass matrix zeroed'

    if (size(S)==size(D)) then
       call allocate(loc_mat,mesh%mass_h%sparsity)
       ns=>mesh%nh
       allocate(Qmm(mesh%nh%loc,mesh%nh%loc))
    else
        call allocate(loc_mat,mesh%mass_u%sparsity)
        ns=> mesh%nm
        allocate(Qmm(mesh%nm%loc,mesh%nm%loc))
    end if
    call zero(loc_mat)

    pc_type = PCSOR  ! Preconditioner context
    ksp_type = KSPGMRES ! Krylov subspace context

    do i=1,n_layers

    !zero RHS
    rhs=0.0



    do ele = 1, N_elements

       m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
         
       
       if (size(S)==size(D)) then
          s_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
       else
          s_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       end if

       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)

       call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nh, &
            dn_t=dnm_t, dm_t = dnh_t, detwei = detwei)
       call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nu, &
            dm_t = dnu_t)

       ! 1/32D [ tr (grad u ^T . grad u) + div u *div u]

       gradhlocgi(1,:)=-matmul(bottom(h_ele),dnh_t(:,:,1))
       gradhlocgi(2,:)=-matmul(bottom(h_ele),dnh_t(:,:,2))

       do j=i+1,n_layers
          gradhlocgi(1,:)=gradhlocgi(1,:)&
               +matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,1))
          gradhlocgi(2,:)=gradhlocgi(2,:)+&
               matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,2))
       end do

       do j=1,n_layers
          gradulocgi(j,1,1,:)=matmul(u(2*(j-1)*n_vels+u_ele),dnu_t(:,:,1))
          gradulocgi(j,2,1,:)=matmul(u(2*(j-1)*n_vels+u_ele),dnu_t(:,:,2))
          gradulocgi(j,1,2,:)=matmul(u((2*(j-1)+1)*n_vels+u_ele),dnu_t(:,:,1))
          gradulocgi(j,2,2,:)=matmul(u((2*(j-1)+1)*n_vels+u_ele),dnu_t(:,:,2))
          ulocgi(j,1,:)=matmul(u(2*(j-1)*n_vels+u_ele),mesh%nu%n)
          ulocgi(j,2,:)=matmul(u((2*(j-1)+1)*n_vels+u_ele),mesh%nu%n)
          djlocgi(j,:)=matmul(d((j-1)*n_dens+u_ele),mesh%nh%n)
          graddlocgi(j,1,:)=matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,1))
          graddlocgi(j,2,:)=matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,2))
       end do
     


       rhs(h_ele)=shape_rhs(mesh%nh,&
           1.0/2.0*djlocgi(i,:)*(gradulocgi(i,1,1,:)*gradulocgi(i,1,1,:)&
           +gradulocgi(i,1,2,:)*gradulocgi(i,2,1,:)&
           +gradulocgi(i,2,1,:)*gradulocgi(i,1,2,:)&
           +gradulocgi(i,2,2,:)*gradulocgi(i,2,2,:)&
           +(gradulocgi(i,1,1,:)+gradulocgi(i,2,2,:))&
           *(gradulocgi(i,1,1,:)+gradulocgi(i,2,2,:)))&
           *detwei)

       !u .(u . grad grad h)

       rhs(s_ele)=rhs(s_ele)-shape_rhs(ns,&
            (gradulocgi(i,1,1,:)+gradulocgi(i,2,2,:))&
            *(ulocgi(i,1,:)*gradhlocgi(1,:)&
            +ulocgi(i,2,:)*gradhlocgi(2,:))*detwei)

       rhs(h_ele)=rhs(s_ele)-dshape_dot_vector_rhs(ns%dn,&
       ulocgi(i,:,:),(ulocgi(i,1,:)*gradhlocgi(1,:)&
       +ulocgi(i,2,:)*gradhlocgi(2,:))*detwei)
       
       rhs(s_ele)=rhs(s_ele)-shape_rhs(ns,&
            ((ulocgi(i,1,:)*gradulocgi(i,1,1,:)+&
            ulocgi(i,2,:)*gradulocgi(i,2,1,:))*gradhlocgi(1,:)&
            +(ulocgi(i,1,:)*gradulocgi(i,1,2,:)+&
            ulocgi(i,2,:)*gradulocgi(i,2,2,:))*gradhlocgi(2,:))&
            *detwei)

       do j=1,n_layers

          ! div (div D_j u_j u_j )
          vlocgi=ulocgi(j,:,:)


            rhs(s_ele)=rhs(s_ele)-dshape_dot_vector_rhs(ns%dn,&
                 vlocgi,&
                  djlocgi(j,:)*(gradulocgi(j,1,1,:)+gradulocgi(j,2,2,:))&
               +(ulocgi(j,1,:)*graddlocgi(j,1,:)&
               +(ulocgi(j,1,:)*graddlocgi(j,1,:)))*detwei)

            vlocgi(1,:)=ulocgi(j,1,:)*gradulocgi(j,1,1,:)&
                 +ulocgi(j,2,:)*gradulocgi(j,2,1,:)
             vlocgi(2,:)=ulocgi(j,1,:)*gradulocgi(j,1,2,:)&
                 +ulocgi(j,2,:)*gradulocgi(j,2,2,:)

             rhs(s_ele)=rhs(s_ele)-dshape_dot_vector_rhs(ns%dn,&
                  vlocgi,djlocgi(j,:)*detwei)

          end do
          
            if(i==1) then
             Qmm=shape_shape(ns,ns,detwei)
                 call addto(loc_mat,s_ele,s_ele,Qmm)
              end if
          end do

          call gallopede_solve(s((i-1)*N_dens+1:i*N_dens), &
                  loc_mat,RHS,ksp_type,pc_type,1e-30,20000)

    end do

    call deallocate(loc_mat)

  end subroutine get_S

subroutine get_Kf(D,u,bottom,rho,mesh,Kf)
    implicit none
    real, dimension(:),intent(in):: u
    real, dimension(:), intent(in) :: D
    real, dimension(:), intent(in) :: rho
    real, dimension(:), intent(in) :: bottom
    type(dg_mesh), intent(inout) :: mesh
    real, dimension(:), intent(out) :: Kf

    !locals

    integer :: ele, i,j,k,jj
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele,n_ele,cX_ele
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(n_moms,2) :: prhs
    real, dimension(2,mesh%nh%loc) :: drhs, ele_X2

    real, dimension(mesh%nu%loc) :: u1l,u2l
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(n_layers,mesh%nu%ngi) :: Djlocgi
    real, dimension(mesh%nu%ngi) :: detwei
    real, dimension(n_layers,2,2,mesh%nu%ngi) :: gradulocgi
    real, dimension(n_layers,2,mesh%nu%ngi) :: gradDlocgi,ulocgi
    real, dimension(mesh%nh%loc) :: Xloc
    real, dimension(mesh%nh%loc,0:n_layers) ::  Hloc

    real, dimension(2,mesh%nm%loc,mesh%nh%loc):: dQmh
    real, dimension(mesh%nm%loc,mesh%nu%loc):: Quu
    real, dimension(mesh%nm%loc,mesh%nm%loc):: Qmm
    real, dimension(2,mesh%nm%loc,mesh%nm%loc):: dQmm
    type(csr_matrix) :: mom_eqn
    KSPType :: ksp_type
    PCType  :: pc_type

    call allocate(mom_eqn,mesh%sparse_m%sparsity)
    call zero(mom_eqn)
    ewrite(1,*) 'subroutine get_barotropic_pressure'

!    call zero (mesh%sparse_m)
!    call zero (mesh%mass_cX)

!    if (all(mesh%sparse_m%val==0)) print*, 'Mass matrix zeroed'

    pc_type = PCSOR  ! Preconditioner context
    ksp_type = KSPGMRES ! Krylov subspace context

    !zero RHS



    prhs=0.0

    do ele = 1, N_elements

       m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       cX_ele=>mesh%EVList_cX((ELE-1)*3+1:ELE*3)
       n_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
         
       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)

       call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nh, &
            dn_t=dnm_t, dm_t = dnh_t, detwei = detwei)
       call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nu, &
            dm_t = dnu_t)


       ! 1/3 D [ tr (grad u ^T . grad u) + div u *div u]

       do j=1,n_layers
          ulocgi(j,1,:)=matmul(u(2*(j-1)*n_vels+u_ele),mesh%nu%n)
          ulocgi(j,2,:)=matmul(u((2*(j-1)+1)*n_vels+u_ele),mesh%nu%n)
          djlocgi(j,:)=matmul(d((j-1)*n_dens+u_ele),mesh%nh%n)
          graddlocgi(j,1,:)=matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,1))
          graddlocgi(j,2,:)=matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,2))
          gradulocgi(j,1,1,:)=matmul(u(2*(j-1)*n_vels+u_ele),dnu_t(:,:,1))
          gradulocgi(j,1,2,:)=matmul(u(2*(j-1)*n_vels+u_ele),dnu_t(:,:,2))
          gradulocgi(j,2,1,:)=matmul(u((2*(j-1)+1)*n_vels+u_ele),dnu_t(:,:,1))
          gradulocgi(j,2,2,:)=matmul(u((2*(j-1)+1)*n_vels+u_ele),dnu_t(:,:,2))
       end do
     

      do i=1,n_layers

!       Quu=dshape_dot_vector_shape(dnm_t,&
!            ulocgi(i,:,:),mesh%nu,djlocgi(i,:)*detwei)
              
       prhs(m_ele,1)=prhs(m_ele,1)&
            +shape_rhs(mesh%nm,detwei*(&
            djlocgi(i,:)*ulocgi(i,1,:)*&
            (gradulocgi(i,1,1,:)+gradulocgi(i,2,2,:))&
            +ulocgi(i,1,:)*ulocgi(i,1,:)*graddlocgi(i,1,:)&
            +ulocgi(i,1,:)*ulocgi(i,2,:)*graddlocgi(i,2,:)&
            +djlocgi(i,:)*ulocgi(i,1,:)*gradulocgi(i,1,1,:)&
            +djlocgi(i,:)*ulocgi(i,2,:)*gradulocgi(i,1,2,:)))
!            +matmul(Quu,u(2*(i-1)*n_vels+u_ele))
        prhs(m_ele,2)=prhs(m_ele,2)&
             +shape_rhs(mesh%nm,detwei*(&
            djlocgi(i,:)*ulocgi(i,2,:)*&
            (gradulocgi(i,1,1,:)+gradulocgi(i,2,2,:))&
            +ulocgi(i,2,:)*ulocgi(i,1,:)*graddlocgi(i,1,:)&
            +ulocgi(i,2,:)*ulocgi(i,2,:)*graddlocgi(i,2,:)&
            +djlocgi(i,:)*ulocgi(i,1,:)*gradulocgi(i,2,1,:)&
            +djlocgi(i,:)*ulocgi(i,2,:)*gradulocgi(i,2,2,:)))
!           +matmul(Quu,u((2*(i-1)+1)*n_vels+u_ele))

             end do
       
             Qmm=shape_shape(mesh%nm,mesh%nm,detwei)
             call addto(mom_eqn,m_ele,m_ele,Qmm)

          end do

       call gallopede_solve(Kf(1:n_moms), &
                  mom_eqn,PRHS(:,1),ksp_type,pc_type,1e-30,20000)
        call gallopede_solve(Kf(n_moms+1:2*n_moms), &
                  mom_eqn,PRHS(:,2),ksp_type,pc_type,1e-30,20000)

      call deallocate(mom_eqn)

  end subroutine get_Kf



    subroutine get_nonlinear_non_layer2(D,u,bottom,rho,mesh,p,np,n1,n2,&
         nn1,nn2,h0)
    implicit none
    real, dimension(:),intent(in):: u
    real, dimension(:), intent(in) :: D
    real, dimension(:), intent(in) :: rho
    real, dimension(:), intent(in) :: bottom
    type(dg_mesh), intent(inout) :: mesh
    real, dimension(:), intent(out), optional :: n1,n2,nn1,nn2
    real, dimension(:), intent(in), optional :: h0

    real, dimension(:) :: p,np

    !locals

    integer :: ele, i,j,k,jj
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele,n_ele,cX_ele
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(n_moms) :: rhs1, rhs2, nrhs1,nrhs2
    real, dimension(n_pres) :: nprhs
    real, dimension(n_dens) :: hh
    real, dimension(2,mesh%nh%loc) :: drhs, ele_X2

    real, dimension(mesh%nu%loc) :: u1l,u2l
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(n_layers,mesh%nu%ngi) :: Djlocgi, Ajlocgi, Cjlocgi,Dj0locgi
    real, dimension(mesh%nu%ngi) :: detwei, lamlocgi
    real, dimension(2,2,mesh%nu%ngi) :: gradulocgi
    real, dimension(2,0:n_layers,mesh%nu%ngi) :: gradhlocgi
    real, dimension(2,mesh%nu%ngi) :: gradlamlocgi
    real, dimension(2,n_layers,mesh%nu%ngi) :: gradDjlocgi,ujlocgi
    real, dimension(n_dens*n_layers) :: divCu, udotgradC
    real, dimension(mesh%nu%ngi) :: udotgradClocgi, divCulocgi
    real, dimension(2,mesh%nh%ngi) :: gradd2locgi
    real, dimension(mesh%nh%loc) :: Xloc
    real, dimension(mesh%nh%loc,0:n_layers) ::  Hloc

    real, dimension(2,mesh%nm%loc,mesh%nh%loc):: dQmh
    real, dimension(mesh%nm%loc,mesh%nm%loc):: Qmm
    real, dimension(2,mesh%nm%loc,mesh%nm%loc):: dQmm
    KSPType :: ksp_type
    PCType  :: pc_type

    assert(size(D)==n_dens*n_layers)
    assert(size(rho)==n_layers)
    assert(size(bottom)==n_dens)

    ewrite(1,*) 'subroutine get_nonlinear_non_layer2'

!    call zero (mesh%sparse)
!    call zero (mesh%mass_cX)

!    if (all(mesh%sparse_m%val==0)) print*, 'Mass matrix zeroed'

    pc_type = PCSOR  ! Preconditioner context
    ksp_type = KSPGMRES ! Krylov subspace context

    if (present(h0)) then
       hh=h0
    else
       do i=1,n_dens
          hh(i)=sum(D( (/ (i+(j-1)*n_dens,j=1,n_layers) /)))
       end do
    end if

    do i=1,n_layers

       !zero RHS
       nprhs=0.0


       do ele = 1, N_elements

          m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
          u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
          cX_ele=>mesh%EVList_cX((ELE-1)*3+1:ELE*3)
          n_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
         
          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)

          call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nh, &
               dn_t=dnm_t, dm_t = dnh_t, detwei = detwei)
          call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nu, &
               dm_t = dnu_t)


     
          graddjlocgi=0.0
             

          do j=1,n_layers
             Dj0locgi(j,:)=sum(D((j-1)*n_dens+1:j*n_dens))/n_dens
             Djlocgi(j,:)=matmul(D((j-1)*n_dens+h_ele),mesh%nh%n)
             gradDjlocgi(1,j,:)=matmul(D((j-1)*n_dens+h_ele),dnh_t(:,:,1))
             gradDjlocgi(2,j,:)=matmul(D((j-1)*n_dens+h_ele),dnh_t(:,:,2))
          end do

          do jj=1,mesh%nh%loc
                Hloc(jj,0)=-bottom(h_ele(jj))+hh(h_ele(jj))
          end do

          do j=1,n_layers
             do jj=1,mesh%nh%loc
                Hloc(jj,j)=-bottom(h_ele(jj))
                do k=j+1,n_layers                     
                   Hloc(jj,j)=Hloc(jj,j)+D((k-1)*n_dens+h_ele(jj))
                end do
             end do
          end do

          do j=0,n_layers
             gradhlocgi(1,j,:)=matmul(Hloc(:,j),dnh_t(:,:,1))
             gradhlocgi(2,j,:)=matmul(Hloc(:,j),dnh_t(:,:,2))
          end do


          do j=1,n_layers
             u1l=u((2*(j-1))*n_vels+u_ele)
             u2l=u((2*(j-1)+1)*n_vels+u_ele)
             ujlocgi(1,j,:)=matmul(u1l,mesh%nu%n)
             ujlocgi(2,j,:)=matmul(u2l,mesh%nu%n)
             Ajlocgi(j,:)=matmul(u1l,dnu_t(:,:,1))+&
                  matmul(u2l,dnu_t(:,:,2))
          end do
          do j=1,n_layers
             Cjlocgi(j,:)=djlocgi(j,:)*(&
                  djlocgi(j,:)*Ajlocgi(j,:)/2.0&
                  -ujlocgi(1,j,:)*gradhlocgi(1,j,:)&
                  -ujlocgi(2,j,:)*gradhlocgi(2,j,:)&
                  )

             do k=j+1,n_layers
                Cjlocgi(j,:)=cjlocgi(j,:)+djlocgi(j,:)*(&
                     +Djlocgi(k,:)*Ajlocgi(k,:)&
                     +ujlocgi(1,k,:)*graddjlocgi(1,k,:)&
                     +ujlocgi(2,k,:)*graddjlocgi(2,k,:))
             end do

          end do

          if(.true.) then

             if (.true.) then

                ! 0.5*(DA+B)**2=0.5(-u_i.grad h_i+sum_j=i^N div D_j u_j ) **2             
                             
                nprhs(cX_ele)=nprhs(cX_ele)&
                     +shape_rhs(mesh%nm,0.5*rho(i)*detwei*(&
                     -ujlocgi(1,i,:)*gradhlocgi(1,i-1,:)&
                     -ujlocgi(2,i,:)*gradhlocgi(2,i-1,:)&
                     +sum(ujlocgi(1,i:n_layers,:)*graddjlocgi(1,i:n_layers,:)&
                     +ujlocgi(2,i:n_layers,:)*graddjlocgi(2,i:n_layers,:),1)&
                     +sum(Djlocgi(i:n_layers,:)*Ajlocgi(i:n_layers,:),1)&
                     )**2)

             end if

             if (.true.) then

                !   sum_j  rho_j div ((u_j-u_i)C_j)+C_jdiv u_i
                
                do j=1,i-1

                   nprhs(cX_ele)=nprhs(cX_ele)&
                        -dshape_dot_vector_rhs(dnm_t,&
                         ujlocgi(:,j,:)- (1.0-IMPLICIT_FACTOR)*ujlocgi(:,i,:),&
                         rho(j)*detwei*cjlocgi(j,:))&
                        +shape_rhs(mesh%nm,detwei*rho(j)*Ajlocgi(i,:)&
                         *cjlocgi(j,:)*(1.0-IMPLICIT_FACTOR))
                        
                end do
             end if

              ! u.grad (D^2 F)/D = (2 F u.grad D  + D u.grad F)
              ! do as wFu.grad D -DFAw -DFu.grad w

             nprhs(cX_ele)=nprhs(cX_ele)+IMPLICIT_FACTOR*(&
                  shape_rhs(mesh%nm,rho(i)*detwei*(djlocgi(i,:)*Ajlocgi(i,:)/3.0&
                  +(-sum(ujlocgi(:,i,:)*gradhlocgi(:,i,:),1)&
                  +sum(sum(ujlocgi(:,i+1:n_layers,:)&
                  *graddjlocgi(:,i+1:n_layers,:),2),1)&
                  +sum(Djlocgi(i+1:n_layers,:)*Ajlocgi(i+1:n_layers,:),1))/2.0)&
                  *(sum(ujlocgi(:,i,:)*graddjlocgi(:,i,:),1)&
                  -djlocgi(i,:)*Ajlocgi(i,:)))&
                  -dshape_dot_vector_rhs(dnm_t,ujlocgi(:,i,:),rho(i)*detwei&
                  *djlocgi(i,:)*(djlocgi(i,:)*Ajlocgi(i,:)/3.0&
                  +(-sum(ujlocgi(:,i,:)*gradhlocgi(:,i,:),1)&
                  +sum(sum(ujlocgi(:,i+1:n_layers,:)&
                  *graddjlocgi(:,i+1:n_layers,:),2),1)&
                  +sum(Djlocgi(i+1:n_layers,:)*Ajlocgi(i+1:n_layers,:),1))/2.0)))

               
         !   G_iu_i.grad h_i+1 can be dones straight off
             nprhs(cX_ele)=nprhs(cX_ele)+IMPLICIT_FACTOR*(&
                  shape_rhs(mesh%nm,rho(i)*detwei*&
                  (djlocgi(i,:)*Ajlocgi(i,:)/2.0&
                  +(-sum(ujlocgi(:,i,:)*gradhlocgi(:,i,:),1)&
                  +sum(sum(ujlocgi(:,i+1:n_layers,:)&
                  *graddjlocgi(:,i+1:n_layers,:),2),1)&
                  +sum(Djlocgi(i+1:n_layers,:)*Ajlocgi(i+1:n_layers,:),1)))&
                  *sum(ujlocgi(:,i,:)*gradhlocgi(:,i,:),1)))
          end if

          if(i==1) then

             Qmm=shape_shape(mesh%nm,mesh%nm,detwei)
             call addto(mesh%sparse_m,m_ele,m_ele,Qmm)
             call addto(mesh%mass_cX,cX_ele,cX_ele,Qmm)
             Qmm=dshape_dot_dshape(dnm_t,dnm_t,(50.0)**2*detwei)
             call addto(mesh%mass_cX,cX_ele,cX_ele,Qmm)
!                 print*, i, ele, sum(Qmm), sum(mesh%mass_h%val)
          end if


          if (RIGID_LID) then

             nprhs(cX_ele)=nprhs(cX_ele)&
              -(1.0/6.0)*shape_rhs(mesh%nm,&
                  rho(1)*detwei*(djlocgi(1,:)*Ajlocgi(1,:))&
                  *(djlocgi(1,:)*Ajlocgi(1,:)))
!                  -(1.0/3.0)*shape_rhs(mesh%nm,&
!                  rho(1)*detwei*(djlocgi(1,:)*Ajlocgi(1,:))&
!                  *(djlocgi(1,:)*Ajlocgi(i,:)))&
!                  +(1.0/3.0)*shape_rhs(mesh%nm,&
!              rho(1)*detwei*sum((ujlocgi(:,i,:)-ujlocgi(:,1,:))*graddjlocgi(:,1,:),1))&
!              +(1.0/3.0)*dshape_dot_vector_rhs(dnh_t,ujlocgi(:,i,:)-ujlocgi(:,1,:),&
!              rho(1)*detwei*djlocgi(1,:)*djlocgi(1,:)*Ajlocgi(1,:))

                  
          end if

       end do

       

       print*, 'prhs=', sum(nprhs)
       call gallopede_solve(np((i-1)*N_pres+1:i*N_pres), &
            mesh%mass_cX,NPRHS,ksp_type,pc_type,1e-30,20000)
       
    end do



    p=0.0

    if(.true.) then
       

       ! rho_i*0.5*u_i.u_i-rho_i g h_i -sum_j=1^i-1 rho_j gD_j

       if (RIGID_LID) then


             do j=1,n_layers 
             p((j-1)*n_dens+1:j*n_dens)=&
                  p((j-1)*n_dens+1:j*n_dens)&
                  +rho(j)*g0*bottom
             if (.true.) then
                p((j-1)*n_dens+1:j*n_dens)=&
                     p((j-1)*n_dens+1:j*n_dens)&
                     +(0.5-IMPLICIT_FACTOR)&
                     *(rho(j)*u(2*(j-1)*n_vels+1:2*(j-1)*n_vels+n_vels)**2&
                     +rho(1)*u(1:n_vels)**2&
                     -2.0*rho(1)*u(1:n_vels)*u(2*(j-1)*n_vels+1:2*(j-1)*n_vels+n_vels))&
                     +(0.5-IMPLICIT_FACTOR)&
                     *(rho(j)*u(2*(j-1)*n_vels+n_vels+1:2*j*n_vels)**2&
                      +rho(1)*u(1+n_vels:2*n_vels)**2&
                     -2.0*rho(1)*u(1+n_vels:2*n_vels)&
                          *u((2*(j-1)+1)*n_vels+1:2*j*n_vels))
             end if
             do k=2,n_layers
                p((j-1)*n_dens+1:j*n_dens)=&
                     p((j-1)*n_dens+1:j*n_dens)&
                      -g0*(rho(j)-rho(1))*D((k-1)*n_dens+1:k*n_Dens)
             end do
             do k=2,j-1
                p((j-1)*n_dens+1:j*n_dens)=&
                     p((j-1)*n_dens+1:j*n_dens)&
                     -g0*(rho(k)-rho(j))*D((k-1)*n_dens+1:k*n_Dens)
             end do
          end do

       else
          
       do j=1,n_layers 
             p((j-1)*n_dens+1:j*n_dens)=&
                  p((j-1)*n_dens+1:j*n_dens)&
                  +rho(j)*g0*bottom
             if (.true.) then
                p((j-1)*n_dens+1:j*n_dens)=&
                     p((j-1)*n_dens+1:j*n_dens)&
                     +(0.5-IMPLICIT_FACTOR)*rho(j)*&
                     (u(2*(j-1)*n_vels+1:2*(j-1)*n_vels+n_vels)**2)&
                     +(0.5-IMPLICIT_FACTOR)*rho(j)*&
                     (u(2*(j-1)*n_vels+n_vels+1:2*j*n_vels)**2)
             end if
             p((j-1)*n_dens+1:j*n_dens)=&
                     p((j-1)*n_dens+1:j*n_dens)&
                     -g0*rho(j)*hh
             do k=1,j-1
                p((j-1)*n_dens+1:j*n_dens)=&
                     p((j-1)*n_dens+1:j*n_dens)&
                     -g0*(rho(k)-rho(j))*D((k-1)*n_dens+1:k*n_Dens)
             end do
          end do
       end if

    end if



       do i=1,n_layers

          ! Normalize pressure in each layer

          p((i-1)*N_dens+1:i*N_dens)=p((i-1)*N_dens+1:i*N_dens)-&
               p((i-1)*n_dens+1)
          np((i-1)*N_pres+1:i*N_pres)=np((i-1)*N_pres+1:i*N_pres)-&
               np((i-1)*n_pres+1)
!               sum(p((i-1)*N_dens+1:i*N_dens))/n_dens

       end do

       if (.false.) then

           call zero (mesh%mass_cX)

           do i=1,n_layers

             

             nrhs2=0.0
             
             do ele = 1, N_elements
                
                m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
                u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
                h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
                X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
                cX_ele=>mesh%EVList_cX((ELE-1)*3+1:ELE*3)
                n_ele=>mesh%EVList_h((ELE-1)*mesh%nh%LOC+1:ELE*mesh%Nh%LOC)

                ele_X(1,:)=mesh%X(X_ele)
                ele_X(2,:)=mesh%Y(X_ele)
             
                call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nh, &
                     dn_t=dnm_t, dm_t = dnh_t, detwei = detwei)
                call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nu, &
                     dm_t = dnu_t)



                if (i==1) then 
                   Qmm=shape_shape(mesh%nm,mesh%nm,detwei)
                   call addto(mesh%mass_cX,cX_ele,cX_ele,Qmm)
                   Qmm=dshape_dot_dshape(dnm_t,dnm_t,(100.0)**2*detwei)
                   call addto(mesh%mass_cX,cX_ele,cX_ele,Qmm)
                end if
                
                Qmm=shape_shape(mesh%nm,mesh%nm,detwei)

                nprhs(cX_ele)=nprhs(cX_ele)+&
                     matmul(Qmm,np(cX_ele+(i-1)*n_pres))

             end do

             call gallopede_solve(np((i-1)*N_pres+1:i*N_pres), &
                  mesh%mass_cX,NPRHS,ksp_type,pc_type,1e-20,20000)
          end do
             
       end if

       if (present(n1)) then

          do i=1,n_layers

             
             rhs1=0.0
             rhs2=0.0
             nrhs1=0.0
             nrhs2=0.0
             

             ! get (diagnostic) forcings
             
             do ele = 1, N_elements
                
                m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
                u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
                h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
                X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
                cX_ele=>mesh%EVList_cX((ELE-1)*3+1:ELE*3)
                n_ele=>mesh%EVList_h((ELE-1)*mesh%nh%LOC+1:ELE*mesh%Nh%LOC)
                
             ele_X(1,:)=mesh%X(X_ele)
             ele_X(2,:)=mesh%Y(X_ele)
             
             call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nh, &
                  dn_t=dnm_t, dm_t = dnh_t, detwei = detwei)
             call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nu, &
                  dm_t = dnu_t)

             dQmh=shape_dshape(mesh%nm,dnh_t,detwei)
             rhs1(m_ele)=rhs1(m_ele)+matmul(dQmh(1,:,:),p((i-1)*n_dens+n_ele))
             rhs2(m_ele)=rhs2(m_ele)+matmul(dQmh(2,:,:),p((i-1)*n_dens+n_ele))
             dQmm=shape_dshape(mesh%nm,dnm_t,detwei)
             nrhs1(m_ele)=nrhs1(m_ele)+matmul(dQmm(1,:,:),np((i-1)*n_pres+cX_ele))
             nrhs2(m_ele)=nrhs2(m_ele)+matmul(dQmm(2,:,:),np((i-1)*n_pres+cX_ele))

          end do
         
          call gallopede_solve(n1((i-1)*N_moms+1:i*N_moms), &
               mesh%Sparse_m,rhs1,ksp_type,pc_type,1e-20,20000) 
          call gallopede_solve(n2((i-1)*N_moms+1:i*N_moms), &
               mesh%sparse_m,rhs2,ksp_type,pc_type,1e-20,20000)
          call gallopede_solve(nn1((i-1)*N_moms+1:i*N_moms), &
               mesh%sparse_m,nrhs1,ksp_type,pc_type,1e-20,20000) 
          call gallopede_solve(nn2((i-1)*N_moms+1:i*N_moms), &
               mesh%sparse_m,nrhs2,ksp_type,pc_type,1e-20,20000)

      end do

      
   end if

   if (.not. NONHYDROSTATIC_PRESSURE_FLAG) np=0.0
   if (.not. HYDROSTATIC_PRESSURE_FLAG) p=0.0


    ewrite(1,*) 'end subroutine get_nonlinear_non_layer'
            
  end subroutine get_nonlinear_non_layer2
  
  subroutine getA(u,mesh,A)
    implicit none
    
    real, dimension(N_vels*N_Layers*2), intent(in) :: u
    type(dg_mesh), intent(in) :: mesh
    real, dimension(N_dens*N_layers), intent(out) :: A
       
    !locals
    integer :: ele, i,iloc
    real, dimension(:), allocatable :: RHS
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele
    real, dimension(2,mesh%nx%loc) :: ele_X
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nu%ngi) :: divulocgi,detwei
    real, dimension(mesh%nu%ngi) :: u1locgi,u2locgi
    type(csr_matrix) :: denmat
    real, dimension(mesh%nh%loc,mesh%nh%loc) :: Qhh


    KSPType :: ksp_type
    PCType  :: pc_type

    call allocate(denmat,mesh%Mass_h%sparsity)
    call zero(denmat)


    ewrite(1,*) 'subroutine geta'
    


    pc_type = PCSOR  ! Preconditioner context
    ksp_type = KSPCG ! Krylov subspace context

    allocate( RHS(N_dens) )

    open(233,file='A.dat')

    do i = 1, N_Layers
       RHS = 0.
       do ele = 1, N_elements
          m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
          u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)

          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)

          call transform_to_physical(ele_X, n=mesh%nx, m = mesh%nu, &
               dm_t = dnu_t)
          call transform_to_physical(ele_X, n=mesh%nx, m = mesh%nh, &
               dm_t = dnh_t, detwei = detwei)

                divulocgi = matmul(transpose(dnu_t(:,:,1)), &
                     u((i-1)*2*N_vels+u_ele)) + &
                     matmul(transpose(dnu_t(:,:,2)), &
                     u(((i-1)*2+1)*N_vels+u_ele))
                u1locgi = matmul(transpose(mesh%nu%n), &
                     u((i-1)*2*N_vels+u_ele))
                u2locgi =  matmul(transpose(mesh%nu%n), &
                     u(((i-1)*2+1)*N_vels+u_ele))
                forall (iloc=1:mesh%nh%loc)
                   rhs(h_ele(iloc))=rhs(h_ele(iloc))-&
                        dot_product(dnh_t(iloc,:,1)*u1locgi&
                        +dnh_t(iloc,:,2)*u2locgi,detwei)
                end forall
                if(i==1) then
                   Qhh=shape_shape(mesh%nh,mesh%nh,detwei)
                   call addto(denmat,h_ele,h_ele,Qhh)
                end if
       end do

       ewrite(1,*) 'Getting A'
       call gallopede_solve(A((i-1)*N_dens+1:N_dens), &
            denmat,RHS,ksp_type,pc_type,1e-30,2000)


       do ele = 1, N_elements
          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          write(233,*) A((i-1)*N_dens+h_ele)
       end do

    end do

    close(233)

    deallocate( RHS )
    call deallocate(denmat)

    ewrite(1,*) 'end subroutine getA'

  end subroutine getA

  subroutine get_u_squared(u,mesh,u_squared)
    implicit none
    
    real, dimension(N_vels*N_Layers*2), intent(in) :: u
    type(dg_mesh), intent(in) :: mesh
    real, dimension(N_dens*N_layers), intent(out) :: u_squared

    !locals
    integer :: ele, i
    real, dimension(:), allocatable :: RHS
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nu%ngi) :: u_squared_locgi,detwei

    KSPType :: ksp_type
    PCType  :: pc_type

    ewrite(1,*) 'subroutine get u_squared'
    
    pc_type = PCSOR  ! Preconditioner context
    ksp_type = KSPCG ! Krylov subspace context

    allocate( RHS(N_dens) )

    do i = 1, N_Layers

       if (n_vels .ne. n_dens) then

          RHS = 0.

          do ele = 1, N_elements
             m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
             u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
             h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
             X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)

             ele_X(1,:)=mesh%X(X_ele)
             ele_X(2,:)=mesh%Y(X_ele)

             call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nu, &
                  dm_t = dnu_t, detwei = detwei)
             
             u_squared_locgi = matmul(transpose(mesh%nu%n),u(((i-1)*2)*N_vels+u_ele)**2)&
                  +matmul(transpose(mesh%nu%n),u(((i-1)*2+1)*N_vels+u_ele)**2)
             
             RHS(h_ele) = RHS(h_ele) + shape_rhs(mesh%nh,detwei*u_squared_locgi)
          end do
       
          call gallopede_solve(u_squared((i-1)*N_dens+1:i*N_dens), &
               mesh%Mass_h,RHS,ksp_type,pc_type,1e-30,2000)

       else

          u_squared((i-1)*N_dens+1:i*n_dens)=sum(&
               u(((i-1)*2)*N_vels+1:((i-1)*2)*N_vels+n_vels)**2&
               + u(((i-1)*2+1)*N_vels+1:((i-1)*2+1)*N_vels+n_vels)**2)


       end if
    end do

    deallocate( RHS )

    ewrite(1,*) 'end subroutine getA'
  end subroutine get_u_squared

  subroutine getB(u,bottom,D,mesh,B)
    implicit none
    real, dimension(N_vels*N_Layers*2), target, intent(in) :: u
    real, dimension(N_dens*N_Layers), intent(in) :: D
    type(dg_mesh), intent(in) :: mesh
    real, dimension(N_dens*N_layers), intent(out) :: B
    real, dimension(N_verts), intent(in) :: bottom

    !locals
    integer :: ele, i,j
    real, dimension(:), allocatable :: RHS
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
    real, dimensioN(N_Layers,mesh%nu%ngi) :: divulocgi, ulocgi
    real, dimension(2,mesh%nu%ngi) :: gradbotlocgi
    real, dimension(N_Layers,2,mesh%nu%ngi) :: graddlocgi
    real, dimension(mesh%nu%ngi) :: dlocgi
    real, dimension(mesh%nu%ngi) :: detwei

    KSPType :: ksp_type
    PCType  :: pc_type

    ewrite(1,*) 'subroutine getb'
    
    pc_type = PCSOR  ! Preconditioner context
    ksp_type = KSPCG ! Krylov subspace context
    open(233, file='B.dat')

    allocate( RHS(N_dens) )

    do i = 1, N_Layers
       RHS = 0.
       do ele = 1, N_elements
          m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
          u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)

          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)

          call transform_to_physical(ele_X, mesh%nm, m = mesh%nu, &
               dn_t = dnm_t, dm_t = dnu_t, detwei = detwei)
          call transform_to_physical(ele_X, mesh%nm, m = mesh%nh, &
               dm_t = dnh_t)

          if(i<N_Layers) then
             do j = 1,N_Layers
                divulocgi(j,:) = &
                     matmul(transpose(dnu_t(:,:,1)),&
                     u((j-1)*2*N_vels + u_ele)) + &
                     matmul(transpose(dnu_t(:,:,2)),&
                     u(((j-1)*2+1)*N_vels + u_ele))
                graddlocgi(j,1,:) = matmul(transpose(dnh_t(:,:,1)), &
                     D((j-1)*N_dens + h_ele))
                graddlocgi(j,2,:) = matmul(transpose(dnh_t(:,:,2)), &
                     D((j-1)*N_dens + h_ele))
             end do
          end if

          gradbotlocgi(1,:) = matmul(transpose(dnm_t(:,:,1)), &
                  bottom(X_ele))
          gradbotlocgi(2,:) = matmul(transpose(dnm_t(:,:,2)), &
                  bottom(X_ele))

          ulocgi(1,:) = matmul(transpose(mesh%nu%n),&
               u((i-1)*2*N_vels + u_ele))
          ulocgi(1,:) = matmul(transpose(mesh%nu%n),&
               u(((i-1)*2+1)*N_vels + u_ele))
          RHS(h_ele) = RHS(h_ele) + &
               shape_rhs(mesh%nh,detwei*sum(ulocgi*gradbotlocgi,1))

          if(i<N_Layers) then
             do j = i+1, N_Layers
                RHS(h_ele) = RHS(h_ele) - &
                     shape_rhs(mesh%nh,detwei*sum(ulocgi*graddlocgi(j,:,:),1))
             end do
          end if

          if(i<N_Layers) then
             do j = i+1, N_Layers
                ulocgi(1,:) = &
                     matmul(transpose(mesh%nu%n),u((j-1)*2*N_vels + u_ele))
                ulocgi(2,:) = &
                     matmul(transpose(mesh%nu%n),u(((j-1)*2+1)*N_vels + u_ele))
                dlocgi(:) = matmul(transpose(mesh%nh%n), &
                     D((j-1)*N_dens + h_ele))

                rhs(h_ele)= rhs(h_ele) -dshape_dot_vector_rhs(dnh_t,&
                     ulocgi,detwei*dlocgi)
                
             end do
          end if
       end do
       call gallopede_solve(B((i-1)*N_dens+1:N_dens), &
            mesh%Mass_h,RHS,ksp_type,pc_type,1e-30,2000)

       do ele = 1, N_elements
          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          write(233,*) B((i-1)*N_dens+h_ele)
       end do

    end do

    ewrite(1,*) 'end subroutine getb'
    close(233)

  end subroutine getB

subroutine getdivCu2(D,bottom,u,mesh,divcu)
    implicit none    
    real, dimension(N_vels*N_Layers*2), target, intent(in) :: u
    type(dg_mesh), intent(inout) :: mesh
    real, dimension(N_dens*N_layers), intent(in) :: D
    real, dimension(n_dens) :: bottom 
    real, dimension(N_dens*N_Layers), intent(out) :: divcu
    
    !locals
    integer :: ele, i,j, iloc
    real, dimension(n_dens) :: RHS
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
    real, dimension(mesh%nu%loc) :: u1,u2
    real, dimension(mesh%nh%loc) :: Cloc
    real, dimension(n_layers,2,mesh%nh%ngi) :: graddlocgi
    real, dimension(2,mesh%nh%ngi) :: gradbotlocgi
    real, dimension(mesh%nh%ngi) :: detwei
    real, dimension(n_layers,mesh%nh%ngi) :: dlocgi,divulocgi
    real, dimension(2,n_layers,mesh%nh%ngi) :: ulocgi
    real, dimension(mesh%nh%loc,mesh%nh%loc) :: Qhh
    
    KSPType :: ksp_type
    PCType  :: pc_type

    ewrite(1,*) 'subroutine get divcu'

    pc_type = PCSOR  ! Preconditioner context
    ksp_type = KSPCG ! Krylov subspace context
    call zero(mesh%mass_h)

    do i = 1, N_Layers
       RHS = 0.
       do ele = 1, N_elements
          m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
          u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)

          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)

          call transform_to_physical(ele_X, mesh%nm, m = mesh%nu, &
               dn_t = dnm_t, dm_t = dnu_t, detwei = detwei)
          call transform_to_physical(ele_X, mesh%nm, m = mesh%nh, &
               dm_t = dnh_t)

          do j=1,n_layers
             u1 = u((j-1)*2*N_vels + u_ele)
             u2 = u(((j-1)*2+1)*N_vels + u_ele)
 
             divulocgi(j,:) = matmul(u1,dnu_t(:,:,1))+&
                  matmul(u2,dnu_t(:,:,2))
             ulocgi(1,j,:)=matmul(u1,mesh%nu%n)
             ulocgi(2,j,:)=matmul(u2,mesh%nu%n)
             dlocgi(j,:) = matmul(d((j-1)*n_dens+h_ele),mesh%nh%n)
       
          end do

          if(i<N_Layers) then
             do j = 1,N_Layers
                graddlocgi(j,1,:) = &
                     matmul(D((j-1)*N_dens + h_ele),dnh_t(:,:,1))
                graddlocgi(j,2,:) = &
                     matmul(D((j-1)*N_dens + h_ele),dnh_t(:,:,2))
             end do
          end if

          gradbotlocgi(1,:) = matmul(bottom(h_ele),dnh_t(:,:,1))
          gradbotlocgi(2,:) = matmul(bottom(h_ele),dnh_t(:,:,2))

!           <w Div(D_i u_i.grad b>

          RHS(h_ele) = RHS(h_ele) - &
               dshape_dot_vector_rhs(dnh_t,ulocgi(:,i,:),&
               detwei*dlocgi(i,:)*sum(ulocgi(:,i,:)*gradbotlocgi(:,:),1)) 

          if(i<N_Layers) then
             do j = i+1, N_Layers

!              -<w Div(D_i u_i.grad \sum_j={i+1}^{N}D_j))>

                RHS(h_ele) = RHS(h_ele) + &
                     dshape_dot_vector_rhs(dnh_t,ulocgi(:,i,:),&
                     detwei*dlocgi(i,:)*&
                     sum(ulocgi(:,i,:)*graddlocgi(j,:,:),1))

!              <w, Div( D_i \sum_j={i+1}^{N}div (D_j u_j)>

                RHS(h_ele) = RHS(h_ele) - &
                     dshape_dot_vector_rhs(dnh_t,ulocgi(:,i,:),&
                     detwei*dlocgi(i,:)* &
                     divulocgi(j,:)*dlocgi(j,:))
    
                RHS(h_ele) = RHS(h_ele) - &
                     dshape_dot_vector_rhs(dnh_t,ulocgi(:,i,:),&
                     detwei*dlocgi(i,:)* &
                     sum(ulocgi(:,j,:)*graddlocgi(j,:,:),1))
             end do
          end if


!               <w ,div (D_i u_i D_i div u_i/2)>

          rhs(h_ele)=rhs(h_ele)-dshape_dot_vector_rhs(dnh_t,&
               ulocgi(:,i,:),0.5*detwei*dlocgi(i,:)*dlocgi(i,:)*divulocgi(i,:))
          
          if (i == 1) then

             Qhh= shape_shape(mesh%nh,mesh%nh,detwei)
             call addto(mesh%mass_h,h_ele,h_ele,Qhh)

             end if

       end do

       call gallopede_solve(divCu((i-1)*N_dens+1:i*N_dens), &
            mesh%Mass_h,RHS,ksp_type,pc_type,1.d-31,10000)       

    end do

  end subroutine getdivCu2

  subroutine getdivCu(D,A,B,u,mesh,divcu)
    implicit none    
    real, dimension(N_vels*N_Layers*2), target, intent(in) :: u
    type(dg_mesh), intent(in) :: mesh
    real, dimension(N_dens*N_layers), intent(in) :: B,D,A
    real, dimension(N_dens*N_Layers), intent(out) :: divcu
    
    !locals
    integer :: ele, i,j, iloc
    real, dimension(:), allocatable :: RHS
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
    real, dimension(mesh%nu%loc) :: u1,u2
    real, dimension(mesh%nh%loc) :: Cloc
    real, dimension(2,mesh%nh%ngi) :: gradclocgi
    real, dimension(mesh%nh%ngi) :: divulocgi,detwei
    real, dimension(mesh%nh%ngi) :: u1locgi,u2locgi,clocgi
    
    KSPType :: ksp_type
    PCType  :: pc_type
    
    open(233,file='divCu.dat')

    ewrite(1,*) 'subroutine get divcu'

    pc_type = PCSOR  ! Preconditioner context
    ksp_type = KSPCG ! Krylov subspace context

    allocate( RHS(N_dens) )

    do i = 1, N_Layers
       RHS = 0.
       do ele = 1, N_elements
          m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
          u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)

          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)

          call transform_to_physical(ele_X, mesh%nm, m = mesh%nu, &
               dn_t = dnm_t, dm_t = dnu_t, detwei = detwei)
          call transform_to_physical(ele_X, mesh%nm, m = mesh%nh, &
               dm_t = dnh_t)

          Cloc = D((i-1)*N_dens + h_ele) * ( &
               D((i-1)*N_dens + h_ele) * &
               A((i-1)*N_dens + h_ele) / 2.0 + &
               B((i-1)*N_dens + h_ele) )

          u1 = u((i-1)*2*N_vels + u_ele)
          u2 = u(((i-1)*2+1)*N_vels + u_ele)
          gradclocgi(1,:) = matmul(transpose(dnh_t(:,:,1)),Cloc)
          gradclocgi(2,:) = matmul(transpose(dnh_t(:,:,2)),Cloc)
          divulocgi = matmul(transpose(dnu_t(:,:,1)),u1) + &
               matmul(transpose(dnu_t(:,:,2)),u2)
          u1locgi = matmul(transpose(mesh%nu%n),u1)
          u2locgi = matmul(transpose(mesh%nu%n),u2)

          forall (iloc=1:mesh%nh%loc)
             rhs(h_ele(iloc))=rhs(h_ele(iloc))-&
                  dot_product(dnh_t(iloc,:,1)*u1locgi&
                  +dnh_t(iloc,:,2)*u2locgi,detwei*Clocgi)
          end forall

       end do

       call gallopede_solve(divCu((i-1)*N_dens+1:N_dens), &
            mesh%Mass_h,RHS,ksp_type,pc_type,1e-31,10000)       

       do ele = 1, N_elements
          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          write(233,*) divCu((i-1)*N_dens+h_ele)
       end do

    end do

    close(233)
        ewrite(1,*) 'end subroutine getdivcu'

  end subroutine getdivCu

 subroutine getudotgradc(rho,D,A,B,u,mesh,udotgradc)
    implicit none
    real, dimension(N_layers), intent(in) :: rho
    real, dimension(N_vels*N_Layers*2), target, intent(in) :: u
    type(dg_mesh), intent(in) :: mesh
    real, dimension(N_dens*N_layers), intent(in) :: D, A,B

    real, dimension(N_dens*N_Layers), intent(out) :: udotgradc
    
    !locals
    integer :: ele, i,j
    real, dimension(:), allocatable :: RHS
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
    real, dimension(N_Layers,2,mesh%nu%ngi) :: graddlocgi
    real, dimension(mesh%nu%loc) :: u1,u2
    real, dimension(mesh%nh%loc) :: Cloc
    real, dimension(2,mesh%nh%ngi) :: gradclocgi
    real, dimension(mesh%nh%ngi) :: divulocgi,detwei
    real, dimension(mesh%nh%ngi) :: u1locgi,u2locgi,clocgi
   
    
    KSPType :: ksp_type
    PCType  :: pc_type

    ewrite(1,*) 'subroutine get_udotgradc'
    
    pc_type = PCSOR  ! Preconditioner context
    ksp_type = KSPCG ! Krylov subspace context


    open(233,file='udotgradc.dat')
    allocate( RHS(N_dens) )

    do i = 1, N_Layers
       RHS = 0.
       do ele = 1, N_elements
          m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
          u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)

          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)

          u1 = u((i-1)*2*N_vels + u_ele)
          u2 = u(((i-1)*2+1)*N_vels + u_ele)

          call transform_to_physical(ele_X, mesh%nm, m = mesh%nu, &
               dn_t = dnm_t, dm_t = dnu_t, detwei = detwei)
          call transform_to_physical(ele_X, mesh%nm, m = mesh%nh, &
               dm_t = dnh_t)

          Cloc = 0.
          
          do j=1,i-1

          
          end do

          gradclocgi(1,:) = matmul(transpose(dnh_t(:,:,1)),Cloc)
          gradclocgi(2,:) = matmul(transpose(dnh_t(:,:,2)),Cloc)
          u1locgi = matmul(transpose(mesh%nu%n),u1)
          u2locgi = matmul(transpose(mesh%nu%n),u2)

          RHS(h_ele) = RHS(h_ele) + shape_rhs( &
               mesh%nh,detwei* &
               (gradclocgi(1,:)*u1locgi + gradclocgi(2,:)*u2locgi))
       end do

       call gallopede_solve(udotgradc((i-1)*N_dens+1:N_dens), &
            mesh%Mass_h,RHS,ksp_type,pc_type,1e-20,2000)       

       do ele = 1, N_elements
          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          write(233,*) udotgradc((i-1)*N_dens+h_ele)
       end do

    end do


    close(233)
    ewrite(1,*) 'end subroutine getudotgradc'

  end subroutine getudotgradc


  subroutine getudotgradc2(rho,D,bottom,u,mesh,udotgradc)
    implicit none
    real, dimension(N_layers), intent(in) :: rho
    real, dimension(N_vels*N_Layers*2), target, intent(in) :: u
    type(dg_mesh), intent(inout) :: mesh
    real, dimension(N_dens*N_layers), intent(in) :: D
    real, dimension(n_dens), intent(in) :: bottom
    real, dimension(N_dens*N_Layers), intent(out) :: udotgradc
    
    !locals
    integer :: ele, i,j,k
    real, dimension(n_dens):: RHS
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
    real, dimension(N_Layers,2,mesh%nu%ngi) :: graddlocgi
    real, dimension(mesh%nu%loc) :: u1,u2
    real, dimension(2,mesh%nh%ngi) :: gradbotlocgi
    real, dimension(mesh%nh%ngi) :: detwei
    real, dimension(n_layers,mesh%nh%ngi) :: dlocgi, divulocgi
    real, dimension(n_layers,2,mesh%nh%ngi) :: ulocgi

    real, dimension(mesh%nh%loc,mesh%nh%loc) :: Qhh

    KSPType :: ksp_type
    PCType  :: pc_type

    ewrite(1,*) 'subroutine get_udotgradc2'
    
    pc_type = PCSOR  ! Preconditioner context
    ksp_type = KSPCG ! Krylov subspace context

    call zero(mesh%mass_h)

    do i = 1, N_Layers
       RHS = 0.
       do ele = 1, N_elements
          m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
          u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)

          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)
          
          call transform_to_physical(ele_X, mesh%nm, m = mesh%nu, &
               dn_t = dnm_t, dm_t = dnu_t, detwei = detwei)
          call transform_to_physical(ele_X, mesh%nm, m = mesh%nh, &
               dm_t = dnh_t)
          
          do j=1,n_layers
             u1 = u((j-1)*2*N_vels + u_ele)
             u2 = u(((j-1)*2+1)*N_vels + u_ele)
             
                divulocgi(j,:) = matmul(u1,dnu_t(:,:,1))+&
                     matmul(u2,dnu_t(:,:,2))
                ulocgi(j,1,:)=matmul(u1,mesh%nu%n)
                ulocgi(j,2,:)=matmul(u2,mesh%nu%n)
                Dlocgi(j,:)=matmul(D((j-1)*n_dens+h_ele),mesh%nh%n)
             end do
             
             
             if(i<N_Layers) then
                do j = 1,N_Layers
                   graddlocgi(j,1,:) = matmul(transpose(dnh_t(:,:,1)), &
                        D((j-1)*N_dens + h_ele))
                   graddlocgi(j,2,:) = matmul(transpose(dnh_t(:,:,2)), &
                        D((j-1)*N_dens + h_ele))
                end do
             end if
             
             gradbotlocgi(1,:) = matmul(transpose(dnh_t(:,:,1)), &
                  bottom(h_ele))
             gradbotlocgi(2,:) = matmul(transpose(dnh_t(:,:,2)), &
                  bottom(h_ele))
             
          do j=1,i-1
             
             RHS(h_ele) = RHS(h_ele) - &
                  dshape_dot_vector_rhs(dnh_t,ulocgi(i,:,:),&
                  detwei*rho(j)*dlocgi(j,:)*ulocgi(j,1,:)*gradbotlocgi(1,:))
             RHS(h_ele) = RHS(h_ele) - &
                  dshape_dot_vector_rhs(dnh_t,ulocgi(:,i,:),&
                  detwei*rho(j)*dlocgi(j,:)*ulocgi(j,2,:)*gradbotlocgi(2,:))
             
             if(i<N_Layers) then
                do k = j+1, N_Layers
                   RHS(h_ele) = RHS(h_ele) + &
                        dshape_dot_vector_rhs(dnh_t,ulocgi(i,:,:),&
                        detwei*rho(j)*dlocgi(j,:)*&
                        ulocgi(j,2,:)*graddlocgi(k,1,:))
                   RHS(h_ele) = RHS(h_ele) + &
                        dshape_dot_vector_rhs(dnh_t,ulocgi(i,:,:),&
                        detwei*rho(j)*dlocgi(j,:)*&
                        ulocgi(j,1,:)*graddlocgi(k,2,:))
                end do
             end if
             
             if(i<N_Layers) then
                do k = j+1, N_Layers
                   RHS(h_ele) = RHS(h_ele) - &
                        dshape_dot_vector_rhs(dnh_t,ulocgi(i,:,:),&
                        detwei*rho(j)*dlocgi(j,:)* &
                        divulocgi(k,:)*dlocgi(k,:))
                   RHS(h_ele) = RHS(h_ele) - &
                        dshape_dot_vector_rhs(dnh_t,ulocgi(i,:,:),&
                        detwei*rho(j)*dlocgi(j,:)* &
                        (ulocgi(k,1,:)*graddlocgi(k,1,:)&
                        +ulocgi(k,2,:)*graddlocgi(k,2,:)))
                end do
             end if
             
             
             rhs(h_ele)=rhs(h_ele)-dshape_dot_vector_rhs(dnh_t,&
                  ulocgi(i,:,:),0.5*detwei*rho(j)*&
                  (dlocgi(j,:)**2.0)*divulocgi(j,:))
             
             rhs(h_ele)=rhs(h_ele)-shape_rhs(mesh%nh,detwei*rho(j)*&
                  sum(ulocgi(j,:,:)*gradbotlocgi,1)*divulocgi(i,:))
             
             do k=j+1,n_layers
                rhs(h_ele)=rhs(h_ele)+shape_rhs(mesh%nh,detwei*rho(j)*&
                     sum(ulocgi(j,:,:)*graddlocgi(k,:,:),1)*divulocgi(i,:))
             end do
             
             
             if(i<N_Layers) then
                do k = j+1, N_Layers
                   RHS(h_ele) = RHS(h_ele) - &
                        shape_rhs(mesh%nh,&
                        detwei*rho(j)*dlocgi(j,:)* &
                        divulocgi(k,:)*dlocgi(k,:)*divulocgi(i,:))
                   RHS(h_ele) = RHS(h_ele) - &
                        shape_rhs(mesh%nh,&
                        detwei*rho(j)*dlocgi(j,:)* &
                        (ulocgi(k,1,:)*graddlocgi(k,1,:)&
                        +ulocgi(k,2,:)*graddlocgi(k,2,:))*divulocgi(i,:))
                end do
             end if
             
             rhs(h_ele)=rhs(h_ele)-shape_rhs(mesh%nh,&
                  0.5*detwei*rho(j)*&
                  (dlocgi(j,:)**2.0)*divulocgi(j,:)*divulocgi(i,:))

          end do
          
          if (i == 1) then
             
             Qhh= shape_shape(mesh%nh,mesh%nh,detwei)
             call addto(mesh%mass_h,h_ele,h_ele,Qhh)

          end if

       end do
       call gallopede_solve(udotgradc((i-1)*N_dens+1:i*N_dens), &
            mesh%Mass_h,RHS,ksp_type,pc_type,1.0d-31,10000)       
  
    end do
    
    
    ewrite(1,*) 'end subroutine getudotgradc2'
    
  end subroutine getudotgradc2
  
  subroutine get_hydrostatic(rho,u,bottom,D,mesh,P)
    implicit none
    
    real, dimension(N_vels*N_Layers*2), target, intent(in) :: u
    type(dg_mesh), intent(in) :: mesh
    real, dimension(N_dens), intent(in) :: bottom
    real, dimension(N_layers) :: rho
    real, dimension(N_dens*N_layers), intent(in) :: D
    real, dimension(N_dens*N_Layers), intent(out) , target :: P
    
    !locals
    integer :: ele, i,j
    real, dimension(:), allocatable :: RHS
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(mesh%nh%ngi) :: detwei
    real, dimension(N_layers,mesh%nh%ngi) :: dlocgi
    real, dimension(mesh%nh%ngi) :: u1locgi,u2locgi,blocgi
    real, dimension(:), pointer :: PP 
    real, dimension(:), allocatable :: u_squared

    KSPType :: ksp_type
    PCType  :: pc_type

    ewrite(1,*) 'subroutine get_hydrostatic'
    
    pc_type = PCSOR  ! Preconditioner context
    ksp_type = KSPCG ! Krylov subspace context

    allocate( RHS(N_dens) )
    allocate(u_squared(N_Layers*N_dens) )
    call get_u_squared(u,mesh,u_squared)

    P = 0.

    do i = 1, N_Layers
       RHS = 0.
       do ele = 1, N_elements
          m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
          u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)

          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)

          call transform_to_physical(ele_X,mesh%nm,detwei = detwei)
    
          u1locgi = &
               matmul(transpose(mesh%nu%n),u((i-1)*2*N_vels + u_ele))
          u2locgi = &
               matmul(transpose(mesh%nu%n),u(((i-1)*2+1)*N_vels + u_ele))
          do j = 1, N_Layers
             dlocgi(j,:) = &
                  matmul(transpose(mesh%nh%n),D((j-1)*N_dens + h_ele))
          end do
          blocgi = matmul(transpose(mesh%nm%n),bottom(X_ele))

             RHS(h_ele) = RHS(h_ele) + shape_rhs(mesh%nh, &
                  detwei*rho(i)*g0*blocgi)

          end do

          PP => P((i-1)*N_dens+1:i*N_dens)
          call gallopede_solve(PP,mesh%Mass_h,&
               RHS,ksp_type,pc_type,1e-20,2000)   
          

          if(NONHYDROSTATIC_PRESSURE_FLAG) then

             P((i-1)*N_dens+1:i*N_dens)=P((i-1)*N_dens+1:i*N_dens)-&
                  rho(i)*u_squared((i-1)*N_dens+1:i*N_dens)
          

          end if

          do j = i, N_Layers
             P((i-1)*N_dens+1:i*N_dens)=P((i-1)*N_dens+1:i*N_dens)-&
                  rho(i)*g0*D((j-1)*N_dens+1:j*N_dens)
          end do
          do j = 1, i-1
             P((i-1)*N_dens+1:i*N_dens)=P((i-1)*N_dens+1:i*N_dens)-&
                  rho(j)*g0*D((j-1)*N_dens+1:j*N_dens)
          end do
       end do

    ewrite(1,*) 'End subroutine get_hydrostatic'

  end subroutine get_hydrostatic

  subroutine curl_free_m(psimat,m1,m2,prhs,D,mesh)
    implicit none
    real, intent(out), dimension(n_dens) :: prhs  
    type(csr_matrix), intent(inout) :: psimat
    real, intent(in), dimension(n_moms) :: m1, m2 
    real, intent(in), dimension(N_dens) :: D
    type(dg_mesh), intent(in) :: mesh

    integer :: ele
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele
    integer, dimension(:), pointer :: X_ele_2, h_ele_2, m_ele_2, u_ele_2

    real, dimension(mesh%nu%ngi) :: dlocgi
    real, dimension(2,mesh%nu%ngi) :: gradlocgi
    real, dimension(mesh%nm%ngi) :: m1locgi, m2locgi, curlmlocgi

    integer, dimension(:), allocatable :: neigh
    integer, dimension(:), pointer :: bdyu, bdyu_2, bdyh, bdyh_2,bdym,bdym_2
    integer :: bdy_i, ele_2, gi, nod
    real, dimension(2,mesh%nm%loc) :: ele_X, ele_X_2
    real, dimension(2,mesh%nm_f%loc) :: ele_Xf, ele_Xf_2
    real :: kmat, kmat11, kmat12, kmat21, kmat22
    integer :: ni , iloc ,jloc
    real :: detwei(mesh%nu%ngi), dlocgi_f(mesh%nh_f%ngi)
    real :: detwei_2(mesh%nu%ngi), detwei_f(mesh%nu_f%ngi)
    real :: normal(2,mesh%nu_f%ngi)

    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t,dnh_t_2
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t,dnm_t_2
    real, dimension(mesh%nm_f%loc,mesh%nm_f%ngi,2) :: dnm_ft

    integer, dimension(mesh%nh_f%loc), target :: surface_h_lno
    integer, dimension(mesh%nu_f%loc), target :: surface_u_lno
    integer, dimension(mesh%nm_f%loc), target :: surface_m_lno


    real, dimension(mesh%nh%loc,mesh%nh%loc) :: Qmm  
    real, dimension(mesh%nh_f%loc,mesh%nh_f%loc) :: Qmm_F  

    ewrite(2,*)("subroutine curl_free_m")

    ewrite(2,*)("zeroing stuff")
    call zero(psimat)

    prhs = 0.

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

       call transform_to_physical(ele_X, mesh%nm, m = mesh%nu,&
            dn_t = dnm_t, dm_t = dnu_t, detwei = detwei)
       call transform_to_physical(ele_X, mesh%nm, m = mesh%nh,&
            dm_t = dnh_t)

       !volume integrals

       Qmm=-dshape_dot_dshape(dnh_t,dnh_t,detwei)
       call addto(psimat,h_ele,h_ele,Qmm)

       dlocgi = matmul(transpose(mesh%nh%n),d(h_ele))
       gradlocgi(1,:) = matmul(transpose(dnh_t(:,:,1)),d(h_ele))
       gradlocgi(2,:) = matmul(transpose(dnh_t(:,:,2)),d(h_ele))
       m1locgi = matmul(transpose(mesh%nm%n),m1(m_ele))
       m2locgi = matmul(transpose(mesh%nm%n),m2(m_ele))
       curlmlocgi = matmul(transpose(dnm_t(:,:,1)),m2(m_ele))&
            -matmul(transpose(dnm_t(:,:,2)),m1(m_ele))

       prhs(h_ele)=shape_rhs(mesh%nh,detwei*&
            curlmlocgi/dlocgi-&
            (gradlocgi(1,:)*m2locgi-gradlocgi(2,:)*m1locgi)/dlocgi**2)

       
        !ewrite(2,*)("surface integrals")

       !===================================================================

!       if (size(neigh)/=row_length(mesh%bdy_list,ele)) then
!          ewrite(2,*)("reallocating neigh")
!          deallocate(neigh)
!          allocate(neigh(row_length(mesh%bdy_list,ele)))
!       end if
!
!       neigh=row_m(mesh%bdy_list,ele)

       !surface integrals
!       neighbourloop2: do ni=1,size(neigh)
!          ele_2=neigh(ni)

          ! check for external bdy
!         if (ele_2==0) then
!             if(.false.) then
!                cycle neighbourloop2
!             else
!                surface_h_lno = boundary_local_num(offnods(ni,:), &
!                     mesh%nh%numbering)
!                bdyh => surface_h_lno
!                surface_u_lno = boundary_local_num(offnods(ni,:), &
!                     mesh%nu%numbering)
!                bdyu => surface_u_lno
!                surface_m_lno = boundary_local_num(offnods(ni,:), &
!                     mesh%nm%numbering)
!                bdym => surface_m_lno
!             end if
!          else

!           u_ele_2=>mesh%EVList_u((ELE_2-1)*mesh%nu%LOC+1:ELE_2*mesh%nu%LOC)
!           m_ele_2=>mesh%EVList_m((ELE_2-1)*mesh%nm%LOC+1:ELE_2*mesh%nm%LOC)
!           h_ele_2=>mesh%EVList_h((ELE_2-1)*mesh%Nh%LOC+1:ELE_2*mesh%Nh%LOC)
!            X_ele_2=>mesh%EVList_X((ELE_2-1)*3+1:ELE_2*3)

!             bdy_i=ival(mesh%bdy_list, ele, ele_2)
!             bdyu=> mesh%bdy_nu_lno( &
!                  (bdy_i-1)*mesh%nu_f%loc+1:bdy_i*mesh%nu_f%loc)
!             bdyh=>mesh%bdy_nh_lno( &
!                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)
!             bdym=>mesh%bdy_nm_lno( &
!                  (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)

!             bdy_i=ival(mesh%bdy_list, ele_2, ele)
!             bdyu_2=> mesh%bdy_nu_lno( &
!                  (bdy_i-1)*mesh%nu_f%loc+1:bdy_i*mesh%nu_f%loc)
!             bdyh_2=>mesh%bdy_nh_lno( &
!                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)
!             bdym_2=>mesh%bdy_nm_lno( &
!                  (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)

             ! Locations of local vertices.
!             ele_X_2(1,:)=mesh%X(X_ele_2)
!             ele_X_2(2,:)=mesh%Y(X_ele_2)

!         end if

          ! Locations of bdy vertices.

!          ele_Xf=ele_X(:,bdym)

          ! Change of coordinates on face
!          call transform_bdy_to_physical(ele_X, ele_Xf, mesh%nm, mesh%nm_f, &
!               dn_ft=dnm_ft, detwei_f = detwei_f,normal = normal)

          !==========================================================

!           Qmm_f=dshape_dot_vector_shape(dnm_ft,normal,mesh%nm,detwei)
       
!           call addto(mommat,m_ele,m_ele,0.5*transpose(Qmm))
!           call addto(mommat,m_ele,m_ele_2(bdym_2),0.5*transpose(Qmm))

!!        end do neighbourloop2

     end do ele_loop

     deallocate( neigh )

     ewrite(2,*)("END subroutine curl_free_m")

   end subroutine curl_free_m


subroutine assemble_q_equation(psimat,prhs,m1,m2,D,Dold,u1,u2,mesh)
  implicit none
  real, intent(out), dimension(n_dens) :: prhs  
    type(csr_matrix), intent(inout) :: psimat
    real, intent(in), dimension(n_moms) :: m1, m2 
    real, intent(in), dimension(N_dens) :: D, Dold
    real, intent(in), dimension(n_vels) :: u1, u2 
    type(dg_mesh), intent(in) :: mesh

    integer :: ele
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele
    integer, dimension(:), pointer :: X_ele_2,u_ele_2, h_ele_2,m_ele_2
  
    real, dimension(2,mesh%nu%ngi) :: ulocgi
    real, dimension(2,mesh%nu_f%ngi) :: ulocgi_f
    real, dimension(mesh%nh%ngi) :: dlocgi
    real, dimension(mesh%nh%ngi) :: dlocgi_f
    real, dimension(2,mesh%nu%ngi) :: gradlocgi
    real, dimension(2,mesh%nu_f%ngi) :: gradlocgi_f
    real, dimension(mesh%nm%ngi) :: m1locgi, m2locgi, curlmlocgi
    real, dimension(mesh%nm_f%ngi) :: m1locgi_f, m2locgi_f, curlmlocgi_f
    real, dimension(mesh%nh%loc) :: curl_M
    real, dimension(mesh%nh_f%loc) :: curl_M_f


    integer, dimension(:), allocatable :: neigh
    integer, dimension(:), pointer :: bdyu, bdyu_2, bdyh, bdyh_2,bdym,bdym_2
    integer :: bdy_i, ele_2, gi, nod
    real, dimension(2,mesh%nm%loc) :: ele_X, ele_X_2
    real, dimension(2,mesh%nm_f%loc) :: ele_Xf, ele_Xf_2
    real :: kmat, kmat11, kmat12, kmat21, kmat22
    integer :: ni , iloc ,jloc
    real :: detwei(mesh%nu%ngi)
    real :: detwei_2(mesh%nu%ngi), detwei_f(mesh%nu_f%ngi)
    real :: normal(2,mesh%nu_f%ngi)

    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t,dnh_t_2
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t,dnm_t_2
    real, dimension(mesh%nm_f%loc,mesh%nm_f%ngi,2) :: dnm_ft
    real, dimension(mesh%nh_f%loc,mesh%nh_f%ngi,2) :: dnh_ft


    integer, dimension(mesh%nh_f%loc), target :: surface_h_lno
    integer, dimension(mesh%nu_f%loc), target :: surface_u_lno
    integer, dimension(mesh%nm_f%loc), target :: surface_m_lno


    real, dimension(mesh%nh%loc,mesh%nh%loc) :: Qmm  
    real, dimension(mesh%nh_f%loc,mesh%nh_f%loc) :: Qmm_f  

    ewrite(2,*)("subroutine curl_free_m")

    ewrite(2,*)("zeroing stuff")
    call zero(psimat)

    prhs = 0.

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

       call transform_to_physical(ele_X, mesh%nm, m = mesh%nu, dn_t = dnm_t, &
            dm_t = dnu_t, detwei = detwei)
       call transform_to_physical(ele_X, mesh%nm, m = mesh%nh, dm_t = dnh_t)

       !volume integrals

       dlocgi = matmul(transpose(mesh%nh%n),Dold(h_ele))
       gradlocgi(1,:) = matmul(transpose(dnh_t(:,:,1)),d(h_ele))
       gradlocgi(2,:) = matmul(transpose(dnh_t(:,:,2)),d(h_ele))
       m1locgi = matmul(transpose(mesh%nm%n),m1(m_ele))
       m2locgi = matmul(transpose(mesh%nm%n),m2(m_ele))
       curlmlocgi = matmul(transpose(dnm_t(:,:,1)),m2(m_ele))&
            -matmul(transpose(dnm_t(:,:,2)),m1(m_ele))


       Qmm=-shape_shape(mesh%nh,mesh%nh,detwei)
       call addto(psimat,h_ele,h_ele,Qmm)
       curl_M=shape_rhs(mesh%nh,detwei*&
            curlmlocgi/dlocgi-&
            (gradlocgi(1,:)*m2locgi-gradlocgi(2,:)*m1locgi)/dlocgi**2)
       prhs(h_ele)=curl_M

       ulocgi(1,:)=matmul(transpose(mesh%nu%n),u1(u_ele))
       ulocgi(2,:)=matmul(transpose(mesh%nu%n),u1(u_ele))

       Qmm=-dshape_dot_vector_shape(dnh_t,ulocgi,mesh%nh, detwei)

!       call addto(psimat,h_ele,h_ele,dt*theta_in*Qmm)
       prhs(h_ele)=prhs-dt*(1)*&
            matmul(Qmm,curl_M)
       
        !ewrite(2,*)("surface integrals")

       !===================================================================

!       if (size(neigh)/=row_length(mesh%bdy_list,ele)) then
!          ewrite(2,*)("reallocating neigh")
!          deallocate(neigh)
!          allocate(neigh(row_length(mesh%bdy_list,ele)))
!      end if

!       neigh=row_m(mesh%bdy_list,ele)

       !surface integrals
!       neighbourloop2: do ni=1,size(neigh)
!          ele_2=neigh(ni)

          ! check for external bdy
!          if (ele_2==0) then
!             if(.false.) then
!                cycle neighbourloop2
!             else
!                surface_h_lno = boundary_local_num(offnods(ni,:), &
!                     mesh%nh%numbering)
!                bdyh => surface_h_lno
!                surface_u_lno = boundary_local_num(offnods(ni,:), &
!                     mesh%nu%numbering)
!                bdyu => surface_u_lno
!                surface_m_lno = boundary_local_num(offnods(ni,:), &
!                     mesh%nm%numbering)
!                bdym => surface_m_lno
!            end if
!          else

!           u_ele_2=>mesh%EVList_u((ELE_2-1)*mesh%nu%LOC+1:ELE_2*mesh%nu%LOC)
!           m_ele_2=>mesh%EVList_m((ELE_2-1)*mesh%nm%LOC+1:ELE_2*mesh%nm%LOC)
!           h_ele_2=>mesh%EVList_h((ELE_2-1)*mesh%Nh%LOC+1:ELE_2*mesh%Nh%LOC)
!           X_ele_2=>mesh%EVList_X((ELE_2-1)*3+1:ELE_2*3)

!             bdy_i=ival(mesh%bdy_list, ele, ele_2)
!             bdyu=> mesh%bdy_nu_lno( &
!                  (bdy_i-1)*mesh%nu_f%loc+1:bdy_i*mesh%nu_f%loc)
!             bdyh=>mesh%bdy_nh_lno( &
!                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)
!             bdym=>mesh%bdy_nm_lno( &
!                 (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)

!             bdy_i=ival(mesh%bdy_list, ele_2, ele)
!             bdyu_2=> mesh%bdy_nu_lno( &
!                  (bdy_i-1)*mesh%nu_f%loc+1:bdy_i*mesh%nu_f%loc)
!             bdyh_2=>mesh%bdy_nh_lno( &
!                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)
!             bdym_2=>mesh%bdy_nm_lno( &
!                 (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)

             ! Locations of local vertices.
!             ele_X_2(1,:)=mesh%X(X_ele_2)
!             ele_X_2(2,:)=mesh%Y(X_ele_2)

!          end if

          ! Locations of bdy vertices.

!          ele_Xf=ele_X(:,bdym)

          ! Change of coordinates on face
 !         call transform_bdy_to_physical(ele_X, ele_Xf, mesh%nm, mesh%nm_f, &
 !              dn_ft=dnm_ft, detwei_f = detwei_f,normal = normal)
 !         call transform_bdy_to_physical(ele_X, ele_Xf, mesh%nu, mesh%nu_f)

          !==========================================================

 !         if (ele_2==0) then

 !         dlocgi_f = matmul(transpose(mesh%nh_f%n),&
 !              0.5*(dold(h_ele(bdyu))+d(h_ele(bdyu))))
 !         gradlocgi_f(1,:) = matmul(transpose(dnh_ft(:,:,1)),&
 !              0.5*(dold(h_ele(bdyu))+d(h_ele(bdyu))))
 !         gradlocgi_f(2,:) = matmul(transpose(dnh_ft(:,:,2)),&
 !              0.5*(dold(h_ele(bdyu))+d(h_ele(bdyu))))
 !         m1locgi_f = matmul(transpose(mesh%nm_f%n),m1(m_ele(bdyu)))
 !         m2locgi_f = matmul(transpose(mesh%nm_f%n),m2(m_ele(bdyu)))
 !         curlmlocgi_f = matmul(transpose(dnm_ft(:,:,1)),m2(m_ele(bdyu)))&
 !              -matmul(transpose(dnm_ft(:,:,2)),m1(m_ele(bdyu)))
 !         curl_M_f=shape_rhs(mesh%nm_f,detwei_f*&
 !              curlmlocgi_f/dlocgi_f-&
 !              (gradlocgi_f(1,:)*m2locgi_f-gradlocgi_f(2,:)*&
 !              m1locgi_f)/dlocgi_f**2)


 !         ulocgi(1,:)=matmul(transpose(mesh%nu_f%n),u1(u_ele(bdyu)))
 !         ulocgi(2,:)=matmul(transpose(mesh%nu_f%n),u1(u_ele(bdyu)))

 !         Qmm_f=shape_shape(mesh%nm_f,mesh%nm_f,detwei_f*&
 !              sum(normal*ulocgi,1))

 !         call addto(mommat,m_ele(bdym),m_ele(bdym), &
 !              0.5*dt*theta_in*Qmm_f)
 !         call addto(mommat,m_ele(bdym),m_ele_2(bdym_2), &
 !              0.5*dt*theta_in*Qmm_f)
 !         prhs(m_ele(bdym))=prhs+0.5*dt*(1-theta_in)*&
 !              matmul(transpose(Qmm_f),curl_M_f)
          
 !         dlocgi_f = matmul(transpose(mesh%nh_f%n),d(h_ele_2(bdyh_2)))
 !         gradlocgi_f(1,:) = matmul(transpose(dnh_ft(:,:,1)),&
 !              d(h_ele_2(bdyu_2)))
 !         gradlocgi_f(2,:) = matmul(transpose(dnh_ft(:,:,2)),&
 !              d(h_ele_2(bdyu_2)))
 !         m1locgi_f = matmul(transpose(mesh%nm_f%n),m1(m_ele_2(bdym_2)))
 !         m2locgi_f = matmul(transpose(mesh%nm_f%n),m2(m_ele_2(bdym_2)))
 !         curlmlocgi_f = matmul(transpose(dnm_ft(:,:,1)),m2(m_ele_2(bdym_2)))&
 !              -matmul(transpose(dnm_ft(:,:,2)),m1(m_ele_2(bdym_2)))
 !         curl_M_f=shape_rhs(mesh%nm_f,detwei_f*&
 !              curlmlocgi_f/dlocgi_f-&
 !              (gradlocgi_f(1,:)*m2locgi_f-gradlocgi_f(2,:)*&
 !              m1locgi_f)/dlocgi_f**2)

 !         prhs(m_ele(bdym))=prhs+0.5*dt*(1-theta_in)*&
 !              matmul(transpose(Qmm_f),curl_M_f)

 !      end if

!       end do neighbourloop2

    end do ele_loop

     deallocate( neigh )

     ewrite(2,*)("END subroutine assemble_q_equation")

   end subroutine assemble_q_equation
   
   subroutine block_assemble_MLCM_momentum_equation(mommat, &
       rhs1, rhs2, m1, m2, &
       u1, u2, D, Dold, P, bottom, rho, mesh)
    implicit none
    real, intent(out), dimension(n_moms) :: rhs1, rhs2
    real, intent(in), dimension(n_vels) :: u1, u2
    real, intent(in), dimension(n_verts) :: bottom
    real, intent(in), dimension(N_Layers) :: rho
    real, intent(in), dimension(n_moms) :: m1, m2
    real, intent(in), dimension(N_dens) :: D, Dold
    real, intent(in), dimension(N_dens) :: p
    type(block_csr_matrix), intent(inout) :: mommat
    type(dg_mesh), intent(in) :: mesh

    !locals
    integer :: ele
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele
    integer, dimension(:), pointer :: u_ele_2,h_ele_2,X_ele_2,m_ele_2
    real, dimension(2,2,mesh%nu%ngi) :: gradulocgi
    real, dimension(mesh%nu%ngi) :: plocgi,dlocgi
    real, dimension(2,mesh%nu%ngi) :: ulocgi
    real, dimension(mesh%nu%ngi,2) :: ulocgim
    real, dimension(2,mesh%nu%ngi) :: graddlocgi
    real, dimension(mesh%nh_f%ngi) :: u1locgi_f, u2locgi_f,plocgi_f
    real, dimension(2,mesh%nh_f%ngi) :: ulocgi_f

    ! List of neighbours of current element.
    integer, dimension(:), allocatable :: neigh
    integer, dimension(:), pointer :: bdyu, bdyu_2, bdyh, bdyh_2,bdym,bdym_2
    integer :: bdy_i, ele_2, gi, nod
    real, dimension(2,mesh%nm%loc) :: ele_X, ele_X_2
    real, dimension(2,mesh%nm_f%loc) :: ele_Xf, ele_Xf_2
    real :: kmat, kmat11, kmat12, kmat21, kmat22
    integer :: ni , iloc ,jloc
    real :: detwei(mesh%nu%ngi), dlocgi_f(mesh%nh_f%ngi)
    real :: detwei_2(mesh%nu%ngi), detwei_f(mesh%nu_f%ngi)
    real :: normal(2,mesh%nu_f%ngi)
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t,dnu_t_2
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t,dnh_t_2
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t,dnm_t_2

    integer, dimension(mesh%nh_f%loc), target :: surface_h_lno
    integer, dimension(mesh%nu_f%loc), target :: surface_u_lno
    integer, dimension(mesh%nm_f%loc), target :: surface_m_lno

    real :: avg_cst,theta_in
    real :: upwinder(mesh%nu_f%ngi)

    !pointer stuff
    type(layer), dimension(N_Layers) :: Dl
    real, dimension(2,mesh%nm%loc) :: dRHS
    real, dimension(2,mesh%nm_f%loc) :: dRHS_f

    real, dimension(2,2,mesh%nm%loc,mesh%nm%loc) :: TQmm
    real, dimension(mesh%nm%loc,mesh%nm%loc) :: Qmm    
    real, dimension(mesh%nm_f%loc,mesh%nm_f%loc) :: Qmm_f
    real, dimension(2,2,mesh%nm_f%loc,mesh%nm_f%loc) :: TQmm_f

    ewrite(2,*)("subroutine assemble_MLCM_momentum_equation")

    ewrite(2,*)("zeroing stuff")
    call zero(mommat)

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

       call transform_to_physical(ele_X, mesh%nm, m = mesh%nu,&
            dn_t = dnm_t, dm_t = dnu_t, detwei = detwei)
       call transform_to_physical(ele_X, mesh%nm, m = mesh%nh,&
            dm_t = dnh_t)
      

       ulocgi(1,:) = matmul(transpose(dnu_t(:,:,1)),u1(u_ele))
       ulocgi(2,:) = matmul(transpose(dnu_t(:,:,2)),u2(u_ele))

       ulocgim(:,1) = matmul(transpose(mesh%nu%n),u1(u_ele))
       ulocgim(:,2) = matmul(transpose(mesh%nu%n),u2(u_ele))

    ! Volume integrals

          dlocgi = matmul(transpose(mesh%nh%n),Dold(h_ele))
          Qmm = shape_shape(mesh%nm,mesh%nm,detwei/dlocgi)
          rhs1(m_ele) = rhs1(m_ele) + matmul(Qmm,m1(m_ele))
          rhs2(m_ele) = rhs2(m_ele) + matmul(Qmm,m2(m_ele))
          
          dlocgi = matmul(transpose(mesh%nh%n),D(h_ele))
          Qmm = shape_shape(mesh%nm,mesh%nm,detwei/dlocgi)

          call addto(mommat,1,1,m_ele,m_ele,Qmm)
          call addto(mommat,2,2,m_ele,m_ele,Qmm)

          ulocgim(:,1) = matmul(transpose(mesh%nu%n),u1(u_ele))
          ulocgim(:,2) = matmul(transpose(mesh%nu%n),u2(u_ele))

          if(NONLINEAR_FLAG) then

             ! Nonlinear term
             ! -dw u.m

             dlocgi = matmul(transpose(mesh%nh%n),0.5*(D(h_ele)+Dold(h_ele)))
             
             TQmm=-dshape_outer_vector_shape(dnm_t,transpose(ulocgim),mesh%nm,&
                  detwei/dlocgi)
                  
          
             rhs1(m_ele)=rhs1(m_ele) - (1-theta_in)*dt*&
                  (matmul(TQmm(1,1,:,:),m1(m_ele))+&
                  matmul(TQmm(1,2,:,:),m2(m_ele)))
             rhs2(m_ele)=rhs2(m_ele) - (1-theta_in)*dt*&
                  (matmul(TQmm(2,1,:,:),m1(m_ele))+&
                  matmul(TQmm(2,2,:,:),m2(m_ele)))
             
             dlocgi = matmul(transpose(mesh%nh%n),D(h_ele))

             TQmm=-dshape_outer_vector_shape(dnm_t,transpose(ulocgim),mesh%nm,&
                 detwei/dlocgi)

             call addto(mommat,1,1,m_ele,m_ele,dt*theta_in*TQmm(1,1,:,:))
             call addto(mommat,1,2,m_ele,m_ele,dt*theta_in*TQmm(1,2,:,:))
             call addto(mommat,2,1,m_ele,m_ele,dt*theta_in*TQmm(2,1,:,:))
             call addto(mommat,2,2,m_ele,m_ele,dt*theta_in*TQmm(2,2,:,:))

          end if

          !Pressure term 

   
          plocgi = matmul(transpose(mesh%nh%n),P(h_ele))

          drhs = dshape_rhs(dnm_t,detwei*plocgi)
       
          rhs1(m_ele) = rhs1(m_ele) - dt*drhs(1,:)
          rhs2(m_ele) = rhs2(m_ele) - dt*drhs(2,:)

          !ewrite(2,*)("surface integrals")

       !===================================================================

       if (size(neigh)/=row_length(mesh%bdy_list,ele)) then
          ewrite(2,*)("reallocating neigh")
          deallocate(neigh)
          allocate(neigh(row_length(mesh%bdy_list,ele)))
       end if

       neigh=row_m(mesh%bdy_list,ele)

       !surface integrals
       neighbourloop3: do ni=1,size(neigh)
          ele_2=neigh(ni)

          ! check for external bdy
          if (ele_2==0) then
             if(.false.) then
                cycle neighbourloop3
             else
                avg_cst = 1.0
                surface_h_lno = boundary_local_num(offnods(ni,:), &
                     mesh%nh%numbering)
                bdyh => surface_h_lno
                surface_u_lno = boundary_local_num(offnods(ni,:), &
                     mesh%nu%numbering)
                bdyu => surface_u_lno
                surface_m_lno = boundary_local_num(offnods(ni,:), &
                     mesh%nm%numbering)
                bdym => surface_m_lno
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
             bdyh=>mesh%bdy_nh_lno( &
                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)
             bdym=>mesh%bdy_nm_lno( &
                  (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)

             bdy_i=ival(mesh%bdy_list, ele_2, ele)
             bdyu_2=> mesh%bdy_nu_lno( &
                  (bdy_i-1)*mesh%nu_f%loc+1:bdy_i*mesh%nu_f%loc)
             bdyh_2=>mesh%bdy_nh_lno( &
                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)
             bdym_2=>mesh%bdy_nm_lno( &
                  (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)

             ! Locations of local vertices.
             ele_X_2(1,:)=mesh%X(X_ele_2)
             ele_X_2(2,:)=mesh%Y(X_ele_2)

          end if

          ! Locations of bdy vertices.

          ele_Xf=ele_X(:,bdym)

          ! Change of coordinates on face
          call transform_bdy_to_physical(ele_X, ele_Xf, mesh%nm, mesh%nm_f, &
               detwei_f = detwei_f,normal = normal)

          !non-linear surface integral
          
          ! m.u  w.n
          if(NONLINEAR_FLAG) then


             if(ele_2.ne.0) then
                
                ulocgi_f(1,:) = matmul(transpose(mesh%nu_f%n), &
                     u1(u_ele(bdyu)))
                ulocgi_f(2,:)= matmul(transpose(mesh%nu_f%n), &
                     u2(u_ele(bdyu)))
                dlocgi_f = matmul(transpose(mesh%nh_f%n),&
                     0.5*(d(h_ele(bdyh))+Dold(h_ele(bdyh))))

                TQmm_f = shape_shape_vector_outer_vector(mesh%nm_f,&
                     mesh%nm_f ,detwei_f/dlocgi_f,&
                     normal,ulocgi_f)

                if(NOUPWIND_FLAG) then
                   upwinder=0.5
                else
                   upwinder=0;
                   where(sum(normal*ulocgi_f,1) > 0) upwinder=1
                end if

!                call addto(mommat,1,1,m_ele(bdym),m_ele(bdym), &
!                     dt*theta_in*0.5*TQmm_f(1,1,:,:))
!                call addto(mommat,1,2,m_ele(bdym),m_ele(bdym), &
!                     dt*theta_in*0.5*TQmm_f(1,2,:,:))
!                call addto(mommat,2,1,m_ele(bdym),m_ele(bdym), &
!                     dt*theta_in*0.5*TQmm_f(2,1,:,:))
!                call addto(mommat,2,2,m_ele(bdym),m_ele(bdym), &
!                     dt*theta_in*0.5*TQmm_f(2,2,:,:))

                rhs1(m_ele(bdym)) = rhs1(m_ele(bdym)) - &
                     dt*(1)*upwinder*(matmul(TQmm_f(1,1,:,:),&
                     m1(m_ele(bdym)))+&
                     matmul(TQmm_f(1,2,:,:),m2(m_ele(bdym))))
                rhs2(m_ele(bdym)) = rhs2(m_ele(bdym)) - &
                     dt*(1)*upwinder*(matmul(TQmm_f(2,1,:,:),&
                     m1(m_ele(bdym)))+&
                     matmul(TQmm_f(2,2,:,:),m2(m_ele(bdym))))
               

                !outside face

!                ulocgi_f(1,:) = matmul(transpose(mesh%nu_f%n), &
!                     u1(u_ele_2(bdyu_2)))
!                ulocgi_f(2,:)= matmul(transpose(mesh%nu_f%n), &
!                     u2(u_ele_2(bdyu_2)))
!                 dlocgi_f = matmul(transpose(mesh%nh_f%n),&
!                      0.5*(d(h_ele(bdyh))+Dold(h_ele(bdyh))))
!                TQmm_f = shape_shape_vector_outer_vector(mesh%nm_f,mesh%nm_f,&
!                     detwei_f/dlocgi_f,normal,ulocgi_f)
!                call addto(mommat,1,1,m_ele(bdym),m_ele(bdym), &
!                     dt*theta_in*0.5*TQmm_f(1,1,:,:))
!                call addto(mommat,1,2,m_ele(bdym),m_ele(bdym), &
!                     dt*theta_in*0.5*TQmm_f(1,2,:,:))
!                call addto(mommat,2,1,m_ele(bdym),m_ele(bdym), &
!                     dt*theta_in*0.5*TQmm_f(2,1,:,:))
!                call addto(mommat,2,2,m_ele(bdym),m_ele(bdym), &
!                     dt*theta_in*0.5*TQmm_f(2,2,:,:))

                rhs1(m_ele(bdym)) = rhs1(m_ele(bdym)) - &
                     dt*(1)*(1-upwinder)*(matmul(TQmm_f(1,1,:,:),&
                     m1(m_ele_2(bdym_2)))+&
                     matmul(TQmm_f(1,2,:,:),m2(m_ele_2(bdym_2))))
                rhs2(m_ele(bdym)) = rhs2(m_ele(bdym)) - &
                     dt*(1)*(1-upwinder)*(matmul(TQmm_f(2,1,:,:),&
                     m1(m_ele_2(bdym_2)))+&
                     matmul(TQmm_f(2,2,:,:),m2(m_ele_2(bdym_2))))
                
             end if

          end if

          !pressure surface integral
          if(.true.) then

             !Pn
             !Pressure is continuous so only do inside face
             plocgi_f = matmul(transpose(mesh%nh_f%n),P(h_ele(bdyh)))

            dRHS_f = -shape_vector_rhs(mesh%nm_f,normal, &
                 detwei_f*plocgi_f)

             RHS1(m_ele(bdym)) = RHS1(m_ele(bdym)) - dt*dRHS_f(1,:)
             RHS2(m_ele(bdym)) = RHS2(m_ele(bdym)) - dt*dRHS_f(2,:)

          end if
       end do neighbourloop3

    end do ele_loop
    deallocate( neigh )

    ewrite(2,*)("END subroutine assemble_MLCM_momentum_equation")

  end subroutine block_assemble_MLCM_momentum_equation

subroutine block_assemble_MLCM_momentum_equation_2(mommat, &
       rhs,m1,m2, phi, &
       u1, u2, P,D, bottom, rho, mesh)
    implicit none
    real, intent(out), dimension(n_dens) :: rhs
    real, intent(in), dimension(n_vels) :: u1, u2
    real, intent(in), dimension(n_moms) :: m1, m2
    real, intent(in), dimension(n_verts) :: bottom
    real, intent(in), dimension(N_Layers) :: rho
    real, intent(in), dimension(n_dens) :: phi
    real, intent(in), dimension(N_dens) :: p, D
    type(csr_matrix), intent(inout) :: mommat
    type(dg_mesh), intent(in) :: mesh

    !locals
    integer :: ele
    integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele
    integer, dimension(:), pointer :: u_ele_2,h_ele_2,X_ele_2,m_ele_2
    real, dimension(2,2,mesh%nu%ngi) :: gradulocgi
    real, dimension(mesh%nu%ngi) :: plocgi,dlocgi
    real, dimension(2,mesh%nu%ngi) :: gradplocgi
    real, dimension(2,mesh%nu%ngi) :: ulocgi
    real, dimension(2,mesh%nu%ngi) :: mlocgi
    real, dimension(2,mesh%nu%ngi) :: graddlocgi
    real, dimension(2,2,mesh%nu%ngi) :: gradmlocgi
    real, dimension(2,mesh%nu%ngi) :: rhslocgi

    ! List of neighbours of current element.
    
    real, dimension(2,mesh%nm%loc) :: ele_X, ele_X_2
    real, dimension(2,mesh%nm_f%loc) :: ele_Xf, ele_Xf_2
    real :: kmat, kmat11, kmat12, kmat21, kmat22
    integer :: ni , iloc ,jloc
    real :: detwei(mesh%nu%ngi)
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t

    !pointer stuff
    type(layer), dimension(N_Layers) :: Dl
    real, dimension(mesh%nh%loc) :: dRHS

    real, dimension(2,2,mesh%nh%loc,mesh%nh%loc) :: Qdhdh 
    real, dimension(mesh%nh%loc,mesh%nh%loc) :: Qhh  

    ewrite(1,*)("subroutine assemble_MLCM_momentum_equation_2")

    ewrite(1,*)("zeroing stuff")
    call zero(mommat)

    rhs = 0.

    ele_loop: do ele = 1, n_elements

       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)

       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)

       call transform_to_physical(ele_X, mesh%nm, m = mesh%nu,&
            dn_t = dnm_t, dm_t = dnu_t, detwei = detwei)
       call transform_to_physical(ele_X, mesh%nm, m = mesh%nh,&
            dm_t = dnh_t)
      
    ! Volume integrals

          Qhh = -dshape_dot_dshape(dnh_t,dnh_t,detwei)
!          rhs(h_ele) = rhs(h_ele) + matmul(Qhh,phi(h_ele))

          call addto(mommat,h_ele,h_ele,Qhh)

          if(NONLINEAR_FLAG) then

             ulocgi(1,:) = matmul(transpose(mesh%nu%n),u1(u_ele))
             ulocgi(2,:) = matmul(transpose(mesh%nu%n),u2(u_ele)) 
             mlocgi(1,:) = matmul(transpose(mesh%nm%n),m1(m_ele))
             mlocgi(2,:) = matmul(transpose(mesh%nm%n),m2(m_ele)) 
             dlocgi= matmul(transpose(mesh%nh%n),D(h_ele))

             gradulocgi(1,1,:) = matmul(transpose(dnu_t(:,:,1)),u1(u_ele))
             gradulocgi(1,2,:) = matmul(transpose(dnu_t(:,:,2)),u1(u_ele))
             gradulocgi(2,1,:) = matmul(transpose(dnu_t(:,:,1)),u2(u_ele))
             gradulocgi(2,2,:) = matmul(transpose(dnu_t(:,:,2)),u2(u_ele))
             gradmlocgi(1,1,:) = matmul(transpose(dnm_t(:,:,1)),m1(m_ele))
             gradmlocgi(1,2,:) = matmul(transpose(dnm_t(:,:,2)),m1(m_ele))
             gradmlocgi(2,1,:) = matmul(transpose(dnm_t(:,:,1)),m2(m_ele))
             gradmlocgi(2,2,:) = matmul(transpose(dnm_t(:,:,2)),m2(m_ele))
             graddlocgi(1,:) = matmul(transpose(dnh_t(:,:,1)),D(h_ele))

             rhslocgi(1,:)=&
                  (ulocgi(2,:)*(gradmlocgi(2,1,:)-gradmlocgi(1,2,:))/dlocgi&
                  +mlocgi(2,:)*(gradulocgi(2,1,:)-gradulocgi(1,2,:))/dlocgi&
                  +sum(ulocgi*gradmlocgi(1,:,:),1)/dlocgi&
                  +sum(mlocgi*gradulocgi(1,:,:),1)/dlocgi&
                  -graddlocgi(1,:)*sum(ulocgi*mlocgi,1)/(dlocgi**2))

             rhslocgi(2,:)=&
                  (-ulocgi(1,:)*(gradmlocgi(2,1,:)-gradmlocgi(1,2,:))/dlocgi&
                  -mlocgi(1,:)*(gradulocgi(2,1,:)-gradulocgi(1,2,:))/dlocgi&
                  +sum(ulocgi*gradmlocgi(2,:,:),1)/dlocgi&
                  +sum(mlocgi*gradulocgi(2,:,:),1)/dlocgi&
                  -graddlocgi(2,:)*sum(ulocgi*mlocgi,1)/(dlocgi**2))

                
             forall (iloc=1:mesh%nh%loc)
                drhs(iloc)=-dot_product(sum(dnh_t(iloc,:,:)*&
                     transpose(rhslocgi),2),detwei)
             end forall
            
             rhs(h_ele)=rhs(h_ele)-dt*drhs


!             call addto(mommat,h_ele,h_ele,dt*theta_in*Qhh(:,:))
!              rhs(h_ele)=rhs(h_ele) - (1-theta_in)*dt*&
!                  (matmul(Qhh(:,:),phi(h_ele))+&

          end if

          !Pressure term 

          gradplocgi(1,:) = matmul(transpose(dnh_t(:,:,1)),P(h_ele))
          gradplocgi(2,:) = matmul(transpose(dnh_t(:,:,2)),P(h_ele))


          forall (iloc=1:mesh%nh%loc)
             drhs(iloc)=-dot_product(sum(dnh_t(iloc,:,:)*&
                  transpose(gradplocgi),2),detwei)
          end forall

          rhs(h_ele)=rhs(h_ele)+dt*drhs     

       end do ele_loop

    ewrite(1,*)("END subroutine assemble_MLCM_momentum_equation_2")

  end subroutine block_assemble_MLCM_momentum_equation_2

subroutine div_m_layer(mommat,rhs1,rhs2,m1,m2,psi,D,mesh)
  implicit none
type(csr_matrix), intent(inout) ::mommat
real, intent(out), dimension(n_moms):: rhs1, rhs2
real, intent(inout), dimension(n_moms) :: m1, m2
real, intent(in), dimension(n_dens):: psi
real, intent(in), dimension(n_dens):: D
type(dg_mesh), intent(in) :: mesh

integer :: ele
integer, dimension(:), pointer :: X_ele,u_ele, h_ele,m_ele
real, dimension(mesh%nh%ngi) :: dlocgi
real, dimension(2, mesh%nm%ngi) ::gradpsilocgi, nonmlocgi
real, dimension(2,mesh%nm%loc) :: ele_X

real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
real :: detwei(mesh%nu%ngi)

 real, dimension(mesh%nm%loc,mesh%nm%loc) :: Qmm  


 ewrite(2,*)("subroutine curl_free_m")

    ewrite(2,*)("zeroing stuff")
    call zero(mommat)

    rhs1 = 0.
    rhs2 = 0.

ele_loop: do ele = 1, n_elements

       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)

       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)

       call transform_to_physical(ele_X, mesh%nm, m = mesh%nu, dn_t = dnm_t, &
            dm_t = dnu_t, detwei = detwei)
       call transform_to_physical(ele_X, mesh%nm, m = mesh%nh, dm_t = dnh_t)

       !volume integrals


       dlocgi = matmul(transpose(mesh%nh%n),d(h_ele))
       gradpsilocgi(1,:)=matmul(transpose(dnh_t(:,:,1)),psi(h_ele))
       gradpsilocgi(2,:)=matmul(transpose(dnh_t(:,:,2)),psi(h_ele))
       
       Qmm=shape_shape(mesh%nm,mesh%nm,detwei/dlocgi)

       call addto(mommat,m_ele,m_ele,Qmm)
       rhs1(m_ele)=matmul(Qmm,m1(m_ele))
       rhs2(m_ele)=matmul(Qmm,m2(m_ele))
          
 !      drhs(m_ele)=dshape_rhs(dnm_t,detwei*psilocgi)
       

       rhs1(m_ele)=rhs1(m_ele)+shape_rhs(mesh%nm,detwei*gradpsilocgi(1,:))
       rhs2(m_ele)=rhs2(m_ele)+shape_rhs(mesh%nm,detwei*gradpsilocgi(2,:))

    end do ele_loop

  end subroutine div_m_layer

end module MLCM_momentum_equation
