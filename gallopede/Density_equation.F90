
#include "fdebug.h"
module density_equation

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
!  use gallopede_solvers
  use data_structures
  use fldebug
  use FETools

  implicit none

  !=============================================================
  !TODO
  !weak flux conditions for inlets
  !=============================================================

  public solve_density_equation, density_smoother, load_quad_d, energy_find, &
       solve_interface_equation, get_density_relation

  private
  
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
  
contains
  
  subroutine solve_density_equation(D,u1,u2,Mesh)
    
    implicit none
    !code to solve the advection diffusion equation
    
    real, dimension(n_dens), intent(inout) :: D
    real, dimension(n_vels), intent(in) :: u1, u2
    type(dg_mesh) :: Mesh

    !local variables
    type(csr_matrix) :: den_eqn_mat
    integer :: ele, iloc, jloc, gi, globi, globj, ele_2, iloc2,jloc2
    integer, dimension(:), pointer :: u_ele, h_ele, X_ele
    real, dimension(2,mesh%nx%loc) :: ele_X
    real :: kmat, dkmat,h_supg,ab_u
    real, dimension(n_dens) :: rhs, newD
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nh%ngi):: detwei
    real, dimension(mesh%nu%ngi):: divulocgi
  

    real, dimension(mesh%nh%loc,mesh%nh%loc) :: Qhh,Qhh2
    real, dimension(2,mesh%nu%ngi) :: ulocgi
    real, dimension(2,mesh%nh%loc) :: ele_Xh


    integer :: nits 
    real :: totmass,vol, tau
    logical, parameter :: mass_sqr = .false., SUPG_FLAG= .false.
    
    KSPType :: ksp_type
    PCType :: pc_type

    assert( size(D)  == n_dens )
    assert( size(u1) == n_vels )
    assert( size(u2) == n_vels )

!    where (abs(u1)>1e3*abs(u2)) u2=0
!    where (abs(u2)>1e3*abs(u1)) u1=0
    
    ewrite(1,*)( "subroutine solve_density_equation" )

    ewrite(2,*) shape(mesh%nu%dn)

    !assemble matrix and RHS

    call allocate(den_eqn_mat,mesh%Mass_h%sparsity)

    nits = 0

    newD = D

    totmass = get_totmass(mesh,D)
    ewrite(2,*)(totmass)

 !   nonlinear_loop: do
 !      if(nits==den_maxnits) exit

       rhs = 0.0

       

       call zero(den_eqn_mat)

       nits = nits + 1

       element_loop: do ele = 1, N_elements
          u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)

          !ele_Xh(:,(/1,3,6/))=ele_X
          !ele_Xh(:,(/2,4,5/))=0.5*(ele_X(:,(/1,3,2/))+ele_X(:,(/2,1,3/)))

          call transform_to_physical(ele_X, n=mesh%nx, m = mesh%nh, &
               dm_t = dnh_t, detwei = detwei)
          call transform_to_physical(ele_X, n=mesh%nx, m = mesh%nu, &
               dm_t = dnu_t)

          !construct u at Gauss points

          ulocgi(1,:) = matmul(u1(u_ele),mesh%nu%n)
          ulocgi(2,:) = matmul(u2(u_ele),mesh%nu%n)
          divulocgi =   matmul(u1(u_ele),dnu_t(:,:,1))&
               +matmul(u2(u_ele),dnu_t(:,:,2))
          ab_u=sqrt(sum(u1(u_ele)**2+u2(u_ele)**2)/mesh%nu%loc)
!          where (abs(ab_u)<1e-10) ab_u=1e-10
          h_supg=0.0
!          do iloc=1,mesh%nu%loc
!          h_supg=h_supg+abs(sum(transpose(ulocgi)*dnu_t(iloc,:,:)))/mesh%nh%loc
!          end do
!          h_supg=2.0*ab_u/h_supg
!          if (sqrt(sum(detwei)/8.0)*(3.0/mesh%nh%loc)*&
!               sum(ab_u/mesh%nu%ngi)> 6*kappa) then
!             tau=0.5*sqrt(sum(detwei)/8.0)*(3.0/mesh%nh%loc)
!          if (ab_u*h_supg/2&
!              .ge. 6.0*kappa) then
!             tau=1.0/(2.0*sqrt(2.0*sum(ulocgi**2/(2.0*mesh%nu%ngi))))
!             tau=h_supg/(2.0*ab_u)
!             tau=dt/2.0
!          else
!             tau=h_supg**2/(12*kappa)
!          end if
          !matrix contributions


          forall (iloc=1:mesh%nh%loc,jloc=1:mesh%nh%loc)
             Qhh(iloc,jloc)=sum(sum(dnh_t(iloc,:,:)*transpose(ulocgi),2)&
                  *sum(dnh_t(jloc,:,:)*transpose(ulocgi),2)*detwei)
          end forall

          h_supg=sqrt(sum(dshape_dot_vector_shape(dnh_t,ulocgi,&
               mesh%nh,detwei)**2)&
               /sum(Qhh**2))

          tau=1.0/sqrt(1.0/h_supg**2 + 4.0/dt**2)
          if (ab_u>0.0)&
               tau=1.0/sqrt(1.0/h_supg**2 + 4.0/dt**2&
               +kappa**2/(h_supg*ab_u**2)**2)
          

          !bulk integrals
          Qhh = shape_shape(mesh%nh,mesh%nh,detwei)
            call addto(den_eqn_mat,h_ele,h_ele,Qhh)
          rhs(h_ele) = rhs(h_ele) + matmul(Qhh,D(h_ele))

          !Advection 

          Qhh=-dshape_dot_vector_shape(dnh_t,ulocgi,mesh%nh,detwei)
!          Qhh=shape_shape(mesh%nh,mesh%nh,detwei*divulocgi)

          rhs(h_ele) = rhs(h_ele) &
               - dt*(1.0-theta)*matmul(Qhh,D(h_ele))
            call addto(den_eqn_mat,h_ele,h_ele,dt*theta*Qhh)




            forall(iloc=1:mesh%nh%loc,jloc=1:mesh%nh%loc)
               Qhh(iloc,jloc)=sum(mesh%nh%n(iloc,:)&
                    *(dnh_t(jloc,:,1)*ulocgi(1,:)&
                    +dnh_t(jloc,:,2)*ulocgi(2,:))&
                    *detwei)
            end forall

!            rhs(h_ele) = rhs(h_ele) &
!                 - dt*(1.0-theta)*matmul(Qhh,D(h_ele))
!            call addto(den_eqn_mat,h_ele,h_ele,dt*theta*Qhh)

!            newD(h_ele)= 0.001*cos(2*3.1415927*(ele_Xh(2,:)-0.1*t))

!            Diffusion

             Qhh=dshape_dot_dshape(dnh_t,dnh_t,detwei*kappa)

             rhs(h_ele) = rhs(h_ele) &
               - dt*(1.0-theta)*matmul(Qhh,D(h_ele))
             call addto(den_eqn_mat,h_ele,h_ele,dt*theta*Qhh)

             if (SUPG_flag) then
                 Qhh = dshape_dot_vector_shape(dnh_t,ulocgi,mesh%nh,detwei)&
                      +shape_shape(mesh%nh,mesh%nh,detwei*divulocgi)
                 Call addto(den_eqn_mat,h_ele,h_ele,tau*Qhh)
                 rhs(h_ele) = rhs(h_ele) + tau*matmul(Qhh,D(h_ele))

                 forall (iloc=1:mesh%nh%loc,jloc=1:mesh%nh%loc)
                    Qhh(iloc,jloc)=sum(sum(dnh_t(iloc,:,:)*transpose(ulocgi),2)&
                      *sum(transpose(dnh_t(jloc,:,:))*ulocgi,1)*detwei)&
                      +sum(mesh%nh%n(iloc,:)*sum(transpose(dnh_t(jloc,:,:))&
                      *ulocgi,1)*detwei*divulocgi)
                 end forall

                 rhs(h_ele) = rhs(h_ele) &
                      - tau*dt*(1.0-theta)*matmul(Qhh,D(h_ele))
                 call addto(den_eqn_mat,h_ele,h_ele,tau*dt*theta*Qhh)
                 
                 Qhh = dshape_dot_vector_shape(dnh_t,ulocgi,mesh%nh,&
                      detwei*divulocgi)&
                      +shape_shape(mesh%nh,mesh%nh,detwei*divulocgi**2)        
      
                 rhs(h_ele) = rhs(h_ele) &
                      - tau*dt*(1.0-theta)*matmul(Qhh,D(h_ele))
                 call addto(den_eqn_mat,h_ele,h_ele,tau*dt*theta*Qhh)

              end if

       end do element_loop

       !solve density equation

       ewrite(1,*)("Solving density equation")

       ksp_type = KSPGMRES
       pc_type = PCSOR

       call gallopede_solve(newD, den_eqn_mat, rhs,&
            ksp_type,pc_type, 1.0d-20, 30000)
   
!    end do nonlinear_loop

    D = newD

    call deallocate( den_eqn_mat )

    ewrite(1,*) ("end subroutine solve_density_equation")

  end subroutine solve_density_equation

 subroutine solve_interface_equation(h,u1,u2,Mesh,theta_in,bottom)
    
    implicit none
    !code to solve the interface transport equation
    
    real, dimension(n_dens*n_layers), intent(inout) :: h
    real, dimension(n_vels*n_layers), intent(in) :: u1, u2
    real, dimension(n_dens), intent(in) :: bottom 
    type(dg_mesh) :: Mesh

    !local variables
    type(block_csr_matrix) :: den_eqn_mat
    integer :: ele, iloc, jloc, gi, globi, globj, ele_2, iloc2,jloc2
    integer, dimension(:), pointer :: u_ele, h_ele, X_ele
    real, dimension(2,mesh%nx%loc) :: ele_X
    real :: kmat, dkmat,theta_in
    real, dimension(n_dens*n_layers) :: rhs, newh
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nh%ngi):: detwei
    real, dimension(mesh%nu%ngi):: divulocgi
  

    real, dimension(mesh%nh%loc,mesh%nh%loc) :: Qhh,Qhh2
    real, dimension(2,mesh%nh%loc,mesh%nh%loc) :: Quh
    real, dimension(2,mesh%nu%ngi) :: ulocgi
    real, dimension(2,mesh%nh%loc) :: ele_Xh


    integer :: nits,j,k,top_layers
    real :: totmass,vol, tau
    logical, parameter :: mass_sqr = .false., SUPG_FLAG= .false.
    
    KSPType :: ksp_type
    PCType :: pc_type

    assert( size(h)  == n_dens )
    assert( size(u1) == n_vels )
    assert( size(u2) == n_vels )

!    where (abs(u1)>1e3*abs(u2)) u2=0
!    where (abs(u2)>1e3*abs(u1)) u1=0
    
    ewrite(1,*)( "subroutine solve_interface_equation" )

    ewrite(2,*) shape(mesh%nu%dn)

    !assemble matrix and RHS

   call allocate(den_eqn_mat,mesh%Mass_h%sparsity, (/n_layers,n_layers/))

    nits = 0

    newh = h

    totmass = get_totmass(mesh,h)
    ewrite(2,*)(totmass)

 !   nonlinear_loop: do
 !      if(nits==den_maxnits) exit

       rhs = 0.0

       

       call zero(den_eqn_mat)

       nits = nits + 1

       element_loop: do ele = 1, N_elements
          u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)

          !ele_Xh(:,(/1,3,6/))=ele_X
          !ele_Xh(:,(/2,4,5/))=0.5*(ele_X(:,(/1,3,2/))+ele_X(:,(/2,1,3/)))

          call transform_to_physical(ele_X, n=mesh%nx, m = mesh%nh, &
               dm_t = dnh_t, detwei = detwei)
          call transform_to_physical(ele_X, n=mesh%nx, m = mesh%nu, &
               dm_t = dnu_t)


          tau=1.0/2.0*sqrt(sum(detwei)/8.0)*(3.0/mesh%nh%loc)



          divulocgi =   matmul(u1(u_ele),dnu_t(:,:,1))&
               +matmul(u2(u_ele),dnu_t(:,:,2))

          !matrix contributions

          !bulk integrals

          
          Qhh = shape_shape(mesh%nh,mesh%nh,detwei)
          DO J=1,N_LAYERS
            call addto(den_eqn_mat,j,j,h_ele,h_ele,Qhh)
            rhs(h_ele+(j-1)*n_dens) = rhs(h_ele+(j-1)*n_dens)&
                 + matmul(Qhh,h(h_ele+(j-1)*n_dens))
         end do

          !Advection 

         

         
         top_layers=1


         do j=top_layers,n_layers
           if (j<n_layers) then
            do k=j,N_layers-1
          !construct u at Gauss points

          ulocgi(1,:) = matmul(u1(u_ele+(k-1)*N_vels),mesh%nu%n)
          ulocgi(2,:) = matmul(u2(u_ele+(k-1)*N_vels),mesh%nu%n)
          Qhh=-dshape_dot_vector_shape(dnh_t,ulocgi,mesh%nh,detwei)

              rhs(h_ele+(j-1)*n_dens) = rhs(h_ele+(j-1)*n_dens) &
               - (1.0-theta_in)*dt*&
               matmul(Qhh,h(h_ele+(k-1)*N_dens)-h(h_ele+k*n_dens))
              call addto(den_eqn_mat,j,k,h_ele,h_ele,dt*theta_in*Qhh)
              call addto(den_eqn_mat,j,k+1,h_ele,h_ele,-dt*theta_in*Qhh)
           end do

          ulocgi(1,:) = matmul(u1(u_ele+(n_layers-1)*N_vels),mesh%nu%n)
          ulocgi(2,:) = matmul(u2(u_ele+(n_layers-1)*N_vels),mesh%nu%n)

          Qhh=-dshape_dot_vector_shape(dnh_t,ulocgi,mesh%nh,detwei)
               rhs(h_ele+(j-1)*n_dens) = rhs(h_ele+(j-1)*n_dens) &
               - (1.0-theta)*dt*matmul(Qhh,h(h_ele+(n_layers-1)*N_dens))
               rhs(h_ele+(j-1)*n_dens) = rhs(h_ele+(j-1)*n_dens) &
                    -dt*matmul(Qhh,bottom(h_ele))
              call addto(den_eqn_mat,j,n_layers,h_ele,h_ele,dt*theta_in*Qhh)

            else
               ulocgi(1,:) = matmul(u1(u_ele+(j-1)*N_vels),mesh%nu%n)
               ulocgi(2,:) = matmul(u2(u_ele+(j-1)*N_vels),mesh%nu%n)

               Qhh=-dshape_dot_vector_shape(dnh_t,ulocgi,mesh%nh,detwei)
               rhs(h_ele+(j-1)*n_dens) = rhs(h_ele+(j-1)*n_dens) &
                    - (1.0-theta_in)*dt*matmul(Qhh,h(h_ele+(j-1)*N_dens))
               rhs(h_ele+(j-1)*n_dens) = rhs(h_ele+(j-1)*n_dens) &
                    -dt*matmul(Qhh,bottom(h_ele))
              call addto(den_eqn_mat,n_layers,n_layers,h_ele,h_ele,dt*theta_in*Qhh)
           end if
         end do

             Qhh=dshape_dot_dshape(dnh_t,dnh_t,detwei*kappa)

             rhs(h_ele) = rhs(h_ele) &
               - dt*(1.0-theta_in)*matmul(Qhh,h(h_ele))
             call addto(den_eqn_mat,1,1,h_ele,h_ele,dt*theta_in*Qhh)

             do k=2,n_layers
                rhs(h_ele+(k-1)*n_dens) = rhs(h_ele+(k-1)*n_dens) &
                     - 0.001*dt*(1.0-theta_in)*matmul(Qhh,h(h_ele+(k-1)*n_dens))
                call addto(den_eqn_mat,k,k,h_ele,h_ele,0.001*dt*theta_in*Qhh) 
             end do


             if (SUPG_flag) then
                 Qhh = dshape_dot_vector_shape(dnh_t,ulocgi,mesh%nh,&
                      tau*detwei)
                 Call addto(den_eqn_mat,1,1,h_ele,h_ele,Qhh)
                 rhs(h_ele) = rhs(h_ele) + matmul(Qhh,h(h_ele))

                 forall (iloc=1:mesh%nh%loc,jloc=1:mesh%nh%loc)
                    Qhh(iloc,jloc)=sum(sum(dnh_t(iloc,:,:)*transpose(ulocgi),2)&
                         *sum(dnh_t(jloc,:,:)*transpose(ulocgi),2)*tau*detwei)
                 end forall
                 
                 rhs(h_ele) = rhs(h_ele) &
                      - dt*(1.0-theta_in)*matmul(Qhh,h(h_ele))
                 call addto(den_eqn_mat,1,1,h_ele,h_ele,dt*theta_in*Qhh)
                 
                 Qhh = dshape_dot_vector_shape(dnh_t,ulocgi,mesh%nh,&
                      tau*detwei*divulocgi)

                 rhs(h_ele) = rhs(h_ele) &
                      - dt*(1.0-theta_in)*matmul(Qhh,h(h_ele))
                 call addto(den_eqn_mat,1,1,h_ele,h_ele,dt*theta_in*Qhh)

              end if

       end do element_loop

       !solve density equation

       ewrite(1,*)("Solving density equation")

       ksp_type = KSPGMRES
       pc_type = PCSOR

       call gallopede_block_solve(newh, den_eqn_mat, rhs,&
            ksp_type,pc_type, 1.0d-31, 30000)
   
!    end do nonlinear_loop

    h = newh

    call deallocate( den_eqn_mat )

    ewrite(1,*) ("end subroutine solve_external_density")

  end subroutine solve_interface_equation




  function get_totmass(mesh,D) result(totmass)
    type(dg_mesh), intent(in) :: mesh
    real, dimension(n_dens), intent(in) :: D
    real :: totmass
    
    !locals
    integer :: ele, iloc, globi
    real, dimension(mesh%nh%ngi) :: detwei
    integer, dimension(:), pointer :: u_ele,X_ele,h_ele
    real, dimension(2,mesh%nx%loc) :: ele_X
    
    totmass = 0.
    do ele = 1, N_elements
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)
       call transform_to_physical(ele_X, n=mesh%nx, detwei = detwei)
       do iloc = 1, mesh%nh%loc
          globi = h_ele(iloc)
          totmass = totmass + sum(detwei * mesh%nh%n(iloc,:)) * D(globi)
       end do
    end do

  end function get_totmass

  subroutine density_smoother(D,mass,lump_mass)
    real, dimension(n_dens*n_layers), intent(inout) :: D
    type(csr_matrix) :: mass
    real :: lump_mass(n_dens)
    !  locals
    integer :: i
    real    :: rhs(n_dens)

    KSPType :: ksp_type
    PCType :: pc_type

    print*, (minval(lump_mass))

    ksp_type = KSPGMRES
    pc_type = PCSOR
    
    do i=0,N_layers-1
       rhs=dot_product(lump_mass,D(i*N_Dens+1:(i+1)*N_Dens))
       call gallopede_solve(D(i*N_Dens+1:(i+1)*N_Dens)&
            ,mass, rhs, ksp_type, pc_type, 1.0e-14, 5000)
       
    end do

  end subroutine density_smoother

 subroutine get_density_relation(D,mesh,rel_d)
    real, dimension(n_dens*n_layers), intent(in) :: D
    type(dg_mesh) :: mesh
    real, dimension(n_dens), intent(out):: rel_d
    !  locals
    integer :: i, ele
    real    :: rhs(n_dens)

    integer, dimension(:), pointer :: X_ele, h_ele
    real, dimension(2,mesh%nm%loc) :: ele_X
    type(csr_matrix)::den_mat

    real, dimension(mesh%nu%ngi) :: detwei, d1locgi
    real, dimension(mesh%nh%loc,mesh%nh%loc):: Qhh

    KSPType :: ksp_type
    PCType :: pc_type

    ksp_type = KSPGMRES
    pc_type = PCSOR
   
    call allocate(den_mat,mesh%Mass_h%sparsity)

    call zero(den_mat)
    rhs=0.0

       do ele=1,n_elements

          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
          
          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)
          
          call transform_to_physical(ele_X, n=mesh%nm, m = mesh%nh, &
               detwei = detwei)
                
          d1locgi=matmul(D(h_ele),mesh%nh%n)
                
          rhs(h_ele)=rhs(h_ele)+shape_rhs(mesh%nh,detwei)

          Qhh=shape_shape(mesh%nh,mesh%nh,detwei*d1locgi)
          
          call addto(den_mat,h_ele,h_ele,Qhh)
          
       end do
       
       call gallopede_solve(rel_d,&
            den_mat, rhs, ksp_type, pc_type, 1.0e-30, 10000)
       
    
  end subroutine get_density_relation


  subroutine load_quad_d(D_out,bottom,u,mesh)
    implicit none
    real, dimension(n_dens*n_layers), intent(out):: D_out
    real, dimension(n_dens), intent(out):: bottom
    real, dimension(2*n_layers*n_vels), intent(out):: u
    type(dg_mesh) :: mesh
    real ,dimension(mesh%nh%loc) :: data 
    real ,dimension(mesh%nu%loc) :: u_data 
    integer :: io1
    ! locals
    integer :: ele, i
    integer, pointer, dimension(:) ::h_ele, X_ele, u_ele
    real, dimension (2,3) :: ele_X
    real, dimension(n_dens*n_layers) :: D_layer
    

    KSPType :: ksp_type
    PCType  :: pc_type

    ewrite(1,*) 'subroutine load_quad_d'
    
    pc_type = PCSOR  ! Preconditioner context
    ksp_type = KSPCG ! Krylov subspace context

    open(unit=2502, file='DQ_initial.dat', status='old', &
         iostat=io1)
    if(io1.ne.0) then
       write (0,*) 'Could not open file DQ_initial.dat for reading'
       stop
    end if
    open(unit=2503, file='bottomQ.dat', status='old', &
         iostat=io1)
    if(io1.ne.0) then
       write (0,*) 'Could not open file bottomQ.dat for reading'
       stop
    end if

    open(unit=2504, file='v_initial.dat', status='old', &
         iostat=io1)
    if(io1.ne.0) then
       write (0,*) 'Could not open file v_initial.dat for reading'
       stop
    end if

    do i=1,n_layers
       do ele=1,n_elements
          h_ele=>mesh%EVList_h((ELE-1)*mesh%nh%loc+1:ELE*mesh%nh%loc)
          u_ele=>mesh%EVList_u((ELE-1)*mesh%nu%loc+1:ELE*mesh%nu%loc)
          X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)
          read(unit=2502, iostat=io1, fmt=*) data
          d_layer(h_ele+(i-1)*n_dens)=data
          read(unit=2504, iostat=io1, fmt=*) u_data
          u(u_ele+2*(i-1)*n_vels)=u_data
          if(i==1) then
             read(unit=2503, iostat=io1, fmt=*) data
             bottom(h_ele)=data
          end if
       end do
       do ele=1,n_elements
          u_ele=>mesh%EVList_u((ELE-1)*mesh%nu%loc+1:ELE*mesh%nu%loc)
          read(unit=2504, iostat=io1, fmt=*) u_data
          u(u_ele+2*(i-1)*n_vels+n_vels)=u_data
       end do
    end do

    

    D_out=d_layer

    close(2502)
    close(2503)
    close(2504)

  end subroutine load_quad_d


  subroutine energy_find(D,u,bottom,mesh,rho)   
    real, dimension(n_layers*n_dens), intent(in) :: D
    real, dimension(n_dens), intent(in) :: bottom
    real, dimension(2*n_layers*n_vels), intent(in) :: u
    real, dimension(n_layers), intent(in) :: rho
    type(dg_mesh) :: Mesh
    integer, dimension(:), pointer :: u_ele, h_ele, X_ele
    real, dimension(2,mesh%nx%loc) :: ele_X
    real, dimension(mesh%nh%ngi):: detwei
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(2,mesh%nu%ngi) :: ulocgi
     real, dimension(mesh%nu%ngi) :: divulocgi
    real, dimension(mesh%nh%ngi) :: dlocgi, hlocgi, blocgi, djlocgi
     real, dimension(2,mesh%nh%ngi) :: gradhlocgi
    real :: energy, totmass(n_layers)
    integer :: ele, i, j, ios

      energy=0.0
      totmass=0.0


    element_loop: do ele = 1, N_elements
       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)

       call transform_to_physical(ele_X, n=mesh%nx, m = mesh%nh, &
            dm_t = dnh_t, detwei = detwei)
       call transform_to_physical(ele_X, n=mesh%nx, m = mesh%nu, &
            dm_t = dnu_t)

       do i=1,n_layers
          ulocgi(1,:)=matmul(u(2*(i-1)*n_vels+u_ele),mesh%nu%n)
          ulocgi(2,:)=matmul(u(2*(i-1)*n_vels+n_vels+u_ele),mesh%nu%n)
          hlocgi=matmul(bottom(h_ele),mesh%nh%n)
          dlocgi=matmul(d((i-1)*n_dens+h_ele),mesh%nh%n)
          gradhlocgi(1,:)=matmul(bottom(h_ele),dnh_t(:,:,1))
          gradhlocgi(2,:)=matmul(bottom(h_ele),dnh_t(:,:,2))
          divulocgi=matmul(u(2*(i-1)*n_vels+u_ele),dnu_t(:,:,1))+&
               matmul(u(2*(i-1)*n_vels+n_vels+u_ele),dnu_t(:,:,2))
          do j=i+1,n_layers
             hlocgi=hlocgi-matmul(d((j-1)*n_dens+h_ele),mesh%nh%n)
             gradhlocgi(1,:)=gradhlocgi(1,:)&
                  -matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,1))
             gradhlocgi(2,:)=gradhlocgi(2,:)&
                  -matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,2))
          end do
        
          energy=energy+0.5*sum(rho(i)*dlocgi*sum(ulocgi*ulocgi,1)*detwei)

          energy=energy+0.5*sum(g0*rho(i)*dlocgi*dlocgi*detwei)
          energy=energy-sum(g0*rho(i)*dlocgi*hlocgi*detwei)

          blocgi=sum(ulocgi*gradhlocgi,1)
          do j=i+1,n_layers
             ulocgi(1,:)=matmul(u(2*(j-1)*n_vels+u_ele),mesh%nu%n)
             ulocgi(2,:)=matmul(u(2*(j-1)*n_vels+n_vels+u_ele),mesh%nu%n)
             djlocgi=matmul(d((j-1)*n_dens+h_ele),mesh%nh%n)
             blocgi=blocgi+ulocgi(1,:)*matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,1))&
                  +djlocgi*matmul(u(2*(j-1)*n_vels+u_ele),dnu_t(:,:,1))
             blocgi=blocgi+ulocgi(2,:)*matmul(d((j-1)*n_dens+h_ele),dnh_t(:,:,2))&
                  +djlocgi*matmul(u(2*(j-1)*n_vels+n_vels+u_ele),dnu_t(:,:,2))
          end do

          energy=energy+0.5*sum(rho(i)*dlocgi*dlocgi*dlocgi**divulocgi*&
               divulocgi*detwei/6.0)
          energy=energy+0.5*sum(rho(i)*dlocgi*dlocgi*divulocgi*blocgi*detwei)
          energy=energy+0.5*sum(rho(i)*dlocgi*blocgi*blocgi*detwei)
          totmass(i)=totmass(i)+sum(detwei*dlocgi)

       end do
    end do element_loop

    if (t==0) then
       open(2550, file='energy.log', status='replace', iostat=ios)
    else
       open(2550, file='energy.log',position='append',&
            status='old', iostat=ios)
    end if

    if (ios .ne. 0) then
       ewrite(1,*) 'Cannont open file energy.log, terminating'
       stop
    end if

    if (t==0) then
       open(2551, file='mass.log', status='replace', iostat=ios)
    else
       open(2551, file='mass.log',position='append',&
            status='old', iostat=ios)
    end if

    if (ios .ne. 0) then
       ewrite(1,*) 'Cannont open file mass.log, terminating'
       stop
    end if

    write(2550,*) t, energy
    write(2551,*) t, totmass

    close(2550)
    close(2551)
       
  end subroutine energy_find


end module density_equation
