#include "fdebug.h"

module MLCM_operator
  
  use elements
  use sparse_tools
  use global_parameters_gallopede
  use transform_elements
  use data_structures
  use fetools
!  use gallopede_solvers
  use mesh_tools
  use multilayer_tools
  use fldebug
  use density_equation
  
  implicit none

  public :: assemble_MLCM_operator, get_vels_layer, mini_out, get_moms_layer,&
       pressure_operator_solve
       
  private
  
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"

  integer :: get_vel_count = 0
      
contains
  
  subroutine assemble_MLCM_operator( D,m,u,mesh, &
       bottom, rho, &
       mommat, &
       rhs1,rhs2)
    implicit none
    type(block_csr_matrix), intent(inout) :: mommat
    real, dimension(:), intent(out) :: rhs1, rhs2
    real, dimension(:), intent(in) :: m
    real, dimension(:), intent(in) :: u
    real, dimension(:), intent(in) :: D
    type(dg_mesh) :: mesh
    real, intent(in), dimension(N_Layers) :: rho
    real, dimension(n_dens), intent(in) :: bottom
    
    !locals
    integer :: ele,i,j,k,iloc,jloc,iX,im,jm,n_min
    real, dimension(n_dens) :: rel_d
    real, dimension(N_layers,mesh%nh%ngi) :: dlocgi,d0locgi
    real, dimension(2,N_layers,mesh%nu%ngi) :: graddlocgi, gradhlocgi
    real, dimension(2,mesh%nu%ngi) :: gradblocgi, v1locgi,v2locgi
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nu%ngi) :: detwei, rel_dlocgi



    ! Local element information
    integer, dimension(:), pointer :: u_ele, X_ele, h_ele, m_ele
    integer, dimension(mesh%nh_f%loc), target :: surface_h_lno
    integer, dimension(mesh%nu_f%loc), target :: surface_u_lno
    integer, dimension(mesh%nm_f%loc), target :: surface_m_lno
    integer, dimension(:), pointer :: u_ele_2=>null(),h_ele_2=>null()
    integer, dimension(:), pointer :: x_ele_2=>null(),m_ele_2=>null()
    real, dimension(n_layers,mesh%nu_f%ngi) :: djlocgi_f
    real, dimension(2,mesh%nm%loc) :: ele_X_2
    real, dimension(2,mesh%nm_f%loc) :: ele_Xf, ele_Xf_2
    integer, dimension(:), allocatable :: neigh
    integer, dimension(:), pointer :: bdyu=>null(), bdyu_2=>null()
    integer, dimension(:), pointer :: bdyh=>null(), bdyh_2=>null()
    integer, dimension(:), pointer :: bdym=>null(), bdym_2=>null()
    integer :: bdy_i, ele_2, ni, gi, nod, face
    real :: detwei_f(mesh%nu_f%ngi)
    real :: normal(2,mesh%nu_f%ngi), normal2(2,mesh%nu_f%ngi)
    real :: avg_cst, ell
    

    real, dimension(mesh%nu%loc) :: lumped_mass
    real, allocatable, dimension(:,:) :: Qu, Qu_f
    real, allocatable, dimension(:,:) :: Qum
    real, allocatable, dimension(:,:,:,:) :: QT, QT_T, QT_f,QT_Tf
    real, allocatable, dimension(:,:,:) :: QuT
    real, allocatable, dimension(:,:,:,:) :: QTum

    real, parameter :: mu=1.e6

    integer :: uloc,mloc, uloc_f,mloc_f, n_vl, n_mm
    type(element_type), pointer :: nm,nu,nm_f,nu_f
    logical :: log1
    real, dimension(mesh%nu%ngi) :: xlocgi, ylocgi


    ewrite(1,*)("subroutine assemble_MLCM_operator")
    print*, size(u),size(rhs1),n_layers*n_vels
    print*, sum(u(2*n_layers*n_vels+1:size(u)))


    if (size(u)==2*n_vels*n_layers) then
       print*, 'u is continuous', size(u), 2*n_vels*n_layers
       uloc=mesh%nu%loc
       uloc_f=mesh%nu_f%loc
       nu=>mesh%nu
       nu_f=>mesh%nu_f
       n_vl=n_vels
    else
       print*, 'u is discontinuous', size(u), 2*n_moms*n_layers
       uloc=mesh%nm%loc
       uloc_f=mesh%nm_f%loc
       nu=>mesh%nm
       nu_f=>mesh%nm_f
       n_vl=n_moms
    end if
    if (size(m)==2*n_vels*n_layers) then
       print*, 'm is continuous'
       mloc=mesh%nu%loc
       mloc_f=mesh%nu_f%loc
       nm=>mesh%nu
       nm_f=>mesh%nu_f
       n_mm=n_vels
    else
       print*, 'm is discontinuous'
       mloc=mesh%nm%loc
       mloc_f=mesh%nm_f%loc
       nm=>mesh%nm
       n_mm=n_moms
       nm_f=>mesh%nm_f
    end if
       allocate(Qu(uloc,uloc), QT(2,2,uloc,uloc),&
       QT_T(2,2,uloc,uloc),  QuT(2,uloc,uloc),&
       QTum(2,2,uloc,mloc), Qum(uloc,mloc)) 
       allocate(Qu_f(uloc_f,uloc_f),QT_f(2,2,uloc_f,uloc_f),&
            QT_Tf(2,2,uloc_f,uloc_f))

!    assert(maxval(mommat%val)<huge(0.0))

    call zero(mommat)
!    call get_density_relation(D,mesh,rel_d)


    log1=mommat%blocks(1)==2*N_Layers
    assert(log1)
    log1=mommat%blocks(2)==2*N_Layers
    assert(log1)

    ewrite(2,*)("ele loop")

    rhs1=0.
    rhs2=0.


    if (RIGID_LID) then
       n_min=2
    else 
       n_min=1
    end if

     ! Initial (and for meshes with only one shape of element, final) size
    ! of neigh.
    ewrite(2,*)("Allocating memory for neigh");
    allocate(neigh(row_length(mesh%bdy_list,1)))

    ele_loop: do ele = 1, n_elements

       call point_set(u,u_ele,mesh,ele)
       call point_set(m,m_ele,mesh,ele)
!       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
!       m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
!       allocate(h_ele(3))
!       h_ele=mesh%EVList_h((ELE-1)*mesh%Nh%LOC+(/1,3,6/))
       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)

       call transform_to_physical(ele_X, mesh%nm, m = mesh%nh, &
            dn_t = dnm_t, dm_t = dnh_t, &
            detwei = detwei)
       call transform_to_physical(ele_X, mesh%nm, m = mesh%nu, &
            dm_t = dnu_t)

!       do i=1,uloc
!          lumped_mass(i)=sum(row_val_ptr(mesh%mass_u,u_ele(i)))
!       end do

       !get dlocgi
       dlocgi = 0.

!       rel_Dlocgi=matmul(rel_D(h_ele),mesh%nh%n)

       do i = 1, N_Layers
          dlocgi(i,:) = matmul(D((i-1)*n_dens+h_ele),mesh%nh%n)
           D0locgi(i,:)=sum(D((I-1)*n_dens+1:I*n_dens))/n_dens
       end do

       assert(maxval(dlocgi)<huge(0.0))

       !get gradblocgi
       gradblocgi = 0.
       do iX = 1,2
          gradblocgi(iX,:) = matmul(bottom(h_ele),dnh_t(:,:,iX))
       end do

       assert(maxval(gradblocgi)<huge(0.0))

       !get graddlocgi
       graddlocgi = 0.
       forall(i = 1:N_Layers, iX=1:2)
          graddlocgi(iX,i,:) = matmul(D((i-1)*n_dens+h_ele),dnh_t(:,:,iX))
       end forall

       assert(maxval(graddlocgi)<huge(0.0))

       !get gradhlocgi
       gradhlocgi = 0.
       do i = 1, N_Layers
           gradhlocgi(:,i,:) = gradblocgi(:,:)
           do j = i+1,n_layers
          gradhlocgi(:,i,:) = gradhlocgi(:,i,:) - &
               graddlocgi(:,j,:)
          end do
       end do

       assert(maxval(gradhlocgi)<huge(0.0))

       !RIGHT HAND SIDE

       !get local mass matrix for RHS
 

 

       do i = 1,N_Layers
      Qum=shape_shape(nu,nm,detwei*dlocgi(i,:))
          Qu=shape_shape(nu,nu,&
               detwei *dlocgi(i,:))
!      Qum=shape_shape(mesh%nu,mesh%nm,detwei)

      assert(maxval(Qum)<huge(0.0))

      rhs1((i-1)*n_vl+u_ele) = rhs1((i-1)*n_vl+u_ele) + &
           matmul(Qum,m((i-1)*2*n_moms+m_ele))
      rhs2((i-1)*n_vl+u_ele) = rhs2((i-1)*n_vl+u_ele) + &
           matmul(Qum,m((i-1)*2*n_mm+n_mm+m_ele))
      end do

       !LEFT HAND SIDE

       !horizontal kinetic energy

       do i = 1, N_Layers
          Qu=shape_shape(nu,nu,&
               detwei *dlocgi(i,:)*rho(i))
!               detwei*rho(i))
          assert(maxval(Qu)<huge(0.0))

          call addto(mommat,2*(i-1)+1,2*(i-1)+1, &
               u_ele,u_ele,Qu)
          call addto(mommat,2*(i-1)+2,2*(i-1)+2, &
               u_ele,u_ele,Qu)
       end do

       !vertical kinetic energy
       
       ! sum_i \rho_i <w_i.grad h_i,D_i u_i.grad h_i>

       if(.not. GN_FLAG) then
          !gradient-free terms
          do i = n_min, N_Layers
             QT=shape_shape_vector_outer_vector( &
                  nu,nu,detwei*dlocgi(i,:)*rho(i), &
                  gradhlocgi(:,i,:),gradhlocgi(:,i,:))
             assert(maxval(QT)<huge(0.0))
!             call tensor_addto(mommat,i,i,u_ele,u_ele,QT)

             forall (im=1:2,jm=1:2)
                QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
             end forall

             call tensor_addto(mommat,i,i,u_ele,u_ele,0.5*(QT_T+QT))
          end do
       end if

       if(.not. GN_FLAG) then

          ! sum_i rho_i<w_i.grad h_i, D_i sum div (D_ju_j) >
          ! sum_i rho_i<u_i.grad h_i, D_i sum div (D_jw_j) >

          do i = n_min, N_Layers
             if(i<N_Layers) then
                do j = i+1, N_Layers

                   QT=shape_shape_vector_outer_vector( &
                        nu,nu,detwei*dlocgi(i,:)*rho(i), &
                        gradhlocgi(:,i,:),graddlocgi(:,j,:))
!                   call tensor_addto(mommat,i,j,u_ele,u_ele,QT)
!                   call tensor_addto(mommat,j,i,u_ele,u_ele,QT,.true.)

                   QT=QT+shape_vector_outer_dshape( &
                        nu,gradhlocgi(:,i,:),nu%dn,&
                        detwei*dlocgi(i,:)*dlocgi(j,:)*rho(i))

                     forall (im=1:2,jm=1:2)
                        QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
                     end forall

                   call tensor_addto(mommat,i,j,u_ele,u_ele,QT)
                   call tensor_addto(mommat,j,i,u_ele,u_ele,QT_T)
                end do
             end if
          end do

       end if

       if(.not. GN_FLAG) then
          ! sum_i<D_i sum_j div(D_jw_j), sum_k div(D_ku_k) >rho_i
          do i = n_min, N_Layers
             if(i<N_Layers) then
                do j = i+1, N_Layers
                   do k = j, N_Layers

                      QT=shape_shape_vector_outer_vector( &
                           nu,nu,detwei*dlocgi(i,:)*rho(i), &
                           graddlocgi(:,j,:),graddlocgi(:,k,:))

                      forall (im=1:2,jm=1:2)
                         QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
                      end forall

                      call tensor_addto(mommat,j,k,u_ele,u_ele,QT)
                      if (k>j) &
                           call tensor_addto(mommat,k,j,u_ele,u_ele,QT_T)


                      QT=shape_vector_outer_dshape( &
                           nu,graddlocgi(:,j,:),nu%dn,&
                           detwei*dlocgi(i,:)*dlocgi(k,:)*rho(i))

!                      call tensor_addto(mommat,j,k,u_ele,u_ele,QT)
!                      call tensor_addto(mommat,k,j,u_ele,u_ele,QT,.true.)
!                      call tensor_addto(mommat,j,k,u_ele,u_ele,(QT_T+QT))
!                      call tensor_addto(mommat,k,j,u_ele,u_ele,QT_T)
                      

                      QT=dshape_outer_vector_shape( &
                           nu%dn,graddlocgi(:,k,:),nu,&
                           detwei*dlocgi(i,:)*dlocgi(j,:)*rho(i))

                      forall (im=1:2,jm=1:2)
                         QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
                      end forall

                      call tensor_addto(mommat,j,k,u_ele,u_ele,QT+QT_T)
                      if (k>j) &
                           call tensor_addto(mommat,k,j,u_ele,u_ele,QT+QT_T)

                      QT =dshape_outer_dshape(nu%dn,nu%dn, &
                           detwei*dlocgi(i,:)*dlocgi(j,:)*dlocgi(k,:)*rho(i))

                      forall (im=1:2,jm=1:2)
                         QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
                      end forall


                      call tensor_addto(mommat,j,k,u_ele,u_ele,QT)
                      if (k>j) &
                           call tensor_addto(mommat,k,j,u_ele,u_ele,QT_T)

                   end do
                end do
             end if
          end do
       end if

          ! sum_i< D_iD_i div w_i,u_i.grad h_i + sum_j div(D_ju_j)>rho_i/2
          ! sum_i< w_i.grad h_i + sum_j div(D_jw_j), D_iD_i div u_i>rho_i/2

       if(.not. GN_FLAG) then
          do i = n_min, N_Layers

!             print*, i, gradhlocgi(:,i,:)

             QT=dshape_outer_vector_shape( &
                  nu%dn,gradhlocgi(:,i,:),nu, &
                  detwei*dlocgi(i,:)*dlocgi(i,:)*rho(i)/2.0)

             forall (im=1:2,jm=1:2)
                QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
             end forall

             call tensor_addto(mommat,i,i,u_ele,u_ele,QT_T+QT)

   

             if(i<N_Layers) then

                do j = i+1, N_Layers

                   QT =dshape_outer_vector_shape( &
                        nu%dn,graddlocgi(:,j,:),nu, &
                        detwei*dlocgi(i,:)*dlocgi(i,:)*rho(i)/2.0)
!                  call tensor_addto(mommat,i,j,u_ele,u_ele,QT)
!                  call tensor_addto(mommat,j,i,u_ele,u_ele,QT,.true.)

                   QT = QT+dshape_outer_dshape(nu%dn,nu%dn, &
                        detwei*dlocgi(i,:)*dlocgi(i,:)*dlocgi(j,:)* &
                        rho(i)/2.0)

                   forall (im=1:2,jm=1:2)
                      QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
                   end forall

                   call tensor_addto(mommat,i,j,u_ele,u_ele,QT)
                   call tensor_addto(mommat,j,i,u_ele,u_ele,QT_T)

                end do

             end if
          end do

       end if

          ! sum_i< D_iD_iD_i div w_i, div u_i>rho_i/3
          do i = n_min, N_Layers

              if(.true.) then

             QT = dshape_outer_dshape(nu%dn,nu%dn, &
!                detwei*dlocgi(i,:)*&
!                  detwei*dlocgi(i,:)*dlocgi(i,:)*dlocgi(i,:)* &
!                  detwei*dlocgi(i,:)*d0locgi(i,:)*d0locgi(i,:)* &
              detwei*dlocgi(i,:)*dlocgi(i,:)*&
                  rho(i)/3.0)

             forall (im=1:2,jm=1:2)
                QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
             end forall


             call tensor_addto(mommat,i,i,u_ele,u_ele,0.5*(QT+QT_T))
!             call tensor_addto(mommat,i,i,u_ele,u_ele,0.5*QT)
!             call tensor_addto(mommat,i,i,u_ele,u_ele,0.5*QT,.true.)


          end if

             QU = dshape_dot_dshape(nu%dn,nu%dn, &
                  detwei*rho(i)*alpha*alpha)
             

             call addto(mommat,2*(i-1)+1,2*(i-1)+1, &
                  u_ele,u_ele,Qu)
             call addto(mommat,2*(i-1)+2,2*(i-1)+2, &
                  u_ele,u_ele,Qu)

          end do

          if (RIGID_LID) then
             
             do i=n_min,n_layers

                do j=n_min,n_layers

                if (.False.) then
             QT = dshape_outer_dshape(nu%dn,nu%dn, &
                  detwei*dlocgi(1,:)*dlocgi(i,:)*dlocgi(j,:))/3.0
             forall (im=1:2,jm=1:2)
                QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
             end forall
             call tensor_addto(mommat,i,j,u_ele,u_ele,0.5*(QT+QT_T))
   
             v1locgi(1,:)= dlocgi(1,:)*graddlocgi(1,i,:)&
                  - dlocgi(i,:)*graddlocgi(1,1,:)
             v1locgi(2,:)= dlocgi(1,:)*graddlocgi(2,i,:)&
                  - dlocgi(i,:)*graddlocgi(2,1,:)
             v2locgi(1,:)= dlocgi(1,:)*graddlocgi(1,j,:)&
                  - dlocgi(j,:)*graddlocgi(1,1,:)
             v2locgi(2,:)= dlocgi(1,:)*graddlocgi(2,j,:)&
                  - dlocgi(j,:)*graddlocgi(2,1,:)

             QT= dshape_outer_vector_shape(nu%dn,&
                  v1locgi,&
                  nu,detwei*dlocgi(j,:))/3.0
             forall (im=1:2,jm=1:2)
                QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
             end forall

             call tensor_addto(mommat,i,j,u_ele,u_ele,QT_T+QT)

             QT= shape_shape_vector_outer_vector(nu,&
                  nu, detwei*rel_dlocgi,&
                 v1locgi,v2locgi)/3.0
             forall (im=1:2,jm=1:2)
                QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
             end forall

             call tensor_addto(mommat,i,j,u_ele,u_ele,0.5*(QT_T+QT))

          end if

       end do

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
       if (size(u)==2*n_moms*n_layers) then
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
                     mesh%nm%numbering)
                bdyu => surface_m_lno
                if (size(m)==2*n_vels*n_layers) then
                   surface_m_lno = boundary_local_num(offnods(ni,:), &
                        mesh%nu%numbering)
                   bdym => surface_u_lno
                else
                   surface_m_lno = boundary_local_num(offnods(ni,:), &
                        mesh%nm%numbering)
                   bdym => surface_m_lno
                end if
             end if
          else

             avg_cst = 0.5
             call point_set(u,u_ele_2,mesh,ele_2)
             call point_set(m,m_ele_2,mesh,ele_2)
!             u_ele_2=>mesh%EVList_u((ELE_2-1)*mesh%nu%LOC+1:ELE_2*mesh%nu%LOC)
!             m_ele_2=>mesh%EVList_m((ELE_2-1)*mesh%nm%LOC+1:ELE_2*mesh%nm%LOC)
             h_ele_2=>mesh%EVList_h((ELE_2-1)*mesh%Nh%LOC+1:ELE_2*mesh%Nh%LOC)
             X_ele_2=>mesh%EVList_X((ELE_2-1)*3+1:ELE_2*3)



             bdy_i=ival(mesh%bdy_list, ele, ele_2)
             bdyh=>mesh%bdy_nh_lno( &
                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)
             if (size(m)==2*n_vels*n_layers) then
                bdym=> mesh%bdy_nu_lno( &
                  (bdy_i-1)*mesh%nu_f%loc+1:bdy_i*mesh%nu_f%loc)
             else
                bdym=>mesh%bdy_nm_lno( &
                  (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)
             end if
             bdyu=>mesh%bdy_nm_lno(&
                  (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)

             bdy_i=ival(mesh%bdy_list, ele_2, ele)
             bdyu_2=> mesh%bdy_nm_lno( &
                  (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)
             bdyh_2=>mesh%bdy_nh_lno( &
                  (bdy_i-1)*mesh%nh_f%loc+1:bdy_i*mesh%nh_f%loc)
             if (size(m)==2*n_vels*n_layers) then
                bdym_2=> mesh%bdy_nu_lno( &
                     (bdy_i-1)*mesh%nu_f%loc+1:bdy_i*mesh%nu_f%loc)
             else
                bdym_2=>mesh%bdy_nm_lno( &
                     (bdy_i-1)*mesh%nm_f%loc+1:bdy_i*mesh%nm_f%loc)
             end if


             ! Locations of local vertices.
             ele_X_2(1,:)=mesh%X(X_ele_2)
             ele_X_2(2,:)=mesh%Y(X_ele_2)

          ! Locations of bdy vertices.

          ele_Xf=ele_X(:,bdyu)
          ele_Xf_2=ele_X_2(:,bdyu_2)

          do i = 1, N_Layers
             djlocgi_f(i,:) = matmul(D((i-1)*n_dens+h_ele(bdyh)),mesh%nh_f%n)
          end do

          ! Change of coordinates on face
          call transform_bdy_to_physical(ele_X, ele_Xf, nu, nu_f, &
               detwei_f = detwei_f,normal = normal)
          call transform_bdy_to_physical(ele_X_2, ele_Xf_2, nu,nu_f, &
               normal = normal2)

          assert(all(mod(ele_Xf,maxval(mesh%X))==mod(ele_XF_2,maxval(mesh%X))))

          !==========================================================
          !Helmholtz terms
          face=get_face(bdyu)

          do i=n_min,n_layers

          Qu_f=-dshape_dot_vector_shape(nu%dn_s(bdyu,:,face,:),&
               normal,nu_f,alpha*detwei_f)
          QT_f=-shape_vector_outer_dshape(nu_f,&
               normal,nu%dn_s(bdyu,:,face,:),&
!               djlocgi_f(i,:)*djlocgi_f(i,:)*&
               djlocgi_f(i,:)*detwei_f*rho(i))
           forall (im=1:2,jm=1:2)
              QT_Tf(im,jm,:,:)=transpose(QT_f(jm,im,:,:))
           end forall

         call tensor_addto(mommat,i,i,u_ele(bdyu),u_ele(bdyu),0.5*(QT_f+QT_Tf))
!         call tensor_addto(mommat,i,i,u_ele_2(bdyu_2),u_ele(bdyu),0.5*(QT_Tf))
!          call tensor_addto(mommat,i,i,u_ele(bdyu),u_ele_2(bdyu_2),0.5*(QT_f))

          if (.not. GN_FLAG) then
             do j = i+1, N_Layers
!                QT_F=-shape_vector_outer_vector_shape(nu_f,&
!                     normal,gradhlocgi_f(:,i,:),nu_f,djlocgi_f(i,:)*&
!                     djlocgi_f(j,:)*rho(i)*detwei_f)
             end do
          end if

          call addto(mommat,2*(i-1)+1,2*(i-1)+1,&
               u_ele(bdyu),u_ele(bdyu),0.5*(Qu_f+transpose(Qu_f)))
           call addto(mommat,2*(i-1)+1,2*(i-1)+1,&
               u_ele(bdyu),u_ele_2(bdyu_2),0.5*(Qu_f))
           call addto(mommat,2*(i-1)+1,2*(i-1)+1,&
               u_ele_2(bdyu_2),u_ele(bdyu),0.5*(transpose(Qu_f)))
          call addto(mommat,2*(i-1)+2,2*(i-1)+2,&
               u_ele(bdyu),u_ele(bdyu),0.5*(Qu_f+transpose(Qu_f)))
          call addto(mommat,2*(i-1)+2,2*(i-1)+2,&
               u_ele(bdyu),u_ele_2(bdyu_2),0.5*(Qu_f))
          call addto(mommat,2*(i-1)+2,2*(i-1)+2,&
               u_ele_2(bdyu_2),u_ele(bdyu),0.5*(transpose(Qu_f)))

          !Penalty term
          ell=1.0/sqrt(sum((ele_Xf(:,1)-ele_Xf(:,2))**2))
          Qu_f=shape_shape(nu_f,nu_f,mu*ell*detwei_f)

          call addto(mommat,2*(i-1)+1,2*(i-1)+1,&
               u_ele(bdyu),u_ele(bdyu),0.25*(Qu_f))
           call addto(mommat,2*(i-1)+1,2*(i-1)+1,&
               u_ele_2(bdyu_2),u_ele(bdyu),-0.25*(Qu_f))
           call addto(mommat,2*(i-1)+2,2*(i-1)+2,&
               u_ele(bdyu),u_ele(bdyu),0.25*(Qu_f))
           call addto(mommat,2*(i-1)+2,2*(i-1)+2,&
               u_ele_2(bdyu_2),u_ele(bdyu),-0.25*(Qu_f))

           end do


          end if
       end do neighbourloop2

    end if

end do ele_loop

    ewrite(1,*)("end subroutine assemble_MLCM_operator")

  end subroutine assemble_MLCM_operator

  subroutine get_vels_layer(Mesh, bcs, D, bottom, rho, u, m, u0,pout)
    type(dg_mesh), intent(in) :: Mesh
    real, dimension(n_layers*n_dens), intent(in) :: D
    real, dimension(n_dens), intent(in) :: bottom
    real, dimension(n_layers), intent(in) :: rho
    type(bc_info), intent(in) :: bcs
    real, dimension(:), intent(inout) :: u
    real, dimension(:), intent(in), optional :: u0
    real, dimension(:), intent(in) :: m
    logical, optional :: pout

    !locals
    PetscErrorCode :: ierr
    KSPType :: ksp_type
    PCType :: pc_type

    real, dimension(size(u)):: urhs, u_check
    type(block_csr_matrix) :: mommat
    integer :: i, j, u_it

    assert(size(D) == n_layers*n_dens)
    assert(size(bottom) == n_dens)
    assert(size(rho) == n_layers)
    if (present(u0)) then
        assert(size(u0) == size(u))
    end if

     if (constant_u .and. (.not. present(pout))) then
        call constant_u_set(u,d,mesh)
     else

    ewrite(1,*)("subroutine get_vels_layers")



    get_vel_count = get_vel_count + 1

    urhs = 0.

    if (size(u)==2*n_vels*n_layers) then
       call allocate(mommat,mesh%Mass_u%sparsity,&
            (/2 * N_layers, 2 * N_layers /))
    else
        call allocate(mommat,mesh%mass_h%sparsity,&
            (/2 * N_layers, 2 * N_layers /))
     end if
    
    assert(maxval(D)<huge(0.0))


    if(present(u0)) then
       u = u0
    else
       u = 0.
    end if

    call assemble_MLCM_operator(D,m,u,mesh,bottom,rho, &
         mommat=mommat, &
         rhs1=urhs(1:size(urhs)/2),&
         rhs2=urhs(size(urhs)/2+1:size(urhs)))

    ewrite(2,*) huge(0.0)
!    assert(maxval(mommat%val)<huge(0.0))

    ewrite(2,*)("solving GN operator equation")
    !solve GN equation
    ksp_type = KSPCG
    pc_type = PCICC

!   ksp_type = KSPGMRES
!    pc_type = PCSOR


!    ewrite(1,*) 'maxval(mommat)', maxval(mommat%val)
    

!    ewrite(1,*) 'maxval(mommat)=', maxval(mommat%val)
!    assert(maxval(mommat%val)<huge(0.0))

    print*, bcs%N_tangents

!    call    mini_out(urhs(1:n_layers*n_vels),mesh%EVList_u,&
!         'urhs1.dat',mesh%nu%loc,n_layers)
!    call    mini_out(urhs(n_layers*n_vels+1:2*n_layers*n_vels),mesh%EVList_u,&
!         'urhs2.dat',mesh%nu%loc,n_layers)

    if(bcs%n_tangents .ne. 0) then
       print*, 'boundary solve'
       call gallopede_block_solve_lift(u,&
            urhs(1:size(urhs)/2),&
            urhs(size(urhs)/2+1:size(urhs)), &
            mommat,&
            bcs,ksp_type, pc_type,&
            rerrtol=1.d-20,abserrtol=1.d-12,noi= 10000,zero=.false.)
    else
       print*, "No boundaries"
       call swap_rhs(urhs)
          u_check=u
          call gallopede_block_solve(u,mommat,&
               urhs,ksp_type, pc_type,&
               errtol=1.d-32,noi= 30000,zero=.false.)
!           call u_smooth(u,m,d,mesh)
       end if
    call deallocate( mommat )

    if (RIGID_LID) call get_u1_rigid_lid(D,u,mesh)
 end if

    ewrite(1,*)("END subroutine get_vels_layers")

  end subroutine get_vels_layer

  subroutine assemble_forward_MLCM_operator( D,m,u,mesh, &
       bottom, rho, &
       mommat1,mommat2, &
       rhs1,rhs2)
    implicit none
    type(csr_matrix), dimension (n_layers), intent(inout) :: mommat1
    type(csr_matrix), intent(inout) ::mommat2
    real, dimension(n_layers*n_dens), intent(out) :: rhs1, rhs2
    real, dimension(2*n_layers*n_moms), intent(in) :: m
    real, dimension(2*n_layers*n_vels), intent(in) :: u
    real, dimension(n_layers*n_dens), intent(in) :: D
    type(dg_mesh) :: mesh
    real, intent(in), dimension(N_Layers) :: rho
    real, dimension(n_dens), intent(in) :: bottom
    
    !locals
    integer :: ele,i,j,k,iloc,jloc,iX,im,jm
    real, dimension(N_layers,mesh%nh%ngi) :: dlocgi
    real, dimension(2,N_layers,mesh%nu%ngi) :: graddlocgi, gradhlocgi
    real, dimension(2,mesh%nu%ngi) :: gradblocgi
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nu%ngi) :: detwei

    ! Local element information
    integer, dimension(:), pointer :: u_ele, X_ele, h_ele, m_ele

    real, dimension(mesh%nu%loc) :: lumped_mass
    real, dimension(mesh%nh%loc,mesh%nu%loc) :: Qu
    real, dimension(mesh%nh%loc,mesh%nh%loc) :: Qhh
    real, dimension(mesh%nm%loc,mesh%nm%loc) :: Qmm
    real, dimension(2,2,mesh%nh%loc,mesh%nu%loc) :: QT, QT_T
 
    logical :: log1
    real, dimension(mesh%nu%ngi) :: xlocgi, ylocgi

    ewrite(2,*)("subroutine assemble_forward_MLCM_operator")

!    assert(maxval(mommat%val)<huge(0.0))
    assert(size(D)==n_layers*n_dens)
    assert(size(m)==2*n_layers*n_moms)
    assert(size(rhs1)==n_layers*n_dens)
    assert(size(rhs2)==n_layers*n_dens)
    assert(size(bottom)==n_dens)


    do i=1,n_layers
       call zero(mommat1(i))
    end do
    call zero(mommat2)

    ewrite(2,*)("ele loop")

    rhs1=0.
    rhs2=0.

    ele_loop: do ele = 1, n_elements

       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)

       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)

       call transform_to_physical(ele_X, mesh%nm, m = mesh%nh, &
            dn_t = dnm_t, dm_t = dnh_t, &
            detwei = detwei)
       call transform_to_physical(ele_X, mesh%nm, m = mesh%nu, &
            dm_t = dnu_t)

       !get dlocgi
       dlocgi = 0.
       do i = 1, N_Layers
          dlocgi(i,:) = matmul(D((i-1)*n_dens+h_ele),mesh%nh%n)
       end do

       assert(maxval(dlocgi)<huge(0.0))

       !get gradblocgi
       gradblocgi = 0.
       do iX = 1,2
          gradblocgi(iX,:) = matmul(bottom(h_ele),dnh_t(:,:,iX))
       end do

       assert(maxval(gradblocgi)<huge(0.0))

       !get graddlocgi
       graddlocgi = 0.
       forall(i = 1:N_Layers, iX=1:2)
          graddlocgi(iX,i,:) = matmul(D((i-1)*n_dens+h_ele),dnh_t(:,:,iX))
       end forall

       assert(maxval(graddlocgi)<huge(0.0))

       !get gradhlocgi
       gradhlocgi = 0.
       do i = 1, N_Layers
          gradhlocgi(:,i,:) = gradblocgi(:,:) - &
               sum(graddlocgi(:,i+1:N_layers,:),2)
       end do

       assert(maxval(gradhlocgi)<huge(0.0))

       !LHS HAND SIDE

       !get local mass matrix


       do i=1,n_layers
          Qhh=shape_shape(mesh%nh,mesh%nh,detwei*dlocgi(i,:))
          call addto(mommat1(i),h_ele,h_ele,Qhh)
       end do
       Qmm=shape_shape(mesh%nm,mesh%nm,detwei)
       call addto(mommat2,m_ele,m_ele,Qmm)


       !RHS HAND SIDE

       !horizontal kinetic energy

       do i = 1, N_Layers
          Qu=shape_shape(mesh%nh,mesh%nu,&
               detwei *dlocgi(i,:)*rho(i))
          assert(maxval(Qu)<huge(0.0))

          rhs1((i-1)*n_dens+h_ele)=&
               rhs1((i-1)*n_dens+h_ele)+&
               matmul(Qu,u(2*(i-1)*n_vels+u_ele))
          rhs2((i-1)*n_dens+h_ele)=&
               rhs2((i-1)*n_dens+h_ele)+&
               matmul(Qu,u(2*(i-1)*n_vels+n_vels+u_ele))

       end do

       !vertical kinetic energy
       
       ! sum_i \rho_i <w_i.grad h_i,D_i u_i.grad h_i>

       if(.not. GN_FLAG) then
          !gradient-free terms
          do i = 1, N_Layers
             QT=shape_shape_vector_outer_vector( &
                  mesh%nh,mesh%nu,detwei*dlocgi(i,:)*rho(i), &
                  gradhlocgi(:,i,:),gradhlocgi(:,i,:))
             assert(maxval(QT)<huge(0.0))

             rhs1((i-1)*n_dens+h_ele)=&
                  rhs1((i-1)*n_dens+h_ele)+&
                  matmul(QT(1,1,:,:),u(2*(i-1)*n_vels+u_ele))
             rhs1((i-1)*n_dens+h_ele)=&
                  rhs1((i-1)*n_dens+h_ele)+&
                  matmul(QT(1,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))
             rhs2((i-1)*n_dens+h_ele)=&
                  rhs2((i-1)*n_dens+h_ele)+&
                  matmul(QT(2,1,:,:),u(2*(i-1)*n_vels+u_ele))
             rhs2((i-1)*n_dens+h_ele)=&
                  rhs2((i-1)*n_dens+h_ele)+&
                  matmul(QT(2,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))
          end do
       end if

       if(.not. GN_FLAG) then

          ! sum_i rho_i<w_i.grad h_i, D_i sum div (D_ju_j) >
          ! sum_i rho_i<u_i.grad h_i, D_i sum div (D_jw_j) >

          do i = 1, N_Layers
             if(i<N_Layers) then
                do j = i+1, N_Layers

                   QT=shape_shape_vector_outer_vector( &
                        mesh%nh,mesh%nu,detwei*dlocgi(i,:)*rho(i), &
                        gradhlocgi(:,i,:),graddlocgi(:,j,:))

                   rhs1((i-1)*n_dens+h_ele)=&
                        rhs1((i-1)*n_dens+h_ele)+&
                        matmul(QT(1,1,:,:),u(2*(j-1)*n_vels+u_ele))
                   rhs1((i-1)*n_dens+h_ele)=&
                        rhs1((i-1)*n_dens+h_ele)+&
                        matmul(QT(1,2,:,:),u(2*(j-1)*n_vels+n_vels+u_ele))
                   rhs2((i-1)*n_dens+h_ele)=&
                        rhs2((i-1)*n_dens+h_ele)+&
                        matmul(QT(2,1,:,:),u(2*(j-1)*n_vels+u_ele))
                   rhs2((i-1)*n_dens+h_ele)=&
                        rhs2((i-1)*n_dens+h_ele)+&
                        matmul(QT(2,2,:,:),u(2*(j-1)*n_vels+n_vels+u_ele))

                   QT=shape_shape_vector_outer_vector( &
                        mesh%nh,mesh%nu,detwei*dlocgi(i,:)*rho(i), &
                        graddlocgi(:,j,:),gradhlocgi(:,i,:))


                   rhs1((j-1)*n_dens+h_ele)=&
                        rhs1((j-1)*n_dens+h_ele)+&
                        matmul(QT(1,1,:,:),u(2*(i-1)*n_vels+u_ele))
                   rhs1((j-1)*n_dens+h_ele)=&
                        rhs1((j-1)*n_dens+h_ele)+&
                        matmul(QT(1,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))
                   rhs2((j-1)*n_dens+h_ele)=&
                        rhs2((j-1)*n_dens+h_ele)+&
                        matmul(QT(2,1,:,:),u(2*(i-1)*n_vels+u_ele))
                   rhs2((j-1)*n_dens+h_ele)=&
                        rhs2((j-1)*n_dens+h_ele)+&
                        matmul(QT(2,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))


                   QT=shape_vector_outer_dshape( &
                        mesh%nh,gradhlocgi(:,i,:),dnu_t,&
                        detwei*dlocgi(i,:)*dlocgi(j,:)*rho(i))

                   rhs1(2*(i-1)*n_dens+h_ele)=&
                        rhs1(2*(i-1)*n_dens+h_ele)+&
                        matmul(QT(1,1,:,:),u(2*(j-1)*n_vels+u_ele))
                   rhs1((i-1)*n_dens+h_ele)=&
                        rhs1((i-1)*n_dens+h_ele)+&
                        matmul(QT(1,2,:,:),u(2*(j-1)*n_vels+n_vels+u_ele))
                   rhs2((i-1)*n_dens+h_ele)=&
                        rhs2((i-1)*n_dens+h_ele)+&
                        matmul(QT(2,1,:,:),u(2*(j-1)*n_vels+u_ele))
                   rhs2((i-1)*n_dens+h_ele)=&
                        rhs2((i-1)*n_dens+h_ele)+&
                        matmul(QT(2,2,:,:),u(2*(j-1)*n_vels+n_vels+u_ele))

                   QT=dshape_outer_vector_shape( &
                        dnh_t,gradhlocgi(:,i,:),mesh%nu,&
                        detwei*dlocgi(j,:)*dlocgi(i,:)*rho(i))

                   rhs1((j-1)*n_dens+h_ele)=&
                        rhs1((j-1)*n_dens+h_ele)+&
                        matmul(QT(1,1,:,:),u(2*(i-1)*n_vels+u_ele))
                   rhs1((j-1)*n_dens+h_ele)=&
                        rhs1((j-1)*n_dens+h_ele)+&
                        matmul(QT(1,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))
                   rhs2((j-1)*n_dens+h_ele)=&
                        rhs2((j-1)*n_dens+h_ele)+&
                        matmul(QT(2,1,:,:),u(2*(i-1)*n_vels+u_ele))
                   rhs2((j-1)*n_dens+h_ele)=&
                        rhs2((j-1)*n_dens+h_ele)+&
                        matmul(QT(2,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))

                end do
             end if
          end do

       end if

       if(.not. GN_FLAG) then
          ! sum_i<D_i sum_j div(D_jw_j), sum_k div(D_ku_k) >rho_i
          do i = 1, N_Layers
             if(i<N_Layers) then
                do j = i+1, N_Layers
                   do k = i+1, N_Layers

                      QT=shape_shape_vector_outer_vector( &
                           mesh%nh,mesh%nu,detwei*dlocgi(i,:)*rho(i), &
                           graddlocgi(:,j,:),graddlocgi(:,k,:))


                   rhs1((j-1)*n_dens+h_ele)=&
                        rhs1((j-1)*n_dens+h_ele)+&
                        matmul(QT(1,1,:,:),u(2*(k-1)*n_vels+u_ele))
                   rhs1((j-1)*n_dens+h_ele)=&
                        rhs1((j-1)*n_dens+h_ele)+&
                        matmul(QT(1,2,:,:),u(2*(k-1)*n_vels+n_vels+u_ele))
                   rhs2((j-1)*n_dens+h_ele)=&
                        rhs2((j-1)*n_dens+h_ele)+&
                        matmul(QT(2,1,:,:),u(2*(k-1)*n_vels+u_ele))
                   rhs2((j-1)*n_dens+h_ele)=&
                        rhs2((j-1)*n_dens+h_ele)+&
                        matmul(QT(2,2,:,:),u(2*(k-1)*n_vels+n_vels+u_ele))


                      QT=shape_vector_outer_dshape( &
                           mesh%nh,graddlocgi(:,j,:),dnu_t,&
                           detwei*dlocgi(i,:)*dlocgi(k,:)*rho(i))

                   rhs1((j-1)*n_dens+h_ele)=&
                        rhs1((j-1)*n_dens+h_ele)+&
                        matmul(QT(1,1,:,:),u(2*(k-1)*n_vels+u_ele))
                   rhs1((j-1)*n_dens+h_ele)=&
                        rhs1((j-1)*n_dens+h_ele)+&
                        matmul(QT(1,2,:,:),u(2*(k-1)*n_vels+n_vels+u_ele))
                   rhs2((j-1)*n_dens+h_ele)=&
                        rhs2((j-1)*n_dens+h_ele)+&
                        matmul(QT(2,1,:,:),u(2*(k-1)*n_vels+u_ele))
                   rhs2((j-1)*n_dens+h_ele)=&
                        rhs2((j-1)*n_dens+h_ele)+&
                        matmul(QT(2,2,:,:),u(2*(k-1)*n_vels+n_vels+u_ele))

                   QT=dshape_outer_vector_shape( &
                        dnh_t,graddlocgi(:,k,:),mesh%nu,&
                        detwei*dlocgi(i,:)*dlocgi(j,:)*rho(i))

                   rhs1((k-1)*n_dens+h_ele)=&
                        rhs1((k-1)*n_dens+h_ele)+&
                        matmul(QT(1,1,:,:),u(2*(j-1)*n_vels+u_ele))
                   rhs1((k-1)*n_dens+h_ele)=&
                        rhs1((k-1)*n_dens+h_ele)+&
                        matmul(QT(1,2,:,:),u(2*(j-1)*n_vels+n_vels+u_ele))
                   rhs2((k-1)*n_dens+h_ele)=&
                        rhs2((k-1)*n_dens+h_ele)+&
                        matmul(QT(2,1,:,:),u(2*(j-1)*n_vels+u_ele))
                   rhs2((k-1)*n_dens+h_ele)=&
                        rhs2((k-1)*n_dens+h_ele)+&
                        matmul(QT(2,2,:,:),u(2*(j-1)*n_vels+n_vels+u_ele))

                   QT =dshape_outer_dshape(dnh_t,dnu_t, &
                        detwei*dlocgi(i,:)*dlocgi(j,:)*dlocgi(k,:)*rho(i))


                   rhs1((j-1)*n_dens+h_ele)=&
                        rhs1((j-1)*n_dens+h_ele)+&
                        matmul(QT(1,1,:,:),u(2*(k-1)*n_vels+u_ele))
                   rhs1((j-1)*n_dens+h_ele)=&
                        rhs1((j-1)*n_dens+h_ele)+&
                        matmul(QT(1,2,:,:),u(2*(k-1)*n_vels+n_vels+u_ele))
                   rhs2((j-1)*n_dens+h_ele)=&
                        rhs2((j-1)*n_dens+h_ele)+&
                        matmul(QT(2,1,:,:),u(2*(k-1)*n_vels+u_ele))
                   rhs2((j-1)*n_dens+h_ele)=&
                        rhs2((j-1)*n_dens+h_ele)+&
                        matmul(QT(2,2,:,:),u(2*(k-1)*n_vels+n_vels+u_ele))

                   end do
                end do
             end if
          end do
       end if

          ! sum_i< D_iD_i div w_i,u_i.grad h_i + sum_j div(D_ju_j)>rho_i/2
          ! sum_i< w_i.grad h_i + sum_j div(D_jw_j), D_iD_i div u_i>rho_i/2

       if(.not. GN_FLAG) then
          do i = 1, N_Layers

            QT=dshape_outer_vector_shape( &
                  dnh_t,gradhlocgi(:,i,:),mesh%nu, &
                  detwei*dlocgi(i,:)*dlocgi(i,:)*rho(i)/2.0)

             rhs1((i-1)*n_dens+h_ele)=&
                  rhs1((i-1)*n_dens+h_ele)+&
                  matmul(QT(1,1,:,:),u(2*(i-1)*n_vels+u_ele))
             rhs1((i-1)*n_dens+h_ele)=&
                  rhs1((i-1)*n_dens+h_ele)+&
                  matmul(QT(1,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))
             rhs2((i-1)*n_dens+h_ele)=&
                  rhs2((i-1)*n_dens+h_ele)+&
                  matmul(QT(2,1,:,:),u(2*(i-1)*n_vels+u_ele))
             rhs2((i-1)*n_dens+h_ele)=&
                  rhs2((i-1)*n_dens+h_ele)+&
                  matmul(QT(2,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))

             QT=shape_vector_outer_dshape( &
                  mesh%nh,gradhlocgi(:,i,:),dnu_t, &
                  detwei*dlocgi(i,:)*dlocgi(i,:)*rho(i)/2.0)

             rhs1((i-1)*n_dens+h_ele)=&
                  rhs1((i-1)*n_dens+h_ele)+&
                  matmul(QT(1,1,:,:),u(2*(i-1)*n_vels+u_ele))
             rhs1((i-1)*n_dens+h_ele)=&
                  rhs1((i-1)*n_dens+h_ele)+&
                  matmul(QT(1,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))
             rhs2((i-1)*n_dens+h_ele)=&
                  rhs2((i-1)*n_dens+h_ele)+&
                  matmul(QT(2,1,:,:),u(2*(i-1)*n_vels+u_ele))
             rhs2((i-1)*n_dens+h_ele)=&
                  rhs2((i-1)*n_dens+h_ele)+&
                  matmul(QT(2,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))

             if(i<N_Layers) then

                do j = i+1, N_Layers

                   QT =dshape_outer_vector_shape( &
                        dnh_t,graddlocgi(:,j,:),mesh%nu, &
                        detwei*dlocgi(i,:)*dlocgi(i,:)*rho(i)/2.0)

                   rhs1((i-1)*n_dens+h_ele)=&
                        rhs1((i-1)*n_dens+h_ele)+&
                        matmul(QT(1,1,:,:),u(2*(j-1)*n_vels+u_ele))
                   rhs1((i-1)*n_dens+h_ele)=&
                        rhs1((i-1)*n_dens+h_ele)+&
                        matmul(QT(1,2,:,:),u(2*(j-1)*n_vels+n_vels+u_ele))
                   rhs2((i-1)*n_dens+h_ele)=&
                        rhs2((i-1)*n_dens+h_ele)+&
                        matmul(QT(2,1,:,:),u(2*(j-1)*n_vels+u_ele))
                   rhs2((i-1)*n_dens+h_ele)=&
                        rhs2((i-1)*n_dens+h_ele)+&
                        matmul(QT(2,2,:,:),u(2*(j-1)*n_vels+n_vels+u_ele))


                   QT =shape_vector_outer_dshape( &
                        mesh%nh,graddlocgi(:,j,:),dnu_t, &
                        detwei*dlocgi(i,:)*dlocgi(i,:)*rho(i)/2.0)

                   rhs1((j-1)*n_dens+h_ele)=&
                        rhs1((j-1)*n_dens+h_ele)+&
                        matmul(QT(1,1,:,:),u(2*(i-1)*n_vels+u_ele))
                   rhs1((j-1)*n_dens+h_ele)=&
                        rhs1((j-1)*n_dens+h_ele)+&
                        matmul(QT(1,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))
                   rhs2((j-1)*n_dens+h_ele)=&
                        rhs2((j-1)*n_dens+h_ele)+&
                        matmul(QT(2,1,:,:),u(2*(i-1)*n_vels+u_ele))
                   rhs2((j-1)*n_dens+h_ele)=&
                        rhs2((j-1)*n_dens+h_ele)+&
                        matmul(QT(2,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))
              

                   QT = dshape_outer_dshape(dnh_t,dnu_t, &
                        detwei*dlocgi(i,:)*dlocgi(i,:)*dlocgi(j,:)* &
                        rho(i)/2.0)

                   rhs1((i-1)*n_dens+h_ele)=&
                        rhs1((i-1)*n_dens+h_ele)+&
                        matmul(QT(1,1,:,:),u(2*(j-1)*n_vels+u_ele))
                   rhs1((i-1)*n_dens+h_ele)=&
                        rhs1((i-1)*n_dens+h_ele)+&
                        matmul(QT(1,2,:,:),u(2*(j-1)*n_vels+n_vels+u_ele))
                   rhs2((i-1)*n_dens+h_ele)=&
                        rhs2((i-1)*n_dens+h_ele)+&
                        matmul(QT(2,1,:,:),u(2*(j-1)*n_vels+u_ele))
                   rhs2((i-1)*n_dens+h_ele)=&
                        rhs2((i-1)*n_dens+h_ele)+&
                        matmul(QT(2,2,:,:),u(2*(j-1)*n_vels+n_vels+u_ele))



                   rhs1((j-1)*n_dens+h_ele)=&
                        rhs1((j-1)*n_dens+h_ele)+&
                        matmul(QT(1,1,:,:),u(2*(i-1)*n_vels+u_ele))
                   rhs1((j-1)*n_dens+h_ele)=&
                        rhs1((j-1)*n_dens+h_ele)+&
                        matmul(QT(1,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))
                   rhs2((j-1)*n_dens+h_ele)=&
                        rhs2((j-1)*n_dens+h_ele)+&
                        matmul(QT(2,1,:,:),u(2*(i-1)*n_vels+u_ele))
                   rhs2((j-1)*n_dens+h_ele)=&
                        rhs2((j-1)*n_dens+h_ele)+&
                        matmul(QT(2,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))
               
                end do

             end if
          end do

       end if

       if(.true.) then
          ! sum_i< D_iD_iD_i div w_i, div u_i>rho_i/3
          do i = 1, N_Layers

             QT = dshape_outer_dshape(dnh_t,dnu_t, &
                  detwei*dlocgi(i,:)*dlocgi(i,:)*dlocgi(i,:)* &
                  rho(i)/3.0)

             rhs1((i-1)*n_dens+h_ele)=&
                  rhs1((i-1)*n_dens+h_ele)+&
                  matmul(QT(1,1,:,:),u(2*(i-1)*n_vels+u_ele))
             rhs1((i-1)*n_dens+h_ele)=&
                  rhs1((i-1)*n_dens+h_ele)+&
                  matmul(QT(1,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))
             rhs2((i-1)*n_dens+h_ele)=&
                  rhs2((i-1)*n_dens+h_ele)+&
                  matmul(QT(2,1,:,:),u(2*(i-1)*n_vels+u_ele))
             rhs2((i-1)*n_dens+h_ele)=&
                  rhs2((i-1)*n_dens+h_ele)+&
                  matmul(QT(2,2,:,:),u(2*(i-1)*n_vels+n_vels+u_ele))


             Qu = dshape_dot_dshape(dnh_t,dnu_t, &
                  detwei*rho(i)*alpha*alpha)
 

             rhs1((i-1)*n_dens+h_ele)=&
                  rhs1((i-1)*n_dens+h_ele)+&
                  matmul(Qu,u(2*(i-1)*n_vels+u_ele))
             rhs2((i-1)*n_dens+h_ele)=&
                  rhs2((i-1)*n_dens+h_ele)+&
                  matmul(Qu,u(2*(i-1)*n_vels+n_vels+u_ele))
          end do

       end if

    end do ele_loop

    ewrite(2,*)("end subroutine assemble_MLCM_operator")

  end subroutine assemble_forward_MLCM_operator

  subroutine h_to_m_rhs(mcont,mesh,mrhs,mloc)
    real, dimension(n_dens), intent(in) :: mcont
    type(dg_mesh) :: mesh
    real, dimension(n_moms), intent(out) :: mrhs,mloc
    integer :: ele, d_pos(3)
    integer, dimension(:), pointer :: u_ele, X_ele, h_ele, m_ele
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(mesh%nu%ngi) :: detwei, mcontlocgi

    mrhs=0.0
    
    if (mesh%nh%loc == 3) then 
       d_pos =(/ 1,2,3 /)
    else
       d_pos = (/1,3,6/)
    end if

    ele_loop: do ele = 1, n_elements
       m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)

       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)

       call transform_to_physical(ele_X, mesh%nm, m = mesh%nh, &
            detwei = detwei)

       mcontlocgi=matmul(mcont(h_ele),mesh%nh%n)

       mrhs(m_ele)=mrhs(m_ele)&
            +shape_rhs(mesh%nm,detwei*mcontlocgi)

       mloc(m_ele)=mcont(h_ele(d_pos))

    end do ele_loop

  end subroutine h_to_m_rhs

  subroutine get_moms_layer(Mesh, bcs, D, bottom, rho, u, m)
    type(dg_mesh), intent(in) :: Mesh
    real, dimension(n_layers*n_dens), intent(in), target :: D
    real, dimension(n_dens), intent(in) :: bottom
    real, dimension(n_layers), intent(in) :: rho
    type(bc_info), intent(in) :: bcs
    real, dimension(2*n_vels*n_layers), intent(in), target :: u
    real, dimension(2*n_layers*n_moms), intent(inout), target :: m

    !locals
    real, dimension(n_dens*n_layers):: urhs1,urhs2
    type(csr_matrix) :: mommat2
    type(csr_matrix), dimension (n_layers) :: mommat1
    real, dimension(n_dens) :: mcont
    real, dimension(n_moms) :: mrhs
    integer :: i, j, k

    PetscErrorCode :: ierr
    KSPType :: ksp_type
    PCType :: pc_type

    ksp_type = KSPGMRES
    pc_type = PCICC

    do i=1,n_layers
       call allocate(mommat1(i),mesh%mass_h%sparsity)
    end do
    call allocate(mommat2,mesh%mass_h%sparsity)

    call assemble_forward_MLCM_operator(D,m,u,mesh,bottom,rho, &
         mommat1=mommat1, mommat2=mommat2, &
!         rhs1=urhs((/((i+2*j*n_moms, i=1,n_moms),j=0,n_layers-1)/)),&
!         rhs2=urhs((/((i+(2*j+1)*n_moms, i=1,n_moms),j=0,n_layers-1)/)))
          rhs1=urhs1,&
         rhs2=urhs2)
    print*, "Max urhs1=", maxval(urhs1)
    print*, "Max urhs2=", maxval(urhs2)

    do i=1,n_layers
         call gallopede_solve(mcont,&
              mommat1(i),&
              urhs1((i-1)*n_dens+1:(i-1)*n_dens+n_dens),&
              ksp_type, pc_type,&
              errtol=1.d-32,noi= 30000)
         call h_to_m_rhs(mcont,mesh,mrhs,&
              m(2*(i-1)*n_moms+1:2*(i-1)*n_moms+n_moms))
         call gallopede_solve(m(2*(i-1)*n_moms+1:2*(i-1)*n_moms+n_moms),&
              mommat2,&
              mrhs,&
              ksp_type, pc_type,&
              errtol=1.d-32,noi= 30000)
         call gallopede_solve(mcont,&
              mommat1(i),&
              urhs2((i-1)*n_dens+1:(i-1)*n_dens+n_dens),&
              ksp_type, pc_type,&
              errtol=1.d-32,noi= 30000)
         call h_to_m_rhs(mcont,mesh,mrhs,&
              m(2*(i-1)*n_moms+n_moms+1:2*(i-1)*n_moms+n_moms+n_moms))
         call gallopede_solve(m(2*(i-1)*n_moms+n_moms+1:2*(i-1)*n_moms+n_moms+n_moms),&
              mommat2,&
              mrhs,&
              ksp_type, pc_type,&
              errtol=1.d-32,noi= 30000)
   end do

    do i=1,n_layers
       call deallocate( mommat1(i) )
    end do
    call deallocate( mommat2 )

  end subroutine get_moms_layer


  subroutine tensor_addto(mommat,i,j,ele1,ele2,QT,switch)
    type(block_csr_matrix), intent(inout) :: mommat
    integer, intent(in) :: i,j
    integer, dimension(:), intent(in) :: ele1,ele2
    real, dimension(:,:,:,:), intent(in) :: QT
    logical, optional :: switch

    if (present(switch)) then

       call addto(mommat,2*(i-1)+1,2*(j-1)+1, &
            ele1,ele2,transpose(QT(1,1,:,:)))
       call addto(mommat,2*(i-1)+1,2*(j-1)+2, &
            ele1,ele2,transpose(QT(2,1,:,:)))
       call addto(mommat,2*(i-1)+2,2*(j-1)+1, &
            ele1,ele2,transpose(QT(1,2,:,:)))
       call addto(mommat,2*(i-1)+2,2*(j-1)+2, &
            ele1,ele2,transpose(QT(2,2,:,:)))

    else

       call addto(mommat,2*(i-1)+1,2*(j-1)+1, &
            ele1,ele2,QT(1,1,:,:))
       call addto(mommat,2*(i-1)+1,2*(j-1)+2, &
            ele1,ele2,QT(1,2,:,:))
       call addto(mommat,2*(i-1)+2,2*(j-1)+1, &
            ele1,ele2,QT(2,1,:,:))
       call addto(mommat,2*(i-1)+2,2*(j-1)+2, &
            ele1,ele2,QT(2,2,:,:))

    end if

  end subroutine tensor_addto

  subroutine u_smooth(u,m,d,mesh)

    real, dimension(:), intent(inout) :: u
    real, dimension(:), intent(in) :: m,d
    type(dg_mesh), intent(in) :: mesh
    real, dimension(2*n_vels) :: rhs
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nu%ngi) :: detwei, dlocgi
    real :: d_m

    type(block_csr_matrix) :: mommat

    ! Local element information
    integer, dimension(:), pointer :: u_ele, X_ele, h_ele, m_ele
    integer :: i, ele


    real, dimension(mesh%nu%loc,mesh%nu%loc) :: Qu
    real, dimension(mesh%nu%loc,mesh%nm%loc) :: Qum
    real, dimension(mesh%nu%loc,mesh%nh%loc) :: Qud
    real, dimension(2,2,mesh%nu%loc,mesh%nu%loc) :: QT

         !locals
    PetscErrorCode :: ierr
    KSPType :: ksp_type
    PCType :: pc_type

    ksp_type = KSPCG
    pc_type = PCICC

     call allocate(mommat,mesh%Mass_u%sparsity, (/2,2/))

       call zero(mommat)

    ewrite(2,*)("ele loop")

    do i=0, n_layers-1


       d_m=1   !0.*sum(d(i*n_dens+1:(i+1)*n_dens))/n_dens
       rhs=0.

       ele_loop: do ele = 1, n_elements
          
          u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
          m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
          X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
          h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
          
          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)

          call transform_to_physical(ele_X, mesh%nm, m = mesh%nh, &
               dn_t = dnm_t, dm_t = dnh_t, &
               detwei = detwei)
          call transform_to_physical(ele_X, mesh%nm, m = mesh%nu, &
               dm_t = dnu_t)

          Qu=shape_shape(mesh%nu,mesh%nu,detwei)
          
          !rhs(u_ele)=rhs(u_ele)+matmul(Qu,u(2*i*n_vels+u_ele))
          !rhs(n_vels+u_ele)=rhs(n_vels+u_ele)&
          !          +matmul(Qu,u(2*i*n_vels+n_vels+u_ele))


          dlocgi=0!*matmul(d(i*n_dens+h_ele)-d_m,mesh%nh%n)

           Qum=shape_shape(mesh%nu,mesh%nm,detwei*(dlocgi+d_m))
          rhs(u_ele)=rhs(u_ele)+matmul(Qum,m(2*i*n_moms+m_ele))
          rhs(n_vels+u_ele)=rhs(n_vels+u_ele)&
               +matmul(Qum,m(2*i*n_moms+n_moms+m_ele))



          Qud=shape_shape(mesh%nu,mesh%nh,detwei)
!          rhs(u_ele)=rhs(u_ele)+matmul(Qud,d(n_dens+h_ele))
!          rhs(n_vels+u_ele)=rhs(n_vels+u_ele)&
!               +matmul(Qud,d(h_ele))

          if (i==0) then

             call addto(mommat,1,1,u_ele,u_ele,Qu)
             call addto(mommat,2,2,u_ele,u_ele,Qu)
            
          
!          Qu=dshape_dot_dshape(dnu_t,dnu_t,detwei*(1.0/30.0))
             QT=dshape_outer_dshape(dnu_t,dnu_t,d_m**3*detwei/3.0)

          call addto(mommat,1,1, u_ele,u_ele,QT(1,1,:,:))
          call addto(mommat,2,1, u_ele,u_ele,QT(2,1,:,:))
          call addto(mommat,1,2, u_ele,u_ele,QT(1,2,:,:))
          call addto(mommat,2,2, u_ele,u_ele,QT(2,2,:,:))


          QT=dshape_outer_dshape(dnu_t,dnu_t,d_m**2*dlocgi*detwei)

          call addto(mommat,1,1, u_ele,u_ele,QT(1,1,:,:))
          call addto(mommat,2,1, u_ele,u_ele,QT(2,1,:,:))
          call addto(mommat,1,2, u_ele,u_ele,QT(1,2,:,:))
          call addto(mommat,2,2, u_ele,u_ele,QT(2,2,:,:))

          QT=dshape_outer_dshape(dnu_t,dnu_t,d_m*dlocgi**2*detwei)

          call addto(mommat,1,1, u_ele,u_ele,QT(1,1,:,:))
          call addto(mommat,2,1, u_ele,u_ele,QT(2,1,:,:))
          call addto(mommat,1,2, u_ele,u_ele,QT(1,2,:,:))
          call addto(mommat,2,2, u_ele,u_ele,QT(2,2,:,:))

          QT=dshape_outer_dshape(dnu_t,dnu_t,dlocgi**3*detwei/3.0)

          call addto(mommat,1,1, u_ele,u_ele,QT(1,1,:,:))
          call addto(mommat,2,1, u_ele,u_ele,QT(2,1,:,:))
          call addto(mommat,1,2, u_ele,u_ele,QT(1,2,:,:))
          call addto(mommat,2,2, u_ele,u_ele,QT(2,2,:,:))
         

       end if

    end do ele_loop

       call gallopede_block_solve(u(2*i*n_vels+1:2*(i+1)*n_vels),mommat,&
            rhs,ksp_type, pc_type,&
               errtol=1.d-32,noi= 30000)

    end do

    call deallocate(mommat)

  end subroutine u_smooth



 subroutine mini_out(field,mesh_EV,f_name,n_locs,levels)
   implicit none
   real, dimension(:), intent(in) :: field
     integer, intent(in) :: levels
     integer, intent(in) :: n_locs
     integer, dimension(:), pointer :: mesh_EV
     character(len=*), intent(in) :: f_name
     


     integer ele, i
     integer, dimension(:), pointer :: h_ele

     open(233,file=f_name)
     
     do i=0,levels-1
        do ele=1,n_elements
           h_ele=>mesh_EV((ELE-1)*n_locs+1:ELE*n_locs)

           write(233,*) field(i*size(field)/levels+h_ele)

        end do
     end do

     close(233)

   end subroutine mini_out

   subroutine swap_rhs(RHS)
     real, dimension(:), intent(inout) :: RHS

     real, dimension(size(RHS)) :: work
     integer :: i, N_U_LAYERS, n_vls


     N_U_LAYERS=n_layers
     n_vls=size(rhs)/(2*n_layers)
     do i=0,n_layers-1
        print*, i, n_vls, size(work)
        work(2*i*n_vls+1:2*i*n_vls+n_vls)=&
             RHS(i*n_vls+1:i*n_vls+n_vls)
        work((2*i+1)*n_vls+1:(2*i+1)*n_vls+n_vls)=&
             RHS((N_U_LAYERS+i)*n_vls+1:(N_U_LAYERS+i)*n_vls+n_vls)
     end do

     RHS=work

   end subroutine swap_rhs

   subroutine constant_u_set(u,D,mesh)

 type(dg_mesh), intent(in) :: Mesh
    real, dimension(:), intent(out) :: u

    real, dimension(:), intent(in) :: D
     

    !locals
    integer :: ele,i,j
    real, dimension(n_layers) :: dd

      type(csr_matrix) :: den_mat

    ! Local element information
    integer, dimension(:), pointer :: u_ele, X_ele,m_ele,h_ele
    real, dimension(2,mesh%nu%loc) :: ele_X
    real, dimension(2,mesh%nx%loc) :: ele_X2

    real, dimension(n_vels) :: rhs
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nu%ngi) :: detwei, dlocgi

    real, dimension(mesh%nu%loc,mesh%nu%loc) :: Qu
    real, dimension(mesh%nu%loc,mesh%nm%loc) :: Qum
    real, dimension(mesh%nu%loc,mesh%nh%loc) :: Qud
    real, dimension(2,2,mesh%nu%loc,mesh%nu%loc) :: QT

         !locals
    PetscErrorCode :: ierr
    KSPType :: ksp_type
    PCType :: pc_type

    ksp_type = KSPCG
    pc_type = PCICC

    dd=(/ 75.0, 275.0 /)
    call allocate(den_mat,mesh%mass_u%sparsity)



    do i=1,n_layers

    call zero(den_mat) 
    rhs=0.0

    do ele=1,n_elements
       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)

       ele_X2(1,:)=mesh%X(X_ele)
       ele_X2(2,:)=mesh%Y(X_ele)

       call transform_to_physical(ele_X2, mesh%nm, m = mesh%nh, &
            dn_t = dnm_t, dm_t = dnh_t, &
            detwei = detwei)
       call transform_to_physical(ele_X2, mesh%nm, m = mesh%nu, &
            dm_t = dnu_t)


       if(mesh%nu%loc == 3) then

          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)

       else

       ele_X(1,(/ 1,3,6 /))=mesh%X(X_ele)
       ele_X(2,(/ 1,3,6 /))=mesh%Y(X_ele)
       
       ele_X(:,2)=0.5*(ele_X(:,1)+ele_X(:,3))
       ele_X(:,4)=0.5*(ele_X(:,1)+ele_X(:,6))
       ele_X(:,5)=0.5*(ele_X(:,3)+ele_X(:,6))

    end if
    


       u(2*(i-1)*n_vels+u_ele)=1.0
 !      u(2*(i-1)*n_vels+u_ele)=0.1*sin(2*3.14159*ele_X(1,:)/1000.0)
 !      u(2*(i-1)*n_vels+u_ele)=1.5*(1.0-dd(i)/D((i-1)*n_dens+u_ele))
       u((2*(i-1)+1)*n_vels+u_ele)=0.0


       dlocgi(:) = matmul(D((i-1)*n_dens+h_ele),mesh%nh%n)

       Qu=shape_shape(mesh%nu,mesh%nu,detwei*dlocgi)

       call addto(den_mat,u_ele,u_ele,Qu)

       dlocgi(:) = matmul(D((i-1)*n_dens+h_ele)-dd(i),mesh%nh%n)

       rhs(u_ele)=rhs(u_ele)+shape_rhs(mesh%nu,&
            1.5*(dlocgi*detwei))

    end do

        call gallopede_solve(u(2*(i-1)*n_vels+1:(2*i-1)*n_vels),den_mat,&
            rhs,ksp_type, pc_type,&
               errtol=1.d-24,noi= 30000)
     end do

   call deallocate(den_mat)

end subroutine constant_u_set

subroutine assemble_CG_MLCM_operator(D,m,u,mesh, &
       bottom, rho, &
       mommat, &
       rhs1,rhs2)
    implicit none
    type(block_csr_matrix), intent(inout) :: mommat
    real, dimension(:), intent(out) :: rhs1, rhs2
    real, dimension(:), intent(in) :: m
    real, dimension(:), intent(in) :: u
    real, dimension(:), intent(in) :: D
    type(dg_mesh) :: mesh
    real, intent(in), dimension(:) :: rho
    real, dimension(:), intent(in) :: bottom
    
    !locals
    integer :: ele,i,j,k,iloc,jloc,iX,im,jm
    real, dimension(N_layers,mesh%nh%ngi) :: dlocgi,d0locgi
    real, dimension(2,N_layers,mesh%nu%ngi) :: graddlocgi, gradhlocgi
    real, dimension(2,mesh%nu%ngi) :: gradblocgi,gradh0locgi
    real, dimension(2,mesh%nm%loc) :: ele_X
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t
    real, dimension(mesh%nm%loc,mesh%nm%ngi,2) :: dnm_t
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
    real, dimension(mesh%nu%ngi) :: detwei

    ! Local element information
    integer, dimension(:), pointer :: u_ele, X_ele, h_ele, m_ele

    real, dimension(mesh%nu%loc) :: lumped_mass
    real, dimension(mesh%nu%loc,mesh%nu%loc) :: Qu
    real, dimension(mesh%nu%loc,mesh%nh%loc) :: Qum
    real, dimension(2,2,mesh%nu%loc,mesh%nu%loc) :: QT, QT_T
    real, dimension(2,mesh%nu%loc,mesh%nu%loc) :: QuT
    real, dimension(2,2,mesh%nu%loc,mesh%nm%loc) :: QTum

    logical :: log1
    real, dimension(mesh%nu%ngi) :: xlocgi, ylocgi


    ewrite(1,*)("subroutine assemble_MLCM_operator")


!    assert(maxval(mommat%val)<huge(0.0))
    assert(size(D)==n_layers*n_dens)
    assert(size(m)==2*n_layers*n_moms)
    assert(size(bottom)==n_dens)

    call zero(mommat)

    log1=mommat%blocks(1)==2*N_Layers
    assert(log1)
    log1=mommat%blocks(2)==2*N_Layers
    assert(log1)

    ewrite(2,*)("ele loop")

    rhs1=0.
    rhs2=0.

    ele_loop: do ele = 1, n_elements

       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       m_ele=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
!       allocate(h_ele(3))
!       h_ele=mesh%EVList_h((ELE-1)*mesh%Nh%LOC+(/1,3,6/))
       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)

       call transform_to_physical(ele_X, mesh%nm, m = mesh%nx, &
            dn_t = dnm_t, dm_t = dnh_t, &
            detwei = detwei)
       call transform_to_physical(ele_X, mesh%nm, m = mesh%nu, &
            dm_t = dnu_t)

       do i=1,mesh%nu%loc
          lumped_mass(i)=sum(row_val_ptr(mesh%mass_u,u_ele(i)))
       end do

       !get dlocgi
       dlocgi = 0.
       do i = 1, N_Layers
          dlocgi(i,:) = matmul(D((i-1)*n_dens+h_ele),mesh%nh%n)
           D0locgi(i,:)=sum(D((I-1)*n_dens+1:I*n_dens))/n_dens
       end do

       assert(maxval(dlocgi)<huge(0.0))

       !get gradblocgi
       gradblocgi = 0.
       do iX = 1,2
          gradblocgi(iX,:) = matmul(bottom(h_ele),dnh_t(:,:,iX))
       end do

       assert(maxval(gradblocgi)<huge(0.0))

       !get graddlocgi
       graddlocgi = 0.
       forall(i = 1:N_Layers, iX=1:2)
          graddlocgi(iX,i,:) = matmul(D((i-1)*n_dens+h_ele),dnh_t(:,:,iX))
       end forall

       assert(maxval(graddlocgi)<huge(0.0))

       !get gradhlocgi
       gradhlocgi = 0.
       do i = 1, N_Layers
          gradhlocgi(:,i,:) = gradblocgi(:,:) - &
               sum(graddlocgi(:,i+1:N_layers,:),2)
       end do

       gradh0locgi(:,:) = -gradblocgi(:,:) + &
               sum(graddlocgi(:,i:N_layers,:),2)

       assert(maxval(gradhlocgi)<huge(0.0))

       !RIGHT HAND SIDE

       !get local mass matrix for RHS
 

 

       do i = 1,N_Layers
!      Qum=shape_shape(mesh%nu,mesh%nm,detwei*dlocgi(i,:))
      Qum=shape_shape(mesh%nu,mesh%nh,detwei)

      assert(maxval(Qum)<huge(0.0))

          rhs1((i-1)*n_vels+u_ele) = rhs1((i-1)*n_vels+u_ele) + &
               matmul(Qum,rho(i)*m((i-1)*2*n_dens+h_ele))
          rhs2((i-1)*n_vels+u_ele) = rhs2((i-1)*n_vels+u_ele) + &
               matmul(Qum,rho(i)*m((i-1)*2*n_dens+n_dens+h_ele))

          
          rhs1((i-1)*n_vels+u_ele) = rhs1((i-1)*n_vels+u_ele) - &
               shape_rhs(mesh%nu,rho(i)*g0*dlocgi(i,:)*gradh0locgi(1,:)*detwei)

          rhs2((i-1)*n_vels+u_ele) = rhs2((i-1)*n_vels+u_ele) - &
               shape_rhs(mesh%nu,rho(i)*g0*dlocgi(i,:)*gradh0locgi(2,:)*detwei)
          

       end do

       if (RIGID_LID) then
          rhs1(n_layers*n_vels+u_ele)=0.0
          rhs2(n_layers*n_vels+u_ele)=0.0
       end if


       !LEFT HAND SIDE

       !horizontal kinetic energy

       do i = 1, N_Layers
          Qu=shape_shape(mesh%nu,mesh%nu,&
               detwei *dlocgi(i,:)*rho(i))
!               detwei*rho(i))
          assert(maxval(Qu)<huge(0.0))

          call addto(mommat,2*(i-1)+1,2*(i-1)+1, &
               u_ele,u_ele,Qu)
          call addto(mommat,2*(i-1)+2,2*(i-1)+2, &
               u_ele,u_ele,Qu)
       end do

       !vertical kinetic energy
       
       ! sum_i \rho_i <w_i.grad h_i,D_i u_i.grad h_i>

       if(.not. GN_FLAG) then
          !gradient-free terms
          do i = 1, N_Layers
             QT=shape_shape_vector_outer_vector( &
                  mesh%nu,mesh%nu,detwei*dlocgi(i,:)*rho(i), &
                  gradhlocgi(:,i,:),gradhlocgi(:,i,:))
             assert(maxval(QT)<huge(0.0))
!             call tensor_addto(mommat,i,i,u_ele,u_ele,QT)

             forall (im=1:2,jm=1:2)
                QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
             end forall

             call tensor_addto(mommat,i,i,u_ele,u_ele,0.5*(QT_T+QT))
          end do
       end if

       if(.not. GN_FLAG) then

          ! sum_i rho_i<w_i.grad h_i, D_i sum div (D_ju_j) >
          ! sum_i rho_i<u_i.grad h_i, D_i sum div (D_jw_j) >

          do i = 1, N_Layers
             if(i<N_Layers) then
                do j = i+1, N_Layers

                   QT=shape_shape_vector_outer_vector( &
                        mesh%nu,mesh%nu,detwei*dlocgi(i,:)*rho(i), &
                        gradhlocgi(:,i,:),graddlocgi(:,j,:))
!                   call tensor_addto(mommat,i,j,u_ele,u_ele,QT)
!                   call tensor_addto(mommat,j,i,u_ele,u_ele,QT,.true.)

                   QT=QT+shape_vector_outer_dshape( &
                        mesh%nu,gradhlocgi(:,i,:),dnu_t,&
                        detwei*dlocgi(i,:)*dlocgi(j,:)*rho(i))

                     forall (im=1:2,jm=1:2)
                        QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
                     end forall

                   call tensor_addto(mommat,i,j,u_ele,u_ele,QT)
                   call tensor_addto(mommat,j,i,u_ele,u_ele,QT_T)
                end do
             end if
          end do

       end if

       if(.not. GN_FLAG) then
          ! sum_i<D_i sum_j div(D_jw_j), sum_k div(D_ku_k) >rho_i
          do i = 1, N_Layers
             if(i<N_Layers) then
                do j = i+1, N_Layers
                   do k = j, N_Layers

                      QT=shape_shape_vector_outer_vector( &
                           mesh%nu,mesh%nu,detwei*dlocgi(i,:)*rho(i), &
                           graddlocgi(:,j,:),graddlocgi(:,k,:))

                      forall (im=1:2,jm=1:2)
                         QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
                      end forall

                      call tensor_addto(mommat,j,k,u_ele,u_ele,QT)
                      if (k>j) &
                           call tensor_addto(mommat,k,j,u_ele,u_ele,QT_T)


                      QT=shape_vector_outer_dshape( &
                           mesh%nu,graddlocgi(:,j,:),dnu_t,&
                           detwei*dlocgi(i,:)*dlocgi(k,:)*rho(i))

!                      call tensor_addto(mommat,j,k,u_ele,u_ele,QT)
!                      call tensor_addto(mommat,k,j,u_ele,u_ele,QT,.true.)
!                      call tensor_addto(mommat,j,k,u_ele,u_ele,(QT_T+QT))
!                      call tensor_addto(mommat,k,j,u_ele,u_ele,QT_T)
                      

                      QT=dshape_outer_vector_shape( &
                           dnu_t,graddlocgi(:,k,:),mesh%nu,&
                           detwei*dlocgi(i,:)*dlocgi(j,:)*rho(i))

                      forall (im=1:2,jm=1:2)
                         QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
                      end forall

                      call tensor_addto(mommat,j,k,u_ele,u_ele,QT+QT_T)
                      if (k>j) &
                           call tensor_addto(mommat,k,j,u_ele,u_ele,QT+QT_T)

                      QT =dshape_outer_dshape(dnu_t,dnu_t, &
                           detwei*dlocgi(i,:)*dlocgi(j,:)*dlocgi(k,:)*rho(i))

                      forall (im=1:2,jm=1:2)
                         QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
                      end forall


                      call tensor_addto(mommat,j,k,u_ele,u_ele,QT)
                      if (k>j) &
                           call tensor_addto(mommat,k,j,u_ele,u_ele,QT_T)

                   end do
                end do
             end if
          end do
       end if

          ! sum_i< D_iD_i div w_i,u_i.grad h_i + sum_j div(D_ju_j)>rho_i/2
          ! sum_i< w_i.grad h_i + sum_j div(D_jw_j), D_iD_i div u_i>rho_i/2

       if(.not. GN_FLAG) then
          do i = 1, N_Layers

             QT=dshape_outer_vector_shape( &
                  dnu_t,gradhlocgi(:,i,:),mesh%nu, &
                  detwei*dlocgi(i,:)*dlocgi(i,:)*rho(i)/2.0)

             forall (im=1:2,jm=1:2)
                QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
             end forall

             call tensor_addto(mommat,i,i,u_ele,u_ele,QT_T+QT)

   

             if(i<N_Layers) then

                do j = i+1, N_Layers

                   QT =dshape_outer_vector_shape( &
                        dnu_t,graddlocgi(:,j,:),mesh%nu, &
                        detwei*dlocgi(i,:)*dlocgi(i,:)*rho(i)/2.0)
!                  call tensor_addto(mommat,i,j,u_ele,u_ele,QT)
!                  call tensor_addto(mommat,j,i,u_ele,u_ele,QT,.true.)

                   QT = QT+dshape_outer_dshape(dnu_t,dnu_t, &
                        detwei*dlocgi(i,:)*dlocgi(i,:)*dlocgi(j,:)* &
                        rho(i)/2.0)

                   forall (im=1:2,jm=1:2)
                      QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
                   end forall

                   call tensor_addto(mommat,i,j,u_ele,u_ele,QT)
                   call tensor_addto(mommat,j,i,u_ele,u_ele,QT_T)

                end do

             end if
          end do

       end if

       if(.true.) then
          ! sum_i< D_iD_iD_i div w_i, div u_i>rho_i/3
          do i = 1, N_Layers

             QT = dshape_outer_dshape(dnu_t,dnu_t, &
!                  detwei*1.0*dlocgi(i,:))
                  detwei*dlocgi(i,:)*dlocgi(i,:)*dlocgi(i,:)* &
!                  detwei*dlocgi(i,:)*d0locgi(i,:)*d0locgi(i,:)* &
!              detwei*dlocgi(i,:)*dlocgi(i,:)*&
                  rho(i)/3.0)

             forall (im=1:2,jm=1:2)
                QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
             end forall

             call tensor_addto(mommat,i,i,u_ele,u_ele,0.5*(QT+QT_T))
!             call tensor_addto(mommat,i,i,u_ele,u_ele,0.5*QT)
!             call tensor_addto(mommat,i,i,u_ele,u_ele,0.5*QT,.true.)


             QU = dshape_dot_dshape(dnu_t,dnu_t, &
                  detwei*rho(i)*alpha*alpha)
             

             call addto(mommat,2*(i-1)+1,2*(i-1)+1, &
                  u_ele,u_ele,Qu)
             call addto(mommat,2*(i-1)+2,2*(i-1)+2, &
                  u_ele,u_ele,Qu)

          end do

          if (RIGID_LID) then

             !<div w_R sum div Du>
             !<div Dw_i div u_r>
             
             do i=1,n_layers


                if (.false.) then
             QT = dshape_outer_dshape(dnu_t,dnu_t, &
                  detwei*dlocgi(i,:))
             forall (im=1:2,jm=1:2)
                QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
             end forall
             call tensor_addto(mommat,i,n_layers+1,u_ele,u_ele,QT_T)
             call tensor_addto(mommat,n_layers+1,i,u_ele,u_ele,QT)
             QT= dshape_outer_vector_shape(dnu_t,&
                  graddlocgi(:,i,:),mesh%nu,detwei)
             forall (im=1:2,jm=1:2)
                QT_T(im,jm,:,:)=transpose(QT(jm,im,:,:))
             end forall

             call tensor_addto(mommat,i,n_layers+1,u_ele,u_ele,QT_T)
             call tensor_addto(mommat,n_layers+1,i,u_ele,u_ele,QT)

             end if

             if (.true.) then

 !            QuT=-shape_dshape(mesh%nu,dnu_t,rho(i)*dlocgi(i,:)*detwei)             
            QuT=dshape_shape(dnu_t,mesh%nu,rho(i)*dlocgi(i,:)*detwei)             
             call addto(mommat,2*(n_layers)+1,2*(i-1)+1, &
               u_ele,u_ele,QuT(1,:,:))
             call addto(mommat,2*(n_layers)+1,2*(i-1)+2, &
               u_ele,u_ele,QuT(2,:,:))
             call addto(mommat,2*(i-1)+1,2*(n_layers)+1, &
               u_ele,u_ele,transpose(QuT(1,:,:)))
             call addto(mommat,2*(i-1)+2,2*(n_layers)+1, &
               u_ele,u_ele,transpose(QuT(2,:,:)))

             Qu=0.0
             Qu(1,1)=1.0e64

             call addto(mommat,2*(n_layers)+1,2*(n_layers)+1, &
               u_ele,u_ele,Qu(:,:))

!             QuT=-shape_shape_vector(mesh%nu,mesh%nu,&
!                  rho(i)*detwei,graddlocgi(:,i,:))
!             
!              call addto(mommat,2*(n_layers)+1,2*(i-1)+1, &
!               u_ele,u_ele,QuT(1,:,:))
!             call addto(mommat,2*(n_layers)+1,2*(i-1)+2, &
!               u_ele,u_ele,QuT(2,:,:))
!             call addto(mommat,2*(i-1)+1,2*(n_layers)+1, &
!               u_ele,u_ele,transpose(QuT(1,:,:)))
!             call addto(mommat,2*(i-1)+2,2*(n_layers)+1, &
!               u_ele,u_ele,transpose(QuT(2,:,:)))

             Qu=shape_shape(mesh%nu,mesh%nu,detwei)
             
              call addto(mommat,2*(n_layers)+2,2*(n_layers)+2, &
               u_ele,u_ele,Qu)
             
              QU = -dshape_dot_dshape(dnu_t,dnu_t, &
                  detwei*100.0)

              call addto(mommat,2*(n_layers)+1,2*(n_layers)+1, &
                   u_ele,u_ele,Qu)

           end if

           if (.false.) then

             QuT=dshape_shape(dnu_t,mesh%nu,rho(i)*dlocgi(i,:)*detwei)             
             call addto(mommat,2*(n_layers)+1,2*(i-1)+1, &
               u_ele,u_ele,QuT(1,:,:))
             call addto(mommat,2*(n_layers)+1,2*(i-1)+2, &
               u_ele,u_ele,QuT(2,:,:))
                
             Qu=-shape_shape(mesh%nu,mesh%nu,detwei)
             
             call addto(mommat,2*(i-1)+1,2*(n_layers)+1, &
               u_ele,u_ele,Qu)
             call addto(mommat,2*(i-1)+2,2*(n_layers)+2, &
               u_ele,u_ele,Qu)
             
             QuT=dshape_shape(dnu_t,mesh%nu,rho(i)*detwei)  

             call addto(mommat,2*(n_layers)+2,2*(i-1)+1, &
               u_ele,u_ele,-QuT(2,:,:))
             call addto(mommat,2*(n_layers)+2,2*(i-1)+2, &
               u_ele,u_ele,QuT(1,:,:))

           end if


             end do
          
          end if

       end if

!       deallocate(h_ele)
    end do ele_loop

    ewrite(2,*)("end subroutine assemble_MLCM_operator")


  end subroutine assemble_CG_MLCM_operator

 subroutine pressure_operator_solve(Mesh, D, bottom, rho, u_out,&
      m, bcs, u0)
    type(dg_mesh), intent(in) :: Mesh
    real, dimension(:), intent(in) :: D
    real, dimension(:), intent(in) :: bottom
    real, dimension(:), intent(in) :: rho
    type(bc_info), intent(in) :: bcs
    real, dimension(:), intent(inout) :: u_out
    real, dimension(:), intent(in), optional :: u0
    real, dimension(:), intent(in) :: m

    !locals
    PetscErrorCode :: ierr
    KSPType :: ksp_type
    PCType :: pc_type

    real, dimension(2*n_vels*n_layers):: urhs,u
    type(block_csr_matrix) :: mommat
    integer :: i, j, u_it,k

    ewrite(1,*)("subroutine pressure operator solve")



    get_vel_count = get_vel_count + 1

    urhs = 0.

     call allocate(mommat,mesh%mass_u%sparsity,&
            (/2 * N_layers, 2 * N_layers /))   
    assert(maxval(D)<huge(0.0))


    if(present(u0)) then
       u = u0
    else
       u = 0.
    end if

    call assemble_CG_MLCM_operator(D,m,u,mesh,bottom,rho, &
         mommat=mommat, &
         rhs1=urhs(1:n_vels*n_layers),&
         rhs2=urhs(n_vels*n_layers+1:2*n_vels*n_layers))

    ewrite(2,*) huge(0.0)
!    assert(maxval(mommat%val)<huge(0.0))

    ewrite(2,*)("solving GN operator equation")
    !solve GN equation
!    ksp_type = KSPCG
!    pc_type = PCICC

if (RIGID_LID) then
    ksp_type=KSPCG
    pc_type=PCICC
else
   ksp_type = KSPCG
    pc_type = PCICC
end if


!    ewrite(1,*) 'maxval(mommat)', maxval(mommat%val)
    

!    ewrite(1,*) 'maxval(mommat)=', maxval(mommat%val)
!    assert(maxval(mommat%val)<huge(0.0))

    print*, bcs%N_tangents

!    call    mini_out(urhs(1:n_layers*n_vels),mesh%EVList_u,&
!         'urhs1.dat',mesh%nu%loc,n_layers)
!    call    mini_out(urhs(n_layers*n_vels+1:2*n_layers*n_vels),mesh%EVList_u,&
!         'urhs2.dat',mesh%nu%loc,n_layers)

    if(bcs%n_tangents .ne. 0) then
       print*, 'boundary solve'
       call gallopede_block_solve_lift(u,&
            urhs(1:size(urhs)/2),&
            urhs(size(urhs)/2+1:size(urhs)), &
            mommat,&
            bcs,ksp_type, pc_type,&
            rerrtol=1.d-20,abserrtol=1.d-12,noi= 10000,zero=.false.)
    else
       print*, "No boundaries"
       call swap_rhs(urhs)
       call gallopede_block_solve(u,mommat,&
               urhs,ksp_type, pc_type,&
               errtol=1.d-32,noi= 30000,zero=.false.)
!           call u_smooth(u,m,d,mesh)
       end if

!       u_out=0.0
!       urhs=0.0

       do i=1,n_layers
!          do j=1,n_layers
!             do k=1,n_vels
!             urhs(2*(i-1)*n_vels+k)=&
!                  urhs(2*(i-1)*n_vels+k)&
!                  +sum(row_val_ptr(mommat,2*(i-1)+1,2*(j-1)+1,k)&
!                  *u(2*(j-1)*n_vels+row_m_ptr(mommat,k))&
!                  +row_val_ptr(mommat,2*(i-1)+1,2*(j-1)+2,k)&
!                  *u((2*(j-1)+1)*n_vels+row_m_ptr(mommat,k)))

!             urhs((2*(i-1)+1)*n_vels+k)=&
!                  urhs((2*(i-1)+1)*n_vels+k)&
!                  +sum(row_val_ptr(mommat,2*(i-1)+2,2*(j-1)+1,k)&
!                  *u(2*(j-1)*n_vels+row_m_ptr(mommat,k))&
!                  +row_val_ptr(mommat,2*(i-1)+2,2*(j-1)+2,k)&
!                  *u((2*(j-1)+1)*n_vels+row_m_ptr(mommat,k)))
!             end do
!          end do

!          call gallopede_solve(u_out(1:n_vels),mesh%mass_h,&
!               urhs(1:n_vels),ksp_type, pc_type,&
!               errtol=1.d-32,noi= 30000)          
!         call gallopede_solve(u_out(n_vels+1:2*n_vels),mesh%mass_h,&
!               urhs(n_vels+1:2*n_vels),ksp_type, pc_type,&
!               errtol=1.d-32,noi= 30000) 


             u_out(2*(i-1)*n_vels+1:(2*(i-1)+1)*n_vels)=&
                  d((i-1)*n_dens+1:i*n_dens)*u(2*(i-1)&
                  *n_vels+1:(2*(i-1)+1)*n_vels)
!                  -m(2*(i-1)*n_vels+1:(2*(i-1)+1)*n_vels)/rho(i)
!                  -u_out(2*(i-1)*n_vels+1:(2*(i-1)+1)*n_vels)/rho(i)&
!            +d((i-1)*n_dens+1:i*n_dens)*u(2*(i-1)*n_vels+1:(2*(i-1)+1)*n_vels)
!             +m(2*(i-1)*n_vels+1:(2*(i-1)+1)*n_vels)

               u_out((2*(i-1)+1)*n_vels+1:2*i*n_vels)=&
                    d((i-1)*n_dens+1:i*n_dens)&
                    *u(2*(i-1)*n_vels+1:(2*(i-1)+1)*n_vels)
!                    -m((2*(i-1)+1)*n_vels+1:2*i*n_vels)/rho(i)
!                  -U_out((2*(i-1)+1)*n_vels+1:2*i*n_vels)/rho(i)&
!             +d((i-1)*n_dens+1:i*n_dens)*u((2*(i-1)+1)*n_vels+1:2*i*n_vels)
!             +m((2*(i-1)+1)*n_vels+1:2*n_vels)
          end do

          

       call deallocate(mommat)

    ewrite(2,*)("END subroutine get_vels_layers")

  end subroutine pressure_operator_solve

  subroutine get_u1_rigid_lid(D,u,mesh)

   real, dimension(2*n_vels*n_layers), intent(inout):: u
   real, dimension(n_dens*n_layers), intent(in) :: D
   type(dg_mesh) :: mesh

   integer :: j

   real, dimension(n_dens) :: rel_D


   call get_density_relation(D,mesh,rel_D)

   u(1:2*n_vels)=0.0

   do j=2,n_layers

      u(1:n_vels)=u(1:n_vels)-rel_d*d((j-1)*n_dens+1:j*n_dens)&
           *u(2*(j-1)*n_vels+1:(2*(j-1)+1)*n_vels)
      u(1+n_vels:2*n_vels)=u(1+n_vels:2*n_vels)&
           -rel_d*d((j-1)*n_dens+1:j*n_dens)&
           *u((2*(j-1)+1)*n_vels+1:2*(j)*n_vels)

   end do

    end subroutine get_u1_rigid_lid

    subroutine point_set(variable,point_list,mesh,ele)
      real, dimension(:), intent(in) :: variable
      integer, dimension(:), pointer :: point_list
      type(dg_mesh) :: mesh
      integer, intent(in) :: ele

    if (size(variable)==2*n_vels*n_layers) then
      point_list=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
    elseif (size(variable)==2*n_moms*n_layers) then
       point_list=>mesh%EVList_m((ELE-1)*mesh%Nm%LOC+1:ELE*mesh%Nm%LOC)
    elseif (size(variable)==n_dens*n_layers) then
       point_list=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)
    elseif (size(variable)==n_verts) then
       point_list=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
    end if


    end subroutine point_set

    function get_face(m_pnts)
      integer :: get_face, face
      integer, intent(in) :: m_pnts(2)
      select case(m_pnts(1))
      case(1)
         face=3
      case(2)
         face=1
      case(3)
         face=2
      end select
      get_face=face
    end function get_face

end module MLCM_operator
