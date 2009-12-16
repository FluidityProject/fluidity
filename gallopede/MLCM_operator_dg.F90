#include "fdebug.h"

module MLCM_operator
  
  use elements
  use sparse_tools
  use quadrature
  use global_numbering
  use shape_functions
  use global_parameters_gallopede
  use adjacency_lists
  use transform_elements
  use dgtools
  use text_io
  use data_structures
  use fetools
  use vector_tools
  use gallopede_solvers
  use vtk_io
  use mesh_tools
  use multilayer_tools
  use fldebug
  
  implicit none

  public :: assemble_MLCM_operator, get_vels_layer

  private
  
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
      
contains
  
  subroutine assemble_MLCM_operator( D, mesh, &
       bottom, rho, &
       mommat, &
       rhs1,rhs2, &
       m1,m2, &
       u1,u2)
    implicit none
    type(block_csr_matrix), intent(out), optional :: mommat
    type(layer), dimension(N_Layers), intent(inout), optional :: rhs1, rhs2, m1, m2
    type(layer), dimension(N_Layers), intent(in) :: D
    type(dg_mesh) :: mesh
    type(layer), dimension(n_layers), intent(in), optional :: u1,u2
    real, dimension(N_verts), intent(in) :: bottom
    real, dimension(N_layers), intent(in) :: rho

    !locals
    integer :: ele, globi, globj, iloc, jloc, gi,ni,ele_2,bdy_i,bcnt
    integer :: i,j,k, iX, jX
    real :: kmat, kmat11, kmat12, kmat21, kmat22
    real, dimension(N_layers,mesh%nh%loc) :: dloc
    real, dimension(N_layers,mesh%nh%ngi) :: dlocgi
    real, dimension(N_layers,mesh%nh_f%ngi) :: dlocgi_f
    real, dimension(2,N_layers,mesh%nu%ngi) :: gradblocgi
    real, dimension(2,3) :: ele_X, ele_X_2
    real, dimension(2,2) :: ele_Xf, ele_Xf_2
    real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t, dnu_t_2
    real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t, dnh_t_2
    integer, dimension(:), pointer :: bdyu, bdyu_2, bdyh, bdyh_2
    real, dimension(mesh%nu%ngi) :: detwei, detwei_2 
    real, dimension(mesh%nu_f%ngi) :: detwei_f
    real, dimension(2,mesh%nu_f%ngi) :: normal

    ! List of neighbours of current element.
    integer, dimension(:), allocatable :: neigh

    ! Local element information
    integer, dimension(:), pointer :: u_ele, u_ele_2, X_ele, X_ele_2, h_ele

    real :: val1, llength


    real, dimension(3,3) :: bcmat

    ! Local node number map for big NC element.
    integer, dimension(9) :: local_glno

    ! Local integration matrices for big NC element.
    real, dimension(N_Layers,mesh%nh%loc,mesh%nh%loc) :: Dmat
    real, dimension(2,mesh%nh%loc,mesh%nu%loc+3*mesh%nu_f%loc) :: B
    real, dimension(2,N_Layers,N_Layers,mesh%nh%loc,mesh%nu%loc+3*mesh%nu_f%loc) :: C
    real, dimension(2,2,N_Layers,N_Layers,mesh%nu%loc+3*mesh%nu_f%loc,mesh%nu%loc+3*mesh%nu_f%loc) :: BQB
    real, dimension(mesh%nh%loc,mesh%nh%loc) :: Qd
    real, dimension(mesh%nu%loc,mesh%nu%loc) :: Qu
    real, dimension(2,mesh%nh_f%loc,mesh%nu_f%loc) :: Q_surf
    logical :: log1

    ewrite(2,*)("subroutine assemble_MLCM_operator")

    log1=mommat%blocks(1)==2*N_Layers
    assert(log1)
    log1=mommat%blocks(2)==2*N_Layers
    assert(log1)

    ewrite(2,*)("Allocating memory for neigh");
    allocate(neigh(row_length(mesh%bdy_list,1)))

    if(present(rhs1)) then
       assert(present(rhs2))
       !rhs1 = 0.
       !rhs2 = 0.
       log1 = present(u1).or.present(m1)
       assert( log1 )
    end if

    if(present(m1)) then
       assert(present(m2))
       assert(.not.present(u1))
       assert(.not.present(u2))
    end if

    if(present(u1)) then
       assert(present(u2))
       assert(.not.present(m1))
       assert(.not.present(m2))
    end if

    ewrite(2,*)("ele loop")
    ele_loop: do ele = 1, n_elements

       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)

       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)
       call transform_to_physical(ele_X, mesh%nu, m = mesh%nh, &
            dn_t = dnu_t, dm_t = dnh_t, &
            detwei = detwei)

       !get dlocgi
       
       dloc = 0.
       dlocgi = 0.
       do i = 1, N_Layers
          dloc(i,:) = D(i)%l(h_ele)
          do iloc = 1, mesh%nh%loc
             dlocgi(i,:) = dlocgi(i,:) + &
                  mesh%nh%n(iloc,:) * dloc(i,iloc)
          end do
       end do

       !get gradblocgi
       gradblocgi = 0.
       do i = 1, N_Layers
          do iloc = 1,mesh%nu%loc
             do iX = 1,2
                gradblocgi(iX,i,:) = gradblocgi(iX,i,:) + &
                     dnu_t(iloc,:,iX) * bottom(X_ele(iloc))
             end do
          end do
       end do

       !volume integrals
       if(present(rhs1).and.present(m1)) then
          
          !get local mass matrix
          Qu=shape_shape(mesh%nu,mesh%nu,detwei)
          
          iloc_loop: do iloc = 1, mesh%nu%loc
             globi = u_ele(iloc)
             jloc_loop: do jloc = 1, mesh%nu%loc
                globj = u_ele(jloc)

                forall(i = 1:N_Layers)
                   rhs1(i)%l(globi) = rhs1(i)%l(globi) + &
                        Qu(iloc,jloc)*m1(i)%l(globj)
                   rhs2(i)%l(globi) = rhs2(i)%l(globi) + &
                        Qu(iloc,jloc)*m2(i)%l(globj)
                end forall
             end do jloc_loop
          end do iloc_loop
       end if

       !get scaled mass matrix

       do i = 1, N_Layers
          Qu=shape_shape(mesh%nu,mesh%nu,detwei * dloc(i,:) )
          
          iloc_dloop: do iloc = 1, mesh%nu%loc
             globi = u_ele(iloc)
             jloc_dloop: do jloc = 1, mesh%nu%loc
                globj = u_ele(jloc)
                   
                if(present(mommat)) then
                   call addto(mommat,2*(i-1)+1,2*(i-1)+1, &
                        globi,globj,rho(i)*Qu(iloc,jloc))
                   call addto(mommat,2*(i-1)+2,2*(i-1)+2, &
                        globi,globj,rho(i)*Qu(iloc,jloc))
                end if
                if(present(rhs1).and.present(u1)) then
                   rhs1(i)%l(globi) = rhs1(i)%l(globi) + &
                        rho(i)*Qu(iloc,jloc)*u1(i)%l(globj)
                   rhs2(i)%l(globi) = rhs2(i)%l(globi) + &
                        rho(i)*Qu(iloc,jloc)*u2(i)%l(globj)
                end if
             end do jloc_dloop
          end do iloc_dloop
          
       end do

       QD(:,:) = shape_shape(mesh%nh,mesh%nh,detwei)

       ! First part of local numbering is for this element.
       local_glno(1:mesh%nu%loc)=u_ele
       local_glno(mesh%nu%loc+1:)=0

       !------------------------------------------------------------------
       ! Element internal integral.
       !------------------------------------------------------------------

       B=0.0
       C=0.0

       B(:,:,1:mesh%nu%loc) = -dshape_shape(dnh_t, mesh%nu, detwei)
       do i = 1, N_Layers
          if(i < N_Layers) then
             do j = i+1, N_Layers
                C(i,j,:,:,1:mesh%nu%loc) = - &
                     dshape_shape(dnh_t, mesh%nu, &
                     detwei * dlocgi(j,:) )
             end do
          end if
       end do

       !effect of bottom boundary
       do i = 1, N_Layers
          C(i,i,:,:,1:mesh%nu%loc) = &
             shape_shape_vector(mesh%nh,mesh%nu, &
             detwei, gradblocgi(:,i,:) )
       end do

       if (size(neigh)/=row_length(mesh%bdy_list,ele)) then
          ewrite(2,*)("reallocating neigh")
          deallocate(neigh)
          allocate(neigh(row_length(mesh%bdy_list,ele)))
       end if

       neigh=row_m(mesh%bdy_list,ele)

       !loop over neighbours of ele

       neighbourloop: do ni=1,size(neigh)
          ele_2=neigh(ni)

          if (ele_2==0) cycle neighbourloop

          bdy_i=ival(mesh%bdy_list, ele, ele_2)
          bdyu=>mesh%bdy_nu_lno((bdy_i-1)*mesh%nu_f%loc+1: &
               bdy_i*mesh%nu_f%loc)
          bdyh=>mesh%bdy_nh_lno((bdy_i-1)*mesh%nh_f%loc+1: &
               bdy_i*mesh%nh_f%loc)

          bdy_i=ival(mesh%bdy_list, ele_2, ele)
          bdyu_2=> mesh%bdy_nu_lno((bdy_i-1)*mesh%nu_f%loc+1: &
               bdy_i*mesh%nu_f%loc)
          bdyh_2=> mesh%bdy_nh_lno((bdy_i-1)*mesh%nh_f%loc+1: &
               bdy_i*mesh%nh_f%loc)

          u_ele_2=> &
               mesh%EVList_u((ELE_2-1)*mesh%nu%LOC+1: &
               ELE_2*mesh%nu%LOC)
          !end if

          ! Locations of bdy vertices.
          ele_Xf=ele_X(:,bdyu)

          ! Change of coordinates on face
          call transform_bdy_to_physical(ele_X, ele_Xf, &
               mesh%nu, mesh%nu_f, &
               detw_f = detwei_f,normal = normal)

          if(ele_2.ne.0) then
             ! Values of local node map
             local_glno(mesh%nu%loc+(ni-1)*mesh%nu_f%loc+1:3+ni*mesh%nu_f%loc)= &
                  u_ele_2(bdyu_2)
          end if

          dlocgi_f = 0.
          do i = 1, N_Layers
             do iloc = 1, mesh%nh_f%loc
                dlocgi_f(i,:) = dlocgi_f(i,:) + &
                     mesh%nh_f%n(iloc,:) * D(i)%l(h_ele(bdyh(iloc)))
             end do
          end do

          !==============================================================
          if(ele_2.ne.0) then

             Q_surf=shape_shape_vector(mesh%nh_f,mesh%nu_f, &
                  detwei_f, normal)
             
             B(:,bdyh,bdyu)=B(:,bdyh,bdyu)+0.5*Q_surf
             B(:,bdyh,mesh%nu%loc+(ni-1)*mesh%nu_f%loc+1:mesh%nu%loc+ni*mesh%nu_f%loc)= &
                  +0.5*Q_surf

             do i = 1, N_Layers
                if(i<N_Layers) then
                   do j = i+1, N_Layers
                      Q_surf=shape_shape_vector(mesh%nh_f,mesh%nu_f, &
                           detwei_f * dlocgi_f(j,:) , normal)
                      C(i,j,:,bdyh,bdyu)=C(i,j,:,bdyh,bdyu)+0.5*Q_surf
                      C(i,j,:,bdyh,mesh%nu%loc+(ni-1)*mesh%nu_f%loc+1:mesh%nu%loc+ni*mesh%nu_f%loc)= &
                      C(i,j,:,bdyh,mesh%nu%loc+(ni-1)*mesh%nu_f%loc+1:mesh%nu%loc+ni*mesh%nu_f%loc) &
                           +0.5*Q_surf
                   end do
                end if
             end do
          end if
          !===============================================================
       end do neighbourloop

       BQB=0.0

       !BB part
       dmat = 0.
       forall(iloc=1:mesh%nh%loc, i=1:N_Layers) 
          dmat(i,iloc,iloc) = rho(i) * dloc(i,iloc)**3/3.0
       end forall
       forall(i=1:2, iX = 1:2, jX = 1:2)
          BQB(i,i,iX,jX,:,:)=BQB(i,i,iX,jX,:,:) &
               +matmul(transpose(matmul(Qd, &
               matmul(dmat(i,:,:),B(iX,:,:)))),B(jX,:,:))
       end forall
       
       !BC part
       dmat = 0.
       forall(iloc=1:mesh%nh%loc, i=1:N_Layers) 
          dmat(i,iloc,iloc) = rho(i) * dloc(i,iloc)**2/2.0
       end forall
       forall(i=1:2, j=1:2, iX = 1:2, jX = 1:2)
          BQB(i,j,iX,jX,:,:)=BQB(i,j,iX,jX,:,:) &
               +matmul(transpose(matmul(Qd, &
               matmul(dmat(i,:,:),C(j,i,iX,:,:)))),B(jX,:,:))
       end forall
       
       !CB part
       forall(i=1:2, j=1:2, iX = 1:2, jX = 1:2)
          BQB(i,j,iX,jX,:,:)=BQB(i,j,iX,jX,:,:) &
               +matmul(transpose(matmul(Qd, &
               matmul(dmat(i,:,:),B(iX,:,:)))),C(i,j,jX,:,:))
       end forall
       
       !CC part
       dmat = 0.
       forall(iloc=1:mesh%nh%loc, i=1:N_Layers) 
          dmat(i,iloc,iloc) = rho(i) * dloc(i,iloc)
       end forall
       do i = 1,N_layers
          do j = 1,N_Layers
             do k = 1,N_Layers
                do iX = 1,2
                   do jX = 1,2
                      BQB(i,j,iX,jX,:,:)=BQB(i,j,iX,jX,:,:) &
                           +matmul(transpose(matmul(Qd, &
                           matmul(dmat(k,:,:),C(k,i,iX,:,:)))),C(k,j,jX,:,:))
                   end do
                end do
             end do
          end do
       end do

       ! Put the contribution for this element into the matrix.
       iloop: do iloc=1,9
          globi=local_glno(iloc)

          ! Exclude boundaries for some methods
          if (globi==0) cycle iloop

          jloop: do jloc=1,9
             globj=local_glno(jloc)

             ! Exclude boundaries.
             if (globj==0) cycle jloop

             if(present(mommat)) then
                ! Insert value in the matrix.
                do i = 1, N_Layers
                   do j = 1, N_Layers
                      call addto(mommat, (i-1)*2 + 1, (j-1)*2 + 1, &
                           globi, globj, &
                           BQB(i,j,1,1,iloc,jloc))
                      call addto(mommat, (i-1)*2 + 1, (j-1)*2 + 2, &
                           globi, globj, &
                           BQB(i,j,1,2,iloc,jloc))
                      call addto(mommat, (i-1)*2 + 2, (j-1)*2 + 1, &
                           globi, globj, &
                           BQB(i,j,2,1,iloc,jloc))
                      call addto(mommat, (i-1)*2 + 2, (j-1)*2 + 2, &
                           globi, globj, &
                           BQB(i,j,2,2,iloc,jloc))
                   end do
                end do
             end if
             if(present(rhs1).and.present(u1)) then
                do i = 1, N_Layers
                   do j = 1, N_Layers
                      rhs1(i)%l(globi) = rhs1(i)%l(globi) + &
                           BQB(i,j,1,1,iloc,jloc) * u1(j)%l(globj) + &
                           BQB(i,j,1,2,iloc,jloc) * u2(j)%l(globj)
                      rhs2(i)%l(globi) = rhs2(i)%l(globi) + &
                           BQB(i,j,2,1,iloc,jloc) * u1(j)%l(globj) + &
                           BQB(i,j,2,2,iloc,jloc) * u2(j)%l(globj)
                   end do
                end do
             end if

          end do jloop

       end do iloop

    end do ele_loop

    ewrite(2,*)("end subroutine assemble_MLCM_operator")

  end subroutine assemble_MLCM_operator

  subroutine add_penalty_term_layers(mom11,mom12,mom21,mom22,u1,u2,mesh)
    type(csr_matrix), intent(inout) :: mom11, mom12, mom21, mom22
    real, dimension(n_vels), intent(inout) :: u1, u2
    type(dg_mesh) :: mesh

    !locals
    integer :: ele, globi, globj, iloc, jloc, gi,ni,ele_2,bdy_i,bcnt
    integer :: i,j, globiE, globiF, globjE, globjF
    real :: kmat, kmat11, kmat12, kmat21, kmat22
    real, dimension(2,3) :: ele_X, ele_X_2
    real, dimension(2,2) :: ele_Xf
    integer, dimension(:), pointer :: bdy, bdy_2
    real, dimension(mesh%nu_f%ngi) :: detwei_f, udotnE,udotnF
    real :: udotn 
    real, dimension(2,mesh%nu_f%ngi) :: normal


    ! List of neighbours of current element.
    integer, dimension(:), allocatable :: neigh

    ! Local element information
    integer, dimension(:), pointer :: u_ele, u_ele_2, X_ele, X_ele_2

    real, dimension(mesh%nu_f%loc,mesh%nu_f%loc) :: Q_penE, Q_penF

    real :: val1, llength, sgn

    ewrite(2,*)("Subroutine add_penalty_term_layers")

    ewrite(3,*)(eta)

    if(eta>0.0) then

       ewrite(2,*)("Allocating memory for neigh");
       allocate(neigh(row_length(mesh%bdy_list,1)))

       ele_loop: do ele = 1, n_elements

          u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
          X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)

          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)

          if (size(neigh)/=row_length(mesh%bdy_list,ele)) then
             deallocate(neigh)
             allocate(neigh(row_length(mesh%bdy_list,ele)))
          end if

          neigh=row_m(mesh%bdy_list,ele)

          !loop over neighbours of ele

          neighbourloop: do ni=1,size(neigh)
             ele_2=neigh(ni)

             if (ele_2==0) cycle neighbourloop

             u_ele_2=> &
                  mesh%EVList_u((ELE_2-1)*mesh%nu%LOC+1:ELE_2*mesh%nu%LOC)
             X_ele_2=> mesh%EVList_X((ELE_2-1)*3+1:ELE_2*3)

             bdy_i=ival(mesh%bdy_list, ele, ele_2)
             bdy=>mesh%bdy_nu_lno((bdy_i-1)*2+1:bdy_i*2)

             bdy_i=ival(mesh%bdy_list, ele_2, ele)
             bdy_2=> mesh%bdy_nu_lno((bdy_i-1)*2+1:bdy_i*2)

             ! Locations of bdy vertices.
             ele_Xf=ele_X(:,bdy)

             ! Change of coordinates in second element.
             call transform_bdy_to_physical(ele_X, ele_Xf, &
                  mesh%nu, mesh%nu_f, &
                  detw_f = detwei_f,normal = normal)

             if(.true.) then
                !get udotn
                udotnE = 0.
                do iloc = 1,mesh%nu_f%loc
                   globi = u_ele(bdy(iloc))
                   udotnE = udotnE + &
                        mesh%nu_f%n(iloc,:) * ( &
                        u1(globi) * normal(1,:) + &
                        u2(globi) * normal(2,:) )
                end do

                udotnF = 0.
                do iloc = 1,mesh%nu_f%loc

                   globi = u_ele(bdy(iloc))
                   udotnF = udotnF + &
                        mesh%nu_f%n(iloc,:) * ( &
                        u1(globi) * normal(1,:) + &
                        u2(globi) * normal(2,:) )

                end do
             end if

             Q_penE = shape_shape(mesh%nu_f,mesh%nu_f, &
                  detwei_f * abs(udotnE))
             Q_penF = shape_shape(mesh%nu_f,mesh%nu_f, &
                  detwei_f * abs(udotnF))

             do iloc = 1,mesh%nu_f%loc

                globiE = u_ele(bdy(iloc))
                globiF = u_ele_2(bdy_2(iloc))

                do jloc = 1,mesh%nu_f%loc

                   globjE = u_ele(bdy(jloc))
                   globjF = u_ele_2(bdy_2(jloc))

                   call addto(mom11,globiE,globjE, &
                        eta*Q_penE(iloc,jloc))
                   call addto(mom22,globiE,globjE, &
                        eta*Q_penE(iloc,jloc))
                   call addto(mom11,globiE,globjF, &
                        -eta*Q_penF(iloc,jloc))
                   call addto(mom22,globiE,globjF, &
                        -eta*Q_penF(iloc,jloc))
                end do
             end do

          end do neighbourloop
       end do ele_loop

    end if

    ewrite(2,*)("END Subroutine add_penalty_term_layers")

  end subroutine add_penalty_term_layers

  subroutine get_vels_layer(Mesh, D, bottom, rho, u, m, u0)
    type(dg_mesh), intent(in) :: Mesh
    real, dimension(N_Layers*N_dens), intent(in), target :: D
    real, dimension(N_verts), intent(in) :: bottom
    real, dimension(N_Layers), intent(in) :: rho
    real, dimension(N_vels*2*N_Layers), intent(out) :: u
    real, dimension(n_vels*2*N_Layers), intent(in), optional :: u0
    real, dimension(n_vels*2*N_Layers), intent(in), target :: m

    !locals
    PetscErrorCode :: ierr
    KSPType :: ksp_type
    PCType :: pc_type

    real, dimension(:), allocatable, target :: urhs
    type(layer), dimension(:), pointer :: urhs1, urhs2, m1, m2, Dl
    type(block_csr_matrix) :: mommat
    integer :: i, j

    ewrite(2,*)("subroutine get_vels_layers")

    allocate( urhs1(N_layers) )
    allocate( urhs2(N_layers) )
    allocate( m1(N_layers) )
    allocate( m2(N_layers) )
    allocate( Dl(N_layers) )
    
    ewrite(2,*)("allocating memory")
    allocate( urhs(N_vels*2*N_Layers) )
    do i = 0, N_Layers-1
       urhs1(i+1)%l => urhs(1 + i*2*N_vels:N_vels + i*2*N_vels)
       urhs2(i+1)%l => urhs(N_vels+1 + i*2*N_vels:N_vels*2 + i*2*N_vels)
       m1(i+1)%l => m(1 + i*2*N_vels:N_vels + i*2*N_vels)
       m2(i+1)%l => m(N_vels + 1 + i*2*N_vels:2*N_vels + i*2*N_vels)
       Dl(i+1)%l => D(1 + i*N_dens:N_dens + i*N_dens)
    end do

    urhs = 0.

    ewrite(2,*)("cloning Mass")
    ewrite(3,*)(N_Layers)
    mommat = block_clone(mesh%Mass_u, (/2 * N_layers, 2 * N_layers /))

    call zero(mommat)
    ewrite(3,*)(sum(mesh%X))

    call assemble_MLCM_operator(Dl, mesh,bottom,rho, &
         mommat=mommat, &
         rhs1=urhs1,rhs2=urhs2, &
         m1=m1,m2=m2)

    !call add_penalty_term(mom11,mom12,mom21,mom22,u1,u2,mesh)

    ewrite(2,*)("solving GN operator equation")
    !solve GN equation
    ksp_type = KSPCG
    pc_type = PCICC

    if(present(u0)) then
       u = u0
    else
       u = 0.
    end if

    ewrite(3,*)(sum(URHS))

    !check out patricks solver code
    call gallopede_block_solve(u, mommat, urhs, &
         ksp_type, pc_type, 1.0e-10, 500)

    call deallocate( mommat )
    
    deallocate( urhs ) 
    deallocate( urhs1 )
    deallocate( urhs2 )
    deallocate( m1 )
    deallocate( m2 )
    deallocate( Dl )

    ewrite(2,*)("END subroutine get_vels_layers")

  end subroutine get_vels_layer

end module MLCM_operator
