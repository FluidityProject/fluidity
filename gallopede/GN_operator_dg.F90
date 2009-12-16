#include "fdebug.h"

module GN_operator
  
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
  use fldebug
  
  implicit none

  public :: assemble_GN_alpha_operator, assemble_GN_D_operator, add_penalty_term, getmass, get_vels, get_momentum

  private
  
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
      
contains
  
  subroutine assemble_GN_alpha_operator(mom11,mom12,mom21,mom22, &
       rhs1,rhs2,D,m1,m2, &
       mesh, &
       u1,u2)
    implicit none
    type(csr_matrix), intent(out) :: mom11, mom12, mom21, mom22
    real, dimension(n_vels), intent(inout) :: rhs1, rhs2, m1, m2
    real, dimension(n_verts), intent(in) :: D
    type(dg_mesh) :: mesh
    real, dimension(n_vels), intent(in), optional :: u1,u2

    !locals
    integer :: ele, globi, globj, iloc, jloc, gi,ni,ele_2,bdy_i,bcnt
    integer :: i,j
    real :: kmat, kmat11, kmat12, kmat21, kmat22
    real, dimension(mesh%nu%ngi) :: dloc
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
    integer, dimension(:), pointer :: u_ele, u_ele_2, X_ele, X_ele_2

    real :: val1, llength

    real, dimension(3,3) :: bcmat

    ! Local node number map for big NC element.
    integer, dimension(9) :: local_glno

    ! Local integration matrices for big NC element.
    real, dimension(2,mesh%nh%loc,mesh%nu%loc+3*mesh%nu_f%loc) :: B
    real, dimension(2,2,mesh%nu%loc+3*mesh%nu_f%loc,mesh%nu%loc+3*mesh%nu_f%loc) :: BQB
    real, dimension(mesh%nh%loc,mesh%nh%loc) :: Qd
    real, dimension(mesh%nu%loc,mesh%nu%loc) :: Qu
    real, dimension(2,mesh%nh_f%loc,mesh%nu_f%loc) :: Q_surf

    real :: sgn

    ewrite(2,*)("subroutine assemble_GN_alpha_operator")

    ewrite(2,*)("Allocating memory for neigh");
    allocate(neigh(row_length(mesh%bdy_list,1)))

    call zero(mom11)
    call zero(mom12)
    call zero(mom21)
    call zero(mom22)

    rhs1 = 0.
    rhs2 = 0.

    ele_loop: do ele = 1, n_elements

       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)

       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)
       call transform_to_physical(ele_X, mesh%nu, m = mesh%nh, &
            dn_t = dnu_t, dm_t = dnh_t, &
            detwei = detwei)

       !volume integrals
       !get local mass matrix
       Qu=shape_shape(mesh%nu,mesh%nu,detwei)

       iloc_loop: do iloc = 1, mesh%nu%loc
          globi = u_ele(iloc)
          jloc_loop: do jloc = 1, mesh%nu%loc
             globj = u_ele(jloc)

             rhs1(globi) = rhs1(globi) + Qu(iloc,jloc)*m1(globj)
             rhs2(globi) = rhs2(globi) + Qu(iloc,jloc)*m2(globj)

             call addto(mom11,globi,globj,Qu(iloc,jloc))
             call addto(mom22,globi,globj,Qu(iloc,jloc))

          end do jloc_loop
       end do iloc_loop

       if(.true.) then

          QD = shape_shape(mesh%nh,mesh%nh,detwei)
          call invert(QD)

          ! First part of local numbering is for this element.
          local_glno(1:mesh%nu%loc)=u_ele
          local_glno(mesh%nu%loc+1:)=0

          !------------------------------------------------------------------
          ! Element internal integral.
          !------------------------------------------------------------------

          B=0.0

          B(:,:,1:mesh%nu%loc)= dshape_shape(dnh_t, mesh%nu, detwei)

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
             !if (ele_2==0) then
             !   bdy => offnods(ni,:)
             !else
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

             Q_surf=shape_shape_vector(mesh%nh_f,mesh%nu_f, &
                  detwei_f, normal)

             if(ele_2.ne.0) then
                ! Values of local node map
                local_glno(mesh%nu%loc+(ni-1)*mesh%nu_f%loc+1:3+ni*mesh%nu_f%loc)= &
                     u_ele_2(bdyu_2)
             end if

             sgn = -1
             if(normal(1,1)*gvec(1) + normal(1,2)*gvec(2) >0.0) sgn = 1

             !==============================================================
             if(ele_2.ne.0) then
                B(:,bdyh,bdyu)=B(:,bdyh,bdyu)-0.5*Q_surf
                B(:,bdyh,mesh%nu%loc+(ni-1)*mesh%nu_f%loc+1:mesh%nu%loc+ni*mesh%nu_f%loc)= &
                     -0.5*Q_surf
             end if
             !===============================================================
          end do neighbourloop

          BQB=0.0

          forall(i = 1:2, j = 1:2)
             BQB(i,j,:,:)=BQB(i,j,:,:) &
                  +matmul(transpose(matmul(Qd,B(i,:,:))),B(j,:,:))
          end forall

          ! Put the contribution for this element into the matrix.
          iloop: do iloc=1,9
             globi=local_glno(iloc)

             ! Exclude boundaries for some methods
             if (globi==0) cycle iloop

             jloop: do jloc=1,9
                globj=local_glno(jloc)

                ! Exclude boundaries.
                if (globj==0) cycle jloop

                ! Insert value in the matrix.
                if(.false.) then
                   call addto(mom11, globi, globj,alpha1*alpha1* &
                        (BQB(1,1,iloc,jloc) + BQB(2,2,iloc,jloc)))
                   call addto(mom22, globi, globj,alpha1*alpha1* &
                        (BQB(2,2,iloc,jloc) + BQB(1,1,iloc,jloc)))
                else

                   call addto(mom11, globi, globj,alpha1*alpha1* &
                        BQB(1,1,iloc,jloc))
                   call addto(mom12, globi, globj,alpha1*alpha1* &
                        BQB(1,2,iloc,jloc))
                   call addto(mom21, globi, globj,alpha1*alpha1* &
                        BQB(2,1,iloc,jloc))
                   call addto(mom22, globi, globj,alpha1*alpha1* &
                        BQB(2,2,iloc,jloc))
                end if
             end do jloop

          end do iloop

       end if

    end do ele_loop

    ewrite(2,*)("end subroutine assemble_GN_operator")

  end subroutine assemble_GN_alpha_operator

  subroutine assemble_GN_D_operator( D, mesh, &
       mom11,mom12,mom21,mom22, &
       rhs1,rhs2, &
       m1,m2, &
       u1,u2)
    implicit none
    type(csr_matrix), intent(out), optional :: mom11, mom12, mom21, mom22
    real, dimension(n_vels), intent(inout), optional :: rhs1, rhs2, m1, m2
    real, dimension(n_dens), intent(in) :: D
    type(dg_mesh) :: mesh
    real, dimension(n_vels), intent(in), optional :: u1,u2

    !locals
    integer :: ele, globi, globj, iloc, jloc, gi,ni,ele_2,bdy_i,bcnt
    integer :: i,j
    real :: kmat, kmat11, kmat12, kmat21, kmat22
    real, dimension(mesh%nu%ngi) :: dloc
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
    real, dimension(2,mesh%nh%loc,mesh%nu%loc+3*mesh%nu_f%loc) :: B
    real, dimension(2,2,mesh%nu%loc+3*mesh%nu_f%loc,mesh%nu%loc+3*mesh%nu_f%loc) :: BQB
    real, dimension(mesh%nh%loc,mesh%nh%loc) :: Qd
    real, dimension(mesh%nu%loc,mesh%nu%loc) :: Qu
    real, dimension(2,mesh%nh_f%loc,mesh%nu_f%loc) :: Q_surf

    ewrite(2,*)("subroutine assemble_GN_D_operator")

    ewrite(2,*)("Allocating memory for neigh");
    allocate(neigh(row_length(mesh%bdy_list,1)))

    if(present(mom11)) then
       assert(present(mom12))
       assert(present(mom21))
       assert(present(mom22))
       !call zero(mom11)
       !call zero(mom12)
       !call zero(mom21)
       !call zero(mom22)
    end if

    if(present(rhs1)) then
       assert(present(rhs2))
       !rhs1 = 0.
       !rhs2 = 0.
       assert(present(u1).or.present(m1))
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

    ele_loop: do ele = 1, n_elements

       u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)

       ele_X(1,:)=mesh%X(X_ele)
       ele_X(2,:)=mesh%Y(X_ele)
       call transform_to_physical(ele_X, mesh%nu, m = mesh%nh, &
            dn_t = dnu_t, dm_t = dnh_t, &
            detwei = detwei)

       !get dloc
       dloc = matmul(transpose(mesh%nh%n),D(h_ele))

       !volume integrals

       if(present(rhs1).and.present(m1)) then
          
          !get local mass matrix
          Qu=shape_shape(mesh%nu,mesh%nu,detwei)
          
          iloc_loop: do iloc = 1, mesh%nu%loc
             globi = u_ele(iloc)
             jloc_loop: do jloc = 1, mesh%nu%loc
                globj = u_ele(jloc)
                
                rhs1(globi) = rhs1(globi) + Qu(iloc,jloc)*m1(globj)
                rhs2(globi) = rhs2(globi) + Qu(iloc,jloc)*m2(globj)
                
             end do jloc_loop
          end do iloc_loop
       end if

       !get scaled mass matrix
       Qu=shape_shape(mesh%nu,mesh%nu,detwei * dloc )
       
       iloc_dloop: do iloc = 1, mesh%nu%loc
          globi = u_ele(iloc)
          jloc_dloop: do jloc = 1, mesh%nu%loc
             globj = u_ele(jloc)
             
             if(present(mom11)) then
                call addto(mom11,globi,globj,Qu(iloc,jloc))
                call addto(mom22,globi,globj,Qu(iloc,jloc))
             end if
             if(present(rhs1).and.present(u1)) then
                rhs1(globi) = rhs1(globi) + Qu(iloc,jloc)*u1(globj)
                rhs2(globi) = rhs2(globi) + Qu(iloc,jloc)*u2(globj)
             end if
          end do jloc_dloop
       end do iloc_dloop
    
       QD = shape_shape(mesh%nh,mesh%nh,3.0 * detwei / (Dloc**3) )
       call invert(QD)

       ! First part of local numbering is for this element.
       local_glno(1:mesh%nu%loc)=u_ele
       local_glno(mesh%nu%loc+1:)=0

       !------------------------------------------------------------------
       ! Element internal integral.
       !------------------------------------------------------------------

       B=0.0

       B(:,:,1:mesh%nu%loc)= dshape_shape(dnh_t, mesh%nu, detwei)

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
          !if (ele_2==0) then
          !   bdy => offnods(ni,:)
          !else
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

          Q_surf=shape_shape_vector(mesh%nh_f,mesh%nu_f, &
               detwei_f, normal)

          if(ele_2.ne.0) then
             ! Values of local node map
             local_glno(mesh%nu%loc+(ni-1)*mesh%nu_f%loc+1:3+ni*mesh%nu_f%loc)= &
                  u_ele_2(bdyu_2)
          end if

          !==============================================================
          if(ele_2.ne.0) then
             B(:,bdyh,bdyu)=B(:,bdyh,bdyu)-0.5*Q_surf
             B(:,bdyh,mesh%nu%loc+(ni-1)*mesh%nu_f%loc+1:mesh%nu%loc+ni*mesh%nu_f%loc)= &
                  -0.5*Q_surf
          end if
          !===============================================================
       end do neighbourloop

       BQB=0.0

       forall(i = 1:2, j = 1:2)
          BQB(i,j,:,:)=BQB(i,j,:,:) &
               +matmul(transpose(matmul(Qd,B(i,:,:))),B(j,:,:))
       end forall

       ! Put the contribution for this element into the matrix.
       iloop: do iloc=1,9
          globi=local_glno(iloc)

          ! Exclude boundaries for some methods
          if (globi==0) cycle iloop

          jloop: do jloc=1,9
             globj=local_glno(jloc)

             ! Exclude boundaries.
             if (globj==0) cycle jloop

             if(present(mom11)) then
                ! Insert value in the matrix.
                call addto(mom11, globi, globj,BQB(1,1,iloc,jloc))
                call addto(mom12, globi, globj,BQB(1,2,iloc,jloc))
                call addto(mom21, globi, globj,BQB(2,1,iloc,jloc))
                call addto(mom22, globi, globj,BQB(2,2,iloc,jloc))
             end if
             if(present(rhs1).and.present(u1)) then
                rhs1(globi) = rhs1(globi) + BQB(1,1,iloc,jloc) * u1(globj) + &
                     BQB(1,2,iloc,jloc) * u2(globj)
                rhs2(globi) = rhs2(globi) + BQB(2,1,iloc,jloc) * u1(globj) + &
                     BQB(2,2,iloc,jloc) * u2(globj)                
             end if

          end do jloop

       end do iloop

    end do ele_loop

    ewrite(2,*)("end subroutine assemble_GN_D_operator")

  end subroutine assemble_GN_D_operator

  subroutine add_penalty_term(mom11,mom12,mom21,mom22,u1,u2,mesh)
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

    ewrite(2,*)("Subroutine add_penalty_term")

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

    ewrite(2,*)("END Subroutine add_penalty_term")

  end subroutine add_penalty_term

  subroutine getmass(mass,n1,evlist1,n2,evlist2)
    type(csr_matrix), intent(inout) :: mass
    type(element_type), intent(in) :: n1, n2
    integer, dimension(:), target :: evlist1, evlist2

    !locals
    integer :: ele, iloc, jloc, gi, globi, globj,totele
    integer, dimension(:), pointer :: vlist1,vlist2
    real :: kmat

    ASSERT(n1%ngi==n2%ngi)
    totele = size(evlist1)/n1%loc

    call zero(mass)

    do ele = 1, totele
       vlist1 => evlist1((ele-1)*n1%loc+1:ele*n1%loc)
       vlist2 => evlist2((ele-1)*n2%loc+1:ele*n2%loc)

       do iloc = 1,n1%loc
          globi = vlist1(iloc)
          do jloc = 1,n2%loc
             globj = vlist2(jloc)
             kmat = 0.
             do gi = 1,n1%ngi
                kmat = kmat + n1%n(iloc,gi)*n2%n(jloc,gi)
             end do

             call addto(mass,globi,globj,kmat)

          end do
       end do
    end do

  end subroutine getmass

  subroutine get_vels(Mesh, D, u, m, u0)
    type(dg_mesh), intent(in) :: Mesh
    real, dimension(n_dens), intent(in) :: D
    real, dimension(n_vels*2), intent(out) :: u
    real, dimension(n_vels*2), intent(in), optional :: u0
    real, dimension(n_vels*2), intent(in), target :: m

    !locals
    PetscErrorCode :: ierr
    KSPType :: ksp_type
    PCType :: pc_type

    real, dimension(:), allocatable, target :: urhs
    real, dimension(:), pointer :: urhs1, urhs2, m1, m2
    type(csr_matrix) :: mom11,mom12,mom21,mom22
    type(block_csr_matrix) :: mommat

    ewrite(2,*)("subroutine get_vels")

    ewrite(3,*)(sum(m))
    ewrite(3,*)(sum(D))

    ewrite(2,*)("allocating memory")
    allocate( urhs(N_vels*2) )
    urhs1 => urhs(1:N_vels)
    urhs2 => urhs(N_vels+1:N_vels*2)

    urhs = 0.

    m1 => m(1:N_vels)
    m2 => m(N_vels+1:N_vels*2)

    ewrite(2,*)("cloning Mass")
    mommat = block_clone(mesh%Mass_u, (/2, 2/))
    mom11 = block(mommat,1,1)
    mom12 = block(mommat,1,2)
    mom21 = block(mommat,2,1)
    mom22 = block(mommat,2,2)

    call assemble_GN_D_operator( D, mesh, &
         mom11=mom11,mom12=mom12,mom21=mom21,mom22=mom22, &
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

    call gallopede_block_solve(u, mommat, urhs, &
         ksp_type, pc_type, 1.0e-10, 500)

    call deallocate( mommat )
    deallocate( urhs )

    ewrite(2,*)("END subroutine get_vels")

  end subroutine get_vels

  subroutine get_momentum(Mesh, D, u1, u2, m1, m2)
    type(dg_mesh), intent(in) :: Mesh
    real, dimension(n_dens), intent(in) :: D
    real, dimension(n_vels), intent(in) :: u1, u2
    real, dimension(n_vels), intent(out) :: m1, m2

    !locals
    integer :: ele
    integer, dimension(:), pointer :: u_ele, X_ele
    real, dimension(2,3) :: ele_X
    real, dimension(mesh%nu%ngi) :: detwei
    real, dimension(mesh%nu%loc,mesh%nu%loc) :: Qu
    real, dimension(:), allocatable :: rhs1, rhs2

    allocate( rhs1(n_vels) )
    allocate( rhs2(n_vels) )

    m1 = 0.
    m2 = 0.

    rhs1 = 0.
    rhs2 = 0.

    call assemble_GN_D_operator( D, mesh, &
         rhs1=rhs1,rhs2=rhs2, &
         u1=u1,u2=u2)

    if(.true.) then
       ele_loop: do ele = 1, n_elements
          
          u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
          X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
          
          ele_X(1,:)=mesh%X(X_ele)
          ele_X(2,:)=mesh%Y(X_ele)
          call transform_to_physical(ele_X, mesh%nu, detwei = detwei)
          
          Qu=shape_shape(mesh%nu,mesh%nu,detwei)
          call invert(Qu)
          
          m1(u_ele) = matmul(Qu,rhs1(u_ele))
          m2(u_ele) = matmul(Qu,rhs2(u_ele))
          
       end do ele_loop
       
       else 
          m1 = rhs1
          m2 = rhs2
       end if

    deallocate( rhs1 )
    deallocate( rhs2 )

  end subroutine get_momentum

end module GN_operator
