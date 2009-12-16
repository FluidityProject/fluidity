#include "fdebug.h"

module helmholtz

  !a P1nc solver for helmholtz equation

  !coded by Colin Cotter and David Ham May-June 2006
  !for last modification see CVS
  
  !links against libdfluidity
  use elements
  use sparse_tools
  use quadrature
  use global_numbering
  use shape_functions
  use global_parameters_gallopede
  use adjacency_lists
  use transform_elements
  use dgtools
  use FEtools
  use Solvers
  use fldebug

  implicit none

contains
  
  subroutine Assembl_helmholtz_eqn(Helm,X,Y,EVList_h,EVList_u, &
       m1,m2,u1_rhs,u2_rhs, &
       nu,nh,nh_f, &
       b_seg_list,b_seg_nh_lno,nc,alpha1)

    implicit none

    !type(block_csr_matrix), intent(out)::Helm
    type(csr_matrix), intent(out)::Helm
    type(csr_matrix), intent(in)::b_seg_list
    real, intent(in)::X(:),Y(:)
    integer, dimension(:), target, intent(in)::EVList_h,EVList_u
    real, intent(in)::m1(:),m2(:)
    real, intent(out)::u1_rhs(:),u2_rhs(:)
    type(element_type) :: nu, nh, nh_f   ! Shape Functions
    integer, target, intent(in) :: b_seg_nh_lno(entries(b_seg_list)*nh_f%loc)
    logical, intent(in) :: nc
    real, intent(in) :: alpha1

    !local variables

    integer::ele,iloc,jloc,gi,globi,globj,ni,ele_2,b_seg_i
    real, dimension(2,nh%loc) :: ele_X, ele_X_2
    real, dimension(2,nh_f%loc) :: ele_Xf, ele_Xf_2
    real ::  dnu_t(nu%loc,nu%ngi,2), dnu_t_2(nu%loc,nu%ngi,2)
    real, dimension(nh%loc,nh%ngi,2) :: dnh_t
    integer, dimension(:), pointer :: b_seg, b_seg_2
    real :: detwei(nh%ngi), detwei_2(nh%ngi), detwei_f(nh_f%ngi)
    integer :: cj,i,j,test
    real :: normal(nh_f%ngi,2)
    real :: kmat

    ! List of neighbours of current element.
    integer, dimension(:), allocatable :: neigh

    ! Local element information
    integer, dimension(:), pointer :: u_ele, u_ele_2, h_ele, h_ele_2

    real :: val1, llength

    real, dimension(2,3,2,2) :: bcmatcg
    real, dimension(3,3) :: bcmatnc
    real, dimension(3,3) :: cg2ncmat
    real, dimension(3,3) :: nc2cgmat

    ! Local node number map for big NC element.
    integer, dimension(9) :: local_glno

    ! Local integration matrices for big NC element.
    real, dimension(2,3,9) :: B
    real, dimension(9,9) :: BQB
    real, dimension(3,3) :: Q 
    real, dimension(2,2,2) :: Q_surf_cg
    real, dimension(2,3,3) :: Q_surf_nc

    MSG("subroutine Assembl_helmholtz_eqn");

    nc2cgmat(:,1) = (/0.0,0.5,0.5/) 
    nc2cgmat(:,2) = (/0.5,0.0,0.5/)
    nc2cgmat(:,3) = (/0.5,0.5,0.0/)

    cg2ncmat(:,1) = (/-1.0, 1.0,  1.0/)
    cg2ncmat(:,2) = (/1.0, -1.0,  1.0/)
    cg2ncmat(:,3) = (/1.0,  1.0, -1.0/)

    ! Initial (and for meshes with only one shape of element, final) size
    ! of neigh.

    MSG("Allocating memory for neigh");
    allocate(neigh(row_length(b_seg_list,1)))

    MSG("Zeroing Helm");
    call zero(Helm)

    MSG("Zeroing u1_rhs u2_rhs");
    u1_rhs = 0.
    u2_rhs = 0.

    MSG("Elements loop");
    element_loop: do ele = 1, N_Elements
       u_ele=>EVList_u((ELE-1)*Nu%LOC+1:ELE*Nu%LOC)
       h_ele=>EVList_h((ELE-1)*Nh%LOC+1:ELE*Nh%LOC)

       ! Locations of local vertices.
       ele_X(1,:)=X(h_ele)
       ele_X(2,:)=Y(h_ele)

       ! Transform derivatives and weights into physical space.
       if(nc) then
          call transform_to_physical(ele_X, nh, m = nu, dm_t = dnu_t, &
               detwei = detwei)
       else
          call transform_to_physical(ele_X, nh, dn_t = dnh_t, detwei = detwei)
       end if

       ! Construct layer depth at the gauss points
       !hloc = 0.
       !do iloc = 1, nh%loc
       !   do gi = 1,nh%ngi
       !      hloc(gi) = hloc(gi) + nh%n(iloc,gi)*h(h_ele(iloc))
       !   end do
       !end do

       !-------------------------------------------------------------------
       ! Element internal integral.
       !-------------------------------------------------------------------

       if(nc) then

          !----------------------------------------------------------------
          ! Start of nc
          !----------------------------------------------------------------

          iloc_loop_nc: do iloc=1,nu%loc
             globi=u_ele(iloc)

             jloc_loop_nc: do jloc=1,nu%loc
                globj=u_ele(jloc)

                kmat = 0.

                gi_loop_nc: do gi=1,nu%ngi

                   !matrix contributions

                   kmat = kmat + detwei(gi)*nu%n(iloc,gi)*nu%n(jloc,gi)

                   !right-hand side
                   u1_rhs(globi) = u1_rhs(globi) + detwei(gi)* &
                        nu%n(iloc,gi)*nu%n(jloc,gi)*m1(globj)

                   u2_rhs(globi) = u2_rhs(globi) + detwei(gi)* &
                        nu%n(iloc,gi)*nu%n(jloc,gi)*m2(globj)

                end do gi_loop_nc

                ! Insert kmat in the matrix.

                call addto(Helm,globi,globj, kmat)

             end do jloc_loop_nc
          end do iloc_loop_nc

          !get local mass matrix
          Q=shape_shape(nu,nu,detwei)

          do iloc = 1,3
                Q(iloc,iloc) = 1.0/Q(iloc,iloc)
          end do

          ! First part of local numbering is for this element.
          local_glno(1:3)=u_ele
          local_glno(4:)=0

          !-------------------------------------------------------------------
          ! Element internal integral.
          !-------------------------------------------------------------------

          B=0.0
          B(:,:,1:3)= shape_dshape(nu, dnu_t, detwei)

          if (size(neigh)/=row_length(b_seg_list,ele)) then
             deallocate(neigh)
             allocate(neigh(row_length(b_seg_list,ele)))
          end if

          neigh=row_m(b_seg_list,ele)

          !loop over neighbours of ele

          neighbourloop: do ni=1,size(neigh)
             ele_2=neigh(ni)

             ! Skip external b_segs (neumann bcs)
             if (ele_2==0) cycle neighbourloop

             u_ele_2=>EVList_u((ELE_2-1)*nu%LOC+1:ELE_2*nu%LOC)
             h_ele_2=>EVList_h((ELE_2-1)*Nh%LOC+1:ELE_2*Nh%LOC)
             !x_ele_2=>EVList_h((ELE_2-1)*Nh%LOC+1:ELE_2*Nh%LOC)

             b_seg_i=ival(b_seg_list, ele, ele_2)
             b_seg=>b_seg_nh_lno((b_seg_i-1)*2+1:b_seg_i*2)

             b_seg_i=ival(b_seg_list, ele_2, ele)
             b_seg_2=>b_seg_nh_lno((b_seg_i-1)*2+1:b_seg_i*2)

             ! Locations of local vertices.
             ele_X_2(1,:)=X(h_ele_2)
             ele_X_2(2,:)=Y(h_ele_2)

             ! Locations of b_seg vertices.
             ele_Xf=ele_X(:,b_seg)

             ! Change of coordinates in second element.
             call transform_to_physical(ele_X_2, nh, m = nu, dm_t = dnu_t_2, &
                  detwei = detwei_2)
             call transform_face_to_physical(ele_X, ele_Xf, nh, nh_f, &
                  detw_f = detwei_f,normal = normal)

             ! Values of local node map
             local_glno(3+(ni-1)*2+1:3+ni*2)= u_ele_2(b_seg_2)

             Q_surf_cg=shape_shape_vector(nh_f,nh_f,detwei_f, normal)

             ! Interior face part of average.
             do i = 1,2
                Q_surf_nc(i,:,:)=matmul(cg2ncmat(:,b_seg), &
                     matmul(Q_surf_cg(i,:,:), &
                  cg2ncmat(b_seg,:)))
             end do
             B(:,b_seg,b_seg)=B(:,b_seg,b_seg)-0.5*Q_surf_nc(:,b_seg,b_seg)

             ! Exterior face part of average
             do i = 1,2
                Q_surf_nc(i,:,:)=matmul(cg2ncmat(:,b_seg), &
                     matmul(Q_surf_cg(i,:,:), &
                     cg2ncmat(b_seg_2,:)))
             end do
             B(:,b_seg,3+(ni-1)*2+1:3+ni*2)=+0.5*Q_surf_nc(:,b_seg,b_seg_2)

          end do neighbourloop

          BQB=0.0

          do i=1,2
             BQB=BQB+matmul(matmul(transpose(B(i,:,:)),Q),B(i,:,:))
          end do

          ! Put the contribution for this element into the matrix.
          iloop: do iloc=1,9
             globi=local_glno(iloc)

             ! Exclude boundaries.
             if (globi==0) cycle iloop

             jloop: do jloc=1,9
                globj=local_glno(jloc)

                ! Exclude boundaries.
                if (globj==0) cycle jloop

                ! Insert value in the matrix.
                call addto(Helm, globi, globj,alpha1*alpha1*BQB(iloc,jloc))
             end do jloop

          end do iloop

          !=================================================================
          !end of nc
          !=================================================================

       else

          iloc_loop: do iloc=1,nh%loc
             globi=h_ele(iloc)

             jloc_loop: do jloc=1,nh%loc
                globj=h_ele(jloc)

                kmat = 0.

                gi_loop: do gi=1,nh%ngi

                   !matrix contributions

                   kmat = kmat + detwei(gi)*(nh%n(iloc,gi)*nh%n(jloc,gi) + &
                        alpha1*alpha1*dot_product(dnh_t(iloc,gi,1:2), &
                        dnh_t(jloc,gi,1:2)) )

                   !right-hand side
                   u1_rhs(globi) = u1_rhs(globi) + detwei(gi)* &
                        nh%n(iloc,gi)*nh%n(jloc,gi)*m1(globj)

                   u2_rhs(globi) = u2_rhs(globi) + detwei(gi)* &
                        nh%n(iloc,gi)*nh%n(jloc,gi)*m2(globj)

                end do gi_loop

                ! Insert kmat in the matrix.

                call addto(Helm,globi,globj, kmat)

             end do jloc_loop
          end do iloc_loop

       end if

    end do element_loop

    deallocate( neigh )
    MSG("END subroutine Assembl_helmholtz_eqn");

  end subroutine Assembl_helmholtz_eqn

  subroutine solve_helmholtz_eqn(x, A, b, preerr, prenoi)
    use parallel_tools
    ! routine to encapsulate the call to the solver.
    ! Solves Ax=b with initial condition for x passed in
    real, dimension(:), intent(inout) :: X
    type(csr_matrix), intent(in) :: A
    real, dimension(:), intent(in) :: b
    real, intent(in) :: preerr
    integer, intent(in) :: prenoi

    ! Variables used by the conjugate gradient solver.
    integer :: imatst, timm, tmisou, kits, para, tnoit1
    real :: trelax, teror1
    logical :: momsym
    real, dimension(:,:), allocatable :: workv
    real, dimension(:), allocatable :: work1d
    real :: rdum(0)
    integer :: iidum(0)
    integer :: nonods

    imatst=0
    timm=-5
    tmisou=1
    trelax=1.0
    momsym=.false.
    teror1=preerr
    tnoit1=prenoi
    iidum=0
    allocate(workv(size(A%val), 5))
    allocate(work1d(4*size(A%val)))

    IF(IsParallel()) THEN
       para = 1
    else
       para=0
       nonods = size(X)
    end if
    MSG("going inot SOLCG------")
    call pminmx(b,size(b),'******difvec  ')

    call solcg(X, b, size(X), size(X), nonods, .true.,  &
         A%val, A%findrm, A%colm, &
         size(A%colm), size(A%colm), &
         0, kits)

    MSG("just out of SOLCG")
    CHECK(kits)

  end subroutine solve_helmholtz_eqn

end module helmholtz
