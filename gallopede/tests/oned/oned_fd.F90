! One dimensional Finite Difference code for MLCM equation set
! Written by James Percival

! TO DO : FIX alternative BOUNDARY CONDITIONS
!         Proper IO


module MLCM_1d_globals
implicit none

integer, parameter :: mult=5
real :: t=0.0, dt=1.0,t_max=120000.0,kappa=1000.0, c=0.0
real :: dx=10.0*mult, g=10.0, t_out=20.0, alpha=0.0, alpha2=0.0
real :: kappa2=0.0,implicit_factor=0.00, alpha3=00.0
!integer, parameter :: n_extra=1000, n_vels=((1199-1)/2+n_extra)
!integer, parameter :: n_extra=500, n_vels=((550)+n_extra)
!integer, parameter :: n_extra=100, n_vels=((3400)/(2*mult)+n_extra)
integer, parameter :: n_extra=00, n_vels=((3400)/(2*mult)+n_extra)
integer, parameter :: n_dens=n_vels, n_layers=2
integer :: out=0
logical :: implicit_time_step=.false.
logical :: RIGID_LID=.false.

end module MLCM_1d_globals

module MLCM_1d_routines
use MLCM_1d_globals
implicit none

contains

 subroutine Initialize_Petsc()
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscviewer.h"
    PetscErrorCode :: ierr
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    CHKERRQ(ierr) 
  end subroutine Initialize_Petsc

  
  subroutine time_step(rho,D,m,u,b)
    real, dimension(n_layers), intent(in) :: rho
    real, dimension(n_dens,n_layers), intent(inout) :: D
    real, dimension(n_vels,n_layers), intent(out) :: u
    real, dimension(n_vels,n_layers), intent(inout) :: m
    real, dimension(n_dens) :: b
    real, dimension(4) :: k,rk
    real, dimension(n_vels,n_layers,4) :: fD
    real, dimension(n_vels,n_layers,4) :: fm
    real, dimension(n_dens,n_layers) ::Dw,Du
    real, dimension(n_vels,n_layers) :: mw

    integer :: i

    mw=m
    Dw=D

    fd=0.0; fm=0.0

! Picard Iteration "Backward Euler"


if (implicit_time_step) then
    do i=1,20
       call get_u(rho,mw,Dw,b,u)
       call get_D_flux(Dw,u,fd(:,:,1))
       call get_m_flux(rho,mw,u,Dw,b,fm(:,:,1))    
       
       Dw=d+dt*fD(:,:,1)
       mw=m+dt*fm(:,:,1)
    end do

    D=Dw
    m=mw

 end if

! Runge-Kutta

if (.not. implicit_time_step) then
    k=(/ 0.5 , 0.5 , 1.0, 0.0 /) 
    rk=(/1.0,2.0,2.0,1.0/)/6.0

    do i=1,4


       call get_u(rho,mw,Dw,b,u)
!       u(:,1)=b/400
!       u(:,2)=b/400
       call get_D_flux(Dw,u,fd(:,:,i))
       call get_m_flux(rho,mw,u,Dw,b,fm(:,:,i))    

       Dw=d+dt*k(i)*fD(:,:,i)
       mw=m+dt*k(i)*fm(:,:,i)

       print*, 'Sum flux m=', sum(abs(fm(:,:,i)))
       print*, 'Sum flux D=', sum(abs(fd(:,:,i)))
       print*, 'Sum u=', sum(abs(u)) 

    end do

    forall (i=1:4)
       fd(:,:,i)=fd(:,:,i)*rk(i)
       fm(:,:,i)=fm(:,:,i)*rk(i)
    end forall
    D=D+dt*sum(fd,3)
    m=m+dt*sum(fm,3)
end if


!     call get_u(rho,m,D,b,u)

       print*, 'Sum flux m=', sum(abs(m(:,:)))
       print*, 'Sum flux D=', sum(abs(sum(fd(:,:,:),3)))
       print*, 'Sum u=', sum(abs(u)) 

  end subroutine time_step

  subroutine get_u(rho,m,D,b,u)
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscviewer.h"

    real, dimension(n_layers), intent(in) :: rho
    real, dimension(n_dens), intent(in) :: b
    real, dimension(n_vels,n_layers), intent(in) :: m
    real, dimension(n_dens,n_layers), intent(in) :: D
    real, dimension(n_vels,n_layers), intent(inout) :: u
    real, dimension(n_vels*n_layers) :: temp
    integer, dimension(3) :: loc
    real, dimension(n_vels*n_layers,n_layers,3) :: mat
    integer :: i,j,k,kk


    KSPType :: ksp_type=KSPGMRES
    PCType :: pc_type=PCSOR
    PetscErrorCode :: ierr 
    Vec :: x, bvec
    Mat :: A
    PetscTruth:: pstrue
    PetscInt :: n,np(n_vels*n_layers), iterations
    KSP :: ksp ! Krylov subspace context
    PC  :: pc  ! Preconditioner context
    KSPConvergedReason :: reason

    n=n_vels*n_layers



    print*, "Creating vector x"
    call VecCreateSeq(PETSC_COMM_WORLD, n, x,ierr)
    CHKERRQ(ierr)
    print*, "Creating vector b"  
    call VecCreateSeq(PETSC_COMM_WORLD, n, bvec,ierr)

!    print*, "setting sizes"
!    call VecSetSizes(bvec, n, n, ierr)
!    call VecSetSizes(x, n, n, ierr)
    
    ! Zero the answer.
    print*, "Zeroing answer"
      call VecZeroEntries(x, ierr)
      do i=1,n_layers
         temp((i-1)*n_vels+1:i*n_vels)=u(:,i)
      end do
      call VecSetValues(x, n, (/ (i,i=0,n_vels*n_layers-1)/),&
           temp, INSERT_VALUES, ierr)
      print*, "Setting values in vecs"
      do i=1,n_layers
         temp((i-1)*n_vels+1:i*n_vels)=m(:,i)
      end do
      call VecSetValues(bvec, n, (/ (i,i=0,n_vels*n_layers-1)/),&
           temp, INSERT_VALUES, ierr)
      call VecAssemblyBegin(bvec, ierr)
      call VecAssemblyEnd(bvec, ierr)

      !Get the actual matrix content

      mat=0.0
      call get_matrix(rho,D,b,mat)
!      mat=0.0
!      do j=1,n_layers
!         mat((j-1)*n_vels+1:j*n_vels,j,2)=1.0
!      end do
      print*, "construct the matrix"

      np=3*n_layers

      call MatCreateSeqAIJ(MPI_COMM_SELF,n,&
           n,0,np, A, ierr)
      call MatZeroEntries(A, ierr)
      do i=0,n_vels-1
         do j=0,n_layers-1
            do k=0,n_layers-1
               loc(1:3)=i+(/-1,0,1/)
               where (loc .lt. 0) loc =loc +n_vels
               where (loc .ge. n_vels) loc = loc -n_vels
               call MatSetValues(A, 1, j*n_vels+i,3, &
                    loc+k*n_vels, &
                    mat(j*n_vels+i+1,k+1,:), INSERT_VALUES, ierr)
            end do
         end do
      end do

      print*, "assemble the matrix"

      call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

      print*, "Set up the solver context"

      call KSPCreate(MPI_COMM_SELF, ksp, ierr); CHKERRQ(ierr)
      call KSPSetInitialGuessNonZero(ksp,PETSC_TRUE,ierr)
      call KSPSetType(ksp, ksp_type, ierr); CHKERRQ(ierr)
      call KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN, ierr); 
      CHKERRQ(ierr)
      call KSPGetPC(ksp, pc, ierr); CHKERRQ(ierr)
      call PCSetType(pc, pc_type, ierr); CHKERRQ(ierr)
      call KSPSetTolerances(ksp,1.0e-16, 1.0e-20,     &
           PETSC_DEFAULT_DOUBLE_PRECISION,30000,ierr); CHKERRQ(ierr)

      print*,  "Solve the system"


      call KSPSetMonitor(ksp,KSPDefaultMonitor,PETSC_NULL_OBJECT,&
           PETSC_NULL_FUNCTION,ierr)
      call KSPSolve(ksp, bvec, x, ierr)
      call KSPGetConvergedReason(ksp, reason, ierr)
      call KSPGetIterationNumber(ksp, iterations, ierr)

      print*, reason
      print*, iterations

      if(reason<0) then
         if (reason==-3) then
            print*, "failed to solve matrix in given iterations"
            stop
         else
            print*, "failed to solve matrix : diverged"
            stop
         end if
      end if

      call VecGetValues(x,n,(/ (i, i=0,n_vels*n_layers-1)/),&
           temp, ierr)
      do i=1,n_layers
         u(:,i)=temp((i-1)*n_vels+1:i*n_vels)
      end do

      print*, "clean up"

      call VecDestroy(x, ierr)
      call VecDestroy(bvec, ierr)
      call MatDestroy(A, ierr)
      call KSPDestroy(ksp, ierr)

    end subroutine get_u

  subroutine smoother(xin,param)
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscviewer.h"

    real, dimension(n_vels,n_layers), intent(inout) :: xin
    real, intent(in) :: param
    real, dimension(n_vels*n_layers) :: temp
    integer, dimension(3) :: loc
    integer :: i,j,k,kk


    KSPType :: ksp_type=KSPGMRES
    PCType :: pc_type=PCSOR
    PetscErrorCode :: ierr 
    Vec :: x, bvec
    Mat :: A
    PetscTruth:: pstrue
    PetscInt :: n,np(n_vels*n_layers), iterations
    KSP :: ksp ! Krylov subspace context
    PC  :: pc  ! Preconditioner context
    KSPConvergedReason :: reason

    n=n_vels*n_layers



    print*, "Creating vector x"
    call VecCreateSeq(PETSC_COMM_WORLD, n, x,ierr)
    CHKERRQ(ierr)
    print*, "Creating vector b" 
    call VecCreateSeq(PETSC_COMM_WORLD, n, bvec,ierr)

!    print*, "setting sizes"
!    call VecSetSizes(bvec, n, n, ierr)
!    call VecSetSizes(x, n, n, ierr)
    
    ! Zero the answer.
    print*, "Zeroing answer"
      call VecZeroEntries(x, ierr)
      do i=1,n_layers
         temp((i-1)*n_vels+1:i*n_vels)=xin(:,i)
      end do
      call VecSetValues(x, n, (/ (i,i=0,n_vels*n_layers-1)/),&
           temp, INSERT_VALUES, ierr)
      print*, "Setting values in vecs"
      do i=1,n_layers
         temp((i-1)*n_vels+1:i*n_vels)=xin(:,i)
      end do
      call VecSetValues(bvec, n, (/ (i,i=0,n_vels*n_layers-1)/),&
           temp, INSERT_VALUES, ierr)
      call VecAssemblyBegin(bvec, ierr)
      call VecAssemblyEnd(bvec, ierr)

      print*, "construct the matrix"

      np=3*n_layers

      call MatCreateSeqAIJ(MPI_COMM_SELF,n,&
           n,0,np, A, ierr)
      call MatZeroEntries(A, ierr)
      do i=0,n_vels-1
         do j=0,n_layers-1
               loc(1:3)=i+(/-1,0,1/)
               where (loc .lt. 0) loc =loc +n_vels
               where (loc .ge. n_vels) loc = loc -n_vels
               call MatSetValues(A, 1, j*n_vels+i,3, &
                    loc+j*n_vels, &
                    (/-param/(dx**2),1.0+2.0*param/(dx**2),-param/(dx**2)/),&
                    INSERT_VALUES, ierr)
         end do
      end do

      print*, "assemble the matrix"

      call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

      print*, "Set up the solver context"

      call KSPCreate(MPI_COMM_SELF, ksp, ierr); CHKERRQ(ierr)
      call KSPSetInitialGuessNonZero(ksp,PETSC_TRUE,ierr)
      call KSPSetType(ksp, ksp_type, ierr); CHKERRQ(ierr)
      call KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN, ierr); 
      CHKERRQ(ierr)
      call KSPGetPC(ksp, pc, ierr); CHKERRQ(ierr)
      call PCSetType(pc, pc_type, ierr); CHKERRQ(ierr)
      call KSPSetTolerances(ksp,1.0e-16, 1.0e-15,     &
           PETSC_DEFAULT_DOUBLE_PRECISION,30000,ierr); CHKERRQ(ierr)

      print*,  "Solve the system"


      call KSPSetMonitor(ksp,KSPDefaultMonitor,PETSC_NULL_OBJECT,&
           PETSC_NULL_FUNCTION,ierr)
      call KSPSolve(ksp, bvec, x, ierr)
      call KSPGetConvergedReason(ksp, reason, ierr)
      call KSPGetIterationNumber(ksp, iterations, ierr)

      print*, reason
      print*, iterations

      if(reason<0) then
         if (reason==-3) then
            print*, "failed to solve matrix in given iterations"
            stop
         else
            print*, "failed to solve matrix : diverged"
            stop
         end if
      end if

      call VecGetValues(x,n,(/ (i, i=0,n_vels*n_layers-1)/),&
           temp, ierr)
      do i=1,n_layers
         xin(:,i)=temp((i-1)*n_vels+1:i*n_vels)
      end do

      print*, "clean up"

      call VecDestroy(x, ierr)
      call VecDestroy(bvec, ierr)
      call MatDestroy(A, ierr)
      call KSPDestroy(ksp, ierr)


    end subroutine smoother

  subroutine get_D_flux(Din,uin,f)
    real, intent(in) :: Din(n_dens,n_layers), uin(n_vels,n_layers)
    real, intent(out) :: f(n_dens,n_layers)
    integer :: i,j
    real, dimension(0:n_vels+1,n_layers) :: phi, diff1, diff2
    real, dimension(0:n_dens+2,n_layers)::D
    real, dimension(0:n_dens+2)::h
    real, dimension(0:n_vels+1,n_layers)::u
    real, dimension(n_dens,n_layers):: Outing

     ! enforce periodicity

     D(0,:)=Din(n_dens,:);D(1:n_dens,:)=Din;D(n_dens+1:n_dens+2,:)=Din(1:2,:)  
     u(0,:)=uin(n_vels,:);u(1:n_vels,:)=uin;u(n_vels+1,:)=uin(1,:)  

        h(:)=sum(D,2)
    ! Uses QUICK upwinding method for (uD)_x term

     phi(:,:)=0

    forall(i=1:n_vels)
       where (u(i,:)-c .gt. 0) 
          phi(i,:)=(u(i,:)-c)*(0.5*(D(i+1,:)+D(i,:))&
               -0.125*(D(i+1,:)-2.0*D(i,:)+D(i-1,:))&
               )
       elsewhere (u(i,:)-c .lt. 0)
          phi(i,:)=(u(i,:)-c)*(0.5*(D(i+1,:)+D(i,:))&
               -0.125*(D(i+2,:)-2.0*D(i+1,:)+D(i,:))&
               )
       end where
end forall

forall (i=1:n_dens)
  diff1(i,:)=(D(i+1,:)-2.0*D(i,:)+D(i-1,:))/(dx**2)
end forall
diff1(0,:)=diff1(n_dens,:)
diff1(n_dens+1,:)=diff1(1,:)

forall (i=1:n_dens)
  diff2(i,:)=(Diff1(i+1,:)-2.0*Diff1(i,:)+Diff1(i-1,:))/(dx**2)
end forall
diff2(0,:)=diff2(n_dens,:)
diff2(n_dens+1,:)=diff2(1,:)

       ! hack for top level
if (.false.) then
do j=2,n_layers
   do i=1,n_vels
       if ((u(i,j)-u(i,1)).gt.0) then
        phi(i,1)=phi(i,1)+(u(i,j)-u(i,1))*(0.5*(D(i+1,j)+D(i,j))&
               -0.125*(D(i+1,j)-2.0*D(i,j)+D(i-1,j))&
               )
       elseif ((u(i,j)-u(i,1)) .lt. 0) then
          phi(i,1)=phi(i,1)+(u(i,j)-u(i,1))*(0.5*(D(i+1,j)+D(i,j))&
               -0.125*(D(i+2,j)-2.0*D(i+1,j)+D(i,j))&
               )
       end if

    end do

 end do

end if

    phi(0,:)=phi(n_vels,:)

 
    forall(i=1:n_dens)
    f(i,1)=-(phi(i,1)-phi(i-1,1))/dx+sum(kappa*diff1(i,:))-kappa2*diff2(i,1)
       f(i,2:n_layers)=-(phi(i,2:n_layers)-phi(i-1,2:n_layers))/dx+0.0*kappa*diff1(i,2:n_layers)-kappa2*diff2(i,2:n_layers)
    end forall


    if (t==20) then

    forall(i=1:n_dens)
       outing(i,:)=D(i+1,:)*(u(i+1,:)-u(i,:))/(0.5*dx)
    end forall

    open(55,file='init_nDdu.dat')
    write(55,*) outing
    close(55)
    
     forall(i=1:n_dens)
       outing(i,:)=u(i,:)*(D(i+1,:)-D(i,:))/(0.5*dx)
    end forall

    open(55,file='init_nudD.dat')
    write(55,*) outing
    close(55)
       
    end if

  end subroutine get_D_flux

  subroutine get_m_flux(rho,mn,uin,Din,b,flux)
    real, intent(in) :: rho(n_layers)
    real, intent(in) :: Din(n_dens,n_layers), uin(n_vels,n_layers)
    real, intent(in) :: mn(n_vels,n_layers), b(n_dens)
    real, intent(out) :: flux(n_vels,n_layers)
    integer :: i
    real, dimension(0:n_dens+1,n_layers) :: phi
    real, dimension(n_dens+1,n_layers) :: ld,ld_g
    real, dimension(0:n_dens+1,n_layers) :: ld_m
    real, dimension(-1:n_vels+1,n_layers)::m
    real, dimension(0:n_vels+1,n_layers)::u
    real, dimension(0:n_dens+1,n_layers)::D
    real, dimension(n_dens,n_layers) :: outing
    character(len=10) :: write_buffer

     ! enforce periodicity

     m(-1:0,:)=mn(n_vels-1:n_vels,:);m(1:n_vels,:)=mn;m(n_vels+1,:)=mn(1,:)  
     u(0,:)=uin(n_vels,:);u(1:n_vels,:)=uin;u(n_vels+1,:)=uin(1,:)
     D(0,:)=Din(n_dens,:);D(1:n_dens,:)=Din;D(n_dens+1,:)=Din(1,:)

    !Get pressure forcings
    ld=0.0
    call get_ld(ld,ld_m,ld_g,rho,uin,Din,b)
    
    ! Single (mu)_x term, thanks nice form (m= M/d remember) 

    forall(i=1:n_dens)
!          phi(i,:)= ((u(i,:))-c)*((m(i,:))&
       where ((1.0-implicit_factor)*(0.5*(u(i,:)+u(i-1,:))-c) .gt. 0)
          phi(i,:)= (0.5*(u(i,:)+u(i-1,:))-c)*(0.5*(m(i,:)+m(i-1,:))&
               -0.125*(m(i,:)-2.0*m(i-1,:)+m(i-2,:))&
               )
       elsewhere
!          phi(i,:)=(0.5*(u(i,:)+u(i-1,:))-c)*(0.5*(m(i,:)+m(i-1,:))&
          phi(i,:)= (0.5*(u(i,:)+u(i-1,:))-c)*(0.5*(m(i,:)+m(i-1,:))&
               -0.125*(m(i+1,:)-2.0*m(i,:)+m(i-1,:))&
               )
       end where
    end forall
    phi(0,:)=phi(n_dens,:)
    phi(n_dens+1,:)=phi(1,:)

     forall(i=1:n_dens)
        flux(i,:)=(ld(i+1,:)-ld(i,:))/dx
     end forall

     if (.true.) then
        call smoother(flux,alpha2)
     end if

    forall(i=1:n_dens)
       flux(i,:)=flux(i,:)-(1.0-implicit_factor)*(phi(i+1,:)-phi(i,:))/(dx)&
!            +(ld(i+1,:)-ld(i,:))/dx&
            +alpha3*(m(i+1,:)-2.0*m(i,:)+m(i-1,:))/(dx*dx)
!            +(ld_m(i+1,:)-ld_m(i-1,:))/(2.0*dx)
!            +ld_g(i,:)
    end forall

    ! "Tidal forcing term"


    if (.false.) then
       flux=flux+50.0*1000.0*(2.0*3.1415/(12.0*3600.0))*&
            sin(2.0*3.1415*t/(12.0*3600.0))&
            *(0.5*(1/D(2:n_dens+1,:)+1/D(1:n_dens,:)))
    end if
    ! All done


    if (t==20) then

       write(write_buffer,'(i4.4)') out
    open(88,file='init_ndld')
    write(88,*) (ld(2:n_dens,:)-ld(1:n_dens-1,:))/dx
    close(88)
    open(88,file='init_np')
    do i=1,n_layers
       ld_m(1:n_dens+1,i)=ld(:,i)-sum(ld(:,i))/(n_dens+1)
    end do
    write(88,*) ld_m(1:n_dens+1,:)
    close(88)
     open(88,file='init_ndlnh')
    write(88,*) (ld(2:n_dens,:)-ld(1:n_dens-1,:))/dx-ld_g(1:n_dens-1,:)
    close(88)
    open(88,file='init_ndlh')
    write(88,*) ld_g(:,:)
    close(88)
    forall (i=1:n_dens)
       outing(i,:)=u(i,:)*(m(i+1,:)-m(i-1,:))/(2.0*dx)
    end forall
    open(88,file='init_nudm.dat')
    write(88,*) outing
    close(88)
    forall (i=1:n_dens)
       outing(i,:)=m(i,:)*(u(i+1,:)-u(i-1,:))/(2.0*dx)
    end forall
    open(88,file='init_nmdu.dat')
    write(88,*) outing
    close(88)


end if

  end subroutine get_m_flux

  subroutine get_ld(ld,ld_m,ld_g,rho,uin,Din,bin)
    real, dimension(n_vels+1,n_layers), intent(out) :: ld,ld_g
    real, dimension(0:n_vels+1,n_layers), intent(out) :: ld_m
    real, dimension(n_layers), intent(in) :: rho
    real, dimension(n_vels,n_layers), intent(in) :: uin
    real, dimension(n_dens,n_layers), intent(in) :: Din
    real, dimension(n_dens), intent(in) :: bin
    integer :: i,j,k, nl_min
    real, dimension(n_layers) :: rho_local
    real, dimension(0:n_dens+1,n_layers) :: A,B_D,B_DD
   real, dimension(n_dens,n_layers) ::  p2
    real, dimension(0:n_vels+1,n_layers) :: B_u,h_x,h_xx
    real, dimension(0:n_vels+1,n_layers) ::p1
     !local Periodic versions
     real, dimension(0:n_dens+1):: b
     real, dimension(0:n_dens+2,n_layers)::D
     real, dimension(0:n_vels+1,n_layers)::u,h
     real, dimension(-1:n_vels+1,n_layers):: C_u
     real :: rho_bt
 
     if (RIGID_LID) then
        nl_min=2
     else 
        nl_min=1
     end if

     b(0)=bin(n_dens);b(1:n_dens)=bin;b(n_dens+1)=bin(1)
     D(0,:)=Din(n_dens,:);D(1:n_dens,:)=Din;D(n_dens+1:n_dens+2,:)=Din(1:2,:)  
     u(0,:)=uin(n_vels,:);u(1:n_vels,:)=uin;u(n_vels+1,:)=uin(1,:)  
    
! Aim to produce entire forcing term, (dl/dD) on D points. Write in terms of D,
! A= u_x and B= u(b-sum D)_x +sum (Du)_x. Make these in right positions, A on D
! points, B_D the B terms on D points and B_u the B terms on u points.

   forall (i=1:n_dens)
       where ((u(i,:)) .gt. 0)
          C_u(i,:)= (u(i,:))*(0.5*(D(i+1,:)+D(i,:))&
               -0.125*(D(i+1,:)-2.0*D(i,:)+D(i-1,:))&
               )
       elsewhere
          c_u(i,:)=(u(i,:))*(0.5*(D(i+1,:)+D(i,:))&
               -0.125*(D(i+2,:)-2.0*D(i+1,:)+D(i,:))&
               )
       end where
      
   end forall
   c_u(n_dens+1,:)=C_u(1,:)
   c_u(-1:0,:)=c_u(n_dens-1:n_dens,:)
    forall (i=1:n_dens,j=1:n_layers)
       A(i,j)=(u(i,j)-u(i-1,j))/dx

!       B_D(i,j)=0.5*(u(i,j)+u(i-1,j))*sum(b(i+1)-b(i-1)&
!            -D(i+1,j+1:n_layers)+D(i-1,j+1:n_layers))/(2.0*dx)&
!            +sum((c_u(i,j+1:n_layers)-c_u(i-1,j+1:n_layers)))/dx
!       B_DD(i,j)=sum((c_u(i,j+1:n_layers)-c_u(i-1,j+1:n_layers)))/dx
       B_D(i,j)=sum(D(i,j+1:n_layers)*&
            (u(i,j+1:n_layers)-u(i-1,j+1:n_layers))/dx)
    end forall
    forall (i=0:n_vels,j=1:n_layers)
       B_u(i,j)=u(i,j)*(b(i+1)-b(i)&
            -sum(D(i+1,j+1:n_layers))+sum(D(i,j+1:n_layers)))/dx&
            +sum(u(i,j+1:n_layers)&
            *(D(i+1,j+1:n_layers)-D(i,j+1:n_layers))/dx)
    end forall
    b_u(N_dens+1,:)=b_u(1,:)


    A(0,:)=A(n_dens,:); A(n_dens+1,:)=A(1,:)
    B_D(0,:)=B_D(n_dens,:); B_D(n_dens+1,:)=B_D(1,:)
    B_DD(0,:)=B_DD(n_dens,:); B_DD(n_dens+1,:)=B_DD(1,:)

!    print*, 'A=', sum(abs(A))
!    print*, 'B_u=', sum(abs(B_u))
!    print*, 'B_d=', sum(abs(B_d))
!    print*, 'u=', sum(abs(u))

 

! Initialize
    ld=0.0
    ld_m=0.0
    ld_g=0.0


! Only things on u points first. Only u^2/2 + B_u^2/2
    forall (i=0:n_vels,j=1:n_layers)
       p1(i,j)=0.5*((1.0-2.0*implicit_factor)*u(i,j)**2+(B_u(i,j))**2)
    end forall

if (t==20) then
   open(88,file='init_nu2')
   write(88,*) u**2
   close(88)
end if

    forall (i=1:n_dens,j=1:n_layers)
!       ld(i,j)=0.5*rho(j)*(0.5*(u(i,j)+u(i-1,j)))**2
       ld(i,j)=0.5*rho(j)*(p1(i,j)+p1(i-1,j))
!       ld_m(i,j)=rho(j)*p1(i,j)
!        ld_m(i,j)=0.5*(p1(i,j)+p1(i-1,j))
    end forall



! Now for a mixed term. (DA+B_D)*(B_u)
     forall (i=0:n_vels,j=1:n_layers)
        p1(i,j)=B_u(i,j)
     end forall

if (t==20) then
   open(88,file='init_nDAB')
   write(88,*) 0.5*(D(1:n_dens,:)*A(1:n_dens,:)+B_D(1:n_dens,:)&
        +0.5*(B_u(1:n_dens,:)+B_u(0:n_dens-1,:)))**2
   close(88)
end if

     forall (i=1:n_dens,j=1:n_layers)
       ld(i,j)=ld(i,j)+0.5*rho(j)*(D(i,j)*A(i,j)+B_D(i,j))&
            *(p1(i,j)+p1(i-1,j))
!       ld_m(i,j)=0.5*rho(j)*(D(i,j)*A(i,j)+B_D(i,j))&
!            *(p1(i,j)+p1(i-1,j))
    end forall

! get the rest of the D point stuff out of the way.

      forall (j=1:n_layers)
         rho_local(j+1:n_layers)=rho(j)
         rho_local(1:j)=rho(1:j)
         if (RIGID_LID) then
            rho_bt=rho(j)
            else rho_bt=rho(j)-rho(1)
         forall (i=1:n_dens)
            h(i,j)=g*rho(j)*(b(i))&
                 -g*sum(rho_bt*D(i,nl_min:n_layers))&
                 -g*sum((rho(nl_min:j-1)-rho(j))*D(i,nl_min:j-1))
!                 +0.5*rho(j)*(D(i,j)*A(i,j)+B_D(i,j))**2
             p2(i,j)=0.5*rho(j)*(D(i,j)*A(i,j)+B_D(i,j))**2 
         end forall
      end forall
      h(n_dens+1,:)=h(1,:)
      
      forall (i=1:n_dens,j=1:n_layers)
!         ld(i,j)=ld(i,j)+p1(i,j)
         ld(i,j)=ld(i,j)+h(i,j)+p2(i,j)
         ld_g(i,j)=(h(i+1,j)-h(i,j))/(dx)&
              +0.5*rho(j)*(u(i+1,j)**2-u(i-1,j)**2)/(2.0*dx)
!         ld_m(i,j)=0.5*rho(j)*(D(i,j)*A(i,j)+B_D(i,j))**2
      end forall

! Terms from implicit formulation :
! 1/D(D^3A/3+D^2B/2)_x+(DA/2+B)uh_x

  
 forall (i=1:n_dens,j=1:n_layers)
    p2(i,j)=(rho(j)/D(i,j))*(&
         (D(i+1,j)**3*A(i+1,j)-D(i-1,j)**3*A(i-1,j))/(6.0*dx)&
         +(D(i+1,j)**2*B_d(i+1,j)-D(i-1,j)**2*B_d(i-1,j))/(4.0*dx)&
         )&
         +rho(j)*D(i,j)*(B_u(i,j)-B_u(i-1,j))/(2.0*dx)
!         +rho(j)*D(i,j)*(((u(i,j)-u(i-1,j))/dx)*&
!         (b(i+1)-b(i-1)-sum(D(i+1,j+1:n_layers))+sum(D(i-1,j+1:n_layers)))&
!              /(4.0*dx)&
!         +sum(((u(i,j+1:n_layers)-u(i-1,j+1:n_layers))/dx)*&
!         (D(i+1,j+1:n_layers)-D(i-1,j+1:n_layers)))&
!              /(4.0*dx)&
!         +0.5*(u(i,j)+u(i-1,j))*(b(i+1)-2.0*b(i)+b(i-1)&
!         +sum(-D(i+1,j+1:n_layers)+2.0*D(i,j+1:n_layers)&
!         -D(i-1,j+1:n_layers)))/(2.0*(dx**2))&
!         +0.5*sum((u(i,j+1:n_layers)+u(i-1,j+1:n_layers))*(&
!         D(i+1,j+1:n_layers)-2.0*D(i,j+1:n_layers)&
!         +D(i-1,j+1:n_layers)))/(2.0*(dx**2))&
!         )

    p1(i,j)=rho(j)*B_u(i,j)*(D(i+1,j)-D(i,j))/(dx)
 end forall
 p1(0,:)=p1(n_dens,:)

     forall (i=1:n_dens,j=1:n_layers)
        ld(i,j)=ld(i,j)+implicit_factor*(p2(i,j)&
             +0.5*(p1(i,j)+p1(i-1,j))&
             )
      end forall
      forall (i=0:n_dens,j=1:n_layers)
         p1(i,j)=rho(j)*u(i,j)*(b(i+1)-b(i)&
              -sum(D(i+1,j+1:n_layers))+sum(D(i,j+1:n_layers)))/dx
      end forall

      forall (i=1:n_dens,j=1:n_layers)
         ld(i,j)=ld(i,j)+implicit_factor*(&
             (D(i,j)*A(i,j)/2.0+B_d(i,j))*0.5*(p1(i,j)+p1(i-1,j))&
             +0.5*(B_u(i,j)*p1(i,j)+B_u(i-1,j)*p1(i-1,j)))
      end forall


! Only Cu terms remain. This may get messy. Div(Cu) first.Splits into two
! like the previous term. First  (C_D.u_x+ D.(B_u u)_x) 
! C_d is D(DA/2+B_d)

      forall (i=1:n_dens,j=1:n_layers)
         p2(i,j)=sum(rho(1:j-1)*&
              (0.5*D(i,1:j-1)*D(i,1:j-1)*A(i,1:j-1)&
              +D(i,1:j-1)*B_D(i,1:j-1))*A(i,1:j-1)&
          +D(i,1:j-1)*(B_u(i,1:j-1)*u(i,1:j-1)-B_u(i-1,1:j-1)*u(i-1,1:j-1))/dx&
                   )
      end forall
      
      forall (i=1:n_dens,j=1:n_layers)
!        ld(i,j)=ld(i,j)+p2(i,j)
!        ld_m(i,j)=sum(rho(1:j-1)*&
!            D(i,1:j-1)*(B_u(i,1:j-1)*u(i,1:j-1)-B_u(i-1,1:j-1)*u(i-1,1:j-1))/dx&
!                   )
      end forall
      
! Now for u.(C_D)_x + u.B_u.D_x. All on u points
     forall (i=0:n_dens,j=1:n_layers)
         p1(i,j)=sum(u(i,1:j-1)*rho(1:j-1)*(&
              ((0.5*D(i+1,1:j-1)**2*A(i+1,1:j-1)+D(i+1,1:j-1)*(B_D(i+1,1:j-1)))&
              -(0.5*D(i,1:j-1)**2*A(i,1:j-1)+D(i,1:j-1)*B_D(i,1:j-1)))/dx&
             +B_u(i,1:j-1)*(D(i+1,1:j-1)-D(i,1:j-1))/dx&
              ))
      end forall

! Average onto D points

       forall (i=1:n_dens,j=1:n_layers)
!        ld(i,j)=ld(i,j)+0.5*(p1(i,j)+p1(i-1,j))
      end forall

! Try a mixed approach?

     forall (i=1:n_dens,j=1:n_layers)
!         p2(i,j)=sum(0.5*B_D(i,1:j-1)*(u(i,1:j-1)+u(i-1,1:j-1))*rho(1:j-1)*(&
!              (D(i+1,1:j-1)-D(i-1,1:j-1))/(2.0*dx)&
!              ))
      end forall

    forall (i=1:n_dens,j=1:n_layers)
 !       ld(i,j)=ld(i,j)+p2(i,j)
      end forall

    forall (i=0:n_dens,j=1:n_layers)
!         p1(i,j)=sum(u(i,1:j-1)*rho(1:j-1)*(&
!              ((0.5*D(i+1,1:j-1)**2*A(i+1,1:j-1)+D(i+1,1:j-1)*(B_D(i+1,1:j-1)))&
!              -(0.5*D(i,1:j-1)**2*A(i,1:j-1)+D(i,1:j-1)*B_D(i,1:j-1)))/dx&
!              ))
      end forall


      !Now for u.gradC. Now only have points on u, unfortunately.

      ! First u.(C_D)_x +u.B_u.D_x
      forall (i=0:n_dens,j=1:n_layers)
         p1(i,j)=u(i,j)*sum(rho(1:j-1)*&
              (0.5*D(i+1,1:j-1)**2*A(i+1,1:j-1)&
              -0.5*D(i,1:j-1)**2*A(i,1:j-1)&
              +D(i+1,1:j-1)*B_D(i+1,1:j-1)&
              -D(i,1:j-1)*B_D(i,1:j-1))/dx)&
              +u(i,j)*sum(rho(1:j-1)*B_u(i,1:j-1)&
              *(D(i+1,1:j-1)-D(i,1:j-1))/dx)
      end forall

      forall (i=1:n_dens,j=1:n_layers)
!         ld(i,j)=ld(i,j)-0.5*(p1(i,j)+p1(i-1,j))
      end forall

      ! Finally u.D.(B_u)_x. Must average to u or D points. D doesn't work.
      ! We try u.
      p1=0.0
      do i=1,n_dens
         do j=1,n_layers
            do k=1,j-1
               p1(i,j)=p1(i,j)+u(i,j)&
                    *0.5*(rho(k)*(D(i,k)+D(i+1,k))&
                    *(B_u(i+1,k)-B_u(i-1,k)))/(2.0*dx)
            end do
         end do
      end do
      p1(0,:)=p1(n_dens,:)

      forall (i=1:n_dens,j=1:n_layers)
!         ld_m(i,j)=ld_m(i,j)-p1(i,j)
!         ld(i,j)=ld(i,j)-0.5*(p1(i,j)+p1(i-1,j))
!         ld_m(i,j)=-0.5*(p1(i,j)+p1(i-1,j))
      end forall

! Look at using (u_j-u_i).grad C_j + C_j div u_j
      p1=0.0
      c_u=0.0
      do i=1,n_dens
         do j=1,n_layers
            do k=1,j-1
               C_u(i,j)=C_u(i,j)+(0.5*(D(i+1,k)**2*A(i+1,k)-D(i,k)**2*A(i,k))&
                    +D(i+1,k)*B_D(i+1,k)-D(i,k)*B_D(i,k))/dx
            end do
         end do
      end do
      C_u(-1:0,:)=C_u(n_dens-1:n_dens,:)
      C_u(n_dens+1,:)=C_u(1,:)

      do i=1,n_dens
         do j=1,n_layers
            do k=1,j-1

               if (0.5*(u(i,k)+u(i-1,k)&
                    -(1.0-implicit_factor)*(u(i,j)+u(i-1,j))) .gt. 0) then
                  p1(i,j)=p1(i,j)+rho(k)*0.5*(u(i,k)+u(i-1,k)&
                       -(1.0-implicit_factor)*(u(i,j)+u(i-1,j)))*&
                  (0.5*(C_u(i,k)+C_u(i-1,k))&
!               -0.125*(C_u(i,k)-2.0*C_u(i-1,k)+C_u(i-2,k))&
               )
               else
                 p1(i,j)=p1(i,j)+rho(k)*0.5*(u(i,k)+u(i-1,k)&
                      -(1.0-implicit_factor)*(u(i,j)+u(i-1,j)))*&
                  (0.5*(C_u(i,k)+C_u(i-1,k))&
!              -0.125*(C_u(i+1,k)-2.0*C_u(i,k)+C_u(i-1,k))&
               )
              end if
!               p1(i,j)=p1(i,j)+rho(k)*(u(i,k)-u(i,j))*&
!                    (D(i+1,k)**2*A(i+1,k)-D(i,k)**2*A(i,k)&
!                    +D(i+1,k)*B_D(i+1,k)-D(i,k)*B_D(i,k)&
!                    +B_u(i,k)*(D(i+1,k)-D(i,k)))/(dx)&
!                    +rho(k)*(u(i,k)-u(i,j))*&
!                   0.5*(D(i+1,k)+D(i,k))*(B_u(i+1,k)-B_u(i-1,k))/(2.0*dx)
            end do
         end do
      end do
      p1(0,:)=p1(n_dens,:)

      p2=0.0
      c_u=0.0
      do i=1,n_dens
         do j=1,n_layers
            do k=1,j-1
               C_u(i,j)=C_u(i,j)&
                    +0.5*(D(i,k)+D(i+1,k))*(B_u(i+1,k)-B_u(i-1,k))/(2.0*dx)
            end do
         end do
      end do
      C_u(-1:0,:)=C_u(n_dens-1:n_dens,:)
      C_u(n_dens+1,:)=C_u(1,:)
   
      do i=1,n_dens
         do j=1,n_layers
            do k=nl_min,j-1

               if (0.5*(u(i,k)+u(i-1,k)&
                    -(1.0-implicit_factor)*(u(i,j)+u(i-1,j))) .gt. 0) then
                  p2(i,j)=p2(i,j)+rho(k)*0.5*(u(i,k)+u(i-1,k)&
                       -(1.0-implicit_factor)*(u(i,j)+u(i-1,j)))*&
                  (0.5*(C_u(i,k)+C_u(i-1,k))&
!               -0.125*(C_u(i,k)-2.0*C_u(i-1,k)+C_u(i-2,k))&
               )
               else
                 p2(i,j)=p2(i,j)+rho(k)*0.5*(u(i,k)+u(i-1,k)&
                      -(1.0-implicit_factor)*(u(i,j)+u(i-1,j)))*&
                  (0.5*(C_u(i,k)+C_u(i-1,k))&
!              -0.125*(C_u(i+1,k)-2.0*C_u(i,k)+C_u(i-1,k))&
               )
              end if
!               p2(i,j)=p2(i,j)+rho(k)*0.5*(u(i,k)+u(i-1,k)-u(i,j)-u(i-1,j))*&
!                    (D(i,k)*(B_u(i,k)-B_u(i-1,k)))/(dx)&
!                   + D(i,k)*(0.5*(u(i,k)+u(i-1,k))*(b(i+1)-2.0*b(i)+b(i-1))&
!                    +sum(0.5*(u(i,k+1:N_layers)+u(i-1,k+1:n_layers)&
!                    -u(i,k)-u(i-1,k))*(D(i+1,k+1:n_layers)&
!                    -2.0*D(i,k+1:n_layers)+D(i-1,k+1:n_layers))))/(dx**2)
            end do
         end do
      end do

      forall (i=1:n_dens,j=nl_min:n_layers)
         ld(i,j)=ld(i,j)+p2(i,j)+0.5*(p1(i,j)+p1(i-1,j))
         ld_m(i,j)=p2(i,j)+0.5*(p1(i,j)+p1(i-1,j))
      end forall


      p2=0.0
      do i=1,n_dens
         do j=1,n_layers
            do k=nl_min,j-1
               p2(i,j)=p2(i,j)+rho(k)*A(i,k)*&
                    (D(i,k)**2*A(i,k)+D(i,k)*B_D(i,k))
            end do
         end do
      end do

      forall (i=1:n_dens,j=nl_min:n_layers)
        ld(i,j)=ld(i,j)+p2(i,j)
        ld_m(i,j)=ld_m(i,j)+p2(i,j)
!        ld_m(i,j)=p2(i,j)
      end forall

      if (RIGID_LID) then
         
         forall (i=1:n_vels,j=nl_min:n_layers)
            p1(i,j)=0.5*((1.0-2.0*implicit_factor)*u(i,1)**2&
                 -u(i,j)*u(i,1)+(u(i,j)-u(i,1))*&
                 (D(i+1,1)**3*A(i+1,1)-D(i,1)**3*A(i,1))&
                 /(dx*0.5*(D(i,1)+D(i+1,1)))
         end forall
         forall (i=0:n_vels+1,j=nl_min:n_layers)
            p2(i,j)=0.5*D(i,1)*D(i,1)*A(i,1)*A(i,1)
         end forall

      end if

      ! Have now got forcing. Hopefully.

      ld(n_dens+1,:)=ld(1,:)
      ld_g(n_dens+1,:)=ld_g(1,:)
      ld_m(n_dens+1,:)=ld_m(1,:)
      ld_m(0,:)=ld_m(n_dens,:)

      if (t==20) then
         open(88,file='init_nnl')
         write(88,*) ld_m
         close(88)
      end if

   end subroutine get_ld

   subroutine get_matrix(rho,Din,bin,mat)
     real, dimension(n_layers), intent(in) :: rho
     real, dimension(n_dens), intent(in) :: bin
     real, dimension(n_dens,n_layers), intent(in) :: Din
     real, dimension(n_vels*n_layers,n_layers,-1:1) :: mat
     integer :: i,j,k,kk


     !local Periodic versions
     real, dimension(0:n_dens+2):: b
     real, dimension(0:n_dens+2,n_layers)::D
     real, dimension(0:n_vels+1,n_layers)::h_x,h_xx

     ! Subrountine makes the (block triangular) matrix to be inverted in
     ! the momentum velocity solve

     b(0)=bin(n_dens);b(1:n_dens)=bin;b(n_dens+1:n_dens+2)=bin(1:2)
     D(0,:)=Din(n_dens,:);D(1:n_dens,:)=Din;D(n_dens+1:n_dens+2,:)=Din(1:2,:)  
 
     mat=0.0

     if( .true.) then

     ! first term is rho*u
     forall (i=1:n_vels,j=1:n_layers)
        mat((j-1)*n_vels+i,j,0)=rho(j)
     end forall

     end if

     ! next -1/D*(D^3A/3)_x. Turns into -D^2u_xx/3 -D.D_x.u_x

     if (.true.) then

     forall (i=1:n_vels,j=1:n_layers)   
        mat((j-1)*n_vels+i,j,0)=mat((j-1)*n_vels+i,j,0)+&
             2.0*rho(j)*((0.5*(D(i+1,j)**2+D(i,j)**2)))/(3.0*dx**2)
        mat((j-1)*n_vels+i,j,1)=&
             -rho(j)*((0.5*(D(i+1,j)**2+D(i,j)**2)))/(3.0*dx**2)&
             -rho(j)*0.5*(D(i+1,j)+D(i,j))*(D(i+1,j)-D(i,j))/(2.0*dx**2)
        mat((j-1)*n_vels+i,j,-1)=&
             -rho(j)*((0.5*(D(i+1,j)**2+D(i,j)**2)))/(3.0*dx**2)&
             +rho(j)*0.5*(D(i+1,j)+D(i,j))*(D(i+1,j)-D(i,j))/(2.0*dx**2)
     end forall


     end if

     ! Worth making h_x before we start

     forall (i=1:n_vels,j=1:n_layers)
        h_x(i,j)=(b(i+1)-b(i)&
             -(sum(D(i+1,j+1:n_layers))-sum(D(i,j+1:n_layers))))/dx
     end forall
     h_x(0,:)=h_x(n_vels,:);h_x(n_vels+1,:)=h_x(1,:)

   ! Also make h_xx

    forall (i=1:n_vels,j=1:n_layers)
       h_xx(i,j)=(h_x(i+1,j)-h_x(i-1,j))/(2.0*dx)
    end forall
    h_xx(0,:)=h_xx(n_vels,:);h_xx(n_vels+1,:)=h_xx(1,:)



     ! now -1/D*(D^2*B/2)_x. Two sets of terms :
     ! -D_x.(u.(b-sum D)_x +sum [ D.u_x+u.D_x])
     ! -D.(u_x.(b-sum D)_x +u.(b- sum D)_xx + sum [D.u_xx +2D_x.u_x+u.D_xx])/2
     ! Do the terms on u_i (un-summed) first.

    if (.true.) then

     forall (i=1:n_vels,j=1:n_layers)

        mat((j-1)*n_vels+i,j,0)=mat((j-1)*n_vels+i,j,0)&
             -rho(j)*(D(i+1,j)-D(i,j))*h_x(i,j)/dx&
            -rho(j)*0.5*(D(i+1,j)+D(i,j))*h_xx(i,j)/2.0
      

        mat(i+(j-1)*n_vels,j,1)=mat(i+(j-1)*n_vels,j,1)&
             -rho(j)*0.5*(D(i+1,j)+D(i,j))&
             *h_x(i,j)/(2.0*2.0*dx)

        mat((j-1)*n_vels+i,j,-1)=mat((j-1)*n_vels+i,j,-1)&
             +rho(j)*0.5*(D(i+1,j)+D(i,j))*&
             h_x(i,j)/(2.0*2.0*dx)
     end forall

     end if


 ! -D_x.(u.(b-sum D)_x +sum [ D.u_x+u.D_x])
 ! -D.(u_x.(b-sum D)_x +u.(b- sum D)_xx + sum [D.u_xx +2D_x.u_x+u.D_xx])/2

     if (.true.) then
      ! Now for the interface terms
     forall (i=1:n_vels,j=1:n_layers)
        forall (k=j+1:n_layers)
           mat((j-1)*n_vels+i,k,0)=mat((j-1)*n_vels+i,k,0)-&
                (D(i+1,k)-D(i,k))*(rho(j)*(D(i+1,j)-D(i,j)))&
                /(dx**2)&
                -(D(i+2,k)-D(i+1,k)-D(i,k)+D(i-1,k))/(2.0*dx**2)*&
                (rho(j)*0.5*(D(i+1,j)+D(i,j)))/2.0&
                +2.0*0.5*(D(i+1,k)+D(i,k))*(rho(j)*&
                0.5*(D(i+1,j)+D(i,j)))/(2.0*dx**2)
           mat((j-1)*n_vels+i,k,1)=mat((j-1)*n_vels+i,k,1)-&
                0.5*(D(i+1,k)+D(i,k))*(rho(j)*&
                (D(i+1,j)-D(i,j))/dx)/(2.0*dx)&
                -2.0*(D(i+1,k)-D(i,k))/dx*(rho(j)*&
                0.5*(D(i+1,j)+D(i,j)))/(2.0*2.0*dx)&
                -0.5*(D(i+1,k)+D(i,k))*(rho(j)*&
                0.5*(D(i+1,j)+D(i,j)))/(2.0*dx**2)
           mat((j-1)*n_vels+i,k,-1)=mat((j-1)*n_vels+i,k,-1)+&
                0.5*(D(i+1,k)+D(i,k))*(rho(j)*&
                (D(i+1,j)-D(i,j))/dx)/(2.0*dx)&
                +2.0*(D(i+1,k)-D(i,k))/dx*(rho(j)*&
                0.5*(D(i+1,j)+D(i,j)))/(2.0*2.0*dx)&
                -0.5*(D(i+1,k)+D(i,k))*(rho(j)*&
                0.5*(D(i+1,j)+D(i,j)))/(2.0*dx**2)
        end forall
     end forall

     end if

     ! Onto (DA/2+u.(b- sum D)_x+sum[D.u_x+u.D_x]).sum.(b-D)_x.
     ! Again non summed terms first

     if (.true.) then
        
     forall (i=1:n_vels,j=1:n_layers)
        mat((j-1)*n_vels+i,j,0)=mat((j-1)*n_vels+i,j,0)+&
             rho(j)*h_x(i,j)*h_x(i,j)
        mat((j-1)*n_vels+i,j,1)=mat((j-1)*n_vels+i,j,1)+&
             rho(j)*0.5*(D(i+1,j)+D(i,j))*h_x(i,j)/(4.0*dx)
        mat((j-1)*n_vels+i,j,-1)=mat((j-1)*n_vels+i,j,-1)-&
             rho(j)*0.5*(D(i+1,j)+D(i,j))*h_x(i,j)/(4.0*dx)
     end forall

     end if

     if (.true.) then

     ! Now the sum terms

     forall (i=1:n_vels,j=1:n_layers)
        forall (k=j+1:n_layers)
           mat((j-1)*n_vels+i,k,0)=mat((j-1)*n_vels+i,k,0)+&
                (D(i+1,k)-D(i,k))/dx*(rho(j)*h_x(i,j))
           mat((j-1)*n_vels+i,k,1)=mat((j-1)*n_vels+i,k,1)&
                +0.5*(D(i+1,k)+D(i,k))*(rho(j)*h_x(i,j))/(2.0*dx)
           mat((j-1)*n_vels+i,k,-1)=mat((j-1)*n_vels+i,k,-1)&
                -0.5*(D(i+1,k)+D(i,k))*(rho(j)*h_x(i,j))/(2.0*dx)
        end forall
     end forall

     end if

     ! Last term, -sum (rho(j)*D *(D u_x/2+u. h_x +sum [D.u_x+u.D_x]))_x
     ! Explicit terms are :
     ! -sum rho(j) [D.D_x.u_x+D^2.u_xx/2 +D.u_x. h_x + D.u. h_xx
     ! +D_x.u.h_x +D_x. sum [D.u_x+u.D_x]
     ! + D. sum [D.u_xx+2D_x.u_x+u.D_xx] ]

       ! level terms first

     if (.true.) then

     forall (i=1:n_vels,j=1:n_layers)
        forall (k=1:j-1)
           mat((j-1)*n_vels+i,k,0)=mat((j-1)*n_vels+i,k,0)+rho(k)*(&
                2.0*0.5*(D(i+1,k)**2+D(i,k)**2)/(2.0*dx**2)&
                -0.5*(D(i+1,k)+D(i,k))*h_xx(i,k)&
                -(D(i+1,k)-D(i,k))*h_x(i,k)/dx)
           mat((j-1)*n_vels+i,k,1)=mat((j-1)*n_vels+i,k,1)+rho(k)*(&
                -0.5*(D(i+1,k)+D(i,k))*&
                (D(i+1,k)-D(i,k))/(2.0*dx**2)&
                -0.5*(D(i+1,k)**2+D(i,k)**2)/(2.0*dx**2)&
                -0.5*(D(i+1,k)+D(i,k))*h_x(i,k)/(2.0*dx))
           mat((j-1)*n_vels+i,k,-1)=mat((j-1)*n_vels+i,k,-1)+rho(k)*(&
                +0.5*(D(i+1,k)+D(i,k))*&
                (D(i+1,k)-D(i,k))/(2.0*dx**2)&
                -0.5*(D(i+1,k)**2+D(i,k)**2)/(2.0*dx**2)&
                +0.5*(D(i+1,k)+D(i,k))*h_x(i,k)/(2.0*dx))
        end forall
     end forall

     end if
 
     if (.true.) then
     ! Finally terms inside the double sum
     ! -sum rho(j) [D.D_x.u_x+D^2.u_xx/2 +D.u_x. h_x + D.u. h_xx
     ! +D_x.u.h_x +D_x. sum [D.u_x+u.D_x]
     ! + D. sum [D.u_xx+2D_x.u_x+u.D_xx] ]

     do i=1,n_vels
        do j=1,n_layers
        do k=1,j-1
           do kk=k+1,n_layers
              mat((j-1)*n_vels+i,kk,0)=mat((j-1)*n_vels+i,kk,0)&
                   -(D(i+1,kk)-D(i,kk))*(rho(k)*&
                   (D(i+1,k)-D(i,k)))/dx**2&
                   -(D(i+2,kk)-D(i+1,kk)-D(i,kk)+D(i-1,kk))*(rho(k)*&
                   (0.5*(D(i+1,k)+D(i,k))))/(2.0*dx**2)&
                   +2.0*0.5*(D(i+1,kk)+D(i,kk))*(rho(k)*&
                   0.5*(D(i+1,k)+D(i,k)))/(dx**2)
              mat((j-1)*n_vels+i,kk,1)=mat((j-1)*n_vels+i,kk,1)&
                   -0.5*(D(i+1,kk)+D(i,kk))*(rho(k)*&
                   (D(i+1,k)-D(i,k)))/(2.0*dx**2)&
                   -0.5*(D(i+1,kk)+D(i,kk))*(rho(k)&
                   *0.5*(D(i+1,k)+D(i,k)))/dx**2&
                -2.0*(D(i+1,kk)-D(i,kk))*(rho(k)&
                *0.5*(D(i+1,k)+D(i,k)))/(2.0*dx**2)
              mat((j-1)*n_vels+i,kk,-1)=mat((j-1)*n_vels+i,kk,-1)&
                   +0.5*(D(i+1,kk)+D(i,kk))*(rho(k)*&
                   (D(i+1,k)-D(i,k)))/(2.0*dx**2)&
                   -0.5*(D(i+1,kk)+D(i,kk))*(rho(k)&
                   *0.5*(D(i+1,k)+D(i,k)))/dx**2&
                   +2.0*(D(i+1,kk)-D(i,kk))*(rho(k)&
                   *0.5*(D(i+1,k)+D(i,k)))/(2.0*dx**2)  
           end do
        end do
     end do
  end do

  end if




  if (.true.) then

     do i=1,n_vels
        do j=1,n_layers
              mat((j-1)*n_vels+i,j,0)=mat((j-1)*n_vels+i,j,0)+2.0*(alpha/dx)**2
             mat((j-1)*n_vels+i,j,1)=mat((j-1)*n_vels+i,j,1)-(alpha/dx)**2
              mat((j-1)*n_vels+i,j,-1)=mat((j-1)*n_vels+i,j,-1)-(alpha/dx)**2
           end do
        end do
     end if


     ! Matrix made (modulo boundary conditions)

end subroutine get_matrix

end module MLCM_1d_routines

module MLCM_1d_IO
use MLCM_1d_routines
use MLCM_1d_globals
implicit none
contains

  subroutine read_data(rho,D,b,m,u)
    real, dimension(n_layers), intent(out) :: rho
    real, dimension(n_dens), intent(out) :: b
    real, dimension(n_dens,n_layers), intent(out) :: D
    real, dimension(n_vels,n_layers), intent(out) ::m
    real, dimension(n_vels,n_layers), intent(out) ::u
    real, dimension(2*mult*(n_vels-n_extra)+1,n_layers) ::buffer
    real ::  f
    integer :: i,j, read_bit=33
    character(len=8) :: read_buffer
    real, dimension(n_layers) ::dd

!    dd=(/50.0,950.0/)
    dd=(/75.0,275.0/)    

if (read_bit==1) then
!   out=187
   t=out*t_out
   write(read_buffer,'(i4.4)') out
   open(17,file='MLCM_D'//trim(read_buffer))
   read(17,*) D
   read(17,*) b
   close(17)
   open(17, file='MLCM_m'//trim(read_buffer))
   read(17,*) m
   close(17)

   rho=(/1020.0,1025.0/)

else if (read_bit==2) then

   rho=(/1023.0,1026.0/)
   dd=(/75.0,275.0/)    
   c=1.5
!   open(17, file='MLCM_Din.dat')
   open(17, file='ham_d.dat')
   print* , 'Size of buffer=', size(buffer)
   read(17,*) buffer
   close(17)
!   open(17, file='MLCM_uin.dat')
!   read(17,*) u
!   close(17)
!   open(17, file='MLCM_min.dat')
!   read(17,*) m
!   close(17)

!   c=sqrt(2.0*(100.0*900.0)/(1000.0)*10.0*(rho(2)-rho(1))/rho(1))
!   c=sqrt(5.1*(dd(1)*dd(2))/(dd(1)+dd(2))*10.0*(rho(2)-rho(1))/rho(1))
! c=sqrt(5.1*(50.0*950.0)/(1000.0)*10.0*(rho(2)-rho(1))/rho(1))


   do j=1,n_layers
      D(1:n_extra/2,j)=buffer(1,j)
      D(n_dens-n_extra/2:n_dens,j)=buffer(1,j)    
   end do
   D(n_extra/2+1:n_dens-n_extra/2,:)=&
        buffer(1:mult*(2*(n_dens-n_extra)-1):2*mult,:) 
   print*, c
   do j=1,n_layers
      u(:,j)=c*(1.0-dd(j)/buffer(mult+1,j))
   end do
   do i=n_extra/2+1,n_vels-n_extra/2
      do j=1,n_layers
      u(i,j)=c*(1.0-dd(j)/buffer(mult*(2*(i-1-n_extra/2)+1)+1,j))
      end do
   end do

!   u=0.0

   b=sum(D(1,:))

   m=0.0
!   call utom(c,rho,buffer,b,u,m)
   call utom(c,rho,D,b,u,m)
   if( n_extra>200) then
   do j=1,n_layers
      do i=1,n_extra/2+50
         m(i,j)=m(i,j)*exp(-(n_extra/2+50.0-i)/50.0)
         m(n_dens +1-i,j)=m(n_dens+1-i,j)*exp(-(n_extra/2+50.0-i)/50.0)
      end do
!      m(n_dens-500:n_dens,j)=m(n_dens-501,j)
   end do
   end if

!      call smoother(m,40.0)

else if (read_bit==3) then

   rho=(/1000.0,1020.0/)
  
   open(17, file='MLCM_Dcut.dat')
   read(17,*) D(1:n_vels-n_extra,:)
   close(17)
   open(17, file='MLCM_ucut.dat')
   read(17,*) u(1:n_vels-n_extra,:)
   close(17)
!   open(17, file='MLCM_min.dat')
!   read(17,*) m
!   close(17)

!   c=sqrt(2.0*(100.0*900.0)/(1000.0)*10.0*(rho(2)-rho(1))/rho(1))
!   c=sqrt(5.1*(dd(1)*dd(2))/(dd(1)+dd(2))*10.0*(rho(2)-rho(1))/rho(1))
 c=sqrt(5.1*(50.0*950.0)/(1000.0)*10.0*(rho(2)-rho(1))/rho(1))


 print*, D(1,:)
 print*, D(n_vels-n_extra,:)

   do j=1,n_layers
      do i=1,n_extra
         D(n_vels-n_extra+i,j)=(D(1,j)*(i-1)&
              +D(n_vels-n_extra,j)*(n_extra-i))/(n_extra-1)
         u(n_vels-n_extra+i,j)=(u(1,j)*(i-1)&
              +u(n_vels-n_extra,j)*(n_extra-i))/(n_extra-1)
      end do    
   end do

!   u=0.0

   b=1000.0

   m=0.0
!   call utom(c,rho,buffer,b,u,m)
   call utom(c,rho,D,b,u,m)

      call smoother(m,40.0)

else if (read_bit==10) then

   rho=(/1023.0,1026.0/)
   dd=(/75.0,275.0/)    
   c=1.5
   open(17, file='ham_d.dat')
   print* , 'Size of buffer=', size(buffer)
   read(17,*) buffer
   close(17)

   D(1:(n_dens-n_extra)/2,:)=buffer(1:(n_dens-n_extra)-1:2,:)
   D((n_dens-n_extra)/2+1:n_dens-n_extra,:)=buffer((n_dens-n_extra)-1:1:-2,:)
   do i=1,n_layers
      D(n_dens-n_extra+1:n_dens,i)=buffer(1,i)
   end do

   print*, c
   do i=1,n_vels/2
      do j=1,n_layers
      u(i,j)=c*(1.0-dd(j)/buffer(2*(i-n_extra/2),j))
      u(n_vels/2+i,j)=-c*(1.0-dd(j)/buffer(2*(n_vels/2+1-i),j))
      end do
   end do

!   u=0.0

   b=sum(D(1,:))

   m=0.0
!   call utom(c,rho,buffer,b,u,m)
   call utom(c,rho,D,b,u,m)

else if (read_bit==11) then

   rho=(/1023.0,1026.0/)
   dd=(/75.0,275.0/)    
   c=1.5
   open(17, file='ham_d.dat')
   print* , 'Size of buffer=', size(buffer)
   read(17,*) buffer
   close(17)

   D(1:n_dens-n_extra,:)=buffer(1:mult*(2*(n_dens-n_extra)-1):2*mult,:)
   do i=1,n_layers
   D(n_dens-n_extra+1:n_dens,i)=buffer(1,i)
   end do


   print*, c
   do i=1,n_vels-n_extra
      do j=1,n_layers
      u(i,j)=c*(1.0-dd(j)/buffer(mult*(2*i-1)+1,j))
      end do
   end do
   do i=1,n_layers
   u(n_vels-n_extra+1:n_vels,i)=c*(1.0-dd(i)/buffer(mult+1,i))
   end do


   

!   u=0.0

   b=sum(D(1,:))

   do i=60,n_dens
      if (i<=(n_dens-60)/2+60) then
         f=50.0*(-1.0-tanh(1.0*(i-60-2*(n_dens-60)/16)/10.0))
      else
!         f=50.0*tanh((3*n_dens/4-i)/2.0)
         f=50.0*(-1.0-tanh(1.0*(14*(n_dens-60)/16-i+60)/10.0))
      end if
      b(i)=b(i)-f
      D(i,2)=D(i,2)-f
   end do


   m=0.0
!   call utom(c,rho,buffer,b,u,m)
   call utom(c,rho,D,b,u,m)


else if (read_bit==12) then

   rho=(/1023.0,1026.0/)
   dd=(/75.0,275.0/)    
   c=1.4
   open(17, file='ham_d.dat')
   print* , 'Size of buffer=', size(buffer)
   read(17,*) buffer
   close(17)

   D(1:(n_dens-n_extra)/2,:)=buffer(1:mult*((n_dens-n_extra)-1):2*mult,:)




   print*, c
   do i=0,n_vels/2-1
      do j=1,n_layers
      u(i+1,j)=c*(1.0-dd(j)/buffer(2*mult*(i-n_extra/2)+mult+1,j))
      end do
   end do

   c=1.6
   open(17, file='ham_d2.dat')
   print* , 'Size of buffer=', size(buffer)
   read(17,*) buffer
   close(17)


   D((n_dens-n_extra)/2+1:n_dens-n_extra,:)=&
        buffer(1:mult*((n_dens-n_extra)-1):2*mult,:)
   do i=1,n_layers
      D(n_dens-n_extra+1:n_dens,i)=buffer(1,i)
   end do

   do i=0,n_vels/2-1
      do j=1,n_layers
      u(n_vels/2+i+1,j)=c*(1.0-dd(j)/buffer(2*mult*(i)+mult+1,j))
      end do
   end do

!   u=0.0

   b=sum(D(1,:))

   m=0.0
!   call utom(c,rho,buffer,b,u,m)
   call utom(c,rho,D,b,u,m)

else
    do i=1,size(rho)
       rho(i)=1015.0+5.0*i
    end do

    do j=1,size(D,2)
       do i=1,size(D,1)
          D(i,j)=50.0
          if( j==3) D(i,j)=900.0
       end do
    end do

   do i=1,size(D,1)
      f=100*exp(-(dx*(i-1-size(D,1)/2)/1.0e2)**2)
!      D(i,size(D,2))=D(i,size(D,2))-f
      b(i)=sum(D(i,:))
      D(i,1)=D(i,1)+0.1*f
       D(i,2)=D(i,2)-0.1*f
       m(i,:)=0.0
    end do

    m=0.0
end if
!    m=100.00

c=0.0


  end subroutine read_data

  subroutine utom(c,rho,Din,botin,uin,mout)
    real, intent(in) :: c
    real, dimension(:), intent(in) :: rho, botin
    real, dimension(:,:), intent(in) :: Din
    real, dimension(:,:), intent(inout) :: uin, mout
    integer :: i,j,k
    real, dimension(size(mout,1),n_layers) :: m
    real, dimension(0:n_dens+2):: bot
    real, dimension(0:size(Din,1)+2,n_layers):: D! , A, F, B, u
    real, dimension(0:size(uin,1)+2,n_layers):: u
    real, dimension(n_layers*size(uin,1),n_layers,-1:1) :: mat
    real, dimension(2) ::dd
    integer :: nv, nd

    nv=size(uin,1)
    nd=size(Din,1)

!    dd=(/100.0,900.0/)
!    dd=(/50.0,950.0/)
    dd=(/75.0,275.0/)

    bot(0)=botin(n_dens);bot(1:n_dens)=botin;bot(n_dens+1:n_dens+2)=botin(1:2)  
    u(0,:)=uin(size(uin,1),:);u(1:size(uin,1),:)=uin;
   u(size(uin,1)+1:size(uin,1)+2,:)=uin(1:2,:)  
    D(0,:)=Din(size(Din,1),:);D(1:size(Din,1),:)=Din;
 D(size(Din,1)+1:size(Din,1)+2,:)=Din(1:2,:)  
    

!    call  get_matrix(rho,Din(1:size(Din,1):2,:),botin,mat)
    call  get_matrix(rho,Din,botin,mat)
    m=0.0
    mout=0.0
    do j=1,n_layers
       do i=1,size(mout,1)
          do k=1,n_layers
             mout(i,j)=mout(i,j)+mat(n_vels*(j-1)+i,k,0)*u(i,k)&
             +mat(nv*(j-1)+i,k,-1)*u(i-1,k)&
             +mat(nv*(j-1)+i,k,1)*u(i+1,k)

          end do

!!$             m(i,j)=rho(j)*(&
!!$             c*(D(2*i,j)-dd(j))&
!!$                  -c*dd(j)*( -(D(2*i+1,j)-D(2*i-1,j))**2.0/((3.0/2.0)*dx**2)&
!!$              +D(2*i,j)*(D(2*i+1,j)-2.0*D(2*i,j)+D(2*i-1,j))/((3.0/4.0)*dx**2))&
!!$                  -c*dd(j)*(&
!!$                 (sum(D(2*i+1,j+1:n_layers))-sum(D(2*i-1,j+1:n_layers)))&
!!$                  *(D(2*i+1,j)-D(2*i-1,j))/(2.0*(dx**2.0))&
!!$                  +D(2*i,j)*(sum(D(2*i+1,j+1:n_layers))&
!!$                  -2.0*sum(D(2*i,j+1:n_layers))+sum(D(2*i-1,j+1:n_layers)))&
!!$                  /((2.0/4.0)*(dx**2)))&
!!$            -c*dd(j)*(&
!!$            (D(2*i+1,j)-D(2*i-1,j))/(2.0*dx)&
!!$             +sum(D(2*i+1,j+1:n_layers)-D(2*i-1,j+1:n_layers))/dx&
!!$             )*(sum(D(2*i+1,j+1:n_layers)-D(2*i-1,j+1:n_layers)))/dx&
!!$                  )
!!$                  
!!$             m(i,j)=m(i,j)/D(2*i,j)
!!$
!!$             do k=1,j-1
!!$                m(i,j)=m(i,j)-rho(k)*c*dd(k)*(&
!!$                     (d(2*i+1,k)-2.0*d(2*i,k)+D(2*i-1,k))/(dx**2/4.0)/2.0&
!!$                     +sum((d(2*i+1,k+1:N_layers)-2.0*d(2*i,k+1:n_layers)&
!!$                     +D(2*i-1,k+1:N_Layers))/(dx**2/4.0)))
!!$             end do
       end do
    end do
    
!    mout=m
!    mout=mout-m
                
   end subroutine utom

  subroutine output_data(out,u,m,D,b)
    integer, intent(in) :: out
    real, dimension(:,:), intent(in) :: u,m,D
    real, dimension(:), intent(in) ::  b
    character(len=10) :: write_buffer

    write(write_buffer,'(i4.4)') out
    open(88,file='MLCM_D'//trim(write_buffer))
    write(88,*) D
    write(88,*) b
    close(88)
    open(88,file='MLCM_u'//trim(write_buffer))
    write(88,*) u
    close(88)
    open(88,file='MLCM_m'//trim(write_buffer))
    write(88,*) m
    close(88)
       
  end subroutine output_data
   



end module MLCM_1d_IO

program MLCM_1d
use MLCM_1d_globals
use MLCM_1d_routines
use MLCM_1d_IO
implicit none

real, allocatable, dimension(:,:) :: m, u, D 
real, allocatable, dimension(:) :: rho,b

! Must make subroutine call to initialise PETSc
call initialize_petsc()

allocate(D(n_dens,n_layers))
allocate(b(n_dens))
allocate(u(n_vels,n_layers))
allocate(m(n_vels,n_layers))
allocate(rho(n_layers))

! Get initial conditions from text-file
call read_data(rho,D,b,m,u)
call output_data(out,u,m,D,b)
call get_u(rho,m,D,b,u)
call output_data(out,u,m,D,b)
out=out+1

! Main timestepping loop (uses RK4 method)

! Note that m here is actually m/D fom the PRL paper

do while (t .le. t_max)

   
   call time_step(rho,D,m,u,b)
!   alpha3=0.999*alpha3
   t=t+dt

   if (t-out*t_out .ge. 0) then
      call output_data(out,u,m,D,b)
      out=out+1
   end if

end do

end program MLCM_1d



    
