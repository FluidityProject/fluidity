module analysis_module
contains
subroutine analysis6(A, E, S, d, ndim, nrens, nrobs, verbose)
! Computes the analysed ensemble for A using the square root formulation
! in algorithm A.2 from Evensen 2004.  I.e. the sophisitcated inversion 
! algorithm is used, which appears to be more stable with large m.
! This algorith does not use the full R but assumes (nrens-1) R = E E^T and E is 
! specified.

   use m_multa
   implicit none
   integer, intent(in) :: ndim             ! dimension of model state
   integer, intent(in) :: nrens            ! number of ensemble members
   integer, intent(in) :: nrobs            ! number of observations
   
   real, intent(inout) :: A(ndim,nrens)    ! ensemble matrix
   real, intent(in)    :: S(nrobs,nrens)   ! matrix holding HA' 
   real, intent(in)    :: d(nrobs)         ! vector holding d-HA
   real, intent(in)    :: E(nrobs,nrens)   ! matrix holding observation perturbations
   logical, intent(in) :: verbose


   real ave(ndim)       ! ensemble mean

   real S0(nrobs,nrens)

   real, allocatable :: U0(:,:),sig0(:),VT0(:,:)
   real, allocatable :: U1(:,:),sig1(:),VT1(:,:)
   real, allocatable :: U2(:,:),sig2(:),VT2(:,:)
   real, allocatable :: X0(:,:),X1(:,:),X2(:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real, allocatable, dimension(:)   :: work,isigma

   real X3(nrens,nrens)
   real, allocatable :: y1(:)
   real, allocatable :: y2(:)
   real y3(nrobs)
   real y4(nrens) 

   real sigsum,sigsum1
   integer ierr,nrsigma,i,j,lwork, nrmin
   integer iblkmax

   nrmin=min(nrobs,nrens)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subtract mean from ensemble
   ave(:)=A(:,1)
   do i=2,nrens
      ave(:)=ave(:)+A(:,i)
   enddo
   ave=(1.0/real(nrens))*ave
   do i=1,nrens
      A(:,i)=A(:,i)-ave(:)
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute SVD of S=HA'  ->  U0, sig0
   S0=S
   allocate (U0(nrobs,nrmin)  )
   allocate (sig0(nrmin))
   allocate (VT0(1,1))
   lwork=2*max(3*nrens+nrobs,5*nrens)
   allocate(work(lwork))
   sig0=0.0
!$OMP CRITICAL
   call dgesvd('S', 'N', nrobs, nrens, S0, nrobs, sig0, U0, nrobs, VT0, nrens, work, lwork, ierr)
!$OMP END CRITICAL
   deallocate(work,VT0)
   if (ierr /= 0) then
      print *,'analysis6: ierr from call dgesvd 0= ',ierr; stop
   endif

   sigsum=sum( sig0(1:nrmin) )
   sigsum1=0.0
! Significant eigenvalues.
   nrsigma=0
   do i=1,nrmin                       
      if (sigsum1/sigsum < 0.9999) then
         nrsigma=nrsigma+1
         sigsum1=sigsum1+sig0(i)
      else
         sig0(i:nrmin)=0.0
         exit
      endif
   enddo

   if (verbose) then
      write(*,'(a,i5,g13.5)') ' dominant sing. values and share ',nrsigma,sigsum1/sigsum
      write(*,'(5g11.3)')sig0
   endif

   do i=1,nrsigma
       sig0(i) = 1.0/sig0(i)
   enddo
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute X0=sig0^{*T} U0^T E      

! X0= U0^T R:q
 
   allocate(X0(nrmin,nrens))
   call dgemm('t', 'n', nrmin, nrens, nrobs, 1.0d0, U0, nrobs, E, nrobs, 0.0d0, X0, nrmin)
   !X0=matmul(transpose(U0),E)

   do j=1,nrens
   do i=1,nrmin
      X0(i,j)=sig0(i)*X0(i,j)
   enddo
   enddo
 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute singular value decomposition  of X0(nrmin,nrens)
   allocate (U1(nrmin,nrmin)  )
   allocate (sig1(nrmin))
   allocate (VT1(1,1))
   allocate (VT0(1,1))
   lwork=2*max(3*nrens+nrobs,5*nrens)
   allocate(work(lwork))
   sig1=0.0
!$OMP CRITICAL
   call dgesvd('S', 'N', nrmin, nrens, X0, nrmin, sig1, U1, nrmin, VT0, nrens, work, lwork, ierr)
!$OMP END CRITICAL
   deallocate(work,VT1)
   if (ierr /= 0) then
      print *,'analysis6: ierr from call dgesvd 1= ',ierr; stop
   endif

   do i=1,nrmin
      sig1(i)=1.0/(1.0+sig1(i)**2)
   enddo
  
   deallocate (X0,VT0)
 !print *,'sig1:',nrmin
!print '(5g11.3)',sig1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! X1 = U0 * sig0^{-1} * U1
   do j=1,nrmin
   do i=1,nrmin
      U1(i,j)=sig0(i)*U1(i,j)
   enddo
   enddo


   allocate(X1(nrobs,nrmin))
   call dgemm('n','n',nrobs,nrmin,nrmin, 1.0,U0,nrobs, U1,nrmin, 0.0,X1,nrobs)
  ! X1=matmul(U0,U1)
   deallocate (U0,sig0,U1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update mean
   allocate(y1(1:nrmin))
   call dgemv('t',nrobs,nrmin,1.0,X1,nrobs,d,1,0.0,y1 ,1)
   allocate(y2(1:nrmin))
   y2=sig1*y1  
   call dgemv('n',nrobs,nrmin,1.0,X1,nrobs,y2,1,0.0,y3 ,1)
   call dgemv('t',nrobs,nrens,1.0,S ,nrobs,y3,1,0.0,y4 ,1)
   call dgemv('n',ndim ,nrens,1.0,A ,ndim ,y4,1,1.0,ave,1)
   deallocate(y1,y2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! X2 = (I+sig1^2)^{-0.5} * X1^T * S  
   allocate(X2(nrmin,nrens))
   call dgemm('t','n',nrmin,nrens,nrobs,1.0,X1,nrobs, S,nrobs, 0.0,X2,nrmin)
  !X2=matmul(transpose(X1),S)

   do j=1,nrens
   do i=1,nrmin
      X2(i,j)=sqrt(sig1(i))*X2(i,j)
   enddo
   enddo

   deallocate (sig1,X1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!
! SVD of X2
   lwork=2*max(3*nrens+nrens,5*nrens)
   allocate (U2(nrmin,nrmin), sig2(nrmin), VT2(nrens,nrens), work(lwork))
   sig2=0.0
!$OMP CRITICAL
   call dgesvd('N', 'A', nrmin, nrens, X2, nrmin, sig2, U2, nrmin, VT2, nrens, work, lwork, ierr)
!$OMP END CRITICAL
   deallocate(work)
   if (ierr /= 0) then
      print *,'ierr from call dgesvd 2 = ',ierr
      stop
   endif

   allocate(isigma(nrmin))
   isigma=1.0
   do i=1,nrmin
      if ( sig2(i) > 1.0 ) print *,'WARNING (analysis 6): sig2 > 1',i,sig2(i)
      isigma(i)=sqrt( max(1.0-sig2(i)**2,0.0) )
   enddo

   do j=1,nrens
      X3(:,j)=VT2(j,:)
   enddo


   do j=1,nrmin
      X3(:,j)=X3(:,j)*isigma(j)
   enddo

   !print '(a)','A5: sig2: '
   !print '(5g11.3)',sig2(1:nrmin)


! Final ensemble perturbations
   iblkmax=min(ndim,200)
   call  multa(A, X3, ndim, nrens, iblkmax )

   
   do i=1,nrens
      A(:,i)=A(:,i)+ave(:)
   enddo

   deallocate (U2,sig2,VT2,isigma)
 
end subroutine analysis6

subroutine analysis2(A, D, R, S, ndim, nrens, nrobs, verbose)
! Computes the analysed ensemble for A 
   use m_multa
   implicit none
   integer, intent(in) :: ndim             ! dimension of model state
   integer, intent(in) :: nrens            ! number of ensemble members
   integer, intent(in) :: nrobs            ! number of observations
   
   real, intent(inout) :: A(ndim,nrens)   ! ensemble matrix
   real, intent(in)    :: D(nrobs,nrens)  ! vector holding observation innovations
   real, intent(in)    :: S(nrobs,nrens)  ! matrix holding HA' 
   real, intent(inout) :: R(nrobs,nrobs)  ! Error covariance matrix for observations
   logical, intent(in) :: verbose


   real, allocatable, dimension(:,:) :: X1,X2,U,X4,Reps,V
   real X3(nrobs,nrens)
   real, allocatable, dimension(:)   :: sig,work

   real sigsum,sigsum1,oneobs(1,1)
   integer ierr,nrsigma,i,j,lwork, m
   integer iblkmax
   character(len=2) tag2



   if (nrobs > 1) then
!      R=float(nrens)*R+matmul(S,transpose(S))
      call dgemm('n','t', nrobs, nrobs, nrens, &
                 1.0, S, nrobs, &
                      S, nrobs, &
        float(nrens), R, nrobs)



      allocate (U(nrobs,nrobs)  )
      allocate (V(nrobs,nrobs)  )
      allocate (sig(nrobs)  )
      lwork=2*max(3*nrobs+nrobs,5*nrobs)
      allocate(work(lwork))
      sig=0.0
!$OMP CRITICAL
      call dgesvd('A', 'A', nrobs, nrobs, R, nrobs, sig, U, nrobs, V, nrobs, work, lwork, ierr)
!$OMP END CRITICAL
      deallocate(work)
      if (ierr /= 0) then
         print *,'ierr from call dgesvd= ',ierr
         stop
      endif

!      open(10,file='sigma.dat')
!         do i=1,nrobs
!            write(10,'(i5,g12.3)')i,sig(i)
!         enddo
!      close(10)

      sigsum=sum( sig(1:nrobs) )
      sigsum1=0.0
   ! Significant eigenvalues.
      nrsigma=0
      do i=1,nrobs                 ! singular values are in descending order
         if (sigsum1/sigsum < 0.999) then
            nrsigma=nrsigma+1
            sigsum1=sigsum1+sig(i)
            sig(i) = 1.0/sig(i)
         else
            sig(i:nrobs)=0.0
            exit
         endif
      enddo

      if (verbose) then
         write(*,'(a,i5,g13.5)') ' dominant sing. values and share ',nrsigma,sigsum1/sigsum
         write(*,'(8g12.3)')1./sig(1:nrsigma), sig(nrsigma+1:nrobs)
      endif

      allocate (X1(nrobs,nrobs))
      do i=1,nrobs
      do j=1,nrobs
         X1(i,j)=sig(i)*U(j,i)
      enddo
      enddo
      deallocate(sig)

      allocate (X2(nrobs,nrens))
!     X2=matmul(X1,D)
      call dgemm('n','n',nrobs,nrens,nrobs,1.0,X1,nrobs,D ,nrobs,0.0,X2,nrobs)
      deallocate(X1) 

!     X3=matmul(V,X2)
      call dgemm('t','n',nrobs,nrens,nrobs,1.0,V ,nrobs,X2,nrobs,0.0,X3,nrobs)
      deallocate(V)
      deallocate(X2)

   else
      oneobs=matmul(S,transpose(S))+R
      print *,'oneobs: ',oneobs(1,1)
      X3=D/oneobs(1,1)
   endif

   if (2_8*ndim*nrobs < 1_8*nrens*(nrobs+ndim)) then
!    Code for few observations ( m<nN/(2n-N) )
      if (verbose) print * ,'analysis: Representer approach is used'
      allocate (Reps(ndim,nrobs))

!    Reps=matmul(A,transpose(S))
      call dgemm('n','t',ndim,nrobs,nrens,1.0,A,ndim,S,nrobs,0.0,Reps,ndim)
!      call printreps(Reps,nrobs)

!    A=A+matmul(Reps,X3)
      call dgemm('n','n',ndim,nrens,nrobs,1.0,Reps,ndim,X3,nrobs,1.0,A,ndim)
      deallocate(Reps)

      tag2(1:2)='X3'
!      open(10,file='X.uf',form='unformatted')
!         write(10)tag2,nrens,nrobs,X3,S
!      close(10)

   else
      if (verbose) print * ,'analysis: X5 appraoch is used'
      allocate(X4(nrens,nrens))
!      X4=matmul(transpose(S),X3)
      call dgemm('t','n',nrens,nrens,nrobs,1.0,S,nrobs,X3,nrobs,0.0,X4,nrens)
      do i=1,nrens
         X4(i,i)=X4(i,i)+1.0  ! Addition with Af
      enddo

      iblkmax=min(ndim,200)
      call  multa(A, X4, ndim, nrens, iblkmax )

      tag2(1:2)='X5'
!      open(10,file='X.uf',form='unformatted')
!         write(10)tag2,nrens,X4
!      close(10)

!      open(10,file='X5.dat')
!         write(10,*)'TITLE = "X5"'
!         write(10,*)'VARIABLES = "iens_f" "iens_a" "X5"'
!         write(10,'(2(a,i5),a)')' ZONE  F=BLOCK, I=',nrens,', J=',nrens,', K=1'
!         write(10,'(30I5)')((i,i=1,nrens),j=1,nrens)
!         write(10,'(30I5)')((j,i=1,nrens),j=1,nrens)
!         write(10,'(10(1x,e12.5))')((X4(i,j),i=1,nrens),j=1,nrens)
!      close(10)

!      open(10,file='X5col.dat')
!         write(10,'(a)')'These should all be ones'
!         do j=1,nrens
!            write(10,'(i5,f10.4)')j,sum(X4(:,j))  
!         enddo
!      close(10)

!      open(10,file='X5row.dat')
!         do j=1,nrens
!            write(10,'(i5,f10.4)')j,sum(X4(j,:))/float(nrens)
!          enddo
!      close(10)

      deallocate(X4)
   endif 


end subroutine analysis2

subroutine analysis(A, D, E, S, ndim, nrens, nrobs, verbose)
! Computes the analysed ensemble for A 
   use m_multa
   implicit none
   integer, intent(in) :: ndim             ! dimension of model state
   integer, intent(in) :: nrens            ! number of ensemble members
   integer, intent(in) :: nrobs            ! number of observations
   
   real, intent(inout) :: A(ndim,nrens)   ! ensemble matrix
   real, intent(in)    :: D(nrobs,nrens)  ! matrix holding innovations
   real, intent(in)    :: S(nrobs,nrens)  ! matrix holding HA' 
   real, intent(in)    :: E(nrobs,nrens)  ! matrix holding observation perturbations

   logical, intent(in) :: verbose

   real, allocatable, dimension(:,:) :: X1,X2,U,X4,Reps
   real ES(nrobs,nrens), X3(nrobs,nrens), V(nrens,nrens)
   real, allocatable, dimension(:)   :: sig,work

! Diagnostic tests
   integer, parameter :: target= 3*22+1  ! TEM(01) after u,v,d
!   integer, parameter :: target= 4*22+1  ! SAL(1) after u,v,d,t
!   integer, parameter :: target= 12  ! Uc(12) 
   real, dimension(nrobs) :: Kweights, Y1
   real, allocatable :: Y2(:)

   real sigsum,sigsum1
   integer ierr,nrsigma,i,j,lwork, nrmin
   integer iblkmax

   logical global   ! for printout purpose only
   character*9 id   ! information : nature of data + julian day

   nrmin=min(nrobs,nrens)

   ES=S+E

! Compute SVD of HA'+E  ->  U, sig
   allocate (U(nrobs,nrmin)  )
   allocate (sig(nrmin)  )
   lwork=2*max(3*nrens+nrobs,5*nrens)
   allocate(work(lwork))

   sig=0.0
!$OMP CRITICAL
   call dgesvd('S', 'N', nrobs, nrens, ES, nrobs, sig, U, nrobs, V, nrens, work, lwork, ierr)
!$OMP END CRITICAL

   deallocate(work)
   if (ierr /= 0) then
      print *,'ierr from call dgesvd= ',ierr
      stop
   endif
   do i=1,nrmin
      sig(i)=sig(i)**2
   enddo

   sigsum=sum( sig(1:nrmin) )
   sigsum1=0.0
! Significant eigenvalues.
   nrsigma=0
   do i=1,nrmin                 ! singular values are in descending order
      if (sigsum1/sigsum < 0.999) then
         nrsigma=nrsigma+1
         sigsum1=sigsum1+sig(i)
         sig(i) = 1.0/sig(i)
      else
         sig(i:nrmin)=0.0
         exit
      endif
   enddo
   if (verbose) then
      open(10,file='sigma.dat')
      write(10,'(a,i5,g13.5)') 'nb dominant sing. values and share ', &
           nrsigma,sigsum1/sigsum
      write(10,'(10f8.3)')1./sig(1:nrsigma), sig(nrsigma+1:nrmin)
      close(10)
   endif



   allocate (X1(nrmin,nrobs))
   do i=1,nrmin
   do j=1,nrobs
      X1(i,j)=sig(i)*U(j,i)
   enddo
   enddo
   deallocate(sig)

! Diagnostics tests
   if (verbose) then
      Y1=matmul(A(target,:),transpose(S))
      allocate(Y2(nrmin))
      Y2=matmul(Y1,U)
      Kweights=matmul(Y2,X1)
      deallocate(Y2)
      print *,'Forecast, TEM(01) member 1 =', A(target,1)
      print*, 'Kriging weights '
      print '(10f8.3)', Kweights
      print*, 'Innovation mem1 '
      print '(10f8.3)', D(:,1)
   endif

   allocate (X2(nrmin,nrens))
!   X2=matmul(X1,D)
   call dgemm('n','n',nrmin,nrens,nrobs,1.0,X1,nrmin, &
              D,nrobs,0.0,X2,nrmin)
   deallocate(X1) 

!   X3=matmul(U,X2)
   call dgemm('n','n',nrobs,nrens,nrmin,1.0,U,nrobs,             &
              X2,nrmin,0.0,X3,nrobs)
   deallocate(U)
   deallocate(X2)


   if (2_8*ndim*nrobs < 1_8*nrens*(nrobs+ndim)) then
!    Code for few observations ( m<nN/(2n-N) )
      if (verbose) print* ,' Representers computed'
      allocate (Reps(ndim,nrobs))    ! Representers
!      Reps=matmul(A,transpose(S))
      call dgemm('n','t',ndim,nrobs,nrens,1.0,A,ndim,S,nrobs,0.0,Reps,ndim)
!      A=A+matmul(Reps,X3)
      call dgemm('n','n',ndim,nrens,nrobs,1.0,Reps,ndim,X3,nrobs,1.0,A,ndim)
      deallocate(Reps)
   else
      if (verbose) print* ,' nrl-combination X5 computed'
      allocate(X4(nrens,nrens))
!      X4=matmul(transpose(S),X3)
      call dgemm('t','n',nrens,nrens,nrobs,1.0,S,nrobs,X3,nrobs,0.0,X4,nrens)
      do i=1,nrens
         X4(i,i)=X4(i,i)+1.0  ! Addition with Af
      enddo

      iblkmax=min(ndim,200)
      call  multa(A, X4, ndim, nrens, iblkmax )

! Output of global analysis nrl-combination matrix X5  ('not really linear')
!      open(10,file='assimilation.in')
!         read(10,'(l1)') global
!      close(10)

!      if(global.or.verbose) then
!        open(10,file='infile.data')
!           read(10,'(a)') id
!           read(10,'(a)') id
!           read(10,'(a)') id
!           read(10,'(a9)') id
!        close(10)

!        print *, ' Writing X5 in X5.uf'
        open(10,file='X5.uf',status='replace',form='unformatted')
!        print *, ' Writing X5 in X5_'//trim(id)//'.uf'
!        open(10,file='X5_'//trim(id)//'.uf',status='replace',form='unformatted')
        write(10) ((X4(i,j),i=1,nrens),j=1,nrens)
        close(10)

!        open(10,file='X5'//trim(id)//'.dat')
        open(10,file='X5.dat')
        write(10,*)'TITLE = "X5"'
        write(10,*)'VARIABLES = "iens_f" "iens_a" "X5"'
        write(10,'(2(a,i5),a)')' ZONE  F=BLOCK, I=',nrens,', J=',nrens,', K=1'
        write(10,'(30I5)')((i,i=1,nrens),j=1,nrens)
        write(10,'(30I5)')((j,i=1,nrens),j=1,nrens)
        write(10,900)((X4(i,j),i=1,nrens),j=1,nrens)
        close(10)
900 format(10(1x,e12.5))
!      endif
!later: Also compute determinant of X5 by eigendecomposition

      open(10,file='X5col.dat')
        do j=1,nrens
           write(10,'(i5,f10.4)')j,sum(X4(:,j))  ! should be one E(Xa)=E(Xf)
        enddo
      close(10)
      open(10,file='X5row.dat')
        do j=1,nrens
!       write(10,'(i5,g12.4)')j,sum(X4(j,:))
           write(10,'(i5,f10.4)')j,sum(X4(j,:))/float(nrens)
        enddo
      close(10)

      deallocate(X4)
   endif 

   if(verbose) print *,'Analysis, TEM(01), member 1 =', A(target,1)

end subroutine analysis
end module analysis_module
