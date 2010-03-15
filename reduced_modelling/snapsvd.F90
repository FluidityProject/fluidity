module snapsvd_module
#include "fdebug.h"
  implicit none

contains 
  subroutine snapsvd(m,n,snapmatrix,nsvd_total,nsvd,usv,svdval)

    !   
    !     This code shows how to use ARPACK to find a few of the
    !     largest singular values(sigma) and corresponding right singular 
    !     vectors (v) for the the matrix A by solving the symmetric problem:
    !          
    !                        (A'*A)*v = sigma*v
    ! 
    !     where A is an m by n real matrix.
    !
    !     This code may be easily modified to estimate the 2-norm
    !     condition number  largest(sigma)/smallest(sigma) by setting
    !     which = 'BE' below.  This will ask for a few of the smallest
    !     and a few of the largest singular values simultaneously.
    !     The condition number could then be estimated by taking
    !     the ratio of the largest and smallest singular values.
    !
    !     This formulation is appropriate when  m  .ge.  n.
    !     Reverse the roles of A and A' in the case that  m .le. n.
    !
    !     The main points illustrated here are 
    !
    !        1) How to declare sufficient memory to find NEV 
    !           largest singular values of A .  
    !
    !        2) Illustration of the reverse communication interface 
    !           needed to utilize the top level ARPACK routine DSAUPD 
    !           that computes the quantities needed to construct
    !           the desired singular values and vectors(if requested).
    !
    !        3) How to extract the desired singular values and vectors
    !           using the ARPACK routine DSEUPD.
    !
    !        4) How to construct the left singular vectors U from the 
    !           right singular vectors V to obtain the decomposition
    !
    !                        A*V = U*S
    !
    !           where S = diag(sigma_1, sigma_2, ..., sigma_k).
    !
    !     The only thing that must be supplied in order to use this
    !     routine on your problem is to change the array dimensions 
    !     appropriately, to specify WHICH singular values you want to 
    !     compute and to supply a the matrix-vector products 
    !
    !                         w <-  Ax
    !                         y <-  A'w
    !
    !     in place of the calls  to AV( ) and ATV( ) respectively below.  
    !
    !     Further documentation is available in the header of DSAUPD
    !     which may be found in the SRC directory.
    !
    !     This codes implements
    !
    !\Example-1
    !     ... Suppose we want to solve A'A*v = sigma*v in regular mode,
    !         where A is derived from the simplest finite difference 
    !         discretization of the 2-dimensional kernel  K(s,t)dt  where
    !
    !                 K(s,t) =  s(t-1)   if 0 .le. s .le. t .le. 1,
    !                           t(s-1)   if 0 .le. t .lt. s .le. 1. 
    !
    !         See subroutines AV  and ATV for details.
    !     ... OP = A'*A  and  B = I.
    !     ... Assume "call av (n,x,y)" computes y = A*x
    !     ... Assume "call atv (n,y,w)" computes w = A'*y
    !     ... Assume exact shifts are used
    !     ...
    !
    !\BeginLib
    !
    !\Routines called:
    !     dsaupd  ARPACK reverse communication interface routine.
    !     dseupd  ARPACK routine that returns Ritz values and (optionally)
    !             Ritz vectors.
    !     dnrm2   Level 1 BLAS that computes the norm of a vector.
    !     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
    !     dscal   Level 1 BLAS thst computes x <- x*alpha.
    !     dcopy   Level 1 BLAS thst computes y <- x.
    !
    !-----------------------------------------------------------------------
    !
    !     %------------------------------------------------------%
    !     | Storage Declarations:                                |
    !     |                                                      |
    !     | It is assumed that A is M by N with M .ge. N.        |
    !     |                                                      |
    !     | The maximum dimensions for all arrays are            |
    !     | set here to accommodate a problem size of            |
    !     | M .le. MAXM  and  N .le. MAXN                        |
    !     |                                                      |
    !     | The NEV right singular vectors will be computed in   |
    !     | the N by NCV array V.                                |
    !     |                                                      |
    !     | The NEV left singular vectors will be computed in    |
    !     | the M by NEV array U.                                |
    !     |                                                      |
    !     | NEV is the number of singular values requested.      |
    !     |     See specifications for ARPACK usage below.       |
    !     |                                                      |
    !     | NCV is the largest number of basis vectors that will |
    !     |     be used in the Implicitly Restarted Arnoldi      |
    !     |     Process.  Work per major iteration is            |
    !     |     proportional to N*NCV*NCV.                       |
    !     |                                                      |
    !     | You must set:                                        |
    !     |                                                      |
    !     | MAXM:   Maximum number of rows of the A allowed.     |
    !     | MAXN:   Maximum number of columns of the A allowed.  |
    !     | MAXNEV: Maximum NEV allowed                          |
    !     | MAXNCV: Maximum NCV allowed                          |
    !     %------------------------------------------------------%
    !

    use FLDebug

    !include 'paramnew.h'

    integer   maxm, maxn, maxnev, maxncv, ldv, ldu, m, n, i
    integer   ndigit,logfil, msgets, msaitr, msapps, msaupd, msaup2, mseigt, mseupd

    !        parameter (maxm = istate,maxn=nrsnapshots,maxnev=nrsnapshots-1,    &
    !      maxncv=nrsnapshots, ldu = maxm, ldv=maxn )
    !     %--------------%
    !     | Local Arrays |
    !     %--------------%
    !
    !Double precision v(ldv,maxncv), u(ldu, maxnev),            & 
    !                 workl(maxncv*(maxncv+8)), workd(3*maxn),  &
    !                 s(maxncv,2), resid(maxn), ax(maxm)
    !logical          select(maxncv)
    integer          iparam(11), ipntr(11)

    REAL, ALLOCATABLE, DIMENSION(:,:):: v,u,s
    REAL, ALLOCATABLE, DIMENSION(:)::   workl,workd,resid,ax
    logical, ALLOCATABLE, DIMENSION(:)::   select

    !
    !     %---------------%
    !     | Local Scalars |
    !     %---------------%
    !
    character        bmat*1, which*2
    integer          ido, nev, ncv, lworkl, info, ierr,   &
         j, ishfts, maxitr, mode1, nconv
    logical          rvec
    Double precision toll, sigma, temp
    !
    !     %------------%
    !     | Parameters |
    !     %------------%
    !
    Double precision one, zero
    parameter        (one = 1.0D+0, zero = 0.0D+0)


    !ccccc declare input
    double precision snapmatrix(m,n)
    integer nsvd,nsvd_total
    !cccccccc  declare output
    double precision usv(m,nsvd),svdval(nsvd) ! left svd, svdval
    !
    !  
    !     %-----------------------------%
    !     | BLAS & LAPACK routines used |
    !     %-----------------------------%
    !
    Double precision dnrm2

    !
    !  %-------------%
    !  | local       |
    !  %-------------%

    REAL TOTAL_E,PARTIAL_E,PERCENT_E

    external         dnrm2, daxpy, dcopy, dscal
    !
    !     %-----------------------%
    !     | Executable Statements |
    !     %-----------------------%
    !
    !     %-------------------------------------------------%
    !     | The following include statement and assignments |
    !     | initiate trace output from the internal         |
    !     | actions of ARPACK.  See debug.doc in the        |
    !     | DOCUMENTS directory for usage.  Initially, the  |
    !     | most useful information will be a breakdown of  |
    !     | time spent in the various stages of computation |
    !     | given by setting msaupd = 1.                    |
    !     %-------------------------------------------------%
    !
    !      include 'debug1.h'
#ifdef HAVE_LIBARPACK

    maxm = m
    maxn=n
    maxnev=n-1
    maxncv=n
    ldu = maxm
    ldv=maxn

    ALLOCATE(v(ldv,maxncv))
    ALLOCATE(u(ldu,maxnev))
    ALLOCATE(s(maxncv,2))
    ALLOCATE(workl(maxncv*(maxncv+8)))
    ALLOCATE(workd(3*maxn))
    ALLOCATE(resid(maxn))
    ALLOCATE(ax(maxm))
    ALLOCATE(select(maxncv))

    ewrite(3,*) 'in snapsvd.F'
    ewrite(3,*) M,maxm
    !      maxm=m
    !      ldu = maxm
    ndigit = -3
    logfil = 6
    msgets = 0
    msaitr = 0 
    msapps = 0
    msaupd = 1
    msaup2 = 0
    mseigt = 0
    mseupd = 0
    !
    !     %-------------------------------------------------%
    !     | The following sets dimensions for this problem. |
    !     %-------------------------------------------------%
    !
    !      if (n .ne. nrsnapshots) then
    !       write(*,*) 'dimension error in svd routine',n,nrsnapshots
    !       stop
    !      end if
    !
    !     %------------------------------------------------%
    !     | Specifications for ARPACK usage are set        | 
    !     | below:                                         |
    !     |                                                |
    !     |    1) NEV = 4 asks for 4 singular values to be |  
    !     |       computed.                                | 
    !     |                                                |
    !     |    2) NCV = 20 sets the length of the Arnoldi  |
    !     |       factorization                            |
    !     |                                                |
    !     |    3) This is a standard problem               |
    !     |         (indicated by bmat  = 'I')             |
    !     |                                                |
    !     |    4) Ask for the NEV singular values of       |
    !     |       largest magnitude                        |
    !     |         (indicated by which = 'LM')            |
    !     |       See documentation in DSAUPD for the      |
    !     |       other options SM, BE.                    | 
    !     |                                                |
    !     | Note: NEV and NCV must satisfy the following   |
    !     |       conditions:                              |
    !     |                 NEV <= MAXNEV,                 |
    !     |             NEV + 1 <= NCV <= MAXNCV           |
    !     %------------------------------------------------%
    !
    nev   = nsvd_total
    ncv   = n
    bmat  = 'I'
    which = 'LM'
    !
    if ( n .gt. maxn ) then
       ewrite(3,*)  ' ERROR with _SVD: N is greater than MAXN '
       go to 9000
    else if ( m .gt. maxm ) then
       ewrite(3,*)  ' ERROR with _SVD: M is greater than MAXM ',M,MAXM
       go to 9000
    else if ( nev .gt. maxnev ) then
       ewrite(3,*)  ' ERROR with _SVD: NEV is greater than MAXNEV '
       go to 9000
    else if ( ncv .gt. maxncv ) then
       ewrite(3,*)  ' ERROR with _SVD: NCV is greater than MAXNCV '
       go to 9000
    end if
    !
    !     %-----------------------------------------------------%
    !     | Specification of stopping rules and initial         |
    !     | conditions before calling DSAUPD                    |
    !     |                                                     |
    !     |           abs(sigmaC - sigmaT) < TOL*abs(sigmaC)    |
    !     |               computed   true                       |
    !     |                                                     |
    !     |      If TOL .le. 0,  then TOL <- macheps            |
    !     |              (machine precision) is used.           |
    !     |                                                     |
    !     | IDO  is the REVERSE COMMUNICATION parameter         |
    !     |      used to specify actions to be taken on return  |
    !     |      from DSAUPD. (See usage below.)                |
    !     |                                                     |
    !     |      It MUST initially be set to 0 before the first |
    !     |      call to DSAUPD.                                | 
    !     |                                                     |
    !     | INFO on entry specifies starting vector information |
    !     |      and on return indicates error codes            |
    !     |                                                     |
    !     |      Initially, setting INFO=0 indicates that a     | 
    !     |      random starting vector is requested to         |
    !     |      start the ARNOLDI iteration.  Setting INFO to  |
    !     |      a nonzero value on the initial call is used    |
    !     |      if you want to specify your own starting       |
    !     |      vector (This vector must be placed in RESID.)  | 
    !     |                                                     |
    !     | The work array WORKL is used in DSAUPD as           | 
    !     | workspace.  Its dimension LWORKL is set as          |
    !     | illustrated below.                                  |
    !     %-----------------------------------------------------%
    !
    lworkl = ncv*(ncv+8)
    toll = 0.0d-3 ! 1.d-3 
    info = 0
    ido = 0
    !
    !     %---------------------------------------------------%
    !     | Specification of Algorithm Mode:                  |
    !     |                                                   |
    !     | This program uses the exact shift strategy        |
    !     | (indicated by setting IPARAM(1) = 1.)             |
    !     | IPARAM(3) specifies the maximum number of Arnoldi |
    !     | iterations allowed.  Mode 1 of DSAUPD is used     |
    !     | (IPARAM(7) = 1). All these options can be changed |
    !     | by the user. For details see the documentation in |
    !     | DSAUPD.                                           |
    !     %---------------------------------------------------%
    !
    ishfts = 1
    maxitr = n
    mode1 = 1
    !
    iparam(1) = ishfts
    !                
    iparam(3) = maxitr
    !                  
    iparam(7) = mode1
    !
    !     %------------------------------------------------%
    !     | M A I N   L O O P (Reverse communication loop) |
    !     %------------------------------------------------%
    !
10  continue
    !
    !        %---------------------------------------------%
    !        | Repeatedly call the routine DSAUPD and take | 
    !        | actions indicated by parameter IDO until    |
    !        | either convergence is indicated or maxitr   |
    !        | has been exceeded.                          |
    !        %---------------------------------------------%
    !
    call dsaupd ( ido, bmat, n, which, nev, toll, resid,      &
         ncv, v, ldv, iparam, ipntr, workd, workl,   &
         lworkl, info )
    !
    if (ido .eq. -1 .or. ido .eq. 1) then
       !
       !           %---------------------------------------%
       !           | Perform matrix vector multiplications |
       !           |              w <--- A*x       (av())  |
       !           |              y <--- A'*w      (atv()) |
       !           | The user should supply his/her own    |
       !           | matrix vector multiplication routines |
       !           | here that takes workd(ipntr(1)) as    |
       !           | the input, and returns the result in  |
       !           | workd(ipntr(2)).                      |
       !           %---------------------------------------%
       !
       call av (m,n,snapmatrix,workd(ipntr(1)),ax) 
       call atv (m,n,snapmatrix,ax,workd(ipntr(2)))
       !
       !           %-----------------------------------------%
       !           | L O O P   B A C K to call DSAUPD again. |
       !           %-----------------------------------------%
       !
       go to 10
       !
    end if
    !
    !     %----------------------------------------%
    !     | Either we have convergence or there is |
    !     | an error.                              |
    !     %----------------------------------------%
    !
    if ( info .lt. 0 ) then
       !
       !        %--------------------------%
       !        | Error message. Check the |
       !        | documentation in DSAUPD. |
       !        %--------------------------%
       !
       ewrite(3,*)  ' '
       ewrite(3,*)  ' Error with _saupd, info = ', info
       ewrite(3,*)  ' Check documentation in _saupd '
       ewrite(3,*)  ' '
       !
    else 
       !
       !        %--------------------------------------------%
       !        | No fatal errors occurred.                  |
       !        | Post-Process using DSEUPD.                 |
       !        |                                            |
       !        | Computed singular values may be extracted. |  
       !        |                                            |
       !        | Singular vectors may also be computed now  |
       !        | if desired.  (indicated by rvec = .true.)  | 
       !        |                                            |
       !        | The routine DSEUPD now called to do this   |
       !        | post processing                            | 
       !        %--------------------------------------------%
       !           
       rvec = .true.
       !
       call dseupd ( rvec, 'All', select, s, v, ldv, sigma,   &
            bmat, n, which, nev, toll, resid, ncv, v, ldv,    &
            iparam, ipntr, workd, workl, lworkl, ierr )
       !
       !        %-----------------------------------------------%
       !        | Singular values are returned in the first     |
       !        | column of the two dimensional array S         |
       !        | and the corresponding right singular vectors  | 
       !        | are returned in the first NEV columns of the  |
       !        | two dimensional array V as requested here.    |
       !        %-----------------------------------------------%
       !

       ewrite(3,*) ipntr(7)
       ewrite(3,*) ipntr(7)+ncv
       ewrite(3,*) (workl(ipntr(7)+ncv+j-1),j=1,ncv)

       TOTAL_E=0.0
       DO J=1,NCV
          IF(workl(ipntr(7)+ncv+j-1).ge.0.0) THEN
             TOTAL_E=TOTAL_E+sqrt(workl(ipntr(7)+ncv+j-1))
          ENDIF
       ENDDO

       if ( ierr .ne. 0) then
          !
          !           %------------------------------------%
          !           | Error condition:                   |
          !           | Check the documentation of DSEUPD. |
          !           %------------------------------------%
          !
          ewrite(3,*)  ' '
          ewrite(3,*)  ' Error with _seupd, info = ', ierr
          ewrite(3,*)  ' Check the documentation of _seupd. '
          ewrite(3,*)  ' '
          !
       else
          !
          nconv =  iparam(5)
          ewrite(3,*) 'nconv',nconv
          PARTIAL_E=0.0
          do 20 j=1, nconv
             s(j,1) = sqrt(s(j,1))
             if(j.le.nsvd) then
                PARTIAL_E=PARTIAL_E+s(j,1)
             endif
             !
             !              %-----------------------------%
             !              | Compute the left singular   |
             !              | vectors from the formula    |
             !              |                             |
             !              |     u = Av/sigma            |
             !              |                             |
             !              | u should have norm 1 so     |
             !              | divide by norm(Av) instead. |
             !              %-----------------------------%
             !
             call av(m, n, snapmatrix, v(1,j), ax)
             call dcopy(m, ax, 1, u(1,j), 1)
             temp = one/dnrm2(m, u(1,j), 1)
             call dscal(m, temp, u(1,j), 1)
             !
             !              %---------------------------%
             !              |                           |
             !              | Compute the residual norm |
             !              |                           |
             !              |   ||  A*v - sigma*u ||    |
             !              |                           |
             !              | for the NCONV accurately  |
             !              | computed singular values  |
             !              | and vectors.  (iparam(5)  |
             !              | indicates how many are    |
             !              | accurate to the requested |
             !              | tolerance).               |
             !              | Store the result in 2nd   |
             !              | column of array S.        |
             !              %---------------------------%
             !
             call daxpy(m, -s(j,1), u(1,j), 1, ax, 1)
             s(j,2) = dnrm2(m, ax, 1)
             !
20        end do

          PERCENT_E = 100.0*PARTIAL_E/TOTAL_E
          ewrite(3,*)  PARTIAL_E,TOTAL_E
          EWRITE(3,*) 'The energry captured(%)',PERCENT_E
          !
          !           %-------------------------------%
          !           | Display computed residuals    |
          !           %-------------------------------%
          call dmout(6, nconv, 2, s, maxncv, -6,           &
               'Singular values and direct residuals')
       end if
       !
       !        %------------------------------------------%
       !        | Print additional convergence information |
       !        %------------------------------------------%
       !
       if ( info .eq. 1) then
          ewrite(3,*)  ' '
          ewrite(3,*)  ' Maximum number of iterations reached.'
          ewrite(3,*)  ' '
       else if ( info .eq. 3) then
          ewrite(3,*)  ' ' 
          ewrite(3,*)  ' No shifts could be applied during implicit',   &
               ' Arnoldi update, try increasing NCV.'
          ewrite(3,*)  ' '
       end if
       !
       ewrite(3,*)  ' '
       ewrite(3,*)  ' _SVD '
       ewrite(3,*)  ' ==== '
       ewrite(3,*)  ' '
       ewrite(3,*)  ' Size of the matrix is ', n
       ewrite(3,*)  ' The number of Ritz values requested is ', nev
       ewrite(3,*)  ' The number of Arnoldi vectors generated',     &
            ' (NCV) is ', ncv
       ewrite(3,*)  ' What portion of the spectrum: ', which
       ewrite(3,*)  ' The number of converged Ritz values is ',     &
            nconv 
       ewrite(3,*)  ' The number of Implicit Arnoldi update',       &
            ' iterations taken is ', iparam(3)
       ewrite(3,*)  ' The number of OP*x is ', iparam(9)
       ewrite(3,*)  ' The convergence criterion is ', toll
       ewrite(3,*)  ' '
       !
    end if
    !
    !     %-------------------------%
    !     | Done with program dsvd. |
    !     %-------------------------%
    !
9000 continue

    if (nconv .ne. nsvd_total) then
       write(*,*) 'convergence error in svd computation'
       stop
    end if

    do j = nsvd,1,-1
       svdval(nconv+1-j)= s(j,1)  ! decreasing order
    enddo


    do j = nsvd,1,-1        
       do i=1,m       
          usv(i,nconv+1-j) = u(i,j)
       enddo
    enddo
    ewrite(3,*) 'the end of snapsvd'
1000 format (1X, E24.16, 1X, E24.16)
2000 format (100(1X,E24.16)) 

    DEALLOCATE(v)
    DEALLOCATE(u)
    DEALLOCATE(s)
    DEALLOCATE(workl)
    DEALLOCATE(workd)
    DEALLOCATE(resid)
    DEALLOCATE(ax)
    DEALLOCATE(select)

#else
    FLAbort("No ARPACK support - cannot run reduced model.")
#endif
    return
  end subroutine snapsvd


  !*******************H1***********************
  subroutine snapsvd_H1(m,n,epsilon,snapmatrix,snapmatrix_dx,snapmatrix_dy,snapmatrix_dz,nsvd_total,nsvd,usv,svdval,D3)

    !   
    !     This code shows how to use ARPACK to find a few of the
    !     largest singular values(sigma) and corresponding right singular 
    !     vectors (v) for the the matrix A by solving the symmetric problem:
    !          
    !                        (A'*A)*v = sigma*v
    ! 
    !     where A is an m by n real matrix.
    !
    !     This code may be easily modified to estimate the 2-norm
    !     condition number  largest(sigma)/smallest(sigma) by setting
    !     which = 'BE' below.  This will ask for a few of the smallest
    !     and a few of the largest singular values simultaneously.
    !     The condition number could then be estimated by taking
    !     the ratio of the largest and smallest singular values.
    !
    !     This formulation is appropriate when  m  .ge.  n.
    !     Reverse the roles of A and A' in the case that  m .le. n.
    !
    !     The main points illustrated here are 
    !
    !        1) How to declare sufficient memory to find NEV 
    !           largest singular values of A .  
    !
    !        2) Illustration of the reverse communication interface 
    !           needed to utilize the top level ARPACK routine DSAUPD 
    !           that computes the quantities needed to construct
    !           the desired singular values and vectors(if requested).
    !
    !        3) How to extract the desired singular values and vectors
    !           using the ARPACK routine DSEUPD.
    !
    !        4) How to construct the left singular vectors U from the 
    !           right singular vectors V to obtain the decomposition
    !
    !                        A*V = U*S
    !
    !           where S = diag(sigma_1, sigma_2, ..., sigma_k).
    !
    !     The only thing that must be supplied in order to use this
    !     routine on your problem is to change the array dimensions 
    !     appropriately, to specify WHICH singular values you want to 
    !     compute and to supply a the matrix-vector products 
    !
    !                         w <-  Ax
    !                         y <-  A'w
    !
    !     in place of the calls  to AV( ) and ATV( ) respectively below.  
    !
    !     Further documentation is available in the header of DSAUPD
    !     which may be found in the SRC directory.
    !
    !     This codes implements
    !
    !\Example-1
    !     ... Suppose we want to solve A'A*v = sigma*v in regular mode,
    !         where A is derived from the simplest finite difference 
    !         discretization of the 2-dimensional kernel  K(s,t)dt  where
    !
    !                 K(s,t) =  s(t-1)   if 0 .le. s .le. t .le. 1,
    !                           t(s-1)   if 0 .le. t .lt. s .le. 1. 
    !
    !         See subroutines AV  and ATV for details.
    !     ... OP = A'*A  and  B = I.
    !     ... Assume "call av (n,x,y)" computes y = A*x
    !     ... Assume "call atv (n,y,w)" computes w = A'*y
    !     ... Assume exact shifts are used
    !     ...
    !
    !\BeginLib
    !
    !\Routines called:
    !     dsaupd  ARPACK reverse communication interface routine.
    !     dseupd  ARPACK routine that returns Ritz values and (optionally)
    !             Ritz vectors.
    !     dnrm2   Level 1 BLAS that computes the norm of a vector.
    !     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
    !     dscal   Level 1 BLAS thst computes x <- x*alpha.
    !     dcopy   Level 1 BLAS thst computes y <- x.
    !
    !-----------------------------------------------------------------------
    !
    !     %------------------------------------------------------%
    !     | Storage Declarations:                                |
    !     |                                                      |
    !     | It is assumed that A is M by N with M .ge. N.        |
    !     |                                                      |
    !     | The maximum dimensions for all arrays are            |
    !     | set here to accommodate a problem size of            |
    !     | M .le. MAXM  and  N .le. MAXN                        |
    !     |                                                      |
    !     | The NEV right singular vectors will be computed in   |
    !     | the N by NCV array V.                                |
    !     |                                                      |
    !     | The NEV left singular vectors will be computed in    |
    !     | the M by NEV array U.                                |
    !     |                                                      |
    !     | NEV is the number of singular values requested.      |
    !     |     See specifications for ARPACK usage below.       |
    !     |                                                      |
    !     | NCV is the largest number of basis vectors that will |
    !     |     be used in the Implicitly Restarted Arnoldi      |
    !     |     Process.  Work per major iteration is            |
    !     |     proportional to N*NCV*NCV.                       |
    !     |                                                      |
    !     | You must set:                                        |
    !     |                                                      |
    !     | MAXM:   Maximum number of rows of the A allowed.     |
    !     | MAXN:   Maximum number of columns of the A allowed.  |
    !     | MAXNEV: Maximum NEV allowed                          |
    !     | MAXNCV: Maximum NCV allowed                          |
    !     %------------------------------------------------------%
    !

    use FLDebug

    !IMPLICIT NONE
    !include 'paramnew.h'

    integer   maxm, maxn, maxnev, maxncv, ldv, ldu,m,n,i
    integer   ndigit,logfil, msgets, msaitr, msapps, msaupd, msaup2, mseigt, mseupd

    !        parameter (maxm = istate,maxn=nrsnapshots,maxnev=nrsnapshots-1,    &
    !      maxncv=nrsnapshots, ldu = maxm, ldv=maxn )
    !     %--------------%
    !     | Local Arrays |
    !     %--------------%
    !
    !Double precision v(ldv,maxncv), u(ldu, maxnev),            & 
    !                 workl(maxncv*(maxncv+8)), workd(3*maxn),  &
    !                 s(maxncv,2), resid(maxn), ax(maxm)
    !logical          select(maxncv)
    integer          iparam(11), ipntr(11)

    REAL, ALLOCATABLE, DIMENSION(:,:):: v,u,s
    REAL, ALLOCATABLE, DIMENSION(:)::   workl,workd,resid,ax
    REAL, ALLOCATABLE, DIMENSION(:)::   ay,ay_dx,ay_dy,ay_dz
    logical, ALLOCATABLE, DIMENSION(:)::   select
    logical ::D3

    !
    !     %---------------%
    !     | Local Scalars |
    !     %---------------%
    !
    character        bmat*1, which*2
    integer          ido, nev, ncv, lworkl, info, ierr,   &
         j, ishfts, maxitr, mode1, nconv
    logical          rvec
    Double precision toll, sigma, temp
    !
    !     %------------%
    !     | Parameters |
    !     %------------%
    !
    Double precision one, zero
    parameter        (one = 1.0D+0, zero = 0.0D+0)


    !ccccc declare input
    double precision snapmatrix(m,n)
    double precision snapmatrix_dx(m,n),snapmatrix_dy(m,n),snapmatrix_dz(m,n)
    real epsilon
    integer nsvd,nsvd_total
    !cccccccc  declare output
    double precision usv(m,nsvd),svdval(nsvd) ! left svd, svdval
    !
    !  
    !     %-----------------------------%
    !     | BLAS & LAPACK routines used |
    !     %-----------------------------%
    !
    Double precision dnrm2

    !
    !  %-------------%
    !  | local       |
    !  %-------------%

    REAL TOTAL_E,PARTIAL_E,PERCENT_E

    external         dnrm2, daxpy, dcopy, dscal
    !
    !     %-----------------------%
    !     | Executable Statements |
    !     %-----------------------%
    !
    !     %-------------------------------------------------%
    !     | The following include statement and assignments |
    !     | initiate trace output from the internal         |
    !     | actions of ARPACK.  See debug.doc in the        |
    !     | DOCUMENTS directory for usage.  Initially, the  |
    !     | most useful information will be a breakdown of  |
    !     | time spent in the various stages of computation |
    !     | given by setting msaupd = 1.                    |
    !     %-------------------------------------------------%
    !
    !      include 'debug1.h'
#ifdef HAVE_LIBARPACK

    maxm = m
    maxn=n
    maxnev=n
    maxncv=n
    ldu = maxm
    ldv=maxn

    ALLOCATE(v(ldv,maxncv))
    ALLOCATE(u(ldu, maxnev))
    ALLOCATE(s(maxncv,2))
    ALLOCATE(workl(maxncv*(maxncv+8)))
    ALLOCATE(workd(3*maxn))
    ALLOCATE(resid(maxn))
    ALLOCATE(ax(maxm))
    ALLOCATE(select(maxncv))
    ALLOCATE(ay(maxn))
    ALLOCATE(ay_dx(maxn))
    ALLOCATE(ay_dy(maxn))
    ALLOCATE(ay_dz(maxn))

    ewrite(3,*) 'in snapsvd.F'
    ewrite(3,*) M,maxm
    !      maxm=m
    !      ldu = maxm
    ndigit = -3
    logfil = 6
    msgets = 0
    msaitr = 0 
    msapps = 0
    msaupd = 1
    msaup2 = 0
    mseigt = 0
    mseupd = 0
    ay=0.0
    ay_dx=0.0
    ay_dy=0.0
    ay_dz=0.0
    !
    !     %-------------------------------------------------%
    !     | The following sets dimensions for this problem. |
    !     %-------------------------------------------------%
    !
    !      if (n .ne. nrsnapshots) then
    !       write(*,*) 'dimension error in svd routine',n,nrsnapshots
    !       stop
    !      end if
    !
    !     %------------------------------------------------%
    !     | Specifications for ARPACK usage are set        | 
    !     | below:                                         |
    !     |                                                |
    !     |    1) NEV = 4 asks for 4 singular values to be |  
    !     |       computed.                                | 
    !     |                                                |
    !     |    2) NCV = 20 sets the length of the Arnoldi  |
    !     |       factorization                            |
    !     |                                                |
    !     |    3) This is a standard problem               |
    !     |         (indicated by bmat  = 'I')             |
    !     |                                                |
    !     |    4) Ask for the NEV singular values of       |
    !     |       largest magnitude                        |
    !     |         (indicated by which = 'LM')            |
    !     |       See documentation in DSAUPD for the      |
    !     |       other options SM, BE.                    | 
    !     |                                                |
    !     | Note: NEV and NCV must satisfy the following   |
    !     |       conditions:                              |
    !     |                 NEV <= MAXNEV,                 |
    !     |             NEV + 1 <= NCV <= MAXNCV           |
    !     %------------------------------------------------%
    !
    nev   = nsvd_total
    ncv   = n
    bmat  = 'I'
    which = 'LM'
    !
    if ( n .gt. maxn ) then
       ewrite(3,*)  ' ERROR with _SVD: N is greater than MAXN '
       go to 9000
    else if ( m .gt. maxm ) then
       ewrite(3,*)  ' ERROR with _SVD: M is greater than MAXM ',M,MAXM
       go to 9000
    else if ( nev .gt. maxnev ) then
       ewrite(3,*)  ' ERROR with _SVD: NEV is greater than MAXNEV '
       go to 9000
    else if ( ncv .gt. maxncv ) then
       ewrite(3,*)  ' ERROR with _SVD: NCV is greater than MAXNCV '
       go to 9000
    end if
    !
    !     %-----------------------------------------------------%
    !     | Specification of stopping rules and initial         |
    !     | conditions before calling DSAUPD                    |
    !     |                                                     |
    !     |           abs(sigmaC - sigmaT) < TOL*abs(sigmaC)    |
    !     |               computed   true                       |
    !     |                                                     |
    !     |      If TOL .le. 0,  then TOL <- macheps            |
    !     |              (machine precision) is used.           |
    !     |                                                     |
    !     | IDO  is the REVERSE COMMUNICATION parameter         |
    !     |      used to specify actions to be taken on return  |
    !     |      from DSAUPD. (See usage below.)                |
    !     |                                                     |
    !     |      It MUST initially be set to 0 before the first |
    !     |      call to DSAUPD.                                | 
    !     |                                                     |
    !     | INFO on entry specifies starting vector information |
    !     |      and on return indicates error codes            |
    !     |                                                     |
    !     |      Initially, setting INFO=0 indicates that a     | 
    !     |      random starting vector is requested to         |
    !     |      start the ARNOLDI iteration.  Setting INFO to  |
    !     |      a nonzero value on the initial call is used    |
    !     |      if you want to specify your own starting       |
    !     |      vector (This vector must be placed in RESID.)  | 
    !     |                                                     |
    !     | The work array WORKL is used in DSAUPD as           | 
    !     | workspace.  Its dimension LWORKL is set as          |
    !     | illustrated below.                                  |
    !     %-----------------------------------------------------%
    !
    lworkl = ncv*(ncv+8)
    toll = 0.0d-3 ! 1.d-3 
    info = 0
    ido = 0
    !
    !     %---------------------------------------------------%
    !     | Specification of Algorithm Mode:                  |
    !     |                                                   |
    !     | This program uses the exact shift strategy        |
    !     | (indicated by setting IPARAM(1) = 1.)             |
    !     | IPARAM(3) specifies the maximum number of Arnoldi |
    !     | iterations allowed.  Mode 1 of DSAUPD is used     |
    !     | (IPARAM(7) = 1). All these options can be changed |
    !     | by the user. For details see the documentation in |
    !     | DSAUPD.                                           |
    !     %---------------------------------------------------%
    !
    ishfts = 1
    maxitr = n
    mode1 = 1
    !
    iparam(1) = ishfts
    !                
    iparam(3) = maxitr
    !                  
    iparam(7) = mode1
    !
    !     %------------------------------------------------%
    !     | M A I N   L O O P (Reverse communication loop) |
    !     %------------------------------------------------%
    !
10  continue
    !
    !        %---------------------------------------------%
    !        | Repeatedly call the routine DSAUPD and take | 
    !        | actions indicated by parameter IDO until    |
    !        | either convergence is indicated or maxitr   |
    !        | has been exceeded.                          |
    !        %---------------------------------------------%
    !
    call dsaupd ( ido, bmat, n, which, nev, toll, resid,      &
         ncv, v, ldv, iparam, ipntr, workd, workl,   &
         lworkl, info )
    !
    if (ido .eq. -1 .or. ido .eq. 1) then
       !
       !           %---------------------------------------%
       !           | Perform matrix vector multiplications |
       !           |              w <--- A*x       (av())  |
       !           |              y <--- A'*w      (atv()) |
       !           | The user should supply his/her own    |
       !           | matrix vector multiplication routines |
       !           | here that takes workd(ipntr(1)) as    |
       !           | the input, and returns the result in  |
       !           | workd(ipntr(2)).                      |
       !           %---------------------------------------%
       !
       call av (m,n,snapmatrix,workd(ipntr(1)),ax(1:m)) 
       call atv (m,n,snapmatrix,ax,ay)

       call av (m,n,snapmatrix_dx,workd(ipntr(1)),ax(1:m)) 
       call atv (m,n,snapmatrix_dx,ax,ay_dx)

       call av (m,n,snapmatrix_dy,workd(ipntr(1)),ax(1:m)) 
       call atv (m,n,snapmatrix_dy,ax,ay_dy)

       if(.false.) then
          call av (m,n,snapmatrix_dz,workd(ipntr(1)),ax(1:m)) 
          call atv (m,n,snapmatrix_dz,ax,ay_dz)
       endif
       open(1,file='left.dat')
       write(1,*) ay(1:n)
       close(1)
       open(1,file='left_dx.dat')
       write(1,*) ay_dx(1:n)
       close(1)
       open(1,file='left_dy.dat')
       write(1,*) ay_dy(1:n)
       close(1)
       open(1,file='left_dz.dat')
       write(1,*) ay_dz(1:n)
       close(1)
       !     stop 1
       workd(ipntr(2):ipntr(2)+n-1)=ay(1:n)+epsilon*(ay_dx(1:n)+ay_dy(1:n)+ay_dz(1:n))
       !
       !           %-----------------------------------------%
       !           | L O O P   B A C K to call DSAUPD again. |
       !           %-----------------------------------------%
       !
       go to 10
       !
    end if
    !
    !     %----------------------------------------%
    !     | Either we have convergence or there is |
    !     | an error.                              |
    !     %----------------------------------------%
    !
    if ( info .lt. 0 ) then
       !
       !        %--------------------------%
       !        | Error message. Check the |
       !        | documentation in DSAUPD. |
       !        %--------------------------%
       !
       ewrite(3,*)  ' '
       ewrite(3,*)  ' Error with _saupd, info = ', info
       ewrite(3,*)  ' Check documentation in _saupd '
       ewrite(3,*)  ' '
       !
    else 
       !
       !        %--------------------------------------------%
       !        | No fatal errors occurred.                  |
       !        | Post-Process using DSEUPD.                 |
       !        |                                            |
       !        | Computed singular values may be extracted. |  
       !        |                                            |
       !        | Singular vectors may also be computed now  |
       !        | if desired.  (indicated by rvec = .true.)  | 
       !        |                                            |
       !        | The routine DSEUPD now called to do this   |
       !        | post processing                            | 
       !        %--------------------------------------------%
       !           
       rvec = .true.
       !
       call dseupd ( rvec, 'All', select, s, v, ldv, sigma,   &
            bmat, n, which, nev, toll, resid, ncv, v, ldv,    &
            iparam, ipntr, workd, workl, lworkl, ierr )
       !
       !        %-----------------------------------------------%
       !        | Singular values are returned in the first     |
       !        | column of the two dimensional array S         |
       !        | and the corresponding right singular vectors  | 
       !        | are returned in the first NEV columns of the  |
       !        | two dimensional array V as requested here.    |
       !        %-----------------------------------------------%
       !

       ewrite(3,*) ipntr(7)
       ewrite(3,*) ipntr(7)+ncv
       ewrite(3,*) (workl(ipntr(7)+ncv+j-1),j=1,ncv)

       TOTAL_E=0.0
       DO J=1,NCV
          IF(workl(ipntr(7)+ncv+j-1).ge.0.0) THEN
             TOTAL_E=TOTAL_E+sqrt(workl(ipntr(7)+ncv+j-1))
          ENDIF
       ENDDO
       if ( ierr .ne. 0) then
          !
          !           %------------------------------------%
          !           | Error condition:                   |
          !           | Check the documentation of DSEUPD. |
          !           %------------------------------------%
          !
          ewrite(3,*)  ' '
          ewrite(3,*)  ' Error with _seupd, info = ', ierr
          ewrite(3,*)  ' Check the documentation of _seupd. '
          ewrite(3,*)  ' '
          !
       else
          !
          nconv =  iparam(5)
          ewrite(3,*) 'nconv',nconv
          PARTIAL_E=0.0
          do 20 j=1, nconv
             s(j,1) = sqrt(s(j,1))
             if(j.le.nsvd) then
                PARTIAL_E=PARTIAL_E+s(j,1)
             endif
             !
             !              %-----------------------------%
             !              | Compute the left singular   |
             !              | vectors from the formula    |
             !              |                             |
             !              |     u = Av/sigma            |
             !              |                             |
             !              | u should have norm 1 so     |
             !              | divide by norm(Av) instead. |
             !              %-----------------------------%
             !
             call av(m, n, snapmatrix, v(1,j), ax)
             call dcopy(m, ax, 1, u(1,j), 1)
             temp = one/dnrm2(m, u(1,j), 1)
             call dscal(m, temp, u(1,j), 1)
             !
             !              %---------------------------%
             !              |                           |
             !              | Compute the residual norm |
             !              |                           |
             !              |   ||  A*v - sigma*u ||    |
             !              |                           |
             !              | for the NCONV accurately  |
             !              | computed singular values  |
             !              | and vectors.  (iparam(5)  |
             !              | indicates how many are    |
             !              | accurate to the requested |
             !              | tolerance).               |
             !              | Store the result in 2nd   |
             !              | column of array S.        |
             !              %---------------------------%
             !
             call daxpy(m, -s(j,1), u(1,j), 1, ax, 1)
             s(j,2) = dnrm2(m, ax, 1)
             !
20        end do

          PERCENT_E = 100.0*PARTIAL_E/TOTAL_E
          ewrite(3,*)  PARTIAL_E,TOTAL_E
          EWRITE(3,*) 'The energry captured(%)',PERCENT_E
          !
          !           %-------------------------------%
          !           | Display computed residuals    |
          !           %-------------------------------%
          call dmout(6, nconv, 2, s, maxncv, -6,           &
               'Singular values and direct residuals')
       end if
       !
       !        %------------------------------------------%
       !        | Print additional convergence information |
       !        %------------------------------------------%
       !
       if ( info .eq. 1) then
          ewrite(3,*)  ' '
          ewrite(3,*)  ' Maximum number of iterations reached.'
          ewrite(3,*)  ' '
       else if ( info .eq. 3) then
          ewrite(3,*)  ' ' 
          ewrite(3,*)  ' No shifts could be applied during implicit',   &
               ' Arnoldi update, try increasing NCV.'
          ewrite(3,*)  ' '
       end if
       !
       ewrite(3,*)  ' '
       ewrite(3,*)  ' _SVD '
       ewrite(3,*)  ' ==== '
       ewrite(3,*)  ' '
       ewrite(3,*)  ' Size of the matrix is ', n
       ewrite(3,*)  ' The number of Ritz values requested is ', nev
       ewrite(3,*)  ' The number of Arnoldi vectors generated',     &
            ' (NCV) is ', ncv
       ewrite(3,*)  ' What portion of the spectrum: ', which
       ewrite(3,*)  ' The number of converged Ritz values is ',     &
            nconv 
       ewrite(3,*)  ' The number of Implicit Arnoldi update',       &
            ' iterations taken is ', iparam(3)
       ewrite(3,*)  ' The number of OP*x is ', iparam(9)
       ewrite(3,*)  ' The convergence criterion is ', toll
       ewrite(3,*)  ' '
       !
    end if
    !
    !     %-------------------------%
    !     | Done with program dsvd. |
    !     %-------------------------%
    !
9000 continue

    if (nconv .ne. nsvd_total) then
       write(*,*) 'convergence error in svd computation'
       stop
    end if

    do j = nsvd,1,-1
       svdval(nconv+1-j)= s(j,1)  ! decreasing order
    enddo


    do j = nsvd,1,-1        
       do i=1,m       
          usv(i,nconv+1-j) = u(i,j)
       enddo
    enddo
    ewrite(3,*) 'the end of snapsvd'
1000 format (1X, E24.16, 1X, E24.16)
2000 format (100(1X,E24.16)) 

    DEALLOCATE(v)
    DEALLOCATE(u)
    DEALLOCATE(s)
    DEALLOCATE(workl)
    DEALLOCATE(workd)
    DEALLOCATE(resid)
    DEALLOCATE(ax)
    DEALLOCATE(select)
    DEALLOCATE(ay)
    DEALLOCATE(ay_dx)
    DEALLOCATE(ay_dy)
    DEALLOCATE(ay_dz)

#else
    FLAbort("No ARPACK support - cannot run reduced model.")
#endif
    return
  end subroutine snapsvd_H1
  !********************


  ! 
  ! ------------------------------------------------------------------
  !     matrix vector subroutines
  !
  !-------------------------------------------------------------------
  !
  subroutine av (m, n, snapmatrix, x, w)
    !
    !     computes  w <- A*x

    !      include 'paramnew.h'

    double precision x(n), w(m), psum
    integer i,j,m,n
    double precision snapmatrix(m,n) 


    !      if (n .ne. nrsnapshots) then
    !         write(*,*)'matrix dimension error: no of snapshots'
    !         stop                     
    !      end if

    !      if(m .ne. istate) then 
    !         write(*,*)'matrix dimension error: state dimension'
    !         stop                     
    !      end if

    do i=1,m
       psum = 0.0d0
       do j=1,n
          psum = psum + snapmatrix(i,j)*x(j)
       enddo
       w(i) = psum
    enddo

    return
  end subroutine av
  !
  !-------------------------------------------------------------------
  !
  subroutine atv (m, n, snapmatrix, w, y)
    !
    !     computes  y <- A'*w

    !      include 'paramnew.h'
    integer         m, n
    double precision w(m), y(n),psum

    integer i,j
    double precision snapmatrix(m,n) 

    !      if (n .ne. nrsnapshots .or. m .ne. istate) then
    !         write(*,*)'transpose matrix dimension error'
    !         stop                     
    !      end if

    do j=1,n
       psum = 0.0d0
       do i=1,m
          psum = psum + snapmatrix(i,j)*w(i)
       enddo
       y(j) = psum
    enddo

    return
  end subroutine atv

end module snapsvd_module
