C     ******************************************************************
C     ******************************************************************

      subroutine moresor(n,g,bnnz,blin,bcol,bval,bdiag,delta,sigma1,
     +sigma2,eps,maxit,l,pd,p,chcnt,memfail,msinfo)
      
      implicit none

C     SCALAR ARGUMENTS
      logical memfail,pd
      integer bnnz,chcnt,maxit,msinfo,n
      double precision delta,eps,l,sigma1,sigma2

C     ARRAY ARGUMENTS
      integer bcol(bnnz),bdiag(n),blin(bnnz)
      double precision bval(bnnz),g(n),p(n)

C     Solves the problem
C
C     minimize    psi(w) = 1/2 w^TBw + g^Tw 
C     subject to  ||w|| <= delta
C
C     Using the method described in "Computing a trust region step", 
C     by More and Sorensen.

C     msinfo:
C
C     0: both g and B are null
C     1: first convergence criterion is satisfied
C     2: second convergence criterion is satisfied
C     3: third convergence criterion is satisfied
C     5: maximum allowed number of iterations is achieved

#include "dim.par"
#include "machconst.com"
#include "outtyp.com"

C     LOCAL SCALARS
      integer col,dum,i,idx,iter,lin,lssinfo
      double precision b1n,d,delta2,geucn,ll,llant,ls,lsant,lu,luant,
     +        peucn,peucn2,ptz,rpeucn2,rzeucn2,sgn,tau,teucn2,tmp,ueucn2

C     LOCAL ARRAYS
      integer wi(nmax)
      double precision t(nmax),wd(nmax),z(nmax)

      delta2 = delta**2

      msinfo  = 0

      memfail = .false.

      call lssini(.false.,.false.,.true.)

C     step 1: initialize ls (lower bound on l) with max{-bii}, where bii 
C     are the elements of the diagonal of B
      
      ls = -bval(bdiag(1))
      
      do i = 2,n
          lin = bdiag(i)
          ls  = max( ls, -bval(lin) )
      end do
      
C     Calculate ||B||1, B sparse
      
      do i = 1,n
          wd(i) = 0.0d0
      end do     
      
      do i = 1,bnnz
          lin = blin(i)
          col = bcol(i)
          
          if ( ( lin .le. n ) .and. ( col .le. n ) ) then
              wd(col) = wd(col) + abs( bval(i) )
              if ( lin .ne. col ) then
                  wd(lin) = wd(lin) + abs( bval(i) )
              end if
          end if
      end do
      
      b1n = wd(1)
      do i = 2,n
          b1n = max( b1n, wd(i) )
      end do
      
C     step 2: initialize ll (lower bound on l) with 
C             max{0, ls, ||g||/delta - ||B||1}, where ||B||1 is the 
C             1-norm of the matrix B

      geucn = 0.0d0
      do i = 1,n
          geucn = geucn + g(i)**2
      end do
      geucn = sqrt( geucn )
      
      ll = (geucn / delta) - b1n 
      ll = max( 0.0d0, ll )
      ll = max( ls, ll )

C     step 3: initialize lu (upper bound on l) with ||g||/delta + ||B||1

      lu = (geucn / delta) + b1n

C     If the matrix is null, there is nothing to be done

      if ( ( abs( ll ) .le. macheps23 ) .and.
     +     ( abs( lu ) .le. macheps23 ) .and. 
     +     ( abs( ls ) .le. macheps23 ) ) then
          msinfo = 0
          go to 21
      end if
      
C     step 4: initialize iteration counter

      iter = 1

      call lssana(n,bnnz,blin,bcol,bval,bdiag,lssinfo)
      
      if ( lssinfo .eq. 6 ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          memfail = .true.
          return

      else if ( lssinfo .eq. 7 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          memfail = .true.
          return

      else if ( lssinfo .eq. 8 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          memfail = .true.
          return

      end if

C     step 5: safeguard of l (ensures that l is bigger than ll)
      
 5    continue

      l = max( l, ll )

C     step 6: safeguard of l (ensures that l is smaller than lu)

      l = min( l, lu )

C     step 7: safeguard of l

      if ( ( l .le. ls + macheps23 * max( abs( ls ), 1.0d0 ) )
     +     .and. ( iter .ne. 1 ) ) then
          l = max( 1.0d-3 * lu, sqrt( ll*lu ) )
      end if

C     step 8: try to use the Cholesky decomposition: (B +lI) = R^TR.
C             If the decomposition is successfull, R is stored in the 
C             upper triangular portion of B (including the diagonal) and 
C             pd is set to true.
C             If the decomposition fails, d and idx are set as explained 
C             before, pd is set to false, and the Euclidian-norm of u is 
C             calculated (see explanation of variable ueucn2)

      do i = 1,n
          wd(i) = l
      end do

      call lssfac(n,bnnz,blin,bcol,bval,bdiag,wd,idx,d,dum,lssinfo)
      chcnt = chcnt + 1
      
      if ( lssinfo .eq. 6 ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          memfail = .true.
          return

      else if ( lssinfo .eq. 7 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          memfail = .true.
          return

      else if ( lssinfo .eq. 8 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          memfail = .true.
          return

      end if

      if ( ( lssinfo .eq. 1 ) .or. ( lssinfo .eq. 2 ) ) then
          pd = .false.
      else ! lssinfo .eq. 0
          pd = .true.
      end if

C     In this case (B + lI) is not positive definite, and d and idx are 
C     calculated. Because p cannot be calculated (it is not possible to 
C     solve the system using Cholesky factorization), the values of l, 
C     ll and ls are updated, the iteration counter is increased and a
C     new iteration is started

      if ( .not. pd ) then

C         Print information (current iteration)

          if ( iprintinn .ge. 5 ) then
              write(*, 1000) iter
              write(*, 1010) ls,ll,lu,l
              write(*, 1070)
              
              write(file10_unit,1000) iter
              write(file10_unit,1010) ls,ll,lu,l
              write(file10_unit,1070) 
          end if
          
          llant = ll
          luant = lu
          lsant = ls
          
          if ( lu-l .le. macheps23 * max( abs( lu ),1.0d0 ) ) then
              lu = lu + macheps23 * max( abs( lu ),1.0d0 )
              l  = lu
          else
              call scalcu(n,bnnz,blin,bcol,bval,bdiag,l,idx,p,ueucn2,wd,
     +        memfail)
              chcnt = chcnt + 1
              
              if ( memfail ) then
                  go to 22
              end if
              
              ll = max( l, ll )
              ls = max( l + (d / ueucn2), ls )
              ll = max( ll, ls )
              l  = ls
          end if
          
          iter = iter + 1
          
C         Test whether the number of iterations is exhausted

          if ( iter .gt. maxit ) then
              msinfo = 5
              
              if ( iprintinn .ge. 5 ) then
                  write(*, 1090) 
                  write(file10_unit,1090) 
              end if
              
              go to 22
          end if
          
          go to 5
      end if
      
C     step 9: solve R^TRp = -g for p and calculate the squared 
C             Euclidian-norm of p
      
      do i = 1,n
          p(i) = -g(i)
      end do
      
      call lsssol(n,p)

C     Euclidian-norm of Rp = p^T R^TRp = p^T (-g) = - p^Tg
      
      rpeucn2 = 0.0d0
      do i = 1,n
          rpeucn2 = rpeucn2 - p(i) * g(i)
      end do
      
      peucn2 = 0.0d0
      do i = 1,n
          peucn2 = peucn2 + p(i)**2
      end do
      peucn = sqrt( peucn2 )
      
C     step 10: calculate z and tau, where tau * z is the approximation 
C              of the eigenvector associated with the smallest 
C              eigenvalue of B

      if ( peucn .lt. delta ) then

C        Calculate z

          call scalcz(n,wi,wd,z)
          
C         Calculate z Euclidian-norm

          tmp = 0.0d0
          do i = 1,n
              tmp = tmp + z(i)**2
          end do
          tmp = sqrt( tmp )
          
C         Divide z by its norm

          do i = 1,n
              z(i) = z(i) / tmp
          end do
          
C         Calculate the squared Euclidian-norm of the product Rz.         
C         Note that z^T R^T Rz = z^T (B + lI) z
          
          do i = 1,n
              wd(i) = 0.0d0
          end do
          
          do i = 1,bnnz
              lin = blin(i)
              col = bcol(i)
              
              if ( lin .eq. col ) then
                  wd(lin) = wd(lin) + (bval(i) + l) * z(col)
              else
                  wd(lin) = wd(lin) + bval(i) * z(col)
                  wd(col) = wd(col) + bval(i) * z(lin)
              end if
          end do
          
          rzeucn2 = 0.0d0
          do i = 1,n
              rzeucn2 = rzeucn2 + z(i) * wd(i)
          end do
          
C         Calculate tau

          ptz = 0.0d0
          do i = 1,n
              ptz = ptz + p(i) * z(i)
          end do
          
          if ( ptz .lt. 0.0d0 ) then
              sgn = -1.0d0
          else
              sgn =  1.0d0
          end if
          
          tmp = delta2 - peucn2
          
          tau = (ptz**2) + tmp
          tau = tmp / (ptz + sgn * sqrt( tau ))
          
      end if

C     Print informations (current iteration)

      if ( iprintinn .ge. 5 ) then
          write(*, 1000) iter
          write(*, 1010) ls,ll,lu,l
          write(*, 1020) peucn
          
          write(file10_unit,1000) iter
          write(file10_unit,1010) ls,ll,lu,l
          write(file10_unit,1020) peucn
          
          if ( peucn .lt. delta ) then
              write(*, 1030) sqrt( peucn2 + 2 * tau * ptz + tau ** 2 )
              write(file10_unit,1030) sqrt( peucn2 + 2 * tau * ptz + tau ** 2 )
          else
              write(*, 1030) peucn
              write(file10_unit,1030) peucn
          end if
          
          write(*, 1050) delta
          write(file10_unit,1050) delta
      end if
      
C     steps 11 and 12: update ll, lu and ls

      llant = ll
      luant = lu
      lsant = ls
      
      if ( peucn .lt. delta ) then
          lu = min( l, lu )
          ls = max( l - rzeucn2, ls )
      else
          ll = max( l, ll )
      end if

C     step 13: update ls when B + lI is not positive definite.
C              This was done right after the Cholesky decomposition 
C              failure

C     step 14: update ll

      ll = max( ll, ls )
      
C     step 15: convergence test

      if ( ( abs( l ) .le. eps ) .and. ( peucn .le. delta ) ) then
          msinfo = 1
          go to 21
      end if

C     step 16: second convergence test

      if ( abs( delta - peucn ) .le. sigma1*delta ) then
          msinfo = 2
      end if

C     step 17: convergence test for the hard case

      tmp = rpeucn2 + l * delta2
      tmp = max( sigma2, tmp )
      tmp = tmp * sigma1 * (2 - sigma1)
      
      if ( ( peucn .lt. delta ) .and. ( (rzeucn2*(tau**2) .le. tmp ) 
     +     .or. (lu - ll .le. eps) ) ) then     
          
          msinfo = msinfo + 3
          go to 20
      end if
      
      if ( msinfo .eq. 2 ) then
          go to 21
      end if
      
C     step 21: Calculate l to be used in the next iteration

      if ( ( abs( geucn ) .gt. eps ) .or. 
     +     ( ( l .le. ll + macheps23 * max( abs( ll ), 1.0d0 ) )
     +     .and. ( ls .le. ll ) ) ) then 
          
C         Solve R^T t = p for t and calculate t squared Euclidian-norm
         
          do i = 1,n
              t(i) = p(i)
          end do
          
          call lsssoltr('T',n,t)
          
          teucn2 = 0.0d0
          do i = 1,n
              teucn2 = teucn2 + t(i)**2
          end do
          
C         Update l using Newton's method update
        
          l = l + (peucn2 / teucn2) * ((peucn - delta) / delta)
      else
          l = ls
      end if
      
C     step 22: update iteration counter

      iter = iter + 1

C     Test whether the number of iterations is exhausted

      if ( iter .gt. maxit ) then
          msinfo = 5
          
          if ( iprintinn .ge. 5 ) then
              write(*, 1090) 
              write(file10_unit,1090) 
          end if
          
          go to 22
      end if      
      
C     step 23: start a new iteration

      if ( ( abs( llant-ll ) .le. 
     +     macheps23 * max( abs( ll ), 1.0d0 ) ) .and. 
     +     ( abs( luant-lu ) .le.
     +     macheps23 * max( abs( lu ), 1.0d0 ) ) .and. 
     +     ( abs( lsant-ls ) .le. 
     +     macheps23 * max( abs( ls ), 1.0d0 ) ) ) then
          
          ll = ll + macheps23 * max( abs( ll ), 1.0d0 )
          lu = lu - macheps23 * max( abs( lu ), 1.0d0 )
      end if
      
      go to 5

C     steps 18, 19 and 20:

C     The solution is given by p in 3 cases:
C     - if the first convergence criterion is satisfied;
C     - if only the second convergence criterion is satisfied;
C     - if both the second and the third convergence criteria are 
C       satisfied, but the squared Euclidian-norm of R(tau*z) is 
C       strictly bigger than l*(delta2 - peucn2)

C     The solution is given by p + tau*z when:
C     - just the third convergence criterion is satisfied;
C     - both the second and the third convergence criteria are 
C       satisfied, but the squared Euclidian-norm of R(tau*z) is smaller 
C       or equal to l*(delta2 - peucn2)

 20   tmp = (rzeucn2 * (tau**2)) - l * (delta2 - peucn2)
      
      if ( ( msinfo .eq. 3 ) .or. 
     +     ( ( msinfo .eq. 5 ) .and. ( tmp .le. 0.0d0 ) ) ) then
          
          peucn2 = 0.0d0
          do i = 1,n
              p(i) = p(i) + tau * z(i)
              peucn2 = peucn2 + p(i)**2
          end do
          peucn = sqrt( peucn2 )
          
          msinfo = 3
      else
          msinfo = 2
      end if
      
C     Print informations

 21   if ( iprintinn .ge. 5 ) then
          write(*, 1060) iter
          write(*, 1010) ls,ll,lu,l
          write(*, 1080) 
          write(*, 1040) peucn
          
          write(file10_unit,1060) iter
          write(file10_unit,1010) ls,ll,lu,l
          write(file10_unit,1080) 
          write(file10_unit,1040) peucn
      end if
      
 22   continue
      
C     Non-executable statements

 1000 format(/,10X,'More-Sorensen iteration: ',I7)
 1010 format(  10X,'ls = ',1P,D11.4,
     +         10X,'ll = ',1P,D11.4,
     +         10X,'lu = ',1P,D11.4,
     +         10X,'l  = ',1P,D11.4)
 1020 format(  10X,'Euclidian-norm of p: ',1P,D7.1)
 1030 format(  10X,'Euclidian-norm of p + z: ',1P,D7.1)
 1040 format(  10X,'Euclidian-norm of step: ',1P,D7.1)
 1050 format(  10X,'delta: ',1P,D7.1)
 1060 format(  10X,'Number of iterations: ',I7)
 1070 format(  10X,'Matrix is not positive definite!')

 1080 format(  10X,'Flag of More-Sorensen: ',
     +             'convergence criterion satisfied.')
 1090 format(  10X,'Flag of More-Sorensen: ',
     +             'maximum number of iterations achieved.')

      end

C     ******************************************************************
C     ******************************************************************
      
      subroutine scalcu(n,annz,alin,acol,aval,adiag,l,idx,u,ueucn2,wd,
     +memfail)
            
      implicit none

C     SCALAR ARGUMENTS
      logical memfail
      integer annz,idx,n
      double precision l,ueucn2
      
C     ARRAY ARGUMENTS
      integer acol(annz),adiag(n),alin(annz)
      double precision aval(annz),u(idx),wd(idx-1)
      
C     Solves a sparse linear system of the form (A + lI + ekek^Td)u = 0,
C     where A + lI is a definite positive matrix in R^{k x k}.
C     u is a vector of k positions with u(k) = 1.

C     LOCAL SCALARS
      integer col,i,lin,lssinfo
      
      if ( idx .eq. 1 ) then
          u(idx) = 1.0d0
          ueucn2 = 1.0d0
          return
      end if
      
      do i = 1,idx
          u(i) = 0.0d0
      end do
      
C     Permute columns and rows of A

      call lsspermind(annz,alin,acol)
      
C     Eliminate columns that have index greater than idx and define u
C     as idx-th column of A + lI
      
      do i = 1,annz
          col = acol(i)
          lin = alin(i)
          
          if ( ( col .eq. idx ) .and. ( lin .lt. idx ) ) then
              u(lin) = u(lin) - aval(i)
          else if ( ( lin .eq. idx ) .and. ( col .lt. idx ) ) then
              u(col) = u(col) - aval(i)
          end if
          
          if ( ( col .ge. idx ) .or. ( lin .ge. idx ) ) then
              acol(i) = col + n
              alin(i) = lin + n
          end if
      end do
      
C     Solve system (A + lI)x = u

      call lssini(.false.,.true.,.false.)
      
      call lsspermvec(n,adiag)
      
      do i = 1,idx-1
          wd(i) = l
      end do
      
      call lssafsol(idx-1,annz,alin,acol,aval,adiag,wd,u,lssinfo)
      
      if ( lssinfo .eq. 6 ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          memfail = .true.
          return

      else if ( lssinfo .eq. 7 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          memfail = .true.
          return

      else if ( lssinfo .eq. 8 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          memfail = .true.
          return

      end if
      
      call lssunpermvec(n,adiag)
      
      call lssini(.false.,.false.,.true.)
      
C     Undo modifications in acol and alin
      
      do i = 1,annz
          col = acol(i)
          
          if ( col .gt. n ) then
              acol(i) = acol(i) - n
              alin(i) = alin(i) - n
          end if
      end do
      
      call lssunpermind(n,annz,alin,acol)
      
      u(idx) = 1.0d0
      
      ueucn2 = 0.0d0
      do i = 1,idx
          ueucn2 = ueucn2 + u(i)**2
      end do
      
      end
      
C     ******************************************************************
C     ******************************************************************      
      
      subroutine scalcz(n,rowind,rowval,z)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      integer rowind(n)
      double precision rowval(n),z(n)
      
C     Sparse implementation of the technique presented by Cline, Moler,
C     Stewart and Wilkinson to estimate the condition number of a upper
C     triangular matrix.

#include "machconst.com"

C     LOCAL SCALARS
      integer col,i,k,rownnz
      double precision acs,acz,d,ek,rki,s,sm,w,wk,wkm
      double precision lssgetd

C     EXTERNAL FUNCTIONS
      external lssgetd

      call lsssetrow(n)
      
      do i = 1,n
          z(i) = 0.0d0
      end do
      
      acz = 0.0d0
      acs = 1.0d0
      ek  = 1.0d0
      
      do k = 1,n
          d = lssgetd(k)
          
          if ( abs( acs*z(k) ) .gt. macheps23 ) then
              ek = dsign(ek,-acs*z(k))
          end if
          
          if ( abs( ek - acs*z(k) ) .gt. abs( d ) ) then
              s   = abs( d ) / abs( ek - acs*z(k) ) 
              acs = acs * s
              ek  = s * ek
          end if
          
          wk  =  ek - acs * z(k)
          wkm = -ek - acs * z(k)
          s   = abs( wk )
          sm  = abs( wkm )
          
          if ( d .eq. 0.0d0 ) then 
              wk  = 1.0d0
              wkm = 1.0d0
          else
              wk  = wk / d
              wkm = wkm / d
          end if
          
          if ( k .eq. n ) then
              go to 10
          end if
          
          call lssgetrow(n,k,rownnz,rowind,rowval)
          
          sm = sm + acs * acz
          
          do i = 1,rownnz
              
              col = rowind(i)
              if ( col .gt. k ) then
                  rki    = rowval(i)
                  sm     = sm - abs( acs * z(col) )
                  sm     = sm + abs( acs * z(col) + wkm * rki )
                  acz    = acz - abs( z(col) )
                  acz    = acz + abs( z(col) + (wk/acs) * rki )
                  z(col) = z(col) + (wk/acs) * rki
              end if
          end do
          
          s = s + acs * acz
          
          if ( s .lt. 
     +         sm - macheps23 * max( abs( sm ), 1.0d0 ) ) then
              w = wkm - wk
              wk = wkm
              
              do i = 1,rownnz
                  
                  col = rowind(i)
                  if ( col .gt. k ) then
                      rki    = rowval(i)
                      acz    = acz - abs( z(col) )
                      acz    = acz + abs( z(col) + (w/acs) * rki )
                      z(col) = z(col) + (w/acs) * rki
                  end if
              end do
              
          end if
          acz = acz - abs( z(k+1) )
 10       z(k) = wk/acs
      end do
      
C     Divide z by its 1-norm to avoid overflow
      
      s = 0.0d0
      do i = 1,n
          s = s + abs( z(i) )
      end do
      
      do i = 1,n
          z(i) = z(i) / s
      end do
      
C     Solve Rz = y

      call lsssoltr(' ',n,z)

      end

C     ******************************************************************
C     ******************************************************************

C     moresor:
C
C     Method to minimize sparse quadratic functions subjected to 
C     ||w|| <= delta
C
C     minimize    psi(w) = 1/2 w^TBw + g^Tw 
C     subject to  ||w|| <= delta
C
C     Method described in "Computing a trust region step", by More and 
C     Sorensen.
C
C     The main ideia of this method is to find a positive scalar \mslamb 
C     that is a zero of the function
C
C     phi(\mslamb) = 1/||p|| - 1/delta,
C      
C     where p is the solution of the linear system
C
C     (B + \mslamb I)p = -g.
C
C     Note that symmetric matrix B, vector g and positive real number 
C     delta are presented in the minimization problem above. I is the 
C     identity matrix.
C
C     The method used to find the zero of that function is basically the 
C     Newton method to find roots.
C
C     On Entry
C
C     n        integer
C              dimension
C
C     g        double precision g(n)
C              vector used to define the quadratic function
C
C     bnnz     integer
C              number of nonzero elements of B
C
C     blin     integer blin(bnnz)
C              row indices of nonzero elements of B
C
C     bcol     integer bcol(bnnz)
C              column indices of nonzero elements of B
C
C     bval     double precision bval(bnnz)
C              nonzero elements of B, that is,
C              B(blin(i),bcol(i)) = bval(i). Since B is symmetric, just
C              one element B(i,j) or B(j,i) must appear in its sparse
C              representation. If more than one pair corresponding to
C              the same position of B appears in the sparse
C              representation, the multiple entries will be summed.
C              If any blin(i) or bcol(i) is out of range, the entry will
C              be ignored
C     
C     bdiag    integer bdiag(n)
C              indices of diagonal elements of B in blin, bcol and bval
C
C     delta    double precision
C              trust-region radius
C     
C     sigma1   double precision
C              allowed error for convergence criteria 1 and 2
C
C     sigma2   double precision
C              allowed error for convergence criteria 3 (hard case)
C
C     eps      double precision
C              allowed error
C
C     maxit    integer
C              maximum number of allowed iterations
C     
C     l        double precision
C              initial value for mslamb
C
C     On Return
C      
C     l        double precision
C              value that gives p as a solution to the minimization 
C              problem, because p is also solution to 
C              (B + l I)p = -g
C
C     pd       logical
C              set to true if the last Cholesky decomposition is 
C              successfull
C
C     p        double precision p(n)
C              solution to problem
C              minimize     psi(w)
C              subjected to ||w|| <= delta
C
C     chcnt    integer
C              number of Cholesky decompositions
C
C     memfail  logical
C              true iff linear solver failed because of lack of memory
C
C     msinfo   integer
C              stores which convergence criteria was satisfied:
C
C              0 = both g and B are null;
C
C              1 = first convergence criterion is satisfied;
C
C              2 = second convergence criterion is satisfied;
C
C              3 = third convergence criterion is satisfied
C
C              5 = maximum allowed number of iterations is achieved.

C     ******************************************************************
C     ******************************************************************

C     scalcu:
C
C     Solve a sparse linear system of the form (A + lI + ekek^Td)u = 0,
C     where A + lI is a definite positive matrix in R^{k x k}.
C     u is a vector of k positions with u(k) = 1.
C
C     On Entry
C      
C     n        integer
C              dimension of A
C
C     annz     integer
C              number of nonzero elements of A
C      
C     alin     integer alin(annz)
C              row indices of nonzero elements of A
C
C     acol     integer acol(annz)
C              column indices of nonzero elements of A
C
C     aval     integer aval(annz)
C              nonzero elements of A, that is,
C              A(alin(i),acol(i)) = aval(i). Since A is symmetric, just
C              one element A(i,j) or A(j,i) must appear in its sparse
C              representation. If more than one pair corresponding to
C              the same position of A appears in the sparse
C              representation, the multiple entries will be summed.
C              If any alin(i) or acol(i) is out of range, the entry will
C              be ignored
C     
C     adiag    integer adiag(n)
C              indices of diagonal elements of A in alin, acol and aval
C
C     l        double precision
C              used to compute A + lI
C      
C     idx      integer
C              index k 
C
C     On Return
C      
C     u        double precision u(n)
C              system solution
C
C     ueucn2   double precision
C              u squared Euclidian-norm
C
C     memfail  logical
C              true iff linear solver failed because of lack of memory

C     ******************************************************************
C     ******************************************************************      
      
C     scalcz:
C
C     Sparse implementation of the technique presented by Cline, Moler,
C     Stewart and Wilkinson to estimate the condition number of a matrix.
C     This technique is used by the More-Sorensen method to calculate an
C     approximation to the eigenvector associated to the smallest 
C     eigenvalue of a matrix B (\mslamb_1). 
C     In this technique, when \mslamb approaches -\mslamb_1, \|Rz\| 
C     approaches 0. This insures that z is an approximation to the 
C     wanted eigenvector.
C     Basically, it solves R^Ty = e, choosing e(k) as 1 or -1 (whatever 
C     gives maximum local growth of y). Then, it solves Rz = y.
C     Note that R is a sparse matrix given by P^T D^0.5 L^T P (obtained
C     applying subroutine MA27BD).
C      
C     On Entry
C
C     n        integer      
C              dimension of R
C      
C     rowind   integer rowind(n)
C     rowval   double precision rowval(n)
C              working arrays
C      
C     On Return
C      
C     z        double precision z(n)
C              approximation of the eigenvector of B associated to
C              \mslamb_1
