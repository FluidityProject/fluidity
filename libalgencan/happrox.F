C     *****************************************************************
C     *****************************************************************

      subroutine comphapp(n,m,rho,equatn)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision rho(m)

C     This subroutine computes an approximation H of the Hessian of the
C     Augmented Lagrangian following a very simple idea: "discard the 
C     second order terms and then correct the remaining matrix in order 
C     to satisfy a secant equation".
C
C     Hence, H takes the form
C
C     H = B + S + rho A^t A, 
C
C     where S is the spectral correction of (rho A^t A) and B is 
C     the BFGS correction of (S + rho A^t A). More specifically,
C
C     S = hlspg I,
C
C     where
C
C     hlspg = max(lspgmi, min(lspgma, s^t (y - rho A^t A s) / s^t s)),
C
C     D = S + rho A^t A,
C
C     and
C
C     B = [ y y ^t / ( y^t s ) ] - [ D s ( D s )^t / ( s^t D s ) ].
C
C     Note that this subroutine does not compute matrix H explicitly,
C     but computes some quantities that will be used later to compute 
C     the product of H by a vector p.
C
C     The quantities computed by this subroutine are:
C
C     (a) hlspg = s^t (y - rho A^t A s) / (s^t s)
C
C     (b) hds = D s = ( hlspg I + rho A^t A ) s, and
C
C     (c) hstds = <s,hds>.

#include "dim.par"
#include "algconst.par"
#include "graddat.com"
#include "sydat.com"
#include "happdat.com"

C     LOCAL SCALARS
      integer i,j
      double precision ats

C     ------------------------------------------------------------------
C     Compute hds = rho A^t A s
C     ------------------------------------------------------------------

      do i = 1,n
          hds(i) = 0.0d0
      end do

      do j = 1,m

          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

C             COMPUTE THE INNER PRODUCT <a,s>
              ats = 0.0d0
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  ats = ats + jcval(i) * s(jcvar(i))
              end do

              ats = ats * rho(j)

C             ADD rho * ats * a
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  hds(jcvar(i)) = hds(jcvar(i)) + ats * jcval(i) 
              end do

          end if

      end do

      hstds = 0.0d0
      do i = 1,n
          hstds = hstds + s(i) * hds(i)
      end do

C     ------------------------------------------------------------------
C     Compute hlspg = s^t (y - rho A^t A s) / (s^t s)
C     ------------------------------------------------------------------

      if ( sty - hstds .le. 0.0d0 ) then
          hlspg = lspgmi
      else
          hlspg = max( lspgmi, min( (sty - hstds) / sts, lspgma ) )
      end if

      do i = 1,n
          hds(i) = hds(i) + hlspg * s(i)
      end do

      hstds = hstds + hlspg * sts

      end

C     *****************************************************************
C     *****************************************************************

      subroutine applyhapp(n,m,rho,equatn,goth,p,hp)

      implicit none

C     SCALAR ARGUMENTS
      logical goth
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision hp(n),p(n),rho(m)

C     This subroutine computes the product of the matrix computed by 
C     subroutine comphapp times vector p.

#include "dim.par"
#include "machconst.com"
#include "graddat.com"
#include "sydat.com"
#include "itetyp.com"
#include "happdat.com"

C     LOCAL SCALARS
      integer i,j
      double precision atp,c1,c2,ptds,pty

C     ------------------------------------------------------------------
C     Compute Hessian approximation
C     ------------------------------------------------------------------

      if ( .not. goth ) then

          goth = .true.
          call comphapp(n,m,rho,equatn)

      end if

C     ------------------------------------------------------------------
C     Compute ( hlspg I ) p
C     ------------------------------------------------------------------

      do i = 1,n
          hp(i) = hlspg * p(i)
      end do

C     ------------------------------------------------------------------
C     Add ( rho A^T A ) p
C     ------------------------------------------------------------------

      do j = 1,m

          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

C             COMPUTE THE INNER PRODUCT <a,p>
              atp = 0.0d0
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  atp = atp + jcval(i) * p(jcvar(i))
              end do

              atp = atp * rho(j)

C             ADD rho * atp * a
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  hp(jcvar(i)) = hp(jcvar(i)) + atp * jcval(i)
              end do

          end if

      end do

C     ------------------------------------------------------------------
C     Add B p, 
C     where B = [ y y ^t / ( y^t s ) ] - [ D s ( D s )^t / ( s^t D s ) ]
C     ------------------------------------------------------------------

      if ( sameface .and. sty .gt. macheps12 * seucn * yeucn ) then

          pty = 0.0d0
          ptds = 0.0d0
          do i = 1,n
              pty = pty + p(i) * y(i)
              ptds = ptds + p(i) * hds(i)
          end do

          c1 = pty / sty
          c2 = ptds / hstds
          do i = 1,n
              hp(i) = hp(i) + c1 * y(i) - c2 * hds(i)
          end do

      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine comphpre(n,m,rho,equatn)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision rho(m)

C     Consider the preconditioner 
C
C         P = Q + E + diag(rho A^t A) 
C
C     for matrix 
C
C         H = B + S + rho A^t A, 
C
C     where E is the spectral correction of diag(rho A^t A) and Q is 
C     the BFGS correction of (E + diag(rho A^t A)), while S and B are
C     the spectral and BFGS corrections of matrix (rho A^t A), 
C     respectively. 
C
C     This subroutine computes:
C
C     (a) pdiag = diag(rho A^t A),
C
C     (b) plspg such that E = plspg I, 
C
C     (c) psmdy = s - D^-1 y, where D = E + diag(rho A^t A), and
C
C     (d) the inner product psmdty = <psmdy,y>.
C
C     These quantities will be used later, in subroutine applyp, to
C     compute z = P^{-1} r.

#include "dim.par"
#include "machconst.com"
#include "algconst.par"
#include "graddat.com"
#include "sydat.com"
#include "itetyp.com"
#include "hpredat.com"

C     LOCAL SCALARS
      integer i,j
      double precision sttmp

C     ------------------------------------------------------------------
C     Compute diag( rho A^t A )
C     ------------------------------------------------------------------

      do i = 1,n
          pdiag(i) = 0.0d0
      end do

      do j = 1,m
          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  pdiag(jcvar(i)) = 
     +            pdiag(jcvar(i)) + rho(j) * jcval(i) ** 2
              end do
          end if
      end do

C     ------------------------------------------------------------------
C     Compute plspg = s^t (y - diag( rho A^t A ) s) / (s^t s)
C     ------------------------------------------------------------------

      sttmp = 0.0d0
      do i = 1,n
          sttmp = sttmp + pdiag(i) * s(i) ** 2
      end do

      if ( sty - sttmp .le. 0.0d0 ) then
          plspg = lspgmi
      else
          plspg = max( lspgmi, min( ( sty - sttmp ) / sts, lspgma ) )
      end if

C     ------------------------------------------------------------------
C     Compute the BFGS correction Q of ( E + diag( rho A^t A ) )
C
C     Q = [ (s - D^-1 y) s^t + s (s - D^-1 y)^t ] / s^t y -
C         [ <s - D^-1 y, y> s s^t ] / (s^t y)^2,
C
C     where D = ( E + diag( rho A^t A ) )
C     ------------------------------------------------------------------
 
      if ( sameface .and. sty .gt. macheps12 * seucn * yeucn ) then

          psmdyty = 0.0d0
          do i = 1,n
              psmdy(i) = s(i) - y(i) / ( plspg + pdiag(i) )
              psmdyty = psmdyty + psmdy(i) * y(i)
          end do

      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine applyhpre(n,m,rho,equatn,gotp,r,z)

      implicit none

C     SCALAR ARGUMENTS
      logical gotp
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision r(n),rho(m),z(n)

C     This subroutine computes the product of the inverse of the matrix
C     computed by subroutine comphpre times vector r, i.e., z = P^{-1} r.

#include "dim.par"
#include "machconst.com"
#include "sydat.com"
#include "itetyp.com"
#include "hpredat.com"

C     LOCAL SCALARS
      integer i
      double precision c1,c2,psmdytr,str

C     ------------------------------------------------------------------
C     Compute P
C     ------------------------------------------------------------------

      if ( .not. gotp ) then

          gotp = .true.
          call comphpre(n,m,rho,equatn)

      end if

C     ------------------------------------------------------------------
C     Compute ( E + diag( rho A^T A ) )^{-1} r
C     ------------------------------------------------------------------

      do i = 1,n
          z(i) = r(i) / ( plspg + pdiag(i) )
      end do

C     ------------------------------------------------------------------
C     Add Q^{-1} r, where
C
C     Q^{-1} = [ (s - D^-1 y) s^t + s (s - D^-1 y)^t ] / s^t y -
C              [ <s - D^-1 y, y> s s^t ] / (s^t y)^2
C
C     and D = ( E + diag( rho A^T A ) )
C     ------------------------------------------------------------------

      if ( sameface .and. sty .gt. macheps12 * seucn * yeucn ) then

          str = 0.0d0
          psmdytr = 0.0d0
          do i = 1,n
              str = str + s(i) * r(i)
              psmdytr = psmdytr + psmdy(i) * r(i)
          end do

          c1 = str / sty
          c2 = psmdytr / sty - psmdyty * str / sty ** 2

          do i = 1,n
              z(i) = z(i) + c1 * psmdy(i) + c2 * s(i) 
          end do

      end if

      end
