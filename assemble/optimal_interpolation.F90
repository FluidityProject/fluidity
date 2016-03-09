!! Copyright (C) 2005 Alexander Barth <barth.alexander@gmail.com>
!! Copyright (C) 2006 David Saunders (NASA Ames Research Center)
!! Copyright (C) 2008 Gian Franco Marras (CINECA) <g.marras@cineca.it>
!!
!! This program is free software; you can redistribute it and/or modify it under
!! the terms of the GNU General Public License as published by the Free Software
!! Foundation; either version 3 of the License, or (at your option) any later
!! version.
!!
!! This program is distributed in the hope that it will be useful, but WITHOUT
!! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
!! details.
!!
!! You should have received a copy of the GNU General Public License along with
!! this program; if not, see <http://www.gnu.org/licenses/>.

!! Fortran 90 module for n-dimensional optimal interpolation.
!! Dependencies: LAPACK (dsyev)

#define DIAG_OBS_COVAR
      module optimal_interpolation
!     working precision 
!     4 = simple precision, 8 is double precision

      integer, parameter :: wp = 8


      contains

!     --------------------------------------------------------------------------

      subroutine select_nearest(x,ox,param,m,index,distance)

!     Select the m observations from ox(1:nd,1:on) closest to point x(1:nd).

!     Arguments:

      implicit none
      real(wp),intent(in)  :: x(:),ox(:,:),param(:)
      integer, intent(in)  :: m
      real(wp),intent(out) :: distance(m)
      integer, intent(out) :: index(m)

!     Local variables:

      real(wp) :: d(size(ox,2))
      integer  :: i

!     Execution:

!     Calculate a measure of (squared) distance to each observation:

      do i=1,size(ox,2)
        d(i) = sum(((x - ox(:,i)) * param)**2)
      end do

      call sort(d,m,index)

      distance = d(index)

      end subroutine select_nearest


!     --------------------------------------------------------------------------

      subroutine sort(d,m,pannier)

!     Return the indices of the m smallest elements in d(:).
!     The algorithm is succinctly coded, but would a heap sort be faster?

!     Arguments:

      implicit none
      real(wp), intent(in)  :: d(:)
      integer,  intent(in)  :: m
      integer,  intent(out) :: pannier(m)

      integer :: i,max_pannier(1)

      do i=1,m
        pannier(i) = i
      end do

      max_pannier = maxloc(d(pannier))

      do i=m+1,size(d)
        if (d(i) .lt. d(pannier(max_pannier(1)))) then
          pannier(max_pannier(1)) = i
          max_pannier = maxloc(d(pannier))
        end if
      end do

      end subroutine sort

!     --------------------------------------------------------------------------

      subroutine observation_covariance(ovar,index,R)
      implicit none
      real(wp),intent(in) ::ovar(:)
      integer,intent(in) :: index(:)
#     ifdef DIAG_OBS_COVAR
      real(wp),    intent (out) :: R(size(index))
#     else
      real(wp),    intent (out) :: R(size(index),size(index))
      integer :: i
#     endif


#     ifdef DIAG_OBS_COVAR
      R = ovar(index)
#     else
      R = 0
      do i=1,size(index)
        R(i,i) = ovar(index(i))
      end do
#     endif

      end subroutine observation_covariance

!     --------------------------------------------------------------------------
!     Modified Bessel functions of 2nd kind (order 1)

      function mod_bessel_K1(x)       
      implicit none
      real(wp), intent(in) :: x
      real(wp) :: mod_bessel_K1

      real(8) :: y
      real(8), parameter :: & 
           P1 = 1.0D0, &
           P2 = 0.15443144D0, &
           P3 = -0.67278579D0, &
           P4 = -0.18156897D0, &
           P5 = -0.1919402D-1, &
           P6 = -0.110404D-2, &
           P7 = -0.4686D-4, &
           Q1 = 1.25331414D0, &
           Q2 = 0.23498619D0, &
           Q3 = -0.3655620D-1, &
           Q4 = 0.1504268D-1, &
           Q5 = -0.780353D-2, &
           Q6 = 0.325614D-2, &
           Q7 = -0.68245D-3

      if(x <= 0.) stop 'error x <= 0' 

      if(x <= 2.0) then
         y = x * x * 0.25
         mod_bessel_K1 = (log(x/2.0)*mod_bessel_I1(x))+(1.0/x)*(P1+y*(P2+y*(P3+ &
              y*(P4+y*(P5+y*(P6+y*P7))))))
      else
         y = 2.0 / x
         mod_bessel_K1 = (exp(-x)/sqrt(x))*(Q1+y*(Q2+y*(Q3+   &
              y*(Q4+y*(Q5+y*(Q6+y*Q7))))))
      endif

      end function

!     --------------------------------------------------------------------------
!     Modified Bessel functions of 1st kind (order 1)

      function mod_bessel_I1(x)
      implicit none
      real(wp) :: x
      real(wp) :: mod_bessel_I1

      real(8) :: y, ax
      real(8), parameter :: & 
           P1 = 0.5D0, &
           P2 = 0.87890594D0, &
           P3 = 0.51498869D0, & 
           P4 = 0.15084934D0, &
           P5 = 0.2658733D-1, &
           P6 = 0.301532D-2, &
           P7 = 0.32411D-3, &
           Q1 = 0.39894228D0, &
           Q2 = -0.3988024D-1, &
           Q3 = -0.362018D-2, &
           Q4 = 0.163801D-2, &
           Q5 = -0.1031555D-1, &
           Q6 = 0.2282967D-1, &
           Q7 = -0.2895312D-1, &
           Q8 = 0.1787654D-1, &
           Q9 = -0.420059D-2

      if(abs(x) < 3.75) then
         y = x*x / (3.75*3.75)
         mod_bessel_I1 = x*(P1+y*(P2+y*(P3+y*(P4+y*(P5+y*(P6+y*P7))))))
      else
         ax = abs(x)
         y = 3.75 / ax
         mod_bessel_I1 = (exp(ax)/sqrt(ax))*(Q1+y*(Q2+y*(Q3+         &
              y*(Q4+y*(Q5+y*(Q6+y*(Q7+y*(Q8+y*Q9))))))))
         if (x < 0.) mod_bessel_I1 = - mod_bessel_I1
      endif

      end function

!     --------------------------------------------------------------------------

      function  background_covariance_diva(x1,x2,param) result(c)
      implicit none
      real(wp),intent(in) :: x1(:),x2(:),param(:)
      real(wp) :: c

      real(wp) :: d(size(x1))
      real(wp) :: rn 

      d = (x1 - x2)*param
      rn = sqrt(sum(d**2))

      if (rn == 0) then
        c = 1
      else
        c = rn * mod_bessel_K1(rn)
      end if

      end function background_covariance_diva


!     --------------------------------------------------------------------------

      function  background_covariance_gaussian(x1,x2,param) result(c)
      implicit none
      real(wp),intent(in) :: x1(:),x2(:),param(:)
      real(wp) :: c

      real(wp) :: d(size(x1))

      d = (x1 - x2)*param
      c = exp(-sum(d**2))

      end function background_covariance_gaussian

!     --------------------------------------------------------------------------


      function  background_covariance(x1,x2,param) result(c)
      implicit none
      real(wp),intent(in) :: x1(:),x2(:),param(:)
      real(wp) :: c

#     if defined(BACKGROUND_COVARIANCE_DIVA)
      c = background_covariance_diva(x1,x2,param)
#     else
      c = background_covariance_gaussian(x1,x2,param)
#     endif

      end function background_covariance




!     --------------------------------------------------------------------------

      subroutine pinv (A, tolerance, work, D)

!     Compute pseudo-inverse A+ of symmetric A in factored form U D+ U', where
!     U overwrites A and D is diagonal matrix D+.

!     Saunders notes: Working with the factors of the pseudo-inverse is
!                     preferable to multiplying them together (more stable,
!                     less arithmetic).  Also, the choice of tolerance is
!                     not straightforward.  If A is noisy, try 1.e-2 * || A ||,
!                     else machine-eps. * || A ||  (1 norms).
!     Arguments:

      real(wp), intent (inout) :: A(:,:)  ! Upper triangle input; orthogonal U out
      real(wp), intent (in)    :: tolerance
      real(wp), intent (out)   :: work(:), D(:)

!     Local variables:

      integer :: i, info, N

!     Execution:

      N = size (A,1)

!     Eigendecomposition/SVD of symmetric A:

      call dsyev ('V', 'U', N, A, N, D, work, size (work), info)

!     Diagonal factor D+ of pseudo-inverse:

      do i = 1, N
         if (D(i) > tolerance) then
            D(i) = 1. / D(i)
         else
            D(i) = 0.
         end if
      end do


      end subroutine pinv

!     --------------------------------------------------------------------------

      function pinv_workspace(N) result(lwork)

!     Determine the workspace needed for dsyev.

      implicit none
      integer,intent(in) :: N
      integer :: lwork

      integer :: info
      real(wp) :: dummy,rwork

      call dsyev('V','U', N, dummy,N, dummy, rwork, -1, info )
      lwork = ceiling(rwork)

      end function

!     --------------------------------------------------------------------------

      subroutine optiminterp(ox,of,ovar,param,m,gx,gf,gvar)

!     Main optimal interpolation routine

      implicit none
      real(wp), intent(in)  :: ox(:,:), of(:,:)  ! Observations
      real(wp), intent(in)  :: ovar(:)           ! Observation error variances
      real(wp), intent(in)  :: param(:)          ! inverse of correlation lengths
      integer,  intent(in)  :: m                 ! # nearest observations used
      real(wp), intent(in)  :: gx(:,:)           ! Target grid coords.
      real(wp), intent(out) :: gf(:,:), gvar(:)  ! Interpolated fields.
                                                 ! and error variances
!     Local variables:

      real(wp) :: PH(m), PHiA(m),A(m,m),D(m)

#ifdef DIAG_OBS_COVAR
      real(wp) :: R(m)
#else
      real(wp) :: R(m,m)
#endif

      integer  :: gn,nf,index(m)
      real(wp) :: distance(m)

      integer  :: i,j1,j2,j,lwork
      real(wp) :: tolerance = 1e-5

#     ifdef VERBOSE
      integer  :: percentage_done
#     endif
#     ifdef STATIC_WORKSPACE
      real(wp) :: work((m+2)*m)
#     else
      real(wp), allocatable :: work(:)
#     endif

!     Execution:
      gn = size(gx,2)
      nf = size (of, 1)  ! # functions at each observation point

#     ifndef STATIC_WORKSPACE
!     query and allocate workspace for pseudo-inverse
      lwork = pinv_workspace(m)
#     endif

!$omp parallel  default(none) private(work,i,PH,PHiA,index,distance,j1,j2,R,D,A) &
!$omp shared(gf,gn,gvar,gx,lwork,m,of,ovar,ox,param,tolerance)

#     ifndef STATIC_WORKSPACE
      allocate(work(lwork))
#     endif
      
#     ifdef VERBOSE
      percentage_done = 0
#     endif

!$omp do 
      do i=1,gn
print*, '1111111111111111111111',gn,i
#       ifdef VERBOSE
        if (percentage_done .ne. int(100.*real(i)/real(gn))) then
          percentage_done = int(100.*real(i)/real(gn))
          write(6,*) 'done ',percentage_done
        end if
#       endif
print*, '22222222222222222222222222222'

!       get the indexes of the nearest observations
        call select_nearest(gx(:,i),ox,param,m,index,distance)
     
!       form compute the error covariance matrix of the observation 
print*, '33333333333333333333333333333333'

        call observation_covariance(ovar,index,R)
print*, '444444444444444444444444444444444444'

!       form the error covariance matrix background field

        do j2 = 1, m
          ! Upper triangle only

           do j1 = 1, j2  
             A(j1,j2) = &
                 background_covariance(ox(:,index(j1)),ox(:,index(j2)),param)
           end do

           PH(j2) = background_covariance(gx(:,i),ox(:,index(j2)),param)
        end do

print*, '555555555555555555555555555555555555'
!       covariance matrix of the innovation

#ifdef DIAG_OBS_COVAR
        do j = 1, m
           A(j,j) = A(j,j) + R(j)
        end do
#else
        A  = A + R
#endif

!       pseudo inverse of the covariance matrix of the innovation

        call pinv(A,tolerance,work,D)

print*, '66666666666666666666666666666666666'
        PHiA = matmul (A, D * matmul (PH, A))

!       compute the analysis for all fields

        gf(:,i) = matmul(of(:,index),PHiA)

!       compute the error variance of the analysis

        gvar(i) = 1. - dot_product(PHiA,PH)

print*, '777777777777777777777777777777777777'
      end do
!$omp end do

#     ifndef STATIC_WORKSPACE
      deallocate(work)
#     endif
 
!$omp end parallel 


      end subroutine optiminterp

#ifdef FORTRAN_2003_INTEROP

      subroutine optiminterp_gw(n,gn,on,nparam,ox,of,ovar,
     &                           param,m,gx,gf,gvar) bind(c)
      USE ISO_C_BINDING
      implicit none
      integer(C_INT) :: m,n,gn,on,nparam
      real(C_DOUBLE) :: gx(n,gn),ox(n,on),of(on),ovar(on),param(nparam)
      real(C_DOUBLE) :: gf(gn),gvar(gn)


      call optiminterp(ox,of,ovar,param,m,gx,gf,gvar)
      
      end subroutine optiminterp_gw


#endif

      end module optimal_interpolation
