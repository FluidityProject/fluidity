
#include "fdebug.h"

 module smolyak
     
      use rbf_interp
      use spud 
      ！use fluids_module
   
     implicit none
     
      
    public ::  lagsum
    integer :: lagsum
    integer :: total_nd
    integer :: nsvd
contains

subroutine tensor_product_main (timestep,podnum, pod_coef, total_timestep , sparse_max )

 implicit none 
  integer, intent(in) :: total_timestep,podnum ,timestep,sparse_max
  real, intent(inout) :: pod_coef(podnum)  
  integer :: m,nd  
  real ( kind = 8 ) int_error  
  ! real ( kind = 8 ) x1d(n1d)  
  integer :: nsvd
  integer ::  i,k,j 
  real , dimension(:,:), allocatable ::  pod_coef_all_obv, pod_coef_all, pod_coef_all_lag_tmp   
  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:) 
  integer ( kind = 4 ), allocatable :: n_1d(:) 
  integer  ni
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:) !ZD(ND), the function evaluated at the points XD.
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)  
  real ( kind = 8 ), allocatable :: fe(:,:) ! exact values of function produced from fluidity
   real ( kind = 8 ), allocatable :: fi(:) 
  integer :: inteval,total_nd
  call get_option(&
            '/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
  m = podnum
  inteval=2
  nd= total_timestep/inteval  
   ni = 1 
  allocate(pod_coef_all_obv(total_timestep,podnum))
  allocate(pod_coef_all(total_timestep,podnum)) 
  
  
  allocate ( n_1d(1:m) ) 
  allocate ( a(1:m) )
  allocate ( b(1:m) )
 
  allocate ( ze(1:ni) )
  allocate ( zi(1:ni))
 
  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00
  
  !  Define the interpolation evaluation information.
  open(40,file='total', ACTION='read') 
  read(40,*) total_nd
  close(40)
  allocate ( zd(1:total_nd)) 
  allocate (fe(1:m,1:total_nd))
  allocate (xd(m,total_nd))
  allocate ( xi(1:m,1:ni) )
  allocate ( fi(1:ni) ) 
  open(41,file='coef_pod_all_fvalue', ACTION='READ')       
  read(41,*)((fe(i,j),i=1,m),j=1,total_nd) 
  close(41)
  open(41,file='xd', ACTION='READ')  
  read(41,*)((xd(i,j),i=1,m),j=1,total_nd) 
  close(41)
  n_1d(:)=3
  do k=1,podnum
   xi(:,1)=pod_coef(:)
   zd(:)=fe(k,:) 
  call lagrange_interp_nd_value_test2 ( m, n_1d, a, b, total_nd,xd,zd, ni, xi, zi )
   pod_coef(k)=zi(ni)
  enddo

  deallocate ( a )
  deallocate ( b )
  deallocate ( n_1d )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( zd )
  deallocate ( ze )
  deallocate ( zi )
  deallocate (pod_coef_all_obv)
  deallocate (pod_coef_all)
  deallocate (fi )
  return 
  end 


subroutine smolyak_main (timestep,podnum, pod_coef, total_timestep , sparse_max )

!*****************************************************************************80
!
!! TEST01: sequence of sparse interpolants to an M-dimensional function.
!
!  Discussion:
!
!    We have functions that can generate a Lagrange interpolant to data
!    in M dimensions, with specified order or level in each dimension.
!
!    We use the Lagrange function as the inner evaluator for a sparse
!    grid procedure. 
!
!    The procedure computes sparse interpolants of levels 0 to SPARSE_MAX
!    to a given function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error,app_error1
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  real ( kind = 8 ), allocatable :: fi(:) 
  logical more
  integer ( kind = 4 ) nd,nd_lq
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
 ! real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ) r0
  real ( kind = 8 ), allocatable :: xd(:,:),xd_lq(:,:),pod_coef_all_obv(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:),fd(:)
  real ( kind = 8 ), allocatable :: ze(:),temp_coef(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)
  real ( kind = 8 ), allocatable :: fe(:,:) !exact values of function produced from fluidity
  integer :: m,i,j,k,ll
  integer, intent(in) :: timestep,total_timestep,podnum
  real, intent(inout) :: pod_coef(podnum) 
  real ( kind = 8 ) w_lq(total_timestep)
 !real ( kind = 8 ) fd(nd_lq)
  integer :: n1d
  real :: mean, positive,negtive,num_pos,num_neg
  real minin(podnum),maxin(podnum)
  allocate  (xd_lq(podnum,total_timestep))
   
  write ( *, '(a)' ) ''
  !write ( *, '(a)' ) 'TEST01:'
 ! write ( *, '(a)' ) '  Sparse interpolation for a function f(x) of M-dimensional argument.'
 ! write ( *, '(a)' ) '  Use a sequence of sparse grids of levels 0 through SPARSE_MAX.'
 ! write ( *, '(a)' ) '  Invoke a general Lagrange interpolant function to do this.'
 ! write ( *, '(a)' ) ''
 ! write ( *, '(a)' ) '  Compare the exact function and the interpolants at a grid of points.'
 ! write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The "order" is the sum of the orders of all the product grids'
  write ( *, '(a)' ) '  used to make a particular sparse grid.'
!
!  User input.
!
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i4)' ) '  Maximum sparse grid level = ', sparse_max

  call get_option(&
            '/reduced_model/pod_basis_formation/pod_basis_count', nsvd) 
  m = podnum
  nd= total_timestep 
  !nd=160
  nd_lq=200
  n1d=3
  ! podnum=3*nsvd
  allocate(pod_coef_all_obv(total_timestep,podnum))
  allocate(fd(nd_lq))
  allocate(temp_coef(m))
  !allocate(pod_coef_all(total_timestep,podnum))
  print *, 'dd'
  open(unit=61,file='coef_pod_all_obv')
  read(61,*)((pod_coef_all_obv(j,k),k=1,podnum),j=1,total_timestep)
  close(61)

     Do i=1,m   ! cannot use minval because some of the vector is null.          
          minin(i)=pod_coef_all_obv(1,i)
          maxin(i)=pod_coef_all_obv(1,i)      
         do j=1,total_timestep
           if(minin(i)>pod_coef_all_obv(j,i)) then
              minin(i)=pod_coef_all_obv(j,i)          
           endif
           if(maxin(i)<pod_coef_all_obv(j,i)) then
              maxin(i)=pod_coef_all_obv(j,i)
           endif
         enddo
        ENDDO 
   
  positive=0
  negtive=0
  mean=0
  num_pos=0
  num_neg=0
  do j=1,m
    do k=1,nd
      if(pod_coef_all_obv(k,j)>0) then 
         positive=positive+pod_coef_all_obv(k,j)
         num_pos=num_pos+1
      else 
         negtive=negtive+pod_coef_all_obv(k,j)
         num_neg=num_neg+1
      endif
      mean=mean+pod_coef_all_obv(k,j)
    enddo
  enddo
    mean=mean/(num_pos+num_neg)
    positive=positive/num_pos
    negtive=negtive/num_neg
  do j=1,m
   do k=1,nd
      xd_lq(j,k)=pod_coef_all_obv(k,j) 
   enddo
  enddo
  
  m=podnum
  allocate ( ind(1:m) )  ! level of each dimension, so has the size of m
  
  !  Define the region.
  !
  allocate ( a(1:m) )
  allocate ( b(1:m) )

  a(1:m) = minin(1:m)!positive!0.0D+00
  b(1:m) = maxin(1:m)!negtive!2.0D+00

  !  Define the interpolation evaluation information.
  ! open(40,file='total', ACTION='read') 
  ! read(40,*) total_nd
  ! close(40)
  allocate (fe(1:m,1:total_nd))
  ! open(41,file='coef_pod_all_fvalue', ACTION='READ')       
  ! read(41,*)((fe(i,j),i=1,m),j=1,total_nd) 
  ! close(41)
 
   ni = 1
   allocate ( xi(1:m,1:ni) )
   allocate ( fi(1:ni) ) 
   xi(:,1)=pod_coef(:)
   allocate ( ze(1:ni) )
  !call f_sinr ( m, ni, xi, ze ) 
  ! Compute a sequence of sparse grid interpolants of increasing level.
  
   allocate ( zpi(1:ni) )
   allocate ( zi(1:ni) )
   
  do k=1, podnum
    sparse_min = 0
    lagsum=0
      do j=1,nd_lq
           fd(j)=pod_coef_all_obv(j+1,k)  ! target function value is the l^th next timestep's pod_coefficient 
      enddo 
  do l_max = sparse_min, sparse_max

    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )  
     !    Output, integer ( kind = 4 ) C(0:L_MAX), the coefficients for objects 
     !    at sublevels 0 through L_MAX.
     !    Output, integer ( kind = 4 ) W(0:L_MAX), the number of objects at 
     !    sublevels 0 through L_MAX.

    zi(1:ni) = 0.0D+00
    nd_total = 0

    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.

      do 
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
       
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )
         
        allocate ( zd(1:nd))
        !call f_sinr ( m, nd, xd, zd )
        !get the values of interpolation points 
        !zd(1:nd)=fe(k,(nd_total+1):(nd_total + nd))
        do i=1, nd 
        if(nd_total .eq. 0) then
           open(10,file='xd', ACTION='WRITE')
         !stop 111
         else
           open(10,file='xd', position='append',ACTION='WRITE')
         endif
        write(10,*)  xd(:,i)
        close(10)
        
        !  call fluids()
          open(10,file='pod_coef_smolyak')
        !  read(10,*) temp_coef(1:m) 
          close(10)
        !  zd(i)=temp_coef(k)
        enddo
        
       
       !  r0 = ( b(1) - a(1) ) / real ( 3, kind = 8 )
       ! call rbf_weight ( m, nd_lq, xd_lq, r0, phi1, fd, w_lq ) 
       ! call rbf_interp_nd ( m, nd_lq, xd_lq, r0, phi1, w_lq, nd, xd, zd ) 
         
        !  Use the grid to evaluate the interpolant. 
         call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi )
!
!  Weighted the interpolant values and add to the sparse grid interpolant.
!
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)

        deallocate ( xd )
        deallocate ( zd )

        if ( .not. more ) then
          exit
        end if
      end do
    end do
   pod_coef(k)=zi(1)

!  Compare sparse interpolant and exact function at interpolation points.
!
   !  app_error = r8vec_norm_affine ( ni, zi, ze ) / real ( ni, kind = 8 )
   !  app_error1 = r8vec_norm_affine ( ni, zi, zd ) / real ( ni, kind = 8 )
   !  write ( *, '(2x,i2,2x,i8,2x,e8.2)' ) l, nd_total, app_error1 

    deallocate ( c )
    deallocate ( w )

  end do  ! do l_max = sparse_min, sparse_max
  enddo  
      
  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi ) 
  deallocate (fe)
  deallocate(xd_lq)
  deallocate(pod_coef_all_obv)
  deallocate(fd)
  deallocate(temp_coef)
  return
end
 


subroutine main_call()

!*****************************************************************************80
!
!! MAIN is the main program for SPARSE_INTERP_ND_PRB.
!
!  Discussion:
!
!    SPARSE_INTERP_ND_PRB tests the SPARSE_INTERP_ND library.
!
 
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) sparse_max

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_INTERP_ND_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the SPARSE_INTERP_ND library.'
  write ( *, '(a)' ) '  The R8LIB library is also required.'

  m = 1
  sparse_max = 9
 !  call test01 ( m, sparse_max )

  m = 2
  sparse_max = 9
  ! call test01 ( m, sparse_max )
 
  m = 3
  sparse_max = 9
  !call test01 ( m, sparse_max )

  m = 6
  sparse_max = 7
 !  call test01 ( m, sparse_max )
 
   m = 2
  sparse_max = 7
 ! call testfluidity(m, sparse_max )

!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_INTERP_ND_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end



subroutine test01 ( m, sparse_max )

!*****************************************************************************80
!
!! TEST01: sequence of sparse interpolants to an M-dimensional function.
!
!  Discussion:
!
!    We have functions that can generate a Lagrange interpolant to data
!    in M dimensions, with specified order or level in each dimension.
!
!    We use the Lagrange function as the inner evaluator for a sparse
!    grid procedure. 
!
!    The procedure computes sparse interpolants of levels 0 to SPARSE_MAX
!    to a given function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
 ! real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Sparse interpolation for a function f(x) of M-dimensional argument.'
  write ( *, '(a)' ) '  Use a sequence of sparse grids of levels 0 through SPARSE_MAX.'
  write ( *, '(a)' ) '  Invoke a general Lagrange interpolant function to do this.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Compare the exact function and the interpolants at a grid of points.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The "order" is the sum of the orders of all the product grids'
  write ( *, '(a)' ) '  used to make a particular sparse grid.'
!
!  User input.
!
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i4)' ) '  Maximum sparse grid level = ', sparse_max

  allocate ( ind(1:m) )
!
!  Define the region.
!
  allocate ( a(1:m) )
  allocate ( b(1:m) )

  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00
!
!  Define the interpolation evaluation information.
!
  ni = 100
  seed = 123456789
  allocate ( xi(m,ni) )
  call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )

  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni

  allocate ( ze(1:ni) )
  call f_sinr ( m, ni, xi, ze )
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''

  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )

  sparse_min = 0

  do l_max = sparse_min, sparse_max

    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0

    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.

      do
!
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )
        allocate ( zd(1:nd) )
        call f_sinr ( m, nd, xd, zd )
!
!  Use the grid to evaluate the interpolant.
!
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi )
!
!  Weighted the interpolant values and add to the sparse grid interpolant.
!
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)

        deallocate ( xd )
        deallocate ( zd )

        if ( .not. more ) then
          exit
        end if

      end do

    end do
!
!  Compare sparse interpolant and exact function at interpolation points.
!
    app_error = r8vec_norm_affine ( ni, zi, ze ) / real ( ni, kind = 8 )

    write ( *, '(2x,i2,2x,i8,2x,e8.2)' ) l, nd_total, app_error

    deallocate ( c )
    deallocate ( w )

  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )

  return
end
subroutine test02 ( m, sparse_max )

!*****************************************************************************80
!
!! TEST01: sequence of sparse interpolants to an M-dimensional function.
!
!  Discussion:
!
!    We have functions that can generate a Lagrange interpolant to data
!    in M dimensions, with specified order or level in each dimension.
!
!    We use the Lagrange function as the inner evaluator for a sparse
!    grid procedure. 
!
!    The procedure computes sparse interpolants of levels 0 to SPARSE_MAX
!    to a given function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
  !real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)
  integer :: j,k
  integer :: total_timestep, num_pod
  real :: pod_coef_all_obv(2001, 72)
  total_timestep=2001
  num_pod=72
  
  open(unit=61,file='coef_pod_all_obv')
  read(61,*)((pod_coef_all_obv(j,k),k=1,num_pod),j=1,total_timestep)
  close(61)
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Sparse interpolation for a function f(x) of M-dimensional argument.'
  write ( *, '(a)' ) '  Use a sequence of sparse grids of levels 0 through SPARSE_MAX.'
  write ( *, '(a)' ) '  Invoke a general Lagrange interpolant function to do this.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Compare the exact function and the interpolants at a grid of points.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The "order" is the sum of the orders of all the product grids'
  write ( *, '(a)' ) '  used to make a particular sparse grid.'
!
!  User input.
!
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i4)' ) '  Maximum sparse grid level = ', sparse_max

  allocate ( ind(1:m) )
!
!  Define the region.
!
  allocate ( a(1:m) )
  allocate ( b(1:m) )

  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00
!
!  Define the interpolation evaluation information.
!
  ni = 100
  seed = 123456789
  allocate ( xi(m,ni) )
 

  do j=1,m
   do k=1,ni
   xi(j,k)=pod_coef_all_obv(k,j)
  enddo
  enddo
  call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )

  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni

  allocate ( ze(1:ni) )
  call f_sinr ( m, ni, xi, ze )
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''

  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )

  sparse_min = 0

  do l_max = sparse_min, sparse_max

    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0

    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.

      do
!
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )
        allocate ( zd(1:nd) )
        call f_sinr ( m, nd, xd, zd )
!
!  Use the grid to evaluate the interpolant.
!
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi )
!
!  Weighted the interpolant values and add to the sparse grid interpolant.
!
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)

        deallocate ( xd )
        deallocate ( zd )

        if ( .not. more ) then
          exit
        end if

      end do

    end do
!
!  Compare sparse interpolant and exact function at interpolation points.
!
    app_error = r8vec_norm_affine ( ni, zi, ze ) / real ( ni, kind = 8 )

    write ( *, '(2x,i2,2x,i8,2x,e8.2)' ) l, nd_total, app_error

    deallocate ( c )
    deallocate ( w )

  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )

  return
end
 

subroutine f_sinr ( m, n, x, z )

!*****************************************************************************80
!
!  F_SINR is a scalar function of an M-dimensional argument, to be interpolated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the points.
!
!    Output, real ( kind = 8 ) Z(N), the value of the function at each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) r(n)
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) z(n)
      
  r(1:n) = sqrt ( sum ( x(1:m,1:n)**2, dim = 1 ) )
  z(1:n) = sin ( r(1:n) )

  return
end
subroutine f_sinr_test ( m, n, x, z )

!*****************************************************************************80
!
!! F_SINR is a scalar function of an M-dimensional argument, to be interpolated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the points.
!
!    Output, real ( kind = 8 ) Z(N), the value of the function at each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) r(n)
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) z(n)
      
  !r(1:n) = sqrt ( sum ( x(1:m,1:n)**2, dim = 1 ) )
 ! z(1:n) = sin ( r(1:n) )
   
  integer :: j,k
  integer :: total_timestep, num_pod
  real :: pod_coef_all_obv(2001, 72)
  total_timestep=2001
  num_pod=72
  
  open(unit=61,file='coef_pod_all_obv')
  read(61,*)((pod_coef_all_obv(j,k),k=1,num_pod),j=1,total_timestep)
  close(61)
  z(1:n)=pod_coef_all_obv(1:n,1)
   
  return
end

subroutine cc_compute_points ( n, points )

!*****************************************************************************80
!
!! CC_COMPUTE_POINTS: abscissas of a Clenshaw Curtis rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [-1,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) POINTS(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) points(n)

  if ( n < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CC_COMPUTE_POINTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop

  else if ( n == 1 ) then

    points(1) = 0.0D+00

  else

    do i = 1, n
      points(i) = cos ( real ( n - i, kind = 8 ) * pi &
                      / real ( n - 1, kind = 8 ) )
    end do

    points(1) = -1.0D+00
    if ( mod ( n, 2 ) == 1 ) then
      points((n+1)/2) = 0.0D+00
    end if
    points(n) = +1.0D+00

  end if

  return
end
subroutine comp_next ( n, k, a, more, h, t )

!*****************************************************************************80
!
!! COMP_NEXT computes the compositions of the integer N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K parts is an ordered sequence
!    of K nonnegative integers which sum to N.  The compositions (1,2,1)
!    and (1,1,2) are considered to be distinct.
!
!    The routine computes one composition on each call until there are no more.
!    For instance, one composition of 6 into 3 parts is
!    3+2+1, another would be 6+0+0.
!
!    On the first call to this routine, set MORE = FALSE.  The routine
!    will compute the first element in the sequence of compositions, and
!    return it, as well as setting MORE = TRUE.  If more compositions
!    are desired, call again, and again.  Each time, the routine will
!    return with a new composition.
!
!    However, when the LAST composition in the sequence is computed
!    and returned, the routine will reset MORE to FALSE, signaling that
!    the end of the sequence has been reached.
!
!    This routine originally used a SAVE statement to maintain the
!    variables H and T.  I have decided that it is safer
!    to pass these variables as arguments, even though the user should
!    never alter them.  This allows this routine to safely shuffle
!    between several ongoing calculations.
!
!
!    There are 28 compositions of 6 into three parts.  This routine will
!    produce those compositions in the following order:
!
!     I         A
!     -     ---------
!     1     6   0   0
!     2     5   1   0
!     3     4   2   0
!     4     3   3   0
!     5     2   4   0
!     6     1   5   0
!     7     0   6   0
!     8     5   0   1
!     9     4   1   1
!    10     3   2   1
!    11     2   3   1
!    12     1   4   1
!    13     0   5   1
!    14     4   0   2
!    15     3   1   2
!    16     2   2   2
!    17     1   3   2
!    18     0   4   2
!    19     3   0   3
!    20     2   1   3
!    21     1   2   3
!    22     0   3   3
!    23     2   0   4
!    24     1   1   4
!    25     0   2   4
!    26     1   0   5
!    27     0   1   5
!    28     0   0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 July 2008
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer whose compositions are desired.
!
!    Input, integer ( kind = 4 ) K, the number of parts in the composition.
!
!    Input/output, integer ( kind = 4 ) A(K), the parts of the composition.
!
!    Input/output, logical MORE, set by the user to start the computation,
!    and by the routine to terminate it.
!
!    Input/output, integer ( kind = 4 ) H, T, two internal parameters needed
!    for the computation.  The user should allocate space for these in the
!    calling program, include them in the calling sequence, but never alter
!    them!
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) h
  logical more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t
!
!  The first computation.
!
  if ( .not. more ) then

    t = n
    h = 0
    a(1) = n
    a(2:k) = 0
!
!  The next computation.
!
  else
!
!  If the first entry A(1) is positive, then set H to zero,
!  so that when we increment H, it points to A(1); we will decrement A(1) by 1
!  and increment A(2).
!
    if ( 1 < t ) then
      h = 0
    end if
!
!  Otherwise, A(1) is 0.  Then by H + 1 is the entry we incremented last time.
!  Set H = H + 1, zero A(H), adding all but one of its value to A(1),
!  and incrementing A(H+1) by 1.
!
    h = h + 1
    t = a(h)
    a(h) = 0
    a(1) = t - 1
    a(h+1) = a(h+1) + 1

  end if
!
!  This is the last element of the sequence if all the
!  items are in the last slot.
!
  more = ( a(k) /= n )

  return
end


 
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
subroutine lagrange_basis_1d ( nd, xd, ni, xi, lb ) 

!*****************************************************************************80
!
!! LAGRANGE_BASIS_1D evaluates a 1D Lagrange basis.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(ND), the interpolation nodes.
!
!    Input, integer ( kind = 4 ) NI, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XI(NI), the evaluation points.
!
!    Output, real ( kind = 8 ) LB(NI,ND), the value, at the I-th point XI, 
!    of the Jth basis function.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) lb(ni,nd)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  
  do i = 1, ni
    do j = 1, nd
      lb(i,j) = product ( ( xi(i) - xd(1:j-1)  ) / ( xd(j) - xd(1:j-1)  ) ) &
              * product ( ( xi(i) - xd(j+1:nd) ) / ( xd(j) - xd(j+1:nd) ) )
    end do
  end do

  return
end
subroutine lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_GRID2 sets an M-dimensional Lagrange interpolant grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) IND(M), the index or level of the 1D rule 
!    to be used in each dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
!
!    Output, real ( kind = 8 ) XD(M,ND), the points at which data was sampled.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(m)
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: x_1d(:)
  real ( kind = 8 ) xd(m,nd)
!
!  Compute the data points.
!
  xd(1:m,1:nd) = 0.0D+00
  do i = 1, m
    call order_from_level_135 ( ind(i), n )
    allocate ( x_1d(1:n) )
    call cc_compute_points ( n, x_1d )
    x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
                          + ( 1.0D+00 + x_1d(1:n) ) * b(i) )
    call r8vec_direct_product ( i, n, x_1d, m, nd, xd )
    deallocate ( x_1d )
  end do
  return
end
subroutine lagrange_interp_nd_size2 ( m, ind, nd )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_SIZE2 sizes an M-dimensional Lagrange interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) IND(M), the index or level of the 1D rule 
!    to be used in each dimension.
!
!    Output, integer ( kind = 4 ) ND, the number of points in the product grid.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(m)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nd
!
!  Determine the number of data points.
!
  nd = 1
  do i = 1, m
    call order_from_level_135 ( ind(i), n )
    nd = nd * n
  end do

  return
end
subroutine lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zi )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_VALUE2 evaluates an ND Lagrange interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) IND(M), the index or level of the 1D rule 
!    to be used in each dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
!
!    Input, real ( kind = 8 ) ZD(ND), the function evaluated at the points XD.
!
!    Input, integer ( kind = 4 ) NI, the number of points at which the 
!    interpolant is to be evaluated.
!
!    Input, real ( kind = 8 ) XI(M,NI), the points at which the interpolant 
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) ZI(NI), the interpolant evaluated at the 
!    points XI.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer   ni

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: value(:)
  real ( kind = 8 ) w(nd)
  real ( kind = 8 ), allocatable :: x_1d(:)
  real ( kind = 8 ) xi(m,ni)
  real ( kind = 8 ) zd(nd)
  real ( kind = 8 ) zi(ni)

  do j = 1, ni

    w(1:nd) = 1.0D+00

    do i = 1, m
      call order_from_level_135 ( ind(i), n )
      allocate ( x_1d(1:n) )
      allocate ( value(1:n) )
      call cc_compute_points ( n, x_1d )
      x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
                            + ( 1.0D+00 + x_1d(1:n) ) * b(i) )
      call lagrange_basis_1d ( n, x_1d, 1, xi(i,j), value )
      call r8vec_direct_product2 ( i, n, value, m, nd, w )
      deallocate ( value )
      deallocate ( x_1d )
    end do

    zi(j) = dot_product ( w, zd )

  end do

  return
end

subroutine lagrange_interp_nd_value ( m, n_1d, a, b, nd, zd, ni, xi, zi )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_VALUE evaluates an ND Lagrange interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N_1D(M), the order of the 1D rule to be used
!    in each dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
!
!    Input, real ( kind = 8 ) ZD(ND), the function evaluated at the points XD.
!
!    Input, integer ( kind = 4 ) NI, the number of points at which the 
!    interpolant is to be evaluated.
!
!    Input, real ( kind = 8 ) XI(M,NI), the points at which the interpolant is 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) ZI(NI), the interpolant evaluated at the 
!    points XI.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_1d(m)
  real ( kind = 8 ), allocatable :: value(:)
  real ( kind = 8 ) w(nd)
  real ( kind = 8 ), allocatable :: x_1d(:)
  real ( kind = 8 ) xi(m,ni)
  real ( kind = 8 ) zd(nd)
  real ( kind = 8 ) zi(ni)

  do j = 1, ni

    w(1:nd) = 1.0D+00

    do i = 1, m
      n = n_1d(i)
      allocate ( x_1d(1:n) )
      allocate ( value(1:n) )
      call cc_compute_points ( n, x_1d )
      x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
                            + ( 1.0D+00 + x_1d(1:n) ) * b(i) )
      call lagrange_basis_1d ( n, x_1d, 1, xi(i,j), value )
      call r8vec_direct_product2 ( i, n, value, m, nd, w )
      deallocate ( value )
      deallocate ( x_1d )
    end do

    zi(j) = dot_product ( w, zd )

  end do

  return
end
subroutine lagrange_interp_nd_value_test2 ( m, n_1d, a, b, nd, xd,zd, ni, xi, zi )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_VALUE evaluates an ND Lagrange interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N_1D(M), the order of the 1D rule to be used
!    in each dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
!
!    Input, real ( kind = 8 ) ZD(ND), the function evaluated at the points XD.
!
!    Input, integer ( kind = 4 ) NI, the number of points at which the 
!    interpolant is to be evaluated.
!
!    Input, real ( kind = 8 ) XI(M,NI), the points at which the interpolant is 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) ZI(NI), the interpolant evaluated at the 
!    points XI.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j 
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_1d(m)
  real ( kind = 8 ), allocatable :: value(:)
  real ( kind = 8 ) w(nd)
  real ( kind = 8 ), allocatable :: x_1d(:)
  real ( kind = 8 ) xi(m,ni),xd(m,nd)
  real ( kind = 8 ) zd(nd)
  real ( kind = 8 ) zi(ni)
  integer ::k
    integer:: total_timestep, num_pod
   
  do j = 1, ni

    w(1:nd) = 1.0D+00

    do i = 1, m
      n = n_1d(i)
      allocate ( x_1d(1:n) )
      allocate ( value(1:n) )
      call cc_compute_points ( n, x_1d )
      x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
                            + ( 1.0D+00 + x_1d(1:n) ) * b(i) )
    !  call lagrange_basis_1d ( n, x_1d, 1, xi(i,j), value )
      call lagrange_basis_1d ( n, xd(i,:), 1, xi(i,j), value ) 
      ! print *, 'value', value
      
      call r8vec_direct_product2 ( i, n, value, m, nd, w )
      !print *, 'value', value
      deallocate ( value )
      deallocate ( x_1d )
    end do
    !print *, 'weight', w
    !print *, 'zd', zd
    zi(j) = dot_product ( w, zd )
    print *, 'zj', zi(j)
  end do

  return
end
subroutine order_from_level_135 ( l, n )

!*****************************************************************************80
!
!! ORDER_FROM_LEVEL_135 evaluates the 135 level-to-order relationship.
!
!  Discussion:
!
!    Clenshaw Curtis rules, and some others, often use the following
!    scheme:
!
!    L: 0  1  2  3   4   5
!    N: 1  3  5  9  17  33 ... 2^L+1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, the level, which should be 0 or greater.
!
!    Output, integer ( kind = 4 ) N, the order.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  if ( l < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ORDER_FROM_LEVEL_135 - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of L!'
    stop
  else if ( l == 0 ) then
    n = 1
  else
    n = ( 2 ** l ) + 1
  end if

  return
end
subroutine smolyak_coefficients ( l_max, m, c, w )

!*****************************************************************************80
!
!! SMOLYAK_COEFFICIENTS returns the Smolyak coefficients and counts.
!
!  Discussion:
!
!    The Smolyak sparse interpolant can be written as:
!
!      A(L,M)(X) = sum ( L-M+1 <= |L| <= L_max ) 
!        C(|L|) * g(l1)(x1) * g(l2)(x2) * ... * g(lm)(xm).
!
!    where:
!
!    * L=(l1,l2,...,lm) is a vector of M nonnegative integers;
!    * |L| is the sum of the entries of L;
!    * X=(x1,x2,...,xm) is an M-dimensional point in a product space;
!    * g(i)(xj) is the i-th 1-d interpolation function in dimension j;
!
!    Note that:
!
!    * W(|L|) will represent the number of distinct interpolants for which
!      the sublevel, or sum of the L vector entries, is |L|;
!
!    * the coefficients C and counts W will be zero for sublevels 
!      0 through L_MAX - M (and MATLAB indices 1 through L_MAX-M+1).
!
!    * it will be the case that W' * C = 1, essentially because the interpolant
!      to the identity function must be the identity function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L_MAX, the (maximum) level.
!    0 <= L_MAX.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!    1 <= M.
!
!    Output, integer ( kind = 4 ) C(0:L_MAX), the coefficients for objects 
!    at sublevels 0 through L_MAX.
!
!    Output, integer ( kind = 4 ) W(0:L_MAX), the number of objects at 
!    sublevels 0 through L_MAX.
!
  implicit none

  integer ( kind = 4 ) l_max

  integer ( kind = 4 ) c(0:l_max)
  !integer ( kind = 4 ) i4_choose
  !integer ( kind = 4 ) i4_mop
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_min
  integer ( kind = 4 ) m
  integer ( kind = 4 ) w(0:l_max)

  l_min = max ( l_max - m + 1, 0 )

  c(0:l_min-1) = 0
  do l = l_min, l_max
    c(l) = i4_mop ( l_max - l ) * i4_choose ( m - 1, l_max - l )
  end do

  w(0:l_min-1) = 0
  do l = l_min, l_max
    w(l) = i4_choose ( l + m - 1, m - 1 )
  end do



  contains
   function i4_mop ( i ) result (i4mop)

!*****************************************************************************80
!
!! I4_MOP returns the I-th power of -1 as an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the power of -1.
!
!    Output, integer ( kind = 4 ) I4_MOP, the I-th power of -1.
!
  !implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4mop

  if ( mod ( i, 2 ) == 0 ) then
    i4mop = 1
  else
    i4mop = -1
  end if
 
  end function i4_mop

  function i4_choose ( n, k ) result(i4choose)

!*****************************************************************************80
!
!! I4_CHOOSE computes the binomial coefficient C(N,K) as an I4.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, are the values of N and K.
!
!    Output, integer ( kind = 4 ) I4_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4choose
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0

  else if ( mn == 0 ) then

    value = 1

  else

    mx = max ( k, n - k )
    value = mx + 1

    do i = 2, mn
      value = ( value * ( mx + i ) ) / i
    end do

  end if

  i4choose = value

  return
end
  
end 

subroutine r8mat_uniform_abvec ( m, n, a, b, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_ABVEC returns a scaled pseudorandom R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    A(I) <= R(I,J) <= B(I)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_UNIFORM_ABVEC - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if
     ! print *, seed
      r(i,j) = a(i) + ( b(i) - a(i) ) * real ( seed, kind = 8 ) &
        * 4.656612875D-10

    end do
  end do

  return
end

function r8vec_norm_affine1 ( n, v0, v1 )

!*****************************************************************************80
!
!! R8VEC_NORM_AFFINE returns the affine norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The affine vector L2 norm is defined as:
!
!      R8VEC_NORM_AFFINE(V0,V1)
!        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the vectors.
!
!    Input, real ( kind = 8 ) V0(N), the base vector.
!
!    Input, real ( kind = 8 ) V1(N), the vector whose affine norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM_AFFINE, the L2 norm of V1-V0.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) r8vec_norm_affine1
  real ( kind = 8 ) v0(n)
  real ( kind = 8 ) v1(n)

  r8vec_norm_affine1 = sqrt ( sum ( ( v0(1:n) - v1(1:n) )**2 ) )

  return
end

subroutine r8vec_direct_product22 ( factor_index, factor_order, factor_value, &
  factor_num, point_num, w )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out the task involving the weights W.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      W(1:4) = ( 2, 3, 5, 7 )
!
!    Rule 2:
!      Order = 3
!      W(1:3) = ( 11, 13, 17 )
!
!    Rule 3:
!      Order = 2
!      W(1:2) = ( 19, 23 )
!
!    Product Rule:
!      Order = 24
!      W(1:24) =
!        ( 2 * 11 * 19 )
!        ( 3 * 11 * 19 )
!        ( 4 * 11 * 19 )
!        ( 7 * 11 * 19 )
!        ( 2 * 13 * 19 )
!        ( 3 * 13 * 19 )
!        ( 5 * 13 * 19 )
!        ( 7 * 13 * 19 )
!        ( 2 * 17 * 19 )
!        ( 3 * 17 * 19 )
!        ( 5 * 17 * 19 )
!        ( 7 * 17 * 19 )
!        ( 2 * 11 * 23 )
!        ( 3 * 11 * 23 )
!        ( 5 * 11 * 23 )
!        ( 7 * 11 * 23 )
!        ( 2 * 13 * 23 )
!        ( 3 * 13 * 23 )
!        ( 5 * 13 * 23 )
!        ( 7 * 13 * 23 )
!        ( 2 * 17 * 23 )
!        ( 3 * 17 * 23 )
!        ( 5 * 17 * 23 )
!        ( 7 * 17 * 23 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, real ( kind = 8 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, real ( kind = 8 ) W(POINT_NUM), the elements of the
!    direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) START, the first location of a block of values
!    to set.
!
!    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values
!    to set.
!
!    Local, integer ( kind = 4 ) SKIP, the distance from the current value 
!    of START to the next location of a block of values to set.
!
!    Local, integer ( kind = 4 ) REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  real ( kind = 8 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  real ( kind = 8 ) w(point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    w(1:point_num) = 1.0D+00
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      w(start:start+contig-1) = w(start:start+contig-1) * factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end

subroutine r8vec_direct_product1 ( factor_index, factor_order, factor_value, &
  factor_num, point_num, x )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out that task for the abscissas X.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      X(1:4) = ( 1, 2, 3, 4 )
!
!    Rule 2:
!      Order = 3
!      X(1:3) = ( 10, 20, 30 )
!
!    Rule 3:
!      Order = 2
!      X(1:2) = ( 100, 200 )
!
!    Product Rule:
!      Order = 24
!      X(1:24) =
!        ( 1, 10, 100 )
!        ( 2, 10, 100 )
!        ( 3, 10, 100 )
!        ( 4, 10, 100 )
!        ( 1, 20, 100 )
!        ( 2, 20, 100 )
!        ( 3, 20, 100 )
!        ( 4, 20, 100 )
!        ( 1, 30, 100 )
!        ( 2, 30, 100 )
!        ( 3, 30, 100 )
!        ( 4, 30, 100 )
!        ( 1, 10, 200 )
!        ( 2, 10, 200 )
!        ( 3, 10, 200 )
!        ( 4, 10, 200 )
!        ( 1, 20, 200 )
!        ( 2, 20, 200 )
!        ( 3, 20, 200 )
!        ( 4, 20, 200 )
!        ( 1, 30, 200 )
!        ( 2, 30, 200 )
!        ( 3, 30, 200 )
!        ( 4, 30, 200 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, real ( kind = 8 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, real ( kind = 8 ) X(FACTOR_NUM,POINT_NUM), the elements of
!    the direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) START, the first location of a block of 
!    values to set.
!
!    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values 
!    to set.
!
!    Local, integer ( kind = 4 ) SKIP, the distance from the current value 
!    of START to the next location of a block of values to set.
!
!    Local, integer ( kind = 4 ) REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  real ( kind = 8 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  real ( kind = 8 ) x(factor_num,point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    x(1:factor_num,1:point_num) = 0.0D+00
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      x(factor_index,start:start+contig-1) = factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end

subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end 

end module smolyak
