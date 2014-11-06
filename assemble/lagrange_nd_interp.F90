#include "fdebug.h"

module lag_nd_interp
use spud
 contains 

 subroutine lag_nd_interp_main (podnum, l, total_timestep ) 
  implicit none 
  integer, intent(in) :: total_timestep,podnum, l
  !real, intent(inout) :: pod_coef(podnum)  
  integer :: m,nd  
  real ( kind = 8 ) int_error  
  ! real ( kind = 8 ) x1d(n1d)  
  integer :: nsvd
  integer ::  i,k,j 
  real , dimension(:,:), allocatable ::  pod_coef_all_obv, pod_coef_all, pod_coef_all_lag_tmp   
  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:) 
  integer ( kind = 4 ), allocatable :: n_1d(:) 
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)  
  integer :: inteval
  call get_option(&
            '/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
  m = podnum
  inteval=2
  nd= total_timestep/inteval  
  ni = total_timestep!nd
  allocate(pod_coef_all_obv(total_timestep,podnum))
  allocate(pod_coef_all(total_timestep,podnum)) 
  allocate(xd(m,nd))
  allocate ( zd(1:nd)) 
  allocate ( n_1d(1:m) ) 
  allocate ( a(1:m) )
  allocate ( b(1:m) )
  allocate ( xi(1:m,1:ni) )
  allocate ( ze(1:ni) )
  allocate ( zi(1:ni))
  n_1d(1:m) =  total_timestep/inteval 
  open(unit=61,file='coef_pod_all_obv')
  read(61,*)((pod_coef_all_obv(j,k),k=1,podnum),j=1,total_timestep)
  close(61) 
  open(unit=81,file='coef_pod_all')
  read(81,*)((pod_coef_all(j,k),k=1,podnum),j=1,total_timestep)
  close(81)  
  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00  
   do j=1,m 
   do k=1,nd
     xd(j,k)=pod_coef_all(k*inteval,j)!k
    ! print *, 'xd(j,k)', xd(j,k) 
   enddo 
   enddo 
  
 ! call f_sinr_test_ann ( m, nd, xd, zd )
   do k=1,nd 
     zd(k)=pod_coef_all(k*inteval,l)
   enddo
  !zd(1:nd)=pod_coef_all(1:nd,l)
  write ( *, '(a)' ) ''
 ! write ( *, '(a)' ) '         Zinterp          Zexact      Error'
  write ( *, '(a)' ) '' 
  seed = 123456789 
  call r8mat_uniform_01 ( m, ni, seed, xi )
   do j=1,m
   do k=1,ni
     xi(j,k)=pod_coef_all(k,j)
   enddo
   enddo
  
  !   call f_sinr_test_ann ( m, nd, xd, ze)
  ze(1:ni)=pod_coef_all(1:ni,l) 
  call lagrange_interp_nd_value_test2 ( m, n_1d, a, b, nd, xd,zd, ni, xi, zi )

  ! do j = 1,  ni
    ! write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) &
  !    zi(j), ze(j), abs ( zi(j) - ze(j)) 
  ! end do
     if(l==1) then
       open(40,file='coef_pod_all_lag_tmp')
        write(40,*)(zi(:))
       close(40) 
     else
        open(101,file='coef_pod_all_lag_tmp', position='append',ACTION='WRITE')
        write(101,*)(zi(i),i=1,ni)
        close(101) 
     !  stop 111
    endif

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
 
  return
end
 


subroutine   testscatterdata_ann ( )
 implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  external f_sinr
 ! integer ( kind = 4 ) j,k
  integer ( kind = 4 ) m
  integer ( kind = 4 ), allocatable :: n_1d(:)
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real :: t,t1
  integer :: j,k
  integer :: total_timestep, num_pod
  real :: pod_coef_all_obv(2001, 72)
  total_timestep=2001
  num_pod=72
  
  open(unit=61,file='coef_pod_all_obv')
  read(61,*)((pod_coef_all_obv(j,k),k=1,num_pod),j=1,total_timestep)
  close(61) 
    
   !z(1:n)=pod_coef_all_obv(1:n,1) 
  m = num_pod

  allocate ( n_1d(1:m) )
  n_1d(1:m) = 1000
  nd=n_1d(1)
  allocate ( a(1:m) )
  allocate ( b(1:m) )
  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

   !call lagrange_interp_nd_size ( m, n_1d, nd )

  allocate ( xd(1:m,1:nd) )
  allocate ( zd(1:nd) )
   
   do j=1,m 
   do k=1,nd
     xd(j,k)=pod_coef_all_obv(k,j)!k
    ! print *, 'xd(j,k)', xd(j,k) 
   enddo 
   enddo 
  !call lagrange_interp_nd_grid ( m, n_1d, a, b, nd, xd )
  call f_sinr_test_ann ( m, nd, xd, zd )
  
!  Evaluate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '         Zinterp          Zexact      Error'
  write ( *, '(a)' ) ''

  ni = nd
  seed = 123456789
  allocate ( xi(1:m,1:ni) )
  call r8mat_uniform_01 ( m, ni, seed, xi )
  xi(:,:)=xd(:,:)
  allocate ( ze(1:ni) )
    call f_sinr_test_ann ( m, nd, xd, ze)
  ! call f_sinr ( m, nd, xd, zd )

  allocate ( zi(1:ni))
  call lagrange_interp_nd_value_test2 ( m, n_1d, a, b, nd, xd,zd, ni, xi, zi )

  do j = 1,  ni
    write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      zi(j), ze(j), abs ( zi(j) - ze(j)) 
  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( n_1d )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( zd )
  deallocate ( ze )
  deallocate ( zi )
 
  return
end

 
subroutine f_obv ( m, n, x, z )

 
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


subroutine f_sinr ( m, n, x, z )

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
    z(1:n) = sin ( r(1:n) ) !original
   
  return
end

subroutine f_sinr_test ( m, n, x, z ) 
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n !number of points

  real ( kind = 8 ) r(n)
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) z(n)
  real :: t
  
  integer :: j,k
  integer :: total_timestep, num_pod
  real :: pod_coef_all_obv(2001, 72)
  total_timestep=2001
  num_pod=72
  
  open(unit=61,file='coef_pod_all_obv')
  read(61,*)((pod_coef_all_obv(j,k),k=1,num_pod),j=1,total_timestep)
  close(61) 
    r(1:n) = sqrt ( sum ( x(1:m,1:n)**2, dim = 1 ) )
    z(1:n) = sin ( r(1:n) ) !original
   ! z(1:n)=pod_coef_all_obv(1:n,1)
    
    t=0
    do j=1, n
      z(j)=t  !!test data
     t=t+0.5
    enddo
  return
end

subroutine f_sinr_test_ann ( m, n, x, z ) 
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n !number of points

  real ( kind = 8 ) r(n)
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) z(n)
  real :: t
  
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
subroutine f_poly352 ( m, n, x, z )

!*****************************************************************************80
!
!! F_POLY253 is a scalar function of a 3-dimensional argument, to be interpolated.
!
!  Discussion:
!
!    The polynomial is of maximum degrees 3, 5, and 2, in the three variables.
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
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the points.
!
!    Output, real ( kind = 8 ) Z(N,1), the value of the function at each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) z(n)

  z(1:n) = 1.0 + x(1,1:n) ** 2 * x(2,1:n) ** 5 * x(3,1:n) ** 2 &
    + x(1,1:n) * x(2,1:n) ** 2 * x(3,1:n) + x(1,1:n) ** 3

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
       !print *, 'xddddxiiii', xi(i),xd(i)
      lb(i,j) = product ( ( xi(i) - xd(1:j-1)  ) / ( xd(j) - xd(1:j-1)  ) ) &
              * product ( ( xi(i) - xd(j+1:nd) ) / ( xd(j) - xd(j+1:nd) ) )
    end do
  end do

  return
end
subroutine lagrange_interp_nd_grid ( m, n_1d, a, b, nd, xd )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_GRID sets an M-dimensional Lagrange interpolant grid.
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
!    Input, integer ( kind = 4 ) N_1D(M), the order of the 1D rule to be used
!    in each dimension.
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
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_1d(m)
  real ( kind = 8 ), allocatable :: x_1d(:)
  real ( kind = 8 ) xd(m,nd)
!
!  Compute the data points.
!
  xd(1:m,1:nd) = 0.0D+00
  do i = 1, m
    n = n_1d(i)
    allocate ( x_1d(1:n) )
    call cc_compute_points ( n, x_1d )
    x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
                          + ( 1.0D+00 + x_1d(1:n) ) * b(i) )
    call r8vec_direct_product ( i, n, x_1d, m, nd, xd )
    deallocate ( x_1d )
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
    call order_from_level_135 ( ind(i), n ) ! n = ( 2 ** l ) + 1
    allocate ( x_1d(1:n) )
    call cc_compute_points ( n, x_1d )
    x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
                          + ( 1.0D+00 + x_1d(1:n) ) * b(i) )
    call r8vec_direct_product ( i, n, x_1d, m, nd, xd )
    deallocate ( x_1d )
  end do

  return
end
subroutine lagrange_interp_nd_size ( m, n_1d, nd )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_SIZE sizes an M-dimensional Lagrange interpolant.
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
!    Output, integer ( kind = 4 ) ND, the number of points in the product grid.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) n_1d(m)
  integer ( kind = 4 ) nd
!
!  Determine the number of data points.
!
  nd = product ( n_1d(1:m) )

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
    call order_from_level_135 ( ind(i), n ) ! n = ( 2 ** l ) + 1
    nd = nd * n
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
   real :: pod_coef_all_obv(2001, 72)
  
  total_timestep=2001
  num_pod=72
  
 ! open(unit=61,file='coef_pod_all_obv')
!  read(61,*)((pod_coef_all_obv(j,k),k=1,num_pod),j=1,total_timestep)
 ! close(61)
  
   do j=1,m
   do k=1,nd
    !  xd(j,k)=k 
   enddo
   enddo
 !   xi(:,:)=xd(:,:)
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

subroutine lagrange_interp_nd_value_test ( m, n_1d, a, b, nd, zd, ni, xi, zi )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_VALUE evaluates an ND Lagrange interpolant.
!
!  Licensing:
!
 
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
  integer ( kind = 4 ) i,k
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_1d(m)
  real ( kind = 8 ), allocatable :: value(:)
  real ( kind = 8 ) w(nd)
  real ( kind = 8 ), allocatable :: x_1d(:)
  real ( kind = 8 ) xi(m,ni),xd(m,nd)
  real ( kind = 8 ) zd(nd)
  real ( kind = 8 ) zi(ni)
  integer:: total_timestep, num_pod
   real :: pod_coef_all_obv(2001, 72)
  
  total_timestep=2001
  num_pod=72
  
  open(unit=61,file='coef_pod_all_obv')
  read(61,*)((pod_coef_all_obv(j,k),k=1,num_pod),j=1,total_timestep)
  close(61)
  
   do j=1,m
   do k=1,nd
     xd(j,k)=k !pod_coef_all_obv(k,j)
   enddo
   enddo
   xi=xd
  ! print *, '1111xd(i,:)xi(i,:)', xd(:,:),xi(:,:) 
   !stop 111
  do j = 1, ni

    w(1:nd) = 1.0D+00

    do i = 1, m
      n = n_1d(i)
      allocate ( x_1d(1:n) )
     ! allocate ( value(1:n) )
       allocate ( value(1:nd) )
      call cc_compute_points ( n, x_1d )
      x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
                            + ( 1.0D+00 + x_1d(1:n) ) * b(i) ) 
     ! print *, 'before lagrange_basis_1d'
      print *, 'xd(i,:)xi(i,:)', xd(i,:),xi(i,:)
     ! call lagrange_basis_1d ( nd, xd(i,:), 1, xi(i,j), value)  
     call lagrange_basis_1d ( n, x_1d, 1, xi(i,j), value )
      print *, 'value', value
     ! print *, 'after lagrange_basis_1d'
      call r8vec_direct_product2 ( i, n, value, m, nd, w )
     ! print *, 'after r8vec_direct_product2'
       deallocate ( value )
     !  print *, 'after end'
      deallocate ( x_1d )
    end do
    print *, 'w ze before', w,zd
    zi(j) = dot_product ( w, zd )
    print *, 'w ze after', w,zd
  end do
 
  return
end

subroutine lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zi )

!*****************************************************************************80
!
!  LAGRANGE_INTERP_ND_VALUE2 evaluates an ND Lagrange interpolant.
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
  integer ( kind = 4 ) ni

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
    print *, 'e ze before',  w, zd
    zi(j) = dot_product ( w, zd )
    print *, 'e ze after',  w, zd
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
!    28 September 2012
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
    write ( *, '(a)' ) ''
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
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end

function r8vec_norm_affine ( n, v0, v1 )

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

  real ( kind = 8 ) r8vec_norm_affine
  real ( kind = 8 ) v0(n)
  real ( kind = 8 ) v1(n)

  r8vec_norm_affine = sqrt ( sum ( ( v0(1:n) - v1(1:n) )**2 ) )

  return
end

subroutine r8vec_direct_product ( factor_index, factor_order, factor_value, &
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
subroutine r8vec_direct_product2 ( factor_index, factor_order, factor_value, &
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
!        ( 5 * 11 * 19 )
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

 
 
subroutine testscatterdata_1d ( )
 implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  external f_sinr
  integer ( kind = 4 ) j,k
  integer ( kind = 4 ) m
  integer ( kind = 4 ), allocatable :: n_1d(:)
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real :: t,t1
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Interpolate in 2D, using orders.'
  write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_GRID sets the interpolant.'
  write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_VALUE evaluates it.'

  m = 1

  allocate ( n_1d(1:m) )
  n_1d(1:m) = 9
  nd=9
  allocate ( a(1:m) )
  allocate ( b(1:m) )
  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

   !call lagrange_interp_nd_size ( m, n_1d, nd )

  allocate ( xd(1:m,1:nd) )
  allocate ( zd(1:nd) )
   t=0 
   do j=1,m
      t1=1
   do k=1,nd
     xd(j,k)=k
     print *, 'xd(j,k)', xd(j,k)
     t1=t1+1
   enddo
     t=t+1
   enddo
  
  !call lagrange_interp_nd_grid ( m, n_1d, a, b, nd, xd )
  call f_sinr_test ( m, nd, xd, zd )
 !  call f_sinr ( m, nd, xd, zd )
!
!  Evaluate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '         Zinterp          Zexact      Error'
  write ( *, '(a)' ) ''

  ni = nd
  seed = 123456789
  allocate ( xi(1:m,1:ni) )
  call r8mat_uniform_01 ( m, ni, seed, xi )
  xi(:,:)=xd(:,:)
  allocate ( ze(1:ni) )
    call f_sinr_test ( m, nd, xd, ze)
  ! call f_sinr ( m, nd, xd, zd )

  allocate ( zi(1:ni))
  call lagrange_interp_nd_value_test2 ( m, n_1d, a, b, nd, xd,zd, ni, xi, zi )

  do j = 1,  ni
    write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      zi(j), ze(j), abs ( zi(j) - ze(j)) 
  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( n_1d )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( zd )
  deallocate ( ze )
  deallocate ( zi )
 
  return
end


subroutine testscatterdata_nd ( )
 implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: b(:)
  external f_sinr
  integer ( kind = 4 ) j,k
  integer ( kind = 4 ) m
  integer ( kind = 4 ), allocatable :: n_1d(:)
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real :: t,t1
   
  m = 133

  allocate ( n_1d(1:m) )
  n_1d(1:m) = 163
  nd=n_1d(1)
  allocate ( a(1:m) )
  allocate ( b(1:m) )
  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00

   !call lagrange_interp_nd_size ( m, n_1d, nd )

  allocate ( xd(1:m,1:nd) )
  allocate ( zd(1:nd) )
   
   do j=1,m 
   do k=1,nd
     xd(j,k)=k
     print *, 'xd(j,k)', xd(j,k) 
   enddo 
   enddo
   
  call f_sinr_test ( m, nd, xd, zd )
 !  call f_sinr ( m, nd, xd, zd )
!
!  Evaluate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '         Zinterp          Zexact      Error'
  write ( *, '(a)' ) ''

  ni = nd
  seed = 123456789
  allocate ( xi(1:m,1:ni) )
  call r8mat_uniform_01 ( m, ni, seed, xi )
  xi(:,:)=xd(:,:)
  allocate ( ze(1:ni) )
    call f_sinr_test ( m, nd, xd, ze)
  ! call f_sinr ( m, nd, xd, zd )

  allocate ( zi(1:ni))
  call lagrange_interp_nd_value_test2 ( m, n_1d, a, b, nd, xd,zd, ni, xi, zi )

  do j = 1,  ni
    write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      zi(j), ze(j), abs ( zi(j) - ze(j)) 
  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( n_1d )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( zd )
  deallocate ( ze )
  deallocate ( zi )
 
  return
end

end module lag_nd_interp
