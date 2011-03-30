program main
  implicit none
  integer :: ndim
  integer, parameter :: nphase = 2
  real, dimension( : ), allocatable :: field
  integer :: u_nonods, u_nods


  open( 5, file = 'filedim', status = 'unknown' )
  read( 5, * ) ndim
  close( 5 )

  allocate( field( ndim ))
  u_nonods = ndim / 2

  field = 0.

  do u_nods = 1, ndim

     if( u_nods <= u_nonods )then
        field( u_nods ) = 981. * ( 1.05 - .71 ) * 0.02 ! -981. * 1.05 ! -981. * ( 1.05 - .71 )
     elseif( u_nods > u_nonods ) then
        field( u_nods ) = 0.
     end if

     print*, field( u_nods )

  end do

end program main
