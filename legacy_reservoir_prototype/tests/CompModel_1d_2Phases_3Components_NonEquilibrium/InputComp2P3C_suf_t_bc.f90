program main
  implicit none
  integer :: ndim
  real, dimension( : ), allocatable :: field
  integer :: idim, iphase
  integer, parameter :: nphase = 2

  open( 5, file = 'filedim', status = 'unknown' )
  read( 5, * ) ndim
  close( 5 )

  allocate( field( ndim ))

  field = 0.
  field( 1 ) = 0.5

  do iphase = 2, nphase
     field( 1 + ( iphase - 1 ) * 2 ) = 0.5
     field( 2 + ( iphase - 1 ) * 2 ) = 0.
  end do

  do idim = 1, ndim
     print*, field( idim )
  end do


end program main


