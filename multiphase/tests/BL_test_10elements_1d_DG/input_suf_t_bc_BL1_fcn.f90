program main
  implicit none
  integer :: ndim
  integer, parameter :: nphase = 2
  real, dimension( : ), allocatable :: field
  integer :: i, ndim2

  open( 5, file = 'filedim', status = 'unknown' )
  read( 5, * ) ndim
  close( 5 )

  allocate( field( ndim ))

  !field = 0.
  field( 1 ) = .5
  field( 2 ) = 0.

  do i = 2, nphase
     field( 1 + ( i - 1 ) * 2 ) = .5
     field( 2 + ( i - 1 ) * 2 ) = 0.
  end do

  do i = 1, ndim
     print*, field( i )
  end do

end program main
