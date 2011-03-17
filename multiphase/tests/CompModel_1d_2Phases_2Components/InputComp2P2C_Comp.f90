
program main
  implicit none
  integer :: ndim
  integer, parameter :: nphase = 2, ncomp = 2
  real, dimension( : ), allocatable :: field !( cv_pha_nonods * ncomp )
  integer :: i, ndim2

  open( 5, file = 'filedim', status = 'unknown' )
  read( 5, * ) ndim
  close( 5 )

  allocate( field( ndim ))

  field = 0.

  ndim2 = ndim / ncomp

  do i = ndim2 + 1, ndim, 1
     field( i ) = 1.
  end do

  do i = 1, ndim
     print*, field( i )
  end do

end program main
