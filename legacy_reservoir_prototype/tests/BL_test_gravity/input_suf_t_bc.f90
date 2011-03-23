
program main
  implicit none
  integer :: ndim
  real, dimension( : ), allocatable :: field
  integer :: idim

  open( 5, file = 'filedim', status = 'unknown' )
  read( 5, * ) ndim
  close( 5 )

  allocate( field( ndim1 ))

  field = 0.
  field( 1 ) = 1.

  do idim = 1, ndim
     print*, field( idim )
  end do


end program main




