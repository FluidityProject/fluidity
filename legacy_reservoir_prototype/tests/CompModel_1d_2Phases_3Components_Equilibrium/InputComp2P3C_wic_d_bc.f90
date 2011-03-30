program main
  implicit none
  integer :: ndim
  integer, dimension( : ), allocatable :: field
  integer :: i

  open( 5, file = 'filedim', status = 'unknown' )
  read( 5, * ) ndim
  close( 5 )

  allocate( field( ndim ))

  field = 0
  field( 1 ) = 1

  do i = 1, ndim
     print*, field( i )
  end do

end program main
