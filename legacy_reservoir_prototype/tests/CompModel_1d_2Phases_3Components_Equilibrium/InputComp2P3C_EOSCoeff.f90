program main
  implicit none
  integer :: ndim1, ndim2
  real, dimension( :, : ), allocatable :: field
  integer :: idim1, idim2

  open( 5, file = 'filedim', status = 'unknown' )
  read( 5, * ) ndim1, ndim2
  close( 5 )

  allocate( field( ndim1, ndim2 ))

  field = 0
  field( 1 : ndim1, 1 ) = 1.

  do idim1 = 1, ndim1
     do idim2 = 1, ndim2
        print*, field( idim1, idim2 )
     end do
  end do

end program main
