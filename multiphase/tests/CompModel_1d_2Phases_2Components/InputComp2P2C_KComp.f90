program main
  implicit none
  integer :: ndim1, ndim2, ndim3
  real, dimension( :, :, : ), allocatable :: field
  integer :: idim1, idim2, idim3

  open( 5, file = 'filedim', status = 'unknown' )
  read( 5, * ) ndim1, ndim2, ndim3
  close( 5 )

  allocate( field( ndim1, ndim2, ndim3 ))
  field = 0.

  do idim1 = 1, ndim1
     do idim2 = 1, ndim2
        do idim3 = 1, ndim3
           if( idim1 == 1 )then ! CO2
              if( idim2 /= idim3 ) field( idim1, idim2, idim3 ) = 1. / 4.5
           elseif( idim1 == 2 )then ! Brine
              if( idim2 /= idim3 ) field( idim1, idim2, idim3 ) = 1. / 0.2     
           end if
        end do
     end do
  end do

  do idim1 = 1, ndim1
     do idim2 = 1, ndim2
        do idim3 = 1, ndim3
           print*, field( idim1, idim2, idim3 )
        end do
     end do
  end do

end program main
