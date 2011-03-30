
program main
  implicit none
  integer :: ncomp, nphase
  real, dimension( :, :, : ), allocatable :: K_Comp
  integer :: iphase, jphase, icomp, dummy
  ! K_Comp based upon data on Table 6 of Analytical Methods in 
  ! Improved Oil Recovery -- Tara LaForce


  open( 5, file = 'filedim', status = 'unknown' )
  read( 5, * ) ncomp, nphase, dummy
  close( 5 )

  allocate( K_Comp( ncomp, nphase, nphase ))

  K_Comp = 0.

  ! Component 1:
  K_Comp( 1, 1, 2 ) = 3.2
  K_Comp( 1, 2, 1 ) = 3.2

  ! Component 2:
  K_Comp( 2, 1, 2 ) = 1.5
  K_Comp( 2, 2, 1 ) = 1.5

  ! Component 3:
  K_Comp( 3, 1, 2 ) = 0.1
  K_Comp( 3, 2, 1 ) = 0.1

  do icomp = 1, ncomp
     do iphase = 1, nphase
        do jphase = 1, nphase
           print*, K_Comp( icomp, iphase, jphase )
        end do
     end do
  end do


end program main
