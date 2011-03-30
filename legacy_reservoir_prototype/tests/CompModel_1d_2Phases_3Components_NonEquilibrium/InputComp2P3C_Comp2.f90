
program main
  implicit none
  integer :: cv_nonods, nphase, ncomp, dummy
  real, dimension( : ), allocatable :: comp
  ! Local variables
  integer :: icomp, iphase, i
  ! K_Comp based upon data on Table 6 of Analytical Methods in 
  ! Improved Oil Recovery -- Tara LaForce
  ! Assuming Comp_1(Phase Liq)


  open( 5, file = 'filedim', status = 'unknown' )
  read( 5, * ) cv_nonods, nphase, ncomp, dummy
  close( 5 )

  allocate( comp( cv_nonods * nphase * ncomp ))

  comp = 0.

  ! Component 1, Phase 2:
  icomp = 1
  iphase = 2
  comp(( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1 : &
       ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods ) = 0.10

  ! Component 2, Phase 2:
  icomp = 2
  iphase = 2
  comp(( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1 : &
       ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods ) = 0.25

  ! Component 3, Phase 2:
  icomp = 3
  iphase = 2
  comp(( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1 : &
       ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods ) = 0.65

  do i = 1, cv_nonods * nphase * ncomp
     print*, comp( i )
  end do

end program main
