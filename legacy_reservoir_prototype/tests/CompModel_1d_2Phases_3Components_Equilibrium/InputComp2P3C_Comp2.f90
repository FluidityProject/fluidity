
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

  ! Component 1, Phase 1:
  icomp = 1
  iphase = 1
     comp(( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1 : &
          ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods ) = 0.688172

  ! Component 1, Phase 2:
  iphase = 2
     comp(( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1 : &
          ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods ) = 0.215054

  ! Component 2, Phase 1:
  icomp = 2
  iphase = 1
     comp(( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1 : &
          ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods ) = 0.250000

  ! Component 2, Phase 2:
  iphase = 2
     comp(( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1 : &
          ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods ) = 0.166666

  ! Component 3, Phase 1:
  icomp = 3
  iphase = 1
     comp(( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1 : &
          ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods ) = 0.618280E-01

  ! Component 3, Phase 2:
  iphase = 2
     comp(( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1 : &
          ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods ) = 0.618280

  do i = 1, cv_nonods * nphase * ncomp
     print*, comp( i )
  end do

end program main
