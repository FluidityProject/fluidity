
program main
  implicit none
  integer :: cv_nonods, nphase, ncomp, dummy
  real, dimension( : ), allocatable :: comp
  ! Local variables
  integer :: icomp, iphase, i


  open( 5, file = 'filedim', status = 'unknown' )
  read( 5, * ) cv_nonods, nphase, ncomp, dummy
  close( 5 )

  allocate( comp( cv_nonods * nphase * ncomp ))

  comp = 0.

  ! Component 2, Phase 2: ( initially, just liquid brine )
  icomp = 2
  iphase = 2
     comp(( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1 : &
          ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods ) = 1.

  do i = 1, cv_nonods * nphase * ncomp
     print*, comp( i )
  end do

end program main
