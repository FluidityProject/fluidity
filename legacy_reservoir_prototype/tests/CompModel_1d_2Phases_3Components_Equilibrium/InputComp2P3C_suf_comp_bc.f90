
program main
  implicit none
  integer :: stotel, cv_snloc, nphase, ncomp
  real, dimension( : ), allocatable ::  suf_comp_bc
  ! Local variables
  integer :: icomp, iphase, i


  open( 5, file = 'filedim', status = 'unknown' )
  read( 5, * ) stotel, cv_snloc, nphase, ncomp
  close( 5 )

  allocate( suf_comp_bc( stotel * cv_snloc * nphase * ncomp ))

  suf_comp_bc = 0.

  ! Assuming all components are injected under
  ! equilibrium conditions

  ! Component 1, Phase 1:
  icomp = 1
  iphase = 1
  suf_comp_bc( stotel * cv_snloc * nphase * ( icomp - 1 ) + &
       stotel * cv_snloc * (iphase - 1 ) + 1 ) = 0.688172

  ! Component 1, Phase 2:
  iphase = 2
  suf_comp_bc( stotel * cv_snloc * nphase * ( icomp - 1 ) + &
       stotel * cv_snloc * (iphase - 1 ) + 1 ) = 0.215054

  ! Component 2, Phase 1:
  icomp = 2
  iphase = 1
  suf_comp_bc( stotel * cv_snloc * nphase * ( icomp - 1 ) + &
       stotel * cv_snloc * (iphase - 1 ) + 1 ) = 0.250000

  ! Component 2, Phase 2:
  iphase = 2
  suf_comp_bc( stotel * cv_snloc * nphase * ( icomp - 1 ) + &
       stotel * cv_snloc * (iphase - 1 ) + 1 ) = 0.166667

  ! Component 3, Phase 1:
  icomp = 3
  iphase = 1
  suf_comp_bc( stotel * cv_snloc * nphase * ( icomp - 1 ) + &
       stotel * cv_snloc * (iphase - 1 ) + 1 ) = 0.618280E-01

  ! Component 3, Phase 2:
  iphase = 2
  suf_comp_bc( stotel * cv_snloc * nphase * ( icomp - 1 ) + &
       stotel * cv_snloc * (iphase - 1 ) + 1 ) = 0.618280


  do i = 1, stotel * cv_snloc * nphase * ncomp
     print*, suf_comp_bc( i )
  end do

end program main
