
program main
  implicit none
  integer :: stotel, cv_snloc, nphase, ncomp
  real, dimension( : ), allocatable ::  suf_comp_bc
  ! Local variables
  integer :: icomp, iphase, i
  ! Injection of components based upon data on Table 9 of Analytical Methods in 
  ! Improved Oil Recovery -- Tara LaForce


  open( 5, file = 'filedim', status = 'unknown' )
  read( 5, * ) stotel, cv_snloc, nphase, ncomp
  close( 5 )

  allocate( suf_comp_bc( stotel * cv_snloc * nphase * ncomp ))

  suf_comp_bc = 0.

  iphase = 1
  ! Component 1, Phase 1: ( injecting )
  icomp = 1
  suf_comp_bc( stotel * cv_snloc * nphase * ( icomp - 1 ) + &
       stotel * cv_snloc * (iphase - 1 ) + 1 ) = 0.787

  ! Component 2, Phase 1: ( injecting )
  icomp = 2
  suf_comp_bc( stotel * cv_snloc * nphase * ( icomp - 1 ) + &
       stotel * cv_snloc * (iphase - 1 ) + 1 ) = 0.213


  do i = 1, stotel * cv_snloc * nphase * ncomp
     print*, suf_comp_bc( i )
  end do

end program main
