program suf_t_bc_funct
  implicit none
  ! Allocating  SUF_T_BC( STOTEL * CV_SNLOC * NPHASE )
  integer, parameter :: stotel= 2, cv_snloc = 1, nphase = 2
  real, dimension( stotel * cv_snloc * nphase ) :: suf_t_bc
  ! Local
  integer :: iphase, i

  suf_t_bc( 1 ) = 0.5
  suf_t_bc( 2 ) = 0.
  do iphase = 2, nphase
     suf_t_bc( 1 + ( iphase - 1 ) * 2 ) = 0.5
     suf_t_bc( 2 + ( iphase - 1 ) * 2 ) = 0.
  end do

  print*, stotel * cv_snloc * nphase
  do i = 1, stotel * cv_snloc * nphase
     print*, suf_t_bc( i )
  end do

  return
end program suf_t_bc_funct

