
program input_suf_u_bc
  implicit none
  integer, parameter :: stotel = 2, u_snloc = 3, nphase = 2
  real, dimension( stotel * u_snloc * nphase ) :: field
  integer :: i
  real :: value

  field = 0.
  field( 1 ) = 1.

  do i = 1, stotel * nphase * u_snloc
    print*, field( i )
  end do

end program input_suf_u_bc
