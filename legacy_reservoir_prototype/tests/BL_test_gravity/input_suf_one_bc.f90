
program input_suf_one_bc
  implicit none
  integer, parameter :: stotel = 2, cv_snloc = 1, nphase = 1
  real, dimension( stotel * cv_snloc * nphase ) :: field
  integer :: i
  real :: value

  field = 0
  field( 1 ) = 1.
  field( 2 ) = 1.

  do i = 1, stotel * cv_snloc * nphase
    print*, field( i )
  end do

end program input_suf_one_bc
