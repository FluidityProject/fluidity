
program input_wic_t_bc
  implicit none
  integer, parameter :: stotel = 2, nphase = 1
  integer, dimension( stotel * nphase ) :: field
  integer :: i
  real :: value

  field = 0
  field( 1 ) = 1

  do i = 1, stotel * nphase
    print*, field( i )
  end do

end program input_wic_t_bc
