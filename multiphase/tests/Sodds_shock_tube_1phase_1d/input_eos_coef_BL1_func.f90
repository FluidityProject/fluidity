
program input_eos_coef_funct
  implicit none
  ! Allocating EOS_COEFS( NPHASE, NCOEF )
  integer, parameter :: nphase = 2, ncoef = 10
  real, dimension( nphase, ncoef ) :: eos_coefs
  ! Local
  integer :: i, j

 ! print*, nphase, ncoef
  eos_coefs = 0.
  eos_coefs( 1 : nphase, 1 ) = 1.
  eos_coefs( 1, 2 ) = 0.01
 
  do i = 1, nphase
     do j = 1, ncoef
        print*, eos_coefs( i, j )
     end do
  end do

  return
end program input_eos_coef_funct
