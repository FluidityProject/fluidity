
program input_perm_funct
  implicit none
  ! Allocating  PERM( TOTELE, NDIM, NDIM )
  integer, parameter :: totele = 50, ndim = 1
  real, dimension( totele, ndim, ndim ) :: perm
  ! Local
  integer :: i, j, k
  real, parameter :: perm_value = 1.

  print*, totele, ndim, ndim

  perm = 0.
  perm( 1 : totele, 1, 1 ) = perm_value

  do i = 1, totele
     do j = 1, ndim
        do k = 1, ndim
           print*, perm( i, j, k )
        end do
     end do
  end do

  return
end program input_perm_funct
