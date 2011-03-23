
program input_satura_funct
  implicit none
  ! Allocating  SATURA( CV_PHA_NONODS )
  integer, parameter :: cv_nloc = 3, totele = 10, nphase = 2
  integer, parameter :: cv_nonods = ( cv_nloc - 1 ) * totele + 1, cv_pha_nonods = cv_nonods * nphase
  real, dimension( cv_pha_nonods ) :: satura
  ! Local
  integer :: i

 ! print*, cv_pha_nonods
  do i = 1, cv_pha_nonods
     if( i <= cv_nonods ) then
        satura( i ) = 0.
     else
        satura( i ) = 1.
     endif
     print*, satura( i )
  end do

  return
end program input_satura_funct

