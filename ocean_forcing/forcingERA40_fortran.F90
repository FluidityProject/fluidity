

! wrap up some of the coords module
subroutine rotate_wind(longitude, latitude, u, v, r3u, r3v, r3w)
  use Coordinates
  implicit none
  real, intent(in)::longitude, latitude, r3u, r3v, r3w
  real, intent(out)::u,v
  

  call rotate2ll(longitude, latitude, r3u, r3v, r3w, u, v)

end subroutine


