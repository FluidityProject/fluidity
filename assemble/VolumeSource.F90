#include "fdebug.h"

module VolumeSource
  use fldebug
  use Futils
  use state_module
  use Fields
  implicit none

  logical :: initialised = .false.
  logical :: have_volumesource = .false.
  real :: originX,originY,Radius,Bottom,horiz_ramp,vert_ramp,strength
  logical :: GOT_volumesource_IT
  integer :: IT_volumesource
contains  

  subroutine Volumesource_initialise()
    integer :: unit, io
    namelist/volumeflux/originx,originy,radius,bottom,horiz_ramp,vert_ramp,strength

    if(.not.initialised) then
       ewrite(2,*) 'Reading in file'
       unit=free_unit()
       open(unit=unit, file="volumesource.dat", action="read",&
            & status="old", iostat=io)
       if(io == 0) then
          ewrite(2,*) 'Reading data'
          read(unit, nml=volumeflux)
          initialised = .true.
          have_volumesource = .true.
       end if
       initialised = .true.
       ewrite(2,*) 'volumeflux',originx,originy,radius,bottom,horiz_ramp,vert_ramp,strength
    end if

    ewrite(2,*) 'have_volumesource',have_volumesource
    
  end subroutine Volumesource_initialise

  subroutine set_volumesource(X,Y,Z,sourceT,state)
    real, dimension(:), intent(in) :: X,Y,Z
    real, dimension(:), intent(out) :: sourceT
    type(state_type), intent(inout) :: state

    integer :: I
    real :: R
    real :: source
    integer :: stat
    type(scalar_field), pointer :: source_field

    if(size(X).ne.size(sourceT)) then
       FLAbort('Only works for xnonod == nonods at present')
    end if

    if(have_volumesource) then
       do I = 1, size(X)
          source = 0.
          if(Z(I)>bottom-vert_ramp) then
             R = sqrt((X(I)-originX)**2 + (Y(I)-originY)**2)
             if(R<radius+horiz_ramp) then
                source = strength
                if(R>radius) then
                   source = source*(1.0 - (R-radius)/horiz_ramp)
                end if
                if(Z(I)<bottom) then
                   source = source*(1.0 - (bottom-Z(I))/vert_ramp)
                end if
             end if
          end if
          sourceT(I) = source
       end do

       !copy into field
       source_field => extract_scalar_field(state,'VolumeTemperatureSource',stat)
       if(stat==0) then
          source_field%val = sourceT
       else
          FLAbort("couldnt find field")
       end if
    end if

  end subroutine set_volumesource

end module VolumeSource
