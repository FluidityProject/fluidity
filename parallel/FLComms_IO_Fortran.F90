!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    C.Pain@Imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module flcomms_io

  use fldebug

  implicit none
  
  private
  
  public :: read_halos, write_halos
  
  interface
    function cread_halos(basename, basename_len)
      integer, intent(in) :: basename_len
      character(len = basename_len), intent(in) :: basename
      integer :: cread_halos
    end function cread_halos
  
    function cwrite_halos(basename, basename_len)
      integer, intent(in) :: basename_len
      character(len = basename_len), intent(in) :: basename
      integer :: cwrite_halos
    end function cwrite_halos
  end interface
  
contains

  subroutine read_halos(basename, stat)
    !!< Read the halos from the supplied halo file with the supplied base name.
    !!< The file read from will be:
    !!<   basename_[process number].halo

    character(len = *), intent(in) :: basename
    integer, optional, intent(out) :: stat
    
    integer :: lstat
    
    ewrite(1, *) "In read_halos"
    
    if(present(stat)) then
      stat = 0
    end if
    
    lstat = cread_halos(trim(basename), len_trim(basename))
    if(lstat /= 0) then
      if(present(stat)) then
        stat = lstat
      else
        FLAbort("Failed to read halos from file with base name " // trim(basename))
      end if
      ewrite(1, *) "Exiting read_halos after an error"
      return
    end if
    
    ewrite(1, *) "Exiting read_halos"
    
  end subroutine read_halos

  subroutine write_halos(basename, stat)
    !!< Write out the current halos to a halo file with the supplied base name.
    !!< The file written to will be:
    !!<   basename_[process number].halo
    
    character(len = *), intent(in) :: basename
    integer, optional, intent(out) :: stat
    
    integer :: lstat
    
    ewrite(1, *) "In write_halos"
    
    if(present(stat)) then
      stat = 0
    end if
    
    lstat = cwrite_halos(trim(basename), len_trim(basename))
    if(lstat /= 0) then
      if(present(stat)) then
        stat = lstat
      else
        FLAbort("Failed to write out halos to file with base name " // trim(basename))
      end if
      ewrite(1, *) "Exiting write_halos after an error"
      return
    end if
    
    ewrite(1, *) "Exiting write_halos"
    
  end subroutine write_halos

end module flcomms_io
