! Copyright (C) 2006 Imperial College London and others.
!   
! Please see the AUTHORS file in the main source directory for a full list
! of copyright holders.
!
! Prof. C Pain
! Applied Modelling and Computation Group
! Department of Earth Science and Engineering
! Imperial College London
!
! amcgsoftware@imperial.ac.uk
!  
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
! USA

module SampleNetCDF
  use FLDebug
  implicit none
contains
  
  subroutine SampleNetCDF_Open(name, id)
    character(len=*), intent(in)::name
    integer, intent(out)::id
    
    call samplenetcdf_open_c(name, len_trim(name), id)
  end subroutine SampleNetCDF_Open

  subroutine SampleNetCDF_SetVariable(id, varname)
    character(len=*), intent(in)::varname
    integer, intent(out)::id
    
    call samplenetcdf_setvariable_c(id, varname, len_trim(varname))
  end subroutine SampleNetCDF_SetVariable

  subroutine SampleNetCDF_GetValue(id, longitude, latitude, val)
    integer, intent(out)::id
    real, intent(out)::longitude, latitude
    real, intent(out)::val

    call samplenetcdf_getvalue_c(id, longitude, latitude, val)
  end subroutine SampleNetCDF_GetValue

  subroutine SampleNetCDF_Close(id)
    integer, intent(out)::id
    
    call samplenetcdf_close_c(id)
  end subroutine SampleNetCDF_Close
  
end module SampleNetCDF
