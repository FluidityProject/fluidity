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

module fluxes
  use FLDebug
  use global_parameters
  use, intrinsic :: iso_c_binding
  
  implicit none
  
  interface
     subroutine fluxes_AddFieldOfInterest_c(scalar) bind(c)
       use, intrinsic :: iso_c_binding
       character(c_char), intent(in) :: scalar
     end subroutine fluxes_AddFieldOfInterest_c

     subroutine fluxes_ClearFields_c() bind(c)
       use, intrinsic :: iso_c_binding
     end subroutine fluxes_ClearFields_c

     subroutine fluxes_GetScalars_c(longitude, latitude, scalars) bind(c)
       use, intrinsic :: iso_c_binding
       real(c_double), intent(in) :: longitude
       real(c_double), intent(in) :: latitude
       real(c_double), dimension(*), intent(out) :: scalars
     end subroutine fluxes_GetScalars_c
     
     subroutine fluxes_GetScalar_c(name, longitude, latitude, scalar) bind(c)
       use, intrinsic :: iso_c_binding
       character(c_char), intent(in) :: name
       real(c_double), intent(in) :: longitude
       real(c_double), intent(in) :: latitude
       real(c_double), intent(out) :: scalar
     end subroutine fluxes_GetScalar_c
     
     subroutine fluxes_RegisterDataFile_c(filename) bind(c)
       use, intrinsic :: iso_c_binding
       character(c_char), intent(in) :: filename
     end subroutine fluxes_RegisterDataFile_c
     
     subroutine fluxes_SetSimulationTimeUnits_c(units) bind(c)
       use, intrinsic :: iso_c_binding
       character(c_char), intent(in) :: units
     end subroutine fluxes_SetSimulationTimeUnits_c
     
     subroutine fluxes_SetTimeSeconds_c(time) bind(c)
       use, intrinsic :: iso_c_binding
       real(c_double), intent(in) :: time
     end subroutine fluxes_SetTimeSeconds_c
  end interface
  
contains
  
  subroutine fluxes_AddFieldOfInterest(scalar)
    character(len=*), intent(in) :: scalar
    call fluxes_AddFieldOfInterest_c(trim(scalar)//c_null_char)
  end subroutine fluxes_AddFieldOfInterest
  
  subroutine fluxes_ClearFields()
    call fluxes_ClearFields_c()
  end subroutine fluxes_ClearFields
  
  subroutine fluxes_GetScalars(longitude, latitude, scalars)
    real, intent(in) :: longitude, latitude
    real, dimension(:), intent(out) :: scalars
    call fluxes_GetScalars_c(longitude, latitude, scalars)
  end subroutine fluxes_GetScalars
  
  subroutine fluxes_GetScalar(name, longitude, latitude, scalar)
    character(len=*), intent(in) :: name
    real, intent(in) :: longitude, latitude
    real, intent(out) :: scalar
    call fluxes_GetScalar_c(trim(name)//c_null_char, longitude, latitude, scalar)
  end subroutine fluxes_GetScalar
  
  subroutine fluxes_RegisterDataFile(filename)
    character(len=*), intent(in) :: filename
    call fluxes_RegisterDataFile_c(trim(filename)//c_null_char)
  end subroutine fluxes_RegisterDataFile
  
  subroutine fluxes_SetSimulationTimeunits(units)
    character(len=*), intent(in) :: units
    call fluxes_SetSimulationTimeunits_c(trim(units)//c_null_char)
  end subroutine fluxes_SetSimulationTimeunits
  
  subroutine fluxes_SetTimeSeconds(time)
    real(real_8), intent(in) :: time
    call fluxes_SetTimeSeconds_c(time)
  end subroutine fluxes_SetTimeSeconds
  
end module fluxes
