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

module climatology
  use FLDebug
  use global_parameters
  use, intrinsic :: iso_c_binding
  
  implicit none
  
  interface
     subroutine climatology_SetSimulationTimeunits_c(units) bind(c)
       use, intrinsic :: iso_c_binding
       character(c_char), intent(in) :: units
     end subroutine climatology_SetSimulationTimeunits_c

     subroutine climatology_SetClimatology_c(filename) bind(c)
       use, intrinsic :: iso_c_binding
       character(c_char), intent(in) :: filename
     end subroutine climatology_SetClimatology_c
     
     subroutine climatology_SetTimeSeconds_c(time) bind(c)
       use, intrinsic :: iso_c_binding
       real(c_double), intent(in) :: time
     end subroutine climatology_SetTimeSeconds_c
     
     subroutine climatology_GetValue_c(name, x, y, z, value) bind(c)
       use, intrinsic :: iso_c_binding
       character(c_char), intent(in) :: name
       real(c_double), intent(in) :: x, y, z
       real(c_double), intent(out) :: value
     end subroutine climatology_GetValue_c
     
     subroutine climatology_GetSurfaceValue_c(name, x, y, z, value) bind(c)
       use, intrinsic :: iso_c_binding
       character(c_char), intent(in) :: name
       real(c_double), intent(in) :: x, y, z
       real(c_double), intent(out) :: value
     end subroutine climatology_GetSurfaceValue_c
  end interface
  
contains
  
  subroutine climatology_SetSimulationTimeunits(units)
    character(len=*), intent(in) :: units
    call climatology_SetSimulationTimeunits_c(trim(units)//c_null_char)
  end subroutine climatology_SetSimulationTimeunits
  
  subroutine climatology_SetClimatology(filename)
    character(len=*), intent(in) :: filename
    call climatology_SetClimatology_c(trim(filename)//c_null_char)
  end subroutine climatology_SetClimatology

  subroutine climatology_SetTimeSeconds(time)
    real(real_8), intent(in) :: time
    call climatology_SetTimeSeconds_c(time)
  end subroutine climatology_SetTimeSeconds

  subroutine climatology_GetValue(name, x, y, z, value)
    character(len=*), intent(in) :: name
    real, intent(in) :: x, y, z
    real, intent(out) :: value
    call climatology_GetValue_c(trim(name)//c_null_char, x, y, z, value)
  end subroutine climatology_GetValue

  subroutine climatology_GetSurfaceValue(name, xyz, value)
    character(len=*), intent(in) :: name
    real, dimension(1:3), intent(in) :: xyz
    real, intent(out) :: value
    call climatology_GetSurfaceValue_c(trim(name)//c_null_char, xyz(1), xyz(2), xyz(3), value)
  end subroutine climatology_GetSurfaceValue
  
end module climatology
