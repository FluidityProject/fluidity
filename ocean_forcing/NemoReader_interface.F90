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

module nemo_v2
  use FLDebug
  use global_parameters
  use, intrinsic :: iso_c_binding
  
  implicit none
  
  interface
     subroutine nemo_v2_AddFieldOfInterest_c(scalar) bind(c)
       use, intrinsic :: iso_c_binding
       character(c_char), intent(in) :: scalar
     end subroutine nemo_v2_AddFieldOfInterest_c

     subroutine nemo_v2_ClearFields_c() bind(c)
       use, intrinsic :: iso_c_binding
     end subroutine nemo_v2_ClearFields_c

     subroutine nemo_v2_GetScalars_c(longitude, latitude, p_depth, scalars) bind(c)
       use, intrinsic :: iso_c_binding
       real(c_double), intent(in) :: longitude, latitude, p_depth
       real(c_double), dimension(*), intent(out) :: scalars
     end subroutine nemo_v2_GetScalars_c
     
     subroutine nemo_v2_RegisterDataFile_c(filename) bind(c)
       use, intrinsic :: iso_c_binding
       character(c_char), intent(in) :: filename
     end subroutine nemo_v2_RegisterDataFile_c
     
     subroutine nemo_v2_SetSimulationTimeUnits_c(units) bind(c)
       use, intrinsic :: iso_c_binding
       character(c_char), intent(in) :: units
     end subroutine nemo_v2_SetSimulationTimeUnits_c
     
     subroutine nemo_v2_SetTimeSeconds_c(time) bind(c)
       use, intrinsic :: iso_c_binding
       real(c_double), intent(in) :: time
     end subroutine nemo_v2_SetTimeSeconds_c

     subroutine get_nemo_variables_c(time, X, Y, Z, DEPTH, Te, Sa, U, V, W, SSH, NNodes) bind(c)
       use, intrinsic :: iso_c_binding
       real(c_double), intent(in) :: time
       real(c_double), dimension(*), intent(in) :: X, Y, Z, DEPTH
       real(c_double), dimension(*), intent(out) :: Te, Sa, U, V, W, SSH
       integer(c_int), intent(in) :: NNodes
     end subroutine get_nemo_variables_c
  end interface
  
contains
  
  subroutine nemo_v2_AddFieldOfInterest(scalar)
    character(len=*), intent(in) :: scalar
    call nemo_v2_AddFieldOfInterest_c(trim(scalar)//c_null_char)
  end subroutine nemo_v2_AddFieldOfInterest
  
  subroutine nemo_v2_ClearFields()
    call nemo_v2_ClearFields_c()
  end subroutine nemo_v2_ClearFields
  
  subroutine nemo_v2_GetScalars(longitude, latitude, p_depth, scalars)
    real, intent(in) :: longitude, latitude, p_depth
    real, dimension(:), intent(out) :: scalars
    call nemo_v2_GetScalars_c(longitude, latitude, p_depth, scalars)
  end subroutine nemo_v2_GetScalars
  
  subroutine nemo_v2_RegisterDataFile(filename)
    character(len=*), intent(in) :: filename
    call nemo_v2_RegisterDataFile_c(trim(filename)//c_null_char)
  end subroutine nemo_v2_RegisterDataFile
  
  subroutine nemo_v2_SetSimulationTimeunits(units)
    character(len=*), intent(in) :: units
    call nemo_v2_SetSimulationTimeunits_c(trim(units)//c_null_char)
  end subroutine nemo_v2_SetSimulationTimeunits
  
  subroutine nemo_v2_SetTimeSeconds(time)
    real(real_8), intent(in) :: time
    call nemo_v2_SetTimeSeconds_c(time)
  end subroutine nemo_v2_SetTimeSeconds
  
  subroutine get_nemo_variables(time, X, Y, Z, DEPTH, Te, Sa, U, V, W, SSH, NNodes)
    use, intrinsic :: iso_c_binding
    real(real_8), intent(in) :: time
    real, dimension(:), intent(in) :: X, Y, Z, DEPTH
    real, dimension(:), intent(out) :: Te, Sa, U, V, W, SSH
    integer, intent(in) :: NNodes
    call get_nemo_variables_c(time, X, Y, Z, DEPTH, Te, Sa, U, V, W, SSH, NNodes)
  end subroutine get_nemo_variables
  
end module nemo_v2
