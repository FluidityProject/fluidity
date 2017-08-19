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
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
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
module cv_fields
  !!< Module containing general tools for discretising Control Volume problems.
  use spud
  use fldebug
  use global_parameters, only: FIELD_NAME_LEN
  use fields
  use state_module
  use cvtools
  use diagnostic_fields

  implicit none

  private
  public :: cv_disc_get_cfl_no

contains

  subroutine cv_disc_get_cfl_no(field_option_path, state, mesh, cfl_no, density_option_path)
    !!< This subroutine returns the cfl number (if needed) 
    !!< based on the option paths provided.
    !!< If it is not required it returns a constant field.
    !!< this subroutine allocates a field... please deallocate after!

    !! option paths used to determine if a cfl no is needed for field
    character(len=*), dimension(:), intent(in) :: field_option_path
    ! bucket full of fields
    type(state_type), intent(inout) :: state
    ! mesh to allocate the cfl_no on
    type(mesh_type), intent(inout) :: mesh
    !! a scalar field for the cfl no.... still to be allocated
    type(scalar_field), intent(inout) :: cfl_no
    !! option path used to determine if a cfl no is needed for density
    character(len=*), dimension(:), optional, intent(in) :: density_option_path

    character(len=FIELD_NAME_LEN) :: cfl_type
    integer :: f, cfl_stat, stat
    ! somewhere to put strings temporarily
    character(len=FIELD_NAME_LEN) :: tmpstring
    integer :: nfields

    nfields = size(field_option_path)
    if(present(density_option_path)) then
      ! this allows multiple grouped fields to be evaluated at once 
      ! (mostly used for coupled_cv)
      assert(nfields==size(density_option_path))
    end if

    ! hmmm, do we need a cfl no.?
    ! there are three reasons we might...
    ! the field or density discretisation might need it
    ! or we might need it to subcycle
    ! at the moment these all have to use the same definition
    ! of the courant number
    cfl_type="start"
    cfl_stat=1
    
    do f = 1, nfields
      ! check to see if the field discretisation requires a courant number
      call get_option(trim(complete_cv_field_path(field_option_path(f)))//&
                            "/face_value[0]/courant_number[0]/name", &
                            tmpstring, stat)
      if(stat==0) then
        if(trim(cfl_type)=="start") then
          cfl_type=tmpstring
          cfl_stat=stat
        elseif(trim(cfl_type)/=trim(tmpstring)) then
          ewrite(-1,*) "Attempting to discretise two fields using different courant numbers."
          FLExit("This is not currently supported.")
        end if
      end if

      if(present(density_option_path)) then
        ! check to see if the density discretisation requires a courant number
        call get_option(trim(complete_cv_field_path(density_option_path(f)))//&
                              "/face_value[0]/courant_number[0]/name", &
                              tmpstring, stat)
        if(stat==0) then
          if(trim(cfl_type)=="start") then
            cfl_type=tmpstring
            cfl_stat=stat
          elseif(trim(cfl_type)/=trim(tmpstring)) then
            ewrite(-1,*) "Attempting to discretise two fields using different courant numbers."
            FLExit("This is not currently supported.")
          end if
        end if
      end if

      ! check to see if we need the courant number to subcycle with
      call get_option(trim(field_option_path(f))//&
                            "/prognostic/temporal_discretisation&
                            &/control_volumes/maximum_courant_number_per_subcycle&
                            &/courant_number[0]/name", &
                            tmpstring, stat)
      if(stat==0) then
        if(trim(cfl_type)=="start") then
          cfl_type=tmpstring
          cfl_stat=stat
        elseif(trim(cfl_type)/=trim(tmpstring)) then
          ewrite(-1,*) "Attempting to discretise face values "//&
                      "using a "//trim(cfl_type)//" courant number"
          ewrite(-1,*) "and to subcycle "//&
                      "using a "//trim(tmpstring)//" courant number."
          FLExit("This is not currently supported.")
        end if
      end if

      ! check to see if we need the courant number for the limiter
      call get_option(trim(complete_cv_field_path(field_option_path(f)))//&
                            "/face_value[0]/limit_face_value/limiter[0]&
                            &/courant_number[0]/name", &
                            tmpstring, stat)
      if(stat==0) then
        if(trim(cfl_type)=="start") then
          cfl_type=tmpstring
          cfl_stat=stat
        elseif(trim(cfl_type)/=trim(tmpstring)) then
          ewrite(-1,*) "Attempting to discretise face values or subcycle"//&
                      "using a "//trim(cfl_type)//" courant number"
          ewrite(-1,*) "and to limit "//&
                      "using a "//trim(tmpstring)//" courant number."
          FLExit("This is not currently supported.")
        end if
      end if
      
    end do

    if (cfl_stat==0) then
      ! otherwise we want to calculate a node centred field of the cfl number
      call allocate(cfl_no, mesh, "CourantNumber")
      call calculate_diagnostic_variable(state, trim(cfl_type), cfl_no, &
          &option_path=trim(complete_cv_field_path(field_option_path(1)))//"/face_value[0]/courant_number[0]")
    else
      ! if we don't need a cfl number then just set it all to 1
      call allocate(cfl_no, mesh, "CourantNumber", field_type=FIELD_TYPE_CONSTANT)
      call set(cfl_no, 1.0)
    end if

  end subroutine cv_disc_get_cfl_no

end module cv_fields
