!    Copyright (C) 2008 Imperial College London and others.
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

module LERM
  use fldebug
  use spud
  use global_parameters, only: PYTHON_FUNC_LEN, OPTION_PATH_LEN
  use detector_data_types

implicit none

  private

  public :: LERM_get_diatom_fgroup, LERM_update_living_diatom

  interface

    subroutine LERM_update_living_diatom(vars, env, dt) bind(c, name='updateLivingDiatom_c')
      use :: iso_c_binding
      implicit none
      real(c_float), dimension(11), intent(inout) :: vars
      real(c_float), dimension(5), intent(inout) :: env
      real(c_float), intent(in) :: dt
    end subroutine LERM_update_living_diatom

  end interface

contains

  subroutine LERM_get_diatom_fgroup(fgroup, fg_path)
    type(functional_group), intent(inout) :: fgroup
    character(len=*), intent(in) :: fg_path

    character(len=OPTION_PATH_LEN) :: stage_buffer
    integer :: i, stage_count

    ! Record FG name and Agents field path
    call get_option(trim(fg_path)//"/name", fgroup%name)
    if (have_option(trim(fg_path)//"/scalar_field::Agents")) then
       fgroup%agents_field_path = trim(fg_path)//"/scalar_field::Agents"
    else
       FLExit("No Agents field specified for functional group "//trim(fgroup%name))
    end if

    ! Record all associated stage names
    stage_count = option_count(trim(fg_path)//"/stages/stage")
    if (stage_count>0) then
       allocate(fgroup%stage_names%ptr(stage_count))
       do i=1, stage_count
          write(stage_buffer, "(a,i0,a)") trim(fg_path)//"/stages/stage[",i-1,"]"
          call get_option(trim(stage_buffer)//"/name", fgroup%stage_names%ptr(i))
       end do
    end if

    allocate(fgroup%variables(11))
    fgroup%variables(1)%name = "Stage"
    fgroup%variables(1)%field_type = 0

    fgroup%variables(2)%name = "C_New"
    fgroup%variables(2)%field_type = 1
    fgroup%variables(2)%field_name = "DiatomBiomass"
    fgroup%variables(2)%field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"

    fgroup%variables(3)%name = "C_Old"
    fgroup%variables(3)%field_type = 0

    fgroup%variables(4)%name = "Carbon_Pool"
    fgroup%variables(4)%field_type = 1
    fgroup%variables(4)%field_name = "DiatomParticulateCarbon"
    fgroup%variables(4)%field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"

    fgroup%variables(5)%name = "Chlorophyll_Pool"
    fgroup%variables(5)%field_type = 1
    fgroup%variables(5)%field_name = "DiatomParticulateChlorophyll"
    fgroup%variables(5)%field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"

    fgroup%variables(6)%name = "Nitrate_Ing"
    fgroup%variables(6)%field_type = 2
    fgroup%variables(6)%field_name = "DiatomRequestNitrate"
    fgroup%variables(6)%field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"
    fgroup%variables(6)%depletion_field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"
    fgroup%variables(6)%chemfield = "DissolvedNitrate"

    fgroup%variables(7)%name = "Nitrate_Pool"
    fgroup%variables(7)%field_type = 1
    fgroup%variables(7)%field_name = "DiatomParticulateNitrate"
    fgroup%variables(7)%field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"

    fgroup%variables(8)%name = "Silicate_Ing"
    fgroup%variables(8)%field_type = 2
    fgroup%variables(8)%field_name = "DiatomRequestSilicate"
    fgroup%variables(8)%field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"
    fgroup%variables(8)%depletion_field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"
    fgroup%variables(8)%chemfield = "DissolvedSilicate"

    fgroup%variables(9)%name = "Silicate_Pool"
    fgroup%variables(9)%field_type = 1
    fgroup%variables(9)%field_name = "DiatomParticulateSilicate"
    fgroup%variables(9)%field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"

    fgroup%variables(10)%name = "Ammonium_Ing"
    fgroup%variables(10)%field_type = 2
    fgroup%variables(10)%field_name = "DiatomRequestAmmonium"
    fgroup%variables(10)%field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"
    fgroup%variables(10)%depletion_field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"
    fgroup%variables(10)%chemfield = "DissolvedAmmonium"

    fgroup%variables(11)%name = "Ammonium_Pool"
    fgroup%variables(11)%field_type = 1
    fgroup%variables(11)%field_name = "DiatomParticulateAmmonium"
    fgroup%variables(11)%field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"

  end subroutine LERM_get_diatom_fgroup

end module LERM
