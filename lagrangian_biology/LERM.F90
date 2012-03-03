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
  use fields
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
  use state_module
  use iso_c_binding
  use detector_data_types
  use Profiler

implicit none

  private

  public :: LERM_get_diatom_fgroup, LERM_get_diatom_env_fields
  public :: LERM_initialise_living_diatom
  public :: LERM_update_living_diatom, LERM_update_dead_diatom

  interface

    subroutine lerm_updateLivingDiatom(vars, n_vars, env, n_env, dt) bind(c, name='updateLivingDiatom_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: n_vars, n_env
      real(c_double), dimension(n_vars), intent(inout) :: vars
      real(c_double), dimension(n_env), intent(inout) :: env
      real(c_double), intent(in) :: dt
    end subroutine lerm_updateLivingDiatom

    subroutine lerm_updateDeadDiatom(vars, n_vars, env, n_env, dt) bind(c, name='updateDeadDiatom_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: n_vars, n_env
      real(c_double), dimension(n_vars), intent(inout) :: vars
      real(c_double), dimension(n_env), intent(inout) :: env
      real(c_double), intent(in) :: dt
    end subroutine lerm_updateDeadDiatom

  end interface

contains

  subroutine LERM_update_living_diatom(agent, agent_list, state, dt)
    type(detector_type), pointer, intent(in) :: agent
    type(detector_linked_list), intent(in) :: agent_list
    type(state_type), intent(inout) :: state
    real, intent(in) :: dt

    real, dimension(size(agent_list%env_field_name)) :: env_field_values
    type(scalar_field), pointer :: env_field
    integer :: i

    call profiler_tic(trim(agent_list%name)//"::lerm_biology_update")

    do i=1, size(agent_list%env_field_name)
       env_field=>extract_scalar_field(state,trim(agent_list%env_field_name(i)))
       env_field_values(i)=eval_field(agent%element,env_field,agent%local_coords)
    end do

    call lerm_updateLivingDiatom(agent%biology, size(agent%biology), env_field_values, size(env_field_values), dt) 

    call profiler_toc(trim(agent_list%name)//"::lerm_biology_update")

  end subroutine LERM_update_living_diatom

  subroutine LERM_update_dead_diatom(agent, agent_list, state, dt)
    type(detector_type), pointer, intent(in) :: agent
    type(detector_linked_list), intent(in) :: agent_list
    type(state_type), intent(inout) :: state
    real, intent(in) :: dt

    real, dimension(size(agent_list%env_field_name)) :: env_field_values
    type(scalar_field), pointer :: env_field
    integer :: i

    call profiler_tic(trim(agent_list%name)//"::lerm_biology_update")

    do i=1, size(agent_list%env_field_name)
       env_field=>extract_scalar_field(state,trim(agent_list%env_field_name(i)))
       env_field_values(i)=eval_field(agent%element,env_field,agent%local_coords)
    end do

    call lerm_updateDeadDiatom(agent%biology, size(agent%biology), env_field_values, size(env_field_values), dt) 

    call profiler_toc(trim(agent_list%name)//"::lerm_biology_update")

  end subroutine LERM_update_dead_diatom

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

    allocate(fgroup%variables(13))
    fgroup%variables(1)%name = "Stage"
    fgroup%variables(1)%field_type = 0

    fgroup%variables(2)%name = "C_New"
    fgroup%variables(2)%field_type = 1
    fgroup%variables(2)%field_name = "DiatomEnsembleSize"
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
    fgroup%variables(6)%pool_index = 7

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
    fgroup%variables(8)%pool_index = 10

    fgroup%variables(9)%name = "Silicate_Rel"
    fgroup%variables(9)%field_type = 3
    fgroup%variables(9)%field_name = "DiatomReleaseSilicate"
    fgroup%variables(9)%field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"
    fgroup%variables(9)%chemfield = "DissolvedSilicate"
    fgroup%variables(9)%pool_index = 10

    fgroup%variables(10)%name = "Silicate_Pool"
    fgroup%variables(10)%field_type = 1
    fgroup%variables(10)%field_name = "DiatomParticulateSilicate"
    fgroup%variables(10)%field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"

    fgroup%variables(11)%name = "Ammonium_Ing"
    fgroup%variables(11)%field_type = 2
    fgroup%variables(11)%field_name = "DiatomRequestAmmonium"
    fgroup%variables(11)%field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"
    fgroup%variables(11)%depletion_field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"
    fgroup%variables(11)%chemfield = "DissolvedAmmonium"
    fgroup%variables(11)%pool_index = 13

    fgroup%variables(12)%name = "Ammonium_Rel"
    fgroup%variables(12)%field_type = 3
    fgroup%variables(12)%field_name = "DiatomReleaseAmmonium"
    fgroup%variables(12)%field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"
    fgroup%variables(12)%chemfield = "DissolvedAmmonium"
    fgroup%variables(12)%pool_index = 13

    fgroup%variables(13)%name = "Ammonium_Pool"
    fgroup%variables(13)%field_type = 1
    fgroup%variables(13)%field_name = "DiatomParticulateAmmonium"
    fgroup%variables(13)%field_path = trim(fg_path)//"/variables_lerm/scalar_field::Particulate"

  end subroutine LERM_get_diatom_fgroup

  subroutine LERM_get_diatom_env_fields(env_fields)
    character(len=FIELD_NAME_LEN), dimension(:), allocatable, intent(inout) :: env_fields

    allocate(env_fields(5))
    env_fields(1) = "DissolvedAmmonium"
    env_fields(2) = "DissolvedNitrate"
    env_fields(3) = "DissolvedSilicate"
    env_fields(4) = "Temperature"
    env_fields(5) = "Irradiance"

  end subroutine LERM_get_diatom_env_fields

  subroutine LERM_initialise_living_diatom(vars)
    real, dimension(13), intent(inout) :: vars

    vars(1) = 0.0;
    vars(2) = 50000.0;
    vars(3) = 50000.0;
    vars(4) = 1.5E-8;
    vars(5) = 2.7E-9;
    vars(6) = 0.0;
    vars(7) = 0.0;
    vars(8) = 0.0;
    vars(9) = 0.0;
    vars(10) = 1.05E-9;
    vars(11) = 0.0;
    vars(12) = 0.0;
    vars(13) = 2.2E-9;
  end subroutine LERM_initialise_living_diatom

end module LERM
