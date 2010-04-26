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
!    C.Pain@Imperial.ac.uk
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

module turbine

  use spud
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use futils

  implicit none

  private
  public :: turbine_check_options

contains

  subroutine turbine_check_options
    character(len=OPTION_PATH_LEN):: turbine_path, turbine_name, bc_name
    integer :: notur, i, j
    logical :: have_dirichlet_model, have_dgflux_model

    ! Don't check turbine configuration if it's not included in the model!
    if (.not.have_option("/turbine_model")) return

    have_dirichlet_model=.false.
    have_dgflux_model=.false.

    ! loop through turbines
    notur = option_count("/turbine_model/turbine")
    do i=0, notur-1
       turbine_path="/turbine_model/turbine["//int2str(i)//"]"
       if (have_option(trim(turbine_path)//"/dirichlet")) then
           have_dirichlet_model=.true.
           ! The specified b.c.'s  in the turbine model must be dirichlet boundary conditions with normal_component.
           do j=1,2
              call get_option("/turbine_model/turbine["//int2str(i)//"]/dirichlet/boundary_condition_name_"//int2str(j)//"/name", bc_name)
              if (.not. have_option("/material_phase[0]/vector_field::Velocity/prognostic/boundary_conditions::"//trim(bc_name)//"/type::dirichlet/align_bc_with_surface/normal_component")) then
                 call get_option("/turbine_model/turbine["//int2str(i)//"]/name", turbine_name)
                 FLExit("Error while checking the options for turbine '"//trim(turbine_name)//"': Turbine model boundary has to be dirichlet boundary conditions with a normal_component.")
              end if
           end do
       elseif (have_option(trim(turbine_path)//"/dgflux")) then
            have_dgflux_model=.true.
            FLExit("DGFlux turbine not supported yet. Use Dirichlet turbine")
       else
          FLAbort("Unknown turbine model specified!")
       end if
    end do

    if (have_dirichlet_model) then
       ! We need the FreeSurface field.
       if (.not.have_option("/material_phase[0]/scalar_field::FreeSurface")) then
          FLExit("Turbine modelling requires FreeSurface to be activated.")
       end if
    end if

  end subroutine turbine_check_options

end module turbine
