!    Copyright (C) 2007 Imperial College London and others.
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
module initialise_ocean_forcing_module
  use FLDebug
  use spud
  use global_parameters, only: OPTION_PATH_LEN
  use climatology
  use fluxes
  use nemo_v2
  use nemo_states_module

  implicit none

  private

  public initialise_ocean_forcing_readers

contains

  subroutine initialise_ocean_forcing_readers
    !! Sets up the readers for ocean conditions (climatology, boundary
    !! conditions, and forcing).

    character(len=OPTION_PATH_LEN) :: option, dataset

    ewrite(1,*) "In initialise_ocean_forcing_readers"

    if(have_option("/timestepping/current_time/time_units")) then
       call get_option("/timestepping/current_time/time_units/date", option)

       call fluxes_setsimulationtimeunits(trim(option))
       call climatology_setsimulationtimeunits(trim(option))
       call NEMO_v2_setsimulationtimeunits(trim(option))
    end if

    if(have_option("/environmental_data/climatology/file_name")) then
       call get_option("/environmental_data/climatology/file_name", option)

       call climatology_setclimatology(trim(option))
    end if

    if(have_option("/ocean_forcing/bulk_formulae/input_file")) then
       call get_option("/ocean_forcing/bulk_formulae/input_file/file_name", option)

       dataset = "ERA40" ! default
       if(have_option("/ocean_forcing/bulk_formulae/file_type")) then
          call get_option("/ocean_forcing/bulk_formulae/file_type/filetype/name", dataset)
       end if
       
       if (dataset == "ERA40") then
          call fluxes_registerdatafile(trim(option))
          !                field from NetCDF file     dex |   Physical meaning
          call fluxes_addfieldofinterest("u10")   !   0   | 10 metre U wind component
          call fluxes_addfieldofinterest("v10")   !   1   | 10 metre V wind component
          call fluxes_addfieldofinterest("ssrd")  !   2   | Surface solar radiation
          call fluxes_addfieldofinterest("strd")  !   3   | Surface thermal radiation 
          call fluxes_addfieldofinterest("ro")    !   4   | Runoff
          call fluxes_addfieldofinterest("tp")    !   5   | Total precipitation
          call fluxes_addfieldofinterest("d2m")   !   6   | Dew point temp at 2m
          call fluxes_addfieldofinterest("t2m")   !   7   | Air temp at 2m 
          call fluxes_addfieldofinterest("msl")   !   8   | Mean sea level pressure
       else
          FLExit("unsupported bulk formula input file type. Choose ERA40")
       end if
    end if

    if(have_option("/ocean_biology/lagrangian_ensemble/hyperlight")) then
       call fluxes_addfieldofinterest("tcc")      !       | Total cloud cover 
    end if

    if(have_option("/ocean_forcing/external_data_boundary_conditions")) then
       call get_option("/ocean_forcing/external_data_boundary_conditions/input_file/file_name", option)

       ewrite(2,*)"Registering external data forcing file: " // trim(option)
       call NEMO_v2_RegisterDataFile(option)

       call NEMO_v2_AddFieldOfInterest("temperature")   !   0   | Sea temperature
       call NEMO_v2_AddFieldOfInterest("salinity")      !   1   | Salinity
       call NEMO_v2_AddFieldOfInterest("u")             !   2   | Azimuthal velocity
       call NEMO_v2_AddFieldOfInterest("v")             !   3   | Meridional velocity
       call NEMO_v2_AddFieldOfInterest("ssh")           !   4   | Sea surface height
    end if

    ewrite(1, *) "Exiting initialise_ocean_forcing_readers"
  end subroutine initialise_ocean_forcing_readers

end module initialise_ocean_forcing_module
