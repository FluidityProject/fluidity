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
#include "fdebug.h"

module hyperlight
  use spud
  use global_parameters, only:   OPTION_PATH_LEN
  use fields
  use state_module
  use fluxes

  implicit none
contains

  subroutine set_irradiance_from_hyperlight(state)
    type(state_type), intent(inout) :: state
    type(vector_field),pointer  :: coord
    type(scalar_field),pointer :: chlorophyll, irradiance_field
    character(len=OPTION_PATH_LEN) :: field_name
    character(len=1024) :: time_units
    integer :: node, date, f
    real :: x, y, z, lat_long(2), chl, bf_chl, cdom, wind, wind_u, wind_v, cloud, time
    real :: irradiance, lambda, timestep, scalar, scalars(9), euphotic_ratio

    ewrite(1,*) "Running Hyperlight"
    call hyperlight_reset()
    ewrite(2,*) "Hyperlight: setting default parameters"

    ! passing current and start time to Hyperlight
    call get_option("/timestepping/current_time/time_units/date", time_units)
    call get_option("/timestepping/current_time", time)
    call hyperlight_set_date_time(time, trim(time_units))

    ! passing lat long (single position only, for now...)
    if (have_option("/ocean_forcing/bulk_formulae/position/single_location")) then
       call hyperlight_set_single_position(1)
       call get_option("/ocean_forcing/bulk_formulae/position", lat_long)
       call hyperlight_set_coords(lat_long(1), lat_long(2))
    else
       ewrite(-1,*) "Hyperlight currently requires ocean forcing with single_location!"
       FLExit("Hyperlight: Incorrect latitude and longitude specified.")
    end if

    ! passing BF and CDOM ratio as specified
    call get_option("/ocean_biology/lagrangian_ensemble/hyperlight/BF_chl", bf_chl)
    call get_option("/ocean_biology/lagrangian_ensemble/hyperlight/CDOM", cdom)

    ! getting windspeed and cloudcover
    if (have_option("/ocean_biology/lagrangian_ensemble/hyperlight/CloudCover")) then
       call get_option("/ocean_biology/lagrangian_ensemble/hyperlight/CloudCover", cloud)
    else
       call fluxes_getscalar("tcc", lat_long(2), lat_long(1), cloud)
    end if
    if (have_option("/ocean_biology/lagrangian_ensemble/hyperlight/WindSpeed")) then
       call get_option("/ocean_biology/lagrangian_ensemble/hyperlight/WindSpeed", wind)
    else
       call fluxes_getscalar("u10", lat_long(2), lat_long(1), wind_u)
       call fluxes_getscalar("v10", lat_long(2), lat_long(1), wind_v)
       wind =  sqrt(wind_u**2.0 + wind_v**2.0);
    end if
    call hyperlight_set_params(bf_chl, cdom, wind, cloud)

    ! Performance parameter: set the percentage of surface irradiance after which Hyperlight stops computing
    if (have_option("/ocean_biology/lagrangian_ensemble/hyperlight/EuphoticRatio")) then
       call get_option("/ocean_biology/lagrangian_ensemble/hyperlight/EuphoticRatio", euphotic_ratio)
       call hyperlight_set_euphotic_ratio(euphotic_ratio)
    else
       call hyperlight_set_euphotic_ratio(0.0)
    end if

    ! set nodes on Hyperlight grid
    ewrite(2,*) "Hyperlight: creating grid"
    coord=>extract_vector_field(state, "Coordinate")
    chlorophyll=>extract_scalar_field(state, "Chlorophyll")
    do node=1,node_count(coord)
       chl = node_val(chlorophyll, node)
       call hyperlight_set_node(node_val(coord, node), chl)
    end do

    ! run it
    ewrite(2,*) "Hyperlight: run()"
    call hyperlight_run()

    ! copy Hyperlight result to irradiance fields
    ewrite(2,*) "Hyperlight: setting irradiance fields"
    frequency_field_loop: do f=0,35
       lambda = 350.0 + (f * 10.0)
       field_name="Irradiance_"//int2str(NINT(lambda))
       irradiance_field=>extract_scalar_field(state, field_name)
       node_loop: do node=1,node_count(coord)
         call hyperlight_query_node(node_val(coord, node), lambda, irradiance)
         call set(irradiance_field, node, irradiance)
       end do node_loop
    end do frequency_field_loop

  end subroutine set_irradiance_from_hyperlight

  subroutine hyperlight_init()
    call hyperlight_grid_init
  end subroutine hyperlight_init

  subroutine hyperlight_set_node(coord, chl)
    real, intent(in) :: chl
    real, dimension(:), intent(in) :: coord
    real :: z

    assert(size(coord)==3)
    z = coord(3)
    if (z .lt. 0.0) then
        z = -coord(3)
    end if
    call hyperlight_grid_set_node(coord(1), coord(2), z, chl)
  end subroutine hyperlight_set_node

  subroutine hyperlight_query_node(coord, lambda, irradiance)
    real, intent(in) :: lambda
    real, dimension(:), intent(in) :: coord
    real, intent(out) :: irradiance
    real :: z

    assert(size(coord)==3)
    z = coord(3)
    if (z .lt. 0.0) then
        z = -coord(3)
    end if
    call hyperlight_grid_query_node(coord(1), coord(2), z, lambda, irradiance)
  end subroutine hyperlight_query_node

end module hyperlight
