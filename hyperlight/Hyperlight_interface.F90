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
  use state_module
  use fields
  use global_parameters, only:   OPTION_PATH_LEN
  implicit none
contains

  subroutine set_irradiance_from_hyperlight(state)
    type(state_type), intent(inout) :: state
    type(vector_field) :: coord
    type(scalar_field) :: chlorophyll
    type(scalar_field) :: irradiance_field
    character(len=OPTION_PATH_LEN) :: field_name
    integer :: node, date, f
    real :: x, y, z, chl, time, lat, long, bf_chl, cdom, wind, cloud
    real :: irradiance, lambda

    ewrite(1,*) "Running Hyperlight"
    call hyperlight_reset()
    ewrite(2,*) "Hyperlight: setting default parameters"

    !!! some default configuration, for now ...
    date = 19970114
    time = 17.0
    call hyperlight_set_date_time(date, time)

    lat = 31.6960
    long = -64.1650
    call hyperlight_set_coords(lat, long)

    bf_chl = 0.0185
    cdom = 0.0
    wind =5.0
    cloud = 0.3
    call hyperlight_set_params(bf_chl, cdom, wind, cloud)
    !!!


    ! set nodes on Hyperlight grid
    ewrite(2,*) "Hyperlight: creating grid"
    coord=extract_vector_field(state, "Coordinate")
    chlorophyll=extract_scalar_field(state, "Chlorophyll")
    do node=1,node_count(coord)
       chl = node_val(chlorophyll, node)
       call hyperlight_set_node(node_val(coord, node), chl)
    end do

    ! run it
    ewrite(2,*) "Hyperlight: run()"
    call hyperlight_run()

    ! copy Hyperlight result to irradiance fields
    ewrite(2,*) "Hyperlight: setting irradiance fields"
    do node=1,node_count(coord)
       frequency_field_loop: do f=0,35
         lambda = 350.0 + (f * 10.0)
         call hyperlight_query_node(node_val(coord, node), lambda, irradiance)

         field_name="Irradiance_"//int2str(NINT(lambda))
         irradiance_field=extract_scalar_field(state, field_name)
         call set(irradiance_field, node, irradiance)
       end do frequency_field_loop
    end do

  end subroutine set_irradiance_from_hyperlight

  subroutine hyperlight_init()
    call hyperlight_grid_init_c
  end subroutine hyperlight_init

  subroutine hyperlight_reset()
    call hyperlight_grid_reset_c
  end subroutine hyperlight_reset
  
  subroutine hyperlight_set_date_time(date, time_gmt)
    integer, intent(in)::date
    real, intent(in)::time_gmt
    call hyperlight_grid_set_date_time_c(date, time_gmt)
  end subroutine hyperlight_set_date_time

  subroutine hyperlight_set_coords(latitude, longitude)
    real, intent(in)::latitude, longitude
    call hyperlight_grid_set_coords_c(latitude, longitude)
  end subroutine hyperlight_set_coords

  subroutine hyperlight_set_params(bf_chl, cdom_ratio, windspeed, cloudcover)
    real, intent(in)::bf_chl, cdom_ratio, windspeed, cloudcover
    call hyperlight_grid_set_params_c(bf_chl, cdom_ratio, windspeed, cloudcover)
  end subroutine hyperlight_set_params

  subroutine hyperlight_set_node(coord, chl)
    real, intent(in) :: chl
    real, dimension(:), intent(in) :: coord
    real :: z

    assert(size(coord)==3)
    z = coord(3)
    if (z .lt. 0.0) then
        z = -coord(3)
    end if
    call hyperlight_grid_set_node_c(coord(1), coord(2), z, chl)
  end subroutine hyperlight_set_node

  subroutine hyperlight_run()
    call hyperlight_grid_run_c()
  end subroutine hyperlight_run

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
    call hyperlight_grid_query_node_c(coord(1), coord(2), z, lambda, irradiance)
  end subroutine hyperlight_query_node

end module hyperlight
