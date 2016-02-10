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

module OceanSurfaceForcing
  use global_parameters, only : OPTION_PATH_LEN, pi
  use elements
  use spud
  use transform_elements
  use fetools
  use fields
  use state_module
  use boundary_conditions
  use coordinates
  
  implicit none

  private
  public:: wind_forcing

contains

  subroutine wind_forcing(state, rhs)
    !!< Implements wind forcing from new options.
    type(state_type), intent(in):: state
    type(vector_field), intent(inout):: rhs
    
    type(vector_field), pointer:: velocity, positions, wind_surface_field
    type(scalar_field), pointer:: wind_drag_coefficient
    type(element_type):: faceshape
    character(len=OPTION_PATH_LEN):: bctype
    real, dimension(:,:), allocatable:: llpos_at_quads, ll_wind_at_quads
    real, dimension(:,:), allocatable:: wind_at_quads
    real, dimension(:), allocatable:: detwei, C_D, unorm
    real rho_air
    logical apply_wind_formula, on_sphere
    integer, dimension(:), pointer:: surface_element_list
    integer sngi, wdim, nobcs
    integer i, j, k, sele
    logical:: parallel_dg
    
    ewrite(1,*) 'Inside wind_forcing'
    ewrite_minmax(rhs)
   
    velocity => extract_vector_field(state, "Velocity")
    positions => extract_vector_field(state, "Coordinate")
   
    parallel_dg=continuity(velocity)<0 .and. IsParallel()
    on_sphere=have_option('/geometry/spherical_earth')
    faceshape=face_shape(velocity, 1)
    sngi=face_ngi(velocity, 1)
    if (on_sphere) then
      wdim=3 ! dimension of the wind
      assert(velocity%dim==3)
    else
      wdim=velocity%dim-1 ! dimension of the wind
    end if
    
    allocate(detwei(1:sngi), C_D(1:sngi), unorm(1:sngi), &
      wind_at_quads(1:wdim,1:sngi))
    if (on_sphere) then
      allocate( llpos_at_quads(1:2,1:sngi), ll_wind_at_quads(1:2,1:sngi) )
    end if
    
    nobcs = get_boundary_condition_count(velocity)
    do i=1, nobcs
      call get_boundary_condition(velocity, i, type=bctype, &
          surface_element_list=surface_element_list)
      if (bctype=='wind_forcing') then  
        wind_surface_field => extract_surface_field(velocity, i, "WindSurfaceField")
        apply_wind_formula=has_scalar_surface_field(velocity, i, "WindDragCoefficient")
        if (apply_wind_formula) then
           wind_drag_coefficient => extract_scalar_surface_field(velocity, &
              i, "WindDragCoefficient")
           call get_option(trim(velocity%option_path)// &
              '/prognostic/boundary_conditions['//int2str(i-1)//']&
              &/type[0]/wind_velocity/density_air', rho_air)
        end if
              
        do j=1, size(surface_element_list) 
          sele=surface_element_list(j) ! face/surface element nr.
          
          if (parallel_dg) then
            if (.not. element_owned(velocity, face_ele(velocity,sele))) cycle
          end if
          
          ! compute integration weights detwei
          call transform_facet_to_physical(positions, sele, detwei)
                    
          ! compute wind_at_quads: wind velocity at quadr. points
          ! OR (if .not. apply_wind_formula): wind stress at quadr. points
          if (on_sphere) then
            call LongitudeLatitude(face_val_at_quad(positions,sele), &
               llpos_at_quads(1,:), llpos_at_quads(2,:))
            ll_wind_at_quads=ele_val_at_quad(wind_surface_field, j)
            call ll2r3_rotate( &
                llpos_at_quads(1,:), llpos_at_quads(2,:), &
                ll_wind_at_quads(1,:), ll_wind_at_quads(2,:), &
                wind_at_quads(1,:), wind_at_quads(2,:), wind_at_quads(3,:))
          else
          
            wind_at_quads=ele_val_at_quad(wind_surface_field, j)
          end if


          if (apply_wind_formula) then
            
            ! wind_at_quads is actual wind velocity
            
            ! make wind_at_quads relative velocity (specified wind-sea surface velocity)
            ! at quadrature points
            do k=1, wdim
              wind_at_quads(k,:)=wind_at_quads(k,:)- &
                  face_val_at_quad(velocity, sele, dim=k)
            end do
            ! drag coefficient:
            C_D=ele_val_at_quad(wind_drag_coefficient, j)
          
            ! compute its norm, sum over dim=1 is sum over components
            unorm=sqrt(sum(wind_at_quads**2, dim=1))

            ! multiply at each gauss point:
            detwei=detwei*C_D*unorm*rho_air
          end if
          
          do k=1, wdim
            ! add surface forcing in rhs of momentum equation:
            call addto(rhs, k, face_global_nodes(velocity, sele), &
              shape_rhs( faceshape, wind_at_quads(k,:)*detwei ))
          end do
        end do
      end if
    end do
      
    deallocate(detwei, C_D, unorm, wind_at_quads)
    if (on_sphere) then
      deallocate(llpos_at_quads, ll_wind_at_quads)
    end if
    
    ewrite_minmax(rhs)
   
  end subroutine wind_forcing

end module OceanSurfaceForcing
