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

module foam_drainage
  ! This module contains the options used for Foam drainage
  use fldebug
  use spud
  use fields
  use state_module
  use field_options
  use cv_upwind_values
  use state_fields_module
  use diagnostic_fields_matrices

  implicit none

  private
  public :: calculate_drainage_source_absor

contains

  subroutine calculate_drainage_source_absor(state)

    type(state_type), intent(inout) :: state

    integer :: i, stat
    type(vector_field), pointer :: K1, foamvel, liquidvelocity, liqcontentvel, source, absor
    type(scalar_field), pointer :: p, liquidcontent, K2

    type(scalar_field) :: p_remap, rho_remap, K2_remap
    type(vector_field) :: foamvel_remap, K1_remap
    type(vector_field), pointer :: x

    real :: atmospheric_pressure
    real, allocatable, dimension(:) :: absor_val

          x => extract_vector_field(state, "Coordinate")

          ! K1 is a vector based on the properties of the liquid within the Plateau borders: (0.0, -density*gravity/3*drag_coefficient*viscosity, 0.0)
          K1 => extract_vector_field(state,'VelocityDrainageK1')

          ! K2 is a scalar based on the properties of the liquid within the Plateau borders: sqrt(sqrt(3)-pi/2)*surface_tension/(6*drag_coefficient*viscosity)
          K2 => extract_scalar_field(state,'VelocityDrainageK2')

          p => extract_scalar_field(state,'Pressure')

          foamvel => extract_vector_field(state,'FoamVelocity')


          if (have_option("/material_phase[0]/vector_field::FoamLiquidContentVelocity/diagnostic")) then
            liqcontentvel => extract_vector_field(state,'FoamLiquidContentVelocity')
            liquidcontent => extract_scalar_field(state,'Density')
            liquidvelocity => extract_vector_field(state,'Velocity')
            call zero(liqcontentvel)
            call allocate(rho_remap, liquidvelocity%mesh, 'RemappedLiquidContent')
            call remap_field(liquidcontent, rho_remap) 
            do i=1, node_count(liqcontentvel)
                 call set(liqcontentvel, i, ( node_val(rho_remap, i)*node_val(liquidvelocity, i)  ) )
            enddo
            call deallocate(rho_remap)
          else
            ewrite (0,*)"WARNING: You should have a FoamLiquidContentVelocity field if you want to obtain the volumetric liquid flow at the boundaries"
          endif

          liquidvelocity => extract_vector_field(state,'Velocity')

          source => extract_vector_field(state,'VelocitySource')
          if (have_option(trim(source%option_path)//'/diagnostic')) then
            call get_option(trim(p%option_path)//'/prognostic/atmospheric_pressure', &
                            atmospheric_pressure, default=0.0)            
            call zero(source)

            call allocate(p_remap, source%mesh, 'RemappedPressure')
            call allocate(K1_remap, source%dim, source%mesh, 'RemappedK1')
            call allocate(K2_remap, source%mesh, 'RemappedK2')
            call allocate(foamvel_remap, source%dim, source%mesh, 'RemappedFoamVelocity')
            call remap_field(p, p_remap)            
            call remap_field(K1, K1_remap)
            call remap_field(K2, K2_remap)
            call remap_field(foamvel, foamvel_remap)

          do i=1, node_count(source)

            call set(source, i, ((-node_val(K1_remap, i)/node_val(K2_remap, i))*(node_val(p_remap, i)+atmospheric_pressure )**1.5 + (node_val(foamvel_remap,i)/node_val(K2_remap, i))*(node_val(p_remap, i)+atmospheric_pressure )**0.5))

          enddo

            call deallocate(p_remap)
            call deallocate(K1_remap)
            call deallocate(K2_remap)
            call deallocate(foamvel_remap)

          endif

          absor => extract_vector_field(state,'VelocityAbsorption')
          if (have_option(trim(absor%option_path)//'/diagnostic')) then
             call get_option(trim(p%option_path)//'/prognostic/atmospheric_pressure', &
                               atmospheric_pressure, default=3.6E-8)
             call zero(absor)

            call allocate(p_remap, absor%mesh, 'RemappedPressure')
            call allocate(K2_remap, source%mesh, 'RemappedK2')
            call remap_field(p, p_remap)
            call remap_field(K2, K2_remap)            

            allocate (absor_val(absor%dim))  

             do i=1, node_count(absor)
                absor_val = (1/node_val(K2_remap, i))*(node_val(p_remap, i)+atmospheric_pressure )**0.5
                call set(absor, i, absor_val )
             enddo

            deallocate(absor_val)

            call deallocate(p_remap)
            call deallocate(K2_remap)

          endif

  end subroutine calculate_drainage_source_absor

end module foam_drainage




