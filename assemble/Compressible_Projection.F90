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

module compressible_projection
  use fldebug
  use state_module
  use sparse_tools
  use spud
  use fields
  use sparse_matrices_fields
  use multimaterial_module, only: extract_prognostic_pressure, &
                                  calculate_material_eos
  use global_parameters, only: new_options, OPTION_PATH_LEN
  use legacy_cv_shape_functions
  use fefields, only: compute_lumped_mass
  implicit none 

  ! Buffer for output messages.
  character(len=255), private :: message

  private
  public assemble_compressible_projection

contains

  subroutine assemble_compressible_projection(state, cmc, rhs, dt, cmcget)

    ! inputs:
    ! bucket full of fields
    type(state_type), dimension(:), intent(inout) :: state

    type(csr_matrix), intent(inout) :: cmc
    type(scalar_field), intent(inout) :: rhs

    real, intent(in) :: dt
    logical, intent(in) :: cmcget

    ! local:
    integer :: i, stat, norm_stat
    character(len=OPTION_PATH_LEN) :: pressure_option_path
    character(len=FIELD_NAME_LEN) :: normalisation_field

    type(scalar_field) :: materialpressure, materialdrhodp, normdensity, &
                          normolddensity, normmatdrhodpp, normdrhodp
    type(scalar_field), pointer :: normalisation, &
                                   volumefraction, oldvolumefraction, materialdensity, oldmaterialdensity
    type(scalar_field), pointer :: dummy_ones

    type(scalar_field), pointer :: pressure
    type(vector_field), pointer :: positions
    type(scalar_field) :: lumped_mass, tempfield

    ! Cause the pressure warnings to only happen once.
    logical, save :: pressure_warned=.false.

    real :: atmospheric_pressure

    ewrite(1,*) 'Entering assemble_compressible_projection'

    pressure_option_path=""
    pressure=>extract_prognostic_pressure(state, stat=stat)
    if(stat==0) then
       pressure_option_path=trim(pressure%option_path)
    else if((stat==1).and.(new_options).and.(.not.pressure_warned)) then
       ewrite(0,*) "Warning: No prognostic pressure found in state."
       ewrite(0,*) "Strange that you've got here without a prognostic press&
            &ure."
       pressure_warned=.true.
    else if((stat==2).and.(new_options).and.(.not.pressure_warned)) then
       ewrite(0,*) "Warning: Multiple prognostic pressures found."
       ewrite(0,*) "Not sure if new options are compatible with this yet."
       pressure_warned=.true.
    end if
   
    if(have_option(trim(pressure_option_path)//"/prognostic/scheme/use_compressible_projection_method")) THEN

      ! only do all this if we need to make cmc (otherwise we'd be adding repeatedly)
      if(cmcget) then

        positions=>extract_vector_field(state(1), "Coordinate")
        call allocate(lumped_mass, pressure%mesh, "LumpedMassField")
        call allocate(tempfield, pressure%mesh, "TemporaryAssemblyField")
        call compute_lumped_mass(positions, lumped_mass)

        allocate(dummy_ones)
        call allocate(dummy_ones, pressure%mesh, "DummyOnesField")
        call set(dummy_ones, 1.0)

        call zero(rhs)

        call get_option(trim(pressure_option_path)//'/prognostic/atmospheric_pressure', &
                        atmospheric_pressure, default=0.0)

        call get_option(trim(pressure_option_path)//"/prognostic/scheme/use_compressible_projection_method/normalisation/name", &
                        normalisation_field, stat=norm_stat)

        call allocate(materialpressure, pressure%mesh, 'MaterialEOSPressure')
        call allocate(materialdrhodp, pressure%mesh, 'DerivativeMaterialdensityWRTBulkPressure')

        call allocate(normdensity, pressure%mesh, 'NormalisedMaterialDensity')
        call allocate(normolddensity, pressure%mesh, 'NormalisedOldMaterialDensity')
        call allocate(normmatdrhodpp, pressure%mesh, 'NormalisedMaterialPressure')
        call allocate(normdrhodp, pressure%mesh, 'NormalisedDrhodp')

        normdensity%val = 0.0
        normolddensity%val = 0.0
        normmatdrhodpp%val = 0.0
        normdrhodp%val=0.0

        do i = 1,size(state)

          materialpressure%val=0.0
          materialdrhodp%val=0.0

          call calculate_material_eos(state(i), materialpressure=materialpressure, materialdrhodp=materialdrhodp)

          volumefraction=>extract_scalar_field(state(i),'MaterialVolumeFraction', stat=stat)
          if(stat==0) then
            oldvolumefraction=>extract_scalar_field(state(i),'OldMaterialVolumeFraction')
            materialdensity=>extract_scalar_field(state(i),'MaterialDensity')
            oldmaterialdensity=>extract_scalar_field(state(i),'OldMaterialDensity')

            if(norm_stat==0) then
              normalisation=>extract_scalar_field(state(i), trim(normalisation_field))
            else
              normalisation=>dummy_ones
            end if

            normdensity%val = normdensity%val &
                              + materialdensity%val*volumefraction%val/ &
                                normalisation%val
            normolddensity%val = normolddensity%val &
                                + oldmaterialdensity%val*oldvolumefraction%val/ &
                                  normalisation%val
            normmatdrhodpp%val = normmatdrhodpp%val &
                                  + materialpressure%val*materialdrhodp%val*volumefraction%val/ &
                                    normalisation%val
            normdrhodp%val = normdrhodp%val &
                              + materialdrhodp%val*volumefraction%val/ &
                                normalisation%val
          endif

        end do

        call zero(tempfield)
        tempfield%val = (1./(dt*dt))*lumped_mass%val*normdrhodp%val

        call addto_diag(cmc, tempfield)

        rhs%val = (1./dt)*lumped_mass%val* &
                          ( &
                            normolddensity%val &
                          - normdensity%val &
                          ) &
               +(1./dt)*lumped_mass%val* &
                          ( &
                            normmatdrhodpp%val &
                          - normdrhodp%val*(pressure%val+atmospheric_pressure) &
                          )

        call deallocate(normdensity)
        call deallocate(normolddensity)
        call deallocate(normmatdrhodpp)
        call deallocate(normdrhodp)

        call deallocate(materialpressure)
        call deallocate(materialdrhodp)

        call deallocate(lumped_mass)
        call deallocate(tempfield)
        call deallocate(dummy_ones)
        deallocate(dummy_ones)

      end if

    end if

  end subroutine assemble_compressible_projection

end module compressible_projection

