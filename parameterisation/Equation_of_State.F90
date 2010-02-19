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

module equation_of_state
  !!< This module contains functions used to evaluate the equation of state.
  use fldebug
  use fields
  use state_module
  use global_parameters, only: OPTION_PATH_LEN
  use spud
  implicit none
  
  private
  public :: calculate_perturbation_density, mcD_J_W_F2002, &
            compressible_eos, compressible_material_eos

contains

  subroutine calculate_perturbation_density(state, density, reference_density)
    !!< Calculates the perturbation density (i.e. the reference density is already subtracted)
    !!< of a state with equation_of_state fluids/linear or 
    !!< fluids/ocean_pade_approximation.
    type(state_type), intent(in):: state
    type(scalar_field), intent(inout) :: density
    real, intent(out), optional :: reference_density
    
    type(vector_field), pointer:: u
    type(scalar_field), pointer:: T, S, oldT, oldS, topdis
    type(scalar_field) DeltaT, DeltaS, remapT, remapS, fluidconcentration,&
         & sedimentdensity
    character(len=OPTION_PATH_LEN) option_path, dep_option_path, class_name
    logical, dimension(:), allocatable:: done
    logical include_depth_below
    real T0, S0, gamma, rho_0, salt, temp, dist, dens, theta
    integer, dimension(:), pointer:: density_nodes
    integer ele, i, node
    
    ewrite(1,*) 'In calculate_perturbation_density'
    
    u => extract_vector_field(state, "Velocity")
    
    call zero(density)
    
    option_path='/material_phase::'//trim(state%name)//'/equation_of_state/fluids'
    
    call get_option(trim(u%option_path)//'/prognostic/temporal_discretisation/relaxation', &
                    theta, default = 1.0)
    
    rho_0 = 0.0
    
    if (have_option(trim(option_path)//'/linear')) then
    
       option_path=trim(option_path)//'/linear'
       
       if (have_option(trim(option_path)//'/temperature_dependency')) then
          dep_option_path=trim(option_path)//'/temperature_dependency'
          call get_option(trim(dep_option_path)//'/reference_temperature', T0)
          call get_option(trim(dep_option_path)//'/thermal_expansion_coefficient', gamma)
          T => extract_scalar_field(state, "Temperature")
          oldT => extract_scalar_field(state, "OldTemperature")
          call allocate(deltaT, density%mesh, "DeltaT")
          call allocate(remapT, density%mesh, "RemapT")
          
          ! deltaT=theta*T+(1-theta)*oldT-T0
          call remap_field(T, remapT)
          call set(deltaT, remapT)
          call scale(deltaT, theta)
          
          call remap_field(oldT, remapT)
          call addto(deltaT, remapT, 1.0-theta)
          call addto(deltaT, -T0)
          ! density=density-gamma*deltaT
          call addto(density, deltaT, scale=-gamma)
          call deallocate(deltaT)
          call deallocate(remapT)
       end if
       
       if (have_option(trim(option_path)//'/salinity_dependency')) then
          dep_option_path=trim(option_path)//'/salinity_dependency'
          call get_option(trim(dep_option_path)//'/reference_salinity', S0)
          call get_option(trim(dep_option_path)//'/saline_contraction_coefficient', gamma)
          S => extract_scalar_field(state, "Salinity")
          oldS => extract_scalar_field(state, "OldSalinity")
          call allocate(deltaS, density%mesh, "DeltaS")
          call allocate(remapS, density%mesh, "RemapS")
          
          ! deltaS=theta*S+(1-theta)*oldS-S0
          call remap_field(S, remapS)
          call set(deltaS, remapS)
          call scale(deltaS, theta)
          
          call remap_field(oldS, remapS)
          call addto(deltaS, remapS, 1.0-theta)
          call addto(deltaS, -S0)
          ! density=density+gamma*deltaS
          call addto(density, deltaS, scale=gamma)
          call deallocate(deltaS)
          call deallocate(remapS)
       end if
       
       call get_option(trim(option_path)//'/reference_density', rho_0)
       call scale(density, rho_0)
       
    elseif (have_option(trim(option_path)//'/ocean_pade_approximation')) then
      
       option_path=trim(option_path)//'/ocean_pade_approximation'
       
       include_depth_below=have_option(trim(option_path)//'/include_depth_below_surface')
       
       T => extract_scalar_field(state, "Temperature")
       oldT => extract_scalar_field(state, "OldTemperature")
       S => extract_scalar_field(state, "Salinity")
       oldS => extract_scalar_field(state, "OldSalinity")
       if (include_depth_below) then
          topdis => extract_scalar_field(state, "DistanceToTop")
       endif
      
       allocate( done(1:node_count(density)) )
       done=.false.
       
       do ele=1, element_count(density)

          density_nodes => ele_nodes(density, ele)

          do i=1,size(density_nodes)
             node=density_nodes(i)
             ! In the continuous case ensure we only do each calculation once.
             if (done(node)) cycle
             done(node)=.true.
            
             salt=theta*node_val(S, node)+(1-theta)*node_val(oldS, node)
             temp=node_val(T, node)+(1-theta)*node_val(oldT, node)
             if (include_depth_below) then
                dist=node_val(topdis, node)
             else
                dist=0.0
             end if            
               
             call mcD_J_W_F2002(dens,temp,salt,dist)
             call addto(density, node, dens)
          end do
           
       end do
         
       ! reference density is assumed 1 for the pade approximation
       rho_0=1.0

    end if
    
    if (have_option('/material_phase::'//trim(state%name)//'/sediment'))&
         & then 

       call allocate(deltaS, density%mesh, "DeltaS")
       call allocate(remapS, density%mesh, "RemapS")
       call allocate(sedimentdensity, density%mesh, "SedimentDensity")
       ! fluidconcentration is 1-sedimentconcentration.
       call allocate(fluidconcentration, density%mesh, &
            "FluidConcentration")
       call zero(sedimentdensity)
       call set(fluidconcentration, 1.0)

       do i=1,option_count('/material_phase::'//trim(state%name)&
            //'/sediment/sediment_class')
          option_path='/material_phase::'//trim(state%name)//&
               '/sediment/sediment_class['//int2str(i-1)//"]"
          
          call get_option(trim(option_path)//"/name", class_name)

          call get_option(trim(option_path)//'/density', gamma)
          ! Note that 1000.0 is a hack. We actually need to properly
          ! account for the reference density of seawater.
          gamma=rho_0*(gamma/1000.0) - rho_0

          S => extract_scalar_field(state, &
               "SedimentConcentration"//trim(class_name))
          oldS => extract_scalar_field(state, &
               "OldSedimentConcentration"//trim(class_name))
          
          ! deltaS=theta*S+(1-theta)*oldS-S0
          call remap_field(S, remapS)
          call set(deltaS, remapS)
          call scale(deltaS, theta)
          
          call remap_field(oldS, remapS)
          call addto(deltaS, remapS, 1.0-theta)
          ! density=density+gamma*deltaS
          call addto(sedimentdensity, deltaS, scale=gamma)
          call addto(fluidconcentration, deltaS, scale=-1.0)

       end do
       
       call scale(density,fluidconcentration)
       call addto(density,sedimentdensity)

       call deallocate(deltaS)
       call deallocate(remapS)
       call deallocate(fluidconcentration)
       call deallocate(sedimentdensity)
       
    end if

    if(present(reference_density)) then
      reference_density = rho_0
    end if
    
  end subroutine calculate_perturbation_density

  subroutine mcD_J_W_F2002(density,T,Salinity,distance_to_top)
    !!<  function to evaluate density from the 2002 McDougall, Jackett,
    !!< Wright and Feistel equation of state using Pade approximation.  
    real, intent(out) :: density
    real, intent(in) :: T,Salinity,distance_to_top
   
    real :: p,p1,p2,S

    ! Salinity can be negitive because it's numerically diffused,
    ! some regions may be initialised with zero salinity, and
    ! therefore undershoot may occur.
    S = max(Salinity, 0.0)

    ! calculate pressure in decibars from hydrostatic pressure
    ! using reference density 1000 kg m^-2

    p = 9.81*1000.0*distance_to_top*1.0e-4

    !     evaluate top and bottom of Pade approximant
      
    p1 = 9.99843699e2 &
         + 7.35212840*T - 5.45928211e-2*(T**2) + 3.98476704e-4*(T**3) &
         + 2.96938239*S - 7.23268813e-3*S*T + 2.12382341e-3*(S**2) &
         + 1.04004591e-2*p + 1.03970529e-7*p*(T**2) &
         + 5.18761880e-6*p*S - 3.24041825e-8*(p**2) &
         - 1.23869360e-11*(p**2)*(t**2)

    p2 = 1.0 &
         + 7.28606739e-3*T - 4.60835542e-5*(T**2) + 3.68390573e-7*(T**3) &
         + 1.80809186e-10*(T**4) &
         + 2.14691708e-3*S - 9.27062484e-6*S*T - 1.78343643e-10*S*(T**3) &
         + 4.76534122e-6*(S**1.5) + 1.63410736e-9*(S**1.5)*(T**2) &
         + 5.30848875e-6*p -3.03175128e-16*(p**2)*(t**3) &
         - 1.27934137e-17*(p**3)*T
    
    ! calculate the resulting density
    
    density = p1/p2
    
    ! the perturbation density
    
    density = (density-1000.0)/1000.0
    
  end subroutine mcD_J_W_F2002
  
  subroutine compressible_eos(state,density,&
                                    pressure,drhodp)

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout), optional :: density, &
                                             pressure, drhodp

    !locals
    integer :: stat, gstat, cstat
    type(scalar_field), pointer :: pressure_local, energy_local, density_local
    character(len=4000) :: thismaterial_phase, eos_path
    real :: reference_density, ratio_specific_heats
    real :: bulk_sound_speed_squared, atmospheric_pressure
    type(scalar_field) :: drhodp_local, energy_remap, pressure_remap, density_remap
    logical :: incompressible

    ewrite(1,*) 'Entering compressible_eos'

    if (present(density)) then
      call allocate(drhodp_local, density%mesh, 'Localdrhop')
    else if (present(pressure)) then
      call allocate(drhodp_local, pressure%mesh, 'Localdrhop')
    else if (present(drhodp)) then
      call allocate(drhodp_local, drhodp%mesh, 'Localdrhop')
    else
      FLAbort("No point in being in here if you don't want anything out.")
    end if

    thismaterial_phase = '/material_phase::'//trim(state%name)
    eos_path = trim(thismaterial_phase)//'/equation_of_state'

    if(have_option(trim(eos_path)//'/compressible')) then

      if(have_option(trim(eos_path)//'/compressible/miegrunneisen')) then
        call get_option(trim(eos_path)//'/compressible/miegrunneisen/reference_density', &
                        reference_density, default=0.0)
        
        call get_option(trim(eos_path)//'/compressible/miegrunneisen/ratio_specific_heats', &
                        ratio_specific_heats, stat=gstat)
        if(gstat/=0) then
          ratio_specific_heats=1.0
        end if
        
        call get_option(trim(eos_path)//'/compressible/miegrunneisen/bulk_sound_speed_squared', &
                        bulk_sound_speed_squared, stat=cstat)
        if(cstat/=0) then
          bulk_sound_speed_squared=0.0
        end if
        
        incompressible = ((gstat/=0).and.(cstat/=0))
        if(incompressible) then
          ewrite(0,*) "Selected compressible eos but not specified a bulk_sound_speed_squared or a ratio_specific_heats."
        end if
        
        call zero(drhodp_local)
        
        if(.not.incompressible) then
          energy_local=>extract_scalar_field(state,'InternalEnergy',stat=stat)
          ! drhodp = 1.0/( bulk_sound_speed_squared + (ratio_specific_heats - 1.0)*energy )
          if((stat==0).and.(gstat==0)) then   ! we have an internal energy field and we want to use it
            call allocate(energy_remap, drhodp_local%mesh, 'RemappedInternalEnergy')
            call remap_field(energy_local, energy_remap)
            
            call addto(drhodp_local, energy_remap, (ratio_specific_heats-1.0))
            
            call deallocate(energy_remap)
          end if
          call addto(drhodp_local, bulk_sound_speed_squared)
          call invert(drhodp_local)
        end if

        if(present(density)) then
          ! calculate the density
          ! density may equal density in state depending on how this
          ! subroutine is called
          if(incompressible) then
            ! density = reference_density
            call set(density, reference_density)
          else
            pressure_local=>extract_scalar_field(state,'Pressure',stat=stat)
            if (stat==0) then
              assert(density%mesh==drhodp_local%mesh)
            
              ! density = drhodp_local*(pressure_local + atmospheric_pressure
              !                       + bulk_sound_speed_squared*reference_density)
              call get_option(trim(pressure_local%option_path)//'/prognostic/atmospheric_pressure', &
                              atmospheric_pressure, default=0.0)
              
              call allocate(pressure_remap, drhodp_local%mesh, "RemappedPressure")
              call remap_field(pressure_local, pressure_remap)
              
              call set(density, reference_density*bulk_sound_speed_squared + atmospheric_pressure)
              call addto(density, pressure_remap)
              call scale(density, drhodp_local)
              
              call deallocate(pressure_remap)
            else
              FLExit('No Pressure in material_phase::'//trim(state%name))
            end if
          end if
        end if

        if(present(pressure)) then
          if(incompressible) then
            ! pressure is unrelated to density in this case
            call zero(pressure)
          else
            ! calculate the pressure using the eos and the calculated (probably prognostic)
            ! density
            density_local=>extract_scalar_field(state,'Density',stat=stat)
            if (stat==0) then
              assert(pressure%mesh==drhodp_local%mesh)
              
              ! pressure = density_local/drhodp_local &
              !          - bulk_sound_speed_squared*reference_density
              
              call allocate(density_remap, drhodp_local%mesh, "RemappedDensity")
              call remap_field(density_local, density_remap)
              
              call set(pressure, drhodp_local)
              call invert(pressure)
              call scale(pressure, density_remap)
              call addto(pressure, -bulk_sound_speed_squared*reference_density)
              
              call deallocate(density_remap)
            else
              FLExit('No Density in material_phase::'//trim(state%name))
            end if
          end if
        end if

!       elseif(have_option(trim(eos_path)//'/compressible/foam')) then
      ! perhaps place a new foam eos here?
      
!       else
      ! place other compressible eos here

      end if

    end if

    if(present(density)) then
      ewrite_minmax(density%val)
    end if

    if(present(pressure)) then
      ewrite_minmax(pressure%val)
    end if

    if(present(drhodp)) then
      if((cstat/=0).and.(gstat/=0)) then
        ! pressure is unrelated to density in this case
        call zero(drhodp)
      else
        assert(drhodp%mesh==drhodp_local%mesh)
        call set(drhodp, drhodp_local)
      end if
      ewrite_minmax(drhodp%val)
    end if

    call deallocate(drhodp_local)

  end subroutine compressible_eos
  
  subroutine compressible_material_eos(state,materialdensity,&
                                    materialpressure,materialdrhodp)

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout), optional :: materialdensity, &
                                             materialpressure, materialdrhodp

    !locals
    integer :: stat, gstat, cstat
    type(scalar_field), pointer :: pressure, materialenergy, materialdensity_local
    character(len=4000) :: thismaterial_phase, eos_path
    real :: reference_density, ratio_specific_heats
    real :: bulk_sound_speed_squared, atmospheric_pressure
    type(scalar_field) :: drhodp

    ewrite(1,*) 'Entering compressible_material_eos'

    if (present(materialdensity)) then
      call allocate(drhodp, materialdensity%mesh, 'Gradient of density wrt pressure')
    else if (present(materialpressure)) then
      call allocate(drhodp, materialpressure%mesh, 'Gradient of density wrt pressure')
    else if (present(materialdrhodp)) then
      call allocate(drhodp, materialdrhodp%mesh, 'Gradient of density wrt pressure')
    else
      FLAbort("No point in being in here if you don't want anything out.")
    end if

    thismaterial_phase = '/material_phase::'//trim(state%name)
    eos_path = trim(thismaterial_phase)//'/equation_of_state'

    if(have_option(trim(eos_path)//'/compressible')) then

      if(have_option(trim(eos_path)//'/compressible/miegrunneisen')) then
        call get_option(trim(eos_path)//'/compressible/miegrunneisen/reference_density', &
                        reference_density, default=0.0)
        call get_option(trim(eos_path)//'/compressible/miegrunneisen/ratio_specific_heats', &
                        ratio_specific_heats, stat=gstat)
        if(gstat/=0) then
          ratio_specific_heats=1.0
        end if
        call get_option(trim(eos_path)//'/compressible/miegrunneisen/bulk_sound_speed_squared', &
                        bulk_sound_speed_squared, stat = cstat)
        if(cstat/=0) then
          bulk_sound_speed_squared=0.0
        end if
        if((gstat/=0).and.(cstat/=0)) then
          FLExit("Must set either a bulk_sound_speed_squared or a ratio_specific_heats.")
        end if
        materialenergy=>extract_scalar_field(state,'MaterialInternalEnergy',stat=stat)
        if(stat==0) then   ! we have an internal energy field
          drhodp%val = 1.0/( bulk_sound_speed_squared + (ratio_specific_heats - 1.0)*materialenergy%val )
        else               ! we don't have an internal energy field
          call set(drhodp, 1.0/bulk_sound_speed_squared)
        end if

        if(present(materialdensity)) then
          ! calculate the materialdensity
          ! materialdensity can equal materialdensity in state depending on how this
          ! subroutine is called
          pressure=>extract_scalar_field(state,'Pressure',stat=stat)
          if (stat==0) then
            call get_option(trim(pressure%option_path)//'/prognostic/atmospheric_pressure', &
                            atmospheric_pressure, default=0.0)
            materialdensity%val = drhodp%val*(pressure%val + atmospheric_pressure &
                                   + bulk_sound_speed_squared*reference_density)
          else
            FLExit('No Pressure in material_phase::'//trim(state%name))
          end if
        end if

        if(present(materialpressure)) then
          ! calculate the materialpressure using the eos and the calculated (probably prognostic)
          ! materialdensity
          ! materialpressure /= bulk pressure
          materialdensity_local=>extract_scalar_field(state,'MaterialDensity',stat=stat)
          if (stat==0) then
            materialpressure%val = materialdensity_local%val/drhodp%val &
                                  - bulk_sound_speed_squared*reference_density
          else
            FLExit('No MaterialDensity in material_phase::'//trim(state%name))
          end if
        end if

!       else
!       ! place other compressible material eos here

      end if

!     else
!     ! an incompressible option?

    end if      

    if(present(materialdensity)) then
      ewrite_minmax(materialdensity%val)
    end if

    if(present(materialpressure)) then
      ewrite_minmax(materialpressure%val)
    end if

    if(present(materialdrhodp)) then
      materialdrhodp%val=drhodp%val
      ewrite_minmax(materialdrhodp%val)
    end if

    call deallocate(drhodp)

  end subroutine compressible_material_eos
  
end module equation_of_state
