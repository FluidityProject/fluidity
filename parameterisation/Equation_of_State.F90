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

module equation_of_state
  !!< This module contains functions used to evaluate the equation of state.
  use fldebug
  use fields
  use state_module
  use global_parameters, only: OPTION_PATH_LEN
  use sediment, only: get_n_sediment_fields, get_sediment_item
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
    character(len=OPTION_PATH_LEN) option_path, dep_option_path, sediment_field_name, class_name, sfield_name
    logical, dimension(:), allocatable:: done
    logical include_depth_below
    real T0, S0, gamma, rho_0, salt, temp, dist, dens, theta
    integer, dimension(:), pointer:: density_nodes
    integer ele, i, node, n_sediment_fields, f
    
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

       if (have_option(trim(option_path)//'/generic_scalar_field_dependency')) then
          do f = 1, option_count(trim(option_path)//'/generic_scalar_field_dependency')
             dep_option_path=trim(option_path)//'/generic_scalar_field_dependency['//int2str(f-1)//']'
             call get_option(trim(dep_option_path)//'/name', sfield_name)
             call get_option(trim(dep_option_path)//'/reference_value', T0)
             call get_option(trim(dep_option_path)//'/expansion_coefficient', gamma)
             T => extract_scalar_field(state, trim(sfield_name))
             oldT => extract_scalar_field(state, "Old"//trim(sfield_name))
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
          end do
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
       call zero(sedimentdensity)

       n_sediment_fields = get_n_sediment_fields()

       do i=1,n_sediment_fields
       
          call get_sediment_item(state, i, S)

          call get_sediment_item(state, i, 'submerged_specific_gravity', gamma)
          
          gamma = gamma * rho_0

          oldS => extract_scalar_field(state, &
               "Old"//trim(S%name))
          
          ! deltaS=theta*S+(1-theta)*oldS-S0
          call remap_field(S, remapS)
          call set(deltaS, remapS)
          call scale(deltaS, theta)
          
          call remap_field(oldS, remapS)
          call addto(deltaS, remapS, 1.0-theta)
          ! density=density+gamma*deltaS
          call addto(sedimentdensity, deltaS, scale=gamma)

       end do
       
       call addto(density,sedimentdensity)

       call deallocate(deltaS)
       call deallocate(remapS)
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
  
  subroutine compressible_eos(state, density, pressure, drhodp)

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout), optional :: density, pressure, drhodp

    character(len=OPTION_PATH_LEN) :: eos_path
    type(scalar_field) :: drhodp_local

    ewrite(1,*) 'Entering compressible_eos'

    if (present(drhodp)) then
      drhodp_local=drhodp
      if (present(density)) then
        assert(drhodp%mesh==density%mesh)
      end if
      if (present(pressure)) then
        assert(drhodp%mesh==pressure%mesh)
      end if
    else if (present(density)) then
      call allocate(drhodp_local, density%mesh, 'Localdrhop')
    else if (present(pressure)) then
      call allocate(drhodp_local, pressure%mesh, 'Localdrhop')
    else
      FLAbort("No point in being in here if you don't want anything out.")
    end if

    eos_path = trim(state%option_path)//'/equation_of_state'

    if(have_option(trim(eos_path)//'/compressible')) then
      
      ! each of the following compressible_eos_XXX() routines should always calculate drhodp
      ! (zero if density does not depend on pressure) and calculate density and
      ! pressure if present
      
      if(have_option(trim(eos_path)//'/compressible/stiffened_gas')) then
        
        ! standard stiffened gas eos
        
        call compressible_eos_stiffened_gas(state, eos_path, drhodp_local, &
            density=density, pressure=pressure)
      
      else if(have_option(trim(eos_path)//'/compressible/giraldo')) then
        
        ! Eq. of state commonly used in atmospheric applications. See
        ! Giraldo et. al., J. Comp. Phys., vol. 227 (2008), 3849-3877. 
        ! density= P_0/(R*T)*(P/P_0)^((R+c_v)/c_p)
        
        call compressible_eos_giraldo(state, eos_path, drhodp_local, &
            density=density, pressure=pressure)


      elseif(have_option(trim(eos_path)//'/compressible/foam')) then
        
        ! eos used in foam modelling
        
        call compressible_eos_foam(state, eos_path, drhodp_local, &
            density=density, pressure=pressure)

      end if
      
    else
    
      ! I presume we dont' actually want to be here
      FLAbort('Gone into compressible_eos without having equation_of_state/compressible')

    end if

    if(present(density)) then
      ewrite_minmax(density)
    end if

    if(present(pressure)) then
      ewrite_minmax(pressure)
    end if

    if(present(drhodp)) then      
      ewrite_minmax(drhodp)
    else
      call deallocate(drhodp_local)
    end if

  end subroutine compressible_eos
    
  subroutine compressible_eos_stiffened_gas(state, eos_path, drhodp, &
    density, pressure)
    ! Standard stiffened gas equation
    type(state_type), intent(inout) :: state
    character(len=*), intent(in):: eos_path
    type(scalar_field), intent(inout) :: drhodp
    type(scalar_field), intent(inout), optional :: density, pressure
    
    !locals
    integer :: stat, gstat, cstat
    type(scalar_field), pointer :: pressure_local, energy_local, density_local
    real :: reference_density, ratio_specific_heats
    real :: bulk_sound_speed_squared, atmospheric_pressure
    type(scalar_field) :: energy_remap, pressure_remap, density_remap
    logical :: incompressible
    
    call get_option(trim(eos_path)//'/compressible/stiffened_gas/reference_density', &
                        reference_density, default=0.0)
        
    call get_option(trim(eos_path)//'/compressible/stiffened_gas/ratio_specific_heats', &
                    ratio_specific_heats, stat=gstat)
    if(gstat/=0) then
      ratio_specific_heats=1.0
    end if
    
    call get_option(trim(eos_path)//'/compressible/stiffened_gas/bulk_sound_speed_squared', &
                    bulk_sound_speed_squared, stat=cstat)
    if(cstat/=0) then
      bulk_sound_speed_squared=0.0
    end if
    
    incompressible = ((gstat/=0).and.(cstat/=0))
    if(incompressible) then
      ewrite(0,*) "Selected compressible eos but not specified a bulk_sound_speed_squared or a ratio_specific_heats."
    end if
    
    call zero(drhodp)
    
    if(.not.incompressible) then
      energy_local=>extract_scalar_field(state,'InternalEnergy',stat=stat)
      ! drhodp = 1.0/( bulk_sound_speed_squared + (ratio_specific_heats - 1.0)*energy )
      if((stat==0).and.(gstat==0)) then   ! we have an internal energy field and we want to use it
        call allocate(energy_remap, drhodp%mesh, 'RemappedInternalEnergy')
        call remap_field(energy_local, energy_remap)
        
        call addto(drhodp, energy_remap, (ratio_specific_heats-1.0))
        
        call deallocate(energy_remap)
      end if
      call addto(drhodp, bulk_sound_speed_squared)
      call invert(drhodp)
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
          assert(density%mesh==drhodp%mesh)
        
          ! density = drhodp*(pressure_local + atmospheric_pressure
          !                  + bulk_sound_speed_squared*reference_density)
          call get_option(trim(pressure_local%option_path)//'/prognostic/atmospheric_pressure', &
                          atmospheric_pressure, default=0.0)
          
          call allocate(pressure_remap, drhodp%mesh, "RemappedPressure")
          call remap_field(pressure_local, pressure_remap)
          
          call set(density, reference_density*bulk_sound_speed_squared + atmospheric_pressure)
          call addto(density, pressure_remap)
          call scale(density, drhodp)
          
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
          assert(pressure%mesh==drhodp%mesh)
          
          ! pressure = density_local/drhodp &
          !          - bulk_sound_speed_squared*reference_density
          
          call allocate(density_remap, drhodp%mesh, "RemappedDensity")
          call remap_field(density_local, density_remap)
          
          call set(pressure, drhodp)
          call invert(pressure)
          call scale(pressure, density_remap)
          call addto(pressure, -bulk_sound_speed_squared*reference_density)
          
          call deallocate(density_remap)
        else
          FLExit('No Density in material_phase::'//trim(state%name))
        end if
      end if
    end if

  end subroutine compressible_eos_stiffened_gas
  
  subroutine compressible_eos_giraldo(state, eos_path, drhodp, &
    density, pressure)
    ! Eq. of state commonly used in atmospheric applications. See
    ! Giraldo et. al., J. Comp. Phys., vol. 227 (2008), 3849-3877. 
    ! density= P_0/(R*T)*(P/P_0)^((R+c_v)/c_p)
    type(state_type), intent(inout) :: state
    character(len=*), intent(in):: eos_path
    type(scalar_field), intent(inout) :: drhodp
    type(scalar_field), intent(inout), optional :: density, pressure
      
    ! locals
    integer :: stat, gstat, cstat, pstat, tstat
    type(scalar_field), pointer :: pressure_local, density_local, temperature_local
    real :: reference_density, p_0, c_p, c_v
    real :: drhodp_node, power
    real :: R
    type(scalar_field) :: pressure_remap, density_remap, temperature_remap
    logical :: incompressible
    integer :: node
    
    call get_option(trim(eos_path)//'/compressible/giraldo/reference_pressure', &
                    p_0, default=1.0e5)
    
    call get_option(trim(eos_path)//'/compressible/giraldo/C_P', &
                    c_p, stat=gstat)
    if(gstat/=0) then
      c_p=1.0
    end if
    
    call get_option(trim(eos_path)//'/compressible/giraldo/C_V', &
                    c_v, stat=cstat)
    if(cstat/=0) then
      c_v=1.0
    end if

    R=c_p-c_v
    
    incompressible = ((gstat/=0).or.(cstat/=0))
    if(incompressible) then
      ewrite(0,*) "Selected compressible eos but not specified either C_P or C_V."
    end if
    
    call zero(drhodp)
    
    if(.not.incompressible) then
      pressure_local=>extract_scalar_field(state,'Pressure',stat=pstat)
      temperature_local=>extract_scalar_field(state,'Temperature',stat=tstat)
      if ((pstat==0).and.(tstat==0)) then
        ! drhodp = ((R+c_v)/c_p)*1.0/( R*T) * (P/P_0)^((R+c_v-c_p)/c_p)
        call allocate(pressure_remap, drhodp%mesh, 'RemappedPressure')
        call remap_field(pressure_local, pressure_remap)
        call allocate(temperature_remap, drhodp%mesh, 'RemappedTemperature')
        call remap_field(temperature_local, temperature_remap)

        power=(R+c_v-c_p)/c_p
        do node=1,node_count(drhodp)
          drhodp_node=((c_v+R)/c_p)*1.0/(R*node_val(temperature_remap,node))*(node_val(pressure_remap,node)/p_0)**(power)
          call set(drhodp, node, drhodp_node)
        end do
        
        call deallocate(temperature_remap)
      else
        FLExit('No Pressure or temperature in material_phase::'//trim(state%name))
      endif
    end if

    if(present(density)) then
      ! calculate the density
      ! density may equal density in state depending on how this
      ! subroutine is called
      if(incompressible) then
        ! density = reference_density
        call set(density, reference_density)
      else
        assert(density%mesh==drhodp%mesh)
        call set(density, pressure_remap)
        call scale(density, drhodp)
        call scale(density, 1.0/(1.0+power))
          
        call deallocate(pressure_remap)
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
          assert(pressure%mesh==drhodp%mesh)
          
          ! pressure = density_local/drhodp
          
          call allocate(density_remap, drhodp%mesh, "RemappedDensity")
          call remap_field(density_local, density_remap)
          
          call set(pressure, drhodp)
          call invert(pressure)
          call scale(pressure, density_remap)
          call scale(pressure, (1.0+power))
          
          call deallocate(density_remap)
        else
          FLExit('No Density in material_phase::'//trim(state%name))
        end if
      end if
    end if

  end subroutine compressible_eos_giraldo

  subroutine compressible_eos_foam(state, eos_path, drhodp, &
    density, pressure)
    ! Foam EoS Used with compressible simulations of liquid drainage in foams.
    ! It describes the liquid content in the foam as the product of the  Plateau 
    ! border cross sectional area and the local Plateau  border length per unit volume (lambda).
    type(state_type), intent(inout) :: state
    character(len=*), intent(in):: eos_path
    type(scalar_field), intent(inout) :: drhodp
    type(scalar_field), intent(inout), optional :: density, pressure

    ! locals
    integer :: pstat, dstat
    type(scalar_field), pointer :: pressure_local, density_local, drainagelambda_local
    real :: atmospheric_pressure
    type(scalar_field) :: pressure_remap, density_remap, drainagelambda_remap

    call zero(drhodp)

    pressure_local => extract_scalar_field(state,'Pressure', stat=pstat)

    drainagelambda_local => extract_scalar_field(state,'DrainageLambda')

    call allocate(drainagelambda_remap, drhodp%mesh, 'RemappedDrainageLambda')
    call remap_field(drainagelambda_local, drainagelambda_remap)

    call addto(drhodp, drainagelambda_remap)

    call deallocate(drainagelambda_remap)

    if(present(density)) then
      if (pstat==0) then
        assert(density%mesh==drhodp%mesh)

        call get_option(trim(pressure_local%option_path)//'/prognostic/atmospheric_pressure', &
                        atmospheric_pressure, default=0.0)

        call allocate(pressure_remap, drhodp%mesh, "RemappedPressure")
        call remap_field(pressure_local, pressure_remap)

        call set(density, atmospheric_pressure)
        call addto(density, pressure_remap)
        call scale(density, drhodp)

        call deallocate(pressure_remap)
      else
        FLExit('No Pressure in material_phase::'//trim(state%name))
      end if
    end if

    if(present(pressure)) then
      density_local=>extract_scalar_field(state,'Density',stat=dstat)
      if (dstat==0) then
        assert(pressure%mesh==drhodp%mesh)

        call get_option(trim(pressure_local%option_path)//'/prognostic/atmospheric_pressure', &
                        atmospheric_pressure, default=0.0)

        call allocate(density_remap, drhodp%mesh, "RemappedDensity")
        call remap_field(density_local, density_remap)

        call set(pressure, drhodp)
        call invert(pressure)
        call scale(pressure, density_remap)

        call deallocate(density_remap)
      else
        FLExit('No Density in material_phase::'//trim(state%name))
      end if
    end if
        
  end subroutine compressible_eos_foam
        
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

      if(have_option(trim(eos_path)//'/compressible/stiffened_gas')) then
        call get_option(trim(eos_path)//'/compressible/stiffened_gas/reference_density', &
                        reference_density, default=0.0)
        call get_option(trim(eos_path)//'/compressible/stiffened_gas/ratio_specific_heats', &
                        ratio_specific_heats, stat=gstat)
        if(gstat/=0) then
          ratio_specific_heats=1.0
        end if
        call get_option(trim(eos_path)//'/compressible/stiffened_gas/bulk_sound_speed_squared', &
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
      ewrite_minmax(materialdensity)
    end if

    if(present(materialpressure)) then
      ewrite_minmax(materialpressure)
    end if

    if(present(materialdrhodp)) then
      materialdrhodp%val=drhodp%val
      ewrite_minmax(materialdrhodp)
    end if

    call deallocate(drhodp)

  end subroutine compressible_material_eos
  
end module equation_of_state
