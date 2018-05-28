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

module diagnostic_fields_wrapper
  !!< A module to link to diagnostic variable calculations.

  use fldebug
  use global_parameters, only: FIELD_NAME_LEN, timestep
  use futils
  use spud
  use parallel_tools
  use fetools
  use fields
  use sparse_matrices_fields
  use state_module
  use field_derivatives
  use field_options, only: do_not_recalculate
  use diagnostic_fields, only: calculate_diagnostic_variable
  use equation_of_state
  use multiphase_module
  use diagnostic_fields_matrices
  use multimaterial_module, only: calculate_material_mass, &
       calculate_bulk_material_pressure, &
       calculate_sum_material_volume_fractions, &
       calculate_material_volume
  use tidal_module, only: calculate_diagnostic_equilibrium_pressure
  use free_surface_module, only: calculate_diagnostic_free_surface, &
       calculate_diagnostic_wettingdrying_alpha
  use vorticity_diagnostics
  use momentum_diagnostic_fields
  use sediment_diagnostics
  use dqmom
  use geostrophic_pressure
  
  implicit none

  private
  
  public :: calculate_diagnostic_variables

contains

  subroutine calculate_diagnostic_variables(state, exclude_nonrecalculated)
    !!< Updates diagnostic fields in the supplied states.

    type(state_type), dimension(:) :: state
    logical, intent(in), optional :: exclude_nonrecalculated

    integer :: i,stat
    type(scalar_field), pointer :: s_field
    type(vector_field), pointer :: v_field
    logical :: diagnostic

    ! An array of submaterials of the current phase in state(istate).
    type(state_type), dimension(:), pointer :: submaterials
    
    ewrite(1, *) "In calculate_diagnostic_variables"
 
    do i = 1, size(state)

       ! start of fields that can be called through the generic calculate_diagnostic_variable
       ! interface, i.e. - those that only require things available in f90modules

       s_field => extract_scalar_field(state(i), "CFLNumber", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "CFLNumber", &
             & s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "GridReynoldsNumber", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "GridReynoldsNumber", &
             & s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "GridPecletNumber", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "GridPecletNumber", &
             & s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "ControlVolumeCFLNumber", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "ControlVolumeCFLNumber", &
             & s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "DG_CourantNumber", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "DG_CourantNumber", &
             & s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "CVMaterialDensityCFLNumber", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "CVMaterialDensityCFLNumber", &
             & s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "KineticEnergyDensity", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "KineticEnergyDensity", &
             & s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i),  "HorizontalVelocityDivergence", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "HorizontalVelocityDivergence", s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), &
         & "VelocityDivergence", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "VelocityDivergence", s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), &
         & "PerturbationDensity", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           ! this routine returns the density used in the buoyancy term, which we're not really interested in
           ! but it computes the PerturbationDensity as a side effect. Note that this means it will happen twice
           ! as it will be recalculated at the beginning of Momemtum_Equation after the Temperature and Salinity
           ! fields have been solved for.
           call calculate_perturbation_density(state(i), s_field)
         end if
       end if

       ! this diagnostic field depends on PerturbationDensity
       s_field => extract_scalar_field(state(i), "GravitationalPotentialEnergyDensity", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "GravitationalPotentialEnergyDensity", s_field)
         end if
       end if
       
       ! this diagnostic field depends on PerturbationDensity
       s_field => extract_scalar_field(state(i), "IsopycnalCoordinate", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "IsopycnalCoordinate", &
             & s_field)
         end if
       end if
       
       ! Must be calculated after IsopycnalCoordinate
       s_field => extract_scalar_field(state(i), "BackgroundPotentialEnergyDensity", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "BackgroundPotentialEnergyDensity", s_field)
         end if
       end if
       
       v_field => extract_vector_field(state(i), "InnerElementFullVelocity", stat)
       if(stat == 0) then
         if(recalculate(trim(v_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "InnerElementFullVelocity", &
             & v_field)
         end if
       end if

       ! Must be calculated after InnerElementFullVelocity
       v_field => extract_vector_field(state(i), "InnerElementFullVorticity", stat)
       if(stat == 0) then
         if(recalculate(trim(v_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "InnerElementFullVorticity", &
             & v_field)
         end if
       end if

       v_field => extract_vector_field(state(i), "InnerElementVorticity", stat)
       if(stat == 0) then
         if(recalculate(trim(v_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "InnerElementVorticity", &
             & v_field)
         end if
       end if

       v_field => extract_vector_field(state(i), "DgMappedVelocity", stat)
       if(stat == 0) then
         if(recalculate(trim(v_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "DgMappedVelocity", &
             & v_field)
         end if
       end if

       v_field => extract_vector_field(state(i), "DgMappedVorticity", stat)
       if(stat == 0) then
         if(recalculate(trim(v_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "DgMappedVorticity", &
             & v_field)
         end if
       end if
       
       s_field => extract_scalar_field(state(i), "HorizontalStreamFunction", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "HorizontalStreamFunction", s_field)
         end if
       end if
       
       s_field => extract_scalar_field(state(i), "Speed", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "Speed", &
             & s_field)
         end if
       end if
       
       s_field => extract_scalar_field(state(i), "DiffusiveDissipation", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "DiffusiveDissipation", &
             & s_field)
         end if
       end if
       
       s_field => extract_scalar_field(state(i), "RichardsonNumber", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "RichardsonNumber", &
             & s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "StreamFunction", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "StreamFunction", s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "MultiplyConnectedStreamFunction", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "MultiplyConnectedStreamFunction", s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "Time", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "Time", s_field)
         end if
       end if

       v_field => extract_vector_field(state(i), "LinearMomentum", stat)
       if(stat == 0) then
         if(recalculate(trim(v_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "LinearMomentum", v_field)
         end if
       end if

       v_field => extract_vector_field(state(i), "DiagnosticCoordinate", stat)
       if(stat == 0) then
         if(recalculate(trim(v_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "DiagnosticCoordinate", v_field)
         end if
       end if

       v_field => extract_vector_field(state(i), "BedShearStress", stat)  
       if(stat == 0) then  
         if(recalculate(trim(v_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "BedShearStress", v_field)  
         end if
       end if

       v_field => extract_vector_field(state(i), "MaxBedShearStress", stat)  
       if(stat == 0) then  
         if(recalculate(trim(v_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "MaxBedShearStress", v_field)  
         end if
       end if

       s_field => extract_scalar_field(state(i), "GalerkinProjection", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "GalerkinProjection", s_field)
         end if
       end if

       v_field => extract_vector_field(state(i), "GalerkinProjection", stat)
       if(stat == 0) then
         if(recalculate(trim(v_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "GalerkinProjection", v_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "UniversalNumber", stat)
       if(stat == 0) then
          call calculate_diagnostic_variable(state(i), "UniversalNumber", s_field)
       end if

       s_field => extract_scalar_field(state(i), "NodeOwner", stat)
       if(stat == 0) then
          call calculate_diagnostic_variable(state(i), "NodeOwner", s_field)
       end if

       ! end of fields that can be called through the generic calculate_diagnostic_variable
       ! interface, i.e. - those that only require things available in f90modules

       ! start of fields that cannot be called through the generic calculate_diagnostic_variable
       ! interface, i.e. - those that need things from assemble
       s_field => extract_scalar_field(state(i), "ControlVolumeDivergence", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_divergence_cv(state(i), s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "FiniteElementDivergence", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_divergence_fe(state(i), s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "MaterialVolumeFraction", stat)
       if(stat == 0) then
         diagnostic = have_option(trim(s_field%option_path)//"/diagnostic")
         if(diagnostic .and. .not.(aliased(s_field))) then
           if(recalculate(trim(s_field%option_path))) then
             call calculate_sum_material_volume_fractions(state, s_field)
             call scale(s_field, -1.0)
             call addto(s_field, 1.0)
           end if
         end if
       end if
       
       s_field => extract_scalar_field(state(i), "MaterialMass", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_material_mass(state(i), s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "MaterialVolume", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_material_volume(state(i), s_field)
         end if
       end if
       
       s_field => extract_scalar_field(state(i), "MaterialDensity", stat)
       if(stat == 0) then
         diagnostic = have_option(trim(s_field%option_path)//"/diagnostic")
         if(diagnostic .and. .not.(aliased(s_field))) then
           if(recalculate(trim(s_field%option_path))) then
            call calculate_densities(state(i), bulk_density=s_field)
           end if
         end if
       end if

       s_field => extract_scalar_field(state(i), "Density", stat)
       if(stat == 0) then
         diagnostic = have_option(trim(s_field%option_path)//"/diagnostic")
         if(diagnostic .and. .not.(aliased(s_field))) then
           if(recalculate(trim(s_field%option_path))) then
             if(option_count("/material_phase/vector_field::Velocity/prognostic") > 1) then 
               call get_phase_submaterials(state, i, submaterials)
               call calculate_densities(submaterials, bulk_density=s_field)
               deallocate(submaterials)
             else
               call calculate_densities(state, bulk_density=s_field)
             end if
           end if
         end if
       end if

       s_field => extract_scalar_field(state(i), "MaterialEOSDensity", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call compressible_material_eos(state(i), materialdensity=s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "MaterialPressure", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call compressible_material_eos(state(i), materialpressure=s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "BulkMaterialPressure", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_bulk_material_pressure(state, s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "SumMaterialVolumeFractions", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_sum_material_volume_fractions(state, s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "FreeSurface", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_free_surface(state(i), s_field)
         end if
       end if
       
       s_field => extract_scalar_field(state(i), "WettingDryingAlpha", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_wettingdrying_alpha(state(i), s_field)
         end if
       end if

       s_field => extract_scalar_field(state(i), "EquilibriumPressure", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_equilibrium_pressure(state(i), s_field)
         end if
       end if

       v_field => extract_vector_field(state(i), "ControlVolumeDivergenceTransposed", stat)
       if(stat == 0) then
         if(recalculate(trim(v_field%option_path))) then
           call calculate_div_t_cv(state(i), v_field)
         end if
       end if

       v_field => extract_vector_field(state(i), "FiniteElementGradient", stat)
       if(stat == 0) then
         if(recalculate(trim(v_field%option_path))) then
           call calculate_grad_fe(state(i), v_field)
         end if
       end if

       v_field => extract_vector_field(state(i), "FiniteElementDivergenceTransposed", stat)
       if(stat == 0) then
         if(recalculate(trim(v_field%option_path))) then
           call calculate_div_t_fe(state(i), v_field)
         end if
       end if
              
       v_field => extract_vector_field(state(i), "PlanetaryVorticity", stat)
       if(stat == 0) then
         if(recalculate(trim(v_field%option_path))) then
           call calculate_planetary_vorticity(state(i), v_field)
         end if
       end if
       
       v_field => extract_vector_field(state(i), "AbsoluteVorticity", stat)
       if(stat == 0) then
         if(recalculate(trim(v_field%option_path))) then
           call calculate_absolute_vorticity(state(i), v_field)
         end if
       end if
       
       s_field => extract_scalar_field(state(i), "PotentialVorticity", stat)
       if(stat == 0) then
         diagnostic = have_option(trim(s_field%option_path)//"/diagnostic")
         if(diagnostic .and. recalculate(trim(s_field%option_path))) then
           call calculate_potential_vorticity(state(i), s_field)
         end if
       end if
       
       s_field => extract_scalar_field(state(i), "RelativePotentialVorticity", stat)
       if(stat == 0) then
         if(recalculate(trim(s_field%option_path))) then
           call calculate_relative_potential_vorticity(state(i), s_field)
         end if
       end if
       ! End of vorticity diagnostics

       ! Start of sediment diagnostics.
       if (have_option("/material_phase[0]/sediment")) then
          call calculate_sediment_sinking_velocity(state(i))
          call calculate_sediment_active_layer_d50(state(i))
          call calculate_sediment_active_layer_sigma(state(i))
          call calculate_sediment_active_layer_volume_fractions(state(i))
       end if
       ! End of sediment diagnostics.
       
       ! Start of population balance diagnostics.
       call dqmom_calculate_moments(state(i))
       call dqmom_calculate_statistics(state(i))
       ! End of population balance diagnostics.

       ! Multiphase-related diagnostic fields
       s_field => extract_scalar_field(state(i), "PhaseVolumeFraction", stat)
       if(stat == 0) then
         diagnostic = have_option(trim(s_field%option_path)//"/diagnostic")
         if(diagnostic .and. .not.(aliased(s_field))) then
           if(recalculate(trim(s_field%option_path))) then
            call calculate_diagnostic_phase_volume_fraction(state)
           end if
         end if
       end if

       s_field => extract_scalar_field(state(i), "SumVelocityDivergence", stat)
       if(stat == 0) then
         ! Check that we are running a multiphase simulation
         if(option_count("/material_phase/vector_field::Velocity/prognostic") > 1) then 
            diagnostic = have_option(trim(s_field%option_path)//"/diagnostic")
            if(diagnostic .and. .not.(aliased(s_field))) then
               if(recalculate(trim(s_field%option_path))) then
                  call calculate_sum_velocity_divergence(state, s_field)
               end if
            end if
         else
            FLExit("The SumVelocityDivergence field is only used in multiphase simulations.")
         end if
       end if
       
       s_field => extract_scalar_field(state(i), "CompressibleContinuityResidual", stat)
       if(stat == 0) then
         ! Check that we are running a compressible multiphase simulation
         if(option_count("/material_phase/vector_field::Velocity/prognostic") > 1 .and. option_count("/material_phase/equation_of_state/compressible") > 0) then 
            diagnostic = have_option(trim(s_field%option_path)//"/diagnostic")
            if(diagnostic .and. .not.(aliased(s_field))) then
               if(recalculate(trim(s_field%option_path))) then
                  call calculate_compressible_continuity_residual(state, s_field)
               end if
            end if
         else
            FLExit("The CompressibleContinuityResidual field is only used in compressible multiphase simulations.")
         end if
       end if

       ! end of fields that cannot be called through the generic
       ! calculate_diagnostic_variable interface, i.e. - those that need things
       ! higher than femtools in the build
       
       ! the following fields need to be here in case they are taking the difference with
       ! other diagnostic fields
       s_field => extract_scalar_field(state(i), "ScalarAbsoluteDifference", stat)  
       if(stat == 0) then  
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "AbsoluteDifference", s_field)  
         end if
       end if
       
       v_field => extract_vector_field(state(i), "VectorAbsoluteDifference", stat)  
       if(stat == 0) then  
         if(recalculate(trim(v_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "AbsoluteDifference", v_field)  
         end if
       end if

       s_field => extract_scalar_field(state(i), "AbsoluteDifference", stat)  
       if(stat == 0) then  
         if(recalculate(trim(s_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "AbsoluteDifference", s_field)  
         end if
       end if
       
       v_field => extract_vector_field(state(i), "AbsoluteDifference", stat)  
       if(stat == 0) then  
         if(recalculate(trim(v_field%option_path))) then
           call calculate_diagnostic_variable(state(i), "AbsoluteDifference", v_field)  
         end if
       end if

    end do
    
    ewrite(1, *) "Exiting calculate_diagnostic_variables"
    
  contains
  
    logical function recalculate(option_path)
      character(len=*) :: option_path
      
      recalculate = ((.not.present_and_true(exclude_nonrecalculated)).or. &
           (.not.do_not_recalculate(option_path)))
    
    end function recalculate

  end subroutine calculate_diagnostic_variables

end module diagnostic_fields_wrapper
