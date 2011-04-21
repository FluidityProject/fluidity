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

module momentum_diagnostic_fields
  use FLDebug
  use equation_of_state
  use fields
  use fefields
  use state_module
  use spud
  use state_module
  use state_fields_module
  use sparsity_patterns_meshes
  use solvers
  use field_priority_lists
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use multimaterial_module
  use multiphase_module
  use diagnostic_fields_wrapper_new
  implicit none

  interface calculate_densities
    module procedure calculate_densities_single_state, calculate_densities_multiple_states
  end interface

  private
  public :: calculate_momentum_diagnostics, calculate_densities
  public :: calculate_diagnostic_velocity

contains

  subroutine calculate_momentum_diagnostics(state, istate, submaterials, submaterials_istate, u_name)
    !< A subroutine to group together all the diagnostic calculations that
    !< must happen before a momentum solve.
  
    type(state_type), dimension(:), intent(inout) :: state
    integer, intent(in) :: istate
    ! An array of submaterials of the current phase in state(istate).
    type(state_type), dimension(:), intent(inout) :: submaterials
    ! The index of the current phase (i.e. state(istate)) in the submaterials array
    integer :: submaterials_istate
    ! name of the u field, either "Velocity" or "Momentum"
    character(len=*), intent(in) :: u_name
    
    ! Local variables  
    type(scalar_field), pointer :: bulk_density, buoyancy_density, sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield
    
    integer :: stat
    logical :: gravity, diagnostic
    
    ewrite(1,*) 'Entering calculate_momentum_diagnostics'
    
    ! This needs to be done first or none of the following multimaterial algorithms will work...
    call calculate_diagnostic_material_volume_fraction(submaterials)
    call calculate_diagnostic_phase_volume_fraction(state)

    ! Calculate the density according to the eos... do the buoyancy density and the density
    ! at the same time to save computations
    ! don't calculate buoyancy if no gravity
    gravity = have_option("/physical_parameters/gravity")
    bulk_density => extract_scalar_field(submaterials(submaterials_istate), 'Density', stat)
    diagnostic = .false.
    if (stat==0) diagnostic = have_option(trim(bulk_density%option_path)//'/diagnostic')
    if(diagnostic.and.gravity) then
      buoyancy_density => extract_scalar_field(submaterials(submaterials_istate), trim(u_name)//'BuoyancyDensity')
      call calculate_densities(submaterials,&
                                buoyancy_density=buoyancy_density, &
                                bulk_density=bulk_density, &
                                momentum_diagnostic=.true.)
    else if(diagnostic) then
      call calculate_densities(submaterials,&
                                bulk_density=bulk_density, &
                                momentum_diagnostic=.true.)
    else if(gravity) then
      buoyancy_density => extract_scalar_field(submaterials(submaterials_istate), trim(u_name)//'BuoyancyDensity')
      call calculate_densities(submaterials,&
                                buoyancy_density=buoyancy_density, &
                                momentum_diagnostic=.true.)
    end if
    
    vfield => extract_vector_field(submaterials(submaterials_istate), trim(u_name)//"Absorption", stat = stat)
    if(stat == 0) then
      if(have_option(trim(vfield%option_path) // "/diagnostic")) then
        call calculate_diagnostic_variable(submaterials, submaterials_istate, vfield)
      end if
    end if

    vfield => extract_vector_field(submaterials(submaterials_istate), trim(u_name)//"Source", stat = stat)
    if(stat == 0) then
      if(have_option(trim(vfield%option_path) // "/diagnostic")) then
        call calculate_diagnostic_variable(state, istate, vfield)
      end if
    end if

    tfield => extract_tensor_field(submaterials(submaterials_istate),'Viscosity',stat)
    if (stat==0) then
      diagnostic = have_option(trim(tfield%option_path)//'/diagnostic')
      if(diagnostic) then
        call calculate_diagnostic_variable(submaterials, submaterials_istate, tfield)
      end if
    end if

    tfield => extract_tensor_field(submaterials(submaterials_istate), trim(u_name)//'SurfaceTension', stat)
    if(stat==0) then
      diagnostic = have_option(trim(tfield%option_path)//'/diagnostic')
      if(diagnostic) then
        call calculate_surfacetension(submaterials, tfield)
      end if
    end if

    ! diagnostic Pressure (only for compressible) calculated from
    ! Density and InternalEnergie via compressible eos
    sfield => extract_scalar_field(submaterials(submaterials_istate), 'Pressure', stat)
    if(stat==0) then
      diagnostic = have_option(trim(sfield%option_path)//'/diagnostic')
      if(diagnostic) then
        call calculate_diagnostic_pressure(submaterials(submaterials_istate), sfield)
      end if
    end if

    ewrite(1,*) 'Exiting calculate_momentum_diagnostics'
    
  end subroutine calculate_momentum_diagnostics

  subroutine calculate_densities_single_state(state, buoyancy_density, bulk_density, &
                                              momentum_diagnostic)
  
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout), optional, target :: buoyancy_density
    type(scalar_field), intent(inout), optional, target :: bulk_density
    logical, intent(in), optional :: momentum_diagnostic
    
    type(state_type), dimension(1) :: states
  
    states = (/state/)
    call calculate_densities(states, buoyancy_density=buoyancy_density, bulk_density=bulk_density, &
                             momentum_diagnostic=momentum_diagnostic)
    state = states(1)
  
  end subroutine calculate_densities_single_state

  subroutine calculate_densities_multiple_states(state, buoyancy_density, bulk_density, &
                                                 momentum_diagnostic)
  
    type(state_type), dimension(:), intent(inout) :: state
    type(scalar_field), intent(inout), optional, target :: buoyancy_density
    type(scalar_field), intent(inout), optional, target :: bulk_density
    logical, intent(in), optional ::momentum_diagnostic
    
    type(scalar_field) :: eosdensity
    type(scalar_field), pointer :: tmpdensity
    type(scalar_field) :: bulksumvolumefractionsbound
    type(scalar_field) :: buoyancysumvolumefractionsbound
    type(mesh_type), pointer :: mesh
    integer, dimension(size(state)) :: state_order
    logical :: subtract_out_hydrostatic, multimaterial
    character(len=OPTION_PATH_LEN) :: option_path
    integer :: subtract_count, materialvolumefraction_count
    real :: hydrostatic_rho0, reference_density
    integer :: i, stat
    
    logical :: boussinesq
    type(vector_field), pointer :: velocity
    real :: boussinesq_rho0
    
    if(.not.present(buoyancy_density).and..not.present(bulk_density)) then
      ! coding error
      FLAbort("No point calling me if I don't have anything to do.")
    end if
    
    if(present(buoyancy_density)) call zero(buoyancy_density)
    if(present(bulk_density)) call zero(bulk_density)
    
    if(present(bulk_density)) then
      mesh => bulk_density%mesh
    else
      mesh => buoyancy_density%mesh
    end if
    
    multimaterial = .false.
    materialvolumefraction_count = 0
    subtract_count = 0
    if(size(state)>1) then
      do i = 1, size(state)
        if(has_scalar_field(state(i), "MaterialVolumeFraction")) then
          materialvolumefraction_count = materialvolumefraction_count + 1
        end if
        
        option_path='/material_phase::'//trim(state(i)%name)//'/equation_of_state'
        
        subtract_count = subtract_count + &
            option_count(trim(option_path)//'/fluids/linear/subtract_out_hydrostatic_level') + &
            option_count(trim(option_path)//'/fluids/ocean_pade_approximation')
        
      end do
      if(size(state)/=materialvolumefraction_count) then
        FLExit("Multiple material_phases but not all of them have MaterialVolumeFractions.")
      end if
      if(subtract_count>1) then
        FLExit("You can only select one material_phase to use the reference_density from to subtract out the hydrostatic level.")
      end if
      
      multimaterial = .true.
      
      ! allocate a bounding field for the volume fractions
      if(present(bulk_density)) then
        call allocate(bulksumvolumefractionsbound, mesh, "SumMaterialVolumeFractionsBound")
        call set(bulksumvolumefractionsbound, 1.0)
      end if
      if(present(buoyancy_density)) then
        call allocate(buoyancysumvolumefractionsbound, mesh, "SumMaterialVolumeFractionsBound")
        call set(buoyancysumvolumefractionsbound, 1.0)
      end if

      ! get the order in which states should be processed      
      call order_states_priority(state, state_order)

      ! this needs to be done first or none of the following multimaterial algorithms will work...
      call calculate_diagnostic_material_volume_fraction(state)
    else
      assert(size(state_order)==1)
      ! set up a dummy state ordering for the single material case
      state_order(1) = 1      
    end if
    
    boussinesq = .false.
    hydrostatic_rho0 = 0.0
    state_loop: do i = 1, size(state)
  
      option_path='/material_phase::'//trim(state(state_order(i))%name)//'/equation_of_state'
      
      if(have_option(trim(option_path)//'/fluids')) then
        ! we have a fluids eos
      
        subtract_out_hydrostatic = &
          have_option(trim(option_path)//'/fluids/linear/subtract_out_hydrostatic_level') .or. &
          have_option(trim(option_path)//'/fluids/ocean_pade_approximation')
        
        call allocate(eosdensity, mesh, "LocalPerturbationDensity")
        
        call calculate_perturbation_density(state(state_order(i)), eosdensity, reference_density)
        
        if(multimaterial) then
          ! if multimaterial we have to subtract out a single reference density at the end
          ! rather than one per material so add it in always for now
          call addto(eosdensity, reference_density)
        end if
        
        if(present(buoyancy_density)) then
          if(multimaterial) then
            if(subtract_out_hydrostatic) then
              ! if multimaterial we have to subtract out a single global value at the end
              ! so save it for now
              hydrostatic_rho0 = reference_density
            end if
            call add_scaled_material_property(state(state_order(i)), buoyancy_density, eosdensity, &
                                              sumvolumefractionsbound=buoyancysumvolumefractionsbound, &
                                              momentum_diagnostic=momentum_diagnostic)
          else
            call set(buoyancy_density, eosdensity)
            if(.not.subtract_out_hydrostatic) then
              call addto(buoyancy_density, reference_density)
            end if
          end if
          
          ! find out if the velocity in this state is *the* (i.e. not aliased)
          ! prognostic one and if it's using a Boussinesq equation.
          ! if it is record the rho0 and it will be used later to scale
          ! the buoyancy density
          velocity => extract_vector_field(state(state_order(i)), "Velocity", stat)
          if(stat==0) then
            if(.not.aliased(velocity)) then
              if (have_option(trim(velocity%option_path)//"/prognostic/equation::Boussinesq")) then
                ! have we already found a state where the velocity was using Boussinesq?
                if(boussinesq) then
                  ! uh oh... looks like you're using multiphase... good luck with that...
                  ! everything here at the moment assumes a single prognostic velocity
                  FLExit("Two nonaliased velocities using equation type Boussinesq.  Don't know what to do.")
                end if
                boussinesq=.true.
                boussinesq_rho0 = reference_density
              end if
            end if
          end if
          
        end if
        
        if(present(bulk_density)) then
          if(multimaterial) then
            ! the perturbation density has already had the reference density added to it
            ! if you're multimaterial
            call add_scaled_material_property(state(state_order(i)), bulk_density, eosdensity, &
                                              sumvolumefractionsbound=bulksumvolumefractionsbound, &
                                              momentum_diagnostic=momentum_diagnostic)
          else
            call set(bulk_density, eosdensity)
            call addto(bulk_density, reference_density)
          end if
        end if
        
        call deallocate(eosdensity)
        
      else
        ! we don't have a fluids eos
        
        tmpdensity => extract_scalar_field(state(state_order(i)), "MaterialDensity", stat)
        if(stat==0) then
          if(multimaterial) then
            if(present(buoyancy_density)) then
              call add_scaled_material_property(state(state_order(i)), buoyancy_density, tmpdensity, &
                                                sumvolumefractionsbound=buoyancysumvolumefractionsbound, &
                                                momentum_diagnostic=momentum_diagnostic)
            end if
            if(present(bulk_density)) then
              call add_scaled_material_property(state(state_order(i)), bulk_density, tmpdensity, &
                                                sumvolumefractionsbound=bulksumvolumefractionsbound, &
                                                momentum_diagnostic=momentum_diagnostic)
            end if
          else
            if(present(buoyancy_density)) then
              call remap_field(tmpdensity, buoyancy_density)
            end if
            if(present(bulk_density)) then
              call remap_field(tmpdensity, bulk_density)
            end if
          end if
        else
          if(multimaterial) then
            FLExit("No multimaterial MaterialDensity or fluid eos provided")
          else
            if(have_option(trim(option_path)//'/compressible')) then
              call allocate(eosdensity, mesh, "LocalCompressibleEOSDensity")
            
              call compressible_eos(state(state_order(i)), density=eosdensity)
              
              if(present(bulk_density)) then
                call set(bulk_density, eosdensity)
              end if
              if(present(buoyancy_density)) then
                call set(buoyancy_density, eosdensity)
              end if
            
              call deallocate(eosdensity)
            else
              tmpdensity => extract_scalar_field(state(state_order(i)), "Density", stat)
              if(stat==0) then
                if(present(buoyancy_density)) then
                  call remap_field(tmpdensity, buoyancy_density)
                end if
                if(present(bulk_density)) then
                  call remap_field(tmpdensity, bulk_density)
                end if
              else
                if(present(buoyancy_density)) then
                  FLExit("You haven't provide enough information to set the buoyancy density.")
                end if
                if(present(bulk_density)) then
                  ! coding error... hopefully
                  FLAbort("How on Earth did you get here without a density?!")
                end if
              end if
            end if
          end if
        end if

      end if
        
    end do state_loop
    
    if(present(buoyancy_density)) then
      if(multimaterial) call addto(buoyancy_density, -hydrostatic_rho0)
      
      if(boussinesq) then
        ! the buoyancy density is being used in a Boussinesq eqn
        ! therefore it needs to be scaled by rho0:
        call scale(buoyancy_density, 1./boussinesq_rho0)
      end if
    end if

    if(multimaterial) then
      if(present(buoyancy_density)) then
        call deallocate(buoyancysumvolumefractionsbound)
      end if
      if(present(bulk_density)) then
        call deallocate(bulksumvolumefractionsbound)
      end if
    end if
  
  end subroutine calculate_densities_multiple_states

  subroutine calculate_diagnostic_pressure(state, pressure)
    ! diagnostic Pressure (only for compressible) calculated from
    ! Density and InternalEnergie via compressible eos
    type(state_type), intent(inout):: state
    type(scalar_field), intent(inout):: pressure

    ewrite(1,*) "In calculate_diagnostic_pressure"

    if (have_option(trim(state%option_path)//'/equation_of_state/compressible')) then
      call compressible_eos(state, pressure=pressure)
    else
      FLExit("Diagnostic pressure can only be used in combination with a compressible equation of state.")
    end if

    ewrite_minmax(pressure)

  end subroutine calculate_diagnostic_pressure

  subroutine calculate_diagnostic_velocity(state, u)
    ! diagnostic Velocity, only for compressible with prognostic Momentum
    ! calculated simply by Velocity=Momentum/Density
    type(state_type), intent(inout):: state
    type(vector_field), intent(inout):: u

    type(vector_field), pointer:: positions, momentum
    type(vector_field):: proj_rhs
    type(scalar_field), pointer:: density
    type(scalar_field):: lumped_mass 
    type(csr_matrix):: mass
    type(csr_sparsity), pointer:: mass_sparsity
    integer:: ele, i

    ewrite(1,*) "Inside calculate_diagnostic_velocity"

    positions => extract_vector_field(state, "Coordinate")
    ! there must be a prognostic momentum
    momentum => extract_vector_field(state, "Momentum")
    ! a density is also required
    density => extract_scalar_field(state, "Density")

    if (continuity(u)>=0) then

      call allocate(proj_rhs, u%dim, u%mesh, "DiagnosticVelocityProjectionRHS")
      call zero(proj_rhs)
      do ele=1, element_count(u)
        call calculate_diagnostic_velocity_cg_ele(proj_rhs, ele)
      end do

      if (have_option(trim(u%option_path)//"/diagnostic/solver")) then
        mass_sparsity => get_csr_sparsity_firstorder(state, u%mesh, u%mesh)
        call allocate(mass, mass_sparsity, name="DiagnosticVelocityMass")
        call compute_mass(positions, u%mesh, mass, density=density)
 
        call petsc_solve(u, mass, proj_rhs)

        call deallocate(mass)
      else
        ewrite(1,*) "No solver options specified - lumping mass to compute diagnostic velocity"
        call allocate(lumped_mass, u%mesh, name="DiagnosticVelocityLumpedMass")
        call compute_lumped_mass(positions, lumped_mass, density=density)
        call set(u, proj_rhs)
        do i=1, u%dim
          ewrite_minmax( u%val(i,:) )
        end do
        ewrite_minmax( lumped_mass )
        call inverse_scale(u, lumped_mass)
        call deallocate(lumped_mass)
      end if

      call deallocate(proj_rhs)

    else
      do ele=1, element_count(u)
        call calculate_diagnostic_velocity_dg_ele(ele)
      end do
    end if

    ewrite_minmax( density )
    do i=1, u%dim
      ewrite_minmax( u%val(i,:) )
    end do
    do i=1, momentum%dim
      ewrite_minmax( momentum%val(i,:) )
    end do

    contains

    subroutine calculate_diagnostic_velocity_cg_ele(proj_rhs, ele)
      type(vector_field), intent(inout):: proj_rhs
      integer, intent(in):: ele

      real, dimension(ele_ngi(u, ele)):: detwei

      call transform_to_physical(positions, ele, detwei)
      call addto(proj_rhs, ele_nodes(proj_rhs, ele), &
          shape_vector_rhs(ele_shape(proj_rhs, ele), &
              ele_val_at_quad(momentum, ele), detwei))

    end subroutine calculate_diagnostic_velocity_cg_ele

    subroutine calculate_diagnostic_velocity_dg_ele(ele)
      integer, intent(in):: ele

      real, dimension(ele_loc(u, ele), ele_loc(u, ele)):: local_mass
      real, dimension(ele_loc(u, ele), u%dim):: local_rhs
      real, dimension(ele_ngi(u, ele)):: detwei


      call transform_to_physical(positions, ele, detwei)
      ! this has to be transposed because of the way solve_mulitple (in actually fact
      ! dgesv from lapack) expects multiple rhs's
      local_rhs=transpose(shape_vector_rhs(ele_shape(u, ele),  &
          ele_val_at_quad(momentum, ele), detwei/ele_val_at_quad(density, ele)))

      ! this solve takes the rhs and overwrites it with the solution
      call solve(local_mass, local_rhs)

      ! now transpose it again
      call set(u, ele_nodes(u, ele), transpose(local_rhs))
      
    end subroutine calculate_diagnostic_velocity_dg_ele

  end subroutine calculate_diagnostic_velocity

  subroutine momentum_diagnostics_fields_check_options

    character(len=OPTION_PATH_LEN):: phase_path
    integer:: i
  
    do i=0, option_count('/material_phase')-1
       phase_path = '/material_phase[' // int2str(i) // ']'
       if (have_option(trim(phase_path)//'/scalar_field::Pressure/diagnostic')) then
         if (.not. have_option(trim(phase_path)//'/equation_of_state/compressible')) then
           FLExit("Diagnostic pressure can only be used in combination with a compressible equation of state.")
         end if
         if (have_option(trim(phase_path)//'/scalar_field::MaterialVolumeFraction')) then
           FLExit("Diagnostic pressure currently does not work with multi-material")
         end if

       end if
    end do

  end subroutine momentum_diagnostics_fields_check_options

end module momentum_diagnostic_fields
