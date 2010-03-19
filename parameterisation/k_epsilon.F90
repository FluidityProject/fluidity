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
!    but WITHOUT ANY WARRANTY; without even the implied arranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module k_epsilon
  use quadrature
  use elements
  use field_derivatives
  use fields
  use sparse_matrices_fields
  use state_module
  use spud
  use global_parameters, only:   OPTION_PATH_LEN
  use state_fields_module
  use populate_state_module
  use boundary_conditions
  use fields_manipulation
  !use ieee_arithmetic
  use fetools
  use FLDebug

  implicit none

  private

  ! old turbulent kinetic energy, lengthscale, eddy viscosity, eddy diffusivity
  type(scalar_field), save :: tke_old, ll, EV, ED
  type(scalar_field), save :: surface_kk_values, surface_eps_values
  ! Minimum values of 2 fields to initialise, and empirical constants from flml
  real, save               :: eps_min, k_min, c_mu, c_eps_1, c_eps_2, sigma_eps, sigma_k
  integer, save            :: nNodes
  logical, save            :: initialised, calculate_bcs

  ! The following are the public subroutines
  public :: keps_init, keps_cleanup, keps_tke, keps_eps, keps_eddyvisc, keps_adapt_mesh, keps_check_options

  ! General plan is:
  !  - Init in populate_state.
  !  - If solve is about to do TKE, call keps_tke (which calculates P and sets source/absorption for solve).
  !  - If solve is about to do epsilon, call keps_eps (which fixes TKE surfaces, set source/absorption for solve).
  !  - After keps_eps solve, keps_eddyvisc recalculates the eddy viscosity and adds it to the viscosity field.
  !  - keps_eddyvisc also recalculates the lengthscale.
  !  - keps_adapt_options repopulates the fields after an adapt.
  !  - When done, clean-up.
  !
  ! TurbulentKineticEnergy and TurbulentDissipation need to have higher priority then other prognostic 
  ! fields such as temperature, velocity and salinity for this to work.
  ! Priority is set in preprocessor/Field_Priority_Lists

contains

!----------
!    - check we have the right fields (if not abort)
!    - initialise parameters based on options
!    - allocate space for optional fields, which are module level variables (to save passing them around)
!----------

subroutine keps_init(state)

    type(state_type), intent(inout) :: state
    type(scalar_field), pointer     :: scalarField

    ewrite(1,*)'Now in k_epsilon turbulence model - keps_init'

    ! Do we have the calculate_BCs option? Better have!
    calculate_bcs = &
    have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/calculate_boundaries")

    ! Allocate the temporary, module-level variables
    call keps_allocate_fields(state)

    ! Get the 5 model constants
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/C_mu',          C_mu, default = 0.09)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/C_eps_1',    c_eps_1, default = 1.44)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/C_eps_2',    c_eps_2, default = 1.92)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/sigma_k',     sigma_k, default = 1.0)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/sigma_eps', sigma_eps, default = 1.3)

    ewrite(1,*) "k-epsilon parameters"
    ewrite(1,*) "--------------------------------------------"
    ewrite(1,*) "c_mu: ",     c_mu
    ewrite(1,*) "c_eps_1: ",  c_eps_1
    ewrite(1,*) "c_eps_2: ",  c_eps_2
    ewrite(1,*) "sigma_k: ",  sigma_k
    ewrite(1,*) "sigma_eps: ",sigma_eps
    ewrite(1,*) "Calculating BCs: ", calculate_bcs
    ewrite(1,*) "--------------------------------------------"

    ! Minimum values: hard-coded. Must be something small and nonzero.
    k_min   = 7.6e-6
    eps_min = 1.e-9
    
    ! initialise 2 fields (k, epsilon) with minimum values
    scalarField => extract_scalar_field(state, "TurbulentKineticEnergy")
    call set(scalarField,k_min)
    scalarField => extract_scalar_field(state, "TurbulentDissipation")
    call set(scalarField,eps_min)

    ! initialise other fields
    call set(EV,1.e-6)
    call set(ED,1.e-6)
    call set(ll, min(k_min**1.5/eps_min, 10.0) ) ! need to set a minimum using some sensible function

    ewrite(1,*) "Internal fields initial values (minmax)"
    ewrite(1,*) "--------------------------------------------"
    ewrite_minmax(EV)
    ewrite_minmax(ED)
    ewrite_minmax(ll)
    ewrite(1,*) "--------------------------------------------"

    ! initialise surface
    if (calculate_bcs) then
        initialised = .false.
        ewrite(1,*) "initialised: ", initialised
        ewrite(1,*) "calling keps_init_surfaces"
        call keps_init_surfaces(state)
    end if
end subroutine keps_init

!------------------------------------------------------------------------------

subroutine keps_tke(state)

    type(state_type), intent(inout)   :: state
    type(scalar_field), pointer       :: source, absorption, kk, eps
    type(scalar_field), pointer       :: scalar_surface, scalarField
    type(vector_field), pointer       :: velocity, positions
    type(tensor_field), pointer       :: kk_diff, background_diff, visc
    type(tensor_field)                :: DU_DX  !, DU_DXT
    type(scalar_field), dimension(1)  :: grtemp
    type(scalar_field)                :: utemp
    integer                           :: i, ii, j, stat
    logical, allocatable, dimension(:):: derivs

    ewrite(1,*) "In keps_tke"

    positions  => extract_vector_field(state, "Coordinate")
    velocity   => extract_vector_field(state, "Velocity")
    visc       => extract_tensor_field(state, "Viscosity")
    source     => extract_scalar_field(state, "TurbulentKineticEnergySource")
    absorption => extract_scalar_field(state, "TurbulentKineticEnergyAbsorption")
    kk_diff    => extract_tensor_field(state, "TurbulentKineticEnergyDiffusivity")
    kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps        => extract_scalar_field(state, "TurbulentDissipation")

    ! clip at k_min and eps_min. We can't have negative values now can we.
    do i=1,nNodes
        call set(kk, i, max(kk%val(i), k_min) )
        call set(eps, i, max(eps%val(i), eps_min) )
    end do

    ! now allocate our temporary gradient fields
    call allocate(DU_DX, visc%mesh, "DUDX")
    allocate( derivs(velocity%dim) )
    derivs = .false.

    call zero(DU_DX)
    do i = 1, velocity%dim
        do j = 1, velocity%dim
            utemp  = extract_scalar_field_from_vector_field(velocity, i)
            grtemp = extract_scalar_field_from_tensor_field(DU_DX, i, j)
            if (i==j) derivs(j) = .true.
            call differentiate_field( utemp, positions, derivs, grtemp )
            derivs = .false.
        end do
    end do

    ewrite(1,*) "In kinetic energy production term loop"
    do ii = 1, nNodes
        do i = 1, velocity%dim
            do j = 1, velocity%dim
                call addto(source, ii, EV%val(ii) * ( DU_DX%val(i,j,ii) + DU_DX%val(j,i,ii) ) &
                * DU_DX%val(i,j,ii) )
            end do
        end do
        call set(absorption, ii, eps%val(ii) / kk%val(ii))
    end do

    call deallocate(DU_DX)
    deallocate(derivs)

    ! set diffusivity for KK - copied from Jon's update 28 Feb revision 12719
    call zero(kk_diff)
    background_diff => extract_tensor_field(state, "BackgroundDiffusivity")
    call set(kk_diff, kk_diff%dim, kk_diff%dim, EV, scale=1. / sigma_k)
    call addto(KK_diff, background_diff)

    ! puts the BC boundary values in surface_kk_values (module level variable):
    if (calculate_bcs) then
        call tke_bc(state)   ! (almost) zero Dirichlet BC: sets to k_min
        ! map these onto the actual BCs in kk - like GLS, but is it necessary if no fix_surface_values?
        scalar_surface => extract_surface_field(KK, 'tke_boundary', "value")
        call remap_field(surface_kk_values, scalar_surface)
    end if
    
    ! finally, we need a copy of the old TKE for eps, so grab it before we solve.
    do i=1,nNodes
        call set(tke_old, i, kk%val(i) )
    end do

    ewrite_minmax(source%val(:))
    ewrite_minmax(absorption%val(:))

    ! set source and absorption terms in optional output fields
    ! TKE source and absorption are ready for the solve (see Fluids.F90)
    scalarField => extract_scalar_field(state, "Source1", stat)
    if(stat == 0) then
       call set(scalarField, source)  
    end if
    scalarField => extract_scalar_field(state, "Absorption1", stat)
    if(stat == 0) then
       call set(scalarField, absorption)  
    end if

end subroutine keps_tke

!----------------------------------------------------------------------------------

subroutine keps_eps(state)

    type(state_type), intent(inout)  :: state
    type(scalar_field), pointer      :: source, absorption, kk, eps, kksource
    type(scalar_field), pointer      :: scalarField, scalar_surface
    type(vector_field), pointer      :: positions
    type(tensor_field), pointer      :: eps_diff, background_diff
    real                             :: prod, diss, EpsOverTke
    integer                          :: i, stat, sele, NNodes_sur
    integer, dimension(:), pointer   :: surface_elements, surface_node_list
    character(len=FIELD_NAME_LEN)    :: bc_type

    ewrite(1,*) "In keps_eps"

    positions  => extract_vector_field(state, "Coordinates")
    kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
    kksource   => extract_scalar_field(state, "TurbulentKineticEnergySource")
    eps        => extract_scalar_field(state, "TurbulentDissipation")
    source     => extract_scalar_field(state, "TurbulentDissipationSource")
    absorption => extract_scalar_field(state, "TurbulentDissipationAbsorption")
    eps_diff   => extract_tensor_field(state, "TurbulentDissipationDiffusivity")

    do i=1,nNodes
        ! clip at k_min: it may be < 0 after solve
        call set(kk, i, max(kk%val(i), k_min))
        ! re-construct eps at "old" timestep
        call set(eps, i, max((tke_old%val(i)**1.5) / ll%val(i), eps_min) )
        ! compute RHS terms in epsilon equation
        EpsOverTke = eps%val(i) / tke_old%val(i)
        prod       = c_eps_1 * EpsOverTke * kksource%val(i)   ! kk source term couples equations
        diss       = c_eps_2 * EpsOverTke
        call set(source, i, prod / eps%val(i) )
        call set(absorption, i, diss)
    end do

    ! Set diffusivity for Eps - copied from Jon's update 28 Feb revision 12719
    call zero(eps_diff)
    background_diff => extract_tensor_field(state, "BackgroundDiffusivity")
    call set(eps_diff, eps_diff%dim, eps_diff%dim, ED, scale = 1. / sigma_eps)
    call addto(eps_diff, background_diff)

    ! puts the BC boundary values in surface_eps_values (module level variable):
    if (calculate_bcs) then
        call get_boundary_condition(eps, name='eps_boundary', type=bc_type, &
                                       surface_node_list=surface_node_list, &
                                       surface_element_list=surface_elements)
        ewrite(1,*) "tke bc type: ",bc_type
        NNodes_sur = size(surface_node_list)
        do i = 1, size(surface_elements)
            sele = surface_elements(i)
            call eps_bc(state, kk, eps, positions, sele)
        end do
        ! map these onto the actual BC in eps - JON - NOT REQUIRED IF NOT FIX_SURFACE_VALUES?
        scalar_surface => extract_surface_field(eps, 'eps_boundary', "value")
        call remap_field(surface_eps_values, scalar_surface)
    end if

    ! TurbulentDissipation is now ready for solving (see Fluids.F90)
    ewrite_minmax(source%val(:))
    ewrite_minmax(absorption%val(:))

    ! set source and absorption terms in optional output fields
    scalarField => extract_scalar_field(state, "Source2", stat)
    if(stat == 0) then
        call set(scalarField, source)  
    end if
    scalarField => extract_scalar_field(state, "Absorption2", stat)
    if(stat == 0) then
        call set(scalarField, absorption)  
    end if

end subroutine keps_eps

!----------
! eddyvisc calculates the lengthscale, and then the eddy viscosity.
! Eddy viscosity is placed in the velocity viscosity
! GLS has a "fix surface values" option (see options tree notes)
! that may stabilise some simulations by re-imposing BCs here. I don't use it.
!----------

subroutine keps_eddyvisc(state)

    type(state_type), intent(inout)  :: state
    type(tensor_field), pointer      :: eddy_visc, eddy_diff, viscosity, background_visc, background_diff
    type(scalar_field), pointer      :: kk, eps, scalarField
    integer                          :: i, stat

    ewrite(1,*) "In keps_eddyvisc"

    kk        => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps       => extract_scalar_field(state, "TurbulentDissipation")
    eddy_visc => extract_tensor_field(state, "EddyViscosity")
    eddy_diff => extract_tensor_field(state, "EddyDiffusivity")
    viscosity => extract_tensor_field(state, "Viscosity")

    do i=1,nNodes
        call set(eps, i, max(eps%val(i),eps_min) )   ! clip at eps_min
        ! compute dissipative scale - MAY NEED CHANGING TO INCLUDE LENGTH LIMITER
        call set(ll, i, (kk%val(i))**1.5 / eps%val(i))

    ! calculate wall-normal element mesh size - from weak bcs
    !call compute_inverse_jacobian( ele_val(x, ele), ele_shape(x, ele), invJ )
    !G = matmul(transpose(invJ(:,:,1)), invJ(:,:,1))
    !n(:,1) = normal_bdy(:,1)
    !hb = 1. / sqrt( matmul(matmul(transpose(n), G), n) )
    !h  = hb(1, 1)



        ! calculate viscosity and diffusivity for next step and for use in other fields
        call set( EV, i, C_mu * sqrt(kk%val(i)) * ll%val(i) )       ! momentum
        call set( ED, i, C_mu * sqrt(kk%val(i)) * ll%val(i) )       ! temperature etc.
    end do

    ewrite(1,*) "Set k-epsilon eddy-diffusivity and eddy-viscosity tensors for use in other fields"
    call zero(eddy_visc) ! zero it first as we're using an addto below
    call zero(eddy_diff)
    call set(eddy_visc, eddy_visc%dim, eddy_visc%dim, EV) 
    call set(eddy_diff, eddy_diff%dim, eddy_diff%dim, ED)

    background_visc => extract_tensor_field(state, "BackgroundViscosity")
    call addto(eddy_visc, background_visc)
    background_diff => extract_tensor_field(state, "BackgroundDiffusivity")
    call addto(eddy_diff, background_diff)

    ewrite_minmax(kk)
    ewrite_minmax(eps)
    ewrite_minmax(EV)
    ewrite_minmax(ED)
    ewrite_minmax(ll)

    ewrite(1,*) "Set viscosity field"
    call zero(viscosity)
    call set(viscosity,viscosity%dim,viscosity%dim,EV)
    call addto(viscosity,background_visc)

    ! Set output on optional fields
    scalarField => extract_scalar_field(state, "LengthScale", stat)
    if(stat == 0) then
        call set(scalarField, ll) 
    end if
    scalarField => extract_scalar_field(state, "EddyViscosity", stat)
    if(stat == 0) then
        call set(scalarField, EV)  
    end if

end subroutine keps_eddyvisc

!------------------------------------------------------------------------------------

subroutine keps_cleanup()

    ewrite(1,*) "In keps_cleanup"

    call deallocate(ll)
    call deallocate(EV)
    call deallocate(ED)
    call deallocate(tke_old)
    if (calculate_bcs) then
        call deallocate(surface_kk_values)
        call deallocate(surface_eps_values)
    end if

end subroutine keps_cleanup

!---------
! Needs to be called after an adapt to reset the fields
! and arrays within the module
!----------

subroutine keps_adapt_mesh(state)

    type(state_type), intent(inout) :: state

    ewrite(1,*) "In keps_adapt_mesh"
    call keps_allocate_fields(state)   ! reallocate everything
    if (calculate_bcs) then
        call keps_init_surfaces(state) ! re-do the boundaries
    end if
    ! We need to repopulate the fields internal to this module, post adapt.
    call keps_tke(state)
    call keps_eps(state)
    call keps_eddyvisc(state)

end subroutine keps_adapt_mesh

!---------------------------------------------------------------------------------

subroutine keps_check_options

    ewrite(1,*) "In keps_check_options"

    ! Don't do k-epsilon if it's not included in the model!
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/")) return

    ! checking for required fields
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentKineticEnergy")) then
        FLExit("You need TurbulentKineticEnergy field for k-epsilon")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentDissipation")) then
        FLExit("You need TurbulentDissipation field for k-epsilon")
    end if
    ! check that the diffusivity is on for the two turbulent fields, and are they diagnostic?
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &tensor_field::Diffusivity")) then
        FLExit("You need TurbulentKineticEnergy Diffusivity field for k-epsilon")
    end if    
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &tensor_field::Diffusivity/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentKineticEnergy Diffusivity field set to diagnostic/internal")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &tensor_field::Diffusivity")) then
        FLExit("You need TurbulentDissipation Diffusivity field for k-epsilon")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &tensor_field::Diffusivity/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentDissipation Diffusivity field set to diagnostic/internal")
    end if
    ! source terms...similarly, not sure I need these set in my flml?
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &tensor_field::Source")) then
        FLExit("You need TurbulentKineticEnergy Source field for k-epsilon")
    end if    
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &tensor_field::Source/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentKineticEnergy Source field set to diagnostic/internal")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &tensor_field::Source")) then
        FLExit("You need TurbulentDissipation Source field for k-epsilon")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &tensor_field::Source/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentDissipation Source field set to diagnostic/internal")
    end if
    ! absorption terms
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &tensor_field::Absorption")) then
        FLExit("You need TurbulentKineticEnergy Absorption field for k-epsilon")
    end if    
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &tensor_field::Absorption/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentKineticEnergy Absorption field set to diagnostic/internal")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &tensor_field::Absorption")) then
        FLExit("You need TurbulentDissipation Absorption field for k-epsilon")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &tensor_field::Absorption/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentDissipation Absorption field set to diagnostic/internal")
    end if
    ! background diffusivities also needed?
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &tensor_field::BackgroundDiffusivity/prescribed")) then
        FLExit("You need BackgroundDiffusivity tensor field for k-epsilon")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &tensor_field::BackgroundViscosity/prescribed")) then
        FLExit("You need BackgroundViscosity tensor field for k-epsilon")
    end if
    ! check for velocity
    if (.not.have_option("/material_phase[0]/vector_field::Velocity")) then
        FLExit("You need Velocity field for k-epsilon")
    end if
    ! these fields allow the new diffusivities/viscosities to be used in other calculations
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &tensor_field::EddyViscosity")) then
        FLExit("You need EddyViscosity field for k-epsilon")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &tensor_field::EddyDiffusivity")) then
        FLExit("You need EddyDiffusivity field for k-epsilon")
    end if
    ! check there's a viscosity somewhere
    if (.not.have_option("/material_phase[0]/vector_field::Velocity/prognostic/&
                          &tensor_field::Viscosity/")) then
        FLExit("Need viscosity switched on under the Velocity field for k-epsilon.") 
    end if
    ! check that the user has switched Velocity/viscosity to diagnostic
    if (.not.have_option("/material_phase[0]/vector_field::Velocity/prognostic/&
                          &tensor_field::Viscosity/diagnostic/")) then
        FLExit("You need to switch the viscosity field under Velocity to diagnostic/internal")
    end if

  end subroutine keps_check_options


!------------------------------------------------------------------!
!                       Private subroutines                        !
!------------------------------------------------------------------!

!----------
! initialise the surface meshes used for the BCS. Called at startup.
!----------
subroutine keps_init_surfaces(state)

    type(state_type), intent(in)     :: state
    type(scalar_field), pointer      :: kk, eps
    type(mesh_type), pointer         :: surface_mesh

    ewrite(1,*) "In keps_init_surfaces"

    ! grab hold of some essential fields
    kk  => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps => extract_scalar_field(state, "TurbulentDissipation")

    ! if we're already initialised, then deallocate surface fields to make space for new ones
    ! Not needed now??? See latest GLS
    if (initialised) then
        ewrite(1,*) "initialised: ", initialised
        call deallocate(surface_kk_values)
        call deallocate(surface_eps_values)
        call deallocate(surface_mesh)
    end if

    ! ONLY ONE MESH/ONE CONDITION FOR NOW. Must have the same surface ids in flml!
    ! If domain is open, we also need inflow BCs: be careful applying them!
    ! N.B. JON DOES THIS PURELY FOR THE FIX_SURFACE_VALUES OPTION:
    ! THIS IS AN EXTRA MESH USED TO IMPOSE LENGTHSCALES ETC. ON THE ACTUAL
    ! BOUNDARY IF THE OPTION IS CHOSEN.
    ! BCS_FROM_OPTIONS ALREADY CREATES THE ACTUAL SURFACE MESH,
    ! WHICH SHOULD BE ALL I NEED.

    ewrite(1,*) "allocate space for surface kk and eps values"
    call get_boundary_condition(kk, 'tke_boundary', surface_mesh=surface_mesh)
    call allocate(surface_kk_values, surface_mesh, name="SurfaceValuesTKE")
    call get_boundary_condition(eps, 'eps_boundary', surface_mesh=surface_mesh)
    call allocate(surface_eps_values, surface_mesh, name="SurfaceValuesEps")

    if(node_count(surface_kk_values) /= node_count(surface_eps_values)) then
        FLAbort("kk and eps BCs are on different surface meshes")
    end if

    ewrite(1,*) "Leaving keps_init_surfaces"
    initialised = .true.

end subroutine keps_init_surfaces

!----------
! tke_bc calculates the Dirichlet BCs on the TKE (kk) field.
! BC is applied to the surface ids on which no-slip velocity bcs exist.
!----------

subroutine tke_bc(state)

    type(state_type), intent(in)     :: state
    type(scalar_field), pointer      :: kk
    character(len=FIELD_NAME_LEN)    :: bc_type
    integer                          :: i, NNodes_sur
    integer, dimension(:), pointer   :: surface_node_list


    kk => extract_scalar_field(state, "TurbulentKineticEnergy")
    ewrite(1,*) "tke bc count: ", get_boundary_condition_count(kk)

        call get_boundary_condition(kk, 'tke_boundary', type=bc_type, surface_node_list=surface_node_list)
        ewrite(1,*) "tke bc type: ", bc_type
        NNodes_sur = size(surface_node_list)
        select case(bc_type)
        case("Dirichlet")
            do i=1,NNodes_sur
                call set(surface_kk_values, i, k_min)   ! not zero, that causes NaNs
                call set(kk, surface_node_list(i), surface_kk_values%val(i))
            end do
        case default
            FLAbort('Unknown surface BC for TKE')
        end select

end subroutine tke_bc

!----------
! eps_bc calculates the boundary conditions on the eps field = 2*EV*(d(k^0.5)/dn)^2.
!----------

subroutine eps_bc(state, kk, eps, positions, sele)

    type(state_type), intent(inout)                           :: state
    type(scalar_field), pointer, intent(inout)                :: kk, eps
    type(vector_field), intent(in)                            :: positions
    type(element_type)                                        :: shape_kk, shape_bdy
    integer, intent(in)                                       :: sele
    integer                                                   :: i, j, ele, NNodes_sur
    integer, dimension(:), pointer                            :: surface_elements, surface_node_list
    real                                                      :: eps_bdy
    real, dimension(face_ngi(positions, sele))                :: detwei_bdy
    real, dimension(positions%dim, face_ngi(positions, sele)) :: normal_bdy
    real, dimension(ele_loc(positions, sele), &
                    ele_ngi(positions, sele), positions%dim)  :: dshape_kk
    real, dimension(face_loc(positions, sele))                :: detwei
    real, dimension(face_ngi(positions, sele))                :: rhs
    real, dimension(ele_loc(positions, sele))                 :: dkdn

    ewrite(1,*) "In eps_bc"

        !do i=1,NNodes_sur
            !call set(surface_eps_values, i, eps_min)   ! not zero, that causes NaNs

    ! 05/03/10: Chris Pain suggests a good idea.
    ! Instead of evaluating eps = f(d(k^.5)/dn) at the nodes, do it at the Gauss quadrature points on the surface.
    ! See notes. dk/dn can be written as integral dotted with normal.
    ! Then set integral of (surface_shape_fns *(eps - f(d(k^.5)/dn) ) ) = 0 over the wall-adjacent element surface.
    ! Loop over surface elements to construct a surface mass matrix (not currently used in Fluidity).
    ! Lump the mass matrix to diagonalise it, then simply solve the system for epsilon.
    ! If/when nodal values of epsilon on surface are required, they can be easily extracted.

    ! Form RHS
    do i = 1, size(surface_elements)
        !sele = surface_elements(i)
        ele  = face_ele(positions, sele)   ! element with facet sele
        shape_kk = ele_shape(positions, i)
        ! Get dshape_kk and element quadrature weights: ngi
        call transform_to_physical( positions, ele, &
             shape_kk, dshape=dshape_kk, detwei=detwei )
        ! Get boundary normal and transformed element quadrature weights over surface.
        call transform_facet_to_physical( positions, sele, detwei_f=detwei_bdy, normal=normal_bdy )
        ! dot normal with dshape_kk to get dk/dn: dim * ngi
        dkdn = dshape_dot_vector_rhs( dshape_kk, vector=normal_bdy, detwei=detwei )
        ! Get surface shape functions
        shape_bdy = face_shape( positions, sele )
        ! Construct RHS surface integral
        rhs = shape_rhs( shape_bdy, detwei_bdy ) * 2.0 * kk%val(i) * EV%val(i) * dkdn**2.0
    end do

!.....................................

            ! copy the boundary values onto the mesh using the global node id
            !call set(eps, surface_node_list(i), surface_eps_values%val(i))
        !end do
    !case default
        !FLAbort('Unknown surface BC for TurbulentDissipation')
    !end select

end subroutine eps_bc

!------------------------------------------------------------------------------------------

subroutine gradient_surface_ele(infield, positions, ele, sele, gradient)
    ! Differentiates a field on surfaces of an element.

    type(scalar_field), intent(in)                              :: infield
    type(vector_field), intent(in)                              :: positions
    type(element_type)                                          :: augmented_shape
    type(element_type), pointer                                 :: u_shape, u_f_shape, x_shape
    integer, intent(in)                                         :: sele
    integer                                                     :: i, ele, l_face_number
    logical                                                     :: compute
    real, dimension(positions%dim, &
          face_loc(positions, sele)),         intent(inout)     :: gradient
    real, dimension(face_ngi(positions, sele))                  :: detwei
    real, dimension(positions%dim, &
          positions%dim, ele_ngi(infield, sele))                :: invJ
    real, dimension(positions%dim, &
          positions%dim, face_ngi(infield, sele))               :: invJ_face
    real, dimension(mesh_dim(positions), &
          face_loc(positions, sele), face_loc(positions, sele)) :: tensor1
    real, dimension(face_loc(positions, sele))                  :: tensor2
    real, dimension(face_loc(positions, sele), &
          face_ngi(positions, sele), mesh_dim(positions))       :: vol_dshape_face

    ewrite(1,*) "In gradient_surface_ele: calculated k-epsilon BCs"

    x_shape   => ele_shape(positions, ele)
    u_shape   => ele_shape(infield, ele)
    u_f_shape => face_shape(infield, sele)

    ! don't compute if the field is constant
    compute= (maxval(infield%val) /= minval(infield%val))

    ewrite(1,*) "Compute inverse Jacobian for k-epsilon boundary condition"

    call compute_inverse_jacobian( ele_val(positions, ele), ele_shape(positions, ele), invJ )
    invJ_face = spread(invJ(:, :, 1), 3, size(invJ_face, 3))
    l_face_number = local_face_number(infield, sele)
    augmented_shape = make_element_shape(x_shape%loc, u_shape%dim, &
         u_shape%degree, u_shape%quadrature, quad_s=u_f_shape%quadrature )
    vol_dshape_face = eval_volume_dshape_at_face_quad( augmented_shape, l_face_number, invJ_face )

    ewrite(1,*) "Compute detwei for k-epsilon boundary condition"

    ! Compute detwei
    call transform_facet_to_physical(positions, sele, detwei_f=detwei)
    tensor1 = shape_dshape(face_shape(positions, sele), vol_dshape_face, detwei)
    tensor2 = face_val(infield, sele)
    do i=1,size(tensor1)
        ewrite(1,*) "tensor1 size:",(size(tensor1,i))
        ewrite(1,*) "tensor2 size:",(size(tensor2,i))
    end do
    if (compute) then
        gradient = tensormul(tensor1, tensor2)
    else
        gradient = 0
    end if

    do i=1,size(gradient)
        ewrite(1,*) "gradient tensor size:",(size(gradient,i))
    end do

end subroutine gradient_surface_ele

!------------------------------------------------------------------------------------------

subroutine keps_allocate_fields(state)

    type(state_type), intent(inout) :: state
    type(vector_field), pointer     :: vectorField
    ewrite(1,*) "In keps_allocate_fields"

    vectorField => extract_vector_field(state,"Velocity")
    nNodes = node_count(vectorField)

    ! allocate some space for the fields we need for calculations, but are optional in the model
    call allocate(ll,      vectorField%mesh, "LengthScale")
    call allocate(EV,      vectorField%mesh, "EddyVisc")
    call allocate(ED,      vectorField%mesh, "EddyDiff")
    call allocate(tke_old, vectorField%mesh, "Old_TKE")

end subroutine keps_allocate_fields

end module k_epsilon

