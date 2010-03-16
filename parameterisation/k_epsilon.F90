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
  !use unittest_tools
  use ieee_arithmetic
  use FLDebug

  implicit none

  private

  ! old turbulent kinetic energy, lengthscale, eddy viscosity, eddy diffusivity, production
  type(scalar_field), save :: tke_old, ll, EV, ED, P
  type(scalar_field), save :: surface_kk_values, surface_eps_values
  ! Minimum values of 2 fields to initialise, and empirical constants from flml
  real, save               :: eps_min, k_min, c_mu, c_eps_1, c_eps_2, sigma_eps, sigma_k
  integer, save                        :: nNodes, NNodes_sur
  integer, dimension(:), pointer, save :: surface_node_list, surface_element_list
  logical, save                        :: initialised, calculate_bcs

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
  ! fields such as temperature, velocity and salinity for this to work

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
    call set(P, 0.0)
    call set(EV,1.e-6)
    call set(ED,1.e-6)
    call set(ll, min(k_min**1.5/eps_min, 10.0) ) ! need to set a minimum using some sensible function

    ewrite(1,*) "Internal fields initial values (minmax)"
    ewrite(1,*) "--------------------------------------------"
    ewrite_minmax(P)
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
    type(vector_field), pointer       :: positions, velocity
    type(tensor_field), pointer       :: kk_diff, background_diff, tensorField
    type(tensor_field)                :: DU_DX
    type(scalar_field), dimension(1)  :: grtemp
    type(scalar_field)                :: utemp
    real                              :: diss, duidxj, dujdxi
    integer                           :: i, ii, j, stat
    logical, allocatable, dimension(:):: derivs
    logical                           :: nan

    ewrite(1,*) "In keps_tke"

    positions  => extract_vector_field(state, "Coordinate")
    source     => extract_scalar_field(state, "TurbulentKineticEnergySource")
    absorption => extract_scalar_field(state, "TurbulentKineticEnergyAbsorption")
    kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
    kk_diff    => extract_tensor_field(state, "TurbulentKineticEnergyDiffusivity")
    eps        => extract_scalar_field(state, "TurbulentDissipation")
    ! Temporary fields for calculating gradient:
    velocity   => extract_vector_field(state, "Velocity")
    tensorField=> extract_tensor_field(state, "EddyViscosity")

    ! clip at k_min and eps_min. We can't have negative values noe can we.
    do i=1,nNodes
        call set(kk, i, max(kk%val(i),k_min))
        call set(eps, i, max(eps%val(i),eps_min))
    end do

    ! now allocate our temporary gradient fields
    call allocate(DU_DX, tensorField%mesh, "DUDX")
    if(DU_DX%refcount%count < 0) then
        FLAbort("tensor field refcount is less than zero!")
    end if

    allocate(derivs(velocity%dim))
    derivs = .false.

    call zero(DU_DX)
    do i = 1, velocity%dim
        do j = 1, velocity%dim
            utemp  = extract_scalar_field_from_vector_field(velocity, i)
            grtemp = extract_scalar_field_from_tensor_field(DU_DX, i, j)
            if (i==j) derivs(j) = .true.
            call differentiate_field(utemp, positions, derivs, grtemp )
            derivs = .false.
        end do
    end do

    ewrite(1,*) "In kinetic energy production term loop"
    do ii = 1, nNodes
        do i = 1, velocity%dim
            do j = 1, velocity%dim
                duidxj=DU_DX%val(i,j,ii)
                dujdxi=DU_DX%val(j,i,ii)
                nan=ieee_is_nan(duidxj)
                if(nan) then
                    ewrite(1,*) "duidxj is NaN at node/dim/dim: ", ii, i, j, DU_DX%val(:,:,ii)
                end if
                nan=ieee_is_nan(dujdxi)
                if(nan) then
                    ewrite(1,*) "dujdxi is NaN at node/dim/dim: ", ii, i, j, DU_DX%val(:,:,ii)
                end if
                ! Production term is sum of components of tensor
                call addto(P, ii, EV%val(ii) * DU_DX%val(i,j,ii)+DU_DX%val(j,i,ii) &
                * DU_DX%val(i,j,ii) )
            end do
        end do
        diss = eps%val(ii)    ! dissipation term of k equation is equal to epsilon
        call set(source, ii, P%val(ii) )   ! source term
        call set(absorption, ii, diss / kk%val(ii))    ! absorption term
        if(kk%val(i)==0) then
            ewrite(1,*) "zero at node: ", i
            FLAbort("ERROR: kk field contains a zero!")
        end if
    end do

    call deallocate(DU_DX)
    deallocate(derivs)

    ! set diffusivity for KK - copied from Jon's update 28 Feb revision 12719
    call zero(kk_diff)
    background_diff => extract_tensor_field(state, "BackgroundDiffusivity")
    call set(kk_diff,kk_diff%dim,kk_diff%dim,EV,scale=1./sigma_k)
    call addto(KK_diff,background_diff)

    ! puts the BC boundary values in surface_kk_values (module level variable):
    if (calculate_bcs) then
        call tke_bc(state)
        ! map these onto the actual BCs in kk - like GLS, but is it necessary if no fix_surface_values?
        scalar_surface => extract_surface_field(KK, 'tke_boundary', "value")
        call remap_field(surface_kk_values, scalar_surface)
    end if
    
    ! finally, we need a copy of the old TKE for eps, so grab it before we solve. Clip at k_min!!!
    do i=1,nNodes
        call set(tke_old, i, max(kk%val(i),k_min))
    end do

    ! that's the TKE set up ready for the solve which is the next thing to happen (see Fluids.F90)
    ewrite_minmax(source%val(:))
    ewrite_minmax(absorption%val(:))

    ! set source and absorption terms in optional output fields
    scalarField => extract_scalar_field(state, "Source1", stat)
    if(stat == 0) then
       call set(scalarField,source)  
    end if
    scalarField => extract_scalar_field(state, "Absorption1", stat)
    if(stat == 0) then
       call set(scalarField,absorption)  
    end if
    
    do i=1, NNodes
        nan=ieee_is_nan(kk%val(i))
        if(nan) then
            ewrite(1,*) "NaN at node: ", i
            FLAbort("ERROR: kk field contains NaN!")
        end if
        nan=ieee_is_nan(tke_old%val(i))
        if(nan) then
            ewrite(1,*) "NaN at node: ", i
            FLAbort("ERROR: tke_old field contains NaN!")
        end if
        nan=ieee_is_nan(eps%val(i))
        if(nan) then
            ewrite(1,*) "NaN at node: ", i
            FLAbort("ERROR: eps field contains NaN!")
        end if
        nan=ieee_is_nan(P%val(i))
        if(nan) then
            ewrite(1,*) "NaN at node: ", i
            FLAbort("ERROR: P field contains NaN!")
        end if
        nan=ieee_is_nan(source%val(i))
        if(nan) then
            ewrite(1,*) "NaN at node: ", i
            FLAbort("ERROR: source field contains NaN!")
        end if
        nan=ieee_is_nan(absorption%val(i))
        if(nan) then
            ewrite(1,*) "NaN at node: ", i
            FLAbort("ERROR: absorption field contains NaN!")
        end if
        ! Check for negative fields. These are not allowed.
        if(kk%val(i) < 0.) then
            FLAbort("ERROR: kk field is negative!")
        end if
        if(eps%val(i) < 0.) then
            FLAbort("ERROR: eps field is negative!")
        end if
    end do

    

end subroutine keps_tke

!----------------------------------------------------------------------------------

subroutine keps_eps(state)

    type(state_type), intent(inout)  :: state
    type(scalar_field), pointer      :: source, absorption, kk, eps, scalarField, scalar_surface
    type(tensor_field), pointer      :: eps_diff, background_diff
    real                             :: prod, diss, EpsOverTke
    character(len=FIELD_NAME_LEN)    :: bc_type
    integer                          :: i, stat
    logical                          :: nan

    ewrite(1,*) "In keps_eps"

    kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps        => extract_scalar_field(state, "TurbulentDissipation")
    source     => extract_scalar_field(state, "TurbulentDissipationSource")
    absorption => extract_scalar_field(state, "TurbulentDissipationAbsorption")
    eps_diff   => extract_tensor_field(state, "TurbulentDissipationDiffusivity")

    ! clip at k_min and eps_min
    do i=1,nNodes
        call set(kk, i, max(kk%val(i),k_min))
        call set(eps, i, max(eps%val(i),eps_min))
    end do

    ! re-construct eps at "old" timestep. Hence the reason we don't need to work out eps?
    do i=1,nNodes
        call set(eps, i, max((tke_old%val(i)**1.5)/ll%val(i), eps_min) )
        if(tke_old%val(i)==0) then
            ewrite(1,*) "zero at node: ", i
            FLAbort("ERROR: tke_old field contains a zero!")
        end if
        nan=ieee_is_nan(tke_old%val(i))
        if(nan) then
            ewrite(1,*) "NaN at node: ", i
            FLAbort("ERROR: tke_old field contains NaN!")
        end if
        nan=ieee_is_nan(eps%val(i))
        if(nan) then
            ewrite(1,*) "NaN at node: ", i
            FLAbort("ERROR: Eps field contains NaN!")
        end if
        if(eps%val(i)==0) then
            ewrite(1,*) "zero at node: ", i
            FLAbort("ERROR: eps field contains a zero!")
        end if
        nan=ieee_is_nan(ll%val(i))
        if(nan) then
            ewrite(1,*) "NaN at node: ", i
            FLAbort("ERROR: ll field contains NaN!")
        end if

        ! compute RHS production terms in epsilon equation.
        EpsOverTke = eps%val(i) / tke_old%val(i)
        prod       = c_eps_1 * EpsOverTke * max(P%val(i), 1.e-6)
        diss       = c_eps_2 * EpsOverTke * eps%val(i)
        call set(source, i, prod)
        call set(absorption, i, diss/eps%val(i))
    end do

    ! Set diffusivity for Eps - copied from Jon's update 28 Feb revision 12719
    call zero(eps_diff)
    background_diff => extract_tensor_field(state, "BackgroundDiffusivity")
    call set(eps_diff,eps_diff%dim,eps_diff%dim,ED,scale=1./sigma_eps)
    call addto(eps_diff,background_diff)

    ! puts the BC boundary values in surface_eps_values (module level variable):
    if (calculate_bcs) then
        call eps_bc(state)
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
        call set(scalarField,source)  
    end if
    scalarField => extract_scalar_field(state, "Absorption2", stat)
    if(stat == 0) then
        call set(scalarField,absorption)  
    end if

    ! Check for negative fields. These are not allowed.
    do i=1, NNodes
        if(kk%val(i) < 0.) then
            FLAbort("ERROR: kk field is negative!")
        end if
        if(eps%val(i) < 0.) then
            FLAbort("ERROR: eps field is negative!")
        end if
    end do

end subroutine keps_eps

!----------
! eddyvisc calculates the lengthscale, and then the viscosity.
! Viscosity is placed in the velocity viscosity
!----------

subroutine keps_eddyvisc(state)

    type(state_type), intent(inout)  :: state
    type(tensor_field), pointer      :: eddy_visc, eddy_diff, viscosity, background_visc, background_diff
    type(scalar_field), pointer      :: kk_state, eps_state
    type(scalar_field)               :: kk, eps, kk_copy
    !real                             :: x         ! For GLS-style eddy viscosity
    character(len=FIELD_NAME_LEN)    :: bc_type
    integer                          :: i, dims
    real                             :: tke

    ewrite(1,*) "In keps_eddyvisc"

    kk_state  => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps_state => extract_scalar_field(state, "TurbulentDissipation")
    eddy_visc => extract_tensor_field(state, "EddyViscosity")
    eddy_diff => extract_tensor_field(state, "EddyDiffusivity")
    viscosity => extract_tensor_field(state, "Viscosity")

    ewrite(1,*) "allocating temporary fields"
    call allocate(kk,kk_state%mesh,"TKE")
    call allocate(eps,eps_state%mesh,"PSI")
    call allocate(kk_copy,kk_state%mesh,"TKE_COPY")
    call set(kk,kk_state)
    call set(eps,eps_state)

    ewrite(1,*) "Calling TKE BC"
    ! We need to add a dirichlet BC to the TKE first,
    ! so we get the right value for the rest of this calculation on the surfaces
    ! JON: DO I STILL NEED TO DO THIS IF I DON'T FIX_SURFACE_VALUES?
    call tke_bc(state)

    ! From latest version of GLS: ask John exactly why we do this.
    ewrite(1,*) "setting kk_copy and kk_state"
    call set(kk_copy,kk_state)
    call set(kk_state,kk)
    call eps_bc(state)
    call set(kk_state,kk_copy)

    ewrite(1,*) "In eddy viscosity"

    do i=1,nNodes
        tke = kk%val(i)
        call set(eps, i, tke**1.5/ll%val(i))
        ! clip at eps_min
        call set(eps, i, max(eps%val(i),eps_min))
        ! compute dissipative scale - MAY NEED CHANGING
        call set(ll, i, tke**1.5/eps%val(i))
    end do

    ! calculate viscosity and diffusivity for next step and for use in other fields
    do i=1,nNodes
        ! GLS formula for eddy viscosity:
        !x = sqrt(kk%val(i))*ll%val(i)
        !call set( EV, i, x)

        ! Classic? (Ferziger & Peric) formula for eddy viscosity:
        call set( EV, i, C_mu*tke**2./eps%val(i) )       ! momentum
        call set( ED, i, C_mu*tke**2./eps%val(i) )       ! tracer
    end do

    ewrite(1,*) "Set k-epsilon tensors for use in other fields"

    ! Set the eddy_diffusivity and viscosity tensors for use by other fields
    call zero(eddy_visc) ! zero it first as we're using an addto below
    call zero(eddy_diff)
    call set(eddy_visc,eddy_visc%dim,eddy_visc%dim,EV) 
    call set(eddy_diff,eddy_diff%dim,eddy_diff%dim,ED)

    ewrite(1,*) "Set background viscosity and diffusivity fields in state"

    background_visc => extract_tensor_field(state, "BackgroundViscosity")
    call addto(eddy_visc,background_visc)
    background_diff => extract_tensor_field(state, "BackgroundDiffusivity")
    call addto(eddy_diff,background_diff)

    ewrite_minmax(kk)
    ewrite_minmax(eps)
    ewrite_minmax(EV)
    ewrite_minmax(ED)
    ewrite_minmax(ll)

    do i=1, NNodes
        if(kk%val(i) < 0.) then
            FLAbort("ERROR: kk field is negative!")
        end if
        if(eps%val(i) < 0.) then
            FLAbort("ERROR: eps field is negative!")
        end if
    end do

    call deallocate(kk)
    call deallocate(eps)
    call deallocate(kk_copy)

    ewrite(1,*) "Set viscosity field"

    ! Set viscosity
    call zero(viscosity)
    call set(viscosity,viscosity%dim,viscosity%dim,EV)
    call addto(viscosity,background_visc) 

    ! Set output on optional fields - if the field exists, stick something in it.
    ! We only need to do this to those fields that we have allocated ourselves.
    call keps_output_fields(state)

end subroutine keps_eddyvisc

!------------------------------------------------------------------------------------

subroutine keps_cleanup()

    logical   :: have_surface_mesh

    ewrite(1,*) "In keps_cleanup"

    if(associated(surface_kk_values%mesh%faces)) then  ! surface_mesh is a pointer?
        have_surface_mesh = associated(surface_kk_values%mesh%faces%boundary_ids)
        ewrite(1,*) "surface_kk_values field has the following boundary ids: ", &
                    (surface_kk_values%mesh%faces%boundary_ids)
    else
        have_surface_mesh = .false.
        ewrite(1,*) "surface_kk_values field mesh does not have boundary ids associated with it"
    end if

    call deallocate(ll)
    call deallocate(P)
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
    integer                          :: nbcs

    ewrite(1,*) "In keps_init_surfaces"

    ! grab hold of some essential fields
    kk  => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps => extract_scalar_field(state, "TurbulentDissipation")

    ! if we're already initialised, then deallocate surface fields to make space for new ones
    ! Not needed now! See latest GLS
    if (initialised) then
            ewrite(1,*) "initialised: ", initialised
        call deallocate(surface_kk_values)
        call deallocate(surface_eps_values)
        call deallocate(surface_mesh)
    end if

    ewrite(1,*) "Grab BCs for kk and eps fields"
    ! ONLY ONE MESH/ONE CONDITION FOR NOW. Must have the same surface ids in flml!
    ! N.B. JON DOES THIS PURELY FOR THE FIX_SURFACE_VALUES OPTION:
    ! THIS IS AN EXTRA MESH USED TO IMPOSE LENGTHSCALES ETC. ON THE ACTUAL
    ! BOUNDARY IF THE OPTION IS CHOSEN.
    ! BCS_FROM_OPTIONS ALREADY CREATES THE ACTUAL SURFACE MESH,
    ! WHICH SHOULD BE ALL I NEED.

    do nbcs=1, get_boundary_condition_count(kk)
        call get_boundary_condition(field=kk, n=nbcs, surface_mesh=surface_mesh)
        ! Set module-level counter
        NNodes_sur = node_count(surface_mesh)
    end do

    ewrite(1,*) "allocate space for surface kk values"
    call allocate(surface_kk_values,surface_mesh, name="SurfaceValuesTKE")

    do nbcs=1, get_boundary_condition_count(eps)
        call get_boundary_condition(field=eps, n=nbcs, surface_mesh=surface_mesh)
        if(node_count(surface_mesh) /= NNodes_sur) then
            FLAbort("kk and eps BCs are on different surface meshes")
        end if
    end do

    ewrite(1,*) "allocate space for surface eps values"
    call allocate(surface_eps_values,surface_mesh, name="SurfaceValuesEps")

    ewrite(1,*) "Leaving keps_init_surfaces"
    initialised = .true.

end subroutine keps_init_surfaces

!----------
! tke_bc calculates the Dirichlet BCs on the TKE (kk) field
!----------

subroutine tke_bc(state)

    type(state_type), intent(in)     :: state
    type(scalar_field), pointer      :: kk
    character(len=FIELD_NAME_LEN)    :: bc_type
    integer                          :: i, bcs
    logical                          :: have_surface_mesh

    if(associated(surface_kk_values%mesh%faces)) then  ! surface_mesh is a pointer?
        have_surface_mesh = associated(surface_kk_values%mesh%faces%boundary_ids)
        ewrite(1,*) "surface_kk_values field has the following boundary ids: ", (surface_kk_values%mesh%faces%boundary_ids)
    else
        have_surface_mesh = .false.
        ewrite(1,*) "surface_kk_values field mesh does not have boundary ids associated with it"
    end if

    kk => extract_scalar_field(state, "TurbulentKineticEnergy")
    ewrite(1,*) "tke bc count: ", get_boundary_condition_count(kk)

    do bcs=1, get_boundary_condition_count(kk)
        call get_boundary_condition(field=kk, n=bcs, type=bc_type, surface_node_list=surface_node_list)
        ewrite(1,*) "tke bc type: ",bc_type
        select case(bc_type)
        case("Dirichlet")
            do i=1,NNodes_sur
                call set(surface_kk_values,i,0.0)
                call set(kk,surface_node_list(i),surface_kk_values%val(i))
            end do
        case default
            FLAbort('Unknown surface BC for TKE')
        end select
    end do

end subroutine tke_bc

!----------
! eps_bc calculates the boundary conditions on the eps field = 2*EV*(d(k^0.5)/dn)^2.
!----------

subroutine eps_bc(state)

    type(state_type), intent(inout) :: state
    integer                         :: i, j, ele, sele, bcs
    type(scalar_field), pointer     :: kk, eps
    type(vector_field), pointer     :: positions
    character(len=FIELD_NAME_LEN)   :: bc_type
    real, dimension(1)              :: dkdn
    real                            :: eps_bdy
    logical                         :: have_surface_mesh

    ! 05/03/10: Chris Pain suggests a good idea.
    ! Instead of evaluating eps = f(d(k^.5)/dn) at the nodes, do it at the Gauss quadrature points on the surface.
    ! See notes. dk/dn can be written as integral dotted with normal.
    ! Then set integral of (surface_shape_fns *(eps - f(d(k^.5)/dn) ) ) = 0 over the wall-adjacent element surface.
    ! Loop over surface elements to construct a surface mass matrix (not currently used in Fluidity).
    ! Lump the mass matrix to diagonalise it, then simply solve the system for epsilon.
    ! If/when nodal values of epsilon on surface are required, they can be easily extracted.

    ewrite(1,*) "In eps_bc"

    if(associated(surface_eps_values%mesh%faces)) then  ! surface_mesh is associated with field?
        have_surface_mesh = associated(surface_eps_values%mesh%faces%boundary_ids)
        ewrite(1,*) "surface_eps_values field has the following boundary ids: ", &
                    (surface_eps_values%mesh%faces%boundary_ids)
    else
        have_surface_mesh = .false.
        ewrite(1,*) "surface_eps_values field mesh does not have boundary ids associated with it"
    end if

    ! grab hold of some essential fields
    positions => extract_vector_field(state, "Coordinate")
    kk        => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps       => extract_scalar_field(state, "TurbulentDissipation")
    ewrite(1,*) "eps bc count: ", get_boundary_condition_count(eps)

    ewrite(1,*) "Get boundary condition"
    do bcs=1, get_boundary_condition_count(eps)
        call get_boundary_condition(field=eps, n=bcs, type=bc_type, &
                               surface_node_list=surface_node_list, &
                               surface_element_list=surface_element_list)
        ewrite(1,*) "tke bc type: ",bc_type
        select case(bc_type)
        case("Dirichlet")

            do i=1,NNodes_sur
                call set(surface_eps_values,i,0.0)

!        do i = 1, size(surface_element_list)
!            sele = surface_element_list(i)
!            ele  = face_ele(positions, sele)
!            kk%val = sqrt(kk%val)

            ! Calculate the k gradient w.r.t. the boundary normal at surface element sele:
!            call gradient_surface_ele(kk, positions, ele, sele, dkdn )
!            ewrite(1,*) "dkdn(i):", (i), (dkdn)

!            dkdn = dkdn**2.
!            dkdn = 2. * EV%val * dkdn
!            eps_bdy = 0

!            do j=1,size(dkdn)
!                eps_bdy = eps_bdy + dkdn(j)**2.
!            end do

!            eps_bdy = sqrt(eps_bdy)
!            call set(surface_eps_values, sele, eps_bdy )

            ! copy the boundary values onto the mesh using the global node id
                call set(eps,surface_node_list(i),surface_eps_values%val(i))
            end do

        case default
            FLAbort('Unknown surface BC for TurbulentDissipation')
        end select
    end do

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

    ! allocate some space for the fields we need for calculations, but are optional in the model
    call allocate(ll,      vectorField%mesh, "LengthScale")
    call allocate(P,       vectorField%mesh, "ShearProduction")
    call allocate(EV,      vectorField%mesh, "EddyVisc")
    call allocate(ED,      vectorField%mesh, "EddyDiff")
    call allocate(tke_old, vectorField%mesh, "Old_TKE")

    nNodes = node_count(vectorField)

end subroutine keps_allocate_fields

!-------------------------------------------------------------------------------------------

subroutine keps_output_fields(state)

    type(state_type), intent(in)     :: state
    type(scalar_field), pointer      :: scalarField
    integer                          :: stat

    ewrite(1,*) "In keps_output_fields"

    scalarField => extract_scalar_field(state, "LengthScale", stat)
    if(stat == 0) then
        call set(scalarField,ll) 
    end if
    scalarField => extract_scalar_field(state, "ShearProduction", stat)
    if(stat == 0) then
        call set(scalarField,P) 
    end if
    scalarField => extract_scalar_field(state, "EddyViscosity", stat)
    if(stat == 0) then
       call set(scalarField,EV)  
    end if

end subroutine keps_output_fields

end module k_epsilon

