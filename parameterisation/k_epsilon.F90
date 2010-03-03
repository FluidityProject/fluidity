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
  use boundary_conditions
  use fields_manipulation
  use FLDebug

  implicit none

  private

  ! These variables are the parameters requried by k-epsilon. 
  ! They are all private to prevent tampering
  ! and saved so that we don't need to call keps_init every time we need them.

  ! old turbulent kinetic energy, lengthscale, eddy viscosity, eddy diffusivity, production
  type(scalar_field), save :: tke_old, ll, EV, ED, P
  type(tensor_field), save :: DU_DX ! Gradient tensor. does this need to be saved? Can I initialise it later?
  type(scalar_field), save :: surface_values, surface_kk_values
  ! Minimum values of 2 fields to initialise. And empirical constants from Diamond
  real, save               :: eps_min, k_min, c_mu, c_eps_1, c_eps_2, sigma_eps, sigma_k
  integer, save                        :: nNodes, NNodes_sur
  integer, dimension(:), pointer, save :: surface_nodes, surface_node_list, surface_element_list
  logical, save                        :: initialised, calculate_bcs

  ! The following are the public subroutines
  public :: keps_init, keps_cleanup, keps_tke, keps_eps, keps_eddyvisc, keps_adapt_mesh, keps_check_options

  ! General plan is:
  !  - Init in populate_state
  !  - If solve is about to do TKE, call calc_tke (which calculates P and sets source/absorption for solve)
  !  - If solve is about to do epsilon, call calc_eps (which fixes TKE surfaces, set source/absorption for solve)
  !  - After calc_eps solve, recalculate the viscosity and the lengthscale
  !  - When done, clean-up
  !
  ! TurbulentKineticEnergy and TurbulentDiffusivity need to have higher priority then other prognostic 
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
    type(vector_field), pointer     :: vectorField
    type(tensor_field), pointer     :: tensorField
    integer                         :: stat

    ewrite(1,*)'Now in k_epsilon turbulence model - JB'

    ! Do we have BCs specified for both fields? Better have!
    calculate_bcs = &
    have_option("/material_phase[0]/subgridscale_parameterisations/k_epsilon/scalar_field::TurbulentKineticEnergy/prognostic/boundary_conditions[0]") &
    .and. have_option("/material_phase[0]/subgridscale_parameterisations/k_epsilon/scalar_field::TurbulentDiffusivity/prognostic/boundary_conditions[0]")

    ! Allocate the temporary, module-level variables
    call keps_allocate_fields(state)

    ! populate the 5 model constants
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/C_mu', C_mu, default = 0.09)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/C_eps_1', c_eps_1, default = 1.44)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/C_eps_2', c_eps_2, default = 1.92)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/sigma_k', sigma_k, default = 1.0)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/sigma_eps', sigma_eps, default = 1.3)

    ! Hard-coded by Jon Hill. Must be something small and nonzero.
    k_min   = 7.6e-6
    eps_min = 1.e-12

    ewrite(1,*) "Parameters"
    ewrite(1,*) "--------------------------------------------"
    ewrite(1,*) "c_mu: ",     c_mu
    ewrite(1,*) "c_eps_1: ",  c_eps_1
    ewrite(1,*) "c_eps_2: ",  c_eps_2
    ewrite(1,*) "sigma_k: ",  sigma_k
    ewrite(1,*) "sigma_eps: ",sigma_eps
    ewrite(1,*) "--------------------------------------------"
    
    ! initialise 2 fields (k, epsilon) with minimum values
    scalarField => extract_scalar_field(state, "TurbulentKineticEnergy")
    call set(scalarField,k_min)
    scalarField => extract_scalar_field(state, "TurbulentDiffusivity")
    call set(scalarField,eps_min)

    ! initialise other fields
    call set(P, 0.0)
    scalarField => extract_scalar_field(state, "EddyViscosity")
    call set(scalarField,1e-6)
    scalarField => extract_scalar_field(state, "LengthScale")
    call set(scalarField,k_min**1.5/eps_min)
    call set(ll,k_min**1.5/eps_min)

    ! initialise surface
    if (calculate_bcs) then
        initialised = .false.
        call keps_init_surfaces(state)
    end if
end subroutine keps_init

!------------------------------------------------------------------------------

subroutine keps_tke(state)

    type(state_type), intent(inout)   :: state
    type(scalar_field), pointer       :: source, absorption, kk, eps, utemp
    type(scalar_field), pointer       :: scalar_surface, scalarField
    type(vector_field), pointer       :: positions
    type(tensor_field), pointer       :: kk_diff, tensorField
    real                              :: diss
    integer                           :: i, ii, j, stat
    character(len=FIELD_NAME_LEN)     :: bc_type
    type(scalar_field), dimension(1)  :: grtemp
    type(vector_field)                :: velocity
    logical, allocatable, dimension(:):: derivs

    positions  => extract_vector_field(state, "Coordinate")
    source     => extract_scalar_field(state, "TurbulentKineticEnergySource")
    absorption => extract_scalar_field(state, "TurbulentKineticEnergyAbsorption")
    kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
    kk_diff    => extract_tensor_field(state, "TurbulentKineticEnergyDiffusivity")
    eps        => extract_scalar_field(state, "TurbulentDiffusivity")
    ! Temporary field for calculating gradient field:
    utemp      => extract_scalar_field(state, "Velocity")
    
    ! now allocate our temp fields
    call allocate(DU_DX, tensorField%mesh, "DUDX")
    allocate( derivs(velocity%dim))
    derivs = .false.

    do i = 1, DU_DX%dim
        do j = 1, DU_DX%dim
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
                ! Production term is sum of components of tensor
                call addto(P, ii, node_val(EV, ii) * DU_DX%val(i,j,ii)+DU_DX%val(j,i,ii) &
                * DU_DX%val(i,j,ii)  )
            end do
        end do
        diss = node_val(eps, ii)
        call set(source, ii, P%val(ii) )
        call set(absorption, ii, diss / node_val(kk, ii))
    end do

    deallocate(derivs)

    ! set diffusivity for kk.
    call set(kk_diff,kk_diff%dim,kk_diff%dim,EV,scale=1./sigma_k)
    ! add in background (need to grab from Diamond, rather than hardcode)
    do ii=1,nNodes
        call addto(kk_diff,kk_diff%dim,kk_diff%dim,ii,1.0e-6)
    end do

    ! boundary conditions
    if (calculate_bcs) then
        call get_option("/material_phase[0]/subgridscale_parameterisations/k_epsilon/scalar_field::TurbulentKineticEnergy/prognostic/boundary_conditions[0]/type", bc_type)
        ! puts the BC boundary values in surface_values (module level variable):
        do i=1,NNodes_sur 
            call tke_bc(state, bc_type)
        end do
        ! map these onto the actual BC in kk
        scalar_surface => extract_surface_field(kk, 'tke_boundary', "value")
        call remap_field(surface_values, scalar_surface)
    end if
    
    ! finally, we need a copy of the old TKE for eps, so grab it before we solve
    call set(tke_old, kk)

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

end subroutine keps_tke


!----------------------------------------------------------------------------------

subroutine keps_eps(state)

    type(state_type), intent(inout)  :: state
    type(scalar_field), pointer      :: source, absorption, kk, eps, scalarField, scalar_surface
    type(tensor_field), pointer      :: eps_diff
    real                             :: prod, diss, EpsOverTke
    character(len=FIELD_NAME_LEN)    :: bc_type
    integer                          :: i, stat, sele

    ewrite(1,*) "In calc_eps"

    kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps        => extract_scalar_field(state, "TurbulentDiffusivity")
    source     => extract_scalar_field(state, "TurbulentDiffusivitySource")
    absorption => extract_scalar_field(state, "TurbulentDiffusivityAbsorption")
    eps_diff   => extract_tensor_field(state, "TurbulentDiffusivityDiffusivity")

    ! clip at k_min
    do i=1,nNodes
        call set(kk,i, max(node_val(kk,i),k_min))
        call set(eps,i,max(node_val(eps,i),eps_min))
    end do

    ! re-construct eps at "old" timestep. Hence the reason we don't need to work out eps?
    do i=1,nNodes
        call set(eps, i, (node_val(tke_old,i))**1.5)
        call scale(eps, 1./ll%val(i))
    end do

    ! compute RHS production terms in epsilon equation. Removed buoyancy (see GLS)
    do i=1,nNodes
        EpsOverTke = node_val(eps,i)/node_val(tke_old,i)
        prod       = c_eps_1*EpsOverTke*node_val(P,i)
        diss       = c_eps_2*EpsOverTke*node_val(eps,i)
        call set(source, i, prod)
        call set(absorption, i, diss/node_val(eps,i))
    end do

    ! Set eddy viscosity scaled by constant (sigma_eps) for TurbulentDiffusivity
    call set(eps_diff,eps_diff%dim,eps_diff%dim,EV,scale=1./sigma_eps)
    ! Add in background (need to grab from Diamond rather than hardcode)
    do i=1,nNodes
        call addto(eps_diff,eps_diff%dim,eps_diff%dim,i,1.0e-6)
    end do

    ! boundary conditions
    if (calculate_bcs) then
        call get_option("/material_phase[0]/subgridscale_parameterisations/k_epsilon/scalar_field:: &
                          & TurbulentKineticEnergy/prognostic/boundary_conditions[0]/type",bc_type)
        ! puts the BC boundary values in surface_values (module level variable):
        call eps_bc(state, sele, bc_type)
        ! map these onto the actual BC in eps
        scalar_surface => extract_surface_field(eps, 'eps_boundary', "value")
        call remap_field(surface_values, scalar_surface)
    end if

    ! TurbulentDiffusivity is now ready for solving (see Fluids.F90)
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

end subroutine keps_eps

!----------
! eddyvisc calculates the lengthscale, and then the viscosity.
! Viscosity is placed in the velocity viscosity
!----------

subroutine keps_eddyvisc(state)

    type(state_type), intent(inout)  :: state
    type(scalar_field), pointer      :: kk_state, eps_state
    type(scalar_field)               :: kk, eps, kk_copy
    type(tensor_field), pointer      :: eddy_visc, eddy_diff, viscosity, background_visc, background_diff
    !real                             :: x
    character(len=FIELD_NAME_LEN)    :: bc_type
    integer                          :: i, stat, sele
    real                             :: tke

    kk_state  => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps_state => extract_scalar_field(state, "TurbulentDiffusivity")
    eddy_visc => extract_tensor_field(state, "EddyViscosity",stat)
    eddy_diff => extract_tensor_field(state, "EddyDiffusivity",stat)
    viscosity => extract_tensor_field(state, "Viscosity",stat)

    call allocate(kk,kk_state%mesh,"TKE")
    call allocate(eps,eps_state%mesh,"PSI")
    call allocate(kk_copy,kk_state%mesh,"TKE_COPY")
    call set(kk,kk_state)
    call set(eps,eps_state)
    call get_boundary_condition(eps, 'eps_boundary', type=bc_type)

    ! We need to add a dirichlet BC to the TKE first, so we get the right value for the rest
    ! of this calculation on the surfaces
    do i=1,NNodes_sur
        call tke_bc(state, 'dirichlet')
        ! copy the boundary values onto the mesh using the global node id
        call set(kk,surface_nodes(i),node_val(surface_values,i))
    end do

    ! From latest version of GLS: ask John exactly why we do this.
    call set(kk_copy,kk_state)
    call set(kk_state,kk)
    call eps_bc(state, sele, 'dirichlet')
    call set(kk_state,kk_copy)

write(*,*) "In eddy viscosity" 
    do i=1,nNodes
        tke = node_val(kk,i)
        call set(eps, i, tke**1.5/node_val(ll,i))
        ! clip at eps_min
        call set(eps, i, max(node_val(eps,i),eps_min))
        ! compute dissipative scale - MAY NEED CHANGING
        call set(ll, i, tke**1.5/node_val(eps,i))
    end do

    ! calculate viscosity and diffusivity for next step and for use in other fields
    do i=1,nNodes
        ! GLS formula for eddy viscosity:
        !x = sqrt(node_val(kk,i))*node_val(ll,i)
        !call set( EV, i, x)

        ! Classic? (Ferziger & Peric) formula for eddy viscosity:
        call set( EV, i, C_mu * tke**2. / eps%val(i) ) ! momentum
        ED = EV                                        ! tracer
    end do

    ewrite_minmax(ll)
    ewrite_minmax(kk)
    ewrite_minmax(eps)

    !set the eddy_diffusivity and viscosity tensors for use by other fields
    call zero(eddy_visc) ! zero it first as we're using an addto below
    call zero(eddy_diff)
    call set(eddy_visc,eddy_visc%dim,eddy_visc%dim,EV) 
    call set(eddy_diff,eddy_diff%dim,eddy_diff%dim,ED)

    background_visc => extract_tensor_field(state, "BackgroundViscosity",stat)
    if(stat == 0) then
        call addto(eddy_visc,background_visc)
    endif
    background_diff => extract_tensor_field(state, "BackgroundDiffusivity",stat)
    if(stat == 0) then
        call addto(eddy_diff,background_diff)
    endif

    ! Add background viscosity
    if (stat == 0) then   ! we have a background viscosity - see above
        call zero(viscosity)
        call set(viscosity,viscosity%dim,viscosity%dim,EV)
        call addto(viscosity,background_visc) 
    else                  ! we don't - hard code a reasonable value
        call set(viscosity,viscosity%dim,viscosity%dim,EV)
        do i = 1,nNodes
            call addto(viscosity,viscosity%dim,viscosity%dim,i,1.0e-6) 
        end do
    end if

    ! Set output on optional fields - if the field exists, stick something in it.
    ! We only need to do this to those fields that we have allocated ourselves.
    call keps_output_fields(state)

end subroutine keps_eddyvisc

!------------------------------------------------------------------------------------

subroutine keps_cleanup()

    call deallocate(ll)
    call deallocate(P)
    call deallocate(EV)
    call deallocate(DU_DX)
    call deallocate(tke_old)
    if (calculate_bcs) then
        call deallocate(surface_values)
        call deallocate(surface_kk_values)
    end if

end subroutine keps_cleanup

!---------
! Needs to be called after an adapt to reset the fields
! and arrays within the module
!----------

subroutine keps_adapt_mesh(state)

    type(state_type), intent(inout) :: state

    ewrite(1,*) "In keps_adapt_mesh"
    call keps_allocate_fields(state) ! reallocate everything
    if (calculate_bcs) then
        call keps_init_surfaces(state) ! re-do the boundaries
    end if
    ! We need to repopulate the fields internal to this module, post adapt.
    call keps_tke(state)
    call keps_eps(state)
    call keps_eddyvisc(state)

end subroutine keps_adapt_mesh



subroutine keps_check_options
    
    character(len=FIELD_NAME_LEN) :: buffer
    integer                       :: stat
    real                          :: min_tke

    ! Don't do k-epsilon if it's not included in the model!
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/")) return

    ! checking for required fields
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentKineticEnergy")) then
        FLExit("You need TurbulentKineticEnergy field for k-epsilon")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentDiffusivity")) then
        FLExit("You need TurbulentDiffusivity field for k-epsilon")
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
                          &scalar_field::TurbulentDiffusivity/prognostic/&
                          &tensor_field::Diffusivity")) then
        FLExit("You need TurbulentDiffusivity Diffusivity field for k-epsilon")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentDiffusivity/prognostic/&
                          &tensor_field::Diffusivity/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentDiffusivity Diffusivity field set to diagnostic/internal")
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
                          &scalar_field::TurbulentDiffusivity/prognostic/&
                          &tensor_field::Source")) then
        FLExit("You need TurbulentDiffusivity Source field for k-epsilon")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentDiffusivity/prognostic/&
                          &tensor_field::Source/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentDiffusivity Source field set to diagnostic/internal")
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
                          &scalar_field::TurbulentDiffusivity/prognostic/&
                          &tensor_field::Absorption")) then
        FLExit("You need TurbulentDiffusivity Absorption field for k-epsilon")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &scalar_field::TurbulentDiffusivity/prognostic/&
                          &tensor_field::Absorption/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentDiffusivity Absorption field set to diagnostic/internal")
    end if

    ! background diffusivities also needed?
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &tensor_field::BackgroundDiffusivity/prescribed")) then
        FLExit("You need GLSBackgroundDiffusivity tensor field for k-epsilon")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &tensor_field::BackgroundViscosity/prescribed")) then
        FLExit("You need GLSBackgroundViscosity tensor field for k-epsilon")
    end if

    ! check for velocity
    if (.not.have_option("/material_phase[0]/vector_field::Velocity")) then
        FLExit("You need Velocity field for k-epsilon")
    end if

    ! these two fields allow the new diffusivities/viscosities to be used in
    ! other calculations - we need them!
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

    ! check if priorities have been set - if so warn the user this might screw things up
    if (have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                     &scalar_field::TurbulentKineticEnergy/prognostic/priority")) then
        ewrite(-1,*)("WARNING: Priorities for the k-epsilon fields are set internally. Setting them in the FLML might mess things up")
    end if
    if (have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                     &scalar_field::TurbulentDiffusivity/prognostic/priority")) then
        ewrite(-1,*)("WARNING: Priorities for the k-epsilon fields are set internally. Setting them in the FLML might mess things up")
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
    character(len=OPTION_PATH_LEN)   :: input_mesh_name
    type(scalar_field), pointer      :: kk, bc_field
    type(mesh_type)                  :: input_mesh, boundary_mesh

    ! grab hold of some essential fields
    bc_field => extract_scalar_field(state, "BoundaryConditions")
    kk => extract_scalar_field(state, "TurbulentKineticEnergy")

    ! if we're already initialised, then deallocate surface fields to make space for new ones
    if (initialised) then
        call deallocate(surface_values)
    end if

    ! create a surface mesh to place values onto. ONLY ONE MESH/ONE CONDITION FOR NOW.
    call get_option(trim(kk%option_path)//'/prognostic/mesh/name', input_mesh_name)
    input_mesh = extract_mesh(state, input_mesh_name);
    call get_boundary_condition(bc_field, name='BC', surface_element_list=surface_element_list)
    call create_surface_mesh(boundary_mesh, surface_nodes, input_mesh, surface_element_list, 'SurfaceMesh')
    NNodes_sur = node_count(boundary_mesh)
    call allocate(surface_values, boundary_mesh, name="SurfaceValues")
    call allocate(surface_kk_values,boundary_mesh, name="SurfaceValuesTKE")

    initialised = .true.

end subroutine keps_init_surfaces

!----------
! tke_bc calculates the Dirichlet BCs on the TKE (kk) field
!----------

subroutine tke_bc(state, bc_type)

    type(state_type), intent(in)     :: state
    character(len=*), intent(in)     :: bc_type
    integer                          :: i

    select case(bc_type)
    case("dirichlet")
        call set(surface_values,i,0.0)
    case default
        FLAbort('Unknown surface BC for TKE')
    end select

end subroutine tke_bc

!----------
! eps_bc calculates the boundary conditions on the eps field = 2*EV*(d(k^0.5)/dn)^2.
!----------
subroutine eps_bc(state, sele, bc_type)

    type(state_type), intent(inout)           :: state
    integer                                   :: i, j, sele
    type(scalar_field), pointer               :: kk, eps
    type(vector_field), pointer               :: u, positions
    character(len=*), intent(in)              :: bc_type
    real, dimension(1)                        :: dkdn
    real                                      :: eps_bdy

    ! grab hold of some essential fields
    positions => extract_vector_field(state, "Coordinate")
    kk        => extract_scalar_field(state, "TurbulentKineticEnergy")
    u         => extract_vector_field(state, "Velocity")

    do i=1, get_boundary_condition_count(u)
       call get_boundary_condition(u, i, surface_node_list=surface_nodes, &
            surface_element_list=surface_element_list)
    end do

    select case(bc_type)
    case("dirichlet")
        do i=1, NNodes_sur
            kk%val = sqrt(kk%val)
            ! Calculate the k gradient w.r.t. the boundary normal at node i:
            call gradient_surface_ele(kk, positions, sele, dkdn, state )
            dkdn = dkdn**2.
            dkdn = 2. * EV%val * dkdn
            eps_bdy = 0
            do j=1,size(dkdn)
                eps_bdy = eps_bdy + dkdn(j)**2.
            end do
            eps_bdy = sqrt(eps_bdy)
            call set(surface_values, sele, eps_bdy )
            ! copy the boundary values onto the mesh using the global node id
            call set(eps,surface_nodes(i),node_val(surface_values,i))
        end do
    case default
        FLAbort('Unknown surface BC for TurbulentDiffusivity')
    end select

end subroutine eps_bc

!------------------------------------------------------------------------------------------

subroutine gradient_surface_ele(infield, positions, sele, gradient, state)
    ! Differentiates a field on surfaces of an element.

    type(scalar_field), intent(in)                               :: infield
    type(vector_field), intent(in)                               :: positions
    type(state_type), intent(in)                                 :: state
    integer, intent(in)                                          :: sele
    real, dimension(positions%dim, positions%dim), intent(inout) :: gradient

    type(element_type)                                           :: augmented_shape
    type(element_type), pointer                                  :: u_shape, u_f_shape, x_shape
    integer                                                      :: i, ele, loc, l_face_number
    logical                                                      :: compute
    real, dimension(face_ngi(positions, sele))                   :: detwei
    real, dimension(positions%dim, positions%dim, ele_ngi(infield, sele))        :: invJ
    real, dimension(positions%dim, positions%dim, face_ngi(infield, sele))       :: invJ_face
    real, dimension(mesh_dim(positions), face_loc(positions, sele), face_loc(positions, sele)) :: r
    real, dimension(face_loc(positions, sele), face_ngi(positions, sele), mesh_dim(positions)) :: vol_dshape_face

    ele       =  face_ele(positions, sele)
    x_shape   => ele_shape(positions, ele)
    u_shape   => ele_shape(infield, ele)
    u_f_shape => face_shape(infield, sele)

    ! don't compute if the field is constant
    compute= (maxval(infield%val) /= minval(infield%val))

    call compute_inverse_jacobian( ele_val(positions, ele), ele_shape(positions, ele), invJ )
    invJ_face = spread(invJ(:, :, 1), 3, size(invJ_face, 3))
    l_face_number = local_face_number(infield, sele)
    augmented_shape = make_element_shape(x_shape%loc, u_shape%dim, &
         u_shape%degree, u_shape%quadrature, quad_s=u_f_shape%quadrature )
    vol_dshape_face = eval_volume_dshape_at_face_quad( augmented_shape, l_face_number, invJ_face )

    ! Compute detwei
    call transform_facet_to_physical(positions, sele, detwei_f=detwei)
    r = shape_dshape(face_shape(positions, sele), vol_dshape_face, detwei)

    if (compute) then
        gradient = tensormul(r, ele_val(infield, ele), 3)
    end if

end subroutine gradient_surface_ele

!------------------------------------------------------------------------------------------

subroutine keps_allocate_fields(state)

    type(state_type), intent(inout) :: state
    type(vector_field), pointer     :: vectorField

    vectorField => extract_vector_field(state,"Velocity")

    ! allocate some space for the fields we need for calculations, but are optional in the model
    call allocate(ll,      vectorField%mesh, "LengthScale") ! DIMENSION?
    call allocate(P,       vectorField%mesh, "ShearProduction") ! DIMENSION?
    call allocate(EV,     vectorField%mesh, "EddyViscosity") ! DIMENSION?
    call allocate(tke_old, vectorField%mesh, "Old_TKE") ! DIMENSION?

    nNodes = node_count(vectorField)

end subroutine keps_allocate_fields

!-------------------------------------------------------------------------------------------

subroutine keps_output_fields(state)

    type(state_type), intent(in)     :: state
    type(scalar_field), pointer      :: scalarField
    integer                          :: stat

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

