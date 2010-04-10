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
  use sparse_tools
  use sparse_tools_petsc
  use sparsity_patterns_meshes
  use state_module
  use spud
  use global_parameters, only:   OPTION_PATH_LEN
  use state_fields_module
  use populate_state_module
  use boundary_conditions
  use fields_manipulation
  use ieee_arithmetic
  use fetools
  use vector_tools
  use solvers
  use node_boundary
  use quicksort
  use FLDebug

  implicit none

  private

  ! old turbulent kinetic energy, lengthscale, eddy viscosity, eddy diffusivity
  type(scalar_field), save :: tke_old, ll, EV, ED
  !type(scalar_field), save :: surface_kk_values, surface_eps_values
  ! Minimum values of 2 fields to initialise, and empirical constants from flml
  real, save               :: eps_min, k_min, c_mu, c_eps_1, c_eps_2, sigma_eps, sigma_k
  integer, save            :: nNodes, count_rhs, count_mass, index
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
    k_min   = 1.e-9   ! used to be 7.6e-6
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

end subroutine keps_init

!------------------------------------------------------------------------------

subroutine keps_tke(state)

    type(state_type), intent(inout)   :: state
    type(scalar_field), pointer       :: source, absorption, kk, eps
    type(scalar_field), pointer       :: scalarField
    type(vector_field), pointer       :: velocity, positions
    type(tensor_field), pointer       :: kk_diff, background_diff, visc
    type(tensor_field)                :: DU_DX
    type(scalar_field), dimension(1)  :: grtemp
    type(scalar_field)                :: utemp
    integer                           :: i, ii, j, stat, bc
    logical, allocatable, dimension(:):: derivs
    character(len=FIELD_NAME_LEN)     :: bc_type, bc_name

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

        do bc = 1, get_boundary_condition_count(kk)

            call get_boundary_condition(kk, bc, name=bc_name, type=bc_type)
            ewrite(1,*) "tke bc type: ", bc_type
            ewrite(1,*) "tke bc name: ", bc_name
            ewrite(1,*) "bc count: ", bc

            if (bc_name=='tke_boundary') then

                call tke_bc(state, bc)   ! (almost) zero Dirichlet BC: sets to k_min

            end if

        end do

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
    type(scalar_field), pointer      :: scalarField
    type(vector_field), pointer      :: positions
    type(tensor_field), pointer      :: eps_diff, background_diff
    type(mesh_type), pointer         :: surface_mesh
    type(patch_type)                 :: current_patch
    type(scalar_field)               :: rhs_vector, surface_eps_values, ele_size
    type(petsc_csr_matrix)           :: lumped_mass
    type(csr_sparsity), pointer      :: eps_sparsity

    real                             :: prod, diss, EpsOverTke, &
                                        & h, normal_size_node, ll_limit
    real, dimension(:), pointer      :: x, z
    integer                          :: i, j, bc, stat, ele, sele, nele
    integer, dimension(:), pointer   :: surface_elements, surface_node_list
    character(len=FIELD_NAME_LEN)    :: bc_type, bc_name

    ewrite(1,*) "In keps_eps"

    positions  => extract_vector_field(state, "Coordinate")
    kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
    kksource   => extract_scalar_field(state, "TurbulentKineticEnergySource")
    eps        => extract_scalar_field(state, "TurbulentDissipation")
    source     => extract_scalar_field(state, "TurbulentDissipationSource")
    absorption => extract_scalar_field(state, "TurbulentDissipationAbsorption")
    eps_diff   => extract_tensor_field(state, "TurbulentDissipationDiffusivity")
    x => positions%val(1)%ptr  
    z => positions%val(2)%ptr

    do i = 1, NNodes    ! clip at k_min: it may be < 0 after solve

        call set(kk, i, max(kk%val(i), k_min))    ! re-construct eps at "old" timestep
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

    ! puts the BC boundary values in surface_eps_values:
    if (calculate_bcs) then

        do bc = 1, get_boundary_condition_count(kk)   ! use kk bcs for eps

            call get_boundary_condition(kk, bc, name=bc_name, type=bc_type, &
                 surface_node_list=surface_node_list, surface_mesh=surface_mesh, &
                 surface_element_list=surface_elements)

            ewrite(1,*) "eps bc name: ", bc_name
            ewrite(1,*) "eps bc type: ", bc_type

            if (bc_name /= 'inflow') then   ! nasty way of selecting just the sides. TEMPORARY

                ! create a sparse matrix, rhs and solution vectors
                eps_sparsity => get_csr_sparsity_firstorder( state, surface_mesh, surface_mesh )
                call allocate(lumped_mass, eps_sparsity, (/1,1/), diagonal=.true., name="lumped_mass")
                call allocate(surface_eps_values, surface_mesh, name="SurfaceValuesEps")
                call allocate(rhs_vector, surface_mesh, name="RHS")
                call allocate(ele_size, surface_mesh, name="elementsize")
                call zero(lumped_mass)
                call zero(rhs_vector)
                call zero(surface_eps_values)
                call zero(ele_size)

                ewrite(1,*) "surface field size: ", size( surface_eps_values%val )

                index = 1                             ! initialise array counter

                do i = 1, size(surface_elements)      ! snloc
                    ele  = face_ele(positions, i)     ! index of the element which owns face
                    sele = surface_elements(i)        ! global face id
                    ! Calculate boundary condition at surface element i and add to petsc system
                    call eps_bc(ele, sele, kk, positions, &
                                surface_node_list, rhs_vector, lumped_mass, h)
                    call set(ele_size, i, h)      ! Surface field of wall-normal element sizes
                end do

                ! Solve the linear system
                surface_eps_values%option_path = eps%option_path   ! important: get solver options here.
                call petsc_solve( surface_eps_values, lumped_mass, rhs_vector )

                ! Check for NaNs etc.
                do i = 1, node_count(surface_eps_values)
                    print *, i, surface_node_list(i), x(surface_node_list(i)), &
                             z(surface_node_list(i)), node_val(surface_eps_values,i)
                    if(ieee_is_nan(node_val(rhs_vector,i) ) ) then
                        ewrite(1,*) "NaN in rhs_vector at: ", i
                    end if
                end do

                ! map these onto the actual BC in eps - JON - NOT REQUIRED IF NOT FIX_SURFACE_VALUES?
                !scalar_surface => extract_surface_field(eps, 'eps_boundary', "value")
                !call remap_field( surface_eps_values, scalar_surface )

                do i = 1, node_count(surface_eps_values)

                    ! apply boundary condition to epsilon field
                    call set( eps, surface_node_list(i), surface_eps_values%val(i) )

                    ! for each node in the surface, get the elements and find the dz
                    ! average out for this (i.e. add and divide by nodes per element)
                    current_patch = get_patch_ele(surface_mesh, i)
                    nele = current_patch%count
                    normal_size_node = 0.

                    do j = 1, nele   ! get sizes at adjacent surface elements
                        sele = surface_elements( current_patch%elements(j) )
                        ele = face_ele( positions, sele )
                        
                        normal_size_node = normal_size_node + ele_size%val(j)
                    end do

                    normal_size_node = normal_size_node / nele    ! normalise

                    ! Set lengthscale limit. See Yap (1987) or Craft & Launder (1996)
                    ll_limit = max( ( 0.83 * eps%val(i)**2. / kk%val(i) * &
                    ( kk%val(i)**1.5 / 2.5 / eps%val(i) / normal_size_node - 1. ) * &
                    ( kk%val(i)**1.5 / 2.5 / eps%val(i) / normal_size_node )**2. ), 0. )

                    ! limit the lengthscale at surface
                    ll%val(i) = max( ll_limit, ll%val(i) )
                    ewrite(1,*) "lengthscale limited at node: ", i, surface_node_list(i), ll%val(i)
                    deallocate(current_patch%elements)

                end do

                call deallocate( rhs_vector )
                call deallocate( lumped_mass )
                call deallocate( surface_eps_values )
                call deallocate( ele_size )

                ewrite(1,*) "next surface..."
            end if

        end do

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
    type(vector_field), pointer      :: positions
    type(scalar_field), pointer      :: kk, eps, scalarField
    integer                          :: i, stat
    logical                          :: on_bound

    ewrite(1,*) "In keps_eddyvisc"

    kk        => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps       => extract_scalar_field(state, "TurbulentDissipation")
    positions => extract_vector_field(state, "Coordinate")
    eddy_visc => extract_tensor_field(state, "EddyViscosity")
    eddy_diff => extract_tensor_field(state, "EddyDiffusivity")
    viscosity => extract_tensor_field(state, "Viscosity")

    do i = 1, NNodes

        call set(eps, i, max(eps%val(i),eps_min) )   ! clip epsilon field at eps_min

        on_bound = node_lies_on_boundary( positions%mesh, positions, i )
        if(on_bound .eqv. .true.) then

            ! if lengthscale at surface node > kk**1.5/eps, change it
            if( ll%val(i) > (kk%val(i)**1.5 / eps%val(i)) ) then
                call set(ll, i, (kk%val(i))**1.5 / eps%val(i))
                ewrite(1,*) "lengthscale at surface node is changed: ", i, ll%val(i)
            else
                ewrite(1,*) "lengthscale at surface node is unchanged: ", i, ll%val(i)
            end if

        else   ! set lengthscale at all other nodes to kk**1.5/eps

            call set(ll, i, min( (kk%val(i))**1.5 / eps%val(i), 1.0 ) )
            ewrite(1,*) "lengthscale at volume node is: ", i, ll%val(i)

        end if

        if(ll%val(i) > 2.0 ) then
            ewrite(1,*) "WARNING: lengthscale at node is too big: ", i
        end if

        on_bound = .false.

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

end subroutine keps_cleanup

!---------
! Needs to be called after an adapt to reset the fields
! and arrays within the module
!----------

subroutine keps_adapt_mesh(state)

    type(state_type), intent(inout) :: state

    ewrite(1,*) "In keps_adapt_mesh"
    call keps_allocate_fields(state)   ! reallocate everything
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


subroutine tke_bc(state, bc)

    type(state_type), intent(in)     :: state
    type(scalar_field), pointer      :: kk
    type(scalar_field)               :: surface_kk_values
    type(mesh_type), pointer         :: surface_mesh
    character(len=FIELD_NAME_LEN)    :: bc_type, bc_name
    integer                          :: i
    integer, intent(in)              :: bc
    integer, dimension(:), pointer   :: surface_node_list


    kk => extract_scalar_field(state, "TurbulentKineticEnergy")

    call get_boundary_condition(kk, bc, name=bc_name, type=bc_type, &
         surface_mesh=surface_mesh, surface_node_list=surface_node_list)

    call allocate(surface_kk_values, surface_mesh, name="SurfaceValuesTKE")

    do i = 1, size(surface_node_list)
        call set( surface_kk_values, i, k_min )   ! not zero, that causes NaNs
        call set( kk, surface_node_list(i), surface_kk_values%val(i) )
    end do

    call deallocate(surface_kk_values)

end subroutine tke_bc

!----------
! eps_bc calculates the boundary conditions on the eps field = 2*EV*(d(k^0.5)/dn)^2.
!----------

subroutine eps_bc(ele, sele, kk, positions, surface_node_list, rhs_vector, lumped_mass, h)

    ! 05/03/10: Chris Pain suggests a good idea.
    ! Evaluate eps = f(d(k^.5)/dn) at the Gauss quadrature points on the surface.
    ! See notes. dk/dn can be written as integral dotted with normal.
    ! Then set integral (surface_shape_fns *(eps - f(d(k^.5)/dn) ) ) = 0.
    ! Loop over surface elements to construct a surface mass matrix.
    ! Lump the mass matrix to diagonalise it, then simply solve the system for epsilon.

    type(scalar_field), pointer, intent(inout)    :: kk
    type(vector_field), intent(in)                :: positions
    type(scalar_field), intent(inout)             :: rhs_vector
    type(petsc_csr_matrix), intent(inout)         :: lumped_mass
    integer, intent(in)                           :: ele, sele
    integer, dimension(:), pointer, intent(in)    :: surface_node_list

    type(element_type)                            :: shape_kk, fshape_kk, augmented_shape
    integer                                       :: i, j, snloc, quad, n
    integer, dimension(face_loc(positions, sele)) :: nodes_bdy, new_node_list
    logical                                       :: on_bound
    real, dimension(ele_ngi(kk, ele))             :: detwei
    real, dimension(face_ngi(kk, sele))           :: detwei_bdy, fields_quad, dkdn_contracted
    real, dimension(positions%dim, positions%dim) :: G
    real, dimension(positions%dim, 1)             :: normal
    real, dimension(1,1)                          :: hb
    real                                          :: h     ! wall-normal distance

    real, dimension(positions%dim, positions%dim, ele_ngi(kk, sele))     :: invJ
    real, dimension(positions%dim, positions%dim, face_ngi(kk, sele))    :: invJ_face
    real, dimension(ele_loc(kk, sele), ele_ngi(kk, sele), positions%dim) :: dshape_kk
    real, dimension(ele_loc(kk, ele), face_ngi(kk, sele), positions%dim) :: fdshape_kk
    real, dimension(positions%dim, face_ngi(positions, sele))            :: normal_bdy
    real, dimension(face_loc(positions,sele))                            :: rhs
    real, dimension(face_loc(positions,sele), face_ngi(positions,sele))  :: dkdn
    real, dimension(face_loc(positions,sele), face_loc(positions,sele))  :: mass

    ewrite(1,*) "In eps_bc"

    ! Get ids and shape functions
    quad      = face_ngi(kk, sele)    ! no. of gauss points in surface element
    snloc     = face_loc(kk, sele)    ! no. of nodes on surface element
    shape_kk  = ele_shape(kk, ele)    ! shape functions in element
    fshape_kk = face_shape(kk, sele)  ! shape functions in surface element
    nodes_bdy = face_global_nodes( positions, sele ) ! global node ids in surface element

    do i = 1, size(surface_node_list)
        ! I want the indices in surface_node_list that match nodes_bdy for my addto counter.
        do j = 1, snloc
            if (surface_node_list(i) == nodes_bdy(j)) then
                new_node_list(j) = index
                index = index + 1   ! array counter
            end if
        end do
    end do
    ewrite(1,*) "nodes_bdy: ", nodes_bdy

    if(new_node_list(1) > new_node_list(2)) then   ! inelegant solution to the problem of
        n = size(new_node_list)   ! nodes_bdy having descending indices, so that in above
        new_node_list = new_node_list(n:1:-1)   ! loop 2nd index comes first
        ewrite(1,*) "node list order reversed: ", new_node_list
    else
        ewrite(1,*) "new node list: ", new_node_list
    end if
    index = index - 1

    ! Get dshape_kk and element quadrature weights: ngi
    call transform_to_physical( positions, ele, shape_kk, dshape=dshape_kk, detwei=detwei )

    ! Need inverse Jacobian for augmented shape functions
    call compute_inverse_jacobian( ele_val(positions, ele), shape_kk, invJ )
    invJ_face = spread(invJ(:, :, 1), 3, size(invJ_face, 3))

    ! Transform element shape functions to face
    augmented_shape = make_element_shape(shape_kk%loc, shape_kk%dim, &
                      shape_kk%degree, shape_kk%quadrature, quad_s = fshape_kk%quadrature )

    ! Get dshape on face. nloc x sngi x dim
    fdshape_kk = eval_volume_dshape_at_face_quad( augmented_shape, &
                            local_face_number(kk, sele), invJ_face )

    ! Get boundary normal and transformed element quadrature weights over surface
    call transform_facet_to_physical( positions, sele, detwei_f=detwei_bdy, normal=normal_bdy )

    ! calculate wall-normal element mesh size - from Weak_BCs - NOT NEEDED IF I USE THE FUNCTION BELOW
    G = matmul(transpose(invJ(:,:,1)), invJ(:,:,1))
    normal(:,1) = normal_bdy(:,1)
    hb = 1. / sqrt( matmul(matmul(transpose(normal), G), normal) )
    h  = hb(1, 1)

    do i = 1, size(nodes_bdy)         ! loop over nodes in element: nloc
        ! apparently this doesn't need a global node id
        on_bound = node_lies_on_boundary( kk%mesh, positions, nodes_bdy(i) )
        if(on_bound .eqv. .true.) then                     ! Exclude non-boundary node(s?)
            do j = 1, quad                       ! loop over gauss points in surface: sngi
                ! dot normal with dshape_kk to get shape fn gradient w.r.t. surface normal
                dkdn(i,j) = dot_product( normal_bdy(:,j), fdshape_kk(i,j,:) )
            end do
        end if
        on_bound = .false.
    end do

    do j = 1, quad      ! Sum dshape gradient contributions over snloc
        dkdn_contracted(j) = sum(dkdn(:,j), 1)
        dkdn_contracted(j) = ( dkdn_contracted(j) ) ** 2.0
    end do

    ewrite(1,*) "dk/dn array: ", dkdn
    ewrite(1,*) "contracted dk/dn array: ", dkdn_contracted

    ! get scalar field values at quadrature points: matmul( 1*snloc, snloc*sngi )
    fields_quad = 2.0 * face_val_at_quad(kk, sele) * face_val_at_quad(EV, sele)
    ! Perform rhs surface integral at node sele
    rhs = shape_rhs( fshape_kk, detwei_bdy * dkdn_contracted * fields_quad )
    ! block to add to mass matrix diagonal
    mass = shape_shape( fshape_kk, fshape_kk, detwei_bdy )

    ewrite(1,*) "rhs: ", rhs
    ewrite(1,*) "mass: ", mass

    ! Add contributions to RHS and lumped mass matrix
    call addto( rhs_vector, new_node_list, rhs )
    call addto_diag( lumped_mass, 1, 1, new_node_list, sum(mass, 2) )

end subroutine eps_bc

!------------------------------------------------------------------------------------------

subroutine keps_allocate_fields(state)

    type(state_type), intent(inout) :: state
    type(vector_field), pointer     :: vectorField
    ewrite(1,*) "In keps_allocate_fields"

    vectorField => extract_vector_field(state, "Velocity")
    nNodes = node_count(vectorField)

    ! allocate some space for the fields we need for calculations, but are optional in the model
    call allocate(ll,      vectorField%mesh, "LengthScale")
    call allocate(EV,      vectorField%mesh, "EddyVisc")
    call allocate(ED,      vectorField%mesh, "EddyDiff")
    call allocate(tke_old, vectorField%mesh, "Old_TKE")

end subroutine keps_allocate_fields

!------------------------------------------------------------------------------------------
! not currently used
!------------------------------------------------------------------------------------------

function get_normal_element_size(ele, sele, positions, kk) result (h)

    type(vector_field)                     :: positions
    type(scalar_field)                     :: kk
    integer                                :: ele, sele
    integer                                :: ndim, snloc
    integer, dimension(face_loc(kk, sele)) :: nodes_bdy
    real                                   :: h
    real, dimension(1,1)                   :: hb
    real, dimension(positions%dim,1)       :: n

    real, dimension(positions%dim,positions%dim)                     :: G
    real, dimension(positions%dim, positions%dim, ele_ngi(kk, sele)) :: invJ
    real, dimension(positions%dim, face_ngi(kk, sele))               :: normal_bdy
    real, dimension(face_ngi(kk, sele))                              :: detwei_bdy

    ndim        = positions%dim
    snloc       = face_loc(kk, sele)
    nodes_bdy = face_global_nodes(kk, sele)

    call compute_inverse_jacobian( ele_val(positions, ele), &
            ele_shape(positions, ele), invJ )

    ! Can I do this away from the boundary?
    call transform_facet_to_physical( positions, sele, detwei_f=detwei_bdy, normal=normal_bdy )
    
    ! calculate wall-normal element mesh size
    G = matmul(transpose(invJ(:,:,1)), invJ(:,:,1))
    n(:,1) = normal_bdy(:,1)
    hb = 1. / sqrt( matmul(matmul(transpose(n), G), n) )
    h  = hb(1,1) 
end function get_normal_element_size

end module k_epsilon

