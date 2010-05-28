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

  ! old turbulent kinetic energy, lengthscale, eddy viscosity, production
  type(scalar_field), save :: tke_old, ll, EV, P, tkeovereps
  ! Minimum values of 2 fields to initialise, and empirical constants from flml
  real, save               :: eps_init, k_init, ll_max, fields_min, &
                              c_mu, c_eps_1, c_eps_2, sigma_eps, sigma_k
  integer, save            :: nnodes, count_rhs, count_mass, index
  logical, save            :: calculate_bcs, limit_length

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
!    - initialise parameters based on options
!----------

subroutine keps_init(state)

    type(state_type), intent(inout) :: state
    type(scalar_field), pointer     :: scalarField

    ewrite(1,*)'Now in k_epsilon turbulence model - keps_init'

    ! Do we have the calculate_BCs option?
    calculate_bcs = &
    have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/calculate_boundaries")

    ! Do we have the limit_lengthscale option?
    limit_length = &
    have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/limit_length")

    ! Allocate the temporary, module-level variables
    call keps_allocate_fields(state)

    ! Get the 5 model constants
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/C_mu',          C_mu, default = 0.09)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/C_eps_1',    c_eps_1, default = 1.44)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/C_eps_2',    c_eps_2, default = 1.92)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/sigma_k',    sigma_k, default = 1.0)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/sigma_eps',sigma_eps, default = 1.3)
  
    ! initialise fields (k, epsilon) with minimum values
    scalarField => extract_scalar_field(state, "TurbulentKineticEnergy")
    call get_option(trim(scalarField%option_path)// &
                    "/prognostic/initial_condition::WholeMesh/constant", k_init)
    call set(scalarField, k_init)

    ! initialise other fields
    ll_max = 1.   ! hard coded to size of backward facing step inlet
    fields_min = 1.e-12

    call set( ll, ll_max )
    !eps_init = k_init**1.5 / ll_max   ! make eps_init a dependent variable
    scalarField => extract_scalar_field(state, "TurbulentDissipation")
    call get_option(trim(scalarField%option_path)// &
                    "/prognostic/initial_condition::WholeMesh/constant", eps_init)
    call set(scalarField, eps_init)
    call set(EV, c_mu * k_init**2. / eps_init)

    ewrite(1,*) "k-epsilon parameters"
    ewrite(1,*) "--------------------------------------------"
    ewrite(1,*) "c_mu: ",     c_mu
    ewrite(1,*) "c_eps_1: ",  c_eps_1
    ewrite(1,*) "c_eps_2: ",  c_eps_2
    ewrite(1,*) "sigma_k: ",  sigma_k
    ewrite(1,*) "sigma_eps: ",sigma_eps
    ewrite(1,*) "fields_min: ",fields_min
    ewrite(1,*) "k_init: ",    k_init
    ewrite(1,*) "eps_init: ",  eps_init
    ewrite(1,*) "ll_max: ",   ll_max
    ewrite(1,*) "EV_init: ",   c_mu * k_init**2. / eps_init
    ewrite(1,*) "Calculating BCs: ", calculate_bcs
    ewrite(1,*) "limiting lengthscale on surfaces: ", limit_length
    ewrite(1,*) "--------------------------------------------"

end subroutine keps_init

!------------------------------------------------------------------------------

subroutine keps_tke(state)

    type(state_type), intent(inout)    :: state
    type(scalar_field), pointer        :: source, absorption, kk, eps
    type(vector_field), pointer        :: positions, nu
    type(tensor_field), pointer        :: kk_diff, visc
    type(element_type)                 :: shape_velo
    integer                            :: i, ele, bc
    real, allocatable, dimension(:)    :: detwei
    real, allocatable, dimension(:,:,:):: dshape_velo
    real, allocatable, dimension(:)    :: strain_ngi
    real, allocatable, dimension(:)    :: strain_loc
    real                               :: epsovertke
    integer, dimension(:), pointer   :: surface_elements, surface_node_list
    character(len=FIELD_NAME_LEN)    :: bc_type, bc_name
    type(mesh_type), pointer         :: surface_mesh

    ewrite(1,*) "In keps_tke"

    positions  => extract_vector_field(state, "Coordinate")
    nu         => extract_vector_field(state, "NonlinearVelocity")
    visc       => extract_tensor_field(state, "Viscosity")
    source     => extract_scalar_field(state, "TurbulentKineticEnergySource")
    absorption => extract_scalar_field(state, "TurbulentKineticEnergyAbsorption")
    kk_diff    => extract_tensor_field(state, "TurbulentKineticEnergyDiffusivity")
    kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps        => extract_scalar_field(state, "TurbulentDissipation")
    shape_velo  = ele_shape(nu, ele)

    ! Calculate TKE production using strain rate function
    call zero(source)
    call zero(absorption)

    do ele = 1, ele_count(nu)

        allocate( dshape_velo (ele_loc(nu, ele), ele_ngi(nu, ele), nu%dim) )
        allocate( detwei (ele_ngi(nu, ele) ) )
        allocate( strain_ngi (ele_ngi(nu, ele) ) )
        allocate( strain_loc (ele_loc(nu, ele) ) )

        call transform_to_physical( positions, ele, shape_velo, dshape=dshape_velo, detwei=detwei )

        strain_ngi = double_dot_product(dshape_velo, ele_val(nu, ele) )
        strain_loc = shape_rhs( shape_velo, detwei * strain_ngi )

        ! Sum components of tensor
        call addto(source, ele_nodes(nu, ele), abs(ele_val(EV, ele) * strain_loc) )

        deallocate( dshape_velo, detwei, strain_ngi, strain_loc )

    end do

    ! Calculate TKE absorption
    do i = 1, nnodes
        epsovertke = eps%val(i) / kk%val(i)
        call set(absorption, i, epsovertke)
    end do

    ! puts the BC boundary values in surface_eps_values:
    if (calculate_bcs) then    ! do I need this any more?

        do bc = 1, get_boundary_condition_count(kk)   ! use kk bcs for eps

            call get_boundary_condition(kk, bc, name=bc_name, type=bc_type, &
                 surface_node_list=surface_node_list, surface_mesh=surface_mesh, &
                 surface_element_list=surface_elements)

            if (bc_type == 'k_epsilon') then

                ewrite(1,*) "Calculating k bc: ", bc_name
                ! what goes in here? Lacasse and Lew set a Dirichlet BC using wall stress.
                call tke_bc(state, bc, surface_node_list)

            end if

        end do

    end if

    ! set diffusivity for KK - copied from Jon's update 28 Feb revision 12719
    call zero(kk_diff)
    do i = 1, kk_diff%dim
        call set(kk_diff, i, i, EV, scale=1. / sigma_k)
    end do

    ! finally, we need a copy of the old TKE and source for eps, so grab it before we solve.
    do i=1,nnodes
        call set(tke_old, i, max(kk%val(i), fields_min) )
        call set(P, i, source%val(i))
    end do

    ewrite_minmax(tke_old)
    ewrite_minmax(kk)
    ewrite_minmax(eps)
    ewrite_minmax(EV)
    ewrite_minmax(ll)

end subroutine keps_tke

!----------------------------------------------------------------------------------
! IMPORTANT! USE OLD VALUES OF EVERYTHING FOR CONSISTENCY
!----------------------------------------------------------------------------------

subroutine keps_eps(state)

    type(state_type), intent(inout)  :: state
    type(scalar_field), pointer      :: source, absorption, kk, eps
    type(vector_field), pointer      :: positions
    type(tensor_field), pointer      :: eps_diff
    type(mesh_type), pointer         :: surface_mesh
    type(scalar_field)               :: rhs_vector, surface_eps_values
    type(petsc_csr_matrix)           :: lumped_mass
    type(csr_sparsity), pointer      :: eps_sparsity
    real                             :: prod, diss, epsovertke
    integer                          :: i, j, bc, ele, sele
    integer, dimension(:), pointer   :: surface_elements, surface_node_list
    character(len=FIELD_NAME_LEN)    :: bc_type, bc_name

    ewrite(1,*) "In keps_eps"

    positions  => extract_vector_field(state, "Coordinate")
    kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps        => extract_scalar_field(state, "TurbulentDissipation")
    source     => extract_scalar_field(state, "TurbulentDissipationSource")
    absorption => extract_scalar_field(state, "TurbulentDissipationAbsorption")
    eps_diff   => extract_tensor_field(state, "TurbulentDissipationDiffusivity")

    do i = 1, nnodes

        ! re-construct eps at "old" timestep
        eps%val(i) = max( tke_old%val(i)**1.5 / ll%val(i), fields_min )

        ! compute RHS terms in epsilon equation
        EpsOverTke = eps%val(i) / tke_old%val(i)
        prod       = c_eps_1 * epsovertke * P%val(i)   ! kk source term at old timestep
        diss       = c_eps_2 * epsovertke

        ! Implicit form of equation
        call set(source, i, prod )
        call set(absorption, i, diss )

    end do

    ! Set diffusivity for Eps
    call zero(eps_diff)
    do i = 1, eps_diff%dim
        call set(eps_diff, i, i, EV, scale = 1. / sigma_eps)
    end do

    ! puts the BC boundary values in surface_eps_values:
    if (calculate_bcs) then

        do bc = 1, get_boundary_condition_count(eps)   ! use kk bcs for eps

            call get_boundary_condition(eps, bc, name=bc_name, type=bc_type, &
                 surface_node_list=surface_node_list, surface_mesh=surface_mesh, &
                 surface_element_list=surface_elements)

            if (bc_type == 'k_epsilon') then

                ewrite(1,*) "Calculating eps bc: ", bc_name

                ! create a sparse matrix, rhs and solution vectors
                eps_sparsity => get_csr_sparsity_firstorder( state, surface_mesh, surface_mesh )
                call allocate(lumped_mass, eps_sparsity, (/1,1/), diagonal=.true., name="lumped_mass")
                call allocate(surface_eps_values, surface_mesh, name="SurfaceValuesEps")
                call allocate(rhs_vector, surface_mesh, name="RHS")
                call zero(lumped_mass)
                call zero(rhs_vector)
                call zero(surface_eps_values)

                ewrite(1,*) "surface field size: ", size( surface_eps_values%val )

                index = 1                             ! initialise array counter

                do i = 1, size(surface_elements)      ! snloc
                    ele  = face_ele(positions, i)     ! index of the element which owns face
                    sele = surface_elements(i)        ! global face id

                    ! Calculate boundary condition at surface element i and add to petsc system
                    call eps_bc(ele, sele, kk, positions, &
                                surface_node_list, rhs_vector, lumped_mass)
                end do

                ! Solve the linear system. Important: get solver options here.
                surface_eps_values%option_path = eps%option_path
                call petsc_solve( surface_eps_values, lumped_mass, rhs_vector )

                ! apply boundary condition to epsilon field
                do j = 1, node_count(surface_eps_values)

                    i = surface_node_list(j)
                    eps%val(i) = max(surface_eps_values%val(j), eps_init)

                end do

                call deallocate( rhs_vector )
                call deallocate( lumped_mass )
                call deallocate( surface_eps_values )

            end if

        end do

    end if

    ewrite_minmax(tke_old)
    ewrite_minmax(kk)
    ewrite_minmax(eps)
    ewrite_minmax(EV)
    ewrite_minmax(ll)

end subroutine keps_eps

!----------
! eddyvisc calculates the lengthscale, and then the eddy viscosity.
! Eddy viscosity is placed in the velocity viscosity
! GLS has a "fix surface values" option (see options tree notes)
! that may stabilise some simulations by re-imposing BCs here. I don't use it.
!----------

subroutine keps_eddyvisc(state)

    type(state_type), intent(inout)  :: state
    type(tensor_field), pointer      :: eddy_visc, viscosity, background
    type(vector_field), pointer      :: positions
    type(scalar_field), pointer      :: kk, eps, scalarField
    integer                          :: i, j, stat

    ewrite(1,*) "In keps_eddyvisc"

    kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps        => extract_scalar_field(state, "TurbulentDissipation")
    positions  => extract_vector_field(state, "Coordinate")
    eddy_visc  => extract_tensor_field(state, "EddyViscosity")
    viscosity  => extract_tensor_field(state, "Viscosity")
    background => extract_tensor_field(state, "BackgroundViscosity")

    ! Calculate new eddy viscosity and lengthscale
    call set(viscosity, background)

    do i = 1, nnodes

        call set(kk, i, max(kk%val(i), fields_min) )     ! clip k field at fields_min
        call set(eps, i, max(eps%val(i), fields_min) )   ! clip epsilon field at fields_min

        ! set lengthscale at all other nodes to kk**1.5/eps
        call set(ll, i, min( kk%val(i)**1.5 / eps%val(i), ll_max ) )

        ! calculate viscosity for next step and for use in other fields
        call set( EV, i, C_mu * (kk%val(i))**2. / eps%val(i) )

        ! calculate ratio of fields (a diagnostic)
        call set( tkeovereps, i, kk%val(i) / eps%val(i) )

    end do

    ! Limit the lengthscale on surfaces
    ! IDEA: Limit lengthscale to 2x adaptivity element size tensor?
    if(limit_length) then
        call limit_lengthscale(state)
    end if

    ewrite(1,*) "Set k-epsilon eddy-diffusivity and eddy-viscosity tensors for use in other fields"
    call zero(eddy_visc)    ! zero it first as we're using an addto below

    do i = 1, eddy_visc%dim
        call set(eddy_visc, i, i, EV)   !tensors are isotropic
    end do

    ewrite_minmax(kk)
    ewrite_minmax(eps)
    ewrite_minmax(EV)
    ewrite_minmax(ll)
    ewrite_minmax(tkeovereps)

    viscosity%val = viscosity%val + eddy_visc%val

    do i = 1, viscosity%dim
        do j = 1, viscosity%dim
            if(j<i) cycle
            ewrite_minmax(viscosity%val(i,j,:))
        end do
    end do

    ! Set output on optional fields
    scalarField => extract_scalar_field(state, "LengthScale", stat)
    if(stat == 0) then
        call set(scalarField, ll) 
    end if
    scalarField => extract_scalar_field(state, "ScalarEddyViscosity", stat)
    if(stat == 0) then
        call set(scalarField, EV)  
    end if
    scalarField => extract_scalar_field(state, "TKEOverEpsilon", stat)
    if(stat == 0) then
        call set(scalarField, tkeovereps)  
    end if

end subroutine keps_eddyvisc

!------------------------------------------------------------------------------------

subroutine keps_cleanup()

    ewrite(1,*) "In keps_cleanup"

    call deallocate(ll)
    call deallocate(EV)
    call deallocate(tke_old)
    call deallocate(P)
    call deallocate(tkeovereps)

end subroutine keps_cleanup

!------------------------------------------------------------------------------------!
! Needs to be called after an adapt to reset the fields and arrays within the module !
!------------------------------------------------------------------------------------!

subroutine keps_adapt_mesh(state)

    type(state_type), intent(inout) :: state

    ewrite(1,*) "In keps_adapt_mesh"
    call keps_allocate_fields(state)   ! reallocate everything
    call keps_eddyvisc(state)
    call keps_tke(state)
    !call keps_eps(state)
    !call keps_eddyvisc(state)

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
    ! check that diffusivity is on for the two turbulent fields, and diagnostic
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
    ! source terms
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
    ! check for velocity
    if (.not.have_option("/material_phase[0]/vector_field::Velocity")) then
        FLExit("You need Velocity field for k-epsilon")
    end if
    ! these fields allow the new diffusivities/viscosities to be used in other calculations
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/&
                          &tensor_field::EddyViscosity")) then
        FLExit("You need EddyViscosity field for k-epsilon")
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

subroutine keps_allocate_fields(state)

    type(state_type), intent(inout) :: state
    type(vector_field), pointer     :: vectorField
    ewrite(1,*) "In keps_allocate_fields"

    vectorField => extract_vector_field(state, "Velocity")
    nnodes = node_count(vectorField)

    ! allocate some space for the fields we need for calculations
    call allocate(ll,         vectorField%mesh, "LengthScale")
    call allocate(EV,         vectorField%mesh, "ScalarEddyViscosity")
    call allocate(tke_old,    vectorField%mesh, "Old_TKE")
    call allocate(P,          vectorField%mesh, "Production")
    call allocate(tkeovereps, vectorField%mesh, "TKEoverEpsilon")

end subroutine keps_allocate_fields

!------------------------------------------------------------------------------------------

subroutine tke_bc(state, bc, surface_node_list)

    type(state_type), intent(in)     :: state
    type(scalar_field), pointer      :: kk
    type(mesh_type), pointer         :: surface_mesh
    character(len=FIELD_NAME_LEN)    :: bc_type, bc_name
    integer                          :: i
    integer, intent(in)              :: bc
    integer, dimension(:), pointer, intent(in)   :: surface_node_list

    kk => extract_scalar_field(state, "TurbulentKineticEnergy")

    ! Lacasse and Lew set a Dirichlet BC using wall stress.
    do i = 1, size(surface_node_list)
        call set( kk, surface_node_list(i), fields_min )
    end do

end subroutine tke_bc

!------------------------------------------------------------------------------------------

subroutine eps_bc(ele, sele, kk, positions, surface_node_list, rhs_vector, lumped_mass)

    ! 05/03/10: Chris Pain suggests a good idea.
    ! Evaluate eps = f(d(k^.5)/dn) at the Gauss quadrature points on the surface.
    ! dk/dn can be written as integral dotted with normal.
    ! Then set integral (surface_shape_fns *(eps - 2*EV*(d(k^.5)/dn) ) ) = 0.
    ! Loop over surface elements to construct a surface mass matrix.
    ! Lump the mass matrix to diagonalise it, then simply solve the system for epsilon.

    type(scalar_field), pointer                   :: kk
    type(vector_field), intent(in)                :: positions
    type(scalar_field), intent(inout)             :: rhs_vector
    type(petsc_csr_matrix), intent(inout)         :: lumped_mass
    type(element_type), pointer                   :: shape_kk, fshape_kk
    type(element_type)                            :: augmented_shape
    integer, dimension(:), pointer                :: surface_node_list
    integer, intent(in)                           :: ele, sele
    integer                                       :: i, j, snloc, quad, n
    integer, dimension(face_loc(positions, sele)) :: nodes_bdy, new_node_list
    logical                                       :: on_bound
    real, dimension(ele_ngi(kk, ele))             :: detwei
    real, dimension(face_ngi(kk, sele))           :: detwei_bdy, fields_quad, dkdn_contracted

    real, dimension(positions%dim, positions%dim, ele_ngi(kk, sele))     :: invJ
    real, dimension(positions%dim, positions%dim, face_ngi(kk, sele))    :: invJ_face
    real, dimension(ele_loc(kk, ele), ele_ngi(kk, ele), positions%dim)   :: dshape_kk
    real, dimension(ele_loc(kk, ele), face_ngi(kk, sele), positions%dim) :: fdshape_kk
    real, dimension(positions%dim, face_ngi(positions, sele))            :: normal_bdy
    real, dimension(face_loc(positions,sele))                            :: rhs
    real, dimension(face_loc(positions,sele), face_ngi(positions,sele))  :: dkdn
    real, dimension(face_loc(positions,sele), face_loc(positions,sele))  :: mass

    ! Get ids and shape functions
    quad      = face_ngi(kk, sele)    ! no. of gauss points in surface element
    snloc     = face_loc(kk, sele)    ! no. of nodes on surface element
    shape_kk  => ele_shape(kk, ele)    ! shape functions in element
    fshape_kk => face_shape(kk, sele)  ! shape functions in surface element
    nodes_bdy = face_global_nodes( positions, sele ) ! global node ids in surface element

    do i = 1, size(surface_node_list)
        ! I want the indices in surface_node_list that match nodes_bdy for my addto counter.
        do j = 1, snloc
            if (surface_node_list(i) == nodes_bdy(j)) then
                new_node_list(j) = index
                index = index + 1           ! array counter
            end if
        end do
    end do

    ! inelegant solution to the problem of nodes_bdy having descending
    ! indices, so that in above loop 2nd index comes first
    if(new_node_list(1) > new_node_list(2)) then
        n = size(new_node_list)
        new_node_list = new_node_list(n:1:-1)
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

    call deallocate(augmented_shape)

    ! Get boundary normal and transformed element quadrature weights over surface
    call transform_facet_to_physical( positions, sele, detwei_f=detwei_bdy, normal=normal_bdy )

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

    ! get scalar field values at quadrature points: matmul( 1*snloc, snloc*sngi )
    fields_quad = 2.0 * face_val_at_quad(kk, sele) * face_val_at_quad(EV, sele)

    ! Perform rhs surface integral at node sele
    rhs = shape_rhs( fshape_kk, detwei_bdy * dkdn_contracted * fields_quad )

    ! block to add to mass matrix diagonal
    mass = shape_shape( fshape_kk, fshape_kk, detwei_bdy )

    ! Add contributions to RHS and lumped mass matrix
    call addto( rhs_vector, new_node_list, rhs )
    call addto_diag( lumped_mass, 1, 1, new_node_list, sum(mass, 2) )

end subroutine eps_bc

!------------------------------------------------------------------------------------------

subroutine limit_lengthscale(state)
    
    type(state_type), intent(in)     :: state
    type(vector_field), pointer      :: positions
    type(scalar_field), pointer      :: kk, eps
    integer                          :: i, j, k, bc, nele, ele, sele
    real                             :: h, ll_limit
    type(patch_type)                 :: current_patch
    type(mesh_type), pointer         :: sur_mesh
    character(len=FIELD_NAME_LEN)    :: bc_type, bc_name
    integer, dimension(:), pointer   :: surface_node_list, surface_element_list

    positions => extract_vector_field(state, "Coordinate")
    kk        => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps       => extract_scalar_field(state, "TurbulentDissipation")

    do bc = 1, get_boundary_condition_count(kk)   ! use kk bcs for eps

        call get_boundary_condition(kk, bc, name=bc_name, type=bc_type, &
             surface_node_list=surface_node_list, surface_mesh=sur_mesh, &
             surface_element_list=surface_element_list)

        if (bc_type == 'k_epsilon') then
            ! for each node in the surface, get the elements and find the normal
            ! distance, then average (add and divide by nodes per element)

            do i = 1, size(surface_node_list)

                k = surface_node_list(i)
                current_patch = get_patch_ele(sur_mesh, i)
                nele = current_patch%count
                h = 0.

                do j = 1, nele    ! sum over local elements
                    sele = surface_element_list(current_patch%elements(j))
                    ele = face_ele(positions, sele)
                    h = h + get_normal_element_size(ele, sele, positions, kk)
                end do

                h = h / nele    ! normalise

                ! Set lengthscale limit. See Yap (1987) or Craft & Launder (1996)
                ll_limit = 0.83 * eps%val(k)**2. / kk%val(k) * &
                kk%val(k)**1.5 / 2.5 / eps%val(k) / h - 1. * &
                ( kk%val(k)**1.5 / 2.5 / eps%val(k) / h )**2.

                ! set the UPPER limit for lengthscale at surface
                ll%val(k) = min( max( ll_limit, fields_min), ll_max )
                !print *, ll_limit, ll_max, min( max( ll_limit, fields_min), ll_max )

                deallocate(current_patch%elements)

            end do

        end if

    end do
    
end subroutine limit_lengthscale

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

    ndim      = positions%dim
    snloc     = face_loc(kk, sele)
    nodes_bdy = face_global_nodes(kk, sele)

    call compute_inverse_jacobian( ele_val(positions, ele), ele_shape(positions, ele), invJ )

    ! Can I do this away from the boundary?
    call transform_facet_to_physical( positions, sele, detwei_f=detwei_bdy, normal=normal_bdy )
    
    ! calculate wall-normal element mesh size
    G = matmul(transpose(invJ(:,:,1)), invJ(:,:,1))
    n(:,1) = normal_bdy(:,1)
    hb = 1. / sqrt( matmul(matmul(transpose(n), G), n) )
    h  = hb(1,1)

end function get_normal_element_size

!------------------------------------------------------------------------------------------
! Computes the strain rate for the LES model. Double-dot product results in a scalar:
! t:s = [ t11s11 + t12s21 + t13s31 + t21s12 + t22s22 + t23s32 + t31s13 + t32s23 + t33s33 ]
!------------------------------------------------------------------------------------------

function double_dot_product(du_t, nu)

    ! derivative of velocity shape function (nloc x ngi x dim)
    real, dimension(:,:,:), intent(in):: du_t
    real, dimension(:,:), intent(in):: nu     ! nonlinear velocity (dim x nloc)
    real, dimension( size(du_t,2) )            :: double_dot_product
    real, dimension(size(du_t,3),size(du_t,3)) :: S, T
    integer dim, ngi, nloc
    integer gi, i, j

    nloc=size(du_t,1)
    ngi=size(du_t,2)
    dim=size(du_t,3)

    do gi=1, ngi

       s = matmul( nu, du_t(:,gi,:) )
       T = S
       S = S + transpose(S)
       double_dot_product = 0.

       do i = 1, dim
           do j = 1, dim
               double_dot_product(gi) = double_dot_product(gi) + T(i,j) * S(j,i)
           end do
       end do

    end do

end function double_dot_product

end module k_epsilon

