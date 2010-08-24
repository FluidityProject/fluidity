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
  use field_options
  use state_module
  use spud
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use state_fields_module
  use populate_state_module
  use boundary_conditions
  use boundary_conditions_from_options
  use fields_manipulation
  use surface_integrals
  use fetools
  use vector_tools
  use FLDebug

implicit none

  private

  ! old turbulent kinetic energy, lengthscale, eddy viscosity, production, field ratios
  type(scalar_field), save            :: tke_old, ll, EV, P, tkeovereps, epsovertke
  ! Minimum values of 2 fields to initialise, and empirical constants from flml
  real, save                          :: eps_init, k_init, ll_max, visc, &
                                         c_mu, c_eps_1, c_eps_2, sigma_eps, sigma_k
  real, save                          :: fields_min = 1.e-20
  integer, save                       :: nnodes
  logical, save                       :: limit_length, do_k_bc, do_eps_bc
  character(len=FIELD_NAME_LEN), save :: src_abs

  public :: keps_init, keps_cleanup, keps_tke, keps_eps, keps_eddyvisc, keps_adapt_mesh, keps_check_options

  ! Outline:
  !  - Init in populate_state.
  !  - call keps_tke (which calculates production and sets source/absorption/diffusivity for solve).
  !  - call keps_eps (which sets source/absorption/diffusivity for solve).
  !  - After keps_eps solve, keps_eddyvisc recalculates the eddy viscosity and adds it to the viscosity field.
  !  - Wall functions are added to selected boundaries in keps_bcs and wall_functions.
  !  - keps_adapt_options repopulates the fields after an adapt.
  !  - When done, clean-up.
  !
  ! TurbulentKineticEnergy and TurbulentDissipation need to have higher priority then other prognostic 
  ! fields such as temperature, velocity and salinity for this to work.
  ! Priority is set in preprocessor/Field_Priority_Lists

contains

!----------
! initialise parameters based on options
!----------

subroutine keps_init(state)

    type(state_type), intent(inout) :: state
    type(scalar_field), pointer     :: scalarField
    character(len=FIELD_NAME_LEN)   :: bc_type
    character(len=OPTION_PATH_LEN)  :: keps_path
    integer                         :: i

    ewrite(1,*)'Now in k_epsilon turbulence model - keps_init'
    keps_path = "/material_phase[0]/subgridscale_parameterisations/k-epsilon"

    ! Allocate the temporary, module-level variables
    call keps_allocate_fields(state)

    ! Are source and absorption terms implicit/explicit?
    call get_option(trim(keps_path)//'/source_absorption', src_abs)

    ! Do we have the limit_lengthscale option?
    limit_length = &
    have_option(trim(keps_path)//'/limit_length')
    ! Get maximum lengthscale. This should be a geometric constraint on eddy size.
    call get_option(trim(keps_path)//'/L_max', ll_max)
    call set(ll, ll_max)

    ! Get the 5 model constants
    call get_option(trim(keps_path)//'/C_mu', C_mu, default = 0.09)
    call get_option(trim(keps_path)//'/C_eps_1', c_eps_1, default = 1.44)
    call get_option(trim(keps_path)//'n/C_eps_2', c_eps_2, default = 1.92)
    call get_option(trim(keps_path)//'/sigma_k', sigma_k, default = 1.0)
    call get_option(trim(keps_path)//'n/sigma_eps', sigma_eps, default = 1.3)

    ! initialise k field with minimum values
    scalarField => extract_scalar_field(state, "TurbulentKineticEnergy")
    call get_option(trim(scalarField%option_path)// &
                    "/prognostic/initial_condition::WholeMesh/constant", k_init)
    call set(scalarField, k_init)

    ! Are we calculating k boundary conditions within this module?
    do_k_bc = .false.
    do i = 1, get_boundary_condition_count(scalarField)
       call get_boundary_condition(scalarField, i, type=bc_type)
       if (bc_type=="k_epsilon") then
          do_k_bc =.true.
          exit
       end if
    end do

    ! initialise epsilon field with minimum values
    scalarField => extract_scalar_field(state, "TurbulentDissipation")
    call get_option(trim(scalarField%option_path)// &
                    "/prognostic/initial_condition::WholeMesh/constant", eps_init)
    call set(scalarField, eps_init)

    ! Are we calculating epsilon boundary conditions within this module?
    do i = 1, get_boundary_condition_count(scalarField)
       call get_boundary_condition(scalarField, i, type=bc_type)
       if (bc_type=="k_epsilon") then
          do_eps_bc =.true.
          exit
       end if
    end do

    ! Get background viscosity
    call get_option("/material_phase::Fluid/subgridscale_parameterisations/&
                          &k-epsilon/tensor_field::BackgroundViscosity/prescribed/&
                          &value::WholeMesh/isotropic/constant", visc)
    

    ! initialise eddy viscosity field with minimum values
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
    ewrite(1,*) "background visc: ",   visc
    ewrite(1,*) "limiting lengthscale on surfaces: ", limit_length
    ewrite(1,*) "calculating k bcs within module: ", do_k_bc
    ewrite(1,*) "calculating eps bcs within module: ", do_eps_bc
    ewrite(1,*) "--------------------------------------------"

end subroutine keps_init

!------------------------------------------------------------------------------

subroutine keps_tke(state)

    type(state_type), intent(inout)    :: state
    type(scalar_field), pointer        :: source, absorption, kk, eps
    type(vector_field), pointer        :: positions, nu
    type(tensor_field), pointer        :: background_diff, kk_diff
    type(element_type)                 :: shape_kk
    integer                            :: i, ele
    real, allocatable, dimension(:)    :: detwei, strain_ngi, strain_loc
    real, allocatable, dimension(:,:,:):: dshape_kk
    real                               :: residual

    ewrite(1,*) "In keps_tke"

    positions       => extract_vector_field(state, "Coordinate")
    nu              => extract_vector_field(state, "NonlinearVelocity")
    source          => extract_scalar_field(state, "TurbulentKineticEnergySource")
    absorption      => extract_scalar_field(state, "TurbulentKineticEnergyAbsorption")
    kk_diff         => extract_tensor_field(state, "TurbulentKineticEnergyDiffusivity")
    kk              => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps             => extract_scalar_field(state, "TurbulentDissipation")
    background_diff => extract_tensor_field(state, "BackgroundDiffusivity")

    call zero(source)
    call zero(absorption)

    do ele = 1, ele_count(nu)
        shape_kk =  ele_shape(kk, ele)
        allocate( dshape_kk (ele_loc(nu, ele), ele_ngi(nu, ele), nu%dim) )
        allocate( detwei (ele_ngi(nu, ele) ) )
        allocate( strain_ngi (ele_ngi(nu, ele) ) )
        allocate( strain_loc (ele_loc(nu, ele) ) )

        call transform_to_physical( positions, ele, shape_kk, dshape=dshape_kk, detwei=detwei )

        ! Calculate TKE production using strain rate (double_dot_product) function
        strain_ngi = double_dot_product(dshape_kk, ele_val(nu, ele) )
        strain_loc = shape_rhs( shape_kk, detwei * strain_ngi )

        ! Sum components of tensor. Ensure non-negative.
        call addto(source, ele_nodes(nu, ele), max(ele_val(EV, ele) * strain_loc, 0.) )

        deallocate( dshape_kk, detwei, strain_ngi, strain_loc )
    end do

    ! Calculate TKE source and absorption. Implicit or explicit.
    select case (src_abs)
    case ("explicit")
        do i = 1, nnodes
            call set(absorption, i, eps%val(i) / kk%val(i))
            ! Get a copy of old source, for eps equation
            call set(P,  i, source%val(i))
            ! Get a copy of old k for eps equation, before we solve k.
            call set(tke_old, i, max(kk%val(i), fields_min) )
        end do
    case ("implicit")
        do i = 1, nnodes
            ! Dissipation - production
            residual = eps%val(i) / kk%val(i) - source%val(i)
            ! Puts source into absorption. Ensures positivity of terms.
            call set(source, i, -min(0.0, residual) )
            call set(absorption, i, max(0.0, residual) )
            ! Get a copy of old source, for eps equation
            call set(P,  i, source%val(i))
            ! Get a copy of old k for eps equation, before we solve k.
            call set(tke_old, i, max(kk%val(i), fields_min) )
        end do
    case default
        FLAbort("Invalid implicitness option for k")
    end select

    ! Set diffusivity for k equation.
    call zero(kk_diff)
    do i = 1, kk_diff%dim
        call set(kk_diff, i, i, EV, scale=1. / sigma_k)
        ewrite_minmax(kk_diff%val(i,i,:))
    end do

    ! Add special boundary conditions to k field, if selected.
    call keps_bcs(state, kk)

    ewrite_minmax(tke_old)
    ewrite_minmax(P)
    ewrite_minmax(kk)

end subroutine keps_tke

!----------------------------------------------------------------------------------

subroutine keps_eps(state)

    type(state_type), intent(inout)    :: state
    type(scalar_field), pointer        :: source, absorption, eps
    type(tensor_field), pointer        :: background_diff, eps_diff
    real                               :: residual
    integer                            :: i

    ewrite(1,*) "In keps_eps"
    eps             => extract_scalar_field(state, "TurbulentDissipation")
    source          => extract_scalar_field(state, "TurbulentDissipationSource")
    absorption      => extract_scalar_field(state, "TurbulentDissipationAbsorption")
    eps_diff        => extract_tensor_field(state, "TurbulentDissipationDiffusivity")
    background_diff => extract_tensor_field(state, "BackgroundDiffusivity")

    ! Calculate TKE source and absorption. Implicit or explicit.
    select case (src_abs)
    case ("explicit")
        do i = 1, nnodes
            ! compute RHS terms in epsilon equation
            call set(absorption, i, c_eps_2 * eps%val(i) / tke_old%val(i))
            ! P is kk source term at old timestep
            call set(source, i, c_eps_1 * eps%val(i) / tke_old%val(i) * P%val(i) )
        end do
    case ("implicit")
        do i = 1, nnodes
            ! compute RHS terms in epsilon equation
            residual = c_eps_2 * eps%val(i) / tke_old%val(i) - &
                       c_eps_1 * eps%val(i) / tke_old%val(i) * P%val(i)
            ! Puts source into absorption. Ensures positivity of terms.
            call set(source, i, -min(0.0, residual) )
            call set(absorption, i, max(0.0, residual) )
        end do
    case default
        FLAbort("Invalid implicitness option for eps")
    end select

    ewrite_minmax(tke_old)
    ewrite_minmax(eps)
    ewrite_minmax(P)
    ewrite_minmax(source)
    ewrite_minmax(absorption)

    ! Set diffusivity for Eps
    call zero(eps_diff)
    do i = 1, eps_diff%dim
        call set(eps_diff, i,i, EV, scale=1./sigma_eps)
        ewrite_minmax(eps_diff%val(i,i,:))
    end do

    ! Add special boundary conditions to eps field, if selected.
    call keps_bcs(state, eps)

    ewrite_minmax(source)
    ewrite_minmax(absorption)
    ewrite_minmax(eps)

end subroutine keps_eps

!----------
! eddyvisc calculates the lengthscale, and then the eddy viscosity.
! Eddy viscosity is added to the background viscosity.
!----------

subroutine keps_eddyvisc(state)

    type(state_type), intent(inout)  :: state
    type(tensor_field), pointer      :: eddy_visc, eddy_diff, viscosity, bg_visc, bg_diff
    type(vector_field), pointer      :: positions
    type(scalar_field), pointer      :: kk, eps, scalarField
    integer                          :: i, stat

    ewrite(1,*) "In keps_eddyvisc"
    kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps        => extract_scalar_field(state, "TurbulentDissipation")
    positions  => extract_vector_field(state, "Coordinate")
    eddy_visc  => extract_tensor_field(state, "KEpsEddyViscosity")
    eddy_diff  => extract_tensor_field(state, "KEpsEddyDiffusivity")
    viscosity  => extract_tensor_field(state, "Viscosity")
    bg_visc    => extract_tensor_field(state, "BackgroundViscosity")
    bg_diff    => extract_tensor_field(state, "BackgroundDiffusivity")

    ewrite_minmax(kk)
    ewrite_minmax(eps)
    ewrite_minmax(EV)

    ! Calculate new eddy viscosity and lengthscale
    call set(viscosity, bg_visc)

    do i = 1, nnodes
        call set(kk, i, max(kk%val(i), fields_min) )     ! clip k field at fields_min
        call set(eps, i, max(eps%val(i), fields_min) )   ! clip epsilon field at fields_min

        ! set lengthscale
        call set(ll, i, min( c_mu * kk%val(i)**1.5 / eps%val(i), ll_max ) )

        ! calculate ratios of fields (diagnostic)
        call set( tkeovereps, i, kk%val(i) / eps%val(i) )
        call set( epsovertke, i, 1. / tkeovereps%val(i) )

        ! calculate eddy viscosity
        call set( EV, i, C_mu * (kk%val(i))**2. / eps%val(i) )
    end do

    ! Limit the lengthscale on surfaces
    if(limit_length) then
        call limit_lengthscale(state)
    end if

    ewrite(1,*) "Set k-epsilon eddy-diffusivity and eddy-viscosity tensors"
    call zero(eddy_visc)    ! zero it first as we're using an addto below
    call zero(eddy_diff)

    !tensors are isotropic
    do i = 1, eddy_visc%dim
        call set(eddy_visc, i, i, EV)
        ! Using Schmidt number for k (usually 1.0), in absence of a better idea.
        call set(eddy_diff, i, i, EV, scale=1./sigma_k)
    end do

    call addto(eddy_visc, bg_visc)
    call addto(eddy_diff, bg_diff)

    ! Check components of viscosity
    do i = 1, viscosity%dim
        ewrite_minmax(viscosity%val(i,i,:))
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
    scalarField => extract_scalar_field(state, "EpsilonOverTKE", stat)
    if(stat == 0) then
        call set(scalarField, epsovertke)  
    end if

    ewrite_minmax(kk)
    ewrite_minmax(eps)
    ewrite_minmax(EV)

end subroutine keps_eddyvisc

!------------------------------------------------------------------------------------

subroutine keps_cleanup()

    ewrite(1,*) "In keps_cleanup"

    call deallocate(ll)
    call deallocate(EV)
    call deallocate(tke_old)
    call deallocate(P)
    call deallocate(tkeovereps)
    call deallocate(epsovertke)

end subroutine keps_cleanup

!------------------------------------------------------------------------------------!
! Called after an adapt to reset the fields and arrays within the module             !
!------------------------------------------------------------------------------------!

subroutine keps_adapt_mesh(state)

    type(state_type), intent(inout) :: state

    ewrite(1,*) "In keps_adapt_mesh"
    call keps_allocate_fields(state)
    call keps_eddyvisc(state)
    call keps_tke(state)
    call keps_eps(state)
    call keps_eddyvisc(state)

end subroutine keps_adapt_mesh

!---------------------------------------------------------------------------------

subroutine keps_check_options

    character(len=OPTION_PATH_LEN) :: option_path
    character(len=FIELD_NAME_LEN)  :: kmsh, emsh, vmsh
    integer                        :: dimension

    ewrite(1,*) "In keps_check_options"
    option_path = "/material_phase[0]/subgridscale_parameterisations/k-epsilon"

    ! one dimensional problems not supported
    call get_option("/geometry/dimension/", dimension) 
    if (dimension .eq. 1 .and. have_option(trim(option_path))) then
        FLExit("k-epsilon model is only supported for dimension > 1")
    end if
    ! Don't do k-epsilon if it's not included in the model!
    if (.not.have_option(trim(option_path))) return

    ! checking for required fields
    if (.not.have_option(trim(option_path)//"/scalar_field::TurbulentKineticEnergy")) then
        FLExit("You need TurbulentKineticEnergy field for k-epsilon")
    end if
    if (.not.have_option(trim(option_path)//"/scalar_field::TurbulentDissipation")) then
        FLExit("You need TurbulentDissipation field for k-epsilon")
    end if
    ! Check that TurbulentKineticEnergy and TurbulentDissipation fields are on the same mesh
    call get_option(trim(option_path)//"/scalar_field::TurbulentKineticEnergy/prognostic/mesh", kmsh)
    call get_option(trim(option_path)//"/scalar_field::TurbulentDissipation/prognostic/mesh", emsh)
    call get_option("/material_phase::Fluid/vector_field::Velocity/prognostic/mesh", vmsh)
    if(.not. kmsh==emsh .or. .not. kmsh==vmsh .or. .not. emsh==vmsh) then
        FLExit("You must use the Velocity mesh for TurbulentKineticEnergy and TurbulentDissipation fields")
    end if
    ! check that diffusivity is on for the two turbulent fields, and diagnostic
    if (.not.have_option(trim(option_path)//"/scalar_field::TurbulentKineticEnergy/&
        &prognostic/tensor_field::Diffusivity")) then
        FLExit("You need TurbulentKineticEnergy Diffusivity field for k-epsilon")
    end if    
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &tensor_field::Diffusivity/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentKineticEnergy Diffusivity field set to diagnostic/internal")
    end if
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &tensor_field::Diffusivity")) then
        FLExit("You need TurbulentDissipation Diffusivity field for k-epsilon")
    end if
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &tensor_field::Diffusivity/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentDissipation Diffusivity field set to diagnostic/internal")
    end if
    ! source terms
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &tensor_field::Source")) then
        FLExit("You need TurbulentKineticEnergy Source field for k-epsilon")
    end if    
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &tensor_field::Source/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentKineticEnergy Source field set to diagnostic/internal")
    end if
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &tensor_field::Source")) then
        FLExit("You need TurbulentDissipation Source field for k-epsilon")
    end if
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &tensor_field::Source/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentDissipation Source field set to diagnostic/internal")
    end if
    ! absorption terms
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &tensor_field::Absorption")) then
        FLExit("You need TurbulentKineticEnergy Absorption field for k-epsilon")
    end if    
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &tensor_field::Absorption/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentKineticEnergy Absorption field set to diagnostic/internal")
    end if
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &tensor_field::Absorption")) then
        FLExit("You need TurbulentDissipation Absorption field for k-epsilon")
    end if
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &tensor_field::Absorption/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentDissipation Absorption field set to diagnostic/internal")
    end if
    ! check for velocity
    if (.not.have_option("/material_phase[0]/vector_field::Velocity")) then
        FLExit("You need Velocity field for k-epsilon")
    end if
    ! these fields allow the new diffusivities/viscosities to be used in other calculations
    if (.not.have_option(trim(option_path)//"/&
                          &tensor_field::KEpsEddyViscosity")) then
        FLExit("You need KEpsEddyViscosity field for k-epsilon")
    end if
    if (.not.have_option(trim(option_path)//"/&
                          &tensor_field::KEpsEddyDiffusivity")) then
        FLExit("You need KEpsEddyDiffusivity field for k-epsilon")
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
    ! check that the background viscosity is an isotropic constant
    if (.not.have_option(trim(option_path)//"/tensor_field::BackgroundViscosity/prescribed/&
                          &value::WholeMesh/isotropic/constant")) then
        FLExit("You need to switch the BackgroundViscosity field to isotropic/constant")
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
    call allocate(epsovertke, vectorField%mesh, "EpsilonOverTKE")

end subroutine keps_allocate_fields

!--------------------------------------------------------------------------------!
! This gets and applies locally defined boundary conditions (wall functions)     !
!--------------------------------------------------------------------------------!

subroutine keps_bcs(state, field)

    type(state_type)                        :: state
    type(scalar_field), pointer, intent(inout) :: field   ! k or epsilon
    type(scalar_field), pointer             :: kk, surface_field   ! always need kk
    type(vector_field), pointer             :: positions
    type(scalar_field)                      :: rhs_field
    integer                                 :: i, j, ele, sele, node
    integer, dimension(:), pointer          :: surface_elements, surface_node_list
    character(len=FIELD_NAME_LEN)           :: bc_type, bc_name, wall_fns
    character(len=OPTION_PATH_LEN)          :: bc_path
    integer, allocatable, dimension(:)      :: nodes_bdy
    real, dimension(:), allocatable         :: rhs

    kk        => extract_scalar_field(state, "TurbulentKineticEnergy")
    positions => extract_vector_field(state, "Coordinate")

    do j = 1, get_boundary_condition_count(field)

       call get_boundary_condition(field, j, name=bc_name, type=bc_type, &
            surface_node_list=surface_node_list, surface_element_list=surface_elements)

       if (bc_type == 'k_epsilon') then

          ! Do we have high- or low-Reynolds number options for wall functions?
          bc_path=field%bc%boundary_condition(j)%option_path
          bc_path=trim(bc_path)//"/type"
          call get_option(trim(bc_path)//"/wall_functions", wall_fns)

          ewrite(1,*) "Calculating field BC: ", trim(field%name), &
          trim(bc_name), ', ', trim(bc_type)

          call allocate(rhs_field, kk%mesh, name="RHS")
          call zero(rhs_field)

          ! Surface element loop
          ewrite(1,*) "Entering field bc loop"
          do i = 1, size(surface_elements)
             sele = surface_elements(i)
             ele  = face_ele(field, sele)
             allocate(rhs(face_loc(field, sele) ))
             allocate(nodes_bdy(face_loc(field, sele)))
             nodes_bdy =  face_global_nodes(field, sele)

             ! Calculate wall function
             call keps_wall_function(field, kk, ele, sele, positions, wall_fns, rhs)
             call addto( rhs_field, nodes_bdy, rhs )

             deallocate(rhs); deallocate(nodes_bdy)
          end do

          ! Set values in surface field
          ewrite(1,*) "Applying BC values to surface field"
          if (associated(surface_field)) then
             surface_field => extract_surface_field(field, j, "value")
             do i = 1, size(surface_node_list)
                node = surface_node_list(i)
                call set( surface_field, i, rhs_field%val(node) )
             end do
          else
             ewrite(1,*) "No surface fields associated!"
          end if
          call deallocate(rhs_field)

       end if
    end do

end subroutine keps_bcs

!--------------------------------------------------------------------------------!
! Only used if bc type = k_epsilon for field.                                    !
!--------------------------------------------------------------------------------!

subroutine keps_wall_function(field, kk, ele, sele, positions, wall_fns, rhs)

    type(scalar_field), pointer, intent(in)              :: field, kk
    type(vector_field), pointer, intent(in)              :: positions
    integer, intent(in)                                  :: ele, sele
    character(len=FIELD_NAME_LEN), intent(in)            :: wall_fns
    real, dimension(face_loc(field,sele)), intent(inout) :: rhs
    type(element_type), pointer                          :: shape_field, fshape_field
    type(element_type)                                   :: augmented_shape
    integer                                              :: i, j, snloc, sngi
    real                                                 :: kappa, h
    real, dimension(1,1)                                 :: hb
    real, dimension(ele_ngi(field,ele))                  :: detwei
    real, dimension(face_ngi(field,sele))                :: detwei_bdy,fields_sngi,gradn_c
    real, dimension(positions%dim,1)                     :: n
    real, dimension(positions%dim,positions%dim)         :: G
    real, dimension(positions%dim,positions%dim,ele_ngi(field,sele))       :: invJ
    real, dimension(positions%dim,positions%dim,face_ngi(field,sele))      :: invJ_face
    real, dimension(ele_loc(field,ele),ele_ngi(field,ele),positions%dim)   :: dshape_field
    real, dimension(ele_loc(field,ele),face_ngi(field,sele),positions%dim) :: fdshape_field
    real, dimension(positions%dim,face_ngi(field,sele))                    :: normal_bdy
    real, dimension(face_loc(field,sele),face_ngi(field,sele))             :: gradn


    ! Get boundary normal and transformed element quadrature weights over surface
    call transform_facet_to_physical( positions, sele, detwei_f=detwei_bdy, normal=normal_bdy )

    ! Get ids, lists and shape functions
    sngi      =  face_ngi(field, sele)    ! no. of gauss points in surface element
    snloc     =  face_loc(field, sele)    ! no. of nodes on surface element
    shape_field  => ele_shape(field, ele)    ! shape functions in volume element
    fshape_field => face_shape(field, sele)  ! shape functions in surface element

    ! Get dshape_field and element quadrature weights: ngi
    call transform_to_physical( positions, ele, shape_field, dshape=dshape_field, detwei=detwei )

    call compute_inverse_jacobian( ele_val(positions, ele), shape_field, invJ )
    invJ_face = spread(invJ(:, :, 1), 3, size(invJ_face, 3))

    ! Transform element shape functions to face
    augmented_shape = make_element_shape(shape_field%loc, shape_field%dim, &
                    shape_field%degree, shape_field%quadrature, quad_s=fshape_field%quadrature )

    fdshape_field = eval_volume_dshape_at_face_quad( augmented_shape, &
                           local_face_number(field, sele), invJ_face )

    call deallocate(augmented_shape)

    ! dot normal with fdshape_field to get shape fn gradient w.r.t. surface normal
    do i = 1, snloc
       do j = 1, sngi
          gradn(i,j) = dot_product( normal_bdy(:,j), fdshape_field(i,j,:) )
       end do
    end do

    ! low Re wall functions for k and epsilon: see e.g. Wilcox (1994)
    if(wall_fns=="low_Re") then

       if (field%name=="TurbulentKineticEnergy") then
          ! Set Dirichlet here: what value? Must be tuned to Reynolds number.
          ewrite(1,*) "Low-Reynolds number k bc: set to zero for now."
          fields_sngi = 0.0

          ! Perform rhs surface integration
          rhs = shape_rhs( fshape_field, detwei_bdy * fields_sngi )

       else if (field%name=="TurbulentDissipation") then
          ! Evaluate eps = f(d(k^.5)/dn) at the Gauss quadrature points on the surface.
          ! dk/dn can be written as integral dotted with normal.
          do j = 1, sngi      ! Contract dshape gradient over snloc and square
             gradn_c(j) = sum(gradn(:,j), 1)
             gradn_c(j) = ( gradn_c(j) ) ** 2.0
          end do

          ! get k*viscosity at quadrature points
          fields_sngi = 2.0 * face_val_at_quad(kk, sele) * visc

          ! Perform rhs surface integration
          rhs = shape_rhs( fshape_field, detwei_bdy * gradn_c * fields_sngi )
       end if

    ! high Re shear-stress wall functions for k and epsilon: see e.g. Wilcox (1994)
    else
       ! Evaluate stress = grad U.n at the Gauss quadrature points on the surface.
       do j = 1, sngi      ! Contract dshape gradient over snloc
          gradn_c(j) = sum(gradn(:,j), 1)
       end do

       if (field%name=="TurbulentKineticEnergy") then
          ! calculate boundary value
          fields_sngi = 1./C_mu**0.5 * visc * gradn_c

       else if (field%name=="TurbulentDissipation") then
          ! calculate wall-normal element size
          G = matmul(transpose(invJ(:,:,1)), invJ(:,:,1))
          n(:,1) = normal_bdy(:,1)
          hb = 1. / sqrt( matmul(matmul(transpose(n), G), n) )
          h  = hb(1,1)

          ! Von Karman's constant
          kappa = 0.41

          ! get boundary condition values at quadrature points
          fields_sngi = (visc * gradn_c)**1.5 / kappa / h
       end if

       ! Perform rhs surface integration
       rhs = shape_rhs( fshape_field, detwei_bdy * fields_sngi )
    end if

end subroutine keps_wall_function

!----------------------------------------------------------------------------------!
! Use the lengthscale limiter on surfaces (see Yap (1987), Craft & Launder(1996)). !
!----------------------------------------------------------------------------------!

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

        call get_boundary_condition(eps, bc, name=bc_name, type=bc_type, &
             surface_node_list=surface_node_list, surface_mesh=sur_mesh, &
             surface_element_list=surface_element_list)

        ! for each node in the surface, get the elements and find the normal
        ! distance, then average (add and divide by nodes per element)
        if (bc_type == 'k_epsilon') then

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

                h = h / nele      ! normalise

                ! Set lengthscale limit. See Yap (1987) or Craft & Launder (1996)
                ll_limit = 0.83 * eps%val(k)**2. / kk%val(k) * &
                (kk%val(k)**1.5 / 2.5 / eps%val(k) / h - 1. ) * &
                ( kk%val(k)**1.5 / 2.5 / eps%val(k) / h )**2.

                ! set the UPPER limit for lengthscale at surface
                ll%val(k) = min( max( ll_limit, fields_min), ll_max )

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
    real                                   :: h
    real, dimension(1,1)                   :: hb
    real, dimension(positions%dim,1)       :: n

    real, dimension(positions%dim,positions%dim)                     :: G
    real, dimension(positions%dim, positions%dim, ele_ngi(kk, sele)) :: invJ
    real, dimension(positions%dim, face_ngi(kk, sele))               :: normal_bdy
    real, dimension(face_ngi(kk, sele))                              :: detwei_bdy

    call compute_inverse_jacobian( ele_val(positions, ele), ele_shape(positions, ele), invJ )
    call transform_facet_to_physical( positions, sele, detwei_f=detwei_bdy, normal=normal_bdy )
    
    ! calculate wall-normal element mesh size
    G = matmul(transpose(invJ(:,:,1)), invJ(:,:,1))
    n(:,1) = normal_bdy(:,1)
    hb = 1. / sqrt( matmul(matmul(transpose(n), G), n) )
    h  = hb(1,1)

end function get_normal_element_size

!------------------------------------------------------------------------------------------!
! Computes the strain rate for the LES model. Double-dot product results in a scalar:      !
! t:s = [ t11s11 + t12s21 + t13s31 + t21s12 + t22s22 + t23s32 + t31s13 + t32s23 + t33s33 ] !
!------------------------------------------------------------------------------------------!

function double_dot_product(du_t, nu)

    real, dimension(:,:,:), intent(in)         :: du_t   ! derivative of velocity shape function
    real, dimension(:,:), intent(in)           :: nu     ! nonlinear velocity
    real, dimension( size(du_t,2) )            :: double_dot_product
    real, dimension(size(du_t,3),size(du_t,3)) :: S, T
    integer                                    :: dim, ngi, gi, i, j

    ngi  = size(du_t,2)
    dim  = size(du_t,3)

    do gi=1, ngi

       S = matmul( nu, du_t(:,gi,:) )
       T = S
       !print *, "gi, S: ", gi, S
       S = S + transpose(S)
       double_dot_product(gi) = 0.

       do i = 1, dim
           do j = 1, dim
               double_dot_product(gi) = double_dot_product(gi) + T(i,j) * S(j,i)
           end do
       end do
    end do

end function double_dot_product

end module k_epsilon
