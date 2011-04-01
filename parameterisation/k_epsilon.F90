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
  use boundary_conditions
  use fields_manipulation
  use surface_integrals
  use fetools
  use vector_tools
  use sparsity_patterns_meshes
  use FLDebug

implicit none

  private

  ! old turbulent kinetic energy, lengthscale, eddy viscosity, production, field ratios
  type(scalar_field), save            :: tke_old, ll, tkeovereps, epsovertke
  ! Minimum values of 2 fields to initialise, and empirical constants from flml
  real, save                          :: eps_init, k_init, ll_max, visc, &
                                         c_mu, c_eps_1, c_eps_2, sigma_eps, sigma_k
  real, save                          :: fields_min = 1.e-20
  integer, save                       :: nnodes ! for VELOCITY - not necessarily the same for k & eps!!!
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
    scalarField => extract_scalar_field(state, "ScalarEddyViscosity")
    call set(scalarField, c_mu * k_init**2. / eps_init)

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
    ewrite(1,*) "implicit/explicit source/absorption terms: ",   trim(src_abs)
    ewrite(1,*) "limiting lengthscale on surfaces: ", limit_length
    ewrite(1,*) "calculating k bcs within module: ", do_k_bc
    ewrite(1,*) "calculating eps bcs within module: ", do_eps_bc
    ewrite(1,*) "--------------------------------------------"

end subroutine keps_init

!------------------------------------------------------------------------------

subroutine keps_tke(state)

    type(state_type), intent(inout)    :: state
    type(scalar_field), pointer        :: source_kk, absorption_kk, kk, eps, EV, lumped_mass, strain_rate
    type(scalar_field)                 :: src_rhs, abs_rhs, strain_rhs
    type(vector_field), pointer        :: positions, nu, u
    type(tensor_field), pointer        :: background_diff, kk_diff
    type(element_type), pointer        :: shape_kk
    integer                            :: i, ele
    integer, pointer, dimension(:)     :: nodes_kk
    real                               :: residual
    real, allocatable, dimension(:)    :: detwei, strain_ngi, rhs_addto, lumpedmass
    real, allocatable, dimension(:,:,:):: dshape_kk

    ewrite(1,*) "In keps_tke"

    positions       => extract_vector_field(state, "Coordinate")
    nu              => extract_vector_field(state, "NonlinearVelocity")
    u              => extract_vector_field(state, "Velocity")
    source_kk       => extract_scalar_field(state, "TurbulentKineticEnergySource")
    absorption_kk   => extract_scalar_field(state, "TurbulentKineticEnergyAbsorption")
    kk_diff         => extract_tensor_field(state, "TurbulentKineticEnergyDiffusivity")
    kk              => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps             => extract_scalar_field(state, "TurbulentDissipation")
    !background_diff => extract_tensor_field(state, "BackgroundDiffusivity")
    EV              => extract_scalar_field(state, "ScalarEddyViscosity")
    strain_rate     => extract_scalar_field(state, "StrainRate")
    lumped_mass     => extract_scalar_field(state, "LumpedMass")

    ewrite_minmax(tke_old)
    ewrite_minmax(kk)
    ewrite_minmax(source_kk)
    ewrite_minmax(absorption_kk)

    call zero(source_kk); call zero(absorption_kk); call zero(strain_rate); call zero(lumped_mass)

    ! Set copy of old kk for eps solve
    call set(tke_old, kk)

    ewrite_minmax(tke_old)
    ewrite_minmax(kk)
    ewrite_minmax(source_kk)
    ewrite_minmax(absorption_kk)
    ewrite_minmax(strain_rate)

    call allocate(src_rhs, kk%mesh, name="KKSRCRHS")
    call allocate(abs_rhs, kk%mesh, name="KKABSRHS")
    call allocate(strain_rhs, kk%mesh, name="STRAINRHS")

    call zero(src_rhs); call zero(abs_rhs); call zero(strain_rhs)

    ! Assembly loop
    do ele = 1, ele_count(kk)
        shape_kk => ele_shape(kk, ele)
        nodes_kk => ele_nodes(kk, ele)

        allocate(dshape_kk (size(nodes_kk), ele_ngi(kk, ele), positions%dim))
        allocate(detwei (ele_ngi(kk, ele)))
        allocate(strain_ngi (ele_ngi(kk, ele)))
        allocate(rhs_addto(ele_loc(kk, ele)))
        allocate(lumpedmass(ele_loc(kk, ele)))

        call transform_to_physical( positions, ele, shape_kk, dshape=dshape_kk, detwei=detwei )

        ! Lumped mass diagnostic field:
        lumpedmass = shape_rhs(shape_kk, detwei)
        call addto(lumped_mass, nodes_kk, lumpedmass)

        ! Calculate TKE production at ngi using strain rate (double_dot_product) function
        strain_ngi = double_dot_product(dshape_kk, ele_val(u, ele) )

        ! Strain rate diagnostic field:
        rhs_addto = shape_rhs(shape_kk, detwei*strain_ngi)
        call addto(strain_rhs, nodes_kk, rhs_addto)

        ! Source term:
        rhs_addto = shape_rhs(shape_kk, detwei*strain_ngi*ele_val_at_quad(EV, ele))
        call addto(src_rhs, nodes_kk, rhs_addto)

        ! Absorption term:
        rhs_addto = shape_rhs(shape_kk, detwei*ele_val_at_quad(eps,ele)/ele_val_at_quad(kk,ele))
        call addto(abs_rhs, nodes_kk, rhs_addto(:)/lumpedmass(:))

        deallocate(dshape_kk, detwei, strain_ngi, rhs_addto, lumpedmass)
    end do

    ewrite_minmax(source_kk)
    ewrite_minmax(absorption_kk)

    ! Node loop: set src/abs values at nodes
    do i = 1, node_count(kk)
        call set(strain_rate, i, node_val(strain_rhs,i))
      select case (src_abs)
      case ("explicit")
        call set(source_kk, i, node_val(src_rhs,i))
        call set(absorption_kk, i, node_val(abs_rhs,i))
      case ("implicit")
        residual = node_val(abs_rhs,i) - node_val(src_rhs,i)
        call set(source_kk, i, -min(0.0, residual) )
        call set(absorption_kk, i, max(0.0, residual) )
      case default
        FLAbort("Invalid implicitness option for k")
      end select
    end do

    ewrite_minmax(source_kk)
    ewrite_minmax(absorption_kk)

    call deallocate(src_rhs); call deallocate(abs_rhs)

    ! Set diffusivity for k equation.
    call zero(kk_diff)
    do i = 1, kk_diff%dim(1)
        call set(kk_diff, i, i, EV, scale=1. / sigma_k)
        ewrite_minmax(kk_diff%val(i,i,:))
    end do

    ! Add special boundary conditions to k field, if selected.
    call keps_bcs(state, kk)

    ewrite_minmax(tke_old)
    ewrite_minmax(kk)
    ewrite_minmax(source_kk)
    ewrite_minmax(absorption_kk)
    ewrite_minmax(strain_rate)

end subroutine keps_tke

!----------------------------------------------------------------------------------

subroutine keps_eps(state)

    type(state_type), intent(inout)    :: state
    type(scalar_field), pointer        :: source_eps, source_kk, absorption_eps, eps, EV, lumped_mass
    type(scalar_field)                 :: src_rhs, abs_rhs
    type(vector_field), pointer        :: positions
    type(tensor_field), pointer        :: background_diff, eps_diff
    type(element_type), pointer        :: shape_eps
    integer                            :: i, ele
    real                               :: residual
    real, allocatable, dimension(:)    :: detwei, rhs_addto, lumpedmass
    integer, pointer, dimension(:)     :: nodes_eps

    ewrite(1,*) "In keps_eps"
    eps             => extract_scalar_field(state, "TurbulentDissipation")
    source_eps      => extract_scalar_field(state, "TurbulentDissipationSource")
    source_kk       => extract_scalar_field(state, "TurbulentKineticEnergySource")
    absorption_eps  => extract_scalar_field(state, "TurbulentDissipationAbsorption")
    eps_diff        => extract_tensor_field(state, "TurbulentDissipationDiffusivity")
    !background_diff => extract_tensor_field(state, "BackgroundDiffusivity")
    EV              => extract_scalar_field(state, "ScalarEddyViscosity")
    positions       => extract_vector_field(state, "Coordinate")

    ewrite_minmax(tke_old)
    ewrite_minmax(source_kk)
    ewrite_minmax(eps)
    ewrite_minmax(source_eps)
    ewrite_minmax(absorption_eps)

    call allocate(src_rhs, eps%mesh, name="EPSSRCRHS")
    call allocate(abs_rhs, eps%mesh, name="EPSABSRHS")

    call zero(src_rhs); call zero(abs_rhs)
    call zero(source_eps); call zero(absorption_eps)

    ! Assembly loop
    do ele = 1, ele_count(eps)
        shape_eps => ele_shape(eps, ele)
        nodes_eps => ele_nodes(eps, ele)

        allocate(detwei (ele_ngi(eps, ele)))
        allocate(rhs_addto(ele_loc(eps, ele)))
        allocate(lumpedmass(ele_loc(eps, ele)))

        ! Get transformed quadrature weights and shape function gradients
        call transform_to_physical(positions, ele, detwei=detwei)

        lumpedmass = shape_rhs(shape_eps, detwei)

        ! Source term:
        rhs_addto = shape_rhs(shape_eps, detwei*c_eps_1*ele_val_at_quad(eps,ele)/ &
                              ele_val_at_quad(tke_old,ele)*ele_val_at_quad(source_kk,ele))
        call addto(src_rhs, nodes_eps, rhs_addto(:)/lumpedmass(:))

        ! Absorption term:
        rhs_addto = shape_rhs(shape_eps, detwei*c_eps_2*ele_val_at_quad(eps,ele)/ &
                              ele_val_at_quad(tke_old,ele))

        call addto(abs_rhs, nodes_eps, rhs_addto(:)/lumpedmass(:))

        deallocate(detwei, rhs_addto, lumpedmass)
    end do

    ! Node loop: set src/abs values at nodes
    do i = 1, node_count(eps)
      select case (src_abs)
      case ("explicit")
        call set(source_eps, i, node_val(src_rhs,i))
        call set(absorption_eps, i, node_val(abs_rhs,i))
      case ("implicit")
        residual = node_val(abs_rhs,i) - node_val(src_rhs,i)
        ! Puts source into absorption. Ensures positivity of terms.
        call set(source_eps, i, -min(0.0, residual) )
        call set(absorption_eps, i, max(0.0, residual) )
      case default
        FLAbort("Invalid implicitness option for k")
      end select
    end do

    call deallocate(src_rhs); call deallocate(abs_rhs)

    ! Set diffusivity for Eps
    call zero(eps_diff)
    do i = 1, eps_diff%dim(1)
        call set(eps_diff, i,i, EV, scale=1./sigma_eps)
        ewrite_minmax(eps_diff%val(i,i,:))
    end do

    ! Add special boundary conditions to eps field, if selected.
    call keps_bcs(state, eps)

    ewrite_minmax(tke_old)
    ewrite_minmax(source_kk)
    ewrite_minmax(eps)
    ewrite_minmax(source_eps)
    ewrite_minmax(absorption_eps)

end subroutine keps_eps

!----------
! eddyvisc calculates the lengthscale, and then the eddy viscosity.
! Eddy viscosity is added to the background viscosity.
!----------

subroutine keps_eddyvisc(state)

    type(state_type), intent(inout)  :: state
    type(tensor_field), pointer      :: eddy_visc, eddy_diff, viscosity, diffusivity, bg_visc, bg_diff
    type(vector_field), pointer      :: positions
    type(scalar_field), pointer      :: kk, eps, EV, scalarField, lumped_mass
    type(scalar_field)               :: ev_rhs
    type(element_type), pointer      :: shape_ev
    integer                          :: i, ele, stat
    integer, pointer, dimension(:)   :: nodes_ev
    real, allocatable, dimension(:)  :: detwei, rhs_addto, lumpedmass

    ewrite(1,*) "In keps_eddyvisc"
    kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps        => extract_scalar_field(state, "TurbulentDissipation")
    positions  => extract_vector_field(state, "Coordinate")
    eddy_visc  => extract_tensor_field(state, "KEpsEddyViscosity")
    eddy_diff  => extract_tensor_field(state, "KEpsEddyDiffusivity")
    viscosity  => extract_tensor_field(state, "Viscosity")
    !diffusivity=> extract_tensor_field(state, "Diffusivity")
    bg_visc    => extract_tensor_field(state, "BackgroundViscosity")
    !bg_diff    => extract_tensor_field(state, "BackgroundDiffusivity")
    EV         => extract_scalar_field(state, "ScalarEddyViscosity")

    ewrite_minmax(kk)
    ewrite_minmax(eps)
    ewrite_minmax(EV)

    call allocate(ev_rhs, EV%mesh, name="EVRHS")
    call zero(ev_rhs); call zero(ll); call zero(EV); call zero(tkeovereps); call zero(epsovertke)

    ! Initialise viscosity to background value
    call set(viscosity, bg_visc)
    !call set(diffusivity, bg_diff)

    ! Check components of viscosity
    do i = 1, viscosity%dim(1)
        ewrite_minmax(viscosity%val(i,i,:))
        !ewrite_minmax(diffusivity%val(i,i,:))
    end do

    ! Calculate scalar eddy viscosity by integration over element
    do ele = 1, ele_count(EV)
      nodes_ev => ele_nodes(EV, ele)
      shape_ev =>  ele_shape(EV, ele)
      allocate(detwei (ele_ngi(EV, ele)))
      allocate(rhs_addto (ele_loc(EV, ele)))
      allocate(lumpedmass(ele_loc(EV, ele)))

      call transform_to_physical(positions, ele, detwei=detwei)

      lumpedmass = shape_rhs(shape_ev, detwei)
      rhs_addto = shape_rhs(shape_ev, detwei*C_mu*ele_val_at_quad(kk,ele)**2./ele_val_at_quad(eps,ele))

      call addto(ev_rhs, nodes_ev, rhs_addto(:)/lumpedmass(:))

      deallocate(detwei, rhs_addto, lumpedmass)
    end do

    ! Node loop: set eddy viscosity at nodes
    do i = 1, node_count(EV)
      call set(EV, i, node_val(ev_rhs,i))
      ! Now set diagnostic fields. These do not need to be assembled by integration
      ! because we do not need to know the values at gauss points.
      call set(ll, i, min(c_mu * node_val(kk,i)**1.5 / node_val(eps,i), ll_max))
      call set( tkeovereps, i, node_val(kk,i) / node_val(eps,i) )
      call set( epsovertke, i, 1. / node_val(tkeovereps,i) )
    end do

    call deallocate(ev_rhs)

    ! Limit the lengthscale on surfaces
    if(limit_length) then
      call limit_lengthscale(state)
    end if

    ewrite(1,*) "Set k-epsilon eddy-diffusivity and eddy-viscosity tensors"
    call zero(eddy_visc)    ! zero it first as we're using an addto below
    call zero(eddy_diff)

    ! eddy tensors are isotropic
    do i = 1, eddy_visc%dim(1)
      call set(eddy_visc, i, i, EV)
      ! Using Schmidt number for k (usually 1.0), in absence of a Prandtl number option for temperature.
      call set(eddy_diff, i, i, EV, scale=1./sigma_k)
    end do

    ! Add turbulence model contribution to viscosity field
    call addto(viscosity, eddy_visc)
    !call addto(diffusivity, eddy_diff)

    ! Check components of viscosity
    do i = 1, viscosity%dim(1)
        ewrite_minmax(viscosity%val(i,i,:))
    end do

    ! Set output on optional fields
    scalarField => extract_scalar_field(state, "LengthScale", stat)
    if(stat == 0) then
        call set(scalarField, ll) 
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
    call deallocate(tke_old)
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
    !if (.not.have_option(trim(option_path)//"/&
    !                      &tensor_field::KEpsEddyDiffusivity")) then
    !    FLExit("You need KEpsEddyDiffusivity field for k-epsilon")
    !end if
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
    nnodes = node_count(vectorField) ! for VELOCITY - not necessarily the same for k & eps!!!

    ! allocate some space for the fields we need for calculations
    call allocate(ll,         vectorField%mesh, "LengthScale")
    call allocate(tke_old,    vectorField%mesh, "Old_TKE")
    call allocate(tkeovereps, vectorField%mesh, "TKEoverEpsilon")
    call allocate(epsovertke, vectorField%mesh, "EpsilonOverTKE")

end subroutine keps_allocate_fields

!--------------------------------------------------------------------------------!
! This gets and applies locally defined boundary conditions (wall functions)     !
!--------------------------------------------------------------------------------!

subroutine keps_bcs(state, field)

    type(state_type), intent(inout)            :: state
    type(scalar_field), pointer, intent(inout) :: field   ! k or epsilon
    type(scalar_field), pointer                :: surface_field, EV
    type(vector_field), pointer                :: positions, u
    type(scalar_field)                         :: rhs_field
    integer                                    :: i, j, ele, sele, node
    integer, dimension(:), pointer             :: surface_elements, surface_node_list
    character(len=FIELD_NAME_LEN)              :: bc_type, bc_name, wall_fns
    character(len=OPTION_PATH_LEN)             :: bc_path
    integer, allocatable, dimension(:)         :: nodes_bdy
    real, dimension(:), allocatable            :: rhs

    positions => extract_vector_field(state, "Coordinate")
    u         => extract_vector_field(state, "Velocity")
    EV        => extract_scalar_field(state, "ScalarEddyViscosity")

    do j = 1, get_boundary_condition_count(field)

       call get_boundary_condition(field, j, name=bc_name, type=bc_type, &
            surface_node_list=surface_node_list, surface_element_list=surface_elements)

       if (bc_type == 'k_epsilon') then

          ! Do we have high- or low-Reynolds number options for wall functions?
          bc_path=field%bc%boundary_condition(j)%option_path
          bc_path=trim(bc_path)//"/type"
          call get_option(trim(bc_path)//"/wall_functions", wall_fns)

          ewrite(1,*) "Calculating field BC: ",trim(field%name),trim(bc_name),' ',trim(bc_type),' ',trim(wall_fns)

          call allocate(rhs_field, field%mesh, name="RHS")
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
             call keps_wall_function(state, field, ele, sele, positions, u, EV, wall_fns, rhs)
             call addto(rhs_field, nodes_bdy, rhs)

             deallocate(rhs); deallocate(nodes_bdy)
          end do

          ! Set values in surface field
          surface_field => extract_surface_field(field, j, "value")
          ewrite(1,*) "Applying BC values to surface field: ", surface_field%name, surface_field%mesh%name

          do i = 1, size(surface_node_list)
             node = surface_node_list(i)
             call set( surface_field, i, rhs_field%val(node) )
          end do

          call deallocate(rhs_field)

       end if
    end do

end subroutine keps_bcs

!--------------------------------------------------------------------------------!
! Only used if bc type == k_epsilon for field.                                    !
!--------------------------------------------------------------------------------!

subroutine keps_wall_function(state, field, ele, sele, positions, u, EV, wall_fns, rhs)

    type(state_type), intent(inout)                      :: state
    type(scalar_field), pointer, intent(in)              :: field, EV
    type(vector_field), pointer, intent(in)              :: positions, u
    integer, intent(in)                                  :: ele, sele
    character(len=FIELD_NAME_LEN), intent(in)            :: wall_fns
    real, dimension(face_loc(field,sele)), intent(inout) :: rhs

    type(scalar_field), pointer                          :: kk
    type(tensor_field), pointer                          :: bg_visc
    type(element_type), pointer                          :: shape_field, fshape_field
    type(element_type)                                   :: augmented_shape
    integer                                              :: gi, sngi, snloc
    real                                                 :: kappa, h
    real, dimension(1,1)                                 :: hb
    real, dimension(ele_ngi(field,ele))                  :: detwei
    real, dimension(face_ngi(field,sele))                :: detwei_bdy, bc_value, sfield_quad
    real, dimension(face_loc(field,sele))                :: lumpedmass
    real, dimension(positions%dim,1)                     :: n
    real, dimension(positions%dim,positions%dim)         :: G
    real, dimension(positions%dim,face_ngi(field,sele))                    :: normal_bdy, grad_gi, vfield_quad
    real, dimension(positions%dim,positions%dim,ele_ngi(field,sele))       :: invJ
    real, dimension(positions%dim,positions%dim,face_ngi(field,sele))      :: invJ_face, tfield_quad
    real, dimension(ele_loc(field,ele),ele_ngi(field,ele),positions%dim)   :: dshape_field
    real, dimension(ele_loc(field,ele),face_ngi(field,sele),positions%dim) :: fdshape_field

    ! Get boundary normal and transformed element quadrature weights over surface
    call transform_facet_to_physical( positions, sele, detwei_f=detwei_bdy, normal=normal_bdy )

    ! Get ids, lists and shape functions
    sngi      =  face_ngi(field, sele)    ! no. of gauss points in surface element
    snloc     =  face_loc(field, sele)    ! no. of nodes on surface element
    shape_field  => ele_shape(field, ele)    ! scalar field shape functions in volume element
    fshape_field => face_shape(field, sele)  ! scalar field shape functions in surface element

    ! Get dshape_field and element quadrature weights: ngi
    call transform_to_physical( positions, ele, shape_field, dshape=dshape_field, detwei=detwei, invJ=invJ )

    !call compute_inverse_jacobian( ele_val(positions, ele), shape_x, invJ )
    invJ_face = spread(invJ(:, :, 1), 3, size(invJ_face, 3))

    ! Transform element shape functions to face
    augmented_shape = make_element_shape(ele_vertices(field, ele), shape_field%dim, &
                    shape_field%degree, shape_field%quadrature, quad_s=fshape_field%quadrature)

    ! nloc x sngi x dim
    fdshape_field = eval_volume_dshape_at_face_quad( augmented_shape, &
                           local_face_number(field, sele), invJ_face)

    lumpedmass = shape_rhs(fshape_field, detwei_bdy)
    call deallocate(augmented_shape)

    ! low Re wall functions for k and epsilon: see e.g. Wilcox (1994)
    if(wall_fns=="low_Re") then
       if (field%name=="TurbulentKineticEnergy") then
          rhs = 0.0
       else if (field%name=="TurbulentDissipation") then
          bg_visc => extract_tensor_field(state, "BackgroundViscosity")
          kk      => extract_scalar_field(state, "TurbulentKineticEnergy")
          tfield_quad = face_val_at_quad(bg_visc,sele)
          sfield_quad = face_val_at_quad(kk, sele)

          do gi = 1, sngi
             ! Grad_n operator
             grad_gi(:,gi) = matmul(fdshape_field(:,gi,:), normal_bdy(:,gi))
             ! Grad_n operator applied to k^0.5
             grad_gi(:,gi) = grad_gi(:,gi) * abs(sfield_quad(gi))**0.5
             ! Get magnitude of vector by taking sqrt(grad_n.grad_n). Multiply by viscosity^2.
             bc_value(gi) = norm2(grad_gi(:,gi)) * tfield_quad(1,1,gi) * 2.0
          end do

          rhs = shape_rhs(fshape_field, detwei_bdy*bc_value)
          rhs(:) = rhs(:)/lumpedmass(:)
       end if

    ! high Re shear-stress wall functions for k and epsilon: see e.g. Wilcox (1994), Mathieu p.360
    else if(wall_fns=="high_Re") then
       vfield_quad = face_val_at_quad(u,sele)
       sfield_quad = face_val_at_quad(EV,sele)

       do gi = 1, sngi
          ! Grad_n operator
          grad_gi(:,gi) = matmul(fdshape_field(:,gi,:), normal_bdy(:,gi))
          ! Multiply grad_n operator by velocity vector (dim*sngi)
          grad_gi(:,gi) = grad_gi(:,gi) * vfield_quad(:,gi)
          ! Subtract normal component of velocity, leaving tangent components:
          ! grad_n*U - n*((grad_n*U).n) (dim*sngi)
          grad_gi(:,gi) = grad_gi(:,gi) - (normal_bdy(:,gi) * dot_product(grad_gi(:,gi), normal_bdy(:,gi)))
          ! Get streamwise component by taking sqrt(grad_n.grad_n). Multiply by eddy viscosity.
          bc_value(gi) = norm2(grad_gi(:,gi)) * sfield_quad(gi)
       end do

       if (field%name=="TurbulentKineticEnergy") then
          rhs = shape_rhs(fshape_field, detwei_bdy*bc_value/C_mu**0.5)
          rhs(:) = rhs(:)/lumpedmass(:)
       else if (field%name=="TurbulentDissipation") then
          ! calculate wall-normal element size
          G = matmul(transpose(invJ(:,:,1)), invJ(:,:,1))
          n(:,1) = normal_bdy(:,1)
          hb = 1. / sqrt( matmul(matmul(transpose(n), G), n) )
          h  = hb(1,1)
          ! Von Karman's constant
          kappa = 0.43

          rhs = shape_rhs(fshape_field, detwei_bdy*bc_value**1.5/kappa/h)
          rhs(:) = rhs(:)/lumpedmass(:)
       end if
    else
       FLAbort("Unknown wall function option for k_epsilon boundary conditions!")
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
       !ewrite(2,*) "nu,du_t: ", nu, du_t(:,gi,:)
       S = matmul( nu, du_t(:,gi,:) )
       T = S
       !print *, "gi, S: ", gi, S
       S = S + transpose(S)
       double_dot_product(gi) = 0.

       do i = 1, dim
           do j = 1, dim
               double_dot_product(gi) = double_dot_product(gi) + T(i,j) * S(j,i)
               !ewrite(2,*) "i,j,T_ij,S_ji,ddp: ", i,j,T(i,j),S(j,i),double_dot_product(gi)
           end do
       end do
    end do

end function double_dot_product

end module k_epsilon
