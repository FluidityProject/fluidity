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
  use sparse_tools
  use sparse_tools_petsc
  use petsc_tools
  use sparsity_patterns_meshes
  use parallel_tools
  use halos_base
  use halo_data_types
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
  use FLDebug

implicit none

  private

  ! old turbulent kinetic energy, lengthscale, eddy viscosity, production
  type(scalar_field), save      :: tke_old, ll, EV, P, tkeovereps
  ! Minimum values of 2 fields to initialise, and empirical constants from flml
  real, save                    :: eps_init, k_init, ll_max, &
                                   c_mu, c_eps_1, c_eps_2, sigma_eps, sigma_k
  real, save                    :: fields_min = 1.e-20
  integer, save                 :: nnodes
  logical, save                 :: limit_length, do_k_bc, do_eps_bc
  character(len=FIELD_NAME_LEN) :: src_abs

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
    character(len=FIELD_NAME_LEN)   :: bc_type
    integer                         :: i

    ewrite(1,*)'Now in k_epsilon turbulence model - keps_init'

    ! Allocate the temporary, module-level variables
    call keps_allocate_fields(state)

    ! Are source and absorption terms implicit/explicit?
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/option', src_abs)

    ! Do we have the limit_lengthscale option?
    limit_length = &
    have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/limit_length")
    ! Get maximum lengthscale. This should be a geometric constraint on eddy size.
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/L_max', ll_max)
    call set(ll, ll_max)

    ! Get the 5 model constants
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/C_mu', C_mu, default = 0.09)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/C_eps_1', c_eps_1, default = 1.44)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/C_eps_2', c_eps_2, default = 1.92)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/sigma_k', sigma_k, default = 1.0)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/sigma_eps', sigma_eps, default = 1.3)

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
    !eps_init = k_init**1.5 / ll_max   ! make eps_init a dependent variable
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
    type(element_type)                 :: shape_velo
    integer                            :: i, ele, bc
    real, allocatable, dimension(:)    :: detwei
    real, allocatable, dimension(:,:,:):: dshape_velo
    real, allocatable, dimension(:)    :: strain_ngi
    real, allocatable, dimension(:)    :: strain_loc
    real                               :: residual
    integer, dimension(:), pointer     :: surface_elements, surface_node_list
    character(len=FIELD_NAME_LEN)      :: bc_type, bc_name
    type(mesh_type), pointer           :: surface_mesh

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

        shape_velo =  ele_shape(nu, ele)          ! Velocity element type
        allocate( dshape_velo (ele_loc(nu, ele), ele_ngi(nu, ele), nu%dim) )
        allocate( detwei (ele_ngi(nu, ele) ) )
        allocate( strain_ngi (ele_ngi(nu, ele) ) )
        allocate( strain_loc (ele_loc(nu, ele) ) )

        call transform_to_physical( positions, ele, shape_velo, dshape=dshape_velo, detwei=detwei )

        ! Calculate TKE production using strain rate function
        strain_ngi = double_dot_product(dshape_velo, ele_val(nu, ele) )
        strain_loc = shape_rhs( shape_velo, detwei * strain_ngi )

        ! Sum components of tensor. Ensure non-negative.
        call addto(source, ele_nodes(nu, ele), max(ele_val(EV, ele) * strain_loc, 0.) )

        deallocate( dshape_velo, detwei, strain_ngi, strain_loc )

    end do

    ! Calculate TKE source and absorption. Implicit or explicit.
    select case (src_abs)
    case ("explicit")
        do i = 1, nnodes
            call set(absorption, i, eps%val(i) / kk%val(i))
            ! Get a copy of old source, for eps equation
            call set(P,  i, source%val(i))
            ! Get a copy of old k and source, for eps equation, before we solve k.
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
            ! Get a copy of old k and source, for eps equation, before we solve k.
            call set(tke_old, i, max(kk%val(i), fields_min) )
        end do
    case default
        FLAbort("Invalid implicitness option for k")
    end select

    ! Add special boundary conditions to k field, if selected.
    if (do_k_bc) then
       ewrite(1,*) "Entering k BC calculation"
       do bc = 1, get_boundary_condition_count(kk)

          call get_boundary_condition(kk, bc, name=bc_name, type=bc_type, &
               surface_node_list=surface_node_list, surface_mesh=surface_mesh, &
               surface_element_list=surface_elements)

          ewrite(1,*) "k bc name, type: ", bc_name, bc_type

          if (bc_type == 'k_epsilon') then
             ewrite(1,*) "Calculating k bc: ", bc_name
             ! what goes in here? Lacasse and Lew set a Dirichlet BC using wall stress.
             call tke_bc(state, surface_node_list)
          end if
       end do
    end if

    ! Set diffusivity for k equation.
    call zero(kk_diff)
    call set(kk_diff, kk_diff%dim, kk_diff%dim, EV, scale=1. / sigma_k)
    call addto(kk_diff, background_diff)

    ewrite_minmax(tke_old)
    ewrite_minmax(P)
    ewrite_minmax(kk)

end subroutine keps_tke

!----------------------------------------------------------------------------------

subroutine keps_eps(state)

    type(state_type), intent(inout)    :: state
    type(scalar_field), pointer        :: source, absorption, kk, eps
    type(vector_field), pointer        :: positions
    type(tensor_field), pointer        :: background_diff, eps_diff
    type(scalar_field)                 :: rhs_vector, surface_eps_values
    type(petsc_csr_matrix)             :: lumped_mass
    type(halo_type), pointer           :: halo
    real                               :: residual
    integer                            :: i, j, ele, sele
    integer, dimension(:), pointer     :: surface_elements, surface_node_list
    character(len=FIELD_NAME_LEN)      :: bc_type, bc_name
    integer, dimension(:), allocatable :: dnnz, onnz
    integer                            :: node, nodes, mynodes
    real, dimension(:), allocatable    :: rhs
    real, dimension(:,:), allocatable  :: mass
    integer, dimension(:), allocatable :: nodes_bdy


    ewrite(1,*) "In keps_eps"

    positions       => extract_vector_field(state, "Coordinate")
    kk              => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps             => extract_scalar_field(state, "TurbulentDissipation")
    source          => extract_scalar_field(state, "TurbulentDissipationSource")
    absorption      => extract_scalar_field(state, "TurbulentDissipationAbsorption")
    eps_diff        => extract_tensor_field(state, "TurbulentDissipationDiffusivity")
    background_diff => extract_tensor_field(state, "BackgroundDiffusivity")

    ! Calculate TKE source and absorption. Implicit or explicit.
    select case (src_abs)
    case ("explicit")
        do i = 1, nnodes
            ! re-construct eps at "old" timestep
            eps%val(i) = max( c_mu * tke_old%val(i)**2. / EV%val(i), fields_min )
            ! compute RHS terms in epsilon equation
            call set(absorption, i, c_eps_2 * eps%val(i) / tke_old%val(i))
            ! P is kk source term at old timestep
            call set(source, i, c_eps_1 * eps%val(i) / tke_old%val(i) * P%val(i) )
        end do
    case ("implicit")
        do i = 1, nnodes
            ! re-construct eps at "old" timestep
            eps%val(i) = max( c_mu * tke_old%val(i)**2. / EV%val(i), fields_min )
            ! compute RHS terms in epsilon equation
            residual = c_eps_2 * eps%val(i) / tke_old%val(i) - &
                       c_eps_1 * eps%val(i) / tke_old%val(i) * P%val(i)
            ! Puts source into absorption. Ensures positivity of terms.
            call set(source, i, -min(0.0, residual) )
            call set(absorption, i, max(0.0, residual) )
        end do
    case default
        FLAbort("Invalid implicitness option for k")
    end select

    ewrite_minmax(tke_old)
    ewrite_minmax(eps)
    ewrite_minmax(P)
    ewrite_minmax(source)
    ewrite_minmax(absorption)

    ! Set diffusivity for Eps
    call zero(eps_diff)
    call set(eps_diff, eps_diff%dim, eps_diff%dim, EV, scale=1./sigma_eps)
    call addto(eps_diff, background_diff)

    ! Add special boundary conditions to eps field, if selected.
    if(do_eps_bc) then
       ewrite(1,*) "Entering epsilon BC calculation"
       ! put the boundary values in surface_eps_values.
       nodes=node_count(eps)
       if (associated(eps%mesh%halos)) then
          halo => eps%mesh%halos(1)
          mynodes = halo_nowned_nodes(halo)
       else
          nullify(halo)
          mynodes = nodes
       end if

       ewrite(1,*) "nodes, mynodes: ", nodes, mynodes
       allocate(dnnz(1:mynodes), onnz(1:mynodes))
       dnnz = 1; onnz = 0

       do i = 1, get_boundary_condition_count(eps)
          call get_boundary_condition(eps, i, type=bc_type, &
               surface_node_list=surface_node_list)

          if (bc_type=="k_epsilon") then
             ! Check for duplicated bcs at node
             do j = 1, size(surface_node_list)
                node=surface_node_list(j)
                if (node > mynodes) cycle
                if ( dnnz(node) > 1 ) then
                   FLAbort("Two boundary condition specifications for the same node.")
                end if
                dnnz(node) = 2
             end do
          end if
       end do

       call allocate(lumped_mass, nodes, nodes, &
            dnnz, onnz, (/1,1/), "lumped_mass", halo=halo)
       call allocate(rhs_vector, eps%mesh, name="RHS")
       call allocate(surface_eps_values, eps%mesh, name="SurfaceValuesEps")
       call zero(rhs_vector)
       call zero(surface_eps_values)
       call zero(lumped_mass)

       ! put a 1.0 on the diagonal unless it's a halo node or a surface node
       do i = 1, nodes
          if (i > mynodes) cycle
          if (dnnz(i) == 2) cycle
          call addto(lumped_mass, 1, 1, i, i, 1.0)
       end do

       do i = 1, get_boundary_condition_count(eps)
          call get_boundary_condition(eps, i, type=bc_type, &
                       surface_node_list=surface_node_list, &
                       surface_element_list=surface_elements)

          if (bc_type=="k_epsilon") then
             ! Calculate boundary condition at surface element and add to system
             ewrite(1,*) "Calculating epsilon boundary condition: ", &
                  trim(bc_name), ', ', trim(bc_type), ', ', size(surface_node_list)

             do j = 1, size(surface_elements)       ! local element list
                sele = surface_elements(j)          ! global face id
                ele  = face_ele(positions, sele)    ! index of the element which owns face

                allocate(rhs(face_loc(positions, sele) ))
                allocate(mass( 1:(face_loc(positions, sele)), 1:(face_loc(positions, sele) ) ) )
                allocate(nodes_bdy(face_loc(positions, ele) ))
                nodes_bdy =  face_global_nodes(positions, sele)    ! global node ids

                call eps_bc(ele, sele, kk, positions, nodes_bdy, rhs, mass)

                ! Add contributions to RHS and lumped mass matrix using global numbering
                call addto( rhs_vector, nodes_bdy, rhs )
                call addto_diag( lumped_mass, 1, 1, nodes_bdy, sum(mass, 2) )

                deallocate(rhs); deallocate(mass); deallocate(nodes_bdy)
             end do
          end if
       end do

       call assemble(lumped_mass)

       ! matrixdump for mass matrix
       !ewrite(1,*) "writing lumped mass to file: matrixdump"
       !call dump_matrix("matrixdump", lumped_mass)
       !ewrite(1,*) "writing rhs vector to file: rhs_vector_field2file.dat"
       !call field2file('rhs_vector_field2file.dat', rhs_vector)

       ! Solve the petsc linear system. Important: get solver options here
       surface_eps_values%option_path = eps%option_path
       call petsc_solve( surface_eps_values, lumped_mass, rhs_vector )

       ! Print surface eps values:
       !call field2file('surfaceepsvalues_field2file.dat', surface_eps_values)
       !ewrite(1,*) "writing surface eps values to file: surfaceepsvalues_field2file.dat"

       ! Finally, apply all k-epsilon type boundary conditions to epsilon field
       do j = 1, nodes
           ! Ignore zeros, they are not on the boundaries of interest
           if(surface_eps_values%val(j) == 0) cycle
           eps%val(j) = max(surface_eps_values%val(j), fields_min)
       end do

       call deallocate( rhs_vector )
       call deallocate( lumped_mass )
       call deallocate( surface_eps_values )

    end if

    ewrite_minmax(source)
    ewrite_minmax(absorption)
    ewrite_minmax(eps)

end subroutine keps_eps

!----------
! eddyvisc calculates the lengthscale, and then the eddy viscosity.
! Eddy viscosity is added to the velocity viscosity.
!----------

subroutine keps_eddyvisc(state)

    type(state_type), intent(inout)  :: state
    type(tensor_field), pointer      :: eddy_visc, viscosity, background
    type(vector_field), pointer      :: positions
    type(scalar_field), pointer      :: kk, eps, scalarField
    integer                          :: i, stat

    ewrite(1,*) "In keps_eddyvisc"

    kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps        => extract_scalar_field(state, "TurbulentDissipation")
    positions  => extract_vector_field(state, "Coordinate")
    eddy_visc  => extract_tensor_field(state, "EddyViscosity")
    viscosity  => extract_tensor_field(state, "Viscosity")
    background => extract_tensor_field(state, "BackgroundViscosity")

    ewrite_minmax(kk)
    ewrite_minmax(eps)

    ! Calculate new eddy viscosity and lengthscale
    call set(viscosity, background)

    do i = 1, nnodes

        call set(kk, i, max(kk%val(i), fields_min) )     ! clip k field at fields_min
        call set(eps, i, max(eps%val(i), fields_min) )   ! clip epsilon field at fields_min

        ! set lengthscale
        call set(ll, i, min( c_mu * kk%val(i)**1.5 / eps%val(i), ll_max ) )

        ! calculate ratio of fields (a diagnostic)
        call set( tkeovereps, i, kk%val(i) / eps%val(i) )

        ! calculate eddy viscosity
        call set( EV, i, C_mu * (kk%val(i))**2. / eps%val(i) )
    end do

    ! Limit the lengthscale on surfaces
    if(limit_length) then
        call limit_lengthscale(state)
    end if

    ewrite(1,*) "Set k-epsilon eddy-diffusivity and eddy-viscosity tensors"
    call zero(eddy_visc)    ! zero it first as we're using an addto below

    do i = 1, eddy_visc%dim
        call set(eddy_visc, i, i, EV)   !tensors are isotropic
    end do

    viscosity%val = viscosity%val + eddy_visc%val
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

    ewrite_minmax(kk)
    ewrite_minmax(eps)

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

!-------------------------------------------------------------------------

!subroutine keps_init_surfaces(state)

!    type(state_type), intent(in)     :: state
!    type(scalar_field), pointer      :: kk, eps

!    ewrite(1,*) "Initialising the k-epsilon surfaces required for BCs"

    ! grab hold of some essential field
!    kk => extract_scalar_field(state, "TurbulentKineticEnergy")
!    eps => extract_scalar_field(state, "TurbulentDissipation")

    ! Create fields of surface values for k and eps wherever k_epsilon bc is applied.
!    call get_entire_boundary_condition(kk, (/"k_epsilon"/), kk_boundary_values, kk_bc_list)
!    call get_entire_boundary_condition(eps, (/"k_epsilon"/), eps_boundary_values, eps_bc_list)

!end subroutine keps_init_surfaces

!--------------------------------------------------------------------------------!
! Only used if bc type = k_epsilon for k field. Hard-coded to Dirichlet for now. !
! May be changed to a Dirichlet BC using wall stress (see Lacasse and Lew).      !
!--------------------------------------------------------------------------------!

subroutine tke_bc(state, surface_node_list)

    type(state_type), intent(in)     :: state
    type(scalar_field), pointer      :: kk
    integer                          :: i
    integer, dimension(:), pointer, intent(in)   :: surface_node_list

    kk => extract_scalar_field(state, "TurbulentKineticEnergy")

    do i = 1, size(surface_node_list)
        call set( kk, surface_node_list(i), fields_min )
    end do

end subroutine tke_bc

!-----------------------------------------------------------------------------------!
! Evaluate eps = f(d(k^.5)/dn) at the Gauss quadrature points on the surface.       !
! dk/dn can be written as integral dotted with normal.                              !
! Then set integral (surface_shape_fns *(eps - 2*EV*(d(k^.5)/dn) ) ) = 0.           !
! Loop over surface elements to construct a surface mass matrix.                    !
! Lump the mass matrix to diagonalise it, then simply solve the system for epsilon. !
!-----------------------------------------------------------------------------------!

subroutine eps_bc(ele, sele, kk, positions, nodes_bdy, rhs, mass)

    type(scalar_field), pointer                             :: kk
    type(vector_field), intent(in)                          :: positions
    type(element_type), pointer                             :: shape_kk, fshape_kk
    type(element_type)                                      :: augmented_shape
    integer, intent(in)                                     :: ele, sele
    real, dimension(face_loc(positions,sele)), intent(out)  :: rhs
    real, dimension(face_loc(positions,sele), &
                    face_loc(positions,sele)), intent(out)  :: mass
    integer                                                 :: i, j, snloc, sngi
    integer, dimension(face_loc(positions, ele))            :: nodes_bdy
    real, dimension(ele_ngi(kk, ele))                       :: detwei
    real, dimension(face_ngi(kk, sele))  :: detwei_bdy, fields_sngi, dkdn_contracted

    real, dimension(positions%dim, positions%dim, ele_ngi(kk, sele))     :: invJ
    real, dimension(positions%dim, positions%dim, face_ngi(kk, sele))    :: invJ_face
    real, dimension(ele_loc(kk, ele), ele_ngi(kk, ele), positions%dim)   :: dshape_kk
    real, dimension(ele_loc(kk, ele), face_ngi(kk, sele), positions%dim) :: fdshape_kk
    real, dimension(positions%dim, face_ngi(positions, sele))            :: normal_bdy
    real, dimension(face_loc(positions,sele), face_ngi(positions,sele))  :: dkdn


    ! Get ids, lists and shape functions
    sngi      =  face_ngi(kk, sele)    ! no. of gauss points in surface element
    snloc     =  face_loc(kk, sele)    ! no. of nodes on surface element
    shape_kk  => ele_shape(kk, ele)    ! shape functions in volume element
    fshape_kk => face_shape(kk, sele)  ! shape functions in surface element

    ! Get dshape_kk and element quadrature weights: ngi
    call transform_to_physical( positions, ele, shape_kk, dshape=dshape_kk, detwei=detwei )

    ! Need inverse Jacobian for augmented shape functions
    call compute_inverse_jacobian( ele_val(positions, ele), shape_kk, invJ )
    invJ_face = spread(invJ(:, :, 1), 3, size(invJ_face, 3))

    ! Transform element shape functions to face
    augmented_shape = make_element_shape(shape_kk%loc, shape_kk%dim, &
                      shape_kk%degree, shape_kk%quadrature, quad_s=fshape_kk%quadrature )

    ! Get dshape on face. nloc x sngi x dim
    fdshape_kk = eval_volume_dshape_at_face_quad( augmented_shape, &
                            local_face_number(kk, sele), invJ_face )

    call deallocate(augmented_shape)

    ! Get boundary normal and transformed element quadrature weights over surface
    call transform_facet_to_physical( positions, sele, detwei_f=detwei_bdy, normal=normal_bdy )

    ! dot normal with dshape_kk to get shape fn gradient w.r.t. surface normal
    do i = 1, size(nodes_bdy)
        do j = 1, sngi
            dkdn(i,j) = dot_product( normal_bdy(:,j), fdshape_kk(i,j,:) )
        end do
    end do

    do j = 1, sngi      ! Sum dshape gradient contributions over snloc
        dkdn_contracted(j) = sum(dkdn(:,j), 1)
        dkdn_contracted(j) = ( dkdn_contracted(j) ) ** 2.0
    end do

    ! get scalar field values at quadrature points: matmul( 1*snloc, snloc*sngi )
    fields_sngi = 2.0 * face_val_at_quad(kk, sele) * face_val_at_quad(EV, sele)

    ! Perform rhs surface integral
    rhs = shape_rhs( fshape_kk, detwei_bdy * dkdn_contracted * fields_sngi )

    ! block to add to mass matrix diagonal
    mass = shape_shape( fshape_kk, fshape_kk, detwei_bdy )

end subroutine eps_bc

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
    integer                                    :: dim, ngi, nloc, gi, i, j

    nloc = size(du_t,1)
    ngi  = size(du_t,2)
    dim  = size(du_t,3)

    do gi=1, ngi

       S = matmul( nu, du_t(:,gi,:) )
       T = S
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
