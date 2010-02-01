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
!  use rotated_boundary_conditions

  implicit none

  private

  ! These variables are the parameters requried by k-epsilon. 
  ! They are all private to prevent tampering
  ! and saved so that we don't need to call keps_init every time we need them.

  integer, save            :: nNodes
  ! turbulent kinetic energy, epsilon, lengthscale, eddy viscosity, production
  type(scalar_field), save :: TKE_old, eps, ll, nut, P
  ! deltas for subroutine calc_tke and calc_eps
  type(scalar_field), save :: delta_KK, delta_eps
  real, save               :: eps_min = 1.e-12, k_min
  ! Experimental constants from Diamond
  real, save               :: c_mu, c_eps_1, c_eps_2, sigma_eps, sigma_k, theta

  ! Gradient tensor
  type(tensor_field), save             :: DU_DX
  ! these are the fields and variables for the surface values
  type(scalar_field), save             :: surface_values
  ! these are used to populate the bcs
  type(scalar_field), save             :: surface_KK_values
  ! for the eps BC
  integer, save                        :: NNodes_sur
  integer, dimension(:), pointer, save :: surface_nodes
  integer, dimension(:), pointer, save :: surface_element_list
  logical, save                        :: initialised, calculate_bcs

  ! The following are the public subroutines
  public :: keps_init, keps_cleanup, calc_tke, eddyvisc, calc_eps, keps_init_surfaces

  ! General plan is:
  !  - Init in populate_state
  !  - If solve is about to do TKE, call calc_tke (which calculates P and sets source/absorption for solve)
  !  - If solve is about to do epsilon, call calc_eps (which fixes TKE surfaces, set source/absorption for solve)
  !  - After calc_eps solve, recalculate the viscosity and the lengthscale
  !  - When done, clean-up
  !
  ! TurbulentKineticEnergy and Epsilon need to have higher priority then other prognostic 
  ! fields such as temperature, velocity and salinity for this to work

contains

!----------
! init does the following:
!    - check we have the right fields (if not abort)
!    - initialise  parameters based on options
!    - allocate space for optional fields, which are module level variables (to save passing them around)
!----------
subroutine keps_init(state)

    type(state_type), intent(inout) :: state
    type(scalar_field), pointer     :: scalarField, s_cur,s_old, dtwo_kk
    type(vector_field), pointer     :: vectorField
    type(tensor_field), pointer     :: tensorField
    integer                         :: stat

    ewrite(1,*)'Now in k_epsilon turbulence model - JB'

    ! check we have the minimal amount of fields
    scalarField  => extract_scalar_field(state, "TurbulentKineticEnergy",stat)
    if(stat/=0) FLAbort("Need TurbulentKineticEnergy field")
    tensorField => extract_tensor_field(state, "TurbulentKineticEnergyDiffusivity",stat)
    if(stat/=0) FLAbort("Need TurbulentKineticEnergyDiffusivity field")
    scalarField => extract_scalar_field(state, "Epsilon",stat)
    if(stat/=0) FLAbort("Need Epsilon field")
    tensorField => extract_tensor_field(state, "EpsilonDiffusivity",stat)
    if(stat/=0) FLAbort("Need EpsilonDiffusivity field")
    scalarField  => extract_scalar_field(state, "TurbulentKineticEnergySource",stat)
    if(stat/=0) FLAbort("Need source for TurbulentKineticEnergySource field")
    scalarField  => extract_scalar_field(state, "EpsilonSource",stat)
    if(stat/=0) FLAbort("Need source for EpsilonSource field")
    scalarField  => extract_scalar_field(state, "TurbulentKineticEnergyAbsorption",stat)
    if(stat/=0) FLAbort("Need source for TurbulentKineticEnergyAbsorption field")
    scalarField  => extract_scalar_field(state, "EpsilonAbsorption",stat)
    if(stat/=0) FLAbort("Need source for EpsilonAbsorption field")
    scalarField  => extract_scalar_field(state, "LengthScale",stat)
    if(stat/=0) FLAbort("Need source for LengthScale field")
    tensorField => extract_tensor_field(state, "Viscosity",stat)
    if(stat/=0) FLAbort("Need Viscosity")
    tensorField  => extract_tensor_field(state, "EddyViscosity",stat)
    if(stat/=0) FLAbort("Need EddyViscosity")
    vectorField => extract_vector_field(state, "Velocity",stat)
    if(stat/=0) FLAbort("Need Velocity")

    ! Do we have BCs specified for both fields? Better be!
    calculate_bcs = &
    have_option("/material_phase[0]/subgridscale_parameterisations/k_epsilon/scalar_field::TurbulentKineticEnergy/prognostic/boundary_conditions[0]") &
    .and. have_option("/material_phase[0]/subgridscale_parameterisations/k_epsilon/scalar_field::Epsilon/prognostic/boundary_conditions[0]")
    if (calculate_bcs) then
        FLAbort("Need to specify Dirichlet boundary conditions for TKE and Epsilon")
    end if

    ! populate some useful variables including the 5 model constants
    nNodes = node_count(scalarField)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k_epsilon/scalar_field::TurbulentKineticEnergy/prognostic/temporal_discretisation/theta', theta, default = 1.0)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/C_mu', C_mu, default = 0.09)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/C_eps_1', c_eps_1, default = 1.44)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/C_eps_2', c_eps_2, default = 1.92)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/sigma_k', sigma_k, default = 1.0)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/sigma_eps', sigma_eps, default = 1.3)

    ! allocate some space for the fields we need for calculations, but are optional in the model
    ! we're going to allocate these on the velocity mesh as we need one of these...
    call allocate(ll,         vectorField%mesh, "LengthScale")
    call allocate(P,          vectorField%mesh, "ShearProduction")
    call allocate(nut,        vectorField%mesh, "EddyViscosity")
    call allocate(eps,        vectorField%mesh, "TKE_Dissipation")
!    call allocate(Fwall,      vectorField%mesh, "WallFunction")
    call allocate(DU_DX,      tensorField%mesh, "DUDX")
    call allocate(TKE_old,    vectorField%mesh, "Old_TKE")
    call allocate(delta_KK,   vectorField%mesh, "deltaKK")
    call allocate(dtwo_kk,    vectorField%mesh, "dtwoKK")
    call allocate(delta_eps,  vectorField%mesh, "deltaeps")

    ! Hard-coded by Jon Hill. Must be something small and nonzero.
    k_min = 7.6e-6
    eps_min = 1.e-12

    ewrite(1,*) "Parameters"
    ewrite(1,*) "--------------------------------------------"
    ewrite(1,*) "c_mu: ",c_mu
    ewrite(1,*) "c_eps_1: ",c_eps_1
    ewrite(1,*) "c_eps_2: ",c_eps_2
    ewrite(1,*) "sigma_k: ",sigma_k
    ewrite(1,*) "sigma_eps: ",sigma_eps
    ewrite(1,*) "--------------------------------------------"
    
    ! initialise 2 fields (k, epsilon) with minimum values
    s_cur => extract_scalar_field(state, "TurbulentKineticEnergy")
    s_old => extract_scalar_field(state, "OldTurbulentKineticEnergy")
    call set(s_cur,k_min)
    call set(s_old,k_min)
    s_cur => extract_scalar_field(state, "Epsilon")
    s_old => extract_scalar_field(state, "OldEpsilon")
    call set(s_cur,eps_min)
    call set(s_old,eps_min)

    ! initialise other fields
    call set(P, 0.0)
    s_cur => extract_scalar_field(state, "EddyViscosity")
    s_old => extract_scalar_field(state, "OldEddyViscosity")    
    call set(s_cur,1e-6)
    call set(s_old,1e-6)
    s_cur => extract_scalar_field(state, "LengthScale")
    s_old => extract_scalar_field(state, "OldLengthScale")    
    call set(s_cur,k_min**1.5/eps_min)
    call set(s_old,k_min**1.5/eps_min)
    call set(ll,k_min**1.5/eps_min)

    ! initialise surface
    if (calculate_bcs) then
        initialised = .false.
        call keps_init_surfaces(state)
    end if
end subroutine keps_init


!----------
! Update TKE
!----------
subroutine calc_tke(state)

    type(state_type), intent(inout)   :: state
    type(scalar_field), pointer       :: source, absorption, kk, eps, s_Old, scalarField, s_current, &
                                         & scalar_surface, utemp
    type(vector_field), pointer       :: v_Old, v_current, positions
    type(tensor_field), pointer       :: kk_diff
    real                              :: diss
    integer                           :: i, ii, j, stat
    character(len=FIELD_NAME_LEN)     :: bc_type
    type(scalar_field), dimension(1)  :: grtemp
    type(vector_field)                :: velocity
    logical, allocatable, dimension(:):: derivs

!    call gls_buoyancy(state)
!subroutine gls_buoyancy(state)
  
    ! grab variables required from state - already checked in init, so no need to check here
    positions => extract_vector_field(state, "Coordinate")
    v_Old => extract_vector_field(state, "OldVelocity")
    v_current => extract_vector_field(state, "Velocity")
    source => extract_scalar_field(state, "TurbulentKineticEnergySource")
    absorption  => extract_scalar_field(state, "TurbulentKineticEnergyAbsorption")
    kk => extract_scalar_field(state, "TurbulentKineticEnergy")
    kk_diff => extract_tensor_field(state, "TurbulentKineticEnergyDiffusivity")
eps => extract_scalar_field(state, "Epsilon")
    ! Temporary field for calculating gradient field:
    utemp => extract_scalar_field(state, "Velocity")

    ! Update the velocity field:
    call set(velocity, v_current)
    call scale(velocity, theta)
    call addto(velocity, v_Old, 1.0-theta)
    
    ! now allocate our temp fields
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
                call addto(P, ii, node_val(nut, ii) * DU_DX%val(i,j,ii)+DU_DX%val(j,i,ii) &
                * DU_DX%val(i,j,ii)  )
            end do
        end do
        diss = node_val(eps, ii)
        call set(source, ii, P%val(ii) )
        call set(absorption, ii, diss / node_val(delta_KK, ii))
    end do

    call deallocate(velocity)
    deallocate(derivs)

!end subroutine gls_buoyancy

    ! work out the deltas:
    !delta_KK = theta*kk + (1-theta)*OldKK
    s_Old => extract_scalar_field(state, "OldTurbulentKineticEnergy")
    call set(delta_KK, kk)
    call scale(delta_KK, theta)
    call addto(delta_KK, s_Old, 1.0-theta)
    s_Old => extract_scalar_field(state, "OldEddyViscosity")
    s_current => extract_scalar_field(state, "EddyViscosity")
    call set(nut, s_current)
    call scale(nut, theta)
    call addto(nut, s_Old, 1.0-theta)

    ! set diffusivity for KK.
    call set(kk_diff,kk_diff%dim,kk_diff%dim,nut,scale=1./sigma_k)
    ! add in background (need to grab from Diamond, rather than hardcode)
    do ii=1,nNodes
        call addto(kk_diff,kk_diff%dim,kk_diff%dim,ii,1.0e-6)
    end do

    ! boundary conditions
    if (calculate_bcs) then
        call get_option("/material_phase[0]/subgridscale_parameterisations/k_epsilon/scalar_field::TurbulentKineticEnergy/prognostic/boundary_conditions[0]/type", bc_type)
        ! puts the BC boundary values in surface_values (module level variable):
        call tke_bc(state,bc_type)
        ! map these onto the actual BC in kk
        scalar_surface => extract_surface_field(KK, 'tke_boundary', "value")
        call remap_field(surface_values, scalar_surface)
    end if
    
    ! finally, we need a copy of the old TKE for eps, so grab it before we solve
    call set(TKE_old,delta_KK)

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

end subroutine calc_tke


!----------
! Calculate Epsilon
!----------
subroutine calc_eps(state)

    type(state_type), intent(inout)  :: state 
    
    type(scalar_field), pointer      :: source, absorption, kk, kkOld, eps, s_Old, scalarField, s_current, scalar_surface
    type(tensor_field), pointer      :: eps_visc, eps_diff
    real                             :: prod, diss, EpsOverTke
    character(len=FIELD_NAME_LEN)    :: bc_type
    integer                          :: i, stat

    source => extract_scalar_field(state, "EpsilonSource")
    absorption  => extract_scalar_field(state, "EpsilonAbsorption")
    kk => extract_scalar_field(state, "TurbulentKineticEnergy")
    kkOld => extract_scalar_field(state, "OldTurbulentKineticEnergy")
    eps => extract_scalar_field(state, "Epsilon")
    eps_diff => extract_tensor_field(state, "EpsilonDiffusivity")

    ! work out the delta_KK = theta*kk + (1-theta)*OldKK
    call set(delta_KK, kk)
    call scale(delta_KK, theta)
    call addto(delta_KK, kkOld, 1.0-theta)
    s_Old => extract_scalar_field(state, "OldLengthScale")
    s_current => extract_scalar_field(state, "LengthScale")
    call set(ll, s_current)
    call scale(ll, theta)
    call addto(ll, s_Old, 1.0-theta)
    s_Old => extract_scalar_field(state, "OldShearProduction")
    s_current => extract_scalar_field(state, "ShearProduction")
    call set(P, s_current)
    call scale(P, theta)
    call addto(P, s_Old, 1.0-theta)
    s_Old => extract_scalar_field(state, "OldEddyViscosity")
    s_current => extract_scalar_field(state, "EddyViscosity")
    call set(nut, s_current)
    call scale(nut, theta)
    call addto(nut, s_Old, 1.0-theta)

    !
    if (calculate_bcs) then ! BC for epsilon is d^2k/dn^2
        call get_boundary_condition(kk, 'tke_boundary', type=bc_type)
     
        ! We need to just add a dirichlet BC to the TKE first, so we get the right value for the rest
        ! of this calculation on the surfaces
        ! only do this if we're using a neumann BC on the solve.
        if (bc_type == "neumann") then
    
     
            ! call the bc code, but specify we want dirichlet
            call tke_bc(state, 'dirichlet')
     
            ! copy the values onto the mesh using the global node id
            do i=1,NNodes_sur
                call set(delta_KK,surface_nodes(i),node_val(surface_values,i))
            end do
        end if
    end if

    ! clip at k_min
    do i=1,nNodes
        call set(delta_KK,i, max(node_val(delta_KK,i),k_min))
        ! do we need to set KK to min too? Not sure, so doing it anyway.
        call set(KK,i, max(node_val(KK,i),k_min))
        !call set(TKE_old,i, max(node_val(TKE_old,i),k_min))
    end do

    ! re-construct eps at "old" timestep
    do i=1,nNodes
        ! This is based on delta_KK *before* the TKE solve
        ! Hence the reason we don't need to work out delta_eps - JON: IS THIS STILL THE CASE?
        call set(eps, i, (node_val(TKE_old,i))**1.5)
        call scale(eps, 1./ll%val(i))
    end do
write(*,*) "In Epsilon calculation"
    ! compute RHS
    do i=1,nNodes
        ! compute production terms in eps-equation
        ! P has already had the correct non-linear iteration set
        EpsOverTke  = node_val(eps,i)/node_val(TKE_old,i)
        prod        = c_eps_1*EpsOverTke*node_val(P,i)
        diss        = c_eps_2*EpsOverTke*node_val(eps,i)
        write(*,*) prod,diss,EpsOverTke,node_val(eps,i),node_val(TKE_old,i),node_val(P,i),node_val(eps,i)
        ! removed buoyancy if statement (see GLS)
        call set(source, i, prod)
        call set(absorption, i, diss/node_val(eps,i))
    end do

    ! Set eddy viscosity scaled by constant (sigma_eps) for Epsilon
    call set(eps_visc,eps_visc%dim,eps_visc%dim,nut,scale=1./sigma_eps)
    ! Add in background (need to grab from Diamond rather than hardcode)
    do i=1,nNodes
        call addto(eps_visc,eps_visc%dim,eps_visc%dim,i,1.0e-6)
    end do

    ! boundary conditions
    if (calculate_bcs) then
        call get_option( &
        "/material_phase[0]/subgridscale_parameterisations/k_epsilon/scalar_field::TurbulentKineticEnergy/prognostic/boundary_conditions[0]/type",&
        bc_type)
        ! puts the BC boundary values in surface_values (module level variable):
        call eps_bc(state,bc_type)
        ! map these onto the actual BC in kk
        scalar_surface => extract_surface_field(KK, 'eps_boundary', "value")
        call remap_field(surface_values, scalar_surface)
    end if

    ! Epsilon is now ready for solving (see Fluids.F90)
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

end subroutine calc_eps


!----------
! eddyvisc calculates the lengthscale, and then the viscosity.
! These are placed in the GLS fields ready for other tracer fields to use
! Viscosity is placed in the velocity viscosity
!----------
subroutine eddyvisc(state)

    type(state_type), intent(inout)  :: state
    type(scalar_field), pointer      :: KK, eps, epsOld, kkOld
    type(tensor_field), pointer      :: eddy_visc, viscosity, background_visc
    real                             :: x
    character(len=FIELD_NAME_LEN)    :: bc_type
    integer                          :: i, stat
    real                             :: tke

    KK => extract_scalar_field(state, "TurbulentKineticEnergy")
    KKOld => extract_scalar_field(state, "OldTurbulentKineticEnergy")
    eps => extract_scalar_field(state, "Epsilon")
    epsOld => extract_scalar_field(state, "OldEpsilon")
    eddy_visc => extract_tensor_field(state, "EddyViscosity",stat)
    viscosity => extract_tensor_field(state, "Viscosity",stat) 

    call set(delta_eps, eps)
    call scale(delta_eps, theta)
    call addto(delta_eps, epsOld, 1.0-theta)
    call set(delta_KK, kk)
    call scale(delta_KK, theta)
    call addto(delta_KK, kkOld, 1.0-theta)

    ! Hard-coded in the GLS option "fix_surface_values": I want to do this by default (I think).
    call get_boundary_condition(eps, 'eps_boundary', type=bc_type)
    ! We need to just add a dirichlet BC to the TKE first, so we get the right value for the rest
    ! of this calculation on the surfaces
    ! only do this if we're using a neumann BC on the solve.
    if (bc_type == "neumann") then
        ! call the bc code, but specify we want dirichlet
        call eps_bc(state, bc_type)
        ! copy the values onto the mesh using the global node id
        do i=1,NNodes_sur
            call set(delta_eps,surface_nodes(i),node_val(surface_values,i))
        end do   
    end if

write(*,*) "In diffusivity" 
    do i=1,nNodes
        tke = delta_KK%val(i)
        call set(eps, i, tke**1.5/node_val(ll,i))
        ! clip at eps_min
        call set(eps, i, max(node_val(eps,i),eps_min))
        ! compute dissipative scale - MAY NEED CHANGING
        call set(ll, i, tke**1.5/node_val(eps,i))
    end do

    ! calculate viscosity for next step and for use in other fields
    do i=1,nNodes
        call set( nut, i, C_mu * tke**2. / eps%val(i) )
        write(*,*) node_val(delta_KK,i),node_val(delta_eps,i),node_val(eps,i),node_val(ll,i), x, &
        node_val(nut,i)
    end do

    !set the eddy_diffusivity and viscosity tensors for use by other fields
    call zero(eddy_visc) ! zero it first as we're using an addto below
    call set(eddy_visc,eddy_visc%dim,eddy_visc%dim,nut) 

    background_visc => extract_tensor_field(state, "BackgroundViscosity",stat)
    if(stat == 0) then
        call addto(eddy_visc,background_visc)    
    endif   

    ! set the velocity viscosity to the computed eddy viscosity
    call set(viscosity,viscosity%dim,viscosity%dim,nut) 
    ! Add background
    !TODO: grab from Diamond rather than hardcode
    do i = 1,nNodes
        call addto(viscosity,viscosity%dim,viscosity%dim,i,1.0e-6) 
    end do

    ! Set output on optional fields - if the field exists, stick something in it
    ! We only need to do this to those fields that we haven't pulled from state, but
    ! allocated ourselves
    call output_fields(state)

end subroutine eddyvisc


!----------
! cleanup does...have a guess...go on.
!----------
subroutine keps_cleanup()

    ! deallocate all our variables 
    call deallocate(ll)
    call deallocate(P)
    call deallocate(nut)
    call deallocate(eps)
!    call deallocate(Fwall)
    call deallocate(DU_DX)
    call deallocate(TKE_old)
    call deallocate(delta_KK)
    call deallocate(delta_eps)
    if (calculate_bcs) then
        call deallocate(surface_values)
        call deallocate(surface_kk_values)
    end if

end subroutine keps_cleanup


!---------
! initialise the surface meshes used for the BCS
! Called at startup and after an adapt
!----------
subroutine keps_init_surfaces(state)
    type(state_type), intent(in)     :: state
    character(len=OPTION_PATH_LEN)   :: input_mesh_name
    type(scalar_field), pointer      :: kk, bc_field
    character(len=FIELD_NAME_LEN)    :: bc_type
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
    call deallocate(boundary_mesh)

    initialised = .true.
end subroutine keps_init_surfaces



!------------------------------------------------------------------!
!------------------------------------------------------------------!
!                                                                  !
!                       Private subroutines                        !
!                                                                  !
!------------------------------------------------------------------!
!------------------------------------------------------------------!



!----------
! tke_bc calculates the boundary conditions on the TKE (kk) field
! Boundary can be either Dirichlet or Neumann.
!----------
subroutine tke_bc(state, bc_type)

    type(state_type), intent(in)     :: state
    character(len=*), intent(in)     :: bc_type
    integer                          :: i

    select case(bc_type)
    case("dirichlet")
    do i=1,NNodes_sur
        call set(surface_values,i,0.0)
    end do
    case default
        FLAbort('Unknown surface BC for TKE')
    end select

end subroutine tke_bc

!----------
! eps_bc calculates the boundary conditions on the eps (eps) field
! Hard-coded to be eps = d^2k/dn^2 for now.
!----------
subroutine eps_bc(state, bc_type)

    type(state_type), intent(in)     :: state
    type(vector_field), pointer:: normal
    integer                          :: i
    character(len=FIELD_NAME_LEN)    :: bc_type
    type(scalar_field), pointer      :: kk, kkOld, dtwo_kk
    type(vector_field), pointer      :: u
    character(len=OPTION_PATH_LEN)   :: bc_option_path



    ! grab hold of some essential fields
    kk      => extract_scalar_field(state, "TurbulentKineticEnergy")    
    kkOld   => extract_scalar_field(state, "OldTurbulentKineticEnergy")
    dtwo_kk => extract_scalar_field(state, "SecondDerivativeTKE")
    u       => extract_vector_field(state, "Velocity")
    ! work out the delta_KK = theta*kk + (1-theta)*OldKK
    call set(delta_KK, kk)
    call scale(delta_KK, theta)
    call addto(delta_KK, kkOld, 1.0-theta)
    

    ! HOW DO I GET THE NORMAL TO THE BOUNDARY?
    ! Procedure:
        ! 1. Find boundary elements using something from Fields_Base.
        ! 2. Looping over elements, call:
        ! transform_facet_to_physical_full(X, face, detwei_f, normal) 
        ! (this computes the normal to a given facet) from Module Transform_Elements.
        ! N.B. I want the normals at NODES not elements.
        ! 4. How do I input an arbitrary vector into the differentiate_field function?
        ! 
        ! call remap_field_to_surface(dtwo_kk, surface_kk_values, surface_element_list)
        ! Take second derivative of surface_kk_values and set value equal to that
        ! value = node_val(surface_kk_values,i)
        ! call set(surface_values,i,value)
        ! end do

    ! OR may be able to adapt this routine from Rotated_BCs. Although this could be
    ! overcomplicating things!
    ! Rotated_BCs is not available for scalar fields. Can I do a hack to make it
    ! available? 



    do i=1, get_boundary_condition_count(u)
       call get_boundary_condition(u, i, type=bc_type, &
            surface_node_list=surface_nodes, &
            option_path=bc_option_path)

       if (bc_type=="dirichlet" .and. &
            have_option(trim(bc_option_path)//"/type[0]/k_epsilon_boundary")) then
            normal   => extract_surface_field(u, i, "normal")
       end if
    end do


!    do i=1,NNodes_sur
!        call set(surface_values, i, dtwo_kk%val(i) )
!        call differentiate_field(delta_kk, i, normal, derivatives, dtwo_kk)
!    end do


end subroutine eps_bc


!---------
! Output the optional fields if they exist in state
!---------
subroutine output_fields(state)

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

    scalarField => extract_scalar_field(state, "DissipationEpsilon", stat)
    if(stat == 0) then
        call set(scalarField,eps) 
    end if  

    scalarField => extract_scalar_field(state, "Viscosity", stat)
    if(stat == 0) then
       call set(scalarField,nut)  
    end if

end subroutine output_fields

end module k_epsilon

