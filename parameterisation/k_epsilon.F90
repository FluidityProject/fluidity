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
  use allsorts
  use global_parameters, only:   OPTION_PATH_LEN
  use state_fields_module
  use boundary_conditions
  use fields_manipulation
  use FLDebug

  implicit none

  private

  ! These variables are the parameters requried by k-epsilon. 
  ! They are all private to prevent tampering
  ! and save'd so that we don't need to call keps_init every time some GLS-y
  ! happens.

  real, save               :: sigma_eps, sigma_k, kappa
  integer, save            :: nNodes
  ! Used in production terms
  type(scalar_field), save :: TKE_old, ll
  type(scalar_field), save :: eps, Fwall
  type(scalar_field), save :: nut, P
  ! eddy viscosity, production
  type(scalar_field), save :: delta_KK, delta_eps
  real, save               :: eps_min = 1.e-12, k_min
  real, save               :: c_mu, c_eps_1, c_eps_2, theta

  ! Grad tensors
  type(tensor_field), save             :: DU_DX

  ! these are the fields and variables for the surface values
  type(scalar_field), save             :: surface_values
  ! these are used to populate the bcs
  type(scalar_field), save             :: surface_KK_values
  ! for the eps BC
  integer, save                        :: NNodes_sur
  integer, dimension(:), pointer, save :: surface_nodes
  integer, dimension(:), pointer, save :: surface_element_list
  logical, save                        :: initialised, calculate_bcs, fix_surface_values
  real, allocatable, dimension(:)      :: dzb
  
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

    type(scalar_field), pointer    :: scalarField
    type(vector_field), pointer    :: vectorField
    type(tensor_field), pointer    :: tensorField
    real                           :: N
    type(scalar_field), pointer    :: distanceToBottom, distanceToTop 
    integer                        :: i, stat
    character(len=FIELD_NAME_LEN)  :: option
    type(scalar_field), pointer    :: s_cur,s_old

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
    scalarField  => extract_scalar_field(state, "Viscosity",stat)
    if(stat/=0) FLAbort("Need source for Viscosity field")
    scalarField  => extract_scalar_field(state, "LengthScale",stat)
    if(stat/=0) FLAbort("Need source for LengthScale field")
    tensorField => extract_tensor_field(state, "Viscosity",stat)
    if(stat/=0) FLAbort("Need viscosity")
    tensorField  => extract_tensor_field(state, "EddyViscosity",stat)
    if(stat/=0) FLAbort("Need EddyViscosity")
    vectorField => extract_vector_field(state, "Velocity",stat)
    if(stat/=0) FLAbort("Need velocity")
   
    ! populate some useful variables including the 5 model constants
    nNodes = node_count(scalarField)
    call get_option('/material_phase[0]/subgridscale_parameterisations/k-epsilon/TurbulentKineticEnergy/prognostic/temporal_discretisation/theta', theta, default = 1.0)
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
    call allocate(Fwall,      vectorField%mesh, "WallFunction")
    call allocate(DU_DX,      tensorField%mesh, "DUDX")
    call allocate(TKE_old,    vectorField%mesh, "Old_TKE")
    call allocate(delta_KK,   vectorField%mesh, "deltaKK")
    call allocate(delta_eps,  vectorField%mesh, "deltaeps")

    k_min = 7.6e-6
    eps_min = 1.e-12
    call set( Fwall, 1.0 )

    ewrite(1,*) "Parameters"
    ewrite(1,*) "--------------------------------------------"
    ewrite(1,*) "c_mu: ",c_mu
    ewrite(1,*) "c_eps_1: ",c_eps_1
    ewrite(1,*) "c_eps_2: ",c_eps_2
    ewrite(1,*) "sigma_k: ",sigma_k
    ewrite(1,*) "sigma_eps: ",sigma_eps
    ewrite(1,*) "--------------------------------------------"
    
    ! initialise 2 k-epsilon fields with minimum values
    s_cur => extract_scalar_field(state, "TurbulentKineticEnergy")
    s_old => extract_scalar_field(state, "OldTurbulentKineticEnergy")
    call set(s_cur,k_min)
    call set(s_old,k_min)
    s_cur => extract_scalar_field(state, "Epsilon")
    s_old => extract_scalar_field(state, "OldEpsilon")
    call set(s_cur,eps_min)
    call set(s_old,eps_min)

    ! init other fields
    call set(P, 0.0)
    s_cur => extract_scalar_field(state, "EddyViscosity")
    s_old => extract_scalar_field(state, "OldEddyViscosity")    
    call set(s_cur,1e-6)
    call set(s_old,1e-6)
    s_cur => extract_scalar_field(state, "DissipationEpsilon")
    s_old => extract_scalar_field(state, "OldDissipationEpsilon")    
    call set(s_cur,eps_min)
    call set(s_old,eps_min)
    call set(eps,eps_min)    
    s_cur => extract_scalar_field(state, "LengthScale")
    s_old => extract_scalar_field(state, "OldLengthScale")    
    call set(s_cur,k_min**1.5/eps_min)
    call set(s_old,k_min**1.5/eps_min)
    call set(ll,k_min**1.5/eps_min)

    ! intilise surface - only if we need to though
    if (calculate_bcs) then
        initialised = .false.
        call keps_init_surfaces(state)
        call calculate_dz(state)
    end if


    ! we're all done!
end subroutine keps_init

!----------
! Update TKE
!----------
subroutine calc_tke(state)

    type(state_type), intent(inout)  :: state 
    
    type(scalar_field), pointer      :: source, absorption, kk, eps, s_Old, scalarField, s_current
    type(tensor_field), pointer      :: kk_diff
    real                             :: prod, diss
    integer                          :: i, ii, j, stat
    character(len=FIELD_NAME_LEN)    :: bc_type
    type(scalar_field), pointer      :: scalar_surface, utemp
    type(scalar_field), dimension(1)      :: grtemp
    type(vector_field), pointer      :: positions, gradient
!   Formerly variables in gls_buoyancy:
    type(vector_field)               :: velocity
    type(vector_field), pointer      :: v_Old, v_current
    logical, allocatable, dimension(:) :: derivs

    
!    call gls_buoyancy(state)

!subroutine gls_buoyancy(state)
  
    ! grab variables required from state - already checked in init, so no need to check here
    positions => extract_vector_field(state, "Coordinate")
    v_Old => extract_vector_field(state, "OldVelocity")
    v_current => extract_vector_field(state, "Velocity")
    source => extract_scalar_field(state, "TurbulentKineticEnergySource")
    absorption  => extract_scalar_field(state, "TurbulentKineticEnergyAbsorption")
    
    call set(velocity, v_current)
    call scale(velocity, theta)
    call addto(velocity, v_Old, 1.0-theta)
    
    ! now allocate our temp fields
    allocate( derivs(velocity%dim))
    derivs = .false.

    do i = 1, DU_DX%dim
        utemp => extract_scalar_field(state, "Velocity")
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
                ! tensor to be multiplied by eddy viscosity
                call addto(P, ii, node_val(nut, ii) * DU_DX%val(i,j,ii)+DU_DX%val(j,i,ii) &
                * DU_DX%val(i,j,ii)  )
                ! Sum components of tensor to calculate production term
            end do
        end do
        diss = node_val(eps, ii)
        call set(source, ii, P%val(ii) )
        call set(absorption, ii, diss / node_val(delta_KK, ii))
    end do

    call deallocate(velocity)
    deallocate(derivs)

!end subroutine gls_buoyancy

    kk => extract_scalar_field(state, "TurbulentKineticEnergy")
    kk_diff => extract_tensor_field(state, "TurbulentKineticEnergyDiffusivity")


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
    s_Old => extract_scalar_field(state, "OldDissipationEpsilon")
    s_current => extract_scalar_field(state, "DissipationEpsilon")
    call set(eps, s_current)
    call scale(eps, theta)
    call addto(eps, s_Old, 1.0-theta)


    ! set diffusivity for KK
    call set(kk_diff,kk_diff%dim,kk_diff%dim,nut,scale=1./sigma_k)
    ! add in background (need to grab from Diamond, rather than hardcode)
    do ii=1,nNodes
        call addto(kk_diff,kk_diff%dim,kk_diff%dim,ii,1.0e-6)
    end do

    ! boundary conditions
    if (calculate_bcs) then
        call get_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/calculate_boundaries/", bc_type)
        call tke_bc(state,bc_type)
        ! above puts the BC boundary values in top_surface_values and bottom_surface_values module level variables
        ! map these onto the actual BCs in kk
!        scalar_surface => extract_surface_field(KK, 'tke_boundary', "value")
!        call remap_field(bottom_surface_values, scalar_surface)
!        scalar_surface => extract_surface_field(eps, 'td_boundary', "value")
!        call remap_field(top_surface_values, scalar_surface)
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
    
    type(scalar_field), pointer      :: source, absorption, kk, kkOld, eps, s_Old, scalarField, s_current
    type(tensor_field), pointer      :: eps_visc, eps_diff
    real                             :: prod,diss,EpsOverTke
    character(len=FIELD_NAME_LEN)    :: bc_type
    integer                          :: i, stat
    type(scalar_field), pointer      :: scalar_surface

    source => extract_scalar_field(state, "EpsilonSource")
    absorption  => extract_scalar_field(state, "EpsilonAbsorption")
    kk => extract_scalar_field(state, "TurbulentKineticEnergy")
    kkOld => extract_scalar_field(state, "OldTurbulentKineticEnergy")
    eps => extract_scalar_field(state, "Epsilon")
    eps_diff => extract_tensor_field(state, "EpsilonDiffusivity")
!    length => extract_scalar_field(state, "LengthScale")

    ! work out the delta_KK = theta*kk + (1-theta)*OldKK
    call set(delta_KK, kk)
    call scale(delta_KK, theta)
    call addto(delta_KK, kkOld, 1.0-theta)
    s_Old => extract_scalar_field(state, "OldLengthScale")
    s_current => extract_scalar_field(state, "LengthScale")
    call set(ll, s_current)
    call scale(ll, theta)
    call addto(ll, s_Old, 1.0-theta)
    s_Old => extract_scalar_field(state, "OldDissipationEpsilon")
    s_current => extract_scalar_field(state, "DissipationEpsilon")
    call set(eps, s_current)
    call scale(eps, theta)
    call addto(eps, s_Old, 1.0-theta)
    s_Old => extract_scalar_field(state, "OldShearProduction")
    s_current => extract_scalar_field(state, "ShearProduction")
    call set(P, s_current)
    call scale(P, theta)
    call addto(P, s_Old, 1.0-theta)
    s_Old => extract_scalar_field(state, "OldViscosity")
    s_current => extract_scalar_field(state, "Viscosity")
    call set(nut, s_current)
    call scale(nut, theta)
    call addto(nut, s_Old, 1.0-theta)
    

    if (fix_surface_values) then
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
write(*,*) "In Eps"
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

    ! Set eddy viscosity for Eps
    call set(eps_visc,eps_visc%dim,eps_visc%dim,nut,scale=1./sigma_eps)
    ! Add in background (need to grab from Diamond rather than hardcode)
    do i=1,nNodes
        call addto(eps_visc,eps_visc%dim,eps_visc%dim,i,1.0e-6)
    end do

    ! boundary conditions
    if (calculate_bcs) then
        call get_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon/calculate_boundaries/", bc_type)
        call eps_bc(state,bc_type)
        ! above puts the BC boundary values in top_surface_values and bottom_surface_values module level variables
        ! map these onto the actual BCs in eps
!        scalar_surface => extract_surface_field(eps, 'eps_boundary', "value")
!        call remap_field(surface_values, scalar_surface)
    end if


    ! Eps is now ready for solving (see Fluids.F90)
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
! diffusivity fixes the top/bottom boundaries of Eps
! then calulates the lengthscale, and then uses those to calculate the 
! viscosity
! These are placed in the GLS fields ready for other tracer fields to use
! Viscosity is placed in the velocity viscosity
!----------
subroutine eddyvisc(state)

    type(state_type), intent(inout)  :: state 
    
    type(scalar_field), pointer      :: KK, eps, epsOld, kkOld, s_current, s_Old
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


    if (fix_surface_values) then
        call get_boundary_condition(eps, 'eps_boundary', type=bc_type)
    
        ! We need to just add a dirichlet BC to the TKE first, so we get the right value for the rest
        ! of this calculation on the surfaces
        ! only do this if we're using a neumann BC on the solve.
        if (bc_type == "neumann") then
    
            ! call the bc code, but specify we want dirichlet
            call eps_bc(state, 'dirichlet')
    
            ! copy the values onto the mesh using the global node id
            do i=1,NNodes_sur
                call set(delta_eps,surface_nodes(i),node_val(surface_values,i))
            end do   
        end if
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

        !write(*,*) node_val(delta_KK,i),node_val(delta_eps,i),node_val(eps,i),node_val(ll,i), x, &
        ! & node_val(nut,i)

    end do

    !set the eddy_diffusivity and viscosity tensors for use by other fields
    call zero(eddy_visc) ! zero it first as we're using an addto below
    call set(eddy_visc,eddy_visc%dim,eddy_visc%dim,nut) 

    background_visc => extract_tensor_field(state, "BackgroundViscosity",stat)
    if(stat == 0) then
        call addto(eddy_visc,background_visc)    
    endif   

    ! add the viscosity to the velocity viscosity
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
    call deallocate(Fwall)
    call deallocate(DU_DX)
    call deallocate(TKE_old)
    call deallocate(delta_KK)
    call deallocate(delta_eps)
    if (calculate_bcs) then
        call deallocate(surface_values)
        call deallocate(surface_kk_values)
        deallocate(dzb)
    end if

end subroutine keps_cleanup


!---------
! initialise the surface meshes used for the BCS
! Called at startup and after an adapt
!----------
subroutine keps_init_surfaces(state)
    type(state_type), intent(in)     :: state  

    type(scalar_field), pointer      :: distanceToTop, distanceToBottom
    character(len=OPTION_PATH_LEN)   :: input_mesh_name
    type(scalar_field), pointer      :: kk
    type(mesh_type)                  :: input_mesh, ocean_mesh
   
    ! grab hold of some essential field
    distanceToTop => extract_scalar_field(state, "DistanceToTop")
    distanceToBottom => extract_scalar_field(state, "DistanceToBottom")
    kk => extract_scalar_field(state, "TurbulentKineticEnergy")


    ! if we're already initialised, then deallocate surface fields to make space for new ones
    if (initialised) then
        call deallocate(surface_values)
        deallocate(dzb)
    end if

    ! create a surface mesh to place values onto. This is for the top surface
    call get_option(trim(kk%option_path)//'/prognostic/mesh/name', input_mesh_name)
    input_mesh = extract_mesh(state, input_mesh_name);
!    call get_boundary_condition(distanceToTop, name='top', surface_element_list=surface_element_list)
!    call create_surface_mesh(ocean_mesh, surface_nodes, input_mesh, surface_element_list, 'OceanSurface')
!    NNodes_sur = node_count(ocean_mesh)
!    call allocate(surface_values, ocean_mesh, name="surface")
!    call allocate(surface_kk_values,ocean_mesh, name="surface_tke")
!    call deallocate(ocean_mesh)
    ! bottom
!    call get_boundary_condition(distanceToBottom, name='bottom', surface_element_list=bottom_surface_element_list)
!    call create_surface_mesh(ocean_mesh, bottom_surface_nodes, input_mesh, bottom_surface_element_list, 'OceanBottom')

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

    type(vector_field), pointer      :: positions 
!    real                             :: gravity_magnitude
    integer                          :: i
!    real, allocatable, dimension(:)  :: z0s, z0b, u_taus_squared, u_taub_squared
 
    ! grab hold of some essential field

    positions => extract_vector_field(state, "Coordinate") 



end subroutine tke_bc

!----------
! eps_bc calculates the boundary conditions on the eps (eps) field
! Boundary can be either Dirichlet or Neumann.
!----------
subroutine eps_bc(state, bc_type)

    type(state_type), intent(in)     :: state  
    character(len=*), intent(in)     :: bc_type

    type(vector_field), pointer      :: positions 
    integer                          :: i
!    real, allocatable, dimension(:)  :: z0s, z0b, u_taus_squared, u_taub_squared
    type(scalar_field), pointer      :: kk, kkOld
    real                             :: value
 
    ! grab hold of some essential field

    positions => extract_vector_field(state, "Coordinate") 
    kk => extract_scalar_field(state, "TurbulentKineticEnergy")    
    kkOld => extract_scalar_field(state, "OldTurbulentKineticEnergy")

    ! work out the delta_KK = theta*kk + (1-theta)*OldKK
    call set(delta_KK, kk)
    call scale(delta_KK, theta)
    call addto(delta_KK, kkOld, 1.0-theta)



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

    scalarField => extract_scalar_field(state, "WallFunction", stat)
    if(stat == 0) then
        call set(scalarField,FWall) 
    end if           

    scalarField => extract_scalar_field(state, "Viscosity", stat)
    if(stat == 0) then
       call set(scalarField,nut)  
    end if


end subroutine output_fields


end module k_epsilon

