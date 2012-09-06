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
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN, timestep, current_time
  use state_fields_module
  use boundary_conditions
  use fields_manipulation
  use surface_integrals
  use fetools
  use vector_tools
  use sparsity_patterns_meshes
  use FLDebug
  use vtk_interfaces
  use solvers

implicit none

  private

  ! locally allocatad fields
  real, save     :: fields_min = 1.0e-11
  logical, save  :: low_Re = .false.                     

  public :: keps_advdif_diagnostics, keps_momentum_diagnostics, keps_bcs, &
       & k_epsilon_check_options, tensor_inner_product

  ! Outline:
  !  - call diagnostics to obtain source terms and calculate eddy viscosity
  !  - after each scalar field solve recalculates the eddy viscosity
  !  - wall functions are added to selected boundaries in keps_bcs and wall_functions

contains

subroutine keps_advdif_diagnostics(state)

  type(state_type), intent(inout) :: state
  
  call keps_damping_functions(state, advdif=.true.)
  call keps_eddyvisc(state, advdif=.true.)
  call keps_diffusion(state)
  call keps_tracer_diffusion(state)
  call keps_calculate_rhs(state)

end subroutine keps_advdif_diagnostics

subroutine keps_momentum_diagnostics(state)

  type(state_type), intent(inout) :: state
  
  call keps_damping_functions(state, advdif=.false.)
  call keps_eddyvisc(state, advdif=.false.)
  call keps_momentum_source(state)

end subroutine keps_momentum_diagnostics

!--------------------------------------------------------------------------------!

subroutine keps_damping_functions(state, advdif)

  type(state_type), intent(in) :: state
  logical, intent(in) :: advdif
  
  type(scalar_field), pointer :: f_1, f_2, f_mu, y, dummydensity, density
  type(scalar_field) :: k, eps
  type(tensor_field), pointer :: bg_visc
  integer :: node, stat
  real :: rhs, Re_T, R_y, fields_max
  character(len=FIELD_NAME_LEN) :: equation_type
  character(len=OPTION_PATH_LEN) :: option_path

  ewrite(1,*) 'in keps_damping_functions'

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-epsilon/'

  f_1 => extract_scalar_field(state, "f_1")
  f_2 => extract_scalar_field(state, "f_2")
  f_mu => extract_scalar_field(state, "f_mu")

  ! initialise low_Re damping functions
  call set(f_1, 1.0)
  call set(f_2, 1.0)
  call set(f_mu, 1.0)

  call get_option(trim(state%option_path)// &
       & "/subgridscale_parameterisations/k-epsilon/max_damping_value", fields_max) 

  ! Low Reynolds damping functions
  ! Check for low reynolds boundary condition and calculate damping functions
  ! Lam-Bremhorst model (Wilcox 1998 - Turbulence modelling for CFD)
  if (low_Re .or. &
       have_option(trim(option_path)//"debugging_options/enable_lowRe_damping")) then
  
     bg_visc => extract_tensor_field(state, "BackgroundViscosity")
     y => extract_scalar_field(state, "DistanceToWall", stat = stat)
     if (stat /= 0) then
        FLAbort("I need the distance to the wall - enable a DistanceToWall field")
     end if

     call time_averaged_value(state, k, 'TurbulentKineticEnergy', advdif, option_path)
     call time_averaged_value(state, eps, 'TurbulentDissipation', advdif, option_path)

     allocate(dummydensity)
     call allocate(dummydensity, f_1%mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
     call set(dummydensity, 1.0)
     dummydensity%option_path = ""

     ! Depending on the equation type, extract the density or set it to some dummy field allocated above
     call get_option(trim(state%option_path)//&
          "/vector_field::Velocity/prognostic/equation[0]/name", equation_type)
     select case(equation_type)
     case("LinearMomentum")
        density=>extract_scalar_field(state, "Density")
     case("Boussinesq")
        density=>dummydensity
     case("Drainage")
        density=>dummydensity
     case default
        ! developer error... out of sync options input and code
        FLAbort("Unknown equation type for velocity")
     end select

     node_loop: do node = 1, node_count(k)

        if (node_val(k,node) <= fields_min .or. &
             & node_val(y,node) <= fields_min .or. &
             & node_val(bg_visc,1,1,node) <= fields_min .or. &
             & node_val(density,node) <= fields_min .or. &
             & node_val(eps,node) <= fields_min) then
           call set(f_mu, node, 0.0)
           call set(f_1, node, 0.0)
           call set(f_2, node, 0.0)

           cycle node_loop
        end if

        Re_T = (node_val(density,node) * node_val(k,node)**2.0) / &
             (node_val(eps,node) * node_val(bg_visc,1,1,node))
        R_y = (node_val(density,node) * node_val(k,node)**0.5 * node_val(y,node)) / &
             node_val(bg_visc,1,1,node)
        
        rhs = (- exp(- 0.0165*R_y) + 1.0)**2.0 * (20.5/Re_T + 1.0)
        if (rhs > 1.0) then
           call set(f_mu, node, 1.0)
           call set(f_1, node, 1.0)
           call set(f_2, node, 1.0)
           
           cycle node_loop
        end if
        call set(f_mu, node, min(rhs,fields_max))

        rhs = (0.05/node_val(f_mu,node))**3.0 + 1.0
        call set(f_1, node, min(rhs,fields_max))

        rhs = -exp(- Re_T**2.0) + 1.0
        call set(f_2, node, min(rhs,fields_max))       
        
     end do node_loop

     call deallocate(k)
     call deallocate(eps)
     call deallocate(dummydensity)
     deallocate(dummydensity)

  end if

end subroutine keps_damping_functions
    
!------------------------------------------------------------------------------!

subroutine keps_calculate_rhs(state)

  type(state_type), intent(inout) :: state

  type(scalar_field), dimension(3) :: src_abs_terms
  type(scalar_field), dimension(2) :: fields
  type(scalar_field), pointer :: src, abs, f_1, f_2, debug
  type(scalar_field) :: src_to_abs
  type(vector_field), pointer :: x, u, g
  type(scalar_field), pointer :: dummydensity, density, buoyancy_density, scalar_eddy_visc
  integer :: i, ele, term, stat
  real :: g_magnitude, c_eps_1, c_eps_2, sigma_p
  logical :: have_buoyancy_turbulence = .true., lump_mass
  character(len=OPTION_PATH_LEN) :: option_path 
  character(len=FIELD_NAME_LEN), dimension(2) :: field_names
  character(len=FIELD_NAME_LEN) :: equation_type, implementation

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-epsilon/'

  if (.not. have_option(trim(option_path))) then 
     return
  end if

  ewrite(1,*) 'In calculate k-epsilon rhs'

  ! get model constants
  call get_option(trim(option_path)//'/C_eps_1', c_eps_1, default = 1.44)
  call get_option(trim(option_path)//'/C_eps_2', c_eps_2, default = 1.92)
  call get_option(trim(option_path)//'/sigma_p', sigma_p, default = 1.0)
  
  ! get field data
  x => extract_vector_field(state, "Coordinate")
  u => extract_vector_field(state, "NonlinearVelocity")
  scalar_eddy_visc => extract_scalar_field(state, "ScalarEddyViscosity")
  f_1 => extract_scalar_field(state, "f_1")
  f_2 => extract_scalar_field(state, "f_2")
  g => extract_vector_field(state, "GravityDirection", stat)
  call get_option('/physical_parameters/gravity/magnitude', g_magnitude, stat)
  if (stat /= 0) then
     have_buoyancy_turbulence = .false.
  end if

  allocate(dummydensity)
  call allocate(dummydensity, X%mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
  call set(dummydensity, 1.0)
  dummydensity%option_path = ""
  
  ! Depending on the equation type, extract the density or set it to some dummy field allocated above
  call get_option(trim(state%option_path)//&
       "/vector_field::Velocity/prognostic/equation[0]/name", equation_type)
  select case(equation_type)
     case("LinearMomentum")
        density=>extract_scalar_field(state, "Density")
     case("Boussinesq")
        density=>dummydensity
     case("Drainage")
        density=>dummydensity
     case default
        ! developer error... out of sync options input and code
        FLAbort("Unknown equation type for velocity")
  end select
  
  if(have_buoyancy_turbulence) then
     buoyancy_density => extract_scalar_field(state, 'VelocityBuoyancyDensity')
  end if

  field_names(1) = 'TurbulentKineticEnergy'
  field_names(2) = 'TurbulentDissipation'

  field_loop: do i = 1, 2
     if (have_option(trim(option_path)//'scalar_field::'// &
          trim(field_names(i))//'/prescribed')) then
        cycle
     end if

     !-----------------------------------------------------------------------------------
     
     ! Setup
     src => extract_scalar_field(state, trim(field_names(i))//"Source")
     abs => extract_scalar_field(state, trim(field_names(i))//"Absorption")

     call time_averaged_value(state, fields(1), trim(field_names(i)), .true., option_path)
     call time_averaged_value(state, fields(2), trim(field_names(3-i)), .true., option_path)

     call allocate(src_abs_terms(1), fields(1)%mesh, name="production_term")
     call allocate(src_abs_terms(2), fields(1)%mesh, name="destruction_term")
     call allocate(src_abs_terms(3), fields(1)%mesh, name="buoyancy_term")
     call zero(src_abs_terms(1)); call zero(src_abs_terms(2)); call zero(src_abs_terms(3))
     call zero(src); call zero(abs)

     !-----------------------------------------------------------------------------------

     ! Assembly loop
     do ele = 1, ele_count(fields(1))
        call assemble_rhs_ele(src_abs_terms, fields(i), fields(3-i), scalar_eddy_visc, u, &
             density, buoyancy_density, have_buoyancy_turbulence, g, g_magnitude, x, &
             c_eps_1, c_eps_2, sigma_p, f_1, f_2, ele, i)
     end do

     ! For non-DG we apply inverse mass globally
     if(continuity(fields(1))>=0) then
        lump_mass = have_option(trim(option_path)//'mass_lumping/lump_mass')
        do term = 1, 3
           call solve_cg_inv_mass(state, src_abs_terms(term), lump_mass, option_path)           
        end do
     end if
     !-----------------------------------------------------------------------------------

     ! Source disabling for debugging purposes
     do term = 1, 3
        if(have_option(trim(option_path)//'debugging_options/disable_'//&
             trim(src_abs_terms(term)%name))) then
           call zero(src_abs_terms(term))
        end if
     end do    
     !-----------------------------------------------------------------------------------

     ! Produce debugging output
     do term = 1, 3
        debug => extract_scalar_field(state, &
          trim(field_names(i))//"_"//trim(src_abs_terms(term)%name), stat)
        if (stat == 0) then
           call set(debug, src_abs_terms(term))
        end if
     end do
     !-----------------------------------------------------------------------------------
     
     ! Implement terms as source or absorbtion
     do term = 1, 3
        call get_option(trim(option_path)//&
             'time_discretisation/source_term_implementation/'//&
             trim(src_abs_terms(term)%name), implementation)
        select case(implementation)
        case("source")
           call addto(src, src_abs_terms(term))
        case("absorbtion")
           call allocate(src_to_abs, fields(1)%mesh, name='SourceToAbsorbtion')
           call set(src_to_abs, fields(1))
           where (src_to_abs%val >= fields_min)
              src_to_abs%val=1./src_to_abs%val
           elsewhere
              src_to_abs%val=1./fields_min
           end where
           call scale(src_abs_terms(term), src_to_abs)
           call addto(abs, src_abs_terms(term), -1.0)
           call deallocate(src_to_abs)
        case default
           ! developer error... out of sync options input and code
           FLAbort("Unknown implementation type for k-epsilon source terms")
        end select
     end do
     !-----------------------------------------------------------------------------------

     ! This allows user-specified source and absorption terms, so that an MMS test can be
     ! set up.
     debug => extract_scalar_field(state, &
             trim(field_names(i))//"PrescribedSource", stat)
     if (stat == 0) then
        call addto(src, debug)
     end if
     !-----------------------------------------------------------------------------------

     ! Deallocate fields
     do term = 1, 3
        call deallocate(src_abs_terms(term))
     end do
     call deallocate(fields(1))
     call deallocate(fields(2))

  end do field_loop
  
  call deallocate(dummydensity)
  deallocate(dummydensity)

end subroutine keps_calculate_rhs
    
!------------------------------------------------------------------------------!

subroutine assemble_rhs_ele(src_abs_terms, k, eps, scalar_eddy_visc, u, density, &
     buoyancy_density, have_buoyancy_turbulence, g, g_magnitude, &
     X, c_eps_1, c_eps_2, sigma_p, f_1, f_2, ele, field_id)

  type(scalar_field), dimension(3), intent(inout) :: src_abs_terms
  type(scalar_field), intent(in) :: k, eps, scalar_eddy_visc, f_1, f_2
  type(vector_field), intent(in) :: X, u, g
  type(scalar_field), intent(in) :: density, buoyancy_density
  real, intent(in) :: g_magnitude, c_eps_1, c_eps_2, sigma_p
  logical, intent(in) :: have_buoyancy_turbulence
  integer, intent(in) :: ele, field_id

  real, dimension(ele_loc(k, ele), ele_ngi(k, ele), x%dim) :: dshape
  real, dimension(ele_ngi(k, ele)) :: detwei, rhs, scalar_eddy_visc_ele, k_ele, eps_ele
  real, dimension(3, ele_loc(k, ele)) :: rhs_addto
  integer, dimension(ele_loc(k, ele)) :: nodes
  real, dimension(ele_loc(k, ele), ele_loc(k, ele)) :: invmass
  real, dimension(u%dim, u%dim, ele_ngi(k, ele)) :: reynolds_stress, grad_u
  type(element_type), pointer :: shape
  integer :: term, ngi, dim, gi, i
  
  ! For buoyancy turbulence stuff
  real, dimension(u%dim, ele_ngi(u, ele))  :: vector, u_quad, g_quad
  real :: u_z, u_xy
  real, dimension(ele_ngi(u, ele)) :: scalar, c_eps_3
  type(element_type), pointer :: shape_density
  real, dimension(:, :, :), allocatable :: dshape_density

  shape => ele_shape(k, ele)
  nodes = ele_nodes(k, ele)

  call transform_to_physical( X, ele, shape, dshape=dshape, detwei=detwei )

  ! get bounded values of k and epsilon for source terms
  ! this doesn't change the field values of k and epsilon
  k_ele = ele_val_at_quad(k,ele)
  eps_ele = ele_val_at_quad(eps, ele)
  ngi = ele_ngi(u, ele)
  do gi = 1, ngi
     k_ele(gi) = max(k_ele(gi), fields_min)
     eps_ele(gi) = max(eps_ele(gi), fields_min)
  end do

  ! Compute Reynolds stress
  grad_u = ele_grad_at_quad(u, ele, dshape)
  scalar_eddy_visc_ele = ele_val_at_quad(scalar_eddy_visc, ele)
  dim = u%dim
  do gi = 1, ngi
     reynolds_stress(:,:,gi) = scalar_eddy_visc_ele(gi)*(grad_u(:,:,gi) + transpose(grad_u(:,:,gi)))
  end do
  do i = 1, dim
     reynolds_stress(i,i,:) = reynolds_stress(i,i,:) - (2./3.)*k_ele*ele_val_at_quad(density, ele)
  end do

  ! Compute P
  rhs = tensor_inner_product(reynolds_stress, grad_u)
  if (field_id==2) then
     rhs = rhs*c_eps_1*ele_val_at_quad(f_1,ele)*eps_ele/k_ele
  end if
  rhs_addto(1,:) = shape_rhs(shape, detwei*rhs)

  ! A:
  rhs = -1.0*eps_ele*ele_val_at_quad(density, ele)
  if (field_id==2) then
     rhs = rhs*c_eps_2*ele_val_at_quad(f_2,ele)*eps_ele/k_ele
  end if
  rhs_addto(2,:) = shape_rhs(shape, detwei*rhs)

  ! Gk:  
  ! Calculate buoyancy turbulence term and add to addto array
  if(have_buoyancy_turbulence) then    
  
    ! calculate scalar and vector components of the source term    
    allocate(dshape_density(ele_loc(buoyancy_density, ele), ele_ngi(buoyancy_density, ele), X%dim))
    if(.not.(buoyancy_density%mesh == k%mesh)) then
       shape_density => ele_shape(buoyancy_density, ele)
       call transform_to_physical( X, ele, shape_density, dshape=dshape_density ) 
    else
       dshape_density = dshape
    end if
     
    scalar = -1.0*g_magnitude*ele_val_at_quad(scalar_eddy_visc, ele)/(sigma_p*ele_val_at_quad(density,ele))
    vector = ele_val_at_quad(g, ele)*ele_grad_at_quad(buoyancy_density, ele, dshape_density)
    
    ! multiply vector component by scalar and sum across dimensions - note that the
    ! vector part has been multiplied by the gravitational direction so that it is
    ! zero everywhere apart from in this direction.
    do gi = 1, ngi
       scalar(gi) = sum(scalar(gi) * vector(:, gi))
    end do
   
    if (field_id == 2) then
       ! calculate c_eps_3 = tanh(v/u)
       g_quad = ele_val_at_quad(g, ele)
       u_quad = ele_val_at_quad(u, ele)
       do gi = 1, ngi
          ! get components of velocity in direction of gravity and in other directions
          u_z = dot_product(g_quad(:, gi), u_quad(:, gi))
          u_xy = (norm2(u_quad(:, gi))**2.0 - u_z**2.0)**0.5
          if (u_xy > fields_min) then
             c_eps_3(gi) = tanh(u_z/u_xy) 
          else
             c_eps_3(gi) = 1.0
          end if
       end do     
       scalar = scalar*c_eps_1*ele_val_at_quad(f_1,ele)*c_eps_3*eps_ele/k_ele
    end if

    ! multiply by determinate weights, integrate and assign to rhs
    rhs_addto(3,:) = shape_rhs(shape, scalar * detwei)
    
    deallocate(dshape_density)
    
  else
    ! No buoyancy term, so set this part of the array to zero.
    rhs_addto(3,:) = 0.0    
  end if

  ! In the DG case we apply the inverse mass locally.
  if(continuity(k)<0) then
     invmass = inverse(shape_shape(shape, shape, detwei))
     do term = 1, 3
        rhs_addto(term,:) = matmul(rhs_addto(term,:), invmass)
     end do
  end if

  do term = 1, 3
     call addto(src_abs_terms(term), nodes, rhs_addto(term,:))
  end do

end subroutine assemble_rhs_ele

!------------------------------------------------------------------------------!
! calculate inner product for 2xN matrices dim,dim,N                           !
!------------------------------------------------------------------------------!
function tensor_inner_product(A, B)

  real, dimension(:,:,:), intent(in) :: A, B
  
  real, dimension(size(A,1), size(A,2), size(A,3)) :: C
  real, dimension(size(A,3)) :: tensor_inner_product
  integer :: i
  
  C = A*B
  do i = 1, size(A,3)
     tensor_inner_product(i) = sum(C(:,:,i))
  end do

end function tensor_inner_product

!----------
! eddyvisc calculates the lengthscale and the eddy viscosity
! Eddy viscosity is added to the background viscosity.
!----------
subroutine keps_eddyvisc(state, advdif)

  type(state_type), intent(inout)  :: state
  logical, intent(in) :: advdif

  type(tensor_field), pointer      :: eddy_visc, viscosity, bg_visc
  type(vector_field), pointer      :: x, u
  type(scalar_field)               :: kk, eps
  type(scalar_field), pointer      :: scalar_eddy_visc, ll, f_mu, density, dummydensity
  type(scalar_field)               :: ev_rhs
  integer                          :: i, j, ele, stat
  
  ! Options grabbed from the options tree
  real                             :: c_mu
  character(len=OPTION_PATH_LEN)   :: option_path
  logical                          :: lump_mass, have_visc = .true.
  character(len=FIELD_NAME_LEN)    :: equation_type

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-epsilon/'

  if (.not. have_option(trim(option_path))) then 
     return
  end if

  ewrite(1,*) "In keps_eddyvisc"

  ! Get model constant
  call get_option(trim(option_path)//'/C_mu', c_mu, default = 0.09)
  
  ! Get field data
  call time_averaged_value(state, kk, "TurbulentKineticEnergy", advdif, option_path)
  call time_averaged_value(state, eps, "TurbulentDissipation", advdif, option_path)
  x  => extract_vector_field(state, "Coordinate")
  u          => extract_vector_field(state, "NonlinearVelocity")
  eddy_visc  => extract_tensor_field(state, "EddyViscosity")
  f_mu       => extract_scalar_field(state, "f_mu")
  bg_visc    => extract_tensor_field(state, "BackgroundViscosity")
  scalar_eddy_visc         => extract_scalar_field(state, "ScalarEddyViscosity")
  ll         => extract_scalar_field(state, "LengthScale")
  viscosity  => extract_tensor_field(state, "Viscosity", stat)
  if (stat /= 0) then
     have_visc = .false.
  end if

  ewrite_minmax(kk)
  ewrite_minmax(eps)
  ewrite_minmax(scalar_eddy_visc)
  
  allocate(dummydensity)
  call allocate(dummydensity, X%mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
  call set(dummydensity, 1.0)
  dummydensity%option_path = ""
  
  ! Depending on the equation type, extract the density or set it to some dummy field allocated above
  call get_option(trim(state%option_path)//&
       "/vector_field::Velocity/prognostic/equation[0]/name", equation_type)
  select case(equation_type)
     case("LinearMomentum")
        density=>extract_scalar_field(state, "Density")
     case("Boussinesq")
        density=>dummydensity
     case("Drainage")
        density=>dummydensity
     case default
        ! developer error... out of sync options input and code
        FLAbort("Unknown equation type for velocity")
  end select

  call allocate(ev_rhs, scalar_eddy_visc%mesh, name="EVRHS")
  call zero(ev_rhs)

  ! Initialise viscosity to background value
  if (have_visc) then
     call set(viscosity, bg_visc)
  end if
  
  ! Compute the length scale diagnostic field here.
  do i = 1, node_count(scalar_eddy_visc)
     call set(ll, i, max(node_val(kk,i), fields_min)**1.5 / max(node_val(eps,i), fields_min))
  end do

  ! Calculate scalar eddy viscosity by integration over element
  do ele = 1, ele_count(scalar_eddy_visc)
     call keps_eddyvisc_ele(ele, X, kk, eps, scalar_eddy_visc, f_mu, density, ev_rhs)
  end do

  ! For non-DG we apply inverse mass globally
  if(continuity(scalar_eddy_visc)>=0) then
     lump_mass = have_option(trim(option_path)//'mass_lumping/lump_mass')
     call solve_cg_inv_mass(state, ev_rhs, lump_mass, option_path)  
  end if
  
  ! Allow for prescribed eddy-viscosity
  if (.not. have_option(trim(option_path)//'/scalar_field::ScalarEddyViscosity/prescribed')) then
     call set(scalar_eddy_visc, ev_rhs)
  end if

  call deallocate(ev_rhs)
  call deallocate(kk)
  call deallocate(eps)
  
  call deallocate(dummydensity)
  deallocate(dummydensity)

  ewrite(2,*) "Setting k-epsilon eddy-viscosity tensor"
  call zero(eddy_visc)

  ! this is skipped if zero_eddy_viscosity is set - this is the easiest way to
  ! disable feedback from the k-epsilon model back into the rest of the model
  if (.not. have_option(trim(option_path)//'debugging_options/zero_reynolds_stress_tensor')) then
     do i = 1, eddy_visc%dim(1)
        do j = 1, eddy_visc%dim(1)
           call set(eddy_visc, i, j, scalar_eddy_visc)
        end do
     end do
  end if

  ! Add turbulence model contribution to viscosity field
  if (have_visc) then
     call addto(viscosity, eddy_visc)
  end if

  ewrite_minmax(eddy_visc)
  ewrite_minmax(viscosity)
  ewrite_minmax(scalar_eddy_visc)
  ewrite_minmax(ll)  
  
  contains
  
   subroutine keps_eddyvisc_ele(ele, X, kk, eps, scalar_eddy_visc, f_mu, density, ev_rhs)
   
      type(vector_field), intent(in)   :: x
      type(scalar_field), intent(in)   :: kk, eps, scalar_eddy_visc, f_mu, density
      type(scalar_field), intent(inout):: ev_rhs
      integer, intent(in)              :: ele
      
      type(element_type), pointer      :: shape_ev
      integer, pointer, dimension(:)   :: nodes_ev
      real, dimension(ele_ngi(scalar_eddy_visc, ele)) :: detwei
      real, dimension(ele_loc(scalar_eddy_visc, ele)) :: rhs_addto
      real, dimension(ele_loc(scalar_eddy_visc, ele), ele_loc(scalar_eddy_visc, ele)) :: invmass
      real, dimension(ele_ngi(kk, ele)) :: kk_at_quad, eps_at_quad
      
   
      nodes_ev => ele_nodes(scalar_eddy_visc, ele)
      shape_ev =>  ele_shape(scalar_eddy_visc, ele)
      
      ! Get detwei
      call transform_to_physical(X, ele, detwei=detwei)
      
      ! Get the k and epsilon values at the Gauss points
      kk_at_quad = ele_val_at_quad(kk,ele)
      eps_at_quad = ele_val_at_quad(eps,ele)
      
      ! Clip the field values at the Gauss points.
      ! Note 1: This isn't a permanent change directly to the field itself,
      ! only to the values used in the computation of the eddy viscosity.
      ! Note 2: Can't allow negative/zero epsilon or k.
      ! Note 3: Here we assume all fields have the same number of
      ! Gauss points per element.
      where (kk_at_quad < fields_min)
         kk_at_quad = fields_min
      end where
      where (eps_at_quad < fields_min)
         eps_at_quad = fields_min
      end where
      
      ! Compute the eddy viscosity
      rhs_addto = shape_rhs(shape_ev, detwei*C_mu*ele_val_at_quad(density,ele)*&
                     ele_val_at_quad(f_mu,ele)*(kk_at_quad**2.0)/eps_at_quad)
            
      ! In the DG case we will apply the inverse mass locally.
      if(continuity(scalar_eddy_visc)<0) then
         invmass = inverse(shape_shape(shape_ev, shape_ev, detwei))
         rhs_addto = matmul(rhs_addto, invmass)
      end if
      
      ! Add the element's contribution to the nodes of ev_rhs
      call addto(ev_rhs, nodes_ev, rhs_addto)    
   
   end subroutine keps_eddyvisc_ele

end subroutine keps_eddyvisc

!---------------------------------------------------------------------------------

subroutine keps_diffusion(state)

  ! calculates k and epsilon field diffusivities
  type(state_type), intent(inout)   :: state

  type(tensor_field), pointer :: diff, bg_visc, eddy_visc
  real :: sigma_k, sigma_eps
  integer :: i, j
  character(len=OPTION_PATH_LEN) :: option_path 

  ewrite(1,*) 'in keps_diffusion'
  
  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-epsilon/'

  eddy_visc => extract_tensor_field(state, "EddyViscosity")
  bg_visc => extract_tensor_field(state, "BackgroundViscosity")
  call get_option(trim(option_path)//'/sigma_k', sigma_k, default = 1.0)
  call get_option(trim(option_path)//'/sigma_eps', sigma_eps, default = 1.3)

  ! Set diffusivity
  diff => extract_tensor_field(state, "TurbulentKineticEnergyDiffusivity")
  call zero(diff)
  do i = 1, node_count(diff)
     do j = 1, diff%dim(1)
        call addto(diff, j, j, i, node_val(bg_visc, j, j, i))
        call addto(diff, j, j, i, node_val(eddy_visc, j, j, i) / sigma_k)
     end do
  end do
  diff => extract_tensor_field(state, "TurbulentDissipationDiffusivity")
  call zero(diff)
  do i = 1, node_count(diff)
     do j = 1, diff%dim(1)
        call addto(diff, j, j, i, node_val(bg_visc, j, j, i))
        call addto(diff, j, j, i, node_val(eddy_visc, j, j, i) / sigma_eps)
     end do
  end do

end subroutine keps_diffusion

!---------------------------------------------------------------------------------

subroutine keps_tracer_diffusion(state)

  ! calculates scalar field diffusivity based upon eddy viscosity and background
  !  diffusivity
  type(state_type), intent(inout)   :: state

  type(tensor_field), pointer       :: t_field
  integer                           :: i_field, i, stat
  real                              :: sigma_p, local_background_diffusivity
  type(scalar_field)                :: local_background_diffusivity_field
  type(scalar_field), pointer       :: scalar_eddy_viscosity, s_field
  type(tensor_field), pointer       :: global_background_diffusivity
  type(tensor_field)                :: background_diffusivity

  ewrite(1,*) 'In keps_tracer_diffusion'

  do i_field = 1, scalar_field_count(state)
     s_field => extract_scalar_field(state, i_field)

     if (have_option(trim(s_field%option_path)//&
          '/prognostic/subgridscale_parameterisation::k-epsilon')) then

        ewrite(1,*) 'Calculating turbulent diffusivity for field: ', s_field%name
        
        ! check options
        if (.not.(have_option(trim(state%option_path)//'/subgridscale_parameterisations/k-epsilon')))&
             & then
           FLExit('you must have /subgridscale_parameterisations/k-epsilon to be able to calculate diffusivity based upon the k-epsilon model')
        end if

        t_field => extract_tensor_field(state, trim(s_field%name)//'Diffusivity', stat=stat) 
        if (stat /= 0) then
           FLExit('you must have a Diffusivity field to be able to calculate diffusivity based upon the k-epsilon model')        
        else if (.not. have_option(trim(t_field%option_path)//"/diagnostic/algorithm::Internal")) then
           FLExit('you must have a diagnostic Diffusivity field with algorithm::Internal to be able to calculate diffusivity based upon the k-epsilon model')  
        end if

        ! get sigma_p number
        call get_option(trim(state%option_path)//'/subgridscale_parameterisations/k-epsilon/sigma_p', sigma_p)

        ! allocate and zero required fields
        call allocate(background_diffusivity, t_field%mesh, name="background_diffusivity")
        call zero(background_diffusivity)
        call allocate(local_background_diffusivity_field, t_field%mesh, &
             name="local_background_diffusivity_field")
        call zero(local_background_diffusivity_field)

        ! set background_diffusivity (local takes precendence over global)
        call get_option(trim(s_field%option_path)//&
             '/prognostic/subgridscale_parameterisation::k-epsilon/background_diffusivity', &
             local_background_diffusivity, stat=stat)
        if (stat == 0) then 
           ! set local isotropic background diffusivity
           call addto(local_background_diffusivity_field, local_background_diffusivity)
           do i = 1, background_diffusivity%dim(1)
              call set(background_diffusivity, i, i, local_background_diffusivity_field)
           end do
        else
           global_background_diffusivity => extract_tensor_field(state, 'BackgroundDiffusivity', stat=stat)
           if (stat == 0) then 
              call set(background_diffusivity, global_background_diffusivity)
           end if
        end if

        ! get eddy viscosity
        scalar_eddy_viscosity => extract_scalar_field(state, 'ScalarEddyViscosity', stat)

        call zero(t_field)
        call addto(t_field, background_diffusivity)
        do i = 1, t_field%dim(1)
           call addto(t_field, i, i, scalar_eddy_viscosity, 1.0/sigma_p)
        end do

        call deallocate(background_diffusivity)
        call deallocate(local_background_diffusivity_field)

     end if
  end do

end subroutine keps_tracer_diffusion

!--------------------------------------------------------------------------------!
! calculates the reynolds stress tensor correction term 
! - grad(2/3 k rho delta(ij)) 
! this is added to prescribed momentum source fields or sets diagnostic source 
! fields 
!--------------------------------------------------------------------------------!
subroutine keps_momentum_source(state)

  type(state_type), intent(inout)  :: state
  
  type(scalar_field), pointer :: dummydensity, density
  type(scalar_field) :: k
  type(vector_field), pointer :: source, x, u
  type(vector_field) :: rhs
  integer :: i
  logical :: lump_mass
  character(len=OPTION_PATH_LEN) :: option_path 
  character(len=FIELD_NAME_LEN)    :: equation_type

  ewrite(1,*) 'In keps_momentum_source'

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-epsilon/'

  ! exit if there is no k-epsilon model or if the velocity field is prescribed (hence no source term)
  if (.not.have_option(trim(option_path)) .or. &
      have_option(trim(state%option_path)//"/vector_field::Velocity/prescribed")) then 
     return
  end if

  source => extract_vector_field(state, "VelocitySource")
  call zero(source)

  ! exit after zero'ing field if zero_reynolds_stress_tensor is disabled
  ! zero'd because this is what the calling subroutine expects - the calling routine handles reapplying
  ! prescribed source terms
  if (have_option(trim(option_path)//'debugging_options/zero_reynolds_stress_tensor')) then 
     return
  end if
  
  x => extract_vector_field(state, "Coordinate")
  u => extract_vector_field(state, "NonlinearVelocity")
  call time_averaged_value(state, k, "TurbulentKineticEnergy", .false., option_path)
  
  allocate(dummydensity)
  call allocate(dummydensity, x%mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
  call set(dummydensity, 1.0)
  dummydensity%option_path = ""
  
  ! Depending on the equation type, extract the density or set it to some dummy field allocated above
  call get_option(trim(state%option_path)//&
       "/vector_field::Velocity/prognostic/equation[0]/name", equation_type)
  select case(equation_type)
     case("LinearMomentum")
        density=>extract_scalar_field(state, "Density")
     case("Boussinesq")
        density=>dummydensity
     case("Drainage")
        density=>dummydensity
     case default
        ! developer error... out of sync options input and code
        FLAbort("Unknown equation type for velocity")
  end select

  call allocate(rhs, source%dim, source%mesh, name='TempSource')
  call zero(rhs)
  do i = 1, ele_count(k)
     call keps_momentum_source_ele()
  end do

  ! For non-DG we apply inverse mass globally
  if(continuity(k)>=0) then
     lump_mass = have_option(trim(option_path)//&
          'mass_lumping/lump_mass')
     call solve_cg_inv_mass_vector(state, rhs, lump_mass, option_path)  
     call addto(source, rhs)
  end if

  ewrite_minmax(source)
  
  call deallocate(rhs)
  call deallocate(k)
  call deallocate(dummydensity)
  deallocate(dummydensity)

 contains
    
  subroutine keps_momentum_source_ele()
      
    real, dimension(ele_loc(k, i), ele_ngi(k, i), x%dim) :: dshape_k
    real, dimension(ele_loc(density, i), ele_ngi(density, i), x%dim) :: dshape_density
    real, dimension(ele_ngi(k, i)) :: detwei
    integer, dimension(ele_loc(k, i)) :: nodes
    real, dimension(ele_loc(k, i), ele_loc(k, i)) :: invmass
    type(element_type), pointer :: shape_k, shape_density
    real, dimension(x%dim, ele_loc(k, i)) :: rhs_addto
    real, dimension(x%dim, ele_ngi(k, i)) :: grad_k
    real, dimension(x%dim, ele_ngi(density, i)) :: grad_density

    shape_k => ele_shape(k, i)
    nodes = ele_nodes(source, i)

    call transform_to_physical( x, i, shape_k, dshape=dshape_k, detwei=detwei )

    if(.not.(density%mesh == k%mesh)) then
       shape_density => ele_shape(density, i)
       call transform_to_physical( x, i, shape_density, dshape=dshape_density ) 
    else
       dshape_density = dshape_k
    end if
     
    grad_k = ele_grad_at_quad(k, i, dshape_k)
    grad_density = ele_grad_at_quad(density, i, dshape_density)
    ! IMPORTANT: This gets added to the VelocitySource term. In Momentum_CG/DG.F90
    ! the VelocitySource gets multiplied by the Density field, so we have
    ! divided through by Density here in order to cancel this out
    ! and get the desired result.
    rhs_addto = shape_vector_rhs(shape_k, -(2./3.)*grad_k, detwei) + &
                shape_vector_rhs(shape_k, -(2./3.)*grad_density, detwei*ele_val_at_quad(k,i)/ele_val_at_quad(density,i)) 
      
    ! In the DG case we apply the inverse mass locally.
    if(continuity(k)<0) then
       invmass = inverse(shape_shape(shape_k, shape_k, detwei))
       rhs_addto = matmul(rhs_addto, invmass)
       call addto(source, nodes, rhs_addto)  
    else
       call addto(rhs, nodes, rhs_addto) 
    end if

  end subroutine keps_momentum_source_ele

end subroutine keps_momentum_source

!--------------------------------------------------------------------------------!
! This gets and applies locally defined boundary conditions (wall functions)     !
!--------------------------------------------------------------------------------!

subroutine keps_bcs(state)

  type(state_type), intent(in)               :: state
  type(scalar_field), pointer                :: field1, field2    ! k or epsilon
  type(scalar_field), pointer                :: f_1, f_2, f_mu
  type(scalar_field), pointer                :: surface_field, scalar_eddy_visc
  type(scalar_field), pointer                :: density, dummydensity
  type(vector_field), pointer                :: X, u
  type(tensor_field), pointer                :: bg_visc
  type(scalar_field)                         :: rhs_field, surface_values, masslump
  type(mesh_type), pointer                   :: surface_mesh
  integer                                    :: i, j, sele, index, nbcs, stat
  integer, dimension(:), pointer             :: surface_elements, surface_node_list
  character(len=FIELD_NAME_LEN)              :: bc_type, bc_name, wall_fns
  character(len=OPTION_PATH_LEN)             :: bc_path, bc_path_i, option_path 
  real                                       :: c_mu
  character(len=FIELD_NAME_LEN)              :: equation_type

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-epsilon/'

  ewrite(2,*) "In keps_bcs"

  X => extract_vector_field(state, "Coordinate")
  u                 => extract_vector_field(state, "Velocity")
  scalar_eddy_visc  => extract_scalar_field(state, "ScalarEddyViscosity")
  bg_visc           => extract_tensor_field(state, "BackgroundViscosity")
  f_1               => extract_scalar_field(state, "f_1")
  f_2               => extract_scalar_field(state, "f_2")
  f_mu              => extract_scalar_field(state, "f_mu")

  allocate(dummydensity)
  call allocate(dummydensity, X%mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
  call set(dummydensity, 1.0)
  dummydensity%option_path = ""
  
  ! Depending on the equation type, extract the density or set it to some dummy field allocated above
  call get_option(trim(u%option_path)//"/prognostic/equation[0]/name", equation_type)
  select case(equation_type)
     case("LinearMomentum")
        density=>extract_scalar_field(state, "Density")
     case("Boussinesq")
        density=>dummydensity
     case("Drainage")
        density=>dummydensity
     case default
        ! developer error... out of sync options input and code
        FLAbort("Unknown equation type for velocity")
  end select

  call get_option(trim(option_path)//"C_mu", c_mu, default = 0.09)

  field_loop: do index=1,2

     if(index==1) then
        field1 => extract_scalar_field(state, "TurbulentKineticEnergy")
        field2 => null()
     else
        field1 => extract_scalar_field(state, "TurbulentDissipation")
        field2 => extract_scalar_field(state, "TurbulentKineticEnergy")
     end if

     bc_path=trim(field1%option_path)//'/prognostic/boundary_conditions'
     nbcs=option_count(trim(bc_path))

     ! Loop over boundary conditions for field1
     boundary_conditions: do i=0, nbcs-1

        bc_path_i=trim(bc_path)//"["//int2str(i)//"]"

        ! Get name and type of boundary condition
        call get_option(trim(bc_path_i)//"/name", bc_name)
        call get_option(trim(bc_path_i)//"/type[0]/name", bc_type)
        ! Do we have high- or low-Reynolds-number wall functions?
        call get_option(trim(bc_path_i)//"/type::k_epsilon/", wall_fns, stat=stat)

        if (trim(bc_type)=="k_epsilon" .and. wall_fns=="low_Re") then
           
           ! lowRe BC's are just zero Dirichlet or Neumann - damping functions get calculated in 
           ! keps_calc_rhs
           low_Re = .true.

        else if (trim(bc_type)=="k_epsilon" .and. wall_fns=="high_Re") then
           
           ! Get bc by name. Get type just to make sure it's now dirichlet
           call get_boundary_condition(field1, name=bc_name, type=bc_type, surface_node_list=surface_node_list, &
                surface_element_list=surface_elements, surface_mesh=surface_mesh)
           
           ewrite(2,*) "Calculating highRe k-epsilon BC: ",trim(field1%name),' ',trim(bc_name),' ',trim(bc_type),' ',trim(wall_fns)

           ! Get surface field already created in bcs_from_options
           surface_field => extract_surface_field(field1, bc_name=bc_name, name="value")
           call zero(surface_field)

           call allocate(surface_values, surface_mesh, name="surfacevalues")
           call allocate(rhs_field, field1%mesh, name="rhs")
           call zero(surface_values); call zero(rhs_field)

           if(continuity(field1)>=0) then
              call allocate(masslump, field1%mesh, 'Masslump')
              call zero(masslump)
           end if

           ! Calculate high-Re wall function
           do sele = 1, ele_count(surface_mesh)
              call keps_wall_function(field1, X, u, masslump, scalar_eddy_visc, density, sele, &
                   index, c_mu, rhs_field)
           end do
           
           ! In the continuous case we globally apply inverse mass
           if(continuity(field1)>=0) then
             where (masslump%val/=0.0)
               masslump%val=1./masslump%val
             end where
             call scale(rhs_field, masslump)
             call deallocate(masslump)
           end if

           ! Put values onto surface mesh
           call remap_field_to_surface(rhs_field, surface_values, surface_elements)
           ewrite_minmax(rhs_field)
           do j = 1, size(surface_node_list)
              call set(surface_field, j, node_val(surface_values, j))
           end do

           call deallocate(surface_values); call deallocate(rhs_field)

        end if
     end do boundary_conditions
  end do field_loop

  call deallocate(dummydensity)
  deallocate(dummydensity)

end subroutine keps_bcs

!--------------------------------------------------------------------------------!
! Only used if bc type == k_epsilon for field and high_Re.                       !
!--------------------------------------------------------------------------------!
subroutine keps_wall_function(field1,X,U,masslump,scalar_eddy_visc,density,sele,index,c_mu,rhs_field)

  type(scalar_field), pointer, intent(in)              :: field1, scalar_eddy_visc, density
  type(vector_field), pointer, intent(in)              :: X, U
  type(scalar_field), intent(inout)                    :: masslump, rhs_field
  integer, intent(in)                                  :: sele, index
  type(element_type), pointer                          :: shape, f_shape
  type(element_type)                                   :: augmented_shape
  integer                                              :: i, j, gi, ele, dim
  integer, dimension(face_loc(rhs_field, sele))        :: nodes_bdy
  real                                                 :: kappa, h, c_mu
  real, dimension(1,1)                                 :: hb
  real, dimension(face_ngi(field1,sele))               :: detwei, tau_wall_at_quad, visc_at_quad, density_at_quad
  real, dimension(X%dim,1)                             :: n
  real, dimension(X%dim,X%dim)                         :: G
  real, dimension(face_loc(rhs_field, sele))           :: rhs
  real, dimension(X%dim,face_ngi(field1,sele))         :: normal, normal_shear_at_quad
  real, dimension(X%dim, X%dim, face_ngi(field1, sele)):: f_invJ, grad_U_at_quad, shear_at_quad
  real, dimension(face_loc(field1, sele), face_loc(field1, sele))                     :: fmass
  real, dimension(X%dim, X%dim, ele_ngi(field1, face_ele(field1, sele)))              :: invJ
  real, dimension(ele_loc(X, face_ele(field1, sele)), face_ngi(field1, sele), X%dim)  :: ele_dshape_at_face_quad

  nodes_bdy = face_global_nodes(rhs_field, sele) ! nodes in rhs_field
  ele       = face_ele(field1, sele) ! ele number for volume mesh
  dim       = mesh_dim(X) ! field dimension 

  ! get shape functions
  f_shape => face_shape(field1, sele)  
  shape   => ele_shape(field1, ele)
  
  ! generate shape functions that include quadrature points on the face required
  ! check that the shape does not already have these first
  assert(shape%degree == 1)
  if(associated(shape%dn_s)) then
    augmented_shape = shape
    call incref(augmented_shape)
  else
     augmented_shape = make_element_shape(shape%loc, shape%dim, shape%degree, shape%quadrature, quad_s=f_shape%quadrature)
  end if
    
  ! assumes that the jacobian is the same for all quadrature points
  ! this is not valid for spheres!
  call compute_inverse_jacobian(ele_val(X, ele), ele_shape(X, ele), invj = invJ)
  assert(ele_numbering_family(shape) == FAMILY_SIMPLEX)
  f_invJ = spread(invJ(:, :, 1), 3, size(f_invJ, 3))
   
  call transform_facet_to_physical(X, sele, detwei_f = detwei, normal = normal)

  ! Evaluate the volume element shape function derivatives at the surface
  ! element quadrature points
  ele_dshape_at_face_quad = eval_volume_dshape_at_face_quad(augmented_shape, &
    & local_face_number(X, sele), f_invJ)

  ! Calculate grad of U at the surface element quadrature points
  do i=1, dim
    do j=1, dim
      grad_U_at_quad(i, j, :) = matmul(ele_val(U, j, ele), ele_dshape_at_face_quad(:,:,i))
    end do
  end do

  density_at_quad = face_val_at_quad(density, sele)
  visc_at_quad = face_val_at_quad(scalar_eddy_visc, sele)
  do gi = 1, face_ngi(X, sele)
    ! Multiply by visosity
    shear_at_quad(:,:,gi) = grad_U_at_quad(:,:,gi)*visc_at_quad(gi)
    ! Multiply by surface normal (dim,sgi) to obtain shear in direction normal
    ! to surface - transpose (because fluidity stores data in row-major order??)
    normal_shear_at_quad(:,gi) = matmul(shear_at_quad(:,:,gi), normal(:,gi))
    ! Wall shear stress tau_wall/rho = nu_T*du/dn.
    ! Get streamwise velocity gradient by taking sqrt(grad_n.grad_n).
    ! Eddy viscosity is dynamic not kinematic so divide by density.
    tau_wall_at_quad(gi) = norm2(normal_shear_at_quad(:,gi))
  end do

  if (index==1) then
    ! k = tau_wall/rho/c_mu**0.5
    rhs = shape_rhs(f_shape, tau_wall_at_quad/density_at_quad/c_mu**0.5*detwei)
  else if (index==2) then
    ! calculate wall-normal element size
    G = matmul(transpose(invJ(:,:,1)), invJ(:,:,1))
    n(:,1) = normal(:,1)
    hb = 1. / sqrt( matmul(matmul(transpose(n), G), n) )
    h  = hb(1,1)
    ! Von Karman's constant
    kappa = 0.43
    ! epsilon = (tau_wall/rho)**1.5/kappa/h
    rhs = shape_rhs(f_shape, (tau_wall_at_quad/density_at_quad)**1.5/kappa/h*detwei)  
  end if

  fmass = shape_shape(f_shape, f_shape, detwei)
  ! In the CG case we need to calculate a global lumped mass
  if(continuity(field1)>=0) then
    call addto(masslump, nodes_bdy, sum(fmass,1))
  else ! In the DG case we will apply the inverse mass locally.
    do i = 1,dim
    rhs = matmul(inverse(fmass), rhs)
    end do
  end if

  ! add to rhs field
  call addto(rhs_field, nodes_bdy, rhs)

  call deallocate(augmented_shape)

end subroutine keps_wall_function

!---------------------------------------------------------------------------------

subroutine time_averaged_value(state, A, field_name, advdif, option_path)
  
  type(state_type), intent(in) :: state
  type(scalar_field), intent(inout) :: A
  character(len=*), intent(in) :: field_name
  logical, intent(in) :: advdif    ! advdif or mom - whether to use old or iterated values
  character(len=OPTION_PATH_LEN), intent(in) :: option_path 

  real :: theta
  type(scalar_field), pointer :: old, iterated
  
  call get_option(trim(option_path)//'time_discretisation/theta', theta)

  old => extract_scalar_field(state, "Old"//trim(field_name))
  if (advdif) then
     iterated => extract_scalar_field(state, "Iterated"//trim(field_name))
  else
     iterated => extract_scalar_field(state, trim(field_name))
  end if

  call allocate(A, old%mesh, name="Nonlinear"//trim(field_name))
  call zero(A)
  call addto(A, old, 1.0-theta)
  call addto(A, iterated, theta)

  ewrite_minmax(old)
  ewrite_minmax(iterated)
  ewrite_minmax(A)

end subroutine time_averaged_value

!---------------------------------------------------------------------------------

subroutine solve_cg_inv_mass(state, A, lump, option_path)
  
  type(state_type), intent(inout) :: state
  type(scalar_field), intent(inout) :: A
  logical, intent(in) :: lump
  character(len=OPTION_PATH_LEN), intent(in) :: option_path 

  type(scalar_field), pointer :: lumped_mass
  type(csr_matrix), pointer :: mass_matrix
  type(scalar_field) :: inv_lumped_mass, x
  
  if (lump) then
     call allocate(inv_lumped_mass, A%mesh)
     lumped_mass => get_lumped_mass(state, A%mesh)
     call invert(lumped_mass, inv_lumped_mass)
     call scale(A, inv_lumped_mass)
     call deallocate(inv_lumped_mass)
  else
     call allocate(x, A%mesh)
     mass_matrix => get_mass_matrix(state, A%mesh)
     call petsc_solve(x, mass_matrix, A, &
          trim(option_path)//&
          'mass_lumping/solve_using_mass_matrix/')
     call set(A, x)
     call deallocate(x)
  end if

end subroutine solve_cg_inv_mass

!---------------------------------------------------------------------------------

subroutine solve_cg_inv_mass_vector(state, A, lump, option_path)
  
  type(state_type), intent(inout) :: state
  type(vector_field), intent(inout) :: A
  logical, intent(in) :: lump
  character(len=OPTION_PATH_LEN), intent(in) :: option_path 

  type(scalar_field), pointer :: lumped_mass
  type(csr_matrix), pointer :: mass_matrix
  type(scalar_field) :: inv_lumped_mass
  type(vector_field) :: x
  
  if (lump) then
     call allocate(inv_lumped_mass, A%mesh)
     lumped_mass => get_lumped_mass(state, A%mesh)
     call invert(lumped_mass, inv_lumped_mass)
     call scale(A, inv_lumped_mass)
     call deallocate(inv_lumped_mass)
  else
     call allocate(x, A%dim, A%mesh)
     mass_matrix => get_mass_matrix(state, A%mesh)
     call petsc_solve(x, mass_matrix, A, &
          trim(option_path)//&
          'mass_lumping/solve_using_mass_matrix/')
     call set(A, x)
     call deallocate(x)
  end if

end subroutine solve_cg_inv_mass_vector

!---------------------------------------------------------------------------------

subroutine k_epsilon_check_options

  character(len=OPTION_PATH_LEN) :: option_path
  character(len=FIELD_NAME_LEN)  :: kmsh, emsh, vmsh
  integer                        :: dimension, stat, n_phases, istate

  ewrite(1,*) "In keps_check_options"

  n_phases = option_count("/material_phase")
  
  if(option_count("/material_phase/subgridscale_parameterisations/k-epsilon") > 1) then
     FLExit("The k-epsilon model can only be applied to a single-phase.")
  end if

  do istate = 0, n_phases-1

     option_path = "/material_phase["//int2str(istate)//"]/subgridscale_parameterisations/k-epsilon"
     
     ! one dimensional problems not supported
     call get_option("/geometry/dimension/", dimension) 
     if (dimension .eq. 1 .and. have_option(trim(option_path))) then
        FLExit("k-epsilon model is only supported for dimension > 1")
     end if
     ! Don't do k-epsilon if it's not included in the model!
     if (.not.have_option(trim(option_path))) return

     ! checking for required fields
     if (have_option(trim(option_path)//"/scalar_field::TurbulentKineticEnergy/prognostic")) then
        ! diffusivity is on and diagnostic
        if (.not.have_option(trim(option_path)//"/scalar_field::TurbulentKineticEnergy"//&
             &"/prognostic/tensor_field::Diffusivity")) then
           FLExit("You need TurbulentKineticEnergy Diffusivity field for k-epsilon")
        end if
        if (.not.have_option(trim(option_path)//&
             &"/scalar_field::TurbulentKineticEnergy/prognostic/"//&
             &"/tensor_field::Diffusivity/diagnostic/algorithm::Internal")) then
           FLExit("You need TurbulentKineticEnergy Diffusivity field set to diagnostic/internal")
        end if
        ! source terms
        if (.not.have_option(trim(option_path)//&
             &"/scalar_field::TurbulentKineticEnergy/prognostic"//&
             &"/scalar_field::Source")) then
           FLExit("You need TurbulentKineticEnergy Source field for k-epsilon")
        end if
        if (.not. have_option(trim(option_path)//&
             &"/scalar_field::TurbulentKineticEnergy/prognostic"//&
             &"/scalar_field::Source/diagnostic/algorithm::Internal")) then
           FLExit("You need TurbulentKineticEnergy Source field set to diagnostic/internal")
        end if
        ! absorption terms
        if (.not.have_option(trim(option_path)//&
             &"/scalar_field::TurbulentKineticEnergy/prognostic"//&
             &"/scalar_field::Absorption")) then
           FLExit("You need TurbulentKineticEnergy Absorption field for k-epsilon")
        end if
        if (.not.have_option(trim(option_path)//&
             &"/scalar_field::TurbulentKineticEnergy/prognostic"//&
             &"/scalar_field::Absorption/diagnostic/algorithm::Internal")) then
           FLExit("You need TurbulentKineticEnergy Absorption field set to diagnostic/internal")
        end if
     else if (have_option(trim(option_path)// &
          "/scalar_field::TurbulentKineticEnergy/prescribed")) then
        ewrite(0,*) "WARNING: TurbulentKineticEnergy field is prescribed"
     else
        FLExit("You need prognostic/prescribed TurbulentKineticEnergy field for k-epsilon")
     end if
     if (have_option(trim(option_path)//"/scalar_field::TurbulentDissipation/prognostic")) then
        ! diffusivity is on and diagnostic
        if (.not.have_option(trim(option_path)//"/scalar_field::TurbulentDissipation"//&
             &"/prognostic/tensor_field::Diffusivity")) then
           FLExit("You need TurbulentDissipation Diffusivity field for k-epsilon")
        end if
        if (.not.have_option(trim(option_path)//&
             &"/scalar_field::TurbulentDissipation/prognostic/"//&
             &"/tensor_field::Diffusivity/diagnostic/algorithm::Internal")) then
           FLExit("You need TurbulentDissipation Diffusivity field set to diagnostic/internal")
        end if
        ! source terms
        if (.not.have_option(trim(option_path)//&
             &"/scalar_field::TurbulentDissipation/prognostic"//&
             &"/scalar_field::Source")) then
           FLExit("You need TurbulentDissipation Source field for k-epsilon")
        end if
        if (.not. have_option(trim(option_path)//&
             &"/scalar_field::TurbulentDissipation/prognostic"//&
             &"/scalar_field::Source/diagnostic/algorithm::Internal")) then
           FLExit("You need TurbulentDissipation Source field set to diagnostic/internal")
        end if
        ! absorption terms
        if (.not.have_option(trim(option_path)//&
             &"/scalar_field::TurbulentDissipation/prognostic"//&
             &"/scalar_field::Absorption")) then
           FLExit("You need TurbulentDissipation Absorption field for k-epsilon")
        end if
        if (.not.have_option(trim(option_path)//&
             &"/scalar_field::TurbulentDissipation/prognostic"//&
             &"/scalar_field::Absorption/diagnostic/algorithm::Internal")) then
           FLExit("You need TurbulentDissipation Absorption field set to diagnostic/internal")
        end if
     else if (have_option(trim(option_path)// &
          "/scalar_field::TurbulentDissipation/prescribed")) then
        ewrite(0,*) "WARNING: TurbulentDissipation field is prescribed"
     else
        FLExit("You need prognostic/prescribed TurbulentDissipation field for k-epsilon")
     end if

     ! Check that TurbulentKineticEnergy and TurbulentDissipation fields are on the same
     !  mesh as the velocity
     call get_option(trim(option_path)//&
          &"/scalar_field::TurbulentKineticEnergy/prognostic/mesh/name", kmsh, stat)
     if (stat /= 0) then
        call get_option(trim(option_path)//&
             &"/scalar_field::TurbulentKineticEnergy/prescribed/mesh/name", kmsh,&
             & stat)
     end if
     call get_option(trim(option_path)//&
          &"/scalar_field::TurbulentDissipation/prognostic/mesh/name", emsh, stat)
     if (stat /= 0) then
        call get_option(trim(option_path)//&
             &"/scalar_field::TurbulentDissipation/prescribed/mesh/name", emsh,&
             & stat)
     end if
     call get_option("/material_phase["//int2str(istate)//"]/vector_field::Velocity/prognostic/mesh/name", vmsh,&
          & stat)
     if (stat /= 0) then
        call get_option("/material_phase["//int2str(istate)//"]/vector_field::Velocity/prescribed/mesh/name", vmsh,&
             & stat)
        if (stat /= 0) then
           FLExit("You must use a prognostic or prescribed Velocity field")
        end if
     end if
     if(.not. kmsh==emsh .or. .not. kmsh==vmsh .or. .not. emsh==vmsh) then
        FLExit("You must use the Velocity mesh for TurbulentKineticEnergy and TurbulentDissipation fields")
     end if

     ! Velocity field options
     if (.not.have_option("/material_phase["//int2str(istate)//"]/vector_field::Velocity/prognostic"//&
          "/tensor_field::Viscosity/") .and. &
          .not.have_option("/material_phase["//int2str(istate)//"]/vector_field::Velocity/prescribed")) then
        FLExit("Need viscosity switched on under the Velocity field for k-epsilon.") 
     end if
     ! check that the user has switched Velocity/viscosity to diagnostic
     if (have_option("/material_phase["//int2str(istate)//"]/vector_field::Velocity/prognostic") .and. &
          .not.have_option("/material_phase["//int2str(istate)//"]/vector_field::Velocity/prognostic"//&
          "/tensor_field::Viscosity/diagnostic/")) then
        FLExit("You need to switch the viscosity field under Velocity to diagnostic/internal")
     end if
     ! check that the user has enabled a Velocity Source field
     if (have_option("/material_phase["//int2str(istate)//"]/vector_field::Velocity/prognostic") .and. &
          .not.have_option("/material_phase["//int2str(istate)//"]/vector_field::Velocity/prognostic"//&
          &"/vector_field::Source/")) then
        FLExit("A velocity source field is required for the reynolds stress adjustment (-2/3 k delta(ij))")
     end if

     ! Check ScalarEddyViscosity is diagnostic
     if (have_option(trim(option_path)//'/scalar_field::ScalarEddyViscosity/prescribed')) then
        ewrite(0,*) "WARNING: ScalarEddyViscosity field is prescribed"
     end if

  end do

end subroutine k_epsilon_check_options

end module k_epsilon
