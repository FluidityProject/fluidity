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

module k_omega
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN, timestep, current_time
  use fldebug
  use futils, only: int2str
  use vector_tools
  use quadrature
  use spud
  use sparse_tools
  use elements
  use fetools
  use parallel_fields
  use fields
  use state_module
  use boundary_conditions
  use vtk_interfaces
  use field_derivatives
  use field_options
  use sparsity_patterns_meshes
  use state_fields_module
  use surface_integrals
  use solvers
  use smoothing_module 
  use fields_manipulation
  use fetools

implicit none

  private

  ! locally allocatad fields
  real, save     :: fields_min = 1.0e-11 !A! remove
  !A!logical, save  :: low_Re = .false.

  public :: komega_advdif_diagnostics, komega_momentum_diagnostics, & !A! komega_bcs removed
       & k_omega_check_options, tensor_inner_product

  ! Outline:
  !  - call diagnostics to obtain source terms and calculate eddy viscosity
  !  - after each scalar field solve recalculates the eddy viscosity

contains

subroutine komega_advdif_diagnostics(state)            !A! called in Fluids.F90

  type(state_type), intent(inout) :: state
  
  call komega_blending_functions(state, advdif=.true.) !A! blending functions for k-omega SST !A!
  call komega_eddyvisc(state, advdif=.true.)           !A! eddy viscosity
  call komega_diffusion(state)                         !A! diff coeff
  call komega_tracer_diffusion(state)                  !A! ???
  call komega_calculate_rhs(state)                     !A! source/sink terms on RHS of adv-diff eqs

end subroutine komega_advdif_diagnostics

subroutine komega_momentum_diagnostics(state)          !A! called in Momentum_Diagnostic_Fields.F90

  type(state_type), intent(inout) :: state
  
  call komega_eddyvisc(state, advdif=.false.)          !A! eddy viscosity

end subroutine komega_momentum_diagnostics
 
!--------------------------------------------------------------------------------!

subroutine komega_blending_functions(state, advdif) !A! used for k-omega SST only

  type(state_type), intent(in)   :: state
  logical, intent(in)            :: advdif
  type(vector_field), pointer    :: X
  type(scalar_field), pointer    :: F_1, F_2, y, dummydensity, density
  type(scalar_field)             :: k, omega
!  type(scalar_field), target     :: dummydensity
  type(vector_field), pointer    :: grad_k, grad_omg
  type(scalar_field), pointer    :: CD_komg
  type(tensor_field), pointer    :: bg_visc
  integer                        :: stat, ele, node
  real                           :: beta_star, sigma_omg2
  character(len=FIELD_NAME_LEN)  :: equation_type
  character(len=OPTION_PATH_LEN) :: option_path

  real                           :: a_3, a_4, a_5, arg_1, arg_2
  real                           :: CD_komg_val, F_1_val, F_2_val
  logical                        :: have_SST

  ewrite(1,*) 'in komega_blending_functions'

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-omega/'

  have_SST = have_option(trim(state%option_path)//"/subgridscale_parameterisations/k-omega/k-omega_SST")
  ewrite(1,*) 'k-omega SST? ', have_SST

  if(have_SST) then
    F_1 => extract_scalar_field(state, "F_1")
    F_2 => extract_scalar_field(state, "F_2")
    grad_k   => extract_vector_field(state, "Grad_k")
    grad_omg => extract_vector_field(state, "Grad_omg")
    CD_komg  => extract_scalar_field(state, "CD_komg")

    !A! initialise blending functions
    call set(F_1, 1.0)
    call set(F_2, 1.0)

    X       => extract_vector_field(state, "Coordinate")
    bg_visc => extract_tensor_field(state, "BackgroundViscosity")
    y       => extract_scalar_field(state, "DistanceToWall", stat = stat)

    call time_averaged_value(state, k, 'TurbulentKineticEnergy', advdif, option_path)
    call time_averaged_value(state, omega, 'TurbulentFrequency', advdif, option_path)

    call get_option(trim(option_path)//'/Beta_Star', beta_star, default = 0.09)
    call get_option(trim(option_path)//'/Sigma_Omg2', sigma_omg2, default = 0.856)

    !A! set dummy density to 1.0
    allocate(dummydensity)
    call allocate(dummydensity, F_1%mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
    call set(dummydensity, 1.0)
    dummydensity%option_path = ""
  
    ! Depending on the equation type, extract the density or set it to some dummy field allocated above
    call get_option(trim(state%option_path)//&
         "/vector_field::Velocity/prognostic/equation[0]/name", equation_type, stat=stat)
    if (stat /= 0) then
      density=>dummydensity
    else
      select case(equation_type)
        case("LinearMomentum")
          density=>extract_scalar_field(state, "Density")
        case("Boussinesq") !A! this is our current case!
          density=>dummydensity
        case("Drainage")
          density=>dummydensity
        case default
          ! developer error... out of sync options input and code
          FLAbort("Unknown equation type for velocity")
      end select
    end if

    call grad(k, X, grad_k)                       !A! grad(k)
    call grad(omega, X, grad_omg)                 !A! grad(omega)
    call inner_product(CD_komg, grad_k, grad_omg) !A! grad(k).grad(omega)

    node_loop: do node = 1, node_count(k)

      CD_komg_val = (2.0*node_val(density,node)*sigma_omg2*node_val(CD_komg,node)) / &
                    node_val(omega,node)

      a_3 = sqrt(node_val(k,node)) / (beta_star*node_val(omega,node)*node_val(y,node))
      a_4 = ( 500.0*node_val(bg_visc,1,1,node) )/( node_val(omega,node)*node_val(y,node)**2.0 )
      a_5 = ( 4.0*node_val(density,node)*sigma_omg2*node_val(k,node) ) / &
            ( max(node_val(CD_komg,node),1.0e-10)*node_val(y,node)**2.0 ) !A! Menter 1994 (Menter 2003 => e-10)
  
      arg_1 = min( max(a_3,a_4) , a_5 )
      arg_2 = max( 2.0*a_3, a_4 )

      F_1_val = tanh(arg_1**4.0)
      F_2_val = tanh(arg_2**2.0)

      ! set values of blending functions
      call set(CD_komg, node, CD_komg_val)
      call set(F_1, node, F_1_val)
      call set(F_2, node, F_2_val)

    end do node_loop

    call deallocate(k)
    call deallocate(omega)
    call deallocate(dummydensity)
    deallocate(dummydensity)
  end if

end subroutine komega_blending_functions

!------------------------------------------------------------------------------!

subroutine komega_calculate_rhs(state)

  type(state_type), intent(inout)  :: state

  type(scalar_field), dimension(2) :: src_abs_terms !A! 3 to 2
  type(scalar_field), dimension(2) :: fields
  type(scalar_field), pointer      :: src, abs, debug !A! f_1, f_2
  type(scalar_field)               :: src_to_abs
  type(vector_field), pointer      :: x, u, g
  type(scalar_field), pointer      :: dummydensity, density, scalar_eddy_visc !A! buoyancy_density
  integer :: i, ele, term, stat
  real    :: g_magnitude, alpha, beta, beta_star !A! c_eps_1, c_eps_2, sigma_p
  logical :: lump_mass !A! have_buoyancy_turbulence = .true.
  character(len=OPTION_PATH_LEN)              :: option_path 
  character(len=FIELD_NAME_LEN), dimension(2) :: field_names
  character(len=FIELD_NAME_LEN)               :: equation_type, implementation

  !A! ADM correction terms
  real    :: C_T, C_omega, beta_p, beta_d
  integer :: region_id_near_disc
  integer :: ndiscs
  integer, dimension(:), allocatable :: region_id_disc
  type(scalar_field), pointer :: NodeOnTurbine

  !A! SST
  logical :: have_SST
  real    :: gama_1, gama_2, beta_1, beta_2
  type(scalar_field), pointer :: F_1, F_2, CD_komg
  real                        :: sigma_omg2

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-omega/'

  if (.not. have_option(trim(option_path))) then 
     return
  end if

  ewrite(1,*) 'In calculate k-omega rhs' !!

  have_SST = have_option(trim(state%option_path)//"/subgridscale_parameterisations/k-omega/k-omega_SST")
  ewrite(1,*) 'k-omega SST? ', have_SST

  ! get model constants
  call get_option(trim(option_path)//'/Alpha', alpha, default = 5.0/9.0) !A! will 5/9 work?
  call get_option(trim(option_path)//'/Beta', beta, default = 0.075)
  call get_option(trim(option_path)//'/Beta_Star', beta_star, default = 0.09)

  !A! get ADM correction term constants
  call get_option('/ADM/C_T', C_T, default = 0.0)
  call get_option('/ADM/C_Omega', C_omega, default = 4.0) ! Rados et al (2008)
  call get_option('/ADM/Beta_p', beta_p, default = 0.05)  ! Rethore et al (2009)
  call get_option('/ADM/Beta_d', beta_d, default = 1.50)  ! Rethore et al (2009)

  call get_option('/ADM/NearDiscRegionID', region_id_near_disc, default = 902)

  call get_option('/ADM/NDiscs', ndiscs)
  allocate(region_id_disc(ndiscs))
  region_id_disc = 0
  call get_option('/ADM/DiscRegionID', region_id_disc)

  ! get field data
  x => extract_vector_field(state, "Coordinate")
  u => extract_vector_field(state, "NonlinearVelocity")
  scalar_eddy_visc => extract_scalar_field(state, "ScalarEddyViscosity")
  NodeOnTurbine => extract_scalar_field(state, "NodeOnTurbine") !A!

  !A! SST terms
  if(have_SST) then
    F_1 => extract_scalar_field(state, "F_1")
    F_2 => extract_scalar_field(state, "F_2")
    CD_komg => extract_scalar_field(state, "CD_komg")
    call get_option(trim(option_path)//'/Beta_1', beta_1, default = 0.0750) !A! same as beta
    call get_option(trim(option_path)//'/Beta_2', beta_2, default = 0.0828)
    call get_option(trim(option_path)//'/Gamma_1', gama_1, default = 0.5532) !A! c.f. 5/9
    call get_option(trim(option_path)//'/Gamma_2', gama_2, default = 0.4403)
    call get_option(trim(option_path)//'/Sigma_omg2', sigma_omg2, default = 0.85616)
  end if

!A!  g => extract_vector_field(state, "GravityDirection", stat)
!A!  call get_option('/physical_parameters/gravity/magnitude', g_magnitude, stat)
!A!  if (stat /= 0) then
!A!     have_buoyancy_turbulence = .false. !!! interesting
!A!  end if

  allocate(dummydensity)
  call allocate(dummydensity, X%mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
  call set(dummydensity, 1.0)
  dummydensity%option_path = ""
  
  ! Depending on the equation type, extract the density or set it to some dummy field allocated above
  call get_option(trim(state%option_path)//&
       "/vector_field::Velocity/prognostic/equation[0]/name", equation_type, stat=stat)
  if (stat /= 0) then
    density=>dummydensity
  else
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
  end if
  
!A!  if(have_buoyancy_turbulence) then
!A!     buoyancy_density => extract_scalar_field(state, 'VelocityBuoyancyDensity')
!A!  else
!A!     buoyancy_density => dummydensity
!A!  end if

  field_names(1) = 'TurbulentKineticEnergy'
  field_names(2) = 'TurbulentFrequency'

  !A! start with k and then do omega, why?
  field_loop: do i = 1, 2
     if (have_option(trim(option_path)//'scalar_field::'// &
          trim(field_names(i))//'/prescribed')) then !A! for debugging purposes
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
!A!     call allocate(src_abs_terms(3), fields(1)%mesh, name="buoyancy_term")
     call zero(src_abs_terms(1)); call zero(src_abs_terms(2)) !A! ; call zero(src_abs_terms(3))
     call zero(src); call zero(abs)

     !-----------------------------------------------------------------------------------

     !A! Assembly loop
     do ele = 1, ele_count(fields(1))

        call assemble_rhs_ele(src_abs_terms, fields(i), fields(3-i), scalar_eddy_visc, u, &
             density, g, g_magnitude, x, &
             alpha, beta, beta_star, ele, i, &
             C_T, C_omega, beta_p, beta_d, ndiscs, region_id_disc, region_id_near_disc, NodeOnTurbine, &
             have_SST, gama_1, gama_2, beta_1, beta_2, F_1, F_2, CD_komg) !A! buoyancy_density, have_buoyancy_turbulence
     end do

     ! For non-DG we apply inverse mass globally
     if(continuity(fields(1))>=0) then
        lump_mass = have_option(trim(option_path)//'mass_terms/lump_mass')
        do term = 1, 2 !A! 3 to 2
           call solve_cg_inv_mass(state, src_abs_terms(term), lump_mass, option_path)           
        end do
     end if
     !-----------------------------------------------------------------------------------

     ! Source disabling for debugging purposes
     do term = 1, 2 !!! 3 to 2
        if(have_option(trim(option_path)//'debugging_options/disable_'//&
             trim(src_abs_terms(term)%name))) then
           call zero(src_abs_terms(term))
        end if
     end do    
     !-----------------------------------------------------------------------------------

     ! Produce debugging output
     do term = 1, 2 !A! 3 to 2
        debug => extract_scalar_field(state, &
          trim(field_names(i))//"_"//trim(src_abs_terms(term)%name), stat)
        if (stat == 0) then
           call set(debug, src_abs_terms(term))
        end if
     end do
     !-----------------------------------------------------------------------------------
     
     ! Implement terms as source or absorbtion
     do term = 1, 2 !A! 3 to 2
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
           FLAbort("Unknown implementation type for k-omega source terms")
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
     do term = 1, 2 !A! 3 to 2
        call deallocate(src_abs_terms(term))
     end do
     call deallocate(fields(1))
     call deallocate(fields(2))

  end do field_loop

  call deallocate(dummydensity)
  deallocate(dummydensity)

end subroutine komega_calculate_rhs
    
!------------------------------------------------------------------------------!

subroutine assemble_rhs_ele(src_abs_terms, k, omega, scalar_eddy_visc, u, density, &
     g, g_magnitude, &
     X, alpha, beta, beta_star, ele, field_id, &
     C_T, C_omega, beta_p, beta_d, ndiscs, region_id_disc, region_id_near_disc, NodeOnTurbine, &
     have_SST, gama_1, gama_2, beta_1, beta_2, F_1, F_2, CD_komg) !A! have_buoyancy_turbulence, buoyancy_density

  type(scalar_field), dimension(2), intent(inout) :: src_abs_terms !A! 3 to 2
  type(scalar_field), intent(in) :: k, omega, scalar_eddy_visc
  type(vector_field), intent(in) :: X, u, g
  type(scalar_field), intent(in) :: density !A! buoyancy_density
  real, intent(in)               :: g_magnitude, alpha, beta, beta_star
!A!  logical, intent(in) :: have_buoyancy_turbulence
  integer, intent(in)            :: ele, field_id

  real, dimension(ele_loc(k, ele), ele_ngi(k, ele), x%dim) :: dshape
  real, dimension(ele_ngi(k, ele))    :: detwei, rhs, scalar_eddy_visc_ele, k_ele, omega_ele
  real, dimension(2, ele_loc(k, ele)) :: rhs_addto !A! 3 to 2
  integer, dimension(ele_loc(k, ele)) :: nodes
  real, dimension(ele_loc(k, ele), ele_loc(k, ele)) :: invmass
  real, dimension(u%dim, u%dim, ele_ngi(k, ele))    :: reynolds_stress, grad_u
  type(element_type), pointer                       :: shape
  integer :: term, ngi, dim, gi, i

  ! For buoyancy turbulence stuff
!A!  real, dimension(u%dim, ele_ngi(u, ele))  :: vector, u_quad, g_quad
!A!  real :: u_z, u_xy
!A!  real, dimension(ele_ngi(u, ele)) :: scalar, c_eps_3
!A!  type(element_type), pointer :: shape_density
!A!  real, dimension(:, :, :), allocatable :: dshape_density

  !A! For ADM correction terms
  type(scalar_field), intent(in)   :: NodeOnTurbine
  real, intent(in)                 :: C_T, C_omega, beta_p, beta_d
  integer                          :: discs
  integer, intent(in)              :: region_id_near_disc, ndiscs
  integer, dimension(ndiscs), intent(in)              :: region_id_disc
  real                             :: a_fac, C_x, u_ele, volume
  real, dimension(u%dim)           :: vel_integral
  real, dimension(ele_ngi(k, ele)) :: ADM_source

  !A! SST terms
  logical, intent(in)            :: have_SST
  real, intent(in)               :: gama_1, gama_2, beta_1, beta_2
  type(scalar_field), intent(in) :: F_1, F_2, CD_komg

  shape => ele_shape(k, ele)
  nodes = ele_nodes(k, ele)

  ! transforms triangular elements into standard isosceles triangles !?
  call transform_to_physical( X, ele, shape, dshape=dshape, detwei=detwei )

  ! get bounded values of k and omega for source terms
  ! this doesn't change the field values of k and omega
  k_ele = ele_val_at_quad(k,ele) !A! k value at quadrature points
  omega_ele = ele_val_at_quad(omega, ele)
  ngi = ele_ngi(u, ele) !A! ngi equals 11 => No. quadrature points
  do gi = 1, ngi
     k_ele(gi) = max(k_ele(gi), fields_min)
     omega_ele(gi) = max(omega_ele(gi), fields_min)
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

  !A! ADM: get u_ele for elements within the disc
  do discs = 1, ndiscs
    if (X%mesh%region_ids(ele)==region_id_disc(discs)) then
      volume=dot_product(ele_val_at_quad(NodeOnTurbine, ele), detwei)
      vel_integral=matmul(matmul(ele_val(u, ele), u%mesh%shape%n), detwei)
      u_ele = vel_integral(1)/volume
    end if
  end do

  ! Compute P (Production)
  if (field_id==1) then !A! id=1 => k
    rhs = tensor_inner_product(reynolds_stress, grad_u)
    if(have_SST) then !A! Menter k-omega SST (2003):
      rhs = min(rhs,10.0*omega_ele*k_ele*beta_star*ele_val_at_quad(density, ele))
    endif
    do discs = 1, ndiscs
      if (X%mesh%region_ids(ele)==region_id_disc(discs)) then ! Disc production term
       a_fac = 0.5*(1.0-sqrt(1.0-C_T))
       C_x = (4.0*a_fac)/(1.0 - a_fac)
       rhs = rhs +      0.5*C_x*(  beta_p*u_ele**3.0 -  beta_d*u_ele*k_ele) ! Rethore et al
      end if
    end do
  end if

  if (field_id==2) then !A! id=2 => omega
    rhs = tensor_inner_product(reynolds_stress, grad_u)
    if(have_SST) then !A! Menter k-omega SST (2003):
      rhs = min(rhs,10.0*omega_ele*k_ele*beta_star*ele_val_at_quad(density, ele))
    endif
    ADM_Source = (C_omega*rhs*rhs)/(k_ele*k_ele) !ADM!

    if(have_SST) then !A! Menter k-omega SST:

      rhs = rhs*( ele_val_at_quad(F_1, ele)*gama_1 + (1.0-ele_val_at_quad(F_1, ele))*gama_2 ) * &
                ( ele_val_at_quad(density, ele) / scalar_eddy_visc_ele ) + &
                (1.0-ele_val_at_quad(F_1, ele))*ele_val_at_quad(CD_komg, ele)

      do discs = 1, ndiscs
        if ((X%mesh%region_ids(ele)==region_id_disc(discs)) .or. (X%mesh%region_ids(ele)==region_id_near_disc)) then !ADM!
          rhs = rhs + ADM_Source
        endif
      end do

    else !A! Wilcox k-omega:

      rhs = rhs*alpha*(omega_ele/k_ele) !A! since eddy_visc = rho*(k/omega)
      do discs = 1, ndiscs
        if ((X%mesh%region_ids(ele)==region_id_disc(discs)) .or. (X%mesh%region_ids(ele)==region_id_near_disc)) then !ADM!
          rhs = rhs + ADM_Source
        end if
      end do

    endif
  end if
  rhs_addto(1,:) = shape_rhs(shape, detwei*rhs) !A! shape_rhs => integration

  ! Compute A (Absorption)
  rhs = -1.0*omega_ele*k_ele*beta_star*ele_val_at_quad(density, ele) !A! identical for k-omega and k-omega SST
  if (field_id==2) then
    if(have_SST) then
      rhs = -1.0*omega_ele*omega_ele*ele_val_at_quad(density, ele)*&
            (ele_val_at_quad(F_1, ele)*beta_1 + (1.0-ele_val_at_quad(F_1, ele))*beta_2)
    else
      rhs = -1.0*omega_ele*omega_ele*ele_val_at_quad(density, ele)*beta
    endif
  end if
  rhs_addto(2,:) = shape_rhs(shape, detwei*rhs)

  ! Gk:  
  ! Calculate buoyancy turbulence term and add to addto array
  ! No buoyancy term, so set this part of the array to zero. !A! Currently no buoyancy
!A!  rhs_addto(3,:) = 0.0    

  ! In the DG case we apply the inverse mass locally.
  if(continuity(k)<0) then
     invmass = inverse(shape_shape(shape, shape, detwei))
     do term = 1, 2 !A! 3 to 2
        rhs_addto(term,:) = matmul(rhs_addto(term,:), invmass)
     end do
  end if

  do term = 1, 2 !A! 3 to 2
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

!---------------------------------------------------------------------------------------!
! eddyvisc calculates the lengthscale and the eddy viscosity
! Eddy viscosity is added to the background viscosity.
!---------------------------------------------------------------------------------------!
subroutine komega_eddyvisc(state, advdif)

  type(state_type), intent(inout)  :: state
  logical, intent(in) :: advdif

  type(tensor_field), pointer      :: eddy_visc, viscosity, bg_visc
  type(vector_field), pointer      :: x, u
  type(scalar_field)               :: kk, omega ! uses kk instead of k
  type(scalar_field), pointer      :: scalar_eddy_visc, ll, density, dummydensity
  type(scalar_field)               :: ev_rhs
  integer                          :: i, j, ele, stat
  
  ! Options grabbed from the options tree
  character(len=OPTION_PATH_LEN)   :: option_path
  logical                          :: lump_mass, have_visc = .true.
  character(len=FIELD_NAME_LEN)    :: equation_type

  !A! SST terms
  logical                          :: have_SST
  type(scalar_field), pointer      :: F_2
  real                             :: beta_star

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-omega/'

  if (.not. have_option(trim(option_path))) then 
     return
  end if

  ewrite(1,*) "In komega_eddyvisc"

  have_SST = have_option(trim(state%option_path)//"/subgridscale_parameterisations/k-omega/k-omega_SST")
  ewrite(1,*) 'k-omega SST? ', have_SST

  ! Get field data
  call time_averaged_value(state, kk, "TurbulentKineticEnergy", advdif, option_path)
  call time_averaged_value(state, omega, "TurbulentFrequency", advdif, option_path)

  x                => extract_vector_field(state, "Coordinate")
  u                => extract_vector_field(state, "NonlinearVelocity")
  eddy_visc        => extract_tensor_field(state, "EddyViscosity")
  bg_visc          => extract_tensor_field(state, "BackgroundViscosity")
  scalar_eddy_visc => extract_scalar_field(state, "ScalarEddyViscosity")
  ll               => extract_scalar_field(state, "LengthScale")
  viscosity        => extract_tensor_field(state, "Viscosity", stat)
  if (stat /= 0) then
     have_visc = .false.
  end if

  call get_option(trim(option_path)//'/Beta_Star', beta_star, default = 0.09) !!

  !A! SST terms
  if(have_SST) then
    F_2       => extract_scalar_field(state, "F_2")
  end if

  ewrite_minmax(kk)
  ewrite_minmax(omega)
  ewrite_minmax(scalar_eddy_visc)
  
  allocate(dummydensity)
  call allocate(dummydensity, X%mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
  call set(dummydensity, 1.0)
  dummydensity%option_path = ""
  
  ! Depending on the equation type, extract the density or set it to some dummy field allocated above
  call get_option(trim(state%option_path)//&
       "/vector_field::Velocity/prognostic/equation[0]/name", equation_type, stat=stat)
  if (stat /= 0) then
    density=>dummydensity
  else
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
  end if

  call allocate(ev_rhs, scalar_eddy_visc%mesh, name="EVRHS")
  call zero(ev_rhs)

  ! Initialise viscosity to background value
  if (have_visc) then
     call set(viscosity, bg_visc)
  end if
  
  ! Compute the length scale diagnostic field here.
  do i = 1, node_count(scalar_eddy_visc)
     call set(ll, i, max(node_val(kk,i), fields_min)**0.5 / max(node_val(omega,i)*beta_star, fields_min))
  end do

  ! Calculate scalar eddy viscosity by integration over element
  do ele = 1, ele_count(scalar_eddy_visc)
     call komega_eddyvisc_ele(ele, X, kk, omega, scalar_eddy_visc, density, ev_rhs, &
                              F_2, u, have_SST)
  end do

  ! For non-DG we apply inverse mass globally
  if(continuity(scalar_eddy_visc)>=0) then
     lump_mass = have_option(trim(option_path)//'mass_terms/lump_mass')
     call solve_cg_inv_mass(state, ev_rhs, lump_mass, option_path)  
  end if
  
  ! Allow for prescribed eddy-viscosity
  if (.not. have_option(trim(option_path)//'/scalar_field::ScalarEddyViscosity/prescribed')) then
     call set(scalar_eddy_visc, ev_rhs)
  end if

  call deallocate(ev_rhs)
  call deallocate(kk)
  call deallocate(omega)
  
  call deallocate(dummydensity)
  deallocate(dummydensity)

  ewrite(2,*) "Setting k-omega eddy-viscosity tensor"
  call zero(eddy_visc)

  ! this is skipped if zero_eddy_viscosity is set - this is the easiest way to
  ! disable feedback from the k-omega model back into the rest of the model
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
  
   subroutine komega_eddyvisc_ele(ele, X, kk, omega, scalar_eddy_visc, density, ev_rhs, &
                                  F_2, u, have_SST)
   
      type(vector_field), intent(in)    :: x
      type(scalar_field), intent(in)    :: kk, omega, scalar_eddy_visc, density
      type(scalar_field), intent(inout) :: ev_rhs
      integer, intent(in)               :: ele
      
      type(element_type), pointer       :: shape_ev
      integer, pointer, dimension(:)    :: nodes_ev
      real, dimension(ele_ngi(scalar_eddy_visc, ele)) :: detwei
      real, dimension(ele_loc(scalar_eddy_visc, ele)) :: rhs_addto
      real, dimension(ele_loc(scalar_eddy_visc, ele), ele_loc(scalar_eddy_visc, ele)) :: invmass
      real, dimension(ele_ngi(kk, ele)) :: kk_at_quad, omega_at_quad

      !A! SST terms
      type(scalar_field), intent(in)    :: F_2
      real                              :: a_1
      logical, intent(in)               :: have_SST

      !A! invariant strain
      integer                                                    :: ngi, gi
      type(vector_field), intent(in)                             :: u
      type(element_type), pointer                                :: shape
      real, dimension(ele_loc(kk, ele), ele_ngi(kk, ele), x%dim) :: dshape
      real, dimension(u%dim, u%dim, ele_ngi(kk, ele))            :: grad_u, rate_of_strain
      real, dimension(ele_ngi(scalar_eddy_visc, ele))            :: strain_invariant

      ngi = ele_ngi(u, ele) !A! ngi equals 11 => No. quadrature points per element

      shape => ele_shape(kk, ele)
      call transform_to_physical( X, ele, shape, dshape=dshape, detwei=detwei )

      grad_u = ele_grad_at_quad(u, ele, dshape)
      do gi = 1, ngi
        rate_of_strain(:,:,gi) = 0.5*( grad_u(:,:,gi) + transpose(grad_u(:,:,gi)) )
      end do

      strain_invariant = sqrt( 2.0 * tensor_inner_product(rate_of_strain, rate_of_strain) )

      nodes_ev => ele_nodes(scalar_eddy_visc, ele)
      shape_ev =>  ele_shape(scalar_eddy_visc, ele)

      ! Get the k and omega values at the Gauss points
      kk_at_quad = ele_val_at_quad(kk,ele)
      omega_at_quad = ele_val_at_quad(omega,ele)
      
      ! Clip the field values at the Gauss points.
      ! Note 1: This isn't a permanent change directly to the field itself,
      ! only to the values used in the computation of the eddy viscosity.
      ! Note 2: Can't allow negative/zero omega or k.
      ! Note 3: Here we assume all fields have the same number of
      ! Gauss points per element.
      where (kk_at_quad < fields_min)
         kk_at_quad = fields_min
      end where
      where (omega_at_quad < fields_min)
         omega_at_quad = fields_min
      end where

      !A! Set eddy_visc depending on have_SST:
      if(have_SST) then
        call get_option(trim(option_path)//'/a_1', a_1, default = 0.31)
        rhs_addto = shape_rhs(shape_ev, detwei*ele_val_at_quad(density,ele)*&
                    a_1*kk_at_quad/&
                    max(a_1*omega_at_quad , strain_invariant*ele_val_at_quad(F_2,ele) ))
      else
        rhs_addto = shape_rhs(shape_ev, detwei*ele_val_at_quad(density,ele)*&
                    (kk_at_quad/omega_at_quad))
      endif

      ! In the DG case we will apply the inverse mass locally.
      if(continuity(scalar_eddy_visc)<0) then
         invmass = inverse(shape_shape(shape_ev, shape_ev, detwei))
         rhs_addto = matmul(rhs_addto, invmass)
      end if

      ! Add the element's contribution to the nodes of ev_rhs
      call addto(ev_rhs, nodes_ev, rhs_addto)    


   end subroutine komega_eddyvisc_ele

end subroutine komega_eddyvisc

!---------------------------------------------------------------------------------

subroutine komega_diffusion(state)

  ! calculates k and omega field diffusivities
  type(state_type), intent(inout)   :: state
  type(tensor_field), pointer       :: diff, bg_visc, eddy_visc
  real                              :: sigma, sigma_star
  integer                           :: i, j
  character(len=OPTION_PATH_LEN)    :: option_path

  !A! SST terms
  real                              :: sigma_k1, sigma_k2, sigma_omg1, sigma_omg2
  real                              :: sigma_k, sigma_omg
  type(scalar_field), pointer       :: F_1
  logical                           :: have_SST

  ewrite(1,*) 'in komega_diffusion'

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-omega/'

  have_SST = have_option(trim(state%option_path)//"/subgridscale_parameterisations/k-omega/k-omega_SST")
  ewrite(1,*) 'k-omega SST? ', have_SST

!  if (have_option(trim(option_path)//"/scalar_field::TurbulentKineticEnergy/prognostic")) then
!    if (.not.have_option(trim(option_path)//"/scalar_field::TurbulentKineticEnergy"//&
!         &"/prognostic/tensor_field::Diffusivity")) then
!      ewrite(1,*) 'NoDiffusivity'
!    else

  eddy_visc => extract_tensor_field(state, "EddyViscosity")
  bg_visc   => extract_tensor_field(state, "BackgroundViscosity")

  call get_option(trim(option_path)//'/Sigma', sigma, default = 0.5)
  call get_option(trim(option_path)//'/Sigma_Star', sigma_star, default = 0.5)

  !A! SST terms
  if(have_SST) then
    F_1       => extract_scalar_field(state, "F_1")
    call get_option(trim(option_path)//'/Sigma_k1', sigma_k1, default = 0.85034) !A! Wilcox: 0.50
    call get_option(trim(option_path)//'/Sigma_k2', sigma_k2, default = 1.0)
    call get_option(trim(option_path)//'/Sigma_omg1', sigma_omg1, default = 0.5)
    call get_option(trim(option_path)//'/Sigma_omg2', sigma_omg2, default = 0.85616)
  end if

  ! Set diffusivity
  diff => extract_tensor_field(state, "TurbulentKineticEnergyDiffusivity")
  call zero(diff)
  do i = 1, node_count(diff)
     do j = 1, diff%dim(1)
        call addto(diff, j, j, i, node_val(bg_visc, j, j, i)) ! why j j i ?
        if (have_SST) then
          sigma_k = F_1%val(i)*sigma_k1 + (1.0-F_1%val(i))*sigma_k2
          call addto(diff, j, j, i, node_val(eddy_visc, j, j, i) * sigma_k)
        else
          call addto(diff, j, j, i, node_val(eddy_visc, j, j, i) * sigma_star)
        endif
     end do
  end do
  diff => extract_tensor_field(state, "TurbulentFrequencyDiffusivity") !!
  call zero(diff)
  do i = 1, node_count(diff)
     do j = 1, diff%dim(1)
        call addto(diff, j, j, i, node_val(bg_visc, j, j, i))
        if (have_SST) then
          sigma_omg = F_1%val(i)*sigma_omg1 + (1.0-F_1%val(i))*sigma_omg2
          call addto(diff, j, j, i, node_val(eddy_visc, j, j, i) * sigma_omg)
        else
          call addto(diff, j, j, i, node_val(eddy_visc, j, j, i) * sigma)
        endif
     end do
  end do

!    end if !A!
!  end if !A!

end subroutine komega_diffusion

!---------------------------------------------------------------------------------

subroutine komega_tracer_diffusion(state)

  ! calculates scalar field diffusivity based upon eddy viscosity and background
  !  diffusivity
  type(state_type), intent(inout)   :: state

  type(tensor_field), pointer       :: t_field
  integer                           :: i_field, i, stat
  real                              :: local_background_diffusivity
  type(scalar_field)                :: local_background_diffusivity_field
  type(scalar_field), pointer       :: scalar_eddy_viscosity, s_field
  type(tensor_field), pointer       :: global_background_diffusivity
  type(tensor_field)                :: background_diffusivity

  ewrite(1,*) 'In komega_tracer_diffusion'

  do i_field = 1, scalar_field_count(state)
     s_field => extract_scalar_field(state, i_field)

     if (have_option(trim(s_field%option_path)//&
          '/prognostic/subgridscale_parameterisation::k-omega')) then

        ewrite(1,*) 'Calculating turbulent diffusivity for field: ', s_field%name
        
        ! check options
        if (.not.(have_option(trim(state%option_path)//'/subgridscale_parameterisations/k-omega')))&
             & then
           FLExit('you must have /subgridscale_parameterisations/k-omega to be able to calculate diffusivity based upon the k-omega model')
        end if

        t_field => extract_tensor_field(state, trim(s_field%name)//'Diffusivity', stat=stat) 
        if (stat /= 0) then
           FLExit('you must have a Diffusivity field to be able to calculate diffusivity based upon the k-omega model')
        else if (.not. have_option(trim(t_field%option_path)//"/diagnostic/algorithm::Internal")) then
           FLExit('you must have a diagnostic Diffusivity field with algorithm::Internal to be able to calculate diffusivity based upon the k-omega model') !!
        end if

        ! get sigma_p number !! what is this? buoyancy term?
        ! call get_option(trim(state%option_path)//'/subgridscale_parameterisations/k-omega/sigma_p', sigma_p)

        ! allocate and zero required fields
        call allocate(background_diffusivity, t_field%mesh, name="background_diffusivity")
        call zero(background_diffusivity)
        call allocate(local_background_diffusivity_field, t_field%mesh, &
             name="local_background_diffusivity_field")
        call zero(local_background_diffusivity_field)

        ! set background_diffusivity (local takes precendence over global)
        call get_option(trim(s_field%option_path)//&
             '/prognostic/subgridscale_parameterisation::k-omega/background_diffusivity', & !!
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
        !do i = 1, t_field%dim(1)
        !   call addto(t_field, i, i, scalar_eddy_viscosity, 1.0/sigma_p)
        !end do

        call deallocate(background_diffusivity)
        call deallocate(local_background_diffusivity_field)

     end if
  end do

end subroutine komega_tracer_diffusion

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
          'mass_terms/use_consistent_mass_matrix/')
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
          'mass_terms/use_consistent_mass_matrix/')
     call set(A, x)
     call deallocate(x)
  end if

end subroutine solve_cg_inv_mass_vector

!---------------------------------------------------------------------------------

subroutine k_omega_check_options

  character(len=OPTION_PATH_LEN) :: option_path
  character(len=FIELD_NAME_LEN)  :: kmsh, emsh, vmsh
  integer                        :: dimension, stat, n_phases, istate

  ewrite(1,*) "In komega_check_options"

  n_phases = option_count("/material_phase")
  
  if(option_count("/material_phase/subgridscale_parameterisations/k-omega") > 1) then
     FLExit("The k-omega model can only be applied to a single-phase.")
  end if

  do istate = 0, n_phases-1

     option_path = "/material_phase["//int2str(istate)//"]/subgridscale_parameterisations/k-omega"
     
     ! one dimensional problems not supported
     call get_option("/geometry/dimension/", dimension) 
     if (dimension .eq. 1 .and. have_option(trim(option_path))) then
        FLExit("k-omega model is only supported for dimension > 1")
     end if
     ! Don't do k-omega if it's not included in the model!
     if (.not.have_option(trim(option_path))) return

     ! checking for required fields
     if (have_option(trim(option_path)//"/scalar_field::TurbulentKineticEnergy/prognostic")) then
        ! diffusivity is on and diagnostic
        if (.not.have_option(trim(option_path)//"/scalar_field::TurbulentKineticEnergy"//&
             &"/prognostic/tensor_field::Diffusivity")) then
           FLExit("You need TurbulentKineticEnergy Diffusivity field for k-omega")
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
           FLExit("You need TurbulentKineticEnergy Source field for k-omega")
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
           FLExit("You need TurbulentKineticEnergy Absorption field for k-omega")
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
        FLExit("You need prognostic/prescribed TurbulentKineticEnergy field for k-omega")
     end if
     if (have_option(trim(option_path)//"/scalar_field::TurbulentFrequency/prognostic")) then
        ! diffusivity is on and diagnostic
        if (.not.have_option(trim(option_path)//"/scalar_field::TurbulentFrequency"//&
             &"/prognostic/tensor_field::Diffusivity")) then
           FLExit("You need TurbulentFrequency Diffusivity field for k-omega")
        end if
        if (.not.have_option(trim(option_path)//&
             &"/scalar_field::TurbulentFrequency/prognostic/"//&
             &"/tensor_field::Diffusivity/diagnostic/algorithm::Internal")) then
           FLExit("You need TurbulentFrequency Diffusivity field set to diagnostic/internal")
        end if
        ! source terms
        if (.not.have_option(trim(option_path)//&
             &"/scalar_field::TurbulentFrequency/prognostic"//&
             &"/scalar_field::Source")) then
           FLExit("You need TurbulentFrequency Source field for k-omega")
        end if
        if (.not. have_option(trim(option_path)//&
             &"/scalar_field::TurbulentFrequency/prognostic"//&
             &"/scalar_field::Source/diagnostic/algorithm::Internal")) then
           FLExit("You need TurbulentFrequency Source field set to diagnostic/internal")
        end if
        ! absorption terms
        if (.not.have_option(trim(option_path)//&
             &"/scalar_field::TurbulentFrequency/prognostic"//&
             &"/scalar_field::Absorption")) then
           FLExit("You need TurbulentFrequency Absorption field for k-omega")
        end if
        if (.not.have_option(trim(option_path)//&
             &"/scalar_field::TurbulentFrequency/prognostic"//&
             &"/scalar_field::Absorption/diagnostic/algorithm::Internal")) then
           FLExit("You need TurbulentFrequency Absorption field set to diagnostic/internal")
        end if
     else if (have_option(trim(option_path)// &
          "/scalar_field::TurbulentFrequency/prescribed")) then
        ewrite(0,*) "WARNING: TurbulentFrequency field is prescribed"
     else
        FLExit("You need prognostic/prescribed TurbulentFrequency field for k-omega")
     end if

     ! Check that TurbulentKineticEnergy and TurbulentFrequency fields are on the same 
     !  mesh as the velocity
     call get_option(trim(option_path)//&
          &"/scalar_field::TurbulentKineticEnergy/prognostic/mesh/name", kmsh, stat)
     if (stat /= 0) then
        call get_option(trim(option_path)//&
             &"/scalar_field::TurbulentKineticEnergy/prescribed/mesh/name", kmsh,&
             & stat)
     end if
     call get_option(trim(option_path)//&
          &"/scalar_field::TurbulentFrequency/prognostic/mesh/name", emsh, stat) 
     if (stat /= 0) then
        call get_option(trim(option_path)//&
             &"/scalar_field::TurbulentFrequency/prescribed/mesh/name", emsh,& 
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
!A!     if(.not. kmsh==emsh .or. .not. kmsh==vmsh .or. .not. emsh==vmsh) then
!A!        FLExit("You must use the Velocity mesh for TurbulentKineticEnergy and TurbulentFrequency fields") 
!A!     end if

     ! Velocity field options
     if (.not.have_option("/material_phase["//int2str(istate)//"]/vector_field::Velocity/prognostic"//&
          "/tensor_field::Viscosity/") .and. &
          .not.have_option("/material_phase["//int2str(istate)//"]/vector_field::Velocity/prescribed")) then
        FLExit("Need viscosity switched on under the Velocity field for k-omega.") 
     end if
     ! check that the user has switched Velocity/viscosity to diagnostic
     if (have_option("/material_phase["//int2str(istate)//"]/vector_field::Velocity/prognostic") .and. &
          .not.have_option("/material_phase["//int2str(istate)//"]/vector_field::Velocity/prognostic"//&
          "/tensor_field::Viscosity/diagnostic/")) then
        FLExit("You need to switch the viscosity field under Velocity to diagnostic/internal")
     end if

     ! Check ScalarEddyViscosity is diagnostic
     if (have_option(trim(option_path)//'/scalar_field::ScalarEddyViscosity/prescribed')) then
        ewrite(0,*) "WARNING: ScalarEddyViscosity field is prescribed"
     end if

  end do

end subroutine k_omega_check_options 

end module k_omega !!
