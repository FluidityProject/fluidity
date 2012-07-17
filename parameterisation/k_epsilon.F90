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
  use initialise_fields_module

implicit none

  private

  ! locally allocatad fields
  real, save                          :: fields_min = 1.e-10

  public :: keps_diagnostics, keps_eddyvisc, keps_bcs, keps_adapt_mesh,&
       & k_epsilon_check_options, keps_momentum_source

  ! Outline:
  !  - call diagnostics to obtain source terms and calculate eddy viscosity
  !  - after each scalar field solve recalculates the eddy viscosity
  !  - wall functions are added to selected boundaries in keps_bcs and wall_functions
  !  - keps_adapt_options repopulates the fields after an adapt

contains

subroutine keps_diagnostics(state)

  type(state_type), intent(inout) :: state

  call keps_eddyvisc(state)
  call keps_calculate_rhs(state)

end subroutine keps_diagnostics

subroutine keps_calculate_rhs(state)

  type(state_type), intent(inout) :: state

  type(scalar_field), dimension(3) :: src_abs_terms
  type(scalar_field_pointer), dimension(2) :: fields
  type(scalar_field), pointer :: src, abs, EV, f_1, f_2, debug
  type(scalar_field_pointer), allocatable, dimension(:) :: buoyant_fields
  type(scalar_field) :: src_to_abs
  type(vector_field), pointer :: positions, u, g
  type(tensor_field), pointer :: diff
  integer :: i, node, ele, term, stat
  real :: g_magnitude, c_eps_1, c_eps_2, sigma_eps, sigma_k
  real, allocatable, dimension(:) :: beta, delta_t
  logical :: prescribed, gravity = .true., have_buoyant_fields =.false., lump_mass
  character(len=OPTION_PATH_LEN) :: option_path 
  character(len=FIELD_NAME_LEN), dimension(2) :: field_names

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-epsilon/'

  if (.not. have_option(trim(option_path))) then 
     return
  end if

  ewrite(1,*) 'In calculate k-epsilon rhs'

  ! get model constants
  call get_option(trim(option_path)//'/C_eps_1', c_eps_1, default = 1.44)
  call get_option(trim(option_path)//'/C_eps_2', c_eps_2, default = 1.92)
  call get_option(trim(option_path)//'/sigma_k', sigma_k, default = 1.0)
  call get_option(trim(option_path)//'/sigma_eps', sigma_eps, default = 1.3)
  
  ! get field data
  positions    => extract_vector_field(state, "Coordinate")
  u            => extract_vector_field(state, "Velocity")
  EV           => extract_scalar_field(state, "ScalarEddyViscosity")
  f_1          => extract_scalar_field(state, "f_1")
  f_2          => extract_scalar_field(state, "f_2")
  g            => extract_vector_field(state, "GravityDirection", stat)
  call get_option('physical_parameters/gravity/magnitude', g_magnitude, stat)
  if (stat /= 0) then
     gravity = .false.
  end if
  call get_scalar_field_buoyancy_data(state, buoyant_fields, beta, delta_t)
  if (allocated(buoyant_fields)) then
     have_buoyant_fields = .true.
  end if

  field_names(1) = 'TurbulentKineticEnergy'
  field_names(2) = 'TurbulentDissipation'

  field_loop: do i = 1, 2
     !-----------------------------------------------------------------------------------
     
     ! Setup

     src            => extract_scalar_field(state, trim(field_names(i))//"Source")
     abs            => extract_scalar_field(state, trim(field_names(i))//"Absorption")
     fields(1)%ptr  => extract_scalar_field(state, trim(field_names(i)))
     fields(2)%ptr  => extract_scalar_field(state, trim(field_names(3-i)))

     call allocate(src_abs_terms(1), fields(1)%ptr%mesh, name="Production")
     call allocate(src_abs_terms(2), fields(1)%ptr%mesh, name="Destruction")
     call allocate(src_abs_terms(3), fields(1)%ptr%mesh, name="Buoyancy")
     call zero(src_abs_terms(1)); call zero(src_abs_terms(2)); call zero(src_abs_terms(3))
     call zero(src); call zero(abs)
     !-----------------------------------------------------------------------------------

     ! Assembly loop
     
     do ele = 1, ele_count(fields(1)%ptr)
        call assemble_rhs_ele(src_abs_terms, fields(i)%ptr, fields(3-i)%ptr, EV, u, buoyant_fields, &
             & have_buoyant_fields, beta, delta_t, g, g_magnitude, gravity, positions,&
             & c_eps_1, c_eps_2, sigma_k, sigma_eps, f_1, f_2, ele, i)
     end do

     ! For non-DG we apply inverse mass globally
     if(continuity(fields(1)%ptr)>=0) then
        lump_mass = have_option(trim(option_path)//'mass_lumping_in_diagnostics/lump_mass')
        do term = 1, 3
           call solve_cg_inv_mass(state, src_abs_terms(term), lump_mass, option_path)           
        end do
     end if
     !-----------------------------------------------------------------------------------

     ! Produce debugging output

     debug => extract_scalar_field(state, &
          trim(field_names(i))//"Production", stat)
     if (stat == 0) then
        call set(debug, src_abs_terms(1))
     end if
     debug => extract_scalar_field(state, &
          trim(field_names(i))//"Destruction", stat)
     if (stat == 0) then
        call set(debug, src_abs_terms(2))
     end if
     debug => extract_scalar_field(state, &
          trim(field_names(i))//"BuoyancyTerm", stat)
     if (stat == 0) then
        call set(debug, src_abs_terms(3))
     end if
     !-----------------------------------------------------------------------------------
     
     ! This allows user-specified implicit/explicit rhs terms

     ewrite(2,*) "Calculating k source and absorption"
     ! Set implicit/explicit source/absorbtion terms
     if(have_option(trim(option_path)//'implicit_source')) then
        call allocate(src_to_abs, fields(1)%ptr%mesh, name='SourceToAbsorbtion')
        call set(src_to_abs, fields(1)%ptr)
        where (src_to_abs%val >= fields_min)
           src_to_abs%val=1./src_to_abs%val
        elsewhere
           src_to_abs%val=1./fields_min
        end where
        call scale(src_abs_terms(1), src_to_abs)
        call addto(abs, src_abs_terms(1), -1.0)
        call deallocate(src_to_abs)
     else
        call addto(src, src_abs_terms(1))
     end if
     if(have_option(trim(option_path)//'implicit_buoyancy')) then
        call allocate(src_to_abs, fields(1)%ptr%mesh, name='SourceToAbsorbtion')
        call set(src_to_abs, fields(1)%ptr)
        where (src_to_abs%val >= fields_min)
           src_to_abs%val = 1./src_to_abs%val
        elsewhere
           src_to_abs%val = 1./fields_min
        end where
        call scale(src_abs_terms(3), src_to_abs)
        call addto(abs, src_abs_terms(3), -1.0)
        call deallocate(src_to_abs)
     else
        call addto(src, src_abs_terms(3))
     end if
     if(have_option(trim(option_path)//'explicit_absorbtion')) then
        call scale(src_abs_terms(2), fields(1)%ptr)
        call addto(src, src_abs_terms(2), -1.0)
     else
        call addto(abs, src_abs_terms(2))
     end if
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

  end do field_loop

  call deallocate_scalar_field_buoyancy_data(buoyant_fields, beta, delta_t)

  ! Set diffusivity
  diff => extract_tensor_field(state, trim(field_names(1))//"Diffusivity")
  do i = 1, diff%dim(1)
     call set(diff, i, i, EV, scale=1. / sigma_k)
  end do
  diff => extract_tensor_field(state, trim(field_names(2))//"Diffusivity")
  do i = 1, diff%dim(1)
     call set(diff, i, i, EV, scale=1. / sigma_eps)
  end do

end subroutine keps_calculate_rhs
    
!------------------------------------------------------------------------------!
! calculate the source and absorbtion terms                                    !
!------------------------------------------------------------------------------!
subroutine assemble_rhs_ele(src_abs_terms, k, eps, EV, u, buoyant_fields, &
     have_buoyant_fields, beta, delta_t, g, g_magnitude, gravity, positions, &
     c_eps_1, c_eps_2, sigma_k, sigma_eps, f_1, f_2, ele, field_id)

  type(scalar_field), dimension(3), intent(inout) :: src_abs_terms
  type(scalar_field), intent(in) :: k, eps, EV, f_1, f_2
  type(scalar_field_pointer), dimension(:) :: buoyant_fields
  type(vector_field), intent(in) :: positions, u, g
  real, dimension(:) :: beta, delta_t
  real, intent(in) :: g_magnitude, c_eps_1, c_eps_2, sigma_k, sigma_eps
  logical, intent(in) :: gravity, have_buoyant_fields
  integer, intent(in) :: ele, field_id

  real, dimension(ele_loc(k, ele), ele_ngi(k, ele), positions%dim) :: dshape
  real, dimension(ele_ngi(k, ele)) :: detwei, strain_ngi, inv_k, rhs
  real, dimension(3, ele_loc(k, ele)) :: rhs_addto
  integer, dimension(ele_loc(k, ele)) :: nodes
  real, dimension(ele_loc(k, ele), ele_loc(k, ele)) :: invmass
  real, dimension(u%dim, u%dim, ele_ngi(k, ele)) :: rhs_tensor
  type(element_type), pointer :: shape
  integer :: i_loc, i_field, term

  shape => ele_shape(k, ele)
  nodes = ele_nodes(k, ele)

  call transform_to_physical( positions, ele, shape, dshape=dshape, detwei=detwei )

  ! require inverse of k field for some terms
  inv_k = ele_val_at_quad(k,ele)
  where (inv_k >= fields_min)
     inv_k = 1.0/inv_k
  elsewhere
     inv_k = 1.0/fields_min
  end where

  ! P:
  rhs_tensor = reynolds_stresses(u, EV, dshape, ele)
  rhs = tensor_inner_product(rhs_tensor, ele_grad_at_quad(u, ele, dshape))
  if (field_id==2) then
     rhs = rhs*c_eps_1*ele_val_at_quad(f_1,ele)*ele_val_at_quad(eps,ele)*inv_k
  end if
  rhs_addto(1,:) = shape_rhs(shape, detwei*rhs)

  ! A:
  rhs = ele_val_at_quad(eps,ele)*inv_k
  if (field_id==2) then
     rhs = rhs*c_eps_2*ele_val_at_quad(f_2,ele)
  end if
  rhs_addto(2,:) = shape_rhs(shape, detwei*rhs)

  ! Gk:
  ! zero buoyancy addto array
  do i_loc = 1, ele_loc(k, ele)
     rhs_addto(3,i_loc) = 0.0
  end do
  ! check buoyant fields are possible and present
  if (have_buoyant_fields .and. gravity) then 
     ! loop through buoyant fields, calculate source term and add to addto array
     do i_field = 1, size(buoyant_fields, 1)
        call calculate_buoyancy_term(rhs_addto, EV, u, buoyant_fields(i_field)%ptr, &
             & beta(i_field), delta_t(i_field), g, g_magnitude, positions, ele, detwei,&
             & shape, field_id, c_eps_1, f_1)
     end do
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
! calculate reynolds stresses   = EV*symmetric gradient                        !
!------------------------------------------------------------------------------!
function reynolds_stresses(u, EV, dshape, ele)

  type(vector_field), intent(in) :: u
  type(scalar_field), intent(in) :: EV
  integer, intent(in) :: ele
  real, dimension(:, :, :), intent(in) :: dshape
  
  real, dimension(u%dim, u%dim, ele_ngi(u, ele)) :: reynolds_stresses
  real, dimension(u%dim, u%dim, ele_ngi(u, ele)) :: grad_u
  real, dimension(ele_ngi(u, ele)) :: EV_ele
  integer :: ngi, dim, gi, i

  grad_u = ele_grad_at_quad(u, ele, dshape)
  EV_ele = ele_val_at_quad(EV, ele)
  
  ngi = ele_ngi(u, ele)
  dim = u%dim

  do gi = 1, ngi
     reynolds_stresses(:,:,gi) = EV_ele(gi)*(grad_u(:,:,gi) + transpose(grad_u(:,:,gi)))
  end do

end function reynolds_stresses

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
    
!------------------------------------------------------------------------------!
! calculate the buoyancy source term for each buoyant scalar field             !
!------------------------------------------------------------------------------!
subroutine calculate_buoyancy_term(rhs_addto, EV, u, buoyant_field, beta,&
     & delta_t, g, g_magnitude, positions, ele, detwei, shape, field_id, c_eps_1, f_1)

  real, intent(inout), dimension(:,:) :: rhs_addto
  type(scalar_field), intent(in) :: EV, buoyant_field, f_1 
  type(vector_field), intent(in) :: positions, u, g
  real, intent(in), dimension(:) :: detwei
  real, intent(in) :: g_magnitude, c_eps_1, beta, delta_t
  type(element_type), intent(in), pointer :: shape
  integer, intent(in) :: ele, field_id

  real, dimension(u%dim, ele_ngi(u, ele))  :: vector, u_z, u_xy
  real, dimension(ele_ngi(u, ele)) :: scalar, c_eps_3
  real, dimension(ele_loc(buoyant_field, ele),ele_ngi(buoyant_field, ele),positions%dim) :: dshape_s
  type(element_type), pointer :: shape_s

  integer :: i_gi, i_loc, i_field

  ! get dshape for scalar field so that we can obtain gradients 
  shape_s => ele_shape(buoyant_field, ele)
  call transform_to_physical( positions, ele, shape_s, dshape=dshape_s )        

  ! calculate scalar and vector components of the source term
  vector = ele_val_at_quad(g, ele)*ele_grad_at_quad(buoyant_field,&
       & ele, dshape_s)
  scalar = -1.0*beta*g_magnitude*ele_val_at_quad(EV, ele)/delta_t

  ! multiply vector component by scalar and sum across dimensions - note that the
  ! vector part has been multiplied by the gravitational direction so the it is
  ! zero everywhere apart from in this direction
  do i_gi = 1, ele_ngi(u, ele)
     scalar(i_gi) = sum(scalar(i_gi) * vector(:, i_gi))
  end do
  
  if (field_id==2) then  
     ! get components of velocity in direction of gravity and in other directions
     u_z = abs(ele_val_at_quad(g, ele)) * ele_val_at_quad(u, ele)
     u_xy = ele_val_at_quad(u, ele) - u_z
     ! calculate c_eps_3 = tanh(v/u)
     do i_gi = 1, ele_ngi(u, ele)
        if (norm2(u_xy(:, i_gi)) /= 0.0) then
           c_eps_3(i_gi) = tanh(norm2(u_z(:, i_gi))/norm2(u_xy(:, i_gi))) 
        else
           c_eps_3(i_gi) = 1.0
        end if
     end do
     scalar = scalar*c_eps_1*ele_val_at_quad(f_1,ele)*c_eps_3
  end if

  ! multiply by determinate weights, integrate and add to rhs
  rhs_addto(3,:) = rhs_addto(3,:) + shape_rhs(shape, scalar * detwei)     

end subroutine calculate_buoyancy_term

!----------
! eddyvisc calculates the lengthscale, and then the eddy viscosity.!
! Eddy viscosity is added to the background viscosity.
!----------

subroutine keps_eddyvisc(state)

  type(state_type), intent(inout)  :: state
  type(tensor_field), pointer      :: eddy_visc, viscosity, bg_visc
  type(vector_field), pointer      :: positions
  type(scalar_field), pointer      :: kk, eps, EV, ll, f_mu
  type(scalar_field)               :: ev_rhs
  type(element_type), pointer      :: shape_ev
  integer                          :: i, ele, stat
  integer, pointer, dimension(:)   :: nodes_ev
  real, allocatable, dimension(:)  :: detwei, rhs_addto
  real, allocatable, dimension(:,:):: invmass
  real                             :: c_mu, lmax
  character(len=OPTION_PATH_LEN)   :: option_path
  logical                          :: lump_mass, have_visc = .true.

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-epsilon/'

  if (.not. have_option(trim(option_path))) then 
     return
  end if

  ewrite(1,*) "In keps_eddyvisc"

  ! get model constants
  call get_option(trim(option_path)//'/lengthscale_limit', lmax)
  call get_option(trim(option_path)//'/C_mu', c_mu, default = 0.09)
  
  ! get field data
  kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
  eps        => extract_scalar_field(state, "TurbulentDissipation")
  positions  => extract_vector_field(state, "Coordinate")
  eddy_visc  => extract_tensor_field(state, "EddyViscosity")
  f_mu       => extract_scalar_field(state, "f_mu")
  bg_visc    => extract_tensor_field(state, "BackgroundViscosity")
  EV         => extract_scalar_field(state, "ScalarEddyViscosity")
  ll         => extract_scalar_field(state, "LengthScale")
  viscosity  => extract_tensor_field(state, "Viscosity", stat)
  if (stat /= 0) then
     have_visc = .false.
  end if

  ewrite_minmax(kk)
  ewrite_minmax(eps)
  ewrite_minmax(EV)

  call allocate(ev_rhs, EV%mesh, name="EVRHS")
  call zero(ev_rhs)

  ! Initialise viscosity to background value
  if (have_visc) then
     call set(viscosity, bg_visc)
  end if

  !Clip fields: can't allow negative/zero epsilon or k
  do i = 1, node_count(EV)
     call set(kk, i, max(node_val(kk,i), fields_min))
     call set(eps, i, max(node_val(eps,i), fields_min))
     ! Limit lengthscale to prevent instablilities.
     call set(ll, i, min(node_val(kk,i)**1.5 / node_val(eps,i), lmax))
  end do

  ! Calculate scalar eddy viscosity by integration over element
  do ele = 1, ele_count(EV)
     nodes_ev => ele_nodes(EV, ele)
     shape_ev =>  ele_shape(EV, ele)
     allocate(detwei (ele_ngi(EV, ele)))
     allocate(rhs_addto (ele_loc(EV, ele)))
     allocate(invmass (ele_loc(EV, ele), ele_loc(EV, ele)))
     call transform_to_physical(positions, ele, detwei=detwei)
     rhs_addto = shape_rhs(shape_ev, detwei*C_mu*ele_val_at_quad(f_mu,ele)*ele_val_at_quad(kk,ele)**0.5*ele_val_at_quad(ll,ele))
     ! In the DG case we will apply the inverse mass locally.
     if(continuity(EV)<0) then
        invmass = inverse(shape_shape(shape_ev, shape_ev, detwei))
        rhs_addto = matmul(rhs_addto, invmass)
     end if
     call addto(ev_rhs, nodes_ev, rhs_addto)
     deallocate(detwei, rhs_addto, invmass)
  end do

  ! For non-DG we apply inverse mass globally
  if(continuity(EV)>=0) then
     lump_mass = have_option(trim(option_path)//'mass_lumping_in_diagnostics/lump_mass')
     call solve_cg_inv_mass(state, ev_rhs, lump_mass, option_path)  
  end if

  call set(EV, ev_rhs)
  call deallocate(ev_rhs)

  ewrite(2,*) "Setting k-epsilon eddy-viscosity tensor"
  call zero(eddy_visc)

  ! eddy viscosity tensor is isotropic
  ! this is skipped if zero_eddy_viscosity is set - this is the easiest way to
  ! disable feedback from the k-epsilon model back into the rest of the model
  if (.not. have_option(trim(option_path)//'debugging_options/zero_reynolds_stress_tensor')) then
     do i = 1, eddy_visc%dim(1)
        call set(eddy_visc, i, i, EV)
     end do
  end if

  ! Add turbulence model contribution to viscosity field
  if (have_visc) then
     call addto(viscosity, eddy_visc)
  end if

  ewrite_minmax(eddy_visc)
  ewrite_minmax(viscosity)
  ewrite_minmax(kk)
  ewrite_minmax(eps)
  ewrite_minmax(EV)
  ewrite_minmax(ll)

end subroutine keps_eddyvisc

!--------------------------------------------------------------------------------!
! calculates the reynolds stress tensor correction term 
! - grad(2/3 k delta(ij)) 
! this is added to prescribed momentum source fields or sets diagnostic source 
! fields 
!--------------------------------------------------------------------------------!
subroutine keps_momentum_source(state)

  type(state_type), intent(inout)  :: state
  
  type(scalar_field), pointer :: k, lumped_mass
  type(vector_field), pointer :: source, x
  type(vector_field) :: prescribed_source
  integer :: i
  logical :: prescribed, lump_mass
  character(len=OPTION_PATH_LEN) :: option_path 

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-epsilon/'

  if (.not.have_option(trim(option_path)) .or. &
      have_option(trim(option_path)//'debugging_options/zero_reynolds_stress_tensor') .or. &
      have_option("/material_phase[0]/vector_field::Velocity/prescribed")) then 
     return
  end if

  k         => extract_scalar_field(state, "TurbulentKineticEnergy")
  source    => extract_vector_field(state, "VelocitySource")
  x         => extract_vector_field(state, "Coordinate")

  call zero(source)
  do i = 1, ele_count(k)
     call keps_momentum_source_ele()
  end do

  ! For non-DG we apply inverse mass globally
  if(continuity(k)>=0) then
     lump_mass = have_option(trim(option_path)//&
          'mass_lumping_in_diagnostics/lump_mass')
     call solve_cg_inv_mass_vector(state, source, lump_mass, option_path)  
  end if

  ! Allow for prescribed momentum source
  prescribed = (have_option(trim(source%option_path)//'/prescribed/'))
  if(prescribed) then
     call allocate(prescribed_source, source%dim, source%mesh, name='PrescribedSource')
     call initialise_field_over_regions(prescribed_source, &
          trim(source%option_path)//'/prescribed/value', &
          x)
     call addto(source, prescribed_source)
     call deallocate(prescribed_source)
  end if

contains
    
  subroutine keps_momentum_source_ele()
      
    real, dimension(ele_loc(k, i), ele_ngi(k, i), x%dim) :: dshape
    real, dimension(ele_ngi(k, i)) :: detwei
    integer, dimension(ele_loc(k, i)) :: nodes
    real, dimension(ele_loc(k, i), ele_loc(k, i)) :: invmass
    type(element_type), pointer :: shape
    real, dimension(ele_ngi(source, i)) :: rhs
    real, dimension(x%dim, ele_ngi(k, i)) :: grad_k
    integer :: j

    shape => ele_shape(k, i)
    nodes = ele_nodes(source, i)

    call transform_to_physical( x, i, shape, dshape=dshape, detwei=detwei )

    grad_k = ele_grad_at_quad(k, i, dshape)
    do j = 1, x%dim
       rhs = shape_rhs(shape, -(2./3.)*detwei*grad_k(j,:))
       
       ! In the DG case we apply the inverse mass locally.
       if(continuity(k)<0) then
          invmass = inverse(shape_shape(shape, shape, detwei))
          rhs = matmul(rhs, invmass)
       end if

       call addto(source, j, nodes, rhs)  
    end do

  end subroutine keps_momentum_source_ele

end subroutine keps_momentum_source

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
          'mass_lumping_in_diagnostics/solve_using_mass_matrix/')
     call set(A, x)
     call deallocate(x)
  end if

end subroutine solve_cg_inv_mass

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
          'mass_lumping_in_diagnostics/solve_using_mass_matrix/')
     call set(A, x)
     call deallocate(x)
  end if

end subroutine solve_cg_inv_mass_vector

!--------------------------------------------------------------------------------!
! This gets and applies locally defined boundary conditions (wall functions)     !
!--------------------------------------------------------------------------------!

subroutine keps_bcs(state)

  type(state_type), intent(in)               :: state
  type(scalar_field), pointer                :: field1, field2    ! k or epsilon
  type(scalar_field), pointer                :: f_1, f_2, f_mu, y
  type(scalar_field), pointer                :: surface_field, EV
  type(vector_field), pointer                :: positions, u
  type(tensor_field), pointer                :: bg_visc
  type(scalar_field)                         :: rhs_field, surface_values
  type(mesh_type), pointer                   :: surface_mesh
  integer                                    :: i, j, ele, sele, index, nbcs, stat, node
  integer, dimension(:), pointer             :: surface_elements, surface_node_list
  character(len=FIELD_NAME_LEN)              :: bc_type, bc_name, wall_fns
  character(len=OPTION_PATH_LEN)             :: bc_path, bc_path_i, option_path 
  real                                       :: c_mu

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-epsilon/'

  ewrite(2,*) "In keps_bcs"

  positions => extract_vector_field(state, "Coordinate")
  u         => extract_vector_field(state, "Velocity")
  EV        => extract_scalar_field(state, "ScalarEddyViscosity")
  bg_visc   => extract_tensor_field(state, "BackgroundViscosity")
  f_1       => extract_scalar_field(state, "f_1")
  f_2       => extract_scalar_field(state, "f_2")
  f_mu      => extract_scalar_field(state, "f_mu")

  ! initialise low_Re damping functions
  call set(f_1, 1.0)
  call set(f_2, 1.0)
  call set(f_mu, 1.0)

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

        if (trim(bc_type) .eq. "k_epsilon") then
           ! Get bc by name. Get type just to make sure it's now dirichlet
           call get_boundary_condition(field1, name=bc_name, type=bc_type, surface_node_list=surface_node_list, &
                surface_element_list=surface_elements, surface_mesh=surface_mesh)

           ! Do we have high- or low-Reynolds-number wall functions?
           call get_option(trim(bc_path_i)//"/type::k_epsilon/", wall_fns)

           ewrite(2,*) "Calculating field BC: ",trim(field1%name),' ',trim(bc_name),' ',trim(bc_type),' ',trim(wall_fns)

           ! Get surface field already created in bcs_from_options
           surface_field => extract_surface_field(field1, bc_name=bc_name, name="value")
           call zero(surface_field)

           call allocate(surface_values, surface_mesh, name="surfacevalues")
           call allocate(rhs_field, field1%mesh, name="rhs")
           call zero(surface_values); call zero(rhs_field)

           do j = 1, ele_count(surface_mesh)
              sele = surface_elements(j)
              ele  = face_ele(rhs_field, sele)

              ! Calculate wall function
              call keps_wall_function(field1,field2,positions,u,bg_visc,EV,ele,sele,index,wall_fns,c_mu,rhs_field)
           end do

           ! Put values onto surface mesh
           call remap_field_to_surface(rhs_field, surface_values, surface_elements)
           ewrite_minmax(rhs_field)
           do j = 1, size(surface_node_list)
              call set(surface_field, j, node_val(surface_values, j))
           end do

           call deallocate(surface_values); call deallocate(rhs_field)

           ! Check for low reynolds boundary condition and calculate damping functions
           ! Lam-Bremhorst model (Wilcox 1998 - Turbulence modelling for CFD)
           if (wall_fns=="low_Re" .and. index==2) then
              y => extract_scalar_field(state, "DistanceToWall", stat = stat)
              if (stat /= 0) then
                 FLAbort("I need the distance to the wall - enable a DistanceToWall field")
              end if
              do node = 1, node_count(field1)
                 call keps_damping_functions(state,field2,field1,f_1,f_2,f_mu,y,bg_visc,node)
              end do
           end if

        end if
     end do boundary_conditions
  end do field_loop

end subroutine keps_bcs

!--------------------------------------------------------------------------------!
! Only used if bc type == k_epsilon for field and low_Re                         !
!--------------------------------------------------------------------------------!

subroutine keps_damping_functions(state,k,eps,f_1,f_2,f_mu,y,bg_visc,node)

  type(state_type), intent(in) :: state
  type(scalar_field), intent(inout) :: f_1, f_2, f_mu
  type(scalar_field), intent(in) :: k, eps, y
  type(tensor_field), intent(in) :: bg_visc
  integer, intent(in) :: node

  real :: rhs, Re_T, R_y, fields_max

  call get_option(trim(state%option_path)// &
       & "/subgridscale_parameterisations/k-epsilon/max_damping_value", fields_max) 

  if ((node_val(k,node) .eq. 0.0) .or. &
       & (node_val(y,node) .eq. 0.0) .or. &
       & (node_val(bg_visc,1,1,node) .eq. 0.0) .or. &
       & (node_val(eps,node) .eq. 0.0)) then
     call set(f_mu, node, 0.0)
     call set(f_1, node, 0.0)
     call set(f_2, node, 0.0)
     return
  end if

  if (node_val(bg_visc,1,1,node) /= 0.0) then
     if (node_val(eps,node) /= 0.0) then
        Re_T = node_val(k,node)**2.0 / (node_val(eps,node) * node_val(bg_visc,1,1,node))
     else 
        Re_T = 1e5
     end if
     R_y = node_val(k,node)**0.5 * node_val(y,node) / node_val(bg_visc,1,1,node)
  else
     Re_T = 1e5
     R_y = 1e5
  end if

  rhs = (- exp(- 0.0165*R_y) + 1.0)**2.0 * (20.5/Re_T + 1.0)
  if (rhs > 1.0) then
     call set(f_mu, node, 1.0)
     call set(f_1, node, 1.0)
     call set(f_2, node, 1.0)
     return
  end if
  call set(f_mu, node, min(rhs,fields_max))

  rhs = (0.05/node_val(f_mu,node))**3.0 + 1.0
  call set(f_1, node, min(rhs,fields_max))

  rhs = -exp(- Re_T**2.0) + 1.0
  call set(f_2, node, min(rhs,fields_max))

end subroutine keps_damping_functions

!--------------------------------------------------------------------------------!
! Only used if bc type == k_epsilon for field.                                   !
!--------------------------------------------------------------------------------!

subroutine keps_wall_function(field1,field2,positions,u,bg_visc,EV,ele,sele,index,wall_fns,c_mu,rhs_field)

  type(scalar_field), pointer, intent(in)              :: field1, field2, EV
  type(vector_field), pointer, intent(in)              :: positions, u
  type(tensor_field), pointer, intent(in)              :: bg_visc
  integer, intent(in)                                  :: ele, sele, index
  character(len=FIELD_NAME_LEN), intent(in)            :: wall_fns
  type(scalar_field), intent(inout)                    :: rhs_field

  type(element_type), pointer                          :: shape, fshape
  integer                                              :: i, j, gi, sgi, sloc
  real                                                 :: kappa, h, c_mu
  real, dimension(1,1)                                 :: hb
  real, dimension(ele_ngi(field1,ele))                 :: detwei
  real, dimension(face_ngi(field1,sele))               :: detwei_bdy, ustar, q_sgin, visc_sgi
  real, dimension(face_loc(field1,sele))               :: lumpedfmass
  real, dimension(ele_loc(field1,ele))                 :: sqrt_k
  real, dimension(positions%dim,1)                     :: n
  real, dimension(positions%dim,positions%dim)         :: G
  real, dimension(positions%dim,ele_ngi(field1,ele))   :: grad_k
  real, dimension(positions%dim,face_ngi(field1,sele)) :: normal_bdy, q_sgi, qq_sgin
  real, dimension(positions%dim,ele_loc(field1,ele))   :: q
  real, dimension(positions%dim,face_loc(field1,ele))  :: q_s
  real, dimension(face_loc(field1,ele),face_loc(field1,ele))   :: fmass
  real, dimension(ele_loc(field1,ele),ele_loc(field1,ele))     :: invmass
  real, dimension(positions%dim,positions%dim,ele_loc(field1,sele))       :: qq
  real, dimension(positions%dim,positions%dim,ele_ngi(field1,sele))       :: grad_u, invJ
  real, dimension(positions%dim,positions%dim,face_loc(field1,sele))      :: qq_s
  real, dimension(positions%dim,positions%dim,face_ngi(field1,sele))      :: bg_visc_sgi, qq_sgi
  real, dimension(ele_loc(field1,ele),ele_ngi(field1,ele),positions%dim)  :: dshape
  real, dimension(face_loc(rhs_field, sele))           :: rhs
  integer, dimension(face_loc(rhs_field, sele))        :: nodes_bdy

  ! Get ids, lists and shape functions
  sgi      =  face_ngi(field1, sele)    ! no. of gauss points in surface element
  sloc     =  face_loc(field1, sele)    ! no. of nodes on surface element
  shape    => ele_shape(field1, ele)    ! scalar field shape functions in volume element
  fshape   => face_shape(field1, sele)  ! scalar field shape functions in surface element
  nodes_bdy = face_global_nodes(rhs_field, sele) ! nodes in rhs_field

  ! zero rhs field
  rhs = 0.0

  ! Get shape fn gradients, element/face quadrature weights, and surface normal
  call transform_to_physical( positions, ele, shape, dshape=dshape, detwei=detwei, invJ=invJ )
  call transform_facet_to_physical( positions, sele, detwei_f=detwei_bdy, normal=normal_bdy )

  invmass = shape_shape(shape, shape, detwei)
  call invert(invmass)

  fmass = shape_shape(fshape, fshape, detwei_bdy)  !*density_gi*vfrac_gi to be safe?
  lumpedfmass = sum(fmass, 2)

  ! low Re wall functions for k and epsilon: see e.g. Wilcox (1994)
  if(wall_fns=="low_Re") then
     if (index==1) then
        rhs = 0.0
     else if (index==2) then
        bg_visc_sgi = face_val_at_quad(bg_visc,sele)
        visc_sgi = bg_visc_sgi(1,1,:)

        ! grad(k**0.5) (dim, ngi)
        sqrt_k = ele_val(field2, ele)
        sqrt_k = sqrt(abs(sqrt_k))
        ! Lifted from ele_grad_at_quad function:
        do i=1, positions%dim
           grad_k(i,:) = matmul(sqrt_k, dshape(:,:,i))
        end do

        ! grad(k**0.5) at ele_nodes (dim,loc)
        q = shape_vector_rhs(shape, grad_k, detwei)

        q = matmul(q,invmass)
        ! Pick surface nodes (dim,sloc)
        q_s = q(:,face_local_nodes(field1,sele))
        q_sgi = matmul(q_s, fshape%n)

        !dot with surface normal
        do gi = 1, sgi
           q_sgin(gi) = dot_product(q_sgi(:,gi),normal_bdy(:,gi))
        end do

        ! integral of 2*nu*(grad(k**0.5))**2 (sloc)
        rhs = shape_rhs(fshape, detwei_bdy*q_sgin**2.0*visc_sgi*2.0)
        rhs = rhs/lumpedfmass
     end if

     ! high Re shear-stress wall functions for k and epsilon: see e.g. Wilcox (1994), Mathieu p.360
  else if(wall_fns=="high_Re") then
     visc_sgi = face_val_at_quad(EV,sele)
     grad_u = ele_grad_at_quad(u, ele, dshape)

     ! grad(U) at ele_nodes (dim,dim,loc)
     qq = shape_tensor_rhs(shape, grad_u, detwei)

     do i=1,u%dim
        do j=1,u%dim
           qq(i,j,:) = matmul(invmass,qq(i,j,:))

           ! Pick surface nodes (dim,dim,sloc)
           qq_s(i,j,:) = qq(i,j,face_local_nodes(field1,sele))

           ! Get values at surface quadrature (dim,dim,sgi)
           qq_sgi(i,j,:) = matmul(qq_s(i,j,:), fshape%n)
        end do
     end do

     do gi = 1, sgi
        !dot with surface normal (dim,sgi)
        qq_sgin(:,gi) = matmul(qq_sgi(:,:,gi),normal_bdy(:,gi))

        ! Subtract normal component of velocity, leaving tangent components:
        qq_sgin(:,gi) = qq_sgin(:,gi)-normal_bdy(:,gi)*dot_product(qq_sgin(:,gi),normal_bdy(:,gi))

        ! Get streamwise component by taking sqrt(grad_n.grad_n). Multiply by eddy viscosity.
        ustar(gi) = norm2(qq_sgin(:,gi)) * visc_sgi(gi)
     end do

     if (index==1) then
        rhs = shape_rhs(fshape, detwei_bdy*ustar/c_mu**0.5)
        rhs = rhs/lumpedfmass
     else if (index==2) then
        ! calculate wall-normal element size
        G = matmul(transpose(invJ(:,:,1)), invJ(:,:,1))
        n(:,1) = normal_bdy(:,1)
        hb = 1. / sqrt( matmul(matmul(transpose(n), G), n) )
        h  = hb(1,1)
        ! Von Karman's constant
        kappa = 0.43

        rhs = shape_rhs(fshape, detwei_bdy*ustar**1.5/kappa/h)
        rhs = rhs/lumpedfmass
     end if
  else
     FLAbort("Unknown wall function option for k_epsilon boundary conditions!")
  end if

  ! Add element contribution to rhs field
  call addto(rhs_field, nodes_bdy, rhs)

end subroutine keps_wall_function

!------------------------------------------------------------------------------------!
! Called after an adapt to reset the fields and arrays within the module             !
!------------------------------------------------------------------------------------!

subroutine keps_adapt_mesh(state)

  type(state_type), intent(inout) :: state

  ewrite(1,*) "In keps_adapt_mesh"

  call keps_eddyvisc(state)
  call keps_calculate_rhs(state)
  call keps_eddyvisc(state)

end subroutine keps_adapt_mesh

!---------------------------------------------------------------------------------

subroutine k_epsilon_check_options

  ! THIS WILL ONLY WORK FOR SINGLE PHASE MODELS
  
  character(len=OPTION_PATH_LEN) :: option_path
  character(len=FIELD_NAME_LEN)  :: kmsh, emsh, vmsh
  integer                        :: dimension, stat

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
  ! Check that TurbulentKineticEnergy and TurbulentDissipation fields are on the same
  !  mesh as the velocity
  call get_option(trim(option_path)//&
       &"/scalar_field::TurbulentKineticEnergy/prognostic/mesh/name", kmsh)
  call get_option(trim(option_path)//&
       &"/scalar_field::TurbulentDissipation/prognostic/mesh/name", emsh)
  call get_option("/material_phase[0]/vector_field::Velocity/prognostic/mesh/name", vmsh,&
       & stat)
  if (stat /= 0) then
     call get_option("/material_phase[0]/vector_field::Velocity/prescribed/mesh/name", vmsh,&
       & stat)
     if (stat /= 0) then
        FLExit("You must use a prognostic or prescribed Velocity field")
     end if
  end if
  if(.not. kmsh==emsh .or. .not. kmsh==vmsh .or. .not. emsh==vmsh) then
     FLExit("You must use the Velocity mesh for TurbulentKineticEnergy and TurbulentDissipation fields")
  end if
  ! check that diffusivity is on for the two turbulent fields, and diagnostic
  if (.not.have_option(trim(option_path)//"/scalar_field::TurbulentKineticEnergy"//&
       &"/prognostic/tensor_field::Diffusivity")) then
     FLExit("You need TurbulentKineticEnergy Diffusivity field for k-epsilon")
  end if
  if (.not.have_option(trim(option_path)//&
       &"/scalar_field::TurbulentKineticEnergy/prognostic/"//&
       &"/tensor_field::Diffusivity/diagnostic/algorithm::Internal")) then
     FLExit("You need TurbulentKineticEnergy Diffusivity field set to diagnostic/internal")
  end if
  if (.not.have_option(trim(option_path)//&
       &"/scalar_field::TurbulentDissipation/prognostic"//&
       &"/tensor_field::Diffusivity")) then
     FLExit("You need TurbulentDissipation Diffusivity field for k-epsilon")
  end if
  if (.not.have_option(trim(option_path)//&
       &"/scalar_field::TurbulentDissipation/prognostic"//&
       &"/tensor_field::Diffusivity/diagnostic/algorithm::Internal")) then
     FLExit("You need TurbulentDissipation Diffusivity field set to diagnostic/internal")
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
       &"/scalar_field::TurbulentKineticEnergy/prognostic"//&
       &"/scalar_field::Absorption")) then
     FLExit("You need TurbulentKineticEnergy Absorption field for k-epsilon")
  end if
  if (.not.have_option(trim(option_path)//&
       &"/scalar_field::TurbulentKineticEnergy/prognostic"//&
       &"/scalar_field::Absorption/diagnostic/algorithm::Internal")) then
     FLExit("You need TurbulentKineticEnergy Absorption field set to diagnostic/internal")
  end if
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
  ! Velocity field options
  if (.not.have_option("/material_phase[0]/vector_field::Velocity/prognostic"//&
       "/tensor_field::Viscosity/") .and. &
       .not.have_option("/material_phase[0]/vector_field::Velocity/prescribed")) then
     FLExit("Need viscosity switched on under the Velocity field for k-epsilon.") 
     ! check that the user has switched Velocity/viscosity to diagnostic
     if (.not.have_option("/material_phase[0]/vector_field::Velocity/prognostic"//&
          &"/tensor_field::Viscosity/diagnostic/")) then
        FLExit("You need to switch the viscosity field under Velocity to diagnostic/internal")
     end if
     ! check that the user has enabled a Velocity Source field
     if (.not.have_option("/material_phase[0]/vector_field::Velocity/prognostic"//&
          &"/vector_field::Source/")) then
        FLExit("A velocity source field is required for the reynolds stress adjustment (-2/3 k delta(ij))")
     end if
  end if

end subroutine k_epsilon_check_options

!----------------------------------------------------------------------------------------! 
! Collect information required from scalar fields in order to calculate buoyancy source  !
! terms                                                                                  !
!----------------------------------------------------------------------------------------!
subroutine get_scalar_field_buoyancy_data(state, buoyant_fields, beta, delta_t)
  
  type(state_type), intent(in)                                      :: state
  type(scalar_field_pointer), allocatable, dimension(:), intent(out):: buoyant_fields
  real, allocatable, dimension(:), intent(out)                      :: beta, delta_t
  type(scalar_field_pointer), allocatable, dimension(:)             :: temp_buoyant_fields
  real, allocatable, dimension(:)                                   :: temp_beta, temp_delta_t
  type(scalar_field), pointer                                       :: field
  integer                                                           :: n_fields, i_field, stat

  ! Get number of scalar fields that are children of this state
  n_fields = scalar_field_count(state)

  ! Loop over scalar fields and copy required information to arrays if buoyancy_effects
  !  are selected for the field
  scalar_field_loop: do i_field = 1, n_fields

     field => extract_scalar_field(state, i_field)

     if (have_option(trim(field%option_path)//&
          &'/prognostic/subgridscale_parameterisation::k-epsilon/buoyancy_effects')) then

        ! determine allocation requirements, resize array and store required values
        if (allocated(buoyant_fields)) then

           allocate(temp_buoyant_fields(size(buoyant_fields, 1)))
           allocate(temp_beta(size(beta, 1)))
           allocate(temp_delta_t(size(delta_t, 1)))

           temp_buoyant_fields = buoyant_fields
           temp_beta = beta
           temp_delta_t = delta_t

           deallocate(buoyant_fields)
           deallocate(beta)
           deallocate(delta_t)

           allocate(buoyant_fields(size(temp_buoyant_fields, 1) + 1))
           allocate(beta(size(temp_beta, 1) + 1))
           allocate(delta_t(size(temp_delta_t, 1) + 1))

           buoyant_fields(1:size(temp_buoyant_fields,1)) = temp_buoyant_fields
           beta(1:size(temp_beta, 1)) = temp_beta
           delta_t(1:size(temp_delta_t, 1)) = temp_delta_t

           buoyant_fields(ubound(buoyant_fields,1))%ptr => extract_scalar_field(state, i_field)
           call get_option(trim(field%option_path)//&
                &'/prognostic/subgridscale_parameterisation::k-epsilon/buoyancy_effects/beta', &
                & beta(ubound(beta,1))) 
           call get_option(trim(field%option_path)//&
                &'/prognostic/subgridscale_parameterisation::k-epsilon/prandtl_schmidt_number', &
                & delta_t(ubound(delta_t,1)), stat) 

           deallocate(temp_buoyant_fields)
           deallocate(temp_beta)
           deallocate(temp_delta_t)

        else

           allocate(buoyant_fields(1))
           allocate(beta(1))
           allocate(delta_t(1))

           buoyant_fields(1)%ptr => extract_scalar_field(state, i_field)
           call get_option(trim(field%option_path)//&
                &'/prognostic/subgridscale_parameterisation::k-epsilon/buoyancy_effects/beta', &
                & beta(1)) 
           call get_option(trim(field%option_path)//&
                &'/prognostic/subgridscale_parameterisation::k-epsilon/prandtl_schmidt_number', &
                & delta_t(1))         

        end if

     end if

  end do scalar_field_loop

end subroutine get_scalar_field_buoyancy_data

!----------------------------------------------------------------------------------------! 
! Deallocates scalar buoyancy information                                                !
!----------------------------------------------------------------------------------------!
subroutine deallocate_scalar_field_buoyancy_data(buoyant_fields, beta, delta_t)
  
  type(scalar_field_pointer), allocatable, dimension(:), intent(inout):: buoyant_fields
  real, allocatable, dimension(:), intent(inout)                      :: beta, delta_t

  ! deallocate buoyancy data
  if (allocated(buoyant_fields)) then
     deallocate(buoyant_fields)
     deallocate(beta)
     deallocate(delta_t)
  end if

end subroutine deallocate_scalar_field_buoyancy_data

end module k_epsilon
