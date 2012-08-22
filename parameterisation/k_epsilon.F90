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
  real, save                          :: fields_min = 1.0e-10

  public :: keps_diagnostics, keps_eddyvisc, keps_bcs, keps_adapt_mesh,&
       & k_epsilon_check_options, keps_momentum_source, tensor_inner_product

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
  type(scalar_field) :: src_to_abs
  type(vector_field), pointer :: positions, u, g
  type(scalar_field), pointer :: dummydensity, density, buoyancy_density
  type(tensor_field), pointer :: diff, visc
  integer :: i, j, node, ele, term, stat
  real :: g_magnitude, c_eps_1, c_eps_2, sigma_eps, sigma_k, prandtl_schmidt_number
  logical :: prescribed, gravity = .true., have_buoyancy_turbulence = .true., lump_mass
  character(len=OPTION_PATH_LEN) :: option_path 
  character(len=FIELD_NAME_LEN), dimension(2) :: field_names
  character(len=FIELD_NAME_LEN) :: equation_type

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
  call get_option(trim(option_path)//'/prandtl_schmidt_number', prandtl_schmidt_number, default = 1.0)
  
  ! get field data
  positions    => extract_vector_field(state, "Coordinate")
  u            => extract_vector_field(state, "Velocity")
  EV           => extract_scalar_field(state, "ScalarEddyViscosity")
  visc         => extract_tensor_field(state, "BackgroundViscosity")
  f_1          => extract_scalar_field(state, "f_1")
  f_2          => extract_scalar_field(state, "f_2")
  g            => extract_vector_field(state, "GravityDirection", stat)
  call get_option('physical_parameters/gravity/magnitude', g_magnitude, stat)
  if (stat /= 0) then
     gravity = .false.
     have_buoyancy_turbulence = .false.
  end if
  
  allocate(dummydensity)
  call allocate(dummydensity, positions%mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
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
  
  if(have_buoyancy_turbulence) then
     buoyancy_density => extract_scalar_field(state, "VelocityBuoyancyDensity")
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
        call assemble_rhs_ele(src_abs_terms, fields(i)%ptr, fields(3-i)%ptr, EV, u, equation_type, &
             & density, buoyancy_density, have_buoyancy_turbulence, g, g_magnitude, gravity, positions, &
             & c_eps_1, c_eps_2, sigma_k, sigma_eps, prandtl_schmidt_number, f_1, f_2, ele, i)
     end do

     ! For non-DG we apply inverse mass globally
     if(continuity(fields(1)%ptr)>=0) then
        lump_mass = have_option(trim(option_path)//'mass_lumping_in_diagnostics/lump_mass')
        do term = 1, 3
           call solve_cg_inv_mass(state, src_abs_terms(term), lump_mass, option_path)           
        end do
     end if
     !-----------------------------------------------------------------------------------

     ! Source disabling for debugging purposes
     if(have_option(trim(option_path)//'debugging_options/disable_production')) then
        call zero(src_abs_terms(1))
     end if
     if(have_option(trim(option_path)//'debugging_options/disable_destruction')) then
        call zero(src_abs_terms(2))
     end if
     if(have_option(trim(option_path)//'debugging_options/disable_buoyancy')) then
        call zero(src_abs_terms(3))
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
  
  call deallocate(dummydensity)
  deallocate(dummydensity)

  ! Set diffusivity
  diff => extract_tensor_field(state, trim(field_names(1))//"Diffusivity")
  call zero(diff)
  do i = 1, node_count(diff)
     do j = 1, diff%dim(1)
        call addto(diff, j, j, i, node_val(visc, j, j, i))
        call addto(diff, j, j, i, node_val(EV, i) / sigma_k)
     end do
  end do
  diff => extract_tensor_field(state, trim(field_names(2))//"Diffusivity")
  call zero(diff)
  do i = 1, node_count(diff)
     do j = 1, diff%dim(1)
        call addto(diff, j, j, i, node_val(visc, j, j, i))
        call addto(diff, j, j, i, node_val(EV, i) / sigma_k)
     end do
  end do

end subroutine keps_calculate_rhs
    
!------------------------------------------------------------------------------!
! calculate the source and absorbtion terms                                    !
!------------------------------------------------------------------------------!
subroutine assemble_rhs_ele(src_abs_terms, k, eps, EV, u, equation_type, density, &
     buoyancy_density, have_buoyancy_turbulence, g, g_magnitude, gravity, &
     positions, c_eps_1, c_eps_2, sigma_k, sigma_eps, prandtl_schmidt_number, f_1, f_2, ele, field_id)

  type(scalar_field), dimension(3), intent(inout) :: src_abs_terms
  type(scalar_field), intent(in) :: k, eps, EV, f_1, f_2
  type(vector_field), intent(in) :: positions, u, g
  character(len=FIELD_NAME_LEN), intent(in) :: equation_type
  type(scalar_field), intent(in) :: density, buoyancy_density
  real, intent(in) :: g_magnitude, c_eps_1, c_eps_2, sigma_k, sigma_eps, prandtl_schmidt_number
  logical, intent(in) :: gravity, have_buoyancy_turbulence
  integer, intent(in) :: ele, field_id

  real, dimension(ele_loc(k, ele), ele_ngi(k, ele), positions%dim) :: dshape
  real, dimension(ele_ngi(k, ele)) :: detwei, strain_ngi, inv_k, rhs, EV_ele, k_ele
  real, dimension(3, ele_loc(k, ele)) :: rhs_addto
  integer, dimension(ele_loc(k, ele)) :: nodes
  real, dimension(ele_loc(k, ele), ele_loc(k, ele)) :: invmass
  real, dimension(u%dim, u%dim, ele_ngi(k, ele)) :: reynolds_stress, grad_u
  type(element_type), pointer :: shape
  integer :: i_loc, term, ngi, dim, gi, i
  
  ! For buoyancy turbulence stuff
  real, dimension(u%dim, ele_ngi(u, ele))  :: vector, u_z, u_xy
  real, dimension(ele_ngi(u, ele)) :: scalar, c_eps_3
  real, dimension(ele_loc(buoyancy_density, ele),ele_ngi(buoyancy_density, ele),positions%dim) :: dshape_density
  type(element_type), pointer :: shape_density

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

  ! Compute Reynolds stress
  grad_u = ele_grad_at_quad(u, ele, dshape)
  EV_ele = ele_val_at_quad(EV, ele)
  k_ele = ele_val_at_quad(k, ele)
  ngi = ele_ngi(u, ele)
  dim = u%dim
  do gi = 1, ngi
     reynolds_stress(:,:,gi) = EV_ele(gi)*(grad_u(:,:,gi) + transpose(grad_u(:,:,gi)))
  end do
  do i = 1, dim
     reynolds_stress(i,i,:) = reynolds_stress(i,i,:) - (2./3.)*k_ele*ele_val_at_quad(density, ele)
  end do
  ! Compute P
  rhs = tensor_inner_product(reynolds_stress, grad_u)
  if (field_id==2) then
     rhs = rhs*c_eps_1*ele_val_at_quad(f_1,ele)*ele_val_at_quad(eps,ele)*inv_k
  end if
  rhs_addto(1,:) = shape_rhs(shape, detwei*rhs)

  ! A:
  rhs = ele_val_at_quad(eps,ele)*ele_val_at_quad(density, ele)*inv_k
  if (field_id==2) then
     rhs = rhs*c_eps_2*ele_val_at_quad(f_2,ele)
  end if
  rhs_addto(2,:) = shape_rhs(shape, detwei*rhs)

  ! Gk:  
  ! Calculate buoyancy turbulence term and add to addto array
  if(have_buoyancy_turbulence) then    
  
    ! calculate scalar and vector components of the source term    
    if(.not.(buoyancy_density%mesh == k%mesh)) then
       shape_density => ele_shape(buoyancy_density, ele)
       call transform_to_physical( positions, ele, shape_density, dshape=dshape_density ) 
    else
       dshape_density = dshape
    end if
     
    scalar = -1.0*g_magnitude*ele_val_at_quad(EV, ele)/(prandtl_schmidt_number*ele_val_at_quad(density,ele))
    vector = ele_val_at_quad(g, ele)*ele_grad_at_quad(buoyancy_density, ele, dshape_density)
    
    ! multiply vector component by scalar and sum across dimensions - note that the
    ! vector part has been multiplied by the gravitational direction so that it is
    ! zero everywhere apart from in this direction.
    do gi = 1, ngi
       scalar(gi) = sum(scalar(gi) * vector(:, gi))
    end do
   
    if (field_id == 2) then
       ! get components of velocity in direction of gravity and in other directions
       u_z = abs(ele_val_at_quad(g, ele)) * ele_val_at_quad(u, ele)
       u_xy = ele_val_at_quad(u, ele) - u_z
       ! calculate c_eps_3 = tanh(v/u)
       do gi = 1, ngi
          if (norm2(u_xy(:, gi)) > fields_min) then
             c_eps_3(gi) = tanh(norm2(u_z(:, gi))/norm2(u_xy(:, gi))) 
          else
             c_eps_3(gi) = 1.0
          end if
       end do     
       scalar = scalar*c_eps_1*ele_val_at_quad(f_1,ele)*c_eps_3*ele_val_at_quad(eps,ele)*inv_k
    end if

    ! multiply by determinate weights, integrate and assign to rhs
    rhs_addto(3,:) = shape_rhs(shape, scalar * detwei)
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
! eddyvisc calculates the lengthscale, and then the eddy viscosity!
! Eddy viscosity is added to the background viscosity.
!----------
subroutine keps_eddyvisc(state)

  type(state_type), intent(inout)  :: state
  type(tensor_field), pointer      :: eddy_visc, viscosity, bg_visc
  type(vector_field), pointer      :: positions, u
  type(scalar_field), pointer      :: kk, eps, EV, ll, f_mu, density, dummydensity
  type(scalar_field)               :: ev_rhs
  integer                          :: i, j, ele, stat
  
  ! Options grabbed from the options tree
  real                             :: c_mu, lmax
  character(len=OPTION_PATH_LEN)   :: option_path
  logical                          :: lump_mass, limit_length_scale, have_visc = .true.
  character(len=FIELD_NAME_LEN)    :: equation_type

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-epsilon/'

  if (.not. have_option(trim(option_path))) then 
     return
  end if

  ewrite(1,*) "In keps_eddyvisc"

  ! Get model constants
  limit_length_scale = have_option(trim(option_path)//'/limit_length_scale')
  if(limit_length_scale) then
     ! Limit length scale to prevent instablilities.
     call get_option(trim(option_path)//'/limit_length_scale', lmax)
  end if
  call get_option(trim(option_path)//'/C_mu', c_mu, default = 0.09)
  
  ! Get field data
  kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
  eps        => extract_scalar_field(state, "TurbulentDissipation")
  positions  => extract_vector_field(state, "Coordinate")
  u          => extract_vector_field(state, "Velocity")
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
  
  allocate(dummydensity)
  call allocate(dummydensity, positions%mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
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

  call allocate(ev_rhs, EV%mesh, name="EVRHS")
  call zero(ev_rhs)

  ! Initialise viscosity to background value
  if (have_visc) then
     call set(viscosity, bg_visc)
  end if
  
  ! Compute the length scale diagnostic field here.
  do i = 1, node_count(EV)
     ! Also clip the k and epsilon fields at the nodes.
     call set(kk, i, max(node_val(kk,i), fields_min))
     call set(eps, i, max(node_val(eps,i), fields_min))
     
     ! Limit lengthscale to prevent instablilities.
     if(limit_length_scale) then
        call set(ll, i, min(node_val(kk,i)**1.5 / node_val(eps,i), lmax))
     else
        call set(ll, i, node_val(kk,i)**1.5 / node_val(eps,i))
     end if
  end do

  ! Calculate scalar eddy viscosity by integration over element
  do ele = 1, ele_count(EV)
     call keps_eddyvisc_ele(ele, positions, kk, eps, EV, ll, f_mu, density, ev_rhs)
  end do

  ! For non-DG we apply inverse mass globally
  if(continuity(EV)>=0) then
     lump_mass = have_option(trim(option_path)//'mass_lumping_in_diagnostics/lump_mass')
     call solve_cg_inv_mass(state, ev_rhs, lump_mass, option_path)  
  end if
  
  ! Allow for prescribed eddy-viscosity
  if (.not. have_option(trim(option_path)//'/scalar_field::ScalarEddyViscosity/prescribed')) then
     call set(EV, ev_rhs)
  end if

  call deallocate(ev_rhs)
  
  call deallocate(dummydensity)
  deallocate(dummydensity)

  ewrite(2,*) "Setting k-epsilon eddy-viscosity tensor"
  call zero(eddy_visc)

  ! this is skipped if zero_eddy_viscosity is set - this is the easiest way to
  ! disable feedback from the k-epsilon model back into the rest of the model
  if (.not. have_option(trim(option_path)//'debugging_options/zero_reynolds_stress_tensor')) then
     do i = 1, eddy_visc%dim(1)
        do j = 1, eddy_visc%dim(1)
           call set(eddy_visc, i, j, EV)
        end do
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
  
  
  contains
  
   subroutine keps_eddyvisc_ele(ele, positions, kk, eps, EV, ll, f_mu, density, ev_rhs)
   
      type(vector_field), pointer      :: positions
      type(scalar_field), pointer      :: kk, eps, EV, ll, f_mu, density
      type(scalar_field), intent(inout) :: ev_rhs
      integer, intent(in)              :: ele
      
      
      type(element_type), pointer      :: shape_ev
      integer, pointer, dimension(:)   :: nodes_ev
      integer                          :: i 
      real, dimension(ele_ngi(EV, ele)) :: detwei
      real, dimension(ele_loc(EV, ele)) :: rhs_addto
      real, dimension(ele_loc(EV, ele), ele_loc(EV, ele)) :: invmass
      real, dimension(ele_ngi(kk, ele)) :: kk_at_quad, eps_at_quad, ll_at_quad
      
   
      nodes_ev => ele_nodes(EV, ele)
      shape_ev =>  ele_shape(EV, ele)
      
      ! Get detwei
      call transform_to_physical(positions, ele, detwei=detwei)
      
      ! Get the k, epsilon and the length scale values at the Gauss points
      kk_at_quad = ele_val_at_quad(kk,ele)
      eps_at_quad = ele_val_at_quad(eps,ele)
      if(limit_length_scale) then
         ll_at_quad = ele_val_at_quad(ll,ele)
      end if
      
      ! Clip the field values at the Gauss points.
      ! Note 1: This isn't a permanent change directly to the field itself,
      ! only to the values used in the computation of the eddy viscosity.
      ! Note 2: Can't allow negative/zero epsilon or k.
      ! Note 3: Here we assume all fields have the same number of
      ! Gauss points per element.
      do i = 1, ele_ngi(kk, ele)
         ! k
         if(kk_at_quad(i) < fields_min) then
            kk_at_quad(i) = fields_min
         end if
         ! epsilon
         if(eps_at_quad(i) < fields_min) then
            eps_at_quad(i) = fields_min
         end if
         ! Limit lengthscale to prevent instablilities.
         if(limit_length_scale) then
            if(ll_at_quad(i) > lmax) then
               ll_at_quad(i) = lmax
            end if 
         end if
      end do
      
      ! Compute the eddy viscosity
      if(limit_length_scale) then
         rhs_addto = shape_rhs(shape_ev, detwei*C_mu*ele_val_at_quad(f_mu,ele)* &
                     (kk_at_quad**0.5)*ll_at_quad)
      else
         rhs_addto = shape_rhs(shape_ev, detwei*C_mu*ele_val_at_quad(density,ele)*&
                     ele_val_at_quad(f_mu,ele)*(kk_at_quad**2.0)/eps_at_quad)
      end if
            
      ! In the DG case we will apply the inverse mass locally.
      if(continuity(EV)<0) then
         invmass = inverse(shape_shape(shape_ev, shape_ev, detwei))
         rhs_addto = matmul(rhs_addto, invmass)
      end if
      
      ! Add the element's contribution to the nodes of ev_rhs
      call addto(ev_rhs, nodes_ev, rhs_addto)    
   
   end subroutine keps_eddyvisc_ele

end subroutine keps_eddyvisc

!--------------------------------------------------------------------------------!
! calculates the reynolds stress tensor correction term 
! - grad(2/3 k rho delta(ij)) 
! this is added to prescribed momentum source fields or sets diagnostic source 
! fields 
!--------------------------------------------------------------------------------!
subroutine keps_momentum_source(state)

  type(state_type), intent(inout)  :: state
  
  type(scalar_field), pointer :: k, lumped_mass, dummydensity, density
  type(vector_field), pointer :: source, x, u
  type(vector_field) :: rhs
  type(vector_field) :: prescribed_source
  integer :: i
  logical :: prescribed, lump_mass
  character(len=OPTION_PATH_LEN) :: option_path 
  character(len=FIELD_NAME_LEN)    :: equation_type

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-epsilon/'

  if (.not.have_option(trim(option_path)) .or. &
      have_option(trim(option_path)//'debugging_options/zero_reynolds_stress_tensor') .or. &
      have_option(trim(state%option_path)//"/vector_field::Velocity/prescribed")) then 
     return
  end if

  k         => extract_scalar_field(state, "TurbulentKineticEnergy")
  source    => extract_vector_field(state, "VelocitySource")
  x         => extract_vector_field(state, "Coordinate")
  u         => extract_vector_field(state, "Velocity")
  
  allocate(dummydensity)
  call allocate(dummydensity, x%mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
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
  
  ! Allow for prescribed momentum source
  prescribed = (have_option(trim(source%option_path)//'/prescribed/'))
  if(.not.prescribed) then
     call zero(source)
  end if

  call allocate(rhs, source%dim, source%mesh, name='TempSource')
  call zero(rhs)
  do i = 1, ele_count(k)
     call keps_momentum_source_ele()
  end do

  ! For non-DG we apply inverse mass globally
  if(continuity(k)>=0) then
     lump_mass = have_option(trim(option_path)//&
          'mass_lumping_in_diagnostics/lump_mass')
     call solve_cg_inv_mass_vector(state, rhs, lump_mass, option_path)  
     call addto(source, rhs)
  end if
  
  call deallocate(rhs)
  
  call deallocate(dummydensity)
  deallocate(dummydensity)

  ! This code isn't needed because the adjustment to the field is done before the non
  ! -linear iteration loop
  ! ! Allow for prescribed momentum source
  ! prescribed = (have_option(trim(source%option_path)//'/prescribed/'))
  ! if(prescribed) then
  !    call allocate(prescribed_source, source%dim, source%mesh, name='PrescribedSource')
  !    call initialise_field_over_regions(prescribed_source, &
  !         trim(source%option_path)//'/prescribed/value', &
  !         x)
  !    call addto(source, prescribed_source)
  !    call deallocate(prescribed_source)
  ! end if

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
  type(scalar_field), pointer                :: density, dummydensity
  type(vector_field), pointer                :: positions, u
  type(tensor_field), pointer                :: bg_visc
  type(scalar_field)                         :: rhs_field, surface_values
  type(mesh_type), pointer                   :: surface_mesh
  integer                                    :: i, j, ele, sele, index, nbcs, stat, node
  integer, dimension(:), pointer             :: surface_elements, surface_node_list
  character(len=FIELD_NAME_LEN)              :: bc_type, bc_name, wall_fns
  character(len=OPTION_PATH_LEN)             :: bc_path, bc_path_i, option_path 
  real                                       :: c_mu
  character(len=FIELD_NAME_LEN)              :: equation_type

  option_path = trim(state%option_path)//'/subgridscale_parameterisations/k-epsilon/'

  ewrite(2,*) "In keps_bcs"

  positions => extract_vector_field(state, "Coordinate")
  u         => extract_vector_field(state, "Velocity")
  EV        => extract_scalar_field(state, "ScalarEddyViscosity")
  bg_visc   => extract_tensor_field(state, "BackgroundViscosity")
  f_1       => extract_scalar_field(state, "f_1")
  f_2       => extract_scalar_field(state, "f_2")
  f_mu      => extract_scalar_field(state, "f_mu")

  allocate(dummydensity)
  call allocate(dummydensity, positions%mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
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

           if(wall_fns=="high_Re") then
              call allocate(surface_values, surface_mesh, name="surfacevalues")
              call allocate(rhs_field, field1%mesh, name="rhs")
              call zero(surface_values); call zero(rhs_field)

              do j = 1, ele_count(surface_mesh)
                 sele = surface_elements(j)
                 ele  = face_ele(rhs_field, sele)
                 
                 ! Calculate wall function
                 call keps_wall_function(field1,field2,positions,u,bg_visc,EV,density,ele,sele,index,c_mu,rhs_field)
              end do
              
              ! Put values onto surface mesh
              call remap_field_to_surface(rhs_field, surface_values, surface_elements)
              ewrite_minmax(rhs_field)
              do j = 1, size(surface_node_list)
                 call set(surface_field, j, node_val(surface_values, j))
              end do

              call deallocate(surface_values); call deallocate(rhs_field)
           end if           

           ! assume low Re wall functions if .not. high_Re: see e.g. Wilcox (1994)
           ! k = 0; dEps/dy = 0
           ! Now check for low reynolds boundary condition and calculate damping functions
           ! Lam-Bremhorst model (Wilcox 1998 - Turbulence modelling for CFD)
           if (wall_fns=="low_Re" .and. index==2) then
              y => extract_scalar_field(state, "DistanceToWall", stat = stat)
              if (stat /= 0) then
                 FLAbort("I need the distance to the wall - enable a DistanceToWall field")
              end if
              do node = 1, node_count(field1)
                 call keps_damping_functions(state,field2,field1,f_1,f_2,f_mu,y,bg_visc,density,node)
              end do
           end if

        end if
     end do boundary_conditions
  end do field_loop

  ! This is for the lowRe MMS test. This was complicated as values are limited actually
  ! near the wall so we need a way to test this code away from directly next to the
  ! wall. This was the only way I could see of enabling the damping functions without
  ! actually having to have a low_Re boundary. See tests/mms_rans_p2p1_keps_lowRe
  if (have_option(trim(option_path)//"/debugging_options/enable_lowRe_damping")) then
     field1 => extract_scalar_field(state, "TurbulentDissipation")
     field2 => extract_scalar_field(state, "TurbulentKineticEnergy")
     y => extract_scalar_field(state, "DistanceToWall", stat = stat)
     if (stat /= 0) then
        FLAbort("I need the distance to the wall - enable a DistanceToWall field")
     end if
     do node = 1, node_count(field1)
        call keps_damping_functions(state,field2,field1,f_1,f_2,f_mu,y,bg_visc,density,node)
     end do
  end if

  call deallocate(dummydensity)
  deallocate(dummydensity)

end subroutine keps_bcs

!--------------------------------------------------------------------------------!
! Only used if bc type == k_epsilon for field and low_Re                         !
!--------------------------------------------------------------------------------!

subroutine keps_damping_functions(state,k,eps,f_1,f_2,f_mu,y,bg_visc,density,node)

  type(state_type), intent(in) :: state
  type(scalar_field), intent(inout) :: f_1, f_2, f_mu
  type(scalar_field), intent(in) :: k, eps, y, density
  type(tensor_field), intent(in) :: bg_visc
  integer, intent(in) :: node

  real :: rhs, Re_T, R_y, fields_max

  call get_option(trim(state%option_path)// &
       & "/subgridscale_parameterisations/k-epsilon/max_damping_value", fields_max) 

  if ((node_val(k,node) .eq. 0.0) .or. &
       & (node_val(y,node) .eq. 0.0) .or. &
       & (node_val(bg_visc,1,1,node) .eq. 0.0) .or. &
       & (node_val(density,node) .eq. 0.0) .or. &
       & (node_val(eps,node) .eq. 0.0)) then
     call set(f_mu, node, 0.0)
     call set(f_1, node, 0.0)
     call set(f_2, node, 0.0)
     return
  end if

  if (node_val(bg_visc,1,1,node) /= 0.0) then
     if (node_val(eps,node) /= 0.0) then
        Re_T = (node_val(density,node) * node_val(k,node)**2.0) / (node_val(eps,node) * node_val(bg_visc,1,1,node))
     else 
        Re_T = 1e5
     end if
     R_y = (node_val(density,node) * node_val(k,node)**0.5 * node_val(y,node)) / node_val(bg_visc,1,1,node)
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

subroutine keps_wall_function(field1,field2,positions,u,bg_visc,EV,density,ele,sele,index,c_mu,rhs_field)

  type(scalar_field), pointer, intent(in)              :: field1, field2, EV, density
  type(vector_field), pointer, intent(in)              :: positions, u
  type(tensor_field), pointer, intent(in)              :: bg_visc
  integer, intent(in)                                  :: ele, sele, index
  type(scalar_field), intent(inout)                    :: rhs_field

  type(element_type), pointer                          :: shape, fshape
  integer                                              :: i, j, gi, sgi, sloc
  real                                                 :: kappa, h, c_mu
  real, dimension(1,1)                                 :: hb
  real, dimension(ele_ngi(field1,ele))                 :: detwei
  real, dimension(face_ngi(field1,sele))               :: detwei_bdy, ustar, q_sgin, visc_sgi, density_sgi
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

  ! high Re shear-stress wall functions for k and epsilon: see e.g. Wilcox (1994), Mathieu p.360
  visc_sgi = face_val_at_quad(EV,sele)
  density_sgi = face_val_at_quad(density,sele)
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
     ! Note: we divide by density here to write the BC in terms of a dynamic viscosity,
     ! rather than a kinematic one.
     ustar(gi) = norm2(qq_sgin(:,gi)) * visc_sgi(gi) / density_sgi(gi)
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

subroutine k_epsilon_check_options(state)

  type(state_type) :: state  
  character(len=OPTION_PATH_LEN) :: option_path
  character(len=FIELD_NAME_LEN)  :: kmsh, emsh, vmsh
  integer                        :: dimension, stat

  ewrite(1,*) "In keps_check_options"
  option_path = trim(state%option_path)//"/subgridscale_parameterisations/k-epsilon"

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
  call get_option(trim(state%option_path)//"/vector_field::Velocity/prognostic/mesh/name", vmsh,&
       & stat)
  if (stat /= 0) then
     call get_option(trim(state%option_path)//"/vector_field::Velocity/prescribed/mesh/name", vmsh,&
       & stat)
     if (stat /= 0) then
        FLExit("You must use a prognostic or prescribed Velocity field")
     end if
  end if
  if(.not. kmsh==emsh .or. .not. kmsh==vmsh .or. .not. emsh==vmsh) then
     FLExit("You must use the Velocity mesh for TurbulentKineticEnergy and TurbulentDissipation fields")
  end if

  ! Velocity field options
  if (.not.have_option(trim(state%option_path)//"/vector_field::Velocity/prognostic"//&
       "/tensor_field::Viscosity/") .and. &
       .not.have_option(trim(state%option_path)//"/vector_field::Velocity/prescribed")) then
     FLExit("Need viscosity switched on under the Velocity field for k-epsilon.") 
  end if
  ! check that the user has switched Velocity/viscosity to diagnostic
  if (have_option(trim(state%option_path)//"/vector_field::Velocity/prognostic") .and. &
       .not.have_option(trim(state%option_path)//"/vector_field::Velocity/prognostic"//&
       "/tensor_field::Viscosity/diagnostic/")) then
     FLExit("You need to switch the viscosity field under Velocity to diagnostic/internal")
  end if
  ! check that the user has enabled a Velocity Source field
  if (have_option(trim(state%option_path)//"/vector_field::Velocity/prognostic") .and. &
       .not.have_option(trim(state%option_path)//"/vector_field::Velocity/prognostic"//&
       &"/vector_field::Source/")) then
     FLExit("A velocity source field is required for the reynolds stress adjustment (-2/3 k delta(ij))")
  end if

  ! Check ScalarEddyViscosity is diagnostic
  if (have_option(trim(option_path)//'/scalar_field::ScalarEddyViscosity/prescribed')) then
     ewrite(0,*) "WARNING: ScalarEddyViscosity field is prescribed"
  end if

end subroutine k_epsilon_check_options

end module k_epsilon
