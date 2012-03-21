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
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module copy_outof_into_state
  !! This module enables the multiphase prototype code to interact with state
  !! First copying everything required from state to the old variable space
  !! Second copying what is needed back to state for output 

  use fldebug
  use state_module
  use fields
  use fields_base
  use field_options
  use spud
  use populate_state_module
  use diagnostic_variables
  use diagnostic_fields
  use diagnostic_fields_wrapper
  use global_parameters, only: option_path_len
  use diagnostic_fields_wrapper_new

  use element_numbering
  use boundary_conditions
  
  use quicksort


  implicit none
  
  private
  
  public :: copy_outof_state, copy_into_state
  
  contains
  
    subroutine copy_outof_state(state, &
         nonlinear_iterations, nonlinear_iteration_tolerance, &
                                ! Begin here all the variables from read_scalar
         nphases, ncomps, totele, ndim, nlev, &
         u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, mat_nloc, &
         cv_snloc, u_snloc, p_snloc, stotel, &
         cv_ndgln, u_ndgln, p_ndgln, x_ndgln, xu_ndgln, mat_ndgln, &
         cv_sndgln, p_sndgln, u_sndgln, &
         ncoef, nuabs_coefs, &
         u_ele_type, p_ele_type, mat_ele_type, cv_ele_type, &
         cv_sele_type, u_sele_type, &
         nits, nits_internal, ndpset, noit_dim, &
         nits_flux_lim_volfra, nits_flux_lim_comp, &
         t_disopt, u_disopt, v_disopt, t_dg_vel_int_opt, &
         u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
         capil_pres_opt, ncapil_pres_coef, comp_diffusion_opt, ncomp_diff_coef, &
         patmos, p_ini, t_beta, v_beta, t_theta, v_theta, &
         u_theta, domain_length, &
         lump_eqns, volfra_use_theta_flux, volfra_get_theta_flux, &
         comp_use_theta_flux, comp_get_theta_flux, &
         t_use_theta_flux, t_get_theta_flux, &
                                ! Now the variables from read_all
         in_ele_upwind, dg_ele_upwind, &
         Mobility, alpha_beta, &
         KComp_Sigmoid, Comp_Sum2One, &
         wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc, &
         wic_comp_bc, uabs_option, eos_option, cp_option, &
         suf_vol_bc, suf_d_bc, suf_cpd_bc, suf_t_bc, suf_p_bc, &
         suf_u_bc, suf_v_bc, suf_w_bc, suf_one_bc, suf_comp_bc, &
         suf_u_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob1, suf_v_bc_rob2, &
         suf_w_bc_rob1, suf_w_bc_rob2, suf_t_bc_rob1, suf_t_bc_rob2, &
         suf_vol_bc_rob1, suf_vol_bc_rob2, &
         suf_comp_bc_rob1, suf_comp_bc_rob2, &
         x, y, z, xu, yu, zu, nu, nv, nw, ug, vg, wg, &
         u_source, t_source, v_source, comp_source, &
         u, v, w, &
         den, satura, volfra, comp, t, p, cv_p, volfra_pore, &
         cv_one, Viscosity, & 
                                ! cv_one corresponds to density in single-phase advection problem
         uabs_coefs, &
         eos_coefs, cp_coefs, &
         comp_diff_coef, capil_pres_coef, &
         u_abs_stab, u_absorb, comp_absorb, &
         t_absorb, v_absorb, &
         perm, K_Comp, &
         comp_diffusion, &
                                ! Now adding other things which we have taken inside this routine to define
         cv_nonods, p_nonods, u_nonods, x_nonods, xu_nonods, mat_nonods, &
         have_temperature_fields, nits_flux_lim_t)

      type(state_type), dimension(:), intent(inout) :: state

      logical, intent(inout) :: have_temperature_fields

      ! coordinate, velocity and pressure meshes
      type(mesh_type), pointer :: cmesh => null()
      type(mesh_type), pointer :: vmesh => null()
      type(mesh_type), pointer :: pmesh  => null()

      type(scalar_field), pointer :: porosity, permeability, density, pressure, &
           phasevolumefraction, pvf_source, &
           scalarfield, scalarfield_source, &
           componentmassfraction, &
           phasevolumefraction_bc, density_bc, pressure_bc, scalarfield_bc

      type(vector_field), pointer :: velocity, velocity_source
      type(vector_field), pointer :: velocity_bc

      type(tensor_field), pointer :: viscosity_ph1, viscosity_ph2, t_permeability

      integer :: nonlinear_iterations, stat, nstates

      real :: nonlinear_iteration_tolerance

      character(len=OPTION_PATH_LEN) :: field_name, material_phase_name, option_path

      !! temporary variables only needed for interfacing purposes

      type(vector_field), pointer :: positions
      type(vector_field) :: cv_positions

      integer :: i, j, k, l, cv_nonods, p_nonods, mat_nonods, &
           x_nonods, xu_nonods, u_nonods, cv_nod, u_nod, xu_nod, x_nod, mat_nod, &
           shared_nodes, u_nloc2

      real :: coord_min, coord_max, eos_value

      real :: const, dx
      real, dimension(:), allocatable :: const_vec
      real, dimension(:,:), allocatable :: const_array


      integer, dimension(:), allocatable :: cv_ndgln, u_ndgln, p_ndgln, x_ndgln, xu_ndgln, mat_ndgln, &
           cv_sndgln, p_sndgln, u_sndgln
      integer, dimension(:), pointer :: element_nodes

      real, dimension(:), allocatable :: initial_constant_velocity

      logical :: is_isotropic, is_symmetric, is_diagonal, ndgln_switch

      !! Variables needed by the prototype code
      !! and therefore needing to be pulled out of state or
      !! derived from information in state:

      ! Scalars (from read_scalar())
      integer :: nphases, ncomps, totele, ndim, nlev, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, mat_nloc, &
           cv_snloc, u_snloc, p_snloc, stotel, &
           ncoef, nuabs_coefs, &
           u_ele_type, p_ele_type, mat_ele_type, cv_ele_type, &
           cv_sele_type, u_sele_type

      integer :: nits, nits_internal, ndpset, noit_dim, &
           nits_flux_lim_volfra, nits_flux_lim_comp, nits_flux_lim_t

      integer :: t_disopt, u_disopt, v_disopt, t_dg_vel_int_opt, &
           u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt

      integer :: capil_pres_opt, ncapil_pres_coef, comp_diffusion_opt, ncomp_diff_coef

      real :: patmos, p_ini, t_beta, v_beta, t_theta, v_theta, &
           u_theta, domain_length

      logical :: lump_eqns, volfra_use_theta_flux, volfra_get_theta_flux, &
           comp_use_theta_flux, comp_get_theta_flux, t_use_theta_flux, t_get_theta_flux


      ! Others (from read_all())

      integer :: in_ele_upwind, dg_ele_upwind

      real :: Mobility, alpha_beta

      logical :: KComp_Sigmoid, Comp_Sum2One

      integer, dimension( : ), allocatable :: wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc, &
           wic_comp_bc, uabs_option, eos_option, cp_option

      real, dimension( : ), allocatable :: suf_vol_bc, suf_d_bc, suf_cpd_bc, suf_t_bc, suf_p_bc, &
           suf_u_bc, suf_v_bc, suf_w_bc, suf_one_bc, suf_comp_bc, &
           suf_u_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob1, suf_v_bc_rob2, &
           suf_w_bc_rob1, suf_w_bc_rob2, suf_t_bc_rob1, suf_t_bc_rob2, &
           suf_vol_bc_rob1, suf_vol_bc_rob2, &
           suf_comp_bc_rob1, suf_comp_bc_rob2, &
           x, y, z, xu, yu, zu, nu, nv, nw, ug, vg, wg, &
           u_source, t_source, v_source, comp_source, &
           u, v, w, &
           den, satura, volfra, comp, t, p, cv_p, volfra_pore, &
           cv_one, Viscosity

      real, dimension( :, : ), allocatable :: uabs_coefs
      real, dimension( :, : ), allocatable :: eos_coefs, cp_coefs

      real, dimension( : , : , : ), allocatable :: comp_diff_coef, capil_pres_coef, &
           u_abs_stab, u_absorb, comp_absorb, &
           t_absorb, v_absorb, &
           perm, K_Comp

      real, dimension( : , : , : , : ), allocatable :: comp_diffusion

      integer :: Velocity_BC_Type, Pressure_BC_Type, density_bc_type, &
           component_bc_type, pvf_bc_type, Temperature_bc_type, shape_option(2), kk

      real :: component_suf_bc

      integer, dimension(:), allocatable :: Velocity_SufID_BC, Pressure_SufID_BC, density_sufid_bc, &
           component_sufid_bc, pvf_sufid_bc, Temperature_sufid_bc

      ! Gravity terms to be linked with u_source
      logical :: have_gravity
      real :: gravity_magnitude, delta_den, grm
      type(vector_field) :: gravity_direction

      integer :: nobcs, nodi

      character(len=option_path_len) :: vel_element_type
      integer, dimension( : ), allocatable :: index
      real, dimension( : ), allocatable :: x_temp
      logical :: is_overlapping   

      real :: val

      ewrite(3,*) 'In copy_outof_state'

      nlev = 3

      nstates = option_count("/material_phase")

      ewrite(3,*) ' nstates:', nstates

      ! Assume there are the same number of components in each phase (will need to check this eventually)
      ncomps=0
      do i=1,nstates
         if (have_option("/material_phase[" // int2str(i-1) // "]/is_multiphase_component")) then
            ncomps=ncomps+1
         end if
      end do

      nphases = nstates - ncomps
      assert(nphases > 0)
      ewrite(3,*)'nstates, ncomps, nphases:', nstates, ncomps, nphases

      ! Get the coordinate, velocity and pressure mesh (which must exist)
      cmesh => extract_mesh(state, "CoordinateMesh", stat=stat)

      if (stat /= 0) then
         FLExit('Require mesh with name CoordinateMesh')
      end if

      vmesh => extract_mesh(state, "VelocityMesh", stat=stat)

      if (stat /= 0) then
         FLExit('Require mesh with name VelocityMesh')
      end if

      pmesh => extract_mesh(state, "PressureMesh", stat=stat)

      if (stat /= 0) then
         FLExit('Require mesh with name PressureMesh')
      end if

      positions => extract_vector_field( state, "Coordinate" )
      cv_positions = get_coordinate_field( state(1), pmesh )
      totele = ele_count( cmesh )
      stotel = surface_element_count( cmesh )


      call get_option("/geometry/dimension",ndim)

      if( ndim == 1 ) then
         ! This is a strictly 1d property, so will have a suitably 1d method for finding it in state!
         coord_min=1.0e9
         coord_max=-1.0e-9
         do i=1,node_count(positions)
            coord_min=min(coord_min,positions%val(X_,i))
            coord_max=max(coord_max,positions%val(X_,i))
         end do
         domain_length = coord_max - coord_min
         ! Still need dx for the old way of setting up x (and y, z)
         dx = domain_length/real(totele)
      end if

      ! get the vel element type.
      call get_option('/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           vel_element_type)

      select case(ndim)
      case(1)
         xu_nloc = 2
         cv_nloc = pmesh%shape%degree + 1
         shared_nodes = 1
      case(2)
         xu_nloc = 3
         cv_nloc = tr(pmesh%shape%degree+1)
         select case(pmesh%shape%degree)
         case(1)
            shared_nodes=2
         case(2)
            shared_nodes=3
         end select
      case(3)
         xu_nloc = 4
         cv_nloc = te(pmesh%shape%degree+1)
         select case(pmesh%shape%degree)
         case(1)
            shared_nodes=3
         case(2)
            shared_nodes=6
         end select
      end select
      x_nloc = cv_nloc
      p_nloc = cv_nloc
      mat_nloc = cv_nloc
      u_nloc = xu_nloc    

      is_overlapping = .false.
      if ( trim( vel_element_type ) == 'overlapping' ) is_overlapping = .true. 

      if( is_overlapping ) u_nloc = cv_nloc * xu_nloc

      ewrite(3,*)'ndim, cv_nloc, x_nloc, xu_nloc, u_nloc:', ndim, cv_nloc, x_nloc, xu_nloc, u_nloc 

      if (pmesh%continuity>=0) then
         ! Continuous pressure mesh
         cv_nonods = ( cv_nloc - shared_nodes ) * totele + shared_nodes
      else
         ! Discontinuous pressure mesh
         cv_nonods = cv_nloc * totele
      endif

      ! u_snloc is 1 for the single phase advection case:
      if (nstates==1) then
         u_snloc = 1
         xu_nloc = u_nloc
         p_nonods = u_nloc * totele ! = u_nonods
      else
         u_snloc = 3
         p_nonods = cv_nonods
      endif

      Select Case( ndim )
      case( 1 )
         cv_snloc = 1
      case( 2 )
         u_snloc = face_loc( vmesh, 1 )
         cv_snloc = face_loc( pmesh, 1 )
         ! cv_snloc = 3 ! quadratic triangle
      case default; FLAbort("Incorrect number of dimension - ndim")
      end Select
      p_snloc = 2
      ! x_snloc = 1
      mat_nonods = mat_nloc * totele
      if( is_overlapping ) u_snloc = u_snloc * cv_nloc

      ewrite(3,*)'cv_nonods, u_snloc, cv_snloc:', cv_nonods, u_snloc, cv_snloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !! GLOBAL NUMBERING SECTION
      !!
      !! nb. This could all be moved to a separate subrtn now, but it needed to be moved
      !! as it's needed by the field setup bits below
      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!
      !!  x, y, z are ordered in the global numbering order, and not in spatial order
      !!  This seems to be consistent with other places in the code (cv-adv-dif, etc)
!!!!!!
      !! Switch on (true) or off (false) using Fluidity node numbering
      ndgln_switch = .true.
      if( ndim == 1 ) ndgln_switch = .false.

      xu_nonods = max( ( xu_nloc - 1 ) * totele + 1, totele )

      allocate( cv_ndgln( totele * cv_nloc ) )
      allocate( p_ndgln( totele * p_nloc ) )
      allocate( u_ndgln( totele * u_nloc ) )
      allocate( xu_ndgln( totele * xu_nloc ) )
      allocate( x_ndgln( totele * cv_nloc ) )
      allocate( mat_ndgln( totele * mat_nloc ) )

      allocate( xu( xu_nonods ) )
      allocate( yu( xu_nonods ) )
      allocate( zu( xu_nonods ) )
      xu = 0.
      yu = 0.
      zu = 0.

      u_nonods = u_nloc * totele

      allocate( nu( u_nonods * nphases ))
      allocate( nv( u_nonods * nphases ))
      allocate( nw( u_nonods * nphases ))
      allocate( ug( u_nonods * nphases ))
      allocate( vg( u_nonods * nphases ))
      allocate( wg( u_nonods * nphases ))
      nu = 1.
      nv = 0.
      nw = 0.
      if ( ndim > 1 ) nv = 1.
      if ( ndim > 2 ) nw = 1.
      ug = 0.
      vg = 0.
      wg = 0.


      !! global numbering to allow proper copying of fields
      cv_nod = 0
      u_nod = 0
      xu_nod = 0
      x_nod = 0
      mat_nod = 0
      pressure => extract_scalar_field( state(1), "Pressure")
      velocity => extract_vector_field( state(1), "Velocity")
      Loop_OverElements1: do k = 1,totele
         if ( ndgln_switch ) then
            element_nodes => ele_nodes( pressure, k )
            do j = 1, size( element_nodes )
               cv_ndgln( ( k - 1 ) * size( element_nodes ) + j ) = element_nodes( j )  
               !ewrite(3,*)'ele, nod, cv_ndgln:', &
               !     k, j, cv_ndgln( ( k - 1 ) * size( element_nodes ) + j )          
            end do
         else
            do j = 1, cv_nloc
               cv_nod = cv_nod + 1
               cv_ndgln( ( k - 1 ) * cv_nloc + j ) = cv_nod
            end do
            if( cv_nonods /= totele * cv_nloc ) cv_nod = cv_nod - 1
         end if

         !! mat_ndgln on the pressure mesh but always discontinuous
         if ( ndgln_switch ) then
            element_nodes => ele_nodes( pressure, k )
            do j = 1, size( element_nodes )
               mat_ndgln( ( k - 1 ) * size( element_nodes ) + j ) = element_nodes( j )  
            end do
         else
            do j = 1, mat_nloc
               mat_nod = mat_nod + 1
               mat_ndgln( ( k - 1 ) * mat_nloc + j ) = mat_nod
               !ewrite(3,*)'ele, nod, mat_ndgln:', k, j, mat_nod
            end do
         end if

         !if( ndgln_switch .and. ndim>1 ) then
         !   element_nodes => ele_nodes( velocity, k )
         !   do j = 1, size( element_nodes )
         !      u_ndgln( ( k - 1 ) * size( element_nodes ) + j ) = element_nodes( j )  
         !      ewrite(3,*)'ele, nod, u_ndgln:', k, j, u_ndgln( ( k - 1 ) * size( element_nodes ) + j )
         !   end do
         !else
         do j = 1, u_nloc
            u_nod = u_nod + 1
            u_ndgln(( k - 1 ) * u_nloc + j ) = u_nod
            ewrite(3,*)'ele, nod, u_ndgln:', k, j, u_nod
         end do
         !end if

         if( ndgln_switch )then
            element_nodes => ele_nodes( positions, k )
            !    do j = 1, size( element_nodes )
            do j = 1, xu_nloc
               xu_nod = xu_nod + 1
               xu_ndgln( ( k - 1 ) * xu_nloc + j ) = element_nodes( j )
               !       xu( xu_nod ) = positions % val( 1, xu_nod )
               if( .true. ) then
                  xu( element_nodes( j ) ) = positions % val( 1, element_nodes( j ) )
                  if( ndim > 1 ) yu( element_nodes( j ) ) = positions % val( 2, element_nodes( j ) )
                  if( ndim > 2 ) zu( element_nodes( j ) ) = positions % val( 3, element_nodes( j ) )
               else
                  xu( xu_nod ) = positions % val( 1, element_nodes( j ) )
                  if ( ndim > 1 ) yu( xu_nod ) = positions %val( 2, element_nodes( j ) )
                  if ( ndim > 2 ) zu( xu_nod ) = positions %val( 3, element_nodes( j ) ) 
               end if
               !ewrite(3,*) 'ele, xu_iloc, xu_ndgln :', k, j, xu_ndgln( ( k - 1 ) * xu_nloc + j )
               !ewrite(3,*) 'xu_nod, xu, yu:', element_nodes( j ), xu( element_nodes( j ) ), yu( element_nodes( j ) )
            end do
            if( .not. .true. )then
               if( xu_nloc /= 1 ) xu_nod = xu_nod - 1
            end if
         else
            do j = 1, xu_nloc
               xu_nod = xu_nod + 1
               xu_ndgln((k-1)*xu_nloc+j) = xu_nod
               xu( xu_nod ) = positions%val(1,xu_nod)
               if (ndim>1) yu(xu_nod) = positions%val(2,xu_nod)
               if (ndim>2) zu(xu_nod) = positions%val(3,xu_nod)
            end do
            if( xu_nloc /= 1 ) xu_nod = xu_nod - 1
         end if
      end do Loop_OverElements1

      p_ndgln = cv_ndgln

      if ( ndgln_switch) then
         x_ndgln = cv_ndgln

         ewrite(3,*)' '
         do k = 1, totele
            do j = 1, x_nloc
               ewrite(3,*)'ele, loc, x_ndgln:', k, j, x_ndgln( ( k - 1 ) * x_nloc + j )
            end do
         end do
      end if

      !x_nonods = max( ( x_nloc - 1 ) * totele + 1, totele )
      !allocate( x( x_nonods ) )
      !allocate( y( x_nonods ) )
      !allocate( z( x_nonods ) )
      !x = 0.
      !y = 0.
      !z = 0.

      if ( ndgln_switch) then
         x_nonods = size( positions % val( 1, : ) )
         allocate( x( x_nonods ) )
         allocate( y( x_nonods ) )
         allocate( z( x_nonods ) )
         x = 0.
         y = 0.
         z = 0.
         allocate( index( 1 : x_nonods ) )
         allocate( x_temp( 1 : x_nonods ))

         x = positions % val( 1, : )
         if( ndim > 1 ) y = positions % val( 2, : )
         if( ndim > 2 ) z = positions % val( 3, : )
         do k = 1, totele
            element_nodes => ele_nodes( pressure, k )
            do j = 1, size( element_nodes )
               x_ndgln( ( k - 1 ) * size( element_nodes ) + j ) = element_nodes( j )
               ewrite(3,*)'ele, loc, x_ndgln2:', k, j, element_nodes( j )
            end do
         end do

         if ( .false. ) then
            x = cv_positions % val( 1, : )
            x_temp = x
            index = 0
            call qsort( x_temp, index )
            x_nod = 0
            do j = 1, x_nonods 
               x_nod = x_nod + 1
               x( x_nod ) = x_temp( index( j ))
            end do
            if ( ndim > 1 ) then
               y = cv_positions % val( 2, : )
               x_temp = y
               index = 0
               call qsort( x_temp, index )
               x_nod = 0
               do j = 1, x_nonods 
                  x_nod = x_nod + 1
                  y( x_nod ) = x_temp( index( j ))
               end do
            endif
            if ( ndim > 2 ) then
               z = cv_positions % val( 3, : )
               x_temp = z
               index = 0
               call qsort( x_temp, index )
               x_nod = 0
               do j = 1, x_nonods 
                  x_nod = x_nod + 1
                  z( x_nod ) = x_temp( index( j ))
               end do
            end if
            x_nod=0
            do k = 1, totele
               do j = 1, x_nloc
                  x_nod = x_nod + 1
                  x_ndgln( ( k - 1 ) * x_nloc + j ) = x_nod
               end do
               if( x_nloc /= 1 ) x_nod = x_nod - 1
            end do
         end if
      else ! Do it the old way (just for 1D)
         x_nonods = max( ( x_nloc - 1 ) * totele + 1, totele )
         allocate( x( x_nonods ) )
         allocate( y( x_nonods ) )
         allocate( z( x_nonods ) )
         x = 0.
         y = 0.
         z = 0.
         if( ndim == 1 ) then
            x_nod=0
            x = 0
            do i=1,totele
               do j = 1, x_nloc
                  x_nod = x_nod + 1
                  x_ndgln( ( i - 1 ) * x_nloc + j ) = x_nod
                  if ( x_nloc == 1 ) then
                     x( x_nod ) = ( real( i - 1 ) + 0.5 ) * dx 
                  else
                     x( x_nod ) = real( i - 1 ) * dx + real( j - 1 ) * dx / real ( x_nloc - 1 )
                  end if
               end do
               if( x_nloc /= 1 ) x_nod = x_nod - 1
            end do
            ewrite(3,*) 'x:', size( x ), x
            ewrite(3,*) 'x_ndgln2: ', size( x_ndgln ), x_ndgln
         end if
      end if

      !-
      !- u_ele_type = cv_ele_type = p_ele_type will flag the dimension and 
      !- type of element:
      !- = 1 or 2: 1D (linear and quadratic, respectively)
      !- = 3 or 4: triangle (linear or quadratic, respectively)
      !- = 5 or 6: quadrilateral (bi-linear or tri-linear, respectively)
      !- = 7 or 8: tetrahedron (linear or quadratic, respectively)
      !- = 9 or 10: hexahedron (bi-linear or tri-linear, respectively)
      !-
      !- In order to determine if the model utilises a overlapping 
      !- formulation, just use in the appropriate subroutine:
      !- is_overlapping = .false.
      !- if ( trim( vel_element_type ) == 'overlapping' ) is_overlapping = .true. 
      !-
      !- We need a mechanism to flag this throughout the code, for testing
      !- the current 1D cases. For now:
      !-
      Select Case( ndim )
      case( 1 )
         p_ele_type = 2
         cv_ele_type = 2
         u_ele_type = 2
      case( 2 ) ! Hack ... just linear triangle
         p_ele_type = 3
         cv_ele_type = 3
         u_ele_type = 3
      case default; FLAbort("Incorrect number of dimension - ndim")
      end Select
      !   if( .true. ) then
      !      p_ele_type = 2
      !      cv_ele_type = 2
      !      u_ele_type = 2
      !   else
      !      p_ele_type = 1
      !      cv_ele_type = 1
      !      u_ele_type = 1
      !   end if

      !- 
      !- The follow will disapear soon as they are never used:
      !- mat_ele_type, u_sele_type, cv_sele_type
      mat_ele_type = 1 ; u_sele_type = 1 ; cv_sele_type = 1

      allocate( cv_sndgln( stotel * cv_snloc ) )
      allocate( p_sndgln( stotel * p_snloc ) )
      allocate( u_sndgln( stotel * u_snloc ) )

      Mapping_Surfaces: Select Case( ndim ) ! THIS NEED TO BE CHANGED
      case( 1 )
         cv_sndgln( 1 ) = 1
         cv_sndgln( 2 ) = cv_ndgln( size( cv_ndgln ) )
         p_sndgln( 1 ) = 1
         p_sndgln( 2 ) = p_ndgln( size( p_ndgln ) )
         ! if( u_ele_type == 2 ) then
         if( is_overlapping ) then
            u_nloc2 = u_nloc / cv_nloc
            do j = 1, cv_nloc
               u_sndgln( j ) = 1 + ( j-1 ) * u_nloc2
               u_sndgln( j + cv_nloc ) = u_nonods - u_nloc + j * u_nloc2
            end do
         else
            u_sndgln( 1 ) = 1
            u_sndgln( 2 ) = u_nonods
         end if
      case( 2 : 3 )
         ewrite(3,*)'stotel, cv_snloc:', stotel, cv_snloc
         call getsndgln( pmesh, cv_sndgln )
         do j = 1, stotel
            ewrite(3,*)'cv_sndgln:', j, ( cv_sndgln( ( j - 1 ) * cv_snloc + k ), k = 1, cv_snloc )
         end do

         ewrite(3,*)'stotel, u_snloc:', stotel, u_snloc
         call getsndgln( vmesh, u_sndgln )
         do j = 1, stotel
            ewrite(3,*)'u_sndgln:', j, ( u_sndgln( ( j - 1 ) * u_snloc + k ), k = 1, u_snloc )
         end do
      case default
         FLAbort("Don't have surface global node numbers for this dimension, i.e., > 3")
      end select Mapping_Surfaces

      !! EoS things are going to be done very differently
      !! Currently the default value of ncoef is 10 so I'm going
      !! to put that here until we have a more concrete idea of what
      !! we're dealing with
      ncoef = 10

      !! I don't understand why absorption needs a coefficient, or what it is
      !! but it's equal to 1 in all the test cases
      nuabs_coefs = 1

      ewrite(3,*) ' Getting iteration info'
      call get_option( '/timestepping/nonlinear_iterations', nonlinear_iterations, &
           default = 3 )
      nits = nonlinear_iterations

      if( have_option( '/timestepping/nonlinear_iterations/tolerance' )) then
         call get_option( '/timestepping/nonlinear_iterations/tolerance', &
              nonlinear_iteration_tolerance )
      else
         nonlinear_iteration_tolerance = 1.e-6
      end if

      ! This one is only for compositional problems
      call get_option('/material_phase[0]/scalar_field::component1/&
           &prognostic/temporal_discretisation/control_volumes/number_advection_iterations',&
           &nits_internal, default=1)

      ! a mystery
      ndpset = 0
      noit_dim = 5

      ! Get the number of advection CV face value iterations.
      call get_option("/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/" // &
           "temporal_discretisation/control_volumes/number_advection_iterations", nits_flux_lim_volfra, default=3)

      call get_option("/material_phase[" // int2str(nphases) // "]/scalar_field::ComponentMassFractionPhase1/" // &
           "temporal_discretisation/control_volumes/number_advection_iterations", nits_flux_lim_comp, default=3)

      call get_option("/material_phase[0]/scalar_field::Temperature/prognostic/" // &
           "temporal_discretisation/control_volumes/number_advection_iterations", nits_flux_lim_t, default=3)

      !! disopt options: going to need to change the schema I think
      !       =0      1st order in space          Theta=specified    UNIVERSAL
      !       =1      1st order in space          Theta=non-linear   UNIVERSAL
      !       =2      Trapezoidal rule in space   Theta=specified    UNIVERSAL
      !       =2 if isotropic limiter then FEM-quadratic & stratification adjust. Theta=non-linear 
      !       =3      Trapezoidal rule in space   Theta=non-linear   UNIVERSAL
      !       =4      Finite elements in space    Theta=specified    UNIVERSAL
      !       =5      Finite elements in space    Theta=non-linear   UNIVERSAL
      !       =6      Finite elements in space    Theta=specified    NONE
      !       =7      Finite elements in space    Theta=non-linear   NONE
      !       =8      Finite elements in space    Theta=specified    DOWNWIND+
      !       =9      Finite elements in space    Theta=non-linear   DOWNWIND+

      if (have_option('/material_phase[0]/scalar_field::Temperature/prognostic/' // &
           'spatial_discretisation/control_volumes/face_value::FiniteElement/limit_face_value')) then
         t_disopt = 8 
      else
         t_disopt = 1
      endif

      ! hard wired option
      u_disopt = 1

      if (have_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/control_volumes/face_value::FirstOrderUpwind')) then
         v_disopt = 0 ! Unless theta=non_linear ???
      elseif (have_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/control_volumes/face_value::Trapezoidal')) then
         v_disopt = 2 ! Unless theta=non_linear ???
      elseif (have_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/control_volumes/face_value::FiniteElement/do_not_limit_face_value')) then
         v_disopt = 6 ! Unless theta=non_linear ???
      elseif (have_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/control_volumes/face_value::FiniteElement/limit_face_value')) then
         v_disopt = 8 !! Unless all the other options, but need to be able to get 8 here
      endif

      t_dg_vel_int_opt = 0
      u_dg_vel_int_opt = 4 ! Not used -- it can be deleted
      v_dg_vel_int_opt = 4
      w_dg_vel_int_opt = 0 ! Not used -- it can be deleted

      ewrite(3,*) ' Getting capillary pressure options'
      if (have_option('/porous_media/multiphase_parameters/cp_A')) then
         ! 1 is the only option available at the moment...
         capil_pres_opt = 1
         ! and even this one doesn't depend on any coefficients!
         ncapil_pres_coef = 0
      else
         capil_pres_opt = 0
         ncapil_pres_coef = 0
      end if

      !! These are not currently used in any of the code
      comp_diffusion_opt = 0
      ncomp_diff_coef = 0

      call get_option('/material_phase[0]/scalar_field::Pressure/prognostic/&
           &atmospheric_pressure', patmos, default=0.0)
      call get_option('/material_phase[0]/scalar_field::Pressure/prognostic/&
           &initial_condition::WholeMesh/constant',p_ini, default=0.0)

      call get_option('/material_phase[0]/scalar_field::Temperature/prognostic/' // &
           'spatial_discretisation/conservative_advection', t_beta, default=0.5)

      call get_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/conservative_advection', v_beta)

      call get_option('/material_phase[0]/scalar_field::Temperature/prognostic/' // &
           'temporal_discretisation/theta', t_theta, default=0.5)

      call get_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'temporal_discretisation/theta', v_theta)
      call get_option('/material_phase[0]/vector_field::Velocity/prognostic/' // &
           'temporal_discretisation/theta', u_theta)

      !! I'm going to get this one from the PhaseVolumeFraction scalar_field
      lump_eqns = have_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/continuous_galerkin/mass_terms/lump_mass_matrix')

      volfra_use_theta_flux = .false.
      volfra_get_theta_flux = .true.
      comp_use_theta_flux = .false.
      comp_get_theta_flux = .true.
      t_use_theta_flux = .false.
      t_get_theta_flux = .true.

      ewrite(3,*) 'Finished stuff from read_scalar'

      ! IN/DG_ELE_UPWIND are options for optimisation of upwinding across faces in the overlapping
      ! formulation. The data structure and options for this formulation need to be added later. 
      in_ele_upwind = 3
      dg_ele_upwind = 3

      ! Calculating Mobility
      ewrite(3,*) "going to get viscosities"

      ! Also have to allocate and initialise Viscosity
      allocate(Viscosity( cv_nonods * nphases ))
      Viscosity = 0.

      viscosity_ph1 => extract_tensor_field(state(1), "Viscosity")

      ! why is mobility needed as a scalar variable?
      if( have_option( "/physical_parameters/mobility" ))then
         call get_option( "/physical_parameters/mobility", Mobility )
         Viscosity( 1 : cv_nonods ) = Viscosity_Ph1%val( 1, 1, 1 )
      elseif ( have_option( "/material_phase[1]/vector_field::Velocity/prognostic/" // &
           "tensor_field::Viscosity/prescribed/value::WholeMesh/" // &
           "isotropic")) then
         Viscosity( 1 : cv_nonods ) = Viscosity_Ph1%val( 1, 1, 1 )
         viscosity_ph2 => extract_tensor_field(state(2), "Viscosity")
         Mobility =  Viscosity_Ph2%val( 1, 1, 1 ) / Viscosity_Ph1%val( 1, 1, 1 )         
      elseif( nphases == 1 ) then
         Mobility = 0. 
      end if

      ewrite(3,*) "Got viscosities & mobility"

!!!
!!! Options below are for the multi-component flow model
!!! 
      call get_option("/material_phase[" // int2str(nstates-ncomps) // &
           "]/is_multiphase_component/alpha_beta", alpha_beta, default=1.0)
      KComp_Sigmoid = have_option("/material_phase[" // int2str(nstates-ncomps) // &
           "]/is_multiphase_component/KComp_Sigmoid")
      Comp_Sum2One = have_option("/material_phase[" // int2str(nstates-ncomps) // &
           "]/is_multiphase_component/Comp_Sum2One")

!!!
!!! Porosity and Permeability: it WILL be necessary to change the permeability as it
!!! is defined in the PC as a tensor with dimension ( totele, ndim, ndim )
!!!
!!! Assuming for now that porosity is constant across an element
      porosity => extract_scalar_field(state, "Porosity")

      allocate(volfra_pore(totele))
      pore_ele_loop: do k = 1,element_count(porosity)
         element_nodes => ele_nodes( porosity, k )

         volfra_pore(k) = porosity%val( element_nodes( 1 ))
      end do pore_ele_loop

      ! Assuming for now that permeability is constant across an element
      if (have_option("/porous_media/scalar_field::Permeability")) then
         permeability => extract_scalar_field(state(1), "Permeability")
         allocate(perm(totele, ndim, ndim))
         perm = 0.
         perm_ele_loop: do k = 1,element_count(permeability)
            element_nodes => ele_nodes(permeability,k)
            !perm(k, 1:ndim, 1:ndim)=permeability%val(element_nodes(1))

            ! fyi this is still wrong but will have to do for now
            forall(i=1:ndim)perm(k, i, i)=permeability%val(element_nodes(1))
         end do perm_ele_loop

      elseif (have_option("/porous_media/tensor_field::Permeability")) then

         option_path="/porous_media/tensor_field::Permeability"
         allocate(perm(totele, ndim, ndim))
         perm = 0.
         if (have_option(trim(option_path)//"/prescribed")) then
            option_path=trim(option_path)//"/prescribed/value[0]"
         else
            FLAbort("Permeability field not prescribed - sort this out")
         endif
         is_isotropic=have_option(trim(option_path)//"/isotropic")
         is_diagonal=have_option(trim(option_path)//"/diagonal")
         is_symmetric=have_option(trim(option_path)//"/anisotropic_symmetric")

         if(is_isotropic) then

            option_path=trim(option_path)//"/isotropic"

            if(have_option(trim(option_path)//"/constant")) then
               call get_option(trim(option_path)//"/constant", const)
               do i=1, ndim
                  perm(1:totele,i,i)=const
               end do
            else if(have_option(trim(option_path)//"/python")) then
               t_permeability => extract_tensor_field(state(1), "Permeability")
               do k = 1, element_count(t_permeability)
                  element_nodes => ele_nodes(t_permeability,k)
                  do i=1,ndim
                     perm(k,i,i) = t_permeability%val(i, i, element_nodes(1))
                  end do
               end do
            else if (have_option(trim(option_path)//"/generic_function")) then
               FLExit("Generic functions are obsolete. Use a Python function.")
            else
               FLExit("Incorrect initial condition for field")
            end if

         else if(is_diagonal) then

            option_path=trim(option_path)//"/diagonal"

            if(have_option(trim(option_path)//"/constant")) then
               allocate(const_vec(ndim))
               const_vec=0.
               call get_option(trim(option_path)//"/constant", const_vec)
               do i=1, ndim
                  perm(1:totele,i,i)=const_vec(i)
               end do
               deallocate(const_vec)
            else if(have_option(trim(option_path)//"/python")) then
               t_permeability => extract_tensor_field(state(1), "Permeability")
               do i=1,ndim
                  do j=1,element_count(t_permeability)
                     perm(j,i,i) = t_permeability%val(i, i, j)
                  end do
               end do
            else if (have_option(trim(option_path)//"/generic_function")) then
               FLExit("Generic functions are obsolete. Use a Python function.")
            else
               FLExit("Incorrect initial condition for field")
            end if

         else

            ! Set path
            if(is_symmetric) then
               option_path=trim(option_path)//"/anisotropic_symmetric"
            else
               option_path=trim(option_path)//"/anisotropic_asymmetric"
            end if

            if(have_option(trim(option_path)//"/constant")) then
               allocate(const_array(ndim, ndim))
               const_array=0.
               call get_option(trim(option_path)//"/constant", const_array)
               do i=1,ndim
                  do j=1,ndim
                     perm(1:totele, i, j) = const_array(i, j)
                  end do
               end do
               deallocate(const_array)
            else if (have_option(trim(option_path)//"/python")) then
               t_permeability => extract_tensor_field(state(1), "Permeability")
               do k = 1, element_count(t_permeability)
                  element_nodes => ele_nodes(t_permeability,k)
                  do i=1,ndim
                     do j=1,ndim
                        perm(k,i,j) = t_permeability%val(i, j, element_nodes(1))
                     end do
                  end do
               end do

            else if (have_option(trim(option_path)//"/generic_function")) then
               FLExit("Generic functions are obsolete. Use a Python function.")
            else
               FLExit("Incorrect initial condition for field")
            end if
         end if

      end if

!!!
!!! WIC_X_BC (in which X = D, U, V, W, P, T, COMP and VOL) controls the boundary conditions
!!! type applied. == 1 (Dirichlet), = 2 (Robin), = 3 (Newman) 
!!!
!!!
!!! Components Boundary Conditions


      allocate( wic_comp_bc( stotel * nphases ))
      allocate( suf_comp_bc( stotel * cv_snloc * nphases * ncomps ))
      wic_comp_bc = 0
      suf_comp_bc = 0.

      Loop_Component_BC: do i=nphases, nphases+ncomps-1
         do j=1,nphases
            if( have_option("/material_phase[" // int2str(i) // "]/scalar_field::ComponentMassFractionPhase" // &
                 int2str(j) //"/prognostic/" // "boundary_conditions[0]/type::dirichlet" )) then

               shape_option=option_shape("/material_phase[" // int2str(i) // &
                    "]/scalar_field::ComponentMassFractionPhase" // &
                    int2str(j) //"/prognostic/boundary_conditions[0]/surface_ids")

               allocate(component_sufid_bc(1:shape_option(1)))
               Component_BC_Type = 1
!               nobcs = get_boundary_condition_count(  )

               call get_option( "/material_phase[" // int2str(i) // "]/scalar_field::ComponentMassFractionPhase" // & 
                    int2str(j) //"/prognostic/" // "boundary_conditions[0]/surface_ids", component_sufid_bc )

               call get_option( "/material_phase[" // int2str(i) // "]/scalar_field::ComponentMassFractionPhase" // &
                    int2str(j) //"/prognostic/" // "boundary_conditions[0]/type::dirichlet/constant", Component_Suf_BC )

               do k=1,shape_option(1)
                  wic_comp_bc( component_sufid_bc(k) + nphases*(j-1) ) = Component_BC_Type
                  suf_comp_bc( component_sufid_bc(k) + nphases*(j-1) + nphases*ncomps*(i-nphases) ) = Component_Suf_BC
               enddo
               deallocate(Component_sufid_bc)
            endif

         end do

      end do Loop_Component_BC

      allocate( suf_cpd_bc( stotel * cv_snloc * nphases ))
      suf_cpd_bc = 0.      

      ! Extra allocation to stop the code dying when also running everything through the old io
      allocate( suf_one_bc( stotel * cv_snloc * nphases ))
      suf_one_bc=0.

      ewrite(3,*) "Done with boundary conditions"

!!!
!!! Robin Boundary conditions need to be added at a later stage
!!!
      allocate( suf_u_bc_rob1( stotel * u_snloc * nphases ))
      allocate( suf_v_bc_rob1( stotel * u_snloc * nphases ))
      allocate( suf_w_bc_rob1( stotel * u_snloc * nphases ))
      allocate( suf_u_bc_rob2( stotel * u_snloc * nphases ))
      allocate( suf_v_bc_rob2( stotel * u_snloc * nphases ))
      allocate( suf_w_bc_rob2( stotel * u_snloc * nphases ))
      allocate( suf_t_bc_rob1( stotel * nphases ))
      allocate( suf_t_bc_rob2( stotel * nphases ))
      allocate( suf_vol_bc_rob1( stotel * nphases ))
      allocate( suf_vol_bc_rob2( stotel * nphases ))
      allocate( suf_comp_bc_rob1( stotel * nphases ))
      allocate( suf_comp_bc_rob2( stotel * nphases ))
      suf_u_bc_rob1 = 0.
      suf_u_bc_rob2 = 0.
      suf_v_bc_rob1 = 0.
      suf_v_bc_rob2 = 0.
      suf_w_bc_rob1 = 0.
      suf_w_bc_rob2 = 0.
      suf_t_bc_rob1 = 0.
      suf_t_bc_rob2 = 0.
      suf_vol_bc_rob1 = 0.
      suf_vol_bc_rob2 = 0.
      suf_comp_bc_rob1 = 0.
      suf_comp_bc_rob2 = 0.

!!!
!!!  End of Boundary Conditions Section
!!!

!!!
!!!  Density
!!!
      ewrite(3,*) "going to get density..."

      Loop_Density: do i = 1, nphases
         density => extract_scalar_field( state( i ), "Density")
         if ( .not. allocated( den )) allocate( den( nphases * node_count( density )))

         do j = node_count( density ), 1, -1
            den(( i - 1 ) * node_count( density ) + j ) = density%val( j )
            ! This will make sure that the fields in the PC *does not* contain any boundary conditions
            ! elements. This need to be changed along with the future data structure to take into 
            ! account the overlapping formulation.
            if( j == 1 ) den(( i - 1 ) * node_count( density ) + j ) = &
                 den(( i - 1 ) * node_count( density ) + j + 1)     
         end do

         if( .not. ( allocated( wic_d_bc ) .and. allocated( suf_d_bc ))) then 
            allocate( wic_d_bc( stotel * nphases ))
            allocate( suf_d_bc( stotel * p_snloc * nphases ))
            wic_d_bc = 0
            suf_d_bc = 0.
         end if



      Conditional_Density_BC: if( have_option( "/material_phase[" // int2str(i-1) // &
              "]/scalar_field::Density/" // &
              "prognostic/boundary_conditions[0]/type::dirichlet" )) then

            shape_option = option_shape( "/material_phase[" // int2str(i-1) // "]/scalar_field::" // &
                 "Density/prognostic/boundary_conditions[0]/surface_ids" )
            allocate(  density_sufid_bc( 1 : shape_option( 1 )))

            pvf_bc_type = 1
            nobcs = get_boundary_condition_count( density )


           do k = 1, nobcs
             density_bc => extract_surface_field( density, k, "value" )
              
             Density_BC_Type = 1

             call get_option( '/material_phase[' // int2str(i-1) // &
                 ']/scalar_field::Density/prognostic/' // &
                 'boundary_conditions[0]/surface_ids', Density_SufID_BC )


             ! Specify the boundary condition type.
             do j=1,shape_option(1)
               wic_d_bc( density_sufid_bc( j ) + ( i - 1 ) * nphases ) = density_bc_type
             enddo
      
             do j = 1, stotel
              if(pmesh%faces%boundary_ids(j) .eq. density_sufid_bc(k)) then
               do kk = 1, p_snloc
                  suf_d_bc( ( i - 1 ) * stotel*p_snloc + (j-1)*p_snloc + kk ) = density_bc%val( k )
               end do
              end if 
             end do
      
           end do ! End of BC loop


          deallocate(density_sufid_bc)

         endif Conditional_Density_BC

      enddo Loop_Density

!!!
!!! Pressure
!!!
      ewrite(3,*) "pressure..."
      pressure => extract_scalar_field( state( 1 ), "Pressure" )
      allocate( p( node_count( pressure )))

      ! This will make sure that the fields in the PC *does not* contain any boundary conditions
      ! elements. This need to be changed along with the future data structure to take into 
      ! account the overlapping formulation.
      do i= node_count( pressure ), 1, -1
         p( i ) = pressure%val( i )
         if( i == 1 ) p ( i ) = p( i + 1 )
      enddo

      !! Control-volume pressure used for shock tube initialisation
      allocate( cv_p( node_count( pressure )))
      cv_p = p

      

      if( .not. ( allocated( wic_p_bc ) .and. allocated( suf_p_bc ))) then
          allocate( wic_p_bc( stotel * nphases ))
          allocate( suf_p_bc( stotel * p_snloc * nphases ))
          wic_p_bc = 0
          suf_p_bc = 0.
      end if

      
      i = 1 ! The first phase is set, and the others alias to this phase for pressure.   
      Conditional_Pressure_BC: if( have_option( '/material_phase[' // int2str(i-1) // &
            ']/scalar_field::Pressure/prognostic/' // &
            'boundary_conditions[0]/type::dirichlet' )) then

           shape_option=option_shape('/material_phase[' // int2str(i-1) // &
               ']/scalar_field::Pressure/' // &
               'prognostic/boundary_conditions[0]/surface_ids')

           allocate( pressure_sufid_bc( 1 : shape_option( 1 )))


           nobcs = get_boundary_condition_count( pressure )
           do k = 1, nobcs
             pressure_bc => extract_surface_field( pressure, k, "value" )

             Pressure_BC_Type = 1
             call get_option( '/material_phase[' // int2str(i-1) // &
               ']/scalar_field::Pressure/prognostic/' // &
               'boundary_conditions[0]/surface_ids', Pressure_SufID_BC )

             ! Specify the boundary condition type.
             do j = 1, shape_option( 1 )
               wic_p_bc( pressure_sufid_bc( j ) + ( i - 1 ) * nphases ) = pressure_bc_type
             enddo
      
             do j = 1, stotel
              if(pmesh%faces%boundary_ids(j) .eq. pressure_sufid_bc(k)) then
               do kk = 1, p_snloc
                  suf_p_bc( ( i - 1 ) * stotel*p_snloc + (j-1)*p_snloc + kk ) = pressure_bc%val( k )
               end do
              end if 
             end do
      
           end do ! End of BC loop
             
         end if Conditional_Pressure_BC
         
         
        ! Copy this to the other phases.
        if(nphases > 1) then
         do i = 2, nphases
          wic_p_bc( pressure_sufid_bc( 1 ) + ( i - 1 ) * nphases ) = Pressure_BC_Type
          do j = 1, stotel
           do k = 1, p_snloc
             suf_p_bc( ( i - 1 ) * stotel + (j-1)*p_snloc + k ) =  suf_p_bc( ( 1 - 1 ) * stotel + (j-1)*p_snloc + k )
           end do
          end do    
         end do 
        end if
     
      
!!!
!!! Volume Fraction (or Saturation) and associated boundary conditions:
!!!
      ewrite(3,*) "phasevolumefraction..."

      Loop_VolumeFraction: do i = 1, nphases

         phasevolumefraction => extract_scalar_field(state(i), "PhaseVolumeFraction")

         if ( .not. allocated( satura )) allocate( satura( nphases * node_count( phasevolumefraction )))

         ! This will make sure that the fields in the PC *does not* contain any boundary conditions
         ! elements. This need to be changed along with the future data structure to take into 
         ! account the overlapping formulation.
         !do j = node_count( phasevolumefraction ), 1, -1
         !   satura( ( i - 1 ) * node_count( phasevolumefraction ) + j ) = phasevolumefraction%val( j )
         !   if( j == 1 ) satura( ( i - 1 ) * node_count( phasevolumefraction ) + j ) = &
         !                satura( ( i - 1 ) * node_count( phasevolumefraction ) + j + 1)
         !enddo

         ! this assumes a constant initial condition for now
         ! but does not take into account any boundary conditions
         ! and is valid for multi-d.
         call get_option('/material_phase[' // int2str(i-1) // &
              ']/scalar_field::PhaseVolumeFraction/prognostic/' // &
              'initial_condition[0]/constant', val)

         satura( (i-1)*node_count(phasevolumefraction)+1 : &
              &   i   *node_count(phasevolumefraction) ) = val

         if (.not. allocated(wic_vol_bc)) then
            allocate( wic_vol_bc( stotel * nphases ))
            wic_vol_bc = 0.
         endif

         if (.not. allocated(suf_vol_bc)) then
            allocate( suf_vol_bc( stotel * 1 * nphases ))
            suf_vol_bc = 0.
         endif



         Conditional_VolumeFraction_BC: if( have_option( "/material_phase[" // int2str(i-1) // &
              "]/scalar_field::PhaseVolumeFraction/" // &
              "prognostic/boundary_conditions[0]/type::dirichlet" )) then

            shape_option = option_shape( "/material_phase[" // int2str(i-1) // "]/scalar_field::" // &
                 "PhaseVolumeFraction/prognostic/boundary_conditions[0]/surface_ids" )
            allocate( pvf_sufid_bc( 1 : shape_option( 1 )))

            pvf_bc_type = 1
            nobcs = get_boundary_condition_count( phasevolumefraction )


           do k = 1, nobcs
             phasevolumefraction_bc => extract_surface_field( phasevolumefraction, k, "value" )
              
             pvf_bc_type = 1
             call get_option( "/material_phase[" // int2str(i-1) //"]/scalar_field::" // &
                 "PhaseVolumeFraction/prognostic/boundary_conditions[0]/surface_ids", pvf_sufid_bc )

             ! Specify the boundary condition type.
             do j = 1, shape_option( 1 )
               wic_vol_bc( pvf_sufid_bc(j) + ( i - 1 ) * nphases ) = pvf_bc_type
             enddo
         
      
             do j = 1, stotel
              if(pmesh%faces%boundary_ids(j) .eq. pvf_sufid_bc(k)) then
               do kk = 1, p_snloc
                  suf_vol_bc( ( i - 1 ) * stotel*p_snloc + (j-1)*p_snloc + kk ) = phasevolumefraction_bc%val( k )
               end do
              end if 
             end do
      
           end do ! End of BC loop


          deallocate(pvf_sufid_bc)

         endif Conditional_VolumeFraction_BC

      enddo Loop_VolumeFraction

      ! Also need to allocate and initialise volfra
      allocate(volfra( cv_nonods * nphases ))
      volfra=0.

!!!
!!! Velocity and associated boundary conditions:
!!!
      ewrite(3,*) "velocity..."

      Loop_Velocity: do i = 1, nphases

         velocity => extract_vector_field(state(i), "Velocity")
         option_path = "/material_phase["//int2str(i-1)//"]/vector_field::Velocity"
         if (.not. allocated(u)) then
            allocate(u(nphases*size(velocity%val(1, :))))
            allocate(v(nphases*size(velocity%val(1, :))))
            allocate(w(nphases*size(velocity%val(1, :))))
            u=0.
            v=0.
            w=0.
         endif

         ! facility to initialise velocity field with the initial value from diamond
         ! as long as it's a constant value
         if (.not. allocated(initial_constant_velocity)) allocate(initial_constant_velocity(ndim))
         initial_constant_velocity=0.
         call get_option(trim(option_path)//"/prognostic/initial_condition::WholeMesh/constant", &
              initial_constant_velocity, stat)

         if (stat==0) then
            u((i-1)*size(velocity%val(1, :))+1:i*size(velocity%val(1, :)))=initial_constant_velocity(1)
            if (ndim>1) v((i-1)*size(velocity%val(1, :))+1:i*size(velocity%val(1, :)))=initial_constant_velocity(2)
            if (ndim>2) w((i-1)*size(velocity%val(1, :))+1:i*size(velocity%val(1, :)))=initial_constant_velocity(3)
         endif

         ! This will make sure that the fields in the PC *does not* contain any boundary conditions
         ! elements. This need to be changed along with the future data structure to take into 
         ! account the overlapping formulation.
         do j = size(velocity%val(1, :)), 1, -1
            u((i-1)*size(velocity%val(1, :))+j)=velocity%val(X_, j)
            if( j == 1 ) u( ( i - 1 ) * size(velocity%val(1, :)) + j ) = &
                 u( ( i - 1 ) * size(velocity%val(1, :)) + j + 1)
            if (ndim>1) then
               v((i-1)*size(velocity%val(1, :))+j)=velocity%val(Y_, j)
               if( j == 1 ) v( ( i - 1 ) * size(velocity%val(1, :)) + j ) = &
                    v( ( i - 1 ) * size(velocity%val(1, :)) + j + 1)
            endif
            if (ndim>2) then
               w((i-1)*size(velocity%val(1, :))+j)=velocity%val(Z_, j)
               if( j == 1 ) w( ( i - 1 ) * size(velocity%val(1, :)) + j ) = &
                    w( ( i - 1 ) * size(velocity%val(1, :)) + j + 1)
            endif
         enddo

         if (.not.allocated(wic_u_bc)) then
            allocate( wic_u_bc( stotel * nphases ))
            wic_u_bc = 0
         endif
         if (.not.allocated(suf_u_bc)) then
            allocate( suf_u_bc( stotel * u_snloc * nphases ))
            suf_u_bc = 0.
         endif
         if (.not.allocated(suf_v_bc)) then
            allocate( suf_v_bc( stotel * u_snloc * nphases ))
            suf_v_bc = 0.
         endif
         if (.not.allocated(suf_w_bc)) then
            allocate( suf_w_bc( stotel * u_snloc * nphases ))
            suf_w_bc = 0.
         endif
  

         Conditional_Velocity_BC: if( have_option( '/material_phase[' // int2str(i-1) // &
              ']/vector_field::Velocity/prognostic/' // &
              'boundary_conditions[0]/type::dirichlet' )) then

            shape_option=option_shape("/material_phase[" // int2str(i-1) // "]/vector_field::Velocity/&
                 &prognostic/boundary_conditions[0]/surface_ids")
            if( .not. allocated( velocity_sufid_bc ))allocate(velocity_sufid_bc(1:shape_option(1)))

            Velocity_BC_Type = 1
            nobcs = get_boundary_condition_count( velocity )
            
            do k = 1, nobcs
               velocity_bc => extract_surface_field( velocity, k, "value" )
           
            
                call get_option( "/material_phase[" // int2str(i-1) // "]/vector_field::Velocity/" // &
                 "prognostic/boundary_conditions[0]/surface_ids", Velocity_SufID_BC )
           
!
                 ! Specify the boundary condition type.
                 do j = 1, shape_option( 1 )
                   wic_u_bc( velocity_sufid_bc(j) + ( i - 1 ) * nphases ) = Velocity_BC_Type
                 enddo
      
             do j = 1, stotel
!              print*, vmesh%faces%boundary_ids(j) , velocity_sufid_bc(k)
              if(vmesh%faces%boundary_ids(j) .eq. velocity_sufid_bc(k)) then
               do kk = 1, p_snloc
                  suf_u_bc( ( i - 1 ) * stotel*cv_snloc + (j-1)*cv_snloc + kk ) = velocity_bc%val( 1, k )
                  if (velocity%dim>1) suf_v_bc( ( i - 1 ) * stotel*cv_snloc + (j-1)*cv_snloc + kk ) = velocity_bc%val( 2, k )
                  if (velocity%dim>2) suf_w_bc( ( i - 1 ) * stotel*cv_snloc + (j-1)*cv_snloc + kk ) = velocity_bc%val( 3, k )
!                  print*, suf_u_bc( ( i - 1 ) * stotel*cv_snloc + (j-1)*cv_snloc + kk ), j
               end do
              end if 
           
            end do
            
           end do ! End of BC Loop 
            
         endif Conditional_Velocity_BC
         
!         print*, suf_u_bc, wic_u_bc


      enddo Loop_Velocity


      allocate(uabs_option(nphases))
      allocate(eos_option(nphases))
      allocate(cp_option(nphases))
      do i=1,nphases
         option_path = "/material_phase[" // int2str(i-1) // "]/multiphase_properties/relperm_type"
         if (have_option(trim(option_path)//"/Corey")) uabs_option(i)=3
         if (have_option(trim(option_path)//"/Corey/boost_at_zero_saturation")) uabs_option(i)=5

         option_path = "/material_phase[" // int2str(i-1) // "]/equation_of_state"
         if (have_option(trim(option_path)//"/incompressible/linear")) then
            eos_option(i) = 2
         elseif (have_option(trim(option_path)//"/compressible/stiffened_gas")) then
            eos_option(i) = 1
         elseif (have_option(trim(option_path)//"/compressible/exponential_oil_gas")) then
            eos_option(i) = 3
         else
            FLAbort("Unknown EoS option for phase "// int2str(i))
         endif
         cp_option(i) = 1
      enddo

      !! uabs_coefs is currently only used in rel perm options which aren't
      !! selected in any of the test cases, ie the 'standard polynomial
      !! representation of relative permeability'
      allocate(uabs_coefs(nphases, nuabs_coefs))
      allocate(eos_coefs(nphases, ncoef))
      eos_coefs=0.
      !! Capillary pressure isn't used at all at the moment
      allocate(cp_coefs(nphases, nphases))
      do i=1,nphases
         uabs_coefs(i,1) = 1.
         if (eos_option(i)==2) then
            if (have_option("/material_phase[" // int2str(i-1) // "]/equation_of_state/incompressible/linear/all_equal")) then
               call get_option("/material_phase[" // int2str(i-1) // "]/equation_of_state/incompressible/linear/all_equal", eos_value)
               eos_coefs(i,1) = eos_value
            else
               call get_option("/material_phase[" // int2str(i-1) // "]/equation_of_state/incompressible/linear/specify_all", eos_coefs(i, :))
            endif
         elseif (eos_option(i)==1) then
            call get_option("/material_phase[" // int2str(i-1) // "]/equation_of_state/compressible/stiffened_gas/eos_option1", eos_coefs(i, 1))
            call get_option("/material_phase[" // int2str(i-1) // "]/equation_of_state/compressible/stiffened_gas/eos_option2", eos_coefs(i, 2))
            eos_coefs(i, 3:ncoef) = 0.
         endif

      enddo
      cp_coefs = 1.

      ewrite(3,*) 'Getting source terms -- gravity '
      ! Gravity is associated with the u_source term
      call get_option( "/physical_parameters/gravity/magnitude", gravity_magnitude, stat )
      have_gravity = ( stat == 0 )

      if( have_gravity ) then
         gravity_direction = extract_vector_field(state(1), 'GravityDirection', stat )
         ! Normalise direction vector
         grm=0
         do i=1,ndim
            grm=grm + gravity_direction%val(i,1)**2
         end do
         do i=1,ndim
            gravity_direction%val(i,:) = gravity_direction%val(i,:)/sqrt(grm)
         end do
      end if

      !!if( have_gravity ) then
      !!   if( have_option( '/physical_parameters/gravity/vector_field::' // &
      !!        'GravityDirection/prescribed/value::WholeMesh' ))then
      !!      call get_option( '/physical_parameters/gravity/vector_field::' // &
      !!           'GravityDirection/prescribed/value::WholeMesh/constant', &
      !!           gravity_direction, stat )
      !!   end if
      !!end if

      ewrite(3, *)"Getting source terms -- velocity "
      allocate( u_source( ndim * u_nonods * nphases ))
      u_source = 0.
      Conditional_VelocitySource: if( have_option( '/material_phase[0]/vector_field::Velocity/' // &
           'prognostic/vector_field::Source' )) then 
         ! This is still not working as the length of node_count(velocity_source) =
         ! node_count(velocity) /= u_nonods
         do i=1,nphases
            velocity_source => extract_vector_field(state(i), "VelocitySource", stat)
            if (stat==0) then
               do j=1,node_count(velocity_source)
                  u_source((i-1)*node_count(velocity_source)+j)=velocity_source%val(X_, j)
                  if (ndim>1) u_source(u_nonods*nphases + (i-1)*node_count(velocity_source)+j) = velocity_source%val(Y_, j)
                  if (ndim>2) u_source(2*u_nonods*nphases + (i-1)*node_count(velocity_source)+j) = velocity_source%val(Z_, j)                  
               enddo
            else
               u_source = 0.
            endif
         enddo
      end if Conditional_VelocitySource

      if (have_gravity) then
         do i = 1, nphases - 1, 1
            delta_den = 0.
            do j = 2, node_count( density )
               delta_den = delta_den + ( den((i)*node_count(density)+j) - &
                    den((i-1)*node_count(density)+j) ) / &
                    real( ( node_count( density ) - 1 ) * max( 1, ( nphases - 1 )))
            end do
            do j = 1, u_nonods
               do k = 1, ndim
                  u_source( ( k - 1 ) * u_nonods * nphases + ( i - 1 ) * u_nonods + j  ) = &
                       u_source( ( k - 1 ) * u_nonods * nphases + ( i - 1 ) * u_nonods + j  ) + &
                       delta_den * gravity_magnitude * gravity_direction%val(k,1) * min( 0.027, & 
                       domain_length /  real( totele * cv_nloc ))
               end do
            end do
         end do
      end if

      ewrite(3,*) 'Getting PVF Source'
      do i=1,nphases
         pvf_source => extract_scalar_field(state(i), "PhaseVolumeFractionSource", stat)
         if (.not.allocated(v_source)) allocate(v_source(cv_nonods*nphases))
         if (stat==0) then
            do j=1,node_count(pvf_source)
               v_source((i-1)*node_count(pvf_source)+j)=pvf_source%val(j)
            enddo
         else
            v_source = 0.
         endif
      enddo

      ewrite(3,*) 'Getting temperature field'

      allocate( t( cv_nonods * nphases ))
      t = 0.0

      Conditional_ExtraScalarField: if( have_option( "/material_phase[0]/" // &
           "scalar_field::Temperature" ))then

         have_temperature_fields = .true.

         Loop_Temperature: do i = 1, nphases
            scalarfield => extract_scalar_field(state(i), "Temperature")

            t_ele_loop: do k = 1,element_count(scalarfield)

               element_nodes => ele_nodes(scalarfield,k)

               t_node_loop: do j = 1,size(element_nodes)

                  t((cv_ndgln((k-1)*size(element_nodes)+j)) + (i-1)*node_count(scalarfield)) = &
                       scalarfield%val(element_nodes(j))

               end do t_node_loop

            end do t_ele_loop

            scalarfield_source => extract_scalar_field( state( i ), trim( field_name ) // &
                 "Source", stat)

            if ( .not. allocated( t_source )) allocate( t_source( cv_nonods * nphases ))
            if ( stat == 0 ) then
               do j = 1, node_count( scalarfield_source )
                  t_source(( i - 1 ) * node_count( scalarfield_source) + j ) = &
                       scalarfield_source%val( j )
               enddo
            else
               t_source = 0.
            endif

            if( ( .not. allocated( wic_t_bc )) .and. ( .not. allocated( suf_t_bc ))) then
               allocate( wic_t_bc( stotel * nphases ))
               allocate( suf_t_bc( stotel * cv_nloc * nphases ))
               wic_t_bc = 0.
               suf_t_bc = 0.
            endif



            Conditional_Temperature_BC: if( have_option( "/material_phase[" // int2str(i-1) // &
                 "]/scalar_field::Temperature/prognostic/boundary_conditions[0]/" // &
                 "type::dirichlet" )) then

               shape_option=option_shape("/material_phase[" // int2str(i-1) // &
                    "]/scalar_field::Temperature/prognostic/boundary_conditions[0]/" // &
                    "surface_ids")

               if( .not. allocated( Temperature_sufid_bc )) &
                    allocate( Temperature_sufid_bc( 1 : shape_option( 1 )))
               Temperature_bc_type = 1

               call get_option( "/material_phase[" // int2str(i-1) // "]/scalar_field::" // &
                    "Temperature/prognostic/" // &
                    "boundary_conditions[0]/surface_ids", Temperature_sufid_bc )

               do j=1,shape_option(1)
                  wic_t_bc( Temperature_sufid_bc(j) + (i-1)*nphases ) = Temperature_bc_type
               enddo

               nobcs = get_boundary_condition_count( scalarfield )


               do k = 1, nobcs
                
                scalarfield_bc => extract_surface_field( scalarfield, k, "value" )
              
                 ! Specify the boundary condition type.
                do j=1,shape_option(1)
                 wic_t_bc( Temperature_sufid_bc( j ) + ( i - 1 ) * nphases ) = Temperature_bc_type
                enddo
       
                do j = 1, stotel
                 if(vmesh%faces%boundary_ids(j) .eq. density_sufid_bc(k)) then
                  do kk = 1, p_snloc
                   suf_t_bc( ( i - 1 ) * stotel*cv_snloc + (j-1)*cv_snloc + kk ) = scalarfield_bc%val( k )
                  end do
                 end if 
               end do
      
              end do ! End of BC loop

              deallocate(Temperature_sufid_bc)

            endif Conditional_Temperature_BC

         end do Loop_Temperature

      else 

         have_temperature_fields = .false.

      end if Conditional_ExtraScalarField

      ewrite(3,*) 'Getting component source'
      allocate(comp_source(cv_nonods*nphases))
      comp_source=0.
      ! comp is stored in the order
      !   comp1 phase1
      !   comp1 phase2
      !   comp2 phase1
      !   comp2 phase2
      !   etc
      ! Component order doesn't really matter but the material_phase names
      ! should match up
      if (ncomps>0) then
         allocate(comp(cv_nonods*nphases*ncomps))
         ! Assume for now that components have been inserted in state AFTER all the phases
         do i=nstates-ncomps+1,nstates
            do j=1,option_count("/material_phase[" // int2str(i-1) //"]/scalar_field")
               componentmassfraction => extract_scalar_field(state(i), j)
               if (componentmassfraction%name(1:21) == "ComponentMassFraction") then
                  call get_option("/material_phase[" // int2str(i-1) //"]/scalar_field[" // int2str(j-1) //&
                       &"]/material_phase_name", material_phase_name)
                  do k=1,nphases
                     ewrite(3,*) 'mp_name, state_name', trim(material_phase_name), state(k)%name
                     if (trim(material_phase_name) == trim(state(k)%name)) then
                        do l=1,node_count(componentmassfraction)
                           comp( ((i-(1+nphases))*nphases+(k-1))*node_count(componentmassfraction)+l ) = componentmassfraction%val(l)
                        enddo
                     endif
                  enddo
               endif
            enddo
         enddo
      endif

      ewrite(3,*) 'Getting capillary pressure options and absorptions'
      if (ncapil_pres_coef>0) then
         allocate(capil_pres_coef(ncapil_pres_coef,nphases,nphases))
         capil_pres_coef=0.
      endif
      ! Not sure what this one does yet
      allocate(u_abs_stab(cv_nloc*totele, ndim*nphases, ndim*nphases))
      u_abs_stab=0.
      ! or this one
      allocate(u_absorb(cv_nloc*totele, ndim*nphases, ndim*nphases))
      u_absorb=0.
      allocate(comp_absorb(cv_nonods, nphases, nphases))
      comp_absorb=0.
      allocate( t_absorb( cv_nonods, nphases, nphases ))
      t_absorb=0.
      allocate( v_absorb( cv_nonods, nphases, nphases ))
      v_absorb=0.
      allocate( k_comp( ncomps, nphases, nphases ))
      k_comp=0.
      if (KComp_Sigmoid) then
         do i=1, ncomps
            call get_option('material_phase['// int2str(i+nphases-1) //']/is_multiphase_component/' // &
                 'KComp_Sigmoid/k_comp', k_comp(i, 1, 1))
            k_comp(i, 1:nphases, 1:nphases) = k_comp(i, 1, 1)
         end do
      end if

      allocate( comp_diffusion( cv_nloc*totele, ndim, ndim, nphases ))
      comp_diffusion=0.0
      allocate( comp_diff_coef( ncomps, ncomp_diff_coef, nphases ))
      comp_diff_coef=0.0
      allocate( cv_one( nphases * cv_nonods ))
      cv_one = 1.0

      ewrite(3,*) "Leaving copy_outof_state"
      
!      stop 999

    end subroutine copy_outof_state

    subroutine copy_into_state(state, &
         proto_saturations, &
         proto_temperatures, &
         proto_pressure, &
         proto_velocity_u, &
         proto_velocity_v, &
         proto_velocity_w, &                               
         proto_densities, &
         proto_components, &
         ncomp, &
         nphase, &
         cv_ndgln, &
         p_ndgln, &
         u_ndgln, &
         ndim)

      !!< Copy prototype solution arrays into fluidity state array for output

      type(state_type), dimension(:), intent(inout) :: state
      real, dimension(:), intent(in) :: proto_saturations
      real, dimension(:), intent(in) :: proto_temperatures      
      real, dimension(:), intent(in) :: proto_pressure
      real, dimension(:), intent(in) :: proto_velocity_u
      real, dimension(:), intent(in) :: proto_velocity_v
      real, dimension(:), intent(in) :: proto_velocity_w      
      real, dimension(:), intent(in) :: proto_densities
      real, dimension(:), intent(in) :: proto_components
!      real, dimension(:,:,:), intent(in) :: velocity_dg
      integer, intent(in) :: ncomp      
      integer, intent(in) :: nphase
      integer, dimension(:), intent(in) :: cv_ndgln
      integer, dimension(:), intent(in) :: p_ndgln
      integer, dimension(:), intent(in) :: u_ndgln
      integer, intent(in) :: ndim

      ! local variables        
      integer :: stat
      integer :: i,j,k,p,ele,jloc
      integer :: number_nodes
      integer :: nstates
      integer :: nloc
      integer, dimension(:), pointer :: element_nodes => null()
      type(scalar_field), pointer :: phasevolumefraction => null()
      type(scalar_field), pointer :: phasetemperature => null()      
      type(scalar_field), pointer :: pressure => null()
      type(vector_field), pointer :: velocity => null()
      type(scalar_field), pointer :: density => null()
      type(scalar_field), pointer :: componentmassfraction => null()
      character(len=option_path_len) :: material_phase_name

      ewrite(3,*) "In copy_into_state"

      assert(size(state) >= nphase)

      phase_loop: do p = 1,nphase

         phasevolumefraction => extract_scalar_field(state(p), "PhaseVolumeFraction", stat=stat)

         if (stat /= 0) then 

            ewrite(1,*) 'Issue in prototype interface for phase ',p

            FLAbort('Failed to extract phase volume fraction from state in copy_into_state')

         end if

         number_nodes = node_count(phasevolumefraction)

         volf_ele_loop: do i = 1,element_count(phasevolumefraction)

            element_nodes => ele_nodes(phasevolumefraction,i)

            nloc = size(element_nodes)

            volf_node_loop: do j = 1,nloc

               call set(phasevolumefraction, &
                    element_nodes(j), &
                    proto_saturations((cv_ndgln((i-1)*nloc+j)) + (p-1)*number_nodes))

            end do volf_node_loop

         end do volf_ele_loop

         phasetemperature => extract_scalar_field(state(p), "Temperature", stat=stat)

         found_temp: if (stat == 0) then 

            ewrite(1,*) 'In copy in to state and found temperature for phase ',p

            number_nodes = node_count(phasetemperature)

            temp_ele_loop: do i = 1,element_count(phasetemperature)

               element_nodes => ele_nodes(phasetemperature,i)

               nloc = size(element_nodes)

               temp_node_loop: do j = 1,nloc

                  call set(phasetemperature, &
                       element_nodes(j), &
                       proto_temperatures((cv_ndgln((i-1)*nloc+j)) + (p-1)*number_nodes))

               end do temp_node_loop

            end do temp_ele_loop

         end if found_temp

         density => extract_scalar_field(state(p), "Density", stat=stat)

         if (stat /= 0) then 

            ewrite(1,*) 'Issue in prototype interface for phase ',p

            FLAbort('Failed to extract phase density from state in copy_into_state')

         end if

         number_nodes = node_count(density)

         den_ele_loop: do i = 1,element_count(density)

            element_nodes => ele_nodes(density,i)

            nloc = size(element_nodes)

            den_node_loop: do j = 1,nloc

               call set(density, &
                    element_nodes(j), &
                    proto_densities((cv_ndgln((i-1)*nloc+j)) + (p-1)*number_nodes))

            end do den_node_loop

         end do den_ele_loop

         velocity => extract_vector_field(state(p), "Velocity", stat=stat)

         if (stat /= 0) then 

            ewrite(1,*) 'Issue in prototype interface for phase ',p

            FLAbort('Failed to extract phase velocity from state in copy_into_state')

         end if

         number_nodes = node_count(velocity)

         vel_ele_loop: do i = 1,element_count(velocity)

            element_nodes => ele_nodes(velocity,i)

            nloc = size(element_nodes)

            vel_node_loop: do j = 1,nloc

               ! set u
               call set(velocity, &
                    1, &
                    element_nodes(j), &
                    proto_velocity_u((u_ndgln((i-1)*nloc+j)) + (p-1)*number_nodes))

               ! set v
               if (ndim > 1) call set(velocity, &
                    2, &
                    element_nodes(j), &
                    proto_velocity_v((u_ndgln((i-1)*nloc+j)) + (p-1)*number_nodes))

               ! set w
               if (ndim > 2) call set(velocity, &
                    3, &
                    element_nodes(j), &
                    proto_velocity_w((u_ndgln((i-1)*nloc+j)) + (p-1)*number_nodes))

            end do vel_node_loop

         end do vel_ele_loop

      end do phase_loop

      ! comp is stored in the order
      !   comp1 phase1
      !   comp1 phase2
      !   comp2 phase1
      !   comp2 phase2
      !   etc
      have_comp: if (ncomp > 0) then

         ! there is a state for each phase and each component
         nstates = size(state)

         ! Assume for now that components have been inserted in state AFTER all the phases
         comp_loop: do i = nstates-ncomp+1,nstates

            ! inspect each scalar field of this state   
            sfield_loop: do j = 1,option_count("/material_phase[" // int2str(i-1) //"]/scalar_field")

               ! extract each scalar field of this state
               componentmassfraction => extract_scalar_field(state(i), j)

               ! determine if this scalar field has the magic name for components   
               is_compfrac: if (componentmassfraction%name(1:21) == "ComponentMassFraction") then

                  ! find the phase this component fraction is associated with   
                  call get_option("/material_phase[" // int2str(i-1) //"]/scalar_field[" // int2str(j-1) //&
                       &"]/material_phase_name", material_phase_name)

                  ! find the phase index   
                  phase_index_loop: do k = 1,nphase

                     phase_name_check: if (trim(material_phase_name) == state(k)%name) then

                        number_nodes = node_count(componentmassfraction)

                        comp_ele_loop: do ele = 1,element_count(componentmassfraction)

                           element_nodes => ele_nodes(componentmassfraction,ele)

                           nloc = size(element_nodes)

                           comp_node_loop: do jloc= 1,nloc

                              call set(componentmassfraction, &
                                   element_nodes(jloc), &
                                   proto_components((cv_ndgln((ele-1)*nloc+jloc)) + ((i-(1+nphase))*nphase+(k-1))*number_nodes))

                           end do comp_node_loop

                        end do comp_ele_loop

                        cycle sfield_loop

                     end if phase_name_check

                  end do phase_index_loop

               end if is_compfrac

            end do sfield_loop

         end do comp_loop

      end if have_comp

      pressure => extract_scalar_field(state(1), "Pressure")    

      press_ele_loop: do i = 1,ele_count(pressure)

         element_nodes => ele_nodes(pressure,i)

         nloc = size(element_nodes)

         press_node_loop: do j = 1,nloc

            call set(pressure, &
                 element_nodes(j), &
                 proto_pressure(p_ndgln((i-1)*nloc+j)))

         end do press_node_loop

      end do press_ele_loop

      ewrite(3,*) "Leaving copy_into_state"

      return

    end subroutine copy_into_state

  end module copy_outof_into_state
