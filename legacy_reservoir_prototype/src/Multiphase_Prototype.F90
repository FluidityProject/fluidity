
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

module mp_prototype
  ! Prototype modules
  use input_var
  use multiphase_mom_press_volf
  use multiphase_field_advection
  use spact
  use printout
  
  ! Transition modules
  use copy_outof_into_state
  
  ! New modules
  use fldebug
  use state_module


    implicit none

  contains

    subroutine multiphase_prototype(state, dt, current_time, finish_time, &
                                    nonlinear_iterations, nonlinear_iteration_tolerance)

      implicit none
      
      !! New variable declaration
      
      type(state_type), dimension(:), pointer :: state
      
      integer :: nonlinear_iterations  !! equal to nits in prototype code
      
      real :: dt, current_time, finish_time
      real :: nonlinear_iteration_tolerance

      
      !! Old variable declaration

      integer :: problem, nphase, ncomp, totele, ndim, nlev, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, mat_nloc, &
           cv_snloc, u_snloc, p_snloc, x_snloc, stotel, &
           ncoef, nuabs_coefs, &
           u_ele_type, p_ele_type, mat_ele_type, cv_ele_type, &
           cv_sele_type, u_sele_type

      integer :: ntime, ntime_dump, nits, nits_internal, ndpset, noit_dim, &
           nits_flux_lim_volfra, nits_flux_lim_comp

      real :: patmos, p_ini, t_ini, t_beta, v_beta, t_theta, v_theta, &
           u_theta, domain_length, Mobility

      integer :: t_disopt, u_disopt, v_disopt, t_dg_vel_int_opt, &
           u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
           nopt_vel_upwind_coefs

      logical :: lump_eqns, volfra_use_theta_flux, volfra_get_theta_flux, &
           comp_use_theta_flux, comp_get_theta_flux, KComp_Sigmoid, Comp_Sum2One

      real :: volfra_error, volfra_relax, volfra_relax_diag, volfra_relax_row,  & 
           scalar_error, scalar_relax, scalar_relax_diag, scalar_relax_row, & 
           global_error, global_relax, global_relax_diag, global_relax_row, & 
           velocity_error, velocity_relax, velocity_relax_diag, velocity_relax_row, & 
           pressure_error, pressure_relax, pressure_relax_diag, pressure_relax_row, &
           mass_matrix_error, mass_matrix_relax, mass_matrix_relax_diag, mass_matrix_relax_row 

      integer :: volfra_relax_number_iterations, scalar_relax_number_iterations, &
           global_relax_number_iterations,  velocity_relax_number_iterations, &
           pressure_relax_number_iterations, mass_matrix_relax_number_iterations, &
           in_ele_upwind, dg_ele_upwind

      integer :: comp_diffusion_opt, ncomp_diff_coef, capil_pres_opt, ncapil_pres_coef
      real, dimension( : , : , : ), allocatable :: comp_diff_coef, capil_pres_coef
      real, dimension( : , : , : , : ), allocatable :: comp_diffusion

      ! Derived data
      integer :: cv_nonods, u_nonods, p_nonods, cv_pha_nonods, u_pha_nonods, mat_nonods, &
           x_nonods, xu_nonods, nlenmcy,ncp_coefs
      integer, dimension( : ), allocatable :: wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc, &
           wic_comp_bc
      real, dimension( : ), allocatable :: suf_vol_bc, suf_d_bc, suf_cpd_bc, suf_t_bc, suf_p_bc, &
           suf_u_bc, suf_v_bc, suf_w_bc, suf_one_bc, suf_comp_bc, &
           suf_u_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob1, suf_v_bc_rob2, &
           suf_w_bc_rob1, suf_w_bc_rob2, suf_t_bc_rob1, suf_t_bc_rob2, &
           suf_vol_bc_rob1, suf_vol_bc_rob2, &
           suf_comp_bc_rob1, suf_comp_bc_rob2, &
           opt_vel_upwind_coefs

      real, dimension( : ), allocatable :: sat_error_relax2_noit, t_error_relax2_noit, gl_error_relax2_noit, &
           u_error_relax2_noit, p_error_relax2_noit, mass_error_relax2_noit

      real, dimension( : ), allocatable :: x, y, z, xu, yu, zu, nu, nv, nw, ug, vg, wg

      integer, dimension( : ), allocatable :: uabs_option
      real, dimension( :, : ), allocatable :: uabs_coefs
      real, dimension( :, :, : ), allocatable :: u_abs_stab, u_absorb, comp_absorb
      real, dimension( : ), allocatable :: u_source, t_source, v_source, comp_source
      real, dimension( :, :, : ), allocatable :: t_absorb, v_absorb

      real, dimension( :, :, :, : ), allocatable :: udiffusion, tdiffusion
      real, dimension( : ), allocatable :: u, v, w

      real, dimension( : ), allocatable :: den, satura, volfra, comp, t, p, cv_p, volfra_pore, &
           cv_one, Viscosity
      real, dimension( :, :, : ), allocatable :: perm

      integer, dimension( : ), allocatable :: eos_option, cp_option
      real, dimension( :, : ), allocatable :: eos_coefs, cp_coefs
      real, dimension( : , : , : ), allocatable :: K_Comp
      real :: alpha_beta

      ! Working vectors
      real, dimension( : ), allocatable :: uold, vold, wold, nuold, nvold, nwold, &
           denold, uden, udenold, deriv, saturaold, compold, volold, told, cpden, cpdenold, pold, cv_pold
      integer, dimension( : ), allocatable :: u_ndgln, xu_ndgln, cv_ndgln, x_ndgln, p_ndgln, &
           mat_ndgln, u_sndgln, cv_sndgln, x_sndgln, p_sndgln

      ! For matrix sparsity
      integer :: ncoldgm_pha
      integer, dimension( : ), allocatable :: finmcy, midmcy, midacv, finacv, finele, midele, &
           centc, findc, midcmc, findcmc, centct, findct, findgm_pha, middgm_pha, findm, midm
      integer, dimension( : ), allocatable :: colacv, colele, colct, colc, coldgm_pha, colmcy, colcmc, colm
      integer :: mxnele, mx_nct, mx_nc, mx_ncolcmc, mx_ncoldgm_pha, mx_ncolmcy, mx_ncolacv, mx_ncolm
      integer :: ncolacv, ncolmcy, ncolele, ncolct, ncolc, ncolcmc, ncolm
      integer :: nkcomp, iphase, icomp

      integer :: option_debug
      integer, parameter :: unit_input = 5, unit_debug = 101


      open( unit_input, file = 'input.dat', status = 'unknown' )
      open( unit_debug, file = 'mirror_int_data.dat', status = 'unknown' )
      !open( 357, file = 'flog.dat', status = 'unknown')
      !open( 357, file = '/dev/null', status = 'unknown')

      ewrite(3,*) 'In multiphase_prototype'

!!!!!!!!
!!!!!   Major insert required here to pull everything required out of state
!!!!!    - for this we need to list everything that's read in and not derived
!!!!!    - although some derived stuff might be available through state as well
!!!!!!!!
      call copy_outof_state(state, dt, current_time, finish_time, &
                            nonlinear_iterations, nonlinear_iteration_tolerance)

      call read_scalar( unit_input, option_debug, problem, nphase, ncomp, totele, ndim, nlev, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, &
           cv_snloc, u_snloc, p_snloc, x_snloc, stotel, &
           ncoef, nuabs_coefs, &
           u_ele_type, p_ele_type, mat_ele_type, cv_ele_type, &
           cv_sele_type, u_sele_type, ntime, ntime_dump, nits, nits_internal, noit_dim, &
           nits_flux_lim_volfra, nits_flux_lim_comp, &
           ndpset, &
           dt, patmos, p_ini, t_ini, &
           t_beta, v_beta, t_theta, v_theta, u_theta, &
           t_disopt, u_disopt, v_disopt, t_dg_vel_int_opt, &
           u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
           lump_eqns, & 
           volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux, &
           capil_pres_opt, ncapil_pres_coef, comp_diffusion_opt, ncomp_diff_coef, &
           domain_length )

      if( option_debug == 357 )then
         open( 357, file = 'flog.dat', status = 'unknown')
      else
         open( 357, file = '/dev/null', status = 'unknown')
      endif

      if( u_snloc < 0 ) u_snloc = 1 * nlev
      mat_nloc = cv_nloc

      select case( problem )
      case( -1 ); ! CV-Adv DG
         cv_nonods = cv_nloc * totele
      case( 2 ); ! BL DG
         cv_nonods = cv_nloc * totele
      case default ;
         cv_nonods = ( cv_nloc - 1 ) * totele + 1
      end select

      u_nonods = u_nloc * totele
      p_nonods = cv_nonods
      mat_nonods = mat_nloc * totele
      cv_pha_nonods = cv_nonods * nphase
      u_pha_nonods = u_nonods * nphase
      ncp_coefs = nphase
      x_nonods = max(( x_nloc - 1 ) * totele + 1, totele )
      xu_nonods = max(( xu_nloc - 1 ) * totele + 1, totele )
      nlenmcy = u_pha_nonods + cv_nonods
      nopt_vel_upwind_coefs = mat_nonods * nphase * ndim * ndim * 2

      if( ndpset < 0 ) ndpset = cv_nonods

      ! Variables in which the dimensions depend upon input data
      allocate( wic_vol_bc( stotel * nphase ))
      allocate( wic_d_bc( stotel * nphase ))
      allocate( wic_u_bc( stotel * nphase ))
      allocate( wic_p_bc( stotel * nphase ))
      allocate( wic_t_bc( stotel * nphase ))
      allocate( wic_comp_bc( stotel * nphase ))

      allocate( suf_vol_bc( stotel * cv_snloc * nphase ))
      allocate( suf_d_bc( stotel * cv_snloc * nphase ))
      allocate( suf_cpd_bc( stotel * cv_snloc * nphase ))
      allocate( suf_t_bc( stotel * cv_snloc * nphase ))

      allocate( suf_p_bc( stotel * p_snloc * nphase ))
      allocate( suf_u_bc( stotel * u_snloc * nphase ))
      allocate( suf_v_bc( stotel * u_snloc * nphase ))
      allocate( suf_w_bc( stotel * u_snloc * nphase ))
      allocate( suf_one_bc( stotel * cv_snloc * nphase ))
      allocate( suf_comp_bc( stotel * cv_snloc * nphase * ncomp ))

      allocate( suf_u_bc_rob1( stotel * u_snloc * nphase ))
      allocate( suf_v_bc_rob1( stotel * u_snloc * nphase ))
      allocate( suf_w_bc_rob1( stotel * u_snloc * nphase ))
      allocate( suf_u_bc_rob2( stotel * u_snloc * nphase ))
      allocate( suf_v_bc_rob2( stotel * u_snloc * nphase ))
      allocate( suf_w_bc_rob2( stotel * u_snloc * nphase ))
      allocate( suf_t_bc_rob1( stotel * cv_snloc * nphase ))
      allocate( suf_t_bc_rob2( stotel * cv_snloc * nphase ))
      allocate( suf_vol_bc_rob1( stotel * cv_snloc * nphase ))
      allocate( suf_vol_bc_rob2( stotel * cv_snloc * nphase ))
      allocate( suf_comp_bc_rob1( stotel * cv_snloc * nphase ))
      allocate( suf_comp_bc_rob2( stotel * cv_snloc * nphase ))

      allocate( opt_vel_upwind_coefs( nopt_vel_upwind_coefs ))

      allocate( x( x_nonods ))
      allocate( y( x_nonods ))
      allocate( z( x_nonods ))
      allocate( xu( xu_nonods ))
      allocate( yu( xu_nonods ))
      allocate( zu( xu_nonods ))
      allocate( nu( u_nonods * nphase ))
      allocate( nv( u_nonods * nphase ))
      allocate( nw( u_nonods * nphase ))
      allocate( ug( u_nonods * nphase ))
      allocate( vg( u_nonods * nphase ))
      allocate( wg( u_nonods * nphase ))

      allocate( uabs_option( nphase ))
      allocate( uabs_coefs( nphase, nuabs_coefs ))
      allocate( u_abs_stab( mat_nonods, ndim * nphase, ndim * nphase ))
      allocate( u_absorb( mat_nonods, ndim * nphase, ndim * nphase ))
      allocate( t_absorb( cv_nonods, nphase, nphase ))
      allocate( v_absorb( cv_nonods, nphase, nphase ))
      allocate( comp_absorb( cv_nonods, nphase, nphase ))
      allocate( u_source( u_pha_nonods ))
      allocate( t_source( cv_pha_nonods ))
      allocate( v_source( cv_pha_nonods ))
      allocate( comp_source( cv_pha_nonods ))

      allocate( udiffusion( mat_nonods, ndim, ndim, nphase ))
      allocate( tdiffusion( mat_nonods, ndim, ndim, nphase ))

      ! Now initialising remaining prognostic and diagnostic velocity fields
      allocate( u( u_pha_nonods ))
      allocate( v( u_pha_nonods ))
      allocate( w( u_pha_nonods ))

      ! Now initialising remaining prognostic and diagnostic scalar fields
      allocate( den( cv_pha_nonods ))
      allocate( satura( cv_pha_nonods ))
      allocate( comp( cv_pha_nonods * ncomp ))
      allocate( volfra( cv_pha_nonods ))
      allocate( t( cv_pha_nonods ))
      allocate( cv_one( cv_pha_nonods ))
      allocate( Viscosity( cv_pha_nonods ))
      allocate( p( cv_nonods ))
      allocate( cv_p( cv_nonods ))

      allocate( volfra_pore( totele ))

      allocate( perm( totele, ndim, ndim ))

      ! Input of EOS and Absorption coefficients
      allocate( eos_option( nphase ))
      allocate( cp_option( nphase ))
      allocate( eos_coefs( nphase, ncoef ))
      allocate( cp_coefs( nphase, ncp_coefs ))

      allocate( capil_pres_coef( ncapil_pres_coef, nphase, nphase ))
      allocate( comp_diffusion( mat_nonods, ndim, ndim, nphase ))
      allocate( comp_diff_coef( ncomp, ncomp_diff_coef, nphase ))

      nkcomp = Combination( ncomp, 2 )
      allocate( K_Comp( ncomp, nphase, nphase ))

      call read_all( unit_input, nphase, ncomp, totele, ndim, &
           cv_snloc, u_snloc, p_snloc, stotel, &
           ncoef, nuabs_coefs, ncp_coefs, & 
           cv_nonods, u_nonods, &
           mat_nonods, x_nonods, xu_nonods, &
           nopt_vel_upwind_coefs, &
           wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc, wic_comp_bc, & 
           suf_vol_bc, suf_d_bc, suf_cpd_bc, suf_t_bc, suf_p_bc, &
           suf_u_bc, suf_v_bc, suf_w_bc, suf_one_bc, suf_comp_bc, &
           suf_u_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob1, suf_v_bc_rob2, &
           suf_w_bc_rob1, suf_w_bc_rob2, suf_t_bc_rob1, suf_t_bc_rob2, &
           suf_vol_bc_rob1, suf_vol_bc_rob2, &
           suf_comp_bc_rob1, suf_comp_bc_rob2, &
           opt_vel_upwind_coefs, &
           volfra_error, volfra_relax, volfra_relax_diag, volfra_relax_row, volfra_relax_number_iterations, & 
           scalar_error, scalar_relax, scalar_relax_diag, scalar_relax_row, scalar_relax_number_iterations, & 
           global_error, global_relax, global_relax_diag, global_relax_row, global_relax_number_iterations, & 
           velocity_error, velocity_relax, velocity_relax_diag, velocity_relax_row, velocity_relax_number_iterations, & 
           pressure_error, pressure_relax, pressure_relax_diag, pressure_relax_row, pressure_relax_number_iterations, & 
           mass_matrix_error, mass_matrix_relax, mass_matrix_relax_diag, mass_matrix_relax_row, mass_matrix_relax_number_iterations, &
           in_ele_upwind, dg_ele_upwind, &
           x, y, z, xu, yu, zu, nu, nv, nw, ug, vg, wg, &
           uabs_option, uabs_coefs, u_abs_stab, Mobility, &
           u_absorb, t_absorb, v_absorb, comp_absorb, &
           u_source, t_source, v_source, comp_source, udiffusion, tdiffusion, &
           ncomp_diff_coef, comp_diffusion, comp_diff_coef, &
           ncapil_pres_coef, capil_pres_coef, & 
           u, v, w, &
           den, satura, comp, volfra, t, cv_one, p, cv_p, volfra_pore, Viscosity, perm, &
           KComp_Sigmoid, K_Comp, Comp_Sum2One, alpha_beta, &
           eos_option, cp_option, eos_coefs, cp_coefs )


      close( unit_input )


      ! Allocating solver parameters into the following arrays (dimension = 4)
      allocate( sat_error_relax2_noit( noit_dim ))
      allocate( t_error_relax2_noit( noit_dim ))
      allocate( gl_error_relax2_noit( noit_dim ))
      allocate( u_error_relax2_noit( noit_dim ))
      allocate( p_error_relax2_noit( noit_dim ))
      allocate( mass_error_relax2_noit( noit_dim ))

      sat_error_relax2_noit( 1 ) = volfra_error
      sat_error_relax2_noit( 2 ) = volfra_relax
      sat_error_relax2_noit( 3 ) = volfra_relax_diag
      sat_error_relax2_noit( 4 ) = volfra_relax_row
      sat_error_relax2_noit( 5 ) = real( volfra_relax_number_iterations )

      t_error_relax2_noit( 1 ) = scalar_error
      t_error_relax2_noit( 2 ) = scalar_relax
      t_error_relax2_noit( 3 ) = scalar_relax_diag
      t_error_relax2_noit( 4 ) = scalar_relax_row
      t_error_relax2_noit( 5 ) = real( scalar_relax_number_iterations )

      gl_error_relax2_noit( 1 ) = global_error
      gl_error_relax2_noit( 2 ) = global_relax
      gl_error_relax2_noit( 3 ) = global_relax_diag
      gl_error_relax2_noit( 4 ) = global_relax_row
      gl_error_relax2_noit( 5 ) = real( global_relax_number_iterations )

      u_error_relax2_noit( 1 ) = velocity_error
      u_error_relax2_noit( 2 ) = velocity_relax
      u_error_relax2_noit( 3 ) = velocity_relax_diag
      u_error_relax2_noit( 4 ) = velocity_relax_row
      u_error_relax2_noit( 5 ) = real( velocity_relax_number_iterations )

      p_error_relax2_noit( 1 ) = pressure_error
      p_error_relax2_noit( 2 ) = pressure_relax
      p_error_relax2_noit( 3 ) = pressure_relax_diag
      p_error_relax2_noit( 4 ) = pressure_relax_row
      p_error_relax2_noit( 5 ) = real( pressure_relax_number_iterations )

      mass_error_relax2_noit( 1 ) = mass_matrix_error 
      mass_error_relax2_noit( 2 ) = mass_matrix_relax
      mass_error_relax2_noit( 3 ) = mass_matrix_relax_diag 
      mass_error_relax2_noit( 4 ) = mass_matrix_relax_row 
      mass_error_relax2_noit( 5 ) = real( mass_matrix_relax_number_iterations )

      ! Set up scalar and vector fields for previous time-step
      allocate( uold( u_pha_nonods ))
      allocate( vold( u_pha_nonods ))
      allocate( wold( u_pha_nonods ))
      allocate( nuold( u_pha_nonods ))
      allocate( nvold( u_pha_nonods ))
      allocate( nwold( u_pha_nonods ))
      allocate( denold( cv_pha_nonods ))
      allocate( uden( cv_pha_nonods ))
      allocate( udenold( cv_pha_nonods ))
      allocate( deriv( cv_pha_nonods ))
      allocate( saturaold( cv_pha_nonods ))
      allocate( compold( cv_pha_nonods * ncomp ))
      allocate( volold( cv_pha_nonods ))
      allocate( told( cv_pha_nonods ))
      allocate( cpden( cv_pha_nonods ))
      allocate( cpdenold( cv_pha_nonods ))
      allocate( pold( cv_nonods ))
      allocate( cv_pold( cv_nonods ))

      uold = 0.
      vold = 0.
      wold = 0.
      denold = 1.
      uden = 0.
      udenold = 0.
      deriv = 0.
      told = 0.
      cpden = 0.
      pold = 0.
      nuold = 0.
      nvold = 0.
      nwold = 0.

      ! Set up Global node number for velocity and scalar fields
      allocate( u_ndgln( totele * u_nloc ))
      allocate( xu_ndgln( totele * xu_nloc ))
      allocate( cv_ndgln( totele * cv_nloc ))
      allocate( x_ndgln( totele * cv_nloc ))
      allocate( p_ndgln( totele * p_nloc ))
      allocate( mat_ndgln( totele * mat_nloc ))
      allocate( u_sndgln( stotel * u_snloc ))
      allocate( cv_sndgln( stotel * cv_snloc ))
      allocate( x_sndgln( stotel * cv_snloc ))
      allocate( p_sndgln( stotel * p_snloc ))

      u_ndgln = 0
      xu_ndgln = 0
      cv_ndgln = 0
      x_ndgln = 0
      p_ndgln = 0
      mat_ndgln =0
      u_sndgln = 0
      cv_sndgln =0
      x_sndgln = 0
      p_sndgln = 0

      call allocating_global_nodes( ndim, totele, domain_length, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, mat_nloc, &
           cv_snloc, u_snloc, p_snloc, stotel, &
           cv_nonods, u_nonods, x_nonods, xu_nonods, &
           u_ele_type, cv_ele_type, &
           x, xu,  &
           u_ndgln, xu_ndgln, cv_ndgln, x_ndgln, p_ndgln, &
           mat_ndgln, u_sndgln, cv_sndgln, p_sndgln )

      ! Initialising T and Told. This should be properly initialised through a generic function
      call initialise_scalar_fields( &
           problem, ndim, nphase, totele, domain_length, &
           x_nloc, cv_nloc, x_nonods, cv_nonods,  &
           x_ndgln, cv_ndgln, &
           x, told, t )

      close( unit_input )

      ! Mirroring Input dat

      if( .true. ) call mirror_data( unit_debug, problem, nphase, ncomp, totele, ndim, nlev, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, &
           cv_snloc,  p_snloc, stotel, &
           ncoef, nuabs_coefs, &
           u_ele_type, p_ele_type, mat_ele_type, cv_ele_type, &
           cv_sele_type, u_sele_type, &
           ntime, nits, ndpset, &
           dt, patmos, p_ini, t_ini, &
           t_beta, v_beta, t_theta, v_theta, u_theta, &
           t_disopt, u_disopt, v_disopt, t_dg_vel_int_opt, &
           u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
           domain_length, u_snloc, mat_nloc, cv_nonods, u_nonods, &
           p_nonods, mat_nonods, ncp_coefs, x_nonods, xu_nonods, &
           nlenmcy, &
           nopt_vel_upwind_coefs, &
           u_ndgln, xu_ndgln, cv_ndgln, x_ndgln, p_ndgln, &
           mat_ndgln, u_sndgln, cv_sndgln, x_sndgln, p_sndgln, &
           wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc, & 
           suf_vol_bc, suf_d_bc, suf_cpd_bc, suf_t_bc, suf_p_bc, &
           suf_u_bc, &
           suf_u_bc_rob1, suf_u_bc_rob2, &
           opt_vel_upwind_coefs, &
           x, xu, nu, ug, &
           uabs_option, u_abs_stab, u_absorb, &
           u_source, &
           u,  &
           den, satura, comp, p, cv_p, volfra_pore, perm )


      ! Sparsity patterns
      allocate( finmcy( nlenmcy + 1 )) ! Force balance plus cty multi-phase eqns
      allocate( midmcy( nlenmcy ))
      allocate( midacv( cv_pha_nonods )) ! CV multi-phase eqns 
      allocate( finacv( cv_pha_nonods + 1 ))
      allocate( finele( totele + 1 )) ! Element connectivity
      allocate( midele( totele ))
      allocate( centc( u_nonods )) ! C sparsity operating on pressure in force balance
      allocate( findc( u_nonods + 1 ))
      allocate( midcmc( cv_nonods )) ! pressure matrix for projection method
      allocate( findcmc( cv_nonods + 1 ))
      allocate( centct( cv_nonods )) ! CT sparsity - global cty eqn.
      allocate( findct( cv_nonods + 1 ))
      allocate( findgm_pha( u_pha_nonods + 1 )) ! Force balance sparsity
      allocate( middgm_pha( u_pha_nonods ))
      allocate( findm( cv_nonods + 1 )) ! Sparsity for the CV-FEM
      allocate( midm( cv_nonods ))


      ! Defining lengths
      mxnele = ( 2 * ndim + 1 ) * totele
      mx_nct = cv_nonods * ( 2 * u_nloc + 1 ) * ndim * nphase
      mx_nc = mx_nct  

      ! select case( problem )
      ! case( -2 ); ! CV-Adv (Cty)
      !    mx_ncolcmc = ( 2 * cv_nloc + 1 ) * cv_nonods
      ! case( -1 ); ! CV-Adv (DG)
      !    mx_ncolcmc = ( 2 * cv_nloc + 1 ) * cv_nonods
      ! case( 0 ); ! CV-Adv (Std)
      !    mx_ncolcmc = ( 2 * cv_nloc + 3 ) * cv_nonods
      ! case( 1, 2); ! BL-test1
      mx_ncolcmc = ( 2 * ( cv_nloc + 2 ) + 1 ) * cv_nonods
      ! end select

      mx_ncoldgm_pha = mxnele * ( u_nloc * ndim )**2 * nphase + totele * ( u_nloc * ndim * nphase )**2
      mx_ncolmcy = mx_ncoldgm_pha + mx_nct + mx_nc + mx_ncolcmc
      mx_ncolacv = ( 2 * ndim + 1 ) * cv_nonods * nphase + cv_nonods * ( nphase - 1 ) * nphase
      mx_ncolm = mxnele * cv_nloc * cv_nloc

      allocate( colmcy( mx_ncolmcy ))
      allocate( colacv( mx_ncolacv ))
      allocate( colele( mxnele ))
      allocate( colct( mx_nct ))
      allocate( colc( mx_nc ))
      allocate( coldgm_pha( mx_ncoldgm_pha ))
      allocate( colcmc( mx_ncolcmc ))
      allocate( colm( mx_ncolm ))

      call get_spars_pats( &
           ndim, u_nonods * nphase, cv_nonods * nphase, &
           u_nonods, cv_nonods, &
           u_nloc, cv_nloc, nphase, totele, u_ndgln, &
           mx_ncolacv, ncolacv, finacv, colacv, midacv, & ! CV multi-phase eqns (e.g. vol frac, temp)
           nlenmcy, mx_ncolmcy, ncolmcy, finmcy, colmcy, midmcy, & ! Force balance plus cty multi-phase eqns
           mxnele, ncolele, midele, finele, colele, & ! Element connectivity 
           mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, & ! Force balance sparsity  
           mx_nct, ncolct, findct, colct, & ! CT sparsity - global cty eqn
           mx_nc, ncolc, findc, colc, & ! C sparsity operating on pressure in force balance
           mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, & ! pressure matrix for projection method
           mx_ncolm, ncolm, findm, colm, midm, u_ele_type )

      close( unit_debug )

      if( .true. ) call check_sparsity( &
           u_nonods * nphase, cv_nonods * nphase, &
           u_nonods, cv_nonods, totele, &
           mx_ncolacv, ncolacv, finacv, colacv, midacv, & ! CV multi-phase eqns (e.g. vol frac, temp)
           nlenmcy, mx_ncolmcy, ncolmcy, finmcy, colmcy, midmcy, & ! Force balance plus cty multi-phase eqns
           mxnele, ncolele, midele, finele, colele, & ! Element connectivity 
           mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, & ! Force balance sparsity  
           mx_nct, ncolct, findct, colct, & ! CT sparsity - global cty eqn
           mx_nc, ncolc, findc, colc, & ! C sparsity operating on pressure in force balance
           mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, & ! pressure matrix for projection method
           mx_ncolm, ncolm, findm, colm, midm ) ! CV-FEM matrix

      Select Case( problem )

      Case( : 0 ) ; ! CV-Adv test case: -2( Cty ), -1( DG ), 0( Std )

         call solve_multiphase_field_advection( problem, nphase, ncomp, totele, ndim, &
              u_nloc, xu_nloc, cv_nloc, x_nloc, mat_nloc, &
              cv_snloc, u_snloc, stotel, &
              domain_length, &
                                ! Element types
              u_ele_type, p_ele_type, cv_ele_type, &
              cv_sele_type, u_sele_type, &
                                ! Total time loop and initialisation parameters
              ntime, ntime_dump, nits, &
              nits_flux_lim_volfra, nits_flux_lim_comp, &
              dt,  &
                                ! Discretisation parameters
              t_beta, t_theta, t_disopt, t_dg_vel_int_opt, lump_eqns, &
              volfra_use_theta_flux, volfra_get_theta_flux, &
              opt_vel_upwind_coefs, nopt_vel_upwind_coefs, &
              noit_dim, &
              sat_error_relax2_noit, t_error_relax2_noit, gl_error_relax2_noit, &
              u_error_relax2_noit, p_error_relax2_noit, mass_error_relax2_noit, &
              in_ele_upwind, dg_ele_upwind, &
                                ! Total nodes for different meshes
              cv_nonods, u_nonods, mat_nonods, x_nonods, &
              u_ndgln, xu_ndgln, cv_ndgln, x_ndgln, &
              mat_ndgln, u_sndgln, cv_sndgln, &
                                ! Boundary conditions and surface elements
              wic_d_bc, wic_u_bc, wic_t_bc, &
              suf_d_bc, suf_t_bc, suf_u_bc, suf_v_bc, suf_w_bc, &
              suf_t_bc_rob1, suf_t_bc_rob2, &
                                ! Positions and grid velocities
              x, y, z, nu, nv, nw, ug, vg, wg, &
                                ! Absorption and source terms and coefficients
              t_absorb, t_source, &
                                ! Diffusion parameters
              tdiffusion, &
                                ! Scalar fields
              t, p, cv_one, volfra_pore, deriv, &
              told, nuold, nvold, nwold, &
                                ! Matrices sparsity
              mx_ncolacv, ncolacv, finacv, colacv, midacv, & ! CV multi-phase eqns (e.g. vol frac, temp)
              mx_nct, ncolct, findct, colct, & ! CT sparsity - global cty eqn
              ncolm, findm, colm, midm, & ! CV-FEM matrix
              mxnele, ncolele, finele, colele ) ! Element connectivity 

      Case( 1,2 ); 

         call solve_multiphase_mom_press_volf( nphase, ncomp, totele, ndim, &
                                ! Nodes et misc
              u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, mat_nloc, &
              cv_snloc, u_snloc, p_snloc, stotel, &
              ncoef, nuabs_coefs, &
                                ! Element types
              u_ele_type, p_ele_type, cv_ele_type, &
              cv_sele_type, u_sele_type, &
                                ! Total time loop and initialisation parameters
              ntime, ntime_dump, nits, nits_internal, &
              nits_flux_lim_volfra, nits_flux_lim_comp, & 
              ndpset, &
              dt, &
                                ! Discretisation parameters
              v_beta, v_theta, &
              v_disopt, &
              v_dg_vel_int_opt, &
              t_beta, t_theta, t_disopt, t_dg_vel_int_opt, lump_eqns, &
              volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux, &
              opt_vel_upwind_coefs, nopt_vel_upwind_coefs, &
              noit_dim, &
              sat_error_relax2_noit, t_error_relax2_noit, gl_error_relax2_noit, &
              u_error_relax2_noit, p_error_relax2_noit, mass_error_relax2_noit, &
              in_ele_upwind, dg_ele_upwind, &
                                ! Total nodes for different meshes
              domain_length, &
              cv_nonods, u_nonods, cv_pha_nonods, u_pha_nonods, mat_nonods, &
              x_nonods, xu_nonods, &
              u_ndgln, xu_ndgln, cv_ndgln, x_ndgln, p_ndgln, mat_ndgln, &
              u_sndgln, cv_sndgln, p_sndgln, &
                                ! Boundary conditions and surface elements
              wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc, wic_comp_bc, &
              suf_vol_bc, suf_d_bc, suf_p_bc, suf_t_bc, &
              suf_u_bc, suf_v_bc, suf_w_bc, suf_comp_bc, &
              suf_u_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob1, suf_v_bc_rob2, &
              suf_w_bc_rob1, suf_w_bc_rob2, suf_comp_bc_rob1, suf_comp_bc_rob2, &
              suf_t_bc_rob1, suf_t_bc_rob2, suf_vol_bc_rob1, suf_vol_bc_rob2, &
                                ! Positions and grid velocities
              x, y, z, nu, nv, nw, ug, vg, wg, &
                                ! Absorption and source terms and coefficients
              uabs_option, uabs_coefs, u_abs_stab, Mobility, &
              u_absorb, v_absorb, comp_absorb, &
              u_source, v_source, comp_source, &
              t_absorb, t_source, Comp_Sum2One, &
                                ! Diffusion parameters
              udiffusion, &
              tdiffusion, &
              comp_diffusion_opt, ncomp_diff_coef, comp_diffusion, comp_diff_coef, &
                                ! Velocities and scalar fields
              u, v, w, &
              den, satura, comp, t, p, cv_p, cv_one, volfra_pore, Viscosity, &
              perm, &
              uold, vold, wold, denold, saturaold, compold, uden, udenold, deriv, &
              told, pold, cv_pold, nuold, nvold, nwold, &
                                ! EOS terms
              eos_option, &
              eos_coefs, &
              KComp_Sigmoid, K_Comp, alpha_beta, &
                                ! Capillarity pressure terms
              capil_pres_opt, ncapil_pres_coef, capil_pres_coef, &
                                ! Matrices sparsity
              mx_ncolacv, ncolacv, finacv, colacv, midacv, & ! CV multi-phase eqns (e.g. vol frac, temp)
              nlenmcy, mx_ncolmcy, ncolmcy, finmcy, colmcy, midmcy, & ! Force balance plus cty multi-phase eqns
              mxnele, ncolele, finele, colele, & ! Element connectivity 
              mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, & ! Force balance sparsity  
              mx_nct, ncolct, findct, colct, & ! CT sparsity - global cty eqn
              mx_nc, ncolc, findc, colc, & ! C sparsity operating on pressure in force balance
              mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, & ! pressure matrix for projection method
              mx_ncolm, ncolm, findm, colm, midm ) ! CV-FEM matrix

      end Select
      
      call copy_into_state(state, satura, p)

      ewrite(3,*) 'Leaving multiphase_prototype'
      close( 357 )

      return
    end subroutine multiphase_prototype

  end module mp_prototype

