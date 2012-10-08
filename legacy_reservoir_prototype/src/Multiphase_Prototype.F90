
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
  use multiphase_mom_press_volf
  use spact
  use printout
  
  ! Transition modules
  use copy_outof_into_state
  
  ! New modules
  use fldebug
  use state_module


    implicit none

  contains

    subroutine multiphase_prototype(state, dt, &
         nonlinear_iterations, nonlinear_iteration_tolerance, &
         dump_no)

      !! New variable declaration

      type(state_type), dimension(:), pointer :: state

      integer :: dump_no, nonlinear_iterations  !! equal to nits in prototype code

      real :: dt 
      real :: nonlinear_iteration_tolerance


      !! Old variable declaration

      integer :: nphase, ncomp, totele, ndim, nlev, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, mat_nloc, &
           cv_snloc, u_snloc, p_snloc, x_snloc, stotel, &
           ncoef, nuabs_coefs, &
           u_ele_type, p_ele_type, mat_ele_type, cv_ele_type, &
           cv_sele_type, u_sele_type

      integer :: ntime_dump, nits, nits_internal, ndpset, noit_dim, &
           nits_flux_lim_volfra, nits_flux_lim_comp, nits_flux_lim_t

      real :: patmos, p_ini, t_beta, v_beta, t_theta, v_theta, &
           u_theta, domain_length, Mobility

      integer :: t_disopt, u_disopt, v_disopt, t_dg_vel_int_opt, &
           u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
           nopt_vel_upwind_coefs

      logical :: lump_eqns, volfra_use_theta_flux, volfra_get_theta_flux, &
           comp_use_theta_flux, comp_get_theta_flux, KComp_Sigmoid, Comp_Sum2One, &
           t_use_theta_flux, t_get_theta_flux

      integer :: in_ele_upwind, dg_ele_upwind

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
      integer :: nkcomp, mx_nface_p1

      integer, parameter :: new_unit_debug = 304

      logical :: have_temperature_fields, scale_momentum_by_volume_fraction = .false.

      open( new_unit_debug, file = 'mirror_new.dat', status = 'unknown' )

      ewrite(3,*) 'In multiphase_prototype'

      call copy_outof_state(state, &
           nonlinear_iterations, nonlinear_iteration_tolerance, &
                                ! Begin here all the variables from read_scalar
           nphase, ncomp, totele, ndim, nlev, &
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
           uabs_coefs, &
           eos_coefs, cp_coefs, scale_momentum_by_volume_fraction,  &
           comp_diff_coef, capil_pres_coef, &
           u_abs_stab, u_absorb, comp_absorb, &
           t_absorb, v_absorb, &
           perm, K_Comp, &
           comp_diffusion, &
                                ! Now adding other things which we have taken inside this routine to define
           cv_nonods, p_nonods, u_nonods, x_nonods, xu_nonods, mat_nonods, &
           have_temperature_fields, nits_flux_lim_t)

      ! Test ground for sorting out memory problems

      ! Going to move to here a load of random things
      nopt_vel_upwind_coefs = mat_nonods * nphase * ndim * ndim * 2

      ! This should really be in the copy routine, but it isn't used
      ! anyway
      allocate( opt_vel_upwind_coefs( nopt_vel_upwind_coefs ))
      opt_vel_upwind_coefs = 0.

      ewrite(3,*) 'cv_ele_type', cv_ele_type
      ewrite(3,*) 'cv_sele_type', cv_sele_type
      ewrite(3,*) 'u_sele_type', u_sele_type

      ! Variables in which the dimensions depend upon input data
      allocate( udiffusion( mat_nonods, ndim, ndim, nphase ))
      allocate( tdiffusion( mat_nonods, ndim, ndim, nphase ))
      ! These can later be added into the schema, but as we are not solving for diffusion now,
      ! these can be done latter. For now, we are just allocating memory and initialise them
      udiffusion = 0.
      tdiffusion = 0.

      nkcomp = Combination( ncomp, 2 )

      ! This should really be in the copy routine, but it isn't used
      ! anyway
      !      allocate( opt_vel_upwind_coefs( nopt_vel_upwind_coefs ))
      opt_vel_upwind_coefs = 0.
      nopt_vel_upwind_coefs = mat_nonods * nphase * ndim * ndim * 2

      if( u_snloc < 0 ) u_snloc = 1 * nlev
      cv_pha_nonods = cv_nonods * nphase
      u_pha_nonods = u_nonods * nphase * ndim
      ncp_coefs = nphase
      nlenmcy = u_pha_nonods + cv_nonods

      ! Set up scalar and vector fields for previous time-step
      allocate( uold( u_nonods * nphase  ))
      allocate( vold( u_nonods * nphase ))
      allocate( wold( u_nonods * nphase ))
      allocate( nuold( u_nonods * nphase ))
      allocate( nvold( u_nonods * nphase ))
      allocate( nwold( u_nonods * nphase ))
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
      nuold = 0.
      nvold = 0.
      nwold = 0.
      denold = 1.
      uden = 0.
      udenold = 0.
      deriv = 0.
      saturaold = 0.
      compold = 0.
      volold = 0.
      told = 0.
      cpden = 0.
      cpdenold = 0.
      pold = 0.
      cv_pold = 0.
      told = t

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
      mx_nface_p1 = 2 * ndim + 1
      mxnele = mx_nface_p1 * totele
      if(cv_nonods==cv_nloc*totele) then ! is DG for pressure...
         mx_nct = totele * u_nloc * cv_nloc * mx_nface_p1 * ndim * nphase * totele
      else
         mx_nct = totele * u_nloc * cv_nloc * ndim * nphase
      endif
      mx_nc = mx_nct  

      ! Assume the DG representation requires the most stroage...
      mx_ncolcmc = mx_nface_p1**3 * cv_nloc * cv_nloc * totele  

      mx_ncoldgm_pha = ( mxnele + totele ) * ( u_nloc * ndim * nphase)**2
      mx_ncolmcy = mx_ncoldgm_pha + mx_nct + mx_nc + mx_ncolcmc
      mx_ncolacv = 3 * mx_nface_p1 * cv_nonods * nphase + cv_nonods * ( nphase - 1 ) * nphase

      if(cv_nonods==cv_nloc*totele) then ! is DG for pressure
         mx_ncolm = mxnele * cv_nloc * cv_nloc * nphase * totele
      else
         mx_ncolm = mxnele * cv_nloc * cv_nloc * nphase
      endif

      allocate( colmcy( mx_ncolmcy ) )
      allocate( colacv( mx_ncolacv ) )
      allocate( colele( mxnele ) )
      allocate( colct( mx_nct ) )
      allocate( colc( mx_nc ) )
      allocate( coldgm_pha( mx_ncoldgm_pha ) )
      allocate( colcmc( mx_ncolcmc ) )
      allocate( colm( mx_ncolm ) )

      call get_spars_pats( &
           ndim, nphase, totele, u_nonods * nphase * ndim, cv_nonods * nphase, &
           u_nonods, cv_nonods, x_nonods, &
           cv_ele_type, u_ele_type, &
           u_nloc, cv_nloc, x_nloc, xu_nloc, mat_nloc, &
           u_snloc, cv_snloc, x_snloc, &
           u_ndgln, cv_ndgln, x_ndgln, xu_ndgln, &
                                ! CV multi-phase eqns (e.g. vol frac, temp)
           mx_ncolacv, ncolacv, finacv, colacv, midacv, &
                                ! Force balance plus cty multi-phase eqns
           nlenmcy, mx_ncolmcy, ncolmcy, finmcy, colmcy, midmcy, &
                                ! Element connectivity
           mxnele, ncolele, midele, finele, colele, &
                                ! Force balance sparsity
           mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, &
                                ! CT sparsity - global cty eqn
           mx_nct, ncolct, findct, colct, &
                                ! C sparsity operating on pressure in force balance
           mx_nc, ncolc, findc, colc, &
                                ! pressure matrix for projection method
           mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, &
                                ! CV-FEM matrix
           mx_ncolm, ncolm, findm, colm, midm, mx_nface_p1 )

      !call check_sparsity( &
      !     u_nonods * nphase * ndim, cv_nonods * nphase, &
      !     u_nonods, cv_nonods, totele, &
      !     mx_ncolacv, ncolacv, finacv, colacv, midacv, & ! CV multi-phase eqns (e.g. vol frac, temp)
      !     nlenmcy, mx_ncolmcy, ncolmcy, finmcy, colmcy, midmcy, & ! Force balance plus cty multi-phase eqns
      !     mxnele, ncolele, midele, finele, colele, & ! Element connectivity 
      !     mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, & ! Force balance sparsity  
      !     mx_nct, ncolct, findct, colct, & ! CT sparsity - global cty eqn
      !     mx_nc, ncolc, findc, colc, & ! C sparsity operating on pressure in force balance
      !     mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, & ! pressure matrix for projection method
      !     mx_ncolm, ncolm, findm, colm, midm ) ! CV-FEM matrix

      !     ewrite(3,*) 'mat_nloc, mat_nonods: ', mat_nloc, mat_nonods
      !     ewrite(3,*) 'mat_ndgln: ', mat_ndgln

      call solve_multiphase_mom_press_volf( state, nphase, ncomp, totele, ndim, &
                                ! Nodes et misc
           u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, mat_nloc, &
           cv_snloc, u_snloc, p_snloc, stotel, &
                                ! Element types
           u_ele_type, p_ele_type, cv_ele_type, &
           cv_sele_type, u_sele_type, &
                                ! Total time loop and initialisation parameters
           ntime_dump, nits, nits_internal, dump_no, &
           nits_flux_lim_volfra, nits_flux_lim_comp, & 
           ndpset, &
                                ! Discretisation parameters
           v_beta, v_theta, &
           v_disopt, &
           v_dg_vel_int_opt, &
           t_beta, t_theta, t_disopt, t_dg_vel_int_opt, lump_eqns, &
           volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux, &
           opt_vel_upwind_coefs, nopt_vel_upwind_coefs, &
           noit_dim, &
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
           x, y, z, xu, yu, zu, nu, nv, nw, ug, vg, wg, &
                                ! Absorption and source terms and coefficients
           u_abs_stab, Mobility, &
           u_absorb, v_absorb, comp_absorb, &
           u_source, v_source, comp_source, &
           t_absorb, t_source, &
                                ! Diffusion parameters
           udiffusion, &
           tdiffusion, &
           comp_diffusion_opt, ncomp_diff_coef, comp_diffusion, comp_diff_coef, &
                                ! Velocities and scalar fields
           u, v, w, &
           den, satura, comp, t, p, cv_p, volfra_pore, &
           perm, &
           uold, vold, wold, denold, saturaold, compold, uden, udenold, deriv, &
           told, pold, cv_pold, nuold, nvold, nwold, &
                                ! EOS terms
           K_Comp, alpha_beta, &
                                ! Matrices sparsity
           mx_ncolacv, ncolacv, finacv, colacv, midacv, & ! CV multi-phase eqns (e.g. vol frac, temp)
           nlenmcy, mx_ncolmcy, ncolmcy, finmcy, colmcy, midmcy, & ! Force balance plus cty multi-phase eqns
           mxnele, ncolele, finele, colele, & ! Element connectivity 
           mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, & ! Force balance sparsity  
           mx_nct, ncolct, findct, colct, & ! CT sparsity - global cty eqn
           mx_nc, ncolc, findc, colc, & ! C sparsity operating on pressure in force balance
           mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, & ! pressure matrix for projection method
           mx_ncolm, ncolm, findm, colm, midm, & ! CV-FEM matrix
           have_temperature_fields, scale_momentum_by_volume_fraction, cv_one, nits_flux_lim_t, t_use_theta_flux, &
           t_get_theta_flux) 

      ewrite(3,*) 'Leaving multiphase_prototype'

    end subroutine multiphase_prototype

    integer function Combination( n, r )
      ! This function performs the combinatorial:
      ! C(n,r) = n! / ( (n-r)! r! )
      integer :: n, r

      Combination = Permut( n ) / ( Permut( n - r ) * Permut( r ))

    end function Combination

    integer function Permut( n )
      ! This function performs probabilistic permutation:
      ! P(n) = n!
      integer :: n
      integer :: i

      permut = 1
      do i = 1, n
         permut = permut * i
      end do

    end function permut

  end module mp_prototype

