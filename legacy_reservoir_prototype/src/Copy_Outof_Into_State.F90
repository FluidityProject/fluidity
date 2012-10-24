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
!    but WITHOUT ANY WARRANTY; without seven the implied warranty of
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
    use field_options
    use spud
    use populate_state_module
    use diagnostic_variables
    use diagnostic_fields
    use diagnostic_fields_wrapper
    use global_parameters, only: option_path_len
    use diagnostic_fields_wrapper_new

    use element_numbering
    use shape_functions
    use fefields
    use boundary_conditions
    use printout

    use quicksort


    implicit none

    private

    public :: copy_outof_state, copy_into_state

   interface Get_Ndgln
     module procedure Get_Scalar_Ndgln, Get_Vector_Ndgln, Get_Mesh_Ndgln
   end interface Get_Ndgln

   interface Get_SNdgln
     module procedure Get_Scalar_SNdgln, Get_Vector_SNdgln
   end interface Get_SNdgln

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
         eos_coefs, cp_coefs,scale_momentum_by_volume_fraction, &
         comp_diff_coef, capil_pres_coef, &
         u_abs_stab, u_absorb, comp_absorb, &
         t_absorb, v_absorb, &
         perm, K_Comp, &
         comp_diffusion, &
                                ! Now adding other things which we have taken inside this routine to define
         cv_nonods, p_nonods, u_nonods, x_nonods, xu_nonods, mat_nonods, &
         have_temperature_fields, nits_flux_lim_t)

      type(state_type), dimension(:), intent(inout) :: state

      logical, intent(inout) :: have_temperature_fields, scale_momentum_by_volume_fraction 

      ! coordinate, velocity and pressure meshes
      type(mesh_type), pointer :: cmesh => null()
      type(mesh_type), pointer :: vmesh => null()
      type(mesh_type), pointer :: pmesh  => null()

      type(scalar_field), pointer :: porosity, permeability, &
           density, pressure, pvf_source, field

      type(vector_field), pointer :: velocity, velocity_source
      type(vector_field), pointer :: velocity_bc

      type(mesh_type), pointer :: velocity_cg_mesh, pressure_cg_mesh
      type(vector_field) :: velocity_cg

      type(tensor_field), pointer :: viscosity_ph1, viscosity_ph2, t_permeability

      integer :: nonlinear_iterations, stat, nstates

      real :: nonlinear_iteration_tolerance

      character(len=OPTION_PATH_LEN) :: option_path

      !! temporary variables only needed for interfacing purposes
      type(vector_field), pointer :: positions

      integer :: i, j, k, l, cv_nonods, p_nonods, mat_nonods, &
           x_nonods, x_nonods_p1, xu_nonods, u_nonods, cv_nod, u_nod, xu_nod, x_nod, mat_nod, &
           shared_nodes, u_nloc2

      real :: coord_min, coord_max, eos_value

      real :: const, dx
      real, dimension(:), allocatable :: const_vec
      real, dimension(:,:), allocatable :: const_array

      integer, dimension(:), allocatable :: cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, &
           xu_ndgln, mat_ndgln, cv_sndgln, p_sndgln, u_sndgln, u_sndgln2, x_sndgln
      integer, dimension(:), pointer :: element_nodes

      real, dimension(:), allocatable :: initial_constant_velocity

      logical :: is_isotropic, is_symmetric, is_diagonal, ndgln_switch

      !! Variables needed by the prototype code
      !! and therefore needing to be pulled out of state or
      !! derived from information in state:

      ! Scalars (from read_scalar())
      integer :: nphases, ncomps, totele, ndim, nlev, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           cv_snloc, u_snloc, p_snloc, stotel, &
           ncoef, ncoef_max, nuabs_coefs, &
           u_ele_type, p_ele_type, mat_ele_type, cv_ele_type, &
           cv_sele_type, u_sele_type, x_snloc

      integer :: u_snloc2, ilev, ele, sele, u_siloc, u_siloc2, inod_remain

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
           t_absorb, v_absorb, vel_absorb, &
           perm, K_Comp

      real, dimension( : , : , : , : ), allocatable :: comp_diffusion

      integer :: Velocity_BC_Type, shape_option(2), kk

      integer, dimension(:), allocatable :: Velocity_SufID_BC

      ! Gravity terms to be linked with u_source
      logical :: have_gravity
      real :: gravity_magnitude, delta_den, grm
      type(vector_field) :: gravity_direction

      integer :: nobcs

      character(len=option_path_len) :: vel_element_type
      integer, dimension( : ), allocatable :: index
      real, dimension( : ), allocatable :: x_temp
      logical :: is_overlapping

      real :: val

      ! velocity bcs
      integer :: count
      integer, dimension(:), allocatable :: face_nodes

      ewrite(3,*) 'In copy_outof_state'

      ! Getting all primary scalars necessary to the model
      call Get_Primary_Scalars( state, &         
           nphases, nstates, ncomps, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, dx, &
           is_overlapping )

      ewrite(3,*) 'nphases, nstates, ncomps, totele, ndim, stotel:', &
           nphases, nstates, ncomps, totele, ndim, stotel
      ewrite(3,*) 'u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, mat_nloc:', &
           u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, mat_nloc
      ewrite(3,*) 'x_snloc, cv_snloc, u_snloc, p_snloc:', &
           x_snloc, cv_snloc, u_snloc, p_snloc
      ewrite(3,*) 'cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, p_nonods :', &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, p_nonods
      ewrite(3,*) 'is_overlapping????:',  is_overlapping

      positions => extract_vector_field( state( 1 ), "Coordinate" )
      allocate( x_ndgln_p1( totele * x_nloc_p1 ) ) ; x_ndgln_p1 = 0
      ewrite(3,*)'X_NDGLN_P1:'
      call Get_Ndgln( x_ndgln_p1, positions )
      ewrite(3,*) '   '


      ! Allocating Volume-based Global Node Numbers:
      ! Positions/Coordinates
      ! NB This is a continuous version of the Pressure mesh
      pressure_cg_mesh => extract_mesh( state( 1 ), "PressureMesh_Continuous" )
      allocate( x_ndgln( totele * x_nloc ) ) ; x_ndgln = 0
      ewrite(3,*)'X_NDGLN:'
      call Get_Ndgln( x_ndgln, pressure_cg_mesh )
      ewrite(3,*) '   '

      ! Pressure and Control Volume
      pressure => extract_scalar_field( state( 1 ), "Pressure" )
      allocate( cv_ndgln( totele * cv_nloc ) ) ; cv_ndgln = 0
      ewrite(3,*)'CV_NDGLN'
      call Get_Ndgln( cv_ndgln, pressure )
      allocate( p_ndgln( totele * p_nloc ) ) ; p_ndgln = 0
      allocate( mat_ndgln( totele * mat_nloc ) ) ; mat_ndgln = 0
      p_ndgln = cv_ndgln
      mat_ndgln = cv_ndgln
      ewrite(3,*) '   '

      ! Velocities
      velocity => extract_vector_field( state( 1 ), "Velocity" )
      allocate( u_ndgln( totele * u_nloc ) ) ;  u_ndgln = 0
      ewrite(3,*)'U_NDGLN'
      call Get_Ndgln( u_ndgln, velocity, is_overlapping, cv_nloc )
      ewrite(3,*) '   '

      ! Velocity in the continuous space
      velocity_cg_mesh => extract_mesh( state( 1 ), "VelocityMesh_Continuous" )
      allocate( xu_ndgln( totele * xu_nloc ) ) ;  xu_ndgln = 0
      ewrite(3,*)'XU_NDGLN'
      call Get_Ndgln( xu_ndgln, velocity_cg_mesh )
      ewrite(3,*) '   '

      ! Allocating Surface-based Global Node Numbers:
      ! Control Volumes
      allocate( cv_sndgln( stotel * cv_snloc ) ) ; cv_sndgln = 0
      ewrite(3,*)'CV_SNGLN'
      call Get_SNdgln( cv_sndgln, pressure )
      ewrite(3,*) '   '

      ! Pressure
      allocate( p_sndgln( stotel * p_snloc ) ) ; p_sndgln = 0
      p_sndgln = cv_sndgln

      ! Velocities
      if ( is_overlapping ) then
         u_snloc2 = u_snloc / cv_nloc
         u_nloc2 = u_nloc / cv_nloc
      else
         u_snloc2 = u_snloc
         u_nloc2 = u_nloc
      end if

      allocate( u_sndgln2( stotel * u_snloc2 ) ) ; u_sndgln2 = 0
      ewrite(3,*)'U_SNGLN'
      call Get_SNdgln( u_sndgln2, velocity )
      allocate( u_sndgln( stotel * u_snloc ) ) ; u_sndgln = 0

      if ( is_overlapping ) then
         ! Convert u_sndgln2 to overlapping u_sndgln
         do sele = 1, stotel
            do u_siloc2 = 1, u_snloc2
               do ilev = 1, cv_nloc
                  u_siloc = ( ilev - 1 ) * u_snloc2 + u_siloc2
                  ele = int( ( u_sndgln2 ( ( sele - 1 ) * u_snloc2 + &
                       u_siloc2 ) - 1 ) / u_nloc2) + 1
                  inod_remain = u_sndgln2( ( sele - 1 ) * u_snloc2 + &
                       u_siloc2 ) - ( ele - 1 ) * u_nloc2
                  u_sndgln( ( sele - 1 ) * u_snloc + u_siloc ) = &
                       ( ele - 1 ) * u_nloc + inod_remain + &
                       ( ilev - 1 ) * u_nloc2
               end do
            end do
         end do
      else
         u_sndgln=u_sndgln2
      end if
      deallocate( u_sndgln2 )

      ewrite(3,*) ' Final u_sndgln: '
      do sele = 1, stotel
         ewrite(3,*) sele, ( u_sndgln( (sele - 1 ) * u_snloc + k ), k = 1, u_snloc )
      end do
      ewrite(3,*) '  '


      ! Velocity on the coordinate mesh (CG)
      allocate( xu( xu_nonods )) ; xu = 0.
      allocate( yu( xu_nonods )) ; yu = 0.
      allocate( zu( xu_nonods )) ; zu = 0.

      ewrite(3,*) ' X/Y/ZU:'

      call allocate( velocity_cg, ndim, velocity_cg_mesh, "Velocity_CG_Coordinates" )
      velocity_cg % val( :, : )= 0
      call project_field( positions, velocity_cg, positions )

      !! Global numbering to allow proper copying of fields
      do k = 1, totele
         element_nodes => ele_nodes( velocity_cg_mesh, k )
         do j = 1, xu_nloc
            xu( element_nodes( j ) ) = velocity_cg % val( 1, element_nodes( j ) )
            if( ndim > 1 ) yu( element_nodes( j ) ) = velocity_cg % val( 2, element_nodes( j ) )
            if( ndim > 2 ) zu( element_nodes( j ) ) = velocity_cg % val( 3, element_nodes( j ) )
            ewrite(3,*)'ele, ndgln, xu, yu:', k, xu_ndgln( ( k - 1 ) * xu_nloc + j ), &
                 xu( xu_ndgln( ( k - 1 ) * xu_nloc + j ) ), &
                 yu( xu_ndgln( ( k - 1 ) * xu_nloc + j ))
         end do
         ewrite(3,*) ' '
      end do
      call deallocate( velocity_cg )

      ! Coordinate mesh
      allocate( x( x_nonods ) ) ; x = 0.
      allocate( y( x_nonods ) ) ; y = 0.
      allocate( z( x_nonods ) ) ; z = 0.

      call xp1_2_xp2( ndim, totele, x_nloc, x_nloc_p1, x_nonods_p1, x_nonods, &
           x_ndgln_p1, x_ndgln, positions, &
           x, y, z )

      do k = 1, totele
         do j = 1, x_nloc
            ewrite(3,*)'ele, ndgln, x, y:', k, x_ndgln( ( k - 1 ) * x_nloc + j ), &
                 x( x_ndgln( ( k - 1 ) * x_nloc + j ) ), &
                 y( x_ndgln( ( k - 1 ) * x_nloc + j ))
         end do
      end do

      ! Defining element type 
      call Get_Ele_Type(x_nloc, cv_ele_type, p_ele_type, u_ele_type)
      ewrite(3,*)'cv_ele_type:', cv_ele_type

      ! Defining  solvers options
      call Get_Solvers_Options( nphases, &
           nits, nits_internal, &
           ndpset, nits_flux_lim_volfra, nits_flux_lim_comp, &
           nits_flux_lim_t, t_disopt, v_disopt, t_beta, v_beta, &
           t_theta, v_theta, u_theta, &
           nonlinear_iteration_tolerance, lump_eqns )

      ! Initial pressure -- this can be deleted later on.
      call get_option('/material_phase[0]/scalar_field::Pressure/prognostic/' // &
           'atmospheric_pressure', patmos, default=0.0)
      call get_option('/material_phase[0]/scalar_field::Pressure/prognostic/' // &
           'initial_condition::WholeMesh/constant',p_ini, default=0.0)


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

      ! IN/DG_ELE_UPWIND are options for optimisation of upwinding across faces in the overlapping
      ! formulation. The data structure and options for this formulation need to be added later. 
      in_ele_upwind =  3
      dg_ele_upwind = 3

      ! Variables *ENTIRELY HARD-WIRED* that *MUST BE REPLACED* in the future
      t_dg_vel_int_opt = 1
      u_dg_vel_int_opt = 4 ! Not used -- it can be deleted
      v_dg_vel_int_opt = 4
      if (.not.is_overlapping) v_dg_vel_int_opt = 1
      w_dg_vel_int_opt = 0 ! Not used -- it can be deleted
      comp_diffusion_opt = 0
      ncomp_diff_coef = 0
      volfra_use_theta_flux = .false.
      volfra_get_theta_flux = .true.
      comp_use_theta_flux = .false.
      comp_get_theta_flux = .true.
      t_use_theta_flux = .false.
      t_get_theta_flux = .true.
      nlev = 3
      mat_ele_type = 1 ; u_sele_type = 1 ; cv_sele_type = 1


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
      allocate(volfra_pore(totele)) ; volfra_pore = 0.
      pore_ele_loop: do k = 1,element_count(porosity)
         element_nodes => ele_nodes( porosity, k )
         volfra_pore(k) = porosity%val( element_nodes( 1 ))
      end do pore_ele_loop

      ! Assuming for now that permeability is constant across an element
      if (have_option("/porous_media/scalar_field::Permeability")) then

         permeability => extract_scalar_field(state(1), "Permeability")
         allocate(perm(totele, ndim, ndim)) ; perm = 0.
         perm_ele_loop: do k = 1,element_count(permeability)
            element_nodes => ele_nodes(permeability,k)
            forall(i=1:ndim)perm(k, i, i) = permeability%val(element_nodes(1))
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

      else

         FLAbort("Unknown permeability option")

      end if


!!!
!!! Robin Boundary conditions need to be added at a later stage
!!!
      allocate( suf_u_bc_rob1( stotel * u_snloc * nphases ))
      allocate( suf_v_bc_rob1( stotel * u_snloc * nphases ))
      allocate( suf_w_bc_rob1( stotel * u_snloc * nphases ))
      allocate( suf_u_bc_rob2( stotel * u_snloc * nphases ))
      allocate( suf_v_bc_rob2( stotel * u_snloc * nphases ))
      allocate( suf_w_bc_rob2( stotel * u_snloc * nphases ))
      allocate( suf_t_bc_rob1( stotel * cv_snloc * nphases ))
      allocate( suf_t_bc_rob2( stotel * cv_snloc * nphases ))
      allocate( suf_vol_bc_rob1( stotel * cv_snloc * nphases ))
      allocate( suf_vol_bc_rob2( stotel * cv_snloc * nphases ))
      allocate( suf_comp_bc_rob1( stotel * cv_snloc * nphases ))
      allocate( suf_comp_bc_rob2( stotel * cv_snloc * nphases ))
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
!!! Initialise fields and specify BCs for all fields
!!! BC_Type :: == 1 (Dirichlet), == 2 (Robin), == 3 (Newman) 
!!!

!!!
!!! Components
!!!
      ewrite(3,*) "Allocating Components:"

      allocate( comp( cv_nonods * nphases * ncomps )) ; comp = 0.
      allocate( comp_source( cv_nonods * nphases )) ; comp_source = 0.
      allocate( wic_comp_bc( stotel * nphases )) ; wic_comp_bc = 0
      allocate( suf_comp_bc( stotel * cv_snloc * nphases * ncomps )) ; suf_comp_bc = 0.

      Loop_Components: do i = nphases + 1, nphases + ncomps ! Component loop
         Loop_Phases_Components: do j = 1, nphases ! Phase loop

            field => extract_scalar_field( state( i ), &
                 "ComponentMassFractionPhase" // int2str( j ) )

            k = ( i - ( nphases + 1 ) ) * nphases * cv_nonods + &
                 ( j - 1 ) * cv_nonods
            kk = ( i - ( nphases + 1 ) ) * nphases * stotel * cv_snloc + &
                 ( j - 1 ) * stotel * cv_snloc


ewrite(3,*)'@@:', node_count( field )
ewrite(3,*)'-->:',k + 1, k + node_count( field ), kk + 1, kk + stotel * cv_snloc

            call Get_CompositionFields_Outof_State( state, nphases, i, j, field, &
                 comp( k + 1 : k + cv_nonods ), wic_comp_bc, &
                 kk + 1, kk + stotel * cv_snloc, &
                 suf_comp_bc( kk + 1 : kk + stotel * cv_snloc ), &
                 field_prot_source = &
                 comp_source( ( j - 1 ) * cv_nonods + 1 : &
                 ( j - 1 ) * cv_nonods + cv_nonods ) )

         end do Loop_Phases_Components
      end do Loop_Components
      !  ewrite(3,*)'comp:', comp
      !  ewrite(3,*)'wic_comp_bc:', wic_comp_bc
      !  ewrite(3,*)'suf_comp_bc:', suf_comp_bc
      !  ewrite(3,*)'comp_source:', comp_source


      ewrite(3,*)'comp:', cv_nonods, nphases, ncomps
      do i = 1, ncomps
         do j = 1, nphases
            ewrite(3,*) i, j, comp( ( i - 1 ) * nphases * cv_nonods + ( j - 1 ) * cv_nonods + 1 : &
                 ( i - 1 ) * nphases * cv_nonods + ( j - 1 ) * cv_nonods + cv_nonods )
         end do
      end do

      ewrite(3,*)'suf_comp_bc:', stotel, cv_snloc
      do i = 1, ncomps
         do j = 1, nphases
            ewrite(3,*) i, j, suf_comp_bc( ( i - 1 ) * nphases * stotel * cv_snloc + &
                 ( j - 1 ) * stotel * cv_snloc + 1 : &
                 ( i - 1 ) * nphases * stotel * cv_snloc + &
                 ( j - 1 ) * stotel * cv_snloc + stotel * cv_snloc )
            ewrite(3,*)' '
         end do
      end do

!!!  Density
!!!
      ewrite(3,*) "Densities"

      Loop_Density: do i = 1, nphases
         field => extract_scalar_field( state( i ), "Density")
         if ( .not. allocated( den )) then
            allocate( den( nphases * node_count( field ))) ; den = 0.
         endif
         if ( .not. allocated( wic_d_bc )) then
            allocate( wic_d_bc( stotel * nphases ) ) ; wic_d_bc = 0
         end if
         if ( .not. allocated( suf_d_bc )) then 
            allocate( suf_d_bc( stotel * p_snloc * nphases ) ) ; suf_d_bc = 0.
         end if
         k = ( i - 1 ) * node_count( field )
         call Get_ScalarFields_Outof_State( state, i, field, &
              den(  k + 1 : k + node_count( field ) ), wic_d_bc, suf_d_bc )
      end do Loop_Density

!!!
!!! Pressure
!!!
      ewrite(3,*) "pressure"

      field => extract_scalar_field( state( 1 ), "Pressure" )
      allocate( p( node_count( field ))) ; p=0.

      if( .not. ( allocated( wic_p_bc ) .and. allocated( suf_p_bc ))) then
         allocate( wic_p_bc( stotel * nphases )) ; wic_p_bc = 0
         allocate( suf_p_bc( stotel * p_snloc * nphases )) ; suf_p_bc = 0.
      end if

      call Get_ScalarFields_Outof_State( state, 1, field, &
           p, wic_p_bc, suf_p_bc )

      ! Control-volume pressure
      allocate( cv_p( node_count( field )))
      cv_p = p

      ! Copy this to the other phases
      if(nphases > 1) then
         do i = 2, nphases
            wic_p_bc( ( i - 1 ) * stotel + 1 : i * stotel ) = wic_p_bc( 1 : stotel )
            suf_p_bc( ( i - 1 ) * stotel * p_snloc + 1 : i * stotel * p_snloc ) = suf_p_bc( 1 : stotel * p_snloc)
         end do
      end if

!!!
!!! Volume Fraction (or Saturation)
!!!
      ewrite(3,*) "phasevolumefraction"

      Loop_VolumeFraction: do i = 1, nphases

         field => extract_scalar_field(state(i), "PhaseVolumeFraction")
         if (.not. allocated(satura)) then
            allocate( satura( nphases * node_count( field ))) ; satura=0.
         end if
         if (.not. allocated(wic_vol_bc)) then
            allocate( wic_vol_bc( stotel * nphases ))
            wic_vol_bc = 0
         endif
         if (.not. allocated(suf_vol_bc)) then
            allocate( suf_vol_bc( stotel * cv_snloc * nphases ))
            suf_vol_bc = 0.
         endif
         if(.not. allocated( v_source ) ) then
            allocate( v_source( cv_nonods * nphases ) ) ; v_source = 0.
         end if
         k = ( i - 1 ) * node_count( field )
         call Get_ScalarFields_Outof_State( state, i, field, &
              satura( k + 1 : k + node_count( field ) ), wic_vol_bc, suf_vol_bc, &
              field_prot_source = v_source )

      enddo Loop_VolumeFraction

      ! Also need to allocate and initialise volfra
      allocate(volfra( cv_nonods * nphases ))
      volfra = 0.

!!!
!!! Velocity and associated boundary conditions:
!!!
      ewrite(3,*) "velocity"

      if (.not. allocated(u)) then
         allocate( u( nphases * u_nonods )) ; u = 0.
         allocate( v( nphases * u_nonods )) ; v = 0.
         allocate( w( nphases * u_nonods )) ; w = 0. 
      end if

      if (.not.allocated(wic_u_bc)) then
         allocate( wic_u_bc( stotel * nphases ))
         wic_u_bc = 0
      end if
      if (.not.allocated(suf_u_bc)) then
         allocate( suf_u_bc( stotel * u_snloc * nphases ))
         suf_u_bc = 0.
      end if
      if (.not.allocated(suf_v_bc)) then
         allocate( suf_v_bc( stotel * u_snloc * nphases ))
         suf_v_bc = 0.
      end if
      if (.not.allocated(suf_w_bc)) then
         allocate( suf_w_bc( stotel * u_snloc * nphases ))
         suf_w_bc = 0.
      end if

      allocate( nu(nphases * u_nonods )) ; nu = 0.
      allocate( nv(nphases * u_nonods )) ; nv = 0.
      allocate( nw(nphases * u_nonods )) ; nw = 0.

      allocate( ug(nphases * u_nonods )) ; ug = 0.
      allocate( vg(nphases * u_nonods )) ; vg = 0.
      allocate( wg(nphases * u_nonods )) ; wg = 0.


      allocate( u_absorb( cv_nloc * totele, ndim * nphases, &
           ndim * nphases ) ) ; u_absorb = 0.
      allocate( vel_absorb( u_nloc * totele, ndim * nphases, &
           ndim * nphases ) ) ; vel_absorb = 0.

      Loop_Velocity: do i = 1, nphases

         velocity => extract_vector_field(state(i), "Velocity")
         call Get_VectorFields_Outof_State( state, i, velocity, &
              u, v, w, nu, nv, nw, &
              wic_u_bc, suf_u_bc, suf_v_bc, suf_w_bc, u_source, vel_absorb, is_overlapping )

      end do Loop_Velocity    

      allocate( comp_absorb( cv_nonods, nphases, nphases ) ) 
      allocate( t_absorb( cv_nonods, nphases, nphases ) ) 
      allocate( v_absorb( cv_nonods, nphases, nphases ) ) 
      comp_absorb=0. ; t_absorb=0. ; v_absorb=0.

      allocate(uabs_option(nphases))
      allocate(eos_option(nphases))
      allocate(cp_option(nphases))
      do i=1,nphases
         option_path = "/material_phase[" // int2str(i-1) // "]/multiphase_properties/relperm_type"
         if (have_option(trim(option_path)//"/Corey")) uabs_option(i)=3
         if (have_option(trim(option_path)//"/Corey/boost_at_zero_saturation")) uabs_option(i)=5

!!$         option_path = "/material_phase[" // int2str(i-1) // "]/equation_of_state"
!!$         if (have_option(trim(option_path)//"/incompressible/linear")) then
!!$            eos_option(i) = 2
!!$         elseif (have_option(trim(option_path)//"/compressible/stiffened_gas")) then
!!$            eos_option(i) = 1
!!$         elseif (have_option(trim(option_path)//"/compressible/exponential_oil_gas")) then
!!$            eos_option(i) = 3
!!$         elseif (have_option(trim(option_path)//"/compressible/linear_in_pressure")) then
!!$            ! do nothing here
!!$         elseif (have_option(trim(option_path)//"/compressible/exponential_in_pressure")) then
!!$            ! do nothing here
!!$         else
!!$            FLAbort("Unknown EoS option for phase "// int2str(i))
!!$         endif

         cp_option(i) = 1
         option_path = "/material_phase[" // int2str(i-1) // "]/scale_momentum_by_volume_fraction"
         if(have_option(trim(option_path))) scale_momentum_by_volume_fraction = .true.
      enddo


      !! uabs_coefs is currently only used in rel perm options which aren't
      !! selected in any of the test cases, ie the 'standard polynomial
      !! representation of relative permeability'
      nuabs_coefs = 1
      ncoef = 10
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
               call get_option("/material_phase[" // int2str(i-1) // "]/equation_of_state/incompressible/linear/specify_all/coefficients", eos_coefs(i, :))
            endif
         elseif (eos_option(i)==1) then
            call get_option("/material_phase[" // int2str(i-1) // "]/equation_of_state/compressible/stiffened_gas/eos_option1", eos_coefs(i, 1))
            call get_option("/material_phase[" // int2str(i-1) // "]/equation_of_state/compressible/stiffened_gas/eos_option2", eos_coefs(i, 2))
            eos_coefs(i, 3:ncoef) = 0.
         endif

      end do
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
            gravity_direction%val(i,:) = gravity_direction%val(i,:) / sqrt( grm )
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

      !density => extract_scalar_field(state(1), "Density")
      !if (have_gravity) then
      !   do i = 1, nphases - 1
      !      delta_den = 0.
      !      do j = 2, node_count( density )
      !         delta_den = delta_den + ( den((i)*node_count(density)+j) - &
      !              den((i-1)*node_count(density)+j) ) / &
      !              real( ( node_count( density ) - 1 ) * max( 1, ( nphases - 1 )))
      !      end do
      !      do j = 1, u_nonods
      !         do k = 1, ndim
      !            u_source( ( k - 1 ) * u_nonods * nphases + ( i - 1 ) * u_nonods + j  ) = &
      !                 u_source( ( k - 1 ) * u_nonods * nphases + ( i - 1 ) * u_nonods + j  ) + &
      !                 delta_den * gravity_magnitude * gravity_direction%val(k,1) * min( 0.027, & 
      !                 domain_length /  real( totele * cv_nloc ))
      !         end do
      !      end do
      !   end do
      !end if


!!!
!!!  Temperature
!!!

      ewrite(3,*) 'Getting temperature field'

      allocate( t( cv_nonods * nphases ))
      t = 0.
      have_temperature_fields = .false.

      Conditional_ExtraScalarField: if( have_option( "/material_phase[0]/" // &
           "scalar_field::Temperature" ))then

         have_temperature_fields = .true.

         Loop_Temperature: do i = 1, nphases

            field => extract_scalar_field(state(i), "Temperature")

            if ( .not. allocated( t_source )) then
               allocate( t_source( cv_nonods * nphases )) ; t_source = 0.
            end if
            if ( (.not. allocated( wic_t_bc )) .and. (.not. allocated( suf_t_bc ))) then
               allocate( wic_t_bc( stotel * nphases )) ; wic_t_bc = 0
               allocate( suf_t_bc( stotel * cv_nloc * nphases )) ; suf_t_bc = 0.
            endif
            k = ( i - 1 ) * node_count( field )
            call Get_ScalarFields_Outof_State( state, i, field, &
                 t( k + 1 : k + node_count( field ) ), wic_t_bc, suf_t_bc, &
                 field_prot_source=t_source( k + 1 : k + node_count( field ) ))

         end do Loop_Temperature

      end if Conditional_ExtraScalarField

      ewrite(3,*) 'Getting capillary pressure options and absorptions'
      if (ncapil_pres_coef>0) then
         allocate(capil_pres_coef(ncapil_pres_coef,nphases,nphases))
         capil_pres_coef=0.
      endif
      ! Not sure what this one does yet
      allocate(u_abs_stab(cv_nloc*totele, ndim*nphases, ndim*nphases))
      u_abs_stab=0.
      ! or this one
      allocate( k_comp( ncomps, nphases, nphases ))
      k_comp=0.
      if (KComp_Sigmoid) then
         do i=1, ncomps
            call get_option('material_phase['// int2str(i+nphases-1) //']/is_multiphase_component/' // &
                 'KComp_Sigmoid/K_Comp', k_comp(i, 1, 1))
            k_comp(i, 1:nphases, 1:nphases) = k_comp(i, 1, 1)
         end do
      end if

      allocate( comp_diffusion( cv_nloc*totele, ndim, ndim, nphases ))
      comp_diffusion=0.0
      allocate( comp_diff_coef( ncomps, ncomp_diff_coef, nphases ))
      comp_diff_coef=0.0
      allocate( cv_one( nphases * cv_nonods ))
      cv_one = 1.0

      if( .false. ) then

         allocate( x_sndgln( stotel * cv_snloc ) ) ; x_sndgln = 0
         call Get_SNdgln( x_sndgln, positions )

         call test_bc( ndim, nphases, &
              u_nonods, cv_nonods, x_nonods, &
              u_snloc, p_snloc, cv_snloc, stotel, u_sndgln, p_sndgln, cv_sndgln, x_sndgln, &
                                ! For force balance eqn:
              suf_u_bc, suf_v_bc, suf_w_bc, suf_p_bc, &
              suf_u_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob1, suf_v_bc_rob2, &
              suf_w_bc_rob1, suf_w_bc_rob2, &
              wic_u_bc, wic_p_bc, &
                                ! For cty eqn: 
              suf_vol_bc, suf_d_bc, &
              suf_vol_bc_rob1, suf_vol_bc_rob2, &
              wic_vol_bc, wic_d_bc, &
              x, y, z )

         deallocate( x_sndgln )

      end if

      ewrite(3,*) "Leaving copy_outof_state"

      return

    end subroutine Copy_outof_state

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
      integer :: nloc, nlev
      integer, dimension(:), pointer :: element_nodes => null()
      type(scalar_field), pointer :: phasevolumefraction => null()
      type(scalar_field), pointer :: phasetemperature => null()      
      type(scalar_field), pointer :: pressure => null()
      type(vector_field), pointer :: velocity => null()
      type(scalar_field), pointer :: density => null()
      type(scalar_field), pointer :: componentmassfraction => null()
      character(len=option_path_len) :: material_phase_name, vel_element_type
      logical :: is_overlapping
      real, dimension(:), allocatable :: proto_velocity_u_tmp, proto_velocity_v_tmp, proto_velocity_w_tmp
      ewrite(3,*) "In copy_into_state"

      assert(size(state) >= nphase)

      ! Deal with overlapping velocity...
      ! Get the vel element type.
      call get_option('/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           vel_element_type)
      is_overlapping = .false.
      if ( trim( vel_element_type ) == 'overlapping' ) is_overlapping = .true. 

      velocity => extract_vector_field(state(1), "Velocity", stat=stat)
      if (stat /= 0) then 
         FLAbort('Failed to extract phase velocity from state in copy_into_state')
      end if
      number_nodes = node_count( velocity )

      allocate( proto_velocity_u_tmp( number_nodes * nphase ) ) ; proto_velocity_u_tmp = 0.
      allocate( proto_velocity_v_tmp( number_nodes * nphase ) ) ; proto_velocity_v_tmp = 0.
      allocate( proto_velocity_w_tmp( number_nodes * nphase ) ) ; proto_velocity_w_tmp = 0.

      if ( is_overlapping ) then

         ! in case of overlapping elements just take
         ! the average of the various levels
         !
         ! this means that this field cannot be used 
         ! for any cv-fem calculations
         !
         ! ONLY for visualisation purposes

         pressure => extract_scalar_field(state(1), "Pressure")
         nlev = ele_loc( pressure, 1 )
         nloc = ele_loc( velocity, 1 )

         do p = 1, nphase
            do i = 1, element_count( velocity )
               do j = 1, nlev
                  proto_velocity_u_tmp ( (p-1)*number_nodes + (i-1)*nloc + 1 : (p-1)*number_nodes + i*nloc ) = &
                       proto_velocity_u_tmp ( (p-1)*number_nodes + (i-1)*nloc + 1 : (p-1)*number_nodes + i*nloc ) + &
                       proto_velocity_u ( (p-1)*number_nodes*nlev + (i-1)*nloc*nlev + (j-1)*nloc+ 1 : & 
                       &                        (p-1)*number_nodes*nlev + (i-1)*nloc*nlev + j*nloc )

                  proto_velocity_v_tmp ( (p-1)*number_nodes + (i-1)*nloc + 1 : (p-1)*number_nodes + i*nloc ) = &
                       proto_velocity_v_tmp ( (p-1)*number_nodes + (i-1)*nloc + 1 : (p-1)*number_nodes + i*nloc ) + &
                       proto_velocity_v ( (p-1)*number_nodes*nlev + (i-1)*nloc*nlev + (j-1)*nloc+ 1 : & 
                       &                        (p-1)*number_nodes*nlev + (i-1)*nloc*nlev + j*nloc )

                  proto_velocity_w_tmp ( (p-1)*number_nodes + (i-1)*nloc + 1 : (p-1)*number_nodes + i*nloc ) = &
                       proto_velocity_w_tmp ( (p-1)*number_nodes + (i-1)*nloc + 1 : (p-1)*number_nodes + i*nloc ) + &
                       proto_velocity_w ( (p-1)*number_nodes*nlev + (i-1)*nloc*nlev + (j-1)*nloc+ 1 : & 
                       &                        (p-1)*number_nodes*nlev + (i-1)*nloc*nlev + j*nloc )
               end do
            end do
         end do

         proto_velocity_u_tmp = proto_velocity_u_tmp / nlev
         proto_velocity_v_tmp = proto_velocity_v_tmp / nlev
         proto_velocity_w_tmp = proto_velocity_w_tmp / nlev

      else
         proto_velocity_u_tmp = proto_velocity_u
         proto_velocity_v_tmp = proto_velocity_v
         proto_velocity_w_tmp = proto_velocity_w
      end if

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
                    proto_velocity_u_tmp( (i-1)*nloc + j + (p-1)*number_nodes) )

               ! set v
               if (ndim > 1) call set(velocity, &
                    2, &
                    element_nodes(j), &
                    proto_velocity_v_tmp( (i-1)*nloc + j + (p-1)*number_nodes) )

               ! set w
               if (ndim > 2) call set(velocity, &
                    3, &
                    element_nodes(j), &
                    proto_velocity_w_tmp( (i-1)*nloc + j + (p-1)*number_nodes) )

            end do vel_node_loop

         end do vel_ele_loop

      end do phase_loop

      deallocate( proto_velocity_u_tmp, proto_velocity_v_tmp, proto_velocity_w_tmp )

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


    subroutine Get_Primary_Scalars( state, &         
         nphases, nstates, ncomps, totele, ndim, stotel, &
         u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
         x_snloc, cv_snloc, u_snloc, p_snloc, &
         cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, dx, &
         is_overlapping )
      !! This subroutine extracts all primary variables associated with the mesh from state,
      !! and associated them with the variables used in the MultiFluids model.
      implicit none
      type(state_type), dimension(:), intent( in ) :: state
      integer, intent( inout ) :: nphases, nstates, ncomps, totele, ndim, &
           stotel, u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, cv_nonods, mat_nonods, &
           u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods
      real, intent( inout ) :: dx
      logical, intent( inout ) :: is_overlapping

      ! Local variables
      character( len = option_path_len ) :: vel_element_type
      type( vector_field ), pointer :: positions, velocity
      type( scalar_field ), pointer :: pressure
      type( mesh_type ), pointer :: velocity_cg_mesh, pressure_cg_mesh
      integer :: i, degree

      ewrite(3,*)' In Get_Primary_Scalars'

      ! Defining dimension and nstates
      call get_option( "/geometry/dimension", ndim )
      nstates = option_count( "/material_phase" )

      ! Assume there are the same number of components in each phase
      ! (will need to check this eventually)
      ncomps = 0
      do i = 1, nstates
         if( have_option( "/material_phase[" // int2str(i-1) // &
              "]/is_multiphase_component" ) ) then
            ncomps = ncomps + 1
         end if
      end do
      nphases = nstates - ncomps
      assert( nphases > 0 ) ! Check if there is more than 0 phases

      ! Get the vel element type.
      call get_option('/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           vel_element_type)
      is_overlapping = .false.
      if ( trim( vel_element_type ) == 'overlapping' ) is_overlapping = .true. 

      positions => extract_vector_field( state, "Coordinate" )
      pressure_cg_mesh => extract_mesh( state, "PressureMesh_Continuous" )

      ! Defining number of elements and surface elements, coordinates, locs and snlocs
      totele = ele_count( positions )
      stotel = surface_element_count( positions )

      ! Coordinates
      x_nloc_p1 = ele_loc( positions, 1 )
      x_nloc = ele_loc( pressure_cg_mesh, 1 )
      x_snloc = face_loc( pressure_cg_mesh, 1 )
      x_nonods_p1 = node_count( positions )
      x_nonods = node_count( pressure_cg_mesh )

      ! Pressure, Control Volumes and Materials
      pressure => extract_scalar_field( state, "Pressure" )
      p_nloc = ele_loc( pressure, 1 )
      p_snloc = face_loc( pressure, 1 )
      p_nonods = node_count( pressure )
      cv_nloc = p_nloc
      cv_snloc = p_snloc
      cv_nonods = p_nonods
      mat_nloc = cv_nloc
      mat_nonods = mat_nloc * totele

      ! Velocities and velocities (DG) associated with the continuous space (CG)
      velocity => extract_vector_field( state, "Velocity" )
      u_nloc = ele_loc( velocity, 1 )
      u_snloc = face_loc( velocity, 1 )
      u_nonods = node_count( velocity )

      ! Get the continuous space of the velocity field
      velocity_cg_mesh => extract_mesh( state, "VelocityMesh_Continuous" )
      xu_nloc = ele_loc( velocity_cg_mesh, 1 )
      xu_nonods = max(( xu_nloc - 1 ) * totele + 1, totele )

      ! Take care of overlapping elements
      if( ( is_overlapping ) .and. ( ndim > 1 ) ) u_nonods = u_nonods * cv_nloc
      if( ( is_overlapping ) .and. ( ndim > 1 ) ) u_nloc = u_nloc * cv_nloc
      if( is_overlapping ) u_snloc = u_snloc * cv_nloc

      ! Used just for 1D:
      dx = maxval( positions%val(1,:) ) - minval( positions%val(1,:) )

      return
    end subroutine Get_Primary_Scalars

    subroutine Get_Scalar_Ndgln( ndgln, field, is_overlapping, cv_nloc )
      implicit none
      type( scalar_field ), intent( in ) :: field
      integer, intent( in ), optional :: cv_nloc
      logical, intent( in ), optional :: is_overlapping
      integer, dimension( : ), intent( inout ) :: ndgln
      ! Local variables
      integer, dimension( : ), pointer :: nloc
      integer :: ele, iloc

      do ele = 1, ele_count( field )
         nloc => ele_nodes( field, ele )
         do iloc = 1, ele_loc( field, ele )
            ndgln( ( ele - 1 ) * ele_loc( field, ele ) + iloc ) =  nloc( iloc )
            ewrite(3,*)'ele, iloc, ndgln:', ele, iloc, &
                 ndgln( ( ele - 1 ) * ele_loc( field, ele ) + iloc )
         end do
      end do

      return
    end subroutine Get_Scalar_Ndgln

    subroutine Get_Vector_Ndgln( ndgln, field, is_overlapping, cv_nloc )
      implicit none
      type( vector_field ), intent( in ) :: field
      integer, intent( in ), optional :: cv_nloc
      logical, intent( in ), optional :: is_overlapping
      integer, dimension( : ), intent( inout ) :: ndgln
      ! Local variables
      integer, dimension( : ), pointer :: nloc
      integer :: ele, iloc, count, cv_nloc2, ndim
      logical :: is_overlapping2

      
      call get_option( "/geometry/dimension", ndim )
      is_overlapping2 = .false.
      if ( present( is_overlapping ) ) is_overlapping2 = is_overlapping
      cv_nloc2 = 1
      if ( ( is_overlapping2 ) .and. ( ndim > 1 ) ) cv_nloc2 = cv_nloc

      count = 0
      do ele = 1, ele_count( field )
         nloc => ele_nodes( field, ele )
         do iloc = 1, ele_loc( field, 1 ) * cv_nloc2
            if( is_overlapping2 ) then
               count = count + 1
               ndgln( ( ele - 1 ) * ele_loc( field, 1 ) * cv_nloc2 + iloc ) = count
            else
               ndgln( ( ele - 1 ) * ele_loc( field, 1 ) * cv_nloc2  + iloc ) =  nloc( iloc )
            end if
            ewrite(3,*)'ele, iloc, ndgln:', ele, iloc, &
                 ndgln( ( ele - 1 ) * ele_loc( field, 1 ) * cv_nloc2  + iloc )
         end do
      end do

      return
    end subroutine Get_Vector_Ndgln


   subroutine Get_Mesh_Ndgln( ndgln, mesh, is_overlapping, cv_nloc )
      implicit none
      type( mesh_type ), intent( in ) :: mesh
      integer, intent( in ), optional :: cv_nloc
      logical, intent( in ), optional :: is_overlapping
      integer, dimension( : ), intent( inout ) :: ndgln
      ! Local variables
      integer, dimension( : ), pointer :: nloc
      integer :: ele, iloc

      do ele = 1, ele_count( mesh )
         nloc => ele_nodes( mesh, ele )
         do iloc = 1, ele_loc( mesh, ele )
            ndgln( ( ele - 1 ) * ele_loc( mesh, ele ) + iloc ) =  nloc( iloc )
            ewrite(3,*)'ele, iloc, ndgln:', ele, iloc, &
                 ndgln( ( ele - 1 ) * ele_loc( mesh, ele ) + iloc )
         end do
      end do

      return
    end subroutine Get_Mesh_Ndgln


    subroutine Get_Scalar_SNdgln( sndgln, field, is_overlapping, cv_nloc  )
      implicit none
      type( scalar_field ), intent( in ) :: field
      integer, dimension( : ), intent( inout ) :: sndgln
      integer, intent( in ), optional :: cv_nloc
      logical, intent( in ), optional :: is_overlapping
      ! Local variables
      integer, dimension( : ), allocatable :: snloc
      integer :: sele, iloc

      allocate( snloc( face_loc( field, 1 ) ) )
      do sele = 1, surface_element_count( field )
         snloc = face_global_nodes( field, sele )
         do iloc = 1, face_loc( field, sele )
            sndgln( ( sele - 1 ) * face_loc( field, sele ) + iloc ) =  snloc( iloc )
            ewrite(3,*)'sele, iloc, sndgln:', sele, iloc, &
                 sndgln( ( sele - 1 ) * face_loc( field, sele ) + iloc )
         end do
      end do

      deallocate( snloc )

      return
    end subroutine Get_Scalar_SNdgln

    subroutine Get_Vector_SNdgln( sndgln, field, is_overlapping, cv_nloc  )
      implicit none
      type( vector_field ), intent( in ) :: field
      integer, dimension( : ), intent( inout ) :: sndgln
      integer, intent( in ), optional :: cv_nloc
      logical, intent( in ), optional :: is_overlapping
      ! Local variables
      integer, dimension( : ), allocatable :: snloc
      integer :: sele, iloc

      allocate( snloc( face_loc( field, 1 ) ) )
      do sele = 1, surface_element_count( field )
         snloc = face_global_nodes( field, sele )
         do iloc = 1, face_loc( field, sele )
            sndgln( ( sele - 1 ) * face_loc( field, sele ) + iloc ) =  snloc( iloc )
            ewrite(3,*)'sele, iloc, sndgln:', sele, iloc, &
                 sndgln( ( sele - 1 ) * face_loc( field, sele ) + iloc )
         end do
      end do

      deallocate( snloc )

      return
    end subroutine Get_Vector_SNdgln

    subroutine Get_CompositionFields_Outof_State( state, nphases, icomp, iphase, field, &
         field_prot, wic_bc, &
         kprime, kprime2, &
         suf_bc, &
         field_prot_source, field_prot_absorption )
      implicit none
      type( state_type ), dimension( : ), intent( inout ) :: state
      integer, intent( in ) :: nphases, icomp, iphase
      type( scalar_field ), pointer :: field, field_prot_bc
      real, dimension( : ), intent( inout ) :: field_prot
      real, dimension( : ), intent( inout ), optional :: field_prot_source, field_prot_absorption
      integer, dimension( : ), intent( inout ) :: wic_bc
      integer, intent( in ) :: kprime, kprime2
      real, dimension( kprime : kprime2 ), intent( inout ) :: suf_bc
      ! Local variables
      type( mesh_type ), pointer :: pmesh, cmesh
      type(scalar_field), pointer :: pressure, field_source, field_absorption
      type( scalar_field ) :: dummy
      type( vector_field ), pointer :: positions
      integer, dimension( : ), allocatable :: sufid_bc, face_nodes
      integer :: shape_option( 2 )
      character( len = option_path_len ) :: option_path, field_name
      logical :: have_source, have_absorption
      integer :: nstates, stotel, nobcs, bc_type, i, j, k, kk, sele, stat, snloc
      real :: initial_constant
      character( len = 8192 ) :: func

      field_name = trim( field % name )
      positions => extract_vector_field( state( 1 ), "Coordinate" )
      pressure => extract_scalar_field( state( 1 ), "Pressure" )
      pmesh => extract_mesh( state, "PressureMesh" )
      cmesh => extract_mesh( state, "CoordinateMesh" )

      field_source => extract_scalar_field( state( iphase ), field_name // "Source", stat )
      have_source = ( stat == 0 )
      field_absorption => extract_scalar_field( state( iphase ), field_name // "Absorption", stat )
      have_absorption = ( stat == 0 )

      snloc = face_loc( pressure, 1 )
      nstates = option_count("/material_phase")
      stotel = surface_element_count( cmesh )

      option_path = "/material_phase[" // int2str( icomp - 1 ) // &
           "]/scalar_field::ComponentMassFractionPhase" // &
           int2str( iphase )

      Conditional_Composition_MassFraction: if ( have_option( trim( option_path ) // &
           "/prognostic/initial_condition::WholeMesh/constant" ) ) then

         call get_option( trim( option_path ) // &
              "/prognostic/initial_condition::WholeMesh/constant", initial_constant )
         field_prot = initial_constant

      elseif( have_option( trim( option_path ) // &
           "/prognostic/initial_condition::WholeMesh/python") ) then

         call get_option( trim( option_path ) // &
              "/prognostic/initial_condition::WholeMesh/python", func )

         call allocate( dummy, field % mesh, "dummy" )
         call get_option( "/timestepping/current_time", current_time )
         call set_from_python_function( dummy, trim( func ), positions, current_time )
         field_prot = dummy % val
         call deallocate( dummy )

      else
         ewrite(-1,*) "No initial condition for field::", trim( field_name )
         FLAbort( " Check initial conditions " )

      end if Conditional_Composition_MassFraction


      Conditional_Composition_BC: if ( have_option( trim( option_path ) // &
           "/prognostic/boundary_conditions[0]/type::dirichlet" )) then

         bc_type = 1
         nobcs = get_boundary_condition_count( field )

         Loop_Over_BC: do k = 1, nobcs
            field_prot_bc => extract_surface_field( field, k, "value" )
            shape_option = option_shape( trim( option_path ) // &
                 "/prognostic/boundary_conditions[" // &
                 int2str( k - 1 ) // "]/surface_ids" )
            allocate( sufid_bc( 1 : shape_option( 1 ) ) )

            call get_option( trim( option_path ) // &
                 "/prognostic/boundary_conditions[" // &
                 int2str( k - 1 ) // "]/surface_ids", sufid_bc )

            allocate( face_nodes( face_loc( field, 1 ) ) )
            sele = 1
            do j = 1, stotel
               if( any ( sufid_bc == pmesh % faces % boundary_ids( j ) ) ) then
                  wic_bc( j + ( iphase - 1 ) * stotel ) = bc_type
                  face_nodes = ele_nodes( field_prot_bc, sele )
                  do kk = 1, snloc
                    suf_bc( ( icomp - ( nphases + 1 ) ) * nphases * stotel * snloc + &
                       ( iphase - 1 ) * stotel * snloc + ( j - 1 ) * snloc + kk ) = &
                          field_prot_bc % val( face_nodes( 1 ) )
                  end do
                  sele = sele + 1
               end if
            end do

            deallocate( face_nodes )
            deallocate( sufid_bc )

         end do Loop_Over_BC ! End of BC loop

      end if Conditional_Composition_BC


      if ( have_source )  then
         do j = 1, node_count( field_source )
            field_prot_source( ( iphase - 1 ) * node_count( field_source ) + j ) = &
                 field_source % val( j )
         end do
      end if

      if ( have_absorption ) then
         do j = 1, node_count( field_absorption )
            field_prot_absorption( ( iphase - 1 ) * node_count( field_absorption ) + j ) = &
                 field_absorption % val( j )
         end do
      end if

      return
    end subroutine Get_CompositionFields_Outof_State

    subroutine Get_ScalarFields_Outof_State( state, iphase, field, &
         field_prot, wic_bc, suf_bc, field_prot_source, field_prot_absorption )
      implicit none
      type( state_type ), dimension( : ), intent( inout ) :: state
      integer, intent( in ) :: iphase
      type( scalar_field ), pointer :: field, field_prot_bc
      real, dimension( : ), intent( inout ) :: field_prot
      real, dimension( : ), intent( inout ), optional :: field_prot_source, field_prot_absorption
      integer, dimension( : ), intent( inout ) :: wic_bc
      real, dimension( : ), intent( inout ) :: suf_bc

      ! Local variables
      type( mesh_type ), pointer :: pmesh, cmesh
      type(scalar_field), pointer :: pressure, field_source, field_absorption
      type(scalar_field) :: dummy
      type(vector_field), pointer :: positions
      integer, dimension(:), allocatable :: sufid_bc
      character( len = option_path_len ) :: option_path, field_name
      integer :: stotel, nobcs, bc_type, i, j, k, kk, sele
      integer :: nstates, nphases, ncomps, snloc, stat
      integer :: shape_option(2)
      real :: initial_constant
      logical :: have_source, have_absorption
      integer, dimension(:), allocatable :: face_nodes
      character( len = 8192 ) :: func

      field_name = trim( field % name )

      field_source => extract_scalar_field( state( iphase ), field_name // "Source", stat )
      have_source = (stat == 0)

      field_absorption => extract_scalar_field( state( iphase ), field_name // "Absorption", stat )
      have_absorption = (stat == 0)

      pressure => extract_scalar_field(state(1), "Pressure")
      ! this is basically p_snloc
      snloc = face_loc( pressure, 1 )

      pmesh => extract_mesh( state, "PressureMesh" )
      cmesh => extract_mesh( state, "CoordinateMesh" )

      positions => extract_vector_field(state(1), "Coordinate")

      stotel = surface_element_count( cmesh )
      nstates = option_count("/material_phase")
      ncomps = 0
      do i = 1, nstates
         if (have_option("/material_phase[" // int2str(i-1) // "]/is_multiphase_component")) then
            ncomps=ncomps+1
         end if
      end do
      nphases = nstates - ncomps

      option_path = "/material_phase["//int2str(iphase-1)//"]/scalar_field::"//trim(field_name)

      if ( have_option( trim(option_path)//"/prognostic/initial_condition::WholeMesh/constant") ) then

         call get_option(trim(option_path)//"/prognostic/initial_condition::WholeMesh/constant", &
              initial_constant)
         field_prot = initial_constant

      else if ( have_option( trim(option_path)//"/prognostic/initial_condition::WholeMesh/python") ) then


         call get_option( trim(option_path)//"/prognostic/initial_condition::WholeMesh/python", func )

         call allocate( dummy, field%mesh, "dummy" )


         call get_option("/timestepping/current_time", current_time)
         call set_from_python_function(dummy, trim(func), positions, current_time)
         field_prot = dummy%val

         call deallocate(dummy)

      else

         ewrite(-1, *) "No initial condition for field::", trim(field_name)
         FLAbort("Check initial conditions")

      end if

      Conditional_Field_BC: if( have_option( trim( option_path ) // &
           "/prognostic/boundary_conditions[0]/type::dirichlet" )) then

         BC_Type = 1
         nobcs = get_boundary_condition_count( field )

         do k = 1, nobcs
            field_prot_bc => extract_surface_field( field, k, "value" )

            shape_option = option_shape( trim( option_path ) // &
                 "/prognostic/boundary_conditions["// int2str(k-1)//"]/surface_ids" )

            allocate( sufid_bc( 1 : shape_option( 1 ) ) )

            call get_option( '/material_phase[' // int2str(iphase-1) // &
                 ']/scalar_field::' // trim(field_name) // '/prognostic/' // &
                 'boundary_conditions['//int2str(k-1)//']/surface_ids', SufID_BC )

            allocate( face_nodes( face_loc(field,1) ) )

            sele=1
            do j = 1, stotel
               if( any ( SufID_BC == pmesh%faces%boundary_ids(j) ) ) then

                  wic_bc( j + ( iphase - 1 ) * stotel ) = BC_Type

                  face_nodes = ele_nodes(field_prot_bc, sele)
                  do kk = 1, snloc
                     suf_bc( ( iphase - 1 ) * stotel * snloc + ( j - 1 ) * snloc + kk ) = &
                          field_prot_bc%val(face_nodes( 1 ))
                  end do

                  sele=sele+1
               end if
            end do

            deallocate(face_nodes)
            deallocate( sufid_bc )

         end do ! End of BC loop

      end if Conditional_Field_BC

      if (have_source) then
         do j = 1, node_count( field_source )
            field_prot_source( ( iphase - 1 ) * node_count( field_source ) + j ) = &
                 field_source%val(j)
         end do
      end if

      if (have_absorption) then
         do j = 1, node_count( field_absorption )
            field_prot_absorption( ( iphase - 1 ) * node_count( field_absorption ) + j ) = &
                 field_absorption%val(j)
         end do
      end if

      return
    end subroutine Get_ScalarFields_Outof_State

    subroutine Get_VectorFields_Outof_State( state, iphase, field, &
         field_u_prot, field_v_prot, field_w_prot, field_nu_prot, field_nv_prot, field_nw_prot, &
         wic_bc, suf_u_bc, suf_v_bc, suf_w_bc, field_prot_source, field_prot_absorption, is_overlapping )
      implicit none
      type( state_type ), dimension( : ), intent( inout ) :: state
      integer, intent( in ) :: iphase
      logical, intent( in ) :: is_overlapping
      type(vector_field), pointer :: field, field_prot_bc
      real, dimension( : ), intent( inout ) :: field_u_prot, field_v_prot, field_w_prot, &
           field_nu_prot, field_nv_prot, field_nw_prot
      real, dimension( : ), intent( inout ), optional :: field_prot_source
      real, dimension( : , :, : ), intent( inout ), optional :: field_prot_absorption
      integer, dimension( : ), intent( inout ) :: wic_bc
      real, dimension( : ), intent( inout ) :: suf_u_bc, suf_v_bc, suf_w_bc

      ! Local variables
      type( mesh_type ), pointer :: pmesh, cmesh
      type(vector_field) :: dummy
      type(vector_field), pointer :: positions
      type(scalar_field), pointer :: pressure, field_source, field_absorption
      integer, dimension(:), allocatable :: sufid_bc, face_nodes
      character( len = option_path_len ) :: option_path, field_name
      integer :: ndim, nlev, stotel, snloc, snloc2, nonods, nobcs, bc_type, j, k, kk, l, &
           shape_option( 2 ), count, nloc, nloc2, ele, ilev
      real, dimension( : ), allocatable :: initial_constant
      logical :: have_source, have_absorption
      character(len=8192) :: func

      pmesh => extract_mesh(state, "PressureMesh" )
      cmesh => extract_mesh(state, "CoordinateMesh" )

      positions => extract_vector_field(state(1), "Coordinate")
      ndim = field % dim

      nlev = 1
      if ( is_overlapping ) nlev = ele_loc(pmesh, 1)

      nonods = node_count( field ) * nlev

      stotel = surface_element_count( cmesh )

      snloc2 = face_loc( field, 1 )
      snloc = snloc2 * nlev
      nloc2 = ele_loc( field, 1 )
      nloc = nloc2 * nlev

      field_name = trim( field % name )

      option_path = "/material_phase["//int2str(iphase-1)//"]/vector_field::"//trim(field_name)
      if ( have_option( trim(option_path)//"/prognostic/initial_condition::WholeMesh/constant") ) then

         allocate(initial_constant(ndim)) ; initial_constant=0.
         call get_option(trim(option_path)//"/prognostic/initial_condition::WholeMesh/constant", &
              initial_constant)

         field_u_prot( (iphase-1) * nonods + 1 : iphase * nonods ) = initial_constant(1)
         if (ndim>1) field_v_prot( (iphase-1) * nonods + 1 : iphase * nonods ) = initial_constant(2)
         if (ndim>2) field_w_prot( (iphase-1) * nonods + 1 : iphase * nonods ) = initial_constant(3)

         deallocate( initial_constant )

      else if ( have_option( trim(option_path)//"/prognostic/initial_condition::WholeMesh/python") ) then

         call get_option( trim(option_path)//"/prognostic/initial_condition::WholeMesh/python", func )

         call allocate( dummy, field%dim, field%mesh, "dummy" )
         call get_option("/timestepping/current_time", current_time)
         call set_from_python_function(dummy, trim(func), positions, current_time)

         do ele = 1, ele_count( field )
            do ilev = 1, nlev
               field_u_prot( (iphase-1) * nonods + (ele-1) * nloc + (ilev-1) * nloc2 + 1 : &
                    &           (iphase-1) * nonods + (ele-1) * nloc + (ilev-1) * nloc2 + nloc2 ) = &
                    &    ele_val( dummy, 1, ele )

               if (ndim>1) &
                    field_v_prot( (iphase-1) * nonods + (ele-1) * nloc + (ilev-1) * nloc2 + 1 : &
                    &                (iphase-1) * nonods + (ele-1) * nloc + (ilev-1) * nloc2 + nloc2 ) = &
                    &         ele_val( dummy, 2, ele )

               if (ndim>2) &
                    field_v_prot( (iphase-1) * nonods + (ele-1) * nloc + (ilev-1) * nloc2 + 1 : &
                    &                (iphase-1) * nonods + (ele-1) * nloc + (ilev-1) * nloc2 + nloc2 ) = &
                    &          ele_val( dummy, 3, ele )
            end do
         end do

         call deallocate(dummy)

      else

         ewrite(-1, *) "No initial condition for field::", trim(field_name)
         FLAbort("Check initial conditions")

      end if

      ! set nu to u
      field_nu_prot((iphase-1) * nonods+1 : iphase * nonods) = &
           &      field_u_prot((iphase-1) * nonods+1 : iphase * nonods)            
      if (ndim>1) field_nv_prot((iphase-1) * nonods+1 : iphase * nonods) = &
           &      field_v_prot((iphase-1) * nonods+1 : iphase * nonods)
      if (ndim>2) field_nw_prot((iphase-1) * nonods+1 : iphase * nonods) = &
           &      field_w_prot((iphase-1) * nonods+1 : iphase * nonods)

      Conditional_Field_BC: if( have_option( '/material_phase[' // int2str(iphase-1) // &
           ']/vector_field::' // trim(field_name) // ' /prognostic/' // &
           'boundary_conditions[0]/type::dirichlet' )) then

         BC_Type = 1
         nobcs = get_boundary_condition_count( field )

         do k = 1, nobcs
            field_prot_bc => extract_surface_field( field, k, "value" )

            shape_option = option_shape('/material_phase[' // int2str(iphase-1) // &
                 ']/vector_field::' // trim(field_name) //&
                 '/prognostic/boundary_conditions[' //int2str(k-1)// ']/surface_ids' )

            allocate(sufid_bc(1:shape_option(1)))

            call get_option( '/material_phase[' // int2str(iphase-1) // &
                 ']/vector_field::' // trim(field_name) // &
                 '/prognostic/boundary_conditions[' // int2str(k-1) // &
                 ']/surface_ids', SufID_BC )

            allocate(face_nodes(face_loc(field, 1)))
            face_nodes=(/( l, l=1, snloc2)/)

            do j = 1, stotel

               if( any ( sufid_bc == field%mesh%faces%boundary_ids(j) ) ) then 

                  wic_bc( j  + ( iphase - 1 ) * stotel ) = BC_Type

                  count = 1
                  do kk = 1, snloc
                     suf_u_bc( ( iphase - 1 ) * stotel * snloc + (j-1) * snloc + kk ) = &
                          field_prot_bc%val( 1, face_nodes(count) )
                     if (ndim>1) suf_v_bc( ( iphase - 1 ) * stotel * snloc + & 
                          (j-1) * snloc + kk ) = field_prot_bc%val( 2, face_nodes(count) )
                     if (ndim>2) suf_w_bc( ( iphase - 1 ) * stotel * snloc + &
                          (j-1) * snloc + kk ) = field_prot_bc%val( 3, face_nodes(count) )

                     count = count + 1
                     if ( mod(kk, snloc2)==0. ) count=1
                  end do

                  face_nodes = face_nodes + snloc2

               end if
            end do

            deallocate(face_nodes)
            deallocate(sufid_bc)

         end do ! End of BC Loop 

      endif Conditional_Field_BC

      return
    end subroutine Get_VectorFields_Outof_State


    subroutine Get_Ele_Type(x_nloc, cv_ele_type, p_ele_type, u_ele_type)

      !-
      !- u_ele_type = cv_ele_type = p_ele_type will flag the dimension and 
      !- type of element:
      !- = 1 or 2: 1D (linear and quadratic, respectively)
      !- = 3 or 4: triangle (linear or quadratic, respectively)
      !- = 5 or 6: quadrilateral (bi-linear or tri-linear, respectively)
      !- = 7 or 8: tetrahedron (linear or quadratic, respectively)
      !- = 9 or 10: hexahedron (bi-linear or tri-linear, respectively)
      !-

      implicit none

      integer, intent(in) :: x_nloc
      integer, intent( inout ) :: cv_ele_type, p_ele_type, u_ele_type

      ! Local variables
      integer :: ndim, degree

      call get_option("/geometry/dimension", ndim)

      call get_option( &
           "/geometry/mesh::PressureMesh/from_mesh/mesh_shape/polynomial_degree", &
           degree)

      Select Case( ndim )
      case( 1 ) ! ndim

         Select Case( degree )

         case( 1 ) ! degree

            ! ndim=1; p=1
            cv_ele_type = 1

         case( 2 ) ! degree

            ! ndim=1; p=2
            cv_ele_type = 2

         case default; FLAbort("Degree error")

         end Select ! degree

      case( 2 ) ! ndim

         Select Case( degree )

         case( 1 ) ! degree

            Select Case( x_nloc )

            case( 3 ) ! x_nloc

               ! ndim=2; p=1; x_nloc=3
               ! linear triangle
               cv_ele_type = 3

            case( 4 ) ! x_nloc

               ! ndim=2; p=1; x_nloc=4
               ! bilinear quad
               cv_ele_type = 5

            case default; FLAbort("X_nloc error")

            end Select ! x_nloc

         case( 2 ) ! degree

            Select Case( x_nloc )

            case( 6 ) ! x_nloc

               ! ndim=2; p=2; x_nloc=3
               ! quadratic triangle
               cv_ele_type = 4

            case( 10 ) ! x_nloc

               ! ndim=2; p=2; x_nloc=4
               ! bi-quadratic quad
               cv_ele_type = 6

            case default; FLAbort("X_nloc error")

            end Select ! x_nloc

         case default; FLAbort("Degree error")

         end Select ! degree

      case( 3 ) ! ndim

         Select Case( degree )

         case( 1 ) ! degree

            Select Case( x_nloc )

            case( 4 ) ! x_nloc

               ! ndim=3; p=1; x_nloc=4
               ! linear tets
               cv_ele_type = 7

            case( 8 ) ! x_nloc

               ! ndim=3; p=1; x_nloc=8
               ! tri-linear hex
               cv_ele_type = 9

            case default; FLAbort("X_nloc error")

            end Select ! x_nloc

         case( 2 ) ! degree

            Select Case( x_nloc )

            case( 10 ) ! x_nloc

               ! ndim=3; p=2; x_nloc=4
               ! quadratic tet
               cv_ele_type = 8

            case( 27 ) ! x_nloc

               ! ndim=3; p=2; x_nloc=8
               ! bilinear quad
               cv_ele_type = 10

            case default; FLAbort("X_nloc error")

            end Select ! x_nloc

         case default; FLAbort("Degree error")

         end Select ! degree

      end Select ! ndim

      p_ele_type = cv_ele_type ; u_ele_type = cv_ele_type

      return
    end subroutine Get_Ele_Type


    subroutine Get_Solvers_Options( nphases, &
         nits, nits_internal, &
         ndpset, nits_flux_lim_volfra, nits_flux_lim_comp, &
         nits_flux_lim_t, t_disopt, v_disopt, t_beta, v_beta, &
         t_theta, v_theta, u_theta, &
         nonlinear_iteration_tolerance, lump_eqns  )
      ! This subroutine extract all solvers information from diamons
      implicit none
      integer, intent( in ) :: nphases
      integer, intent( inout ) :: nits, nits_internal, ndpset, &
           nits_flux_lim_volfra, nits_flux_lim_comp, &
           nits_flux_lim_t, t_disopt, v_disopt
      real, intent( inout ) :: t_beta, v_beta, t_theta, v_theta, u_theta, &
           nonlinear_iteration_tolerance
      logical, intent( inout ) :: lump_eqns

      call get_option( '/timestepping/nonlinear_iterations', nits, &
           default = 3 )

      if( have_option( '/timestepping/nonlinear_iterations/tolerance' )) then
         call get_option( '/timestepping/nonlinear_iterations/tolerance', &
              nonlinear_iteration_tolerance )
      else
         nonlinear_iteration_tolerance = 1.e-6
      end if

      ! This one is only for compositional problems
      call get_option('/material_phase[0]/scalar_field::component1/' // &
           'prognostic/temporal_discretisation/control_volumes/' // &
           'number_advection_iterations', nits_internal, default=1)

      ! Reference pressure node to be set
      call get_option('/material_phase[0]/scalar_field::Pressure/' // &
           'prognostic/spatial_discretisation/reference_node', ndpset, default = 0 )


      ! Get the number of advection CV face value iterations.
      call get_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/' // &
           'prognostic/temporal_discretisation/control_volumes/' // &
           'number_advection_iterations', nits_flux_lim_volfra, default=3 )

      call get_option('/material_phase[' // int2str(nphases) // &
           ']/scalar_field::ComponentMassFractionPhase1/' // &
           'temporal_discretisation/control_volumes/' // &
           'number_advection_iterations', nits_flux_lim_comp, default=3 )

      call get_option('/material_phase[0]/scalar_field::Temperature/prognostic/' // &
           'temporal_discretisation/control_volumes/number_advection_iterations', &
           nits_flux_lim_t, default=3 )

      ! Still hard-wired: MUST BE REVIEWED LATER
      t_disopt = 1

      if ( have_option('/material_phase[0]/scalar_field::Temperature/prognostic/' // &
           'spatial_discretisation/control_volumes/face_value::FiniteElement/' // &
           'limit_face_value/limiter::Extrema' ) ) then
         t_disopt = 9

      else
         if ( have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/' // &
              'spatial_discretisation/control_volumes/face_value::FiniteElement/' // &
              'limit_face_value' ) ) t_disopt = 5
      end if

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
      !       =8      Finite elements in space    Theta=specified    DOWNWIND+INTERFACE TRACKING
      !       =9      Finite elements in space    Theta=non-linear   DOWNWIND+INTERFACE TRACKING

      ! Add this variables into spud and replace them in cv_assemb and multi_dyncore
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
         if( have_option( '/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
              'spatial_discretisation/control_volumes/face_value::FiniteElement/' // &
              'limit_face_value/limiter::Sweby' ) ) then
            v_disopt = 5 !! Unless all the other options, but need to be able to get 8 here
         elseif( have_option( '/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
              'spatial_discretisation/control_volumes/face_value::FiniteElement/' // &
              'limit_face_value/limiter::Extrema' ) ) then
            v_disopt = 9
         else
            v_disopt = 8
         end if
      endif

      call get_option('/material_phase[0]/scalar_field::Temperature/prognostic/' // &
           'spatial_discretisation/conservative_advection', t_beta, default=0.0)
      call get_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/conservative_advection', v_beta )
      call get_option('/material_phase[0]/scalar_field::Temperature/prognostic/' // &
           'temporal_discretisation/theta', t_theta, default=1.)

      call get_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'temporal_discretisation/theta', v_theta)

      call get_option('/material_phase[0]/vector_field::Velocity/prognostic/' // &
           'temporal_discretisation/theta', u_theta)

      !! I'm going to get this one from the PhaseVolumeFraction scalar_field
      lump_eqns = have_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/' // &
           'prognostic/spatial_discretisation/continuous_galerkin/' // &
           'mass_terms/lump_mass_matrix')

      return
    end subroutine Get_Solvers_Options

    subroutine xp1_2_xp2( ndim, totele, x_nloc_p2, x_nloc_p1, x_nonods_p1, x_nonods_p2, &
         x_ndgln_p1, x_ndgln_p2, positions, &
         x, y, z )
      ! This subrt maps the coordinate P1 mesh into a P2 mesh. 
      implicit none
      integer, intent( in ) :: ndim, totele, x_nloc_p2, x_nloc_p1, x_nonods_p1, x_nonods_p2
      integer, dimension( totele * x_nloc_p1 ), intent( in ) :: x_ndgln_p1
      integer, dimension( totele * x_nloc_p2 ), intent( in ) :: x_ndgln_p2
      real, dimension( x_nonods_p2 ), intent( inout ) :: x, y, z
      type(vector_field), intent(in) :: positions

      ! Local variables
      real, dimension( x_nonods_p1 ) :: x_p1, y_p1, z_p1
      integer, dimension( : ), allocatable :: iloclist_p2
      real, dimension( : ), allocatable :: x2, y2, z2
      integer :: ele, iloc, inod
      real :: xnod1, xnod2, ynod1, ynod2, xtemp, ytemp

      allocate( iloclist_p2 ( x_nloc_p2 ) ) ; iloclist_p2 = 0
      allocate( x2( x_nloc_p2 ) ) ; x2 = 0.
      allocate( y2( x_nloc_p2 ) ) ; y2 = 0.
      allocate( z2( x_nloc_p2 ) ) ; z2 = 0.

      if( ndim == 2 ) then
         iloclist_p2 = (/ 1, 4, 2, 5, 6, 3 /) 
      elseif( ndim == 3 ) then
         iloclist_p2 = (/ 1, 5, 2, 6, 7, 3, 8, 9, 10, 4 /)
      else
         iloclist_p2 = (/ 1, 3, 2 /)           
         !         FLAbort( "Still need to be done" )
      end if

      Conditional_Pn: if ( x_nloc_p2 == 3 .or. x_nloc_p2 == 4 ) then

         if( ( x_nloc_p2 == 3 ) .and. ( ndim == 1 ) ) then ! 1D quadratic
            x_p1 = positions % val( 1, : )
            do ele = 1, totele
               x2 = 0.
               do iloc = 1, x_nloc_p1
                  x2( iloc ) = x_p1( x_ndgln_p1( ( ele - 1 ) * x_nloc_p1 + iloc )) 
               end do
               do iloc = 1, x_nloc_p1 - 1
                  xnod1 = x2( iloc )
                  xnod2 = x2( iloc + 1 )
                  x2( x_nloc_p1 + iloc ) = 0.5 * (  xnod1  +  xnod2  )              
               end do
               do iloc = 1, x_nloc_p2
                  inod = x_ndgln_p2( ( ele - 1 ) * x_nloc_p2 + iloc )
                  x( inod ) = x2( iloclist_p2( iloc ) )
               end do
            end do

         else

            x = positions % val( 1, : )
            if( ndim > 1 ) y = positions % val( 2, : )
            if( ndim > 2 ) z = positions % val( 3, : )
         end if

      else if ( x_nloc_p2 == 6 .or. x_nloc_p2 == 10 ) then

         x_p1 = positions % val( 1, : )
         y_p1 = positions % val( 2, : )
         if (ndim == 3)  z_p1 = positions % val( 3, : )

         if ( ( x_nloc_p2 == 6 ) .and. ( ndim == 2 ) ) then ! 2D P2 Tri

            do ele = 1, totele

               x2 = 0. ; y2 = 0.
               do iloc = 1, x_nloc_p1
                  x2( iloc ) = x_p1( x_ndgln_p1( ( ele - 1 ) * x_nloc_p1 + iloc )) 
                  y2( iloc ) = y_p1( x_ndgln_p1( ( ele - 1 ) * x_nloc_p1 + iloc ))
               end do

               do iloc = 1, x_nloc_p1
                  if( iloc < x_nloc_p1 ) then
                     xnod1 = x2( iloc )      ; ynod1 = y2( iloc )
                     xnod2 = x2( iloc + 1 ); ynod2 = y2( iloc + 1 )
                  else
                     xnod1 = x2( iloc ) ; ynod1 = y2( iloc )
                     xnod2 = x2( 1 )    ; ynod2 = y2 ( 1 )
                  end if
                  x2( x_nloc_p1 + iloc ) = 0.5 * (  xnod1  +  xnod2  )
                  y2( x_nloc_p1 + iloc ) = 0.5 * (  ynod1  +  ynod2  )
               end do

               xtemp = x2( 5 ) ; ytemp = y2( 5 )
               x2( 5 ) = x2( 6 ) ; y2( 5 ) = y2( 6 )
               x2( 6 ) = xtemp ; y2( 6 ) = ytemp

               do iloc = 1, x_nloc_p2
                  inod = x_ndgln_p2( ( ele - 1 ) * x_nloc_p2 + iloc )
                  x( inod ) = x2( iloclist_p2( iloc ) )
                  y( inod ) = y2( iloclist_p2( iloc ) )
               end do

            end do

         else ! Quadratic Tets

            do ele = 1, totele

               x2 = 0. ; y2 = 0. ; z2 = 0.
               do iloc = 1, x_nloc_p1
                  x2( iloc ) = x_p1( x_ndgln_p1( ( ele - 1 ) * x_nloc_p1 + iloc )) 
                  y2( iloc ) = y_p1( x_ndgln_p1( ( ele - 1 ) * x_nloc_p1 + iloc ))
                  z2( iloc ) = z_p1( x_ndgln_p1( ( ele - 1 ) * x_nloc_p1 + iloc ))
               end do

               x2( 5 ) = 0.5 * (x2(1) + x2(2) )
               y2( 5 ) = 0.5 * (y2(1) + y2(2) )
               z2( 5 ) = 0.5 * (z2(1) + z2(2) )

               x2( 6 ) = 0.5 * (x2(1) + x2(3) )
               y2( 6 ) = 0.5 * (y2(1) + y2(3) )
               z2( 6 ) = 0.5 * (z2(1) + z2(3) )

               x2( 7 ) = 0.5 * (x2(2) + x2(3) )
               y2( 7 ) = 0.5 * (y2(2) + y2(3) )
               z2( 7 ) = 0.5 * (z2(2) + z2(3) )

               x2( 8 ) = 0.5 * (x2(1) + x2(4) )
               y2( 8 ) = 0.5 * (y2(1) + y2(4) )
               z2( 8 ) = 0.5 * (z2(1) + z2(4) )

               x2( 9 ) = 0.5 * (x2(2) + x2(4) )
               y2( 9 ) = 0.5 * (y2(2) + y2(4) )
               z2( 9 ) = 0.5 * (z2(2) + z2(4) )

               x2( 10 ) = 0.5 * (x2(3) + x2(4) )
               y2( 10 ) = 0.5 * (y2(3) + y2(4) )
               z2( 10 ) = 0.5 * (z2(3) + z2(4) )


               do iloc = 1, x_nloc_p2
                  inod = x_ndgln_p2( ( ele - 1 ) * x_nloc_p2 + iloc )
                  x( inod ) = x2( iloclist_p2( iloc ) )
                  y( inod ) = y2( iloclist_p2( iloc ) )
                  z( inod ) = z2( iloclist_p2( iloc ) )
               end do

            end do

         end if
      end if Conditional_Pn


      deallocate( iloclist_p2 )
      deallocate( x2 )
      deallocate( y2 )
      deallocate( z2 )


      return
    end subroutine xp1_2_xp2

  end module copy_outof_into_state
