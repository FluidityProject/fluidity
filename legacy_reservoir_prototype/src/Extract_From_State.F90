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

  module Copy_Outof_State
    !! This module enables the multiphase prototype code to interact with state by
    !! copying everything required from state to the MP-space.

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
    use futils, only: int2str
    !use printout
    !use quicksort


    implicit none

    private

    public :: Get_Primary_Scalars, Compute_Node_Global_Numbers, Extracting_MeshDependentFields_From_State, &
         Get_ScalarFields_Outof_State, Get_CompositionFields_Outof_State, Get_VectorFields_Outof_State, &
         Extract_Position_Field, xp1_2_xp2, Get_Ele_Type, Get_Discretisation_Options

    interface Get_Ndgln
       module procedure Get_Scalar_Ndgln, Get_Vector_Ndgln, Get_Mesh_Ndgln
    end interface Get_Ndgln

    interface Get_SNdgln
       module procedure Get_Scalar_SNdgln, Get_Vector_SNdgln
    end interface Get_SNdgln

  contains


    subroutine Get_Primary_Scalars( state, &         
         nphase, nstate, ncomp, totele, ndim, stotel, &
         u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
         x_snloc, cv_snloc, u_snloc, p_snloc, &
         cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, dx, &
         is_overlapping )
!!$ This subroutine extracts all primary variables associated with the mesh from state,
!!$ and associated them with the variables used in the MultiFluids model.
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      integer, intent( inout ) :: nphase, nstate, ncomp, totele, ndim, stotel
      integer, intent( inout ), optional :: u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, &
           mat_nloc, x_snloc, cv_snloc, u_snloc, p_snloc, cv_nonods, mat_nonods, u_nonods, &
           xu_nonods, x_nonods, x_nonods_p1, p_nonods
      real, intent( inout ), optional :: dx
      logical, intent( inout ), optional :: is_overlapping

!!$ Local variables
      character( len = option_path_len ) :: vel_element_type
      type( vector_field ), pointer :: positions, velocity
      type( scalar_field ), pointer :: pressure
      type( mesh_type ), pointer :: velocity_cg_mesh, pressure_cg_mesh
      integer :: i, degree

      ewrite(3,*)' In Get_Primary_Scalars'

!!$ Defining dimension and nstate
      call get_option( '/geometry/dimension', ndim )
      nstate = option_count( '/material_phase' )

!!$ Assume there are the same number of components in each phase (will need to check this eventually)
      ncomp = 0
      do i = 1, nstate
         if( have_option( '/material_phase[' // int2str(i-1) // &
              ']/is_multiphase_component' ) ) then
            ncomp = ncomp + 1
         end if
      end do
      nphase = nstate - ncomp
      assert( nphase > 0 ) ! Check if there is more than 0 phases

!!$ Get the vel element type.
      if( present( is_overlapping ) ) then
         call get_option('/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
              vel_element_type )
         is_overlapping = .false.
         if ( trim( vel_element_type ) == 'overlapping' ) is_overlapping = .true. 
      end if

      positions => extract_vector_field( state, 'Coordinate' )
      pressure_cg_mesh => extract_mesh( state, 'PressureMesh_Continuous' )

!!$ Defining number of elements and surface elements, coordinates, locs and snlocs
      totele = ele_count( positions )
      stotel = surface_element_count( positions )

!!$ Coordinates
      if( present( x_nloc_p1 ) ) x_nloc_p1 = ele_loc( positions, 1 )
      if( present( x_nloc ) ) x_nloc = ele_loc( pressure_cg_mesh, 1 )
      if( present( x_snloc ) ) x_snloc = face_loc( pressure_cg_mesh, 1 )
      if( present( x_nonods_p1 ) ) x_nonods_p1 = node_count( positions )
      if( present( x_nonods ) ) x_nonods = node_count( pressure_cg_mesh )

!!$ Pressure, Control Volumes and Materials
      pressure => extract_scalar_field( state, 'Pressure' )
      if( present( p_nloc ) ) p_nloc = ele_loc( pressure, 1 )
      if( present( p_snloc ) ) p_snloc = face_loc( pressure, 1 )
      if( present( p_nonods ) ) p_nonods = node_count( pressure )
      if( present( cv_nloc ) ) cv_nloc = p_nloc
      if( present( cv_snloc ) ) cv_snloc = p_snloc
      if( present( cv_nonods ) ) cv_nonods = p_nonods
      if( present( mat_nloc ) ) mat_nloc = cv_nloc
      if( present( mat_nonods ) ) mat_nonods = mat_nloc * totele

!!$ Velocities and velocities (DG) associated with the continuous space (CG)
      velocity => extract_vector_field( state, 'Velocity' )
      if( present( u_nloc ) ) u_nloc = ele_loc( velocity, 1 )
      if( present( u_snloc ) ) u_snloc = face_loc( velocity, 1 )
      if( present( u_nonods ) ) u_nonods = node_count( velocity )

!!$ Get the continuous space of the velocity field
      velocity_cg_mesh => extract_mesh( state, 'VelocityMesh_Continuous' )
      if( present( xu_nloc ) ) xu_nloc = ele_loc( velocity_cg_mesh, 1 )
      if( present( xu_nonods ) ) xu_nonods = max(( xu_nloc - 1 ) * totele + 1, totele )

!!$ Take care of overlapping elements
      if( present( is_overlapping ) ) then
         if( ( is_overlapping ) .and. ( ndim > 1 ) ) u_nonods = u_nonods * cv_nloc
         if( ( is_overlapping ) .and. ( ndim > 1 ) ) u_nloc = u_nloc * cv_nloc
         if( is_overlapping ) u_snloc = u_snloc * cv_nloc 
      end if

!!$ Used just for 1D:
      if( present( dx ) ) dx = maxval( positions % val( 1, : ) ) - minval( positions % val( 1, : ) )

      return
    end subroutine Get_Primary_Scalars


    subroutine Compute_Node_Global_Numbers( state, &
         is_overlapping, totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
         cv_snloc, p_snloc, u_snloc, &
         cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln, &
         cv_sndgln, p_sndgln, u_sndgln )
!!$ This subroutine calculates the global node numbers requested to operates in the MP-space.
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      type( vector_field ), pointer :: positions, velocity
      type( mesh_type ), pointer :: pressure_cg_mesh, velocity_cg_mesh
      type( scalar_field ), pointer :: pressure
      logical, intent( in ) :: is_overlapping
      integer, intent( in ) :: totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
           cv_snloc, p_snloc, u_snloc
      integer, dimension( : ) :: cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln, &
           cv_sndgln, p_sndgln, u_sndgln
!!$ Local variables 
      integer, dimension( : ), allocatable :: u_sndgln2
      integer :: u_nloc2, u_snloc2, sele, u_siloc2, ilev, u_siloc, ele, inod_remain

!!$ Linear mesh coordinate
      positions => extract_vector_field( state( 1 ), 'Coordinate' )
      call Get_Ndgln( x_ndgln_p1, positions )

!!$ Positions/Coordinates
      pressure_cg_mesh => extract_mesh( state( 1 ), 'PressureMesh_Continuous' )
      call Get_Ndgln( x_ndgln, pressure_cg_mesh )

!!$ Pressure, control volume and material
      pressure => extract_scalar_field( state( 1 ), 'Pressure' )
      call Get_Ndgln( cv_ndgln, pressure )
      p_ndgln = cv_ndgln ; mat_ndgln = cv_ndgln

!!$ Velocities
      velocity => extract_vector_field( state( 1 ), 'Velocity' )
      call Get_Ndgln( u_ndgln, velocity, is_overlapping, cv_nloc )

!!$ Velocity in the continuous space
      velocity_cg_mesh => extract_mesh( state( 1 ), 'VelocityMesh_Continuous' )
      call Get_Ndgln( xu_ndgln, velocity_cg_mesh )

!!$ Surface-based global node numbers for control volumes and pressure
      call Get_SNdgln( cv_sndgln, pressure )
      p_sndgln = cv_sndgln

!!$ Velocities
      if ( is_overlapping ) then
         u_snloc2 = u_snloc / cv_nloc
         u_nloc2 = u_nloc / cv_nloc
      else
         u_snloc2 = u_snloc
         u_nloc2 = u_nloc
      end if

      allocate( u_sndgln2( stotel * u_snloc2 ) ) ; u_sndgln2 = 0
      call Get_SNdgln( u_sndgln2, velocity )

      if ( is_overlapping ) then ! Convert u_sndgln2 to overlapping u_sndgln
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
         u_sndgln = u_sndgln2
      end if

      deallocate( u_sndgln2 )

      return
    end subroutine Compute_Node_Global_Numbers


    subroutine Get_Ele_Type( x_nloc, cv_ele_type, p_ele_type, u_ele_type, &
         mat_ele_type, u_sele_type, cv_sele_type )
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
      integer, intent( inout ), optional :: mat_ele_type, u_sele_type, cv_sele_type
!!$ Local variables
      integer :: ndim, degree

      call get_option( '/geometry/dimension', ndim) 

      call get_option( &
           '/geometry/mesh::PressureMesh/from_mesh/mesh_shape/polynomial_degree', &
           degree )

      Select Case( ndim )
      case( 1 ) ! ndim

         Select Case( degree )

         case( 1 ) ! degree

!!$ ndim=1; p=1
            cv_ele_type = 1

         case( 2 ) ! degree

            ! ndim=1; p=2
            cv_ele_type = 2

         case default; FLAbort('Degree error')

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

            case default; FLAbort('X_nloc error')

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

            case default; FLAbort('X_nloc error')

            end Select ! x_nloc

         case default; FLAbort('Degree error')

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

            case default; FLAbort('X_nloc error')

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

            case default; FLAbort('X_nloc error')

            end Select ! x_nloc

         case default; FLAbort('Degree error')

         end Select ! degree

      end Select ! ndim

      p_ele_type = cv_ele_type ; u_ele_type = cv_ele_type

!!$ The following options are hardcoded and need to be either deleted from the code tree or
!!$ added into the schema.
      if( present( mat_ele_type ) ) mat_ele_type = 1 
      if( present( u_sele_type ) ) u_sele_type = 1 
      if( present( cv_sele_type ) ) cv_sele_type = 1

      return
    end subroutine Get_Ele_Type



    subroutine Get_Discretisation_Options( state, &
         is_overlapping, &
         t_disopt, v_disopt, t_beta, v_beta, t_theta, v_theta, u_theta, &
         t_dg_vel_int_opt, u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
         comp_diffusion_opt, ncomp_diff_coef, in_ele_upwind, dg_ele_upwind, &
         nits_flux_lim_t, nits_flux_lim_volfra, nits_flux_lim_comp, &
         volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux, &
         t_use_theta_flux, t_get_theta_flux, scale_momentum_by_volume_fraction )
!!$ This subroutine extract all discretisation options from the schema
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      logical, intent( in ) :: is_overlapping
      integer, intent( inout ) :: t_disopt, v_disopt
      real, intent( inout ) :: t_beta, v_beta, t_theta, v_theta, u_theta
      integer, intent( inout ) :: t_dg_vel_int_opt, u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
           comp_diffusion_opt, ncomp_diff_coef, in_ele_upwind, dg_ele_upwind, &
           nits_flux_lim_t, nits_flux_lim_volfra, nits_flux_lim_comp
      logical, intent( inout ) :: volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, &
           comp_get_theta_flux, t_use_theta_flux, t_get_theta_flux, scale_momentum_by_volume_fraction

!!$ Local variables:
      integer :: nphase, nstate, ncomp, totele, ndim, stotel, iphase
      character( len = option_path_len ) :: option_path, option_path2, option_path3

!!$ DISOPT Options:
!!$ =0      1st order in space          Theta=specified    UNIVERSAL
!!$ =1      1st order in space          Theta=non-linear   UNIVERSAL
!!$ =2      Trapezoidal rule in space   Theta=specified    UNIVERSAL
!!$ =2      if isotropic limiter then FEM-quadratic & stratification adjust. Theta=non-linear 
!!$ =3      Trapezoidal rule in space   Theta=non-linear   UNIVERSAL
!!$ =4      Finite elements in space    Theta=specified    UNIVERSAL
!!$ =5      Finite elements in space    Theta=non-linear   UNIVERSAL
!!$ =6      Finite elements in space    Theta=specified    NONE
!!$ =7      Finite elements in space    Theta=non-linear   NONE
!!$ =8      Finite elements in space    Theta=specified    DOWNWIND+INTERFACE TRACKING
!!$ =9      Finite elements in space    Theta=non-linear   DOWNWIND+INTERFACE TRACKING

!!$ Extracting primary scalars as local variables:
      call Get_Primary_Scalars( state, &         
           nphase, nstate, ncomp, totele, ndim, stotel )

!!$ Solving Advection Field: Temperature
      option_path = '/material_phase[0]/scalar_field::Temperature'
      option_path2 = trim( option_path ) //  '/prognostic/spatial_discretisation'
      option_path3 = trim( option_path ) //  '/prognostic/temporal_discretisation/control_volumes/number_advection_iterations'
      t_disopt = 1

      call get_option( trim( option_path3 ), nits_flux_lim_t, default = 3 )

      Conditional_TDISOPT: if( have_option( trim( option_path2 ) ) ) then
         if( have_option( trim( option_path2 ) // '/control_volumes/face_value::FiniteElement/limit_face_value/' // &
              'limiter::Extrema' ) ) then
            t_disopt = 9
         else
            if( have_option( trim( option_path2 ) // '/control_volumes/face_value::FiniteElement/limit_face_value' ) ) &
                 t_disopt = 5
         end if
      end if Conditional_TDISOPT

      call get_option( trim( option_path2 ) // '/conservative_advection', t_beta, default = 0.0 )
      call get_option( '/material_phase[0]/scalar_field::Temperature/prognostic/temporal_discretisation/theta', &
           t_theta, default = 1. )

!!$ Solving Advection Field: Volume fraction
      option_path = '/material_phase[0]/scalar_field::PhaseVolumeFraction'
      option_path2 = trim( option_path ) // '/prognostic/spatial_discretisation/control_volumes/face_value'
      option_path3 = trim( option_path ) // '/prognostic/temporal_discretisation/control_volumes/number_advection_iterations'
      v_disopt = 8

      call get_option( trim( option_path3 ), nits_flux_lim_volfra, default = 3 )

      Conditional_VDISOPT: if( have_option( trim( option_path ) ) ) then
         if( have_option( trim( option_path2 ) // '::FirstOrderUpwind' ) ) v_disopt = 0
         if( have_option( trim( option_path2 ) // '::Trapezoidal' ) ) v_disopt = 2
         if( have_option( trim( option_path2 ) // 'FiniteElement/do_not_limit_face_value' ) ) v_disopt = 6
         if( have_option( trim( option_path2 ) // 'FiniteElement/limit_face_value/limiter::Sweby' ) ) v_disopt = 5
         if( have_option( trim( option_path2 ) // 'FiniteElement/limit_face_value/limiter::Extrema' ) ) v_disopt = 9
      end if Conditional_VDISOPT

      call get_option( trim( option_path ) // '/prognostic/spatial_discretisation/conservative_advection', v_beta )
      call get_option( trim( option_path ) // '/prognostic/temporal_discretisation/theta', v_theta )

!!$ Solving Velocity Field
      call get_option( '/material_phase[0]/vector_field::Velocity/prognostic/temporal_discretisation/theta', u_theta )

!!$ Solving Component Field
      option_path3 = '/material_phase[' // int2str( nphase ) // ']/scalar_field::ComponentMassFractionPhase1/' // &
           'temporal_discretisation/control_volumes/number_advection_iterations'
      call get_option( trim( option_path3 ), nits_flux_lim_comp, default = 3 )

!!$ Scaling factor for the momentum equation
      scale_momentum_by_volume_fraction = .false.
      do iphase = 1, nphase
         option_path = '/material_phase[' // int2str( iphase - 1 ) // ']/scale_momentum_by_volume_fraction'
         if( have_option( trim( option_path ) ) ) scale_momentum_by_volume_fraction = .true.
      end do

!!$ Options below are hardcoded and need to be added into the schema
      t_dg_vel_int_opt = 1 ; u_dg_vel_int_opt = 4 ; v_dg_vel_int_opt = 4 ; w_dg_vel_int_opt = 0
      if( .not. is_overlapping ) v_dg_vel_int_opt = 1
      comp_diffusion_opt = 0 ; ncomp_diff_coef = 0
      volfra_use_theta_flux = .false. ; volfra_get_theta_flux = .true.
      comp_use_theta_flux = .false.   ; comp_get_theta_flux = .true.
      t_use_theta_flux = .false.      ; t_get_theta_flux = .true.

!!$ IN/DG_ELE_UPWIND are options for optimisation of upwinding across faces in the overlapping
!!$ formulation. The data structure and options for this formulation need to be added later. 
      in_ele_upwind =  3 ; dg_ele_upwind = 3

      return
    end subroutine Get_Discretisation_Options



    subroutine Extract_Position_Field( state, &
         xu, yu, zu )
!!$ This subroutine extracts the spatial coordinates fields from state-space and copy them into
!!$ MP-space
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      type( mesh_type ), pointer :: velocity_cg_mesh
      type( vector_field ), pointer :: positions
      type( vector_field ) :: velocity_cg
      real, dimension( : ) :: xu, yu, zu
!!$ Local variables
      integer :: ndim, totele, xu_nloc, xu_nonods, ele, iloc
      integer, dimension( : ), pointer :: element_nodes

      call get_option( '/geometry/dimension', ndim )
      positions => extract_vector_field( state, 'Coordinate' )
      totele = ele_count( positions )

!!$ Velocity in the continuous space
      velocity_cg_mesh => extract_mesh( state( 1 ), 'VelocityMesh_Continuous' )
      call allocate( velocity_cg, ndim, velocity_cg_mesh, 'Velocity_CG_Coordinates' )
      velocity_cg % val( :, : )= 0
      call project_field( positions, velocity_cg, positions )
      xu_nloc = ele_loc( velocity_cg_mesh, 1 )
      xu_nonods = max( ( xu_nloc - 1 ) * totele + 1, totele )

      Loop_Elements: do ele = 1, totele
         element_nodes => ele_nodes( velocity_cg_mesh, ele )
         Loop_Local_Nodes: do iloc = 1, xu_nloc
            xu( element_nodes( iloc ) ) = velocity_cg % val( 1, element_nodes( iloc ) )
            if( ndim > 1 )  yu( element_nodes( iloc ) ) = velocity_cg % val( 2, element_nodes( iloc ) )
            if( ndim > 2 )  zu( element_nodes( iloc ) ) = velocity_cg % val( 3, element_nodes( iloc ) )
         end do Loop_Local_Nodes
      end do Loop_Elements

      call deallocate( velocity_cg )

      return
    end subroutine Extract_Position_Field


    subroutine xp1_2_xp2( state, &
         x_nloc_p2, x_nloc_p1, x_nonods_p1, x_nonods_p2, &
         x_ndgln_p1, x_ndgln_p2, &
         x, y, z )
      ! This subrt maps the coordinate P1 mesh into a P2 mesh. 
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      type( vector_field ), pointer :: positions
      integer, intent( in ) :: x_nloc_p2, x_nloc_p1, x_nonods_p1, x_nonods_p2
      integer, dimension( : ), intent( in ) :: x_ndgln_p1
      integer, dimension( : ), intent( in ) :: x_ndgln_p2
      real, dimension( : ), intent( inout ) :: x, y, z

      ! Local variables
      real, dimension( x_nonods_p1 ) :: x_p1, y_p1, z_p1
      integer, dimension( : ), allocatable :: iloclist_p2
      real, dimension( : ), allocatable :: x2, y2, z2
      integer :: ndim, totele, ele, iloc, inod
      real :: xnod1, xnod2, ynod1, ynod2, xtemp, ytemp

      call get_option( '/geometry/dimension', ndim )
      positions => extract_vector_field( state( 1 ), 'Coordinate' )
      totele = ele_count( positions )

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

      deallocate( iloclist_p2, x2, y2, z2 )

      return
    end subroutine xp1_2_xp2


    subroutine Extracting_MeshDependentFields_From_State( state, &
         xu, yu, zu, x, y, z, &
         PhaseVolumeFraction, PhaseVolumeFraction_BC_Spatial, PhaseVolumeFraction_BC, PhaseVolumeFraction_Source, &
         Pressure_CV, Pressure_FEM, Pressure_FEM_BC_Spatial, Pressure_FEM_BC, &
         Density, Density_BC_Spatial, Density_BC, &
         Component, Component_BC_Spatial, Component_BC, Component_Source, &
         Velocity_U, Velocity_V, Velocity_W, Velocity_NU, Velocity_NV, Velocity_NW, &
         Velocity_U_BC_Spatial, Velocity_U_BC, Velocity_V_BC, Velocity_W_BC, Velocity_U_Source, Velocity_Absorption, &
         Temperature, Temperature_BC_Spatial, Temperature_BC, Temperature_Source, &
         Porosity, Permeability  )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      integer, dimension( : ), intent( inout ) :: PhaseVolumeFraction_BC_Spatial, Pressure_FEM_BC_Spatial, &
           Density_BC_Spatial, Component_BC_Spatial, Velocity_U_BC_Spatial, Temperature_BC_Spatial
      real, dimension( : ), intent( inout ) :: xu, yu, zu, x, y, z, &
           PhaseVolumeFraction, PhaseVolumeFraction_BC, PhaseVolumeFraction_Source, &
           Pressure_CV, Pressure_FEM, Pressure_FEM_BC, &
           Density, Density_BC, &
           Component, Component_BC, Component_Source, &
           Velocity_U, Velocity_V, Velocity_W, Velocity_NU, Velocity_NV, Velocity_NW, &
           Velocity_U_BC, Velocity_V_BC, Velocity_W_BC, Velocity_U_Source, &
           Temperature, Temperature_BC, Temperature_Source, &
           Porosity
      real, dimension( :, :, : ), intent( inout ) :: Velocity_Absorption, Permeability

!!$ Local variables
      type( scalar_field ), pointer :: field
      type( vector_field ), pointer :: vectorfield, grav_vectorfield
      type( tensor_field ), pointer :: tensorfield
      integer, dimension( : ), pointer :: element_nodes
      character( len = option_path_len ) :: option_path
      integer :: nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, &
           cv_ele_type, p_ele_type, u_ele_type, &
           stat, istate, iphase, jphase, icomp, ele, idim, jdim, knod, knod2
      integer, dimension( : ), allocatable :: cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, &
           xu_ndgln, mat_ndgln, cv_sndgln, p_sndgln, u_sndgln
      real :: dx
      logical :: is_overlapping, is_isotropic, is_diagonal, is_symmetric, have_gravity
      real, dimension( :, : ), allocatable :: constant

!!$ Extracting spatial resolution
      call Get_Primary_Scalars( state, &         
           nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, dx, &
           is_overlapping )

!!$ Calculating Global Node Numbers
      allocate( x_ndgln_p1( totele * x_nloc_p1 ), x_ndgln( totele * x_nloc ), cv_ndgln( totele * cv_nloc ), &
           p_ndgln( totele * p_nloc ), mat_ndgln( totele * mat_nloc ), u_ndgln( totele * u_nloc ), &
           xu_ndgln( totele * xu_nloc ), cv_sndgln( stotel * cv_snloc ), p_sndgln( stotel * p_snloc ), &
           u_sndgln( stotel * u_snloc ) )

      x_ndgln_p1 = 0 ; x_ndgln = 0 ; cv_ndgln = 0 ; p_ndgln = 0 ; mat_ndgln = 0 ; u_ndgln = 0 ; &
           xu_ndgln = 0 ; cv_sndgln = 0 ; p_sndgln = 0 ; u_sndgln = 0

      call Compute_Node_Global_Numbers( state, &
           is_overlapping, totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
           cv_snloc, p_snloc, u_snloc, &
           cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln, &
           cv_sndgln, p_sndgln, u_sndgln )

      call Get_Ele_Type( x_nloc, cv_ele_type, p_ele_type, u_ele_type )

      xu = 0. ; yu = 0. ; zu = 0.
      call Extract_Position_Field( state, &
           xu, yu, zu )

      x = 0. ; y = 0. ; z = 0.
      call xp1_2_xp2( state, &
           x_nloc, x_nloc_p1, x_nonods_p1, x_nonods, &
           x_ndgln_p1, x_ndgln, &
           x, y, z )

!!$
!!$ Extracting Volume Fraction (or Saturation) Field:
!!$
      Loop_VolumeFraction: do iphase = 1, nphase
         field => extract_scalar_field( state( iphase ), 'PhaseVolumeFraction' )
         knod = ( iphase - 1 ) * node_count( field )
         call Get_ScalarFields_Outof_State( state, iphase, field, &
              PhaseVolumeFraction( knod + 1 : knod + node_count( field ) ), PhaseVolumeFraction_BC_Spatial, &
              PhaseVolumeFraction_BC, &
              field_prot_source = PhaseVolumeFraction_Source )
      end do Loop_VolumeFraction

!!$
!!$ Extracting Pressure Field:
!!$
      field => extract_scalar_field( state( 1 ), 'Pressure' )
      call Get_ScalarFields_Outof_State( state, 1, field, &
           Pressure_FEM, Pressure_FEM_BC_Spatial, Pressure_FEM_BC )
      Pressure_CV = Pressure_FEM
      if( nphase > 1 ) then ! Copy this to the other phases
         do iphase = 2, nphase
            Pressure_FEM_BC_Spatial( ( iphase - 1 ) * stotel + 1 : iphase * stotel ) = &
                 Pressure_FEM_BC_Spatial( 1 : stotel )
            Pressure_FEM_BC( ( iphase - 1 ) * stotel * p_snloc + 1 : iphase * stotel * p_snloc ) = &
                 Pressure_FEM_BC( 1 : stotel * p_snloc )
         end do
      end if

!!$
!!$ Extracting Density Field:
!!$
      Loop_Density: do iphase = 1, nphase
         field => extract_scalar_field( state( iphase ), 'Density' )
         knod = ( iphase - 1 ) * node_count( field )
         call Get_ScalarFields_Outof_State( state, iphase, field, &
              Density(  knod + 1 : knod + node_count( field ) ), Density_BC_Spatial, Density_BC )
      end do Loop_Density

!!$
!!$ Extracting Components Field:
!!$
      Loop_Components: do icomp = nphase + 1, nphase + ncomp ! Component loop
         Loop_Phases_Components: do iphase = 1, nphase ! Phase loop

            field => extract_scalar_field( state( icomp ), 'ComponentMassFractionPhase' // int2str( iphase ) )

            knod = ( icomp - ( nphase + 1 ) ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods
            knod2 = ( icomp - ( nphase + 1 ) ) * nphase * stotel * cv_snloc + &
                 ( iphase - 1 ) * stotel * cv_snloc

            call Get_CompositionFields_Outof_State( state, nphase, icomp, iphase, field, &
                 Component( knod + 1 : knod + cv_nonods ), Component_BC_Spatial, &
                 knod2 + 1, knod2 + stotel * cv_snloc, &
                 Component_BC( knod2 + 1 : knod2 + stotel * cv_snloc ), &
                 field_prot_source = Component_Source( ( iphase - 1 ) * cv_nonods + 1 : &
                 ( iphase - 1 ) * cv_nonods + cv_nonods ) )

         end do Loop_Phases_Components
      end do Loop_Components

!!$
!!$ Extracting Velocity Field:
!!$
      Loop_Velocity: do iphase = 1, nphase
         vectorfield => extract_vector_field( state( iphase ), 'Velocity' )
         call Get_VectorFields_Outof_State( state, iphase, vectorfield, &
              Velocity_U, Velocity_V, Velocity_W, Velocity_NU, Velocity_NV, Velocity_NW, &
              Velocity_U_BC_Spatial, Velocity_U_BC, Velocity_V_BC, Velocity_W_BC, &
              Velocity_U_Source, Velocity_Absorption, is_overlapping )
      end do Loop_Velocity

!!$
!!$ Extracting Temperature Field:
!!$
      do istate = 1, nstate
         Conditional_Temperature: if( have_option( '/material_phase[' // int2str( istate - 1 ) // &
              ']/scalar_field::Temperature' ) ) then
            field => extract_scalar_field( state( istate ), 'Temperature' )
            knod = ( istate - 1 ) * node_count( field )
            call Get_ScalarFields_Outof_State( state, istate, field, &
                 Temperature( knod + 1 : knod + node_count( field ) ), &
                 Temperature_BC_Spatial, Temperature_BC, &
                 field_prot_source = Temperature_Source( knod + 1 : knod + node_count( field ) ) )          
         end if Conditional_Temperature
      end do


!!$
!!$ Extracting Porosity field (assuming for now that porosity is constant
!!$  across an element):
!!$
      field => extract_scalar_field( state, 'Porosity' )
      do ele = 1, element_count( field )
         element_nodes => ele_nodes( field, ele )
         Porosity( ele ) = field % val ( element_nodes( 1 ) )
      end do

!!$
!!$ Extracting Permeability Field:
!!$
      Conditional_PermeabilityField: if( have_option( '/porous_media/scalar_field::Permeability' ) ) then
         field => extract_scalar_field( state( 1 ), 'Permeability' )
         do ele = 1, element_count( field ) 
            element_nodes => ele_nodes( field, ele )
            forall( idim = 1 : ndim ) Permeability( ele, idim, idim ) = field % val( element_nodes( 1 ) )
         end do

      elseif( have_option( '/porous_media/tensor_field::Permeability' ) ) then
         option_path = '/porous_media/tensor_field::Permeability'
         if( have_option( trim( option_path ) // '/prescribed' ) ) then
            option_path = trim( option_path ) // '/prescribed/value[0]'
         else
            FLAbort( 'Permeability field not defined - sort this out' )
         endif
         is_isotropic = have_option( trim( option_path ) // '/isotropic' )
         is_diagonal = have_option( trim( option_path ) // '/diagonal' )
         is_symmetric = have_option( trim( option_path ) // '/anisotropic_symmetric' )

         Conditional_TensorPermeability: if( is_isotropic )then
            option_path = trim( option_path ) // '/isotropic'

            if( have_option( trim( option_path ) // '/constant')) then
               allocate( constant( 1, 1 ) ) ; constant = 0.
               call get_option( trim( option_path ) // '/constant', constant( 1, 1 ) )
               do idim = 1, ndim
                  Permeability( 1 : totele, idim, idim ) = constant( 1, 1 )
               end do
               deallocate( constant )

            elseif( have_option( trim( option_path ) // '/python' ) ) then
               tensorfield => extract_tensor_field( state( 1 ), 'Permeability' )
               do ele = 1, element_count( tensorfield )
                  element_nodes => ele_nodes( tensorfield, ele )
                  do idim = 1, ndim
                     Permeability( ele, idim, idim ) = &
                          tensorfield % val(idim, idim, element_nodes( 1 ) )
                  end do
               end do

            else
               FLExit( 'Incorrect initial condition for field' )
            end if

         elseif( is_diagonal )then ! Conditional_TensorPermeability
            option_path = trim( option_path ) // '/diagonal'

            if( have_option( trim( option_path ) // '/constant' ) )then
               allocate(constant( ndim, 1 ) ) ; constant = 0.
               call get_option( trim( option_path ) // '/constant', constant )
               do idim = 1, ndim
                  Permeability( 1 : totele, idim, idim ) = constant( idim, 1 )
               end do
               deallocate( constant )

            elseif( have_option( trim( option_path ) // '/python' ) )then
               tensorfield => extract_tensor_field( state( 1 ), 'Permeability' )
               do idim = 1, ndim
                  do ele = 1, element_count( tensorfield )
                     Permeability( ele, idim, idim ) = tensorfield % val( idim, idim,  ele )
                  end do
               end do

            else
               FLExit( 'Incorrect initial condition for field' )
            end if

         else
            if( is_symmetric ) then
               option_path = trim( option_path ) // '/anisotropic_symmetric'
            else
               option_path = trim( option_path ) // '/anisotropic_asymmetric'
            end if

            Conditional_Perm_Anisotropic: if( have_option( trim( option_path ) // '/constant' ) ) then
               allocate( constant( ndim, ndim ) ) ; constant = 0.
               call get_option( trim( option_path ) // '/constant', constant )
               do idim = 1, ndim
                  do jdim = 1, ndim
                     Permeability( 1 : totele, idim, jdim ) = constant( idim, jdim )
                  end do
               end do
               deallocate( constant )

            elseif( have_option( trim( option_path ) // '/python' ) ) then
               tensorfield => extract_tensor_field( state( 1 ), 'Permeability' )
               do ele = 1, element_count( tensorfield )
                  element_nodes => ele_nodes( tensorfield, ele )
                  do idim = 1, ndim
                     do jdim = 1, ndim
                        Permeability( ele, idim, jdim ) = &
                             tensorfield % val( idim, jdim, element_nodes( 1 ) )
                     end do
                  end do
               end do

            else 
               FLExit( 'Incorrect initial condition for field ')
            end if Conditional_Perm_Anisotropic

         end if Conditional_TensorPermeability

      end if Conditional_PermeabilityField


      deallocate( cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, &
           xu_ndgln, mat_ndgln, cv_sndgln, p_sndgln, u_sndgln )

      return
    end subroutine Extracting_MeshDependentFields_From_State




    subroutine Get_ScalarFields_Outof_State( state, iphase, field, &
         field_prot, wic_bc, suf_bc, field_prot_source, field_prot_absorption )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
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
      integer :: nstate, nphase, ncomp, snloc, stat
      integer :: shape_option(2)
      real :: initial_constant
      logical :: have_source, have_absorption
      integer, dimension(:), allocatable :: face_nodes
      character( len = 8192 ) :: func

      field_name = trim( field % name )

      field_source => extract_scalar_field( state( iphase ), field_name // 'Source', stat )
      have_source = (stat == 0)

      field_absorption => extract_scalar_field( state( iphase ), field_name // 'Absorption', stat )
      have_absorption = (stat == 0)

      pressure => extract_scalar_field(state(1), 'Pressure')
      ! this is basically p_snloc
      snloc = face_loc( pressure, 1 )

      pmesh => extract_mesh( state, 'PressureMesh' )
      cmesh => extract_mesh( state, 'CoordinateMesh' )

      positions => extract_vector_field(state(1), 'Coordinate')

      stotel = surface_element_count( cmesh )
      nstate = option_count('/material_phase')
      ncomp = 0
      do i = 1, nstate
         if (have_option('/material_phase[' // int2str(i-1) // ']/is_multiphase_component')) then
            ncomp=ncomp+1
         end if
      end do
      nphase = nstate - ncomp

      option_path = '/material_phase['//int2str(iphase-1)//']/scalar_field::'//trim(field_name)

      if ( have_option( trim(option_path)//'/prognostic/initial_condition::WholeMesh/constant') ) then

         call get_option(trim(option_path)//'/prognostic/initial_condition::WholeMesh/constant', &
              initial_constant)
         field_prot = initial_constant

      else if ( have_option( trim(option_path)//'/prognostic/initial_condition::WholeMesh/python') ) then


         call get_option( trim(option_path)//'/prognostic/initial_condition::WholeMesh/python', func )

         call allocate( dummy, field%mesh, 'dummy' )


         call get_option('/timestepping/current_time', current_time)
         call set_from_python_function(dummy, trim(func), positions, current_time)
         field_prot = dummy%val

         call deallocate(dummy)

      else

         ewrite(-1, *) 'No initial condition for field::', trim(field_name)
         FLAbort('Check initial conditions')

      end if

      Conditional_Field_BC: if( have_option( trim( option_path ) // &
           '/prognostic/boundary_conditions[0]/type::dirichlet' )) then

         BC_Type = 1
         nobcs = get_boundary_condition_count( field )

         do k = 1, nobcs
            field_prot_bc => extract_surface_field( field, k, 'value' )

            shape_option = option_shape( trim( option_path ) // &
                 '/prognostic/boundary_conditions['// int2str(k-1)//']/surface_ids' )

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


    subroutine Get_CompositionFields_Outof_State( state, nphase, icomp, iphase, field, &
         field_prot, wic_bc, &
         kprime, kprime2, &
         suf_bc, &
         field_prot_source, field_prot_absorption )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      integer, intent( in ) :: nphase, icomp, iphase
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
      integer :: nstate, stotel, nobcs, bc_type, i, j, k, kk, sele, stat, snloc
      real :: initial_constant
      character( len = 8192 ) :: func

      field_name = trim( field % name )
      positions => extract_vector_field( state( 1 ), 'Coordinate' )
      pressure => extract_scalar_field( state( 1 ), 'Pressure' )
      pmesh => extract_mesh( state, 'PressureMesh' )
      cmesh => extract_mesh( state, 'CoordinateMesh' )

      field_source => extract_scalar_field( state( iphase ), field_name // 'Source', stat )
      have_source = ( stat == 0 )
      field_absorption => extract_scalar_field( state( iphase ), field_name // 'Absorption', stat )
      have_absorption = ( stat == 0 )

      snloc = face_loc( pressure, 1 )
      nstate = option_count('/material_phase')
      stotel = surface_element_count( cmesh )

      option_path = '/material_phase[' // int2str( icomp - 1 ) // &
           ']/scalar_field::ComponentMassFractionPhase' // &
           int2str( iphase )

      Conditional_Composition_MassFraction: if ( have_option( trim( option_path ) // &
           '/prognostic/initial_condition::WholeMesh/constant' ) ) then

         call get_option( trim( option_path ) // &
              '/prognostic/initial_condition::WholeMesh/constant', initial_constant )
         field_prot = initial_constant

      elseif( have_option( trim( option_path ) // &
           '/prognostic/initial_condition::WholeMesh/python') ) then

         call get_option( trim( option_path ) // &
              '/prognostic/initial_condition::WholeMesh/python', func )

         call allocate( dummy, field % mesh, 'dummy' )
         call get_option( '/timestepping/current_time', current_time )
         call set_from_python_function( dummy, trim( func ), positions, current_time )
         field_prot = dummy % val
         call deallocate( dummy )

      else
         ewrite(-1,*) 'No initial condition for field::', trim( field_name )
         FLAbort( ' Check initial conditions ' )

      end if Conditional_Composition_MassFraction


      Conditional_Composition_BC: if ( have_option( trim( option_path ) // &
           '/prognostic/boundary_conditions[0]/type::dirichlet' )) then

         bc_type = 1
         nobcs = get_boundary_condition_count( field )

         Loop_Over_BC: do k = 1, nobcs
            field_prot_bc => extract_surface_field( field, k, 'value' )
            shape_option = option_shape( trim( option_path ) // &
                 '/prognostic/boundary_conditions[' // &
                 int2str( k - 1 ) // ']/surface_ids' )
            allocate( sufid_bc( 1 : shape_option( 1 ) ) )

            call get_option( trim( option_path ) // &
                 '/prognostic/boundary_conditions[' // &
                 int2str( k - 1 ) // ']/surface_ids', sufid_bc )

            allocate( face_nodes( face_loc( field, 1 ) ) )
            sele = 1
            do j = 1, stotel
               if( any ( sufid_bc == pmesh % faces % boundary_ids( j ) ) ) then
                  wic_bc( j + ( iphase - 1 ) * stotel ) = bc_type
                  face_nodes = ele_nodes( field_prot_bc, sele )
                  do kk = 1, snloc
                     suf_bc( ( icomp - ( nphase + 1 ) ) * nphase * stotel * snloc + &
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



    subroutine Get_VectorFields_Outof_State( state, iphase, field, &
         field_u_prot, field_v_prot, field_w_prot, field_nu_prot, field_nv_prot, field_nw_prot, &
         wic_bc, suf_u_bc, suf_v_bc, suf_w_bc, field_prot_source, field_prot_absorption, is_overlapping )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
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
      integer :: ndim, stotel, snloc, snloc2, nonods, nobcs, bc_type, j, k, kk, l, &
           shape_option( 2 ), count
      real, dimension( : ), allocatable :: initial_constant
      logical :: have_source, have_absorption
      character(len=8192) :: func

      pmesh => extract_mesh(state, 'PressureMesh' )
      cmesh => extract_mesh(state, 'CoordinateMesh' )

      positions => extract_vector_field(state(1), 'Coordinate')

      ndim = field % dim
      stotel = surface_element_count( cmesh )
      snloc2 = face_loc(field, 1)
      snloc = snloc2
      if ( is_overlapping ) snloc = snloc2 * ele_loc(pmesh, 1)

      nonods = node_count( field )
      field_name = trim( field % name )

      option_path = '/material_phase['//int2str(iphase-1)//']/vector_field::'//trim(field_name)
      if ( have_option( trim(option_path)//'/prognostic/initial_condition::WholeMesh/constant') ) then

         allocate(initial_constant(ndim)) ; initial_constant=0.
         call get_option(trim(option_path)//'/prognostic/initial_condition::WholeMesh/constant', &
              initial_constant)

         field_u_prot((iphase-1) * nonods+1 : iphase * nonods)=initial_constant(1)
         if (ndim>1) field_v_prot((iphase-1) * nonods+1 : iphase * nonods)=initial_constant(2)
         if (ndim>2) field_w_prot((iphase-1) * nonods+1 : iphase * nonods)=initial_constant(3)

         deallocate( initial_constant )

      else if ( have_option( trim(option_path)//'/prognostic/initial_condition::WholeMesh/python') ) then

         call get_option( trim(option_path)//'/prognostic/initial_condition::WholeMesh/python', func )

         call allocate( dummy, field%dim, field%mesh, 'dummy' )
         call get_option('/timestepping/current_time', current_time)
         call set_from_python_function(dummy, trim(func), positions, current_time)

         field_u_prot((iphase-1) * nonods+1 : iphase * nonods)=dummy%val(1,:)
         if (ndim>1) field_v_prot((iphase-1) * nonods+1 : iphase * nonods)=dummy%val(2,:)
         if (ndim>2) field_w_prot((iphase-1) * nonods+1 : iphase * nonods)=dummy%val(3,:)

         call deallocate(dummy)

      else

         ewrite(-1, *) 'No initial condition for field::', trim(field_name)
         FLAbort('Check initial conditions')

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
            field_prot_bc => extract_surface_field( field, k, 'value' )

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





!!$
!!$ Module Get_Ndgln Interfaces

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
!!$            ewrite(3,*)'ele, iloc, ndgln:', ele, iloc, &
!!$                 ndgln( ( ele - 1 ) * ele_loc( field, ele ) + iloc )
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


      call get_option( '/geometry/dimension', ndim )
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
!!$            ewrite(3,*)'ele, iloc, ndgln:', ele, iloc, &
!!$                 ndgln( ( ele - 1 ) * ele_loc( field, 1 ) * cv_nloc2  + iloc )
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
!!$            ewrite(3,*)'ele, iloc, ndgln:', ele, iloc, &
!!$                 ndgln( ( ele - 1 ) * ele_loc( mesh, ele ) + iloc )
         end do
      end do

      return
    end subroutine Get_Mesh_Ndgln

!!$
!!$ Module Get_SNdgln Interfaces

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
!!$            ewrite(3,*)'sele, iloc, sndgln:', sele, iloc, &
!!$                 sndgln( ( sele - 1 ) * face_loc( field, sele ) + iloc )
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
!!$            ewrite(3,*)'sele, iloc, sndgln:', sele, iloc, &
!!$                 sndgln( ( sele - 1 ) * face_loc( field, sele ) + iloc )
         end do
      end do

      deallocate( snloc )

      return
    end subroutine Get_Vector_SNdgln


  end module Copy_Outof_State

