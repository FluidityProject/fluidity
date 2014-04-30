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
    use global_parameters, only: option_path_len, is_overlapping, is_compact_overlapping
    use diagnostic_fields_wrapper_new
    use element_numbering
    use shape_functions
    use fefields
    use boundary_conditions
    use futils, only: int2str
    use boundary_conditions_from_options
 use memory_diagnostics
 use initialise_fields_module, only: initialise_field_over_regions

    !use printout
    !use quicksort


    implicit none

    private

    public :: Get_Primary_Scalars, Compute_Node_Global_Numbers, Extracting_MeshDependentFields_From_State, &
         Extract_TensorFields_Outof_State, Extract_Position_Field, xp1_2_xp2, Get_Ele_Type, Get_Discretisation_Options, &
         print_from_state, update_boundary_conditions, pack_multistate, finalise_multistate, get_ndglno, Adaptive_NonLinear,&
         get_var_from_packed_state, as_vector

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
         cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, dx )
!!$ This subroutine extracts all primary variables associated with the mesh from state,
!!$ and associated them with the variables used in the MultiFluids model.
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      integer, intent( inout ) :: nphase, nstate, ncomp, totele, ndim, stotel
      integer, intent( inout ), optional :: u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, &
           mat_nloc, x_snloc, cv_snloc, u_snloc, p_snloc, cv_nonods, mat_nonods, u_nonods, &
           xu_nonods, x_nonods, x_nonods_p1, p_nonods
      real, intent( inout ), optional :: dx

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
      call get_option('/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           vel_element_type )
      is_overlapping = .false.
      if ( trim( vel_element_type ) == 'overlapping' ) is_overlapping = .true.

      is_compact_overlapping = .false.
      if ( trim( vel_element_type ) == 'is_compact_overlapping' ) is_compact_overlapping = .true.

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
      if( present( u_nonods ) ) then
         if( ( is_overlapping ) .and. ( ndim > 1 ) ) u_nonods = u_nonods * cv_nloc
         if( ( is_overlapping ) .and. ( ndim > 1 ) ) u_nloc = u_nloc * cv_nloc
         if( is_overlapping ) u_snloc = u_snloc * cv_nloc 
      end if

!!$ Used just for 1D:
      if( present( dx ) ) dx = maxval( positions % val( 1, : ) ) - minval( positions % val( 1, : ) )

      return
    end subroutine Get_Primary_Scalars

    function get_ndglno(mesh) result(ndglno)
          type(mesh_type) :: mesh
          integer, dimension(:), pointer  ::  ndglno

          ndglno=> mesh%ndglno
        end function get_ndglno


    subroutine Compute_Node_Global_Numbers( state, &
         totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
         cv_snloc, p_snloc, u_snloc, &
         cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln, &
         cv_sndgln, p_sndgln, u_sndgln )
!!$ This subroutine calculates the global node numbers requested to operates in the MP-space.
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      type( vector_field ), pointer :: positions, velocity
      type( mesh_type ), pointer :: pressure_cg_mesh, velocity_cg_mesh
      type( scalar_field ), pointer :: pressure
      integer, intent( in ) :: totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
           cv_snloc, p_snloc, u_snloc
      integer, dimension( : ), pointer  :: cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln
      integer, dimension( : ) ::  cv_sndgln, p_sndgln, u_sndgln
!!$ Local variables 
      integer, dimension( : ), allocatable :: u_sndgln2
      integer :: u_nloc2, u_snloc2, sele, u_siloc2, ilev, u_siloc, ele, inod_remain, i

      x_ndgln_p1=>get_ndglno(extract_mesh(state(1),"CoordinateMesh"))
      x_ndgln=>get_ndglno(extract_mesh(state(1),"PressureMesh_Continuous"))
      cv_ndgln=>get_ndglno(extract_mesh(state(1),"PressureMesh"))
      p_ndgln=>get_ndglno(extract_mesh(state(1),"PressureMesh"))
      mat_ndgln=>get_ndglno(extract_mesh(state(1),"PressureMesh_Discontinuous"))
      u_ndgln=>get_ndglno(extract_mesh(state(1),"InternalVelocityMesh"))
      xu_ndgln=>get_ndglno(extract_mesh(state(1),"VelocityMesh_Continuous"))

!!$ Linear mesh coordinate
      positions => extract_vector_field( state( 1 ), 'Coordinate' )
!!$      call Get_Ndgln( x_ndgln_p1, positions )
!!$
!!$ Positions/Coordinates
      pressure_cg_mesh => extract_mesh( state( 1 ), 'PressureMesh_Continuous' )
!!$      call Get_Ndgln( x_ndgln, pressure_cg_mesh )
!!$
!!$ Pressure, control volume and material
      pressure => extract_scalar_field( state( 1 ), 'Pressure' )
!!$      call Get_Ndgln( cv_ndgln, pressure )
!!$      p_ndgln = cv_ndgln
!!$      mat_ndgln = (/ (i, i = 1, totele * cv_nloc ) /)
!!$
!!$ Velocities
      velocity => extract_vector_field( state( 1 ), 'Velocity' )
!!$      call Get_Ndgln( u_ndgln, velocity, cv_nloc )
!!$
!!$ Velocity in the continuous space
      velocity_cg_mesh => extract_mesh( state( 1 ), 'VelocityMesh_Continuous' )
!!$      call Get_Ndgln( xu_ndgln, velocity_cg_mesh )

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
         t_disopt, v_disopt, t_beta, v_beta, t_theta, v_theta, u_theta, &
         t_dg_vel_int_opt, u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
         comp_diffusion_opt, ncomp_diff_coef, in_ele_upwind, dg_ele_upwind, &
         nits_flux_lim_t, nits_flux_lim_volfra, nits_flux_lim_comp, &
         volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux, &
         t_use_theta_flux, t_get_theta_flux, scale_momentum_by_volume_fraction )
!!$ This subroutine extract all discretisation options from the schema
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
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
              'limiter::CompressiveAdvection' ) ) then
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
         if( have_option( trim( option_path2 ) // '::FiniteElement/do_not_limit_face_value' ) ) v_disopt = 6
         if( have_option( trim( option_path2 ) // '::FiniteElement/limit_face_value/limiter::Sweby' ) ) v_disopt = 5
         if( have_option( trim( option_path2 ) // '::FiniteElement/limit_face_value/limiter::CompressiveAdvection' ) ) v_disopt = 9
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
      comp_use_theta_flux = .false. ; comp_get_theta_flux = .true.
      t_use_theta_flux = .false. ; t_get_theta_flux = .true.

!!$ IN/DG_ELE_UPWIND are options for optimisation of upwinding across faces in the overlapping
!!$ formulation. The data structure and options for this formulation need to be added later. 
      in_ele_upwind = 3 ; dg_ele_upwind = 3

      ! simplify things...
      !in_ele_upwind = 1 ; dg_ele_upwind = 1
      !u_dg_vel_int_opt = 1
      !v_dg_vel_int_opt = 1
      !v_disopt = 0

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


    subroutine Extracting_MeshDependentFields_From_State( state, packed_state, initialised, &
         PhaseVolumeFraction, PhaseVolumeFraction_Source, &
         Component, Component_Source, &
         Velocity_U_Source, Velocity_Absorption, &
         Temperature, Temperature_Source,  &
         Permeability  )
      implicit none
      type( state_type ), dimension( : ), intent( inout ) :: state
      type( state_type ), intent( inout ) :: packed_state

      logical, intent( in ) :: initialised
      real, dimension( : ), intent( inout ) :: &
           PhaseVolumeFraction, PhaseVolumeFraction_Source, &
           Component, Component_Source, &
           Velocity_U_Source, &
           Temperature, Temperature_Source
      real, dimension( :, :, : ), intent( inout ) :: Velocity_Absorption, Permeability

!!$ Local variables
      type( scalar_field ), pointer :: scalarfield
      type( vector_field ), pointer :: vectorfield, grav_vectorfield, x_all
      type( tensor_field ), pointer :: tensorfield
      integer, dimension( : ), pointer :: element_nodes
      character( len = option_path_len ) :: option_path
      integer :: nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, &
           cv_ele_type, p_ele_type, u_ele_type, &
           stat, istate, iphase, jphase, icomp, ele, idim, jdim, knod, knod2
      integer, dimension( : ), pointer :: cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, &
           xu_ndgln, mat_ndgln, cv_sndgln, p_sndgln, u_sndgln
      real :: dx
      logical :: is_overlapping, is_isotropic, is_diagonal, is_symmetric, have_gravity
      real, dimension( :, : ), allocatable :: constant
      real, dimension( : ), allocatable :: x, y, z, dummy

!!$ Extracting spatial resolution
      call Get_Primary_Scalars( state, &         
           nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, dx )

!!$ Calculating Global Node Numbers
      allocate( cv_sndgln( stotel * cv_snloc ), p_sndgln( stotel * p_snloc ), &
           u_sndgln( stotel * u_snloc ), dummy( cv_nonods ) )

!      x_ndgln_p1 = 0 ; x_ndgln = 0 ; cv_ndgln = 0 ; p_ndgln = 0 ; mat_ndgln = 0
 !     u_ndgln = 0 ;  xu_ndgln = 0 ; 
      cv_sndgln = 0 ; p_sndgln = 0 ; u_sndgln = 0

      call Compute_Node_Global_Numbers( state, &
           totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
           cv_snloc, p_snloc, u_snloc, &
           cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln, &
           cv_sndgln, p_sndgln, u_sndgln )

      call Get_Ele_Type( x_nloc, cv_ele_type, p_ele_type, u_ele_type )



      allocate( x( x_nonods ) , y( x_nonods ), z( x_nonods ) )

      x = 0. ; y = 0. ; z = 0.
      call xp1_2_xp2( state, &
           x_nloc, x_nloc_p1, x_nonods_p1, x_nonods, &
           x_ndgln_p1, x_ndgln, &
           x, y, z )

      x_all => extract_vector_field( packed_state, "PressureCoordinate" )
      do idim = 1, ndim
          if( idim ==1 ) x_all % val( idim, : ) = x
          if( idim ==2 )   x_all % val( idim, : ) = y
          if( idim ==3 )   x_all % val( idim, : ) = z
      end do

!!$
!!$ Extracting Volume Fraction (or Saturation) Field:
!!$
      Loop_VolumeFraction: do iphase = 1, nphase
         scalarfield => extract_scalar_field( state( iphase ), 'PhaseVolumeFraction' )
         knod = ( iphase - 1 ) * node_count( scalarfield )
         call Get_ScalarFields_Outof_State( state, initialised, iphase, scalarfield, &
              PhaseVolumeFraction( knod + 1 : knod + node_count( scalarfield ) ), &
              field_prot_source = PhaseVolumeFraction_Source )
      end do Loop_VolumeFraction

!!$
!!$ Extracting Pressure Field:
!!$
      scalarfield => extract_scalar_field( state( 1 ), 'Pressure' )
      call Get_ScalarFields_Outof_State( state, initialised, 1, scalarfield, &
           dummy )

!!$
!!$ Extracting Density Field:
!!$
      Loop_Density: do iphase = 1, nphase
         scalarfield => extract_scalar_field( state( iphase ), 'Density' )
         !knod = ( iphase - 1 ) * node_count( scalarfield )
         call Get_ScalarFields_Outof_State( state, initialised, iphase, scalarfield, &
              dummy)
      end do Loop_Density

!!$
!!$ Extracting Components Field:
!!$
      Loop_Components: do icomp = nphase + 1, nphase + ncomp ! Component loop
         Loop_Phases_Components: do iphase = 1, nphase ! Phase loop

            scalarfield => extract_scalar_field( state( icomp ), 'ComponentMassFractionPhase' // int2str( iphase ) )

            knod = ( icomp - ( nphase + 1 ) ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods
            knod2 = ( icomp - ( nphase + 1 ) ) * nphase * stotel * cv_snloc + &
                 ( iphase - 1 ) * stotel * cv_snloc

            call Get_CompositionFields_Outof_State( state, initialised, nphase, icomp, iphase, scalarfield, &
                 Component( knod + 1 : knod + cv_nonods ), &
                 field_prot_source = Component_Source( ( iphase - 1 ) * cv_nonods + 1 : &
                 ( iphase - 1 ) * cv_nonods + cv_nonods ) )

         end do Loop_Phases_Components
      end do Loop_Components

!!$
!!$ Extracting Velocity Field:
!!$
      Loop_Velocity: do iphase = 1, nphase
         vectorfield => extract_vector_field( state( iphase ), 'Velocity' )
         call Get_VectorFields_Outof_State( state, initialised, iphase, vectorfield, &
              field_prot_source=Velocity_U_Source, field_prot_absorption=Velocity_Absorption )
      end do Loop_Velocity

!!$
!!$ Extracting Temperature Field:
!!$
      do iphase = 1, nphase
         Conditional_Temperature: if( have_option( '/material_phase[' // int2str( iphase - 1 ) // &
              ']/scalar_field::Temperature' ) ) then
            scalarfield => extract_scalar_field( state( iphase ), 'Temperature' )
            knod = ( iphase - 1 ) * node_count( scalarfield )
            call Get_ScalarFields_Outof_State( state, initialised, iphase, scalarfield, &
                 Temperature( knod + 1 : knod + node_count( scalarfield ) ), &
                 field_prot_source = Temperature_Source( knod + 1 : knod + node_count( scalarfield ) ) )
         end if Conditional_Temperature
      end do

!!$
!!$ Extracting Permeability Field:
!!$
      Permeability = 0.
      Conditional_PermeabilityField: if( have_option( '/porous_media/scalar_field::Permeability' ) ) then

         scalarfield => extract_scalar_field( state( 1 ), 'Permeability' )
         do ele = 1, element_count( scalarfield ) 
            element_nodes => ele_nodes( scalarfield, ele )
            forall( idim = 1 : ndim ) Permeability( ele, idim, idim ) = scalarfield % val( element_nodes( 1 ) )
         end do

      elseif( have_option( '/porous_media/tensor_field::Permeability' ) ) then

         tensorfield => extract_tensor_field( state( 1 ), 'Permeability' )
         option_path =  '/porous_media/tensor_field::Permeability'
         call Extract_TensorFields_Outof_State( state, 1, &
              tensorfield, option_path, &
              Permeability )

      elseif( have_option( '/porous_media/vector_field::Permeability' ) ) then
         FLAbort( 'Permeability Vector Field is not defined yet.' )

      end if Conditional_PermeabilityField

      deallocate( cv_sndgln, p_sndgln, u_sndgln, dummy )

      return
    end subroutine Extracting_MeshDependentFields_From_State




    subroutine Get_ScalarFields_Outof_State( state, initialised, iphase, field, &
         field_prot, wic_bc, suf_bc, field_prot_source, field_prot_absorption, suf_bc_rob1, suf_bc_rob2 )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      logical, intent( in ) :: initialised
      integer, intent( in ) :: iphase
      type( scalar_field ), pointer :: field, field_prot_bc, field_prot_bc1, field_prot_bc2
      real, dimension( : ), intent( inout ) :: field_prot
      real, dimension( : ), intent( inout ), optional :: field_prot_source, field_prot_absorption, suf_bc_rob1, suf_bc_rob2
      integer, dimension( : ), intent( inout ), optional :: wic_bc
      real, dimension( : ), intent( inout ), optional :: suf_bc

      ! Local variables
      type( mesh_type ), pointer :: pmesh, cmesh
      type(scalar_field), pointer :: pressure, field_source, field_absorption
      type(scalar_field) :: dummy
      type(vector_field), pointer :: positions
      integer, dimension(:), allocatable :: sufid_bc
      character( len = option_path_len ) :: option_path, option_path2, field_name
      integer :: stotel, nobcs, bc_type, i, j, k, kk, sele
      integer :: nstate, nphase, ncomp, snloc, stat
      integer :: shape_option(2)
      real :: initial_constant
      logical :: have_source, have_absorption
      integer, dimension(:), allocatable :: face_nodes
      character( len = 8192 ) :: func

      field_name = trim( field % name )
      have_source = .false. ; have_absorption = .false.

      Conditional_SourceField: if( present( field_prot_source ) ) then
         field_source => extract_scalar_field( state( iphase ), trim(field_name) // 'Source', stat )
         have_source = ( stat == 0 )

         if ( have_source ) then
            do j = 1, node_count( field_source )
               field_prot_source( ( iphase - 1 ) * node_count( field_source ) + j ) = &
                    field_source % val( j )
            end do
         end if
      end if Conditional_SourceField

      Conditional_AbsorptionField: if( present( field_prot_absorption ) ) then
         field_absorption => extract_scalar_field( state( iphase ), trim(field_name) // 'Absorption', stat )
         have_absorption = ( stat == 0 )
         if ( have_absorption ) then
            do j = 1, node_count( field_absorption )
               field_prot_absorption( ( iphase - 1 ) * node_count( field_absorption ) + j ) = &
                    field_absorption % val( j )
            end do
         end if
      end if Conditional_AbsorptionField

      pressure => extract_scalar_field( state( 1 ), 'Pressure' )
      pmesh => extract_mesh( state, 'PressureMesh' )
      cmesh => extract_mesh( state, 'CoordinateMesh' )
      positions => extract_vector_field( state( 1 ), 'Coordinate' )

      snloc = face_loc( pressure, 1 ) ; stotel = surface_element_count( cmesh ) ; &
           nstate = option_count( '/material_phase' )

      ncomp = 0
      do i = 1, nstate
         if( have_option( '/material_phase[' // int2str( i - 1) // ']/is_multiphase_component' ) )then
            ncomp = ncomp + 1
         end if
      end do
      nphase = nstate - ncomp
      option_path = '/material_phase['//int2str( iphase - 1 )//']/scalar_field::'//trim( field_name )
      ewrite(3,*)'option_path:', trim( option_path )


!!$ This will need to be ammended later on to take into account python functions that impose 
!!$ time-dependent field changes
      Conditional_InitialisationFromFLML: if( initialised ) then ! Extracting from state after initialisation 
         field_prot = field % val
!!$         field_prot( 1 : node_count( field ) ) = field % val
!!$         field_prot( ( iphase - 1 ) * node_count( field ) + 1 : iphase * node_count( field ) ) = &
!!$              field % val

      else !Initialisation before adapt
         if( have_option( trim( option_path ) // '/prognostic/initial_condition::WholeMesh/constant' ) )then
            call get_option(trim( option_path ) // '/prognostic/initial_condition::WholeMesh/constant', &
                 initial_constant )
            field_prot = initial_constant

         elseif( have_option( trim( option_path ) // '/prognostic/initial_condition::WholeMesh/python ') )then
            call get_option( trim( option_path ) // '/prognostic/initial_condition::WholeMesh/python', func )
            call allocate( dummy, field % mesh, 'dummy' )
            call get_option('/timestepping/current_time', current_time)
            call set_from_python_function(dummy, trim(func), positions, current_time)
            field_prot = dummy % val
            call deallocate( dummy )
         elseif( have_option( trim( option_path ) // '/prognostic/initial_condition/from_file')) then
            field_prot = field % val

         else if (have_option( trim( option_path ) // '/prognostic/initial_condition') )then
         call allocate( dummy, field % mesh, 'dummy' )
         call get_option('/timestepping/current_time', current_time)
         call initialise_field_over_regions(dummy, trim( option_path ) // '/prognostic/initial_condition', positions, current_time)
         field_prot = dummy%val
         call deallocate( dummy )

         else
            ewrite(-1,*) 'No initial condition for field::', trim( field_name )
            FLAbort( 'Check initial conditions' )
         end if
      end if Conditional_InitialisationFromFLML

!!$ Boundary conditions
      if (present( wic_bc ) ) then
      
      option_path2 = trim( option_path ) // '/prognostic/boundary_conditions['
      nobcs = get_boundary_condition_count( field )
      Loop_BC: do k = 1, nobcs

         option_path = trim( option_path2 ) // int2str( k - 1 ) // ']/surface_ids'
         shape_option = option_shape( trim( option_path ) )
         allocate( SufID_BC( 1 : shape_option( 1 ) ) )
         call get_option( trim( option_path ), SufID_BC )
         allocate( face_nodes( face_loc( field, 1) ) )

         option_path = trim( option_path2 ) // int2str( k - 1 ) // ']/'

         Conditional_Field_BC: if( have_option( trim( option_path ) // 'type::dirichlet' ) ) then

            BC_Type = 1
            field_prot_bc => extract_surface_field( field, k, 'value' )

            sele = 1
            do j = 1, stotel
               if( any ( SufID_BC == pmesh % faces % boundary_ids( j ) ) ) then
                  wic_bc( j + ( iphase - 1 ) * stotel ) = BC_Type
                  face_nodes = ele_nodes( field_prot_bc, sele )
                  do kk = 1, snloc
                     suf_bc( ( iphase - 1 ) * stotel * snloc + ( j - 1 ) * snloc + kk ) = &
                          field_prot_bc % val( face_nodes( kk ) )
                  end do
                  sele = sele + 1
               end if
            end do

         else if( have_option( trim( option_path ) // 'type::robin' ) ) then

            BC_Type = 2

            do j = 1, stotel
               if( any ( SufID_BC == pmesh % faces % boundary_ids( j ) ) ) then
                  wic_bc( j + ( iphase - 1 ) * stotel ) = BC_Type
               end if
            end do

            ! calculate this later on...
            suf_bc_rob1 = 0.
            suf_bc_rob2 = 0.

         end if Conditional_Field_BC

         deallocate( face_nodes, sufid_bc )

      end do Loop_BC

   end if

      return
    end subroutine Get_ScalarFields_Outof_State


    subroutine Get_CompositionFields_Outof_State( state, initialised, nphase, icomp, iphase, field, &
         field_prot, wic_bc, &
         kprime, kprime2, &
         suf_bc, &
         field_prot_source, field_prot_absorption )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      logical, intent( in ) :: initialised
      integer, intent( in ) :: nphase, icomp, iphase
      type( scalar_field ), pointer :: field, field_prot_bc
      real, dimension( : ), intent( inout ) :: field_prot
      real, dimension( : ), intent( inout ), optional :: field_prot_source, field_prot_absorption
      integer, dimension( : ), intent( inout ), optional :: wic_bc
      integer, intent( in ), optional  :: kprime, kprime2
      real, dimension(  :  ), intent( inout ), optional  :: suf_bc
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

      Conditional_InitialisedFromFLML: if( initialised ) then
         field_prot = field % val
      else
         !option_path = '/material_phase[' // int2str( icomp - 1 ) // &
         !     ']/scalar_field::ComponentMassFractionPhase' // &
         !     int2str( iphase )

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
         elseif( have_option( trim( option_path ) // '/prognostic/initial_condition/from_file')) then
            field_prot = field % val
         else
            ewrite(-1,*) 'No initial condition for field::', trim( field_name )
            FLAbort( ' Check initial conditions ' )

         end if Conditional_Composition_MassFraction

      end if Conditional_InitialisedFromFLML
      if ( present(wic_bc)) then  

      
      Conditional_Composition_BC: if ( have_option( trim( option_path ) // &
           '/prognostic/boundary_conditions[0]/type::dirichlet' )) then

         BC_Type = 1
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
                          field_prot_bc % val( face_nodes( kk ) )
                  end do
                  sele = sele + 1
               end if
            end do

            deallocate( face_nodes )
            deallocate( sufid_bc )

         end do Loop_Over_BC ! End of BC loop

      end if Conditional_Composition_BC

   end if


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



    subroutine Get_VectorFields_Outof_State( state, initialised, iphase, field, &
         wic_bc, wic_momu_bc, suf_u_bc, suf_v_bc, suf_w_bc, &
         suf_momu_bc, suf_momv_bc, suf_momw_bc, &
         field_prot_source, field_prot_absorption )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      logical, intent( in ) :: initialised
      integer, intent( in ) :: iphase
      type( vector_field ), pointer :: field, field_prot_bc
 !     real, dimension( : ), intent( inout ) :: field_u_prot, field_v_prot, field_w_prot, &
 !          field_nu_prot, field_nv_prot, field_nw_prot
      real, dimension( : ), intent( inout ), optional :: field_prot_source
      real, dimension( : , :, : ), intent( inout ), optional :: field_prot_absorption
      integer, dimension( : ), intent( inout ), optional :: wic_bc, wic_momu_bc 
      real, dimension( : ), intent( inout ), optional :: suf_u_bc, suf_v_bc, suf_w_bc
      real, dimension( : ), intent( inout ), optional :: suf_momu_bc, suf_momv_bc, suf_momw_bc

      ! Local variables
      type( mesh_type ), pointer :: pmesh, cmesh
      type(vector_field) :: dummy
      type(vector_field), pointer :: positions, field_source, field_absorption, chunk
      type(scalar_field), pointer :: pressure
      integer, dimension(:), allocatable :: sufid_bc, face_nodes
      character( len = option_path_len ) :: option_path, option_path2, field_name, bct
      integer :: ndim, stotel, snloc, snloc2, nonods, nobcs, bc_type, i,j, k, kk, l, &
           shape_option( 2 ), count, u_nonods, idim, stat,nlev,nloc, ele
      real, dimension( : ), allocatable :: initial_constant
      integer, dimension(:), allocatable :: order
      logical :: have_source, have_absorption
      character(len=8192) :: func

      pmesh => extract_mesh(state, 'PressureMesh' )
      cmesh => extract_mesh(state, 'CoordinateMesh' )
      positions => extract_vector_field( state( 1 ), 'Coordinate' )

      ndim = field % dim
      stotel = surface_element_count( cmesh )
      snloc2 = face_loc( field, 1)
      snloc = snloc2
      if ( is_overlapping ) snloc = snloc2 * ele_loc( pmesh, 1)
      nonods = node_count( field )
      field_name = trim( field % name )
      u_nonods = nonods
      nlev=1
      if ( is_overlapping ) then
         u_nonods = nonods * ele_loc( pmesh, 1)
         nlev=ele_loc( pmesh, 1)
      end if
      


      have_absorption = .false.
      Conditional_AbsorptionField: if( present( field_prot_absorption ) ) then
         field_absorption => extract_vector_field( state( iphase ), trim(field_name) // 'Absorption', stat )
         option_path = '/material_phase[' // int2str( iphase - 1 ) // ']/vector_field::' // trim( field_name ) // &
              '/prognostic/vector_field::Absorption/diagnostic/algorithm::vector_python_diagnostic'
         have_absorption =  have_option( trim(option_path) )
         if ( have_absorption ) then
            do idim = 1, ndim
               field_prot_absorption( :, idim + (iphase-1)*ndim, idim + (iphase-1)*ndim ) =  &
                    field_absorption % val( idim, : )
            end do
         else
            do idim = 1, ndim
               field_prot_absorption( :, idim + (iphase-1)*ndim, idim + (iphase-1)*ndim ) = 0.0 
            end do
         end if
      end if Conditional_AbsorptionField

      if (present(wic_bc)) then

      option_path = '/material_phase[' // int2str( iphase - 1 )// ']/vector_field::' // trim( field_name )
      option_path2 = trim( option_path ) // '/prognostic/boundary_conditions['

      nobcs = get_boundary_condition_count( field )
      Loop_BC: do k = 1, nobcs

         field_prot_bc => extract_surface_field( field, k, 'value' )

         option_path = trim( option_path2 ) // int2str( k - 1 ) // ']/surface_ids'
         shape_option = option_shape( trim( option_path ) )
         allocate( SufID_BC( 1 : shape_option( 1 ) ) )
         call get_option( trim( option_path ), SufID_BC )
         allocate( face_nodes( face_loc( field, 1) ) )

         option_path = trim( option_path2 ) // int2str( k - 1 ) // ']/'

         Conditional_Field_BC: if( have_option( trim( option_path ) // 'type::dirichlet' ) ) then

            BC_Type = 1

            face_nodes = (/ ( l, l = 1, snloc2 ) /)

            do j = 1, stotel
               if( any ( sufid_bc == field % mesh % faces % boundary_ids( j ) ) ) then 
                  wic_bc( j  + ( iphase - 1 ) * stotel ) = BC_Type
                  count = 1
                  do kk = 1, snloc
                     suf_u_bc( ( iphase - 1 ) * stotel * snloc + ( j - 1 ) * snloc + kk ) = &
                          field_prot_bc % val( 1, face_nodes( count ) )
                     if( ndim > 1 ) suf_v_bc( ( iphase - 1 ) * stotel * snloc + & 
                          ( j - 1 ) * snloc + kk ) = field_prot_bc % val( 2, face_nodes( count ) )
                     if( ndim > 2 ) suf_w_bc( ( iphase - 1 ) * stotel * snloc + &
                          ( j - 1 ) * snloc + kk ) = field_prot_bc % val( 3, face_nodes( count ) )
                     count = count + 1
                     if ( mod( kk, snloc2 ) == 0. ) count = 1
                  end do
                  face_nodes = face_nodes + snloc2
               end if
            end do

         else if( have_option( trim( option_path ) // 'type::momentum' ) ) then

            call get_option( trim( option_path ) // 'type::momentum/boundary/', bct )
            if ( trim( bct ) == "incoming") then
               BC_Type = 1
            else if ( trim( bct ) == "open") then
               BC_Type = 5
            else
               FLAbort( 'Wrong momentum boundary condition in diamond.' )
            end if

            face_nodes = (/ ( l, l = 1, snloc2 ) /)
            do j = 1, stotel
               if( any ( sufid_bc == field % mesh % faces % boundary_ids( j ) ) ) then 
                  wic_momu_bc( j  + ( iphase - 1 ) * stotel ) = BC_Type
                  count = 1
                  do kk = 1, snloc
                     suf_momu_bc( ( iphase - 1 ) * stotel * snloc + ( j - 1 ) * snloc + kk ) = &
                          field_prot_bc % val( 1, face_nodes( count ) )
                     if( ndim > 1 ) suf_momv_bc( ( iphase - 1 ) * stotel * snloc + & 
                          ( j - 1 ) * snloc + kk ) = field_prot_bc % val( 2, face_nodes( count ) )
                     if( ndim > 2 ) suf_momw_bc( ( iphase - 1 ) * stotel * snloc + &
                          ( j - 1 ) * snloc + kk ) = field_prot_bc % val( 3, face_nodes( count ) )
                     count = count + 1
                     if ( mod( kk, snloc2 ) == 0. ) count = 1
                  end do
                  face_nodes = face_nodes + snloc2
               end if
            end do

         end if Conditional_Field_BC

         deallocate( face_nodes, SufID_BC )

      end do Loop_BC

   end if

      return
    end subroutine Get_VectorFields_Outof_State



    subroutine Extract_TensorFields_Outof_State( state, istate_field, &
         field, field_path, &
         field_prot_tensor, &
         GlobalNodeNumber )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      integer, intent( in ) :: istate_field
      type( tensor_field ), pointer :: field
      character( len = option_path_len ), intent( in ) :: field_path
      real, dimension( :, :, : ), intent( inout ) :: field_prot_tensor
      integer, dimension( : ), intent( in ), optional :: GlobalNodeNumber
!!$ Local variables
      type( scalar_field ), pointer :: pressure
      type( vector_field ), pointer :: positions
      type( tensor_field ), pointer :: tensorfield
      integer :: nstate, ndim, totele, istate, ele, idim, jdim, iloc, inod
      character( len = option_path_len ) :: field_name, option_path, option_permeability, option_viscosity
      logical :: compute_viscosity, compute_permeability, is_isotropic, is_diagonal, is_anisotropic
      real, dimension( :, : ), allocatable :: constant
      integer, dimension( : ), pointer :: element_nodes

      call get_option( '/geometry/dimension', ndim )
      nstate = option_count( '/material_phase' )
      positions => extract_vector_field( state, 'Coordinate' )
      pressure => extract_scalar_field( state, 'Pressure' )

      totele = ele_count( positions ) ; field_name = trim( field % name )
      compute_viscosity = .false. ; compute_permeability = .false.

!!$ Defining logicals for the field associated with the tensor. For now, it is set up for either 
!!$ Permeability( totele, ndim, ndim ) and Viscosity( mat_nonods, ndim, ndim, nphase ). We may need to
!!$ extend the subroutine to also with an arbitrary tensor shape.
      Conditional_WhichField: if( trim( field_path ) == '/material_phase[' // int2str( istate_field - 1 ) // &
           ']/vector_field::Velocity/prognostic/tensor_field::' // trim( field_name ) )then
         compute_viscosity = .true.
      elseif( trim( field_path ) == '/porous_media/tensor_field::Permeability' ) then 
         compute_permeability = .true.
      else
         FLAbort( 'Tensor Field was not defined in the schema - sort this out' )
      end if Conditional_WhichField


!!$ Now looping over state
      LoopOverState: do istate = 1, nstate

         if( istate /= istate_field ) cycle LoopOverState        

!!$ Defining flags:
         if( compute_permeability ) then
            option_permeability = '/porous_media/tensor_field::Permeability/prescribed/value[0]'
            if( have_option( trim( option_permeability ) ) ) then
               option_path = trim( option_permeability )
            else
               FLAbort( 'Option path for the Permeability tensor is incomplete.' )
            end if
         end if

         if( compute_viscosity ) then
            option_viscosity = '/material_phase[' // int2str( istate - 1 ) // ']/vector_field::' // &
                 'Velocity/prognostic/tensor_field::' // trim( field_name ) // '/prescribed/value[0]'
            if( have_option( trim( option_viscosity ) ) ) then
               option_path = trim( option_viscosity )
            else
               FLAbort( 'Option path for the Viscosity tensor is incomplete.' )
            end if
         end if
!!$
!!$         Conditional_InitialisationFromFLML:if ( initialised ) then
!!$               tensorfield => extract_tensor_field( state( istate ), trim( field_name ) )
!!$            do ele = 1, element_count( tensorfield )
!!$               element_nodes => ele_nodes( tensorfield, ele )
!!$               do idim = 1, ndim
!!$                  do jdim = 1, ndim 
!!$                     if( compute_permeability )then
!!$                        field_prot_tensor( ele, idim, jdim ) = field % val( idim, jdim, element_nodes( 1 ) )
!!$                     elseif( compute_viscosity )then
!!$                         do iloc = 1, ele_loc( pressure, 1 )
!!$                              inod = GlobalNodeNumber( ( ele - 1 ) * node_count( pressure ) + iloc )
!!$                              field_prot_tensor( inod, idim, jdim ) = &
!!$                                   tensorfield % val( idim, jdim, element_nodes( 1 ) )                               
!!$                         end do
!!$                     else
!!$                        FLAbort( 'Incorrect path for the tensor field' )
!!$                     end if
!!$                  end do
!!$               end do
!!$            end do
!!$         end if Conditional_InitialisationFromFLML
!!$
!!$

         is_isotropic = have_option( trim( option_path ) // '/isotropic' )
         is_diagonal = have_option( trim( option_path ) // '/diagonal' )
         is_anisotropic = have_option( trim( option_path ) // '/anisotropic_symmetric' ) .or. &
              have_option( trim( option_path ) // '/anisotropic_asymmetric' )
!!$ Endof Defining flags

!!$ Isotropic Tensor:
         Conditional_Tensor: if ( is_isotropic ) then
            option_path = trim( option_path ) // '/isotropic'

            if( have_option( trim( option_path ) // '/constant')) then
               allocate( constant( 1, 1 ) ) ; constant = 0.
               call get_option( trim( option_path ) // '/constant', constant( 1, 1 ) )
               do idim = 1, ndim
                  field_prot_tensor( : , idim, idim ) = constant( 1, 1 )
               end do
               deallocate( constant )

            elseif( have_option( trim( option_path ) // '/python' ) ) then
               tensorfield => extract_tensor_field( state( istate ), trim( field_name ) )
               do ele = 1, element_count( tensorfield )
                  element_nodes => ele_nodes( tensorfield, ele )
                  do idim = 1, ndim
                     if( compute_permeability ) then
                        field_prot_tensor( ele, idim, idim ) = &
                             tensorfield % val( idim, idim, element_nodes( 1 ) )
                     elseif( compute_viscosity ) then
                        do iloc = 1, ele_loc( pressure, 1 )
                           inod = GlobalNodeNumber( ( ele - 1 ) * node_count( pressure ) + iloc )
                           field_prot_tensor( inod , idim, idim ) =   &
                                tensorfield % val( idim, idim, element_nodes( 1 ) )
                        end do
                     else
                        FLAbort( 'Option path for the tensor field is incomplete.' )
                     end if
                  end do
               end do

            else
               FLExit( 'Incorrect initial condition for field' )
            end if
!!$ Endof Isotropic Tensor

!!$ Diagonal Tensor:
         elseif( is_diagonal )then ! If_Conditional_Tensor
            option_path = trim( option_path ) // '/diagonal'

            if( have_option( trim( option_path ) // '/constant' ) ) then
               allocate(constant( ndim, 1 ) ) ; constant = 0.
               call get_option( trim( option_path ) // '/constant', constant )
               do idim = 1, ndim
                  field_prot_tensor( : , idim, idim ) = constant( idim, 1 )
               end do
               deallocate( constant )

            elseif( have_option( trim( option_path ) // '/python' ) )then
               tensorfield => extract_tensor_field( state( istate ), trim( field_name ) )
               element_nodes => ele_nodes( tensorfield, ele )
               do idim = 1, ndim
                  do ele = 1, element_count( tensorfield )
                     if( compute_permeability ) then
                        field_prot_tensor( ele, idim, idim ) = tensorfield % val( idim, idim,  ele )
                     elseif( compute_viscosity ) then
                        do iloc = 1, ele_loc( pressure, 1 )
                           inod = GlobalNodeNumber( ( ele - 1 ) * node_count( pressure ) + iloc )
                           field_prot_tensor( inod , idim, idim ) =   &
                                tensorfield % val( idim, idim, element_nodes( 1 ) )
                        end do
                     else
                        FLAbort( 'Option path for the tensor field is incomplete.' )
                     end if
                  end do
               end do

            else
               FLExit( 'Incorrect initial condition for field' )
            end if
!!$ Endof Diagonal Tensor

!!$ Anisotropic Tensor:
         elseif( is_anisotropic )then  ! If_Conditional_Tensor
            if( have_option( trim( option_path ) // '/anisotropic_symmetric' ) ) then
               option_path = trim( option_path ) // '/anisotropic_symmetric'
            elseif( have_option( trim( option_path ) // '/anisotropic_asymmetric' ) ) then
               option_path = trim( option_path ) // '/anisotropic_asymmetric'
            else
               FLAbort( 'Option path for the tensor field is incomplete.' )
            end if

            if( have_option( trim( option_path ) // '/constant' ) ) then
               allocate( constant( ndim, ndim ) ) ; constant = 0.
               call get_option( trim( option_path ) // '/constant', constant )
               do idim = 1, ndim
                  do jdim = 1, ndim
                     field_prot_tensor( : , idim, jdim ) = constant( idim, jdim )
                  end do
               end do

            elseif( have_option( trim( option_path ) // '/python' ) ) then
               tensorfield => extract_tensor_field( state( istate ), trim( field_name ) )
               do ele = 1, element_count( tensorfield )
                  element_nodes => ele_nodes( tensorfield, ele )
                  do idim = 1, ndim
                     do jdim = 1, ndim
                        if( compute_permeability ) then
                           field_prot_tensor( ele, idim, jdim ) = &
                                tensorfield % val( idim, jdim, element_nodes( 1 ) )
                        elseif( compute_viscosity ) then
                           do iloc = 1, ele_loc( pressure, 1 )
                              inod = GlobalNodeNumber( ( ele - 1 ) * node_count( pressure ) + iloc )
                              field_prot_tensor( inod, idim, jdim ) = &
                                   tensorfield % val( idim, jdim, element_nodes( 1 ) )
                           end do
                        else 
                           FLAbort( 'Option path for the tensor field is incomplete.' )
                        end if
                     end do
                  end do
               end do

            else 
               FLExit( 'Incorrect initial condition for field' )

            end if
!!$ Endof Anisotropic Tensor

         else
            FLExit( 'Incorrect initial condition for field' )

         end if Conditional_Tensor

      end do LoopOverState

      return
    end subroutine Extract_TensorFields_Outof_State


!!$
!!$ Module Get_Ndgln Interfaces

    subroutine Get_Scalar_Ndgln( ndgln, field, cv_nloc )
      implicit none
      type( scalar_field ), intent( in ) :: field
      integer, intent( in ), optional :: cv_nloc
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

    subroutine Get_Vector_Ndgln( ndgln, field, cv_nloc )
      implicit none
      type( vector_field ), intent( in ) :: field
      integer, intent( in ), optional :: cv_nloc
      integer, dimension( : ), intent( inout ) :: ndgln
      ! Local variables
      integer, dimension( : ), pointer :: nloc
      integer :: ele, iloc, count, cv_nloc2, ndim

      call get_option( '/geometry/dimension', ndim )
      cv_nloc2 = 1

      if ( trim( field % name ) == 'Velocity' ) then
         if ( is_overlapping .and. ( ndim > 1 ) ) cv_nloc2 = cv_nloc
      end if

      count = 0
      do ele = 1, ele_count( field )
         nloc => ele_nodes( field, ele )
         do iloc = 1, ele_loc( field, 1 ) * cv_nloc2
            if( is_overlapping ) then
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


    subroutine Get_Mesh_Ndgln( ndgln, mesh, cv_nloc )
      implicit none
      type( mesh_type ), intent( in ) :: mesh
      integer, intent( in ), optional :: cv_nloc
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

    subroutine Get_Scalar_SNdgln( sndgln, field, cv_nloc  )
      implicit none
      type( scalar_field ), intent( in ) :: field
      integer, dimension( : ), intent( inout ) :: sndgln
      integer, intent( in ), optional :: cv_nloc
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

    subroutine Get_Vector_SNdgln( sndgln, field, cv_nloc  )
      implicit none
      type( vector_field ), intent( in ) :: field
      integer, dimension( : ), intent( inout ) :: sndgln
      integer, intent( in ), optional :: cv_nloc
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




    subroutine print_from_state( state, field_prot )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      real, dimension( : ), intent( in ) :: field_prot
      !
      type( scalar_field ), pointer :: field
      integer :: nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, knod, istate
      real :: dx
      logical :: initialised
      integer, dimension( : ), pointer :: cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln
      integer, dimension( : ), allocatable ::     cv_sndgln, p_sndgln, u_sndgln, Temperature_BC_Spatial
      real, dimension( : ), allocatable :: Temperature, Temperature_BC 

      call Get_Primary_Scalars( state, &         
           nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, dx )

!!$ Calculating Global Node Numbers
      allocate( cv_sndgln( stotel * cv_snloc ), p_sndgln( stotel * p_snloc ), &
           u_sndgln( stotel * u_snloc ) )

  !    x_ndgln_p1 = 0 ; x_ndgln = 0 ; cv_ndgln = 0 ; p_ndgln = 0 ; mat_ndgln = 0 ; u_ndgln = 0 ; xu_ndgln = 0 ; &
           cv_sndgln = 0 ; p_sndgln = 0 ; u_sndgln = 0

      allocate( temperature( nphase * cv_nonods ), temperature_bc_spatial( nphase * stotel ), &
           temperature_bc( stotel * cv_snloc * nphase ) )

      call Compute_Node_Global_Numbers( state, &
           totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
           cv_snloc, p_snloc, u_snloc, &
           cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln, &
           cv_sndgln, p_sndgln, u_sndgln )

      initialised = .true.
      do istate = 1, nstate
         Conditional_Temperature: if( have_option( '/material_phase[' // int2str( istate - 1 ) // &
              ']/scalar_field::Temperature' ) ) then
            field => extract_scalar_field( state( istate ), 'Temperature' )
            knod = ( istate - 1 ) * node_count( field )
            call Get_ScalarFields_Outof_State( state, initialised, istate, field, &
                 Temperature( knod + 1 : knod + node_count( field ) ), &
                 Temperature_BC_Spatial, Temperature_BC ) !, &
!!$                 field_prot_source = Temperature_Source( knod + 1 : knod + node_count( field ) ) )
         end if Conditional_Temperature
      end do

      ewrite(3,*)'::temperature::', norm2( temperature )
      do istate = 1, cv_nonods
         ewrite(3,*) istate, field_prot( istate ), temperature( istate ), field%val(istate)
      end do

      deallocate( cv_sndgln, p_sndgln, u_sndgln, Temperature_BC_Spatial, Temperature, Temperature_BC ) 

      return
    end subroutine print_from_state


    subroutine update_boundary_conditions( state, stotel, cv_snloc, nphase, & 
         &                                 suf_t_bc, suf_t_bc_rob1, suf_t_bc_rob2 )
      implicit none
      type( state_type ), dimension( : ), intent( in ) :: state
      integer, intent( in ) :: stotel, cv_snloc, nphase
      real, dimension( stotel * cv_snloc * nphase ), intent( inout ) :: suf_t_bc, suf_t_bc_rob1, suf_t_bc_rob2
      !
      character( len = option_path_len ) :: option_path, option_path2, field_name, name
      integer :: shape_option(2), iphase, nobcs, kk, k, j , sele, stat
      integer, dimension( : ), allocatable :: SufID_BC
      integer, dimension( : ), pointer:: surface_element_list, face_nodes

      type( scalar_field ), pointer :: field, field_prot_bc, field_prot1, field_prot2
      type( scalar_field ), pointer :: field_prot_bc1, field_prot_bc2
      type( scalar_field ) :: field_prot_bc1f, field_prot_bc2f
      type( mesh_type ), pointer :: pmesh, surface_mesh

      suf_t_bc = 0. ; suf_t_bc_rob1 = 0. ; suf_t_bc_rob2 = 0. 

      call set_boundary_conditions_values( state, shift_time = .true. )
     
      pmesh => extract_mesh( state, 'PressureMesh' )

      do iphase = 1, nphase

         field_name = 'Temperature'
         field => extract_scalar_field( state( iphase ), trim( field_name ) )

         option_path = '/material_phase['//int2str( iphase - 1 )//']/scalar_field::'//trim( field_name )
         option_path2 = trim( option_path ) // '/prognostic/boundary_conditions['

         nobcs = get_boundary_condition_count( field )

         Loop_BC: do k = 1, nobcs

            option_path = trim( option_path2 ) // int2str( k - 1 ) // ']/surface_ids'
            shape_option = option_shape( trim( option_path ) )
            allocate( SufID_BC( 1 : shape_option( 1 ) ) )
            call get_option( trim( option_path ), SufID_BC )

            option_path = trim( option_path2 ) // int2str( k - 1 ) // ']/'

            Conditional_Field_BC: if( have_option( trim( option_path ) // 'type::dirichlet' ) ) then

               !BC_Type = 1
               field_prot_bc => extract_surface_field( field, k, 'value' )

               sele = 1
               do j = 1, stotel
                  if( any ( SufID_BC == pmesh % faces % boundary_ids( j ) ) ) then
                     !wic_bc( j + ( iphase - 1 ) * stotel ) = BC_Type
                     face_nodes => ele_nodes( field_prot_bc, sele )
                     do kk = 1, cv_snloc
                        suf_t_bc( ( iphase - 1 ) * stotel * cv_snloc + ( j - 1 ) * cv_snloc + kk ) = &
                             field_prot_bc % val( face_nodes( 1 ) )
                     end do
                     sele = sele + 1
                  end if
               end do

            else if( have_option( trim( option_path ) // 'type::robin' ) ) then

               !BC_Type = 2
               if( have_option( trim( option_path ) // 'type::robin/order_zero_coefficient/from_field' ) ) then
             
                  call get_boundary_condition( field, k, surface_mesh = surface_mesh, &
                       surface_element_list = surface_element_list )

                  call allocate( field_prot_bc1f, surface_mesh, "Robin1" )
                  call allocate( field_prot_bc2f, surface_mesh, "Robin2" )

                  call get_option( trim( option_path ) // "type::robin/order_zero_coefficient/from_field/name", name )
                  field_prot1 => extract_scalar_field( state( iphase ), name, stat )
                  if(stat /= 0) FLExit( "Could not extract parent field 1. Check options file?" )

                  call get_option( trim( option_path ) // "type::robin/order_one_coefficient/from_field/name", name )
                  field_prot2 => extract_scalar_field( state( iphase ), name, stat )
                  if(stat /= 0) FLExit( "Could not extract parent field 2. Check options file?" )

                  call remap_field_to_surface( field_prot1, field_prot_bc1f, surface_element_list )
                  call remap_field_to_surface( field_prot2, field_prot_bc2f, surface_element_list )

                  ! copy back memory
                  sele = 1
                  do j = 1, stotel
                     if( any ( SufID_BC == pmesh % faces % boundary_ids( j ) ) ) then
                        !wic_bc( j + ( iphase - 1 ) * stotel ) = BC_Type
                        face_nodes => ele_nodes( field_prot_bc1f, sele )
                        do kk = 1, cv_snloc
                           suf_t_bc_rob1( ( iphase - 1 ) * stotel * cv_snloc + ( j - 1 ) * cv_snloc + kk ) = &
                                field_prot_bc1f % val( face_nodes( 1 ) )
                           suf_t_bc_rob2( ( iphase - 1 ) * stotel * cv_snloc + ( j - 1 ) * cv_snloc + kk ) = &
                                field_prot_bc2f % val( face_nodes( 1 ) )
                        end do
                        sele = sele + 1
                     end if
                  end do

                  call deallocate( field_prot_bc1f )
                  call deallocate( field_prot_bc2f )

               else

                  field_prot_bc1 => extract_surface_field( field, k, 'order_zero_coefficient' )
                  field_prot_bc2 => extract_surface_field( field, k, 'order_one_coefficient' )

                  sele = 1
                  do j = 1, stotel
                     if( any ( SufID_BC == pmesh % faces % boundary_ids( j ) ) ) then
                        !wic_bc( j + ( iphase - 1 ) * stotel ) = BC_Type
                        face_nodes => ele_nodes( field_prot_bc1, sele )
                        do kk = 1, cv_snloc
                           suf_t_bc_rob1( ( iphase - 1 ) * stotel * cv_snloc + ( j - 1 ) * cv_snloc + kk ) = &
                                field_prot_bc1 % val( face_nodes( 1 ) )
                           suf_t_bc_rob2( ( iphase - 1 ) * stotel * cv_snloc + ( j - 1 ) * cv_snloc + kk ) = &
                                field_prot_bc2 % val( face_nodes( 1 ) )
                        end do
                        sele = sele + 1
                     end if
                  end do

               end if

            end if Conditional_Field_BC

            deallocate( SufID_BC )

         end do Loop_BC

      end do

    end subroutine update_boundary_conditions

    subroutine pack_multistate(state, packed_state,&
         multiphase_state, multicomponent_state, pmulti_state ) 
      
      type(state_type), dimension(:), intent(inout):: state
      type(state_type), dimension(:), intent(inout), pointer :: &
           multiphase_state, multicomponent_state
      type(state_type) ::  packed_state 
    
      type(state_type), dimension(:,:), pointer, optional :: pmulti_state

      type(state_type), dimension(:,:), pointer :: multi_state
      
      integer :: i,nphase,ncomp,ndim, stat, iphase, icomp, idim
      
      type(scalar_field), pointer :: pressure, p2,sfield
      type(vector_field), pointer :: velocity, position, vfield
      
      type(scalar_field) :: porosity
      type(vector_field) :: p_position, u_position, m_position
      type(tensor_field) :: permeability
      type(mesh_type) :: ovmesh,lmesh,nvmesh, element_mesh
      type(element_type) :: overlapping_shape, vel_shape, element_shape
      character( len = option_path_len ) :: vel_element_type

      logical :: has_density, has_phase_volume_fraction

      ncomp=option_count('/material_phase/is_multiphase_component')
      nphase=size(state)-ncomp
      position=>extract_vector_field(state(1),"Coordinate")
      ndim=mesh_dim(position)

      if (has_scalar_field(state(1),"Porosity")) then
         sfield=>extract_scalar_field(state(1),"Porosity")
         element_mesh=sfield%mesh
      else
         element_shape=make_element_shape(position%mesh%shape,degree=0)
         element_mesh=make_mesh(position%mesh,element_shape,&
              continuity=-1,name="ElementMesh")
      end if

      call get_option('/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           vel_element_type )

      allocate(multiphase_state(nphase))
      allocate(multicomponent_state(ncomp))
      allocate(multi_state(max(1,ncomp),nphase))

      pressure=>extract_scalar_field(state(1),"Pressure")
      call insert(packed_state,pressure,"Pressure")
      call insert(packed_state,pressure%mesh,"PressureMesh")

      call add_new_memory(packed_state,pressure,"FEPressure")
      call add_new_memory(packed_state,pressure,"OldFEPressure")
      call add_new_memory(packed_state,pressure,"CVPressure")
      call add_new_memory(packed_state,pressure,"OldCVPressure")

      p2=>extract_scalar_field(packed_state,"FEPressure")
      call set( p2, pressure )

      p2=>extract_scalar_field(packed_state,"CVPressure")
      call set( p2, pressure )


      call insert_sfield(packed_state,"FEDensity",1,nphase)

      call insert_sfield(packed_state,"Density",1,nphase)
      call insert_sfield(packed_state,"DensityHeatCapacity",1,nphase)


      if (option_count("/material_phase/scalar_field::Temperature")>0) then
         call insert_sfield(packed_state,"Temperature",1,nphase,&
	         add_source=.true.,add_absorption=.true.)
          call insert_sfield(packed_state,"FETemperature",1,nphase)
      end if
      call insert_sfield(packed_state,"PhaseVolumeFraction",1,nphase,&
           add_source=.true.)
      call insert_sfield(packed_state,"FEPhaseVolumeFraction",1,nphase)

      velocity=>extract_vector_field(state(1),"Velocity")
      call insert(packed_state,velocity%mesh,"VelocityMesh")
      ovmesh=make_mesh(position%mesh,&
           shape=velocity%mesh%shape,&
           continuity=1,name="VelocityMesh_Continuous")
      call insert(packed_state,ovmesh,"VelocityMesh_Continuous")
      if ( .not. has_mesh(state(1),"VelocityMesh_Continuous") ) &
           call insert(state(1),ovmesh,"VelocityMesh_Continuous")
      call allocate(u_position,ndim,ovmesh,"VelocityCoordinate")
      ovmesh=make_mesh(position%mesh,&
           shape=pressure%mesh%shape,&
           continuity=1,name="PressureMesh_Continuous")
      call insert(packed_state,ovmesh,"PressureMesh_Continuous")
      if ( .not. has_mesh(state(1),"PressureMesh_Continuous") ) &
           call insert(state(1),ovmesh,"PressureMesh_Continuous")

      call allocate(p_position,ndim,ovmesh,"PressureCoordinate")
      ovmesh=make_mesh(position%mesh,&
           shape=pressure%mesh%shape,&
           continuity=-11,name="PressureMesh_Discontinuous")
      call insert(packed_state,ovmesh,"PressureMesh_Discontinuous")
      if ( .not. has_mesh(state(1),"PressureMesh_Discontinuous") ) &
           call insert(state(1),ovmesh,"PressureMesh_Discontinuous")
      call allocate(m_position,ndim,ovmesh,"MaterialCoordinate")    

      call remap_field( position, u_position )
      call remap_field( position, p_position )
      call remap_field( position, m_position )
  
      call insert(packed_state,position,"Coordinate")
      call insert(packed_state,u_position,"VelocityCoordinate")
      call insert(packed_state,p_position,"PressureCoordinate")
      call insert(packed_state,m_position,"MaterialCoordinate")
      call deallocate(p_position)
      call deallocate(u_position)
      call deallocate(m_position)

      has_density=has_scalar_field(state(1),"Density")
      has_phase_volume_fraction=has_scalar_field(state(1),"PhaseVolumeFraction")

      if (trim(vel_element_type)=='overlapping') then
         !overlapping
         ovmesh=make_mesh(position%mesh,&
              shape=velocity%mesh%shape,& 
              continuity=velocity%mesh%continuity,&
              overlapping_shape=pressure%mesh%shape)
         overlapping_shape=pressure%mesh%shape

         do i=1,size(state)
            if (has_vector_field(state(i),"Velocity")) then

               velocity=>extract_vector_field(state(i),"Velocity",stat)
               velocity%mesh%shape%numbering=>&
                    get_overlapping_version(velocity%mesh%shape%numbering)
               velocity%mesh%overlapping_shape=overlapping_shape
            end if
         end do
         call insert_vfield(packed_state,"Velocity",ovmesh,add_source=.true.)
         call insert_vfield(packed_state,"NonlinearVelocity",ovmesh)
         call insert(state(1),ovmesh,"InternalVelocityMesh")
         call insert(packed_state,ovmesh,"InternalVelocityMesh")
      else
         call insert(packed_state,velocity%mesh,"InternalVelocityMesh")
         call insert_vfield(packed_state,"Velocity",add_source=.true.)
         call insert_vfield(packed_state,"NonlinearVelocity")
         call insert(state(1),velocity%mesh,"InternalVelocityMesh")
      end if

      call unpack_multiphase(packed_state,multiphase_state)  
      if (ncomp>0) then
         call insert_sfield(packed_state,"ComponentDensity",ncomp,nphase)
         call insert_sfield(packed_state,"ComponentMassFraction",ncomp,&
              nphase,add_source=.true.)

         call insert_sfield(packed_state,"FEComponentDensity",ncomp,nphase)
         call insert_sfield(packed_state,"FEComponentMassFraction",ncomp,nphase)

         call unpack_multicomponent(packed_state,multicomponent_state)  
      end if

      iphase=1
      icomp=1
      do i=1,size(state)
         if(have_option(trim(state(i)%option_path)&
              //'/is_multiphase_component')) then
            velocity=>extract_vector_field(state(i),"Velocity",stat)
            if (stat==0) velocity%wrapped=.true.
            velocity=>extract_vector_field(state(i),"OldVelocity",stat)
            if (stat==0) velocity%wrapped=.true.
            velocity=>extract_vector_field(state(i),"IteratedVelocity",stat)
            if (stat==0) velocity%wrapped=.true.


            call unpack_component_sfield(state(i),packed_state,"FEComponentDensity",icomp)
            !call unpack_component_sfield(state(i),packed_state,"FEOldComponentDensity",icomp)
            call unpack_component_sfield(state(i),packed_state,"FEComponentMassFraction",icomp)
            !call unpack_component_sfield(state(i),packed_state,"FEOldComponentMassFraction",icomp)


            call unpack_component_sfield(state(i),packed_state,"OldComponentDensity",icomp,prefix='Old')
            call unpack_component_sfield(state(i),packed_state,"OldComponentMassFraction",icomp,prefix='Old')
            call unpack_component_sfield(state(i),packed_state,"ComponentDensity",icomp)
            call unpack_component_sfield(state(i),packed_state,"ComponentMassFraction",icomp)

            call insert_components_in_multi_state(multi_state(icomp,:),state(i))

            icomp=icomp+1
            cycle
         else

            
            if (has_density) then
               call unpack_sfield(state(i),packed_state,"IteratedDensity",1,iphase,&
                    check_paired(extract_scalar_field(state(i),"Density"),&
                    extract_scalar_field(state(i),"IteratedDensity")))
               call unpack_sfield(state(i),packed_state,"OldDensity",1,iphase,&
                    check_paired(extract_scalar_field(state(i),"Density"),&
                    extract_scalar_field(state(i),"OldDensity")))
               call unpack_sfield(state(i),packed_state,"Density",1,iphase)
               call insert(multi_state(1,iphase), extract_scalar_field(state(i),"Density"),"Density")
            end if

            if(have_option(trim(state(i)%option_path)&
                 //'/scalar_field::Temperature')) then
               call unpack_sfield(state(i),packed_state,"Temperature",1,iphase)
               call unpack_sfield(state(i),packed_state,"OldTemperature",1,iphase)
               call unpack_sfield(state(i),packed_state,"IteratedTemperature",1,iphase)
               call unpack_sfield(state(i),packed_state,"TemperatureSource",1,iphase)
               call unpack_sfield(state(i),packed_state,"TemperatureAbsorption",1,iphase)
               call insert(multi_state(1,iphase), extract_scalar_field(state(i),"Temperature"),"Temperature")
            end if

            
            if(has_phase_volume_fraction) then
               call unpack_sfield(state(i),packed_state,"IteratedPhaseVolumeFraction",1,iphase,&
                    check_paired(extract_scalar_field(state(i),"PhaseVolumeFraction"),&
                    extract_scalar_field(state(i),"IteratedPhaseVolumeFraction")))
               call unpack_sfield(state(i),packed_state,"OldPhaseVolumeFraction",1,iphase,&
                    check_paired(extract_scalar_field(state(i),"PhaseVolumeFraction"),&
                    extract_scalar_field(state(i),"OldPhaseVolumeFraction")))
               call unpack_sfield(state(i),packed_state,"PhaseVolumeFraction",1,iphase)
               call unpack_sfield(state(i),packed_state,"PhaseVolumeFractionSource",1,iphase)
               call insert(multi_state(1,iphase), extract_scalar_field(state(i),"PhaseVolumeFraction"),"PhaseVolumeFraction")
            end if

            call unpack_vfield(state(i),packed_state,"IteratedVelocity",iphase,&
                 check_vpaired(extract_vector_field(state(i),"Velocity"),&
                 extract_vector_field(state(i),"IteratedVelocity")))
            call unpack_vfield(state(i),packed_state,"OldVelocity",iphase,&
                 check_vpaired(extract_vector_field(state(i),"Velocity"),&
                 extract_vector_field(state(i),"OldVelocity")))
            call unpack_vfield(state(i),packed_state,"Velocity",iphase)
            call insert(multi_state(1,iphase), extract_vector_field(state(i),"Velocity"),"Velocity")

            iphase=iphase+1
         end if
      end do

      if (option_count("/material_phase/scalar_field::Temperature")>0) call allocate_multiphase_scalar_bcs(packed_state,multi_state,"Temperature")
      call allocate_multiphase_scalar_bcs(packed_state,multi_state,"Density")
      call allocate_multiphase_scalar_bcs(packed_state,multi_state,"PhaseVolumeFraction")
      call allocate_multiphase_vector_bcs(packed_state,multi_state,"Velocity")


      call allocate_multicomponent_scalar_bcs(multicomponent_state,multi_state,"ComponentMassFraction")
      call allocate_multicomponent_scalar_bcs(multicomponent_state,multi_state,"ComponentDensity")

      !! // How to add density as a dependant of olddensity,  as an example
      !! // Suppose we have a previous declaration

      !! type(tensor_field), pointer :: old_density, density

      !! old_density=>extract_tensor_field(packed_state,"PackedOldDensity")
      !! density=>extract_tensor_field(packed_state,"PackedDensity")
      !! call add_dependant_field(old_density,density)

      !! // if we now make a  call 

      !! call mark_as_updated(density)

      !! // then is_updated(density) is now .true. (no other changes)
      !! // If we then make a  call 
      

      !! call mark_as_updated(old_density)
      !! // then is_updated(old_density) is now .true. and is_updated(density) is .false.


      if (present(pmulti_state)) then
         pmulti_state=>multi_state
      else
         deallocate(multi_state)
      end if


      if (has_scalar_field(state(1),"Porosity")) then
         call insert(packed_state,extract_scalar_field(state(1),"Porosity"),"Porosity")
      else
         call allocate(porosity,element_mesh,"Porosity")
         call set(porosity,1.0)
         call insert(packed_state,porosity,"Porosity")
         call deallocate(porosity)
      end if

      if (has_scalar_field(state(1),"Permeability")) then
         call allocate(permeability,element_mesh,"Permeability",&
              dim=[mesh_dim(position),mesh_dim(position)])
         call zero(permeability)
         sfield=>extract_scalar_field(state(1),"Permeability")
         do idim=1,mesh_dim(position)
            call set(permeability,idim,idim,sfield)
         end do
         call insert(packed_state,permeability,"Permeability")
         call deallocate(permeability)
      else if (has_vector_field(state(1),"Permeability")) then
         call allocate(permeability,element_mesh,"Permeability",&
              dim=[mesh_dim(position),mesh_dim(position)])
         call zero(permeability)
         vfield=>extract_vector_field(state(1),"Permeability")
         call set(permeability,vfield)
         call insert(packed_state,permeability,"Permeability")
         call deallocate(permeability)
      else if (has_vector_field(state(1),"Permeability")) then
         call insert(packed_state,extract_tensor_field(state(1),"Permeability"),"Permeability")
      else
         call allocate(permeability,element_mesh,"Permeability")
         call zero(permeability)
         call insert(packed_state,permeability,"Permeability")
         call deallocate(permeability)
      end if

      contains

        subroutine insert_components_in_multi_state(ms,s)
          type(state_type), dimension(:) :: ms
          type(state_type) :: s

          integer :: i
          character(len=FIELD_NAME_LEN) :: full_name, short_name

          do i=1,size(ms)
             full_name="ComponentMassFractionPhase"//int2str(i)
             short_name="ComponentMassFraction"
             if (has_scalar_field(s,trim(full_name))) then
                call insert(ms(i),extract_scalar_field(s,trim(full_name)),short_name)
             end if
          end do

          do i=1,size(ms)
             full_name="ComponentDensityPhase"//int2str(i)
             short_name="ComponentDensity"
             if (has_scalar_field(s,trim(full_name))) then
                call insert(ms(i),extract_scalar_field(s,trim(full_name)),short_name)
             end if
          end do

        end subroutine insert_components_in_multi_state

        subroutine unpack_multicomponent(mstate,mcstate)
          type(state_type) :: mstate
          type(state_type), dimension(:) :: mcstate

          integer :: icomp

          type(tensor_field), pointer :: component, density
          type(tensor_field) :: vfield


          component=>extract_tensor_field(mstate,"PackedComponentMassFraction")
          density=>extract_tensor_field(mstate,"PackedComponentDensity")

          do icomp=1,ncomp
             call allocate(vfield,component%mesh,"ComponentMassFraction",field_type=FIELD_TYPE_DEFERRED,dim=[1,nphase])
             vfield%option_path=component%option_path
             deallocate(vfield%val)
             vfield%val=>component%val(icomp:icomp,:,:)
             vfield%wrapped=.true.
             call insert(mcstate(icomp),vfield,vfield%name)
             call deallocate(vfield)


             call allocate(vfield,component%mesh,"ComponentDensity",field_type=FiELD_TYPE_DEFERRED,dim=[1,nphase])
             vfield%option_path=component%option_path
             deallocate(vfield%val)
             vfield%val=>component%val(icomp:icomp,:,:)
             vfield%wrapped=.true.
             call insert(mcstate(icomp),vfield,vfield%name)
             call deallocate(vfield)
             
          end do

        end subroutine unpack_multicomponent
          

        subroutine unpack_multiphase(mstate,mpstate)
          type(state_type) :: mstate
          type(state_type), dimension(:) :: mpstate
          integer :: index, iphase

          type(tensor_field), pointer :: tfield
          type(tensor_field) :: mp_tfield


          do index=1,size(mstate%tensor_fields)
             tfield=>extract_tensor_field(mstate,index)
             if (tfield%name(:6)=="Packed") then
                do iphase=1,nphase
                   call allocate(mp_tfield,tfield%mesh,tfield%name(7:),field_type=FiELD_TYPE_DEFERRED,dim=[tfield%dim(1),1])

                   mp_tfield%val=>tfield%val(:,iphase:iphase,:)
                   mp_tfield%updated=>tfield%updated
                   mp_tfield%wrapped=.true.
                   call insert(mpstate(iphase),mp_tfield,mp_tfield%name)
                   call deallocate(mp_tfield)
                end do
             end if
          end do

        END subroutine unpack_multiphase
        

        subroutine add_new_memory(mstate,sfield,name)
          type(state_type) :: mstate
          type(scalar_field) :: sfield
          character (len=*) :: name
          type(scalar_field) :: sfield2

          call allocate(sfield2,sfield%mesh,name)
          call insert(mstate,sfield2,name)
          call deallocate(sfield2)

        end subroutine add_new_memory
          

        subroutine insert_sfield(mstate,name,ncomp,nphase,nmesh,&
             add_source, add_absorption)
          type(state_type), intent(inout) :: mstate
          character(len=*), intent(in) :: name
          type(mesh_type),optional, target :: nmesh
          integer, intent(in) :: ncomp,nphase
          logical, optional, intent(in) :: add_source, add_absorption

          type(scalar_field), pointer :: nfield
          type(mesh_type), pointer :: lmesh
          type(tensor_field) :: mfield
          integer :: stat
          logical :: ladd_source, ladd_absorption

          if (present(add_source)) then
             ladd_source=add_source
          else
             ladd_source=.false.
          end if

          if (present(add_absorption)) then
             ladd_absorption=add_absorption
          else
             ladd_absorption=.false.
          end if

          nfield=>extract_scalar_field(state(1),name,stat)
          if (stat/=0) then
             nfield=>extract_scalar_field(state(1),"Pressure",stat)
          end if

          if (present(nmesh)) then
             lmesh=> nmesh
          else
             lmesh=>nfield%mesh
          end if

          call allocate(mfield,lmesh,"Packed"//name,dim=[ncomp,nphase])
          call insert(mstate,mfield,"Packed"//name)
          call deallocate(mfield)
          call allocate(mfield,lmesh,"PackedOld"//name,dim=[ncomp,nphase])
          call insert(mstate,mfield,"PackedOld"//name)
          call deallocate(mfield)
          call allocate(mfield,lmesh,"PackedIterated"//name,dim=[ncomp,nphase])
          call insert(mstate,mfield,"PackedIterated"//name)
          call deallocate(mfield)
          
          if (ladd_source) then
             call allocate(mfield,lmesh,"Packed"//trim(name)//"Source",&
                  dim=[ncomp,nphase])
             call zero(mfield)
             call insert(mstate,mfield,"Packed"//trim(name)//"Source")
             call deallocate(mfield)
          end if

          if (ladd_absorption) then
             call allocate(mfield,lmesh,"Packed"//trim(name)//"Absorption",&
                  dim=[ncomp,nphase])
             call zero(mfield)
             call insert(mstate,mfield,"Packed"//trim(name)//"Absorption")
             call deallocate(mfield)
          end if

        end subroutine insert_sfield

        subroutine insert_vfield(mstate,name,nmesh, add_source, add_absorption)
          type(state_type), intent(inout) :: mstate
          character(len=*), intent(in) :: name
          type(mesh_type), optional, target :: nmesh
          logical, optional, intent(in) :: add_source, add_absorption   

          type(vector_field), pointer :: nfield
          type(mesh_type), pointer :: lmesh
          type(tensor_field)  :: mfield
          logical :: ladd_source, ladd_absorption

          nfield=>extract_vector_field(state(1),name)
          if (present(nmesh)) then
             lmesh=> nmesh
          else
             lmesh=>nfield%mesh
          end if

          if (present(add_source)) then
             ladd_source=add_source
          else
             ladd_source=.false.
          end if

          if (present(add_absorption)) then
             ladd_absorption=add_absorption
          else
             ladd_absorption=.false.
          end if

          
          call allocate(mfield,lmesh,"Packed"//name,dim=[ndim,nphase])
          call insert(mstate,mfield,"Packed"//name)
          call deallocate(mfield)
          call allocate(mfield,lmesh,"PackedOld"//name,dim=[ndim,nphase])
          call insert(mstate,mfield,"PackedOld"//name)
          call deallocate(mfield)
          call allocate(mfield,lmesh,"PackedIterated"//name,dim=[ndim,nphase])
          call insert(mstate,mfield,"PackedIterated"//name)
          call deallocate(mfield)

          if (ladd_source) then
             call allocate(mfield,lmesh,"Packed"//trim(name)//"Source",&
                  dim=[ndim,nphase])
             call zero(mfield)
             call insert(mstate,mfield,"Packed"//trim(name)//"Source")
             call deallocate(mfield)
          end if

          if (ladd_absorption) then
             call allocate(mfield,lmesh,"Packed"//trim(name)//"Absorption",&
                  dim=[ndim,nphase])
             call zero(mfield)
             call insert(mstate,mfield,"Packed"//trim(name)//"Absorption")
             call deallocate(mfield)
          end if

        end subroutine insert_vfield

        subroutine unpack_sfield(nstate,mstate,name,icomp, iphase,free)
          type(state_type), intent(inout) :: nstate, mstate
          character(len=*) :: name
          integer :: icomp,iphase, stat
          logical, optional :: free 

          type(scalar_field), pointer :: nfield
          type(tensor_field), pointer :: mfield
          logical lfree


          if (present(free)) then
             lfree=free
          else
             lfree=.true.
          end if

          mfield=>extract_tensor_field(mstate,"Packed"//name)

          nfield=>extract_scalar_field(nstate,name,stat)
          if (stat==0) then
             mfield%val(icomp,iphase,:)=nfield%val(:)
             if (icomp==1 .and. iphase == 1) then
                mfield%option_path=nfield%option_path
             end if
             if(lfree) deallocate(nfield%val)
             nfield%val=>mfield%val(icomp,iphase,:)
             nfield%val_stride=ncomp*nphase
             nfield%wrapped=.true.
          end if          

        end subroutine unpack_sfield

        subroutine unpack_vfield(nstate,mstate,name,iphase,free)
          type(state_type), intent(inout) :: nstate, mstate
          character(len=*) :: name
          integer :: iphase, stat
          logical, optional :: free 

          
          type(vector_field), pointer :: nfield
          type(tensor_field), pointer :: mfield
          logical lfree

          if (present(free)) then
             lfree=free
          else
             lfree=.true.
          end if

          nfield=>extract_vector_field(nstate,name,stat)
          mfield=>extract_tensor_field(mstate,"Packed"//name)

          if (stat==0) then
             mfield%val(:,iphase,:)=nfield%val(:,1:node_count(mfield))
             if (icomp==1 .and. iphase == 1) then
                mfield%option_path=nfield%option_path
             end if
             if (lfree) then
#ifdef HAVE_MEMORY_STATS
                call register_deallocation("vector_field", "real", &
                     size(nfield%val), name=nfield%name)
#endif
                deallocate(nfield%val)
             end if
             nfield%val=>mfield%val(:,iphase,:)
             nfield%wrapped=.true.
          end if


        end subroutine unpack_vfield

        subroutine unpack_component_sfield(st,mst,name,ic,prefix)
          type(state_type) :: st, mst
          character (len=*) :: name
          integer :: ic
          integer :: ip
          character (len=*), optional :: prefix


          type(scalar_field), pointer :: nfield,pnfield
          type(tensor_field), pointer :: mfield
          logical :: free
          

          mfield=>extract_tensor_field(mst,"Packed"//name)

          do ip=1,nphase
             if ( has_scalar_field(st,name//"Phase"//int2str(ip)) ) then

                nfield=>extract_scalar_field(st,name//"Phase"//int2str(ip))
                if (present(prefix)) then
                   pnfield=>extract_scalar_field(st,name(len(prefix)+1:)//"Phase"//int2str(ip))
                   free=check_paired(pnfield,nfield)
                else
                   free=.true.
                end if

                
               
                
                mfield%val(ic,ip,:)=nfield%val(:)
                if (ic==1) then
                   mfield%option_path=nfield%option_path
                end if
                if (free) deallocate(nfield%val)
                nfield%val=>mfield%val(ic,ip,:)
                nfield%val_stride=ncomp*nphase
                nfield%wrapped=.true.
             end if
          end do
                
        end subroutine unpack_component_sfield

        subroutine allocate_multiphase_scalar_bcs(s,ms,name)

          type(state_type), intent(inout) :: s
          type(state_type), dimension(:,:), intent(inout) :: ms
          character( len=*) :: name

          type(tensor_field), pointer :: mfield
          type(scalar_field), pointer :: sfield
          type(tensor_boundary_condition) :: tbc
          type(scalar_boundary_condition), pointer ::  bc
          type(tensor_boundary_condition), dimension(:), pointer :: temp

          integer :: iphase,icomp,stat, n, nbc

          nbc=0
          mfield=>extract_tensor_field(s,"Packed"//name)

          allocate(mfield%bc)
          
          do iphase=1,mfield%dim(2)
             do icomp=1,mfield%dim(1)
                sfield=>extract_scalar_field( ms(icomp,iphase),trim(name),stat)
                if (stat/= 0 ) cycle

                allocate(tbc%applies(mfield%dim(1),mfield%dim(2)))
                tbc%applies= .false.
                tbc%applies(icomp,iphase)=.true.
               
                do n=1,get_boundary_condition_count(sfield)

                  bc=>sfield%bc%boundary_condition(n)
                  
                  tbc%name=bc%name
                  tbc%type=bc%type
                  tbc%surface_element_list=>bc%surface_element_list
                  tbc%surface_node_list=>bc%surface_node_list
                  tbc%surface_mesh=>bc%surface_mesh
                  call incref(tbc%surface_mesh)
                  tbc%scalar_surface_fields=>bc%surface_fields
                  
                  nbc=nbc+1
                  if (nbc>1) then
                     temp=>mfield%bc%boundary_condition
                     allocate(mfield%bc%boundary_condition(nbc))
                     mfield%bc%boundary_condition(1:nbc-1)=temp
                     deallocate(temp)
                     mfield%bc%boundary_condition(nbc)=tbc
                  else
                     allocate(mfield%bc%boundary_condition(nbc))
                     mfield%bc%boundary_condition(nbc)=tbc
                  end if
               end do
            end do
         end do
               

       end subroutine allocate_multiphase_scalar_bcs


        subroutine allocate_multiphase_vector_bcs(s,ms,name)

          type(state_type), intent(inout) :: s
          type(state_type), dimension(:,:), intent(inout) :: ms
          character( len=*) :: name

          type(tensor_field), pointer :: mfield
          type(vector_field), pointer :: vfield
          type(tensor_boundary_condition) :: tbc
          type(vector_boundary_condition), pointer ::  bc
          type(tensor_boundary_condition), dimension(:), pointer :: temp

          integer :: iphase, stat, n, nbc

          nbc=0
          mfield=>extract_tensor_field(s,"Packed"//name)

          allocate(mfield%bc)
          
          do iphase=1,mfield%dim(2)
             vfield=>extract_vector_field( ms(1,iphase),trim(name),stat)
             if (stat/= 0 ) cycle

             allocate(tbc%applies(mfield%dim(1),mfield%dim(2)))
             tbc%applies= .false.
             tbc%applies(:,iphase)=.true.
               
             do n=1,get_boundary_condition_count(vfield)

                bc=>vfield%bc%boundary_condition(n)
                  
                tbc%name=bc%name
                tbc%type=bc%type
                tbc%surface_element_list=>bc%surface_element_list
                tbc%surface_node_list=>bc%surface_node_list
                tbc%surface_mesh=>bc%surface_mesh
                call incref(tbc%surface_mesh)
                tbc%vector_surface_fields=>bc%surface_fields
                  
                nbc=nbc+1
                if (nbc>1) then
                   temp=>mfield%bc%boundary_condition
                   allocate(mfield%bc%boundary_condition(nbc))
                   mfield%bc%boundary_condition(1:nbc-1)=temp
                   deallocate(temp)
                   mfield%bc%boundary_condition(nbc)=tbc
                else
                   allocate(mfield%bc%boundary_condition(nbc))
                   mfield%bc%boundary_condition(nbc)=tbc
                end if
             end do
          end do

       end subroutine allocate_multiphase_vector_bcs


subroutine allocate_multicomponent_scalar_bcs(s,ms,name)

          type(state_type), dimension(:), intent(inout) :: s
          type(state_type), dimension(:,:), intent(inout) :: ms
          character( len=*) :: name

          type(tensor_field), pointer :: mfield
          type(scalar_field), pointer :: sfield
          type(tensor_boundary_condition) :: tbc
          type(scalar_boundary_condition), pointer ::  bc
          type(tensor_boundary_condition), dimension(:), pointer :: temp

          integer :: iphase,icomp,stat, n, nbc


          do icomp=1,size(s)
             nbc=0
             mfield=>extract_tensor_field(s(icomp),name)
             allocate(mfield%bc)
          
             do iphase=1,mfield%dim(2)
                sfield=>extract_scalar_field( ms(icomp,iphase),trim(name),stat)

                if (stat/= 0 ) cycle

                allocate(tbc%applies(mfield%dim(1),mfield%dim(2)))
                tbc%applies= .false.
                tbc%applies(icomp,iphase)=.true.
                
                do n=1,get_boundary_condition_count(sfield)
                   
                   bc=>sfield%bc%boundary_condition(n)
                   
                  tbc%name=bc%name
                  tbc%type=bc%type
                  tbc%surface_element_list=>bc%surface_element_list
                  tbc%surface_node_list=>bc%surface_node_list
                  tbc%surface_mesh=>bc%surface_mesh
                  call incref(tbc%surface_mesh)
                  tbc%scalar_surface_fields=>bc%surface_fields
                  
                  nbc=nbc+1
                  if (nbc>1) then
                     temp=>mfield%bc%boundary_condition
                     allocate(mfield%bc%boundary_condition(nbc))
                     mfield%bc%boundary_condition(1:nbc-1)=temp
                     deallocate(temp)
                     mfield%bc%boundary_condition(nbc)=tbc
                  else
                     allocate(mfield%bc%boundary_condition(nbc))
                     mfield%bc%boundary_condition(nbc)=tbc
                  end if
               end do
            end do
         end do
               

       end subroutine allocate_multicomponent_scalar_bcs

       function check_paired(sfield1,sfield2) result(unpaired)
         type(scalar_field) :: sfield1,sfield2
         logical unpaired

         unpaired = .not. associated(sfield2%val,sfield1%val)

       end function check_paired

       function check_vpaired(vfield1,vfield2) result(unpaired)
         type(vector_field) :: vfield1,vfield2
         logical unpaired

         unpaired = .not. associated(vfield2%val,vfield1%val)

       end function check_vpaired


      end subroutine pack_multistate

      function wrap_as_tensor(field) result(tfield)

        type(scalar_field), intent(inout) :: field  
        type(tensor_field), pointer :: tfield

        allocate(tfield)
        call allocate(tfield,field%mesh,name=field%name,dim=[1,1])

        tfield%val(1,1,:)=field%val
        deallocate(tfield%updated)
        tfield%updated=field%updated
        deallocate(field%val)
        deallocate(field%updated)
        field%val=>tfield%val(1,1,:)
        field%updated=>tfield%updated
        tfield%option_path=field%option_path

      end function wrap_as_tensor

      function as_vector(tfield,dim,slice) result(vfield)

        type(tensor_field), intent(inout) :: tfield  
        integer, intent(in) :: dim
        integer, intent(in), optional :: slice
        
        type(vector_field)  :: vfield 
        integer :: lslice

        if (present(slice)) then
           lslice=slice
        else
           lslice=1
        end if

        vfield%name=tfield%name
        vfield%mesh=tfield%mesh
        vfield%option_path=tfield%option_path
        vfield%dim=tfield%dim(dim)
        select case(dim)
        case(1)
           vfield%val=>tfield%val(:,lslice,:)
        case(2)
           vfield%val=>tfield%val(lslice,:,:)
        end select
        vfield%wrapped=.true.

      end function as_vector


      subroutine finalise_multistate(packed_state,multiphase_state,&
           multicomponent_state)

        type(state_type) :: packed_state
        type(state_type), dimension(:), pointer :: multiphase_state, multicomponent_state


        call deallocate(multiphase_state)
        deallocate(multiphase_state)
        call deallocate(multicomponent_state)
        deallocate(multicomponent_state)
        call deallocate(packed_state)

      end subroutine finalise_multistate



    subroutine add_dependant_fields_to_tensor_from_state(infield,state,scalar_field_names,&
         vector_field_names,tensor_field_names)

      !  Convenience subroutine to add a bunch of fields as dependants of infield in one call.

      type(tensor_field) :: infield
      type(state_type) :: state
      character (len=*) , dimension(:), optional :: scalar_field_names, vector_field_names, &
           tensor_field_names

      integer :: i
      type(scalar_field), pointer :: sfield
      type(vector_field), pointer :: vfield
      type(tensor_field), pointer :: tfield

      if (present(scalar_field_names)) then
         do i=1,size(scalar_field_names)
            sfield=>extract_scalar_field(state,trim(scalar_field_names(i)))
            call add_dependant_field(infield,sfield)
         end do
      end if

      if (present(vector_field_names)) then
         do i=1,size(vector_field_names)
            vfield=>extract_vector_field(state,trim(vector_field_names(i)))
            call add_dependant_field(infield,vfield)
         end do
      end if

      if (present(tensor_field_names)) then
         do i=1,size(tensor_field_names)
            tfield=>extract_tensor_field(state,tensor_field_names(i))
            call add_dependant_field(infield,tfield)
         end do
      end if

    end subroutine add_dependant_fields_to_tensor_from_state


    subroutine Adaptive_NonLinear(packed_state, reference_field, its,&
        Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs,order)
        !This subroutine either store variables before the nonlinear timeloop starts, or checks
        !how the nonlinear iterations are going and depending on that increase the timestep
        !or decreases the timestep and repeats that timestep
        Implicit none
        type(state_type), intent(inout) :: packed_state!, backup_state
        real, dimension(:,:,:), allocatable, intent(inout) :: reference_field
        logical, intent(inout) :: Repeat_time_step, ExitNonLinearLoop
        logical, intent(in) :: nonLinearAdaptTs
        integer, intent(in) :: its, order
        !Local variables
        real :: dt
        real, parameter :: check_sat_threshold = 1d-6
        real, dimension(:), pointer :: pressure
        real, dimension(:,:), pointer :: phasevolumefraction
        real, dimension(:,:,:), pointer :: velocity

!        real, dimension(:,:,:), pointer :: tVar, tVar_it
!        real, dimension(:,:), pointer :: vVar, vVar_it
!        real, dimension(:,:), pointer :: sVar, sVar_it

        !Variables for automatic non-linear iterations
        real :: tolerance_between_non_linear, initial_dt, min_ts, max_ts, increase_ts_switch, decrease_ts_switch
        !Variables for adaptive time stepping based on non-linear iterations
        real :: increaseFactor, decreaseFactor, ts_ref_val, s_gc, s_or, acctim
        integer :: variable_selection, i, NonLinearIteration

        !First of all, check if the user wants to do something
        call get_option( '/timestepping/nonlinear_iterations/nonlinear_iterations_automatic', tolerance_between_non_linear, default = -1. )
        if (tolerance_between_non_linear<0) return
        call get_option( '/timestepping/nonlinear_iterations', NonLinearIteration, default = 3 )
          !Get data from diamond. Despite this is slow, as it is done in the outest loop, it should not affect the performance.
         !Variable to check how good nonlinear iterations are going 1 (Pressure), 2 (Velocity), 3 (Saturation)
        call get_option( '/timestepping/nonlinear_iterations/nonlinear_iterations_automatic/adaptive_timestep_nonlinear', &
        variable_selection, default = 3)
        call get_option( '/timestepping/nonlinear_iterations/nonlinear_iterations_automatic/adaptive_timestep_nonlinear/increase_factor', &
        increaseFactor, default = 1.2 )
        call get_option( '/timestepping/nonlinear_iterations/nonlinear_iterations_automatic/adaptive_timestep_nonlinear/decrease_factor', &
        decreaseFactor, default = 2. )
        call get_option( '/timestepping/nonlinear_iterations/nonlinear_iterations_automatic/adaptive_timestep_nonlinear/max_timestep', &
        max_ts, default = huge(min_ts) )
        call get_option( '/timestepping/nonlinear_iterations/nonlinear_iterations_automatic/adaptive_timestep_nonlinear/min_timestep', &
        min_ts, default = -1. )
        call get_option( '/timestepping/nonlinear_iterations/nonlinear_iterations_automatic/adaptive_timestep_nonlinear/increase_ts_switch', &
        increase_ts_switch, default = 1d-3 )
        call get_option( '/timestepping/nonlinear_iterations/nonlinear_iterations_automatic/adaptive_timestep_nonlinear/decrease_ts_switch', &
        decrease_ts_switch, default = 1d-1 )
        !Get irresidual water saturation and irreducible oil to repeat a timestep if the saturation goes out of these values
        call get_option("/material_phase[0]/multiphase_properties/immobile_fraction", &
        s_gc, default=0.0)
        call get_option("/material_phase[1]/multiphase_properties/immobile_fraction", &
        s_or, default=0.0)
        !Get time step
        call get_option( '/timestepping/timestep', initial_dt )
        dt = initial_dt
        !By default the minimum time-steps is ten orders smaller than the initial timestep
        if(min_ts<0) min_ts = initial_dt * 1d-10


        select case (order)
            case (1)!Store or get from backup
                !If we do not have adaptive time stepping then there is nothing to backup
                if (.not.nonLinearAdaptTs) return
                !we either store the data or we recover it if repeting a timestep
                !Procedure to repeat time-steps
                !If  Repeat_time_step then we recover values, else we store them
                call copy_packed_new_to_iterated(packed_state, Repeat_time_step)
            case (2)!Calculate and store reference_field
                if (its==1) then
                    !Store variable to check afterwards
                    call get_var_from_packed_state(packed_state, velocity = velocity, pressure = pressure,&
                    phasevolumefraction = phasevolumefraction)
                    if (allocated(reference_field)) deallocate(reference_field)
                    select case (variable_selection)
                        case (1)
                            allocate (reference_field(1,1,size(pressure,1) ))
                            reference_field(1,1, :) = pressure
                        case (2)
                            allocate (reference_field(size(velocity,1),size(velocity,2),size(velocity,3) ))
                            reference_field(:,:,:) = velocity
                        case default
                            allocate (reference_field(1,size(phasevolumefraction,1),size(phasevolumefraction,2) ))
                            reference_field(1,:,:) = phasevolumefraction
                        end select
                end if

            case default!Check how is the process going on and decide
                 !If Automatic_NonLinerIterations then we compare the variation of the a property from one time step to the next one
                ExitNonLinearLoop = .false.
                Repeat_time_step = .false.
                if (its > 1 ) then

                    call get_var_from_packed_state(packed_state, velocity = velocity, pressure = pressure,&
                    phasevolumefraction = phasevolumefraction)

                    select case (variable_selection)
                        case (1)
                            ts_ref_val = maxval(abs(reference_field(1,1,:)-pressure))
                        case (2)
                             ts_ref_val = maxval(abs(reference_field-velocity))
                        case default
                            ts_ref_val = maxval(abs(reference_field(1,:,:)-phasevolumefraction))
                    end select
                    !If only non-linear iterations
                    if (.not.nonLinearAdaptTs) then
                        !Automatic non-linear iteration checking
                        if (ts_ref_val < tolerance_between_non_linear)  ExitNonLinearLoop = .true.
                        return
                    end if
                    !Check that saturation is between bounds, works for two phases only!!
                    if (have_option(  '/timestepping/nonlinear_iterations/&
                            &nonlinear_iterations_automatic/adaptive_timestep_nonlinear/keep_sat_bounds') )  then
                        Repeat_time_step = (maxval(s_gc-phasevolumefraction(1,:))>check_sat_threshold&
                        .or.maxval(s_or-phasevolumefraction(2,:))>check_sat_threshold)
                    end if

                     !Increase Ts section
                    if ((ts_ref_val < increase_ts_switch .and.dt*increaseFactor<max_ts).and..not.Repeat_time_step) then
                        call get_option( '/timestepping/timestep', dt )
                        dt = dt * increaseFactor
                        call set_option( '/timestepping/timestep', dt )
                        print *, "Time step increased to:", dt
                        ExitNonLinearLoop = .true.
                        return
                    end if

!                    !Exit loop section
!                    if ((ts_ref_val < tolerance_between_non_linear).and..not.Repeat_time_step) then
!                        ExitNonLinearLoop = .true.
!                        return
!                    end if

                    !Decrease Ts section only if we have done at least the 70% of the  nonLinearIterations
                    if ((ts_ref_val > decrease_ts_switch.or.repeat_time_step) &
                    .and.its>=int(0.70*NonLinearIteration)) then

                        if ( dt / decreaseFactor < min_ts) then
                            !Do not decrease
                            Repeat_time_step = .false.
                            ExitNonLinearLoop = .true.
                            deallocate(reference_field)
                            return
                        end if

                        !Decrease time step, reset the time and repeat!
                        call get_option( '/timestepping/timestep', dt )
                        call get_option( '/timestepping/current_time', acctim )
                        acctim = acctim - dt
                        call set_option( '/timestepping/current_time', acctim )
                        dt = dt / decreaseFactor
                        call set_option( '/timestepping/timestep', dt )
                        print *, "Time step decreased to:", dt
                        Repeat_time_step = .true.
                        ExitNonLinearLoop = .true.
                    end if
                end if

        end select

    end subroutine Adaptive_NonLinear


   subroutine copy_packed_new_to_iterated(packed_state, viceversa)
    !Values from packed_state are stored in iterated unless viceversa is true, in that case
    !the iterated values are moved to the new values
     type(state_type), intent(inout) :: packed_state
     logical, intent(in) :: viceversa

     type(scalar_field), pointer :: sfield, nsfield
     type(vector_field), pointer :: vfield, nvfield
     type(tensor_field), pointer :: tfield, ntfield

     integer :: i

     do i=1,size(packed_state%scalar_fields)
        sfield=>packed_state%scalar_fields(i)%ptr
        if (sfield%name(1:14)=="PackedIterated") then
            nsfield=>extract_scalar_field(packed_state,"Packed"//sfield%name(15:))
            if (viceversa) then
                nsfield%val = sfield%val
            else
                sfield%val=nsfield%val
            end if
        end if
    end do

    do i=1,size(packed_state%vector_fields)
        vfield=>packed_state%vector_fields(i)%ptr
         if (vfield%name(1:14)=="PackedIterated") then
             nvfield=>extract_vector_field(packed_state,"Packed"//vfield%name(15:))
             if (viceversa) then
                 nvfield%val = vfield%val
             else
                 vfield%val=nvfield%val
             end if
         end if
     end do

     do i=1,size(packed_state%tensor_fields)
         tfield=>packed_state%tensor_fields(i)%ptr
         if (tfield%name(1:14)=="PackedIterated") then
             ntfield=>extract_tensor_field(packed_state,"Packed"//tfield%name(15:))
             if (viceversa) then
                 ntfield%val  = tfield%val
             else
                 tfield%val=ntfield%val
             end if
         end if
     end do

     sfield=>extract_scalar_field(packed_state,"OldFEPressure")
     nsfield=>extract_scalar_field(packed_state,"FEPressure")
     sfield%val=nsfield%val

     sfield=>extract_scalar_field(packed_state,"OldCVPressure")
     nsfield=>extract_scalar_field(packed_state,"CVPressure")
     sfield%val=nsfield%val

   end subroutine copy_packed_new_to_iterated

    subroutine get_var_from_packed_state(packed_state,FEDensity,&
    OldFEDensity,IteratedFEDensity,Density,OldDensity,IteratedDensity,PhaseVolumeFraction,&
    OldPhaseVolumeFraction,IteratedPhaseVolumeFraction, Velocity, OldVelocity, IteratedVelocity, &
    FEPhaseVolumeFraction, OldFEPhaseVolumeFraction, IteratedFEPhaseVolumeFraction,&
    NonlinearVelocity, OldNonlinearVelocity,IteratedNonlinearVelocity, ComponentDensity, &
    OldComponentDensity, IteratedComponentDensity,ComponentMassFraction, OldComponentMassFraction,&
    Temperature,OldTemperature, IteratedTemperature,FETemperature, OldFETemperature, IteratedFETemperature,&
    IteratedComponentMassFraction, FEComponentDensity, OldFEComponentDensity, IteratedFEComponentDensity,&
    FEComponentMassFraction, OldFEComponentMassFraction, IteratedFEComponentMassFraction,&
    Pressure,FEPressure, OldFEPressure, CVPressure,OldCVPressure,&
    Coordinate, VelocityCoordinate,PressureCoordinate,MaterialCoordinate  )
        !This subroutine returns a pointer to the desired values of a variable stored in packed state
        !All the input variables (but packed_stated) are pointers following the structure of the *_ALL variables
        !and also all of them are optional, hence you can obtaine whichever you want
        !######################EXAMPLE OF USAGE OF THIS SUBROUTINE:#####################################
        !If we want to get the velocity and the phasevolumefraction one should proceed this way:
        !Define variables:
        !real, dimension(:,:,:), pointer :: Velocity_pointer
        !real, dimension(:,:), pointer :: PhaseVolumeFraction_pointer
        !Assign the pointers
        !call get_var_from_packed_state(packed_state, Velocity = Velocity_pointer, PhaseVolumeFraction = PhaseVolumeFraction_pointer)
        !
        ! In this way we only have to introduce the name of the variables we want to get from packed_state
        !########################################################################################
        implicit none
        type(state_type), intent(inout) :: packed_state
        real, optional, dimension(:,:,:), pointer :: Velocity, OldVelocity, IteratedVelocity, NonlinearVelocity, OldNonlinearVelocity,&
        IteratedNonlinearVelocity,ComponentDensity, OldComponentDensity, IteratedComponentDensity,ComponentMassFraction, OldComponentMassFraction,&
        IteratedComponentMassFraction, FEComponentDensity, OldFEComponentDensity, IteratedFEComponentDensity, FEComponentMassFraction, &
        OldFEComponentMassFraction, IteratedFEComponentMassFraction
        real, optional, dimension(:,:), pointer :: FEDensity, OldFEDensity, IteratedFEDensity, Density,&
        OldDensity,IteratedDensity,PhaseVolumeFraction,OldPhaseVolumeFraction,IteratedPhaseVolumeFraction,&
        Temperature, OldTemperature, IteratedTemperature, FETemperature, OldFETemperature, IteratedFETemperature,&
        Coordinate, VelocityCoordinate,PressureCoordinate,MaterialCoordinate, &
        FEPhaseVolumeFraction, OldFEPhaseVolumeFraction, IteratedFEPhaseVolumeFraction
        real, optional, dimension(:), pointer ::Pressure,FEPressure, OldFEPressure, CVPressure,OldCVPressure
        !Local variables
        type(scalar_field), pointer :: sfield
        type(vector_field), pointer :: vfield
        type(tensor_field), pointer :: tfield

        !Scalar stored
        if (present(Pressure)) then
            sfield => extract_scalar_field( packed_state, "Pressure" )
            Pressure =>  sfield%val(:)
        end if
        if (present(FEPressure)) then
            sfield => extract_scalar_field( packed_state, "FEPressure" )
            FEPressure =>  sfield%val(:)
        end if
        if (present(OldFEPressure)) then
            sfield => extract_scalar_field( packed_state, "OldFEPressure" )
            OldFEPressure =>  sfield%val(:)
        end if
        if (present(CVPressure)) then
            sfield => extract_scalar_field( packed_state, "CVPressure" )
            CVPressure =>  sfield%val(:)
        end if
        if (present(OldCVPressure)) then
            sfield => extract_scalar_field( packed_state, "OldCVPressure" )
            OldCVPressure =>  sfield%val(:)
        end if

        !Vectors stored
        if (present(Coordinate)) then
            vfield => extract_vector_field( packed_state, "Coordinate" )
            Coordinate =>  vfield%val(:,:)
        end if
        if (present(VelocityCoordinate)) then
            vfield => extract_vector_field( packed_state, "VelocityCoordinate" )
            VelocityCoordinate =>  vfield%val(:,:)
        end if
        if (present(PressureCoordinate)) then
            vfield => extract_vector_field( packed_state, "PressureCoordinate" )
            PressureCoordinate =>  vfield%val(:,:)
        end if
        if (present(MaterialCoordinate)) then
            vfield => extract_vector_field( packed_state, "MaterialCoordinate" )
            MaterialCoordinate =>  vfield%val(:,:)
        end if

        !Tensors stored
        if (present(FEDensity)) then
            tfield => extract_tensor_field( packed_state, "PackedFEDensity" )
            FEDensity =>  tfield%val(1,:,:)
        end if
        if (present(OldFEDensity)) then
            tfield => extract_tensor_field( packed_state, "PackedOldFEDensity" )
            OldFEDensity =>  tfield%val(1,:,:)
        end if
        if (present(IteratedFEDensity)) then
            tfield => extract_tensor_field( packed_state, "PackedIteratedFEDensity" )
            IteratedFEDensity => tfield%val(1,:,:)
        end if
        if (present(Density)) then
            tfield => extract_tensor_field( packed_state, "PackedDensity" )
            Density => tfield%val(1,:,:)
        end if
        if (present(OldDensity)) then
            tfield => extract_tensor_field( packed_state, "PackedOldDensity" )
            OldDensity => tfield%val(1,:,:)
        end if
        if (present(IteratedDensity)) then
            tfield => extract_tensor_field( packed_state, "PackedIteratedDensity" )
            IteratedDensity =>  tfield%val(1,:,:)
        end if
        if (present(PhaseVolumeFraction)) then
            tfield => extract_tensor_field( packed_state, "PackedPhaseVolumeFraction" )
            PhaseVolumeFraction =>  tfield%val(1,:,:)
        end if
        if (present(OldPhaseVolumeFraction)) then
            tfield => extract_tensor_field( packed_state, "PackedOldPhaseVolumeFraction" )
            OldPhaseVolumeFraction =>  tfield%val(1,:,:)
        end if
        if (present(IteratedPhaseVolumeFraction)) then
            tfield => extract_tensor_field( packed_state, "PackedIteratedPhaseVolumeFraction" )
            IteratedPhaseVolumeFraction =>  tfield%val(1,:,:)
        end if
        if (present(FEPhaseVolumeFraction)) then
            tfield => extract_tensor_field( packed_state, "PackedFEPhaseVolumeFraction" )
            FEPhaseVolumeFraction =>  tfield%val(1,:,:)
        end if
        if (present(OldFEPhaseVolumeFraction)) then
            tfield => extract_tensor_field( packed_state, "PackedOldFEPhaseVolumeFraction" )
            OldFEPhaseVolumeFraction =>  tfield%val(1,:,:)
        end if
        if (present(IteratedFEPhaseVolumeFraction)) then
            tfield => extract_tensor_field( packed_state, "PackedIteratedFEPhaseVolumeFraction" )
            IteratedFEPhaseVolumeFraction =>  tfield%val(1,:,:)
        end if
        if (present(Temperature)) then
            tfield => extract_tensor_field( packed_state, "PackedTemperature" )
            Temperature =>  tfield%val(1,:,:)
        end if
        if (present(OldTemperature)) then
            tfield => extract_tensor_field( packed_state, "PackedOldTemperature" )
            OldTemperature =>  tfield%val(1,:,:)
        end if
        if (present(IteratedTemperature)) then
            tfield => extract_tensor_field( packed_state, "PackedIteratedTemperature" )
            IteratedTemperature =>  tfield%val(1,:,:)
        end if
         if (present(FETemperature)) then
            tfield => extract_tensor_field( packed_state, "PackedFETemperature" )
            FETemperature =>  tfield%val(1,:,:)
        end if
         if (present(OldFETemperature)) then
            tfield => extract_tensor_field( packed_state, "PackedOldFETemperature" )
            OldFETemperature =>  tfield%val(1,:,:)
        end if
        if (present(IteratedFETemperature)) then
            tfield => extract_tensor_field( packed_state, "PackedIteratedFETemperature" )
            IteratedFETemperature =>  tfield%val(1,:,:)
        end if
        if (present(Velocity)) then
            tfield => extract_tensor_field( packed_state, "PackedVelocity" )
            Velocity => tfield%val(:,:,:)
        end if
        if (present(OldVelocity)) then
            tfield => extract_tensor_field( packed_state, "PackedOldVelocity" )
            OldVelocity => tfield%val(:,:,:)
        end if
        if (present(IteratedVelocity)) then
            tfield => extract_tensor_field( packed_state, "PackedIteratedVelocity" )
            IteratedVelocity => tfield%val(:,:,:)
        end if
        if (present(NonlinearVelocity)) then
            tfield => extract_tensor_field( packed_state, "PackedNonlinearVelocity" )
            NonlinearVelocity => tfield%val(:,:,:)
        end if
        if (present(OldNonlinearVelocity)) then
            tfield => extract_tensor_field( packed_state, "PackedOldNonlinearVelocity" )
            OldNonlinearVelocity => tfield%val(:,:,:)
        end if
        if (present(IteratedNonlinearVelocity)) then
            tfield => extract_tensor_field( packed_state, "PackedIteratedNonlinearVelocity" )
            IteratedNonlinearVelocity => tfield%val(:,:,:)
        end if
        if (present(ComponentDensity)) then
            tfield => extract_tensor_field( packed_state, "PackedComponentDensity" )
            ComponentDensity => tfield%val(:,:,:)
        end if
        if (present(OldComponentDensity)) then
            tfield => extract_tensor_field( packed_state, "PackedOldComponentDensity" )
            OldComponentDensity => tfield%val(:,:,:)
        end if
        if (present(IteratedComponentDensity)) then
            tfield => extract_tensor_field( packed_state, "PackedIteratedComponentDensity" )
            IteratedComponentDensity => tfield%val(:,:,:)
        end if
        if (present(ComponentMassFraction)) then
            tfield => extract_tensor_field( packed_state, "PackedComponentMassFraction" )
            ComponentMassFraction => tfield%val(:,:,:)
        end if
        if (present(OldComponentMassFraction)) then
            tfield => extract_tensor_field( packed_state, "PackedOldComponentMassFraction" )
            OldComponentMassFraction => tfield%val(:,:,:)
        end if
        if (present(IteratedComponentMassFraction)) then
            tfield => extract_tensor_field( packed_state, "PackedIteratedComponentMassFraction" )
            IteratedComponentMassFraction => tfield%val(:,:,:)
        end if
        if (present(FEComponentDensity)) then
            tfield => extract_tensor_field( packed_state, "PackedFEComponentDensity" )
            FEComponentDensity => tfield%val(:,:,:)
        end if
        if (present(OldFEComponentDensity)) then
            tfield => extract_tensor_field( packed_state, "PackedOldFEComponentDensity" )
            OldFEComponentDensity => tfield%val(:,:,:)
        end if
        if (present(IteratedFEComponentDensity)) then
            tfield => extract_tensor_field( packed_state, "PackedIteratedFEComponentDensity" )
            IteratedFEComponentDensity => tfield%val(:,:,:)
        end if
        if (present(FEComponentMassFraction)) then
            tfield => extract_tensor_field( packed_state, "PackedFEComponentMassFraction" )
            FEComponentMassFraction => tfield%val(:,:,:)
        end if
        if (present(OldFEComponentMassFraction)) then
            tfield => extract_tensor_field( packed_state, "PackedOldFEComponentMassFraction" )
            OldFEComponentMassFraction => tfield%val(:,:,:)
        end if
        if (present(IteratedFEComponentMassFraction)) then
            tfield => extract_tensor_field( packed_state, "PackedIteratedFEComponentMassFraction" )
            IteratedFEComponentMassFraction => tfield%val(:,:,:)
        end if



    end subroutine get_var_from_packed_state


  end module Copy_Outof_State

