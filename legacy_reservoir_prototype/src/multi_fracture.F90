
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

module multiphase_fractures

  use quadrature
  use elements
  use sparse_tools
  use fields
  use fefields
  use state_module
  use copy_outof_state
  use spud
  use global_parameters, only: option_path_len, field_name_len
  use futils, only: int2str
  use solvers
  use implicit_solids
  use FLDebug

  implicit none

#ifdef USING_FEMDEM
  ! variable name convention:
  ! _r : ring, _v : volume, _vc : volume coarse

  interface
     subroutine y2d_allocate_femdem( string, &
          &     nodes_r, elements_r, edges_r, &
          &     nodes_v, elements_v, edges_v )
       character( len = * ), intent( in ) :: string
       integer, intent( out ) :: nodes_r, elements_r, edges_r, &
            &                    nodes_v, elements_v, edges_v
     end subroutine y2d_allocate_femdem
  end interface

  interface
     subroutine y2d_populate_femdem( ele1_r, ele2_r, ele3_r, &
          &                          face1_r, face2_r, x_r, y_r, &
          &                          p1, p2, p3, p4, &
          &                          ele1_v, ele2_v, ele3_v, &
          &                          face1_v, face2_v, x_v, y_v )
       integer, dimension( * ), intent( out ) :: ele1_r, ele2_r, ele3_r, &
            &                                    ele1_v, ele2_v, ele3_v, &
            &                                    face1_r, face2_r, face1_v, face2_v
       real, dimension( * ), intent( out ) :: x_r, y_r, x_v, y_v, p1, p2, p3, p4
     end subroutine y2d_populate_femdem
  end interface

  interface
     subroutine y2dfemdem( string, dt, rho, p, u, v, u_s, v_s )
       character( len = * ), intent( in ) :: string
       real, intent( in ) :: dt
       real, dimension( * ), intent( in ) :: rho, p, u, v
       real, dimension( * ), intent( out ) :: u_s, v_s
     end subroutine y2dfemdem
  end interface
#endif

  type( vector_field ), save :: positions_r, positions_v, positions_vc
  type( tensor_field ), save :: permeability_r
  character( len = FIELD_NAME_LEN ), save :: femdem_mesh_name


  private
  public :: femdem, blasting

contains

  subroutine femdem( states, totele, cv_nonods, u_nonods, ndim, nphase, cv_nloc, &
       &             cv_ndgln, dt, rho, p, u, v, absorption, perm, porosity )

    implicit none

    integer, intent( in ) :: totele, cv_nonods, u_nonods, nphase, ndim, cv_nloc
    integer, dimension( totele * cv_nloc ), intent( in ) :: cv_ndgln
    real, intent( in ) :: dt
    real, dimension( nphase * cv_nonods ), intent( in ) :: rho
    real, dimension( nphase * u_nonods ), intent( in ) :: u, v
    real, dimension( cv_nonods ), intent( in ) :: p
    type( state_type ), dimension( : ), intent( in ) :: states

    real, dimension( totele * cv_nloc, ndim * nphase, ndim * nphase ), intent( inout ) :: absorption
    real, dimension( totele, ndim, ndim ), intent( inout ) :: perm
    real, dimension( totele ), intent( inout ) :: porosity

    real, dimension( : ), allocatable :: vf, rho_r, p_r, u_r, v_r, u_s, v_s
    integer :: r_nonods

    ! read in ring and solid volume meshes
    ! and simplify the volume mesh
    call initialise_femdem

    ! calculate volume fraction using the coarse mesh
    allocate( vf( totele ) ) ; vf = 0.
    call calculate_volume_fraction( states, totele, vf )

    ! calculate absorption coefficient
    call calculate_absorption( totele, cv_nloc, cv_nonods, ndim, &
         &                     nphase, cv_ndgln, rho, vf, dt, absorption )

    ! calculate porosity
    call calculate_phi( states, totele, vf, porosity )

    ! calculate permeability and porosity
    call calculate_perm( states, totele, ndim, vf, perm )


    r_nonods = node_count( positions_r )
    allocate( rho_r( r_nonods), p_r( r_nonods ), &
         &    u_r( r_nonods ), v_r( r_nonods ), u_s( r_nonods ), v_s( r_nonods ) )
    call interpolate_fields( states, nphase, u_nonods, cv_nonods, r_nonods, &
         &                   cv_nloc, totele, cv_ndgln, &
         &                   rho, p, u, v, rho_r, p_r, u_r, v_r )
    call y2dfemdem( trim( femdem_mesh_name ) // char( 0 ), dt, rho_r, p_r, u_r, v_r, u_s, v_s )

    ! deallocate
    call deallocate_femdem
    deallocate( vf )
    deallocate( rho_r, p_r, u_r, v_r, u_s, v_s )

    return
  end subroutine femdem



  subroutine blasting( states, totele, cv_nonods, u_nonods, ndim, nphase, cv_nloc, &
       &               cv_ndgln, dt, rho, p, u, v, absorption, perm, porosity )

    implicit none

    integer, intent( in ) :: totele, cv_nonods, u_nonods, nphase, ndim, cv_nloc
    integer, dimension( totele * cv_nloc ), intent( in ) :: cv_ndgln
    real, intent( in ) :: dt
    real, dimension( nphase * cv_nonods ), intent( in ) :: rho
    real, dimension( nphase * u_nonods ), intent( in ) :: u, v
    real, dimension( cv_nonods ), intent( in ) :: p
    type( state_type ), dimension( : ), intent( in ) :: states

    real, dimension( totele * cv_nloc, ndim * nphase, ndim * nphase ), intent( inout ) :: absorption
    real, dimension( totele, ndim, ndim ), intent( inout ) :: perm
    real, dimension( totele ), intent( inout ) :: porosity

    real, dimension( : ), allocatable :: rho_r, p_r, u_r, v_r, u_s, v_s 
    integer :: r_nonods

    ! read in ring and solid volume meshes
    ! and simplify the volume mesh
    call initialise_femdem

    r_nonods = node_count( positions_r )
    allocate( rho_r( r_nonods ), p_r( r_nonods ), &
         &    u_r( r_nonods ), v_r( r_nonods ), u_s( r_nonods ), v_s( r_nonods ) )
    call interpolate_fields( states, nphase, u_nonods, cv_nonods, r_nonods, &
         &                   cv_nloc, totele, cv_ndgln, &
         &                   rho, p, u, v, rho_r, p_r, u_r, v_r )
    call y2dfemdem( trim( femdem_mesh_name ) // char( 0 ), dt, rho_r, p_r, u_r, v_r, u_s, v_s )

    ! deallocate
    call deallocate_femdem
    deallocate( rho_r, p_r, u_r, v_r, u_s, v_s )

    return
  end subroutine blasting

  !----------------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------

  subroutine initialise_femdem

    implicit none

    integer :: i, loc, sloc
    integer :: ndim, nodes_r, elements_r, edges_r, &
         &     nodes_v, elements_v, edges_v
    integer, dimension( : ), allocatable :: ele1_r, ele2_r, ele3_r, &
         &                                  ele1_v, ele2_v, ele3_v, &
         &                                  face1_r, face2_r, face1_v, face2_v
    real, dimension( : ), allocatable :: x_r, y_r, x_v, y_v, p1, p2, p3, p4 
    type( quadrature_type ) :: quad
    type( element_type ) :: shape
    integer, dimension( : ), allocatable :: sndglno_r, boundary_ids_r, &
         &                                  sndglno_v, boundary_ids_v
    integer :: quad_degree, poly_degree, continuity
    type( mesh_type ) :: mesh_r, mesh_v, mesh_r_p0

    ewrite(3,*) "inside initialise_femdem"

    !call get_option( "/femdem/name", femdem_mesh_name )
    femdem_mesh_name = "rock.msh"

    call get_option( "/geometry/quadrature/degree", quad_degree )
    call get_option( "/geometry/dimension", ndim )

    if ( ndim == 2 ) then
       loc = 3 ; sloc= 2
    else if ( ndim == 3 ) then
       loc = 4 ; sloc= 3
       FLAbort( "Fracture modelling is supported for 2D only." )
    end if

    call y2d_allocate_femdem( trim( femdem_mesh_name ) // char( 0 ), &
         nodes_r, elements_r, edges_r, nodes_v, elements_v, edges_v )

    ewrite(3,*) "nodes_r, elements_r, edges_r, nodes_v, elements_v, edges_v", &
                 nodes_r, elements_r, edges_r, nodes_v, elements_v, edges_v

    allocate( ele1_r( elements_r ), ele2_r( elements_r ), ele3_r( elements_r ) )
    allocate( face1_r( edges_r ), face2_r( edges_r ) )

    allocate( ele1_v( elements_v ), ele2_v( elements_v ), ele3_v( elements_v ) )
    allocate( face1_v( edges_v ), face2_v( edges_v ) )

    allocate( x_r( nodes_r ), y_r( nodes_r ) )
    allocate( x_v( nodes_v ), y_v( nodes_v ) )

    allocate( p1( elements_r ), p2( elements_r ) )
    allocate( p3( elements_r ), p4( elements_r ) )

    call y2d_populate_femdem( ele1_r, ele2_r, ele3_r, &
         face1_r, face2_r, x_r, y_r, p1, p2, p3, p4, &
         ele1_v, ele2_v, ele3_v, face1_v, face2_v, x_v, y_v )

    quad = make_quadrature( loc, ndim, degree = quad_degree )
    shape = make_element_shape( loc, ndim, 1, quad )

    ! create the ring mesh
    call allocate( mesh_r, nodes_r, elements_r, shape, name="CoordinateMesh" )
    call allocate( positions_r, ndim, mesh_r, name="Coordinate" )

    positions_r%val( 1, : ) = x_r
    positions_r%val( 2, : ) = y_r

    do i = 1, elements_r
       positions_r%mesh%ndglno( (i-1)*loc+1 : i*loc ) = &
            (/ ele1_r(i)+1, ele2_r(i)+1, ele3_r(i)+1 /)
    end do

    allocate( sndglno_r( edges_r * sloc ) ) ; sndglno_r = 0
    allocate( boundary_ids_r( edges_r ) ) ; boundary_ids_r = 666

    do i = 1, edges_r
       sndglno_r( (i-1)*sloc+1 : i*sloc ) = &
            (/ face1_r(i)+1, face2_r(i)+1 /)
    end do

    call add_faces( positions_r%mesh, &
         sndgln = sndglno_r, &
         boundary_ids = boundary_ids_r )

    positions_r%dim = ndim

    deallocate( boundary_ids_r, sndglno_r )
    call deallocate( mesh_r )

    ! create the volume mesh
    call allocate( mesh_v, nodes_v, elements_v, shape, name="CoordinateMesh" )
    call allocate( positions_v, ndim, mesh_v, name="Coordinate" )
    positions_v%val( 1, : ) = x_v
    positions_v%val( 2, : ) = y_v

    do i = 1, elements_v
       positions_v%mesh%ndglno( (i-1)*loc+1 : i*loc ) = &
            (/ ele1_v(i)+1, ele2_v(i)+1, ele3_v(i)+1 /)
    end do

    allocate( sndglno_v( edges_v * sloc ) ) ; sndglno_v = 0
    allocate( boundary_ids_v( edges_v ) ) ; boundary_ids_v = 666

    do i = 1, edges_v
       sndglno_v( (i-1)*sloc+1 : i*sloc ) = &
            (/ face1_v(i)+1, face2_v(i)+1 /)
    end do

    call add_faces( positions_v%mesh, &
         sndgln = sndglno_v, &
         boundary_ids = boundary_ids_v )

    positions_v%dim = ndim

    deallocate( boundary_ids_v, sndglno_v )
    call deallocate( mesh_v )

    call deallocate_element( shape )

    ! coarsen the volume mesh
    call coarsen_mesh_2d( positions_v, positions_vc )

    ! create the ring p0 mesh
    poly_degree = 0 ; continuity = -1
    shape = make_element_shape( loc, ndim, poly_degree, quad )
    mesh_r_p0 = make_mesh( positions_r%mesh, shape, continuity, "P0DG" )

    ! store ring permeability
    call allocate( permeability_r, mesh_r_p0, name = "Permeability" )
    call zero( permeability_r )

    call set_all( permeability_r, 1, 1, p1 )
    call set_all( permeability_r, 1, 2, p2 )
    call set_all( permeability_r, 2, 1, p3 )
    call set_all( permeability_r, 2, 2, p4 )

    ! deallocate
    call deallocate( mesh_r_p0 )
    call deallocate_element( shape )
    call deallocate( quad )

    deallocate( ele1_r, ele2_r, ele3_r, face1_r, face2_r, &
         &      ele1_v, ele2_v, ele3_v, face1_v, face2_v, &
         &      x_r, y_r, x_v, y_v, p1, p2, p3, p4 )

    ewrite(3,*) "leaving initialise_femdem"

    return
  end subroutine initialise_femdem

  !----------------------------------------------------------------------------------------------------------

  subroutine calculate_volume_fraction( states, totele, vf )

    implicit none

    integer, intent( in ) :: totele
    type( state_type ), dimension( : ), intent( in ) :: states
    real, dimension( totele ), intent( inout ) :: vf

    type( state_type ) :: alg_ext, alg_fl
    type( mesh_type ), pointer :: fl_mesh, p0_fl_mesh
    type( vector_field ) :: fl_positions
    type( scalar_field ) :: volume_fraction

    ewrite(3,*) "inside calculate_volume_fraction"

    call insert( alg_ext, positions_vc % mesh, "Mesh" )
    call insert( alg_ext, positions_vc, "Coordinate" )

    fl_mesh => extract_mesh( states( 1 ), "CoordinateMesh" )
    fl_positions = extract_vector_field( states(1), "Coordinate" )

    call insert( alg_fl, fl_mesh, "Mesh" )
    call insert( alg_fl, fl_positions, "Coordinate" )

    p0_fl_mesh => extract_mesh( states( 1 ), "P0DG" )

    ! volume fraction, i.e. porosity...
    call allocate( volume_fraction, p0_fl_mesh, "VolumeFraction" )
    call zero( volume_fraction )

    ewrite(3,*) "...interpolating"

    ! interpolate
    call interpolation_galerkin_femdem( alg_ext, alg_fl, field = volume_fraction )

    vf = volume_fraction % val

    ! ensure vf is between 0. and 1.
    call bound_volume_fraction( vf )

    ! now deallocate
    call deallocate( fl_positions )
    call deallocate( volume_fraction )
    call deallocate( alg_fl )
    call deallocate( alg_ext )

    ewrite(3,*) "leaving calculate_volume_fraction"

    return
  end subroutine calculate_volume_fraction

  !----------------------------------------------------------------------------------------------------------

  subroutine calculate_absorption( totele, cv_nloc, cv_nonods, ndim, &
       &                           nphase, cv_ndgln, rho, vf, dt, absorption )

    implicit none

    integer, intent( in ) :: totele, cv_nloc, cv_nonods, ndim, nphase
    integer, dimension( totele*cv_nloc ), intent( in ) :: cv_ndgln
    real, intent( in ) :: dt
    real, dimension( totele ), intent( in ) :: vf
    real, dimension( nphase * cv_nonods), intent( in ) :: rho

    real, dimension( totele * cv_nloc, ndim * nphase, ndim * nphase ), intent( inout ) :: absorption

    integer :: ele, iloc, mi, ci, iphase, idim, idx
    real :: sigma

    ewrite(3,*) "inside calculate_absorption"

    do ele = 1, totele
       do iloc = 1, cv_nloc

          mi = ( ele - 1 ) * cv_nloc + iloc
          ci = cv_ndgln( ( ele - 1 ) * cv_nloc + iloc )

          do iphase = 1, nphase
             sigma = vf( ele ) * rho( ci + ( iphase - 1 ) * cv_nonods ) / dt
             do idim = 1, ndim
                idx = idim + ( iphase - 1 ) * ndim
                absorption( mi, idx, idx ) = absorption( mi, idx, idx ) + sigma
             end do
          end do

       end do
    end do

    ewrite(3,*) "leaving calculate_absorption"

    return
  end subroutine calculate_absorption

  !----------------------------------------------------------------------------------------------------------

  subroutine calculate_phi( states, totele, vf, porosity )

    implicit none

    type( state_type ), dimension( : ), intent( in ) :: states
    integer, intent( in ) :: totele
    real, dimension( totele ), intent( in ) :: vf

    real, dimension( totele ), intent( inout ) :: porosity

    real :: phi_rock
    type( scalar_field ), pointer :: phi_matrix
    integer :: ele

    ewrite(3,*) "inside calculate_phi"

    ! rock mass porosity
    phi_rock = 0.
    ! matrix (background) porosity
    phi_matrix => extract_scalar_field( states( 1 ), "Porosity" )

    porosity = ( 1. - vf ) * phi_matrix % val + vf * phi_rock 

    ! bound...
    do ele = 1, totele
       porosity( ele ) = max( 0., min( 1., porosity( ele ) ) )
    end do

    ewrite(3,*) "leaving calculate_phi"

    return
  end subroutine calculate_phi

  !----------------------------------------------------------------------------------------------------------

  subroutine calculate_perm( states, totele, ndim, vf, perm )

    implicit none

    integer, intent( in ) :: totele, ndim
    type( state_type ), dimension( : ), intent( in ) :: states
    real, dimension( totele ), intent( in ) :: vf

    real, dimension( totele, ndim, ndim ), intent( inout ) :: perm

    type( state_type ) :: alg_ext, alg_fl
    type( mesh_type ), pointer :: fl_mesh, p0_fl_mesh
    type( vector_field ) :: fl_positions
    type( scalar_field ) :: rvf

    type( tensor_field ), pointer :: permeability_bg_t
    type( scalar_field ), pointer :: permeability_bg_s

    type( scalar_field ) :: field_fl_p11, field_fl_p12, field_fl_p21, field_fl_p22, &
         &                  field_ext_p11, field_ext_p12, field_ext_p21, field_ext_p22

    real :: ring, solid, bg
    real, dimension( :, : ), allocatable :: permeability_v
    real, dimension( :, :, : ), allocatable :: permeability_bg, permeability_rl 

    character( len = OPTION_PATH_LEN ) :: &
         path = "/tmp/galerkin_projection/continuous"
    integer :: ele, idim, jdim

    real, parameter :: tol = 1.e-10

    ewrite(3,*) "inside calculate_perm"

    call insert( alg_ext, positions_r%mesh, "Mesh" )
    call insert( alg_ext, positions_r, "Coordinate" )

    fl_mesh => extract_mesh( states( 1 ), "CoordinateMesh" )
    fl_positions = extract_vector_field( states( 1 ), "Coordinate" )

    call insert( alg_fl, fl_mesh, "Mesh" )
    call insert( alg_fl, fl_positions, "Coordinate" )

    call set_solver_options( path, &
         ksptype = "cg", &
         pctype = "mg", &
         rtol = 1.e-10, &
         atol = 0., &
         max_its = 10000 )

    path = "/tmp"

    p0_fl_mesh => extract_mesh( states( 1 ), "P0DG" )

    ewrite(3,*) "...generating fluids state"

    ! this is the permeability on the fluidity mesh
    call allocate( field_fl_p11, p0_fl_mesh, "Permeability11" )
    call zero( field_fl_p11 )
    call insert( alg_fl, field_fl_p11, "Permeability11" )
    field_fl_p11 % option_path = path

    call allocate( field_fl_p12, p0_fl_mesh, "Permeability12" )
    call zero( field_fl_p12 )
    call insert( alg_fl, field_fl_p12, "Permeability12" )
    field_fl_p12 % option_path = path

    call allocate( field_fl_p21, p0_fl_mesh, "Permeability21" )
    call zero( field_fl_p21 )
    call insert( alg_fl, field_fl_p21, "Permeability21" )
    field_fl_p21 % option_path = path

    call allocate( field_fl_p22, p0_fl_mesh, "Permeability22" )
    call zero( field_fl_p22 )
    call insert( alg_fl, field_fl_p22, "Permeability22" )
    field_fl_p22 % option_path = path

    ewrite(3,*) "...generating femdem/ring state"

    ! this is the permeability on the solid mesh
    field_ext_p11 = extract_scalar_field( permeability_r, 1, 1 )
    call insert( alg_ext, field_ext_p11, "Permeability11" )

    field_ext_p12 = extract_scalar_field( permeability_r, 1, 2 )
    call insert( alg_ext, field_ext_p12, "Permeability12" )

    field_ext_p21 = extract_scalar_field( permeability_r, 2, 1 )
    call insert( alg_ext, field_ext_p21, "Permeability21" )

    field_ext_p22 = extract_scalar_field( permeability_r, 2, 2 )
    call insert( alg_ext, field_ext_p22, "Permeability22" )

    ! volume fraction - this is the ring
    call allocate( rvf, p0_fl_mesh, "VolumeFraction" )
    call zero( rvf )

    ewrite(3,*) "...interpolating"

    ! interpolate
    call interpolation_galerkin_femdem( alg_ext, alg_fl, field = rvf )

    ! bound ring volume fraction
    call bound_volume_fraction( rvf % val )

    allocate( permeability_rl(totele, ndim, ndim) ) ; permeability_rl = 0.
    permeability_rl(:, 1, 1) = field_fl_p11 % val
    permeability_rl(:, 1, 2) = field_fl_p12 % val
    permeability_rl(:, 2, 1) = field_fl_p21 % val
    permeability_rl(:, 2, 2) = field_fl_p22 % val

    ! background and solid permeabilities
    allocate( permeability_bg(totele, ndim, ndim) ) ; permeability_bg = 0.
    if( have_option( "/porous_media/scalar_field::Permeability" ) ) then

       permeability_bg_s => extract_scalar_field( states( 1 ), "Permeability" )
       do ele = 1, element_count( permeability_bg_s ) 
          forall( idim = 1 : ndim ) permeability_bg( ele, idim, idim ) = permeability_bg_s % val( ele )
       end do

    elseif( have_option( "/porous_media/tensor_field::Permeability" ) ) then

       permeability_bg_t => extract_tensor_field( states( 1 ), "Permeability" )
       path = "/porous_media/tensor_field::Permeability"
       call Extract_TensorFields_Outof_State( states, 1, &
            permeability_bg_t, path, permeability_bg )
    end if

    allocate( permeability_v(ndim, ndim) ) ; permeability_v = 0.
    forall( idim = 1:ndim ) permeability_v( idim, idim ) = 0.3 

    do ele = 1, totele

       ! solid, ring and background volume fractions
       solid = vf( ele )
       ring = rvf % val( ele )
       ! if the ring permeability is zero we are not in a fracture
       ! so permeability is a combination of the solid and matrix 
       ! (background) permeabilities
       if ( permeability_rl( ele, 1, 1 ) <= tol ) ring = 0.
       bg = max( 0., 1. - solid - ring )

       forall( idim = 1:ndim, jdim = 1:ndim ) &
            perm( ele, idim, jdim ) = bg * permeability_bg( ele, idim, jdim ) &
            &                       + solid * permeability_v( idim, jdim ) &
            &                       + permeability_rl( ele, idim, jdim )

    end do

    deallocate( permeability_bg, permeability_v, permeability_rl )

    call deallocate( fl_positions )
    call deallocate( rvf )

    call deallocate( field_fl_p22 )
    call deallocate( field_fl_p21 )
    call deallocate( field_fl_p12 )
    call deallocate( field_fl_p11 )

    call deallocate( alg_fl )
    call deallocate( alg_ext )

    ewrite(3,*) "leaving calculate_perm"

    return
  end subroutine calculate_perm

  !----------------------------------------------------------------------------------------------------------


  subroutine interpolate_fields( states, nphase, u_nonods, cv_nonods, r_nonods, &
       &                         cv_nloc, totele, cv_ndgln, &
       &                         rho, p, u, v, rho_r, p_r, u_r, v_r )

    implicit none

    type( state_type ), dimension( : ), intent( in ) :: states
    integer, intent( in ) :: nphase, u_nonods, cv_nonods, r_nonods, cv_nloc 
    integer, dimension( totele * cv_nloc ), intent( in ) :: cv_ndgln
    real, dimension( nphase * u_nonods ), intent( in ) :: u, v
  
    real, dimension( cv_nonods ), intent( in ) :: rho
    real, dimension( cv_nonods ), intent( in ) :: p
    real, dimension( r_nonods ), intent( inout ) :: u_r, v_r, p_r, rho_r

    type( mesh_type ), pointer :: fl_mesh, u_mesh
    type( scalar_field ) :: field_fl_rho, field_fl_p, field_fl_u, field_fl_v, &
         &                  field_ext_rho, field_ext_p, field_ext_u, field_ext_v, &
         &                  u_dg, v_dg
    type( vector_field ) :: fl_positions
    type( scalar_field ), pointer :: pressure
    type( vector_field ), pointer :: velocity
    type( state_type ) :: alg_ext, alg_fl
    real, dimension( : ), allocatable :: u_tmp, v_tmp
    integer, dimension( : ), pointer :: fl_ele_nodes
    integer :: ele, totele, u_nloc, nlev, i, j, k, number_nodes
    character(len=option_path_len) :: vel_element_type
    logical :: is_overlapping
    character( len = OPTION_PATH_LEN ) :: &
         path = "/tmp/galerkin_projection/continuous"

    u_r = 0. ; v_r = 0. ; p_r = 0. ; rho_r = 0.

    fl_mesh => extract_mesh( states( 1 ), "CoordinateMesh" )
    fl_positions = extract_vector_field( states( 1 ), "Coordinate" )

    u_mesh => extract_mesh( states( 1 ), "VelocityMesh" )

    totele = ele_count( fl_mesh )

    call set_solver_options( path, &
         ksptype = "cg", &
         pctype = "mg", &
         rtol = 1.e-10, &
         atol = 0., &
         max_its = 10000 )

    path = "/tmp"

    call allocate( field_fl_rho, fl_mesh, "Density" )
    call zero( field_fl_rho )

    call allocate( field_fl_p, fl_mesh, "Pressure" )
    call zero( field_fl_p )

    call allocate( field_fl_u, fl_mesh, "Velocity1" )
    call zero( field_fl_u )

    call allocate( field_fl_v, fl_mesh, "Velocity2" )
    call zero( field_fl_v )

    ! deal with pressure
    if ( cv_nloc == 6  ) then
       ! linearise pressure for p2
       do ele = 1, totele
          fl_ele_nodes => ele_nodes( fl_mesh, ele )
          field_fl_p % val( fl_ele_nodes( 1 ) ) = p( cv_ndgln( ( ele - 1 ) * cv_nloc + 1 ) )
          field_fl_p % val( fl_ele_nodes( 2 ) ) = p( cv_ndgln( ( ele - 1 ) * cv_nloc + 3 ) )
          field_fl_p % val( fl_ele_nodes( 3 ) ) = p( cv_ndgln( ( ele - 1 ) * cv_nloc + 6 ) )
       end do
    else
       ! just copy memory for p1
       field_fl_p % val = p
    end if

    ! deal with density
    if ( cv_nloc == 6  ) then
       ! linearise density for p2
       do ele = 1, totele
          fl_ele_nodes => ele_nodes( fl_mesh, ele )
          field_fl_rho % val( fl_ele_nodes( 1 ) ) = rho( cv_ndgln( ( ele - 1 ) * cv_nloc + 1 ) )
          field_fl_rho % val( fl_ele_nodes( 2 ) ) = rho( cv_ndgln( ( ele - 1 ) * cv_nloc + 3 ) )
          field_fl_rho % val( fl_ele_nodes( 3 ) ) = rho( cv_ndgln( ( ele - 1 ) * cv_nloc + 6 ) )
       end do
    else
       ! just copy memory for p1
       field_fl_rho % val = rho
    end if

    ! deal with velocity  
    velocity => extract_vector_field( states( 1 ), "Velocity" )
    number_nodes = node_count( velocity )
    u_nloc = ele_loc( velocity, 1 )

    allocate( u_tmp( number_nodes * nphase ) ) ; u_tmp = 0.
    allocate( v_tmp( number_nodes * nphase ) ) ; v_tmp = 0.

    call get_option('/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
         vel_element_type)
    is_overlapping = .false.
    if ( trim( vel_element_type ) == "overlapping" ) is_overlapping = .true.

    if ( is_overlapping ) then
       pressure => extract_scalar_field( states( 1 ), "Pressure" )
       nlev = ele_loc( pressure, 1 )

       do k = 1, nphase
          do i = 1, totele
             do j = 1, nlev
                u_tmp ( (k-1)*number_nodes + (i-1)*u_nloc + 1 : (k-1)*number_nodes + i*u_nloc ) = &
                     u_tmp ( (k-1)*number_nodes + (i-1)*u_nloc + 1 : (k-1)*number_nodes + i*u_nloc ) + &
                     u ( (k-1)*number_nodes*nlev + (i-1)*u_nloc*nlev + (j-1)*u_nloc+ 1 : & 
                     (k-1)*number_nodes*nlev + (i-1)*u_nloc*nlev + j*u_nloc )
                v_tmp ( (k-1)*number_nodes + (i-1)*u_nloc + 1 : (k-1)*number_nodes + i*u_nloc ) = &
                     v_tmp ( (k-1)*number_nodes + (i-1)*u_nloc + 1 : (k-1)*number_nodes + i*u_nloc ) + &
                     v ( (k-1)*number_nodes*nlev + (i-1)*u_nloc*nlev + (j-1)*u_nloc+ 1 : & 
                     (k-1)*number_nodes*nlev + (i-1)*u_nloc*nlev + j*u_nloc )
             end do
          end do
       end do

       u_tmp = u_tmp / nlev
       v_tmp = v_tmp / nlev

    else
       u_tmp = u
       v_tmp = v
    end if

    call allocate( u_dg, u_mesh, "u_dg" )
    call zero( u_dg )
    call allocate( v_dg, u_mesh, "v_dg" )
    call zero( v_dg )
    u_dg % val = u_tmp ; v_dg % val = v_tmp 

    call project_field( u_dg, field_fl_u, fl_positions )
    call project_field( v_dg, field_fl_v, fl_positions )

    ! fluidity state
    call insert( alg_fl, fl_mesh, "Mesh" )
    call insert( alg_fl, fl_positions, "Coordinate" )

    call insert(alg_fl, field_fl_rho, "Density")
    call insert(alg_fl, field_fl_p, "Pressure")
    call insert(alg_fl, field_fl_u, "Velocity1")
    call insert(alg_fl, field_fl_v, "Velocity2")

    ! ring state
    call insert( alg_ext, positions_r%mesh, "Mesh" )
    call insert( alg_ext, positions_r, "Coordinate" )

    call allocate( field_ext_rho, fl_mesh, "Density" )
    call zero( field_ext_rho )
    field_ext_rho % option_path = path
    call insert( alg_ext, field_ext_rho, "Density" )

    call allocate( field_ext_p, fl_mesh, "Pressure" )
    call zero( field_ext_p )
    field_ext_p % option_path = path
    call insert( alg_ext, field_ext_p, "Pressure" )

    call allocate( field_ext_u, fl_mesh, "Velocity1" )
    call zero( field_ext_u )
    field_ext_u % option_path = path
    call insert( alg_ext, field_ext_u, "Velocity1" )

    call allocate( field_ext_v, fl_mesh, "Velocity2" )
    call zero( field_ext_v )
    field_ext_v % option_path = path
    call insert( alg_ext, field_ext_v, "Velocity2" )

    ewrite(3,*) "...interpolating"

    ! interpolate
    call interpolation_galerkin_femdem( alg_fl, alg_ext, femdem_out = .true. )

    ! copy memory
    p_r = field_ext_p % val
    u_r = field_ext_u % val
    v_r = field_ext_v % val

    call deallocate( fl_positions )

    call deallocate( field_fl_rho )
    call deallocate( field_fl_p )
    call deallocate( field_fl_u )
    call deallocate( field_fl_v )

    call deallocate( field_ext_rho )
    call deallocate( field_ext_p )
    call deallocate( field_ext_u )
    call deallocate( field_ext_v )

    call deallocate( alg_fl )
    call deallocate( alg_ext )

    call deallocate( u_dg )
    call deallocate( v_dg )

    deallocate( u_tmp, v_tmp )

    return
  end subroutine interpolate_fields


  !----------------------------------------------------------------------------------------------------------

  subroutine coarsen_mesh_2d( f, c )

    implicit none

    type( vector_field ), intent( in ) :: f
    type( vector_field ), intent( out ) :: c

    integer :: snloc, nloc, nele_s, nele, ele, nnodes_s, &
         &     ele2, siloc, sjloc, iloc, jloc, st, mv, &
         &     edge, nedge, bcs, nbcs, count, ne, i, j, &
         &     ndim, quad_degree
    integer, dimension( : ), allocatable :: ndglno, sndglno, boundary_ids, &
         &                                  ndglno_s, tmp, ele_nodes, &
         &                                  ele2_nodes, bc, n, nodes_s
    real :: area
    real, dimension( : ), allocatable :: X_s, Y_s
    logical :: on_the_wall, delete_ele2
    type( quadrature_type ) :: quad
    type( element_type ) :: shape
    type( mesh_type ) :: mesh

    logical, parameter :: coarsen_mesh = .false.

    nedge = 3 ! triangles
    nloc  = 3 ! 3 nodes per element
    snloc = 2 ! 2 nodes per edge
    ndim  = 2

    nele = element_count( f )
    nbcs = surface_element_count( f )

    if ( coarsen_mesh ) then

       allocate( ndglno( nloc * nele ), sndglno( snloc * nbcs ) ) 
       allocate( ele_nodes( nloc ), ele2_nodes( nloc ) )
       allocate( n( snloc ), bc(snloc ) )

       ndglno = f % mesh % ndglno
       call getsndgln( f % mesh, sndglno)

       nele_s = nele ! number of elements in the simplified mesh

       do ele = 1, nele

          ele_nodes = ndglno( ( ele - 1 ) * nloc + 1 : ele * nloc )

          ! check if this element has been deleted
          if ( all( ele_nodes < 0 ) ) cycle

          do edge = 1, nedge

             n( 1 ) = ele_nodes( edge )
             if ( edge < nedge ) then
                n( 2 ) = ele_nodes( edge+1 )
             else
                n( 2 ) = ele_nodes( 1 )
             end if

             on_the_wall = .false.
             do bcs = 1, nbcs
                bc = sndglno( (bcs-1)*snloc+1 : bcs*snloc )

                ! on_the_wall = .true. if all nodes on the boundary
                if ( all( n == bc ) ) then
                   on_the_wall = .true.
                   exit
                end if
                ! if one node is on the wall 
                ! figure out which one it is
                st = -666
                if ( any( n == bc ) .and. .not.on_the_wall ) then
                   do siloc = 1, snloc
                      do sjloc = 1, snloc
                         if ( n( siloc ) == bc( sjloc ) ) st = n( siloc )
                         if ( n( siloc ) /= bc( sjloc ) ) mv = n( siloc )
                      end do
                   end do
                end if
             end do

             if ( any( sndglno == mv ) ) on_the_wall = .true.

             if ( .not.on_the_wall ) then

                if ( st < 0 ) then
                   st = n( 1 )
                   mv = n( 2 )
                end if

                do ele2 = 1, nele

                   ele2_nodes = ndglno( ( ele2 - 1 ) * nloc + 1 : ele2 * nloc )

                   ! check if this element has been deleted
                   if ( all( ele2_nodes < 0 ) ) cycle

                   do iloc = 1, nloc
                      if ( ele2_nodes( iloc ) == mv ) &
                           ndglno( ( ele2 - 1 ) * nloc + iloc ) = st
                   end do

                   ! update local memory
                   ele2_nodes = ndglno( ( ele2 - 1 ) * nloc + 1 : ele2 * nloc )

                   ! figure out if we need to delete ele2
                   delete_ele2 = .false.
                   if ( ele2 == ele ) then
                      delete_ele2 = .true.
                   else
                      do iloc = 1, nloc
                         do jloc = iloc+1, nloc
                            if ( ele2_nodes( iloc ) == ele2_nodes( jloc ) ) &
                                 delete_ele2 = .true.
                         end do
                      end do
                   end if

                   if ( delete_ele2 ) then
                      ndglno( ( ele2 - 1 ) * nloc + 1 : ele2 * nloc ) = -666
                      nele_s = nele_s - 1
                   end if

                end do ! ele2, nele
             end if ! .not.on_the_wall
          end do ! edge, nedge
       end do ! ele, nele

       ! create simplified mesh

       allocate( ndglno_s( nele_s * nloc ), tmp( nele_s * nloc ) )

       count = 1
       do ele = 1, nele
          ele_nodes = ndglno( ( ele - 1 ) * nloc + 1 : ele * nloc )
          if ( all( ele_nodes < 0) ) cycle
          ndglno_s( ( count - 1 ) * nloc + 1 : count * nloc ) = ele_nodes
          count = count +1
       end do

       ! make sure elements aren't inside out
       do ele = 1, nele_s
          ele_nodes = ndglno_s( ( ele - 1 ) * nloc + 1 : ele * nloc )
          area = triangle_area( &
               f % val( 1, ele_nodes( 1 ) ), f % val( 2, ele_nodes( 1 ) ), &
               f % val( 1, ele_nodes( 2 ) ), f % val( 2, ele_nodes( 2 ) ), &
               f % val( 1, ele_nodes( 3 ) ), f % val( 2, ele_nodes( 3 ) ) )
          if ( area < 0. ) then
             ! swap 2nd and 3rd nodes
             ndglno_s( ( ele - 1 ) * nloc + 2 : ele * nloc ) = (/ ele_nodes( 3 ), ele_nodes( 2 ) /) 
          end if
       end do

       ! make sure we've recovered all the elements
       assert( nele_s == count-1 )

       ! figure out which nodes are still in the mesh and count them
       tmp = ndglno_s
       call delete_duplicates( tmp, nnodes_s )
       allocate( nodes_s( nnodes_s ) ) ; nodes_s = tmp( 1 : nnodes_s )
       ! sort the node numbers
       call ibubble( nodes_s )

       ! re-numbering

       allocate( x_s( nnodes_s ), y_s( nnodes_s ) )

       ! deal with nodes
       do i = 1, nnodes_s
          j = nodes_s( i )
          x_s( i ) = f % val( 1, j )
          y_s( i ) = f % val( 2, j )
       end do

       ! deal with ndgln
       do ele = 1, nele_s
          ele_nodes = ndglno_s( ( ele - 1 ) * nloc + 1 : ele * nloc )
          do iloc = 1, nloc
             ! figure out new node number
             do i = 1, nnodes_s
                if ( nodes_s( i ) == ele_nodes( iloc ) ) then
                   ne = i
                   exit
                end if
             end do
             ! amend ndglno
             ndglno_s( ( ele - 1 ) * nloc + iloc ) = ne
          end do
       end do

       ! deal with sndgln
       do bcs = 1, nbcs
          bc = sndglno( ( bcs - 1 ) * snloc + 1 : bcs * snloc )
          do siloc = 1, snloc
             ! figure out new node number
             do i = 1, nnodes_s
                if ( nodes_s( i ) == bc( siloc ) ) then
                   ne = i
                   exit
                end if
             end do
             ! amend sndglno
             sndglno( ( bcs - 1 ) * snloc + siloc ) = ne
          end do
       end do

       call get_option( "/geometry/quadrature/degree", quad_degree )

       quad = make_quadrature( nloc, ndim, degree = quad_degree )
       shape = make_element_shape( nloc, ndim, 1, quad )

       ! create the coarse mesh
       call allocate( mesh, nnodes_s, nele_s, shape, name="CoordinateMesh" )
       call allocate( c, ndim, mesh, name="Coordinate" )

       c % val( 1, : ) = x_s
       c % val( 2, : ) = y_s

       c % mesh % ndglno = ndglno_s

       allocate( boundary_ids( nbcs ) ) ; boundary_ids = 666
       call add_faces( c % mesh, &
            sndgln = sndglno, &
            boundary_ids = boundary_ids )

       c % dim = ndim

       deallocate( boundary_ids )

       call deallocate( mesh )
       call deallocate_element( shape )
       call deallocate( quad )

       deallocate( ndglno_s, tmp )
       deallocate( x_s, y_s )

       deallocate( ndglno, sndglno )
       deallocate( ele_nodes, ele2_nodes, n, bc )

    else

       call get_option( "/geometry/quadrature/degree", quad_degree )

       quad = make_quadrature( nloc, ndim, degree = quad_degree )
       shape = make_element_shape( nloc, ndim, 1, quad )

       call allocate( mesh, node_count( f ), element_count( f ), shape, name="CoordinateMesh" )
       call allocate( c, ndim, mesh, name="Coordinate" )
       c % val = f % val

       c % mesh % ndglno = f % mesh % ndglno

       allocate( sndglno( snloc*nbcs ) ) ; sndglno = 0
       call getsndgln( f % mesh, sndglno)

       allocate( boundary_ids( nbcs ) ) ; boundary_ids = 666
       call add_faces( c % mesh, &
            sndgln = sndglno, &
            boundary_ids = boundary_ids )

       c % dim = ndim

       deallocate( sndglno, boundary_ids )

       call deallocate( mesh )
       call deallocate_element( shape )
       call deallocate( quad )

    end if

    return
  end subroutine coarsen_mesh_2d

  !----------------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------

  subroutine ibubble( a )

    implicit none

    integer, dimension( : ), intent( inout ) :: a
    integer :: n, i, j, p, q

    n = size( a )

    do i = 1, n
       do j = n, i + 1, -1
          p = a( j - 1 ) 
          q = a( j )
          if ( p > q ) then
              a( j - 1 ) = q
              a( j ) = p
          end if
       end do
    end do

    return
  end subroutine ibubble

  subroutine delete_duplicates( a, count )

    implicit none

    integer, dimension( : ), intent( inout ) :: a
    integer, intent( out ) :: count
    integer :: n, i, j
    integer, dimension(:), allocatable :: tmp

    n = size( a )

    allocate( tmp( n ) ) ; tmp = a

    count = 1
    tmp( 1 ) = tmp( 1 )

    outer: do i = 2, n
       do j = 1, count
          ! Found a match so start looking again
          if ( tmp( j ) == a( i ) ) cycle outer
       end do
       ! No match found so add it to the output
       count = count + 1
       tmp( count ) = a( i )

    end do outer

    a = tmp

    deallocate( tmp )

    return
  end subroutine delete_duplicates

  subroutine bound_volume_fraction ( v, v_min, v_max )

    implicit none

    real, dimension( : ), intent( inout ) :: v
    real, intent( in ), optional :: v_min, v_max
    real :: vmin, vmax
    integer :: i

    vmin = 0. ; vmax = 1.
    if ( present( v_min ) ) vmin = v_min
    if ( present( v_max ) ) vmax = v_max

    do i = 1, size( v )
       v( i ) = max( vmin, min( vmax, v( i ) ) )
    end do

    return
  end subroutine bound_volume_fraction

  subroutine deallocate_femdem

    implicit none

    call deallocate( positions_r )
    call deallocate( positions_v )
    call deallocate( positions_vc )
    call deallocate( permeability_r )

    return
  end subroutine deallocate_femdem

  real function triangle_area( x1, y1, x2, y2, x3, y3 )

    implicit none

    real :: x1, y1, x2, y2, x3, y3

    triangle_area = 0.5 * ( ( x2 * y3 - y2 * x3 ) - x1 * ( y3 - y2 ) + y1 * ( x3 - x2 ) )

    return
  end function triangle_area

end module multiphase_fractures
