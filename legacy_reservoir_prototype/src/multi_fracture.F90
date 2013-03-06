
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
    use state_module

    use spud
    use global_parameters, only: option_path_len, field_name_len
    use futils, only: int2str
    use solvers
    use implicit_solids
    use FLDebug

    implicit none

#ifdef USING_FEMDEM
    interface
       subroutine y2d_allocate_femdem( string, nodes, elements, edges )
         character( len = * ), intent( in ) :: string
         integer, intent( out ) :: nodes, elements, edges
       end subroutine y2d_allocate_femdem
    end interface

    interface
       subroutine y2d_populate_femdem( ele1, ele2 ,ele3, &
            face1, face2, xs, ys, p1, p2, p3, p4 )
         integer, dimension( * ), intent( out ) :: ele1, ele2, ele3
         integer, dimension( * ), intent( out ) :: face1, face2
         real, dimension( * ), intent( out ) :: xs, ys, p1, p2, p3, p4
       end subroutine y2d_populate_femdem
    end interface

    interface 
       subroutine y2dfemdem( xs, ys, p1, p2, p3, p4, p )
         real, dimension( * ), intent( in ) :: p
         real, dimension( * ), intent( out ) :: xs, ys, p1, p2, p3, p4
       end subroutine y2dfemdem
    end interface
#endif

    type( vector_field ), save :: fracture_positions
    type( tensor_field ), save :: permeability

    private
    public :: fractures

  contains

    subroutine fractures( states, totele, ndim, perm )

      integer, intent( in ) :: totele, ndim
      type( state_type ), dimension( : ), intent( in ) :: states
      real, dimension( totele, ndim, ndim ), intent( inout ) :: perm

      ewrite(3,*) 'inside fractures'

      ! define: 1. fracture_positions and 2. permeability
      call initialise_fractures( states )

      ! interpolate conservatively permeability from
      ! the fracture mesh to the prototype mesh
      ! that is interpolate "permeability" to "perm"
      call interpolate_fractures( states, totele, ndim, perm )

      ewrite(3,*) 'leaving fractures'

      return
    end subroutine fractures

    !----------------------------------------------------------------------------------------------------------

    subroutine initialise_fractures( states )

      type( state_type ), dimension( : ), intent( in ) :: states

      character( len = FIELD_NAME_LEN ) :: fracture_mesh_name
      integer :: i, loc, sloc
      integer :: ndim, nodes, elements, edges
      integer, dimension( : ), allocatable :: ele1, ele2, ele3
      integer, dimension( : ), allocatable :: face1, face2
      real, dimension( : ), allocatable :: xs, ys, p1, p2, p3, p4
      type( quadrature_type ) :: quad
      type( element_type ) :: shape
      integer, dimension( : ), allocatable :: sndglno, boundary_ids
      integer :: quad_degree
      type( mesh_type ) :: mesh

      ewrite(3,*) 'inside initialise_fractures'

      call get_option( '/porous_media/Permeability_from_femdem/name', fracture_mesh_name )
      call get_option( "/geometry/quadrature/degree", quad_degree )
      call get_option( "/geometry/dimension", ndim )

      if ( ndim ==2 ) then
         loc = 3
         sloc= 2
      else if ( ndim == 3 ) then
         loc = 4
         sloc= 3
      else
         FLAbort("Fracture modelling is supported for 2 and 3D only.")
      end if

      call y2d_allocate_femdem( trim( fracture_mesh_name )//char(0), &
           nodes, elements, edges )

      allocate( ele1( elements ) ) ; allocate( ele2( elements ) )
      allocate( ele3( elements ) )
      allocate( face1( edges ) ) ; allocate( face2( edges ) )

      allocate( xs( nodes ) ) ; allocate( ys( nodes ) )

      allocate( p1( elements ) ) ; allocate( p2( elements ) )
      allocate( p3( elements ) ) ; allocate( p4( elements ) )

      call y2d_populate_femdem( ele1, ele2, ele3, &
           face1, face2, xs, ys, p1, p2, p3, p4 )
  
      quad = make_quadrature( loc, ndim, degree = quad_degree )
      shape = make_element_shape( loc, ndim, 1, quad )

      call allocate( mesh, nodes, elements, shape, name="CoordinateMesh" )
      call allocate( fracture_positions, ndim, mesh, name="Coordinate" )

      fracture_positions%val( 1, : ) = xs
      fracture_positions%val( 2, : ) = ys

      do i = 1, elements
         fracture_positions%mesh%ndglno( (i-1)*loc+1 : i*loc ) = &
              (/ ele1(i)+1, ele2(i)+1, ele3(i)+1 /)
      end do

      allocate( sndglno( edges * sloc ) ) ; sndglno = 0
      allocate( boundary_ids( edges) ) ; boundary_ids = 666

      do i = 1, edges
         sndglno( (i-1)*sloc+1 : i*sloc ) = &
              (/ face1(i)+1, face2(i)+1 /)
      end do

      call add_faces( fracture_positions%mesh, &
           sndgln = sndglno, &
           boundary_ids = boundary_ids )

      fracture_positions%dim = ndim

      deallocate( sndglno, boundary_ids )
      call deallocate( mesh )
      call deallocate_element( shape )
      call deallocate( quad )

      ! store permeability
      call allocate( permeability, fracture_positions%mesh, name = "Permeability" )
      call zero( permeability )

      call set_all( permeability, 1, 1, p1 )
      call set_all( permeability, 1, 2, p2 )
      call set_all( permeability, 2, 1, p3 )
      call set_all( permeability, 2, 2, p4 )

      deallocate( ele1, ele2, ele3, face1, face2, xs, ys, p1, p2, p3, p4 )

      ewrite(3,*) 'leaving initialise_fractures'

      return
    end subroutine initialise_fractures


    subroutine interpolate_fractures( states, totele, ndim, perm )

      integer, intent( in ) :: totele, ndim
      type( state_type ), dimension( : ), intent( in ) :: states
      real, dimension( totele, ndim, ndim ), intent( inout ) :: perm

      type( state_type ) :: alg_ext, alg_fl
      type( mesh_type ), pointer :: fl_mesh, p0_fl_mesh
      type( vector_field ) :: fl_positions
      type( scalar_field ) :: dummy

      type( scalar_field ), pointer :: pressure

      type( scalar_field ) :: field_fl_p11, field_fl_p12, field_fl_p21, field_fl_p22
      type( scalar_field ) :: field_ext_p11, field_ext_p12, field_ext_p21, field_ext_p22

      character( len = OPTION_PATH_LEN ) :: &
           path = "/tmp/galerkin_projection/continuous"

      ewrite(3,*) 'inside interpolate_fractures'

      pressure => extract_scalar_field( states(1), "Pressure" )

      call insert( alg_ext, fracture_positions%mesh, "Mesh" )
      call insert( alg_ext, fracture_positions, "Coordinate" )

      fl_mesh => extract_mesh( states(1), "CoordinateMesh" )
      fl_positions = extract_vector_field( states(1), "Coordinate" )

      call insert( alg_fl, fl_mesh, "Mesh" )
      call insert( alg_fl, fl_positions, "Coordinate" )

      call set_solver_options(path, &
           ksptype = "cg", &
           pctype = "mg", &
           rtol = 1.e-10, &
           atol = 0., &
           max_its = 10000)

      path = "/tmp"

      p0_fl_mesh => extract_mesh( states(1), "P0Mesh" )

      ewrite(3,*) '...generating multiphase state'

      ! this is the permeability on the fluidity mesh
      call allocate( field_fl_p11, p0_fl_mesh, "Permeability11" )
      call zero( field_fl_p11 )
      call insert( alg_fl, field_fl_p11, "Permeability11" )
      field_fl_p11%option_path = path

      call allocate( field_fl_p12, p0_fl_mesh, "Permeability12" )
      call zero( field_fl_p12 )
      call insert( alg_fl, field_fl_p12, "Permeability12" )
      field_fl_p12%option_path = path

      call allocate( field_fl_p21, p0_fl_mesh, "Permeability21" )
      call zero( field_fl_p21 )
      call insert( alg_fl, field_fl_p21, "Permeability21" )
      field_fl_p21%option_path = path

      call allocate( field_fl_p22, p0_fl_mesh, "Permeability22" )
      call zero( field_fl_p22 )
      call insert( alg_fl, field_fl_p22, "Permeability22" )
      field_fl_p22%option_path = path

      ewrite(3,*) '...generating femdem state'

      ! this is the permeability on the solid mesh
      field_ext_p11 = extract_scalar_field( permeability, 1, 1 )
      call insert( alg_ext, field_ext_p11, "Permeability11" )

      field_ext_p12 = extract_scalar_field( permeability, 1, 2 )
      call insert( alg_ext, field_ext_p12, "Permeability12" )

      field_ext_p21 = extract_scalar_field( permeability, 2, 1 )
      call insert( alg_ext, field_ext_p21, "Permeability21" )

      field_ext_p22 = extract_scalar_field( permeability, 2, 2 )
      call insert( alg_ext, field_ext_p22, "Permeability22" )

      ! dummy field needed for the interpolation
      call allocate( dummy, pressure%mesh, "Dummy" )
      call zero( dummy )

      ewrite(3,*) '...interpolating'

      ! interpolate
      call interpolation_galerkin_femdem( alg_ext, alg_fl, field = dummy )

      ! copy memory to prototype
      perm = 0.
      perm( :, 1, 1 ) = field_fl_p11 % val
      perm( :, 1, 2 ) = field_fl_p12 % val
      perm( :, 2, 1 ) = field_fl_p21 % val
      perm( :, 2, 2 ) = field_fl_p22 % val

      ! now deallocate
      call deallocate( dummy )

      call deallocate( field_ext_p22 )
      call deallocate( field_ext_p21 )
      call deallocate( field_ext_p12 )
      call deallocate( field_ext_p11 )

      call deallocate( field_fl_p22 )
      call deallocate( field_fl_p21 )
      call deallocate( field_fl_p12 )
      call deallocate( field_fl_p11 )

      call deallocate( alg_fl )
      call deallocate( alg_ext )

      ewrite(3,*) 'leaving interpolate_fractures'

      return
    end subroutine interpolate_fractures


    subroutine deallocate_fractures

      call deallocate( fracture_positions )
      call deallocate( permeability )

      return

    end subroutine deallocate_fractures


  end module multiphase_fractures
