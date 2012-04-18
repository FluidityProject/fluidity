  !
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
  
  module pascal_tetrahedra

    use fldebug

    integer, parameter :: nlayer = 2, totele = 8, tet_totele = 4, tri_totele = 3, &
         quad_cv_nloc = 8, no_of_nodes_in_faces = 7, no_faces = 4, nside = 3
    integer, parameter :: no_elements = 8 ! this is hard-coded for now
    real, parameter :: h_scale = 2.0396489026555056 ! Scaling factor to give a unity volume of the local tetrahedron

  contains


    integer function tet_nodes( nlayer_prev, nlayer_now )
      implicit none
      integer :: nlayer_prev, nlayer_now
      integer :: ilayer

      tet_nodes = 0
      do ilayer = nlayer_prev, nlayer_now
         tet_nodes = tet_nodes + ( ilayer + 1 ) * ( ilayer + 2 ) / 2
      end do

      return
    end function tet_nodes

    integer function compute_permutation( n )
      implicit none
      integer :: n
      integer :: i

      compute_permutation = 1
      do i = n, 1, 1
         compute_permutation = compute_permutation * i 
      end do

      return
    end function compute_permutation


    subroutine numbering_tet_elements( x_nonods, triangle_totele, nodeplustetnodes, &
         x_ndgln_big, &
         x, y, z, lx, ly, lz, fem_nod, &
         x_ndgln, x_ndgln2 ) 
      implicit none
      integer, intent( in ) :: x_nonods, triangle_totele
      integer, intent( inout ) :: nodeplustetnodes
      integer, dimension( x_nonods ), intent( inout ) :: x_ndgln_big
      real, dimension( x_nonods ), intent( inout ) :: x, y, z, lx, ly, lz
      integer, dimension( x_nonods ), intent( inout ) :: fem_nod
      integer, dimension( x_nonods ), intent( inout ) :: x_ndgln, x_ndgln2
      ! Local variables
      integer :: ifem ! counter for fem_nod
      integer :: nnods, ele, ele2, ilayer, ilayer0, iloc, iloc2, iside, inod, nod, plusnodes
      logical :: repetitive_nodes

      lx = 0; ly = 0; lz = 0

      ! Setting-up unity volume tetrahedron
      ilayer = nlayer ! Basis of the large tetrahedron

      lx( tet_nodes( 0, ilayer - 1 ) + 1 + ilayer * 0 ) = 0. ! Point node 5
      ly( tet_nodes( 0, ilayer - 1 ) + 1 + ilayer * 0 ) = 0.
      lz( tet_nodes( 0, ilayer - 1 ) + 1 + ilayer * 0 ) = 0.

      lx( tet_nodes( 0, ilayer - 1 ) + 1 + ilayer * 1 ) = 1. * h_scale ! Point node 7
      ly( tet_nodes( 0, ilayer - 1 ) + 1 + ilayer * 1 ) = 0.
      lz( tet_nodes( 0, ilayer - 1 ) + 1 + ilayer * 1 ) = 0.

      lx( tet_nodes( 0, ilayer - 1 ) + 1 + ilayer * 2 ) = 0.5 * h_scale  ! Point node 9
      ly( tet_nodes( 0, ilayer - 1 ) + 1 + ilayer * 2 ) = sqrt( 3. / 4. ) * h_scale
      lz( tet_nodes( 0, ilayer - 1 ) + 1 + ilayer * 2 ) = 0.

      ilayer = 0 ! Top of the large tetrahedron (Point node 1 )
      ilayer0 = nlayer
      lx( tet_nodes( 0, ilayer - 1 ) + 1 + ilayer * 0 ) = 0.5 * &
           ( lx( tet_nodes( 0, ilayer0 - 1 ) + 1 ) + &
           lx( tet_nodes( 0, ilayer0 - 1 ) + 1 + ilayer0 ) ) !* h_scale
      ly( tet_nodes( 0, ilayer - 1 ) + 1 + ilayer * 0 ) = 0.5 * &
           ( ly( tet_nodes( 0, ilayer0 - 1 ) + 1 ) + & 
           ly( tet_nodes( 0, ilayer0 - 1 ) + 1 + ilayer0 * 2 ) ) !* h_scale 
      lz( tet_nodes( 0, ilayer - 1 ) + 1 + ilayer * 0 ) = &
           sqrt( 2. / 3. ) * h_scale

      !ewrite(3,*)  tet_nodes( 0, ilayer - 1 ) + 1 + ilayer * 0, tet_nodes( 0, ilayer0 - 1 ) + 
      !1,tet_nodes( 0, ilayer0 - 1 ) + 1 + ilayer0
      !ewrite(3,*)  lx(1), lx( 5), lx(7), lx(9)
      !stop 986

      fem_nod = 0! FEM_NOD contains the local triangle/tet FEM node
      z = 0.

!!!
!!! Numbering point nodes of the Big Tetrahedron.
!!! 
      ifem = 0
      Loop_Layers0: do ilayer = 0, nlayer, nlayer
         ! Number of nodes in the top (node point 1, ilayer = 0) and the basis (node points
         ! 2, 3, 4, ilayer = nlayer )
         nnods = max( 1, 3 + 3 * ( ilayer - 1 ) ) 
         do iloc = 1, nnods,  max( 1, ilayer )
            ifem = ifem + 1
            x( tet_nodes( 0, ilayer - 1 ) + iloc ) = lx( tet_nodes( 0, ilayer - 1 ) + iloc )
            y( tet_nodes( 0, ilayer - 1 ) + iloc ) = ly( tet_nodes( 0, ilayer - 1 ) + iloc )
            z( tet_nodes( 0, ilayer - 1 ) + iloc ) = lz( tet_nodes( 0, ilayer - 1 ) + iloc )
            fem_nod( tet_nodes( 0, ilayer - 1 ) + iloc ) = ifem
            iloc2 = tet_nodes( 0, ilayer - 1 ) + iloc
            ewrite(3,*)  ilayer, iloc2, x( iloc2 ), y( iloc2 ), z( iloc2 )
         end do
      end do Loop_Layers0

      ewrite(3,*)  '  '
!!!
!!! Numbering the remaining node points in the *BASIS* of the large tetrahedron,
!!! i.e., ilayer =nlayer
!!!
      ilayer = nlayer        
      ewrite(3,*)'X-Coord:'
      call Mapping_QTet_TriangleSide( .true., ilayer, x_nonods, &
           x, fem_nod, ifem )
      ewrite(3,*)'Y-Coord:'
      call Mapping_QTet_TriangleSide( .false., ilayer, x_nonods, &
           y, fem_nod, ifem )
      ewrite(3,*)'Z-Coord:'
      call Mapping_QTet_TriangleSide( .false., ilayer, x_nonods, &
           z, fem_nod, ifem )
!!!
!!! Now for the node points in the centre of the Pascal tetrahedron at 
!!! ilayer = nlayer (i.e., basis of the tetrahedron).
!!!
      ewrite(3,*)'Centre:'
      call Mapping_QTet_TriangleCentre( ilayer, x_nonods, x, y, z )

!!!
!!! Numbering node points at the intermediate layers
!!!
      ilayer = 0
      Loop_Layers10: do ilayer = ( nlayer - 1 ), 1, -1

         ! Node at the vertices of each layer
         ewrite(3,*)'X/Y/Z-Coord Vertices, ilayer:', ilayer 
         call Mapping_QTet_TriangleVertices( .true., ilayer, x_nonods, &
              x, y, z, fem_nod, ifem ) 

         ! Node at the side of each layer
         ewrite(3,*)'X-Coord Sides, ilayer:', ilayer 
         call Mapping_QTet_TriangleSide( .true., ilayer, x_nonods, &
              x, fem_nod, ifem ) 
         ewrite(3,*)'Y-Coord Sides, ilayer:', ilayer 
         call Mapping_QTet_TriangleSide( .false., ilayer, x_nonods, &
              y, fem_nod, ifem )
         ewrite(3,*)'Z-Coord Sides, ilayer:', ilayer 
         call Mapping_QTet_TriangleSide( .false., ilayer, x_nonods, &
              z, fem_nod, ifem )

         ! Now the node points in the central region of the intermediate layers
         ewrite(3,*)'Nodes in the centre of the triangle at layer', ilayer
         call Mapping_QTet_TriangleCentre( ilayer, x_nonods, x, y, z )

      end do Loop_Layers10

      ! Sanity check on fem_nod (just checking the inventory of fem_nod -  
      ! if there is any node index been repeated)
      if( QTet_CheckRepetitiveNodes( x_nonods, fem_nod ) ) stop 901

!!! From here, everything is hard-coded, we need an efficient algorithm that generates
!!! tetrahedra between layers. Here just the 8 tetrahedra from ilayers = [ 0, 2 ]

      ! Numbering x_ndgln
      call Numbering_QTets( x_nonods, x_ndgln )

!!! Now creating triangles from the small tetrahedra stored in x_ndgln and adding the 
!!! control volume nodes. x_ndgln2 now holds node points for all tetrahedra and the 
!!! associated coordinates.
      call QTet_CV( x_nonods, triangle_totele, nodeplustetnodes, &
           x_ndgln, x_ndgln2, x, y, z )

!!! Now setting up CV hexahedra from the list of node points stored in x_ndgln2. The
!!! new numbering is stored in x_ndgln_big
      do ele = 1, totele
         call Mapping_LinTet2Quad( ele, x_nonods, &
              x_ndgln2, x_ndgln_big )
      end do

      do ele2 = 1, 32
         do iloc = 1, 8
            if( x_ndgln_big( ( ele2 - 1 ) * quad_cv_nloc + iloc ) <= 10 ) then
               ewrite(3,*)'nod, fem_nod:', x_ndgln_big( ( ele2 - 1 ) * quad_cv_nloc + iloc ), &
                    fem_nod( x_ndgln_big( ( ele2 - 1 ) * quad_cv_nloc + iloc ) )
            end if
         end do
      end do
      !stop 987

      !x_nonods = maxval(x_ndgln_big)

      return
    end subroutine numbering_tet_elements


    subroutine Mapping_QTet_TriangleVertices( mapfemnode, ilayer, x_nonods, &
         x, y, z, fem_nod, ifem )
      ! This subrt maps the node in the vertices of triangle at *ilayer*. If MapFEMnode
      ! is true then femnode is updated
      implicit none
      logical, intent( in ) :: mapfemnode
      integer, intent( in ) :: ilayer, x_nonods
      real, dimension( x_nonods ), intent( inout ) :: x, y, z
      integer, dimension( x_nonods ), intent( inout ) :: fem_nod
      integer, intent( inout ) :: ifem
      ! Local variables
      integer :: ilayer0, node, nodeplus1, nodeplus2, iside
      real :: dx, dy, dz

      ilayer0 = ilayer + 1
      ! Now node points at the vertices of the triagles (from the bottom to the top)
      ! X-Coordinate
      node = tet_nodes( 0, ilayer - 1 ) + 1
      nodeplus1 = tet_nodes( 0, ilayer - 1 ) + 1 + ilayer
      nodeplus2 = tet_nodes( 0, ilayer - 1 ) + 1 + 2 * ilayer

      Conditional_Layer: if( ilayer == 1 ) then   
         x( node ) = 1. / 4. * &
              ( x( tet_nodes( 0, ilayer ) + 1 + ilayer0 ) - x( tet_nodes( 0, ilayer ) + 1 ) ) + &
              x( tet_nodes( 0, ilayer ) + 1 )
         x( nodeplus1 ) = 3. / 4. * &
              ( x( tet_nodes( 0, ilayer ) + 1 + ilayer0 ) - x( tet_nodes( 0, ilayer ) + 1 ) ) + &
              x( tet_nodes( 0, ilayer ) + 1 )
         x( nodeplus2 ) = x( tet_nodes( 0, ilayer ) + 1 + 2 * ilayer0 )

         y( node ) = 1. / 4. * &
              ( y( tet_nodes( 0, ilayer ) + 1 + 2 * ilayer0 ) -  &
              y( tet_nodes( 0, ilayer ) + 1 ) ) + y( tet_nodes( 0, ilayer ) + 1 )
         y( nodeplus1 ) = 1. / 4. * &
              ( y( tet_nodes( 0, ilayer ) + 1 + 2 * ilayer0 ) -  &
              y( tet_nodes( 0, ilayer ) + 1 ) ) + y( tet_nodes( 0, ilayer ) + 1 )
         y( nodeplus2 ) = 3. / 4. * &
              ( y( tet_nodes( 0, ilayer ) + 1 + 2 * ilayer0 ) -  &
              y( tet_nodes( 0, ilayer ) + 1 ) ) + y( tet_nodes( 0, ilayer ) + 1 )

      else
         x( node ) = 1. / real( ilayer0 ) * &
              ( x( tet_nodes( 0, ilayer ) + 1 + ilayer0 ) - x( tet_nodes( 0, ilayer ) + 1 ) ) + &
              x( tet_nodes( 0, ilayer ) + 1 )
         x( nodeplus1 ) = real( ilayer0 - 1 ) / real( ilayer0 ) * &
              ( x( tet_nodes( 0, ilayer ) + 1 + ilayer0 ) - x( tet_nodes( 0, ilayer ) + 1 ) ) + &
              x( tet_nodes( 0, ilayer ) + 1 )
         x( nodeplus2 ) = x( tet_nodes( 0, ilayer ) + 1 + 2 * ilayer0 )


         y( node ) = 1. / real( ilayer0 ) * &
              ( y( tet_nodes( 0, ilayer ) + 1 + 2 * ilayer0 ) -  &
              y( tet_nodes( 0, ilayer ) + 1 ) ) + y( tet_nodes( 0, ilayer ) + 1 )
         y( nodeplus1 ) = 1. / real( ilayer0 ) * &
              ( y( tet_nodes( 0, ilayer ) + 1 + 2 * ilayer0 ) -  &
              y( tet_nodes( 0, ilayer ) + 1 ) ) + y( tet_nodes( 0, ilayer ) + 1 )
         y( nodeplus2 ) = real( ilayer0 - 1 ) / real( ilayer0 ) * &
              ( y( tet_nodes( 0, ilayer ) + 1 + 2 * ilayer0 ) -  &
              y( tet_nodes( 0, ilayer ) + 1 ) ) + y( tet_nodes( 0, ilayer ) + 1 )

      end if Conditional_Layer
      z( node ) = QTet_Height( ilayer, x_nonods, z )
      z( nodeplus1 ) = QTet_Height( ilayer, x_nonods, z )
      z( nodeplus2 ) = QTet_Height( ilayer, x_nonods, z )

      if( mapfemnode ) then
         do iside = 1, nside
            ifem = ifem + 1
            fem_nod( tet_nodes( 0, ilayer - 1 ) + 1 + ( iside - 1 ) * ilayer ) = ifem
         end do
      end if

      ewrite(3,*)  node, x( node ), y( node ), z( node )
      ewrite(3,*)  nodeplus1, x( nodeplus1 ), y( nodeplus1 ), z( nodeplus1 )
      ewrite(3,*)  nodeplus2, x( nodeplus2 ), y( nodeplus2 ), z( nodeplus2 )


      return
    end subroutine Mapping_QTet_TriangleVertices

    subroutine Mapping_QTet_TriangleSide( mapfemnode, ilayer, x_nonods, &
         x, fem_nod, ifem )
      ! This subrt maps the node in the main side of the triangle at *ilayer*. If MapFEMnode
      ! is true then femnode is updated
      implicit none
      logical, intent( in ) :: mapfemnode
      integer, intent( in ) :: ilayer, x_nonods
      real, dimension( x_nonods ), intent( inout ) :: x
      integer, dimension( x_nonods ), intent( inout ) :: fem_nod
      integer, intent( inout ) :: ifem
      ! Local variables
      integer :: iside, iloc, node, nodeplus, nodeminus
      real :: dx

      Loop_Side0: do iside = 1, nside
         Loop_ILOC0: do iloc = 2, ilayer
            node = tet_nodes( 0, ilayer - 1 ) + ( iside - 1 ) * ilayer + iloc
            if( iside /= 3 ) then
               nodeplus = tet_nodes( 0, ilayer - 1 ) + 1 + iside * ilayer 
               nodeminus = tet_nodes( 0, ilayer - 1 ) + 1 + ( iside - 1 ) * ilayer 
               dx = real( iloc - 1 ) / real( ilayer ) * ( x( nodeplus ) - x( nodeminus ) )
               x( node ) = x( nodeminus ) + dx
            else
               nodeplus = tet_nodes( 0, ilayer - 1 ) + 1
               nodeminus = tet_nodes( 0, ilayer - 1 ) + 1 + ( iside - 1 ) * ilayer
               dx = real( iloc - 1 ) / real( ilayer ) * ( x( nodeplus ) - x( nodeminus ) )
               x( node ) = x( nodeminus ) + dx
            end if
            if( mapfemnode ) then
               ifem = ifem + 1
               fem_nod( node ) = ifem
            end if
            !ewrite(3, '(a10,i3,i3,g12.5)')'SideNodes:', ilayer, node, x( node )
         end do Loop_ILOC0
      end do Loop_Side0
      !ewrite(3,*)  '++'

      return
    end subroutine Mapping_QTet_TriangleSide

    subroutine Mapping_QTet_TriangleCentre( ilayer, x_nonods, x, y, z )
      ! This subrt maps the node in the centre of the triangle at *ilayer*.  There are
      ! tet_nodes( ilayer, ilayer ) - ( 3 + ( ilayer - 1 ) * nside ) nodes in the
      ! centre region
      implicit none
      integer, intent( in ) :: ilayer, x_nonods
      real, dimension( x_nonods ), intent( inout ) :: x, y, z
      ! Local variables
      integer :: node, iloc, ilayer0, nnods, nod, inod

      nnods = tet_nodes( ilayer, ilayer )
      nod = nnods - ( 3 + ( ilayer - 1 ) * nside ) ! Computing the number of centre node points

      Conditional_Centre_Nodes: if( nod > 0 ) then
         node = tet_nodes( 0, ilayer - 1 ) + ( 3 + ( ilayer - 1 ) * nside )
         inod = 0
         Loop_Centre_Nodes0: do ilayer0 = 1, ilayer - 2
            Loop_Centre_Nodes1: do iloc = 1, ilayer0
               inod = inod + 1
               ! Node Point - X coordinates
               x( node + inod ) = real( iloc ) / real ( ilayer0 + 1 ) * &
                    x( tet_nodes( 0, ilayer - 1 ) + 1 + ilayer ) 
               ! Node Point - Y coordinates
               y( node + inod ) = &
                    y( tet_nodes( 0, ilayer - 1 ) + 1 + 2 * ilayer + 1 + ilayer0 )
               ! Node Point - Z coordinates 
               z( node + inod ) = z( tet_nodes( 0, ilayer - 1 ) + 1 )
               ewrite(3,'(i3,1x,i3,3g12.5)') ilayer, node + inod , x( node + inod ), y( node + inod ), z( node + inod )
            end do Loop_Centre_Nodes1
         end do Loop_Centre_Nodes0
      end if Conditional_Centre_Nodes

      return
    end subroutine Mapping_QTet_TriangleCentre

    subroutine Numbering_QTets( x_nonods, x_ndgln )
      ! This (hard-coded) subrt allocates x_ndgln for a given set of tetrahedra (8)
      ! from the large tetrahedron. 
      implicit none
      integer, intent( in ) :: x_nonods
      integer, dimension( x_nonods ), intent( inout ) :: x_ndgln
      ! Local variables
      integer :: ele, ilayer

      ele = 1 ! Element 1 ( 1 2 3 4 )
      ilayer = 0
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 1 ) = tet_nodes( 0, ilayer ) ! 1
      ilayer = 1
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 2 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 0 * ilayer ! 2
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 3 ) = & 
           tet_nodes( 0, ilayer - 1 ) + 1 + 1 * ilayer ! 3
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 4 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 2 * ilayer ! 4
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 5 ) = tet_nodes( 0, nlayer ) + ele

      ele = 2 ! Element 2 ( 2 5 6 10 )
      ilayer = 1
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 1 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 0 * ilayer ! 2
      ilayer = 2
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 2 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 0 * ilayer ! 5
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 3 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 1 + 0 * ilayer ! 6
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 4 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 2 * ilayer + 1 ! 10
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 5 ) = tet_nodes( 0, nlayer ) + ele

      ele = 3 ! Element 3 ( 2 3 6 10 )
      ilayer = 1
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 1 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 0 * ilayer ! 2
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 2 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 1 * ilayer ! 3
      ilayer = 2
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 3 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 1 + 0 * ilayer ! 6
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 4 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 2 * ilayer + 1 ! 10
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 5 ) = tet_nodes( 0, nlayer ) + ele

      ele = 4 ! Element 4 ( 3 6 7 10 )
      ilayer = 1
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 1 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 1 * ilayer  ! 3
      ilayer = 2
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 2 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 1 + 0 * ilayer ! 6
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 3 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 1 * ilayer ! 7
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 4 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 2 * ilayer + 1 ! 10
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 5 ) = tet_nodes( 0, nlayer ) + ele

      ele = 5 ! Element 5 ( 3 7 8 10 )
      ilayer = 1
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 1 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 1 * ilayer ! 3
      ilayer = 2
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 2 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 1 * ilayer ! 7
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 3 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 1 * ilayer + 1 ! 8
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 4 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 2 * ilayer + 1 ! 10
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 5 ) = tet_nodes( 0, nlayer ) + ele

      ele = 6 ! Element 6 ( 3 4 8 10 )
      ilayer = 1
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 1 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 1 * ilayer ! 3
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 2 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 2 * ilayer ! 4
      ilayer = 2
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 3 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 1 * ilayer + 1 ! 8
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 4 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 2 * ilayer + 1! 10
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 5 ) = tet_nodes( 0, nlayer ) + ele

      ele = 7 ! Element 7 ( 4 8 9 10 )
      ilayer = 1
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 1 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 2 * ilayer ! 4
      ilayer = 2
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 2 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 1 * ilayer + 1 ! 8
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 3 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 2 * ilayer ! 9
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 4 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 2 * ilayer + 1 ! 10
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 5 ) = tet_nodes( 0, nlayer ) + ele

      ele = 8 ! Element 8 ( 2 3 4 10 )
      ilayer = 1
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 1 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 0 * ilayer ! 2 
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 2 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 1 * ilayer ! 3
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 3 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 2 * ilayer ! 4
      ilayer = 2
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 4 ) = &
           tet_nodes( 0, ilayer - 1 ) + 1 + 2 * ilayer + 1 ! 10
      x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 5 ) = tet_nodes( 0, nlayer ) + ele

      return
    end subroutine Numbering_QTets

    subroutine QTet_CV( x_nonods, triangle_totele, nodeplustetnodes, &
         x_ndgln, x_ndgln2, x, y, z )
      ! This subrt maps and allocates the CV nodes
      implicit none
      integer, intent( in ) :: x_nonods, triangle_totele
      integer, intent( inout ) :: nodeplustetnodes
      integer, dimension( x_nonods ), intent( inout ) :: x_ndgln, x_ndgln2
      real, dimension( x_nonods ), intent( inout ) :: x, y, z
      ! Local variables
      integer, dimension( 12 ), save :: inc_tri = &
           (/ 1, 2, 3, 1, 2, 4, 1, 3, 4, 2, 3, 4 /)
      integer :: ele, ele2, iloc, iloc2, node, node1, node2, node3, node4

      ewrite(3,*)'in QTet_CV'

      ! Computing node point in the centre of the tetrahedra
      Loop_Elements_Tets0: do ele = 1, totele
         node = x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 5 )
         x( node ) = 0.
         y( node ) = 0.
         z( node ) = 0.
         Loop_Iloc0: do iloc = 1, tet_totele 
            x( node ) = x( node ) + 0.25 * x( x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + iloc ) )
            y( node ) = y( node ) + 0.25 * y( x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + iloc ) )
            z( node ) = z( node ) + 0.25 * z( x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + iloc ) )
         end do Loop_Iloc0
      end do Loop_Elements_Tets0

      ! Allocating 4 triangles within each of the TOTELE tetrahedra
      nodeplustetnodes = tet_nodes( 0, nlayer )
      Loop_Elements_Tets1: do ele = 1, totele 
         Loop_Elements_Triangle: do ele2 =  1, tet_totele
            Loop_Iloc_Triangle: do iloc2 = 1, tri_totele
               node = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + &
                    ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )  
               x_ndgln2( node + iloc2 ) = &
                    x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + inc_tri( ( ele2 - 1 ) * tri_totele + iloc2 ) )
               ewrite(3,*)  ele, ele2, iloc2 + node, x_ndgln2( node + iloc2 )
            end do Loop_Iloc_Triangle
         end do Loop_Elements_Triangle
         nodeplustetnodes = nodeplustetnodes + 1
         x_ndgln2( ele * tet_totele * no_of_nodes_in_faces + ele ) = &
              x_ndgln( ( ele - 1 ) * ( tet_totele + 1 ) + 5 )
         ewrite(3,*)  '->:', ele, ele * tet_totele * no_of_nodes_in_faces + ele, &
              x_ndgln2( ele * tet_totele * no_of_nodes_in_faces +  ele )
      end do Loop_Elements_Tets1


      ! Now mapping the CV node points
      ! nodeplustetnodes = tet_nodes( 0, nlayer )
      Loop_Elements_Tets2: do ele = 1, totele
         ewrite(3,*)'Ele / Totele:', ele, totele
         call Mapping_QTet_CV( ele, x_nonods, nodeplustetnodes, &
              x_ndgln2, x, y, z )
      end do Loop_Elements_Tets2

      ! Renumbering and eliminating repetitive numbering
      call Renumbering_QTet_CV( x_nonods, x_ndgln2, &
           x, y, z )
      if( .true. ) call checking_QTet_CV( x_nonods, x_ndgln2 )

      return
    end subroutine QTet_CV


    subroutine Mapping_QTet_CV( ele, x_nonods, nodeplustetnodes, &
         x_ndgln2, x, y, z )
      ! This subroutine maps the CV's of each element ELE and allocates 
      ! node-points ( x_ndgln ) and coordinates (X/Y/Z)
      implicit none
      integer, intent( in ) :: ele, x_nonods
      integer, intent( inout ) :: nodeplustetnodes
      integer, dimension( x_nonods ), intent( inout ) :: x_ndgln2
      real, dimension( x_nonods ), intent( inout ) :: x, y, z
      ! Local variables
      integer :: ele2, node, iloc, iloc2, inode, nodeplus, nodeminus, nodeplus1

      ewrite(3,*)'in Mapping_QTet_CV'

      Loop_Elements_Triangle: do ele2 = 1, tet_totele
         inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
              ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 ) 
         Loop_Sides: do iloc = 1, tri_totele
            nodeplustetnodes = nodeplustetnodes + 1
            Conditional_Sides: if( iloc /= tri_totele ) then
               nodeplus = x_ndgln2( inode + iloc + 1 )
               nodeminus = x_ndgln2( inode + iloc )
               x_ndgln2( inode + 3 + iloc ) = nodeplustetnodes
               node = x_ndgln2( inode + 3 + iloc )
               x( node ) = 0.5 * ( x( nodeplus ) + x( nodeminus ) )
               y( node ) = 0.5 * ( y( nodeplus ) + y( nodeminus ) )
               z( node ) = 0.5 * ( z( nodeplus ) + z( nodeminus ) )
               ewrite(3,*) ele2, iloc, nodeminus, nodeplus, node
            else
               nodeplus = x_ndgln2( inode + iloc - 2 )
               nodeminus = x_ndgln2( inode + iloc )
               x_ndgln2( inode + 3 + iloc ) = nodeplustetnodes
               node = x_ndgln2( inode + 3 + iloc )
               x( node ) = 0.5 * ( x( nodeplus ) + x( nodeminus ) )
               y( node ) = 0.5 * ( y( nodeplus ) + y( nodeminus ) )
               z( node ) = 0.5 * ( z( nodeplus ) + z( nodeminus ) )
               ewrite(3,*) ele2, iloc, nodeminus, nodeplus, node
            endif Conditional_Sides
            !ewrite(3,*) node, x( node ), y( node ), z( node )
            !ewrite(3,*) ele2, iloc, node
         end do Loop_Sides
         iloc = 4 ! Centre node point
         nodeplustetnodes = nodeplustetnodes + 1
         x_ndgln2( inode + 3 + iloc ) = nodeplustetnodes
         node = x_ndgln2( inode + 3 + iloc )
         nodeminus = x_ndgln2( inode + 1 )
         nodeplus =  x_ndgln2( inode + 2 )
         nodeplus1 = x_ndgln2( inode + 3 )
         x( node ) = 1. / 3. * ( x( nodeplus1 ) + x( nodeplus ) + x( nodeminus ) )
         y( node ) = 1. / 3. * ( y( nodeplus1 ) + y( nodeplus ) + y( nodeminus ) )
         z( node ) = 1. / 3. * ( z( nodeplus1 ) + z( nodeplus ) + z( nodeminus ) )
         ! ewrite(3,*) node, x( node ), y( node ), z( node )
         !ewrite(3,*) ele2, iloc, node
         ewrite(3,*) ele2, iloc, nodeminus, nodeplus, nodeplus1, node
         !stop 87
      end do Loop_Elements_Triangle

      return
    end subroutine Mapping_QTet_CV



    subroutine Renumbering_QTet_CV( x_nonods, x_ndgln2, &
         x, y, z )
      implicit none
      integer, intent( in ) :: x_nonods
      integer, dimension( x_nonods ), intent( inout ) :: x_ndgln2
      real, dimension( x_nonods ), intent( in ) :: x, y, z
      ! Local Variables
      real, parameter :: toler = 1.e-10
      integer :: ele, ele2, ntotel, node1, node2, icount
      integer, dimension( : ), allocatable :: x_temp
      logical :: ltest

      ewrite(3,*)'Renumbering and elliminating repeating nodes'

      allocate( x_temp( x_nonods ) )
      x_temp = 0
      x_temp = x_ndgln2

      ntotel = totele * tet_totele * no_of_nodes_in_faces + totele
      icount = 0

      Loop_Scan1: do ele = 1, ntotel - 1
         node1 = x_ndgln2( ele )
         Loop_Scan2: do ele2 = ele + 1, ntotel
            node2 = x_ndgln2( ele2 )
            Conditional_Coord: if( ( abs( x( node1 ) - x( node2 ) ) <= toler ) .and. &
                 ( abs(y( node1 ) - y( node2 ) ) <= toler ) .and. &
                 ( abs( z( node1 ) - z( node2 ) ) <= toler ) )then
               icount = icount + 1
               x_ndgln2( ele2 ) = x_ndgln2( ele )
               ewrite(3,*)  ele, x( node1 ), y( node1 ), z( node1 )
               ewrite(3,*)  ele2, x( node2 ), y( node2 ), z( node2 )
               ewrite(3,*)  ' ', icount
               !  ewrite(3,*)  ele, ele2, x_temp( ele ), x_ndgln2( ele2 ) 
            end if Conditional_Coord
         end do Loop_Scan2
         !  ewrite(3,*)  x_temp( ele ), x_ndgln2( ele )
      end do Loop_Scan1

      ewrite(3,*)  ' '
      ewrite(3,*)  'comparing the x_ndgln2 with the old one '
      ele2 = 0
      do ele = 1, ntotel
         ewrite(3,*) ele, x( x_temp( ele ) ), y( x_temp( ele ) ), z( x_temp( ele ) )
         ewrite(3,*) ele, x( x_ndgln2( ele ) ), y( x_ndgln2( ele ) ), z( x_ndgln2( ele ) )
         ltest = x_temp( ele ) == x_ndgln2( ele )
         if ( .not. ltest ) then
            !     ele2 = ele2 + 1
            ewrite(3,*), ele
         end if
      end do

      return
    end subroutine Renumbering_QTet_CV

    subroutine checking_QTet_CV( x_nonods, x_ndgln2 )
      implicit none
      integer, intent( in ) :: x_nonods 
      integer, dimension( x_nonods ), intent( in ) :: x_ndgln2
      ! Local variables
      integer :: ele, ele2, inode, iloc, node1, node2, node3, node4

      ewrite(3,*)'Checking QTet_CV nodes after renumbering'
      do ele = 1, totele
         do ele2 = 1, tet_totele
            inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
                 ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )
            do iloc = 1, tri_totele
               if( iloc /= tri_totele ) then
                  node1 = x_ndgln2( inode + iloc + 1 )
                  node2 = x_ndgln2( inode + iloc )
                  node3 = x_ndgln2( inode + 3 + iloc )
                  ewrite(3,*) ele, ele2, iloc, node1, node2, node3
               else
                  node1 = x_ndgln2( inode + iloc - 2 )
                  node2 = x_ndgln2( inode + iloc )
                  node3 = x_ndgln2( inode + 3 + iloc )
                  ewrite(3,*) ele, ele2, iloc, node1, node2, node3
               end if
            end do
            iloc = 4
            node1 = x_ndgln2( inode + 1 )
            node2 = x_ndgln2( inode + 2 )
            node3 = x_ndgln2( inode + 3 )
            node4 = x_ndgln2( inode + 3 + iloc )
            ewrite(3,*) ele, ele2, iloc, node1, node2, node3, node4
         end do
         ewrite(3,*) ele, x_ndgln2( ele * tet_totele * no_of_nodes_in_faces + ele )
      end do

      return
    end subroutine checking_QTet_CV


    subroutine Mapping_LinTet2Quad( ele, x_nonods, &
         x_ndgln2, x_ndgln )
      ! This subrt remaps the node points from a tetrahedron (15) into
      ! a hexahedron. Here we are applying the template node points allocation
      ! used in the linear tetrahedron subrt VOL_CV_TET_SHAPE for each of the
      ! tetrahedra at the large tetrahedron.
      implicit none
      integer, intent( in ) :: ele, x_nonods
      integer, dimension( x_nonods ), intent( in ) :: x_ndgln2
      integer, dimension( x_nonods ), intent( inout ) :: x_ndgln
      ! Local Variables
      integer :: elequad, ele2, inode, node1, node2, node3, node4, &
           node5, node6, node7, node8

!!!
!!! Quad-Element 1:
!!!
      ele2 = 1 
      inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
           ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )
      node1 = x_ndgln2( inode + 2 ) ! Node 1 --> 2
      node2 = x_ndgln2( inode + 5 ) ! Node 2 --> 20
      node5 = x_ndgln2( inode + 7 ) ! Node 8 --> 22
      node6 = x_ndgln2( inode + 4 ) ! Node 9 --> 19

      ele2 = 2 
      inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
           ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )
      node3 = x_ndgln2( inode + 5 ) ! Node 4 --> 24
      node8 = x_ndgln2( inode + 7 ) ! Node 12 --> 26

      ele2 = 4
      inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
           ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )
      node4 = x_ndgln2( inode + 7 ) ! Node 5 --> 34

      node7 = x_ndgln2( ele * tet_totele * no_of_nodes_in_faces + ele ) ! Node 11 --> 11

      elequad = 1 + ( ele - 1 ) * 4
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 1 ) = node1 ! 1 
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 2 ) = node2 ! 2  
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 3 ) = node3 ! 4  
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 4 ) = node4 ! 5  

      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 5 ) = node6 ! 9  
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 6 ) = node5 ! 8  
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 7 ) = node7 ! 11 
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 8 ) = node8 ! 12 

!!!
!!! Quad-Element 2:
!!!
      ele2 = 1
      inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
           ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )
      node1 = x_ndgln2( inode + 5 ) ! Node 2 --> 20
      node2 = x_ndgln2( inode + 3 ) ! Node 3 --> 3 
      node5 = x_ndgln2( inode + 7 ) ! Node 8 --> 22
      node6 = x_ndgln2( inode + 6 ) ! Node 10 --> 21

      ele2 = 2
      inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
           ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )
      node7 = x_ndgln2( inode + 7 ) ! Node 12 --> 26

      ele2 = 3
      inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
           ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )
      node4 = x_ndgln2( inode + 5 ) ! Node 6 --> 28
      node8 = x_ndgln2( inode + 7 ) ! Node 13 --> 30

      ele2 = 4
      inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
           ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )
      node3 = x_ndgln2( inode + 7 )  ! Node 5 --> 34

      elequad = 2 + ( ele - 1 ) * 4
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 1 ) = node1 ! 2
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 2 ) = node2 ! 3
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 3 ) = node3 ! 5
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 4 ) = node4 ! 6

      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 5 ) = node5 ! 8
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 6 ) = node6 ! 10
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 7 ) = node7 ! 12
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 8 ) = node8 ! 13

!!!
!!! Quad-Element 3:
!!!
      ele2 = 1
      inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
           ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )

      ele2 = 2
      inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
           ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )
      node1 = x_ndgln2( inode + 5 ) ! Node 4 --> 24 
      node3 = x_ndgln2( inode + 3 ) ! Node 7 --> 4
      node6 = x_ndgln2( inode + 7 ) ! Node 12 --> 26
      node7 = x_ndgln2( inode + 6 ) ! Node 14 --> 25

      ele2 = 3
      inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
           ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )
      node4 = x_ndgln2( inode + 5 ) ! Node 6 --> 28
      node8 = x_ndgln2( inode + 7 ) ! Node 13 --> 30

      ele2 = 4
      inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
           ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )
      node2 = x_ndgln2( inode + 7 ) ! Node 5 --> 34

      node5 = x_ndgln2( ele * tet_totele * no_of_nodes_in_faces + ele ) ! Node 11 --> 11


      elequad = 3 + ( ele - 1 ) * 4
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 1 ) = node1 ! 4
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 2 ) = node2 ! 5
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 3 ) = node3 ! 7
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 4 ) = node4 ! 6

      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 5 ) = node5 ! 11
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 6 ) = node6 ! 12
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 7 ) = node7 ! 14
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 8 ) = node8 ! 13

      ! Quad-Element 4:

      ele2 = 1
      inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
           ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )
      node1 = x_ndgln2( inode + 4 ) ! Node 9 --> 19
      node2 = x_ndgln2( inode + 7 ) ! Node 8 --> 22
      node5 = x_ndgln2( inode + 1 ) ! Node 15 --> 1
      node6 = x_ndgln2( inode + 6 ) ! Node 10 --> 21

      ele2 = 2
      inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
           ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )
      node4 = x_ndgln2( inode + 7 ) ! Node 12 --> 26
      node7 = x_ndgln2( inode + 6 ) ! Node 14 --> 25

      ele2 = 3
      inode = ( ele - 1 ) * tet_totele * no_of_nodes_in_faces + & 
           ( ele2 - 1 ) * no_of_nodes_in_faces + ( ele - 1 )
      node8 = x_ndgln2( inode + 7 ) ! Node 13 --> 30

      node3 = x_ndgln2( ele * tet_totele * no_of_nodes_in_faces + ele ) ! Node 11 --> 11


      elequad = 4 + ( ele - 1 ) * 4
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 1 ) = node1 ! 9
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 2 ) = node2 ! 8
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 3 ) = node3 ! 11
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 4 ) = node4 ! 12

      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 5 ) = node5 ! 15
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 6 ) = node6 ! 10
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 7 ) = node7 ! 14
      x_ndgln( ( elequad - 1 ) * quad_cv_nloc + 8 ) = node8 ! 13

      ! Printing the hexahedra node points
      ewrite(3,*) 'Printing node points of the hexaedra from tetrahedron', ele 
      do ele2 = 1, 4
         ewrite(3,*) ele2
         ewrite(3,'(8i3)') ( x_ndgln( ( ele2 - 1 ) * quad_cv_nloc + inode ), inode = 1, 8 )
      end do

      return
    end subroutine Mapping_LinTet2Quad



    real function QTet_Height( ilayer, x_nonods, z )
      implicit none
      integer, intent( in ) :: ilayer, x_nonods
      real, dimension( x_nonods ), intent( in ) :: z

      QTet_Height = 0.
      QTet_Height = real ( nlayer - ilayer ) / real( nlayer ) * &
           ( z( tet_nodes( 0, 0 ) ) + z( tet_nodes( 0, nlayer ) ) )

      return
    end function QTet_Height

    logical function QTet_CheckRepetitiveNodes( n, x )
      implicit none
      integer, intent( in ) :: n
      integer, dimension( n ), intent( in ) :: x
      ! Local variables
      integer :: i, j, temp

      ewrite(3,*) ' checking repetitive nodes'
      QTet_CheckRepetitiveNodes = .false.
      do i = 1, n - 1
         if( x( i ) /= 0 ) then
            temp = x( i )
            do j = i + 1 , n
               QTet_CheckRepetitiveNodes = delta_diff( temp, x( j ) )
               if( QTet_CheckRepetitiveNodes ) then
                  ewrite(3,*)  i, j, x( i ), x( j )
                  return
               end if
            end do
         end if

      end do

      return

    end function QTet_CheckRepetitiveNodes

    logical function delta_diff( a, b )
      ! If abs( a-b ) < delta, a == b and the function returns .true.
      implicit none
      integer :: a, b
      real, parameter :: delta = 1.e-6
      real :: c

      delta_diff = .false.
      if( abs( real( a - b ) ) < delta ) delta_diff = .true.

      return
    end function delta_diff

  end module pascal_tetrahedra
