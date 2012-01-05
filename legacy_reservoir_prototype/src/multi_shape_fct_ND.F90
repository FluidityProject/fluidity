  
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
  
  
!!!!==============================================!!!!
!!!!     SHAPE FUNCTIONS SUBRTS FOR MULTI-D       !!!!
!!!! (Quadrilaterals, Triangles, Hexaedra and     !!!!
!!!!             Tetrahedra)                      !!!!
!!!!==============================================!!!!
  

  module shape_functions_Linear_Quadratic

    use fldebug
!    use shape_functions
  use state_module
  use spud
  use global_parameters, only: option_path_len

  contains

    subroutine re2dn4( lowqua, ngi, nloc, mloc, &
         m, weight, n, nlx, nly, &
         sngi, snloc, sweigh, sn, snlx )
      ! This subrt computes shape functions M and N and their derivatives
      ! at the Gauss points.
      ! NB: We may need to define surface elements for p and (u,v,w) 
      implicit none
      logical, intent( in ) :: lowqua
      integer, intent( in ) :: ngi, nloc, mloc
      real, dimension( mloc, ngi ), intent( inout ) :: m
      real, dimension( ngi ), intent( inout ) :: weight
      real, dimension( nloc, ngi ), intent( inout ) :: n, nlx, nly
      integer, intent( in ) :: sngi, snloc
      real, dimension( sngi ), intent( inout ) :: sweigh
      real, dimension( snloc, sngi ), intent( inout ) :: sn, snlx
      ! Local variables:
      integer, parameter :: nl = 16, nlp = 4, npq = 2
      real, dimension( : ), allocatable :: lx, ly, lxp, lyp, weit
      integer :: q, p, gpoi, corn, ndgi, i
      real :: posi
      logical :: getndp

      ! Allocating memory
      allocate( lx( nl ) )
      allocate( ly( nl ) )
      allocate( lxp( nlp ) )
      allocate( lyp( nlp ) )

      ! LX/YP(I) are the local X and Y coordinates of nodal point I
      lxp( 1 ) = -1
      lyp( 1 ) = -1
      lxp( 2 ) = 1
      lyp( 2 ) = -1
      lxp( 3 ) = -1
      lyp( 3 ) = 1
      lxp( 4 ) = 1
      lyp( 4 ) = 1

      Conditional_NGI: if( ngi == 4 ) then
         posi = 1. / sqrt( 3. )
         lx( 1 ) = -posi
         ly( 1 ) = -posi
         lx( 2 ) = posi
         ly( 2 ) = posi
         Loop_Q: do q = 1, npq
            Loop_P: do p = 1, npq
               Loop_Corn: do corn = 1, nlp
                  gpoi = ( q - 1 ) * 2 + p
                  if( mloc == 1 ) m( 1, gpoi ) = 1.
                  weight( gpoi ) = 1.
                  n( corn, gpoi ) = 0.25 * ( 1. + lxp( corn ) * lx( p ) ) * &
                       ( 1. + lyp( corn ) * ly( q ) )
                  nlx( corn, gpoi ) = 0.25 * lxp( corn ) * ( 1. + lyp( corn ) * ly( q ))
                  nly( corn, gpoi ) = 0.25 * lyp( corn ) * ( 1. + lxp( corn ) * lx( q ))
               end do Loop_Corn
            end do Loop_P
         end do Loop_Q

         if( ( sngi > 1 ) .and. ( snloc > 1 ) ) then ! Surface shape functions
            Loop_P2: do p = 1, npq
               Loop_Corn2: do corn = 1, npq
                  gpoi = p
                  sn( corn, gpoi ) = 0.5 * ( 1. + lxp( corn ) * lx( p ) )
                  snlx( corn, gpoi ) = 0.5 * lxp( corn )
                  sweigh( gpoi ) = 0.
               end do Loop_Corn2
            end do Loop_P2
         end if

      else ! If ngi =/ 4
         ndgi = int( sqrt( ngi + 0.1 ) + 0.1 )
         getndp = .false.
         call lagrot( weit, lx, ndgi, getndp )   
         ly( 1 : ndgi ) = lx( 1 : ndgi )
         Loop_Q3: do q = 1, ndgi
            Loop_P3: do p = 1, ndgi
               Loop_Corn3: do corn = 1, nlp
                  gpoi = ( q - 1 ) * ndgi + p
                  if( mloc == 1 ) m( 1, gpoi ) = 1.
                  weight( gpoi ) = weit( p ) * weit( q )
                  n( corn, gpoi ) = 0.25 * ( 1. + lxp( corn ) * lxp( p ) * &
                       ( 1. + lyp( corn ) * ly( q ) ) )
                  nlx( corn, gpoi ) = 0.25 * lxp( corn ) * &
                       ( 1. + lyp( corn ) * ly( q ) )
                  nly( corn, gpoi ) = 0.25 * lyp( corn ) * &
                       ( 1. + lxp( corn ) * lx( q ) )
               end do Loop_Corn3
            end do Loop_P3
         end do Loop_Q3

         if( sngi > 0 ) then
            getndp = .false.
            call lagrot( weit, lx, sngi, getndp )
            Loop_P4: do p = 1, sngi
               Loop_Corn4: do corn = 1, npq
                  gpoi = p
                  sn( corn, gpoi ) = 0.5 * ( 1. + lxp( corn ) * lx( p ) )
                  snlx( corn, gpoi ) = 0.5 * lxp( corn )
                  sweigh( gpoi ) = weit( p )
               end do Loop_Corn4
            end do Loop_P4
         end if

      end if Conditional_NGI

      if( mloc == nloc ) then
         do i = 1, nlp
            do corn = 1, nlp
               m( corn, i ) = n( corn, i )
            end do
         end do
      end if

      ! Deallocating
      deallocate( lx )
      deallocate( ly )
      deallocate( lxp )
      deallocate( lyp )

      return
    end subroutine re2dn4

    subroutine re3dn8( lowqua, ngi, nloc, mloc, &
         m, weight, n, nlx, nly, nlz, &
         sngi, snloc, sweigh, sn, snlx, snly ) 
      ! This subrt. computes the shape functions M and N and their
      ! derivatives at the Gauss points for 3D.
      ! If LOWQUA, then use one point quadrature else use 8 point quadrature.
      !NB.: LX/YP(I) are the local X/Y coordinates of nodal point I.
      implicit none
      logical, intent( in ) :: lowqua
      integer, intent( in ) :: ngi, nloc, mloc
      real, dimension( mloc, ngi ), intent( inout ) :: m
      real, dimension( ngi ), intent( inout ) :: weight
      real, dimension( nloc, ngi ), intent( inout ) :: n, nlx, nly, nlz
      integer, intent( in ) :: sngi, snloc
      real, dimension( sngi ), intent( inout ) :: sweigh
      real, dimension( snloc, sngi ), intent( inout ) :: sn, snlx, snly
      ! Local variables:
      integer, parameter :: nl = 4, nlp = 8, npq = 2
      real, dimension( : ), allocatable :: lx, ly, lz, lxp, lyp, lzp, weit
      integer :: p, q, ir, corn, gpoi, ngi1d, gi
      real :: posi

      ! Allocating memory
      allocate( lx( nl ) )
      allocate( ly( nl ) )
      allocate( lz( nl ) )
      allocate( lxp( nlp ) )
      allocate( lyp( nlp ) )
      allocate( lzp( nlp ) )
      allocate( weit( nl ) )


      if( sngi >= 1 ) then ! Surface integrals
         call re2dn4( lowqua, sngi, snloc, 0, &
              m, sweigh, sn, snlx, snly, &
              0, 0, sweigh, sn, snlx )
      end if

      ! Nodal Point 1
      lxp( 1 ) = -1
      lyp( 1 ) = -1
      lzp( 1 ) = -1
      ! Nodal Point 2
      lxp( 2 ) = 1
      lyp( 2 ) = -1
      lzp( 2 ) = -1
      ! Nodal Point 3
      lxp( 3 ) = -1
      lyp( 3 ) = 1
      lzp( 3 ) = -1
      ! Nodal Point 4
      lxp( 4 ) = 1
      lyp( 4 ) = 1
      lzp( 4 ) = -1
      ! Nodal Point 5
      lxp( 5 ) = -1
      lyp( 5 ) = -1
      lzp( 5 ) = 1
      ! Nodal Point 6
      lxp( 6 ) = 1
      lyp( 6 ) = -1
      lzp( 6 ) = 1
      ! Nodal Point 7
      lxp( 7 ) = -1
      lyp( 7 ) = 1
      lzp( 7 ) = 1
      ! Nodal Point 8 
      lxp( 8 ) = 1
      lyp( 8 ) = 1
      lzp( 8 ) = 1

      Conditional_NGI: Select Case( NGI )

      case( 8 )

         posi = 1. / sqrt( 3. )
         lx( 1 ) = -posi
         ly( 1 ) = -posi
         lz( 1 ) = -posi
         lx( 2 ) = posi
         ly( 2 ) = posi
         lz( 2 ) = posi
         Loop_P1: do p = 1, npq
            Loop_Q1: do q = 1, npq
               Loop_IR1: do ir = 1, npq
                  Loop_Corn1: do corn = 1, nlp
                     gpoi = ( ( p - 1 ) * 2 + 1 ) + ( ir - 1 ) * 4
                     if( mloc > 0 ) m( 1, gpoi ) = 1.
                     ! Weight
                     weight( gpoi ) = 1.
                     ! N
                     n( corn, gpoi ) = 0.125 * ( 1. + lxp( corn ) * lx( p ) ) * &
                          ( 1. + lyp( corn ) * ly( q ) ) * ( 1. + lzp( corn ) * &
                          lz( ir ) )
                     ! x-derivative
                     nlx( corn, gpoi ) = 0.125 * lxp( corn ) * ( 1. + &
                          lyp( corn ) * ly( q ) ) * ( 1. + lzp( corn ) * lz( ir ) ) 
                     ! y-derivative
                     nly( corn, gpoi ) = 0.125 * lyp( corn ) * ( 1. + lxp( corn ) * &
                          lx( p ) ) * ( 1. + lzp( corn ) * lz( ir ) )
                     ! z-derivative
                     nlz( corn, gpoi ) = 0.125 * lzp( corn ) * ( 1. + lxp( corn ) * &
                          lx( p ) ) * ( 1. + lyp( corn ) * ly( q ) )
                  end do Loop_Corn1
               end do Loop_IR1
            end do Loop_Q1
         end do Loop_P1

      case default

         ngi1d = int( ( ngi + 0.1 ) ** ( 1. / 3. ) + 0.1 )
         do gi = 1, ngi1d
            weit( gi ) = rgptwe( gi, ngi1d, .true. )
            lx( gi ) = rgptwe( gi, ngi1d, .false. )
            ly( gi ) = lx( gi )
            lz( gi ) = lx( gi )
         end do
         Loop_P2: do p = 1, ngi1d
            Loop_Q2: do q = 1, ngi1d
               Loop_IR2: do ir = 1, ngi1d
                  Loop_Corn2: do corn = 1, nlp
                     gpoi = ( ( p - 1 ) * ngi1d + q ) + ( ir - 1 ) * ngi1d **2
                     if( mloc > 0 ) m( 1, gpoi ) = 1.
                     ! Weight
                     weight( gpoi ) = weit( p ) * weit( q ) * weit( ir )
                     ! N
                     n( corn, gpoi ) = 0.125 * ( 1. + lxp( corn ) * lx( p ) ) * &
                          ( 1. + lyp( corn ) * ly( q ) ) * ( 1. + lzp( corn ) * lz( ir ) )
                     ! x-derivative
                     nlx( corn, gpoi ) = 0.125 * lxp( corn ) * ( 1. + lyp( corn ) * &
                          ly( q ) ) * ( 1. + lzp( corn ) * lz( ir ) )
                     ! y-derivative
                     nly( corn, gpoi ) = 0.125 * lyp( corn ) * ( 1. + lxp( corn ) * &
                          lx( p ) ) * ( 1. + lzp( corn ) * lz( ir ) )
                     ! z-derivative
                     nlz( corn, gpoi ) = 0.125 * lzp( corn ) * ( 1. + lxp( corn ) * &
                          lx( p ) ) * ( 1. + lyp( corn ) * ly( q ) )
                  end do Loop_Corn2
               end do Loop_IR2
            end do Loop_Q2
         end do Loop_P2

      end Select Conditional_NGI

      ! Deallocating
      deallocate( lx )
      deallocate( ly )
      deallocate( lz )
      deallocate( lxp )
      deallocate( lyp )
      deallocate( lzp )
      deallocate( weit )

      return
    end subroutine re3dn8


    subroutine re2dn9( lowqua, ngi, nloc, mloc, &
         m, weight, n, nlx, nly )
      ! Quadratic variation (2D) for velocity -- 9 node brick element.
      ! Linear variation (2D) for pressure -- 4 node brick element.
      ! NB.: We may need to define surface elements for p and (u,v,w). 
      !
      ! This is for the 2-D 27node element, which is number as follows
      !   |  /Z
      !  Y| /
      !   |/
      !   +-------X
      ! For Z = -1 ...
      !   7   8   9
      !   4   5   6
      !   1   2   3
      ! For M the shape functions have local node numbers
      ! For Z = -1 ...
      !    4       3
      !
      !    1       2
      implicit none
      logical, intent( in ) :: lowqua 
      integer, intent( in ) :: ngi, nloc, mloc
      real, dimension( mloc, ngi ), intent( inout ) :: m
      real, dimension( ngi ), intent( inout ) :: weight
      real, dimension( nloc, ngi ), intent( inout ) :: n, nlx, nly
      ! Local variables:
      integer, parameter :: nl = 3, nlp = 9, npq = 4
      real, dimension( : ), allocatable :: lx, ly, lz, lxp, lyp, weit, &
           xn, yn, zn, dxn, dyn, dzn
      integer :: nquad, p, q, corn, gpoi, ilx, ily, nj
      real :: posi

      ! Allocating memory
      allocate( lx( nl ) )
      allocate( ly( nl ) )
      allocate( lz( nl ) )
      allocate( lxp( nlp ) )
      allocate( lyp( nlp ) )
      allocate( weit( nl ) )
      allocate( xn( nl ) )
      allocate( yn( nl ) )
      allocate( zn( nl ) )
      allocate( dxn( nl ) )
      allocate( dyn( nl ) )
      allocate( dzn( nl ) )

      Conditional_LOWQUA: if( lowqua ) then
         posi = 1. / sqrt( 3. )
         lx( 1 ) = -posi
         ly( 1 ) = -posi
         lz( 1 ) = -posi
         lx( 2 ) = posi
         ly( 2 ) = posi
         lz( 2 ) = posi
         weit( 1 ) = 1.
         weit( 2 ) = 1.
         weit( 3 ) = 1.
         nquad = 2
      else
         posi = 0.774596669241483
         lx( 1 ) = -posi
         ly( 1 ) = -posi
         lx( 2 ) = 0.
         ly( 2 ) = 0.
         lx( 3 ) = posi
         ly( 3 ) = posi
         weit( 1 ) = 0.555555555555556
         weit( 2 ) = 0.888888888888889
         weit( 3 ) = 0.555555555555556
         nquad = 3
      end if Conditional_LOWQUA

      ! Nodal Point 1
      lxp( 1 ) = -1
      lyp( 1 ) = -1
      ! Nodal Point 2
      lxp( 2 ) = 1
      lyp( 2 ) = -1
      ! Nodal Point 3
      lxp( 3 ) = 1
      lyp( 3 ) = 1
      ! Nodal Point 4
      lxp( 4 ) = -1
      lyp( 4 ) = 1

      ! Compute M
      Loop_P1: do p = 1, nquad
         Loop_Q1: do q = 1, nquad
            Loop_Corn1: do corn = 1, npq
               gpoi = ( p - 1 ) * nquad + q
               m( corn, gpoi ) = 0.25 * ( 1. + lxp( corn ) * lx( p ) ) * &
                    ( 1 + lyp( corn ) * ly( q ) )
            end do Loop_Corn1
         end do Loop_Q1
      end do Loop_P1

      ! Compute N
      Loop_ILX1: do ilx = 1, nl
         Loop_ILY1: do ily = 1, nl
            nj = ilx + ( ily - 1 ) * 3
            lxp( nj ) = real( ilx - 2 )
            lyp( nj ) = real( ily - 2 )
         end do Loop_ILY1
      end do Loop_ILX1

      Loop_P2: do p = 1, nquad
         Loop_Q2: do q = 1, nquad
            gpoi = ( p - 1 ) * nquad + q
            ! Weight
            weight( gpoi ) = weit( p ) * weit( q )
            ! XN
            xn( 1 ) = 0.5 * lx( p ) * ( lx( p ) - 1. )
            xn( 2 ) = 1. - lx( p ) * lx( p )
            xn( 3 ) = 0.5 * lx( p ) * ( lx( p ) + 1. )
            ! DXDN
            dxn( 1 ) = 0.5 * ( 2. * lx( p ) - 1. )
            dxn( 2 ) = -2. * lx( p )
            dxn( 3 ) = 0.5 * ( 2. * lx( p ) + 1. )
            ! YN
            yn( 1 ) = 0.5 * ly( q ) * ( ly( q ) - 1. )
            yn( 2 ) = 1. - ly( q ) * ly( q )
            yn( 3 ) = 0.5 * ly( q ) * ( ly( q ) + 1. )
            ! DYDN
            dyn( 1 ) = 0.5 * ( 2. * ly( q ) - 1. )
            dyn( 2 ) = -2. * ly( q )
            dyn( 3 ) = 0.5 * ( 2. * ly( q ) + 1. )
            ! N, NLX and NLY
            Loop_ILX2: do ilx = 1, nl
               Loop_ILY2: do ily = 1, nl
                  nj = ilx + ( ily - 1 ) * 3
                  n( nj, gpoi ) = xn( ilx ) * yn( ily )
                  nlx( nj, gpoi ) = dxn( ilx ) * yn( ily )
                  nly( nj, gpoi ) = xn( ilx ) * dyn( ily )
               end do Loop_ILY2
            end do Loop_ILX2
         end do Loop_Q2
      end do Loop_P2

      ! Deallocating
      deallocate( lx )
      deallocate( ly )
      deallocate( lz )
      deallocate( lxp )
      deallocate( lyp )
      deallocate( weit )
      deallocate( xn )
      deallocate( yn )
      deallocate( zn )
      deallocate( dxn )
      deallocate( dyn )
      deallocate( dzn )

      return
    end subroutine re2dn9

    subroutine re3d27( lowqua, ngi, nloc, mloc, &
         m, weight, n, nlx, nly, nlz )
      ! Quadratic variation (3D) for velocity -- 27 node brick element.
      ! Linear variation (3D) for pressure -- 8 node brick element.
      ! NB.: We may need to define surface elements for p and (u,v,w). 
      !
      ! This is for the 2-D 27node element, which is number as follows
      !   |  /Z
      !  Y| /
      !   |/
      !   +-------X
      ! For Z = -1 ...
      !   7   8   9
      !   4   5   6
      !   1   2   3
      ! For Z = 0 ...
      !   16   17   18
      !   13   14   15
      !   10   11   12
      ! For Z = 1 ...
      !   25    26   27
      !   22    23   24
      !   19    20   21
      ! For M the shape functions have local node numbers
      ! For Z = -1 ...
      !    4       3
      !
      !    1       2
      ! For Z = 1 ...
      !    8       7
      !
      !    5       6
      implicit none
      logical, intent( in ) :: lowqua 
      integer, intent( in ) :: ngi, nloc, mloc
      real, dimension( mloc, ngi ), intent( inout ) :: m
      real, dimension( ngi ), intent( inout ) :: weight
      real, dimension( nloc, ngi ), intent( inout ) :: n, nlx, nly, nlz
      ! Local variables:
      integer, parameter :: nl = 3, nlp = 27, npq = 4
      real, dimension( : ), allocatable :: lx, ly, lz, lxp, lyp, lzp, weit, &
           xn, yn, zn, dxn, dyn, dzn
      integer :: nquad, p, q, ir, corn, gpoi, ilx, ily, ilz, nj
      real :: posi

      ! Allocating memory
      allocate( lx( nl ) )
      allocate( ly( nl ) )
      allocate( lz( nl ) )
      allocate( lxp( nlp ) )
      allocate( lyp( nlp ) )
      allocate( lzp( nlp ) )
      allocate( weit( nl ) )
      allocate( xn( nl ) )
      allocate( yn( nl ) )
      allocate( zn( nl ) )
      allocate( dxn( nl ) )
      allocate( dyn( nl ) )
      allocate( dzn( nl ) )

      Conditional_LOWQUA: if( lowqua ) then
         posi = 1. / sqrt( 3. )
         lx( 1 ) = -posi
         ly( 1 ) = -posi
         lz( 1 ) = -posi
         lx( 2 ) = posi
         ly( 2 ) = posi
         lz( 2 ) = posi
         weit( 1 ) = 1.
         weit( 2 ) = 1.
         weit( 3 ) = 1.
         nquad = 2
      else
         posi = 0.774596669241483
         lx( 1 ) = -posi
         ly( 1 ) = -posi
         lz( 1 ) = -posi
         lx( 2 ) = 0.
         ly( 2 ) = 0.
         lz( 2 ) = 0.
         lx( 3 ) = posi
         ly( 3 ) = posi
         lz( 3 ) = posi
         weit( 1 ) = 0.555555555555556
         weit( 2 ) = 0.888888888888889
         weit( 3 ) = 0.555555555555556
         nquad = 3
      end if Conditional_LOWQUA

      ! Nodal Point 1
      lxp( 1 ) = -1
      lyp( 1 ) = -1
      lzp( 1 ) = -1
      ! Nodal Point 2
      lxp( 2 ) = 1
      lyp( 2 ) = -1
      lzp( 2 ) = -1
      ! Nodal Point 3
      lxp( 3 ) = 1
      lyp( 3 ) = 1
      lzp( 3 ) = -1
      ! Nodal Point 4
      lxp( 4 ) = -1
      lyp( 4 ) = 1
      lzp( 4 ) = -1
      ! Nodal Point 5
      lxp( 5 ) = -1
      lyp( 5 ) = -1
      lzp( 5 ) = 1
      ! Nodal Point 6
      lxp( 6 ) = 1
      lyp( 6 ) = -1
      lzp( 6 ) = 1
      ! Nodal Point 7
      lxp( 7 ) = 1
      lyp( 7 ) = 1
      lzp( 7 ) = 1
      ! Nodal Point 8
      lxp( 8 ) = -1
      lyp( 8 ) = 1
      lzp( 8 ) = 1

      ! Compute M
      Loop_P1: do p = 1, nquad
         Loop_Q1: do q = 1, nquad
            Loop_IR1: do ir = 1, nquad
               Loop_Corn1: do corn = 1, npq
                  gpoi = ( p - 1 ) * nquad + q + ( ir - 1 ) * nquad * nquad
                  m( corn, gpoi ) = 0.125 * ( 1. + lxp( corn ) * lx( p ) ) * &
                       ( 1 + lyp( corn ) * ly( q ) * ( 1. + lzp( corn ) * lz( ir ) ) ) 
               end do Loop_Corn1
            end do Loop_IR1
         end do Loop_Q1
      end do Loop_P1


      ! Compute N
      Loop_ILX1: do ilx = 1, nl
         Loop_ILY1: do ily = 1, nl
            Loop_ILZ1: do ilz = 1, nl
               nj = ilx + ( ily - 1 ) * 3 + ( ilz - 1 ) * 9
               lxp( nj ) = real( ilx - 2 )
               lyp( nj ) = real( ily - 2 )
               lzp( nj ) = real( ilz - 2 )
            end do Loop_ILZ1
         end do Loop_ILY1
      end do Loop_ILX1

      Loop_P2: do p = 1, nquad
         Loop_Q2: do q = 1, nquad
            Loop_IR2: do ir = 1, nquad
               gpoi = ( p - 1 ) * nquad + q + ( ir - 1 ) * nquad * nquad
               ! Weight
               weight( gpoi ) = weit( p ) * weit( q ) * weit( ir )
               ! XN
               xn( 1 ) = 0.5 * lx( p ) * ( lx( p ) - 1. )
               xn( 2 ) = 1. - lx( p ) * lx( p )
               xn( 3 ) = 0.5 * lx( p ) * ( lx( p ) + 1. )
               ! DXDN
               dxn( 1 ) = 0.5 * ( 2. * lx( p ) - 1. )
               dxn( 2 ) = -2. * lx( p )
               dxn( 3 ) = 0.5 * ( 2. * lx( p ) + 1. )
               ! YN
               yn( 1 ) = 0.5 * ly( q ) * ( ly( q ) - 1. )
               yn( 2 ) = 1. - ly( q ) * ly( q )
               yn( 3 ) = 0.5 * ly( q ) * ( ly( q ) + 1. )
               ! DYDN
               dyn( 1 ) = 0.5 * ( 2. * ly( q ) - 1. )
               dyn( 2 ) = -2. * ly( q )
               dyn( 3 ) = 0.5 * ( 2. * ly( q ) + 1. )
               ! ZN
               zn( 1 ) = 0.5 * lz( ir ) * ( lz( ir ) - 1. )
               zn( 2 ) = 1. - lz( ir ) * lz( ir )
               zn( 3 ) = 0.5 * lz( q ) * ( ly( q ) + 1. )
               ! DZDN
               dzn( 1 ) = 0.5 * ( 2. * lz( ir ) - 1. )
               dzn( 2 ) = -2. * lz( ir )
               dzn( 3 ) = 0.5 * ( 2. * lz( ir ) + 1. )
               ! N, NLX and NLY
               Loop_ILX2: do ilx = 1, nl
                  Loop_ILY2: do ily = 1, nl
                     Loop_ILZ2: do ilz = 1, nl
                        nj = ilx + ( ily - 1 ) * 3 + ( ilz - 1 ) * 9
                        n( nj, gpoi ) = xn( ilx ) * yn( ily ) * zn( ilz )
                        nlx( nj, gpoi ) = dxn( ilx ) * yn( ily ) * zn( ilz )
                        nly( nj, gpoi ) = xn( ilx ) * dyn( ily ) * zn( ilz )
                        nlz( nj, gpoi ) = xn( ilx ) * yn( ily ) * dzn( ilz )
                     end do Loop_ILZ2
                  end do Loop_ILY2
               end do Loop_ILX2
            end do Loop_IR2
         end do Loop_Q2
      end do Loop_P2


      ! Deallocating
      deallocate( lx )
      deallocate( ly )
      deallocate( lz )
      deallocate( lxp )
      deallocate( lyp )
      deallocate( weit )
      deallocate( xn )
      deallocate( yn )
      deallocate( zn )
      deallocate( dxn )
      deallocate( dyn )
      deallocate( dzn )

      return
    end subroutine re3d27 


    subroutine retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
         cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface )
      implicit none
      integer, intent( in ) :: ndim, cv_ele_type, cv_nloc, u_nloc
      integer, intent( inout ) :: cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface
      ! Local variables
      character( len = option_path_len ) :: overlapping_path


      Conditional_EleType: Select Case( cv_ele_type )

      case( 1, 2 ) ! 1D
         Conditional_CV_NLOC1D: Select Case( cv_nloc )
         case( 1 )
            cv_ngi = 1
            scvngi = 2
         case( 2 )
            cv_ngi = 12
            scvngi = 3
         case( 3 )
            cv_ngi = 12
            scvngi = 4
            if( ( u_nloc == 4 ) .or. ( u_nloc == 5 )) cv_ngi = 18
         case default; FLExit(" Invalid integer for cv_nloc ")
         end Select Conditional_CV_NLOC1D
         sbcvngi = 1
         nface = 2

      case( 3, 4 ) ! Triangles
         Conditional_CV_NLOC2D_Tri: Select Case( cv_nloc )
         case( 3 ) ! Linear Triangle
            cv_ngi = 12
            scvngi = 3
            sbcvngi = 2 
         case( 6 ) ! Quadratic Triangle
            cv_ngi = 81
            scvngi = 18
            sbcvngi = 6 
         case default; FLExit(" Invalid integer for cv_nloc ")
         end Select Conditional_CV_NLOC2D_Tri
         nface = 3

      case( 5, 6 ) ! Quads
         Conditional_CV_NLOC2D_Quad: Select Case( cv_nloc )
         case( 4 ) ! Bi-linear Quad
            cv_ngi = 16
            scvngi = 2
            sbcvngi = 8 
         case( 9 ) ! Bi-quad Quad
            cv_ngi = 36
            scvngi = 24
            sbcvngi = 6
         case default; FLExit(" Invalid integer for cv_nloc ")
         end Select Conditional_CV_NLOC2D_Quad
         nface = 4

      case( 7, 8 ) ! Tetrahedra
         Conditional_CV_NLOC3D_Tet: Select Case( cv_nloc )
         case( 4 ) ! Linear 
            cv_ngi = 32
            scvngi = 6
            sbcvngi = 3
         case( 10 ) ! Quadratic (double check this)
            cv_ngi = 128
            scvngi = 96
            sbcvngi = 54
         case default; FLExit(" Invalid integer for cv_nloc ")
         end Select Conditional_CV_NLOC3D_Tet
         nface = 4

      case( 9, 10 )
         Conditional_CV_NLOC3D_Hex: Select Case( cv_nloc )
         case( 8 ) ! Tri-linear Hex
            cv_ngi = 64
            scvngi = 12
            sbcvngi = 4
         case( 27 ) ! Tri-quad Hex
            cv_ngi = 216
            scvngi = 216
            sbcvngi = 36
         case default; FLExit(" Invalid integer for cv_nloc ")
         end Select Conditional_CV_NLOC3D_Hex
         nface = 6

      case default; FLExit( " Invalid integer for cv_ele_type " )

      end Select Conditional_EleType

      if( cv_ele_type > 2 ) scvngi = scvngi + nface * sbcvngi
      cv_ngi_short = cv_ngi

      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) cv_ngi = cv_ngi * cv_nloc


      return
    end subroutine retrieve_ngi


    subroutine lagrot( weit, quadpos, ndgi, getndp )
      ! This computes the weight and points for standard Gaussian quadrature.
      ! If (GETNDP == T) then get the position of the nodes and neglect the weights.
      implicit none
      integer, intent( in ) :: ndgi
      real, dimension( ndgi ), intent( inout ) :: weit, quadpos
      logical, intent( in ) :: getndp
      ! Local variables
      logical :: weight
      integer :: ig

      ewrite(3,*) 'In LAGROT'

      Conditional_GETNDP: if( .not. getndp ) then
         weight = .true.
         do ig = 1, ndgi
            weit( ig ) = rgptwe( ig, ndgi, weight )
         end do

         weight = .false.
         do ig = 1, ndgi
            quadpos( ig ) = rgptwe( ig, ndgi, weight )
         end do

      else
         if( ndgi == 1 ) then
            quadpos( ndgi ) = 0.
         else
            do ig = 1, ndgi
               quadpos( ig ) = -1. + 2. * real( ig - 1 ) / real( ndgi - 1 )
            end do
         end if

      end if Conditional_GETNDP

      ewrite(3,*) 'Leaving LAGROT'

      return
    end subroutine lagrot




    SUBROUTINE DETNLXR( ELE, X,Y,Z, XONDGL, TOTELE, NONODS, NLOC, NGI, &
         N, NLX, NLY, NLZ, WEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
         NX, NY, NZ) 
      IMPLICIT NONE
      INTEGER, intent( in ) :: ELE, TOTELE, NONODS, NLOC, NGI
      INTEGER, DIMENSION( TOTELE * NLOC ) :: XONDGL
      REAL, DIMENSION( NONODS ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( NLOC, NGI ), intent( in ) :: N, NLX, NLY, NLZ 
      REAL, DIMENSION( NGI ), intent( in ) :: WEIGHT
      REAL, DIMENSION( NGI ), intent( inout ) :: DETWEI, RA
      REAL, intent( inout ) :: VOLUME
      LOGICAL, intent( in ) :: D1, D3, DCYL
      REAL, DIMENSION( NLOC, NGI ), intent( inout ) :: NX, NY, NZ
      ! Local variables
      REAL, PARAMETER :: PIE = 3.141592654
      REAL :: AGI, BGI, CGI, DGI, EGI, FGI, GGI, HGI, KGI, A11, A12, A13, A21, &
           A22, A23, A31, A32, A33, DETJ, TWOPIE, RGI
      INTEGER :: GI, L, IGLX
      !
      VOLUME = 0.
      !
      IF(D3) THEN
         do  GI=1,NGI! Was loop 331
            !
            AGI=0.
            BGI=0.
            CGI=0.
            !
            DGI=0.
            EGI=0.
            FGI=0.
            !
            GGI=0.
            HGI=0.
            KGI=0.
            !
            do  L=1,NLOC! Was loop 79
               IGLX=XONDGL((ELE-1)*NLOC+L)
               ! NB R0 does not appear here although the z-coord might be Z+R0. 
               AGI=AGI+NLX(L,GI)*X(IGLX) 
               BGI=BGI+NLX(L,GI)*Y(IGLX) 
               CGI=CGI+NLX(L,GI)*Z(IGLX) 
               !
               DGI=DGI+NLY(L,GI)*X(IGLX) 
               EGI=EGI+NLY(L,GI)*Y(IGLX) 
               FGI=FGI+NLY(L,GI)*Z(IGLX) 
               !
               GGI=GGI+NLZ(L,GI)*X(IGLX) 
               HGI=HGI+NLZ(L,GI)*Y(IGLX) 
               KGI=KGI+NLZ(L,GI)*Z(IGLX)
            end do ! Was loop 79
            !
            DETJ=AGI*(EGI*KGI-FGI*HGI)&
                 -BGI*(DGI*KGI-FGI*GGI)&
                 +CGI*(DGI*HGI-EGI*GGI)
            DETWEI(GI)=DETJ*WEIGHT(GI)
            RA(GI)=1.0
            VOLUME=VOLUME+DETWEI(GI)
            ! For coefficient in the inverse mat of the jacobian. 
            A11= (EGI*KGI-FGI*HGI) /DETJ
            A21=-(DGI*KGI-FGI*GGI) /DETJ
            A31= (DGI*HGI-EGI*GGI) /DETJ
            !
            A12=-(BGI*KGI-CGI*HGI) /DETJ
            A22= (AGI*KGI-CGI*GGI) /DETJ
            A32=-(AGI*HGI-BGI*GGI) /DETJ
            !
            A13= (BGI*FGI-CGI*EGI) /DETJ
            A23=-(AGI*FGI-CGI*DGI) /DETJ
            A33= (AGI*EGI-BGI*DGI) /DETJ
            do  L=1,NLOC! Was loop 373
               NX(L,GI)= A11*NLX(L,GI)+A12*NLY(L,GI)+A13*NLZ(L,GI)
               NY(L,GI)= A21*NLX(L,GI)+A22*NLY(L,GI)+A23*NLZ(L,GI)
               NZ(L,GI)= A31*NLX(L,GI)+A32*NLY(L,GI)+A33*NLZ(L,GI)
            end do ! Was loop 373
            !
         end do ! Was loop 331
         ! IF(D3) THEN...
      ELSE IF(.NOT.D1) THEN
         TWOPIE=1.0 
         IF(DCYL) TWOPIE=2.*PIE
         do  GI=1,NGI! Was loop 1331
            !
            RGI=0.
            !
            AGI=0.
            BGI=0.
            CGI=0.
            DGI=0.
            !
            do  L=1,NLOC! Was loop 179
               IGLX=XONDGL((ELE-1)*NLOC+L)

               AGI=AGI + NLX(L,GI)*X(IGLX) 
               BGI=BGI + NLX(L,GI)*Y(IGLX) 
               CGI=CGI + NLY(L,GI)*X(IGLX) 
               DGI=DGI + NLY(L,GI)*Y(IGLX) 
               !
               RGI=RGI+N(L,GI)*Y(IGLX)
            end do ! Was loop 179
            !
            IF(.NOT.DCYL) RGI=1.0
            !
            DETJ= AGI*DGI-BGI*CGI 
            RA(GI)=RGI
            DETWEI(GI)=TWOPIE*RGI*DETJ*WEIGHT(GI)
            VOLUME=VOLUME+DETWEI(GI)
            !
            do L=1,NLOC
               NX(L,GI)=(DGI*NLX(L,GI)-BGI*NLY(L,GI))/DETJ
               NY(L,GI)=(-CGI*NLX(L,GI)+AGI*NLY(L,GI))/DETJ
               NZ(L,GI)=0.0
            END DO
            !
         end do ! Was loop 1331
         ! ENDOF IF(D3) THEN ELSE...
      ELSE 
         ! For 1D...
         do  GI = 1, NGI
            !
            AGI = 0.
            !
            do  L = 1, NLOC
               IGLX = XONDGL(( ELE - 1 ) * NLOC + L )
               AGI = AGI + NLX( L, GI ) * X( IGLX ) 
            end do
            !
            DETJ = AGI 
            DETWEI( GI ) = DETJ * WEIGHT( GI )
            VOLUME = VOLUME + DETWEI( GI )
            !
            do L = 1, NLOC
               NX( L, GI ) = NLX( L, GI ) / DETJ
               NY( L, GI ) = 0.0
               NZ( L, GI ) = 0.0
            END DO
            !
         end do
         ! ENDOF IF(D3) THEN ELSE...
      ENDIF
      !
      RETURN
    END SUBROUTINE DETNLXR

    real function lagran( diff, lx, inod, ndnod, nodpos )
      implicit none
      ! This return the Lagrange poly assocaited with node INOD at point LX
      ! If DIFF then send back the value of this poly differentiated. 
      logical :: diff
      real :: lx
      integer :: inod, ndnod
      real, dimension( 0 : ndnod ) :: nodpos
      ! Local variables
      integer :: n, k, i, j
      real :: denomi, over, over1

      n = ndnod - 1
      k = inod - 1

      denomi = 1.
      do i = 0, k - 1
         denomi = denomi * ( nodpos( k ) - nodpos( i ) )
      end do
      do i = k + 1, n
         denomi = denomi * ( nodpos( k ) - nodpos( i ) )
      end do

      Conditional_DIFF1: if( .not. diff ) then
         over = 1. 
         do i = 0, k - 1
            over = over * ( lx - nodpos( i ) )
         end do
         do i = k + 1, n
            over = over * ( lx - nodpos( i ) )
         end do
         lagran = over / denomi

      else ! Conditional_DIFF1
         over = 0.
         do j = 0, n
            if( j /= k ) then
               over1 = 1.
               do i = 0, k - 1
                  if( j /= i ) over1 = over1 * ( lx - nodpos( i ) )
               end do
               do i = k + 1, n
                  if( j /= i ) over1 = over1 * ( lx - nodpos( i ) )
               end do
               over = over + over1
            end if
         end do
         lagran = over / denomi

      end if Conditional_DIFF1

      return
    end function lagran



    REAL FUNCTION RGPTWE( IG, ND, WEIGHT )

      IMPLICIT NONE

      ! If WEIGHT is TRUE in function RGPTWE then return the Gauss-pt weight else 
      ! return the Gauss-pt.  There are ND Gauss points -- we are looking for either 
      ! the weight or the x-coord of the IG'th Gauss point. 

      INTEGER :: IG, ND
      LOGICAL :: WEIGHT

      Loop_Weight: IF( WEIGHT ) THEN ! Gauss points weights

         SELECT CASE( ND )

         CASE( 1 ) ; RGPTWE = 2.0

         CASE( 2 ) ; RGPTWE = 1.0

         CASE( 3 ) 
            SELECT CASE( IG )
            CASE( 1, 3 ) ; RGPTWE = 0.555555555555556
            CASE( 2 )    ; RGPTWE = 0.888888888888889
            END SELECT

         CASE( 4 )
            SELECT CASE( IG )
            CASE( 1, 4 ) ; RGPTWE = 0.347854845137454
            CASE( 2, 3 ) ; RGPTWE = 0.652145154862546
            END SELECT

         CASE( 5 )
            SELECT CASE( IG )
            CASE( 1, 5 ) ; RGPTWE = 0.236926885056189
            CASE( 2, 4 ) ; RGPTWE = 0.478628670499366
            CASE( 3 )    ; RGPTWE = 0.568888888888889
            END SELECT

         CASE( 6 )
            SELECT CASE( IG )
            CASE( 1, 6 ) ; RGPTWE = 0.171324492379170
            CASE( 2, 5 ) ; RGPTWE = 0.360761573048139
            CASE( 3, 4 ) ; RGPTWE = 0.467913934572691
            END SELECT

         CASE( 7 )
            SELECT CASE( IG )
            CASE( 1, 7 ) ; RGPTWE = 0.129484966168870
            CASE( 2, 6 ) ; RGPTWE = 0.279705391489277
            CASE( 3, 5 ) ; RGPTWE = 0.381830050505119
            CASE( 4 )    ; RGPTWE = 0.417959183673469
            END SELECT

         CASE( 8 )
            SELECT CASE( IG )
            CASE( 1, 8 ) ; RGPTWE = 0.101228536290376
            CASE( 2, 7 ) ; RGPTWE = 0.222381034453374
            CASE( 3, 6 ) ; RGPTWE = 0.313706645877877
            CASE( 4, 5 ) ; RGPTWE = 0.362683783378362
            END SELECT

         CASE( 9 )
            SELECT CASE( IG )
            CASE( 1, 9 ) ; RGPTWE = 0.081274388361574
            CASE( 2, 8 ) ; RGPTWE = 0.180648160694857
            CASE( 3, 7 ) ; RGPTWE = 0.260610696402935
            CASE( 4, 6 ) ; RGPTWE = 0.312347077040003
            CASE( 5 )    ; RGPTWE = 0.330239355001260
            END SELECT

         CASE( 10 )
            SELECT CASE( IG )
            CASE( 1, 10 ) ; RGPTWE = 0.066671344308688
            CASE( 2, 9 )  ; RGPTWE = 0.149451349150581
            CASE( 3, 8 )  ; RGPTWE = 0.219086362515982
            CASE( 4, 7 )  ; RGPTWE = 0.269266719309996
            CASE( 5, 6 )  ; RGPTWE = 0.295524224714753
            END SELECT

         CASE DEFAULT; FLExit( "In Lagrot subrt - wrong integer for NDGI " )

         END SELECT

      ELSE

         SELECT CASE( ND ) ! Gauss points

         CASE( 1 ) ; RGPTWE = 0.0

         CASE( 2 ) ; RGPTWE = 0.577350269189626

         CASE( 3 )
            SELECT CASE( IG )
            CASE( 1, 3 ) ; RGPTWE = 0.774596669241483
            CASE( 2 )    ; RGPTWE = 0.0
            END SELECT

         CASE( 4 )
            SELECT CASE( IG )
            CASE( 1, 4 ) ; RGPTWE = 0.861136311594953
            CASE( 2, 3 ) ; RGPTWE = 0.339981043584856
            END SELECT

         CASE( 5 )
            SELECT CASE( IG )
            CASE( 1, 5 ) ; RGPTWE = 0.906179845938664
            CASE( 2, 4 ) ; RGPTWE = 0.538469310105683
            CASE( 3 )    ; RGPTWE = 0.0
            END SELECT

         CASE( 6 )
            SELECT CASE( IG )
            CASE( 1, 6 ) ; RGPTWE = 0.932469514203152
            CASE( 2, 5 ) ; RGPTWE = 0.661209386466265
            CASE( 3, 4 ) ; RGPTWE = 0.238619186083197
            END SELECT

         CASE( 7 )
            SELECT CASE( IG )
            CASE( 1, 7 ) ; RGPTWE = 0.949107912342759
            CASE( 2, 6 ) ; RGPTWE = 0.741531185599394
            CASE( 3, 5 ) ; RGPTWE = 0.405845151377397
            CASE( 4 )    ; RGPTWE = 0.0
            END SELECT

         CASE( 8 )
            SELECT CASE( IG )
            CASE( 1, 8 ) ; RGPTWE = 0.960289856497536
            CASE( 2, 7 ) ; RGPTWE = 0.796666477413627
            CASE( 3, 6 ) ; RGPTWE = 0.525532409916329
            CASE( 4, 5 ) ; RGPTWE = 0.183434642495650
            END SELECT

         CASE( 9 )
            SELECT CASE( IG )
            CASE( 1, 9 ) ; RGPTWE = 0.968160239507626
            CASE( 2, 8 ) ; RGPTWE = 0.836031107326636
            CASE( 3, 7 ) ; RGPTWE = 0.613371432700590
            CASE( 4, 6 ) ; RGPTWE = 0.324253423403809
            CASE( 5)     ; RGPTWE = 0.0
            END SELECT

         CASE( 10 )
            SELECT CASE( IG )
            CASE( 1, 10 ) ; RGPTWE = 0.973906528517172
            CASE( 2, 9 )  ; RGPTWE = 0.865063366688985
            CASE( 3, 8 )  ; RGPTWE = 0.679409568299024
            CASE( 4, 7 )  ; RGPTWE = 0.433395394129247
            CASE( 5, 6 )  ; RGPTWE = 0.148874338981631
            END SELECT

         CASE DEFAULT; FLExit( "In Lagrot subrt - wrong integer for NDGI " )

         END SELECT

         IF( IG <= INT( ( ND / 2 ) + 0.1 )) RGPTWE = -RGPTWE

      END IF Loop_Weight

      RETURN

    END FUNCTION RGPTWE


  end module shape_functions_Linear_Quadratic

  
  module shape_functions_NDim

    use fldebug
    use shape_functions_Linear_Quadratic

  contains

    subroutine quad_1d_shape( cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, n, nlx, un, unlx )
      ! For quadatic elements. Shape functions associated with volume integration 
      ! using both CV basis functions CVN as well as FEM basis functions N (and 
      !its derivatives NLX, NLY, NLZ)
      implicit none
      integer, intent( in ) :: cv_ngi, cv_nloc, u_nloc
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: cvn
      real, dimension( cv_ngi ), intent( inout ) :: cvweigh
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: n, nlx
      real, dimension( u_nloc, cv_ngi ), intent( inout ) :: un, unlx
      ! Local Variables
      integer, parameter :: three = 3
      real, dimension( : ), allocatable :: lx, wei, rdummy, xi_min, xi_max, &
           cv_nodpos, u_nodpos
      integer :: gpoi, icv, p, nquad, iloc
      logical :: getndp, diff, ndiff
      real :: vol_cv, lx_tran, lxgp

      ewrite( 3, * ) 'In QUAD_1D_SHAPE'

      nquad = cv_ngi / cv_nloc

      ! Allocating memory
      allocate( lx( nquad ) )
      allocate( wei( nquad ) )
      allocate( rdummy( max( cv_ngi, u_nloc, cv_nloc ) ) )
      allocate( xi_min( nquad ) )
      allocate( xi_max( nquad ) )
      allocate( cv_nodpos( cv_nloc ) )
      allocate( u_nodpos( u_nloc ) )

      if( ( cv_ngi /= 18 ) .and. ( cv_ngi /= 15 ) .and. ( cv_ngi /= 12 ) .and. &
           ( cv_ngi /= 9 ) .and. ( cv_ngi /= 6 ) .and. ( cv_ngi /= 1 ) ) &
           FLAbort("Incorrect number of quadrature points")

      Case_CV_NLOC: Select Case( cv_nloc )
      case( 4 );
         ! Node 1
         xi_min( 1 ) = -1.
         xi_max( 1 ) = -1. + 1. / 3.
         ! Node 2
         xi_min( 2 ) = -1. + 1. / 3.
         xi_max( 2 ) =  0.
         ! Node 3
         xi_min( 3 ) =  0.
         xi_max( 3 ) =  1. - 1. / 3.
         ! Node 4
         xi_min( 4 ) =  1. - 1. / 3.
         xi_max( 4 ) =  1.
      case( 3 );
         if( .true. ) then
            ! Node 1
            xi_min( 1 ) = -1.
            xi_max( 1 ) = -0.5
            ! Node 2
            xi_min( 2 ) = -0.5
            xi_max( 2 ) =  0.5
            ! Node 3
            xi_min( 3 ) =  0.5
            xi_max( 3 ) =  1.
         else
            ! Node 1
            xi_min( 1 ) = -1.
            xi_max( 1 ) = -1. + 2. / 3.
            ! Node 2
            xi_min( 2 ) = -1. + 2. / 3. 
            xi_max( 2 ) =  1. - 2. / 3. 
            ! Node 3
            xi_min( 3 ) =  1. - 2. / 3.  
            xi_max( 3 ) =  1. 
         end if
      case( 2 );
         ! Node 1
         xi_min( 1 ) = -1.
         xi_max( 1 ) =  0.
         ! Node 2
         xi_min( 2 ) =  0.
         xi_max( 2 ) =  1.
      case( 1 )
         ! Node 1
         xi_min( 1 ) = -1.
         xi_max( 1 ) =  1.
      end Select Case_CV_NLOC

      diff = .true.
      ndiff = .false.
      ! Find the roots of the quadrature points and nodes
      ! also get the weights.
      getndp = .true.

      !     Compute standard Gauss quadrature: weights and points
      call lagrot( rdummy, cv_nodpos, cv_nloc, getndp )
      call lagrot( rdummy, u_nodpos, u_nloc, getndp )

      getndp = .false.

      !     Compute standard Gauss quadrature: weights and points
      call lagrot( wei, lx, nquad, getndp )

      ewrite(3,*)'cv_nloc, u_nloc, nquad, cv_ngi: ', cv_nloc, u_nloc, nquad, cv_ngi
      ewrite(3,*)'wei: ', wei
      ewrite(3,*)'lx: ', lx
      ewrite(3,*)'cv_nodpos: ', cv_nodpos
      ewrite(3,*)'u_nodpos: ', u_nodpos

      gpoi = 0
      cvn = 0.
      Loop_ICV1: do icv = 1, cv_nloc
         vol_cv = xi_max( icv ) - xi_min( icv )
         Loop_P1: do p = 1, nquad
            gpoi = gpoi + 1
            cvn( icv, gpoi ) = 1. ! Mapping to a new local coordinate system
            lx_tran = 0.5 * ( xi_max( icv ) + xi_min( icv ) ) + &
                 0.5 * ( xi_max( icv ) + xi_min( icv ) ) * lx( p )
            lxgp = lx_tran
            cvweigh( gpoi ) = wei( p ) * vol_cv / 2.
            Loop_ILOC1: do iloc = 1, cv_nloc
               n( iloc, gpoi ) = lagran( ndiff, lxgp, iloc, cv_nloc, cv_nodpos )
               nlx( iloc, gpoi ) = lagran( diff, lxgp, iloc, cv_nloc, cv_nodpos )
            end do Loop_ILOC1
            Loop_ILOC2: do iloc = 1, u_nloc
               un( iloc, gpoi ) = lagran( ndiff, lxgp, iloc, u_nloc, u_nodpos )
               unlx( iloc, gpoi ) = lagran( diff, lxgp, iloc, u_nloc, u_nodpos )
            end do Loop_ILOC2
         end do Loop_P1
      end do Loop_ICV1

      deallocate( lx )
      deallocate( wei )
      deallocate( rdummy )
      deallocate( xi_min )
      deallocate( xi_max )
      deallocate( cv_nodpos )
      deallocate( u_nodpos )

      ewrite(3,*) 'Leaving QUAD_1D_SHAPE'

      return
    end subroutine quad_1d_shape


    subroutine quad_nd_shape( ndim, cv_ele_type, cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, &
         n, nlx, nly, nlz, &
         un, unlx, unly, unlz )
      ! For quadratic elements: Shape functions associated with volume integration 
      ! using both CV (CVN) and FEM (N and its derivatives NLX/Y/Z) basis functions.

      implicit none
      integer, intent( in ) :: ndim, cv_ele_type, cv_ngi, cv_nloc, u_nloc
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: cvn
      real, dimension( cv_ngi ), intent( inout ) :: cvweigh 
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: n, nlx, nly, nlz
      real, dimension( u_nloc, cv_ngi ), intent( inout ) :: un, unlx, unly, unlz

      ! Local variables
      real, dimension( :, : ), allocatable :: cvn_dum, cvn_1d_dum, n_1d, &
           nlx_1d, un_1d, unlx_1d, cvn_1d
      real, dimension( : ), allocatable :: cvweigh_dum, cvweigh_1d_dum, cvweigh_1d
      integer :: cv_ngi_1d, cv_nloc_1d, u_nloc_1d

      ewrite(3,*) 'In MD_Shapes subrt: quad_nd_shape'

      Conditional_Dimensionality: Select Case( ndim )
      case( 1 )
         cv_ngi_1d = cv_ngi
         cv_nloc_1d = cv_nloc
         u_nloc_1d = u_nloc

      case( 2 )
         cv_ngi_1d = int( sqrt( 1.e-3 + cv_ngi ))
         cv_nloc_1d = int( sqrt( 1.e-3 + cv_nloc ))
         u_nloc_1d = int( sqrt( 1.e-3 + u_nloc ))

      case( 3 )
         cv_ngi_1d = int( ( 1.e-3 + cv_ngi ) ** 1./3. )
         cv_nloc_1d = int( ( 1.e-3 + cv_nloc ) ** 1./3. )
         u_nloc_1d = int( ( 1.e-3 + u_nloc ) ** 1./3. )

      case default; FLExit( " Invalid integer for NDIM " )

      end Select Conditional_Dimensionality

      ! Allocating memory
      allocate( cvweigh_1d_dum( cv_ngi_1d ) )
      allocate( n_1d( cv_nloc_1d, cv_ngi_1d ) )
      allocate( nlx_1d( cv_nloc_1d, cv_ngi_1d ) )
      allocate( un_1d( u_nloc_1d, cv_ngi_1d ) )
      allocate( unlx_1d( u_nloc_1d, cv_ngi_1d ) )
      allocate( cvn_1d( cv_nloc_1d, cv_ngi_1d ) )

      ! Computing Shape Functions for 1D:
      call quad_1d_shape( cv_ngi_1d, cv_nloc_1d, u_nloc_1d, cvn_1d, cvweigh_1d, &
           n_1d, nlx_1d, un_1d, unlx_1d )

      ! Computing Shape Functions for N:
      call quad_nd_shape_N( cv_ele_type, ndim, cv_ngi, cv_nloc, cvn, cvweigh, &
           n, nlx, nly, nlz, &
           cv_ngi_1d, cv_nloc_1d, cvn_1d, cvweigh_1d, n_1d, nlx_1d )

      ! Allocating memory for UN calculation:
      allocate( cvn_dum( cv_nloc, cv_ngi ) )
      allocate( cvweigh_dum( cv_ngi ) )
      allocate( cvn_1d_dum( cv_nloc_1d, cv_ngi_1d ) )
      allocate( cvweigh_1d_dum( cv_ngi_1d ) )

      ! Computing Shape Functions for UN:
      call quad_nd_shape_N( cv_ele_type, ndim, cv_ngi, u_nloc, cvn_dum, cvweigh_dum, &
           un, unlx, unly, unlz, &
           cv_ngi_1d, u_nloc_1d, cvn_1d_dum, cvweigh_1d_dum, un_1d, unlx_1d )

      deallocate( cvn_dum )
      deallocate( cvn_1d_dum )
      deallocate( n_1d )
      deallocate( nlx_1d )
      deallocate( un_1d )
      deallocate( unlx_1d )
      deallocate( cvweigh_dum )
      deallocate( cvweigh_1d_dum )
      deallocate( cvweigh_1d )
      deallocate( cvn_1d ) 

      return
    end subroutine quad_nd_shape


    subroutine quad_nd_shape_N( cv_ele_type, ndim, cv_ngi, cv_nloc, cvn, cvweigh, &
         n, nlx, nly, nlz, &
         cv_ngi_1d, cv_nloc_1d, cvn_1d, cvweigh_1d, n_1d, nlx_1d )
      ! For quadatic elements -- shape functions associated with volume 
      ! integration using both CV basis functions CVN as well as FEM basis
      ! functions N (and its derivatives NLX, NLY, NLZ)
      implicit none
      integer, intent( in ) :: cv_ele_type, ndim, cv_ngi, cv_nloc
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: cvn
      real, dimension( cv_ngi ), intent( inout ) :: cvweigh
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: n, nlx, nly, nlz
      integer, intent( in ) :: cv_ngi_1d, cv_nloc_1d
      real, dimension( cv_nloc_1d, cv_ngi_1d ), intent( in ) :: cvn_1d
      real, dimension( cv_ngi_1d ), intent( in ) :: cvweigh_1d
      real, dimension( cv_nloc_1d, cv_ngi_1d ), intent( in ) :: n_1d, nlx_1d
      ! Local variables
      integer :: cv_iloc_1d, cv_jloc_1d, cv_kloc_1d, cv_iloc, cv_jloc, cv_kloc, cv_igi, &
           cv_igi_1d, cv_jgi_1d, cv_kgi_1d

      Conditional_Dimensionality: Select Case( ndim )
      case( 1 )
         cvn = cvn_1d
         n = n_1d
         nlx = nlx_1d
         nly = 0.
         nlz = 0.
         cvweigh = cvweigh_1d

      case( 2 )
         cv_iloc = 0
         nlz = 0.
         Loop_CVILOC_1D_2D: do cv_iloc_1d = 1, cv_nloc_1d
            Loop_CVJLOC_1D_2D: do cv_jloc_1d = 1, cv_nloc_1d
               cv_iloc = cv_iloc + 1
               cv_igi = 0
               do cv_igi_1d = 1, cv_nloc_1d
                  do cv_jgi_1d = 1, cv_nloc_1d
                     cv_igi = cv_igi + 1
                     cvn( cv_iloc, cv_igi ) =  &
                          cvn_1d( cv_iloc_1d, cv_igi_1d ) * cvn_1d( cv_jloc_1d, cv_jgi_1d )
                     n( cv_iloc, cv_igi ) =  &
                          n_1d( cv_iloc_1d, cv_igi_1d ) * n_1d( cv_jloc_1d, cv_jgi_1d )
                     nlx( cv_iloc, cv_igi ) =  &
                          nlx_1d( cv_iloc_1d, cv_igi_1d ) * n_1d( cv_jloc_1d, cv_jgi_1d )
                     nly( cv_iloc, cv_igi ) =  &
                          n_1d( cv_iloc_1d, cv_igi_1d ) * nlx_1d( cv_jloc_1d, cv_jgi_1d )
                     cvweigh( cv_igi ) = cvweigh_1d( cv_igi_1d ) * cvweigh_1d( cv_jgi_1d )
                  end do
               end do
            end do Loop_CVJLOC_1D_2D
         end do Loop_CVILOC_1D_2D

      case( 3 )
         cv_iloc = 0
         Loop_CVILOC_1D_3D: do cv_iloc_1d = 1, cv_nloc_1d
            Loop_CVJLOC_1D_3D: do cv_jloc_1d = 1, cv_nloc_1d
               Loop_CVKLOC_1D_3D: do cv_kloc_1d = 1, cv_nloc_1d
                  cv_iloc = cv_iloc + 1
                  cv_igi = 0
                  do cv_igi_1d = 1, cv_nloc_1d
                     do cv_jgi_1d = 1, cv_nloc_1d
                        do cv_kgi_1d = 1, cv_nloc_1d
                           cv_igi = cv_igi + 1
                           cvn( cv_iloc, cv_igi ) =  cvn_1d( cv_iloc_1d, cv_igi_1d ) * &
                                cvn_1d( cv_jloc_1d, cv_jgi_1d ) * cvn_1d( cv_kloc_1d, cv_kgi_1d )
                           n( cv_iloc, cv_igi ) =  n_1d( cv_iloc_1d, cv_igi_1d ) * &
                                n_1d( cv_jloc_1d, cv_jgi_1d ) * n_1d( cv_kloc_1d, cv_kgi_1d )
                           nlx( cv_iloc, cv_igi ) =  nlx_1d( cv_iloc_1d, cv_igi_1d ) * &
                                n_1d( cv_jloc_1d, cv_jgi_1d ) * n_1d( cv_kloc_1d, cv_kgi_1d )
                           nly( cv_iloc, cv_igi ) =  n_1d( cv_iloc_1d, cv_igi_1d ) * &
                                nlx_1d( cv_jloc_1d, cv_jgi_1d ) * n_1d( cv_kloc_1d, cv_kgi_1d )
                           nlz( cv_iloc, cv_igi ) =  n_1d( cv_iloc_1d, cv_igi_1d ) * &
                                n_1d( cv_jloc_1d, cv_jgi_1d ) * nlx_1d( cv_kloc_1d, cv_kgi_1d )
                           cvweigh( cv_igi ) = cvweigh_1d( cv_igi_1d ) * & 
                                cvweigh_1d( cv_jgi_1d ) * cvweigh_1d( cv_kgi_1d )
                        end do
                     end do
                  end do
               end do Loop_CVKLOC_1D_3D
            end do Loop_CVJLOC_1D_3D
         end do Loop_CVILOC_1D_3D

      case default; FLExit( " Invalid integer for NDIM " )

      end Select Conditional_Dimensionality

      return
    end subroutine quad_nd_shape_N


    subroutine vol_cv_tri_shape( cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, &
         n, nlx, nly, nlz, un, unlx, unly, unlz )
      ! Compute shape functions N, UN etc for linear trianles. Shape functions 
      ! associated with volume integration using both CV basis functions CVN, as 
      ! well as FEM basis functions N (and its derivatives NLX, NLY, NLZ). Also 
      ! for velocity basis functions UN, UNLX, UNLY, UNLZ.
      implicit none
      integer, intent( in ) :: cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: cvn
      real, dimension( cv_ngi ), intent( inout ) :: cvweigh
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: n, nlx, nly, nlz
      real, dimension( u_nloc, cv_ngi ), intent( inout ) :: un, unlx, unly, unlz
      ! Local variables
      integer, parameter :: x_nonods = 7, totele = 3, quad_cv_nloc = 4, quad_u_loc_dummy = 1
      real, parameter :: h_scale = 2.70240031 !scaling factor to give a unity area of the local triangle
      integer, dimension( : ), allocatable :: x_ndgln, fem_nod
      real, dimension( : ), allocatable :: lx, ly, lz, x, y, z
      integer :: ele, cv_iloc, quad_cv_ngi
      logical, parameter :: d1 = .false., dcyl = .false., d3 = .false.

      ! Allocating memory
      allocate( lx( totele ) )
      allocate( ly( totele ) )
      allocate( lz( totele ) )
      allocate( x( x_nonods ))
      allocate( y( x_nonods ))
      allocate( z( x_nonods ))
      allocate( fem_nod( x_nonods ) )

      fem_nod = 0

      ! Setting-up unity area triangle
      lx( 1 ) = 0.
      ly( 1 ) = 0.

      lx( 2 ) = 1. * h_scale
      ly( 2 ) = 0.

      lx( 3 ) = .5 * h_scale
      ly( 3 ) = sqrt( 3./4. ) * h_scale

      ! Remmaping
      x( 1 ) = lx( 1 )
      y( 1 ) = ly( 1 )
      fem_nod( 1 ) = 1


      x( 2 ) = 0.5 * ( lx( 1 ) + lx( 2 ) )
      y( 2 ) = 0.5 * ( ly( 1 ) + ly( 2 ) )


      x( 3 ) = lx( 2 )
      y( 3 ) = ly( 2 )
      fem_nod( 3 ) = 2


      x( 4 ) = 0.5 * ( lx( 3 ) + lx( 1 ) )
      y( 4 ) = 0.5 * ( ly( 3 ) + ly( 1 ) )


      x( 5 ) = 1./3. * ( lx( 1 ) + lx( 2 ) + lx( 3 ) )
      y( 5 ) = 1./3. * ( ly( 1 ) + ly( 2 ) + ly( 3 ) )


      x( 6 ) = 0.5 * ( lx( 2 ) + lx( 3 ) )
      y( 6 ) = 0.5 * ( ly( 2 ) + ly( 3 ) )


      x( 7 ) = lx( 3 )
      y( 7 ) = ly( 3 )
      fem_nod( 7 ) = 3

      z = 0.

      ! Computing sudo elements
      allocate( x_ndgln( totele * quad_cv_nloc ) )
      ele = 1
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 1
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 2
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 4
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 5

      ele = 2
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 2
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 3
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 5
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 6

      ele = 3
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 4
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 5
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 7
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 6


      ! Compute the shape functions using these quadrilaterals/hexs:
      ! For pressure:
      call shape_tri_tet( cv_ele_type, ndim, totele, cv_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, &
           n, nlx, nly, nlz, cvweigh )

      ! For velocities:
      call shape_tri_tet( cv_ele_type, ndim, totele, u_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, &
           un, unlx, unly, unlz, cvweigh )

      ! Compute CVN (i.e., CV basis function):
      quad_cv_ngi = cv_ngi / totele
      call calc_cvn_tri_tet( cv_ele_type, totele, cv_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, fem_nod, cvn )

      deallocate( lx )
      deallocate( ly )
      deallocate( lz )
      deallocate( x )
      deallocate( y ) 
      deallocate( z ) 
      deallocate( x_ndgln )
      deallocate( fem_nod )

      return
    end subroutine vol_cv_tri_shape



    subroutine vol_cv_qtri_shape( cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, &
         n, nlx, nly, nlz, un, unlx, unly, unlz )
      ! Compute shape functions N, UN etc for quadratic trianles. Shape functions 
      ! associated with volume integration using both CV basis functions CVN, as 
      ! well as FEM basis functions N (and its derivatives NLX, NLY, NLZ). Also 
      ! for velocity basis functions UN, UNLX, UNLY, UNLZ.
      implicit none
      integer, intent( in ) :: cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: cvn
      real, dimension( cv_ngi ), intent( inout ) :: cvweigh
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: n, nlx, nly, nlz
      real, dimension( u_nloc, cv_ngi ), intent( inout ) :: un, unlx, unly, unlz
      ! Local variables
      integer, parameter :: x_nonods = 16, totele = 9, quad_cv_nloc = 4, quad_u_loc_dummy = 1
      real, parameter :: h_scale = 2.70240031 !scaling factor to give a unity area of the local triangle
      integer, dimension( : ), allocatable :: x_ndgln, fem_nod
      real, dimension( : ), allocatable :: lx, ly, lz, x, y, z
      integer :: ele, cv_iloc, quad_cv_ngi
      logical, parameter :: d1 = .false., dcyl = .false., d3 = .true.

      ! Allocating memory
      allocate( lx( totele ) )
      allocate( ly( totele ) )
      allocate( lz( totele ) )
      allocate( x( x_nonods ))
      allocate( y( x_nonods ))
      allocate( z( x_nonods ))
      allocate( fem_nod( x_nonods ) )

      fem_nod = 0

      ! Setting-up unity area triangle
      lx( 1 ) = 0.
      ly( 1 ) = 0.

      lx( 2 ) = 1. * h_scale
      ly( 2 ) = 0.

      lx( 3 ) = .5 * h_scale
      ly( 3 ) = sqrt( 3./4. ) * h_scale

      ! Remmaping
      ! Node point 1
      x( 1 ) = lx( 1 )
      y( 1 ) = ly( 1 )
      fem_nod( 1 ) = 1
      ! Node point 2
      x( 2 ) = 0.25 * lx( 2 ) + 3. / 4. * lx( 1 ) 
      y( 2 ) = 0.25 * ly( 2 ) + 3. / 4. * ly( 2 ) 
      ! Node point 3
      x( 3 ) = 0.5 * ( lx( 1 ) + lx( 2 ) )
      y( 3 ) = 0.5 * ( ly( 1 ) + ly( 2 ) )
      fem_nod( 3 ) = 2
      ! Node point 4
      x( 4 ) = 0.25 * lx( 1 ) + 3. / 4. * lx( 2 ) 
      y( 4 ) = 0.25 * lx( 1 ) + 3. / 4. * lx( 2 ) 
      ! Node point 5
      x( 5 ) = lx( 2 )
      y( 5 ) = ly( 2 )
      fem_nod( 5 ) = 3
      ! Node point 11
      x( 11 ) = 0.5 * ( lx( 1 ) + lx( 3 ) )
      y( 11 ) = 0.5 * ( ly( 1 ) + ly( 3 ) )
      fem_nod( 11 ) = 4
      ! Node point 12
      x( 12 ) = 0.5 * ( lx( 2 ) + lx( 3 ) )
      y( 12 ) = 0.5 * ( ly( 2 ) + ly( 3 ) )
      fem_nod( 12 ) = 5
      ! Node point 16
      x( 16 ) = lx( 3 )
      y( 16 ) = ly( 3 )
      fem_nod( 16 ) = 6
      ! Node point 6
      x( 6 ) = 0.5 * ( x( 1 ) + x( 11 ) )
      y( 6 ) = 0.5 * ( y( 1 ) + y( 11 ) )
      ! Node point 7
      x( 7 ) = 1. / 3. * ( x( 1 ) + x( 3 ) + x( 11 ) )
      y( 7 ) = 1. / 3. * ( y( 1 ) + y( 3 ) + y( 11 ) )
      ! Node point 8
      x( 8 ) = 1. / 3. * ( x( 3 ) + x( 11 ) + x( 12 ) )
      y( 8 ) = 1. / 3. * ( y( 3 ) + y( 11 ) + y( 12 ) )
      ! Node point 9
      x( 9 ) = 1. / 3. * ( x( 3 ) + x( 5 ) + x( 12 ) )
      y( 9 ) = 1. / 3. * ( y( 3 ) + y( 5 ) + y( 12 ) )
      ! Node point 10
      x( 10 ) = 0.5 * ( x( 5 ) + x( 12 ) )
      y( 10 ) = 0.5 * ( y( 5 ) + y( 12 ) )
      ! Node point 13
      x( 13 ) = 1. / 3. * ( x( 11 ) + x( 12 ) + x( 16 ) )
      y( 13 ) = 1. / 3. * ( y( 11 ) + y( 12 ) + y( 16 ) )
      ! Node point 14
      x( 14 ) = 0.5 * ( x( 11 ) + x( 16 ) )
      y( 14 ) = 0.5 * ( y( 11 ) + y( 16 ) )
      ! Node point 15
      x( 15 ) = 0.5 * ( x( 12 ) + x( 16 ) )
      y( 15 ) = 0.5 * ( y( 12 ) + y( 16 ) )

      z = 0.

      ! Computing sudo elements
      allocate( x_ndgln( totele * quad_cv_nloc ) )
      ele = 1 
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 1 
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 2
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 6
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 7

      ele = 2
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 2
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 3
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 7
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 8

      ele = 3
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 3
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 4
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 8
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 9

      ele = 4
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 4
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 5
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 9
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 10

      ele = 5
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 6
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 7
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 11
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 8

      ele = 6
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 8
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 9
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 12
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 10

      ele = 7
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 11
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 8
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 14
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 13

      ele = 8
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 8
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 12
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 13
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 15

      ele = 9
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 14
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 13
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 16
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 15


      ! Compute the shape functions using these quadrilaterals/hexs:
      ! For pressure:
      call shape_tri_tet( cv_ele_type, ndim, totele, cv_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, &
           n, nlx, nly, nlz, cvweigh )

      ! For velocities:
      call shape_tri_tet( cv_ele_type, ndim, totele, u_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, &
           un, unlx, unly, unlz, cvweigh )

      ! Compute CVN (i.e., CV basis function):
      quad_cv_ngi = cv_ngi / totele
      call calc_cvn_tri_tet( cv_ele_type, totele, cv_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, fem_nod, cvn )


      deallocate( lx )
      deallocate( ly )
      deallocate( lz )
      deallocate( x )
      deallocate( y ) 
      deallocate( z ) 
      deallocate( x_ndgln )
      deallocate( fem_nod )

      return
    end subroutine vol_cv_qtri_shape


    subroutine vol_cv_tet_shape( cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, &
         n, nlx, nly, nlz, un, unlx, unly, unlz )
      ! Compute the shape functions N, UN, etc for linear tetrahedra. 
      ! The shape functions are associated with volume integration using both CV basis 
      ! functions CVN as well as FEM basis functions N (and its derivatives NLX, NLY, NLZ)
      ! also for velocity basis functions UN, UNLX, UNLY, UNLZ
      implicit none
      integer, intent( in ) :: cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: cvn
      real, dimension( cv_ngi ), intent( inout ) :: cvweigh
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: n, nlx, nly, nlz
      real, dimension( u_nloc, cv_ngi ), intent( inout ) :: un, unlx, unly, unlz
      ! Local variables
      integer, parameter :: x_nonods = 15, totele = 4, quad_cv_nloc = 8
      real, parameter :: h_scale = 2.04 !scaling factor to give a unity area of the local triangle
      integer, dimension( : ), allocatable :: x_ndgln, fem_nod
      real, dimension( : ), allocatable :: lx, ly, lz, x, y, z
      integer :: ele, cv_iloc, quad_cv_ngi

      ! Allocating memory
      allocate( lx( totele ) )
      allocate( ly( totele ) )
      allocate( lz( totele ) )
      allocate( x( x_nonods ))
      allocate( y( x_nonods ))
      allocate( z( x_nonods ))
      allocate( fem_nod( x_nonods ) )

      fem_nod = 0! FEM_NOD contains the local triangle/tet FEM node

      ! Setting-up unity volume tetrahedron
      lx( 1 ) = 0.
      ly( 1 ) = 0.
      lz( 1 ) = 0.

      lx( 2 ) = 1. * h_scale
      ly( 2 ) = 0.
      lz( 2 ) = 0.

      lx( 3 ) = .5 * h_scale
      ly( 3 ) = sqrt( 3. / 4. ) * h_scale
      lz( 3 ) = 0.

      lx( 4 ) = 0.5 * ( lx( 1 ) + lx( 2 ) ) * h_scale
      ly( 4 ) = 0.5 * ( ly( 1 ) + ly( 3 ) ) * h_scale 
      lz( 4 ) = sqrt( 2. / 3. ) * h_scale

      ! Remmaping
      z = 0.
      ! Node point 1
      x( 1 ) = lx( 1 )
      y( 1 ) = ly( 1 )
      fem_nod( 1 ) = 1
      ! Node point 2
      x( 2 ) = 0.5 * ( lx( 1 ) + lx( 2 ) )
      y( 2 ) = 0.5 * ( ly( 1 ) + ly( 2 ) )
      ! Node point 3
      x( 3 ) = lx( 2 )
      y( 3 ) = ly( 2 )
      fem_nod( 3 ) = 2
      ! Node point 4
      x( 4 ) = 0.5 * ( lx( 3 ) + lx( 1 ) )
      y( 4 ) = 0.5 * ( ly( 3 ) + ly( 1 ) )
      ! Node point 5
      x( 5 ) = 1./3. * ( lx( 1 ) + lx( 2 ) + lx( 3 ) )
      y( 5 ) = 1./3. * ( ly( 1 ) + ly( 2 ) + ly( 3 ) )
      ! Node point 6
      x( 6 ) = 0.5 * ( lx( 2 ) + lx( 3 ) )
      y( 6 ) = 0.5 * ( ly( 2 ) + ly( 3 ) )
      ! Node point 7
      x( 7 ) = lx( 3 )
      y( 7 ) = ly( 3 )
      fem_nod( 7 ) = 3
      ! Node point 15
      x( 15 ) = lx( 4 )
      y( 15 ) = lx( 4 )
      z( 15 ) = lz( 4 )
      fem_nod( 15 ) = 4
      ! Node point 8
      x( 8 ) = 1./3. * ( x( 1 ) + x( 3 ) + x( 15 ) )
      y( 8 ) = 1./3. * ( y( 1 ) + y( 3 ) + y( 15 ) )
      z( 8 ) = 1./3. * ( z( 1 ) + z( 3 ) + z( 15 ) )
      ! Node point 9
      x( 9 ) = 0.5 * ( x( 1 ) + x( 15 ) )
      y( 9 ) = 0.5 * ( y( 1 ) + y( 15 ) )
      z( 9 ) = 0.5 * ( z( 1 ) + z( 15 ) )
      ! Node point 10
      x( 10 ) = 0.5 * ( x( 3 ) + x( 15 ) )
      y( 10 ) = 0.5 * ( y( 3 ) + y( 15 ) )
      z( 10 ) = 0.5 * ( z( 3 ) + z( 15 ) )
      ! Node point 11
      x( 11 ) = 0.25 * ( lx( 1 ) + lx( 2 ) + lx( 3 ) + lx( 4 ) )
      y( 11 ) = 0.25 * ( ly( 1 ) + ly( 2 ) + ly( 3 ) + ly( 4 ) )
      z( 11 ) = 0.25 * ( lz( 1 ) + lz( 2 ) + lz( 3 ) + lz( 4 ) )
      ! Node point 12
      x( 12 ) = 1./3. * ( x( 1 ) + x( 7 ) + x( 15 ) )
      y( 12 ) = 1./3. * ( y( 1 ) + y( 7 ) + y( 15 ) )
      z( 12 ) = 1./3. * ( z( 1 ) + z( 7 ) + z( 15 ) )
      ! Node point 13
      x( 13 ) = 1./3. * ( x( 3 ) + x( 7 ) + x( 15 ) )
      y( 13 ) = 1./3. * ( y( 3 ) + y( 7 ) + y( 15 ) )
      z( 13 ) = 1./3. * ( z( 3 ) + z( 7 ) + z( 15 ) )
      ! Node point 14
      x( 14 ) = 0.5 * ( x( 7 ) + x( 15 ) )
      y( 14 ) = 0.5 * ( y( 7 ) + y( 15 ) )
      z( 14 ) = 0.5 * ( z( 7 ) + z( 15 ) )

      ! Computing sudo elements
      allocate( x_ndgln( totele * quad_cv_nloc ) )
      ele = 1
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 1
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 2
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 4
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 5

      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 5 ) = 9
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 6 ) = 8
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 7 ) = 11
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 8 ) = 12

      ele = 2
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 2
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 3
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 5
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 6

      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 5 ) = 8
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 6 ) = 10
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 7 ) = 12
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 8 ) = 13

      ele = 3
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 4
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 5
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 7
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 6

      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 5 ) = 11
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 6 ) = 12
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 7 ) = 14
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 8 ) = 13

      ele = 4
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 9
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 8
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 11
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 12

      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 5 ) = 15
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 6 ) = 10
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 7 ) = 14
      x_ndgln( ( ele - 1 ) * quad_cv_nloc + 8 ) = 13

      ! Compute the shape functions using these quadrilaterals/hexs:
      ! For pressure:
      call shape_tri_tet( cv_ele_type, ndim, totele, cv_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, &
           n, nlx, nly, nlz, cvweigh )

      ! For velocities:
      call shape_tri_tet( cv_ele_type, ndim, totele, u_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, &
           un, unlx, unly, unlz, cvweigh )

      ! Compute CVN (i.e., CV basis function):
      quad_cv_ngi = cv_ngi / totele
      call calc_cvn_tri_tet( cv_ele_type, totele, cv_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, fem_nod, cvn )

      deallocate( lx )
      deallocate( ly )
      deallocate( lz )
      deallocate( x )
      deallocate( y ) 
      deallocate( z ) 
      deallocate( x_ndgln )
      deallocate( fem_nod )

      return
    end subroutine vol_cv_tet_shape


    subroutine vol_cv_qtet_shape( cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, &
         n, nlx, nly, nlz, un, unlx, unly, unlz )
      ! Compute shape functions N, UN etc for quadratic tetrahedra. Shape functions 
      ! associated with volume integration using both CV basis functions CVN, as 
      ! well as FEM basis functions N (and its derivatives NLX, NLY, NLZ). Also 
      ! for velocity basis functions UN, UNLX, UNLY, UNLZ.
      use pascal_tetrahedra
      implicit none
      integer, intent( in ) :: cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: cvn
      real, dimension( cv_ngi ), intent( inout ) :: cvweigh
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: n, nlx, nly, nlz
      real, dimension( u_nloc, cv_ngi ), intent( inout ) :: un, unlx, unly, unlz
      ! Local variables
      integer, dimension( : ), allocatable :: x_ndgln, x_ndgln2, x_ndgln_big, tet_ndgln
      integer :: nonods, nod, nnods, ilayer, ilayer0, inod, x_nonods, &
           triangle_totele, iloc, ifem, nodeplustetnodes, quad_cv_ngi
      real, dimension( : ), allocatable :: lx, ly, lz, x, y, z
      integer, dimension( : ), allocatable :: fem_nod


      x_nonods = no_of_nodes_in_faces * no_faces * totele + totele
      triangle_totele = totele * tet_totele

      ! Allocating memory
      allocate( x( x_nonods ) )
      allocate( y( x_nonods ) )
      allocate( z( x_nonods ) )
      allocate( lx( x_nonods ) )
      allocate( ly( x_nonods ) )
      allocate( lz( x_nonods ) )
      allocate( x_ndgln( x_nonods ) )
      allocate( x_ndgln2( x_nonods ) )
      allocate( x_ndgln_big( x_nonods ) )
      allocate( fem_nod( x_nonods ) )

!!!
!!!   CALL FOR QUAD_TETS (or 8 LINEAR TETS)
!!!
      call numbering_tet_elements( x_nonods, triangle_totele, nodeplustetnodes, &
           x_ndgln, &
           x, y, z, lx, ly, lz, fem_nod, &
           x_ndgln_big, x_ndgln2 )  


      ! Compute the shape functions using these quadrilaterals/hexs:
      ! For pressure:
      call shape_tri_tet( cv_ele_type, ndim, totele, cv_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, &
           n, nlx, nly, nlz, cvweigh )

      ! For velocities:
      call shape_tri_tet( cv_ele_type, ndim, totele, u_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, &
           un, unlx, unly, unlz, cvweigh )

      ! Compute CVN (i.e., CV basis function):
      quad_cv_ngi = cv_ngi / totele
      call calc_cvn_tri_tet( cv_ele_type, totele, cv_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, fem_nod, cvn )

      deallocate( lx )
      deallocate( ly )
      deallocate( lz )
      deallocate( x )
      deallocate( y ) 
      deallocate( z ) 
      deallocate( x_ndgln )
      deallocate( x_ndgln2 )
      deallocate( x_ndgln_big )
      deallocate( fem_nod )

      return
    end subroutine vol_cv_qtet_shape



    subroutine shape_tri_tet( cv_ele_type, ndim, totele, cv_nloc, cv_ngi, x_nonods, &
         quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, &
         n, nlx, nly, nlz, cvweigh )
      use shape_functions
      implicit none
      integer, intent( in ) :: cv_ele_type, ndim, totele, cv_nloc, cv_ngi, &
           x_nonods, quad_cv_nloc
      integer, dimension( totele * quad_cv_nloc ), intent( in ) :: x_ndgln
      real, dimension( x_nonods ), intent( in ) :: x, y, z
      real, dimension( quad_cv_nloc ), intent( in ) :: lx, ly, lz ! corner nodes
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: n, nlx, nly, nlz
      real, dimension( cv_ngi ), intent( inout ) ::  cvweigh
      ! Local variables
      logical :: d1, dcyl, d3, lowqua
      integer :: ele, quad_cv_ngi, quad_cv_gi, cv_gi, nwicel, xnod, quad_cv_iloc, &
           quad_u_loc_dummy, mloc, dummy_sngi, dummy_snloc, dummy_smloc
      real :: xgi, ygi, zgi, volume
      real, dimension( : ), allocatable :: quad_l1, quad_l2, quad_l3, quad_l4, &
           quad_cvweight, detwei, ra, rdummy
      real, dimension( :, : ), allocatable :: quad_n, quad_nlx, quad_nly, quad_nlz, &
           quad_nx, quad_ny, quad_nz

      d3 = ( ndim == 3 )
      d1 = .false.
      dcyl = .false. 

      Conditional_Dimensionality1: if( d3 ) then
         quad_cv_ngi = 8
         nwicel = 3
         if( cv_nloc == 10 ) then ! Quadratic hexs
            quad_cv_ngi = 27
            nwicel = 4
         end if
      else
         quad_cv_ngi = 4
         nwicel = 1
         if( cv_nloc == 6 ) then ! Quadratic quads
            quad_cv_ngi = 9
            nwicel = 2
         end if
      end if Conditional_Dimensionality1

      ! Consistency check:
      if( cv_ngi /= totele * quad_cv_ngi ) then
         FLExit( "Wrong number for CV_NGI" )
      end if

      quad_u_loc_dummy = 1
      mloc = 1
      dummy_sngi = 1
      dummy_snloc = 1
      dummy_smloc = 1
      lowqua = .false.

      ! Allocating memory
      allocate( quad_l1( cv_ngi ) )
      allocate( quad_l2( cv_ngi ) )
      allocate( quad_l3( cv_ngi ) )
      allocate( quad_l4( cv_ngi ) )
      allocate( quad_cvweight( quad_cv_ngi ) )
      allocate( detwei( quad_cv_ngi ) )
      allocate( ra( quad_cv_ngi ) )
      allocate( quad_n( quad_cv_nloc, quad_cv_ngi ) )
      allocate( quad_nlx( quad_cv_nloc, quad_cv_ngi ) )
      allocate( quad_nly( quad_cv_nloc, quad_cv_ngi ) )
      allocate( quad_nlz( quad_cv_nloc, quad_cv_ngi ) )
      allocate( quad_nx( quad_cv_nloc, quad_cv_ngi ) )
      allocate( quad_ny( quad_cv_nloc, quad_cv_ngi ) )
      allocate( quad_nz( quad_cv_nloc, quad_cv_ngi ) )
      allocate( rdummy( 10000 ) )

      ! Now we need to compute QUAD_NLX/Y/Z - get the hex or quad
      ! shape functions quad_n etc.
      call shape_l_q_quad( lowqua, quad_cv_ngi, quad_cv_nloc, mloc, &
           dummy_sngi, dummy_snloc, dummy_smloc, rdummy, rdummy, rdummy, rdummy, &
           quad_cvweight, quad_n, quad_nlx, quad_nly, quad_nlz, &
           rdummy, rdummy, rdummy, rdummy, rdummy, rdummy, rdummy, &
           nwicel, d3 )   

      Loop_Elements: do ele = 1, totele ! Calculate DETWEI,RA,NX,NY,NZ for element ELE

         call detnlxr( ele, x, y, z, x_ndgln, totele, x_nonods, quad_cv_nloc, quad_cv_ngi, &
              quad_n, quad_nlx, quad_nly, quad_nlz, quad_cvweight, &
              detwei, ra, volume, d1, d3, dcyl, &       
              quad_nx, quad_ny, quad_nz )

         Loop_NGI: do quad_cv_gi = 1, quad_cv_ngi ! Determine the quadrature points and weights
            cv_gi = ( ele - 1 ) * quad_cv_ngi + quad_cv_gi
            cvweigh( cv_ngi ) = detwei( quad_cv_gi )
            xgi = 0.
            ygi = 0.
            zgi = 0.

            do quad_cv_iloc = 1, quad_cv_nloc
               xnod = x_ndgln( ( ele - 1 ) * quad_cv_nloc + quad_cv_iloc )
               xgi =  xgi + quad_n( quad_cv_iloc, quad_cv_gi ) * x( xnod )
               ygi =  ygi + quad_n( quad_cv_iloc, quad_cv_gi ) * y( xnod )
               if( d3 ) zgi =  zgi + quad_n( quad_cv_iloc, quad_cv_gi ) * z( xnod )
            end do

            if( d3 ) then ! local coords for the triangle
               quad_l1( cv_gi ) = volume_quad_map( 1, xgi, ygi, zgi, lx, ly, lz )
               quad_l2( cv_gi ) = volume_quad_map( 2, xgi, ygi, zgi, lx, ly, lz )
               quad_l3( cv_gi ) = volume_quad_map( 3, xgi, ygi, zgi, lx, ly, lz )
               quad_l4( cv_gi ) = volume_quad_map( 4, xgi, ygi, zgi, lx, ly, lz )
            else
               quad_l1( cv_gi ) = area_quad_map( 1, xgi, ygi, lx, ly )
               quad_l2( cv_gi ) = area_quad_map( 2, xgi, ygi, lx, ly )
               quad_l3( cv_gi ) = area_quad_map( 3, xgi, ygi, lx, ly )
            end if

         end do Loop_NGI

      end do Loop_Elements

      ! Now determine the basis functions and derivatives at the 
      ! quadrature pts quad_L1, quad_L2, quad_L3, quad_L4, etc
      call shatri( quad_l1, quad_l2, quad_l3, quad_l4, rdummy, d3, &
           cv_nloc, cv_ngi, n, nlx, nly, nlz )


      deallocate( quad_l1 )
      deallocate( quad_l2 )
      deallocate( quad_l3 )
      deallocate( quad_l4 )
      deallocate( quad_cvweight )
      deallocate( detwei )
      deallocate( ra )
      deallocate( quad_n )
      deallocate( quad_nlx )
      deallocate( quad_nly )
      deallocate( quad_nlz )
      deallocate( quad_nx )
      deallocate( quad_ny )
      deallocate( quad_nz )
      deallocate( rdummy )

      return
    end subroutine shape_tri_tet


    subroutine shatri( l1, l2, l3, l4, weight, d3, &
         nloc, ngi, &
         n, nlx, nly, nlz )
      implicit none
      integer, intent( in ) :: ngi, nloc
      real, dimension( ngi ), intent( in ) :: l1, l2, l3, l4, weight
      logical, intent( in ) :: d3
      real, dimension( nloc, ngi ), intent( inout ) :: n, nlx, nly, nlz
      ! Local variables
      integer :: gi, ndim, cv_ele_type_dummy, u_nloc_dummy
      real, dimension( :, : ), allocatable :: cvn_dummy, un_dummy, unlx_dummy, &
           unly_dummy, unlz_dummy
      real, dimension( : ), allocatable :: cvweigh_dummy


      Conditional_Dimensionality: if( d3 ) then ! Assume a triangle

         Conditional_NLOC: Select Case( nloc )
         case( 6, 7 )
            do gi = 1, ngi           
               n( 1, gi ) = ( 2. * l1( gi ) - 1. ) * l1( gi )            
               n( 2, gi ) = ( 2. * l2( gi ) - 1. ) * l2( gi )             
               n( 3, gi ) = ( 2. * l3( gi ) - 1. ) * l3( gi )             
               n( 4, gi ) = 4. * l1( gi ) * l2( gi )         
               n( 5, gi ) = 4. * l2( gi ) * l3( gi )         
               n( 6, gi ) = 4. * l1( gi ) * l3( gi )
               ! x-derivative (nb. l1 + l2 + l3 + l4 = 1 )                
               nlx( 1, gi ) = 4. * l1( gi ) - 1.              
               nlx( 2, gi ) = 0.           
               nlx( 3, gi ) = -4. * ( 1. - l2( gi ) ) + 4. * l1( gi ) + 1.                            
               nlx( 4, gi ) = 4 * l2( gi )            
               nlx( 5, gi ) = - 4 * l2( gi )          
               nlx( 6, gi ) = 4. * ( 1. - l2( gi ) ) -8. * l1( gi )            
               ! y derivative
               nly( 1, gi ) = 0.
               nly( 2, gi ) = 4. * l2( gi ) - 1.
               nly( 3, gi ) = -4. * ( 1. - l1( gi ) ) + 4. * l2( gi ) + 1. 
               nly( 4, gi ) = 4. * l1( gi )
               nly( 5, gi ) = 4. * ( 1. - l1( gi ) ) - 8. * l2( gi )
               nly( 6, gi ) = -4 * ( l1( gi ) )   
               if( nloc == 7 ) then  ! Bubble function
                  n( 7, gi ) = l1( gi ) * l2( gi ) * l3( gi )
                  nlx( 7, gi ) = l2( gi ) * ( 1. - l2( gi ) )  &
                       -2. * l1( gi ) * l2( gi )
                  nly( 7, gi ) = l1( gi ) * ( 1. - l1( gi ) )  &
                       -2. * l1( gi ) * l2( gi )
               end if
            end do

         case( 3, 4 )
            do gi = 1, ngi 
               n( 1, gi ) = l1( gi )         
               n( 2, gi ) = l2( gi )         
               n( 3, gi ) = l3( gi )
               ! x-derivative            
               nlx( 1, gi ) = 1.                
               nlx( 2, gi ) = 0.             
               nlx( 3, gi ) = -1.
               ! y-derivative 
               nly( 1, gi ) = 0.          
               nly( 2, gi ) = 1.         
               nly( 3, gi ) = -1.
               ! nloc = 4
               if( nloc == 4 ) then
                  n( 4, gi ) = l1( gi ) * l2( gi ) * l3( gi )                
                  nlx( 4, gi ) = l2( gi ) * ( 1. - l2( gi )) - 2. * l1( gi ) * l2( gi )            
                  nly( 4, gi ) = l1( gi ) * ( 1. - l1( gi )) - 2. * l1( gi ) * l2( gi )              
               end if
            end do

         case( 1 )
            do gi = 1, ngi 
               n( 1, gi ) = 1.
               nlx( 1, gi ) = 0.
               nly( 1, gi ) = 0.
            end do

         case default; FLExit( " Wrong number of NLOC " )
         end Select Conditional_NLOC

      else ! Assume it is a tetrahedron 
         Conditional_NLOC2: Select Case( nloc )
         case( 10, 11 ) ! 5 points quadrature
            !call quadtetshapes( n, nlx, nly, nlz, nloc, ngi )
            allocate( cvn_dummy( nloc, ngi ) )
            allocate( un_dummy( nloc, ngi ) )
            allocate( unlx_dummy( nloc, ngi ) )
            allocate( unly_dummy( nloc, ngi ) )
            allocate( unlz_dummy( nloc, ngi ) )
            allocate( cvweigh_dummy( ngi ) )
            u_nloc_dummy = nloc
            cv_ele_type_dummy = 2
            ndim = 3
            call quad_nd_shape( ndim, cv_ele_type_dummy, ngi, nloc, &
                 u_nloc_dummy, cvn_dummy, cvweigh_dummy, &
                 n, nlx, nly, nlz, &
                 un_dummy, unlx_dummy, unly_dummy, unlz_dummy )
            !call quad_nd_shape( ndim, cv_ele_type_dummy, cv_ngi, cv_nloc, u_nloc_dummy, cvn_dummy, cvweigh_dummy, &
!                 n, nlx, nly, nlz, &
!                 un_dummy, unlx_dummy, unly_dummy, unlz_dummy )
            if( nloc == 11 ) then ! Bubble function
               n( 11, : ) = l1 * l2 * l3 * l4
               nlx( 11, : ) = l2 * l3 * ( -l2 - l3 + 1. ) - 2. * l1 * l2 * l3
               nly( 11, : ) = l1 * l3 * ( -l1 - l3 + 1. ) - 2. * l1 * l2 * l3
               nlz( 11, : ) = l1 * l2 * ( -l1 - l2 + 1. ) - 2. * l1 * l2 * l3
            end if

         case( 4, 5 )
            do gi = 1, ngi
               n( 1, gi ) = l1( gi )
               n( 2, gi ) = l2( gi )
               n( 3, gi ) = l3( gi )
               n( 4, gi ) = l4( gi )
               ! x-derivative       
               nlx( 1, gi ) = 1.         
               nlx( 2, gi ) = 0.       
               nlx( 3, gi ) = 0.     
               nlx( 4, gi ) = -1.
               ! y-derivative 
               nly( 1, gi ) = 0.
               nly( 2, gi ) = 1.
               nly( 3, gi ) = 0.
               nly( 4, gi ) = -1.
               ! z-derivative 
               nlz( 1, gi ) = 0.
               nlz( 2, gi ) = 0.
               nlz( 3, gi ) = 1.
               nlz( 4, gi ) = -1.
               if( nloc == 5 ) then ! Bubble function
                  n( 5, gi ) = l1( gi ) * l2( gi ) * l3( gi ) * l4( gi ) 
                  nlx( 5, gi ) = l2( gi ) * l3( gi ) * ( 1. - l2( gi ) - l3( gi ))  &
                       -2. * ( l1( gi ) * l2( gi ) * l3( gi ) )
                  nly( 5, gi ) = l1( gi ) * l3( gi ) * ( 1. - l1( gi ) - l3( gi ))  &
                       -2. * ( l1( gi ) * l2( gi ) * l3( gi ) )
                  nlz( 5, gi ) = l1( gi ) * l2( gi ) * ( 1. - l2( gi ) - l2( gi ))  &
                       -2. * ( l1( gi ) * l2( gi ) * l3( gi ) )
               endif
            end do

         case( 1 )
            do gi = 1, ngi
               n( 1, gi ) = 1.
               nlx( 1, gi ) = 0.
               nly( 1, gi ) = 0.
               nlz( 1, gi ) = 0.
            end do

         case default; FLExit( "Wrong number for NLOC " )

         end Select Conditional_NLOC2

      end if Conditional_Dimensionality

      deallocate( cvn_dummy )
      deallocate( un_dummy )
      deallocate( unlx_dummy )
      deallocate( unly_dummy )
      deallocate( unlz_dummy )
      deallocate( cvweigh_dummy )
      return
    end subroutine shatri

    subroutine calc_cvn_tri_tet( cv_ele_type, totele, cv_nloc, cv_ngi, x_nonods, &
         quad_cv_nloc, x_ndgln, fem_nod, cvn )
      ! Compute CVN (CV basis function)
      implicit none
      integer, intent( in ) :: cv_ele_type, totele, cv_nloc, cv_ngi, &
           x_nonods, quad_cv_nloc
      integer, dimension( totele * quad_cv_nloc ), intent( in ) :: x_ndgln
      integer, dimension( x_nonods ), intent( in ) :: fem_nod
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: cvn
      ! Local variables
      integer :: ele, nod, quad_cv_iloc, xnod, quad_cv_gi, cv_gi

      cvn = 0.
      Loop_Elements: do ele = 1, totele
         nod = 0
         do quad_cv_iloc = 1, quad_cv_nloc
            xnod = x_ndgln( ( ele - 1 ) * quad_cv_nloc + quad_cv_iloc )
            if( fem_nod( xnod ) /= 0 ) nod = fem_nod( xnod )
         end do

         if( nod == 0 ) FLExit(" Problem with CVN calculation " )

         do quad_cv_gi = 1, cv_ngi
            cv_gi = ( ele - 1 ) * cv_ngi + quad_cv_gi
            cvn( nod, cv_gi ) = 1.
         end do

      end do Loop_Elements

      return
    end subroutine calc_cvn_tri_tet


    subroutine shape_l_q_quad( lowqua, ngi, nloc, mloc, &
         sngi, snloc, smloc, m, mlx, mly, mlz, &
         weight, n, nlx, nly, nlz, &
         sweigh, sn, snlx, snly, sm, smlx, smly, &
         nwicel, d3 )  
      use shape_functions_Linear_Quadratic
      ! This subrt computes shape functions. For now, let's just 
      ! define for one element type.
      ! NB: N may overwrite M if we are not solving for pressure.
      implicit none
      logical, intent( in ) :: lowqua
      integer, intent( in ) :: ngi, nloc, mloc, sngi, snloc, smloc
      real, dimension( mloc, ngi ), intent( inout ) :: m, mlx, mly, mlz
      real, dimension( ngi ), intent( inout ) :: weight
      real, dimension( nloc, ngi ), intent( inout ) :: n, nlx, nly, nlz
      real, dimension( sngi ), intent( inout ) :: sweigh
      real, dimension( snloc, sngi ), intent( inout ) :: sn, snlx, snly
      real, dimension( smloc, sngi ), intent( inout ) :: sm, smlx, smly
      integer, intent( in ) :: nwicel
      logical, intent( in ) :: d3
      ! Local variables
      integer :: ipoly, iqadra

      Select Case( nwicel )

      case default; FLExit( "Wrong option for NWICEL " )

      case( 1 )
         if( .not. d3 ) then
            call re2dn4( lowqua, ngi, nloc, mloc, &
                 m, weight, n, nlx, nly, &
                 sngi, snloc, sweigh, sn, snlx )
         else
            call re3dn8( lowqua, ngi, nloc, mloc, &
                 m, weight, n, nlx, nly, nlz, &
                 sngi, snloc, sweigh, sn, snlx, snly ) 
         end if

      case( 3 )
         if( .not. d3 ) then
            call re2dn9( lowqua, ngi, nloc, mloc, &
                 m, weight, n, nlx, nly )
         else ! Lagrange 27 nodes 3D element - bilinear pressure
            call re3d27( lowqua, ngi, nloc, mloc, &
                 m, weight, n, nlx, nly, nlz )
         end if

      end Select

      return
    end subroutine shape_l_q_quad


!!!-                           -!!!
!!!-     And the functions:    -!!!
!!!-                           -!!!


    real function volume_quad_map( cv_iloc, xgi, ygi, zgi, lx, ly, lz )
      ! Compute the cv_iloc^{th} shape function value at point (xgi, ygi, zgi)
      implicit none
      integer :: cv_iloc
      real :: xgi, ygi, zgi
      integer, parameter :: n = 3
      real, dimension( n ) :: lx, ly, lz
      ! Local variables
      integer :: ipt
      real, dimension( :, : ), allocatable :: lxyz
      real, dimension( : ), allocatable :: xyzgi
      real :: loc_x_coord, loc_y_coord, loc_z_coord, loc_zz_coord

      allocate( lxyz( n + 1, n ) )
      allocate( xyzgi( n ) )

      do ipt = 1, n + 1
         lxyz( ipt, 1 ) = lx( ipt ) 
         lxyz( ipt, 2 ) = ly( ipt )
         lxyz( ipt, 3 ) = lz( ipt )
      end do

      xyzgi( 1 ) = xgi
      xyzgi( 2 ) = ygi
      xyzgi( 3 ) = zgi

      ! Local coordinates:
      loc_x_coord = tet_vol( lxyz( 2, : ), lxyz( 3, : ), lxyz( 4, : ), xyzgi )
      loc_y_coord = tet_vol( lxyz( 1, : ), lxyz( 3, : ), lxyz( 4, : ), xyzgi )
      loc_z_coord = tet_vol( lxyz( 1, : ), lxyz( 2, : ), lxyz( 4, : ), xyzgi )
      loc_zz_coord = 1. - ( loc_x_coord + loc_y_coord + loc_z_coord )

      Select Case( cv_iloc )
      case( 1 ); volume_quad_map = loc_x_coord
      case( 2 ); volume_quad_map = loc_y_coord
      case( 3 ); volume_quad_map = loc_z_coord
      case( 4 ); volume_quad_map = loc_zz_coord
      case default; FLExit( " Wrong number of cv_iloc " )
      end Select

      deallocate( lxyz )
      deallocate( xyzgi )

      return
    end function volume_quad_map


    real function area_quad_map( cv_iloc, xgi, ygi, lx, ly )
      implicit none
      integer :: cv_iloc
      real :: xgi, ygi
      integer, parameter :: n = 3
      real, dimension( n ) :: lx, ly
      ! Local variables
      real :: loc_x_coord, loc_y_coord, loc_z_coord

      loc_x_coord = triareaf( lx( 2 ), ly( 2 ), lx( 3 ), ly( 3 ), xgi, ygi )
      loc_y_coord = triareaf( lx( 1 ), ly( 1 ), lx( 3 ), ly( 3 ), xgi, ygi )
      loc_z_coord = 1. - ( loc_x_coord + loc_y_coord )

      Select Case( cv_iloc )
      case( 1 ); area_quad_map = loc_x_coord
      case( 2 ); area_quad_map = loc_y_coord
      case( 3 ); area_quad_map = loc_z_coord
      case default; FLExit( " Wrong number of cv_iloc " )
      end Select

      return
    end function area_quad_map

    real function tet_vol( a, b, c, d ) 
      implicit none
      integer, parameter :: n = 3
      real, dimension( n ) :: a, b, c, d
      ! Local variables
      real, dimension( : ), allocatable :: am, bm, cm, cp

      allocate( am( n ))
      allocate( bm( n ))
      allocate( cm( n ))
      allocate( cp( n ))

      am = a - d
      bm = b - d
      cm = c - d
      call CrossProduct( n, cp, bm, cm )

      tet_vol = 1. / 6. * abs( dot_product( am, cp ) )

      deallocate( am )
      deallocate( bm )
      deallocate( cm )

      return
    end function tet_vol

    real function triareaf( x1, y1, x2, y2, x3, y3 )
      implicit none
      real :: x1, y1, x2, y2, x3, y3

      triareaf = 0.5 * abs( ( x2 * y3 - y2 * x3 ) - x1 * ( y3 - y2 ) + y1 * ( x3 - x2 ) )

      return
    end function triareaf

    subroutine CrossProduct( n, cp, a, b )
      implicit none
      integer, intent( in ) :: n
      real, dimension( n ), intent( inout ) :: cp
      real, dimension( n ), intent( in ) :: a, b

      cp( 1 ) = a( 2 ) * b( 3 ) - a( 3 ) * b( 2 )
      cp( 2 ) = a( 3 ) * b( 1 ) - a( 1 ) * b( 3 )
      cp( 3 ) = a( 1 ) * b( 2 ) - a( 2 ) * b( 1 )

      return
    end subroutine CrossProduct

  end module shape_functions_NDim

