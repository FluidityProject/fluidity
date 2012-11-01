  
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
!!!!      SHAPE FUNCTIONS SUBRTS FOR MULTI-D      !!!!
!!!!   (Quadrilaterals, Triangles, Hexaedra and   !!!!
!!!!                 Tetrahedra)                  !!!!
!!!!==============================================!!!!
  

  module shape_functions_Linear_Quadratic

    use fldebug
    use state_module
    use spud
    use global_parameters, only: option_path_len

  contains

    subroutine re2dn4( lowqua, ngi, ngi_l, nloc, mloc, &
         m, weight, n, nlx, nly, &
         sngi, snloc, sweigh, sn, snlx, &
         l1, l2 )
      ! This subrt computes shape functions M and N and their derivatives
      ! at the Gauss points.
      ! NB: We may need to define surface elements for p and (u,v,w) 
      implicit none
      logical, intent( in ) :: lowqua
      integer, intent( in ) :: ngi, ngi_l, nloc, mloc
      real, dimension( mloc, ngi ), intent( inout ) :: m
      real, dimension( ngi ), intent( inout ) :: weight
      real, dimension( nloc, ngi ), intent( inout ) :: n, nlx, nly
      integer, intent( in ) :: sngi, snloc
      real, dimension( sngi ), intent( inout ) :: sweigh
      real, dimension( snloc, sngi ), intent( inout ) :: sn, snlx
      real, dimension( ngi_l ), intent( in ) :: l1, l2
      ! Local variables:
      integer, parameter :: nl = 16, nlp = 4, npq = 2
      real, dimension( : ), allocatable :: lx, ly, lxp, lyp, weit
      integer :: q, p, gpoi, corn, ndgi, i, nd_quad_1d
      real :: posi, rlx, rly
      logical :: getndp

      ewrite(3,*)' In re2dn4 subrt. '
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

      Conditional_NGI: if ( ( ngi == 4 ) .or. ( ngi == 1 ) ) then
         posi = 1. / sqrt( 3. )
         lx( 1 ) = -posi
         ly( 1 ) = -posi
         lx( 2 ) = posi
         ly( 2 ) = posi
         nd_quad_1d = 2
         if ( ngi == 1 ) then
            lx( 1 ) = 0. ; ly( 1 ) = 0. ; nd_quad_1d = 1
         end if
         Loop_Q: do q = 1, nd_quad_1d
            Loop_P: do p = 1, nd_quad_1d
               Loop_Corn: do corn = 1, nlp
                  gpoi = ( q - 1 ) * nd_quad_1d + p
                  rlx = lx( p )
                  rly = ly( q )
                  if( ngi_l /= 0 ) then
                     rlx = l1( gpoi )
                     rly = l2( gpoi )
                  endif
                  if( mloc == 1 ) m( 1, gpoi ) = 1.
                  weight( gpoi ) = 1.
                  if( ngi == 1 ) weight( gpoi ) = 4. ! It's a convolution of two weights of 2.
                  n( corn, gpoi ) = 0.25 * ( 1. + lxp( corn ) * rlx ) * &
                       ( 1. + lyp( corn ) * rly )
                  nlx( corn, gpoi ) = 0.25 * lxp( corn ) * ( 1. + lyp( corn ) * rly )
                  nly( corn, gpoi ) = 0.25 * lyp( corn ) * ( 1. + lxp( corn ) * rlx )
               end do Loop_Corn
            end do Loop_P
         end do Loop_Q

         nd_quad_1d = 2
         if ( sngi == 1 ) then
            lx( 1 ) = 0. ; ly( 1 ) = 0. ; nd_quad_1d = 1
         end if

         if( ( sngi == 1 ) .and. ( snloc > 1 ) ) then ! Surface shape functions
            Loop_P2: do p = 1, nd_quad_1d
               Loop_Corn2: do corn = 1, 2
                  gpoi = p
                  sn( corn, gpoi ) = 0.5 * ( 1. + lxp( corn ) * lx( p ) )
                  ewrite(3,*)'sn::::',corn, gpoi,sn( corn, gpoi )
                  snlx( corn, gpoi ) = 0.5 * lxp( corn )
                  ewrite(3,*)'snlx::::',corn, gpoi,snlx( corn, gpoi )
                  sweigh( gpoi ) = 1.
                  if( sngi == 1 ) sweigh( gpoi ) = 2.
               end do Loop_Corn2
            end do Loop_P2
         end if

      else ! If ngi =/ 4
         ndgi = int( sqrt( ngi + 0.1 ) + 0.1 )
         getndp = .false.
         allocate( weit( ndgi ) )
         call lagrot( weit, lx, ndgi, getndp )   
         ly( 1 : ndgi ) = lx( 1 : ndgi )
         Loop_Q3: do q = 1, ndgi
            Loop_P3: do p = 1, ndgi
               Loop_Corn3: do corn = 1, nlp
                  gpoi = ( q - 1 ) * ndgi + p
                  rlx = lx( p )
                  rly = ly( q )
                  if( ngi_l /= 0 ) then
                     rlx = l1( gpoi )
                     rly = l2( gpoi )
                  endif
                  if( mloc == 1 ) m( 1, gpoi ) = 1.
                  weight( gpoi ) = weit( p ) * weit( q )
                  ! n( corn, gpoi ) = 0.25 * ( 1. + lxp( corn ) * lxp( p ) * &
                  n( corn, gpoi ) = 0.25 * ( 1. + lxp( corn ) * rlx * &
                       ( 1. + lyp( corn ) * rly ) )
                  nlx( corn, gpoi ) = 0.25 * lxp( corn ) * &
                       ( 1. + lyp( corn ) * rly )
                  nly( corn, gpoi ) = 0.25 * lyp( corn ) * &
                       ( 1. + lxp( corn ) * rlx )
               end do Loop_Corn3
            end do Loop_P3
         end do Loop_Q3
         deallocate( weit )

         if( sngi > 0 ) then
            getndp = .false.
            allocate( weit( sngi ) ) 
            call lagrot( weit, lx, sngi, getndp )
            Loop_P4: do p = 1, sngi
               Loop_Corn4: do corn = 1, npq
                  gpoi = p
                  sn( corn, gpoi ) = 0.5 * ( 1. + lxp( corn ) * lx( p ) )
                  snlx( corn, gpoi ) = 0.5 * lxp( corn )
                  sweigh( gpoi ) = weit( p )
               end do Loop_Corn4
            end do Loop_P4
            deallocate( weit )
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

    subroutine re3dn8( lowqua, ngi, ngi_l, nloc, mloc, &
         m, weight, n, nlx, nly, nlz, &
         sngi, snloc, sweigh, sn, snlx, snly, &
         l1, l2, l3 ) 
      ! This subrt. computes the shape functions M and N and their
      ! derivatives at the Gauss points for 3D.
      ! If LOWQUA, then use one point quadrature else use 8 point quadrature.
      !NB.: LX/YP(I) are the local X/Y coordinates of nodal point I.
      implicit none
      logical, intent( in ) :: lowqua
      integer, intent( in ) :: ngi, ngi_l, nloc, mloc
      real, dimension( mloc, ngi ), intent( inout ) :: m
      real, dimension( ngi ), intent( inout ) :: weight
      real, dimension( nloc, ngi ), intent( inout ) :: n, nlx, nly, nlz
      integer, intent( in ) :: sngi, snloc
      real, dimension( sngi ), intent( inout ) :: sweigh
      real, dimension( snloc, sngi ), intent( inout ) :: sn, snlx, snly
      real, dimension( ngi_l ), intent( in ) :: l1, l2, l3
      ! Local variables:
      integer, parameter :: nl = 4, nlp = 8, npq = 2
      real, dimension( : ), allocatable :: lx, ly, lz, lxp, lyp, lzp, weit, rdum2
      integer :: p, q, ir, corn, gpoi, ngi1d, gi
      real :: posi, rdum, rlx, rly, rlz

      ewrite(3,*)' In re3dn8 subrt. '
      ! Allocating memory
      allocate( lx( nl ) )
      allocate( ly( nl ) )
      allocate( lz( nl ) )
      allocate( lxp( nlp ) )
      allocate( lyp( nlp ) )
      allocate( lzp( nlp ) )
      allocate( weit( nl ) )
      allocate( rdum2( ngi_l ) )

      if( sngi >= 1 ) then ! Surface integrals
         call re2dn4( lowqua, sngi, 0, snloc, 0, &
              m, sweigh, sn, snlx, snly, &
              0, 0, sweigh, sn, snlx, &
              rdum2, rdum2 )
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
                     gpoi = ( q - 1 ) * npq + p  + ( ir - 1 ) * npq * npq
                     rlx = lx( p )
                     rly = ly( q )
                     rlz = lz( ir )
                     if( ngi_l /= 0 ) then
                        rlx = l1( gpoi )
                        rly = l2( gpoi )
                        rlz = l3( gpoi )
                     endif
                     if( mloc > 0 ) m( 1, gpoi ) = 1.
                     ! Weight
                     weight( gpoi ) = 1.
                     ! N
                     n( corn, gpoi ) = 0.125 * ( 1. + lxp( corn ) * rlx ) * &
                          ( 1. + lyp( corn ) * rly ) * ( 1. + lzp( corn ) * &
                          rlz )
                     ! x-derivative
                     nlx( corn, gpoi ) = 0.125 * lxp( corn ) * ( 1. + &
                          lyp( corn ) * rly ) * ( 1. + lzp( corn ) * rlz ) 
                     ! y-derivative
                     nly( corn, gpoi ) = 0.125 * lyp( corn ) * ( 1. + lxp( corn ) * &
                          rlx ) * ( 1. + lzp( corn ) * rlz )
                     ! z-derivative
                     nlz( corn, gpoi ) = 0.125 * lzp( corn ) * ( 1. + lxp( corn ) * &
                          rlx ) * ( 1. + lyp( corn ) * rly )
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
                     rlx = lx( p )
                     rly = ly( q )
                     rlz = lz( ir )
                     if( ngi_l /= 0 ) then
                        rlx = l1( gpoi )
                        rly = l2( gpoi )
                        rlz = l3( gpoi )
                     endif
                     if( mloc > 0 ) m( 1, gpoi ) = 1.
                     ! Weight
                     weight( gpoi ) = weit( p ) * weit( q ) * weit( ir )
                     ! N
                     n( corn, gpoi ) = 0.125 * ( 1. + lxp( corn ) * rlx ) * &
                          ( 1. + lyp( corn ) * rly ) * ( 1. + lzp( corn ) * rlz )
                     ! x-derivative
                     nlx( corn, gpoi ) = 0.125 * lxp( corn ) * ( 1. + lyp( corn ) * &
                          rly ) * ( 1. + lzp( corn ) * rlz )
                     ! y-derivative
                     nly( corn, gpoi ) = 0.125 * lyp( corn ) * ( 1. + lxp( corn ) * &
                          rlx ) * ( 1. + lzp( corn ) * rlz )
                     ! z-derivative
                     nlz( corn, gpoi ) = 0.125 * lzp( corn ) * ( 1. + lxp( corn ) * &
                          rlx ) * ( 1. + lyp( corn ) * rly )
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


    subroutine re2dn9( lowqua, ngi, ngi_l, nloc, mloc, &
         m, weight, n, nlx, nly, &
         l1, l2 )
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
      integer, intent( in ) :: ngi, ngi_l, nloc, mloc
      real, dimension( mloc, ngi ), intent( inout ) :: m
      real, dimension( ngi ), intent( inout ) :: weight
      real, dimension( nloc, ngi ), intent( inout ) :: n, nlx, nly
      real, dimension( ngi_l ), intent( in ) :: l1, l2
      ! Local variables:
      integer, parameter :: nl = 3, nlp = 9, npq = 4
      real, dimension( : ), allocatable :: lx, ly, lz, lxp, lyp, weit, &
           xn, yn, zn, dxn, dyn, dzn
      integer :: nquad, p, q, corn, gpoi, ilx, ily, nj
      real :: posi, rlx, rly

      ewrite(3,*) 'In re2dn9'

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

      Conditional_LOWQUA:if( ngi == 1 ) then
         lx( 1 ) = 0.0
         ly( 1 ) = 0.0
         lz( 1 ) = 0.0
         weit( 1 ) = 2.
         nquad = 1 
      else if( lowqua .or. (ngi == 4) ) then
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
               rlx = lx( p )
               rly = ly( q )
               if( ngi_l /= 0 ) then
                  rlx = l1( gpoi )
                  rly = l2( gpoi )
               endif
               m( corn, gpoi ) = 0.25 * ( 1. + lxp( corn ) * rlx ) * &
                    ( 1 + lyp( corn ) * rly )
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
            rlx = lx( p )
            rly = ly( q )
            if( ngi_l /= 0 ) then
               rlx = l1( gpoi )
               rly = l2( gpoi )
            endif
            ! Weight
            weight( gpoi ) = weit( p ) * weit( q )
            ! XN
            xn( 1 ) = 0.5 * rlx * ( rlx - 1. )
            xn( 2 ) = 1. - rlx * rlx
            xn( 3 ) = 0.5 * rlx * ( rlx + 1. )
            ! DXDN
            dxn( 1 ) = 0.5 * ( 2. * rlx - 1. )
            dxn( 2 ) = -2. * rlx
            dxn( 3 ) = 0.5 * ( 2. * rlx + 1. )
            ! YN
            yn( 1 ) = 0.5 * rly * ( rly - 1. )
            yn( 2 ) = 1. - rly * rly
            yn( 3 ) = 0.5 * rly * ( rly + 1. )
            ! DYDN
            dyn( 1 ) = 0.5 * ( 2. * rly - 1. )
            dyn( 2 ) = -2. * rly
            dyn( 3 ) = 0.5 * ( 2. * rly + 1. )
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

    subroutine re3d27( lowqua, ngi, ngi_l, nloc, mloc, &
         m, weight, n, nlx, nly, nlz, &
         l1, l2, l3 )
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
      integer, intent( in ) :: ngi, ngi_l, nloc, mloc
      real, dimension( mloc, ngi ), intent( inout ) :: m
      real, dimension( ngi ), intent( inout ) :: weight
      real, dimension( nloc, ngi ), intent( inout ) :: n, nlx, nly, nlz
      real, dimension( ngi_l ), intent( in ) :: l1, l2, l3
      ! Local variables:
      integer, parameter :: nl = 3, nlp = 27, npq = 4
      real, dimension( : ), allocatable :: lx, ly, lz, lxp, lyp, lzp, weit, &
           xn, yn, zn, dxn, dyn, dzn
      integer :: nquad, p, q, ir, corn, gpoi, ilx, ily, ilz, nj
      real :: posi

      ewrite(3,*)'In re3d27'

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

      Conditional_LOWQUA: if((ngi == 1) ) then
         if(ngi_l /= 0) then
            lx( 1 ) = l1(1)
            ly( 1 ) = l2(1)
            lz( 1 ) = l3(1)
         else
            lx( 1 ) = 0.0
            ly( 1 ) = 0.0
            lz( 1 ) = 0.0
         endif
         ! this should be 2, but the
         ! 1pt quadrature gets the 
         ! volume wrong so adjust it...
         weit( 1 ) = 2.0
         nquad=1
      else if( lowqua .or. (ngi == 8) ) then
         posi = 1. / sqrt( 3. )
         if( ngi_l /= 0 ) then
            lx( 1 ) = l1(1)
            ly( 1 ) = l2(1)
            lz( 1 ) = l3(1)
            lx( 2 ) = l1(2)
            ly( 2 ) = l2(2)
            lz( 2 ) = l3(2)
         else
            lx( 1 ) = -posi
            ly( 1 ) = -posi
            lz( 1 ) = -posi
            lx( 2 ) = posi
            ly( 2 ) = posi
            lz( 2 ) = posi
         end if
         weit( 1 ) = 1.
         weit( 2 ) = 1.
         weit( 3 ) = 1.
         nquad = 2
      else
         posi = 0.774596669241483
         if( ngi_l /= 0 ) then
            lx( 1 ) = l1(1)
            ly( 1 ) = l2(1)
            lz( 1 ) = l3(1)
            lx( 2 ) = l1(2)
            ly( 2 ) = l2(2)
            lz( 2 ) = l3(2)
            lx( 3 ) = l1(3)
            ly( 3 ) = l2(3)
            lz( 3 ) = l3(3)
         else
            lx( 1 ) = -posi
            ly( 1 ) = -posi
            lz( 1 ) = -posi
            lx( 2 ) = 0.
            ly( 2 ) = 0.
            lz( 2 ) = 0.
            lx( 3 ) = posi
            ly( 3 ) = posi
            lz( 3 ) = posi
         end if
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

      Loop_P2: do ir = 1, nquad
         Loop_Q2: do q = 1, nquad
            Loop_IR2: do p = 1, nquad
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
               zn( 3 ) = 0.5 * lz( ir ) * ( lz( ir ) + 1. )
               ! DZDN
               dzn( 1 ) = 0.5 * ( 2. * lz( ir ) - 1. )
               dzn( 2 ) = -2. * lz( ir )
               dzn( 3 ) = 0.5 * ( 2. * lz( ir ) + 1. )
               ! N, NLX and NLY
               Loop_ILX2: do ilz = 1, nl
                  Loop_ILY2: do ily = 1, nl
                     Loop_ILZ2: do ilx = 1, nl
                        nj = ilx + ( ily - 1 ) * 3 + ( ilz - 1 ) * 9
                        n( nj, gpoi ) = xn( ilx ) * yn( ily ) * zn( ilz )
                        !ewrite(3,*)'nj, gpoi,n( nj, gpoi ),xn( ilx ), yn( ily ),zn( ilz ):', &
                        !     nj, gpoi,n( nj, gpoi ),xn( ilx ), yn( ily ),zn( ilz )
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
         cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )
      implicit none
      integer, intent( in ) :: ndim, cv_ele_type, cv_nloc, u_nloc
      ! If QUAD_OVER_WHOLE_ELE=.true. then dont divide element into CV's to form quadrature.
      logical, intent( in ) :: QUAD_OVER_WHOLE_ELE
      integer, intent( inout ) :: cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface
      ! Local variables
      ! volume_order & surface_order are the volume and surface 
      ! order of the integration used within each sub-quad/hex for CV approach.
      ! Default value -ve or 0 value is always 1pt quadrature (=1). 
      integer, PARAMETER :: volume_order=1
      !      integer, PARAMETER :: volume_order=2
      integer, PARAMETER :: surface_order=1
      !      integer, PARAMETER :: surface_order=2
      ! whole_ele_volume_order & whole_ele_surface_order are the volume and surface 
      ! order of the integration used within each sub-quad/hex for QUAD_OVER_WHOLE_ELE=.true.
      ! -ve or 0 and take on the default value. 
      integer, PARAMETER :: whole_ele_volume_order=0
      !      integer, PARAMETER :: whole_ele_volume_order=1
      !      integer, PARAMETER :: whole_ele_volume_order=2
      integer, PARAMETER :: whole_ele_surface_order=0
      !      integer, PARAMETER :: whole_ele_surface_order=1
      !      integer, PARAMETER :: whole_ele_surface_order=2
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
            if(QUAD_OVER_WHOLE_ELE) then
               cv_ngi = 3
               sbcvngi = 2 
               scvngi = 2
               if (whole_ele_volume_order==1) cv_ngi = 1
               if (whole_ele_surface_order==1) sbcvngi = 1
               if (whole_ele_surface_order==1) scvngi = 1
               if (whole_ele_volume_order==2) cv_ngi = 3
               if (whole_ele_surface_order==2) sbcvngi = 2
               if (whole_ele_surface_order==2) scvngi = 2
            else
               if (volume_order==1) cv_ngi = 3
               if (volume_order==2) cv_ngi = 3*4

               if (surface_order==1) scvngi = 3
               if (surface_order==1) sbcvngi = 2
               if (surface_order==2) scvngi = 3*2
               if (surface_order==2) sbcvngi = 2*2
            endif
         case( 6 ) ! Quadratic Triangle
            if(QUAD_OVER_WHOLE_ELE) then
               cv_ngi = 7
               sbcvngi = 3
               scvngi = 3
               if (whole_ele_volume_order==1) cv_ngi = 1
               if (whole_ele_surface_order==1) sbcvngi = 1
               if (whole_ele_surface_order==1) scvngi = 1
               if (whole_ele_volume_order==2) cv_ngi = 3
               if (whole_ele_surface_order==2) sbcvngi = 2
               if (whole_ele_surface_order==2) scvngi = 2
               if (whole_ele_volume_order==3) cv_ngi = 7
               if (whole_ele_surface_order==3) sbcvngi = 3
               if (whole_ele_surface_order==3) scvngi = 3
            else
               if (volume_order==1) cv_ngi = 12
               if (volume_order==2) cv_ngi = 12*4

               if (surface_order==1) scvngi = 12
               if (surface_order==1) sbcvngi = 4 
               if (surface_order==2) scvngi = 12*2
               if (surface_order==2) sbcvngi = 4*2
            endif
         case default; FLExit(" Invalid integer for cv_nloc ")
         end Select Conditional_CV_NLOC2D_Tri
         nface = 3

      case( 5, 6 ) ! Quads
         Conditional_CV_NLOC2D_Quad: Select Case( cv_nloc )
         case( 4 ) ! Bi-linear Quad
            if(QUAD_OVER_WHOLE_ELE) then
               cv_ngi = 4
               scvngi = 2
               sbcvngi = 2 
               if (whole_ele_volume_order==1) cv_ngi = 1
               if (whole_ele_surface_order==1) sbcvngi = 1
               if (whole_ele_surface_order==1) scvngi = 1
               if (whole_ele_volume_order==2) cv_ngi = 4
               if (whole_ele_surface_order==2) sbcvngi = 2
               if (whole_ele_surface_order==2) scvngi = 2
               if (whole_ele_volume_order==3) cv_ngi = 9
               if (whole_ele_surface_order==3) sbcvngi = 3
               if (whole_ele_surface_order==3) scvngi = 3
            else
               if (volume_order==1) cv_ngi = 4
               if (volume_order==2) cv_ngi = 16

               if (surface_order==1) scvngi = 4
               if (surface_order==2) scvngi = 4*2
               if (surface_order==1) sbcvngi = 2
               if (surface_order==2) sbcvngi = 4 
            endif
         case( 9 ) ! Bi-quad Quad
            if(QUAD_OVER_WHOLE_ELE) then
               cv_ngi = 9
               sbcvngi = 3
               scvngi = 3
               if (whole_ele_volume_order==1) cv_ngi = 1
               if (whole_ele_surface_order==1) sbcvngi = 1
               if (whole_ele_surface_order==1) scvngi = 1
               if (whole_ele_volume_order==2) cv_ngi = 4
               if (whole_ele_surface_order==2) sbcvngi = 2
               if (whole_ele_surface_order==2) scvngi = 2
               if (whole_ele_volume_order==3) cv_ngi = 9
               if (whole_ele_surface_order==3) sbcvngi = 3
               if (whole_ele_surface_order==3) scvngi = 3
            else
               if (volume_order==1) cv_ngi = 16
               if (volume_order==2) cv_ngi = 16*4
               if (volume_order==3) cv_ngi = 16*9

               if (surface_order==1) scvngi = 16
               if (surface_order==1) sbcvngi = 4
               if (surface_order==2) scvngi = 16*2
               if (surface_order==2) sbcvngi = 4*2
               if (surface_order==3) scvngi = 16*3
               if (surface_order==3) sbcvngi = 4*3
            endif
         case default; FLExit(" Invalid integer for cv_nloc ")
         end Select Conditional_CV_NLOC2D_Quad
         nface = 4

      case( 7, 8 ) ! Tetrahedra
         Conditional_CV_NLOC3D_Tet: Select Case( cv_nloc )
         case( 4 ) ! Linear 
            if(QUAD_OVER_WHOLE_ELE) then
               cv_ngi = 4
               sbcvngi = 3
               scvngi = 3
               if (whole_ele_volume_order==1) cv_ngi = 1
               if (whole_ele_surface_order==1) sbcvngi = 1
               if (whole_ele_surface_order==1) scvngi = 1
               if (whole_ele_volume_order==2) cv_ngi = 4
               if (whole_ele_surface_order==2) sbcvngi = 3
               if (whole_ele_surface_order==2) scvngi = 3
               if (whole_ele_volume_order==3) cv_ngi = 11
               if (whole_ele_surface_order==3) sbcvngi = 7
               if (whole_ele_surface_order==3) scvngi = 7
            else
               if (volume_order==1) cv_ngi = 4
               if (volume_order==2) cv_ngi = 4*8

               if (surface_order==1) scvngi = 6
               if (surface_order==1) sbcvngi = 3
               if (surface_order==2) scvngi = 6*4
               if (surface_order==2) sbcvngi = 3*4
            endif
         case( 10 ) ! Quadratic
            if(QUAD_OVER_WHOLE_ELE) then
               cv_ngi = 11
               sbcvngi = 7
               scvngi = 7
               if (whole_ele_volume_order==1) cv_ngi = 1
               if (whole_ele_surface_order==1) sbcvngi = 1
               if (whole_ele_surface_order==1) scvngi = 1
               if (whole_ele_volume_order==2) cv_ngi = 4
               if (whole_ele_surface_order==2) sbcvngi = 3
               if (whole_ele_surface_order==2) scvngi = 3
               if (whole_ele_volume_order==3) cv_ngi = 11
               if (whole_ele_surface_order==3) sbcvngi = 7
               if (whole_ele_surface_order==3) scvngi = 7
            else

               if (volume_order==1) cv_ngi=8*4*1 ! (1x1x1)
               if (volume_order==2) cv_ngi=8*4*8 ! (2x2x2)
               if (volume_order==3) cv_ngi=8*4*27 ! (3x3x3)

               if (surface_order==1) scvngi = 48 ! 6x8x1 (cv_faces x tets x sngi)
               if (surface_order==1) sbcvngi = 12 ! 1x12 (sngi x cv_faces)
               if (surface_order==2) scvngi = 48*4 ! 6x8x4 (cv_faces x hexs x sngi)
               if (surface_order==2) sbcvngi = 12*4 ! 4x12 (sngi x cv_faces)

            endif
         case default; FLExit(" Invalid integer for cv_nloc ")
         end Select Conditional_CV_NLOC3D_Tet
         nface = 4

      case( 9, 10 )
         Conditional_CV_NLOC3D_Hex: Select Case( cv_nloc )
         case( 8 ) ! Tri-linear Hex
            if(QUAD_OVER_WHOLE_ELE) then
               cv_ngi = 8
               sbcvngi = 4
               scvngi = 4
               if (whole_ele_volume_order==1) cv_ngi = 1
               if (whole_ele_surface_order==1) sbcvngi = 1
               if (whole_ele_surface_order==1) scvngi = 1
               if (whole_ele_volume_order==2) cv_ngi = 8
               if (whole_ele_surface_order==2) sbcvngi = 4
               if (whole_ele_surface_order==2) scvngi = 4
               if (whole_ele_volume_order==3) cv_ngi = 27
               if (whole_ele_surface_order==3) sbcvngi = 9
               if (whole_ele_surface_order==3) scvngi = 9
            else
               if (volume_order==1) cv_ngi = 8
               if (volume_order==2) cv_ngi = 8*8

               if (surface_order==1) scvngi = 12
               if (surface_order==1) sbcvngi = 4
               if (surface_order==2) scvngi = 12*4
               if (surface_order==2) sbcvngi = 4*4
            endif
         case( 27 ) ! Tri-quad Hex
            if(QUAD_OVER_WHOLE_ELE) then
               cv_ngi = 27
               sbcvngi = 9
               scvngi = 9
               if (whole_ele_volume_order==1) cv_ngi = 1
               if (whole_ele_surface_order==1) sbcvngi = 1
               if (whole_ele_surface_order==1) scvngi = 1
               if (whole_ele_volume_order==2) cv_ngi = 8
               if (whole_ele_surface_order==2) sbcvngi = 4
               if (whole_ele_surface_order==2) scvngi = 4
               if (whole_ele_volume_order==3) cv_ngi = 27
               if (whole_ele_surface_order==3) sbcvngi = 9
               if (whole_ele_surface_order==3) scvngi = 9
            else
               if (volume_order==1) cv_ngi = 64
               if (volume_order==2) cv_ngi = 64*8
               if (volume_order==3) cv_ngi = 64*27

               if (surface_order==1) scvngi = 12*8
               if (surface_order==1) sbcvngi = 16
               if (surface_order==2) scvngi = 12*8*4
               if (surface_order==2) sbcvngi = 16*4
               if (surface_order==3) scvngi = 12*8*9
               if (surface_order==3) sbcvngi = 16*9
            endif
         case default; FLExit(" Invalid integer for cv_nloc ")
         end Select Conditional_CV_NLOC3D_Hex
         nface = 6

      case default; FLExit( " Invalid integer for cv_ele_type " )

      end Select Conditional_EleType

      if(.not.QUAD_OVER_WHOLE_ELE) then
         if( cv_ele_type > 2 ) scvngi = scvngi + nface * sbcvngi
      endif
      cv_ngi_short = cv_ngi

      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) cv_ngi = cv_ngi * cv_nloc

      !         ewrite(3,*)'cv_ele_type,cv_ngicv_ngi_short=', &
      !                  cv_ele_type,cv_ngi,cv_ngi_short
      !         ewrite(3,*)'scvngi,sbcvngi=',scvngi,sbcvngi
      !         stop 393
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

      weit = 0. ; quadpos = 0.

      ewrite(3,*) 'In LAGROT'
      ewrite(3,*) 'ndgi, getndp:', ndgi, getndp
      ewrite(3,*) 'weit:', weit( 1 : ndgi )
      ewrite(3,*) 'quadpos:', quadpos( 1 : ndgi )

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
      INTEGER :: GI, L, IGLX, ii
      logical, save :: first = .true.
      real, save :: rsum, rsumabs
      !
      VOLUME = 0.
      !
      IF(D3) THEN
         if( first ) then
            rsum = 0. ; rsumabs = 0.
            first = .false.
         end if
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
               ewrite(3,*)'xndgln, x, nl:', &
                    iglx, l, x(iglx), y(iglx), z(iglx), NLX(L,GI), NLY(L,GI), NLZ(L,GI)
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
            DETWEI(GI)=ABS(DETJ)*WEIGHT(GI)
            RA(GI)=1.0
            VOLUME=VOLUME+DETWEI(GI)
            ewrite(3,*)'gi, detj, weight(gi)', gi, detj, weight(gi)
            rsum = rsum + detj
            rsumabs = rsumabs + abs( detj )
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
         !ewrite(3,*)'ele, sum(detj), sum(abs(detj)):', ele, rsum, rsumabs
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
            DETWEI(GI)=TWOPIE*RGI*ABS(DETJ)*WEIGHT(GI)
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
            DETWEI( GI ) = ABS(DETJ) * WEIGHT( GI )
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
      real, dimension( 0 : ndnod - 1 ) :: nodpos
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

      ewrite(3,*)'In RGPTWE', ig, nd, weight

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
         ewrite(3,*)'in rgptwe, nd:', nd
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

      !  RETURN

    END FUNCTION RGPTWE



    subroutine quad_basis_funs_1d(sngi, snloc,  &
         sweigh, sn, snlx )

      ! determine the 1d shape functions sn and its local derivative slnx. 
      implicit none
      integer, intent( in ) :: sngi, snloc
      real, dimension( snloc, sngi ), intent( inout ) :: sn, snlx
      real, dimension( sngi ), intent( inout ) :: sweigh
      ! local variables...
      integer :: iloc, gpoi
      real :: lxgp
      logical :: diff, ndiff, getndp
      real, dimension( : ), allocatable :: quad_pts, snodpos, rdummy

      allocate(quad_pts(sngi))
      allocate(snodpos(snloc))
      allocate(rdummy(snloc))

      ! Find the roots of the quadrature points and nodes
      ! also get the weights.
      diff = .true.
      ndiff = .false.

      getndp = .true.
      call lagrot( rdummy, snodpos, snloc, getndp )

      getndp = .false.
      !     Compute standard Gauss quadrature: weights and points
      call lagrot( sweigh, quad_pts, sngi, getndp )

      do iloc = 1, snloc
         do gpoi=1,sngi
            lxgp=quad_pts(gpoi)
            sn( iloc, gpoi )   = lagran( ndiff, lxgp, iloc, snloc, snodpos )
            snlx( iloc, gpoi ) = lagran( diff,  lxgp, iloc, snloc, snodpos )
         end do
      end do
      deallocate(quad_pts)
      deallocate(snodpos)
      deallocate(rdummy)
      return
    end subroutine quad_basis_funs_1d



  end module shape_functions_Linear_Quadratic

  
  module shape_functions_NDim

    use fldebug
    use shape_functions_Linear_Quadratic

  contains

    subroutine quad_1d_shape( cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, n, nlx, un, unlx )
      ! For quadratic elements. Shape functions associated with volume integration 
      ! using both CV basis functions CVN as well as FEM basis functions N (and 
      ! its derivatives NLX, NLY, NLZ)
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
      ewrite( 3, * ) 'cv_ngi, cv_nloc', cv_ngi, cv_nloc

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
      ewrite(3,*) 'cv_ngi, cv_nloc, u_nloc::', &
           cv_ngi, cv_nloc, u_nloc

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
         cv_ngi_1d = int( ( 1.e-3 + cv_ngi ) ** (1./3.) )
         cv_nloc_1d = int( ( 1.e-3 + cv_nloc ) ** (1./3.) )
         u_nloc_1d = int( ( 1.e-3 + u_nloc ) ** (1./3.) )

      case default; FLExit( " Invalid integer for NDIM " )

      end Select Conditional_Dimensionality

      ! Allocating memory
      allocate( cvweigh_1d( cv_ngi_1d ) )
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


    subroutine vol_cv_tri_tet_shape( cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, &
         n, nlx, nly, nlz, &
         un, unlx, unly, unlz )
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
      integer, dimension( : ), allocatable :: x_ndgln, fem_nod, x_ndgln_ideal
      real, dimension( : ), allocatable :: lx, ly, lz, x, y, z, cvweigh_dummy, &
           x_ideal, y_ideal, z_ideal
      integer, parameter :: max_totele = 10000, max_x_nonods = 10000
      logical :: d1, dcyl, d3
      integer :: ele, quad_cv_ngi, quad_cv_nloc, totele, x_nonods, &
           cv_gj, cv_gk, cv_iloc, cv_gi, totele_sub
      real :: rsum

      ewrite(3,*)'In vol_cv_tri_tet_shape'

      d1 = ( ndim == 1 )
      dcyl = .false.
      d3 = ( ndim == 3 )

      if( d3 ) then
         Select Case( cv_nloc )
         case( 4 ) ;  quad_cv_nloc = 8  ! Linear tetrahedron
         case( 10 ) ; quad_cv_nloc = 27 ! Quadratic tetrahedron
         case default; FLExit( "Wrong integer for CV_NLOC" )
         end Select
      else
         Select Case( cv_nloc )
         case( 3 ) ;  quad_cv_nloc = 4  ! Linear triangle
         case( 6 ) ; quad_cv_nloc = 9 ! Quadratic triangle
         case default; FLExit( "Wrong integer for CV_NLOC" )
         end Select
      endif

      ! Allocating memory
      allocate( lx( max_x_nonods ) ) ; lx = 0.
      allocate( ly( max_x_nonods ) ) ; ly = 0.
      allocate( lz( max_x_nonods ) ) ; lz = 0.
      allocate( x( max_x_nonods )) ; x = 0.
      allocate( y( max_x_nonods )) ; y =0.
      allocate( z( max_x_nonods )) ; z = 0.
      allocate( fem_nod( max_x_nonods ) ) ; fem_nod = 0
      allocate( x_ndgln( max_totele * quad_cv_nloc ) ) ; x_ndgln = 0
      allocate( cvweigh_dummy( cv_ngi )) ; cvweigh_dummy = 0.
      allocate( x_ideal( max_x_nonods ) ) ; x_ideal = 0.
      allocate( y_ideal( max_x_nonods ) ) ; y_ideal = 0.
      allocate( z_ideal( max_x_nonods ) ) ; z_ideal = 0.
      allocate( x_ndgln_ideal( max_x_nonods ) ) ; x_ndgln_ideal = 0

      ! Get the x_ndgln for the nodes of the triangle or tet or hex/quad super-elements:
      x = 0. ; y = 0. ; z = 0. ; lx = 0. ; ly = 0. ; lz = 0. ; fem_nod = 0 ; x_ndgln = 0
      call Compute_XNDGLN_TriTetQuadHex( cv_ele_type, &
           max_totele, max_x_nonods, quad_cv_nloc, &
           totele, x_nonods, &
           x_ndgln, lx, ly, lz, x, y, z, fem_nod, &
           x_ideal, y_ideal, z_ideal, x_ndgln_ideal )

      ! Compute the shape functions using these quadrilaterals/hexs:
      ! For pressure:
      call shape_tri_tet( cv_ele_type, cv_nloc, &
           cv_ele_type, ndim, totele, cv_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, &
           n, nlx, nly, nlz, cvweigh )

      ! Compute cvn: 
      call Calc_CVN_TriTetQuadHex( cv_ele_type, totele, cv_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, fem_nod, cvn )

      ewrite(3,*)'cvweigh:',cvweigh
      do cv_iloc=1,cv_nloc
         rsum=0.0
         do cv_gi=1,cv_ngi
            rsum=rsum+cvn(cv_iloc,cv_gi)*cvweigh(cv_gi)
         end do
         ewrite(3,*)'cv_iloc,rsum:',cv_iloc,rsum
      end do
      !stop 2922
      if(cv_nloc==10) then
        ewrite(3,*)'cv_nloc=',cv_nloc
         totele_sub=8
         !call test_quad_tet( cv_nloc, cv_ngi, cvn, n, nlx, nly, nlz, &
         call test_quad_tet( cv_nloc, cv_ngi, cvn, n, nlx, nly, nlz, &
              cvweigh, x_ideal, y_ideal, z_ideal, cv_nloc, x_ndgln_ideal, 1) 
      endif

      ! And for velocities:
      if(u_nloc==1) then ! constant basis function throughout element...
         un=1.0
         unlx=0.0
         unly=0.0
         unlz=0.0
      else
         call shape_tri_tet( cv_ele_type, cv_nloc, &
              cv_ele_type, ndim, totele, u_nloc, cv_ngi, x_nonods, &
              quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, &
              un, unlx, unly, unlz, cvweigh_dummy )
      endif

      deallocate( lx )
      deallocate( ly )
      deallocate( lz )
      deallocate( x )
      deallocate( y ) 
      deallocate( z ) 
      deallocate( x_ndgln )
      deallocate( fem_nod )
      deallocate( cvweigh_dummy )

      return
    end subroutine vol_cv_tri_tet_shape



     subroutine test_quad_tet( cv_nloc, cv_ngi, cvn, n, nlx, nly, nlz, &
                       cvweight, x, y, z, x_nonods, x_ndgln2, totele )
! test the volumes of idealised triangle 
      implicit none
      integer, intent( in ) :: cv_nloc, cv_ngi, x_nonods, totele
      real, dimension( x_nonods ), intent( in ) :: x, y, z
      integer, dimension( totele * cv_nloc ), intent( in ) :: x_ndgln2
      real, dimension( cv_ngi ), intent( in ) :: cvweight
      real, dimension( cv_nloc,cv_ngi ), intent( in ) :: cvn
      REAL, DIMENSION( CV_NLOC, CV_NGI ), intent( in ) :: N, NLX, NLY, NLZ 
! local variables...
      integer, dimension( : ), allocatable :: x_ndgln
      real, dimension( : ), allocatable :: DETWEI,RA
      real, dimension( :,: ), allocatable :: NX,NY,NZ
      INTEGER :: NDIM, ele, cv_iloc, cv_gi
      LOGICAL :: D1,D3,DCYL
      REAL :: VOLUME, rsum, rsum2

      ALLOCATE( DETWEI( CV_NGI )) 
      ALLOCATE( x_ndgln( CV_Nloc )) 
      do cv_iloc=1,cv_nloc
          x_ndgln( CV_iloc )=cv_iloc 
      end do
      ALLOCATE( RA( CV_NGI ))
      ALLOCATE( NX( CV_NLOC, CV_NGI ))
      ALLOCATE( NY( CV_NLOC, CV_NGI ))
      ALLOCATE( NZ( CV_NLOC, CV_NGI ))



      ndim=3
      D1 = ( NDIM == 1 )
      D3 = ( NDIM == 3 )
      DCYL = .FALSE. 

      RSUM=0.0
      Loop_Elements: DO ELE = 1, TOTELE

         ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
         CALL DETNLXR( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, CV_NLOC, CV_NGI, &
              N, NLX, NLY, NLZ, CVWEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
              NX, NY, NZ ) 
         EWRITE(3,*)'ele, VOLUME=',ele, VOLUME
         ewrite(3,*)'detwei:',detwei
         EWRITE(3,*)'sum of detwei:',sum(detwei)
         ewrite(3,*)'nx=',nx
         ewrite(3,*)'nlx=',nlx
         RSUM=RSUM+VOLUME
         
         do cv_iloc=1,cv_nloc
            rsum2=0.0
            do cv_gi=1,cv_ngi
               rsum2=rsum2+cvn(cv_iloc,cv_gi)*detwei(cv_gi)
               ewrite(3,*)'cv_gi,cvn(cv_iloc,cv_gi),detwei(cv_gi),CVWEIGHT(cv_gi):', &
                        cv_gi,cvn(cv_iloc,cv_gi),detwei(cv_gi),CVWEIGHT(cv_gi)
            end do
            ewrite(3,*)'cv_iloc,its vol:',cv_iloc,rsum2
         end do
      END DO Loop_Elements

         do cv_iloc=1,cv_nloc
            ewrite(3,*)'cv_iloc,nod,X, Y, Z:',cv_iloc,X_NDGLN(cv_iloc), &
                    X(X_NDGLN(cv_iloc)), Y(X_NDGLN(cv_iloc)), Z(X_NDGLN(cv_iloc))
         end do

      EWRITE(3,*)'VOLUME OF THE DOMAIN(SHOULD BE 1):',RSUM
    
      !STOP 2992

      return
     end subroutine test_quad_tet



    subroutine Compute_XNDGLN_TriTetQuadHex( cv_ele_type, &
         max_totele, max_x_nonods, quad_cv_nloc, &
         totele, x_nonods, &
         x_ndgln, lx, ly, lz, x, y, z, fem_nod, &
         x_ideal, y_ideal, z_ideal, x_ndgln_ideal  )
      ! Get the x_ndgln for the nodes of triangles or tetrahedra
      implicit none
      integer, intent( in ) :: cv_ele_type, max_totele, max_x_nonods, quad_cv_nloc
      integer, intent( inout ) :: totele, x_nonods
      real, dimension( max_x_nonods ), intent( inout ) :: lx, ly, lz, x, y, z
      integer, dimension( max_x_nonods ), intent( inout ) :: fem_nod
      integer, dimension( max_totele * quad_cv_nloc ), intent( inout ) :: x_ndgln
      real, dimension( max_x_nonods ), intent( inout ) :: x_ideal, y_ideal, z_ideal
      integer, dimension( max_x_nonods ), intent( inout ) :: x_ndgln_ideal
      ! Local variables
      integer, dimension( : ), allocatable :: x_ndgln2, x_ndgln_big
      real, dimension( : ), allocatable :: x2, y2, z2
      integer :: quad_u_loc_dummy, ele, quad_cv_iloc, quad_cv_siloc, cv_iloc, &
           x_nonods2, npoly, npoly_x, npoly_y, npoly_z, nelex, neley, nelez, &
           elex, eley, elez, ip, jp, kp, xcount, xnod,  xnod2, quad_cv_nloc2, &
           x_nloc
      real :: h_scale, intvalsx, intvalsy, intvalsz, xstar, xfini, ystar, yfini, &
           zstar, zfini, dis2
      logical :: d1, d3, found

      fem_nod = 0
      quad_u_loc_dummy = 1

      Conditional_ElementTypes: Select Case( cv_ele_type )
      case( 1, 2, 5, 6, 9, 10 ) ! Quadrilaterals and Hexahedra
         d1 = ( ( cv_ele_type == 1 ) .or. ( cv_ele_type == 2 ) )
         d3 = ( ( cv_ele_type == 9 ) .or. ( cv_ele_type == 10 ) )
         npoly_y = 1
         npoly_z = 1
         npoly = 2 ! linear
         if( ( cv_ele_type == 2 ) .or. ( cv_ele_type == 6 ) .or. &
              ( cv_ele_type == 10 ) ) npoly = 3
         npoly_x = npoly

         if( .not. d1 ) npoly_y = npoly
         if( d3 ) npoly_z = npoly
         nelex = npoly_x ; neley = npoly_y ; nelez = npoly_z
         totele = nelex * neley * nelez
         quad_cv_nloc2 = npoly_x * npoly_y * npoly_z
         x_nonods2 = totele * quad_cv_nloc2

         allocate( x2( x_nonods2 ) )
         allocate( y2( x_nonods2 ) )
         allocate( z2( x_nonods2 ) )

         lx( 1 ) = -1.0
         if( .not. d1 ) ly( 1 ) = -1.0
         if( d3 ) lz( 1 ) = -1.0

         lx( 2 ) =  1.0
         if( .not. d1 ) ly( 2 ) = -1.0
         if( d3 ) lz( 2 ) = -1.0

         if( x_nonods >= 3 ) then
            lx( 3 ) = -1.0
            if( .not. d1 ) ly( 3 ) =  1.0
            if( d3 ) lz( 3 ) = -1.0
         endif

         if( x_nonods >= 4 ) then
            lx( 4 ) =  1.0
            if( .not. d1 ) ly( 4 ) =  1.0
            if( d3 ) lz( 4 ) = -1.0
         endif

         if( d3 ) then
            lx( 5 ) = -1.0
            ly( 5 ) = -1.0
            lz( 5 ) =  1.0

            lx( 6 ) =  1.0
            ly( 6 ) = -1.0
            lz( 6 ) =  1.0

            lx( 7 ) = -1.0
            ly( 7 ) =  1.0
            lz( 7 ) =  1.0

            lx( 8 ) =  1.0
            ly( 8 ) =  1.0
            lz( 8 ) =  1.0
         end if


         ! Divide each element into sub-elements for each CV:
         intvalsx = max( 2 * ( nelex - 2 ) + 2, 1 )
         intvalsy = max( 2 * ( neley - 2 ) + 2, 1 )
         intvalsz = max( 2 * ( nelez - 2 ) + 2, 1 )

         Loop_Poly2d_2_1x: do elex = 1, nelex
            Loop_Poly2d_2_2y: do eley = 1, neley
               Loop_Poly2d_2_3z: do elez = 1, nelez
                  ele= ( elez - 1 ) * neley * nelex + ( eley - 1 ) * nelex + elex 
                  xstar = max( -1.0, -1.0 - 2.0 / real( intvalsx ) + 2.0 * &
                       intvalsx * ( elex - 1 ) / real( nelex ) )
                  xfini = min( 1.0, -1.0 - 2.0 / real( intvalsx ) + 2.0 * intvalsx * &
                       ( elex ) / real( nelex ) )

                  if( .not. d1 ) then
                     ystar = max( -1.0, -1.0 - 2.0 / real( intvalsy ) + &
                          2.0 * intvalsy * ( eley - 1 ) / real( neley ) )
                     yfini = min( 1.0, -1.0 - 2.0 / real( intvalsy ) + &
                          2.0 * intvalsy * ( eley ) / real( neley ) )
                  end if
                  if( d3 ) then
                     zstar = max( -1.0, -1.0 - 2.0 / real( intvalsz ) + &
                          2.0 * intvalsz * ( elez - 1 ) / real( nelez ) )
                     zfini = min( 1.0, -1.0 - 2.0 / real( intvalsz ) + &
                          2.0 * intvalsz * ( elez ) / real( nelez ) )
                  end if

                  ! Place a mesh across each sub-element
                  Loop_Poly2d_2_1: do ip = 1, npoly_x
                     Loop_Poly2d_2_2: do jp = 1, npoly_y
                        Loop_Poly2d_2_3: do kp = 1, npoly_z

                           quad_cv_iloc = ( kp - 1 ) * npoly_y * npoly_x + &
                                ( jp - 1 ) * npoly_x + ip 
                           cv_iloc = ( ele - 1 ) * quad_cv_nloc2 + quad_cv_iloc
                           x2 ( cv_iloc ) = xstar + ( xfini - xstar ) * real( ip - 1 ) / &
                                real( npoly_x - 1 )
                           if( .not. d1 ) y2( cv_iloc ) = ystar + ( yfini - ystar ) * &
                                real( jp - 1 ) / real( npoly_y - 1 )
                           if( d3 ) z2( cv_iloc ) = zstar + ( zfini - zstar ) * &
                                real( kp - 1 ) / real( npoly_y - 1 )
                           quad_cv_siloc = quad_cv_siloc + 1
                           x_ndgln( cv_iloc ) = &
                                cv_iloc
                        end do Loop_Poly2d_2_3
                     end do Loop_Poly2d_2_2
                  end do Loop_Poly2d_2_1

               end do Loop_Poly2d_2_3z
            end do Loop_Poly2d_2_2y
         end do Loop_Poly2d_2_1x

         ! Remove repetition of nodal points 
         xcount = 0
         do xnod = 1, x_nonods2
            found = .false.
            do xnod2 = 1, xnod - 1
               dis2 = abs( x2( xnod ) - x2( xnod2 ) )
               if( .not. d1 ) dis2 = dis2 + abs( y2( xnod ) - y2( xnod2 ) )
               if( d3 ) dis2 = dis2 + abs( z2( xnod ) - z2( xnod2 ) ) 
               if( dis2 < 1.e-4 ) found = .true.
            end do
            if( .not. found ) then
               xcount = xcount + 1
               x( xcount) = x2( xnod )
               y( xcount) = y2( xnod )
               z( xcount) = z2( xnod )
               x_ndgln( xnod ) = xcount 
            endif
         end do
         x_nonods = xcount

         deallocate( x2 )
         deallocate( y2 )
         deallocate( z2 )

      case( 3 ) ! Linear Triangles 

         totele = 3
         x_nonods = 7
         h_scale =  1.5196713713031851 ! Scaling factor to give a unity area of the local triangle

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

      case( 4 ) ! Quadratic Triangles 
         x_nonods = max_x_nonods ! x_nonods will be assessed and updated in Make_QTri
         totele = 12
         x = 0. ; y = 0. ; z = 0.
         x_ndgln = 0
         call Make_QTri( totele, quad_cv_nloc, max_x_nonods, x_nonods, &
              x_ndgln, lx, ly, x, y, fem_nod )

         !! Just debugging local numbering:
         ewrite(3,*)'Just out of Make_QTri'
         ewrite(3,*)'cv_ele_type, totele, x_nonods, quad_cv_nloc :', &
              cv_ele_type, totele, x_nonods, quad_cv_nloc
         ewrite(3,*)'fem_nod:', ( fem_nod( ele ), ele = 1, 6 )
         ewrite(3,*)'lx:', ( lx( ele ), ele = 1, 3 )
         ewrite(3,*)'ly:', ( ly( ele ), ele = 1, 3 )
         ewrite(3,*)'x_ndgln:'
         do ele = 1, totele
            ewrite(3,*) ele, ( x_ndgln( ( ele - 1 ) * quad_cv_nloc + cv_iloc ), &
                 cv_iloc = 1, quad_cv_nloc )
         end do

         ewrite(3,*)'X / Y / Z'
         do ele = 1, totele
            do cv_iloc = 1, quad_cv_nloc
               xnod = x_ndgln( ( ele - 1 ) * quad_cv_nloc + cv_iloc )
               ewrite(3,*) ele, cv_iloc, xnod, x( xnod ), y( xnod ), z ( xnod )
            end do
         end do

      case( 7 ) ! Linear Tetrahedra
         x_nonods = 15
         totele = 4
         h_scale = 2.0396489026555056 ! Scaling factor to give a unity volume of the local tetrahedron

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

         lx( 4 ) = 0.5 * ( lx( 1 ) + lx( 2 ) )
         ly( 4 ) = 0.5 * ( ly( 1 ) + ly( 3 ) ) 
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
         y( 15 ) = ly( 4 )
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
         ele = 1

         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 9
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 8
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 12
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 11

         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 5 ) = 1
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 6 ) = 2
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 7 ) = 4
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 8 ) = 5

         ele = 2

         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 8
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 10
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 11
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 13

         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 5 ) = 2
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 6 ) = 3
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 7 ) = 5
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 8 ) = 6

         ele = 3

         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 12
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 11
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 14
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 13

         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 5 ) = 4
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 6 ) = 5
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 7 ) = 7
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 8 ) = 6

         ele = 4

         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 1 ) = 15
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 2 ) = 10
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 3 ) = 14
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 4 ) = 13

         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 5 ) = 9
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 6 ) = 8
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 7 ) = 12
         x_ndgln( ( ele - 1 ) * quad_cv_nloc + 8 ) = 11

      case( 8 ) ! Quadratic Tetrahedra 

         totele = 8 ;  x_nonods = max_x_nonods ; x_nloc = 4
         call  Make_QTets( totele, quad_cv_nloc, x_nloc, max_x_nonods, x_nonods, &
              x_ndgln, lx, ly, lz, x, y, z, fem_nod, &
              x_ideal, y_ideal, z_ideal, x_ndgln_ideal )

         !         call get_x_ndgln_qtet( max_x_nonods, max_totele, quad_cv_nloc, &
         !              x_nonods, totele, fem_nod, x_ndgln, &
         !              x, y, z, lx, ly, lz )

      case default; FLExit( "Wrong integer for CV_ELE_TYPE" )
      end Select Conditional_ElementTypes

      ewrite(3,*) ' Leaving Compute_XNDGLN_TriTetQuadHex'

      return
    end subroutine Compute_XNDGLN_TriTetQuadHex

    subroutine get_x_ndgln_qtet( max_x_nonods2, max_totele2, quad_cv_nloc2, &
         x_nonods2, totele2, fem_nod, x_ndgln_return, &
         x, y, z, lx, ly, lz )
      use pascal_tetrahedra
      implicit none
      integer, intent( in ) :: max_x_nonods2, max_totele2, quad_cv_nloc2
      integer, intent( inout ) :: x_nonods2, totele2
      integer, dimension( max_x_nonods2 ), intent( inout ) :: fem_nod
      integer, dimension( max_totele2 * quad_cv_nloc2 ), intent( inout ) :: x_ndgln_return
      real, dimension( max_x_nonods2 ), intent( inout ) :: lx, ly, lz, x, y, z
      ! Local variables
      integer, dimension( : ), allocatable :: x_ndgln_big, x_ndgln2
      integer :: triangle_totele, nodeplustetnodes

      totele2 = totele * tet_totele
      triangle_totele = totele * tet_totele
      x_nonods2 = no_of_nodes_in_faces * no_faces * triangle_totele + triangle_totele
      allocate( x_ndgln2( x_nonods2 ) )
      allocate( x_ndgln_big( x_nonods2 ) )

      ! Sanity check
      if( x_nonods2 > max_totele2 * quad_cv_nloc2 ) &
           FLExit( "Dimension for x_ndgln is wrong" )
      
      call numbering_tet_elements( x_nonods2, triangle_totele, nodeplustetnodes, &
           x_ndgln_return, &
           x, y, z, lx, ly, lz, fem_nod, &
           x_ndgln_big, x_ndgln2 ) 

      deallocate( x_ndgln2 )
      deallocate( x_ndgln_big )

      return
    end subroutine get_x_ndgln_qtet



    subroutine suf_cv_tri_tet_shape( cv_ele_type, ndim, scvngi, cv_nloc, u_nloc, scvfeweigh, &
         scvfen, scvfenlx, scvfenly, scvfenlz, scvfenslx, scvfensly,  &
         sufen, sufenlx, sufenly, sufenlz, sufenslx, sufensly, &
         cv_neiloc, cvfem_neiloc, ufem_neiloc )
      !        sn, snlx, snly, snlz, sufnlx, sufnly,  &
      !        sun, sunlx, sunly, sunlz, sufunlx, sufunly  )
      ! Compute shape functions N, UN etc for linear triangles. Shape functions 
      ! associated with volume integration using both CV basis functions CVN, as 
      ! well as FEM basis functions SN (and its derivatives SNLX, SNLY, SNLZ). Also 
      ! for velocity basis functions SUN, SUNLX, SUNLY, SUNLZ.
      ! Also the derivatives along the CV faces: sufnlx, sufnly, sufunlx, sufunly  
      implicit none
      integer, intent( in ) :: cv_ele_type, ndim, scvngi, cv_nloc, u_nloc
      real, dimension( scvngi ), intent( inout ) :: scvfeweigh
      real, dimension( cv_nloc, scvngi ), intent( inout ) :: scvfen, scvfenlx, scvfenly, &
           scvfenlz, scvfenslx, scvfensly
      real, dimension( u_nloc, scvngi ), intent( inout ) :: sufen, sufenlx, sufenly, sufenlz, &
           sufenslx, sufensly
      integer, dimension( cv_nloc, scvngi ), intent( inout ) :: cv_neiloc, cvfem_neiloc
      integer, dimension( u_nloc, scvngi ), intent( inout ) :: ufem_neiloc
      ! Local variables
      integer, parameter :: max_totele = 1000, max_x_nonods = 10000
      integer, dimension( : ), allocatable :: x_ndgln, fem_nod, x_ndgln_ideal
      integer, dimension( :,: ), allocatable :: cv_neiloc_cells_dummy
      real, dimension( : ), allocatable :: lx, ly, lz, x, y, z, scvfeweigh_dummy, &
           x_ideal, y_ideal, z_ideal
      logical :: d1, dcyl, d3
      integer :: x_nonods, totele, ele, cv_iloc, quad_cv_ngi, quad_cv_nloc, inod

      ewrite(3,*)' In vol_cv_tri_shape'

      d1 = ( ndim == 1 )
      dcyl = .false.
      d3 = ( ndim == 3 )

      if( d3 ) then
         Select Case( cv_nloc )
         case( 4 ) ;  quad_cv_nloc = 8  ! Linear tetrahedron
         case( 10 ) ; quad_cv_nloc = 27 ! Quadratic tetrahedron
         case default; FLExit( "Wrong integer for CV_NLOC" )
         end Select
      else
         Select Case( cv_nloc )
         case( 3 ) ;  quad_cv_nloc = 4  ! Linear triangle
         case( 6 ) ; quad_cv_nloc = 9 ! Quadratic triangle
         case default; FLExit( "Wrong integer for CV_NLOC" )
         end Select
      endif

      ! Allocating memory
      allocate( lx( max_x_nonods ) )
      allocate( ly( max_x_nonods ) )
      allocate( lz( max_x_nonods ) )
      allocate( x( max_x_nonods ))
      allocate( y( max_x_nonods ))
      allocate( z( max_x_nonods ))
      allocate( x_ndgln( max_totele * quad_cv_nloc ) )
      allocate( fem_nod( max_x_nonods ) )
      allocate( scvfeweigh_dummy( scvngi ) )
      allocate( cv_neiloc_cells_dummy( cv_nloc, scvngi ) )
      allocate( x_ideal( max_x_nonods ) ) ; x_ideal = 0.
      allocate( y_ideal( max_x_nonods ) ) ; y_ideal = 0.
      allocate( z_ideal( max_x_nonods ) ) ; z_ideal = 0.
      allocate( x_ndgln_ideal( max_x_nonods ) ) ; x_ndgln_ideal = 0

      ! Get the x_ndgln for the nodes of triangles, tetrahedra, quadrilaterals or hexahedra
      ! super-elements:
      x = 0. ; y = 0. ; z = 0. ; lx = 0. ; ly = 0. ; lz = 0. ; fem_nod = 0 ; x_ndgln = 0
      call Compute_XNDGLN_TriTetQuadHex( cv_ele_type, &
           max_totele, max_x_nonods, quad_cv_nloc, &
           totele, x_nonods, &
           x_ndgln, lx, ly, lz, x, y, z, fem_nod, &
           x_ideal, y_ideal, z_ideal, x_ndgln_ideal )

      ! Compute the shape functions using these quadrilaterals and hexahedra:
      ! For pressure:
      call Compute_SurfaceShapeFunctions_Triangle_Tetrahedron( &
           cv_nloc, cv_ele_type, &
           cv_ele_type, ndim, totele, cv_nloc, scvngi, x_nonods, &
           quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, fem_nod, &
           scvfen, scvfenlx, scvfenly, scvfenlz, scvfenslx, scvfensly, &
           scvfeweigh, cv_neiloc, cvfem_neiloc )

      ewrite(3,*)'Shape Functions for scalar fields -- SCVFEN'
      call PrintOutFunMat( cv_nloc, scvngi, scvfen )

      ewrite(3,*)'Shape Functions for scalar fields -- SCVFENLX'
      call PrintOutFunMat( cv_nloc, scvngi, scvfenlx )

      ewrite(3,*)'Shape Functions for scalar fields -- SCVFENLY'
      call PrintOutFunMat( cv_nloc, scvngi, scvfenly )

      ewrite(3,*)'Shape Functions for scalar fields -- SCVFENLZ'
      call PrintOutFunMat( cv_nloc, scvngi, scvfenlz )

      ewrite(3,*)'Shape Functions for scalar fields -- SCVFENSLX'
      call PrintOutFunMat( cv_nloc, scvngi, scvfenslx )

      ewrite(3,*)'Shape Functions for scalar fields -- SCVFENSLY'
      call PrintOutFunMat( cv_nloc, scvngi, scvfensly )

      ewrite(3,*)'Shape Functions for scalar fields -- SCVFEWEIGH'
      ewrite(3,*) ( scvfeweigh( cv_iloc ), cv_iloc = 1, scvngi )

      ! And for velocities:
      if( u_nloc == 1 ) then ! a constant basis function 
         sufen = 1.0 
         sufenlx = 0.0 
         sufenly = 0.0
         sufenlz = 0.0 
         sufenslx = 0.0
         sufensly = 0.0
      else
         call Compute_SurfaceShapeFunctions_Triangle_Tetrahedron( &
              cv_nloc, cv_ele_type, &
              cv_ele_type, ndim, totele, u_nloc, scvngi, x_nonods, &
              quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, fem_nod, &
              sufen, sufenlx, sufenly, sufenlz, sufenslx, sufensly, &
              scvfeweigh_dummy, cv_neiloc_cells_dummy, ufem_neiloc )
      endif

      ewrite(3,*)'Shape Functions for velocity fields -- SUFEN'
      call PrintOutFunMat( u_nloc, scvngi, sufen )

      ewrite(3,*)'Shape Functions for velocity fields -- SUFENLX'
      call PrintOutFunMat( u_nloc, scvngi, sufenlx )

      ewrite(3,*)'Shape Functions for velocity fields -- SUFENLY'
      call PrintOutFunMat( u_nloc, scvngi, sufenly )

      ewrite(3,*)'Shape Functions for velocity fields -- SUFENLZ'
      call PrintOutFunMat( u_nloc, scvngi, sufenlz )

      ewrite(3,*)'Shape Functions for velocity fields -- SUFENSLX'
      call PrintOutFunMat( u_nloc, scvngi, sufenslx )

      ewrite(3,*)'Shape Functions for velocity fields -- SUFENSLY'
      call PrintOutFunMat( u_nloc, scvngi, sufensly )

      ewrite(3,*)'Shape Functions for velocity fields -- SCVFEWEIGH'
      ewrite(3,*) ( scvfeweigh( cv_iloc ), cv_iloc = 1, scvngi )
   
      deallocate( lx )
      deallocate( ly )
      deallocate( lz )
      deallocate( x )
      deallocate( y ) 
      deallocate( z ) 
      deallocate( x_ndgln )
      deallocate( fem_nod )
      deallocate( scvfeweigh_dummy )
      deallocate( cv_neiloc_cells_dummy )

      return
    end subroutine suf_cv_tri_tet_shape


    subroutine Compute_SurfaceShapeFunctions_Triangle_Tetrahedron( &
         cv_nloc_cells, cv_ele_type_cells, &
         cv_ele_type, ndim, totele, cv_nloc, scvngi, x_nonods, &
         quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, fem_nod, &
         sn, snlx, snly, snlz, sufnlx, sufnly, &
         scvweigh, cv_neiloc_cells, cvfem_neiloc )
      ! this subroutine calculates shape functions sn, snlx, snly, snlz, sufnlx, sufnly, 
      ! their weights scvweigh and local connectivity cv_neiloc_cells, cvfem_neiloc 
      ! on the boundaries of the cv_nloc_cells control volumes. 
      ! The control volume types are defines using cv_ele_type_cells. 
      ! cv_neiloc_cells is associated with the CV cells 
      ! cvfem_neiloc is associated with the FEM basis SN etc.
      use shape_functions
      implicit none
      integer, intent( in ) :: cv_nloc_cells, cv_ele_type_cells
      integer, intent( in ) :: cv_ele_type, ndim, totele, cv_nloc, scvngi, &
           x_nonods, quad_cv_nloc
      integer, dimension( totele * quad_cv_nloc ), intent( in ) :: x_ndgln
      integer, dimension( x_nonods ), intent( in ) :: fem_nod
      integer, dimension( cv_nloc_cells, scvngi ), intent( inout ) :: cv_neiloc_cells
      integer, dimension( cv_nloc, scvngi ), intent( inout ) :: cvfem_neiloc
      real, dimension( x_nonods ), intent( in ) :: x, y, z
      real, dimension( quad_cv_nloc ), intent( in ) :: lx, ly, lz
      real, dimension( cv_nloc, scvngi ), intent( inout ) :: sn, snlx, snly, snlz, &
           sufnlx, sufnly
      real, dimension( scvngi ), intent( inout ) ::  scvweigh
      ! Local variables
      logical, dimension( : ), allocatable :: remove_ig_pt
      integer, dimension( : ), allocatable :: x_sndgln, next_to_cv_iloc_gi, &
           next_to_cv_iloc_gi2 
      integer, dimension( :, : ), allocatable :: loc_2nd_lev, cv_neiloc_cells_2, &
           cv_neiloc_cells_3
      real, dimension( : ), allocatable :: l1_2, l2_2, l3_2, l4_2, &
           quad_cvweight, detwei, ra, rdummy, &
           normxn, normyn, normzn, quad_scvweight, quad_sdetwei, sra, &
           loc_coord_nod_l1, loc_coord_nod_l2, loc_coord_nod_l3, loc_coord_nod_l4, &
           gl_quad_l1, gl_quad_l2, gl_quad_l3, gl_quad_l4, gl_quad_scvweigh, &
           xsl, ysl, zsl, scvweigh_2, l1, l2, l3, l4
      real, dimension( :, : ), allocatable :: quad_n, quad_nlx, quad_nly, quad_nlz, &
           quad_nx, quad_ny, quad_nz, sn_i_xj, quad_sn, quad_snlx, quad_snly, &
           quad_snx, quad_sny, quad_sm, quad_smlx, quad_smly, &
           suf_quad_sn, suf_quad_snlx, suf_quad_snly, &
           gl_quad_sn, gl_quad_snlx, gl_quad_snly, gl_quad_snlz, &
           sn_2, snlx_2, snly_2, snlz_2, &
           suf_snlx_2, suf_snly_2, suf_snlx, suf_snly, rdummy2
      logical :: d1, dcyl, d3, d2, lowqua, found, found_fem_nod, &
           zer_l1, zer_l2, zer_l3, zer_l4, tri_tet
      integer :: ele, quad_cv_ngi, quad_cv_gi, cv_gi, nwicel, xnod, quad_cv_iloc, &
           quad_u_loc_dummy, mloc, dummy_sngi, dummy_snloc, dummy_smloc, &
           ip, jp, kp, sele, iface, stotel, nface, npoly, quad_cv_snloc, quad_cv_sngi, &
           cv_iloc, cv_jloc, cv_iloc_belong, cv_ngi_2, cv_sgi, cv_sgj, &
           quad_cv_sgi, quad_cv_siloc, cv_sgk, cv_sngi_2, quad_cv_sjloc, quad_sgi, &
           xnodi, xnodj, nodi, nodj, cv_iloc_cells, cv_jloc_cells, npoly_ngi, icount
      real :: xgi, ygi, zgi, volume, sarea, normx, normy, normz, d2_quad
      real :: half_side_length

      ewrite(3,*)'Compute_SurfaceShapeFunctions_Triangle_Tetrahedron'
      ewrite(3,*)'scvngi=',scvngi
      ewrite(3,*)'lx:', lx( 1 : quad_cv_nloc )
      ewrite(3,*)'ly:', ly( 1 : quad_cv_nloc )
      ewrite(3,*)'lz:', lz( 1 : quad_cv_nloc )

      sn = 0.0
      snlx = 0.0
      snly = 0.0
      snlz = 0.0
      sufnlx = 0.0
      sufnly = 0.0

      d3 = ( ndim == 3 )
      d2 = ( ndim == 2 )
      d1 = ( ndim == 1 )
      dcyl = .false. 

      ! Compute some dummy variables
      ! npoly=2 for linear 2 node elements, =3 for quadratic 3 node elements in 1D. 
      call dummy_tri_tet( d1, d3, quad_cv_ngi, quad_cv_nloc, &
           dummy_sngi, dummy_snloc, nwicel, & 
           cv_nloc_cells, scvngi, totele, quad_u_loc_dummy, &
           mloc, dummy_smloc, lowqua, npoly, npoly_ngi )

      quad_cv_snloc = npoly ** ( ndim - 1 )
      quad_cv_sngi = max(1,npoly_ngi) ** ( ndim - 1 )

      nface = 1
      if( d2 ) nface = 4
      if( d3 ) nface = 6
      stotel = nface
      allocate( x_sndgln( stotel * quad_cv_snloc ) )
      allocate( loc_2nd_lev( stotel, quad_cv_snloc ) )

      ! Allocating memory
      allocate( quad_cvweight( quad_cv_ngi ) ) ; quad_cvweight = 0.
      allocate( detwei( quad_cv_ngi ) )
      allocate( ra( quad_cv_ngi ) )
      allocate( quad_n( quad_cv_nloc, quad_cv_ngi ) ) ; quad_n = 0.
      allocate( quad_nlx( quad_cv_nloc, quad_cv_ngi ) ) ; quad_nlx = 0.
      allocate( quad_nly( quad_cv_nloc, quad_cv_ngi ) ) ; quad_nly = 0.
      allocate( quad_nlz( quad_cv_nloc, quad_cv_ngi ) ) ; quad_nlz = 0.
      allocate( quad_nx( quad_cv_nloc, quad_cv_ngi ) ) ; quad_nx = 0.
      allocate( quad_ny( quad_cv_nloc, quad_cv_ngi ) ) ; quad_ny = 0.
      allocate( quad_nz( quad_cv_nloc, quad_cv_ngi ) ) ; quad_nz = 0.
      allocate( rdummy( 10000 ) ) ; rdummy = 0.
      allocate( rdummy2( 100, 100 ) ) ; rdummy2 = 0.

      ! For surfaces
      allocate( quad_scvweight( quad_cv_sngi ) ) ; quad_scvweight = 0.
      allocate( quad_sdetwei( quad_cv_sngi ) ) ; quad_sdetwei = 0.
      allocate( sra( quad_cv_sngi ) ) ; sra = 0.
      allocate( quad_sn( quad_cv_snloc, quad_cv_sngi ) ) ; quad_sn = 0.
      allocate( quad_snlx( quad_cv_snloc, quad_cv_sngi ) ) ; quad_snlx = 0.
      allocate( quad_snly( quad_cv_snloc, quad_cv_sngi ) ) ; quad_snly = 0.
      allocate( quad_snx( quad_cv_snloc, quad_cv_sngi ) ) ; quad_snx = 0.
      allocate( quad_sny( quad_cv_snloc, quad_cv_sngi ) ) ; quad_sny = 0.
      allocate( quad_sm( quad_cv_nloc, quad_cv_sngi ) ) ; quad_sm = 0.
      allocate( quad_smlx( quad_cv_snloc, quad_cv_sngi ) ) ; quad_smlx = 0.
      allocate( quad_smly( quad_cv_snloc, quad_cv_sngi ) ) ; quad_smly = 0.

      allocate( loc_coord_nod_l1( x_nonods ) ) ; loc_coord_nod_l1 = 0.
      allocate( loc_coord_nod_l2( x_nonods ) ) ; loc_coord_nod_l2 = 0.
      allocate( loc_coord_nod_l3( x_nonods ) ) ; loc_coord_nod_l3 = 0.
      allocate( loc_coord_nod_l4( x_nonods ) ) ; loc_coord_nod_l4 = 0.
      allocate( sn_i_xj( cv_nloc, x_nonods ) ) ; sn_i_xj = 0.

      allocate( suf_quad_sn( cv_nloc, stotel * quad_cv_sngi * totele ) ) ; suf_quad_sn = 0.
      allocate( suf_quad_snlx( cv_nloc, stotel * quad_cv_sngi * totele ) ) ; suf_quad_snlx = 0.
      allocate( suf_quad_snly( cv_nloc, stotel * quad_cv_sngi * totele ) ) ; suf_quad_snly = 0.

      allocate( gl_quad_l1( stotel * quad_cv_sngi * totele ) ) ; gl_quad_l1 = 0.
      allocate( gl_quad_l2( stotel * quad_cv_sngi * totele ) ) ; gl_quad_l2 = 0.
      allocate( gl_quad_l3( stotel * quad_cv_sngi * totele ) ) ; gl_quad_l3 = 0.
      allocate( gl_quad_l4( stotel * quad_cv_sngi * totele ) ) ; gl_quad_l4 = 0.

      allocate( gl_quad_sn( cv_nloc, stotel * quad_cv_sngi * totele ) )
      allocate( gl_quad_snlx( cv_nloc, stotel * quad_cv_sngi * totele ) )
      allocate( gl_quad_snly( cv_nloc, stotel * quad_cv_sngi * totele ) )
      allocate( gl_quad_snlz( cv_nloc, stotel * quad_cv_sngi * totele ) )
      allocate( gl_quad_scvweigh( stotel * quad_cv_sngi * totele ) )

      allocate( remove_ig_pt( stotel * quad_cv_sngi * totele ) )
      allocate( next_to_cv_iloc_gi( stotel * quad_cv_sngi * totele ) )
      allocate( next_to_cv_iloc_gi2( stotel * quad_cv_sngi * totele ) )
      next_to_cv_iloc_gi = 0 ; next_to_cv_iloc_gi2 = 0

      allocate( normxn( quad_cv_sngi ) ) 
      allocate( normyn( quad_cv_sngi ) ) 
      allocate( normzn( quad_cv_sngi ) )

      allocate( xsl( quad_cv_snloc ) ) 
      allocate( ysl( quad_cv_snloc ) ) 
      allocate( zsl( quad_cv_snloc ) ) 

      allocate( l1( stotel * quad_cv_sngi * totele ) ) ; l1 = 0.
      allocate( l2( stotel * quad_cv_sngi * totele ) ) ; l2 = 0.
      allocate( l3( stotel * quad_cv_sngi * totele ) ) ; l3 = 0.
      allocate( l4( stotel * quad_cv_sngi * totele ) ) ; l4 = 0.

      normx = 1. ; normy = 0. ; normz = 0.

      ! Now we need to compute QUAD_NLX/Y/Z - get the hex or quad
      ! shape functions quad_n etc.
      ewrite(3,*) 'before shape_l_q_quad lowqua, quad_cv_ngi, mloc:', &
           lowqua, quad_cv_ngi, mloc
      ewrite(3,*) 'dummy_sngi, dummy_snloc, dummy_smloc, quad_cv_ngi:', &
           dummy_sngi, dummy_snloc, dummy_smloc, quad_cv_ngi
       ewrite(3,*)'scvngi,quad_cv_sngi=',scvngi,quad_cv_sngi
!        stop 331

      ! Work out local coords of the nodes
      loc_coord_nod_l1 = 0. ; loc_coord_nod_l2 = 0. ; loc_coord_nod_l3 = 0. ; &
           loc_coord_nod_l4 = 0.
      do xnod = 1, x_nonods
         if( d3 ) then 
            loc_coord_nod_l1( xnod ) = &
                 volume_quad_map( 1, x( xnod ), y( xnod ), z( xnod ), lx(1:4), ly(1:4), lz(1:4) )
            loc_coord_nod_l2( xnod ) = &
                 volume_quad_map( 2, x( xnod ), y( xnod ), z( xnod ), lx(1:4), ly(1:4), lz(1:4) )
            loc_coord_nod_l3( xnod ) = &
                 volume_quad_map( 3, x( xnod ), y( xnod ), z( xnod ), lx(1:4), ly(1:4), lz(1:4) )
            loc_coord_nod_l4( xnod ) = &
                 volume_quad_map( 4, x( xnod ), y( xnod ), z( xnod ), lx(1:4), ly(1:4), lz(1:4) )
         else if( .not. d1 ) then ! 2D
            loc_coord_nod_l1( xnod ) = &
                 area_quad_map( 1, x( xnod ), y( xnod ), lx, ly )
            loc_coord_nod_l2( xnod ) = &
                 area_quad_map( 2, x( xnod ), y( xnod ), lx, ly )
            loc_coord_nod_l3( xnod ) = &
                 area_quad_map( 3, x( xnod ), y( xnod ), lx, ly )
         else 
            loc_coord_nod_l1( xnod ) = 1.0
         end if
      end do

      ! Get the shape functions on lines (in 2D) and quadrilateria surfaces in 3D: 
      call shape_l_q_quad( lowqua, quad_cv_ngi, quad_cv_nloc, mloc, &
           quad_cv_sngi, quad_cv_snloc, dummy_smloc, rdummy, rdummy, rdummy, rdummy, &
           quad_cvweight, quad_n, quad_nlx, quad_nly, quad_nlz, &
           quad_scvweight, quad_sn, quad_snlx, quad_snly, quad_sm, quad_smlx, quad_smly, &
           nwicel, d3 )
      ewrite(3,*)'quad_sn:', quad_sn
      ewrite(3,*)'quad_snlx:', quad_snlx
      ! Checking the output from shape_l_q_quad
      ewrite(3,*)'quad_cvweight:', ( quad_cvweight( xnod ), xnod = 1, quad_cv_ngi )
      ewrite(3,*)'quad_scvweight:', ( quad_scvweight( xnod ), xnod = 1, quad_cv_sngi )
      do xnod = 1, quad_cv_nloc
         ewrite(3,*)'quad_n:', xnod, ( quad_n( xnod, ele ), ele = 1, quad_cv_ngi )
         ewrite(3,*)'quad_nlx:', xnod, ( quad_nlx( xnod, ele ), ele = 1, quad_cv_ngi )
         ewrite(3,*)'quad_nly:', xnod, ( quad_nly( xnod, ele ), ele = 1, quad_cv_ngi )
         ewrite(3,*)''
      end do
      do xnod = 1, quad_cv_snloc
         ewrite(3,*)'quad_sn:', xnod, ( quad_sn( xnod, ele ), ele = 1, quad_cv_sngi )
         ewrite(3,*)'quad_snlx:', xnod, ( quad_snlx( xnod, ele ), ele = 1, quad_cv_sngi )
         ewrite(3,*)'quad_snly:', xnod, ( quad_snly( xnod, ele ), ele = 1, quad_cv_sngi )
         ewrite(3,*)''
      end do

      ! Now determine the basis functions at the 
      ! node pts
      tri_tet = .true.
      call shatri_hex( loc_coord_nod_l1, loc_coord_nod_l2, loc_coord_nod_l3, &
           loc_coord_nod_l4, rdummy, d3, &
           cv_nloc, x_nonods, &
           sn_i_xj, rdummy2, rdummy2, rdummy2, &
           tri_tet  )
      ! Chacking the output from shatri
      do xnod = 1, x_nonods
         ewrite(3,*)'loc_coord_nod_l1/4, sum:', loc_coord_nod_l1(xnod), &
              loc_coord_nod_l2(xnod), loc_coord_nod_l3(xnod), loc_coord_nod_l4(xnod)
      end do
      do xnod = 1, cv_nloc
         ewrite(3,*)'sn_i_xj:', xnod, ( sn_i_xj( xnod, ele ), ele = 1, x_nonods )
      end do
      !stop 2929

      Loop_Elements: do ele = 1, totele ! Calculate SDETWEI,RA,SNX,SNY,SNZ for element ELE
         ! What is the fem node belonging to this element (CV_ILOC):
         cv_iloc_belong = 0
         do quad_cv_iloc = 1, quad_cv_nloc
            xnod = x_ndgln( ( ele - 1 ) * quad_cv_nloc + quad_cv_iloc )
            if( fem_nod( xnod ) /= 0 ) cv_iloc_belong = xnod

            ewrite(3,*) ele, quad_cv_iloc, &
                 ( ele - 1 ) * quad_cv_nloc + quad_cv_iloc, &
                 xnod, fem_nod(xnod), &
                 x(xnod), y(xnod),z(xnod)
         end do

         Loop_SurfaceElements: do sele = 1, stotel ! Extract surface nodes
            iface = sele
            quad_cv_siloc = 0
            Conditional_Dimension: if( d1 ) then
               quad_cv_iloc = 1
               found = .true.
               if( found ) then
                  quad_cv_siloc = 1
                  x_sndgln( ( sele - 1 ) * quad_cv_snloc + quad_cv_siloc ) = &
                       x_ndgln( ( ele - 1 ) * quad_cv_nloc + quad_cv_iloc )
                  loc_2nd_lev( sele, quad_cv_siloc ) = quad_cv_iloc
               end if

            elseif( .not. d3 ) then ! 2D
               Loop_Poly2d_1_1: do ip = 1, npoly
                  Loop_Poly2d_1_2: do jp = 1, npoly
                     quad_cv_iloc = ( jp - 1 ) * npoly + ip 
                     found = .false.
                     Select Case( iface )
                     case( 1 )
                        if( ip == 1 ) found = .true.
                     case( 2 )  
                        if( ip == npoly ) found = .true.
                     case( 3 )
                        if( jp == 1 ) found = .true.
                     case( 4 )
                        if( jp == npoly ) found = .true.
                     case default; FLExit( "Wrong integer for IFACE" )
                     end Select
                     if( found ) then
                        quad_cv_siloc = quad_cv_siloc + 1
                        x_sndgln( ( sele - 1 ) * quad_cv_snloc + quad_cv_siloc ) = &
                             x_ndgln( ( ele - 1 ) * quad_cv_nloc + quad_cv_iloc )
                        loc_2nd_lev( sele, quad_cv_siloc ) = quad_cv_iloc
                        ewrite(3, *) 'ele,sele,quad_cv_snloc,quad_cv_nloc,quad_cv_siloc,x_sndgln:', &
                             ele,sele,quad_cv_snloc,quad_cv_nloc,quad_cv_siloc,&
                             x_sndgln( ( sele - 1 ) * quad_cv_snloc + quad_cv_siloc )
                     end if
                  end do Loop_Poly2d_1_2
               end do Loop_Poly2d_1_1

            else ! 3D
               Loop_Poly3d_2_1: do ip = 1, npoly
                  Loop_Poly3d_2_2: do jp = 1, npoly
                     Loop_Poly3d_2_3: do kp = 1, npoly
                        quad_cv_iloc = ( kp - 1 ) * npoly * npoly + &
                             ( jp - 1 ) * npoly + ip
                        found = .false.
                        Select Case( iface )
                        case( 1 )
                           if( ip == 1 ) found = .true.
                        case( 2 )  
                           if( ip == npoly ) found = .true.
                        case( 3 )
                           if( jp == 1 ) found = .true.
                        case( 4 )
                           if( jp == npoly ) found = .true.
                        case( 5 )
                           if( kp == 1 ) found = .true.
                        case( 6 )
                           if( kp == npoly ) found = .true.
                        case default; FLExit( "Wrong integer for IFACE" )
                        end Select
                        if( found ) then
                           quad_cv_siloc = quad_cv_siloc + 1
                           x_sndgln( ( sele - 1 ) * quad_cv_snloc + quad_cv_siloc ) = &
                                x_ndgln( ( ele - 1 ) * quad_cv_nloc + quad_cv_iloc )
                           loc_2nd_lev( sele, quad_cv_siloc ) = quad_cv_iloc
                           ewrite(3, *) 'ele,sele,quad_cv_snloc,quad_cv_nloc,quad_cv_siloc,x_sndgln:', &
                                ele,sele,quad_cv_snloc,quad_cv_nloc,quad_cv_siloc,&
                                x_sndgln( ( sele - 1 ) * quad_cv_snloc + quad_cv_siloc )
                        end if
                     end do Loop_Poly3d_2_3
                  end do Loop_Poly3d_2_2
               end do Loop_Poly3d_2_1
            endif Conditional_Dimension

            ! Now compute the determinant -- Quad_Sdetwei
            do quad_cv_siloc = 1, quad_cv_snloc
               xnod = x_sndgln( ( sele - 1 ) * quad_cv_snloc + quad_cv_siloc )
               xsl( quad_cv_siloc ) = x( xnod )
               ysl( quad_cv_siloc ) = y( xnod )
               if( d3 ) zsl( quad_cv_siloc ) = z( xnod )
               !ewrite(3,*)'x/y/zsl:', ele, sele, quad_cv_siloc, xsl( quad_cv_siloc ), &
               !     ysl( quad_cv_siloc ), zsl( quad_cv_siloc )
            end do

            call dgsdetnxloc2( quad_cv_snloc, quad_cv_sngi, &
                 xsl, ysl, zsl, &
                 quad_sn, quad_snlx, quad_snly, quad_scvweight, quad_sdetwei, sarea, &
                 ( ndim == 1 ), ( ndim == 3 ), ( ndim == -2 ), &
                 normxn, normyn, normzn, &
                 normx, normy, normz )
            !ewrite(3,*)'normx/y/z:', normx, normy, normz
            !ewrite(3,*)'quad_sdetwei:', ( quad_sdetwei( xnod ), xnod = 1, quad_cv_sngi )
            !do xnod = 1, quad_cv_sngi
            !  ewrite(3,*)'normx/y/zn:', xnod, normxn(xnod), normyn(xnod), normzn(xnod)
            !end do

            ! Take out quadrature points that are inside a CV: 
            ! NB: fem_nod(xnod) = 0 if not a local finite element node ELSE = local node no.
            ! This is determined by looking to see if any of the face nodes 
            ! have a local fem_nod. 
            found_fem_nod = .false.
            do quad_cv_sjloc = 1, quad_cv_snloc
               xnod = x_sndgln( ( sele - 1 ) * quad_cv_snloc + quad_cv_sjloc )
               if( fem_nod( xnod ) /= 0 ) found_fem_nod = .true.
               !ewrite(3,*)'xnod, fem_nod:', xnod, fem_nod( xnod ), found_fem_nod
            end do

            ! Keep record of the nodes that are neighbours to this surface quadrature pt
            do quad_cv_sgi = 1, quad_cv_sngi 
               quad_sgi = stotel * quad_cv_sngi * ( ele - 1 ) + quad_cv_sngi * &
                    ( sele - 1 ) + quad_cv_sgi
               next_to_cv_iloc_gi( quad_sgi ) = cv_iloc_belong
            end do

            ! Determine the derivative along the surface wrt to the basis functions SN. 
            Loop_QUAD_CV_SGI: do quad_cv_sgi = 1, quad_cv_sngi 
               cv_sgi = stotel * quad_cv_sngi * ( ele - 1 ) + quad_cv_sngi * &
                    ( sele - 1 ) + quad_cv_sgi
               if( found_fem_nod ) then
                  remove_ig_pt( cv_sgi ) = .true.
               else
                  remove_ig_pt( cv_sgi ) = .false.
               end if

               Loop_CV_ILOC: do cv_iloc = 1, cv_nloc
                  suf_quad_sn( cv_iloc, cv_sgi ) = 0.0
                  suf_quad_snlx( cv_iloc, cv_sgi ) = 0.0
                  suf_quad_snly( cv_iloc, cv_sgi ) = 0.0
                  Loop_QUAD_CV_SJLOC: do quad_cv_sjloc = 1, quad_cv_snloc
                     xnod = x_sndgln( ( sele - 1 ) * quad_cv_snloc + quad_cv_sjloc )
                     suf_quad_sn( cv_iloc, cv_sgi ) = suf_quad_sn( cv_iloc, cv_sgi ) + &
                          quad_sn( quad_cv_sjloc, quad_cv_sgi ) * sn_i_xj( cv_iloc, xnod )
                     suf_quad_snlx( cv_iloc, cv_sgi )= suf_quad_snlx( cv_iloc, cv_sgi ) + &
                          quad_snlx( quad_cv_sjloc, quad_cv_sgi ) * sn_i_xj( cv_iloc, xnod )
                     if( d3 ) suf_quad_snly( cv_iloc, cv_sgi ) =  &
                          suf_quad_snly( cv_iloc, cv_sgi ) + &
                          quad_snly( quad_cv_sjloc, quad_cv_sgi ) * sn_i_xj( cv_iloc, xnod )
                  end do Loop_QUAD_CV_SJLOC
               end do Loop_CV_ILOC

            end do Loop_QUAD_CV_SGI

            ewrite(3,*)'ele, sele, totele, quad_cv_sngi, quad_cv_snloc:', &
                 ele, sele, totele, quad_cv_sngi, quad_cv_snloc

            ! Determine the quadrature points and weights
            Loop_NGI: do quad_cv_sgi = 1, quad_cv_sngi 
               cv_sgi = ( ele - 1 ) * quad_cv_sngi * stotel + ( sele - 1 ) * &
                    quad_cv_sngi + quad_cv_sgi
               ! we do not multiply the weight by the determinant as we are operating in 
               ! local coordinates...
               ! gl_quad_scvweigh( cv_sgi ) = quad_sdetwei( quad_cv_sgi )
               gl_quad_scvweigh( cv_sgi ) = quad_scvweight( quad_cv_sgi )

               ewrite(3,*) ''

               xgi = 0. ; ygi = 0. ; zgi = 0.
               do quad_cv_siloc = 1, quad_cv_snloc
                  xnod = x_sndgln( ( sele - 1 ) * quad_cv_snloc + quad_cv_siloc )
                  xgi =  xgi + quad_sn( quad_cv_siloc, quad_cv_sgi ) * x( xnod )
                  ygi =  ygi + quad_sn( quad_cv_siloc, quad_cv_sgi ) * y( xnod )
                  if( d3 ) zgi =  zgi + quad_sn( quad_cv_siloc, quad_cv_sgi ) * z( xnod )

                  !ewrite(3,*)  quad_cv_siloc, x( xnod ), y( xnod ), z( xnod ), '::', &
                  !     xnod

               end do

               !ewrite(3,*) x_ndgln( ( ele - 1 ) * quad_cv_snloc + 1:&
               !     ( ele - 1 ) * quad_cv_snloc + quad_cv_nloc )
               !ewrite(3,*)  '------>> quad_cv_gi, xgi,  ygi, zgi', &
               !     quad_cv_sgi, cv_sgi, xgi,  ygi, zgi

               ! ewrite(3,*) lx(1:4)
               ! ewrite(3,*) ly(1:4)
               ! ewrite(3,*) lz(1:4)
               ! ewrite(3,*) volume_quad_map( 1, sum(lx(1:4))/4., sum(ly(1:4))/4., sum(lz(1:4))/4., lx, ly, lz )
               ! ewrite(3,*) volume_quad_map( 2, sum(lx(1:4))/4., sum(ly(1:4))/4., sum(lz(1:4))/4., lx, ly, lz )
               ! ewrite(3,*) volume_quad_map( 3, sum(lx(1:4))/4., sum(ly(1:4))/4., sum(lz(1:4))/4., lx, ly, lz )
               ! ewrite(3,*) volume_quad_map( 4, sum(lx(1:4))/4., sum(ly(1:4))/4., sum(lz(1:4))/4., lx, ly, lz )

               if( d3 ) then ! local coords for the tetrahedron
                  gl_quad_l1( cv_sgi ) = volume_quad_map( 1, xgi, ygi, zgi, lx(1:4), ly(1:4), lz(1:4) )
                  gl_quad_l2( cv_sgi ) = volume_quad_map( 2, xgi, ygi, zgi, lx(1:4), ly(1:4), lz(1:4) )
                  gl_quad_l3( cv_sgi ) = volume_quad_map( 3, xgi, ygi, zgi, lx(1:4), ly(1:4), lz(1:4) )
                  gl_quad_l4( cv_sgi ) = volume_quad_map( 4, xgi, ygi, zgi, lx(1:4), ly(1:4), lz(1:4) )
               else if(.not.d1) then ! 2D
                  gl_quad_l1( cv_sgi ) = area_quad_map( 1, xgi, ygi, lx, ly )
                  gl_quad_l2( cv_sgi ) = area_quad_map( 2, xgi, ygi, lx, ly )
                  gl_quad_l3( cv_sgi ) = area_quad_map( 3, xgi, ygi, lx, ly )
               else 
                  gl_quad_l1( cv_sgi ) = 1.0
               end if
            end do Loop_NGI

         end do Loop_SurfaceElements

      end do Loop_Elements

      ! Now determine the basis functions and derivatives at the 
      ! quadrature pts quad_L1, quad_L2, quad_L3, quad_L4, etc
      tri_tet = .true.
      call shatri_hex( gl_quad_l1, gl_quad_l2, gl_quad_l3, gl_quad_l4, rdummy, d3, &
           cv_nloc, stotel * quad_cv_sngi * totele, &
           gl_quad_sn, gl_quad_snlx, gl_quad_snly, gl_quad_snlz, &
           tri_tet )
      ewrite(3,*)'gl_quad_l1:', ( gl_quad_l1( cv_sgi ), cv_sgi = 1, stotel * quad_cv_sngi * totele )
      ewrite(3,*)'gl_quad_l2:', ( gl_quad_l2( cv_sgi ), cv_sgi = 1, stotel * quad_cv_sngi * totele )
      ewrite(3,*)'gl_quad_l3:', ( gl_quad_l3( cv_sgi ), cv_sgi = 1, stotel * quad_cv_sngi * totele )
      ewrite(3,*)'gl_quad_l4:', ( gl_quad_l4( cv_sgi ), cv_sgi = 1, stotel * quad_cv_sngi * totele )
      ewrite(3,*)' '

      do xnod = 1, cv_nloc
         ewrite(3,*)'gl_quad_sn:', xnod, ( gl_quad_sn( xnod, ele ), &
              ele = 1, stotel * quad_cv_sngi * totele )
         ewrite(3,*)'gl_quad_snlx:', xnod, ( gl_quad_snlx( xnod, ele ), &
              ele = 1, stotel * quad_cv_sngi * totele )
         ewrite(3,*)'gl_quad_snly:', xnod, ( gl_quad_snly( xnod, ele ), &
              ele = 1, stotel * quad_cv_sngi * totele )
      end do
      ! 
      ! Find shared quadrature points to see what is on the other side 
      do cv_sgi = 1, stotel * quad_cv_sngi * totele
         found = .false.
         do cv_sgj = 1, stotel * quad_cv_sngi * totele
            if( cv_sgi /= cv_sgj ) then
               d2_quad = abs( gl_quad_l1( cv_sgi ) - gl_quad_l1( cv_sgj ) ) &
                    + abs( gl_quad_l2( cv_sgi ) - gl_quad_l2( cv_sgj ) ) &
                    + abs( gl_quad_l3( cv_sgi ) - gl_quad_l3( cv_sgj ) )
               if( d3 ) d2_quad = d2_quad  &
                    + abs( gl_quad_l4( cv_sgi ) - gl_quad_l4( cv_sgj ) )
               if( d2_quad < 1.0e-5 ) then ! same pt
                  found = .true. 
                  ! ewrite(3,*) cv_sgi, cv_sgj, next_to_cv_iloc_gi( cv_sgi ), &
                  ! next_to_cv_iloc_gi( cv_sgj ), fem_nod( next_to_cv_iloc_gi( cv_sgj ) )
                  ! next_to_cv_iloc_gi2( cv_sgi ) contains the node on the other side.
                  next_to_cv_iloc_gi2( cv_sgi ) = next_to_cv_iloc_gi( cv_sgj )
               endif
            endif
         end do
      end do

      ! adjust remove_ig_pt to always set to false if on the boundary of the element: 
      do cv_sgi = 1, stotel * quad_cv_sngi * totele
         if( next_to_cv_iloc_gi2( cv_sgi ) == 0 ) then
            remove_ig_pt( cv_sgi )= .false.
         endif
      end do

      ewrite(3,*)'stotel, quad_cv_sngi, totele:',stotel, quad_cv_sngi, totele
      ewrite(3,*)'stotel * quad_cv_sngi * totele:',stotel * quad_cv_sngi * totele

      allocate( sn_2( cv_nloc, stotel * quad_cv_sngi * totele ) ) ; sn_2 = 0.
      allocate( suf_snlx_2( cv_nloc, stotel * quad_cv_sngi * totele ) ) ; suf_snlx_2 = 0.
      allocate( suf_snly_2( cv_nloc, stotel * quad_cv_sngi * totele ) ) ; suf_snly_2 = 0.
      allocate( snlx_2( cv_nloc, stotel * quad_cv_sngi * totele ) ) ; snlx_2 = 0.
      allocate( snly_2( cv_nloc, stotel * quad_cv_sngi * totele ) ) ; snly_2 = 0.
      allocate( snlz_2( cv_nloc, stotel * quad_cv_sngi * totele ) ) ; snlz_2 = 0.
      allocate( scvweigh_2( stotel * quad_cv_sngi * totele ) )
      allocate( cv_neiloc_cells_2( cv_nloc_cells, stotel * quad_cv_sngi * totele ) )
      allocate( cv_neiloc_cells_3( cv_nloc_cells, stotel * quad_cv_sngi * totele ) )
      allocate( l1_2( stotel * quad_cv_sngi * totele ) ) ; l1_2 = 0.
      allocate( l2_2( stotel * quad_cv_sngi * totele ) ) ; l2_2 = 0.
      allocate( l3_2( stotel * quad_cv_sngi * totele ) ) ; l3_2 = 0.
      allocate( l4_2( stotel * quad_cv_sngi * totele ) ) ; l4_2 = 0.

      cv_neiloc_cells_3 = 0
      ! Define cv_neiloc_cells_3 from next_to_cv_iloc_gi2; next_to_cv_iloc_gi: 
      Loop_Neiloc3: do cv_sgi = 1, stotel * quad_cv_sngi * totele
         Conditional_Neiloc3: if( .not. remove_ig_pt( cv_sgi ) ) then
            xnodi = next_to_cv_iloc_gi( cv_sgi )
            xnodj = next_to_cv_iloc_gi2( cv_sgi )
            if ( xnodi == 0 ) then
               nodi = 0
            else
               nodi = fem_nod( xnodi )
            end if
            if ( xnodj == 0 ) then
               nodj = 0
            else
               nodj = fem_nod( xnodj )
            end if

            if ( ( nodi /= 0 ) .and. ( nodj /= 0 ) ) then
               cv_neiloc_cells_3( nodi, cv_sgi ) = nodj
               ! If there is no node on the other side of the quadrature pt then set -1
            elseif ( nodi /= 0 ) then
               cv_neiloc_cells_3( nodi, cv_sgi ) = -1
            elseif ( nodj /= 0 ) then
               cv_neiloc_cells_3( nodj, cv_sgi ) = -1
            endif

         endif Conditional_Neiloc3
      end do Loop_Neiloc3

      ! Take out quadrature points that are inside a CV: 
      ! NB fem_nod(xnod) = 0 if not a local finite element node ELSE = local node no.
      ! This is determined by looking to see if any of the face nodes 
      ! have a local fem node. 

      ewrite(3,*) 'suf_quad_snlx:', suf_quad_snlx
      ewrite(3,*) 'suf_quad_snly:', suf_quad_snly
      ewrite(3,*) 'remove_ig_pt:', remove_ig_pt
      ewrite(3,*) 'fem_nod:', fem_nod
      ewrite(3,*) 'stotel, quad_cv_sngi, totele:', stotel, &
           quad_cv_sngi, totele, stotel * quad_cv_sngi * totele

      cv_sgk = 0 ; cv_neiloc_cells_2 = 0
      do cv_sgi = 1, stotel * quad_cv_sngi * totele
         if( .not. remove_ig_pt( cv_sgi ) ) then
            cv_sgk = cv_sgk + 1
            ! this is also gl_quad_sn( :, cv_sgi ):
            sn_2( :, cv_sgk ) = suf_quad_sn( :, cv_sgi )
            suf_snlx_2( :, cv_sgk ) = suf_quad_snlx( :, cv_sgi )
            suf_snly_2( :, cv_sgk ) = suf_quad_snly( :, cv_sgi )
            snlx_2( :, cv_sgk ) = gl_quad_snlx( :, cv_sgi )
            snly_2( :, cv_sgk ) = gl_quad_snly( :, cv_sgi )
            snlz_2( :, cv_sgk ) = gl_quad_snlz( :, cv_sgi )
            scvweigh_2( cv_sgk ) = gl_quad_scvweigh( cv_sgi )
            cv_neiloc_cells_2( :, cv_sgk ) = cv_neiloc_cells_3( :, cv_sgi )
            l1_2( cv_sgk ) = gl_quad_l1( cv_sgi )
            l2_2( cv_sgk ) = gl_quad_l2( cv_sgi )
            l3_2( cv_sgk ) = gl_quad_l3( cv_sgi )
            if( d3 ) l4_2( cv_sgk ) = gl_quad_l4( cv_sgi )
            ewrite(3,*)'cv_sgk, l1/2/3/4:', cv_sgk, l1_2( cv_sgk ), &
                 l2_2( cv_sgk ), l3_2( cv_sgk ), l4_2( cv_sgk )
         end if
      end do

      ewrite(3,*) 'suf_snlx_2:', suf_snlx_2

      cv_sngi_2 = cv_sgk
      ewrite(3,*) 'quad_cv_ngi, cv_sngi_2:', quad_cv_sngi, cv_sngi_2

      ewrite(3,*) 'd3, stotel * quad_cv_sngi * totele, cv_sngi_2:', &
           d3, stotel * quad_cv_sngi * totele, cv_sngi_2
!       stop 83

      ! Take out repetition of quadrature points
      cv_sgk = 0 ; cv_neiloc_cells = 0 ; l1 = 0. ; l2 = 0. ; l3 = 0. ; l4=0.
      Loop_CV_SGI: do cv_sgi = 1, cv_sngi_2
         found = .false.
         Loop_CV_SGJ: do cv_sgj = 1, cv_sgi - 1
            d2_quad = abs( l1_2( cv_sgi ) - l1_2( cv_sgj ) ) &
                 + abs( l2_2( cv_sgi ) - l2_2( cv_sgj ) ) &
                 + abs( l3_2( cv_sgi ) - l3_2( cv_sgj ) ) 
            if( d3 ) d2_quad = d2_quad &
                 + abs( l4_2( cv_sgi ) - l4_2( cv_sgj ) ) 
            if( d2_quad < 1.0e-5 ) found = .true. ! same pt -- remove repetition
         end do Loop_CV_SGJ

         if( .not. found ) then
            cv_sgk = cv_sgk + 1

            ewrite(3,*)cv_sgi, cv_sgk,scvngi,cv_sngi_2, sn_2( 1:3, cv_sgi )

            sn( :, cv_sgk ) = sn_2( :, cv_sgi )
            if( ndim.ge.2 ) sufnlx( :, cv_sgk ) = suf_snlx_2( :, cv_sgi )
            if( ndim.ge.3 ) sufnly( :, cv_sgk ) = suf_snly_2( :, cv_sgi )
            snlx( :, cv_sgk ) = snlx_2( :, cv_sgi )
            if( ndim.ge.2 ) snly( :, cv_sgk ) = snly_2( :, cv_sgi )
            if( ndim.ge.3 ) snlz( :, cv_sgk ) = snlz_2( :, cv_sgi )
            scvweigh( cv_sgk ) = scvweigh_2( cv_sgi )
            l1( cv_sgk ) = l1_2( cv_sgi )
            l2( cv_sgk ) = l2_2( cv_sgi )
            l3( cv_sgk ) = l3_2( cv_sgi )
            if( d3 ) l4( cv_sgk ) = l4_2( cv_sgi )
            cv_neiloc_cells( :, cv_sgk ) = cv_neiloc_cells_2( :, cv_sgi )    
         endif
      end do Loop_CV_SGI

      ewrite(3,*) 'sufnlx( : , 1 : scvngi ):', sufnlx( :, 1 : scvngi )

      if( scvngi /= cv_sgk ) then
         ewrite(3,*) 'scvngi, cv_sgk:', scvngi, cv_sgk
         FLExit( "SCVNGI /= CV_SGK " )
      endif

      ewrite(3,*)'l1:', l1
      ewrite(3,*)'l2:', l2
      ewrite(3,*)'l3:', l3
      ewrite(3,*)'l4:', l4
      ewrite(3,*)'sn:', sn
      ewrite(3,*)'snlx:', snlx
      ewrite(3,*)'snly:', snly
      ewrite(3,*)'snlz:', snlz

      ! Remapping over the common quadrature points across CV
      do cv_iloc_cells = 1, cv_nloc_cells
         do cv_sgi = 1, scvngi
            if ( cv_neiloc_cells( cv_iloc_cells, cv_sgi ) > 0 ) then
               cv_jloc_cells = cv_neiloc_cells( cv_iloc_cells, cv_sgi )
               cv_neiloc_cells( cv_jloc_cells, cv_sgi ) = cv_iloc_cells
            end if
         end do
      end do

      icount=0
      cvfem_neiloc( :, : ) = 0
      ! calculate cvfem_neiloc from local coords l1-4:
      Loop_SurfaceQuadrature: do cv_sgi = 1, scvngi
         zer_l1 = ( abs( l1( cv_sgi ) ) < 1.0e-4 )
         zer_l2 = ( abs( l2( cv_sgi ) ) < 1.0e-4 )
         zer_l3 = ( abs( l3( cv_sgi ) ) < 1.0e-4 )
         if ( d3 ) zer_l4 = ( abs( l4( cv_sgi )) < 1.0e-4 )

!         if(zer_l4 ) icount=icount+1
         if( zer_l1 .or. zer_l2 .or. zer_l3 .or. zer_l4 )  icount=icount+1

         ewrite(3,*) 'gi, l1/2/3/4:', cv_sgi, abs( l1( cv_sgi ) ), abs( l2( cv_sgi ) ), &
              abs( l3( cv_sgi ) ), abs( l4( cv_sgi ) )

         Conditional_Neiloc: if ( d3 ) then
            ewrite(3,*)'cv_sgi, zer_l1/2/3/4:', cv_sgi, (zer_l1.or.zer_l2.or.zer_l3.or.zer_l4)
            if ( zer_l1 .or. zer_l2 .or. zer_l3 .or. zer_l4 ) then
               ! on the surface of the element: 
               if(cv_nloc==10) then
                 if(zer_l1) then
! face 1:
                    cvfem_neiloc( 6, cv_sgi )=-1
                    cvfem_neiloc( 5, cv_sgi )=-1
                    cvfem_neiloc( 3, cv_sgi )=-1
                    cvfem_neiloc( 9, cv_sgi )=-1
                    cvfem_neiloc( 8, cv_sgi )=-1
                    cvfem_neiloc( 10, cv_sgi )=-1
                 else if(zer_l2) then
! face 2:
                    cvfem_neiloc( 1, cv_sgi )=-1
                    cvfem_neiloc( 4, cv_sgi )=-1
                    cvfem_neiloc( 6, cv_sgi )=-1
                    cvfem_neiloc( 7, cv_sgi )=-1
                    cvfem_neiloc( 9, cv_sgi )=-1
                    cvfem_neiloc( 10, cv_sgi )=-1
                 else if(zer_l3) then
! face 3:
                    cvfem_neiloc( 1, cv_sgi )=-1
                    cvfem_neiloc( 2, cv_sgi )=-1
                    cvfem_neiloc( 3, cv_sgi )=-1
                    cvfem_neiloc( 7, cv_sgi )=-1
                    cvfem_neiloc( 8, cv_sgi )=-1
                    cvfem_neiloc( 10, cv_sgi )=-1
                 else if(zer_l4) then
! face 4:
                    cvfem_neiloc( 1, cv_sgi )=-1
                    cvfem_neiloc( 2, cv_sgi )=-1
                    cvfem_neiloc( 3, cv_sgi )=-1
                    cvfem_neiloc( 4, cv_sgi )=-1
                    cvfem_neiloc( 5, cv_sgi )=-1
                    cvfem_neiloc( 6, cv_sgi )=-1
                 endif
               else
               do cv_iloc = 1, cv_nloc
                  ewrite(3,*)'iloc, sn',cv_iloc, abs(sn( cv_iloc, cv_sgi )), &
                       abs ( sn( cv_iloc, cv_sgi )) > 1.e-4
                  if ( abs ( sn( cv_iloc, cv_sgi )) > 1.e-4 ) &
                       cvfem_neiloc( cv_iloc, cv_sgi ) = -1

               end do
               endif
            endif
         else
            ewrite(3,*)'cv_sgi, zer_l1/2/3:', cv_sgi, (zer_l1.or.zer_l2.or.zer_l3)
            if ( zer_l1 .or. zer_l2 .or. zer_l3 ) then
               ! on the surface of the element: 
               if(cv_nloc==6) then ! quadratic triangle
                 if(zer_l1) then
! face 1:
                    cvfem_neiloc( 3, cv_sgi )=-1
                    cvfem_neiloc( 5, cv_sgi )=-1
                    cvfem_neiloc( 6, cv_sgi )=-1
                 else if(zer_l2) then
! face 2:
                    cvfem_neiloc( 1, cv_sgi )=-1
                    cvfem_neiloc( 4, cv_sgi )=-1
                    cvfem_neiloc( 6, cv_sgi )=-1
                 else if(zer_l3) then
! face 3:
                    cvfem_neiloc( 1, cv_sgi )=-1
                    cvfem_neiloc( 2, cv_sgi )=-1
                    cvfem_neiloc( 3, cv_sgi )=-1
                 endif
               else
                  do cv_iloc = 1, cv_nloc
                     ewrite(3,*)'iloc, sn',cv_iloc, abs(sn( cv_iloc, cv_sgi )), &
                          abs ( sn( cv_iloc, cv_sgi )) > 1.e-4
                     if ( abs( sn( cv_iloc, cv_sgi )) > 1.e-4 ) &
                          cvfem_neiloc( cv_iloc, cv_sgi ) = -1
                     ewrite(3,*)'iloc, sn',cv_iloc, abs(sn( cv_iloc, cv_sgi )), abs ( sn( cv_iloc, cv_sgi )) > 1.e-4
                  end do
               endif
            endif
         endif Conditional_Neiloc

      end do Loop_SurfaceQuadrature

      ewrite(3,*)'icount::',icount

      do quad_cv_siloc = 1, cv_nloc
         do cv_sgi = 1, scvngi
            ewrite(3,*)'iloc, gi, cvfem_neiloc::', &
                 quad_cv_siloc, cv_sgi, cvfem_neiloc( quad_cv_siloc, cv_sgi )
         end do
      end do

      do cv_sgi = 1, scvngi
         ewrite(3,*)'cv_sgi, cvfem_neiloc::', &
              cv_sgi, cvfem_neiloc( :, cv_sgi )
      end do

      ! Calculate cvfem_neiloc from local coordinates l1-4: (hard-wired for linear traingles)
      if( .false. ) then
         cvfem_neiloc = 0

         cv_iloc = 1 ! CV 1
         cvfem_neiloc( cv_iloc, 1 ) = -1
         cvfem_neiloc( cv_iloc, 5 ) = -1
         cvfem_neiloc( cv_iloc, 3 ) = -1
         cvfem_neiloc( cv_iloc, 8 ) = -1

         cv_iloc = 2 ! CV 2
         cvfem_neiloc( cv_iloc, 1 ) = -1
         cvfem_neiloc( cv_iloc, 5 ) = -1
         cvfem_neiloc( cv_iloc, 6 ) = -1
         cvfem_neiloc( cv_iloc, 9 ) = -1

         cv_iloc = 3 ! CV 3
         cvfem_neiloc( cv_iloc, 6 ) = -1
         cvfem_neiloc( cv_iloc, 9 ) = -1
         cvfem_neiloc( cv_iloc, 8 ) = -1
         cvfem_neiloc( cv_iloc, 3 ) = -1
      end if

      ewrite(3,*) 'scvweigh:',scvweigh
      ewrite(3,*) ' '
      do cv_sgi = 1, scvngi
         ewrite(3,*) 'cv_sgi=',cv_sgi
         ewrite(3,*) 'sn( :, cv_sgi )    =', sn( :, cv_sgi )
      end do
      ewrite(3,*) ' '
      do cv_sgi = 1, scvngi
         ewrite(3,*) 'cv_sgi=',cv_sgi
         ewrite(3,*) 'sufnlx( :, cv_sgi )=',sufnlx( :, cv_sgi )
      end do
      ewrite(3,*) ' '
      do cv_sgi = 1, scvngi
         ewrite(3,*) 'cv_sgi=',cv_sgi
         ewrite(3,*) 'snlx( :, cv_sgi ):', snlx( :, cv_sgi )
      end do
      ewrite(3,*) ' '
      do cv_sgi = 1, scvngi
         ewrite(3,*) 'cv_sgi=',cv_sgi
         ewrite(3,*) 'snly( :, cv_sgi ):', snly( :, cv_sgi )
      end do
!      stop 281

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
      deallocate( rdummy2 )
      deallocate( quad_scvweight )
      deallocate( quad_sdetwei )
      deallocate( sra )
      deallocate( quad_sn )
      deallocate( quad_snlx )
      deallocate( quad_snly )
      deallocate( quad_snx )
      deallocate( quad_sny )
      deallocate( quad_sm )
      deallocate( quad_smlx )
      deallocate( quad_smly )
      deallocate( loc_coord_nod_l1 )
      deallocate( loc_coord_nod_l2 )
      deallocate( loc_coord_nod_l3 )
      deallocate( loc_coord_nod_l4 )
      deallocate( sn_i_xj )
      deallocate( suf_quad_sn )
      deallocate( suf_quad_snlx )
      deallocate( suf_quad_snly )
      deallocate( gl_quad_l1 )
      deallocate( gl_quad_l2 )
      deallocate( gl_quad_l3 )
      deallocate( gl_quad_l4 )
      deallocate( gl_quad_sn )
      deallocate( gl_quad_snlx )
      deallocate( gl_quad_snly )
      deallocate( gl_quad_snlz )
      deallocate( gl_quad_scvweigh )
      deallocate( remove_ig_pt )
      deallocate( next_to_cv_iloc_gi )
      deallocate( next_to_cv_iloc_gi2 )

      deallocate( normxn ) 
      deallocate( normyn ) 
      deallocate( normzn ) 
      deallocate( xsl )  
      deallocate( ysl )
      deallocate( zsl )   
      deallocate( suf_snlx_2 )
      deallocate( suf_snly_2 )
      deallocate( sn_2 )
      deallocate( snlx_2 )
      deallocate( snly_2 )
      deallocate( snlz_2 )
      deallocate( scvweigh_2 )
      deallocate( cv_neiloc_cells_2 )
      deallocate( cv_neiloc_cells_3 )
      deallocate( l1_2 )
      deallocate( l2_2 )
      deallocate( l3_2 )
      deallocate( l4_2 )
      deallocate( l1 )
      deallocate( l2 )
      deallocate( l3 )
      deallocate( l4 )

      return
    end subroutine Compute_SurfaceShapeFunctions_Triangle_Tetrahedron


    subroutine dummy_tri_tet( d1, d3, quad_cv_ngi, quad_cv_nloc, &
         dummy_sngi, dummy_snloc, nwicel, & 
         cv_nloc, cv_sngi, totele, quad_u_loc_dummy, &
         mloc, dummy_smloc, lowqua, npoly, npoly_ngi )
      implicit none
      ! Compute some local variables for suf_shape_tri_tet
      logical, intent( in ) :: d1, d3
      integer, intent( inout ) :: quad_cv_ngi
      integer, intent( in ) :: quad_cv_nloc
      integer, intent( inout ) :: dummy_sngi, dummy_snloc, nwicel
      integer, intent( in ) :: cv_nloc, cv_sngi, totele 
      integer, intent( inout ) :: quad_u_loc_dummy, mloc, dummy_smloc
      logical, intent( inout ) :: lowqua
      integer, intent( inout ) :: npoly, npoly_ngi

      ewrite(3,*)'In dummy_tri_tet'

      quad_u_loc_dummy = 1
      mloc = 1
      dummy_smloc = 1
      lowqua = .false.

      Conditional_Dimensionality1: if( d3 ) then
         quad_cv_ngi = 8
         dummy_sngi = 4
         dummy_snloc = 4
         nwicel = 1
         if( cv_nloc == 10 ) then ! Quadratic hexs
            quad_cv_ngi = 27
            nwicel = 3
            dummy_sngi = 9
            dummy_snloc = 9
            mloc = 8
         end if
      else
         quad_cv_ngi = 4
         dummy_sngi = 2
         dummy_snloc = 2
         nwicel = 1
         if ( cv_nloc == 6 ) then ! Quadratic quads
            quad_cv_ngi = 9
            dummy_sngi = 3
            dummy_snloc = 3
            nwicel = 3
            mloc = 4
         end if
      end if Conditional_Dimensionality1

      ! Consistency check:
      if( .false. )then
         if( cv_sngi /= totele * quad_cv_ngi ) &
              FLExit( "Wrong number for CV_NGI" )
      end if

      npoly = 2 ! Linear default
      if( d1 ) then
         npoly = 1
      elseif ( .not. d3 ) then ! 2D
         if( quad_cv_nloc == 9 ) npoly = 3 ! Quadratic
      else ! 3D
         if( quad_cv_nloc == 27 ) npoly = 3 ! Quadratic
      endif
      if( quad_cv_nloc == 1 ) npoly = 1

      ! for the quadrature pts:
      npoly_ngi = npoly
      if (( .not. d1 ).and. (.not. d3) ) then ! 2D
         if( quad_cv_nloc == 4 ) then
            if(cv_sngi == 9*3 ) npoly_ngi = 3 ! 3 pt
            if(cv_sngi == 9*2 ) npoly_ngi = 2 ! 2 pt
            if(cv_sngi == 9*1 ) npoly_ngi = 1 ! 1 pt
         else if( quad_cv_nloc == 9 ) then
            if(cv_sngi == 24*3 ) npoly_ngi = 3 ! 3 pt
            if(cv_sngi == 24*2 ) npoly_ngi = 2 ! 2 pt
            if(cv_sngi == 24*1 ) npoly_ngi = 1 ! 1 pt
         endif 
      else ! 3D
         if( quad_cv_nloc == 8 ) then
            npoly_ngi = 3 ! Linear
            if(cv_sngi == 18*9 ) npoly_ngi = 3 ! 3*3 pt
            if(cv_sngi == 18*4 ) npoly_ngi = 2 ! 2*2 pt
            if(cv_sngi == 18*1 ) npoly_ngi = 1 ! 1 pt
         else if( quad_cv_nloc == 27 ) then
            npoly_ngi = 3 ! Quadratic
            if(cv_sngi == 96*9 ) npoly_ngi = 3 ! 3*3 pt
            if(cv_sngi == 96*4 ) npoly_ngi = 2 ! 2*2 pt
            if(cv_sngi == 96*1 ) npoly_ngi = 1 ! 1 pt
         endif 
      endif
      if( quad_cv_nloc == 1 ) npoly_ngi = 1

      if( d1 .and. ( npoly /= 1 )) then
         ewrite(3,*) 'npoly is wrong -- it should be 1 instead of ', npoly
         FLExit( "Wrong number for NPOLY" )
      endif

      return
    end subroutine dummy_tri_tet



    subroutine shape_tri_tet( cv_ele_type_cells, cv_nloc_cells, &
         cv_ele_type, ndim, totele, cv_nloc, cv_ngi, x_nonods, &
         quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, &
         n, nlx, nly, nlz, cvweigh )

      ! Determine the volume shape functions n, nlx, nly, nlz 
      ! and weights cvweigh for the cv_nloc_cells CV cells. 

      use shape_functions
      implicit none
      integer, intent( in ) :: cv_ele_type, ndim, totele, cv_nloc, cv_ngi, &
           x_nonods, quad_cv_nloc, cv_ele_type_cells, cv_nloc_cells
      integer, dimension( totele * quad_cv_nloc ), intent( in ) :: x_ndgln
      real, dimension( x_nonods ), intent( in ) :: x, y, z
      real, dimension( quad_cv_nloc ), intent( in ) :: lx, ly, lz
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: n, nlx, nly, nlz
      real, dimension( cv_ngi ), intent( inout ) ::  cvweigh
      ! Local variables
      logical :: d1, dcyl, d3, lowqua
      integer :: ele, quad_cv_ngi, quad_cv_gi, cv_gi, nwicel, xnod, quad_cv_iloc, &
           quad_u_loc_dummy, mloc, dummy_sngi, dummy_snloc, dummy_smloc, &
           id,jd,kd
      real :: xgi, ygi, zgi, volume, rsum, xstar, xfini, ystar, yfini 
      real :: zstar, zfini,xcoord,ycoord,zcoord
      real :: rsum1,rsum2,rsum3,rsum4,rsum5,rsum6,rsum7,rsum8
      integer :: nod1,nod2,nod3, nod4,nod5,nod6,nod7,nod8
      real, dimension( : ), allocatable :: quad_l1, quad_l2, quad_l3, quad_l4, &
           quad_cvweight, detwei, ra, rdummy, &
           x_temp, y_temp, z_temp
      real, dimension( :, : ), allocatable :: quad_n, quad_nlx, quad_nly, quad_nlz, &
           quad_nx, quad_ny, quad_nz
      integer, dimension( : ), allocatable :: x_ndgln_temp, nod_pt


      ewrite(3,*)'In shape_tri_tet'

      d3 = ( ndim == 3 )
      d1 = ( ndim == 1 )
      dcyl = .false. 

      Conditional_Dimensionality1: if( d3 ) then
         quad_cv_ngi = 8
         dummy_sngi = 4
         dummy_snloc = 4
         nwicel = 1 ! was 3
         quad_u_loc_dummy = 1
         mloc = 1
         dummy_smloc = 1
         lowqua = .false.
         if( cv_nloc_cells == 4 ) then ! Linear hexs
            quad_cv_ngi =  8  !27
            if(cv_ngi==4) quad_cv_ngi = 1 ! one pt quadrature...
            if(cv_ngi==4*8) quad_cv_ngi = 8 ! 2 pt quadrature...
            if(cv_ngi==4*27) quad_cv_ngi = 27 ! 3 pt quadrature...
         endif
         if( cv_nloc_cells == 10 ) then ! Quadratic hexs
            quad_cv_ngi =  8  !27
            if(cv_ngi==32) quad_cv_ngi = 1 ! one pt quadrature...
            if(cv_ngi==32*8) quad_cv_ngi = 8 ! 2 pt quadrature...
            if(cv_ngi==32*27) quad_cv_ngi = 27 ! 3 pt quadrature...
            nwicel = 3
            dummy_sngi = 9
            dummy_snloc = 9
            mloc = 8
         end if
      else
         quad_cv_ngi = 4
         dummy_sngi = 2
         dummy_snloc = 2
         nwicel = 1
         quad_u_loc_dummy = 1
         mloc = 1
         dummy_smloc = 1
         lowqua = .false.
         if ( cv_nloc_cells == 3 ) then !  linear
            !!quad_cv_ngi = 9
            if(cv_ngi==3) quad_cv_ngi = 1 ! one pt quadrature...
            if(cv_ngi==3*4) quad_cv_ngi = 4 ! 2 pt quadrature...
            if(cv_ngi==3*9) quad_cv_ngi = 9 ! 3 pt quadrature...
         endif

         if ( cv_nloc_cells == 6 ) then ! Quadratic quads
            !!quad_cv_ngi = 9
            if(cv_ngi==12) quad_cv_ngi = 1 ! one pt quadrature...
            if(cv_ngi==12*4) quad_cv_ngi = 4 ! 2 pt quadrature...
            if(cv_ngi==12*9) quad_cv_ngi = 9  ! 3 pt quadrature...
            dummy_sngi = 3
            dummy_snloc = 3
            nwicel = 3
            mloc = 4
         end if
      end if Conditional_Dimensionality1

      if ( cv_ngi == 9 ) quad_cv_ngi = 3

      ! Consistency check:
      if( cv_ngi /= totele * quad_cv_ngi ) then
         ewrite(3,*)'cv_ngi, totele, quad_cv_ngi,totele * quad_cv_ngi:', &
                  cv_ngi, totele, quad_cv_ngi,totele * quad_cv_ngi
         FLExit( "Wrong number for CV_NGI" )
      end if

      ! Allocating memory
      allocate( quad_l1( cv_ngi ) ) ; quad_l1 = 0.
      allocate( quad_l2( cv_ngi ) ) ; quad_l2 = 0.
      allocate( quad_l3( cv_ngi ) ) ; quad_l3 = 0.
      allocate( quad_l4( cv_ngi ) ) ; quad_l4 = 0.
      allocate( quad_cvweight( quad_cv_ngi ) ) ; quad_cvweight = 0.
      allocate( detwei( quad_cv_ngi ) ) ; detwei = 0.
      allocate( ra( quad_cv_ngi ) ) ; ra = 0.
      allocate( quad_n( quad_cv_nloc, quad_cv_ngi ) ) ;  quad_n = 0.
      allocate( quad_nlx( quad_cv_nloc, quad_cv_ngi ) ) ; quad_nlx = 0.
      allocate( quad_nly( quad_cv_nloc, quad_cv_ngi ) ) ; quad_nly = 0.
      allocate( quad_nlz( quad_cv_nloc, quad_cv_ngi ) ) ; quad_nlz = 0.
      allocate( quad_nx( quad_cv_nloc, quad_cv_ngi ) ) ; quad_nx = 0.
      allocate( quad_ny( quad_cv_nloc, quad_cv_ngi ) ) ; quad_ny = 0.
      allocate( quad_nz( quad_cv_nloc, quad_cv_ngi ) ) ; quad_nz = 0.
      allocate( rdummy( 10000 ) ) ;  rdummy = 0.

      ewrite(3,*)'Just b4 shape_l_q_quad from shape_tri_tet'
      ewrite(3,*) 'totele, x_nonods, cv_nloc_cells, cv_nloc, cv_ngi:', &
           totele, x_nonods, cv_nloc_cells, cv_nloc, cv_ngi
      ewrite(3,*) 'quad_cv_nloc, quad_cv_ngi:', &
           quad_cv_nloc, quad_cv_ngi
      ewrite(3,*)'quad_cv_ngi, lowqua, mloc, nwicel:', &
           quad_cv_ngi, lowqua, mloc, nwicel
!       stop 3821

      ! Now we need to compute QUAD_NLX/Y/Z - get the hex or quad
      ! shape functions quad_n etc.
      call shape_l_q_quad( lowqua, quad_cv_ngi, quad_cv_nloc, mloc, &
           dummy_sngi, dummy_snloc, dummy_smloc, rdummy, rdummy, rdummy, rdummy, &
           quad_cvweight, quad_n, quad_nlx, quad_nly, quad_nlz, &
           rdummy, rdummy, rdummy, rdummy, rdummy, rdummy, rdummy, &
           nwicel, d3 )   

      ewrite(3,*)'quad_n:',quad_n
      ewrite(3,*)'quad_nlx:',quad_nlx
      ewrite(3,*)'quad_nly:',quad_nly
      ewrite(3,*)'quad_nlz:',quad_nlz

      !        stop 2921

      if(.false.) then
         xnod = x_ndgln( 1 )
         ewrite(3,*)'node1, x,y,z:',xnod,x(xnod),y(xnod),z(xnod)
         xnod = x_ndgln( 3 )
         ewrite(3,*)'node3, x,y,z:',xnod,x(xnod),y(xnod),z(xnod)
         xnod = x_ndgln( 7 )
         ewrite(3,*)'node7, x,y,z:',xnod,x(xnod),y(xnod),z(xnod)
         xnod = x_ndgln( 9 )
         ewrite(3,*)'node9, x,y,z:',xnod,x(xnod),y(xnod),z(xnod)

         xnod = x_ndgln( 19 )
         ewrite(3,*)'node19, x,y,z:',xnod,x(xnod),y(xnod),z(xnod)
         xnod = x_ndgln( 21 )
         ewrite(3,*)'node21, x,y,z:',xnod,x(xnod),y(xnod),z(xnod)
         xnod = x_ndgln( 25 )
         ewrite(3,*)'node25, x,y,z:',xnod,x(xnod),y(xnod),z(xnod)
         xnod = x_ndgln( 27 )
         ewrite(3,*)'node27, x,y,z:',xnod,x(xnod),y(xnod),z(xnod)
         stop 3923
      endif

      Loop_Elements: do ele = 1, totele ! Calculate DETWEI,RA,NX,NY,NZ for element ELE


         if(.false.) then
            ewrite(3,*) '+++++++ELE::', ele
            allocate(x_temp(27))
            allocate(y_temp(27))
            allocate(z_temp(27))
            allocate(x_ndgln_temp(27))
            xstar=-0.
            xfini=+1.
            ystar=-0.
            yfini=+1.
            zstar=-0.
            zfini=+1.
            ! redfine hex nodes
            do id=1,3
               do jd=1,3
                  do kd=1,3
                     quad_cv_iloc=(kd-1)*9 +(jd-1)*3+ id
                     xcoord= xstar + (real(id-1)/2.) * (xfini-xstar)
                     ycoord= ystar + (real(jd-1)/2.) * (yfini-ystar)
                     zcoord= zstar + (real(kd-1)/2.) * (zfini-zstar)
                     x_ndgln_temp(quad_cv_iloc)=quad_cv_iloc
                     x_temp(quad_cv_iloc)=xcoord
                     y_temp(quad_cv_iloc)=ycoord
                     z_temp(quad_cv_iloc)=zcoord
                  end do
               end do
            end do


            call detnlxr( ele, x_temp, y_temp, z_temp, x_ndgln_temp, 1, 27, quad_cv_nloc, quad_cv_ngi, &
                 quad_n, quad_nlx, quad_nly, quad_nlz, quad_cvweight, &
                 detwei, ra, volume, d1, d3, dcyl, &       
                 quad_nx, quad_ny, quad_nz )

         else if (.false.) then
 
            allocate(x_temp(x_nonods))
            allocate(y_temp(x_nonods))
            allocate(z_temp(x_nonods))
            x_temp=x
            y_temp=y
            z_temp=z
            allocate(x_ndgln_temp(27))
            x_ndgln_temp(1:27)=x_ndgln(1:27)

            if(.true.) then
               xstar=-0.
               xfini=+1.
               ystar=-0.
               yfini=+1.
               zstar=-0.
               zfini=+1.
               ! redfine hex nodes
               do id=1,3
                  do jd=1,3
                     do kd=1,3
                        quad_cv_iloc=(kd-1)*9 +(jd-1)*3+ id
                        xcoord= xstar + (real(id-1)/2.) * (xfini-xstar)
                        ycoord= ystar + (real(jd-1)/2.) * (yfini-ystar)
                        zcoord= zstar + (real(kd-1)/2.) * (zfini-zstar)
                        x_ndgln_temp(quad_cv_iloc)=quad_cv_iloc
                        x_temp(quad_cv_iloc)=xcoord
                        y_temp(quad_cv_iloc)=ycoord
                        z_temp(quad_cv_iloc)=zcoord
                     end do
                  end do
               end do
               x_temp(1)=0.5
               !       x_temp(5)=0.5
            endif

            allocate(nod_pt(27))
            !  nod1 to nod8 are the corners of the hex:
            nod1=x_ndgln_temp((ele-1)*quad_cv_nloc+1) 
            nod2=x_ndgln_temp((ele-1)*quad_cv_nloc+3) 
            nod3=x_ndgln_temp((ele-1)*quad_cv_nloc+7) 
            nod4=x_ndgln_temp((ele-1)*quad_cv_nloc+9) 

            nod5=x_ndgln_temp((ele-1)*quad_cv_nloc+19) 
            nod6=x_ndgln_temp((ele-1)*quad_cv_nloc+21) 
            nod7=x_ndgln_temp((ele-1)*quad_cv_nloc+25) 
            nod8=x_ndgln_temp((ele-1)*quad_cv_nloc+27) 
            nod_pt(1:27)=x_ndgln_temp((ele-1)*quad_cv_nloc+1:(ele-1)*quad_cv_nloc+27)
            ! redefine hex:
            ! 1st level:
            x_temp(nod_pt(2))=0.5*(x_temp(nod1)+x_temp(nod2))
            x_temp(nod_pt(4))=0.5*(x_temp(nod1)+x_temp(nod3))
            x_temp(nod_pt(6))=0.5*(x_temp(nod2)+x_temp(nod4))
            x_temp(nod_pt(8))=0.5*(x_temp(nod3)+x_temp(nod4))
            x_temp(nod_pt(5))=0.25*( x_temp(nod1)+x_temp(nod2) +x_temp(nod3)+x_temp(nod4) )

            ! 3rd level:
            x_temp(nod_pt(2+18))=0.5*(x_temp(nod5)+x(nod6))
            x_temp(nod_pt(4+18))=0.5*(x_temp(nod5)+x_temp(nod7))
            x_temp(nod_pt(6+18))=0.5*(x_temp(nod6)+x_temp(nod8))
            x_temp(nod_pt(8+18))=0.5*(x_temp(nod7)+x_temp(nod8))
            x_temp(nod_pt(5+18))=0.25*( x_temp(nod5)+x_temp(nod6) +x_temp(nod7)+x_temp(nod8) )

            ! 2nd level:
            do id=1,9 ! average of the 2 levels...
               x_temp(nod_pt(id+9))=0.5*(x_temp(nod_pt(id))+x_temp(nod_pt(id+18)))
            end do

            ! y:
            ! 1st level:
            y_temp(nod_pt(2))=0.5*(y_temp(nod1)+y_temp(nod2))
            y_temp(nod_pt(4))=0.5*(y_temp(nod1)+y_temp(nod3))
            y_temp(nod_pt(6))=0.5*(y_temp(nod2)+y_temp(nod4))
            y_temp(nod_pt(8))=0.5*(y_temp(nod3)+y_temp(nod4))
            y_temp(nod_pt(5))=0.25*( y_temp(nod1)+y_temp(nod2) +y_temp(nod3)+y_temp(nod4) )

            ! 3rd level:
            y_temp(nod_pt(2+18))=0.5*(y_temp(nod5)+y_temp(nod6))
            y_temp(nod_pt(4+18))=0.5*(y_temp(nod5)+y_temp(nod7))
            y_temp(nod_pt(6+18))=0.5*(y_temp(nod6)+y_temp(nod8))
            y_temp(nod_pt(8+18))=0.5*(y_temp(nod7)+y_temp(nod8))
            y_temp(nod_pt(5+18))=0.25*( y_temp(nod5)+y_temp(nod6) +y_temp(nod7)+y_temp(nod8) )

            ! 2nd level:
            do id=1,9 ! average of the 2 levels...
               y_temp(nod_pt(id+9))=0.5*(y_temp(nod_pt(id))+y_temp(nod_pt(id+18)))
            end do

            ! z: 
            if(.false.) then
               z_temp=0.0
               z_temp(nod_pt(2+18))=0.1
               z_temp(nod_pt(4+18))=0.1
               z_temp(nod_pt(6+18))=0.1
               z_temp(nod_pt(8+18))=0.1
               z_temp(nod_pt(5+18))=0.1
            else
               ! 1st level:
               z_temp(nod_pt(2))=0.5*(z_temp(nod1)+z_temp(nod2))
               z_temp(nod_pt(4))=0.5*(z_temp(nod1)+z_temp(nod3))
               z_temp(nod_pt(6))=0.5*(z_temp(nod2)+z_temp(nod4))
               z_temp(nod_pt(8))=0.5*(z_temp(nod3)+z_temp(nod4))
               z_temp(nod_pt(5))=0.25*( z_temp(nod1)+z_temp(nod2) +z_temp(nod3)+z_temp(nod4) )

               ! 3rd level:
               z_temp(nod_pt(2+18))=0.5*(z_temp(nod5)+z_temp(nod6))
               z_temp(nod_pt(4+18))=0.5*(z_temp(nod5)+z_temp(nod7))
               z_temp(nod_pt(6+18))=0.5*(z_temp(nod6)+z_temp(nod8))
               z_temp(nod_pt(8+18))=0.5*(z_temp(nod7)+z_temp(nod8))
               z_temp(nod_pt(5+18))=0.25*( z_temp(nod5)+z_temp(nod6) +z_temp(nod7)+z_temp(nod8) )
            endif

            ! 2nd level:
            do id=1,9 ! average of the 2 levels...
               z_temp(nod_pt(id+9))=0.5*(z_temp(nod_pt(id))+z_temp(nod_pt(id+18)))
            end do


            ewrite(3,*)'d1, d3, dcyl:',d1, d3, dcyl
            !         call detnlxr( ele, x, y, z, x_ndgln, totele, x_nonods, quad_cv_nloc, quad_cv_ngi, &
            call detnlxr( ele, x_temp, y_temp, z_temp, x_ndgln_temp, 1, x_nonods, quad_cv_nloc, quad_cv_ngi, &
                 quad_n, quad_nlx, quad_nly, quad_nlz, quad_cvweight, &
                 detwei, ra, volume, d1, d3, dcyl, &       
                 quad_nx, quad_ny, quad_nz )

            !  nod1 to nod8 are the corners of the hex:
            nod1=x_ndgln_temp((ele-1)*quad_cv_nloc+1) 
            nod2=x_ndgln_temp((ele-1)*quad_cv_nloc+3) 
            nod3=x_ndgln_temp((ele-1)*quad_cv_nloc+7) 
            nod4=x_ndgln_temp((ele-1)*quad_cv_nloc+9) 

            nod5=x_ndgln_temp((ele-1)*quad_cv_nloc+19) 
            nod6=x_ndgln_temp((ele-1)*quad_cv_nloc+21) 
            nod7=x_ndgln_temp((ele-1)*quad_cv_nloc+25) 
            nod8=x_ndgln_temp((ele-1)*quad_cv_nloc+27) 

            ! calculate the volume of the hex: 

            if(.true.) then
               rsum1=Volume_TetHex( .false., &
                    x_temp(nod1), x_temp(nod2), x_temp(nod3), x_temp(nod7), &
                    y_temp(nod1), y_temp(nod2), y_temp(nod3), y_temp(nod7), &
                    z_temp(nod1), z_temp(nod2), z_temp(nod3), z_temp(nod7) )

               rsum2=Volume_TetHex( .false., &
                    x_temp(nod5), x_temp(nod6), x_temp(nod7), x_temp(nod1), &
                    y_temp(nod5), y_temp(nod6), y_temp(nod7), y_temp(nod1), &
                    z_temp(nod5), z_temp(nod6), z_temp(nod7), z_temp(nod1) )

               rsum3=Volume_TetHex( .false., &
                    x_temp(nod6), x_temp(nod7), x_temp(nod2), x_temp(nod1), &
                    y_temp(nod6), y_temp(nod7), y_temp(nod2), y_temp(nod1), &
                    z_temp(nod6), z_temp(nod7), z_temp(nod2), z_temp(nod1) )

               rsum4=Volume_TetHex( .false., &
                    x_temp(nod2), x_temp(nod3), x_temp(nod4), x_temp(nod8), &
                    y_temp(nod2), y_temp(nod3), y_temp(nod4), y_temp(nod8), &
                    z_temp(nod2), z_temp(nod3), z_temp(nod4), z_temp(nod8) )

               rsum5=Volume_TetHex( .false., &
                    x_temp(nod6), x_temp(nod7), x_temp(nod8), x_temp(nod2), &
                    y_temp(nod6), y_temp(nod7), y_temp(nod8), y_temp(nod2), &
                    z_temp(nod6), z_temp(nod7), z_temp(nod8), z_temp(nod2) )

               rsum6=Volume_TetHex( .false., &
                    x_temp(nod7), x_temp(nod3), x_temp(nod2), x_temp(nod8), &
                    y_temp(nod7), y_temp(nod3), y_temp(nod2), y_temp(nod8), &
                    z_temp(nod7), z_temp(nod3), z_temp(nod2), z_temp(nod8) )
            else
               rsum1=Volume_TetHex( .false., &
                    x(nod1), x(nod2), x(nod3), x(nod7), &
                    y(nod1), y(nod2), y(nod3), y(nod7), &
                    z(nod1), z(nod2), z(nod3), z(nod7) )

               rsum2=Volume_TetHex( .false., &
                    x(nod5), x(nod6), x(nod7), x(nod1), &
                    y(nod5), y(nod6), y(nod7), y(nod1), &
                    z(nod5), z(nod6), z(nod7), z(nod1) )

               rsum3=Volume_TetHex( .false., &
                    x(nod6), x(nod7), x(nod2), x(nod1), &
                    y(nod6), y(nod7), y(nod2), y(nod1), &
                    z(nod6), z(nod7), z(nod2), z(nod1) )

               rsum4=Volume_TetHex( .false., &
                    x(nod2), x(nod3), x(nod4), x(nod8), &
                    y(nod2), y(nod3), y(nod4), y(nod8), &
                    z(nod2), z(nod3), z(nod4), z(nod8) )

               rsum5=Volume_TetHex( .false., &
                    x(nod6), x(nod7), x(nod8), x(nod2), &
                    y(nod6), y(nod7), y(nod8), y(nod2), &
                    z(nod6), z(nod7), z(nod8), z(nod2) )

               rsum6=Volume_TetHex( .false., &
                    x(nod7), x(nod3), x(nod2), x(nod8), &
                    y(nod7), y(nod3), y(nod2), y(nod8), &
                    z(nod7), z(nod3), z(nod2), z(nod8) )
            endif

            ewrite(3,*)'rsum1,rsum2,rsum3,rsum4,rsum5,rsum6:', &
                 rsum1,rsum2,rsum3,rsum4,rsum5,rsum6

            ewrite(3,*) 'should be 1/32=', 1./32.
            ewrite(3,*)'rsum1+rsum2+rsum3+rsum4+rsum5+rsum6:', &
                 abs(rsum1) + abs(rsum2) + abs(rsum3) + &
                 abs(rsum4) + abs(rsum5) + abs(rsum6)

            ewrite(3,*)'detwei, volume:',detwei, volume

            !stop 2921


         endif


         call detnlxr( ele, x, y, z, x_ndgln, totele, x_nonods, quad_cv_nloc, quad_cv_ngi, &
              quad_n, quad_nlx, quad_nly, quad_nlz, quad_cvweight, &
              detwei, ra, volume, d1, d3, dcyl, &       
              quad_nx, quad_ny, quad_nz )
         
         ewrite(3,*)'quad_cv_ngi=',quad_cv_ngi
         ewrite(3,*)'detwei for ele=:', ele, detwei
! adjust the volume so that we get the volume correct: 
         if(d3.and.(quad_cv_ngi==1)) detwei=detwei*(2.0310311939315246**3/8.0)
!         stop 3821

         do quad_cv_iloc = 1, quad_cv_nloc
            ewrite(3,*)' quad_cv_iloc: ', quad_cv_iloc

            ewrite(3,*)'quad_nlx:', ( quad_nlx( quad_cv_iloc, quad_cv_gi ), &
                 quad_cv_gi = 1, quad_cv_ngi )

            ewrite(3,*)'quad_nx:', ( quad_nx( quad_cv_iloc, quad_cv_gi ), &
                 quad_cv_gi = 1, quad_cv_ngi )
            ewrite(3,*)'quad_ny:', ( quad_ny( quad_cv_iloc, quad_cv_gi ), &
                 quad_cv_gi = 1, quad_cv_ngi )
            ewrite(3,*)'quad_nz:', ( quad_nz( quad_cv_iloc, quad_cv_gi ), &
                 quad_cv_gi = 1, quad_cv_ngi )
            ewrite(3,*)'detwei:', ( detwei( quad_cv_gi ), quad_cv_gi = 1, quad_cv_ngi )
         end do

         ewrite(3,*) 'vol:', sum(detwei)

!         stop 3838

         Loop_NGI: do quad_cv_gi = 1, quad_cv_ngi ! Determine the quadrature points and weights
            cv_gi = ( ele - 1 ) * quad_cv_ngi + quad_cv_gi
            ! the weights need to sum to 0.5 in 2D triangles and 1./6. in 3D tetrahedra
            if( ndim == 1)  cvweigh( cv_gi ) = detwei( quad_cv_gi )
            if( ndim == 2 ) cvweigh( cv_gi ) = 0.5 * detwei( quad_cv_gi )
            if( ndim == 3 ) cvweigh( cv_gi ) = ( 1. / 6. ) * detwei( quad_cv_gi )
            xgi = 0.
            ygi = 0.
            zgi = 0.
            do quad_cv_iloc = 1, quad_cv_nloc
               xnod = x_ndgln( ( ele - 1 ) * quad_cv_nloc + quad_cv_iloc )
               xgi = xgi + quad_n( quad_cv_iloc, quad_cv_gi ) * x( xnod )
               ygi = ygi + quad_n( quad_cv_iloc, quad_cv_gi ) * y( xnod )
               if( d3 ) zgi = zgi + quad_n( quad_cv_iloc, quad_cv_gi ) * z( xnod )
            end do

            !ewrite(3,*)  'quad_cv_gi, xgi,  ygi, zgi', quad_cv_gi, xgi,  ygi, zgi
            !ewrite(3,*) lx(1:4)
            !ewrite(3,*) ly(1:4)
            !ewrite(3,*) lz(1:4)
            !ewrite(3,*) volume_quad_map( 1, sum(lx(1:4))/4., sum(ly(1:4))/4., sum(lz(1:4))/4., lx, ly, lz )
            !ewrite(3,*) volume_quad_map( 2, sum(lx(1:4))/4., sum(ly(1:4))/4., sum(lz(1:4))/4., lx, ly, lz )
            !ewrite(3,*) volume_quad_map( 3, sum(lx(1:4))/4., sum(ly(1:4))/4., sum(lz(1:4))/4., lx, ly, lz )
            !ewrite(3,*) volume_quad_map( 4, sum(lx(1:4))/4., sum(ly(1:4))/4., sum(lz(1:4))/4., lx, ly, lz )

            if( d3 ) then ! local coords for the triangle
               quad_l1( cv_gi ) = volume_quad_map( 1, xgi, ygi, zgi, lx(1:4), ly(1:4), lz(1:4) )
               quad_l2( cv_gi ) = volume_quad_map( 2, xgi, ygi, zgi, lx(1:4), ly(1:4), lz(1:4) )
               quad_l3( cv_gi ) = volume_quad_map( 3, xgi, ygi, zgi, lx(1:4), ly(1:4), lz(1:4) )
               quad_l4( cv_gi ) = volume_quad_map( 4, xgi, ygi, zgi, lx(1:4), ly(1:4), lz(1:4) )
            else
               quad_l1( cv_gi ) = area_quad_map( 1, xgi, ygi, lx( 1 : 3 ), ly( 1 : 3 ) )
               quad_l2( cv_gi ) = area_quad_map( 2, xgi, ygi, lx( 1 : 3 ), ly( 1 : 3 ) )
               quad_l3( cv_gi ) = area_quad_map( 3, xgi, ygi, lx( 1 : 3 ), ly( 1 : 3 ) )
            end if

         end do Loop_NGI
      end do Loop_Elements

      ! Now determine the basis functions and derivatives at the 
      ! quadrature pts quad_L1, quad_L2, quad_L3, quad_L4, etc
      call shatri_hex( quad_l1, quad_l2, quad_l3, quad_l4, rdummy, d3, &
           cv_nloc, cv_ngi, n, nlx, nly, nlz, &
           .true. )
      ewrite(3,*)'cvweigh:',cvweigh
      ewrite(3,*)'nlx(1,:):',nlx(1,:)
      rsum=-1.e+5
      do cv_gi = 1, cv_ngi
         rsum = max(rsum,nlx(1,cv_gi))
      end do
      ewrite(3,*)'max(nlx(1,:)):',rsum
      rsum=+1.e+5
      do cv_gi = 1, cv_ngi
         rsum = min(rsum,nlx(1,cv_gi))
      end do
      ewrite(3,*)'min(nlx(1,:)):',rsum

      ewrite(3,*)'sum(nlx(1,:)):',sum(nlx(1,:))
      rsum = 0.0
      do cv_gi = 1, cv_ngi
         rsum = rsum + cvweigh( cv_gi )
      end do
      ewrite(3,*)'rsum:', rsum
      ewrite(3,*)'sum(cvweigh):',sum(cvweigh)
      !stop 2921

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



    subroutine shatri_hex( l1, l2, l3, l4, weight, d3, &
         nloc, ngi, &
         n, nlx, nly, nlz, &
         tri_tet )
      implicit none
      integer, intent( in ) :: nloc, ngi
      logical, intent( in ) :: tri_tet
      real, dimension( ngi ), intent( in ) :: l1, l2, l3, l4
      real, dimension( ngi ), intent( inout ) :: weight
      logical, intent( in ) :: d3
      real, dimension( nloc, ngi ), intent( inout ) :: n, nlx, nly, nlz
      ! Local variables
      logical :: lowqua
      integer :: nwicel, mloc, snloc, sngi
      real, dimension( : ), allocatable :: rdum, sweigh
      real, dimension( :, : ), allocatable ::  m, sn, snlx, snly

      ewrite(3,*)'In shatri_hex'

      Conditional_Dimension: if( tri_tet ) then ! traingles and tets

         call shatri( l1, l2, l3, l4, weight, d3, &
              nloc, ngi, &
              n, nlx, nly, nlz )

      else
         ! Get the shape functions on lines (in 2D) and quadrilateral surfaces in 3D 
         if( .not. d3 ) then ! 2D:
            if( nloc == 4 ) nwicel = 1
            if( nloc == 9 ) nwicel = 3
         else
            if( nloc == 8 ) nwicel = 1
            if( nloc == 27 ) nwicel = 3
         endif
         lowqua = .false. 
         mloc = 0
         snloc = 0
         sngi = 0

         Select Case( nwicel )

         case default; FLExit( "Wrong option for NWICEL " )

         case( 1 )
            allocate( rdum( 0 ) )
            allocate( m( mloc, ngi ) )
            allocate( sweigh( sngi ) )
            allocate( sn( snloc, sngi ) )
            allocate( snlx( snloc, sngi ) )
            allocate( snly( snloc, sngi ) )
            if( .not. d3 ) then
               call re2dn4( lowqua, ngi, 0, nloc, mloc, &
                    m, weight, n, nlx, nly, &
                    sngi, snloc, sweigh, sn, snlx, &
                    rdum, rdum )
            else
               call re3dn8( lowqua, ngi, 0, nloc, mloc, &
                    m, weight, n, nlx, nly, nlz, &
                    sngi, snloc, sweigh, sn, snlx, snly, &
                    rdum, rdum, rdum )
            end if
            deallocate( rdum )
            deallocate( m )
            deallocate( sweigh )
            deallocate( sn )
            deallocate( snlx )
            deallocate( snly )

         case( 3 )
            if( .not. d3 ) then
               call re2dn9( lowqua, ngi, 0, nloc, mloc, &
                    m, weight, n, nlx, nly, &
                    rdum, rdum )
            else ! Lagrange 27 nodes 3D element - bilinear pressure
               call re3d27( lowqua, ngi, 0, nloc, mloc, &
                    m, weight, n, nlx, nly, nlz, &
                    rdum, rdum, rdum )
            end if

         end Select

      endif Conditional_Dimension

      return
    end subroutine shatri_hex


    subroutine shatri( l1, l2, l3, l4, weight, d3, &
         nloc, ngi, &
         n, nlx, nly, nlz )
      implicit none
      integer, intent( in ) :: nloc, ngi
      real, dimension( ngi ), intent( in ) :: l1, l2, l3, l4, weight
      logical, intent( in ) :: d3
      real, dimension( nloc, ngi ), intent( inout ) :: n, nlx, nly, nlz
      ! Local variables
      logical :: base_order
      integer :: gi, ndim, cv_ele_type_dummy, u_nloc_dummy
      real, dimension( :, : ), allocatable :: cvn_dummy, un_dummy, unlx_dummy, &
           unly_dummy, unlz_dummy
      real, dimension( : ), allocatable :: cvweigh_dummy

      ewrite(3,*)'In shatri d3,nloc=',d3,nloc

      Conditional_Dimensionality: if( .not. d3 ) then ! Assume a triangle

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

            base_order=.true.
            if(base_order) then
               ! order so that the 1st nodes are on the base...
              call base_order_tri(n,nloc,ngi)
              call base_order_tri(nlx,nloc,ngi)
              call base_order_tri(nly,nloc,ngi)
            endif

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
            ewrite(3,*)'going into SHATRIold'
            call SHATRIold(L1, L2, L3, L4, cvweigh_dummy, (ndim==3), &
     &               NLOC,NGI,  &
     &               N,NLX,NLY,NLZ)

            base_order=.true.
            if(base_order) then
               ! order so that the 1st nodes are on the base...
              call base_order_tet(n,nloc,ngi)
              call base_order_tet(nlx,nloc,ngi)
              call base_order_tet(nly,nloc,ngi)
              call base_order_tet(nlz,nloc,ngi)
            endif

            deallocate( cvn_dummy )
            deallocate( un_dummy )
            deallocate( unlx_dummy )
            deallocate( unly_dummy )
            deallocate( unlz_dummy )
            deallocate( cvweigh_dummy )

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

      return
    end subroutine shatri



    subroutine base_order_tri(n,nloc,ngi)
      ! order so that the 1st nodes are on the base for a 
      ! quadratic triangle...
      implicit none
      integer, intent( in ) :: nloc, ngi
      real, dimension( nloc, ngi ), intent( inout ) :: n
      ! local variables...
      real, dimension( :, : ), allocatable :: rn
      integer, dimension( : ), allocatable :: old2new
      integer :: iloc

      allocate(rn(nloc,ngi))
      allocate(old2new(nloc))
      rn=n

      old2new(1)=1
      old2new(2)=4
      old2new(3)=2
      old2new(4)=6
      old2new(5)=5
      old2new(6)=3

    do iloc=1,nloc
      n(iloc,:)=rn(old2new(iloc),:)
    end do
      return
    end subroutine base_order_tri



    subroutine base_order_tet(n,nloc,ngi)
      ! order so that the 1st nodes are on the base of tet for 
      ! a quadratic tet...
      implicit none
      integer, intent( in ) :: nloc, ngi
      real, dimension( nloc, ngi ), intent( inout ) :: n
      ! local variables...
      real, dimension( :, : ), allocatable :: rn
      integer, dimension( : ), allocatable :: old2new
      integer :: iloc

      allocate(rn(nloc,ngi))
      allocate(old2new(nloc))
      rn=n

      old2new(1)=1
      old2new(2)=2
      old2new(3)=3
      old2new(4)=6
      old2new(5)=4
      old2new(6)=5

      old2new(7)=7
      old2new(8)=8
      old2new(9)=9
      old2new(10)=10

    do iloc=1,nloc
      n(iloc,:)=rn(old2new(iloc),:)
    end do
      return
    end subroutine base_order_tet




    subroutine Calc_CVN_TriTetQuadHex( cv_ele_type, totele, cv_nloc, cv_ngi, x_nonods, &
         quad_cv_nloc, x_ndgln, fem_nod, cvn )
      ! Compute CVN (CV basis function) for triangles, tetrahedra, quadrilaterals and hexahedra
      implicit none
      integer, intent( in ) :: cv_ele_type, totele, cv_nloc, cv_ngi, &
           x_nonods, quad_cv_nloc
      integer, dimension( totele * quad_cv_nloc ), intent( in ) :: x_ndgln
      integer, dimension( x_nonods ), intent( in ) :: fem_nod
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: cvn
      ! Local variables
      integer :: ele, nod, quad_cv_iloc, xnod, quad_cv_gi, cv_gi, &
           quad_cv_ngi

      ewrite(3,*)'In Calc_CVN_TriTetQuadHex'

      cvn = 0.
      Loop_Elements: do ele = 1, totele
         nod = 0
         do quad_cv_iloc = 1, quad_cv_nloc
            xnod = x_ndgln( ( ele - 1 ) * quad_cv_nloc + quad_cv_iloc )
            if( fem_nod( xnod ) /= 0 ) nod = fem_nod( xnod )
            ewrite(3,*) 'ele, xnod, fem_nod, nod:', ele, xnod, fem_nod( xnod ), nod
         end do

         if( nod == 0 ) FLExit(" Problem with CVN calculation " )

            quad_cv_ngi = cv_ngi / totele
            do quad_cv_gi = 1, quad_cv_ngi
               cv_gi = ( ele - 1 ) * quad_cv_ngi + quad_cv_gi
               cvn( nod, cv_gi ) = 1.
            end do

      end do Loop_Elements

      return
    end subroutine Calc_CVN_TriTetQuadHex


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
      integer :: ipoly, iqadra, iloc
      real, dimension( : ), allocatable :: rdum 

      ewrite(3,*)' In shape_l_q_quad', nwicel, d3

      allocate( rdum( 1 ) )

      Select Case( nwicel )

      case default; FLExit( "Wrong option for NWICEL " )

      case( 1 )
         if( .not. d3 ) then
            call re2dn4( lowqua, ngi, 0, nloc, mloc, &
                 m, weight, n, nlx, nly, &
                 sngi, snloc, sweigh, sn, snlx, &
                 rdum, rdum )
         else
            call re3dn8( lowqua, ngi, 0, nloc, mloc, &
                 m, weight, n, nlx, nly, nlz, &
                 sngi, snloc, sweigh, sn, snlx, snly, &
                 rdum, rdum, rdum ) 
         end if

      case( 3 )
         if( .not. d3 ) then
            call re2dn9( lowqua, ngi, 0, nloc, mloc, &
                 m, weight, n, nlx, nly, &
                 rdum, rdum )

            call quad_basis_funs_1d(sngi, snloc,  &
                 sweigh, sn, snlx )

         else ! Lagrange 27 nodes 3D element - bilinear pressure
            call re3d27( lowqua, ngi, 0, nloc, mloc, &
                 m, weight, n, nlx, nly, nlz, &
                 rdum, rdum, rdum )

            call re2dn9( lowqua, sngi, 0, snloc, mloc, &
                 m, sweigh, sn, snlx, snly, &
                 rdum, rdum )
         end if

      end Select

      deallocate( rdum )

      return
    end subroutine shape_l_q_quad

!!!-                                      -!!!
!!!-     And the functions:    -!!!
!!!-                                      -!!!


    real function volume_quad_map( cv_iloc, xgi, ygi, zgi, lx, ly, lz )
      ! Compute the cv_iloc^{th} shape function value at point (xgi, ygi, zgi)
      implicit none
      integer :: cv_iloc
      real :: xgi, ygi, zgi
      real, dimension( : ) :: lx, ly, lz
      ! Local variables
      real :: lxyz( size(lx), 3 ), xyzgi( 3 )
      real :: loc_x_coord, loc_y_coord, loc_z_coord, loc_zz_coord

      lxyz( :, 1 ) = lx( : ) 
      lxyz( :, 2 ) = ly( : )
      lxyz( :, 3 ) = lz( : )

      xyzgi( : ) = (/ xgi, ygi, zgi /)

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

    subroutine PrintOutFunMat( n, m, a )
      implicit none
      integer, intent( in ) :: n, m
      real, dimension( n, m ), intent( in ) :: a
      ! Local variables
      integer :: in, im

      do in = 1, n
         ewrite(3,*) in, ( a( in, im ), im = 1, m )
      end do

      return
    end subroutine PrintOutFunMat


    SUBROUTINE DGSDETNXLOC2( SNLOC, SNGI, &
         XSL, YSL, ZSL, &
         SN, SNLX, SNLY, SWEIGH, SDETWE, SAREA, &
         D1, D3, DCYL, &
         NORMXN, NORMYN, NORMZN, &
         NORMX, NORMY, NORMZ )
      IMPLICIT NONE

      INTEGER, intent( in ) :: SNLOC, SNGI
      REAL, DIMENSION( SNLOC ), intent( in ) :: XSL, YSL, ZSL
      REAL, DIMENSION( SNLOC, SNGI ), intent( in ) :: SN, SNLX, SNLY
      REAL, DIMENSION( SNGI ), intent( in ) :: SWEIGH
      REAL, DIMENSION( SNGI ), intent( inout ) :: SDETWE 
      REAL, intent( inout ) ::  SAREA
      LOGICAL, intent( in ) ::  D1,D3,DCYL
      REAL, DIMENSION( SNGI ), intent( inout ) :: NORMXN, NORMYN, NORMZN
      REAL, intent( in ) :: NORMX, NORMY, NORMZ
      ! Local variables
      real, parameter :: pi = 3.141592654
      INTEGER :: GI, SL, IGLX
      REAL :: DXDLX, DXDLY, DYDLX, DYDLY, DZDLX, DZDLY
      REAL :: A, B, C, DETJ, RGI, TWOPI

      SAREA=0.

      IF(D3) THEN
         DO GI=1,SNGI

            DXDLX=0.
            DXDLY=0.
            DYDLX=0.
            DYDLY=0.
            DZDLX=0.
            DZDLY=0.

            DO SL=1,SNLOC
               DXDLX=DXDLX + SNLX(SL,GI)*XSL(SL)
               DXDLY=DXDLY + SNLY(SL,GI)*XSL(SL)
               DYDLX=DYDLX + SNLX(SL,GI)*YSL(SL)
               DYDLY=DYDLY + SNLY(SL,GI)*YSL(SL)
               DZDLX=DZDLX + SNLX(SL,GI)*ZSL(SL)
               DZDLY=DZDLY + SNLY(SL,GI)*ZSL(SL)
            END DO
            A = DYDLX*DZDLY - DYDLY*DZDLX
            B = DXDLX*DZDLY - DXDLY*DZDLX
            C = DXDLX*DYDLY - DXDLY*DYDLX

            DETJ=SQRT( A**2 + B**2 + C**2)
            SDETWE(GI)=DETJ*SWEIGH(GI)
            SAREA=SAREA+SDETWE(GI)

            ! Calculate the normal at the Gauss pts...
            ! Perform x-product. N=T1 x T2
            CALL NORMGI(NORMXN(GI),NORMYN(GI),NORMZN(GI), &
                 DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY, &
                 NORMX,NORMY,NORMZ)
         END DO
      ELSE IF(.NOT.D1) THEN
         TWOPI=1.0
         IF(DCYL) TWOPI=2.*PI

         DO GI=1,SNGI
            RGI=0.
            DXDLX=0.
            DXDLY=0.
            DYDLX=0.
            DYDLY=0.
            DZDLX=0.
            ! DZDLY=1 is to calculate the normal.
            DZDLY=1.
            DO SL=1,SNLOC
               DXDLX=DXDLX + SNLX(SL,GI)*XSL(SL)
               DYDLX=DYDLX + SNLX(SL,GI)*YSL(SL)
               RGI=RGI+SN(SL,GI)*YSL(SL)
               !RGI=RGI+SN(SL,GI)*YSL(IGLX)
            END DO
            IF(.NOT.DCYL) RGI=1.0
            DETJ=SQRT( DXDLX**2 + DYDLX**2 )
            SDETWE(GI)=TWOPI*RGI*DETJ*SWEIGH(GI)
            SAREA=SAREA+SDETWE(GI)
            CALL NORMGI(NORMXN(GI),NORMYN(GI),NORMZN(GI), &
                 DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY, &
                 NORMX,NORMY,NORMZ)
         END DO

      ELSE ! For 1D...
         DO GI = 1, SNGI
            DXDLX = 0.
            DO SL = 1, SNLOC
               DXDLX = DXDLX + SNLX( SL, GI ) * XSL( SL )
            END DO
            SDETWE( GI ) = SWEIGH( GI )
            SAREA = SAREA + SDETWE( GI )
            NORMXN( GI ) = NORMX
            NORMYN( GI ) = 0.0
            NORMZN( GI ) = 0.0
         END DO

      ENDIF

      RETURN

    END SUBROUTINE DGSDETNXLOC2

    SUBROUTINE NORMGI( NORMXN, NORMYN, NORMZN, &
         DXDLX, DYDLX, DZDLX, DXDLY, DYDLY, DZDLY, &
         NORMX, NORMY, NORMZ) 
      ! Calculate the normal at the Gauss pts
      ! Perform x-product. N=T1 x T2
      implicit none
      REAL, intent( inout ) :: NORMXN, NORMYN, NORMZN
      REAL, intent( in )    :: DXDLX, DYDLX, DZDLX, DXDLY, DYDLY, DZDLY
      REAL, intent( in )    :: NORMX, NORMY, NORMZ
      ! Local variables
      REAL :: RN, SIRN

      CALL XPROD( NORMXN, NORMYN, NORMZN, &
           DXDLX, DYDLX, DZDLX, &
           DXDLY, DYDLY, DZDLY )

      RN = SQRT( NORMXN**2 + NORMYN**2 + NORMZN**2 )

      SIRN = SIGN( 1.0 / RN, NORMXN * NORMX + NORMYN * NORMY + NORMZN * NORMZ )

      NORMXN = SIRN * NORMXN
      NORMYN = SIRN * NORMYN
      NORMZN = SIRN * NORMZN

      RETURN

    END SUBROUTINE NORMGI


    SUBROUTINE XPROD( AX, AY, AZ, &
         BX, BY, BZ, &
         CX, CY, CZ )
      implicit none
      REAL, intent( inout ) :: AX, AY, AZ
      REAL, intent( in )    :: BX, BY, BZ, CX, CY, CZ

      ! Perform x-product. a=b x c 
      AX =    BY * CZ - BZ * CY
      AY = -( BX * CZ - BZ * CX )
      AZ =    BX * CY - BY * CX

      RETURN
    END subroutine XPROD


!!!!

    subroutine Make_QTri( totele, x_nloc, max_x_nonods, x_nonods, &
         x_ndgln, lx, ly, x, y, fem_nod )
      implicit none
      integer, intent( in ) :: totele, x_nloc, max_x_nonods
      integer, intent( inout ) :: x_nonods
      integer, dimension( totele * x_nloc ), intent( inout ) :: x_ndgln
      real, dimension( max_x_nonods ), intent( inout ) :: lx, ly, x, y
      integer, dimension( max_x_nonods ), intent( inout ) :: fem_nod
      ! Local variables
      ! Scaling factor to give a unity area of the local triangle
      real, parameter :: h_scale =  1.5196713713031851
      integer, parameter :: x_nloc_big = 7, totele_big = 4, &
           nonods_big = x_nloc_big * totele_big
      logical, dimension( : ), allocatable :: node_belong
      integer, dimension( : ), allocatable :: x_ndgln_big, old2new
      real, dimension( : ), allocatable :: x_big, y_big
      integer :: ele_big, increment_ele_big, ele, x_iloc, elebig_ele, iloc, xnod, &
           ele_ref, mx_x_nonods, inod, count_nod

      ewrite(3,*)' In Make_QTri'

      ! Setting-up unity area triangle
      lx( 1 ) = 0.
      ly( 1 ) = 0.

      lx( 2 ) = 1. * h_scale
      ly( 2 ) = 0.

      lx( 3 ) = .5 * h_scale
      ly( 3 ) = sqrt( 3./4. ) * h_scale

      ! Remaping Big Triangle
      ! Node point 1
      x( 1 ) = lx( 1 )
      y( 1 ) = ly( 1 )
      fem_nod( 1 ) = 1
      ! Node point 3
      x( 3 ) = lx( 2 )
      y( 3 ) = ly( 2 )
      fem_nod( 3 ) = 3
      ! Node point 6
      x( 6 ) = lx( 3 )
      y( 6 ) = ly( 3 )
      fem_nod( 6 ) = 6

      ! Remaping Mid Triangle
      ! Node point 2
      x( 2 ) = 0.5 * ( lx( 1 ) + lx( 2 ) )
      y( 2 ) = 0.5 * ( ly( 1 ) + ly( 2 ) )
      fem_nod( 2 ) = 2
      ! Node point 4
      x( 4 ) = 0.5 * ( lx( 1 ) + lx( 3 ) )
      y( 4 ) = 0.5 * ( ly( 1 ) + ly( 3 ) )
      fem_nod( 4 ) = 4
      ! Node point 5
      x( 5 ) = 0.5 * ( lx( 2 ) + lx( 3 ) )
      y( 5 ) = 0.5 * ( ly( 2 ) + ly( 3 ) )
      fem_nod( 5 ) = 5

      mx_x_nonods=x_nonods

      allocate( x_ndgln_big( x_nloc_big * totele_big ) )
      x_ndgln_big = 0

      allocate( x_big( x_nonods * x_nloc ) ) ;  x_big = 0.
      allocate( y_big( x_nonods * x_nloc ) ) ;  y_big = 0.
      x_big( 1 : x_nonods ) = x ( 1 : x_nonods )
      y_big( 1 : x_nonods ) = y ( 1 : x_nonods )

      ele_big = 1
      x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 1 ) = 1
      x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 2 ) = 2
      x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 3 ) = 4

      ele_big = 2
      x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 1 ) = 2
      x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 2 ) = 4
      x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 3 ) = 5

      ele_big = 3
      x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 1 ) = 2
      x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 2 ) = 3
      x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 3 ) = 5

      ele_big = 4
      x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 1 ) = 4
      x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 2 ) = 5
      x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 3 ) = 6

      increment_ele_big = 6
      ele_ref = 1
      do ele_big = 1, 4
         call Computing_Small_QTriangles( ele_big, x_nloc_big, totele_big, &
              x_nonods * x_nloc, increment_ele_big, x_ndgln_big, &
              x_big, y_big ) 

         call ReMaping_Fields_QTriangles( ele_big, x_nloc_big, totele_big, &
              x_nonods * x_nloc, x_ndgln_big, x_big, y_big, &
              x_nonods, x_nloc, totele, x_ndgln, ele_ref, x, y ) 
      end do

      ! Just eliminating repetitive nodes in the 4 nodes pts of the quads
      call Eliminating_Repetitive_Nodes( totele, x_nloc, x_nonods, .false., &
           x_ndgln, x, y ) 

! the elements are anti-clockwise correct to make the local numbering correct...
      do ele = 1, totele
        xnod=x_ndgln( ( ele - 1 ) * x_nloc + 3 )
        x_ndgln( ( ele - 1 ) * x_nloc + 3 )=x_ndgln( ( ele - 1 ) * x_nloc + 4 )
        x_ndgln( ( ele - 1 ) * x_nloc + 4 )=xnod
      end do
! these elements are inside out so correct them...
      ele=3
      xnod=x_ndgln( ( ele - 1 ) * x_nloc + 1 )
      x_ndgln( ( ele - 1 ) * x_nloc + 1 )=x_ndgln( ( ele - 1 ) * x_nloc + 2 )
      x_ndgln( ( ele - 1 ) * x_nloc + 2 )=xnod
      xnod=x_ndgln( ( ele - 1 ) * x_nloc + 3 )
      x_ndgln( ( ele - 1 ) * x_nloc + 3 )=x_ndgln( ( ele - 1 ) * x_nloc + 4 )
      x_ndgln( ( ele - 1 ) * x_nloc + 4 )=xnod
      ele=4
      xnod=x_ndgln( ( ele - 1 ) * x_nloc + 1 )
      x_ndgln( ( ele - 1 ) * x_nloc + 1 )=x_ndgln( ( ele - 1 ) * x_nloc + 2 )
      x_ndgln( ( ele - 1 ) * x_nloc + 2 )=xnod
      xnod=x_ndgln( ( ele - 1 ) * x_nloc + 3 )
      x_ndgln( ( ele - 1 ) * x_nloc + 3 )=x_ndgln( ( ele - 1 ) * x_nloc + 4 )
      x_ndgln( ( ele - 1 ) * x_nloc + 4 )=xnod
      ele=5
      xnod=x_ndgln( ( ele - 1 ) * x_nloc + 1 )
      x_ndgln( ( ele - 1 ) * x_nloc + 1 )=x_ndgln( ( ele - 1 ) * x_nloc + 2 )
      x_ndgln( ( ele - 1 ) * x_nloc + 2 )=xnod
      xnod=x_ndgln( ( ele - 1 ) * x_nloc + 3 )
      x_ndgln( ( ele - 1 ) * x_nloc + 3 )=x_ndgln( ( ele - 1 ) * x_nloc + 4 )
      x_ndgln( ( ele - 1 ) * x_nloc + 4 )=xnod
      ele=9
      xnod=x_ndgln( ( ele - 1 ) * x_nloc + 1 )
      x_ndgln( ( ele - 1 ) * x_nloc + 1 )=x_ndgln( ( ele - 1 ) * x_nloc + 2 )
      x_ndgln( ( ele - 1 ) * x_nloc + 2 )=xnod
      xnod=x_ndgln( ( ele - 1 ) * x_nloc + 3 )
      x_ndgln( ( ele - 1 ) * x_nloc + 3 )=x_ndgln( ( ele - 1 ) * x_nloc + 4 )
      x_ndgln( ( ele - 1 ) * x_nloc + 4 )=xnod
      ele=12
      xnod=x_ndgln( ( ele - 1 ) * x_nloc + 1 )
      x_ndgln( ( ele - 1 ) * x_nloc + 1 )=x_ndgln( ( ele - 1 ) * x_nloc + 2 )
      x_ndgln( ( ele - 1 ) * x_nloc + 2 )=xnod
      xnod=x_ndgln( ( ele - 1 ) * x_nloc + 3 )
      x_ndgln( ( ele - 1 ) * x_nloc + 3 )=x_ndgln( ( ele - 1 ) * x_nloc + 4 )
      x_ndgln( ( ele - 1 ) * x_nloc + 4 )=xnod

! At this pt we need to shift the node numbers down so there are only 19 nodes: 
      allocate(node_belong(x_nonods)) ; node_belong=.false.
      allocate(old2new(x_nonods)) ; old2new=0
      do ele=1,totele
        do x_iloc=1,4
          inod=x_ndgln( ( ele - 1 ) * x_nloc + x_iloc )
          node_belong(inod)=.true.
        end do
      end do
!
      count_nod=0
      do inod=1,x_nonods
        if(node_belong(inod)) then
           count_nod=count_nod+1
           old2new(inod)=count_nod           
           x(count_nod)=x(inod)
           y(count_nod)=y(inod)
        endif
      end do

      do ele = 1, totele 
         do iloc = 1, 4
            inod = x_ndgln( ( ele - 1 ) * x_nloc + iloc )
            x_ndgln( ( ele - 1 ) * x_nloc + iloc )=old2new(inod)
         end do
      end do
      x_nonods=count_nod


      ewrite(3,*)'xndgln0:'
      do ele = 1, totele
         ewrite(3,*) ele, ( x_ndgln( ( ele - 1 ) * x_nloc + x_iloc ), x_iloc = 1, 4 )
         ewrite(3,*) ele, ( x(x_ndgln( ( ele - 1 ) * x_nloc + x_iloc )), x_iloc = 1, 4 )
         ewrite(3,*) ele, ( y(x_ndgln( ( ele - 1 ) * x_nloc + x_iloc )), x_iloc = 1, 4 )
      end do

      ! Extra Nodes in each Quad (now 9 / quad )
      call Adding_Extra_Parametric_Nodes( totele, x_nloc, mx_x_nonods, &
           x_ndgln, x, y )

      mx_x_nonods=maxval(x_ndgln)
      call Eliminating_Repetitive_Nodes_all( totele, x_nloc, x_nonods, mx_x_nonods,  &
           x_ndgln, x, y )

      x_nonods = maxval( x_ndgln )

      deallocate( x_ndgln_big )
      deallocate( x_big )
      deallocate( y_big )

      return
    end subroutine Make_QTri

    subroutine Computing_Small_QTriangles( ele_big, x_nloc_big, totele_big, &
         x_nonods_big, increment_ele_big, x_ndgln_big,  &
         x_big, y_big ) 
      implicit none
      integer, intent( in ) :: ele_big, x_nloc_big, totele_big, x_nonods_big
      integer, intent( inout ) :: increment_ele_big
      integer, dimension( x_nloc_big * totele_big ), intent( inout ) :: x_ndgln_big
      real, dimension( x_nonods_big ), intent( inout ) :: x_big, y_big
      ! Local variables
      integer :: iloc, xnod1, xnod2, xnod3, iloc2, x_nloc_big2

      ewrite(3,*)' In Computing_Small_QTriangles '

      x_nloc_big2 = x_nloc_big / 2
      iloc2 = x_nloc_big2
      do iloc = 1, x_nloc_big2
         if( iloc < x_nloc_big2 ) then
            xnod1 = x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + iloc )
            xnod2 = x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + iloc + 1 )
         else
            xnod1 = x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + iloc )
            xnod2 = x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 1 )
         end if
         iloc2 = iloc2 + 1
         increment_ele_big = increment_ele_big + 1
         x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + iloc2 ) = &
              increment_ele_big
         x_big( increment_ele_big ) = 0.5 * ( x_big( xnod1 ) + x_big( xnod2 ) )
         y_big( increment_ele_big ) = 0.5 * ( y_big( xnod1 ) + y_big( xnod2 ) )
      end do

      increment_ele_big = increment_ele_big + 1
      iloc2 = iloc2 + 1
      x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + iloc2 ) = &
           increment_ele_big
      xnod1 = x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 1 )
      xnod2 = x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 2 )
      xnod3 = x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 3 )  
      x_big( increment_ele_big ) = 1. / 3. * ( x_big( xnod1 ) + &
           x_big( xnod2 ) + x_big( xnod3 ) )
      y_big( increment_ele_big ) = 1. / 3. * ( y_big( xnod1 ) + &
           y_big( xnod2 ) + y_big( xnod3 ) )

      return
    end subroutine Computing_Small_QTriangles

    subroutine ReMaping_Fields_QTriangles( ele_big, x_nloc_big, totele_big, &
         x_nonods_big, x_ndgln_big, x_big, y_big, &
         x_nonods, x_nloc, totele, x_ndgln, ele_ref, x, y )
      implicit none
      integer, intent( in ) :: ele_big, x_nloc_big, totele_big, x_nonods_big
      integer, dimension( x_nloc_big * totele_big ), intent( in ) :: x_ndgln_big
      real, dimension( x_nonods_big ), intent( in ) :: x_big, y_big
      integer, intent( in ) :: x_nonods, x_nloc, totele
      integer, dimension( x_nloc * totele ), intent( inout ) :: x_ndgln
      integer, intent( inout ) :: ele_ref
      real, dimension( x_nonods ), intent( inout ) :: x, y
      ! Local variables
      integer :: xnod, iloc, x_nloc_big2, ele_ref2, ele

      ewrite(3,*)' In ReMaping_Fields_QTriangles'

      x_nloc_big2 = x_nloc_big / 2

      ! Element 1
      ele = ele_ref
      ele_ref2 = ele_ref
      iloc = 1
      x_ndgln( ( ele - 1 ) * x_nloc + iloc ) = &
           x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 1 )

      iloc = 2
      x_ndgln( ( ele - 1 ) * x_nloc + iloc ) = &
           x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + x_nloc_big2 + 1 )

      iloc = 3
      x_ndgln( ( ele - 1 ) * x_nloc + iloc ) = &
           x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + x_nloc_big2 + 4 )

      iloc = 4
      x_ndgln( ( ele - 1 ) * x_nloc + iloc ) = &
           x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + x_nloc_big2 + 3 )

      ! Element 2
      ele = ele + 1
      iloc = 1
      x_ndgln( ( ele - 1 ) * x_nloc + iloc ) = &
           x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + x_nloc_big2 + 1 )

      iloc = 2
      x_ndgln( ( ele - 1 ) * x_nloc + iloc ) = &
           x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 2 )

      iloc = 3
      x_ndgln( ( ele - 1 ) * x_nloc + iloc ) = &
           x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + x_nloc_big2 + 2 )

      iloc = 4
      x_ndgln( ( ele - 1 ) * x_nloc + iloc ) = &
           x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + x_nloc_big2 + 4 )

      ! Element 3
      ele = ele + 1
      iloc = 1
      x_ndgln( ( ele - 1 ) * x_nloc + iloc ) = &
           x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + x_nloc_big2 + 3 )

      iloc = 2
      x_ndgln( ( ele - 1 ) * x_nloc + iloc ) = &
           x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + 3 )

      iloc = 3
      x_ndgln( ( ele - 1 ) * x_nloc + iloc ) = &
           x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + x_nloc_big2 + 2 )

      iloc = 4
      x_ndgln( ( ele - 1 ) * x_nloc + iloc ) = &
           x_ndgln_big( ( ele_big - 1 ) * x_nloc_big + x_nloc_big2 + 4 )

      ele_ref = ele_ref + 3

      do ele = ele_ref2, ele_ref2 + 2
         do iloc = 1, 4 ! Just the initial quad node points
            xnod = x_ndgln( ( ele - 1 ) * x_nloc + iloc )
            x( xnod ) = x_big( xnod )
            y( xnod ) = y_big( xnod )
         end do
      end do

      return
    end subroutine ReMaping_Fields_QTriangles


    subroutine Eliminating_Repetitive_Nodes( totele, x_nloc, x_nonods, over_all, &
         x_ndgln, x, y )
      implicit none
      integer, intent( in ) :: totele, x_nloc
      integer, intent( in ) :: x_nonods
      logical, intent( in ) :: over_all
      integer, dimension( totele * x_nloc ), intent( inout ) :: x_ndgln
      real, dimension( x_nonods ), intent( in ) :: x, y
      ! Local variables
      real, parameter :: toler = 1.e-5
      integer :: ele, ele2, iloc, jloc, inod, jnod, jnod2, x_nloc2, isum, iref
      integer :: count_nod
      logical :: found
      integer, dimension( : ), allocatable :: x_ndgln2,new2old,old2new
      real, dimension( : ), allocatable :: x2,y2

      ewrite(3,*) 'In Eliminating_Repetitive_Nodes'


    if(.true.) then

      x_nloc2 = x_nloc
      if( .not. over_all)  x_nloc2 = 4

      do ele = 1, totele - 1
         do iloc = 1, x_nloc2
            inod = x_ndgln( ( ele - 1 ) * x_nloc + iloc )
            do ele2 = ele + 1, totele
               do jloc = 1, x_nloc2
                  jnod = x_ndgln( ( ele2 - 1 ) * x_nloc + jloc )
                  if( ( abs( x( inod ) - x( jnod ) ) <= toler ) .and. &
                       ( abs( y( inod ) - y( jnod ) ) <= toler ) ) then
                     x_ndgln( ( ele2 - 1 ) * x_nloc + jloc ) = &
                          x_ndgln( ( ele - 1 ) * x_nloc + iloc )
                  end if
               end do
            end do
         end do
      end do

      if( over_all )then
         isum = 1 
         do iloc = 2, totele * x_nloc
            found = .false. ; jnod2 = 6666
            inod = x_ndgln( iloc )
            do jloc = 1, iloc - 1
               jnod = x_ndgln( jloc )
               if( inod == jnod ) then
                  found = .true.
                  jnod2 = jnod
               end if
            end do
            if( .not. found ) isum = isum + 1    
            ewrite(3,*)'::', iloc, inod, jnod2, isum
         end do
      end if

    endif

      return
    end subroutine Eliminating_Repetitive_Nodes


    subroutine Eliminating_Repetitive_Nodes_all( totele, x_nloc, x_nonods, &
         mx_x_nonods, x_ndgln, x, y, z )
      implicit none
      integer, intent( in ) :: totele, x_nloc, mx_x_nonods
      integer, intent( inout ) :: x_nonods
      integer, dimension( totele * x_nloc ), intent( inout ) :: x_ndgln
      real, dimension( mx_x_nonods ), intent( inout ) :: x, y
      real, dimension( mx_x_nonods ), intent( inout ), optional :: z
      ! Local variables
      real, parameter :: toler = 1.e-10
      integer :: ele, ele2, iloc, jloc, inod, jnod, jnod2, x_nloc2, isum, iref
      integer :: count_nod, jnod_found
      integer, dimension( : ), allocatable :: x_ndgln2, new2old, old2new
      real, dimension( : ), allocatable :: x2, y2, z2

      ewrite(3,*) 'In Eliminating_Repetitive_Nodes_all'
      ewrite(3,*) 'x_nonods, mx_x_nonods:', maxval(x_ndgln), mx_x_nonods

      x_nonods = maxval( x_ndgln ) 
      allocate( new2old( x_nonods ) ) ; new2old = 0
      allocate( old2new( x_nonods ) ) ; old2new = 0
      allocate( x2( x_nonods ) )
      allocate( y2( x_nonods ) )
      allocate( z2( x_nonods ) )

      x2(1:x_nonods)=x(1:x_nonods)
      y2(1:x_nonods)=y(1:x_nonods)
      if (present(z)) then
         z2(1:x_nonods)=z(1:x_nonods)
      else
         z2=0.
      end if

      x=0. ; y=0. 
      if (present(z)) z=0.

      count_nod=0
      do inod=1,x_nonods
         jnod_found=0
         do jnod=1,count_nod
            if( ( abs( x2( inod ) - x( jnod ) ) <= toler ) .and. &
                 ( abs( y2( inod ) - y( jnod ) ) <= toler ) ) then
               if( present( z ) ) then
                  if( ( abs( z2( inod ) - z( jnod ) ) <= toler ) ) jnod_found=jnod
               else
                  jnod_found=jnod
               end if
            end if
         end do
         if(jnod_found==0) then
            count_nod=count_nod+1
            new2old(count_nod)=inod
            old2new(inod)=count_nod
            x(count_nod)=x2(inod)
            y(count_nod)=y2(inod)
            if (present(z)) z(count_nod)=z2(inod)
         else
            old2new(inod)=jnod_found
         end if
      end do

      do ele = 1, totele 
         do iloc = 1, x_nloc
            inod = x_ndgln( ( ele - 1 ) * x_nloc + iloc )
            x_ndgln( ( ele - 1 ) * x_nloc + iloc )=old2new(inod)
         end do
      end do
      x_nonods = count_nod

      return
    end subroutine Eliminating_Repetitive_Nodes_all


    subroutine Adding_Extra_Parametric_Nodes( totele, x_nloc, mx_x_nonods, &
         x_ndgln, x, y )
      implicit none
      integer, intent( in ) :: totele, x_nloc, mx_x_nonods
      integer, dimension( totele * x_nloc ), intent( inout ) :: x_ndgln
      real, dimension( mx_x_nonods ), intent( inout ) :: x, y
      ! Local variables
      integer, dimension( : ), allocatable :: x_ndgln2, loclist
      integer :: ele, iloc, iloc2, x_loc_ref, xnod1, xnod2, xnod3, xnod4, npoly, inod
      integer :: iloc_list(4),jloc_list(4),iiloc,ii,jloc,x_iloc
      real :: rsumx, rsumy

      ewrite(3,*) 'In Adding_Extra_Parametric_Nodes'

      x_loc_ref = maxval( x_ndgln ) 
      ewrite(3,*)' x_loc_ref:', x_loc_ref
      do ele = 1, totele
         ewrite(3,*)'x_ndgln:', ele, &
              ( x_ndgln( ( ele - 1 ) * x_nloc + iloc ) , iloc = 1, x_nloc )
      end do

      if( .true. ) then
         iloc_list(1)=1
         jloc_list(1)=2
         iloc_list(2)=3
         jloc_list(2)=4
         iloc_list(3)=1
         jloc_list(3)=3
         iloc_list(4)=2
         jloc_list(4)=4
         do ele = 1, totele
            iiloc=4

            do ii=1,4
               iloc=iloc_list(ii)
               jloc=jloc_list(ii)
               xnod1 = x_ndgln( ( ele - 1 ) * x_nloc + iloc )
               xnod2 = x_ndgln( ( ele - 1 ) * x_nloc + jloc )
               x_loc_ref = x_loc_ref + 1
               iiloc=iiloc+1
               x_ndgln( ( ele - 1 ) * x_nloc + iiloc ) = x_loc_ref
               x( x_loc_ref ) = 0.5 * ( x( xnod1 ) + x( xnod2 ) ) 
               y( x_loc_ref ) = 0.5 * ( y( xnod1 ) + y( xnod2 ) ) 
            end do

            xnod1 = x_ndgln( ( ele - 1 ) * x_nloc + 1 )
            xnod2 = x_ndgln( ( ele - 1 ) * x_nloc + 2 )
            xnod3 = x_ndgln( ( ele - 1 ) * x_nloc + 3 )
            xnod4 = x_ndgln( ( ele - 1 ) * x_nloc + 4 )
            x_loc_ref = x_loc_ref + 1
            iiloc=iiloc+1
            x_ndgln( ( ele - 1 ) * x_nloc + iiloc ) = x_loc_ref
            x( x_loc_ref ) = 0.25 * ( x( xnod1 ) + x( xnod2 ) + x( xnod3 ) + x( xnod4 )) 
            y( x_loc_ref ) = 0.25 * ( y( xnod1 ) + y( xnod2 ) + y( xnod3 ) + y( xnod4 ) ) 
         end do
      end if ! false

      if( .false. ) then
         do ele = 1, totele
            iloc2 = 4
            do iloc = 1, 4
               iloc2 = iloc2 + 1
               if ( iloc < 4 ) then
                  xnod1 = x_ndgln( ( ele - 1 ) * x_nloc + iloc )
                  xnod2 = x_ndgln( ( ele - 1 ) * x_nloc + iloc + 1 )
                  x_ndgln( ( ele - 1 ) * x_nloc + iloc2 ) = x_loc_ref
               else
                  xnod1 = x_ndgln( ( ele - 1 ) * x_nloc + iloc )
                  xnod2 = x_ndgln( ( ele - 1 ) * x_nloc + 1 )
                  x_ndgln( ( ele - 1 ) * x_nloc + iloc2 ) = x_loc_ref
               end if
               xnod3 = x_loc_ref
               x( xnod3 ) = 0.5 * ( x( xnod1 ) + x( xnod2 ) ) 
               y( xnod3 ) = 0.5 * ( y( xnod1 ) + y( xnod2 ) ) 
               x_loc_ref = x_loc_ref + 1
            end do

            x_ndgln( ( ele - 1 ) * x_nloc + x_nloc ) = x_loc_ref
            xnod4 = x_loc_ref
            rsumx = 0. ; rsumy = 0.
            do iloc = 1, x_nloc - 1
               xnod1 = x_ndgln( ( ele - 1 ) * x_nloc + iloc )
               rsumx = rsumx + 1. / real( x_nloc - 1 ) * x( xnod1 )
               rsumy = rsumy + 1. / real( x_nloc - 1 ) * y( xnod1 )
            end do
            x( xnod4 ) = rsumx 
            y( xnod4 ) = rsumy

            x_loc_ref = x_loc_ref + 1
         end do
      end if

      allocate( x_ndgln2( totele * x_nloc ) ) ; x_ndgln2 = 0
      allocate( loclist( x_nloc ) ) ; loclist = 0
      x_ndgln2 = x_ndgln ; x_ndgln = 0 ; inod = 0 ; jloc = 0

      loclist = (/ 1, 5, 2, 7, 9, 8, 3, 6, 4 /)

      do ele = 1, totele
         do iloc = 1, x_nloc
            jloc = loclist( iloc )
            x_ndgln( ( ele - 1 ) * x_nloc + iloc ) = &
                 x_ndgln2( ( ele - 1 ) * x_nloc + jloc )  
         end do
      end do

      ewrite(3,*) '-'
      ewrite(3,*) ''

      do ele = 1, totele
         ewrite(3,*)'x_ndgln:', ele, &
              ( x_ndgln( ( ele - 1 ) * x_nloc + iloc ) , iloc = 1, x_nloc )
      end do

      deallocate( x_ndgln2 )
      deallocate( loclist )

      return
    end subroutine Adding_Extra_Parametric_Nodes

    subroutine Adding_Extra_Parametric_Nodes_Tet( totele, x_nloc, mx_x_nonods, &
         x_ndgln, x, y, z )
      implicit none
      integer, intent( in ) :: totele, x_nloc, mx_x_nonods
      integer, dimension( totele * x_nloc ), intent( inout ) :: x_ndgln
      real, dimension( mx_x_nonods ), intent( inout ) :: x, y, z
      ! Local variables
      integer, dimension( : ), allocatable :: x_ndgln2, loclist
      integer :: ele, iloc, iloc2, x_loc_ref
      integer :: xnod1, xnod2, xnod3, xnod4
      integer :: xnod5, xnod6, xnod7, xnod8 
      integer :: npoly, inod
      integer :: iloc_list(4),jloc_list(4),iiloc,jloc, ii, jj
      real :: rsumx, rsumy

      ewrite(3,*) 'In Adding_Extra_Parametric_Nodes'

      x_loc_ref = maxval( x_ndgln ) 

      if( .true. ) then
         iloc_list(1)=1
         jloc_list(1)=1
         iloc_list(2)=1
         jloc_list(2)=1
         iloc_list(3)=1
         jloc_list(3)=1
         iloc_list(4)=1
         jloc_list(4)=1

         do ele = 1, totele ! loop over hexes
            iiloc=8

            do ii = 1, 6 ! loop over faces

               do jj = 1, 4 ! loop over edges

                  iloc=1
                  jloc=1
                  xnod1 = x_ndgln( ( ele - 1 ) * x_nloc + iloc )
                  xnod2 = x_ndgln( ( ele - 1 ) * x_nloc + jloc )

                  x_loc_ref = x_loc_ref + 1
                  iiloc=iiloc+1

                  ! add a parametric node on the edge
                  x_ndgln( ( ele - 1 ) * x_nloc + iiloc ) = x_loc_ref
                  x( x_loc_ref ) = 0.5 * ( x( xnod1 ) + x( xnod2 ) ) 
                  y( x_loc_ref ) = 0.5 * ( y( xnod1 ) + y( xnod2 ) ) 
                  z( x_loc_ref ) = 0.5 * ( z( xnod1 ) + z( xnod2 ) ) 
               end do

               xnod1 = x_ndgln( ( ele - 1 ) * x_nloc + 1 ) ! this 1 2 3 4 is wrong
               xnod2 = x_ndgln( ( ele - 1 ) * x_nloc + 2 )
               xnod3 = x_ndgln( ( ele - 1 ) * x_nloc + 3 )
               xnod4 = x_ndgln( ( ele - 1 ) * x_nloc + 4 )

               x_loc_ref = x_loc_ref + 1
               iiloc=iiloc+1

               ! add a parametric node on the face
               x_ndgln( ( ele - 1 ) * x_nloc + iiloc ) = x_loc_ref
               x( x_loc_ref ) = 0.25 * ( x( xnod1 ) + x( xnod2 ) + x( xnod3 ) + x( xnod4 ) ) 
               y( x_loc_ref ) = 0.25 * ( y( xnod1 ) + y( xnod2 ) + y( xnod3 ) + y( xnod4 ) )
               z( x_loc_ref ) = 0.25 * ( z( xnod1 ) + z( xnod2 ) + z( xnod3 ) + z( xnod4 ) )
            end do ! loop over faces
            
            xnod1 = x_ndgln( ( ele - 1 ) * x_nloc + 1 ) ! this 1 2 3 4... is wrong
            xnod2 = x_ndgln( ( ele - 1 ) * x_nloc + 2 )
            xnod3 = x_ndgln( ( ele - 1 ) * x_nloc + 3 )
            xnod4 = x_ndgln( ( ele - 1 ) * x_nloc + 4 )
            xnod5 = x_ndgln( ( ele - 1 ) * x_nloc + 5 )
            xnod6 = x_ndgln( ( ele - 1 ) * x_nloc + 6 )
            xnod7 = x_ndgln( ( ele - 1 ) * x_nloc + 7 )
            xnod8 = x_ndgln( ( ele - 1 ) * x_nloc + 8 )

            x_loc_ref = x_loc_ref + 1
            iiloc=iiloc+1

            ! add a parametric node on the volume 
            x_ndgln( ( ele - 1 ) * x_nloc + iiloc ) = x_loc_ref
            x( x_loc_ref ) = 0.125 * ( x( xnod1 ) + x( xnod2 ) + x( xnod3 ) + x( xnod4 ) + &
                 x( xnod5 ) + x( xnod6 ) + x( xnod7 ) + x( xnod8 ) )
            y( x_loc_ref ) = 0.125 * ( y( xnod1 ) + y( xnod2 ) + y( xnod3 ) + y( xnod4 ) + &
                 y( xnod5 ) + y( xnod6 ) + y( xnod7 ) + y( xnod8 ) )
            z( x_loc_ref ) = 0.125 * ( z( xnod1 ) + z( xnod2 ) + z( xnod3 ) + z( xnod4 ) + &
                 z( xnod5 ) + z( xnod6 ) + z( xnod7 ) + z( xnod8 ) )

         end do ! loop over hexes

      end if

      allocate( x_ndgln2( totele * x_nloc ) ) ; x_ndgln2 = 0
      allocate( loclist( x_nloc ) ) ; loclist = 0
      x_ndgln2 = x_ndgln ; x_ndgln = 0 ; inod = 0 ; jloc = 0

      loclist = (/ 1, 5, 2, 6, 7, 8, 3, 9, 4 /) ! this is wrong

      do ele = 1, totele
         do iloc = 1, x_nloc
            jloc = loclist( iloc )
            x_ndgln( ( ele - 1 ) * x_nloc + iloc ) = &
                 x_ndgln2( ( ele - 1 ) * x_nloc + jloc )  
         end do
      end do

      do ele = 1, totele
         ewrite(3,*)'x_ndgln:', ele, &
              ( x_ndgln( ( ele - 1 ) * x_nloc + iloc ) , iloc = 1, x_nloc )
      end do

      deallocate( x_ndgln2 )
      deallocate( loclist )

      return
    end subroutine Adding_Extra_Parametric_Nodes_Tet


!!!!!!!!
!!!!!!!! Making and Numbering Quadratic Tetrahedron and 8 Hexahedra 
!!!!!!!!


    subroutine Make_QTets( totele, quad_cv_nloc, x_nloc, max_x_nonods, x_nonods, &
         x_ndgln_real, lx, ly, lz, x, y, z, fem_nod, &
         xp2, yp2, zp2, x_ndgln_p2 )
      ! This subrt creates the local coordinates and node points for:
      ! (a) quadratic tetrahedra of unit volume and (b) 27 points of
      ! the 8 hexahedra  within the 8 linear tetrahedra.
      ! FEM_NOD is the local numbering of the FEM representation of
      ! the unit volume quadratic tetrahedron.
      implicit none
      integer, intent( inout ) :: totele, x_nonods
      integer, intent( in ) :: quad_cv_nloc, x_nloc, max_x_nonods
      integer, dimension( max_x_nonods ), intent( inout ) :: x_ndgln_real
      real, dimension( max_x_nonods ), intent( inout ) :: lx, ly, lz, x, y, z
      integer, dimension( max_x_nonods ), intent( inout ) :: fem_nod
      real, dimension( max_x_nonods ), intent( inout ) :: xp2, yp2, zp2
      integer, dimension( max_x_nonods ), intent( inout ) :: x_ndgln_p2
      ! Local variables
      real, parameter :: h_scale = 2.0396489026555056
      integer, parameter :: number_of_hexs = 4,  number_of_nodes = 27
      integer :: ele, iloc, istart, ifinish
      integer, dimension( : ), allocatable :: x_ndgln, iloclist
      real, dimension( : ), allocatable :: Volume_P1
      real :: Volume_P2, Volume_P1_Tets

      ewrite(3,*)' In Make_QTets'

!!!
!!! Setting up unity area quadratic tetrahedron
!!!
      ! Level 2 (basis of the tetrahedron )
      lx( 1 ) = 0.
      ly( 1 ) = 0.
      lz( 1 ) = 0.

      lx( 2 ) = 1. * h_scale
      ly( 2 ) = 0.
      lz( 2 ) = 0.

      lx( 3 ) = .5 * h_scale
      ly( 3 ) = sqrt( 3. / 4. ) * h_scale
      lz( 3 ) = 0.

      lx( 4 ) = 0.5 * ( lx( 1 ) + lx( 2 ) )
      ly( 4 ) = 0.5 * ( ly( 1 ) + ly( 3 ) )
      lz( 4 ) = sqrt( 2. / 3. ) * h_scale

      ewrite(3,*)'lx:', lx( 1 : 4 )
      ewrite(3,*)'ly:', ly( 1 : 4 )
      ewrite(3,*)'lz:', lz( 1 : 4 )

      ! Computing X / Y / Z and FEM_NOD for the tetrahedra
      xp2 = 0. ; yp2 = 0. ; zp2 = 0.

      xp2( 1 ) = 0.
      yp2( 1 ) = 0.
      zp2( 1 ) = 0.

      xp2( 3 ) = 1. * h_scale
      yp2( 3 ) = 0.
      zp2( 3 ) = 0.

      xp2( 2 ) = 0.5 * ( xp2( 1 ) + xp2( 3 ) )
      yp2( 2 ) = 0.5 * ( yp2( 1 ) + yp2( 3 ) )
      zp2( 2 ) = 0.5 * ( zp2( 1 ) + zp2( 3 ) )

      xp2( 6 ) = 0.5 * h_scale
      yp2( 6 ) = sqrt( 3. / 4. ) * h_scale
      zp2( 6 ) = 0.

      xp2( 5 ) = 0.5 * ( xp2( 3 ) + xp2( 6 ) )
      yp2( 5 ) = 0.5 * ( yp2( 3 ) + yp2( 6 ) )
      zp2( 5 ) = 0.5 * ( zp2( 3 ) + zp2( 6 ) ) 

      xp2( 4 ) = 0.5 * ( xp2( 1 ) + xp2( 6 ) )
      yp2( 4 ) = 0.5 * ( yp2( 1 ) + yp2( 6 ) )
      zp2( 4 ) = 0.5 * ( zp2( 1 ) + zp2( 6 ) )

      ! Level 0 (top of the tetrahedron )
      xp2( 10 ) = 0.5 * ( xp2( 1 ) + xp2( 3 ) )
      yp2( 10 ) = 0.5 * ( yp2( 1 ) + yp2( 6 ) )
      zp2( 10 ) = sqrt( 2. / 3. ) * h_scale

      ! Level 1 (Intermediate level)
      xp2( 7 ) = 0.5 * ( xp2( 1 ) + xp2( 10 ) )
      yp2( 7 ) = 0.5 * ( yp2( 1 ) + yp2( 10 ) )
      zp2( 7 ) = 0.5 * ( zp2( 1 ) + zp2( 10 ) )

      xp2( 8 ) = 0.5 * ( xp2( 3 ) + xp2( 10 ) )
      yp2( 8 ) = 0.5 * ( yp2( 3 ) + yp2( 10 ) )
      zp2( 8 ) = 0.5 * ( zp2( 3 ) + zp2( 10 ) )

      xp2( 9 ) = 0.5 * ( xp2( 6 ) + xp2( 10 ) )
      yp2( 9 ) = 0.5 * ( yp2( 6 ) + yp2( 10 ) )
      zp2( 9 ) = 0.5 * ( zp2( 6 ) + zp2( 10 ) )

      Volume_P2 = Volume_TetHex( .false., &
           xp2( 1 ), xp2( 3 ), xp2( 6 ), xp2( 10 ), &
           yp2( 1 ), yp2( 3 ), yp2( 6 ), yp2( 10 ), &
           zp2( 1 ), zp2( 3 ), zp2( 6 ), zp2( 10 ) )
      ewrite(3,*)'Volume of P2 Tets:', Volume_P2

      ! Defining FEM_NODs
      do iloc = 1, 10
         fem_nod( iloc ) = iloc
      end do

      ewrite(3,*)' X/Y/Z for the P2 Tetrahedron:'
      do iloc = 1, 10
         ewrite(3,*) iloc, xp2( iloc ), yp2( iloc ), zp2( iloc )
      end do

!!!
!!! Defining linear tetrahedra (8 within the quadratic tetrahedron)
!!!
      allocate( x_ndgln( max_x_nonods ) ) ; x_ndgln = 0
      !allocate( x_ndgln_p2( max_x_nonods ) ) ; x_ndgln_p2 = 0
      allocate( Volume_P1( 8 ) ) ; Volume_P1 = 0.

      ele = 1
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 1 ) = 7
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 2 ) = 8
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 3 ) = 9
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 4 ) = 10

      Volume_P1( ele ) = Volume_TetHex( .false., &
           xp2( 7 ), xp2( 8 ), xp2( 9 ), xp2( 10 ), &
           yp2( 7 ), yp2( 8 ), yp2( 9 ), yp2( 10 ), &
           zp2( 7 ), zp2( 8 ), zp2( 9 ), zp2( 10 ) )
      ewrite(3,*)'Tet, Volume of P1 Tets:', ele, Volume_P1( ele )

      ele = 2
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 1 ) = 1
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 2 ) = 2
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 3 ) = 4
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 4 ) = 7

      Volume_P1( ele ) = Volume_TetHex( .false., &
           xp2( 1 ), xp2( 2 ), xp2( 4 ), xp2( 7 ), &
           yp2( 1 ), yp2( 2 ), yp2( 4 ), yp2( 7 ), &
           zp2( 1 ), zp2( 2 ), zp2( 4 ), zp2( 7 ) )
      ewrite(3,*)'Tet, Volume of P1 Tets:', ele, Volume_P1( ele )

      ele = 3
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 1 ) = 2
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 2 ) = 7
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 3 ) = 8
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 4 ) = 4

      Volume_P1( ele ) = Volume_TetHex( .false., &
           xp2( 2 ), xp2( 7 ), xp2( 8 ), xp2( 4 ), &
           yp2( 2 ), yp2( 7 ), yp2( 8 ), yp2( 4 ), &
           zp2( 2 ), zp2( 7 ), zp2( 8 ), zp2( 4 ) )
      ewrite(3,*)'Tet, Volume of P1 Tets:', ele, Volume_P1( ele ) 


      ele = 4
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 1 ) = 2
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 2 ) = 3
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 3 ) = 4
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 4 ) = 8

      Volume_P1( ele ) = Volume_TetHex( .false., &
           xp2( 2 ), xp2( 3 ), xp2( 4 ), xp2( 8 ), &
           yp2( 2 ), yp2( 3 ), yp2( 4 ), yp2( 8 ), &
           zp2( 2 ), zp2( 3 ), zp2( 4 ), zp2( 8 ) )
      ewrite(3,*)'Tet, Volume of P1 Tets:', ele, Volume_P1( ele )

      ele = 5
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 1 ) = 3
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 2 ) = 5
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 3 ) = 4
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 4 ) = 8

      Volume_P1( ele ) = Volume_TetHex( .false., &
           xp2( 3 ), xp2( 5 ), xp2( 4 ), xp2( 8 ), &
           yp2( 3 ), yp2( 5 ), yp2( 4 ), yp2( 8 ), &
           zp2( 3 ), zp2( 5 ), zp2( 4 ), zp2( 8 ) )
      ewrite(3,*)'Tet, Volume of P1 Tets:', ele, Volume_P1( ele )

      ele = 6
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 1 ) = 4
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 2 ) = 5
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 3 ) = 9
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 4 ) = 8

      Volume_P1( ele ) = Volume_TetHex( .false., &
           xp2( 4 ), xp2( 5 ), xp2( 9 ), xp2( 8 ), &
           yp2( 4 ), yp2( 5 ), yp2( 9 ), yp2( 8 ), &
           zp2( 4 ), zp2( 5 ), zp2( 9 ), zp2( 8 ) )
      ewrite(3,*)'Tet, Volume of P1 Tets:', ele, Volume_P1( ele )


      ele = 7
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 1 ) = 5
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 2 ) = 6
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 3 ) = 4
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 4 ) = 9

      Volume_P1( ele ) = Volume_TetHex( .false., &
           xp2( 5 ), xp2( 6 ), xp2( 4 ), xp2( 9 ), &
           yp2( 5 ), yp2( 6 ), yp2( 4 ), yp2( 9 ), &
           zp2( 5 ), zp2( 6 ), zp2( 4 ), zp2( 9 ) )
      ewrite(3,*)'Tet, Volume of P1 Tets:', ele, Volume_P1( ele )

      ele = 8
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 1 ) = 7
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 2 ) = 9
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 3 ) = 8
      x_ndgln_p2( ( ele - 1 ) * x_nloc + 4 ) = 4

      Volume_P1( ele ) = Volume_TetHex( .false., &
           xp2( 7 ), xp2( 9 ), xp2( 8 ), xp2( 4 ), &
           yp2( 7 ), yp2( 9 ), yp2( 8 ), yp2( 4 ), &
           zp2( 7 ), zp2( 9 ), zp2( 8 ), zp2( 4 ) )
      ewrite(3,*)'Tet, Volume of P1 Tets:', ele, Volume_P1( ele )

      Volume_P1_Tets = 0.
      do ele = 1, totele
         Volume_P1_Tets = Volume_P1_Tets + Volume_P1( ele ) 
      end do

      ewrite(3,*)' Total Volume of P1 Tets:', Volume_P1_Tets

      if( abs( Volume_P1_Tets - Volume_P2 ) >= 1.e-7 ) then
           FLAbort( "Volumes of P2 and the sum of 8 P1s dont match " )
        end if

      do ele = 1, totele
         call Make_Linear_Tetrahedron( ele, quad_cv_nloc, x_nloc, x_nonods, &
              number_of_hexs, &
              xp2, yp2, zp2, x, y, z, &
              x_ndgln_p2, x_ndgln )
      end do

      call Make_BiLinear_Hexahedra( totele, number_of_hexs, quad_cv_nloc, x_nonods, &
           x, y, z, x_ndgln )

      ! Now building the final X_NDGLN with the numbering consistent with the 
      ! remanining of the model

      allocate( iloclist( quad_cv_nloc ) ) ; iloclist = 0
      iloclist = &
           (/  5, 14,  6, 17, 18, 15,  7, 16,  8, &  ! Level 1
           19, 23, 20, 26, 27, 24, 21, 25, 22, & ! Level 2
           1,  9,  2, 12, 13, 10,  3, 11,  4  /)       ! Level 3

      totele = totele * number_of_hexs ! From now, TOTELE is actual total number of hexahedra
      do ele = 1, totele
         do iloc = 1, quad_cv_nloc
            x_ndgln_real( ( ele - 1 ) * quad_cv_nloc + iloc ) = &
                 x_ndgln( ( ele - 1 ) * quad_cv_nloc + iloclist( iloc ) )
         end do
      end do
      x_nonods = maxval( x_ndgln_real )
      ewrite(3,*) 'just before leaving Make_QTets'
      ewrite(3,*) 'x_nonods:', x_nonods

      ele = 30+2
      ewrite(3,*) x_ndgln_real( ( ele - 1 ) * quad_cv_nloc + 1 :  ( ele - 1 ) * quad_cv_nloc + 27 )
      do iloc = 1, quad_cv_nloc
         ewrite(3,*) x_ndgln_real( ( ele - 1 ) * quad_cv_nloc + iloc ), &
              x(x_ndgln_real( ( ele - 1 ) * quad_cv_nloc + iloc )), &
              y(x_ndgln_real( ( ele - 1 ) * quad_cv_nloc + iloc )), &
              z(x_ndgln_real( ( ele - 1 ) * quad_cv_nloc + iloc ))
      end do
      !stop 999

      deallocate( x_ndgln )
      deallocate( iloclist )

      return
    end subroutine Make_QTets


    subroutine Make_Linear_Tetrahedron( ele, quad_cv_nloc, x_nloc, x_nonods, &
         number_of_hexs, &
         xp2, yp2, zp2, x, y, z, &
         x_ndgln_p2, x_ndgln )
      implicit none
      integer, intent( in ) :: ele, quad_cv_nloc, x_nloc, x_nonods, &
           &                           number_of_hexs
      real, dimension( 10 ), intent( in ) :: xp2, yp2, zp2
      real, dimension( x_nonods ), intent( inout ) :: x, y, z
      integer, dimension( x_nonods ), intent( in ) :: x_ndgln_p2
      integer, dimension( x_nonods ), intent( inout ) :: x_ndgln

      ! Local variables
      real, dimension( 4 ) :: lx, ly, lz
      integer :: ele_hex, hex_numbering, jloc, iloc, kloc, piloc
      integer, dimension( : ), allocatable :: inod

      ewrite(3,*)' In Make_Linear_Tetrahedron'

      do iloc = 1, x_nloc
         lx( iloc ) = xp2( x_ndgln_p2( ( ele - 1 ) * x_nloc + iloc ) )
         ly( iloc ) = yp2( x_ndgln_p2( ( ele - 1 ) * x_nloc + iloc ) )
         lz( iloc ) = zp2( x_ndgln_p2( ( ele - 1 ) * x_nloc + iloc ) )

         ewrite(3,*) iloc, x_ndgln_p2( ( ele - 1 ) * x_nloc + iloc ) 
      end do

      ewrite(3,*)'lx:', lx( 1 : 4 )
      ewrite(3,*)'ly:', ly( 1 : 4 )
      ewrite(3,*)'lz:', lz( 1 : 4 )

      ! Creating the numbering for the P1 tets and hexs
      allocate( inod( 11 ) ) ; inod = 0
      do iloc = 1, 11
         inod( iloc ) = ( ele - 1 ) * quad_cv_nloc * number_of_hexs + 10 + iloc
      end do

      do iloc = 1, 11
         ewrite(3,*)'x_ndgln(), inod:',  &
              ( ele - 1 ) * quad_cv_nloc * number_of_hexs +iloc, &
              inod( iloc )
      end do

!!!
!!! Remmaping
!!!

      ! Node point 1
      !x( inod + 1 ) = lx( 1 )
      !y( inod + 1 ) = ly( 1 )
      ! Node point 2
      x( inod( 1 ) ) = 0.5 * ( lx( 1 ) + lx( 2 ) )
      y( inod( 1 ) ) = 0.5 * ( ly( 1 ) + ly( 2 ) )
      z( inod( 1 ) ) = 0.5 * ( lz( 1 ) + lz( 2 ) )
      ! Node point 3
      !x( inod + 3 ) = lx( 2 )
      !y( inod + 3 ) = ly( 2 )
      ! Node point 4
      x( inod( 2 ) ) = 0.5 * ( lx( 3 ) + lx( 1 ) )
      y( inod( 2 ) ) = 0.5 * ( ly( 3 ) + ly( 1 ) )
      z( inod( 2 ) ) = 0.5 * ( lz( 3 ) + lz( 1 ) )
      ! Node point 5
      x( inod( 3 ) ) = 1./3. * ( lx( 1 ) + lx( 2 ) + lx( 3 ) )
      y( inod( 3 ) ) = 1./3. * ( ly( 1 ) + ly( 2 ) + ly( 3 ) )
      z( inod( 3 ) ) = 1./3. * ( lz( 1 ) + lz( 2 ) + lz( 3 ) )
      ! Node point 6
      x( inod( 4 ) ) = 0.5 * ( lx( 2 ) + lx( 3 ) )
      y( inod( 4 ) ) = 0.5 * ( ly( 2 ) + ly( 3 ) )
      z( inod( 4 ) ) = 0.5 * ( lz( 2 ) + lz( 3 ) )
      ! Node point 7
      !x( inod + 7 ) = lx( 3 )
      !y( inod + 7 ) = ly( 3 )
      ! Node point 15
      !x( inod + 15 ) = lx( 4 )
      !y( inod + 15 ) = ly( 4 )
      ! Node point 8
      x( inod( 5 ) ) = 1./3. * ( lx( 1 ) + lx( 2 ) + lx( 4 ) )
      y( inod( 5 ) ) = 1./3. * ( ly( 1 ) + ly( 2 ) + ly( 4 ) )
      z( inod( 5 ) ) = 1./3. * ( lz( 1 ) + lz( 2 ) + lz( 4 ) )
      ! Node point 9
      x( inod( 6 ) ) = 0.5 * ( lx( 1 ) + lx( 4 ) )
      y( inod( 6 ) ) = 0.5 * ( ly( 1 ) + ly( 4 ) )
      z( inod( 6 ) ) = 0.5 * ( lz( 1 ) + lz( 4 ) )
      ! Node point 10
      x( inod( 7 ) ) = 0.5 * ( lx( 2 ) + lx( 4 ) )
      y( inod( 7 ) ) = 0.5 * ( ly( 2 ) + ly( 4 ) )
      z( inod( 7 ) ) = 0.5 * ( lz( 2 ) + lz( 4 ) )
      ! Node point 11
      x( inod( 8 ) ) = 0.25 * ( lx( 1 ) + lx( 2 ) + lx( 3 ) + lx( 4 ) )
      y( inod( 8 ) ) = 0.25 * ( ly( 1 ) + ly( 2 ) + ly( 3 ) + ly( 4 ) )
      z( inod( 8 ) ) = 0.25 * ( lz( 1 ) + lz( 2 ) + lz( 3 ) + lz( 4 ) )
      ! Node point 12
      x( inod( 9 ) ) = 1./3. * ( lx( 1 ) + lx( 3 ) + lx( 4 ) )
      y( inod( 9 ) ) = 1./3. * ( ly( 1 ) + ly( 3 ) + ly( 4 ) )
      z( inod( 9 ) ) = 1./3. * ( lz( 1 ) + lz( 3 ) + lz( 4 ) )
      ! Node point 13
      x( inod( 10 ) ) = 1./3. * ( lx( 2 ) + lx( 3 ) + lx( 4 ) )
      y( inod( 10 ) ) = 1./3. * ( ly( 2 ) + ly( 3 ) + ly( 4 ) )
      z( inod( 10 ) ) = 1./3. * ( lz( 2 ) + lz( 3 ) + lz( 4 ) )
      ! Node point 14
      x( inod( 11 ) ) = 0.5 * ( lx( 3 ) + lx( 4 ) )
      y( inod( 11 ) ) = 0.5 * ( ly( 3 ) + ly( 4 ) )
      z( inod( 11 ) ) = 0.5 * ( lz( 3 ) + lz( 4 ) )

!!!
!!! Computing pseudo elements
!!!

      hex_numbering = ( ele - 1 ) * number_of_hexs

      ele_hex = hex_numbering + 1 ! Hex 1
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 1 ) = inod( 6 ) ! Node 9
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 2 ) = inod( 5 ) ! Node 8
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 3 ) = inod( 9 ) ! Node 12
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 4 ) = inod( 8 ) ! Node 11
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 5 ) = &
           x_ndgln_p2( ( ele - 1 ) * x_nloc + 1 )                          ! Node 1
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 6 ) = inod( 1 ) ! Node 2
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 7 ) = inod( 2 ) ! Node 4
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 8 ) = inod( 3 ) ! Node 5

      x( x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 5 )  ) = &
           xp2( x_ndgln_p2( ( ele - 1 ) * x_nloc + 1 )   )
      y( x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 5 )  ) = &
           yp2( x_ndgln_p2( ( ele - 1 ) * x_nloc + 1 )   )
      z( x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 5 )  ) = &
           zp2( x_ndgln_p2( ( ele - 1 ) * x_nloc + 1 )   )
      
      ele_hex = hex_numbering + 2 ! Hex 2
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 1 ) = inod( 5 )  ! Node 8
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 2 ) = inod( 7 )  ! Node 10
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 3 ) = inod( 8 )  ! Node 11
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 4 ) = inod( 10 )! Node 13
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 5 ) = inod( 1 )  ! Node 2
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 6 ) = &
           x_ndgln_p2( ( ele - 1 ) * x_nloc + 2 )                           ! Node 3
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 7 ) = inod( 3 )  ! Node 5
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 8 ) = inod( 4 )  ! Node 6

      x( x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 6 )  ) = &
           xp2( x_ndgln_p2( ( ele - 1 ) * x_nloc + 2 )   )
      y( x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 6 )  ) = &
           yp2( x_ndgln_p2( ( ele - 1 ) * x_nloc + 2 )   )
      z( x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 6 )  ) = &
           zp2( x_ndgln_p2( ( ele - 1 ) * x_nloc + 2 )   )

      ele_hex = hex_numbering + 3 ! Hex 3
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 1 ) = inod( 9 )   ! Node 12
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 2 ) = inod( 8 )   ! Node 11
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 3 ) = inod( 11 ) ! Node 14
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 4 ) = inod( 10 ) ! Node 13
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 5 ) = inod( 2 )   ! Node 4
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 6 ) = inod( 3 )   ! Node 5
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 7 ) = &
           x_ndgln_p2( ( ele - 1 ) * x_nloc + 3 )                            ! Node 7
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 8 ) = inod( 4 )   ! Node 6

      x( x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 7 )  ) = &
           xp2( x_ndgln_p2( ( ele - 1 ) * x_nloc + 3 )   )
      y( x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 7 )  ) = &
           yp2( x_ndgln_p2( ( ele - 1 ) * x_nloc + 3 )   )
      z( x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 7 )  ) = &
           zp2( x_ndgln_p2( ( ele - 1 ) * x_nloc + 3 )   )

      ele_hex = hex_numbering + 4 ! Hex 4
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 1 ) = &
           x_ndgln_p2( ( ele - 1 ) * x_nloc + 4 )                            ! Node 15
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 2 ) = inod( 7 )   ! Node 10
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 3 ) = inod( 11 ) ! Node 14
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 4 ) = inod( 10 ) ! Node 13
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 5 ) = inod( 6 )   ! Node 9
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 6 ) = inod( 5 )   ! Node 8
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 7 ) = inod( 9 )   ! Node 12
      x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 8 ) = inod( 8 )   ! Node 11

      x( x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 1 )  ) = &
           xp2( x_ndgln_p2( ( ele - 1 ) * x_nloc + 4 )   )
      y( x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 1 )  ) = &
           yp2( x_ndgln_p2( ( ele - 1 ) * x_nloc + 4 )   )
      z( x_ndgln( ( ele_hex - 1 ) * quad_cv_nloc + 1 )  ) = &
           zp2( x_ndgln_p2( ( ele - 1 ) * x_nloc + 4 )   )

      ewrite(3,*)'Checking the initial numbering - X_NDGLN, for the hexs:'
      ewrite(3,*) ' ele, ele_hex, xndgln()=, xndgln:'
      do ele_hex = 1, 4
         kloc = ( ele - 1 ) * number_of_hexs + ele_hex
         do iloc = 1, 8
            ewrite(3,*) ele, ele_hex, ( kloc - 1 ) * quad_cv_nloc + iloc, &
                 x_ndgln( ( kloc - 1 ) * quad_cv_nloc + iloc ), &
                 x( x_ndgln(( kloc - 1 ) * quad_cv_nloc + iloc) ), &
                 y( x_ndgln(( kloc - 1 ) * quad_cv_nloc + iloc) ), &
                 z( x_ndgln(( kloc - 1 ) * quad_cv_nloc + iloc) )
         end do
      end do

      deallocate( inod )

      return
    end subroutine Make_Linear_Tetrahedron


  subroutine Make_BiLinear_Hexahedra( totele, number_of_hexs, quad_cv_nloc, x_nonods, &
       x, y, z, x_ndgln )
    implicit none
    integer, intent( in ) :: totele, number_of_hexs, quad_cv_nloc
    integer, intent( inout ) :: x_nonods
    real, dimension( x_nonods ), intent( inout ) :: x, y, z
    integer, dimension( x_nonods ), intent( inout ) :: x_ndgln
    ! Local variables
    integer :: ele, ele_hex, iloc, kloc, ele_hex2

    Loop_Tets: do ele = 1, totele
       Loop_Hexs: do ele_hex = 1, number_of_hexs
          call Adding_Parametric_Nodes_Hex( ele, ele_hex, totele, number_of_hexs, &
               quad_cv_nloc, x_nonods, &
               x, y, z, x_ndgln )
       end do Loop_Hexs
    end do Loop_Tets

    call Eliminating_Repetitive_Nodes_all( totele * number_of_hexs, quad_cv_nloc, x_nonods, &
         x_nonods, x_ndgln, x, y, z )

    return
  end subroutine Make_BiLinear_Hexahedra
 

  subroutine Adding_Parametric_Nodes_Hex( ele, ele_hex, totele, number_of_hexs, &
       quad_cv_nloc, x_nonods, &
       x, y, z, x_ndgln )
    implicit none
    integer, intent( in ) :: ele, ele_hex, totele, number_of_hexs, quad_cv_nloc, x_nonods
    real, dimension( x_nonods ), intent( inout ) :: x, y, z
    integer, dimension( x_nonods ), intent( inout ) :: x_ndgln
    ! Local variables
    integer, parameter :: jloc2 = 8
    integer :: inod2, jnod, knod, iloc, ele_hex2, kloc
    integer, dimension( : ), allocatable :: inod

    ewrite(3,*)' In Adding_Parametric_Nodes_Hex '

    ele_hex2 = ( ele - 1 ) * number_of_hexs  + ele_hex

    ewrite(3,*)' Hexahedron', ele_hex, 'is based on the following nodes:'
    do iloc = 1, 8
       ewrite(3,*) ( ele_hex2 - 1 ) * quad_cv_nloc + iloc , &
            x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + iloc )
    end do

    ewrite(3,*)' Master node points for the Hex:'
    do iloc = 1, 8
       jnod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + iloc )
       ewrite(3,*) jnod, x( jnod ), y( jnod ), z( jnod )
    end do

    inod2 = maxval( x_ndgln )

    allocate( inod ( 19 ) ) ; inod = 0
    do iloc = 1, 19 ! Extra parametric nodes 
       x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + jloc2 + iloc ) = inod2 + iloc
       inod( iloc ) = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + jloc2 + iloc )
       ewrite(3,*)'inod:', iloc, ( ele_hex2 - 1 ) * quad_cv_nloc + jloc2 + iloc, inod( iloc )
    end do

    ewrite(3,*)'Numbering / Address for the Hex with parametric nodes'
    ewrite(3,*) x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + jloc2 + 1 : &
         ( ele_hex2 - 1 ) * quad_cv_nloc + jloc2 + 19 )

!!!
!!! Level 1
!!!
    ! Node 9
    jnod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 1 )
    knod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 2 )

    x( inod( 1 ) ) = 0.5 * ( x( jnod ) + x( knod ) )
    y( inod( 1 ) ) = 0.5 * ( y( jnod ) + y( knod ) )
    z( inod( 1 ) ) = 0.5 * ( z( jnod ) + z( knod ) )

    ! Node 10
    jnod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 2 )
    knod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 4 )

    x( inod( 2 ) ) = 0.5 * ( x( jnod ) + x( knod ) )
    y( inod( 2 ) ) = 0.5 * ( y( jnod ) + y( knod ) )
    z( inod( 2 ) ) = 0.5 * ( z( jnod ) + z( knod ) )

    ! Node 11
    jnod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 3 )
    knod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 4 )

    x( inod( 3 ) ) = 0.5 * ( x( jnod ) + x( knod ) )
    y( inod( 3 ) ) = 0.5 * ( y( jnod ) + y( knod ) )
    z( inod( 3 ) ) = 0.5 * ( z( jnod ) + z( knod ) )

    ! Node 12
    jnod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 1 )
    knod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 3 )

    x( inod( 4 ) ) = 0.5 * ( x( jnod ) + x( knod ) )
    y( inod( 4 ) ) = 0.5 * ( y( jnod ) + y( knod ) )
    z( inod( 4 ) ) = 0.5 * ( z( jnod ) + z( knod ) )

    ! Node 13
    x( inod( 5 ) ) = 0.25 * ( sum( x ( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 1  : &
         ( ele_hex2 - 1 ) * quad_cv_nloc + 4 ) ) ) )
    y( inod( 5 ) ) = 0.25 * ( sum( y ( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 1  : &
         ( ele_hex2 - 1 ) * quad_cv_nloc + 4 ) ) ) )
    z( inod( 5 ) ) = 0.25 * ( sum( z ( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 1  : &
         ( ele_hex2 - 1 ) * quad_cv_nloc + 4 ) ) ) )

!!!
!!! Level 3
!!!
    ! Node 14
    jnod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 5 )
    knod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 6 )

    x( inod( 6 ) ) = 0.5 * ( x( jnod ) + x( knod ) )
    y( inod( 6 ) ) = 0.5 * ( y( jnod ) + y( knod ) )
    z( inod( 6 ) ) = 0.5 * ( z( jnod ) + z( knod ) )

    ! Node 15
    jnod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 6 )
    knod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 8 )

    x( inod( 7 ) ) = 0.5 * ( x( jnod ) + x( knod ) )
    y( inod( 7 ) ) = 0.5 * ( y( jnod ) + y( knod ) )
    z( inod( 7 ) ) = 0.5 * ( z( jnod ) + z( knod ) )

    ! Node 16
    jnod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 7 )
    knod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 8 )

    x( inod( 8 ) ) = 0.5 * ( x( jnod ) + x( knod ) )
    y( inod( 8 ) ) = 0.5 * ( y( jnod ) + y( knod ) )
    z( inod( 8 ) ) = 0.5 * ( z( jnod ) + z( knod ) )

    ! Node 17
    jnod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 5 )
    knod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 7 )

    x( inod( 9 ) ) = 0.5 * ( x( jnod ) + x( knod ) )
    y( inod( 9 ) ) = 0.5 * ( y( jnod ) + y( knod ) )
    z( inod( 9 ) ) = 0.5 * ( z( jnod ) + z( knod ) )

    ! Node 18
    x( inod( 10 ) ) = 0.25 * ( sum( x ( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 5  : &
         ( ele_hex2 - 1 ) * quad_cv_nloc + 8 ) ) ) )
    y( inod( 10 ) ) = 0.25 * ( sum( y ( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 5  : &
         ( ele_hex2 - 1 ) * quad_cv_nloc + 8 ) ) ) )
    z( inod( 10 ) ) = 0.25 * ( sum( z ( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 5  : &
         ( ele_hex2 - 1 ) * quad_cv_nloc + 8 ) ) ) )

!!!
!!! Level 2
!!!
    ! Node 19
    x( inod( 11 ) ) = 0.5 * ( x( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 1 ) ) + &
         x( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 5 ) ) )
    y( inod( 11 ) ) = 0.5 * ( y( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 1 ) ) + &
         y( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 5 ) ) )
    z( inod( 11 ) ) = 0.5 * ( z( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 1 ) ) + &
         z( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 5 ) ) )

    ! Node 20
    x( inod( 12 ) ) = 0.5 * ( x( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 2 ) ) + &
         x( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 6 ) ) )
    y( inod( 12 ) ) = 0.5 * ( y( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 2 ) ) + &
         y( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 6 ) ) )
    z( inod( 12 ) ) = 0.5 * ( z( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 2 ) ) + &
         z( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 6 ) ) )

    ! Node 21
    x( inod( 13 ) ) = 0.5 * ( x( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 3 ) ) + &
         x( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 7 ) ) )
    y( inod( 13 ) ) = 0.5 * ( y( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 3 ) ) + &
         y( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 7 ) ) )
    z( inod( 13 ) ) = 0.5 * ( z( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 3 ) ) + &
         z( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 7 ) ) )

    ! Node 22 
    x( inod( 14 ) ) = 0.5 * ( x( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 4 ) ) + &
         x( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 8 ) ) )
    y( inod( 14 ) ) = 0.5 * ( y( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 4 ) ) + &
         y( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 8 ) ) )
    z( inod( 14 ) ) = 0.5 * ( z( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 4 ) ) + &
         z( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 8 ) ) )

    ! Node 23
    jnod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 19 )
    knod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 20 )

    x( inod( 15 ) ) = 0.5 * ( x( jnod ) + x( knod ) )
    y( inod( 15 ) ) = 0.5 * ( y( jnod ) + y( knod ) )
    z( inod( 15 ) ) = 0.5 * ( z( jnod ) + z( knod ) )

    ! Node 24
    jnod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 20 )
    knod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 22 )

    x( inod( 16 ) ) = 0.5 * ( x( jnod ) + x( knod ) )
    y( inod( 16 ) ) = 0.5 * ( y( jnod ) + y( knod ) )
    z( inod( 16 ) ) = 0.5 * ( z( jnod ) + z( knod ) )

    ! Node 25
    jnod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 21 )
    knod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 22 )

    x( inod( 17 ) ) = 0.5 * ( x( jnod ) + x( knod ) )
    y( inod( 17 ) ) = 0.5 * ( y( jnod ) + y( knod ) )
    z( inod( 17 ) ) = 0.5 * ( z( jnod ) + z( knod ) )

    ! Node 26
    jnod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 19 )
    knod = x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 21 )

    x( inod( 18 ) ) = 0.5 * ( x( jnod ) + x( knod ) )
    y( inod( 18 ) ) = 0.5 * ( y( jnod ) + y( knod ) )
    z( inod( 18 ) ) = 0.5 * ( z( jnod ) + z( knod ) )

    ! Node 27
    x( inod( 19 ) ) = 0.25 * ( sum( x ( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 19 : &
         ( ele_hex2 - 1 ) * quad_cv_nloc + 22 ) ) ) )
    y( inod( 19 ) ) = 0.25 * ( sum( y ( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 19 : &
         ( ele_hex2 - 1 ) * quad_cv_nloc + 22 ) ) ) )
    z( inod( 19 ) ) = 0.25 * ( sum( z ( x_ndgln( ( ele_hex2 - 1 ) * quad_cv_nloc + 19 : &
         ( ele_hex2 - 1 ) * quad_cv_nloc + 22 ) ) ) )

    ewrite(3,*)'X / Y / Z of parametric nodes:'
    do iloc = 1, 19
       ewrite(3,*) inod( iloc ), x( inod ( iloc ) ), y( inod ( iloc ) ), z( inod ( iloc ) )
    end do

    ewrite(3,*)'Checking the initial numbering - X_NDGLN22, for the hexs:'
    ewrite(3,*) ' ele, ele_hex, xndgln()=, xndgln:'
    do ele_hex2 = 1, ele_hex 
       kloc = ( ele - 1 ) * number_of_hexs + ele_hex2
       do iloc = 1, 27
          ewrite(3,*) ele, ele_hex2, ( kloc - 1 ) * quad_cv_nloc + iloc, &
               x_ndgln( ( kloc - 1 ) * quad_cv_nloc + iloc ) 
       end do
    end do

    deallocate( inod ) 

    return
  end subroutine Adding_Parametric_Nodes_Hex


!!!!!!!!
!!!!!!!!
!!!!!!!!


   SUBROUTINE SHAPE_one_ele2(&
         ndim, cv_ele_type, &
         cv_ngi, cv_nloc, u_nloc,  &
                                ! Volume shape functions
         cvweight, cvfen, cvfenlx, cvfenly, cvfenlz, &
         ufen, ufenlx, ufenly, ufenlz, &
                                ! Surface of each CV shape functions
         sbcvngi,  &
         sbcvfen, sbcvfenslx, sbcvfensly, sbcvfeweigh, &
         sbufen, sbufenslx, sbufensly, &
                                ! Surface element shape funcs
         nface, &
         cv_sloclist, u_sloclist, cv_snloc, u_snloc ) 
      ! This subrt defines the sub-control volume and FEM shape functions.
      ! Shape functions associated with volume integration using both CV basis 
      ! functions CVN as well as FEM basis functions CVFEN (and its derivatives 
      ! CVFENLX, CVFENLY, CVFENLZ)
      implicit none
      integer, intent( in ) :: ndim, cv_ele_type, cv_ngi, cv_nloc, u_nloc
      real, dimension( cv_ngi ), intent( inout ) :: cvweight
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: cvfen, cvfenlx, cvfenly, cvfenlz
      real, dimension( u_nloc, cv_ngi ), intent( inout ) :: ufen, ufenlx, ufenly, ufenlz
      integer, intent( in ) :: sbcvngi
      integer, intent( in ) :: cv_snloc, u_snloc
      real, dimension( cv_snloc, sbcvngi ), intent( inout ) :: sbcvfen, sbcvfenslx, sbcvfensly
      real, dimension( sbcvngi ), intent( inout ) :: sbcvfeweigh
      real, dimension( u_snloc, sbcvngi ), intent( inout ) :: sbufen, sbufenslx, sbufensly
      integer, intent( in ) :: nface
      integer, dimension( nface, cv_snloc ), intent( inout ) :: cv_sloclist
      integer, dimension( nface, u_snloc ), intent( inout ) :: u_sloclist

! local variables...
      real, dimension( :, : ), allocatable :: M,MLX,MLY,MLZ, sm,SMLX,SMLY
      integer :: MLOC,SMLOC,NWICEL
      logical :: LOWQUA,d3
      REAL :: RUB(1000)

      ewrite(3,*)'just inside SHAPE_one_ele' 
!      stop 7299

      LOWQUA=.false.
      MLOC=1
      SMLOC=1
      ALLOCATE(M(MLOC,CV_NGI))
      ALLOCATE(MLX(MLOC,CV_NGI))
      ALLOCATE(MLY(MLOC,CV_NGI))
      ALLOCATE(MLZ(MLOC,CV_NGI))
      ALLOCATE(SM(SMLOC,sbcvngi))
      ALLOCATE(SMLX(SMLOC,sbcvngi))
      ALLOCATE(SMLY(SMLOC,sbcvngi))
      d3=(ndim==3)
 
! for pressure...
      Conditional_Dimensionality1: if( d3 ) then
         nwicel = 2
         if( cv_nloc == 8 ) then ! Linear hex
            nwicel = 1
         else if( cv_nloc == 27 ) then ! Quadratic hex
            nwicel = 3
         else if( cv_nloc == 4 ) then ! Linear tets
            nwicel = 4
         else if( cv_nloc == 10 ) then ! Quadratic tets
            nwicel = 5
         end if
      else
         nwicel = 2
         if( cv_nloc == 4 ) then ! Linear hex
            nwicel = 1
         else if( cv_nloc == 9 ) then ! Quadratic hex
            nwicel = 3
         else if ( cv_nloc == 3 ) then ! Linear tets
            nwicel = 4
         else if ( cv_nloc == 6 ) then ! Quadratic tets
            nwicel = 5
         end if
      end if Conditional_Dimensionality1
            
! for pressure...
           CALL SHAPE(LOWQUA,cv_NGI,cv_NLOC,MLOC, sbcvngi,cv_SNLOC,SMLOC,  &
        M,MLX,MLY,MLZ,cvWEIGHT,cvfen, cvfenlx, cvfenly, cvfenlz ,         & 
        sbcvfeweigh,sbcvfen, sbcvfenslx, sbcvfensly, SM,SMLX,SMLY,            &
        NWICEL,ndim==3)


! for velocity...
      Conditional_Dimensionality2: if( d3 ) then
         nwicel = 2
         if( u_nloc == 8 ) then ! Linear hex
            nwicel = 1
         else if( u_nloc == 27 ) then ! Quadratic hex
            nwicel = 3
         else if( u_nloc == 4 ) then ! Linear tets
            nwicel = 4
         else if( u_nloc == 10 ) then ! Quadratic tets
            nwicel = 5
         end if
      else
         nwicel = 2
         if( u_nloc == 4 ) then ! Linear hex
            nwicel = 1
         else if( u_nloc == 9 ) then ! Quadratic hex
            nwicel = 3
         else if ( u_nloc == 3 ) then ! Linear tets
            nwicel = 4
         else if ( u_nloc == 6 ) then ! Quadratic tets
            nwicel = 5
         end if
      end if Conditional_Dimensionality2
! for velocity...
           CALL SHAPE(LOWQUA,cv_NGI,u_NLOC,MLOC, sbcvngi,u_SNLOC,SMLOC,    &
        M,MLX,MLY,MLZ,RUB,ufen, ufenlx, ufenly, ufenlz,              & 
        RUB,sbufen, sbufenslx, sbufensly, SM,SMLX,SMLY,               &
        NWICEL,ndim==3)
      ! Determine CV_SLOCLIST & U_SLOCLIST
      CALL DETERMIN_SLOCLIST( CV_SLOCLIST, CV_NLOC, CV_SNLOC, NFACE, &
           NDIM, CV_ELE_TYPE )
      IF( U_SNLOC == 1 ) THEN
         U_SLOCLIST( 1, 1 ) = 1
         U_SLOCLIST( 2, 1 ) = U_NLOC
      ELSE
         CALL DETERMIN_SLOCLIST( U_SLOCLIST, U_NLOC, U_SNLOC, NFACE, &
              NDIM, CV_ELE_TYPE )
      ENDIF

      IF(NDIM.LT.3) THEN
         CVfenlz=0.0
         ufenlz=0.0
         sbcvfensly=0.0
         sbufensly=0.0
      ENDIF
      IF(NDIM.LT.2) THEN
         CVfenly=0.0
         ufenly=0.0
         sbcvfenslx=0.0
         sbufenslx=0.0
      ENDIF

      ewrite(3,*)'ndim, cv_ele_type,cv_ngi, cv_nloc, u_nloc:', &
               ndim, cv_ele_type,cv_ngi, cv_nloc, u_nloc
      ewrite(3,*)'cvweight:',cvweight
      ewrite(3,*)'cvfen:',cvfen
      ewrite(3,*)'cvfenlx:',cvfenlx
      ewrite(3,*)'cvfenly:',cvfenly
      ewrite(3,*)'cvfenlz:',cvfenlz 
      ewrite(3,*)'ufen:',ufen
      ewrite(3,*)'ufenlx:',ufenlx
      ewrite(3,*)'ufenly:',ufenly
      ewrite(3,*)'ufenlz:',ufenlz
      ewrite(3,*)'sbcvngi=',sbcvngi
      ewrite(3,*)'sbcvfen:',sbcvfen
      ewrite(3,*)'sbcvfenslx:',sbcvfenslx
      ewrite(3,*)'sbcvfensly:',sbcvfensly
      ewrite(3,*)'sbcvfeweigh:',sbcvfeweigh
                 
      ewrite(3,*)'sbufen:',sbufen
      ewrite(3,*)'sbufenslx:',sbufenslx
      ewrite(3,*)'sbufensly:',sbufensly

      ewrite(3,*)'nface:',nface
      ewrite(3,*)'cv_sloclist:', cv_sloclist
      ewrite(3,*)'u_sloclist:',u_sloclist
      ewrite(3,*)'cv_snloc, u_snloc:',cv_snloc, u_snloc
!      stop 145

      END SUBROUTINE SHAPE_one_ele2




   SUBROUTINE SHAPE(LOWQUA,NGI,NLOC,MLOC, SNGI,SNLOC,SMLOC,    &
        M,MLX,MLY,MLZ,WEIGHT,N,NLX,NLY,NLZ,                    & 
        SWEIGH,SN,SNLX,SNLY, SM,SMLX,SMLY,                     &
        NWICEL,D3)
     LOGICAL, INTENT(IN)::LOWQUA
     INTEGER, INTENT(IN)::NGI,NLOC,MLOC,SNGI,SNLOC,SMLOC
     REAL, INTENT(OUT)::M(MLOC,NGI),MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
     REAL, INTENT(OUT)::WEIGHT(NGI)
     REAL, INTENT(OUT)::N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
     REAL, INTENT(OUT)::SWEIGH(SNGI)
     REAL, INTENT(OUT)::SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
     REAL, INTENT(OUT)::SM(SMLOC,SNGI),SMLX(SMLOC,SNGI),SMLY(SMLOC,SNGI)
     INTEGER, INTENT(IN)::NWICEL
     LOGICAL, INTENT(IN)::D3
     
     INTEGER IPOLY,IQADRA,gi,gj,ggi,i,j,ii
     
     ewrite(3,*)'inside shape LOWQUA,NGI,NLOC,MLOC, SNGI,SNLOC,SMLOC:', &
                           LOWQUA,NGI,NLOC,MLOC, SNGI,SNLOC,SMLOC
     ewrite(3,*)'NWICEL,d3:',NWICEL,d3
     IF(NWICEL.EQ.1) THEN
        IF(.NOT.D3) THEN
           CALL RE2DN4(LOWQUA,NGI,0,NLOC,MLOC, &
                M,WEIGHT,N,NLX,NLY,          &
                SNGI,SNLOC,SWEIGH,SN,SNLX, &
                m,m)
        ELSE
           CALL RE3DN8(LOWQUA,NGI,0,NLOC,MLOC, &
                M,WEIGHT,N,NLX,NLY,NLZ,      &
                SNGI,SNLOC,SWEIGH,SN,SNLX,SNLY, &
                m,m,m)
        ENDIF
     ENDIF
     
     IF(NWICEL.EQ.2) THEN
        ewrite(3,*)'option not avaialble'
        stop 3832
!        IF(.NOT.D3) THEN
!           CALL RE2DN8(LOWQUA,NGI,NLOC,MLOC, &
!                M,WEIGHT,N,NLX,NLY )
!        ELSE
!           ! SERENDIPITY 20 NODE 3-D ELEMENT -BILINEAR PRESSURE
!           CALL RE3D20(LOWQUA,NGI,NLOC,MLOC, &
!                M,WEIGHT,N,NLX,NLY,NLZ )
!        ENDIF
     ENDIF
     
     IF(NWICEL.EQ.3) THEN
        IF(.NOT.D3) THEN
           CALL RE2DN9(LOWQUA,NGI,0,NLOC,MLOC, &
                M,WEIGHT,N,NLX,NLY, &
                m,m)
           sweigh=0.0
           sn=0.0
           snlx=0.0
           do gi=1,3
             do gj=1,3
               ggi=(gj-1)*3+gi
               sweigh(gi)=sweigh(gi)+WEIGHT(ggi)
               do i=1,3
                 do j=1,3
                   ii=(j-1)*3+i
                   sn(i,gi)=sn(i,gi)+n(ii,ggi)
                   snlx(i,gi)=snlx(i,gi)+nlx(ii,ggi)
                 end do
               end do
             end do
           end do
        ELSE
           ! LAGRANGE 27 NODE 3-D ELEMENT -BILINEAR PRESSURE
           CALL RE3D27(LOWQUA,NGI,0,NLOC,MLOC, &
                M,WEIGHT,N,NLX,NLY,NLZ, &
                m,m,m)
           CALL RE2DN9(LOWQUA,SNGI,0,SNLOC,MLOC, &
                M,SWEIGH,SN,SNLX,SNLY, &
                m,m)
        ENDIF
     ENDIF
     
     IF((NWICEL.EQ.4).or.(NWICEL.EQ.5)) THEN 
! works for linear or quadratic triangles or tets...
           CALL TR2or3DQU(NGI,NLOC,MLOC,        &
                M,MLX,MLY,MLZ,                  &
                WEIGHT,N,NLX,NLY,NLZ,           &
                SNGI,SNLOC,SWEIGH,SN,SNLX,SNLY, &
                SMLOC,                          &
                SM,SMLX,SMLY,D3)
!       ewrite(3,*)'weight:',weight
!      STOP 3321
     ENDIF
     
     IF(NWICEL.GE.100) THEN
        ! A Spectal element using Legendra, Lagrange or Chebichef polynomials. 
        CALL SPECTR(NGI,NLOC,MLOC, &
             M,WEIGHT,N,NLX,NLY,NLZ,D3,.NOT.D3, IPOLY,IQADRA)
     ENDIF
     
   END SUBROUTINE SHAPE
   



   SUBROUTINE TR2DQU(NGI,NLOC,MLOC,&
     &      M,MLX,MLY,&
     &      WEIGHT,N,NLX,NLY, &
     &      SNGI,SNLOC,SWEIGH,SN,SNLX,&
     &      SMLOC,&
     &      SM,SMLX )
!      This subroutine defines the shape functions M and N and their
!      derivatives at the Gauss points for quadratic elements. 
! For 3-D FLOW. 
      INTEGER NLOC,MLOC,NGI,SNGI,SNLOC,SMLOC
      REAL M(MLOC,NGI),MLX(MLOC,NGI),MLY(MLOC,NGI)
      REAL WEIGHT(NGI)
      REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI)
      REAL SWEIGH(SNGI)
      REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI)
      REAL SM(SMLOC,SNGI),SMLX(SMLOC,SNGI)
! Local variables...
      REAL POSI,TLY
      INTEGER IPOLY,IQADRA
      REAL L1(20),L2(20),L3(20)
      REAL RUB(10)
      LOGICAL DD3
! NB LXP(I) AND LYP(I) ARE THE LOCAL X AND Y COORDS OF NODAL POINT I
       
      ewrite(3,*) 'HERE 1 MLOC,NLOC,NGI=',MLOC,NLOC,NGI
      ewrite(3,*) 'HERE 2'
!
! Get the quadrature positions and weights for TRIANGLES...
     DD3=.FALSE.
     CALL TRIQUAold(L1, L2, L3, L3, WEIGHT, DD3,NGI)

! Work out the shape functions and there derivatives...
        CALL SHATRIold(L1, L2, L3, L3, WEIGHT, DD3,&
     &              NLOC,NGI,&
     &              N,NLX,NLY,NLY) 
        CALL SHATRIold(L1, L2, L3, L3, WEIGHT, DD3,&
     &              MLOC,NGI,&
     &              M,MLX,MLY,MLY) 


      IF(SNGI.GT.0) THEN
        ewrite(3,*)'for surfaces SNGI,SNLOC,smloc:',SNGI,SNLOC,smloc
! IQADRA=1 corresponds to Gaussian quadrature.
         IQADRA=1
! IPOLY=1 is for Lagrange polynomials.
         IPOLY=1

          ewrite(3,*)'for sn:'
         CALL SPECTR(SNGI,SNLOC,0,&
     &   RUB,SWEIGH,SN,SNLX,SNLX,SNLX,.FALSE.,.FALSE., IPOLY,IQADRA)

          ewrite(3,*)'for sm:'
         CALL SPECTR(SNGI,SMLOC,0,&
     &   RUB,SWEIGH,SM,SMLX,SMLX,SMLX,.FALSE.,.FALSE., IPOLY,IQADRA)
      ENDIF
      
   END subroutine tr2dqu
     


   SUBROUTINE TR2or3DQU(NGI,NLOC,MLOC,&
        &     M,MLX,MLY,MLZ,&
        &     WEIGHT,N,NLX,NLY,NLZ,&
        &     SNGI,SNLOC,SWEIGH,SN,SNLX,SNLY,&
        &     SMLOC,&
        &     SM,SMLX,SMLY,D3)
     !     This subroutine defines the shape functions M and N and their
     !     derivatives at the Gauss points for quadratic elements. 
     !     For 3-D FLOW. 
     INTEGER SNGI,SNLOC,SMLOC,NLOC,MLOC,NGI
     REAL SWEIGH(SNGI)
     REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
     REAL SM(SMLOC,SNGI),SMLX(SMLOC,SNGI),SMLY(SMLOC,SNGI)
     REAL M(MLOC,NGI),MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
     REAL WEIGHT(NGI)
     REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
     LOGICAL D3
     ! Local variables...
     REAL RUB(500)
     LOGICAL DD3,base_order
     REAL L1(50),L2(50),L3(50),L4(50)
     integer IQADRA,IPOLY
     ! NB LXP(I) AND LYP(I) ARE THE LOCAL X AND Y COORDS OF NODAL POINT I

     ewrite(3,*) 'HERE 1 MLOC,NLOC,NGI=',MLOC,NLOC,NGI
     ewrite(3,*) 'HERE 2'

     ! Get the quadrature positions and weights for TRIANGLES or TETS...
     DD3=D3
     CALL TRIQUAold(L1, L2, L3, L4, WEIGHT, DD3,NGI)
     !ewrite(3,*)'l1:',l1(1:ngi)
     !ewrite(3,*)'l2:',l2(1:ngi)
     !ewrite(3,*)'l3:',l3(1:ngi)
     !if(d3) ewrite(3,*)'l4:',l4(1:ngi)
     !ewrite(3,*)'weight:',weight
     !stop 2821

     ! Work out the shape functions and there derivatives...
     CALL SHATRIold(L1, L2, L3, L4, WEIGHT, DD3,&
          &              NLOC,NGI,&
          &              N,NLX,NLY,NLZ)
     ! re-arrange ordering for quadratic elements...
     if(d3) then
        if(nloc==10) then
           base_order=.true.
           if(base_order) then
              ! order so that the 1st nodes are on the base...
              call base_order_tet(n,nloc,ngi)
              call base_order_tet(nlx,nloc,ngi)
              call base_order_tet(nly,nloc,ngi)
              call base_order_tet(nlz,nloc,ngi)
           endif
        endif
     else
        if(nloc==6) then
           base_order=.true.
           if(base_order) then
              ! order so that the 1st nodes are on the base...
              call base_order_tri(n,nloc,ngi)
              call base_order_tri(nlx,nloc,ngi)
              call base_order_tri(nly,nloc,ngi)
           endif
        endif
     endif
     ewrite(3,*)'n::',n
     ewrite(3,*)'nlx::',nlx
     ewrite(3,*)'nly::',nly
     CALL SHATRIold(L1, L2, L3, L4, WEIGHT, DD3,&
          MLOC,NGI,&
          M,MLX,MLY,MLZ) 

     IF(SNGI.GT.0) THEN

        IF(D3) THEN
           DD3=.FALSE.
           CALL TRIQUAold(L1, L2, L3, L4, SWEIGH, DD3,SNGI)

           ! Work out the shape functions and there derivatives...
           CALL SHATRIold(L1, L2, L3, L4, SWEIGH, DD3,&
                              SNLOC,SNGI,&
                              SN,SNLX,SNLY,RUB) 
           if(snloc==6) then
              base_order=.true.
              if(base_order) then
                 ! order so that the 1st nodes are on the base...
                 call base_order_tri(sn,snloc,sngi)
                 call base_order_tri(snlx,snloc,sngi)
                 call base_order_tri(snly,snloc,sngi)
              endif
           endif
           CALL SHATRIold(L1, L2, L3, L4, SWEIGH, DD3,&
                              SMLOC,NGI,&
                              SM,SMLX,SMLY,RUB) 
        ELSE
           ewrite(3,*)'for surfaces SNGI,SNLOC,smloc:',SNGI,SNLOC,smloc
           ! IQADRA=1 corresponds to Gaussian quadrature.
           IQADRA=1
           ! IPOLY=1 is for Lagrange polynomials.
           IPOLY=1

           ewrite(3,*)'for sn IPOLY,IQADRA,SNGI,SNLOC:', &
                IPOLY,IQADRA,SNGI,SNLOC
           CALL SPECTR(SNGI,SNLOC,0,&
                   RUB,SWEIGH,SN,SNLX,RUB,RUB,.FALSE.,.FALSE., IPOLY,IQADRA)
           ewrite(3,*)'+++for sn SWEIGH:',SWEIGH

           if(.false.) then
              ewrite(3,*)'for sm:'
              CALL SPECTR(SNGI,SMLOC,0,&
                      RUB,SWEIGH,SM,SMLX,RUB,RUB,.FALSE.,.FALSE., IPOLY,IQADRA)
           endif
        ENDIF

     ENDIF
     ! the weights need to sum to 0.5 in 2D triangles and 1./6. in 3D
     IF(D3) THEN
        !        WEIGHT=(1./6.)*WEIGHT
        WEIGHT=1.*WEIGHT
        SWEIGH=0.5*SWEIGH
     ELSE ! 2d...
        WEIGHT=0.5*WEIGHT
     ENDIF

   end subroutine tr2or3dqu


      
   SUBROUTINE TR2D(LOWQUA,NGI,NLOC,MLOC,&
         M,WEIGHT,N,NLX,NLY, &
          SNGI,SNLOC,SWEIGH,SN,SNLX )
!     This subroutine defines the shape functions M and N and their
!     derivatives at the Gauss points
!     For 3-D FLOW. 
      INTEGER NGI,NLOC,MLOC
      REAL M(MLOC,NGI),WEIGHT(NGI)
      REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI)
      INTEGER SNGI,SNLOC
      REAL SWEIGH(SNGI)
      REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI)
      REAL POSI,TLY
      INTEGER P,Q,CORN,GPOI,ILOC,JLOC,GI
      LOGICAL LOWQUA,GETNDP
      REAL WEIT(20),LX(20),LXP(2)
      INTEGER I
! NB LXP(I) AND LYP(I) ARE THE LOCAL X AND Y COORDS OF NODAL POINT I
       
      ewrite(3,*) 'HERE 1 MLOC,NLOC,NGI=',MLOC,NLOC,NGI
      ewrite(3,*) 'HERE 2'

      IF((NLOC.NE.3).OR.(NGI.NE.3)) THEN
          ewrite(3,*)'PROBLEM IN TR2D NLOC,NGI:',NLOC,NGI
          stop 282
      ENDIF

! This is for one point quadrature. 
      do  I=1,3
         NLX(I,1)=1.
         NLY(I,1)=1.
      end do

      WEIGHT(1)=1.

      N(1,1)=1./3.
      N(2,1)=1./3.
      N(3,1)=1./3.

! This is for 3 point quadrature. 
      do  GI=1,3
         NLX(1,GI)=1.
         NLY(1,GI)=0.

         NLX(2,GI)=0.
         NLY(2,GI)=1.

         NLX(3,GI)=-1.
         NLY(3,GI)=-1.
      end do


      N(1,1)=0.5
      N(2,1)=0.5
      N(3,1)=0.

      N(1,2)=0.
      N(2,2)=0.5
      N(3,2)=0.5

      N(1,3)=0.5
      N(2,3)=0.
      N(3,3)=0.5

      WEIGHT(1)=1./3.
      WEIGHT(2)=1./3.
      WEIGHT(3)=1./3.

      IF(MLOC.EQ.1) THEN
         M(1,1)=1.0
         M(1,2)=1.0
         M(1,3)=1.0
      ENDIF
      IF(MLOC.EQ.NLOC) M = N

      IF(SNGI.GT.0) THEN
         LXP(1)=-1
         LXP(2)= 1
         GETNDP=.FALSE.
         CALL LAGROT(WEIT,LX,SNGI,GETNDP)
         do  P=1,SNGI
            do  CORN=1,2! Was loop 327
               GPOI=P
               SN(CORN,GPOI)=0.5*(1.+LXP(CORN)*LX(P))
               SNLX(CORN,GPOI)=0.5*LXP(CORN)
               SWEIGH(GPOI)=WEIT(P)
            end do
         end do 
      ENDIF
      
   END subroutine tr2d




   SUBROUTINE TR3D(LOWQUA,NGI,NLOC,MLOC,&
          M,WEIGHT,N,NLX,NLY,NLZ,&
          SNGI,SNLOC,SWEIGH,SN,SNLX,SNLY)
!     This subroutine defines the shape functions M and N and their
!     derivatives at the Gauss points
!     For 3-D FLOW.
      INTEGER NGI,NLOC,MLOC
      REAL ALPHA,BETA
      PARAMETER(ALPHA=0.58541020,BETA=0.13819660)
      INTEGER SNGI,SNLOC
      REAL SWEIGH(SNGI)
      REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
      REAL M(MLOC,NGI),WEIGHT(NGI)
      REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
      REAL RUB(20)
      INTEGER P,Q,CORN,GPOI,ILOC,JLOC,GI
      LOGICAL LOWQUA
      INTEGER I
! NB LXP(I) AND LYP(I) ARE THE LOCAL X AND Y COORDS OF NODAL POINT I
       
      ewrite(3,*) 'HERE 1 MLOC,NLOC,NGI=',MLOC,NLOC,NGI
      ewrite(3,*) 'HERE 2'

      IF((NLOC.NE.4).OR.(NGI.NE.4)) THEN
         ewrite(3,*) 'PROBLEM IN TR3D'
         STOP 201
      ENDIF

! This is for one point. 

! This is for 4 point quadrature. 
      do  GI=1,NGI
         NLX(1,GI)=1.
         NLY(1,GI)=0.
         NLZ(1,GI)=0.

         NLX(2,GI)=0.
         NLY(2,GI)=1.
         NLZ(2,GI)=0.

         NLX(3,GI)=0.
         NLY(3,GI)=0.
         NLZ(3,GI)=1.

         NLX(4,GI)=-1.
         NLY(4,GI)=-1.
         NLZ(4,GI)=-1.
      end do 

      do I=1,4
         do GI=1,4
            N(I,GI)=BETA
         END DO
      END DO

      do I=1,4
         N(I,I)=ALPHA
         WEIGHT(I)=0.25
         IF(MLOC.EQ.1) M(1,I)=1.0
      END DO
      IF(MLOC.EQ.NLOC) M = N

      IF(SNGI.GT.0) THEN
         CALL TR2D(.FALSE.,SNGI,SNLOC,SNLOC,&
           RUB,SWEIGH,SN,SNLX,SNLY, &
           0,0,RUB,RUB,RUB )
      ENDIF 

      do I=1,NGI
         WEIGHT(I)=WEIGHT(I)/6.
      END DO
      do I=1,SNGI
         SWEIGH(I)=SWEIGH(I)/2.
      END DO
      
   end subroutine tr3d
      

!
!
   SUBROUTINE SHATRIold(L1, L2, L3, L4, WEIGHT, D3, &
        NLOC,NGI,  &
        N,NLX,NLY,NLZ) 
     ! Work out the shape functions and there derivatives...
     IMPLICIT NONE
     INTEGER NLOC,NGI
     LOGICAL D3
     REAL L1(NGI), L2(NGI), L3(NGI), L4(NGI)
     REAL WEIGHT(NGI)
     REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
     ! Local variables...
     INTEGER GI
     !
     IF(.NOT.D3) THEN
        ! Assume a triangle...
        IF((NLOC.EQ.6).OR.(NLOC.EQ.7)) THEN
           DO 10 GI=1,NGI
              N(1,GI)=(2.*L1(GI)-1.)*L1(GI)
              N(2,GI)=(2.*L2(GI)-1.)*L2(GI)
              N(3,GI)=(2.*L3(GI)-1.)*L3(GI)
              !
              N(4,GI)=4.*L1(GI)*L2(GI)
              N(5,GI)=4.*L2(GI)*L3(GI)
              N(6,GI)=4.*L1(GI)*L3(GI)

              !
              ! nb L1+L2+L3+L4=1
              ! x-derivative...
              NLX(1,GI)=4.*L1(GI)-1.
              NLX(2,GI)=0.
              NLX(3,GI)=-4.*(1.-L2(GI))+4.*L1(GI) + 1.
              !
              NLX(4,GI)=4.*L2(GI)
              NLX(5,GI)=-4.*L2(GI)
              NLX(6,GI)=4.*(1.-L2(GI))-8.*L1(GI)
              !
              ! y-derivative...
              NLY(1,GI)=0.
              NLY(2,GI)=4.*L2(GI)-1.0
              NLY(3,GI)=-4.*(1.-L1(GI))+4.*L2(GI) + 1.
              !
              NLY(4,GI)=4.*L1(GI)
              NLY(5,GI)=4.*(1.-L1(GI))-8.*L2(GI)
              NLY(6,GI)=-4.*L1(GI)
              IF(NLOC.EQ.7) THEN
                 ! Bubble function...
                 N(7,GI)  =L1(GI)*L2(GI)*L3(GI)
                 NLX(7,GI)=L2(GI)*(1.-L2(GI))-2.*L1(GI)*L2(GI)
                 NLY(7,GI)=L1(GI)*(1.-L1(GI))-2.*L1(GI)*L2(GI)
              ENDIF
10            CONTINUE
              ! ENDOF IF(NLOC.EQ.6) THEN...
           ENDIF
           !
           IF((NLOC.EQ.3).OR.(NLOC.EQ.4)) THEN
              DO 20 GI=1,NGI
                 N(1,GI)=L1(GI)
                 N(2,GI)=L2(GI)
                 N(3,GI)=L3(GI)
                 !
                 NLX(1,GI)=1.0
                 NLX(2,GI)=0.0
                 NLX(3,GI)=-1.0
                 !
                 NLY(1,GI)=0.0
                 NLY(2,GI)=1.0
                 NLY(3,GI)=-1.0
                 IF(NLOC.EQ.4) THEN
                    ! Bubble function...
                    N(4,GI)  =L1(GI)*L2(GI)*L3(GI)
                    NLX(4,GI)=L2(GI)*(1.-L2(GI))-2.*L1(GI)*L2(GI)
                    NLY(4,GI)=L1(GI)*(1.-L1(GI))-2.*L1(GI)*L2(GI)
                 ENDIF
20               CONTINUE
              ENDIF
              !
              IF(NLOC.EQ.1) THEN
                 DO 30 GI=1,NGI
                    N(1,GI)=1.0
                    NLX(1,GI)=0.0
                    NLY(1,GI)=0.0
30                  CONTINUE
                 ENDIF
                 !
                 ! ENDOF IF(.NOT.D3) THEN
              ENDIF
              !
              !
              IF(D3) THEN
                 ! Assume a tet...
                 ! This is for 5 point quadrature. 
                 IF((NLOC.EQ.10).OR.(NLOC.EQ.11)) THEN
                    DO 40 GI=1,NGI
                       !         ewrite(3,*)'gi,L1(GI),L2(GI),L3(GI),L4(GI):',gi,L1(GI),L2(GI),L3(GI),L4(GI)
                       N(1,GI)=(2.*L1(GI)-1.)*L1(GI)
                       N(3,GI)=(2.*L2(GI)-1.)*L2(GI)
                       N(5,GI)=(2.*L3(GI)-1.)*L3(GI)
                       N(10,GI)=(2.*L4(GI)-1.)*L4(GI)

                       !if(L1(GI).gt.-1.93) ewrite(3,*)'gi,L1(GI), L2(GI), L3(GI), L4(GI),N(1,GI):', &
                       !                            gi,L1(GI), L2(GI), L3(GI), L4(GI),N(1,GI)
                       !
                       !
                       N(2,GI)=4.*L1(GI)*L2(GI)
                       N(6,GI)=4.*L1(GI)*L3(GI)
                       N(7,GI)=4.*L1(GI)*L4(GI)
                       !
                       N(4,GI) =4.*L2(GI)*L3(GI)
                       N(9,GI) =4.*L3(GI)*L4(GI)
                       N(8,GI)=4.*L2(GI)*L4(GI)
                       ! nb L1+L2+L3+L4=1
                       ! x-derivative...
                       NLX(1,GI)=4.*L1(GI)-1.
                       NLX(3,GI)=0.
                       NLX(5,GI)=0.
                       NLX(10,GI)=-4.*(1.-L2(GI)-L3(GI))+4.*L1(GI) + 1.
                       !if(L1(GI).gt.-1.93) ewrite(3,*)'Nlx(1,GI):', &
                       !     Nlx(1,GI)
                       !
                       NLX(2,GI)=4.*L2(GI)
                       NLX(6,GI)=4.*L3(GI)
                       NLX(7,GI)=4.*(L4(GI)-L1(GI))
                       !
                       NLX(4,GI) =0.
                       NLX(9,GI) =-4.*L3(GI)
                       NLX(8,GI)=-4.*L2(GI)
                       !
                       ! y-derivative...
                       NLY(1,GI)=0.
                       NLY(3,GI)=4.*L2(GI)-1.0
                       NLY(5,GI)=0.
                       NLY(10,GI)=-4.*(1.-L1(GI)-L3(GI))+4.*L2(GI) + 1.
                       !
                       NLY(2,GI)=4.*L1(GI)
                       NLY(6,GI)=0.
                       NLY(7,GI)=-4.*L1(GI)
                       !
                       NLY(4,GI) =4.*L3(GI)
                       NLY(9,GI) =-4.*L3(GI)
                       NLY(8,GI)=4.*(1-L1(GI)-L3(GI))-8.*L2(GI)
                       !
                       ! z-derivative...
                       NLZ(1,GI)=0.
                       NLZ(3,GI)=0.
                       NLZ(5,GI)=4.*L3(GI)-1.
                       NLZ(10,GI)=-4.*(1.-L1(GI)-L2(GI))+4.*L3(GI) + 1.
                       !
                       NLZ(2,GI)=0.
                       NLZ(6,GI)=4.*L1(GI)
                       NLZ(7,GI)=-4.*L1(GI)
                       !
                       NLZ(4,GI) =4.*L2(GI)
                       NLZ(9,GI) =4.*(1.-L1(GI)-L2(GI))-8.*L3(GI)
                       NLZ(8,GI)=-4.*L2(GI)
                       IF(NLOC.EQ.11) THEN
                          ! Bubble function...
                          N(11,GI)  =L1(GI)*L2(GI)*L3(GI)*L4(GI)
                          NLX(11,GI)=L2(GI)*L3(GI)*(1.-L2(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
                          NLY(11,GI)=L1(GI)*L3(GI)*(1.-L1(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
                          NLZ(11,GI)=L1(GI)*L2(GI)*(1.-L1(GI)-L2(GI))-2.*L1(GI)*L2(GI)*L3(GI)
                       ENDIF
                       !
40                     CONTINUE
                       ! ENDOF IF(NLOC.EQ.10) THEN...
                    ENDIF 
                    !
                    IF((NLOC.EQ.4).OR.(NLOC.EQ.5)) THEN
                       DO 50 GI=1,NGI
                          N(1,GI)=L1(GI)
                          N(2,GI)=L2(GI)
                          N(3,GI)=L3(GI)
                          N(4,GI)=L4(GI)
                          !
                          NLX(1,GI)=1.0
                          NLX(2,GI)=0
                          NLX(3,GI)=0
                          NLX(4,GI)=-1.0
                          !
                          NLY(1,GI)=0.0
                          NLY(2,GI)=1.0
                          NLY(3,GI)=0.0
                          NLY(4,GI)=-1.0
                          !
                          NLZ(1,GI)=0.0
                          NLZ(2,GI)=0.0
                          NLZ(3,GI)=1.0
                          NLZ(4,GI)=-1.0
                          IF(NLOC.EQ.5) THEN
                             ! Bubble function...
                             N(5,GI)  =L1(GI)*L2(GI)*L3(GI)*L4(GI)
                             NLX(5,GI)=L2(GI)*L3(GI)*(1.-L2(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
                             NLY(5,GI)=L1(GI)*L3(GI)*(1.-L1(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
                             NLZ(5,GI)=L1(GI)*L2(GI)*(1.-L1(GI)-L2(GI))-2.*L1(GI)*L2(GI)*L3(GI)
                          ENDIF
50                        CONTINUE 
                       ENDIF
                       !
                       IF(NLOC.EQ.1) THEN
                          DO 60 GI=1,NGI
                             N(1,GI)=1.0
                             NLX(1,GI)=0.0
                             NLY(1,GI)=0.0
                             NLZ(1,GI)=0.0
60                           CONTINUE 
                          ENDIF
                          !
                          ! ENDOF IF(D3) THEN...
                       ENDIF
                       !
                       RETURN
                     END SUBROUTINE SHATRIold
!
!
!
!
        SUBROUTINE TRIQUAold(L1, L2, L3, L4, WEIGHT, D3,NGI)
! This sub calculates the local corrds L1, L2, L3, L4 and 
! weights at the quadrature points. 
! If D3 it does this for 3Dtetrahedra elements else 
! triangular elements.
        IMPLICIT NONE
        INTEGER NGI
        LOGICAL D3
        REAL L1(NGI), L2(NGI), L3(NGI), L4(NGI), WEIGHT(NGI)
! Local variables...
        REAL ALPHA,BETA
        REAL ALPHA1,BETA1
        REAL ALPHA2,BETA2
        INTEGER I
!
        IF(D3) THEN
! this is for a tetrahedra element...
! This is for one point. 
          IF(NGI.EQ.1) THEN
! Degree of precision is 1
             DO I=1,NGI
               L1(I)=0.25
               L2(I)=0.25
               L3(I)=0.25
               L4(I)=0.25
               WEIGHT(I)=1.0
             END DO
           ENDIF
       
           IF(NGI.EQ.4) THEN
! Degree of precision is 2
              ALPHA=0.58541020
              BETA=0.13819660
              DO I=1,NGI
                 L1(I)=BETA
                 L2(I)=BETA
                 L3(I)=BETA
                 L4(I)=BETA
                 WEIGHT(I)=0.25
              END DO
              L1(1)=ALPHA
              L2(2)=ALPHA
              L3(3)=ALPHA
              L4(4)=ALPHA
           ENDIF
       
           IF(NGI.EQ.5) THEN
! Degree of precision is 3
             L1(1)=0.25
             L2(1)=0.25
             L3(1)=0.25
             L4(1)=0.25
             WEIGHT(1)=-4./5.
!
              DO I=2,NGI
                 L1(I)=1./6.
                 L2(I)=1./6.
                 L3(I)=1./6.
                 L4(I)=1./6.
                 WEIGHT(I)=9./20.
              END DO
             L1(2)=0.5
             L2(3)=0.5
             L3(4)=0.5
             L4(5)=0.5
       ENDIF
!
       IF(NGI.EQ.11) THEN
! Degree of precision is 4
         ALPHA=(1.+SQRT(5./14.))/4.0
         BETA =(1.-SQRT(5./14.))/4.0
         I=1
               L1(I)=0.25
               L2(I)=0.25
               L3(I)=0.25
               WEIGHT(I)=-6.*74.0/5625.0
       DO I=2,5
               L1(I)=1./14.
               L2(I)=1./14.
               L3(I)=1./14.
               WEIGHT(I)=6.*343./45000.
       END DO
               L1(2)=11./14.
               L2(3)=11./14.
               L3(4)=11./14.
       DO I=6,11
               L1(I)=ALPHA
               L2(I)=ALPHA
               L3(I)=ALPHA
               WEIGHT(I)=6.*56.0/2250.0
       END DO
               L3(6)=BETA
               L2(7)=BETA
               L2(8)=BETA
               L3(8)=BETA
               L1(9)=BETA
               L1(10)=BETA
               L3(10)=BETA
               L1(11)=BETA
               L2(11)=BETA
! ENDOF IF(NGI.EQ.11) THEN...
       ENDIF
         DO I=1,NGI
                 L4(I)=1.0-L1(I)-L2(I)-L3(I)
         END DO
! Now multiply by 1/6. to get weigts correct...
         DO I=1,NGI
           WEIGHT(I)=WEIGHT(I)/6.
         END DO
! ENDOF IF(D3) THEN...
        ENDIF
!
        IF(.NOT.D3) THEN
! 2-D TRAINGULAR ELEMENTS...
          IF(NGI.EQ.1) THEN
! LINEAR
            I=1
                  L1(I)=1./3.
                  L2(I)=1./3.
                  WEIGHT(I)=1.0
          ENDIF
!
          IF(NGI.EQ.3) THEN
! QUADRASTIC
            DO I=1,NGI
                    L1(I)=0.5
                    L2(I)=0.5
                    WEIGHT(I)=1.0/3.0
            END DO
                    L1(2)=0.0
                    L2(3)=0.0
           ENDIF
!
          IF(NGI.EQ.4) THEN
! CUBIC
            I=1
                  L1(I)=1./3.
                  L2(I)=1./3.
                  WEIGHT(I)=-27./48.
            DO I=2,NGI
                    L1(I)=0.2
                    L2(I)=0.2
                    WEIGHT(I)=25./48.
            END DO
                    L1(1)=0.6
                    L2(2)=0.6
           ENDIF
!
          IF(NGI.EQ.7) THEN
! QUNTIC
         ALPHA1=0.0597158717
         BETA1 =0.4701420641
         ALPHA2=0.7974269853
         BETA2 =0.1012865073
            I=1
                    L1(I)=1./3.
                    L2(I)=1./3.
                    WEIGHT(I)=0.225
            DO I=2,4
                    L1(I)=BETA1
                    L2(I)=BETA1
                    WEIGHT(I)=0.1323941527
            END DO
                  L1(2)=ALPHA1
                  L2(4)=ALPHA1
            DO I=5,7
                    L1(I)=BETA2
                    L2(I)=BETA2
                    WEIGHT(I)=0.1259391805
            END DO
                  L1(5)=ALPHA2
                  L2(6)=ALPHA2
! ENDOF IF(NGI.EQ.7) THEN...
          ENDIF
!
          DO I=1,NGI
                  L3(I)=1.0-L1(I)-L2(I)
          END DO
! ENDOF IF(.NOT.D3) THEN...
        ENDIF
!
        RETURN 
        END subroutine TRIQUAold
!
!
!
!
       SUBROUTINE SPECTR(NGI,NLOC,MLOC,&
     &      M,WEIGHT,N,NLX,NLY,NLZ,D3,D2, IPOLY,IQADRA  )
       IMPLICIT NONE
       INTEGER NGI,NLOC,MLOC,IPOLY,IQADRA
       REAL M(MLOC,NGI),WEIGHT(NGI)
       REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI)
        REAL NLZ(NLOC,NGI)
         REAL RGPTWE
         REAL WEIT(30),NODPOS(30),QUAPOS(30)
         INTEGER GPOI
         LOGICAL DIFF,NDIFF,D3,D2
       INTEGER NDGI,NDNOD,NMDNOD,IGR,IGQ,IGP,KNOD,JNOD,INOD,ILOC
       REAL LXGP,LYGP,LZGP
! This subroutine defines a spectal element.
! IPOLY defines the element type and IQADRA the quadrature. 
! In 2-D the spectral local node numbering is as..
! 7 8 9
! 4 5 6
! 1 2 3 
! For 3-D...
! lz=-1
! 3 4
! 1 2
! and for lz=1
! 7 8
! 5 6
       
         ewrite(3,*)'inside SPECTR IPOLY,IQADRA', IPOLY,IQADRA
!
         DIFF=.TRUE.
         NDIFF=.FALSE.
         IF(D3) THEN
         NDGI  =INT((NGI**(1./3.))+0.1)
         NDNOD =INT((NLOC**(1./3.))+0.1)
         NMDNOD=INT((MLOC**(1./3.))+0.1)
!
! Find the roots of the quadrature points and nodes
! also get the weights. 
       ewrite(3,*)'about to go into inside GTROOT IPOLY,IQADRA',IPOLY,IQADRA
        CALL GTROOT(IPOLY,IQADRA,WEIT,NODPOS,QUAPOS,NDGI,NDNOD)
       ewrite(3,*)'outside GTROOT'
      do  IGR=1,NDGI! Was loop 101
      do  IGQ=1,NDGI! Was loop 101
      do  IGP=1,NDGI! Was loop 101
         GPOI=IGP + (IGQ-1)*NDGI + (IGR-1)*NDGI*NDGI
!
!           WEIGHT(GPOI)
!     &        =RGPTWE(IGP,NDGI,.TRUE.)*RGPTWE(IGQ,NDGI,.TRUE.)
!     &        *RGPTWE(IGR,NDGI,.TRUE.)
           WEIGHT(GPOI)=WEIT(IGP)*WEIT(IGQ)*WEIT(IGR)
!
           LXGP=QUAPOS(IGP)
           LYGP=QUAPOS(IGQ)
           LZGP=QUAPOS(IGR)
! NB If TRUE in function RGPTWE then return the Gauss-pt weight 
! else return the Gauss-pt. 
!
      do  KNOD=1,NDNOD! Was loop 20
      do  JNOD=1,NDNOD! Was loop 20
      do  INOD=1,NDNOD! Was loop 20
           ILOC=INOD + (JNOD-1)*NDNOD + (KNOD-1)*NDNOD*NDNOD
!
             N(ILOC,GPOI)&
     &        =SPECFU(NDIFF,LXGP,INOD,NDNOD,IPOLY,NODPOS)&
     &        *SPECFU(NDIFF,LYGP,JNOD,NDNOD,IPOLY,NODPOS)&
     &        *SPECFU(NDIFF,LZGP,KNOD,NDNOD,IPOLY,NODPOS)
!
             NLX(ILOC,GPOI)&
     &        =SPECFU(DIFF,LXGP, INOD,NDNOD,IPOLY,NODPOS)&
     &        *SPECFU(NDIFF,LYGP,JNOD,NDNOD,IPOLY,NODPOS)&
     &        *SPECFU(NDIFF,LZGP,KNOD,NDNOD,IPOLY,NODPOS)
!
             NLY(ILOC,GPOI)&
     &        =SPECFU(NDIFF,LXGP,INOD,NDNOD,IPOLY,NODPOS)&
     &        *SPECFU(DIFF, LYGP,JNOD,NDNOD,IPOLY,NODPOS)&
     &        *SPECFU(NDIFF,LZGP,KNOD,NDNOD,IPOLY,NODPOS)
!
             NLZ(ILOC,GPOI)&
     &        =SPECFU(NDIFF,LXGP,INOD,NDNOD,IPOLY,NODPOS)&
     &        *SPECFU(NDIFF,LYGP,JNOD,NDNOD,IPOLY,NODPOS)&
     &        *SPECFU(DIFF,LZGP, KNOD,NDNOD,IPOLY,NODPOS)
!
      end do ! Was loop 20
      end do ! Was loop 20
      end do ! Was loop 20
      end do ! Was loop 101
      end do ! Was loop 101
      end do ! Was loop 101
!
!
! Find the roots of the quadrature points and nodes
! also get the weights. 
       ewrite(3,*)'2about to go into inside GTROOT IPOLY,IQADRA',IPOLY,IQADRA
        CALL GTROOT(IPOLY,IQADRA,WEIT,NODPOS,QUAPOS,NDGI,NMDNOD)
       ewrite(3,*)'2out of GTROOT'
      do  IGR=1,NDGI! Was loop 102
      do  IGQ=1,NDGI! Was loop 102
      do  IGP=1,NDGI! Was loop 102
         GPOI=IGP + (IGQ-1)*NDGI + (IGR-1)*NDGI*NDGI
!
           LXGP=QUAPOS(IGP)
           LYGP=QUAPOS(IGQ)
           LZGP=QUAPOS(IGR)

      do  KNOD=1,NMDNOD! Was loop 30
      do  JNOD=1,NMDNOD! Was loop 30
      do  INOD=1,NMDNOD! Was loop 30
           ILOC=INOD + (JNOD-1)*NMDNOD + (KNOD-1)*NMDNOD*NMDNOD
!
             M(ILOC,GPOI)&
     &        =SPECFU(NDIFF,LXGP,INOD,NMDNOD,IPOLY,NODPOS)&
     &        *SPECFU(NDIFF,LYGP,JNOD,NMDNOD,IPOLY,NODPOS)&
     &        *SPECFU(NDIFF,LZGP,KNOD,NMDNOD,IPOLY,NODPOS)
!
      end do ! Was loop 30
      end do ! Was loop 30
      end do ! Was loop 30
      end do ! Was loop 102
      end do ! Was loop 102
      end do ! Was loop 102
         ENDIF
!
         IF(D2) THEN
         NDGI  =INT((NGI**(1./2.))+0.1)
         NDNOD =INT((NLOC**(1./2.))+0.1)
         NMDNOD=INT((MLOC**(1./2.))+0.1)
!
! Find the roots of the quadrature points and nodes
! also get the weights. 
        CALL GTROOT(IPOLY,IQADRA,WEIT,NODPOS,QUAPOS,NDGI,NDNOD)
      do  IGQ=1,NDGI! Was loop 10
      do  IGP=1,NDGI! Was loop 10
         GPOI=IGP + (IGQ-1)*NDGI
!
           WEIGHT(GPOI)=WEIT(IGP)*WEIT(IGQ)
!
           LXGP=QUAPOS(IGP)
           LYGP=QUAPOS(IGQ)
! NB If TRUE in function RGPTWE then return the Gauss-pt weight 
! else return the Gauss-pt. 
!
      do  JNOD=1,NDNOD! Was loop 120
      do  INOD=1,NDNOD! Was loop 120
           ILOC=INOD + (JNOD-1)*NDNOD 
!
             N(ILOC,GPOI)&
     &        =SPECFU(NDIFF,LXGP,INOD,NDNOD,IPOLY,NODPOS)&
     &        *SPECFU(NDIFF,LYGP,JNOD,NDNOD,IPOLY,NODPOS)
!
             NLX(ILOC,GPOI)&
     &        =SPECFU(DIFF, LXGP,INOD,NDNOD,IPOLY,NODPOS)&
     &        *SPECFU(NDIFF,LYGP,JNOD,NDNOD,IPOLY,NODPOS)
!
             NLY(ILOC,GPOI)&
     &        =SPECFU(NDIFF,LXGP,INOD,NDNOD,IPOLY,NODPOS)&
     &        *SPECFU(DIFF, LYGP,JNOD,NDNOD,IPOLY,NODPOS)
!
      end do ! Was loop 120
      end do ! Was loop 120
      end do ! Was loop 10
      end do ! Was loop 10
!
! Find the roots of the quadrature points and nodes
! also get the weights. 
        CALL GTROOT(IPOLY,IQADRA,WEIT,NODPOS,QUAPOS,NDGI,NMDNOD)
      do  IGQ=1,NDGI! Was loop 11
      do  IGP=1,NDGI! Was loop 11
         GPOI=IGP + (IGQ-1)*NDGI 
           LXGP=QUAPOS(IGP)
           LYGP=QUAPOS(IGQ)
      do  JNOD=1,NMDNOD! Was loop 130
      do  INOD=1,NMDNOD! Was loop 130
           ILOC=INOD + (JNOD-1)*NMDNOD 
!
             M(ILOC,GPOI)&
     &        =SPECFU(NDIFF,LXGP,INOD,NMDNOD,IPOLY,NODPOS)&
     &        *SPECFU(NDIFF,LYGP,JNOD,NMDNOD,IPOLY,NODPOS)
!
      end do ! Was loop 130
      end do ! Was loop 130
!
      end do ! Was loop 11
      end do ! Was loop 11
         ENDIF
!
!
         IF((.NOT.D2).AND.(.NOT.D3)) THEN
         NDGI  =NGI
         NDNOD =NLOC
         NMDNOD=MLOC
!
! Find the roots of the quadrature points and nodes
! also get the weights. 
        CALL GTROOT(IPOLY,IQADRA,WEIT,NODPOS,QUAPOS,NDGI,NDNOD)
        ewrite(3,*)'NDGI,NDNOD,NLOC:',NDGI,NDNOD,NLOC
        ewrite(3,*)'WEIT(1:ndgi):',WEIT(1:ndgi)
        ewrite(3,*)'NODPOS(1:ndnod):',NODPOS(1:ndnod)
        ewrite(3,*)'QUAPOS(1:ndgi):',QUAPOS(1:ndgi)
      do  IGP=1,NDGI! Was loop 1000
         GPOI=IGP 
!
           WEIGHT(GPOI)=WEIT(IGP)
!
           LXGP=QUAPOS(IGP)
! NB If TRUE in function RGPTWE then return the Gauss-pt weight 
! else return the Gauss-pt. 
!
      do  INOD=1,NDNOD! Was loop 12000
           ILOC=INOD 
!
             N(ILOC,GPOI)&
     &        =SPECFU(NDIFF,LXGP,INOD,NDNOD,IPOLY,NODPOS)
!
             NLX(ILOC,GPOI)&
     &        =SPECFU(DIFF, LXGP,INOD,NDNOD,IPOLY,NODPOS)
         ewrite(3,*)'ILOC,GPOI,N(ILOC,GPOI),NLX(ILOC,GPOI):', &
                  ILOC,GPOI,N(ILOC,GPOI),NLX(ILOC,GPOI)
!
      end do ! Was loop 12000
      end do ! Was loop 1000
         ewrite(3,*)'n WEIGHT:',WEIGHT
!
! Find the roots of the quadrature points and nodes
! also get the weights. 
       ewrite(3,*)'this is for m which we dont care about:'
        CALL GTROOT(IPOLY,IQADRA,WEIT,NODPOS,QUAPOS,NDGI,NMDNOD)
      do  IGP=1,NDGI! Was loop 1100
         GPOI=IGP 
           LXGP=QUAPOS(IGP)
      do  INOD=1,NMDNOD! Was loop 13000
           ILOC=INOD
!
             M(ILOC,GPOI)&
     &        =SPECFU(NDIFF,LXGP,INOD,NMDNOD,IPOLY,NODPOS)
!
      end do ! Was loop 13000
!
      end do ! Was loop 1100
       ewrite(3,*)'...finished this is for m which we dont care about:'
         ENDIF
         END SUBROUTINE SPECTR



      
      SUBROUTINE GTROOT(IPOLY,IQADRA,WEIT,NODPOS,QUAPOS,NDGI,NDNOD)
      IMPLICIT NONE
      INTEGER IPOLY,IQADRA,NDGI,NDNOD
      REAL WEIT(NDGI),QUAPOS(NDGI),NODPOS(NDNOD)
      LOGICAL GETNDP
!     This sub returns the weights WEIT the quadrature points QUAPOS and 
!     the node points NODPOS. 
!     NODAL POISTIONS ******
!     NB if GETNDP then find the nodal positions
      GETNDP=.TRUE.
!     Compute standard Lagrange nodal points
      IF(IQADRA.EQ.1) CALL LAGROT(NODPOS,NODPOS,NDNOD,GETNDP)
!     Compute Chebyshev-Gauss-Lobatto nodal points.
      IF(IQADRA.EQ.2) CALL CHEROT(NODPOS,NODPOS,NDNOD,GETNDP)
      IF(IQADRA.EQ.3) CALL CHEROT(NODPOS,NODPOS,NDNOD,GETNDP)
!     Compute Legendre-Gauss-Lobatto nodal points.
      IF(IQADRA.EQ.4) CALL LEGROT(NODPOS,NODPOS,NDNOD,GETNDP)
!     
!     QUADRATURE************
      GETNDP=.FALSE.
!     Compute standard Gauss quadrature. weits and points
      IF(IQADRA.EQ.1) CALL LAGROT(WEIT,QUAPOS,NDGI,GETNDP)
!     Compute Chebyshev-Gauss-Lobatto quadrature.
      IF(IQADRA.EQ.2) CALL CHEROT(WEIT,QUAPOS,NDGI,GETNDP)
      IF(IQADRA.EQ.3) CALL CHEROT(WEIT,QUAPOS,NDGI,GETNDP)
!     Compute Legendre-Gauss-Lobatto quadrature.
      IF(IQADRA.EQ.4) CALL LEGROT(WEIT,QUAPOS,NDGI,GETNDP)
      END SUBROUTINE GTROOT
      



      REAL FUNCTION SPECFU(DIFF,LXGP,INOD,NDNOD,IPOLY,NODPOS)
      LOGICAL DIFF
      INTEGER INOD,NDNOD,IPOLY
      REAL LXGP,NODPOS(NDNOD)
!     INOD contains the node at which the polynomial is associated with
!     LXGP is the position at which the polynomial is to be avaluated.\
!     If(DIFF) then find the D poly/DX. 
!     
      IF(IPOLY.EQ.1) SPECFU=LAGRAN(DIFF,LXGP,INOD,NDNOD,NODPOS)
!     
      IF(IPOLY.EQ.2) SPECFU=CHEBY1(DIFF,LXGP,INOD,NDNOD,NODPOS)
!     
      IF(IPOLY.EQ.3) SPECFU=CHEBY2(DIFF,LXGP,INOD,NDNOD,NODPOS)
!     
      IF(IPOLY.EQ.4) SPECFU=LEGEND(DIFF,LXGP,INOD,NDNOD,NODPOS)
      END FUNCTION SPECFU



      SUBROUTINE CHEROT(WEIT,QUAPOS,NDGI,GETNDP)
      IMPLICIT NONE
      INTEGER NDGI
      REAL PIE
      PARAMETER(PIE=3.141569254)
!     This computes the weight and points for Chebyshev-Gauss-Lobatto quadrature.
!     See page 67 of:Spectral Methods in Fluid Dynamics, C.Canuto
!     IF(GETNDP) then get the POSITION OF THE NODES 
!     AND DONT BOTHER WITH THE WEITS.
      REAL WEIT(NDGI),QUAPOS(NDGI)
      LOGICAL GETNDP
      INTEGER IG,J
!     
      IF(.NOT.GETNDP) THEN
!     THE WEIGHTS...
         WEIT(1)       =PIE/(2.*(NDGI-1))
         WEIT(NDGI)=PIE/(2.*(NDGI-1))
      do IG=2,NDGI-1
            WEIT(IG)   =PIE/(NDGI-1)
         END DO
      ENDIF
!     
!     The quad points...
      do IG=1,NDGI
         J=IG-1
         QUAPOS(IG)=COS(PIE*REAL(J)/REAL(NDGI-1))
      END DO
      END SUBROUTINE CHEROT

      SUBROUTINE LEGROT(WEIT,QUAPOS,NDGI,GETNDP)
      IMPLICIT NONE
!     This computes the weight and points for Chebyshev-Gauss-Lobatto quadrature.
!     See page 69 of:Spectral Methods in Fluid Dynamics, C.Canuto
!     IF(GETNDP) then get the POSITION OF THE NODES 
!     AND DONT BOTHER WITH THE WEITS.
      INTEGER NDGI
      REAL WEIT(NDGI),QUAPOS(NDGI)
      LOGICAL GETNDP
      INTEGER N,IG
!     Work out the root's i.e the quad positions first.
      CALL LROOTS(QUAPOS,NDGI)
      IF(.NOT.GETNDP) THEN
!     THE WEIGHTS...
         N=NDGI-1
      do IG=1,NDGI
            WEIT(IG)=2./(REAL(N*(N+1))*PLEGEN(QUAPOS(IG),NDGI-1)**2)
         END DO
      ENDIF
      END SUBROUTINE LEGROT

  real function PLEGEN(LX,K)
      REAL LX
      REAL R
      INTEGER K,L
      R=0.
      DO L=0,INT(K/2)
         R=R+((-1)**L)*binomial_coefficient(K,L)*binomial_coefficient(2*K-2*L,K)*LX**(K-2*L)
      END DO
      PLEGEN=R/REAL(2**K)
  end function plegen

   
  real function binomial_coefficient(K,L)
    !<! Calculate binomial coefficients
    
      INTEGER K,L
      binomial_coefficient=factorial(K)/(factorial(L)*factorial(K-L))

  end function binomial_coefficient

      SUBROUTINE LROOTS(QUAPOS,NDGI)
      IMPLICIT NONE
      INTEGER NDGI
      REAL QUAPOS(NDGI)
      REAL ALPHA,BETA,RKEEP
      INTEGER N,I
!     This sub works out the Gauss-Lobatto-Legendre roots.
      ALPHA = 0.
      BETA   = 0.
!     
      N=NDGI-1
      CALL JACOBL(N,ALPHA,BETA,QUAPOS)
!     Now reverse ordering. 
      do I=1,INT(0.1+ NDGI/2)
         RKEEP=QUAPOS(I)
         QUAPOS(I)=QUAPOS(NDGI+1-I)
         QUAPOS(NDGI+1-I)=RKEEP
      END DO
      END SUBROUTINE LROOTS

      REAL FUNCTION CHEBY1(DIFF,LX,INOD,NDNOD,NODPOS)
      IMPLICIT NONE
      INTEGER NDNOD,INOD
      REAL NODPOS(NDNOD),LX,RNX
      LOGICAL DIFF,DIFF2
!     If DIFF then returns the spectral function DIFFERENTIATED W.R.T X 
!     associated.
!     This function returns the spectral function associated 
!     with node INOD at POINT LX 
!     NDNOD=no of nodes in 1-D.
!     NDGI=no of Gauss pts in 1-D.
!     NB The nodes are at the points COS(pie*J/2.) j=0,..,ndgi-1
      REAL CJ,CN,HX,XPT
      INTEGER NX,N
      
      NX=NDNOD-1
      RNX=REAL(NX)
      DIFF2=.FALSE.
!     
      CJ=1.
      IF((INOD.EQ.1).OR.(INOD.EQ.NDNOD)) CJ=2.
      HX=0.
!     
      CN=2.
      N=0
      XPT=NODPOS(INOD)
      HX=HX + ( (2./RNX)/(CJ*CN) )*TCHEB(N,XPT,.FALSE.,DIFF2)&
     &     *TCHEB(N,LX,DIFF,DIFF2)
      CN=1.
      do  N=1,NX-1! Was loop 10
         HX=HX + ( (2./RNX)/(CJ*CN) )*TCHEB(N,XPT,.FALSE.,DIFF2)&
     &        *TCHEB(N,LX,DIFF,DIFF2)
      end do ! Was loop 10
      CN=2.
      N=NX
      HX=HX + ( (2./RNX)/(CJ*CN))*TCHEB(N,XPT,.FALSE.,DIFF2)&
     &     *TCHEB(N,LX,DIFF,DIFF2)
      CHEBY1=HX
      END FUNCTION CHEBY1
      


      
      REAL FUNCTION CHEBY2(DIFF,LX,INOD,NDNOD,NODPOS)
      IMPLICIT NONE
      INTEGER NDNOD,INOD
      REAL NODPOS(NDNOD)
      REAL LX,R,RR,RCONST
      LOGICAL DIFF
!     If DIFF then returns the spectral function DIFFERENTIATED W.R.T X 
!     associated.
!     This function returns the spectral function associated 
!     with node INOD at POINT LX 
!     NDNOD=no of nodes in 1-D.
!     NDGI=no of Gauss pts in 1-D.
!     NB The nodes are at the points COS(pie*J/2.) j=0,..,ndgi-1
      REAL CJ,CN,HX,XPT
      IF(.NOT.DIFF) THEN
         IF(ABS(LX-NODPOS(INOD)).LT.1.E-10) THEN
            CHEBY2=1.
         ELSE
            R=(-1.)**(INOD)*(1.-LX**2)&
     &       *TCHEB(NDNOD-1,LX,.TRUE.,.FALSE.)
             CJ=1.
             IF((INOD.EQ.1).OR.(INOD.EQ.NDNOD)) CJ=2.
               R=R/(CJ*(NDNOD-1)**2 *(LX-NODPOS(INOD)) )
               CHEBY2=R
             ENDIF
           ELSE
             IF(ABS(LX-NODPOS(INOD)).LT.1.E-10) THEN
               IF(ABS(LX+1.).LT.1.E-10) THEN
                 CHEBY2= (2.*(NDNOD-1)**2+1)/6.
               ELSE
                 IF(ABS(LX-1.).LT.1.E-10) THEN
                   CHEBY2=-(2.*(NDNOD-1)**2+1)/6.
                 ELSE
                   CHEBY2=-LX/(2.*(1.-LX**2))
                 ENDIF
               ENDIF
             ELSE
               R=(-1.)**(INOD)*(1.-LX**2)&
     &          *TCHEB(NDNOD-1,LX,.FALSE.,.TRUE.)
               CJ=1.
               IF((INOD.EQ.1).OR.(INOD.EQ.NDNOD)) CJ=2.
               R=R/(CJ*(NDNOD-1)**2 *(LX-NODPOS(INOD)) )
               RR=(-1.)**(INOD)*TCHEB(NDNOD-1,LX,.TRUE.,.FALSE.)&
     &           /(CJ*REAL(NDNOD-1)**2)
               RCONST=-(2*LX+(1.-LX**2)/(LX-NODPOS(INOD) ) )&
     &              /(LX-NODPOS(INOD))
               RR=RR*RCONST
               CHEBY2=RR
             ENDIF
           ENDIF            
         END FUNCTION CHEBY2



      REAL FUNCTION TCHEB(N,XPT,DIFF,DIFF2)
!      use math_utilities
      IMPLICIT NONE
      LOGICAL DIFF,DIFF2
      REAL DDT,DT,T,XPT,TTEMP,TM1,DTM1,DTEMP,R,RR
      INTEGER K,L,N
      INTEGER NI
! If DIFF then return the n'th Chebyshef polynomial 
! differentiated w.r.t x.
! If DIFF2 then form the 2'nd derivative. 
! This sub returns the value of the K'th Chebyshef polynomial at a 
! point XPT.
!
! This formula can be found in:Spectral Methods for Fluid Dynamics, page 66
        K=N
!
       IF((.NOT.DIFF).AND.(.NOT.DIFF2)) THEN
         IF(N.EQ.0) T=1.
         IF(N.EQ.1) T=XPT
         IF(N.GT.1) THEN
          TM1=1.
          T=XPT
      do  NI=2,N! Was loop 110
            TTEMP=T
            T=2.*XPT*TTEMP - TM1
            TM1=TTEMP
      end do ! Was loop 110
         ENDIF
         TCHEB=T
        ENDIF
!
       IF(DIFF.AND.(.NOT.DIFF2)) THEN
! This part forms the differential w.r.t x.
         IF(N.EQ.0) DT=0.
         IF(N.EQ.1) THEN 
           T=XPT
           DT=1.
         ENDIF
         IF(N.EQ.2) THEN 
           T=2.*XPT*XPT -1
           DT=4.*XPT 
         ENDIF
         IF(N.GT.2) THEN
           TM1=XPT
           DTM1=1.
           T=2.*XPT*XPT -1
           DT=4.*XPT 
      do  NI=2,N! Was loop 10
            TTEMP=T
            T=2.*XPT*TTEMP - TM1
            TM1=TTEMP
!
            DTEMP=DT
            DT=2.*XPT*DTEMP+2.*TTEMP-DTM1
            DTM1=DTEMP
      end do ! Was loop 10
         ENDIF
         TCHEB=DT
         ENDIF
!
         IF(DIFF2) THEN
            IF(N.LE.1) THEN
               DDT=0.
            ELSE
              R=0.
      do  L=0,INT(0.1+K/2)! Was loop 50
                RR=(-1.)**K*factorial(K-L-1)/(factorial(L)*factorial(K-2*L))
                RR=RR*2**(K-2*L)
! The following is for 2'nd derivative. 
                RR=RR*(K-2.*L)*(K-2.*L-1)*XPT**(K-2*L-2)
                R=R+RR
      end do ! Was loop 50
              DDT=R*REAL(K)*0.5
            ENDIF
           TCHEB=DDT
         ENDIF
         END FUNCTION TCHEB


    recursive function factorial(n) result(f)
      ! Calculate n!
      integer :: f
      integer, intent(in) :: n

      if (n==0) then
         f=1
      else
         f=n*factorial(n-1)
      end if

    end function factorial



      REAL FUNCTION LEGEND(DIFF,LX,INOD,NDNOD,NODPOS)
      INTEGER INOD,NDNOD
      REAL NODPOS(NDNOD)
      REAL LX
      LOGICAL DIFF
!     If DIFF then returns the spectral function DIFFERENTIATED W.R.T X 
!     associated.
!     This function returns the spectral function associated 
!     with node INOD at POINT LX 
!     NDNOD=no of nodes in 1-D.
!     NDGI=no of Gauss pts in 1-D.
!     NB The nodes are at the points COS(pie*J/2.) j=0,..,ndgi-1
      REAL CJ,CN,HX,XPT
      LEGEND=1.
      END FUNCTION LEGEND


    SUBROUTINE DETERMIN_SLOCLIST( CV_SLOCLIST, CV_NLOC, CV_SNLOC, NFACE,  &
         ndim, cv_ele_type )
      ! determine CV_SLOCLIST
      IMPLICIT NONE
      INTEGER, intent( in ) :: CV_NLOC, CV_SNLOC, NFACE, ndim, cv_ele_type
      INTEGER, DIMENSION( NFACE, CV_SNLOC ), intent( inout ) :: CV_SLOCLIST
      ! Local variables
      INTEGER :: IFACE, quad_cv_siloc, quad_cv_iloc, NPOLY, IP, JP, KP
      LOGICAL :: FOUND

      CONDITIONAL_DIMENSION: IF(NDIM==1) THEN
         ! 1d: 
         IF(NFACE.NE.2) THEN
            EWRITE(3,*) 'NFACE not correct NFACE=',NFACE
            STOP 4332
         END IF
         CV_SLOCLIST(1,1)=1
         CV_SLOCLIST(2,1)=CV_NLOC
      ELSE IF(NDIM==2) THEN
         ! 2d: 
         ! linear triangle: 
         IF(CV_NLOC==3) THEN
            IF(NFACE/=3) THEN
               EWRITE(3,*) 'NFACE not correct NFACE=',NFACE
               STOP 4333
            END IF
            CV_SLOCLIST(1,1)=1
            CV_SLOCLIST(1,2)=2
            CV_SLOCLIST(2,1)=1
            CV_SLOCLIST(2,2)=3
            CV_SLOCLIST(3,1)=2
            CV_SLOCLIST(3,2)=3
            ! linear quad: 
         ELSE IF(CV_NLOC==4) THEN
            IF(NFACE/=4) THEN
               EWRITE(3,*) 'NFACE not correct NFACE=',NFACE
               STOP 4334
            END IF
            CV_SLOCLIST(1,1)=1
            CV_SLOCLIST(1,2)=3
            CV_SLOCLIST(2,1)=2
            CV_SLOCLIST(2,2)=4
            CV_SLOCLIST(3,1)=1
            CV_SLOCLIST(3,2)=2
            CV_SLOCLIST(4,1)=1
            CV_SLOCLIST(4,2)=2
            ! quadratic triangle: 
         ELSE IF(CV_NLOC==6) THEN
            IF(NFACE/=3) THEN
               EWRITE(3,*) 'NFACE not correct NFACE=',NFACE
               STOP 4335
            END IF
            CV_SLOCLIST(1,1)=1
            CV_SLOCLIST(1,2)=2
            CV_SLOCLIST(1,3)=3
            CV_SLOCLIST(2,1)=1
            CV_SLOCLIST(2,2)=4
            CV_SLOCLIST(2,3)=6
            CV_SLOCLIST(3,1)=3
            CV_SLOCLIST(3,2)=5
            CV_SLOCLIST(3,3)=6
            ! quadratic quad: 
         ELSE IF(CV_NLOC==9) THEN
            IF(NFACE/=4) THEN
               EWRITE(3,*) 'NFACE not correct NFACE=',NFACE
               STOP 4336
            ENDIF
            CV_SLOCLIST(1,1)=1
            CV_SLOCLIST(1,2)=4
            CV_SLOCLIST(1,3)=7
            CV_SLOCLIST(2,1)=3
            CV_SLOCLIST(2,2)=6
            CV_SLOCLIST(2,3)=9
            CV_SLOCLIST(3,1)=1
            CV_SLOCLIST(3,2)=2
            CV_SLOCLIST(3,3)=3
            CV_SLOCLIST(4,1)=7
            CV_SLOCLIST(4,2)=8
            CV_SLOCLIST(4,3)=9

         END IF

      ELSE IF(NDIM==3) THEN
         ! 3d: 
         ! linear triangle: 
         IF(CV_NLOC==4) THEN
            IF(NFACE/=4) THEN
               EWRITE(3,*) 'NFACE not correct NFACE=',NFACE
               STOP 4337
            END IF
            CV_SLOCLIST(1,1)=1
            CV_SLOCLIST(1,2)=2
            CV_SLOCLIST(1,3)=3
            CV_SLOCLIST(2,1)=1
            CV_SLOCLIST(2,2)=4
            CV_SLOCLIST(2,3)=2
            CV_SLOCLIST(3,1)=1
            CV_SLOCLIST(3,2)=3
            CV_SLOCLIST(3,3)=4
            CV_SLOCLIST(4,1)=2
            CV_SLOCLIST(4,2)=3
            CV_SLOCLIST(4,3)=4
            ! quadratic triangle: 
         ELSE IF(CV_NLOC==10) THEN
            IF(NFACE/=4) THEN
               EWRITE(3,*) 'NFACE not correct NFACE=',NFACE
               STOP 4338
            END IF
            CV_SLOCLIST(1,1)=1
            CV_SLOCLIST(1,2)=2
            CV_SLOCLIST(1,3)=3
            CV_SLOCLIST(1,4)=4
            CV_SLOCLIST(1,5)=5
            CV_SLOCLIST(1,6)=6

            CV_SLOCLIST(2,1)=1
            CV_SLOCLIST(2,2)=4
            CV_SLOCLIST(2,3)=6
            CV_SLOCLIST(2,4)=7
            CV_SLOCLIST(2,5)=9
            CV_SLOCLIST(2,6)=10

            CV_SLOCLIST(3,1)=6
            CV_SLOCLIST(3,2)=5
            CV_SLOCLIST(3,3)=3
            CV_SLOCLIST(3,4)=9
            CV_SLOCLIST(3,5)=8
            CV_SLOCLIST(3,6)=10

            CV_SLOCLIST(4,1)=3
            CV_SLOCLIST(4,2)=2
            CV_SLOCLIST(4,3)=1
            CV_SLOCLIST(4,4)=8
            CV_SLOCLIST(4,5)=7
            CV_SLOCLIST(4,6)=10
            ! general hex: 
         ELSE 
            IF(NFACE/=6) THEN
               EWRITE(3,*) 'NFACE not correct NFACE=',NFACE
               STOP 4339
            ENDIF
            npoly=INT(REAL(CV_NLOC+0.01)**0.333333)
            do iface=1,nface
               quad_cv_siloc =0
               Loop_Poly2d_2_1: do ip = 1, npoly
                  Loop_Poly2d_2_2: do jp = 1, npoly
                     Loop_Poly2d_2_3: do kp = 1, npoly
                        quad_cv_iloc = ( kp - 1 ) * npoly * npoly + &
                             ( jp - 1 ) * npoly + ip 
                        found = .false.
                        Select Case( iface )
                        case( 1 )
                           if( ip == 1 ) found = .true.
                        case( 2 )  
                           if( ip == npoly ) found = .true.
                        case( 3 )
                           if( jp == 1 ) found = .true.
                        case( 4 )
                           if( jp == npoly ) found = .true.
                        case( 5 )
                           if( kp == 1 ) found = .true.
                        case( 6 )
                           if( kp == npoly ) found = .true.
                        case default; FLExit( "Wrong integer for IFACE" )
                        end Select
                        if( found ) then
                           quad_cv_siloc = quad_cv_siloc + 1
                           CV_SLOCLIST(iface,quad_cv_siloc)=quad_cv_iloc

                        end if
                     end do Loop_Poly2d_2_3
                  end do Loop_Poly2d_2_2
               end do Loop_Poly2d_2_1
            end do

         END IF
      END IF CONDITIONAL_DIMENSION

      RETURN
    END SUBROUTINE DETERMIN_SLOCLIST



      SUBROUTINE JACOBL(N,ALPHA,BETA,XJAC)
      IMPLICIT NONE
!     COMPUTES THE GAUSS-LOBATTO COLLOCATION POINTS FOR JACOBI POLYNOMIALS
!     
!     N:       DEGREE OF APPROXIMATION
!     ALPHA:   PARAMETER IN JACOBI WEIGHT
!     BETA:    PARAMETER IN JACOBI WEIGHT
!     
!     XJAC:    OUTPUT ARRAY WITH THE GAUSS-LOBATTO ROOTS
!     THEY ARE ORDERED FROM LARGEST (+1.0) TO SMALLEST (-1.0)
!     
      INTEGER N
      REAL ALPHA,BETA
!      IMPLICIT REAL(A-H,O-Z)
!      REAL XJAC(1)
      REAL XJAC(N+1)
      REAL ALP,BET,RV
      REAL PNP1P,PDNP1P,PNP,PDNP,PNM1P,PDNM1,PNP1M,PDNP1M,PNM,PDNM,PNM1M
      REAL DET,RP,RM,A,B,DTH,CD,SD,CS,SS,X,PNP1,PDNP1,PN,PDN,PNM1,POLY
      REAL PDER,RECSUM,DELX,CSSAVE
      INTEGER NP,NH,J,K,JM,I,NPP
      COMMON /JACPAR/ALP,BET,RV
      INTEGER KSTOP
      DATA KSTOP/10/
      REAL EPS
      DATA EPS/1.0E-12/
      ALP = ALPHA
      BET =BETA
      RV = 1 + ALP
      NP = N+1
!
!  COMPUTE THE PARAMETERS IN THE POLYNOMIAL WHOSE ROOTS ARE DESIRED
!
      CALL JACOBF(NP,PNP1P,PDNP1P,PNP,PDNP,PNM1P,PDNM1,1.0)
      CALL JACOBF(NP,PNP1M,PDNP1M,PNM,PDNM,PNM1M,PDNM1,-1.0)
      DET = PNP*PNM1M-PNM*PNM1P
      RP = -PNP1P
      RM = -PNP1M
      A = (RP*PNM1M-RM*PNM1P)/DET
      B = (RM*PNP-RP*PNM)/DET
!
      XJAC(1) = 1.0
      NH = (N+1)/2
!
!  SET-UP RECURSION RELATION FOR INITIAL GUESS FOR THE ROOTS
!
      DTH = 3.14159265/(2*N+1)
      CD = COS(2.*DTH)
      SD = SIN(2.*DTH)
      CS = COS(DTH)
      SS = SIN(DTH)
!
!  COMPUTE THE FIRST HALF OF THE ROOTS BY POLYNOMIAL DEFLATION
!
      do  J=2,NH! Was loop 39
         X = CS
      do  K=1,KSTOP! Was loop 29
            CALL JACOBF(NP,PNP1,PDNP1,PN,PDN,PNM1,PDNM1,X)
            POLY = PNP1+A*PN+B*PNM1
            PDER = PDNP1+A*PDN+B*PDNM1
            RECSUM = 0.0
            JM = J-1
      do  I=1,JM! Was loop 27
               RECSUM = RECSUM+1.0/(X-XJAC(I))
      end do ! Was loop 27
28          CONTINUE
            DELX = -POLY/(PDER-RECSUM*POLY)
            X = X+DELX
            IF(ABS(DELX) .LT. EPS) GO TO 30
      end do ! Was loop 29
30       CONTINUE
         XJAC(J) = X
         CSSAVE = CS*CD-SS*SD
         SS = CS*SD+SS*CD
         CS = CSSAVE
      end do ! Was loop 39
      XJAC(NP) = -1.0
      NPP = N+2
!
! USE SYMMETRY FOR SECOND HALF OF THE ROOTS
!
      do  I=2,NH! Was loop 49
         XJAC(NPP-I) = -XJAC(I)
      end do ! Was loop 49
      IF(N .NE. 2*(N/2)) RETURN
      XJAC(NH+1) = 0.0
      RETURN
      END SUBROUTINE JACOBL
      
      
      SUBROUTINE JACOBF(N,POLY,PDER,POLYM1,PDERM1,POLYM2,PDERM2,X)
      IMPLICIT NONE
!     
!     COMPUTES THE JACOBI POLYNOMIAL (POLY) AND ITS DERIVATIVE
!     (PDER) OF DEGREE  N  AT  X
!     
      INTEGER N
      REAL APB,POLY,PDER,POLYM1,PDERM1,POLYM2,PDERM2,X
!     IMPLICIT REAL(A-H,O-Z)
      COMMON /JACPAR/ALP,BET,RV
      REAL ALP,BET,RV,POLYLST,PDERLST,A1,A2,B3,A3,A4
      REAL POLYN,PDERN,PSAVE,PDSAVE
      INTEGER K
      APB = ALP+BET
      POLY = 1.0
      PDER = 0.0
      IF(N .EQ. 0) RETURN
      POLYLST = POLY
      PDERLST = PDER
      POLY = RV * X
      PDER = RV
      IF(N .EQ. 1) RETURN
      do K=2,N
         A1 = 2.*K*(K+APB)*(2.*K+APB-2.)
         A2 = (2.*K+APB-1.)*(ALP**2-BET**2)
         B3 = (2.*K+APB-2.)
         A3 = B3*(B3+1.)*(B3+2.)
         A4 = 2.*(K+ALP-1)*(K+BET-1.)*(2.*K+APB)
         POLYN = ((A2+A3*X)*POLY-A4*POLYLST)*A1
         PDERN = ((A2+A3*X)*PDER-A4*PDERLST+A3*POLY)*A1
         PSAVE = POLYLST
         PDSAVE = PDERLST
         POLYLST = POLY
         POLY = POLYN
         PDERLST = PDER
         PDER = PDERN
      END DO
      POLYM1 = POLYLST
      PDERM1 = PDERLST
      POLYM2 = PSAVE
      PDERM2 = PDSAVE
      RETURN
      END SUBROUTINE JACOBF




  real function Volume_TetHex( hexs, &
       x1, x2, x3, x4, &
       y1, y2, y3, y4, &
       z1, z2, z3, z4 )
    implicit none
    logical :: hexs
    real :: x1, x2, x3, x4, &
       y1, y2, y3, y4, &
       z1, z2, z3, z4 
    integer :: iloc

    Volume_TetHex = &
         
         ( x4 - x1 ) * ( &
         ( y2 - y1 ) * ( z3 - z1 ) - ( z2 - z1 ) * ( y3 - y1 ) &
         ) + &
         
         ( y4 - y1 ) * ( &
         ( z2 - z1 ) * ( x3 - x1 ) - ( x2 - x1 ) * ( z3 - z1 ) &
         ) + &
         
         ( z4 - z1 ) * ( &
         ( x2 - x1 ) * ( y3 - y1 ) - ( y2 - y1 ) * ( x3 - x1 ) &
         )

    if( .not. hexs ) Volume_TetHex = Volume_TetHex / 6.


    return
  end function Volume_TetHex


  end module shape_functions_NDim

