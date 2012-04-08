  
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

      Conditional_LOWQUA: if( lowqua .or. (ngi == 4) ) then
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

      Conditional_LOWQUA: if( lowqua .or. (ngi == 8) ) then
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
         cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )
      implicit none
      integer, intent( in ) :: ndim, cv_ele_type, cv_nloc, u_nloc
! If QUAD_OVER_WHOLE_ELE=.true. then dont divide element into CV's to form quadrature.
      logical, intent( in ) :: QUAD_OVER_WHOLE_ELE
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
            cv_ngi = 48 ! 36
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
      integer, dimension( : ), allocatable :: x_ndgln, fem_nod
      real, dimension( : ), allocatable :: lx, ly, lz, x, y, z, cvweigh_dummy
      integer, parameter :: max_totele = 1000, max_x_nonods = 1000
      logical :: d1, dcyl, d3
      integer :: ele, quad_cv_ngi, quad_cv_nloc, totele, x_nonods, &
           cv_gj,cv_gk,cv_iloc

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
      allocate( lx( max_x_nonods ) )
      allocate( ly( max_x_nonods ) )
      allocate( lz( max_x_nonods ) )
      allocate( x( max_x_nonods ))
      allocate( y( max_x_nonods ))
      allocate( z( max_x_nonods ))
      allocate( fem_nod( max_x_nonods ) )
      allocate( x_ndgln( max_totele * quad_cv_nloc ) )
      allocate(cvweigh_dummy(cv_ngi)) 

      ! Get the x_ndgln for the nodes of the triangle or tet or hex/quad super-elements:
      call Compute_XNDGLN_TriTetQuadHex( cv_ele_type, &
           max_totele, max_x_nonods, quad_cv_nloc, &
           totele, x_nonods, &
           x_ndgln, lx, ly, lz, x, y, z, fem_nod )

      ! Compute the shape functions using these quadrilaterals/hexs:
      ! For pressure:
      call shape_tri_tet( cv_ele_type, cv_nloc, &
           cv_ele_type, ndim, totele, cv_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, x, y, z, lx, ly, lz, &
           n, nlx, nly, nlz, cvweigh )

      ! Compute cvn: 
      call Calc_CVN_TriTetQuadHex( cv_ele_type, totele, cv_nloc, cv_ngi, x_nonods, &
           quad_cv_nloc, x_ndgln, fem_nod, cvn )

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



    subroutine Compute_XNDGLN_TriTetQuadHex( cv_ele_type, &
         max_totele, max_x_nonods, quad_cv_nloc, &
         totele, x_nonods, &
         x_ndgln, lx, ly, lz, x, y, z, fem_nod )
      ! Get the x_ndgln for the nodes of triangles or tetrahedra
      implicit none
      integer, intent( in ) :: cv_ele_type, max_totele, max_x_nonods, quad_cv_nloc
      integer, intent( inout ) :: totele, x_nonods
      real, dimension( max_x_nonods ), intent( inout ) :: lx, ly, lz, x, y, z
      integer, dimension( max_x_nonods ), intent( inout ) :: fem_nod
      integer, dimension( max_totele * quad_cv_nloc ), intent( inout ) :: x_ndgln
      ! Local variables
      integer, dimension( : ), allocatable :: x_ndgln2, x_ndgln_big
      real, dimension( : ), allocatable :: x2, y2, z2
      integer :: quad_u_loc_dummy, ele, quad_cv_iloc, quad_cv_siloc, cv_iloc, &
           x_nonods2, npoly, npoly_x, npoly_y, npoly_z, nelex, neley, nelez, &
           elex, eley, elez, ip, jp, kp, xcount, xnod,  xnod2, quad_cv_nloc2
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
         !h_scale = 2.70240031 ! Scaling factor to give a unity area of the local triangle
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

         x_nonods = 16
         totele = 9
         !h_scale = 2.70240031 ! Scaling factor to give a unity area of the local triangle
         h_scale =  1.5196713713031851 ! Scaling factor to give a unity area of the local triangle

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

      case( 7 ) ! Linear Tetrahedra
         x_nonods = 15
         totele = 4
         h_scale = 2.04 ! Scaling factor to give a unity area of the local tetrahedron

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

      case( 8 ) ! Quadratic Tetrahedra 
         call get_x_ndgln_qtet( max_x_nonods, max_totele, quad_cv_nloc, &
              x_nonods, totele, fem_nod, x_ndgln, &
              x, y, z, lx, ly, lz )

      case default; FLExit( "Wrong integer for CV_ELE_TYPE" )
      end Select Conditional_ElementTypes

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

      totele2 = totele
      x_nonods2 =  no_of_nodes_in_faces * no_faces * totele2 + totele2
      triangle_totele = totele2 * tet_totele
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
      integer, parameter :: max_totele = 1000, max_x_nonods = 1000
      integer, dimension( : ), allocatable :: x_ndgln, fem_nod
      integer, dimension( :,: ), allocatable :: cv_neiloc_cells_dummy
      real, dimension( : ), allocatable :: lx, ly, lz, x, y, z, scvfeweigh_dummy
      logical :: d1, dcyl, d3
      integer :: x_nonods, totele, ele, cv_iloc, quad_cv_ngi, quad_cv_nloc

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
      allocate( cv_neiloc_cells_dummy( cv_nloc,scvngi ) )

      ! Get the x_ndgln for the nodes of triangles, tetrahedra, quadrilaterals or hexahedra
      ! super-elements:
      x = 0. ; y = 0. ; z = 0. ; lx = 0. ; ly = 0. ; lz = 0.
      call Compute_XNDGLN_TriTetQuadHex( cv_ele_type, &
           max_totele, max_x_nonods, quad_cv_nloc, &
           totele, x_nonods, &
           x_ndgln, lx, ly, lz, x, y, z, fem_nod )

      ewrite(3,*)'cv_ele_type, totele, x_nonods:', cv_ele_type, totele, x_nonods
      ewrite(3,*)'x_ndgln:'
      do ele = 1, totele
         ewrite(3,*)( x_ndgln( ( ele - 1 ) * quad_cv_nloc + cv_iloc ), &
              cv_iloc = 1, quad_cv_nloc )
      end do
      ewrite(3,*)'x:', ( x( cv_iloc ), cv_iloc = 1, x_nonods )
      ewrite(3,*)'y:', ( y( cv_iloc ), cv_iloc = 1, x_nonods )
      ewrite(3,*)'z:', ( z( cv_iloc ), cv_iloc = 1, x_nonods )
      ewrite(3,*)'lx:', ( lx( cv_iloc ), cv_iloc = 1, 3 )
      ewrite(3,*)'ly:', ( ly( cv_iloc ), cv_iloc = 1, 3 )
      ewrite(3,*)'lz:', ( lz( cv_iloc ), cv_iloc = 1, 3 )

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
      if(u_nloc==1) then ! a constant basis function... 
        sufen=1.0 
        sufenlx=0.0 
        sufenly=0.0
        sufenlz=0.0 
        sufenslx=0.0
        sufensly=0.0
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
      integer, dimension( cv_nloc_cells, scvngi ), intent( inout ) :: cvfem_neiloc
      real, dimension( x_nonods ), intent( in ) :: x, y, z
      real, dimension( quad_cv_nloc ), intent( in ) :: lx, ly, lz ! corner nodes
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
           xnodi, xnodj, nodi, nodj, cv_iloc_cells, cv_jloc_cells
      real :: xgi, ygi, zgi, volume, sarea, normx, normy, normz, d2_quad
      real :: half_side_length

      ewrite(3,*)'In shape_tri_tet s'

      sn=0.0
      snlx=0.0
      snly=0.0
      snlz=0.0
      sufnlx=0.0
      sufnly=0.0

      d3 = ( ndim == 3 )
      d2 = ( ndim == 2 )
      d1 = .false.
      dcyl = .false. 

      ! Compute some dummy variables
      ! npoly=2 for linear 2 node elements, =3 for quadratic 3 node elements in 1D. 
      call dummy_tri_tet( d1, d3, quad_cv_ngi, quad_cv_nloc, &
           dummy_sngi, dummy_snloc, nwicel, & 
           cv_nloc_cells, scvngi, totele, quad_u_loc_dummy, &
           mloc, dummy_smloc, lowqua, npoly )

      quad_cv_snloc = npoly ** ( ndim - 1 )
      !      quad_cv_sngi = quad_cv_snloc
      quad_cv_sngi = max(1,npoly-1) ** ( ndim - 1 )

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

      allocate( l1( stotel * quad_cv_sngi * totele ) )
      allocate( l2( stotel * quad_cv_sngi * totele ) )
      allocate( l3( stotel * quad_cv_sngi * totele ) )
      allocate( l4( stotel * quad_cv_sngi * totele ) )

      normx = 1. ; normy = 0. ; normz = 0.

      ! Now we need to compute QUAD_NLX/Y/Z - get the hex or quad
      ! shape functions quad_n etc.
      ewrite(3,*) 'before shape_l_q_quad lowqua, quad_cv_ngi, mloc:', &
           lowqua, quad_cv_ngi, mloc
      ewrite(3,*) 'dummy_sngi, dummy_snloc, dummy_smloc, quad_cv_ngi:', &
           dummy_sngi, dummy_snloc, dummy_smloc, quad_cv_ngi

      ! Work out local coords of the nodes
      loc_coord_nod_l1 = 0. ; loc_coord_nod_l2 = 0. ; loc_coord_nod_l3 = 0. ; &
           loc_coord_nod_l4 = 0.
      do xnod = 1, x_nonods
         if( d3 ) then 
            loc_coord_nod_l1( xnod ) = &
                 volume_quad_map( 1, x( xnod ), y( xnod ), z( xnod ), lx, ly, lz )
            loc_coord_nod_l2( xnod ) = &
                 volume_quad_map( 2, x( xnod ), y( xnod ), z( xnod ), lx, ly, lz )
            loc_coord_nod_l3( xnod ) = &
                 volume_quad_map( 3, x( xnod ), y( xnod ), z( xnod ), lx, ly, lz )
            loc_coord_nod_l4( xnod ) = &
                 volume_quad_map( 4, x( xnod ), y( xnod ), z( xnod ), lx, ly, lz )
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
         ewrite(3,*)'loc_coord_nod_l1/3x:',loc_coord_nod_l1(xnod), loc_coord_nod_l2(xnod), &
              loc_coord_nod_l3(xnod)
      end do
      do xnod = 1, cv_nloc
         ewrite(3,*)'sn_i_xj:', xnod, ( sn_i_xj( xnod, ele ), ele = 1, x_nonods )
      end do

      Loop_Elements: do ele = 1, totele ! Calculate SDETWEI,RA,SNX,SNY,SNZ for element ELE
         ! What is the fem node belonging to this element (CV_ILOC):
         cv_iloc_belong = 0
         do quad_cv_iloc = 1, quad_cv_nloc
            xnod = x_ndgln( ( ele - 1 ) * quad_cv_nloc + quad_cv_iloc )
            if( fem_nod( xnod ) /= 0 ) cv_iloc_belong = xnod 
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
                     end if
                  end do Loop_Poly2d_1_2
               end do Loop_Poly2d_1_1

            else ! 3D
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
                           x_sndgln( ( sele - 1 ) * quad_cv_snloc + quad_cv_siloc ) = &
                                x_ndgln( ( ele - 1 ) * quad_cv_nloc + quad_cv_iloc )
                           loc_2nd_lev( sele, quad_cv_siloc ) = quad_cv_iloc
                        end if
                     end do Loop_Poly2d_2_3
                  end do Loop_Poly2d_2_2
               end do Loop_Poly2d_2_1
            endif Conditional_Dimension

            ! Now compute the determinant -- Quad_Sdetwei
            do quad_cv_siloc = 1, quad_cv_snloc
               xnod = x_sndgln( ( sele - 1 ) * quad_cv_snloc + quad_cv_siloc )
               xsl( quad_cv_siloc ) = x( xnod )
               ysl( quad_cv_siloc ) = y( xnod )
               if( d3 ) zsl( quad_cv_siloc ) = z( xnod )
               ewrite(3,*)'x/ysl:', ele, sele, quad_cv_siloc, xsl( quad_cv_siloc ), &
                    ysl( quad_cv_siloc )
            end do

            call dgsdetnxloc2( quad_cv_snloc, quad_cv_sngi, &
                 xsl, ysl, zsl, &
                 quad_sn, quad_snlx, quad_snly, quad_scvweight, quad_sdetwei, sarea, &
                 ( ndim == 1 ), ( ndim == 3 ), ( ndim == -2 ), &
                 normxn, normyn, normzn, &
                 normx, normy, normz )
            ewrite(3,*)'normx/y/z:', normx, normy, normz
            ewrite(3,*)'quad_sdetwei:', ( quad_sdetwei( xnod ), xnod = 1, quad_cv_sngi )
            do xnod = 1, quad_cv_sngi
               ewrite(3,*)'normx/y/zn:', normxn(xnod), normyn(xnod), normyn(xnod)
            end do

            ! Take out quadrature points that are inside a CV: 
            ! NB: fem_nod(xnod) = 0 if not a local finite element node ELSE = local node no.
            ! This is determined by looking to see if any of the face nodes 
            ! have a local fem_nod. 
            found_fem_nod = .false.
            do quad_cv_sjloc = 1, quad_cv_snloc
               xnod = x_sndgln( ( sele - 1 ) * quad_cv_snloc + quad_cv_sjloc )
               if( fem_nod( xnod ) /= 0 ) found_fem_nod = .true.
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

            ewrite(3,*)'ele, totele, quad_cv_sngi:', &
                 ele, totele, quad_cv_sngi

            ! Determine the quadrature points and weights
            Loop_NGI: do quad_cv_sgi = 1, quad_cv_sngi 
               cv_sgi = ( ele - 1 ) * quad_cv_sngi * stotel + ( sele - 1 ) * &
                    quad_cv_sngi + quad_cv_sgi
! we do not multiply the weight by the determinant as we are operating in 
! local coordinates...
!               gl_quad_scvweigh( cv_sgi ) = quad_sdetwei( quad_cv_sgi )
               gl_quad_scvweigh( cv_sgi ) = quad_scvweight( quad_cv_sgi )

               xgi = 0. ; ygi = 0. ; zgi = 0.
               do quad_cv_siloc = 1, quad_cv_snloc
                  xnod = x_sndgln( ( sele - 1 ) * quad_cv_snloc + quad_cv_siloc )
                  xgi =  xgi + quad_sn( quad_cv_siloc, quad_cv_sgi ) * x( xnod )
                  ygi =  ygi + quad_sn( quad_cv_siloc, quad_cv_sgi ) * y( xnod )
                  if( d3 ) zgi =  zgi + quad_sn( quad_cv_iloc, quad_cv_sgi ) * z( xnod )
               end do

               if( d3 ) then ! local coords for the tetrahedron
                  gl_quad_l1( cv_sgi ) = volume_quad_map( 1, xgi, ygi, zgi, lx, ly, lz )
                  gl_quad_l2( cv_sgi ) = volume_quad_map( 2, xgi, ygi, zgi, lx, ly, lz )
                  gl_quad_l3( cv_sgi ) = volume_quad_map( 3, xgi, ygi, zgi, lx, ly, lz )
                  gl_quad_l4( cv_sgi ) = volume_quad_map( 4, xgi, ygi, zgi, lx, ly, lz )
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
              ele = stotel * 1, quad_cv_sngi * totele )
         ewrite(3,*)'gl_quad_snly:', xnod, ( gl_quad_snly( xnod, ele ), &
              ele = stotel * 1, quad_cv_sngi * totele )
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
                  ewrite(3,*) cv_sgi, cv_sgj, next_to_cv_iloc_gi( cv_sgi ), &
                       next_to_cv_iloc_gi( cv_sgj ), fem_nod( next_to_cv_iloc_gi( cv_sgj ) )
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

      allocate( sn_2( cv_nloc, stotel * quad_cv_sngi * totele ) ) ; sn_2 = 0.
      allocate( suf_snlx_2( cv_nloc, stotel * quad_cv_sngi * totele ) ) ; suf_snlx_2 = 0.
      allocate( suf_snly_2( cv_nloc, stotel * quad_cv_sngi * totele ) ) ; suf_snly_2 = 0.
      allocate( snlx_2( cv_nloc, stotel * quad_cv_sngi * totele ) ) ; snlx_2 = 0.
      allocate( snly_2( cv_nloc, stotel * quad_cv_sngi * totele ) ) ; snly_2 = 0.
      allocate( snlz_2( cv_nloc, stotel * quad_cv_sngi * totele ) ) ; snlz_2 = 0.
      allocate( scvweigh_2( stotel * quad_cv_sngi * totele ) )
      allocate( cv_neiloc_cells_2( cv_nloc_cells, stotel * quad_cv_sngi * totele ) )
      allocate( cv_neiloc_cells_3( cv_nloc_cells, stotel * quad_cv_sngi * totele ) )
      allocate( l1_2( stotel * quad_cv_sngi * totele ) )
      allocate( l2_2( stotel * quad_cv_sngi * totele ) )
      allocate( l3_2( stotel * quad_cv_sngi * totele ) )
      allocate( l4_2( stotel * quad_cv_sngi * totele ) )

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
         end if
      end do

      ewrite(3,*) 'suf_snlx_2:', suf_snlx_2

      cv_sngi_2 = cv_sgk

      ! Take out repetition of quadrature points
      cv_sgk = 0 ; cv_neiloc_cells = 0 ; l1 = 0. ; l2 = 0. ; l3 = 0.
      Loop_CV_SGI: do cv_sgi = 1, cv_sngi_2
         found = .false.
         Loop_CV_SGJ: do cv_sgj = 1, cv_sgi - 1
            d2_quad = abs( l1_2( cv_sgi ) - l1_2( cv_sgj ) ) &
                 + abs( l2_2( cv_sgi ) - l2_2( cv_sgj ) ) &
                 + abs( l3_2( cv_sgi ) - l3_2( cv_sgj ) ) 
            if( d3 ) d2_quad = d2_quad &
                 + abs( l4_2( cv_sgi ) - l4_2( cv_sgj ) ) 
            if( d2_quad < 1.0e-5 ) found = .true. ! same pt -- remove repetation
         end do Loop_CV_SGJ

         if( .not. found ) then
            cv_sgk = cv_sgk + 1
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
            ewrite(3,*) 'cv_sgk, cv_sgi :' ,cv_sgk, cv_sgi 
            ewrite(3,*) 'suf_snlx_2( :, cv_sgi ):', suf_snlx_2( :, cv_sgi )      
            ewrite(3,*) 'sufnlx( :, cv_sgk ):', sufnlx( :, cv_sgk )       
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

      cvfem_neiloc = 0
      ! calculate cvfem_neiloc from local coords l1-4:
      Loop_SurfaceQuadrature: do cv_sgi = 1, scvngi
         zer_l1 = ( abs( l1( cv_sgi ) ) < 1.0e-4 )
         zer_l2 = ( abs( l2( cv_sgi ) ) < 1.0e-4 )
         zer_l3 = ( abs( l3( cv_sgi ) ) < 1.0e-4 )
         if ( d3 ) zer_l4 = ( abs( l4( cv_sgi )) < 1.0e-4 )
         Conditional_Neiloc: if ( d3 ) then
            if ( zer_l1 .or. zer_l2 .or. zer_l3 .or. zer_l4 ) then
               ! on the surface of the element: 
               do cv_iloc = 1, cv_nloc
                  if ( abs ( sn( cv_iloc, cv_sgi )) > 1.e-4 ) &
                       cvfem_neiloc( cv_iloc, cv_sgi ) = -1
               end do
            endif
         else
            ewrite(3,*)'cv_sgi, zer_l1/2/3:', cv_sgi, (zer_l1.or.zer_l2.or.zer_l3) 
            if ( zer_l1 .or. zer_l2 .or. zer_l3 ) then
               ! on the surface of the element: 
               do cv_iloc = 1, cv_nloc
                  if ( abs( sn( cv_iloc, cv_sgi )) > 1.e-4 ) &
                       cvfem_neiloc( cv_iloc, cv_sgi ) = -1
               end do
            endif
         endif Conditional_Neiloc
      end do Loop_SurfaceQuadrature

      do quad_cv_siloc = 1, cv_nloc
         do cv_sgi = 1, scvngi
            ewrite(3,*)'iloc, gi, cvfem_neiloc::', &
                 quad_cv_siloc, cv_sgi, cvfem_neiloc( quad_cv_siloc, cv_sgi )
         end do
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
         mloc, dummy_smloc, lowqua, npoly )
      implicit none
      ! Compute some local variables for suf_shape_tri_tet
      logical, intent( in ) :: d1, d3
      integer, intent( inout ) :: quad_cv_ngi
      integer, intent( in ) :: quad_cv_nloc
      integer, intent( inout ) :: dummy_sngi, dummy_snloc, nwicel
      integer, intent( in ) :: cv_nloc, cv_sngi, totele 
      integer, intent( inout ) :: quad_u_loc_dummy, mloc, dummy_smloc
      logical, intent( inout ) :: lowqua
      integer, intent( inout ) :: npoly

      ewrite(3,*)'In dummy_tri_tet'

      Conditional_Dimensionality1: if( d3 ) then
         quad_cv_ngi = 8
         dummy_sngi = 4
         dummy_snloc = 4
         nwicel = 1 ! was 3
         if( cv_nloc == 10 ) then ! Quadratic hexs
            quad_cv_ngi = 27
            nwicel = 4
            dummy_sngi = 9
            dummy_snloc = 9
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
            nwicel = 2
         end if
      end if Conditional_Dimensionality1

      !      if ( cv_sngi == 9 ) quad_cv_ngi = 3
      ! Consistency check:
      if( .false. )then
      !if( cv_sngi /= totele * quad_cv_ngi ) then
         FLExit( "Wrong number for CV_NGI" )
      end if

      quad_u_loc_dummy = 1
      mloc = 1
      dummy_smloc = 1
      lowqua = .false.

      npoly = 2 ! Linear default
      if( d1 ) then
         npoly = 1
      elseif ( .not. d3 ) then ! 2D
         if( quad_cv_nloc == 9 ) npoly = 3 ! Quadratic
      else ! 3D
         if( quad_cv_nloc == 27 ) npoly = 3 ! Quadratic
      endif

      if( quad_cv_nloc == 1 ) npoly = 1

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
      real, dimension( quad_cv_nloc ), intent( in ) :: lx, ly, lz ! corner nodes
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: n, nlx, nly, nlz
      real, dimension( cv_ngi ), intent( inout ) ::  cvweigh
      ! Local variables
      logical :: d1, dcyl, d3, lowqua
      integer :: ele, quad_cv_ngi, quad_cv_gi, cv_gi, nwicel, xnod, quad_cv_iloc, &
           quad_u_loc_dummy, mloc, dummy_sngi, dummy_snloc, dummy_smloc
      real :: xgi, ygi, zgi, volume, rsum
      real, dimension( : ), allocatable :: quad_l1, quad_l2, quad_l3, quad_l4, &
           quad_cvweight, detwei, ra, rdummy
      real, dimension( :, : ), allocatable :: quad_n, quad_nlx, quad_nly, quad_nlz, &
           quad_nx, quad_ny, quad_nz

      ewrite(3,*)'in shape_tri_tet'

      d3 = ( ndim == 3 )
      d1 = ( ndim == 1 )
      dcyl = .false. 

      Conditional_Dimensionality1: if( d3 ) then
         quad_cv_ngi = 8
         dummy_sngi = 4
         dummy_snloc = 4
         nwicel = 1 ! was 3
         if( cv_nloc_cells == 10 ) then ! Quadratic hexs
            quad_cv_ngi = 27
            nwicel = 4
            dummy_sngi = 9
            dummy_snloc = 9
         end if
      else
         quad_cv_ngi = 4
         dummy_sngi = 2
         dummy_snloc = 2
         nwicel = 1
         if ( cv_nloc_cells == 6 ) then ! Quadratic quads
            quad_cv_ngi = 9
            dummy_sngi = 3
            dummy_snloc = 3
            nwicel = 2
         end if
      end if Conditional_Dimensionality1

      if ( cv_ngi == 9 ) quad_cv_ngi = 3
      ewrite(3,*) 'totele, cv_nloc_cells, cv_nloc, cv_ngi, quad_cv_ngi:', &
           totele, cv_nloc_cells, cv_nloc, cv_ngi, quad_cv_ngi

      ! Consistency check:
      if( cv_ngi /= totele * quad_cv_ngi ) then
         FLExit( "Wrong number for CV_NGI" )
      end if

      quad_u_loc_dummy = 1
      mloc = 1
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
            ewrite(3,*)'ele, totele, quad_cv_gi, quad_cv_ngi, cv_gi, cv_ngi:', ele, totele, quad_cv_gi, quad_cv_ngi, cv_gi, cv_ngi
            ! the weights need to sum to 0.5 in 2D triangles and 1./6. in 3D tets...
            if(ndim==1) cvweigh( cv_gi ) = detwei( quad_cv_gi )
            if(ndim==2) cvweigh( cv_gi ) = 0.5*detwei( quad_cv_gi )
            if(ndim==3) cvweigh( cv_gi ) = (1./6.)*detwei( quad_cv_gi )
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
      call shatri_hex( quad_l1, quad_l2, quad_l3, quad_l4, rdummy, d3, &
           cv_nloc, cv_ngi, n, nlx, nly, nlz, &
           .true. )
      ewrite(3,*)'cvweigh:',cvweigh
      rsum=0.0
      do cv_gi=1,cv_ngi
         rsum=rsum+cvweigh(cv_gi)
      end do
      ewrite(3,*)'rsum:',rsum
!      stop 2821
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
      integer :: gi, ndim, cv_ele_type_dummy, u_nloc_dummy
      real, dimension( :, : ), allocatable :: cvn_dummy, un_dummy, unlx_dummy, &
           unly_dummy, unlz_dummy
      real, dimension( : ), allocatable :: cvweigh_dummy

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

      cvn = 0.
      Loop_Elements: do ele = 1, totele
         nod = 0
         do quad_cv_iloc = 1, quad_cv_nloc
            xnod = x_ndgln( ( ele - 1 ) * quad_cv_nloc + quad_cv_iloc )
            if( fem_nod( xnod ) /= 0 ) nod = fem_nod( xnod )
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
         else ! Lagrange 27 nodes 3D element - bilinear pressure
            call re3d27( lowqua, ngi, 0, nloc, mloc, &
                 m, weight, n, nlx, nly, nlz, &
                 rdum, rdum, rdum )
         end if

      end Select

      deallocate( rdum )

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
      real, dimension( : ) :: lx, ly, lz
      ! Local variables
      integer :: ipt
      real, dimension( :, : ), allocatable :: lxyz
      real :: xyzgi(3)
      real :: loc_x_coord, loc_y_coord, loc_z_coord, loc_zz_coord

      allocate( lxyz( n + 1, n ) )

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
      REAL, DIMENSION( SNGI ), intent( inout ) :: NORMXN, NORMYN ,NORMZN
      REAL, intent( inout ) :: NORMX, NORMY, NORMZ
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
            NORMZ=0.
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


  end module shape_functions_NDim

