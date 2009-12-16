!     Copyright (C) 2006 Imperial College London and others.
!     
!     Please see the AUTHORS file in the main source directory for a full list
!     of copyright holders.
!     
!     Prof. C Pain
!     Applied Modelling and Computation Group
!     Department of Earth Science and Engineering
!     Imperial College London
!     
!     C.Pain@Imperial.ac.uk
!     
!     This library is free software; you can redistribute it and/or
!     modify it under the terms of the GNU Lesser General Public
!     License as published by the Free Software Foundation,
!     version 2.1 of the License.
!     
!     This library is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!     Lesser General Public License for more details.
!     
!     You should have received a copy of the GNU Lesser General Public
!     License along with this library; if not, write to the Free Software
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!     USA

#include "fdebug.h"

module sdetnx_module

  use fldebug
  use shape_transformations

  implicit none

  private

  public :: sdetnx2

contains

SUBROUTINE SDETNX(ELE,LOCNOD, &
             XONDGL, TOTELE,XNONOD,NLOC,SNLOC,SNGI,  &
             X,Y,Z, &
             SN,SNLX,SNLY, SWEIGH, DETWEI,SAREA, D3,DCYL,  &
             NORMXN,NORMYN,NORMZN, &
             NORMX,NORMY,NORMZ) 
       INTEGER TOTELE,XNONOD,NLOC,SNLOC,SNGI
       INTEGER ELE,LOCNOD(SNLOC),XONDGL(TOTELE*NLOC)
          REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
          REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI) 
          REAL SWEIGH(SNGI) 
       REAL DETWEI(SNGI)
       REAL SAREA
       LOGICAL D3,DCYL 
       REAL NORMXN(SNGI),NORMYN(SNGI),NORMZN(SNGI)
       REAL NORMX,NORMY,NORMZ
       ! Local variables...
       REAL PIE
       PARAMETER(PIE=3.141592654) 
       INTEGER GI,SL,L,IGLX
       REAL DXDLX,DXDLY,DYDLX,DYDLY,DZDLX,DZDLY
       REAL A,B,C,DETJ,RGI,TWOPIE
!
         SAREA=0.
!

         IF(D3) THEN
       DO GI=1,SNGI
!
       DXDLX=0.
       DXDLY=0.
       DYDLX=0.
       DYDLY=0.
       DZDLX=0.
       DZDLY=0.
!
       DO SL=1,SNLOC
           L=LOCNOD(SL) 
           IGLX=XONDGL((ELE-1)*NLOC+L)
         DXDLX=DXDLX + SNLX(SL,GI)*X(IGLX)
         DXDLY=DXDLY + SNLY(SL,GI)*X(IGLX) 
         DYDLX=DYDLX + SNLX(SL,GI)*Y(IGLX) 
         DYDLY=DYDLY + SNLY(SL,GI)*Y(IGLX) 
         DZDLX=DZDLX + SNLX(SL,GI)*Z(IGLX) 
         DZDLY=DZDLY + SNLY(SL,GI)*Z(IGLX) 
        ENDDO
          A = DYDLX*DZDLY - DYDLY*DZDLX
          B = DXDLX*DZDLY - DXDLY*DZDLX
          C = DXDLX*DYDLY - DXDLY*DYDLX
!
         DETJ=SQRT( A**2 + B**2 + C**2)
         DETWEI(GI)=DETJ*SWEIGH(GI)
           SAREA=SAREA+DETWEI(GI)
!
! Calculate the normal at the Gauss pts...
!  TANX1=DXDLX,TANY1=DYDLX,TANZ1=DZDLX,    TANX2=DXDLY,TANY2=DYDLY,TANZ2=DZDLY
! Perform x-product. N=T1 x T2
           CALL NORMGI(NORMXN(GI),NORMYN(GI),NORMZN(GI),  &
          DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY, &
          NORMX,NORMY,NORMZ) 
!
        ENDDO
! IF(D3) THEN...
         ELSE
! ENDOF IF(D3) THEN ELSE...
         TWOPIE=1.0 
         IF(DCYL) TWOPIE=2.*PIE
!
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
             L=LOCNOD(SL) 
           IGLX=XONDGL((ELE-1)*NLOC+L) 
           DXDLX=DXDLX + SNLX(SL,GI)*X(IGLX) 
           DYDLX=DYDLX + SNLX(SL,GI)*Y(IGLX) 
           RGI=RGI+SN(SL,GI)*Y(IGLX)
!             ewrite(3,*) 'IGLX,X(IGLX),Y(IGLX):',IGLX,X(IGLX),Y(IGLX)
        ENDDO
         IF(.NOT.DCYL) RGI=1.0 
         DETJ=SQRT( DXDLX**2 + DYDLX**2 )
         DETWEI(GI)=TWOPIE*RGI*DETJ*SWEIGH(GI)
           SAREA=SAREA+DETWEI(GI) 
!           ewrite(3,*) 'DETJ,DETWEI(GI),SWEIGH(GI),SAREA:',
!     &              DETJ,DETWEI(GI),SWEIGH(GI),SAREA
!           ewrite(3,*) 'DETWEI(GI),TWOPIE,RGI,DETJ,SWEIGH(GI):',
!     &              DETWEI(GI),TWOPIE,RGI,DETJ,SWEIGH(GI)
!
! Calculate the normal at the Gauss pts...
!  TANX1=DXDLX,TANY1=DYDLX,TANZ1=DZDLX,    TANX2=DXDLY,TANY2=DYDLY,TANZ2=DZDLY
! Perform x-product. N=T1 x T2
           NORMZ=0.
           CALL NORMGI(NORMXN(GI),NORMYN(GI),NORMZN(GI),  &
          DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY, &
          NORMX,NORMY,NORMZ) 
!
         ENDDO
         ENDIF
         RETURN
         END subroutine SDETNX
!
!
!
!
       SUBROUTINE SDETNX2(ELE,SNDGLN, &
             STOTEL,XNONOD,SNLOC,SNGI, & 
             X,Y,Z, &
             SN,SNLX,SNLY, SWEIGH, DETWEI,SAREA, D3,DCYL,  &
             NORMXN,NORMYN,NORMZN, &
             NORMX,NORMY,NORMZ) 
       INTEGER STOTEL,XNONOD,SNLOC,SNGI
       INTEGER ELE,SNDGLN(STOTEL*SNLOC)
          REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
          REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI) 
          REAL SWEIGH(SNGI) 
       REAL DETWEI(SNGI)
       REAL SAREA
       LOGICAL D3,DCYL 
       REAL NORMXN(SNGI),NORMYN(SNGI),NORMZN(SNGI)
       REAL NORMX,NORMY,NORMZ
! Local variables...
       REAL PIE
       PARAMETER(PIE=3.141592654) 
       INTEGER GI,SL,IGLX
       REAL DXDLX,DXDLY,DYDLX,DYDLY,DZDLX,DZDLY
       REAL A,B,C,DETJ,RGI,TWOPIE
!
         SAREA=0.
!

         IF(D3) THEN
       DO GI=1,SNGI
!
       DXDLX=0.
       DXDLY=0.
       DYDLX=0.
       DYDLY=0.
       DZDLX=0.
       DZDLY=0.
!
       DO SL=1,SNLOC
           IGLX=SNDGLN((ELE-1)*SNLOC+SL)
         DXDLX=DXDLX + SNLX(SL,GI)*X(IGLX)
         DXDLY=DXDLY + SNLY(SL,GI)*X(IGLX) 
         DYDLX=DYDLX + SNLX(SL,GI)*Y(IGLX) 
         DYDLY=DYDLY + SNLY(SL,GI)*Y(IGLX) 
         DZDLX=DZDLX + SNLX(SL,GI)*Z(IGLX) 
         DZDLY=DZDLY + SNLY(SL,GI)*Z(IGLX) 
       ENDDO
          A = DYDLX*DZDLY - DYDLY*DZDLX
          B = DXDLX*DZDLY - DXDLY*DZDLX
          C = DXDLX*DYDLY - DXDLY*DYDLX
!
         DETJ=SQRT( A**2 + B**2 + C**2)
         DETWEI(GI)=DETJ*SWEIGH(GI)
           SAREA=SAREA+DETWEI(GI)
!
! Calculate the normal at the Gauss pts...
!  TANX1=DXDLX,TANY1=DYDLX,TANZ1=DZDLX,    TANX2=DXDLY,TANY2=DYDLY,TANZ2=DZDLY
! Perform x-product. N=T1 x T2
         CALL NORMGI(NORMXN(GI),NORMYN(GI),NORMZN(GI), &
          DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY, &
          NORMX,NORMY,NORMZ)  
!
        ENDDO
! IF(D3) THEN...
         ELSE
! ENDOF IF(D3) THEN ELSE...
         TWOPIE=1.0 
         IF(DCYL) TWOPIE=2.*PIE
!
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
           IGLX=SNDGLN((ELE-1)*SNLOC+SL)
           DXDLX=DXDLX + SNLX(SL,GI)*X(IGLX) 
           DYDLX=DYDLX + SNLX(SL,GI)*Y(IGLX) 
           RGI=RGI+SN(SL,GI)*Y(IGLX)
!             ewrite(3,*) 'IGLX,X(IGLX),Y(IGLX):',IGLX,X(IGLX),Y(IGLX)
         ENDDO
         IF(.NOT.DCYL) RGI=1.0 
         DETJ=SQRT( DXDLX**2 + DYDLX**2 )
         DETWEI(GI)=TWOPIE*RGI*DETJ*SWEIGH(GI)
           SAREA=SAREA+DETWEI(GI) 
!           ewrite(3,*) 'DETJ,DETWEI(GI),SWEIGH(GI),SAREA:',
!     &              DETJ,DETWEI(GI),SWEIGH(GI),SAREA
!           ewrite(3,*) 'DETWEI(GI),TWOPIE,RGI,DETJ,SWEIGH(GI):',
!     &              DETWEI(GI),TWOPIE,RGI,DETJ,SWEIGH(GI)
!
! Calculate the normal at the Gauss pts...
!  TANX1=DXDLX,TANY1=DYDLX,TANZ1=DZDLX,    TANX2=DXDLY,TANY2=DYDLY,TANZ2=DZDLY
! Perform x-product. N=T1 x T2
           NORMZ=0.
       CALL NORMGI(NORMXN(GI),NORMYN(GI),NORMZN(GI), &
          DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY, &
          NORMX,NORMY,NORMZ) 
!
         ENDDO
         ENDIF
         RETURN
         END subroutine SDETNX2

end module sdetnx_module
