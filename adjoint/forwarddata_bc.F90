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
!    C.Pain@Imperial.ac.uk
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
SUBROUTINE FORWARDDATA_BC(BCU1,BCU2,BCV1,BCV2,BCW1,BCW2,   &
     NOBCU,NOBCV,NOBCW,D3,ITSTIME,MAXREC,MRECORD,          &
     NOCVA,CVA,IMEM_CVA,                                   &
     CVA_X,CVA_Y,CVA_Z,                                    &
     X,Y,Z,NONODS,NPHASE)
   use FLDebug
  IMPLICIT NONE
  INTEGER ::NOBCU,NOBCV,NOBCW
  INTEGER ::NONODS,NPHASE
  INTEGER ::BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
  REAL    ::BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
  REAL	  ::X(NONODS*NPHASE),Y(NONODS*NPHASE),Z(NONODS*NPHASE)
  LOGICAL ::D3
  INTEGER ::NOCVA
  REAL	  ::CVA(NOCVA),CVA_X(NOCVA),CVA_Y(NOCVA),CVA_Z(NOCVA)
  INTEGER ::IMEM_CVA(NOCVA)
  INTEGER ::ITSTIME,MAXREC
  INTEGER ::MRECORD(MAXREC)

  INTEGER ::K1,K2,II

  ewrite(3,*) 'NOBCU,NOBCV,NOBCW',NOBCU,NOBCV,NOBCW
  K1 = MRECORD(ITSTIME)
  K2 = MRECORD(ITSTIME)+NOBCU-1
  CVA(K1:K2) = BCU1
  IMEM_CVA(K1:K2) = BCU2
!  DO II =1,NOBCU
!     CVA_X(K1+II-1) = X(BCU2(II))
!     CVA_Y(K1+II-1) = Y(BCU2(II))
!     CVA_Z(K1+II-1) = Z(BCU2(II))
!     if ( (bcu1(ii)-0.0).lt.1.0E-6) then
!        ewrite(3,*) 'here bcu=0',X(BCU2(II)),Y(BCU2(II)),Z(BCU2(II)),BCU1(II)
!     endif
!  ENDDO

  !	goto 10
  K1 = K2+1
  K2 = K2+NOBCV
  CVA(K1:K2) = BCV1
  IMEM_CVA(K1:K2) = BCV2
!  DO II =1,NOBCV
!     CVA_X(K1+II-1) = X(BCV2(II))
!     CVA_Y(K1+II-1) = Y(BCV2(II))
!     CVA_Z(K1+II-1) = Z(BCV2(II))
!  ENDDO

  IF(D3) THEN
     !	  ewrite(3,*) 'BCW',NOBCW,BCW1
     K1 = K2+1
     K2 = K2+NOBCW
     CVA(K1:K2) = BCW1
     IMEM_CVA(K1:K2) = BCW2
     do ii = 1,nobcw
        if(bcw2(ii).eq.0) then 
          ewrite(3,*) 'bcw2,ii',imem_cva(ii),ii
          stop 1234
        endif
     enddo
 !    DO II =1,NOBCW
 !    CVA_X(K1+II-1) = X(BCW2(II))
 !    CVA_Y(K1+II-1) = Y(BCW2(II))
 !    CVA_Z(K1+II-1) = Z(BCW2(II))
 !    ENDDO
  ENDIF

  !10	continue	
!  MRECORD(ITSTIME+1) = K2+1
  ewrite(3,*) 'MRECORD in forward_bc',MRECORD(ITSTIME+1)
END SUBROUTINE FORWARDDATA_BC


!==================================


SUBROUTINE FORWARDDATA_BCH(BCU1,BCU2,BCV1,BCV2,BCW1,BCW2,   &
     NOBCU,NOBCV,NOBCW,D3,ITSTIME,MAXREC,MRECORD,          &
     NOCVA,CVA,IMEM_CVA,                                   &
     CVA_X,CVA_Y,CVA_Z,                                    &
     X,Y,Z,NONODS,NPHASE,                                  &
     NOBCH,BCH2,BCH1,HYDSAP,NHYDSP,GOTFRE)
   use FLDebug
  IMPLICIT NONE
  INTEGER ::NOBCU,NOBCV,NOBCW
  INTEGER ::NONODS,NPHASE
  INTEGER ::BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
  REAL    ::BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
  REAL	  ::X(NONODS*NPHASE),Y(NONODS*NPHASE),Z(NONODS*NPHASE)
  LOGICAL ::D3
  INTEGER ::NOCVA
  REAL	  ::CVA(NOCVA),CVA_X(NOCVA),CVA_Y(NOCVA),CVA_Z(NOCVA)
  INTEGER ::IMEM_CVA(NOCVA)
  INTEGER ::ITSTIME,MAXREC
  INTEGER ::MRECORD(MAXREC)

  ! free surface
  INTEGER  ::NOBCH
  REAL     ::BCH1(NOBCH)
  INTEGER  ::BCH2(NOBCH)
  INTEGER  ::NHYDSP
  REAL     ::HYDSAP(NHYDSP)
  LOGICAL  ::GOTFRE

  INTEGER ::K1,K2,II,NOD

  ewrite(3,*) 'NOBCU,NOBCV,NOBCW,NOBCH',NOBCU,NOBCV,NOBCW,NOBCH

  IF(GOTFRE) THEN

    K1 = MRECORD(ITSTIME)
    K2 = MRECORD(ITSTIME)+NOBCH-1

     CVA(K1:K2) = BCH1
     IMEM_CVA(K1:K2) = BCH2
!      DO  II = K1,K2
!        CVA(II) = BCH1(II-K1+1)
!        IMEM_CVA(II) = BCH2(II-K1+1)
!      ENDDO

  ELSE
    K1 = MRECORD(ITSTIME)
    K2 = MRECORD(ITSTIME)+NOBCU-1
    CVA(K1:K2) = BCU1
    IMEM_CVA(K1:K2) = BCU2
!    DO II =1,NOBCU
!     CVA_X(K1+II-1) = X(BCU2(II))
!     CVA_Y(K1+II-1) = Y(BCU2(II))
!     CVA_Z(K1+II-1) = Z(BCU2(II))
!     if ( (bcu1(ii)-0.0).lt.1.0E-6) then
!        ewrite(3,*) 'here bcu=0',X(BCU2(II)),Y(BCU2(II)),Z(BCU2(II)),BCU1(II)
!     endif
!    ENDDO

  !	goto 10
    K1 = K2+1
    K2 = K2+NOBCV
    CVA(K1:K2) = BCV1
    IMEM_CVA(K1:K2) = BCV2
!    DO II =1,NOBCV
!     CVA_X(K1+II-1) = X(BCV2(II))
!     CVA_Y(K1+II-1) = Y(BCV2(II))
!     CVA_Z(K1+II-1) = Z(BCV2(II))
!    ENDDO

    IF(D3) THEN
     !	  ewrite(3,*) 'BCW',NOBCW,BCW1
       K1 = K2+1
       K2 = K2+NOBCW
       CVA(K1:K2) = BCW1
       IMEM_CVA(K1:K2) = BCW2
       do ii = 1,nobcw
        if(bcw2(ii).eq.0) then 
          ewrite(3,*) 'bcw2,ii',imem_cva(ii),ii
          stop 1234
        endif
       enddo
 !     DO II =1,NOBCW
 !     CVA_X(K1+II-1) = X(BCW2(II))
 !     CVA_Y(K1+II-1) = Y(BCW2(II))
 !     CVA_Z(K1+II-1) = Z(BCW2(II))
 !     ENDDO
    ENDIF

  ENDIF
  !10	continue	
!  MRECORD(ITSTIME+1) = K2+1
  ewrite(3,*)  'k1,k2',k1,k2
  ewrite(3,*) 'MRECORD in forward_bc',MRECORD(ITSTIME)
END SUBROUTINE FORWARDDATA_BCH

