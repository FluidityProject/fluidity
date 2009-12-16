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

#include "confdefs.h"
#include "fdebug.h"
module mesh_connections
use FLDebug
use Element_Numbering, only: ilink2, silink2
use shape_transformations
implicit none
contains

    SUBROUTINE GETFINELE4(TOTELE,NLOC,SNLOC,NDGLNO, COLELE,NONODS)
   ! This sub caluculates COLELE the element
   ! connectivitiy list in order of faces.
   INTEGER TOTELE,NLOC,SNLOC,NONODS
   INTEGER NDGLNO(TOTELE*NLOC)
   INTEGER COLELE(TOTELE*5)
! LOCAL VARIABLES...
       INTEGER NFACE
       PARAMETER(NFACE=4)
       INTEGER LOCLIST(NFACE,3),CFACE(5),EFACE(5)
   INTEGER ELE,ILOC,JLOC,ILOC2,NOD,INOD,JNOD,COUNT,ELE2,I,HIT
   INTEGER IFACE,IFACE2,KFACE
   INTEGER NCOLEL
   LOGICAL FOUND

   INTEGER, ALLOCATABLE, DIMENSION(:)::FINELE
   INTEGER, ALLOCATABLE, DIMENSION(:)::FINTRAN
   INTEGER, ALLOCATABLE, DIMENSION(:)::COLTRAN
   INTEGER, ALLOCATABLE, DIMENSION(:)::ICOUNT
   ALLOCATE(FINELE(TOTELE+1))
   ALLOCATE(FINTRAN(NONODS+1))
   ALLOCATE(COLTRAN(TOTELE*5))
   ALLOCATE(ICOUNT(MAX(NONODS,TOTELE)))

!   print*, 'now in getfinele4'
!   print*, 'totele, nloc,snloc, nonods',totele, nloc,snloc, nonods
!   print*, 'ndglno(1:totele*nloc)', ndglno(1), ndglno(totele*nloc)
!   print*, 'max-min ndglno', maxval(ndglno), minval(ndglno)
!   print *,'just inside getfinele'
!   pause 


   ICOUNT(1:NONODS) = 0
   DO ELE=1,TOTELE
     DO ILOC=1,NLOC
        NOD=NDGLNO((ELE-1)*NLOC+ILOC)
        ICOUNT(NOD)=ICOUNT(NOD)+1
      END DO
   END DO

   FINTRAN(1)=1
   DO NOD=1,NONODS
     FINTRAN(NOD+1)=FINTRAN(NOD)+ICOUNT(NOD)
   END DO

   ICOUNT(1:NONODS) = 0
   DO ELE=1,TOTELE
     DO ILOC=1,NLOC
        NOD=NDGLNO((ELE-1)*NLOC+ILOC)
        COLTRAN(FINTRAN(NOD)+ICOUNT(NOD))=ELE
        ICOUNT(NOD)=ICOUNT(NOD)+1
      END DO
   END DO

   ICOUNT(1:TOTELE) = 0
   COLELE(1:TOTELE*5) = 0
   DO ELE=1,TOTELE
     DO ILOC=1,NLOC
        NOD=NDGLNO((ELE-1)*NLOC+ILOC)
        DO COUNT=FINTRAN(NOD),FINTRAN(NOD+1)-1
          ELE2=COLTRAN(COUNT)
! Put ELE2 into list FINELE,COLELE
          FOUND=.FALSE.
          DO I=1,ICOUNT(ELE)
            IF(COLELE((ELE-1)*5+I).EQ.ELE2) FOUND=.TRUE.
          END DO
          IF(.NOT.FOUND) THEN
! Do elements ELE and ELE2 share at least 3nodes?
            HIT=0
            DO ILOC2=1,NLOC
              INOD=NDGLNO((ELE-1)*NLOC+ILOC2)
              DO JLOC=1,NLOC
                JNOD=NDGLNO((ELE2-1)*NLOC+JLOC)
                IF(INOD.EQ.JNOD) HIT=HIT+1
              END DO
            END DO
            IF(HIT.GE.SNLOC) THEN
               ICOUNT(ELE)=ICOUNT(ELE)+1
               COLELE((ELE-1)*5+ICOUNT(ELE))=ELE2
            ENDIF
          ENDIF
        END DO
      END DO
   END DO



       IFACE=1
       LOCLIST(IFACE,1)=1
       LOCLIST(IFACE,2)=2
       LOCLIST(IFACE,3)=3
       IFACE=2
       LOCLIST(IFACE,1)=1
       LOCLIST(IFACE,2)=2
       LOCLIST(IFACE,3)=4
       IFACE=3
       LOCLIST(IFACE,1)=3
       LOCLIST(IFACE,2)=2
       LOCLIST(IFACE,3)=4
       IFACE=4
       LOCLIST(IFACE,1)=1
       LOCLIST(IFACE,2)=3
       LOCLIST(IFACE,3)=4

!   ele=3
!     DO I=1,5
!! Which face is it?
!       ELE2=COLELE((ELE-1)*5+I)
!       EWRITE(3,*) 'I=',I
!       do iloc=1,nloc
!         INOD=NDGLNO((ELE-1)*NLOC+ILOC)
!              DO JLOC=1,NLOC
!                JNOD=NDGLNO((ELE2-1)*NLOC+JLOC)
!                if(inod.eq.jnod) then
!                   ewrite(3,*) 'ele2,jnod:',ele2,jnod
!                endif
!              END DO
!       END DO
!      END DO
!       IFACE=3
!      ewrite(3,*) 'NDGLNO((ELE-1)*NLOC+LOCLIST(IFACE,1)):',NDGLNO((ELE-1)*NLOC+LOCLIST(IFACE,1))
!      ewrite(3,*) 'NDGLNO((ELE-1)*NLOC+LOCLIST(IFACE,2)):',NDGLNO((ELE-1)*NLOC+LOCLIST(IFACE,2))
!      ewrite(3,*) 'NDGLNO((ELE-1)*NLOC+LOCLIST(IFACE,3)):',NDGLNO((ELE-1)*NLOC+LOCLIST(IFACE,3))
!
!     ewrite(3,*) 'before reordering'

!     ewrite(3,*) 'COLELE:',(colele(i),i=1,15)
! Put the central element at the end of list COLELE((ELE-1)*5+5)
   DO ELE=1,TOTELE
     DO IFACE=1,NFACE
! Which face is it?
       ELE2=COLELE((ELE-1)*5+IFACE)
       IF(ELE2.EQ.ELE) THEN
! swop over
         COLELE((ELE-1)*5+IFACE)=COLELE((ELE-1)*5+5)
         COLELE((ELE-1)*5+5)=ELE
       ENDIF
     END DO
   END DO


! Now put COLELE in order of faces...
   DO ELE=1,TOTELE
     DO IFACE2=1,NFACE
! Which face is it?
       ELE2=COLELE((ELE-1)*5+IFACE2)
!       if(ele.eq.3) then
!         ewrite(3,*) 'ele2:',ele2
!       endif
       IF(ELE2.NE.0) THEN
          KFACE=0
          DO IFACE=1,NFACE
            HIT=0
            DO ILOC=1,NLOC
              DO I=1,3
                 IF(NDGLNO((ELE-1)*NLOC+LOCLIST(IFACE,I)) &
                 .EQ.NDGLNO((ELE2-1)*NLOC+ILOC)) HIT=HIT+1
              END DO
            END DO
            IF(HIT.EQ.3) KFACE=IFACE
          END DO
          CFACE(IFACE2)=KFACE
          EFACE(IFACE2)=ELE2
        ELSE
          CFACE(IFACE2)=0
          EFACE(IFACE2)=0
        ENDIF
      END DO

! Put in order of faces...(leave the last 5th entry)
     DO IFACE2=1,NFACE
       COLELE((ELE-1)*5+IFACE2)=0
       DO IFACE=1,NFACE
         IF(CFACE(IFACE).EQ.IFACE2) COLELE((ELE-1)*5+IFACE2)=EFACE(IFACE)
       END DO
!       if(ele.eq.3) then
!          ewrite(3,*) 'colele((ELE-1)*5+IFACE2):',colele((ELE-1)*5+IFACE2)
!       endif
     END DO

   END DO

!      STOP 39321
!     ewrite(3,*) 'COLELE:',(colele(i),i=1,15)
!     stop 393

!    print*, 'at the end of getfinele4'
!    pause
   RETURN
 END SUBROUTINE GETFINELE4

 SUBROUTINE MLOCELE2ELEFACE(MLOCOTHILOC,MLOCSILOC2ILOC,SNLOC,NLOC, &
                            SMLOC,MLOC,OTHILOC,SILOC2ILOC)
   ! This subroutine forms MLOCOTHILOC and MLOCSILOC2ILOC, the pointer to the
   ! local node number on the other side of the face and the
   ! pointer of the local node number given a surface node number.
   INTEGER SNLOC,NLOC,SMLOC,MLOC
   INTEGER MLOCOTHILOC(MLOC),MLOCSILOC2ILOC(SMLOC)
   INTEGER OTHILOC(NLOC),SILOC2ILOC(SNLOC)
   ! Local variables...
   INTEGER ILOC,ILOC2,JLOC,JLOC2,SKLOC,KLOC,KLOC2
! iNTEGER FUNCTION...
   INTEGER SILOC,SJLOC
   !

! Calculate MLOCSILOC2ILOC
   DO SILOC=1,SNLOC
     ILOC=SILOC2ILOC(SILOC)
     ILOC2=OTHILOC(ILOC)
     DO SJLOC=1,SNLOC
       JLOC=SILOC2ILOC(SJLOC)
       JLOC2=OTHILOC(JLOC)
       SKLOC=SILINK2(SILOC,SJLOC)
       KLOC=ILINK2(ILOC,JLOC)
       KLOC2=ILINK2(ILOC2,JLOC2)
       MLOCSILOC2ILOC(SKLOC)=KLOC
       MLOCOTHILOC(KLOC)=KLOC2
     END DO
   END DO
   RETURN
 END SUBROUTINE MLOCELE2ELEFACE

SUBROUTINE ELE2ELEFACE(OTHILOC,SILOC2ILOC,SNLOC,NLOC,TOTELE, &
      XONDGL,NDGLNO,ELE,ELE2)
   ! This subroutine forms OTHILOC from SILOC2ILOC, the pointer to the
   ! local node number on the other side of the face and the
   ! pointer of the local node number given a surface node number.
   INTEGER SNLOC,NLOC,TOTELE,ELE,ELE2
   INTEGER OTHILOC(NLOC),SILOC2ILOC(SNLOC)
   INTEGER XONDGL(TOTELE*NLOC),NDGLNO(TOTELE*NLOC)
   ! Local variables...
   INTEGER ICOUNT,SILOC,ILOC,ILOC2,XNODI,XNODI2
   INTEGER NODI,NODI2
   !
   OTHILOC(1:NLOC) = 0
   DO SILOC=1,SNLOC
      ILOC=SILOC2ILOC(SILOC)
      XNODI=XONDGL((ELE-1)*NLOC+ILOC)
      NODI =NDGLNO((ELE-1)*NLOC+ILOC)
      DO ILOC2=1,NLOC
         XNODI2=XONDGL((ELE2-1)*NLOC+ILOC2)
         NODI2 =NDGLNO((ELE2-1)*NLOC+ILOC2)
         IF(NODI.EQ.NODI2) THEN
            OTHILOC(ILOC)=ILOC2
         ENDIF
      END DO
   END DO
!   IF(ELE.EQ.ELE2) THEN
!     ewrite(3,*)' ele2eleface wrongly called ele,ele2:',ele,ele2
!     stop 3282
!   ENDIF
!   IF(ELE2.eq.0) THEN
!     ewrite(3,*)' ele2eleface wrongly called ele2:',ele2
!     stop 3283
!   ENDIF
 END SUBROUTINE ELE2ELEFACE

end module mesh_connections
