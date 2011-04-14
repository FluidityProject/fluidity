C Copyright (C) 2006 Imperial College London and others.
C 
C Please see the AUTHORS file in the main source directory for a full list
C of copyright holders.
C 
C Adrian Umpleby
C Applied Modelling and Computation Group
C Department of Earth Science and Engineering
C Imperial College London
C 
C adrian@Imperial.ac.uk
C 
C This library is free software; you can redistribute it and/or
C modify it under the terms of the GNU Lesser General Public
C License as published by the Free Software Foundation; either
C version 2.1 of the License.
C 
C This library is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C Lesser General Public License for more details.
C 
C You should have received a copy of the GNU Lesser General Public
C License along with this library; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
C USA
      SUBROUTINE EDGELS( BIGLST, NODLST, IPTR, IELM,
     :                   ELPTRS, EDPTRS, NDPTRS, CURADD,
     :                   SURFAC, INTRNL )
C-----------------------------------------------------------------------
C
C - This subroutine finds a list of elements which surround an edge,
C - returning the list in ELPTRS and the number of them in CURADD.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
C
      INCLUDE 'blktim.i'
C
      INTEGER ELPTRS(3,MAXEDE), EDPTRS(3,MAXEDE+1), NDPTRS(MAXEDE+1),
     :        IPTR, IELM, CURADD, SURFAC, INTRNL
C
      INTEGER EDGLST(6), CONELM(4), NXT, LST, INEWEL, IFLAG, INEWED,
     :        CURCHK, I, CURELM, ED1, ED2, ED3, INOD, NODS(4), I1, I2,
     :        ELM1, ELM2, ENDE, IREG, PRVREG, FSTREG, INFO(6), LSTELE,
     :        LSTUED, LSTLED, CUREDG, VALINF, UPPELE, LOWELE, OPSTED,
     :        CHKEL1, CHKEL2, UPPED1, UPPED2, LOWED1, LOWED2, SWCHEL,
     :        SWUPED, SWLWED, UPPEND, LOWEND, j
C
      LOGICAL GETFLG, ISGEOM, ISINSD, goterr
C
      CALL ADAPT_SECOND( TIMSTT )
C
      goterr = .false.
c      if( iptr .eq. 32454 ) goterr = .true.
c
      if( goterr ) then
        print*,' '
        print*,'Into EDGELS... ',ierr
        print*,'   IPTR,IELM: ',IPTR,IELM
        print*,'   NODES: ',biglst(1,iptr),biglst(2,iptr)
      end if
C
      CURADD = 0
      CURCHK = 1
      SURFAC = 0
      INTRNL = 0
C
      IF( IERR .NE. 0 ) RETURN
C
      INFO(1) = 6861915
      INFO(2) = 4953702
      INFO(3) = 6103629
      INFO(4) = 1878962
      INFO(5) = 2735377
      INFO(6) = 3872931
C
      ISGEOM = GETFLG( BIGLST, NODLST, IPTR, 2 )
      ISINSD = GETFLG( BIGLST, NODLST, IPTR, 3 ) .AND. .NOT. ISGEOM
C
      CURELM = IELM
      LSTELE = 0
      LSTUED = 0
      LSTLED = 0
C - just to get rid of compile warnings...
      UPPEND = 0
      SWLWED = 0
      SWUPED = 0
      PRVREG = 0
      LOWEND = 0
      FSTREG = 0
C
  5   CONTINUE
C
      CALL ELMINF( BIGLST, NODLST, CURELM, CONELM, EDGLST,
     :             INEWEL, IFLAG, IREG, NXT, LST )
C
      if( ierr .ne. 0 ) then
         print*,'*** EDGELS: Got error from ELMINF (1)'
         return
      end if
C
      CUREDG = 1
 10   IF( EDGLST(CUREDG) .NE. IPTR ) THEN
         CUREDG = CUREDG + 1
         IF( CUREDG .LT. 7 ) GOTO 10
         PRINT*,'*** EDGELS: FAILED TO FIND EDGE IN ELEMENT!'
         print*,curelm,iptr
         print*,conelm
         print*,edglst
         ierr = -3
         return
      END IF
C
      VALINF = INFO(CUREDG)
      UPPELE = MOD(VALINF,4) + 1
      VALINF = VALINF / 4
      LOWELE = MOD(VALINF,4) + 1
      VALINF = VALINF / 4
      OPSTED = MOD(VALINF,8)
      VALINF = VALINF / 8
      CHKEL1 = MOD(VALINF,4) + 1
      VALINF = VALINF / 4
      UPPED1 = MOD(VALINF,8)
      VALINF = VALINF / 8
      LOWED1 = MOD(VALINF,8)
      VALINF = VALINF / 8
      CHKEL2 = MOD(VALINF,4) + 1
      VALINF = VALINF / 4
      UPPED2 = MOD(VALINF,8)
      VALINF = VALINF / 8
      LOWED2 = MOD(VALINF,8)
      if( valinf .ne. lowed2 ) then
         print*,'*** EDGELS: INFO FAILED TO BE CONSISTENT?!?'
         print*,'*** EDGELS: BAD EDGE NUMBER FROM INFO'
         print*,curelm,iptr,ielm,valinf
         print*,conelm
         print*,edglst
         print*,curedg,opsted,info(curedg)
         print*,uppele,lowele,chkel1,chkel2
         print*,upped1,lowed1,upped2,lowed2
         ierr = -4
         return
      end if
C
      if( opsted .gt. 6 .or. upped1 .gt. 6 .or. lowed1 .gt. 6 .or.
     :    upped2 .gt. 6 .or. lowed2 .gt. 6 ) then
         print*,'*** EDGELS: BAD EDGE NUMBER FROM INFO'
         print*,curelm,iptr,ielm
         print*,conelm
         print*,edglst
         print*,curedg,opsted,info(curedg)
         print*,uppele,lowele,chkel1,chkel2
         print*,upped1,lowed1,upped2,lowed2
         ierr = -4
         return
      end if
C
      CHKEL1 = CONELM(CHKEL1)
      CHKEL2 = CONELM(CHKEL2)
C
      UPPELE = CONELM(UPPELE)
      LOWELE = CONELM(LOWELE)
c
      conelm(1) = chkel1
      conelm(2) = chkel2
      conelm(3) = uppele
      conelm(4) = lowele
      do i = 1, 4
         do j = i+1, 4
            if( conelm(i).eq.conelm(j) .and. conelm(i).gt.0 ) then
               print*,'*** EDGELS: Got same face twice'
               print*,ielm,iptr,curelm,curadd
               print*,conelm
               print*,edglst
               ierr = -554
               return
            end if
         end do
      end do
C
      OPSTED = EDGLST(OPSTED)
      UPPED1 = EDGLST(UPPED1)
      LOWED1 = EDGLST(LOWED1)
      UPPED2 = EDGLST(UPPED2)
      LOWED2 = EDGLST(LOWED2)
c
      edglst(1) = iptr
      edglst(2) = opsted
      edglst(3) = upped1
      edglst(4) = lowed1
      edglst(5) = upped2
      edglst(6) = lowed2
      do i = 1, 6
         do j = i+1, 6
            if( edglst(i) .eq. edglst(j) ) then
               print*,'*** EDGELS: Got same edge twice'
               print*,ielm,iptr,curelm,curadd
               print*,conelm
               print*,edglst
               ierr = -553
               return
            end if
         end do
      end do
C
      CURADD = CURADD + 1
      ELPTRS(1,CURADD) = CURELM
      EDPTRS(3,CURADD) = OPSTED
C
      IF( LSTELE .EQ. 0 ) THEN
         INTRNL = IREG
         FSTREG = IREG
         PRVREG = IREG
         IF( CHKEL1 .LT. 0 ) THEN
            LSTELE = CHKEL1
            LSTUED = UPPED1
            LSTLED = LOWED1
            SURFAC = -CHKEL1
            UPPEND = LSTUED
            LOWEND = LSTLED
         ELSE
            LSTELE = CHKEL2
            LSTUED = UPPED2
            LSTLED = LOWED2
            SWCHEL = LSTELE
            SWUPED = LSTUED
            SWLWED = LSTLED
            IF( CHKEL2 .LT. 0 ) THEN
               SURFAC = -CHKEL2
               UPPEND = UPPED2
               LOWEND = LOWED2
            END IF
         END IF
      ELSE IF( IREG .NE. PRVREG ) THEN
         IF( SURFAC .GT. 0 ) THEN
            INTRNL = -1
         ELSE IF( INTRNL .GE. 0 ) THEN
            INTRNL = -CURADD
            SURFAC = -CURADD
         ELSE IF( IREG .NE. FSTREG ) THEN
            INTRNL = -1
         END IF
      ELSE IF( INTRNL .LT. -1 .AND. SURFAC .LT. 0 .AND.
     :         IREG .NE. FSTREG ) THEN
         SURFAC = SURFAC - 1
      END IF
C
      PRVREG = IREG
C
      IF( CHKEL1 .EQ. LSTELE ) THEN
         CHKEL1 = CHKEL2
         CHKEL2 = LSTELE
         I = UPPED1
         UPPED1 = UPPED2
         UPPED2 = I
         I = LOWED1
         LOWED1 = LOWED2
         LOWED2 = I
      END IF
C
      IF( LSTUED .EQ. LOWED2 ) THEN
         ED1 = LOWED2
         LOWED2 = UPPED2
         UPPED2 = ED1
         ED1 = LOWED1
         LOWED1 = UPPED1
         UPPED1 = ED1
         ED1 = UPPELE
         UPPELE = LOWELE
         LOWELE = ED1
      ELSE IF( LSTLED .NE. LOWED2 .OR. LSTUED .NE. UPPED2 ) THEN
         PRINT*,'*** EDGELS: Did not find last edge!'
         print*,ielm,iptr,curelm,surfac,curadd
         print*,chkel1,chkel2,lstele
         print*,upped1,upped2,lstued
         print*,lowed1,lowed2,lstled
         print*,(elptrs(1,i),i=1,curadd)
         print*,(elptrs(2,i),i=1,curadd)
         print*,(elptrs(3,i),i=1,curadd)
         print*,(edptrs(1,i),i=1,curadd)
         print*,(edptrs(2,i),i=1,curadd)
         print*,(edptrs(3,i),i=1,curadd)
         ierr = -553
         return
      END IF
C
      ELPTRS(2,CURADD) = UPPELE
      ELPTRS(3,CURADD) = LOWELE
      LSTUED = UPPED1
      LSTLED = LOWED1
      EDPTRS(1,CURADD) = UPPED1
      EDPTRS(2,CURADD) = LOWED1
C
      IF( CHKEL1 .LT. 0 ) THEN
C
         IF( SURFAC .LE. 0 ) THEN
C
            IF( ISINSD ) THEN
               PRINT*,'*** EDGELS: INCONSISTENCY IN EDGE FLAGS'
               print*,curelm,iptr,ielm,curadd
               print*,conelm
               print*,edglst
               IERR = -555
               RETURN
            END IF
C
            IF( INTRNL .NE. IREG ) INTRNL = -1
            SURFAC = -CHKEL1
C
            DO CURCHK = 1, CURADD/2
               ED1 = ELPTRS(1,CURCHK)
               ELPTRS(1,CURCHK) = ELPTRS(1,CURADD-CURCHK+1)
               ELPTRS(1,CURADD-CURCHK+1) = ED1
               ED1 = ELPTRS(2,CURCHK)
               ELPTRS(2,CURCHK) = ELPTRS(2,CURADD-CURCHK+1)
               ELPTRS(2,CURADD-CURCHK+1) = ED1
               ED1 = ELPTRS(3,CURCHK)
               ELPTRS(3,CURCHK) = ELPTRS(3,CURADD-CURCHK+1)
               ELPTRS(3,CURADD-CURCHK+1) = ED1
               ED1 = EDPTRS(3,CURCHK)
               EDPTRS(3,CURCHK) = EDPTRS(3,CURADD-CURCHK+1)
               EDPTRS(3,CURADD-CURCHK+1) = ED1
            END DO
C
            DO CURCHK = 1, (CURADD-1)/2
               ED1 = EDPTRS(1,CURCHK)
               ED2 = EDPTRS(2,CURCHK)
               EDPTRS(1,CURCHK) = EDPTRS(1,CURADD-CURCHK)
               EDPTRS(2,CURCHK) = EDPTRS(2,CURADD-CURCHK)
               EDPTRS(1,CURADD-CURCHK) = ED1
               EDPTRS(2,CURADD-CURCHK) = ED2
            END DO
C
            UPPEND = EDPTRS(1,CURADD)
            LOWEND = EDPTRS(2,CURADD)
            EDPTRS(1,CURADD) = SWUPED
            EDPTRS(2,CURADD) = SWLWED
            LSTUED = SWUPED
            LSTLED = SWLWED
C
            LSTELE = IELM
            CURELM = SWCHEL
            IF( SWCHEL .GT. 0 ) GOTO 5
C
            PRINT*,'*** EDGELS: Got a surface after a switch?'
            print*,ielm,iptr,swchel,surfac
            ierr = -551
            return
C
         END IF
C
         EDPTRS(1,CURADD+1) = UPPEND
         EDPTRS(2,CURADD+1) = LOWEND
C
         GOTO 400
C
      END IF
C
      LSTELE = CURELM
      CURELM = CHKEL1
C
      IF( CURELM .NE. IELM ) GOTO 5
C
      IF( SURFAC .GT. 0 ) THEN
         print*,'*** EDGELS: Got first elem again with surface'
         print*,ielm,iptr,curelm,surfac
         ierr = -552
         return
      END IF
C
 400  CONTINUE
C
      CALL EDGINF( BIGLST, NODLST, IPTR, NODS, INEWED, IFLAG, NXT, LST )
C
      DO CURCHK = 1, CURADD
         CALL CMNNOD( BIGLST, NODLST, EDPTRS(1,CURCHK),
     :                EDPTRS(2,CURCHK), INOD, I1, I2 )
         IF( I2 .EQ. NODS(1) ) THEN
            ED1 = EDPTRS(2,CURCHK)
            EDPTRS(2,CURCHK) = EDPTRS(1,CURCHK)
            EDPTRS(1,CURCHK) = ED1
            ED1 = ELPTRS(2,CURCHK)
            ELPTRS(2,CURCHK) = ELPTRS(3,CURCHK)
            ELPTRS(3,CURCHK) = ED1
         END IF
         IF( INOD .EQ. 0 ) THEN
            PRINT*,'*** EDGELS: FAILED TO FIND COMMON NODE!'
            PRINT*,'  edges: ',edptrs(1,curchk),edptrs(2,curchk)
            inod = edptrs(1,curchk)
            print*,'  nods1: ',biglst(1,inod),biglst(2,inod)
            inod = edptrs(2,curchk)
            print*,'  nods2: ',biglst(1,inod),biglst(2,inod)
            IERR = -10201
            RETURN
         END IF
         NDPTRS(CURCHK) = INOD
      END DO
C
      IF( SURFAC .GT. 0 ) THEN
         CALL CMNNOD( BIGLST, NODLST, EDPTRS(1,CURADD+1),
     :                EDPTRS(2,CURADD+1), INOD, I1, I2 )
         IF( I2 .EQ. NODS(1) ) THEN
            ED1 = EDPTRS(2,CURADD+1)
            EDPTRS(2,CURADD+1) = EDPTRS(1,CURADD+1)
            EDPTRS(1,CURADD+1) = ED1
         END IF
         IF( INOD .EQ. 0 ) THEN
            PRINT*,'*** EDGELS: FAILED TO FIND COMMON SURFACE NODE!'
            PRINT*,'  edges: ',edptrs(1,curadd+1),edptrs(2,curadd+1)
            inod = edptrs(1,curadd+1)
            print*,'  nods1: ',biglst(1,inod),biglst(2,inod)
            inod = edptrs(2,curadd+1)
            print*,'  nods2: ',biglst(1,inod),biglst(2,inod)
            IERR = -10202
            RETURN
         END IF
         NDPTRS(CURADD+1) = INOD
      END IF
C
 900  if( goterr ) then
       print*,'Listing info: ',ielm,iptr,surfac,intrnl
       print*,(elptrs(1,i),i=1,curadd)
       print*,(elptrs(2,i),i=1,curadd)
       print*,(elptrs(3,i),i=1,curadd)
       curchk = curadd
       if( surfac .gt. 0 ) curchk = curchk + 1
       print*,(edptrs(1,i),i=1,curchk)
       print*,(edptrs(2,i),i=1,curchk)
       print*,(edptrs(3,i),i=1,curadd)
      end if
c
      RETURN
      END
C

