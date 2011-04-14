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
      SUBROUTINE CNCTEL( BIGLST, NODLST, ELM1, ELM2 )
C-----------------------------------------------------------------------
C
C - This subroutine adds ELM2 to the list of connected elements for ELM1,
C - and vice versa. If they already exist in the lists then it does nothing.
C
C - ELM2 can be sent down as a negative surface value, in which case
C - this subroutine will connect ELM1 to that surface. If element ELM1
C - is already connected to that surface, then this subroutine attempts
C - to make another face of ELM1 connect to the surface. If there is no
C - unconnected face (with zero value) then an error will occur.
C - i.e. you can send down ELM2 as a preconnected element to ELM1 and
C - there will be no action taken (as above), BUT you cannot expect the
C - same for the surfaces. This is because it is possible for an element
C - to have more than one face on a surface but not more than one face
C - connected to only one other element.
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
      INTEGER ELM1, ELM2, i
c
c      logical debug
c
c      debug = .false.
c      if( elm1 .eq. 57095 .or. elm2 .eq. 57095 ) debug = .true.
c
c      if( debug ) then
c        print*,'Into CNCTEL: ',elm1,elm2,ierr
c        print*,'   conelm for 1: ',(biglst(i,elm1),i=1,4)
c        print*,'   conelm for 2: ',(biglst(i,elm1),i=1,4)
c      end if
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( ELM1 .LE. 0 ) THEN
         PRINT*,'*** CNCTEL: ELM1 IS NOT A VALID POINTER'
         print*,elm1,elm2
         IERR = -1
         RETURN
      END IF
C
      IF( BIGLST(4,ELM1) .LT. 0 ) THEN
         PRINT*,'*** CNCTEL: ELM1 IS NOT AN ELEMENT'
         PRINT*,ELM1,elm2,BIGLST(4,ELM1)
         IERR = -2
         RETURN
      END IF
C
      IF( ELM2 .EQ. 0 ) THEN
C
         PRINT*,'*** CNCTEL: ELM2 IS NOT A VALID POINTER'
         print*,elm1,elm2
         IERR = -1
         RETURN
C
      ELSE IF( ELM2 .GT. 0 ) THEN
C
         IF( BIGLST(4,ELM2) .LT. 0 ) THEN
            PRINT*,'*** CNCTEL: ELM2 IS NOT AN ELEMENT'
            PRINT*,elm1,ELM2,BIGLST(4,ELM2)
            IERR = -2
            RETURN
         END IF
C
         IF( BIGLST(1,ELM1) .EQ. ELM2 ) GOTO 100
         IF( BIGLST(2,ELM1) .EQ. ELM2 ) GOTO 100
         IF( BIGLST(3,ELM1) .EQ. ELM2 ) GOTO 100
         IF( BIGLST(4,ELM1)/2 .EQ. ELM2 ) GOTO 100
C
      ELSE IF( BIGLST(4,ELM1) .EQ. 1 ) THEN
C
         IERR = 0
C
         IF( BIGLST(1,ELM1) .GT. 0 ) THEN
            BIGLST(4,ELM1) = BIGLST(1,ELM1)*2 + 1
            BIGLST(1,ELM1) = ELM2
         ELSE IF( BIGLST(2,ELM1) .GT. 0 ) THEN
            BIGLST(4,ELM1) = BIGLST(2,ELM1)*2 + 1
            BIGLST(2,ELM1) = ELM2
         ELSE IF( BIGLST(3,ELM1) .GT. 0 ) THEN
            BIGLST(4,ELM1) = BIGLST(3,ELM1)*2 + 1
            BIGLST(3,ELM1) = ELM2
         ELSE
            PRINT*,'*** CNCTEL: ELM1 CANNOT HAVE ALL SIDES ON SURFACES'
            IERR = -98
         END IF
C
         GOTO 200
C
      END IF
C
      IERR = 0
C
      IF( BIGLST(1,ELM1) .EQ. 0 ) THEN
         BIGLST(1,ELM1) = ELM2
      ELSE IF( BIGLST(2,ELM1) .EQ. 0 ) THEN
         BIGLST(2,ELM1) = ELM2
      ELSE IF( BIGLST(3,ELM1) .EQ. 0 ) THEN
         BIGLST(3,ELM1) = ELM2
      ELSE IF( BIGLST(4,ELM1)/2 .EQ. 0 ) THEN
         BIGLST(4,ELM1) = ELM2*2 + 1
      ELSE
         PRINT*,'*** CNCTEL: ELM1 IS ALREADY FULLY CONNECTED'
         print*,elm2,elm1,(biglst(i,elm1),i=1,4)
         IERR = -25
      END IF
C
 100  CONTINUE
C
      IF( ELM2 .LT. 0 ) GOTO 200
C
      IF( BIGLST(1,ELM2) .EQ. ELM1 ) GOTO 200
      IF( BIGLST(2,ELM2) .EQ. ELM1 ) GOTO 200
      IF( BIGLST(3,ELM2) .EQ. ELM1 ) GOTO 200
      IF( BIGLST(4,ELM2)/2 .EQ. ELM1 ) GOTO 200
C
      IF( BIGLST(1,ELM2) .EQ. 0 ) THEN
         BIGLST(1,ELM2) = ELM1
      ELSE IF( BIGLST(2,ELM2) .EQ. 0 ) THEN
         BIGLST(2,ELM2) = ELM1
      ELSE IF( BIGLST(3,ELM2) .EQ. 0 ) THEN
         BIGLST(3,ELM2) = ELM1
      ELSE IF( BIGLST(4,ELM2) .EQ. 1 ) THEN
         BIGLST(4,ELM2) = ELM1*2 + 1
      ELSE
         PRINT*,'*** CNCTEL: ELM2 IS ALREADY FULLY CONNECTED'
         print*,elm1,elm2,(biglst(i,elm2),i=1,4)
         IERR = -25
      END IF
C
 200  CONTINUE
c
c      if( debug ) then
c        print*,'Leaving CNCTEL: ',elm1,elm2,ierr
c        print*,'   conelm for 1: ',(biglst(i,elm1),i=1,4)
c        print*,'   conelm for 2: ',(biglst(i,elm1),i=1,4)
c      end if
C
      RETURN
      END
C
