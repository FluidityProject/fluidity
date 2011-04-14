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
      SUBROUTINE GETHED( HEDFIL, MSHFIL, MSHBIN,
     :                   NNOD, NELM, NGEOM, NCON )
C  ---------------------------------------------------------------------
C
C
C
C  ---------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER HEDFIL, MSHFIL, NNOD, NELM, NCON, NGEOM
C
      LOGICAL MSHBIN
C
      INTEGER MF, MFO, KH, I1, I2, I3, I4, I5, I6, I7, I8, I9
C
      CHARACTER*128 TITLE, FILNAM, HOME
C
      CALL PXFGETENV( 'EVENT_HOME', 10, HOME, KH, I1 )
C
      READ( HEDFIL, '(A)' ) TITLE
      READ( HEDFIL, * )
C
      READ( HEDFIL, '(10X,2I5)' ) NNOD, NGEOM
      READ( HEDFIL, * )
C
      IF( MOD( NNOD, 10 ) .EQ. 3 ) READ( HEDFIL, * )
      READ( HEDFIL, * )
C
      READ( HEDFIL, '(10X,A)' ) FILNAM
      MF  = 128
C
  5   MFO = MF
      CALL STRPQT( FILNAM, MF )
      MF = INDEX( FILNAM, ' ' ) - 1
      IF( MF .NE. MFO ) GOTO 5
C
      IF( MF .GT. 1 ) THEN
         IF( FILNAM(1:1) .NE. '/' ) THEN
            FILNAM = HOME(1:KH)//'/data/mesh/'//FILNAM(1:MF)
            MF = INDEX( FILNAM, ' ' ) - 1
         END IF
         OPEN( UNIT=MSHFIL, FILE=FILNAM(1:MF), STATUS='OLD',
     :         FORM='FORMATTED' )
      ELSE
         MSHFIL = HEDFIL
      END IF
C
      MSHBIN = .FALSE.
      READ( MSHFIL, '(25X,3I5)', ERR=10 ) NELM, NNOD, NCON
      GOTO 20
C
  10  MSHBIN = .TRUE.
      IF( MSHFIL .EQ. HEDFIL ) THEN
         WRITE(0,*) 'ERROR READING FROM EVENT DATA FILE'
         NGEOM = -1
         RETURN
      END IF
C
      OPEN( UNIT=MSHFIL, FILE=FILNAM(1:MF), STATUS='OLD',
     :      FORM='UNFORMATTED' )
C
      READ( MSHFIL ) I1, I2, I3, NELM, NNOD, NCON, I4, I5, I6,
     :               I7, I8, I9
C
  20  CONTINUE
C
C
      RETURN
      END
C
