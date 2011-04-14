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

C> This subroutine makes a guess at the sizes required for the real and
C> integer arrays in ADPTVY, returning the guesses in SZRL and SZINT.
C>
C> It assumes that NNOD, NELM and TOTFRE (the number of fields) are known,
C> plus at least one of XPCTEL or XPCTND (expected elements or nodes).
C>
C> An error results in negative return values for SZRL and SZINT.
C>
C> Note that only the first three (NNOD,NELM,TOTFRE) can be passed down as
C> a constant, since any of the others which are given non-positive values
C> will be returned as a different value (as outlined below). 
C>
C> If one of XPCTEL or XPCTND is not provided (i.e. <=0) then it is worked
C> out based upon the ratio of NELM to NNOD in the current mesh (which may
C> or may not be a good estimate of what's to come...)
C>
C> If SZENLS is provided (i.e. >0) it is used, else assumed to be 4*NELM.
C>
C> If XPCTSE (expected no. of surf.els) is provided (i.e. >0) it is used,
C> otherwise a guess is made based on the ratio of NSELM to NELM - note
C> that if NSELM is not provided (i.e. <=0) then XPCTSE is taken to be
C> half of XPCTEL (a quite high upper bound).
C>
C>
C> The calculation of the required sizes is performed using the following
C> rough guidelines...
C>
C> Integer array must be large enough to contain the following...
C>
C> - node-element list base pointer for old mesh: NNOD+1
C> - node-element list for old mesh: SZENLS
C> - element connectivity list for the old mesh: 4*NELM (in 3D)
C> - the adapting edge/element info: >30*NWNELM, probably >35*NWNELM
C> - the new elementwise region info: NWNELM
C> - the new volume element-node list base pointer: NWNELM+1
C> - the new volume element-node list: 4*NWNELM (in 3D)
C> - the new surface element-node list base pointer: NWNSEL+1
C> - the new surface element-node list: 3*NWNSEL (in 3D)
C>
C> Real array must be large enough to contain the following...
C>
C> - the adapting node info: >17*NWNNOD, probably >20*NWNNOD
C> - the new node coordinates: 3*NWNNOD
C> - the interpolated field: TOTFRE*NWNNOD
C> - the interpolated metrics (if METRIC is TRUE): 9*NWNNOD

C> @param[in] NNOD Number of mesh nodes.
C> @param[in] NELM Number of mesh elements.
C> @param[in] SZENLS Size of the element-node list.
C> @param[in] NSELM Number of mesh surface elements.
C> @param[in] TOTFRE Number of degrees of freedom.
C> @param[in] XPCTEL Expected number of elements after an adapt.
C> @param[in] XPCTND Expected number of nodes after an adapt.
C> @param[in] XPCTSE Expected number of surface elements after an adapt.
C> @param[in] METRIC Is the metric to be interpolated?
C> @param[out] SZINT Size of the integer array required by adaptivity.
C> @param[out] SZRL Size of the real array required by adaptivity.

      SUBROUTINE ADAPTMEM( NNOD, NELM, SZENLS, NSELM, TOTFRE,
     &                     XPCTEL, XPCTND, XPCTSE,
     &                     METRIC, SZINT, SZRL )
      IMPLICIT NONE
C
      INTEGER NNOD, NELM, SZENLS, NSELM, TOTFRE, XPCTEL, XPCTND,
     &        XPCTSE, SZINT, SZRL
C
      LOGICAL METRIC
C
C - use this long int to deal with any problems with large numbers in the
C - multiplication below..
C
      INTEGER*8 P
      INTEGER MAXEL, MAXND, MAXSE
      REAL FACTOR
      INCLUDE 'blkerr.i'
C
      SZRL  = -1
      SZINT = -1
C
C - make a guess about the expected no. of surface elements
C
      IF( XPCTEL .LE. 0 ) THEN
         IF( XPCTND .LE. 0 ) THEN
            WRITE(0,*)
     &         '*** ADAPTMEM: Need at least one of XPCTEL or XPCTND > 0'
            WRITE(0,*) 'XPCTEL, XPCTND: ',XPCTEL,XPCTND
            RETURN
         END IF
         P = 1
         P = P*NELM*XPCTND/NNOD
         XPCTEL = P
      END IF
C
      IF( NSELM .LE. 0 ) THEN
         XPCTSE = XPCTEL/2
      ELSE IF( XPCTSE .LE. 0 ) THEN
         P = 1
         P = P*NSELM*XPCTEL/NELM
         XPCTSE = P
      END IF
C
      IF( XPCTND .LE. 0 ) THEN
         P = 1
         P = P*NNOD*XPCTEL/NELM
         XPCTND = P
      END IF
C
      IF( SZENLS .LE. 0 ) THEN
         SZENLS = 4*NELM
      END IF
C
C - make sure we take the maximum nodes, etc, and build in some leeway
C
      MAXEL = MAX(NELM,XPCTEL)*4./3.
      MAXND = MAX(NNOD,XPCTND)*4./3.
      MAXSE = MAX(NSELM,XPCTSE)*4./3.
C
      SZINT = NNOD+1 + SZENLS + 4*NELM + 41*MAXEL+1 + 4*MAXSE+1
C
      SZRL = (23+TOTFRE)*MAXND
C
      IF( METRIC ) SZRL = SZRL + 9*MAXND

      RETURN
      END SUBROUTINE ADAPTMEM
C
