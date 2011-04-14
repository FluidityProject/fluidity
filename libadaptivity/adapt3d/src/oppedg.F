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
      SUBROUTINE OPPEDG( BIGLST, NODLST, IEDG, IELM, OPEDG, NXT, LST )
C-----------------------------------------------------------------------
C
C
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
      INTEGER IEDG, IELM, OPEDG, NXT, LST
C
      INTEGER EDGLST(6), CONELS(4), NODS(4), IFLAG, INEW, NXT1, LST1,
     :        IREG, I
C
      LOGICAL chkcns
C
      IF( IERR .NE. 0 ) RETURN
C
ccc      chkcns = .false.
C
      CALL ELMINF( BIGLST, NODLST, IELM, CONELS, EDGLST,
     :             INEW, IFLAG, IREG, NXT, LST )
C
cc      IF( IERR .NE. 0 ) THEN
cc        print*,'*** OPPEDG: GOT ERROR FROM ELMINF'
cc         print*,'   edge,element: ',iedg,ielm
cc         RETURN
cc      END IF
C
      if( edglst(1) .eq. iedg ) then
         opedg = edglst(5)
      else if( edglst(2) .eq. iedg ) then
         opedg = edglst(6)
      else if( edglst(3) .eq. iedg ) then
         opedg = edglst(4)
      else if( edglst(4) .eq. iedg ) then
         opedg = edglst(3)
      else if( edglst(5) .eq. iedg ) then
         opedg = edglst(1)
      else if( edglst(6) .eq. iedg ) then
         opedg = edglst(2)
ccc
ccc      CALL EDGINF( BIGLST, NODLST, IEDG, NODS,
ccc     :             INEW, IFLAG, NXT1, LST1 )
C
cc      IF( IERR .NE. 0 ) THEN
cc         print*,'*** OPPEDG: GOT ERROR FROM EDGINF'
cc         print*,'   edge,element: ',iedg,ielm
cc         print*,'   conelm: ',conels
cc         print*,'   edglst: ',edglst
cc         RETURN
cc      END IF
C
ccc      I = 0
C
ccc  10  IF( I .LT. 6 ) THEN
ccc         I = I + 1
ccc         OPEDG = EDGLST(I)
ccc         IF( BIGLST(1,OPEDG) .EQ. NODS(1) ) GOTO 10
ccc         IF( BIGLST(2,OPEDG) .EQ. NODS(1) ) GOTO 10
ccc         IF( BIGLST(1,OPEDG) .EQ. NODS(2) ) GOTO 10
ccc         IF( BIGLST(2,OPEDG) .EQ. NODS(2) ) GOTO 10
c         if( iedg.eq.1 .and. i.ne.5 ) then
c            print*,'+++ OPPEDG: inconsistency: ',iedg,i,opedg
c         else if( iedg.eq.2 .and. i.ne.6 ) then
c            print*,'+++ OPPEDG: inconsistency: ',iedg,i,opedg
c         else if( iedg.eq.3 .and. i.ne.4 ) then
c            print*,'+++ OPPEDG: inconsistency: ',iedg,i,opedg
c         else if( iedg.eq.4 .and. i.ne.3 ) then
c            print*,'+++ OPPEDG: inconsistency: ',iedg,i,opedg
c         else if( iedg.eq.5 .and. i.ne.1 ) then
c            print*,'+++ OPPEDG: inconsistency: ',iedg,i,opedg
c         else if( iedg.eq.6 .and. i.ne.2 ) then
c            print*,'+++ OPPEDG: inconsistency: ',iedg,i,opedg
c         endif 
      ELSE
         OPEDG = 0
         PRINT*,'*** OPPEDG: COULD NOT FIND OPPOSITE EDGE'
         PRINT*,'  Edge: ',iedg
         do i = 1, 6
            print*,edglst(i),biglst(1,edglst(i)),biglst(2,edglst(i))
         end do
         IERR = -150
      END IF
c
c - this does a consistency check to see if the rest of the
c - element's edges have the nodes of the supplied and opposite edges
c
ccc      if( chkcns ) then
c
ccc        CALL EDGINF( BIGLST, NODLST, OPEDG, NODS(3),
ccc     :               INEW, IFLAG, NXT1, LST1 )
c
ccc        do i = 1, 6
c
ccc          inew  = biglst(1,edglst(i))
ccc          iflag = biglst(2,edglst(i))
c
ccc          if( inew .ne. nods(1) .and.
ccc     :        inew .ne. nods(2) .and.
ccc     :        inew .ne. nods(3) .and.
ccc     :        inew .ne. nods(4) ) then
ccc            print*,'*** OPPEDG: EDGE HAS NODE NOT IN ELEMENT!'
ccc            print*,'   elm: ',ielm,ireg,i
ccc            print*,'   nodes: ',nods
ccc            print*,'   conel: ',conels
ccc            print*,'   edges: ',edglst
ccc            print*,'   ednds: ',inew,iflag,edglst(i)
ccc            ierr = -11112
ccc            return
ccc          end if
c
ccc          if( iflag .ne. nods(1) .and.
ccc     :        iflag .ne. nods(2) .and.
ccc     :        iflag .ne. nods(3) .and.
ccc     :        iflag .ne. nods(4) ) then
ccc            print*,'*** OPPEDG: EDGE HAS NODE NOT IN ELEMENT!'
ccc            print*,'   elm: ',ielm,ireg,i
ccc            print*,'   nodes: ',nods
ccc            print*,'   conel: ',conels
ccc            print*,'   edges: ',edglst
ccc            print*,'   ednds: ',inew,iflag,edglst(i)
ccc            ierr = -11112
ccc            return
ccc          end if
c
ccc        end do
c
ccc      end if
C
      RETURN
      END
C
