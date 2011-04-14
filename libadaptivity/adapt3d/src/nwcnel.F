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
      SUBROUTINE NWCNEL( BIGLST, NODLST, ELM1, ELM2,
     :                   IFAC, IED1, IED2, IED3 )
C-----------------------------------------------------------------------
C
C - This subroutine adds ELM2 to the list of connected elements for ELM1,
C - and vice versa. If they already exist in the lists (in the correct
C - positions) then it does nothing.
C
C - ELM2 can be sent down as a negative surface value, in which case
C - this subroutine will connect ELM1 to that surface.
C
C - IFAC is the face of ELM1 that ELM2 should connect to. A zero value
C - will cause the edges (IED1,IED2,IED3) to be used to find the face,
C - otherwise the edges are ignored. If IFAC is zero, and the edges are
c - zero, then it will try to find the common edges between ELM1 and ELM2
C - and use those to find the appropriate faces to connect up.
C - If the face (IFAC) and edges (IED1,IED2,IED3) are all zero, and ELM2
C - is a surface (i.e. negative) then ... what...???
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
      INTEGER ELM1, ELM2, IFAC, IED1, IED2, IED3
C
      INTEGER I, WCHFAC, IFC2, NXT, NXT2, J, ieo1, ieo2, ieo3, ifco
C
      ieo1 = ied1
      ieo2 = ied2
      ieo3 = ied3
      ifco = ifac
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
         WRITE(0,*) '*** CNCTEL: ELM1 IS NOT A VALID POINTER'
         WRITE(0,*) elm1,elm2
         IERR = -1
         RETURN
      END IF
C
      IF( BIGLST(4,ELM1) .LT. 0 ) THEN
         WRITE(0,*) '*** CNCTEL: ELM1 IS NOT AN ELEMENT'
         WRITE(0,*) ELM1,elm2,BIGLST(4,ELM1)
         IERR = -2
         RETURN
      END IF
C
      IF( ELM2 .EQ. 0 ) THEN
C
         WRITE(0,*) '*** CNCTEL: ELM2 IS NOT A VALID POINTER'
         WRITE(0,*) elm1,elm2
         IERR = -1
         RETURN
C
      ELSE IF( ELM2 .GT. 0 ) THEN
C
         IF( BIGLST(4,ELM2) .LT. 0 ) THEN
            WRITE(0,*) '*** CNCTEL: ELM2 IS NOT AN ELEMENT'
            WRITE(0,*) elm1,ELM2,BIGLST(4,ELM2)
            IERR = -2
            RETURN
         END IF
C
c         IF( BIGLST(1,ELM1) .EQ. ELM2 ) GOTO 100
c         IF( BIGLST(2,ELM1) .EQ. ELM2 ) GOTO 100
c         IF( BIGLST(3,ELM1) .EQ. ELM2 ) GOTO 100
c         IF( BIGLST(4,ELM1)/2 .EQ. ELM2 ) GOTO 100
C
      END IF
C
      IF( IFAC .LT. 0 ) THEN
C
         WRITE(0,*) '*** CNCTEL: INVALID FACE: ',IFAC
         WRITE(0,*) ELM1,ELM2
         IERR = -2013
         RETURN
C
      ELSE IF( IFAC .GT. 4 ) THEN
C
         WRITE(0,*) '*** CNCTEL: INVALID FACE: ',IFAC
         WRITE(0,*) ELM1,ELM2
         IERR = -2013
         RETURN
C
      ELSE IF( IFAC .EQ. 0 ) THEN
C
C - no face given, so have to find it (from edges, if given)
C
         IF( IED1 .EQ. 0 ) THEN
C
C - no edges given either!
C
            NXT = BIGLST(NXTBIG,ELM1)
            NXT2 = BIGLST(NXTBIG,NXT)
C
            IF( ELM2 .GT. 0 ) THEN
C - find edges common to ELM1 and ELM2
               CALL CMNEDS( BIGLST, NODLST, ELM1, ELM2,
     :                      IED1, IED2, IED3 )
C
               IF( IERR .NE. 0 ) THEN
                  WRITE(0,*) '*** NWCNEL: Got error from CMNEDS (1)'
                  RETURN
               END IF
C
               IF( IED1.EQ.0 .OR. IED2.EQ.0 .OR. IED3.EQ.0 ) THEN
                  WRITE(0,*) '*** NWCNEL: Lost edge after CMNEDS!'
                  WRITE(0,*) elm1,elm2
                  WRITE(0,*) ied1,ied2,ied3
                  WRITE(0,*) (biglst(j,elm1),j=1,4)
                  WRITE(0,*) (biglst(j,nxt),j=1,4),
     :                       (biglst(j,nxt2),j=1,2)
                  if(elm2.gt.0) then
                    WRITE(0,*) (biglst(j,elm2),j=1,4)
                    nxt = biglst(nxtbig,elm2)
                    nxt2 = biglst(nxtbig,nxt)
                    WRITE(0,*) (biglst(j,nxt),j=1,4),
     :                     (biglst(j,nxt2),j=1,2)
                  endif
                  IERR = -6547
                  return
               END IF
C
            ELSE
c
c - what to do here????????????
c
c - ELM2 is a surface, but we do not have any info about which face
c - of ELM1 should have this surface...
C
C - just look for face entry which contains zero (it should only get
C - into this situation from SHREDG, and then there should be only one
C - entry in the faces that contains zero...)
C
               IF(BIGLST(1,ELM1) .EQ. 0 .OR.
     :            BIGLST(1,ELM1) .EQ. ELM2 ) THEN
                  IFAC = 1
                  IED1 = BIGLST(1,NXT)
                  IED2 = BIGLST(2,NXT)
                  IED3 = BIGLST(3,NXT)
               ELSE IF(BIGLST(2,ELM1) .EQ. 0 .OR.
     :                 BIGLST(2,ELM1) .EQ. ELM2 ) THEN
                  IFAC = 2
                  IED1 = BIGLST(1,NXT)
                  IED2 = BIGLST(4,NXT)
                  IED3 = BIGLST(2,NXT2)
               ELSE IF(BIGLST(3,ELM1) .EQ. 0 .OR.
     :                 BIGLST(3,ELM1) .EQ. ELM2 ) THEN
                  IFAC = 3
                  IED1 = BIGLST(3,NXT)
                  IED2 = BIGLST(1,NXT2)
                  IED3 = BIGLST(2,NXT2)
               ELSE IF(BIGLST(4,ELM1) .LT. 2 .OR.
     :                 BIGLST(4,ELM1)/2 .EQ. ELM2 ) THEN
                  IFAC = 4
                  IED1 = BIGLST(2,NXT)
                  IED2 = BIGLST(4,NXT)
                  IED3 = BIGLST(1,NXT2)
               ELSE
                  WRITE(0,*) '*** NWCNEL: Could not find empty face!'
                  WRITE(0,*) elm1,elm2
                  WRITE(0,*) ied1,ied2,ied3
                  WRITE(0,*) (biglst(i,elm1),i=1,4)
                  WRITE(0,*) (biglst(i,nxt),i=1,4),
     :                       (biglst(i,nxt2),i=1,2)
                  if(elm2.gt.0) then
                    WRITE(0,*) (biglst(i,elm2),i=1,4)
                    nxt = biglst(nxtbig,elm2)
                    nxt2 = biglst(nxtbig,nxt)
                    WRITE(0,*) (biglst(i,nxt),i=1,4),
     :                         (biglst(i,nxt2),i=1,2)
                  endif
                  IERR = -6544
                  RETURN
               END IF
C
C - check that was the only empty face
C
               IF(IFAC .LT. 3) THEN
C
                  DO I = IFAC+1, 3
                     IF(BIGLST(I,ELM1) .EQ. 0 ) THEN
                        WRITE(0,*)
     :                      '*** NWCNEL: Found another empty face!'
                        WRITE(0,*) elm1,elm2
                        WRITE(0,*) ied1,ied2,ied3
                        WRITE(0,*) (biglst(j,elm1),j=1,4)
                        WRITE(0,*) (biglst(j,nxt),j=1,4),
     :                         (biglst(j,nxt2),j=1,2)
                        if(elm2.gt.0) then
                          WRITE(0,*) (biglst(j,elm2),j=1,4)
                          nxt = biglst(nxtbig,elm2)
                          nxt2 = biglst(nxtbig,nxt)
                          WRITE(0,*) (biglst(j,nxt),j=1,4),
     :                               (biglst(j,nxt2),j=1,2)
                        endif
                        IERR = -6545
                        RETURN
                     END IF
                  END DO
C
                  IF(BIGLST(4,ELM1) .LT. 2) THEN
                     WRITE(0,*) '*** NWCNEL: Found another empty face!'
                     WRITE(0,*) elm1,elm2
                     WRITE(0,*) ied1,ied2,ied3
                     WRITE(0,*) (biglst(j,elm1),j=1,4)
                     WRITE(0,*) (biglst(j,nxt),j=1,4),
     :                          (biglst(j,nxt2),j=1,2)
                     if(elm2.gt.0) then
                       WRITE(0,*) (biglst(j,elm2),j=1,4)
                       nxt = biglst(nxtbig,elm2)
                       nxt2 = biglst(nxtbig,nxt)
                       WRITE(0,*) (biglst(j,nxt),j=1,4),
     :                            (biglst(j,nxt2),j=1,2)
                     endif
                     IERR = -6545
                     RETURN
                  END IF
C
               ELSE IF( IFAC .LT. 4 ) THEN
C
                  IF(BIGLST(4,ELM1) .LT. 2) THEN
                     WRITE(0,*) '*** NWCNEL: Found another empty face!'
                     WRITE(0,*) elm1,elm2
                     WRITE(0,*) ied1,ied2,ied3
                     WRITE(0,*) (biglst(j,elm1),j=1,4)
                     WRITE(0,*) (biglst(i,nxt),i=1,4),
     :                          (biglst(i,nxt2),i=1,2)
                     if(elm2.gt.0) then
                        WRITE(0,*) (biglst(j,elm2),j=1,4)
                        nxt = biglst(nxtbig,elm2)
                        nxt2 = biglst(nxtbig,nxt)
                        WRITE(0,*) (biglst(i,nxt),i=1,4),
     :                             (biglst(i,nxt2),i=1,2)
                     endif
                     IERR = -6545
                     RETURN
                  END IF
C
               END IF
C
               IF( IED1.EQ.0 .OR. IED2.EQ.0 .OR. IED3.EQ.0 ) THEN
                  WRITE(0,*) '*** NWCNEL: Lost edge after face search!!'
                  WRITE(0,*) elm1,elm2
                  WRITE(0,*) ied1,ied2,ied3
                  WRITE(0,*) (biglst(j,elm1),j=1,4)
                  WRITE(0,*) (biglst(j,nxt),j=1,4),
     :                       (biglst(j,nxt2),j=1,2)
                  if(elm2.gt.0) then
                    WRITE(0,*) (biglst(j,elm2),j=1,4)
                    nxt = biglst(nxtbig,elm2)
                    nxt2 = biglst(nxtbig,nxt)
                    WRITE(0,*) (biglst(j,nxt),j=1,4),
     :                         (biglst(j,nxt2),j=1,2)
                  endif
                  IERR = -6548
                  return
               END IF
C
            END IF
C
         END IF
C
C - find out which face of ELM1 has edges IED1, IED2 & ED3
C
         IFAC = WCHFAC( BIGLST, NODLST, ELM1, IED1, IED2, IED3 )
C
         IF( IERR .NE. 0 ) THEN
            WRITE(0,*) '*** NWCNEL: Got error from WCHFAC (1)'
            WRITE(0,*) 'Elms, orig fac: ',elm1,elm2,ifco
            WRITE(0,*) 'Orig eds: ',ieo1,ieo2,ieo3
            if( elm2 .gt. 0 ) then
               WRITE(0,*) 'conel2: ',(biglst(i,elm2),i=1,3),
     :                               biglst(4,elm2)/2
               nxt  = biglst(nxtbig,elm2)
               nxt2 = biglst(nxtbig,nxt)
               WRITE(0,*) 'edgls2: ',(biglst(i,nxt),i=1,4),
     :                               biglst(1,nxt2),biglst(2,nxt2)
            end if
            RETURN
         END IF
C
      END IF
C
      IF( ELM2 .GT. 0 ) THEN
C
C - connecting ELM1 to a real element (i.e. ELM2 is not a surface)
C
         IF( IED1 .EQ. 0 ) THEN
C - edges not given, so find which edges are common to ELM1 and ELM2
            CALL CMNEDS( BIGLST, NODLST, ELM1, ELM2,
     :                   IED1, IED2, IED3 )
            IF( IERR .NE. 0 ) THEN
               WRITE(0,*) '*** NWCNEL: Got error from CMNEDS (2)'
               RETURN
            END IF
            IF( IED3 .EQ. 0 ) THEN
               WRITE(0,*)
     :            '*** NWCNEL: not found all 3 edges in CMNEDS (2)'
               WRITE(0,*) elm1,elm2,ifac,ifc2
               WRITE(0,*) ied1,ied2,ied3
               WRITE(0,*) (biglst(i,elm1),i=1,4)
               nxt = biglst(nxtbig,elm1)
               nxt2 = biglst(nxtbig,nxt)
               WRITE(0,*) (biglst(i,nxt),i=1,4),(biglst(i,nxt2),i=1,2)
               WRITE(0,*) (biglst(i,elm2),i=1,4)
               nxt = biglst(nxtbig,elm2)
               nxt2 = biglst(nxtbig,nxt)
               WRITE(0,*) (biglst(i,nxt),i=1,4),(biglst(i,nxt2),i=1,2)
               IERR = -6549
               RETURN
            END IF
C - quick consistency check that these edges do correspond to face IFAC
            IFC2 = WCHFAC( BIGLST, NODLST, ELM1, IED1, IED2, IED3 )
            IF( IERR .NE. 0 ) THEN
               WRITE(0,*) '*** NWCNEL: Got error from WCHFAC (2)'
               RETURN
            END IF
            IF( IFC2 .NE. IFAC ) THEN
               WRITE(0,*) '*** NWCNEL: Inconsistent face and edges'
               WRITE(0,*) 'els, fcs: ',elm1,elm2,ifac,ifc2
               WRITE(0,*) 'common eds: ',ied1,ied2,ied3
               WRITE(0,*) 'conel1: ',(biglst(i,elm1),i=1,3),
     :                               biglst(4,elm1)/2
               nxt = biglst(nxtbig,elm1)
               nxt2 = biglst(nxtbig,nxt)
               WRITE(0,*) 'edgls1: ',(biglst(i,nxt),i=1,4),
     :                               (biglst(i,nxt2),i=1,2)
               nxt2 = biglst(1,nxt2)
               nxt = biglst(1,nxt)
               WRITE(0,*) 'nodls1: ',(biglst(i,nxt),i=1,2),
     :                               (biglst(i,nxt2),i=1,2)
               WRITE(0,*) 'conel2: ',(biglst(i,elm2),i=1,3),
     :                               biglst(4,elm2)/2
               nxt = biglst(nxtbig,elm2)
               nxt2 = biglst(nxtbig,nxt)
               WRITE(0,*) 'edgls2: ',(biglst(i,nxt),i=1,4),
     :                               (biglst(i,nxt2),i=1,2)
               nxt2 = biglst(1,nxt2)
               nxt = biglst(1,nxt)
               WRITE(0,*) 'nodls2: ',(biglst(i,nxt),i=1,2),
     :                               (biglst(i,nxt2),i=1,2)
               IERR = -6546
               RETURN
            END IF
         END IF
C
C - find which face of ELM2 has edges IED1, IED2 & IED3
C
         IFC2 = WCHFAC( BIGLST, NODLST, ELM2, IED1, IED2, IED3 )
C
         IF( IERR .NE. 0 ) THEN
            WRITE(0,*) '*** NWCNEL: Got error from WCHFAC (3)'
            RETURN
         END IF
C
         IF( IFC2 .EQ. 4 ) THEN
C - special case of face 4 (of ELM2)
            IF( BIGLST(4,ELM2) .NE. 1 ) THEN
C - face 4 (of ELM2) already connected -error if not to ELM1
               IF( BIGLST(4,ELM2) .NE. ELM1*2+1 ) THEN
                  IERR = -6543
                  WRITE(0,*)
     :               '*** NWCNEL: ELEMENT 2 FACE ALREADY CONNECTED!'
                  WRITE(0,*) elm1,elm2,ifac,ifc2
                  WRITE(0,*) ied1,ied2,ied3
                  WRITE(0,*) (biglst(j,elm1),j=1,4)
                  nxt = biglst(nxtbig,elm1)
                  nxt2 = biglst(nxtbig,nxt)
                  WRITE(0,*) (biglst(j,nxt),j=1,4),
     :                       (biglst(j,nxt2),j=1,2)
                  if(elm2.gt.0) then
                    WRITE(0,*) (biglst(j,elm2),j=1,4)
                    nxt = biglst(nxtbig,elm2)
                    nxt2 = biglst(nxtbig,nxt)
                    WRITE(0,*) (biglst(j,nxt),j=1,4),
     :                         (biglst(j,nxt2),j=1,2)
                  endif
                  RETURN
               END IF
            ELSE
C - connect face 4 (of ELM2) to ELM1
               BIGLST(4,ELM2) = ELM1*2 + 1
            END IF
         ELSE IF( BIGLST(IFC2,ELM2) .NE. 0 ) THEN
C - face IFC2 (of ELM2) is already connected -error if not to ELM1
            IF( BIGLST(IFC2,ELM2) .NE. ELM1 ) THEN
               IERR = -6543
               WRITE(0,*)
     :             '*** NWCNEL: ELEMENT 2 FACE ALREADY CONNECTED!'
               WRITE(0,*) elm1,elm2,ifac,ifc2
               WRITE(0,*) ied1,ied2,ied3
               WRITE(0,*) (biglst(i,elm1),i=1,4)
               if(elm2.gt.0) then
                 WRITE(0,*) (biglst(i,elm2),i=1,4)
               endif
               RETURN
            END IF
         ELSE
C - connect face IFC2 (of ELM2) to ELM1
            BIGLST(IFC2,ELM2) = ELM1
         END IF
C
      END IF
C
      IF( IFAC .EQ. 4 ) THEN
C
C - deal with special case of face 4 (of ELM1)
C
         IF( ELM2 .LT. 0 ) THEN
c
c - what to do here????????????
C
c - ELM2 is a surface, but the face of ELM1 requested for connection is
c - number 4, which we can't do without swapping things around...
C
C - just have to swap with another face (3) and hope for the best...
C
C - check that face four (of ELM1) is empty
C
            IF(BIGLST(4,ELM1) .NE. 1) THEN
               IERR = -6543
               WRITE(0,*)
     :             '*** NWCNEL: ELEMENT 1 FACE ALREADY CONNECTED!'
               WRITE(0,*) elm1,elm2,ifac,ifc2
               WRITE(0,*) (biglst(i,elm1),i=1,4)
               if(elm2.gt.0) WRITE(0,*) (biglst(i,elm2),i=1,4)
               RETURN
            END IF
C
            IF( BIGLST(3,ELM1) .GE. 0 ) THEN
C - swap the faces (3 and 4 of ELM1) 
               I = BIGLST(3,ELM1)
               BIGLST(3,ELM1) = ELM2
               BIGLST(4,ELM1) = 1 + I*2
C - also swap edges 2/4 with 3/6 (of ELM1) to keep it consistent 
               NXT = BIGLST(NXTBIG,ELM1)
               I = BIGLST(2,NXT)
               BIGLST(2,NXT) = BIGLST(3,NXT)
               BIGLST(3,NXT) = I
               I = BIGLST(4,NXT)
               NXT2 = BIGLST(NXTBIG,NXT)
               BIGLST(4,NXT) = BIGLST(2,NXT2)
               BIGLST(2,NXT2) = I
            ELSE IF( BIGLST(2,ELM1) .GE. 0 ) THEN
C - swap the faces (2 and 4 of ELM1) 
               I = BIGLST(2,ELM1)
               BIGLST(2,ELM1) = ELM2
               BIGLST(4,ELM1) = 1 + I*2
C - also swap edges 1/5 with 2/6 (of ELM1) to keep it consistent 
               NXT = BIGLST(NXTBIG,ELM1)
               I = BIGLST(1,NXT)
               BIGLST(1,NXT) = BIGLST(2,NXT)
               BIGLST(2,NXT) = I
               NXT2 = BIGLST(NXTBIG,NXT)
               I = BIGLST(1,NXT2)
               BIGLST(1,NXT2) = BIGLST(2,NXT2)
               BIGLST(2,NXT2) = I
            ELSE IF( BIGLST(1,ELM1) .GE. 0 ) THEN
C - swap the faces (1 and 4 of ELM1) 
               I = BIGLST(1,ELM1)
               BIGLST(1,ELM1) = ELM2
               BIGLST(4,ELM1) = 1 + I*2
C - also swap edges 1/3 with 4/5 (of ELM1) to keep it consistent 
               NXT = BIGLST(NXTBIG,ELM1)
               I = BIGLST(1,NXT)
               BIGLST(1,NXT) = BIGLST(4,NXT)
               BIGLST(4,NXT) = I
               I = BIGLST(3,NXT)
               NXT2 = BIGLST(NXTBIG,NXT)
               BIGLST(3,NXT) = BIGLST(1,NXT2)
               BIGLST(1,NXT2) = I
            ELSE
               WRITE(0,*)
     :           '*** MWCNEL: CANNOT HAVE ALL SIDES ON SURFACES'
               WRITE(0,*) elm1,elm2,ifac,ifc2
               WRITE(0,*) (biglst(i,elm1),i=1,4)
               IERR = -98
               RETURN
            END IF
C
         ELSE IF( BIGLST(4,ELM1) .NE. 1 ) THEN
C - face 4 (of ELM1) already connected -error if not to ElM2
            IF( BIGLST(4,ELM1) .NE. ELM2*2+1 ) THEN
               IERR = -6543
               WRITE(0,*)
     :             '*** NWCNEL: ELEMENT 1 FACE ALREADY CONNECTED!'
               WRITE(0,*) elm1,elm2,ifac,ifc2
               WRITE(0,*) (biglst(i,elm1),i=1,4)
               if(elm2.gt.0) WRITE(0,*) (biglst(i,elm2),i=1,4)
               RETURN
            END IF
         ELSE
C - connect face 4 (of ELM1) to ELM2
            BIGLST(4,ELM1) = ELM2*2 + 1
         END IF
C
      ELSE IF( BIGLST(IFAC,ELM1) .NE. 0 ) THEN
C - face IFAC (of ELM1) already connected -error if not to ELM2
         IF( BIGLST(IFAC,ELM1) .NE. ELM2 ) THEN
            IERR = -6543
            WRITE(0,*) '*** NWCNEL: ELEMENT 1 FACE ALREADY CONNECTED!'
            WRITE(0,*) elm1,elm2,ifac,ifc2
            WRITE(0,*) (biglst(i,elm1),i=1,4)
            if(elm2.gt.0) WRITE(0,*) (biglst(i,elm2),i=1,4)
         END IF
      ELSE
C - connect face IFAC (of ELM1) to ELM2
         BIGLST(IFAC,ELM1) = ELM2
      END IF
C
      RETURN
      END
C
