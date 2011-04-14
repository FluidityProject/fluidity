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
      SUBROUTINE MSHOUT( ERRMES, FIELDS, NFIELD, FLDNMS, FLDNAM, NFREE,
     :                   NODX, NODY, NODZ, NNOD, NGEOM,
     :                   ENLBAS, ENLIST, NELM, SZENLS,
     :                   MSHNAM, LENFIL, COUNT, NOASK,
     :                   OUTMSH, OUTFLD, OUTERR, OUTALL,
     :                   ELMREG, PUTREG, ELMAT, PUTMAT, ELMSOR, PUTSOR,
     :                   NEWNOD, PUTNWN, DOMNOD, PUTDOM )
C-----------------------------------------------------------------------
C
C - This subroutine dumps the mesh to a file, the fields to another file,
C - and the error measure to another.
C
C - The filename is sent as MSHNAM (with LENFIL the length <110 ).
C - Extensions are added to this filename for each of the files dumped:
C - the first extension is integer COUNT converted into a string of
C - three characters (with leading zeros - i.e. five is 005).
C - The next extension is '.mesh' for the mesh file, '.field' for the fields
C - and '.error' for the error measure, unless OUTALL is true, in which
C - case the mesh and field files are combined into a file with '.all'
C - extension.
C
C - For example, if MSHNAM is "output-file" (with LENFIL = 12), and
C - COUNT is seven, then the files will be "output-file.007.mesh",
C - "output-file.007.field", and "output-file.007.error".
C - If OUTALL is true then the '.mesh' and '.field' files will be replaced
C - by a single '.all' file.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER NNOD, NELM, SZENLS, ELMREG(NELM), ENLBAS(NELM+1),
     :        ENLIST(SZENLS), LENFIL, COUNT, NFIELD, NFREE, NGEOM,
     :        ELMAT(NELM), ELMSOR(NELM), NEWNOD(NNOD), DOMNOD(NNOD)
C
      REAL ERRMES(NNOD), NODX(NNOD), NODY(NNOD), NODZ(NNOD),
     :     FIELDS(NFREE,NNOD,NFIELD)
C
      CHARACTER*128 MSHNAM
      CHARACTER*20 FLDNAM(NFIELD)
C
      LOGICAL NOASK, OUTMSH, OUTFLD, OUTERR, OUTALL, PUTMAT, PUTREG,
     :        PUTNWN, PUTDOM, PUTSOR, FLDNMS
C-----------------------------------------------------------------------
C
C - ERRMES is the array containing the nodwise error measure
C - FIELDS is the array containing the NFIELD fields with NFREE freedoms
C - NFIELD is the number of fields
C - FLDNMS is used to switch on writing of the field names
C - FLDNAM contains a 20 character text description for each field
C - NFREE is the number of freedoms for each field
C - NODX,NODY,NODZ are the node coordinates
C - NNOD is the number of nodes
C - NGEOM is the geometry ( 1,2,3: 1d  4,5: 2d  6: 3d )
C - ELMREG is the element region/material id for each element
C - ENLBAS is the element node list base pointer
C - ENLIST is the element node list
C - NELM is the number of elements
C - SZENLS is the size of the element node list
C - MSHNAM is a string with the base filename (to which extensions are added)
C - LENFIL is the length of the base filename (not including extensions)
C - COUNT is an integer which is used as an extension to the base filename
C - NOASK is a logical which inhibits the 'Do you want to...?' questions
C - OUTMSH is a logical which switches dumping of the mesh file on/off
C - OUTFLD is a logical which switches dumping of the field file on/off
C - OUTERR is a logical which switches dumping of the error file on/off
C - OUTALL is a logical which overrides OUTMSH and OUTFLD to create a
C          combined mesh/field file with a great deal more flexibility
C          than both the mesh and field files, since they make use of...
C - PUTREG is sent TRUE if ELMREG is to be dumped in the '.all' file
C - PUTMAT is sent TRUE if ELMAT  is to be dumped in the '.all' file
C - PUTSOR is sent TRUE if ELMSOR is to be dumped in the '.all' file
C - PUTNWN is sent TRUE if NEWNOD is to be dumped in the '.all' file
C - PUTDOM is sent TRUE if DOMNOD is to be dumped in the '.all' file
C
C-----------------------------------------------------------------------
      INTEGER I, OUTFIL, NCNSTR, IELM, INOD, NREG, ILN, NLN,
     :        IBAS, NODS(20), J, K, ERRFIL, FLDFIL, IFLAGS, IFLEN
C
      LOGICAL DMPMSH, DMPERR, DMPFLD, RSPNSE
C
      CHARACTER CNTCHR*3, NUMLST*10
C
      OUTFIL = 16
      ERRFIL = 16
      FLDFIL = 16
C
      IFLEN  = 20
C
      NUMLST = '0123456789'
C
      IBAS = COUNT
      I = MOD(IBAS,10)
      IBAS = ( IBAS - I )/10
      J = MOD(IBAS,10)
      IBAS = ( IBAS - J )/10
      K = MOD(IBAS,10)
      IBAS = ( IBAS - K )/10
      I = I + 1
      J = J + 1
      K = K + 1
      CNTCHR = NUMLST(K:K)//NUMLST(J:J)//NUMLST(I:I)
C
      IF( .NOT. OUTALL ) THEN
C
       DMPMSH = OUTMSH
       IF( DMPMSH .AND. .NOT. NOASK ) THEN
          DMPMSH = RSPNSE( 'Do you want to dump the mesh?' ,29 )
       END IF
C
       IF( DMPMSH ) THEN
C
          NREG = -100000000
          DO IELM = 1, NELM
             NREG = MAX( NREG, ELMREG(IELM) )
          END DO
C 
          MSHNAM = MSHNAM(1:LENFIL)//'.'//CNTCHR//'.mesh'
          OPEN( UNIT = OUTFIL, FILE = MSHNAM(1:LENFIL+9),
     :          STATUS = 'UNKNOWN', FORM = 'UNFORMATTED' )
C     :          STATUS = 'UNKNOWN', FORM = 'FORMATTED' )
C 
          NCNSTR = 0
C 
C          WRITE( OUTFIL, 990 ) NELM, NNOD, NCNSTR
          WRITE( OUTFIL ) 1, NREG, 0, NELM, NNOD, NCNSTR, 0, 0, 0,
     :                    1, 1, 0
C 
          DO IELM = 1, NELM
             IBAS = ENLBAS(IELM)
             NLN  = ENLBAS(IELM+1) - IBAS
             DO ILN = 1, NLN
                NODS(ILN) = ENLIST(IBAS+ILN)
             END DO
C             WRITE( OUTFIL, 991 ) IELM,1,1,ELMREG(IELM),0,NODS
             WRITE( OUTFIL ) IELM,1,1,ELMREG(IELM),0,NLN,
     :                       (NODS(ILN),ILN=1,NLN)
          END DO
C 
          DO INOD = 1, NNOD
C             WRITE( OUTFIL, 992 ) INOD,NODX(INOD),NODY(INOD),NODZ(INOD),
C     :                            INOD, 1, 1, 1
             WRITE( OUTFIL ) NODX(INOD),NODY(INOD),NODZ(INOD),
     :                       INOD, 1, 1, 1
          END DO
C
          PRINT*,'MESH WRITTEN TO FILE ',MSHNAM(1:LENFIL+9)
C
          CLOSE( OUTFIL )
C
       END IF
C
       DMPFLD = OUTFLD
       IF( DMPFLD .AND. .NOT. NOASK ) THEN
          DMPFLD = RSPNSE( 'Do you want to dump the fields?' ,30 )
       END IF
C
       IF( DMPFLD ) THEN
C
          MSHNAM = MSHNAM(1:LENFIL)//'.'//CNTCHR//'.field'
          OPEN( UNIT = FLDFIL, FILE = MSHNAM(1:LENFIL+10),
     :          STATUS = 'UNKNOWN', FORM = 'UNFORMATTED' )
C
          WRITE( FLDFIL ) NNOD, NFREE, NFIELD, 0.0, 0, 0.0, 1.0
C
          WRITE( FLDFIL ) NODX
          WRITE( FLDFIL ) NODY
          WRITE( FLDFIL ) NODZ
C
          WRITE( FLDFIL ) FIELDS
C
          CLOSE( FLDFIL )
C
          PRINT*,'FIELDS WRITTEN TO FILE ',MSHNAM(1:LENFIL+10)
C
       END IF
C
       DMPERR = OUTERR
       IF( DMPERR .AND. .NOT. NOASK ) THEN
          DMPERR = RSPNSE( 'Do you want to dump the error?' ,30 )
       END IF
C
       IF( DMPERR ) THEN
C
          MSHNAM = MSHNAM(1:LENFIL)//'.'//CNTCHR//'.error'
          OPEN( UNIT = ERRFIL, FILE = MSHNAM(1:LENFIL+10),
     :          STATUS = 'UNKNOWN', FORM = 'UNFORMATTED' )
C
          WRITE( ERRFIL ) NNOD, 1, 1, 0.0, 0, 0.0, 1.0
C
          WRITE( ERRFIL ) NODX
          WRITE( ERRFIL ) NODY
          WRITE( ERRFIL ) NODZ
          WRITE( ERRFIL ) ERRMES
C
          CLOSE( ERRFIL )
C
          PRINT*,'ERROR MEASURE WRITTEN TO FILE ',MSHNAM(1:LENFIL+10)
C
       END IF
C
C
      ELSE
C
       IF( NOASK ) THEN
          DMPFLD = .TRUE.
       ELSE IF( OUTERR ) THEN
          DMPFLD = RSPNSE( 'Do you want to dump everything?' ,31 )
       ELSE
          DMPFLD = RSPNSE( 'Do you want to dump the mesh/field?' ,35 )
       END IF
C
       IF( DMPFLD ) THEN
C
          IFLAGS = 131072 + 65536 + 4096
          IF( FLDNMS ) IFLAGS = IFLAGS + 8192
          IF( PUTSOR ) IFLAGS = IFLAGS + 128
          IF( PUTREG ) IFLAGS = IFLAGS + 64
          IF( PUTMAT ) IFLAGS = IFLAGS + 32
          IF( PUTNWN ) IFLAGS = IFLAGS + 2
          IF( PUTDOM ) IFLAGS = IFLAGS + 1 
C
          MSHNAM = MSHNAM(1:LENFIL)//'.'//CNTCHR//'.all'
          OPEN( UNIT = FLDFIL, FILE = MSHNAM(1:LENFIL+8),
     :          STATUS = 'UNKNOWN', FORM = 'UNFORMATTED' )
C
          WRITE( FLDFIL ) IFLAGS, NGEOM, NFREE, 0, NNOD, NFREE, NELM,
     :                    0, 0, 0, NFIELD, 0, 0.0, 0.0, 1.0
C
          WRITE( FLDFIL ) ENLBAS
C
          WRITE( FLDFIL ) ENLIST
C
          IF( FLDNMS ) THEN
             DO I = 1, NFIELD
                WRITE( FLDFIL ) IFLEN, FLDNAM(I)
             END DO
          END IF
C
          IF( PUTSOR ) WRITE( FLDFIL ) ELMSOR
          IF( PUTREG ) WRITE( FLDFIL ) ELMREG
          IF( PUTMAT ) WRITE( FLDFIL ) ELMAT
          IF( PUTNWN ) WRITE( FLDFIL ) NEWNOD
          IF( PUTDOM ) WRITE( FLDFIL ) DOMNOD
C
          WRITE( FLDFIL ) NODX
          IF( NGEOM .GT. 3 ) THEN
             WRITE( FLDFIL ) NODY
             IF( NGEOM .GT. 5 ) WRITE( FLDFIL ) NODZ
          END IF
C
          WRITE( FLDFIL ) FIELDS
C
          IF( OUTERR ) THEN
C
           WRITE( FLDFIL ) 0, NGEOM, 1, 0, NNOD, 1, NELM, 0,
     :                     0, 0, 1, 0, 0.0, 0.0, 1.0
C
           WRITE( FLDFIL ) ERRMES
C
          END IF
C
          CLOSE( FLDFIL )
C
          IF( OUTERR ) THEN
             PRINT*,'MESH, FIELDS, ERROR WRITTEN TO FILE ',
     :                                     MSHNAM(1:LENFIL+10)
          ELSE
             PRINT*,'MESH AND FIELDS WRITTEN TO FILE ',
     :                                     MSHNAM(1:LENFIL+10)
          END IF
C
       END IF
C
      END IF
      MSHNAM = MSHNAM(1:LENFIL)
C
 990  FORMAT( 'SIZES                    ', 3I5 )
 991  FORMAT( 'ELEM ', 9I5 )
 992  FORMAT( 'NODE ', I5, 3(1PE15.7), 4I5 )
      RETURN
      END
C
