C ======================================================================
      Subroutine saveMani(
C ======================================================================
     &      nP, nF, nE, nPv, nFv, nEv,
     &      XYP, IPF, IPE, IPV, IFV, IEV, lbE,
     &      ParCrv, iFnc,
     &      fName)
C ======================================================================
C  Routine saves the mesh in the single file fName.
C  The must have extension .ani.
C ======================================================================
C     Integer nP, nF, nE, nPv, nFv, nEv
      real  XYP(2, *)
      Integer IPE(3, *), IPF(4, *), IPV(*), IFV(*), IEV(*)
      Integer lbE(*)

      real  ParCrv(2, *)
      Integer iFnc(*)

      Character*(*) fName

c local variables
      Character*30  fNameExt

C ======================================================================
      i = 1
      Do while( fName(i:i+3) .NE. '.ani')
         i = i + 1
      End do

      fNameExt = fName(1:i) // 'ani'

      Write(*,'(/,A,A)') 'Saving mesh ', fNameExt

c ... count curved edges
      nC = 0
      Do n = 1, nF
         If(IPF(3, n).NE.0) nC = nC + 1
      End do


c ... save the mesh header
      kps = 10
      kfs = kps + nP + 2
      kes = kfs + nF + 2
      kcs = kes + nE + 2

      kfp = kcs + nC + 2
      kff = kfp + nPv + 2
      kfe = kff + nFv + 2

      kis = kfe + nEv + 2

      Open(10, file=fNameExt, status='UNKNOWN')
      Write(10, '(3(A,I7),A)') 'T points:        ', 
     &          nP, ' (lines ', kps, ' - ', kps + nP - 1, ')'
      Write(10, '(3(A,I7),A)') 'T edges:         ',
     &          nF, ' (lines ', kfs, ' - ', kfs + nF - 1, ')'
      Write(10, '(3(A,I7),A)') 'T elements:      ',
     &          nE, ' (lines ', kes, ' - ', kes + nE - 1, ')'

      If(nC.NE.0) Then
         Write(10, '(3(A,I7),A)') 'T curved edges:  ',
     &             nC, ' (lines ', kcs, ' - ', kcs + nC - 1, ')'
      Else
         Write(10, '(A,I7)') 'T curved faces:  ', nC
      End if

      If(nPv.NE.0) Then
         Write(10, '(3(A,I7),A)') 'T fixed points:  ',
     &             nPv, ' (lines ', kfp, ' - ', kfp + nPv - 1, ')'
      Else
         Write(10, '(A,I7)') 'T fixed points:  ', nPv
      End if

      If(nFv.NE.0) Then
         Write(10, '(3(A,I7),A)') 'T fixed edges:   ', 
     &             nFv, ' (lines ', kff, ' - ', kff + nFv - 1, ')'
      Else
         Write(10, '(A,I7)') 'T fixed edges:   ', nFv
      End if

      If(nEv.NE.0) Then
         Write(10, '(3(A,I7),A)') 'T fixed elements:', 
     &             nEv, ' (lines ', kfe, ' - ', kfe + nEv - 1, ')'
      Else
         Write(10, '(A,I7)') 'T fixed elements:', nEv 
      End if


c ... save the mesh
      Write(10,*)
      Write(10,*) nP, ' # of nodes'
      Do n = 1, nP
         Write(10,*) (XYP(j, n), j = 1, 2)
      End do

      Write(10,*)
      Write(10,*) nF, ' # of edges'
      k = 0
      Do n = 1, nF
         Write(10,*) (IPF(j, n), j = 1, 3), k, IPF(4, n)
      End do

      Write(10,*)
      Write(10,*) nE, ' # of elements'
      Do n = 1, nE
         Write(10,*) (IPE(j, n), j = 1, 3), lbE(n)
      End do


c ... saving curvilinear edges
      Write(10,*)
      Write(10,*) nC, ' # of curvilinear faces'
      Do n = 1, nC
         Write(10,*) (parCrv(j, n), j = 1, 2), iFnc(n)
      End do


c ... saving the fixed mesh points, faces and elements
      Write(10,*)
      Write(10,*) nPv, ' # number of fixed points'
      Do n = 1, nPv
         Write(10,*) IPV(n)
      End do


      Write(10,*)
      Write(10,*) nFv, ' # number of fixed faces'
      Do n = 1, nFv
         Write(10,*) IFV(n)
      End do


      Write(10,*)
      Write(10,*) nEv, ' # number of fixed elements'
      Do n = 1, nEv
         Write(10,*) IEV(n)
      End do

      Close(10)

      Return
      End



C ======================================================================
      Subroutine saveS(nP, Sol, fName)
C ======================================================================
C group (Q)
      real  Sol(*)
 
      Character*(*) fName
 
C group (Local variables)
      Character*30 fNameExt
 
C ======================================================================
      i = 1
      Do while( fName(i:i+3) .NE. '.sol')
         i = i + 1
      End do

      fNameExt = fName(1:i) // 'sol'
      Write(*,'(A,A)') 'Saving solution ', fNameExt

C save the solution associated to the mesh
      Open(10, file=fNameExt, status='UNKNOWN')
      Write(10,*) nP
 
      Write(10,*)
      Do n = 1, nP
         Write(10,*) Sol(n)
      End do
      Close(10)
 
      Return
      End



C ======================================================================
      Subroutine saveMgmv(nP, nE, XYP, IPE, fName)
C ======================================================================
C Routine saves mesh in the GMV file. 
C
C *** Remarks:
C        1. The size of the working memory is nE
C ======================================================================
C group (M)
      real  XYP(2, *)
      Integer IPE(3, *)
      Character*(*) fName

C group (Local variables)
      real       z
      Character*30 fNameExt

C ======================================================================
      z = 0D0

      fNameExt = fName // '.gmv'
      Open(10, file=fNameExt, status='UNKNOWN')

      Write(10, '(A)') 'gmvinput ascii'
      Write(10, *)
      Write(10, *) 'nodev ', nP
      Do n = 1, nP
         Write(10, *) (XYP(i, n), i = 1, 2), z
      End do

c ... save faces
      Write(10, '(A,I10)') 'cells ', nE
      Do n = 1, nE
         Write(10, *) 'general 1'
         Write(10, *) '3 ', (IPE(i, n), i = 1, 3)
      End do

      Write(10, '(A)') 'endgmv'

      Close(10)

      Return
      End



