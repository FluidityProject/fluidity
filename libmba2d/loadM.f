C ================================================================
      Subroutine loadMone(
C ================================================================
c group (M)
     &      nP,  MaxP,  nF,  MaxF,  nE,  MaxE,
     &      nPv, MaxPV, nFv, MaxFV, nEv, MaxEV,
     &      XYP, IPF, IPE, IPV, IFV, IEV, lbE,
     &      ParCrv, iFnc,
     &      fName)
C ==========================================================
C Routines read the input mesh from single file fName.sol.
C ==========================================================
C group (M)
C     Integer MaxP, MaxF, MaxE
      real  XYP(2, *)
      Integer IPE(3, *), IPF(4, *), IPV(*), IFV(*), IEV(*) 
      Integer lbE(*)

      real  ParCrv(2, *)
      Integer iFnc(*)

      Character*(*) fName

c group (Local variables)
      Logical flagP, flagF, flagE, flagC, flagPV, flagFV, flagEV
      Character*30  fNameExt

C ==========================================================
      i = 1
      Do while( fName(i:i+3) .NE. '.ani')
         i = i + 1
      End do
      fNameExt = fName(1:i) // 'ani'

      Write(*,'(A,A)') 'Loading mesh ', fNameExt

      Open(10, file=fNameExt, status='OLD', ERR=1000)
      Read(10,*) flagP 
      Read(10,*) flagF
      Read(10,*) flagE
      Read(10,*) flagC
      Read(10,*) flagPV
      Read(10,*) flagFV
      Read(10,*) flagEV

      Read(10,*)
      Read(10,*) nP
      If(nP.GT.MaxP) Call errMes(1003, 'loadMone',
     &                          'local parameter MaxP is small')
      Do n = 1, nP
         Read(10,*) (XYP(j, n), j = 1, 2)
      End do


c ... reading the boundary edges: 
c     indices of the edge ends,                  1 2
c     number of respective line of curved edges  3
c     dummy integer                              4
c     label of the edge                          5
      Read(10,*)
      Read(10,*) nF
      If(nF.GT.MaxF) Call errMes(1004, 'loadMone',
     &                          'local parameter MaxF is small')
      Do n = 1, nF
         Read(10,*) (IPF(j, n), j = 1, 3), k, IPF(4, n)
      End do


      Read(10,*)
      Read(10,*) nE
      If(nE.GT.MaxE) Call errMes(1006, 'loadMone',
     &                          'local parameter MaxE is small')
      Do n = 1, nE
         Read(10,*) (IPE(j, n), j = 1, 3), lbE(n)
      End do


c ... reading parametrization of the curved boundary edges:
c     parametrization of the edge ends                      1 2
c     reference to the function of the parametrization      3
      Read(10,*)
      Read(10,*) nCrv
      If(nCrv.GT.MaxF) Call errMes(1004, 'loadMone',
     &                            'local parameter MaxF is small')
      Do n = 1, nCrv
         Read(10,*) (parCrv(j, n), j = 1, 2), iFnc(n)
      End do


c ... reading the fixed mesh points, faces and elements
      Read(10,*)
      Read(10,*) nPv
      If(nPv.GT.MaxPv) Call errMes(1008, 'loadM', 
     &                            'local parameter MaxPV is small')
      Do n = 1, nPv
         Read(10,*) IPV(n)
      End do


      Read(10,*)
      Read(10,*) nFv
      If(nFv.GT.MaxFV) Call errMes(1008, 'loadM', 
     &                            'local parameter MaxFV is small')
      Do n = 1, nFv
         Read(10,*) IFV(n)
      End do


      Read(10,*)
      Read(10,*) nEv
      If(nEv.GT.MaxEV) Call errMes(1008, 'loadM', 
     &                            'local parameter MaxEV is small')
      Do n = 1, nEv
         Read(10,*) IEV(n)
      End do

      Close(10)

      Return

 1000 Continue
      Call errMes(4001, 'loadMone', 'Input file name is wrong ')
      Return
      End



C ==========================================================
      Subroutine loadS(nP, Sol, fName)
C ==========================================================
C Routine reads the solution from file fName.sol.
C ==========================================================
      real        Sol(*)
      Character*(*) fName

C group (Local variables)
      Character*30  fNameExt

C ==========================================================
      i = 1
      Do while( fName(i:i+3) .NE. '.sol')
         i = i + 1
      End do
      fNameExt = fName(1:i) // 'sol'

      Write(*,'(A,A)') 'Loading solution ', fNameExt

      Open(10, file=fNameExt, status='UNKNOWN')
      Read(10,*) nPw
      If(nPw.NE.nP) Call errMes(4003, 'loadS', 
     &                   'mesh & data files are incompartible')

      Read(10,*)
      Do n = 1, nP
         Read(10,*) Sol(n)
      End do
      Close(10)

      Return
      End



C ==========================================================
      Subroutine loadM(
C ==========================================================
c group (M)
     &      nP,  MaxP,  nF,  MaxF,  nE,  MaxE,
     &      nPv, MaxPV, nFv, MaxFV, nEv, MaxEV,
     &      XYP, IPF, IPE, IPV, IFV, IEV, lbE,
     &      ParCrv, iFnc,
     &      fName)
C ==========================================================
C Routines read the input mesh from files fName.* 
C 
C *** Remark: no extension is required.
C ==========================================================
C group (M)
C     Integer MaxP, MaxF, MaxE
      real  XYP(2, *)
      Integer IPE(3, *), IPF(4, *), IPV(*), IFV(*), IEV(*) 
      Integer lbE(*)

      real  ParCrv(2, *)
      Integer iFnc(*)

      Character*(*) fName

C group (Local variables)
      Character*30 fNameExt

C ==========================================================
      Write(*,'(A,A)') 'Loading mesh ', fName

c ... reading coordinates of the nodes
      fNameExt = fName // '.vrt'
      Open(10, file=fNameExt, status='OLD', ERR=1000)
      Read(10,*) nP
      If(nP.GT.MaxP) Call errMes(1003, 'loadM', 
     &                          'local parameter MaxP is small')

      Read(10,*)
      Do n = 1, nP
         Read(10,*) (XYP(j, n), j = 1, 2)
      End do
      Close(10)


c ... reading the connectivity table
c     indices of the vertices,   1 2 3
c     label of the triangle      4
      fNameExt = fName // '.tri'
      Open(10, file=fNameExt, status='OLD', ERR=1000)
      Read(10,*) nE
      If(nE.GT.MaxE) Call errMes(1006, 'loadM', 
     &                          'local parameter MaxE is small')

      Read(10,*)
      Do n = 1, nE
         Read(10,*) (IPE(j, n), j = 1, 3), lbE(n)
      End do
      Close(10)


c ... reading the boundary edges: 
c     indices of the edge ends,                  1 2
c     number of respective line in fname//.crv,  3
c     dummy integer                              4
c     label of the edge                          5
      nF = 0
      fNameExt = fName // '.bnd'
      Open(10, file=fNameExt, status='OLD', ERR=100)
      Read(10,*) nF
      If(nF.GT.MaxF) Call errMes(1004, 'loadM', 
     &                          'local parameter MaxF is small')

      Read(10,*)
      Do n = 1, nF
         Read(10,*) (IPF(j, n), j = 1, 3), k, IPF(4, n)
      End do
      Close(10)


c ... reading parametrization of the curved boundary edges:
c     parametrization of the edge ends                      1 2
c     reference to the function of the parametrization      3
 100  nCrv = 0 
      fNameExt = fName // '.crv'
      Open(10, file=fNameExt, status='OLD', ERR=200)
      Read(10,*) nCrv

      Read(10,*)
      Do n = 1, nCrv
         Read(10,*) (parCrv(j, n), j = 1, 2), iFnc(n)
      End do
      Close(10)


c ... reading fixed (steady) nodes of the mesh
 200  nPv = 0
      fNameExt = fName // '.fix'
      Open(10, file=fNameExt, status='OLD', ERR=300)
      Read(10,*) nPv
      If(nPv.GT.MaxPv) Call errMes(1008, 'loadM', 
     &                            'local parameter MaxPV is small')

      Read(10,*)
      Do n = 1, nPv
         Read(10,*) IPV(n)
      End do
      Close(10)


c ... reading fixed (steady) nodes of the mesh
 300  nFv = 0
      fNameExt = fName // '.ffv'
      Open(10, file=fNameExt, status='OLD', ERR=400)
      Read(10,*) nFv
      If(nFv.GT.MaxFV) Call errMes(1008, 'loadM', 
     &                            'local parameter MaxFV is small')

      Read(10,*)
      Do n = 1, nFv
         Read(10,*) IFV(n)
      End do
      Close(10)


c ... reading fixed (steady) edges of the mesh
 400  nEv = 0
      fNameExt = fName // '.fev'
      Open(10, file=fNameExt, status='OLD', ERR=9000)
      Read(10,*) nEv
      If(nEv.GT.MaxEV) Call errMes(1008, 'loadM', 
     &                            'local parameter MaxEV is small')

      Read(10,*)
      Do n = 1, nEv
         Read(10,*) IEV(n)
      End do
      Close(10)

 1000 Call errMes(4001, 'loadM', 'missing files')

 9000 Return
      End


C ================================================================
      Subroutine loadMgmv(
C ================================================================
c group (M)
     &      nP,  MaxP,  nF,  MaxF,  nE,  MaxE,
     &      nPv, MaxPV, nFv, MaxFV, nEv, MaxEV,
     &      XYP, IPF, IPE, IPV, IFV, IEV, lbE,
     &      ParCrv, iFnc,
     &      fName)
C ==========================================================
C Routines read the input mesh from single file fName.gmv.
C ==========================================================
C group (M)
C     Integer MaxP, MaxF, MaxE
      real  XYP(2, *)
      Integer IPE(3, *), IPF(4, *), IPV(*), IFV(*), IEV(*) 
      Integer lbE(*)

      real  ParCrv(2, *)
      Integer iFnc(*)

      Character*(*) fName

c group (Local variables)
      Character*30  fNameExt, keyword

C ==========================================================
      fNameExt = fName // '.gmv'

      Write(*,'(A,A)') 'Loading mesh ', fNameExt

      Open(10, file=fNameExt, status='OLD', ERR=1000)
      Read(10,*) 
      Read(10,*) 
      Read(10,*) 

      Read(10,*) keyword, nP
      If(nP.GT.MaxP) Call errMes(1003, 'loadMgmv',
     &                          'local parameter MaxP is small')
      Do n = 1, nP
         Read(10,*) (XYP(j, n), j = 1, 2)
      End do


c ... reading the boundary edges: 
      nF = 0
      If(nF.GT.MaxF) Call errMes(1004, 'loadMgmv',
     &                          'local parameter MaxF is small')


c ... reading the elements: 
      Read(10,*) keyword, nE
      If(nE.GT.MaxE) Call errMes(1006, 'loadMgmv',
     &                          'local parameter MaxE is small')
      Do n = 1, nE
         Read(10,*) keyword, m, (IPE(j, n), j = 1, 3)
      End do


c ... reading parametrization of the curved boundary edges:
      nCrv = 0
      If(nCrv.GT.MaxF) Call errMes(1004, 'loadMgmv',
     &                            'local parameter MaxF is small')


c ... reading the fixed mesh points, faces and elements
      nPv = 0 
      nFv = 0 
      nEv = 0 

      Close(10)

      Return

 1000 Continue
      Call errMes(4001, 'loadMgmv', 'Input file name is wrong ')

      Return
      End

