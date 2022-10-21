      module mba3d_F2E
C
      use mba3d_auxSE
      use mba3d_auxSF
      use mba3d_auxSP
      use mba3d_error
      use mba3d_list
      use mba3d_makQ
      use mba3d_update
C
      contains
C
C ================================================================
      Subroutine F2E(
C ================================================================
c group (M)
     &           iwF, iwE,
     &           nE, MaxE, XYP, IPE,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           L1E, L2E, nL2, nStep,
     &           IHolE,
     &           status,
C group (Q)
     &           HesP, rQuality, detG, qE,
C group (S)
     &           lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &           flag)
C ================================================================
      include 'makS.fd'
      include 'color.fd'
      include 'operat.fd'
      include 'shutF.fd'
      include 'status.fd'
C ================================================================
C group (M)
      Integer IPE(5, *)
      Real*8  XYP(3, *)

      Real*8  hStar

      Integer L1E(2, *), L2E(*), nStep(4)
      Integer IHolE(*)

      Integer ICP(*), IEP(*)
      Integer IFE(4, *), IEE(4, *)

      Integer status

C group (Q)
      Real*8  HesP(6, *), rQuality
      Real*8  detG(*), qE(*)

C group (S)
      Integer iFu(*), iEu(*), IPFu(4, *), IPEu(5, *)
      Real*8  qEu(*)

C group (Flag)
      Logical flag

C ================================================================
C group (Local functions)

C group (Local variables)
      Integer iFs(MaxS), iEs(MaxS), IPFs(4, MaxS), IPEs(5, MaxS)
      Real*8  qEs(MaxS)

      Integer ip(5)
      Real*8  v
      Logical flagSHUT, flagTM, flag1, flag2, flag3

C ================================================================
      flag = .FALSE.

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1

      iF = IFE(iwF, iwE)

      i1 = iwF
      i2 = ip(i1 + 1)
      i3 = ip(i2 + 1)

      iP1 = IPE(i1, iwE)
      iP2 = IPE(i2, iwE)
      iP3 = IPE(i3, iwE)

      iE1 = iwE
      iE2 = IEE(iwF, iE1)


C ... checking the case when face -> edge operation is impossible
      If(iE2.EQ.0) Goto 1000
      If(iF.NE.0)  Goto 1000

      i4 = ip(i3 + 1)
      iPa = IPE(i4, iwE)

      Do j1 = 1, 4
         j2 = ip(j1 + 1)
         j3 = ip(j2 + 1)
         If(check33(i1, i2, i3, j1, j2, j3)) Then
            j4 = ip(j3 + 1)
            iPb = IPE(j4, iE2)
         End if
      End do


C ... checking the case when face -> edge operation is impossible
      flagSHUT = lOpenTriag
      Call shutF(XYP(1, iPa), XYP(1, iPb),
     &           XYP(1, iP1), XYP(1, iP2), XYP(1, iP3),
     &           flagSHUT)
      If(.NOT.flagSHUT) Goto 1000


C ... number of inverted elements
      flagTM = ifXnode(status, ANITangledMesh)
      if(flagTM) Then
         nBad = 0
         If(qEs(iE1).LE.0D0) nBad = nBad + 1
         If(qEs(iE2).LE.0D0) nBad = nBad + 1
      End if


C ... saving the initial superelement structure
      Call copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)



C ... making a virtual evaluation of the quality
      Call findSE(lE, iEs, iE1, nE1)
      Call calQE(
     &     HesP(1, iPa), detG(iPa), XYP(1, iPa),
     &     HesP(1, iPb), detG(iPb), XYP(1, iPb),
     &     HesP(1, iP1), detG(iP1), XYP(1, iP1),
     &     HesP(1, iP3), detG(iP3), XYP(1, iP3),
     &     hStar, qEs(nE1), v)

      If(qEs(nE1).LE.rQuality) Goto 1000

      Call findSE(lE, iEs, iE2, nE2)
      Call calQE(
     &     HesP(1, iPa), detG(iPa), XYP(1, iPa),
     &     HesP(1, iPb), detG(iPb), XYP(1, iPb),
     &     HesP(1, iP2), detG(iP2), XYP(1, iP2),
     &     HesP(1, iP3), detG(iP3), XYP(1, iP3),
     &     hStar, qEs(nE2), v)

      If(qEs(nE2).LE.rQuality) Goto 1000

      lE = lE + 1
      If(lE.GT.MaxS) Goto 9000

      Call calQE(
     &     HesP(1, iPa), detG(iPa), XYP(1, iPa),
     &     HesP(1, iPb), detG(iPb), XYP(1, iPb),
     &     HesP(1, iP1), detG(iP1), XYP(1, iP1),
     &     HesP(1, iP2), detG(iP2), XYP(1, iP2),
     &     hStar, qEs(lE), v)

      If(qEs(lE).LE.rQuality) Goto 1000


C ... checking for boundary elements
      If(ifXnode(status, ANIForbidBoundaryElements) .AND.
     &  ifXnode(ICP(iPa), jBnode) .AND.
     &  ifXnode(ICP(iPb), jBnode)) Then
        flag1 = ifXnode(ICP(iP1), jBnode)  
        flag2 = ifXnode(ICP(iP2), jBnode) 
        If(flag1 .AND. flag2) Goto 1000

        flag3 = ifXnode(ICP(iP3), jBnode)
        If(flag1 .AND. flag3) Goto 1000
        If(flag2 .AND. flag3) Goto 1000
      End if


C ... creating new elements
      ide = IPEs(5, nE1)

      IPEs(1, nE1) = iPa
      IPEs(2, nE1) = iPb
      IPEs(3, nE1) = iP1
      IPEs(4, nE1) = iP3
      IPEs(5, nE1) = ide

      IPEs(1, nE2) = iPa
      IPEs(2, nE2) = iPb
      IPEs(3, nE2) = iP2
      IPEs(4, nE2) = iP3
      IPEs(5, nE2) = ide

      IPEs(1, lE)  = iPa
      IPEs(2, lE)  = iPb
      IPEs(3, lE)  = iP1
      IPEs(4, lE)  = iP2
      IPEs(5, lE)  = ide


C ... checking for inverted elements
      If(flagTM) Then
         Call updQb(nE1, lE, iEs, XYP, IPEs, qEs)
         Call updQb(nE2, lE, iEs, XYP, IPEs, qEs)
         Call updQb(lE,  lE, iEs, XYP, IPEs, qEs)

         mBad = 0
         If(qEs(nE1).LE.0D0) mBad = mBad + 1
         If(qEs(nE2).LE.0D0) mBad = mBad + 1
         If(qEs(lE) .LE.0D0) mBad = mBad + 1

         If(mBad.GE.nBad) Goto 1000
      End if


C ... updating the grid
      flag = .TRUE.

      Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(nE1), qEs(nE1))
      Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(nE2), qEs(nE2))

      Call eleDel(iEs(nE1), IPE, IEE)
      Call eleDel(iEs(nE2), IPE, IEE)

      Call eleAdd(nE, MaxE, IHolE)
      Call lstAdd(nE, L1E, nL2, L2E, nStep, IHolE,
     &            qE, qEs(lE), iEs(lE))
      Call eleDel(iEs(lE), IPE, IEE)

      Call eleUpd(nE1, IEP, IPE, IFE,  IEE,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs)
      Call eleUpd(nE2, IEP, IPE, IFE,  IEE,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs)
      Call eleUpd(lE,  IEP, IPE, IFE,  IEE,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs)

 1000 Return
 9000 Call errMes(1007, 'F2E', 'local parameter MaxS is small')
      End Subroutine F2E
C
      End Module mba3d_F2E
