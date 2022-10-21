      Module mba3d_deletP
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
      Subroutine deletP(
C ================================================================
c group (M)
     &           iwP, iwE,
     &           nP, nF, MaxF, nE, MaxE,
     &           XYP, IPF, IPE,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           L1E, L2E, nL2, nStep,
     &           IHolP, IHolF, IHolE,
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
      include 'magic.fd'
      include 'status.fd'
C ================================================================
C group (M)
      Integer IPF(4, *), IPE(5, *)
      Real*8  XYP(3, *)

      Real*8  hStar

      Integer L1E(2, *), L2E(*), nStep(4)
      Integer IHolP(*), IHolF(*), IHolE(*)

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
      Logical flag, flagTM

C ================================================================
C group (Local functions)

C group (Local variables)
      Integer iFs(MaxS), iEs(MaxS)
      Integer IPFs(4, MaxS), IPEs(5, MaxS)
      Real*8  qEs(MaxS)

      Integer ip(5), iref(4), iPt(3), nPt(3), iEd(1)

      Integer iP1s(MaxS), iP2s(MaxS), ICP2s(MaxS)
      Integer iDPs(2, MaxS), iNPs(2, MaxS)
      Integer iDFs(MaxS), iDEs(MaxS), iOs(MaxS)
      Integer IPSs(3, MaxS), IESs(MaxS)

      Real*8  oldVolume, newVolume, v
      Logical flagBNDs, flagOrient, flagFACE, flagFACE2, flagEDGE
      Logical l1

C ================================================================
      flag = .FALSE.

      Call copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1


      iref(1) = 1
      iref(2) = 2
      iref(3) = 3
      iref(4) = 1


      iPd = IPE(iwP, iwE)


C ... checking the case when deleting is impossible
      If(ifXnode(ICP(iPd), jVnode)) Goto 1000
      If(ifXnode(ICP(iPd), jTnode)) Goto 1000


C ... analyzing the superelement entries
      Call makSF(iPd, lF, iFs, lE, iEs, IPEs, IPF, IFE, MaxS,
     &           lP1, iP1s, lP2, iP2s, ICP2s,
     &           lS, IPSs, IESs, lDF, iDFs, lDE, iDEs, 
     &           flagFACE, flagEDGE)
C ... intersection of two interfaces at the point has more than 3 edges
      If(flagEDGE .AND. lP1.GE.3) Goto 1000

      If(lP2 + 1.GE.nP) Goto 1000


      lDP = lP2
      Do n = 1, lDP
         iDPs(1, n) = iPd
         iDPs(2, n) = iP2s(n)
      End do


      lFold = lF

      Do n = 1, lDF
         nFt = iDFs(n)
         iFs(nFt) = -iFs(nFt)
      End do


      Do n = 1, lDE
         nEt = iDEs(n)
         iEs(nEt) = -iEs(nEt)
      End do


C ... checking the number of inverted elements
      flagTM = ifXnode(status, ANITangledMesh)
      nBad = 0
      if(flagTM) Then
         Do n = 1, lDE
            nEt = iDEs(n)
            If(qEs(nEt).LE.0D0) nBad = nBad + 1
         End do
      End if

      flagTM = flagTM .AND. nBad.GT.0


C ... making a virtual evaluation of the quality
      oldVolume = 0D0
      Do 20 n = 1, lDE
         nEd = iDEs(n)

         iP1 = IPEs(1, nEd)
         iP2 = IPEs(2, nEd)
         iP3 = IPEs(3, nEd)
         iP4 = IPEs(4, nEd)

         v = calVol(XYP(1, iP1), XYP(1, iP2),
     &              XYP(1, iP3), XYP(1, iP4))

         oldVolume = oldVolume + dabs(v)
         iOs(n) = dsign(1D0, v)
 20   Continue


      Do 80 n = 1, lP1
         iPa = iP1s(n)

c  ...  checking for the orientation
         Do 40 k = 1, lS
            iP1 = IPSs(1, k)
            iP2 = IPSs(2, k)
            iP3 = IPSs(3, k)

            iEd(1) = iabs(iEs(iDEs(k)))
            If(check13(iPa, iP1, iP2, iP3)) Then
c  ...  checking in details this face
               Do 30 m = 1, lS
                  If(m.EQ.k) Goto 30

                  icnt = 0
                  Do i = 1, 3
                     If(IPSs(i, k).NE.iPa) Then
                        If(check13(IPSs(i, k), IPSs(1, m), IPSs(2, m),
     &                             IPSs(3, m))) icnt = icnt + 1
                     End if
                  End do

                  If(icnt.EQ.2) Goto 40
 30            Continue

               Goto 80
            End if

            Call chkSO(iPd, XYP(1, iPa), XYP, IPE,
     &                 1, iEd, iOs(k), flagOrient)
            If(.NOT.flagOrient) Goto 80
 40      Continue


c  ...  trying for the quality
         lEadd = 0D0 
         newVolume = 0D0
         Do 50 k = 1, lS
            iP1 = iPSs(1, k)
            iP2 = iPSs(2, k)
            iP3 = iPSs(3, k)

            If(check13(iPa, iP1, iP2, iP3)) Goto 50

            lEadd = lEadd + 1
            lEnew = lE + lEadd
            If(lEnew.GT.MaxS) Goto 9000
            iEs(lEnew) = 0

            Call calQE(
     &           HesP(1, iPa), detG(iPa), XYP(1, iPa),
     &           HesP(1, iP1), detG(iP1), XYP(1, iP1),
     &           HesP(1, iP2), detG(iP2), XYP(1, iP2),
     &           HesP(1, iP3), detG(iP3), XYP(1, iP3),
     &           hStar, qEs(lEnew), v)

            If(qEs(lEnew).LE.rQuality) Goto 80

            IPEs(1, lEnew) = iPa
            IPEs(2, lEnew) = iP1
            IPEs(3, lEnew) = iP2
            IPEs(4, lEnew) = iP3
            IPEs(5, lEnew) = IPEs(5, IESs(k))

            newVolume = newVolume + dabs(v)
 50      Continue

         lEold = lE
         lE = lEnew
         Goto 100
 80   Continue


      Goto 1000


 100  Continue
C ... checking for the global volume
      If(.NOT.flagTM .AND.
     &   dabs(oldVolume - newVolume).GT.volPREC*oldVolume) Then
C        Write(*, 5000) oldVolume, newVolume
         Goto 1000
      End if


C ... checking for surrounding points
      lNP = 0
      Do n = 1, lP2
         iPb = iP2s(n)

         If(iPa.NE.iPb) Then
            lNP = lNP + 1
            If(lNP.GT.MaxS) Goto 9000

            iNPs(1, lNP) = iPa
            iNPs(2, lNP) = iPb
         End if
      End do


      nArms = 0
      If(ifXnode(status, ANIUse2ArmRule)) Then
         If(ifXnode(ICP(iPd), jInode)) Then
            nArms = 2
         Else If(ifXnode(ICP(iPd), jBnode)) Then
            nArms = 1
         End if
      End if


      If(nArms.GT.0) Then
         Call chkSPf(nArms, iPd, iDELET, ICP, IEP, IPE, IEE, lP, iOs)
         Call chkSPb(2, lDP, iDPs, lNP, iNPs, iDELET,
     &               ICP, IEP, IPE, IEE, lP, iOs, iPbad, flagBNDs)
         If(flagBNDs) Goto 1000
      End if


c ... checking for boundary elements
      If(ifXnode(status, ANIForbidBoundaryElements)) Then
         Do 800 n = lEold + 1, lE
            Do i = 1, 4
               iPu = IPEs(i, n)
               If(ifXnode(ICP(iPu), jInode)) Goto 800
            End do
            Goto 1000
 800     Continue
      End if


C ... checking for inverted elements
      If(flagTM) Then
         Do n = 1, lE
            If(iEs(n).GE.0) Then
               Call updQb(n, lE, iEs, XYP, IPEs, qEs)
            End if
         End do

         mBad = 0
         Do n = lEold, lE
            If(iEs(n).GE.0 .AND. qEs(n).LE.0D0) mBad = mBad + 1
         End do

         If(mBad.GE.nBad) Goto 1000
      End if


C ... analyzing curvilinear and plane faces
      If(flagFACE) Then
         Do n = lEold + 1, lE
            Do 200 i1 = 1, 4
               i2 = ip(i1 + 1)
               i3 = ip(i2 + 1)

               iPt(1) = IPEs(i1, n)
               iPt(2) = IPEs(i2, n)
               iPt(3) = IPEs(i3, n)


               Do i = 1, 3
                  Call findSE(lP2, iP2s, iPt(i), nPt(i))
                  If(nPt(i).LE.0) Goto 200
               End do

               ICF1s = min(ICP2s(nPt(1)), ICP2s(nPt(2)))
               ICF1s = min(ICF1s,         ICP2s(nPt(3)))

               ICF2s = max(ICP2s(nPt(1)), ICP2s(nPt(2)))
               ICF2s = max(ICF2s,         ICP2s(nPt(3)))


               flagFACE2 = .FALSE.
               If(ICF1s.GT.0) Then
                  If(ICF1s.EQ.ICF2s) flagFACE2 = .TRUE.
                  If(lP1.EQ.2) Then
                     icnt = 0
                     Do j1 = 1, 3
                        If(iPt(j1).EQ.iP1s(1) .OR.
     &                     iPt(j1).EQ.iP1s(2)) Then
                           icnt = icnt + 1

                           j2 = iref(j1 + 1)
                           j3 = iref(j2 + 1)
                           l1 = ICP2s(nPt(j2)) .EQ. ICP2s(nPt(j3))
                        End if
                     End do

                     If(icnt.EQ.2) flagFACE2 = .TRUE.
                     If(icnt.EQ.1 .AND. l1) flagFACE2 = .TRUE.
                  End if
               End if

               If(flagFACE2) Then
c                 Do m = lFold + 1, lF
                  Do m = 1, lF
                     If(check33(iPt(1), iPt(2), iPt(3), IPFs(1, m),
     &                  IPFs(2, m), IPFs(3, m))) Goto 200
                  End do

                  lF = lF + 1
                  If(lE.GT.MaxS) Goto 9000

                  Call facAdd(iFs(lF), nF, MaxF, IHolF)

                  IPFs(1, lF) = iPt(1)
                  IPFs(2, lF) = iPt(2)
                  IPFs(3, lF) = iPt(3)
                  IPFs(4, lF) = ICF1s
               End if
 200        Continue
         End do
      End if


C ... updating the grid
      flag = .TRUE.

      Call pntDel(iPd, nP, ICP, IHolP)

      If(flagFACE) Then
         Do n = 1, lFold
            iFt = iFs(n)
            If(iFt.LE.0) Then
               iFt = -iFt
               Call facDel(iFt, nF, IPF, IHolF)
            End if
         End do


         Do n = lFold + 1, lF
            Call facUpd(n, IPF, iFs, IPFs)
         End do
      End if


      Do n = 1, lE
         iEt = iEs(n)
         If(iEt.EQ.0) Then
            Call eleAdd(nE, MaxE, IHolE)
            Call lstAdd(nE, L1E, nL2, L2E, nStep, IHolE,
     &                  qE, qEs(n), iEs(n))
            Call eleDel(iEs(n), IPE, IEE)
         Else If(iEt.LT.0) Then
            iEt = -iEt
            Call lstDel(nE, L1E, nL2, L2E, nStep, IHolE, qE, iEt)
            Call eleDel(iEt, IPE, IEE)
         End if
      End do


      Do n = 1, lE
         If(iEs(n).GT.0) Then
            Call eleUpd(n,  IEP, IPE, IFE,  IEE,
     &                  lF, lE,  iFs, iEs, IPFs, IPEs)
         End if
      End do


 5000 Format('Warning in deletP: bad volumes =', 2E16.8)

 1000 Return
 9000 Call errMes(1007, 'deletP', 'local parameter MaxS is small')
      End Subroutine deletP
C
      End Module mba3d_deletP
