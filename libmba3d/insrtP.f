      Module mba3d_insrtP
C
      use mba3d_auxSE
      use mba3d_auxSP
      use mba3d_auxSR
      use mba3d_error
      use mba3d_lintrp3D
      use mba3d_list
      use mba3d_makM
      use mba3d_makQ
      use mba3d_update
C
      contains
C
C ================================================================
      Subroutine insrtP(
C ================================================================
C group (M)
     &           iwR, iwE,
     &           nP, MaxP, nF, MaxF, nE, MaxE,
     &           XYP, IPF, IPE,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           L1E, L2E, nL2, nStep,
     &           IHolP, IHolF, IHolE,
     &           status,
C group (Q)
     &           HesP, rQuality, detG, qE,
     &           MetricFunction, flagAnalytic,
C group (S)
     &           lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &           nPw, nEw, XYPw, HesPw, IPEw,
     &           miLINTRP, mrLINTRP, iSE, rSE, iControl,
     &           flag)
C ================================================================
      include 'makS.fd'
      include 'color.fd'
      include 'operat.fd'
      include 'status.fd'
C ================================================================
C   Inserting a point in the middle of mesh edge. The new is
C   projected to a curvilinear surface or edge if necessary.
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
      Logical flagAnalytic
 
      Integer  MetricFunction
      EXTERNAL MetricFunction

C group (S)
      Integer iFu(*), iEu(*), IPFu(4, *), IPEu(5, *)
      Real*8  qEu(*)

C group (W)
      Integer IPEw(4, *), iSE(*)
      Real*8  XYPw(3, *), HesPw(6, *)
      Real*8  rSE(*)

C group (Flag)
      Logical flag, flagTM

C ================================================================
C group (Functions)

C group (Local variables)
      Integer iFs(MaxS), iEs(MaxS), ICFs(MaxS)
      Integer IPFs(4, MaxS), IPEs(5, MaxS)
      Real*8  qEs(MaxS), HesPs(6), detGs, XYPs(3)
      Real*8  XYPo(3)

      Integer ip(5)
      Integer iRs(3, MaxS), iPR(2, 6)
      Integer iPs(MaxS), iDPs(2, 1), iNPs(2, 1)
      Real*8  v1, v2

      Logical flagBNDs, flagLoop, flagEDGE
      Logical flagFBE

C ================================================================
      DATA iPR /1,2, 1,3, 1,4, 2,3, 2,4, 3,4/
C ================================================================

      flag = .FALSE.

      Call copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1


      iPa = IPE(iPR(1, iwR), iwE)
      iPb = IPE(iPR(2, iwR), iwE)

C ... checking the case when insertion is impossible
      Call makSR(iPa, iPb, lE, iEs, IPEs, lR, iRs, flagEDGE)
      If(.NOT.flagEDGE) Goto 1000

      flagLoop = .TRUE.
      If(iRs(3, lR).NE.iRs(2, 1)) flagLoop = .FALSE.

      lPEnd = lR
      If(.NOT.flagLoop) lPend = lR + 1


C ... checking the number of inverted elements
      flagTM = ifXnode(status, ANITangledMesh)
      if(flagTM) Then
         nBad = 0
         Do n = 1, lR
            If(qEs(iRs(1, n)).LE.0D0) nBad = nBad + 1
         End do
         If(nBad.GT.0) Goto 1000
      End if


C ... making a virtual evaluation of the quality
      Call clrSR(iPa, iPb, ICP, IPF, IFE, lF, iFs, lE, iEs, ICRab)
      ICPs = ICRab
      If(ifXnode(ICPs, jTnode)) Goto 1000
      If(ifXnode(ICP(iPa), jTnode) .AND. 
     &   ifXnode(ICP(iPb), jTnode)) Goto 1000


      Do i = 1, 3
         XYPs(i) = (XYP(i, iPa) + XYP(i, iPb)) / 2
      End do


      If(ifXnode(ICRab, jInode)) Then
         Do i = 1, 6
            HesPs(i) = (HesP(i, iPa) + HesP(i, iPb)) / 2
         End do

         LDH = 6
         nXY = 1
         If( .NOT.flagAnalytic ) Then
            Call LINTRP3D(nEw, IPEw, nPw, XYPw, LDH, HesPw, nXY, XYPs,
     &                    HesPs, iSE, miLINTRP, rSE, mrLINTRP, iControl)
         Else
            Call scaleBack(XYPs, XYPo)
            Call iniQ_analytic(1, XYPo, MetricFunction, HesPs)
         End if
      Else
         lP = lPend + 2
         Do n = 1, lPend
            If(n.LT.lPEnd) Then
               iPs(n) = iRs(2, n)
            Else
               iPs(n) = iRs(3, n - 1)
            End if
         End do

         iPs(lP - 1) = iPa
         iPs(lP) = iPb
         Call HesBnd(lP, iPs, ICP, HesP, HesPs)
      End if

      Call calDet(HesPs, detGs)


      flagFBE = ifXnode(status, ANIForbidBoundaryElements)
      flagFBE = flagFBE .AND. ifXnode(ICPs, jBnode)
      flagFBE = flagFBE .AND. 
     &        (ifXnode(ICP(iPa), jBnode) .OR. 
     &         ifXnode(ICP(iPb), jBnode)) 


      lEold = lE
      Do n = 1, lR
         iP1 = iRs(2, n)
         iP2 = iRs(3, n)

c  ...  checking for boundary elements
         If(flagFBE) Then
            If(ifXnode(ICP(iP1), jBnode) .AND. 
     &         ifXnode(ICP(iP2), jBnode)) Goto 1000 
         End if

         nE1 = iRs(1, n)
         Call calQE(
     &        HesP(1, iPa), detG(iPa), XYP(1, iPa),
     &        HesP(1, iP1), detG(iP1), XYP(1, iP1),
     &        HesP(1, iP2), detG(iP2), XYP(1, iP2),
     &        HesPs, detGs, XYPs,
     &        hStar, qEs(nE1), v1)

         If(qEs(nE1).LE.rQuality) Goto 1000

         lE = lE + 1
         If(lE.GT.MaxS) Goto 9000
         iEs(lE) = 0

         Call calQE(
     &        HesP(1, iPb), detG(iPb), XYP(1, iPb),
     &        HesP(1, iP1), detG(iP1), XYP(1, iP1),
     &        HesP(1, iP2), detG(iP2), XYP(1, iP2),
     &        HesPs, detGs, XYPs,
     &        hStar, qEs(lE), v2)

         If(qEs(lE).LE.rQuality) Goto 1000
      End do


C ... checking for surrounding points
      flagBNDs = .FALSE.
      If(ifXnode(status, ANIUse2ArmRule)) Then
         Do n = 1, lR
            iPt = iRs(2, n)
            If(ifXnode(ICP(iPt), jBnode)) flagBNDs = .TRUE.
         End do
      End if


      If(flagBNDs) Then
         iDPs(1, 1) = iPa
         iDPs(2, 1) = iPb

         If(ifXnode(ICP(iPa), jBnode) .AND. 
     &      ifXnode(ICP(iPb), jBnode)) Then
            Call chkSPf(1, iPa, iINSRT, ICP, IEP, IPE, IEE, lP, iPs)
            Call chkSPb(2, 1, iDPs, 0, iNPs, iINSRT,
     &                  ICP, IEP, IPE, IEE, lP, iPs, iPbad, flagBNDs)
            If(flagBNDs) Goto 1000

            Call chkSPf(1, iPb, iINSRT, ICP, IEP, IPE, IEE, lP, iPs)
            Call chkSPb(2, 1, iDPs, 0, iNPs, iINSRT,
     &                  ICP, IEP, IPE, IEE, lP, iPs, iPbad, flagBNDs)
            If(flagBNDs) Goto 1000
         End if
      End if


c  ...   checking for the new boundary point
      If(ifXnode(ICPs, jBnode)) Then
         If(ifXnode(status, ANIUse2ArmRule)) Then
            flagBNDs = .TRUE.
            Do n = 1, lR
               iPs(n) = iRs(2, n)
               If(ifXnode(ICP(iPt), jInode)) flagBNDs = .FALSE.
            End do
         Else
            flagBNDs = .FALSE.
         End if

         If(flagBNDs) Then
            iPs(lR + 1) = iPa
            iPs(lR + 2) = iPb

            lP = lR + 2
            Call chkSPb(1, 1, iDPs, 0, iNPs, iINSRT,
     &                  ICP, IEP, IPE, IEE, lP, iPs, iPbad, flagBNDs)
            If(flagBNDs) Goto 1000
         End if


c  ...  checking for boundary elements
         If(flagFBE .AND.
     &      ifXnode(ICP(iPa), jBnode) .AND. 
     &      ifXnode(ICP(iPb), jBnode)) Then
            Do 100 n = 1, lR
               iP1 = iRs(2, n)
               iP2 = iRs(3, n)

               If(ifXnode(ICP(iP1), jInode)) Goto 100
               If(ifXnode(ICP(iP2), jInode)) Goto 100

               Goto 1000
 100        Continue
         End if
      End if


C ... analyzing curvilinear and plane faces
      Call pntAdd(iPc, nP, MaxP, ICP,  XYP,  HesP,  detG, IHolP,
     &                           ICPs, XYPs, HesPs, detGs)

      lFold = lF

      Do n = 1, lPEnd
         If(n.LT.lPEnd) Then
            iP1 = iRs(2, n)
         Else
            iP1 = iRs(3, n - 1)
         End if

         Do m = 1, lFold
            If(check33(iPa, iPb, iP1,
     &                 IPFs(1, m), IPFs(2, m), IPFs(3, m))) Then
               nF1 = m
               iCLRf = IPF(4, iFs(nF1))

               IPFs(1, nF1) = iPa
               IPFs(2, nF1) = iP1
               IPFs(3, nF1) = iPc
               IPFs(4, nF1) = iCLRf


               lF = lF + 1
               IPFs(1, lF) = iPb
               IPFs(2, lF) = iP1
               IPFs(3, lF) = iPc
               IPFs(4, lF) = iCLRf

               ICFs(lF) = nF1

               Call facAdd(iFs(lF), nF, MaxF, IHolF)
            End if
         End do
      End do


C ... updating the grid
      flag = .TRUE.

      Do n = lFold + 1, lF
         nF1 = ICFs(n)

         Call facUpd(nF1, IPF, iFs, IPFs)
         Call facUpd(n,   IPF, iFs, IPFs)
      End do


      Do n = 1, lR
         iP1 = iRs(2, n)
         iP2 = iRs(3, n)

         nE1 = iRs(1, n)
         ide = IPEs(5, nE1)

         IPEs(1, nE1) = iPa
         IPEs(2, nE1) = iP1
         IPEs(3, nE1) = iP2
         IPEs(4, nE1) = iPc
         IPEs(5, nE1) = ide

         nE2 = lEold + n
         IPEs(1, nE2) = iPb
         IPEs(2, nE2) = iP1
         IPEs(3, nE2) = iP2
         IPEs(4, nE2) = iPc
         IPEs(5, nE2) = ide

         Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(nE1), qEs(nE1))
         Call eleDel(iEs(nE1), IPE, IEE)

         Call eleAdd(nE, MaxE, IHolE)
         Call lstAdd(nE, L1E, nL2, L2E, nStep, IHolE,
     &               qE, qEs(nE2), iEs(nE2))
         Call eleDel(iEs(nE2), IPE, IEE)
      End do


      Do n = 1, lR
         nE1 = iRs(1, n)
         Call eleUpd(nE1, IEP, IPE, IFE,  IEE,
     &               lF,  lE,  iFs, iEs, IPFs, IPEs)

         nE2 = lEold + n
         Call eleUpd(nE2, IEP, IPE, IFE,  IEE,
     &               lF,  lE,  iFs, iEs, IPFs, IPEs)
      End do

 1000 Return
 9000 Call errMes(1007, 'insrtP', 'local parameter MaxS is small')
      End Subroutine insrtP
C
      End Module mba3d_insrtP
