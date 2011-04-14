C ====================================================
      Subroutine insrtP(
C ====================================================
c group (M)
     &            iwF, iwE,
     &            nP, MaxP, nF, MaxF, nE, MaxE,
     &            XYP, IPF, IPE, lbE,
     &            calCrv, parCrv, iFnc,
     &            hStar,
     &            ICP, IEP, IFE, IEE,
     &            L1E, L2E, nL2, nStep,
     &            IHolP, IHolF, IHolE,
     &            status,
c group (CRV)
     &            L1Et, L2Et, tE,
     &            nL2t, nStept, nEt,
     &            nCrvFnc, LFnc, ILt,
C group (Q)
     &            HesP, Quality, rQuality, detG, qE,
     &            MetricFunction, flagAnalytic,
C group (S)
     &            lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &            nPw, nEw, XYPw, HesPw, IPEw,
     &            iSE, rSE,
     &            flag)
C ====================================================
      include 'makS.fd'
      include 'colors.fd'
      include 'status.fd'
      include 'operat.fd'
C ====================================================
C Routine realizes one of the mesh operations: inserts
C a point at the middle of edge iwF of element iwE.
C
C *** DATA FLOW CHART
C
C  -> check if the edge can be split
C  -> check that 2-arm rule will be preserved
C  -> insert a virtual point 
C  -> project the point onto curved boundary
C
C  -> virtual evaluation of the element quality
C
C  -> check that no traingles were tangled (for curved edge)
C  -> check that "smoothed" element quality has incresed
C
C  -> update the quality of new mesh elements
C  -> update the list of curvilinear edges
C  -> update mesh cross-references for the new elements
C
C ====================================================
C group (M)
      Integer IPF(4, *), IPE(3, *), lbE(*)
      real  XYP(2, *)

      EXTERNAL calCrv
      Integer iFnc(*)
      real  parCrv(2, *), hStar

      Integer L1E(2, *), L2E(*), nStep(4)
      Integer IHolP(*), IHolF(*), IHolE(*)

      Integer ICP(*), IEP(*)
      Integer IFE(3, *), IEE(3, *)

      Integer status

C group (CRV)
      real  tE(*)
      Integer L1Et(2, *), L2Et(*)
      Integer nL2t(*), nStept(4, *), nEt(*)
      Integer LFnc(*), ILt(*)

C group (Q)
      real  HesP(3, *), Quality, rQuality
      real  detG(*), qE(*)
      Logical flagAnalytic
 
      Integer  MetricFunction
      EXTERNAL MetricFunction

C group (S)
      Integer iFu(*), iEu(*), IPFu(2, *), IPEu(3, *)
      real  qEu(*)

C group (W)
      Integer IPEw(3, *), iSE(*)
      real  XYPw(2, *), HesPw(3, *)
      real  rSE(*)

C group (Flag)
      Logical flag

C group (Local variables)
      Integer iFs(MaxS), iEs(MaxS), IPFs(2, MaxS), IPEs(3, MaxS)
      real  prjXYPs(2), XYPs(2), HesPs(3), detGs, qEs(MaxS)
      real  XYPt(2), XYPo(2)

      Integer iPs(MaxS), ip(4)
      real  calVol, vol1, vol2
      real  par(6), t1, t2, tc
      Logical flagBNDs, flagFirst

      Integer minClr
      Logical ifXnode

C ====================================================
      flag = .FALSE.

      Call copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      i1 = iwF
      i2 = ip(i1 + 1)

      iP1 = IPE(i1, iwE)
      iP2 = IPE(i2, iwE)


c ... checking for faces which can not be split
      ICPs = minClr(ICP(iP1), ICP(iP2))
      If(ifXnode(ICPs, jTnode)) goto 1000


      iE1 = iwE
      iE2 = IEE(iwF, iE1)

      i3 = ip(i2 + 1)
      iPa = IPE(i3, iwE)

      If(iE2.NE.0) Then
         Do i = 1, 3
            If(IPE(i, iE2).NE.iP1 .AND. IPE(i, iE2).NE.iP2) Then
               i4 = i
               iPb = IPE(i, iE2)
            End if
         End do
      End if


C ... checking for surrounding points
      iF = IFE(iwF, iwE)
      If(iF.EQ.0) Then
         ICPs = jInode
      Else
         lP = 2
         iPs(1) = iP1
         iPs(2) = iP2

         If(ifXnode(status, ANIUse2ArmRule)) Then
            Call chkSPb(iP1, iP2, 0, 0, iINSRT,
     &                  ICP, IEP, IPE, IEE, lP, iPs, flagBNDs)
            If(flagBNDs) goto 1000
         End if

         If(iE2.NE.0) Then
            ICPs = jInode + jSnode
         Else
            ICPs = jBnode + jSnode
         End if
      End if


C ... skipping the inverted elements
      If(ifXnode(status, ANIUntangleMesh)) Then
         If(qE(iE1).LE.0D0) goto 1000
         If(iE2.GT.0) Then
            If(qE(iE2).LE.0D0) goto 1000
         End if
      End if


C ... making the point for inserting
      iFNCs = 0
      Do i = 1, 3
         HesPs(i) = (HesP(i, iP1) + HesP(i, iP2)) / 2
      End do

      If(ifXnode(ICPs, jSnode)) Then
         Call infoF(iF, iP1, iP2, iF1, iF2, iPc, iPd,
     &              par, IPF, parCrv, lF, iFs)

         iCRVs = IPF(3, iF)
         If(iCRVs.NE.0) Then
            t1 = par(3)
            t2 = par(4)
            tc = (t1 + t2) / 2

            iFNCs = iFNC(iCRVs)
            Call aniCrv(tc, XYPs, iFNCs, calCrv)

            Call findSE(nCrvFnc, LFnc, iFNCs, k)
            ir = ILt(k)
            Call prjCrv(XYPs, prjXYPs, iFNCs, tc, calCrv,
     &                  L1Et(1, ir), L2Et(ir), nL2t(k), nStept(1, k),
     &                  nEt(k), tE(ir))

            If( .NOT.flagAnalytic ) Then
            Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1,
     &                  prjXYPs, HesPs, iSE, rSE, .FALSE.)
            Else
               Call scaleBack(prjXYPs, XYPo)
               Call iniQ_analytic(1, XYPo, MetricFunction, HesPs)
            End if

            Do i = 1, 2
               XYPt(i) = (XYP(i, iP1) + XYP(i, iP2)) / 2
            End do
         Else
            Do i = 1, 2
               XYPs(i) = (XYP(i, iP1) + XYP(i, iP2)) / 2
            End do

            If( .NOT.flagAnalytic ) Then
            Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1,
     &                  XYPs, HesPs, iSE, rSE, .FALSE.)
            Else
               Call scaleBack(XYPs, XYPo)
               Call iniQ_analytic(1, XYPo, MetricFunction, HesPs)
            End if
         End if
      Else
         iCRVs = 0
         Do i = 1, 2
            XYPs(i) = (XYP(i, iP1) + XYP(i, iP2)) / 2
         End do

         If( .NOT.flagAnalytic ) Then
         Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1,
     &               XYPs, HesPs, iSE, rSE, .FALSE.)
         Else
            Call scaleBack(XYPs, XYPo)
            Call iniQ_analytic(1, XYPo, MetricFunction, HesPs)
         End if
      End if


      Call calDet(HesPs, detGs)


C ... making a virtual evaluation of the quality
      lEadd = 0
      Call findSE(lE, iEs, iE1, nE1)
      Call calQE(
     &     HesP(1, iP1), XYP(1, iP1),
     &     HesP(1, iPa), XYP(1, iPa),
     &     HesPs,        XYPs,
     &     hStar, qEs(nE1))

      If(qEs(nE1).LE.rQuality) goto 1000

      lE = lE + 1
      Call calQE(
     &     HesP(1, iP2), XYP(1, iP2),
     &     HesP(1, iPa), XYP(1, iPa),
     &     HesPs,        XYPs,
     &     hStar, qEs(lE))

      If(qEs(lE).LE.rQuality) goto 1000


C  ...  checking for the orientation
      If(iCRVs.NE.0) Then
         vol1 = calVol(XYP(1, iPa), XYP(1, iP1), XYPs)
         vol2 = calVol(XYP(1, iPa), XYP(1, iP1), XYPt)

         If(vol1 * vol2.LE.0D0) goto 1000

         vol1 = calVol(XYP(1, iPa), XYP(1, iP2), XYPs)
         vol2 = calVol(XYP(1, iPa), XYP(1, iP2), XYPt)

         If(vol1 * vol2.LE.0D0) goto 1000
      End if


      If(iE2.NE.0) Then
         lEadd = 1
         Call findSE(lE, iEs, iE2, nE2)
         Call calQE(
     &        HesP(1, iP1), XYP(1, iP1),
     &        HesP(1, iPb), XYP(1, iPb),
     &        HesPs,        XYPs,
     &        hStar, qEs(nE2))

         If(qEs(nE2).LE.rQuality) goto 1000

         lE = lE + 1
         Call calQE(
     &        HesP(1, iP2), XYP(1, iP2),
     &        HesP(1, iPb), XYP(1, iPb),
     &        HesPs,        XYPs,
     &        hStar, qEs(lE))

         If(qEs(lE).LE.rQuality) goto 1000


C  ...  checking for the orientation
         If(iCRVs.NE.0) Then
            vol1 = calVol(XYP(1, iPb), XYP(1, iP1), XYPs)
            vol2 = calVol(XYP(1, iPb), XYP(1, iP1), XYPt)

            If(vol1 * vol2.LE.0D0) goto 1000

            vol1 = calVol(XYP(1, iPb), XYP(1, iP2), XYPs)
            vol2 = calVol(XYP(1, iPb), XYP(1, iP2), XYPt)

            If(vol1 * vol2.LE.0D0) goto 1000
         End if
      End if


C ... checking for the triangles orientation
      If(iF.NE.0 .AND. iCRVs.NE.0) Then
         vol1 = calVol(XYP(1, iPa), XYPs, XYP(1, iP2))
         vol2 = calVol(XYP(1, iPa), XYPs, XYP(1, iP1))
         If(vol1 * vol2.GE.0D0) goto 1000

         If(iE2.NE.0) Then
            vol1 = calVol(XYP(1, iPb), XYPs, XYP(1, iP2))
            vol2 = calVol(XYP(1, iPb), XYPs, XYP(1, iP1))
            If(vol1 * vol2.GE.0D0) goto 1000
         End if
      End if



C ... updating the grid
      flag = .TRUE.

      Call pntAdd(iPc, nP, MaxP, ICP,  XYP,  HesP,  detG, IHolP,
     &                           ICPs, XYPs, HesPs, detGs)

      IPEs(1, nE1) = iP1
      IPEs(2, nE1) = iPc
      IPEs(3, nE1) = iPa

      IPEs(1, lE - lEadd) = iP2
      IPEs(2, lE - lEadd) = iPc
      IPEs(3, lE - lEadd) = iPa

      If(ifXnode(status, ANISmoothMesh)) Then
         flagFirst = .TRUE.
         Call updQE(XYP, lE, iEs, IPEs,
     &              HesP, rQuality, detG, hStar, qEs, flagFirst)
         If(.NOT.flagFirst) Then
            Call pntDel(iPc, nP, ICP, IHolP)
            flag = .FALSE.

            goto 1000
         End if
      End if 


      If(iF.NE.0) Then
         Call facAdd(iF2, nF, MaxF, IHolF)

         Call findSE(lF, iFs, iF, nF1)
         IPFs(1, nF1) = iP1
         IPFs(2, nF1) = iPc

         lF = lF + 1
         iFs(lF) = iF2
         IPFs(1, lF) = iPc
         IPFs(2, lF) = iP2

         iCRVs = IPF(3, iF)
         iFNCs = iFnc(iF)
         iBNDs = iPF(4, iF)
         Call facUpd(nF1, IPF, parCrv, iFnc,
     &               iFs, IPFs, iCRVs, iFNCs, iBNDs, par(3), tc)

         If(iCRVs.NE.0) iCRVs = iF2
         Call facUpd(lF,  IPF, parCrv, iFnc,
     &               iFs, IPFs, iCRVs, iFNCs, iBNDs, tc, par(4))
      End if


      Call eleDel(iEs(nE1), IPE, IEE)

      If(iE2.NE.0) Then
         IPEs(1, nE2) = iP1
         IPEs(2, nE2) = iPc
         IPEs(3, nE2) = iPb

         IPEs(1, lE) = iP2
         IPEs(2, lE) = iPc
         IPEs(3, lE) = iPb

         Call eleDel(iEs(nE2), IPE, IEE)
      End if

      Do n = lE - lEadd, lE
         Call eleAdd(nE, MaxE, IHolE)
         Call lstAdd(nE, L1E, nL2, L2E, nStep, IHolE,
     &               qE, qEs(n), iEs(n))
         Call eleDel(iEs(n), IPE, IEE)
      End do

      Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(nE1), qEs(nE1))
      If(iE2.NE.0) Then
         Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(nE2), qEs(nE2))
      End if

      Call eleUpd(nE1, IEP, IPE, IFE, IEE,
     &            lF, lE, iFs, iEs, IPFs, IPEs)
      Call eleUpd(lE, IEP, IPE, IFE, IEE,
     &            lF, lE, iFs, iEs, IPFs, IPEs)

      lbE(iEs(lE - lEadd)) = lbE(iE1)

      If(iE2.NE.0) Then
         Call eleUpd(nE2, IEP, IPE, IFE, IEE,
     &               lF, lE, iFs, iEs, IPFs, IPEs)
         Call eleUpd(lE - 1, IEP, IPE, IFE, IEE,
     &               lF, lE, iFs, iEs, IPFs, IPEs)

         lbE(iEs(lE)) = lbE(iE2)
      End if

 1000 Return
      End


