C =======================================================
      Subroutine clpsF1(
C =======================================================
c group (M)
     &            iwF, iwE,
     &            nP, nF, nE,
     &            XYP, IPF, IPE, 
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
     &            HesP, rQuality, detG, qE,
     &            MetricFunction, flagAnalytic,
C group (S)
     &            lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &            nPw, nEw, XYPw, HesPw, IPEw,
     &            iSE, rSE, w1, w2,
     &            flag)
C =======================================================
      include 'makS.fd'
      include 'colors.fd'
      include 'status.fd'
      include 'operat.fd'
C =======================================================
C Routine realizes one of the mesh operations: collapses
C edge iwF of element iwE to its middle point.
C
C *** DATA FLOW CHART
C
C  -> simple check if the edge can be collapsed
C  -> collect information about the edge
c
c  -> define a point to collapse the edge
c  ----> mid-point of the edge 
c  ----> projection onto the curved boundary
c  ----> one of the terminal points of the edge 
c
c  -> virtual evaluation of the element quality
c
c  -> check that no boundary triangles were created
c  -> check that 2-arm rule helds
c  -> check that no traingles were tangled
c
c  -> update the quality of new mesh elements
c  -> update the list of curvilinear edges
c  -> update mesh cross-references for the new elements  
C =======================================================
C group (M)
      Integer IPF(4, *), IPE(3, *)
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
      real  HesP(3, *), rQuality
      real  detG(*), qE(*)
      Logical flagAnalytic
 
      Integer  MetricFunction
      EXTERNAL MetricFunction

C group (S)
      Integer iFu(*), iEu(*), IPFu(3, *), IPEu(3, *)
      real  qEu(*)

C group (W)
      Integer IPEw(3, *), iSE(*)
      real  XYPw(2, *), HesPw(3, *)
      real  rSE(*), w1, w2

C group (Flag)
      Logical flag

C group (Local variables)
      Integer iFs(MaxS), iEs(MaxS), IPFs(2, MaxS), IPEs(3, MaxS)
      real  prjXYPs(2), XYPs(2), HesPs(3), detGs, qEs(MaxS)
      real  XYPo(2)

      Integer ip(4), iPs(MaxS), iOs(MaxS)
      real  par(6), t1, t2, tc
      Logical flagBNDs, flagOrient, flagFirst, flagTM
      
      Logical chkTangled, ifXnode
      Integer minClr

C =======================================================
      flag = .FALSE.

      Call copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      i1 = iwF
      i2 = ip(i1 + 1)
      i3 = ip(i2 + 1)

      iP1 = IPE(i1, iwE)
      iP2 = IPE(i2, iwE)
      iP3 = IPE(i3, iwE)

      iF = IFE(iwF, iwE)

      iE1 = iwE
      iE2 = IEE(iwF, iE1)

C ... checking for faces which can not be collapsed
      ICP1 = ICP(iP1)
      ICP2 = ICP(iP2)
      ICP3 = ICP(iP3)


      ICPs = minClr(ICP1, ICP2)
      If(ifXnode(ICPs, jVnode)) goto 1000
      If(ifXnode(ICPs, jTnode)) goto 1000

      If(ifXnode(ICPs, jSnode) .AND. iF.EQ.0) goto 1000

      ICPt = maxClr(ICP1, ICP2)
      If(ifXnode(ICPt, jVnode) .AND. ifXnode(ICPt, jTnode)) Goto 1000


C ... gathering information about the face
      Do i = 1, 6
         par(i) = 0D0
      End do

      If(iF.NE.0) Then
         Call infoF(iF, iP1, iP2, iF1, iF2, iPc, iPd,
     &              par, IPF, parCrv, lF, iFs)

         If(iPc.LE.0 .OR. iPd.LE.0) goto 1000
      End if


C ... checking for inverted elements
      flagTM = ifXnode(status, ANIUntangleMesh)

      if(flagTM) Then
         nBad = 0
         Do n = 1, lE
            If(qEs(n).LE.0D0) nBad = nBad + 1
         End do
      End if

      if(flagTM) Then
         flagTM = nBad.GT.0
      End if
      

C ... finding a point in which we collapse the edge
      iFNCs = 0
      If(ICP1.EQ.jInode .AND. ICP2.EQ.jInode) Then
         ICPs = jInode
         Do i = 1, 2
            XYPs(i) = XYP(i, iP1) * w1 + XYP(i, iP2) * w2
         End do

         Do i = 1, 3
            HesPs(i) = HesP(i, iP1) * w1 + HesP(i, iP2) * w2
         End do

         If(.NOT.flagAnalytic) Then
         Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1, 
     &               XYPs, HesPs, iSE, rSE, .FALSE.)
         Else
            Call scaleBack(XYPs, XYPo)
            Call iniQ_analytic(1, XYPo, MetricFunction, HesPs)
         End if

         Call calDet(HesPs, detGs)
       Else If(ifXnode(ICP1, jSnode) .AND. ICP2.EQ.jInode .OR.
     &         ifXnode(ICP1, jTnode) .AND. ICP2.EQ.jInode .OR.
     &         ifXnode(ICP1, jTnode) .AND. ifXnode(ICPs, jBnode) .OR.
     &         ifXnode(ICP1, jVnode)) Then
         ICPs = ICP1
         Do i = 1, 2
            XYPs(i) = XYP(i, iP1)
         End do

         Do i = 1, 3
            HesPs(i) = HesP(i, iP1)
         End do
         detGs = detG(iP1)

         t1 = par(2)
         t2 = par(3)
       Else If(ifXnode(ICP2, jSnode) .AND. ICP1.EQ.jInode .OR.
     &         ifXnode(ICP2, jTnode) .AND. ICP1.EQ.jInode .OR.
     &         ifXnode(ICP2, jTnode) .AND. ifXnode(ICPs, jBnode) .OR.
     &         ifXnode(ICP2, jVnode)) Then
         ICPs = ICP2
         Do i = 1, 2
            XYPs(i) = XYP(i, iP2)
         End do

         Do i = 1, 3
            HesPs(i) = HesP(i, iP2)
         End do
         detGs = detG(iP2)

         t1 = par(4)
         t2 = par(5)

C ...    changing order of points to be consistent with the previous case
         i = iP2
         iP2 = iP1
         iP1 = i
      Else If(ifXnode(ICPs, jSnode)) Then
c!       ICPs = minClr(ICP1, ICP2)
         iCRVs = IPF(3, iF)
         If(iCRVs.NE.0) Then
            t1 = par(3)
            t2 = par(4)
            tc = (t1 + t2) / 2

            iFNCs = iFnc(iCRVs)
            Call aniCrv(tc, XYPs, iFNCs, calCrv)

            t1 = tc
            t2 = tc

            Call findSE(nCrvFnc, LFnc, iFNCs, k)
            ir = ILt(k)
            Call prjCrv(XYPs, prjXYPs, iFNCs, tc, calCrv,
     &                  L1Et(1, ir), L2Et(ir), nL2t(k), nStept(1, k),
     &                  nEt(k), tE(ir))
  
            If(.NOT.flagAnalytic) Then
            Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1, 
     &                  prjXYPs, HesPs, iSE, rSE, .FALSE.)
         Else
               Call scaleBack(prjXYPs, XYPo)
               Call iniQ_analytic(1, XYPo, MetricFunction, HesPs)
            End if
         Else
            Do i = 1, 2
               XYPs(i) = (XYP(i, iP1) + XYP(i, iP2)) / 2
            End do

            If(.NOT.flagAnalytic) Then
            Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1, 
     &                  XYPs, HesPs, iSE, rSE, .FALSE.)
            Else
               Call scaleBack(XYPs, XYPo)
               Call iniQ_analytic(1, XYPo, MetricFunction, HesPs)
            End if
         End if

         Call calDet(HesPs, detGs)
      Else
         goto 1000
      End if


C ... virtual evaluation of the superelement quality
      Do 10 n = 1, lE
         iE = iEs(n)
         If(iE.EQ.iE1 .OR. iE.EQ.iE2) Then
            iEs(n) = -iEs(n)
            goto 10
         End if

         Do i1 = 1, 3
            If(IPEs(i1, n).EQ.iP1 .OR. IPEs(i1, n).EQ.iP2) Then
               i2 = ip(i1 + 1)
               i3 = ip(i2 + 1)

               IPEs(i1, n) = iP1
               iPa = IPEs(i2, n)
               iPb = IPEs(i3, n)

               Call calQE(
     &              HesP(1, iPa), XYP(1, iPa),
     &              HesP(1, iPb), XYP(1, iPb),
     &              HesPs,        XYPs,
     &              hStar, qEs(n))

               If((1.05*qEs(n)).LE.rQuality) goto 1000
               goto 10
            End if
         End do

         iEs(n) = 0
 10   Continue


C ... checking for boundary triangles
      If(ifXnode(status, ANIForbidBoundaryElements)) Then
         Do 20 n = 1, lE
            iE = iEs(n)
            If(iE.LE.0) goto 20

            Do i = 1, 3
               iPt = IPEs(i, n)
               If(ifXnode(ICP(iPt), jInode)) goto 20
            End do

            goto 1000
 20      Continue
      End if


C ... checking for surrounding points (not ICP1 but ICP(iP1))
      If(ifXnode(status, ANIUse2ArmRule)) Then
         If(ifXnode(ICP(iP1), jBnode) .AND. ICP(iP2).EQ.jInode) Then
            Call chkSPf(iP1, iP2, ICP, IEP, IPE, IEE, lP, iPs)
            Call chkSPb(iP1, iP2, 0, 0, iCLPS,
     &                  ICP, IEP, IPE, IEE, lP, iPs, flagBNDs)
            If(flagBNDs) goto 1000
         End if
      End if


C ... checking for orientation of triangles
      Call calSO(XYP, IPE, lE, iEs, iOs)
      Call chkSO(iP1, iP2, XYPs, XYP, IPE, lE, iEs, iOs, flagOrient)
      If(.NOT.flagOrient) goto 1000


C ... checking for inverted elements
      If(flagTM) Then
         Do n = 1, lE
            If(iEs(n).GE.0 .AND. qEs(n).GT.0D0) Then
               Call updQb(n, lE, iEs, XYP, IPEs, qEs)
            End if
         End do

         mBad = 0
         Do n = 1, lE
            If(iEs(n).GE.0 .AND. qEs(n).LE.0D0) mBad = mBad + 1
         End do

         If(mBad.GE.nBad) goto 1000

C  ...  colapsing may result in topologically wrong mesh
        flagTM = chkTangled(lE, iEs, IPEs)
        If(flagTM) goto 1000
      End if


c ... updating the quality
      If(ifXnode(status, ANISmoothMesh)) Then
         flagFirst = .TRUE.
         Call updQE(XYP, lE, iEs, IPEs,
     &              HesP, rQuality, detG, hStar, qEs, flagFirst)
         If(.NOT.flagFirst) goto 1000
      End if 


C ... updating the grid
      flag = .TRUE.

      If(ifXnode(ICP1, jTnode) .or. ifXnode(ICP2, jTnode)) then
        call addXnode(ICPs, jTnode)
      end if

      Call pntUpd(iP1, ICP,  XYP,  HesP,  detG,
     &                 ICPs, XYPs, HesPs, detGs)
      Call pntDel(iP2, nP, ICP, IHolP)


      If(ifXnode(ICP1, jSnode) .AND. ifXnode(ICP2, jSnode)) Then
         Call findSE(lF, iFs, iF, nFd)
         Call facDel(iF, nF, IPF, iFnc, IHolF)
         iFs(nFd) = -iFs(nFd)

         Call findSE(lF, iFs, iF1, nF1)
         Call findSE(lF, iFs, iF2, nF2)

         IPFs(1, nF1) = iPc
         IPFs(2, nF1) = iP1

         IPFs(1, nF2) = iP1
         IPFs(2, nF2) = iPd

         iCRVs = IPF(3, iF1)
         iFNCs = iFnc(iF1)
         iBNDs = IPF(4, iF1)

         Call facUpd(nF1, IPF, parCrv, iFnc,
     &               iFs, IPFs, iCRVs, iFNCs, iBNDs, par(1), t1)

         iCRVs = IPF(3, iF2)
         iFNCs = iFnc(iF2)
         iBNDs = IPF(4, iF2)
         Call facUpd(nF2, IPF, parCrv, iFnc,
     &               iFs, IPFs, iCRVs, iFNCs, iBNDs, t2, par(6))
      End if


      Do n = 1, lE
         iEt = iEs(n)
         If(iEt.LT.0) Then
            iEt = -iEt
            Call lstDel(nE, L1E, nL2, L2E, nStep, IHolE, qE, iEt)
            Call eleDel(iEt, IPE, IEE)
         Else If(iEt.GT.0) Then
            Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEt, qEs(n))

            Call eleUpd(n, IEP, IPE, IFE, IEE,
     &                  lF, lE, iFs, iEs, IPFs, IPEs)
         End if
      End do

1000  Return
      End



C =======================================================
      Subroutine clpsF2(
C =======================================================
c group (M)
     &            iwP, iwE,
     &            nP, nF, nE,
     &            XYP, IPF, IPE, 
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
     &            HesP, rQuality, detG, qE,
     &            MetricFunction, flagAnalytic,
C group (S)
     &            lEu, iEu, 
C group (W)
     &            nPw, nEw, XYPw, HesPw, IPEw,
     &            iSE, rSE,
     &            flag)
C =======================================================
      include 'makS.fd'
      include 'colors.fd'
      include 'status.fd'
C =======================================================
C Routine realizes one of the mesh operations: collapses
C an edge ending at point iwP of element iwE. It uses the
C main routine clpsF1.
C =======================================================
C group (M)
      Integer IPF(4, *), IPE(3, *)
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
      real  HesP(3, *), rQuality
      real  detG(*), qE(*)
      Logical flagAnalytic
 
      Integer  MetricFunction
      EXTERNAL MetricFunction

C group (S)
      Integer iEu(*)

C group (W)
      Integer IPEw(3, *), iSE(*)
      real  XYPw(2, *), HesPw(3, *)
      real  rSE(*)

C group (Flag)
      Logical flag

C group (Local variables)
      Integer iFs(MaxS), iEs(MaxS), IPFs(2, MaxS), IPEs(3, MaxS)
      real  qEs(MaxS)

      Integer ip(4)
      real  d, w1, w2

      Logical ifXnode

C =======================================================
      flag = .FALSE.

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

C ... checking points which can not be involved
      iP1 = IPE(iwP, iwE)
      If(ifXnode(ICP(iP1), jSnode)) goto 1000
      If(ifXnode(ICP(iP1), jTnode)) goto 1000

      i2 = ip(iwP + 1)
      i3 = ip(i2 + 1)
      iPa = IPE(i2, iwE)
      iPb = IPE(i3, iwE)

      Do 10 k = 1, lEu
         kE = iEu(k)
         If(kE.EQ.iwE) goto 10

         Do i = 1, 3
            If(IPE(i, kE).EQ.iP1) Then
               i1 = i
               goto 5
            End if
         End do
         goto 10

 5       i2 = ip(i1 + 1)
         i3 = ip(i2 + 1)

         iP2 = IPE(i2, kE)
         iP3 = IPE(i3, kE)

         d = (XYP(1, iP2) - XYP(1, iP1)) *
     &       (XYP(2, iP3) - XYP(2, iP1)) -
     &       (XYP(1, iP3) - XYP(1, iP1)) *
     &       (XYP(2, iP2) - XYP(2, iP1))

         If(d.LT.0D0) Then
            iwF = min(i1, i3)
            itF = max(i1, i3)
            If(iP3.EQ.iPa .OR. iP3.EQ.iPb) goto 10
         Else
            iwF = min(i1, i2)
            itF = max(i1, i2)
            If(iP2.EQ.iPa .OR. iP2.EQ.iPb) goto 10
         End if
         If(iwF.EQ.1 .AND. itF.EQ.3) iwF = 3

         If(iwF.EQ.i1) Then
            w1 = 0D0
            w2 = 1D0
         Else
            w1 = 1D0
            w2 = 0D0
         End if

         Call makSE(kE,  IEP, IPF, IPE, IFE,  IEE,  qE, MaxS,
     &              lFs, lEs, iFs, iEs, IPFs, IPEs, qEs,
     &              status)

         Call clpsF1(
c group (M)
     &            iwF, kE,
     &            nP, nF, nE,
     &            XYP, IPF, IPE, 
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
     &            HesP, rQuality, detG, qE,
     &            MetricFunction, flagAnalytic,
C group (S)
     &            lFs, lEs, iFs, iEs, IPFs, IPEs, qEs,
C group (W)
     &            nPw, nEw, XYPw, HesPw, IPEw,
     &            iSE, rSE, w1, w2,
     &            flag)
        If(flag) goto 1000
 10   Continue

 1000 Return
      End


