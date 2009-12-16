      Module mba3d_clpsR
C
      use mba3d_auxSE
      use mba3d_lintrp3D
      use mba3d_list
      use mba3d_makM
      use mba3d_makQ
      use mba3d_update
C
      contains
C
C ================================================================
      Subroutine clpsR(
C ================================================================
c group (M)
     &           iwR, iwE,
     &           nP, nE, XYP, IPE,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           L1E, L2E, nL2, nStep,
     &           IHolP, IHolE,
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
      include 'shutF.fd'
      include 'status.fd'
C ================================================================
C group (M)
      Integer IPE(5, *)
      Real*8  XYP(3, *)

      Real*8  hStar

      Integer L1E(2, *), L2E(*), nStep(4)
      Integer IHolP(*), IHolE(*)

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
      Logical flag, flagOrient

C ================================================================
C group (Local functions)

C group (Local variables)
      Integer iFs(MaxS), iEs(MaxS), IPFs(4, MaxS), IPEs(5, MaxS)
      Integer iPf(MaxS)
      Real*8  qEs(MaxS)

      Integer ip(5), ipr(2, 6), iDs(2), iNs(2)
      Real*8  XYPs(3), HesPs(6), detGs
      Real*8  XYPo(3), w1, w2, v
      Logical flagTM, flagBNDs

      DATA    ipr/1,2, 1,3, 1,4, 2,3, 2,4, 3,4/

C ================================================================
      flag = .FALSE.

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1

      i1 = ipr(1, iwR)
      i2 = ipr(2, iwR)

      iP1 = IPE(i1, iwE)
      iP2 = IPE(i2, iwE)


C ... check the case when the edge can not be collapsed
      ICP1 = ICP(iP1)
      ICP2 = ICP(iP2)

      ICPs = minClr(ICP1, ICP2)
      If(ifXnode(ICPs, jVnode)) Goto 1000
      If(ifXnode(ICPs, jTnode)) Goto 1000


C ... add miscaleneous restrictions including missing algorithms
      If(ifXnode(ICP1, jSnode) .OR. ifXnode(ICP2, jSnode)) Goto 1000


C ... number of inverted elements (the code is missing)
      flagTM = ifXnode(status, ANITangledMesh)
      if(flagTM) Then
         Goto 1000  

C        nBad = 0
C        If(qEs(iE1).LE.0D0) nBad = nBad + 1
C        If(qEs(iE2).LE.0D0) nBad = nBad + 1
      End if


C ... finding a point in which we collapse the edge
      w1 = 0.5D0
      w2 = 0.5D0

      iFNCs = 0
      If(ICP1.EQ.jInode .AND. ICP2.EQ.jInode) Then
         ICPs = jInode
         Do i = 1, 3
            XYPs(i) = XYP(i, iP1) * w1 + XYP(i, iP2) * w2
         End do

         Do i = 1, 6
            HesPs(i) = HesP(i, iP1) * w1 + HesP(i, iP2) * w2
         End do

         If(.NOT.flagAnalytic) Then
            LDH = 6
            nXY = 1
            Call LINTRP3D(nEw, IPEw, nPw, XYPw, LDH, HesPw, nXY, XYPs,
     &                    HesPs, iSE, miLINTRP, rSE, mrLINTRP, iControl)
         Else
            Call scaleBack(XYPs, XYPo)
            Call iniQ_analytic(1, XYPo, MetricFunction, HesPs)
         End if

 
         Call calDet(HesPs, detGs)

      Else If(ifXnode(ICP1, jSnode) .AND. ICP2.EQ.jInode .OR.
     &        ifXnode(ICP1, jTnode) .AND. ICP2.EQ.jInode .OR.
     &        ifXnode(ICP1, jVnode)) Then
         ICPs = ICP1
         Do i = 1, 3
            XYPs(i) = XYP(i, iP1)
         End do

         Do i = 1, 6
            HesPs(i) = HesP(i, iP1)
         End do
         detGs = detG(iP1)

      Else If(ifXnode(ICP2, jSnode) .AND. ICP1.EQ.jInode .OR.
     &        ifXnode(ICP2, jTnode) .AND. ICP1.EQ.jInode .OR.
     &        ifXnode(ICP2, jVnode)) Then
         ICPs = ICP2
         Do i = 1, 3
            XYPs(i) = XYP(i, iP2)
         End do

         Do i = 1, 6
            HesPs(i) = HesP(i, iP2)
         End do
         detGs = detG(iP2)

C ...    changing order of points to be consistent with the previous case
         i = iP2
         iP2 = iP1
         iP1 = i

      Else
         Goto 1000
      End if

 
C ... saving the initial superelement structure
      Call copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)


C ... making a virtual evaluation of the quality
      Do 10 n = 1, lE
         iE = iEs(n)
         
         If(check1j(iP1, IPE(1, iE)) .AND. 
     &      check1j(iP2, IPE(1, iE))) Then
            iEs(n) = -iEs(n)
            Goto 10
         End if

         Do i1 = 1, 4
            If(IPEs(i1, n).EQ.iP1 .OR. IPEs(i1, n).EQ.iP2) Then
               i2 = ip(i1 + 1)
               i3 = ip(i2 + 1)
               i4 = ip(i3 + 1)

               IPEs(i1, n) = iP1
               iPa = IPEs(i2, n)
               iPb = IPEs(i3, n)
               iPc = IPEs(i4, n)

               Call calQE(
     &              HesPs, detGs, XYPs,
     &              HesP(1, iPa), detG(iPa), XYP(1, iPa),
     &              HesP(1, iPb), detG(iPb), XYP(1, iPb),
     &              HesP(1, iPc), detG(iPc), XYP(1, iPc),
     &              hStar, qEs(n), v)

               If(qEs(n).LE.rQuality) Goto 1000

c  ...  check orientation of neigboors
               iEt = IEE(i1, iE)
               If(iEt.GT.0) Then
                  Call chkSOtet(iP1, iP2, iPa, iPb, XYPs, XYP, 
     &                          IPE(1, iEt), v, flagOrient) 
                  If(.NOT.flagOrient) Goto 1000
               End if

               iEt = IEE(i3, iE)
               If(iEt.GT.0) Then
                  Call chkSOtet(iP1, iP2, iPb, iPc, XYPs, XYP, 
     &                          IPE(1, iEt), v, flagOrient) 
                  If(.NOT.flagOrient) Goto 1000
               End if

               iEt = IEE(i4, iE)
               If(iEt.GT.0) Then
                  Call chkSOtet(iP1, iP2, iPc, iPa, XYPs, XYP, 
     &                          IPE(1, iEt), v, flagOrient) 
                  If(.NOT.flagOrient) Goto 1000
               End if

               Goto 10
            End if
         End do

         iEs(n) = 0
 10   Continue


C ... checking for boundary elements
      If(ifXnode(status, ANIForbidBoundaryElements)) Then
         Do 20 n = 1, lE
            iE = iEs(n)
            If(iE.LE.0) Goto 20

            Do i = 1, 4
               iPt = IPEs(i, n)
               If(ifXnode(ICP(iPt), jInode)) Goto 20
            End do

            Goto 1000
 20      Continue
      End if


C ... checking for surrounding points (not ICP1 but ICP(iP1))
c ... We need to study only pathes which were ending at iP2.
c ... Thus, each boundary point connected to iP2 should have 
c ... another interior point neighboor.
      If(ifXnode(status, ANIUse2ArmRule)) Then
         If(ifXnode(ICP(iP1), jBnode) .AND. ICP(iP2).EQ.jInode) Then
            iDs(1) = iP2
            iDs(2) = iP1
            Call chkSPf(2, iP2, iCLPSR, ICP, IEP, IPE, IEE, lP, iPf)
            Call chkSPb(2, 1, iDs, 0, iNs, iCLPSR,
     &                  ICP, IEP, IPE, IEE, lP, iPf, iPbad, flagBNDs)
            If(flagBNDs) Goto 1000
         End if
      End if


C ... checking for inverted elements
      If(flagTM) Then
c        Do n = 1, lE
c           If(iEs(n).GE.0 .AND. qEs(n).GT.0D0) Then
c              Call updQb(n, lE, iEs, XYP, IPEs, qEs)
c           End if
c        End do
c
c        mBad = 0
c        Do n = 1, lE
c           If(iEs(n).GE.0 .AND. qEs(n).LE.0D0) mBad = mBad + 1
c        End do
c
c        If(mBad.GE.nBad) Goto 1000
      End if


C ... updating the grid
      flag = .TRUE.

      Call pntUpd(iP1, ICP,  XYP,  HesP,  detG,
     &                 ICPs, XYPs, HesPs, detGs)
      Call pntDel(iP2, nP, ICP, IHolP)

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


 1000 Return
 9000 Call errMes(1007, 'clpsR', 'local parameter MaxS is small')
      End Subroutine clpsR



C ================================================================
      Subroutine chkSOtet(iP1, iP2, iPa, iPb, 
     &                    XYPs, XYP, IPE, v1, flag) 
C ================================================================
C Routine checks orientation of {iP1, iPa, iPb, *} and IPE. 
C
C flag = .TRUE. if the tets are oriented correctly. 
C flag = .TRUE. if the second tet has both iP1 and iP2
C ================================================================
      Integer iP1, iP2, iPa, iPb, IPE(4)
      Real*8  XYPs(3), XYP(3, *), v1
      Logical flag

      Real*8  v2

C ================================================================
      Do i = 1, 4
         iPt = IPE(i)

         If(.NOT.check13(iPt, iP1, iPa, iPb) .AND. iPt.NE.iP2) Then
            v2 = calVol(XYPs, XYP(1, iPa), XYP(1, iPb), XYP(1, iPt))
            flag = (v1 * v2).LT.0D0 
            Goto 9000
         End if
      End do

      flag = .TRUE.
      
 9000 Return
      End Subroutine chkSOtet
C
      End Module mba3d_clpsR
