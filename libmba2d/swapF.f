C ========================================================
      Subroutine swapF(
C ========================================================
c group (M)
     &            iwF, iwE,
     &            nE, XYP, IPE,
     &            hStar,
     &            ICP, IEP, IFE, IEE,
     &            L1E, L2E, nL2, nStep,
     &            status,
C group (Q)
     &            HesP, rQuality, detG, qE,
C group (S)
     &            lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &            flag)
C ========================================================
      include 'makS.fd'
      include 'colors.fd'
      include 'status.fd'
      include 'operat.fd'
C ========================================================
C group (M)
      Integer IPE(3, *)
      real  XYP(2, *)

      real  hStar

      Integer L1E(2, *), L2E(*), nStep(4)

      Integer ICP(*), IEP(*)
      Integer IFE(3, *), IEE(3, *)

      Integer status

C group (Q)
      real  HesP(3, *), rQuality
      real  detG(*), qE(*)

C group (S)
      Integer iFu(*), iEu(*), IPFu(2, *), IPEu(3, *)
      real  qEu(*)

C group (Flag)
      Logical flag

C group (Local variables)
      Integer iFs(MaxS), iEs(MaxS), IPFs(2, MaxS), IPEs(3, MaxS)
      real  qEs(MaxS)

      Integer ip(4), iPs(MaxS)
      real  t1, t2
      Logical flagBNDs, flagFirst

      Integer minClr
      Logical ifXnode

C ========================================================
      flag = .FALSE.

      Call copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      iF = IFE(iwF, iwE)

      i1 = iwF
      i2 = ip(i1 + 1)

      iP1 = IPE(i1, iwE)
      iP2 = IPE(i2, iwE)

      iE1 = iwE
      iE2 = IEE(iwF, iE1)


c ... checking the case when swapping is impossible
      If(iE2.EQ.0) goto 1000
      If(iF.NE.0)  goto 1000

      ICPs = minClr(ICP(iP1), ICP(iP2))
      If(ifXnode(ICPs, jTnode)) goto 1000


      i3 = ip(i2 + 1)
      iPa = IPE(i3, iwE)

      Do i = 1, 3
         If(IPE(i, iE2).NE.iP1 .AND. IPE(i, iE2).NE.iP2) Then
            i4 = i
            iPb = IPE(i, iE2)
         End if
      End do


C ... checking for boundary triangles
      If(ifXnode(status, ANIForbidBoundaryElements)) Then
         If(ifXnode(ICP(iPa), jInode)) goto 100
         If(ifXnode(ICP(iPb), jInode)) goto 100

         If(ifXnode(ICP(iP1), jBnode)) goto 1000
         If(ifXnode(ICP(iP2), jBnode)) goto 1000

 100     Continue
      End if


c ... checking the case when swapping is impossible
      t1 =(XYP(1, iPa) - XYP(1, iP1)) * (XYP(2, iPb) - XYP(2, iP1)) -
     &    (XYP(2, iPa) - XYP(2, iP1)) * (XYP(1, iPb) - XYP(1, iP1))

      t2 =(XYP(1, iPa) - XYP(1, iP2)) * (XYP(2, iPb) - XYP(2, iP2)) -
     &    (XYP(2, iPa) - XYP(2, iP2)) * (XYP(1, iPb) - XYP(1, iP2))

      If(t1 * t2.GE.0D0) goto 1000


C ... skipping the inverted elements
      if(ifXnode(status, ANIUntangleMesh)) Then
         If(qE(iE1).LE.0D0) goto 1000
         If(qE(iE2).LE.0D0) goto 1000
      End if


c ... checking for surrounding points
      If(ifXnode(status, ANIUse2ArmRule)) Then
         If(ifXnode(ICP(iPa), jBnode) .OR. 
     &      ifXnode(ICP(iPb), jBnode)) Then
            Call chkSPf(iP1, iP2, ICP, IEP, IPE, IEE, lP, iPs)
            Call chkSPb(iP1, iP2, iPa, iPb, iSWAP,
     &                  ICP, IEP, IPE, IEE, lP, iPs, flagBNDs)
            If(flagBNDs) goto 1000
         End if
      End if


C ... making a virtual evaluation of the quality
      Call findSE(lE, iEs, iE1, nE1)
      Call calQE(
     &     HesP(1, iP1), XYP(1, iP1),
     &     HesP(1, iPa), XYP(1, iPa),
     &     HesP(1, iPb), XYP(1, iPb),
     &     hStar, qEs(nE1))

      If(qEs(nE1).LE.rQuality) goto 1000

      Call findSE(lE, iEs, iE2, nE2)
      Call calQE(
     &     HesP(1, iP2), XYP(1, iP2),
     &     HesP(1, iPa), XYP(1, iPa),
     &     HesP(1, iPb), XYP(1, iPb),
     &     hStar, qEs(nE2))

      If(qEs(nE2).LE.rQuality) goto 1000


C ... updating the grid
      flag = .TRUE.

      IPEs(1, nE1) = iP1
      IPEs(2, nE1) = iPa
      IPEs(3, nE1) = iPb

      IPEs(1, nE2) = iP2
      IPEs(2, nE2) = iPa
      IPEs(3, nE2) = iPb


c ... updating the quality
      If(ifXnode(status, ANISmoothMesh)) Then
         flagFirst = .TRUE.
         Call updQE(XYP, lE, iEs, IPEs,
     &              HesP, rQuality, detG, hStar, qEs, flagFirst)
         If(.NOT.flagFirst) Then
            flag = .FALSE.
            goto 1000
         End if
      End if 


      Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(nE1), qEs(nE1))
      Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(nE2), qEs(nE2))

      Call eleDel(iEs(nE1), IPE, IEE)
      Call eleDel(iEs(nE2), IPE, IEE)

      Call eleUpd(nE1, IEP, IPE, IFE, IEE,
     &            lF, lE, iFs, iEs, IPFs, IPEs)
      Call eleUpd(nE2, IEP, IPE, IFE, IEE,
     &            lF, lE, iFs, iEs, IPFs, IPEs)

 1000 Return
      End


