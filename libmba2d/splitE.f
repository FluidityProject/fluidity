C ====================================================
      Subroutine splitE(
C ====================================================
c group (M)
     &            iwE,
     &            nP, MaxP, nF, MaxF, nE, MaxE,
     &            XYP, IPF, IPE, lbE,
     &            parCrv, iFnc,
     &            hStar,
     &            ICP, IEP, IFE, IEE,
     &            L1E, L2E, nL2, nStep,
     &            IHolP, IHolF, IHolE,
     &            status,
C group (Q)
     &            HesP, Quality, rQuality,
     &            detG, qE,
C group (S)
     &            lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &            flag)
C ====================================================
      include 'makS.fd'
      include 'colors.fd'
      include 'operat.fd'
      include 'status.fd'
C ====================================================
C Routine splits triangle into 3 triangles by inserting
C one interior point at the center of mass.
C
C *** Remarks:
C        1. The mesh quality may drop down after the
C           splitting.
C ====================================================
C group (M)
      Integer IPF(4, *), IPE(3, *), lbE(*)
      real  XYP(2, *)

      Integer iFnc(*)
      real  parCrv(2, *), hStar

      Integer L1E(2, *), L2E(*), nStep(4)
      Integer IHolP(*), IHolF(*), IHolE(*)

      Integer ICP(*), IEP(*)
      Integer IFE(3, *), IEE(3, *)

      Integer status

C group (Q)
      real  HesP(3, *), Quality, rQuality
      real  detG(*), qE(*)

C group (S)
      Integer iFu(*), iEu(*), IPFu(2, *), IPEu(3, *)
      real  qEu(*)

C group (W)
      Logical flag

C group (Local variables)
      Integer iFs(MaxS), iEs(MaxS), IPFs(2, MaxS), IPEs(3, MaxS)
      real  XYPs(2), HesPs(3), detGs, qEs(MaxS)

      Integer ip(4)
      Integer minClr
      Logical flagFirst, ifXnode

C ====================================================
      flag = .FALSE.

      Call copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      Call findSE(lE, iEs, iwE, nE1)

      iP1 = IPE(1, iwE)
      iP2 = IPE(2, iwE)
      iP3 = IPE(3, iwE)


c ... check for faces which can not be split
      ICPs = minClr(ICP(iP1), ICP(iP2))
      ICPs = minClr(ICPs, ICP(iP3))
      If(ifXnode(ICPs, jTnode)) goto 1000


c ... create the interior point
      ICPs = jInode
      Do i = 1, 2
         XYPs(i) = (XYP(i, iP1) + XYP(i, iP2) + XYP(i, iP3)) / 3
      End do

      Do i = 1, 3
         HesPs(i) = (HesP(i, iP1) + HesP(i, iP2) + HesP(i, iP3)) / 3
      End do

      Call calDet(HesPs, detGs)

      Call pntAdd(iPs, nP, MaxP, ICP,  XYP,  HesP,  detG, IHolP,
     &                           ICPs, XYPs, HesPs, detGs)

c ... create 3 elements
      lEold = lE
      Do i1 = 1, 3
         i2 = ip(i1 + 1)

         iP1 = IPE(i1, iwE)
         iP2 = IPE(i2, iwE)

         If(i1.EQ.1) Then
            kE = nE1
         Else
            lE = lE + 1
            If(lE.GT.MaxS) goto 9000
            kE = lE
         End if

c        Call calQE(
c    &        HesP(1, iP1), detG(iP1), XYP(1, iP1),
c    &        HesP(1, iPa), detG(iPa), XYP(1, iPa),
c    &        HesPs, detGs, XYPs,
c    &        hStar, qEs(kE))

         IPEs(1, kE) = iP1 
         IPEs(2, kE) = iP2
         IPEs(3, kE) = iPs 
      End do


c ... update the quality
c     If(ifXnode(status, ANISmoothMesh)) Then
c        flagFirst = .FALSE.
c        Call updQE(XYP, lE, iEs, IPEs,
c    &              HesP, rQuality, detG, hStar, qEs, flagFirst)
c     End if 


C ... update the grid
      flag = .TRUE.

C!!!  next line simulates lstAdd
      qE(iEs(nE1)) = qEs(nE1)
      Call eleDel(iEs(nE1), IPE, IEE)

      Do n = lEold + 1, lE
         Call eleAdd(nE, MaxE, IHolE)
C!!!     3 next lines simulate lstAdd
         nE = nE + 1
         iEs(n) = nE 
         qE(iEs(n)) = qEs(n)
         Call eleDel(iEs(n), IPE, IEE)
      End do

      Call eleUpd(nE1, IEP, IPE, IFE, IEE,
     &            lF, lE, iFs, iEs, IPFs, IPEs)
      Do n = lEold + 1, lE
         Call eleUpd(n, IEP, IPE, IFE, IEE,
     &               lF, lE, iFs, iEs, IPFs, IPEs)

         lbE(iEs(n)) = lbE(iwE)
      End do

 1000 Return

 9000 Call errMes(1007, 'splitE', 'local parameter MaxS is small')
      End


