      Module mba3d_splitE
C
      use mba3d_error
      use mba3d_makQ
      use mba3d_update
C
      contains
C
C ================================================================
      Subroutine splitE(
C ================================================================
C group (M)
     &           iwE,
     &           nP, MaxP, nE, MaxE,
     &           XYP, IPE,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           IHolP, IHolE,
C group (Q)
     &           HesP, detG, qE,
C group (S)
     &           lF, lE, iFs, iEs, IPFs, IPEs, qEs,
C group (W)
     &           flag)
C ================================================================
      include 'makS.fd'
      include 'color.fd'
      include 'operat.fd'
C ================================================================
C The operation does not change the quality list, since the last
C does not work properly when the quality is descreased.
C
C The metric in the new point is interpolated, since 8-tree is not
C yet available. The latter comes from an attempt to save memory.
C ================================================================
C group (M)
      Integer IPE(5, *)
      Real*8  XYP(3, *)

      Real*8  hStar

      Integer IHolP(*), IHolE(*)

      Integer ICP(*), IEP(*)
      Integer IFE(4, *), IEE(4, *)

C group (Q)
      Real*8  HesP(6, *), detG(*), qE(*)

C group (S)
      Integer iFs(*), iEs(*), IPFs(4, *), IPEs(5, *)
      Real*8  qEs(*)

C group (Flag)
      Logical flag
C ================================================================
C group (Local variables)
      Integer ip(5)
      Real*8  HesPs(6), detGs, XYPs(3), v

C ================================================================
      flag = .TRUE.

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1


      iE1 = iwE
      Call findSE(lE, iEs, iE1, nE1)

      iP1 = IPEs(1, nE1)
      iP2 = IPEs(2, nE1)
      iP3 = IPEs(3, nE1)
      iP4 = IPEs(4, nE1)


c ... creating an inner point
      ICPs = jInode
      Do i = 1, 3
         XYPs(i) = (XYP(i, iP1) + XYP(i, iP2)
     &            + XYP(i, iP3) + XYP(i, iP4)) / 4
      End do

      Do i = 1, 6
         HesPs(i) = (HesP(i, iP1) + HesP(i, iP2)
     &             + HesP(i, iP3) + HesP(i, iP4)) / 4
      End do

c!!!  LDH = 6
c!!!  nXY = 1
c!!!  Call LINTRP3D(nEw, IPEw, nPw, XYPw, LDH, HesPw, nXY, XYPs,
c!!! &              HesPs, iSE, miLINTRP, rSE, mrLINTRP, iControl)

      Call calDet(HesPs, detGs)

      Call pntAdd(iPs, nP, MaxP, ICP,  XYP,  HesP,  detG, IHolP,
     &                           ICPs, XYPs, HesPs, detGs)

c ... creating 4 elements
      lEold = lE
      Do i1 = 1, 4
         i2 = ip(i1 + 1)
         i3 = ip(i2 + 1)

         iP1 = IPE(i1, iwE)
         iP2 = IPE(i2, iwE)
         iP3 = IPE(i3, iwE)

         If(i1.EQ.1) Then
            kE = nE1
         Else
            lE = lE + 1
            If(lE.GT.MaxS) Goto 9000
            kE = lE
         End if

         Call calQE(
     &        HesP(1, iP1), detG(iP1), XYP(1, iP1),
     &        HesP(1, iP2), detG(iP2), XYP(1, iP2),
     &        HesP(1, iP3), detG(iP3), XYP(1, iP3),
     &        HesPs, detGs, XYPs,
     &        hStar, qEs(kE), v)
         
         IPEs(1, kE) = iP1 
         IPEs(2, kE) = iP2
         IPEs(3, kE) = iP3 
         IPEs(4, kE) = iPs 
         IPEs(5, kE) = IPE(5, iE1)
      End do


C ... updating the grid
C!!!  Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(nE1), qEs(nE1))
C!!!  next line simulates lstAdd
      qE(iEs(nE1)) = qEs(nE1)
      Call eleDel(iEs(nE1), IPE, IEE)
                                                                                
      Do n = lEold + 1, lE
         Call eleAdd(nE, MaxE, IHolE)
C!!!     Call lstAdd(nE, L1E, nL2, L2E, nStep, IHolE,
C!!! &               qE, qEs(n), iEs(n))
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
      End do

      Return
 9000 Call errMes(1007, 'splitE', 'local parameter MaxS is small')
      End Subroutine splitE
C
      End Module mba3d_splitE
