C ======================================================
      Subroutine moveP(
C ======================================================
c group (M)
     &            iwP, iwE,
     &            nE,
     &            XYP, IPF, IPE,
     &            calCrv, parCrv, iFnc,
     &            hStar,
     &            ICP, IEP, IEE,
     &            L1E, L2E, nL2, nStep,
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
     &            iSE, rSE,
     &            icnt, rMove, flag)
C ======================================================
      include 'makS.fd'
      include 'colors.fd'
      include 'status.fd'
      include 'magic.fd'
C ======================================================
C Routine moves a vertex of triangle iwE to increase the
C quality of mesh element in the superelement associated
C with the vertex.
C ======================================================
C group (M)
      Integer IPF(4, *), IPE(3, *)
      real  XYP(2, *)

      EXTERNAL calCrv
      Integer iFnc(*)
      real  parCrv(2, *), hStar

      Integer L1E(2, *), L2E(*), nStep(4)

      Integer ICP(*), IEP(*)
      Integer IEE(3, *)

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
      Integer iFu(*), iEu(*), IPFu(2, *), IPEu(3, *)
      real  qEu(*)

C group (W)
      Integer IPEw(3, *), iSE(*)
      real  XYPw(2, *), HesPw(3, *)
      real  rSE(*), rMove

C group (Flag)
      Logical flag

C group (Local variables)
      Integer iFs(MaxS), iEs(MaxS), IPFs(2, MaxS), IPEs(3, MaxS)
      real  XYPs(2, 2), HesPs(3), detGs, qEs(MaxS)
      real  prjXYPs(2, 2)

C ... for nonlinear minimization procedure
      real  NLnFnc, ZZ(2), U(2), U1(2)

      Integer ip(4), iCRVs(2), iFNCs(2)
      real  distSP, calVol, calEdge, heit, v, d
      real  par(4), q(3), qMin, hMin, hMax, xStep, yStep, tStep
      real  tc(2), t1, t2, t3

      Logical ifXnode, flagTM

C ======================================================
      flag = .FALSE.

      Call copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      iP1 = IPE(iwP, iwE)
      ICPs = ICP(iP1)


C ... checking the case when moving is impossible
      If(ifXnode(ICPs, jVnode)) goto 1000
      If(ifXnode(ICPs, jTnode)) goto 1000

      Do i = 1, 3
         HesPs(i) = HesP(i, iP1)
      End do


c ... removing elements from superelement
      heit = 1D12
      Do 10 n = 1, lE
         iE = iEs(n)

         Do i1 = 1, 3
            If(IPE(i1, iE).EQ.iP1) Then
               i2 = ip(i1 + 1) 
               i3 = ip(i2 + 1)

               iP2 = IPE(i2, iE)
               iP3 = IPE(i3, iE)

               v = calVol( XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))
               d = calEdge(XYP(1, iP2), XYP(1, iP3))
               heit = min(heit, 2D0 * abs(v) / d)

               goto 10
            End if
         End do

         If(.NOT.ifXnode(status, ANISmoothMesh)) iEs(n) = 0
 10   Continue

      
C ... checking for inverted elements
      flagTM = ifXnode(status, ANIUntangleMesh) 
      if(flagTM) Then
         nBad = 0
         Do n = 1, lE
            If(qEs(n).LE.0D0) nBad = nBad + 1
         End do

         If(heit.LE.0D0) goto 1000
      End if

      if(flagTM) Then
         flagTM = nBad.GT.0
      End if


C ... computing the gradient
      xStep = distSP(iP1,  1D0, 0D0, lE, IPEs, XYP)
      tStep = distSP(iP1, -1D0, 0D0, lE, IPEs, XYP)
      If(tStep.GT.xStep) xStep = -tStep
      xStep = DiscreteGrad * min(abs(xStep), heit) * sign(1.0, xStep)

      yStep = distSP(iP1, 0D0,  1D0, lE, IPEs, XYP)
      tStep = distSP(iP1, 0D0, -1D0, lE, IPEs, XYP)
      If(tStep.GT.yStep) yStep = -tStep
      yStep = DiscreteGrad * min(abs(yStep), heit) * sign(1.0, yStep)


      If(ICPs.EQ.jInode) Then
         XYPs(1, 1) = XYP(1, iP1) + xStep
         XYPs(2, 1) = XYP(2, iP1)
         XYPs(1, 2) = XYP(1, iP1)
         XYPs(2, 2) = XYP(2, iP1) + yStep

         Do i = 1, 2
            iCRVs(i) = 0
            iFNCs(i) = 0
         End do
      Else if(ifXnode(ICPs, jSnode)) Then
         Call infoP(iP1, iPa, iPb, iF1, iF2, par,
     &              IPF, parCrv, lF, iFs)
         Do i = 1, 2
            XYPs(i, 1) = XYP(i, iP1) +
     &                   DiscreteGrad * (XYP(i, iPa) - XYP(i, iP1))
            XYPs(i, 2) = XYP(i, iP1) +
     &                   DiscreteGrad * (XYP(i, iPb) - XYP(i, iP1))
         End do

         iCRVs(1) = iPF(3, iF1)
         iCRVs(2) = iPF(3, iF2)

         iFNCs(1) = iFnc(iF1)
         iFNCs(2) = iFnc(iF2)

         tc(1) = par(2) + DiscreteGrad * (par(1) - par(2))
         tc(2) = par(3) + DiscreteGrad * (par(4) - par(3))

         Do i = 1, 2
            Call findSE(nCrvFnc, LFnc, iFNCs(i), k)
            If(k.GT.0) Then
               ir = ILt(k)
               Call prjCrv(XYPs(1, i), prjXYPs(1, i), iFNCs(i), tc(i),
     &                                                          calCrv,
     &                     L1Et(1, ir), L2Et(ir), nL2t(k), 
     &                     nStept(1, k), nEt(k), tE(ir))
               Do j = 1, 2
                  XYPs(j, i) = prjXYPs(j, i)
               End do
            End if
         End do
      End if


c ... calculating the search direction
      q(1) = 1D0 - rQuality
      Do i = 1, 2
         q(i + 1) = NLnFnc(
C group (F)
     &       2, XYPs(1, i),
C group (ANI)
     &       XYP, IPE, IEE, HesP, hStar, status,
     &       lE, iEs, prjXYPs, IPEs, HesPs, detGs, qEs,
     &       nPw, nEw, XYPw, HesPw, IPEw,
     &       MetricFunction, flagAnalytic,
     &       iSE, rSE, iP1, 0, calCrv,
     &       L1Et, L2Et, tE,
     &       nL2t, nStept, nEt, nCrvFnc, LFnc, ILt)


C  ...   updating the spoiled values of qEs
         Do n = 1, lE
            qEs(n) = qEu(n)
         End do
      End do


      nU = 2
      hMin = 0D0
      iref = 1

      Do i = 1, 2
         U(i) = XYP(i, iP1)
      End do
      If(ICPs.EQ.jInode) Then
         ZZ(1) = (q(1) - q(2)) / xStep
         ZZ(2) = (q(1) - q(3)) / yStep
         hMax = distSP(iP1, ZZ(1), ZZ(2), lE, IPEs, XYP)

c  ...   check that the direction is not zero
         If(ZZ(1) ** 2 + ZZ(2) ** 2.LE.0D0) goto 1000
         If(hMax.LE.0D0)                    goto 1000

      Else if(ifXnode(ICPs, jSnode)) Then
         If(q(2).LE.q(1)) Then
            If(iCRVs(1).NE.0) Then
               nU = 1

               iFmove = iF1
               iPend = iPa
               hMin = min(par(1), par(2))
               hMax = max(par(1), par(2))
               iref = 1

               t1 = par(2)
               t2 = par(1)
               ZZ(1) = t2 - t1

               U(1) = t1
            Else
               Do i = 1, 2
                  ZZ(i) = XYP(i, iPa) - XYP(i, iP1)
               End do
               hMax = sqrt(ZZ(1) ** 2 + ZZ(2) ** 2)
            End if
         Else If(q(3).LE.q(1)) Then
            If(iCRVs(2).NE.0) Then
               nU = 1

               iFmove = iF2
               iPend = iPb
               hMin = min(par(3), par(4))
               hMax = max(par(3), par(4))
               iref = 2

               t1 = par(3)
               t2 = par(4)
               ZZ(1) = t2 - t1

               U(1) = t1
            Else
               Do i = 1, 2
                  ZZ(i) = XYP(i, iPb) - XYP(i, iP1)
               End do
               hMax = sqrt(ZZ(1) ** 2 + ZZ(2) ** 2)
            End if
         Else
            goto 1000
         End if
      End if

      hMax = hMax - 0.2 * (hMax - hMin)
      If(flagTM) hMax = 2 * hMax


C ... minimizing the functional
      qMin = q(1)

      Call minim(
C group(F)
     &     nU, U, ZZ, hMin, hMax, qMin, U1,
     &     icnt, rMove, flag,
C group (ANI)
     &     XYP, IPE, IEE, HesP, hStar, status,
     &     lE, iEs, XYPs, IPEs, detGs, HesPs, qEs,
     &     nPw, nEw, XYPw, HesPw, IPEw, 
     &     MetricFunction, flagAnalytic,
     &     iSE, rSE, iP1, iFNCs(iref), calCrv,
     &     L1Et, L2Et, tE,
     &     nL2t, nStept, nEt, nCrvFnc, LFnc, ILt,
C group (Tangle)
     &     flagTM, nBad)
      If(rMove.EQ.0D0) goto 1000


C ... analysing information from the previous routine
      iref = 1
      If(.NOT.flag) Then
         If(nU.EQ.1) goto 1000

         If(q(2).LT.q(1)) Then
            iref = 1
         Else If(q(3).LT.q(1)) Then
            iref = 2
         Else
            goto 1000
         End if

         qMin = NLnFnc(
C group (F)
     &       2, XYPs(1, iref),
C group (ANI)
     &       XYP, IPE, IEE, HesP, hStar, status,
     &       lE, iEs, prjXYPs, IPEs, HesPs, detGs, qEs,
     &       nPw, nEw, XYPw, HesPw, IPEw,
     &       MetricFunction, flagAnalytic,
     &       iSE, rSE, iP1, 0, calCrv,
     &       L1Et, L2Et, tE,
     &       nL2t, nStept, nEt, nCrvFnc, LFnc, ILt)
      End if


C ... updating the grid
      flag = .TRUE.
      Call pntUpd(iP1, ICP,  XYP,  HesP,  detG,
     &                 ICPs, XYPs(1, iref), HesPs, detGs)

C ... updating for inverted elements
      If(flagTM) Then
         Do n = 1, lE
            If(iEs(n).GT.0) Then
               Call updQb(n, lE, iEs, XYP, IPEs, qEs)
            End if
         End do
      End if

      If(ifXnode(ICPs, jSnode)) Then
         Call findSE(lF, iFs, iF1, nF1)
         Call findSE(lF, iFs, iF2, nF2)

         IPFs(1, nF1) = iPa
         IPFs(2, nF1) = iP1

         IPFs(1, nF2) = iP1
         IPFs(2, nF2) = iPb

         iBNDs = IPF(4, iF1)

         t1 = par(1)
         t3 = U(1)
         t2 = par(4)

         Call facUpd(nF1, IPF, parCrv, iFnc,
     &        iFs, IPFs, iCRVs(1), iFNCs(1), iBNDs, t1, t3)

         Call facUpd(nF2, IPF, parCrv, iFnc,
     &        iFs, IPFs, iCRVs(2), iFNCs(2), iBNDs, t3, t2)
      End if

      Do 20 n = 1, lE
         If(iEs(n).LE.0) goto 20
         Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(n), qEs(n))
 20   Continue

 1000 Return
      End



C ======================================================
      real Function distSP(iPo, nx, ny, lE, IPEs, XYP)
C ======================================================
      Integer IPEs(3, *)
      real  XYP(2, *), nx, ny

C group (Local variables)
      Integer iref(4)
      real  ox, oy, ax, ay, bx, by, dt, dta, dtb, dtc, dx, dy
      real  cx, cy, aNorm, bNorm, bet, gam
      Logical flag

C ======================================================
      flag = .FALSE.

      iref(1) = 1
      iref(2) = 2
      iref(3) = 3
      iref(4) = 1

      distSP = 1D24
      Do 20 n = 1, lE
         Do i1 = 1, 3
            If(IPEs(i1, n).EQ.iPo) Then
               i2 = iref(i1 + 1)
               i3 = iref(i2 + 1)

               iPa = IPEs(i2, n)
               iPb = IPEs(i3, n)

               ox = XYP(1, iPo)
               oy = XYP(2, iPo)

               ax = XYP(1, iPa) - ox
               ay = XYP(2, iPa) - oy

               bx = XYP(1, iPb) - ox
               by = XYP(2, iPb) - oy
               goto 10
            End if
         End do
         goto 20

c  ...  calculating the bisectris between (ax, ay) and (bx, by)
 10      aNorm = sqrt(ax ** 2 + ay ** 2)
         bNorm = sqrt(bx ** 2 + by ** 2)
         cx = ax / aNorm + bx / bNorm
         cy = ay / aNorm + by / bNorm

         dt  = ax * cy - cx * ay
         dtc = nx * cy - cx * ny
         If(dtc / dt.LT.0D0) goto 15

         dta = ax * ny - nx * ay
         If(dta / dt.GE.0D0) goto 17

 15      dt = -dt
         If(dtc / dt.LT.0D0) goto 20

         dtb = bx * ny - nx * by
         If(dtb / dt.LT.0D0) goto 20

 17      Continue
         If(abs(ny).LE.1D-12) Then
            dy = oy

            If(abs(ay - by).LE.1D-12) goto 20
            dx = ((bx + ox) * ay - (ax + ox) * by) / (ay - by)
         Else
            dt  = nx / ny
            bet = bx - ax + dt * (ay - by)
            gam = (by + oy) * ax - (ay  + oy) * bx +
     &            dt * oy * (by - ay)

            dy = -gam / bet
            dx = ox + dt * (dy - oy)
         End if

         distSP = min(distSP, (dx - ox) ** 2 + (dy - oy) ** 2)
         flag = .TRUE.
 20   Continue
      distSP = sqrt(distSP)

      If(.NOT.flag) distSP = -distSP
      Return
      End
