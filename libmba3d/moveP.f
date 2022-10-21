      Module mba3d_moveP
C
      use mba3d_auxSF
      use mba3d_auxSP
      use mba3d_list
      use mba3d_makQ
      use mba3d_minim
      use mba3d_nlnfnc
      use mba3d_update
C
      contains
C
C ================================================================
      Subroutine moveP(
C ================================================================
c group (M)
     &           iwP, iwE,
     &           nE, XYP, IPF, IPE,
     &           hStar,
     &           ICP, IFE, 
     &           L1E, L2E, nL2, nStep,
     &           status,
C group (Q)
     &           HesP, rQuality, detG, qE,
     &           MetricFunction, flagAnalytic,
C group (S)
C!!! &           lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &           lFu, lEu, iFu, iEu, IPFs, IPEs, qEu,
C group (W)
     &           nPw, nEw, XYPw, HesPw, IPEw,
     &           miLINTRP, mrLINTRP, iSE, rSE, iControl, rMove,
     &           flag)
C ================================================================
C Remark: the routine was optimized (reference !!!):
C         IPFu => IPFs  since IPFu is not changed
C         IPEu => IPEs  since IPEu is not changed
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
      Integer ICP(*), IFE(4, *)

      Integer status

C group (Q)
      Real*8  HesP(6, *), rQuality
      Real*8  detG(*), qE(*)
      Logical flagAnalytic

      Integer  MetricFunction
      EXTERNAL MetricFunction

C group (S)
C!!!  Integer iFu(*), iEu(*), IPFu(4, *), IPEu(5, *)
      Integer iFu(*), iEu(*), IPFs(4, *), IPEs(5, *)
      Real*8  qEu(*)

C group (W)
      Integer IPEw(4, *), iSE(*)
      Real*8  XYPw(3, *), HesPw(6, *)
      Real*8  rSE(*), rMove

C group (Flag)
      Logical flag

C ================================================================
C group (Local Functions)

C group (Local variables)
      Integer ip(5), iref(4), iPt(3)
      Integer iPs(MaxS), iRs(MaxS), iFs(MaxS), iCs(MaxS), iEs(MaxS)

      Real*8  XYPt(3), XYPs(3, 3)
      Real*8  HesPs(6), detGs, qEs(MaxS)

      Logical flagTM

C ... for nonlinear minimization procedure
      Real*8  Z(3, 3), U(3), U1(3)

      Real*8  q(4), qMin, hMin, hMax
      Real*8  xStep(3), t1, t2
      Real*8  heit, v, s

C ================================================================
      flag = .FALSE.

      Call copySEmove(lFu, lEu, iFu, iEu, qEu,
     &                lF,  lE,  iFs, iEs, qEs)

      iref(1) = 1
      iref(2) = 2
      iref(3) = 3
      iref(4) = 1

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1

      iP1 = IPE(iwP, iwE)
      ICPs = ICP(iP1)


C ... checking the case when moving is impossible
      If(ifXnode(ICPs, jVnode)) Goto 1000
      If(ifXnode(ICPs, jTnode)) Goto 1000


C ... checking the number of inverted elements
      flagTM = ifXnode(status, ANITangledMesh)
      nBad = 0
      if(flagTM) Then
         Do n = 1, lE
            If(qEs(n).LE.0D0) nBad = nBad + 1
         End do
      End if

      flagTM = flagTM .AND. nBad.GT.0


      heit = 1D12
      Do 10 n = 1, lE
         iE = iEs(n)

         Do i1 = 1, 4
            If(IPE(i1, iE).EQ.iP1) Then
               i2 = ip(i1 + 1)
               i3 = ip(i2 + 1)
               i4 = ip(i3 + 1)

               iP2 = IPE(i2, iE)
               iP3 = IPE(i3, iE)
               iP4 = IPE(i4, iE)

               v = calVol(XYP(1, iP1), XYP(1, iP2),
     &                    XYP(1, iP3), XYP(1, iP4))
               s = calSqr(XYP(1, iP2), XYP(1, iP3), XYP(1, iP4))
               heit = min(heit, 6D0 * dabs(v) / s)

               Goto 10
            End if
         End do

         If(.NOT.flagTM) iEs(n) = -iE
 10   Continue


C ... calculating a search direction
      lC1 = 0

      If(ifXnode(ICPs, jRnode)) Then
         nZ = 1
         Do 20 n = 1, lF
            iPt(1) = IPFs(1, n)
            iPt(2) = IPFs(2, n)
            iPt(3) = IPFs(3, n)

            Do i1 = 1, 3
               If(iP1.EQ.iPt(i1)) Then
                  lC1 = lC1 + 1
                  iCs(lC1) = n
                  iPs(lC1) = -1
                  iRs(lC1) = IPF(4, iFs(n))

                  i2 = i1
                  Do j = 1, 2
                     i2 = iref(i2 + 1)

                     Call clrSR(iP1, iPt(i2), ICP, IPF, IFE,
     &                          lF, iFs, lE, iEs, ICRab)

                     If(ICRab.EQ.ICPs) Then
                        iPs(lC1) = iPt(i2)

                        Do i = 1, 3
                           Z(i, 1) = XYP(i, iPt(i2)) - XYP(i, iP1)
                        End do

                        If(iPs(lC1).EQ.iPs(1)) Goto 20
                     End if
                  End do
               End if
            End do
 20      Continue

      Else If(ifXnode(ICPs, jSnode)) Then
         nZ = 2
         Do n = 1, lF
            iPt(1) = IPFs(1, n)
            iPt(2) = IPFs(2, n)
            iPt(3) = IPFs(3, n)

            Do i1 = 1, 3
               If(iP1.EQ.iPt(i1)) Then
                  i2 = iref(i1 + 1)
                  i3 = iref(i2 + 1)

                  Do i = 1, 3
                     Z(i, 1) = XYP(i, iPt(i2)) - XYP(i, iP1)
                     Z(i, 2) = XYP(i, iPt(i3)) - XYP(i, iP1)
                  End do

                  Call VecMul(Z(1, 1), Z(1, 2), Z(1, 3))
                  Call VecMul(Z(1, 1), Z(1, 3), Z(1, 2))
               End if
            End do
         End do

      Else If(ifXnode(ICPs, jInode)) Then
         nZ = 3
         Do i = 1, 3
            Do j = 1, 3
               Z(i, j) = 0D0
            End do
            Z(i, i) = 1D0
         End do
      End if


C ... calculating the terminal points in the search directions
      Do i = 1, nZ
         If(nZ.EQ.1) Then
            xStep(i) = heit

            Do j = 1, 3
               XYPs(j, i) = XYP(j, iP1) + xStep(i) * Z(j, i)
            End do
         Else If(nZ.EQ.2) Then
            t1 = distSR(iP1,  Z(1,i),  Z(2,i),  Z(3,i),
     &                  lF, iFs, IPFs, XYP)
            t2 = distSR(iP1, -Z(1,i), -Z(2,i), -Z(3,i),
     &                  lF, iFs, IPFs, XYP)

            If(t2.GT.t1) t1 = -t2
            xStep(i) = 2D-2 * min(dabs(t1), heit) * dsign(1D0, t1)

            Do j = 1, 3
               XYPs(j, i) = XYP(j, iP1) + xStep(i) * Z(j, i)
            End do

         Else If(nZ.EQ.3) Then
            t1 = distSF(iP1,  Z(1,i),  Z(2,i),  Z(3,i),
     &                  lE, iEs, IPEs, XYP)
            t2 = distSF(iP1, -Z(1,i), -Z(2,i), -Z(3,i),
     &                 lE, iEs, IPEs, XYP)

C  ...  the case when shut is along a face
            If(t1.LT.0D0 .AND. t2.LT.0D0) Goto 1000

            If(t2.GT.t1) t1 = -t2
            xStep(i) = 2D-2 * min(dabs(t1), heit) * dsign(1D0, t1)

            Do j = 1, 3
               XYPs(j, i) = XYP(j, iP1)
            End do

            XYPs(i, i) = XYP(i, iP1) + xStep(i)
         End if
      End do



      q(1) = 1D0 - rQuality
      Do i = 1, nZ
         q(i + 1) = NLnFnc(
C group (F)
     &       3, XYPs(1, i),
C group (ANI)
     &       XYP, HesP, detG, hStar,
     &       iP1, lE, iEs, XYPt, ICPs, IPEs, HesPs, detGs, qEs,
     &       nPw, nEw, XYPw, HesPw, IPEw,
     &       MetricFunction, flagAnalytic,
     &       miLINTRP, mrLINTRP, iSE, rSE, iControl)

C  ...   updating the spoiled values of qEs
         Do n = 1, lE
            qEs(n) = qEu(n)
         End do
      End do



      nU = 3
      Do i = 1, 3
         U(i) = XYP(i, iP1)
      End do


      hMin = 0D0
      If(ifXnode(ICPs, jRnode)) Then
         If(q(2).GT.q(1)) Then
            Do i = 1, 3
               Z(i, 1) = -Z(i, 1)
            End do
         End if

         hMax = heit
      Else If(ifXnode(ICPs, jSnode)) Then
         t1 = (q(1) - q(2)) / xStep(1)
         t2 = (q(1) - q(3)) / xStep(2)

         Do i = 1, 3
            Z(i, 1) = t1 * Z(i, 1) + t2 * Z(i, 2)
         End do

         hMax = distSR(iP1, Z(1,1), Z(2,1), Z(3,1),
     &                 lF, iFs, IPFs, XYP)
         If(hMax.LT.0D0) Goto 1000

         hMax = heit

      Else If(ifXnode(ICPs, jInode)) Then
         Do i = 1, 3
            Z(i, 1) = (q(1) - q(i + 1)) / xStep(i)
         End do

         hMax = distSF(iP1, Z(1,1), Z(2,1), Z(3,1), lE, iEs, IPEs, XYP)
         If(hMax.LT.0D0) Goto 1000
      End if


      hMax = hMax - 0.2 * (hMax - hMin)
      If(flagTM) hMax = 2 * hMax 


      qMin = q(1)
      Call minim(
C group(F)
     &     nU, U, Z, hMin, hMax, qMin, U1,
     &     icnt, rMove, .FALSE., flag,
C group (ANI)
     &     XYP, IPE, detG, HesP, hStar,
     &     iP1, lE, iEs, XYPs, ICPs, IPEs, detGs, HesPs, qEs,
     &     nPw, nEw, XYPw, HesPw, IPEw,
     &     MetricFunction, flagAnalytic,
     &     miLINTRP, mrLINTRP, iSE, rSE, iControl,
     &     flagTM, nBad)
       If(rMove.EQ.0D0) Goto 1000


C ... analysing information from the previous routine
      m = 1
      If(.NOT.flag .AND. ICPs.EQ.jInode) Then
         If(q(2).LT.q(1)) Then
            m = 1
         Else If(q(3).LT.q(1)) Then
            m = 2
         Else If(q(4).LT.q(1)) Then
            m = 3
         Else
            Goto 1000
         End if


         qMin = NLnFnc(
C group (F)
     &       3, XYPs(1, m),
C group (ANI)
     &       XYP, HesP, detG, hStar,
     &       iP1, lE, iEs, XYPt, ICPs, IPEs, HesPs, detGs, qEs,
     &       nPw, nEw, XYPw, HesPw, IPEw,
     &       MetricFunction, flagAnalytic,
     &       miLINTRP, mrLINTRP, iSE, rSE, iControl)
      Else If(.NOT.flag) Then
         Goto 1000
      End if



C ... analyzing curvilinear and plane faces
c     Call copySQ(0, qEs, XYP(1, iP1), HesPs, detGs,
c    &               qEt, XYPt,        HesPt, detGt)

      Call pntUpd(iP1,  ICP,  XYP,  HesP,  detG,
     &            ICPs, XYPs(1, m), HesPs, detGs)


C ... checking for inverted elements
      If(flagTM) Then
         Do n = 1, lE
            If(iEs(n).GE.0) Then
               Call updQb(n, lE, iEs, XYP, IPEs, qEs)
            End if
         End do
      End if


C ... updating the grid
      flag = .TRUE.

      Do 200 n = 1, lE
         If(iEs(n).LE.0) Goto 200
         Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(n), qEs(n))
 200  Continue

 1000 Return
      End Subroutine moveP



C ================================================================
      Real*8 Function distSR(iPo, nx, ny, nz, lF, iFs, IPFs, XYP)
C ================================================================
      Real*8    PREC
      Parameter(PREC = 1D-10)

      Integer iFs(*), IPFs(4, *)
      Real*8  XYP(3, *), nx, ny, nz

C ================================================================
C group (Local variables)
      Integer iref(4)
      Real*8  XYC(3), XYD(3), VecA(3), VecB(3)
      Real*8  v1, v2

C ================================================================
      iref(1) = 1
      iref(2) = 2
      iref(3) = 3
      iref(4) = 1

      distSR = 1D12
      Do 20 n = 1, lF
         If(iFs(n).LE.0) Goto 20
         Do i1 = 1, 3
            If(IPFs(i1, n).EQ.iPo) Then
               i2 = iref(i1 + 1)
               i3 = iref(i2 + 1)

               iPa = IPFs(i2, n)
               iPb = IPFs(i3, n)

               XYC(1) = XYP(1, iPo) + nx
               XYC(2) = XYP(2, iPo) + ny
               XYC(3) = XYP(3, iPo) + nz

               v1 = calVol(XYP(1, iPo), XYP(1, iPa), XYP(1, iPb), XYC)
               If(dabs(v1).GT.PREC) Goto 20


               Do i = 1, 3
                  VecA(i) = XYP(i, iPa) - XYP(i, iPo)
                  VecB(i) = XYP(i, iPb) - XYP(i, iPo)
               End do

               Call VecMul(VecA, VecB, XYD)

               Do i = 1, 3
                  XYD(i) = XYP(i, iPo) + XYD(i)
               End do

               v1 = calVol(XYP(1, iPo), XYD, XYC, XYP(1, iPa))
               v2 = calVol(XYP(1, iPo), XYD, XYC, XYP(1, iPb))
               If(v1 * v2.GT.0D0) Goto 20

               distSR = 1D0
               Goto 1000
            End if
         End do
 20   Continue

      distSR = -distSR
 1000 Return
      End Function distSR



C ================================================================
      Real*8 Function distSF(iPo, nx, ny, nz, lE, iEs, IPEs, XYP)
C ================================================================
      include 'shutF.fd'
C ================================================================
      Integer iEs(*), IPEs(5, *)
      Real*8  XYP(3, *), nx, ny, nz

C group (Local variables)
      Integer iref(5)
      Real*8  XYD(3), d, h
      Real*8  v, v1, v2, v3
      Logical flagSHUT

C ================================================================
      iref(1) = 1
      iref(2) = 2
      iref(3) = 3
      iref(4) = 4
      iref(5) = 1

      d = dsqrt(nx ** 2 + ny ** 2 + nz ** 2)

      distSF = 1D12
      Do 20 n = 1, lE
         If(iEs(n).LE.0) Goto 20
         Do i1 = 1, 4
            If(IPEs(i1, n).EQ.iPo) Then
               i2 = iref(i1 + 1)
               i3 = iref(i2 + 1)
               i4 = iref(i3 + 1)

               iPa = IPEs(i2, n)
               iPb = IPEs(i3, n)
               iPc = IPEs(i4, n)

               XYD(1) = XYP(1, iPo) + nx
               XYD(2) = XYP(2, iPo) + ny
               XYD(3) = XYP(3, iPo) + nz

               flagSHUT = lOpenTriag
               Call shutF(XYP(1, iPo), XYD,
     &                    XYP(1, iPa), XYP(1, iPb), XYP(1, iPc),
     &                    flagSHUT)
               If(.NOT.flagSHUT) Goto 20

               v  = calVol(XYP(1, iPo), XYP(1, iPa),
     &                     XYP(1, iPb), XYP(1, iPc))
               v1 = calVol(XYP(1, iPo), XYD, XYP(1, iPa), XYP(1, iPb))
               v2 = calVol(XYP(1, iPo), XYD, XYP(1, iPa), XYP(1, iPc))
               v3 = calVol(XYP(1, iPo), XYD, XYP(1, iPb), XYP(1, iPc))

               h = dabs(v) / (dabs(v1) + dabs(v2) + dabs(v3))

c  ...  checking for the orientation of points O and D
               XYD(1) = XYP(1, iPo) + 2 * h * nx
               XYD(2) = XYP(2, iPo) + 2 * h * ny
               XYD(3) = XYP(3, iPo) + 2 * h * nz

               v1 = calVol(XYD, XYP(1, iPa), XYP(1, iPb), XYP(1, iPc))
               If(v * v1.LT.0D0) Then
                  distSF = h * d
                  Goto 1000
               End if
            End if
         End do
 20   Continue

      distSF = -distSF
 1000 Return
      End Function distSF
C
      End Module mba3d_moveP
