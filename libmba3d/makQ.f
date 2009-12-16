      Module mba3d_makQ
C
      use mba3d_auxSE
      use mba3d_auxSP
      use mba3d_calVol
C
      contains
C
C ================================================================
      Subroutine makQ(
C ================================================================
     &     nLoop,
c group (M)
     &     nP, nE, XYP, IPE, nEv, IEV,
     &     nEStar, hStar,
c group (Q)
     &     HesP, detG, qE)
C ================================================================
      include 'magic.fd'
      include 'cubature.fd'
C ================================================================
C Quality computation for mesh elements.
C
C Pre-conditions:  1. connectivity structure {IPE(4, *), XYP(3, *)}
C                  2. tensor metric field HesP(6, *)
C
C Post-conditions: 1. determinant of the tensor metric, detG(*)
C                  2. quality of mesh elements, qE(*) 
C
C Remark: nLoop is used in parallel and simulation codes.
C         It should equal to 1 in sequential code.
C ================================================================
C group (M)
      Real*8  XYP(3, *)
      Integer IPE(4, *), IEV(*)

      Real*8  hStar

C group (Q)
      Real*8  HesP(6, *)
      Real*8  detG(*), qE(*)

C group (Local variables)
      Real*8  d1, d2, d3, d4, dsum
      Real*8  VStar, Vn

C ================================================================
      Do n = 1, nP
         Call calDet(HesP(1, n), detG(n))
      End do

      If(nLoop.EQ.1) Then
         VStar = 0D0
         Do n = 1, nE
            i1 = IPE(1, n)
            i2 = IPE(2, n)
            i3 = IPE(3, n)
            i4 = IPE(4, n)

            d1 = detG(i1)
            d2 = detG(i2)
            d3 = detG(i3)
            d4 = detG(i4)
           
            dsum = dsqrt(d1 * T2B + (d2 + d3 + d4) * T2A) 
     &           + dsqrt(d2 * T2B + (d3 + d4 + d1) * T2A) 
     &           + dsqrt(d3 * T2B + (d4 + d1 + d2) * T2A) 
     &           + dsqrt(d4 * T2B + (d1 + d2 + d3) * T2A) 

            Vn = calVol(XYP(1, i1), XYP(1, i2),
     &                  XYP(1, i3), XYP(1, i4))
            VStar = VStar + dabs(Vn) * dsum / 4
         End do
         hStar = (VStar / nEStar * magicNumber *
     &                    12D0 / dsqrt(2D0)) ** 3.333D-1
      End if

      Do n = 1, nE
         i1 = IPE(1, n)
         i2 = IPE(2, n)
         i3 = IPE(3, n)
         i4 = IPE(4, n)

         Call calQE(
     &        HesP(1, i1), detG(i1), XYP(1, i1),
     &        HesP(1, i2), detG(i2), XYP(1, i2),
     &        HesP(1, i3), detG(i3), XYP(1, i3),
     &        HesP(1, i4), detG(i4), XYP(1, i4),
     &        hStar, qE(n), Vn)
      End do


c ... set quality of fixed elements to 1
      Do n = 1, nEv
         qE(IEV(n)) = 1D0
      End do 

      Return
      End Subroutine makQ



C ================================================================
      Subroutine updQa(n, XYP, IPE, IEE, qE)
C ================================================================
C Initial quality modification for tangled elements and 
C their closest (face-) neighboors.
C ================================================================
      Real*8  XYP(3, *), qE(*)
      Integer IPE(5, *), IEE(4, *)

C (Local variables)
      Integer ip(5)
      Real*8  v1, v2

C ================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1

      Do 100 i1 = 1, 4
         iE = IEE(i1, n)
         If(iE.LE.0) Goto 100

         i2 = ip(i1 + 1)
         i3 = ip(i2 + 1)

         iP1 = IPE(i1, n)
         iP2 = IPE(i2, n)
         iP3 = IPE(i3, n)

         Do j1 = 1, 4
            j2 = ip(j1 + 1)
            j3 = ip(j2 + 1)

            jP1 = IPE(j1, iE)
            jP2 = IPE(j2, iE)
            jP3 = IPE(j3, iE)

            If(check33(iP1, iP2, iP3, jP1, jP2, jP3)) Then
               i4  = ip(i3 + 1)
               iP4 = IPE(i4, n)

               j4  = ip(j3 + 1)
               jP4 = IPE(j4, iE)

               v1 = calVol(XYP(1, iP1), XYP(1, iP2),
     &                     XYP(1, iP3), XYP(1, iP4))
               v2 = calVol(XYP(1, iP1), XYP(1, iP2),
     &                     XYP(1, iP3), XYP(1, jP4))

               If(v1 * v2.GE.0D0) Then
                  qE(n)  = -dabs(qE(n))
                  qE(iE) = -dabs(qE(iE))
               End if
               Goto 100
            End if
         End do
 100  Continue

      Return
      End Subroutine updQa


C ================================================================
      Subroutine updQb(nEs, lE, iEs, XYP, IPEs, qEs)
C ================================================================
C Dynamic quality modification for tangled elements inside
C a super-element.
C
C Remark: non-efficient, time-consuming, but robust algorithm.
C ================================================================
      Real*8  XYP(3, *), qEs(*)
      Integer iEs(*), IPEs(5, *)

C group (Local variables)
      Integer ip(5)
      Real*8  v1, v2

C ================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1

      Do 100 i1 = 1, 4
         i2 = ip(i1 + 1)
         i3 = ip(i2 + 1)

         iP1 = IPEs(i1, nEs)
         iP2 = IPEs(i2, nEs)
         iP3 = IPEs(i3, nEs)

         Do 20 k = 1, lE
            If(iEs(k).LT.0)  Goto 20
            If(k.EQ.nEs)     Goto 20

            Do j1 = 1, 4
               j2 = ip(j1 + 1)
               j3 = ip(j2 + 1)

               jP1 = IPEs(j1, k)
               jP2 = IPEs(j2, k)
               jP3 = IPEs(j3, k)

               If(check33(iP1, iP2, iP3, jP1, jP2, jP3)) Then
                  i4  = ip(i3 + 1)
                  iP4 = IPEs(i4, nEs)

                  j4  = ip(j3 + 1)
                  jP4 = IPEs(j4, k)

                  v1 = calVol(XYP(1, iP1), XYP(1, iP2),
     &                        XYP(1, iP3), XYP(1, iP4))
                  v2 = calVol(XYP(1, iP1), XYP(1, iP2),
     &                        XYP(1, iP3), XYP(1, jP4))

                  If(v1 * v2.GE.0D0) Then
                     qEs(nEs) = -dabs(qEs(nEs))
                     qEs(k)   = -dabs(qEs(k))
                  End if

                  Goto 100
               End if
            End do
 20      Continue
 100  Continue

      Return
      End Subroutine updQb


C ================================================================
      Real*8 Function calSqr(xy1, xy2, xy3)
C ================================================================
      Real*8 xy1(3), xy2(3), xy3(3)

C Local variables
      Real*8 ax, ay, az, bx, by, bz

      ax = xy1(1) - xy3(1)
      ay = xy1(2) - xy3(2)
      az = xy1(3) - xy3(3)

      bx = xy2(1) - xy3(1)
      by = xy2(2) - xy3(2)
      bz = xy2(3) - xy3(3)

      calSqr = 5D-1 * dsqrt((ay * bz - az * by) ** 2 +
     &                      (ax * bz - az * bx) ** 2 +
     &                      (ax * by - ay * bx) ** 2)
      Return
      End Function calSqr


C ================================================================
      Subroutine calDet(HesP, detG)
C ================================================================
      include 'magic.fd'
C ================================================================
C Routine computes the determinant of the metric HesP.
C
C Remark: robustness of the overall code was increased by
C         replacing errMes() with the computation of |H|.
C ================================================================
      Real*8 HesP(6), detG

c Local arrays for LAPACK
      Real*8  A(3, 3), E(3), rW(15)
      Integer info

C ================================================================
      detG = HesP(1) * (HesP(2) * HesP(3) - HesP(5) ** 2) -
     &       HesP(4) * (HesP(4) * HesP(3) - HesP(5) * HesP(6)) +
     &       HesP(6) * (HesP(4) * HesP(5) - HesP(6) * HesP(2))

      If(detG.LE.0D0) Then
        A(1, 1) = HesP(1)
        A(2, 2) = HesP(2)
        A(3, 3) = HesP(3)

        A(1, 2) = HesP(4)
        A(2, 3) = HesP(5)
        A(1, 3) = HesP(6)

        Call dsyev('V', 'U', 3, A, 3, E, rW, 15, info)
        If(info.NE.0) Call errMes(3011, 'calDet',
     &                    'Error in LAPACK routine dsyev')

        E(1) = dabs(E(1))
        E(2) = dabs(E(2))
        E(3) = dabs(E(3))

        E(1) = max( E(1), E(3) * AniRatio )

        Do i = 1, 6
           HesP(i) = 0D0
        End do

        Do i = 1, 3
          HesP(1) = HesP(1) + E(i) * A(1, i) ** 2 
          HesP(2) = HesP(2) + E(i) * A(2, i) ** 2
          HesP(3) = HesP(3) + E(i) * A(3, i) ** 2

          HesP(4) = HesP(4) + E(i) * A(1, i) * A(2, i) 
          HesP(5) = HesP(5) + E(i) * A(2, i) * A(3, i) 
          HesP(6) = HesP(6) + E(i) * A(1, i) * A(3, i) 
        End do

        detG = HesP(1) * (HesP(2) * HesP(3) - HesP(5) ** 2) -
     &         HesP(4) * (HesP(4) * HesP(3) - HesP(5) * HesP(6)) +
     &         HesP(6) * (HesP(4) * HesP(5) - HesP(6) * HesP(2))
      End if

      Return
      End Subroutine calDet


C ================================================================
      Subroutine calQE(
C ================================================================
     &      Hes1, det1, xy1,
     &      Hes2, det2, xy2,
     &      Hes3, det3, xy3,
     &      Hes4, det4, xy4,
     &      hStar, qE, volume)
C ================================================================
C Computing quality of tetrahedron {xy1, ..., xy4} assuming that
C the metric field is linear.
C ================================================================
      include 'cubature.fd'
C ================================================================
      Real*8 Hes1(6), det1, xy1(3)
      Real*8 Hes2(6), det2, xy2(3)
      Real*8 Hes3(6), det3, xy3(3)
      Real*8 Hes4(6), det4, xy4(3)
      Real*8 hStar, qE, volume

C group (Local variables)
      Real*8 HesAvg(6), d1, d2, d3, d4, dsum
      Real*8 Pk, Vk
      Real*8 F, x1, y1, z1

C group (Function)
      F(x1) = (x1 * (2D0 - x1)) ** 5

C ================================================================
      Pk = 0D0
      Do n = 1, 6
         If(n.EQ.1) Then
            x1 = xy1(1) - xy4(1)
            y1 = xy1(2) - xy4(2)
            z1 = xy1(3) - xy4(3)
            Do i = 1, 6
               HesAvg(i) = (Hes1(i) + Hes4(i)) / 2
            End do
         Else If(n.EQ.2) Then
            x1 = xy2(1) - xy4(1)
            y1 = xy2(2) - xy4(2)
            z1 = xy2(3) - xy4(3)
            Do i = 1, 6
               HesAvg(i) = (Hes2(i) + Hes4(i)) / 2
            End do
         Else If(n.EQ.3) Then
            x1 = xy3(1) - xy4(1)
            y1 = xy3(2) - xy4(2)
            z1 = xy3(3) - xy4(3)
            Do i = 1, 6
               HesAvg(i) = (Hes3(i) + Hes4(i)) / 2
            End do
         Else If(n.EQ.4) Then
            x1 = xy1(1) - xy3(1)
            y1 = xy1(2) - xy3(2)
            z1 = xy1(3) - xy3(3)
            Do i = 1, 6
               HesAvg(i) = (Hes1(i) + Hes3(i)) / 2
            End do
         Else If(n.EQ.5) Then
            x1 = xy2(1) - xy3(1)
            y1 = xy2(2) - xy3(2)
            z1 = xy2(3) - xy3(3)
            Do i = 1, 6
               HesAvg(i) = (Hes2(i) + Hes3(i)) / 2
            End do
         Else If(n.EQ.6) Then
            x1 = xy1(1) - xy2(1)
            y1 = xy1(2) - xy2(2)
            z1 = xy1(3) - xy2(3)
            Do i = 1, 6
               HesAvg(i) = (Hes1(i) + Hes2(i)) / 2
            End do
         End if
         Pk = Pk + dsqrt(HesAvg(1) * x1 * x1 +
     &                   HesAvg(2) * y1 * y1 +
     &                   HesAvg(3) * z1 * z1 +
     &               2 * HesAvg(4) * x1 * y1 +
     &               2 * HesAvg(5) * y1 * z1 +
     &               2 * HesAvg(6) * x1 * z1)
      End do

      d1 = det1
      d2 = det2 
      d3 = det3 
      d4 = det4 
           
      dsum = dsqrt(d1 * T2B + (d2 + d3 + d4) * T2A) 
     &     + dsqrt(d2 * T2B + (d3 + d4 + d1) * T2A) 
     &     + dsqrt(d3 * T2B + (d4 + d1 + d2) * T2A) 
     &     + dsqrt(d4 * T2B + (d1 + d2 + d3) * T2A) 

      volume = calVol(xy1, xy2, xy3, xy4)
      Vk = volume * dsum / 4


      x1 = Pk / (6 * hStar)
      x1 = min(x1, 1D0 / x1)

      qE = 1832.8208D0 * dabs(Vk) / (Pk ** 3) * F(x1)

      Return
      End Subroutine calQE



C ================================================================
      Subroutine calQF(
C ================================================================
     &      Hes1, xy1, Hes2, xy2, Hes3, xy3, Hes4, xy4,
     &      hStar, iF, iR)
C ================================================================
      Real*8  Hes1(6), xy1(3)
      Real*8  Hes2(6), xy2(3)
      Real*8  Hes3(6), xy3(3)
      Real*8  Hes4(6), xy4(3)
      Real*8  hStar

      Integer iF(4), iR(6)

C group (Local variables)
      Real*8  HesAvg(6)
      Real*8  qF(4), qR(6)
      Real*8  F, x1, y1, z1

C group (Function)
      F(x1) = (x1 * (2D0 - x1)) ** 5

C ================================================================
      Do n = 1, 6
         If(n.EQ.1) Then
            x1 = xy1(1) - xy2(1)
            y1 = xy1(2) - xy2(2)
            z1 = xy1(3) - xy2(3)
            Do i = 1, 6
               HesAvg(i) = (Hes1(i) + Hes2(i)) / 2
            End do
         Else If(n.EQ.2) Then
            x1 = xy1(1) - xy3(1)
            y1 = xy1(2) - xy3(2)
            z1 = xy1(3) - xy3(3)
            Do i = 1, 6
               HesAvg(i) = (Hes1(i) + Hes3(i)) / 2
            End do
         Else If(n.EQ.3) Then
            x1 = xy1(1) - xy4(1)
            y1 = xy1(2) - xy4(2)
            z1 = xy1(3) - xy4(3)
            Do i = 1, 6
               HesAvg(i) = (Hes1(i) + Hes4(i)) / 2
            End do
         Else If(n.EQ.4) Then
            x1 = xy2(1) - xy3(1)
            y1 = xy2(2) - xy3(2)
            z1 = xy2(3) - xy3(3)
            Do i = 1, 6
               HesAvg(i) = (Hes2(i) + Hes3(i)) / 2
            End do
         Else If(n.EQ.5) Then
            x1 = xy2(1) - xy4(1)
            y1 = xy2(2) - xy4(2)
            z1 = xy2(3) - xy4(3)
            Do i = 1, 6
               HesAvg(i) = (Hes2(i) + Hes4(i)) / 2
            End do
         Else If(n.EQ.6) Then
            x1 = xy3(1) - xy4(1)
            y1 = xy3(2) - xy4(2)
            z1 = xy3(3) - xy4(3)
            Do i = 1, 6
               HesAvg(i) = (Hes3(i) + Hes4(i)) / 2
            End do
         End if

         x1 = dsqrt(HesAvg(1) * x1 * x1 +
     &              HesAvg(2) * y1 * y1 +
     &              HesAvg(3) * z1 * z1 +
     &          2 * HesAvg(4) * x1 * y1 +
     &          2 * HesAvg(5) * y1 * z1 +
     &          2 * HesAvg(6) * x1 * z1) / hStar
         x1 = min(x1, 1D0 / x1)

         iR(n) = n
         qR(n) = F(x1)
      End do


      qF(1) = min(qR(1), qR(2))
      qF(1) = min(qF(1), qR(4))

      qF(2) = min(qR(4), qR(5))
      qF(2) = min(qF(2), qR(6))

      qF(3) = min(qR(2), qR(3))
      qF(3) = min(qF(3), qR(6))

      qF(4) = min(qR(1), qR(3))
      qF(4) = min(qF(4), qR(5))

      Do i = 1, 4
         iF(i) = i
      End do


      Do i = 1, 3
         iMin = i
         Do j = i + 1, 4
            If(qF(j).LT.qF(iMin)) iMin = j
         End do

         x1 = qF(i)
         qF(i) = qF(iMin)
         qF(iMin) = x1

         k = iF(i)
         iF(i) = iF(iMin)
         iF(iMin) = k
      End do


      Do i = 1, 5
         iMin = i
         Do j = i + 1, 6
            If(qR(j).LT.qR(iMin)) iMin = j
         End do

         x1 = qR(i)
         qR(i) = qR(iMin)
         qR(iMin) = x1

         k = iR(i)
         iR(i) = iR(iMin)
         iR(iMin) = k
      End do

      Return
      End Subroutine calQF



C ================================================================
      Subroutine HesBnd(lP, iPs, ICP, HesP, HesPs)
C ================================================================
      include 'color.fd'
C ================================================================
      Integer iPs(*), ICP(*)
      Real*8  HesP(6, *), HesPs(6)

C (Local variables)
      Real*8  hesB(6), hesI(6)

C ================================================================
      nI = 0
      nB = 0

      Do i = 1, 6
         hesI(i) = 0D0
         hesB(i) = 0D0
      End do

      Do n = 1, lP
         iP = iPs(n)
         If(ifXnode(ICP(iP), jInode)) Then
            nI = nI + 1
            Do i = 1, 6
               hesI(i) = hesI(i) + HesP(i, iP)
            End do
         Else
            nB = nB + 1
            Do i = 1, 6
               hesB(i) = hesB(i) + HesP(i, iP)
            End do
         End if
      End do

      If(nI.GT.0) Then
         Do i = 1, 6
            HesPs(i) = hesI(i) / nI
         End do
      Else If(nB.GT.0) Then
         Do i = 1, 6
            HesPs(i) = hesB(i) / nB
         End do
      Else
         Call errMes(6001, 'HesBnd', 'system error')
      End if

      Return
      End Subroutine HesBnd


C ==========================================================
      Subroutine iniQ_analytic(nP, XYP, MetricFunction, HesP)
C ==========================================================
C  Three Fortran routines below create a metric field which
C  is 3x3 variable positive definite diagonal tensor HesP,
C
C             M11   M12   M13 
C      HesP = M12   M22   M23 
C             M13   M23   M33
C
C  where Mij = Mij(x, y, z).
C
C  The tensor element are enumerated in the following order:
C  HesP_{11}, HesP_{22}, HesP_{33}, HesP_{12}, HesP_{23}, HesP_{13}
C
C ==========================================================
      Real*8   XYP(3, *)
      Real*8   HesP(6, *)

      Integer  MetricFunction
      EXTERNAL MetricFunction

      Real*8   x, y, z, Metric(3, 3)
C ==========================================================
      Do n = 1, nP
         x = XYP(1, n)
         y = XYP(2, n)
         z = XYP(3, n)

         i = MetricFunction(x, y, z, Metric)

         HesP(1, n) = Metric(1,1)
         HesP(2, n) = Metric(2,2)
         HesP(3, n) = Metric(3,3)
         HesP(4, n) = Metric(1,2)
         HesP(5, n) = Metric(2,3)
         HesP(6, n) = Metric(1,3)
      End do

      Return
      End Subroutine iniQ_analytic


C ================================================================
      Real*8 Function avgQ(nE, qE, L1E, L2E)
C ================================================================
      Real*8  qE(*)
      Integer L2E(*), L1E(2, *) 

      avgQ = 0D0

      iE = L2E(1)
      Do n = 1, nE
         avgQ = avgQ + qE(iE)
         iE = L1E(2, iE)
      End do

      avgQ = avgQ / nE

      Return
      End Function avgQ
C
      End Module mba3d_makQ
