C ================================================================
C  Changes since makQ_orig.f:
C  1. constant metric is replaced by linear metric
C  2. use LAPACK's routine dsyev instead of Linpack's dsvdc
C  3. added routine computing average quality
C ================================================================

C ================================================================
      Subroutine makQ(
C ================================================================
C Routine computes hStar and Quality of elements in the initial 
c mesh. 
C ================================================================
c group (M)
     &     nP, nE, XYP, IPE, IEE, nEV, IEV,
     &     nEStar, hStar, status,
c group (Q)
     &     HesP, detG, qE)
C ================================================================
      include 'status.fd'
C ================================================================
C group (M)
C     Integer MaxP, MaxF, MaxE, nEStar
      real  XYP(2, *)
      Integer IPE(3, *), IEE(3, *), IEV(*)
      Integer status

      real  hStar

C group (Q)
      real  HesP(3, *)
      real  detG(*), qE(*)

C group (Local variables)
      Integer iref(4)
      real  HesAvg(3), detAvg, dsum, calVol
      real  VStar, Vn, qEt
      Logical ifXnode

C ================================================================
      iref(1) = 1
      iref(2) = 2
      iref(3) = 3
      iref(4) = 1

      Do n = 1, nP
         Call calDet(HesP(1, n), detG(n))
      End do

      VStar = 0D0
      Do n = 1, nE
         dsum = 0D0

         Do i1 = 1, 3 
            i2 = iref(i1 + 1)
            iP1 = IPE(i1, n)
            iP2 = IPE(i2, n)

            Do i = 1, 3
               HesAvg(i) = (HesP(i, iP1) + HesP(i, iP2)) / 2
            End do

            Call calDet(HesAvg, detAvg)

            dsum = dsum + sqrt(detAvg)
         End do

         i3 = iref(i2 + 1)
         iP3 = IPE(i3, n)

         Vn = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))
         VStar = VStar + abs(Vn) * dsum / 3
      End do

      hStar = sqrt(VStar / nEStar * 4D0 / sqrt(3D0))

      Do n = 1, nE
         i1 = IPE(1, n)
         i2 = IPE(2, n)
         i3 = IPE(3, n)

         Call calQE(
     &        HesP(1, i1), XYP(1, i1),
     &        HesP(1, i2), XYP(1, i2),
     &        HesP(1, i3), XYP(1, i3),
     &        hStar, qE(n))

         If(ifXnode(status, ANISmoothMesh)) Then
            Do i = 1, 3
               k = IEE(i, n)
               If(k.GT.0) Then
                  j1 = IPE(1, k)
                  j2 = IPE(2, k)
                  j3 = IPE(3, k)

                  Call calQE(
     &                 HesP(1, j1), XYP(1, i1),
     &                 HesP(1, j2), XYP(1, i2),
     &                 HesP(1, j3), XYP(1, i3),
     &                 hStar, qEt)

                  qE(n) = min(qE(n), qEt)
               End if
            End do
         End if
      End do


c ... set quality of fixed elements to 1
      Do n = 1, nEv
         qE(IEV(n)) = 1D0
      End do

      Return
      End



C ================================================================
      real Function calVol(xy1, xy2, xy3)
C ================================================================
C Routine computes triangle area
C ================================================================
      real xy1(2), xy2(2), xy3(2)

      calVol = ((xy1(1) - xy3(1)) * (xy2(2)-xy3(2)) -
     &          (xy1(2) - xy3(2)) * (xy2(1)-xy3(1))) / 2
      Return
      End



C ================================================================
      Subroutine calDet(HesP, detG)
C ================================================================
      include 'magic.fd'
C ================================================================
C Routine computes determinat det(H) of matrix H where
C       | HesP(1)  HesP(3) |
C   H = |                  |
C       | HesP(3)  HesP(2) |  
C ================================================================
      real  HesP(3), detG

C Local arrays for LAPACK
      real  A(2, 2), E(2), rW(10)
      Integer info

C ================================================================
      detG = HesP(1) * HesP(2) - HesP(3) ** 2

      If(detG.LE.0D0) Then
        A(1, 1) = HesP(1)
        A(2, 2) = HesP(2)
        A(1, 2) = HesP(3)

        Call dsyev('V', 'U', 2, A, 2, E, rW, 10, info)
        If(info.NE.0) Call errMes(3011, 'calDet', 
     &                    'Error in Lapack routine dsyev')

        E(1) = abs(E(1))
        E(2) = abs(E(2))

        E(1) = max( E(1), E(2) * AniRatio )

        If(E(2).EQ.0D0) Then
           E(1) = AniEigenvalue
           E(2) = AniEigenvalue
        End if

        HesP(1) = E(1) * A(1, 1) ** 2 + E(2) * A(1, 2) ** 2
        HesP(2) = E(1) * A(2, 1) ** 2 + E(2) * A(2, 2) ** 2
        HesP(3) = E(1) * A(1, 1) * A(2, 1) + E(2) * A(1, 2) * A(2, 2)

        detG = HesP(1) * HesP(2) - HesP(3) ** 2
      End if

      Return
      End



C ================================================================
      Subroutine calQE(
C ================================================================
     &      Hes1, xy1, Hes2, xy2, Hes3, xy3,
     &      hStar, qE)
C ================================================================
C Routines computes quality of triangle {xy1, xy2, xy3} assuming
C that the metric is linear.
C
C *** Remarks:
C        1. round-off errors require to use |Lk|
C        2. moved division by 2 in HesAvg in calculation of Lk 
C ================================================================
      real Hes1(3), xy1(2)
      real Hes2(3), xy2(2)
      real Hes3(3), xy3(2)
      real hStar, qE

C group (Local variables)
      real HesAvg(3), det12, det13, det23, dsum
      real Pk, Vk, Lk
      real F, x1, y1, x2, y2

C group (FUnction)
      F(x1) = (x1 * (2D0 - x1)) ** 3

C ================================================================
      Do i = 1, 3
         HesAvg(i) = Hes2(i) + Hes3(i)
      End do
      Call calDet(HesAvg, det23)

      x1 = xy3(1) - xy2(1)
      y1 = xy3(2) - xy2(2)
      Lk = abs(HesAvg(1) * x1 * x1 + HesAvg(2) * y1 * y1 +
     &                            2 * HesAvg(3) * x1 * y1) / 2
      Pk = sqrt(Lk)

      Do i = 1, 3
         HesAvg(i) = Hes1(i) + Hes2(i)
      End do
      Call calDet(HesAvg, det12)

      x1 = xy1(1) - xy2(1)
      y1 = xy1(2) - xy2(2)
      Lk = abs(HesAvg(1) * x1 * x1 + HesAvg(2) * y1 * y1 +
     &                            2 * HesAvg(3) * x1 * y1) / 2
      Pk = Pk + sqrt(Lk)

      Do i = 1, 3
         HesAvg(i) = Hes1(i) + Hes3(i)
      End do
      Call calDet(HesAvg, det13)

      x2 = xy1(1) - xy3(1)
      y2 = xy1(2) - xy3(2)
      Lk = abs(HesAvg(1) * x2 * x2 + HesAvg(2) * y2 * y2 +
     &                            2 * HesAvg(3) * x2 * y2) / 2
      Pk = Pk + sqrt(Lk)

      dsum = (sqrt(det12) + sqrt(det23) + sqrt(det13)) / 2

c     Vk = abs(x1 * y2 - y1 * x2) * 5D-1 * dsum / 3
c     qE = 20.784619D0 * Vk / (Pk ** 2) * F(x1)

      Vk = abs(x1 * y2 - y1 * x2) * dsum

      x1 = Pk / (3 * hStar)
      If(x1.GT.1D0) x1 = 1D0 / x1

      qE = 3.46410316666D0 * Vk / (Pk ** 2) * F(x1)

      Return
      End



C ================================================================
      Subroutine calQF(
C ================================================================
     &      Hes1, det1, xy1,
     &      Hes2, det2, xy2,
     &      Hes3, det3, xy3,
     &      hStar, iw)
C ================================================================
C Routine orders mesh edges in assending of their quality 
C ================================================================
      real  Hes1(3), det1, xy1(2)
      real  Hes2(3), det2, xy2(2)
      real  Hes3(3), det3, xy3(2)
      real  hStar

      Integer iw(3)

C group (Local variables)
      real  HesMax(3), detMax
      real  qF(3)
      real  F, x1, y1

C group (Function)
      F(x1) = (x1 * (2D0 - x1)) ** 3

C ================================================================
      detMax = det1
      Do i = 1, 3
         HesMax(i) = Hes1(i)
      End do

      If(det2.GT.detMax) Then
         detMax = det2
         Do i = 1, 3
            HesMax(i) = Hes2(i)
         End do
      End if

      If(det3.GT.detMax) Then
         detMax = det3
         Do i = 1, 3
            HesMax(i) = Hes3(i)
         End do
      End if

      x1 = xy1(1) - xy2(1)
      y1 = xy1(2) - xy2(2)
      x1 = sqrt(HesMax(1) * x1 * x1 + HesMax(2) * y1 * y1 +
     &                             2 * HesMax(3) * x1 * y1) / hStar
      x1 = min(x1, 1D0 / x1)
      qF(1) = F(x1)


      x1 = xy3(1) - xy2(1)
      y1 = xy3(2) - xy2(2)
      x1 = sqrt(HesMax(1) * x1 * x1 + HesMax(2) * y1 * y1 +
     &                             2 * HesMax(3) * x1 * y1) / hStar
      x1 = min(x1, 1D0 / x1)
      qF(2) = F(x1)

      x1 = xy1(1) - xy3(1)
      y1 = xy1(2) - xy3(2)
      x1 = sqrt(HesMax(1) * x1 * x1 + HesMax(2) * y1 * y1 +
     &                             2 * HesMax(3) * x1 * y1) / hStar
      x1 = min(x1, 1D0 / x1)
      qF(3) = F(x1)

      iw(1) = 1
      iw(2) = 2
      iw(3) = 3

      Do i = 1, 2
         iMin = i
         Do j = i + 1, 3
            If(qF(j).LT.qF(iMin)) iMin = j
         End do

         x1 = qF(i)
         qF(i) = qF(iMin)
         qF(iMin) = x1

         k = iw(i)
         iw(i) = iw(iMin)
         iw(iMin) = k
      End do

      Return
      End



C ================================================================
      Subroutine updQE(XYP, lE, iEs, IPEs, 
     &                 HesP, rQuality, detG, hStar, qEs, flag)
C ================================================================
C Routine is used for mesh smoothing (status & ANISmoothMesh = 1)
C
C *** Remark: the routine is obsolete
C ================================================================
      Integer IPEs(3, *), iEs(*)
      real  XYP(2, *),  HesP(3, *)
      real  rQuality,  detG(*), hStar, qEs(*)
      Logical flag

C group (Local variables)
      Integer iref(4)
      real  qEt
      Logical flagFirst

C ==========================================================
      iref(1) = 1
      iref(2) = 2
      iref(3) = 3
      iref(4) = 1

      flagFirst = flag
      flag = .FALSE.

      Do 50 n = 1, lE
         If(iEs(n).LE.0) goto 50

         Do 30 i1 = 1, 3
            i2 = iref(i1 + 1)

            iP1 = IPEs(i1, n)
            iP2 = IPEs(i2, n)

            Do 20 k = n + 1, lE
               If(iEs(k).LE.0)  goto 20

               Do j1 = 1, 3
                  j2 = iref(j1 + 1)

                  jP1 = IPEs(j1, k)
                  jP2 = IPEs(j2, k)

                  If(iP1.EQ.jP1 .AND. iP2.EQ.jP2 .OR.
     &               iP1.EQ.jP2 .AND. iP2.EQ.jP1) Then

                     i3 = iref(i2 + 1)
                     iP3 = IPEs(i3, n)

                     j3 = iref(j2 + 1)
                     jP3 = IPEs(j3, k)

                     Call calQE(
     &                    HesP(1, jP1), XYP(1, iP1),
     &                    HesP(1, jP2), XYP(1, iP2),
     &                    HesP(1, jP3), XYP(1, iP3),
     &                    hStar, qEt)

                     If(qEt.LT.rQuality .AND. flagFirst) goto 1000
                     qEs(n) = min(qEs(n), qEt)

                     Call calQE(
     &                    HesP(1, iP1), XYP(1, jP1),
     &                    HesP(1, iP2), XYP(1, jP2),
     &                    HesP(1, iP3), XYP(1, jP3),
     &                    hStar, qEt)

                     If(qEt.LT.rQuality. AND. flagFirst) goto 1000
                     qEs(k) = min(qEs(k), qEt)

                     goto 30
                  End if
               End do
 20         Continue
 30      Continue
 50   Continue

      flag = .TRUE.

 1000 Return
      End



C ================================================================
      real Function avgQ(nE, qE, L1E, L2E)
C ================================================================
      real  qE(*)
      Integer L2E(*), L1E(2, *) 

      avgQ = 0D0

      iE = L2E(1)
      Do n = 1, nE
         avgQ = avgQ + qE(iE)
         iE = L1E(2, iE)
      End do

      avgQ = avgQ / nE

      Return
      End

