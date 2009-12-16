      Module mba3d_nlnfnc
C
      use mba3d_auxSP
      use mba3d_lintrp3D
      use mba3d_makM
      use mba3d_makQ
C
      contains
C
C ================================================================
      Double Precision Function NLnFnc(
C ================================================================
C group (F)
     &       nU, U,
C group (ANI)
     &       XYP, HesP, detG, hStar,
     &       iPs, lE, iEs, XYPs, ICPs, IPEs, HesPs, detGs, qEs,
     &       nPw, nEw, XYPw, HesPw, IPEw,
     &       MetricFunction, flagAnalytic,
     &       miLINTRP, mrLINTRP, iSE, rSE, iControl)
C ================================================================
      include 'color.fd'
C ================================================================
C group (F)
      Real*8 U(*)

C group (ANI)
      Integer iEs(*), IPEs(5, *), IPEw(4, *), iSE(*)
      Real*8  XYP(3, *), HesP(6, *), detG(*), hStar
      Real*8  HesPs(*),  detGs,  qEs(*), XYPs(*)
      Real*8  XYPw(3, *), HesPw(6, *), rSE(*)
      Logical flagAnalytic

      Integer  MetricFunction
      EXTERNAL MetricFunction

C group (Local variables)
      Integer ip(5)
      Real*8  XYPo(3), v

C ================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1


      Do i = 1, 3
         XYPs(i) = U(i)
      End do

      If(ifXnode(ICPs, jInode)) Then
         Do i = 1, 6
            HesPs(i) = HesP(i, iPs)
         End do

         LDH = 6
         nXY = 1
         If( .NOT.flagAnalytic ) Then
            Call LINTRP3D(nEw, IPEw, nPw, XYPw, LDH, HesPw, nXY, XYPs,
     &                    HesPs, iSE, miLINTRP, rSE, mrLINTRP, iControl)
         Else
            Call scaleBack(XYPs, XYPo)
            Call iniQ_analytic(1, XYPo, MetricFunction, HesPs)
         End if

         Call calDet(HesPs, detGs)
      Else
         Do i = 1, 6
            HesPs(i) = HesP(i, iPs)
         End do

         detGs = detG(iPs)
      End if


      NLnFnc = 1D0
      Do 10 n = 1, lE
         iE = iEs(n)
         If(iE.LE.0) Goto 10

         Do i1 = 1, 4
            If(IPEs(i1, n).EQ.iPs) Then
               i2 = ip(i1 + 1)
               i3 = ip(i2 + 1)
               i4 = ip(i3 + 1)

               iPa = IPEs(i2, n)
               iPb = IPEs(i3, n)
               iPc = IPEs(i4, n)

               Call calQE(
     &              HesPs, detGs, XYPs,
     &              HesP(1, iPa), detG(iPa), XYP(1, iPa),
     &              HesP(1, iPb), detG(iPb), XYP(1, iPb),
     &              HesP(1, iPc), detG(iPc), XYP(1, iPc),
     &              hStar, qEs(n), v)
               Goto 5
            End if
         End do
 5       NLnFnc = min(NLnFnc, qEs(n))
 10   Continue

      NLnFnc = 1D0 - NLnFnc
      Return
      End Function NlnFnc
C
      End Module mba3d_nlnfnc
