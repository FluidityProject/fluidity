C ==========================================================
      real Function NLnFnc(
C ==========================================================
C group (F)
     &       nU, U,
C group (ANI)
     &       XYP, IPE, IEE, HesP, hStar, status,
     &       lE, iEs, XYPs, IPEs, HesPs, detGs, qEs,
     &       nPw, nEw, XYPw, HesPw, IPEw,
     &       MetricFunction, flagAnalytic,
     &       iSE, rSE, iP1, iFNCs, calCrv,
     &       L1Et, L2Et, tE,
     &       nL2t, nStept, nEt, nCrvFnc, LFnc, ILt)
C ==========================================================
      include 'status.fd'
C ==========================================================
C group (F)
      real U(*)

C group (ANI)
      Integer IPE(3, *), IEE(3, *)
      Integer iEs(*), IPEs(3, *), IPEw(3, *), iSE(*), status
      real  XYP(2, *), HesP(3, *), hStar

      real  XYPs(2),   HesPs(*),  detGs,  qEs(*)
      real  XYPw(2, *), HesPw(3, *), rSE(*)
      real  tE(*)

      EXTERNAL calCrv

      Integer L1Et(2, *), L2Et(*)
      Integer nL2t(*), nStept(4, *), nEt(*)
      Integer LFnc(*), ILt(*)

      Logical flagAnalytic

      Integer  MetricFunction
      EXTERNAL MetricFunction

C group (Local variables)
      Integer ip(4)
      real  prjXYPs(2), XYPo(2), tc, qEt
      Logical check1j, ifXnode

C ==========================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      If(nU.EQ.2) Then
         Do i = 1, 2
            XYPs(i) = U(i)
         End do
      Else if(nU.EQ.1) Then
         tc = U(1)

         Call aniCrv(tc, XYPs, iFNCs, calCrv)
      End if

      Call findSE(nCrvFnc, LFnc, iFNCs, k)
      If(k.GT.0) Then
         ir = ILt(k)
         Call prjCrv(XYPs, prjXYPs, iFNCs, tc, calCrv,
     &               L1Et(1, ir), L2Et(ir), nL2t(k), nStept(1, k),
     &               nEt(k), tE(ir))

         If(.NOT.flagAnalytic) Then
         Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1, 
     &               prjXYPs, HesPs, iSE, rSE, .FALSE.)
      Else 
            Call scaleBack(prjXYPs, XYPo)
            Call iniQ_analytic(1, XYPo, MetricFunction, HesPs)
         End if
      Else 
         If(.NOT.flagAnalytic) Then
         Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1, 
     &               XYPs, HesPs, iSE, rSE, .FALSE.)
         Else
            Call scaleBack(XYPs, XYPo)
            Call iniQ_analytic(1, XYPo, MetricFunction, HesPs)
         End if
      End if

      Call calDet(HesPs, detGs)


      Do 50 n = 1, lE
         iE = iEs(n)
         If(iE.LE.0) goto 50

         Do i1 = 1, 3
            If(IPEs(i1, n).EQ.iP1) Then
               i2 = ip(i1 + 1)
               i3 = ip(i2 + 1)

               iPa = IPEs(i2, n)
               iPb = IPEs(i3, n)

               Call calQE(
     &              HesPs,        XYPs,
     &              HesP(1, iPa), XYP(1, iPa),
     &              HesP(1, iPb), XYP(1, iPb),
     &              hStar, qEs(n))

c  ...  updating the quality
               If(ifXnode(status, ANISmoothMesh)) Then
                  Do 10 j1 = 1, 3
                     jE = IEE(j1, iE)
                     If(jE.LE.0) goto 10

                     jP1 = IPE(1, jE)
                     jP2 = IPE(2, jE)
                     jP3 = IPE(3, jE)

                     If(check1j(iP1, IPE(1, jE))) Then 
                        if(jP2.EQ.iP1) jP2 = jP1 
                        if(jP3.EQ.iP1) jP3 = jP1 

                        Call calQE(
     &                       HesPs,        XYPs,
     &                       HesP(1, jP2), XYP(1, iPa),
     &                       HesP(1, jP3), XYP(1, iPb),
     &                       hStar, qEt)

                        qEs(n) = min(qEs(n), qEt)
                     Else
                        Call calQE(
     &                       HesP(1, jP1), XYPs,
     &                       HesP(1, jP2), XYP(1, iPa),
     &                       HesP(1, jP3), XYP(1, iPb),
     &                       hStar, qEt)

                        qEs(n) = min(qEs(n), qEt)

                        Call calQE(
     &                       HesPs,        XYP(1, jP1),
     &                       HesP(1, iPa), XYP(1, jP2),
     &                       HesP(1, iPb), XYP(1, jP3),
     &                       hStar, qEt)

                        Call findSE(lE, iEs, jE, m)
                        If(m.LE.0) Call errMes(6002, 
     &                                  'NLnFnc', 'System error') 
                        qEs(m) = min(qEs(m), qEt)
                     End if
 10               Continue
               End if
               
               goto 50
            End if
         End do
 50   Continue

      NLnFnc = 1D0
      Do n = 1, lE
         NLnFnc = min(NLnFnc, qEs(n))
      End do
      NLnFnc = 1D0 - NLnFnc
      Return
      End
