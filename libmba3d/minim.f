      Module mba3d_minim
C
      use mba3d_auxSE
      use mba3d_makQ
      use mba3d_nlnfnc
C
      contains
C
C ================================================================
      Subroutine minim(
C ================================================================
C group(F)
     &       nU, U, ZZ, Lamda1, Lamda2, fMin, U1,
     &       icnt, rMove, flagIni, flagResult,
C group (ANI)
     &       XYP, IPE, detG, HesP, hStar,
     &       iPs, lE, iEs, XYPs, ICPs, IPEs, detGs, HesPs, qEs,
     &       nPw, nEw, XYPw, HesPw, IPEw,
     &       MetricFunction, flagAnalytic,
     &       miLINTRP, mrLINTRP, iSE, rSE, iControl,
     &       flagTM, nBad)
C ================================================================
      include 'makS.fd'
      include 'status.fd'
C ================================================================
C MINIMIZATION OF THE NONLINEAR FUNCTIONAL FUNC(U)
C IN GIVEN DIRECTION  ZZ
C ================================================================
C group (F)
      Real*8  U(*), ZZ(*), U1(*)
      Real*8  Lamda1, Lamda2, fMin, rMove
      Logical flagIni, flagResult

C group (ANI)
      Integer IPE(5, *), iEs(*), IPEs(5, *), IPEw(4, *), iSE(*)
      Real*8  XYP(3, *),  HesP(6, *), detG(*), hStar
      Real*8  XYPs(3, *), HesPs(*),   detGs, qEs(*)
      Real*8  XYPw(3, *), HesPw(6, *), rSE(*)
      Logical flagAnalytic

      Integer  MetricFunction
      EXTERNAL MetricFunction

      Logical flagTM
      Integer nBad

C ================================================================
C group (Local variables)
      Real*8  Lamda, dLamda, step, x, f0, f1
      Real*8  aStep

      Real*8  XYPu(3), XYPt(3), qEt(MaxS), HesPt(6), detGt

      Integer iOs(MaxS)
      Logical flagOrient

C ================================================================
      flagResult = .FALSE.

      x = 0
      Do i = 1, nU
         x = x + ZZ(i) ** 2
      End do
      x = dsqrt(x)
      If(x.EQ.0D0) Goto 9000

      Do i = 1, nU
         ZZ(i) = ZZ(i) / x
      End do

      x = Lamda1
      Lamda1 = min(x, Lamda2)
      Lamda2 = max(x, Lamda2)

      icnt = 0
      step = 5D-3
      aStep = 0D0
      dLamda = Lamda2 - Lamda1

      If(flagIni) Then
         f0 = NLnFnc(nU, U,
C group (ANI)
     &        XYP, HesP, detG, hStar,
     &        iPs, lE, iEs, XYPs, ICPs, IPEs, HesPs, detGs, qEs,
     &        nPw, nEw, XYPw, HesPw, IPEw,
     &        MetricFunction, flagAnalytic,
     &        miLINTRP, mrLINTRP, iSE, rSE, iControl)
      Else
         f0 = fMin
      End if

      Call copySQ(lE, qEs, XYPs, HesPs, detGs,
     &                qEt, XYPt, HesPt, detGt)


 10   Lamda = step  * dLamda
      If(Lamda.GT.(1D0 - aStep - step) * dLamda) Then
         f1 = f0 + 1D0
      Else
         Do i = 1, nU
            U1(i) = U(i) + Lamda * ZZ(i)
         End do


         icnt = icnt + 1
         f1 = NLnFnc(nU, U1,
C group (ANI)
     &        XYP, HesP, detG, hStar,
     &        iPs, lE, iEs, XYPs, ICPs, IPEs, HesPs, detGs, qEs,
     &        nPw, nEw, XYPw, HesPw, IPEw,
     &        MetricFunction, flagAnalytic,
     &        miLINTRP, mrLINTRP, iSE, rSE, iControl)


C  ...  checking for inverted elements
         If(flagTM) Then
            Do i = 1, 3
               XYPu(i) = XYP(i, iPs)
               XYP(i, iPs) = XYPs(i, 1)
            End do

            Do n = 1, lE
               If(iEs(n).GE.0) Then
                  Call updQb(n, lE, iEs, XYP, IPEs, qEs)
               End if
            End do

            Do i = 1, 3
               XYP(i, iPs) = XYPu(i)
            End do
 
            mBad = 0
            Do n = 1, lE
               If(qEs(n).LE.0D0) mBad = mBad + 1
            End do
            If(mBad.GT.nBad) f1 = f0 + 1D0

C  ...  checking for the tetrahedrons orientation
         Else    
            Call calSO(XYP, IPE, lE, iEs, iOs)
            Call chkSO(iPs, XYPs, XYP, IPE, lE, iEs, iOs, flagOrient)
            If(.NOT.flagOrient) f1 = f0 + 1D0
         End if
      End if

      If(f1.LT.f0) Then
         f0 = f1
         Do i = 1, nU
            U(i) = U1(i)
         End do

         flagResult = .TRUE.

         Call copySQ(lE, qEs, XYPs, HesPs, detGs,
     &                   qEt, XYPt, HesPt, detGt)

         aStep = aStep + step
         step = 4D0 * step
      Else
         step = step / 4D0
         If(step.LT.1D-3) Then
            fMin = f0
            rMove = aStep * dLamda / Lamda2
            Goto 1000
         End if
      End if
      Goto 10

 1000 Continue
      Call copySQ(lE, qEt, XYPt, HesPt, detGt,
     &                qEs, XYPs, HesPs, detGs)

 9000 Return
      End Subroutine minim
C
      End module mba3d_minim
