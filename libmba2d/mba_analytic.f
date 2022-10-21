C ==========================================================
      Subroutine mbaAnalytic(
C ==========================================================
c group (M)
     &      nP, MaxP, nF, MaxF, nE, MaxE, nPv,
     &      XYP, IPF, IPE, IPV,
     &      CrvFunction, ParCrv, iFnc,
     &      nEStar, 
c group (D)
     &      nFv, nEv, IFV, IEV, lbE,
     &      flagAuto, status,
c group (Q)
     &      MaxSkipE, MaxQItr,
     &      MetricFunction, Quality, rQuality, 
c group (W)
     &      MaxWr, MaxWi, rW, iW,
     &      iPrint, iERR)
C ==========================================================
      include 'lintrp.fd'
      include 'status.fd'
C ==========================================================
C  VARIABLES & PARAMETER are described in mba_nodal.f
C
C  MetricFunction - integer function created by the user (see 
C                   example in file forlibmba.f)
C
C    Integer Function MetricFunction(x, y, z, Metric)
C
C  This routine creates a metric at the given point (x,y, z). The
C  metric is a 2x2 positive definite symmetric tensor:
C
C                M11   M12
C      Metric =     
C                M12   M22
C
C  Only the upper triangular part of array Metric must be defined.
C
C
C *** Authors: K. Lipnikov (lipnikov@hotmail.com)
C              Y. Vassilevski (vasilevs@dodo.inm.ras.ru)
C ==========================================================
C group (M)
C     Integer MaxP, MaxF, MaxE, MaxFnc
      real  XYP(2, *)
      Integer IPE(3, *), IPF(4, *), IPV(*)

      EXTERNAL CrvFunction
      real   ParCrv(2, *)
      Integer  iFnc(*)

C group (D)
      Integer nFv, nEv
      Integer IFV(*), IEV(*), lbE(*)
      
      Logical flagAuto
      Integer status

C group (Q)
      Integer  MaxSkipE, MaxQItr
      real   Quality, rQuality
      Logical  flagAnalytic

      Integer  MetricFunction
      EXTERNAL MetricFunction

C group (W)
      real  rW(*)
      Integer iW(*)

C group (Local variables)
      real  hStar

C ==========================================================
C group (Common blocks)
      Integer iDomBnd, iMatBnd
      Common /aniBND/ iDomBnd, iMatBnd
 
      real  refXYP(2), scaXYP(2)
      Common /rescale/refXYP, scaXYP

C ==========================================================
      iERR = 0

      Do i = 1, 2
         refXYP(i) = 0D0
         scaXYP(i) = 1D0
      End do


c ... refine initial mesh when nE is very small
c ... it increases robustness of the code 
      Do while(nE < nEStar / 15 .AND. nE.LE.500 .AND. nEv+nFv.EQ.0)
         iIFE = 1
         iiW  = iIFE + 3 * nE
         nWi  = iiW  + 3 * nE + nP 
         If(nWi.GT.MaxWi) goto 100

         iSol = 1
         nWr  = iSol + MaxP
         If(nWr.GT.MaxWr) goto 100
 
         If(iPrint.GE.1) Write(*,5001) nP, nE

         Do i = 1, nP
            rW(iSol + i - 1) = 0D0
         End do 

         Call uniformRefinement(
     &        nP, MaxP, nF, MaxF, nE, MaxE,
     &        XYP, IPF, IPE, lbE,
     &        CrvFunction, ParCrv, iFnc, iW(iIFE),
     &        rW(iSol), 1, iW(iiW), MaxWi)
      End do 


 100  miLINTRP = 10 * nP + 3 * nE + 6
      mrLINTRP =  4 * nP + MaxH + 4

      inEt = 1
      inStept = inEt + MaxF
      inL2t = inStept + 4 * MaxF
      iLFnc = inL2t + MaxF
      iILt  = iLFnc + MaxF
      iL1Et = iILt + MaxF
      iL2Et = iL1Et + 2 * MaxF
      iIHolP = iL2Et + 2 * MaxF
      iIHolF = iIHolP + MaxP
      iIHolE = iIHolF + MaxF
      iICP = iIHolE + MaxE
      iIEP = iICP + MaxP
      iIFE = iIEP + MaxP
      iIEE = iIFE + 3 * MaxE
      iL1E = iIEE + 3 * MaxE
      iL2E = iL1E + 2 * MaxE
      iIPEw = iL2E + 2 * MaxE
      iiSE  = iIPEw + 3 * nE
      iIEPw = iiSE + miLINTRP
c ... we need twice less memory for backReferences
      inEPw = iIEPw + max(6 * nE, 4 * MaxF)
      nWi   = inEPw + max(3 * MaxP, 2 * MaxF)


      iHesP = 1
      itE = iHesP + 3 * MaxP
      idG = itE + MaxF
      iqE = idG + MaxP
      iHesPw = iqE + MaxE
      iXYPw = iHesPw + 3 * nP
      irSE = iXYPw + 2 * nP
      nWr  = irSE + max(mrLINTRP, max(nE, MaxF))


      iW(1) = nWi
      iW(2) = nWr
      If(nWi.GT.MaxWi) Then
         iERR = 1001
         goto 1000
      End if

      If(nWr.GT.MaxWr) Then
         iERR = 1002
         goto 1000
      End if


      Do n = 1, nWr
         rW(n) = 0D0
      End do

      Do n = 1, nWi
         iW(n) = 0
      End do


c ... compute the analytic metric
      Call iniQ_analytic(nP, XYP, MetricFunction, rW(iHesP))


c ... scale geometry to unit cube
      Call scale2Square(nP, XYP, .TRUE.)


c ... print Ani2D header
      If(iPrint.GE.1) Write(*, 5004) Quality, nEStar, MaxQItr


c ... set up default status
      Call setStatus(flagAuto, status, iPrint)


c ... call the main module
      flagAnalytic = .TRUE.
      Call ani2(
c group (M)
     &      nP, MaxP, nF, MaxF, nE, MaxE, nPv,
     &      XYP, IPF, IPE, IPV,
     &      CrvFunction, ParCrv, iFnc,
     &      nEStar, hStar,
     &      iW(iICP), iW(iIEP),
     &      iW(iIFE), iW(iIEE),
     &      iW(iL1E), iW(iL2E),
     &      iW(iIHolP), iW(iIHolF), iW(iIHolE),
     &      iW(iIEPw), iW(inEPw),
     &      iW(iIPEw), iW(iiSE),
c group (Dev)
     &      nFv, nEv, IFV, IEV, lbE,
     &      flagAuto, status,
c group (CRV)
     &      iW(iL1Et), iW(iL2Et), rW(itE),
     &      iW(inL2t), iW(inStept), iW(inEt),
     &      iW(iLFnc), iW(iILt),
c group (Q)
     &      MaxSkipE, MaxQItr, nQItr,
     &      rW(iHesP), Quality, rQuality,
     &      rW(idG), rW(iqE), rW(iXYPw), rW(iHesPw), rW(irSE),
     &      MetricFunction, flagAnalytic,
c group (ERR)
     &      iPrint, iERR)


c ... rescale geometry back
      Call scale2Square(nP, XYP, .FALSE.)


c ... returning sadditional information
      Do n = 1, nP
         iW(n) = iW(iICP + n - 1)
      End do
      iW(nP + 1) = nQItr

      rW(1) = 0D0
      rW(2) = rQuality
      rW(3) = hStar


 1000 If(iERR.EQ.0 .OR. iERR.EQ.1000) Return
      Call errMes(iERR, 'mbaMetric', 
     &            'See error.f for error description')

      Return

 5001 Format('Refining mesh:', I6, ' pts  and', I7, ' elements')

 5004 Format(/,
     &    'STONE FLOWER! (1997-2007), version 2.0', /,
     &    'Target: Quality', F5.2, ' with', I8, 
     &    ' triangles for at most', I8, ' iterations') 
      End



C ==========================================================
      Subroutine iniQ_analytic(nP, XYP, MetricFunction, HesP)
C ==========================================================
C  Three Fortran routines below create a metric field which
C  is 2x2 variable positive definite symmetric tensor HesP,
C             F(x,y)  H(x,y)
C      HesP =
C             H(x,y)  G(x,y)
C ==========================================================
C group (M)
      real   XYP(2, *)

C group (Q)
      Integer  MetricFunction
      EXTERNAL MetricFunction

      real   HesP(3, *)

C group (Local variables)
      real   x, y, Metric(2, 2)

C ==========================================================
      real  refXYP(2), scaXYP(2)
      Common /rescale/refXYP, scaXYP

C ==========================================================
      Do n = 1, nP
         x = refXYP(1) + XYP(1, n) * scaXYP(1)
         y = refXYP(2) + XYP(2, n) * scaXYP(2)

         i = MetricFunction(x, y, Metric)

         HesP(1, n) = Metric(1, 1)
         HesP(2, n) = Metric(2, 2)
         HesP(3, n) = Metric(1, 2)
      End do

      Return
      End




