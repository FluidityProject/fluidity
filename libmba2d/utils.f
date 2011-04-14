C ================================================================
      Subroutine uniformRefinement(
C ================================================================
     &      nP, MaxP, nF, MaxF, nE, MaxE,  
     &      XYP, IPF, IPE, lbE,
     &      CrvFunction, ParCrv, iFnc, IRE,
     &      F, LDF, iW, MaxWi)
C ================================================================
C Routine refines the input mesh and interpolates (linearly) the
C nodal function F(LDF, *).
C 
C *** Remarks:
C        1. The size of working memory is 3 * nE + nP
C ================================================================
      real  XYP(2, *)
      Integer IPE(3, *), IPF(4, *), lbE(*)

      EXTERNAL CrvFunction
      real  ParCrv(2, *)
      Integer iFnc(*)

      Integer IRE(3, *), iW(*)
      Integer LDF
      real  F(LDF, *)

C (Local variables)
      Integer iref(4), iEt(8)
      real  t1, t2, t3, s
      Logical cmpE, tangled, check22, flagE

      DATA    iref /1,2,3,1/

      real  refXYP(2), scaXYP(2)
      Common /rescale/refXYP, scaXYP

C ==========================================================
      If(scaXYP(1).LE.0D0) Then
         Do i = 1, 2
            refXYP(i) = 0D0
            scaXYP(i) = 1D0
         End do
      End if

      inEP = 1
      iIEP = inEP + nP
      iEnd = iIEP + 3 * nE

      If(iEnd.GT.MaxWi) Call errMes(1001, 'uniformRefinement',
     &                             'not enough working memory')


c ... compute map E -> F
      Call listE2R(nP, nFtot, nE, IPE, IRE, iW(inEP), iW(iIEP))

      nPo = nP
      nFo = nF
      nEo = nE

      nP = nPo + nFtot
      nF = 2 * nFo
      nE = 4 * nEo

      If(nP.GT.MaxP) Call errMes(1003, 'uniformRefinement',
     &                          'local parameter MaxP is small')
      If(nF.GT.MaxF) Call errMes(1004, 'uniformRefinement',
     &                          'local parameter MaxF is small')
      If(nE.GT.MaxE) Call errMes(1006, 'uniformRefinement',
     &                          'local parameter MaxE is small')


c ... split edges
      Do 10 n = 1, nFo
         iP1 = IPF(1, n)
         iP2 = IPF(2, n)

         flagE = cmpE(iP1, iP2, iW(iIEP), iW(inEP), 0, iE)

         Do i1 = 1, 3
            i2 = iref(i1 + 1)

            jP1 = IPE(i1, iE)
            jP2 = IPE(i2, iE)

            If(check22(iP1, iP2, jP1, jP2)) Then
               kP1 = nPo + IRE(i1, iE)

               IPF(1, nFo + n) = iP1
               IPF(2, nFo + n) = kP1

               IPF(3, nFo + n) = 0
               IPF(4, nFo + n) = IPF(4, n)

               IPF(1, n) = kP1
               goto 10
            End if
         End do
  10  Continue


c ... split elements
      kE = nEo 
      Do n = 1, nEo
         Do i1 = 1, 3
            i2 = iref(i1 + 1)
            i3 = iref(i2 + 1)

            iP1 = IPE(i1, n)
            iP2 = IPE(i2, n)

            jP1 = nPo + IRE(i1, n)
            jP3 = nPo + IRE(i3, n)

            kE = kE + 1
            IPE(1, kE) = iP1
            IPE(2, kE) = jP1
            IPE(3, kE) = jP3

            lbE(kE) = lbE(n)

            Do i = 1, 2
               XYP(i, jP1) = (XYP(i, iP1) + XYP(i, iP2)) / 2
            End do

            Do i = 1, LDF
               F(i, jP1) = (F(i, iP1) + F(i, iP2)) / 2
            End do
         End do

         Do i = 1, 3
            IPE(i, n) = nPo + IRE(i, n)
         End do
      End do


c ... split curved edges
      kC = 0
      Do n = 1, nFo
         If(iFnc(n).GT.0) kC = kC + 1
      End do

      Do n = 1, nFo
         nC = IPF(3, n) 

         If(nC.GT.0) Then
            s = 1.0D0

            iloop = 0
 20         iloop = iloop + 1
            If(iloop.GT.3) Call errMes(6004, 'uniformRefinement', 
     &                          'Mesh is tangled after refinement')

            t1 = ParCrv(1, nC)
            t2 = ParCrv(2, nC)

            s = s / 2
            t3 = s * t1 + (1 - s) * t2
            ParCrv(1, nC) = t3

            kP1 = IPF(1, n)
            Call aniCrv(t3, XYP(1, kP1), iFnc(nC), CrvFunction)

c  ...  check for inverted elements    
            iP1 = IPF(1, nFo + n) 
            iP2 = IPF(2, n) 

            nEt = 0
            iE1 = 0
            Do k = 1, 2
               If(cmpE(iP1, iP2, iW(iIEP), iW(inEP), iE1, iE2)) Then
                  Do i = 1, 3 
                     iEt(nEt + i) = nEo + 3 * (iE2 - 1) + i 
                  End do
                  iEt(nEt + 4) = iE2
                  nEt = nEt + 4
               End if

               iE1 = iE2
            End do

            Do i = 1, nEt
               Do j = max(i + 1, 4), nEt
                  If(tangled(iEt(i), iEt(j), XYP, IPE)) goto 20
               End do
            End do

c  ...  add new curved edge
            kC = kC + 1
            IPF(3, nFo + n) = kC

            iFnc(kC) = iFnc(nC)
            ParCrv(1, kC) = t1
            ParCrv(2, kC) = t3
         End if
      End do

      Return
      End



C ================================================================
      Subroutine orientBoundary(nP, nF, nE, XYP, IPF, IPE, iW, MaxWi)
C ================================================================
C Routine orients the external boundary of a given mesh in such 
C a way that the domain is on the left of an edge. In orwer words, 
C IPF(1, *) and IPF(2, *) are flipped if neccessary.
C
C *** Remarks:
C        1. The size of working memory is 3 * nE + 2 * nF + nP
C ================================================================
      real  XYP(2, *)
      Integer IPE(3, *), IPF(4, *), iW(*)

C (Local variables)
      real  v, calVol
      Integer iref(4)
      DATA    iref /1,2,3,1/

C ================================================================
      inEP = 1
      iIEP = inEP + nP
      iIFE = iIEP + 2 * nF
      iEnd = iIFE + 3 * nE

      If(iEnd.GE.MaxWi) Call errMes(1001, 'orientBoundary',
     &                             'not enough working memory')

c ... compute map E -> F
      Call listE2F(nP, nF, nE, IPF, IPE, iW(iIFE), iW(inEP), iW(iIEP))

      Do n = 1, nE
         Do i1 = 1, 3
            ife = iW(iIFE + 3 * (n - 1) + i1 - 1)
            If(ife.GT.0) Then 
               i2 = iref(i1 + 1)
               i3 = iref(i2 + 1)
 
               iP1 = IPF(1, ife)  
               iP2 = IPF(2, ife)  
               iP3 = IPE(i3, n)  

               v = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))
               If(v.GT.0D0) Call swapii(IPF(1, ife), IPF(2, ife))
            End if
         End do
      End do      

      Return
      End



C ================================================================
      Subroutine listE2R(nP, nR, nE, IPE, IRE, nEP, IEP)
C ================================================================
C  The routine creates connectivity lists E->R for mesh edges 
C
C  *** Remarks:
C         1. Working memory is nEP(nP), IEP(3 * nE)
C ================================================================
      Integer IPE(3, *), IRE(3, *)
      Integer IEP(*), nEP(*)

C ================================================================
C group (Local variables)
      Integer iref(4)
      Logical cmpE, check22

      DATA    iref /1,2,3,1/

C ================================================================
      Call backReferences(nP, nE, 3, 3, IPE, nEP, IEP)

      Do n = 1, nE
         Do i = 1, 3
            IRE(i, n) = 0
         End do
      End do

      nR = 0
      Do n = 1, nE
         Do 10 i1 = 1, 3
            If(IRE(i1, n).EQ.0) Then
               nR = nR + 1
               IRE(i1, n) = nR

               i2 = iref(i1 + 1)

               iP1 = IPE(i1, n)
               iP2 = IPE(i2, n)

               If(cmpE(iP1, iP2, IEP, nEP, n, iE2)) Then
                  Do j1 = 1, 3
                     j2 = iref(j1 + 1)

                     jP1 = IPE(j1, iE2)
                     jP2 = IPE(j2, iE2)
                     If(check22(iP1, iP2, jP1, jP2)) Then
                        IRE(j1, iE2) = nR
                        goto 10
                     End if
                  End do
               End if
            End if
 10      Continue
      End do

      Return
      End



C ================================================================
      Subroutine listR2R(nP, nR, nE, MaxL, IPE, nRR, IRR, iW)
C ================================================================
C  The routine creates connectivity lists R->R for mesh edges 
C
C  *** Remarks:
C         1. iW(*) - working memory of size 9 * nE
C ================================================================
      Integer IPE(3, *), nRR(*), IRR(*)
      Integer iW(*)

C ================================================================
      iIRE = 1
      inEP = iIRE + 3 * nE
      iIEP = inEP + nP
      iEnd = iIEP + 3 * nE 
      
      Call listE2R(nP, nR, nE, IPE, iW(iIRE), iW(inEP), iW(iIEP))

      inER = inEP 
      iIER = inER + nR
      iEnd = iIER + 3 * nE
      Call backReferences(nR, nE, 3,3, iW(iIRE), iW(inER), iW(iIER))

      nL = 0 

      i2 = 0
      Do n = 1, nR
         i1 = i2 + 1
         i2 = iW(inER + n - 1)

         Do m = i1, i2
            iE = iW(iIER + m - 1)

            Do 100 j = 1, 3
               iRt = iW(iIRE + 3 * (iE - 1) + j - 1)

               If(iRt.EQ.n .AND. m.GT.i1) goto 100

               nL = nL + 1
               If(nL.GT.MaxL) 
     &            Call errMes(2011, 'listR2R',
     &                       'user parameter MaxL is small')

               IRR(nL) = iRt
 100        Continue
         End do

         nRR(n) = nL
      End do

      Return
      End



C ================================================================
      Subroutine listR2P(nP, nR, nE, MaxR, IPE, IPR, nEP, IEP)
C ================================================================
C  The routine creates connectivity lists R->P for mesh edges 
C
C  *** Remarks:
C         1. Working memory is nEP(nP), IEP(3 * nE)
C ================================================================
      Integer IPE(3, *), IPR(2, *)
      Integer IEP(*), nEP(*)

C ================================================================
C group (Local variables)
      Integer iref(4)
      Logical cmpE, flag

      DATA    iref /1,2,3,1/

C ================================================================
      Call backReferences(nP, nE, 3, 3, IPE, nEP, IEP)

      nR = 0
      Do n = 1, nE

         Do 10 i1 = 1, 3
            i2 = iref(i1 + 1)

            iP1 = IPE(i1, n)
            iP2 = IPE(i2, n)

            flag = cmpE(iP1, iP2, IEP, nEP, n, iE2)

            If(iE2.EQ.0 .OR. iE2.GT.n) Then
               nR = nR + 1
               If(nR.GT.MaxR) Then
                  Call errMes(2011, 'listR2P',
     &                       'user parameter MaxR is small')
               End if

               IPR(1, nR) = iP1
               IPR(2, nR) = iP2

               goto 10
            End if
 10      Continue
      End do

      Return
      End



C ================================================================
      Subroutine listP2P(nP, nE, MaxList, IPE, nPP, IPP, iW)
C ================================================================
C  The routine creates connectivity lists P->P for mesh points.
C
C  *** Remarks:
C         1. iW(*) - working memory of size 3 * nE + nP
C ================================================================
      Integer IPE(3, *), nPP(*), IPP(*)
      Integer iW(*)

C ================================================================
      inEP = 0
      iIEP = inEP + nP

      Call backReferences(nP, nE, 3,3, IPE, iW(inEP + 1), iW(iIEP + 1))
 
c ... main algorithm: array nEP is overloaded inside
      nL = 0

      i2 = 0
      Do n = 1, nP
         nLo = nL

         i1 = i2 + 1
         i2 = iW(inEP + n)

         Do m = i1, i2
            iE = iW(iIEP + m)

            Do j = 1, 3
               iPt = IPE(j, iE)

               If(iW(inEP + iPt).GT.0) Then
                  nL = nL + 1
                  If(nL.GT.MaxList) Call errMes(2011, 'listP2P',
     &                             'user parameter MaxList is small')

                  IPP(nL) = iPt

                  iW(inEP + iPt) = -iW(inEP + iPt)
               End if
            End do
         End do

         nPP(n) = nL

c  ...  recovering values of array nEP
         Do m = nLo + 1, nL
            iPt = IPP(m)
            iW(inEP + iPt) = -iW(inEP + iPt)
         End do
      End do

      Return
      End



C ================================================================
      Subroutine listE2F(nP, nF, nE, IPF, IPE, IFE, nEP, IEP)
C ================================================================
C  The routine computes connectivity list E->F for BOUNDARY edges
C
C  *** Remarks:
C         1. Working memory is nEP(nP), IEP(2 * nF)
C ================================================================
      Integer IPF(4, *), IPE(3, *), IFE(3, *)
      Integer IEP(*), nEP(*)

C ================================================================
C group (Local variables)
      Integer iref(4)
      Logical cmpE

      DATA    iref /1,2,3,1/

C ================================================================
      Call backReferences(nP, nF, 2, 4, IPF, nEP, IEP)

      Do n = 1, nE
         Do i1 = 1, 3
            IFE(i1, n) = 0

            i2 = iref(i1 + 1)

            iP1 = IPE(i1, n)
            iP2 = IPE(i2, n)

            If(cmpE(iP1, iP2, IEP, nEP, 0, iF)) Then
               IFE(i1, n) = iF
            End if
         End do
      End do

      Return
      End



C ================================================================
      Subroutine listE2E(nP, nE, IPE, IEE, nEP, IEP)
C ================================================================
C  The routine computes connectivity lists E->E for neighboring
C  triangles. 
C
C  *** Remarks:
C         1. Working memory is nEP(nP), IEP(3 * nE)
C ================================================================
      Integer IPE(3, *), IEE(3, *), IEP(*), nEP(*)

C group (Local variables)
      Integer iref(4)
      Logical cmpE

      DATA    iref /1,2,3,1/

C ================================================================
      Call backReferences(nP, nE, 3, 3, IPE, nEP, IEP)

      Do n = 1, nE
         Do i1 = 1, 3
            i2 = iref(i1 + 1)

            ip1 = IPE(i1, n)
            ip2 = IPE(i2, n)

            IEE(i1, n) = 0
            If(cmpE(ip1, ip2, IEP, nEP, n, iE2)) Then
               IEE(i1, n) = iE2
            End if
         End do
      End do

      Return
      End



C ================================================================
      Subroutine listConv(
     &           nP, nR, nE, nEP, IEP, IRE, nX, MaxX, nRP, IRP, iW)
C ================================================================
C  The routine convolutes unstructured maps X->Y and Y->Z to get 
C  the map X->Z. For examples, if X means points (P), Y means
C  elements (E), and Z means edges (R), we get the map from a point
C  to all edges in the elements having this point.  
C
C  Only the first map X->Y is structured
C
C  *** Remarks:
C         1. Working memory is iW(nR)
C ================================================================
      Integer nP, nR, nE, MaxX
      Integer nEP(*), IEP(*), IRE(3, *), nRP(*), IRP(*), iW(*)
C ================================================================
      Do n = 1, nR
         iW(n) = 1
      End do

      nX = 0
      i2 = 0

      Do n = 1, nP
         nX0 = nX

         i1 = i2 + 1
         i2 = nEP(n)

         Do i = i1, i2
            iE = IEP(i)

            Do j = 1, 3
               iR = IRE(j, iE)  
               If(iW(iR).GT.0) Then
                  nX = nX + 1
                  If(nX.GT.MaxX) Call errMes(2011, 'listConv',
     &                               'user parameter MaxX is small')
                  IRP(nX) = iR

                  iW(iR) = -iW(iR)
               End if
            End do
         End do 

         nRP(n) = nX
        
c ...    restore array iW
         Do i = nX0 + 1, nX
            iR = IRP(i)
            iW(iR) = -iW(iR)
         End do 
      End do      

      Return
      End



C ================================================================
      Subroutine reverseMap(nP, nR, nRP, IRP, nPR, IPR)
C ================================================================
C Routine creates map R->P reverse to the map P->R.
C ================================================================
      Integer nRP(*), IRP(*), nPR(*), IPR(*)

      Do n = 1, nR
         nPR(n) = 0
      End do

      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = nRP(n)

         Do i = i1, i2
            iR = IRP(i)
            nPR(iR) = nPR(iR) + 1
         End do
      End do

      Do n = 2, nR
         nPR(n) = nPR(n) + nPR(n - 1)
      End do

      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = nRP(n)

         Do i = i1, i2
            iR = IRP(i)

            iref = nPR(iR)
            IPR(iref) = n

            nPR(iR) = iref - 1
         End do
      End do

      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = nRP(n)

         Do i = i1, i2
            iR = IRP(i)
            nPR(iR) = nPR(iR) + 1
         End do
      End do

      Return
      End



C ================================================================
      Subroutine backReferences(nP, nE, L, M, IPE, nEP, IEP)
C ================================================================
C Routine creates map P->E reverse to the map E->P.
C     nEP(P) - nEP(P-1) = number of elements having common 
C                         point P.
C     IPE([nEP(P-1) + 1 : nEP(P)]) = list of elements having
C                         common point P.  
C ================================================================
      Integer IPE(M, *), nEP(*), IEP(*)

      Do n = 1, nP
         nEP(n) = 0
      End do

      Do n = 1, nE
         Do i = 1, L
            i1 = IPE(i, n)
            nEP(i1) = nEP(i1) + 1
         End do
      End do

      Do n = 2, nP
         nEP(n) = nEP(n) + nEP(n - 1)
      End do

      Do n = 1, nE
         Do i = 1, L
            i1 = IPE(i, n)

            iref = nEP(i1)
            IEP(iref) = n

            nEP(i1) = iref - 1
         End do
      End do

      Do n = 1, nE
         Do i = 1, L
            i1 = IPE(i, n)
            nEP(i1) = nEP(i1) + 1
         End do
      End do

      Return
      End



