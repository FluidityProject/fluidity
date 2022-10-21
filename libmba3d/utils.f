      Module mba3d_utils
C
      use mba3d_auxSE
      use mba3d_auxSF
      use mba3d_auxSP
      use mba3d_cmpE
      use mba3d_error
C
      contains
C
C ================================================================
      Subroutine listP2P(nP, nE, MaxList, IPE, nPP, IPP, iW)
C ================================================================
C  The routine creates connectivity lists P->P for mesh points.
C
C  *** Remarks:
C         1. iW(*) - working memory of size 4 * nE + nP
C ================================================================
      Integer IPE(4, *), nPP(*), IPP(*)
      Integer iW(*)

C ================================================================
      inEP = 0
      iIEP = inEP + nP

      Call backReferences(nP, nE, 4,4, IPE, iW(inEP + 1), iW(iIEP + 1))
 
c ... main algorithm: array nEP is overloaded inside
      nL = 0

      i2 = 0
      Do n = 1, nP
         nLo = nL

         i1 = i2 + 1
         i2 = iW(inEP + n)

         Do m = i1, i2
            iE = iW(iIEP + m)

            Do j = 1, 4
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
      End Subroutine listP2P


 
C ================================================================
      Subroutine listR2P(nP, nR, nE, MaxR, IPE, IPR, iW)
C ================================================================
C  The routine reates connectivity list R -> P for mesh edges. 
C  The algorithm has a linear arithmetical complexity. 
C
C  Parameters:
C      nP, nR, nE - the number of points, edges and elements
C      MaxR       - the maximal number of edges
C
C      IPE(4, nE) - connectivity list of tetrahedra: E -> P
C      IPR(2, nR) - connectivity list of edges:      R -> P
C
C      iW(*) - working memory of size 4 * nE + nP
C ================================================================
      Integer IPE(4, *), IPR(2, *)
      Integer iW(*)

C ================================================================
      inEP = 0
      iIEP = inEP + nP

      Call backReferences(nP, nE, 4, 4, IPE, iW(inEP + 1), iW(iIEP + 1))


c ... main algorithm: array nEP is overloaded inside
      nR = 0
   
      i2 = 0
      Do n = 1, nP
         nRo = nR

         i1 = i2 + 1
         i2 = iW(inEP + n)

         Do m = i1, i2
            iE = iW(iIEP + m)

            Do j = 1, 4
               iPt = IPE(j, iE)

               If(iPt.GT.n .AND. iW(inEP + iPt).GT.0) Then
                  nR = nR + 1
                  If(nR.GT.MaxR) Then
                     Call errMes(2011, 'listR2P', 
     &                          'user parameter MaxR is small')
                  End if

                  IPR(1, nR) = n
                  IPR(2, nR) = iPt

                  iW(inEP + iPt) = -iW(inEP + iPt)
               End if
            End do
         End do

c  ...  recovering values of array nEP
         Do m = nRo + 1, nR
            iPt = IPR(2, m)
            iW(inEP + iPt) = -iW(inEP + iPt)
         End do
      End do

      Return
      End Subroutine listR2P



C ================================================================
      Subroutine listR2R(nP, nR, nE, MaxL, IPE, nRR, IRR, iW, iERR)
C ================================================================
C  The routine creates connectivity lists R->R for mesh edges.
C  Routine returns 0 upon successful completion.
C
C  *** Remarks:
C         1. iW(*) - working memory of size 12 * nE + nR
C ================================================================
      Integer IPE(4, *), nRR(*), IRR(*)
      Integer iW(*)

C ================================================================
      iERR = 0

      iIRE = 1
      inEP = iIRE + 6 * nE
      iIEP = inEP + nP
      iEnd = iIEP + 4 * nE

      Call listE2R(nP, nR, nE, IPE, iW(iIRE), iW(inEP), iW(iIEP))

      inER = inEP
      iIER = inER + nR
      iEnd = iIER + 6 * nE
      Call backReferences(nR, nE, 6,6, iW(iIRE), iW(inER), iW(iIER))

      nL = 0

      i2 = 0
      Do n = 1, nR
         nLo = nL

         i1 = i2 + 1
         i2 = iW(inER + n - 1)

         Do m = i1, i2
            iE = iW(iIER + m - 1)

            Do 100 j = 1, 6
               iRt = iW(iIRE + 6 * (iE - 1) + j - 1)

               Do k = nLo + 1, nL
                  If(iRt.EQ.IRR(k)) Goto 100
               End do

               nL = nL + 1
               If(nL.GT.MaxL) Then
                  iERR = n
                  Goto 9000
               End if

               IRR(nL) = iRt
 100        Continue
         End do

         nRR(n) = nL
      End do

 9000 Return
      End Subroutine listR2R



C ================================================================
      Subroutine listE2R(nP, nR, nE, IPE, IRE, nEP, IEP)
C ================================================================
C  The routine computes connectivity lists for mesh edges
C
C  nEP(*) - working array of size nP
C  IEP(*) - working array of size 4*nE
C ================================================================
      Integer IPE(4, *), IRE(6, *)
      Integer IEP(*), nEP(*)
C ================================================================
C group (Local variables)
      Integer ipr(2, 6)

      DATA ipr /1,2, 1,3, 1,4, 2,3, 2,4, 3,4/

C ================================================================
      Call backReferences(nP, nE, 4, 4, IPE, nEP, IEP)

      Do n = 1, nE
         Do i = 1, 6
            IRE(i, n) = 0
         End do
      End do

      nR = 0
      Do n = 1, nE
         Do 20 i = 1, 6
            If(IRE(i, n).EQ.0) Then
               nR = nR + 1
               IRE(i, n) = nR

               ip1 = IPE(ipr(1, i), n)
               ip2 = IPE(ipr(2, i), n)
               ip3 = max(ip1, ip2)

               m1 = 1
               If(ip3.GT.1) m1 = nEP(ip3 - 1) + 1

               Do m = m1, nEP(ip3)
                  iE = IEP(m)
                  Do j = 1, 6
                     jp1 = IPE(ipr(1, j), iE)
                     jp2 = IPE(ipr(2, j), iE)

                     If(check22(ip1, ip2, jp1, jp2)) Then
                        If(IRE(j, iE).EQ.0) IRE(j, iE) = nR
                     End if
                  End do
               End do
            End if
 20      Continue
      End do

      Return
      End Subroutine listE2R



C ================================================================
      Subroutine listE2F(nP, nF, nE, IPE, IFE, nEP, IEP)
C ================================================================
C  The routine creates connectivity list E->F for mesh elements
C
C  *** Remarks:
C         1. Working memory is nEP(nP), IEP(4 * nE)
C ================================================================
      Integer IPE(4, *), IFE(4, *)
      Integer IEP(*), nEP(*)

C ================================================================
C group (Local variables)
      Integer iref(5)

      DATA    iref /1,2,3,4,1/

C ================================================================
      Call backReferences(nP, nE, 4, 4, IPE, nEP, IEP)

      Do n = 1, nE
         Do i = 1, 4
            IFE(i, n) = 0
         End do
      End do

      nF = 0
      Do n = 1, nE
         Do 10 i1 = 1, 4
            If(IFE(i1, n).EQ.0) Then
               nF = nF + 1
               IFE(i1, n) = nF

               i2 = iref(i1 + 1)
               i3 = iref(i2 + 1)

               iP1 = IPE(i1, n)
               iP2 = IPE(i2, n)
               iP3 = IPE(i3, n)

               If(cmpE(iP1, iP2, iP3, IEP, nEP, n, iE2)) Then
                  Do j1 = 1, 4
                     j2 = iref(j1 + 1)
                     j3 = iref(j2 + 1)

                     jP1 = IPE(j1, iE2)
                     jP2 = IPE(j2, iE2)
                     jP3 = IPE(j3, iE2)
                     If(check33(iP1, iP2, iP3, jP1, jP2, jP3)) Then
                        IFE(j1, iE2) = nF
                        Goto 10
                     End if
                  End do
               End if
            End if
 10      Continue
      End do

      Return
      End Subroutine listE2F



C ================================================================
      Subroutine listE2Fb(nP, nFb, nE, IPF, IPE, IFE, nEP, IEP)
C ================================================================
C  The routine computes connectivity list E->F for BOUNDARY faces
C
C  *** Remarks:
C         1. Working memory is nEP(nP), IEP(3 * nFb)
C ================================================================
      Integer IPF(3, *), IPE(4, *), IFE(4, *)
      Integer IEP(*), nEP(*)

C ================================================================
C group (Local variables)
      Integer iref(5)

      DATA    iref /1,2,3,4,1/

C ================================================================
      Call backReferences(nP, nFb, 3, 3, IPF, nEP, IEP)

      Do n = 1, nE
         Do i1 = 1, 4
            IFE(i1, n) = 0

            i2 = iref(i1 + 1)
            i3 = iref(i2 + 1)

            iP1 = IPE(i1, n)
            iP2 = IPE(i2, n)
            iP3 = IPE(i3, n)

            If(cmpE(iP1, iP2, iP3, IEP, nEP, 0, iF)) Then
               IFE(i1, n) = iF
            End if
         End do
      End do

      Return
      End Subroutine listE2Fb



C ================================================================
      Subroutine listConv(
     &           nP, nR, nE, nEP, IEP, L, IRE, 
     &           nX, MaxX, nRP, IRP, iW, iERR)
C ================================================================
C  The routine convolutes unstructured map X->Y and structured map 
C  Y->Z to get the map X->Z. For examples, if X means points (P), 
C  Y means elements (E), and Z means edges (R), we get the map from 
C  a point to all edges in the elements having this point.  
C
C  Routine returns 0 upon successful completion.
C
C  *** Remarks:
C         1. Working memory is iW(nR)
C ================================================================
      Integer nP, nR, nE, L, MaxX
      Integer nEP(*), IEP(*), IRE(L, *), nRP(*), IRP(*), iW(*)
C ================================================================
      iERR = 0

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

            Do j = 1, L
               iR = IRE(j, iE)  
               If(iW(iR).GT.0) Then
                  nX = nX + 1
                  If(nX.GT.MaxX) Then
                     iERR = n
                     Goto 9000
                  End if

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

 9000 Return
      End Subroutine listConv



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
      End Subroutine reverseMap



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
      End Subroutine backReferences



C ================================================================
      Subroutine dualNormals(nP, nR, nE, IPE, IPR, XYP, NRM, iW)
C ================================================================
C  The routine computes normals to surfaces of a dual mesh
C  which separate points of a primary mesh.
C
C  Parameters:
C      nP, nR, nE - the number of points, edges and elements in the
C                   primary mesh
C
C      IPE(4, nE) - connectivity list of tetrahedra: E -> P
C      IPR(2, nR) - connectivity list of edges:      R -> P
C
C      XYP(3, nP) - Cartesian coordinates of mesh points
C      NRM(3, nR) - array of unit vectors normal to surfaces of a
C                   dual mesh separating point of the primary mesh
C
C      iW(*) - working memory of size 4 * nE + nP
C ================================================================
      Integer IPE(4, *), IPR(2, *)
      Integer iW(*)

      Real*8  XYP(3, *), NRM(3, *)

C ================================================================
C (Local variables)
      Real*8  XYA(3), XYB(3), XYE(3), XYC(3), XYD(3)
      Real*8  XYM(3), XYN(3), XYR(3)
      Real*8  s
 
C ================================================================
      inEP = 0
      iIEP = inEP + nP

      Call backReferences(nP, nE, 4, 4, IPE, iW(inEP + 1), iW(iIEP + 1))

      Do n = 1, nR
         iP1 = IPR(1, n)
         iP2 = IPR(2, n)

         Do i = 1, 3
            NRM(i, n) = 0D0
            XYE(i) = (XYP(i, iP1) + XYP(i, iP2)) / 2

            XYR(i) = XYP(i, iP2) - XYP(i, iP1)
         End do

         i1 = 1
         If(iP1.GT.1) i1 = iW(inEP + iP1 - 1) + 1

         i2 = iW(inEP + iP1)

         s = 0D0
         Do 10 m = i1, i2
            iE = iW(iIEP + m)

            icnt = 0
            Do j = 1, 4
               iPt = IPE(j, iE)
               If(iPt.NE.iP1 .AND. iPt.NE.iP2) Then
                  icnt = icnt + 1
                  If(icnt.EQ.1) iP3 = iPt
                  If(icnt.EQ.2) iP4 = iPt
               End if
            End do

            If(icnt.NE.2) Goto 10

c   ...   computing vertices of tet differ from iP1 & iP2  
            Do i = 1, 3
               XYC(i) = (XYP(i, iP1) + XYP(i, iP2) 
     &                 + XYP(i, iP3) + XYP(i, iP4)) / 4

               XYA(i) = (XYP(i, iP1) + XYP(i, iP2) 
     &                               + XYP(i, iP3)) / 3 - XYC(i)

               XYB(i) = (XYP(i, iP1) + XYP(i, iP2) 
     &                               + XYP(i, iP4)) / 3 - XYC(i)

               XYD(i) = XYE(i) - XYC(i)
            End do

            Call VecMul(XYD, XYA, XYM)
            Call VecMul(XYD, XYB, XYN)

            s = s + calNorm(XYM) + calNorm(XYN)
            
c    ...    orienting the normal vector from iP1 to iP2
            If(DotMul(XYR, XYM).LT.0D0) Then
               Do i = 1, 3
                  NRM(i, n) = NRM(i, n) - XYM(i)
               End do
            Else
               Do i = 1, 3
                  NRM(i, n) = NRM(i, n) + XYM(i)
               End do
            End if

            If(DotMul(XYR, XYN).LT.0D0) Then
               Do i = 1, 3
                  NRM(i, n) = NRM(i, n) - XYN(i)
               End do
            Else
               Do i = 1, 3
                  NRM(i, n) = NRM(i, n) + XYN(i)
               End do
            End if
 10      Continue

c  ...  rescaling the area of the interface separating iP1 & iP2
c        s = s / calNorm(NRM(1, n))
c        Do i = 1, 3
c           NRM(i, n) = NRM(i, n) * s
c        End do
      End do

      Return
      End Subroutine dualNormals



C ================================================================
      Subroutine addBoundaryFaces(
C ================================================================
     &      nP, nF, MaxF, nE,  
     &      XYP, IPF, IPE, lbF, lbE, 
     &      iW)
C ================================================================
C     iW(*) - working memory of size 2 * nP + 3 * nF + 4 * nE
C ================================================================
      include 'makS.fd'
C ================================================================
      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *), lbF(*), lbE(*)

      Integer iW(*)

C ================================================================
C group (Functions)

C group (Local variables)
      Integer ip(5)
      Logical flagE, flagF

C ================================================================
      iERR = 0
 
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1

      iIFP = 1
      iIEP = iIFP + 3 * nF 
      inFP = iIEP + 4 * nE
      inEP = inFP + nP 

c ... checking that [iVface, MaxS] is clear
      Do n = 1, nF
         If(lbF(n).GE.iVface) Call errMes(1011, 'addBoundaryFaces',
     &                  'reserved boundary identificator is used')
      End do


C ... creating an auxiliary structure
      Call backReferences(nP, nF, 3, 3, IPF, iW(inFP), iW(iIFP))
      Call backReferences(nP, nE, 4, 4, IPE, iW(inEP), iW(iIEP))


C ... creating material and missing boundaries
      k = 0
      kMax = MaxS - iMface
      Do n = 1, nE
         Do i1 = 1, 4
            i2 = ip(i1 + 1)
            i3 = ip(i2 + 1)

            ip1 = IPE(i1, n)
            ip2 = IPE(i2, n)
            ip3 = IPE(i3, n)

            flagF = cmpE(ip1, ip2, ip3, iW(iIFP), iW(inFP), 0, iFt)
            flagE = cmpE(ip1, ip2, ip3, iW(iIEP), iW(inEP), n, iEt)

            If(.NOT.flagF .AND. .NOT.flagE) Then
               nF = nF + 1
               If(nF.GT.MaxF) Call errMes(1004, 'addBoundaryFaces',
     &                            'local parameter MaxF is small')

               IPF(1, nF) = IPE(i1, n)
               IPF(2, nF) = IPE(i2, n)
               IPF(3, nF) = IPE(i3, n)
               lbF(nF) = iVface
            End if
         End do
      End do

      Return
      End Subroutine addBoundaryFaces



C ================================================================
      Subroutine addMaterialFaces(
C ================================================================
     &      nP, nF, MaxF, nE,  
     &      XYP, IPF, IPE, lbF, lbE, 
     &      iW)
C ================================================================
C     iW(*) - working memory of size 2 * nP + 3 * nF + 4 * nE
C ================================================================
      include 'makS.fd'
C ================================================================
      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *), lbF(*), lbE(*)

      Integer iW(*)

C ================================================================
C group (Functions)

C group (Local variables)
      Integer ip(5)
      Logical flagE, flagF

      Integer Mlist(2, MaxS-iMface+1)

C ================================================================
      iERR = 0
 
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1

      iIFP = 1
      iIEP = iIFP + 3 * nF 
      inFP = iIEP + 4 * nE
      inEP = inFP + nP 

c ... checking that [iVface, MaxS] is clear
      Do n = 1, nF
         If(lbF(n).GE.iVface) Call errMes(1011, 'addMaterialFaces',
     &                  'reserved boundary identificator is used')
      End do


C ... creating an auxiliary structure
      Call backReferences(nP, nF, 3, 3, IPF, iW(inFP), iW(iIFP))
      Call backReferences(nP, nE, 4, 4, IPE, iW(inEP), iW(iIEP))


C ... creating material and missing boundaries
      k = 0
      kMax = MaxS - iMface
      Do n = 1, nE
         Do i1 = 1, 4
            i2 = ip(i1 + 1)
            i3 = ip(i2 + 1)

            ip1 = IPE(i1, n)
            ip2 = IPE(i2, n)
            ip3 = IPE(i3, n)

            flagF = cmpE(ip1, ip2, ip3, iW(iIFP), iW(inFP), 0, iFt)
            flagE = cmpE(ip1, ip2, ip3, iW(iIEP), iW(inEP), n, iEt)

            If(.NOT.flagF .AND. flagE .AND. iEt.GT.n) Then
               mat1 = lbE(n)
               mat2 = lbE(iEt)
               If(mat1.NE.mat2) Then
c  ...  searching for this pair in the list of material interfaces
                  Do i = 1, k
                     If((Mlist(1, i).EQ.mat1 .AND.
     &                   Mlist(2, i).EQ.mat2) .OR.
     &                  (Mlist(1, i).EQ.mat2 .AND.
     &                   Mlist(2, i).EQ.mat1)) Then
                        iCface = iMface + i - 1
                        Goto 1
                     End if
                  End do

c  ...  making the new material interface
                  k = k + 1
                  If(k.GT.kMax) Call errMes(1010, 'addMaterialFaces',
     &                        'not enough memory for material faces')
                  Mlist(1, k) = mat1
                  Mlist(2, k) = mat2
                  iCface = iMface + k - 1

 1                nF = nF + 1
                  If(nF.GT.MaxF) Call errMes(1004, 'addMaterialFaces',
     &                               'local parameter MaxF is small')

                  IPF(1, nF) = IPE(i1, n)
                  IPF(2, nF) = IPE(i2, n)
                  IPF(3, nF) = IPE(i3, n)
                  lbF(nF) = iCface
               End if
            End if
         End do
      End do

      Return
      End Subroutine addMaterialFaces



C ================================================================
      Subroutine global2local(
C ================================================================
     &           myID, ICE,
c group (Mg)
     &           nP, nF, MaxF, nE, 
     &           XYP, IPF, IPE, lbF, lbE,
     &           nPv, nFv, nEv, IPV, IFV, IEV,
c group (Ml)
     &           nPl, nFl, nEl, 
     &           XYPl, IPFl, IPEl, lbFl, lbEl,
     &           nPvl, nFvl, nEvl, IPVl, IFVl, IEVl, 
c group (I)
     &           nFvi, IPPl, IFFl, 
c group (W)
     &           IPw, IFw, IEw)
C ================================================================
C  The subroutine extracts a submesh from the global mesh using
C  tetrahedra colored by myID color in array ICE. 
C
C  The fixed triangles are placed at the beginning of corresponding 
C  list. The interfaces triangles are created with utility 
C  addBoundaryFaces and are placed right after the fixed surface 
C  triangles. The number of these triangles equals to nFvi. Note 
C  that these triangles are not added to the list of fixed triangles. 
C
C  The intereface points are also placed at the beginning of list.
C  The user don't need to include them in the list of fixed trianges. 
C  It will be done automatically when the interface triangles will 
C  be added to the list of fixed triangles.
C
C  IPFl(MaxF) - since we need to add interface triangles, nFl may be
C               bigger than nF.
C
C  IPPl(nPl)  - references to a global enumeration necessary to
C               glue submeshes in one global mesh:
C                  0 - point is located outside the DD interfaces;
C                 >0 - splitted points from different subdomains
C                      have the same value of IPPl.
C
C IFFl(nFl)   - references to global enumeration of surface triangles:
C                  0 - the auxiliary interface triangle
C                 >0 - splitted triangles from different meshes have
C                      the same value of IFFl.
C
C  In order to use the information inside IPPl, the interface 
C  triangles should not be modified and be always at the beginning
C  of the corresponding list. The user has to add the interface 
C  triangles (nFvi triangles) to the list of fixed triangles.
C
C  Working memory: IPw(2*nP), IFw(nF), IEw(2*nP + 3*nF + 4*nE)
C
C  Remarks: The number of interface triangles is correct if and 
C  only if the global mesh is complete, i.e. there are no missing 
C  surface triangles. If necessary, the user has to use utility 
C  addBoundaryFaces before splitting the global mesh.
C
C ================================================================
      include 'makS.fd'
      include 'color.fd'
      include 'magic.fd'
C ================================================================
C group (Mg)
      Integer myID, ICE(*)

      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *), lbF(*), lbE(*)
      Integer IPV(*), IFV(*), IEV(*)

C group (Ml)
      Real*8  XYPl(3, *)
      Integer IPEl(4, *), IPFl(3, *), lbFl(*), lbEl(*)
      Integer IPVl(*), IFVl(*), IEVl(*)

C group (I)
      Integer IPPl(*), IFFl(*)

C group (W)
      Integer IPw(*), IFw(*), IEw(*)

c group (Local variables)
      Logical flag 

C ================================================================
      iP1 = 1
      iP2 = iP1 + nP

      Call makTnode(nP, nP, nE, IPw(iP2), IPE, ICE)

c ... marking points of the subdomain with myID
      Do n = 1, nP
         IPw(n) = 0
      End do

      nEl = 0
      Do n = 1, nE
         If(ICE(n).EQ.myID) Then
            nEl = nEl + 1
            IEw(n) = nEl

            Do i = 1, 4
               iP = IPE(i, n)
               IPw(iP) = iP

               IPEl(i, nEl) = iP
            End do

            lbEl(nEl) = lbE(n)
         End if
      End do

      nEvl = 0
      Do n = 1, nEv
         iE = IEV(n)
         If(ICE(iE).EQ.myID) Then
            nEvl = nEvl + 1
            IEVl(nEvl) = IEw(iE)
         End if
      End do


c ... coping points
      nPl = 0
      Do m = 1, 2
         Do n = 1, nP
            iPt = iP2 + n - 1
            flag = (m.EQ.1 .AND. IPw(iPt).GT.0) .OR. 
     &             (m.EQ.2 .AND. IPw(iPt).EQ.0)  

            If(IPw(n).NE.0 .AND. flag) Then
               nPl = nPl + 1
               IPw(n) = nPl

               Do i = 1, 3
                  XYPl(i, nPl) = XYP(i, n)
               End do
                                                                                                                           
               IPPl(nPl) = IPw(iP2 + n - 1)
            End if
         End do
      End do


      nPvl = 0
      Do n = 1, nPv
         If(IPw(IPV(n)).NE.0) Then
            nPvl = nPvl + 1
            IPVl(nPvl) = IPw(IPV(n))
         End if
      End do


C ... coping faces (fixed faces)
      Do n = 1, MaxF
         IFFl(n) = 0
      End do

      Do n = 1, nF
         IFw(n) = 0
      End do

      nFl = 0
      Do 100 n = 1, nFv
         iF = IFV(n)
         Do i = 1, 3
            If(IPw(IPF(i, iF)).EQ.0) Goto 100
         End do

         nFl = nFl + 1
         IFw(iF) = nFl
         Do i = 1, 3
            IPFl(i, nFl) = IPw(IPF(i, iF))
         End do

         lbFl(nFl) = lbF(iF)
         IFFl(nFl) = iF
 100  Continue


c ... converting some surface triangles to interface triangles
      Do n = 1, nP
         IPw(iP2 + n) = IPw(n)
      End do

      Do n = 1, nE
         If(ICE(n).NE.myID) Then 
            Do i = 1, 4
               iP1 = IPE(i, n)
               IPw(iP2 + iP1) = -IPw(iP2 + iP1)
            End do     
         End if
      End do
      
      mFvi = 0
      Do 200 n = 1, nF
         If(IFw(n).GT.0) Goto 200
         Do i = 1, 3
            If(IPw(iP2 + IPF(i, n)).GE.0) Goto 200
         End do

         mFvi = mFvi + 1

         nFl = nFl + 1
         IFw(n) = nFl
         Do i = 1, 3
            IPFl(i, nFl) = IPw(IPF(i, n))
         End do

         lbFl(nFl) = lbF(n)
         IFFl(nFl) = n
 200  Continue
      nFvt = nFl


c ... coping the rest of faces
      Do 300 n = 1, nF
         If(IFw(n).GT.0) Goto 300
         Do i = 1, 3
            If(IPw(IPF(i, n)).EQ.0) Goto 300
         End do

         nFl = nFl + 1
         IFw(n) = nFl
         Do i = 1, 3
            IPFl(i, nFl) = IPw(IPF(i, n))
         End do

         lbFl(nFl) = lbF(n)
         IFFl(nFl) = n
 300  Continue

      nFvl = 0
      Do n = 1, nFv
         iF = IFV(n)
         If(IFw(iF).NE.0) Then
            nFvl = nFvl + 1
            IFVl(nFvl) = IFw(iF)
         End if
      End do


c ... computing local coordinates for elements
      Do n = 1, nEl
         Do i = 1, 4
            IPEl(i, n) = IPw(IPEl(i, n))
         End do
      End do


c ... computing the interface triangles
      nFlo = nFl
      Call addBoundaryFaces(
     &      nPl, nFl, MaxF, nEl,  
     &      XYPl, IPFl, IPEl, lbFl, lbEl, 
     &      IEw)

      nFvi = nFl - nFlo


c ... moving the fixed and interface triangles
      m = nFvl + mFvi + 1  
      Do n = nFl, nFlo + 1, -1
         Do i = 1, 3
            Call swapii(IPFl(i, n), IPFl(i, m))
         End do

         Call swapii(lbFl(n), lbFl(m))
         Call swapii(IFFl(n), IFFl(m))
         m = m + 1
         If(m.GT.nFlo) Goto 500
      End do


c ... changing color of interface triangles (if possible)
 500  nFvi = nFvi + mFvi

      icfree = iVface
      Do 600 ic = 1, iVface - 1
         Do n = 1, nFl
            If(lbFl(n).EQ.ic) Goto 600
         End do
         icfree = ic
         Goto 700
 600  Continue

 700  Continue
      Do n = 1, nFl
         If(lbFl(n).EQ.iVface) lbFl(n) = icfree
      End do

      Return
      End Subroutine global2local



C ================================================================
      Subroutine local2global(
C ================================================================
     &           nMeshes, 
c group (Mg)
     &           MaxP, MaxF, MaxE, 
     &           nP, nF, nE, 
     &           XYP, IPF, IPE, lbF, lbE,
     &           nPv, nFv, nEv, IPV, IFV, IEV,
c group (Ml)
     &           nPl, nFl, nEl, 
     &           XYPl, IPFl, IPEl, lbFl, lbEl,
     &           nPvl, nFvl, nEvl, IPVl, IFVl, IEVl,
c group (I)
     &           IPPl, IFFl, 
c group (W)
     &           IPs, IPw, IFw, IEw)
C ================================================================
C  The subroutine gathers submeshes into a global mesh using
C  references to global enumeration of mesh points kept in array 
C  IPPl(*) and the global enumeration of faces in array IFFl(*).
C
C  The meshes are kept as continuous lists of data. For instance,
C  array XYP(1:3, 1:nP1) keeps coordinates of mesh points in the 
C  first grid. Continuation of the array, XYP(1:3, nP1 : nP1 + nP2),
C  keeps coordinates of the mesh points of the second grid. And
C  so on. The array having such a structure are:
C
C    XYPl(3, *) - mesh coordinates
C    IPFl(3, *) - connectivity list for triangles
C    IPEl(4, *) - connectivity list for tetrahedra
C
C    IPVl(*) - list of fixed points
C    IFVl(*) - list of fixed triangles
C    IEVl(*) - list of fixed tetrahedra  
C
C    IPPl(*) - references to global enumeration of mesh point
C    IFFl(*) - references to global enumeration of mesh triangles 
C
C
C  The necessary information to extract submesh data are kept in
C  the following arrays:
C
C    nPl(nMeshes)  - numbers of points
C    nFl(nMeshes)  - numbers of faces
C    nEl(nMeshes)  - numbers of tetrahedra 
C
C    nPvl(nMeshes) - numbers of fixed points
C    nFvl(nMeshes) - numbers of fixed triangles
C    nEvl(nMeshes) - numbers of fixed tetrahedra 
C
C  Working memory: IPs(MaxP), IPw(MaxP), IFw(MaxF), IEw(4*MaxE)
C 
C  Remark: Array IPFl is destroyed in the algorithm.
C ================================================================
C group (Mg)
      Real*8  XYP(3, *)
      Integer IPF(3, *), IPE(4, *) 
      Integer lbF(*), lbE(*)

      Integer IPV(*), IFV(*), IEV(*)

C group (Ml)
      Integer nPl(*), nFl(*), nEl(*)

      Real*8  XYPl(3, *)
      Integer IPFl(3, *), IPEl(4, *)
      Integer lbFl(*), lbEl(*)

      Integer nPvl(*), nFvl(*), nEvl(*), IPVl(*), IFVl(*), IEVl(*)

C group (I)
      Integer IPPl(*), IFFl(*)

C group (W)
      Integer IPs(*), IPw(*), IFw(*), IEw(*)

c group (Local variables)

C ================================================================
      mP = 0
      mF = 0
      mE = 0     
      Do n = 1, nMeshes
         mP = mP + nPl(n)
         mF = mF + nFl(n)
         mE = mE + nEl(n)
      End do
      If(mP.GT.MaxP) Call errMes(1003, 'local2global',
     &                   'local parameter MaxP is small')
      If(mF.GT.MaxF) Call errMes(1004, 'local2global',
     &                   'local parameter MaxF is small')
      If(mE.GT.MaxE) Call errMes(1006, 'local2global',
     &                   'local parameter MaxE is small')

      Do n = 1, mP
         IPs(n) = 0
         IPw(n) = 0
      End do

c ... counting points on interfeices 
      kP = 0
      Do n = 1, mP
         If(IPPl(n).NE.0) Then
            If(IPs(IPPl(n)).EQ.0) Then
               kP = kP + 1
               IPs(IPPl(n)) = kP
            End if
         End if
      End do

c ... counting the rest of points
      Do n = 1, mP
         If(IPPl(n).EQ.0) Then
            kP = kP + 1
            IPw(n) = kP
         Else
            IPw(n) = IPs(IPPl(n))
         End if
      End do


c ... making points
      nP = kP
      Do n = 1, mP
         Do i = 1, 3
            XYP(i, IPw(n)) = XYPl(i, n)
         End do
      End do


c ... making vertices
      Do n = 1, nP
         IPs(n) = 0
      End do

      mP = 0
      mV = 0
      Do k = 1, nMeshes
         Do n = 1, nPvl(k)
            mV = mV + 1
            IPs(IPw(mP + IPVl(mV))) = 1
         End do
         mP = mP + nPl(k)
      End do

      nPv = 0
      Do n = 1, nP
         If(IPs(n).NE.0) Then
            nPv = nPv + 1
            IPV(nPv) = n
         End if
      End do


c ... making faces
      mP = 0
      mF = 0
      Do k = 1, nMeshes
         Do n = 1, nFl(k)
            mF = mF + 1

            IFw(mF) = 0

            Do i = 1, 3
               IPFl(i, mF) = IPw(mP + IPFl(i, mF))
            End do
         End do
         mP = mP + nPl(k)
      End do

      Call backReferences(nP, mF, 3, 3, IPFl, IPs, IEw)

      Do n = 1, mF
         If(IFw(n).EQ.0) Then
            ip1 = IPFl(1, n)
            ip2 = IPFl(2, n)
            ip3 = IPFl(3, n)

            If(cmpE(ip1, ip2, ip3, IEw, IPs, n, iF)) Then
               If(iF.GT.n) IFw(n) = iF
            End if
         End if
      End do

      nF = 0
      Do n = 1, mF
         If(IFw(n).EQ.0 .AND. IFFl(n).GT.0) Then
            nF = nF + 1
            IFw(n) = nF

            Do i = 1, 3
               IPF(i, nF) = IPFl(i, n)
            End do
            lbF(nF) = lbFl(n)
         Else
            IFw(n) = 0
         End if
      End do


c ... making fixed faces
      Do n = 1, nF
         IEw(n) = 0
      End do

      mF = 0
      mV = 0
      Do k = 1, nMeshes
         Do n = 1, nFvl(k)
            mV = mV + 1
            IEw(IFw(mF + IFVl(mV))) = 1
         End do
         mF = mF + nFl(k)
      End do

      nFv = 0
      Do n = 1, nF
         If(IEw(n).NE.0) Then
            nFv = nFv + 1
            IFV(nFv) = n
         End if
      End do


c ... making elements
      mP = 0
      nE = 0
      Do k = 1, nMeshes
         Do n = 1, nEl(k)
            nE = nE + 1
            
            Do i = 1, 4
               IPE(i, nE) = IPw(mP + IPEl(i, nE))
            End do
            lbE(nE) = lbEl(nE)
         End do
         mP = mP + nPl(k)
      End do


c ... making fixed elements
      nEv = 0

      mE = 0
      mV = 0
      Do k = 1, nMeshes
         Do n = 1, nEvl(k)
            nEv = nEv + 1
            IEV(nEv) = mE + IEVl(mV + n)
         End do
         mE = mE + nEl(k)
         mV = mV + nEvl(k)
      End do

      Return
      End Subroutine local2global



C ================================================================
      Subroutine makTnode(nP, MaxP, nE, IPP, IPE, ICE)
C ================================================================
C   The T-nodes are marked. The original connectivity
C   list IPE is used. 
C ================================================================
      Integer IPP(*), IPE(4, *), ICE(*)

      Do n = 1, MaxP
         IPP(n) = 0
      End do

      Do n = 1, nE
         ict = ICE(n) + 1
         Do i = 1, 4
            iP = IPE(i, n)
            icn = IPP(iP)
            If(icn.EQ.0) Then
               IPP(iP) = ict
            Else If(icn.NE.ict) Then
               IPP(iP) = -1
            End if
         End do
      End do


      Do n = 1, nP
         If(IPP(n).EQ.-1) Then
            IPP(n) = n
         Else
            IPP(n) = 0
         End if
      End do

      Return
      End Subroutine makTnode
C
      End Module mba3d_utils
