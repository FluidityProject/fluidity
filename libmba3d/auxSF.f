      Module mba3d_auxSF
C
      use mba3d_auxSE
      use mba3d_auxSR
      use mba3d_error
      use mba3d_findSE
      use mba3d_makQ
C
      contains
C
C ================================================================
      Subroutine makSF(
C ================================================================
     &           iP, lF, iFs, lE, iEs, IPEs, IPF, IFE, MaxS,
     &           lP1, iP1s, lP2, iP2s, ICP2s,
     &           lS, IPSs, IESs, lDF, iDFs, lDE, iDEs, 
     &           flagFACE, flagEDGE)
C ================================================================
C Compute:
C   (a) list of faces IPSs(3, lS) around point iP is created.
C
C   (b) list of surface points connected with point iP 
C
C   (c) edge case : list of points on an edge passing through iP
C       face case : list of points on the face passing through iP
C ================================================================
      Integer iFs(*), iEs(*)
      Integer IPEs(5, *), IPF(4, *), IFE(4, *)
      Integer iP1s(*), iP2s(*), ICP2s(*)
      Integer IPSs(3, *), IESs(*), iDFs(*), iDEs(*)
      Logical flagFACE, flagEDGE

C group (Local variables)
      Integer iref(5), iPt(3)
      Logical flag

C ================================================================
      iref(1) = 1
      iref(2) = 2
      iref(3) = 3
      iref(4) = 4
      iref(5) = 1


      lP1 = 0
      lP2 = 0

      lS  = 0

      lDF = 0
      lDE = 0


      flagFACE = .FALSE.
      flagEDGE = .FALSE.

      Do n = 1, lE
         flag = .FALSE.
         Do 10 i1 = 1, 4
            If(IPEs(i1, n).EQ.iP) Then
               flag = .TRUE.

               lDE = lDE + 1
               iDEs(lDE) = n

               i2 = i1
               Do i = 1, 3
                  i2 = iref(i2 + 1)
                  iPt(i) = IPEs(i2, n)

                  Call findSE(lP2, iP2s, iPt(i), nPt)

                  If(nPt.EQ.0) Then
                     lP2 = lP2 + 1
                     If(lP2.GT.MaxS) Goto 1000

                     iP2s(lP2) = iPt(i)
                     ICP2s(lP2) = 0
                  End if
               End do


               Do m = 1, lS
                  If(check33(iPt(1), iPt(2), iPt(3),
     &               IPSs(1, m), IPSs(2, m), IPSs(3, m))) Goto 10
               End do

               lS = lS + 1
               If(lS.GT.MaxS) Goto 1000

               Do i = 1, 3
                  IPSs(i, lS) = iPt(i)
               End do
               IESs(lS) = n

               Goto 15
            End if
 10     Continue


 15     If(flag) Then
           Do 20 i1 = 1, 4
              iFt = IFE(i1, iEs(n))

              If(iFt.NE.0) Then
                 i2 = iref(i1 + 1)
                 i3 = iref(i2 + 1)

                 iPt(1) = IPEs(i1, n)
                 iPt(2) = IPEs(i2, n)
                 iPt(3) = IPEs(i3, n)

                 If(check13(iP, iPt(1), iPt(2), iPt(3))) Then
                    iCLRs = IPF(4, iFt)
                    If(iCLRs.LT.0) Call errMes(6001, 'auxSF',
     &                                        'system error')

                    Do i = 1, 3
                       Call findSE(lP2, iP2s, iPt(i), nPt)
                       If(nPt.GT.0) Then
                          If(ICP2s(nPt).EQ.0) Then
                             flagFACE = .TRUE.

                             ICP2s(nPt) = iCLRs
                          Else If(ICP2s(nPt).NE.iCLRs) Then
                             flagEDGE = .TRUE.

                             Call findSE(lP1, iP1s, iPt(i), mPt)
                             If(mPt.EQ.0) Then
                                lP1 = lP1 + 1
                                If(lP1.GT.MaxS) Goto 1000

                                iP1s(lP1) = iPt(i)
                             End if

                             ICP2s(nPt) = ICP2s(nPt) + iCLRs
                           End if
                       End if
                    End do

                    Call findSE(lF, iFs, iFt, nFt)
                    Call findSE(lDF, iDFs, nFt, k)
                    If(k.LE.0) Then
                       lDF = lDF + 1
                       If(lDF.GT.MaxS) Goto 1000

                       iDFs(lDF) = nFt
                    End if
                 End if
              End if
 20        Continue
        End if
      End do


      If(flagEDGE) Then
         Goto 9000
      Else If(flagFACE) Then
         Do n = 1, lP2
            If(ICP2s(n).NE.0) Then
               lP1 = lP1 + 1
               iP1s(lP1) = iP2s(n)
            End if
         End do
      Else
         lP1 = lP2
         Do n = 1, lP1
            iP1s(n) = iP2s(n)
         End do
      End if

 9000 Return
 1000 Call errMes(1007, 'makSF', 'local variable MaxS is small')
      End Subroutine makSF



C ================================================================
      Subroutine shutF(ZZ1, ZZ2, xyp1, xyp2, xyp3, flag)
C ================================================================
C Routine returns .TRUE. if the ray [ZZ1, ZZ2) intersects
C triangle {xyz1, xyz2, xyz3}. If flag = .TRUE. the triangle is
C assumed to be the open domain.
C ================================================================
      Real*8  ZZ1(3),  ZZ2(3)
      Real*8  xyp1(3), xyp2(3), xyp3(3)
      Logical flag

C group (Local variables)
      Real*8  v1, v2
      Logical flagOPEN

C ================================================================
      flagOPEN = flag
      flag = .FALSE.


      v1 = calVol(ZZ1, ZZ2, xyp3, xyp1)
      v2 = calVol(ZZ1, ZZ2, xyp3, xyp2)

      If(flagOPEN) Then
         If(v1 * v2.GE.0D0) Goto 1000
      Else
         If(v1 * v2.GT.0D0) Goto 1000
      End if


      v1 = calVol(ZZ1, ZZ2, xyp2, xyp1)
      v2 = calVol(ZZ1, ZZ2, xyp2, xyp3)

      If(flagOPEN) Then
         If(v1 * v2.GE.0D0) Goto 1000
      Else
         If(v1 * v2.GT.0D0) Goto 1000
      End if


      v1 = calVol(ZZ1, ZZ2, xyp1, xyp2)
      v2 = calVol(ZZ1, ZZ2, xyp1, xyp3)

      If(flagOPEN) Then
         If(v1 * v2.GE.0D0) Goto 1000
      Else
         If(v1 * v2.GT.0D0) Goto 1000
      End if

      flag = .TRUE.
 1000 Return
      End Subroutine shutF



C ================================================================
      Real*8 Function projectF(xyp, xy1, xy2, xy3, xyt)
C ================================================================
C Routine finds a point xyt lying inside triangular face 
C {xy1, xy2, xy3} which is the closest point to xyp. It returns
C the distance to this point.
C
C *** Remarks:
C        1. We minimize the distance ||a * x21 + b * x31 - p||^2
C           w.r.t. a and b. 
C ================================================================
      Real*8  xyp(3), xy1(3), xy2(3), xy3(3), xyt(3)

c group (Local variables)
      Real*8  v1(3), v2(3), vp(3)
      Real*8  m11, m12, m22, f1, f2, det, a, b 

C ================================================================
      Do i = 1, 3
         v1(i) = xy2(i) - xy1(i)
         v2(i) = xy3(i) - xy1(i)
         vp(i) = xyp(i) - xy1(i)
      End do

      m11 = DotMul(v1, v1)
      m12 = DotMul(v1, v2)
      m22 = DotMul(v2, v2)

      f1 = DotMul(vp, v1)
      f2 = DotMul(vp, v2)


c ... solve a linear system (Wronskian is always non-zero)
      det = m22 * m11 - m12 * m12

      a = (m22 * f1 - m12 * f2) / det
      b = (m11 * f2 - m12 * f1) / det

c ... analyze few cases (seach along b, a and a + b = 1)
      If(a.LT.0D0) Then
         a = 0D0
         b = max(0D0, f2 / m22)
         b = min(1D0, b)
      Else If(b.LT.0D0) Then
         b = 0D0
         a = max(0D0, f1 / m11)
         a = min(1D0, a)
      Else If(a + b.GT.1D0) Then
         a = max(0D0, (f1 - m12) / (m11 - m12))
         a = min(1D0, a)
         b = 1D0 - a
      End if

 
c ... compute point the closest point
      Do i = 1, 3
         xyt(i) = xy1(i) + a * v1(i) + b * v2(i)
      End do
      
      projectF = calEdge(xyp, xyt)

      Return
      End Function projectF



C ================================================================
      Subroutine bariCoords(xyp, xy1, xy2, xy3, b1, b2, b3)
C ================================================================
C Routine computes barycentric coordinates of point xyp lying
C inside triangle {xy1, xy2, xy3}. They coordinates are always 
C positive and sum to 1 upto round-off errors.
C ================================================================
      Real*8  xyp(3), xy1(3), xy2(3), xy3(3)
      Real*8  b1, b2, b3

c group (Local variables)
      Real*8  a(3), b(3), c(3), art, ar1

C ================================================================
c ... square of area of the triangle
      Do i = 1, 3
         a(i) = xy2(i) - xy1(i)
         b(i) = xy3(i) - xy1(i)
      End do

      Call VecMul(a, b, c)
      art = c(1) ** 2 + c(2) ** 2 + c(3) ** 2


c ... first baricentric coordinate
      Do i = 1, 3
         a(i) = xy2(i) - xyp(i)
         b(i) = xy3(i) - xyp(i)
      End do

      Call VecMul(a, b, c)
      ar1 = c(1) ** 2 + c(2) ** 2 + c(3) ** 2

      b1 = dsqrt(ar1 / art)


c ... second baricentric coordinate
      Do i = 1, 3
         a(i) = xy1(i) - xyp(i)
      End do

      Call VecMul(a, b, c)
      ar1 = c(1) ** 2 + c(2) ** 2 + c(3) ** 2

      b2 = dsqrt(ar1 / art)


c ... third baricentric coordinate
      b3 = 1 - b2 - b1
      If(b1 + b2.GT.1D0) Then
         b1 = min(1D0, b1)
         b2 = max(0D0, 1 - b1) 
         b3 = 0D0
      End if

      Return
      End Subroutine bariCoords



C ================================================================
      Real*8 Function angle2Faces( XYA, XYB, XYC, XYD )
C ================================================================
C Cosine of angle between two faces {A,B,C} and {A,B,D}.
C
C Remark: the angle is equal to -1 for flat faces.
C ================================================================
      Real*8  XYA(3), XYB(3), XYC(3), XYD(3)
      Real*8  V1(3), V2(3), V3(3), XYN(3), XYM(3)
      Real*8  d1, d2

      Do i = 1, 3
         V1(i) = XYB(i) - XYA(i)
         V2(i) = XYC(i) - XYB(i)
         V3(i) = XYD(i) - XYB(i)
      End do

      Call VecMul(V1, V2, XYN)
      Call VecMul(V1, V3, XYM)

      d1 = XYN(1) ** 2 + XYN(2) ** 2 + XYN(3) ** 2
      d2 = XYM(1) ** 2 + XYM(2) ** 2 + XYM(3) ** 2

      angle2Faces = (XYN(1) * XYM(1) +
     &               XYN(2) * XYM(2) +
     &               XYN(3) * XYM(3)) / dsqrt(d1 * d2)

      Return
      End Function angle2Faces



C ================================================================
      Subroutine VecMul(a, b, c)
C ================================================================
      Real*8 a(3), b(3), c(3)

      c(1) = a(2) * b(3) - b(2) * a(3)
      c(2) = b(1) * a(3) - a(1) * b(3)
      c(3) = a(1) * b(2) - b(1) * a(2)
      Return
      End Subroutine VecMul



C ================================================================
      Real*8 Function DotMul(a, b)
C ================================================================
      Real*8 a(3), b(3)

      DotMul = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)
      Return
      End Function DotMul
C
      End Module mba3d_auxSF
