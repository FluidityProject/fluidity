      Module mba3d_auxSE
C
      use mba3d_auxSP
      use mba3d_calVol
      use mba3d_error
C
      contains
C ================================================================
      Subroutine makSE(
C ================================================================
     &           iE, ICP, IEP, IPF, IPE, IFE, IEE, qE, MaxS,
     &           lF, lE,  iFs, iEs, IPFs, IPEs, qEs)
C ================================================================
C A superelement around element iE is created.
C
C Remark: The element iE and his face-heighboors should be included 
C         independently to prevent a situation where vertices of iE
C         give a multy-connected (point-based) superelement.
C         In addition, we repeat algorithm whenever the "face closure"
C         of the superlement touches element iE.
C
C Remark: The 5th column of the IPE is temporary overloaded here.
C         The negative values of IPE(5, *) are used for a quick
C         search in the list iEs.
C
C Remark: Four elements of ICP are temporary overloded here.
C         The negative values of ICP(*) are used for a quick
C         search for vertices of iE.
C
C ================================================================
c group (M)
      Integer ICP(*), IEP(*), IPF(4, *)
      Integer IPE(5, *), IFE(4, *), IEE(4, *)
      Real*8  qE(*)

c group (S)
      Integer iFs(*), iEs(*), IPFs(4, *), IPEs(5, *)
      Real*8  qEs(*)

C ================================================================
      kE = 1
      iEs(1) = iE
      Do i = 1, 4
         iEt = IEE(i, iE)
         If(iEt.GT.0) Then
            kE = kE + 1
            iEs(kE) = iEt
         End if

         ICP(IPE(i, iE)) = -ICP(IPE(i, iE))
      End do

      nMaxS = MaxS - 1
      Do i = 1, 4
         iP = IPE(i, iE)
         Call makSP(iP, IEP, IPE, IEE, nMaxS, kEadd, iEs(kE + 1))

         kE = kE + kEadd
         nMaxS = nMaxS - kEadd
      End do

      lE = 0
      Do 5 n = 1, kE
         iEt = iEs(n)
         If(IPE(5, iEt).LT.0) Goto 5

         lE = lE + 1
         iEs(lE) = iEt

         Do i = 1, 5
            IPEs(i, lE) = IPE(i, iEt)
         End do

         IPE(5, iEt) = -IPE(5, iEt)
  5   Continue


C ... 3D extention: adding of elements adajcent by a face
      nMaxS = MaxS

      k1 = 1
      kE = lE
      ks = 1

 10   m = 0
      Do n = k1, kE, ks
         If(lE.GE.nMaxS - 4) Goto 1000

         Do 15 i = 1, 4
            iEt = IEE(i, iEs(n))
            If(iEt.NE.0) Then
               If(IPE(5, iEt).LT.0) Goto 15

               lE = lE + 1
               iEs(lE) = iEt

               Do j = 1, 5
                  IPEs(j, lE) = IPE(j, iEt)
               End do

               IPE(5, iEt) = -IPE(5, iEt)

c  ...  checking for intersection with tetrahedron iE
               Do j = 1, 4
                  If(ICP(IPE(j, iEt)).LT.0) Then
                     m = m + 1
                     iEs(MaxS - m + 1) = iEt 

                     nMaxS = nMaxS - 1
                     Goto 15
                  End if
               End do 
            End if
 15      Continue
      End do

      If(m.GT.0) Then
         k1 = MaxS
         kE = MaxS - m + 1
         ks = -1

         Goto 10
      End if


c ... restoring the 5th colunm of IPE
      Do n = 1, lE
         IPE(5, iEs(n)) = -IPE(5, iEs(n))
      End do


c ... restoring the ICP
      Do i = 1, 4
         ICP(IPE(i, iE)) = -ICP(IPE(i, iE))
      End do


      lF = 0
      Do n = 1, lE
         kE = iEs(n)
         qEs(n) = qE(kE)

         Do 20 i = 1, 4
            iF = IFE(i, kE)
            If(iF.NE.0) Then
               Do m = 1, lF
                  If(iF.EQ.iFs(m)) Goto 20
               End do

               lF = lF + 1
               iFs(lF) = iF

               Do j = 1, 4
                  IPFs(j, lF) = IPF(j, iF)
               End do
            End if
 20      Continue
      End do

      Return
 1000 Call errMes(1007, 'makSE', 'local variable MaxS is small')
      End Subroutine makSE



C ================================================================
      Subroutine calSO(XYP, IPE, lE, iEs, iOs)
C ================================================================
C Oriented volumes of tetrahedra are saved in array iOs.
C ================================================================
      Real*8  XYP(3, *)
      Integer IPE(5, *), iEs(*), iOs(*)

      Real*8  d

      Do 10 n = 1, lE
         iE = iEs(n)
         If(iE.LE.0) Goto 10

         iP1 = IPE(1, iE)
         iP2 = IPE(2, iE)
         iP3 = IPE(3, iE)
         iP4 = IPE(4, iE)

         d = calVol(XYP(1, iP1), XYP(1, iP2),
     &              XYP(1, iP3), XYP(1, iP4))
         iOs(n) = dsign(1D0, d)
 10   Continue
      Return
      End Subroutine calSO



C ================================================================
      Subroutine chkSO(iP, XYPs, XYP, IPE, lE, iEs, iOs, flag)
C ================================================================
C flag = TRUE if orientation of new tetrahedra coincides with
C             the orientation given in array iOs
C ================================================================
      Real*8  XYPs(3), XYP(3, *)
      Integer IPE(5, *), iEs(*), iOs(*)
      Logical flag

      Real*8  d


      flag = .FALSE.
      Do 10 n = 1, lE
         iE = iEs(n)
         If(iE.LE.0) Goto 10

         iP1 = IPE(1, iE)
         iP2 = IPE(2, iE)
         iP3 = IPE(3, iE)
         iP4 = IPE(4, iE)

         If(iP1.EQ. iP) Then
            d = calVol(XYPs, XYP(1, iP2),
     &                 XYP(1, iP3), XYP(1, iP4))
         Else If(iP2.EQ.iP) Then
            d = calVol(XYP(1, iP1), XYPs,
     &                 XYP(1, iP3), XYP(1, iP4))
         Else If(iP3.EQ.iP) Then
            d = calVol(XYP(1, iP1), XYP(1, iP2),
     &                 XYPs, XYP(1, iP4))
         Else If(iP4.EQ.iP) Then
            d = calVol(XYP(1, iP1), XYP(1, iP2),
     &                 XYP(1, iP3), XYPs)
         End if

         If(d.EQ.0D0) Goto 1000

         iOt = dsign(1D0, d)
         If(iOt.NE.iOs(n)) Goto 1000
 10   Continue

      flag = .TRUE.
 1000 Return
      End Subroutine chkSO



C ================================================================
      Subroutine copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &                  lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)
C ================================================================
C Copy superelement data (u) into superelement data (s).
C ================================================================
      Integer iFu(*), iEu(*), IPFu(4, *), IPEu(5, *)
      Integer iFs(*), iEs(*), IPFs(4, *), IPEs(5, *)
      Real*8  qEu(*), qEs(*)

      lF = lFu
      Do n = 1, lF
         iFs(n) = iFu(n)

         Do i = 1, 4
            IPFs(i, n) = IPFu(i, n)
         End do
      End do

      lE = lEu
      Do n = 1, lE
         iEs(n) = iEu(n)
         qEs(n) = qEu(n)

         Do i = 1, 5
            IPEs(i, n) = IPEu(i, n)
         End do
      End do
      Return
      End Subroutine copySE



C ================================================================
      Subroutine copySEmove(lFu, lEu, iFu, iEu, qEu,
     &                      lF,  lE,  iFs, iEs, qEs)
C ================================================================
C A version of copySE optimized for local operation MOVE.
C ================================================================
      Integer iFu(*), iEu(*)
      Integer iFs(*), iEs(*)
      Real*8  qEu(*), qEs(*)

      lF = lFu
      Do n = 1, lF
         iFs(n) = iFu(n)
      End do

      lE = lEu
      Do n = 1, lE
         iEs(n) = iEu(n)
         qEs(n) = qEu(n)
      End do
      Return
      End Subroutine copySEmove



C ================================================================
      Subroutine copySQ(lE, qEu, XYPu, HesPu, detGu,
     &                      qEs, XYPs, HesPs, detGs)
C ================================================================
C Routine copies data from u-arrays to s-arrays.
C ================================================================
      Real*8  qEu(*), XYPu(3), HesPu(6), detGu
      Real*8  qEs(*), XYPs(3), HesPs(6), detGs

      Do n = 1, lE
         qEs(n) = qEu(n)
      End do

      Do i = 1, 3
         XYPs(i) = XYPu(i)
      End do

      Do i = 1, 6
         HesPs(i) = HesPu(i)
      End do

      detGs = detGu
      Return
      End Subroutine copySQ



C ================================================================
      Logical Function check1j(i, j)
C ================================================================
C check1j = TRUE if i belongs to the set j(4)}.
C ================================================================
      Integer i, j(4)

      check1j = .TRUE.
      If(i.EQ.j(1)) Goto 1000
      If(i.EQ.j(2)) Goto 1000
      If(i.EQ.j(3)) Goto 1000
      If(i.EQ.j(4)) Goto 1000

      check1j = .FALSE.
 1000 Return
      End Function check1j



C ================================================================
      Logical Function check2j(i1, i2, j)
C ================================================================
C check2j = TRUE if {i1, i2} is a subset of j(3).
C ================================================================
      Integer j(3)

      check2j = .FALSE.
      If(i1.NE.j(1) .AND. i1.NE.j(2) .AND. i1.NE.j(3)) Goto 1000
      If(i2.NE.j(1) .AND. i2.NE.j(2) .AND. i2.NE.j(3)) Goto 1000

      check2j = .TRUE.
 1000 Return
      End Function check2j



C ================================================================
      Logical Function check22(i1, i2, j1, j2)
C ================================================================
C check22 = TRUE if pair {i1, i2} coinsides with {j1, j2}.
C ================================================================
      check22 = .FALSE.
      If(i1.NE.j1 .AND. i1.NE.j2) Goto 1000
      If(i2.NE.j1 .AND. i2.NE.j2) Goto 1000

      check22 = .TRUE.
 1000 Return
      End Function check22



C ================================================================
      Logical Function check3j(i1, i2, i3, j)
C ================================================================
C check3j = TRUE if {i1, i2, i3} is a subset of j(4).
C ================================================================
      Integer j(4)

      If(i1.EQ.j(1) .OR. i1.EQ.j(2) .OR.
     &   i1.EQ.j(3) .OR. i1.EQ.j(4)) Then
         If(i2.EQ.j(1) .OR. i2.EQ.j(2) .OR.
     &      i2.EQ.j(3) .OR. i2.EQ.j(4)) Then
            If(i3.EQ.j(1) .OR. i3.EQ.j(2) .OR.
     &         i3.EQ.j(3) .OR. i3.EQ.j(4)) Then
               check3j = .TRUE.
               Goto 1000
            End if
         End if
      End if

      check3j = .FALSE.
 1000 Return
      End Function check3j



C ================================================================
      Logical Function check33(i1, i2, i3, j1, j2, j3)
C ================================================================
C check33 = TRUE if triple {i1, i2, i3} coinsides with {j1, j2, j3}
C ================================================================
      If(i1.EQ.j1 .OR. i1.EQ.j2 .OR. i1.EQ.j3) Then
         If(i2.EQ.j1 .OR. i2.EQ.j2 .OR. i2.EQ.j3) Then
            If(i3.EQ.j1 .OR. i3.EQ.j2 .OR. i3.EQ.j3) Then
               check33 = .TRUE.
               Goto 1000
            End if
         End if
      End if

      check33 = .FALSE.
 1000 Return
      End Function check33



C ================================================================
      Logical Function check13(i1, j1, j2, j3)
C ================================================================
C check13 = TRUE if i1 belongs to the set {j1, j2, j3}.
C ================================================================
      check13 = .TRUE.
      If(i1.EQ.j1) Goto 1000
      If(i1.EQ.j2) Goto 1000
      If(i1.EQ.j3) Goto 1000

      check13 = .FALSE.
 1000 Return
      End Function check13



C ================================================================
      Logical Function check14(i1, j1, j2, j3, j4)
C ================================================================
C check14 = TRUE if i1 belongs to the set {j1, j2, j3, j4}.
C ================================================================
      check14 = .TRUE.
      If(i1.EQ.j1) Goto 1000
      If(i1.EQ.j2) Goto 1000
      If(i1.EQ.j3) Goto 1000
      If(i1.EQ.j4) Goto 1000

      check14 = .FALSE.
 1000 Return
      End Function check14



C ================================================================
      Subroutine swapdd(d1, d2)
C ================================================================
C Routine swaps two real*8 numbers.
C ================================================================
      Real*8  d1, d2, d

      d  = d1
      d1 = d2
      d2 = d

      Return
      End Subroutine swapdd



C ================================================================
      Subroutine swapii(i1, i2)
C ================================================================
C Routine swaps two integer numbers.
C ================================================================
      i  = i1
      i1 = i2
      i2 = i

      Return
      End Subroutine swapii
C
      End Module mba3d_auxSE
