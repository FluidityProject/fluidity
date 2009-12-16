      Module mba3d_update
C
      use mba3d_auxSE
      use mba3d_error
C
      contains
C
C ==========================================================
C PONITS ... POINT ... POINTS ... POINT
C ==========================================================
      Subroutine pntAdd(iP, nP, MaxP, ICP,  XYP,  HesP,  detG, IHolP,
     &                                ICPs, XYPs, HesPs, detGs)
C ==========================================================
      Integer ICP(*), IHolP(*)
      Real*8  XYP(3, *),  HesP(6, *),  detG(*)
      Real*8  XYPs(3, *), HesPs(6, *), detGs

      nHolP = IHolP(1)
      nP = nP + 1
      If(nP.GT.MaxP) Call errMes(1003, 'pntAdd',
     &                   'local parameter MaxP is small')

      If(nHolP.EQ.0) Then
         iP = nP
      Else
         iP = IHolP(nHolP + 1)
         nHolP = nHolP - 1
         IHolP(1) = nHolP
      End if

      Call pntUpd(iP, ICP, XYP, HesP, detG, ICPs, XYPs, HesPs, detGs)

      Return
      End Subroutine pntAdd


C ==========================================================
      Subroutine pntUpd(iP, ICP,  XYP,  HesP,  detG,
     &                      ICPs, XYPs, HesPs, detGs)
C ==========================================================
      Integer ICP(*)
      Real*8  XYPs(3), XYP(3, *)
      Real*8  HesPs(6), HesP(6, *)
      Real*8  detGs, detG(*)

      ICP(iP) = ICPs

      Do i = 1, 3
         XYP(i, iP) = XYPs(i)
      End do

      Do i = 1, 6
         HesP(i, iP) = HesPs(i)
      End do

      detG(iP) = detGs
      Return
      End Subroutine pntUpd



C ==========================================================
      Subroutine pntDel(iP, nP, ICP, IHolP)
C ==========================================================
      Integer ICP(*), IHolP(*)

      nHolP = IHolP(1)
      nHolP = nHolP + 1
      IHolP(nHolP + 1) = iP

      ICP(iP) = -ICP(iP)

      IHolP(1) = nHolP
      nP = nP - 1
      Return
      End Subroutine pntDel


C ==========================================================
C FACES ... FACE ... FACES ... FACE
C ==========================================================
      Subroutine facAdd(iF, nF, MaxF, IHolF)
C ==========================================================
      Integer IHolF(*)

      nHolF = IHolF(1)
      nF = nF + 1
      If(nF.GT.MaxF) Call errMes(1004, 'facAdd',
     &                   'local parameter MaxF is small')

      If(nHolF.EQ.0) Then
         iF = nF
      Else
         iF = IHolF(nHolF + 1)
         nHolF = nHolF - 1
         IHolF(1) = nHolF
      End if

      Return
      End Subroutine facAdd


C ==========================================================
      Subroutine facUpd(nFs, IPF, iFs, IPFs)
C ==========================================================
      Integer IPF(4, *)
      Integer iFs(*), IPFs(4, *)

      iF = iFs(nFs)
      Do i = 1, 4
         IPF(i, iF) = iPFs(i, nFs)
      End do

      Return
      End Subroutine facUpd


C ==========================================================
      Subroutine facDel(iF, nF, IPF, IHolF)
C ==========================================================
      Integer IPF(4, *), IHolF(*)

      nHolF = IHolF(1)
      nHolF = nHolF + 1
      IHolF(nHolF + 1) = iF

      Do i = 1, 4
         IPF(i, iF) = -iabs(IPF(i, iF))
      End do

      IHolF(1) = nHolF
      nF = nF - 1
      Return
      End Subroutine facDel


C ==========================================================
C ELEMENTS ... ELEMENT ... ELEMENTS ... ELEMENT
C ==========================================================
      Subroutine eleAdd(nE, MaxE, IHolE)
C ==========================================================
      Integer IHolE(*)

      If(nE.GE.MaxE .AND. IHolE(1).EQ.0)
     &     Call errMes(1006, 'eleAdd',
     &                'local parameter MaxE is small')

      Return
      End Subroutine eleAdd


C ==========================================================
      Subroutine eleUpd(nEs, IEP, IPE, IFE,  IEE,
     &                  lF,  lE,  iFs, iEs, IPFs, IPEs)
C ==========================================================
      include 'makS.fd'
      include 'color.fd'

      Integer IEP(*), IPE(5, *), IFE(4, *), IEE(4, *)
      Integer IPFs(4, *), IPEs(5, *)
      Integer iFs(*), iEs(*)

C group (Local variables)
      Integer iEu(MaxS), iref(5)

C ==========================================================
      iref(1) = 1
      iref(2) = 2
      iref(3) = 3
      iref(4) = 4
      iref(5) = 1

      iE = iEs(nEs)

      Do i = 1, 4
         IPE(i, iE) = IPEs(i, nEs)
         IEP(IPE(i, iE)) = iE
         IFE(i, iE) = 0
      End do
      IPE(5, iE) = IPEs(5, nEs)

      Do 30 i1 = 1, 4
         i2 = iref(i1 + 1)
         i3 = iref(i2 + 1)

         iP1 = IPEs(i1, nEs)
         iP2 = IPEs(i2, nEs)
         iP3 = IPEs(i3, nEs)

         Do n = 1, lF
            If(iFs(n).GT.0) Then
               jP1 = IPFs(1, n)
               jP2 = IPFs(2, n)
               jP3 = IPFs(3, n)

               If(check33(iP1, iP2, iP3, jP1, jP2, jP3)) Then
                  IFE(i1, iE) = iFs(n)
                  Goto 10
               End if
            End if
         End do

 10      If(IFE(i1, iE).NE.0) Then
C  ...  analyzing the detailed structure of this face 
C  ...  the face may be interior and with boundary points
            Call makSP(iP1, IEP, IPE, IEE, MaxS, kE, iEu)
            Do 15 k = 1, kE
               iEt = iEu(k)
               Do n = 1, lE
                  If(iEt.EQ.iabs(iEs(n))) Goto 15
               End do
               Do j1 = 1, 4
                  j2 = iref(j1 + 1)
                  j3 = iref(j2 + 1)
                  
                  jP1 = IPE(j1, iEt)
                  jP2 = IPE(j2, iEt)
                  jP3 = IPE(j3, iEt)
                  If(check33(iP1, iP2, iP3, jP1, jP2, jP3)) Goto 30
               End do
 15         Continue

            IEE(i1, iE) = 0
         End if

         Do 20 k = 1, lE
            kE = iEs(k)
            If(kE.LE.0)  Goto 20
            If(.NOT.check3j(iP1, iP2, iP3, IPEs(1, k))) Goto 20
            If(kE.EQ.iE) Goto 20

            Do j1 = 1, 4
               j2 = iref(j1 + 1)
               j3 = iref(j2 + 1)

               jP1 = IPEs(j1, k)
               jP2 = IPEs(j2, k)
               jP3 = IPEs(j3, k)

               If(check33(iP1, iP2, iP3, jP1, jP2, jP3)) Then
                  IEE(i1, iE) = kE
                  IEE(j1, kE) = iE

                  IFE(j1, kE) = IFE(i1, iE)
                  Goto 30
               End if
            End do
 20      Continue
 30   Continue
      Return
      End Subroutine eleUpd


C ==========================================================
      Subroutine eleDel(iE, IPE, IEE)
C ==========================================================
      Integer IPE(5, *), IEE(4, *)

      IPE(1, iE) = 0

      Do i = 1, 4
         iEt = IEE(i, iE)
         If(iEt.NE.0) Then
            Do j = 1, 4
               If(IEE(j, iEt).EQ.iE) IEE(j, iEt) = 0
            End do
            IEE(i, iE) = 0
         End if
      End do
      Return
      End Subroutine eleDel
C
      End Module mba3d_update
