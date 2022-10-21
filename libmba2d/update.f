C ==========================================================
C PONITS ... POINT ... POINTS ... POINT
C ==========================================================
      Subroutine pntAdd(iP, nP, MaxP, ICP,  XYP,  HesP,  detG, IHolP,
     &                                ICPs, XYPs, HesPs, detGs)
      Integer ICP(*), IHolP(*)
      real  XYP(2, *),  HesP(3, *),  detG(*)
      real  XYPs(2),    HesPs(3),    detGs

      nHolP = IHolP(1)
      nP = nP + 1
      If(nP.GT.MaxP) goto 1000

      If(nHolP.EQ.0) Then
         iP = nP
      Else
         iP = IHolP(nHolP + 1)
         nHolP = nHolP - 1
         IHolP(1) = nHolP
      End if

      Call pntUpd(iP, ICP, XYP, HesP, detG, ICPs, XYPs, HesPs, detGs)

      Return

 1000 Continue
      Call errMes(1003, 'pntAdd', 
     &           'local variable MaxP is too small')
      End


C ==========================================================
      Subroutine pntUpd(iP, ICP,  XYP,  HesP,  detG,
     &                      ICPs, XYPs, HesPs, detGs)
C ==========================================================
      Integer ICP(*)
      real  XYPs(2), XYP(2, *)
      real  HesPs(3), HesP(3, *)
      real  detGs, detG(*)

      ICP(iP) = ICPs

      Do i = 1, 2
         XYP(i, iP) = XYPs(i)
      End do

      Do i = 1, 3
         HesP(i, iP) = HesPs(i)
      End do

      detG(iP) = detGs
      Return
      End



C ==========================================================
      Subroutine pntDel(iP, nP, ICP, IHolP)
C ==========================================================
      Integer ICP(*), IHolP(*)

      nHolP = IHolP(1)
      nHolP = nHolP + 1
      IHolP(nHolP + 1) = iP

      ICP(iP) = 0

      IHolP(1) = nHolP
      nP = nP - 1
      Return
      End



C ==========================================================
C FACES ... FACE ... FACES ... FACE
C ==========================================================
      Subroutine facAdd(iF, nF, MaxF, IHolF)
      Integer IHolF(*)

      nHolF = IHolF(1)
      nF = nF + 1
      If(nF.GT.MaxF) goto 1000

      If(nHolF.EQ.0) Then
         iF = nF
      Else
         iF = IHolF(nHolF + 1)
         nHolF = nHolF - 1
         IHolF(1) = nHolF
      End if

      Return

 1000 Continue
      Call errMes(1004, 'facAdd', 
     &           'local variable MaxF is too small')
      End


C ==========================================================
      Subroutine facUpd(nFs, IPF, parCrv, iFnc,
     &                  iFs, IPFs, iCRVs, iFNCs, iBNDs, t1, t2)
C ==========================================================
      Integer IPF(4, *), iFnc(*)
      Integer iFs(*), IPFs(2, *)
      real  parCrv(2, *), t1, t2

      iF = iFs(nFs)
      Do i = 1, 2
         IPF(i, iF) = iPFs(i, nFs)
      End do

      IPF(3, iF) = iCRVs
      iFnc(iF) = iFNCs

      If(iCRVs.NE.0) Then
         parCrv(1, iF) = t1
         parCrv(2, iF) = t2
      Else
         parCrv(1, iF) = 0D0
         parCrv(2, iF) = 0D0
      End if

      IPF(4, iF) = iBNDs
      Return
      End



C ==========================================================
      Subroutine facDel(iF, nF, IPF, iFnc, IHolF)
C ==========================================================
      Integer IPF(4, *), iFnc(*), IHolF(*)

      nHolF = IHolF(1)
      nHolF = nHolF + 1
      IHolF(nHolF + 1) = iF

      Do i = 1, 4
         IPF(i, iF) = 0
      End do
      iFnc(iF) = 0 

      IHolF(1) = nHolF
      nF = nF - 1
      Return
      End



C ==========================================================
C ELEMENTS ... ELEMENT ... ELEMENTS ... ELEMENT
C ==========================================================
      Subroutine eleAdd(nE, MaxE, IHolE)
      Integer IHolE(*)

      If(nE.GE.MaxE .AND. IHolE(1).EQ.0) goto 1000

      Return

 1000 Continue
      Call errMes(1006, 'eleAdd', 
     &           'local variable MaxE is too small')
      End



C ==========================================================
      Subroutine eleUpd(nEs, IEP, IPE, IFE, IEE,
     &                  lF, lE, iFs, iEs, IPFs, IPEs)
C ==========================================================
      include 'makS.fd'
C ==========================================================
      Integer IEP(*), IPE(3, *), IFE(3, *), IEE(3, *)
      Integer IPFs(2, *), IPEs(3, *)
      Integer iFs(*), iEs(*)

C group (Local variables)
      Integer iEu(MaxS)
      Integer iref(4)

      Logical check22

C ==========================================================
      kE = 0
      nMaxS = MaxS
      iref(1) = 1
      iref(2) = 2
      iref(3) = 3
      iref(4) = 1

      iE = iEs(nEs)

      Do i = 1, 3
         IPE(i, iE) = IPEs(i, nEs)
         IEP(IPE(i, iE)) = iE
         IFE(i, iE) = 0
      End do

      Do 30 i1 = 1, 3
         i2 = iref(i1 + 1)

         iP1 = IPEs(i1, nEs)
         iP2 = IPEs(i2, nEs)

         Do n = 1, lF
            If(iFs(n).GT.0) Then
               jP1 = IPFs(1, n)
               jP2 = IPFs(2, n)

               If(iP1.EQ.jP1 .AND. iP2.EQ.jP2 .OR.
     &            iP1.EQ.jP2 .AND. iP2.EQ.jP1) Then
                  IFE(i1, iE) = iFs(n)
                  goto 10
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
                  If(iEt.EQ.iabs(iEs(n))) goto 15
               End do
               Do j1 = 1, 3
                  j2 = iref(j1 + 1)
                  
                  jP1 = IPE(j1, iEt)
                  jP2 = IPE(j2, iEt)
                  If(check22(iP1, iP2, jP1, jP2)) goto 30
               End do
 15         Continue

            IEE(i1, iE) = 0
         End if
         
         Do 20 k = 1, lE
            kE = iEs(k)
            If(kE.LE.0)  goto 20
            If(kE.EQ.iE) goto 20

            Do j1 = 1, 3
               j2 = iref(j1 + 1)

               jP1 = IPEs(j1, k)
               jP2 = IPEs(j2, k)

               If(iP1.EQ.jP1 .AND. iP2.EQ.jP2 .OR.
     &            iP1.EQ.jP2 .AND. iP2.EQ.jP1) Then
                  IEE(i1, iE) = kE
                  IEE(j1, kE) = iE
                  goto 30
               End if
            End do
 20      Continue
 30   Continue

      Return
      End



C ==========================================================
      Subroutine eleDel(iE, IPE, IEE)
C ==========================================================
      Integer IPE(3, *), IEE(3, *)

      IPE(1, iE) = 0

      Do i = 1, 3
         iEt = IEE(i, iE)
         If(iEt.NE.0) Then 
            Do j = 1, 3   
               If(IEE(j, iEt).EQ.iE) IEE(j, iEt) = 0
            End do
            IEE(i, iE) = 0
         End if 
      End do
      Return
      End

