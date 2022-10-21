C =========================================================
      Subroutine infoP(iP1, iPa, iPb, iF1, iF2, par,
     &                 IPF, parCrv, lf, iFs)
C =========================================================
      Integer IPF(4, *), iFs(*)
      real  parCrv(2, *), par(*)
      Logical flag1

C =========================================================
      flag1 = .FALSE.
      Do 10 n = 1, lF
         iFc = iFs(n)
         If(iFc.LE.0) goto 10

         iPl = IPF(1, iFc)
         iPr = IPF(2, iFc)

         If(iPl.EQ.iP1 .AND. flag1) Then
            i4  = 2
            iF2 = iFc
            iPb = iPr
         End if

         If(iPr.EQ.iP1 .AND. flag1) Then
            i4  = 1
            iF2 = iFc
            iPb = iPl
         End if

         If(iPl.EQ.iP1 .AND. .NOT.flag1) Then
            i1  = 2
            iF1 = iFc
            iPa = iPr
            flag1 = .TRUE.
         End if

         If(iPr.EQ.iP1 .AND. .NOT.flag1) Then
            i1  = 1
            iF1 = iFc
            iPa = iPl
            flag1 = .TRUE.
         End if
 10   Continue
      i2 = 3 - i1
      i3 = 3 - i4

      par(1) = parCrv(i1, iF1)
      par(2) = parCrv(i2, iF1)
      par(3) = parCrv(i3, iF2)
      par(4) = parCrv(i4, iF2)

      Return
      End


C =========================================================
      Subroutine infoF(iF, iP1, iP2, iF1, iF2, iPc, iPd,
     &            par, IPF, parCrv, lF, iFs)
C =========================================================
      implicit none
      Integer IPF(4, *), iFs(*)
      real  parCrv(2, *), par(*)
      integer :: iF, iP1, iP2, iF1, iF2, iPc, iPd, lF

      integer :: i1, i2, i3, i4, i5, i6
      integer :: iFc, iBnd, iPl, iPr, n

      i1 = 0; i2 = 0; i3 = 0; i4 = 0; i5 = 0; i6 = 0

C =========================================================
      iPc = 0
      iPd = 0

      iF1 = 0
      iF2 = 0

      iBnd = IPF(4, iF)

      Do 10 n = 1, lF
         iFc = iFs(n)
         If(iFc.LE.0) goto 10
         If(IPF(4, iFc).NE.iBnd) goto 10

         iPl = IPF(1, iFc)
         iPr = IPF(2, iFc)

         If(iPl.EQ.iP1 .AND. iPr.NE.iP2) Then
            i1  = 2
            iF1 = iFc
            iPc = iPr
         End if

         If(iPr.EQ.iP1 .AND. iPl.NE.iP2) Then
            i1  = 1
            iF1 = iFc
            iPc = iPl
         End if

         If(iPl.EQ.iP1 .AND. iPr.EQ.iP2) i3 = 1
         If(iPr.EQ.iP1 .AND. iPl.EQ.iP2) i3 = 2

         If(iPl.EQ.iP2 .AND. iPr.NE.iP1) Then
            i6  = 2
            iF2 = iFc
            iPd = iPr
         End if

         If(iPr.EQ.iP2 .AND. iPl.NE.iP1) Then
            i6  = 1
            iF2 = iFc
            iPd = iPl
         End if
 10   Continue
      i2 = 3 - i1
      i4 = 3 - i3
      i5 = 3 - i6


      If(iF1.NE.0) Then 
         par(1) = parCrv(i1, iF1)
         par(2) = parCrv(i2, iF1)
      End if

      par(3) = parCrv(i3, iF)
      par(4) = parCrv(i4, iF)

      If(iF2.NE.0) Then
         par(5) = parCrv(i5, iF2)
         par(6) = parCrv(i6, iF2)
      End if

      Return
      End


C =========================================================
      Subroutine makSE(iE, IEP, IPF, IPE, IFE, IEE, qE, MaxS,
     &                 lF, lE, iFs, iEs, IPFs, IPEs, qEs,
     &                 status)
C =========================================================
C Remark: The 1st column of the IPE is temporary overloaded.
C         The negative values of IPE(1, *) are used for a 
C         quick search in the list iEs.
C =========================================================
      include 'status.fd'
C =========================================================
c group (M)
      Integer IEP(*), IPF(4, *)
      Integer IPE(3, *), IFE(3, *), IEE(3, *)
      real  qE(*)

c group (S)
      Integer iFs(*), iEs(*), IPFs(2, *), IPEs(3, *)
      real  qEs(*)

c group (Control)
      Integer status

c group (Local variables)
      Logical ifXnode

C =========================================================
      kE = 0
      nMaxS = MaxS

      Do i = 1, 3
         iP = IPE(i, iE)
         Call makSP(iP, IEP, IPE, IEE, nMaxS, kEadd, iEs(kE + 1))

         kE = kE + kEadd
         nMaxS = nMaxS - kEadd
      End do

      lE = 0
      Do 10 n = 1, kE
         iEt = iEs(n)
         If(IPE(1, iEt).LT.0) goto 10 

         lE = lE + 1
         iEs(lE) = iEt

         Do i = 1, 3
            IPEs(i, lE) = IPE(i, iEt)
         End do

         IPE(1, iEt) = -IPE(1, iEt)
 10   Continue


c ... adding elements adjacent by an edge (for local topological operations)
      If(ifXnode(status, ANISmoothMesh) .OR.
     &   ifXnode(status, ANIUntangleMesh)) Then
         kE = lE
         Do n = 1, kE
            iEt = iEs(n)
            If(lE.GT.MaxS - 3) goto 1000

            Do 40 i = 1, 3
               jEt = IEE(i, iEt)
               If(jEt.LE.0) goto 40
               If(IPE(1, jEt).LT.0) goto 40 

               lE = lE + 1
               iEs(lE) = jEt

               Do k = 1, 3
                  IPEs(k, lE) = IPE(k, jEt)
               End do

               IPE(1, jEt) = -IPE(1, jEt)
  40        Continue
         End do
      End if


c ... restoring the 1st colunm of IPE
      Do n = 1, lE
         IPE(1, iEs(n)) = -IPE(1, iEs(n))
      End do


c ... adding surface edges
      lF = 0
      Do n = 1, lE
         kE = iEs(n)
         qEs(n) = qE(kE)
         If(IPEs(1, n).EQ.0) Call errMes(6001, 
     &                            'makSE', 'System error') 

         Do 20 i = 1, 3
            iF = IFE(i, kE)
            If(iF.NE.0) Then
               Do m = 1, lF 
                  If(iF.EQ.iFs(m)) goto 20
               End do

               lF = lF + 1
               iFs(lF) = iF

               Do j = 1, 2
                  IPFs(j, lF) = IPF(j, iF)
               End do
            End if
 20      Continue
      End do

      Return
 1000 Continue
      Call errMes(1007, 'makSE', 
     &           'local variable MaxS is too small')
      End


C =========================================================
      Subroutine makSP(iP, IEP, IPE, IEE, MaxS, nS, iSE)
C =========================================================
C Remark: the 3rd column of IPE is overloaded to avoid 
C         seachinf th ine list iSE.
C =========================================================
c group (M)
      Integer IEP(*)
      Integer IPE(3, *), IEE(3, *)

C group (S)
      Integer iSE(*)

C group (Local variables)
      Logical repeat

C =========================================================
      nS = 1
      iSE(nS) = IEP(iP)

      iE = iSE(1)
      IPE(3, iE) = -IPE(3, iE)

 1    repeat = .FALSE.
      n1 = nS
      n2 = max(1, nS - 1)

      Do 4 n = n2, n1
         If(nS.GE.MaxS-3) goto 1000

         Do 2 i = 1, 3
            iE = IEE(i, iSE(n))
            If(iE.EQ.0) goto 2
            If(IPE(3, iE).LT.0) goto 2

            If(iP.EQ.IPE(1, iE) .OR.
     &         iP.EQ.IPE(2, iE) .OR.
     &         iP.EQ.IPE(3, iE)) Then
c              Do k = nS, 1, -1
c                 If(iE.EQ.iSE(k)) goto 2
c              End do

               repeat = .TRUE.
               nS = nS + 1

               iSE(nS) = iE
               IPE(3, iE) = -IPE(3, iE)
            End if
 2       Continue
 4    Continue

      If(repeat) goto 1

c ... restoring the overloaded values
      Do k = 1, nS 
         iE = iSE(k)
         IPE(3, iE) = -IPE(3, iE)
      End do

      Return
 1000 Continue
      Call errMes(1007, 'makSP', 
     &           'local variable MaxS is too small')
      End



C =========================================================
      Subroutine chkSPf(iP1, iP2, ICP, IEP, IPE, IEE, lPf, iPf)
C =========================================================
      include 'makS.fd'
      include 'colors.fd'
C =========================================================
c group (M)
      Integer ICP(*), IEP(*)
      Integer IPE(3, *), IEE(3, *)

C group (S)
      Integer iPf(*)

C group (Local variables)
      Integer iEs(MaxS)
      Logical ifXnode

C =========================================================
      lPf = 0
      Do nStep = 1, 2
         If(nStep.EQ.1) Then
            m2 = 2
            iPf(1) = iP1
            iPf(2) = iP2
         Else
            m2 = lPf
         End if

         Do m = 1, m2
            Call makSP(iPf(m), IEP, IPE, IEE, MaxS, lE, iEs)
            Do n = 1, lE
               iE = iEs(n)
               Do 10 i = 1, 3
                  iPt = IPE(i, iE)
                  If(ifXnode(ICP(iPt), jBnode)) Then
                     Call findSE(lPf, iPf, iPt, nPt)
                     If(nPt.NE.0) goto 10

                     lPf = lPf + 1
                     If(lPf.GT.MaxS) goto 1000
                     iPf(lPf) = iPt
                  End if
 10            Continue
            End do
         End do
      End do
      Return

 1000 Continue
      Call errMes(1007, 'chkSPf', 
     &           'local variable MaxS is too small')
      End


C =========================================================
      Subroutine chkSPb(iD1, iD2, iN1, iN2, iOPERAT,
     &                    ICP, IEP, IPE, IEE, lPf, iPf, flag)
C =========================================================
      include 'makS.fd'
      include 'colors.fd'
      include 'operat.fd'
C =========================================================
c group (M)
      Integer ICP(*), IEP(*)
      Integer IPE(3, *), IEE(3, *)

C group (S)
      Integer iPf(*)
      Logical flag

C group (Local variables)
      Integer iPb(MaxS), iEs(MaxS)
      Logical ifXnode

C =========================================================
      flag = .TRUE.
      Do 10 l = 1, lPf
         lPb = 1
         iPb(1) = iPf(l)

         Do nStep = 1, 2
            m2 = lPb

            Do m = 1, m2
               iPc = iPb(m)
               Call makSP(iPc, IEP, IPE, IEE, MaxS, lE, iEs)
               Do n = 1, lE
                  iE = iEs(n)
                  Do 5 i = 1, 3
                     iPt = IPE(i, iE)

                     If(iOPERAT.EQ.iCLPS) Then
                        If(iPt.EQ.iD2) goto 1
                     Else If(iOPERAT.EQ.iSWAP .OR.
     &                       iOPERAT.EQ.iINSRT) Then
                        If(iPc.EQ.iD1 .AND. iPt.EQ.iD2) goto 5
                        If(iPc.EQ.iD2 .AND. iPt.EQ.iD1) goto 5
                     End if

                     If(ifXnode(ICP(iPt), jInode)) goto 10

 1                   If(nStep.EQ.1) Then
                        If(iOPERAT.EQ.iCLPS) Then
                           If(iPt.EQ.iD1) Then
                              lPb = lPb + 2
                              If(lPb.GT.MaxS) goto 1000

                              iPb(lPb - 1) = iD1
                              iPb(lPb) = iD2
                           Else If(iPt.EQ.iD2) Then
                              lPb = lPb + 2
                              If(lPb.GT.MaxS) goto 1000
                              
                              iPb(lPb - 1) = iD2
                              iPb(lPb) = iD1
                           End if
                        Else If(iOPERAT.EQ.iSWAP) Then
                           If(iPt.EQ.iN1) Then
                              If(ifXnode(ICP(iN2), jInode)) goto 10
                           Else If(iPt.EQ.iN2) Then
                              If(ifXnode(ICP(iN1), jInode)) goto 10
                           End if
                        End if

                        Call findSE(lPb, iPb, iPt, nPt)
                        If(nPt.NE.0) goto 5

                        lPb = lPb + 1
                        If(lPb.GT.MaxS) goto 1000
                        iPb(lPb) = iPt
                     End if
 5                Continue
               End do
            End do
         End do
         Return
 10   Continue

      flag = .FALSE.
      Return

 1000 Continue
      Call errMes(1007, 'chkSPb', 
     &           'local variable MaxS is too small')
      End


C =========================================================
      Subroutine calSO(XYP, IPE, lE, iEs, iOs)
C =========================================================
      real  XYP(2, *)
      Integer IPE(3, *), iEs(*), iOs(*)

      real  d

C =========================================================
      Do 10 n = 1, lE
         iE = iEs(n)
         If(iE.LE.0) goto 10

         iP1 = IPE(1, iE)
         iP2 = IPE(2, iE)
         iP3 = IPE(3, iE)

         d = (XYP(1, iP2) - XYP(1, iP1)) *
     &       (XYP(2, iP3) - XYP(2, iP1)) -
     &       (XYP(1, iP3) - XYP(1, iP1)) *
     &       (XYP(2, iP2) - XYP(2, iP1))
         iOs(n) = sign(1.0, d)
 10   Continue
      Return
      End




C =========================================================
      Subroutine chkSO(iP1, iP2, XYPs, XYP, IPE, lE, iEs, iOs, flag)
C =========================================================
      real  XYPs(2), XYP(2, *)
      Integer IPE(3, *), iEs(*), iOs(*)
      Logical flag

      real  x(3), y(3), d

C =========================================================
      Do 10 n = 1, lE
         iE = iEs(n)
         If(iE.LE.0) goto 10

         Do i = 1, 3
            iPt = IPE(i, iE)
            If(iPt.EQ. iP1 .OR. iPt.EQ.iP2) Then
               x(i) = XYPs(1)
               y(i) = XYPs(2)
            Else
               x(i) = XYP(1, iPt)
               y(i) = XYP(2, iPt)
            End if
         End do

         d = (x(2) - x(1)) * (y(3) - y(1)) -
     &       (x(3) - x(1)) * (y(2) - y(1))

         If(d.EQ.0D0) goto 1000

         iOt = sign(1.0, d)
         If(iOt.NE.iOs(n)) Then
            flag = .FALSE.
            goto 1000
         End if
 10   Continue

      flag = .TRUE.
 1000 Return
      End



C =========================================================
C     Subroutine findSP(lP, iPs, iP, nPs)
C     Subroutine findSF(lF, iFs, iF, nFs)
      Subroutine findSE(lE, iEs, iE, nEs)
C =========================================================
C Search for a position of number iE in array iEs(lE).
C Zero position is returned if search had failed.
C ================================================================
      Integer iEs(*)

      nEs = 0
      Do n = 1, lE
         If(iEs(n).EQ.iE) Then
            nEs = n
            goto 1000
         End if
      End do
 1000 Return
      End


C =========================================================
      Subroutine copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &                  lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)
C =========================================================
      Integer iFu(*), iEu(*), IPFu(2, *), IPEu(3, *)
      Integer iFs(*), iEs(*), IPFs(2, *), IPEs(3, *)
      real  qEu(*), qEs(*)

C =========================================================
      lF = lFu
      Do n = 1, lF
         iFs(n) = iFu(n)

         Do i = 1, 2
            IPFs(i, n) = IPFu(i, n)
         End do
      End do

      lE = lEu
      Do n = 1, lE
         iEs(n) = iEu(n)
         qEs(n) = qEu(n)

         Do i = 1, 3
            IPEs(i, n) = IPEu(i, n)
         End do
      End do
      Return
      End


C =========================================================
      Subroutine copySQ(lE, qEu, XYPu, HesPu, detGu,
     &                      qEs, XYPs, HesPs, detGs)
C =========================================================
      real  qEu(*), XYPu(2), HesPu(3), detGu
      real  qEs(*), XYPs(2), HesPs(3), detGs

C =========================================================
      Do n = 1, lE
         qEs(n) = qEu(n)
      End do

      Do i = 1, 2
         XYPs(i) = XYPu(i)
      End do

      Do i = 1, 3
         HesPs(i) = HesPu(i)
      End do

      detGs = detGu
      Return
      End


