      Module mba3d_swapR
C
      use mba3d_auxSE
      use mba3d_auxSF
      use mba3d_auxSP
      use mba3d_error
      use mba3d_list
      use mba3d_makM
      use mba3d_makQ
      use mba3d_update
C
      contains
C
C ================================================================
      Subroutine swapR(
C ================================================================
c group (M)
     &           iwR, iwE,
     &           nF, MaxF, nE, MaxE,
     &           XYP, IPF, IPE,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           L1E, L2E, nL2, nStep,
     &           IHolF, IHolE,
     &           status,
C group (Q)
     &           HesP, rQuality, detG, qE,
C group (S)
     &           lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &           flag)
C ================================================================
      include 'makS.fd'
      include 'color.fd'
      include 'operat.fd'
      include 'shutF.fd'
      include 'magic.fd'
      include 'status.fd'
C ================================================================
C group (M)
      Integer IPF(4, *), IPE(5, *)
      Real*8  XYP(3, *)

      Real*8  hStar

      Integer L1E(2, *), L2E(*), nStep(4)
      Integer IHolF(*), IHolE(*)

      Integer ICP(*), IEP(*)
      Integer IFE(4, *), IEE(4, *)

      Integer status

C group (Q)
      Real*8  HesP(6, *), rQuality
      Real*8  detG(*), qE(*)

C group (S)
      Integer iFu(*), iEu(*), IPFu(4, *), IPEu(5, *)
      Real*8  qEu(*)

C group (Flag)
      Logical flag, flagTM

C ================================================================
C group (Functions)

C group (Local variables)
      Integer iFs(MaxS), iEs(MaxS)
      Integer IPFs(4, MaxS), IPEs(5, MaxS), ICPs(MaxS)
      Real*8  qEs(MaxS)

      Integer ip(5), iref(4), iPR(2, 6), iPt(3)
      Integer iDPs(2, MaxS), iNPs(2, MaxS)
      Integer iRs(3, MaxS), iPs(MaxS)
      Integer iLev(MaxS), iDel(MaxS)

      Real*8  XYPs(3), XYPd(3), v1, v2, c1
      Real*8  oldVolume, newVolume, vEs(MaxS)

      Logical flagBNDs, flagLOOP, flagSHUT, flagFACE, flagEDGE
      Logical l1, l2

C ================================================================
      DATA iPR /1,2, 1,3, 1,4, 2,3, 2,4, 3,4/
C ================================================================

      flag = .FALSE.

      Call copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1

      iref(1) = 1
      iref(2) = 2
      iref(3) = 3
      iref(4) = 1

      iPa = IPE(iPR(1, iwR), iwE)
      iPb = IPE(iPR(2, iwR), iwE)


C ... checking the case when swapping is impossible
      If(ifXnode(ICP(iPa), jTnode) .AND. 
     &   ifXnode(ICP(iPb), jTnode)) Goto 1000

      Call clrSR(iPa, iPb, ICP, IPF, IFE, lF, iFs, lE, iEs, ICRab)
      If(ifXnode(ICRab, jRnode)) Goto 1000



C ... checking for surrounding points
      Call makSR(iPa, iPb, lE, iEs, IPEs, lR, iRs, flagEDGE)
      If(.NOT.flagEDGE) Goto 1000

      If(ifXnode(status, ANIUse2ArmRule) .OR.
     &   ifXnode(status, ANIForbidBoundaryElements)) Then
         flagBNDs = .TRUE.
         Do n = 1, lR
            If(ifXnode(ICP(iRs(2, n)), jInode)) flagBNDs = .FALSE.
         End do
      Else
         flagBNDs = .FALSE.
      End if

      If(flagBNDs) Then
         iDPs(1, 1) = iPa
         iDPs(2, 1) = iPb

         If(ifXnode(ICP(iPa), jBnode) .AND. 
     &      ifXnode(ICP(iPb), jInode)) Then
            Call chkSPf(1, iPa, iSWAP, ICP, IEP, IPE, IEE, lP, iPs)
            Call chkSPb(2, 1, iDPs, 0, iNPs, iSWAP,
     &                  ICP, IEP, IPE, IEE, lP, iPs, iPbad, flagBNDs)
            If(flagBNDs) Goto 1000

         Else If(ifXnode(ICP(iPb), jBnode) .AND. 
     &           ifXnode(ICP(iPa), jInode)) Then
            Call chkSPf(1, iPb, iSWAP, ICP, IEP, IPE, IEE, lP, iPs)
            Call chkSPb(2, 1, iDPs, 0, iNPs, iSWAP,
     &                  ICP, IEP, IPE, IEE, lP, iPs, iPbad, flagBNDs)
            If(flagBNDs) Goto 1000

         Else If(ifXnode(ICP(iPb), jBnode) .AND. 
     &           ifXnode(ICP(iPa), jBnode)) Then
            lP = 2
            iPs(1) = iPa
            iPs(2) = iPb
            Call chkSPb(2, 1, iDPs, 0, iNPs, iSWAP,
     &                  ICP, IEP, IPE, IEE, lP, iPs, iPbad, flagBNDs)
            If(flagBNDs) Goto 1000
         End if
      End if


      flagLOOP = .TRUE.
      If(iRs(3, lR).NE.iRs(2, 1)) flagLOOP = .FALSE.


      lPend = lR
      Do n = 1, lR
         iPs(n) = iRs(2, n)
      End do
      If(.NOT.flagLOOP) Then
         lPend = lR + 1
         iPs(lPend) = iRs(3, lR)
      End if


C ... search for boundary points adjacent to points A and B
      flagFACE = .FALSE.
      icnt = 0
      Do n = 1, lF
         iPt(1) = IPFs(1, n)
         iPt(2) = IPFs(2, n)
         iPt(3) = IPFs(3, n)

         l1 = check13(iPa, iPt(1), iPt(2), iPt(3))
         l2 = check13(iPb, iPt(1), iPt(2), iPt(3))

         If(l1 .AND. l2) Then
            Do i = 1, 3
               If(iPt(i).NE.iPa .AND. iPt(i).NE.iPb) Then
                  flagFACE = .TRUE.
                  icnt = icnt + 1

                  If(icnt.EQ.1) iPc = iPt(i)
                  If(icnt.EQ.2) iPd = iPt(i)

                  iCLRf = IPFs(4, n)
                  iFs(n) = -iFs(n)
               End if
            End do
         End if
      End do


C ... only two faces are allowed to have the same edge
      If(flagFACE .AND. icnt.NE.2) Goto 1000

C ... checking piece-wise plane material boundaries
      If(flagFACE) Then
         c1 = angle2Faces(XYP(1, iPa), XYP(1, iPb),
     &                    XYP(1, iPc), XYP(1, iPd))
         If(c1.GT.angPREC ** 2 - 1D0) Goto 1000
      End if


C ... checking the number of inverted elements
      flagTM = ifXnode(status, ANITangledMesh)
      if(flagTM) Then
         nBad = 0
         Do n = 1, lR
            If(qEs(iRs(1, n)).LE.0D0) nBad = nBad + 1
         End do
         If(nBad.GT.0) Goto 1000
      End if


C ... coloring the superelement by element colors
C ... array iLev is overloaded
      Do n = 1, lR + 1
         ICPs(n) = 0
      End do

      Do n = 1, lR
         nEt = iRs(1, n)
         ide = IPEs(5, nEt)

         If(ICPs(n).NE.ide) ICPs(n) = ICPs(n) + ide
         If(ICPs(n + 1).NE.ide) ICPs(n + 1) = ICPs(n + 1) + ide
      End do

      If(flagLOOP .AND. ICPs(1).NE.ICPs(lR + 1)) 
     &   ICPs(1) = ICPs(1) + ICPs(lR + 1)


C ... checking that we have at most 2 different colors
      Do n = 1, lPend
         iLev(n) = ICPs(n)
      End do

      ic = countColors(lPend, iLev) 
      If(ic.GT.3) Goto 1000


C ... making a virtual evaluation of the quality
      lEold = lE
      oldVolume = 0D0
      Do n = 1, lR
         nEt = iRs(1, n)

         iP1 = IPEs(1, nEt)
         iP2 = IPEs(2, nEt)
         iP3 = IPEs(3, nEt)
         iP4 = IPEs(4, nEt)

         oldVolume = oldVolume +
     &               dabs(calVol(XYP(1, iP1), XYP(1, iP2),
     &                           XYP(1, iP3), XYP(1, iP4)))
      End do


C ... starting the main loop
      idx = 1
      nLev = lPend - 2
      iLev(1) = 0

      nItrSwap = 0

 100  iLev(idx) = iLev(idx) + 1
      Do 400 n1 = iLev(idx), lPend
         nItrSwap = nItrSwap + 1
         If(nItrSwap.GT.MaxItrSwap) Goto 1000

         If(iPs(n1).LE.0) Goto 400

         iP1 = iPs(n1)
         Do m = n1 + 1, lPend
            If(iPs(m).GT.0) Then
              m1 = m
              Goto 200
            End if
         End do
         Goto 400

 200     iP2 = iPs(m1)
         Do k = m1 + 1, lPend
            If(iPs(k).GT.0) Then
              k1 = k
              Goto 300
            End if
         End do
         Goto 400

 300     iP3 = iPs(k1)


         If(idx.LT.nLev) Then
            flagSHUT = lCloseTriag
            Call shutF(XYP(1, iP1), XYP(1, iP3),
     &                 XYP(1, iPa), XYP(1, iPb), XYP(1, iP2),
     &                 flagSHUT)
            If(.NOT.flagSHUT) Goto 400
         Else
            If(m1.NE.n1 + 1) Then
               i = iP2
               iP2 = iP3
               iP3 = i
            End if

            Do i = 1, 3
               XYPs(i) = (XYP(i, iP3) + XYP(i, iP1)) / 2
               XYPd(i) = (XYP(i, iPa) + XYP(i, iPb)) / 2
            End do

            v1 = calVol(XYP(1, iPa), XYP(1, iP1), XYP(1, iP2), XYPs)
            v2 = calVol(XYP(1, iPa), XYP(1, iP1), XYP(1, iP2), XYPd)

            If(v1 * v2.LE.0D0) Goto 400

            v1 = calVol(XYP(1, iPb), XYP(1, iP1), XYP(1, iP2), XYPs)
            v2 = calVol(XYP(1, iPb), XYP(1, iP1), XYP(1, iP2), XYPd)

            If(v1 * v2.LE.0D0) Goto 400
         End if


         nE1 = iRs(1, m1)
         Call calQE(
     &        HesP(1, iPa), detG(iPa), XYP(1, iPa),
     &        HesP(1, iP1), detG(iP1), XYP(1, iP1),
     &        HesP(1, iP2), detG(iP2), XYP(1, iP2),
     &        HesP(1, iP3), detG(iP3), XYP(1, iP3),
     &        hStar, qEs(nE1), vEs(nE1))

         If(qEs(nE1).LE.rQuality) Goto 400

         lE = lE + 1
         If(lE.GT.MaxS) Goto 9000
         iEs(lE) = 0

         Call calQE(
     &        HesP(1, iPb), detG(iPb), XYP(1, iPb),
     &        HesP(1, iP1), detG(iP1), XYP(1, iP1),
     &        HesP(1, iP2), detG(iP2), XYP(1, iP2),
     &        HesP(1, iP3), detG(iP3), XYP(1, iP3),
     &        hStar, qEs(lE), vEs(lE))

         If(qEs(lE).LE.rQuality) Then
            lE = lE - 1
            Goto 400
         End if

c  ...  go to the next level
         ide = min(ICPs(n1), ICPs(m1))
         ide = min(ide,      ICPs(k1))
c        ide = IPEs(5, nE1)  keep the mistaken choice

         IPEs(1, nE1) = iPa
         IPEs(2, nE1) = iP1
         IPEs(3, nE1) = iP2
         IPEs(4, nE1) = iP3
         IPEs(5, nE1) = ide

         IPEs(1, lE)  = iPb
         IPEs(2, lE)  = iP1
         IPEs(3, lE)  = iP2
         IPEs(4, lE)  = iP3
         IPEs(5, lE)  = ide

         iNPs(1, idx) = iPs(n1)
         iNPs(2, idx) = iPs(k1)

         iDel(idx) = m1
         iPs(m1) = -iPs(m1)

         If(idx.EQ.nLev) Then
            newVolume = 0D0
            Do n = lEold + 1, lE
               newVolume = newVolume + dabs(vEs(n))
            End do

            Do n = 1, nLev
               nEt = iRs(1, iDel(n))
               newVolume = newVolume + dabs(vEs(nEt))
            End do

C ... checking only in the case of plane surfaces
            If(dabs(oldVolume - newVolume).GT.volPREC*oldVolume) Then
C              Write(*, 5000) oldVolume, newVolume
               Goto 400
            End if

            Goto 500
         Else If(idx.EQ.nLev - 1 .AND. flagFACE) Then
            Do n = 1, nLev - 1
               Do j = 1, 2
                  If(iNPs(j, n).EQ.iPc .AND.
     &               iNPs(3 - j, n).EQ.iPd) Then
                     Goto 350
                  End if
               End do
            End do

            Goto 400
         End if

 350     idx = idx + 1
         iLev(idx) = 0
         Goto 100
 400  Continue


      idx = idx - 1
      If(idx.EQ.0) Goto 1000

      m1 = iDel(idx)
      iPs(m1) = -iPs(m1)

      lE = lE - 1
      Goto 100



C ... analyzing curvilinear and plane faces
 500  lFold = lF

      If(flagFACE) Then
         lF = lF + 2
         If(lF.GT.MaxS) Goto 9000

         Do m = 1, lFold
            If(check33(iPa, iPc, iPd, IPFs(1, m),
     &                 IPFs(2, m), IPFs(3, m))) Goto 1000
         End do

         IPFs(1, lF - 1) = iPa
         IPFs(2, lF - 1) = iPc
         IPFs(3, lF - 1) = iPd
         IPFs(4, lF - 1) = iCLRf

         Do m = 1, lFold
            If(check33(iPb, iPc, iPd, IPFs(1, m),
     &                 IPFs(2, m), IPFs(3, m))) Goto 1000
         End do

         IPFs(1, lF) = iPb
         IPFs(2, lF) = iPc
         IPFs(3, lF) = iPd
         IPFs(4, lF) = iCLRf

         Call facAdd(iFs(lF - 1), nF, MaxF, IHolF)
         Call facAdd(iFs(lF), nF, MaxF, IHolF)
      End if


c ... checking for boundary elements
      If(ifXnode(status, ANIForbidBoundaryElements)) Then 
         Do 700 n = 2, lPend - 1
            nEt = iRs(1, n)
            Do i = 1, 4
               iPu = IPEs(i, nEt)
               If(ifXnode(ICP(iPu), jInode)) Goto 700
            End do
            Goto 2000
 700     Continue

         Do 800 n = lEold + 1, lE
            Do i = 1, 4
               iPu = IPEs(i, n)
               If(ifXnode(ICP(iPu), jInode)) Goto 800
            End do
            Goto 2000
 800     Continue
      End if


C ... updating the grid
      flag = .TRUE.

      If(flagFACE) Then
         Do n = 1, lFold
            iFt = iFs(n)
            If(iFt.LE.0) Then
               iFt = -iFt
               Call facDel(iFt, nF, IPF, IHolF)
            End if
         End do


         Call facUpd(lF - 1, IPF, iFs, IPFs)
         Call facUpd(lF, IPF, iFs, IPFs)
      End if


      nEt = iRs(1, 1)
      iEt = iEs(nEt)
      iEs(nEt) = -iEt
      Call lstDel(nE, L1E, nL2, L2E, nStep, IHolE, qE, iEt)
      Call eleDel(iEt, IPE, IEE)

      If(flagLOOP) Then
         nEt = iRs(1, lR)
         iEt = iEs(nEt)
         iEs(nEt) = -iEt
         Call lstDel(nE, L1E, nL2, L2E, nStep, IHolE, qE, iEt)
         Call eleDel(iEt, IPE, IEE)
      End if


      Do n = 2, lPend - 1
         nEt = iRs(1, n)
         iEt = iEs(nEt)
         Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEt, qEs(nEt))
         Call eleDel(iEt, IPE, IEE)
      End do

      Do n = lEold + 1, lE
         Call eleAdd(nE, MaxE, IHolE)
         Call lstAdd(nE, L1E, nL2, L2E, nStep, IHolE,
     &               qE, qEs(n), iEs(n))
         Call eleDel(iEs(n), IPE, IEE)
      End do


      Do n = 2, lPend - 1
         nEt = iRs(1, n)
         Call eleUpd(nEt, IEP, IPE, IFE,  IEE,
     &               lF,  lE,  iFs, iEs, IPFs, IPEs)
      End do

      Do n = lEold + 1, lE
         Call eleUpd(n,  IEP, IPE, IFE,  IEE,
     &               lF, lE,  iFs, iEs, IPFs, IPEs)
      End do


 5000 Format('Warning in swapR: bad volumes =', 2E16.8)

 1000 Return

c ... delete wrong faces
 2000 Do n = lFold + 1, lF
         Call facDel(iFs(n), nF, IPF, IHolF)
      End do
      Return

 9000 Call errMes(1007, 'swapR', 'local parameter MaxS is small')
      End Subroutine swapR
C
      End Module mba3d_swapR
