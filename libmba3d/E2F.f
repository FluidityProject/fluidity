      Module mba3d_E2F
C
      use mba3d_auxSE
      use mba3d_auxSF
      use mba3d_auxSP
      use mba3d_auxSR
      use mba3d_error
      use mba3d_list
      use mba3d_makQ
      use mba3d_update
C
      contains
C
C ================================================================
      Subroutine E2F(
C ================================================================
c group (M)
     &           iwR, iwE,
     &           nF, MaxF, nE, 
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
      Logical flag

C ================================================================
C group (Functions)

C group (Local variables)
      Integer iFs(MaxS), iEs(MaxS)
      Integer IPFs(4, MaxS), IPEs(5, MaxS)
      Real*8  qEs(MaxS)

      Integer iRs(3, MaxS), iPR(2, 6)
      Integer iDs(2, 1), iNs(2, 1), iPs(MaxS)
      Integer iPt(3)
      Real*8  v1, v2, c1

      Logical flagBNDs, flagFACE, flagEDGE
      Logical flagTM, l1, l2

C ================================================================
      DATA iPR /1,2, 1,3, 1,4, 2,3, 2,4, 3,4/
C ================================================================
      flag = .FALSE.

C ... checking operations which do not modify the initial superelement
      iPa = IPE(iPR(1, iwR), iwE)
      iPb = IPE(iPR(2, iwR), iwE)
      Call makSR(iPa, iPb, lEu, iEu, IPEu, lR, iRs, flagEDGE)
      If(.NOT.flagEDGE) Goto 1000


C ... checking the case when edge -> face operation is impossible
      If(lR.NE.3) Goto 1000

      iP1 = iRs(2, 1)
      iP2 = iRs(3, 1)
      iP3 = iRs(3, 2)
      iP4 = iRs(3, 3)

      If(iP1.NE.iP4) Goto 1000

      Call clrSR(iPa, iPb, ICP, IPF, IFE, lFu, iFu, lEu, iEu, ICRab)
      If(ifXnode(ICRab, jRnode)) Goto 1000


C ... checking the orientation
      v1 = calVol(XYP(1, iP1), XYP(1, iP2),
     &            XYP(1, iP3), XYP(1, iPa))
      v2 = calVol(XYP(1, iP1), XYP(1, iP2),
     &            XYP(1, iP3), XYP(1, iPb))
      If(v1 * v2.GE.0D0) Goto 1000


C ... saving the initial superelement structure (it is changed in the sequel)
      Call copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)


C ... number of inverted elements
      flagTM = ifXnode(status, ANITangledMesh)
      if(flagTM) Then
         nBad = 0
         If(qEs(iRs(1, 1)).LE.0D0) nBad = nBad + 1
         If(qEs(iRs(1, 2)).LE.0D0) nBad = nBad + 1
         If(qEs(iRs(1, 3)).LE.0D0) nBad = nBad + 1

         If(nBad.GT.0) Goto 1000
      End if


C ... making a virtual evaluation of the quality
      nE1 = iRs(1, 1)
      nE2 = iRs(1, 2)
      nE3 = iRs(1, 3)
      ide = IPEs(5, nE1)

      Call calQE(
     &     HesP(1, iPa), detG(iPa), XYP(1, iPa),
     &     HesP(1, iP1), detG(iP1), XYP(1, iP1),
     &     HesP(1, iP2), detG(iP2), XYP(1, iP2),
     &     HesP(1, iP3), detG(iP3), XYP(1, iP3),
     &     hStar, qEs(nE1), v1)

      If(qEs(nE1).LE.rQuality) Goto 1000

      Call calQE(
     &     HesP(1, iPb), detG(iPb), XYP(1, iPb),
     &     HesP(1, iP1), detG(iP1), XYP(1, iP1),
     &     HesP(1, iP2), detG(iP2), XYP(1, iP2),
     &     HesP(1, iP3), detG(iP3), XYP(1, iP3),
     &     hStar, qEs(nE2), v2)

      If(qEs(nE2).LE.rQuality) Goto 1000

C ... checking for surrounding points
      If(ifXnode(ICP(iP1), jBnode) .AND.
     &   ifXnode(ICP(iP2), jBnode) .AND.
     &   ifXnode(ICP(iP3), jBnode)) Then

         If(ifXnode(status, ANIUse2ArmRule)) Then
           iDs(1, 1) = iPa
           iDs(2, 1) = iPb

           If(ifXnode(ICP(iPa), jBnode) .AND. 
     &        ifXnode(ICP(iPb), jInode)) Then
              Call chkSPf(1, iPa, iE2F, ICP, IEP, IPE, IEE, lP, iPs)
              Call chkSPb(2, 1, iDs, 0, iNs, iE2F,
     &                    ICP, IEP, IPE, IEE, lP, iPs, iPbad, flagBNDs)
              If(flagBNDs) Goto 1000

           Else If(ifXnode(ICP(iPb), jBnode) .AND. 
     &             ifXnode(ICP(iPa), jInode)) Then
              Call chkSPf(1, iPb, iE2F, ICP, IEP, IPE, IEE, lP, iPs)
              Call chkSPb(2, 1, iDs, 0, iNs, iE2F,
     &                    ICP, IEP, IPE, IEE, lP, iPs, iPbad, flagBNDs)
              If(flagBNDs) Goto 1000

           Else If(ifXnode(ICP(iPb), jBnode) .AND. 
     &             ifXnode(ICP(iPa), jBnode)) Then
              lP = 2
              iPs(1) = iPa
              iPs(2) = iPb
              Call chkSPb(2, 1, iDs, 0, iNs, iE2F,
     &                    ICP, IEP, IPE, IEE, lP, iPs, iPbad, flagBNDs)
              If(flagBNDs) Goto 1000
           End if
         End if
      End if



C ... search for boundary points adjacent to point A and B
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


C ... analyzing curvilinear and plane faces
      lFold = lF

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
         If(ifXnode(ICP(iP1), jInode)) Goto 100
         If(ifXnode(ICP(iP2), jInode)) Goto 100
         If(ifXnode(ICP(iP3), jInode)) Goto 100

         If(ifXnode(ICP(iPa), jBnode)) Goto 1000
         If(ifXnode(ICP(iPb), jBnode)) Goto 1000
      End if


C ... updating the grid
 100  flag = .TRUE.

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


C ... creating new elements
      IPEs(1, nE1) = iPa
      IPEs(2, nE1) = iP1
      IPEs(3, nE1) = iP2
      IPEs(4, nE1) = iP3
      IPEs(5, nE1) = ide

      IPEs(1, nE2) = iPb
      IPEs(2, nE2) = iP1
      IPEs(3, nE2) = iP2
      IPEs(4, nE2) = iP3
      IPEs(5, nE2) = ide


      iE3 = iEs(nE3)
      iEs(nE3) = -iE3

      Call lstDel(nE,  L1E, nL2, L2E, nStep, IHolE, qE, iE3)
      Call eleDel(iE3, IPE, IEE)

      Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(nE1), qEs(nE1))
      Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(nE2), qEs(nE2))

      Call eleDel(iEs(nE1), IPE, IEE)
      Call eleDel(iEs(nE2), IPE, IEE)

      Call eleUpd(nE1, IEP, IPE, IFE,  IEE,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs)
      Call eleUpd(nE2, IEP, IPE, IFE,  IEE,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs)

 1000 Return
 9000 Call errMes(1007, 'E2F', 'local parameter MaxS is small')
      End Subroutine E2F
C
      End Module mba3d_E2F
