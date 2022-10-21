      Module mba3d_ani2
C
      use mba3d_clpsR
      use mba3d_deletP
      use mba3d_dsort
      use mba3d_E2F
      use mba3d_F2E
      use mba3d_auxSE
      use mba3d_error
      use mba3d_insrtP
      use mba3d_lintrp3D
      use mba3d_makM
      use mba3d_moveP
      use mba3d_output
      use mba3d_splitE
      use mba3d_swapR
C
      contains
C
C ================================================================
      Subroutine ani2(
C ================================================================
c group (M)
     &      nP, MaxP, nF, MaxF, nE, MaxE, 
     &      XYP, IPF, IPE, 
     &      nEStar, hStar,
     &      ICP, IPP, IEP, IFE, IEE,
     &      L1E, L2E,
     &      IHolP, IHolF, IHolE,
     &      IEPw, nEPw,
     &      miLINTRP, mrLINTRP, IPEw, iSE,
     &      flagAuto, status,
c group (Dev)
     &      nPv, nFv, nEv, IPV, IFV, IEV, 
c group (Q)
     &      MaxSkipE, MaxQItr, MaxBaskets,
     &      HesP, Quality, rQuality, 
     &      detG, qE, nPw, nEw, XYPw, HesPw, rSE,
     &      MetricFunction, flagAnalytic,
c group (ERR)
     &      flagFILE, chanelOUT, nLines, output,
     &      iPrint, nQItr, iERR)
C ================================================================
      include 'makS.fd'
      include 'color.fd'
      include 'magic.fd'
      include 'status.fd'
      include 'output.fd'
      include 'operat.fd'
C ================================================================
C The core of the package where the algorithm structure is realized
C
C We mention here only essential parameters that differ from ones
C in MeshMetric and MeshSolution.
C
C     IPE(5, MaxE) - columns 1,2,3,4 are the connectivity
C                            list of elements;
C                    column  5 is the element identificator
C
C     IPF(4, MaxF) - columns 1,2,3 are the connectivity
C                            list of faces;
C                    column  4 is the face identificator
C
C     IPEw(4, nE)  - container for keeping original
C                    connectivity list for elements
C ================================================================
C group (M)
C     Integer MaxP, MaxF, MaxE, nPv, nEStar
      Real*8  XYP(3, *)
      Integer IPE(5, *), IPF(4, *), IPV(*)

      Real*8  hStar

      Integer ICP(*), IPP(*), IEP(*)
      Integer IFE(4, *), IEE(4, *)
      Integer L1E(2, *), L2E(*), nStep(4)
      Integer IHolP(*), IHolF(*), IHolE(*)

      Integer IEPw(*), nEPw(*)
      Integer IPEw(4, *), iSE(*)

      Logical flagAuto
      Integer status

C group (Dev)
      Integer IFV(*), IEV(*) 

C group (Q)
C     Integer MaxSkipE, MaxQItr, MaxBaskets
      Real*8  HesP(6, *)
      Real*8  Quality, rQuality
      Real*8  detG(*), qE(*)
      Real*8  XYPw(3, *), HesPw(6, *), rSE(*)
      Logical flagAnalytic

      Integer  MetricFunction
      EXTERNAL MetricFunction

C group (output)
      Logical flagFILE
      Integer chanelOUT, nLines
      Character*(LineLenght) output(*)

C ================================================================
C group (Local functions)

C group (Local variables)
      Integer iwF(4), iwR(6)
      Integer iFu(MaxS), iEu(MaxS)
      Integer IPFu(4, MaxS), IPEu(5, MaxS)
      Real*8  qEu(MaxS)

      Real*8  XYPs(3), HesPs(6), rMove
      Logical flagTest

      Real*8  tm1, tm2

      Real*8  ver, dvo, sao
      Real*8  aer, fvo, fao  
      Real*8  s, sa, rR
      Real*8  aQuality

      Logical flagUAR, flagFBE, flagTM
      Character*(LineLenght) message

C ================================================================
      Real*8  aniXY0(3), aniXYA(3), aniXYB(3)
      Common /aniPlane/aniXY0, aniXYA, aniXYB

C ================================================================
      Call mba3d_seconds(tm1)

      nF2E = 0
      nE2F = 0
      nSWAP  = 0
      nINSRT = 0
      nDELET = 0
      nMOVE  = 0
      nCLPS  = 0

      mF2E = 0
      mE2F = 0
      mSWAP  = 0
      mINSRT = 0
      mDELET = 0
      mMOVE  = 0
      mCLPS  = 0

      flagTM = .TRUE.


c ... Loop initialization
      nQItrAdd = 0
      nQItrBig = 0
      nQItr = 0
      iERR = 0

      nFo = nF
      Call makM(
c group (M)
     &      nP, nF, MaxF, nE, MaxE, 
     &      XYP, IPF, IPE, 
     &      ICP, IPP, IEP, IFE, IEE,
     &      IHolP, IHolF, IHolE,
     &      IEPw, nEPw,
     &      status,
c group (Dev) 
     &      nPv, nFv, nEv, IPV, IFV, IEV,
c group (iERR)
     &      iERR)
      If(iERR.NE.0) Goto 9000

      If(nFo.NE.nF) Then
         If(.NOT.flagAuto) Call errMes(4001, 'ani2', 
     &                                'inconsistent input data')
         If(iPrint.GE.1) Then
            If(flagFILE) Then
               Write(message, 5009) nF - nFo
               Call addOut(nLines, message, output)
            Else
               Write(*, 5009) nF - nFo
            End if
         End if
      End if


c ... gather statistics about elements. We assume that IHolE(1) = 0
      If(iPrint.GE.3) Then
         mE = 0
         Do n = 1, nE
            If(IPE(1, n).GT.0) Then
               mE = mE + 1
               IEPw(mE) = IPE(5, n)
            End if
         End do

         ic = countColors(mE, IEPw)

         Do i = 1, ic
            s = domainVolume(nE, XYP, IPE, IEPw(i))

            If(flagFILE) Then
               Write(message, 5002) i, ic, IEPw(i), s
               Call addOut(nLines, message, output)
            Else
               Write(*, 5002) i, ic, IEPw(i), s
            End if
         End do

         s  = fixedVolume(nEv, XYP, IPE, IEV)
         sa = fixedArea(  nFv, XYP, IPF, IFV)

         If(flagFILE) Then
            Write(message, 5007) nEv, s
            Call addOut(nLines, message, output)
            Write(message, 5008) nEv, sa
            Call addOut(nLines, message, output)
         Else
            Write(*, 5007) nEv, s
            Write(*, 5008) nEv, sa
         End if
      End if


c ... gather statistics about faces. We allow IHolF(1) >= 0
      If(iPrint.GE.4) Then
         mF = 0
         Do n = 1, nF + IHolF(1)
            If(IPF(4, n).GT.0) Then
               mF = mF + 1
               IEPw(mF) = IPF(4, n)
            End if
         End do

         ic = countColors(mF, IEPw)

         Do i = 1, ic
            s = surfaceArea(nF + IHolF(1), XYP, IPF, IEPw(i))

            If(flagFILE) Then
               Write(message, 5003) i, ic, IEPw(i), s
               Call addOut(nLines, message, output)
            Else
               Write(*, 5003) i, ic, IEPw(i), s
            End if
         End do
      End if


c ... check the mesh (level 1)
      Call chkM(
c group (M)
     &          nP, nF, nE,
     &          XYP, IPF, IPE,
     &          ICP, IFE, IEE,
     &          rR, status,
c group (W)
     &          1, nEPw, IEPw)


c ... check the mesh (level 2)
      fvo = fixedVolume(nEv, XYP, IPE, IEV)
      fao = fixedArea(  nFv, XYP, IPF, IFV)


c ... output the statistics 
      If(iPrint.GE.2) Then
         If(flagFILE) Then
            Write(message, 5004) rR, status
            Call addOut(nLines, message, output)
         Else
            Write(*, 5004) rR, status
         End if
      End if


c ... modify the grid to satisfy some restrictios (usually with FEs)
      flagUAR = ifXnode(status, ANIUse2ArmRule)
      flagFBE = ifXnode(status, ANIForbidBoundaryElements)

      If(flagUAR .OR. flagFBE) Then
         Call mba3d_seconds(tm2)
         If(iPrint.GE.2) Then
            rQuality = 1D0 
            Do n = 1, nE
               rQuality = min(rQuality, qE(n))
            End do

            If(flagFILE) Then
               Write(message, 5000) nQItr, rQuality,
     &                              nP, nF, nE, tm2 - tm1
               Call addOut(nLines, message, output)
            Else
               Write(*, 5000) nQItr, rQuality, nP, nF, nE, tm2 - tm1
            End if
         End if
         
         kE = nE
         Do 20 n = 1, kE
            Do i = 1, 4
               iPt = IPE(i, n)
               If(ifXnode(ICP(iPt), jInode)) Goto 20
            End do

            nQItr = nQItr + 1

            Call makSE(n, ICP, IEP, IPF, IPE, IFE, IEE, qE, MaxS,
     &                 lF, lE, iFu, iEu, IPFu, IPEu, qEu)

            Call splitE(
C group (M)
     &           n,
     &           nP, MaxP, nE, MaxE,
     &           XYP, IPE,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           IHolP, IHolE,
C group (Q)
     &           HesP, detG, qE,
C group (S)
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
     &           flagTest)
 20      Continue
      End if


c ... initialize the list
      If(ifXnode(status, ANITangledMesh)) Then
         Do n = 1, nE
            Call updQa(n, XYP, IPE, IEE, qE)
         End do
      End if

      Do n = 1, nE
         rSE(n) = qE(n) 
         L2E(n) = n
      End do
      Call DSORT(rSE, L2E, nE, 2, iERR)
      If(iERR.NE.0) Goto 9000


      nStep(1) = max(100.0, log(real(nEStar)))
      nStep(2) = nStep(1) / 4
      nStep(3) = MaxE
      nStep(4) = 0

      Call lstMak(nEStar, nE, L1E, L2E, nL2, nStep, IHolE)


c ... initialize the 8-tree
      LDH = 6
      nXY = 0
      iControl = 1100 + chanelOUT 

      If( .NOT.flagAnalytic ) Then
         Call LINTRP3D(nEw, IPEw, nPw, XYPw, LDH, HesPw, nXY, XYPs,
     &                 HesPs, iSE, miLINTRP, rSE, mrLINTRP, iControl)
      End if

c ... switch-off initialization
      iControl = iControl + 1000 


c ... output of initial qualities
      If(iPrint.GE.2) 
     &   Call countBadElements(nE, L1E, L2E, qE, Quality, 
     &                         nLines, output, flagFILE)


c ... empty the basket and repeat the main loop
 100  nQItrBig = nQItrBig + 1
      If(nQItrBig.GT.MaxBaskets) Then
         iERR = 1000
         Goto 1000
      End if


      Call mba3d_seconds(tm2)
      If(iPrint.GE.3 .OR.
     &   iPrint.GE.1 .AND. nQItrBig.EQ.1) Then
         aQuality = avgQ(nE, qE, L1E, L2E)

         If(flagFILE) Then
            Write(message, 5010) aQuality, rR, status
            Call addOut(nLines, message, output)

            Write(message, 5000) nQItr, qE(L2E(1)), nP, nF, nE, tm2-tm1
            Call addOut(nLines, message, output)
         Else
            Write(*, 5010) aQuality, rR, status
            Write(*, 5000) nQItr, qE(L2E(1)), nP, nF, nE, tm2-tm1
         End if
      End if
      If(MaxQItr.LE.0) Goto 1100


      nSkipE = 0
 300  nQItr = nQItr + 1
      If(nQItr.GT.MaxQItr + nQItrAdd) Then
         If(nQItrBig.EQ.1) Then
            nQItrAdd = nSkipE
            Goto 100
         End if

         iERR = 1000
         Goto 1000
      End if


c ... check the mesh (level 3)
      If(qE(L2E(1)).GT.0D0 .AND. flagTM) Then
         flagTM = .FALSE.

         Call delXnode(status, ANITangledMesh)

         dvo = domainVolume(nE + IHolE(1), XYP, IPE, 0)
         sao = surfaceArea( nF + IHolF(1), XYP, IPF, 0)
      End if


c ... check for termination signals (dummy function for simulators)   
      If((nQItr / 200) * 200.EQ.nQItr  .AND.  probeAny()) Then
         If(nQItrBig.EQ.1) Then
            nQItrAdd = nSkipE
            Goto 100
         End if

         iERR = 1000
         Goto 1000
      End if


      iwE = L2E(1)
      Do i = 1, nSkipE
         iwE = L1E(2, iwE)
      End do

      If((nQItr / 1000) * 1000.EQ.nQItr) Then
         Call mba3d_seconds(tm2)
         If(iPrint.GE.2) Then
            If(flagFILE) Then
               Write(message, 5000) nQItr, qE(L2E(1)),
     &                              nP, nF, nE, tm2 - tm1
               Call addOut(nLines, message, output)
            Else
               Write(*, 5000) nQItr, qE(L2E(1)),
     &                        nP, nF, nE, tm2 - tm1
            End if
         End if

         m = 1D1 * (tm2 - tm1)
         If(iPrint.GE.2 .AND. (m / 19) * 19.EQ.m .AND. m.GT.0
     &                  .AND. .NOT.flagFILE) Call info
      End if
      If(iPrint.GE.3 .AND. (nQItr / 10000) * 10000.EQ.nQItr) Then
         If(flagFILE) Then
            Write(message,6001) nSkipE, nQItrBig, iwE
            Call addOut(nLines, message, output)
            Write(message,6002)
            Call addOut(nLines, message, output)
            Write(message,6003) 
     &            nF2E, nE2F, nSWAP, nINSRT, nDELET, nMOVE, nCLPS
            Call addOut(nLines, message, output)
            Write(message,6003) 
     &            mF2E, mE2F, mSWAP, mINSRT, mDELET, mMOVE, mCLPS
            Call addOut(nLines, message, output)
         Else
            Write(*, 5001) 
     &            nSkipE, nQItrBig, iwE,
     &            nF2E, nE2F, nSWAP, nINSRT, nDELET, nMOVE, nCLPS,
     &            mF2E, mE2F, mSWAP, mINSRT, mDELET, mMOVE, mCLPS
         End if
      End if

      rQuality = qE(iwE)
      If(rQuality.GT.Quality) Then
         If(nSkipE.EQ.0) Goto 1000

         Goto 100
      End if

      iP1 = IPE(1, iwE)
      iP2 = IPE(2, iwE)
      iP3 = IPE(3, iwE)
      iP4 = IPE(4, iwE)


      Call calQF(HesP(1, iP1), XYP(1, iP1),
     &           HesP(1, iP2), XYP(1, iP2),
     &           HesP(1, iP3), XYP(1, iP3),
     &           HesP(1, iP4), XYP(1, iP4),
     &           hStar, iwF, iwR)

      Call makSE(iwE, ICP, IEP, IPF, IPE, IFE, IEE, qE, MaxS,
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu)


      Do iop = 1, NUMop
C ... analyze edges of the element
        If(iINSRTop.EQ.iop) Then
          Do i = 1, 6
            Call insrtP(
c group (M)
     &           iwR(i), iwE,
     &           nP, MaxP, nF, MaxF, nE, MaxE,
     &           XYP, IPF, IPE,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           L1E, L2E, nL2, nStep,
     &           IHolP, IHolF, IHolE,
     &           status,
C group (Q)
     &           HesP, rQuality, detG, qE,
     &           MetricFunction, flagAnalytic,
C group (S)
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &           nPw, nEw, XYPw, HesPw, IPEw,
     &           miLINTRP, mrLINTRP, iSE, rSE, iControl,
     &           flagTest)
            mINSRT = mINSRT + 1
            If(flagTest) Then
               nINSRT = nINSRT + 1
               Goto 400
            End if
          End do

C ... analyze faces of the element
        Else If(iF2Eop.EQ.iop) Then
          Do i = 1, 4
            Call F2E(
c group (M)
     &           iwF(i), iwE,
     &           nE, MaxE, XYP, IPE,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           L1E, L2E, nL2, nStep,
     &           IHolE,
     &           status,
C group (Q)
     &           HesP, rQuality, detG, qE,
C group (S)
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
c group (W)
     &           flagTest)
            mF2E = mF2E + 1
            If(flagTest) Then
               nF2E = nF2E + 1
               Goto 400
            End if
          End do

C ... analyze edges of the element
        Else If(iE2Fop.EQ.iop) Then
          Do i = 1, 6
            Call E2F(
c group (M)
     &           iwR(i), iwE,
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
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &           flagTest)
            mE2F = mE2F + 1
            If(flagTest) Then
               nE2F = nE2F + 1
               Goto 400
            End if
          End do

C ... analyze edges of the element
        Else If(iCLPSRop.EQ.iop) Then
          Do i = 1, 6
            Call clpsR(
c group (M)
     &           i, iwE,
     &           nP, nE, XYP, IPE,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           L1E, L2E, nL2, nStep,
     &           IHolP, IHolE,
     &           status,
C group (Q)
     &           HesP, rQuality, detG, qE,
     &           MetricFunction, flagAnalytic,
C group (S)
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &           nPw, nEw, XYPw, HesPw, IPEw,
     &           miLINTRP, mrLINTRP, iSE, rSE, iControl,
     &           flagTest)

            mCLPS = mCLPS + 1
            If(flagTest) Then
               nCLPS = nCLPS + 1
               Goto 400
            End if
          End do

c ... swap edge (general method)
        Else If(iSWAPop.EQ.iop) Then
          Do i = 1, 6
            Call swapR(
c group (M)
     &           iwR(i), iwE,
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
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &           flagTest)
            mSWAP = mSWAP + 1
            If(flagTest) Then
               nSWAP = nSWAP + 1
               Goto 400
            End if
          End do 

C ... analyze points of the element
        Else If(iDELETop.EQ.iop) Then
          Do i = 1, 4
            Call deletP(
c group (M)
     &           i, iwE,
     &           nP, nF, MaxF, nE, MaxE,
     &           XYP, IPF, IPE,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           L1E, L2E, nL2, nStep,
     &           IHolP, IHolF, IHolE,
     &           status,
C group (Q)
     &           HesP, rQuality, detG, qE,
C group (S)
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &           flagTest)
            mDELET = mDELET + 1
            If(flagTest) Then
               nDELET = nDELET + 1
               Goto 400
            End if
          End do

c ... moving point
        Else If(iMOVEop.EQ.iop) Then
          Do i = 1, 4
            Call moveP(
c group (M)
     &           i, iwE,
     &           nE, XYP, IPF, IPE,
     &           hStar,
     &           ICP, IFE, 
     &           L1E, L2E, nL2, nStep,
     &           status,
C group (Q)
     &           HesP, rQuality, detG, qE,
     &           MetricFunction, flagAnalytic,
C group (S)
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &           nPw, nEw, XYPw, HesPw, IPEw,
     &           miLINTRP, mrLINTRP, iSE, rSE, iControl, rMove,
     &           flagTest)
            mMOVE = mMOVE + 1
            If(flagTest) Then
               nMOVE = nMOVE + 1
               Goto 400
            End if
          End do
        End if
      End do


C ... does an operation increase the quality sufficiently enough?
 400  itE = L2E(1)
      Do i = 1, nSkipE
         itE = L1E(2, itE)
      End do
      If(itE.EQ.iwE) flagTest = .FALSE.


      If(flagTest) Goto 300


      nSkipE = nSkipE + 1
      If(nSkipE.GT.MaxSkipE .OR. nSkipE.GE.nE) Goto 100
      Goto 300

 1000 rQuality = qE(L2E(1))


      Call mba3d_seconds(tm2)
      If(iPrint.GE.1) Then
         If(flagFILE) Then
            Write(message, 5000) nQItr, qE(L2E(1)),
     &                           nP, nF, nE, tm2 - tm1
            Call addOut(nLines, message, output)
         Else
            Write(*, 5000) nQItr, qE(L2E(1)), nP, nF, nE, tm2 - tm1
         End if
      End if


C ... calculate the number of bad tetrahedrons
 1100 Continue

      If(iPrint.GE.3) Then
         If(flagFILE) Then
            Write(message,6001) nSkipE, nQItrBig, iwE
            Call addOut(nLines, message, output)
            Write(message,6002)
            Call addOut(nLines, message, output)
            Write(message,6003) 
     &            nF2E, nE2F, nSWAP, nINSRT, nDELET, nMOVE, nCLPS
            Call addOut(nLines, message, output)
            Write(message,6003) 
     &            mF2E, mE2F, mSWAP, mINSRT, mDELET, mMOVE, mCLPS
            Call addOut(nLines, message, output)
         Else
            Write(*, 5001) 
     &            nSkipE, nQItrBig, iwE,
     &            nF2E, nE2F, nSWAP, nINSRT, nDELET, nMOVE, nCLPS,
     &            mF2E, mE2F, mSWAP, mINSRT, mDELET, mMOVE, mCLPS
         End if
      End if

      If(iPrint.GE.2) 
     &   Call countBadElements(nE, L1E, L2E, qE, Quality, 
     &                         nLines, output, flagFILE)



C ... update the mesh
      nFold = nF
      Call updM(
c group (M)
     &            nP, nF, nE, 
     &            XYP, IPF, IPE, 
     &            ICP, IPP, IFE, IEE,
     &            IHolP, IHolF, IHolE,
     &            status,
c group (Dev)
     &            nPv, nFv, nEv, IPV, IFV, IEV,
c group (Q)
     &            HesP, qE,
c group (W)
     &            IEPw)


c ... check the mesh (level 1)
      Call chkM(
c group (M)
     &          nP, nF, nE,
     &          XYP, IPF, IPE,
     &          ICP, IFE, IEE,
     &          rR, status,
c group (W)
     &          0, nEPw, IEPw)


c ... print out details
      If(iPrint.GE.1) Then
         averageQ = avgQ(nE, qE, L1E, L2E)
         If(flagFILE) Then
            Write(message, 5010) averageQ, rR, status
            Call addOut(nLines, message, output)
         Else
            Write(*, 5010) averageQ, rR, status
         End if
      End if


      If(iPrint.GE.3) Then
         mE = 0
         Do n = 1, nE
            If(IPE(1, n).GT.0) Then
               mE = mE + 1
               IEPw(mE) = IPE(5, n)
            End if
         End do
                                      
         ic = countColors(mE, IEPw)
                                               
         Do i = 1, ic
            s = domainVolume(nE, XYP, IPE, IEPw(i))
                                    
            If(flagFILE) Then
               Write(message, 5002) i, ic, IEPw(i), s
               Call addOut(nLines, message, output)
            Else
               Write(*, 5002) i, ic, IEPw(i), s
            End if
         End do
      End if


c ... check the mesh (level 2)
      If(nF.GT.0) Then 
         If(iPrint.GE.3) Then
            If(flagFILE) Then
               Write(message, 5006) dvo, sao
               Call addOut(nLines, message, output)
            Else
               Write(*, 5006) dvo, sao
            End if
         End if
         ver = dabs(dvo - domainVolume(nE, XYP, IPE, 0)) / dvo
         aer = dabs(sao - surfaceArea( nF, XYP, IPF, 0)) / sao
         If(ver.GT.volPREC .OR. 
     &      aer.GT.volPREC .AND. nFold.EQ.nF) Call wrnMes(5201, 
     &        'ani2', 'total volume or/and surface is not preserved')

         ver = dabs(fvo - fixedVolume(nEv, XYP, IPE, IEV)) / dvo
         aer = dabs(fao - fixedArea(  nFv, XYP, IPF, IFV)) / sao
         If(ver.GT.volPREC .OR. aer.GT.volPREC) Call errMes(5201, 
     &        'ani2', 'fixed volume or/and surface is not preserved')

      End if

      If(iPrint.GE.2) Then
         If(flagFILE) Then
            Write(message, 5004) rR, status
            Call addOut(nLines, message, output)
         Else
            Write(*, 5004) rR, status
         End if
      End if


c ... display output
 5000 Format('ITRs:', I7, ' Q=', E10.4, '  #P,F,E:', I7,I8,I9,
     &       '  tm=', F6.1, 's')

 5001 Format(/, I4, '  elements in the ', I3, '-th basket,',
     &          I7, '  is the current bad element', /,
     &'Face-Edge Edge-Face  Gen.Swap  Insert  Delete    Move  Collapse'
     & ,/, 2(I9,I10,I10,I8,I8,I8,I10/))

 5002 Format( 'Domain:', I3, '/', I2, ',  clr=', I4, ',   vol=', E12.6)
 5003 Format('Surface:', I3, '/', I2, ',  clr=', I4, ',  area=', E12.6)

 5004 Format('Maximal R/r =', E10.3,
     &       '  (R/r = 3 for equilateral tetrahedron),  status.fd:', I4)

 5006 Format('Domain volume and total boundary area: ', 2E14.6)

 5007 Format( /,'Fixed domain:', I6, '  tets,   vol=', E12.6)
 5008 Format(   'Fixed  area: ', I6, ' faces,  area=', E12.6)

 5009 Format('Warning:', I6, ' new faces have been added') 

 5010 Format('Avg Quality is', E11.4, ',  Maximal R/r =', E11.4,
     &       ',  status.fd:', I5)


c ... parallel output
 6001 Format(I4, '  elements in the ', I3, '-th basket,',
     &       I7, '  is the current bad element')
 6002 Format(
     &'Face-Edge  Edge-Face  Gen.Swap  Insert  Delete   Move  Collapse')
 6003 Format(I10,I12,I10,I8,I8,I7,I9)

 9000 Return
      End Subroutine ani2
C
      End Module mba3d_ani2
