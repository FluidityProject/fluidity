C ==========================================================
      Subroutine ani2(
C ==========================================================
c group (M)
     &      nP, MaxP, nF, MaxF, nE, MaxE, nPv,
     &      XYP, IPF, IPE, IPV,
     &      calCrv, ParCrv, iFnc,
     &      nEStar, hStar,
     &      ICP, IEP, IFE, IEE,
     &      L1E, L2E,
     &      IHolP, IHolF, IHolE,
     &      IEPw, nEPw,
     &      IPEw, iSE,
c group (Dev)
     &      nFv, nEv, IFV, IEV, lbE,
     &      flagAuto, status,
c group (CRV)
     &      L1Et, L2Et, tE,
     &      nL2t, nStept, nEt,
     &      LFnc, ILt,
c group (Q)
     &      MaxSkipE, MaxQItr, nQItr,
     &      HesP, Quality, rQuality,
     &      detG, qE, XYPw, HesPw, rSE,
     &      MetricFunction, flagAnalytic,
c group (ERR)
     &      iPrint, iERR)
C ==========================================================
      include 'makS.fd'
      include 'colors.fd'
      include 'status.fd'
      include 'operat.fd'
      include 'magic.fd'
C ==========================================================
C The main driver routine. 
C
C *** DATA FLOW CHART
C
C  -> put points in increasing order
C  -> dublicate initial metric
C  -> build various mesh structures (cross maps)
C  -> check initial mesh and self-check other mesh structures
C  -> add mesh points to satify 2-arm rule  
C
C  -> compute qualities of mesh elements
C  -> create a list of triangle qualities
C  -> loop of curvilinear boundaries
C  ->  -> created lists of curvilinear edges
C
C  -> create quad-tree
C  -> while( #loop <= #Baskets )
C  -----> begin infinite basket-loop
C  ---------> take an element (E) with the worst quality
C  ---------> order edges of element E
C  ---------> collect a superelement around E
C
C  ---------> loop over edges of E
C  -------------> collapse the edge
C  -------------> insert point in the middle of the edge
C  -------------> swap the edge
C
C  ---------> loop over vertices of E
C  -------------> swap an edge ending at the vertex of E
C  -------------> move the vertex of E
C  
C  ---------> add E to the basket is the operations have failed
C  ---------> terminate basket-loop if the basket is full
C  -----> clean the basket
C
C  -> calculate the distribution of bad elements
C  -> compress the mesh data structure (clean holes)
C  -> check the final mesh
C  -> check the area and the perimetr of the meshed area
C 
C ==========================================================
C group (M)
C     Integer MaxP, MaxF, MaxE, nPv, nEStar
      real  XYP(2, *)
      Integer IPE(3, *), IPF(4, *), IPV(*)

      EXTERNAL calCrv
      real  ParCrv(2, *)
      Integer iFnc(*)

      real  hStar

      Integer ICP(*), IEP(*)
      Integer IFE(3, *), IEE(3, *)
      Integer L1E(2, *), L2E(*), nStep(4)
      Integer IHolP(*), IHolF(*), IHolE(*)

      Integer IEPw(*), nEPw(*)
      Integer IPEw(3, *), iSE(*)

c group (Dev)
      Integer nFv, nEv
      Integer IFV(*), IEV(*), lbE(*)

      Logical flagAuto
      Integer status

C group (CRV)
      real  tE(*)
      Integer L1Et(2, *), L2Et(*)
      Integer nL2t(*), nStept(4, *), nEt(*)
      Integer LFnc(*), ILt(*)

C group (Q)
C     Integer MaxSkipE, MaxQItr
      real  HesP(3, *)
      real  Quality, rQuality
      real  detG(*), qE(*)
      real  XYPw(2, *), HesPw(3, *),  rSE(*)
      Logical flagAnalytic

      Integer  MetricFunction
      EXTERNAL MetricFunction

C group (Local variables)
      Integer iw(3)
      Integer iFu(MaxS), iEu(MaxS)
      Integer IPFu(2, MaxS), IPEu(3, MaxS)
      real  qEu(MaxS)

      real  XYPs(2), HesPs(3), rMove
      Logical flagTM, flagTest

      Integer nCLPS1, nCLPS2, nINSRT, nMOVE, nSWAP
      real  tm1, tm2

      real  domainArea, domainPerimetr, dao, dpo, rR
      real  dao_new, dpo_new
      real  avgQ, aQuality

      Logical ifXnode

      Integer iDummy(1)
      Real    ANItime, tmdata(2)

C ==========================================================
      Integer iDomBnd, iMatBnd
      Common /aniBND/ iDomBnd, iMatBnd
C ==========================================================
      tm1 = ANItime(tmdata)

      nCLPS1 = 0
      nCLPS2 = 0
      nINSRT = 0
      nSWAP  = 0
      nMOVE  = 0
      nNOTHING = 0

      mCLPS1 = 0
      mCLPS2 = 0
      mINSRT = 0
      mSWAP  = 0
      mMOVE  = 0

      flagTM = .TRUE.

      nQItr = 0

c ... put initial data in increasing order
c     Do n = 1, nE
c        Call orderijk(IPE(1, n), IPE(2, n), IPE(3, n))
c     End do


c ... duplicating initial data
      nPw = nP
      Do n = 1, nPw
         Do i = 1, 3
            HesPw(i, n) = HesP(i, n)
         End do

         Do i = 1, 2
            XYPw(i, n) = XYP(i, n)
         End do
      End do

      nEw = nE
      Do n = 1, nEw
         Do i = 1, 3
            IPEw(i, n) = IPE(i, n)
         End do
      End do


c ... building mesh data structure
      nPVo = nPv
      nPo = nP
      nFo = nF
      Call makM(
c group (M)
     &      nP, MaxP, nF, MaxF, nE, MaxE, nPv,
     &      XYP, IPF, IPE, IPV,
     &      parCrv, iFnc,
     &      ICP, IEP, IFE, IEE,
     &      IHolP, IHolF, IHolE,
     &      IEPw, nEPw,
c group (Dev)
     &      nFv, nEv, IFV, IEV, lbE,
     &      status,
c group (iERR)
     &      iERR)
      If(iERR.NE.0) goto 9000


      If(nPVo.NE.nPV) Then
         If(.NOT.flagAuto) Call errMes(4001, 'ani2', 
     &                                'inconsistent input data')
         If(iPrint.GE.1) Write(*, 5007) nPV - nPVo
      End if

      If(nFo.NE.nF) Then
         If(.NOT.flagAuto) Call errMes(4001, 'ani2', 
     &                                'inconsistent input data')
         If(iPrint.GE.1) Write(*, 5008) nF - nFo
      End if


c ... checking (and self-checking) of the initial mesh
      Call chkM(nP, nF, nE,
     &            XYP, IPF, IPE, IFE, IEE, lbE,
     &          calCrv, parCrv, iFnc,
     &            ICP, rR, status)


c ... remove boundary triangles
      If(ifXnode(status, ANIForbidBoundaryElements)) Then
         kE = nE
         Do 20 n = 1, kE
            Do i = 1, 3
               iPt = IPE(i, n)
               If(ifXnode(ICP(iPt), jInode)) goto 20
            End do

            nQItr = nQItr + 1

            Call makSE(n, IEP, IPF, IPE, IFE, IEE, qE, MaxS,
     &                 lF, lE, iFu, iEu, IPFu, IPEu, qEu,
     &                 status)

            Call splitE(
c group (M)
     &           n,
     &           nP, MaxP, nF, MaxF, nE, MaxE,
     &           XYP, IPF, IPE, lbE,
     &           parCrv, iFnc,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           L1E, L2E, nL2, nStep,
     &           IHolP, IHolF, IHolE,
     &           status,
C group (Q)
     &           HesP, Quality, rQuality,
     &           detG, qE,
C group (S)
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &           flagTest)

!            If(.NOT.flagTest) Call errMes(4103, 'ani2.f', 
!     &                       'The input data contradicts status')
 20      Continue
      End if


c ... compute quality of elements
      Call makQ(
c group (M)
     &     nP, nE, XYP, IPE, IEE, nEV, IEV,
     &     nEStar, hStar, status,
c group (Q)
     &     HesP, detG, qE)


c ... initilize the list
      If(ifXnode(status, ANIUntangleMesh)) Then
         Do n = 1, nE
            Call updQa(n, XYP, IPE, IEE, qE)
         End do
      End if

      Do n = 1, nE
         rSE(n) = qE(n)
         L2E(n) = n
      End do
      Call DSORT(rSE, L2E, nE, 2, iERR)
      If(iERR.NE.0) goto 9000

c nStep(1) - typical interval length 
c nStep(2) - rank of interval length is [nStep(1)-nStep(2),nStep(1)+nStep(2)]
c nStep(3) - ipos2 of the second part of L2E
c nStep(4) - output channel in case of debugging (=0 for no debugging)

c     nStep(1) = sqrt(real(nEStar))
c     nStep(2) = 0 
      nStep(1) = log(real(nEStar)) 
      nStep(2) = nStep(1) / 4 
      nStep(3) = MaxE
      nStep(4) = 0
      Call lstMak(nEStar, nE, L1E, L2E, nL2, nStep, IHolE)

      If(iPrint.GE.3) Then
         Call statistics(nE, L1E, L2E, Quality, qE) 
      End if


c ... CURVILINEAR FACES
      Call calCrvFnc(IPF, nF, iFnc, LFnc, nCrvFnc)

      ir = 1
      Do n = 1, nCrvFnc
         ILt(n) = ir
         Call tEMak(tE(ir), nEt(n), MaxF, IPF, nF, parCrv, iFnc,
     &              LFnc(n))

         If(nEt(n).GT.0) Then
            Do i = 1, nEt(n)
               rSE(i) = tE(ir + i - 1)
               L2Et(ir + i - 1) = i
            End do
            Call DSORT(rSE, L2Et(ir), nEt(n), 2, iERR)
            If(iERR.NE.0) goto 9000

            nStept(1, n) = sqrt(real(nEt(n)))
            nStept(2, n) = 0
            nStept(3, n) = MaxF
            nStept(4, n) = 0 
            Call lstMak(nEt(n), nEt(n), L1Et(1, ir), L2Et(ir),
     &                  nL2t(n), nStept(1, n), iDummy)
         End if
         ir = ir + nEt(n)
      End do


c ... initialize the 4-tree
      If( .NOT.flagAnalytic ) Then
         Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 0, XYPs,
     &               HesPs,  iSE, rSE, .TRUE.)
      End if


      tm2 = ANItime(tmdata)
      If(iPrint.GE.1) Then
         aQuality = avgQ(nE, qE, L1E, L2E)

         Write(*, 5006) aQuality, rR, status
         Write(*, 5000) nQItr, qE(L2E(1)), nP, nF, nE, tm2 - tm1
      End if


C ... main loop
      nQItrAdd = 0
      nQItrBig = 0
      iERR = 0

 100  nQItrBig = nQItrBig + 1
      If(nQItrBig.GT.MaxBaskets) Then
         iERR = 1000
         goto 1000
      End if


      tm2 = ANItime(tmdata)

      nSkipE = 0
 300  nQItr = nQItr + 1

      If(nQItr.GT.MaxQItr + nQItrAdd) Then
         If(nQItrBig.GE.1) Then
            nQItrAdd = nSkipE
            goto 100
         End if

         iERR = 1000
         goto 1000
      End if

      iwE = L2E(1)
      Do i = 1, nSkipE
         iwE = L1E(2, iwE)
      End do


c ... check the tangled mesh 
      If(qE(L2E(1)).GT.0D0 .AND. flagTM) Then
         flagTM = .FALSE.

         Call delXnode(status, ANIUntangleMesh)

         dao = domainArea(nE + IHolE(1), XYP, IPE)
         dpo = domainPerimetr(nF + IHolF(1), XYP, IPF)
      End if


c ... output statistics
      If((nQItr / 5000) * 5000.EQ.nQItr) Then
         tm2 = ANItime(tmdata)
         If(iPrint.GE.2)
     &      Write(*, 5000) nQItr, qE(L2E(1)), nP, nF, nE, tm2 - tm1
      End if

      rQuality = qE(iwE)
      If(rQuality.GT.Quality) Then
         If(nSkipE.EQ.0) goto 1000
         goto 100
      End if

      iP1 = IPE(1, iwE)
      iP2 = IPE(2, iwE)
      iP3 = IPE(3, iwE)

      Call calQF(HesP(1, iP1), detG(iP1), XYP(1, iP1),
     &           HesP(1, iP2), detG(iP2), XYP(1, iP2),
     &           HesP(1, iP3), detG(iP3), XYP(1, iP3),
     &           hStar, iw)

      Call makSE(iwE, IEP, IPF, IPE, IFE, IEE, qE, MaxS,
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
     &           status)


      Do i = 1, 3
         If(iCLPSop.EQ.1) Then
            Call clpsF1(
c group (M)
     &           iw(i), iwE,
     &           nP, nF, nE,
     &           XYP, IPF, IPE, 
     &           calCrv, parCrv, iFnc,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           L1E, L2E, nL2, nStep,
     &           IHolP, IHolF, IHolE,
     &           status,
c group (CRV)
     &           L1Et, L2Et, tE,
     &           nL2t, nStept, nEt,
     &           nCrvFnc, LFnc, ILt,
C group (Q)
     &           HesP, rQuality, detG, qE,
     &           MetricFunction, flagAnalytic,
C group (S)
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &           nPw, nEw, XYPw, HesPw, IPEw,
     &           iSE, rSE, 5D-1, 5D-1,
     &           flagTest)
            mCLPS1 = mCLPS1 + 1
            If(flagTest) Then
               nCLPS1 = nCLPS1 + 1
               goto 400
            End if
         End if

         If(iINSRTop.EQ.1) Then
            Call insrtP(
c group (M)
     &           iw(i), iwE,
     &           nP, MaxP, nF, MaxF, nE, MaxE,
     &           XYP, IPF, IPE, lbE,
     &           calCrv, parCrv, iFnc,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           L1E, L2E, nL2, nStep,
     &           IHolP, IHolF, IHolE,
     &           status,
c group (CRV)
     &           L1Et, L2Et, tE,
     &           nL2t, nStept, nEt,
     &           nCrvFnc, LFnc, ILt,
C group (Q)
     &           HesP, Quality, rQuality, detG, qE,
     &           MetricFunction, flagAnalytic,
C group (S)
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &           nPw, nEw, XYPw, HesPw, IPEw,
     &           iSE, rSE,
     &           flagTest)
            mINSRT = mINSRT + 1
            If(flagTest) Then
               nINSRT = nINSRT + 1
               goto 400
            End if
         End if

         If(iSWAPop.EQ.1) Then
            Call swapF(
c group (M)
     &           iw(i), iwE,
     &           nE, XYP, IPE,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           L1E, L2E, nL2, nStep,
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
               goto 400
            End if
         End if
      End do


      Do i = 1, 3
         If(iCLPSop.EQ.1) Then
            Call clpsF2(
c group (M)
     &           i, iwE,
     &           nP, nF, nE,
     &           XYP, IPF, IPE, 
     &           calCrv, parCrv, iFnc,
     &           hStar,
     &           ICP, IEP, IFE, IEE,
     &           L1E, L2E, nL2, nStep,
     &           IHolP, IHolF, IHolE,
     &           status,
c group (CRV)
     &           L1Et, L2Et, tE,
     &           nL2t, nStept, nEt,
     &           nCrvFnc, LFnc, ILt,
C group (Q)
     &           HesP, rQuality, detG, qE,
     &           MetricFunction, flagAnalytic,
C group (S)
     &           lE, iEu, 
C group (W)
     &           nPw, nEw, XYPw, HesPw, IPEw,
     &           iSE, rSE,
     &           flagTest)
            mCLPS2 = mCLPS2 + 1
            If(flagTest) Then
               nCLPS2 = nCLPS2 + 1
               goto 400
            End if
         End if


         If(iMOVEop.EQ.1) Then
            Call moveP(
c group (M)
     &           i, iwE,
     &           nE,
     &           XYP, IPF, IPE,
     &           calCrv, parCrv, iFnc,
     &           hStar,
     &           ICP, IEP, IEE,
     &           L1E, L2E, nL2, nStep,
     &           status,
c group (CRV)
     &           L1Et, L2Et, tE,
     &           nL2t, nStept, nEt,
     &           nCrvFnc, LFnc, ILt,
C group (Q)
     &           HesP, rQuality, detG, qE,
     &           MetricFunction, flagAnalytic,
C group (S)
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
C group (W)
     &           nPw, nEw, XYPw, HesPw, IPEw,
     &           iSE, rSE,
     &           icnt, rMove, flagTest)
            mMOVE = mMOVE + 1
            If(flagTest) Then
               nMOVE = nMOVE + 1
               goto 400
            End if
         End if
      End do

      nNOTHING = nNOTHING + 1

 400  Continue 
      If(flagTest) goto 300


      nSkipE = nSkipE + 1
      If(nSkipE.GT.MaxSkipE .OR. nSkipE.GE.nE) goto 100
      goto 300

 1000 rQuality = qE(L2E(1))

      tm2 = ANItime(tmdata)
      If(iPrint.GE.1) 
     &   Write(*, 5000) nQItr - 1, qE(L2E(1)), nP, nF, nE, tm2 - tm1


C ... calculating the number of bad triangles
      If(iPrint.GE.3) Then
         Write(*,5009) tmdata(1), tmdata(2)

         Write(*,5003) nCLPS1, nINSRT, nSWAP, nCLPS2, nMOVE, nNOTHING,
     &                 mCLPS1, mINSRT, mSWAP, mCLPS2, mMOVE

         Call statistics(nE, L1E, L2E, Quality, qE) 
      End if


C ... remove 'holes' from the final grid
      nFo = nF 
      Call updM(
c group (M)
     &            nP, nF, nE, nPv,
     &            XYP, IPF, IPE, IPV,
     &            parCrv, iFnc,
     &            ICP, IEP, IFE, IEE, lbE,
     &            status,
     &            IHolP, IHolF, IHolE,
     &            qE, nEPw)

c     Call draw(nP, nF, nE, XYP, ICP, IPF, IPE, 'fin.ps')
c     Call draw_Q(nP, nE, XYP, IPE, qE, Quality, 'qE.eps')

      Call chkM(nP, nF, nE,
     &            XYP, IPF, IPE, IFE, IEE, lbE,
     &          calCrv, parCrv, iFnc,
     &            ICP, rR, status)

      If(iPrint.GE.1) Then
         averageQ = avgQ(nE, qE, L1E, L2E)
         Write(*, 5006) averageQ, rR, status
      End if

      If(nCrvFnc.EQ.0 .AND. nF.GT.0) Then 
         dao_new = domainArea(nE, XYP, IPE)
         dpo_new = domainPerimetr(nF, XYP, IPF)

         If(iPrint.GE.2) Write(*, 5005) dao, dpo, dao_new, dpo_new

         dao = abs(dao - dao_new) / dao
         dpo = abs(dpo - dpo_new) / dpo

c         If(dao.GT.1D-8) Call errMes(6005, 'ani2.f', 
c     &                       'Lose of the domain area')

c         If(dpo.GT.1D-8 .AND. nF.EQ.nFo) Call errMes(6005, 'ani2.f', 
c     &                       'Lose of the total boundary lenght')
      End if


 5000 Format('ITRs:', I6, '  Q=', E10.4, '   #P#F#E:', 3I8,
     &     '   tm=', F7.2, 's')

 5001 Format('Cleaning basket #', I3, ' with', I4, ' elements')

 5003 Format(/,
     & 'Collapse I-B   Insert   Swapping   Collapse I-E   Moving  Nothin
     &g',/, I12,I9,I11,I15,I9,I9,/,I12,I9,I11,I15,I9)

 5005 Format('Domain area and total      boundary: ', 2E14.6,/,
     &       'Domain area and user-given boundary: ', 2E14.6)

 5006 Format('Avg Quality  =', E11.4, ',  Maximal R/r =', E11.4,
     &       ',  status.fd:', I5)

 5007 Format('Warning:', I6, ' new fix vertices have been added') 
 5008 Format('Warning:', I6, ' new edges have been added') 

 5009 Format('User time:', F7.2, ',  system time:', F7.2)

 9000 Return
      End


  
C ==========================================================
      Subroutine statistics(nE, L1E, L2E, Quality, qE)
C ==========================================================
      Integer L1E(2, *), L2E(*)
      real  Quality, qE(*)

C (Local variables)
      Integer nBad(11)

C ==========================================================
      Do n = 1, 11
         nBad(n) = 0
      End do

      icnt = 0
      iE = L2E(1)
      Do n = 1, nE
         If(qE(iE).LT.Quality) icnt = icnt + 1
         Do k = 1, 10
            If((k - 1) * 1D-1.LE.qE(iE) .AND.
     &                           qE(iE).LT.k * 1D-1) Then
               nBad(k + 1) = nBad(k + 1) + 1
               goto 500
            End if
         End do
         nBad(1) = nBad(1) + 1

 500     iE = L1E(2, iE)
      End do

      Write(*,5003) (1D-1 * i, i = 0, 10), (nBad(i), i = 1, 11)
      If(icnt.GT.0) Write(*,5004) icnt

 5003 Format(/, 'Distribution of triangles by the quality',/,
     &       F5.1,10F7.1,/, I5,10I7)
 
 5004 Format('Number of bad triangles =', I5) 
      Return
      End

