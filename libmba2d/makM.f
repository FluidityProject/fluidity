C ================================================================
      Subroutine makM(
C ================================================================
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
c group (ERR)
     &      iERR)
C ================================================================
      include 'makS.fd'
      include 'colors.fd'
      include 'status.fd'
C ================================================================
C group (M)
C     Integer MaxP, MaxF, MaxE, nPv
      real  XYP(2, *)
      Integer IPE(3, *), IPF(4, *), IPV(*)

      real  parCrv(2, *)
      Integer iFnc(*)

      Integer ICP(*), IEP(*)
      Integer IFE(3, *), IEE(3, *)
      Integer IHolP(*), IHolF(*), IHolE(*)
      Integer IEPw(*), nEPw(*)

c group (Dev)
      Integer nFv, nEv
      Integer IFV(*), IEV(*), lbE(*)
      Integer status

C group (Local variables)
      Integer Mlist(2, MaxS), Clist(MaxS)
      Integer IPFs(2)

      Integer ip(4) 
      Logical ifXnode, cmpP, cmpE, sharpAngle

C ================================================================
      Integer iDomBnd, iMatBnd
      Common /aniBND/ iDomBnd, iMatBnd
C ==========================================================
      iERR = 0

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      IHolP(1) = 0
      IHolF(1) = 0
      IHolE(1) = 0


C ... update parCrv & iFnc
      Do n = nF, 1, -1
         iCrv = IPF(3, n)
         If(iCrv.GT.n) Then
            iERR = 1101
            goto 1000
         End if
         If(iCrv.NE.0) Then
            IPF(3, n) = n
            iFnc(n) = iFnc(iCrv)
            If(n.NE.iCrv) iFnc(iCrv) = 0

            Do i = 1, 2
               parCrv(i, n) = parCrv(i, iCrv)
            End do
         Else
            iFnc(n) = 0
            Do i = 1, 2
               parCrv(i, n) = 0D0
            End do
         End if
      End do


C ... create an auxiliary structure
      Call backReferences(nP, nE, 3, 3, IPE, nEPw, IEPw)

      Do n = 1, MaxE
         Do i = 1, 3
            IEE(i, n) = 0
            IFE(i, n) = 0
         End do
      End do

C ... create IEE & IEP
      Do n = 1, nE
         Do i1 = 1, 3
            i2 = ip(i1 + 1)

            ip1 = IPE(i1, n)
            ip2 = IPE(i2, n)

            IEP(ip1) = n

            If(cmpE(ip1, ip2, IEPw, nEPw, n, iE2)) Then
               IEE(i1, n) = iE2
            End if
         End do
      End do


c ... creating an auxiliary structure
      Call backReferences(nP, nF, 2, 4, IPF, nEPw, IEPw)

      Do n = 1, nP
         ICP(n) = 0
      End do

C ... create IFE  
      Do n = 1, nE
         Do i1 = 1, 3
            i2 = ip(i1 + 1)

            ip1 = IPE(i1, n)
            ip2 = IPE(i2, n)

            If(cmpE(ip1, ip2, IEPw, nEPw, 0, iF)) Then
               IFE(i1, n) = iF
            End if
         End do
      End do


C ... count existing surfaces
      ilist = 0
      icmax = 0
      Do n = 1, nF
         ic1 = IPF(4, n)
         Call findSE(ilist, Clist, ic1, k)

         If(k.EQ.0) Then
            ilist = ilist + 1
            If(ilist.GT.MaxS) Call errMes(1007, 'makM',
     &                             'local parameter MaxS is small')

c  ...  user-given interface
            Do i = 1, 2
               Mlist(i, ilist) = 0
            End do

            Clist(ilist) = ic1
            icmax = max(ic1, icmax) 
         End if
      End do


C ... create temporary surfaces (material and boundary)
      iDomBnd = icmax + 1
      iMatBnd = iDomBnd + 3 * nE

      icbnd = 0
      icmat = 0

      Do n = 1, nE
         Do 100 i1 = 1, 3
            iFt = IFE(i1, n)
            iEt = IEE(i1, n)

            If(iEt.GT.n .OR. iFt.GT.0) goto 100

            ic1 = lbE(n)
            If(iEt.GT.0) Then
               ic2 = lbE(iEt)
            Else If(iFt.EQ.0) Then
               ic2 = 0
            Else
               ic2 = ic1
            End if

c  ...  order materials such that ic1 > ic2
            If(ic1.LT.ic2) Call swapii(ic1, ic2)

            If(ic1.NE.ic2) Then
               Do i = 1, ilist
                  If(Mlist(1, i).EQ.ic1 .AND.
     &               Mlist(2, i).EQ.ic2) Then
                     iBNDs = Clist(i)
                     goto 10
                  End if
               End do

               If(ic2.EQ.0) Then
                  icbnd = icbnd + 1
                  iBNDs = iDomBnd + icbnd
               Else
                  icmat = icmat + 1
                  iBNDs = iMatBnd + icmat
               End if

               ilist = ilist + 1
               If(ilist.GT.MaxS) Call errMes(1007, 'makM',
     &                                'local parameter MaxS is small')
               Mlist(1, ilist) = ic1
               Mlist(2, ilist) = ic2
               Clist(ilist) = iBNDs

 10            Call facAdd(iF, nF, MaxF, IHolF)

               i2 = ip(i1 + 1)

               IPFs(1) = IPE(i1, n)  
               IPFs(2) = IPE(i2, n)

               Call facUpd(1, IPF, parCrv, iFnc,
     &                     (/iF/), IPFs, 0, 0, iBNDs, 0D0, 0D0)

c  ...  update IFE
               IFE(i1, n) = iF
               If(iEt.GT.0) Then
                  Do j1 = 1, 3
                     If(IEE(j1, iEt).EQ.n) IFE(j1, iEt) = iF
                  End do
               End if
            End if
 100     Continue
      End do
      

C ... create ICP
      Do n = 1, nF
         Do i = 1, 2
            iP1 = IPF(i, n)
            Call addXnode(ICP(iP1), jSnode)
         End do
      End do


      Do n = 1, nE
         Do i1 = 1, 3
            iE = IEE(i1, n)
            If(iE.EQ.0) Then
               i2 = ip(i1 + 1)

               iP1 = IPE(i1, n)
               iP2 = IPE(i2, n)

               Call addXnode(ICP(iP1), jBnode)
               Call addXnode(ICP(iP2), jBnode)
            End if
         End do
      End do

      Do n = 1, nP
         If(.NOT.ifXnode(ICP(n), jBnode)) Call addXnode(ICP(n), jInode)
      End do

      Do n = 1, nPv
         iP1 = IPV(n)
         Call addXnode(ICP(iP1), jVnode)
      End do


c ... first realization of the fixed edges and triangles
      Do n = 1, nFv
         iFt = IFV(n)

         Do i = 1, 2
            iPt = IPF(i, iFt)
            Call addXnode(ICP(iPt), jTnode)
         End do
      End do

      Do n = 1, nEv
         iEt = IEV(n)

         Do i = 1, 3
            iPt = IPE(i, iEt)
            Call addXnode(ICP(iPt), jTnode)
         End do
      End do


c ... color the points (vertices)
      Call backReferences(nP, nF, 2, 4, IPF, nEPw, IEPw)

      Do n = 1, nP
         If(cmpP(n, IEPw, nEPw, IPF)) Call addXnode(ICP(n), jVnode)
      End do


c ... color sharp angles
      Do n = 1, nP
         If(.NOT.ifXnode(ICP(n), jVnode)) Then
            If(sharpAngle(n, IEPw, nEPw, XYP, IPF)) 
     &         Call addXnode(ICP(n), jVnode)
         End if
      End do


c ... color the points (fixed boundary points)
      If(ifXnode(status, ANIFixBoundaryPoints)) Then
         Do n = 1, nE
            Do i1 = 1, 3
               If(IEE(i1, n).EQ.0) Then
                  i2 = ip(i1 + 1)

                  iP1 = IPE(i1, n)
                  iP2 = IPE(i2, n)
                  Call addXnode(ICP(iP1), jVnode)
                  Call addXnode(ICP(iP2), jVnode)
               End if
            End do
         End do
      End if


c ... color the points (fixed boundary edges)
      If(ifXnode(status, ANIFixBoundaryEdges)) Then
         Do n = 1, nE
            Do i1 = 1, 3
               If(IEE(i1, n).EQ.0) Then
                  i2 = ip(i1 + 1)

                  iP1 = IPE(i1, n)
                  iP2 = IPE(i2, n)
                  Call addXnode(ICP(iP1), jTnode)
                  Call addXnode(ICP(iP2), jTnode)
               End if
            End do
         End do
      End if

 1000 Return
      End



C ================================================================
      Subroutine updM(
C ================================================================
c group (M)
     &            nP, nF, nE, nPv,
     &            XYP, IPF, IPE, IPV,
     &            parCrv, iFnc,
     &            ICP, IEP, IFE, IEE, lbE,
     &            status,
     &            IHolP, IHolF, IHolE,
     &            qE, IPw)
C ================================================================
      include 'status.fd'
C ================================================================
C *** Remarks:
C     1. size of working memory is equal to MaxP, IPw(MaxP)
C ================================================================
C group (M)
      Integer IPF(4, *), IPE(3, *), IPV(*)
      real  XYP(2, *)

      Integer iFnc(*)
      real  parCrv(2, *)

      Integer ICP(*), IEP(*), IFE(3, *), IEE(3, *), lbE(*)
      Integer status

      Integer IHolP(*), IHolF(*), IHolE(*)

C group (Q)
      real  qE(*)

C group (W)
      Integer IPw(*)

C group (Local variables)
      Logical ifXnode
 
C ================================================================
      Integer iDomBnd, iMatBnd
      Common /aniBND/ iDomBnd, iMatBnd
C ================================================================
c ... delete references to material or fictitious faces (Id >= iDomBnd)
      If(ifXnode(status, ANIDeleteTemporaryEdges)) Then
         mF = nF + IHolF(1)
         mE = nE + IHolE(1)

         Do n = 1, mE
            If(IPE(1, n).GT.0) Then
               Do i = 1, 3
                  iFt = IFE(i, n)
                  If(iFt.GT.0) Then
                     If(IPF(4, iFt).GE.iDomBnd) IFE(i, n) = 0
                  End if
               End do
            End if
        End do

         Do n = 1, mF
            If(IPF(1, n).GT.0 .AND. IPF(4, n).GE.iDomBnd) Then
               Call facDel(n, nF, IPF, iFnc, IHolF)
            End if
         End do
      End if


      nP = nP + IHolP(1)
      nF = nF + IHolF(1)
      nE = nE + IHolE(1)


c ... fill in holes in points
      Do n = 1, nP
         IPw(n) = 0
      End do

      lP = IHolP(1)
      Do n = 1, lP
         iP = IHolP(n + 1)
         IPw(iP) = -1
      End do

      mP = 0
      Do n = 1, nP
         If(IPw(n).EQ.0) Then
            mP = mP + 1
            IPw(n) = mP

            Do i = 1, 2
               XYP(i, mP) = XYP(i, n)
            End do

            ICP(mP) = ICP(n)
         End if
      End do


      Do n = mP + 1, nP
         ICP(n) = 0
      End do

      Do n = 1, nPv
         IPV(n) = IPw(IPV(n))
      End do

      Do n = 1, nF
         If(IPF(1, n).GT.0) Then
            Do i = 1, 2
               IPF(i, n) = IPw(IPF(i, n))
            End do
         End if
      End do

      Do n = 1, nE
         If(IPE(1, n).GT.0) Then
            Do i = 1, 3
               IPE(i, n) = IPw(IPE(i, n))
            End do
         End if
      End do

      nP = mP


c ... fill in holes in edges
      lF = IHolF(1)
      Do 200 n = 1, lF
         iF = IHolF(n + 1)

         Do m = nF, iF + 1, -1
            If(IPF(1, m).NE.0) Then
               kF = m
               goto 20
            End if
         End do
         goto 200

 20      Do i = 1, 4
            IPF(i, iF) = IPF(i, kF)
         End do

         Do i = 1, 2
            parCrv(i, iF) = parCrv(i, kF)
         End do

         iFnc(iF) = iFnc(kF)

         Do k = 1, nE
            Do i = 1, 3
               If(IFE(i, k).EQ.kF) IFE(i, k) = iF
            End do
         End do

         IPF(1, kF) = 0
 200     nF = nF - 1

      icnt = 0
      Do n = 1, nF
         iCrv = IPF(3, n)
         If(iCrv.NE.0) Then
            icnt = icnt + 1

            IPF(3, n) = icnt

            iFnc(icnt) = iFnc(n)
            Do i = 1, 2
               parCrv(i, icnt) = parCrv(i, n)
            End do
         End if
      End do


c ... fill in holes in elements
      lE = IHolE(1)
      Do 300 n = 1, lE
         iE = IHolE(n + 1)

         Do m = nE, iE + 1, -1
            If(IPE(1, m).NE.0) Then
               kE = m
               goto 30
            End if
         End do
         goto 300

 30      Do i = 1, 3
            IPE(i, iE) = IPE(i, kE)
         End do

         qE(iE) = qE(kE)

         lbE(iE) = lbE(kE)

c  ...   update auxiliary structures
         Do i = 1, 3
            IFE(i, iE) = IFE(i, kE)
            IEE(i, iE) = IEE(i, kE)
            iEt = IEE(i, iE)
            If(iEt.GT.0) Then
              Do j = 1, 3
                 If(IEE(j, iEt).EQ.kE) IEE(j, iEt) = iE
              End do
            End if
         End do


         IPE(1, kE) = 0
 300     nE = nE - 1

      Return
      End



C ================================================================
      Subroutine chkM(
C ================================================================
c group (M)
     &            nP, nF, nE,
     &            XYP, IPF, IPE, IFE, IEE, lbE,
     &            calCrv, parCrv, iFnc,
     &            ICP, rR, status)
C ================================================================
      include 'colors.fd'
      include 'status.fd'
C ================================================================
C Routines checks topology of the input and the oputput meshes.
C ================================================================
C group (M)
      Integer IPF(4, *), IPE(3, *), IFE(3, *), IEE(3, *), lbE(*)
      real  XYP(2, *)

      EXTERNAL calCrv
      Integer iFnc(*)
      real  parCrv(2, *)

      Integer ICP(*), status
      real  rR

C group (Local variables)
      Integer ip(4)
      real  XYPs(2), err, rOut, rIn
      real  sqrEdge, calVol
      Logical ifXnode, check22

C ================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

c ... check color of mesh points
      Do n = 1, nP
         If(ICP(n).LE.0) 
     &      Call errMes(5001, 'chkM', 'wrong point color')

         If(ifXnode(ICP(n), jInode) .AND. 
     &      ifXnode(ICP(n), jBnode)) 
     &      Call errMes(5015, 'chkM', 'color contradiction')
      End do


c ... check labels of mesh faces
      Do n = 1, nF
         If(IPF(1, n).LE.0) 
     &      Call errMes(5002, 'chkM', 'wrong map edge -> points')

         If(IPF(3, n).NE.0 .AND. iFnc(n).LE.0) 
     &      Call errMes(5003, 'chkM', 'wrong Id of curvilinear edge')
      End do


c ... check elements
      Do n = 1, nE
         nTet = n

         Do i = 1, 3
            If(IPE(1, n).LE.0) 
     &         Call errMes(5004, 'chkM', 'wrong map element -> points')
         End do

         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)
         If(iP1.EQ.iP2 .OR. iP1.EQ.iP3 .OR. iP2.EQ.iP3) Then
            iERR = 5006
            goto 500
         End if

         If(lbE(n).LE.0) 
     &      Call errMes(5016, 'chkM', 
     &                 'element identificator (lbE) is not positive')

         iF1 = IFE(1, n)
         iF2 = IFE(2, n)
         iF3 = IFE(3, n)

         If(iF1.EQ.iF2 .AND. iF1.GT.0 .OR.
     &      iF1.EQ.iF3 .AND. iF1.GT.0 .OR.
     &      iF2.EQ.iF3 .AND. iF2.GT.0) Then
            iERR = 5007
            goto 500
         End if


         iE1 = IEE(1, n)
         iE2 = IEE(2, n)
         iE3 = IEE(3, n)

         If(iE1.EQ.iE2 .AND. iE1.NE.0 .OR.
     &      iE1.EQ.iE3 .AND. iE1.NE.0 .OR.
     &      iE2.EQ.iE3 .AND. iE2.NE.0) Then
            iERR = 5008
            goto 500
         End if


         Do 20 i1 = 1, 3
            iF = IFE(i1, n)
            iE = IEE(i1, n)
            If(iF.EQ.0 .AND. iE.EQ.0 .AND. nF.GT.0) Then
               iERR = 5009
               goto 500
            End if

            i2 = ip(i1 + 1)

            iP1 = IPE(i1, n)
            iP2 = IPE(i2, n)

            If(iF.NE.0) Then
               jP1 = IPF(1, iF)
               jP2 = IPF(2, iF)

               If(.NOT.check22(iP1, iP2, jP1, jP2)) Then
                  iERR = 5010
                  goto 500
               End if
            End if

            If(iE.NE.0) Then
               If(iF.NE.0) Then
                  Do j1 = 1, 3
                     If(IFE(j1, iE).EQ.iF) goto 10
                  End do

                  iERR = 5011
                  goto 500
               End if

  10           Do j1 = 1, 3
                  j2 = ip(j1 + 1)

                  jP1 = IPE(j1, iE)
                  jP2 = IPE(j2, iE)

                  If(check22(iP1, iP2, jP1, jP2)) Then
                     If(IEE(j1, iE).NE.n) Then
                        iERR = 5012
                        goto 500
                     End if

                     i3  = ip(i2 + 1)
                     iP3 = IPE(i3, n)

                     j3  = ip(j2 + 1)
                     jP3 = IPE(j3, iE)

                     v1 = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))
                     v2 = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, jP3))

c  ...   check for inverted elements
                     If(v1 * v2.GE.0D0 .AND.
     &                  .NOT.ifXnode(status, ANIUntangleMesh)) Then
                        iERR = 5013
                        goto 500
                     End if
                     goto 20
                  End if
               End do

               iERR = 5014
               goto 500
            End if
 20      Continue
      End do


c ... check the parametrization
      Do n = 1, nF
         If(IPF(3, n).GT.0) Then
            iCrv = IPF(3, n)
            Do i = 1, 2
               iPt = IPF(i, n)
               Call aniCrv(parCrv(i, iCrv), XYPs, iFnc(iCRV), calCrv)
               err = sqrEdge(XYP(1, iPt), XYPs)

               If(err.GE.1D-16)
     &            Call errMes(5006, 'chkM', 'wrong parametrization')
            End do
         End if
      End do   


c ... compute R / r
      rR = 0D0
      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         Call RandR(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), rOut, rIn)

         rR = max(rR, rOut / rIn)
      End do

      Return


 500  Write(*, 5000) nTet, iERR, (IPE(i, nTet), i = 1, 3),
     &                           (IFE(i, nTet), i = 1, 3),
     &                           (IEE(i, nTet), i = 1, 3),
     &                           (ICP(IPE(i, nTet)), i = 1, 3)
      Do k = 1, 3
         iE = IEE(k, nTet)
         If(iE.GT.0) Then
            Write(*, 5000) iE, iERR, (IPE(i, iE), i = 1, 3),
     &                               (IFE(i, iE), i = 1, 3),
     &                               (IEE(i, iE), i = 1, 3),
     &                               (ICP(IPE(i, iE)), i = 1, 3)
         End if
      End do
      Do k = 1, 3
         iF = IFE(k, nTet)
         If(iF.GT.0) Then
            Write(*, 5004) iF, nF, (IPF(i, iF), i = 1, 2),
     &                             (ICP(IPF(i, max(1, iF))), i = 1, 2)
         End if
      End do
      Call errMes(iERR, 'chkM', 'triangles are wrong')

 5000 Format('Error in checking triangle =', I6, '   iERR=', I5, /,
     &       'Points    =', 3I6, /,
     &       'Edges     =', 3I6, /,
     &       'Triangles =', 3I6, /,
     &       'colors    =', 3I6, /)

 5004 Format('Edge ', I6, '  (', I5, ')  of the bad triangle', /,
     &       'Points =', 2I6, /,
     &       'colors =', 2I6, /)

      End



C ================================================================
      Logical Function cmpE(i1, i2, IEP, nEP, iE1, iE2)
C ================================================================
C Search for element iE2, other than iE1, in the intersection of
C sublists associated with i1 and i2. The routine returns value
C .FALSe. and iE2 = 0 when thre is no such element. 
C ================================================================
      Integer IEP(*), nEP(*)

C group (Local variables)
      Integer ib(2), ie(2), ip(2)

      ip(1) = i1
      ip(2) = i2
      Do i = 1, 2
         If(ip(i).EQ.1) Then
            ib(i) = 1
         Else
            ib(i) = nEP(ip(i) - 1) + 1
         End if
         ie(i) = nEP(ip(i))
      End do

      Do 10 i = ib(1), ie(1)
         iE2 = IEP(i)
         If(iE2.EQ.iE1) goto 10
         Do j = ib(2), ie(2)
            If(iE2.EQ.IEP(j)) Then
               cmpE = .TRUE.
               goto 1000
            End if
         End do
 10   Continue

      iE2 = 0
      cmpE = .FALSE.
 1000 Return
      End



C ================================================================
      Logical Function cmpF(i1, IFP, nFP, iF1, iF2)
C ================================================================
C Search of one coincident (different from iF1) in two lists
C ================================================================
      Integer IFP(*), nFP(*)

      ip = i1
      If(ip.EQ.1) Then
         ib = 1
      Else
         ib = nFP(ip - 1) + 1
      End if
      ie = nFP(ip)

      Do 10 i = ib, ie
         iF2 = IFP(i)
         If(iF2.EQ.iF1) goto 10

         cmpF = .TRUE.
         goto 1000
 10   Continue

      cmpF = .FALSE.
 1000 Return
      End



C ================================================================
      Logical Function cmpP(iP, IFP, nFP, IPF)
C ================================================================
C cmpP = TRUE if point iP belongs to two faces 
C with different colors. Otherwise cmpP = FALSE.
C
C Remark: new structure of IFP, nFP
C ================================================================
      Integer IFP(*), nFP(*), IPF(4, *)

      If(iP.EQ.1) Then
         ib = 1
      Else
         ib = nFP(iP - 1) + 1
      End if
      ie = nFP(iP)

      Do i = ib, ie
         ICF1 = IPF(4, IFP(i))
         Do j = i + 1, ie
            ICF2 = IPF(4, IFP(j))
            If(ICF1.NE.ICF2) Then
               cmpP = .TRUE.
               goto 1000
            End if
         End do
      End do

      cmpP = .FALSE.
 1000 Return
      End



C ================================================================
      Logical Function sharpAngle(n, IFP, nFP, XYP, IPF)
C ================================================================
      include 'magic.fd'
C ================================================================
C Routine checks that the angle between two surface edges is
C less than MaxSharpAngle.
C ================================================================
      Integer IFP(*), nFP(*), IPF(4, *)
      real  XYP(2, *)

c group (Local variabels)
      real  DotMul, calNorm
      real  v1(2), v2(2), angle, norms

C ================================================================
      sharpAngle = .FALSE.

      If(n.EQ.1) Then
         ib = 1
      Else
         ib = nFP(n - 1) + 1
      End if
      ie = nFP(n)

      nFt = ie - ib + 1
      If(nFt.LT.2) goto 1000
      If(nFt.GT.2) goto 500


c ... find two neighboring points for point n
      iF1 = IFP(ib)
      iF2 = IFP(ie)

      Do i = 1, 2
         iPt = IPF(i, iF1)
         If(iPt.NE.n) iP1 = iPt

         iPt = IPF(i, iF2)
         If(iPt.NE.n) iP2 = iPt
      End do

      Do i = 1, 2
         v1(i) = XYP(i, iP1) - XYP(i, n)
         v2(i) = XYP(i, iP2) - XYP(i, n)
      End do  

      angle = DotMul(v1, v2)
      norms = calNorm(v1) * calNorm(v2)

      If(angle.LE.MaxSharpAngle*norms) goto 1000
   
  500 sharpAngle = .TRUE.

 1000 Return
      End



C ==============================================================
      Subroutine RandR(XY1, XY2, XY3, rOut, rIn)
C ==============================================================
C Routine computes curcumscribed and inscribed radii for the 
c triangle defined by three vertices.
C ==============================================================
      real  XY1(3), XY2(3), XY3(3)
      real  rOut, rIn, calEdge

      real  a, b, c, s
C ==============================================================
      a = calEdge(XY1, XY2)
      b = calEdge(XY1, XY3)
      c = calEdge(XY2, XY3)

      s = (a + b + c) / 2
 
      rOut = a * b * c / (4 * sqrt(s * (a + b - s) * 
     &                                  (a + c - s) * (b + c - s)))

      rIn = a * b * c / (4 * ROut * s)

      Return
      End



C ================================================================
      real  Function domainArea(nE, XYP, IPE)
C ================================================================
      real  XYP(2, *)
      Integer IPE(3, *)

c (Local variables)
      real s

      s = 0D0
      Do 10 n = 1, nE
         If(IPE(1, n).EQ.0) goto 10

         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         s = s + abs((XYP(1, iP2) - XYP(1, iP1)) *
     &                (XYP(2, iP3) - XYP(2, iP1)) -
     &                (XYP(1, iP3) - XYP(1, iP1)) *
     &                (XYP(2, iP2) - XYP(2, iP1)))
 10   Continue

      domainArea = 5D-1 * s 
      Return
      End



C ================================================================
      real  Function domainPerimetr(nF, XYP, IPF)
C ================================================================
      real  XYP(2, *)
      Integer IPF(4, *)

c (Local variables)
      real s

      s = 0D0
      Do 10 n = 1, nF
         If(IPF(1, n).EQ.0) goto 10

         iP1 = IPF(1, n)
         iP2 = IPF(2, n)

         s = s + sqrt((XYP(1, iP2) - XYP(1, iP1)) ** 2 +
     &                 (XYP(2, iP2) - XYP(2, iP1)) ** 2)
 10   Continue

      domainPerimetr = s 
      Return
      End



