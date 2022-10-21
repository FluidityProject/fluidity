      Module mba3d_refine
C
      use mba3d_auxSE
      use mba3d_cmpE
      use mba3d_error
      use mba3d_makQ
      use mba3d_utils
C
      contains
C
C ================================================================
      Subroutine initializeRefinement(
C ================================================================
     &      nP, nE, XYP, IPE, 
     &      MapMtr, Ref2MapMtr)
C ================================================================
C   Fill mapping matrix for each tetrahedron to make it equilateral
C   Fill reference array to MapMtr
C ================================================================
      Real*8  XYP(3, *), MapMtr(3, 3, *)
      Integer IPE(4, *), Ref2MapMtr(*)

C Local
      Real*8  AA(3,3), A(3,3), B(3,3), XYPet(3,4)
      Integer IPIV(3), info

c Equilateral tet
      XYPet(1,1) = 0d0
      XYPet(2,1) = 0d0
      XYPet(3,1) = 0d0
      XYPet(1,2) = 2*dsqrt(2d0)/3
      XYPet(2,2) = 0d0
      XYPet(3,2) = -4d0/3
      XYPet(1,3) = -dsqrt(2d0)/3
      XYPet(2,3) =  dsqrt(6d0)/3
      XYPet(3,3) = -4d0/3
      XYPet(1,4) = -dsqrt(2d0)/3
      XYPet(2,4) = -dsqrt(6d0)/3
      XYPet(3,4) = -4d0/3

c Find map Ortotet to Equitet
      Do i=1,3
         Do j=1,3
            AA(i,j) = XYPet(i,j+1) - XYPet(i,1)
         End do
      End do

c Loop over tets
      Do n = 1, nE
c ...  find map  Ortotet to Tet
         Do i=1,3
            Do j=1,3
               A(i,j) = XYP(i, IPE(j+1, n)) - XYP(i, IPE(1, n))
            End do
         End do

c ...  find map Tet to Ortotet
         Do i=1,3
            Do j=1,3
               B(i,j) = 0d0
            End do
            B(i,i) = 1d0
         End do

         Call dgesv(3, 3, A, 3, IPIV, B, 3, info)
         If (info.ne.0) Call errMes(1001,'InitializeRefinement',
     &                              'error in dgesv')

c ...  find map Tet to Equitet
         Do i=1,3
            Do j=1,3
               MapMtr(i,j,n) = 0d0
               Do k=1,3
                  MapMtr(i,j,n) = MapMtr(i,j,n) + AA(i,k)*B(k,j)
               End do
            End do
         End do

         Ref2MapMtr(n) = n
      End do

      Return
      End Subroutine initializeRefinement


C ================================================================
C This routine may be called ONLY after initializeRefinement
C ================================================================
      Subroutine uniformRefinement(
C ================================================================
     &      nP, MaxP, nF, MaxF, nE, MaxE,  
     &      XYP, IPF, IPE, lbF, lbE,
     &      MapMtr, Ref2MapMtr,
     &      iW, MaxWi)
C ================================================================
C    MapMtr - data array of length 9*nEcoarse
C    Ref2MapMtr - data array (lengths on input and output are 9*nE)
C    iW(*) - working array of size MaxWi which is at least 
C            nP + 14 * nE 
C ================================================================
      Real*8  XYP(3, *), MapMtr(9, *)
      Integer IPE(4, *), IPF(3, *), lbF(*), lbE(*), Ref2MapMtr(*)

      Integer iW(*)

C (Local variables)
      Integer IPEs(4)
      Real*8  XYPt(3, 4), Qt, minQt(3)

      Logical flagE
      Integer refCrnr(4,4), refOct(4,4,3), loop(4)

C ================================================================
      DATA   refOct /4,5,6,1, 2,3,6,1, 1,3,5,6, 1,2,4,6,
     &               1,3,2,5, 3,6,2,5, 6,4,2,5, 4,1,2,5,
     &               1,2,3,4, 2,6,3,4, 6,5,3,4, 5,1,3,4/
      DATA   refCrnr/1,2,3,1, 1,4,5,2, 2,4,6,3, 3,5,6,4/
      DATA   loop/1, 2, 3, 1/

C ================================================================
      inEP = 1
      iIEP = inEP + nP
      iIRE = iIEP + 4 * nE 
      iIPE = iIRE + 6 * nE
      iEnd = iIPE + 4 * nE


      Call listE2R(nP, nR, nE, IPE, iW(iIRE), iW(inEP), iW(iIEP))


      nPo = nP
      nFo = nF
      nEo = nE

      nP = nPo + nR
      nF = 4 * nFo
      nE = 8 * nEo

      If(nP.GT.MaxP) Call errMes(1003, 'uniformRefinement',
     &                          'local parameter MaxP is small')
      If(nF.GT.MaxF) Call errMes(1004, 'uniformRefinement',
     &                          'local parameter MaxF is small')
      If(nE.GT.MaxE) Call errMes(1006, 'uniformRefinement',
     &                          'local parameter MaxE is small')


c ... mapping P -> E
      Call backReferences(nPo, nEo, 4, 4, IPE, iW(inEP), iW(iIEP))
      k = 0
      Do n = 1, nEo
         Do i = 1, 4
            k = k + 1
            iW(iIPE + k - 1) = IPE(i, n)
         End do
      End do


c ... splitting elements
      kR = 0
      Do n = 1, nEo
         Do i1 = 1, 3
            Do i2 = i1 + 1, 4
               iP1 = IPE(i1, n)
               iP2 = IPE(i2, n)

               kR = kR + 1
               iR = iW(iIRE + kR - 1)

               Do i = 1, 3
                  XYP(i, nPo + iR) = (XYP(i, iP1) + XYP(i, iP2)) / 2
               End do
            End do
         End do

         Do i = 1, 4
            IPEs(i) = IPE(i, n)
         End do

c ... find optimal splitting of inner octahedron (ii=1,2,3)
         Do ii = 1, 3
          minQt(ii)=1d0
          Do i = 1, 4
            Do k = 1, 4
               j = nPo + iW(iIRE + 6 * n + refOct(k,i,ii) - 7)
               jj= Ref2MapMtr(n) 
               if (jj.lt.0) jj = -jj 
               call FixedMap(XYP(1, j), MapMtr(1, jj), XYPt(1, k))
            End do
            call calREG(XYPt(1,1),XYPt(1,2),XYPt(1,3),XYPt(1,4), Qt)
            minQt(ii) = min( Qt, minQt(ii) )
          End do
         End do
         If (minQt(1).ge.minQt(2).and.minQt(1).ge.minQt(3)) Then
             ii = 1
         Else if (minQt(2).ge.minQt(1).and.minQt(2).ge.minQt(3)) Then
             ii = 2
         Else if (minQt(3).ge.minQt(1).and.minQt(3).ge.minQt(2)) Then
             ii = 3
         End if

c ... apply optimal splitting
         Do i = 1, 4
            iE = (i + 3) * nEo + n
            Do k = 1, 4
               IPE(k, iE) = nPo + iW(iIRE + 6 * n + refOct(k,i,ii) - 7)
            End do
            lbE(iE) = lbE(n)
            Ref2MapMtr(iE) = Ref2MapMtr(n)

            iE = (4 - i) * nEo + n
            Do k = 1, 3
               IPE(k, iE) = nPo + iW(iIRE + 6 * n + refCrnr(k, i) - 7)
            End do
            IPE(4, iE) = IPEs(i)
            lbE(iE) = lbE(n)
            Ref2MapMtr(iE) = Ref2MapMtr(n)
         End do
      End do


c ... splitting faces
      Do n = 1, nFo
         iP1 = IPF(1, n)
         iP2 = IPF(2, n)
         iP3 = IPF(3, n)

         flagE = cmpE(iP1, iP2, iP3, iW(iIEP), iW(inEP), 0, iE)

         jF = 3 * nFo + n

         k = 0
         Do i1 = 1, 3
            Do 10 i2 = i1 + 1, 4
               k = k + 1
               iP1 = iW(iIPE + 4 * iE + i1 - 5)
               iP2 = iW(iIPE + 4 * iE + i2 - 5)

               Do j1 = 1, 3
                  j2 = loop(j1 + 1)

                  jP1 = IPF(j1, n)
                  jP2 = IPF(j2, n)

                  If(check22(iP1, iP2, jP1, jP2)) Then
                     IPF(j1, jF) = nPo + iW(iIRE + 6 * iE + k - 7)
                     Goto 10
                  End if
               End do
  10        Continue 
         End do
         lbF(jF) = lbF(n)

         Do i1 = 1, 3
            iF = (3 - i1) * nFo + n 

            i2 = loop(i1 + 1)
            i3 = loop(i2 + 1)

            IPF(1, iF) = IPF(i1, n)
            IPF(2, iF) = IPF(i1, jF)
            IPF(3, iF) = IPF(i3, jF)
            lbF(iF) = lbF(n)
         End do
      End do

      Return
      End Subroutine uniformRefinement



C ================================================================
C This routine may be called ONLY after initializeRefinement
C ================================================================
      Subroutine localRefinement(
C ================================================================
     &      nP, MaxP, nF, MaxF, nE, MaxE,  
     &      XYP, IPF, IPE, lbF, lbE,
     &      MapMtr, Ref2MapMtr, SplitFlag,
     &      iW, MaxWi)
C ================================================================
C  Routine refines uniformly elements marked by SplitFlag.
C  SplitFlag=.true. is ignored in case of negative sign of Ref2MapMtr.
C  Ref2MapMtr is assigned  negative sign at an element if it was split non-uniformly.
C
C  Remarks:
C    MapMtr - data array of length 9*nEcoarse
C    Ref2MapMtr - data array (lengths on input and output are 9*nE)
C    iW(*) - working array of size MaxWi which is at least 
C            nP + 14 * nE + 3*nR + 4*min(MaxF,4*nF)
C ================================================================
      Real*8  XYP(3, *), MapMtr(9, *)
      Integer IPE(4, *), IPF(3, *), lbF(*), lbE(*), Ref2MapMtr(*)
      Logical SplitFlag(*)

      Integer iW(*)

C (Local variables)
      Integer IPEs(6)
      Real*8  XYPt(3, 4), Qt, minQt(3)

      Logical flagE
      Integer refCrnr(4,4), refOct(4,4,3), loop(4), 
     &        facedge(3,4), facount(4), medge(6), opposite(3)

C ================================================================

      DATA   refOct /4,5,6,1, 2,3,6,1, 1,3,5,6, 1,2,4,6,
     &               1,3,2,5, 3,6,2,5, 6,4,2,5, 4,1,2,5,
     &               1,2,3,4, 2,6,3,4, 6,5,3,4, 5,1,3,4/
      DATA   refCrnr/1,2,3,1, 1,4,5,2, 2,4,6,3, 3,5,6,4/
      DATA   loop/1, 2, 3, 1/
      DATA   facedge/1,2,4, 2,3,6, 1,3,5, 4,5,6/

C ================================================================
      iIRE = 1
      iIPE = iIRE + 6 * nE
      inEP = iIPE + 4 * nE 
      iIEP = inEP + nP
      iEnd = iIEP + 4 * nE 


      Call listE2R(nP, nR, nE, IPE, iW(iIRE), iW(inEP), iW(iIEP))

      ilbR   = iEnd
      iendsR = ilbR + nR
      iIPF   = iendsR + 2*nR
      iEnd   = iIPF + 4*min(MaxF,4*nF)

      If(iEnd.GT.MaxWi) Call errMes(1001, 'localRefinement',
     &                             'not enough working memory')

      nPo = nP
      nFo = nF
      nEo = nE


c ... mark new points on edges
      Do n = 1, nR
         iW(ilbR + n - 1) = 0 
      End do

c ... mark mother edges
      Do n = 1, nE
         If( SplitFlag(n).and.Ref2MapMtr(n).GT.0 ) Then
            Do i = 1, 6
               iRt = iW(iIRE + 6 * (n-1) + i - 1)
               iW(ilbR + iRt - 1) = 1
            End do
         End if
      End do

c ... count mother edges
      mR = 0
      Do iRt = 1, nR
         If(iW(ilbR + iRt - 1).GT.0) Then 
            nP = nP + 1
            mR = mR + 1
            iW(ilbR + iRt - 1) = mR
         End if
      End do

c ... assign father nodes to mother edges
      Do n = 1, nE
         k = 0
         Do i1 = 1, 3
            Do i2 = i1 + 1, 4
               iP1 = IPE(i1, n)
               iP2 = IPE(i2, n)

               k = k + 1
               iRt = iW(iIRE + 6 * (n-1) + k - 1)
               iRc = iW(ilbR + iRt - 1)
               If(iRc.GT.0) Then
                  iW(iendsR + (iRc - 1)*2)      = iP1
                  iW(iendsR + (iRc - 1)*2 + 1)  = iP2
               End if
            End do
         End do
      End do


c ... add more points so that each face have 0, 1, or 3 new points
      Do n = 1, nEo
         k = 0

c        count new points for faces of each tet
         Do j1 = 1, 4
            facount(j1) = 0
         End do

         Do i1 = 1, 3          
            Do i2 = i1 + 1, 4  ! loop over edges
               k = k + 1
               iRt = iW(iIRE + 6 * (n-1) + k - 1)
               iRc = iW(ilbR + iRt - 1)

               If(iRc.GT.0) Then
c                 Loop over faces (j1) having k-th edge
                  Do j1 = 1, 4
                     Do j2 = 1, 3
                        If(facedge(j2,j1).eq.k) Then
                           facount(j1) = facount(j1) + 1
                        End if
                     End do
                  End do
               End if
            End do
         End do

c        Add new point to faces with 2 new points
         k = 0
         Do i1 = 1, 3
            Do i2 = i1 + 1, 4 ! loop over edges
               k = k + 1
               iRt = iW(iIRE + 6 * (n-1) + k - 1)
               iRc = iW(ilbR + iRt - 1)

               If(iRc.EQ.0) Then
c                 Loop over faces (j1) having k-th edge
                  Do j1 = 1, 4
                     Do j2 = 1, 3
                        If(facedge(j2,j1).eq.k.and.facount(j1).eq.2)Then
c                        add extra mother edges
                           nP = nP + 1
                           mR = mR + 1
                           iW(ilbR + iRt - 1) = mR
                           If (k.eq.1) Then
                               iP1 = IPE(1,n)
                               iP2 = IPE(2,n)
                           Else If (k.eq.2) Then
                               iP1 = IPE(1,n)
                               iP2 = IPE(3,n)
                           Else If (k.eq.3) Then
                               iP1 = IPE(1,n)
                               iP2 = IPE(4,n)
                           Else If (k.eq.4) Then
                               iP1 = IPE(2,n)
                               iP2 = IPE(3,n)
                           Else If (k.eq.5) Then
                               iP1 = IPE(2,n)
                               iP2 = IPE(4,n)
                           Else If (k.eq.6) Then
                               iP1 = IPE(3,n)
                               iP2 = IPE(4,n)
                           End if

                           iW(iendsR + (mR - 1)*2)     = iP1
                           iW(iendsR + (mR - 1)*2 + 1) = iP2
                           facount(j1) = 3

                           Do k1 = 1, 4
                              Do k2 = 1, 3
                                If(facedge(k2,k1).eq.k.and.k1.ne.j1)Then
                                   facount(k1) = facount(k1) + 1
                                End if
                              End do
                           End do
                           Goto 1
                        End if
                     End do
                  End do
               End if

1              Continue
            End do
         End do
      End do


c ... check for faces with 2 new points
      Do n = 1, nEo
         k = 0

c        count new points for faces of each tet
         Do j1 = 1, 4
            facount(j1) = 0
         End do

         Do i1 = 1, 3
            Do i2 = i1 + 1, 4
               k = k + 1
               iRt = iW(iIRE + 6 * (n-1) + k - 1)
               iRc = iW(ilbR + iRt - 1)

               If(iRc.GT.0) Then
c                 Loop over faces (j1) having k-th edge
                  Do j1 = 1, 4
                     Do j2 = 1, 3
                        If(facedge(j2,j1).eq.k) Then
                           facount(j1) = facount(j1) + 1
                        End if
                     End do
                  End do
               End if
            End do
         End do

         Do j1 = 1, 4
            If(facount(j1).eq.2) Call errMes(1002, 'localRefinement',
     &                          'facount=2')
c If this error happen, it implies that algorithm for removing situation 
c with 2 new points at a face has to be refined!
         End do
      End do



      If(nP.GT.MaxP) Call errMes(1003, 'localRefinement',
     &                          'local parameter MaxP is small')


c ... create new points, the point number is nPo plus the edge number
      Do iRt = 1, nR
         iRc = iW(ilbR + iRt - 1)
         If(iRc.GT.0) Then
            iPt = nPo + iRc
            iP1 = iW(iendsR + (iRc - 1)*2)
            iP2 = iW(iendsR + (iRc - 1)*2 + 1)
            Do i = 1, 3
               XYP(i, iPt) = (XYP(i, iP1) + XYP(i, iP2)) / 2
            End do
         End if
      End do



c ... mapping P -> E
      Call backReferences(nPo, nEo, 4, 4, IPE, iW(inEP), iW(iIEP))
      k = 0
      Do n = 1, nEo
         Do i = 1, 4
            k = k + 1
            iW(iIPE + k - 1) = IPE(i, n)
         End do
      End do


c ... splitting elements
      Do n = 1, nEo
c        count new points
         newP = 0
         Do i = 1, 6
            iRt = iW(iIRE + 6 * (n-1) + i - 1)
            iRc = iW(ilbR + iRt - 1)
            If (iRc.GT.0) Then
                newP = newP + 1
                medge(newP) = iRc ! index of mother edge
            End if
         End do

         Do i = 1, 4
            IPEs(i) = IPE(i, n)
         End do
         IPEs(5) = lbE(n)
         IPEs(6) = Ref2MapMtr(n)

         If(newP.EQ.1) Then
c           tet is split into 2 tets,n and nE+1. iP1,iP2 are endpoints of mother edge
            nE = nE + 1
            If(nE.GT.MaxE) Call errMes(1006, 'localRefinement',
     &                                'local parameter MaxE is small')
            Do k = 1, 4
               IPE(k, nE) = IPEs(k)
            End do

            lbE(n)  = IPEs(5)
            lbE(nE) = IPEs(5)
            Ref2MapMtr(n)  = -iabs( IPEs(6) )
            Ref2MapMtr(nE) = -iabs( IPEs(6) )
            iP1 = iW(iendsR + (medge(1) - 1)*2)
            iP2 = iW(iendsR + (medge(1) - 1)*2+1)

            Do i = 1, 4
               If (iP2.eq.IPE(i, nE)) IPE(i, nE) = nPo + medge(1)
               If (iP1.eq.IPE(i, n )) IPE(i, n ) = nPo + medge(1)
            End do

         Else If(newP.EQ.2) Then
c           tet is split into 4 tets, n, nE+1,nE+2,nE+3. 
            nE = nE + 3
            If(nE.GT.MaxE) Call errMes(1006, 'localRefinement',
     &                                'local parameter MaxE is small')
            Do k = 1, 4
               IPE(k, nE-2) = IPEs(k)
               IPE(k, nE-1) = IPEs(k)
               IPE(k, nE  ) = IPEs(k)
            End do
            lbE(n)    = IPEs(5)
            lbE(nE-2) = IPEs(5)
            lbE(nE-1) = IPEs(5)
            lbE(nE  ) = IPEs(5)
            Ref2MapMtr(n)    = -iabs( IPEs(6) )
            Ref2MapMtr(nE-2) = -iabs( IPEs(6) )
            Ref2MapMtr(nE-1) = -iabs( IPEs(6) )
            Ref2MapMtr(nE  ) = -iabs( IPEs(6) )
c           iP1,iP2,IP3,iP4 are endpoints of mother edge
            iP1 = iW(iendsR + (medge(1) - 1)*2)
            iP2 = iW(iendsR + (medge(1) - 1)*2+1)
            iP3 = iW(iendsR + (medge(2) - 1)*2)
            iP4 = iW(iendsR + (medge(2) - 1)*2+1)

            Do i = 1, 4
               If (iP2.eq.IPE(i, n ))   IPE(i, n   ) = nPo + medge(1)
               If (iP4.eq.IPE(i, n ))   IPE(i, n   ) = nPo + medge(2)
               If (iP2.eq.IPE(i, nE-2)) IPE(i, nE-2) = nPo + medge(1)
               If (iP3.eq.IPE(i, nE-2)) IPE(i, nE-2) = nPo + medge(2)
               If (iP1.eq.IPE(i, nE-1)) IPE(i, nE-1) = nPo + medge(1)
               If (iP4.eq.IPE(i, nE-1)) IPE(i, nE-1) = nPo + medge(2)
               If (iP1.eq.IPE(i, nE  )) IPE(i, nE  ) = nPo + medge(1)
               If (iP3.eq.IPE(i, nE  )) IPE(i, nE  ) = nPo + medge(2)
            End do

         Else If(newP.EQ.3) Then
c           tet is split into 4 tets, n, nE+1,nE+2,nE+3.
            nE = nE + 3
            If(nE.GT.MaxE) Call errMes(1006, 'localRefinement',
     &                                'local parameter MaxE is small')
            Do i = 1, 3
               iP1 = iW(iendsR + (medge(i)-1)*2)
               iP2 = iW(iendsR + (medge(i)-1)*2+1)
               Do k = 1, 4
                  If(iP1.eq.IPEs(k).or.iP2.eq.IPEs(k)) IPEs(k)=0
               End do
               jP1 = iW(iendsR + (medge(loop(i+1))-1)*2)
               jP2 = iW(iendsR + (medge(loop(i+1))-1)*2+1)
               If (iP1.eq.jP1.or.iP1.eq.jP2) opposite(i) = iP1
               If (iP2.eq.jP1.or.iP2.eq.jP2) opposite(i) = iP2
            End do

            Do k = 1, 4
               If (IPEs(k).ne.0) iP0 = IPEs(k)
            End do

            IPE(1, n) = iP0
            IPE(2, n) = nPo + medge(1)
            IPE(3, n) = nPo + medge(2)
            IPE(4, n) = nPo + medge(3)
            lbE(n)    = IPEs(5)
            Ref2MapMtr(n) = -iabs( IPEs(6) )
            Do k = 1, 3
               IPE(1, nE-k+1) = iP0
               IPE(2, nE-k+1) = nPo + medge(k)
               IPE(3, nE-k+1) = nPo + medge(loop(k+1))
               IPE(4, nE-k+1) = opposite(k)
               lbE(nE-k+1) = IPEs(5)
               Ref2MapMtr(nE-k+1)    = -iabs( IPEs(6) )
            End do

         Else If(newP.EQ.4) Then
              Call errMes(1004, 'localRefinement',
     &                          'four new points in tet')
         Else If(newP.EQ.5) Then
              Call errMes(1005, 'localRefinement',
     &                          'five new points in tet')
         Else If(newP.EQ.6) Then
c ... find optimal splitting of inner octahedron (ii=1,2,3)
            Do ii = 1, 3
             minQt(ii)=1d0
             Do i = 1, 4
               Do k = 1, 4
                  iRt = iW(iIRE + 6 * n + refOct(k,i,ii) - 7)
                  iRc = iW(ilbR + iRt - 1)
                  j = nPo + iRc
                  jj= Ref2MapMtr(n)
                  if (jj.lt.0) jj = -jj
                  call FixedMap(XYP(1, j), MapMtr(1, jj), XYPt(1, k))
               End do
               call calREG(XYPt(1,1),XYPt(1,2),XYPt(1,3),XYPt(1,4), Qt)
               minQt(ii) = min( Qt, minQt(ii) )
             End do
            End do
            If (minQt(1).ge.minQt(2).and.minQt(1).ge.minQt(3)) Then
                ii = 1
            Else if (minQt(2).ge.minQt(1).and.minQt(2).ge.minQt(3)) Then
                ii = 2
            Else if (minQt(3).ge.minQt(1).and.minQt(3).ge.minQt(2)) Then
                ii = 3
            End if
c ... apply optimal splitting
            Do i = 1, 4
               nE = nE + 1
c  ...  interior tet
               If(nE.GT.MaxE) Call errMes(1006, 'localRefinement',
     &                                'local parameter MaxE is small')
               Do k = 1, 4
                  iRt = iW(iIRE + 6 * n + refOct(k,i,ii) - 7)
                  iRc = iW(ilbR + iRt - 1)
                  IPE(k, nE) = nPo + iRc
               End do
               lbE(nE) = IPEs(5)
               Ref2MapMtr(nE) = IPEs(6) 

c  ...  vertex-based tet (i<4) or original tet (i=4)
               If(i.EQ.4) Then
                  iE = n
               Else
                  nE = nE + 1
                  If(nE.GT.MaxE) Call errMes(1006, 'localRefinement',
     &                                'local parameter MaxE is small')
                  iE = nE
               End if

               Do k = 1, 3
                  iRt = iW(iIRE + 6 * n + refCrnr(k, i) - 7)
                  iRc = iW(ilbR + iRt - 1)
                  IPE(k, iE) = nPo + iRc
               End do
               IPE(4, iE) = IPEs(i)
               lbE(iE) = IPEs(5)
               Ref2MapMtr(iE) = IPEs(6) 
            End do
         End if
      End do


c ... splitting faces 
      nF = 0
      Do n = 1, nFo
         iP1 = IPF(1, n)
         iP2 = IPF(2, n)
         iP3 = IPF(3, n)

         flagE = cmpE(iP1, iP2, iP3, iW(iIEP), iW(inEP), 0, iE)

c        count new points at face and fill medge
         newP = 0
         Do i = 1, 6
            iRt = iW(iIRE + 6 * (iE-1) + i - 1)
            iRc = iW(ilbR + iRt - 1)
            If (iRc.GT.0) Then
                jP1 = iW(iendsR + (iRc-1)*2)
                jP2 = iW(iendsR + (iRc-1)*2+1)
                If (check13(jP1,iP1,iP2,iP3).and.
     &              check13(jP2,iP1,iP2,iP3)) Then
                    newP = newP + 1
                    medge(newP) = iRc
                End if
            End if
         End do
         
         If(newP.EQ.0) Then
            nF = nF + 1
            If(nF.GT.MaxF) Call errMes(1007, 'localRefinement',
     &                                'local parameter MaxF is small')
            iW(iIPF + (nF-1)*4  ) = iP1
            iW(iIPF + (nF-1)*4+1) = iP2
            iW(iIPF + (nF-1)*4+2) = iP3
            iW(iIPF + (nF-1)*4+3) = lbF(n)

         Else If(newP.EQ.1) Then
            nF = nF + 2
            If(nF.GT.MaxF) Call errMes(1007, 'localRefinement',
     &                                'local parameter MaxF is small')
            jP1 = iW(iendsR + (medge(1)-1)*2)
            jP2 = iW(iendsR + (medge(1)-1)*2+1)
            If(check22(iP1, iP2, jP1, jP2)) Then
               iPa = iP3
               iPb = iP1
               iPc = iP2
            End if
            If(check22(iP1, iP3, jP1, jP2)) Then
               iPa = iP2
               iPb = iP1
               iPc = iP3
            End if
            If(check22(iP2, iP3, jP1, jP2)) Then
               iPa = iP1
               iPb = iP2
               iPc = iP3
            End if
            iW(iIPF + (nF-2)*4  ) = iPa
            iW(iIPF + (nF-2)*4+1) = iPb
            iW(iIPF + (nF-2)*4+2) = nPo + medge(1)
            iW(iIPF + (nF-2)*4+3) = lbF(n) 
            iW(iIPF + (nF-1)*4  ) = iPa
            iW(iIPF + (nF-1)*4+1) = iPc
            iW(iIPF + (nF-1)*4+2) = nPo + medge(1)
            iW(iIPF + (nF-1)*4+3) = lbF(n) 

         Else If(newP.EQ.3) Then
c ... define opposite
            Do i = 1, 3
               iP1 = iW(iendsR + (medge(i)-1)*2)
               iP2 = iW(iendsR + (medge(i)-1)*2+1)
               jP1 = iW(iendsR + (medge(loop(i+1))-1)*2)
               jP2 = iW(iendsR + (medge(loop(i+1))-1)*2+1)
               If (iP1.eq.jP1.or.iP1.eq.jP2) opposite(i) = iP1
               If (iP2.eq.jP1.or.iP2.eq.jP2) opposite(i) = iP2
            End do
c ...
            nF = nF + 4
            If(nF.GT.MaxF) Call errMes(1007, 'localRefinement',
     &                                'local parameter MaxF is small')
c           interior triangle
            iPa = nPo + medge(1)
            iPb = nPo + medge(2)
            iPc = nPo + medge(3)
            iW(iIPF + (nF-4)*4  ) = iPa
            iW(iIPF + (nF-4)*4+1) = iPb
            iW(iIPF + (nF-4)*4+2) = iPc
            iW(iIPF + (nF-4)*4+3) = lbF(n) 
            Do k = 1, 3 ! corner triangles
               iW(iIPF + (nF-k)*4  ) = opposite(k)
               iW(iIPF + (nF-k)*4+1) = nPo + medge(k)
               iW(iIPF + (nF-k)*4+2) = nPo + medge(loop(k+1))
               iW(iIPF + (nF-k)*4+3) = lbF(n) 
            End do
         Else
            Call errMes(1008, 'localRefinement',
     &                  'newP for face  is wrong')
         End if 
      End do

      Do n = 1, nF
         IPF(1,n) = iW(iIPF + (n-1)*4  )
         IPF(2,n) = iW(iIPF + (n-1)*4+1)
         IPF(3,n) = iW(iIPF + (n-1)*4+2)
         lbF(  n) = iW(iIPF + (n-1)*4+3)
      End do

      Return
      End Subroutine localRefinement


C ================================================================
      Subroutine calREG(xy1, xy2, xy3, xy4, qE)
C ================================================================
C Computing regularity quality of tetrahedron {xy1, ..., xy4} 
C ================================================================
      Real*8 xy1(3), xy2(3), xy3(3), xy4(3)
      Real*8 qE

C group (Local variables)
      Real*8 Pk, Vk, volume
      Real*8 x1, y1, z1


C ================================================================
      Pk = 0D0
      Do n = 1, 6
         If(n.EQ.1) Then
            x1 = xy1(1) - xy4(1)
            y1 = xy1(2) - xy4(2)
            z1 = xy1(3) - xy4(3)
         Else If(n.EQ.2) Then
            x1 = xy2(1) - xy4(1)
            y1 = xy2(2) - xy4(2)
            z1 = xy2(3) - xy4(3)
         Else If(n.EQ.3) Then
            x1 = xy3(1) - xy4(1)
            y1 = xy3(2) - xy4(2)
            z1 = xy3(3) - xy4(3)
         Else If(n.EQ.4) Then
            x1 = xy1(1) - xy3(1)
            y1 = xy1(2) - xy3(2)
            z1 = xy1(3) - xy3(3)
         Else If(n.EQ.5) Then
            x1 = xy2(1) - xy3(1)
            y1 = xy2(2) - xy3(2)
            z1 = xy2(3) - xy3(3)
         Else If(n.EQ.6) Then
            x1 = xy1(1) - xy2(1)
            y1 = xy1(2) - xy2(2)
            z1 = xy1(3) - xy2(3)
         End if
         Pk = Pk + dsqrt(x1 * x1 +  y1 * y1 + z1 * z1)
      End do

      volume = calVol(xy1, xy2, xy3, xy4)
      Vk = dabs(volume) 

      qE = 1832.8208D0 * Vk / (Pk ** 3) 

      Return
      End Subroutine calREG

      Subroutine FixedMap(XYP, MapMatrix, XYPt)
      Implicit none
      Real*8   XYP(3), MapMatrix(3,3), XYPt(3)

      Integer i,j

      do i=1,3
       XYPt(i) = 0d0
       do j=1,3
          XYPt(i) = XYPt(i) + MapMatrix(i,j)*XYP(j)
       end do
      end do

      return
      end Subroutine FixedMap
C
      End Module mba3d_refine
