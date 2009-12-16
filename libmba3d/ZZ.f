      Module mba3d_ZZ
C
      use mba3d_error
      use mba3d_findSE
      use mba3d_makQ
      use mba3d_utils
C
      contains
C
C ================================================================
      Subroutine P02P1(
C ================================================================
     &           nP, nF, nR, nE, XYP, IPF, IPE, 
     &           fP0, fP1, 
     &           MaxWr, MaxWi, rW, iW)
C ================================================================
      implicit integer (i)
      include 'makS.fd'
C ================================================================
C  The routine maps a discontinuous piece-wise constant function 
C  with d.o.f in elements onto a continuous piece-wise linear 
C  function with d.o.f. at vertices. We use the ZZ method for that.
C
C  Limitation to the mesh: Each boundary node can be connected with
C                          an interior note by at most 2 mesh edges.
C
C  *** Remarks 
C         1. Extrapolation at boundary nodes adds some diffusion
C            for functions with sharp gradients near boundary.
C ================================================================
C  *** Input
C         nP  - the number of mesh nodes
C         nF  - the number of boundary edges
C         nR  - the number of mesh edges
C         nE  - the number of triangles
C
C         XYP - coordinates of the mesh nodes
C         IPF - map: boundary edge -> vertices 
C         IPE - map: triangle -> vertices
C
C         fP0 - discontinuous piece-wise constant function
C
C  *** Output: 
C         fP1 - continuous piece-wice linear function
C                    
C  *** Work memory: 
C          rW - real*8  array of size MaxWr
C          iW - integer array of size MaxWi
C ================================================================
      Integer  nP, nF, nR, nE

      Real*8   XYP(3, *)
      Integer  IPF(3, *), IPE(4, *)

      Real*8   fP0(*), fP1(*)

      Integer  MaxWr, MaxWi
      Integer  iW(*)
      Real*8   rW(*)

C ... external functions

C ... local variables
      Real*8   xy(3, MaxS), values(MaxS), weights(MaxS)
      Integer  kbuf(MaxS)

      Real*8   vol, s

C ================================================================
C ... distribute memory
      iBnd = 0
      iArm = iBnd + nP
      inEP = iArm + nP
      iIEP = inEP + nP
      iiW  = iIEP + 4 * nE

      iMass = 0
      irW   = iMass + nP

      If(MaxWi.LT.iiW) 
     &   Call errMes(1001, 'P02P1', 'Increase size of MaxWi')

      If(MaxWr.LT.irW) 
     &   Call errMes(1002, 'P02P1', 'Increase size of MaxWr')
      

C ... create maps: vertix -> triangles
      Call backReferences(nP, nE, 4, 4, IPE, iW(inEP+1), iW(iIEP+1))


C ... mark boundary nodes
      Do n = 1, 2 * nP 
         iW(n) = 0
      End do

      Do n = 1, nF
         Do i = 1, 2
            iP = IPF(i, n)
            iW(iBnd + iP) = 1
         End do
      End do


C ... initialize Mass and fP1
      Do n = 1, nP
         rW(iMass + n) = 0D0
         fP1(n) = 0D0
      End do


C ... interpolate in interior points with least squares
      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = iW(inEP + n)

         If(i2 - i1 + 1.GT.MaxS) Call errMes(1007,
     &      'DG2P1', 'Local parameter MaxS is small')

         Do i = i1, i2
            iE = iW(iIEP + i)
   
            iP1 = IPE(1, iE)
            iP2 = IPE(2, iE)
            iP3 = IPE(3, iE)
            iP4 = IPE(4, iE)

            vol = calVol(XYP(1, iP1), XYP(1, iP2), 
     &                                XYP(1, iP3), XYP(1, iP4))
            vol = dabs(vol)

            If(iW(iBnd + n).GT.0) Then
               Do j = 1, 4
                  iP = IPE(j, iE)
                  If(iW(iBnd + iP).EQ.0) 
     &               rW(iMass + n) = rW(iMass + n) + vol / 20
               End do
            Else
               rW(iMass + n) = rW(iMass + n) + vol / 4

               k = i - i1 + 1
               values(k) = fP0(iE)
               Do j = 1, 3
                  xy(j, k) = (XYP(j,iP1) + XYP(j,iP2) 
     &                                   + XYP(j,iP3) + XYP(j,iP4)) / 4
               End do
            End if
         End do

         If(iW(iBnd + n).EQ.0) Then
            k = i2 - i1 + 1
            fP1(n) = LSvalue(k, xy, values, XYP(1, n))
         End if
      End do


C ... extrapolate fP1 at boundary nodes from nearest interior nodes
      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = iW(inEP + n)

         If(iW(iBnd+n).EQ.1) Then
            lbuf = 0

            Do i = i1, i2
               iE = iW(iIEP + i)
   
               iP1 = IPE(1, iE)
               iP2 = IPE(2, iE)
               iP3 = IPE(3, iE)
               iP4 = IPE(4, iE)

               vol = calVol(XYP(1, iP1), XYP(1, iP2), 
     &                                   XYP(1, iP3), XYP(i, iP4))
               vol = dabs(vol)

               Do j = 1, 4
                  iP = IPE(j, iE)

                  Call findSE(lbuf, kbuf, iP, k)

                  If(k.GT.0) Then
                     weights(k) = weights(k) + vol / 20
                  Else If(iW(iBnd + iP).EQ.0) Then
                     lbuf = lbuf + 1
                     kbuf(lbuf) = iP 
                     weights(lbuf) = vol / 20
                  End if
               End do
            End do

            Do k = 1, lbuf
               iP = kbuf(k)
               fP1(n) = fP1(n) + fP1(iP) * weights(k) / rW(iMass + n)
            End do

            If(lbuf.GT.0) iW(iArm + n) = 1
         End if 
      End do
     

c ... remove boundary node from the list
      Do n = 1, nP
         If(iW(iArm + n).EQ.1) iW(iBnd + n) = 0
      End do


c ... the second level of extrapolation
      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = iW(inEP + n)

         If(iW(iBnd + n).EQ.1) Then
            k = 0
            s = 0D0

            Do i = i1, i2
               iE = iW(iIEP + i)
               Do j = 1, 4
                  iP = IPE(j, iE)
                  If(iW(iBnd + iP).EQ.0) Then
                     k = k + 1
                     s = s + fP1(iP)  
                  End if
               End do
            End do

            If(k.EQ.0) 
     &         Call errMes(2011, 'P02P1', 'Mesh violates 2-arm rule')

            fP1(n) = s / k
         End if
      End do

      Return
      End Subroutine P02P1



C ================================================================
      Real*8 function LSvalue(k, xy, values, xy0)
C ================================================================
C  This routine uses least square linear fit to points xy(2,*) and
C  evaluates the value of the linear function at point xy0.
C ================================================================
      implicit none

      Integer  k
      Real*8   xy(3,*), values(*), xy0(3)

c (local variables)
      Real*8   A(4, 4), S(4), work(40)
      Integer  i, j, ipiv(4), info     
      
C ================================================================
      Do i = 1, 4 
         Do j = 1, 4 
            A(i, j) = 0D0
         End do
         S(i) = 0D0
      End do

c ... generate the least squares matrix
      Do i = 1, k
         A(1,1) = A(1,1) + xy(1, i) * xy(1, i)
         A(1,2) = A(1,2) + xy(1, i) * xy(2, i)
         A(1,3) = A(1,3) + xy(1, i) * xy(3, i)
         A(1,4) = A(1,4) + xy(1, i) 

         A(2,2) = A(2,2) + xy(2, i) * xy(2, i)
         A(2,3) = A(2,3) + xy(2, i) * xy(3, i)
         A(2,4) = A(2,4) + xy(2, i) 

         A(3,3) = A(3,3) + xy(3, i) * xy(3, i)
         A(3,4) = A(3,4) + xy(3, i) 
      End do
      A(4,4) = k

c ... generate the RHS
      do i = 1, k
         S(1) = S(1) + xy(1, i) * values(i)
         S(2) = S(2) + xy(2, i) * values(i)
         S(3) = S(3) + xy(3, i) * values(i)
         S(4) = S(4) + values(i)
      end do

      Call dsysv('U', 4, 1, A, 4, ipiv, S, 4, work, 40, info)

      If(info.NE.0) Call errMes(3011, 'LSvalue',
     &                   'Error in Lapack routine dsysv')

c ... evaluate the linear function
      LSvalue = S(1) * xy0(1) + S(2) * xy0(2) + S(3) * xy0(3) + S(4)

      Return
      End function LSvalue
C
      End Module mba3d_ZZ
