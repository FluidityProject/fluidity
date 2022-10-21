C ==========================================================
      Subroutine prjCrv(XY, prjXY, iFNC, t, calCrv,
     &                  L1Et, L2Et, nL2t, nStept, nEt, tE)
C ==========================================================
C Routine computes a point which lies on the boundary of the 
C input polygonal domain and is the closest point to XY. The 
C later is the point lying on new boundary. 
C ==========================================================
      real  XY(2), prjXY(2), t

      EXTERNAL calCrv

      real  tE(*)
      Integer L1Et(2, *), L2Et(*), nStept(4)

      real  XYprev(2), XYnext(2), tnext, tprev
      Integer iPrevL2t, PrevL2, iPrevL1t, PrevL1


      If(iFNC.EQ.0) Then
         prjXY(1) = XY(1)
         prjXY(2) = XY(2)
      Else
         iPrevL2t = PrevL2(nL2t, L2Et, tE, t)
         iPrevL1t = PrevL1(nEt, L1Et, L2Et, iPrevL2t, nStept(1), tE, t)
         tprev = tE(iPrevL1t)
         tnext = tE(L1Et(2, iPrevL1t))

         Call  aniCrv(tprev, XYprev, iFNC, calCrv)
         Call  aniCrv(tnext, XYnext, iFNC, calCrv)
         prjXY(1) = (XYprev(1) * (tnext - t) / (tnext - tprev) +
     &               XYnext(1) * (t - tprev) / (tnext - tprev))
         prjXY(2) = (XYprev(2) * (tnext - t) / (tnext - tprev) +
     &               XYnext(2) * (t - tprev) / (tnext - tprev))
      End if

      Return
      End



C ==========================================================
      Subroutine calCrvFnc(IPF, nF, iFnc, LFnc, nCrvFnc)
C ==========================================================
C Routine computes the number of different functions 
C describing curvilinear boundaries and internal interfaces.
C ==========================================================
      Integer IPF(4, *), iFnc(*), LFnc(*)

      nCrvFnc = 0
      Do n = 1, nF
         iCrv = IPF(3, n)
         If(iCrv.NE.0) Then
            it = iFnc(n)
            Call findSE(nCrvFnc, LFnc, it, k)
            If(k.EQ.0) Then
               nCrvFnc = nCrvFnc + 1
               LFnc(nCrvFnc) = it
            End if
         End if
      End do

      Return
      End


C ==========================================================
      Subroutine tEMak(tE, nEt, MaxF, IPF, nF, parCrv, iFnc,
     &                 iCrvFnc)
C ==========================================================
C Routine creates un-ordered list of parameters for boundary
C part described by iCrvFnc-th function.
C ==========================================================
      real   tE(*), parCrv(2, *)
      Integer  IPF(4, *), iFnc(*)

      real   tmin
      Logical  flag

c ... compute minimal parameter for boundary iFnc
      flag = .TRUE.
      Do i = 1, nF
         If(IPF(3, i).NE.0) Then
            If(iFnc(i).EQ.iCrvFnc) Then
               Do k = 1, 2
                  If(flag) Then
                     tmin = parCrv(k, i)
                     flag = .FALSE.
                  Else
                     tmin = min(tmin, parCrv(k, i))
                  End if
               End do
            End if
         End if
      End do


c ... collect maxima of two parameters in tE 
      ic = 1
      tE(ic) = tmin  
      Do i = 1, nF
         If(IPF(3, i).NE.0) Then
            If(iFnc(i).EQ.iCrvFnc) Then
               ic = ic + 1
c  ...  this should never happen
               If(ic.GT.MaxF) Call errMes(6003, 
     &                             'tEMak', 'System error')
               If(parCrv(1, i).GT.parCrv(2, i)) Then
                  tE(ic) = parCrv(1, i)
               Else
                  tE(ic) = parCrv(2, i)
               End if
            End if
         End if
      End do

      nEt = ic
      Return
      End
