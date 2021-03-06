C ================================================================
C The routines below work with an ordered list of reals.
C The available operations are:
C     (a) intialize the list
C     (b) add new element in the list
C     (c) update an element value
C     (d) exclude an element from the list
C     (e) careful check the list data
C
C *** Remarsk:
C         1. This is the new version of list routines. It requires 
C            double memory for L2E (compared to the old version), 
C            and use iCntl(4) instead of nstep. Besides, it uses 
C            iCntl as the input data to be defined outside:
C
C            iCntl(1) is typical interval length, we recommend sqrt(nEStar)
C            iCntl(2) is rank of interval length belonging to 
C                     [iCntl(1)-iCntl(2), iCntl(1)+iCntl(2)]
C            iCntl(3) is the pointer to the middle of L2E
C            iCntl(4) is the output channel in case of debugging 
C                     (=0 for no debugging)
C
C ================================================================
      Subroutine lstMak(nEStar, nE, L1E, L2E, nL2, iCntl, LHol)
C ================================================================
      Implicit none
      Integer L1E(2,*), L2E(*), LHol(*), nEStar, nE, nL2, iCntl(4)

      Integer n, nstep, ipos2

C ================================================================
      nstep = iCntl(1)
      ipos2 = iCntl(3)

      L1E(1,L2E(1)) = 0
      L1E(2,L2E(nE)) = 0

      Do n = 1, nE - 1
         L1E(2,L2E(n)) = L2E(n + 1)
      End do

      Do n = 2, nE
         L1E(1,L2E(n)) = L2E(n - 1)
      End do

      nL2 = 0
      Do n = 1, nE, nstep
         nL2 = nL2 + 1
         L2E(nL2) = L2E(n)
         L2E(ipos2+nL2) = nstep
      End do
      L2E(ipos2+nL2) = nE - nstep*(nL2-1)


      LHol(1) = 0

      Return
      End



C ================================================================
      Subroutine lstUpd(nE, L1E, nL2, L2E, iCntl, qE, iE,qiE)
C ================================================================
C     Input:  iE,qiE
C     Output: Updated qE and L1E,L2E
C ================================================================
      Implicit none
      Integer           L1E(2,*), L2E(*), nE, nL2, iE, iCntl(4)
      real  qE(*),qiE
      Integer           PrevL2, PrevL1, PrevL2iE

      Integer           nstep,nthr,ipos2
      Integer           i,ia,iPrevL1,iPrevL2,iNextL1,iLast,k,ib

      Logical           flagFP

C ================================================================
      nstep = iCntl(1)
      nthr  = iCntl(2)
      ipos2 = iCntl(3)
C     Chan  = iCntl(4)

      If (L1E(1,iE).eq.0.and.L1E(2,iE).eq.0) then
         Call errMes(5101, 'lstUpd', 'no element in the list')
      End if

      Call fpcheck(qiE, flagFP)

      If(flagFP) Then
         Call errMes(5111, 'lstUpd', 'bad input value(NAN/INF)')
      End if


      iPrevL2 = PrevL2iE( nL2, L2E, qE, iE, L1E )
         

      If (L1E(1,iE).ne.0) then
         L1E(2,L1E(1,iE)) = L1E(2,iE)
      Else
         L1E(1,L1E(2,iE)) = 0
      End if
      If (L1E(2,iE).ne.0) then
         L1E(1,L1E(2,iE)) = L1E(1,iE)
      Else
         L1E(2,L1E(1,iE)) = 0
      End if
      nE = nE - 1


      If (iPrevL2.eq.0) then
       ia = 1
       ib = 1
      Else
       ib = iPrevL2
       If (L2E(iPrevL2).eq.iE) then
         ia = 0
       Else
         ia = 1
       End if
       if (ib.lt.nL2) then
           if (iE.eq.L2E(ib+1)) ib=ib+1
       end if
      End if
      if (L2E(ipos2+ib)-1.ge.nstep-nthr) then
          L2E(ipos2+ib) = L2E(ipos2+ib)-1
          If (L2E(ib).eq.iE) L2E(ib) = L1E(2,L2E(ib))
      else
          Do i = iPrevL2+ia, nL2
             L2E(i) = L1E(2,L2E(i))
          End do
c  ...    the last one may be just positive !!!
          if (L2E(ipos2+nL2) - 1.ge.1) then 
              L2E(ipos2+nL2) = L2E(ipos2+nL2) - 1
          else
              If (L2E(nL2).ne.0) then
                 Call errMes(5102, 'lstUpd', 'L2E(nL2) must be 0 here')
              end if
              nL2 = nL2 - 1
          end if
      end if

      L1E(1,iE) = 0
      L1E(2,iE) = 0


c ... add new element
      iPrevL2 = PrevL2( nL2, L2E, qE, qiE )
      if (iPrevL2.eq.0) then
       iPrevL1 = 0
      else
       iPrevL1 = PrevL1( nE, L1E, L2E, iPrevL2, nstep, qE, qiE )
      end if


      nE = nE + 1
      qE(iE) = qiE

      If(iPrevL1.ne.0) then
         iNextL1 = L1E(2,iPrevL1)
      Else
         iNextL1 = L2E(1)
      End if
      If (iNextL1.ne.0) L1E(1,iNextL1) = iE
      If (iPrevL1.ne.0) L1E(2,iPrevL1) = iE
      L1E(1,iE)      = iPrevL1
      L1E(2,iE)      = iNextL1


      ib = iPrevL2
      if (iPrevL2.eq.0) ib=1

      if (L2E(ipos2+ib)+1.le.nstep+nthr) then
          L2E(ipos2+ib) = L2E(ipos2+ib) + 1
          If (iPrevL2.eq.0) L2E(ib) = L1E(1,L2E(ib))
      else
          Do i = iPrevL2+1, nL2
             L2E(i) = L1E(1,L2E(i))
          End do
          if (L2E(ipos2+nL2) + 1.le.nstep+nthr) then
              L2E(ipos2+nL2) = L2E(ipos2+nL2) + 1
          else
              iLast = L2E(nL2)
              Do while (.true.)
                 k = L1E(2,iLast)
                 If (k.eq.0) then
                    nL2 = nL2 + 1
                    L2E(nL2) = iLast
                    L2E(ipos2+nL2) = 1
                    goto 1
                 End If
                 iLast = k
              End do
          end if
      End if

 1    Continue

      Return  
      End



C ================================================================
      Subroutine lstDel(nE, L1E, nL2, L2E, iCntl, LHol, qE, iE)
C ================================================================
C     Input:  iE
C     Output: Updated L1E,L2E,nE,nL2,LHol
C ================================================================
      Implicit none
      Integer            L1E(2,*),L2E(*),LHol(*),nE,nL2,iE,iCntl(4)
      real   qE(*)
      Integer            PrevL2iE

      Integer            nstep,nthr,ipos2
      Integer            ia,i,iPrevL2,ib

C ================================================================
      nstep = iCntl(1)
      nthr  = iCntl(2)
      ipos2 = iCntl(3)
C     Chan  = iCntl(4)

      If (L1E(1,iE).eq.0.and.L1E(2,iE).eq.0) then
         Call errMes(5103, 'lstDel', 'no element in the list')
      End if

      iPrevL2 = PrevL2iE( nL2, L2E, qE, iE, L1E )

      If (L1E(1,iE).ne.0) then
         L1E(2,L1E(1,iE)) = L1E(2,iE)
      Else
         L1E(1,L1E(2,iE)) = 0
      End if
      If (L1E(2,iE).ne.0) then
         L1E(1,L1E(2,iE)) = L1E(1,iE)
      Else
         L1E(2,L1E(1,iE)) = 0
      End if
      nE = nE - 1

      If (iPrevL2.eq.0) then
       ia = 1
       ib = 1
      Else
       ib = iPrevL2
       If (L2E(iPrevL2).eq.iE) then
         ia = 0
       Else
         ia = 1
       End if
       if (ib.lt.nL2) then
           if (iE.eq.L2E(ib+1)) ib=ib+1
       end if
      End if
      if (L2E(ipos2+ib)-1.ge.nstep-nthr) then
          L2E(ipos2+ib) = L2E(ipos2+ib)-1
          If (L2E(ib).eq.iE) L2E(ib) = L1E(2,L2E(ib))
      else
          Do i = iPrevL2+ia, nL2
             L2E(i) = L1E(2,L2E(i))
          End do
c  ...   the last one may be just positive !!!
          if (L2E(ipos2+nL2) - 1.ge.1) then 
              L2E(ipos2+nL2) = L2E(ipos2+nL2) - 1
          else
              If (L2E(nL2).ne.0) then
                 Call errMes(5104, 'lstDel', 'L2E(nL2) must be 0 here')
              End if
              nL2 = nL2 - 1
          end if
      End if

      L1E(1,iE) = 0
      L1E(2,iE) = 0
      LHol(1) = LHol(1) + 1
      LHol(LHol(1)+1) = iE

      Return  
      End



C ================================================================
      Subroutine lstAdd(nE, L1E, nL2, L2E, iCntl, LHol, qE, qiE, iE)
C ================================================================
C     Input:  qiE,LHol
C     Output: iE, updated L1E,L2E,nE,nL2,LHol,qe
C ================================================================
      Implicit none
      Integer           L1E(2,*), L2E(*), LHol(*),nE,nL2,iE,iCntl(4)
      real  qE(*),qiE
      Integer           PrevL2, PrevL1

      Integer           nstep,nthr,ipos2
      Integer           iPrevL1,iNextL1,iPrevL2,iLast,i,k,ib

      Logical           flagFP
  
C ================================================================
      nstep = iCntl(1)
      nthr  = iCntl(2)
      ipos2 = iCntl(3)
C     Chan  = iCntl(4)

      Call fpcheck(qiE, flagFP)

      If(flagFP) Then
         Call errMes(5113, 'lstAdd', 'bad input value(NAN/INF)')
      End if


      iPrevL2 = PrevL2( nL2, L2E, qE, qiE )
      if (iPrevL2.eq.0) then
       iPrevL1 = 0
      else
       iPrevL1 = PrevL1( nE, L1E, L2E, iPrevL2, nstep, qE, qiE )
      end if

      nE = nE + 1
      If (LHol(1).eq.0) then
         iE = nE
      Else
         iE = LHol(LHol(1)+1)
         LHol(1) = LHol(1) - 1
      End if
      qE(iE) = qiE

      If(iPrevL1.ne.0) then
         iNextL1 = L1E(2,iPrevL1)
      Else
         iNextL1 = L2E(1)
      End if
      If (iNextL1.ne.0) L1E(1,iNextL1) = iE
      If (iPrevL1.ne.0) L1E(2,iPrevL1) = iE
      L1E(1,iE)      = iPrevL1
      L1E(2,iE)      = iNextL1


      ib = iPrevL2
      if (iPrevL2.eq.0) ib=1

      if (L2E(ipos2+ib)+1.le.nstep+nthr) then
          L2E(ipos2+ib) = L2E(ipos2+ib) + 1
          If (iPrevL2.eq.0) L2E(ib) = L1E(1,L2E(ib))
      else
          Do i = iPrevL2+1, nL2
             L2E(i) = L1E(1,L2E(i))
          End do
          if (L2E(ipos2+nL2) + 1.le.nstep+nthr) then
              L2E(ipos2+nL2) = L2E(ipos2+nL2) + 1
          else
              iLast = L2E(nL2)
              Do while (.true.)
                 k = L1E(2,iLast)
                 If (k.eq.0) then
                    nL2 = nL2 + 1
                    L2E(nL2) = iLast
                    L2E(ipos2+nL2) = 1
                    goto 1
                 End If
                 iLast = k
              End do
          end if
      End if

 1    Continue

      Return
      End


C ================================================================
      Integer Function PrevL2(nL2, L2E, qE, qiE)
C ================================================================
      Implicit none
      Integer           L2E(*), nL2
      real  qE(*),qiE

      Integer           i1,i2,i3  

C ================================================================
      If (qiE.lt.qE(L2E(1))) then
          PrevL2 = 0
          Return
      End If
      If (qiE.ge.qE(L2E(nL2))) then
          PrevL2 = nL2
          Return
      End If
      i1 = 1
      i2 = nL2-1
      Do while (.true.)
        i3 = (i1+i2)/2
        If (i1.eq.i2) then
          PrevL2 = i1 
          Return
        End if
        If (i1.eq.i2-1) then
          If (qiE.lt.qE(L2E(i2))) then
              PrevL2 = i1 
          Else
              PrevL2 = i2 
          End if
          Return
        End if
        If (qiE.lt.qE(L2E(i3))) then
            i2 = i3
        Else
            i1 = i3
        End If
      End do

      Return
      End



C ================================================================
      Integer Function PrevL2iE(nL2, L2E, qE, iE, L1E)
C ================================================================
      Implicit none
      Integer L2E(*), nL2, L1E(2,*), iE
      real   qE(*)

      Integer i, ii,j,k, i1,i2,i3
      real   qiE

C ================================================================
      If (L2E(1).eq.iE) then
          PrevL2iE = 0
          Return
      End if
      qiE = qE(iE)


      If (qiE.lt.qE(L2E(1))) then
          goto 55
      End If
      If (qiE.ge.qE(L2E(nL2))) then
          goto 55
      End If
      i1 = 1
      i2 = nL2-1
      Do while (.true.)
        i3 = (i1+i2)/2
        If (i1.eq.i2) then
          i = i1 +1
          goto 22
        End if
        If (i1.eq.i2-1) then
          If (qiE.lt.qE(L2E(i2))) then
              i = i1 +1 
          Else
              i = i2 +1
          End if
          goto 22
        End if
        If (qiE.lt.qE(L2E(i3))) then
            i2 = i3
        Else
            i1 = i3
        End If
      End do

22    continue
      Do ii = i, 2, -1
       j = L2E(ii)
       Do while (.true.)
          if (j.eq.iE) goto 33 
          k = L1E(1,j)
          if (k.eq.L2E(ii-1)) goto 44
          j = k
       End do
33     PrevL2iE = ii - 1
       Return
44     continue
      End do


55    PrevL2iE = nL2

      If (qiE.eq.qE(L2E(nL2))) then
           Do ii = nL2, 2, -1
            j = L2E(ii)
            Do while (.true.)
               if (j.eq.iE) goto 888 
               k = L1E(1,j)
               if (k.eq.L2E(ii-1)) goto 999
               j = k
            End do
888         PrevL2iE = ii - 1
            Return
999         continue
           End do
      End if

      Return
      End

  

C ================================================================
      Integer Function PrevL1(nE, L1E, L2E, iPrevL2, nstep, qE, qiE)
C ================================================================
      Implicit none
      Integer L1E(2,*), L2E(*), nE, iPrevL2, nstep
      real   qE(*),qiE

      Integer j,k

C ================================================================
      if (iPrevL2.le.0) then
         Call errMes(5005, 'PrevL1', 'wrong iPrevL2')
      end if
      j = L2E(iPrevL2)
      Do while (.true.)
         If (qiE.lt.qE(j)) then
            PrevL1 = L1E(1,j)
            Return
         End If
         k = L1E(2,j)
         If (k.eq.0) then
            PrevL1 = j
            Return
         End If
         j = k
      End do

      Return
      End



