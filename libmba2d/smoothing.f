C ================================================================
      Subroutine smoothingS(
C ================================================================
C Routine smoothes the piecewise linear function Sol defined at 
C mesh points. The value of the smoothed function at a mesh point 
C P is equal to the integral average over the superelement 
C associated with point P. This method is known as the Steklov 
c smoothing.
C ================================================================
     &           nP, nE, XYP, IPE, Sol,
     &           MaxWr, MaxWi, rW, iW)
C ================================================================
      real  XYP(2, *), Sol(*)
      Integer IPE(3, *)

      Integer MaxWi, MaxWr
      Integer iW(*)
      real  rW(*)

c group (Local variables)
      real  v, calVol
      Character*80 message

C ================================================================
      inEP = 0
      iIEP = inEP + nP
      iEnd = iIEP + 3 * nE

      If(iEnd.GT.MaxWi) Then
         iERR = 1001
         Write(message,'(A,I10)')
     &        'The approximate size of iW is ', iEnd
         Call errMes(iERR, 'smoothingS', message)
      End if

      iVol  = 0
      iSup = iVol + nE
      iSol = iSup + nP
      iEnd  = iSol + nP

      If(iEnd.GT.MaxWr) Then
         iERR = 1002
         Write(message,'(A,I10)')
     &        'The approximate size of rW is ', iEnd
         Call errMes(iERR, 'smoothingS', message)
      End if


C ... creating an auxiliary structure
      Call backReferences(nP, nE, 3, 3, IPE, iW(inEP + 1), iW(iIEP + 1))


C ... computing volumes of elemnts and superelements
      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         v = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))

         rW(iVol + n) = abs(v)
      End do

      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = iW(inEP + n)

         v = 0D0
         Do i = i1, i2
            iE = iW(iIEP + i)
            v = v + rW(iVol + iE)
         End do
         rW(iSup + n) = v
      End do
      

C ... itegrating the piecewise linear function SOL
      Do n = 1, nP
         rW(iSol + n) = 0D0
      End do

      Do n = 1, nE
         v = 0D0
         Do i = 1, 3
            iP1 = IPE(i, n)
            v = v + Sol(iP1) 
         End do
         v = v * rW(iVol + n) / 3 

         Do i = 1, 3
            iP1 = IPE(i, n)
            rW(iSol + iP1) = rW(iSol + iP1) + v  / rW(iSup + iP1)
         End do
      End do
      
      Do n = 1, nP
         Sol(n) = rW(iSol + n)
      End do

      Return
      End


