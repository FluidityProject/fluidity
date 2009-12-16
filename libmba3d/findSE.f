      Module mba3d_findSE
C
      contains
C
C ================================================================
C     Subroutine findSP(lP, iPs, iP, nPs)
C     Subroutine findSF(lF, iFs, iF, nFs)
      Subroutine findSE(lE, iEs, iE, nEs)
C ================================================================
C Search for index i such that iEs(i) = iE.
C Zero is returned when the search fails.
C ================================================================
      Integer iEs(*)

      nEs = 0
      Do n = 1, lE
         If(iEs(n).EQ.iE) Then
            nEs = n
            Goto 1000
         End if
      End do
 1000 Return
      End Subroutine findSE
C
      End module mba3d_findSE
