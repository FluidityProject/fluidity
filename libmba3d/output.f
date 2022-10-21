      Module mba3d_output
C
      contains
C
C ================================================================
      Subroutine countBadElements(nE, L1E, L2E, qE, Quality, 
     &                            nLines, output, flagFILE)
C ================================================================
      include 'output.fd'
C ================================================================
      Integer L1E(2, *), L2E(*)
      Real*8  qE(*), Quality

      Character*(LineLenght) output(*)
      Logical flagFILE

c group (Local variables)
      Integer nBad(10) 
      Character*(LineLenght) message

C ================================================================
      Do n = 1, 10
         nBad(n) = 0
      End do

      icnt = 0
      iE = L2E(1)
      Do n = 1, nE
         If(qE(iE).LT.Quality) icnt = icnt + 1
         Do k = 1, 10
            If((k - 1) * 1D-1.LE.qE(iE) .AND.
     &           qE(iE).LT.k * 1D-1) Then
               nBad(k) = nBad(k) + 1
               Goto 500
            End if
         End do
 500     iE = L1E(2, iE)
      End do


      If(flagFILE) Then
         Write(message,6005)
         Call addOut(nLines, message, output)
         Write(message,6006) (1D-1 * i, i = 1, 10)
         Call addOut(nLines, message, output)
         Write(message,6008) (nBad(i), i = 1, 10)
         Call addOut(nLines, message, output)
         Write(message,6007) icnt
         Call addOut(nLines, message, output)
      Else
         Write(*, 5005) (1D-1 * i, i = 1, 10),
     &                  (nBad(i), i = 1, 10), icnt
      End if


c ... screen output
 5005 Format(
     & ' *** Distribution of tetrahedrons by the quality',/,
     &  10F8.1,/, 10I8,/, ' *** Number of bad tetrahedrons =', I8,/)

c ... parallel output
 6005 Format(' *** Distribution of tetrahedrons by the quality')
 6006 Format(10F6.1)
 6008 Format(10I6)
 6007 Format(' *** Number of bad tetrahedrons =', I8)

      Return
      End Subroutine countBadElements


C ================================================================
      Subroutine addOut(nLines, message, output)
C ================================================================
      include 'output.fd'
C ================================================================
      Character*(*) message, output(*)

      nLines = min(nLines + 1, MaxLines)
      output(nLines) = message

      Return
      End Subroutine addOut
C
      End Module mba3d_output
