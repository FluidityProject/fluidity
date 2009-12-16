      Module mba3d_cmpE
C
      contains
C
C ================================================================
      Logical Function cmpE(i1, i2, i3, IEP, nEP, iE1, iE2)
C ================================================================
C cmpE = TRUE if iE2 != iE1 and iE2 = {i1, i2, i3, *}  
C ================================================================
      Integer IEP(*), nEP(*)

C group (Local variables)
      Integer ib(3), ie(3), ip(3)

      ip(1) = i1
      ip(2) = i2
      ip(3) = i3
      Do i = 1, 3
         If(ip(i).EQ.1) Then
            ib(i) = 1
         Else
            ib(i) = nEP(ip(i) - 1) + 1
         End if
         ie(i) = nEP(ip(i))
      End do

      Do 10 i = ib(1), ie(1)
         iE2 = IEP(i)
         If(iE2.EQ.iE1) Goto 10
         Do j = ib(2), ie(2)
            If(iE2.EQ.IEP(j)) Then
               Do k = ib(3), ie(3)
                  If(iE2.EQ.IEP(k)) Then
                     cmpE = .TRUE.
                     Goto 1000
                  End if
               End do
               Goto 10
            End if
         End do
 10   Continue

      cmpE = .FALSE.
 1000 Return
      End Function cmpE
C
      End Module mba3d_cmpE
