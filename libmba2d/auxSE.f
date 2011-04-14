C ================================================================
      Logical function tangled(iE1, iE2, XYP, IPE)
C ================================================================
C Routines returns .TRUE. is triangles iE1 and iE2 are inverted.
C ================================================================
      real  XYP(2, *)
      Integer IPE(3, *)
      
      Integer iref(4)
      real  calVol, v1, v2
      Logical check22

      DATA    iref /1,2,3,1/

C ================================================================
      tangled = .FALSE.

      Do 20 i1 = 1, 3
         i2 = iref(i1 + 1)

         iP1 = IPE(i1, iE1)
         iP2 = IPE(i2, iE1)

         Do j1 = 1, 3
            j2 = iref(j1 + 1)

            jP1 = IPE(j1, iE2)
            jP2 = IPE(j2, iE2)

            If(check22(iP1, iP2, jP1, jP2)) Then
               i3  = iref(i2 + 1)
               iP3 = IPE(i3, iE1)

               j3  = iref(j2 + 1)
               jP3 = IPE(j3, iE2)

               v1 = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))
               v2 = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, jP3))

               If(v1 * v2.GE.0D0) Then
                  tangled = .TRUE.
                  goto 9000 
               End if
            End if
         End do
 20   Continue

 9000 Return
      End



C ================================================================
      Subroutine orderijk(i, j, k)
C ================================================================
C This routines sets i, j, k in non-decreasing order
C ================================================================
      If(i.LT.j) Then
         n = j
         j = i
         i = n
      End if

      If(j.LT.k) Then
         n = k
         k = j
         j = n
      End if

      If(i.LT.j) Then
         n = j
         j = i
         i = n
      End if
      Return
      End

