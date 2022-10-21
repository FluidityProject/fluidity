C ================================================================
      Subroutine updQa(n, XYP, IPE, IEE, qE)
C ================================================================
C Initial quality modification for tangled elements and their 
C edge-neighboors.
C ================================================================
      real  XYP(2, *), qE(*)
      Integer IPE(3, *), IEE(3, *)

C (Local variables)
      Integer ip(4)
      real  calVol, v1, v2
      Logical check22

C ================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      Do 100 i1 = 1, 3
         iE = IEE(i1, n)
         If(iE.LE.0) goto 100

         i2 = ip(i1 + 1)
         i3 = ip(i2 + 1)

         iP1 = IPE(i1, n)
         iP2 = IPE(i2, n)

         Do j1 = 1, 3
            j2 = ip(j1 + 1)

            jP1 = IPE(j1, iE)
            jP2 = IPE(j2, iE)

            If(check22(iP1, iP2, jP1, jP2)) Then
               i3  = ip(i2 + 1)
               iP3 = IPE(i3, n)

               j3  = ip(j2 + 1)
               jP3 = IPE(j3, iE)

               v1 = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))
               v2 = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, jP3))

               If(v1 * v2.GE.0D0) Then
                  qE(n)  = -abs(qE(n))
                  qE(iE) = -abs(qE(iE))
               End if
               goto 100
            End if
         End do
 100  Continue

      Return
      End


C ================================================================
      Subroutine updQb(nEs, lE, iEs, XYP, IPEs, qEs)
C ================================================================
C Dynamic quality modification for tangled elements inside
C a super-element.
C
C Remark: non-efficient, time-consuming, but robust algorithm.
C ================================================================
      real  XYP(2, *), qEs(*)
      Integer iEs(*), IPEs(3, *)

C group (Local variables)
      Integer ip(4)
      real  calVol, v1, v2
      Logical check22

C ================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      Do 100 i1 = 1, 3
         i2 = ip(i1 + 1)

         iP1 = IPEs(i1, nEs)
         iP2 = IPEs(i2, nEs)

         Do 20 k = 1, lE
            If(iEs(k).LT.0)  goto 20
            If(k.EQ.nEs)     goto 20

            Do j1 = 1, 3
               j2 = ip(j1 + 1)

               jP1 = IPEs(j1, k)
               jP2 = IPEs(j2, k)

               If(check22(iP1, iP2, jP1, jP2)) Then
                  i3  = ip(i2 + 1)
                  iP3 = IPEs(i3, nEs)

                  j3  = ip(j2 + 1)
                  jP3 = IPEs(j3, k)

                  v1 = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))
                  v2 = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, jP3))

                  If(v1 * v2.GE.0D0) Then
                     qEs(nEs) = -abs(qEs(nEs))
                     qEs(k)   = -abs(qEs(k))
                  End if

                  goto 100
               End if
            End do
 20      Continue
 100  Continue

 1000 Return
      End



C ================================================================
      Logical Function  chkTangled(lE, iEs, IPEs)
C ================================================================
C Local mesh modifications for tangled mesh may result is a 
C topologically wrong mesh (chkTangled = TRUE).
C ================================================================
      Integer iEs(*), IPEs(3, *)

C group (Local variables)
      Integer ip(4), tEdge
      Logical check22

C ================================================================
      chkTangled = .TRUE.

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

c ...  more than two triangles with a common edge
       Do 20 n = 1, lE
          iEt = iEs(n)
          If(iEt.LE.0) goto 20
          
          tEdge = 0
          Do i1 = 1, 3
             i2 = ip(i1 + 1)

             iP1 = IPEs(i1, n)
             iP2 = IPEs(i2, n)

             nEdge = 0
             Do 10 m = 1, lE
                jEt = iEs(m)
                If(m.EQ.n .OR. jEt.LE.0) goto 10

                Do j1 = 1, 3
                   j2 = ip(j1 + 1)

                   jP1 = IPEs(j1, m)
                   jP2 = IPEs(j2, m)
                   If(check22(iP1, iP2, jP1, jP2)) nEdge = nEdge + 1
                End do
 10          Continue            
             If(nEdge.GT.2) goto 1000
             tEdge = tEdge + nEdge
          End do 

c  ...  no neighboors
          If(tEdge.EQ.0) goto 1000
 20    Continue
       

c ...  there is no topological defects
       chkTangled = .FALSE.

 1000  Return
       End

