C ================================================================
      real Function sqrEdge(xy1, xy2)
C ================================================================
C Routine computes sqruare distance between two points
C ================================================================
      real xy1(2), xy2(2)

      sqrEdge = (xy1(1) - xy2(1)) ** 2 + (xy1(2) - xy2(2)) ** 2
      Return
      End


C ================================================================
      real Function calEdge(xy1, xy2)
C ================================================================
C Routine computes distance between two points
C ================================================================
      real xy1(2), xy2(2)

      calEdge = sqrt((xy1(1) - xy2(1)) ** 2 +
     &                (xy1(2) - xy2(2)) ** 2)
      Return
      End



C ================================================================
      real Function calNorm(xy1)
C ================================================================
C Routine calculates norm of vector xy1 
C ================================================================
      real xy1(2)

      calNorm = sqrt(xy1(1) ** 2 + xy1(2) ** 2)
      Return
      End



C ================================================================
      Subroutine extNormal(xy1, xy2, xy3, xyn)
C ================================================================
C Routines compute external normal vector to the edge {xy1, xy2}
C of triangle {xy1, xy2, xy3}
C ================================================================
      real xy1(2), xy2(2), xy3(2), xyn(2)
      real x, y, d

      x = xy2(1) - xy1(1)
      y = xy2(2) - xy1(2)

      d = sqrt(x * x + y * y)
 
      xyn(1) = -y / d
      xyn(2) =  x / d

c ... orientation
      x = xy3(1) - xy1(1)
      y = xy3(2) - xy1(2)

      If( x*xyn(1) + y*xyn(2).GT.0D0 ) Then
         xyn(1) = -xyn(1)
         xyn(2) = -xyn(2)
      End if 

      Return
      End


C ================================================================
      Subroutine calNormal(xy1, xy2, xyn)
C ================================================================
C Routines compute a normal vector to the edge {xy1, xy2}
C ================================================================
      real xy1(2), xy2(2), xyn(2)
      real x, y, d

      x = xy2(1) - xy1(1)
      y = xy2(2) - xy1(2)

      d = sqrt(x * x + y * y)
 
      xyn(1) = -y / d
      xyn(2) =  x / d

      Return
      End



C ================================================================
      real Function DotMul(a, b)
C ================================================================
C Routine computes scalar product of two vectors
C ================================================================
      real a(2), b(2)

      DotMul = a(1) * b(1) + a(2) * b(2)
      Return
      End



C ================================================================
      real Function VecMul(a, b)
C ================================================================
C Routine computes vector product a x b which is a number in 2D.
C ================================================================
      real a(2), b(2)

      VecMul = a(1) * b(2) - a(2) * b(1)
      Return
      End



C ================================================================
      Logical Function check1j(i1, j)
C ================================================================
C check1j = TRUE if i1 belongs to the set j(3).
C ================================================================
      Integer j(3)

      check1j = .TRUE.
      If(i1.EQ.j(1)) goto 1000
      If(i1.EQ.j(2)) goto 1000
      If(i1.EQ.j(3)) goto 1000
 
      check1j = .FALSE.
 1000 Return
      End



C ================================================================
      Logical Function check12(i1, j1, j2)
C ================================================================
C check12 = TRUE if i1 belongs to the set {j1, j2}.
C ================================================================
      check12 = .TRUE.
      If(i1.EQ.j1) goto 1000
      If(i1.EQ.j2) goto 1000

      check12 = .FALSE.
 1000 Return
      End


C ================================================================
      Logical Function check22(i1, i2, j1, j2)
C ================================================================
C check22 = TRUE if pair {i1, i2} coinsides with {j1, j2}.
C ================================================================
      check22 = .FALSE.
      If(i1.NE.j1 .AND. i1.NE.j2) goto 1000
      If(i2.NE.j1 .AND. i2.NE.j2) goto 1000

      check22 = .TRUE.
 1000 Return
      End


C ================================================================
      Subroutine swapii(i1, i2)
C ================================================================
C Routine swaps two integers 
C ================================================================
      i  = i1
      i1 = i2
      i2 = i

      Return
      End



C ================================================================
      Subroutine swapdd(d1, d2)
C ================================================================
C Routine swaps two real numbers
C ================================================================
      real d, d1, d2
      d  = d1
      d1 = d2
      d2 = d

      Return
      End

