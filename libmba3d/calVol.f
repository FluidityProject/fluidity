      Module mba3d_calVol
C
      contains
C
C ================================================================
      Real*8 Function calVol(xy1, xy2, xy3, xy4)
C ================================================================
C The oriented volume of the tetrahedron given by 4 vertices
C ================================================================
      Real*8 xy1(3), xy2(3), xy3(3), xy4(3)

C group (Local variables)
      Real*8 v(3, 3)

      Do i = 1, 3
         v(i, 1) = xy1(i) - xy4(i)
      End do

      Do i = 1, 3
         v(i, 2) = xy2(i) - xy4(i)
      End do

      Do i = 1, 3
         v(i, 3) = xy3(i) - xy4(i)
      End do

      calVol = (v(1, 1) * (v(2, 2) * v(3, 3) - v(3, 2) * v(2, 3)) +
     &          v(2, 1) * (v(3, 2) * v(1, 3) - v(1, 2) * v(3, 3)) +
     &          v(3, 1) * (v(1, 2) * v(2, 3) - v(2, 2) * v(1, 3))) /
     &          6D0
      Return
      End Function calVol
C
      End Module mba3d_calVol
