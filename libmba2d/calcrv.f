C ================================================================
      Subroutine aniCrv(tc, XYP, iFnc, calCrv)
C ================================================================
C The function is used as a buffer between the real function and
C the code. It rescales the Cartesian coordinates.
C ================================================================
      real  tc, XYP(2)
      Integer iFnc
      EXTERNAL calCrv

      real  refXYP(2), scaXYP(2)
      Common /rescale/refXYP, scaXYP

C ================================================================
      Call calCrv(tc, XYP, iFnc)

      Do i = 1, 2
         XYP(i) = (XYP(i) - refXYP(i)) * scaXYP(i)
      End do

      Return
      End



C ================================================================
      Subroutine scale2Square(nP, XYP, flag)
C ================================================================
C Routine scales the model to the square [0.1, 0.9]^2. We allow
C 10% freedom for curved edges.
C ================================================================
      real  XYP(2, *)
      Logical flag

C ================================================================
      real  minXYP(2), maxXYP(2), scale, size

      real  refXYP(2), scaXYP(2)
      Common /rescale/refXYP, scaXYP

C ================================================================
      If(flag) Then
         Do i = 1, 2
            minXYP(i) = XYP(i, 1)
            maxXYP(i) = XYP(i, 1)
         End do

         Do n = 2, nP
            Do i = 1, 2
               minXYP(i) = min(minXYP(i), XYP(i, n))
               maxXYP(i) = max(maxXYP(i), XYP(i, n))
            End do
         End do

c  ...  add 10% for the bounding box
c         Do i = 1, 2
c            size = (maxXYP(i) - minXYP(i)) / 10
c            minXYP(i) = minXYP(i) - size
c            maxXYP(i) = maxXYP(i) + size
c         End do

         Do i = 1, 2
            refXYP(i) = minXYP(i)
            scaXYP(i) = 1D0 / (maxXYP(i) - minXYP(i))
         End do

         scale = min(scaXYP(1), scaXYP(2))

         Do i = 1, 2
            scaXYP(i) = scale
         End do

         Do n = 1, nP
            Do i = 1, 2
               XYP(i, n) = (XYP(i, n) - refXYP(i)) * scaXYP(i)
            End do
         End do
      Else
         Do i = 1, 2
            scaXYP(i) = 1D0 / scaXYP(i)
         End do

         Do n = 1, nP
            Do i = 1, 2
               XYP(i, n) = refXYP(i) + XYP(i, n) * scaXYP(i)
            End do
         End do
      End if

      Return
      End


C ================================================================
      Subroutine scaleBack(XYPi, XYPo)
C ================================================================
C  Routine computes physical coordinates of point XYPi 
C ================================================================
      real  XYPi(2), XYPo(2)

      real  refXYP(2), scaXYP(2)
      Common /rescale/refXYP, scaXYP

C ================================================================
      Do i = 1, 2
         XYPo(i) = refXYP(i) + XYPi(i) / scaXYP(i)
      End do

      Return
      End






