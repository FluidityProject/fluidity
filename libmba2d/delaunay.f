C ================================================================
      Subroutine Delaunay(
C ================================================================
     &           nP, nE, XYP, IPE,
     &           MaxWr, MaxWi, rW, iW)
C ================================================================
      implicit none
      include 'makS.fd'
C ================================================================
C  The routine builds the Delaunay triangulation from the existing
C  triangution by swapping edges in pairs of triangles.
C
C  Data flow:
C      1. calculate map E->E
C      2. loop over triangles
C         2.1 check criterium for sum of opposite angles 
C         2.2 swap triangles
C         2.3 update the map
C
C  *** Remarks 
C ================================================================
      Integer  nP, nE

      real   XYP(2, *)
      Integer  IPE(3, *)

      Integer  MaxWr, MaxWi, iW(*)
      real   rW(*)

      Integer  iref(4)
      Logical  check22, DelonePair, flagREPEAT, flagDELONE

      Integer  i1,i2,i3, j1,j2,j3, k, n
      Integer  iP1,iP2,iP3, jP1,jP2,jP3, iE1,iE2, jE2,jE3
      Integer  inEp, iIEP, iIEE, iEnd, iEmem, jEmem, kEmem
      Integer  nswap, nswapold, nloop
      real   sa1, sa2

      DATA     iref /1,2,3,1/

C ================================================================
      iIEE = 1
      inEP = iIEE + 3 * nE
      iIEP = inEP + nP
      iEnd = iIEP + 3 * nE

      If(iEnd.GT.MaxWi) Call errMes(1001,
     &      'delaunay', 'MaxWi is too small')

      Call listE2E(nP, nE, IPE, iW(iIEE), iW(inEP), iw(iIEP))

c ... initialize loops
      nswapold = nE
      nloop = 0

c ... new loop
 100  flagREPEAT = .FALSE.
      nswap = 0

      Do n = 1, nE
         Do 30 i1 = 1, 3
            iEmem = iIEE + 3 * (n - 1) - 1
            iE1 = iW(iEmem + i1)
            If(iE1.LE.n) goto 30

            i2 = iref(i1 + 1)
            i3 = iref(i2 + 1)

            iP1 = IPE(i1, n) 
            iP2 = IPE(i2, n) 
            iP3 = IPE(i3, n) 

            Do j1 = 1, 3
               j2 = iref(j1 + 1) 

               jP1 = IPE(j1, iE1)
               jP2 = IPE(j2, iE1)
               If(check22(iP1, iP2, jP1, jP2)) Then
                  j3 = iref(j2 + 1) 
                  jP3 = IPE(j3, iE1)
                  goto 10 
               End if
            End do

 10         Continue
            flagDELONE = DelonePair(XYP(1, iP1), XYP(1, iP2), 
     &                              XYP(1, iP3), XYP(1, jP3), sa1)

            If(.NOT.flagDELONE) Then
               flagREPEAT = .TRUE.

               flagDELONE = DelonePair(XYP(1, iP3), XYP(1, jP3), 
     &                                 XYP(1, iP1), XYP(1, jP2), sa2)
               If(flagDELONE .OR. sa2.GT.sa1) Then
                  nswap = nswap + 1
                  iE2 = iW(iEmem + i2) 

                  jEmem = iIEE + 3 * (iE1 - 1) - 1
                  jE2 = iW(jEmem + j2)
                  jE3 = iW(jEmem + j3)
                  If(iP1.NE.jP1) Call swapii(jE3, jE2)

                  IPE(i1, n) = iP3
                  IPE(i2, n) = jP3
                  IPE(i3, n) = iP1

                  IPE(j1, iE1) = iP3
                  IPE(j2, iE1) = jP3
                  IPE(j3, iE1) = iP2

                  iW(iEmem + i2) = jE3
                  If(jE3.GT.0) Then
                     kEmem = iIEE + 3 * (jE3 - 1) - 1 
                     Do k = 1, 3
                        If(iW(kEmem + k).EQ.iE1) Then
                           iW(kEmem + k) = n
                           goto 20
                        End if
                     End do

c                    Call draw_T(nP, 0, nE, XYP, iW, iW, IPE, 'fin.ps')
                     Call errMes(6006, 
     &                   'delaunay.f', 'Cannot build Delaunay mesh')
                  End if

 20               iW(jEmem + j2) = jE2
                  iW(jEmem + j3) = iE2
                  If(iE2.GT.0) Then
                     kEmem = iIEE + 3 * (iE2 - 1) - 1
                     Do k = 1, 3
                        If(iW(kEmem + k).EQ.n) Then
                           iW(kEmem + k) = iE1 
                           goto 30
                        End if
                     End do

c                    Call draw_T(nP, 0, nE, XYP, iW, iW, IPE, 'fin.ps')
                     Call errMes(6006, 
     &                   'delaunay.f', 'Cannot build Delaunay mesh')
                  End if
               End if
            End if
 30      Continue
      End do

      If(flagREPEAT) Then
         nloop = nloop + 1
 
         If(nloop.GT.nE) Then
c          Call draw_T(nP, 0, nE, XYP, iW, iW, IPE, 'fin.ps')
           Call wrnMes(6006, 
     &                'delaunay.f', 'Cannot build Delaunay mesh')
           Return
         End if

         nswapold = nswap
         goto 100
      End if
 
      Return
      End



C ================================================================
      Logical Function DelonePair(xy1, xy2, xy3, xy4, sa)
C ================================================================
C  This roune checks that the pair {1,2,3} and {1,2,4} is the
C  Delaunay pair. It returns .FALSE.otherwise.
C ================================================================
      real  xy1(2), xy2(2), xy3(2), xy4(2)

      real  sa, sb, ca, cb
C ================================================================
      DelonePair = .FALSE.

      ca = 0
      cb = 0

      Do i = 1, 2 
         ca = ca + (xy4(i) - xy1(i)) * (xy4(i) - xy2(i)) 
         cb = cb + (xy3(i) - xy1(i)) * (xy3(i) - xy2(i)) 
      End do
            
c ... Delaunay condition fails
      If(ca.LT.0 .AND. cb.LT.0) goto 9000

c ... Delaunay condition true
      If(ca.GE.0 .AND. cb.GE.0) Then
         DelonePair = .TRUE.
         goto 9000
      End if

c ... the rest of the formula
      sa = (xy4(1) - xy1(1)) * (xy4(2) - xy2(2)) 
     &   - (xy4(1) - xy2(1)) * (xy4(2) - xy1(2))

      sb = (xy3(1) - xy1(1)) * (xy3(2) - xy2(2)) 
     &   - (xy3(1) - xy2(1)) * (xy3(2) - xy1(2)) 

c ... checking when swapping is not possible
c     If(sa * sb.GE.0) goto 9000

      sa = abs(sa) * cb + ca * abs(sb)
      If(sa.GE.0) DelonePair = .TRUE.
        
 9000 Return
      End



C ==============================================================
      Subroutine RandXY(xy1, xy2, xy3, xyc, r)
C ==============================================================
C Routine computes the center and radius of the circle 
c curcumscribed around triangle defined by three vertices.
C ==============================================================
      real  xy1(2), xy2(2), xy3(2)
      real  xyc(2), r

      real  a, b, c, x1, y1, x2, y2, r1, r2
C ==============================================================
      x1 = xy1(1) - xy3(1)
      y1 = xy1(2) - xy3(2)

      x2 = xy2(1) - xy3(1)
      y2 = xy2(2) - xy3(2)

      r1 = x1 ** 2 + y1 ** 2
      r2 = x2 ** 2 + y2 ** 2

      a = x1 * y2 - y1 * x2
      b = r1 * y2 - r2 * y1
      c = r1 * x2 - r2 * x1  

      xyc(1) = b / (2 * a)
      xyc(2) =-c / (2 * a)

      r = sqrt(xyc(1) ** 2 + xyc(2) ** 2)

      xyc(1) = xyc(1) + xy3(1)
      xyc(2) = xyc(2) + xy3(2)

      Return
      End


C ==========================================================
      Subroutine draw_D(nP, nF, nE, XYP, ICP, IPF, IPE, fName)
C ==========================================================
      include 'colors.fd'
C ==========================================================
C Routine converts the Delaunay mesh into a postcript file fName.
C ==========================================================
C group (M)
      real  XYP(2, *)
      Integer ICP(*), IPF(4, *), IPE(3, *)

C group (File)
      Character*(*) fName

C group (Local variables)
      Real  mmTOpt, kx, ky
      Real  x1, y1, x2, y2, x3, y3
      Real  fxmax, fymax

      real  calEdge
      Logical ifXnode

      real  xyc(2), rOut, rad
      Character*30 fNameExt

C ==========================================================
      i = 1
      Do while( fName(i:i+2) .NE. '.ps')
         i = i + 1
      End do
      fNameExt = fName(1:i+2)

      mmTOpt = 72.0 / 25.4
      Pi = 4.0*atan(1.0)
      rdTOdg = 180.0 / Pi
      xBig = 150.0
      yBig = 150.0

      fxmax = XYP(1, 1)
      fxmin = fxmax
      fymax = XYP(2, 1)
      fymin = fymax

      Do n = 2, nP
         fxmin = min(fxmin, real(XYP(1, n)))
         fxmax = max(fxmax, real(XYP(1, n)))
         fymin = min(fymin, real(XYP(2, n)))
         fymax = max(fymax, real(XYP(2, n)))
      End do

      kx = xBig / (fxmax-fxmin) * mmTOpt
      ky = yBig / (fymax-fymin) * mmTOpt

      If(kx.GT.ky) kx = ky

      ibx = (fxmax-fxmin) * kx + 10
      iby = (fymax-fymin) * kx + 10

      Open(30, file = fNameExt, status='UNKNOWN')

      Call headerPS(30, ibx, iby)

      Do n = 1, nP
         x1 = (XYP(1, n) - fxmin) * kx
         y1 = (XYP(2, n) - fymin) * kx
c        Write(30,*) x1,y1, ' m', x1, y1, n, ' intTOtext ctext'
      End do

      Do 10 n = 1, nE
         If(IPE(1, n).EQ.0) goto 10
         x1 = (XYP(1, IPE(1, n)) - fxmin) * kx
         y1 = (XYP(2, IPE(1, n)) - fymin) * kx

         x2 = (XYP(1, IPE(2, n)) - fxmin) * kx
         y2 = (XYP(2, IPE(2, n)) - fymin) * kx

         x3 = (XYP(1, IPE(3, n)) - fxmin) * kx
         y3 = (XYP(2, IPE(3, n)) - fymin) * kx

         Write(30,*) 'newpath'
         Write(30,*) x1,y1, ' m', x2,y2, ' l', x3,y3, ' l'
         Write(30,'(A)') ' closepath cBlack stroke'

c  ...   mark fix points
         rad = fxmax - fxmin
         Do i = 1, 3
            iP1 = IPE(i, n)
            Do j = i + 1, 3
               iP2 = IPE(j, n)

               rad = min(rad, calEdge(XYP(1, iP1), XYP(1, iP2)))
            End do
         End do
         rad = rad * kx / 4
         rad = min(rad, 1D0)

         Do i = 1, 3
            iP1 = IPE(i, n)
            If(ifXnode(ICP(iP1), jVnode)) Then
               x1 = (XYP(1, iP1) - fxmin) * kx
               y1 = (XYP(2, iP1) - fymin) * kx
               Write(30,*) x1,y1, rad, ' cRed fillCircle'
            End if
         End do
 10   Continue


c ... draw boundaries using thick lines
      Do 20 n = 1, nF
         If(IPF(1, n).LE.0) goto 20
         iP1 = IPF(1, n)
         iP2 = IPF(2, n)

         x1 = (XYP(1, iP1) - fxmin) * kx
         y1 = (XYP(2, iP1) - fymin) * kx

         x2 = (XYP(1, iP2) - fxmin) * kx
         y2 = (XYP(2, iP2) - fymin) * kx

         Write(30,*) x1,y1, ' m', x2,y2, ' l '
c        Write(30,*) x1,y1, x2,y2, ' emptyArrow '

         ic = IPF(4, n) / 5
         ic = max(0, IPF(4, n) - 5 * ic)
         Write(30,'(A,I1,A)') 'c', ic, ' stroke'

         x1 = (x1 + x2) / 2
         y1 = (y1 + y2) / 2
c        Write(30,*) x1,y1, ' m', x1, y1, n, ' intTOtext ctext'
 20   Continue


c ... draw the circumcribed circle
      Do 30 n = 1, nE
         If(IPE(1, n).EQ.0) goto 30

         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         Call RandXY(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), xyc, rOut)

         x1  = (xyc(1) - fxmin) * kx
         y1  = (xyc(2) - fymin) * kx
         rad = rOut * kx

         Write(30,*) ' gsave', x1, y1,
     &               ' translate newpath 0 0 ', rad,
     &               ' 0 360 arc closepath cGray stroke grestore'
 30   Continue

      Write(30,'(A)') 'showpage '
      Close(30)

      Return
      End

