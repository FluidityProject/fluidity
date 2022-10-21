C*************************************************************************
      subroutine LINTRP( NT, tri, NV, vrt, LDF, F, NXY, xy, G,
     .                     imem, dmem, flagIni )
C*************************************************************************
      implicit none
      include 'lintrp.fd'
C*************************************************************************
C     Subrouine  LINTRP  restores the values of peicewise linear  func-
C     tion determined in the  nodes  of regular triangulation  anywhere
C     inside the one above.
C
C     Parameters of subroutine in order to locate (without typing):
C       NT         - number of triangles
C       tri(3, NT) - list of triangles ( 3 vertices + label )
C       NV         - number of nodes (vetices) of regular triangulation
C       vrt(2, NV) - coords of vertices
C       LDF        - leading dimension of vector function F and G
C       F(LDF, NV) - the values of vector-function in the nodes of
C                    triangulation
C       NXY        - number of points where the value of function to be
C                    determined
C       xy(2, NXY) - coords of points
C       G(LDF, NXY)- the returned values of vector-function in given
C                    points
C       imem(4 * NQ + 3 * NT + NV + 1)
C                  - auxiliary array of integers
C       dmem(2 * NQ)
C                  - auxiliary array of real numbers
C
C       flagIni    - flag 'Initialisation is needed'
C     NOTE:
C       The constant NQ lies in range   1/3 * NV <  NQ  < 2/3 * NV  and
C       dependes on spacing of triangulaion. To guaranty the robustness
C       of subroutine set NQ to  2/3 * NV. The exact sizes of auxiliary
C       storage been required  by  subroutine are returned in the first
C       elements of correspondent area.
C
C*************************************************************************
        integer NT, NV, LDF, NXY
        integer tri(3, NT)
        real vrt(2, NV), F(LDF, NV)
        real xy(2, NXY), G(LDF, NXY)
        integer imem(*)
        real dmem(*)
        logical flagIni

        real h(MaxH)
        integer i, NQT, QT, REF, ITRI
        logical flagCvxCrvBnd
C*************************************************************************
        if(flagIni) then
          QT = 5
          call INITQT( NQT, imem(QT), h, dmem(MaxH + 1) )
          do i = 1, NV
            call DROPS( NQT, imem(QT), h, dmem(MaxH + 1), vrt, i )
          enddo
C         if (NQT.gt.NV) then
C            write(*,*)'NQ=',NQT, 'NV=',NV
C            stop
C         end if
          REF  = QT + 4 * NQT
          ITRI = REF + NV + 1
        else
          QT  = imem(1)
          NQT = imem(2)
          REF = imem(3)
          ITRI =imem(4)

          do i = 1, MaxH
             h(i) = dmem(i)
          end do
        end if

        flagCvxCrvBnd = .true.
C       flagCvxCrvBnd is flag 
C                 'Nearby a true curved boundary patch the domain is convex'
C      If the domain is polygonal or nearby a true curved boundary
C    patch the domain is concave, and for a given point there is no an
C    element containing it within tolerance PREC=10^{-6} , then the ERROR
C    is assumed. The code could not find an element within PREC tolerance
C    for MaxTrials.
C      If the domain has a curved boundary and nearby a true curved boundary 
C    patch the domain is convex, it is possible that a given
C    point is in the domain but out of the mesh. In this case the tolerance
C    PREC will be relaxed until an element will be found. 

        call RESTORE( NT, tri, NV, vrt, LDF, F, NXY, xy, G,
     .                imem(QT), dmem(MaxH + 1),
     .                imem(REF), imem(ITRI),
     .                h, flagIni, flagCvxCrvBnd )

c        imem(1) = NQT * 4 + 3 * NT + NV + 1
c        dmem(1) = NQT * 2

        if(flagIni) then
          imem(1) = QT
          imem(2) = NQT
          imem(3) = REF
          imem(4) = ITRI

          do i = 1, MaxH
             dmem(i) = h(i)
          end do
        end if
      return
      end



C*************************************************************************
      subroutine INITQT( NQT, QT, h, XYc )
C*************************************************************************
      implicit none
      include 'lintrp.fd'
C*************************************************************************
C Routine initializes a quadtree sructure
C*************************************************************************
        integer NQT, QT(4, *)
        real XYc(2, *)
        real h(*)

        integer k, NEWQT
C*************************************************************************
        h(1) = 0.5
        do k = 2, MaxH
          h(k) = h(k-1) / 2
        enddo
        NQT = 0
        NQT = NEWQT( NQT, QT )
        XYc(1, NQT) = 0.5
        XYc(2, NQT) = 0.5
      return
      end



C*************************************************************************
      integer function NEWQT ( NQT, QT )
C*************************************************************************
C Allocation of new quadtree
C*************************************************************************
      implicit none
        integer NQT
        integer QT(4, *)

        integer i

C*************************************************************************
        NQT = NQT + 1
        do i = 1, 4
          QT(i, NQT) = 0
        enddo
        NEWQT = NQT
      return
      end



C*************************************************************************
      subroutine SETIJ ( XYc, XY, I, J )
C*************************************************************************
C  Definition of quadrant of point XY with respect to the centre of
C  quadtree XYc
C*************************************************************************
      implicit none
        real XYc(2), XY(2)
        integer I, J

C*************************************************************************
        if ( XY(1).GT.XYc(1) ) then
          I = 2
        else
          I = 1
        endif
        if ( XY(2).GT.XYc(2) ) then
          J = 2
        else
          J = 1
        endif
      return
      end



C*************************************************************************
      subroutine DROPS ( NQT, QT, h, XYc, XY, idx )
C*************************************************************************
      implicit none
      include 'lintrp.fd'
C*************************************************************************
C  The arrangement of a new point in quadtree
C*************************************************************************
        integer NQT, idx
        integer QT(2, 2, *)
        real XYc(2, *), XY(2, *)
        real h(*)

        integer L, ip, ptr, i, j, new, NEWQT, dir(2)

C*************************************************************************
        dir(1)=-1
        dir(2)= 1

        ip = 1
        L = 1
        do while ( .TRUE. )
          L = L + 1
          if ( L.GT.MaxH ) 
     &      Call errMes(6201, 'DROPS', 'No memory for quadtree')

          call SETIJ( XYc(1, ip), XY(1, idx), i, j )
          ptr = QT(i, j, ip)
          if ( ptr.EQ.0 ) then
            QT(i, j, ip) = idx
            return
          else if ( ptr.LT.0 ) then
            ip = -ptr
          else if ( ptr.GT.0 ) then
            new = NEWQT( NQT, QT )
            QT(i, j, ip) = -new
            XYc(1, new) = XYc(1, ip) + dir(i)*h(L)
            XYc(2, new) = XYc(2, ip) + dir(j)*h(L)
            call SETIJ( XYc(1, new), XY(1, ptr), i, j )
            QT(i, j, new) = ptr
            ip = new
          endif
        enddo
      return
      end



C*************************************************************************
      real function SQRDST( a, b )
C*************************************************************************
C  Function returns the square of distance between points of plane
C*************************************************************************
      real a(2), b(2)
  
        SQRDST = (a(1) - b(1))**2 + (a(2) - b(2))**2
      return
      end



C*************************************************************************
      subroutine ORDER2( XYc, h, XY, ord, sqrd )
C*************************************************************************
C  Ordering of subquads of quadtree QT by their distantness from
C  the given point XY
C*************************************************************************
      implicit none
        integer ord(4)
        real  XYc(2), XY(2), h, sqrd(4)

        integer i, j
        real dist(2), sqr(2, 2)
        integer ofs

C*************************************************************************
        ofs(i, j) = i + 2 * (j - 1)

        do i = 1, 2
          dist(i)   = abs( XY(i) - XYc(i) )
          sqr(i, 1) = dist(i)**2
          sqr(i, 2) = dim( dist(i), h )**2
        enddo
        call SETIJ( XYc, XY, i, j )

        ord(1) = ofs(i, j)
        sqrd(1) = sqr(1, 1) + sqr(2, 1)
        if ( dist(1).LT.dist(2) ) then
          ord(2) = ofs(3-i, j)
          ord(3) = ofs(i, 3-j)
          sqrd(2) = sqr(1, 1) + sqr(2, 2)
          sqrd(3) = sqr(1, 2) + sqr(2, 1)
        else
          ord(2) = ofs(i, 3-j)
          ord(3) = ofs(3-i, j)
          sqrd(2) = sqr(1, 2) + sqr(2, 1)
          sqrd(3) = sqr(1, 1) + sqr(2, 2)
        endif
        ord(4) = ofs(3-i, 3-j)
        sqrd(4) = sqr(1, 2) + sqr(2, 2)
      return
      end



C*************************************************************************
      integer function NEARST ( QT, XYc, xy, point, h )
C*************************************************************************
      implicit none
      include 'lintrp.fd'
C*************************************************************************
C  This function returns the index of the quadtree point which is
C  the nearest to given one
C*************************************************************************
        integer QT(4, *)
        real XYc(2, *), xy(2, *), point(2)
        real h(*)

        real sqrd(4, MaxH), SQRDST, sqdist, min
        integer L, ip(MaxH), k(MaxH), ord(4, MaxH), ptr

C*************************************************************************
        L   = 1
        ptr = -1
        min = 2.0
        NEARST = 0
        do while ( .TRUE. )
          if ( ptr.LT.0 ) then
            L = L + 1
            if ( L.GT.MaxH ) 
     &        Call errMes(6202, 'NEAREST', 'No memory for quadtree')

            ip(L) = -ptr
            call ORDER2( XYc(1, -ptr), h(L), point, ord(1, L),
     .                   sqrd(1, L) )
            k(L) = 1
          else
            if ( ptr.GT.0 ) then
              sqdist = SQRDST( point, xy(1, ptr) )
              if ( min.GT.sqdist ) then
                min = sqdist
                NEARST = ord(k(L),L) + 4*(ip(L) - 1)
              endif
            endif
            do while ( .TRUE. )
              k(L) = k(L) + 1
              if ( k(L).LE.4 ) then
                 if( sqrd(k(L), L).LT.min ) goto 101
              endif
              L = L - 1
              if ( L.EQ.1 ) goto 102
            enddo
          endif
 101      continue
          ptr = QT(ord(k(L), L), ip(L))
        enddo
 102    continue
      return
      end



C*************************************************************************
      integer function BASETRI(QT,XYc, vrt,NTR,tri, ref,itri, XY,h,PREC)
C*************************************************************************
      implicit none
      include 'lintrp.fd'
C*************************************************************************
C  Following function is used to determine the underlying triangle
C  with relative tolerance reltol
C  for point XY being inside of triangulation
C*************************************************************************
        integer QT(*), tri(3, *), ref(*), itri(*), NTR
        real XYc(2, *), vrt(2, *), XY(2)
        real h(*), PREC  

        integer i, k, n, idx, ip, buf(2, MaxBuf), NEARST
        integer j, m, jj, idx2, iT, tbuf(0:MaxTrSSEl),MaxTrials
        logical ENCLOSE

C*************************************************************************
        BASETRI = 0
        k = 0
        MaxTrials = max(int(NTR * RelativeTrials), 20)

        do while ( k.lt.MaxTrials )
          ip = NEARST( QT, XYc, vrt, XY, h )
          if ( ip.EQ.0 ) goto 101
          tbuf(0) = 0
          idx = QT(ip)
          do i = ref(idx), ref(idx+1)-1
            do jj = 1, tbuf(0)
               if (tbuf(jj).eq.itri(i)) goto 1
            end do
            if ( ENCLOSE( XY, vrt, tri(1, itri(i)), PREC ) ) then
               BASETRI = itri(i)
               goto 101
            endif
            tbuf(0) = tbuf(0) + 1
            if ( tbuf(0).GT.MaxTrSSEl ) 
     &         Call errMes(6203, 'BASETRI', 'No memory for tbuffer')
            tbuf( tbuf(0) ) = itri(i)
1           continue
          enddo

C ...  checking for neighbooring triangles   
          do i = ref(idx), ref(idx+1)-1
            iT = itri(i)
            Do j = 1, 3
               idx2 = tri(j, iT)
               Do m = ref(idx2), ref(idx2+1)-1
                  do jj = 1, tbuf(0)
                    if (tbuf(jj).eq.itri(m)) goto 2
                  end do
                  if ( ENCLOSE( XY, vrt, tri(1, itri(m)), PREC ) )then
                     BASETRI = itri(m)
                     goto 101
                  endif
                  tbuf(0) = tbuf(0) + 1
                  if ( tbuf(0).GT.MaxTrSSEl ) 
     &               Call errMes(6203,'BASETRI','No memory for tbuffer')
                  tbuf( tbuf(0) ) = itri(m)
2                 continue
               End do   
            End do
          enddo
C ...  end of


          k = k + 1
          if ( k.GT.MaxBuf ) 
     &      Call errMes(6203, 'BASETRI', 'No memory for buffer')
          buf(1, k) = ip
          buf(2, k) = idx
          QT(ip) = 0
        enddo
 101    do n = k, 1, -1
          QT(buf(1, n)) = buf(2, n)
        enddo
      return
      end



C*************************************************************************
      subroutine RESTORE(  NTR, tri, NVR, vrt, LDF, f,
     .                     NXY,  xy, g, QT, XYc, ref, itri,
     .                     h, flagIni, flagCvxCrvBnd )
C*************************************************************************
C  Code builds the transposed list of triangle grid and restores
C  the values of function F in the points of interest with aid
C  quadtree structure
C*************************************************************************
      implicit none
        integer NVR, NTR, LDF, NXY
        integer QT(2, 2, *), tri(3, NTR), ref(*), itri(*)
        real XYc(2, *)
        real vrt(2, NVR), xy(2, NXY)
        real f(LDF, NVR), g(LDF, NXY)
        real h(*)
        logical flagIni, flagCvxCrvBnd
c Predefined tolerance for checking whether the point belongs to an element
      real PREC
      parameter( PREC = 1D-6 )
c Predefined number of possible relaxations of the above tolerance (6 orders of 10, here)
      integer    MaxRelax
      parameter( MaxRelax = 50 )


        integer i, j, idx,  ip(3), BASETRI
        real a, b, c, d, PREC1

        real  detf

        integer irel
C*************************************************************************
        if(flagIni) then
          do i = 1, NVR+1
            ref(i) = 0
          enddo
          do i = 1, NTR
            do j = 1, 3
              idx = tri(j, i)
              ref(idx) = ref(idx) + 1
            enddo
          enddo
          ref(1) = ref(1) + 1
          do i = 2, NVR+1
            ref(i) = ref(i-1) + ref(i)
          enddo
          do i = 1, NTR
            do j = 1, 3
              idx = tri(j, i)
              ref(idx) = ref(idx) - 1
              itri(ref(idx)) = i
            enddo
          enddo
        end if

        do i = 1, NXY
          PREC1 = PREC
          irel = 0
11        idx = BASETRI(QT,XYc, vrt,NTR,tri, ref, itri,xy(1,i), h,PREC1)
          if(idx.le.0) then
c            write(*,'(A,2F11.7,A)') 'Point: ', xy(1, i), xy(2, i),
c    &             ' is out of the mesh even approximately. '
C if Curvelinear boundaries exist and nearby the domain is convex, 
C    then probably the node is out of the mesh
C Relaxation of PREC is applied
            if (flagCvxCrvBnd) then
               if (PREC1.lt.5D-2) then
                  PREC1 = PREC1*10
               else
                  PREC1 = PREC1*3.16
               end if
               if (irel.lt.MaxRelax) then
                  irel = irel + 1
                  goto 11
               else
                  Call errMes(6204, 'lintrp2D',
     .                 'Max number of relaxations for PREC is reached')
               end if
C  if no Curvelinear boundaries, then error is assumed or inconsistency
C of parameters in lintrp.fd. REINSTALL constants in lintrp.fd!
            else
               Call errMes(6204, 'lintrp2D',
     .             'Failed to find element within PREC tolerance')
            end if
          end if


          do j = 1, 3
            ip(j) = tri(j, idx)
          enddo
          a = detf(2, 3, i,ip,xy,vrt)
          b = detf(3, 1, i,ip,xy,vrt)
          c = detf(1, 2, i,ip,xy,vrt)
          d = a + b + c

          if(d.eq.0D0) then
             a = 1D0
             d = 1D0
          end if

          do j = 1, LDF
             G(j, i) = (a * F(j, ip(1)) +
     .                  b * F(j, ip(2)) + c * F(j, ip(3))) / d
          end do
        enddo
      return
      end



C*************************************************************************
      real function detf(i,j,ii,ip,xy,vrt)
C*************************************************************************
      implicit none
      integer i,j,ii,ip(*)
      real xy(2,*),vrt(2,*)

        detf = (vrt(1, ip(i)) - xy(1,ii))*(vrt(2, ip(j)) - xy(2,ii))
     &       - (vrt(2, ip(i)) - xy(2,ii))*(vrt(1, ip(j)) - xy(1,ii))
      return
      end



C*************************************************************************
      logical function ENCLOSE ( XY, vrt, TRI, tol )
C*************************************************************************
C  Determines does the point XY belongs to triangle TRI
C*************************************************************************
      implicit none
        real XY(2), vrt(2, *), tol
        integer TRI(3)

        integer i, j, k
        real a, b, c, alph, beta
        real x, y, F, frac

C*************************************************************************
        F( alph, beta ) = a * alph + b * beta + c
        x(i) = vrt(1, TRI(i))
        y(i) = vrt(2, TRI(i))

        k = 3
        j = 2
        do i = 1, 3
          a = y(j) - y(i)
          b = x(i) - x(j)
          c = x(j) * y(i) - x(i) * y(j)
          frac = F( XY(1), XY(2) ) / F( x(k), y(k) )
          if  ( frac.LE.-tol ) then
            ENCLOSE = .FALSE.
            return
          endif
          j = k
          k = i
        enddo
        ENCLOSE = .TRUE.

      return
      end


C ================================================================
      real Function sizeQT(point, imem, dmem)
C ================================================================
      include 'lintrp.fd'
C ================================================================
C  The function returns the size of the quadtree cell which contains
C  the given point with coords point(2).
C
C  PARAMETERS:
C     point(2) - user given point inside the unit square (0,1)^2
C
C     dmem(*)  - real  working memory of LINTRP2D
C     imem(*)  - integer working memory of LINTRP2D
C
C  REMARK 1. This function is the wrapping for function SizeHost.
C
C ================================================================
      real point(2), dmem(*)
      Integer          imem(*)

      real SizeHost

      iQT = imem(1)
      sizeQT = SizeHost(imem(iQT), dmem(MaxH + 1), point, dmem)

      Return
      End



C ================================================================
      real Function SizeHost ( QT, XYc, point, h )
C ================================================================
      implicit none
      include 'lintrp.fd'
C ================================================================
      integer QT(4, *)
      real XYc(2, *), point(2)
      real h(*)

      real sqrd(4, MaxH)
      integer L, ip(MaxH), ord(4, MaxH), ptr

C ================================================================
      L   = 1
      ptr = -1
      SizeHost = 1D0

      do while ( .TRUE. )
         if ( ptr.LT.0 ) then
            L = L + 1
            if ( L.GT.MaxH ) Call errMes(1009, 'SizeHost', 
     &                            'local parameter MaxH is small')

            ip(L) = -ptr
            call ORDER2( XYc(1, -ptr), h(L), point, ord(1, L),
     .                   sqrd(1, L) )
         else
            SizeHost = h(L)
            return
         endif
         ptr = QT(ord(1, L), ip(L))
      enddo

      Call errMes(6204, 'SizeHost', 'Host cell was not found')

      return
      end


