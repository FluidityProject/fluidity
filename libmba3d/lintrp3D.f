      Module mba3d_lintrp3D
C
      use mba3d_error
      use mba3d_makQ
C
      contains
C
C ================================================================
      Subroutine LINTRP3D( 
C ================================================================
     .           NT, tet, NV, vrt, LDF, F, NXYZ, xyz, G,
     .           imem,Nimem, dmem,Ndmem, iControl )
C ================================================================
C  Subrouine  LINTRP3D  interpolates the peicewise linear  function
C  determined at the nodes of 3D triangulation onto a given set of points.
C
C  Parameters of subroutine in order to locate (without typing):
C       NT         - number of tetrahedra
C       tet(4, NT) - list of tetrahedra ( 4 vertices )
C       NV         - number of nodes (vetices) of  triangulation
C       vrt(3, NV) - coords of vertices
C
C       LDF          - leading dimension of vector function F and G
C       F(LDF, NV)   - the values of vector-function in the nodes of
C                      triangulation
C       NXY          - number of points where the value of function to be
C                      determined
C       xyz(3, NXYZ) - coords of points
C       G(LDF, NXYZ) - the returned values of vector-function in given
C                      points
C
C       imem(Nimem) - auxiliary array of integers of length Nimem
C       dmem(Ndmem) - auxiliary array of double precision numbers of 
C                     length Ndmem
C
C       iContol - 4-digit integer representing 3 control parameters
C            1st digit = 1 means "Initialisation is needed"
C                      = 2 means "Initialisation is not needed"
C                      otherwise the input is erroneous
C            2nd digit  = 1 means  "Nearby a true curved boundary patch 
C                                   the domain is convex" 
C                         2 means  "Nearby a true curved boundary patch 
C                                   the domain is not convex"
C                        otherwise  the input is erroneous
C            3,4 digits = 00 means "No output information is provided"
C                  ab   > 00 means "Important warnings are written to 
C                                   the channel (NUNIT=ab )"
C ================================================================
C  REMARK 1.
C    If the domain is polyhedral, and for a given point there is no an 
C    element containing it within tolerance PREC=10^{-6} , then the ERROR 
C    is assumed. The code could not find an element within PREC tolerance 
C    for MaxTrials.
C
C  REMARK 2.
C    If the domain has a curved boundary, it is possible that a given 
C    point is in the domain but out of the mesh. In this case the tolerance 
C    PREC will be relaxed until an element will be found. Yet, the warning 
C    about the point will be issued if NUNIT > 0. That is why for the first 
C    usages of the code NUNIT > 0 are recommended.
C
C ================================================================
      implicit none
      include 'lintrp.fd'
C ================================================================
      integer NT, NV, LDF, NXYZ, Nimem, Ndmem
      integer tet(4, NT)
      double precision vrt(3, NV), F(LDF, NV)
      double precision xyz(3, NXYZ), G(LDF, NXYZ)
      integer imem(*)
      double precision dmem(*)
      integer iControl
      
      double precision h(MaxH)
      integer i, NQT, QT, REF, ITET
      logical flags(2)
      integer NUNIT

C ================================================================
c ... decoding control parameters
      i = iControl/1000
      if ( i .eq. 1 ) then
         flags(1) = .true.
      else if ( i .eq. 2 ) then
         flags(1) = .false.
      else
         Call errMes(6101, 'lintrp3D',
     .              'Wrong 1st digit in iControl')
      end if 

      i = (iControl-1000*(iControl/1000))/100
      if ( i .eq. 1 ) then
         flags(2) = .true.
      else if ( i .eq. 2 ) then
         flags(2) = .false.
      else
         Call errMes(6102, 'lintrp3D',
     .              'Wrong 2nd digit in iControl')
      end if
      NUNIT = mod(iControl,100)


c ... distributing working memory
      if(flags(1)) then
         QT = 5
         call INITQT( NQT, imem(QT), h, dmem(MaxH + 1) )
         do i = 1, NV
            call DROPS( NQT, imem(QT), Nimem-4, h,dmem(MaxH + 1),vrt,i)
         enddo
         REF  = QT + 8 * NQT
         ITET = REF + NV + 1
      else
         QT  = imem(1)
         NQT = imem(2)
         REF = imem(3)
         ITET =imem(4)
         
         do i = 1, MaxH
            h(i) = dmem(i)
         end do
      end if
      
      if (ITET+4*NT.gt.Nimem) then
         Call errMes(1001, 'lintrp3D',
     .              'not enough memory for Integer arrays')
      end if
      if (MAXH+3*NQT.gt.Ndmem) then
         Call  errMes(1002, 'lintrp3D',
     .        'not enough memory for Real*8 arrays')
      end if


c ... calling the main module 
      call RESTORE( NT, tet, NV, vrt, LDF, F, NXYZ, xyz, G,
     .              imem(QT), dmem(MaxH + 1),
     .              imem(REF), imem(ITET),
     .              h, flags, imem(ITET+4*NT), Nimem-ITET-4*NT,
     .              NUNIT ) 


c ... saving memory distribution
      if(flags(1)) then
         imem(1) = QT
         imem(2) = NQT
         imem(3) = REF
         imem(4) = ITET
         
         do i = 1, MaxH
            dmem(i) = h(i)
         end do
      end if
      return
      end Subroutine LINTRP3D



C ================================================================
      Subroutine INITQT( NQT, QT, h, XYZc )
C ================================================================
C     Initializing of octtree sructure
C ================================================================
      implicit none
      include 'lintrp.fd'
C ================================================================
      integer NQT, QT(8, *)
      double precision XYZc(3, *)
      double precision h(*)

      integer k

C ================================================================
      h(1) = 0.5
      do k = 2, MaxH
         h(k) = h(k-1) / 2
      enddo

      NQT = 0
      NQT = NEWQT( NQT, QT )
      XYZc(1, NQT) = 0.5
      XYZc(2, NQT) = 0.5
      XYZc(3, NQT) = 0.5

      return
      end Subroutine INITQT



C ================================================================
      Integer function NEWQT ( NQT, QT )
C ================================================================
C  The function allocates a new octtree
C ================================================================
      implicit none
C ================================================================
      integer NQT
      integer QT(8, *)

      integer i

C ================================================================
      NQT = NQT + 1
      do i = 1, 8
         QT(i, NQT) = 0
      enddo
      NEWQT = NQT
      return
      end function NEWQT



C ================================================================
      Subroutine SETIJ ( XYZc, XYZ, I, J, K )
C ================================================================
C  Definition of octant of point XYZ with respect to the centre of
C  quadtree XYZc.
C ================================================================
      implicit none
C ================================================================
      double precision XYZc(3), XYZ(3)
      integer I, J, K

      if ( XYZ(1).GT.XYZc(1) ) then
         I = 2
      else
         I = 1
      endif

      if ( XYZ(2).GT.XYZc(2) ) then
         J = 2
      else
         J = 1
      endif

      if ( XYZ(3).GT.XYZc(3) ) then
         K = 2
      else
         K = 1
      endif

      return
      end Subroutine SETIJ



C ================================================================
      Subroutine DROPS ( NQT, QT, NQTav, h, XYZc, XYZ, idx )
C ================================================================
C     The arrangement of a new point in octtree
C ================================================================
      implicit none
      include 'lintrp.fd'
C ================================================================
      integer NQT, idx, NQTav
      integer QT(2, 2, 2, *)
      double precision XYZc(3, *), XYZ(3, *)
      double precision h(*)

      integer L, ip, ptr, i, j, k, new, dir(2)

C ================================================================
      dir(1)=-1
      dir(2)= 1

      ip = 1
      L = 1
      do while ( .TRUE. )
         L = L + 1
         if ( L.GT.MaxH ) Goto 1000
         call SETIJ( XYZc(1, ip), XYZ(1, idx), i, j, k )
         ptr = QT(i, j, k, ip)
         if ( ptr.EQ.0 ) then
            QT(i, j, k, ip) = idx
            return
         else if ( ptr.LT.0 ) then
            ip = -ptr
         else if ( ptr.GT.0 ) then
            new = NEWQT( NQT, QT )
            if (NQT.gt.NQTav) then
               Call errMes(1001, 'drops',
     .                     'not enough memory for Integer arrays')
            end if
            QT(i, j, k, ip) = -new
            XYZc(1, new) = XYZc(1, ip) + dir(i)*h(L)
            XYzc(2, new) = XYZc(2, ip) + dir(j)*h(L)
            XYZc(3, new) = XYZc(3, ip) + dir(k)*h(L)
            call SETIJ( XYZc(1, new), XYZ(1, ptr), i, j, k )
            QT(i, j, k, new) = ptr
            ip = new
         endif
      enddo
      return
 1000 Continue
      Call errMes(1009, 'drops', 'local parameter MaxH is small')
      end Subroutine DROPS



C ================================================================
      Double Precision Function SQRDST( a, b )
C ================================================================
C  Function returns the square of distance between points of a space
C ================================================================
      double precision a(3), b(3)

      SQRDST = (a(1) - b(1))**2 + (a(2) - b(2))**2 + (a(3) - b(3))**2

      return
      end Function SQRDST



C ================================================================
      Subroutine ORDER2( XYZc, h, XYZ, ord, sqrd )
C ================================================================
C  Ordering of subocts of octtree QT by their distanse from
C  the given point XYZ
C ================================================================
      implicit none
C ================================================================
      integer ord(8)
      double precision  XYZc(3), XYZ(3), h, sqrd(8)

      integer i, j, k
      double precision dist(3), sqr(3, 2)

C ================================================================
      integer ofs
      ofs(i, j, k) = i + 2 * (j - 1) + 4 * (k - 1)

C ================================================================
      do i = 1, 3
         dist(i)   = DABS( XYZ(i) - XYZc(i) )
         sqr(i, 1) = dist(i)**2
         sqr(i, 2) = DDIM( dist(i), h )**2
      enddo
      call SETIJ( XYZc, XYZ, i, j, k )

      ord(1) = ofs(i, j, k)
      sqrd(1) = sqr(1, 1) + sqr(2, 1) + sqr(3, 1)
      if ( dist(1).LE.dist(2).and.dist(2).LE.dist(3) ) then
         ord(2) = ofs(3-i, j, k)
         ord(3) = ofs(i, 3-j, k)
         ord(4) = ofs(3-i, 3-j, k)
         ord(5) = ofs(i, j, 3-k)
         ord(6) = ofs(3-i, j, 3-k)
         ord(7) = ofs(i, 3-j, 3-k)
         sqrd(2) = sqr(1, 1) + sqr(2, 2) + sqr(3, 2)
         sqrd(3) = sqr(1, 2) + sqr(2, 1) + sqr(3, 2)
         sqrd(4) = sqr(1, 1) + sqr(2, 1) + sqr(3, 2)
         sqrd(5) = sqr(1, 2) + sqr(2, 2) + sqr(3, 1)
         sqrd(6) = sqr(1, 1) + sqr(2, 2) + sqr(3, 1)
         sqrd(7) = sqr(1, 2) + sqr(2, 1) + sqr(3, 1)
         goto 1
      end if
      if ( dist(2).LE.dist(1).and.dist(1).LE.dist(3) ) then
         ord(2) = ofs(i, 3-j, k)
         ord(3) = ofs(3-i, j, k)
         ord(4) = ofs(3-i, 3-j, k)
         ord(5) = ofs(i, j, 3-k)
         ord(6) = ofs(i, 3-j, 3-k)
         ord(7) = ofs(3-i, j, 3-k)
         sqrd(2) = sqr(1, 2) + sqr(2, 1) + sqr(3, 2)
         sqrd(3) = sqr(1, 1) + sqr(2, 2) + sqr(3, 2)
         sqrd(4) = sqr(1, 1) + sqr(2, 1) + sqr(3, 2)
         sqrd(5) = sqr(1, 2) + sqr(2, 2) + sqr(3, 1)
         sqrd(6) = sqr(1, 2) + sqr(2, 1) + sqr(3, 1)
         sqrd(7) = sqr(1, 1) + sqr(2, 2) + sqr(3, 1)
         goto 1
      end if
      if ( dist(1).LE.dist(3).and.dist(3).LE.dist(2) ) then
         ord(2) = ofs(3-i, j, k)
         ord(3) = ofs(i, j, 3-k)
         ord(4) = ofs(3-i, j, 3-k)
         ord(5) = ofs(i, 3-j, k)
         ord(6) = ofs(3-i, 3-j, k)
         ord(7) = ofs(i, 3-j, 3-k)
         sqrd(2) = sqr(1, 1) + sqr(2, 2) + sqr(3, 2)
         sqrd(3) = sqr(1, 2) + sqr(2, 2) + sqr(3, 1)
         sqrd(4) = sqr(1, 1) + sqr(2, 2) + sqr(3, 1)
         sqrd(5) = sqr(1, 2) + sqr(2, 1) + sqr(3, 2)
         sqrd(6) = sqr(1, 1) + sqr(2, 1) + sqr(3, 2)
         sqrd(7) = sqr(1, 2) + sqr(2, 1) + sqr(3, 1)
         goto 1
      end if
      if ( dist(2).LE.dist(3).and.dist(3).LE.dist(1) ) then
         ord(2) = ofs(i, 3-j, k)
         ord(3) = ofs(i, j, 3-k)
         ord(4) = ofs(i, 3-j, 3-k)
         ord(5) = ofs(3-i, j, k)
         ord(6) = ofs(3-i, 3-j, k)
         ord(7) = ofs(3-i, j, 3-k)
         sqrd(2) = sqr(1, 2) + sqr(2, 1) + sqr(3, 2)
         sqrd(3) = sqr(1, 2) + sqr(2, 2) + sqr(3, 1)
         sqrd(4) = sqr(1, 2) + sqr(2, 1) + sqr(3, 1)
         sqrd(5) = sqr(1, 1) + sqr(2, 2) + sqr(3, 2)
         sqrd(6) = sqr(1, 1) + sqr(2, 1) + sqr(3, 2)
         sqrd(7) = sqr(1, 1) + sqr(2, 2) + sqr(3, 1)
         goto 1
      end if
      if ( dist(3).LE.dist(2).and.dist(2).LE.dist(1) ) then
         ord(2) = ofs(i, j, 3-k)
         ord(3) = ofs(i, 3-j, k)
         ord(4) = ofs(i, 3-j, 3-k)
         ord(5) = ofs(3-i, j, k)
         ord(6) = ofs(3-i, j, 3-k)
         ord(7) = ofs(3-i, 3-j, k)
         sqrd(2) = sqr(1, 2) + sqr(2, 2) + sqr(3, 1)
         sqrd(3) = sqr(1, 2) + sqr(2, 1) + sqr(3, 2)
         sqrd(4) = sqr(1, 2) + sqr(2, 1) + sqr(3, 1)
         sqrd(5) = sqr(1, 1) + sqr(2, 2) + sqr(3, 2)
         sqrd(6) = sqr(1, 1) + sqr(2, 2) + sqr(3, 1)
         sqrd(7) = sqr(1, 1) + sqr(2, 1) + sqr(3, 2)
         goto 1
      end if
      if ( dist(3).LE.dist(1).and.dist(1).LE.dist(2) ) then
         ord(2) = ofs(i, j, 3-k)
         ord(3) = ofs(3-i, j, k)
         ord(4) = ofs(3-i, j, 3-k)
         ord(5) = ofs(i, 3-j, k)
         ord(6) = ofs(i, 3-j, 3-k)
         ord(7) = ofs(3-i, 3-j, k)
         sqrd(2) = sqr(1, 2) + sqr(2, 2) + sqr(3, 1)
         sqrd(3) = sqr(1, 1) + sqr(2, 2) + sqr(3, 2)
         sqrd(4) = sqr(1, 1) + sqr(2, 2) + sqr(3, 1)
         sqrd(5) = sqr(1, 2) + sqr(2, 1) + sqr(3, 2)
         sqrd(6) = sqr(1, 2) + sqr(2, 1) + sqr(3, 1)
         sqrd(7) = sqr(1, 1) + sqr(2, 1) + sqr(3, 2)
         goto 1
      end if
 1    ord(8) = ofs(3-i, 3-j, 3-k)
      sqrd(8) = sqr(1, 2) + sqr(2, 2) + sqr(3, 2)

      return
      end Subroutine ORDER2



C ================================================================
      Integer Function NEARST ( QT, XYZc, xyz, point, h )
C ================================================================
C     This function returns the index of the octtree cell which contains
C     the nearest to given point  grid node
C ================================================================
      implicit none
      include 'lintrp.fd'
C ================================================================
      integer QT(8, *)
      double precision XYZc(3, *), xyz(3, *), point(3)
      double precision h(*)

      double precision sqrd(8, MaxH), sqdist, min
      integer L, ip(MaxH), k(MaxH), ord(8, MaxH), ptr

C ================================================================
      L   = 1
      ptr = -1
      min = 2.0
      NEARST = 0

      do while ( .TRUE. )
         if ( ptr.LT.0 ) then
            L = L + 1
            if ( L.GT.MaxH ) Goto 1000
            ip(L) = -ptr
            call ORDER2( XYZc(1, -ptr), h(L), point, ord(1, L),
     .                   sqrd(1, L) )
            k(L) = 1
         else
            if ( ptr.GT.0 ) then
               sqdist = SQRDST( point, xyz(1, ptr) )
               if ( min.GT.sqdist ) then
                  min = sqdist
                  NEARST = ord(k(L),L) + 8*(ip(L) - 1)
               endif
            endif
            do while ( .TRUE. )
               k(L) = k(L) + 1
               if ( k(L).LE.8 ) then
                  if( sqrd(k(L), L).LT.min ) goto 101
               end if
               L = L - 1
               if ( L.EQ.1 ) goto 102
            enddo
         endif
 101     continue
         ptr = QT(ord(k(L), L), ip(L))
      enddo
 102  continue
      return

 1000 Continue
      Call errMes(1009, 'nearest', 'local parameter MaxH is small')
      end Function NEARST



C ================================================================
      Integer Function BASETET( QT, XYZc, vrt, NT, tet, ref, itet, XYZ,
     &                          h, buf, Nbuf, PREC )
C ================================================================
C  The function  determines the underlying tetrahedron
C  for a point XYZ inside of triangulation
C ================================================================
      implicit none
      include 'lintrp.fd'
C ================================================================
      integer NT
      integer QT(*), tet(4, *), ref(*), itet(*)
      double precision XYZc(3, *), vrt(3, *), XYZ(3)
      double precision h(*), PREC
      integer          buf(2,*)
      integer          Nbuf

      integer i, j, k, n, idx, ip, MaxTrials

C ================================================================
      BASETET = 0
      k = 0
      MaxTrials = max(NT * RelativeTrials, 20D0)

      do while ( k.lt.MaxTrials )
         ip = NEARST( QT, XYZc, vrt, XYZ, h )
         if ( ip.EQ.0 ) goto 101
         idx = QT(ip)
         do i = ref(idx), ref(idx+1)-1
            if ( ENCLOSE( XYZ, vrt, tet(1, itet(i)), PREC ) ) then
               BASETET = itet(i)
               goto 101
            endif
         enddo

c Yuri Vassilevski: instead of searching in a neighborhood, take
c                   another cell of the octtree
C ...  checking for neighbooring tetrahedra
c         do i = ref(idx), ref(idx+1)-1
c           iT = itet(i)
c           Do j = 1, 4
c              idx2 = tet(j, iT)
c              Do m = ref(idx2), ref(idx2+1)-1
c                 if ( ENCLOSE(XYZ,vrt,tet(1,itet(m)),PREC) ) then
c                    BASETET = itet(m)
c                    goto 101
c                 endif
c              End do
c           End do
c         enddo
C ...  end of


         k = k + 1
         if ( 2*k.GT.Nbuf ) goto 1000
         buf(1, k) = ip
         buf(2, k) = idx
         QT(ip) = 0
      enddo
 101  do n = k, 1, -1
         QT(buf(1, n)) = buf(2, n)
      enddo
      return

1000  Continue
      Call errMes(1001, 'basetet',
     .           'not enough memory for Integer arrays')
      end Function BASETET



C ================================================================
      Subroutine RESTORE(NT, tet, NV, vrt,  LDF, F,
     .                   NXYZ,  xyz, G, QT, XYZc, ref, itet,
     .                   h, flags, buf, Nbuf, NUNIT )
C ================================================================
C  Code builds the transposed list of tetrahedra grid and restores
C  the values of vector-function F in the points of interest with the aid
C  of the octtree structure
C ================================================================
      implicit none
      include 'lintrp.fd'
C ================================================================
c Predefined tolerance for checking whether the point belongs to an element
      double precision PREC
      parameter( PREC = 1D-6 )

c Predefined number of possible relaxations of the above tolerance (6 orders of 10, here)
      integer    MaxRelax
      parameter( MaxRelax = 7 )

C ================================================================
      integer NV, NT, LDF, NXYZ
      integer QT(2, 2, 2, *), tet(4, NT), ref(*), itet(*)
      double precision XYZc(3, *)
      double precision vrt(3, NV), xyz(3, NXYZ)
      double precision F(LDF, NV), G(LDF, NXYZ)
      double precision h(*)
      logical          flags(2)
      integer          NUNIT
      integer          buf(*)
      integer          Nbuf

C ================================================================
C (Local variables)
      integer          i, j, idx, irel, ip(4)
      double precision a, b, c, d, det, PREC1

C ================================================================
      if(flags(1)) then
         do i = 1, NV+1
            ref(i) = 0
         enddo
         do i = 1, NT
            do j = 1, 4
               idx = tet(j, i)
               ref(idx) = ref(idx) + 1
            enddo
         enddo
         ref(1) = ref(1) + 1
         do i = 2, NV+1
            ref(i) = ref(i-1) + ref(i)
         enddo
         do i = 1, NT
            do j = 1, 4
               idx = tet(j, i)
               ref(idx) = ref(idx) - 1
               itet(ref(idx)) = i
            enddo
         enddo
      end if
      
      Do 10 i = 1, NXYZ
         PREC1 = PREC
         irel = 0
 11      idx = BASETET( QT, XYZc, vrt, NT, tet, ref, itet, xyz(1, i),
     &                  h,  buf, Nbuf, PREC1 )
         if(idx.le.0) then
C if Curvelinear boundaries exist, then probably the node is out of the mesh
C Relaxation of PREC is applied
            if (flags(2)) then
               if (PREC1.lt.5D-2) then
                  PREC1 = PREC1*10
               else
                  PREC1 = PREC1*3.16
               end if
               
               if (irel.lt.MaxRelax) then
                  irel = irel + 1
                  Goto 11
               else
                  Call errMes(6103, 'lintrp3D',
     .                 'Max number of relaxations for PREC is reached')
               end if
C  if no Curvelinear boundaries, then error is assumed or inconsistency
C of parameters in lintrp.fd. REINSTALL constants in lintrp.fd!
            else
               Call errMes(6104, 'lintrp3D',
     .             'Failed to find element within PREC tolerance')
            end if
         end if
         if (PREC1.ne.PREC) then
            if (NUNIT.GT.0) then
               write(NUNIT,21) PREC1
               write(NUNIT,22) xyz(1,i),xyz(2,i),xyz(3,i)
               write(NUNIT,23) PREC
               write(NUNIT,24)
            end if
         end if
         
         do j = 1, 4
            ip(j) = tet(j, idx)
         enddo
         a = calVol(xyz(1,i), vrt(1,ip(2)), vrt(1,ip(3)), vrt(1,ip(4)))
         b = calVol(xyz(1,i), vrt(1,ip(1)), vrt(1,ip(4)), vrt(1,ip(3)))
         c = calVol(xyz(1,i), vrt(1,ip(1)), vrt(1,ip(2)), vrt(1,ip(4)))
         d = calVol(xyz(1,i), vrt(1,ip(1)), vrt(1,ip(3)), vrt(1,ip(2)))
         det = a + b + c + d
         do j = 1, LDF
            G(j, i) = (a * F(j, ip(1)) + b * F(j, ip(2)) +
     .                 c * F(j, ip(3)) + d * F(j, ip(4))) / det
         end do
 10   Continue  

 21   Format('Caution! PREC1 was relaxed: ', F11.7, ' since')
 22   Format('the point ',3F9.5,' is probably out of the mesh')
 23   Format('After MaxTrials failed to find a host within ',F8.6,
     &       ' tolerance') 
 24   Format('NEED TO REINSTALL constants in lintrp.fd')

      return
      end Subroutine RESTORE




C ================================================================
      logical function ENCLOSE ( XYZ, vrt, TET, PREC )
C ================================================================
C  Determines does the point XYZ belong to the tetrahedra
C ================================================================
      implicit none
C ================================================================
      double precision XYZ(3), vrt(3, *), PREC
      integer TET(4)

      integer          Face(3,4)
      double precision frac, Vlm
      integer i

      data Face/2,3,4, 1,4,3, 1,2,4, 1,3,2/

C ================================================================
      Vlm = calVol( vrt(1,TET(1)), vrt(1,TET(Face(1,1))),
     &              vrt(1,TET(Face(2,1))), vrt(1,TET(Face(3,1))))

      do i = 1, 4
         frac = calVol( XYZ(1), vrt(1,TET(Face(1,i))),
     &                          vrt(1,TET(Face(2,i))), 
     &                          vrt(1,TET(Face(3,i)))) /  Vlm
         if ( frac.LE.-PREC ) then
            ENCLOSE = .FALSE.
            return
         endif
      enddo
      ENCLOSE = .TRUE.

      return
      end function ENCLOSE



C ================================================================
      Double Precision Function sizeQT(point, imem, dmem)
C ================================================================
      include 'lintrp.fd'
C ================================================================
C  The function returns the size of the octtree cell which contains
C  the given point with coords point(3). 
C
C  PARAMETERS:
C     point(3) - user given point inside the unit cube (0,1)^3
C
C     dmem(*)  - real*8  working memory of LINTRP3D
C     imem(*)  - integer working memory of LINTRP3D
C
C  REMARK 1. This function is the wrapping for function SizeHost.
C
C ================================================================
      Double Precision point(3), dmem(*)
      Integer          imem(*)

      iQT = imem(1)
      sizeQT = SizeHost(imem(iQT), dmem(MaxH + 1), point, dmem)

      Return
      End Function sizeQT



C ================================================================
      Double Precision Function SizeHost ( QT, XYZc, point, h )
C ================================================================
      implicit none
      include 'lintrp.fd'
C ================================================================
      integer QT(8, *)
      double precision XYZc(3, *), point(3)
      double precision h(*)

      double precision sqrd(8, MaxH)
      integer L, ip(MaxH), ord(8, MaxH), ptr

C ================================================================
      L   = 1
      ptr = -1
      SizeHost = 1D0

      do while ( .TRUE. )
         if ( ptr.LT.0 ) then
            L = L + 1
            if ( L.GT.MaxH ) Goto 1000

            ip(L) = -ptr
            call ORDER2( XYZc(1, -ptr), h(L), point, ord(1, L),
     .                   sqrd(1, L) )
         else
            SizeHost = h(L)
            return
         endif
         ptr = QT(ord(1, L), ip(L))
      enddo

      Call errMes(6105, 'SizeHost', 'host cell was not found')
      return

 1000 Continue
      Call errMes(1009, 'SizeHost', 'local parameter MaxH is small')
      end Function SizeHost
C
      End Module mba3d_lintrp3D
