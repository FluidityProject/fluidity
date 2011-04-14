C Copyright (C) 2006 Imperial College London and others.
C 
C Please see the AUTHORS file in the main source directory for a full list
C of copyright holders.
C 
C Adrian Umpleby
C Applied Modelling and Computation Group
C Department of Earth Science and Engineering
C Imperial College London
C 
C adrian@Imperial.ac.uk
C 
C This library is free software; you can redistribute it and/or
C modify it under the terms of the GNU Lesser General Public
C License as published by the Free Software Foundation; either
C version 2.1 of the License.
C 
C This library is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C Lesser General Public License for more details.
C 
C You should have received a copy of the GNU Lesser General Public
C License along with this library; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
C USA
      SUBROUTINE NODTST( BIGLST, NODLST,
     :                   ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                   SZENLS, SZNELS, NNOD,   NELM,
     :                   NODX,   NODY,   NODZ,   ORGMTX,
     :                   INOD,   IELM,
     :                   FNCDIF, FNCORG, XN, YN, ZN, BSTELM,
     :                   AVEDIF, AVEORG, XA, YA, ZA, BSTELA,
     :                   ALWAYS, DOTOP, MINCHG )
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER INOD, IELM, BSTELM, BSTELA
C
      REAL FNCDIF, FNCORG, AVEDIF, AVEORG, XN, YN, ZN, XA, YA, ZA,
     :     ALWAYS, DOTOP, MINCHG
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkorg.i'
C
      INCLUDE 'blkerr.i'
C
      INCLUDE 'blktim.i'
C
      INTEGER ELPTRS(MAXEDE), EDPTRS(MAXEDE), NDPTRS(MAXEDE), I, J, K,
     :        NUMEDE, NUMELM, NODS(4), INELM, INEW, JNOD, NXT, LST,
     :        ELNDS(4,MAXEDE), ORGELM, ifail, ii, jj, ki, kj,
     :        NTST, ITST, cnt
C
      REAL F0, F1, F2, F3, F4, F5, F6, FT, T1, T2, T3, T4, T5, T6,
     :     X(MAXEDE), Y(MAXEDE), Z(MAXEDE), sx, sy, sz, fadd,
     :     D1, D2, D3, V1(3), V2(3), V3(3), XO, YO, ZO, V(MAXEDE),
     :     XX(10), YY(10), ZZ(10), VOL1, VOL2, VOL3,
     :     EDGFNC, ELMFNC, UNDLEN, RADS(MAXEDE), ORGM(9), A(9),
     :     r0, r1, r2, r3, r4, r5, r6, r7, r8, r9,
     :     x1, y1, z1, x2, y2, z2, dd(maxede), xp(4), yp(4), zp(4),
     :     vol(maxede), tetvol, det
C
      LOGICAL ISINSD, ISSPLT, GTNDFL, GETFLG, LSTEDG, goterr
C
      CALL ADAPT_SECOND( TIMSTT )
C
      FNCORG = 1E+31
      FNCDIF = 1E+30
      BSTELM = 0
      AVEORG = 1E+31
      AVEDIF = 1E+30
      BSTELA = 0
C
c      goterr = ( DEBUG .and. (inod .eq. 573) )
c      goterr = debug .or. ( anisot .eq. 2 )
      goterr = debug
c      goterr = .true.
c      goterr = ( inod .eq. 106 )
C
cc      IF( GTNDFL(INOD,2) ) THEN
cc         PRINT*,'+++ Warning: attempt to relocate geometry node'
cc         return
cc      end if
C
      CALL NDCNCT( BIGLST, NODLST, INOD, IELM,
     :             ELPTRS, EDPTRS, NDPTRS, NUMEDE, NUMELM )
C
      IF( IERR .NE. 0 ) RETURN
C
c      goterr = (inod.eq.4621) .and. (numelm.eq.20) .and. (numede.eq.16)
c      if( inod.eq.30553 .and. numelm.eq.22 .and. numede.eq.13 )
c     :    goterr = .true.
      if( goterr ) then
         print*,' '
         print*,'Into NODTST... ',inod
         print*,'els,eds: ',numelm,numede
         print*,'els: ',(elptrs(i),i=1,numelm)
         print*,'eds: ',(edptrs(i),i=1,numede)
         print*,'nds: ',(ndptrs(i),i=1,numede)
         print*,' '
       end if
c      LSTEDG = EDGON
C
c      EDGON = .TRUE.
C
c      CALL NODEIG( BIGLST, NODLST, INOD, V1, V2, V3, D1, D2, D3 )
C
      XO = NODLST( 1, INOD )
      YO = NODLST( 2, INOD )
      ZO = NODLST( 3, INOD )
      DO I = 1, 9
         ORGM(I) = NODLST(I+6,INOD)
      END DO
      ORGELM = INT(NODLST(16,INOD))
      INELM  = ORGELM
C
      XN = XO
      YN = YO
      ZN = ZO
c
      XA = XO
      YA = YO
      ZA = ZO
c
      if( goterr ) print*,'   Old pos: ',xo,yo,zo
C
      XX(1) = 0.0
      YY(1) = 0.0
      ZZ(1) = 0.0
C
      ISSPLT = GTNDFL( BIGLST, NODLST, INOD, 4 )
      ISINSD = GTNDFL( BIGLST, NODLST, INOD, 3 ) .AND. .NOT. ISSPLT
C
      I = 0
c
      do k = 1, 9
         a(k) = 0.0
      end do
C
      sx = 0.0
      sy = 0.0
      sz = 0.0
c
      IF( ISINSD ) THEN
C
         DO J = 1, NUMEDE
            IF( GETFLG(BIGLST,NODLST,EDPTRS(J),6) ) RETURN
            X(J) = NODLST(1,NDPTRS(J)) - XO
            Y(J) = NODLST(2,NDPTRS(J)) - YO
            Z(J) = NODLST(3,NDPTRS(J)) - ZO
            t1 = edgfnc( BIGLST, NODLST, ndptrs(j), inod, v(j) )
            do k = 1, 9
               a(k) = a(k) + nodlst(k+6,ndptrs(j))
            end do
            sx = sx + nodlst( 7,ndptrs(j))*x(j)
     :              + nodlst(10,ndptrs(j))*y(j)
     :              + nodlst(13,ndptrs(j))*z(j)
            sy = sy + nodlst( 8,ndptrs(j))*x(j)
     :              + nodlst(11,ndptrs(j))*y(j)
     :              + nodlst(14,ndptrs(j))*z(j)
            sz = sz + nodlst( 9,ndptrs(j))*x(j)
     :              + nodlst(12,ndptrs(j))*y(j)
     :              + nodlst(15,ndptrs(j))*z(j)
         END DO
C
         I = NUMEDE
C
      ELSE IF( ISSPLT ) THEN
C
         DO J = 1, NUMEDE
            IF( GETFLG(BIGLST,NODLST,EDPTRS(J),6) ) RETURN
            IF( GETFLG(BIGLST,NODLST,EDPTRS(J),2) ) THEN
               I = I + 1
               X(I) = NODLST(1,NDPTRS(J)) - XO
               Y(I) = NODLST(2,NDPTRS(J)) - YO
               Z(I) = NODLST(3,NDPTRS(J)) - ZO
               t1 = edgfnc( BIGLST, NODLST, ndptrs(j), inod, v(i) )
               do k = 1, 9
                  a(k) = a(k) + nodlst(k+6,ndptrs(j))
               end do
               sx = sx + nodlst( 7,ndptrs(j))*x(i)
     :                 + nodlst(10,ndptrs(j))*y(i)
     :                 + nodlst(13,ndptrs(j))*z(i)
               sy = sy + nodlst( 8,ndptrs(j))*x(i)
     :                 + nodlst(11,ndptrs(j))*y(i)
     :                 + nodlst(14,ndptrs(j))*z(i)
               sz = sz + nodlst( 9,ndptrs(j))*x(i)
     :                 + nodlst(12,ndptrs(j))*y(i)
     :                 + nodlst(15,ndptrs(j))*z(i)
            END IF
         END DO
C
          IF( I .NE. 2 ) THEN
            PRINT*,'*** NODTST: GEOM-SPLIT NODE NOT ON TWO GEOM EDGES'
            IERR = 601
            RETURN
         END IF
C
      ELSE
C
         DO J = 1, NUMEDE
            IF( GETFLG(BIGLST,NODLST,EDPTRS(J),6) ) RETURN
            IF( .NOT. GETFLG(BIGLST,NODLST,EDPTRS(J),3) ) THEN
               I = I + 1
               X(I) = NODLST(1,NDPTRS(J)) - XO
               Y(I) = NODLST(2,NDPTRS(J)) - YO
               Z(I) = NODLST(3,NDPTRS(J)) - ZO
               t1 = edgfnc( BIGLST, NODLST, ndptrs(j), inod, v(i) )
               do k = 1, 9
                  a(k) = a(k) + nodlst(k+6,ndptrs(j))
               end do
               sx = sx + nodlst( 7,ndptrs(j))*x(i)
     :                 + nodlst(10,ndptrs(j))*y(i)
     :                 + nodlst(13,ndptrs(j))*z(i)
               sy = sy + nodlst( 8,ndptrs(j))*x(i)
     :                 + nodlst(11,ndptrs(j))*y(i)
     :                 + nodlst(14,ndptrs(j))*z(i)
               sz = sz + nodlst( 9,ndptrs(j))*x(i)
     :                 + nodlst(12,ndptrs(j))*y(i)
     :                 + nodlst(15,ndptrs(j))*z(i)
            END IF
         END DO
C
         IF( I .LT. 3 ) THEN
            PRINT*,'*** NODTST: SURF NODE ON LESS THAN 3 SURF EDGES'
            IERR = 601
            RETURN
         END IF
C
      END IF
C
c      EDGON = LSTEDG
c
      fadd = 0.0
c
c      do k = 1, 9
c         a(k) = a(k) + nodlst(k+6,inod)*i
c         fadd = max( fadd, abs(a(k)) )
c      end do
c      fadd = fadd/2
c      a(1) = a(1) + fadd
c      a(5) = a(5) + fadd
c      a(9) = a(9) + fadd
      fadd = abs(a(2)) + abs(a(3))
      a(1) = max(a(1),fadd*1.01)
      fadd = abs(a(4)) + abs(a(6))
      a(5) = max(a(5),fadd*1.01)
      fadd = abs(a(7)) + abs(a(8))
      a(9) = max(a(9),fadd*1.01)
c
      do j = 1, i
         sx = sx + nodlst( 7,inod)*x(j)
     :           + nodlst(10,inod)*y(j)
     :           + nodlst(13,inod)*z(j)
         sy = sy + nodlst( 8,inod)*x(j)
     :           + nodlst(11,inod)*y(j)
     :           + nodlst(14,inod)*z(j)
         sz = sz + nodlst( 9,inod)*x(j)
     :           + nodlst(12,inod)*y(j)
     :           + nodlst(15,inod)*z(j)
      end do
c
      if( goterr ) then
         print*,'   Matrix: '
         print*,a(1),a(4),a(7)
         print*,a(2),a(5),a(8)
         print*,a(3),a(6),a(9)
         print*,'   Vector: '
         print*,sx,sy,sz
      end if
c
      if( i .eq. numede ) then
c
         call invrs3( a, det, ifail )
c
      else if( i .eq. 2 ) then
c
         r1 = a(1)*x(1) + a(4)*y(1) + a(7)*z(1)
         r2 = a(2)*x(1) + a(5)*y(1) + a(8)*z(1)
         r3 = a(3)*x(1) + a(6)*y(1) + a(9)*z(1)
         r4 = x(1)*r1 + y(1)*r2 + z(1)*r3
         r1 = x(1)*sx + y(1)*sy + z(1)*sz
         if( goterr ) print*,'   r1,r4: ',r1,r4
c
      else
c
         do ii = 1, i
            dd(ii) = 1.0/sqrt( x(ii)*x(ii) + y(ii)*y(ii) + z(ii)*z(ii) )
         end do
c
         ii = 0
c
         fadd = 1.0
c
  10     if( ii .lt. i-1 ) then
            ii = ii + 1
            r1 = dd(ii)
            jj = ii
  20        if( jj .lt. i ) then
               jj = jj + 1
               r2 = dd(jj)
               r3 = x(ii)*x(jj) + y(ii)*y(jj) + z(ii)*z(jj)
               r3 = r3*r1*r2
               if( abs(r3) .lt. abs(fadd) ) then
                  fadd = r3
                  ki   = ii
                  kj   = jj
                  if( abs(fadd) .lt. 0.6 ) goto 30
               end if
               goto 20
            end if
            goto 10
         end if
c
  30     continue
c
         x1 = x(ki)*dd(ki)
         y1 = y(ki)*dd(ki)
         z1 = z(ki)*dd(ki)
c
         x2 = x(kj)*dd(kj) - fadd*x1
         y2 = y(kj)*dd(kj) - fadd*y1
         z2 = z(kj)*dd(kj) - fadd*z1
c
         r1 = 1.0/sqrt(x2*x2 + y2*y2 + z2*z2)
c
         x2 = x2*r1
         y2 = y2*r1
         z2 = z2*r1
c
         if( goterr ) then
            print*,'   Chosen ',ki,kj,fadd
            print*,x1,y1,z1
            print*,x2,y2,z2
            print*,x1*x2 + y1*y2 + z1*z2
         end if
c
         r1 = a(1)*x1 + a(4)*y1 + a(7)*z1
         r2 = a(2)*x1 + a(5)*y1 + a(8)*z1
         r3 = a(3)*x1 + a(6)*y1 + a(9)*z1
         r4 = r1*x1 + r2*y1 + r3*z1
         r5 = r1*x2 + r2*y2 + r3*z2
         r1 = a(1)*x2 + a(4)*y2 + a(7)*z2
         r2 = a(2)*x2 + a(5)*y2 + a(8)*z2
         r3 = a(3)*x2 + a(6)*y2 + a(9)*z2
         r6 = r1*x2 + r2*y2 + r3*z2
c
c         print*,'Original matrix:'
c         print*,a(1),a(4),a(7)
c         print*,a(2),a(5),a(8)
c         print*,a(3),a(6),a(9)
c         print*,'pos1: ',x(1),y(1),z(1)
c         print*,'pos2: ',x(2),y(2),z(2)
c         print*,'About to invert...',inod,kk
c         print*,r4,r5,r6
         a(1) = r4
         a(2) = r5
         a(3) = r5
         a(4) = r6
         r1 = x1*sx + y1*sy + z1*sz
         r2 = x2*sx + y2*sy + z2*sz
         if( goterr ) then
            print*,'   Constrained matrix: '
            print*,a(1),a(3)
            print*,a(2),a(4)
            print*,'   Constrained vector: '
            print*,r1,r2
         end if
         call invrs2( a, det, ifail )
c         print*,'Done.'
c
      end if
C
c      if( .not. tstdst ) then
c         do j = 1, i
c            v(j) = 1.0
c         end do
c      end if
c
c      T1 = 0.0
c      T2 = 0.0
c      T3 = 0.0
c      T4 = 0.0
c      T5 = 0.0
c      T6 = 0.0
cC
c      DO J = 1, 10
c         XX(J) = 0.0
c         YY(J) = 0.0
c         ZZ(J) = 0.0
c      END DO
C
c      DO J = 1, I
cC
c         XX(1) = XX(1) + X(J)
c         YY(1) = YY(1) + Y(J)
c         ZZ(1) = ZZ(1) + Z(J)
cC
c         XX(2) = XX(2) + X(J)*V(J)
c         YY(2) = YY(2) + Y(J)*V(J)
c         ZZ(2) = ZZ(2) + Z(J)*V(J)
cC
c         XX(3) = XX(3) + X(J)*V(J)*V(J)
c         YY(3) = YY(3) + Y(J)*V(J)*V(J)
c         ZZ(3) = ZZ(3) + Z(J)*V(J)*V(J)
cC
c         XX(4) = XX(4) + X(J)*V(J)*V(J)*V(J)*V(J)
c         YY(4) = YY(4) + Y(J)*V(J)*V(J)*V(J)*V(J)
c         ZZ(4) = ZZ(4) + Z(J)*V(J)*V(J)*V(J)*V(J)
cC
c         XX(5) = XX(5) + X(J)*V(J)*V(J)*V(J)*V(J)*V(J)*V(J)
c         YY(5) = YY(5) + Y(J)*V(J)*V(J)*V(J)*V(J)*V(J)*V(J)
c         ZZ(5) = ZZ(5) + Z(J)*V(J)*V(J)*V(J)*V(J)*V(J)*V(J)
cC
c         XX(6) = XX(6) + X(J)*V(J)*V(J)*V(J)*V(J)*V(J)*V(J)*V(J)*V(J)
c         YY(6) = YY(6) + Y(J)*V(J)*V(J)*V(J)*V(J)*V(J)*V(J)*V(J)*V(J)
c         ZZ(6) = ZZ(6) + Z(J)*V(J)*V(J)*V(J)*V(J)*V(J)*V(J)*V(J)*V(J)
cC
c         T1 = T1 + 1.0
c         T2 = T2 + V(J)
c         T3 = T3 + V(J)*V(J)
c         T4 = T4 + V(J)*V(J)*V(J)*V(J)
c         T5 = T5 + V(J)*V(J)*V(J)*V(J)*V(J)*V(J)
c         T6 = T6 + V(J)*V(J)*V(J)*V(J)*V(J)*V(J)*V(J)*V(J)
cC
c      END DO
cC
c      XX(1) = XX(1)/T1 + XO
c      YY(1) = YY(1)/T1 + YO
c      ZZ(1) = ZZ(1)/T1 + ZO
cC
c      XX(2) = XX(2)/T2 + XO
c      YY(2) = YY(2)/T2 + YO
c      ZZ(2) = ZZ(2)/T2 + ZO
cC
c      XX(3) = XX(3)/T3 + XO
c      YY(3) = YY(3)/T3 + YO
c      ZZ(3) = ZZ(3)/T3 + ZO
cC
c      XX(4) = XX(4)/T4 + XO
c      YY(4) = YY(4)/T4 + YO
c      ZZ(4) = ZZ(4)/T4 + ZO
cC
c      XX(5) = XX(5)/T5 + XO
c      YY(5) = YY(5)/T5 + YO
c      ZZ(5) = ZZ(5)/T5 + ZO
cC
c      XX(6) = XX(6)/T6 + XO
c      YY(6) = YY(6)/T6 + YO
c      ZZ(6) = ZZ(6)/T6 + ZO
C
c      xx(4) = ( xx(4) + xx(3) + xx(2) ) / 3
c      yy(4) = ( yy(4) + yy(3) + yy(2) ) / 3
c      zz(4) = ( zz(4) + zz(3) + zz(2) ) / 3
c
      if( i .eq. numede ) then
         xx(1) = (sx*a(1) + sy*a(4) + sz*a(7))*1.5 + xo
         yy(1) = (sx*a(2) + sy*a(5) + sz*a(8))*1.5 + yo
         zz(1) = (sx*a(3) + sy*a(6) + sz*a(9))*1.5 + zo
      else if( i .eq. 2 ) then
         r1 = r1/r4
         xx(1) = x(1)*r1*1.5 + xo
         yy(1) = y(1)*r1*1.5 + yo
         zz(1) = z(1)*r1*1.5 + zo
      else
         r3 = r1*a(1) + r2*a(3)
         r4 = r1*a(2) + r2*a(4)
         xx(1) = (x1*r3 + x2*r4)*1.5 + xo
         yy(1) = (y1*r3 + y2*r4)*1.5 + yo
         zz(1) = (z1*r3 + z2*r4)*1.5 + zo
      end if
c
      if( goterr ) print*,'   New position: ',xx(1),yy(1),zz(1)
C
      NTST = 1
c
c      INEW = -1
c      CALL ADDND( BIGLST, NODLST, XN, YN, ZN, JNOD, INEW )
C
c      IF( IERR .NE. 0 ) RETURN
c      if( goterr ) print*,'new dif: ',xn-xo,yn-yo,zn-zo
c      if( abs(xn-xo) .lt. 1e-5 .and.
c     :    abs(yn-yo) .lt. 1e-5 .and.
c     :    abs(zn-zo) .lt. 1e-5 ) goto 999
C
c      if( i .eq. numede .or. i .eq. 2 ) then
         fadd = 1e+20
c      else
c         fadd = 0.0
c      end if
c
      AVEORG = 0.0
      F0 = 0.0
C
      DO I = 1, NUMELM
C
         CALL ELNODS( BIGLST, NODLST, ELPTRS(I), ELNDS(1,I), NXT, LST,
     :               .FALSE. )
C
         FT = ELMFNC( BIGLST, NODLST, ELPTRS(I), ELNDS(1,I), ELNDS(2,I),
     :                           ELNDS(3,I), ELNDS(4,I), RADS(I) )
         if( goterr ) print*,'  org elm,rad,fn: ',i,rads(i),ft
         AVEORG = AVEORG + FT
         F0 = MAX( F0, FT )
c
         do j = 1, 4
            xp(j) = nodlst(1,elnds(j,i))
            yp(j) = nodlst(2,elnds(j,i))
            zp(j) = nodlst(3,elnds(j,i))
         end do
c
         vol(i) = tetvol( xp, yp, zp )
C
      END DO
C
      FNCORG = F0
      AVEORG = AVEORG / NUMELM
C
      IF( FNCORG .LE. DOTOP ) GOTO 600
C
      DO 500 ITST = 1, NTST
C
        cnt = 1
c
 100    xx(itst) = xx(itst)*0.6 + xo*0.4
        yy(itst) = yy(itst)*0.6 + yo*0.4
        zz(itst) = zz(itst)*0.6 + zo*0.4
c
        nodlst(1,inod) = xx(ITST)
        nodlst(2,inod) = yy(ITST)
        nodlst(3,inod) = zz(ITST)
c
        do i = 1, numelm
c
          do j = 1, 4
            xp(j) = nodlst(1,elnds(j,i))
            yp(j) = nodlst(2,elnds(j,i))
            zp(j) = nodlst(3,elnds(j,i))
          end do
c
          if( tetvol( xp, yp, zp )/vol(i) .le. 1e-4 ) then
            cnt = cnt + 1
            if( cnt .lt. 6 ) goto 100
            goto 500
          end if
c
        end do
c
        inelm = abs(inelm)
c
        CALL FNDELM( ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :               SZENLS, SZNELS, NNOD,   NELM,
     :               NODX,   NODY,   NODZ,   ORGMTX,
     :               XX(ITST),  YY(ITST),  ZZ(ITST),
     :               INELM,  .TRUE., NODLST(7,INOD) )

        if( goterr ) print*,'Done FNDELM: ',INELM
        if( inelm .lt. 0 ) then
           print*,'---+++ NODTST: Got node outside mesh +++---'
           if( issplt ) then
              print*,'   A Splitter node'
           else if( .not. isinsd ) then
              print*,'   A surface node'
           else
              print*,'   An internal node'
           end if
           print*,xo,yo,zo
           print*,xx(itst),yy(itst),zz(itst)
        end if
c
        nodlst(1,inod)  = xx(ITST)
        nodlst(2,inod)  = yy(ITST)
        nodlst(3,inod)  = zz(ITST)
        NODLST(16,INOD) = FLOAT(ABS(INELM))
C        call fixeig( nodlst(1,inod), nodlst(7,inod), anisot )
C
        if( goterr ) print*,'New element functionals...'
C
        F1 = 0.0
        F2 = 0.0
C
        DO I = 1, NUMELM
c
           FT = ELMFNC( BIGLST, NODLST, 0, ELNDS(1,I), ELNDS(2,I),
     :                  ELNDS(3,I), ELNDS(4,I), VOL2 )
           if( goterr ) print*,'  new elm,rad,fn: ',i,vol2,ft
           F2 = F2 + FT
           F1 = MAX( F1, FT )
C
           IF( VOL2/RADS(I) .LE. 1E-4 ) then
              cnt = cnt + 1
              if( cnt .lt. 6 ) goto 100
              GOTO 500
           end if
C
        END DO
C
        if( goterr ) print*,'Orig,new fn: ',f0,f1
        if( goterr ) print*,'Orig,new av: ',aveorg,f2/numelm
        F1 = F1 - FNCORG
        F2 = F2 / NUMELM - AVEORG
C
        IF( F1 .LT. FNCDIF ) THEN
           FNCDIF = F1
           XN = XX(ITST)
           YN = YY(ITST)
           ZN = ZZ(ITST)
           BSTELM = INELM
           if( INELM .lt. 0 .and. fncdif .lt. 0.0 ) then
              print*,'*** NODTST: Node outside mesh is better!'
              print*,fncorg,aveorg,fncdif,min(f2,avedif)
           end if
        END IF
C
        IF( F2 .LT. AVEDIF .AND. F1 .LT. 0.0 ) THEN
           AVEDIF = F2
           XA = XX(ITST)
           YA = YY(ITST)
           ZA = ZZ(ITST)
           BSTELA = INELM
           if( INELM .lt. 0 .and. fncdif .lt. 0.0 ) then
              print*,'*** NODTST: Node outside mesh is better!'
              print*,fncorg,aveorg,fncdif,avedif
           end if
        END IF
C
        IF( FNCDIF .LT. ALWAYS ) GOTO 600
        IF( AVEDIF .LT. ALWAYS ) GOTO 600
c
        cnt = cnt + 1
        if( cnt .lt. 6 ) goto 100
c
 500  CONTINUE
C
 600  nodlst(1,inod) = XO
      nodlst(2,inod) = YO
      nodlst(3,inod) = ZO
      DO I = 1, 9
         NODLST(I+6,INOD) = ORGM(I)
      END DO
      NODLST(16,INOD) = FLOAT(ORGELM)
c
c      CALL FNDELM( ENLBAS, ENLIST, SZENLS, NELBAS, NELIST, SZNELS,
c     :             NODX,   NODY,   NODZ,   ORGMTX, NNOD,   NELM,
c     :             xo,  yo,  zo,  ORGELM,  NODLST(7,INOD) )
C      call fixeig( nodlst(1,inod), nodlst(7,inod), anisot )
C
c      END DO
C
c      DO I = 1, NUMEDE
C
c         F0 = MAX( F1, EDGFNC( INOD, NDPTRS(I), UNDLEN ) )
c         F2 = MAX( F3, EDGFNC( JNOD, NDPTRS(I), UNDLEN ) )
C
c      END DO
C
c      FNCORG = MAX( F0, F1 )
c      FNCDIF = MAX( F2, F3 ) - FNCORG
c
c      if( tstdst ) fncdif = -1e+20
C
 999  continue
c      CALL RMVND( JNOD )
C
      CALL ADAPT_SECOND( TIMEND )
C
      TIMNOD = TIMNOD + TIMEND - TIMSTT
      TOTNOD = TOTNOD + 1
C
      if( goterr ) print*,'Leaving NODTST: ',FNCORG,FNCDIF
      RETURN
      END
C
