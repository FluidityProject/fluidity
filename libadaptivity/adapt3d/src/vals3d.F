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
      SUBROUTINE VALS3D( X, Y, Z, NPNTS, XP, YP, ZP, VALS )
C
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C - Note that both a quadrilateral surface and a tetrahedron have four points
C - (wedge has six, hexahedron has eight, triangular surface has three).
C - In order for this subroutine to distinguish between them a tetrahedron
C - should be sent down as having 5 points - only the first four are used.
C
C-----------------------------------------------------------------------
C
      INTEGER NPNTS
      REAL X(NPNTS), Y(NPNTS), Z(NPNTS), VALS(NPNTS), XP, YP, ZP
C
C - This is all rather silly...
C
      REAL X12, Y12, Z12, X13, Y13, Z13, X14, Y14, Z14, X15, Y15, Z15,
     :     X16, Y16, Z16, X17, Y17, Z17, X18, Y18, Z18, X23, Y23, Z23,
     :     X24, Y24, Z24, X34, Y34, Z34, X87, Y87, Z87, X86, Y86, Z86,
     :     X85, Y85, Z85, X84, Y84, Z84, X83, Y83, Z83, X82, Y82, Z82,
     :     X64, Y64, Z64, X65, Y65, Z65, X31, Y31, Z31, X62, Y62, Z62,
     :     X63, Y63, Z63
      REAL XP1, XP2, XP3, YP1, YP2, YP3, ZP1, ZP2, ZP3, XP4, YP4, ZP4
      REAL X1P, Y1P, Z1P, X2P, Y2P, Z2P, X6P, Y6P, Z6P, X8P, Y8P, Z8P
      REAL XNRM, YNRM, ZNRM, XN12, YN12, ZN12, XN23, YN23, ZN23,
     :     XN31, YN31, ZN31, A123, A231, A312, B23, B31, B12
      REAL ZN8743, A81357, A31245, A2341, A134, A1234, A124, A12356,
     :     A21357, A18743, A18642, A1456, A18765, A21346, A213, A234,
     :     A28765, A3412, A324, A31256, YN1357, B13, XN8765, XN8743,
     :     A2456, A28743, A58743, A58642, A412, A38642, A3456, A38765,
     :     A42356, A4123, A413, A41256, A41357, A48765, A51234, A51346,
     :     A61234, A6123, A81256, A71256, A61245, A78642, B1357,
     :     B123, B2356, A5123, A61357, A68743, A71234, A81234, B1346,
     :     B1256, B1234, B1245, B412, B234, B341, B34, B24, B456, B8743,
     :     B8642, B8765, XN341, XN1357, ZN2356, XN24, XN1245, XN1234,
     :     XN123, XN13, XN1256, XN1346, XN234, XN2356, XN456, XN34,
     :     XN412, YN13, YN2356, YN234, YN24, YN456, YN34, YN341, YN412,
     :     YN8642, YN8743, YN8765, ZN1256, ZN123, ZN1234, ZN456, ZN8765,
     :     ZN1346, ZN13, ZN1245, ZN1357, ZN24, ZN234, ZN341, ZN34,
     :     ZN412, ZN8642, XN8642, YN1346, YN1234, YN123, YN1256, YN1245
C
      IF( NPNTS .EQ. 3 ) THEN
C
         X12 = X(2) - X(1)
         X23 = X(3) - X(2)
         X31 = X(1) - X(3)
         Y12 = Y(2) - Y(1)
         Y23 = Y(3) - Y(2)
         Y31 = Y(1) - Y(3)
         Z12 = Z(2) - Z(1)
         Z23 = Z(3) - Z(2)
         Z31 = Z(1) - Z(3)
C
         XP1 = X(1) - XP
         XP2 = X(2) - XP
         XP3 = X(3) - XP
         YP1 = Y(1) - YP
         YP2 = Y(2) - YP
         YP3 = Y(3) - YP
         ZP1 = Z(1) - ZP
         ZP2 = Z(2) - ZP
         ZP3 = Z(3) - ZP
C
         XNRM = -Y12*Z31 + Z12*Y31
         YNRM = -Z12*X31 + X12*Z31
         ZNRM = -X12*Y31 + Y12*X31
C
         XN12 = Y12*ZNRM - Z12*YNRM
         YN12 = Z12*XNRM - X12*ZNRM
         ZN12 = X12*YNRM - Y12*XNRM
C
         XN23 = Y23*ZNRM - Z23*YNRM
         YN23 = Z23*XNRM - X23*ZNRM
         ZN23 = X23*YNRM - Y23*XNRM
C
         XN31 = Y31*ZNRM - Z31*YNRM
         YN31 = Z31*XNRM - X31*ZNRM
         ZN31 = X31*YNRM - Y31*XNRM
C
         A123 = XN23*X12 + YN23*Y12 + ZN23*Z12
         A231 = XN31*X23 + YN31*Y23 + ZN31*Z23
         A312 = XN12*X31 + YN12*Y31 + ZN12*Z31
C
         B23  = XN23*XP2 + YN23*YP2 + ZN23*ZP2
         B31  = XN31*XP3 + YN31*YP3 + ZN31*ZP3
         B12  = XN12*XP1 + YN12*YP1 + ZN12*ZP1
C
         VALS(1) = B23/A123
         VALS(2) = B31/A231
         VALS(3) = B12/A312
C
      ELSE IF( NPNTS .EQ. 4 ) THEN
C
C - Remember that points 3 and 4 are wrong way round for anticlockwise quad.
C
         X12 = X(2) - X(1)
         X13 = X(3) - X(1)
         X14 = X(4) - X(1)
         X24 = X(4) - X(2)
         X34 = X(4) - X(3)
         Y12 = Y(2) - Y(1)
         Y13 = Y(3) - Y(1)
         Y14 = Y(4) - Y(1)
         Y24 = Y(4) - Y(2)
         Y34 = Y(4) - Y(3)
         Z12 = Z(2) - Z(1)
         Z13 = Z(3) - Z(1)
         Z14 = Z(4) - Z(1)
         Z24 = Z(4) - Z(2)
         Z34 = Z(4) - Z(3)
C
         X1P = XP - X(1)
         XP4 = X(4) - XP
         Y1P = YP - Y(1)
         YP4 = Y(4) - YP
         Z1P = ZP - Z(1)
         ZP4 = Z(4) - ZP
C
         XNRM = Y12*Z13 - Z12*Y13
         YNRM = Z12*X13 - X12*Z13
         ZNRM = X12*Y13 - Y12*X13
c         print*,'*** Normal: ',xnrm,ynrm,znrm
C
         XN12 = Y12*ZNRM - Z12*YNRM
         YN12 = Z12*XNRM - X12*ZNRM
         ZN12 = X12*YNRM - Y12*XNRM
C
         XN24 = Y24*ZNRM - Z24*YNRM
         YN24 = Z24*XNRM - X24*ZNRM
         ZN24 = X24*YNRM - Y24*XNRM
C
         XN34 = Y34*ZNRM - Z34*YNRM
         YN34 = Z34*XNRM - X34*ZNRM
         ZN34 = X34*YNRM - Y34*XNRM
C
         XN13 = Y13*ZNRM - Z13*YNRM
         YN13 = Z13*XNRM - X13*ZNRM
         ZN13 = X13*YNRM - Y13*XNRM
c         print*,'*** XN12,YN12,ZN12: ',XN12,YN12,ZN12
c         print*,'*** XN13,YN13,ZN13: ',XN13,YN13,ZN13
c         print*,'*** XN24,YN24,ZN24: ',XN24,YN24,ZN24
c         print*,'*** XN34,YN34,ZN34: ',XN34,YN34,ZN34
C
         A312 = XN12*X13 + YN12*Y13 + ZN12*Z13
         A412 = XN12*X14 + YN12*Y14 + ZN12*Z14
C
         B12  = XN12*X1P + YN12*Y1P + ZN12*Z1P
C
         A124 = XN24*X14 + YN24*Y14 + ZN24*Z14
         A324 = XN24*X34 + YN24*Y34 + ZN24*Z34
C
         B24  = XN24*XP4 + YN24*YP4 + ZN24*ZP4
C
         A134 = XN34*X14 + YN34*Y14 + ZN34*Z14
         A234 = XN34*X24 + YN34*Y24 + ZN34*Z24
C
         B34  = XN34*XP4 + YN34*YP4 + ZN34*ZP4
C
         A213 = XN13*X12 + YN13*Y12 + ZN13*Z12
         A413 = XN13*X14 + YN13*Y14 + ZN13*Z14
C
         B13  = XN13*X1P + YN13*Y1P + ZN13*Z1P
C
c         print*,'*** B12,B13: ',b12,b13
c         print*,'*** B24,B34: ',b24,b34
c         print*,'*** A124,A134: ',a124,a134
c         print*,'*** A213,A234: ',a213,a234
c         print*,'*** A312,A324: ',a312,a324
c         print*,'*** A412,A413: ',a412,a413
         VALS(1) = B24*B34/A124/A134
         VALS(2) = B13*B34/A213/A234
         VALS(3) = B12*B24/A312/A324
         VALS(4) = B12*B13/A412/A413
c         vv = vals(1)+vals(2)+vals(3)+vals(4)
c         if( abs(1.0-vv).gt.1e-4 ) then
c            print*,'VALS3D: ',(vals(ii),ii=1,4)
c            print*,'*** SUM IS NOT ONE ***'
c            stop
c         end if
C
      ELSE IF( NPNTS .EQ. 5 ) THEN
C
C - this is actually a tetrahedron which has only four points
C
         X12 = X(2) - X(1)
         X13 = X(3) - X(1)
         X14 = X(4) - X(1)
         X23 = X(3) - X(2)
         X24 = X(4) - X(2)
         Y12 = Y(2) - Y(1)
         Y13 = Y(3) - Y(1)
         Y14 = Y(4) - Y(1)
         Y23 = Y(3) - Y(2)
         Y24 = Y(4) - Y(2)
         Z12 = Z(2) - Z(1)
         Z13 = Z(3) - Z(1)
         Z14 = Z(4) - Z(1)
         Z23 = Z(3) - Z(2)
         Z24 = Z(4) - Z(2)
C
         X1P = XP - X(1)
         X2P = XP - X(2)
         Y1P = YP - Y(1)
         Y2P = YP - Y(2)
         Z1P = ZP - Z(1)
         Z2P = ZP - Z(2)
C
         XN123 = Y12*Z13 - Z12*Y13
         YN123 = Z12*X13 - X12*Z13
         ZN123 = X12*Y13 - Y12*X13
C
         XN412 = Y12*Z14 - Z12*Y14
         YN412 = Z12*X14 - X12*Z14
         ZN412 = X12*Y14 - Y12*X14
C
         XN234 = Y23*Z24 - Z23*Y24
         YN234 = Z23*X24 - X23*Z24
         ZN234 = X23*Y24 - Y23*X24
C
         XN341 = Y13*Z14 - Z13*Y14
         YN341 = Z13*X14 - X13*Z14
         ZN341 = X13*Y14 - Y13*X14
C
         A1234 = XN234*X12 + YN234*Y12 + ZN234*Z12
         A2341 = XN341*X12 + YN341*Y12 + ZN341*Z12
         A3412 = XN412*X13 + YN412*Y13 + ZN412*Z13
         A4123 = XN123*X14 + YN123*Y14 + ZN123*Z14
C
         B234  = XN234*X2P + YN234*Y2P + ZN234*Z2P
         B341  = XN341*X1P + YN341*Y1P + ZN341*Z1P
         B412  = XN412*X1P + YN412*Y1P + ZN412*Z1P
         B123  = XN123*X1P + YN123*Y1P + ZN123*Z1P
C
         VALS(1) = -B234/A1234
         VALS(2) = B341/A2341
         VALS(3) = B412/A3412
         VALS(4) = B123/A4123
C
      ELSE IF( NPNTS .EQ. 6 ) THEN
C
         X1P = XP - X(1)
         X6P = XP - X(6)
         Y1P = YP - Y(1)
         Y6P = YP - Y(6)
         Z1P = ZP - Z(1)
         Z6P = ZP - Z(6)
C
         X12 = X(2) - X(1)
         X13 = X(3) - X(1)
         X14 = X(4) - X(1)
         X15 = X(5) - X(1)
         X16 = X(6) - X(1)
         Y12 = Y(2) - Y(1)
         Y13 = Y(3) - Y(1)
         Y14 = Y(4) - Y(1)
         Y15 = Y(5) - Y(1)
         Y16 = Y(6) - Y(1)
         Z12 = Z(2) - Z(1)
         Z13 = Z(3) - Z(1)
         Z14 = Z(4) - Z(1)
         Z15 = Z(5) - Z(1)
         Z16 = Z(6) - Z(1)
C
         X65 = X(5) - X(6)
         X64 = X(4) - X(6)
         X63 = X(3) - X(6)
         X62 = X(2) - X(6)
         Y65 = Y(5) - Y(6)
         Y64 = Y(4) - Y(6)
         Y63 = Y(3) - Y(6)
         Y62 = Y(2) - Y(6)
         Z65 = Z(5) - Z(6)
         Z64 = Z(4) - Z(6)
         Z63 = Z(3) - Z(6)
         Z62 = Z(2) - Z(6)
C
         XN123  = Y12*Z13 - Z12*Y13
         YN123  = Z12*X13 - X12*Z13
         ZN123  = X12*Y13 - Y12*X13
C
         XN456  = Y64*Z65 - Z64*Y65
         YN456  = Z64*X65 - X64*Z65
         ZN456  = X64*Y65 - Y64*X65
C
         XN1245 = Y12*Z14 - Z12*Y14
         YN1245 = Z12*X14 - X12*Z14
         ZN1245 = X12*Y14 - Y12*X14
C
         XN2356 = Y63*Z65 - Z63*Y65
         YN2356 = Z63*X65 - X63*Z65
         ZN2356 = X63*Y65 - Y63*X65
C
         XN1346 = Y13*Z14 - Z13*Y14
         YN1346 = Z13*X14 - X13*Z14
         ZN1346 = X13*Y14 - Y13*X14
C
         A4123  = XN123*X14 + YN123*Y14 + ZN123*Z14
         A5123  = XN123*X15 + YN123*Y15 + ZN123*Z15
         A6123  = XN123*X16 + YN123*Y16 + ZN123*Z16
C
         B123   = XN123*X1P + YN123*Y1P + ZN123*Z1P
C
         A1456  =-XN456*X16 - YN456*Y16 - ZN456*Z16
         A2456  = XN456*X62 + YN456*Y62 + ZN456*Z62
         A3456  = XN456*X63 + YN456*Y63 + ZN456*Z63
C
         B456   = XN456*X6P + YN456*Y6P + ZN456*Z6P
C
         A31245 = XN1245*X13 + YN1245*Y13 + ZN1245*Z13
         A61245 = XN1245*X16 + YN1245*Y16 + ZN1245*Z16
C
         B1245  = XN1245*X1P + YN1245*Y1P + ZN1245*Z1P
C
         A12356 =-XN2356*X16 - YN2356*Y16 - ZN2356*Z16
         A42356 = XN2356*X64 + YN2356*Y64 + ZN2356*Z64
C
         B2356  = XN2356*X6P + YN2356*Y6P + ZN2356*Z6P
C
         A21346 = XN1346*X12 + YN1346*Y12 + ZN1346*Z12
         A51346 = XN1346*X65 + YN1346*Y65 + ZN1346*Z65
C
         B1346  = XN1346*X6P + YN1346*Y6P + ZN1346*Z6P
C
         VALS(1) = B456*B2356/A1456/A12356
         VALS(2) = B456*B1346/A2456/A21346
         VALS(3) = B456*B1245/A3456/A31245
         VALS(4) = B123*B2356/A4123/A42356
         VALS(5) = B1346*B123/A5123/A51346
         VALS(6) = B1245*B123/A6123/A61245
C
      ELSE IF( NPNTS .EQ. 8 ) THEN
C
         X1P = XP - X(1)
         X8P = XP - X(8)
         Y1P = YP - Y(1)
         Y8P = YP - Y(8)
         Z1P = ZP - Z(1)
         Z8P = ZP - Z(8)
C
         X12 = X(2) - X(1)
         X13 = X(3) - X(1)
         X15 = X(5) - X(1)
         X84 = X(4) - X(8)
         X86 = X(6) - X(8)
         X87 = X(7) - X(8)
C
         Y12 = Y(2) - Y(1)
         Y13 = Y(3) - Y(1)
         Y15 = Y(5) - Y(1)
         Y84 = Y(4) - Y(8)
         Y86 = Y(6) - Y(8)
         Y87 = Y(7) - Y(8)
C
         Z12 = Z(2) - Z(1)
         Z13 = Z(3) - Z(1)
         Z15 = Z(5) - Z(1)
         Z84 = Z(4) - Z(8)
         Z86 = Z(6) - Z(8)
         Z87 = Z(7) - Z(8)
C
         X14 = X(4) - X(1)
         X16 = X(6) - X(1)
         X17 = X(7) - X(1)
         X18 = X(8) - X(1)
         Y14 = Y(4) - Y(1)
         Y16 = Y(6) - Y(1)
         Y17 = Y(7) - Y(1)
         Y18 = Y(8) - Y(1)
         Z14 = Z(4) - Z(1)
         Z16 = Z(6) - Z(1)
         Z17 = Z(7) - Z(1)
         Z18 = Z(8) - Z(1)
C
         X82 = X(2) - X(8)
         X83 = X(3) - X(8)
         X85 = X(5) - X(8)
         Y82 = Y(2) - Y(8)
         Y83 = Y(3) - Y(8)
         Y85 = Y(5) - Y(8)
         Z82 = Z(2) - Z(8)
         Z83 = Z(3) - Z(8)
         Z85 = Z(5) - Z(8)
C
         XN1234 = Y12*Z13 - Y13*Z12
         YN1234 = Z12*X13 - Z13*X12
         ZN1234 = X12*Y13 - X13*Y12
C
         XN1357 = Y13*Z15 - Y15*Z13
         YN1357 = Z13*X15 - Z15*X13
         ZN1357 = X13*Y15 - X15*Y13
C
         XN1256 = Y15*Z12 - Y12*Z15
         YN1256 = Z15*X12 - Z12*X15
         ZN1256 = X15*Y12 - X12*Y15
C
         XN8765 = Y86*Z87 - Y87*Z86
         YN8765 = Z86*X87 - Z87*X86
         ZN8765 = X86*Y87 - X87*Y86
C
         XN8642 = Y84*Z86 - Y86*Z84
         YN8642 = Z84*X86 - Z86*X84
         ZN8642 = X84*Y86 - X86*Y84
C
         XN8743 = Y87*Z84 - Y84*Z87
         YN8743 = Z87*X84 - Z84*X87
         ZN8743 = X87*Y84 - X84*Y87
C
         A51234 = XN1234*X15 + YN1234*Y15 + ZN1234*Z15
         A61234 = XN1234*X16 + YN1234*Y16 + ZN1234*Z16
         A71234 = XN1234*X17 + YN1234*Y17 + ZN1234*Z17
         A81234 = XN1234*X18 + YN1234*Y18 + ZN1234*Z18
C
         B1234  = XN1234*X1P + YN1234*Y1P + ZN1234*Z1P
C
         A21357 = XN1357*X12 + YN1357*Y12 + ZN1357*Z12
         A41357 = XN1357*X14 + YN1357*Y14 + ZN1357*Z14
         A61357 = XN1357*X16 + YN1357*Y16 + ZN1357*Z16
         A81357 = XN1357*X18 + YN1357*Y18 + ZN1357*Z18
C
         B1357  = XN1357*X1P + YN1357*Y1P + ZN1357*Z1P
C
         A31256 = XN1256*X13 + YN1256*Y13 + ZN1256*Z13
         A41256 = XN1256*X14 + YN1256*Y14 + ZN1256*Z14
         A71256 = XN1256*X17 + YN1256*Y17 + ZN1256*Z17
         A81256 = XN1256*X18 + YN1256*Y18 + ZN1256*Z18
C
         B1256  = XN1256*X1P + YN1256*Y1P + ZN1256*Z1P
C
         A18765 =-XN8765*X18 - YN8765*Y18 - ZN8765*Z18
         A28765 = XN8765*X82 + YN8765*Y82 + ZN8765*Z82
         A38765 = XN8765*X83 + YN8765*Y83 + ZN8765*Z83
         A48765 = XN8765*X84 + YN8765*Y84 + ZN8765*Z84
C
         B8765  = XN8765*X8P + YN8765*Y8P + ZN8765*Z8P
C
         A18642 =-XN8642*X18 - YN8642*Y18 - ZN8642*Z18
         A38642 = XN8642*X83 + YN8642*Y83 + ZN8642*Z83
         A58642 = XN8642*X85 + YN8642*Y85 + ZN8642*Z85
         A78642 = XN8642*X87 + YN8642*Y87 + ZN8642*Z87
C
         B8642  = XN8642*X8P + YN8642*Y8P + ZN8642*Z8P
C
         A18743 =-XN8743*X18 - YN8743*Y18 - ZN8743*Z18
         A28743 = XN8743*X82 + YN8743*Y82 + ZN8743*Z82
         A58743 = XN8743*X85 + YN8743*Y85 + ZN8743*Z85
         A68743 = XN8743*X86 + YN8743*Y86 + ZN8743*Z86
C
         B8743  = XN8743*X8P + YN8743*Y8P + ZN8743*Z8P
C
         VALS(1) = B8765*B8642*B8743/A18765/A18642/A18743
         VALS(2) = B8765*B1357*B8743/A28765/A21357/A28743
         VALS(3) = B8765*B8642*B1256/A38765/A38642/A31256
         VALS(4) = B8765*B1357*B1256/A48765/A41357/A41256
         VALS(5) = B1234*B8642*B8743/A51234/A58642/A58743
         VALS(6) = B1234*B1357*B8743/A61234/A61357/A68743
         VALS(7) = B1234*B8642*B1256/A71234/A78642/A71256
         VALS(8) = B1234*B1357*B1256/A81234/A81357/A81256
C
      END IF
C
      RETURN
      END
C
