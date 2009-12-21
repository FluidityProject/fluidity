!     Copyright (C) 2006 Imperial College London and others.
!     
!     Please see the AUTHORS file in the main source directory for a full list
!     of copyright holders.
!     
!     Prof. C Pain
!     Applied Modelling and Computation Group
!     Department of Earth Science and Engineering
!     Imperial College London
!     
!     C.Pain@Imperial.ac.uk
!     
!     This library is free software; you can redistribute it and/or
!     modify it under the terms of the GNU Lesser General Public
!     License as published by the Free Software Foundation,
!     version 2.1 of the License.
!     
!     This library is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!     Lesser General Public License for more details.
!     
!     You should have received a copy of the GNU Lesser General Public
!     License along with this library; if not, write to the Free Software
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!     USA

#include "fdebug.h"

! This module contains a bunch of routines which did not find a home
! elsewhere in a proper module - perhaps because of insane
! cross-dependencies. All routines placed in this module should only
! contain dependencies satisified by other routines within this module.

module AllSorts

  use fldebug
  use parallel_tools
  
  implicit none
  
  private
  
  public :: anofac, detnlx, detnnn, exptol, locfac, sdetnn, tolfun, xprod, &
    & r2norm, ibuble, pminmx, comper, alomem, chkbc, icopy, raddin
  
contains

  SUBROUTINE ANOFAC(I1,I2,I3, X,Y,Z,D0, &
       SAX,SAY,SAZ, ZERO,XNONOD,FX,FY,FZ)
    ! This sub finds the normal given to vectors on the 
    ! surface and accumpulates the results in SAX,SAY,SAZ
    ! Element centre to face centre vector(approx normal) is FX,FY,FZ

    INTEGER I1,I2,I3,XNONOD
    LOGICAL ZERO
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD),D0
    REAL SAX,SAY,SAZ, FX,FY,FZ
    ! Local variables...
    REAL AX,AY,AZ, BX,BY,BZ, CX,CY,CZ
    REAL DD0
    
    ! Ocean scaling...
    DD0=1.
    IF(D0.NE.0.) DD0=D0 
    
    BX=X(I2)-X(I1)
    BY=Y(I2)-Y(I1)
    BZ=(Z(I2)-Z(I1))*DD0
    
    CX=X(I3)-X(I1)
    CY=Y(I3)-Y(I1)
    CZ=(Z(I3)-Z(I1))*DD0
    ! Perform x-product. 
    CALL XPROD(AX,AY,AZ, BX,BY,BZ, CX,CY,CZ)
    IF(FX*AX+FY*AY+FZ*AZ.LT.0) THEN
       ! Reverse direction of vector A - MAKE IT NORMAL to face of element.
       AX=-AX
       AY=-AY
       AZ=-AZ
    ENDIF
    
    IF(ZERO) THEN
       SAX=0.
       SAY=0.
       SAZ=0.
    ENDIF
    SAX=SAX+AX
    SAY=SAY+AY
    SAZ=SAZ+AZ
    RETURN
  END subroutine ANOFAC
  
  SUBROUTINE DETNNN(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI, &
       N,NLX,NLY,NLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL) 
    REAL PIE
    PARAMETER(PIE=3.141592654)
    INTEGER, intent(in):: ELE,TOTELE,XNONOD,NLOC,NGI
    REAL, intent(in):: X(XNONOD),Y(XNONOD),Z(XNONOD)
    INTEGER, intent(in):: XONDGL(TOTELE*NLOC)
    ! Volume shape functions...
    REAL, intent(in):: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL, intent(in):: WEIGHT(NGI)
    REAL, intent(out):: DETWEI(NGI)
    REAL, intent(out):: VOLUME
    LOGICAL, intent(in):: D3,DCYL
    ! Local variables...
    REAL AGI,BGI,CGI, DGI,EGI,FGI, GGI,HGI,KGI
    REAL DETJ,TWOPIE,RGI
    INTEGER GI,L,IGLX
    
    VOLUME=0. 
    
    IF(D3) THEN
       DO GI=1,NGI
          
          AGI=0.
          BGI=0.
          CGI=0.
          
          DGI=0.
          EGI=0.
          FGI=0.
          
          GGI=0.
          HGI=0.
          KGI=0.
          
          DO L=1,NLOC
             IGLX=XONDGL((ELE-1)*NLOC+L)
             ! NB R0 does not appear here although the z-coord might be Z+R0. 
             AGI=AGI+NLX(L,GI)*X(IGLX) 
             BGI=BGI+NLX(L,GI)*Y(IGLX) 
             CGI=CGI+NLX(L,GI)*Z(IGLX) 
             
             DGI=DGI+NLY(L,GI)*X(IGLX) 
             EGI=EGI+NLY(L,GI)*Y(IGLX) 
             FGI=FGI+NLY(L,GI)*Z(IGLX) 
             
             GGI=GGI+NLZ(L,GI)*X(IGLX) 
             HGI=HGI+NLZ(L,GI)*Y(IGLX) 
             KGI=KGI+NLZ(L,GI)*Z(IGLX)
          end do
          
          DETJ=AGI*(EGI*KGI-FGI*HGI)  &
               -BGI*(DGI*KGI-FGI*GGI) &
               +CGI*(DGI*HGI-EGI*GGI)
          DETWEI(GI)=abs(DETJ)*WEIGHT(GI)
          VOLUME=VOLUME+DETWEI(GI)

       end do
    ELSE
       TWOPIE=1.0 
       IF(DCYL) TWOPIE=2.*PIE
       DO GI=1,NGI
          
          RGI=0.
          
          AGI=0.
          BGI=0.
          CGI=0.
          DGI=0.
          
          DO L=1,NLOC
             IGLX=XONDGL((ELE-1)*NLOC+L)
             
             AGI=AGI + NLX(L,GI)*X(IGLX) 
             BGI=BGI + NLX(L,GI)*Y(IGLX) 
             CGI=CGI + NLY(L,GI)*X(IGLX) 
             DGI=DGI + NLY(L,GI)*Y(IGLX) 
             
             RGI=RGI+N(L,GI)*Y(IGLX)
          END DO
          
          IF(.NOT.DCYL) RGI=1.0
          
          DETJ= AGI*DGI-BGI*CGI 
          DETWEI(GI)=TWOPIE*RGI*DETJ*WEIGHT(GI)
          VOLUME=VOLUME+DETWEI(GI)
       end do
    ENDIF
    
    RETURN
  END subroutine DETNNN

  SUBROUTINE DETNLX(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI,  &
       N,NLX,NLY,NLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL, &
       NX,NY,NZ) 
    REAL PIE
    PARAMETER(PIE=3.141592654)
    INTEGER ELE,TOTELE,XNONOD,NLOC,NGI
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    INTEGER XONDGL(TOTELE*NLOC)
    ! Volume shape functions...
    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)
    REAL DETWEI(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    REAL VOLUME
    LOGICAL D3,DCYL
    ! Local variables...
    REAL AGI,BGI,CGI, DGI,EGI,FGI, GGI,HGI,KGI
    REAL A11,A12,A13, A21,A22,A23, A31,A32,A33
    REAL DETJ,TWOPIE,RGI
    INTEGER GI,L,IGLX
    
    VOLUME=0.
    
    IF(D3) THEN
       DO GI=1,NGI
          
          AGI=0.
          BGI=0.
          CGI=0.
          
          DGI=0.
          EGI=0.
          FGI=0.
          
          GGI=0.
          HGI=0.
          KGI=0.
          
          DO L=1,NLOC
             IGLX=XONDGL((ELE-1)*NLOC+L)
             ! NB R0 does not appear here although the z-coord might be Z+R0. 
             AGI=AGI+NLX(L,GI)*X(IGLX) 
             BGI=BGI+NLX(L,GI)*Y(IGLX) 
             CGI=CGI+NLX(L,GI)*Z(IGLX) 
             
             DGI=DGI+NLY(L,GI)*X(IGLX) 
             EGI=EGI+NLY(L,GI)*Y(IGLX) 
             FGI=FGI+NLY(L,GI)*Z(IGLX) 
             
             GGI=GGI+NLZ(L,GI)*X(IGLX) 
             HGI=HGI+NLZ(L,GI)*Y(IGLX) 
             KGI=KGI+NLZ(L,GI)*Z(IGLX)
          end do
          
          DETJ=AGI*(EGI*KGI-FGI*HGI) &
               -BGI*(DGI*KGI-FGI*GGI) &
               +CGI*(DGI*HGI-EGI*GGI)
          DETWEI(GI)=abs(DETJ)*WEIGHT(GI)
          VOLUME=VOLUME+DETWEI(GI)
          ! For coefficient in the inverse mat of the jacobian. 
          A11= (EGI*KGI-FGI*HGI) /DETJ
          A21=-(DGI*KGI-FGI*GGI) /DETJ
          A31= (DGI*HGI-EGI*GGI) /DETJ
          
          A12=-(BGI*KGI-CGI*HGI) /DETJ
          A22= (AGI*KGI-CGI*GGI) /DETJ
          A32=-(AGI*HGI-BGI*GGI) /DETJ
          
          A13= (BGI*FGI-CGI*EGI) /DETJ
          A23=-(AGI*FGI-CGI*DGI) /DETJ
          A33= (AGI*EGI-BGI*DGI) /DETJ
          DO L=1,NLOC
             NX(L,GI)= A11*NLX(L,GI)+A12*NLY(L,GI)+A13*NLZ(L,GI)
             NY(L,GI)= A21*NLX(L,GI)+A22*NLY(L,GI)+A23*NLZ(L,GI)
             NZ(L,GI)= A31*NLX(L,GI)+A32*NLY(L,GI)+A33*NLZ(L,GI)
          end do
          
       end do
    ELSE
       TWOPIE=1.0 
       IF(DCYL) TWOPIE=2.*PIE
       DO GI=1,NGI
          
          RGI=0.
          
          AGI=0.
          BGI=0.
          CGI=0.
          DGI=0.
          
          DO L=1,NLOC
             IGLX=XONDGL((ELE-1)*NLOC+L)
             
             AGI=AGI + NLX(L,GI)*X(IGLX) 
             BGI=BGI + NLX(L,GI)*Y(IGLX) 
             CGI=CGI + NLY(L,GI)*X(IGLX) 
             DGI=DGI + NLY(L,GI)*Y(IGLX) 
             
             RGI=RGI+N(L,GI)*Y(IGLX)
          end do
          
          IF(.NOT.DCYL) RGI=1.0
          
          DETJ= AGI*DGI-BGI*CGI 
          DETWEI(GI)=TWOPIE*RGI*DETJ*WEIGHT(GI)
          VOLUME=VOLUME+DETWEI(GI)
          
          DO L=1,NLOC
             NX(L,GI)=(DGI*NLX(L,GI)-BGI*NLY(L,GI))/DETJ
             NY(L,GI)=(-CGI*NLX(L,GI)+AGI*NLY(L,GI))/DETJ
             NZ(L,GI)=0.0
          END DO
          
       end do
       
    ENDIF
    
    RETURN
  END subroutine DETNLX

  function exptol(x)
    !!< Return the exponential of X, or zero is X is less than - tolex
    !!< (tolext = 40.0, hardcoded).

      real, intent(in) :: x

      real :: exptol

      real :: tolext
      parameter(tolext = 40.0)

      if(x < -tolext) then
        exptol = 0.0
      else
        exptol = exp(x)
      endif

  end function exptol

  SUBROUTINE LOCFAC(IFACE,LOCNOD,SNLOC)
    ! Get local node ers
    ! Put local no.s on faces by looking at faces from outside element. 
    INTEGER IFACE,SNLOC,LOCNOD(SNLOC)
    
    IF(SNLOC.EQ.2) THEN
       ! 2-D
       IF(IFACE.EQ.1) THEN
          ! For face 1 
          LOCNOD(1)=3
          LOCNOD(2)=1
       ENDIF
       IF(IFACE.EQ.2) THEN
          ! For face 2
          LOCNOD(1)=2
          LOCNOD(2)=4
       ENDIF
       IF(IFACE.EQ.3) THEN
          ! For faces 3 
          LOCNOD(1)=1
          LOCNOD(2)=2
       ENDIF
       IF(IFACE.EQ.4) THEN
          ! For faces 4 
          LOCNOD(1)=4 
          LOCNOD(2)=3
       ENDIF
    ELSE
       ! 3-D
       IF(IFACE.EQ.1) THEN
          ! For face 1 
          LOCNOD(1)=3
          LOCNOD(2)=1
          LOCNOD(3)=7
          LOCNOD(4)=5 
       ENDIF
       IF(IFACE.EQ.2) THEN
          ! For face 2
          LOCNOD(1)=2
          LOCNOD(2)=4
          LOCNOD(3)=6
          LOCNOD(4)=8
       ENDIF
       
       IF(IFACE.EQ.3) THEN
          ! For faces 3 
          LOCNOD(1)=1
          LOCNOD(2)=2
          LOCNOD(3)=5
          LOCNOD(4)=6
       ENDIF
       IF(IFACE.EQ.4) THEN
          ! For faces 4 
          LOCNOD(1)=4 
          LOCNOD(2)=3
          LOCNOD(3)=8
          LOCNOD(4)=7 
       ENDIF
       
       IF(IFACE.EQ.5) THEN
          ! For face 5 
          LOCNOD(1)=2
          LOCNOD(2)=1
          LOCNOD(3)=4
          LOCNOD(4)=3
       ENDIF
       IF(IFACE.EQ.6) THEN
          ! For face 6 
          LOCNOD(1)=5
          LOCNOD(2)=6
          LOCNOD(3)=7
          LOCNOD(4)=8
       ENDIF
    ENDIF
    RETURN
  END subroutine LOCFAC

  REAL FUNCTION R2NORM(VEC,FREDOP,PARA)

    real, intent(in)::VEC(:)
    integer, intent(in)::FREDOP, PARA

    real::rsum, psum

#ifdef HAVE_MPI
    include 'mpif.h'
    integer::mierr
#endif
    integer::i

    rsum=0.0
    do I=1,FREDOP
       RSUM=RSUM+VEC(I)**2
    END DO

#ifdef HAVE_MPI
    IF(PARA.EQ.1) then
       call MPI_Allreduce(rsum, psum, 1, getpreal(), MPI_SUM, MPI_COMM_WORLD, mierr)
       assert(mierr == MPI_SUCCESS)
       rsum = psum
    end if
#endif

    R2NORM=SQRT(RSUM)
  END FUNCTION R2NORM

  SUBROUTINE SDETNN(LOCNOD, &
       XNONOD,SNLOC,SNGI, &
       X,Y,Z, &
       SN,SNLX,SNLY, SWEIGH, DETWEI,SAREA, D3,DCYL) 
    INTEGER, intent(in)::XNONOD,SNLOC,SNGI
    INTEGER, intent(in)::LOCNOD(SNLOC)
    REAL, intent(in)::X(XNONOD),Y(XNONOD),Z(XNONOD)
    REAL, intent(in)::SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI) 
    REAL, intent(in)::SWEIGH(SNGI) 
    LOGICAL, intent(in)::D3,DCYL 
    
    REAL, intent(out)::DETWEI(SNGI), SAREA
    
    ! Local variables...
    real, parameter::PIE=3.141592654
    INTEGER GI,SL,L,IGLX
    REAL DXDLX,DXDLY,DYDLX,DYDLY,DZDLX,DZDLY
    REAL A,B,C,DETJ,RGI,TWOPIE
    
    SAREA=0.
    
    IF(D3) THEN
       DO GI=1,SNGI
          
          DXDLX=0.
          DXDLY=0.
          DYDLX=0.
          DYDLY=0.
          DZDLX=0.
          DZDLY=0.
          
          DO SL=1,SNLOC
             L=LOCNOD(SL) 
             IGLX=L
             DXDLX=DXDLX + SNLX(SL,GI)*X(IGLX)
             DXDLY=DXDLY + SNLY(SL,GI)*X(IGLX) 
             DYDLX=DYDLX + SNLX(SL,GI)*Y(IGLX) 
             DYDLY=DYDLY + SNLY(SL,GI)*Y(IGLX) 
             DZDLX=DZDLX + SNLX(SL,GI)*Z(IGLX) 
             DZDLY=DZDLY + SNLY(SL,GI)*Z(IGLX) 
          end do
          A = DYDLX*DZDLY - DYDLY*DZDLX
          B = DXDLX*DZDLY - DXDLY*DZDLX
          C = DXDLX*DYDLY - DXDLY*DYDLX
          
          DETJ=SQRT( A**2 + B**2 + C**2)
          DETWEI(GI)=abs(DETJ)*SWEIGH(GI)
          SAREA=SAREA+DETWEI(GI)

       end do
    ELSE
       TWOPIE=1.0 
       IF(DCYL) TWOPIE=2.*PIE
       
       DO GI=1,SNGI
          RGI=0.
          DXDLX=0.
          DXDLY=0.
          DYDLX=0.
          DYDLY=0.
          DZDLX=0.
          DZDLY=1.
          DO SL=1,SNLOC
             L=LOCNOD(SL) 
             IGLX=L
             DXDLX=DXDLX + SNLX(SL,GI)*X(IGLX) 
             DYDLX=DYDLX + SNLX(SL,GI)*Y(IGLX) 
             RGI=RGI+SN(SL,GI)*Y(IGLX)
          end do
          IF(.NOT.DCYL) RGI=1.0 
          DETJ=SQRT( DXDLX**2 + DYDLX**2 )
          DETWEI(GI)=TWOPIE*RGI*DETJ*SWEIGH(GI)
          SAREA=SAREA+DETWEI(GI) 
       end do
    ENDIF
    RETURN
  END subroutine SDETNN

  function tolfun(value)
!     ------------------------
!     
!     - this function is a tolerance function for a 
!     - value which is used as a denominator.
!     - If the absolute value of VALUE less than 1E-10 
!     - then it returns SIGN(A,B) i.e. the absolute value
!     - of A times the sign of B where A is TOLERANCE
!     - and B is VALUE.
!     
!     -------------------------------
!     - date last modified : 17/03/2002
!     -------------------------------
!     
!     
      REAL TOLFUN, VALUE
!     
!     - local variables
!     
      REAL TOLERANCE
      PARAMETER( TOLERANCE = 1E-10 )
!     
      IF( ABS(VALUE) .LT. TOLERANCE ) THEN
!     
         TOLFUN = SIGN(TOLERANCE,VALUE)
!     
      ELSE
!     
         TOLFUN = VALUE
!     
      END IF
!     
      RETURN

  end function tolfun
  
  SUBROUTINE XPROD(AX,AY,AZ, BX,BY,BZ, CX,CY,CZ)
    REAL AX,AY,AZ, BX,BY,BZ, CX,CY,CZ
    ! Perform x-product. a=b x c 
    AX=  BY*CZ-BZ*CY
    AY=-(BX*CZ-BZ*CX)
    AZ=  BX*CY-BY*CX
    RETURN
  END subroutine XPROD
  
  SUBROUTINE IBUBLE(LIST,NLIST)

    INTEGER NLIST,LIST(NLIST)
    INTEGER I,J,II
    do I=1,NLIST
       do J=2,NLIST
          IF(LIST(J-1).GT.LIST(J)) THEN
             !     SWOP
             II=LIST(J-1)
             LIST(J-1)=LIST(J)
             LIST(J)=II
          ENDIF
       END DO
    END DO
  END SUBROUTINE IBUBLE
  
  SUBROUTINE RIBUBLE(ILIST,RLIST,NLIST)

    INTEGER NLIST
    INTEGER ILIST(NLIST)
    REAL RLIST(NLIST)
    INTEGER I,J,II
    REAL RII
    do I=1,NLIST
       do J=2,NLIST
          IF(RLIST(J-1).GT.RLIST(J)) THEN
             !     SWOP
             RII=RLIST(J-1)
             RLIST(J-1)=RLIST(J)
             RLIST(J)=RII
             
             II=ILIST(J-1)
             ILIST(J-1)=ILIST(J)
             ILIST(J)=II
          ENDIF
       END DO
    END DO
  END SUBROUTINE RIBUBLE
  
  SUBROUTINE COMPER(PARA,TEMPT,T,NONODS,CHANGE)

    INTEGER PARA, NONODS
    REAL TEMPT(NONODS),T(NONODS),CHANGE
    REAL RR
    INTEGER I
    RR=0.
    do  I=1,NONODS! Was loop 10
       RR=MAX(RR,ABS(T(I)-TEMPT(I)) )
    end do ! Was loop 10
    IF(PARA.EQ.1) CALL ALLMAX(RR)
    CHANGE=MAX(RR,CHANGE)
  END SUBROUTINE COMPER
  
  SUBROUTINE PMINMX(VEC,NONODS,STR)

    REAL, parameter::LARGE=1.E+20
    INTEGER, intent(in)::NONODS
    REAL, intent(in)::VEC(NONODS)
    CHARACTER*(*), intent(in)::STR
    
    INTEGER KKK
    real rmax, rmin
    
    RMAX = maxval(vec)
    RMIN = minval(vec)
    
    KKK = INDEX(STR,'  ') - 1
    ewrite(3,*) STR(1:KKK),' MIN&MAX:',RMIN,RMAX
    RETURN
  END SUBROUTINE PMINMX

  SUBROUTINE ALOMEM(VEC, RPT, NONODS, NRMEM)

    integer, intent(out)::vec
    integer, intent(inout)::rpt
    integer, intent(in)::nonods, nrmem
    
    VEC=RPT
    RPT=RPT+NONODS
    IF( NONODS .LT. 0 ) THEN
       FLAbort('CANNOT ASK FOR NEGATIVE AMOUNT!')
    ELSE IF( VEC .GT. NRMEM .OR. VEC .LE. 0 ) THEN
       FLAbort('CURRENT POSITION OUT OF RANGE!')
    ELSE IF( RPT .GT. NRMEM+1) THEN
       FLAbort('NOT ENOUGH SPACE LEFT')
    ENDIF
    
    RETURN
  END SUBROUTINE ALOMEM
  
  subroutine chkbc(NOBCU,BCU2,NONODS)
    INTEGER NONODS,NOBCU, BCU2(NOBCU)
    INTEGER I,NOD
      
    IF(NOBCU>NONODS) THEN
       ewrite(0,*)'4 THE B.C s are out of bounds'
       ewrite(0,*)'duplicates? NOBCU,NONODS:',NOBCU,NONODS
    ENDIF
    do I=1,NOBCU
      NOD=BCU2(I)
      IF((NOD.LT.1).OR.(NOD.GT.NONODS)) THEN
        ewrite(-1,*)'5 THE B.C s are out of bounds NOD=',NOD
        FLAbort("STOP 5")
      ENDIF
    END DO

  end subroutine chkbc

  SUBROUTINE ICOPY(VEC1,VEC2,NONODS)
    INTEGER NONODS, VEC1(NONODS), VEC2(NONODS)
    INTEGER I
    do I=1,NONODS
       VEC1(I)=VEC2(I)
    END DO
  END subroutine icopy


  SUBROUTINE RADDIN(VEC1,VEC2,NONODS)
    INTEGER NONODS,I
!     This sub adds vector VEC2 to VEC1. 
    REAL VEC1(NONODS),VEC2(NONODS)
    do I=1,NONODS
       VEC1(I)=VEC1(I)+VEC2(I)
    END DO

  END SUBROUTINE RADDIN

end module AllSorts
