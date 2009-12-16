!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    C.Pain@Imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
module detnlxr_module
  implicit none
contains

  SUBROUTINE DETNLXR(ELE, X,Y,Z, XONDGL, TOTELE,NONODS,NLOC,NGI, &
       N,NLX,NLY,NLZ, WEIGHT, DETWEI,RA,VOLUME, D3,DCYL, &
       NX,NY,NZ) 
    IMPLICIT NONE
    REAL PIE
    PARAMETER(PIE=3.141592654)
    INTEGER ELE,TOTELE,NONODS,NLOC,NGI
    REAL X(NONODS),Y(NONODS),Z(NONODS)
    INTEGER XONDGL(TOTELE*NLOC)
    ! Volume shape functions...
    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)
    REAL DETWEI(NGI)
    REAL RA(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    REAL VOLUME
    LOGICAL D3,DCYL
    ! Local variables...
    REAL AGI,BGI,CGI, DGI,EGI,FGI, GGI,HGI,KGI
    REAL A11,A12,A13, A21,A22,A23, A31,A32,A33
    REAL DETJ,TWOPIE,RGI
    INTEGER GI,L,IGLX
    !
    VOLUME=0.
    !
    IF(D3) THEN
       do  GI=1,NGI! Was loop 331
          !
          AGI=0.
          BGI=0.
          CGI=0.
          !
          DGI=0.
          EGI=0.
          FGI=0.
          !
          GGI=0.
          HGI=0.
          KGI=0.
          !
          do  L=1,NLOC! Was loop 79
             IGLX=XONDGL((ELE-1)*NLOC+L)
             ! NB R0 does not appear here although the z-coord might be Z+R0. 
             AGI=AGI+NLX(L,GI)*X(IGLX) 
             BGI=BGI+NLX(L,GI)*Y(IGLX) 
             CGI=CGI+NLX(L,GI)*Z(IGLX) 
             !
             DGI=DGI+NLY(L,GI)*X(IGLX) 
             EGI=EGI+NLY(L,GI)*Y(IGLX) 
             FGI=FGI+NLY(L,GI)*Z(IGLX) 
             !
             GGI=GGI+NLZ(L,GI)*X(IGLX) 
             HGI=HGI+NLZ(L,GI)*Y(IGLX) 
             KGI=KGI+NLZ(L,GI)*Z(IGLX)
          end do ! Was loop 79
          !
          DETJ=AGI*(EGI*KGI-FGI*HGI)&
               -BGI*(DGI*KGI-FGI*GGI)&
               +CGI*(DGI*HGI-EGI*GGI)
          DETWEI(GI)=DETJ*WEIGHT(GI)
          RA(GI)=1.0
          VOLUME=VOLUME+DETWEI(GI)
          ! For coefficient in the inverse mat of the jacobian. 
          A11= (EGI*KGI-FGI*HGI) /DETJ
          A21=-(DGI*KGI-FGI*GGI) /DETJ
          A31= (DGI*HGI-EGI*GGI) /DETJ
          !
          A12=-(BGI*KGI-CGI*HGI) /DETJ
          A22= (AGI*KGI-CGI*GGI) /DETJ
          A32=-(AGI*HGI-BGI*GGI) /DETJ
          !
          A13= (BGI*FGI-CGI*EGI) /DETJ
          A23=-(AGI*FGI-CGI*DGI) /DETJ
          A33= (AGI*EGI-BGI*DGI) /DETJ
          do  L=1,NLOC! Was loop 373
             NX(L,GI)= A11*NLX(L,GI)+A12*NLY(L,GI)+A13*NLZ(L,GI)
             NY(L,GI)= A21*NLX(L,GI)+A22*NLY(L,GI)+A23*NLZ(L,GI)
             NZ(L,GI)= A31*NLX(L,GI)+A32*NLY(L,GI)+A33*NLZ(L,GI)
          end do ! Was loop 373
          !
       end do ! Was loop 331
       ! IF(D3) THEN...
    ELSE
       TWOPIE=1.0 
       IF(DCYL) TWOPIE=2.*PIE
       do  GI=1,NGI! Was loop 1331
          !
          RGI=0.
          !
          AGI=0.
          BGI=0.
          CGI=0.
          DGI=0.
          !
          do  L=1,NLOC! Was loop 179
             IGLX=XONDGL((ELE-1)*NLOC+L)

             AGI=AGI + NLX(L,GI)*X(IGLX) 
             BGI=BGI + NLX(L,GI)*Y(IGLX) 
             CGI=CGI + NLY(L,GI)*X(IGLX) 
             DGI=DGI + NLY(L,GI)*Y(IGLX) 
             !
             RGI=RGI+N(L,GI)*Y(IGLX)
          end do ! Was loop 179
          !
          IF(.NOT.DCYL) RGI=1.0
          !
          DETJ= AGI*DGI-BGI*CGI 
          RA(GI)=RGI
          DETWEI(GI)=TWOPIE*RGI*DETJ*WEIGHT(GI)
          VOLUME=VOLUME+DETWEI(GI)
          !
          do L=1,NLOC
             NX(L,GI)=(DGI*NLX(L,GI)-BGI*NLY(L,GI))/DETJ
             NY(L,GI)=(-CGI*NLX(L,GI)+AGI*NLY(L,GI))/DETJ
             NZ(L,GI)=0.0
          END DO
          !
       end do ! Was loop 1331
       ! ENDOF IF(D3) THEN ELSE...
    ENDIF
    !
    RETURN
  END SUBROUTINE DETNLXR
end module detnlxr_module
