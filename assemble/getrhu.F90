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
#include "fdebug.h"

module getrhu_module

  use solvers
  use parallel_tools
  use flcomms_module
  use FLDebug
  use AllSorts
  use mulinv_module

  implicit none

  private

  public :: getrhu

contains

  SUBROUTINE GETRHU(UZAWA,RHS,C1T,C2T,C3T,&
       FINDCT,COLCT,NCT,&
       U,V,W,VECX,VECY,VECZ,  &
       VEC,VECSEC,DT,&
       FINDRM,COLM,CENTRM,NCOLM, NBIGM,&
       NONODS,FREDOP,ILENG,D2HALF,&
       PARA,halo_tag,NNODP,NNODPP,&
       ML, &
       !     Th following is for rotations only ROTAT=.TRUE.
       ROTAT,NNODRO,NRTDR,   NODROT,&
       NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
       DM1,DM2,DM3, RTDR)
    !     This only works for UZAWA=1 and COUPLE=1 *********** 
    !     This sub is called if ROTAT=.true.
    !     If D2HALF then ILENG=3*NONODS 
    !     we can treat as a 3-D flow. 
    !     If couple we will assume the boundary conditions are 
    !     already in GLM WHICH IS THE GLOBAL MATRIX. 
    INTEGER NCOLM, NBIGM, IMM, NNODP
    INTEGER NNODPP, NRTDR
    REAL RELAX
    INTEGER MATSTR
    PARAMETER(MATSTR=0)
    INTEGER FREDOP,NONODS,ILENG,D2HALF
    !     CONMAS=.0 mass is lumped, CONMAS=:1. mass is consistent. 
    !     NB UZAWA use Uzawas method of pressure determination.  
    !     D3 then solution is truely 3-D. 
    !     COUPLE then the momentum equations are coupled. 
    !     D2HALF then we have a two and a half D simulation in (r,z,theta)coords.
    INTEGER UZAWA
    INTEGER NCT,NOITS
    INTEGER D3
    REAL RHS(FREDOP),DT 
    REAL C1T(NCT),C2T(NCT),C3T(NCT)
    !     NCT=4*FREDOP for simple element. 
    REAL U(NONODS),V(NONODS),W(NONODS)
    REAL VECX(NONODS),VECY(NONODS)
    REAL VECZ(NONODS)
    REAL VEC(ILENG),VECSEC(ILENG)
    !     R is only used if CONMAS is on.
    INTEGER I,J,FINDCT(FREDOP+1),COLCT(NCT),POS
    INTEGER FINDRM(NONODS+1),COLM(NCOLM),CENTRM(NONODS)
    REAL ML(NONODS)
    !     For parallel processing ...
    INTEGER  PARA,halo_tag
    !     NB The normal is pointing out of the domain. 
    !     The rotation matrix in 3-D is R=  
    !     T1X   T1Y   T1Z
    !     T2X   T2Y   T2Z
    !     NX    NY    NZ
    !     The rotation matrix in 2-D is R=  
    !     T1X   T1Y   
    !     NX    NY    
    !     
    INTEGER NNODRO
    REAL NX(NNODRO),NY(NNODRO),NZ(NNODRO)
    REAL T1X(NNODRO),T1Y(NNODRO),T1Z(NNODRO)
    REAL T2X(NNODRO),T2Y(NNODRO),T2Z(NNODRO)
    INTEGER NODROT(NNODRO)
    REAL DM1(NONODS),DM2(NONODS),DM3(NONODS)
    REAL RTDR(NRTDR)
    LOGICAL ROTAT
    INTEGER IDIM, KITS
    REAL normC1T, normC2T, normC3T
    !     RTDR contains RT D R   ^T where D contains DMI's and R is the rotation matrix. 
    !     If D3 or D2HALF then ILENG=3*NONODS else =2*NONODS.
    !     This subroutine actually forms the Right Hand Side of Poissons
    !     pressure equation
    !     

    ewrite(1,*)'INSIDE GETRHU', PARA
    
    if (IsParallel()) then
      CALL HALGET(U,NONODS,NONODS,NNODP,halo_tag)
      CALL HALGET(V,NONODS,NONODS,NNODP,halo_tag)
      IF((ILENG.EQ.3*NONODS).AND.(.NOT.(D2HALF.EQ.1))) THEN
         CALL HALGET(W,NONODS,NONODS,NNODP,halo_tag)
      ENDIF
    end if

    IF(ROTAT.AND.(UZAWA.EQ.0)) THEN
       ewrite(3,*)  "CALL MULINV()"
       !     IF(.FALSE.) THEN
       IDIM=ILENG/NONODS
       D3=0
       IF(IDIM.EQ.3) D3=1
       !     Form vecxec=(R M_L R^T+D)^{-1} vecx
       !     NB R R^T=I so we do not need this
       !     NB (R M_L R^T +D)=(M_L + D)
       !     Because M_L =  M_l  0
       !     0    M_l
       !     
       !     ewrite(3,*) 'VECZ:',VECZ 
       !     ewrite(3,*) 'ml:',ml
       CALL MULINV(VECSEC,VECSEC(NONODS+1),&
            VECSEC((IDIM-1)*NONODS+1),ML,ML,ML,&
            VECX,VECY,VECZ,NONODS,D3,&
                                !     Th following is for rotations only ROTAT=.TRUE.
                                !     &         ROTAT,NNODRO,  NODROT,
            ROTAT, 0,  NODROT,&
            NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
            DM1,DM2,DM3 )
    ELSE
       FLAbort("This can never have worked")
    ENDIF

    ewrite(3,*)'INSIDE GETRHU HERE2', NONODS,FREDOP, ILENG
    !      ewrite(3,*)  "GUTS = ", R2NORM(VECSEC, NNODP, 1), R2NORM(VECSEC(NONODS+1), NNODP, 1), R2NORM(VECSEC(2*NONODS+1), NNODP, 1) 
    ewrite(3,*)&
         "GUTS = ", R2NORM(VECSEC, NNODP, 0), &
         R2NORM(VECSEC(NONODS+1:NONODS+nnodp), NNODP, 0), &
         R2NORM(VECSEC(2*NONODS+1:2*NONODS+nnodp), NNODP, 0) 
    ewrite(3,*)&
         "GUTS1 = ", R2NORM(V, NNODP, 0), R2NORM(U, NNODP, 0), R2NORM(W, NNODP, 0) 
    !      ewrite(3,*)  "GUTS1 = ", R2NORM(V, NNODP, 1), R2NORM(U, NNODP, 1), R2NORM(W, NNODP, 1) 
    IF((ILENG.EQ.3*NONODS).AND.(.NOT.(D2HALF.EQ.1))) THEN

       do I=1,FREDOP
          RHS(I)=0.0
          do POS=FINDCT(I),FINDCT(I+1)-1
             J=COLCT(POS)
             RHS(I) = RHS(I)&
                  - C1T(POS)*(VECSEC(J)         + U(J)/DT)&
                  - C2T(POS)*(VECSEC(J+NONODS)  + V(J)/DT)&
                  - C3T(POS)*(VECSEC(J+2*NONODS)+ W(J)/DT)
             !     We can subtract cty from RHS.
          END DO
       END DO

       normC1T = 0.0
       normC2T = 0.0
       normC3T = 0.0
       do I=1,NNODPP
          do POS=FINDCT(I),FINDCT(I+1)-1
             normC1T = normC1T + C1T(POS)*C1T(POS)
             normC2T = normC2T + C2T(POS)*C2T(POS)
             normC3T = normC3T + C3T(POS)*C3T(POS)
          END DO
       END DO
       CALL ALLSUM(normC1T)
       CALL ALLSUM(normC2T)
       CALL ALLSUM(normC3T)
       ewrite(3,*)  "Gore = ", normC1T,normC2T,normC3T
    ELSE
       do I=1,FREDOP
          RHS(I)=0.0
          do POS=FINDCT(I),FINDCT(I+1)-1
             J=COLCT(POS)
             RHS(I) = RHS(I)&
                  - C1T(POS)*(VECSEC(J)       + U(J)/DT)&
                  - C2T(POS)*(VECSEC(J+NONODS)+ V(J)/DT)
             !     We can subtract cty from RHS.
          END DO
       END DO
    ENDIF
  END SUBROUTINE GETRHU
end module getrhu_module
