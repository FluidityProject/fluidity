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
!
module getrhs_module
use FLDebug
use AllSorts
use flcomms_module
use parallel_tools
implicit none

contains

  SUBROUTINE GETRHS(RHS,C1T,C2T,C3T,&
       FINDCT,COLCT,NCT,&
       U,V,W,DM1,DM2,DM3,VECX,VECY,VECZ,  &
       SAVM,CONTRI,&
       DT,&
       ML,&
       FINDRM,COLM,CENTRM,NCOLM, &
       NONODS,FREDOP,&
       CONMAS,D3, &
       PARA,halo_tag,NNODP,NNODPP) 
    ! This sub does not work for UZAWA=1 and COUPLE=1 *********** 
    INTEGER NNODP, NNODPP
    REAL RELAX
    INTEGER CONMAS,D3,NOITS,PRETYM
    INTEGER FREDOP,NONODS
    !       PARAMETER(TOTELE=15504)
    !       PARAMETER(FREDOP=15504)
    !       PARAMETER(NONODS=17280)
    !       PARAMETER(CONMAS=1)
    ! CONMAS=.0 mass is lumped, CONMAS=:1. mass is consistent. 
    ! NB UZAWA use Uzawas method of pressure determination.  
    ! D3 then solution is truely 3-D. 
    ! COUPLE then the momentum equations are coupled. 
    ! D2HALF then we have a two and a half D simulation in (r,z,theta)coords.
    INTEGER NCT,NCOLM
    REAL RHS(FREDOP),ML(NONODS),DT 
    REAL C1T(NCT),C2T(NCT),C3T(NCT)
    ! NCT=4*FREDOP for simple element. 
    REAL DM1(NONODS),DM2(NONODS)
    REAL DM3(NONODS)
    REAL U(NONODS),V(NONODS),W(NONODS)
    REAL VECX(NONODS),VECY(NONODS),VECZ(NONODS)
    REAL SAVM(NONODS)
    REAL CONTRI(NONODS)
    ! R is only used if CONMAS is on.
    INTEGER I,J,FINDCT(FREDOP+1),COLCT(NCT),POS
    INTEGER CENTRM(NONODS),FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER PARA,halo_tag
    LOGICAL MASSYM

    ewrite(3,*)  "SUBROUTINE GETRHS()"
    !
    ewrite(3,*)'IN GETRHS ', D3
    if (IsParallel()) then
      CALL HALGET(U,NONODS,NONODS,NNODP,halo_tag)
      CALL HALGET(V,NONODS,NONODS,NNODP,halo_tag)
      IF(D3==1) THEN
         CALL HALGET(W,NONODS,NONODS,NNODP,halo_tag)
      ENDIF
    end if
    !     This subroutine actually forms the Right Hand Side of Poissons
    !     pressure equation
    do  I=1,FREDOP! Was loop 401
       RHS(I)=0.
    end do ! Was loop 401
    !     
    IF(CONMAS.EQ.0) THEN
       if (IsParallel()) then
          CALL HALGET(U,NONODS,NONODS,NNODP,halo_tag)
          CALL HALGET(V,NONODS,NONODS,NNODP,halo_tag)
       end if

       IF(D3.EQ.0) THEN
          do I=1,FREDOP
             !     DO 400 I=1,NNODPP
             !     The required divergence is placed in RHS in the subroutine HART.
             do POS=FINDCT(I),FINDCT(I+1)-1
                J=COLCT(POS)
                RHS(I) = RHS(I)&
                     - C1T(POS)*(VECX(J)/(ML(J)+DM1(J)) + U(J)/DT)&
                     - C2T(POS)*(VECY(J)/(ML(J)+DM2(J)) + V(J)/DT)
                !     We can subtract cty from RHS.
             END DO
          END DO
          ewrite(3,*) 'R2NORM(RHS,NNODPP,0):',R2NORM(RHS,NNODPP,0)
       ENDIF
       IF(D3.EQ.1) THEN
          if (IsParallel()) then
             CALL HALGET(W,NONODS,NONODS,NNODP,halo_tag)
          end if
          do I=1,FREDOP
             !     DO 410 I=1,NNODPP
             do POS=FINDCT(I),FINDCT(I+1)-1
                J=COLCT(POS)
                RHS(I) = RHS(I)&
                     - C1T(POS)*(VECX(J)/(ML(J)+DM1(J)) + U(J)/DT)&
                     - C2T(POS)*(VECY(J)/(ML(J)+DM2(J)) + V(J)/DT)&
                     - C3T(POS)*(VECZ(J)/(ML(J)+DM3(J)) + W(J)/DT)
             END DO
          END DO
       ENDIF
    ENDIF

    ewrite(1,*)'LEAVING GETRHS'
  END SUBROUTINE GETRHS
    
end module getrhs_module
