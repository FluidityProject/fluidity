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
module ctmult_module
  implicit none
contains
  
  SUBROUTINE CTMULT(RHS,U,V,W, DT, &
       ZERO,NDPSET, FREDOP,NONODS,D3,  &
       C1T,C2T,C3T,FINDCT,COLCT,NCT, &
                                !     Parallel ...
       PARA,halo_tag,&
       halo_tag_p,NNODP,NNODPP )
    use FLDebug
    use AllSorts
    use flcomms_module
    IMPLICIT NONE
    !     This sub works out RHS=RHS-(C1T*U+C2T*V)/DT
    INTEGER NDPSET, NCT, NNODP, NNODPP
    INTEGER FREDOP,NONODS
    REAL DT
    REAL RHS(FREDOP)
    REAL U(NONODS),V(NONODS),W(NONODS)
    INTEGER FINDCT(FREDOP+1),COLCT(NCT)
    REAL C1T(NCT),C2T(NCT),C3T(NCT)
    INTEGER POS
    INTEGER D3
    LOGICAL ZERO
    !     NNODP   no of vel nodes in subdomain excluding halo nodes.
    !     NNODPP no of pressure nodes in subdomain excluding halo nodes.
    INTEGER  PARA,halo_tag, halo_tag_p
    INTEGER I, J

    !     
    IF(ZERO) THEN
       do  I=1,FREDOP! Was loop 10
          RHS(I) = 0.
       end do ! Was loop 10
    ENDIF
    !     
    !     ewrite(3,*)'IN CTMULT R2NORM OF R2NORM(U,NNODP,PARA)=',
    !     &                                   R2NORM(U,NNODP,PARA)
    !     ewrite(3,*)'IN CTMULT R2NORM OF R2NORM(V,NNODP,PARA)=',
    !     &                                   R2NORM(V,NNODP,PARA)
    !     
    !     ewrite(3,*)'R2NORM(C1T,FINDCT(NNODPP+1)-1,PARA)=',
    !     &               R2NORM(C1T,FINDCT(NNODPP+1)-1,PARA)
    !     ewrite(3,*)'R2NORM(C2T,FINDCT(NNODPP+1)-1,PARA)=',
    !     &               R2NORM(C2T,FINDCT(NNODPP+1)-1,PARA)
    !     
    IF(PARA.EQ.1) CALL HALGET(U,NONODS,NONODS,NNODP,halo_tag)
    IF(PARA.EQ.1) CALL HALGET(V,NONODS,NONODS,NNODP,halo_tag)
    IF(D3.EQ.1) THEN 
       IF(PARA.EQ.1) CALL HALGET(W,NONODS,NONODS,NNODP,halo_tag)
    ENDIF
    !     ewrite(3,*)'IN-CTMULT R2NORM OF R2NORM(U,NNODP,PARA)=',
    !     &                                   R2NORM(U,NNODP,PARA)
    !     ewrite(3,*)'IN-CTMULT R2NORM OF R2NORM(V,NNODP,PARA)=',
    !     &                                   R2NORM(V,NNODP,PARA)
    !     
    IF(D3.EQ.0) THEN
       do I=1,FREDOP
          !     DO 400 I=1,NNODPP
          do POS=FINDCT(I),FINDCT(I+1)-1 
             J=COLCT(POS)
             RHS(I) = RHS(I)&
                  - ( C1T(POS)*U(J)&
                  + C2T(POS)*V(J) )/DT
             !     We can subtract cty from RHS.
          END DO
       END DO
    ENDIF
    !     This is for 3-D case. 
    IF(D3.EQ.1) THEN
       do I=1,FREDOP
          !     DO 411 I=1,NNODPP
          do POS=FINDCT(I),FINDCT(I+1)-1
             J=COLCT(POS)
             RHS(I) = RHS(I)&
                  -   (C1T(POS)*U(J)&
                  + C2T(POS)*V(J)&
                  + C3T(POS)*W(J)      )/DT
             !     We can subtract cty from RHS.
          END DO
       END DO
    ENDIF
    !     IF((NDPSET.GE.1).AND.(NDPSET.LE.FREDOP)) THEN
    !     RHS(NDPSET)=0.
    !     ENDIF
    !     
    IF(PARA.EQ.1) CALL HALGET(RHS,FREDOP,FREDOP,NNODPP,halo_tag_p)
    ewrite(3,*)'IN CTMULT R2NORM OF RHS=',R2NORM(RHS,NNODPP,PARA)
  END SUBROUTINE CTMULT
end module ctmult_module
