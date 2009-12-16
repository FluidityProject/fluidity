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

module solve_momentum_cg

  use fldebug
  use global_parameters, only : new_options
  use parallel_tools
  use solmom_module
  use solvers
  use spud

  implicit none

  private

  public :: getvel

contains

  SUBROUTINE GETVEL(U,V,W,NU,NV,NW,IGUESS,FX,FY,FZ,DT, &
       &    DM1,DM2,DM3,BIGM1,BIGM2,BIGM3,&
       &    FINDRM,COLM,CENTRM,NCOLM,NBIGM, &
       &    SAVDIA,PHI,FL,PHIL,&
                                ! This NPHASE is for no of phases. 
       &    ILENG,NPHASE, NONODS,D3,COUPLE, &
       &    GMRESQ,&
       &     PARA,halo_tag,NNODP, &
                                ! Th following is for rotations only ROTAT=.TRUE.
       &         ROTAT,&
       &         velocity_option_path)
    ! This sub solves for velocity.
    ! If COUPLE=1 then SAVDIA,F,PHI are not used. 
    ! If COUPLE=0 then FL,PHIL are not used. 
    ! NPHASE= no of phases.
    INTEGER MATSTR
    PARAMETER(MATSTR=0)
    INTEGER D3,NONODS,NNODP,ILENG,NPHASE,COUPLE
    INTEGER NCOLM,NBIGM
    INTEGER UZAWA,GMRESQ
    REAL U(NONODS*NPHASE),V(NONODS*NPHASE),W(NONODS*NPHASE)
    ! most up to date vels...
    REAL NU(NONODS*NPHASE),NV(NONODS*NPHASE),NW(NONODS*NPHASE)
    LOGICAL IGUESS
    ! If IGUESS then have an initial guess. 
    REAL FX(NONODS*NPHASE),FY(NONODS*NPHASE),FZ(NONODS*NPHASE)
    REAL DM1(NONODS),DM2(NONODS),DM3(NONODS),DT
    REAL BIGM1(NBIGM),BIGM2(NBIGM),BIGM3(NBIGM)
    REAL SAVDIA(NONODS)
    REAL PHI(NONODS)
    REAL FL(ILENG),PHIL(ILENG)
    ! R contains temperory vectors for solver and has length 
    ! NR=MAX(ILENG,NONODS) * MAX(MULPA*2,UZAWA*5)
    ! ILENG=3*NONODS*NPHASE (if 3-D)
    ! ILENG=2*NONODS*NPHASE (if 2-D)
    !
    INTEGER FINDRM(NONODS+1),COLM(NCOLM),CENTRM(NONODS)
    ! For parallel processing ...
    ! NNODP=no of subdomain nodes excluding halo nodes. 
    INTEGER PARA,halo_tag

    !      options tree path to velocity field for solver options
    character(len=*), intent(in):: velocity_option_path
    LOGICAL ROTAT
    ! RTDR contains RT D R   ^T where D contains DMI's 
    ! and R is the rotation matrix. 
    ! Local variables ...
    INTEGER IP,IPLUS,IPLUS2,IPLUS3,KITS,NOD,NODV
    INTEGER I
    REAL RMAXV

    ewrite(3,*)'INSIDE GETVEL'
    ewrite(3,*) 'inside getvel couple,ROTAT:',couple,ROTAT

    IF((COUPLE.EQ.0).AND.(.NOT.ROTAT)) THEN
       ! This is for when the momentum equations are not coupled. 
       ewrite(3,*) 'solving for u: NNODP=',NNODP
       CALL SOLMOM(U,NU,IGUESS,FX,PHI,DT,&
            &      DM1,BIGM1,FINDRM,COLM,CENTRM,NCOLM, &
            &      SAVDIA,&
            &      NONODS,&
            &      GMRESQ,&
            &     PARA,halo_tag,NNODP,&
            &     velocity_option_path)
       !
       ewrite(3,*) 'solving for v: NNODP=',NNODP
       CALL SOLMOM(V,NV,IGUESS,FY,PHI,DT,&
            &      DM2,BIGM2,FINDRM,COLM,CENTRM,NCOLM, &
            &      SAVDIA,&
            &      NONODS,&
            &      GMRESQ,&
            &     PARA,halo_tag,NNODP,&
            &     velocity_option_path)
       !
       IF(D3.EQ.1) THEN
          ewrite(3,*) 'solving for w: NNODP=',NNODP
          CALL SOLMOM(W,NW,IGUESS,FZ,PHI,DT,&
               &      DM3,BIGM3,FINDRM,COLM,CENTRM,NCOLM, &
               &      SAVDIA,&
               &      NONODS,&
               &      GMRESQ,&
               &     PARA,halo_tag,NNODP,&
               &     velocity_option_path)
       ENDIF
    ENDIF
    !
    IF((COUPLE.EQ.1).OR.ROTAT) THEN
       IF(IGUESS) THEN
          IF(ILENG.EQ.3*NONODS*NPHASE) THEN
             do  IP=1,NPHASE! Was loop 614
                IPLUS =(IP-1)*NONODS
                IPLUS3=(IP-1)*3*NONODS
                do  I=1,NONODS! Was loop 614
                   PHIL(I         +IPLUS3)=(NU(I+IPLUS)-U(I+IPLUS))/DT
                   PHIL(I+NONODS  +IPLUS3)=(NV(I+IPLUS)-V(I+IPLUS))/DT
                   PHIL(I+2*NONODS+IPLUS3)=(NW(I+IPLUS)-W(I+IPLUS))/DT
                end do ! Was loop 614
             end do ! Was loop 614
          ELSE
             do  IP=1,NPHASE! Was loop 624
                IPLUS =(IP-1)*NONODS
                IPLUS2=(IP-1)*2*NONODS
                do  I=1,NONODS! Was loop 624
                   PHIL(I       +IPLUS2)= (NU(I+IPLUS)-U(I+IPLUS))/DT
                   PHIL(I+NONODS+IPLUS2)= (NV(I+IPLUS)-V(I+IPLUS))/DT
                end do ! Was loop 624
             end do ! Was loop 624
          ENDIF
       ELSE
          do  I=1,ILENG! Was loop 79
             PHIL(I)=0.
          end do ! Was loop 79
          ! ENDOF IF(IGUESS) THEN ELSE
       ENDIF
       ewrite(3,*)'INSIDE GETVEL HERE2'
       IF(ILENG.EQ.3*NONODS*NPHASE) THEN
          do  IP=1,NPHASE! Was loop 97
             IPLUS =(IP-1)*NONODS
             IPLUS3=(IP-1)*3*NONODS
             do  I=1,NONODS! Was loop 97
                FL(I         +IPLUS3)=FX(I+IPLUS)
                FL(I+NONODS  +IPLUS3)=FY(I+IPLUS)
                FL(I+2*NONODS+IPLUS3)=FZ(I+IPLUS)
             end do ! Was loop 97
          end do ! Was loop 97
       ELSE
          do  IP=1,NPHASE! Was loop 99
             IPLUS =(IP-1)*NONODS
             IPLUS2=(IP-1)*2*NONODS
             do  I=1,NONODS! Was loop 99
                FL(I       +IPLUS2)  = FX(I+IPLUS)
                FL(I+NONODS+IPLUS2)  = FY(I+IPLUS)
             end do ! Was loop 99
          end do ! Was loop 99
       ENDIF
       ewrite(3,*)'INSIDE GETVEL HERE1'
       !             ewrite(3,*) 'befor mulpas MULPA=',MULPA
       !             ewrite(3,*) 'r2norm(FL,ileng),r2norm(BIGM1,NBIGM):',
       !     &                r2norm(FL,ileng),r2norm(BIGM1,NBIGM)
       !             ewrite(3,*) '******************************'
       !             ewrite(3,*) '******************************'

       !
       IF(GMRESQ.NE.0) THEN
          ewrite(3,*)'INSIDE GETVEL just before GMRES1 ROTAT=',&
               &               ROTAT
          CALL GMRES(PHIL,FL,NONODS,ILENG,NNODP,.TRUE.,&
               &     BIGM1,FINDRM,COLM,NCOLM,NBIGM,&
               & halo_tag,KITS, &
               &         option_path=velocity_option_path)

       ENDIF
       !
       ! Work out the vels from the acceleration PHI. 
       IF(ILENG.EQ.3*NONODS*NPHASE) THEN
          do  IP=1,NPHASE! Was loop 610
             IPLUS =(IP-1)*NONODS
             IPLUS3=(IP-1)*3*NONODS
             do  I=1,NONODS! Was loop 610
                U(I+IPLUS)= PHIL(I         +IPLUS3)*DT + U(I+IPLUS)
                V(I+IPLUS)= PHIL(I+NONODS  +IPLUS3)*DT + V(I+IPLUS)
                W(I+IPLUS)= PHIL(I+2*NONODS+IPLUS3)*DT + W(I+IPLUS)
             end do ! Was loop 610
          end do ! Was loop 610
          !         print *,'fl(6441):',fl(6441)
          !         print *,'u(6441):',u(6441)
          !          print *,'bigm1(centrm(6441)):',bigm1(centrm(6441))
          !     stop 921
       ELSE
          do  IP=1,NPHASE! Was loop 620
             IPLUS =(IP-1)*NONODS
             IPLUS2=(IP-1)*2*NONODS
             do  I=1,NONODS! Was loop 620
                U(I+IPLUS)= PHIL(I       +IPLUS2) *DT + U(I+IPLUS)
                V(I+IPLUS)= PHIL(I+NONODS+IPLUS2) *DT + V(I+IPLUS)
             end do ! Was loop 620
          end do ! Was loop 620
       ENDIF
    ENDIF
    !     
    ewrite(3,*)'IN GET VEL U R2NORM=',R2NORM(U,NNODP,PARA)
    ewrite(3,*)'IN GET VEL V R2NORM=',R2NORM(V,NNODP,PARA)
    ewrite(3,*)'IN GET VEL W R2NORM=',R2NORM(W,NNODP,PARA)
    !     Find max V vel nod NODV
    RMAXV=-1.
    do  NOD=1,NONODS! Was loop 334
       IF(ABS(V(NOD)).GT.RMAXV) THEN
          NODV=NOD
          RMAXV=ABS(V(NOD))
       ENDIF
    end do ! Was loop 334
    !
    IF(PARA.EQ.1) CALL ALLMAX(RMAXV)

    RETURN
  END SUBROUTINE GETVEL

end module solve_momentum_cg
