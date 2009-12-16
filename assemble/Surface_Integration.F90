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
module surface_integration

  use FLDebug
  use sdetnx_module
  
  implicit none

contains

  SUBROUTINE TEMSUF( INCBCS, PUTMAT, DT,THETA,&
       &     VEC,T,NONODS,SUFNOD,STOTEL,SNLOC,SNGI,&
       &     D3,DCYL,DSPH, &
       &     SNDGLN,TSNDGL, SALPHE,TAIRE, &
       &     SN,SNLX,SNLY,SWEIGH, &
       &     BIGM,NCOLM,COLM,FINDRM,  CENTRM,&
       &     LUMP,&
       &     X,Y,Z,R0,D0, &
       &     U,V,W)
    !     INCBCS switches on incomming b.c's
    !     SUFNOD=no of surface nodes. 
    !     This sub does the surface integration for a temperature 
    !     equation in any coordinate system. 
    !     We will now attempt to set k{\partial T\over\partial n}=T_{air}+\alpha T. 
    !     - thus defining SALPHE,TAIRE (surface element dependent) 
    !     If PUTMAT then put contribution into the matrix BIGM else do not. 
    !     The {1\over K} term is contained in \alpha (K=thermal conductivity) 
    !     in which T_{air} is the temp of the surroundings just outside the domain. 
    !     DT,THETA are the time step size and the time stepping parameter.
    !     VEC- r.h.s of temp eqn.
    !     T - temperature from previous time level.
    !     STOTEL,SNLOC,SNGI - no of surface elements, nodes per, Gauss pts per.
    !     D3,DCYL,DSPH - coord system. 
    !     SNDGLN - points to nodes (global) from surface elements. 
    !     SN,SNLX,SNLY,SWEIGH -surface element shape functions. 
    !     R0=the minimum distance from ocean bed to centre of earth.
    !     BIGM,NBIGM,COLM,FINDRM - temp matrix. 
    !     MUPTXX,MUPTYY,MUPTZZ - conductivities (associated with different directions)
    !     R0=the minimum distance from ocean bed to centre of earth.
    !     D0=H/r0 where H is the height of ocean above R0 (H not used here)
    !     For dimensional run D0=0.
    !     NB for cylindrical coords (r,z)           corresponds to (y,x)
    !     NB for spherical   coords (\psi,\theta,r) corresponds to (x,y,z)
    !     If LUMP then lump the mass matrix associated with incomming information.
    IMPLICIT NONE
    REAL PIE,INCLIM
    PARAMETER(PIE=3.141592654,INCLIM=1.E+15)
    REAL TWOPIE
    LOGICAL INCBCS
    INTEGER SUFNOD,NONODS,STOTEL,SNLOC,SNGI,NCOLM
    LOGICAL PUTMAT
    REAL DT,THETA
    REAL VEC(NONODS),T(NONODS)
    REAL TAIRE(SUFNOD),SALPHE(SUFNOD)
    REAL X(NONODS),Y(NONODS),Z(NONODS)
    REAL U(NONODS),V(NONODS),W(NONODS)
    REAL BIGM(NCOLM)
    INTEGER FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER CENTRM(NONODS)
    LOGICAL LUMP
    LOGICAL D3,DCYL,DSPH
    INTEGER SNDGLN(STOTEL*SNLOC),TSNDGL(STOTEL*SNLOC)
    REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
    REAL SWEIGH(SNGI)
    REAL SINX(SNGI),ISINX(SNGI),RMEM(SNGI),DETJ(SNGI)
    REAL DETWEI(SNGI),ALPH(SNGI)
    REAL UD(SNGI),VD(SNGI),WD(SNGI)
    REAL NORMXN(SNGI),NORMYN(SNGI),NORMZN(SNGI)
    !     Local variables...
    REAL RGI,THETGI,SINXGI
    REAL A,B,C
    REAL RBIGM1,RBIGM2
    !     
    INTEGER ELE,L,GI,Q,P,GLOBI,GLOBJ,COUNT,LOWER,UPPER
    INTEGER LEGL,LIGL,LJGL
    INTEGER IGL,IGLT,ILOC,JLOC,ITNOD,JTNOD,INUM,INORM
    REAL NN, DXDLX,DXDLY,DYDLX,DYDLY,DZDLX,DZDLY
    REAL NNX,NNY,NNZ,RNN
    REAL R0,D0,DD0
    REAL NORMX,NORMY,NORMZ,NORM,SAREA,MDOTN
    REAL XCS,YCS,ZCS,XCV,YCV,ZCV
    LOGICAL ONEINC,ONENOR

    DD0=1.
    IF(D0.NE.0) DD0=D0
    !     
    TWOPIE=1.0
    IF(DCYL) TWOPIE=2.*PIE
    !     
    do  ELE=1,STOTEL! Was loop 340
       !     ewrite(3,*) 'ELE=',ELE
       !     
       !     
       !     Work out what sort of b.c to to have for surface element...
       ONEINC=.FALSE.
       ONENOR=.FALSE.
       do  L=1,SNLOC! Was loop 479
          ITNOD=TSNDGL((ELE-1)*SNLOC+L)
          IF(SALPHE(ITNOD).GT.INCLIM) THEN
             ONEINC=.TRUE.
          ELSE
             ONENOR=.TRUE.
          ENDIF
       end do ! Was loop 479
       !     
       IF(.NOT.INCBCS) THEN
          ONENOR=.TRUE.
          ONEINC=.FALSE.
       ENDIF
       !     
       IF(ONEINC.AND.(.NOT.ONENOR)) THEN
          !     Calculate an approx normal (normx,normy,normz)...
          CALL SIMPNORM(NORMX,NORMY,NORMZ,D3,&
               &           SNDGLN,STOTEL,SNLOC,NONODS,NONODS,ELE,&
               &           X,Y,Z,&
               &           NCOLM,FINDRM,COLM) 
       ELSE
          NORMX=1.0
          NORMY=0.0
          NORMZ=0.0
       ENDIF
       !     
       CALL SDETNX2(ELE,SNDGLN,&
            &        STOTEL,NONODS,SNLOC,SNGI, &
            &        X,Y,Z,&
            &        SN,SNLX,SNLY, SWEIGH, DETWEI,SAREA, D3,DCYL, &
            &        NORMXN,NORMYN,NORMZN,&
            &        NORMX,NORMY,NORMZ)
       !     
       do  GI=1,SNGI! Was loop 331
          ALPH(GI)=0.
          UD(GI)=0.
          VD(GI)=0.
          WD(GI)=0.
          do  L=1,SNLOC! Was loop 279
             IGL  =SNDGLN((ELE-1)*SNLOC+L)
             ITNOD=TSNDGL((ELE-1)*SNLOC+L)
             ALPH(GI)=ALPH(GI)+SN(L,GI)*SALPHE(ITNOD)
             UD(GI)=UD(GI)+SN(L,GI)*U(IGL)
             VD(GI)=VD(GI)+SN(L,GI)*V(IGL)
             IF(D3) WD(GI)=WD(GI)+SN(L,GI)*W(IGL)
          end do ! Was loop 279
       end do ! Was loop 331
       !     
       !     NB these are surface elements       
       !     
       !     
       do  ILOC=1,SNLOC! Was loop 350
          GLOBI=SNDGLN((ELE-1)*SNLOC+ILOC)
          ITNOD=TSNDGL((ELE-1)*SNLOC+ILOC)
          !     ewrite(3,*) 'iloc,GLOBI,ITNOD:',iloc,GLOBI,ITNOD
          do  JLOC=1,SNLOC! Was loop 350
             GLOBJ=SNDGLN((ELE-1)*SNLOC+JLOC)
             JTNOD=TSNDGL((ELE-1)*SNLOC+JLOC)
             !     ewrite(3,*) 'jloc,globj:',jloc,globj
             !     LJGL=LGL(L)
             !     
             RBIGM1=0.
             RBIGM2=0.
             do  GI=1,SNGI! Was loop 458
                RNN=SN(ILOC,GI)*SN(JLOC,GI)*DETWEI(GI)
                !     
                IF(ONENOR.AND.(.NOT.ONEINC)) THEN
                   !     Robin b.c...
                   RBIGM1=RBIGM1+RNN
                   RBIGM2=RBIGM2+RNN*ALPH(GI)
                ENDIF
                !     
                IF(ONEINC.AND.(.NOT.ONENOR)) THEN
                   !     Add surface integral only for incomming information...
                   MDOTN=UD(GI)*NORMXN(GI)+VD(GI)*NORMYN(GI)+WD(GI)*NORMZN(GI)
                   IF(MDOTN.LT.0.0) THEN
                      RBIGM1=RBIGM1-MDOTN*RNN
                      RBIGM2=RBIGM2+MDOTN*RNN
                   ENDIF
                ENDIF
             end do ! Was loop 458
             !     
             IF(LUMP) THEN
                VEC(GLOBI)=VEC(GLOBI)&
                     &                 +RBIGM1*TAIRE(JTNOD) + RBIGM2*T(GLOBI)  
                COUNT=CENTRM(GLOBI)
                BIGM(COUNT)=BIGM(COUNT) -DT*THETA*RBIGM2
             ELSE
                VEC(GLOBI)=VEC(GLOBI)&
                     &                 +RBIGM1*TAIRE(JTNOD) + RBIGM2*T(GLOBJ)  
                !     &           +RBIGM *(T(GLOBJ)-TAIRN(LJGL)-TAIRE(ELE) )
                !     
                !     Find count. ***************************************
                IF(PUTMAT) THEN
                   LOWER=FINDRM(GLOBI) 
                   UPPER=FINDRM(GLOBI+1)-1
7000               CONTINUE
                   INUM=LOWER+(UPPER-LOWER+1)/2 
                   IF(GLOBJ.GE.COLM(INUM) )  THEN 
                      LOWER=INUM
                   ELSE
                      UPPER=INUM
                   ENDIF
                   IF(UPPER-LOWER.LE.1) THEN
                      IF(GLOBJ.EQ.COLM(LOWER)) THEN
                         COUNT=LOWER
                      ELSE
                         COUNT=UPPER
                      ENDIF
                      GOTO 9000
                   ENDIF
                   GOTO 7000
9000               CONTINUE
                   !     Find count. ***************************************
                   !     
                   !     ewrite(3,*) 'count=',count
                   !     ewrite(3,*) 'GLOBI,PUTMAT=',GLOBI,PUTMAT
                   !     ewrite(3,*) 'FINDRM(GLOBI),FINDRM(GLOBI+1):',
                   !     &              FINDRM(GLOBI),FINDRM(GLOBI+1)
                   BIGM(COUNT)=BIGM(COUNT) -DT*THETA*RBIGM2
                ENDIF
                !     ewrite(3,*) 'here 4'
                !     ENDOF IF(LUMP) THEN ELSE...
             ENDIF
             !     
          end do ! Was loop 350
       end do ! Was loop 350
       !     
    end do ! Was loop 340
    ewrite(3,*) 'finished temsuf'
    !     ewrite(3,*) 'SALPHE',SALPHE
    !     ewrite(3,*) 'nonods,sufnod:',nonods,sufnod
    !     STOP 393
    RETURN
  END SUBROUTINE TEMSUF
  !
  !
  !
  !
  SUBROUTINE SIMPNORM(NORMX,NORMY,NORMZ,D3,&
       &        SNDGLN,STOTEL,SNLOC,XNONOD,NONODS,ELE,&
       &        X,Y,Z,&
       &        NCOLM,FINDRM,COLM)
    ! Calculate an approx normal (normx,normy,normz)...
    IMPLICIT NONE
    INTEGER STOTEL,SNLOC,XNONOD,NONODS,ELE,NCOLM
    INTEGER SNDGLN(STOTEL*SNLOC)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    INTEGER FINDRM(NONODS+1),COLM(NCOLM)
    LOGICAL D3
    REAL NORMX,NORMY,NORMZ
    ! Local variables...
    REAL XCS,YCS,ZCS,XCV,YCV,ZCV,NORM
    INTEGER ICV,L,IGL,L2,IGL2,COUNT,COUNT2,NODJ,NODJ2
    LOGICAL ISURFA,FOUND,SURF,ALLSUF

    XCS=0.
    YCS=0.
    ZCS=0.
    do  L=1,SNLOC! Was loop 178
       IGL  =SNDGLN((ELE-1)*SNLOC+L)
       XCS=XCS+X(IGL)/REAL(SNLOC)
       YCS=YCS+Y(IGL)/REAL(SNLOC)
       IF(D3) ZCS=ZCS+Z(IGL)/REAL(SNLOC)
    end do ! Was loop 178
    !
    XCV=0.
    YCV=0.
    ZCV=0.
    ICV=0
    L=1
    IGL  =SNDGLN((ELE-1)*SNLOC+L)
    do COUNT=FINDRM(IGL),FINDRM(IGL+1)-1
       NODJ=COLM(COUNT)
       ! Make sure its not a surface node of surface element ELE
       SURF=.FALSE.
       do L2=1,SNLOC
          IGL2=SNDGLN((ELE-1)*SNLOC+L2)
          IF(IGL2.EQ.NODJ) SURF=.TRUE.
       END DO
       IF(.NOT.SURF) THEN
          ! make sure NODJ is connected to all surface nodes and is not a surface node for 
          ! surface element.
          ALLSUF=.TRUE.
          do L2=1,SNLOC
             IGL2=SNDGLN((ELE-1)*SNLOC+L2)
             FOUND=.FALSE.
             do COUNT2=FINDRM(NODJ),FINDRM(NODJ+1)-1
                NODJ2=COLM(COUNT2)
                IF(IGL2.EQ.NODJ2) FOUND=.TRUE.
             END DO
             IF(.NOT.FOUND) ALLSUF=.FALSE.
          END DO
          IF(ALLSUF) THEN
             XCV=XCV+X(NODJ)
             YCV=YCV+Y(NODJ)
             IF(D3) ZCV=ZCV+Z(NODJ)
             ICV=ICV+1
          ENDIF
          ! ENDOF IF(.NOT.SURF) THEN...
       ENDIF
    END DO
    !
    XCV=XCV/REAL(ICV)
    YCV=YCV/REAL(ICV)
    ZCV=ZCV/REAL(ICV)
    NORMX=XCS-XCV
    NORMY=YCS-YCV
    NORMZ=ZCS-ZCV
    NORM=MAX(1.E-15, SQRT(NORMX**2 + NORMY**2 + NORMZ**2)  )
    NORMX=NORMX/NORM
    NORMY=NORMY/NORM
    NORMZ=NORMZ/NORM
    RETURN

  end SUBROUTINE SIMPNORM
end module surface_integration
