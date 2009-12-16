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
!
!
!
module advection_diffusion_cg_2d
  use fldebug

  implicit none

contains

  SUBROUTINE HART2D(U,NU,NV,UG,VG,&
      &          SOURCX,X,Y,&
      &          VECX,&
      &          C1T,C2T, BIGM, ML,&
      &          FINDCT,COLCT,NCT,FREDOP,TOTELE, &
      &          FINDRM,COLM,NCOLM,NONODS,CENTRM, &
      &          M,N,NLX,NLY, WEIGHT, NLOC,NGI,MLOC, &
      &          DISOPT,DISOPN,NDISOT,DT,THETA,BETA,LUMP,MAKSYM,&
      &          ABSLUM,SLUMP,&
      &          GETC12,INMAT,MLCENT, &
      &          NDGLNO,PNDGLN, &
      &          SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
      &          DENPT,MUPTXX,MUPTYY, ABSORB)
    ! This subroutine forms the discretised momentum or a field equation. 
    ! The discretization scheme is controlled through DISOPT(DISCRETIZATION OPTION) 
    ! DISOPT .LT. 0 then optimaly upwinded. 
    ! The following are absolute values of DISOPT...
    ! DISOPT=1 - Petrof Galerkin x-t.
    ! DISOPT=2 - Petrof Galerkin x-t(BUT weighted to get correct diffusion at steady state).
    ! DISOPT=3 - Petrof Galerkin (x,y).
    ! DISOPT=4 - Least squares (x,y,t). (symmetric)
    ! DISOPT=5 - Least squares (x,y).
    ! DISOPT=6 - Least squares THETA-method. (symmetric)
    ! DISOPT=7 - Galerkin THETA-method.(symmetric when MAKSYM=.true.)
    ! DISOPT=8 - Space weighting. 
    ! For DISOPT=7 ...
    !    If THETA=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
    !    If THETA=1.0 then bakward Euler is used in NS equns. THETA=2/3 Galerkin. 
    !    Similarly for the temperature equation.
    ! if NDISOT.GE.10 then set advection to zero and treat advection 
    ! using the high res method.
    ! IF MAKSYM then make BIGM symmetrix I.E take advection terms out of matrix.
    ! IF LUMP then lump the mass matrix else dont. 
    ! IF MLCENT then put ML to the diagonal of BIGM matrix. 
    ! BETA controls the conservative discretisation 
    ! BETA=0. is the usual advection form, BETA=1. is divergence form. 
    ! BETA=0.5 is an average of the divergence and advection form.(recommended).
    ! ABSORB =absorbtion of energy etc. 
    ! MUPTXX,MUPTYY =components of diffusivities. 
    ! NU,NV,NW are for the non-linear terms ordeneraly NU=U,NV=V,NW=W.
    ! UG,VG,WG are the grid velocities. 
    ! NDGLNO=element pter for unknowns. 
    ! SONDGL=element pter for materials and sources. 
    ! VONDGL=element pter for advection NU,UG.
    ! XONDGL=element pter for coordinates(x,y,z)
    ! If INMAT then return matricie(s). 
    ! -------------------------------------------------------------------------------------------
    ! If we are solving for another variable like temperature 
    ! or chemical species then NU,NV,NW will be the velocities 
    ! and U=TEMPERATURE OR CHEMICAL SPECIES. 
    ! -------------------------------------------------------------------------------------------
    ! NB SOURCE, DENPT, viscosity are SNONOD in length and have ele-pt SONDGL.
    INTEGER NCT,NCOLM,NLOC,NGI
    REAL RWIND
    PARAMETER(RWIND=0.5)
    INTEGER SNONOD,VNONOD,XNONOD
    INTEGER MLOC,DISOPT,DISOPN,NDISOT,DISCRA
    REAL MUPTXX(SNONOD),MUPTYY(SNONOD)
    REAL DENPT(SNONOD)
    REAL ABSORB(SNONOD)
    REAL DT,THETA,BETA
    INTEGER TOTELE,NONODS,FREDOP
    REAL U(NONODS)
    REAL UG(VNONOD),VG(VNONOD)
    REAL NU(VNONOD),NV(VNONOD)
    REAL MASS,MATRIX,RHSMAT,SOUMAT,VLMGI,VLKGI
    REAL VF1M,VF1,VF2M,VF2,VF3M,VF3
    REAL VMASS
    REAL THETA1,THETA2,THETA3,AHAT1,L1,ALPHA,THTWDT,INDT,TWOTHI
    REAL PE,HOVERQ,PWEIGI
    REAL GALERK,LEASQR,NLEAXT
    !Local area (used to be passed in from advdif
    real UD(NGI),VD(NGI)
    REAL MUGIXX(NGI),MUGIYY(NGI)
    REAL DENGI(NGI),DENGI2(NGI),L1GI(NGI)
    REAL AGI,BGI,CGI,DGI
    REAL DETJ,DETWEI(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI)
    ! HX,HY-characteristic length scales in x,y directions.
    REAL HXGI,HYGI
    REAL AHAT(NGI),GAMMA(NGI),MASCON(NGI)
    REAL BECTWE(NGI)
    !
    REAL VECX(NONODS)
    REAL BIGM(NCOLM)
    REAL ML(NONODS)
    REAL C1T(NCT),C2T(NCT)
    REAL X(XNONOD),Y(XNONOD)
    REAL SOURCX(SNONOD)
    INTEGER NDGLNO(TOTELE*NLOC),PNDGLN(TOTELE*MLOC)
    INTEGER SONDGL(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
    INTEGER XONDGL(TOTELE*NLOC)
    INTEGER ELE,ILOC,JLOC,L
    INTEGER GI,COUNT,I,GLOBI,GLOBJ,GLOBJS,GLOBIS
    INTEGER POS,LOWER,UPPER,PGLOBJ
    INTEGER FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER CENTRM(NONODS)
    INTEGER FINDCT(FREDOP+1),COLCT(NCT)
  
    REAL M(MLOC,NGI)
    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI)
    REAL WEIGHT(NGI)
    LOGICAL GETC12,MLCENT
    LOGICAL LUMP,MAKSYM,UPWIND,OPWIND,XTPETV
    LOGICAL ABSLUM,SLUMP
    LOGICAL INMAT
    REAL ABLUMP,RSLUMP,R,RUD,RVD,R1,R2,RNN,RR
    INTEGER IGLS,IGLV,IGLX,INUM
    REAL RLUMP,RSYM,RCONSI,RGAMMA,GIVOLF
    !
  
    ewrite(3,*)'In HART2D NONODS,TOTELE,DISOPT',&
        &                         NONODS,TOTELE,DISOPT
    !       ewrite(3,*) 'In HART2D NONODS',NONODS
    !
    ! This subroutine forms a contabution to the Right Hand Side
    ! of Poissons pressure equation, as well as  F1 & F2.
    ! THIS SUBROUTINE FORMS THE MATRICES BIGM,BIGT AND C1T,C2T.
    !
    !        ewrite(3,*) 'SOURCX:',SOURCX
    !        ewrite(3,*) 'SOURCY:',SOURCY
    !       ewrite(3,*) 'x=',x
    !       ewrite(3,*) 'y=',y
    !       stop
    !
    do  I=1,NONODS! Was loop 1370
      ML(I)  =0.
      VECX(I)=0.
    end do ! Was loop 1370
    IF(INMAT) THEN
      do  I=1,NCOLM! Was loop 1380
          BIGM(I)=0.
      end do ! Was loop 1380
    ENDIF
    IF(GETC12) THEN
      do  I=1,NCT! Was loop 20
          C1T(I)=0.
          C2T(I)=0.
      end do ! Was loop 20
    ENDIF

    RLUMP=0.
    RSYM=1.
    IF(LUMP) RLUMP=1.
    ! If MAKSYM then force the matrix to be symmetric 
    ! by treating advection explicitly. 
    IF(MAKSYM) RSYM=0.
    RCONSI=1.-RLUMP
    !
    ABLUMP=0.
    RSLUMP=0.
    IF(ABSLUM) ABLUMP=1.
    IF(SLUMP)  RSLUMP=1.
    !
    ! *********************************************
    INDT=1./DT
    TWOTHI=2./3.
    ! The defaults...
    GALERK=0.0
    LEASQR=0.0
    ALPHA =1.0
    L1    =1.0
    AHAT1 =1.0
    THETA1=1.0
    THETA2=1.0
    THETA3=1.0
    UPWIND=.FALSE.
    OPWIND=.FALSE.
    ! OPWIND=optimal upwinding method. 
    XTPETV=.FALSE.
    RGAMMA=1.0
    ! DISOPT=discretization option 
    DISCRA=ABS(DISOPT)
    IF(DISOPT.LT.0) OPWIND=.TRUE.
    IF(DISCRA.EQ.1) THEN
      ! Petrof Galerkin x-t
      XTPETV=.TRUE.
      UPWIND=.TRUE.
    ENDIF
    IF(DISCRA.EQ.2) THEN
      ! Petrof Galerkin x-t(BUT weighted to get correct diffusion at steady state)
      UPWIND=.TRUE.
    ENDIF
    IF(DISCRA.EQ.3) THEN
      ! Petrof Galerkin (x,y)
      AHAT1=0.
      UPWIND=.TRUE.
    ENDIF
    IF(DISCRA.EQ.4) THEN
      ! Least squares (x,y,t)
      THETA1=2./3.
      L1    =6./DT
      AHAT1 =-3.0
      LEASQR=1.0
      RGAMMA=1.0
    ENDIF
    IF(DISCRA.EQ.5) THEN
      ! Least squares (x,y)
      AHAT1=0.
      L1   =0.
    ENDIF
    IF(DISCRA.EQ.6) THEN
      ! Least squares THETA-method
      THETA3=THETA*3./2.
      AHAT1 =0.
      L1    =1.0/(DT*THETA)
      LEASQR=1.0
      RGAMMA=1.0
    ENDIF
    IF(DISCRA.EQ.7) THEN
      ! Galerkin THETA-method. 
      GALERK=1.0
      AHAT1 =-DT*(3.*THETA-2.0)
      L1    =-3.+6.*THETA
      RGAMMA=0.
    ENDIF
    IF(DISCRA.EQ.8) THEN
      ! Space weight method. 
      THETA1=0.
      THETA2=1.
      THETA3=DT*THETA*3./2.
      AHAT1=0.
      UPWIND=.TRUE.
    ENDIF
    !
    do  GI=1,NGI! Was loop 29
      AHAT(GI) =AHAT1
      GAMMA(GI)=RGAMMA
    end do ! Was loop 29
    ! *********************************************
    !
    THTWDT=THETA3*DT*2./3.
    !
    do  ELE=1,TOTELE! Was loop 340
      !
      !         ewrite(3,*)'ele=',ele
      !
      do  GI=1,NGI! Was loop 331
          !
          AGI=0.
          BGI=0.
          CGI=0.
          DGI=0.
          !
          UD(GI)=0.
          VD(GI)=0.
          DENGI(GI)=0.
          MUGIXX(GI)=0.
          MUGIYY(GI)=0.
          !
          GIVOLF=0.
          !
          do  L=1,NLOC! Was loop 79
            IGLS=SONDGL((ELE-1)*NLOC+L)
            IGLV=VONDGL((ELE-1)*NLOC+L)
            IGLX=XONDGL((ELE-1)*NLOC+L)
            !         ewrite(3,*)'IGLS,IGLV,IGLX',igls,IGLV,IGLX
            !         ewrite(3,*)'X(IGLX),Y(IGLX):',X(IGLX),Y(IGLX)
            AGI=AGI + NLX(L,GI)*X(IGLX) 
            BGI=BGI + NLX(L,GI)*Y(IGLX) 
            CGI=CGI + NLY(L,GI)*X(IGLX) 
            DGI=DGI + NLY(L,GI)*Y(IGLX) 
            !
            !
            IF(NDISOT.LT.10) THEN
                UD(GI)=UD(GI) + N(L,GI)*(NU(IGLV)-UG(IGLV))
                VD(GI)=VD(GI) + N(L,GI)*(NV(IGLV)-VG(IGLV))
            ENDIF
            DENGI(GI) =DENGI(GI) + N(L,GI)*DENPT(IGLS)
            ! The diffusivities at Gaiss pts. 
            MUGIXX(GI)=MUGIXX(GI)+ N(L,GI)*MUPTXX(IGLS)
            MUGIYY(GI)=MUGIYY(GI)+ N(L,GI)*MUPTYY(IGLS)
          end do ! Was loop 79
          !
          DETJ      =AGI*DGI-BGI*CGI
          DETWEI(GI)=abs(DETJ)*WEIGHT(GI)
          !
          !         ewrite(3,*)'GIVOLF=',GIVOLF
          ! Introduce grid vels in non-linear terms.
          !
          R=0.
          HXGI=0.
          HYGI=0.
          do  L=1,NLOC! Was loop 373
            NX(L,GI)=(DGI*NLX(L,GI)-BGI*NLY(L,GI))/DETJ
            NY(L,GI)=(-CGI*NLX(L,GI)+AGI*NLY(L,GI))/DETJ
            !
            HXGI=HXGI+ABS(NX(L,GI))
            HYGI=HYGI+ABS(NY(L,GI))
            IGLS=SONDGL((ELE-1)*NLOC+L)
            IGLV=VONDGL((ELE-1)*NLOC+L)
            IF(NDISOT.LT.10) THEN
                R=R&
                    &      +BETA*(NX(L,GI)*NU(IGLV)+NY(L,GI)*NV(IGLV) )*DENGI(GI)
            ENDIF
            R=R&
                  &      +N(L,GI)*ABSORB(IGLS)
          end do ! Was loop 373
          !
          !           BECTWE(GI)=R*DETWEI(GI)*DENGI(GI)
          BECTWE(GI)=R
          DENGI2(GI)=GALERK + (1.0-GALERK)*DENGI(GI)
          L1GI(GI)=L1+LEASQR*BECTWE(GI)
          !
          ! Work out GAMMA at each Gauss pt(OR centre of element)
          IF(UPWIND) THEN 
            HXGI=1./(0.5*HXGI)
            HYGI=1./(0.5*HYGI)
            HXGI=ABS( DGI*UD(GI)-BGI*VD(GI))/DETJ
            HYGI=ABS(-CGI*UD(GI)+AGI*VD(GI))/DETJ
            ! HX,HY are the characteristic length scales in x,y directions. 
            RUD=MAX(ABS(UD(GI)),1E-7)
            RVD=MAX(ABS(VD(GI)),1E-7)
            ! HOVERQ is the characteristic length scale divided by the speed. 
            ! If  XTPETV  then (x,y,t) Petrov Galerkin Method  else (x,y).
            IF(XTPETV) THEN
                !                     HOVERQ=MIN(HXGI/RUD, HYGI/RVD, DT)
                HOVERQ=MIN(2./MAX(HXGI,HYGI,1.E-7), DT)
            ELSE
                !                     HOVERQ=MIN(HXGI/RUD, HYGI/RVD)
                HOVERQ=2./MAX(HXGI,HYGI,1.E-7)
            ENDIF
            IF(OPWIND) THEN 
                ! (Near optimal choice of upwind parameter) PE=grid Peclet no. 
                ! The PE number here does not take into account anisotropic MU.
                PE=0.5*DENGI(GI)*(RUD**2+RVD**2)*HOVERQ &
                    &             / MAX(MUGIXX(GI),1.E-7) 
                ALPHA=1.0-1.0/MAX(1.0,PE)
            ENDIF
            GAMMA(GI)=RWIND*ALPHA*HOVERQ 
            AHAT(GI) =AHAT1*GAMMA(GI)
          ENDIF
          ! NB.  THTWDT=THETA3*DT*2./3.     
          MASCON(GI) = THETA1*2.*AHAT(GI)*INDT*DENGI(GI)&
              &      +THETA2*(AHAT(GI)*BECTWE(GI)+L1GI(GI)*DENGI(GI))&
              &      +THTWDT*L1GI(GI)*BECTWE(GI)       
      end do ! Was loop 331
      !
      do  ILOC=1,NLOC! Was loop 350
          GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
          GLOBIS=SONDGL((ELE-1)*NLOC+ILOC)
          do  JLOC=1,NLOC! Was loop 360
            GLOBJ =NDGLNO((ELE-1)*NLOC+JLOC)
            GLOBJS=SONDGL((ELE-1)*NLOC+JLOC)
            !
            MASS=0.
            MATRIX=0.
            RHSMAT=0.
            SOUMAT=0.     
            ! 
            do  GI=1,NGI! Was loop 458
                R1= UD(GI)*NX(ILOC,GI)+VD(GI)*NY(ILOC,GI)
                R2=(UD(GI)*NX(JLOC,GI)+VD(GI)*NY(JLOC,GI) )*DENGI(GI)
                !
                !  N.B.   BECTWE(GI) = BETA*CTY(GI)*DENGI(GI) .   
                !
                RNN  =N(ILOC,GI)*N(JLOC,GI)
                VLMGI=RNN*DETWEI(GI)
                VLKGI=NX(ILOC,GI)*NX(JLOC,GI)*MUGIXX(GI)&
                    &              +NY(ILOC,GI)*NY(JLOC,GI)*MUGIYY(GI)
                !
                VF1M =AHAT(GI)*(N(ILOC,GI)*R2*RSYM + VLKGI)
                !               VF1  =AHAT(GI)*(N(ILOC,GI)*R2      + VLKGI)
                !
                VF2  =GAMMA(GI)*R1*N(JLOC,GI)  
                VF2M =VF2*DENGI(GI)
                !
                RR =GAMMA(GI)*R1*(R2+BECTWE(GI)*N(JLOC,GI))&
                    &          +L1GI(GI)*N(ILOC,GI)*R2
                VF3M=RR*RSYM + L1GI(GI)*VLKGI
                !               VF3 =RR      + L1GI(GI)*VLKGI
                ! All the mass matrix contributions...       
                MASS = MASS + MASCON(GI)*VLMGI
                ! Matrix contributions (Excluding mass matrix)...
                MATRIX=MATRIX +  (THETA2*(VF1M + VF2M*RSYM)&
                    &               +THTWDT*VF3M)*DETWEI(GI) 
                ! NB.  THTWDT=THETA3*DT*2./3.
                ! RHSMAT contributions.
                !             RHSMAT=RHSMAT +(2.*DT*VF1 + VF3 
                !     &       + (2.*AHAT(GI)*INDT + L1GI(GI))
                !     &          *BECTWE(GI)*RNN )*DETWEI(GI)
                !
                RHSMAT=RHSMAT +( (2.*AHAT(GI)*INDT + L1GI(GI))&
                    &    *(BECTWE(GI)*RNN + N(ILOC,GI)*R2 + VLKGI)&
                    &    +GAMMA(GI)*R1*(BECTWE(GI)*N(JLOC,GI)+R2) )*DETWEI(GI)
                ! Source Material. 
                SOUMAT=SOUMAT + (2.*INDT*AHAT(GI) + L1GI(GI))*VLMGI &
                    &                       + VF2*DETWEI(GI)
            end do ! Was loop 458
            !
            ML(GLOBI)=ML(GLOBI)+MASS
            !          ewrite(3,*) 'mass,soumat:',mass,soumat
            !
            VECX(GLOBI)=VECX(GLOBI) - RHSMAT*U(GLOBJ) &
                  &                +SOUMAT*SOURCX(GLOBJS)*(1.-RSLUMP) &
                  &                +SOUMAT*SOURCX(GLOBIS)*RSLUMP  
            !
            ! we will now place contributions into the matrices BIGM and BIGT
            ! these are can be symm matrices.
            ! Find count. ***************************************
            IF(INMAT) THEN
                ! lump the mass matrix.(LUMP=1. if we are to LUMP the mass matrix)
                BIGM(CENTRM(GLOBI))=BIGM(CENTRM(GLOBI))+MASS*RLUMP
                !
                LOWER=FINDRM(GLOBI) 
                UPPER=FINDRM(GLOBI+1)-1
7000            CONTINUE
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
9000            CONTINUE
                ! Find count. ***************************************
                !
                !             COUNT=POSBM(ELE,(ILOC-1)*8+JLOC)
                !             DO 564 COUNT=FINDRM(GLOBI),FINDRM(GLOBI+1)-1
                !               IF(COLM(COUNT).EQ.GLOBJ) 
                !     &            BIGM(COUNT)=BIGM(COUNT)+MC1 + VLMD*(LUMP-1)
                !564          CONTINUE
                !
                BIGM(COUNT)=BIGM(COUNT)+MATRIX + MASS*RCONSI
                !
                ! end of IF(INMAT) THEN
            ENDIF
            !
          end do ! Was loop 360
          ! have gravity in -ve y direction.
          !         NIDEN=(N(ILOC,1)*DETJ(1)+N(ILOC,2)*DETJ(2)
          !     &        +N(ILOC,3)*DETJ(3)+N(ILOC,4)*DETJ(4)
          !     &        +N(ILOC,5)*DETJ(5)+N(ILOC,6)*DETJ(6)
          !     &        +N(ILOC,7)*DETJ(7)+N(ILOC,8)*DETJ(8))*DEN
          !           VECY(GLOBI)=VECY(GLOBI)-GRAVTY*NIDEN
          !
          ! *** THE FOLLOWING FORMS THE MATRICES C1T AND C2T ***
          IF(GETC12) THEN
            !
            do  JLOC=1,MLOC! Was loop 123
                !
                PGLOBJ=PNDGLN( (ELE-1)*MLOC +JLOC)
                ! Find count. ***************************************
                LOWER=FINDCT(PGLOBJ) 
                UPPER=FINDCT(PGLOBJ+1)-1
7001            CONTINUE
                INUM=LOWER+(UPPER-LOWER+1)/2 
                IF(GLOBI.GE.COLCT(INUM) )  THEN 
                  LOWER=INUM
                ELSE
                  UPPER=INUM
                ENDIF
                IF(UPPER-LOWER.LE.1) THEN
                  IF(GLOBI.EQ.COLCT(LOWER)) THEN
                      POS=LOWER
                  ELSE
                      POS=UPPER
                  ENDIF
                  GOTO 9001
                ENDIF
                GOTO 7001
9001            CONTINUE
                ! Find count. ***************************************
                do  GI=1,NGI! Was loop 123
                  PWEIGI=2.*AHAT(GI)*INDT+L1GI(GI)
                  !          ewrite(3,*) 'PWEIGI=',PWEIGI
                  C1T(POS)=C1T(POS)+PWEIGI*M(JLOC,GI)*NX(ILOC,GI)*DETWEI(GI)
                  C2T(POS)=C2T(POS)+PWEIGI*M(JLOC,GI)*NY(ILOC,GI)*DETWEI(GI)
                  !      C1T(POS)=1.
                end do ! Was loop 123
            end do ! Was loop 123
            !
          ENDIF
      end do ! Was loop 350
    end do ! Was loop 340
    !
    IF(MLCENT) THEN
      do  I=1,NONODS! Was loop 97
          ML(I)=BIGM(CENTRM(I))
      end do ! Was loop 97
    ENDIF
    !
    !
    ewrite(1,*)'exiting HART2D'
    !       ewrite(3,*) 'vecx:',vecx
    !       ewrite(1,*) 'exiting HART2D'
  END SUBROUTINE HART2D

end module advection_diffusion_cg_2d
