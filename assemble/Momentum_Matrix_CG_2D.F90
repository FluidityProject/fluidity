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
module momentum_matrix_cg_2d
  use fldebug
  use limit_gradient
  implicit none

  private
  public :: diff2d

contains

  SUBROUTINE DIFF2D(U,V,NU,NV,UG,VG,&
      &          SOURCX,SOURCY,X,Y,&
      &          VECX,VECY,&
      &          C1T,C2T, BIGM,NBIGM,  ML,&
      &          FINDCT,COLCT,NCT,FREDOP, TOTELE,&
      &          FINDRM,COLM,NCOLM,NONODS,CENTRM, &
      &          M,N,NLX,NLY, WEIGHT, NLOC,NGI,MLOC, &
      &          DISOPT,DISOPN,DT,THETA,BETA,LUMP,MAKSYM,CONVIS,&
      &          GETC12, &
      &          ABSLUM,SLUMP,&
      &          NDGLNO,PNDGLN, &
      &          SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
      &          DENPT,&
      &          MUPTXX,MUPTXY,MUPTYY, XABSOR,YABSOR, &
      &          CVMVEC)
    ! This subroutine forms the discretised momentum or a field equation. 
    ! The discretization scheme is controlled through DISOPT(DISCRETIZATION OPTION) 
    ! All methods use a balancing diffusion term and THETA time stepping. 
    ! The following are absolute values of DISOPT...
    ! DISOPT=1 - balancing diffusion based on (x,y) space.
    ! DISOPT=2 - Laxwendrof balancing diffusion.
    ! DISOPT=3 - (x,y,t) -balancing diffusion. 
    ! DISOPT=4 - No balancing diffusion.
    ! DISOPT=5 - nonlinear streamline and cross stream diffusion.
    ! DISOPT=6 - nonlinear upwind in steapest direction.
    ! DISOPT=7 - nonlinear streamline cross stream diffusion(but restricted)
    !
    ! if DISOPN.ne.0 then set advection to zero and treat advection 
    ! using the high res method.
    ! SORN contains the div q field NOT N_i * div q.  (IS NOT USED NOW)
    ! N.B.  When DTLAX=0.5 *DT we have something similar to lax-wendrof scheme. 
    ! The Petrof Galerkin weighting is controlled through DTLAX. 
    ! IF MAKSYM then make BIGM symmetrix .
    ! IF LUMP then lump the mass matrix else dont. 
    ! BETA controls the conservative discretisation 
    ! BETA=0. is the usual advection form, BETA=1. is divergence form. 
    ! BETA=0.5 is an average of the divergence and advection form.(recommended).
    ! ABSORBX =absorbtion of energy of x-commentum equation. 
    ! CVMVEC=virtual mass vector (is 0 in non 2-phase flow). 
    ! MUPTXX,MUPTYY =components of diffusivities. 
    ! NU,NV,NW are for the non-linear terms ordeneraly NU=U,NV=V,NW=W.
    ! UG,VG,WG are the grid velocities. 
    ! NDGLNO=element pter for unknowns. 
    ! SONDGL=element pter for materials and sources. 
    ! VONDGL=element pter for advection NU,UG.
    ! XONDGL=element pter for coordinates(x,y,z)
    !
    ! If CONVIS then assume the conservation of viscosity due to div q term, 
    ! if div q is not zero then this must be set to .TRUE. 
    ! NU,NV,NW are for the non-linear terms ordeneraly NU=U,NV=V,NW=W.
    ! UG,VG,WG are the grid velocities. 
    ! NB SOURCE, DENPT, viscosity are SNONOD in length and have ele-pt SONDGL.
    ! THE BLOCKS OF THE MATRIX BIGM ARE ARRANGED AS FOLLOWS
    !   1 * 
    !   2 3        *- MATRIX NOT STORED.  
    ! OR FOR NON-BLOCK SYM ...     1  2
    !                              3  4.
    ! 
    INTEGER NBIGM,NCT,NCOLM,NLOC,NGI
    REAL RWIND
    PARAMETER(RWIND=0.5)
    INTEGER SNONOD,VNONOD,XNONOD
    INTEGER MLOC,DISOPT,DISOPN,DISCRA
    REAL DENPT(SNONOD)
    REAL MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTYY(SNONOD)
    REAL XABSOR(SNONOD),YABSOR(SNONOD)
    REAL CVMVEC(SNONOD)
    REAL DT,THETA,BETA
    INTEGER TOTELE,NONODS,FREDOP
    ! If THETA=0.5 then Crank-Nickolson is used in Navier Stokes equ's.
    ! If THETA=1.0 then bakward Euler is used in NS equns.
    ! Similarly for the temperature equation. 
    REAL U(NONODS),V(NONODS)
    REAL MC1,HL,VLND,VLM,VLMD
    REAL UG(VNONOD),VG(VNONOD)
    REAL NU(VNONOD),NV(VNONOD)
    !
    REAL UD(NGI),VD(NGI)
    REAL HOVERQ
    REAL AGI,BGI,CGI,DGI,DETJ,DETWEI(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI)
    REAL BECTWE(NGI)
    REAL DDETWE(NGI)
    REAL GIVOLF,GIVOLX,GIVOLY
    ! HX,HY-characteristic length scales in x,y directions.
    REAL HXGI,HYGI
    REAL GAMMA(NGI)
    REAL MXXDTW(NGI),MXYDTW(NGI),MYYDTW(NGI)
    REAL XABDGI(NGI),YABDGI(NGI)
    REAL MUGIXX,MUGIXY,MUGIYY
    REAL DENGI(NGI)
    !
    REAL VECX(NONODS),VECY(NONODS)
    !       REAL MAXU(NONODS),MINU(NONODS)
    !       REAL MAXV(NONODS),MINV(NONODS)
    REAL BIGM(NBIGM)
    ! This only catters for NBIGM=3*NCOLM at the moment. 
    REAL ML(NONODS)
    REAL C1T(NCT),C2T(NCT)
    REAL X(XNONOD),Y(XNONOD)
    REAL SOURCX(SNONOD),SOURCY(SNONOD)
    REAL VLKXX,VLKXY,VLKYX,VLKYY
    REAL VXVXY,VYVYX
    REAL VABSX,VABSY
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
    !       REAL GAMMAT(NGI),TUD(NGI),TVD(NGI)
    !       REAL LIMU1(NGI),LIMU2(NGI),LIMT1(NGI),LIMT2(NGI)
    REAL GAMMAU(NGI),GAMMAV(NGI)
    REAL UUD(NGI),UVD(NGI), VUD(NGI),VVD(NGI)
    !
    REAL TWOTHI
    REAL ABLUMP
    LOGICAL ABSLUM,SLUMP
    LOGICAL GETC12
    LOGICAL LUMP,MAKSYM,CONVIS,UPWIND,XTPETV
    LOGICAL BLKSYM
    !
    REAL LOVERQ,HL1,HL2 
    REAL R1U,R2U,R1V,R2V,HH2
    LOGICAL NONLIN
    INTEGER IBL11,IBL12,IBL21,IBL22
    INTEGER IGLS,IGLV,IGLX,ICENT,INUM
    REAL RLUMP,RSLUMP,RSYM,RGAMMA,R,R1,R2,RNN
    REAL VXVXX,VYVYY
  
  
    !        DISOPT=4
    !       call rclear(nu,nonods)
    !       call rclear(nv,nonods)
    ewrite(3,*)'In DIFF2D nonods,snonod,xnonod,vnonod',&
        &                         nonods,snonod,xnonod,vnonod
    !       DO 127 I=1,NONODS
    !            ewrite(3,*)'I,DENPT(I),MUPT(I),SOURCX(I):',
    !     &                     I,DENPT(I),MUPT(I),SOURCX(I)
    !127     CONTINUE
  
    ! This subroutine forms a contabution to the Right Hand Side
    ! of Poissons pressure equation, as well as  F1 & F2.
    ! THIS SUBROUTINE FORMS THE MATRICES BIGM,BIGT AND C1T,C2T.
    !
    ! THE STRUCTURE OF THE MATRIX ...
    IF(NBIGM.EQ.4*NCOLM) THEN
      BLKSYM=.FALSE.
    ELSE
      BLKSYM=.TRUE.
    ENDIF
    IF(BLKSYM) THEN
      IBL11=0
      IBL12=1*NCOLM
      IBL21=1*NCOLM
      IBL22=2*NCOLM
    ELSE
      IBL11=0
      IBL12=1*NCOLM
      IBL21=2*NCOLM
      IBL22=3*NCOLM
    ENDIF
    !
    !
    do  I=1,NONODS! Was loop 1370
      ML(I)=0.
      VECX(I)=0.
      VECY(I)=0.
    end do ! Was loop 1370
    do  I=1,NBIGM! Was loop 1380
      BIGM(I)=0.
    end do ! Was loop 1380
    IF(GETC12) THEN
      do  I=1,NCT! Was loop 20
          C1T(I)=0.
          C2T(I)=0.
      end do ! Was loop 20
    ENDIF
    !
    !           ABSLUM=.true.
    !           SLUMP=.TRUE.
    !
    RLUMP =0.
    ABLUMP=0.
    RSLUMP=0.
    RSYM  =1.
    TWOTHI=0.
    IF(LUMP)   RLUMP =1.
    IF(ABSLUM) ABLUMP=1.
    IF(SLUMP)  RSLUMP=1.
    IF(MAKSYM) RSYM  =0.
    IF(CONVIS) TWOTHI=2./3.
    !        ewrite(3,*) 'ABLUMP,RLUMP,RSLUMP,RSYM,THETA:',
    !     &           ABLUMP,RLUMP,RSLUMP,RSYM,THETA
    !
    !
    ! START OF OPTIONS ******************************
    UPWIND=.FALSE.
    NONLIN=.FALSE.
    ! OPWIND=optimal upwinding method. 
    XTPETV=.FALSE.
    RGAMMA=0.0
    ! DISOPT=discretization option 
    DISCRA=ABS(DISOPT)
    !           IF(DISOPT.LT.0) OPWIND=.TRUE.
    !
    UUD(1:NGI) = 0.0
    UVD(1:NGI) = 0.0
    VUD(1:NGI) = 0.0
    VVD(1:NGI) = 0.0
    GAMMAU(1:NGI) = 0.0
    GAMMAV(1:NGI) = 0.0
    !
    IF(DISCRA.EQ.1) THEN
      ! balancing diffusion based on (x,y) space. 
      UPWIND=.TRUE.
    ENDIF
    IF(DISCRA.EQ.2) THEN
      ! Laxwendrof balancing diffusion. 
      UPWIND=.FALSE.
      RGAMMA=0.5*DT
    ENDIF
    IF(DISCRA.EQ.3) THEN
      ! (x,y,t) -balancing diffusion. 
      UPWIND=.TRUE.
      XTPETV=.TRUE.
    ENDIF
    IF(DISCRA.EQ.4) THEN
      ! No balancing diffusion.
    ENDIF
    IF((DISCRA.GE.5).OR.(DISCRA.LE.7)) THEN
      ! nonlinear streamline and cross stream diffusion
      UPWIND=.TRUE.
      NONLIN=.TRUE.
    ENDIF
  
    !
    do  GI=1,NGI! Was loop 29
      GAMMA(GI)=RGAMMA
    end do ! Was loop 29
    ! END OF OPTIONS *******************************
    !
    !        ewrite(3,*)'HERE533'
    !
    do  ELE=1,TOTELE! Was loop 340
      !          ewrite(3,*)'ELE=',ELE
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
          !
          MUGIXX=0.
          MUGIXY=0.
          MUGIYY=0.
          !
          XABDGI(GI)=0.
          YABDGI(GI)=0.
          !
          do  L=1,NLOC! Was loop 79
            IGLS=SONDGL((ELE-1)*NLOC+L)
            IGLV=VONDGL((ELE-1)*NLOC+L)
            IGLX=XONDGL((ELE-1)*NLOC+L)
            AGI=AGI + NLX(L,GI)*X(IGLX) 
            BGI=BGI + NLX(L,GI)*Y(IGLX) 
            CGI=CGI + NLY(L,GI)*X(IGLX) 
            DGI=DGI + NLY(L,GI)*Y(IGLX)
            !        ewrite(3,*)'IGLX,X(IGLX),Y(IGLX):',IGLX,X(IGLX),Y(IGLX)
            !
            IF(DISOPN.EQ.0) THEN
                UD(GI)=UD(GI) + N(L,GI)*(NU(IGLV)-UG(IGLV))
                VD(GI)=VD(GI) + N(L,GI)*(NV(IGLV)-VG(IGLV))
            ENDIF
            DENGI(GI)=DENGI(GI)+ N(L,GI)*DENPT(IGLS)
            !
            MUGIXX=MUGIXX+ N(L,GI)*MUPTXX(IGLS)
            MUGIXY=MUGIXY+ N(L,GI)*MUPTXY(IGLS)
            MUGIYY=MUGIYY+ N(L,GI)*MUPTYY(IGLS)
            !
            !
            XABDGI(GI)=XABDGI(GI) + N(L,GI)*XABSOR(IGLS)
            YABDGI(GI)=YABDGI(GI) + N(L,GI)*YABSOR(IGLS)
            !
          end do ! Was loop 79
          !
          DETJ=AGI*DGI-BGI*CGI
          !        ewrite(3,*)'DETJ,AGI,DGI,BGI,CGI,GIVOLF:',
          !     &                DETJ,AGI,DGI,BGI,CGI,GIVOLF
          DETWEI(GI)=abs(DETJ)*WEIGHT(GI)
          DDETWE(GI)=DENGI(GI)*DETWEI(GI)
          !         MDETWE(GI)=MUGI(GI)*DETWEI(GI)
          MXXDTW(GI)=MUGIXX*DETWEI(GI)
          MXYDTW(GI)=MUGIXY*DETWEI(GI)
          MYYDTW(GI)=MUGIYY*DETWEI(GI)
          !
          ! Introduce grid vels in non-linear terms.
          !
          R=0.
          do  L=1,NLOC! Was loop 373
            NX(L,GI)=(DGI*NLX(L,GI)-BGI*NLY(L,GI))/DETJ
            NY(L,GI)=(-CGI*NLX(L,GI)+AGI*NLY(L,GI))/DETJ
            IGLV=VONDGL((ELE-1)*NLOC+L)
            R=R+NX(L,GI)*NU(IGLV)+NY(L,GI)*NV(IGLV)
          end do ! Was loop 373
          !
          !
          IF(DISOPN.EQ.0) THEN
            BECTWE(GI)=R*BETA*DETWEI(GI)
          ELSE
            BECTWE(GI)=0.0
          ENDIF
          !
          ! Work out GAMMA at each Gauss pt(OR centre of element)
          IF(UPWIND) THEN 
            !
            HXGI=ABS( DGI*UD(GI)-BGI*VD(GI))/DETJ
            HYGI=ABS(-CGI*UD(GI)+AGI*VD(GI))/DETJ
            !
            IF(NONLIN) THEN
                !
                !C SIMILAR TO HXGI,HYGI (HH2 streamline direction)
                HH2=2./MAX(HXGI,HYGI,1.E-7)
                !
                ! CALCULATE LIMIU,LIMIT used for limiting
                ! limit for U...
                CALL LIMGRA(.FALSE.,NGI,NLOC,NONODS,TOTELE, &
                    &             ELE,GI,&
                    &             NDGLNO, UD,VD,VD, U, HH2, &
                    &             N,NX,NY,NY, (DISCRA.EQ.7),&
                    ! 2-d Inverse of Jacobian...
                    &       DETJ, AGI,BGI,&
                    &             CGI,DGI,&
                    ! 3-d inverse of Jacobian...
                    &             AGI,AGI,AGI,&
                    &             AGI,AGI,AGI,&
                    &             AGI,AGI,AGI,&
                    ! THE LIMITING VARIABLEs
                    &             UUD(GI),UVD(GI),UVD(GI),LOVERQ )
                GAMMAU(GI)=RWIND*LOVERQ
  
                ! limit for V...
                CALL LIMGRA(.FALSE.,NGI,NLOC,NONODS,TOTELE, &
                    &             ELE,GI,&
                    &             NDGLNO, UD,VD,VD, V, HH2, &
                    &             N,NX,NY,NY, (DISCRA.EQ.7),&
                    ! 2-d Inverse of Jacobian...
                    &       DETJ, AGI,BGI,&
                    &             CGI,DGI,&
                    ! 3-d inverse of Jacobian...
                    &             AGI,AGI,AGI,&
                    &             AGI,AGI,AGI,&
                    &             AGI,AGI,AGI,&
                    ! THE LIMITING VARIABLEs
                    &             VUD(GI),VVD(GI),VVD(GI),LOVERQ )
                GAMMAV(GI)=RWIND*LOVERQ
                ! ENDOF IF(NONLIN) THEN...
            ENDIF
            !
            ! HX,HY are the characteristic length scales in x,y directions.
            ! If  XTPETV  then (x,y,t) Petrov Galerkin Method  else (x,y).
            IF(XTPETV) THEN
                !                     HOVERQ=MIN(HXGI/RUD, HYGI/RVD, DT)
                HOVERQ=MIN(2./MAX(HXGI,HYGI,1.E-7), DT)
            ELSE
                !                     HOVERQ=MIN(HXGI/RUD, HYGI/RVD)
                HOVERQ=2./MAX(HXGI,HYGI,1.E-7)
            ENDIF
            GAMMA(GI)=RWIND*HOVERQ
            IF(DISCRA.EQ.6) GAMMA(GI)=0.0
          ENDIF
      end do ! Was loop 331
      !
      do  ILOC=1,NLOC! Was loop 350
          GLOBI =NDGLNO((ELE-1)*NLOC+ILOC)
          GLOBIS=SONDGL((ELE-1)*NLOC+ILOC)
          do  JLOC=1,NLOC! Was loop 360
            GLOBJ =NDGLNO((ELE-1)*NLOC+JLOC)
            GLOBJS=SONDGL((ELE-1)*NLOC+JLOC)
            !
            HL1=0.
            HL2=0.
            VLND=0.
            VLM=0.
            VLMD=0.
            !
            VLKXX=0.    
            VLKXY=0.     
            VLKYX=0.    
            VLKYY=0.
            !
            VXVXY=0.
            VYVYX=0.
            !
            VABSX=0.
            VABSY=0.
            !
            do  GI=1,NGI! Was loop 458
                R1=UD(GI)*NX(ILOC,GI)+VD(GI)*NY(ILOC,GI) 
                R2=(UD(GI)*NX(JLOC,GI)&
                    &                      +VD(GI)*NY(JLOC,GI) )*DETWEI(GI)
                !
                R1U=UUD(GI)*NX(ILOC,GI)+UVD(GI)*NY(ILOC,GI) 
                R2U=(UUD(GI)*NX(JLOC,GI)&
                    &                      +UVD(GI)*NY(JLOC,GI) )*DETWEI(GI)
                !
                R1V=VUD(GI)*NX(ILOC,GI)+VVD(GI)*NY(ILOC,GI) 
                R2V=(VUD(GI)*NX(JLOC,GI)&
                    &                      +VVD(GI)*NY(JLOC,GI) )*DETWEI(GI)
                !
                VLND=VLND+DENGI(GI)*N(ILOC,GI)*(R2   &
                    &                     + BECTWE(GI)*N(JLOC,GI)  )
                RNN=N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
                !  N.B.   BECTWE(GI) = BETA*CTY(GI)*DETWEI(GI) .   
                !
                HL1=HL1 + GAMMA(GI) *DENGI(GI)*R1 *R2&
                    &                 + GAMMAU(GI)*DENGI(GI)*R1U*R2U
                ! 
                HL2=HL2 + GAMMA(GI) *DENGI(GI)*R1 *R2&
                    &                 + GAMMAV(GI)*DENGI(GI)*R1V*R2V
                ! 
                VABSX=VABSX + XABDGI(GI)*RNN
                VABSY=VABSY + YABDGI(GI)*RNN
                !
                VLM=VLM + RNN
                VLMD=VLMD+RNN*DENGI(GI)
                !
                VLKXX=VLKXX+NX(ILOC,GI)*NX(JLOC,GI)*MXXDTW(GI)
                VLKXY=VLKXY+NX(ILOC,GI)*NY(JLOC,GI)*MXYDTW(GI)
                VLKYX=VLKYX+NY(ILOC,GI)*NX(JLOC,GI)*MXYDTW(GI)
                VLKYY=VLKYY+NY(ILOC,GI)*NY(JLOC,GI)*MYYDTW(GI)
                ! the cty parts...
                !               VXVXX=VXVXX+NX(ILOC,GI)*NX(JLOC,GI)*MXXDTW(GI)
                VXVXY=VXVXY+NX(ILOC,GI)*NY(JLOC,GI)*MXXDTW(GI)
                !
                VYVYX=VYVYX+NY(ILOC,GI)*NX(JLOC,GI)*MYYDTW(GI)
                !               VYVYY=VYVYY+NY(ILOC,GI)*NY(JLOC,GI)*MYYDTW(GI)
                !
            end do ! Was loop 458
            VXVXX=VLKXX
            VYVYY=VLKXX
            !
            ! Make the lumped absorption consistent with a lumped source. 
            VABSX=VABSX*(1.-ABLUMP) &
                  &           + VLM*XABSOR(GLOBIS)*ABLUMP
            VABSY=VABSY*(1.-ABLUMP) &
                  &           + VLM*YABSOR(GLOBIS)*ABLUMP
            !
            ! lump the mass matrix.(LUMP=1. if we are to LUMP the mass matrix)
            ICENT=CENTRM(GLOBI)
            BIGM(ICENT+IBL11)=BIGM(ICENT+IBL11)+VLMD*RLUMP&
                  &                        + CVMVEC(GLOBIS)*VLM &
                  &                        +DT*THETA*VABSX*ABLUMP
            BIGM(ICENT+IBL22)=BIGM(ICENT+IBL22)+VLMD*RLUMP&
                  &                        + CVMVEC(GLOBIS)*VLM &
                  &                        +DT*THETA*VABSY*ABLUMP
            ! At the moment ABLUMP only works when VABSX=VABSY because we have only one 
            ! lumped mass matrix. 
            ML(GLOBI)=ML(GLOBI) + VLMD + CVMVEC(GLOBIS)*VLM &
                  &                + DT*THETA*VABSX*ABLUMP
            !           ML(GLOBI)=ML(GLOBI) + VLMD + CVMVEC(GLOBIS)*VLM
  
            !
            !         R=VLND+HL
            R=VLND
            VECX(GLOBI)=VECX(GLOBI) &
                  &                -(R+HL1+VABSX*(1.-ABLUMP)+2.*VLKXX+VLKYY)*U(GLOBJ)&
                  &                -VLKYX*V(GLOBJ)&
                  &      +    TWOTHI*(VXVXX*U(GLOBJ)+VXVXY*V(GLOBJ))&
                  &      + VLM*sourcx(globjS)  *(1.-RSLUMP)
            !
            VECY(GLOBI)=VECY(GLOBI)&
                  &                -(R+HL2+VABSY*(1.-ABLUMP)+VLKXX+2.*VLKYY)*V(GLOBJ)&
                  &                -VLKXY*U(GLOBJ)&
                  &      +    TWOTHI*(VYVYX*U(GLOBJ)+VYVYY*V(GLOBJ))&
                  &      + VLM*sourcy(globjS) *(1.-RSLUMP)
            !
            ! Lumped source & absorption. 
            VECX(GLOBI)=VECX(GLOBI)-VABSX*ABLUMP*U(GLOBI)&
                  &                            +VLM*sourcx(GLOBIS) *RSLUMP
            VECY(GLOBI)=VECY(GLOBI)-VABSY*ABLUMP*V(GLOBI)&
                  &                            +VLM*sourcy(GLOBIS) *RSLUMP
            !           
            !
            !
            ! we will now place contributions into the matrices BIGM and BIGT
            ! these are symm matrices so we only store the upper diagonal
            ! parts of these matrices.
            MC1=(VLND*RSYM  )*THETA*DT&
                  &        + VLMD*(1.-RLUMP) 
  
            ! iF SYM make BIGM matrix symmetric. 
            !      
            ! Find count. ***************************************
            LOWER=FINDRM(GLOBI) 
            UPPER=FINDRM(GLOBI+1)-1
  7000       CONTINUE
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
  9000       CONTINUE
            ! Find count. ***************************************
            !
            !             COUNT=POSBM(ELE,(ILOC-1)*8+JLOC)
            !             DO 564 COUNT=FINDRM(GLOBI),FINDRM(GLOBI+1)-1
            !               IF(COLM(COUNT).EQ.GLOBJ) 
            !     &            BIGM(COUNT)=BIGM(COUNT)+MC1 + VLMD*(LUMP-1)
            !564          CONTINUE
            !
            BIGM(COUNT+IBL11)=BIGM(COUNT+IBL11)&
                  &            +MC1 +DT*THETA*(HL1+VABSX*(1.-ABLUMP)+2.*VLKXX+VLKYY&
                  &            - TWOTHI*VXVXX)
            BIGM(COUNT+IBL22)=BIGM(COUNT+IBL22)&
                  &            +MC1 +DT*THETA*(HL2+VABSY*(1.-ABLUMP)+VLKXX+2.*VLKYY&
                  &            - TWOTHI*VYVYY)
            !
            BIGM(COUNT+IBL21)=BIGM(COUNT+IBL21)&
                  &                  +DT*THETA*(VLKXY - TWOTHI*VYVYX)
            !
            IF(.NOT.BLKSYM) THEN
                BIGM(COUNT+IBL12)=BIGM(COUNT+IBL12)&
                    &                  +DT*THETA*(VLKYX - TWOTHI*VXVXY)
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
  7001          CONTINUE
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
  9001          CONTINUE
                ! Find count. ***************************************
                do  GI=1,NGI! Was loop 123
                  C1T(POS)=C1T(POS)+M(JLOC,GI)*NX(ILOC,GI)*DETWEI(GI)
                  C2T(POS)=C2T(POS)+M(JLOC,GI)*NY(ILOC,GI)*DETWEI(GI)
                end do ! Was loop 123
            end do ! Was loop 123
            !
          ENDIF
      end do ! Was loop 350
    end do ! Was loop 340
    !
    !        ikeep=0
    !        rdis=1.e+20
    !        do i=1,nonods
    !          if( (x(i)-31.82)**2+(y(i)-31.82)**2 .lt. rdis) then
    !             rdis=(x(i)-31.82)**2+(y(i)-31.82)**2
    !             ikeep=i
    !          endif
    !        end do
    !        ewrite(3,*) 'sqrt(rdis),ikeep,x(ikeep),y(ikeep):',
    !     &           sqrt(rdis),ikeep,x(ikeep),y(ikeep)
    !        ewrite(3,*) 'sourcx(ikeep),sourcy(ikeep):',
    !     &           sourcx(ikeep),sourcy(ikeep)
    !        ewrite(3,*) 'vecx(ikeep),vecy(ikeep),oth:',
    !     &           vecx(ikeep),vecy(ikeep),(ml(ikeep)/2.66)*sourcy(ikeep)
    !        ewrite(3,*) 'MUPTXX(ikeep),xabsor(ikeep):',
    !     &           MUPTXX(ikeep),xabsor(ikeep)
    !        ewrite(3,*) 'MUPTxy(ikeep),MUPTyy(ikeep),yabsor(ikeep):',
    !     &           MUPTxy(ikeep),MUPTyy(ikeep),yabsor(ikeep)
    !        ewrite(3,*) 'ml(ikeep)=',ml(ikeep)
    !         stop 3636
    !         ewrite(3,*) 'in diff2d cvmvec',cvmvec
    !        CALL PMINMX(MUPTXX,NONODS,'******MUPTXX  ')
    !        CALL PMINMX(MUPTXY,NONODS,'******MUPTXY  ')
    !        CALL PMINMX(MUPTYY,NONODS,'******MUPTYY  ')
    !        CALL PMINMX(XABSOR,NONODS,'******XABSOR  ')
    !        CALL PMINMX(YABSOR,NONODS,'******YABSOR  ')
    !        CALL PMINMX(SOURCX,NONODS,'******SOURCX  ')
    !        CALL PMINMX(SOURCY,NONODS,'******SOURCY  ')
    !        CALL PMINMX(DENPT, NONODS,'******DENPT   ')
    !        ewrite(3,*) 'GETC12=',GETC12
    !        ewrite(3,*) 'R2NORM(C1T,NCT,0):',R2NORM(C1T,NCT,0)
    !        ewrite(3,*) 'R2NORM(C2T,NCT,0):',R2NORM(C2T,NCT,0)
    ewrite(1,*)'EXITING DIFF2D'
    RETURN
  END SUBROUTINE DIFF2D
  !
end module momentum_matrix_cg_2d

!
!
