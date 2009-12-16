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

! A collection of functions and routines moved out of hart3d, to enable hart3d
! to be converted into a module

#include "fdebug.h"

module hart3d_allsorts

  use allsorts
  use FLDebug
  use position_in_matrix
  use shape_module
  use tr2d_module
  use FETools
  use sml
  use detnlxr_module
  implicit none
  
  private
  
  public :: abmatrixmul, abmatrixmulti, anonvdlim, bothanonvdlim, ctdiff, &
    & finptsstore, gaussiloc, getgxyz, getstoreelewei, getsubscale, &
    & imatpts, matdmatinv, matpts, onvdlim, onvdlim_c, philnodele, &
    & gradbound, FUNADP, DEDADP, invxxx, phacmc
  
contains

  REAL FUNCTION DEDADP(NPROPT,VOLSOL,ALFST2,SPRESS)
    !     This function returns the partial 
    !     derivative of silid volume fraction VOLSOL
    !     w.r.t solid pressure 
    !     and the option NPROPT. 
    !     ALFST2= max volume fraction of solids. 
    REAL, parameter::THIRD=1.0/3.0

    INTEGER, intent(in)::NPROPT
    REAL, intent(in)::VOLSOL,ALFST2,SPRESS
    REAL G0 
    !     
    IF(ABS(NPROPT).EQ.2) THEN
       G0=(1.-(VOLSOL/ALFST2))
       !     
       DEDADP=(ALFST2*G0**4/(3.*VOLSOL +ALFST2*G0))/SPRESS
    ELSE
       G0=(1.-(VOLSOL/ALFST2)**THIRD)
       !     
       DEDADP=(G0**2/(THIRD*(VOLSOL/ALFST2)**THIRD +G0))/SPRESS
    ENDIF
    !     
  END FUNCTION DEDADP

  REAL FUNCTION FUNADP(NPROPT,RADISO,&
       &     VOLSOL,ALFST2,SPRESS,KETPRE)
    !     This function returns the pressure function 
    !     given the solid volume fraction VOLSOL 
    !     and the option NPROPT. 
    !     ALFST2= max volume fraction of solids. 
    REAL THIRD
    PARAMETER(THIRD=1./3.)
    INTEGER NPROPT,RADISO
    REAL VOLSOL,ALFST2,SPRESS,KETPRE
    !     DISFIN=inverse of the radial distribution function
    !     
    GO TO (10,20,30,40,50,60) ABS(NPROPT)
    !     For NPROPT=1 ++++++++++++++++++++++++++++++++++++++++
10  CONTINUE
    FUNADP=(1.-(VOLSOL/ALFST2)**THIRD)/SPRESS
    GO TO 1000
    !     For NPROPT=2 ++++++++++++++++++++++++++++++++++++++++
20  CONTINUE
    FUNADP=(1.-VOLSOL/ALFST2)**3/SPRESS
    GO TO 1000
    !     For NPROPT=3 ++++++++++++++++++++++++++++++++++++++++
30  CONTINUE
    FUNADP=MAX(0.1,VOLSOL)**3 *EXP(-ALFST2*VOLSOL)/SPRESS 
    GO TO 1000
    !     For NPROPT=4 ++++++++++++++++++++++++++++++++++++++++
40  CONTINUE
    FUNADP=MIN(   MAX(0.1,VOLSOL)/SPRESS,&
         &     (1.-VOLSOL/0.6)**3/ALFST2  )
    GO TO 1000
    !     For NPROPT=5 ++++++++++++++++++++++++++++++++++++++++
    !     KINETIC ENERGY METHOD (KETPRE=dengas*2.*(1.+RECOE)*KET)...
50  CONTINUE
    !     FUNADP=MAX(0.05,VOLSOL) * DISFIN(VOLSOL,ALFST2)
    !==> changed here as disfun is used for multiphase
    !      FUNADP=MAX(0.01,VOLSOL) *(1./DISFUN(VOLSOL,ALFST2,RADISO))
    !     &     / KETPRE  
    FUNADP=MAX(0.01,VOLSOL) * 1./&
         &        MAX(1.E-5,DISFUNSUB(VOLSOL,ALFST2,RADISO))&
         &     / KETPRE  
    GO TO 1000
    !     For NPROPT=6 (Fosco Gibilera model) +++++++++++++++++
60  CONTINUE
    FUNADP=(1.-VOLSOL/ALFST2)**3/SPRESS
    !     FUNADP=1./DISFG1(VOLSOL,ALFST2,600.) 
    !     ewrite(3,*) 'VOLSOL,ALFST2:',VOLSOL,ALFST2
    !     FUNADP=DISFG1(max(VOLSOL,ALFST2-0.1),ALFST2,-600.) 
    GO TO 1000
    !     
    !     
1000 CONTINUE
  END FUNCTION FUNADP

         SUBROUTINE GETSUBSCALE(TSUB,T,DINVSOUSUBTSTORE,DMATINVCSTORE,&
     &                          NONODS,TOTELE,NLOC,NDGLNO,NSUBTLOC)
! This sub calculates the inner element model TSUB from T
         INTEGER NONODS,TOTELE,NLOC,NSUBTLOC
         INTEGER NDGLNO(TOTELE*NLOC)
         REAL TSUB(TOTELE*NSUBTLOC),T(NONODS)
         REAL DMATINVCSTORE(TOTELE,NSUBTLOC,NLOC)
         REAL DINVSOUSUBTSTORE(TOTELE,NSUBTLOC)
! Local variables...
         INTEGER ELE,ILOC,JLOC,GLOBJ
!         return
      do ELE=1,TOTELE! Was loop 
      do ILOC=1,NSUBTLOC! Was loop 
             TSUB((ELE-1)*NSUBTLOC+ILOC)=DINVSOUSUBTSTORE(ELE,ILOC)
           END DO
           
      do JLOC=1,NLOC! Was loop 
             GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
      do ILOC=1,NSUBTLOC! Was loop 
               TSUB((ELE-1)*NSUBTLOC+ILOC)=TSUB((ELE-1)*NSUBTLOC+ILOC)&
     &                -DMATINVCSTORE(ELE,ILOC,JLOC)*T(GLOBJ)
! DMATINVCSTORE is not a square matrix if NSUBTLOC.ne.NLOC...
             END DO
           END DO
         END DO
         RETURN
         
  end subroutine getsubscale
         
!     
!    
!     
!    
         SUBROUTINE MATDMATINV(DMAT,DMATINV,NLOC)
! calculate DMATINV...
         INTEGER NLOC
         REAL DMAT(NLOC,NLOC),DMATINV(NLOC,NLOC)
! Local variables...
         REAL MAT(NLOC,NLOC),MAT2(NLOC,NLOC),X(NLOC),B(NLOC)
         
         DMATINV=DMAT
         CALL MATINV(DMATINV,NLOC,NLOC,MAT,MAT2,X,B)

         RETURN

  end subroutine matdmatinv

!
!     
!    
!
          SUBROUTINE ABMATRIXMUL(AB, A,NONODS1,NONODS2, &
     &                               B,NONODS3,NONODS4)
! Perform matrix matrix AB=A*B
          INTEGER NONODS1,NONODS2,NONODS3,NONODS4
          REAL AB(NONODS1,NONODS4)
          REAL A(NONODS1,NONODS2),B(NONODS3,NONODS4)
! Logical I,J
          INTEGER I,J,II
!          IF(NONODS2.NE.NONODS3) STOP 8329
      do I=1,NONODS1! Was loop 
      do J=1,NONODS4! Was loop 
              AB(I,J)=0.0
      do II=1,NONODS2! Was loop 
                AB(I,J)=AB(I,J) + A(I,II)*B(II,J)
              END DO
            END DO
          END DO
          RETURN

  end subroutine abmatrixmul

!
!      
!    
!
          SUBROUTINE ABMATRIXMULTI(AB, A,B,NONODS)
! Perform matrix matrix AB=A*B
          INTEGER NONODS
          REAL AB(NONODS,NONODS)
          REAL A(NONODS,NONODS),B(NONODS,NONODS)
! Logical I,J
          INTEGER I,J,II
      do I=1,NONODS! Was loop 
      do J=1,NONODS! Was loop 
              AB(I,J)=0.0
      do II=1,NONODS! Was loop 
                AB(I,J)=AB(I,J) + A(I,II)*B(II,J)
              END DO
            END DO
          END DO
          RETURN

  end subroutine abmatrixmulti
          
      SUBROUTINE CTDIFF(X,Y,Z,&
     &     SCMC,PML,&
     &     SFINCMC,SCOLCMC,NSCMC,FREDOP,TOTELE, &
     &     N,NLX,NLY,NLZ, &
     &     M,MLX,MLY,MLZ, WEIGHT,NGI,NLOC,MLOC, &
     &     PNDGLN,XONDGL,XNONOD,&
     &     DETWEI,&
     &     UDL,VDL,WDL,&
     &     GAMMA,&
     &     MX,MY,MZ,ISPHERE     )
!     This sub calculates the stabalisation diffusion matrix SCMC. 
!     Eventuall KCMC=SCMC^T PML^-1 SCMC
!     PNDGLN=element pter for unknown pressures. 
!     XONDGL=element pter for coordinates(x,y,z)
!     -------------------------------------------------------------------------------------------
!     If we are solving for another variable like temperature 
!     or chemical species then NU,NV,NW will be the velocities 
!     and U=TEMPERATURE OR CHEMICAL SPECIES. 
!     -------------------------------------------------------------------------------------------
!     NB SOURCE, DENPT, viscosity are SNONOD in length and have ele-pt SONDGL.
! ISPHERE indicates the order of the superparamtertic FEM mapping
      INTEGER XNONOD
      INTEGER MLOC,NLOC,NGI
      INTEGER TOTELE,FREDOP,NSCMC,ISPHERE
!     
      REAL SCMC(NSCMC),PML(FREDOP)
      REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
      INTEGER PNDGLN(TOTELE*MLOC)
      INTEGER XONDGL(TOTELE*NLOC)
      INTEGER SFINCMC(FREDOP+1),SCOLCMC(NSCMC)
!     
      REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
      REAL M(MLOC,NGI),MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
      REAL WEIGHT(NGI)
!     
      REAL DETWEI(NGI)
      REAL UDL(NLOC*NLOC,NGI),VDL(NLOC*NLOC,NGI)
      REAL WDL(NLOC*NLOC,NGI)
      REAL GAMMA(NLOC*NLOC,NGI)
      REAL MX(MLOC,NGI),MY(MLOC,NGI),MZ(MLOC,NGI)
!     HX,HY-characteristic length scales in x,y directions.
!     Local variables...
      REAL DIFLIN
      REAL RN,R1,R2
      REAL XC,YC,ZC
!     
      INTEGER ELE,ILOC,JLOC,L,L1,L2,IGLX1,IGLX2,IGLX,ID,NID
      INTEGER GI,GLOBI,GLOBJ
      INTEGER POSMAT
!     
      REAL HXGI,HYGI,HZGI,DETJ
      REAL MATRIX,PMASS
      REAL HOVERQ
      REAL AGI,BGI,CGI,DGI
      REAL EGI,FGI,GGI,HGI,KGI
      REAL A11,A12,A13, A21,A22,A23, A31,A32,A33
      REAL RWIND
      ewrite(1, *) "SUBROUTINE CTDIFF()"
!     
      RWIND =1./REAL(6)
      NID=NLOC*NLOC
!     
      SCMC(1:NSCMC) = 0.0
      PML(1:FREDOP) = 0.0
!     
!     This subroutine forms a contabution to the Right Hand Side
!     of Poissons pressure equation, as well as  F1 & F2.
!     
      do  ELE=1,TOTELE! Was loop 340
!     
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
!     
            end do ! Was loop 79
!     
            DETJ= AGI*(EGI*KGI-FGI*HGI)&
                 -BGI*(DGI*KGI-FGI*GGI)&
                 +CGI*(DGI*HGI-EGI*GGI)
            DETWEI(GI)=abs(DETJ)*WEIGHT(GI)
!     
!     For coefficient in the inverse mat of the jacobian. 
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
!     
            do  L=1,MLOC! Was loop 373
               MX(L,GI)= A11*MLX(L,GI)+A12*MLY(L,GI)+A13*MLZ(L,GI)
               MY(L,GI)= A21*MLX(L,GI)+A22*MLY(L,GI)+A23*MLZ(L,GI)
               MZ(L,GI)= A31*MLX(L,GI)+A32*MLY(L,GI)+A33*MLZ(L,GI)
            end do ! Was loop 373
!     

!     C The first is the old filter term, the second the new one MDP getting
!     c different results and stabiltiy for tidal applications ????
            if(.false.) then
               RWIND =1./REAL(NLOC)
               NID=NLOC
!     **********calculate normalised velocitys across element...          
               XC=0. 
               YC=0. 
               ZC=0.

               do L=1,NLOC! Was loop 
                  IGLX=XONDGL((ELE-1)*NLOC+L)
                  XC=XC+X(IGLX)/REAL(NLOC)
                  YC=YC+Y(IGLX)/REAL(NLOC)
                  ZC=ZC+Z(IGLX)/REAL(NLOC)
               END DO

               do L=1,NLOC! Was loop 
                  IGLX=XONDGL((ELE-1)*NLOC+L) 
                  ID=L
                  UDL(ID,GI)=X(IGLX)-XC
                  VDL(ID,GI)=Y(IGLX)-YC
                  WDL(ID,GI)=Z(IGLX)-ZC
!     
!     Normalise 
                  RN=SQRT(UDL(ID,GI)**2+VDL(ID,GI)**2+WDL(ID,GI)**2)
                  UDL(ID,GI)=UDL(ID,GI)/RN
                  VDL(ID,GI)=VDL(ID,GI)/RN
                  WDL(ID,GI)=WDL(ID,GI)/RN
!     
!     put the length scale in
                  HXGI=ABS(A11*UDL(ID,GI)+A12*VDL(ID,GI)+A13*WDL(ID,GI))
                  HYGI=ABS(A21*UDL(ID,GI)+A22*VDL(ID,GI)+A23*WDL(ID,GI))
                  HZGI=ABS(A31*UDL(ID,GI)+A32*VDL(ID,GI)+A33*WDL(ID,GI))
!     HX,HY are the characteristic length scales in x,y directions. 
                  HOVERQ=2./MAX(HXGI,HYGI,HZGI,1.E-7) 
!     HOVERQ=MIN(2./MAX(HXGI,HYGI,HZGI,1.E-7), DT)
                  GAMMA(ID,GI)=RWIND*HOVERQ 
               END DO 
!     **********calculate normalised velocitys across element...    
            else
               RWIND =1./REAL(6)
               NID=NLOC*NLOC
!     **********calculate normalised velocitys across element...          
!     XC=0. 
!     YC=0. 
!     ZC=0.
!     DO L=1,NLOC
!     IGLX=XONDGL((ELE-1)*NLOC+L)
!     XC=XC+X(IGLX)/REAL(NLOC)
!     YC=YC+Y(IGLX)/REAL(NLOC)
!     ZC=ZC+Z(IGLX)/REAL(NLOC)
!     END DO
               ID=0
               do L1=1,NLOC! Was loop 
                  IGLX1=XONDGL((ELE-1)*NLOC+L1) 
                  do L2=1,NLOC! Was loop 
                     IGLX2=XONDGL((ELE-1)*NLOC+L2) 
                     ID=ID+1
                     if(l1.eq.l2) then
                        UDL(ID,GI)=0.0
                        VDL(ID,GI)=0.0
                        WDL(ID,GI)=0.0
                        GAMMA(ID,GI)=0.0
                     else
                        UDL(ID,GI)=X(IGLX1)-X(IGLX2)
                        VDL(ID,GI)=Y(IGLX1)-Y(IGLX2)
                        WDL(ID,GI)=Z(IGLX1)-Z(IGLX2)
!     
!     Normalise 
                        RN=SQRT(UDL(ID,GI)**2+VDL(ID,GI)**2+WDL(ID,GI)**2)
                        UDL(ID,GI)=UDL(ID,GI)/RN
                        VDL(ID,GI)=VDL(ID,GI)/RN
                        WDL(ID,GI)=WDL(ID,GI)/RN
!     
!     put the length scale in
                        HXGI=ABS(A11*UDL(ID,GI)+A12*VDL(ID,GI)+A13*WDL(ID,GI))
                        HYGI=ABS(A21*UDL(ID,GI)+A22*VDL(ID,GI)+A23*WDL(ID,GI))
                        HZGI=ABS(A31*UDL(ID,GI)+A32*VDL(ID,GI)+A33*WDL(ID,GI))
!     HX,HY are the characteristic length scales in x,y directions. 
                        HOVERQ=RN
                        GAMMA(ID,GI)=RWIND*HOVERQ 
                     endif
                  END DO 
               END DO 
!     **********calculate normalised velocitys across element...            
            endif   
!     
         end do ! Was loop 331
!     
         do  ILOC=1,MLOC! Was loop 350
            GLOBI=PNDGLN((ELE-1)*MLOC+ILOC)
            do  JLOC=1,MLOC! Was loop 360
               GLOBJ=PNDGLN((ELE-1)*MLOC+JLOC)
!     
               MATRIX=0.
               PMASS =0. 
!     
               do  GI=1,NGI! Was loop 458
                  do  ID=1,NID! Was loop 454
                     R1= UDL(ID,GI)*MX(ILOC,GI)&
                        +VDL(ID,GI)*MY(ILOC,GI)&
                        +WDL(ID,GI)*MZ(ILOC,GI) 
                     R2= UDL(ID,GI)*MX(JLOC,GI)&
                        +VDL(ID,GI)*MY(JLOC,GI) &
                        +WDL(ID,GI)*MZ(JLOC,GI)  
!     
!     Matrix contributions...
                     DIFLIN= - GAMMA(ID,GI)*R1*R2*DETWEI(GI)
                     MATRIX=MATRIX + DIFLIN
!     USE THE COMPONENT OF DIFLIN THE X,Y & Z-DIRECTIONS 
!     RESPECTIVELY FOR C1T,C2T,C3T.
                  end do ! Was loop 454
                  PMASS =PMASS +M(ILOC,GI)*M(JLOC,GI)*DETWEI(GI) 
               end do ! Was loop 458
!     
!     we will now place contributions into the matrices SCMC and BIGT
!     these are can be symm matrices.
!     Find POSMAT. ***************************************
!     
               CALL POSINMAT(POSMAT,GLOBI,GLOBJ,&
                             FREDOP,SFINCMC,SCOLCMC,NSCMC)
!     Find POSMAT. ***************************************
!     
               PML(GLOBI)=PML(GLOBI)+PMASS
               SCMC(POSMAT)=SCMC(POSMAT)+MATRIX 
!     
!     
            end do ! Was loop 360
         end do ! Was loop 350
      end do ! Was loop 340

      RETURN 

  end subroutine ctdiff

!     
!     SUBROUTINE IFINPTS HAS MOVED TO LEGACY_advection_diffusion_CV.F90 (CRGW)
!
!     
!     
!     

      SUBROUTINE IMATPTS(MATPSI,MATPSIOLD,NOD,&
     &     PSI,PSIOLD,NONODS,&
     &     PSI2,PSIOLD2,&
     &     NLOC,TOTELE,NDGLNO,&
     &     FINDRM,COLM,NCOLM,&
     &     X1,Y1,Z1,&
     &     X2,Y2,Z2,&
     &     NORMX1,NORMY1,NORMZ1,&
     &     X,Y,Z,&
!     work space...
     &     FINDELE,COLELE,NCOLEL) 
!     This sub calculates the value of PSI that would be at the 
!     other side of the stencil if we had a linear variation and within 
!     a single element.
!     It does this for the interface tracking method...   
      REAL INFINY,RELAX
      INTEGER OINTER
!     OINTER is the limiting option (as it increases the more 
!     compressive the method becomes).
!     OINTER=0  - locally bounded soln
!     OINTER=1  - make solution bounded beteen [0,1]
!     OINTER=2  - allow the soln to violate local boundedness, but still [0,1]
!     OINTER=3  - allow the soln to violate local boundedness, but still [0,1] more compressive.
!     RELAX \in [0.5,1.0] increases compression with smaller values.
      PARAMETER(INFINY=1.E+20,OINTER=1,RELAX=0.5)
      INTEGER NOD,NONODS,NLOC,TOTELE
      REAL MATPSI,MATPSIOLD,PSI(NONODS),PSIOLD(NONODS)
      REAL PSI2,PSIOLD2
      INTEGER NDGLNO(NLOC*TOTELE)
      INTEGER NCOLM,FINDRM(NONODS+1),COLM(NCOLM)
      REAL X1,Y1,Z1,X2,Y2,Z2,NORMX1,NORMY1,NORMZ1
      REAL X(NONODS),Y(NONODS),Z(NONODS)
      INTEGER NCOLEL
      INTEGER FINDELE(NONODS+1),COLELE(NCOLEL)
!     
!     Local variables...
      REAL XC,YC,ZC
      REAL LOCCORDS(4),LOCCORDSK(4)
      INTEGER LOCNODS(4),LOCNODSK(4)
      INTEGER COUNT,ELE,ILOC,KNOD,JNOD,III
      REAL MINCOR,MINCORK,SUM,DIST12
      REAL MINPSI,MINPSIOLD,MAXPSI,MAXPSIOLD,RN
      REAL T1X,   T1Y,   T1Z,   T2X,   T2Y,   T2Z
      REAL REFX2,REFY2,REFZ2,REFX,REFY,REFZ
      REAL VX,VY,VZ
!     
      XC=X1 - 0.01*(X2-X1)
      YC=Y1 - 0.01*(Y2-Y1)
      ZC=Z1 - 0.01*(Z2-Z1)
!     
      do  III=1,2! Was loop 20
!     III=2 is needed if (XC,YC,ZC) is outside the domain
         IF(III.EQ.2) THEN
!     The rotation matrix in 3-D is R=  
!     NX    NY    NZ
!     T1X   T1Y   T1Z
!     T2X   T2Y   T2Z
!     
            VX=X1-X2
            VY=Y1-Y2
            VZ=Z1-Z2
!     
            CALL XPROD(T2X,T2Y,T2Z, NORMX1,NORMY1,NORMZ1, VX,VY,VZ)
!     
            DIST12=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
            RN=SQRT(T2X**2+T2Y**2+T2Z**2)
            IF(RN.LT.(1.E-5)*DIST12) THEN
!     Simply have VX,VY,VZ going in the opposite direction...
               XC=X1 - VX*0.01
               YC=Y1 - VY*0.01
               ZC=Z1 - VZ*0.01
            ELSE
             T2X= T2X/RN
               T2Y= T2Y/RN
               T2Z= T2Z/RN
!     T1=Nx (-T2)
               CALL XPROD(T1X,T1Y,T1Z, NORMX1,NORMY1,NORMZ1, -T2X,-T2Y,-T2Z)
!     
               REFX2= NORMX1*VX + NORMY1*VY + NORMZ1*VZ
               REFY2= T1X   *VX + T1Y   *VY + T1Z   *VZ
               REFZ2= T2X   *VX + T2Y   *VY + T2Z   *VZ
!     Reflect...
               REFX2=-REFX2
!     MAP BACK USING R^T
               REFX = NORMX1*REFX2 + T1X*REFY2 + T2X*REFZ2
               REFY = NORMY1*REFX2 + T1Y*REFY2 + T2Y*REFZ2
               REFZ = NORMZ1*REFX2 + T1Z*REFY2 + T2Z*REFZ2
!     
!     (REFX,REFY,REFZ) is the reflected direction...
               XC=X1 + REFX*0.01
               YC=Y1 + REFY*0.01
               ZC=Z1 + REFZ*0.01
            ENDIF
         ENDIF
         MINCORK=-INFINY
!     
      do  COUNT=FINDELE(NOD),FINDELE(NOD+1)-1! Was loop 10
            ELE=COLELE(COUNT)
!     
            LOCNODS(1)=NDGLNO((ELE-1)*NLOC+1)
            LOCNODS(2)=NDGLNO((ELE-1)*NLOC+2)
            LOCNODS(3)=NDGLNO((ELE-1)*NLOC+3)
            LOCNODS(4)=NDGLNO((ELE-1)*NLOC+4)
!     
!     C Calculate the local coord but with 4th point replaced by INOD...
!     C Find local coords LOCCORDS of point INOD corresponding to these nodes LOCNODS...
            CALL TRILOCCORDS(XC,YC,ZC, &
     &           LOCCORDS(1),LOCCORDS(2),LOCCORDS(3),LOCCORDS(4),&
!     The 4 corners of the tet...
     &           X(LOCNODS(1)),Y(LOCNODS(1)),Z(LOCNODS(1)),&
     &           X(LOCNODS(2)),Y(LOCNODS(2)),Z(LOCNODS(2)),&
     &           X(LOCNODS(3)),Y(LOCNODS(3)),Z(LOCNODS(3)),&
     &           X(LOCNODS(4)),Y(LOCNODS(4)),Z(LOCNODS(4))  )
!     
            MINCOR=MIN(LOCCORDS(1),LOCCORDS(2),&
     &           LOCCORDS(3),LOCCORDS(4))
            IF(MINCOR.GT.MINCORK) THEN 
               MINCORK=MINCOR
      do ILOC=1,NLOC! Was loop 
                  LOCCORDSK(ILOC)=LOCCORDS(ILOC)
                  LOCNODSK(ILOC)=LOCNODS(ILOC)
               END DO
            ENDIF
      end do ! Was loop 10
         IF((MINCORK.GT.-1.E-4).OR.&
     &        (ABS(NORMX1)+ABS(NORMY1)+ABS(NORMZ1).EQ.0.0)) GOTO 1000
      end do ! Was loop 20
 1000 CONTINUE

!     Set all the negative basis to zero and re-normalise 
!     to put on the face of an element...
      SUM=0.0
      do ILOC=1,NLOC! Was loop 
         LOCCORDSK(ILOC)=MAX(0.0,LOCCORDSK(ILOC))
         SUM=SUM+LOCCORDSK(ILOC)
      END DO
      do ILOC=1,NLOC! Was loop 
         LOCCORDSK(ILOC)=LOCCORDSK(ILOC)/SUM
      END DO
      MATPSI=0.
      MATPSIOLD=0.
      do ILOC=1,NLOC! Was loop 
         MATPSI   =MATPSI   +LOCCORDSK(ILOC)*PSI(LOCNODSK(ILOC))
         MATPSIOLD=MATPSIOLD+LOCCORDSK(ILOC)*PSIOLD(LOCNODSK(ILOC))
      END DO
!     Exaduate difference by a factor of 100.
      MATPSI   =PSI(NOD)   +100.*(MATPSI   -PSI(NOD))
      MATPSIOLD=PSIOLD(NOD)+100.*(MATPSIOLD-PSIOLD(NOD))

!     Now correct to make sure that we get a bounded soln...
      IF(OINTER.EQ.0) THEN
         MINPSI   =PSI(NOD)
         MINPSIOLD=PSIOLD(NOD)
         MAXPSI   =PSI(NOD)
         MAXPSIOLD=PSIOLD(NOD)
      do ILOC=1,NLOC! Was loop 
            KNOD=LOCNODSK(ILOC)
            IF(KNOD.NE.NOD) THEN
!     Search around node KNOD for max and min PSI...
      do COUNT=FINDRM(KNOD),FINDRM(KNOD+1)-1! Was loop 
                  JNOD=COLM(COUNT)
                  MINPSI   =MIN(PSI(JNOD),   MINPSI)
                  MINPSIOLD=MIN(PSIOLD(JNOD),MINPSIOLD)
                  MAXPSI   =MAX(PSI(JNOD),   MAXPSI)
                  MAXPSIOLD=MAX(PSIOLD(JNOD),MAXPSIOLD)
               END DO
            ENDIF
         END DO
      ELSE IF(OINTER.EQ.1) THEN
!     Limit the values so they are between [0,1].
         MINPSI=0.0
         MAXPSI=1.0
         MINPSIOLD=0.0
         MAXPSIOLD=1.0
      ELSE IF(OINTER.EQ.2) THEN
!     OINTER=2  - allow the soln to violate local boundedness, but still [0,1]
         MINPSI=0.0
         MAXPSI=1.0
         MINPSIOLD=0.0
         MAXPSIOLD=1.0
!     Amend MATPSI so that it exsadurates difference...
!     MATPSI   = MATPSI    + 2.0*(MATPSI-PSI(NOD))
!     MATPSIOLD= MATPSIOLD + 2.0*(MATPSIOLD-PSIOLD(NOD))
!     IF(ABS(PSI(NOD)-PSI2).LT.0.25) THEN
         IF((PSI(NOD)-PSI2).GT.0.0) THEN
            IF(ABS(PSI(NOD)-MATPSI).LT.ABS(PSI(NOD)-PSI2)) THEN
               MATPSI   = RELAX*MATPSI    + (1.-RELAX)*(2.*PSI(NOD)-PSI2)
            ENDIF
         ENDIF
         IF((PSIOLD(NOD)-PSIOLD2).GT.0.0) THEN
            IF(ABS(PSIOLD(NOD)-MATPSIOLD).LT.ABS(PSIOLD(NOD)-PSIOLD2)) THEN
               MATPSIOLD= RELAX*MATPSIOLD + (1.-RELAX)*(2.*PSIOLD(NOD)-PSIOLD2)
            ENDIF
         ENDIF
      ELSE IF(OINTER.EQ.3) THEN
!     OINTER=2  - allow the soln to violate local boundedness, but still [0,1]
         MINPSI=0.0
         MAXPSI=1.0
         MINPSIOLD=0.0
         MAXPSIOLD=1.0
!     Amend MATPSI so that it exsadurates difference...
!     MATPSI   = MATPSI    + 2.0*(MATPSI-PSI(NOD))
!     MATPSIOLD= MATPSIOLD + 2.0*(MATPSIOLD-PSIOLD(NOD))
         IF((PSI(NOD)-PSI2).GT.0.0) THEN
            MATPSI   = RELAX*MATPSI    + (1.-RELAX)*(2.*PSI(NOD)-PSI2)
         ENDIF
         IF((PSIOLD(NOD)-PSIOLD2).GT.0.0) THEN
            MATPSIOLD= RELAX*MATPSIOLD + (1.-RELAX)*(2.*PSIOLD(NOD)-PSIOLD2)
         ENDIF
!     stop  34 
      ELSE IF(OINTER.EQ.4) THEN
!     OINTER=2  - allow the soln to violate local boundedness, but still [0,1]
         MINPSI=0.0
         MAXPSI=1.0
         MINPSIOLD=0.0
         MAXPSIOLD=1.0
!     Amend MATPSI so that it exsadurates difference...
         IF((PSI(NOD)-PSI2).GT.0.0) THEN
            MATPSI   = MATPSI    + (PSI(NOD)-PSI2)
         ENDIF
         IF((PSIOLD(NOD)-PSIOLD2).GT.0.0) THEN
            MATPSIOLD   = MATPSIOLD    + (PSIOLD(NOD)-PSIOLD2)
         ENDIF
      ELSE IF(OINTER.EQ.5) THEN
!     OINTER=2  - allow the soln to violate local boundedness, but still [0,1]
         MINPSI=0.0
         MAXPSI=1.0
         MINPSIOLD=0.0
         MAXPSIOLD=1.0
!     Amend MATPSI so that it exsadurates difference...
!     MATPSI   = MATPSI    + 2.0*(MATPSI-PSI(NOD))
!     MATPSIOLD= MATPSIOLD + 2.0*(MATPSIOLD-PSIOLD(NOD))
         MATPSI   = MATPSI    + (PSI(NOD)-PSI2)
         MATPSIOLD   = MATPSIOLD    + (PSIOLD(NOD)-PSIOLD2)
!     stop  34 
      ELSE
         FLAbort("3831")
      ENDIF
!     ENDIF
      MATPSI   =MAX(MIN(MATPSI,   MAXPSI),   MINPSI)
      MATPSIOLD=MAX(MIN(MATPSIOLD,MAXPSIOLD),MINPSIOLD)
!     
      RETURN

  end subroutine imatpts

!     
!     SUBROUTINE FINPTS HAS MOVED TO LEGACY_advection_diffusion_CV.F90 (CRGW)
!     
!
!
        SUBROUTINE GETSTOREELEWEI(PSI,PSIOLD,NONODS,NLOC,TOTELE,NDGLNO, &
     &     MATPSI,MATPSIOLD,FINDRM,COLM,NCOLM,BOUND,&
     &     ELEMATPSI,ELEMATWEI)
! use the stored interpolation coeffs to caclulate MATPSI,MATPSIOLD.
!     This sub finds the matrix values MATPSI for a given point on the 
!     stencil 
        REAL FRALINE
        LOGICAL BOUND
        PARAMETER(FRALINE=0.001)
        INTEGER NONODS,NLOC,TOTELE,NDGLNO(TOTELE*NLOC)
        REAL PSI(NONODS),PSIOLD(NONODS)
        INTEGER NCOLM
        INTEGER FINDRM(NONODS+1),COLM(NCOLM)
        REAL MATPSI(NCOLM),MATPSIOLD(NCOLM)
        INTEGER ELEMATPSI(NCOLM)
        REAL ELEMATWEI(NCOLM*NLOC)
!  LOCAL VARIABLES...
        INTEGER NOD,COUNT,ELEWIC,ILOC,INOD
        INTEGER KNOD,COUNT2,JNOD
        REAL RMATPSI,RMATPSIOLD
        REAL, ALLOCATABLE, DIMENSION(:)::MINPSI
        REAL, ALLOCATABLE, DIMENSION(:)::MINPSIOLD
        REAL, ALLOCATABLE, DIMENSION(:)::MAXPSI
        REAL, ALLOCATABLE, DIMENSION(:)::MAXPSIOLD
        
        if(bound) then
          ALLOCATE(MINPSI(TOTELE))
          ALLOCATE(MINPSIOLD(TOTELE))
          ALLOCATE(MAXPSI(TOTELE))
          ALLOCATE(MAXPSIOLD(TOTELE))
        
! find the max and min local to each element...
        ewrite(1, *),'going into minmaxelewic'
          CALL MINMAXELEWIC(PSI,PSIOLD,NONODS,NLOC,TOTELE,NDGLNO, &
     &     FINDRM,COLM,NCOLM,&
     &     MINPSI,MINPSIOLD,MAXPSI,MAXPSIOLD)
        ewrite(1, *),'out of minmaxelewic'
        endif
        
      do NOD=1,NONODS! Was loop 
!         print *,'nod=',nod
      do COUNT=FINDRM(NOD),FINDRM(NOD+1)-1! Was loop 
!            print *,'count=',count
            IF(NOD.NE.COLM(COUNT)) THEN
              ELEWIC=ELEMATPSI(COUNT)
!              print *,'elewic=',elewic
              RMATPSI=0.0
              RMATPSIOLD=0.0
      do ILOC=1,NLOC! Was loop 
                INOD=NDGLNO((ELEWIC-1)*NLOC+ILOC)
                RMATPSI=RMATPSI+ELEMATWEI((COUNT-1)*NLOC+ILOC)*PSI(INOD)
                RMATPSIOLD=RMATPSIOLD+ELEMATWEI((COUNT-1)*NLOC+ILOC)*PSIOLD(INOD)
              END DO
              
              RMATPSI   =PSI(NOD)   +(1./FRALINE)*(RMATPSI   -PSI(NOD))
              RMATPSIOLD=PSIOLD(NOD)+(1./FRALINE)*(RMATPSIOLD-PSIOLD(NOD))

! make locally bounded...
              if(bound) then
                 MATPSI(COUNT)   =MAX(MIN(RMATPSI,   MAXPSI(ELEWIC)),   &
     &                         MINPSI(ELEWIC))
                 MATPSIOLD(COUNT)=MAX(MIN(RMATPSIOLD,MAXPSIOLD(ELEWIC)),&
     &                         MINPSIOLD(ELEWIC))
              else
                 MATPSI(COUNT)   =RMATPSI
                 MATPSIOLD(COUNT)=RMATPSIOLD
              endif
            ENDIF
          END DO
        END DO
        RETURN
        
  end subroutine getstoreelewei
        
!
!   
!
!
        SUBROUTINE MINMAXELEWIC(PSI,PSIOLD,NONODS,NLOC,TOTELE,NDGLNO, &
     &     FINDRM,COLM,NCOLM,&
     &     MINPSI,MINPSIOLD,MAXPSI,MAXPSIOLD)
! This sub calculates the max and min values of PSI in local vacinity of 
! an element. 
        REAL FRALINE
        PARAMETER(FRALINE=0.001)
        INTEGER NONODS,NLOC,TOTELE,NDGLNO(TOTELE*NLOC)
        REAL PSI(NONODS),PSIOLD(NONODS)
        INTEGER NCOLM
        INTEGER FINDRM(NONODS+1),COLM(NCOLM)
        REAL MINPSI(TOTELE),MINPSIOLD(TOTELE),MAXPSI(TOTELE),MAXPSIOLD(TOTELE)
!  LOCAL VARIABLES...
        INTEGER NOD,COUNT,ELEWIC,ILOC,INOD
        INTEGER KNOD,COUNT2,JNOD
        REAL RMATPSI,RMATPSIOLD
        
! find the max and min local to each element...
      do ELEWIC=1,TOTELE! Was loop 
           MINPSI(ELEWIC)   =1.E+20
           MINPSIOLD(ELEWIC)=1.E+20
           MAXPSI(ELEWIC)   =-1.E+20
           MAXPSIOLD(ELEWIC)=-1.E+20
      do ILOC=1,NLOC! Was loop 
              KNOD=NDGLNO((ELEWIC-1)*NLOC+ILOC)
!     Search around node KNOD for max and min PSI...
      do COUNT2=FINDRM(KNOD),FINDRM(KNOD+1)-1! Was loop 
                  JNOD=COLM(COUNT2)
                  MINPSI(ELEWIC)   =MIN(PSI(JNOD),   MINPSI(ELEWIC))
                  MINPSIOLD(ELEWIC)=MIN(PSIOLD(JNOD),MINPSIOLD(ELEWIC))
                  MAXPSI(ELEWIC)   =MAX(PSI(JNOD),   MAXPSI(ELEWIC))
                  MAXPSIOLD(ELEWIC)=MAX(PSIOLD(JNOD),MAXPSIOLD(ELEWIC))
              END DO
           END DO
        END DO
        RETURN
        
  end subroutine minmaxelewic
        
!     
!     
!     
!     
      SUBROUTINE FINPTSSTORE(PSI,PSIOLD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
     &     MATPSI,MATPSIOLD,FINDRM,COLM,NCOLM,&
     &     XNDGLN,XNONOD, &
     &     X,Y,Z,&
     &     N,NLX,NLY,NLZ, WEIGHT,&
!     work space...
     &     FINDELE,COLELE,NCOLEL,&
     &     ELEMATPSI,ELEMATWEI,IGETSTOR,&
     &     BOUND, REFLECT)
!     This sub finds the matrix values MATPSI for a given point on the 
!     stencil 
! IF IGETSTOR=1 then get ELEMATPSI,ELEMATWEI.
      LOGICAL BOUND,REFLECT
! IF REFLECT then use a reflection condition at boundary to 
! do limiting. 
      INTEGER NONODS,NLOC,NGI,TOTELE,NDGLNO(TOTELE*NLOC)
      REAL PSI(NONODS),PSIOLD(NONODS)
      INTEGER NCOLM,NCOLEL
      INTEGER FINDRM(NONODS+1),COLM(NCOLM)
      REAL MATPSI(NCOLM),MATPSIOLD(NCOLM)
      INTEGER XNDGLN(TOTELE*NLOC),XNONOD
      REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
      REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
      REAL WEIGHT(NGI)
!     work space...
      INTEGER FINDELE(NONODS+1),COLELE(NCOLEL)
      INTEGER IGETSTOR
      INTEGER ELEMATPSI(NCOLM*IGETSTOR)
      REAL ELEMATWEI(NCOLM*NLOC*IGETSTOR)
! ELEWIC is the element to do interpolation from
! LOCCORDSK contains the weights. 
!     Local variables...
      INTEGER NOD,COUNT,NODI,NODJ,ILOC,GI,ELE
      INTEGER ELEWIC,XNOD,XNODJ
      REAL LOCCORDSK(NLOC)
      REAL NORMX1,NORMY1,NORMZ1
      REAL VOLUME,INVH,LENG
!     work space...
      REAL DETWEI(NGI),RA(NGI)
      REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
      REAL, ALLOCATABLE, DIMENSION(:)::NORMX
      REAL, ALLOCATABLE, DIMENSION(:)::NORMY
      REAL, ALLOCATABLE, DIMENSION(:)::NORMZ
      REAL, ALLOCATABLE, DIMENSION(:)::MLUM
        REAL, ALLOCATABLE, DIMENSION(:)::MINPSI
        REAL, ALLOCATABLE, DIMENSION(:)::MINPSIOLD
        REAL, ALLOCATABLE, DIMENSION(:)::MAXPSI
        REAL, ALLOCATABLE, DIMENSION(:)::MAXPSIOLD
        INTEGER, ALLOCATABLE, DIMENSION(:)::NOD2XNOD
        
        NORMX1=0.0
        NORMY1=0.0
        NORMZ1=0.0
        IF(REFLECT) THEN
!     calculate normals...********************
        ALLOCATE(NORMX(NONODS))
        ALLOCATE(NORMY(NONODS))
        ALLOCATE(NORMZ(NONODS))
        ALLOCATE(MLUM(NONODS))
      NORMX(1:NONODS) = 0.0
      NORMY(1:NONODS) = 0.0
      NORMZ(1:NONODS) = 0.0
      MLUM(1:NONODS) = 0.0
      do ELE=1,TOTELE! Was loop 
!     Calculate DETWEI,RA,NX,NY,NZ for element ELE
         CALL DETNLXR(ELE, X,Y,Z, XNDGLN, TOTELE,XNONOD,NLOC,NGI, &
     &        N,NLX,NLY,NLZ, WEIGHT, DETWEI,RA,VOLUME, .TRUE.,.FALSE., &
     &        NX,NY,NZ) 
!     
      do ILOC=1,NLOC! Was loop 
            NODI=NDGLNO((ELE-1)*NLOC+ILOC)
      do GI=1,NGI! Was loop 
               NORMX(NODI)=NORMX(NODI)+NX(ILOC,GI)*DETWEI(GI)
               NORMY(NODI)=NORMY(NODI)+NY(ILOC,GI)*DETWEI(GI)
               NORMZ(NODI)=NORMZ(NODI)+NZ(ILOC,GI)*DETWEI(GI)
               MLUM(NODI) =MLUM(NODI) +N(ILOC,GI) *DETWEI(GI)
            END DO
         END DO
      END DO
!     Renormalise
      do NODI=1,NONODS! Was loop 
         INVH=(ABS(NORMX(NODI))+ABS(NORMY(NODI))+ABS(NORMZ(NODI)))&
     &        /MLUM(NODI)
         IF(INVH.GT.1.E-5) THEN
            LENG=SQRT(NORMX(NODI)**2+NORMY(NODI)**2+NORMZ(NODI)**2)
            NORMX(NODI)=NORMX(NODI)/LENG
            NORMY(NODI)=NORMY(NODI)/LENG
            NORMZ(NODI)=NORMZ(NODI)/LENG
         ELSE
            NORMX(NODI)=0.0
            NORMY(NODI)=0.0
            NORMZ(NODI)=0.0
         ENDIF
      END DO
! In parallel need to distribute NORMX,NORMY,NORMZ
!     calculate normals...************************      
        ENDIF
        
        ALLOCATE(MINPSI(TOTELE))
        ALLOCATE(MINPSIOLD(TOTELE))
        ALLOCATE(MAXPSI(TOTELE))
        ALLOCATE(MAXPSIOLD(TOTELE))
        
        IF(BOUND) THEN
! find the max and min local to each element...
        CALL MINMAXELEWIC(PSI,PSIOLD,NONODS,NLOC,TOTELE,NDGLNO, &
     &     FINDRM,COLM,NCOLM,&
     &     MINPSI,MINPSIOLD,MAXPSI,MAXPSIOLD)
        ENDIF
      
!     
!     Calculate node element list. - CALL TO PHILNODELE MOVED TO CONSTRUCT_ADVECTION_DIFFUSION_CV

      ALLOCATE(NOD2XNOD(NONODS))
      do ELE=1,TOTELE! Was loop 
      do ILOC=1,NLOC! Was loop 
          NOD =NDGLNO((ELE-1)*NLOC+ILOC)
          XNOD=XNDGLN((ELE-1)*NLOC+ILOC)
          NOD2XNOD(NOD)=XNOD
        END DO
      END DO
!     
      do  NOD=1,NONODS! Was loop 10
         XNOD=NOD2XNOD(NOD)
!     
      do  COUNT=FINDRM(NOD),FINDRM(NOD+1)-1! Was loop 20
            NODJ=COLM(COUNT)
            XNODJ=NOD2XNOD(NODJ)
            MATPSI(COUNT)=0.
!     
            IF(NOD.NE.NODJ) THEN
               IF(REFLECT) THEN
                 NORMX1=NORMX(NOD)
                 NORMY1=NORMY(NOD)
                 NORMZ1=NORMZ(NOD)
               ENDIF
               CALL MATPTSSTORE(MATPSI(COUNT),MATPSIOLD(COUNT),NOD,&
     &              PSI,PSIOLD,NONODS,XNONOD,&
     &              NLOC,TOTELE,XNDGLN,NDGLNO,&
     &              FINDRM,COLM,NCOLM,&
     &              X(XNOD),Y(XNOD),Z(XNOD),&
     &              X(XNODJ),Y(XNODJ),Z(XNODJ),&
     &              NORMX1,NORMY1,NORMZ1,&
     &              X,Y,Z,&
!     work space...
     &              FINDELE,COLELE,NCOLEL, &
     &              MINPSI,MINPSIOLD,MAXPSI,MAXPSIOLD, &
     &              ELEWIC,LOCCORDSK,BOUND,REFLECT)
               IF(IGETSTOR.EQ.1) THEN
                  ELEMATPSI(COUNT)=ELEWIC
      do ILOC=1,NLOC! Was loop 
                    ELEMATWEI((COUNT-1)*NLOC+ILOC)=LOCCORDSK(ILOC)
                  END DO
               ENDIF
            ENDIF
      end do ! Was loop 20
      end do ! Was loop 10
      RETURN

  end subroutine finptsstore

!     
!     
!     
!     
      SUBROUTINE MATPTSSTORE(MATPSI,MATPSIOLD,NOD,&
     &     PSI,PSIOLD,NONODS,XNONOD,&
     &     NLOC,TOTELE,XNDGLN,NDGLNO,&
     &     FINDRM,COLM,NCOLM,&
     &     X1,Y1,Z1,&
     &     X2,Y2,Z2,&
     &     NORMX1,NORMY1,NORMZ1,&
     &     X,Y,Z,&
!     work space...
     &     FINDELE,COLELE,NCOLEL,&
     &     MINPSI,MINPSIOLD,MAXPSI,MAXPSIOLD,  &
     &     ELEWIC,LOCCORDSK,BOUND,REFLECT)
!     This sub calculates the value of PSI that would be at the 
!     other side of the stencil if we had a linear variation and within 
!     a single element.     
! IF BOUND then make locally bounded.
      REAL INFINY,FRALINE
      LOGICAL REFLECT
! IF REFLECT then use a reflection condition at boundary to 
! do limiting. 
      PARAMETER(INFINY=1.E+20,FRALINE=0.001)
      LOGICAL BOUND
      INTEGER NOD,NONODS,XNONOD,NLOC,TOTELE
      REAL MATPSI,MATPSIOLD,PSI(NONODS),PSIOLD(NONODS)
      INTEGER XNDGLN(NLOC*TOTELE),NDGLNO(NLOC*TOTELE)
      INTEGER NCOLM,FINDRM(NONODS+1),COLM(NCOLM)
      REAL X1,Y1,Z1,X2,Y2,Z2,NORMX1,NORMY1,NORMZ1
      REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
      INTEGER NCOLEL
      INTEGER FINDELE(NONODS+1),COLELE(NCOLEL)
      REAL MINPSI(TOTELE),MINPSIOLD(TOTELE)
      REAL MAXPSI(TOTELE),MAXPSIOLD(TOTELE)
      INTEGER ELEWIC
      REAL LOCCORDSK(NLOC)
!     
!     Local variables...
      REAL XC,YC,ZC
      REAL LOCCORDS(4)
      INTEGER LOCNODS(4),LOCNODSK(4)
      INTEGER NLOCNODS(4),NLOCNODSK(4)
      INTEGER COUNT,ELE,ILOC,KNOD,JNOD
      REAL MINCOR,MINCORK,SUM
      REAL VX,VY,VZ,T2X,T2Y,T2Z,T1X,T1Y,T1Z,DIST12,RN
      REAL REFX,REFY,REFZ,REFX2,REFY2,REFZ2
!     
      XC=X1 - FRALINE*(X2-X1)
      YC=Y1 - FRALINE*(Y2-Y1)
      ZC=Z1 - FRALINE*(Z2-Z1)
      
      IF(REFLECT) THEN
      IF(ABS(NORMX1)+ABS(NORMY1)+ABS(NORMZ1).NE.0.0) THEN
!  if (XC,YC,ZC) is outside the domain
!     The rotation matrix in 3-D is R=  
!     NX    NY    NZ
!     T1X   T1Y   T1Z
!     T2X   T2Y   T2Z
!     
            VX=X1-X2
            VY=Y1-Y2
            VZ=Z1-Z2
!     
            CALL XPROD(T2X,T2Y,T2Z, NORMX1,NORMY1,NORMZ1, VX,VY,VZ)
!     
            DIST12=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
            RN=SQRT(T2X**2+T2Y**2+T2Z**2)
            IF(RN.LT.(1.E-5)*DIST12) THEN
!     Simply have VX,VY,VZ going in the opposite direction...
               XC=X1 - VX*FRALINE
               YC=Y1 - VY*FRALINE
               ZC=Z1 - VZ*FRALINE
            ELSE
               T2X= T2X/RN
               T2Y= T2Y/RN
               T2Z= T2Z/RN
!     T1=Nx (-T2)
               CALL XPROD(T1X,T1Y,T1Z, NORMX1,NORMY1,NORMZ1, -T2X,-T2Y,-T2Z)
!     
               REFX2= NORMX1*VX + NORMY1*VY + NORMZ1*VZ
               REFY2= T1X   *VX + T1Y   *VY + T1Z   *VZ
               REFZ2= T2X   *VX + T2Y   *VY + T2Z   *VZ
!     Reflect...
               REFX2=-REFX2
!     MAP BACK USING R^T
               REFX = NORMX1*REFX2 + T1X*REFY2 + T2X*REFZ2
               REFY = NORMY1*REFX2 + T1Y*REFY2 + T2Y*REFZ2
               REFZ = NORMZ1*REFX2 + T1Z*REFY2 + T2Z*REFZ2
!     
!     (REFX,REFY,REFZ) is the reflected direction...
               XC=X1 + REFX*FRALINE
               YC=Y1 + REFY*FRALINE
               ZC=Z1 + REFZ*FRALINE
            ENDIF
      ENDIF
      ENDIF
!     
      MINCORK=-INFINY
!     
      do  COUNT=FINDELE(NOD),FINDELE(NOD+1)-1! Was loop 10
         ELE=COLELE(COUNT)
!     
         NLOCNODS(1)=NDGLNO((ELE-1)*NLOC+1)
         NLOCNODS(2)=NDGLNO((ELE-1)*NLOC+2)
         NLOCNODS(3)=NDGLNO((ELE-1)*NLOC+3)
         NLOCNODS(4)=NDGLNO((ELE-1)*NLOC+4)
!     
         LOCNODS(1)=XNDGLN((ELE-1)*NLOC+1)
         LOCNODS(2)=XNDGLN((ELE-1)*NLOC+2)
         LOCNODS(3)=XNDGLN((ELE-1)*NLOC+3)
         LOCNODS(4)=XNDGLN((ELE-1)*NLOC+4)
!     
!     C Calculate the local coord but with 4th point replaced by INOD...
!     C Find local coords LOCCORDS of point INOD corresponding to these nodes LOCNODS...
         CALL TRILOCCORDS(XC,YC,ZC, &
     &        LOCCORDS(1),LOCCORDS(2),LOCCORDS(3),LOCCORDS(4),&
!     The 4 corners of the tet...
     &        X(LOCNODS(1)),Y(LOCNODS(1)),Z(LOCNODS(1)),&
     &        X(LOCNODS(2)),Y(LOCNODS(2)),Z(LOCNODS(2)),&
     &        X(LOCNODS(3)),Y(LOCNODS(3)),Z(LOCNODS(3)),&
     &        X(LOCNODS(4)),Y(LOCNODS(4)),Z(LOCNODS(4))  )
!     
         MINCOR=MIN(LOCCORDS(1),LOCCORDS(2),&
     &        LOCCORDS(3),LOCCORDS(4))
         IF(MINCOR.GT.MINCORK) THEN 
            MINCORK=MINCOR
      do ILOC=1,NLOC! Was loop 
               LOCCORDSK(ILOC)=LOCCORDS(ILOC)
               LOCNODSK(ILOC)=LOCNODS(ILOC)
               NLOCNODSK(ILOC)=NLOCNODS(ILOC)
            END DO
            ELEWIC=ELE
         ENDIF
      end do ! Was loop 10

!     Set all the negative basis to zero and re-normalise 
!     to put on the face of an element...
      SUM=0.0
      do ILOC=1,NLOC! Was loop 
         LOCCORDSK(ILOC)=MAX(0.0,LOCCORDSK(ILOC))
         SUM=SUM+LOCCORDSK(ILOC)
      END DO
      do ILOC=1,NLOC! Was loop 
         LOCCORDSK(ILOC)=LOCCORDSK(ILOC)/SUM
      END DO
      MATPSI=0.
      MATPSIOLD=0.
!      XC=0.0
!      YC=0.0
!      ZC=0.0
      do ILOC=1,NLOC! Was loop 
         MATPSI   =MATPSI   +LOCCORDSK(ILOC)*PSI(NLOCNODSK(ILOC))
         MATPSIOLD=MATPSIOLD+LOCCORDSK(ILOC)*PSIOLD(NLOCNODSK(ILOC))
!         XC=XC+LOCCORDSK(ILOC)*X(LOCNODSK(ILOC))
!         YC=YC+LOCCORDSK(ILOC)*Y(LOCNODSK(ILOC))
!         ZC=ZC+LOCCORDSK(ILOC)*Z(LOCNODSK(ILOC))
      END DO
!     Exaduate difference by a factor of 100.
      MATPSI   =PSI(NOD)   +(1./FRALINE)*(MATPSI   -PSI(NOD))
      MATPSIOLD=PSIOLD(NOD)+(1./FRALINE)*(MATPSIOLD-PSIOLD(NOD))

!     Now correct to make sure that we get a bounded soln...
      IF(BOUND) THEN
        MATPSI   =MAX(MIN(MATPSI,   MAXPSI(ELEWIC)),   MINPSI(ELEWIC))
        MATPSIOLD=MAX(MIN(MATPSIOLD,MAXPSIOLD(ELEWIC)),MINPSIOLD(ELEWIC))
      ENDIF
!     
      RETURN

  end subroutine matptsstore

!     
!     
!     
!     
      SUBROUTINE MATPTS(MATPSI,MATPSIOLD,NOD,&
     &     PSI,PSIOLD,NONODS,&
     &     NLOC,TOTELE,NDGLNO,&
     &     FINDRM,COLM,NCOLM,&
     &     X1,Y1,Z1,&
     &     X2,Y2,Z2,&
     &     X,Y,Z,&
!     work space...
     &     FINDELE,COLELE,NCOLEL) 
!     This sub calculates the value of PSI that would be at the 
!     other side of the stencil if we had a linear variation and within 
!     a single element.     
      REAL INFINY,FRALINE
      PARAMETER(INFINY=1.E+20,FRALINE=0.001)
      INTEGER NOD,NONODS,NLOC,TOTELE
      REAL MATPSI,MATPSIOLD,PSI(NONODS),PSIOLD(NONODS)
      INTEGER NDGLNO(NLOC*TOTELE)
      INTEGER NCOLM,FINDRM(NONODS+1),COLM(NCOLM)
      REAL X1,Y1,Z1,X2,Y2,Z2
      REAL X(NONODS),Y(NONODS),Z(NONODS)
      INTEGER NCOLEL
      INTEGER FINDELE(NONODS+1),COLELE(NCOLEL)
!     
!     Local variables...
      REAL XC,YC,ZC
      REAL LOCCORDS(4),LOCCORDSK(4)
      INTEGER LOCNODS(4),LOCNODSK(4)
      INTEGER COUNT,ELE,ILOC,KNOD,JNOD
      REAL MINCOR,MINCORK,SUM
      REAL MINPSI,MINPSIOLD,MAXPSI,MAXPSIOLD
!     
      XC=X1 - FRALINE*(X2-X1)
      YC=Y1 - FRALINE*(Y2-Y1)
      ZC=Z1 - FRALINE*(Z2-Z1)
!     
      MINCORK=-INFINY
!     
      do  COUNT=FINDELE(NOD),FINDELE(NOD+1)-1! Was loop 10
         ELE=COLELE(COUNT)
!     
         LOCNODS(1)=NDGLNO((ELE-1)*NLOC+1)
         LOCNODS(2)=NDGLNO((ELE-1)*NLOC+2)
         LOCNODS(3)=NDGLNO((ELE-1)*NLOC+3)
         LOCNODS(4)=NDGLNO((ELE-1)*NLOC+4)
!     
!     C Calculate the local coord but with 4th point replaced by INOD...
!     C Find local coords LOCCORDS of point INOD corresponding to these nodes LOCNODS...
         CALL TRILOCCORDS(XC,YC,ZC, &
     &        LOCCORDS(1),LOCCORDS(2),LOCCORDS(3),LOCCORDS(4),&
!     The 4 corners of the tet...
     &        X(LOCNODS(1)),Y(LOCNODS(1)),Z(LOCNODS(1)),&
     &        X(LOCNODS(2)),Y(LOCNODS(2)),Z(LOCNODS(2)),&
     &        X(LOCNODS(3)),Y(LOCNODS(3)),Z(LOCNODS(3)),&
     &        X(LOCNODS(4)),Y(LOCNODS(4)),Z(LOCNODS(4))  )
!     
         MINCOR=MIN(LOCCORDS(1),LOCCORDS(2),&
     &        LOCCORDS(3),LOCCORDS(4))
         IF(MINCOR.GT.MINCORK) THEN 
            MINCORK=MINCOR
      do ILOC=1,NLOC! Was loop 
               LOCCORDSK(ILOC)=LOCCORDS(ILOC)
               LOCNODSK(ILOC)=LOCNODS(ILOC)
            END DO
         ENDIF
      end do ! Was loop 10

!     Set all the negative basis to zero and re-normalise 
!     to put on the face of an element...
      SUM=0.0
      do ILOC=1,NLOC! Was loop 
         LOCCORDSK(ILOC)=MAX(0.0,LOCCORDSK(ILOC))
         SUM=SUM+LOCCORDSK(ILOC)
      END DO
      do ILOC=1,NLOC! Was loop 
         LOCCORDSK(ILOC)=LOCCORDSK(ILOC)/SUM
      END DO
      MATPSI=0.
      MATPSIOLD=0.
!      XC=0.0
!      YC=0.0
!      ZC=0.0
      do ILOC=1,NLOC! Was loop 
         MATPSI   =MATPSI   +LOCCORDSK(ILOC)*PSI(LOCNODSK(ILOC))
         MATPSIOLD=MATPSIOLD+LOCCORDSK(ILOC)*PSIOLD(LOCNODSK(ILOC))
!         XC=XC+LOCCORDSK(ILOC)*X(LOCNODSK(ILOC))
!         YC=YC+LOCCORDSK(ILOC)*Y(LOCNODSK(ILOC))
!         ZC=ZC+LOCCORDSK(ILOC)*Z(LOCNODSK(ILOC))
      END DO
!     Exaduate difference by a factor of 100.
      MATPSI   =PSI(NOD)   +(1./FRALINE)*(MATPSI   -PSI(NOD))
      MATPSIOLD=PSIOLD(NOD)+(1./FRALINE)*(MATPSIOLD-PSIOLD(NOD))

!     Now correct to make sure that we get a bounded soln...
      MINPSI   =PSI(NOD)
      MINPSIOLD=PSIOLD(NOD)
      MAXPSI   =PSI(NOD)
      MAXPSIOLD=PSIOLD(NOD)
      do ILOC=1,NLOC! Was loop 
         KNOD=LOCNODSK(ILOC)
         IF(KNOD.NE.NOD) THEN
!     Search around node KNOD for max and min PSI...
      do COUNT=FINDRM(KNOD),FINDRM(KNOD+1)-1! Was loop 
               JNOD=COLM(COUNT)
               MINPSI   =MIN(PSI(JNOD),   MINPSI)
               MINPSIOLD=MIN(PSIOLD(JNOD),MINPSIOLD)
               MAXPSI   =MAX(PSI(JNOD),   MAXPSI)
               MAXPSIOLD=MAX(PSIOLD(JNOD),MAXPSIOLD)
            END DO
         ENDIF
      END DO
!     ENDIF
      MATPSI   =MAX(MIN(MATPSI,   MAXPSI),   MINPSI)
      MATPSIOLD=MAX(MIN(MATPSIOLD,MAXPSIOLD),MINPSIOLD)
!     
      RETURN

  end subroutine matpts

!     
!     
!     
!     
      SUBROUTINE PHILNODELE(NONODS,FINDELE,COLELE, &
     &     NCOLEL,MXNCOLEL, &
     &     TOTELE,NLOC,NDGLNO, &
     &     NLIST,INLIST)
      !=================================================================
      ! This sub calculates the node to element list FINDELE,COLELE
      ! 
      ! Note NLIST and INLIST are only used locally but are passed 
      ! down from parent routine where they are dynamically allocated.
      !
      ! INPUTS:
      ! ------
      ! NDGLNO  - List of global node numbers
      !
      ! OUTPUTS: 
      ! -------
      ! COLELE  - This is a list of the element numbers that each node
      !           belongs to.  So it lists all elements for node 1, then
      !           all elements for node 2, and so on...
      ! FINDELE - is the pointer to the place in COLELE that gives the 
      !           first element associated with a given global node
      ! 
      ! Called from subroutines IFINPTS and FINPTS, which are
      ! subroutines of CONSTRUCT_ADVECTION_DIFFUSION_CV
      ! 
      ! Description                                   Programmer      Date
      ! ==================================================================
      ! Original version..................................CCP     Unknown!
      ! Comments and warning when out of bounds added.....GSC   2006-08-11
      !
      !================================================================ 
      INTEGER NONODS,FINDELE(NONODS+1)
      INTEGER NCOLEL,MXNCOLEL,COLELE(MXNCOLEL)
      INTEGER TOTELE,NLOC,NDGLNO(TOTELE*NLOC)
      INTEGER NLIST(NONODS),INLIST(NONODS)
!     Local variables...
      INTEGER NOD,ELE,ILOC,COUNT, INOD
!     
      do NOD=1,NONODS! Was loop 
         NLIST(NOD)=0
         INLIST(NOD)=0
      END DO

      ! NLIST is the number of elements each node belongs to...
      do ELE=1,TOTELE! Was loop 
      do ILOC=1,NLOC! Was loop 
            INOD=NDGLNO((ELE-1)*NLOC+ILOC)
            NLIST(INOD)=NLIST(INOD)+1
         END DO
      END DO

      ! FINDELE is a pointer to the first element 
      ! associated with a given global node (NOD)
      COUNT=0
      do NOD=1,NONODS! Was loop 
         FINDELE(NOD)=COUNT+1
         COUNT=COUNT+NLIST(NOD)
      END DO
      FINDELE(NONODS+1)=COUNT+1
      NCOLEL=COUNT

      ! COLELE is a list of the element numbers each node belongs
      ! to stored in the order of the global nodes...
      ! INLIST is the element number the node belongs to.
      do ELE=1,TOTELE! Was loop 
      do ILOC=1,NLOC! Was loop 
            INOD=NDGLNO((ELE-1)*NLOC+ILOC)
            INLIST(INOD)=INLIST(INOD)+1
            IF (FINDELE(INOD)-1+INLIST(INOD).GT.MXNCOLEL) THEN
               FLAbort('COLELE ARRAY OUT OF BOUNDS--SUB:PHILNODELE')
            ENDIF
            COLELE(FINDELE(INOD)-1+INLIST(INOD))=ELE ! 
         END DO
      END DO
      RETURN

  end subroutine philnodele

!     
!     
!     
!     
      Subroutine TRILOCCORDS(Xp,Yp,Zp, &
     &     N1, N2, N3, N4, &
     &     X1,Y1,Z1, &
     &     X2,Y2,Z2, &
     &     X3,Y3,Z3, &
     &     X4,Y4,Z4  )

      Real Xp, Yp, Zp

      Real N1, N2, N3, N4
      
      Real X1,Y1,Z1
      Real X2,Y2,Z2
      Real X3,Y3,Z3
      Real X4,Y4,Z4

      Real Volume

!     calculate element volume...

      Volume = TetVolume(X1, Y1, Z1, &
     &     X2, Y2, Z2, &
     &     X3, Y3, Z3, &
     &     X4, Y4, Z4)
      
      Volume = Volume /6.0


!     vol coords...

      N1 = TetVolume(Xp, Yp, Zp, &
     &     X2, Y2, Z2, &
     &     X3, Y3, Z3, &
     &     X4, Y4, Z4) 

      N1 = N1/(6.0*Volume)

      

      N2 = TetVolume(X1, Y1, Z1, &
     &     Xp, Yp, Zp, &
     &     X3, Y3, Z3, &
     &     X4, Y4, Z4) 
      
      N2 = N2/(6.0*Volume)
      


      N3 = TetVolume(X1, Y1, Z1, &
     &     X2, Y2, Z2, &
     &     Xp, Yp, Zp, &
     &     X4, Y4, Z4) 

      N3 = N3/(6.0*Volume)

      
      N4 = TetVolume(X1, Y1, Z1, &
     &     X2, Y2, Z2, &
     &     X3, Y3, Z3, &
     &     Xp, Yp, Zp) 

      N4 = N4/(6.0*Volume)


      Return

  end subroutine triloccords

! 
!     
!     
      SUBROUTINE BOTHANONVDLIM(TOTELE,INCOME, PELE,PELEOT,&
     &     FIRORD,NOLIMI, &
     &     NCOLM,FINDRM,COLM,&
     &     OPINTE,&
     &     TDLIM,TDCEN,FVT,ETDNEW,MATPSI,&
     &     TDLIMold,TDCENold,FVTold,ETDNEWold,MATPSIold)
!     This sub calculates the limited face values TDADJ(1...SNGI) 
!     from the central difference face values TDCEN(1...SNGI)
!     using a NVD shceme.  INCOME(1...SNGI)=1 for incomming to element ELE 
!     else =0. 
!     If OPINTE=1 then optimal cross wind method for interface tracking. 
      REAL, PARAMETER::TOLER=1.0E-10
!     LIBETA is the flux limiting parameter. 
      INTEGER TOTELE
      REAL TDLIM,TDCEN,FVT
      REAL TDLIMold,TDCENold,FVTold
      INTEGER OPINTE
!     COURAT=Courant no.
      REAL INCOME
      INTEGER PELE,PELEOT
      REAL ETDNEW(TOTELE),ETDNEWold(TOTELE)
      LOGICAL FIRORD,NOLIMI
      INTEGER NCOLM,FINDRM(TOTELE+1),COLM(NCOLM)
      REAL MATPSI(NCOLM),MATPSIold(NCOLM)
!     TDMAX(PELE)=maximum of the surrounding 6 element values of element PELE.
!     TDMIN(PELE)=minimum of the surrounding 6 element values of element PELE.
!     PELEOT=element at other side of current face. 
!     ELEOT2=element at other side of the element ELEOTH. 
!     ELESID=element next to oposing current face. 
!     The elements are arranged in this order...
!     ELEOT2,ELE,   PELEOT,ELESID. 
!     Local variables...
!     This sub finds the neighbouring elements which are...
!     Suppose that this is the face IFACE. 
!     |
!     V
!     ---------------------------------------------------
!     |   ELEOT2   |   ELEOTH   |   ELE     |   ELESID   |  
!     --------------------------------------------------- 
!     TAIN         THALF       TAOUT
!     --------------------------------------------------- 
!     >TEXTIN
!     TEXOUT<
!     --------------------------------------------------- 
!     Local variables...
      INTEGER COUNT
      INTEGER ROW,COL,JCOUNT
! for NEW varibles...
      REAL UCIN,UCOU 
      REAL TUPWIN,TUPWI2,TDELE,DENOIN,CTILIN
      REAL DENOOU,CTILOU
      REAL FTILIN,FTILOU
! for OLD varibles...
      REAL UCINold,UCOUold
      REAL TUPWINold,TUPWI2old,TDELEold,DENOINold,CTILINold
      REAL DENOOUold,CTILOUold
      REAL FTILINold,FTILOUold
!     
      IF(NOLIMI) THEN
         TDLIM=TDCEN
         TDLIMold=TDCENold
         RETURN
      ENDIF
!     
      IF(FIRORD) THEN
!     Velocity is pointing into element...
         TDLIM= FVT 
!     Velocity is pointing into element...
         TDLIMold= FVTold 
         RETURN
      ENDIF
!     
!      IF(PELEOT.NE.PELE) THEN
            IF(INCOME.EQ.0.0) THEN
!     Calculate outgoing TUPWI2...
! **********outgoing flux
               ROW=PELE
               COL=PELEOT
!     - Calculate the the coln COUNT of coln NODJ...
!               JCOUNT=0
!               DO COUNT = FINDRM(ROW),FINDRM(ROW+1)-1
!                  IF( COL .EQ. COLM(COUNT) )  JCOUNT = COUNT
!               END DO
! Find JCOUNT. ***************************************
         CALL POSINMAT(JCOUNT,ROW,COL,TOTELE,FINDRM,COLM,NCOLM)
               TUPWIN=MATPSI(JCOUNT)
               TUPWI2=TUPWIN
               TUPWINold=MATPSIold(JCOUNT)
               TUPWI2old=TUPWINold
! **NEW Calculate normalisation parameters for out going velocities ********
               TDELE  =ETDNEW(PELEOT) 
               UCOU  =ETDNEW(PELE)
               DENOOU=TDELE-TUPWI2
               IF(ABS(DENOOU).LT.TOLER) DENOOU=SIGN(TOLER,DENOOU)
               CTILOU=(UCOU-TUPWI2)/DENOOU
! ----now limit...
               FTILOU=(TDCEN-TUPWI2)/DENOOU
               TDLIM= &
!     Velocity is going out of element...
     &              +( TUPWI2+BOTHNVDFUN(FTILOU,CTILOU)*DENOOU )
!     &              +(1.0-INCOME)
!     &              *( TUPWI2+BOTHNVDFUN(FTILOU,CTILOU)*DENOOU )
!
! **OLD Calculate normalisation parameters for out going velocities ********
               TDELEold  =ETDNEWold(PELEOT) 
               UCOUold  =ETDNEWold(PELE)
               DENOOUold=TDELEold-TUPWI2old
               IF(ABS(DENOOUold).LT.TOLER) DENOOUold=SIGN(TOLER,DENOOUold)
               CTILOUold=(UCOUold-TUPWI2old)/DENOOUold
! ----now limit...
               FTILOUold=(TDCENold-TUPWI2old)/DENOOUold
               TDLIMold= &
!     Velocity is going out of element...
     &              +( TUPWI2old+BOTHNVDFUN(FTILOUold,CTILOUold)*DENOOUold )
!     &              +(1.0-INCOME)
!     &              *( TUPWI2+BOTHNVDFUN(FTILOU,CTILOU)*DENOOU )
            ELSE
!     Calculate incoming TUPWIN...
! **********incoming flux
               ROW=PELEOT
               COL=PELE
!     - Calculate the the coln COUNT of coln NODJ...
!               JCOUNT=0
!               DO COUNT = FINDRM(ROW),FINDRM(ROW+1)-1
!                  IF( COL .EQ. COLM(COUNT) )  JCOUNT = COUNT
!               END DO
! Find JCOUNT. ***************************************
         CALL POSINMAT(JCOUNT,ROW,COL,TOTELE,FINDRM,COLM,NCOLM)
               TUPWIN=MATPSI(JCOUNT)
               TUPWI2=TUPWIN
               TUPWINold=MATPSIold(JCOUNT)
               TUPWI2old=TUPWINold
! **NEW Calculate normalisation parameters for incomming velocities ********
               TDELE =ETDNEW(PELE) 
               UCIN  =ETDNEW(PELEOT)
               DENOIN=TDELE-TUPWIN 
               IF(ABS(DENOIN).LT.TOLER) DENOIN=SIGN(TOLER,DENOIN)
               CTILIN=(UCIN-TUPWIN)/DENOIN
! ----now limit...
               FTILIN=(TDCEN-TUPWIN)/DENOIN
!     Velocity is pointing into element...
               TDLIM= ( TUPWIN+BOTHNVDFUN(FTILIN,CTILIN)*DENOIN )
!               TDLIM= INCOME
!     &              *( TUPWIN+BOTHNVDFUN(FTILIN,CTILIN)*DENOIN )

! **OLD Calculate normalisation parameters for incomming velocities ********
               TDELEold =ETDNEWold(PELE) 
               UCINold  =ETDNEWold(PELEOT)
               DENOINold=TDELEold-TUPWINold
               IF(ABS(DENOINold).LT.TOLER) DENOINold=SIGN(TOLER,DENOINold)
               CTILINold=(UCINold-TUPWINold)/DENOINold
! ----now limit...
               FTILINold=(TDCENold-TUPWINold)/DENOINold
!     Velocity is pointing into element...
               TDLIMold= ( TUPWINold+BOTHNVDFUN(FTILINold,CTILINold)*DENOINold )
!               TDLIMold= INCOME
!     &              *( TUPWINold+BOTHNVDFUN(FTILINold,CTILINold)*DENOINold )
            ENDIF
!      ELSE
!C     Calculate normalisation parameters for incomming velocities ********
!         TUPWIN=ETDNEW(PELE)
!         UCIN  =ETDNEW(PELE)
!         DENOIN=1.
!         CTILIN=0.
!c     Calculate normalisation parameters for out going velocities ********
!         TUPWI2=ETDNEW(PELE)
!         UCOU  =ETDNEW(PELE)
!         DENOOU=1.
!         CTILOU=0.
!         
!         FTILIN=(TDCEN-TUPWIN)/DENOIN
!         FTILOU=(TDCEN-TUPWI2)/DENOOU
!C     Velocity is pointing into element...
!         TDLIM= INCOME
!     &        *( TUPWIN+BOTHNVDFUN(FTILIN,CTILIN)*DENOIN )
!C     Velocity is going out of element...
!     &        +(1.0-INCOME)
!     &        *( TUPWI2+BOTHNVDFUN(FTILOU,CTILOU)*DENOOU )
!      ENDIF
!     
      RETURN

  end subroutine bothanonvdlim

!         
!     
!
!
      REAL FUNCTION BOTHNVDFUN(UFTILD,UCTILD)
      REAL LIMGRA
      PARAMETER(LIMGRA=2.0)
!     LIMGRA=limit gradient (=2 for TVD) (=3 recommended else where).
!     TDCEN= the suggested value on the face (obtain by a central scheme)
!     
      REAL UFTILD,UCTILD
!     Alter UFTILD so it is in non-oscilatoty region of NVD.
!     Local variables...
      REAL UFT
!     
      UFT=UFTILD
      IF((UCTILD.GT.1.0).OR.(UCTILD.LT.0.0)) THEN  
         UFT=UCTILD
      ELSE
!     Place inside the triangle. 
!     Place above line...
         UFT=MAX(UFT,UCTILD)
!     Place in TVD region and below 1.0...
!     UFT=MIN(UFT,max(min(0.5,3.*UCTILD),2.*UCTILD), 1.0)
!     UFT=MIN(UFT,3.*UCTILD,1.0)
!     UFT=MIN(UFT,UCTILD+8.*UCTILD**2,3.*UCTILD,1.0)
!     UFT=MIN(UFT,4.*UCTILD-4.*UCTILD**2,1.0)
!     UFT=MIN(UFT,4.*UCTILD-3.*UCTILD**2,1.0)
!     UFT=MIN(UFT,0.75*(4.*UCTILD-2.5*UCTILD**2),1.0)
!     UFT=MIN(UFT,0.75*(4.*UCTILD-2.*UCTILD**2),1.0)
!     one of the best...
!     UFT=MIN(UFT,3.*UCTILD-1.5*UCTILD**2,1.0)
!     UFT=MIN(UFT,UCTILD+8.*UCTILD**2,1.0)
!     one of the best...
         UFT=MIN(UFT,LIMGRA*UCTILD,1.0)
!     Place in below 1.0
!     UFT=MIN(UFT,1.0)
!     UFT=MIN(1.0,MAX(0.0,UFT),LIMGRA*UCTILD) 
!     UFT= MIN(1.0,MAX(0.0,UFT),LIMGRA*UCTILD) 
!     
!     IF(UFT.GT.0.2*UCTILD+0.8) UFT=0.2*UCTILD+0.8
      ENDIF
      BOTHNVDFUN=UFT
      RETURN

  end function bothnvdfun

!     
!         
!      
!     


!     
!     SUBROUTINE ADVOCEANMOM HAS MOVED TO LEGACY_Advection_diffusion_CV.F90 (GSC, 14/08/06)
!     SUBROUTINE ADVOCEAN2   HAS MOVED TO LEGACY_Advection_diffusion_CV.F90 (GSC, 14/08/06)
!     SUBROUTINE ADVOCEAN    HAS MOVED TO LEGACY_Advection_diffusion_CV.F90 (GSC, 14/08/06)
!     SUBROUTINE CVCTS2      HAS BEEN REPLACED BY NEW CODE STRUCTURE IN Divergence_Matrix_CV.F90
!     SUBROUTINE CVCTS       HAS BEEN REPLACED BY NEW CODE STRUCTURE IN Divergence_Matrix_CV.F90
!     
      SUBROUTINE ANONVDLIM(TOTELE,COURAT,TDLIM,TDCEN, INCOME,CONVOL1,&
     &     CONVOL2,ETDNEW,MATPSI,FIRORD,NOLIMI,NCOLM,FINDRM,COLM,&
     &     OPINTE,GPSIL,FPSIL,ULTRAC  )
!     =============================================================
!     This sub calculates the limited face value TDLIM from the 
!     central-difference face value TDCEN between two neighbouring
!     nodes CONVOL1 and CONVOL2 using an NVD scheme.  
!     
!     INPUTS:-
!     
!     TOTELE - Total number of control volumes 
!     COURAT - Courant number (v*dt/dx) (-ve if normal limiting)
!     TDCEN  - Guess at field value at the face
!     INCOME - 1. => flux from CONVOL1 to CONVOL2; 0. => vice-versa
!     CONVOL1- Node on one side of the CV face
!     CONVOL2- Node on the other side of the CV face
!     ETDNEW - Field variable array
!     MATPSI - Field value at far-upwind location
!     FIRORD - Flag set to .true. for first-order advection
!     NOLIMI - Flag set to .true. for no limiting
!     NCOLM  - Max. number of columns in matrix
!     FINDRM - List of nonzero components of BIGM for each row
!     COLM   - List of column numbers for components of BIGM
!
!     OPINTE - Options for interface tracking scheme under development
!     GPSIL  - Options for interface tracking scheme under development
!     FPSIL  - Options for interface tracking scheme under development
!     
!     OUTPUTS:-
!     
!     TDLIM  - Limited value of field at the control volume face
!     
!     ------------------------------------------------------------------
!     
!     Notation used in this subroutine:
!     
!     Suppose that this is the face through which T is advected
!                               |
!                               V
!      ---------------------------------------------------
!     |            |  CONVOL1   |   CONVOL2  |            |  
!      --------------------------------------------------- 
!     
!               If flow is in this direction --->
!      ---------------------------------------------------
!     |   TUPWIN   |   TDONOR   |   TDNWIN   |            |  
!      --------------------------------------------------- 
!     
!               <--- If flow is in this direction
!      ---------------------------------------------------
!     |            |   TDNWIN   |   TDONOR   |  TUPWIN    |  
!      --------------------------------------------------- 
!     
!     TUPWIN is the upwind value of the field variable T; it is 
!     determined by extrapolating the field-variable variation in the 
!     upwind direction along a line connecting CONVOL1 and CONVOL2.  This 
!     is done in the subroutine IFINPTS and the values stored in the 
!     matrix MATPSI.
!     
!     This routine first defines the upwind, downwind, and donor values
!     depending on the direction of the flow (INCOMING).  It then
!     calculates the normalised donor value using:
!     
!                                DONOR_VALUE - UPWIND_VALUE
!     NORMALISED_DONOR_VALUE = -----------------------------
!                              DOWNWIND_VALUE - UPWIND_VALUE
!     
!     This is then passed to a Normalised Variable Diagram (NVD)
!     function which defines the normalised value at the face as a
!     simple function of the normalised donor value.  The limited
!     face value is then calculated from:
!     
!     FACE_VALUE = UPWIND_VALUE + NORMALISED_FACE_VALUE
!                                   *(DOWNWIND_VALUE - UPWIND_VALUE)
!     
!     The normalised face value is computed in NVDFUNNEW, using either
!     a Univeral-limiter scheme, or a Hyper-compressive advection
!     scheme.  The compressive advection scheme is used if a positive
!     courant number is input as COURAT; COURAT should be <0. to disable 
!     compressive advection and use the universal limiter instead.
!     
!     This scheme now also has an added option for increasing the
!     compressiveness of the advection scheme!  This should only be used
!     when advecting volume/mass fractions for interface tracking.  To
!     use, set the parameter ULTRAC to true.
!     
!     Called from CONSTRUCT_ADVECTION_DIFFUSION_CV
!     
!     Description                                   Programmer      Date
!     ==================================================================
!     Subroutine re-written and commented..............GSC    2006-03-10
!     ULTRA-Compressive option added...................GSC    2006-03-15
!     
!     Inputs and outputs
      INTEGER TOTELE
      REAL COURAT
      REAL TDLIM,TDCEN
      REAL INCOME
      INTEGER CONVOL1,CONVOL2
      REAL ETDNEW(TOTELE)
      REAL MATPSI(NCOLM)
      LOGICAL FIRORD,NOLIMI
      INTEGER NCOLM,FINDRM(TOTELE+1),COLM(NCOLM)
      INTEGER OPINTE
      REAL GPSIL,FPSIL

!     Local variables...
      INTEGER DONNODE           !The node number of the donor CV
      INTEGER RECNODE           !The node number of the receptor CV
      INTEGER COUNT             !Index of components in BIGM
      INTEGER JCOUNT            !Component number within BIGM (for upwind)
      INTEGER ROW,COL           !Row and column index of BIGM
      REAL TDONOR               !Donor-node value
      REAL TUPWIN,TDNWIN        !Up- and Downwind node values
      REAL TDENOM               !Denominator in normalization
      REAL CTILDE               !Normalized donor value
      REAL FTILDE               !Normalized face value (GUESS)
      REAL TOLER
      PARAMETER(TOLER=1.0E-10) 
      LOGICAL ULTRAC            !Set to true for extremely compressive

!     Return with non-limited face value if no limiting desired
      IF(NOLIMI) THEN
         TDLIM=TDCEN
         RETURN
      ENDIF

!     If flux is from one node (CV) to another
      IF(CONVOL2.NE.CONVOL1) THEN

!        Define which of the two CVs is the donor and which is the 
!        receptor (downwind), depending on direction of flux...
         IF(INCOME.EQ.0.0) THEN !Outgoing
            DONNODE = CONVOL1   !Donor
            RECNODE = CONVOL2   !Receptor
         ELSE                   !Incoming
            DONNODE = CONVOL2   !Donor
            RECNODE = CONVOL1   !Receptor
         ENDIF

!        New test using cross-stream terms in limiting fluxes
         IF (OPINTE.EQ.1) THEN
            
!           Define the upwind, downwind and donor values...
            TDNWIN=ETDNEW(RECNODE)
            TUPWIN=FPSIL
            TDONOR=GPSIL

!        Old 1D-method for defining the upwind, downwind and donor values.
         ELSE

!           Define the row and column of BIGM corresponding to the
!           donor and receptor (downwind) nodes (elements)...
            ROW=DONNODE
            COL=RECNODE

!           Find the component of MATPSI that corresponds to the 
!           ROW and COLumn of the receptor and donor nodes...
            JCOUNT=0
      do COUNT = FINDRM(ROW),FINDRM(ROW+1)-1! Was loop 
               IF( COL .EQ. COLM(COUNT) )  JCOUNT = COUNT
            END DO

!           Use this to define the upwind "node" value, 
!           which is stored in the MATPSI matrix... If it can't be
!           found just use the value at CONVOL1 node...
            IF(JCOUNT.EQ.0) THEN
               TUPWIN=ETDNEW(CONVOL1) 
            ELSE
               TUPWIN=MATPSI(JCOUNT) !Far-upwind value
            ENDIF

!           Define the downwind and donor (central) node values...
            TDNWIN=ETDNEW(RECNODE)
            TDONOR=ETDNEW(DONNODE)

         ENDIF

!        --------------------------------------------------------------
!        New VOF-type approach to limiting.  Adjust the up and downwind
!        values for partially-full control volumes depending on the
!        presumed location of the interface...
         IF (ULTRAC.AND.COURAT.GT.0.) THEN
            IF (TDONOR.LT.MAX(TDNWIN,TUPWIN) & !Partially full control volume
               & .AND.TDONOR.GT.MIN(TDNWIN,TUPWIN)) THEN 
               IF (TDNWIN.LT.TUPWIN) THEN !More full CV is UPWIND
                  TDNWIN = 0.
                  TUPWIN = 1.
               ELSE             !More full CV is DOWNWIND
                  TDNWIN = 1.
                  TUPWIN = 0.
               ENDIF
            ENDIF
         ENDIF
!        -------------------------------------- Added by GSC 2006-03-09

!        Calculate the normalised donor node value...
         TDENOM=TDNWIN-TUPWIN
         IF(ABS(TDENOM).LT.TOLER) TDENOM=SIGN(TOLER,TDENOM)
         CTILDE=(TDONOR-TUPWIN)/TDENOM

!     If the flux is to the same node (via reflection from boundary)
      ELSE

!        Calculate normalisation parameters... 
         TUPWIN=ETDNEW(CONVOL1)
         TDONOR=ETDNEW(CONVOL1)
         TDENOM=1.
         CTILDE=0.

      ENDIF

!     Calculate the limited face-value for 1st-order scheme
      IF(FIRORD) THEN

         TDLIM=TDONOR

!     Calculate the limited face-value for higher-order schemes
      ELSE

!        Make a guess at the normalized face values
         FTILDE=(TDCEN-TUPWIN)/TDENOM

!        Compute the limited face value from the NVD
         TDLIM=TUPWIN+NVDFUNNEW(FTILDE,CTILDE,COURAT)*TDENOM

      ENDIF
      
  end subroutine anonvdlim

      SUBROUTINE ONVDLIM(TOTELE,&
     &     TDLIM,TDCEN, INCOME, PELE,PELEOT,&
     &     ETDNEW,TDMIN,TDMAX,FIRORD,NOLIMI) 
!     This sub calculates the limited face values TDADJ(1...SNGI) 
!     from the central difference face values TDCEN(1...SNGI)
!     using a NVD shceme.  INCOME(1...SNGI)=1 for incomming to element ELE 
!     else =0. 
      REAL, parameter::TOLER=1.0E-10
!     LIBETA is the flux limiting parameter. 
      INTEGER TOTELE
      REAL TDLIM,TDCEN
      REAL INCOME
      INTEGER PELE,PELEOT
      REAL ETDNEW(TOTELE),TDMIN(TOTELE),TDMAX(TOTELE)
      LOGICAL FIRORD,NOLIMI
!     TDMAX(PELE)=maximum of the surrounding 6 element values of element PELE.
!     TDMIN(PELE)=minimum of the surrounding 6 element values of element PELE.
!     PELEOT=element at other side of current face. 
!     ELEOT2=element at other side of the element ELEOTH. 
!     ELESID=element next to oposing current face. 
!     The elements are arranged in this order...
!     ELEOT2,ELE,   PELEOT,ELESID. 
!     Local variables...
!     This sub finds the neighbouring elements which are...
!     Suppose that this is the face IFACE. 
!     |
!     V
!     ---------------------------------------------------
!     |   ELEOT2   |   ELEOTH   |   ELE     |   ELESID   |  
!     --------------------------------------------------- 
!     TAIN         THALF       TAOUT
!     --------------------------------------------------- 
!     >TEXTIN
!     TEXOUT<
!     --------------------------------------------------- 
!     Local variables...
      REAL UCIN,UCOU 
      REAL TUPWIN,TUPWI2,TDELE,DENOIN,CTILIN
      REAL DENOOU,CTILOU
      REAL FTILIN,FTILOU

      IF(NOLIMI) THEN
         TDLIM=TDCEN
         RETURN
      ENDIF
!     
!     IF(PELEOT.NE.0) THEN
      IF(PELEOT.NE.PELE) THEN
         IF(ETDNEW(PELEOT).GT.ETDNEW(PELE)) THEN
            TUPWIN=TDMAX(PELEOT)
            TUPWI2=TDMIN(PELE)
         ELSE
            TUPWIN=TDMIN(PELEOT)
            TUPWI2=TDMAX(PELE) 
         ENDIF
!     Calculate normalisation parameters for incomming velocities ********
         TDELE =ETDNEW(PELE) 
         DENOIN=TDELE-TUPWIN 
         IF(ABS(DENOIN).LT.TOLER) DENOIN=SIGN(TOLER,DENOIN)
         UCIN  =ETDNEW(PELEOT)
         CTILIN=(UCIN-TUPWIN)/DENOIN
!     Calculate normalisation parameters for out going velocities ********
         TDELE  =ETDNEW(PELEOT) 
         DENOOU=TDELE-TUPWI2
         IF(ABS(DENOOU).LT.TOLER) DENOOU=SIGN(TOLER,DENOOU)
         UCOU  =ETDNEW(PELE)
         CTILOU=(UCOU-TUPWI2)/DENOOU
      ELSE
!     Calculate normalisation parameters for incomming velocities ********
         TUPWIN=ETDNEW(PELE)
         UCIN  =ETDNEW(PELE)
         DENOIN=1.
         CTILIN=0.
!     Calculate normalisation parameters for out going velocities ********
         TUPWI2=ETDNEW(PELE)
         UCOU  =ETDNEW(PELE)
         DENOOU=1.
         CTILOU=0.
      ENDIF
!     
!     
      IF(FIRORD) THEN
!     Velocity is pointing into element...
         TDLIM= INCOME*UCIN &
!     Velocity is going out of element...
     &        +(1.0-INCOME)*UCOU
      ELSE
         FTILIN=(TDCEN-TUPWIN)/DENOIN
         FTILOU=(TDCEN-TUPWI2)/DENOOU
!     Velocity is pointing into element...
         TDLIM= INCOME&
     &        *( TUPWIN+NVDFUNNEW(FTILIN,CTILIN,-1.0)*DENOIN )&
!     Velocity is going out of element...
     &        +(1.0-INCOME)&
     &        *( TUPWI2+NVDFUNNEW(FTILOU,CTILOU,-1.0)*DENOOU )
      ENDIF
!     

  end subroutine onvdlim
        
      FUNCTION NVDFUNNEW(  UF, UC, COURAT )
!     ======================================
!     The function computes NVDFUNNEW, the normalised value of the 
!     advected variable on the face of the control volume, based on 
!     the normalised value of the advected variable in the donor CV,
!     UC, and the high-order estimate of the face value UF. 
!
!     NVDFUNNEW is limited so that it is in the non-oscillatory 
!     region of normalised variable diagram (NVD).
!
!     This version can also use an extremely compressive advection 
!     scheme, appropriate for advecting volume/mass fractions in
!     multi-material problems.  This scheme is applied if the courant
!     number (COURAT) is positive.  IF YOU DO NOT WANT COMPRESSIVE
!     ADVECTION: call this function with COURAT=-1.
!
!     ------------------------------------------------------------------
!     INPUTS:
!
!     UC     - Normalised advected field value at donor CV (node)
!     UF     - High-order estimate of normalised value at CV face
!     COURAT - Local courant number (flow_vel*dt/dx)
!              (set to -ve number if not using interface tracking scheme)
!
!     OUTPUTS:
!
!     NVDFUNNEW - Limited normalised advected field value at CV face
!     ------------------------------------------------------------------
!
!     The function is based on the flux limiting function in the 
!     paper: Riemann solvers on 3-D unstructured meshes for 
!     radiation transport, written by Dr C.C. Pain et al.
!
!     See also the paper: The ULTIMATE conservative difference scheme 
!     applied to unsteady 1D advection, Leonard B. P., 1991, 
!     Computer Methods in Applied Mechanics and Engineering, 88, 17-74.
!
!     Called from ANONVDLIM
!
!     Description                                   Programmer      Date
!     ==================================================================
!     Original subroutine .............................CCP    2002-01-21
!     Interface tracking option added..................GSC    2006-03-14
!
      REAL NVDFUNNEW,COURAT
      REAL UC, UF

!     XI is the parameter in equation 38 of the Riemann paper
!     If XI is equal to 2 then this corresponds to a TVD condition 
!     in 1-D, a value of XI equal to 3 has been recommended elsewhere.
      REAL XI
      PARAMETER( XI = 2.0 )
!       PARAMETER( XI = 1.0 )

!     Local variables
      REAL TILDEUF, MAXUF

!     For the region 0<UC<1 on the NVD, define the limiter
      IF( ( UC .GT. 0.0 ) .AND. ( UC .LT. 1.0 ) ) THEN

!        For the interface tracking scheme, use Hyper-C
!        (NOTE: at high courant numbers limit is defined by XI,
!               not the Hyper-C scheme.)
         IF(COURAT.GT.0.0) THEN
            TILDEUF = MIN( 1.0, max(UC/COURAT,XI*UC) )

!        For the normal limiting
         ELSE
            MAXUF = MAX( 0.0, UF )
            TILDEUF = MIN( 1.0, XI*UC, MAXUF )
         ENDIF

!     Outside the region 0<UC<1 on the NVD, use first-order upwinding
      ELSE
         TILDEUF = UC
      ENDIF

      NVDFUNNEW = TILDEUF

  end function nvdfunnew
  
!     
!     
!     
!     
      SUBROUTINE ONVDLIM_C(COURANT,TOTELE, &
           TDLIM,TDCEN, INCOME, PELE,PELEOT, &
           ETDNEW,TDMIN,TDMAX,FIRORD,NOLIMI) 
!     This sub calculates the limited face values TDADJ(1...SNGI) 
!     from the central difference face values TDCEN(1...SNGI)
!     using a NVD shceme.  INCOME(1...SNGI)=1 for incomming to element ELE 
!     else =0. 
      REAL, parameter::TOLER=1.0E-10
!     LIBETA is the flux limiting parameter. 
      INTEGER TOTELE
      REAL TDLIM,TDCEN,COURANT
      REAL INCOME
      INTEGER PELE,PELEOT
      REAL ETDNEW(TOTELE),TDMIN(TOTELE),TDMAX(TOTELE)
      LOGICAL FIRORD,NOLIMI
!     TDMAX(PELE)=maximum of the surrounding 6 element values of element PELE.
!     TDMIN(PELE)=minimum of the surrounding 6 element values of element PELE.
!     PELEOT=element at other side of current face. 
!     ELEOT2=element at other side of the element ELEOTH. 
!     ELESID=element next to oposing current face. 
!     The elements are arranged in this order...
!     ELEOT2,ELE,   PELEOT,ELESID. 
!     Local variables...
!     This sub finds the neighbouring elements which are...
!     Suppose that this is the face IFACE. 
!     |
!     V
!     ---------------------------------------------------
!     |   ELEOT2   |   ELEOTH   |   ELE     |   ELESID   |  
!     --------------------------------------------------- 
!     TAIN         THALF       TAOUT
!     --------------------------------------------------- 
!     >TEXTIN
!     TEXOUT<
!     --------------------------------------------------- 
!     Local variables...
      REAL UCIN,UCOU 
      REAL TUPWIN,TUPWI2,TDELE,DENOIN,CTILIN
      REAL DENOOU,CTILOU
      REAL FTILIN,FTILOU

      IF(NOLIMI) THEN
         TDLIM=TDCEN
         RETURN
      ENDIF
!     
!     IF(PELEOT.NE.0) THEN
      IF(PELEOT.NE.PELE) THEN
         IF(ETDNEW(PELEOT).GT.ETDNEW(PELE)) THEN
            TUPWIN=TDMAX(PELEOT)
            TUPWI2=TDMIN(PELE)
         ELSE
            TUPWIN=TDMIN(PELEOT)
            TUPWI2=TDMAX(PELE) 
         ENDIF
!    Calculate normalisation parameters for incomming velocities ********
         TDELE =ETDNEW(PELE) 
         DENOIN=TDELE-TUPWIN 
         IF(ABS(DENOIN).LT.TOLER) DENOIN=SIGN(TOLER,DENOIN)
         UCIN  =ETDNEW(PELEOT)
         CTILIN=(UCIN-TUPWIN)/DENOIN
!     Calculate normalisation parameters for out going velocities ********
         TDELE  =ETDNEW(PELEOT) 
         DENOOU=TDELE-TUPWI2
         IF(ABS(DENOOU).LT.TOLER) DENOOU=SIGN(TOLER,DENOOU)
         UCOU  =ETDNEW(PELE)
         CTILOU=(UCOU-TUPWI2)/DENOOU
      ELSE
!     Calculate normalisation parameters for incomming velocities ********
         TUPWIN=ETDNEW(PELE)
         UCIN  =ETDNEW(PELE)
         DENOIN=1.
         CTILIN=0.
!     Calculate normalisation parameters for out going velocities ********
         TUPWI2=ETDNEW(PELE)
         UCOU  =ETDNEW(PELE)
         DENOOU=1.
         CTILOU=0.
      ENDIF
!     
!     
      IF(FIRORD) THEN
!     Velocity is pointing into element...
         TDLIM= INCOME*UCIN &
!     Velocity is going out of element...
              +(1.0-INCOME)*UCOU
      ELSE
         FTILIN=(TDCEN-TUPWIN)/DENOIN
         FTILOU=(TDCEN-TUPWI2)/DENOOU
!     Velocity is pointing into element...
         TDLIM= INCOME &
              *( TUPWIN+NVDFUNNEW_C(FTILIN,CTILIN,COURANT)*DENOIN ) &
!     Velocity is going out of element...
              +(1.0-INCOME) &
              *( TUPWI2+NVDFUNNEW_C(FTILOU,CTILOU,COURANT)*DENOOU )
      ENDIF
!     
      RETURN

  end subroutine onvdlim_c

!     
!     
!     
!     
!***********************************************************************
      FUNCTION NVDFUNNEW_C(  UF, UC, COURAT )
!     ======================================
!     The function computes NVDFUNNEW, the normalised value of the 
!     advected variable on the face of the control volume, based on 
!     the normalised value of the advected variable in the donor CV,
!     UC, and the high-order estimate of the face value UF. 
!
!     NVDFUNNEW is limited so that it is in the non-oscillatory 
!     region of normalised variable diagram (NVD).
!
!     This version can also use an extremely compressive advection 
!     scheme, appropriate for advecting volume/mass fractions in
!     multi-material problems.  This scheme is applied if the courant
!     number (COURAT) is positive.  IF YOU DO NOT WANT COMPRESSIVE
!     ADVECTION: call this function with COURAT=-1.
!
!     ------------------------------------------------------------------
!     INPUTS:
!
!     UC     - Normalised advected field value at donor CV (node)
!     UF     - High-order estimate of normalised value at CV face
!     COURAT - Local courant number (flow_vel*dt/dx)
!              (set to -ve number if not using interface tracking scheme)
!
!     OUTPUTS:
!
!     NVDFUNNEW - Limited normalised advected field value at CV face
!     ------------------------------------------------------------------
!
!     The function is based on the flux limiting function in the 
!     paper: Riemann solvers on 3-D unstructured meshes for 
!     radiation transport, written by Dr C.C. Pain et al.
!
!     See also the paper: The ULTIMATE conservative difference scheme 
!     applied to unsteady 1D advection, Leonard B. P., 1991, 
!     Computer Methods in Applied Mechanics and Engineering, 88, 17-74.
!
!     Called from ANONVDLIM
!
!     Description                                   Programmer      Date
!     ==================================================================
!     Original subroutine .............................CCP    2002-01-21
!     Interface tracking option added..................GSC    2006-03-14
!
!************************************************************************
      REAL NVDFUNNEW_C,COURAT
      REAL UC, UF

!     XI is the parameter in equation 38 of the Riemann paper
!     If XI is equal to 2 then this corresponds to a TVD condition 
!     in 1-D, a value of XI equal to 3 has been recommended elsewhere.
      REAL XI
!      PARAMETER( XI = 2.0 )
!      PARAMETER( XI = 5.0 )
      PARAMETER( XI = 3.0 )
!      PARAMETER( XI = 10.0 )
!       PARAMETER( XI = 1.0 )

!     Local variables
      REAL TILDEUF, MAXUF, MINUF

!     For the region 0<UC<1 on the NVD, define the limiter
      IF( ( UC .GT. 0.0 ) .AND. ( UC .LT. 1.0 ) ) THEN
!      IF( ( UC .GT. 0.0 ) .AND. ( UC .LT. 2.0 ) ) THEN

!        For the interface tracking scheme, use Hyper-C
!        (NOTE: at high courant numbers limit is defined by XI,
!               not the Hyper-C scheme.)
            MAXUF = MIN( 1.0, max(UC/max(COURAT,1.e-20),XI*UC) )
            MINUF = 0.0
            TILDEUF = MAX(MIN( UF, MAXUF ), MINUF)

!     Outside the region 0<UC<1 on the NVD, use first-order upwinding
      ELSE
!         TILDEUF = max(min(1.0+2.*uc,uf),UC)
         TILDEUF = UC
      ENDIF

      NVDFUNNEW_C = TILDEUF


  end function nvdfunnew_c
 
      SUBROUTINE GAUSSILOC( FINDGPTS, COLGPTS, NCOLGPTS,&
     &     NEILOC,   NLOC,    SVNGI     )
!     ----------------------------------------------------    
!     
!     - this subroutine calculates FINDGPTS,COLGPTS,NCOLGPTS 
!     - which contains given a local node ILOC the Gauss pts 
!     - that are used to integrate around this local node. 
!     
!     -------------------------------
!     - date last modified : 12/02/2002
!     -------------------------------
!     
!     
      INTEGER NLOC, SVNGI
!     
      INTEGER FINDGPTS(NLOC+1)
!     
!     - We have overestimated the size of COLGPTS.
!     
      INTEGER COLGPTS(NLOC*SVNGI)
      INTEGER NEILOC(NLOC,SVNGI)
!     
      INTEGER NCOLGPTS
!     
      INTEGER ILOC, COUNT, GI
!     
      COUNT = 0
!     
      do ILOC = 1, NLOC! Was loop 
!     
         FINDGPTS( ILOC ) = COUNT+1
!     
      do GI = 1, SVNGI! Was loop 
!     
            IF( NEILOC(ILOC,GI) .NE. 0 ) THEN
!     
               COUNT = COUNT + 1
               COLGPTS(COUNT) = GI
!     
            END IF
!     
         END DO
!     
      END DO
!     
      FINDGPTS( NLOC+1 ) = COUNT+1
      NCOLGPTS = COUNT
       
  end subroutine gaussiloc

!     
!
!     
      SUBROUTINE GRADBOUND(NODC,NODD,XNODC,XNODD,T,TOLD,&
     &     X,Y,Z,&
     &     TMIN,TMAX,TOLDMIN,TOLDMAX, &
     &     DIFXT, DIFYT, DIFZT, &
     &     DIFXTOLD, DIFYTOLD, DIFZTOLD,&
     &     NONODS,XNONOD,&
     &     TMAX2,TMIN2,TOLDMAX2,TOLDMIN2 )
!     This sub alters MAXT,MINT,MAXTOLD,MINTOLD in 
!     order to perform anisotropic limiting.
!     The ordering is NODU, NODC, NODD
!     ^ this is where the 
!     face is for information going from left to right. 
!     TMAX2,TMIN2,TOLDMAX2,TOLDMIN2 are the original max and min. 
      INTEGER NONODS,XNONOD,NODC,NODD,XNODC,XNODD
      REAL T(NONODS),TOLD(NONODS)
      REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
      REAL TMIN(NONODS),TMAX(NONODS)
      REAL TOLDMIN(NONODS),TOLDMAX(NONODS)
      REAL DIFXT(NONODS), DIFYT(NONODS), DIFZT(NONODS)
      REAL DIFXTOLD(NONODS),DIFYTOLD(NONODS),DIFZTOLD(NONODS)
      REAL TMIN2(NONODS),TMAX2(NONODS)
      REAL TOLDMIN2(NONODS),TOLDMAX2(NONODS)
!     Local variables
      REAL UD,VD,WD,UDOTGRAD,UDOTGRADOLD
      REAL NEWT,NEWTOLD
!     
      UD=X(XNODD)-X(XNODC)
      VD=Y(XNODD)-Y(XNODC)
      WD=Z(XNODD)-Z(XNODC)
!     
      UDOTGRAD   =UD*DIFXT(NODC)+ VD*DIFYT(NODC) &
     &     +WD*DIFZT(NODC) 
      UDOTGRADOLD=UD*DIFXTOLD(NODC)+ VD*DIFYTOLD(NODC) &
     &     +WD*DIFZTOLD(NODC)
!     
!     generated values of T from gradient..
      NEWT=T(NODD) - 2.*UDOTGRAD
      NEWTOLD=TOLD(NODD) - 2.*UDOTGRADOLD
!     
      TMIN(NODC)=MIN(T(NODC),T(NODD),NEWT) 
      TMAX(NODC)=MAX(T(NODC),T(NODD),NEWT)
!     
      TOLDMIN(NODC)=MIN(TOLD(NODC),TOLD(NODD),NEWTOLD)
      TOLDMAX(NODC)=MAX(TOLD(NODC),TOLD(NODD),NEWTOLD)
!     Amend to make soln bounded...
      TMIN(NODC)=MAX( TMIN(NODC), TMIN2(NODC) )
      TMAX(NODC)=MIN( TMAX(NODC), TMAX2(NODC) )
!     
      TOLDMIN(NODC)=MAX( TOLDMIN(NODC), TOLDMIN2(NODC) )
      TOLDMAX(NODC)=MIN( TOLDMAX(NODC), TOLDMAX2(NODC) )
!     
      RETURN

  end subroutine gradbound

!     
!     
!     
!     
      SUBROUTINE GETGXYZ(T,TOLD,&
     &     DIFXT, DIFYT, DIFZT, &
     &     DIFXTOLD, DIFYTOLD, DIFZTOLD,&
     &     ML, &
     &     TOTELE, &
     &     X,Y,Z,&
     &     N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI, &
     &     D3,DCYL,&
     &     NDGLNO,NONODS,&
     &     XNDGLN,XNONOD, &
     &     DETWEI,RA,&
     &     NX,NY,NZ   )
!     This sub gets the gradients DIFXT,DIFXOLT ETC 
!     of T and TOLD.
      INTEGER NONODS,XNONOD,NLOC,TOTELE,NGI
      REAL X(XNONOD),Y(XNONOD),Z(XNONOD) 
      REAL T(NONODS),TOLD(NONODS)
      REAL DIFXT(NONODS), DIFYT(NONODS), DIFZT(NONODS)
      REAL DIFXTOLD(NONODS), DIFYTOLD(NONODS), DIFZTOLD(NONODS)
      REAL ML(NONODS) 
      REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
      REAL WEIGHT(NGI) 
      LOGICAL D3,DCYL
      INTEGER NDGLNO(TOTELE*NLOC),XNDGLN(TOTELE*NLOC)
      REAL DETWEI(NGI),RA(NGI)
      REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
      
!     Local variables
      INTEGER ELE,ILOC,JLOC
      INTEGER GLOBI,GLOBJ,GI
      REAL VOLUME,NN,NNX,NNY,NNZ
!     
      DIFXT(1:NONODS) = 0.0
      DIFYT(1:NONODS) = 0.0
      DIFZT(1:NONODS) = 0.0
!     
      DIFXTOLD(1:NONODS) = 0.0
      DIFYTOLD(1:NONODS) = 0.0
      DIFZTOLD(1:NONODS) = 0.0
!     
      ML(1:NONODS) = 0.0
!     
      do  ELE=1,TOTELE! Was loop 340
!     
!     Calculate DETWEI,RA,NX,NY,NZ for element ELE
         CALL DETNLXR(ELE, X,Y,Z, XNDGLN, TOTELE,NONODS,NLOC,NGI, &
     &        N,NLX,NLY,NLZ, WEIGHT, DETWEI,RA,VOLUME, D3,DCYL, &
     &        NX,NY,NZ) 
!     
      do  ILOC=1,NLOC! Was loop 250
          GLOBI =NDGLNO((ELE-1)*NLOC+ILOC)
      do  JLOC=1,NLOC! Was loop 260
               GLOBJ =NDGLNO((ELE-1)*NLOC+JLOC)
               NN=0.
               NNX=0. 
               NNY=0.
               NNZ=0.
      do  GI=1,NGI! Was loop 331
                  NN =NN +N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
!     
                  NNX=NNX+N(ILOC,GI)*NX(JLOC,GI)*DETWEI(GI)
                  NNY=NNY+N(ILOC,GI)*NY(JLOC,GI)*DETWEI(GI)
                  NNZ=NNZ+N(ILOC,GI)*NZ(JLOC,GI)*DETWEI(GI)
      end do ! Was loop 331

               DIFXT(GLOBI)=DIFXT(GLOBI)+NNX*T(GLOBJ) 
               DIFYT(GLOBI)=DIFYT(GLOBI)+NNY*T(GLOBJ)
               DIFZT(GLOBI)=DIFZT(GLOBI)+NNZ*T(GLOBJ)

               DIFXTOLD(GLOBI)=DIFXTOLD(GLOBI)+NNX*TOLD(GLOBJ) 
               DIFYTOLD(GLOBI)=DIFYTOLD(GLOBI)+NNY*TOLD(GLOBJ)
               DIFZTOLD(GLOBI)=DIFZTOLD(GLOBI)+NNZ*TOLD(GLOBJ)
!     
               ML(GLOBI)=ML(GLOBI)+NN
      end do ! Was loop 260
      end do ! Was loop 250
!     
      end do ! Was loop 340
!     
      do GLOBI=1,NONODS! Was loop 
         DIFXT(GLOBI)=DIFXT(GLOBI) /ML(GLOBI)
         DIFYT(GLOBI)=DIFYT(GLOBI) /ML(GLOBI)
         DIFZT(GLOBI)=DIFZT(GLOBI) /ML(GLOBI)
!     
         DIFXTOLD(GLOBI)=DIFXTOLD(GLOBI) /ML(GLOBI) 
         DIFYTOLD(GLOBI)=DIFYTOLD(GLOBI) /ML(GLOBI)
         DIFZTOLD(GLOBI)=DIFZTOLD(GLOBI) /ML(GLOBI)
      END DO
      RETURN

  end subroutine getgxyz
  !
  subroutine INVXXX( &
       AGI,BGI,CGI, DGI,EGI,FGI, GGI,HGI,KGI, &
       A11,A12,A13, A21,A22,A23, A31,A32,A33)
    ! This sub finds the inverse of a 3x3 matrix
    REAL AGI,BGI,CGI, DGI,EGI,FGI, GGI,HGI,KGI
    REAL A11,A12,A13, A21,A22,A23, A31,A32,A33
    ! Local variables...
    REAL DETJ
    !
    DETJ=AGI*(EGI*KGI-FGI*HGI) &
         -BGI*(DGI*KGI-FGI*GGI) &
         +CGI*(DGI*HGI-EGI*GGI)
    !
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
    RETURN
  END subroutine INVXXX
  !
  !
  !
  !
  SUBROUTINE PHACMC(NDPSET,D3,CMC,CMCP,C1T,C2T,C3T, &
       C1TP,C2TP,C3TP,&
       FRAPHA,&
       ML1PHA,ML2PHA,ML3PHA,DM1PHA,DM2PHA,DM3PHA,&
       MPDIAG,NMPDIA,&
       NONODS,FREDOP,NCT,NCMC,&
                                ! new line...
       NCMCB,NPRESS,NPROPT,RADISO,ALFST2,SPRESS,KETPRE,&
       MODELB,&
       FINCMC,COLCMC,MIDCMC,&
       FINDCT,COLCT, &
                                ! RMRTIN is only used if ROTAT it has length =3*NONODS in2-D 
                                ! N.B.  NRMRT=NDIM*NDIM*NPHASE*NPHASE*NONODS
                                ! and 6*NONODS in 3-D.  
       WORK,&
       RMRTIN,NRMRT,GETINV,     &
                                ! Th following is for rotations only ROTAT=.TRUE.
       ROTAT,NNODRO,  NODROT,&
       NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,&
                                ! Multi-phase stuff
       NPHASE,NDVOLF, MASSP,LPRES,COMPRE,DT )
    IMPLICIT NONE
    ! ALFST2= max volume fraction of solids. 
    ! SPRESS=2nd pressure scaling parameter.
    INTEGER NONODS,FREDOP
    INTEGER NPHASE,NDPSET,NPROPT,RADISO,NPRESS
    REAL ALFST2,SPRESS
    REAL KETPRE(FREDOP)
    ! KETPRE=2(1+e)*gan. temp.
    ! If MODELB then use model B without gas pressure term in solid momentum.
    LOGICAL MODELB
    INTEGER NCT,NCMC,NCMCB
    REAL FRAPHA(FREDOP*NPHASE)
    REAL ML1PHA(NONODS*NPHASE),ML2PHA(NONODS*NPHASE)
    REAL ML3PHA(NONODS*NPHASE)
    REAL DM1PHA(NONODS*NPHASE),DM2PHA(NONODS*NPHASE)
    REAL DM3PHA(NONODS*NPHASE)
    INTEGER NMPDIA
    REAL MPDIAG(NMPDIA)
    LOGICAL GETINV
    ! If ROTAT then it is assumed that C1T etc already conbtain the 
    ! rotation matrices via C1T=R^T*C1T
    ! This sub obtains the pressure matrices CMC & CMCP(which is non-sym). 
    REAL C1T(NCT*NPRESS),C2T(NCT*NPRESS),C3T(NCT*NPRESS)
    REAL CMC(NCMCB),CMCP(NCMCB)
    REAL C1TP(NCT*NPHASE),C2TP(NCT*NPHASE),C3TP(NCT*NPHASE)
    ! NDVOLF=volume frac of all phases. 
    ! Put volume fraction of this phase into RMEM(1) & make a no not 
    ! too close to 0.0 if not multi-phase then set volume frac=1.
    ! NB NMPDIA=NPHASE*NPHASE*NONODS*NDIM*NDIM.
    ! NMXBD=maximum value of NDIM*NPHASE
    ! IF GETINV then get the inverse of MPDIAG with b.c's in. 
    INTEGER NMXBD
    PARAMETER(NMXBD=10)
    REAL BDIAG(NMXBD,NMXBD),BDX(NMXBD),BDB(NMXBD)
    REAL BDX2(NMXBD),BDB2(NMXBD)
    REAL BDIA2(NMXBD,NMXBD)
    !
    REAL NDVOLF(NONODS*NPHASE)
    REAL MASSP(FREDOP),LPRES(NPHASE*FREDOP),DT
    LOGICAL COMPRE(NPHASE)
    INTEGER NRMRT
    REAL WORK(NONODS),RMRTIN(NRMRT)
    REAL AGI,BGI,CGI,DGI,EGI,FGI,GGI,HGI,KGI
    REAL DETJ
    INTEGER ROW,START,FINISH,ST2,FI2
    INTEGER FINDCT(FREDOP+1),COLCT(NCT)
    INTEGER FINCMC(FREDOP+1),COLCMC(NCMC)
    INTEGER MIDCMC(FREDOP)
    INTEGER ACR,COL,ST1,FI1,I,J
    LOGICAL D3,LD3
    ! NB The normal is pointing out of the domain. 
    !  The rotation matrix in 3-D is R=  
    !   T1X   T1Y   T1Z
    !   T2X   T2Y   T2Z
    !   NX    NY    NZ
    !  The rotation matrix in 2-D is R=  
    !   T1X   T1Y   
    !   NX    NY     
    ! 
    INTEGER NNODRO
    REAL NX(NNODRO),NY(NNODRO),NZ(NNODRO)
    REAL T1X(NNODRO),T1Y(NNODRO),T1Z(NNODRO)
    REAL T2X(NNODRO),T2Y(NNODRO),T2Z(NNODRO)
    INTEGER NODROT(NNODRO)
    LOGICAL ROTAT
    ! A COLOURIN TECHNIQUE WORK BE MORE EFFICIENT 
    ! HERE(SEE BOTTON OF SUB). 
    ! The following forms the matrix C1T RML C1 + C2T RML C2.
    ! Often there is only one diagonal so that ML=ML1,ML2,ML3. 
    ! If D3 then suppose the flow is 3-D else suppose it is 2-D. 
    ! May have to alter this to make it more efficient. 
    ! If a multi-phase solution then accumulate the 
    ! pressure matrix
    ! Local mem...
    INTEGER NDIM,II,NORD,IM,JM,IP,JP
    INTEGER IM1,IM2,IM3,NOD1,NOD2,NOD3
    INTEGER ICOL,ICOLPL,ICOUNT,IBADD,NOD
    INTEGER count,ict,IPLUS,icent
    REAL CTC,CTCP,PRELAX
    REAL CTC11,CTC12,CTC21,CTC22
    REAL CTCP11,CTCP12,CTCP21,CTCP22
    REAL RGAS
    LOGICAL PDIAG,TEST
    LOGICAL YES,PRESS2
    real rmin,rmax,rmin2,rmax2,rmin3,rmax3
    !          INTEGER INCTY1
    !          SAVE INCTY1
    !          DATA INCTY1 /0/
    !          INCTY1=INCTY1+1
    !

    !
    ewrite(3,*) 'inside phacmc'
    !       stop 266
    ! Do we have two pressures?
    PRESS2=(NCMCB.GT.NCMC)
    !
    !       TEST=.TRUE.
    TEST=.FALSE.
    !
    IF(TEST) THEN
       ! Replace C1TP with C1T etc. 
       do IP=1,NPHASE
          ICT=(IP-1)*NCT
          IPLUS=(IP-1)*FREDOP
          do I=1,FREDOP
             do COUNT=FINDCT(I),FINDCT(I+1)-1
                !                C1TP(COUNT+ICT)=C1T(COUNT)
                !                C2TP(COUNT+ICT)=C2T(COUNT)
                C1TP(COUNT+ICT)=FRAPHA(I+IPLUS)*C1T(COUNT)
                C2TP(COUNT+ICT)=FRAPHA(I+IPLUS)*C2T(COUNT)
             END DO
          END DO
       END DO
    ENDIF
    !
    !
    ewrite(1,*)'JUST INSIDE PHACMC'  
    !        ewrite(3,*)'ML1=',ML1 
    !        ewrite(3,*)'DM1=',DM1
    !
    NDIM=2
    IF(D3) NDIM=3
    NORD=NDIM*NPHASE
    IF(NRMRT.LT.NONODS*NORD*NORD) THEN
       ewrite(3,*) 'NOT ENOUGH MEM FOR NRMRT'
       FLAbort("STOP")
    ENDIF
    !
    !         ewrite(3,*) 'NMPDIA,nonods,nphase:',NMPDIA,nonods,nphase
    do nod=1,-nonods
       do im=1,nord
          ewrite(3,*) 'i,mpdiag:',nod,&
               (mpdiag(nod+((im-1)*nord+jm-1)*nonods),jm=1,nord)
          !
       end do
    end do
    !
    !
    ! NB in 2-D it has block numbering ...
    ! 1 *
    ! 2 3
    ! and in 3-D ...
    ! 1 * *
    ! 2 3 *
    ! 4 5 6       *-not stored
    !
    ! Form M^{-1}
    ewrite(3,*) 'BEFORE PHAINV GETINV=',GETINV
    IF(GETINV) THEN
       ewrite(1,*)'entering PHAINV' 
       CALL PHAINV(D3,&
            DM1PHA,DM2PHA,DM3PHA,&
            MPDIAG,NMPDIA,&
            NONODS,FREDOP,  CMC,NCMC,&
                                ! RMRTIN is only used if ROTAT it has length =3*NONODS in2-D 
                                ! N.B.  NRMRT=NDIM*NDIM*NPHASE*NPHASE*NONODS
                                ! and 6*NONODS in 3-D.  
            WORK,&
            RMRTIN,NRMRT,   &
                                ! More work space..
                                !     &         NMXBD,BDIAG,BDIA2,BDX,BDB, 
            NMXBD,BDIAG,BDX,BDB,    &
                                ! Th following is for rotations only ROTAT=.TRUE.
            ROTAT,NNODRO,  NODROT,&
            NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,&
                                ! Multi-phase stuff
            NPHASE)
       ewrite(3,*)'out of PHAINV' 
    ENDIF
    !
    ewrite(3,*) 'the inver in getcmc:'
    do nod=1,-nonods
       do im=1,nord
          ewrite(3,*) 'i,RMRTIN:',nod,&
               (RMRTIN(nod+((im-1)*nord+jm-1)*nonods),jm=1,nord)
          !
       end do
    end do
    !
    !
    ! NOW FORM CMC & CMCP
    !
    !       CALL RCLEAR(CMC, NCMC)
    !       CALL RCLEAR(CMCP,NCMC)
    !
    do  ROW=1,FREDOP! Was loop 303
       !         ewrite(3,*)'row,ml1(row):',row,ml1(row)
       START=FINCMC(ROW)
       FINISH=FINCMC(ROW+1)-1
       ST1=FINDCT(ROW)
       FI1=FINDCT(ROW+1)-1
       do  ACR=START,FINISH! Was loop 313
          COL=COLCMC(ACR)
          CTC11 =0.
          CTC12 =0.
          CTC21 =0.
          CTC22 =0.
          !
          CTCP11=0.
          CTCP12=0.
          CTCP21=0.
          CTCP22=0.
          ! We wish to mult vector column 'ROW' of C1 by vector
          ! column 'COL' of C1.
          ST2=FINDCT(COL)
          FI2=FINDCT(COL+1)-1
          do  I=ST1,FI1! Was loop 323
             do  J=ST2,FI2! Was loop 333
                IF(COLCT(I).EQ.COLCT(J)) THEN
                   ICOL=COLCT(I)
                   !
                   ! NB RMRTIN contains (R M R^T +D)^{-1}
                   ! Work out BDIAG the inverse at node ICOL.
                   do IM=1,NORD
                      do JM=1,NORD
                         BDIAG(IM,JM)=RMRTIN(ICOL+((IM-1)*NORD+JM-1)*NONODS)
                      END DO
                      !               ewrite(3,*) 'icol,im,bdiag:',icol,im,BDIAG(IM,iM)
                   END DO
                   !
                   do IP=1,NPHASE
                      IM1=(IP-1)*NDIM+1
                      IM2=(IP-1)*NDIM+2
                      IM3=(IP-1)*NDIM+3
                      NOD1=((IM1-1)*NORD+IM1-1)*NONODS  + ICOL
                      NOD2=((IM2-1)*NORD+IM2-1)*NONODS  + ICOL
                      NOD3=((IM3-1)*NORD+IM3-1)*NONODS  + ICOL
                      ICOLPL=ICOL+(IP-1)*NONODS
                      !
                      CTC11=CTC11+C1T(I)*C1T(J)*NDVOLF(ICOLPL)&
                           /( max(ML1PHA(ICOLPL),MPDIAG(NOD1)) +DM1PHA(ICOLPL))&
                                !
                           +C2T(I)*C2T(J)*NDVOLF(ICOLPL)&
                           /( max(ML2PHA(ICOLPL),MPDIAG(NOD2)) +DM2PHA(ICOLPL))
                      !
                      IF(D3) &
                           CTC11=CTC11+C3T(I)*C3T(J)*NDVOLF(ICOLPL)&
                           /( max(ML3PHA(ICOLPL),MPDIAG(NOD3)) +DM3PHA(ICOLPL))
                      !
                      BDB((IP-1)*NDIM+1)=C1T(J)
                      BDB((IP-1)*NDIM+2)=C2T(J)
                      IF(D3)BDB((IP-1)*NDIM+3)=C3T(J)
                      IF(MODELB) THEN
                         IF(IP.EQ.2) THEN
                            BDB((IP-1)*NDIM+1)=0.
                            BDB((IP-1)*NDIM+2)=0.
                            IF(D3)BDB((IP-1)*NDIM+3)=0.
                         ENDIF
                      ENDIF
                   END DO
                   !
                   IF(PRESS2) THEN
                      IP=1
                      BDB2((IP-1)*NDIM+1)=0.0
                      BDB2((IP-1)*NDIM+2)=0.0
                      IF(D3)BDB2((IP-1)*NDIM+3)=0.0
                      IP=2
                      BDB2((IP-1)*NDIM+1)=C1T(J+NCT)
                      BDB2((IP-1)*NDIM+2)=C2T(J+NCT)
                      IF(D3)BDB2((IP-1)*NDIM+3)=C3T(J+NCT)
                      !                 BDB2((IP-1)*NDIM+1)=C1TP(J+NCT)+C1T(J+NCT)
                      !                 BDB2((IP-1)*NDIM+2)=C2TP(J+NCT)+C2T(J+NCT)
                      !           IF(D3)BDB2((IP-1)*NDIM+3)=C3TP(J+NCT)+C3T(J+NCT)
                      IM1=(IP-1)*NDIM+1
                      IM2=(IP-1)*NDIM+2
                      IM3=(IP-1)*NDIM+3
                      NOD1=((IM1-1)*NORD+IM1-1)*NONODS  + ICOL
                      NOD2=((IM2-1)*NORD+IM2-1)*NONODS  + ICOL
                      NOD3=((IM3-1)*NORD+IM3-1)*NONODS  + ICOL
                      ICOLPL=ICOL+(IP-1)*NONODS
                      !
                      CTC12=CTC12+C1T(I)*C1T(J)*NDVOLF(ICOLPL)&
                           /( max(ML1PHA(ICOLPL),MPDIAG(NOD1)) +DM1PHA(ICOLPL))&
                                !
                           +C2T(I)*C2T(J)*NDVOLF(ICOLPL)&
                           /( max(ML2PHA(ICOLPL),MPDIAG(NOD2)) +DM2PHA(ICOLPL))
                      !
                      IF(D3) &
                           CTC12=CTC12+C3T(I)*C3T(J)*NDVOLF(ICOLPL)&
                           /( max(ML3PHA(ICOLPL),MPDIAG(NOD3)) +DM3PHA(ICOLPL))
                      !
                      ! ENDOF IF(PRESS2) THEN...
                   ENDIF
                   !               
                   !
                   IF(PRESS2) THEN
                      do IM=1,NORD
                         BDX(IM)=0.
                         BDX2(IM)=0.
                         do JM=1,NORD
                            BDX(IM) =BDX(IM) +BDIAG(IM,JM)*BDB(JM)
                            BDX2(IM)=BDX2(IM)+BDIAG(IM,JM)*BDB2(JM)
                         END DO
                      END DO
                   ELSE
                      do IM=1,NORD
                         BDX(IM)=0.
                         do JM=1,NORD
                            BDX(IM) =BDX(IM) +BDIAG(IM,JM)*BDB(JM)
                         END DO
                      END DO
                   ENDIF
                   !
                   do IP=1,NPHASE
                      BDB((IP-1)*NDIM+1)=C1TP(I+(IP-1)*NCT)
                      BDB((IP-1)*NDIM+2)=C2TP(I+(IP-1)*NCT)
                      IF(D3)BDB((IP-1)*NDIM+3)=C3TP(I+(IP-1)*NCT)
                   END DO
                   IF(PRESS2) THEN
                      IP=1
                      BDB2((IP-1)*NDIM+1)=0.
                      BDB2((IP-1)*NDIM+2)=0.
                      IF(D3)BDB2((IP-1)*NDIM+3)=0.
                      IP=2
                      BDB2((IP-1)*NDIM+1)=C1TP(I+(IP-1)*NCT)
                      BDB2((IP-1)*NDIM+2)=C2TP(I+(IP-1)*NCT)
                      IF(D3)BDB2((IP-1)*NDIM+3)=C3TP(I+(IP-1)*NCT)
                      !        ENDIF
                      !
                      !         IF(PRESS2) THEN
                      do IM=1,NORD
                         CTCP11=CTCP11+BDB(IM) *BDX(IM)
                         CTCP12=CTCP12+BDB(IM) *BDX2(IM)
                         CTCP21=CTCP21+BDB2(IM)*BDX(IM)
                         CTCP22=CTCP22+BDB2(IM)*BDX2(IM)
                      END DO
                   ELSE
                      do IM=1,NORD
                         CTCP11=CTCP11+BDB(IM)*BDX(IM)
                      END DO
                   ENDIF
                   !
                   !
                ENDIF
             end do ! Was loop 333
          end do ! Was loop 323
          CMCP(ACR)       =CTCP11
          CMC(ACR)        =CTC11
          IF(PRESS2) THEN
             CMCP(ACR+  NCMC)=CTCP12
             CMCP(ACR+2*NCMC)=CTCP21
             CMCP(ACR+3*NCMC)=CTCP22 
             !     CMC(ACR+  NCMC) =CTC12
             !     CMC(ACR+2*NCMC) =CTC12
             CMC(ACR+  NCMC) =0.
             CMC(ACR+2*NCMC) =0.
             CMC(ACR+3*NCMC) =CTC12   
          ENDIF
       end do ! Was loop 313
    end do ! Was loop 303
    !
    do IP=1,NPHASE
       rmax=-1.e+20
       rmin=1.e+20
       rmax2=-1.e+20
       rmin2=1.e+20
       rmax3=-1.e+20
       rmin3=1.e+20
       IF(COMPRE(IP)) THEN
          do I=1,FREDOP
             ICOUNT=MIDCMC(I)
             !            rmax2=max(rmax2,CMC(ICOUNT))
             !            rmin2=min(rmin2,CMC(ICOUNT))
             !            rmax3=max(rmax3,CMCp(ICOUNT))
             !            rmin3=min(rmin3,CMCp(ICOUNT))
             CMC(ICOUNT) =CMC(ICOUNT)  &
                  + MASSP(I)*LPRES(I+(IP-1)*FREDOP)/DT
             CMCP(ICOUNT)=CMCP(ICOUNT) &
                  + MASSP(I)*LPRES(I+(IP-1)*FREDOP)/DT
             !            rmax=max(rmax,MASSP(I)*LPRES(I+(IP-1)*FREDOP)/DT)
             !            rmin=min(rmin,MASSP(I)*LPRES(I+(IP-1)*FREDOP)/DT)
          END DO
       ENDIF
       !
       IF(PRESS2.AND.(IP.EQ.2)) THEN
          ewrite(3,*) 'PUTTING DIAGONAL IN FOR 2ND PRESSUR NPROPT',NPROPT
          ewrite(3,*) 'ALFST2=',ALFST2
          ! Add diagonal term in...
          do I=1,FREDOP
             !            ewrite(3,*) 'i,FRAPHA(FREDOP+I):',i,FRAPHA(FREDOP+I)
             ICOUNT=MIDCMC(I)
             RGAS=MAX(0.0,FRAPHA(FREDOP+I)) 
             !
             if(NPROPT.GT.0) then
                ! CONSISTENTLY DISCRETISED METHOD...
                CMC(ICOUNT+3*NCMC) =CMC(ICOUNT +3*NCMC)  &
                     + max(0.0,FUNADP(NPROPT,RADISO,RGAS,ALFST2,SPRESS,KETPRE(I)))&
                     *MASSP(I)/DT**2
                !
                CMCP(ICOUNT+3*NCMC)=CMCP(ICOUNT+3*NCMC) &
                     + max(0.0,FUNADP(NPROPT,RADISO,RGAS,ALFST2,SPRESS,KETPRE(I)))&
                     *MASSP(I)/DT**2
             else
                ! Pressure wave speed method...
                CMC(ICOUNT+3*NCMC) =CMC(ICOUNT +3*NCMC)  &
                     + max(0.0,DEDADP(NPROPT,RGAS,ALFST2,SPRESS))&
                     *MASSP(I)/DT**2
                !
                CMCP(ICOUNT+3*NCMC)=CMCP(ICOUNT+3*NCMC) &
                     + max(0.0,DEDADP(NPROPT,RGAS,ALFST2,SPRESS))&
                     *MASSP(I)/DT**2
                !
             endif
             !
          END DO
          ewrite(3,*) 'FINISHED PUTTING DIAGONAL IN FOR 2ND PRESSUR'
       END IF
       !
       ! Put the compressible mass matrix into the matrices CMC,CMCP. 
       !        ewrite(3,*) 'COMPRE:',COMPRE
       !        I=FREDOP/2
       !        ewrite(3,*) 'I=',I
       !        DO COUNT=FINCMC(I),FINCMC(I+1)-1
       !          ewrite(3,*) 'COLCMC(COUNT),CMC(COUNT),CMCP(COUNT):',
       !     &             COLCMC(COUNT),CMC(COUNT),CMCP(COUNT)
       !          ewrite(3,*) 'CMC(COUNT+  NCMC),CMCP(COUNT+  NCMC):',
       !     &             CMC(COUNT+  NCMC),CMCP(COUNT+  NCMC)
       !          ewrite(3,*) 'CMC(COUNT+2*NCMC),CMCP(COUNT+2*NCMC):',
       !     &             CMC(COUNT+2*NCMC),CMCP(COUNT+2*NCMC)
       !          ewrite(3,*) 'CMC(COUNT+3*NCMC),CMCP(COUNT+3*NCMC):',
       !     &             CMC(COUNT+3*NCMC),CMCP(COUNT+3*NCMC)
       !        END DO
       !        STOP 2556
       !         ewrite(3,*) 'absorption addition to getcmc rmin,rmax=',
       !     &            rmin,rmax
       !         ewrite(3,*) 'original rmin2,rmax2:',rmin2,rmax2
       !         ewrite(3,*) 'original rmin3,rmax3:',rmin3,rmax3
       !         i=1340
       !         ewrite(3,*) 'i,CMC(MIDCMC(i)),CMCp(MIDCMC(i)):',
       !     &            i,CMC(MIDCMC(i)),CMCp(MIDCMC(i))
       !         ewrite(3,*) 'MASSP(I)*LPRES(I+(IP-1)*FREDOP)/DT',
       !     &            MASSP(I)*LPRES(I+(IP-1)*FREDOP)/DT
       !         i=1360
       !         ewrite(3,*) 'i,CMC(MIDCMC(i)),CMCp(MIDCMC(i)):',
       !     &            i,CMC(MIDCMC(i)),CMCp(MIDCMC(i))
       !         ewrite(3,*) 'MASSP(I)*LPRES(I+(IP-1)*FREDOP)/DT',
       !     &            MASSP(I)*LPRES(I+(IP-1)*FREDOP)/DT
    END DO
    !         stop
    !
    IF((NDPSET.GE.1).AND.(NDPSET.LE.FREDOP)) THEN 
       ! SET PRSEEURE NODE NDPSET TO ZERO.
       !        ewrite(3,*)'CMC(MIDCMC(NDPSET))=',CMC(MIDCMC(NDPSET))
       CMC(MIDCMC(NDPSET)) =INFINITY
       CMCP(MIDCMC(NDPSET))=INFINITY
    ENDIF
    !            stop 3636
    !
    !
    IF(ROTAT) THEN
       !         IF(.false.) THEN
       CMC(1:NCMCB) = CMCP(1:NCMCB)
       do I=1,-FREDOP
          ICENT=MIDCMC(I)
          CMC(ICENT)=1.1*CMC(ICENT)
          IF(NPRESS.GT.1) CMC(ICENT+3*NCMC)=1.1*CMC(ICENT+3*NCMC)
       END DO
       !           IF(NPRESS.GT.1) CALL RCLEAR(CMC(1+NCMC),2*NCMC)
    ENDIF
    !
    ! NB RMRTIN contains (R M R^T +D)^{-1}
    ! Work out BDIAG the inverse at node ICOL.
    ewrite(3,*) 'INVERSE OF MPDIAG:'
    do ICOL=1,-NONODS
       ewrite(3,*) 'NOD=',ICOL
       do IM=1,NORD
          do JM=1,NORD
             BDIAG(IM,JM)=RMRTIN(ICOL+((IM-1)*NORD+JM-1)*NONODS)
          END DO
          ewrite(3,*) (BDIAG(IM,JM),JM=1,NORD)
       END DO
    END DO
    !
    !          ewrite(3,*) 'IN PHACMC LPRES:'
    !          ewrite(3,*) LPRES
    !
    do  I=1,-FREDOP! Was loop 314
       !            ewrite(3,*)'I,COLM:',I,(COLCMC(ICOUNT),
       !     &         ICOUNT=FINCMC(I),FINCMC(I+1)-1)
       ewrite(3,*) 'i,LPRES(I),LPRES(I+FREDOP):',&
            i,LPRES(I),LPRES(I+FREDOP)
       !           ewrite(3,*)'cmc(midcmc(i)+3*ncmc)',
       !     &      cmc(midcmc(i)+3*ncmc),cmcp(midcmc(i)+3*ncmc)
       !            ewrite(3,*)'I,CMC:',I,(CMC(ICOUNT),
       !     &         ICOUNT=FINCMC(I),FINCMC(I+1)-1)
    end do ! Was loop 314
    !          stop 8689

    !
    !
    ! THE COLOURING ALGORITHM IS AS FOLLOWS. ...
    !        DO 10 COLOR=1,MXCOLO
    !             DO 20 I=1,FREDOP
    !                 V(I)=0.
    !                 IF(NODCOL(I).EQ.COLOR) V(I)=1.
    !20          CONTINUE
    !C PERFORM MATRIX VECTOR MULTIPLICATION 
    !C IE   x= CT ML C v 
    !             DO 30 I=1,FREDOP
    !             DO 30 COUNT=FINCMC(I),FINCMC(I+1)-1
    !                 IF(NODCOL(COLCMC(COUNT).EQ.COLOR) CMC(COUNT)=X(I)
    !30          CONTINUE
    !10     CONTINUE
    !       ewrite(3,*)'ML1=',ML1
    ewrite(1,*)'exiting getcmc'
    RETURN
  END SUBROUTINE PHACMC
  !
  !

  !
  ! 
  SUBROUTINE MATINV(A,N,NMAX,MAT,MAT2,X,B)
    ! This sub finds the inverse of the matrix A and puts it back in A. 
    ! MAT, MAT2 & X,B are working vectors. 
    IMPLICIT NONE
    INTEGER N,NMAX
    REAL A(NMAX,NMAX),MAT(N,N),MAT2(N,N),X(N),B(N)
    ! Local variables
    INTEGER ICOL,IM,JM

    if(.true.) then
       ! Solve MAT XL=BX (NB BDIAG is overwritten)
       ICOL=1
       B=0.0
       B(ICOL)=1.0
       do IM=1,N
          do JM=1,N
             MAT(IM,JM)=A(IM,JM)
          END DO
       END DO
       CALL SMLINNGOT(MAT,X,B,N,N,.FALSE.)
       ! X contains the column ICOL of inverse
       do IM=1,N
          MAT2(IM,ICOL)=X(IM)
       END DO
       do ICOL=2,N
          ! Form column ICOL of the inverse. 
          B=0.
          B(ICOL)=1.0
          ! Solve MAT X=B (NB MAT is overwritten).  
          CALL SMLINNGOT(MAT,X,B,N,N,.TRUE.)
          ! X contains the column ICOL of inverse
          do IM=1,N
             MAT2(IM,ICOL)=X(IM)
          END DO
       END DO
    else
       !
       do ICOL=1,N
          !
          ! Form column ICOL of the inverse. 
          do IM=1,N
             B(IM)=0.
             do JM=1,N
                MAT(IM,JM)=A(IM,JM)
             END DO
          END DO
          B(ICOL)=1.0
          ! Solve MAT X=B (NB MAT is overwritten).  
          CALL SMLINN(MAT,X,B,N,N)
          ! X contains the column ICOL of inverse
          do IM=1,N
             MAT2(IM,ICOL)=X(IM)
          END DO
          !
       END DO
    endif
    !
    ! Set A to MAT2
    do IM=1,N
       do JM=1,N
          A(IM,JM)=MAT2(IM,JM)
       END DO
    END DO
    RETURN
  END SUBROUTINE MATINV
  !
  !
  !
  !
  SUBROUTINE PHAINV(D3,&
       DM1PHA,DM2PHA,DM3PHA,&
       MPDIAG,NMPDIA,&
       NONODS,FREDOP,  TEMP,NTEMP,&
                                ! RMRTIN is only used if ROTAT it has length =3*NONODS in2-D 
                                ! N.B.  NRMRT=NDIM*NDIM*NPHASE*NPHASE*NONODS
                                ! and 6*NONODS in 3-D.  
       WORK,&
       RMRTIN,NRMRT,   &
                                ! More work space..
                                !     &         NMXBD,BDIAG,BDIA2,BDX,BDB, 
       NMXBD,BDIAG,BDX,BDB,    &
                                ! Th following is for rotations only ROTAT=.TRUE.
       ROTAT,NNODRO,  NODROT,&
       NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,&
                                ! Multi-phase stuff
       NPHASE)
    ! Form M^{-1} for multi-phase flow.
    IMPLICIT NONE
    INTEGER NONODS,FREDOP
    ! If ROTAT then it is assumed that C1T etc already conbtain the 
    ! rotation matrices via C1T=R^T*C1T
    INTEGER NTEMP
    REAL TEMP(NTEMP)
    ! NDVOLF=volume frac of all phases. 
    ! Put volume fraction of this phase into RMEM(1) & make a no not 
    ! too close to 0.0 if not multi-phase then set volume frac=1.
    ! NB NMPDIA=NPHASE*NPHASE*NONODS*NDIM*NDIM.
    INTEGER NPHASE
    ! NMXBD=maximum value of NDIM*NPHASE
    !        PARAMETER(NMXBD=10)
    INTEGER NMXBD
    REAL BDIAG(NMXBD,NMXBD),BDX(NMXBD),BDB(NMXBD)
    !        REAL BDIA2(NMXBD,NMXBD)
    REAL BDIA2(10,10)
    REAL TEMPM1(10,10),TEMPM2(10,10),TEMPV1(10),TEMPV2(10)
    !
    REAL DM1PHA(NONODS*NPHASE),DM2PHA(NONODS*NPHASE)
    REAL DM3PHA(NONODS*NPHASE)
    INTEGER NMPDIA,NRMRT
    REAL MPDIAG(NMPDIA)
    !       REAL MPDIA2(NMPDIA)
    REAL WORK(NONODS),RMRTIN(NRMRT)
    REAL AGI,BGI,CGI,DGI,EGI,FGI,GGI,HGI,KGI
    REAL DETJ
    REAL RML1,RML2,RML3
    INTEGER ROW,START,FINISH,ST2,FI2
    INTEGER ACR,COL,ST1,FI1,I,J
    INTEGER IM,JM,IP,JP,ID,JD
    INTEGER NOD,II
    INTEGER IPL,JPL,IBADD,NPLUS
    LOGICAL D3
    ! NB The normal is pointing out of the domain. 
    !  The rotation matrix in 3-D is R=  
    !   T1X   T1Y   T1Z
    !   T2X   T2Y   T2Z
    !   NX    NY    NZ
    !  The rotation matrix in 2-D is R=  
    !   T1X   T1Y   
    !   NX    NY     
    ! 
    INTEGER NNODRO
    REAL NX(NNODRO),NY(NNODRO),NZ(NNODRO)
    REAL T1X(NNODRO),T1Y(NNODRO),T1Z(NNODRO)
    REAL T2X(NNODRO),T2Y(NNODRO),T2Z(NNODRO)
    INTEGER NODROT(NNODRO)
    LOGICAL ROTAT
    ! Local memory...
    INTEGER NDIM,NORD


    ! Form M^{-1} for multi-phase flow.

    ewrite(3,*)'inside PHAINV' 
    !
    NDIM=2
    IF(D3) NDIM=3
    NORD=NDIM*NPHASE
    !
    WORK(1:NONODS) = 0.0
    RMRTIN(1:NRMRT) = 0.0
    !
    IF(ROTAT) THEN
       do  II=1,NNODRO! Was loop 30
          NOD=NODROT(II)
          WORK(NOD)=1.0
          !
          !             DO IM=1,NORD
          !               DO JM=1,NORD
          !                 BDIAG(IM,JM)=0.
          !               END DO
          !             END DO
          !
          do IP=1,NPHASE
             do JP=1,NPHASE
                !
                IPL=(IP-1)*NDIM
                JPL=(JP-1)*NDIM
                !
                do ID=1,NDIM
                   do JD=1,NDIM
                      IM=(IP-1)*NDIM+ID
                      JM=(JP-1)*NDIM+JD
                      IBADD=((IM-1)*NORD+JM-1)*NONODS  + NOD
                      BDIAG(IPL+ID,JPL+JD)=MPDIAG(IBADD)
                   END DO
                END DO
                !
                ! Form R *(M_L R^T)
                IF(D3) THEN
                   ! form matrix BDIA2=(M_L R^T)
                   CALL SMLMA3(&
                        BDIA2(IPL+1,JPL+1),BDIA2(IPL+1,JPL+2),BDIA2(IPL+1,JPL+3),&
                        BDIA2(IPL+2,JPL+1),BDIA2(IPL+2,JPL+2),BDIA2(IPL+2,JPL+3),&
                        BDIA2(IPL+3,JPL+1),BDIA2(IPL+3,JPL+2),BDIA2(IPL+3,JPL+3),&
                                ! = (NB. this matrix contains M_l)...
                        BDIAG(IPL+1,JPL+1),BDIAG(IPL+1,JPL+2),BDIAG(IPL+1,JPL+3),&
                        BDIAG(IPL+2,JPL+1),BDIAG(IPL+2,JPL+2),BDIAG(IPL+2,JPL+3),&
                        BDIAG(IPL+3,JPL+1),BDIAG(IPL+3,JPL+2),BDIAG(IPL+3,JPL+3),&
                                ! *
                        T1X(II),   T2X(II),    NX(II),&
                        T1Y(II),   T2Y(II),    NY(II),&
                        T1Z(II),   T2Z(II),    NZ(II))
                   ! form matrix BDIAG=R*BDIA2
                   CALL SMLMA3(&
                        BDIAG(IPL+1,JPL+1),BDIAG(IPL+1,JPL+2),BDIAG(IPL+1,JPL+3),&
                        BDIAG(IPL+2,JPL+1),BDIAG(IPL+2,JPL+2),BDIAG(IPL+2,JPL+3),&
                        BDIAG(IPL+3,JPL+1),BDIAG(IPL+3,JPL+2),BDIAG(IPL+3,JPL+3),&
                                ! =
                        T1X(II),   T1Y(II),   T1Z(II),&
                        T2X(II),   T2Y(II),   T2Z(II),&
                        NX(II),    NY(II),    NZ(II),&
                                ! *
                        BDIA2(IPL+1,JPL+1),BDIA2(IPL+1,JPL+2),BDIA2(IPL+1,JPL+3),&
                        BDIA2(IPL+2,JPL+1),BDIA2(IPL+2,JPL+2),BDIA2(IPL+2,JPL+3),&
                        BDIA2(IPL+3,JPL+1),BDIA2(IPL+3,JPL+2),BDIA2(IPL+3,JPL+3))
                ELSE
                   ! form matrix BDIA2=(M_L R^T)
                   CALL SMLMA2(&
                        BDIA2(IPL+1,JPL+1),BDIA2(IPL+1,JPL+2),&
                        BDIA2(IPL+2,JPL+1),BDIA2(IPL+2,JPL+2),&
                                ! = (NB. this matrix contains M_l)...
                        BDIAG(IPL+1,JPL+1),BDIAG(IPL+1,JPL+2),&
                        BDIAG(IPL+2,JPL+1),BDIAG(IPL+2,JPL+2), &
                                ! *
                        T1X(II),    NX(II),&
                        T1Y(II),    NY(II))  
                   ! form matrix BDIAG=R*BDIA2
                   CALL SMLMA2(&
                        BDIAG(IPL+1,JPL+1),BDIAG(IPL+1,JPL+2),&
                        BDIAG(IPL+2,JPL+1),BDIAG(IPL+2,JPL+2),&
                                ! =
                        T1X (II),  T1Y(II),&
                        NX(II)  ,   NY(II),    &
                                ! *
                        BDIA2(IPL+1,JPL+1),BDIA2(IPL+1,JPL+2),&
                        BDIA2(IPL+2,JPL+1),BDIA2(IPL+2,JPL+2))
                END IF
                ! END OF DO JP=1,NPHASE
             END DO
             ! add in b.c's into BDIAG for phase IP. 
             NPLUS=NOD+(IP-1)*NONODS
             BDIAG(IPL+1,IPL+1)=BDIAG(IPL+1,IPL+1)+DM1PHA(NPLUS)
             BDIAG(IPL+2,IPL+2)=BDIAG(IPL+2,IPL+2)+DM2PHA(NPLUS)
             IF(D3)BDIAG(IPL+3,IPL+3)=BDIAG(IPL+3,IPL+3)+DM3PHA(NPLUS)
             ! END OF DO IP=1,NPHASE
          END DO
          !

          !
          ! add in b.c's into BDIAG
          ! Ivert BDIAG and put result in RMRTIN (WORK is work space)
          ! set values of (RMR^T +D)^{-1} put result backinto BDIAG. 
          CALL MATINV(BDIAG,NORD,NMXBD,TEMPM1,TEMPM2,&
               TEMPV1,TEMPV2)
          !
          do IM=1,NORD
             do JM=1,NORD
                RMRTIN( ((IM-1)*NORD+JM-1)*NONODS +NOD)=BDIAG(IM,JM)
             END DO
          END DO
          !
       end do ! Was loop 30
       ! ENDOF IF(ROTAT) THEN...
    ENDIF
    !
    ! Form the rest of the inverses *********
    !          ewrite(3,*) 'forming the rest of the invers'
    do NOD=1,NONODS
       IF(WORK(NOD).LT.0.5) THEN
          !
          do  IM=1,NORD! Was loop 145
             do  JM=1,NORD! Was loop 145
                IBADD=((IM-1)*NORD+JM-1)*NONODS  + NOD
                BDIAG(IM,JM)=MPDIAG(IBADD)
                !               IF(NOD.EQ.188) THEN
                !                 ewrite(3,*) 'IM,JM,BDIAG:',IM,JM,BDIAG(IM,JM)
                !               ENDIF
             end do ! Was loop 145
          end do ! Was loop 145
          !
          ! Add in the b.c's.  
          do IP=1,NPHASE
             IPL=(IP-1)*NDIM
             NPLUS=NOD+(IP-1)*NONODS
             BDIAG(IPL+1,IPL+1)=BDIAG(IPL+1,IPL+1)+DM1PHA(NPLUS)
             BDIAG(IPL+2,IPL+2)=BDIAG(IPL+2,IPL+2)+DM2PHA(NPLUS)
             IF(D3)BDIAG(IPL+3,IPL+3)=BDIAG(IPL+3,IPL+3)+DM3PHA(NPLUS)
             ! END OF DO IP=1,NPHASE
          END DO
          ! Ivert BDIAG and put result in RMRTIN (WORK is work space)
          ! set values of (RMR^T +D)^{-1} put result backinto BDIAG. 
          CALL MATINV(BDIAG,NORD,NMXBD,TEMPM1,TEMPM2,&
               TEMPV1,TEMPV2)
          !
          do IM=1,NORD
             do JM=1,NORD
                RMRTIN( ((IM-1)*NORD+JM-1)*NONODS +NOD)=BDIAG(IM,JM)
                !               IF(NOD.EQ.188) THEN
                !                 ewrite(3,*) 'IM,JM,ibvBDIAG:',IM,JM,BDIAG(IM,JM)
                !               ENDIF
             END DO
          END DO
          !
       ENDIF
    END DO
    ! Form the rest of the inverses ***-endof
    ewrite(3,*)'finishing PHAINV '
    !
    RETURN
  END SUBROUTINE PHAINV

  REAL FUNCTION DISFUNSUB(RGAS,ALFSTA,RADISO)
    ! This sub forms the radial distribution function for spheres.
    ! RGAS=volume fraction of spheres.
    use FLDebug

    IMPLICIT NONE
    REAL::RGAS,ALFSTA
    REAL,PARAMETER::THIRD=1./3.
    INTEGER::RADISO,RADIS2
    REAL::A,B,YD,Y,G0,RRGAS

    REAL::RGAS2,RGAS3

    RADIS2=RADISO+1

    IF(RADIS2 == 1)THEN ! (RADISO=0)
       RGAS3=0.62*RGAS/ALFSTA
       IF(RGAS3 < 0.6)THEN
          RGAS2 = RGAS3
       ELSE
          RGAS2 = 0.6 + 3.*(RGAS3-0.6) 
       ENDIF

       DISFUNSUB = DISFU9(RGAS2,0.62) 

    ELSEIF(RADIS2 == 2)THEN ! (RADISO=1)
       RGAS3=0.62*RGAS/ALFSTA

       if(rgas3 < 0.6) then
          RGAS2=rgas3
       else
          RGAS2=0.6 + 5.*(RGAS3-0.6) 
       endif

       DISFUNSUB= DISFU9(RGAS2,0.62) 

    ELSEIF(RADIS2 == 3)THEN ! (RADISO=2)
       RGAS3=0.62*RGAS/ALFSTA

       if(rgas3.lt.0.6) then
          RGAS2=rgas3
       else
          !          RGAS2=0.6 + 3.*(RGAS3-0.6) 
          RGAS2=0.6 + 10.*(RGAS3-0.6) 
       endif
       DISFUNSUB= DISFU9(RGAS2,0.62) 

    ELSEIF(RADIS2 == 4)THEN ! (RADISO=3)
       RGAS3=0.62*RGAS/ALFSTA

       RGAS2=rgas3

       DISFUNSUB= DISFU9(RGAS2,0.62) 
    ELSEIF(RADIS2 == 5)THEN ! (RADISO=4)
       RGAS3=0.62*RGAS/ALFSTA

       RGAS2=rgas3

       DISFUNSUB= DISFU9(RGAS2,0.62) +EXP(600.*(RGAS-ALFSTA))

    ELSEIF(RADIS2 == 6)THEN ! (RADISO=5) 
       RGAS3=0.62*RGAS/ALFSTA

       RGAS2=rgas3

       DISFUNSUB= DISFU9(RGAS2,0.62) +EXP(400.*(RGAS-ALFSTA))

    ELSEIF(RADIS2 == 7)THEN ! (RADISO=6) 
       RGAS3=0.62*RGAS/ALFSTA

       RGAS2=rgas3

       DISFUNSUB= DISFU9(RGAS2,0.62) +EXP(200.*(RGAS-ALFSTA)) 

    ELSEIF(RADIS2 == 8)THEN ! (RADISO=7) 
       RGAS3=0.62*RGAS/ALFSTA

       RGAS2=RGAS3

       DISFUNSUB= DISFU9(RGAS2,0.62) +EXP(100.*(RGAS-ALFSTA)) 

    ELSEIF(RADIS2 == 9)THEN ! (RADISO=8) 
       !
       DISFUNSUB= 1./(1.-(RGAS/0.62)**THIRD) +EXP(600.*(RGAS-0.58)) 

    ELSEIF(RADIS2 == 10)THEN ! (RADISO=9) 
       ! Lun and Savage (1986)
       DISFUNSUB=(1.-(RGAS/ALFSTA))**(-2.5*ALFSTA)

    ELSEIF(RADIS2 == 11)THEN ! (RADISO=10) 
       ! Gidaspow (1994) but amended to produce the min of unity.
       DISFUNSUB=max(1.0+3.*RGAS,0.6/(1.-(RGAS/ALFSTA)**THIRD))

    ELSEIF(RADIS2 == 12)THEN ! (RADISO=11) 
       ! Gidaspow (1994) but amended to produce the min of unity.
       DISFUNSUB=max(1.0+3.*RGAS,0.6/(1.-(RGAS/0.62)**THIRD))+&
            EXP(400.*(RGAS-0.58)) 

    ELSEIF(RADIS2 == 13)THEN ! (RADISO=12) 
       ! Gidaspow (1994) but amended to produce the min of unity.
       DISFUNSUB=max(1.0+3.*RGAS,0.6/(1.-(RGAS/ALFSTA)**THIRD))+&
            EXP(400.*(RGAS-0.58)) 

    ELSEIF(RADIS2 == 14)THEN ! (RADISO=13) 

       ! Carnahan and Starling (1969) ...
       DISFUNSUB=(1.-RGAS)**(-1) + 1.5*RGAS*(1.-RGAS)**(-2) +&
            0.5*(RGAS**2)*(1.-RGAS)**(-3) +&
            EXP(100.*(RGAS-0.55)) + &
            EXP(300.*(RGAS-0.58)) 

    ENDIF

    RETURN
  END FUNCTION DISFUNSUB

  REAL FUNCTION DISFU9(RGAS,ALFSTA)
    ! This sub forms the radial distribution function for spheres.
    ! RGAS=volume fraction of spheres.
    use FLDebug

    IMPLICIT NONE
    REAL::RGAS,ALFSTA
    REAL,PARAMETER::THIRD=1./3.,DELTAMOX=10.0
    REAL A,B,YD,Y,G0,RRGAS
    REAL,SAVE::G0HALF,G058
    DATA G0HALF /1.E+10/

    IF(RGAS.GT.0.5)THEN
       IF(G0HALF.GT.1000.0)G0HALF=&
            LOG10( 0.6/(1. - (0.5/ALFSTA)**THIRD) )

       ! for G0=40 at RGAS=0.62
       A=13.3333333  - 8.33333*G0HALF
       B=-6.6666666  + 5.16666*G0HALF

       IF(RGAS.GT.0.58) THEN
          ! Give it a steeper accent after 0.58 vol frac
          RRGAS=0.58
          YD=A*RRGAS+B
          G058=YD
          A=25.0 *1.7 - 25.0*G058
          B=-14.5*1.7 + 15.5*G058
       ENDIF

       YD=A*RGAS+B
       Y=10.0**YD
       DISFU9=Y
    ELSE
       G0=0.6/(1. - (RGAS/ALFSTA)**THIRD)
       DISFU9=G0
    ENDIF

    RETURN
  END FUNCTION DISFU9

end module hart3d_allsorts
