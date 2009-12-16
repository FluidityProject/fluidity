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

module mellor_yamada
  use quadrature
  use elements
  use allsorts
  use FLDebug
  use state_module
  use fields
  use equation_of_state
  use rotated_boundary_conditions_legacy
  use flcomms_module
  use sml
  use detnlxr_module
  use legacy_tensor_fields
  implicit none
  
  private
  public mellor_yamada_turbulence, mellor_yamada_state
  
contains

  SUBROUTINE mellor_yamada_turbulence(state, NONODS,SNONOD,&
       CURTQQ,CURTQQL,&
       TQQ,TQQL,PREVTQQ,PREVTQQL,NU,NV,NW,&
       X,Y,Z,D3,DCYL,&
       TEMP,&
       NTSUB,TSUB,NSALSUB,SALSUB, &
       NSUBNVLOC,NUSUB,NVSUB,NWSUB, &
       DENGAM,DENINI,TEMINI,GRAVTY,&
       COGRAX,BSOUX,BSOUY,BSOUZ,&
       VGRAVX,VGRAVY,VGRAVZ,&
       LES,&
       ! Salt...
       GOTSALT,salinity_option,DENGAM_sal,S0,TSALT, gottopbotdis, &
       ! This is the output...
       TMUXX,TMUXY,TMUXZ, &
       TMUYY,TMUYZ,TMUZZ, &
       SMUXX,SMUXY,SMUXZ, &
       SMUYY,SMUYZ,SMUZZ, &
       TMUXX1,TMUXY1,TMUXZ1, &
       TMUYY1,TMUYZ1,TMUZZ1, &
       TMUXX2,TMUXY2,TMUXZ2, &
       TMUYY2,TMUYZ2,TMUZZ2, &
       MUPTXX,MUPTXY,MUPTXZ,&
       MUPTYY,MUPTYZ,MUPTZZ,&
       TSOURC1,TABSOR1,&
       TSOURC2,TABSOR2,&
       ! for visualising viscosity and diffusivity...
       GOTVISFIELD,TVISVIS,TVISDIF, &
       ! Integer arrays and basis functions
       TOTELE,NLOC,NGI,NDGLNO,XONDGL,&
       N,NLX,NLY,NLZ, WEIGHT,&
       ! BCS...
       NOBCTQQ,BCT1QQ,BCT2QQ,&
       NOBCTLQQ,BCT1LQQ,BCT2LQQ,DT,&
       ! parallel stuff...
       nnodp,para,halo_tag)

    ! THIS SUB DOES MELLAR YAMADA TURBULENCE MODELLING
    ! it calculates all the relevant parameters.
    ! There are two eqns for q^2 (TQQ) and lq^2 (TLQQ)
    ! NU,NV,NW is the velocity.
    ! x,y,z are the coords.
    !
    ! TMUZZ, SMUZZ diffusivity of density and salinity.
    ! MUPTXX and similar components are the viscosity.
    ! TMUZZ1, TMUZZ2 diffusivity for the 2 turbulence eqns.
    ! TSOURC1 Source for  qq.
    ! TSOURC2 Source for  lqq.
    ! TABSOR1 Absorption for  qq.
    ! TABSOR2 Absorption for  lqq.

    type(state_type), intent(inout):: state
    REAL A1,A2,B1,B2,C1,E1,E2,E3
    ! these numbers have been checked again Hedongs code...
    PARAMETER(A1=0.92,B1=16.6,A2=0.74,B2=10.1)
    PARAMETER(C1=0.08,E1=1.8,E2=1.33,E3=1.0)
    LOGICAL HEDONG,HANERT,SIMP_DIF
    PARAMETER(HEDONG=.true.,HANERT=.TRUE.,SIMP_DIF=.false.) 
    ! If SIMP_DIF then use a very simple vertical diffusivity/viscocity. 
    ! min values of q^2 and q^2l
    REAL Q2MIN,Q2LMIN
!    PARAMETER(Q2MIN=5.0E-7,Q2LMIN=1.0E-5)
!    PARAMETER(Q2MIN=5.0E-7,Q2LMIN=1.0E-7)
    PARAMETER(Q2MIN=5.0E-7,Q2LMIN=1.0E-8)
    LOGICAL DEN_GRAD_DG
    PARAMETER(DEN_GRAD_DG=.false.) 
! If DEN_GRAD_DG take into account DG representation when forming 
! gradient of density.

    INTEGER NONODS,SNONOD,TOTELE,NLOC,NGI
    INTEGER NDGLNO(TOTELE*NLOC),XONDGL(TOTELE*NLOC)
    REAL TQQ(NONODS),TQQL(NONODS)
    REAL CURTQQ(NONODS),CURTQQL(NONODS)
    REAL PREVTQQ(NONODS),PREVTQQL(NONODS)
    REAL NU(NONODS),NV(NONODS),NW(NONODS)
    REAL X(NONODS),Y(NONODS),Z(NONODS)
    LOGICAL D3,DCYL
    REAL DENGAM,DENINI,TEMINI,GRAVTY
    LOGICAL COGRAX
    REAL BSOUX,BSOUY,BSOUZ
    REAL VGRAVX(NONODS),VGRAVY(NONODS),VGRAVZ(NONODS)
    LOGICAL LES
    ! IF LES then we are using an LES method.
    ! Salt...
    LOGICAL GOTSALT
    ! distance to top and bottom are available:
    logical gottopbotdis
    INTEGER salinity_option
    REAL DENGAM_sal,S0
    REAL TSALT(NONODS)
    REAL TEMP(NONODS)
    INTEGER NTSUB,NSALSUB
    REAL TSUB(TOTELE*NTSUB),SALSUB(TOTELE*NSALSUB)
      INTEGER NSUBNVLOC
      REAL NUSUB(TOTELE*NSUBNVLOC),NVSUB(TOTELE*NSUBNVLOC)
      REAL NWSUB(TOTELE*NSUBNVLOC)
    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)
    REAL TMUXX(NONODS),TMUXY(NONODS),TMUXZ(NONODS)
    REAL TMUYY(NONODS),TMUYZ(NONODS),TMUZZ(NONODS) 
    REAL SMUXX(NONODS),SMUXY(NONODS),SMUXZ(NONODS)
    REAL SMUYY(NONODS),SMUYZ(NONODS),SMUZZ(NONODS)
    REAL TMUXX1(NONODS),TMUXY1(NONODS),TMUXZ1(NONODS)
    REAL TMUYY1(NONODS),TMUYZ1(NONODS),TMUZZ1(NONODS)
    REAL TMUXX2(NONODS),TMUXY2(NONODS),TMUXZ2(NONODS)
    REAL TMUYY2(NONODS),TMUYZ2(NONODS),TMUZZ2(NONODS)

    REAL MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTXZ(SNONOD)
    REAL MUPTYY(SNONOD),MUPTYZ(SNONOD),MUPTZZ(SNONOD)
    REAL TSOURC1(NONODS),TABSOR1(NONODS)
    REAL TSOURC2(NONODS),TABSOR2(NONODS)
    LOGICAL GOTVISFIELD
    REAL TVISVIS(NONODS),TVISDIF(NONODS)
    INTEGER nnodp,para,halo_tag
    
    ! The distance to the top TOPDIS and bottom BOTDIS of the domain 
    ! in the vertical direction is used here.
    type(scalar_field), pointer:: topdis_field, botdis_field
    real, dimension(:), pointer:: topdis, botdis
! BC'S
    INTEGER NOBCTQQ,NOBCTLQQ
    INTEGER BCT2QQ(NOBCTQQ),BCT2LQQ(NOBCTLQQ)
    REAL BCT1QQ(NOBCTLQQ),BCT1LQQ(NOBCTLQQ)
! Local variables...
    REAL DETWEI(NGI),RA(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    REAL DT, topdis_igl
    real maxdz,mindz
! scalars...
    REAL KM,KH,QQGI,QQLGI,DENDZ,DUDZ,DVDZ,VOLUME,DEN
    REAL Q,L,GH,SH,SM,KQ,ZPOS,ZMIN,ZMAX,W
    REAL LWALLM1,K,FUTVAL,N2,LMAX2
    INTEGER ELE,GI,ILOC,JLOC,IGL,JGL,II,NOD,I,J
    REAL RT1,RT2,RT3,SALT
    REAL DDUDX,DDUDY,DDUDZ,DDVDX,DDVDY,DDVDZ
       REAL RTSUB,rTEMP,RSALSUB,rsalt,MAX_IN_ELE_DZ
       real min_difvis,max_difvis
       real min_l,min_q,min_km
       real max_l,max_q,max_km
       real min_qqgi,min_QQLGI
       real max_qqgi,max_QQLGI
! *********ALOCAT MEEMORY**************
    REAL, ALLOCATABLE, DIMENSION(:)::ML
    REAL, ALLOCATABLE, DIMENSION(:)::VDENDX
    REAL, ALLOCATABLE, DIMENSION(:)::VDENDY
    REAL, ALLOCATABLE, DIMENSION(:)::VDENDZ
    REAL, ALLOCATABLE, DIMENSION(:)::VEDUDX
    REAL, ALLOCATABLE, DIMENSION(:)::VEDUDY
    REAL, ALLOCATABLE, DIMENSION(:)::VEDUDZ
    REAL, ALLOCATABLE, DIMENSION(:)::VEDVDX
    REAL, ALLOCATABLE, DIMENSION(:)::VEDVDY
    REAL, ALLOCATABLE, DIMENSION(:)::VEDVDZ
    REAL, ALLOCATABLE, DIMENSION(:)::VEDWDX
    REAL, ALLOCATABLE, DIMENSION(:)::VEDWDY
    REAL, ALLOCATABLE, DIMENSION(:)::VEDWDZ
    REAL, ALLOCATABLE, DIMENSION(:)::DIFDEN
    REAL, ALLOCATABLE, DIMENSION(:)::DIFVIS
    REAL, ALLOCATABLE, DIMENSION(:)::DIFTUB
    REAL, ALLOCATABLE, DIMENSION(:)::WORK
    REAL, ALLOCATABLE, DIMENSION(:)::TDEN
       
    REAL, ALLOCATABLE, DIMENSION(:)::T1X
    REAL, ALLOCATABLE, DIMENSION(:)::T1Y
    REAL, ALLOCATABLE, DIMENSION(:)::T1Z
       
    REAL, ALLOCATABLE, DIMENSION(:)::T2X
    REAL, ALLOCATABLE, DIMENSION(:)::T2Y
    REAL, ALLOCATABLE, DIMENSION(:)::T2Z
    REAL, ALLOCATABLE, DIMENSION(:)::NORMX
    REAL, ALLOCATABLE, DIMENSION(:)::NORMY
    REAL, ALLOCATABLE, DIMENSION(:)::NORMZ
    REAL, ALLOCATABLE, DIMENSION(:)::VERT_DZ
       REAL, ALLOCATABLE, DIMENSION(:)::TDENSUB
! ALOCATE.......       
    ewrite(3,*) 'about to alocat'
    ALLOCATE(ML(NONODS))
    ewrite(3,*) 'after ml'
    ALLOCATE(VDENDX(NONODS))
    ALLOCATE(VDENDY(NONODS))
    ALLOCATE(VDENDZ(NONODS))
    ALLOCATE(VEDUDX(NONODS))
    ALLOCATE(VEDUDY(NONODS))
    ALLOCATE(VEDUDZ(NONODS))
    ewrite(3,*) 'before vedvdx'
    ALLOCATE(VEDVDX(NONODS))
    ALLOCATE(VEDVDY(NONODS))
    ALLOCATE(VEDVDZ(NONODS))
    ALLOCATE(VEDWDX(NONODS))
    ALLOCATE(VEDWDY(NONODS))
    ALLOCATE(VEDWDZ(NONODS))
    ALLOCATE(DIFDEN(NONODS))
    ALLOCATE(DIFVIS(NONODS))
    ALLOCATE(DIFTUB(NONODS))
    ALLOCATE(WORK(NONODS))
    ewrite(3,*) 'before tden'
    ALLOCATE(TDEN(NONODS))
    ALLOCATE(T1X(NONODS))
    ALLOCATE(T1Y(NONODS))
    ALLOCATE(T1Z(NONODS))
    ewrite(3,*) 'before t2x'
    ALLOCATE(T2X(NONODS))
    ALLOCATE(T2Y(NONODS))
    ALLOCATE(T2Z(NONODS)) 
    ewrite(3,*) 'before normx'
    ALLOCATE(NORMX(NONODS))
    ALLOCATE(NORMY(NONODS))
    ALLOCATE(NORMZ(NONODS))
    ALLOCATE(VERT_DZ(NONODS))
    
    ALLOCATE(TDENSUB(TOTELE*NLOC))
! *********END ALOCATE MEMORY**************
    ewrite(1,*) 'finished allocating here'
    
    if (gottopbotdis) then
      topdis_field => extract_scalar_field(state, "DistanceToTop")
      botdis_field => extract_scalar_field(state, "DistanceToBottom")
      topdis => topdis_field%val
      botdis => botdis_field%val
    else if (salinity_option==2) then
      ewrite(0, *) 'Salinity option 2 needs distance to top'
      ewrite(0, *) 'This is only available in combination with a free surface field (ident -29)'
      FLAbort("Sorry!")
    else if (.not. cograx) then
      ewrite(0, *) 'This coordinate option needs distance to top'
      ewrite(0, *) 'This is only available in combination with a free surface field (ident -29)'
      FLAbort("Sorry!")
    end if
    
    ZMIN=1.E+20
    ZMAX=-1.E+20
    do igl=1,nonods
       ZMIN=MIN(ZMIN,Z(IGL))
       ZMAX=MAX(ZMAX,Z(IGL))
    end do

    ML(1:NONODS) = 0.0
    VDENDX(1:NONODS) = 0.0
    VDENDY(1:NONODS) = 0.0
    VDENDZ(1:NONODS) = 0.0

    VEDUDX(1:NONODS) = 0.0
    VEDUDY(1:NONODS) = 0.0
    VEDUDZ(1:NONODS) = 0.0

    VEDVDX(1:NONODS) = 0.0
    VEDVDY(1:NONODS) = 0.0
    VEDVDZ(1:NONODS) = 0.0

    VEDWDX(1:NONODS) = 0.0
    VEDWDY(1:NONODS) = 0.0
    VEDWDZ(1:NONODS) = 0.0
    
    VERT_DZ(1:NONODS) = 0.0
    
    TDENSUB = 0.0

    ! calculate Rotation matrix R... 
    CALL GETROTGRAV(&
         T1X,T1Y,T1Z, &
         T2X,T2Y,T2Z, &
         NORMX,NORMY,NORMZ, &
         COGRAX,BSOUX,BSOUY,BSOUZ, &
         VGRAVX,VGRAVY,VGRAVZ, NONODS)

    ! The distance to the top TOPDIS and bottom BOTDIS of the domain 
    ! in the vertical direction is used here.

    ! Calculate the density
    DO IGL=1,NONODS
       SALT=0.0
       IF(GOTSALT) SALT=TSALT(IGL)

       if (gottopbotdis) then
         topdis_igl=topdis(igl)
       else
         topdis_igl=0.0
       end if
       call GETDEN_TEMP_SALT( &
            salinity_option,DENINI,DENGAM,TEMINI,DENGAM_sal,S0, &
            den,temp(igl),salt,topdis_igl)
       TDEN(IGL)=DEN
    END DO

        IF((NTSUB.NE.0).OR.(GOTSALT.AND.(NSALSUB.NE.0))) THEN
               DO ELE=1,TOTELE
                 DO ILOC=1,NLOC
                    I=(ELE-1)*NLOC+ILOC
                    NOD=NDGLNO((ELE-1)*NLOC+ILOC)
                    RTSUB=0.0
                    IF(NTSUB.NE.0) THEN
                      RTSUB=TSUB(I)
                    ENDIF
                    rTEMP=Temp(NOD)+RTSUB
                    
                    RSALSUB=0.0
                    IF(GOTSALT.AND.(NSALSUB.NE.0)) THEN
                      RSALSUB=SALSUB(I)
                    ENDIF
                    rsalt=0.0
                    if(gotsalt) rsalt=tsalt(NOD)+RSALSUB
                    
                    if(gottopbotdis) then
                      topdis_igl=topdis(nod)
                    else
                      topdis_igl=0.0
                    end if
                      
                    CALL GETDEN_TEMP_SALT( &
     &                   salinity_option,DENINI,DENGAM,TEMINI,DENGAM_sal &
     &                   ,S0,den,rtemp,rsalt,topdis_igl)
                    TDENSUB(I)=den  - TDEN(NOD)
                 ENDDO
               ENDDO
        ENDIF
 !!!!!       TDENSUB=0.0
    
    maxdz=-1.e+20
    mindz= 1.e+20
        
    element_loop: DO ELE=1,TOTELE

       ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
       CALL DETNLXR(ELE, X,Y,Z, XONDGL, TOTELE,NONODS,NLOC,NGI, &
            N,NLX,NLY,NLZ, WEIGHT, DETWEI,RA,VOLUME, D3,DCYL,&
            NX,NY,NZ) 

       DO ILOC=1,NLOC
          IGL=NDGLNO((ELE-1)*NLOC+ILOC)
          MAX_IN_ELE_DZ=0.0
          DO JLOC=1,NLOC
             JGL=NDGLNO((ELE-1)*NLOC+JLOC)
             J=(ELE-1)*NLOC+JLOC
             ! Assume we are solving for temp. TEMP(IGL)
             DO GI=1,NGI
                DEN=TDEN(JGL)

                VDENDX(IGL)=VDENDX(IGL)+N(ILOC,GI)*DETWEI(GI)*DEN*NX(JLOC,GI) 
                VDENDY(IGL)=VDENDY(IGL)+N(ILOC,GI)*DETWEI(GI)*DEN*NY(JLOC,GI) 
                VDENDZ(IGL)=VDENDZ(IGL)+N(ILOC,GI)*DETWEI(GI)*DEN*NZ(JLOC,GI)
                
! calc. lumped vertical length scale...                
               MAX_IN_ELE_DZ=max(MAX_IN_ELE_DZ, &
                ABS(NORMX(IGL)*NX(JLOC,GI)+NORMY(IGL)*NY(JLOC,GI)+NORMZ(IGL)*NZ(JLOC,GI)) )
 
               IF(DEN_GRAD_DG) THEN
! Take into account DG representation when forming gradient of density.
                VDENDX(IGL)=VDENDX(IGL)-NX(ILOC,GI)*DETWEI(GI)*TDENSUB(J)*N(JLOC,GI)
                VDENDY(IGL)=VDENDY(IGL)-NY(ILOC,GI)*DETWEI(GI)*TDENSUB(J)*N(JLOC,GI)
                VDENDZ(IGL)=VDENDZ(IGL)-NZ(ILOC,GI)*DETWEI(GI)*TDENSUB(J)*N(JLOC,GI)
               ENDIF

                VEDUDX(IGL)=VEDUDX(IGL)+N(ILOC,GI)*DETWEI(GI)*NU(JGL)*NX(JLOC,GI)
                VEDUDY(IGL)=VEDUDY(IGL)+N(ILOC,GI)*DETWEI(GI)*NU(JGL)*NY(JLOC,GI)
                VEDUDZ(IGL)=VEDUDZ(IGL)+N(ILOC,GI)*DETWEI(GI)*NU(JGL)*NZ(JLOC,GI)

                VEDVDX(IGL)=VEDVDX(IGL)+N(ILOC,GI)*DETWEI(GI)*NV(JGL)*NX(JLOC,GI)
                VEDVDY(IGL)=VEDVDY(IGL)+N(ILOC,GI)*DETWEI(GI)*NV(JGL)*NY(JLOC,GI)
                VEDVDZ(IGL)=VEDVDZ(IGL)+N(ILOC,GI)*DETWEI(GI)*NV(JGL)*NZ(JLOC,GI)

                VEDWDX(IGL)=VEDWDX(IGL)+N(ILOC,GI)*DETWEI(GI)*NW(JGL)*NX(JLOC,GI)
                VEDWDY(IGL)=VEDWDY(IGL)+N(ILOC,GI)*DETWEI(GI)*NW(JGL)*NY(JLOC,GI)
                VEDWDZ(IGL)=VEDWDZ(IGL)+N(ILOC,GI)*DETWEI(GI)*NW(JGL)*NZ(JLOC,GI)

              IF(NSUBNVLOC.NE.0) THEN
                VEDUDX(IGL)=VEDUDX(IGL)-NX(ILOC,GI)*DETWEI(GI)*NUSUB(J)*N(JLOC,GI)
                VEDUDY(IGL)=VEDUDY(IGL)-NY(ILOC,GI)*DETWEI(GI)*NUSUB(J)*N(JLOC,GI)
                VEDUDZ(IGL)=VEDUDZ(IGL)-NZ(ILOC,GI)*DETWEI(GI)*NUSUB(J)*N(JLOC,GI)

                VEDVDX(IGL)=VEDVDX(IGL)-NX(ILOC,GI)*DETWEI(GI)*NVSUB(J)*N(JLOC,GI)
                VEDVDY(IGL)=VEDVDY(IGL)-NY(ILOC,GI)*DETWEI(GI)*NVSUB(J)*N(JLOC,GI)
                VEDVDZ(IGL)=VEDVDZ(IGL)-NZ(ILOC,GI)*DETWEI(GI)*NVSUB(J)*N(JLOC,GI)

                VEDWDX(IGL)=VEDWDX(IGL)-NX(ILOC,GI)*DETWEI(GI)*NWSUB(J)*N(JLOC,GI)
                VEDWDY(IGL)=VEDWDY(IGL)-NY(ILOC,GI)*DETWEI(GI)*NWSUB(J)*N(JLOC,GI)
                VEDWDZ(IGL)=VEDWDZ(IGL)-NZ(ILOC,GI)*DETWEI(GI)*NWSUB(J)*N(JLOC,GI)
              ENDIF

                ML(IGL)=ML(IGL)+N(ILOC,GI)*DETWEI(GI)*N(JLOC,GI)
             END DO
          END DO
               maxdz=max(maxdz,MAX_IN_ELE_DZ)
               mindz=min(mindz,MAX_IN_ELE_DZ)
! calc. lumped vertical length scale...           
               DO GI=1,NGI     
                 VERT_DZ(IGL)=VERT_DZ(IGL)+N(ILOC,GI)*DETWEI(GI)*MAX_IN_ELE_DZ 
               END DO
       END DO

    END DO element_loop

    DO IGL=1,NONODS
       VDENDX(IGL)=VDENDX(IGL)/ML(IGL)
       VDENDY(IGL)=VDENDY(IGL)/ML(IGL)
       VDENDZ(IGL)=VDENDZ(IGL)/ML(IGL)

       VEDUDX(IGL)=VEDUDX(IGL)/ML(IGL)
       VEDUDY(IGL)=VEDUDY(IGL)/ML(IGL)
       VEDUDZ(IGL)=VEDUDZ(IGL)/ML(IGL)

       VEDVDX(IGL)=VEDVDX(IGL)/ML(IGL)
       VEDVDY(IGL)=VEDVDY(IGL)/ML(IGL)
       VEDVDZ(IGL)=VEDVDZ(IGL)/ML(IGL)

       VEDWDX(IGL)=VEDWDX(IGL)/ML(IGL)
       VEDWDY(IGL)=VEDWDY(IGL)/ML(IGL)
       VEDWDZ(IGL)=VEDWDZ(IGL)/ML(IGL)  
! NB the vertical element size is the inverse                 
       VERT_DZ(IGL)=ML(IGL)/VERT_DZ(IGL)
    END DO

    ewrite(2,*) 'DENINI,DENGAM,TEMINI:',DENINI,DENGAM,TEMINI
    CALL PMINMX(vDENDX,nonods,'******vDENDX  ')
    CALL PMINMX(vDENDY,nonods,'******vDENDY  ')
    CALL PMINMX(vDENDZ,nonods,'******vDENDZ  ')

    CALL PMINMX(VEDUDZ,nonods,'******VEDUDZ  ')
    CALL PMINMX(VEDVDZ,nonods,'******VEDVDZ  ')
    CALL PMINMX(ML,nonods,'******ML  ')
    
       min_l=1.e+20
       min_q=1.e+20
       min_km=1.e+20
       
       max_l=-1.e+20
       max_q=-1.e+20
       max_km=-1.e+20
       
       min_qqgi=1.e+20
       min_QQLGI=1.e+20
       
       max_qqgi=-1.e+20
       max_QQLGI=-1.e+20

    ! Make sure we are dealing with density and not temp
    DO IGL=1,NONODS
       DENDZ=NORMX(IGL)*VDENDX(IGL)+NORMY(IGL)*VDENDY(IGL)+NORMZ(IGL)*VDENDZ(IGL)
       !  The rotation matrix in 3-D is R=  
       !   T1X   T1Y   T1Z
       !   T2X   T2Y   T2Z
       !   NORMX NORMY NORMZ
       ! USING U'=R U we can differentiate U' assuming R is constant
       DDUDX=VEDUDX(IGL)*T1X(IGL)  + VEDVDX(IGL)*T1Y(IGL) + VEDWDX(IGL)*T1Z(IGL)
       DDVDX=VEDUDX(IGL)*T2X(IGL)  + VEDVDX(IGL)*T2Y(IGL) + VEDWDX(IGL)*T2Z(IGL)

       DDUDY=VEDUDY(IGL)*T1X(IGL)  + VEDVDY(IGL)*T1Y(IGL) + VEDWDY(IGL)*T1Z(IGL)
       DDVDY=VEDUDY(IGL)*T2X(IGL)  + VEDVDY(IGL)*T2Y(IGL) + VEDWDY(IGL)*T2Z(IGL)

       DDUDZ=VEDUDZ(IGL)*T1X(IGL)  + VEDVDZ(IGL)*T1Y(IGL) + VEDWDZ(IGL)*T1Z(IGL)
       DDVDZ=VEDUDZ(IGL)*T2Y(IGL)  + VEDVDZ(IGL)*T2Y(IGL) + VEDWDZ(IGL)*T2Z(IGL)

       ! Now use the eqn like that used to calculate DENDZ
       DUDZ=NORMX(IGL)*DDUDX +NORMY(IGL)*DDUDY +NORMZ(IGL)*DDUDZ
       DVDZ=NORMX(IGL)*DDVDX +NORMY(IGL)*DDVDY +NORMZ(IGL)*DDVDZ

       QQGI=MAX(Q2MIN,TQQ(IGL))
       QQLGI=MAX(Q2LMIN,TQQL(IGL))
       ZPOS=Z(IGL)
       
       min_qqgi=min(min_qqgi,qqgi)
       min_QQLGI=min(min_QQLGI,QQLGI)
       
       max_qqgi=max(max_qqgi,qqgi)
       max_QQLGI=max(max_QQLGI,QQLGI)

       Q=SQRT(MAX(QQGI,0.0))
       L=max(0.0,QQLGI)/MAX(1.E-10,QQGI)
!!! look at gradient of density and comment out next 4 lines        
       IF(HANERT) THEN
          N2=-GRAVTY*DENDZ/DENINI
          LMAX2=(0.28*Q*Q)/MAX(1.E-10,N2)
          L=MIN(L,SQRT(LMAX2))
       ENDIF
       
       CURTQQ(IGL)=MAX(Q2MIN,CURTQQ(IGL))
       CURTQQL(IGL)=MAX(Q2LMIN,CURTQQL(IGL))

       TQQ(IGL)=MAX(Q2MIN,TQQ(IGL))
       TQQL(IGL)=MAX(Q2LMIN,TQQL(IGL))
       PREVTQQ(IGL)=MAX(Q2MIN,PREVTQQ(IGL))
       PREVTQQL(IGL)=MAX(Q2LMIN,PREVTQQL(IGL))

       ! This is the Richardson number...
       IF(HEDONG) THEN
          ! This is Hedongs approach 
          GH=((L*L*GRAVTY)/MAX(1.E-10,QQGI*DENINI))*DENDZ
          ! Limit values to avoid problems (Hedong experience)
          GH=MAX(MIN(GH,0.0233),-0.28)
          SM=(0.4275-3.354*GH)/((1.-34.676*GH)*(1.-6.127*GH))
          IF(HANERT) SM=(0.393-3.085*GH)/(1.0-40.803*GH+212.469*GH**2)
          SH=0.494/(1.-34.676*GH) 
       ELSE
          GH=((L*L*GRAVTY)/MAX(1.E-10,QQGI*DENINI))*DENDZ
          GH=MAX(MIN(GH,0.023),-0.28)
          SH=A2*(1.-6.*A1/B1) / TOLFUN(1.-(3.*A2*B2+18.*A1*A2)*GH)
          SM=(SH*(18.*A1*2.+9.*A1*A2)*GH + A1*(1.-3.*C1-6.*A1/B1)) &
               / TOLFUN(1.-9.*A1*A2*GH)
       ENDIF
       ! These are the vertical diffusivities for temp...
       KM=MAX(Q*L*SM,0.0)
       KH=MAX(Q*L*SH,0.0)
       min_l=min(min_l,l)
       min_q=min(min_q,q)
       min_km=min(min_km,km)
       
       max_l=max(max_l,l)
       max_q=max(max_q,q)
       max_km=max(max_km,km)

       ! KQ=vertical diffusivity of qq and qql
       KQ= 0.2*MAX(Q*L,0.0)
       
       ! Diffusion for density and salinity..
       DIFDEN(IGL)=KH
       ! viscosity...
       DIFVIS(IGL)=KM
       ! diffusivity for the 2 turbulence eqns
       DIFTUB(IGL)=KQ

       ! Source and absorption for  qq:
       TSOURC1(IGL)=  &
            ( 2.*KM*(DUDZ**2+DVDZ**2) + (2.*GRAVTY/DENINI)*KH*DENDZ ) 

       TABSOR1(IGL)=2.*Q/(B1*MAX(1.E-10,L)) 
       ! Source and absorption for  lqq:
       ! K is the von Karman constant.
       ! LWALL is the wall proximity function. 
       K=0.4
       if (gottopbotdis) then
         LWALLM1=1./MAX(1.E-10,TOPDIS(IGL)) &
            +1./MAX(1.E-10,BOTDIS(IGL))
         ! Limit L^{-1} based on reasonable values
         !             LWALLM1=MIN(1./(0.001*(ZMAX-ZMIN)), LWALLM1)
         LWALLM1=MIN(1./(0.001*(BOTDIS(IGL)+TOPDIS(IGL))), LWALLM1)
       else
         LWALLM1=1./MAX(1.E-10,ZMAX-ZPOS) &
                           +1./MAX(1.E-10,ZPOS-ZMIN)
       end if

       IF(HEDONG) THEN
          W=1.0+E2*((L*LWALLM1)/K)**2
          TSOURC2(IGL)= &
               ( E1*L*(KM*(DUDZ**2+DVDZ**2)   &
               +GRAVTY*KH*DENDZ/DENINI ) )
  !!!             -1.0*(Q**3)/B1)
          TABSOR2(IGL)=(L*Q/B1)*E2*((LWALLM1/K)**2)+1.0*Q/(B1*MAX(1.E-10,L))
       ELSE
          ! Not Hedongs version
          W=1.+E2*((L/K)*LWALLM1)
          TSOURC2(IGL)= &
               ( E1*L*(KM*(DUDZ**2+DVDZ**2)   &
               +E3*(GRAVTY/DENINI)*KH*DENDZ ) &
               -((Q**3)/B1)*1.0 )
          TABSOR2(IGL)=(Q/B1)*E2*(1./K)*LWALLM1
       ENDIF

    END DO
    
    max_difvis=-1.e+20
    min_difvis= 1.e+20
      DO IGL=1,NONODS
        max_difvis=MAX(max_difvis,DIFVIS(igl))
        min_difvis=MIN(min_difvis,DIFVIS(igl))
      END DO
      
!      print *,'Q2MIN,Q2lMIN:',Q2MIN,Q2lMIN
!      print *,'min_qqgi,min_QQLGI:',min_qqgi,min_QQLGI
!      print *,'max_qqgi,max_QQLGI:',max_qqgi,max_QQLGI
!      print *,'min_l,min_q,min_km:',min_l,min_q,min_km
!      print *,'max_l,max_q,max_km:',max_l,max_q,max_km
!      PRINT *,'min_difvis,max_difvis:',min_difvis,max_difvis
!       DIFVIS=0.0
   
! use a simple method to calculate vertical diffusivity/viscocity...
    IF(SIMP_DIF) THEN
      DO IGL=1,NONODS
        DENDZ=NORMX(IGL)*VDENDX(IGL)+NORMY(IGL)*VDENDY(IGL)+NORMZ(IGL)*VDENDZ(IGL)
        DIFDEN(IGL)=DT*(VERT_DZ(IGL)**2)*GRAVTY*MAX(0.0,DENDZ)
        DIFVIS(igl)=DIFDEN(IGL)*1.0
        DIFTUB(igl)=DIFDEN(IGL)
      END DO
    ENDIF
    CALL PMINMX(DIFVIS,nonods,'******DIFVIS  ')
    CALL PMINMX(VDENDZ,nonods,'******VDENDZ  ')
    CALL PMINMX(VERT_DZ,nonods,'******VERT_DZ  ')
    CALL PMINMX(TDEN,nonods,'******TDEN  ')
    CALL PMINMX(TEMP,nonods,'******TEMP  ') 
!    print *,'gravty,dt,DENGAM,maxdz,mindz:',gravty,dt,DENGAM,maxdz,mindz
!    stop 3832

    ! Communicate DIFDEN,DIFVIS,DIFTUB...
    IF(PARA.EQ.1) THEN
       CALL HALGET(DIFDEN,NONODS,NONODS,NNODP,halo_tag)
       CALL HALGET(DIFVIS,NONODS,NONODS,NNODP,halo_tag)
       CALL HALGET(DIFTUB,NONODS,NONODS,NNODP,halo_tag)
    ENDIF

    IF(COGRAX) THEN
       DO IGL=1,NONODS
          ! Diffusion for density and salinity..
          TMUZZ(IGL)=DIFDEN(IGL)
          ! viscosity...
          MUPTZZ(IGL)=DIFVIS(IGL)
          ! diffusivity for the 2 turbulence eqns
          TMUZZ1(IGL)=DIFTUB(IGL)
          TMUZZ2(IGL)=DIFTUB(IGL)

          ! TEMP CROSS TERMS...         
          TMUXZ(IGL)=0.0
          TMUYZ(IGL)=0.0
          ! QQ CROSS TERMS...         
          TMUXZ1(IGL)=0.0
          TMUYZ1(IGL)=0.0
          ! QQL CROSS TERMS...         
          TMUXZ2(IGL)=0.0
          TMUYZ2(IGL)=0.0
          ! viscosity cross terms...
          IF(.NOT.LES) THEN
             MUPTXZ(IGL)=SQRT(MUPTXX(IGL)*MUPTZZ(IGL))
             MUPTYZ(IGL)=SQRT(MUPTYY(IGL)*MUPTZZ(IGL))
          ENDIF
       END DO
    ELSE
       ! Else assume a SGS ocean...and zero diffusivity except in vertical...
       ! the diffusion matrix should be R^T D R 
       ! NB  X'=R X
       DO IGL=1,NONODS
          IF(.TRUE.) THEN

             CALL SMLMA3( &
                  TMUXX(IGL),TMUXY(IGL),TMUXZ(IGL),&
                  RT1,       TMUYY(IGL),TMUYZ(IGL),&
                  RT2,       RT3,       TMUZZ(IGL),& 
                  ! =
                  T1X(IGL),T2X(IGL),NORMX(IGL), &
                  T1Y(IGL),T2Y(IGL),NORMY(IGL), &
                  T1Z(IGL),T2Z(IGL),NORMZ(IGL), &            
                  ! *          
                  0.0,   0.0,   0.0,  &
                  0.0,   0.0,   0.0,  &
                  DIFDEN(IGL)*NORMX(IGL),DIFDEN(IGL)*NORMY(IGL),DIFDEN(IGL)*NORMZ(IGL)) 

             ! viscosity in tensor form: TENSOR=R DIAG R^T
             CALL SMLMA3( &
                  MUPTXX(IGL),MUPTXY(IGL),MUPTXZ(IGL),&
                  RT1,        MUPTYY(IGL),MUPTYZ(IGL),&
                  RT2,        RT3,        MUPTZZ(IGL),& 
                  ! =
                  T1X(IGL),T2X(IGL),NORMX(IGL), &
                  T1Y(IGL),T2Y(IGL),NORMY(IGL), &
                  T1Z(IGL),T2Z(IGL),NORMZ(IGL), &                
                  ! *          
                  0.0,   0.0,   0.0,  &
                  0.0,   0.0,   0.0,  &
                  DIFVIS(IGL)*NORMX(IGL),DIFVIS(IGL)*NORMY(IGL),DIFVIS(IGL)*NORMZ(IGL)) 

             ! TURBULENCE QQ: TENSOR=R DIAG R^T
             CALL SMLMA3( &
                  TMUXX1(IGL),TMUXY1(IGL),TMUXZ1(IGL),&
                  RT1,        TMUYY1(IGL),TMUYZ1(IGL),&
                  RT2,        RT3,        TMUZZ1(IGL),& 
                  ! =
                  T1X(IGL),T2X(IGL),NORMX(IGL), &
                  T1Y(IGL),T2Y(IGL),NORMY(IGL), &
                  T1Z(IGL),T2Z(IGL),NORMZ(IGL), &                
                  ! *          
                  0.0,   0.0,   0.0,  &
                  0.0,   0.0,   0.0,  &
                  DIFTUB(IGL)*NORMX(IGL),DIFTUB(IGL)*NORMY(IGL),DIFTUB(IGL)*NORMZ(IGL)) 
          ENDIF

       END DO
    ENDIF
    ! QQ and QQL have the same diffusivity
       TMUXX2(1:NONODS) = TMUXX1(1:NONODS)
       TMUXY2(1:NONODS) = TMUXY1(1:NONODS)
       TMUXZ2(1:NONODS) = TMUXZ1(1:NONODS)
       TMUYY2(1:NONODS) = TMUYY1(1:NONODS)
       TMUYZ2(1:NONODS) = TMUYZ1(1:NONODS)
       TMUZZ2(1:NONODS) = TMUZZ1(1:NONODS)

    IF(GOTSALT) THEN
       ! Salt diffusivity
       SMUXX(1:NONODS) = TMUXX(1:NONODS)
       SMUXY(1:NONODS) = TMUXY(1:NONODS)
       SMUXZ(1:NONODS) = TMUXZ(1:NONODS)
       SMUYY(1:NONODS) = TMUYY(1:NONODS)
       SMUYZ(1:NONODS) = TMUYZ(1:NONODS)
       SMUZZ(1:NONODS) = TMUZZ(1:NONODS)
    ENDIF
    
    ! add in the b.c's into the eqns 
    ! NB we are solving for rates of change and the 
    ! b.c's reflect this. Also turbulent variables 
    ! are from the previous time step...
    DO II=1,NOBCTQQ
       NOD=BCT2QQ(II)
       FUTVAL=(B1**(2./3.))*(NU(NOD)**2+NV(NOD)**2+NW(NOD)**2)
 !!      BCT1QQ(II)=(FUTVAL-TQQ(NOD))/DT
       BCT1QQ(II)=(0.0-PREVTQQ(NOD))/DT
    END DO
    DO II=1,NOBCTLQQ
       NOD=BCT2LQQ(II)
       BCT1LQQ(II)=(0.0-PREVTQQL(NOD))/DT
    END DO

    ewrite(2,*) 'gotvisfield=', gotvisfield
    IF(GOTVISFIELD) THEN
       ! Store the vertical viscosity and diffusivity for visualising...
       TVISVIS(1:NONODS) = DIFVIS(1:NONODS)
       TVISDIF(1:NONODS) = DIFDEN(1:NONODS)
    ENDIF
    ewrite(2,*) 'DENINI,DENGAM,TEMINI:',DENINI,DENGAM,TEMINI
    CALL PMINMX(TEMP,nonods,'******TEMP  ')
    CALL PMINMX(TSOURC1,nonods,'******TSOURC1  ')
    CALL PMINMX(TSOURC2,nonods,'******TSOURC2  ')
    CALL PMINMX(TABSOR1,nonods,'******TABSOR1  ')
    CALL PMINMX(MUPTZZ,nonods,'******MUPTZZ  ')
    CALL PMINMX(TMUZZ,nonods,'******TMUZZ  ')
    CALL PMINMX(difden,nonods,'******diSfden  ')
    CALL PMINMX(TMUZZ1,nonods,'******TMUZZ1  ')
    CALL PMINMX(TMUZZ2,nonods,'******TMUZZ2  ')
    CALL PMINMX(TQQ,nonods,'******TQQ  ')
    CALL PMINMX(TQQL,nonods,'******TQQL  ')
    CALL PMINMX(NU,nonods,'******NU  ')
    CALL PMINMX(NV,nonods,'******NV  ')
    CALL PMINMX(NW,nonods,'******NW  ')
    CALL PMINMX(TVISVIS,nonods,'******TVISVIS  ')
    CALL PMINMX(TVISDIF,nonods,'******TVISDIF  ')

! DEALLOCATE.......****************************  
         DEALLOCATE(ML)
         DEALLOCATE(VDENDX)
         DEALLOCATE(VDENDY)
         DEALLOCATE(VDENDZ)
         DEALLOCATE(VEDUDX)
         DEALLOCATE(VEDUDY)
         DEALLOCATE(VEDUDZ)
         DEALLOCATE(VEDVDX)
         DEALLOCATE(VEDVDY)
         DEALLOCATE(VEDVDZ)
         DEALLOCATE(VEDWDX)
         DEALLOCATE(VEDWDY)
         DEALLOCATE(VEDWDZ)
         DEALLOCATE(DIFDEN)
         DEALLOCATE(DIFVIS)
         DEALLOCATE(DIFTUB)
         DEALLOCATE(WORK)
         DEALLOCATE(TDEN)
         DEALLOCATE(T1X)
         DEALLOCATE(T1Y)
         DEALLOCATE(T1Z)
         DEALLOCATE(T2X)
         DEALLOCATE(T2Y)
         DEALLOCATE(T2Z) 
         DEALLOCATE(NORMX)
         DEALLOCATE(NORMY)
         DEALLOCATE(NORMZ)
! *********END DEALOCAT MEMORY**************

 END SUBROUTINE MELLOR_YAMADA_TURBULENCE

 subroutine mellor_yamada_state(state,&
       DENGAM,DENINI,TEMINI,GRAVTY,&
       COGRAX,BSOUX,BSOUY,BSOUZ,&
       LES,&
       ! Salt...
       salinity_option,DENGAM_sal,S0, gottopbotdis, &
       ! BCS...
       NOBCTQQ,BCT1QQ,BCT2QQ,&
       NOBCTLQQ,BCT1LQQ,BCT2LQQ,DT,&
       ! parallel stuff...
       nnodp,para,halo_tag)
   !!< Calculate mellor yamada terms based on state.
   !!< kinetic energy field.
   type(state_type), intent(inout) :: state
   ! The following arguments are legacy code which should one day be
   ! replaced but as they do not explicitly depend on rmem they are being
   ! left for the time being.
   real, intent(in) :: dengam, denini, temini, gravty
   logical, intent(in) :: cograx, les
   real, intent(in) :: BSOUX,BSOUY,BSOUZ
   integer, intent(in) :: salinity_option
   real, intent(in) :: DENGAM_sal,S0
   logical, intent(in) :: gottopbotdis
   integer, intent(in) :: NOBCTQQ, NOBCTLQQ
   real, dimension(:), intent(in) :: BCT1QQ, BCT1LQQ
   integer, dimension(:), intent(in) :: BCT2QQ, BCT2LQQ
   real, intent(in) :: dt
   integer, intent(in) :: nnodp, para, halo_tag

   type(scalar_field) :: oldqq, iteratedqq, qq, qqsource, qqabs
   type(scalar_field) :: oldqql, iteratedqql, qql, qqlsource, qqlabs
   
   type(vector_field) :: X, NU, gravitydirection, NU_subgrid

   type(scalar_field) :: temperature, salinity
   type(tensor_field) :: this_tensor
   logical :: salinity_present, d3

   type(scalar_field) :: temperature_subgrid, salinity_subgrid
   integer :: tsub_present, ssub_present, nusub_present
   real, dimension(:), pointer :: tsub, ssub, nusub, nvsub, nwsub
   real, dimension(:), pointer :: vgravx, vgravy, vgravz

   real, dimension(0), target :: zero

   ! Old style element information
   real, dimension(:, :), allocatable :: n, nlx, nly, nlz
   
   ! Old style storage for diffusivity and viscosity
   real, dimension(:), allocatable ::  &
        TMUXX,TMUXY,TMUXZ, &
        TMUYY,TMUYZ,TMUZZ, &
        SMUXX,SMUXY,SMUXZ, &
        SMUYY,SMUYZ,SMUZZ, &
        TMUXX1,TMUXY1,TMUXZ1, &
        TMUYY1,TMUYZ1,TMUZZ1, &
        TMUXX2,TMUXY2,TMUXZ2, &
        TMUYY2,TMUYZ2,TMUZZ2, &
        MUPTXX,MUPTXY,MUPTXZ,&
        MUPTYY,MUPTYZ,MUPTZZ


   real, dimension(:), allocatable :: TVISVIS,TVISDIF

   integer :: stat, nonods, totele, nloc, ngi

   iteratedqq=extract_scalar_field(state, "IteratedKineticEnergy")
   qq=extract_scalar_field(state, "KineticEnergy")
   oldqq=extract_scalar_field(state, "OldKineticEnergy")
   qqsource=extract_scalar_field(state, "KineticEnergySource")
   qqabs=extract_scalar_field(state, "KineticEnergyAbsorption")

   iteratedqql=extract_scalar_field(state, "IteratedTurbulentLengthScalexKineticEnergy")
   qql=extract_scalar_field(state, "TurbulentLengthScalexKineticEnergy")
   oldqql=extract_scalar_field(state, "OldTurbulentLengthScalexKineticEnergy")
   qqlsource=extract_scalar_field(state, "TurbulentLengthScalexKineticEnergySource")
   qqlabs=extract_scalar_field(state, "TurbulentLengthScalexKineticEnergyAbsorption")

   X=extract_vector_field(state, "Coordinate")
   NU=extract_vector_field(state, "NonlinearVelocity")

   d3=(X%dim==3)

   temperature=extract_scalar_field(state, "Temperature")

   salinity=extract_scalar_field(state, "Salinity", stat)   
   if (stat==0) then
      ! Extra reference to make the deallocate below safe.
      call incref(salinity)
      salinity_present=.true.
   else
      call allocate(salinity, temperature%mesh, "Salinity")
      salinity_present=.false.
   end if

   temperature_subgrid=extract_scalar_field(state,&
        & "TemperatureInnerElement", stat)
   if (stat==0) then
      tsub=>temperature_subgrid%val
      tsub_present=1
   else
      tsub=>zero
      tsub_present=0
   end if

   salinity_subgrid=extract_scalar_field(state,&
        & "SalinityInnerElement", stat)
   if (stat==0) then
      ssub=>salinity_subgrid%val
      ssub_present=1
   else
      ssub=>zero
      ssub_present=0
   end if

   NU_subgrid=extract_vector_field(state, "NonlinearVelocityInnerElement",&
        & stat)
   if (stat==0) then
      nusub=>NU_subgrid%val(1)%ptr
      nvsub=>NU_subgrid%val(2)%ptr
      nwsub=>NU_subgrid%val(3)%ptr
      nusub_present=1
   else
      nusub=>zero
      nvsub=>zero
      nwsub=>zero
      nusub_present=0
   end if

   gravitydirection=extract_vector_field(state,"GravityDirection")
   vgravx=>gravitydirection%val(1)%ptr
   vgravy=>gravitydirection%val(2)%ptr
   vgravz=>gravitydirection%val(3)%ptr

   nonods=node_count(NU)
   totele=element_count(nu)
   ngi=ele_ngi(nu,1)
   nloc=ele_loc(nu,1)
   
   ! These aren't currently supported but need to be here to keep the
   ! mellor_yamada call legal.
   allocate(TVISVIS(nonods),TVISDIF(nonods))

   ! Allocate space for the viscosities and diffusivities and copy them
   ! back afterwards. This is required because the memory layout of
   ! tensor_field%val is not the same as the old code.
   allocate(TMUXX(NONODS),TMUXY(NONODS),TMUXZ(NONODS),&
        TMUYY(NONODS),TMUYZ(NONODS),TMUZZ(NONODS),&
        SMUXX(NONODS),SMUXY(NONODS),SMUXZ(NONODS),&
        SMUYY(NONODS),SMUYZ(NONODS),SMUZZ(NONODS),&
        TMUXX1(NONODS),TMUXY1(NONODS),TMUXZ1(NONODS),&
        TMUYY1(NONODS),TMUYZ1(NONODS),TMUZZ1(NONODS),&
        TMUXX2(NONODS),TMUXY2(NONODS),TMUXZ2(NONODS),&
        TMUYY2(NONODS),TMUYZ2(NONODS),TMUZZ2(NONODS),&
        MUPTXX(NONODs),MUPTXY(NONODs),MUPTXZ(NONODs),&
        MUPTYY(NONODs),MUPTYZ(NONODs),MUPTZZ(NONODs))


   ! Copy the tensor quantities from the tensor fields.
   this_tensor=extract_tensor_field(state, "TemperatureDiffusivity", stat)
   if (stat==0) then
      call copy_tensor_to_legacy(this_tensor, TMUXX,TMUXY,TMUXZ, &
           TMUYY,TMUYZ,TMUZZ)
   end if
      
   this_tensor=extract_tensor_field(state, "SalinityDiffusivity", stat)
   if (stat==0) then
      call copy_tensor_to_legacy(this_tensor, SMUXX,SMUXY,SMUXZ, &
           SMUYY,SMUYZ,SMUZZ)
   end if

   this_tensor=extract_tensor_field(state, "KineticEnergyDiffusivity", stat)
   if (stat==0) then
      call copy_tensor_to_legacy(this_tensor, TMUXX1,TMUXY1,TMUXZ1, &
           TMUYY1,TMUYZ1,TMUZZ1)
   end if
      
   this_tensor=extract_tensor_field(state, &
        "TurbulentLengthScalexKineticEnergyDiffusivity", stat)
   if (stat==0) then
      call copy_tensor_to_legacy(this_tensor, TMUXX2,TMUXY2,TMUXZ2, &
           TMUYY2,TMUYZ2,TMUZZ2)
   end if
      
   this_tensor=extract_tensor_field(state, "Viscosity", stat)
   if (stat==0) then
      call copy_tensor_to_legacy(this_tensor, MUPTXX,MUPTXY,MUPTXZ, &
           MUPTYY,MUPTYZ,MUPTZZ)
   end if

   
   allocate(n(ele_loc(NU%mesh, 1), ele_ngi(NU%mesh, 1)))
   allocate(nlx(ele_loc(NU%mesh, 1), ele_ngi(NU%mesh, 1)))
   allocate(nly(ele_loc(NU%mesh, 1), ele_ngi(NU%mesh, 1)))
   allocate(nlz(ele_loc(NU%mesh, 1), ele_ngi(NU%mesh, 1)))
   call extract_old_element(ele_shape(NU%mesh, 1), n, nlx, nly, nlz)

   call mellor_yamada_turbulence(state, node_count(qq),node_count(qqsource),&
       iteratedqq%val,iteratedqql%val,&
       qq%val,qql%val,oldqq%val,oldqql%val,&
       NU%val(1)%ptr,NU%val(2)%ptr,NU%val(3)%ptr,&
       X%val(1)%ptr,X%val(2)%ptr,X%val(3)%ptr,D3,.false.,&
       temperature%val,&
       tsub_present,tsub,ssub_present,ssub, &
       nusub_present,NUSUB,NVSUB,NWSUB, &
       DENGAM,DENINI,TEMINI,GRAVTY,&
       COGRAX,BSOUX,BSOUY,BSOUZ,&
       VGRAVX,VGRAVY,VGRAVZ,&
       LES,&
       ! Salt...
       salinity_present,salinity_option,DENGAM_sal,S0,salinity%val, gottopbotdis, &
       ! This is the output...
       TMUXX,TMUXY,TMUXZ, &
       TMUYY,TMUYZ,TMUZZ, &
       SMUXX,SMUXY,SMUXZ, &
       SMUYY,SMUYZ,SMUZZ, &
       TMUXX1,TMUXY1,TMUXZ1, &
       TMUYY1,TMUYZ1,TMUZZ1, &
       TMUXX2,TMUXY2,TMUXZ2, &
       TMUYY2,TMUYZ2,TMUZZ2, &
       MUPTXX,MUPTXY,MUPTXZ,&
       MUPTYY,MUPTYZ,MUPTZZ,&
       qqsource%val,qqabs%val,&
       qqlsource%val,qqlabs%val,&
       ! for visualising viscosity and diffusivity...
       .true.,TVISVIS,TVISDIF, &
       ! Integer arrays and basis functions
       TOTELE,NLOC,NGI,NU%mesh%ndglno,X%mesh%ndglno,&
       N,NLX,NLY,NLZ, nu%mesh%shape%quadrature%WEIGHT,&
       ! BCS...
       NOBCTQQ,BCT1QQ,BCT2QQ,&
       NOBCTLQQ,BCT1LQQ,BCT2LQQ,DT,&
       ! parallel stuff...
       nnodp,para,halo_tag)
   
   ! Copy the tensor quantities back to the tensor fields.
   this_tensor=extract_tensor_field(state, "TemperatureDiffusivity", stat)
   if (stat==0) then
      call copy_tensor_from_legacy(this_tensor, TMUXX,TMUXY,TMUXZ, &
           TMUYY,TMUYZ,TMUZZ)
   end if
      
   this_tensor=extract_tensor_field(state, "SalinityDiffusivity", stat)
   if (stat==0) then
      call copy_tensor_from_legacy(this_tensor, SMUXX,SMUXY,SMUXZ, &
           SMUYY,SMUYZ,SMUZZ)
   end if

   this_tensor=extract_tensor_field(state, "KineticEnergyDiffusivity", stat)
   if (stat==0) then
      call copy_tensor_from_legacy(this_tensor, TMUXX1,TMUXY1,TMUXZ1, &
           TMUYY1,TMUYZ1,TMUZZ1)
   end if
      
   this_tensor=extract_tensor_field(state, &
        "TurbulentLengthScalexKineticEnergyDiffusivity", stat)
   if (stat==0) then
      call copy_tensor_from_legacy(this_tensor, TMUXX2,TMUXY2,TMUXZ2, &
           TMUYY2,TMUYZ2,TMUZZ2)
   end if
      
   this_tensor=extract_tensor_field(state, "Viscosity", stat)
   if (stat==0) then
      call copy_tensor_from_legacy(this_tensor, MUPTXX,MUPTXY,MUPTXZ, &
           MUPTYY,MUPTYZ,MUPTZZ)
   end if
   

   ! Drop extra field references.
   call deallocate(salinity)

 end subroutine mellor_yamada_state

end module mellor_yamada
