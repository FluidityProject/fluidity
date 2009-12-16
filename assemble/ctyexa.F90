!    Copyright (C) 2006 Imperial College London and others.
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

module ctyexa_module

  use FLDebug
  use solvers
  use geostrophic
  use sparse_tools
  use spud
  use global_parameters, only: OPTION_PATH_LEN, phase2state_index
  use futils
  use hart3d_allsorts
  use ctmult_module
  use allsorts
  use flcomms_module
  use mulmat_module
  use mulinv_module
  use parallel_tools
  use fields
  use state_module
  use free_surface_module
  use halos

  implicit none

  private
 
  public :: ctyexa, funadp

contains

  SUBROUTINE CTYEXA(U,V,W,P,RHS,DP,&
       &     NDWISP,KCMC,&
       &     CMC,FINCMC,COLCMC,NCMC,&
       &     NCMCB,NPRESS,NPROPT,ALFST2,SPRESS, &
       &     C1T,C2T,C3T,FINDCT,COLCT,NCT, &
       &     ML,&
       &     NBUOY,BOY_ML, &
       &     DM1,DM2,DM3,DT,&
       &     VECL,VECSEC,&
       &     FREDOP,NONODS,ILENG,&
       !     The only time ILENG can not be zero is if ROTAT=true.
       &     CONMAS,D3,&
       &     NDPSET,PARA,halo_tag,&
       &     halo_tag_p,NNODP,NNODPP,&
       !     Th following is for rotations only ROTAT=.TRUE.
       &     ROTAT,NNODRO,NODROT,&
       &     NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
       !     The multiphase stuff...
       &     NPHASE,IPHASE,&
       &     C1TP,C2TP,C3TP,&
       ! FOr traffic 
       &      VOLTRAF,UTRAF,VTRAF,WTRAF,GOTTRAF,&
       !     --------------------------------------START OF ADD BY CRGW 20/03/06
       !           SINGLE PHASE SOLID COMPRESSIBILITY:
       &      MKCOMP, DIVQS, &
       ! Really this is all you need:
       &      state)
    !     
    !     This sub corrects the vels to satisfy discretised cty exactly
    !     using a prejection method. 
    !     DIVQ is the discetised cty vector. 
    INTEGER NCMC
    INTEGER MATSTR
    PARAMETER(MATSTR=0)
    !     ALFST2= max volume fraction of solids. 
    !     SPRESS=2nd pressure scaling parameter.
    !     
    !     NDPSET contains the pressure node to set to zero.
    INTEGER NDPSET

    INTEGER J,COL,NPHASE
    INTEGER CONMAS,COUPLE,D3,D2HALF,UZAWA
    INTEGER FREDOP,NONODS,NCT,ILENG,NCMCB,NPRESS
    INTEGER NPROPT
    REAL ALFST2,SPRESS
    REAL RHS(FREDOP*NPRESS),DP(FREDOP*NPRESS),P(FREDOP*NPRESS)
    REAL U(NONODS*NPHASE),V(NONODS*NPHASE),W(NONODS*NPHASE)
    REAL DM1(NONODS),DM2(NONODS),DM3(NONODS),DT
    logical NDWISP
    REAL KCMC(NCMC)
    REAL CMC(NCMCB)
    REAL ML(NONODS)
    REAL C1T(NCT*NPRESS),C2T(NCT*NPRESS),C3T(NCT*NPRESS)
    REAL C1TP(NCT*NPHASE),C2TP(NCT*NPHASE),C3TP(NCT*NPHASE)
    !     NB GLM already has boundary conditions in. 
    INTEGER FINCMC(FREDOP+1),COLCMC(NCMC)
    INTEGER FINDCT(FREDOP+1),COLCT(NCT)
    REAL VECL(ILENG),VECSEC(ILENG)
    INTEGER NBUOY
    REAL BOY_ML(NBUOY*9*NONODS)
    INTEGER I
    !     For parallel processing ...NB then end of R contains memory for parallel solution
    INTEGER NNODP,NNODPP
    !     NNODP   no of vel nodes in subdomain excluding halo nodes.
    !     NNODPP no of pressure nodes in subdomain excluding halo nodes.
    INTEGER  PARA,halo_tag
    INTEGER halo_tag_p
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
    LOGICAL ROTAT
    !     For multiphase stuff ONLY...
    !     DM1=DM1PHA(I+NONODS*(IPHASE-1)) contains b.c's of IPHASE.
    !     DM2=DM2PHA(I+NONODS*(IPHASE-1)) contains b.c's of IPHASE.
    !     DM3=DM3PHA(I+NONODS*(IPHASE-1)) contains b.c's of IPHASE.
    !     ML=MLPHA(I+NONODS*(IPHASE-1)) contains lumped mass*density of 
    !     current phase IPHASE.
    INTEGER IPHASE
    REAL SUMDIV
    !     INTEGER INCTY2
    !     SAVE INCTY2
    !     DATA INCTY2 /0/
    !     MPDIAG contains the diagonals of the block matrix used instead 
    !     of the mass matrix for the non-sym multi-phase projection method.
    !     INVPHA is the inverseof MPDIAG with b.c's added. 
    LOGICAL FINISH
    INTEGER ITS, NPLOO, IPLOO
    INTEGER IDIM
    !     RTDR contains RT D R   ^T where D contains DMI's and R is the rotation matrix. 
    REAL, ALLOCATABLE, DIMENSION(:)::BOY_INV

    REAL VOLTRAF(nonods),UTRAF(nonods), VTRAF(nonods),WTRAF(nonods)
    INTEGER nod, GOTTRAF

    type(csr_matrix) :: cmc_matrix, prolongator
    type(csr_sparsity) :: cmc_sparsity
    type(halo_type) :: halo
    integer :: stat
    character(len=OPTION_PATH_LEN) pressure_option_path
    integer :: istate
    real, dimension(:), allocatable:: work
    real, dimension(:), allocatable:: worku, workv, workw
    !     --------------------------------------START OF ADD BY CRGW 20/03/06
    !           SINGLE PHASE SOLID COMPRESSIBILITY:
    INTEGER MKCOMP
    REAL DIVQS(FREDOP)
    !     --------------------------------------END OF ADD BY CRGW 20/03/06
    type(state_type), dimension(:), intent(inout) :: state
      
    type(scalar_field), pointer :: pressure

    ewrite(3,*) 'NCMC,NCMCB,NPRESS,NPROPT,ALFST2,SPRESS:',&
         &     NCMC,NCMCB,NPRESS,NPROPT,ALFST2,SPRESS

    istate = phase2state_index(IPHASE)
    ewrite(2,*) 'IPHASE, istate = ', IPHASE, istate

    pressure_option_path='/material_phase['//int2str(istate-1)//&
         & ']/scalar_field::Pressure'

    do ITS=1,1
       ewrite(3,*) 'ITS=',ITS
       !     
       FINISH=.FALSE.
       NPLOO =1
       do IPLOO=1,NPLOO
          if (IsParallel()) then
            CALL HALGET(U,NONODS,NONODS,NNODP,halo_tag)
            CALL HALGET(V,NONODS,NONODS,NNODP,halo_tag)
            IF(D3==1) then
               CALL HALGET(W,NONODS,NONODS,NNODP,halo_tag)
            ENDIF
          end if
          !     
          !     
          !     Form right hand side of pressure correction equation.
          SUMDIV=0.
          RHS=0.
          !     
          !           ----------------------------- CRGW ADDED 20/03/06
          !           IF COMPRESSIBLE MODIFY RHS:
          IF(MKCOMP.GE.1) THEN
             SUMDIV=0.
             EWRITE(2,*) 'ADDING DIVQS TO RHS IN CTYEXA:'
             EWRITE_MINMAX(DIVQS)
             do I=1,FREDOP
                SUMDIV=SUMDIV+DIVQS(I)
                RHS(I)=DIVQS(I)/DT
             END DO
          ENDIF
          !           ----------------------------- CRGW END OF ADD 20/03/06
          IF(NDWISP) THEN
             allocate(work(1:FREDOP))
             !     ADD R=KCMC*P, RHS=RHS-R
             CALL MULMAT(work,P,KCMC,FINCMC,COLCMC,&
                  &              NCMC,NCMC,MATSTR,FREDOP,FREDOP,NNODPP,.false.)
             do I=1,FREDOP
                RHS(I)=RHS(I)-work(i)
             END DO
             deallocate(work)
          ENDIF
          !     
          ewrite(3,*) 'dt:',dt 
          ewrite(3,*)'INSIDE CTYEXA HERE1'
          !     
          IF (GOTTRAF.EQ.1)THEN
            !     RHS=RHS-(C1T (1-voltraf)*U+C2T(1-voltraf)V+C3T(1-voltraf)W)/DT
            !     ewrite(3,*) 'here 121 rhs=',rhs  c comment this out temporarily 12may05

            allocate(worku(1:nonods), workv(1:nonods), workw(1:nonods))
            
            worku=(1.0-voltraf)*u
            workv=(1.0-voltraf)*v
            workw=(1.0-voltraf)*w

            CALL PMINMX(VOLTRAF,NONODS,'******VOLTRAF  ')
            CALL PMINMX(UTRAF,NONODS,'******UTRAF  ')          
            CALL PMINMX(VTRAF,NONODS,'******VTRAF  ')          
            CALL PMINMX(WTRAF,NONODS,'******WTRAF  ')          

            ewrite(3,*)  'CALLING CTMULT 1ST TIME'

            CALL CTMULT(RHS,worku,workv,workw, DT, &
                  &                 .FALSE.,NDPSET, FREDOP,NONODS,D3,  &
                  &              C1T,C2T,C3T,FINDCT,COLCT,NCT, &
                  !     Parallel ...
                  &                 PARA,halo_tag,&
                  &                 halo_tag_p,NNODP,NNODPP )
            

            ewrite(3,*)&
                  & 'In CTYEXA.F: AFTER CALL TO CTMULT and prior ro 2nd Phase R'
            ewrite(3,*)   'GOTTRAF', GOTTRAF,'Time step=', DT


            CALL PMINMX(VOLTRAF,NONODS,'******VOLTRAF  ')
            CALL PMINMX(UTRAF,NONODS,'******UTRAF  ')          
            CALL PMINMX(VTRAF,NONODS,'******VTRAF  ')          
            CALL PMINMX(WTRAF,NONODS,'******WTRAF  ')          
            !     CALL PMINMX(r,NONODS,'******R  ')       

            ewrite(3,*)  'Calculating R for 2nd phase'  

            
            worku=voltraf*u
            workv=voltraf*v
            workw=voltraf*w

            CALL PMINMX(VOLTRAF,NONODS,'******VOLTRAF  ')
            CALL PMINMX(UTRAF,NONODS,'******UTRAF  ')          
            CALL PMINMX(VTRAF,NONODS,'******VTRAF  ')          
            CALL PMINMX(WTRAF,NONODS,'******WTRAF  ')          
            !     CALL PMINMX(r,NONODS,'******R  ')       

            ewrite(3,*) 'In CTYEXA.F - CALL CTMULT for 2nd Phase'

            CALL CTMULT(RHS, worku, workv, workw, DT, &
                  &                 .FALSE.,NDPSET, FREDOP,NONODS,D3,  &
                  &                 C1T,C2T,C3T,FINDCT,COLCT,NCT, &
                  !     Parallel ...
                  &                 PARA,halo_tag,&
                  &              halo_tag_p,NNODP,NNODPP )
            ewrite(3,*)  'IN CTYEXA.F: AFTER 2nd CALL to CTMULT'

            CALL PMINMX(VOLTRAF,NONODS,'******VOLTRAF  ')
            CALL PMINMX(UTRAF,NONODS,'******UTRAF  ')          
            CALL PMINMX(VTRAF,NONODS,'******VTRAF  ')          
            CALL PMINMX(WTRAF,NONODS,'******WTRAF  ')
            
            deallocate(worku,workv,workw)

          ELSE

            !     RHS=RHS-(C1TU+C2TV+C3TW)/DT
            !     ewrite(3,*) 'here 121 rhs=',rhs
            EWRITE(2,*) 'GOING INTO CTMULT:'
            IF(MKCOMP.EQ.3) THEN
                CALL CTMULT(RHS,U,V,W, DT, &
                    &                   .FALSE.,NDPSET, FREDOP,NONODS,D3,  &
                    &                   C1TP,C2TP,C3TP,FINDCT,COLCT,NCT, &
                    !     Parallel ...
                    &                   PARA,halo_tag,&
                    &                   halo_tag_p,NNODP,NNODPP )
                !                 
            ELSE  
                CALL CTMULT(RHS,U,V,W, DT, &
                    &                    .FALSE.,NDPSET, FREDOP,NONODS,D3,  &
                    &                    C1T,C2T,C3T,FINDCT,COLCT,NCT, &
                    !     Parallel ...
                    &                    PARA,halo_tag,&
                    &                    halo_tag_p,NNODP,NNODPP )
            ENDIF
          ENDIF
          !     
          !     set pressure correction DP to zero.
          DP(1:FREDOP*NPRESS) = 0.0
          !     Solve pressure correction equation.
          !     ILENG=0
          UZAWA=0
          COUPLE=0
          D2HALF=0
          
#ifdef HAVE_PETSC
          call import_halo(halo_tag_p, halo, stat)
          
          if (stat==0) then
             cmc_sparsity=wrap(fincmc,&
                  &                     colm=colcmc, row_halo=halo,&
                  &                     column_halo=halo,&
                  &                     name='CMCSparsity')
             cmc_matrix=wrap(cmc_sparsity, val=cmc, name='CMCMatrix')
             call deallocate(halo)
          else
             cmc_sparsity=wrap(fincmc, colm=colcmc, name='CMCSparsity')
             cmc_matrix=wrap(cmc_sparsity, val=cmc, name='CMCMatrix')
          end if
          call deallocate(cmc_sparsity)
          if (have_option(trim(pressure_option_path)// &
             "/prognostic/solver/preconditioner::mg/vertical_lumping")) then
               
             if (.not. has_csr_matrix(state(istate), "VerticalProlongationOperator")) then
                pressure => extract_scalar_field(state(istate), "Pressure")
                prolongator=vertical_prolongator_from_free_surface(state(istate), pressure%mesh)
                call insert(state(istate), prolongator, name="VerticalProlongationOperator")
                call deallocate(prolongator)
             end if
             
             prolongator = extract_csr_matrix(state(istate), &
               "VerticalProlongationOperator")
             call petsc_solve(dp, cmc_matrix, rhs, &
               option_path=pressure_option_path, &
               prolongator=prolongator)
          else
             call petsc_solve(dp, cmc_matrix, rhs,&
                  &                 option_path=pressure_option_path)
          end if
          call deallocate(cmc_matrix)
          ewrite(2,*) "PETSc has solved pressure."
#else
          FLAbort("PETSc needed here.")
#endif

          ewrite_minmax(DP)

          do I=1,FREDOP
            P(I)=P(I)+DP(I)
          END DO
          ewrite_minmax(P)

          ewrite(1,*)'IN CTYEXA HAVE JUST finished pressure solve'

          IF(PARA.EQ.1) CALL HALGET(P,FREDOP,FREDOP*NPRESS,NNODPP,halo_tag_p)
          !     
          IF(ROTAT.AND.(CONMAS.EQ.1)) THEN
             !     
             ewrite(3,*) 'here d1'
             do I=1,ILENG
                VECSEC(I)=0.
             END DO
             IF(D3.EQ.0) THEN
                do J=1,FREDOP
                   !     DO 307 J=1,NNODPP
                   do I=FINDCT(J),FINDCT(J+1)-1 
                      COL=COLCT(I)
                      VECSEC(COL)       =VECSEC(COL)        + C1T(I)*DP(J)
                      VECSEC(COL+NONODS)=VECSEC(COL+NONODS) + C2T(I)*DP(J)
                   END DO
                END DO
             ELSE
                do J=1,FREDOP
                   !     DO 327 J=1,NNODPP
                   do I=FINDCT(J),FINDCT(J+1)-1 
                      COL=COLCT(I)
                      VECSEC(COL)         =VECSEC(COL)          + C1T(I)*DP(J)
                      VECSEC(COL+NONODS)  =VECSEC(COL+NONODS)   + C2T(I)*DP(J)
                      VECSEC(COL+2*NONODS)=VECSEC(COL+2*NONODS) + C3T(I)*DP(J)
                   END DO
                END DO
             ENDIF
                  
             !     Form vecl=(R M_L R^T+D)^{-1} vecsec
             !     NB R R^T=I so this sub is not needed. 
             IDIM=2
             IF(ILENG.EQ.3*NONODS) IDIM=3

             CALL MULINV(VECL,VECL(NONODS+1),&
                     &                 VECL((IDIM-1)*NONODS+1),ML,ML,ML,&
                     &                 VECSEC,VECSEC(NONODS+1),VECSEC((IDIM-1)*NONODS+1),&
                     &                 NONODS,D3,&
                     !     Th following is for rotations only ROTAT=.TRUE.
                     !     &         ROTAT,NNODRO,  NODROT,
                     &                 ROTAT, 0,  NODROT,&
                     &                 NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
                     &                 DM1,DM2,DM3 )

             ewrite(3,*) 'here d4'
             !     
             IF(D3.EQ.0) THEN
                do I=1,NONODS
                   U(I)=U(I) + DT*VECL(I)
                   V(I)=V(I) + DT*VECL(I+NONODS)
                END DO
             ELSE
                do I=1,NONODS
                   U(I)=U(I) + DT*VECL(I)
                   V(I)=V(I) + DT*VECL(I+NONODS)
                   W(I)=W(I) + DT*VECL(I+2*NONODS)
                END DO
             ENDIF
             !     
             !     ENDOF(IF (ROTAT)
          ENDIF
          !     
          !     
          !     IF(.NOT.ROTAT.AND.(CONMAS.EQ.0)) THEN
          IF(CONMAS.EQ.0) THEN
             do i=1,-nonods
                ewrite(3,*)'i,ml(i),dm1(i),dm2(i):',&
                     &                 i,ml(i),dm1(i),dm2(i)
             END DO
             !     
            !     
            ewrite(3,*) 'here d5'
            IF(D3.EQ.0) THEN
                do J=1,FREDOP
                  !     DO 30 J=1,NNODPP
                  do I=FINDCT(J),FINDCT(J+1)-1 
                      COL=COLCT(I)
                      U(COL)=U(COL) + DT*C1T(I)*DP(J)/(ML(COL)+DM1(COL))
                      V(COL)=V(COL) + DT*C2T(I)*DP(J)/(ML(COL)+DM2(COL))
                  END DO
                END DO
            ELSE

                IF(GOTTRAF.EQ.1) THEN
                  ! GOT TRAFFIC...
                  do J=1,FREDOP
                      !     DO 32 J=1,NNODPP
                      do I=FINDCT(J),FINDCT(J+1)-1 
                        COL=COLCT(I)
                        U(COL)=U(COL)+DT*C1T(I)*DP(J)/MAX(1.E-10,(1.-VOLTRAF(COL))*(ML(COL)+DM1(COL)))
                        V(COL)=V(COL)+DT*C2T(I)*DP(J)/MAX(1.E-10,(1.-VOLTRAF(COL))*(ML(COL)+DM2(COL)))
                        W(COL)=W(COL)+DT*C3T(I)*DP(J)/MAX(1.E-10,(1.-VOLTRAF(COL))*(ML(COL)+DM3(COL)))
                      END DO
                  END DO

                  ewrite(3,*)  'AFTER THE TRAFFIC EQUATIONS' 

                  CALL PMINMX(VOLTRAF,nonods,'******VOLTRAF')
                  CALL PMINMX(U,nonods,'******U-VELOCITY')                   
                  CALL PMINMX(V,nonods,'******V-VELOCITY')                   
                  CALL PMINMX(W,nonods,'******W-VELOCITY')  

                ELSE
                
                    IF(NBUOY.NE.0) THEN
                      ALLOCATE(BOY_INV(9*NONODS))
                      DO NOD=1,NONODS
            ! Form A^{-1} 
            ! solve for inverse and put in BOY_INV:
                        CALL INVXXX( &
                &         BOY_ML(NOD)+DM1(NOD),BOY_ML(NOD+NONODS),           BOY_ML(NOD+2*NONODS), &
                &         BOY_ML(NOD+3*NONODS),BOY_ML(NOD+4*NONODS)+DM2(NOD),BOY_ML(NOD+5*NONODS), &
                &         BOY_ML(NOD+6*NONODS),BOY_ML(NOD+7*NONODS),         BOY_ML(NOD+8*NONODS)+DM3(NOD),&
                &         BOY_INV(NOD),         BOY_INV(NOD+NONODS  ),BOY_INV(NOD+2*NONODS),&
                &         BOY_INV(NOD+3*NONODS),BOY_INV(NOD+4*NONODS),BOY_INV(NOD+5*NONODS),&
                &         BOY_INV(NOD+6*NONODS),BOY_INV(NOD+7*NONODS),BOY_INV(NOD+8*NONODS))
                      END DO
                    ENDIF
                    DO J=1,FREDOP
                      DO I=FINDCT(J),FINDCT(J+1)-1 
                          COL=COLCT(I)
                      IF(NBUOY.NE.0) THEN
                      U(COL)=U(COL) + DT*(BOY_INV(COL)*C1T(I)         +BOY_INV(COL+NONODS)*C2T(I) &
                &                       +BOY_INV(COL+2*NONODS)*C3T(I))*DP(J)
                      V(COL)=V(COL) + DT*(BOY_INV(COL+3*NONODS)*C1T(I)+BOY_INV(COL+4*NONODS)*C2T(I)&
                &                       +BOY_INV(COL+5*NONODS)*C3T(I))*DP(J)
                      W(COL)=W(COL) + DT*(BOY_INV(COL+6*NONODS)*C1T(I)+BOY_INV(COL+7*NONODS)*C2T(I)&
                &                       +BOY_INV(COL+8*NONODS)*C3T(I))*DP(J)
                      ELSE
                          U(COL)=U(COL) + DT*C1T(I)*DP(J)/(ML(COL)+DM1(COL))
                          V(COL)=V(COL) + DT*C2T(I)*DP(J)/(ML(COL)+DM2(COL))
                          W(COL)=W(COL) + DT*C3T(I)*DP(J)/(ML(COL)+DM3(COL))
                      ENDIF
                      END DO
                    END DO       
                
                ENDIF

            ENDIF
             !     ENDOF ... IF(CONMAS.EQ.0) THEN
          ENDIF
          !     ewrite(3,*)'U:',U
          !     ewrite(3,*)'V:',V
          !     ewrite(3,*)'W:',W
          !     
          !     
          !     
          IF(PARA.EQ.1) CALL HALGET(U,NONODS,NONODS,NNODP,halo_tag)
          IF(PARA.EQ.1) CALL HALGET(V,NONODS,NONODS,NNODP,halo_tag)
          IF(D3.EQ.1) THEN 
             IF(PARA.EQ.1) CALL HALGET(W,NONODS,NONODS,NNODP,halo_tag)
          ENDIF
          !     Check 
          allocate(work(1:FREDOP))
          IF(NDWISP) THEN
            !     ADD R=KCMC*P 
            CALL MULMAT(work,P,KCMC,FINCMC,COLCMC,&
                  &                 NCMC,NCMC,MATSTR,FREDOP,FREDOP,NNODPP,.false.)
            do I=1,FREDOP
                work(I)=-DT*work(I)
            END DO
          ENDIF

          ! cjc
          ! The following line is not DG safe.
          CALL CALDIV(.not. NDWISP,work,NPHASE,&
              &                 FREDOP,NONODS,D3, &
              &                 U,V,W, &
              &                 C1T,C2T,C3T,FINDCT,COLCT,NCT, &
              !     Parallel ...
              &                 PARA,halo_tag,&
              &                 halo_tag_p,NNODP,NNODPP )
          deallocate(work)
          !     
          !     
          IF(FINISH) RETURN
          !     ENDDO OF IPLOO=1,NPLOO
       END DO
       !     
    END DO

    ewrite(1,*)'just leaving ctyexa:'
    RETURN
  END SUBROUTINE CTYEXA

  SUBROUTINE CALDIV(ZERO,RHS,NPHASE,&
       &     FREDOP,NONODS,D3, &
       &     UPHA,VPHA,WPHA, &
       &     C1TP,C2TP,C3TP,FINDCT,COLCT,NCT, &
       !     Parallel ...
       &     PARA,halo_tag,&
       &     halo_tag_p,NNODP,NNODPP )
    !     This sub calculates div q or sum of div q for multi-phase flow.
    !     then prints 1'st and second largest values of. 
    INTEGER NNODP, NNODPP
    INTEGER NONODS,FREDOP,D3,NPHASE,NCT
    !     Working arrays...
    REAL RHS(FREDOP)
    REAL UPHA(NPHASE*NONODS),VPHA(NPHASE*NONODS)
    REAL WPHA(NPHASE*NONODS) 
    LOGICAL ZERO
    REAL C1TP(NCT*NPHASE),C2TP(NCT*NPHASE),C3TP(NCT*NPHASE)
    INTEGER FINDCT(FREDOP+1),COLCT(NCT)
    INTEGER  PARA,halo_tag,halo_tag_p
    !     Local mem
    REAL DTONE, RMAXDI
    INTEGER IPHA, NPLUS, INCT, I, IIPR, IIPR2
    real, dimension(:), allocatable:: work

    allocate(work(1:FREDOP))
    
    DTONE=1.0

    IF(ZERO) RHS = 0.0

    do IPHA=1,NPHASE 
       NPLUS=1+(IPHA-1)*NONODS
       INCT =1+(IPHA-1)*NCT
       CALL CTMULT(work,UPHA(NPLUS),VPHA(NPLUS),WPHA(NPLUS), &
            &        DTONE, &
            &        .true.,0, FREDOP,NONODS,D3,  &
            &        C1TP(INCT),C2TP(INCT),C3TP(INCT),FINDCT,COLCT,NCT, &
            !     Parallel ...
            &        PARA,halo_tag,&
            &        halo_tag_p,NNODP,NNODPP )
       !     Multiply volume fraction by the result & put into RHS. 
       RHS=RHS+work(1:FREDOP)
    end do
      
    deallocate(work)

    RMAXDI=-1
    iipr=0
    do I=1,FREDOP
       IF(ABS(RHS(I)).GT.RMAXDI) then
          RMAXDI=ABS(RHS(I))
          iipr=i
       endif
    end do
    ewrite(2,*) 'MAX DIVQ in element=',RMAXDI
    ewrite(2,*) 'divq: WHICH OCCURS AT PRES NOD=',iipr
    ewrite(2,*) 'R2NORM(RHS,NNODPP,0):',R2NORM(RHS,NNODPP,0)
    !     
    RMAXDI=-1
    iipr2=0
    do  I=1,FREDOP! Was loop 1209
       IF((ABS(RHS(I)).GT.RMAXDI).AND.(I.NE.iipr)) then
          RMAXDI=ABS(RHS(I))
          iipr2=i
       endif
    end do ! Was loop 1209
    ewrite(3,*) 'SECOND LARGEST sum of DIVQ:'
    ewrite(3,*) 'MAX DIVQ in element=',RMAXDI
    ewrite(3,*) 'divq: WHICH OCCURS AT PRES NOD=',iipr2
    RETURN
  END SUBROUTINE CALDIV
  !     
  !     
  !     
  !     

end module ctyexa_module
