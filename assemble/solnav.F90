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

module module_solnav
    use FLDebug
    use geostrophic
    use Solvers
    use sparse_tools
    use spud
    use ctyexa_module
    use global_parameters, only: OPTION_PATH_LEN, phase2state_index
    use futils
    use solve_momentum_cg
    use AllSorts
    use getrdr_module
    use getrhs_module
    use getrhu_module
    use hart3d_allsorts
    use rotated_boundary_conditions_legacy
    use state_module
    use halos

    implicit none
    
  private
  public solnav

contains

  SUBROUTINE SOLNAV(U,V,W,NU,NV,NW,IGUESS,P,DP,DT,&
       &     NONODS,FREDOP,&
       !     the matrices  - It is assumed that if ROTAT then CIT's & CMC already 
       !     contain rotation matrices ie. 
       &     C1T,C2T,C3T,FINDCT,COLCT,NCT, &
       !     for 2-D flows C1T is put in C3T.  
       &     ML,BIGM1,BIGM2,BIGM3,&
       &     FINDRM,COLM,CENTRM,NCOLM,NBIGM, &
       !     BIGM1 is put in BIGM2,BIGM3 often.  NBIGM is the number of entries in 
       !     the number of entries in the discretised momentum equations. 
       !     if COUPLE (and d3=0 then NBIGM=3*NCOLM and d3=1 then NBIGM=6*NCOLM). 
       &     NDWISP,KCMC,&
       &     CMC,FINCMC,COLCMC,NCMC,&
       &     NCMCB,NPRESS,NPROPT,RADISO,ALFST2,SPRESS,&
       !     some vectors of length NONODS
       &     DM1,DM2,DM3,FX,FY,FZ,&
       !     DMI's are not used when UZAWA=1 and COUPLE=1 AND .NOT.ROTAT
       !     NB COUPLE=1 when ROTAT=.TRUE.
       !     vecs of length FREDOP, DIVQ vec=C1TU+C2TV+C3TW (=0if incompresible flow)
       &     RHS,&
       !     the governing parameters, either 0 or 1. 
       &     UZAWA,CONMAS,COUPLE,D2HALF,D3, &
       &     MULPA,CGSOLQ,GMRESQ,POISON,PROJEC, &
       !     NDPSET contains the pressure node to set to zero. 
       &     NDPSET,PARA,halo_tag,&
       &     halo_tag_p,NNODP,NNODPP,&
       !     The following is for rotations only ROTAT=.TRUE.
       &     ROTAT,NNODRO,NRTDR,   NODROT,&
       &     NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
       &     RTDR,ROTMOM,&
       !     The multiphase stuff...
       &     NPHASE,IPHASE,&
       &     C1TP,C2TP,C3TP, &
       &     NBUOY,BOY_ML,&
       ! FOr traffic -25April05
       &      VOLTRAF,UTRAF,VTRAF,WTRAF,GOTTRAF,&
       !     --------------------------------------START OF ADD BY CRGW 14/03/06
       !           SINGLE PHASE SOLID COMPRESSIBILITY:
       &      MKCOMP, DIVQS, &
       &      state)
    !     --------------------------------------END OF ADD BY CRGW 14/03/06
    !     
    !     IF PRESYM then just the upper triangular part of the 
    !     pressure matrix is stored. 
    !     
    !     THIS SUB SOLVE THE COUPLED MOMENTUM AND CONTINUITY EQUATIONS
    !     FOR U,V,(W) & P. 
    !     There are 3 pressure methods that can be used here 
    !     a projection method(PROJEC=1) a poissons pressure 
    !     equation method(POISON=1) (THESE 2 METHODS CAN BE USED TOGETHER)
    !     and UZAWA=1's pressure method, which in this implementation 
    !     also assumes that the momentum equation matrix(matrices) 
    !     is (are) symmetric. ie) advection is treated explicitly. 
    !     If any vector is not used it will contain vector DM1.
    !     D2HALF=1  ONLY when COUPLE=1 and we have a two and a half D simulation 
    !     in cylindrical coords. 

    INTEGER NBIGM,NDPSET,NNODP,NNODPP
    INTEGER NRTDR
    !     
    INTEGER NONODS,FREDOP
    INTEGER UZAWA,CONMAS,COUPLE,D2HALF,D3
    INTEGER MULPA,CGSOLQ,GMRESQ
    INTEGER POISON,PROJEC,NPHASE
    INTEGER NCT,NCOLM,NCMC,NCMCB,NPRESS,NPROPT,RADISO
    REAL ALFST2,SPRESS
    REAL U(NONODS),V(NONODS),W(NONODS)
    !     most up to date vels...
    REAL NU(NONODS),NV(NONODS),NW(NONODS)
    LOGICAL IGUESS
    REAL P(FREDOP*NPRESS)
    REAL DP(FREDOP*NPRESS)
    REAL DT 
    REAL C1T(NCT*NPRESS),C2T(NCT*NPRESS),C3T(NCT*NPRESS)
    REAL C1TP(NCT),C2TP(NCT),C3TP(NCT)
      INTEGER NBUOY
      REAL BOY_ML(NBUOY*9*NONODS)
    INTEGER FINDCT(FREDOP+1),COLCT(NCT) 
    REAL ML(NONODS)
    REAL BIGM1(NBIGM),BIGM2(NBIGM),BIGM3(NBIGM)
    INTEGER FINDRM(NONODS+1),COLM(NCOLM),CENTRM(NONODS)
    LOGICAL NDWISP
    REAL KCMC(NCMC)
    REAL CMC(NCMC)
    INTEGER FINCMC(FREDOP+1),COLCMC(NCMC)
    REAL DM1(NONODS),DM2(NONODS),DM3(NONODS)
    REAL FX(NONODS),FY(NONODS),FZ(NONODS)
    REAL RHS(FREDOP*NPRESS)
    !     For parallel solution...
    INTEGER PARA,halo_tag
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
    REAL RTDR(NRTDR)
    LOGICAL ROTAT
    !     RTDR contains RT D R   ^T where D contains DMI's and R is the rotation matrix. 
    LOGICAL TRANSP,NTRANS,LD3,ROTMOM,LROTAT
    !     For multiphase stuff ONLY...
    INTEGER IPHASE
    REAL, dimension(:), allocatable:: ukeep, vkeep, wkeep
    !     Extra stuff used in ALPDIS...

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !ccccccccccccccFOR TRAFFIC

    REAL VOLTRAF(nonods),UTRAF(nonods),VTRAF(nonods),WTRAF(nonods) 
    INTEGER GOTTRAF 
    !CCC END OF TRAF

    !     Local memory...
    REAL CHANGE
    INTEGER IEPRE
    INTEGER I,IIIROT,ILENG2
    REAL, ALLOCATABLE, DIMENSION(:)::ML_DIAG

    integer ileng
    real, dimension(:), allocatable:: v1l, v2l, v1, v2, v3
    
    character(len=OPTION_PATH_LEN) pressure_option_path
    character(len=OPTION_PATH_LEN) velocity_option_path
    type(csr_matrix) :: cmc_matrix
    type(csr_sparsity) :: cmc_sparsity
    integer :: istate
    type(halo_type) :: halo
    integer :: stat

    !     --------------------------------------START OF ADD BY CRGW 14/03/06
    !     SINGLE PHASE SOLID COMPRESSIBILITY:
    INTEGER MKCOMP
    REAL DIVQS(FREDOP)
    type(state_type), dimension(:), intent(inout) :: state

    INTEGER, parameter::MATSTR=0

    if (d3==0) then
       ileng=2*nonods
    else
       ileng=3*nonods
    end if
    allocate(v1l(1:ileng), v2l(1:ileng))
    allocate(v1(1:nonods), v2(1:nonods), v3(1:nonods))
    
    ALLOCATE(ML_DIAG(3*NONODS))

    istate = phase2state_index(IPHASE)
    ewrite(2,*) 'IPHASE, istate = ', IPHASE, istate

    ewrite(2,*) 'cjc r2norm(p)', sum(p(1:fredop)*p(1:fredop))
    !     
    !     RETURN IF ************************
    IF(IPHASE.NE.1) RETURN 
    !     
    !     so I'll assume we're in phase 1 then?
    !     anyway option_path to pressure field for solver options in options tree
    pressure_option_path='/material_phase['//int2str(istate-1)// &
         &   ']/scalar_field::Pressure'
    !     and for velocity     
    velocity_option_path='/material_phase['//int2str(istate-1)// &
         &   ']/vector_field::Velocity'

    ewrite(1,*)'JUST INSIDE SOLNAV D3=',D3
    ewrite(3,*)'IPHASE,NPHASE,rotat',&
         &     IPHASE,NPHASE,rotat

    call flcomms_update(halo_tag, FX, 1, 1, 0)
    call flcomms_update(halo_tag, FY, 1, 1, 0)
    IF(D3.EQ.1) THEN
          call flcomms_update(halo_tag, FZ, 1, 1, 0)
    ENDIF
    
    if (NPHASE>1) then
      FLAbort("Multi-phase (NPHASE>1) no longer supported!")
    endif

    !     
    !     
    LD3=.FALSE.
    IF(D3.EQ.1) LD3=.TRUE.
    !     
    IF(ROTAT) THEN
       TRANSP=.FALSE.
       NTRANS=.TRUE.
       !     Form the matrix R^T D R in which D is DM1 etc.
       CALL GETRDR(   NONODS, &
            &        ROTAT,NNODRO,NRTDR,   NODROT,&
            &        NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
            &        DM1,DM2,DM3, RTDR)
       ewrite(3,*)'INSIDE SOLNAV HERE1' 
       !     Rotate velocities q=R q - to rotated system(NB. UPHA=U etc). 
       CALL ROTUVW(U,V,W,&
               &           NONODS,&
               &           NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
               &           NODROT,NNODRO,LD3,NTRANS)
       ewrite(3,*)'INSIDE SOLNAV HERE2' 
    ENDIF
    !     
    !     MAKE SURE ENOUGH STORAGE IS SENT DOWN
    !     Memory needed for pressure determination...
    IIIROT=0
    IF(ROTAT) IIIROT=1 
    !     
    !     Memory needed for pressure CORRECTION determination...
    ILENG2=0
    IF(ROTAT) ILENG2=ILENG
    !     
    ewrite(3,*)'poison,projec,uzawa,mulpa,cgsolq,gmresq',&
         &     poison,projec,uzawa,mulpa,cgsolq,gmresq
    ewrite(3,*)'couple=',couple
    IF((POISON.EQ.1).OR.(UZAWA.EQ.1)) THEN
      !     Find RHS of poissons pressure equation. 
      !     
      IF(ROTAT.OR.((UZAWA.EQ.1).AND.(COUPLE.EQ.1))) THEN
         ewrite(1,*)'GOING INTO GETRHU --4'
         CALL GETRHU(UZAWA,RHS,C1T,C2T,C3T,&
                 &                 FINDCT,COLCT,NCT,&
                 &                 U,V,W,FX,FY,FZ,  &
                 &                 V1L,V2L,DT,&
                 &                 FINDRM,COLM,CENTRM,NCOLM, NCOLM,&
                 &                 NONODS,FREDOP,ILENG,D2HALF,&
                 &                 PARA,halo_tag,NNODP,NNODPP, &
                 &                 ML, &
                 !     The following is for rotations only ROTAT=.TRUE.
                 &                 ROTAT,NNODRO,NRTDR,   NODROT,&
                 &                 NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
                 &                 DM1,DM2,DM3, RTDR)
         ewrite_minmax(RHS)
         ewrite(2,*)  "R2NORM of RHS = ", r2norm(RHS, NNODP, 1), DT
         !     NB FX=VECX etc. 
         !     
         !     .... IF(ROTAT.OR.((UZAWA.EQ.1).AND.(COUPLE.EQ.1))) 
      ELSE
            
          !     
          !     
          ewrite(1,*)'GOING INTO GETRHS --1-8'
          IF(MKCOMP.EQ.3) THEN
             CALL GETRHS(RHS,C1TP,C2TP,C3TP,&
                  &                   FINDCT,COLCT,NCT,&
                  &                   U,V,W,DM1,DM2,DM3,FX,FY,FZ,  &
                  &                   V1,V2,DT,&
                  &                   ML,&
                  &                   FINDRM,COLM,CENTRM,NCOLM, &
                  &                   NONODS,FREDOP,&
                  &                   CONMAS,D3,&
                  &                   PARA,halo_tag,NNODP,NNODPP) 
          ELSE
             CALL GETRHS(RHS,C1T,C2T,C3T,&
                  &                   FINDCT,COLCT,NCT,&
                  &                   U,V,W,DM1,DM2,DM3,FX,FY,FZ,  &
                  &                   V1,V2,DT,&
                  &                   ML,&
                  &                   FINDRM,COLM,CENTRM,NCOLM, &
                  &                   NONODS,FREDOP,&
                  &                   CONMAS,D3,&
                  &                   PARA,halo_tag,NNODP,NNODPP) 
          ENDIF
          ewrite(3,*)'HERE'

          IF(MKCOMP.EQ.3) THEN
             do I = 1,FREDOP
                RHS(I) = RHS(I) + (DIVQS(I)/DT)
             ENDDO
          ENDIF

         !     
         !     .... IF(ROTAT.OR.((UZAWA.EQ.1).AND.(COUPLE.EQ.1))) ELSE...
      ENDIF
      !     
      ewrite_minmax(RHS)
      ewrite_minmax(CMC)
      ewrite_minmax(C1T)
      ewrite_minmax(C2T)
      ewrite_minmax(C3T)
      ewrite_minmax(DIVQS)
      ewrite_minmax(C1TP)
      ewrite_minmax(C2TP)
      ewrite_minmax(C3TP)

      ewrite(3,*)'After GETRHS  ---'

      !     Find the pressure P. 
      IIIROT=0
      IF(ROTAT) IIIROT=1 
      ewrite(3,*)'UZAWA,COUPLE,PARA,ILENG,NONODS,FREDOP:',&
           &           UZAWA,COUPLE,PARA,ILENG,NONODS,FREDOP
      !     
      LROTAT=ROTAT
      ewrite(3,*) '-NCMC,NCMCB,NPRESS,NPROPT,RADISO:',&
           &           NCMC,NCMCB,NPRESS,NPROPT,RADISO
      !     ewrite(3,*) 'CMC(1..100):',(CMC(I),I=1,100)

      ! Ensure that the ndpset node pressure really is zero (if present).
      if (NDPSET>0.and.NDPSET<=FREDOP) then
         P(NDPSET)=0
      end if

      !               IF(PREPRE.GT.1000) THEN
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
          call petsc_solve(p, cmc_matrix, rhs,&
                 &     option_path=pressure_option_path)
          call deallocate(cmc_matrix)
          ewrite(2,*) "PETSc has solved pressure."
#else
          FLAbort("Fluidity no longer works without PETSc")
#endif
    ENDIF  ! poison or uzawa

    CALL GETF12(FX,FY,FZ,P, &
           &           C1T,C2T,C3T,FINDCT,COLCT,NCT, &
           &           NONODS,FREDOP,D3, &
           &           halo_tag_p,NNODPP)

    !     Solve for velocities 

    !     UKEEP is the U vel at previous time level.
    allocate(ukeep(1:nonods), vkeep(1:nonods), wkeep(1:nonods))
    ukeep=u
    vkeep=v
    if (d3/=0) wkeep=w

    CALL GETVEL(U,V,W,NU,NV,NW,IGUESS,FX,FY,FZ,DT, &
       &           DM1,DM2,DM3,BIGM1,BIGM2,BIGM3,&
       &           FINDRM,COLM,CENTRM,NCOLM,NBIGM,&
       &           V1,V3,V1L,V2L,&
       !     This NPHASE is for no of phases. 
       &           ILENG, NPHASE, NONODS,D3,COUPLE,&
       &           GMRESQ,&
       &           PARA,halo_tag,NNODP, &
       !     The following is for rotations only ROTAT=.TRUE.
       &           ROTMOM,&
       &           velocity_option_path)

    CHANGE=0.0
    CALL COMPER(PARA,U,UKEEP,NONODS,CHANGE)
    CALL COMPER(PARA,V,VKEEP,NONODS,CHANGE)
    IF(D3.NE.0)CALL COMPER(PARA,W,WKEEP,NONODS,CHANGE)
    
    deallocate(ukeep, vkeep, wkeep)

    ewrite(1,*)'After GETVEL'
    IF(PROJEC.EQ.1) THEN
      !     Correct the velocity and pressure to satisfy continuity. 
      ewrite(1,*)'about to enter CTYEXA'
      ILENG2=0
      IF(ROTAT) ILENG2=ILENG

      !           ---------------------------------------ADDITION BY CRGW 19/04/06

      CALL CTYEXA(U,V,W,P,RHS,DP,&
              &              NDWISP,KCMC,&
              &              CMC,FINCMC,COLCMC,NCMC,&
              &              NCMCB,NPRESS,NPROPT,ALFST2,SPRESS, &
              &              C1T,C2T,C3T,FINDCT,COLCT,NCT, &
              &              ML,&
              &              NBUOY,BOY_ML, &
              &              DM1,DM2,DM3,DT,&
              &              V1L,V2L,&
              &              FREDOP,NONODS,ILENG2,&
              &              CONMAS,D3,&
              &              NDPSET,&
              &              PARA,halo_tag,&
              &              halo_tag_p,NNODP,NNODPP,&
              !     The following is for rotations only ROTAT=.TRUE.
              &              ROTAT,NNODRO, NODROT,&
              &              NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
              !     The multiphase stuff...
              &              NPHASE,IPHASE,&
              &              C1TP,C2TP,C3TP,&
              !     FOr traffic 
              &              VOLTRAF,UTRAF,VTRAF,WTRAF,GOTTRAF,            &
              !     --------------------------------------START OF ADD BY CRGW 20/03/06
              !           SINGLE PHASE SOLID COMPRESSIBILITY:
              &                 MKCOMP, DIVQS, &
              &                 state)
      ewrite(2,*) 'cjc after ctyexa- size(u)',sum(u(1:nonods)&
              &              *u(1:nonods))

      IEPRE=1

    ENDIF
    
    IF(ROTAT) THEN
       !     ewrite(3,*) 'before rotuvw v(31)=',v(31)
       !     Rotate velocities back q=R^T q - to unrotated system. 
       CALL ROTUVW(U,V,W, &
               &           NONODS,&
               &           NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
               &           NODROT,NNODRO,LD3,TRANSP)
    ENDIF
    
    deallocate(v1,v2,v3,v1l,v2l)
    
  END SUBROUTINE SOLNAV

  SUBROUTINE GETF12(F1,F2,F3,P, &
       C1T,C2T,C3T,FINDCT,COLCT,NCT, &
       NONODS,FREDOP,D3,&
       halo_tag_p,NNODPP)
    use flcomms_module
    IMPLICIT NONE
    INTEGER NCT
    INTEGER NONODS,FREDOP,D3
    REAL F1(NONODS),F2(NONODS),F3(NONODS)
    ! F1,F2,F3 contain the r.h.s. of the momentum equations. 
    REAL P(FREDOP)
    REAL C1T(NCT),C2T(NCT),C3T(NCT)
    INTEGER COL
    INTEGER I,J,FINDCT(FREDOP+1),COLCT(NCT)
    ! For parallel processing ...
    ! NNODPP=no of subdomain nodal pressure nodes excluding halo nodes.
    INTEGER halo_tag_p,NNODPP
    real, dimension(:), allocatable:: work
    ! This subroutine forms F1 and F2 for RHS of momentum equations
    ! Make sure F1=VECX and F2=VECY before going on
    
    !
    if (IsParallel()) then
      allocate(work(1:fredop))
      CALL HALGET(P,FREDOP,FREDOP,NNODPP,&
         halo_tag_p)
      deallocate(work)
    end if

    IF(D3.EQ.0) THEN
      do  J=1,FREDOP! Was loop 440
         !       DO 440 J=1,NNODPP
         do  I=FINDCT(J),FINDCT(J+1)-1! Was loop 440
            COL=COLCT(I)
            F1(COL)=F1(COL)+C1T(I)*P(J)
            F2(COL)=F2(COL)+C2T(I)*P(J)
         end do ! Was loop 440
      end do ! Was loop 440
    ELSE
      do  J=1,FREDOP! Was loop 450
         !          ewrite(3,*) 'j=',j
         !       DO 450 J=1,NNODPP
         do  I=FINDCT(J),FINDCT(J+1)-1! Was loop 450
            COL=COLCT(I)
            !           ewrite(3,*) 'i,c1t(i),c2t(i),c3t(i):',
            !     &              i,c1t(i),c2t(i),c3t(i)
            F1(COL)=F1(COL)+C1T(I)*P(J)
            F2(COL)=F2(COL)+C2T(I)*P(J)
            F3(COL)=F3(COL)+C3T(I)*P(J)
         end do ! Was loop 450
      end do ! Was loop 450
    ENDIF
    
  END SUBROUTINE GETF12

end module module_solnav
