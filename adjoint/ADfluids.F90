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
!    amcgsoftware@imperial.ac.uk
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
!
#include "fdebug.h"
  
  SUBROUTINE ADFLUIDS_BC(filename, filenamelen)
    
    ! ************************************************************************
    ! THE MAIN PROGRAM IS USED TO SIMUALATE 3D FLOW. WE USE THE ADJOINT METHOD
    ! TO ASSIMILATE THE OBSERVED DATA IN ORDER TO IMPROVE THE PREDICTION OF 
    ! THE SIMULATION.THE BASIC MODELS ARE THE NAVIER STOKES, RADIATION(WITH OR 
    ! WITHOUT EVENT) AND OR ADVECTION-DIFFUSION TYPE OF EQUATIONS.
    ! DO 1000 global nonlinaer conjugate gradient Loop
    ! 500 Call FLUIDS_FORWARD--the forward model and calculate the cost function
    !     if start the new search direction
    !       call FLUIDS_ADJOINT--the adjoint model and calculate the gradient
    !     endif
    !     Call NONLINEARCG (including linear search)
    !     if linearsearch isn't convergent then
    !     go to 500 (linear search loop)
    !     endif 
    !     if the nonlineaCG is convergent, then
    !        jump out of the global loop
    !     else 
    !        go to the next global loop staep
    !1000 continue
    ! A list of the extral arguments for the adjoint model
    !FUNCT           :Cost function
    !NOCVA0        :the original number of the control variables in data assimailation
    !NOCVA          :the number of the control variables in data assimailation
    !MAXREC       :the maximum size of the vector which us used to record the control variables
    !CVA                : the voctor storing the control variables
    !IMEM_CVA : the vector storing the relative mesh points of the control variables
    !G                       : the gradient of the cost function
    !MAXGRAD   :the maximum of the the gradient
    !GLOITS          :the maximum iteration number of nonlinear CG
    !LINITS           : the current linear search number
    !CVA_X,CVA_Y,CVA_Z: the postion values of the control varables along the x, y, z directions
    !NOBCU_F,NOBCV_F,NOBCW_F: the number of the boundary for velocities U,V,W
    !NTIME           : the number of the time steps
    !NOSNOD       : the number of positions where we have obervational data
    !NOSTIM        : the number of the time steps when we have the observation data.
    !STIME           : the vector storing the time when we will invert the observational data
    !UEXAC,VEXAC,WEXAC: the velocity from the obervational data in x, y, z directions
    !SX,SY,SZ      : the coordination values x, y, z of the observational data
    !ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z: the spans for Gaussian data srpreading
    !FILNAM_INPUT: the input file name
    !FILNAM_SOURCE------the name of observational data file 
    !FILNAM_F------ the forward data file name which is obained by running GEM 
    !FILNAM_AD------the adjoint data file name which is obained by running GEM 
    !NGLOITS ------ the maximum iteration number of nonlinear CG 
    !D3:  D3=.true. 3D, .false. 2D
    !ADMESH: admesh=.true. mesh adaptivity, .false. nonmesh-adaptivity
    ! ************************************************************************
    !
    use FLDebug
    use signals
    use spud
    use fluids_adjoint_module
    use fluids_forward_module
    use NONLINCG_module
    IMPLICIT NONE
    integer, intent(in) :: filenamelen
    character(len=filenamelen), intent(in) :: filename
    REAL EVESFL(1),EVENTD(1),EVENTT(1),EVENDE(1),EVENTV(1)
    REAL EVECO1(1),EVECO2(1),EVECO3(1),EVECO4(1)
    REAL TEMPR1(1),TEMPR2(1),TEMPR3(1),TEMPR4(1),BURNUP1(1)
    REAL TEMEVE(1)
    REAL::EVEDT=0
    REAL EVEUDT,EVEDTN,ETEMIN,TIME1,TIME2
    REAL EVGDIN,EDENIN,EDENGA,EGAMD2,EGAMD3
    REAL VOILIM,EVTIR1,EVTIR2
    INTEGER EVTIOP,EVSBCY,EVMXOP
    INTEGER VERSEV,NDIMI,NDSIEV(1)
    LOGICAL::EVEFIR=.FALSE.
    LOGICAL GOTGAS
    CHARACTER*240 FILEVE
    INTEGER::IEVENT=0, NENERG=0, NNNODS=0, NNNELE=0, NDELA2=0, NMATEV=0, MXNMTT=1
    INTEGER::MXNDEV=1
    
    INTEGER NRMEM,NIMEM
    REAL, ALLOCATABLE::RMEM(:)
    INTEGER, ALLOCATABLE::IMEM(:)
    
    ! The extral variables for adjoint model, nonlinear conjugate gradient...
    REAL ::FUNCT,MAXGRAD   ! The cost function and maximum gradient
    INTEGER  ::NOCVA        ! Number of the control variables
    INTEGER ::NOCVA0,NOCVA1       ! Intial guess for NOCVA size
    REAL,DIMENSION(:),ALLOCATABLE::G,CVA ! Nonlinear CG
    REAL,DIMENSION(:),ALLOCATABLE::OG,OCVA! Nonlinear CG
    REAL,DIMENSION(:),ALLOCATABLE::CVA_X,CVA_Y,CVA_Z  
    REAL,DIMENSION(:),ALLOCATABLE::OCVA_X,OCVA_Y,OCVA_Z  
    REAL,DIMENSION(:), ALLOCATABLE::OGOLD,OD,OCVAOLD
    INTEGER,DIMENSION(:),ALLOCATABLE::IMEM_CVA,IMEM_OCVA
    REAL,PARAMETER :: ALPHA_AD = -1000.0    ! for the descent gradient method
    INTEGER ::LINITS,GLOITS      ! nonlinear conjugate and linearsearch
    INTEGER ::NGLOITS    ! nonlinear conjugate and linearsearch
    !  REAL    ::GOLDR                      ! linearsearch 
    LOGICAL ::LINCONVER
    REAL,PARAMETER ::CONJCONVER=1E-5    
    INTEGER,PARAMETER ::MAXREC=1000000
    !!   INTEGER,PARAMETER ::MAXREC=5000000
    REAL     ::DT,LTIME
    !  INTEGER,PARAMETER  ::NTIME=INT(LTIME/DT)+1
    INTEGER    ::NTIME
    REAL,PARAMETER     ::GOLDR = 5.0    ! the gold ratio for the linearsearch
    INTEGER,DIMENSION(:),ALLOCATABLE ::NOBCU_F,NOBCV_F,NOBCW_F
    INTEGER,DIMENSION(:),ALLOCATABLE ::NOBCU_F0,NOBCV_F0,NOBCW_F0
    REAL,DIMENSION(:), ALLOCATABLE::GOLD,D,CVAOLD
    INTEGER ::NOSTIM,NOSNOD
    REAL,DIMENSION(:,:),ALLOCATABLE ::UEXAC,VEXAC,WEXAC
    REAL,DIMENSION(:),ALLOCATABLE::SX,SY,SZ,STIME
    REAL ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    ! The variables for the mesh adaptivity....
    LOGICAL  ::ADMESHH,D33
    CHARACTER*480 ::FILNAM_INPUT,FILNAM_SOURCE,FILNAM_F,FILNAM_AD
    
    
    ! Free surface
    INTEGER,DIMENSION(:),ALLOCATABLE ::NOBCH_F
    INTEGER,DIMENSION(:),ALLOCATABLE ::NOBCH_F0
    LOGICAL ::GOTFRE
    
    ! Model covariance
    
    REAL,DIMENSION(:),ALLOCATABLE::gcovari
    REAL,DIMENSION(:),ALLOCATABLE::ogcovari
    REAL      ::F_SMOOTHNESS
    SAVE      ::F_SMOOTHNESS
    
    ! Local ......
    INTEGER ::II,jj,kk,k0,k1                        ! for the loop
    INTEGER ::IO,N
    REAL    ::A,MAXGRAD_AVE,F2,F0,COF
    ! for headland and seamount:  REAL,PARAMETER::COF=1.0
    !**  REAL,PARAMETER::COF=1.0 ! for headland
    INTEGER ::STATUS
    REAL,DIMENSION(:),ALLOCATABLE::hcheck
    REAL FAI1,FAI,xc1
    REAL,PARAMETER ::PIE=3.141592654
    REAL FSZERO
    
    !**check  REAL,PARAMETER ::F0=5.83199752696105 !alpha1=100000,ispex=40000,fix-full
    
    
    !******************************************************************************************************
    !   START ADJOINT MODEL
    
    ! Establish signal handlers.
    !  call initialise_signals
    
    ! Set up the "big" memory
    CALL memory_size(NIMEM, NRMEM)
    ALLOCATE(RMEM(NRMEM))
    ALLOCATE(IMEM(NIMEM))
    
    IF( .not.ALLOCATED(RMEM) ) STOP 101
    IF( .not.ALLOCATED(IMEM) ) STOP 102
    
    
    !         ewrite(3,*) 'running code'
    ewrite(3,*)  'ENTER THE FILE NAME TO READ FROM'
    ewrite(3,*) 'NB IF THIS IS A PARALLEL SIMULATION'
    ewrite(3,*) 'THEN WE MUST ENTER THE NAME OF THE '
    ewrite(3,*) 'FILE WITHOUT ANY SUFIX'
    ewrite(3,*) '240 characters max & if it begins with /'
    ewrite(3,*) 'then assume chomplete path has been entered'
    ewrite(3,*) 'ENTER INPUT FILE NAME'
    ewrite(3,*) ' ****************************** '
    ewrite(3,*)  'FILNAM_SOURCE'
    ewrite(3,*) 'FILNAM_F'
    ewrite(3,*) 'FILNAM_AD'
    ewrite(3,*) 'NGLOITS'
    ewrite(3,*)  'DT,LTIME,NTIME=INT(LTIME/DT)+1'
    ewrite(3,*)  'NOSTIM,NOSNOD'
    ewrite(3,*)  'GOTFRE,ADMESHH,D33'
    ewrite(3,*) ' ****************************** '
    
    READ '(A480)', FILNAM_INPUT
    !  ewrite(3,*) 'ENTER SOURCE FILE NAME',
    !          READ '(A240)', FILNAM_SOURCE
    !  ewrite(3,*) 'ENTER FORWARD FILE NAME',
    !          READ '(A240)', FILNAM_F
    !  ewrite(3,*) 'ENTER ADJOINT FILE NAME',
    !          READ '(A240)', FILNAM_AD
    
    OPEN(3,FILE=FILNAM_INPUT,STATUS='OLD')
    READ(3,*) FILNAM_SOURCE
    ewrite(3,*)  'FILNAM_SOURCE',FILNAM_SOURCE
    READ(3,*) FILNAM_F
    ewrite(3,*) 'FILNAM_F',FILNAM_F
    READ(3,*) FILNAM_AD
    ewrite(3,*) 'FILNAM_AD',FILNAM_AD
    READ(3,*) NGLOITS
    ewrite(3,*) 'NGLOITS',NGLOITS
    READ(3,*) DT,LTIME,NTIME
    ewrite(3,*)  'DT,LTIME,NTIME=INT(LTIME/DT)+1',DT,LTIME,NTIME
    READ(3,*) NOSTIM,NOSNOD
    ewrite(3,*)  'NOSTIM,NOSNOD',NOSTIM,NOSNOD
    READ(3,*) GOTFRE,ADMESHH,D33
    ewrite(3,*)  'GOTFRE,ADMESHH,D33',GOTFRE,ADMESHH,D33
    READ(3,*) F0,COF
    ewrite(3,*)  'F0,COF',F0,COF
    READ(3,*) ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    ewrite(3,*)  'ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z',ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
    READ(3,*) FSZERO
    ewrite(3,*)  FSZERO
    
    CLOSE(3)
    
    
    ALLOCATE(NOBCU_F(NTIME))
    ALLOCATE(NOBCV_F(NTIME))
    IF(D33) ALLOCATE(NOBCW_F(NTIME))
    ALLOCATE(NOBCU_F0(NTIME))
    ALLOCATE(NOBCV_F0(NTIME))
    IF(D33) ALLOCATE(NOBCW_F0(NTIME))
    ALLOCATE(UEXAC(NOSNOD,NOSTIM))
    ALLOCATE(VEXAC(NOSNOD,NOSTIM))
    IF(D33) ALLOCATE(WEXAC(NOSNOD,NOSTIM))
    ALLOCATE(SX(NOSNOD))
    ALLOCATE(SY(NOSNOD))
    IF(D33) ALLOCATE(SZ(NOSNOD))
    ALLOCATE(STIME(NOSTIM))
    
    ! ******test free surface
    ALLOCATE(NOBCH_F(NTIME)) 
    ALLOCATE(NOBCH_F0(NTIME)) 
    NOBCH_F(1:NTIME) = 0
    NOBCH_F0(1:NTIME) = 0
    !****
    
    
    !NULLIFY(CVAOLD)
    !NULLIFY(GOLD)
    !NULLIFY(D)
    
    ! *********************************************
    ! Start the simulation
    
    FUNCT =0.0
    MAXGRAD = 0.0
    !  NOCVA1 =3000000
    !  NOCVA0 =3000000
    NOCVA1 =5000000
    NOCVA0 =5000000
    NOCVA = NOCVA1
    ALLOCATE(CVA(NOCVA))
    ALLOCATE(IMEM_CVA(NOCVA))  
    ALLOCATE(G(NOCVA))  
    ALLOCATE(CVA_X(NOCVA))
    ALLOCATE(CVA_Y(NOCVA))
    ALLOCATE(CVA_Z(NOCVA))
    
    G(1:NOCVA) = 0.0
    CVA(1:NOCVA) = 0.0
    CVA_X(1:NOCVA) = 0.0
    CVA_Y(1:NOCVA) = 0.0
    CVA_Z(1:NOCVA) = 0.0
    ! Model covariance 
    ALLOCATE(gcovari(NOCVA))
    gcovari(1:NOCVA) = 0.0
    
    DO II=1,NOCVA
       IMEM_CVA(II)=0
    ENDDO
    
    DO II =1,NTIME
       NOBCU_F(II)=0
       NOBCV_F(II)=0
       IF(D33) THEN
          NOBCW_F(II)=0
       ENDIF
    ENDDO
    DO II =1,NTIME
       NOBCU_F0(II)=0
       NOBCV_F0(II)=0
       IF(D33) THEN
          NOBCW_F0(II)=0
       ENDIF
    ENDDO
    
    
    
    
    !  CALL EXASOL(UEXAC,VEXAC,WEXAC,SX,SY,SZ,D3,                            &
    !                     STIME,LTIME,DT,NOSTIM,NOSNOD,FILNAM_SOURCE)
    
    ! other........
    OPEN(1,FILE='function.dat',STATUS='REPLACE')
    OPEN(2,FILE='gradient.dat',STATUS='REPLACE')
    OPEN(4,FILE='gradient-ave.dat',STATUS='REPLACE')
    CLOSE(1)
    CLOSE(2)
    CLOSE(4)
    
    ALLOCATE(CVAOLD(NOCVA),STAT=status)
    IF( ALLOCATED(CVAOLD) ) THEN
       CVAOLD(1:NOCVA) = 0.0
    ELSE
       ewrite(-1,*)  "ERROR: Failed to allocate CVAOLD"
       ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
       STOP
    END IF
    ALLOCATE(GOLD(NOCVA),STAT=status)  
    IF( ALLOCATED(GOLD) ) THEN
       GOLD(1:NOCVA) = 0.0
    ELSE
       ewrite(-1,*)  "ERROR: Failed to allocate GOLD"
       ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
       STOP
    END IF
    ALLOCATE(D(NOCVA),STAT=status)
    IF( ALLOCATED(D) ) THEN
       D(1:NOCVA) = 0.0
    ELSE
       ewrite(-1,*)  "ERROR: Failed to allocate D"
       ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
       STOP
    END IF
    
    LINITS = 1
    DO GLOITS=1,NGLOITS
       ewrite(3,*) '*********************'
       ewrite(3,*) 'gloits=',gloits  
       
       !    admesh = .true.
       
       IF(GLOITS.NE.1) THEN
          ! Rerun the forward model, to store X,Y,Z; U,V,W;CVA, CVAOLD,IMEM_CVA .....
          LINITS = 999999
          CALL FORWARDDATA_CONTROL_VARIABLES(CVA,IMEM_CVA,D,            &
               CVA_X,CVA_Y,CVA_Z,NOCVA,D33)
          ! Reallocate the larger memory to CVA, IMEM_CVA, .......
          ! Reallocate CVA
          ALLOCATE(OCVA(NOCVA),STAT=status)
          IF( .NOT.(ALLOCATED(OCVA)) ) THEN
             ewrite(-1,*)  "ERROR: Failed to allocate OCVA"
             ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
             STOP
          END IF
          
          OCVA(1:NOCVA) = CVA(1:NOCVA)
          DEALLOCATE(CVA)
          ALLOCATE(CVA(NOCVA1),STAT=status)
          IF( .NOT.(ALLOCATED(CVA)) ) THEN
             ewrite(-1,*)  "ERROR: Failed to allocate CVA"
             ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
             STOP
          END IF
          CVA(1:NOCVA1) = 0.0
          CVA(1:NOCVA) = OCVA(1:NOCVA)
          DEALLOCATE(OCVA)
          
          ! Reallocate IMEM_CVA
          ALLOCATE(IMEM_OCVA(NOCVA),STAT=status)
          IF( .NOT.(ALLOCATED(IMEM_OCVA)) ) THEN
             ewrite(-1,*)  "ERROR: Failed to allocate IMEM_OCVA"
             ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
             STOP
          END IF
          IMEM_OCVA(1:NOCVA) = IMEM_CVA(1:NOCVA)
          DEALLOCATE(IMEM_CVA)
          ALLOCATE(IMEM_CVA(NOCVA1),STAT=status)
          IF( .NOT.(ALLOCATED(IMEM_CVA)) ) THEN
             ewrite(-1,*)  "ERROR: Failed to allocate IMEM_CVA"
             ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
             STOP
          END IF
          DO II=1,NOCVA1
             IMEM_CVA(II)=0
          ENDDO
          IMEM_CVA(1:NOCVA) = IMEM_OCVA(1:NOCVA)
          DEALLOCATE(IMEM_OCVA)
          
          ! Reallocate G
          ALLOCATE(OG(NOCVA),STAT=status)
          IF( .NOT.(ALLOCATED(OG)) ) THEN
             ewrite(-1,*)  "ERROR: Failed to allocate OG"
             ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
             STOP
          END IF
          OG(1:NOCVA) = G(1:NOCVA)
          DEALLOCATE(G)
          ALLOCATE(G(NOCVA1),STAT=status)
          IF( .NOT.(ALLOCATED(G)) ) THEN
             ewrite(-1,*)  "ERROR: Failed to allocate G"
             ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
             STOP
          END IF
          G(1:NOCVA1) = 0.0
          G(1:NOCVA) = OG(1:NOCVA)
          DEALLOCATE(OG)
          
          ! Reallocate CVA_X
          ALLOCATE(OCVA_X(NOCVA),STAT=status)
          IF( .NOT.(ALLOCATED(OCVA_X)) ) THEN
             ewrite(-1,*)  "ERROR: Failed to allocate OCVA_X"
             ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
             STOP
          END IF
          OCVA_X(1:NOCVA) = CVA_X(1:NOCVA)
          DEALLOCATE(CVA_X)
          ALLOCATE(CVA_X(NOCVA1),STAT=status)
          IF( .NOT.(ALLOCATED(CVA_X)) ) THEN
             ewrite(-1,*)  "ERROR: Failed to allocate CVA_X"
             ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
             STOP
          END IF
          CVA_X(1:NOCVA1) = 0.0
          CVA_X(1:NOCVA) = OCVA_X(1:NOCVA)
          DEALLOCATE(OCVA_X)
          
          ! Reallocate CVA_Y
          ALLOCATE(OCVA_Y(NOCVA),STAT=status)
          IF( .NOT.(ALLOCATED(OCVA_Y)) ) THEN
             ewrite(-1,*)  "ERROR: Failed to allocate OCVA_Y"
             ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
             STOP
          END IF
          OCVA_Y(1:NOCVA) = CVA_Y(1:NOCVA)
          DEALLOCATE(CVA_Y)
          ALLOCATE(CVA_Y(NOCVA1),STAT=status)
          IF( .NOT.(ALLOCATED(CVA_Y)) ) THEN
             ewrite(-1,*)  "ERROR: Failed to allocate CVA_Y"
             ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
             STOP
          END IF
          CVA_Y(1:NOCVA1) = 0.0
          CVA_Y(1:NOCVA) = OCVA_Y(1:NOCVA)
          DEALLOCATE(OCVA_Y)
          
          ! Reallocate CVA_Z
          ALLOCATE(OCVA_Z(NOCVA),STAT=status)
          IF( .NOT.(ALLOCATED(OCVA_Z)) ) THEN
             ewrite(-1,*)  "ERROR: Failed to allocate OCVA_Z"
             ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
             STOP
          END IF
          OCVA_Z(1:NOCVA) = CVA_Z(1:NOCVA)
          DEALLOCATE(CVA_Z)
          ALLOCATE(CVA_Z(NOCVA1),STAT=status)
          IF( .NOT.(ALLOCATED(CVA_Z)) ) THEN
             ewrite(-1,*)  "ERROR: Failed to allocate CVA_Z"
             ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
             STOP
          END IF
          CVA_Z(1:NOCVA1) = 0.0
          CVA_Z(1:NOCVA) = OCVA_Z(1:NOCVA)
          DEALLOCATE(OCVA_Z)
          
          ! Reallocate gradient from Mocel covariance vector
          ALLOCATE(ogcovari(NOCVA),STAT=status)
          IF( .NOT.(ALLOCATED(ogcovari)) ) THEN
             ewrite(-1,*)  "ERROR: Failed to allocate gcovari"
             ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
             STOP
          END IF
          ogcovari(1:NOCVA) = gcovari(1:NOCVA)
          DEALLOCATE(gcovari)
          ALLOCATE(gcovari(NOCVA1),STAT=status)
          IF( .NOT.(ALLOCATED(gcovari)) ) THEN
             ewrite(-1,*)  "ERROR: Failed to allocate gcovari"
             ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
             STOP
          END IF
          gcovari(1:NOCVA1) = 0.0
          gcovari(1:NOCVA) = ogcovari(1:NOCVA)
          DEALLOCATE(ogcovari)
          
          !Reallocate CVAOLD, GOLD, D
          IF(ADMESHH) THEN
             ALLOCATE(OCVAOLD(NOCVA))
             OCVAOLD(1:NOCVA)=CVAOLD(1:NOCVA)
             DEALLOCATE(CVAOLD)
             ALLOCATE(CVAOLD(NOCVA1),STAT=status)
             IF( ALLOCATED(CVAOLD) ) THEN
                CVAOLD(1:NOCVA1) = 0.0
             ELSE
                ewrite(-1,*)  "ERROR: Failed to allocate CVAOLD"
                ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
                STOP
             END IF
             DEALLOCATE(OCVAOLD)
             
             
             ALLOCATE(OGOLD(NOCVA))
             GOLD(1:NOCVA)=OGOLD(1:NOCVA)
             DEALLOCATE(GOLD)
             ALLOCATE(GOLD(NOCVA1),STAT=status)  
             IF( ALLOCATED(GOLD) ) THEN
                GOLD(1:NOCVA1) = 0.0
             ELSE
                ewrite(-1,*)  "ERROR: Failed to allocate GOLD"
                ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
                STOP
             END IF
             DEALLOCATE(OGOLD)
             
             
             ALLOCATE(OD(NOCVA))
             D(1:NOCVA)=OD(1:NOCVA)
             DEALLOCATE(D)
             ALLOCATE(D(NOCVA1),STAT=status)
             IF( ALLOCATED(D) ) THEN
                D(1:NOCVA1) = 0.0
             ELSE
                ewrite(-1,*)  "ERROR: Failed to allocate D"
                ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
                STOP
             END IF
             DEALLOCATE(OD)
             
          ENDIF
          
          
          FILEVE=filename
          CALL cpu_time( TIME1 )
          
          NOCVA = NOCVA1
          
          CALL FLUIDS_FORWARD(IEVENT,EVESFL,EVENTD,EVENTT,EVENDE,EVENTV,       &
               EVECO1,EVECO2,EVECO3,EVECO4,                                    &
               TEMPR1,TEMPR2,TEMPR3,TEMPR4,BURNUP1,                           &
               EVEFIR,EVEDT,EVEUDT,EVEDTN, NENERG, NDELA2,                     &
               NNNODS,NNNELE,                                                  &
               FILEVE,EVGDIN,                                                  &
               EDENIN,ETEMIN, EDENGA,EGAMD2,EGAMD3,                            &
               VOILIM, EVTIR1,EVTIR2, EVTIOP,EVSBCY,EVMXOP,                    &
               NMATEV,MXNMTT,TEMEVE,                                           &
               NDIMI,MXNDEV,NDSIEV, VERSEV,GOTGAS,                             &
               RMEM, NRMEM, IMEM, NIMEM,                                       &
               ! New staff for the data assmilation...
               FUNCT,NOCVA0,NOCVA,MAXREC,CVA,IMEM_CVA,G,MAXGRAD,               &
               GLOITS,LINITS,                                                  &
               CVA_X,CVA_Y,CVA_Z,                                              &
               NOBCU_F,NOBCV_F,NOBCW_F,NTIME,D33,                               &
               NOBCU_F0,NOBCV_F0,NOBCW_F0,                                        &
               ADMESHH,NOSTIM,NOSNOD,STIME,                                     &
               UEXAC,VEXAC,WEXAC,                                              &
               SX,SY,SZ,                                                       &
               ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z,                                &
               FILNAM_INPUT,FILNAM_SOURCE,FILNAM_F,FILNAM_AD,                  &
               CVAOLD,GOLD,D,                                                  &
               NOBCH_F,NOBCH_F0,gcovari,F_SMOOTHNESS,xc1)
          
          CALL cpu_time( TIME2 )
          LINITS = 1
          F2 = FUNCT
          !*******new *****
          FUNCT = ABS(FUNCT-F0)
          !*********
          
          
          ! re-allocate the memory of the previous control variables, gradient
          ! and search direction
          ! Re-allocate exactly the Size of CVA and G
          NOCVA = 0
          
          DO II = 1,NTIME
             IF(GOTFRE) THEN
                NOCVA = NOCVA+NOBCH_F(II)
             ELSE
                NOCVA = NOCVA+NOBCU_F(II)+NOBCV_F(II)
                IF(D33) THEN
                   NOCVA = NOCVA+NOBCW_F(II)
                ENDIF
             ENDIF
          ENDDO
          
          DO II = 1,NTIME
             IF(GOTFRE) THEN
                NOBCH_F0(II) = NOBCH_F(II)
             ELSE
                NOBCU_F0(II) = NOBCU_F(II) 
                NOBCV_F0(II) = NOBCV_F(II) 
                IF(D33) NOBCW_F0(II) = NOBCW_F(II) 
             ENDIF
          ENDDO
          
          ewrite(3,*) 'NOCVA in MAIN',NOCVA
          
          
          ! Reallocate CVA
          ALLOCATE(OCVA(NOCVA),STAT=status)
          OCVA(1:NOCVA) = CVA(1:NOCVA)
          DEALLOCATE(CVA)
          ALLOCATE(CVA(NOCVA),STAT=status)
          CVA(1:NOCVA) = OCVA(1:NOCVA)
          DEALLOCATE(OCVA)
          ! Reallocate IMEM_CVA
          ALLOCATE(IMEM_OCVA(NOCVA),STAT=status)
          IMEM_OCVA(1:NOCVA) = IMEM_CVA(1:NOCVA)
          DEALLOCATE(IMEM_CVA)
          ALLOCATE(IMEM_CVA(NOCVA),STAT=status)
          IMEM_CVA(1:NOCVA) = IMEM_OCVA(1:NOCVA)
          DEALLOCATE(IMEM_OCVA)
          ! Reallocate G
          ALLOCATE(OG(NOCVA),STAT=status)
          OG(1:NOCVA) = G(1:NOCVA)
          DEALLOCATE(G)
          ALLOCATE(G(NOCVA),STAT=status)
          G(1:NOCVA) = OG(1:NOCVA)
          DEALLOCATE(OG)
          
          ! Reallocate CVA_X
          ALLOCATE(OCVA_X(NOCVA),STAT=status)
          OCVA_X(1:NOCVA) = CVA_X(1:NOCVA)
          DEALLOCATE(CVA_X)
          ALLOCATE(CVA_X(NOCVA),STAT=status)
          CVA_X(1:NOCVA) = OCVA_X(1:NOCVA)
          DEALLOCATE(OCVA_X)
          ! Reallocate CVA_Y
          ALLOCATE(OCVA_Y(NOCVA),STAT=status)
          OCVA_Y(1:NOCVA) = CVA_Y(1:NOCVA)
          DEALLOCATE(CVA_Y)
          ALLOCATE(CVA_Y(NOCVA),STAT=status)
          CVA_Y(1:NOCVA) = OCVA_Y(1:NOCVA)
          DEALLOCATE(OCVA_Y)
          ! Reallocate CVA_Z
          ALLOCATE(OCVA_Z(NOCVA),STAT=status)
          OCVA_Z(1:NOCVA) = CVA_Z(1:NOCVA)
          DEALLOCATE(CVA_Z)
          ALLOCATE(CVA_Z(NOCVA),STAT=status)
          CVA_Z(1:NOCVA) = OCVA_Z(1:NOCVA)
          DEALLOCATE(OCVA_Z)
          ! Reallocate gcovari
          ALLOCATE(ogcovari(NOCVA),STAT=status)
          ogcovari(1:NOCVA) = gcovari(1:NOCVA)
          DEALLOCATE(gcovari)
          ALLOCATE(gcovari(NOCVA),STAT=status)
          gcovari(1:NOCVA) = ogcovari(1:NOCVA)
          DEALLOCATE(ogcovari)
          
          !Reallocate CVAOLD, GOLD, D
          IF(ADMESHH) THEN
             
             ALLOCATE(OCVAOLD(NOCVA))
             OCVAOLD(1:NOCVA)=CVAOLD(1:NOCVA)
             DEALLOCATE(CVAOLD)
             ALLOCATE(CVAOLD(NOCVA),STAT=status)
             IF( ALLOCATED(CVAOLD) ) THEN
                CVAOLD(1:NOCVA) = 0.0
             ELSE
                ewrite(-1,*)  "ERROR: Failed to allocate CVAOLD"
                ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
                STOP
             END IF
             CVAOLD(1:NOCVA)=OCVAOLD(1:NOCVA)
             DEALLOCATE(OCVAOLD)
             
             
             ALLOCATE(OGOLD(NOCVA))
             OGOLD(1:NOCVA)=GOLD(1:NOCVA)
             DEALLOCATE(GOLD)
             ALLOCATE(GOLD(NOCVA),STAT=status)  
             IF( ALLOCATED(GOLD) ) THEN
                GOLD(1:NOCVA) = 0.0
             ELSE
                ewrite(-1,*)  "ERROR: Failed to allocate GOLD"
                ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
                STOP
             END IF
             GOLD(1:NOCVA)=OGOLD(1:NOCVA)
             DEALLOCATE(OGOLD)
             
             
             ALLOCATE(OD(NOCVA))
             OD(1:NOCVA)=D(1:NOCVA)
             DEALLOCATE(D)
             ALLOCATE(D(NOCVA),STAT=status)
             IF( ALLOCATED(D) ) THEN
                D(1:NOCVA) = 0.0
             ELSE
                ewrite(-1,*)  "ERROR: Failed to allocate D"
                ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
                STOP
             END IF
             D(1:NOCVA)=OD(1:NOCVA)
             DEALLOCATE(OD)
             
          ENDIF
          
          goto 600
       ENDIF  ! end if (gloits.ne.1)
       
       
       !     IF(GLOITS.NE.1.AND.LINITS.EQ.1) THEN
       !  GOTO 600
       !     ENDIF
       
       k1=1
       !    xc1=0.0
       xc1=-0.1
500    CONTINUE
       
       ! Solve the Navier Stoke eqn. and calculate the cost function..
       !--------------------------------------------------------------
       !This is a stub for F_FLUID without event option.
       
       
       FILEVE=filename
       CALL cpu_time( TIME1 )
       
       ewrite(3,*)  'before calling'
       ewrite(3,*) NTIME,FUNCT,NOCVA0,NOCVA,MAXREC
       ewrite(3,*) '11'
       ewrite(3,*) 'MAXGRAD',MAXGRAD
       ewrite(3,*) 'GLOITS,LINITS=',GLOITS,LINITS
       CALL FLUIDS_FORWARD(IEVENT,EVESFL,EVENTD,EVENTT,EVENDE,EVENTV,       &
            EVECO1,EVECO2,EVECO3,EVECO4,                                    &
            TEMPR1,TEMPR2,TEMPR3,TEMPR4,BURNUP1,                           &
            EVEFIR,EVEDT,EVEUDT,EVEDTN, NENERG, NDELA2,                     &
            NNNODS,NNNELE,                                                  &
            FILEVE,EVGDIN,                                                  &
            EDENIN,ETEMIN, EDENGA,EGAMD2,EGAMD3,                            &
            VOILIM, EVTIR1,EVTIR2, EVTIOP,EVSBCY,EVMXOP,                    &
            NMATEV,MXNMTT,TEMEVE,                                           &
            NDIMI,MXNDEV,NDSIEV, VERSEV,GOTGAS,                             &
            RMEM, NRMEM, IMEM, NIMEM,                                       &
            ! New staff for the data assmilation...
            FUNCT,NOCVA0,NOCVA,MAXREC,CVA,IMEM_CVA,G,MAXGRAD,               &
            GLOITS,LINITS,                                                  &
            CVA_X,CVA_Y,CVA_Z,                                              &
            NOBCU_F,NOBCV_F,NOBCW_F,NTIME,D33,                               &
            NOBCU_F0,NOBCV_F0,NOBCW_F0,                                        &
            ADMESHH,NOSTIM,NOSNOD,STIME,                                     &
            UEXAC,VEXAC,WEXAC,                                              &
            SX,SY,SZ,                                                       &
            ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z,                                &
            FILNAM_INPUT,FILNAM_SOURCE,FILNAM_F,FILNAM_AD,                  &
            CVAOLD,GOLD,D,                                                  &
            NOBCH_F,NOBCH_F0,gcovari,F_SMOOTHNESS,xc1)
       ewrite(3,*) 'funct=',funct
       
       ewrite(3,*) 'out of forward fluids'
       CALL cpu_time( TIME2 )
       
       ewrite(3,*)  'start,end,total cpu time',TIME1,TIME2,TIME2-TIME1
       
       
       F2 = FUNCT
       
       !*******new *****
       FUNCT = ABS(FUNCT-F0)
       !*********
       ewrite(3,*) 'funct=',funct
       !    stop 0
       !
       !
600    CONTINUE
       !     IF(GLOITS.EQ.1.AND.LINITS.EQ.1) THEN
       !        NOCVA = NOCVA0   ! roughly give the number of the control variables
       !     ENDIF
       
       
       IF(GLOITS.EQ.1.AND.LINITS.EQ.1) THEN
          !      IF(GLOITS.EQ.1.AND.LINITS.EQ.1) THEN
          ewrite(3,*) 'NOCVA in MAIN',NOCVA
          ! Re-allocate exactly the Size of CVA and G
          NOCVA = 0
          IF(GOTFRE) THEN
             DO II = 1,NTIME
                NOCVA = NOCVA+NOBCH_F(II)
             ENDDO
          ELSE
             DO II = 1,NTIME
                NOCVA = NOCVA+NOBCU_F(II)+NOBCV_F(II)
                IF(D33) THEN
                   NOCVA = NOCVA+NOBCW_F(II)
                ENDIF
             ENDDO
          ENDIF
          DO II = 1,NTIME
             IF(GOTFRE) THEN
                !free surface
                NOBCH_F0(II) = NOBCH_F(II) 
             ELSE
                NOBCU_F0(II) = NOBCU_F(II) 
                NOBCV_F0(II) = NOBCV_F(II) 
                IF(D33) NOBCW_F0(II) = NOBCW_F(II) 
             ENDIF
          ENDDO
          ewrite(3,*) 'NOBCU_F',NOBCU_F
          ewrite(3,*) 'NOBCV_F',NOBCV_F
          ewrite(3,*) 'NOBCW_F',NOBCW_F
          ewrite(3,*) 'NOCVA in ADmainfl',NOCVA

          open(1,file='cva0.dat',status='replace')
          write(1,*) 'nocva=',nocva
          write(1,*) (cva(ii),ii=1,nocva)
          close(1)
          open(1,file='imem_cva0.dat',status='replace')
          write(1,*) 'nocva=',nocva
          write(1,*) (imem_cva(ii),ii=1,nocva)
          close(1)
          
          ! Reallocate CVA
          ALLOCATE(OCVA(NOCVA),STAT=status)
          OCVA(1:NOCVA) = CVA(1:NOCVA)
          DEALLOCATE(CVA)
          ALLOCATE(CVA(NOCVA),STAT=status)
          CVA(1:NOCVA) = OCVA(1:NOCVA)
          DEALLOCATE(OCVA)
          ! Reallocate IMEM_CVA
          ALLOCATE(IMEM_OCVA(NOCVA),STAT=status)
          IMEM_OCVA(1:NOCVA) = IMEM_CVA(1:NOCVA)
          DEALLOCATE(IMEM_CVA)
          ALLOCATE(IMEM_CVA(NOCVA),STAT=status)
          IMEM_CVA(1:NOCVA) = IMEM_OCVA(1:NOCVA)
          DEALLOCATE(IMEM_OCVA)
          ! Reallocate G
          ALLOCATE(OG(NOCVA),STAT=status)
          OG(1:NOCVA) = G(1:NOCVA)
          DEALLOCATE(G)
          ALLOCATE(G(NOCVA),STAT=status)
          G(1:NOCVA) = OG(1:NOCVA)
          DEALLOCATE(OG)
          
          ! Reallocate CVA_X
          ALLOCATE(OCVA_X(NOCVA),STAT=status)
          OCVA_X(1:NOCVA) = CVA_X(1:NOCVA)
          DEALLOCATE(CVA_X)
          ALLOCATE(CVA_X(NOCVA),STAT=status)
          CVA_X(1:NOCVA) = OCVA_X(1:NOCVA)
          DEALLOCATE(OCVA_X)
          ! Reallocate CVA_Y
          ALLOCATE(OCVA_Y(NOCVA),STAT=status)
          OCVA_Y(1:NOCVA) = CVA_Y(1:NOCVA)
          DEALLOCATE(CVA_Y)
          ALLOCATE(CVA_Y(NOCVA),STAT=status)
          CVA_Y(1:NOCVA) = OCVA_Y(1:NOCVA)
          DEALLOCATE(OCVA_Y)
          ! Reallocate CVA_Z
          ALLOCATE(OCVA_Z(NOCVA),STAT=status)
          OCVA_Z(1:NOCVA) = CVA_Z(1:NOCVA)
          DEALLOCATE(CVA_Z)
          ALLOCATE(CVA_Z(NOCVA),STAT=status)
          CVA_Z(1:NOCVA) = OCVA_Z(1:NOCVA)
          DEALLOCATE(OCVA_Z)
          ! Reallocate gcovari
          ALLOCATE(ogcovari(NOCVA),STAT=status)
          ogcovari(1:NOCVA) = gcovari(1:NOCVA)
          DEALLOCATE(gcovari)
          ALLOCATE(gcovari(NOCVA),STAT=status)
          gcovari(1:NOCVA) = ogcovari(1:NOCVA)
          DEALLOCATE(ogcovari)
          
          !Reallocate CVAOLD, GOLD, D
          IF(ADMESHH) THEN
             
             ALLOCATE(OCVAOLD(NOCVA))
             OCVAOLD(1:NOCVA)=CVAOLD(1:NOCVA)
             DEALLOCATE(CVAOLD)
             ALLOCATE(CVAOLD(NOCVA),STAT=status)
             IF( ALLOCATED(CVAOLD) ) THEN
                CVAOLD(1:NOCVA) = 0.0
             ELSE
                ewrite(-1,*)  "ERROR: Failed to allocate CVAOLD"
                ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
                STOP
             END IF
             CVAOLD(1:NOCVA)=OCVAOLD(1:NOCVA)
             DEALLOCATE(OCVAOLD)
             
             
             ALLOCATE(OGOLD(NOCVA))
             OGOLD(1:NOCVA)=GOLD(1:NOCVA)
             DEALLOCATE(GOLD)
             ALLOCATE(GOLD(NOCVA),STAT=status)  
             IF( ALLOCATED(GOLD) ) THEN
                GOLD(1:NOCVA) = 0.0
             ELSE
                ewrite(-1,*)  "ERROR: Failed to allocate GOLD"
                ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
                STOP
             END IF
             GOLD(1:NOCVA)=OGOLD(1:NOCVA)
             DEALLOCATE(OGOLD)
             
             
             ALLOCATE(OD(NOCVA))
             OD(1:NOCVA)=D(1:NOCVA)
             DEALLOCATE(D)
             ALLOCATE(D(NOCVA),STAT=status)
             IF( ALLOCATED(D) ) THEN
                D(1:NOCVA) = 0.0
             ELSE
                ewrite(-1,*)  "ERROR: Failed to allocate D"
                ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
                STOP
             END IF
             D(1:NOCVA)=OD(1:NOCVA)
             DEALLOCATE(OD)
          ENDIF
       ENDIF
       
       IF(LINITS.EQ.1) THEN
          G(1:NOCVA) = 0.0
          
          ! Solve adjoint eqn, 
          ! and calculate the gradient & save the gradient in the file: gradient_bc.dat ...
          !--------------------------------------------------------------------------
          CALL cpu_time( TIME1 )
          ewrite(3,*) 'Running the adjoint model'
          !    ADMESH=.FALSE. 
          ewrite(3,*) 'admesh before adjoint in main',admeshh
          !     open(1,file='cva2.dat')
          !       write(1,*) 'nocva=',nocva
          !       write(1,*) (cva(ii),ii=1,nocva)
          !     close(1)
          !     open(1,file='imem_cva2.dat')
          !       write(1,*) 'nocva=',nocva
          !       write(1,*) (imem_cva(ii),ii=1,nocva)
          !     close(1)
          
          ewrite(3,*) 'funct=', funct
          !        stop 10
          
          CALL FLUIDS_ADJOINT(IEVENT,EVESFL,EVENTD,EVENTT,EVENDE,EVENTV,       &
               EVECO1,EVECO2,EVECO3,EVECO4,                                    &
               TEMPR1,TEMPR2,TEMPR3,TEMPR4,BURNUP1,                            &
               EVEFIR,EVEDT,EVEUDT,EVEDTN, NENERG, NDELA2,                     &
               NNNODS,NNNELE,                                                  &
               FILEVE,EVGDIN,                                                  &
               EDENIN,ETEMIN, EDENGA,EGAMD2,EGAMD3,                            &
               VOILIM, EVTIR1,EVTIR2, EVTIOP,EVSBCY,EVMXOP,                    &
               NMATEV,MXNMTT,TEMEVE,                                           &
               NDIMI,MXNDEV,NDSIEV, VERSEV,GOTGAS,                             &
               RMEM, NRMEM, IMEM, NIMEM,                                       &
               ! New staff for the data assmilation...
               FUNCT,NOCVA0,NOCVA,MAXREC,CVA,IMEM_CVA,G,MAXGRAD,              &
               GLOITS,LINITS,                                                 &
               CVA_X,CVA_Y,CVA_Z,                                             &
               NOBCU_F,NOBCV_F,NOBCW_F,NTIME,D33,                              &
               NOBCU_F0,NOBCV_F0,NOBCW_F0,                                        &
               ADMESHH,NOSTIM,NOSNOD,STIME,                                    &
               UEXAC,VEXAC,WEXAC,                                             &
               SX,SY,SZ,                                                      &
               ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z,                               &
               FILNAM_INPUT,FILNAM_SOURCE,FILNAM_F,FILNAM_AD,                 &
               CVAOLD,GOLD,D,                                                 &
               NOBCH_F,NOBCH_F0,gcovari,F_SMOOTHNESS)
          
          ewrite(3,*) 'out of adjoint fluids'
          CALL cpu_time( TIME2 )
          ewrite(3,*)  'start,end,total cpu time',TIME1,TIME2,TIME2-TIME1
          
          !      ADMESH=.TRUE.
          ! Adjust the control varaibles by using Nonlinear CG or descent gradient
          !-----------------------------------------------------------------------
          
          ewrite(3,*) 'MAXGRAD',MAXGRAD
          ewrite(3,*) 'funct',funct
          !             stop 11
          
          !***FOR CHECK THE CORRECTNESS of THE ADJOINT MODEL********
          
          if(.false.) then
             
             allocate(hcheck(nocva)) 
             hcheck(1:nocva) = 0.0
             kk=0
             DO II = 1, NTIME
                DO JJ = 1,NOBCH_F(II)
                   HCHECK(Kk+jj) = xc1*sin(2.0*pie*(ii-1)*DT/43200.0)/29.41
                ENDDO
                kk=kk+NOBCH_F(II)
             ENDDO
             
             FAI1 =0.0
             FAI=0.0
             DO II = 1, NOCVA
                FAI1=FAI1+HCHECK(II)*G(II)
             ENDDO
             FAI1 =0.0
             FAI=0.0
             jj=0
             DO II = 1, NOCVA
                if(abs(G(ii)).gt.1.0e-6) then
                   jj=jj+1
                   FAI1=FAI1+HCHECK(II)*G(II)
                endif
             ENDDO
             !    FAI = ABS(FUNCT)/FAI1
             FAI = FUNCT*sqrt(real(jj))/FAI1
             
             open(1,file='gc.dat')
             write(1,*) (g(II),ii=1,nocva)
             close(1)
             
             open(1,file='hc.dat')
             write(1,*) (hcheck(II),ii=1,nocva)
             close(1)
             deallocate(hcheck)
             
             open(1,file='gcorrect.dat',position='append')
             write(1,*) xc1,FAI,FAI1,F2,FUNCT
             close(1)
             
             ewrite(3,*) 'FAI,FAI1,FUNCT',FAI,FAI1,FUNCT
             xc1=xc1*0.1
             if(abs(xc1).le.1.0e-15) then
                stop 2345
             else
                go to 500
             endif

          endif
          !****CHECK THE CORRECTNESS of THE ADJOINT MODEL********
          
          
       ENDIF
       
       !      if(K1.EQ.1) THEN
       !       linits = 1
       !       k1=k1+1
       !       go to 500
       !      endif
       
       ! If using the descent gradient method, then
       
       !   GOTO 700
       
       ! Start the Nonlinear CG 
       !.......................
       
       ! Call the nonlinear conjugate gradient....
       ewrite(3,*) 'NONLINCG'
       IF(LINITS.EQ.1) THEN
          DO II = 1,NOCVA
             !            G(II)=1.0*G(II)
             G(II)=-COF*G(II)
             !                 G(II)=-10.0*G(II)
             !                 G(II)=-100.0*G(II)
          ENDDO
       ENDIF
       IF(GLOITS.eq.1.and.LINITS.EQ.1) THEN
          ewrite(3,*) 'VALFUN=FUNCT',FUNCT
       ENDIF
       
       IF(LINITS.EQ.1.AND.GLOITS.NE.1) THEN
          !     IF(LINITS.EQ.999999.AND.GLOITS.NE.1) THEN
          DO II = 1,NOCVA
             !            G(II)=1.0*G(II)
             GOLD(II)=-COF*GOLD(II)
             !                      GOLD(II)=-10.0*GOLD(II)
             !                 GOLD(II)=-100.0*GOLD(II)
          ENDDO
       ENDIF
       
       
  
       CALL NONLINCG(GLOITS,LINITS,FUNCT,CVA,G,                      &
            GOLDR,LINCONVER,NOCVA,GOLD,D,CVAOLD,FSZERO)

       
       ! Is the linear search convergent?
       IF(.NOT.LINCONVER) THEN
          OPEN(1,FILE='function.dat',POSITION='APPEND')
          WRITE(1,*) GLOITS,LINITS,FUNCT,F2
          CLOSE(1)
          GOTO 500
       ENDIF
       
       MAXGRAD_AVE=0.0
       DO II = 1,NOCVA
          MAXGRAD_AVE=MAXGRAD_AVE+G(II)/COF
       ENDDO
       MAXGRAD_AVE=MAXGRAD_AVE/NOCVA
       
       ! Is the nonlinear conjugate convergent?
       OPEN(1,FILE='function.dat',POSITION='APPEND')
       OPEN(2,FILE='gradient.dat',POSITION='APPEND')
       OPEN(4,FILE='gradient-ave.dat',POSITION='APPEND')
       WRITE(1,*) GLOITS,FUNCT
       WRITE(2,*) GLOITS,MAXGRAD
       WRITE(4,*) GLOITS,MAXGRAD_AVE
       CLOSE(1)
       CLOSE(2)
       CLOSE(4)
       ewrite(3,*) 'GLOITS,FUNCT',GLOITS,FUNCT
       ewrite(3,*) 'GLOITS,MAXGRAD',GLOITS,MAXGRAD
       IF(COF*MAXGRAD.LT.CONJCONVER) GOTO 1500
       
    ENDDO
    
150 FORMAT(1X,I4,1X,F10.3)
1500 ewrite(3,*) 'NONLINEAR CONJUGATE GREDIANT CONVERGENT'
    
    
    DEALLOCATE(RMEM)
    DEALLOCATE(IMEM)
    DEALLOCATE(CVA)
    DEALLOCATE(CVA_X)
    DEALLOCATE(CVA_Y)
    IF(D33) DEALLOCATE(CVA_Z)
    DEALLOCATE(gcovari)
    DEALLOCATE(IMEM_CVA)
    DEALLOCATE(G)  
    DEALLOCATE(CVAOLD)
    DEALLOCATE(GOLD)  
    DEALLOCATE(D)  
    DEALLOCATE(NOBCU_F)
    DEALLOCATE(NOBCV_F)
    IF(D33) DEALLOCATE(NOBCW_F)
    DEALLOCATE(NOBCU_F0)
    DEALLOCATE(NOBCV_F0)
    IF(D33) DEALLOCATE(NOBCW_F0)
    DEALLOCATE(UEXAC)
    DEALLOCATE(VEXAC)
    IF(D33) DEALLOCATE(WEXAC)
    DEALLOCATE(SX)
    DEALLOCATE(SY)
    IF(D33) DEALLOCATE(SZ)
    DEALLOCATE(STIME)
    
    DEALLOCATE(NOBCH_F)
    DEALLOCATE(NOBCH_F0)
    
    ! End the whole simulation
    !************************************************
    
    
    RETURN
  END SUBROUTINE ADFLUIDS_BC

   SUBROUTINE ADFLUIDS_BC_TEST(filename, filenamelen)
     
     ! ************************************************************************
     ! THE MAIN PROGRAM IS USED TO SIMUALATE 3D FLOW. WE USE THE ADJOINT METHOD
     ! TO ASSIMILATE THE OBSERVED DATA IN ORDER TO IMPROVE THE PREDICTION OF 
     ! THE SIMULATION.THE BASIC MODELS ARE THE NAVIER STOKES, RADIATION(WITH OR 
     ! WITHOUT EVENT) AND OR ADVECTION-DIFFUSION TYPE OF EQUATIONS.
     ! DO 1000 global nonlinaer conjugate gradient Loop
     ! 500 Call FLUIDS_FORWARD--the forward model and calculate the cost function
     !     if start the new search direction
     !       call FLUIDS_ADJOINT--the adjoint model and calculate the gradient
     !     endif
     !     Call NONLINEARCG (including linear search)
     !     if linearsearch isn't convergent then
     !     go to 500 (linear search loop)
     !     endif 
     !     if the nonlineaCG is convergent, then
     !        jump out of the global loop
     !     else 
     !        go to the next global loop staep
     !1000 continue
     ! A list of the extral arguments for the adjoint model
     !FUNCT           :Cost function
     !NOCVA0        :the original number of the control variables in data assimailation
     !NOCVA          :the number of the control variables in data assimailation
     !MAXREC       :the maximum size of the vector which us used to record the control variables
     !CVA                : the voctor storing the control variables
     !IMEM_CVA : the vector storing the relative mesh points of the control variables
     !G                       : the gradient of the cost function
     !MAXGRAD   :the maximum of the the gradient
     !GLOITS          :the maximum iteration number of nonlinear CG
     !LINITS           : the current linear search number
     !CVA_X,CVA_Y,CVA_Z: the postion values of the control varables along the x, y, z directions
     !NOBCU_F,NOBCV_F,NOBCW_F: the number of the boundary for velocities U,V,W
     !NTIME           : the number of the time steps
     !NOSNOD       : the number of positions where we have obervational data
     !NOSTIM        : the number of the time steps when we have the observation data.
     !STIME           : the vector storing the time when we will invert the observational data
     !UEXAC,VEXAC,WEXAC: the velocity from the obervational data in x, y, z directions
     !SX,SY,SZ      : the coordination values x, y, z of the observational data
     !ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z: the spans for Gaussian data srpreading
     !FILNAM_INPUT: the input file name
     !FILNAM_SOURCE------the name of observational data file 
     !FILNAM_F------ the forward data file name which is obained by running GEM 
     !FILNAM_AD------the adjoint data file name which is obained by running GEM 
     !NGLOITS ------ the maximum iteration number of nonlinear CG 
     !D3:  D3=.true. 3D, .false. 2D
     !ADMESH: admesh=.true. mesh adaptivity, .false. nonmesh-adaptivity
     ! ************************************************************************
     !
     use FLDebug
     use signals
     use spud
     use fluids_adjoint_module
     use fluids_forward_module
     IMPLICIT NONE
     integer, intent(in) :: filenamelen
     character(len=filenamelen), intent(in) :: filename
     REAL EVESFL(1),EVENTD(1),EVENTT(1),EVENDE(1),EVENTV(1)
     REAL EVECO1(1),EVECO2(1),EVECO3(1),EVECO4(1)
     REAL TEMPR1(1),TEMPR2(1),TEMPR3(1),TEMPR4(1),BURNUP1(1)
     REAL TEMEVE(1)
     REAL::EVEDT=0
     REAL EVEUDT,EVEDTN,ETEMIN,TIME1,TIME2
     REAL EVGDIN,EDENIN,EDENGA,EGAMD2,EGAMD3
     REAL VOILIM,EVTIR1,EVTIR2
     INTEGER EVTIOP,EVSBCY,EVMXOP
     INTEGER VERSEV,NDIMI,NDSIEV(1)
     LOGICAL::EVEFIR=.FALSE.
     LOGICAL GOTGAS
     CHARACTER*240 FILEVE
     INTEGER::IEVENT=0, NENERG=0, NNNODS=0, NNNELE=0, NDELA2=0, NMATEV=0, MXNMTT=1
     INTEGER::MXNDEV=1
     
     INTEGER NRMEM,NIMEM
     REAL, ALLOCATABLE::RMEM(:)
     INTEGER, ALLOCATABLE::IMEM(:)
     
     ! The extral variables for adjoint model, nonlinear conjugate gradient...
     REAL ::FUNCT,MAXGRAD   ! The cost function and maximum gradient
     INTEGER  ::NOCVA        ! Number of the control variables
     INTEGER ::NOCVA0,NOCVA1       ! Intial guess for NOCVA size
     REAL,DIMENSION(:),ALLOCATABLE::G,CVA ! Nonlinear CG
     REAL,DIMENSION(:),ALLOCATABLE::OG,OCVA! Nonlinear CG
     REAL,DIMENSION(:),ALLOCATABLE::CVA_X,CVA_Y,CVA_Z  
     REAL,DIMENSION(:),ALLOCATABLE::OCVA_X,OCVA_Y,OCVA_Z  
     REAL,DIMENSION(:), ALLOCATABLE::OGOLD,OD,OCVAOLD
     INTEGER,DIMENSION(:),ALLOCATABLE::IMEM_CVA,IMEM_OCVA
     REAL,PARAMETER :: ALPHA_AD = -1000.0    ! for the descent gradient method
     INTEGER ::LINITS,GLOITS      ! nonlinear conjugate and linearsearch
     INTEGER ::NGLOITS    ! nonlinear conjugate and linearsearch
     !  REAL    ::GOLDR                      ! linearsearch 
     LOGICAL ::LINCONVER
     REAL,PARAMETER ::CONJCONVER=1E-5    
     INTEGER,PARAMETER ::MAXREC=1000000
     !!   INTEGER,PARAMETER ::MAXREC=5000000
     REAL     ::DT,LTIME
     !  INTEGER,PARAMETER  ::NTIME=INT(LTIME/DT)+1
     INTEGER    ::NTIME
     REAL,PARAMETER     ::GOLDR = 5.0    ! the gold ratio for the linearsearch
     INTEGER,DIMENSION(:),ALLOCATABLE ::NOBCU_F,NOBCV_F,NOBCW_F
     INTEGER,DIMENSION(:),ALLOCATABLE ::NOBCU_F0,NOBCV_F0,NOBCW_F0
     REAL,DIMENSION(:), ALLOCATABLE::GOLD,D,CVAOLD
     INTEGER ::NOSTIM,NOSNOD
     REAL,DIMENSION(:,:),ALLOCATABLE ::UEXAC,VEXAC,WEXAC
     REAL,DIMENSION(:),ALLOCATABLE::SX,SY,SZ,STIME
     REAL ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
     ! The variables for the mesh adaptivity....
     LOGICAL  ::ADMESHH,D33
     CHARACTER*480 ::FILNAM_INPUT,FILNAM_SOURCE,FILNAM_F,FILNAM_AD
     
     
     ! Free surface
     INTEGER,DIMENSION(:),ALLOCATABLE ::NOBCH_F
     INTEGER,DIMENSION(:),ALLOCATABLE ::NOBCH_F0
     LOGICAL ::GOTFRE
     
     ! Model covariance
     
     REAL,DIMENSION(:),ALLOCATABLE::gcovari
     REAL,DIMENSION(:),ALLOCATABLE::ogcovari
     REAL      ::F_SMOOTHNESS
     SAVE      ::F_SMOOTHNESS
     
     ! Local ......
     INTEGER ::II,jj,kk,k0,k1                        ! for the loop
     INTEGER ::IO,N
     REAL    ::A,MAXGRAD_AVE,F2,F0,COF
     INTEGER ::STATUS
     REAL,DIMENSION(:),ALLOCATABLE::hcheck
     REAL FAI1,FAI,xc1
     REAL,PARAMETER ::PIE=3.141592654
     REAL FSZERO
     
     
     !******************************************************************************************************
     !   START ADJOINT MODEL
     
     ! Establish signal handlers.
     !  call initialise_signals
     
     ! Set up the "big" memory
     CALL memory_size(NIMEM, NRMEM)
     ALLOCATE(RMEM(NRMEM))
     ALLOCATE(IMEM(NIMEM))
     
     IF( .not.ALLOCATED(RMEM) ) STOP 101
     IF( .not.ALLOCATED(IMEM) ) STOP 102
     
     
     !         ewrite(3,*) 'running code'
     ewrite(3,*)  'ENTER THE FILE NAME TO READ FROM'
     ewrite(3,*) 'NB IF THIS IS A PARALLEL SIMULATION'
     ewrite(3,*) 'THEN WE MUST ENTER THE NAME OF THE '
     ewrite(3,*) 'FILE WITHOUT ANY SUFIX'
     ewrite(3,*) '240 characters max & if it begins with /'
     ewrite(3,*) 'then assume chomplete path has been entered'
     ewrite(3,*) 'ENTER INPUT FILE NAME'
     ewrite(3,*) ' ****************************** '
     ewrite(3,*)  'FILNAM_SOURCE'
     ewrite(3,*) 'FILNAM_F'
     ewrite(3,*) 'FILNAM_AD'
     ewrite(3,*) 'NGLOITS'
     ewrite(3,*)  'DT,LTIME,NTIME=INT(LTIME/DT)+1'
     ewrite(3,*)  'NOSTIM,NOSNOD'
     ewrite(3,*)  'GOTFRE,ADMESHH,D33'
     ewrite(3,*) ' ****************************** '
     
     READ '(A480)', FILNAM_INPUT
     !  ewrite(3,*) 'ENTER SOURCE FILE NAME',
     !          READ '(A240)', FILNAM_SOURCE
     !  ewrite(3,*) 'ENTER FORWARD FILE NAME',
     !          READ '(A240)', FILNAM_F
     !  ewrite(3,*) 'ENTER ADJOINT FILE NAME',
     !          READ '(A240)', FILNAM_AD
     
     OPEN(3,FILE=FILNAM_INPUT,STATUS='OLD')
     READ(3,*) FILNAM_SOURCE
     ewrite(3,*)  'FILNAM_SOURCE',FILNAM_SOURCE
     READ(3,*) FILNAM_F
     ewrite(3,*) 'FILNAM_F',FILNAM_F
     READ(3,*) FILNAM_AD
     ewrite(3,*) 'FILNAM_AD',FILNAM_AD
     READ(3,*) NGLOITS
     ewrite(3,*) 'NGLOITS',NGLOITS
     READ(3,*) DT,LTIME,NTIME
     ewrite(3,*)  'DT,LTIME,NTIME=INT(LTIME/DT)+1',DT,LTIME,NTIME
     READ(3,*) NOSTIM,NOSNOD
     ewrite(3,*)  'NOSTIM,NOSNOD',NOSTIM,NOSNOD
     READ(3,*) GOTFRE,ADMESHH,D33
     ewrite(3,*)  'GOTFRE,ADMESHH,D33',GOTFRE,ADMESHH,D33
     READ(3,*) F0,COF
     ewrite(3,*)  'F0,COF',F0,COF
     READ(3,*) ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
     ewrite(3,*)  'ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z',ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
     READ(3,*) FSZERO
     ewrite(3,*)  FSZERO
     
     CLOSE(3)
     
     
     ALLOCATE(NOBCU_F(NTIME))
     ALLOCATE(NOBCV_F(NTIME))
     IF(D33) ALLOCATE(NOBCW_F(NTIME))
     ALLOCATE(NOBCU_F0(NTIME))
     ALLOCATE(NOBCV_F0(NTIME))
     IF(D33) ALLOCATE(NOBCW_F0(NTIME))
     ALLOCATE(UEXAC(NOSNOD,NOSTIM))
     ALLOCATE(VEXAC(NOSNOD,NOSTIM))
     IF(D33) ALLOCATE(WEXAC(NOSNOD,NOSTIM))
     ALLOCATE(SX(NOSNOD))
     ALLOCATE(SY(NOSNOD))
     IF(D33) ALLOCATE(SZ(NOSNOD))
     ALLOCATE(STIME(NOSTIM))
     
     ! ******test free surface
     ALLOCATE(NOBCH_F(NTIME)) 
     ALLOCATE(NOBCH_F0(NTIME)) 
     NOBCH_F(1:NTIME) = 0
     NOBCH_F0(1:NTIME) = 0
     !****
     
     
     !NULLIFY(CVAOLD)
     !NULLIFY(GOLD)
     !NULLIFY(D)
     
     ! *********************************************
     ! Start the simulation
     
     FUNCT =0.0
     MAXGRAD = 0.0
     !  NOCVA1 =3000000
     !  NOCVA0 =3000000
     NOCVA1 =5000000
     NOCVA0 =5000000
     NOCVA = NOCVA1
     ALLOCATE(CVA(NOCVA))
     ALLOCATE(IMEM_CVA(NOCVA))  
     ALLOCATE(G(NOCVA))  
     ALLOCATE(CVA_X(NOCVA))
     ALLOCATE(CVA_Y(NOCVA))
     ALLOCATE(CVA_Z(NOCVA))
     
     G(1:NOCVA) = 0.0
     CVA(1:NOCVA) = 0.0
     CVA_X(1:NOCVA) = 0.0
     CVA_Y(1:NOCVA) = 0.0
     CVA_Z(1:NOCVA) = 0.0
     ! Model covariance 
     ALLOCATE(gcovari(NOCVA))
     gcovari(1:NOCVA) = 0.0
     
     DO II=1,NOCVA
        IMEM_CVA(II)=0
     ENDDO
     
     DO II =1,NTIME
        NOBCU_F(II)=0
        NOBCV_F(II)=0
        IF(D33) THEN
           NOBCW_F(II)=0
        ENDIF
     ENDDO
     DO II =1,NTIME
        NOBCU_F0(II)=0
        NOBCV_F0(II)=0
        IF(D33) THEN
           NOBCW_F0(II)=0
        ENDIF
     ENDDO
     
     
     
     
     !  CALL EXASOL(UEXAC,VEXAC,WEXAC,SX,SY,SZ,D3,                            &
     !                     STIME,LTIME,DT,NOSTIM,NOSNOD,FILNAM_SOURCE)
     
     ! other........
     
     ALLOCATE(CVAOLD(NOCVA),STAT=status)
     IF( ALLOCATED(CVAOLD) ) THEN
        CVAOLD(1:NOCVA) = 0.0
     ELSE
        ewrite(-1,*)  "ERROR: Failed to allocate CVAOLD"
        ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
        STOP
     END IF
     ALLOCATE(GOLD(NOCVA),STAT=status)  
     IF( ALLOCATED(GOLD) ) THEN
        GOLD(1:NOCVA) = 0.0
     ELSE
        ewrite(-1,*)  "ERROR: Failed to allocate GOLD"
        ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
        STOP
     END IF
     ALLOCATE(D(NOCVA),STAT=status)
     IF( ALLOCATED(D) ) THEN
        D(1:NOCVA) = 0.0
     ELSE
        ewrite(-1,*)  "ERROR: Failed to allocate D"
        ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
        STOP
     END IF
     
     LINITS = 1
     DO GLOITS=1,NGLOITS
        ewrite(3,*) '*********************'
        ewrite(3,*) 'gloits=',gloits  
        ! Solve the Navier Stoke eqn. and calculate the cost function..
        !--------------------------------------------------------------
        !This is a stub for F_FLUID without event option.
        
        
        FILEVE=filename
        CALL cpu_time( TIME1 )
        
        ewrite(3,*)  'before calling'
        ewrite(3,*) NTIME,FUNCT,NOCVA0,NOCVA,MAXREC
        ewrite(3,*) 'GLOITS,LINITS=',GLOITS,LINITS
        CALL FLUIDS_FORWARD(IEVENT,EVESFL,EVENTD,EVENTT,EVENDE,EVENTV,       &
             EVECO1,EVECO2,EVECO3,EVECO4,                                    &
             TEMPR1,TEMPR2,TEMPR3,TEMPR4,BURNUP1,                           &
             EVEFIR,EVEDT,EVEUDT,EVEDTN, NENERG, NDELA2,                     &
             NNNODS,NNNELE,                                                  &
             FILEVE,EVGDIN,                                                  &
             EDENIN,ETEMIN, EDENGA,EGAMD2,EGAMD3,                            &
             VOILIM, EVTIR1,EVTIR2, EVTIOP,EVSBCY,EVMXOP,                    &
             NMATEV,MXNMTT,TEMEVE,                                           &
             NDIMI,MXNDEV,NDSIEV, VERSEV,GOTGAS,                             &
             RMEM, NRMEM, IMEM, NIMEM,                                       &
             ! New staff for the data assmilation...
             FUNCT,NOCVA0,NOCVA,MAXREC,CVA,IMEM_CVA,G,MAXGRAD,               &
             GLOITS,LINITS,                                                  &
             CVA_X,CVA_Y,CVA_Z,                                              &
             NOBCU_F,NOBCV_F,NOBCW_F,NTIME,D33,                               &
             NOBCU_F0,NOBCV_F0,NOBCW_F0,                                        &
             ADMESHH,NOSTIM,NOSNOD,STIME,                                     &
             UEXAC,VEXAC,WEXAC,                                              &
             SX,SY,SZ,                                                       &
             ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z,                                &
             FILNAM_INPUT,FILNAM_SOURCE,FILNAM_F,FILNAM_AD,                  &
             CVAOLD,GOLD,D,                                                  &
             NOBCH_F,NOBCH_F0,gcovari,F_SMOOTHNESS,xc1)
        ewrite(3,*) 'funct=',funct
        
        ewrite(3,*) 'out of forward fluids'
        CALL cpu_time( TIME2 )
        
        ewrite(3,*)  'start,end,total cpu time',TIME1,TIME2,TIME2-TIME1
        
        ewrite(3,*) 'NOCVA in MAIN',NOCVA
        ! Re-allocate exactly the Size of CVA and G
        NOCVA = 0
        IF(GOTFRE) THEN
           DO II = 1,NTIME
              NOCVA = NOCVA+NOBCH_F(II)
           ENDDO
        ELSE
           DO II = 1,NTIME
              NOCVA = NOCVA+NOBCU_F(II)+NOBCV_F(II)
              IF(D33) THEN
                 NOCVA = NOCVA+NOBCW_F(II)
              ENDIF
           ENDDO
        ENDIF
        DO II = 1,NTIME
           IF(GOTFRE) THEN
              !free surface
              NOBCH_F0(II) = NOBCH_F(II) 
           ELSE
              NOBCU_F0(II) = NOBCU_F(II) 
              NOBCV_F0(II) = NOBCV_F(II) 
              IF(D33) NOBCW_F0(II) = NOBCW_F(II) 
           ENDIF
        ENDDO
        
        ! Reallocate CVA
        ALLOCATE(OCVA(NOCVA),STAT=status)
        OCVA(1:NOCVA) = CVA(1:NOCVA)
        DEALLOCATE(CVA)
        ALLOCATE(CVA(NOCVA),STAT=status)
        CVA(1:NOCVA) = OCVA(1:NOCVA)
        DEALLOCATE(OCVA)
        ewrite(3,*) '1111'
        ! Reallocate IMEM_CVA
        ALLOCATE(IMEM_OCVA(NOCVA),STAT=status)
        IMEM_OCVA(1:NOCVA) = IMEM_CVA(1:NOCVA)
        DEALLOCATE(IMEM_CVA)
        ALLOCATE(IMEM_CVA(NOCVA),STAT=status)
        IMEM_CVA(1:NOCVA) = IMEM_OCVA(1:NOCVA)
        DEALLOCATE(IMEM_OCVA)
        ewrite(3,*) '2222'
        ! Reallocate G
        ALLOCATE(OG(NOCVA),STAT=status)
        OG(1:NOCVA) = G(1:NOCVA)
        DEALLOCATE(G)
        ALLOCATE(G(NOCVA),STAT=status)
        G(1:NOCVA) = OG(1:NOCVA)
        DEALLOCATE(OG)
        ewrite(3,*) '3333'
        
        ! Reallocate CVA_X
        ALLOCATE(OCVA_X(NOCVA),STAT=status)
        OCVA_X(1:NOCVA) = CVA_X(1:NOCVA)
        DEALLOCATE(CVA_X)
        ALLOCATE(CVA_X(NOCVA),STAT=status)
        CVA_X(1:NOCVA) = OCVA_X(1:NOCVA)
        DEALLOCATE(OCVA_X)
        ewrite(3,*) '4444'
        ! Reallocate CVA_Y
        ALLOCATE(OCVA_Y(NOCVA),STAT=status)
        OCVA_Y(1:NOCVA) = CVA_Y(1:NOCVA)
        DEALLOCATE(CVA_Y)
        ALLOCATE(CVA_Y(NOCVA),STAT=status)
        CVA_Y(1:NOCVA) = OCVA_Y(1:NOCVA)
        DEALLOCATE(OCVA_Y)
        ewrite(3,*) '5555'
        ! Reallocate CVA_Z
        ALLOCATE(OCVA_Z(NOCVA),STAT=status)
        OCVA_Z(1:NOCVA) = CVA_Z(1:NOCVA)
        DEALLOCATE(CVA_Z)
        ALLOCATE(CVA_Z(NOCVA),STAT=status)
        CVA_Z(1:NOCVA) = OCVA_Z(1:NOCVA)
        DEALLOCATE(OCVA_Z)
        ewrite(3,*) '6666'
        ! Reallocate gcovari
        ALLOCATE(ogcovari(NOCVA),STAT=status)
        ogcovari(1:NOCVA) = gcovari(1:NOCVA)
        DEALLOCATE(gcovari)
        ALLOCATE(gcovari(NOCVA),STAT=status)
        gcovari(1:NOCVA) = ogcovari(1:NOCVA)
        DEALLOCATE(ogcovari)
        
        !Reallocate CVAOLD, GOLD, D
        
        ALLOCATE(OCVAOLD(NOCVA))
        OCVAOLD(1:NOCVA)=CVAOLD(1:NOCVA)
        DEALLOCATE(CVAOLD)
        ALLOCATE(CVAOLD(NOCVA),STAT=status)
        IF( ALLOCATED(CVAOLD) ) THEN
           CVAOLD(1:NOCVA) = 0.0
        ELSE
           ewrite(-1,*)  "ERROR: Failed to allocate CVAOLD"
           ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
           STOP
        END IF
        CVAOLD(1:NOCVA)=OCVAOLD(1:NOCVA)
        DEALLOCATE(OCVAOLD)
        
        
        ALLOCATE(OGOLD(NOCVA))
        OGOLD(1:NOCVA)=GOLD(1:NOCVA)
        DEALLOCATE(GOLD)
        ALLOCATE(GOLD(NOCVA),STAT=status)  
        IF( ALLOCATED(GOLD) ) THEN
           GOLD(1:NOCVA) = 0.0
        ELSE
           ewrite(-1,*)  "ERROR: Failed to allocate GOLD"
           ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
           STOP
        END IF
        GOLD(1:NOCVA)=OGOLD(1:NOCVA)
        DEALLOCATE(OGOLD)
        
        
        ALLOCATE(OD(NOCVA))
        OD(1:NOCVA)=D(1:NOCVA)
        DEALLOCATE(D)
        ALLOCATE(D(NOCVA),STAT=status)
        IF( ALLOCATED(D) ) THEN
           D(1:NOCVA) = 0.0
        ELSE
           ewrite(-1,*)  "ERROR: Failed to allocate D"
           ewrite(-1,*)  "ERROR: ALLOCATE STAT = ", status
           STOP
        END IF
        D(1:NOCVA)=OD(1:NOCVA)
        DEALLOCATE(OD)
        
        G(1:NOCVA) = 0.0
        
        ! Solve adjoint eqn, 
        ! and calculate the gradient & save the gradient in the file: gradient_bc.dat ...
        !--------------------------------------------------------------------------
        CALL cpu_time( TIME1 )
        ewrite(3,*) 'Running the adjoint model'
        !    ADMESH=.FALSE. 
        ewrite(3,*) 'admesh before adjoint in main',admeshh
        !     open(1,file='cva2.dat')
        !       write(1,*) 'nocva=',nocva
        !       write(1,*) (cva(ii),ii=1,nocva)
        !     close(1)
        !     open(1,file='imem_cva2.dat')
        !       write(1,*) 'nocva=',nocva
        !       write(1,*) (imem_cva(ii),ii=1,nocva)
        !     close(1)
        
        ewrite(3,*) 'funct=', funct
        !        stop 10
        
        CALL FLUIDS_ADJOINT(IEVENT,EVESFL,EVENTD,EVENTT,EVENDE,EVENTV,       &
             EVECO1,EVECO2,EVECO3,EVECO4,                                    &
             TEMPR1,TEMPR2,TEMPR3,TEMPR4,BURNUP1,                            &
             EVEFIR,EVEDT,EVEUDT,EVEDTN, NENERG, NDELA2,                     &
             NNNODS,NNNELE,                                                  &
             FILEVE,EVGDIN,                                                  &
             EDENIN,ETEMIN, EDENGA,EGAMD2,EGAMD3,                            &
             VOILIM, EVTIR1,EVTIR2, EVTIOP,EVSBCY,EVMXOP,                    &
             NMATEV,MXNMTT,TEMEVE,                                           &
             NDIMI,MXNDEV,NDSIEV, VERSEV,GOTGAS,                             &
             RMEM, NRMEM, IMEM, NIMEM,                                       &
             ! New staff for the data assmilation...
             FUNCT,NOCVA0,NOCVA,MAXREC,CVA,IMEM_CVA,G,MAXGRAD,              &
             GLOITS,LINITS,                                                 &
             CVA_X,CVA_Y,CVA_Z,                                             &
             NOBCU_F,NOBCV_F,NOBCW_F,NTIME,D33,                              &
             NOBCU_F0,NOBCV_F0,NOBCW_F0,                                        &
             ADMESHH,NOSTIM,NOSNOD,STIME,                                    &
             UEXAC,VEXAC,WEXAC,                                             &
             SX,SY,SZ,                                                      &
             ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z,                               &
             FILNAM_INPUT,FILNAM_SOURCE,FILNAM_F,FILNAM_AD,                 &
             CVAOLD,GOLD,D,                                                 &
             NOBCH_F,NOBCH_F0,gcovari,F_SMOOTHNESS)
        
        ewrite(3,*) 'out of adjoint fluids'
        CALL cpu_time( TIME2 )
        ewrite(3,*)  'start,end,total cpu time',TIME1,TIME2,TIME2-TIME1
     ENDDO
     
     
     DEALLOCATE(RMEM)
     DEALLOCATE(IMEM)
     DEALLOCATE(CVA)
     DEALLOCATE(CVA_X)
     DEALLOCATE(CVA_Y)
     IF(D33) DEALLOCATE(CVA_Z)
     DEALLOCATE(gcovari)
     DEALLOCATE(IMEM_CVA)
     DEALLOCATE(G)  
     DEALLOCATE(CVAOLD)
     DEALLOCATE(GOLD)  
     DEALLOCATE(D)  
     DEALLOCATE(NOBCU_F)
     DEALLOCATE(NOBCV_F)
     IF(D33) DEALLOCATE(NOBCW_F)
     DEALLOCATE(NOBCU_F0)
     DEALLOCATE(NOBCV_F0)
     IF(D33) DEALLOCATE(NOBCW_F0)
     DEALLOCATE(UEXAC)
     DEALLOCATE(VEXAC)
     IF(D33) DEALLOCATE(WEXAC)
     DEALLOCATE(SX)
     DEALLOCATE(SY)
     IF(D33) DEALLOCATE(SZ)
     DEALLOCATE(STIME)
     
     DEALLOCATE(NOBCH_F)
     DEALLOCATE(NOBCH_F0)
     
     ! End the whole test simulation
     !************************************************
     
     
     RETURN
   END SUBROUTINE ADFLUIDS_BC_TEST
   

   
   
   
   
   SUBROUTINE ADFLUIDS_SENSITIVITY(filename, filenamelen)
     
     ! ************************************************************************
     ! THE MAIN PROGRAM IS USED TO SIMUALATE 3D FLOW. WE USE THE ADJOINT METHOD
     ! TO ASSIMILATE THE OBSERVED DATA IN ORDER TO IMPROVE THE PREDICTION OF 
     ! THE SIMULATION.THE BASIC MODELS ARE THE NAVIER STOKES, RADIATION(WITH OR 
     ! WITHOUT EVENT) AND OR ADVECTION-DIFFUSION TYPE OF EQUATIONS.
     !NTIME : the number of the time steps
     !NOSNOD: the number of positions where we have obervational data
     !NOSTIM: the number of the time steps when we have the observation data.
     !STIME:the vector storing the time when we will invert the observational data
     !UEXAC,VEXAC,WEXAC: the velocity from the obervational data in x, y, z directions
     !SX,SY,SZ      : the coordination values x, y, z of the observational data
     !ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z: the spans for Gaussian data srpreading
     !FILNAM_INPUT: the input file name
     !FILNAM_SOURCE------the name of observational data file 
     !FILNAM_F------ the forward data file name which is obained by running GEM 
     !FILNAM_AD------the adjoint data file name which is obained by running GEM 
     !D3:  D3=.true. 3D, .false. 2D
     !ADMESH: admesh=.true. mesh adaptivity, .false. nonmesh-adaptivity
     ! ************************************************************************
     !
     use FLDebug
     use signals
     use fluids_adjoint_module
     use fluids_forward_module
     IMPLICIT NONE
     integer, intent(in) :: filenamelen
     character(len=filenamelen), intent(in) :: filename
     REAL EVESFL(1),EVENTD(1),EVENTT(1),EVENDE(1),EVENTV(1)
     REAL EVECO1(1),EVECO2(1),EVECO3(1),EVECO4(1)
     REAL TEMPR1(1),TEMPR2(1),TEMPR3(1),TEMPR4(1),BURNUP1(1)
     REAL TEMEVE(1)
     REAL::EVEDT=0
     REAL EVEUDT,EVEDTN,ETEMIN,TIME1,TIME2
     REAL EVGDIN,EDENIN,EDENGA,EGAMD2,EGAMD3
     REAL VOILIM,EVTIR1,EVTIR2
     INTEGER EVTIOP,EVSBCY,EVMXOP
     INTEGER VERSEV,NDIMI,NDSIEV(1)
     LOGICAL::EVEFIR=.FALSE.
     LOGICAL GOTGAS
     CHARACTER*240 FILEVE
     INTEGER::IEVENT=0, NENERG=0, NNNODS=0, NNNELE=0, NDELA2=0, NMATEV=0, MXNMTT=1
     INTEGER::MXNDEV=1
     
     INTEGER NRMEM,NIMEM
     REAL, ALLOCATABLE::RMEM(:)
     INTEGER, ALLOCATABLE::IMEM(:)
     
     !
     ! The extral variables for adjoint model, nonlinear conjugate gradient...
     REAL ::FUNCT,MAXGRAD   ! The cost function and maximum gradient
     INTEGER,parameter       ::NOCVA=100000        ! Number of the control variables
     REAL,DIMENSION(:),ALLOCATABLE::G,CVA ! Nonlinear CG
     INTEGER,DIMENSION(:),ALLOCATABLE::IMEM_CVA,IMEM_OCVA
     REAL,DIMENSION(:),ALLOCATABLE::CVA_X,CVA_Y,CVA_Z  
     INTEGER,DIMENSION(:),ALLOCATABLE ::NOBCU_F,NOBCV_F,NOBCW_F
     INTEGER,DIMENSION(:),ALLOCATABLE ::NOBCU_F0,NOBCV_F0,NOBCW_F0
     REAL          ::DT,LTIME
     INTEGER    ::NTIME,NOCVA0
     INTEGER ::NOSTIM,NOSNOD
     REAL,DIMENSION(:,:),ALLOCATABLE ::UEXAC,VEXAC,WEXAC
     REAL,DIMENSION(:),ALLOCATABLE::SX,SY,SZ,STIME
     REAL ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
     
     ! The variables for the mesh adaptivity....
     LOGICAL       ::ADMESHH,D33
     CHARACTER*480 ::FILNAM_INPUT,FILNAM_SOURCE,FILNAM_F,FILNAM_AD
     REAL,DIMENSION(:),ALLOCATABLE::gcovari
     REAL      ::F_SMOOTHNESS
     SAVE      ::F_SMOOTHNESS
     
     
     ! Free surface
     INTEGER,DIMENSION(:),ALLOCATABLE ::NOBCH_F
     INTEGER,DIMENSION(:),ALLOCATABLE ::NOBCH_F0
     LOGICAL ::GOTFRE
     
     ! Local ......
     INTEGER ::II,jj,kk,k0,k1                        ! for the loop
     INTEGER ::IO,N
     REAL    ::A,MAXGRAD_AVE,F2,F0,COF
     REAL,PARAMETER ::PIE=3.141592654
     REAL FSZERO
     
     
     ! Establish signal handlers.
     call initialise_signals
     
     ! Set up the "big" memory
     CALL memory_size(NIMEM, NRMEM)
     ALLOCATE(RMEM(NRMEM))
     ALLOCATE(IMEM(NIMEM))
     
     IF( .not.ALLOCATED(RMEM) ) STOP 101
     IF( .not.ALLOCATED(IMEM) ) STOP 102
     NOCVA0=10000
     
     !         ewrite(3,*) 'running code'
     ewrite(3,*)  'ENTER THE FILE NAME TO READ FROM'
     ewrite(3,*) 'NB IF THIS IS A PARALLEL SIMULATION'
     ewrite(3,*) 'THEN WE MUST ENTER THE NAME OF THE '
     ewrite(3,*) 'FILE WITHOUT ANY SUFIX'
     ewrite(3,*) '240 characters max & if it begins with /'
     ewrite(3,*) 'then assume chomplete path has been entered'
     ewrite(3,*) 'ENTER INPUT FILE NAME'
     ewrite(3,*) ' ****************************** '
     ewrite(3,*) 'FILNAM_SOURCE'
     ewrite(3,*) 'FILNAM_F'
     ewrite(3,*) 'FILNAM_AD'
     ewrite(3,*) 'DT,LTIME,NTIME=INT(LTIME/DT)+1'
     ewrite(3,*) 'NOSTIM,NOSNOD'
     ewrite(3,*) 'GOTFRE,ADMESHH,D33'
     ewrite(3,*) ' ****************************** '
     
     READ '(A480)', FILNAM_INPUT
     
     OPEN(3,FILE=FILNAM_INPUT,STATUS='OLD')
     READ(3,*) FILNAM_SOURCE
     ewrite(3,*)  'FILNAM_SOURCE',FILNAM_SOURCE
     READ(3,*) FILNAM_F
     ewrite(3,*) 'FILNAM_F',FILNAM_F
     READ(3,*) FILNAM_AD
     ewrite(3,*) 'FILNAM_AD',FILNAM_AD
     READ(3,*) DT,LTIME,NTIME
     ewrite(3,*)  'DT,LTIME,NTIME=INT(LTIME/DT)+1',DT,LTIME,NTIME
     READ(3,*) NOSTIM,NOSNOD
     ewrite(3,*)  'NOSTIM,NOSNOD',NOSTIM,NOSNOD
     READ(3,*) GOTFRE,ADMESHH,D33
     
     CLOSE(3)
     
     
     
     ! *********************************************
     ! Start the simulation
     
     FUNCT =0.0
     MAXGRAD = 0.0
     
     ! Model covariance 
     ALLOCATE(gcovari(NOCVA))
     ALLOCATE(CVA(NOCVA))
     ALLOCATE(IMEM_CVA(NOCVA))       
     ALLOCATE(G(NOCVA))       
     ALLOCATE(CVA_X(NOCVA))
     ALLOCATE(CVA_Y(NOCVA))
     ALLOCATE(CVA_Z(NOCVA))
     
     ! Model covariance 
     ALLOCATE(NOBCU_F(NTIME))
     ALLOCATE(NOBCV_F(NTIME))
     ALLOCATE(NOBCW_F(NTIME))
     ALLOCATE(NOBCU_F0(NTIME))
     ALLOCATE(NOBCV_F0(NTIME))
     ALLOCATE(NOBCW_F0(NTIME))
     ALLOCATE(UEXAC(NOSNOD,NOSTIM))
     ALLOCATE(VEXAC(NOSNOD,NOSTIM))
     ALLOCATE(WEXAC(NOSNOD,NOSTIM))
     ALLOCATE(SX(NOSNOD))
     ALLOCATE(SY(NOSNOD))
     ALLOCATE(SZ(NOSNOD))
     ALLOCATE(STIME(NOSTIM))
     ALLOCATE(NOBCH_F(NTIME)) 
     ALLOCATE(NOBCH_F0(NTIME)) 
     
     
     gcovari=0.0
     CVA=0.0
     IMEM_CVA=0       
     G=0.0       
     CVA_X=0.0
     CVA_Y=0.0
     CVA_Z=0.0
     
     NOBCU_F=0
     NOBCV_F=0
     NOBCW_F=0
     NOBCU_F0=0
     NOBCV_F0=0
     NOBCW_F0=0
     UEXAC=0.0
     VEXAC=0.0
     WEXAC=0.0
     SX=0.0
     SY=0.0
     SZ=0.0
     STIME=0.0
     NOBCH_F=0
     NOBCH_F0=0 
     
     !goto 333
     !*************************************
     ! RUN THE FORWARD MODEL
     !*************************************
     
     ! Solve the Navier Stoke eqn. and calculate the cost function..
     !--------------------------------------------------------------
     !This is a stub for F_FLUID without event option.
     
     FILEVE=filename
     CALL cpu_time( TIME1 )
     ewrite(3,*)  'before running the forward model'
     CALL FLUIDS_FORWARD(IEVENT,EVESFL,EVENTD,EVENTT,EVENDE,EVENTV,     &
          EVECO1,EVECO2,EVECO3,EVECO4,                                    &
          TEMPR1,TEMPR2,TEMPR3,TEMPR4,BURNUP1,                            &
          EVEFIR,EVEDT,EVEUDT,EVEDTN, NENERG, NDELA2,                     &
          NNNODS,NNNELE,                                                  &
          FILEVE,EVGDIN,                                                  &
          EDENIN,ETEMIN, EDENGA,EGAMD2,EGAMD3,                            &
          VOILIM, EVTIR1,EVTIR2, EVTIOP,EVSBCY,EVMXOP,                    &
          NMATEV,MXNMTT,TEMEVE,                                           &
          NDIMI,MXNDEV,NDSIEV, VERSEV,GOTGAS,                             &
          RMEM, NRMEM, IMEM, NIMEM,                                       &
          ! New staff for the data assmilation...
          0.0,NOCVA0,NOCVA,10000000,CVA,IMEM_CVA,G,1.0,                   &
          1,1,                                                            &
          CVA_X,CVA_Y,CVA_Z,                                              &
          NOBCU_F,NOBCV_F,NOBCW_F,NTIME,D33,                              &
          NOBCU_F0,NOBCV_F0,NOBCW_F0,                                     &
          ADMESHH,NOSTIM,NOSNOD,STIME,                                    &
          UEXAC,VEXAC,WEXAC,                                              &
          SX,SY,SZ,                                                       &
          1.0,1.0,1.0,1.0,                                                &
          FILNAM_INPUT,FILNAM_SOURCE,FILNAM_F,FILNAM_AD,                  &
          CVA,G,G,                                                        &
          NOBCH_F,NOBCH_F0,gcovari,0.0,0.01)

     ewrite(3,*) 'out of forward fluids'
     CALL cpu_time( TIME2 )
     
     !        CALL RCLEAR(G,NOCVA)
     G(1:NOCVA)=0.0
     
     
     !*************************************
     !  RUN THE ADJOINT MODEL
     !*************************************
     
     ! Solve adjoint eqn, and calculate the gradient
     ! save the gradient in the file: gradient_bc.dat ...
     !-----------------------------------------------------
     ewrite(3,*) 'Running the adjoint model'
     !     DEALLOCATE(RMEM)
     !     DEALLOCATE(IMEM)
     
     !  ALLOCATE(RMEM(NRMEM))
     !  ALLOCATE(IMEM(NIMEM))
     RMEM=0.0
     IMEM=0
333  continue
     !  call initialise_signals
     
     ! Set up the "big" memory
     !  CALL memory_size(NIMEM, NRMEM)
     !  ALLOCATE(RMEM(NRMEM))
     !  ALLOCATE(IMEM(NIMEM))
     FILEVE=filename
     CALL cpu_time( TIME1 )
     CALL FLUIDS_ADJOINT(IEVENT,EVESFL,EVENTD,EVENTT,EVENDE,EVENTV,     &
          EVECO1,EVECO2,EVECO3,EVECO4,                                    &
          TEMPR1,TEMPR2,TEMPR3,TEMPR4,BURNUP1,                            &
          EVEFIR,EVEDT,EVEUDT,EVEDTN, NENERG, NDELA2,                     &
          NNNODS,NNNELE,                                                  &
          FILEVE,EVGDIN,                                                  &
          EDENIN,ETEMIN, EDENGA,EGAMD2,EGAMD3,                            &
          VOILIM, EVTIR1,EVTIR2, EVTIOP,EVSBCY,EVMXOP,                    &
          NMATEV,MXNMTT,TEMEVE,                                           &
          NDIMI,MXNDEV,NDSIEV, VERSEV,GOTGAS,                             &
          RMEM, NRMEM, IMEM, NIMEM,                                       &
          ! New staff for the data assmilation...
          0.0,NOCVA0,NOCVA,10000000,CVA,IMEM_CVA,G,1.0,                   &
          1,1,                                                            &
          CVA_X,CVA_Y,CVA_Z,                                              &
          NOBCU_F,NOBCV_F,NOBCW_F,NTIME,D33,                              &
          NOBCU_F0,NOBCV_F0,NOBCW_F0,                                     &
          ADMESHH,NOSTIM,NOSNOD,STIME,                                    &
          UEXAC,VEXAC,WEXAC,                                              &
          SX,SY,SZ,                                                       &
          1.0,1.0,1.0,1.0,                                                &
          FILNAM_INPUT,FILNAM_SOURCE,FILNAM_F,FILNAM_AD,                  &
          CVA,G,G,                                                        &
          NOBCH_F,NOBCH_F0,gcovari,0.0)
     
     ewrite(3,*) 'out of adjoint fluids'
     CALL cpu_time( TIME2 )
     ewrite(3,*)  'start,end,total cpu time',TIME1,TIME2,TIME2-TIME1
     
     
     
     DEALLOCATE(RMEM)
     DEALLOCATE(IMEM)
     DEALLOCATE(CVA)
     DEALLOCATE(CVA_X)
     DEALLOCATE(CVA_Y)
     DEALLOCATE(CVA_Z)
     DEALLOCATE(gcovari)
     DEALLOCATE(IMEM_CVA)
     DEALLOCATE(G)       
     DEALLOCATE(NOBCU_F)
     DEALLOCATE(NOBCV_F)
     DEALLOCATE(NOBCW_F)
     DEALLOCATE(NOBCU_F0)
     DEALLOCATE(NOBCV_F0)
     DEALLOCATE(NOBCW_F0)
     DEALLOCATE(UEXAC)
     DEALLOCATE(VEXAC)
     DEALLOCATE(WEXAC)
     DEALLOCATE(SX)
     DEALLOCATE(SY)
     DEALLOCATE(SZ)
     DEALLOCATE(STIME)

     DEALLOCATE(NOBCH_F)
     DEALLOCATE(NOBCH_F0)



     ! End the whole simulation
     !************************************************
     
     
     RETURN
   END SUBROUTINE ADFLUIDS_SENSITIVITY

