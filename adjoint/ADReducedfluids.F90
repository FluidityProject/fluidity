  !    Copyright (C) 2006 Imperial College London and others.      
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
  
#include "fdebug.h"
    SUBROUTINE ADreducedFLUIDS_Initial(filename, filenamelen)
    use FLDebug
    use fluids_module
    use signals
    use AllSorts
    use NONLINCG_module
    IMPLICIT NONE
    
    !
    include '../assemble/paramnew.h'      
    
    
    !------------------
    ! Initial condition
    !------------------
    REAL varincr(istate) 
    REAL varin(istate)
    REAL phi0(0:ntimemax,mxpoi)
    REAL u0(0:ntimemax,mxpoi)
    REAL v0(0:ntimemax,mxpoi)
    REAL w0(0:ntimemax,mxpoi)
    common /fwdtraj/ phi0,u0,v0,w0
    
    !------------------------------------
    !  background variables
    !------------------------------
    integer ibackground
    REAL backterm(istate), bweights(istate)
    
    !---------------------------------
    !     POD variables
    !---------------------------------  
    integer nsvd,nsnap,nsvd_total,nsvd_u,nsvd_phi
    parameter (nsvd = 35)
    parameter (nsvd_u = 35)
    parameter (nsvd_phi = 35)
    parameter (nsvd_total = 3*nsvd_u+nsvd_phi)
    REAL leftsvd(istate,nsvd), svdval(nsvd_total)
    REAL coef_multivar(0:ntimemax,nsvd_total)
    REAL coef(0:ntimemax,nsvd)
    REAL romgrad(nsvd)
    REAL tcoef1(nsvd)
    REAL tcoef_multivar(nsvd_total)
    REAL smean(istate)
    integer nsvd1    
    REAL,DIMENSION(:,:),ALLOCATABLE   :: snapmatrix
    
    !----------------------------------------------------!
    !     Dual-Weighted Method                           !
    !----------------------------------------------------!
    REAL,DIMENSION(:),ALLOCATABLE::  dual_wei_space,FERROR
    LOGICAL IF_dualwei
    !----------------------------------------------------!
    !    Error measurement                         !
    !----------------------------------------------------!
    REAL GLOERR(mxpoi)
    LOGICAL IF_remesh,IF_REMESH_OBS,IF_REMESH_VORTICITY,IF_INVERSION_MESH
    INTEGER remesh
    REAL,DIMENSION(:,:),ALLOCATABLE   :: U0_ADJ,V0_ADJ,W0_ADJ,PHI0_ADJ
    REAL DF
    !----------------------------------------------------!
    !               FOR FLUIIDTY                         !  
    !----------------------------------------------------!
    integer, intent(in) :: filenamelen
    character(len=filenamelen), intent(in) :: filename
    REAL EVESFL(1),EVENTD(1),EVENTT(1),EVENDE(1),EVENTV(1)
    REAL EVECO1(1),EVECO2(1),EVECO3(1),EVECO4(1)
    REAL TEMPR1(1),TEMPR2(1),TEMPR3(1),TEMPR4(1),BURNUP1(1)
    REAL TEMEVE(1)
    REAL EVEDT,EVEUDT,EVEDTN,ETEMIN,TIME1,TIME2
    REAL EVGDIN,EDENIN,EDENGA,EGAMD2,EGAMD3
    REAL VOILIM,EVTIR1,EVTIR2
    INTEGER EVTIOP,EVSBCY,EVMXOP
    INTEGER VERSEV,NDIMI,MXNDEV,NDSIEV(1)
    LOGICAL EVEFIR,GOTGAS
    CHARACTER*240 FILEVE
    INTEGER IEVENT,NENERG,NNNODS,NNNELE,NDELA2,NMATEV,MXNMTT
    INTEGER NRMEM,NIMEM
    REAL, ALLOCATABLE::RMEM(:)
    INTEGER, ALLOCATABLE::IMEM(:)
    
    REAL,DIMENSION(:),ALLOCATABLE     ::    MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ
    REAL,DIMENSION(:),ALLOCATABLE     ::    WEIGHT,DETWEI
    INTEGER,DIMENSION(:),ALLOCATABLE  ::    NDGLNO,XONDGL,PNDGLN,SNDGLN
    REAL,DIMENSION(:),ALLOCATABLE     ::    X,Y,Z,NU,NV,NW,UG,VG,WG
    REAL,DIMENSION(:,:),ALLOCATABLE   ::    N,NLX,NLY,NLZ,M,MLX,MLY,MLZ
    REAL,DIMENSION(:,:),ALLOCATABLE   ::    NX,NY,NZ
    INTEGER,DIMENSION(:),ALLOCATABLE  ::    FINDRM,CENTRM
    !-----
    !  GEOSTROPHIC PRESSURE
    !----------------------
    INTEGER MPGLOC,IFPG,PGNODS
    REAL, ALLOCATABLE, DIMENSION(:,:)::MPG
    REAL, ALLOCATABLE, DIMENSION(:,:)::MPGLX
    REAL, ALLOCATABLE, DIMENSION(:,:)::MPGLY
    REAL, ALLOCATABLE, DIMENSION(:,:)::MPGLZ
    REAL, ALLOCATABLE, DIMENSION(:,:)::MPGX
    REAL, ALLOCATABLE, DIMENSION(:,:)::MPGY
    REAL, ALLOCATABLE, DIMENSION(:,:)::MPGZ
    REAL, ALLOCATABLE, DIMENSION(:,:)  ::PG,PG_POD
    REAL, ALLOCATABLE, DIMENSION(:,:)  ::snapmatrix_PG
    INTEGER, ALLOCATABLE, DIMENSION(:):: PGNDGLNO
    REAL, ALLOCATABLE, DIMENSION(:)  ::smean_PG
    REAL, ALLOCATABLE, DIMENSION(:,:)  ::leftsvd_PG,leftsvd_PG_x,leftsvd_PG_y
    REAL, ALLOCATABLE, DIMENSION(:)  ::varin_PG,varincr_PG,svdval_PG
    REAL, ALLOCATABLE, DIMENSION(:,:)  ::COEF_PG
    REAL, ALLOCATABLE, DIMENSION(:)  ::TCOEF_PG
    INTEGER NPGCOLM0,NPGCOLM
    INTEGER, ALLOCATABLE, DIMENSION(:)::PGFINDRM,PGCOLM,PGCENTRM
    !       INTEGER PGFINDRM(PGNODS+1),PGCOLM(NPGCOLM)
    !       INTEGER PGCENTRM(PGNODS) 
    
    
    CHARACTER(40) FILENAME_PG,FILENAME_MPG,FILENAME1,FILENAME2
    INTEGER K3
    
    ! for coriolis force
    INTEGER GEOBAL
    REAL SCFACTH0
    
    LOGICAL D3,DCYL
    INTEGER NONODS,XNONOD,FREDOP,NPRESS
    INTEGER TOTELE,NLOC,NGI,MLOC,TOTELE0
    INTEGER STOTEL0,STOTEL,SNLOC
    INTEGER NOD,ELE,II,JJ,KK,NSTEP,I,J,K,IP,NF
    REAL DT,LTIME
    
    REAL ftest,ftest2
    real volele,volvol,rnn,VOLUME
    real ml(mxpoi)
    integer iloc,jloc,igl,jgl,itime
    integer NUM
    
    !---------------------------!
    ! obervational data         !
    !---------------------------!
    integer iprojop(mxpoi), idata, ipertdata,iuseobs,iusemean
    INTEGER nsnap_obs,SnapNDT_obs,SnapNDT_obs1
    PARAMETER (SnapNDT_obs=SnapNdT)
    PARAMETER (nsnap_obs=INT(ntimemax/SnapNDT_obs)+1)
    REAL weights(istate)
    common /background/ backterm   
    REAL leftsvd_obs(istate,nsvd), svdval_obs(nsvd_total)
    REAL coef_multivar_obs(0:ntimemax,nsvd_total)
    
    REAL,DIMENSION(:,:),ALLOCATABLE   :: obsv,obsvpert
    INTEGER,DIMENSION(:,:),ALLOCATABLE   ::iflagobs
    REAL,DIMENSION(:,:),ALLOCATABLE   :: snapmatrix_obs
    REAL,DIMENSION(:),ALLOCATABLE     :: smean_obs
    
    ! Wind stress
    !------------
    LOGICAL IFWIND
    
    !--------------------------------------------!
    !     parameters for optimization subroutine !
    !--------------------------------------------!
    integer mxfun, itermax, nmeth, iout, mdim, idev
    integer IFUN,ITER,NFLAG
    REAL pgtol, acc
    REAL w(NSVD*(NSVD+7)/2), w1(5*NSVD+2)
    REAL cost,cost1
    
    CHARACTER*240 tu
    CHARACTER*240 cho
    CHARACTER*240 NAME
    
    ! REAL  time1,time2,cputime,costinit
    ! common /timing/ time1,time2,cputime,costinit   
    !*********************************************
    ! external calcfg
    
    ! For adjoint Model
    REAL, DIMENSION(0:ntimemax,NSVD_TOTAL) :: coef_multivar_ADJ
    LOGICAL ::LINCONVER,IRUNFULLFWD
    INTEGER ::LINITS,GLOITS,GLOITS_POD   ! nonlinear conjugate and linearsearch
    INTEGER,PARAMETER ::NGLOITS=1    ! nonlinear conjugate and linearsearch
    INTEGER ::NOCVA      ! number of the control variables
    REAL    ::VALFUN,FUNCT
    REAL,DIMENSION(:),ALLOCATABLE     :: CVA,CVAOLD
    REAL,DIMENSION(:),ALLOCATABLE     :: G,GOLD
    REAL,DIMENSION(:),ALLOCATABLE     :: D
    REAL,DIMENSION(:),ALLOCATABLE     :: U_ADJ,V_ADJ,W_ADJ
    REAL,DIMENSION(:),ALLOCATABLE     :: varincr_adj
    !
    REAL,PARAMETER ::GOLDR = 2.0    ! the gold ratio for the linearsearch
    REAL,PARAMETER::COF=1.0
    REAL    ::MAXGRAD,MAXGRAD_AVE
    REAL,PARAMETER ::CONJCONVER=1.0E-14    
    !  REAL,PARAMETER ::InitialPert=5.0E-1
    !**  REAL,PARAMETER ::PertCoef=0.1
    REAL,PARAMETER ::PertCoef=0.0
    !  REAL,PARAMETER ::PertCoef=0.0
    REAL,PARAMETER ::InitialPert=0.0
    
    ! error analyssi
    REAL::error_u(mxpoi),error_v(mxpoi),error_w(mxpoi),error_phi(mxpoi)
    REAL,DIMENSION(:,:),ALLOCATABLE  :: leftsvd_obs1,coef_multivar_obs1
    REAL,DIMENSION(:),ALLOCATABLE  :: smean_obs1
    
    
    ! -----rerunforward-------------------------------------------------
    ! 0-run the forward with the exact initial/boundary conditions
    ! 1-run the forward with exact + InitialPert/ exact+PertCoef*exact
    ! 2- inversion of the intial conditions with
    !     excat + InitialPert with/without running the original forward model
    ! 3- inversion of the intial coefficient (for POD bases) with
    !     excatCoef+PerturbationCoef
    INTEGER,PARAMETER ::rerunforward=2  !3
    !------------------------------------------------------------------
    
    !local
    REAL U_num,V_num,W_num,U_obs,V_obs,W_obs
    INTEGER ISNAP,NO_obs
    INTEGER IOS
    LOGICAL IFWRITE,IF_INTERPOL 
    LOGICAL IFSIMPLE_PODada
    LOGICAL IFBG_FUN
    LOGICAL IF_H1
    REAL    epsilon
    
#ifdef HAVE_LIBARPACK
    
    
    !****************************************************!
    !       INPUT and INITILIZE THE VARIABLES            !
    !****************************************************!
    
    write(3,*) 'mx,mxpoi,istate,ntimemax,nrsnapshots,SnapNDT_ob,',mx,mxpoi,istate,ntimemax,nrsnapshots,SnapNDT_obs
    
    IRUNFULLFWD = .false.
    IFWRITE=.false.
    IF_INTERPOL=.true.
    IFSIMPLE_PODada=.false.
    IFBG_FUN=.false.
    IFWIND=.true.
    IFPG=1
    IF_dualwei=.false.
    IF_remesh=.true.
    IF_REMESH_OBS=.false.
    IF_REMESH_VORTICITY=.false.
    IF_INVERSION_MESH=.true.
    remesh=0
    NO_obs=10
    IF_H1=.false.
    ! epsilon = 0.05 !for Re=400
    epsilon = 0.00008 !cp 0.000015
    DF =0.000003
    !* DF = 0.0004
    !** DF = 0.0001
    
    !--------4D-Var setup ---------------------------------
    idata = 0                 ! 1 if realistic data
    ipertdata = 1             ! 1 if perturbed observations
    ibackground = 1           ! 1 if background included in cost
    
    !-------POD setup-------------
    iuseobs = 0               ! 0 if use snapshots of the flow
    iusemean = 1              
    ! ge.1 if use snapmatrix centered around mean, 
    ! 1 use the mean of forward results, 2 use the mean of observtions
    !------------------------------------------------------------------------
    
    
    ! Initise the u0,v0,w0,phi0....
    !--------------------------------
    
    DO K=0,ntimemax
       DO I=1,mxpoi
          u0(k,i)=0.0
          v0(k,i)=0.0
          w0(k,i)=0.0
          phi0(k,i)=0.0
       ENDDO
    ENDDO
    
    ! ADJOINT = .true.
    LINCONVER = .true.
    write(3,*) 'allocte cva,g...'
    
    ! allocate vectors for controls and optimisation
    !-------------------------------------------------
    NOCVA = nsvd_total
    ALLOCATE(CVA(NOCVA))
    ALLOCATE(CVAOLD(NOCVA))
    ALLOCATE(G(NOCVA))
    ALLOCATE(GOLD(NOCVA))
    ALLOCATE(D(NOCVA))
    
    !initialise
    write(3,*) 'initiise cvs.....'
    CVA=0.0
    CVAOLD=0.0
    G=0.0
    GOLD=0.0
    D=0.0
    
    coef_multivar_ADJ =0.0
    COEF_multivar_OBS =0.0
    
    
    !**************************************************!
    !   SETUP OBSERVATIONAL DATA                       !
    !**************************************************!
    ALLOCATE(obsv(0:ntimemax,istate))
    ALLOCATE(iflagobs(0:ntimemax,mxpoi))
    ALLOCATE(smean_obs(istate))
    write(3,*) 'initialise obsv...'
    obsv=0.0
    iflagobs=0
    smean_obs=0.0
    do itime=0,ntimemax
       do ip=1,mxpoi
          if (SnapNDT_obs*int(itime/SnapNDT_obs) .eq. itime) then  
             iflagobs(itime,ip) = 1
          else
             iflagobs(itime,ip) = 0
          end if
       enddo
    enddo
    
8000 CONTINUE
    write(3,*) 'allocte the obsv.....'
    
    ! Dual-Weighted set-up
    !---------------------
    IF(IF_dualwei.and.remesh.ge.1) THEN
       !dual_wei=1.0/dble(nrsnapshots)
       allocate(dual_wei_space(mxpoi))
       allocate(FERROR(9*nonods))
       dual_wei_space =0.0
       FERROR=0.0
       open(1,file='FERROR.dat')
       read(1,*) (FERROR(NOD),NOD=1,NONODS*9)
       close(1)
       
       !  CALL DUALWEIGHT_SPACE(dual_wei_space,FERROR,3,.true.,X,Y,Z,   &
       !                 NONODS,TOTELE,NLOC,NCOLM,NDGLNO,FINDRM,COLM)
       deallocate(FERROR)
    ENDIF
    
    
    
    IF(rerunforward.eq.3.or.rerunforward.eq.2) then
       !---------------------------------------------------------!
       ! Read observational data from files                        !
       !---------------------------------------------------------*!
       IF(remesh.eq.0) THEN
          !goto 110
          ewrite(3,*) 'read observation data'
          open(10, file ='observation.dat')
          do ii=0,ntimemax
             read(10, 2000) (obsv(ii,jj), jj=1,istate)
          enddo
          close(10)
          ewrite(3,*) 'after observation data'
110       continue
          !   nonods=7254
          !        write(name,'(a,i0)') 'observation-newmesh.dat.1'
          !        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
          !        do ii=0,ntimemax
          !           read(1, 2000) (obsv(ii,jj), jj=1,nonods)
          !           read(1, 2000) (obsv(ii,jj), jj=1*mxpoi+1,1*mxpoi+nonods)
          !           read(1, 2000) (obsv(ii,jj), jj=2*mxpoi+1,2*mxpoi+nonods)
          !           read(1, 2000) (obsv(ii,jj), jj=3*mxpoi+1,3*mxpoi+nonods)
          !        enddo
          !        close(1) 
          !#############
          !CALL INTERPL_MESH2(ntimemax,7442,7442,21600,4,4,7442,1,  &
          !   obsv(:,1:7442),obsv(:,1*mxpoi+1:1*mxpoi+7442),obsv(:,2*mxpoi+1:2*mxpoi+7442),obsv(:,3*mxpoi+1:3*mxpoi+7442),U0,V0,W0,PHI0)
          !obsv=0.0
          !obsv(:,1:mxpoi) = U0(:,1:mxpoi)
          !obsv(:,1*mxpoi+1:2*mxpoi) = V0(:,1:mxpoi)
          !obsv(:,2*mxpoi+1:3*mxpoi) = W0(:,1:mxpoi)
          !obsv(:,3*mxpoi+1:4*mxpoi) = PHI0(:,1:mxpoi)
          !        open(10, file ='observation1.dat')
          !        do ii=0,ntimemax
          !           write(10, 2000) (obsv(ii,jj), jj=1,istate)
          !        enddo
          !        close(10)
          !STOP 89000
          !##########
       ELSE IF(remesh.ge.1.or.IF_INVERSION_MESH) THEN
          ewrite(3,*) 'read observation data'
          obsv=0.0
          !          iflagobs=0
          smean_obs=0.0
          write(name,'(a,i0)') 'observation-newmesh.dat.',num
          OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
          do ii=0,ntimemax
             read(1, 2000) (obsv(ii,jj), jj=1,nonods)
             read(1, 2000) (obsv(ii,jj), jj=1*mxpoi+1,1*mxpoi+nonods)
             read(1, 2000) (obsv(ii,jj), jj=2*mxpoi+1,2*mxpoi+nonods)
             read(1, 2000) (obsv(ii,jj), jj=3*mxpoi+1,3*mxpoi+nonods)
          enddo
          close(1) 
          U0=0.0
          V0=0.0
          W0=0.0
          PHI0=0.0
          !        open(1, file ='uvw0-newmesh.dat')
          !        do ii=0,ntimemax
          !           read(1, 2000) (U0(ii,jj), jj=1,nonods)
          !           read(1, 2000) (V0(ii,jj), jj=1,nonods)
          !           read(1, 2000) (W0(ii,jj), jj=1,nonods)
          !           read(1, 2000) (PHI0(ii,jj), jj=1,nonods)
          !        enddo
          !        close(1) 
          
          write(name,'(a,i0)') 'xyz-newmesh.dat.',num
          OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
          x=0.0
          y=0.0
          z=0.0
          read(1, 2000) (X(jj), jj=1,nonods)
          read(1, 2000) (Y(jj), jj=1,nonods)
          read(1, 2000) (Z(jj), jj=1,nonods)
          close(1) 
          
          write(name,'(a,i0)') 'shape-newmesh.dat.',num
          OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
          NDGLNO=0
          XONDGL=0
          PNDGLN=0
          read(1,*) TOTELE,NONODS,NLoc
          read(1, *) (NDGLNO(jj), jj=1,TOTELE*NLOC)
          read(1, *) (XONDGL(jj), jj=1,TOTELE*NLOC)
          read(1, *) (PNDGLN(jj), jj=1,TOTELE*MLOC)
          close(1) 
          ewrite(3,*) 'after observation data'
       ENDIF  !IF(remesh.eq.0) THEN
       
       ALLOCATE(snapmatrix_obs(istate,nsnap_obs))
       
       snapmatrix_obs =0.0
       smean_obs =0.0
       ewrite(3,*) 'before buildup observation matrix'
       IF(IF_dualwei.and.remesh.ge.1) THEN
          call buildsnapmat_dualwei(snapmatrix_obs,smean_obs,1,ntimemax,    &
               obsv,iflagobs,1,SnapNDT_obs,dual_wei_space)
       ELSE
          call buildsnapmat(snapmatrix_obs,smean_obs,1,ntimemax,    &
               obsv,iflagobs,1,SnapNDT_obs)
       ENDIF
       
       open(10, file ='snapmatrix_obs.dat')
       do ii=1,istate
          WRITE(10, 2000) (snapmatrix_obs(ii,jj), jj=1,nsnap_obs)
       enddo
       close(10)
       
       open(10, file ='smean_obs.dat')
       write(10, 2000) (smean_obs(i), i=1,istate)
       close(10) 
       
       leftsvd_obs =0.0
       svdval_obs =0.0
       IF(IF_H1) THEN
          CALL H1_SPACE(epsilon,snapmatrix_obs,smean_obs,istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,  &
               nsnap_obs,LEFTSVD_obs,SVDVAL_obs)
       ELSE
          ewrite(3,*) 'calculation of POD vectors for observational data U'
          call snapsvd(mxpoi,nsnap_obs,snapmatrix_obs(1:mxpoi,:),nsvd,nsvd,leftsvd_obs(1:mxpoi,:),svdval_obs(1:nsvd))
          ewrite(3,*) 'calculation of POD vectors for observational data V'
          call snapsvd(mxpoi,nsnap_obs,snapmatrix_obs(1*mxpoi+1:2*mxpoi,:),nsvd,nsvd,leftsvd_obs(1*mxpoi+1:2*mxpoi,:),svdval_obs(nsvd+1:2*nsvd))
          ewrite(3,*) 'calculation of POD vectors for observational data W'
          call snapsvd(mxpoi,nsnap_obs,snapmatrix_obs(2*mxpoi+1:3*mxpoi,:),nsvd,nsvd,leftsvd_obs(2*mxpoi+1:3*mxpoi,:),svdval_obs(2*nsvd+1:3*nsvd))
          call snapsvd(mxpoi,nsnap_obs,snapmatrix_obs(3*mxpoi+1:4*mxpoi,:),nsvd,nsvd,leftsvd_obs(3*mxpoi+1:4*mxpoi,:),svdval_obs(3*nsvd+1:4*nsvd))
          
       ENDIF
       DEALLOCATE(snapmatrix_obs)
       
       ! project the  observational data onto the POD space
       
       DO II = 0,ntimemax
          call fulltopodproj_multivar(nsvd_total,nsvd,nsvd_u,nsvd_phi,nvar,mxpoi,obsv(II,1:istate)-smean_obs(1:istate), &
               leftsvd_obs,coef_multivar_obs(II,:))
       ENDDO
       
       open(1,file='coef_multivar_obs.dat')
       do II=1,nsvd
          write(1,*) 'uu-coef'
          write(1,*) (coef_multivar_obs(JJ,II),JJ=0,ntimemax)
       enddo
       do II=nsvd+1,2*nsvd
          write(1,*) 'vv-coef'
          write(1,*) (coef_multivar_obs(JJ,II),JJ=0,ntimemax)
       enddo
       do II=2*nsvd+1,3*nsvd
          write(1,*) 'ww-coef'
          write(1,*) (coef_multivar_obs(JJ,II),JJ=0,ntimemax)
       enddo
       close(1)
       
       !        CALL OBS_SNAPMATRIX_SMEAN_SVD(istate,mxpoi,ntimemax,nsvd,nsvd_total,        &
       !           nsnap_obs,SnapNDT_obs,remesh,IF_dualwei,                        &
       !           iflagobs,obsv,smean_obs,                                                &
       !           leftsvd_obs,svdval_obs,coef_multivar_obs)
       
    ENDIF !(rerunforward.eq.3.or.rerunforward.eq.2) then
    
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !  MAIN TASK: Start the inversion procedure       !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    
    DO GLOITS = 1,NGLOITS
       
       GLOITS_POD=1
       LINITS =1
       
       IF(IF_REMESH) THEN
          NUM=remesh+1 
       ELSE
          NUM=GLOITS
       ENDIF
       
       IF(IRUNFULLFWD) THEN !true. i.e. to run the full model and get POD vectors
          
          
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
          !           THE FULL FORWARD MODEL                     !
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
          !******************************************!
          !      1. RUN THE FULL FORWARD MODEL       !
          !      2. Setup observation                !
          !      3. Buildup the Matrix               !
          !******************************************!
          
          
          !******************************************!
          ! 1. RUN THE FULL FORWARD MODDEL           !
          !******************************************!
          
          ! prepare for runnning the model
          !--------------------------------
          
300       CONTINUE
          
          ! Establish signal handlers.
          call initialise_signals
          CALL memory_size(NIMEM, NRMEM)
          write(3,*) 'allocate the memory for running the full froward model'  
          ALLOCATE(RMEM(NRMEM))
          ALLOCATE(IMEM(NIMEM))
          
          ! This is a stub for FLUID without event option.
          IEVENT=0
          NENERG=0                        
          NNNODS=0                
          NNNELE=0
          NDELA2=0
          NMATEV=0
          MXNDEV=1
          MXNMTT=1
          EVEFIR=.FALSE.
          FILEVE=filename
          !--------------------------------------------------
          !goto 600           
          IF(GLOITS.EQ.1) THEN
             !*******!!!!!!!!!! CHANGE BACK!!!!!!!!************
             !        IF(.false.) then
             ! The first nonlinear iteration
             !------------------------------
             
             IF(rerunforward.eq.2) then
                write(3,*) 'mxpoi',mxpoi
                
                ! -------Setup Initial condition---------
                open(1,file='initial0.dat')
                read(1,*) cho
                read(1,*) (VARIN(IP),IP=1,mxpoi)
                read(1,*) 
                read(1,*) (VARIN(mxpoi+IP),IP=1,mxpoi)
                read(1,*) cho
                read(1,*) (VARIN(2*mxpoi+IP),IP=1,mxpoi)
                read(1,*) cho
                read(1,*) (VARIN(3*mxpoi+IP),IP=1,mxpoi)
                close(1)
                
                DO IP =1,mxpoi
                   if(.false.) then
                      U0(0,IP)=VARIN(IP)+InitialPert
                      V0(0,IP)=VARIN(mxpoi+IP)+InitialPert
                      W0(0,IP)=VARIN(2*mxpoi+IP)
                      PHI0(0,IP)=VARIN(3*mxpoi+IP)
                   else
                      U0(0,IP)=VARIN(IP)+PertCoef*VARIN(IP)
                      V0(0,IP)=VARIN(mxpoi+IP)+PertCoef*VARIN(mxpoi+IP)
                      W0(0,IP)=VARIN(2*mxpoi+IP)
                      PHI0(0,IP)=VARIN(3*mxpoi+IP)
                   endif
                ENDDO
                open(1,file='initial0Plusperturbation.dat')
                write(1,*) 'uu'
                write(1,*) (U0(0,IP),IP=1,mxpoi)
                write(1,*) 'vv'
                write(1,*) (V0(0,IP),IP=1,mxpoi)
                write(1,*) 'ww'
                write(1,*) (W0(0,IP),IP=1,mxpoi)
                close(1)
             ENDIF ! end rerunforward.eq.2
          ELSE
             ! gloits.gt.1
             !--------------
             !**************************
             !        U0(0,1:mxpoi)=obsv(80,1:mxpoi)
             !        V0(0,1:mxpoi)=obsv(80,1*mxpoi+1:2*mxpoi)
             !        W0(0,1:mxpoi)=obsv(80,2*mxpoi+1:3*mxpoi)
             !        PHI0(0,1:mxpoi)=obsv(80,3*mxpoi+1:4*mxpoi)
             !        goto 100
             !************************
             
             
             ! -------Initial condition from the optimised one at the
             !        last optimisation procedure---------
             open(1,file='initial_gloits.dat')
             read(1,*) cho
             read(1,*) (U0(0,IP),IP=1,mxpoi)
             read(1,*) cho
             read(1,*) (V0(0,IP),IP=1,mxpoi)
             read(1,*) cho
             read(1,*) (W0(0,IP),IP=1,mxpoi)
             read(1,*) cho
             read(1,*) (PHI0(0,IP),IP=1,mxpoi)
             close(1)
100          continue
             VARIN(1:mxpoi)=U0(0,:)
             VARIN(mxpoi+1:2*mxpoi)=V0(0,:)
             VARIN(2*mxpoi+1:3*mxpoi)=W0(0,:)
             VARIN(3*mxpoi+1:4*mxpoi)=PHI0(0,:)
             
          ENDIF !end gloits.eq.1
600       continue
          !RUN FORWARD MODEL 
          !-------------------------------------------------------------
          
          write(3,*) 'before calling FLUID'
          CALL FLUIDS(IEVENT,EVESFL,EVENTD,EVENTT,EVENDE,EVENTV, &
               EVECO1,EVECO2,EVECO3,EVECO4, &
               TEMPR1,TEMPR2,TEMPR3,TEMPR4,BURNUP1, &
               EVEFIR,EVEDT,EVEUDT,EVEDTN, NENERG, NDELA2, &
               NNNODS,NNNELE, &
               FILEVE,EVGDIN, &
               EDENIN,ETEMIN, EDENGA,EGAMD2,EGAMD3, &
               VOILIM, EVTIR1,EVTIR2, EVTIOP,EVSBCY,EVMXOP, &
               NMATEV,MXNMTT,TEMEVE, &
               NDIMI,MXNDEV,NDSIEV, VERSEV,GOTGAS, &
               RMEM, NRMEM, IMEM, NIMEM)
          
          DEALLOCATE( RMEM )
          DEALLOCATE( IMEM )   
          write(3,*) 'end fwd run'
          write(3,*) 'start,end,total cpu time',TIME1,TIME2,TIME2-TIME1
          
          VARIN(1:MXPOI)=U0(0,:)
          VARIN(MXPOI+1:2*MXPOI)=V0(0,:)
          VARIN(2*MXPOI+1:3*MXPOI)=W0(0,:)
          VARIN(3*mxpoi+1:4*MXPOI)= PHI0(0,:)
          
          ! -------END RUNNING THE FULL FORWARD MODEL-----------------------
          
          write(3,*) 'aa'        
          IF(rerunforward.eq.0) THEN
             
             !************************************!
             ! 2. SETUP THE OBSERVATIONAL DATA    !
             !************************************!
             
             call setobsv(ntimemax,u0,v0,w0,phi0,obsv,iflagobs,SnapNDT_obs)
             
             open(10, file ='observation.dat')
             do ii=0,ntimemax
                WRITE(10, 2000) (obsv(ii,jj), jj=1,istate)
             enddo
             close(10) 
             write(3,*)  'after set up the observation'  
             
             open(10, file ='observation_u.dat')
             do jj=1,mxpoi
                WRITE(10, 2000) (obsv(ii,jj), ii=1,ntimemax)
             enddo
             close(10) 
             
             open(10, file ='observation_v.dat')
             do jj=mxpoi+1,2*mxpoi
                WRITE(10, 2000) (obsv(ii,jj), ii=1,ntimemax)
             enddo
             close(10) 
             
             open(10, file ='observation_w.dat')
             do jj=2*mxpoi+1,3*mxpoi
                WRITE(10, 2000) (obsv(ii,jj), ii=1,ntimemax)
             enddo
             close(10) 
             
             open(10, file ='observation_phi.dat')
             do jj=3*mxpoi+1,4*mxpoi
                WRITE(10, 2000) (obsv(ii,jj), ii=1,ntimemax)
             enddo
             close(10) 
             
             !             open(10, file ='iflagobs.dat')
             !           do ii=0,ntimemax
             !               WRITE(10, 2000) (iflagobs(ii,jj), jj=1,mxpoi)
             !           enddo
             !           close(10) 
             
             write(3,*) 'bb'        
             ! setup observation with perturbation !   
             !-------------------------------------!
             
             !     call pertobsv(ntimemax,obsvpert,iflagobs)
             
             !*************************************************************!
             !     3.save the exact initial condition  for POD model       !     
             !*************************************************************!
             open(1,file='initial0.dat')
             write(1,*) 'uu'
             write(1,*) (VARIN(IP),IP=1,mxpoi)
             write(1,*) 'vv'
             write(1,*) (VARIN(mxpoi+IP),IP=1,mxpoi)
             write(1,*) 'ww'
             write(1,*) (VARIN(2*mxpoi+IP),IP=1,mxpoi)
             write(1,*) 'phi'
             write(1,*) (VARIN(3*mxpoi+IP),IP=1,mxpoi)
             close(1)
             stop 1122
             
          ENDIF  !rerunforward.eq.0
          
          !*************************************************************!
          !            4.snapmatrix construction                        !     
          !*************************************************************!
          write(3,*) 'cc'        
          
          ALLOCATE(snapmatrix(istate,nrsnapshots))
          snapmatrix =0.0
          
          write(3,*) 'dd'        
          if(iusemean.eq.2) smean=smean_obs
          
          write(3,*) 'ee'
          IF(IF_dualwei.and.remesh.ge.1) THEN
             call buildsnapmat_dualwei(snapmatrix,smean,iuseobs,ntimemax,    &
                  obsv,iflagobs,iusemean,SnapNdT,dual_wei_space)
          ELSE        
             call buildsnapmat(snapmatrix,smean,iuseobs,ntimemax,    &
                  obsv,iflagobs,iusemean,SnapNdT)
          ENDIF
          
          
          write(name,'(a,i0)') 'snapmatrix.dat.',GLOITS
          OPEN(10,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
          do ii=1,istate
             WRITE(10, 2000) (snapmatrix(ii,jj), jj=1,nrsnapshots)
          enddo
          close(10)
          
          write(name,'(a,i0)') 'varin.dat.',GLOITS
          OPEN(10,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
          write(10, 2000) (varin(i), i=1,istate)
          close(10)
          
          write(name,'(a,i0)') 'smean.dat.',GLOITS
          OPEN(10,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
          write(10, 2000) (smean(i), i=1,istate)
          close(10)
          
          !        stop 90
       ELSE !IRUNFULLFWD=.false. read the results of the full forwarsd model from the files
          !-------------------------------------------------------
          !*************************************************************!
          !                    snapmatrix construction                  !     
          !*************************************************************!
          
          ALLOCATE(snapmatrix(istate,nrsnapshots))
          snapmatrix =0.0
          smean=0.0
          varin=0.0
          IF(rerunforward.eq.3) then
             smean=smean_obs
             write(3,*) 'varin-'
             varin(:)=obsv(0,:)
          ELSE
             
             IF(IFSIMPLE_PODada) THEN
                write(name,'(a,i0)') 'snapmatrix_simplePODada.dat.',GLOITS
                OPEN(10,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                do ii=1,istate
                   read(10, 2000) (snapmatrix(ii,jj), jj=1,nrsnapshots)
                enddo
                close(10)
                
                write(name,'(a,i0)') 'varin_simplePODada.dat.',GLOITS
                OPEN(10,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                read(10, 2000) (varin(i), i=1,istate)
                close(10)
                
                write(name,'(a,i0)') 'smean_simplePODada.dat.',GLOITS
                OPEN(10,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                read(10, 2000) (smean(i), i=1,istate)
                close(10)
                
             ELSE
                write(name,'(a,i0)') 'snapmatrix.dat.',GLOITS
                OPEN(10,FILE=trim(NAME))
                do ii=1,istate
                   READ(10, 2000) (snapmatrix(ii,jj), jj=1,nrsnapshots)
                enddo
                close(10)
                
                write(name,'(a,i0)') 'varin.dat.',GLOITS
                OPEN(10,FILE=trim(NAME))
                READ(10, 2000) (varin(i), i=1,istate)
                close(10)
                
                write(name,'(a,i0)') 'smean.dat.',GLOITS
                OPEN(10,FILE=trim(NAME))
                READ(10, 2000) (smean(i), i=1,istate)
                close(10)
             ENDIF
             
          ENDIF
       ENDIF ! end IRUNFULLFWD
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       !     POD REDUCED FORWARD AND INVERSION MODELS                        !
       !     1. perform the SVD of snapmatrix                                !
       !     2. run the reduced forward mdoel                                !
       !     3. run the reduced adjoint mdoel                                !
       !     4. optimisation procedure                                       !
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       
       
       
       !**********************************************************************!
       !               1. perform the SVD of snapmatrix                       !
       !                     using ARPACK                                     !
       !**********************************************************************!
       
       IF(rerunforward.eq.3) then
          leftsvd=0.0
          leftsvd=leftsvd_obs
       ELSE
          leftsvd=0.0
          svdval=0.0
          IF(IF_H1) THEN
             CALL  H1_SPACE(epsilon,snapmatrix,smean,istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots,LEFTSVD,SVDVAL)
          ELSE
             write(3,*) 'calculation of POD vectors for U'
             call snapsvd(mxpoi,nrsnapshots,snapmatrix(1:mxpoi,:),nsvd,nsvd,leftsvd(1:mxpoi,:),svdval(1:nsvd))
             write(3,*) 'calculation of POD vectors for V'
             call snapsvd(mxpoi,nrsnapshots,                  &
                  snapmatrix(1*mxpoi+1:2*mxpoi,:),nsvd,nsvd,leftsvd(1*mxpoi+1:2*mxpoi,:),svdval(nsvd+1:2*nsvd))
             write(3,*) 'calculation of POD vectors for W'
             call snapsvd(mxpoi,nrsnapshots,                  &
                  snapmatrix(2*mxpoi+1:3*mxpoi,:),nsvd,nsvd,leftsvd(2*mxpoi+1:3*mxpoi,:),svdval(2*nsvd+1:3*nsvd))
             write(3,*) 'calculation of POD vectors for PHI'
             call snapsvd(mxpoi,nrsnapshots,                  &
                  snapmatrix(3*mxpoi+1:4*mxpoi,:),nsvd,nsvd,leftsvd(3*mxpoi+1:4*mxpoi,:),svdval(3*nsvd+1:4*nsvd))
          ENDIF
          IF(IFSIMPLE_PODada) THEN
             ! for 2D???
             leftsvd(2*mxpoi+1:3*mxpoi,:)=leftsvd_obs(3*mxpoi+1:4*mxpoi,:)
          ENDIF
       ENDIF
       !IF(rerunforward.eq.1.or.rerunforward.eq.2) THEN
       !          smean_obs(1:3*mxpoi) = smean(1:3*mxpoi)
       !!          leftsvd_obs(1:3*mxpoi,:) = leftsvd(1:3*mxpoi,:) 
       !          coef_multivar_obs=0.0
       !          isnap =0
       !          DO II = 0,ntimemax             
       !           if (SnapNDT_obs*int(II/SnapNDT_obs) .eq. II) then
       !              isnap = isnap+1
       !              obsv(II,1:3*mxpoi)=(snapmatrix(1:3*mxpoi,isnap)+smean(1:3*mxpoi) )*nrsnapshots**2
       !              leftsvd_obs(1:3*mxpoi,isnap) = leftsvd(1:3*mxpoi,isnap) 
       !!         call fulltopodproj_multivar(3*nsvd,3*nsvd,nsvd,nsvd,3,mxpoi,obsv(II,1:3*mxpoi),leftsvd_obs,coef_multivar_obs(II,:))
       !          call fulltopodproj_multivar(nsvd_total,nsvd,nsvd_u,   &
       !               nsvd_phi,nvar,mxpoi,obsv(II,1:istate)-smean_obs(1:istate),leftsvd_obs,coef_multivar_obs(II,:))
       !           endif
!          ENDDO
       !ENDIF
       
       
       
       
       DEALLOCATE(snapmatrix)
       
       !write leftsv and sv values to output file
       !-----------------------------------------
       
       
       open(10, file ='leftsvincr.dat')
       do ii=1,istate
          write(10, 2000) (leftsvd(ii,jj), jj=1,nsvd)
       enddo
       close(10)

       ! for the geostrophic pressure
       !------------------------------
       
       IF(IFPG.EQ.1) THEN
          
          
          IF(REMESH.EQ.0) THEN
             !  DO  itime = 0,ntimemax
             !    WRITE(FILENAME_MPG, 20) ITIME
             !20 FORMAT('MPG_data',I5.5)
             !  K3 = index( FILENAME_MPG, ' ' ) - 1
             !   FILENAME2=FILENAME_MPG(1:K3)
             !!   write(3,*)  FILENAME2
             !   OPEN(2,FILE = FILENAME2)
             !!    READ(2,*) PGNODS,TOTELE,NGI,MPGLOC,D3
             !     READ(2,*) PGNODS,TOTELE,MPGLOC,D3
             !  CLOSE(2)
             !  ENDDO
             
             OPEN(1,FILE = 'PG_data00001')
             open(2,file='MPG_data00001')
             READ(2,*) PGNODS,TOTELE,MPGLOC,D3
             READ(2,*) NPGCOLM
             NPGCOLM0=NPGCOLM
             
             NGI=11
             write(3,*) 'PGNODS,TOTELE,MPGLOC,D3',PGNODS,TOTELE,MPGLOC,D3
             IF(.NOT.IRUNFULLFWD.OR.GLOITS.EQ.1) THEN
                ALLOCATE(MPG(MPGLOC,NGI))  
                ALLOCATE(MPGLX(MPGLOC,NGI))
                ALLOCATE(MPGLY(MPGLOC,NGI))
                ALLOCATE(MPGLZ(MPGLOC,NGI))
                !      ALLOCATE(PG(0:ntimemax,PGNODS))
                !      ALLOCATE(PG_POD(1,PGNODS))
                ALLOCATE(PGNDGLNO(2*TOTELE*MPGLOC))
                !      ALLOCATE(snapmatrix_PG(MXPOI+2*TOTELE,nrsnapshots))
                ALLOCATE(smean_PG(MXPOI+2*TOTELE))
                ALLOCATE(leftsvd_PG(MXPOI+2*TOTELE,nsvd))
                ALLOCATE(svdval_PG(nsvd))
                !      ALLOCATE(varin_PG(MXPOI+2*TOTELE))
                !      ALLOCATE(varincr_PG(MXPOI+2*TOTELE))
                ALLOCATE(COEF_PG(0:ntimemax,NSVD_PHI))
                ALLOCATE(TCOEF_PG(NSVD_PHI))
                ALLOCATE(PGFINDRM(MXPOI+2*TOTELE+1))
                ALLOCATE(PGCOLM(2*NPGCOLM0))
                ALLOCATE(PGCENTRM(MXPOI+2*TOTELE))
             ENDIF ! end GLOITS.eq.1
             !      PG=0.0
             !      PG_POD=0.0
             MPG=0.0
             MPGLX=0.0
             MPGLY=0.0
             MPGLZ=0.0
             PGNDGLNO=0
             PGFINDRM=0
             PGCOLM=0
             PGCENTRM=0
             read(2,*) (PGFINDRM(I),I=1,PGNODS+1)
             read(2,*) (PGCOLM(I),I=1,NPGCOLM)
             read(2,*) (PGCENTRM(I),I=1,PGNODS)
             READ(1,*) (PGNDGLNO(I),I=1,TOTELE*MPGLOC)
             close(1)
             close(2)
             
          ENDIF !IF(REMESH.EQ.0) THEN
          !      snapmatrix_PG=0.0
          smean_PG=0.0
          leftsvd_PG=0.0
          svdval_PG=0.0
!      varin_PG=0.0
          !      varincr_PG=0.0
          COEF_PG=0.0
          TCOEF_PG=0.0
          
       ENDIF  ! end IFPG
       
       
       if(.false.) then
          ! test to be orthonormal!!
          !--------------------------
          do ii=1,nsvd
             ftest2=0.0
             do jj=1,nsvd
                if(ii.ne.jj) then
                   ftest=0.0
                   do j=1,mxpoi
                      ftest=ftest+ leftsvd(j,ii)*leftsvd(j,jj)
                   enddo
                endif
                if(abs(ftest).gt.1.0e-18) then
                   write(3,*) ii,jj,ftest
                endif
             enddo
             do j=1,mxpoi
                ftest2=ftest2+ leftsvd(j,ii)*leftsvd(j,ii)
             enddo
             write(3,*) 'ftest2',ftest2
          enddo
          stop
       endif
       
       !-------------------------------------------------------
       
       
       !*****************************************************************!
       !       2. project initial conditions on the POD space              !
       !*****************************************************************!
       
       do i=1,istate
          varincr(i) = varin(i) - smean(i)
       enddo

       write(3,*) 'before project initial conditions on the POD space'
       !  call fulltopodproj(istate,nsvd,varincr,leftsvd,tcoef1)
       call fulltopodproj_multivar(nsvd_total,nsvd,nsvd_u,nsvd_phi,nvar,mxpoi,varincr,leftsvd,tcoef_multivar)
       
       COEF_multivar=0.0
       COEF_multivar(0,:)=TCOEF_multivar(:)
       IF(rerunforward.eq.3) then
          if(.false.) then
             open(1,file='coef_multivar_exact.dat')
             do II=1,nsvd_total
                read(1,*) (coef_multivar_obs(JJ,II),JJ=0,ntimemax)
             enddo
             close(1)
          endif
          
          DO I=1,3*NSVD
             COEF_multivar(0,I)=coef_multivar_obs(0,I)
             COEF_multivar(0,I)=COEF_multivar(0,I)+coef_multivar_obs(0,I)*PertCoef
             !        COEF_multivar(0,I)=0.0
          ENDDO
       ENDIF
    
       
       
       !**************************************************************!
       !               3. RUNNING the reduced model                   !
       !**************************************************************!
       IF(REMESH.EQ.0) THEN
          !READ DATA------------------------------------------  
          
          OPEN(1,file='shape.dat')
          
          READ(1,*) TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC,FREDOP,NPRESS
          READ(1,*) D3,DCYL
          write(3,*)  TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC,FREDOP,NPRESS
          write(3,*)  D3,DCYL
          TOTELE0=TOTELE
          IF(.NOT.IRUNFULLFWD.OR.GLOITS.EQ.1) THEN
             ALLOCATE(N(NLOC,NGI))
             ALLOCATE(NLX(NLOC,NGI))
             ALLOCATE(NLY(NLOC,NGI))
             ALLOCATE(NLZ(NLOC,NGI))
             write(3,*) 'after locating N,NLX..'
             ALLOCATE(M(MLOC,NGI))
             ALLOCATE(MLX(MLOC,NGI))
             ALLOCATE(MLY(MLOC,NGI))
             ALLOCATE(MLZ(MLOC,NGI))
             write(3,*) 'after locating M,MLX..'
             ALLOCATE(MUPTXX(MXPOI))
             ALLOCATE(MUPTXY(MXPOI))
             ALLOCATE(MUPTXZ(MXPOI))
             ALLOCATE(MUPTYY(MXPOI))
             ALLOCATE(MUPTYZ(MXPOI))
             ALLOCATE(MUPTZZ(MXPOI))
             write(3,*) 'after locating the viscosities'
             ALLOCATE(WEIGHT(NGI))
             ALLOCATE(X(MXPOI))
             ALLOCATE(Y(MXPOI))
             ALLOCATE(Z(MXPOI))
             write(3,*) 'after locating the locations,x,y,z'
             II=TOTELE*NLOC
             ALLOCATE(NDGLNO(2*II))
             ALLOCATE(XONDGL(2*II))
             II=TOTELE*MLOC
             ALLOCATE(PNDGLN(2*II))
             
             ALLOCATE(FINDRM(MXPOI+1))
             ALLOCATE(CENTRM(MXPOI))
          ENDIF
          
          N=0.0
          NLX=0.0
          NLY=0.0
          NLZ=0.0
          M=0.0
          MLX=0.0
          MLY=0.0
          MLZ=0.0
          MUPTXX=0.0
          MUPTXY=0.0
          MUPTXZ=0.0
          MUPTYY=0.0
          MUPTYZ=0.0
          MUPTZZ=0.0
          WEIGHT=0.0
          X=0.0
          Y=0.0
          Z=0.0
          NDGLNO=0
          XONDGL=0
          PNDGLN=0
          
          write(3,*) 'staring to read data from FWD'
          DO JJ=1,NGI
             READ(1,*) (M(II,JJ),II=1,MLOC)
             READ(1,*) (MLX(II,JJ),II=1,MLOC)
             READ(1,*) (MLY(II,JJ),II=1,MLOC)
             READ(1,*) (MLZ(II,JJ),II=1,MLOC)
          ENDDO
          
          write(3,*) 'finish reading M,MLX...'
          
          DO JJ=1,NGI
             READ(1,*) (N(II,JJ),II=1,NLOC)
             READ(1,*) (NLX(II,JJ),II=1,NLOC)
             READ(1,*) (NLY(II,JJ),II=1,NLOC)
             READ(1,*) (NLZ(II,JJ),II=1,NLOC)
          ENDDO
          write(3,*) 'finish reading N,NLX...'
          
          READ(1,*) (WEIGHT(JJ),JJ=1,NGI)
          write(3,*) 'finish reading weight'
          READ(1,*) (FINDRM(JJ),JJ=1,NONODS+1)
          READ(1,*) (CENTRM(JJ),JJ=1,NONODS)
          
          CLOSE(1)
          
          OPEN(1,file='XYZ_MUPTXX.dat')
          READ(1,*) DT,LTIME
          FLAbort("Broken - Stephan")
          !READ(1,*) OPTOME,GEOBAL
          READ(1,*) SCFACTH0
          write(3,*)  'DT,LTIME',DT,LTIME
          
          READ(1,*) (MUPTXX(JJ),JJ=1,NONODS)
          READ(1,*) (MUPTXY(JJ),JJ=1,NONODS)
          READ(1,*) (MUPTXZ(JJ),JJ=1,NONODS)
          READ(1,*) (MUPTYY(JJ),JJ=1,NONODS)
          READ(1,*) (MUPTYZ(JJ),JJ=1,NONODS)
          READ(1,*) (MUPTZZ(JJ),JJ=1,NONODS)
          write(3,*) 'finishing reading viscosity'
          READ(1,*) (X(JJ),JJ=1,NONODS)
          READ(1,*) (Y(JJ),JJ=1,NONODS)
          READ(1,*) (Z(JJ),JJ=1,NONODS)
          write(3,*) 'finish reading locations, x,y,z'
          READ(1,*) (NDGLNO(JJ),JJ=1,TOTELE*NLOC)
          READ(1,*) (XONDGL(JJ),JJ=1,TOTELE*NLOC)
          READ(1,*) (PNDGLN(JJ),JJ=1,TOTELE*MLOC)
          
          CLOSE(1)
          !fix it later if taking account into free surface 
          IF(IF_REMESH) THEN 
             OPEN(1,file='surface-shape.dat')
             read(1,*) stotel,snloc
             stotel0=stotel
             II=STOTEL0*SNLOC
             ALLOCATE(SNDGLN(2*II))
             SNDGLN=0
             read(1,*) (sndgln(II),II=1,stotel*snloc)
             CLOSE(1)
          ENDIF
          
          
       ENDIF !IF(REMESH.EQ.0) THEN
       
500    CONTINUE
       
       ! run the reduced model--------------
       !------------------------------------
       
       write(3,*) 'before running the reduced model'
       !###### test gradient#####
       !!open(8,file='FwdGradient.dat')  
       !!open(9,file='FwdGradient1.dat') 
       !!DO I =1,NSVD_TOTAL
       !!   COEF_multivar(0,:)=coef_multivar_obs(0,:)
       !!   COEF_multivar(0,I)=COEF_multivar(0,I)+coef_multivar_obs(0,I)*PertCoef
       !###### test gradient#####
       
       if(IFPG.EQ.1) then
          ALLOCATE(leftsvd_PG_x(MXPOI+2*TOTELE0,nsvd_phi))
          ALLOCATE(leftsvd_PG_y(MXPOI+2*TOTELE0,nsvd_phi))
          leftsvd_PG_x=0.0
          leftsvd_PG_y=0.0
          write(3,*) PGNODS,MXPOI+2*TOTELE
          if(PGNODS.gt.MXPOI+2*TOTELE) stop
          CALL REDUCED_MODEL_DIFFCOEF1_PG1(NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,NONODS,XNONOD,&
               X(1:NONODS),Y(1:NONODS),Z(1:NONODS),    &
               LEFTSVD,SMEAN,COEF_multivar,M,MLX,MLY,MLZ,                                    &
               N,NLX,NLY,NLZ,NDGLNO(1:TOTELE*NLOC),XONDGL(1:TOTELE*NLOC),PNDGLN(1:TOTELE*MLOC),WEIGHT,TOTELE,NLOC,NGI,MLOC,       &
               MUPTXX(1:NONODS),MUPTXY(1:NONODS),MUPTXZ(1:NONODS),MUPTYY(1:NONODS),MUPTYZ(1:NONODS),MUPTZZ(1:NONODS),             &
               NTIMEMAX,DT,D3,DCYL,GEOBAL,SCFACTH0,IFWIND, &
               PGNODS,MPGLOC,MPG,MPGLX,MPGLY,MPGLZ,PGNDGLNO(1:TOTELE*MPGLOC),SMEAN_PG(1:PGNODS),leftsvd_PG_x(1:PGNODS,1:NSVD_PHI), &
               leftsvd_PG_y(1:PGNODS,1:NSVD_PHI),  &
               NPGCOLM,PGFINDRM(1:PGNODS+1),PGCOLM(1:NPGCOLM),PGCENTRM(1:PGNODS),IFWRITE,NUM,LINITS)
          !####
          IF(LINITS.NE.1) THEN
             DEALLOCATE(leftsvd_PG_x)
             DEALLOCATE(leftsvd_PG_y)
          ENDIF
          !!stop 333
       else
          
          CALL REDUCED_MODEL_DIFFCOEF1(NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,NONODS,XNONOD,&
               X(1:NONODS),Y(1:NONODS),Z(1:NONODS),        &
               LEFTSVD,SMEAN,COEF_multivar,M,MLX,MLY,MLZ,                                    &
               N,NLX,NLY,NLZ,NDGLNO(1:TOTELE*NLOC),XONDGL(1:TOTELE*NLOC),PNDGLN(1:TOTELE*MLOC),WEIGHT,TOTELE,NLOC,NGI,MLOC,       &
               MUPTXX(1:NONODS),MUPTXY(1:NONODS),MUPTXZ(1:NONODS),MUPTYY(1:NONODS),MUPTYZ(1:NONODS),MUPTZZ(1:NONODS),   &
               NTIMEMAX,DT,D3,DCYL,GEOBAL,SCFACTH0,IFWIND)
       endif
       
       IF(.false.) THEN
          open(1,file='coef_multivar_exact.dat')
          do II=1,nsvd_total
             write(1,*) (coef_multivar(JJ,II),JJ=0,ntimemax)
          enddo
          close(1)
          !        stop 589
       ENDIF
       write(name,'(a,i0)') 'coef_multivar.dat',LINITS
       OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
       do II=1,nsvd_total
          write(1,*) (coef_multivar(JJ,II),JJ=0,ntimemax)
       enddo
       close(1)
       
       write(name,'(a,i0)') 'coef_multivar.dat_GLOITS_POD',GLOITS_POD
       OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
       do II=1,nsvd_total
          write(1,*) (coef_multivar(JJ,II),JJ=0,ntimemax)
       enddo
       close(1)
       write(3,*) 'quit running the reduced model'
       !**************************************************************
       
       
       
       
       FUNCT =0.0
       NF=0
       
       DO NSTEP = 0,ntimemax
          tcoef_multivar(:)=coef_multivar(NSTEP,:)
          !     write(3,*) 'NSTEP =',NSTEP
          !****************project back from rom to full space
          !     call podtofullproj(istate,nsvd,varincr,leftsvd,tcoef1)
          call podtofullproj_multivar(nsvd_total,nsvd,nsvd_u,nsvd_phi,nvar,mxpoi,varincr,leftsvd,tcoef_multivar(:))
          
          DO IP=1,mxpoi        
             U0(NSTEP,IP)=varincr(IP)+ smean(IP)
             V0(NSTEP,IP)=varincr(mxpoi+IP)+smean(mxpoi+IP)
             W0(NSTEP,IP)=varincr(2*mxpoi+IP)+smean(2*mxpoi+IP)
             PHI0(NSTEP,IP)=varincr(3*mxpoi+IP)+smean(3*mxpoi+IP)
          ENDDO
          
          
          !***************************************!
          ! 4. Calculate the objective function   !
          !***************************************!
          write(3,*) 'calculate the objective function'
          ml=0.0
          CALL CALCULATE_ML(ML,NONODS,XNONOD,xondgl,                           &
               TOTELE,NLOC,NGI,                                    &
               NDGLNO,X,Y,Z,                                       &
               N,NLX,NLY,NLZ,WEIGHT,D3,DCYL)                      
          
          !        CALL FUNCTION_POD(FUNCT,NONODS,LTIME,DT,NTIMEMAX,NSTEP,NSTEP*DT,   &
          !             U0(NSTEP,:),V0(NSTEP,:),W0(NSTEP,:),U0(NSTEP,:),X,Y,Z,        &
          !             TOTELE,NLOC,NGI,                                              &
          !             NDGLNO,XNONOD,xondgl,                                         &
          !             N,NLX,NLY,NLZ,WEIGHT,D3,DCYL)
          
          IF(.false.) THEN
             DO II =1,nonods
                IF(nstep.eq.0) THEN
                   FUNCT =FUNCT+ 0.5*(U0(NSTEP,II)-VARIN(II))**2+                 &
                        0.5*(V0(NSTEP,II)-VARIN(mxpoi+II))**2
                   !                +0.5*(W0(NSTEP,II)-VARIN(2*mxpoi+II))**2
                   !                 write(3,*) 'U0(NSTEP,II),OBSV(NSTEP,II)',U0(NSTEP,II),OBSV(NSTEP,II)
                   !                 write(3,*) 'V0(NSTEP,II),OBSV(NSTEP,1*mxpoi+II)',V0(NSTEP,II),OBSV(NSTEP,1*mxpoi+II)
                   NF=NF+1
                ELSE
                   IF(iflagobs(NSTEP,II).eq.1) THEN
                      FUNCT =FUNCT+ 0.5*(U0(NSTEP,II)-OBSV(NSTEP,II))**2+                 &
                           0.5*(V0(NSTEP,II)-OBSV(NSTEP,1*mxpoi+II))**2
                      !                +0.5*(W0(NSTEP,II)-OBSV(NSTEP,2*mxpoi+II))**2
                      !                 write(3,*) 'U0(NSTEP,II),OBSV(NSTEP,II)',U0(NSTEP,II),OBSV(NSTEP,II)
                      !                 write(3,*) 'V0(NSTEP,II),OBSV(NSTEP,1*mxpoi+II)',V0(NSTEP,II),OBSV(NSTEP,1*mxpoi+II)
                      NF=NF+1
                   ENDIF
                ENDIF
             ENDDO
             
          ELSE
             DO II =1,nonods
                IF(nstep.eq.0) THEN
                   IF(IFBG_FUN) THEN
                      FUNCT =FUNCT+ 0.5*(U0(NSTEP,II)-VARIN(II))**2+                 &
                           0.5*(V0(NSTEP,II)-VARIN(mxpoi+II))**2
                      !                +0.5*(W0(NSTEP,II)-VARIN(2*mxpoi+II))**2
                      !                 write(3,*) 'U0(NSTEP,II),OBSV(NSTEP,II)',U0(NSTEP,II),OBSV(NSTEP,II)
                      !                 write(3,*) 'V0(NSTEP,II),OBSV(NSTEP,1*mxpoi+II)',V0(NSTEP,II),OBSV(NSTEP,1*mxpoi+II)
                   ENDIF
                ELSE
                   
                   !            IF (NO_obs*SnapNDT_obs*int(nstep/(NO_obs*SnapNDT_obs)).eq.nstep) then
                   !            IF (NO_obs*SnapNDT_obs*int(nstep/(NO_obs*SnapNDT_obs)).eq.nstep.and.nstep.ge.int(ntimemax/2)) then
                   IF (NO_obs*SnapNDT_obs*int(nstep/(NO_obs*SnapNDT_obs)).eq.nstep.and.nstep.gt.1) then
                      !            IF (NO_obs*SnapNDT_obs*int(nstep/(NO_obs*SnapNDT_obs))        &
                      !            .eq.nstep.or.nstep.eq.(ntimemax-1)) then  
                      !            IF(iflagobs(NSTEP,II).eq.1) THEN
                      U_num=smean(II)
                      V_num=smean(mxpoi+II)
                      U_obs=smean_obs(II)
                      V_obs=smean_obs(mxpoi+II)
                      DO KK=1,NSVD
                         U_num=U_num+coef_multivar(NSTEP,kk)*LEFTSVD(II,KK)
                         V_num=V_num+coef_multivar(NSTEP,NSVD+kk)*LEFTSVD(mxpoi+II,KK)
                         U_obs=U_obs+coef_multivar_obs(NSTEP,KK)*leftsvd_obs(II,KK)
                         V_obs=V_obs+coef_multivar_obs(NSTEP,NSVD+KK)*leftsvd_obs(mxpoi+II,KK)
                      ENDDO
                      IF(NSTEP.NE.ntimemax) THEN
                         FUNCT =FUNCT+0.5*(U_num-U_obs)**2+0.5*(V_num-V_obs)**2*ML(II)
                         NF=NF+1
                      ENDIF
                      !!              NF=NF+1
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
          
          ! OUTPUT THE DATA FOR PLOTTER
          !----------------------------
          if (20*int(nstep/20) .eq. nstep) then
             Call VTKOUTPUT1P(TOTELE,NONODS,1,NLoc,                               &
                  D3,                                                             &
                  NDGLNO,                                                         &
                  X,Y,Z,                                                          &
                  U0(NSTEP,:),V0(NSTEP,:),W0(NSTEP,:),PHI0(NSTEP,:),U0(NSTEP,:),  &
                  1,                                                              &
                  Nstep, 'Fwdpod ')
             error_u(:)=U0(NSTEP,:)-OBSV(NSTEP,1:mxpoi)
             error_v(:)=V0(NSTEP,:)-OBSV(NSTEP,mxpoi+1:2*mxpoi)
             error_w(:)=W0(NSTEP,:)-OBSV(NSTEP,2*mxpoi+1:3*mxpoi)
             error_phi(:)=PHI0(NSTEP,:)-OBSV(NSTEP,3*mxpoi+1:4*mxpoi)
             
             Call VTKOUTPUT1P(TOTELE,NONODS,1,NLoc,                               &
                  D3,                                                             &
                  NDGLNO,                                                         &
                  X,Y,Z,                                                          &
                  error_u(:),error_v(:),error_w(:),error_phi(:),error_u(:),       &
                  1,                                                              &
                  Nstep, 'Errorpod ')
             
          endif
          
          !              CALL OUTPUT_POD(NSTEP,NONODS,XNONOD,TOTELE0,TOTELE,NLOC,X,Y,Z,NDGLNO)
          write(3,*) 'FUNCT=',FUNCT/(mxpoi),nf
       ENDDO
       
       !  write(3,*) 'funct1=', FUNCT,mxpoi,ntimemax
       FUNCT=FUNCT/(nonods)
       !  FUNCT=FUNCT/(NF-1)
       ! ####### test the gradient#########
       !!  Write(8,*)'DCOEF=',COEF_multivar(0,i),coef_multivar_obs(0,i), &
       !!  COEF_multivar(0,i)-coef_multivar_obs(0,i),FUNCT,FUNCT/(COEF_multivar(0,i)-coef_multivar_obs(0,i))
       !!  write(3,*) 'funct2=',FUNCT
       !!  write(9,*) FUNCT/(COEF_multivar(0,i)-coef_multivar_obs(0,i))
       !!ENDDO
       !!   close(8)
       !!   close(9)
       !   stop 890
       !########
       !***************************
       CALL ERROR_ANALYSIS(istate,ntimemax,mxpoi,nrsnapshots,TOTELE0,TOTELE,NONODS,NLoc,D3,NDGLNO,  &
            X,Y,Z,obsv,u0,v0,w0,phi0,DT,T0,num)
       
       !***************************
       
       
       write(3,*) FUNCT
       IF(remesh.EQ.1) THEN
          write(3,*)  'Completing Reduced order simulation with goal-based error measure'
          goto 3000
       ENDIF
       IF(IF_H1) THEN
          stop 199
       ENDIF
       
       IF(LINITS.EQ.1) THEN
          
          !**************************************************!
          !                 RUN THE ADJOINT MODEL            !
          !**************************************************!
          
          ! Initial condition for running the adjoint model
          !------------------------------------------------
          
          !coef_multivar_obs=coef_multivar
          COEF_multivar_ADJ(0,:)=0.0
          write(3,*) 'before running the reduced adjoint model'
          IF(IFPG.EQ.1) THEN
             IF(IF_REMESH_OBS) THEN
                CALL REDUCED_MODEL_DIFFCOEF1consistent_PG_ADJ(NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,NONODS,XNONOD,&
                     X(1:NONODS),Y(1:NONODS),Z(1:NONODS),  &
                     LEFTSVD,SMEAN,coef_multivar,M,MLX,MLY,MLZ,                                                         &
                     N,NLX,NLY,NLZ,NDGLNO(1:TOTELE*NLOC),XONDGL(1:TOTELE*NLOC),PNDGLN(1:TOTELE*MLOC),WEIGHT,TOTELE,NLOC,NGI,MLOC,  &
                     MUPTXX(1:NONODS),MUPTXY(1:NONODS),MUPTXZ(1:NONODS),MUPTYY(1:NONODS),MUPTYZ(1:NONODS),MUPTZZ(1:NONODS),  &
                     NTIMEMAX,DT,D3,DCYL,GEOBAL,SCFACTH0,coef_multivar_ADJ,    &
                     obsv,iflagobs,istate,coef_multivar_obs,leftsvd_obs,smean_obs,NO_obs*SnapNDT_obs,IFWIND,&
                     PGNODS,MPGLOC,MPG,MPGLX,MPGLY,MPGLZ,PGNDGLNO(1:TOTELE*MPGLOC),SMEAN_PG(1:PGNODS),leftsvd_PG_x(1:PGNODS,1:NSVD_PHI), &
                     leftsvd_PG_y(1:PGNODS,1:NSVD_PHI),  &
                     NPGCOLM,PGFINDRM(1:PGNODS+1),PGCOLM(1:NPGCOLM),PGCENTRM(1:PGNODS),NOCVA,G,IFWRITE,GLOITS,LINITS,IF_REMESH_OBS,IF_REMESH_VORTICITY)
             ELSE
                allocate(coef_multivar_obs1(0:ntimemax,nsvd_total))
                allocate(leftsvd_obs1(istate,nsvd))
                allocate(smean_obs1(istate))
                coef_multivar_obs1=0.0
                leftsvd_obs1=0.0
                smean_obs1=0.0
                CALL REDUCED_MODEL_DIFFCOEF1consistent_PG_ADJ(NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,NONODS,XNONOD,&
                     X(1:NONODS),Y(1:NONODS),Z(1:NONODS),  &
                     LEFTSVD,SMEAN,coef_multivar,M,MLX,MLY,MLZ,                                                         &
                     N,NLX,NLY,NLZ,NDGLNO(1:TOTELE*NLOC),XONDGL(1:TOTELE*NLOC),PNDGLN(1:TOTELE*MLOC),WEIGHT,TOTELE,NLOC,NGI,MLOC,  &
                     MUPTXX(1:NONODS),MUPTXY(1:NONODS),MUPTXZ(1:NONODS),MUPTYY(1:NONODS),MUPTYZ(1:NONODS),MUPTZZ(1:NONODS),  &
                     NTIMEMAX,DT,D3,DCYL,GEOBAL,SCFACTH0,coef_multivar_ADJ,    &
                     obsv,iflagobs,istate,coef_multivar_obs1,leftsvd_obs1,smean_obs1,1,IFWIND,&
                     PGNODS,MPGLOC,MPG,MPGLX,MPGLY,MPGLZ,PGNDGLNO(1:TOTELE*MPGLOC),SMEAN_PG(1:PGNODS),leftsvd_PG_x(1:PGNODS,1:NSVD_PHI), &
                     leftsvd_PG_y(1:PGNODS,1:NSVD_PHI),  &
                     NPGCOLM,PGFINDRM(1:PGNODS+1),PGCOLM(1:NPGCOLM),PGCENTRM(1:PGNODS),NOCVA,G,IFWRITE,NUM,LINITS,IF_REMESH_OBS,IF_REMESH_VORTICITY)
                deallocate(coef_multivar_obs1)
                deallocate(leftsvd_obs1)
                deallocate(smean_obs1)
                
             ENDIF
          ELSE
             CALL REDUCED_MODEL_DIFFCOEF1consistent_ADJ(NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,NONODS,XNONOD,&
                  X(1:NONODS),Y(1:NONODS),Z(1:NONODS),  &
                  LEFTSVD,SMEAN,coef_multivar,M,MLX,MLY,MLZ,                                                         &
                  N,NLX,NLY,NLZ,NDGLNO(1:TOTELE*NLOC),XONDGL(1:TOTELE*NLOC),PNDGLN(1:TOTELE*MLOC),WEIGHT,TOTELE,NLOC,NGI,MLOC,  &
                  MUPTXX(1:NONODS),MUPTXY(1:NONODS),MUPTXZ(1:NONODS),MUPTYY(1:NONODS),MUPTYZ(1:NONODS),MUPTZZ(1:NONODS),   &
                  NTIMEMAX,DT,D3,DCYL,GEOBAL,SCFACTH0,coef_multivar_ADJ,    &
                  obsv,iflagobs,istate,coef_multivar_obs,leftsvd_obs,smean_obs,SnapNDT_obs,IFWIND,NOCVA,G)
          ENDIF
          
          IF(LINITS.EQ.1) THEN
             IF(.NOT.IF_REMESH) THEN
                DEALLOCATE(leftsvd_PG_x)
                DEALLOCATE(leftsvd_PG_y)
             ENDIF
          ENDIF
          
          write(name,'(a,i0)') 'coef_multivar_ADJ.dat',GLOITS_POD
          OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
          do II=1,nsvd_total
             write(1,*) (coef_multivar_ADJ(JJ,II),JJ=0,ntimemax)
          enddo
          close(1)
          open(1,file='gradientcheck.dat')
          do II=1,nsvd_total
             write(1,*) g(ii)
          enddo
          close(1)
          
          ! OUTPUT THE DATA FOR PLOTTER
          !----------------------------
          
          ALLOCATE(U_ADJ(mxpoi))
          ALLOCATE(V_ADJ(mxpoi))
          ALLOCATE(W_ADJ(mxpoi))
          ALLOCATE(varincr_adj(istate))
          
          U_ADJ=0.0
          V_ADJ=0.0
          W_ADJ=0.0
          varincr_adj=0.0
          
          IF(IF_remesh) THEN
             allocate(U0_ADJ(0:ntimemax,mxpoi))
             allocate(V0_ADJ(0:ntimemax,mxpoi))
             allocate(W0_ADJ(0:ntimemax,mxpoi))
             allocate(PHI0_ADJ(0:ntimemax,mxpoi))
             U0_ADJ=0.0
             V0_ADJ=0.0
             W0_ADJ=0.0
             PHI0_ADJ=0.0
          ENDIF
          
          DO NSTEP = 0,ntimemax
             tcoef_multivar(:)=coef_multivar_ADJ(NSTEP,:)
             !     write(3,*) 'NSTEP =',NSTEP
             !****************project back from rom to full space
             call podtofullproj_multivar(nsvd_total,nsvd,nsvd_u,nsvd_phi,nvar,mxpoi,varincr_adj,leftsvd,tcoef_multivar(:))
             
             DO IP=1,mxpoi        
                U_ADJ(IP)=varincr_adj(IP)
                V_ADJ(IP)=varincr_adj(mxpoi+IP)
                W_ADJ(IP)=varincr_adj(2*mxpoi+IP)
             ENDDO
             
             IF(IF_remesh) THEN
                U0_ADJ(NSTEP,:)=U_ADJ(:)
                V0_ADJ(NSTEP,:)=V_ADJ(:)
                W0_ADJ(NSTEP,:)=W_ADJ(:)
                PHI0_ADJ(NSTEP,:)=varincr_adj(3*mxpoi+IP)
             ENDIF
             ! OUTPUT THE DATA FOR PLOTTER
             !----------------------------
             if (20*int(nstep/20) .eq. nstep) then
                Call VTKOUTPUT1P(TOTELE,NONODS,1,NLoc,                               &
                     D3,                                                             &
                     NDGLNO,                                                         &
                     X,Y,Z,                                                          &
                     U_ADJ(:),V_ADJ(:),W_ADJ(:),U_ADJ(:),U_ADJ(:),  &
                     2,                                                              &
                     Nstep, 'Adjpod ')
             endif
             
          ENDDO
          
          DEALLOCATE(U_ADJ)
          DEALLOCATE(V_ADJ)
          DEALLOCATE(W_ADJ)
          DEALLOCATE(varincr_adj)
          
          
          
          IF(IF_remesh) THEN
             write(3,*) 'before goal-based error measure'
             CALL errormeasure_adaptmesh(NSVD_TOTAL,NSVD,NSVD_U,NSVD_PHI,MXPOI,NONODS,XNONOD,FREDOP,NPRESS,&
                  X,Y,Z,                       &
                  U0_ADJ,V0_ADJ,W0_ADJ,PHI0_ADJ,  &
                  obsv(0:ntimemax,1:mxpoi),obsv(0:ntimemax,mxpoi+1:2*mxpoi),obsv(0:ntimemax,2*mxpoi+1:3*mxpoi),obsv(0:ntimemax,3*mxpoi+1:4*mxpoi), &
                  FINDRM,CENTRM,M,MLX,MLY,MLZ,  &
                  N,NLX,NLY,NLZ,NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE0,TOTELE,NLOC,NGI,MLOC,       &
                  MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,                            &
                  NTIMEMAX,DT,D3,DCYL,GEOBAL,SCFACTH0,coef_multivar_ADJ, &
                  obsv,iflagobs,COEF_multivar,coef_multivar_obs,LEFTSVD,SMEAN,leftsvd_obs,smean_obs,SnapNDT_obs,IFWIND,  &
                  PGNODS,MPGLOC,MPG,MPGLX,MPGLY,MPGLZ,PGNDGLNO,SMEAN_PG,leftsvd_PG_x,leftsvd_PG_y,  &
                  NPGCOLM0,NPGCOLM,PGFINDRM,PGCOLM,PGCENTRM,DF,IF_REMESH_OBS,IF_INTERPOL,remesh,  &
                  STOTEL0,STOTEL,SNLOC,SNDGLN,NUM,IF_REMESH_VORTICITY)
             DEALLOCATE(leftsvd_PG_x)
             DEALLOCATE(leftsvd_PG_y)
             
             
             DEALLOCATE(U0_ADJ)
             DEALLOCATE(V0_ADJ)
             DEALLOCATE(W0_ADJ)
             DEALLOCATE(PHI0_ADJ)
             !???????
             IF(IF_INVERSION_MESH) THEN
                IFWRITE=.true.
                IF_REMESH=.false.
                goto 8000
             ELSE
                remesh=remesh+1
                IFWRITE=.true.
                IF_INTERPOL=.true.
                IF(IF_DUALWEI) IF_REMESH=.false.
                goto 8000
             ENDIF
          ENDIF
       ENDIF  !IF (LINITS.EQ.1) THEN
       
       !************************************************!
       !    Start optimization procedure                !
       !************************************************!
       
       ! Method 1:  the Nonlinear CG 
       !-----------------------------
       write(3,*) 'NONLINCG'
       ! Calculation of controls and gradient
       IF(NOCVA.NE.NSVD_TOTAL) THEN
          write(3,*) 'NOCVA.NE.NSVD_TOTAL'
          stop 233
       ENDIF
       CVA(1:3*nsvd) = COEF_multivar(0,1:3*nsvd)
       IF(.FALSE.) THEN
          IF(LINITS.EQ.1) THEN
             G(1:3*nsvd) = coef_multivar_ADJ(NTIMEMAX,1:3*nsvd)
          ENDIF
       ENDIF
       
       open(1,file='g1.dat')
       write(1,*) (g(ii),II=1,nocva)
       close(1)
       
       IF(LINITS.EQ.1) THEN
          DO II = 1,NOCVA
             G(II)=-COF*G(II)
          ENDDO
       ENDIF
       IF(GLOITS.eq.1.and.LINITS.EQ.1) THEN
          write(3,*) 'VALFUN=FUNCT',FUNCT
       ENDIF
       
       IF(LINITS.EQ.1.AND.GLOITS.NE.1) THEN
        DO II = 1,NOCVA
           GOLD(II)=-COF*GOLD(II)
        ENDDO
     ENDIF
     open(1,file='g2.dat')
     write(1,*) (g(ii),II=1,nocva)
     close(1)
     
     CALL NONLINCG(GLOITS_POD,LINITS,FUNCT,CVA,G,                      &
          GOLDR,LINCONVER,NOCVA,GOLD,D,CVAOLD,0.0)
     
     coef_multivar(0,1:3*nsvd)=CVA(1:3*nsvd)
     
     OPEN(1,file='coef1.dat')
     write(1,*) (coef_multivar(0,II),II=1,3*nsvd)
     CLOSE(1)
     ! Is the linear search convergent?
     IF(.NOT.LINCONVER) THEN
        OPEN(1,FILE='function.dat',POSITION='APPEND')
        WRITE(1,*) GLOITS,LINITS,FUNCT
        CLOSE(1)
        GOTO 500
     ENDIF
     
     MAXGRAD=0.0
     MAXGRAD_AVE=0.0
     DO II = 1,NOCVA
        IF(ABS(G(II)/COF).GT.MAXGRAD) MAXGRAD=ABS(G(II)/COF)
        MAXGRAD_AVE=MAXGRAD_AVE+G(II)/COF
     ENDDO
     MAXGRAD_AVE=MAXGRAD_AVE/NOCVA
     
     ! Is the nonlinear conjugate convergent?
     OPEN(1,FILE='function.dat',POSITION='APPEND')
     OPEN(2,FILE='gradient.dat',POSITION='APPEND')
     OPEN(4,FILE='gradient-ave.dat',POSITION='APPEND')
     WRITE(1,*) GLOITS,GLOITS_POD,FUNCT
     WRITE(2,*) GLOITS,GLOITS_POD,MAXGRAD
     WRITE(4,*) GLOITS,GLOITS_POD,MAXGRAD_AVE
     CLOSE(1)
     CLOSE(2)
     CLOSE(4)
     write(3,*) 'GLOITS,FUNCT',GLOITS,FUNCT
     write(3,*) 'GLOITS,MAXGRAD',GLOITS,MAXGRAD
     IF(COF*MAXGRAD.LT.CONJCONVER) THEN
        GOTO 1500
        !     ELSE IF(GLOITS_POD.LT.2) THEN
        !             LINITS=1
        !             GLOITS_POD=GLOITS_POD+1
        !             IFWRITE=.false.
        !             GOTO 500
     ELSE
        IRUNFULLFWD =.true.
        open(1,file='initial_gloits.dat')
        write(1,*) 'uu'
        write(1,*) (U0(0,IP),IP=1,mxpoi)
        write(1,*) 'vv'
        write(1,*) (V0(0,IP),IP=1,mxpoi)
        write(1,*) 'ww'
        write(1,*) (W0(0,IP),IP=1,mxpoi)
        write(1,*) 'phi'
        write(1,*) (PHI0(0,IP),IP=1,mxpoi)
        close(1)
        
        ALLOCATE(snapmatrix(istate,nrsnapshots))
        snapmatrix =0.0
        
        if(iusemean.eq.2) smean=smean_obs
        
        IF(IF_dualwei.and.remesh.ge.1) THEN
           call buildsnapmat_dualwei(snapmatrix,smean,iuseobs,ntimemax,    &
                obsv,iflagobs,iusemean,SnapNdT,dual_wei_space)
        ELSE        
           call buildsnapmat(snapmatrix,smean,iuseobs,ntimemax,    &
                obsv,iflagobs,iusemean,SnapNdT)
        ENDIF
        
        write(name,'(a,i0)') 'snapmatrix_simplePODada.dat.',GLOITS
        OPEN(10,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
        do ii=1,istate
           WRITE(10, 2000) (snapmatrix(ii,jj), jj=1,nrsnapshots)
        enddo
        close(10)
        
        write(name,'(a,i0)') 'varin_simplePODada.dat.',GLOITS
        OPEN(10,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
        write(10, 2000) (varin(i), i=1,istate)
 close(10)
 
 write(name,'(a,i0)') 'smean_simplePODada.dat.',GLOITS
 OPEN(10,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
 write(10, 2000) (smean(i), i=1,istate)
 close(10)
 DEALLOCATE(snapmatrix)
ENDIF

IRUNFULLFWD=.true.
IFWRITE=.TRUE.
ENDDO ! GLOITS = 1,NGLOITS




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!  END the inversion procedure                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

2000 format (20(1X,E24.16)) 
200 format(1X,E24.16)

1500 write(3,*) 'NONLINEAR CONJUGATE GREDIANT CONVERGENT'
3000 CONTINUE
DEALLOCATE(obsv)
DEALLOCATE(iflagobs)
DEALLOCATE(smean_obs)


DEALLOCATE(CVA)
DEALLOCATE(CVAOLD)
DEALLOCATE(G)
DEALLOCATE(GOLD)
DEALLOCATE(D)


DEALLOCATE(N)
DEALLOCATE(NLX)
DEALLOCATE(NLY)
DEALLOCATE(NLZ)
DEALLOCATE(M)
DEALLOCATE(MLX)
DEALLOCATE(MLY)
DEALLOCATE(MLZ)
DEALLOCATE(MUPTXX)
DEALLOCATE(MUPTXY)
DEALLOCATE(MUPTXZ)
DEALLOCATE(MUPTYY)
DEALLOCATE(MUPTYZ)
DEALLOCATE(MUPTZZ)
DEALLOCATE(WEIGHT)
DEALLOCATE(X)
DEALLOCATE(Y)
DEALLOCATE(Z)
DEALLOCATE(NDGLNO)
DEALLOCATE(XONDGL)
DEALLOCATE(PNDGLN)
DEALLOCATE(SNDGLN)
DEALLOCATE(FINDRM)
DEALLOCATE(CENTRM)

IF(IFPG.EQ.1) THEN
   DEALLOCATE(MPG)  
   DEALLOCATE(MPGLX)
   DEALLOCATE(MPGLY)
   DEALLOCATE(MPGLZ)
   !      DEALLOCATE(PG)
   !      DEALLOCATE(PG_POD)
   DEALLOCATE(PGNDGLNO)
   !      DEALLOCATE(snapmatrix_PG)
   DEALLOCATE(smean_PG)
   DEALLOCATE(leftsvd_PG)
   DEALLOCATE(svdval_PG)
   !      DEALLOCATE(varin_PG)
   !      DEALLOCATE(varincr_PG)
   
   if(.true.) then
      DEALLOCATE(PGFINDRM)
      DEALLOCATE(PGCOLM)
      DEALLOCATE(PGCENTRM)
   endif
   IF(IF_dualwei) THEN
      deallocate(dual_wei_space)
   ENDIF

ENDIF

#else
FLExit("No ARPACK support, recompile with ARPACK - cannot run reduced model.")
#endif

RETURN
END SUBROUTINE ADreducedFLUIDS_Initial

!#######################################################      

SUBROUTINE GET_PROJECT_NAME(filename, lenth)
  use FLDebug
  IMPLICIT NONE
  INTEGER, INTENT(IN)::lenth
  CHARACTER(LEN=lenth), INTENT(INOUT)::filename 
  CHARACTER(LEN=4096) buffer
  INTEGER ProcNo,IERR
#ifdef HAVE_MPI
  INCLUDE 'mpif.h'
#endif
  ewrite(3, *) "SUBROUTINE GET_PROJECT_NAME()"
  CALL GetProcNo(PROCNO)
  IF(PROCNO.EQ.1) THEN
     if (len_trim(filename)==0) then
        buffer = ' '
        ewrite(3, *) 'Project name: '
        READ(*, '(A4096)', IOSTAT=IERR) buffer
        IF(LEN_TRIM(buffer).GT.lenth) THEN
           ewrite(3,*) "Filename too long :("
           STOP
        END IF
        filename(1:LEN_TRIM(buffer)) = buffer(1:LEN_TRIM(buffer))
     end if
  END IF
#ifdef HAVE_MPI
  CALL MPI_BCAST(filename, lenth, MPI_CHARACTER, 0, MPI_COMM_WORLD, IERR)
#endif

  RETURN
END SUBROUTINE GET_PROJECT_NAME

SUBROUTINE GetProcNo(ProcNo)
  IMPLICIT NONE
  INTEGER, INTENT(OUT)::ProcNo
#ifdef HAVE_MPI
  INCLUDE 'mpif.h'
  INTEGER IERR
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, PROCNO, IERR)
  PROCNO = PROCNO + 1
#else
  PROCNO = 1
#endif
  RETURN
END SUBROUTINE GetProcNo



SUBROUTINE OBS_SNAPMATRIX_SMEAN_SVD(istate,mxpoi,ntimemax,nsvd,nsvd_total,        &
           nsnap_obs,SnapNDT_obs,remesh,IF_dualwei,                        &
           iflagobs,obsv,smean_obs,                                                &
           leftsvd_obs,svdval_obs,coef_multivar_obs)

INTEGER istate,mxpoi,ntimemax,nsnap_obs,SnapNDT_obs,nsvd,nsvd_total
INTEGER                ::REMESH
LOGICAL                ::IF_dualwei
REAL, INTENT(IN)        :: iflagobs(0:ntimemax,mxpoi)
REAL, INTENT(INOUT)        :: obsv(0:ntimemax,istate)
REAL, INTENT(INOUT)        :: smean_obs(istate)
REAL, INTENT(INOUT)        :: leftsvd_obs(istate,nsvd),svdval_obs(nsvd_total)
REAL, INTENT(INOUT)        :: coef_multivar_obs(0:ntimemax,nsvd_total)
REAL,DIMENSION(:,:),ALLOCATABLE   :: snapmatrix_obs


! LOCAL VARIABLES
INTEGER        :: II,JJ,KK,I,J,K

  !---------------------------------------------------------!
  ! Read observational data from files                        !
  !---------------------------------------------------------*!
  IF(remesh.eq.0) THEN
   write(3,*) 'read observation data'
      open(10, file ='observation.dat')
          do ii=0,ntimemax
           read(10, 2000) (obsv(ii,jj), jj=1,istate)
          enddo
      close(10) 
      open(10, file ='uu.dat')
          do ii=0,ntimemax
           write(10, 2000) (obsv(ii,jj), jj=0*mxpoi+1,mxpoi)
          enddo
      close(10) 
      open(10, file ='vv.dat')
          do ii=0,ntimemax
           write(10, 2000) (obsv(ii,jj), jj=1*mxpoi+1,2*mxpoi)
          enddo
      close(10) 
 
!      open(10, file ='pp.dat')
!          do ii=0,ntimemax
!           read(10, 2000) (obsv(ii,jj), jj=3*mxpoi+1,istate)
!          enddo
!      close(10) 
   write(3,*) 'after observation data'
   ENDIF  !IF(remesh.eq.0) THEN
        
        ALLOCATE(snapmatrix_obs(istate,nsnap_obs))
        snapmatrix_obs =0.0
        smean_obs =0.0
        write(3,*) 'before buildup observation matrix'
        IF(IF_dualwei) THEN
          call buildsnapmat_dualwei(snapmatrix_obs,smean_obs,1,ntimemax,    &
             obsv,iflagobs,1,SnapNDT_obs,dual_wei_space)
        ELSE
          call buildsnapmat(snapmatrix_obs,smean_obs,1,ntimemax,    &
             obsv,iflagobs,1,SnapNDT_obs)
        ENDIF

        open(10, file ='snapmatrix_obs.dat')
        do ii=1,istate
           WRITE(10, 2000) (snapmatrix_obs(ii,jj), jj=1,nsnap_obs)
        enddo
        close(10)

        open(10, file ='uumatrix.dat')
        do ii=1,mxpoi
           WRITE(10, 2000) (snapmatrix_obs(ii,jj), jj=1,nsnap_obs)
        enddo
        close(10)

        open(10, file ='vvmatrix.dat')
        do ii=1*mxpoi+1,2*mxpoi
           WRITE(10, 2000) (snapmatrix_obs(ii,jj), jj=1,nsnap_obs)
        enddo
        close(10)

        open(10, file ='ppmatrix.dat')
        do ii=3*mxpoi+1,istate
           WRITE(10, 2000) (snapmatrix_obs(ii,jj), jj=1,nsnap_obs)
        enddo
        close(10)
 
        open(10, file ='smean_obs.dat')
         write(10, 2000) (smean_obs(i), i=1,istate)
        close(10) 

        leftsvd_obs =0.0
        svdval_obs =0.0
        write(3,*) 'calculation of POD vectors for observational data U'
        call snapsvd(mxpoi,nsnap_obs,snapmatrix_obs(1:mxpoi,:),nsvd,nsvd,leftsvd_obs(1:mxpoi,:),svdval_obs(1:nsvd))
        write(3,*) 'calculation of POD vectors for observational data V'
        call snapsvd(mxpoi,nsnap_obs,snapmatrix_obs(1*mxpoi+1:2*mxpoi,:),nsvd,nsvd,leftsvd_obs(1*mxpoi+1:2*mxpoi,:),svdval_obs(nsvd+1:2*nsvd))
        write(3,*) 'calculation of POD vectors for observational data W'
        call snapsvd(mxpoi,nsnap_obs,snapmatrix_obs(2*mxpoi+1:3*mxpoi,:),nsvd,nsvd,leftsvd_obs(2*mxpoi+1:3*mxpoi,:),svdval_obs(2*nsvd+1:3*nsvd))
        call snapsvd(mxpoi,nsnap_obs,snapmatrix_obs(3*mxpoi+1:4*mxpoi,:),nsvd,nsvd,leftsvd_obs(3*mxpoi+1:4*mxpoi,:),svdval_obs(3*nsvd+1:4*nsvd))
  
        DEALLOCATE(snapmatrix_obs)

        ! cproject the  observational data onto the POD space

        DO II = 0,ntimemax
           call fulltopodproj_multivar(nsvd_total,nsvd,nsvd_u,nsvd_phi,nvar,mxpoi,obsv(II,1:istate)-smean_obs(1:istate), &
                leftsvd_obs,coef_multivar_obs(II,:))
        ENDDO

        open(1,file='coef_multivar_obs.dat')
          do II=1,nsvd
             write(1,*) 'uu-coef'
             write(1,*) (coef_multivar_obs(JJ,II),JJ=0,ntimemax)
          enddo
          do II=nsvd+1,2*nsvd
             write(1,*) 'vv-coef'
             write(1,*) (coef_multivar_obs(JJ,II),JJ=0,ntimemax)
          enddo
          do II=2*nsvd+1,3*nsvd
             write(1,*) 'ww-coef'
             write(1,*) (coef_multivar_obs(JJ,II),JJ=0,ntimemax)
          enddo
        close(1)
2000 format (20(1X,E24.16)) 

END SUBROUTINE OBS_SNAPMATRIX_SMEAN_SVD


SUBROUTINE ERRORMEASURE_ADAPTMESH(NSVD_TOTAL,NSVD,NSVD_U,NSVD_PHI,MXNODS,NONODS,XNONOD,FREDOP,NPRESS,&
     X0,Y0,Z0,                       &
     U0_ADJ,V0_ADJ,W0_ADJ,PHI0_ADJ,  &
     U0_OBS,V0_OBS,W0_OBS,PHI0_OBS,FINDRM0,CENTRM0,   &
     M,MLX,MLY,MLZ,  &
     N,NLX,NLY,NLZ,NDGLNO0,XONDGL0,PNDGLN0,WEIGHT,TOTELE0,TOTELE,NLOC,NGI,MLOC,       &
     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,                            &
     NTIME,DT,D3,DCYL,GEOBAL,SCFACTH0,COEF_ADJ,           &
     obsv,iflagobs,COEF,coef_multivar_obs,LEFTSVD,SMEAN,leftsvd_obs,smean_obs,SnapNDT_obs,IFWIND,  &
     PGNODS,MPGLOC,MPG,MPGLX,MPGLY,MPGLZ,PGNDGLNO,SMEAN_PG,leftsvd_PG_x,leftsvd_PG_y,  &
     NPGCOLM0,NPGCOLM,PGFINDRM0,PGCOLM0,PGCENTRM0,DF,IF_REMESH_OBS,IF_INTERPOL,remesh,     &
     STOTEL0,STOTEL,SNLOC,SNDGLN0,NUM,IF_REMESH_VORTICITY)
         use FLDebug
         use mesh_adaptivity
         use assnav_module
         use tr2d_module
         use sparse_tools
         use AllSorts
  IMPLICIT NONE
    EXTERNAL get_xpctel  
      include '../assemble/paramnew.h'
  INTEGER get_xpctel
  INTEGER, INTENT(IN)  ::NSVD_TOTAL,NSVD,NSVD_U,NSVD_PHI
  INTEGER, INTENT(INOUT)  ::TOTELE0,TOTELE,NONODS,XNONOD,NLOC,MLOC,NGI,FREDOP,NPRESS,MXNODS
  REAL, INTENT(IN)     ::DF
  REAL, INTENT(INOUT),DIMENSION(MXNODS) ::X0,Y0,Z0
  REAL, INTENT(INOUT),DIMENSION(0:NTIME,MXNODS) ::U0_ADJ,V0_ADJ,W0_ADJ,PHI0_ADJ
  REAL, INTENT(INOUT),DIMENSION(0:NTIME,MXNODS) ::U0_OBS,V0_OBS,W0_OBS,PHI0_OBS
  REAL, INTENT(IN),DIMENSION(MXNODS) :: MUPTXX,MUPTXY,MUPTXZ
  REAL, INTENT(IN),DIMENSION(MXNODS) :: MUPTYY,MUPTYZ,MUPTZZ
  REAL, INTENT(INOUT),DIMENSION(NVAR*MXNODS,NSVD) :: LEFTSVD
  REAL, INTENT(INOUT),DIMENSION(NVAR*MXNODS) ::SMEAN
  REAL, INTENT(INOUT),DIMENSION(0:ntime,NSVD_TOTAL) :: COEF
  REAL, INTENT(INOUT),DIMENSION(0:ntime,NSVD_TOTAL) :: COEF_ADJ
  REAL, INTENT(IN)    ::DT
  INTEGER, INTENT(IN) ::NTIME
  LOGICAL, INTENT(IN) ::D3,DCYL
  REAL, INTENT(IN),DIMENSION(NLOC,NGI)  :: N,NLX,NLY,NLZ
  REAL, INTENT(IN),DIMENSION(NLOC,NGI)  :: M,MLX,MLY,MLZ
  REAL, INTENT(IN),DIMENSION(NGI)  ::WEIGHT
  INTEGER, INTENT(INOUT),DIMENSION(2*TOTELE0*NLOC):: NDGLNO0,XONDGL0
  INTEGER, INTENT(INOUT),DIMENSION(2*TOTELE0*MLOC)::PNDGLN0
  INTEGER STOTEL0,STOTEL,SNLOC
  INTEGER, INTENT(INOUT),DIMENSION(2*STOTEL0*SNLOC):: SNDGLN0
 ! coriolos force
  INTEGER, INTENT(IN)  :: GEOBAL
  REAL, INTENT(IN)  :: SCFACTH0
!  INTEGER NOCVA
!  REAL G(NOCVA)

  INTEGER, INTENT(INOUT)  :: FINDRM0(MXNODS+1),CENTRM0(MXNODS)
  INTEGER :: FINDRM,CENTRM,COLM
  INTEGER :: FINDRM1,CENTRM1,COLM1,FINDRM2,CENTRM2,COLM2
  INTEGER NCOLM
  INTEGER ,DIMENSION(:),ALLOCATABLE   ::COLM0
  LOGICAL ::IF_REMESH_OBS,IF_INTERPOL,IF_REMESH_VORTICITY
  ! Wind stress
  !------------
  LOGICAL IFWIND

  ! observational
  integer SnapNDT_obs
  integer,INTENT(IN),DIMENSION(0:ntime,MXNODS) ::iflagobs
  REAL,INTENT(IN),DIMENSION(0:ntime,istate):: obsv
  REAL coef_multivar_obs(0:ntime,NSVD_TOTAL)
  REAL leftsvd_obs(ISTATE,nsvd)
  REAl soux_mean,souy_mean,souz_mean,soux1,souy1,souz1,soux2,souy2,souz2
  REAL smean_obs(istate)
  REAL,DIMENSION(:),ALLOCATABLE   :: Bmat_obs

  ! GEOGRAPHIC PRESSURE
  !--------------------

  INTEGER MPGLOC,PGNODS
  REAL, INTENT(INOUT),DIMENSION(MPGLOC,NGI)  :: MPG,MPGLX,MPGLY,MPGLZ
  INTEGER, INTENT(INOUT),DIMENSION(2*TOTELE0*MPGLOC):: PGNDGLNO
!  REAL, INTENT(INOUT),DIMENSION(MXPOI+2*TOTELE0+1,NSVD_PHI) :: LEFTSVD_PG
  REAL  ::SMEAN_PG(MXPOI+2*TOTELE0)
!  REAL, INTENT(INOUT),DIMENSION(0:ntime1,NSVD_PHI) :: COEF_PG
  REAL,DIMENSION(:,:),ALLOCATABLE   ::MPGX,MPGY,MPGZ
  REAL,DIMENSION(:),ALLOCATABLE     :: WEI_PGx,WEI_PGy,WEI_PGz
  REAL L1(NGI),L2(NGI),L3(NGI),L4(NGI),PGWEIG(NGI) 
  REAL A11(NGI),A12(NGI),A13(NGI)
  REAL A21(NGI),A22(NGI),A23(NGI)
  REAL A31(NGI),A32(NGI),A33(NGI) 
!  REAL,DIMENSION(:,:),ALLOCATABLE   ::LEFTSVD_PG_x,LEFTSVD_PG_y
  REAL   ::LEFTSVD_PG_x(MXPOI+2*TOTELE0,NSVD_PHI),LEFTSVD_PG_y(MXPOI+2*TOTELE0,NSVD_PHI)
  REAL,DIMENSION(:),ALLOCATABLE   ::PGMATRIX_x,PGMATRIX_y,PGMATRIX_mean
  REAL,DIMENSION(:),ALLOCATABLE            ::PGVEC_x,PGVEC_y,PGVEC_mean
!  integer, dimension(:), allocatable:: PGFINDRM, PGCENTRM, PGCOLM
       INTEGER NPGCOLM0,NPGCOLM
       INTEGER PGFINDRM0(MXNODS+2*TOTELE0+1),PGCOLM0(2*NPGCOLM0)
       INTEGER PGCENTRM0(MXNODS+2*TOTELE0) 
       INTEGER PGFINDRM,PGCOLM,PGCENTRM

  
  REAL,DIMENSION(:),ALLOCATABLE   :: R
  INTEGER,DIMENSION(:),ALLOCATABLE   :: IMEM

! This is for the max and min element size tensors...
  INTEGER OPHSAM,NHSAMP,MXHSAM
  REAL CRIUP2,MINCH,MAXCH,MESTP2,MESTP1
  REAL ADATOU,ADATOV,ADATOW,ADATOT
  REAL ADOPTU,ADOPTV,ADOPTW,ADOPTT
  REAL ADWEIU,ADWEIV,ADWEIW,ADWEIT
  INTEGER, PARAMETER:: MESHCO=-40
  INTEGER IDENT(20)
  INTEGER TT(20)

  REAL,DIMENSION(:),ALLOCATABLE   :: XHSAMP,YHSAMP,ZHSAMP
  REAL,DIMENSION(:),ALLOCATABLE   :: HMINXX,HMINXY,HMINXZ
  REAL,DIMENSION(:),ALLOCATABLE   :: HMINYY,HMINYZ,HMINZZ

  REAL,DIMENSION(:),ALLOCATABLE   :: HMAXXX,HMAXXY,HMAXXZ
  REAL,DIMENSION(:),ALLOCATABLE   :: HMAXYY,HMAXYZ,HMAXZZ
  INTEGER remesh

      REAL phi0(0:ntimemax,mxpoi)
      REAL u0(0:ntimemax,mxpoi)
      REAL v0(0:ntimemax,mxpoi)
      REAL w0(0:ntimemax,mxpoi)
      common /fwdtraj/ phi0,u0,v0,w0            

! local variables

INTEGER NIMEM,NR,IPT,RPT
INTEGER X,Y,Z,PYSFUN,SNDGLN,NDGLNO,PNDGLN,XONDGL,NTSOL,PYSFU2
INTEGER II, JJ, KK,ITIME,NOD,ID
INTEGER,PARAMETER ::NDIM=3
INTEGER,PARAMETER ::NDIM2=9
INTEGER NONODS1,TOTELE1,NFIELDS
INTEGER FIELDS1(11),FIELDS2(11)
REAL,DIMENSION(:),ALLOCATABLE   :: RR,RLOCAL,PG
REAL,DIMENSION(:),ALLOCATABLE   :: FERROR,GLOERR,FERROR1,FERROR2
REAL ::GLOERR0,DFtilde
INTEGER:: xpctel,FERROR_POINTER
  INTEGER  ::PARA,halo_tag, halo_tag_p,NNODP,NNODPP
  PARAMETER(PARA=0,halo_tag=1, halo_tag_p=2)
  INTEGER   ::IERROR
INTEGER     :: NCT,COLCT,FINDCT,MIDCMC
INTEGER     :: NCT_PG,COLCT_PG,FINDCT_PG,MIDCMC_PG
LOGICAL NAV
INTEGER MAT1,MAT2,MAT3,MAT4,MAT5 
INTEGER VMAT1,VMAT2,VMAT3,VMAT4
integer, dimension(:), pointer :: boundary_ids, region_ids, coplanar_ids

REAL TT0(MXNODS)

  CHARACTER*240 UNAM,name
  INTEGER NUM
  INTEGER IOS

!  read adaptive paprameters
!-----------------------------
    write(3,*) 'in goal-based error measure'

 ALLOCATE(FERROR(NONODS*3*3))
 ALLOCATE(GLOERR(NONODS))
 ALLOCATE(FERROR1(NONODS*3*3))
 ALLOCATE(FERROR2(NONODS*3*3))
  NAV=.true.
  GLOERR0=1.0
  NTSOL=1
  IDENT=0
  TT=0
  TT0=0.0
  allocate(boundary_ids(stotel))
  boundary_ids=0
  allocate(region_ids(totele))
  region_ids=1
  allocate(coplanar_ids(totele))
  coplanar_ids=0

  write(3,*) 'read from adapt_parameters.dat'

  open(1,file='adapt_parameters.dat')
!     read(1,*) NCOLM
     read(1,*) CRIUP2,MINCH,MAXCH,MESTP1,MESTP2
     read(1,*) OPHSAM,NHSAMP,MXHSAM
     read(1,*) ADATOU,ADATOV,ADATOW,ADATOT
     read(1,*) ADOPTU,ADOPTV,ADOPTW,ADOPTT
     read(1,*) ADWEIU,ADWEIV,ADWEIW,ADWEIT
!  fix it later!!
     read(1,*) IDENT(1)

ALLOCATE(XHSAMP(MXHSAM))
ALLOCATE(YHSAMP(MXHSAM))
ALLOCATE(ZHSAMP(MXHSAM))
ALLOCATE(HMINXX(MXHSAM))
ALLOCATE(HMINXY(MXHSAM))
ALLOCATE(HMINXZ(MXHSAM))
ALLOCATE(HMINYY(MXHSAM))
ALLOCATE(HMINYZ(MXHSAM))
ALLOCATE(HMINZZ(MXHSAM))
ALLOCATE(HMAXXX(MXHSAM))
ALLOCATE(HMAXXY(MXHSAM))
ALLOCATE(HMAXXZ(MXHSAM))
ALLOCATE(HMAXYY(MXHSAM))
ALLOCATE(HMAXYZ(MXHSAM))
ALLOCATE(HMAXZZ(MXHSAM))

XHSAMP=0.0
YHSAMP=0.0
ZHSAMP=0.0
HMINXX=0.0
HMINXY=0.0
HMINXZ=0.0
HMINYY=0.0
HMINYZ=0.0
HMINZZ=0.0
HMAXXX=0.0
HMAXXY=0.0
HMAXXZ=0.0
HMAXYY=0.0
HMAXYZ=0.0
HMAXZZ=0.0

     IF(OPHSAM.NE.0) THEN
        DO ii =1,NHSAMP
           READ(1,*) XHSAMP(II),YHSAMP(II),ZHSAMP(II)
           READ(1,*) HMINXX(II),HMINXY(II),HMINXZ(II),HMINYY(II),HMINYZ(II),HMINZZ(II)
           READ(1,*) HMAXXX(II),HMAXXY(II),HMAXXZ(II),HMAXYY(II),HMAXYZ(II),HMAXZZ(II)
        ENDDO
     ENDIF
  close(1)

  write(3,*) 'after read from adapt_parameters.dat'

!MESTP2 = 100000.0

!HMINXX=0.00002
!HMINYY=0.00002
!HMINZZ=0.00002
!HMAXXX=0.04
!HMAXYY=0.04
!HMAXZZ=0.04

!MESTP2 = 10000.0
!HMINXX=0.001
!HMINYY=0.001
!HMINZZ=0.001
!HMAXXX=0.04
!HMAXYY=0.04
!HMAXZZ=0.04

! allocate R, IMEM
!-----------------
  NIMEM = 200*NONODS
  IPT   = 1
  NR    = 200*NONODS
  RPT   = 1
  ALLOCATE( IMEM(NIMEM) )
  ALLOCATE( R(NR) )
  IMEM=0
  R=0.0

          CALL ALOMEM(FINDRM1, IPT, NONODS+1, NIMEM)
          CALL ALOMEM(CENTRM1, IPT, NONODS,   NIMEM)
          COLM1=IPT          

CALL POSINMC_legacy(TOTELE,NONODS,NLOC,                                &
            IMEM(COLM1),NCOLM,NIMEM-IPT+1,        &
            IMEM(FINDRM1), IMEM(CENTRM1),        &
            NONODS,NLOC,NDGLNO0(1:TOTELE*NLOC),                &
            NDGLNO0(1:TOTELE*NLOC) )

          FINDRM0=0
          CENTRM0=0

ALLOCATE(COLM0(2*NCOLM))
COLM0=0
   colm0(1:NCOLM)=IMEM(COLM1:COLM1+NCOLM-1)
   FINDRM0(1:NONODS+1) = IMEM(FINDRM1:FINDRM1+NONODS+1-1)
   CENTRM0(1:NONODS) = IMEM(CENTRM1:CENTRM1+NONODS-1)

   IPT=IPT+NCOLM

        CALL ALOMEM(FINDCT, IPT, FREDOP+1, NIMEM)
        CALL ALOMEM(MIDCMC, IPT, FREDOP,   NIMEM)
        COLCT=IPT
        CALL POSINMC_legacy(TOTELE,NONODS,NLOC,                               &
                                IMEM(COLCT),NCT,NIMEM-IPT+1,          &
                                IMEM(FINDCT), IMEM(MIDCMC),           &
                                FREDOP,MLOC,PNDGLN0(1:TOTELE*MLOC),   &
                                NDGLNO0(1:TOTELE*NLOC))
        IPT=IPT+NCT

        CALL ALOMEM(FINDCT_PG, IPT, PGNODS+1, NIMEM)
        CALL ALOMEM(MIDCMC_PG, IPT, PGNODS,   NIMEM)
        COLCT_PG=IPT
        CALL POSINMC_legacy(TOTELE,NONODS,NLOC,                               &
                                IMEM(COLCT_PG),NCT_PG,NIMEM-IPT+1,          &
                                IMEM(FINDCT_PG), IMEM(MIDCMC_PG),           &
                                PGNODS,MPGLOC,PGNDGLNO(1:TOTELE*MPGLOC),   &
                                NDGLNO0(1:TOTELE*NLOC))
        IPT=IPT+NCT


!********************************************************
! TASK1: Calcalation of goal-based error measure
!********************************************************
DFtilde=DF/Real(Nonods)
!fix later!!!
TT0=0.0
FERROR = 0.0
FERROR1 =0.0
FERROR2 =0.0

goto 6000

!Method 1a: run the forward model
!--------------------------------
  write(3,*) 'calculate the residual from the forward model'

CALL REDUCED_MODEL_DIFFCOEF1_PG1_errormeasure(NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXNODS,NONODS,XNONOD,fredop,  &
     X0(1:NONODS),Y0(1:NONODS),Z0(1:NONODS),&
     FERROR1,COEF,U0(0:ntime,1:NONODS),V0(0:ntime,1:NONODS),W0(0:ntime,1:NONODS),PHI0(0:ntime,1:NONODS),    &
     U0_ADJ(0:ntime,1:NONODS),V0_ADJ(0:ntime,1:NONODS),W0_ADJ(0:ntime,1:NONODS),PHI0_ADJ(0:ntime,1:NONODS), &
     NCOLM,COLM0(1:NCOLM),FINDRM0(1:NONODS+1),NCT,IMEM(COLCT), IMEM(FINDCT),NCT_PG,IMEM(COLCT_PG), IMEM(FINDCT_PG),  &
     M,MLX,MLY,MLZ,                                    &
     N,NLX,NLY,NLZ,NDGLNO0(1:TOTELE*NLOC),XONDGL0(1:TOTELE*NLOC),PNDGLN0(1:TOTELE*MLOC),WEIGHT,TOTELE,NLOC,NGI,MLOC,      &
     MUPTXX(1:NONODS),MUPTXY(1:NONODS),MUPTXZ(1:NONODS),MUPTYY(1:NONODS),MUPTYZ(1:NONODS),MUPTZZ(1:NONODS),             &
     NTIME,DT,D3,DCYL,GEOBAL,SCFACTH0,IFWIND, &
     PGNODS,MPGLOC,MPG,MPGLX,MPGLY,MPGLZ,PGNDGLNO(1:TOTELE*MPGLOC),SMEAN_PG(1:PGNODS), &
     leftsvd_PG_x(1:PGNODS,1:NSVD_PHI),leftsvd_PG_y(1:PGNODS,1:NSVD_PHI),  &
     NPGCOLM,PGFINDRM0(1:PGNODS+1),PGCOLM0(1:NPGCOLM),PGCENTRM0(1:PGNODS),1.0)

write(name,'(a,i0)') 'FERROR1.dat.',num
OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
  write(1,*) (FERROR1(NOD),NOD=1,NONODS*9)
close(1)

!Method 1b: run the adjoint model
!--------------------------------
  write(3,*) 'calculate the residual from the adjoint model'

open(1,file='U0_ADJ1.dat')
  write(1,*) (U0_ADJ(400,II),II=1,NONODS)
close(1)
CALL REDUCED_MODEL_DIFFCOEF1consistent_PG_ADJ_errormeasure(NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXNODS,NONODS,XNONOD,fredop,&
     X0(1:NONODS),Y0(1:NONODS),Z0(1:NONODS),&
     FERROR2,COEF_ADJ,U0_ADJ(0:ntime,1:NONODS),V0_ADJ(0:ntime,1:NONODS),W0_ADJ(0:ntime,1:NONODS),PHI0_ADJ(0:ntime,1:NONODS), &
     U0(0:ntime,1:NONODS),V0(0:ntime,1:NONODS),W0(0:ntime,1:NONODS),PHI0(0:ntime,1:NONODS),   &
     U0_OBS(0:ntime,1:NONODS),V0_OBS(0:ntime,1:NONODS),W0_OBS(0:ntime,1:NONODS),      &
     NCOLM,COLM0(1:NCOLM),FINDRM0(1:NONODS+1),NCT,IMEM(COLCT), IMEM(FINDCT),NCT_PG,IMEM(COLCT_PG), IMEM(FINDCT_PG),   &
     M,MLX,MLY,MLZ,  &
     N,NLX,NLY,NLZ,NDGLNO0(1:TOTELE*NLOC),XONDGL0(1:TOTELE*NLOC),PNDGLN0(1:TOTELE*MLOC),WEIGHT,TOTELE,NLOC,NGI,MLOC,     &
     MUPTXX(1:NONODS),MUPTXY(1:NONODS),MUPTXZ(1:NONODS),MUPTYY(1:NONODS),MUPTYZ(1:NONODS),MUPTZZ(1:NONODS),       &
     NTIME,DT,D3,DCYL,GEOBAL,SCFACTH0,           &
     obsv,iflagobs,istate,coef_multivar_obs,leftsvd_obs,smean_obs,SnapNDT_obs,IFWIND,  &
     PGNODS,MPGLOC,MPG,MPGLX,MPGLY,MPGLZ,PGNDGLNO(1:TOTELE*MPGLOC),SMEAN_PG(1:PGNODS), &
     leftsvd_PG_x(1:PGNODS,1:NSVD_PHI),leftsvd_PG_y(1:PGNODS,1:NSVD_PHI),  &
     NPGCOLM,PGFINDRM0(1:PGNODS+1),PGCOLM0(1:NPGCOLM),PGCENTRM0(1:PGNODS),1.0,IF_REMESH_OBS,IF_REMESH_VORTICITY)
write(name,'(a,i0)') 'FERROR2.dat.',num
OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
  write(1,*) (FERROR2(NOD),NOD=1,NONODS*9)
close(1)


!--------------------------------

!GLOERR(1:NONODS)=0.5*( GLOERR1(1:NONODS)+GLOERR2(1:NONODS) )

!********************************************************
! TASK2: Calcalation of Hessain matrix
!********************************************************


! call qmesh to calculate the hessian
!--------------------------------------

! this is done bu subroutines REDUCED_MODEL_DIFFCOEF1_PG1_errormeasure and 
! REDUCED_MODEL_DIFFCOEF1consistent_PG_ADJ_errormeasure

!CALL QMESH(.true.,CRIUP2,MINCH,MAXCH,                    &
!     ADATOU,ADATOV,ADATOW,ADATOT,&
!     ADOPTU,ADOPTV,ADOPTW,ADOPTT,&
!     ADWEIU,ADWEIV,ADWEIW,ADWEIT,MXNODS,&
!     PHI0(INT(0.5*NTIME),1:NONODS),U0(INT(0.5*NTIME),1:NONODS),V0(INT(0.5*NTIME),1:NONODS), &
!     W0(INT(0.5*NTIME),1:NONODS),1,TT0(1:NONODS),NTSOL,&
!     NDGLNO0(1:TOTELE*NLOC),PNDGLN0(1:TOTELE*MLOC),TOTELE,NLOC,MLOC,FREDOP,&
!     NAV,D3,DCYL,X0(1:NONODS),Y0(1:NONODS),Z0(1:NONODS),&
!     FINDRM0(1:NONODS+1),COLM0(1:NCOLM),NCOLM,NONODS,CENTRM0(1:NONODS),&
!     FERROR,3,GLOERR0,R,NR, &
!     IMEM,NIMEM,&
!     xpctel,MESTP2,&
!     NGI,WEIGHT,&
!     N,NLX,NLY,NLZ &
!     ,& ! New element size stuff...
!     OPHSAM,NHSAMP,MXHSAM,XHSAMP,YHSAMP,ZHSAMP,&
!     HMINXX,HMINXY,HMINXZ,HMINYY,HMINYZ,HMINZZ,&
!     HMAXXX,HMAXXY,HMAXXZ,HMAXYY,HMAXYZ,HMAXZZ &
!     , & ! parallel stuff...&
!     PARA,halo_tag,NNODP,1.0E10)

!********************************************************
! TASK3: Calcalation of the metric tensor
!********************************************************

 !     DO NOD=1,NONODS
 !           DO II=1,NDIM*NDIM
 !                 FERROR((NOD-1)*NDIM*NDIM+II)= FERROR((NOD-1)*NDIM*NDIM+II)  &
!                        /GLOERR(NOD)
!            END DO
!      END DO

!       DO NOD=1,NONODS*9
!          FERROR(NOD)=0.5*(FERROR1(NOD)+FERROR2(NOD))
!       ENDDO

6000 continue

write(name,'(a,i0)') 'FERROR1.dat.',num
OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
  read(1,*) (FERROR1(NOD),NOD=1,NONODS*9)
close(1)
write(name,'(a,i0)') 'FERROR2.dat.',num
OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
  read(1,*) (FERROR2(NOD),NOD=1,NONODS*9)
close(1)

FERROR1=FERROR1/DF
FERROR2=FERROR2/DF

  NNODP=NONODS
  NNODPP=NNODP

      IF(ADOPTU.EQ.1) THEN
!     Make the error a relative error
         DO ID=1,NDIM2 
            DO NOD=1,NONODS
               FERROR1((NOD-1)*NDIM2+ID)=    &
                   FERROR1((NOD-1)*NDIM2+ID)/MAX(ADATOU,ABS(U0(100,NOD))+ABS(V0(100,NOD))+ABS(W0(100,NOD)) )
               FERROR2((NOD-1)*NDIM2+ID)=    &
                   FERROR2((NOD-1)*NDIM2+ID)/MAX(ADATOU,ABS(U0_ADJ(100,NOD))+ABS(V0_ADJ(100,NOD))+ABS(W0_ADJ(100,NOD)) )
            END DO
         END DO
      ENDIF
!
      IF(ADOPTU.EQ.-1) THEN
!     Mulitpply the error by the magnitude for the field to give preferance to large field values
!     e.g.  around detectors for data assimilation and inversion  ---  MDP
         DO ID=1,NDIM2 
            DO NOD=1,NONODS
               FERROR1((NOD-1)*NDIM2+ID)=    &
                   FERROR1((NOD-1)*NDIM2+ID)/MIN(ADATOU,ABS(U0(100,NOD))+ABS(V0(100,NOD))+ABS(W0(100,NOD)) )
               FERROR2((NOD-1)*NDIM2+ID)=    &
                   FERROR2((NOD-1)*NDIM2+ID)/MIN(ADATOU,ABS(U0_ADJ(100,NOD))+ABS(V0_ADJ(100,NOD))+ABS(W0_ADJ(100,NOD)) )
            END DO
         END DO
      ENDIF
write(3,*) 'remesh',remesh

      DO NOD=1,NNODP

         IF(OPHSAM.EQ.0) THEN
            CALL FLFormMetric1(FERROR1((NOD-1)*9+1),3,1.0,MINCH, MAXCH, MESTP2)
            CALL FLFormMetric1(FERROR2((NOD-1)*9+1),3,1.0,MINCH, MAXCH, MESTP2)
         ELSE
           CALL FLFormMetric2(FERROR1((NOD-1)*9+1),3,1.0,   &
                X0(NOD), Y0(NOD), Z0(NOD),                 &
                OPHSAM,NHSAMP,XHSAMP,YHSAMP,ZHSAMP,        &         
                HMINXX,HMINXY,HMINXZ,HMINYY,HMINYZ,HMINZZ, &                
                HMAXXX,HMAXXY,HMAXXZ,HMAXYY,HMAXYZ,HMAXZZ, &                
                MESTP2)
            CALL FLFormMetric2(FERROR2((NOD-1)*9+1),3,1.0,   &
                X0(NOD), Y0(NOD), Z0(NOD),                 &
                OPHSAM,NHSAMP,XHSAMP,YHSAMP,ZHSAMP,        &         
                HMINXX,HMINXY,HMINXZ,HMINYY,HMINYZ,HMINZZ, &                
                HMAXXX,HMAXXY,HMAXXZ,HMAXYY,HMAXYZ,HMAXZZ, &                
                MESTP2)
         ENDIF 
      END DO

write(3,*) 'after FLflorm'

!       DO NOD=1,NONODS*9
!          FERROR(NOD)=0.5*(FERROR1(NOD)+FERROR2(NOD))
!       ENDDO

!     NB Find the field total METRIC by combining the Hessians...
            CALL ALOMEM(MAT1, RPT, NDIM2, NR)
            CALL ALOMEM(MAT2, RPT, NDIM2, NR)
            CALL ALOMEM(MAT3, RPT, NDIM2, NR)
            CALL ALOMEM(MAT4, RPT, NDIM2, NR)
            CALL ALOMEM(VMAT1, RPT, NDIM, NR)
            CALL ALOMEM(VMAT2, RPT, NDIM, NR)
            CALL ALOMEM(VMAT3, RPT, NDIM, NR)
            CALL ALOMEM(VMAT4, RPT, NDIM, NR)
                CALL JACIMP1(FERROR1,FERROR2,NNODP,NONODS,NDIM,  &
                    R(MAT1),R(MAT2),R(MAT3),R(MAT4),            &
                    R(VMAT1),R(VMAT2),R(VMAT3))

  FERROR= FERROR2
open(1,file='FERROR.dat')
  write(1,*) (FERROR(NOD),NOD=1,NONODS*9)
close(1)
!****************************************************************
! TASK4: Mesh adaptivity based on the goal based error metric
!****************************************************************
  write(3,*) 'find the optimised mesh-adaptive mesh'


  R=0.0
  IMEM=0
  IPT   = 1
  RPT   = 1
  CALL ALOMEM(X,  RPT, MXNODS,    NR)
  CALL ALOMEM(Y,  RPT, MXNODS,    NR)
  CALL ALOMEM(Z,  RPT, MXNODS,    NR)
!  CALL ALOMEM(FERROR_POINTER,  RPT, 9*MXNODS,    NR)
  CALL ALOMEM(TT(1),  RPT, MXNODS,    NR)
  R(X:X+NONODS-1)=X0(1:NONODS)
  R(Y:Y+NONODS-1)=Y0(1:NONODS)
  R(Z:Z+NONODS-1)=Z0(1:NONODS)
  call alomem(FERROR_POINTER, rpt, 9*nonods, nr ) 
! DO II=1,9
!  R(FERROR_POINTER+(II-1)*MXNODS:FERROR_POINTER+(II-1)*MXNODS+NONODS-1)=FERROR(1+(II-1)*NONODS:(II-1)*NONODS+NONODS)
! ENDDO
  R(FERROR_POINTER:FERROR_POINTER+9*nonods-1)=FERROR(1:9*nonods)
!FIX IT LATER
  R(TT(1)+MXNODS-1)=0.0

  call alomem(pysfun, rpt, 3*mxnods, nr) 
  R(PYSFUN:PYSFUN+NONODS-1)=X0(1:NONODS)
  R(PYSFUN+MXNODS:PYSFUN+MXNODS+NONODS-1)=Y0(1:NONODS)
  R(PYSFUN+2*MXNODS:PYSFUN+2*MXNODS+NONODS-1)=Z0(1:NONODS)
  PYSFU2=PYSFUN+(3+1)*MXNODS+1
  CALL ALOMEM(XONDGL,  IPT, 2*TOTELE0*NLOC,    NIMEM)
  CALL ALOMEM(NDGLNO,  IPT, 2*TOTELE0*NLOC,    NIMEM)
  CALL ALOMEM(PNDGLN,  IPT, 2*TOTELE0*MLOC,    NIMEM)
  CALL ALOMEM(SNDGLN,  IPT, 2*STOTEL0*SNLOC,    NIMEM)


  IMEM(XONDGL:XONDGL+TOTELE*NLOC-1)=XONDGL0(1:TOTELE*NLOC)
  IMEM(NDGLNO:NDGLNO+TOTELE*NLOC-1)=NDGLNO0(1:TOTELE*NLOC)
  IMEM(PNDGLN:PNDGLN+TOTELE*MLOC-1)=PNDGLN0(1:TOTELE*NLOC)
  IMEM(SNDGLN:SNDGLN+STOTEL*SNLOC-1)=SNDGLN0(1:STOTEL*SNLOC)
  NONODS1=NONODS
  TOTELE1=TOTELE
write(3,*) 'before adapt_maes'
          CALL adapt_mesh(                                        &
              XONDGL,NDGLNO,SNDGLN,PNDGLN,                        &
              NLOC,SNLOC,MLOC,                                        &
              MXNODS,XNONOD,NONODS,TOTELE,FREDOP,STOTEL,        &
              NNODP,NNODPP,                                        &
              MESTP1,FERROR_POINTER,PYSFUN,PYSFU2,                &
              X,Y,Z,                                                &
              boundary_ids, region_ids, coplanar_ids, &
              R,NR,RPT, IMEM,NIMEM,IPT)

  deallocate(boundary_ids)
  deallocate(region_ids)
  deallocate(coplanar_ids)

write(3,*) 'after adapt_mash'


!**************************************************************************
! TASK5: Interpolate all variables from the original mesh onto the new mesh
!**************************************************************************
IF(IF_INTERPOL) THEN
ALLOCATE(RR(11*NONODS))
RR=0.0
ALLOCATE(RLOCAL(11*NONODS1))
RLOCAL=0.0

DO ITIME=0,ntime
 
  IF(NONODS.NE.FREDOP) THEN
    STOP 34
  ENDIF

      ! set up the current mesh0 system---adjoint mesh
     !................................................
     ! Allocate the memory for RR
     ! RR(1:3*NONODS*NPHASE) for the U,V,W at Mesh1

     NFIELDS = 11
     RLOCAL(1:NONODS1)=U0_OBS(ITIME,1:NONODS1)
     RLOCAL(1*NONODS1+1:2*NONODS1)=V0_OBS(ITIME,1:NONODS1)
     RLOCAL(2*NONODS1+1:3*NONODS1)=W0_OBS(ITIME,1:NONODS1)
     RLOCAL(3*NONODS1+1:4*NONODS1)=PHI0_OBS(ITIME,1:NONODS1)
     RLOCAL(4*NONODS1+1:5*NONODS1)=U0_ADJ(ITIME,1:NONODS1)
     RLOCAL(5*NONODS1+1:6*NONODS1)=V0_ADJ(ITIME,1:NONODS1)
     RLOCAL(6*NONODS1+1:7*NONODS1)=W0_ADJ(ITIME,1:NONODS1)
     RLOCAL(7*NONODS1+1:8*NONODS1)=U0(ITIME,1:NONODS1)
     RLOCAL(8*NONODS1+1:9*NONODS1)=V0(ITIME,1:NONODS1)
     RLOCAL(9*NONODS1+1:10*NONODS1)=W0(ITIME,1:NONODS1)
     RLOCAL(10*NONODS1+1:11*NONODS1)=PHI0(ITIME,1:NONODS1)

     ! the base pointers on mesh for the ajoint model..
     FIELDS2(1) = 1
     FIELDS2(2) = 1*NONODS+1
     FIELDS2(3) = 2*NONODS+1
     FIELDS2(4) = 3*NONODS+1
     FIELDS2(5) = 4*NONODS+1
     FIELDS2(6) = 5*NONODS+1
     FIELDS2(7) = 6*NONODS+1
     FIELDS2(8) = 7*NONODS+1
     FIELDS2(9) = 8*NONODS+1
     FIELDS2(10) = 9*NONODS+1
     FIELDS2(11) = 10*NONODS+1


     ! the base pointers on mesh1 for the forward model..
     !.................................................... 

     FIELDS1(1) = 1
     FIELDS1(2) = 1*NONODS1+1
     FIELDS1(3) = 2*NONODS1+1
     FIELDS1(4) = 3*NONODS1+1
     FIELDS1(5) = 4*NONODS1+1
     FIELDS1(6) = 5*NONODS1+1
     FIELDS1(7) = 6*NONODS1+1
     FIELDS1(8) = 7*NONODS1+1
     FIELDS1(9) = 8*NONODS1+1
     FIELDS1(10) = 9*NONODS1+1
     FIELDS1(11) = 10*NONODS1+1

     ! interpolete NU1,NV1,NW1 into the mesh of the adjoint model
     !............................................................

     CALL FLTetra4toTetra4(NONODS1,TOTELE1,X0,Y0,Z0,NDGLNO0,RLOCAL,FIELDS1,       &
          NFIELDS,NONODS,TOTELE,R(X),R(Y),R(Z),IMEM(NDGLNO),RR,FIELDS2,IERROR)
!     CALL FLTetra4toTetra4(NONODS1,TOTELE1,X0,Y0,Z0,NDGLNO0,RLOCAL,FIELDS1,       &
!          1,NONODS1,TOTELE1,X0,Y0,Z0,NDGLNO0,RR,FIELDS1,IERROR)
          

     U0_OBS(ITIME,1:MXNODS)=0.0
     V0_OBS(ITIME,1:MXNODS)=0.0
     W0_OBS(ITIME,1:MXNODS)=0.0
     PHI0_OBS(ITIME,1:MXNODS)=0.0
     U0_ADJ(ITIME,1:MXNODS)=0.0
     V0_ADJ(ITIME,1:MXNODS)=0.0
     W0_ADJ(ITIME,1:MXNODS)=0.0
     U0(ITIME,1:MXNODS)=0.0
     V0(ITIME,1:MXNODS)=0.0
     W0(ITIME,1:MXNODS)=0.0          
     PHI0(ITIME,1:MXNODS)=0.0          
     U0_OBS(ITIME,1:NONODS)=RR(1:NONODS)
     V0_OBS(ITIME,1:NONODS)=RR(1*NONODS+1:2*NONODS)
     W0_OBS(ITIME,1:NONODS)=RR(2*NONODS+1:3*NONODS)
     PHI0_OBS(ITIME,1:NONODS)=RR(3*NONODS+1:4*NONODS)

     U0_ADJ(ITIME,1:NONODS)=RR(4*NONODS+1:5*NONODS)
     V0_ADJ(ITIME,1:NONODS)=RR(5*NONODS+1:6*NONODS)
     W0_ADJ(ITIME,1:NONODS)=RR(6*NONODS+1:7*NONODS)

     U0(ITIME,1:NONODS)=RR(7*NONODS+1:8*NONODS)
     V0(ITIME,1:NONODS)=RR(8*NONODS+1:9*NONODS)
     W0(ITIME,1:NONODS)=RR(9*NONODS+1:10*NONODS)
     PHI0(ITIME,1:NONODS)=RR(10*NONODS+1:11*NONODS)

     if (20*int(ITIME/20) .eq. ITIME) then
        Call VTKOUTPUT1P(TOTELE,NONODS,1,NLoc,                               &
             D3,                                                             &
             IMEM(NDGLNO),                                                   &
             R(X),R(Y),R(Z),                                                 &
             U0_OBS(ITIME,:),V0_OBS(ITIME,:),W0_OBS(ITIME,:),PHI0_OBS(ITIME,:),U0_OBS(ITIME,:),  &
             1,                                                              &
             ITIME, 'FwdPODR')
     endif
write(3,*) 'TOTELE,NONODS,NLoc',TOTELE,NONODS,NLoc

ENDDO


DEALLOCATE(RR)
DEALLOCATE(RLOCAL)

        write(name,'(a,i0)') 'observation-newmesh.dat.',num
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
        do ii=0,ntimemax
           write(1, 2000) (U0_OBS(ii,jj), jj=1,nonods)
           write(1, 2000) (V0_OBS(ii,jj), jj=1,nonods)
           write(1, 2000) (W0_OBS(ii,jj), jj=1,nonods)
           write(1, 2000) (PHI0_OBS(ii,jj), jj=1,nonods)
        enddo
        close(1) 
!        open(1, file ='uvw0-newmesh.dat')
!        do ii=0,ntimemax
!           write(1, 2000) (U0(ii,jj), jj=1,nonods)
!           write(1, 2000) (V0(ii,jj), jj=1,nonods)
!           write(1, 2000) (W0(ii,jj), jj=1,nonods)
!           write(1, 2000) (PHI0(ii,jj), jj=1,nonods)
!        enddo
!        close(1) 

ELSE
        write(name,'(a,i0)') 'observation-newmesh.dat.',num
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
        do ii=0,ntimemax
           read(1, 2000) (U0_OBS(ii,jj), jj=1,nonods)
           read(1, 2000) (V0_OBS(ii,jj), jj=1,nonods)
           read(1, 2000) (W0_OBS(ii,jj), jj=1,nonods)
           read(1, 2000) (PHI0_OBS(ii,jj), jj=1,nonods)
        enddo
        close(1) 
!        open(1, file ='uvw0-newmesh.dat')
!        do ii=0,ntimemax
!           read(1, 2000) (U0(ii,jj), jj=1,nonods)
!           read(1, 2000) (V0(ii,jj), jj=1,nonods)
!           read(1, 2000) (W0(ii,jj), jj=1,nonods)
!           read(1, 2000) (PHI0(ii,jj), jj=1,nonods)
!        enddo
!        close(1) 
ENDIF


X0=0.0
Y0=0.0
Z0=0.0
XONDGL0=0
NDGLNO0=0
PNDGLN0=0
X0(1:NONODS)=R(X:X+NONODS-1)
Y0(1:NONODS)=R(Y:Y+NONODS-1)
Z0(1:NONODS)=R(Z:Z+NONODS-1)
XONDGL0(1:TOTELE*NLOC)=IMEM(XONDGL:XONDGL+TOTELE*NLOC-1)
NDGLNO0(1:TOTELE*NLOC)=IMEM(NDGLNO:NDGLNO+TOTELE*NLOC-1)
PNDGLN0(1:TOTELE*MLOC)=IMEM(PNDGLN:PNDGLN+TOTELE*MLOC-1)
SNDGLN0(1:STOTEL*SNLOC)=IMEM(SNDGLN:SNDGLN+STOTEL*SNLOC-1)

!****************************************************************
! TASK5: PLOT the obervation on the new mesh
!****************************************************************


DO ITIME=0,ntime
if (20*int(ITIME/20) .eq. ITIME) then
!        Call VTKOUTPUT1P(TOTELE,NONODS,1,NLoc,                               &
!             D3,                                                             &
!             IMEM(NDGLNO:NDGLNO+TOTELE*NLOC-1),                                                         &
!             R(X:X+NONODS-1),R(Y:Y+NONODS-1),R(Z:Z+NONODS-1),                                                          &
!             RR(1:NONODS),RR(1*NONODS+1:2*NONODS),RR(2*NONODS+1:3*NONODS),RR(3*NONODS+1:4*NONODS),RR(1:NONODS),  &
!             1,                                                              &
!             ITIME, 'FwdObs')
        Call VTKOUTPUT1P(TOTELE,NONODS,1,NLoc,                               &
             D3,                                                             &
             NDGLNO0,                                                         &
             X0,Y0,Z0,                                                          &
             U0(ITIME,:),V0(ITIME,:),W0(ITIME,:),PHI0(ITIME,:),U0(ITIME,:),  &
             1,                                                              &
             ITIME, 'Fwdobs ')
endif
ENDDO

!****************************************************************
! TASK6: OUTPUT the new data
!****************************************************************
        write(name,'(a,i0)') 'xyz-newmesh.dat.',num
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
           write(1, 2000) (X0(jj), jj=1,nonods)
           write(1, 2000) (Y0(jj), jj=1,nonods)
           write(1, 2000) (Z0(jj), jj=1,nonods)
        close(1)
           

        write(name,'(a,i0)') 'shape-newmesh.dat.',num
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
           write(1,*) TOTELE,NONODS,NLoc
           write(1, *) (NDGLNO0(jj), jj=1,TOTELE*NLOC)
           write(1, *) (XONDGL0(jj), jj=1,TOTELE*NLOC)
           write(1, *) (PNDGLN0(jj), jj=1,TOTELE*MLOC)
        close(1) 

2000 format (20(1X,E24.16)) 

!****************************************************************
! TASK7: Calculate the colm,FINDRM, CENTRM
!****************************************************************
          CALL ALOMEM(FINDRM2, IPT, NONODS+1, NIMEM)
          CALL ALOMEM(CENTRM2, IPT, NONODS,   NIMEM)
          COLM2= IPT

CALL POSINMC_legacy(TOTELE,NONODS,NLOC,                                &
            IMEM(COLM2),NCOLM,NIMEM-IPT+1,        &
            IMEM(FINDRM2), IMEM(CENTRM2),        &
            NONODS,NLOC,NDGLNO0,                &
            NDGLNO0)
          FINDRM0=0
          CENTRM0=0
DEALLOCATE(COLM0)
ALLOCATE(COLM0(2*NCOLM))
COLM0=0
   colm0(1:NCOLM)=IMEM(COLM2:COLM2+NCOLM-1)
   FINDRM0(1:NONODS+1) = IMEM(FINDRM2:FINDRM2+NONODS+1-1)
   CENTRM0(1:NONODS) = IMEM(CENTRM2:CENTRM2+NONODS-1)

   IPT=IPT+NCOLM

!**************************************************************************
! TASK6: CALCULATION OF PG......
!**************************************************************************

      MPGLOC=10
! Get the quadrature positions and weights for TETS...
         CALL TRIQUA(L1, L2, L3, L4, PGWEIG, D3,NGI)
! Work out the shape functions and there derivatives...
         CALL SHATRI(L1, L2, L3, L4, PGWEIG, D3, &
                 MPGLOC,NGI,                     &     
                 MPG,MPGLX,MPGLY,MPGLZ) 
!
!          CALL ALOMEM(PGNDGLNO, IPT, TOTELE*MPGLOC, NIMEM)
!
! The mid side nodes and corner nodes define...
!recover it later
! find PGNODS and PGNDGLNO
        PGNDGLNO=0
        CALL MIDSNODS(FINDRM0,NONODS,PGNODS,        &
                     CENTRM0,COLM0,NCOLM,                &
                     NDGLNO0,PGNDGLNO,TOTELE,NLOC,MPGLOC) 
! Calculate the matrix sparcity...
          CALL ALOMEM(PGFINDRM, IPT, PGNODS+1, NIMEM)
          CALL ALOMEM(PGCENTRM, IPT, PGNODS,   NIMEM)
          PGCOLM=IPT
          CALL POSINMC_legacy(TOTELE,PGNODS,MLOC,                                &
                                IMEM(PGCOLM),NPGCOLM,NIMEM-IPT+1,        &
                                IMEM(PGFINDRM), IMEM(PGCENTRM),        &
                                PGNODS,MPGLOC,PGNDGLNO,                &
                                PGNDGLNO)
pgcolm0=0
PGFINDRM0=0
PGCENTRM0=0

   pgcolm0(1:NPGCOLM)=IMEM(PGCOLM:PGCOLM+NPGCOLM-1)
   PGFINDRM0(1:PGNODS+1) = IMEM(PGFINDRM:PGFINDRM+PGNODS+1-1)
   PGCENTRM0(1:PGNODS) = IMEM(PGCENTRM:PGCENTRM+PGNODS-1)
   IPT=IPT+NPGCOLM

allocate(pg(pgnods))
pg=0.0

stop 145
CALL IN_OUTPUT_PG(2,PG,PGNODS,TOTELE,MPGLOC,NGI,PGNDGLNO,MPG,MPGLX,MPGLY,MPGLZ,  &
    D3,num,  &
    NPGCOLM,PGFINDRM0,PGCOLM0,PGCENTRM0)
deallocate(pg)
CALL READ_SHAPE_XYZ_MUPTXX2(TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC,FREDOP,NPRESS,        & 
                        M,MLX,MLY,MLZ,N,NLX,NLY,NLZ,WEIGHT,FINDRM0,CENTRM0,                        &
                        D3,DCYL,DT,                                                &
                        GEOBAL,SCFACTH0,                &                        
                        MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,X0,Y0,Z0,NDGLNO0,XONDGL0,PNDGLN0,    &
                        stotel,snloc,sndgln0)


 DEALLOCATE(FERROR)
 DEALLOCATE(GLOERR)
 DEALLOCATE(FERROR1)
 DEALLOCATE(FERROR2)

DEALLOCATE(R)
DEALLOCATE(IMEM)
DEALLOCATE(XHSAMP)
DEALLOCATE(YHSAMP)
DEALLOCATE(ZHSAMP)
DEALLOCATE(HMINXX)
DEALLOCATE(HMINXY)
DEALLOCATE(HMINXZ)
DEALLOCATE(HMINYY)
DEALLOCATE(HMINYZ)
DEALLOCATE(HMINZZ)
DEALLOCATE(HMAXXX)
DEALLOCATE(HMAXXY)
DEALLOCATE(HMAXXZ)
DEALLOCATE(HMAXYY)
DEALLOCATE(HMAXYZ)
DEALLOCATE(HMAXZZ)
DEALLOCATE(COLM0)

END SUBROUTINE ERRORMEASURE_ADAPTMESH
