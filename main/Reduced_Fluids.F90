#include "fdebug.h"

module reduced_fluids_module

  use fldebug
  use fluids_module
  use pod_support
  use reduced_model
  use signals
  use snapsvd_module
  use spud

  implicit none

  private

  public :: reducedfluids

contains

SUBROUTINE Reducedfluids(filename, filenamelen)

  include 'paramnew.h'      

  !------------------
  ! Initial condition
  !------------------
  REAL varincr(istate) 
  REAL varin(istate),advarin(istate)
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
  !---------------------------
  ! obervational data
  !----------------------------
  integer iprojop(mxpoi), idata, ipertdata,iuseobs,iusemean
  integer iflagobs(0:ntimemax,mxpoi)
  REAL obsv(0:ntimemax,istate)
  !  REAL obsvpert(0:ntimemax,istate)
  REAL weights(istate)
  !  common /obsval/ obsv, iflagobs,weights, bweights,ibackground 
  common /background/ backterm   

  ! Wind stress
  !------------
  LOGICAL IFWIND
  !---------------------------------
  !     POD variables
  !---------------------------------  
  integer nsvd,nsnap,nsvd_total,nsvd_u,nsvd_phi
  parameter (nsvd = 35)
  parameter (nsvd_u = 35)
  parameter (nsvd_phi = 35)
  parameter (nsvd_total = 3*nsvd_u+nsvd_phi)
  double precision :: leftsvd(istate,nsvd), svdval(nvar*nsvd)
  REAL coef_multivar(0:ntimemax,nsvd_total), romvar(0:ntimemax,istate)
  REAL coef(0:ntimemax,nsvd)
  REAL romgrad(nsvd)
  REAL tcoef1(nsvd)
  REAL tcoef_multivar(nsvd_total)
  REAL smean(istate)
  double precision :: snapmatrix(istate,nrsnapshots) 
  integer nsvd1    


  !----------------------------------------------------------
  !********************FOR FLUIIDTY**************************
  !----------------------------------------------------------
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

  REAL,DIMENSION(:),ALLOCATABLE     ::    MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ
  REAL,DIMENSION(:),ALLOCATABLE     ::    WEIGHT,DETWEI
  INTEGER,DIMENSION(:),ALLOCATABLE  ::    NDGLNO,XONDGL,PNDGLN
  REAL,DIMENSION(:),ALLOCATABLE     ::    X,Y,Z,NU,NV,NW,UG,VG,WG
  REAL,DIMENSION(:,:),ALLOCATABLE   ::    N,NLX,NLY,NLZ,M,MLX,MLY,MLZ
  REAL,DIMENSION(:,:),ALLOCATABLE   ::    NX,NY,NZ
  INTEGER,DIMENSION(:),ALLOCATABLE::          FINDRM,CENTRM
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
  INTEGER NPGCOLM
  INTEGER, ALLOCATABLE, DIMENSION(:)::PGFINDRM,PGCOLM,PGCENTRM
  !       INTEGER PGFINDRM(PGNODS+1),PGCOLM(NPGCOLM)
  !       INTEGER PGCENTRM(PGNODS) 


  CHARACTER(40) FILENAME_PG,FILENAME_MPG,FILENAME1,FILENAME2
  INTEGER K3
  character(len = 240) :: filenametemp
  ! for coriolis force
  INTEGER GEOBAL
  REAL SCFACTH0

  LOGICAL D3,DCYL
  INTEGER NONODS,XNONOD,FREDOP,NPRESS
  INTEGER TOTELE,NLOC,NGI,MLOC
  INTEGER NOD,ELE,II,JJ,KK,NSTEP,I,J,K,IP
  REAL DT

  REAL ftest,ftest2
  real volele,volvol,rnn,VOLUME
  real ml(mxpoi)
  integer iloc,jloc,igl,jgl,itime

  ! For adjoint Model
  REAL, DIMENSION(0:ntimemax,NSVD_TOTAL) :: COEF_ADJ
  LOGICAL  ::ADJOINT
  INTEGER IOS
  LOGICAL IFWRITE 
! error analyssi
  REAL::error_u(mxpoi),error_v(mxpoi),error_w(mxpoi),error_phi(mxpoi)
 
#ifdef HAVE_LIBARPACK
!%%%%%  if(have_option("/model/fluids/pod")) then
    !****************************************************!
     !   INPUT and INITILIZE THE VARIABLES                !
     !****************************************************!

     EWRITE(3,*) 'mx,mxpoi,istate,ntimemax,nrsnapshots',mx,mxpoi,istate,ntimemax,nrsnapshots
     !        stop

     IFWIND=.true.
     ADJOINT=.false.
     IFPG=1
     ! Initise the u0,v0,w0,phi0....
     !---------------------------

     DO K=0,ntimemax
        DO I=1,mxpoi
           u0(k,i)=0.0
           v0(k,i)=0.0
           w0(k,i)=0.0
           phi0(k,i)=0.0
        ENDDO
     ENDDO

     leftsvd=0.0
     svdval=0.0

     !--------4D-Var setup ---------------------------------
     idata = 0                 ! 1 if realistic data
     ipertdata = 0             ! 1 if perturbed observations
     ibackground = 1           ! 1 if background included in cost

     !-------POD setup-------------
     iuseobs = 0               ! 0 if use snapshots of the flow
     iusemean = 1              ! 1 if use snapmatrix centered around mean
     !------------------------------------------------------------------------

     !******************************************!
     !               RUN THE FULL MODEL         !
     !******************************************!


     ! prepare for runnning the model
     !--------------------------------
     ! Establish signal handlers.

     call initialise_signals


     IF(.true.) THEN !true. i.e. to run the full model and get POD vectors

          FILEVE=filename


        !RUN FORWARD MODEL WITH EXACT INITIAL CONDITIONS
        !-----------------------------------------------
        ewrite(3,*) 'before calling FLUID'
        CALL FLUIDS(FILEVE)

        ewrite(3,*) 'out of fluids'
        ewrite(3,*) 'end fwd run'
        EWRITE(3,*)  'start,end,total cpu time',TIME1,TIME2,TIME2-TIME1
        !          stop
        !**********************************************************************


        call setobsv(ntimemax,u0,v0,w0,phi0,obsv,iflagobs,3)

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
        !  call writeout(ntime,nx,ny,u0,v0,phi0,     &
        !       'uref.dat    ','vref.dat    ','phiref.dat  ')        

        ewrite(3,*)  'after set up the observation'  

        !*************************************************************!
        !          setup observation with perturbation                !    
        !*************************************************************!

        !         do itime = 0,ntimemax
        !           do ip=1,istate
        !           obsvpert(itime,ip) = obsv(itime,ip)
        !           enddo
        !         enddo        


        !     call pertobsv(ntimemax,obsvpert,iflagobs)
        EWRITE(3,*)  'start,end,total cpu time3',TIME1,TIME2,TIME2-TIME1
        EWRITE(3,*) 'mx,mxpoi,istate,ntimemax,nrsnapshots',mx,mxpoi,istate,ntimemax,nrsnapshots

        !*************************************************************!
        !               snapmatrix construction                       !     
        !*************************************************************!
        DO IP=1,mxpoi        
           VARIN(IP)=U0(0,IP)
           VARIN(mxpoi+IP)=V0(0,IP)
           VARIN(2*mxpoi+IP)=W0(0,IP)
           VARIN(3*mxpoi+IP)=PHI0(0,IP)
        ENDDO


        call buildsnapmat(snapmatrix,smean,iuseobs,ntimemax,    &
             obsv,iflagobs,iusemean,SnapNdT)

        ewrite(3,*) 'after building the Matrix'
        open(10, file ='snapmatrix.dat')
        do ii=1,istate
           WRITE(10, 2000) (snapmatrix(ii,jj), jj=1,nrsnapshots)
        enddo
        close(10) 

        open(10, file ='varin.dat')
        write(10, 2000) (varin(i), i=1,istate)
        close(10) 

        open(10, file ='smean.dat')
        write(10, 2000) (smean(i), i=1,istate)
        close(10) 

        open(10, file ='observation.dat')
        do ii=0,ntimemax
           write(10, 2000) (obsv(ii,jj), jj=1,istate)
        enddo
        close(10) 

        !**********************************************************************!
        !              perform the SVD of snapmatrix                           !
        !                       using ARPACK                                   !
        !**********************************************************************!

        !        call snapsvd(istate,nrsnapshots,snapmatrix,nsvd_total,nsvd,leftsvd,svdval)

        ewrite(3,*) 'calculation of POD vectors for U'
        call snapsvd(mxpoi,nrsnapshots,snapmatrix(1:mxpoi,:),nsvd,nsvd,leftsvd(1:mxpoi,:),svdval(1:nsvd))
        ewrite(3,*) 'calculation of POD vectors for V'
        call snapsvd(mxpoi,nrsnapshots,snapmatrix(1*mxpoi+1:2*mxpoi,:),nsvd,nsvd,leftsvd(1*mxpoi+1:2*mxpoi,:),svdval(nsvd+1:2*nsvd))
        ewrite(3,*) 'calculation of POD vectors for W'
        call snapsvd(mxpoi,nrsnapshots,snapmatrix(2*mxpoi+1:3*mxpoi,:),nsvd,nsvd,leftsvd(2*mxpoi+1:3*mxpoi,:),svdval(2*nsvd+1:3*nsvd))
        ewrite(3,*) 'calculation of POD vectors for PHI'
        call snapsvd(mxpoi,nrsnapshots,snapmatrix(3*mxpoi+1:4*mxpoi,:),nsvd,nsvd,leftsvd(3*mxpoi+1:4*mxpoi,:),svdval(3*nsvd+1:4*nsvd))
        ewrite(3,*) 'after caclculating the POD vector'

        open(10, file ='leftsvincr.dat')
        !write leftsv and sv values to output file
        !-----------------------------------------
        do ii=1,istate
           write(10, 2000) (leftsvd(ii,jj), jj=1,nsvd)
        enddo
        close(10) 
        open(10, file ='svdvalincr.dat')
        do i=1,nvar*nsvd
           write(10,200) svdval(i)
        enddo
        close(10)
2000    format (20(1X,E24.16)) 
200     format(1X,E24.16)



     ELSE !.false. read POD vectors from the files
        !-------------------------------------------------------

        !*************************************************************!
        !                snapmatrix construction                      !      
        !*************************************************************!

        !------------------------------------------
        if(.true.) then
           open(10, file ='observation.dat')
           do ii=0,ntimemax
              read(10, 2000) (obsv(ii,jj), jj=1,istate)
           enddo
           close(10) 

           snapmatrix =0.0
           smean =0.0
           ewrite(3,*) 'before buildup matrix'
           call buildsnapmat(snapmatrix,smean,1,ntimemax,    &
                obsv,iflagobs,1,SnapNDT)

           open(10, file ='snapmatrix.dat')
           do ii=1,istate
              WRITE(10, 2000) (snapmatrix(ii,jj), jj=1,nrsnapshots)
           enddo
           close(10)

           open(10, file ='smean.dat')
           write(10, 2000) (smean(i), i=1,istate)
           close(10) 
        endif
        !---------------------------
        ewrite(3,*) '1'
        open(10, file ='snapmatrix.dat')
        do ii=1,istate
           READ(10, 2000) (snapmatrix(ii,jj), jj=1,nrsnapshots)
        enddo
        close(10) 
        ewrite(3,*) '2'
      !  open(10, file ='varin.dat')
      !  read(10, 2000) (varin(i), i=1,istate)
      !  close(10) 
        varin(:)=obsv(0,:)
        ewrite(3,*) '3'
        open(10, file ='smean.dat')
        read(10, 2000) (smean(i), i=1,istate)
        close(10) 
        ewrite(3,*) '4'

        !**********************************************************************!
        !              perform the SVD of snapmatrix                           !
        !                       using ARPACK                                   !
        !**********************************************************************!
        ewrite(3,*) 'calculation of POD vectors for U'
        call snapsvd(mxpoi,nrsnapshots,snapmatrix(1:mxpoi,:),nsvd,nsvd,leftsvd(1:mxpoi,:),svdval(1:nsvd))
        ewrite(3,*) 'calculation of POD vectors for V'
        call snapsvd(mxpoi,nrsnapshots,snapmatrix(1*mxpoi+1:2*mxpoi,:),nsvd,nsvd,leftsvd(1*mxpoi+1:2*mxpoi,:),svdval(nsvd+1:2*nsvd))
        ewrite(3,*) 'calculation of POD vectors for W'
        call snapsvd(mxpoi,nrsnapshots,snapmatrix(2*mxpoi+1:3*mxpoi,:),nsvd,nsvd,leftsvd(2*mxpoi+1:3*mxpoi,:),svdval(2*nsvd+1:3*nsvd))
        ewrite(3,*) 'calculation of POD vectors for PHI'
        call snapsvd(mxpoi,nrsnapshots,snapmatrix(3*mxpoi+1:4*mxpoi,:),nsvd,nsvd,leftsvd(3*mxpoi+1:4*mxpoi,:),svdval(3*nsvd+1:4*nsvd))
        ewrite(3,*) 'after caclculating the POD vector'

        ! for the geostrophic pressure
        !------------------------------
        IF(IFPG.EQ.1) THEN

           DO  itime = 0,0
              WRITE(FILENAME_MPG, 20) ITIME
20            FORMAT('MPG_data',I5.5)
              K3 = index( FILENAME_MPG, ' ' ) - 1
              FILENAME2=FILENAME_MPG(1:K3)
              !   ewrite(3,*)  FILENAME2
              OPEN(2,FILE = FILENAME2)
              !    READ(2,*) PGNODS,TOTELE,NGI,MPGLOC,D3
              READ(2,*) PGNODS,TOTELE,MPGLOC,D3
              CLOSE(2)
           ENDDO

           NGI=11
           ewrite(3,*) 'PGNODS,TOTELE,MPGLOC,D3',PGNODS,TOTELE,MPGLOC,D3
           ALLOCATE(MPG(MPGLOC,NGI))  
           ALLOCATE(MPGLX(MPGLOC,NGI))
           ALLOCATE(MPGLY(MPGLOC,NGI))
           ALLOCATE(MPGLZ(MPGLOC,NGI))
           ALLOCATE(PG(0:ntimemax,PGNODS))
           ALLOCATE(PG_POD(1,PGNODS))
           ALLOCATE(PGNDGLNO(TOTELE*MPGLOC))
           ALLOCATE(snapmatrix_PG(PGNODS,nrsnapshots))
           ALLOCATE(smean_PG(PGNODS))
           ALLOCATE(leftsvd_PG(PGNODS,nsvd))
           ALLOCATE(svdval_PG(nsvd))
           ALLOCATE(varin_PG(PGNODS))
           ALLOCATE(varincr_PG(PGNODS))
           ALLOCATE(COEF_PG(0:ntimemax,NSVD_PHI))
           ALLOCATE(TCOEF_PG(NSVD_PHI))
           MPG=0.0
           MPGLX=0.0
           MPGLY=0.0
           MPGLZ=0.0
           PG=0.0
           PG_POD=0.0
           PGNDGLNO=0
           snapmatrix_PG=0.0
           smean_PG=0.0
           leftsvd_PG=0.0
           svdval_PG=0.0
           varin_PG=0.0
           varincr_PG=0.0
           COEF_PG=0.0
           TCOEF_PG=0.0

              OPEN(1,FILE = 'PG_data00000')
              open(2,file='MPG_data00000')
              READ(2,*) PGNODS,TOTELE,MPGLOC,D3
              READ(2,*) NPGCOLM
              ALLOCATE(PGFINDRM(PGNODS+1))
              ALLOCATE(PGCOLM(NPGCOLM))
              ALLOCATE(PGCENTRM(PGNODS))
              PGFINDRM=0
              PGCOLM=0
              PGCENTRM=0
              read(2,*) (PGFINDRM(I),I=1,PGNODS+1)
              read(2,*) (PGCOLM(I),I=1,NPGCOLM)
              read(2,*) (PGCENTRM(I),I=1,PGNODS)

              READ(1,*) (PGNDGLNO(I),I=1,TOTELE*MPGLOC)
              close(1)
              close(2)

        ENDIF



        !write leftsv and sv values to output file
        !-----------------------------------------


        open(10, file ='leftsvincr.dat')
        do ii=1,istate
           write(10, 2000) (leftsvd(ii,jj), jj=1,nsvd)
        enddo
        close(10)

        open(10, file ='svdvalincr.dat')
        do i=1,nvar*nsvd
           write(10,200) svdval(i)
        enddo
        close(10)

        if(.false.) then
           ! test to be orthonormal!!
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
                    ewrite(3,*) ii,jj,ftest
                 endif
              enddo
              do j=1,mxpoi
                 ftest2=ftest2+ leftsvd(j,ii)*leftsvd(j,ii)
              enddo
              ewrite(3,*) 'ftest2',ftest2
           enddo
           stop
        endif


        !--------------- 
        ewrite(3,*) '111'
        open(10, file ='leftsvincr.dat')
        do i=1,istate
           READ(10, 2000) (leftsvd(i,j), j=1,nsvd)
        enddo
        close(10) 

        !       open(10, file ='svdvalincr.dat')
        !       do i=1,nvar*nsvd
        !          READ(10,200) svdval(i)
        !       enddo
        !       close(10)

        open(10, file ='varin.dat')
        READ(10, 2000) (varin(i), i=1,istate)
        close(10) 

        open(10, file ='smean.dat')
        READ(10, 2000) (smean(i), i=1,istate)
        close(10) 

     ENDIF ! end .true. or false.
     !-------------------------------------------------------

     !*****************************************************************!
     !            project initial conditions on the POD space          !
     !*****************************************************************!

     do i=1,istate
        varincr(i) = varin(i) - smean(i)
     enddo

     ewrite(3,*) 'before project initial conditions on the POD space'
     !  call fulltopodproj(istate,nsvd,varincr,leftsvd,tcoef1)
     call fulltopodproj_multivar(nsvd_total,nsvd,nsvd_u,nsvd_phi,nvar,mxpoi,varincr,leftsvd,tcoef_multivar)

     DO II=0,ntimemax
        DO JJ=1,NSVD_TOTAL
           IF(II.EQ.0) THEN
              COEF_multivar(II,JJ)=TCOEF_multivar(JJ)
           ELSE
              COEF_multivar(II,JJ)=0.0
           ENDIF
        ENDDO
     ENDDO

     ! for geostrophic pressure
     !-------------------------

     IF(IFPG.EQ.1) THEN
        do i=1,PGNODS
           varincr_PG(i) = varin_PG(i) - smean_PG(i)
        enddo
        open(10, file ='varincr_PG.dat')
        write(10, 2000) (varincr_PG(i), i=1,PGNODS)
        close(10) 

        ewrite(3,*) 'before project initial PG on the POD space'
        call fulltopodproj(PGNODS,nsvd_phi,varincr_PG,leftsvd_PG,tcoef_PG)
        DO II=0,ntimemax
           DO JJ=1,NSVD_PHI
              IF(II.EQ.0) THEN
                 COEF_PG(II,JJ)=TCOEF_PG(JJ)
              ELSE
                 COEF_PG(II,JJ)=0.0
              ENDIF
           ENDDO
        ENDDO

     ENDIF


     !**************************************************************!
     !                RUNNING the reduced model                     !
     !**************************************************************!

     !READ DATA------------------------------------------  

     OPEN(1,file='shape.dat')

     READ(1,*) TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC,FREDOP,NPRESS
     READ(1,*) D3,DCYL
     ewrite(3,*)  TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC,FREDOP,NPRESS
     ewrite(3,*)  D3,DCYL

     ALLOCATE(N(NLOC,NGI))
     ALLOCATE(NLX(NLOC,NGI))
     ALLOCATE(NLY(NLOC,NGI))
     ALLOCATE(NLZ(NLOC,NGI))
     ewrite(3,*) 'after locating N,NLX..'
     ALLOCATE(M(MLOC,NGI))
     ALLOCATE(MLX(MLOC,NGI))
     ALLOCATE(MLY(MLOC,NGI))
     ALLOCATE(MLZ(MLOC,NGI))
     ewrite(3,*) 'after locating M,MLX..'
     ALLOCATE(MUPTXX(NONODS))
     ALLOCATE(MUPTXY(NONODS))
     ALLOCATE(MUPTXZ(NONODS))
     ALLOCATE(MUPTYY(NONODS))
     ALLOCATE(MUPTYZ(NONODS))
     ALLOCATE(MUPTZZ(NONODS))
     ewrite(3,*) 'after locating the viscosities'
     ALLOCATE(WEIGHT(NGI))
     ALLOCATE(X(NONODS))
     ALLOCATE(Y(NONODS))
     ALLOCATE(Z(NONODS))
     ewrite(3,*) 'after locating the locations,x,y,z'
     II=TOTELE*NLOC
     ALLOCATE(NDGLNO(II))
     ALLOCATE(XONDGL(II))
     II=TOTELE*MLOC
     ALLOCATE(PNDGLN(II))

     ALLOCATE(FINDRM(NONODS+1))
     ALLOCATE(CENTRM(NONODS))

     ewrite(3,*)  'staring to read data from FWD'
     DO JJ=1,NGI
        READ(1,*) (M(II,JJ),II=1,MLOC)
        READ(1,*) (MLX(II,JJ),II=1,MLOC)
        READ(1,*) (MLY(II,JJ),II=1,MLOC)
        READ(1,*) (MLZ(II,JJ),II=1,MLOC)
     ENDDO

     ewrite(3,*) 'finish reading M,MLX...'

     DO JJ=1,NGI
        READ(1,*) (N(II,JJ),II=1,NLOC)
        READ(1,*) (NLX(II,JJ),II=1,NLOC)
        READ(1,*) (NLY(II,JJ),II=1,NLOC)
        READ(1,*) (NLZ(II,JJ),II=1,NLOC)
     ENDDO
     ewrite(3,*) 'finish reading N,NLX...'

     READ(1,*) (WEIGHT(JJ),JJ=1,NGI)
     ewrite(3,*) 'finish reading weight'
     READ(1,*) (FINDRM(JJ),JJ=1,NONODS+1)
     READ(1,*) (CENTRM(JJ),JJ=1,NONODS)

     CLOSE(1)

     OPEN(1,file='XYZ_MUPTXX.dat')
     READ(1,*) DT
     !READ(1,*) OPTOME,GEOBAL
     FLAbort("This is broken - ask Stephan")
     !READ(1,*) OMEGA,OMEGA1,OMEGA2,OMEGA3,OMEGA4,SCFACTH0
     ewrite(3,*)  'DT',DT

     READ(1,*) (MUPTXX(JJ),JJ=1,NONODS)
     READ(1,*) (MUPTXY(JJ),JJ=1,NONODS)
     READ(1,*) (MUPTXZ(JJ),JJ=1,NONODS)
     READ(1,*) (MUPTYY(JJ),JJ=1,NONODS)
     READ(1,*) (MUPTYZ(JJ),JJ=1,NONODS)
     READ(1,*) (MUPTZZ(JJ),JJ=1,NONODS)
     ewrite(3,*) 'finishing reading viscosity'
     READ(1,*) (X(JJ),JJ=1,NONODS)
     READ(1,*) (Y(JJ),JJ=1,NONODS)
     READ(1,*) (Z(JJ),JJ=1,NONODS)
     ewrite(3,*) 'finish reading locations, x,y,z'
     READ(1,*) (NDGLNO(JJ),JJ=1,TOTELE*NLOC)
     READ(1,*) (XONDGL(JJ),JJ=1,TOTELE*NLOC)
     READ(1,*) (PNDGLN(JJ),JJ=1,TOTELE*MLOC)

     CLOSE(1)




     !------------------------------------
     ! run the reduced model--------------
     !------------------------------------

     ewrite(3,*) 'before running the reduced model'
     if(IFPG.EQ.1) then
        ALLOCATE(leftsvd_PG_x(PGNODS,nsvd_phi))
        ALLOCATE(leftsvd_PG_y(PGNODS,nsvd_phi))
        leftsvd_PG_x=0.0
        leftsvd_PG_y=0.0
        CALL REDUCED_MODEL_DIFFCOEF1_PG1(NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,NONODS,XNONOD,X,Y,Z,                          &
             LEFTSVD,SMEAN,COEF_multivar,M,MLX,MLY,MLZ,                                    &
             N,NLX,NLY,NLZ,NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,       &
             MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,                            &
             NTIMEMAX,DT,D3,DCYL,GEOBAL,SCFACTH0,IFWIND, &
             PGNODS,MPGLOC,MPG,MPGLX,MPGLY,MPGLZ,PGNDGLNO,SMEAN_PG,leftsvd_PG_x,leftsvd_PG_y,  &
             NPGCOLM,PGFINDRM,PGCOLM,PGCENTRM,.true.,1,1)
        !      DEALLOCATE(leftsvd_PG_x)
        !      DEALLOCATE(leftsvd_PG_y)
     else

        CALL REDUCED_MODEL_DIFFCOEF1(NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,NONODS,XNONOD,X,Y,Z,                       &
             LEFTSVD,SMEAN,COEF_multivar,M,MLX,MLY,MLZ,                                    &
             N,NLX,NLY,NLZ,NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,       &
             MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,                            &
             NTIMEMAX,DT,D3,DCYL,GEOBAL,SCFACTH0,IFWIND)
     endif


     ewrite(3,*) 'quit running the reduced model'
     !**************************************************************

     !****************project back from rom to full space
     DO NSTEP = 0,ntimemax
        tcoef_multivar(:)=coef_multivar(NSTEP,:)
        !     ewrite(3,*) 'NSTEP =',NSTEP
        !     call podtofullproj(istate,nsvd,varincr,leftsvd,tcoef1)
        call podtofullproj_multivar(nsvd_total,nsvd,nsvd_u,nsvd_phi,nvar,mxpoi,varincr,leftsvd,tcoef_multivar)

        DO IP=1,mxpoi        
           U0(NSTEP,IP)=varincr(IP)+ smean(IP)
           V0(NSTEP,IP)=varincr(mxpoi+IP)+smean(mxpoi+IP)
           W0(NSTEP,IP)=varincr(2*mxpoi+IP)+smean(2*mxpoi+IP)
           PHI0(NSTEP,IP)=varincr(3*mxpoi+IP)+smean(3*mxpoi+IP)
        ENDDO
        ! for geographic pressure
        !-------------------------

        ewrite(3,*) 'calculating the geostrophic pressure'
        IF(IFPG.EQ.1) THEN
           if(.false.) then
              tcoef_PG(:)=coef_PG(NSTEP,:)
              call podtofullproj(PGNODS,nsvd_phi,varincr_PG,leftsvd_PG,tcoef_PG)
              do ip=1,PGNODS
                 PG_POD(1,IP)=varincr_PG(IP)+smean_PG(IP)
              enddo

           else
              PG_POD(1,:)=SMEAN_PG(:)
              DO JJ=1,NSVD_PHI
                 DO IP=1,PGNODS
                    PG_POD(1,IP)=PG_POD(1,IP)+tcoef_multivar(JJ)*leftsvd_PG_y(IP,JJ)+tcoef_multivar(NSVD+JJ)*leftsvd_PG_x(IP,JJ)
                 ENDDO
              ENDDO
           endif
        ENDIF

        ! OUTPUT THE DATA FOR PLOTTER
        !----------------------------
        ewrite(3,*) 'plot the velocity field'
        !if(NSTEP.LE.40) then
        if (20*int(nstep/20) .eq. nstep) then

           if(40*int(nstep/40) .eq. nstep) then
              if(nstep.eq.0) then
                 open(1,file='PG_time.dat')
                 open(2,file='PG_POD_time.dat')
              endif
              write(1,*) 'nstep=',nstep
              write(2,*) 'nstep=',nstep
              write(1,*) (PG(nstep,IP),IP=1,PGNODS)
              write(2,*) (PG_POD(1,IP),IP=1,PGNODS)
              if(nstep.eq.400) then
                 CLOSE(1)
                 CLOSE(2)
              endif
           endif
     
     filenametemp="Fwdpod"  
           Call VTKOUTPUT1P(TOTELE,NONODS,1,NLoc,                               &
                D3,                                                             &
                NDGLNO,                                                         &
                X,Y,Z,                                                          &
                U0(NSTEP,:),V0(NSTEP,:),W0(NSTEP,:),PHI0(NSTEP,:),U0(NSTEP,:),  &
                1,                                                              &
                Nstep,filenametemp )
        error_u(:)=U0(NSTEP,:)-OBSV(NSTEP,1:mxpoi)
        error_v(:)=V0(NSTEP,:)-OBSV(NSTEP,mxpoi+1:2*mxpoi)
        error_w(:)=W0(NSTEP,:)-OBSV(NSTEP,2*mxpoi+1:3*mxpoi)
        error_phi(:)=PHI0(NSTEP,:)-OBSV(NSTEP,3*mxpoi+1:4*mxpoi)
  filenametemp="Errorpod"
        Call VTKOUTPUT1P(TOTELE,NONODS,1,NLoc,                               &
             D3,                                                             &
             NDGLNO,                                                         &
             X,Y,Z,                                                          &
             error_u(:),error_v(:),error_w(:),error_phi(:),error_u(:),       &
             1,                                                              &
             Nstep, filenametemp)
   endif

        !              CALL OUTPUT_POD(NSTEP,NONODS,XNONOD,TOTELE,NLOC,X,Y,Z,NDGLNO)
     ENDDO


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
     DEALLOCATE(FINDRM)
     DEALLOCATE(CENTRM)

     IF(IFPG.EQ.1) THEN
        DEALLOCATE(MPG)  
        DEALLOCATE(MPGLX)
        DEALLOCATE(MPGLY)
        DEALLOCATE(MPGLZ)
        DEALLOCATE(PG)
        DEALLOCATE(PG_POD)
        DEALLOCATE(PGNDGLNO)
        DEALLOCATE(snapmatrix_PG)
        DEALLOCATE(smean_PG)
        DEALLOCATE(leftsvd_PG)
        DEALLOCATE(svdval_PG)
        DEALLOCATE(varin_PG)
        DEALLOCATE(varincr_PG)

        if(.true.) then
           DEALLOCATE(PGFINDRM)
           DEALLOCATE(PGCOLM)
           DEALLOCATE(PGCENTRM)
           DEALLOCATE(leftsvd_PG_x)
           DEALLOCATE(leftsvd_PG_y)
        endif

     ENDIF

!%%%%  end if
#else
  FLAbort("No ARPACK support - cannot run reduced model.")
#endif
  RETURN
END SUBROUTINE Reducedfluids


SUBROUTINE ERROR_ANALYSIS(istate,ntimemax,mxpoi,nrsnapshots,TOTELE0,TOTELE,NONODS,NLoc,D3,NDGLNO,  &
              X,Y,Z,obsv,u0,v0,w0,phi0,dt,T0,NUM)

  INTEGER ::ntimemax,istate
  INTEGER :: mxpoi,TOTELE0,TOTELE,NONODS,NLoc
  INTEGER :: nrsnapshots
  REAL OBSV(0:ntimemax,istate)
  REAL X(mxpoi),Y(mxpoi),z(mxpoi)
!  REAL snapmatrix(8*mxpoi,nrsnapshots)
  REAL phi0(0:ntimemax,mxpoi)
  REAL u0(0:ntimemax,mxpoi)
  REAL v0(0:ntimemax,mxpoi)
  REAL w0(0:ntimemax,mxpoi)
  LOGICAL D3
  INTEGER NDGLNO(2*TOTELE0*NLoc)
  REAL SXX(0:ntimemax),SYY(0:ntimemax),SXY(0:ntimemax),CORRL(0:ntimemax),RRMS(0:ntimemax),RERROR,AERROR
  REAL RERROR1(0:ntimemax),AERROR1(0:ntimemax)
  REAL AVSXX(0:ntimemax),AVSYY(0:ntimemax)
  REAL SXX2(mxpoi),SYY2(mxpoi),SXY2(mxpoi),CORRL2(mxpoi),RRMS2(mxpoi)
  REAL RERROR2(mxpoi),AERROR2(mxpoi)
  REAL AVSXX2(mxpoi),AVSYY2(mxpoi)
  REAL T,TT,DT,t0
  REAL,PARAMETER:: PI=3.141592654
  REAL,DIMENSION(:),ALLOCATABLE ::num1,num2,f1,T1
  integer :: i,j,k,ii,jj,kk,ip
  CHARACTER*240 UNAM,name
  INTEGER NUM
  INTEGER IOS
  integer :: dump_no
  character(len = 240) :: filename

        AVSXX=0.0
        AVSYY=0.0
        AERROR1=0.0
        RERROR1=0.0
        AERROR2=0.0
        RERROR2=0.0
        AVSXX2=0.0
        AVSYY2=0.0
        
        DO k=0,ntimemax
          DO I = 1,nonods
                AVSXX(k)=AVSXX(k)+sqrt(obsv(k,i)**2+obsv(k,mxpoi+i)**2+obsv(k,2*mxpoi+i)**2)
                AVSYY(k)=AVSYY(k)+sqrt(u0(k,i)**2+v0(k,i)**2+w0(k,i)**2)
                AVSXX2(i)=AVSXX2(i)+sqrt(obsv(k,i)**2+obsv(k,mxpoi+i)**2+obsv(k,2*mxpoi+i)**2)
                AVSYY2(i)=AVSYY2(i)+sqrt(u0(k,i)**2+v0(k,i)**2+w0(k,i)**2)
           ENDDO
        ENDDO
           AVSXX=AVSXX/mxpoi
           AVSYY=AVSYY/mxpoi
           AVSXX2=AVSXX2/(ntimemax+1)
           AVSYY2=AVSYY2/(ntimemax+1)

        SXX=0.0
        SYY=0.0
        SXY=0.0
        CORRL=0.0
        RRMS=0.0
        CORRL2=0.0
        SXX2=0.0
        SYY2=0.0
        SXY2=0.0

        DO k=0,ntimemax

               DO I = 1,nonods
           SXX(k)= SXX(k)+(sqrt(obsv(k,i)**2+obsv(k,mxpoi+i)**2+obsv(k,2*mxpoi+i)**2)-AVSXX(k))**2
           SYY(k)= SYY(k)+(sqrt(u0(k,i)**2+v0(k,i)**2+w0(k,i)**2)-AVSYY(k))**2
           SXY(k)= SXY(k)+(sqrt(obsv(k,i)**2+obsv(k,mxpoi+i)**2+obsv(k,2*mxpoi+i)**2)-AVSXX(k))*   &
                          (sqrt(u0(k,i)**2+v0(k,i)**2+w0(k,i)**2)-AVSYY(k))
           SXX2(I)= SXX2(I)+(sqrt(obsv(k,i)**2+obsv(k,mxpoi+i)**2+obsv(k,2*mxpoi+i)**2)-AVSXX2(I))**2
           SYY2(I)= SYY2(I)+(sqrt(u0(k,i)**2+v0(k,i)**2+w0(k,i)**2)-AVSYY2(I))**2
           SXY2(I)= SXY2(I)+(sqrt(obsv(k,i)**2+obsv(k,mxpoi+i)**2+obsv(k,2*mxpoi+i)**2)-AVSXX2(I))*   &
                          (sqrt(u0(k,i)**2+v0(k,i)**2+w0(k,i)**2)-AVSYY2(I))

           if(ABS(sqrt(u0(k,i)**2+v0(k,i)**2+w0(k,i)**2)).le.1.0E-8) then
              RERROR=0.0
           else
              RERROR=(sqrt( (obsv(k,i)-u0(k,i))**2+(obsv(k,mxpoi+i)-v0(k,i))**2))&
                /sqrt(u0(k,i)**2+v0(k,i)**2)
!              AERROR=sqrt( (obsv(k,i)-u0(k,i))**2+(obsv(k,mxpoi+i)-v0(k,i))**2)
              AERROR=( (obsv(k,i)-u0(k,i))**2+(obsv(k,mxpoi+i)-v0(k,i))**2)
           endif
           RERROR1(k)=RERROR1(k)+RERROR
           AERROR1(k)=AERROR1(k)+AERROR
           RERROR2(I)=RERROR2(I)+RERROR
           AERROR2(I)=AERROR2(I)+AERROR
         ENDDO
          CORRL(k) = ABS(SXY(k))/SQRT(SXX(k)*SYY(k))
      ENDDO

          CORRL2(:) = ABS(SXY2(:))/SQRT(SXX2(:)*SYY2(:))

        RERROR1=RERROR1/nonods
!        AERROR1=AERROR1/nonods
        AERROR1=sqrt(AERROR1/nonods)
        RERROR2=RERROR2/(ntimemax+1)
!        AERROR2=AERROR2/(ntimemax+1)
        AERROR2=sqrt(AERROR2/(ntimemax+1))

           dump_no = 1
           filename = "ErrorAnalysis"
           Call VTKOUTPUT1P(TOTELE,NONODS,1,NLoc,                               &
                D3,                                                             &
                NDGLNO(TOTELE*NLOC),                                                         &
                X(1:nonods),Y(1:nonods),Z(1:nonods),                                                          &
                CORRL2(1:nonods),AERROR2(1:nonods),RERROR2(1:nonods),AERROR2(1:nonods),CORRL2(1:nonods),  &
                1,                                                              &
                dump_no, filename)

        write(name,'(a,i0)') 'CORRL_t.dat.',num
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
          do k=0,ntimemax
            write(1,*) k*DT+T0,CORRL(k)
          enddo
        close(1)
        
        write(name,'(a,i0)') 'RERROR1_t.dat.',num
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
          do k=0,ntimemax
            write(1,*) k*DT+T0,RERROR1(k)
          enddo
        close(1)

       write(name,'(a,i0)') 'AERROR1_t.dat.',num
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
          do k=0,ntimemax
            write(1,*) k*DT+T0,AERROR1(k)
          enddo
        close(1)

END SUBROUTINE ERROR_ANALYSIS

end module reduced_fluids_module
