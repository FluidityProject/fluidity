#include "fdebug.h"

module dynstall


    use fldebug
    use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, PYTHON_FUNC_LEN, pi
    implicit none
    
    real, parameter :: conrad = pi / 180.0 
    real, parameter :: condeg = 180.0 / pi  
    
    type LB_Type
    
    !=====================================
    ! These parameters follow the 
    ! formulation of turbinesFOAM
    logical :: StallFlag = .false.
    logical :: StallFlag_prev = .false. 
    real :: ds,time,time_prev
    real :: X,X_prev,Y,Y_prev,Z,Z_prev
    real :: etaL, etaL_prev
    real :: A1,A2,A3,b1,b2,K1,K2,r,r0,alphaDS0DiffDeg
    real :: D,D_prev,DF,DF_prev
    real :: CV,CV_prev,CNV,CNV_prev
    real :: fprime,fprime_prev,fDoublePrime
    real :: eta,S1,S2,CD0
    real :: Tf, Tv, Tvl
    real :: tau, tau_prev
    integer :: nNewTimes
    real :: fcrit
    real :: CMC
    real :: CNC, CNALpha
    real :: lambdaL,lambdaL_prev
    real :: lambdaM,lambdaM_prev
    real :: TI,E0
    real :: T1,T2
    real :: alphaEquiv
    real :: speed_of_sound
    real :: H,H_prev,J,J_prev
    real :: CNI,CMI,CNP,CN1
    real :: Tp,CN_prev,CNprime
    real :: DP,DP_prev
    real :: alphaPrime,alphaPrime_prev
    real :: DAlpha, Dalpha_prev, TAlpha
    real :: Alpha, Alpha_prev, alphaCrit
    real :: alpha1, alphaDS0, alphaSS, dalphaDS
    real :: Vx,CNF,CT,CN,cmFitExponent,K0,CM
    
    end type LB_Type


    contains

    subroutine dystl_init_LB(lb)
        
        implicit none
        type(LB_Type) :: lb
       

        lb%X=0.0
        lb%X_prev=0.0
        lb%Y=0.0
        lb%Y_prev=0.0
        lb%Z=0.0
        lb%Z_prev=0.0
        lb%A1=0.165
        lb%A2=0.335
        lb%b1=0.14
        lb%b2=0.53
        lb%speed_of_sound=343 ! [m/s]
        lb%A3=0.5
        lb%D=0.0
        lb%D_prev=0.0
        lb%fprime=0.0
        lb%fprime_prev=0.0
        lb%DF=0.0
        lb%DF_prev=0.0
        lb%CV=0.0
        lb%CV_prev=0.0
        lb%CNV=0.0
        lb%CNV_prev=0.0
        lb%eta=0.975
        lb%S1=1e12
        lb%S2=1e12
        lb%CD0=1e12
        lb%StallFlag = .false.
        lb%StallFlag_prev=.false.
        lb%Tp=1.7
        lb%Tf=3.0
        lb%tau=0.0
        lb%tau_prev=0.0
        lb%nNewTimes=0
        lb%fcrit=0.6
        lb%K1=0.0
        lb%K2=0.0
        lb%etaL=0.0
        lb%etaL_prev=0.0
        lb%A3=0.5
        lb%T1=20.0
        lb%T2=4.5
        lb%H=0.0
        lb%H_prev=0.0
        lb%lambdaL=0.0
        lb%lambdaL_prev=0.0
        lb%J=0.0
        lb%J_prev=0.0
        lb%lambdaM=0.0
        lb%lambdaM_prev=0.0
        lb%CMC=0.0
        lb%CMI=0.0
        lb%r=0.0
        lb%TAlpha=6.30
        lb%alphaPrime_prev=0.0
        lb%alphaCrit=17.0
        lb%r=0.0
        lb%r0=0.01
        lb%B1=0.5
        lb%B2=0.2
        lb%E0=0.15
        lb%alphaDS0DiffDeg=3.6
        lb%DAlpha=0.0
        lb%DAlpha_prev=0.0
        lb%Tv=11.0
        lb%Tvl=9.0
        lb%eta=0.975

    end subroutine dystl_init_LB
    
    subroutine calcAlphaEquiv(lb,mach,alpha,chord,Ur,adot)
        implicit none
        ! Calculates the equivalent angle of attack 
        ! after having applied the defeciency functions
        type(LB_Type),intent(inout) :: lb
        real, intent(inout) :: alpha
        real, intent(in) :: mach,chord,adot,ur
        real :: T3,beta

        T3= 1.25*mach
        beta=1-mach**2
        lb%X=lb%X_prev*exp(-beta*lb%ds/lb%T1)+lb%A1*(lb%etaL-lb%etaL_prev)*exp(-beta*lb%ds/(2.0*lb%T1))
        lb%Y=lb%Y_prev*exp(-beta*lb%ds/lb%T2)+lb%A2*(lb%etaL-lb%etaL_prev)*exp(-beta*lb%ds/(2.0*lb%T2))
        lb%Z=lb%Z_prev*exp(-beta*lb%ds/T3)+lb%A1*(lb%etaL-lb%etaL_prev)*exp(-beta*lb%ds/(2.0*T3))
        lb%etaL=alpha + chord/(2.0*Ur)*adot
        lb%alphaEquiv = lb%etaL-lb%X-lb%Y-lb%Z
        if (abs(lb%alphaEquiv)>2*pi) then
            lb%alphaEquiv=mod(lb%alphaEquiv,2*pi)
        endif

    end subroutine calcAlphaEquiv
    
    subroutine calcUnsteady(lb,alpha,mach,chord,Ur,adot,dt)
        implicit none
        type(LB_Type),intent(inout) :: lb
        real,intent(in) :: alpha,chord,Ur,adot,dt,mach

        ! Calculate the circulatory normal force coefficient
        lb%CNC=lb%CNAlpha*lb%alphaEquiv

        ! Calculate the impulsive normal force coefficient
        lb%lambdaL=pi/4.0*(alpha+chord/(4.0*Ur)*adot)
        lb%TI=chord/lb%speed_of_sound*(1.0+3.0*mach)/4.0
        lb%H=lb%H_prev*exp(-dt/lb%TI)+(lb%lambdaL-lb%lambdaL_prev)*exp(-dt/(2.0*lb%TI))
        lb%CNI=4.0/mach*lb%H

        ! Calculate the impulsive moment coefficient
        lb%lambdaM = 3*pi/16*(alpha + chord/(4.0*Ur)*adot)+pi/16*chord/Ur*adot
        lb%J = lb%J_prev*exp(-dt/lb%TI)+(lb%lambdaM-lb%lambdaM_prev)*exp(-dt/(2.0*lb%TI))
        lb%CMI=-4.0/mach*lb%J

        ! Calculate total normal force coefficient
        lb%CNP = lb%CNC + lb%CNI

        ! Apply first-order lag to normal force coefficient
        lb%DP=lb%DP_prev*exp(-lb%ds/lb%Tp)+(lb%CNP-lb%CN_prev)*exp(-lb%ds/(2.0*lb%Tp))
        lb%CNprime=lb%CNP-lb%DP

        ! Calculate lagged angle of attack
        lb%DAlpha=lb%Dalpha_prev*exp(-lb%ds/lb%TAlpha)+(alpha-lb%alpha_Prev)*exp(-lb%ds/(2.0*lb%TAlpha))
        lb%alphaPrime=alpha-lb%DAlpha

        ! Calculate reduced pitch rate
        lb%r=adot*chord/(2.0*ur)

        ! Calculate alphaDS0
        lb%dAlphaDS=lb%alphaDS0DiffDeg/180.0*pi
        lb%alphaDS0 = lb%alphaSS + lb%dalphaDS

        if(abs(lb%CNPrime)>lb%CN1) then
            lb%StallFlag=.true.
        endif
        

    end subroutine calcUnsteady

    subroutine calcSeparated(lb)
        implicit none
        type(LB_type),intent(inout):: lb
        real :: f, m, cmf, cmv

        if (abs(lb%alphaPrime) < lb%alpha1) then
            lb%fprime=1.0-0.4*exp((abs(lb%alphaPrime) -lb%alpha1)/lb%S1) 
        else
            lb%fprime=0.02+0.58*exp((lb%alpha1 -abs(lb%alphaPrime))/lb%S1)
        endif

        if(lb%StallFlag_prev.eqv..false.) then
            lb%tau=0.0
        else
            if(lb%tau.eq.lb%tau_prev) then
                lb%tau=lb%tau_prev+lb%ds
            endif
        endif

        ! Calculate dynamic seperation point
        lb%DF=lb%DF_prev*exp(lb%ds/lb%Tf)+(lb%fprime-lb%fprime_prev)*exp(lb%ds/(2.0*lb%Tf))
        lb%fDoublePrime=lb%fprime-lb%DF

        ! Calculate vortex modulation parameter
        if (lb%tau>=0.and.lb%tau<=lb%Tvl) then
            lb%Vx=sin(pi*lb%tau/(2.0*lb%Tvl))**1.5
        elseif (lb%tau>lb%Tvl) then
            lb%Vx=cos(pi*(lb%tau-lb%Tvl)/lb%Tv)**2
        endif
        
        if (abs(lb%alpha) < lb%alpha_prev) then
            lb%Vx=0.0
        endif

        ! Calculate normal force coefficient including dynamic separation point
        lb%CNF=lb%CNAlpha*lb%alphaEquiv*((1.0+sqrt(lb%fDoublePrime))/2.0)**2

        ! Calculate tangential force coefficient
        lb%CT=lb%eta*lb%CNAlpha*lb%alphaEquiv*lb%alphaEquiv*(sqrt(lb%fDoublePrime)-lb%E0)

        ! Calculate static trailing-edge seperation point
        if (abs(lb%alpha) < lb%alpha1) then
            f=1.0-0.4*exp((abs(lb%alpha)-lb%alpha1)/lb%S1)
        else
            f=0.02+0.58*exp((lb%alpha1-abs(lb%alpha))/lb%S2)
        endif

        ! Evaluate vortex lift contributions
        lb%CNV=lb%B1*(lb%fDoublePrime-f)*lb%Vx

        ! Total normal force coefficient is the combination of that from 
        ! circulatory effects, impulsive effects, dynamic separation, and vortex
        ! lift
        lb%CN=lb%CNF+lb%CNV

        ! Calculate moment coefficient
        m=lb%cmFitExponent
        cmf=(lb%K0+lb%K1*(1-lb%fDoublePrime)+lb%K2*sin(pi*lb%fDoublePrime**m))*lb%CNC
        ! + moment coefficient at Zero lift angle of attack
        cmv = lb%b2*(1.0-cos(pi*lb%tau/lb%Tvl))*lb%CNV
        lb%CM=cmf+cmv+lb%CMI

    end subroutine calcSeparated
    
    subroutine update_LBmodel(lb,time)
        implicit none
        type(LB_type), intent(inout) ::lb
        real,intent(in) :: time

        lb%time_prev=time
        lb%Alpha_prev=lb%Alpha
        
    end subroutine update_LBmodel

    subroutine LeishmanBeddoesCorrect(lb,time,dt,alpha)
        
        implicit none
        type(LB_Type), intent(inout) :: lb
        real, intent(in) :: time, dt, alpha
        
        ! update previous values if time has changed
        if (time.ne.lb%time_prev) then
            lb%nNewTimes=lb%nNewTimes+1
            if (lb%nNewTimes > 1) then
                call update_LBmodel(lb,time)
            end if
        end if


    end subroutine LeishmanBeddoesCorrect

    !subroutine LB_EvalIdealCL(AOA,AOA0,CLa,RefFlag,CLID)

    !! AOA inputs in radians
    !! AOA0 is zero lift AOA
    !! RefFlag defines whether to output reference CL or ideal CL
    !! CLa is reference lift slope (per radian) to be used for reference CL (ideal CLa is 2*pi)
    !implicit none
    !real :: AOA, AOA0, CLa
    !integer :: RefFlag
    !real :: CLID
    !real :: IDS, aID, d1, CLaI, aIDc

    !aID=AOA-AOA0
    !call Force180(aID)
    !! reflect function across axis for abs(aID)>90
    !IDS=1
    !if (aID>pi/2.0) then
    !    aID=pi-aID
    !    IDS=-1.0
    !else if (aID<-pi/2.0) then
    !    aID=-pi-aID
    !    IDS=-1.0
    !end if

    !! If RefFlag is 1, output reference CL, otherwise round off ideal CL at high AOA
    !if (RefFlag==1) then
    !    CLID=IDS*CLa*aID
    !else
    !    ! round off the ideal CL after cutoff AOA
    !    aIDc=30.0*conrad
    !    d1=1.8
    !    CLaI=2.0*pi
    !    if (abs(aID)<aIDc) then
    !        CLID=IDS*CLaI*aID
    !    else
    !        CLID=IDS*(CLaI*(aIDc-1.0/d1*sin(d1*aIDc))+CLaI/d1*sin(d1*aID))
    !    end if
    !end if
    !
    !end subroutine LB_EvalIdealCL

    !subroutine Force180(a)
    !
    !    implicit none

    !    real :: a
    !    ! alpha in radians
    !    if (a>pi) then
    !        a=a-2.0*pi
    !    elseif (a<-pi) then
    !        a=a+2.0*pi
    !    endif
    !end subroutine Force180
    !
    !SUBROUTINE LB_UpdateStates(lb,ds)

    !    implicit none

    !    type(LB_Type),intent(inout) :: lb
    !    real, intent(in) :: ds
    !    integer :: i, nei, j, IsBE
    !    
    !    real :: Tf,TfRef,Tp,TvRef, Tv, TvL

    !    ! Set model parameters. All of these are potentially a function of Mach
    !    ! and are set to low mach values...
    !    Tp=4.0                  ! time constant on LE pressure response to change in CL
    !    TfRef=3.0               ! time constant on TE separation point travel
    !    TvRef=6.3               ! time constant on LE vortex lift indicial function
    !    TvL=11.0                ! Characteristic LE vortex travel time



    !! Eval LE separation state
    !! Note: CLCrit is critical (ideal) CL value for LE separation. This is
    !! approximately equal to the CL that would exist at the angle of attack at
    !! max CL if the CL curve had remained linear
    !if (lb%LESepState==0 .AND. (lb%CLRefLE>lb%CLCritP .OR. lb%CLRefLE<lb%CLCritN)) then
    !    ! In LE separation state
    !    lb%LESepState=1
    !    lb%sLEv=0 ! reset leading edge vortex time counter

    !    ! Set logic state flags (for model diagnosis output)
    !    lb%LB_LogicOutputs(5)=1
    !else if (lb%LESepState==1 .AND. (lb%CLRefLE<lb%CLCritP .AND. lb%CLRefLE>lb%CLCritN)) then
    !    ! Out of LE separation state
    !    lb%LESepState=0
    !    lb%sLEv=0 ! reset leading edge vortex time counter

    !    ! Set logic state flags (for model diagnosis output)
    !    lb%LB_LogicOutputs(5)=2
    !end if

    !! Set time constants based on LE separation state and TE separation point
    !! location. Different time constants for abs(CL) increasing vs. decreasing
    !if (lb%LESepState==1) then
    !    if (lb%sLEv<TvL) then
    !        if (lb%CLRateFlag>0) then
    !            !Tf=3.0*TfRef ! original
    !            Tf=4.0*TfRef
    !            Tv=TvRef

    !            ! Set logic state flags (for model diagnosis output)
    !            lb%LB_LogicOutputs(8)=1
    !        else
    !            Tf=1.0/2.0*TfRef
    !            Tv=1.0/2.0*TvRef

    !            ! Set logic state flags (for model diagnosis output)
    !            lb%LB_LogicOutputs(8)=2
    !        end if

    !        ! Set logic state flags (for model diagnosis output)
    !        lb%LB_LogicOutputs(7)=1
    !    else if (lb%sLEv<2.0*TvL) then
    !        if (lb%CLRateFlag>0) then
    !            ! orig 
    !            !Tf=1.0/3.0*TfRef
    !            !Tv=1.0/4.0*TvRef
    !            Tf=2.0*TfRef
    !            Tv=TvRef

    !            ! Set logic state flags (for model diagnosis output)
    !            lb%LB_LogicOutputs(8)=3
    !        else
    !            Tf=1.0/2.0*TfRef
    !            Tv=1.0/2.0*TvRef

    !            ! Set logic state flags (for model diagnosis output)
    !            lb%LB_LogicOutputs(8)=4                                                        
    !        end if

    !        ! Set logic state flags (for model diagnosis output)
    !        lb%LB_LogicOutputs(7)=2                                                
    !    else
    !        ! orig
    !        !Tf=4.0*TfRef
    !        !Tv=0.9*TvRef
    !        Tf=TfRef
    !        Tv=TvRef

    !        ! Set logic state flags (for model diagnosis output)
    !        lb%LB_LogicOutputs(7)=3
    !    end if

    !    ! Set logic state flags (for model diagnosis output)
    !    lb%LB_LogicOutputs(6)=1
    !else
    !    if (lb%F>0.7) then
    !        Tf=TfRef

    !        ! Set logic state flags (for model diagnosis output)
    !        lb%LB_LogicOutputs(7)=4
    !    else
    !        Tf=2*TfRef

    !        ! Set logic state flags (for model diagnosis output)
    !        lb%LB_LogicOutputs(7)=5
    !    end if
    !    Tv=TvRef

    !    ! Set logic state flags (for model diagnosis output)
    !    lb%LB_LogicOutputs(6)=2
    !end if

    !! update LE vortex time counter if in LE separation state
    !if (lb%LESepState==1) then
    !    lb%sLEv=lb%sLEv+ds

    !    ! Set logic state flags (for model diagnosis output)
    !    lb%LB_LogicOutputs(9)=1
    !end if

    !! Update states, first order lag equations, exponential recursion form (midpoint rule version)
    ! lb%dp=lb%dp*exp(-ds/Tp)+(lb%CLRef-lb%CLRef_Last)*exp(-ds/(2*Tp))
    ! lb%dF=lb%dF*exp(-ds/Tf)+(lb%Fstat-lb%Fstat_Last)*exp(-ds/(2*Tf))
    ! lb%dCNv=lb%dCNv*exp(-ds/Tv)+lb%dcv*exp(-ds/(2*Tv))

    !! update lagged values
    !lb%CLRef_Last=lb%CLRef
    !lb%CLRefLE_Last=lb%CLRefLE
    !lb%Fstat_Last=lb%Fstat 
    !lb%cv_Last=lb%cv

    !End SUBROUTINE LB_UpdateStates

    !SUBROUTINE LB_LogicChecksum(LBCheck)
    !
    !    implicit none
    !    type(LB_type) :: lb
    !    integer :: LBCheck
    !    integer :: Loop

    !    ! Calculates a checksum for the logic states in the LB model

    !    LBCheck=0
    !    do Loop=1,9
    !        LBCheck=LBCheck+lb%LB_LogicOutputs(Loop)*lb%Logic_W(Loop)
    !    end do

    !End SUBROUTINE LB_LogicChecksum


end module
