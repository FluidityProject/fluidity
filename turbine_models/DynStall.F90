#include "fdebug.h"

module dynstall

    use fldebug
    use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, PYTHON_FUNC_LEN, pi
    use airfoils

    implicit none
    
    type LB_Type
    
    !=====================================
    ! These parameters follow the 
    ! formulation of turbinesFOAM
    logical :: StallFlag = .false.
    logical :: StallFlag_prev = .false. 
    logical :: do_calcAlphaEquiv=.false.
    real :: ds,time,time_prev
    real :: X,X_prev,Y,Y_prev
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
    real :: CNI,CMI,CNP,CNP_prev,CN1
    real :: Tp,CN_prev,CNprime
    real :: DP,DP_prev
    real :: alphaPrime,alphaPrime_prev
    real :: DAlpha, Dalpha_prev, TAlpha
    real :: Alpha, Alpha_prev, deltaAlpha,deltaAlpha_prev, alphaCrit
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
        lb%A1=0.3
        lb%A2=0.7
        lb%b1=0.14
        lb%b2=0.53
        lb%alpha=0.0
        lb%DeltaAlpha=0.0
        lb%DeltaAlpha_prev=0.0
        lb%speed_of_sound=343 ! [m/s]
        lb%D=0.0
        lb%D_prev=0.0
        lb%DP=0.0
        lb%DP_prev=0.0
        lb%CNP=0.0
        lb%CNP_prev=0.0
        lb%fprime=0.0
        lb%fprime_prev=0.0
        lb%DF=0.0
        lb%DF_prev=0.0
        lb%CV=0.0
        lb%CV_prev=0.0
        lb%CNV=0.0
        lb%CNV_prev=0.0
        lb%eta=0.95
        lb%S1=1e12
        lb%S2=1e12
        lb%CD0=1e12
        lb%StallFlag_prev=.false.
        lb%Tp=1.7
        lb%Tf=3.0
        lb%Tv=6.0
        lb%Tvl=7.0
        lb%tau=0.0
        lb%tau_prev=0.0
        lb%nNewTimes=0
        lb%fcrit=0.6
        lb%K0=1e-6
        lb%K1=0.0
        lb%K2=0.0
        lb%cmFitExponent=2
        lb%CM=0.0
        
    end subroutine dystl_init_LB
    
    subroutine calcAlphaEquiv(lb,mach)
        implicit none
        ! Calculates the equivalent angle of attack 
        ! after having applied the defeciency functions
        type(LB_Type),intent(inout) :: lb
        real, intent(in) :: mach
        real :: beta

        beta=1-mach**2
        lb%X=lb%X_prev*exp(-beta*lb%ds/lb%T1)+lb%A1*lb%deltaAlpha*exp(-lb%b1*beta*lb%ds/2.0)
        lb%Y=lb%Y_prev*exp(-beta*lb%ds/lb%T2)+lb%A2*lb%deltaAlpha*exp(-lb%b2*beta*lb%ds/2.0)
        lb%alphaEquiv = lb%alpha-lb%X-lb%Y
        if (abs(lb%alphaEquiv)>2*pi) then
            lb%alphaEquiv=mod(lb%alphaEquiv,2*pi)
        endif

    end subroutine calcAlphaEquiv
         
    subroutine calcUnsteady(lb,mach,chord,Ur,dt)
        implicit none
        type(LB_Type),intent(inout) :: lb
        real,intent(in) :: chord,Ur,dt,mach
        real :: kAlpha
        ! Calculate the circulatory normal force coefficient
        lb%CNC=lb%CNAlpha*lb%alphaEquiv

        ! Calculate the impulsive normal force coefficient
        kAlpha=0.75/(1.0-mach+pi*(1.0-mach**2)*mach*mach*(lb%A1*lb%b1+lb%A2*lb%b2))
        lb%TI=chord/lb%speed_of_sound

        lb%D=lb%D_prev*exp(-dt/(kAlpha*lb%TI))+(lb%deltaAlpha-lb%deltaAlpha_prev)*exp(-dt/(2.0*lb%TI))
        lb%CNI=4.0*kalpha*lb%TI/mach*(lb%deltaAlpha/dt-lb%D)

        ! Calculate total normal force coefficient
        lb%CNP = lb%CNC + lb%CNI

        ! Apply first-order lag to normal force coefficient
        lb%DP=lb%DP_prev*exp(-lb%ds/lb%Tp)+(lb%CNP-lb%CN_prev)*exp(-lb%ds/(2.0*lb%Tp))
        lb%CNprime=lb%CNP-lb%DP

        ! Calculate lagged angle of attack
        lb%alphaPrime=lb%CNPrime/lb%CNAlpha

        ! Set the stalled switch
        if(abs(lb%CNPrime)>lb%CN1) then
            lb%StallFlag=.true.
        endif
        
    end subroutine calcUnsteady

    subroutine calcSeparated(lb)
        implicit none
        type(LB_type),intent(inout):: lb
        real :: Tf,Tv, Tst, KN, m, cmf, cpv,cmv

        ! Calculate trailing-edge sepration point

        if (abs(lb%alphaPrime) < lb%alpha1) then
            lb%fprime=1.0-0.3*exp((abs(lb%alphaPrime) -lb%alpha1)/lb%S1) 
        else
            lb%fprime=0.04+0.66*exp((lb%alpha1 -abs(lb%alphaPrime))/lb%S2)
        endif

        Tf=lb%Tf

        ! Modify Tf time constant if necessary
        if (lb%tau > 0.and.lb%tau<=lb%Tvl) then
            Tf=0.5*lb%Tf
        else if (lb%tau>lb%Tvl.and.lb%tau<=2.0*lb%Tvl) then
            Tf=4.0*lb%Tf
        endif

        if(abs(lb%alpha) < abs(lb%Alpha_prev).and.abs(lb%CNPrime)<lb%CN1) then
            Tf=0.5*lb%Tf
        endif

        ! Calculate dynamic seperation point
        lb%DF=lb%DF_prev*exp(lb%ds/lb%Tf)+(lb%fprime-lb%fprime_prev)*exp(lb%ds/(2.0*lb%Tf))
        lb%fDoublePrime=lb%fprime-lb%DF

        if (lb%fDoublePrime <0) then
            lb%fDoublePrime = 0.0
        else if (lb%fDoublePrime >1) then
            lb%fDoublePrime =1.0
        endif

        ! Calculate normal force coefficient including dynamic separation point
        lb%CNF=lb%CNAlpha*lb%alphaEquiv*((1.0+sqrt(lb%fDoublePrime))/2.0)**2+lb%CNI

        ! Calculate tangential force coefficient
        if (lb%fDoublePrime <lb%fcrit) then
        lb%CT=lb%eta*lb%CNAlpha*lb%alphaEquiv*lb%alphaEquiv*(lb%fDoublePrime)**1.5
        else 
        lb%CT=lb%eta*lb%CNAlpha*lb%alphaEquiv*lb%alphaEquiv*sqrt(lb%fDoublePrime)
        endif

        ! Compute vortex shedding process if stalled 
        ! Evaluate vortex tracking time
        if (lb%StallFlag_prev.eqv..false.) then
            lb%tau=0.0
        else
            if (lb%tau.eq.lb%tau_prev) then
                lb%tau=lb%tau_prev+lb%ds 
            endif
        endif
        
        ! Calculate Strouhal number time constant and set tau to zero to 
        ! allow multiple vortex shedding
        Tst=2.0*(1.0-lb%fDoublePrime)/0.19
        if (lb%tau < lb%Tvl.and.(abs(lb%alpha)>abs(lb%alpha_prev))) then
            ! Halve Tv if dAlpha.dr changes sign
            if (lb%deltaAlpha*lb%deltaAlpha_prev<0) then
                Tv=0.5*lb%Tv
            endif
            KN = (1.0+sqrt(lb%fDoublePrime))**2/4.0
            lb%CV=lb%CNC*(1.0-KN)
            lb%CNV=lb%CNV_prev*exp(-lb%ds/Tv)+(lb%CV-lb%CV_prev)*exp(-lb%ds/(2.0*Tv))
        else
            Tv=0.5*lb%Tv
            lb%CV=0.0
            lb%CNV=lb%CNV_prev*exp(-lb%ds/Tv)
        endif

        ! Total normal force coefficient is the combination of that from 
        ! circulatory effects, impulsive effects, dynamic separation, and vortex
        ! lift
        lb%CN=lb%CNF+lb%CNV

        ! Calculate moment coefficient
        m=lb%cmFitExponent
        cmf=(lb%K0+lb%K1*(1-lb%fDoublePrime)+lb%K2*sin(pi*lb%fDoublePrime**m))*lb%CNC
        ! + moment coefficient at Zero lift angle of attack
        cpv =0.20*(1-cos(pi*lb%tau/lb%Tvl))
        cmv = lb%b2*(1.0-cos(pi*lb%tau/lb%Tvl))*lb%CNV
        lb%CM=cmf+cmv

    end subroutine calcSeparated
    
    subroutine update_LBmodel(lb,time)
        implicit none
        type(LB_type), intent(inout) ::lb
        real,intent(in) :: time

        lb%time_prev=time
        lb%Alpha_prev=lb%Alpha
        lb%X_prev=lb%X
        lb%Y_prev=lb%Y
        lb%deltaAlpha_prev = lb%deltaAlpha
        lb%D_prev=lb%D
        lb%DP_prev=lb%DP
        lb%CNP_prev=lb%CNP
        lb%DF_prev=lb%DF
        lb%fprime_prev=lb%fprime
        lb%CV_prev=lb%CV
        lb%CNV_prev=lb%CNV
        lb%StallFlag_prev=lb%StallFlag
        lb%tau_prev=lb%tau

        return

    end subroutine update_LBmodel

    subroutine evalStaticData(airfoil,lb,Re)
        implicit none
        type(LB_Type),intent(inout) :: lb
        type(AirfoilType),intent(in) :: airfoil
        real, intent(in) :: Re
        real :: f,alphaSSP,alphaSSN,CLAlpha,CL0,CD0,CM250
        ! Get static stall angle in radians
        ! Evaluate static coefficient data if it has changed
        call EvalStaticStallParams(airfoil,Re,alphaSSP,alphaSSN,CLAlpha)  
               
        if (lb%alpha>0) then
            lb%AlphaSS = alphaSSP
        else
            lb%AlphaSS = alphaSSN
        endif

        lb%CNALpha=CLAlpha
        f=lb%fcrit 
        lb%alpha1 = lb%alphaSS*0.87
        lb%CN1 = lb%CNAlpha*lb%alpha1*((1.0+sqrt(f))/2.0)**2

        !> Get CL0,CD0,CM250
        call EvalStaticCoeff(Re,airfoil%alzer*condeg,CL0,CD0,CM250,airfoil)
        lb%CD0=CD0


    end subroutine evalStaticData

    subroutine LeishmanBeddoesCorrect(lb,airfoil,time,dt,Urel,chord,alpha,Re,CLdyn,CDdyn,CM25dyn)
        
        implicit none
        type(LB_Type), intent(inout) :: lb
        type(AirfoilType),intent(in) :: airfoil
        real, intent(in) :: time,dt,Urel,chord,alpha,Re
        real, intent(out) :: CLdyn,CDdyn,CM25dyn
        real :: mach

        ! update previous values if time has changed
        if (time.ne.lb%time_prev) then
            lb%nNewTimes=lb%nNewTimes+1
            if (lb%nNewTimes > 1) then
                call update_LBmodel(lb,time)
            end if
        end if

        if (lb%nNewTimes <=1) then
            lb%alpha_Prev=alpha
        endif

        lb%Alpha=alpha
        mach=urel/lb%speed_of_sound 
        lb%deltaAlpha=lb%Alpha-lb%alpha_Prev
        lb%ds=2*Urel*dt/chord

        if (lb%do_calcAlphaEquiv) then
            call calcAlphaEquiv(lb,mach)
        else
            lb%alphaEquiv = lb%alpha
        endif

        call EvalStaticData(airfoil,lb,Re)

        call calcUnsteady(lb,mach,chord,Urel,dt)
        call calcSeparated(lb)

        ! Modify Coefficients
        CLdyn=lb%CN*cos(lb%alpha)+lb%CT*sin(lb%alpha)
        CDdyn=lb%CN*sin(lb%alpha)-lb%CT*sin(lb%alpha)+lb%CD0
        CM25dyn=lb%CM

        return

    end subroutine LeishmanBeddoesCorrect

end module
