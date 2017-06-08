#include "fdebug.h"

module dynstall
    
    use spud
    use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, PYTHON_FUNC_LEN, pi
    use airfoils

    implicit none
    
    type DS_Type 
    !=====================================
    ! These parameters follow the 
    ! formulation of turbinesFOAM

    logical :: StallFlag = .false.
    logical :: StallFlag_prev = .false. 
    logical :: do_calcAlphaEquiv=.false.
    real :: deltaS,time,time_prev
    real :: X,X_prev,Y,Y_prev,Z,Z_prev
    real :: Mach 
    real :: T1,T2,T3
    real :: Tp, Tf, Tv, Tvl ! Time periods
    real :: A1,A2,A3,b1,b2,b3,K1,K2,r,r0,alphaDS0DiffDeg
    real :: alphaEquiv, etaL, etaL_prev
    
    real :: D,D_prev,DF,DF_prev
    real :: CV,CV_prev,CNV,CNV_prev
    real :: fprime,fprime_prev,fDoublePrime
    real :: eta,S1,S2,CD0

    real :: tau, tau_prev
    integer :: nNewTimes
    real :: f ! not needed
    real :: fcrit
    real :: CMC
    real :: CNC, CNALpha
    real :: lambdaL,lambdaL_prev
    real :: lambdaM,lambdaM_prev
    real :: TI,E0

    real :: speed_of_sound
    real :: H,H_prev,J,J_prev
    real :: CNI,CMI,CNP,CNP_prev,CN1
    real :: CN_prev,CNprime
    real :: DP,DP_prev
    real :: alphaPrime,alphaPrime_prev
    real :: DAlpha, Dalpha_prev, TAlpha
    real :: Alpha, Alpha_prev, deltaAlpha,deltaAlpha_prev, alphaCrit
    real :: alpha1, alphaDS0, alphaSS, dalphaDS

    real :: Vx,CNF,CT,CN,cmFitExponent,K0,CM
    
    real, allocatable :: ReList(:), CNAlphaList(:), CD0List(:), CN1List(:), alpha1List(:), S1List(:), S2List(:)

    !=====================================
    end type DS_Type


    contains

    subroutine dystl_init(ds,dynstallpath)
        
        implicit none
        type(DS_Type),intent(inout) :: ds
        character(len=100),intent(in) :: dynstallpath
        character(len=80) :: stallfile
        character(1000) :: ReadLine
        integer :: i, Nstall

        ! READ PARAMETERS
        call get_option(trim(dynstallpath)//"Tp",ds%Tp,default=1.7)
        call get_option(trim(dynstallpath)//"Tf",ds%Tf,default=3.0)
        call get_option(trim(dynstallpath)//"TAlpha",ds%TAlpha,default=6.25)
        call get_option(trim(dynstallpath)//"alphaDS0DiffDeg",ds%alphaDS0DiffDeg,default=3.8)
        call get_option(trim(dynstallpath)//"r0",ds%r0,default=0.01)
        call get_option(trim(dynstallpath)//"Tv",ds%Tv,default=11.0)
        call get_option(trim(dynstallpath)//"Tvl",ds%Tvl,default=8.7)
        call get_option(trim(dynstallpath)//"B1",ds%B1,default=0.5)
        call get_option(trim(dynstallpath)//"eta",ds%eta,default=0.98)
        call get_option(trim(dynstallpath)//"E0",ds%E0,default=0.16)
    
        ! READ parameter or compute them from the Static foil Data
        call get_option(trim(dynstallpath)//"stall_data/file_name",StallFile)

        open(15,file=StallFile)
        ! Read the Number of Blades
        
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) NStall 

        allocate(ds%ReList(Nstall),ds%CNAlphaList(Nstall),ds%CD0List(Nstall),ds%CN1List(Nstall),ds%alpha1List(Nstall),ds%S1List(Nstall),ds%S2List(Nstall))
        ! Read the stations specs
        do i=1,NStall
        
        read(15,'(A)') ReadLine ! Stall parameters ....

        read(ReadLine,*) ds%ReList(i),ds%CNAlphaList(i),ds%CD0List(i),ds%CN1List(i),ds%alpha1List(i),ds%S1List(i),ds%S2List(i)

        end do
        
        close(15)
    
        ewrite(1,*) ds%S2List(:)
        stop
        ! SET ALL OTHER INITIAL CONDITIONS
        ds%X=0.0
        ds%X_prev=0.0
        ds%Y=0.0
        ds%Y_prev=0.0
        ds%Z=0.0
        ds%Z_prev=0.0
 
    end subroutine dystl_init
   
    subroutine calcAlphaEquiv(ds)
        
        implicit none
        ! Calculates the equivalent angle of attack 
        ! after having applied the defeciency functions
        type(DS_Type),intent(inout) :: ds
        real :: beta

        ds%T3=1.25*ds%Mach
        beta=1-ds%Mach**2
        
        ds%X=ds%X_prev*exp(-beta*ds%deltaS/ds%T1)+ds%A1*(ds%etaL-ds%etaL_prev)*exp(-beta*ds%deltaS/(2.0*ds%T1))
        ds%Y=ds%Y_prev*exp(-beta*ds%deltaS/ds%T2)+ds%A2*(ds%etaL-ds%etaL_prev)*exp(-beta*ds%deltaS/(2.0*ds%T2))
        ds%Z=ds%Z_prev*exp(-beta*ds%deltaS/ds%T3)+ds%A3*(ds%etaL-ds%etaL_prev)*exp(-beta*ds%deltaS/(2.0*ds%T3))

        ds%alphaEquiv = ds%alpha-ds%X-ds%Y-ds%Z

        if (abs(ds%alphaEquiv)>2*pi) then
            ds%alphaEquiv=mod(ds%alphaEquiv,2*pi)
        endif

    end subroutine calcAlphaEquiv
    
    subroutine DynstallCorrect(dynstall,airfoil,time,dt,Urel,chord,alpha,Re,CLdyn,CDdyn,CM25dyn)
        
        implicit none
        type(DS_Type), intent(inout) :: dynstall
        type(AirfoilType),intent(in) :: airfoil
        real, intent(in) :: time,dt,Urel,chord,alpha,Re
        real, intent(out) :: CLdyn,CDdyn,CM25dyn
        real :: mach

        ! update previous values if time has changed
        !if (time.ne.dynstall%time_prev) then
        !    dynstall%nNewTimes=dynstall%nNewTimes+1
        !    if (dynstall%nNewTimes > 1) then
        !       call update_DynStall(dynstall,time)
        !    end if
        !end if

        !if (dynstall%nNewTimes <=1) then
        !    dynstall%alpha_Prev=alpha
        !endif

        !dynstall%Alpha=alpha
        !dynstall%mach=urel/dynstall%speed_of_sound 
        !dynstall%deltaAlpha=dynstall%Alpha-dynstall%alpha_Prev
        !dynstall%deltaS=2*Urel*dt/chord

        !if (dynstall%do_calcAlphaEquiv) then
        !    call calcAlphaEquiv(dynstall)
        !else
        !    dynstall%alphaEquiv = dynstall%alpha
        !endif
        !! Evaluate static coefficient data if if has changed from 
        !! the Reynolds number correction
        !
        !call EvalStaticData(airfoil,dynstall,Re)
        
        !call calcUnsteady(lb,mach,chord,Urel,dt)
!        call calcSeparated(lb)
!
!        ! Modify Coefficients
!        CLdyn=lb%CN*cos(lb%alpha)+lb%CT*sin(lb%alpha)
!        CDdyn=lb%CN*sin(lb%alpha)-lb%CT*sin(lb%alpha)+lb%CD0
!        CM25dyn=lb%CM

        !return

    end subroutine DynstallCorrect
     
   ! subroutine calcUnsteady(lb,mach,chord,Ur,dt)
   !     implicit none
   !     type(LB_Type),intent(inout) :: lb
   !     real,intent(in) :: chord,Ur,dt,mach
   !     real :: kAlpha
   !     ! Calculate the circulatory normal force coefficient
   !     lb%CNC=lb%CNAlpha*lb%alphaEquiv

   !     ! Calculate the impulsive normal force coefficient
   !     kAlpha=0.75/(1.0-mach+pi*(1.0-mach**2)*mach*mach*(lb%A1*lb%b1+lb%A2*lb%b2))
   !     lb%TI=chord/lb%speed_of_sound

   !     lb%D=lb%D_prev*exp(-dt/(kAlpha*lb%TI))+(lb%deltaAlpha-lb%deltaAlpha_prev)*exp(-dt/(2.0*lb%TI))
   !     lb%CNI=4.0*kalpha*lb%TI/mach*(lb%deltaAlpha/dt-lb%D)

   !     ! Calculate total normal force coefficient
   !     lb%CNP = lb%CNC + lb%CNI

   !     ! Apply first-order lag to normal force coefficient
   !     lb%DP=lb%DP_prev*exp(-lb%ds/lb%Tp)+(lb%CNP-lb%CN_prev)*exp(-lb%ds/(2.0*lb%Tp))
   !     lb%CNprime=lb%CNP-lb%DP

   !     ! Calculate lagged angle of attack
   !     lb%alphaPrime=lb%CNPrime/lb%CNAlpha

   !     ! Set the stalled switch
   !     if(abs(lb%CNPrime)>lb%CN1) then
   !         lb%StallFlag=.true.
   !     endif
   !     
   ! end subroutine calcUnsteady

   ! subroutine calcSeparated(lb)
   !     implicit none
   !     type(LB_type),intent(inout):: lb
   !     real :: Tf,Tv, Tst, KN, m, cmf, cpv,cmv

   !     ! Calculate trailing-edge sepration point

   !     if (abs(lb%alphaPrime) < lb%alpha1) then
   !         lb%fprime=1.0-0.3*exp((abs(lb%alphaPrime) -lb%alpha1)/lb%S1) 
   !     else
   !         lb%fprime=0.04+0.66*exp((lb%alpha1 -abs(lb%alphaPrime))/lb%S2)
   !     endif

   !     Tf=lb%Tf

   !     ! Modify Tf time constant if necessary
   !     if (lb%tau > 0.and.lb%tau<=lb%Tvl) then
   !         Tf=0.5*lb%Tf
   !     else if (lb%tau>lb%Tvl.and.lb%tau<=2.0*lb%Tvl) then
   !         Tf=4.0*lb%Tf
   !     endif

   !     if(abs(lb%alpha) < abs(lb%Alpha_prev).and.abs(lb%CNPrime)<lb%CN1) then
   !         Tf=0.5*lb%Tf
   !     endif

   !     ! Calculate dynamic seperation point
   !     lb%DF=lb%DF_prev*exp(lb%ds/lb%Tf)+(lb%fprime-lb%fprime_prev)*exp(lb%ds/(2.0*lb%Tf))
   !     lb%fDoublePrime=lb%fprime-lb%DF

   !     if (lb%fDoublePrime <0) then
   !         lb%fDoublePrime = 0.0
   !     else if (lb%fDoublePrime >1) then
   !         lb%fDoublePrime =1.0
   !     endif

   !     ! Calculate normal force coefficient including dynamic separation point
   !     lb%CNF=lb%CNAlpha*lb%alphaEquiv*((1.0+sqrt(lb%fDoublePrime))/2.0)**2+lb%CNI

   !     ! Calculate tangential force coefficient
   !     if (lb%fDoublePrime <lb%fcrit) then
   !     lb%CT=lb%eta*lb%CNAlpha*lb%alphaEquiv*lb%alphaEquiv*(lb%fDoublePrime)**1.5
   !     else 
   !     lb%CT=lb%eta*lb%CNAlpha*lb%alphaEquiv*lb%alphaEquiv*sqrt(lb%fDoublePrime)
   !     endif

   !     ! Compute vortex shedding process if stalled 
   !     ! Evaluate vortex tracking time
   !     if (lb%StallFlag_prev.eqv..false.) then
   !         lb%tau=0.0
   !     else
   !         if (lb%tau.eq.lb%tau_prev) then
   !             lb%tau=lb%tau_prev+lb%ds 
   !         endif
   !     endif
   !     
   !     ! Calculate Strouhal number time constant and set tau to zero to 
   !     ! allow multiple vortex shedding
   !     Tst=2.0*(1.0-lb%fDoublePrime)/0.19
   !     if (lb%tau < lb%Tvl.and.(abs(lb%alpha)>abs(lb%alpha_prev))) then
   !         ! Halve Tv if dAlpha.dr changes sign
   !         if (lb%deltaAlpha*lb%deltaAlpha_prev<0) then
   !             Tv=0.5*lb%Tv
   !         endif
   !         KN = (1.0+sqrt(lb%fDoublePrime))**2/4.0
   !         lb%CV=lb%CNC*(1.0-KN)
   !         lb%CNV=lb%CNV_prev*exp(-lb%ds/Tv)+(lb%CV-lb%CV_prev)*exp(-lb%ds/(2.0*Tv))
   !     else
   !         Tv=0.5*lb%Tv
   !         lb%CV=0.0
   !         lb%CNV=lb%CNV_prev*exp(-lb%ds/Tv)
   !     endif

   !     ! Total normal force coefficient is the combination of that from 
   !     ! circulatory effects, impulsive effects, dynamic separation, and vortex
   !     ! lift
   !     lb%CN=lb%CNF+lb%CNV

   !     ! Calculate moment coefficient
   !     m=lb%cmFitExponent
   !     cmf=(lb%K0+lb%K1*(1-lb%fDoublePrime)+lb%K2*sin(pi*lb%fDoublePrime**m))*lb%CNC
   !     ! + moment coefficient at Zero lift angle of attack
   !     cpv =0.20*(1-cos(pi*lb%tau/lb%Tvl))
   !     cmv = lb%b2*(1.0-cos(pi*lb%tau/lb%Tvl))*lb%CNV
   !     lb%CM=cmf+cmv

   ! end subroutine calcSeparated
   ! 
   ! subroutine update_LBmodel(lb,time)
   !     implicit none
   !     type(LB_type), intent(inout) ::lb
   !     real,intent(in) :: time

   !     lb%time_prev=time
   !     lb%Alpha_prev=lb%Alpha
   !     lb%X_prev=lb%X
   !     lb%Y_prev=lb%Y
   !     lb%deltaAlpha_prev = lb%deltaAlpha
   !     lb%D_prev=lb%D
   !     lb%DP_prev=lb%DP
   !     lb%CNP_prev=lb%CNP
   !     lb%DF_prev=lb%DF
   !     lb%fprime_prev=lb%fprime
   !     lb%CV_prev=lb%CV
   !     lb%CNV_prev=lb%CNV
   !     lb%StallFlag_prev=lb%StallFlag
   !     lb%tau_prev=lb%tau

   !     return

   ! end subroutine update_LBmodel

   ! subroutine evalStaticData(airfoil,lb,Re)
   !     implicit none
   !     type(LB_Type),intent(inout) :: lb
   !     type(AirfoilType),intent(in) :: airfoil
   !     real, intent(in) :: Re
   !     real :: f,alphaSSP,alphaSSN,CLAlpha,CL0,CD0,CM250
   !     ! Get static stall angle in radians
   !     ! Evaluate static coefficient data if it has changed
   !     call EvalStaticStallParams(airfoil,Re,alphaSSP,alphaSSN,CLAlpha)  
   !            
   !     if (lb%alpha>0) then
   !         lb%AlphaSS = alphaSSP
   !     else
   !         lb%AlphaSS = alphaSSN
   !     endif

   !     lb%CNALpha=CLAlpha
   !     f=lb%fcrit 
   !     lb%alpha1 = lb%alphaSS*0.87
   !     lb%CN1 = lb%CNAlpha*lb%alpha1*((1.0+sqrt(f))/2.0)**2

   !     !> Get CL0,CD0,CM250
   !     call EvalStaticCoeff(Re,airfoil%alzer*condeg,CL0,CD0,CM250,airfoil)
   !     lb%CD0=CD0


   ! end subroutine evalStaticData

   ! subroutine LeishmanBeddoesCorrect(lb,airfoil,time,dt,Urel,chord,alpha,Re,CLdyn,CDdyn,CM25dyn)
   !     
   !     implicit none
   !     type(LB_Type), intent(inout) :: lb
   !     type(AirfoilType),intent(in) :: airfoil
   !     real, intent(in) :: time,dt,Urel,chord,alpha,Re
   !     real, intent(out) :: CLdyn,CDdyn,CM25dyn
   !     real :: mach

   !     ! update previous values if time has changed
   !     if (time.ne.lb%time_prev) then
   !         lb%nNewTimes=lb%nNewTimes+1
   !         if (lb%nNewTimes > 1) then
   !             call update_LBmodel(lb,time)
   !         end if
   !     end if

   !     if (lb%nNewTimes <=1) then
   !         lb%alpha_Prev=alpha
   !     endif

   !     lb%Alpha=alpha
   !     mach=urel/lb%speed_of_sound 
   !     lb%deltaAlpha=lb%Alpha-lb%alpha_Prev
   !     lb%ds=2*Urel*dt/chord

   !     if (lb%do_calcAlphaEquiv) then
   !         call calcAlphaEquiv(lb,mach)
   !     else
   !         lb%alphaEquiv = lb%alpha
   !     endif

   !     call EvalStaticData(airfoil,lb,Re)

   !     call calcUnsteady(lb,mach,chord,Urel,dt)
   !     call calcSeparated(lb)

   !     ! Modify Coefficients
   !     CLdyn=lb%CN*cos(lb%alpha)+lb%CT*sin(lb%alpha)
   !     CDdyn=lb%CN*sin(lb%alpha)-lb%CT*sin(lb%alpha)+lb%CD0
   !     CM25dyn=lb%CM

   !     return

   ! end subroutine LeishmanBeddoesCorrect

end module
