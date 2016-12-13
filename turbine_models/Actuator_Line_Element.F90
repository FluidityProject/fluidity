#include "fdebug.h"

module actuator_line_element
  
    use fldebug
    use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, PYTHON_FUNC_LEN, pi
    use futils 
    use airfoils
    use actuator_line_model_utils 
    use dynstall

type ActuatorLineType
    integer :: NElem                    ! Number of Elements of the Blade
    character(len=100):: name           ! Actuator line name
    character(len=100):: geom_file      ! Actuator line file name (is not used for the turbines)

    ! Station parameters
    logical :: FlipN =.false.           ! Flip Normal
    real, allocatable :: QCx(:)         ! Blade quarter-chord line x coordinates at element ends
    real, allocatable :: QCy(:)         ! Blade quarter-chord line y coordinates at element ends
    real, allocatable :: QCz(:)         ! Blade quarter-chord line z coordinates at element ends
    real, allocatable :: tx(:)          ! Blade unit tangent vector (rearward chord line direction) x-comp at element ends
    real, allocatable :: ty(:)          ! Blade unit tangent vector (rearward chord line direction) y-comp at element ends 
    real, allocatable :: tz(:)          ! Blade unit tangent vector (rearward chord line direction) z-comp at element ends  
    real, allocatable :: C(:)           ! Blade chord length at element ends
    real, allocatable :: thick(:)       ! Blade thickness at element ends
    real, allocatable :: pitch(:)       ! Blade station pitch at element ends

    ! Element parameters 
    real, allocatable :: PEx(:)         ! Element centre x coordinates
    real, allocatable :: PEy(:)         ! Element centre y coordinates
    real, allocatable :: PEz(:)         ! Element centre z coordinates
    real, allocatable :: tEx(:)         ! Element unit tangent vector (rearward chord line direction) x-component
    real, allocatable :: tEy(:)         ! Element unit tangent vector (rearward chord line direction) y-component
    real, allocatable :: tEz(:)         ! Element unit tangent vector (rearward chord line direction) z-component
    real, allocatable :: nEx(:)         ! Element unit normal vector x-component
    real, allocatable :: nEy(:)         ! Element unit normal vector y-component
    real, allocatable :: nEz(:)         ! Element unit normal vector z-component
    real, allocatable :: sEx(:)         ! Element unit spanwise vector x-component 
    real, allocatable :: sEy(:)         ! Element unit spanwise vector y-component
    real, allocatable :: sEz(:)         ! Element unit spanwise vector z-component
    real, allocatable :: EC(:)          ! Element chord lenght
    real, allocatable :: EDS(:)         ! Element spanwise distance (length)
    real, allocatable :: EArea(:)       ! Element Area
    real, allocatable :: Eepsilon(:)    ! Element Force Projection Parameter
    real, allocatable :: ERdist(:)      ! Element Distance from the origin 
    real, allocatable :: ETtoC(:)       ! Element thickness to Chord ratio
    
    ! Angle of Attack, Pitch, local Reynolds number and relative velocity
    real, allocatable :: Epitch(:)      ! Element pitch angle
    real, allocatable :: EAOA(:)        ! Element Current angle of Attack (used in added mass terms)
    real, allocatable :: EAOAdot(:)     ! Element AOA rate of change
    real, allocatable :: EUn(:)         ! Element Current normal velocity (used in added mass terms)
    real, allocatable :: EAOA_LAST(:)   ! Element Last angle of Attack (used in added mass terms)
    real, allocatable :: EUn_LAST(:)    ! Element Last normal velocity (used in added mass terms)
    real, allocatable :: ERe(:)         ! Element local Reynolds number
    real, allocatable :: EUr(:)         ! Element local Reynolds number
    
    ! Velocity of the Fluid at the actuator Line Locations
    real, allocatable :: EVx(:)         ! Element Local fluid Velocity in the global x-direction
    real, allocatable :: EVy(:)         ! Element Local fluid Velocity in the global y-direction
    real, allocatable :: EVz(:)         ! Element Local fluid Velocity in the global z-direction
     
    ! Body Velocity of the Actuator Line 
    real, allocatable :: EVbx(:)        ! Element Local body Velocity in the global x-direction
    real, allocatable :: EVby(:)       ! Element Local body Velocity in the global y-direction
    real, allocatable :: EVbz(:)        ! Element Local body Velocity in the global z-direction
    real, allocatable :: EObx(:)        ! Element Local body angular Velocity in the global x-direction
    real, allocatable :: EOby(:)       ! Element Local body angular Velocity in the global y-direction
    real, allocatable :: EObz(:)        ! Element Local body angular Velocity in the global z-direction
    
    ! Element Forces CD, CL CM25 
    real, allocatable :: ECD(:)         ! Element Drag Coefficient
    real, allocatable :: ECL(:)         ! Element Lift Coefficient 
    real, allocatable :: ECM(:)         ! Element Moment Coefficient
    real, allocatable :: ECN(:)         ! Element Normal Force Coefficient
    real, allocatable :: ECT(:)         ! Element Tangential Force Coefficient 
    
    ! Element Forces in the nts direction
    real, allocatable :: EFn(:)         ! Element Force in the normal direction
    real, allocatable :: EFt(:)         ! Element Force in the tangential direction (rearward chord line direction) 
    real, allocatable :: EMS(:)         ! Element Moment 

    ! Element Forces and Torque in the xyz direction
    real, allocatable :: EFx(:)         ! Element Force in the global x-direction
    real, allocatable :: EFy(:)         ! Element Force in the global y-direction
    real, allocatable :: EFz(:)         ! Element Force in the global z-direction
    
    ! Element Airfoil Data
    type(AirfoilType), allocatable :: EAirfoil(:) ! Element Airfoil 
    integer :: NAirfoilData
    type(AirfoilType), allocatable :: AirfoilData(:) ! Element Airfoil 
    type(LB_Type), allocatable :: E_LB_Model(:)   ! Element Leishman-Beddoes Model
    
    ! Forces and Torques on the ActuatorLine 
    real :: Fx     ! Element Force in the global x-direction
    real :: Fy     ! Element Force in the global y-direction
    real :: Fz     ! Element Force in the global z-direction

    real :: Area   ! Effective Airfoil Area
    real :: L 
    ! Degrees of Freedom
    
    real :: COR(3)       ! Center of Rotation
    real :: SpanWise(3)   ! Point of Rotation
   
    ! Unsteady Loading
    logical :: do_added_mass=.false.
    logical :: do_dynamic_stall=.false.
    logical :: do_tip_correction=.false.

    ! Blade/Actuator_line Pitch
    logical :: pitch_control=.false.
    real    :: pitch_start_time
    real    :: pitch_end_time

    ! Harmonic Pitch control parameters
    real    :: angular_pitch_freq
    real    :: pitch_angle_init
    real    :: pitchAmp

end type ActuatorLineType

    contains

    subroutine set_actuatorline_geometry(actuatorline)

    implicit none
    type(ActuatorLineType),intent(inout) :: actuatorline
    real, allocatable :: rR(:),ctoR(:),pitch(:),thick(:)
    real :: SVec(3), theta, origin(3), length 
    integer :: Nstations, iact, Istation, ielem

    ewrite(2,*) 'Entering set_actuatorline_geometry'

    ewrite(2,*) 'Actuatorline Name : ', actuatorline%name 
    ewrite(2,*) '============================='

    call read_actuatorline_geometry(actuatorline%geom_file,length,SVec,rR,ctoR,pitch,thick,Nstations)
    
    call allocate_actuatorline(actuatorline,Nstations)
   
    actuatorline%SpanWise=SVec
    actuatorline%COR=(/0.0, 0.0, 0.0/)
    actuatorline%L=length

    ! The directions of vectors etc are just hacked ...
    actuatorline%Nelem=Nstations-1 
    do istation=1,Nstations
    actuatorline%QCx(istation)=rR(istation)*length*Svec(1)
    actuatorline%QCy(istation)=rR(istation)*length*Svec(2)
    actuatorline%QCz(istation)=rR(istation)*length*Svec(3)      
    
    if(actuatorline%pitch_control) then
         pitch(istation)=actuatorline%pitch_angle_init   
    endif

    actuatorline%tx(istation)= cos(pitch(istation)/180.0*pi)    
    actuatorline%ty(istation)= 0.0
    actuatorline%tz(istation)= -sin(pitch(istation)/180.0*pi)     

    actuatorline%C(istation)=ctoR(istation)*length
    actuatorline%thick(istation)=thick(istation)
    actuatorline%FlipN=.true.
    end do

    call make_actuatorline_geometry(actuatorline)
    
    ! Compute the Area
    do ielem=1,actuatorline%Nelem
    actuatorline%Area=actuatorline%Area+actuatorline%EArea(ielem)
    end do
    
    ! Initial Values for the body linear and angular velocities 
    actuatorline%EVbx(:)=0.0
    actuatorline%EVby(:)=0.0
    actuatorline%EVbz(:)=0.0
    actuatorline%EObx(:)=0.0
    actuatorline%EOby(:)=0.0
    actuatorline%EObz(:)=0.0
    actuatorline%Eepsilon(:)=0.0

    
    ! Populate element Airfoils 
    call populate_blade_airfoils(actuatorline%NElem,actuatorline%NAirfoilData,actuatorline%EAirfoil,actuatorline%AirfoilData,actuatorline%ETtoC)
    
    actuatorline%EAOA_LAST(:)=-666.0
    
    do ielem=1,actuatorline%Nelem
    call dystl_init_LB(actuatorline%E_LB_Model(ielem))
    end do
    
    ewrite(2,*) 'Exiting set_actuatorline_geometry'

    end subroutine set_actuatorline_geometry
    

    subroutine Compute_ActuatorLine_Forces(act_line,visc,dt)
       
    implicit none
    type(ActuatorLineType),intent(inout) :: act_line
    real,intent(in) ::visc,dt
    real :: R(3)
    real :: wRotX,wRotY,wRotZ,Rx,Ry,Rz,ub,vb,wb,u,v,w
    real :: xe,ye,ze,nxe,nye,nze,txe,tye,tze,sxe,sye,sze,ElemArea,ElemChord
    real :: xe5,ye5,ze5,xe75,ye75,ze75,ub5,vb5,wb5,ub75,vb75,wb75,urdn5,urdn75 
    real :: urdn,urdc,wP,ur,alpha,Re,alpha5,alpha75,adotnorm, wPnorm,ur5,ur75
    real :: CL,CD,CN,CT,CLCirc,CM25,MS,FN,FT,FS,FX,Fy,Fz,te, F1, g1
    real :: TRx,TRy,TRz, TRn,TRt, TRs, RotX,RotY,RotZ, dal, relem, ds
    integer :: ielem
  
    ewrite(2,*) 'Entering Compute_Forces '
     
    !===========================================================
    ! Assign global values to local values (to make life easier
    !==========================================================
    do ielem=1,act_line%NElem
    
    xe=act_line%PEX(ielem)
    ye=act_line%PEY(ielem)
    ze=act_line%PEZ(ielem)

    nxe=act_line%nEx(ielem)
    nye=act_line%nEy(ielem)
    nze=act_line%nEz(ielem)
    txe=act_line%tEx(ielem)
    tye=act_line%tEy(ielem)
    tze=act_line%tEz(ielem)
    sxe=act_line%sEx(ielem)
    sye=act_line%sEy(ielem)
    sze=act_line%sEz(ielem)
    ElemArea=act_line%EArea(ielem)
    ElemChord=act_line%EC(ielem) 
    u=act_line%EVx(ielem)
    v=act_line%EVy(ielem)
    w=act_line%EVz(ielem) 

    !=====================================
    ! Solid Body Rotation of the elements
    !=====================================
    ub=act_line%EVbx(ielem)
    vb=act_line%EVby(ielem)
    wb=act_line%EVbz(ielem)
    wRotx=act_line%EObx(ielem)
    wRoty=act_line%EOby(ielem)
    wRotz=act_line%EObz(ielem)

    !==============================================================
    ! Calculate element normal and tangential velocity components. 
    !==============================================================
    urdn=nxe*(u-ub)+nye*(v-vb)+nze*(w-wb)! Normal 
    urdc=txe*(u-ub)+tye*(v-vb)+tze*(w-wb)! Tangential
    wP = sxe*wRotX+sye*wRotY+sze*wRotZ
    
    ur=sqrt(urdn**2.0+urdc**2.0)
    act_line%EUr(ielem)=ur
    
    ! This is the dynamic angle of attack 
    alpha=atan2(urdn,urdc)
    act_line%ERe(ielem) = ur*ElemChord/Visc
   
    wPNorm = wP*ElemChord/(2.0*max(ur,0.001))

    ! Calculate half chord and 75% chord velocites to be used in the pitch rate effects
    if(act_line%do_added_mass) then
    xe5=xe+0.25*ElemChord*txe
    ye5=ye+0.25*ElemChord*tye
    ze5=ze+0.25*ElemChord*tze
    xe75=xe+0.5*ElemChord*txe
    ye75=ye+0.5*ElemChord*tye
    ze75=ze+0.5*ElemChord*tze
    call cross(wRotX,wRotY,wRotZ,xe5,ye5,ze5,ub5,vb5,wb5)
    call cross(wRotX,wRotY,wRotZ,xe75,ye75,ze75,ub75,vb75,wb75)
    urdn5=nxe*(u-ub5)+nye*(v-vb5)+nze*(w-wb5)! Normal 
    ur5=sqrt(urdn5**2+urdc**2)
    alpha5=atan2(urdn5,urdc)
    urdn75=nxe*(u-ub75)+nye*(v-vb75)+nze*(w-wb75)! Normal 
    ur75=sqrt(urdn75**2+urdc**2)
    alpha75=atan2(urdn75,urdc)
    else
    alpha5=alpha
    alpha75=alpha
    wPNorm=0.0
    endif
    
    act_line%EAOA(ielem)=alpha75
    
    if(act_line%EAOA_Last(ielem)<0) then
    dal=0.0
    else
    dal=act_line%EAOA(ielem)-act_line%EAOA_Last(ielem)
    endif
    
    act_line%EAOAdot(ielem)=dal/dt
    
    adotnorm=act_line%EAOAdot(ielem)*ElemChord/(2.0*max(ur,0.001)) ! adot*c/(2*U)
    
    !====================================
    ! Compute the Aerofoil Coefficients
    !====================================
    call compute_aeroCoeffs(act_line%do_dynamic_stall,act_line%do_added_mass,act_line%EAirfoil(ielem),act_line%E_LB_Model(ielem),alpha75,alpha5,act_line%ERe(ielem),wPNorm,adotnorm,CN,CT,CM25,CL,CLCIrc,CD)  
    !===================================
    ! Update Dynamic Stall Model 
    !===================================
    if(act_line%do_dynamic_stall) then
    ds=2.0*ur*dt/ElemChord
    call LB_UpdateStates(act_line%E_LB_MODEL(ielem),ds)
    endif 
    !========================================================
    ! Apply Coeffs to calculate tangential and normal Forces
    !========================================================
    FN=0.5*CN*ElemArea*ur**2.0
    FT=0.5*CT*ElemArea*ur**2.0
    MS=0.5*CM25*ElemChord*ElemArea*ur**2.0
    !===============================================
    ! Compute forces in the X, Y, Z axis and torque  
    !===============================================
    FX=FN*nxe+FT*txe
    FY=FN*nye+FT*tye
    FZ=FN*nze+FT*tze
    !==========================================
    ! Assign the derived types
    !==========================================
    ! Local Load Coefficients
    act_line%ECD(ielem)=CD
    act_line%ECL(ielem)=CL
    act_line%ECM(ielem)=CM25
    act_line%ECN(ielem)=CN
    act_line%ECT(ielem)=CT
    
    ! Local Coordinate-system Forces
    act_line%EFN(ielem)=FN
    act_line%EFT(ielem)=FT
    act_line%EMS(ielem)=MS
    ! Global Forces and Torques
    act_line%EFX(ielem)=FX
    act_line%EFY(ielem)=FY
    act_line%EFZ(ielem)=FZ
    !===============================================
    !! Set the AOA_LAST before exiting the routine
    !===============================================
    act_line%EAOA_LAST(ielem)=alpha75 
    act_line%EUn_last(ielem)=urdn 
    end do 

    ewrite(2,*) 'Exiting Compute_Forces'

    end subroutine compute_Actuatorline_Forces

    subroutine compute_Tower_Forces(tower,visc,time,CL,CD,Str)
        implicit none
        type(ActuatorLineType),intent(inout) :: tower
        real,intent(in) ::visc, time, CL, CD, Str
        real :: R(3),rand(3000)
        real :: xe,ye,ze,nxe,nye,nze,txe,tye,tze,sxe,sye,sze,ElemArea,ElemChord
        real :: u,v,w,ub,vb,wb,urdn,urdc,ur,Diameter,freq,alpha
        real :: CN,CT,CLCirc,CM25,MS,FN,FT,FS,FX,Fy,Fz
        integer :: ielem
    
        ewrite(2,*) 'Entering compute_tower_forces'
        call random_number(rand)
        do ielem=1,tower%NElem
    
            xe=tower%PEX(ielem)
            ye=tower%PEY(ielem)
            ze=tower%PEZ(ielem)
            nxe=tower%nEx(ielem)
            nye=tower%nEy(ielem)
            nze=tower%nEz(ielem)
            txe=tower%tEx(ielem)
            tye=tower%tEy(ielem)
            tze=tower%tEz(ielem)
            sxe=tower%sEx(ielem)
            sye=tower%sEy(ielem)
            sze=tower%sEz(ielem)
            Diameter=tower%EC(ielem) 
            ElemArea=tower%EArea(ielem)
            u=tower%EVx(ielem)
            v=tower%EVy(ielem)
            w=tower%EVz(ielem) 

            ub=0.0
            vb=0.0
            wb=0.0
           
            !==============================================================
            ! Calculate element normal and tangential velocity components. 
            !==============================================================
            urdn=nxe*(u-ub)+nye*(v-vb)+nze*(w-wb)! Normal 
            urdc=txe*(u-ub)+tye*(v-vb)+tze*(w-wb)! Tangential
            ur=sqrt(urdn**2.0+urdc**2.0)
            tower%EUr(ielem)=ur
            alpha=atan2(urdn,urdc)
            tower%EAOA(ielem)=alpha
            tower%ERE(ielem)=ur*Diameter/visc
            freq=Str*ur/max(Diameter,0.0001)
            tower%ECL(ielem)=CL*sin(2.0*freq*pi*time)
            tower%ECL(ielem)=tower%ECL(ielem)*(1.0+0.25*(-1.0+2.0*rand(ielem)))
            tower%ECD(ielem)=CD
            CN=tower%ECL(ielem)*cos(alpha)+tower%ECD(ielem)*sin(alpha)                                   
            CT=-tower%ECL(ielem)*sin(alpha)+tower%ECD(ielem)*cos(alpha) 
            !========================================================
            ! Apply Coeffs to calculate tangential and normal Forces
            !========================================================
            FN=0.5*CN*ElemArea*ur**2.0
            FT=0.5*CT*ElemArea*ur**2.0
            !===============================================
            ! Compute forces in the X, Y, Z axis and torque  
            !===============================================
            FX=FN*nxe+FT*txe
            FY=FN*nye+FT*tye
            FZ=FN*nze+FT*tze
            !==========================================
            ! Assign the derived types
            !==========================================
            tower%ECN(ielem)=CN
            tower%ECT(ielem)=CT
            tower%EFN(ielem)=FN
            tower%EFT(ielem)=FT
            tower%EMS(ielem)=0.0
            tower%EFX(ielem)=FX
            tower%EFY(ielem)=FY
            tower%EFZ(ielem)=FZ
         
        enddo

        ewrite(2,*) 'Exiting compute_tower_forces'

    end subroutine compute_Tower_Forces

    subroutine pitch_actuator_line(act_line)

    implicit none
    type(ActuatorLineType),intent(INOUT) :: act_line
    real :: pitch_angle !Pitch in degrees
    real :: R(3,3), t(3,1)
    integer :: istation, Nstation
    
    !> Change the pitch angle by changing n,t and s unit vectors
    Nstation=act_line%Nelem+1

    do istation=1,Nstation 

    pitch_angle=act_line%pitch(istation)
    !> Define the Pitch Rotation Matrix
    R=reshape((/cos(pitch_angle*pi/180.0),0.0,-sin(pitch_angle*pi/180.0),0.0,1.0,0.0,sin(pitch_angle*pi/180.0),0.0,cos(pitch_angle*pi/180.0)/),(/3,3/))
    t(1,1)=act_line%tx(istation)
    t(2,1)=act_line%ty(istation)
    t(3,1)=act_line%tz(istation)     
    t=matmul(R,t)
    !>Reassign the tangential vector
    act_line%tx(istation)=cos(pitch_angle*pi/180.0)
    act_line%ty(istation)=0.0
    act_line%tz(istation)=-sin(pitch_angle*pi/180.0)
    end do

    call make_actuatorline_geometry(act_line) 

    end subroutine pitch_actuator_line 

    subroutine populate_blade_airfoils(NElem,NData,EAirfoil,AirfoilData,ETtoC)

    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    ! This routine initialises the airfoil struct for the blades
    ! by interpolating from the data
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    implicit none
    integer,intent(IN) :: NElem,NData
    type(AirfoilType),dimension(NElem) :: EAirfoil
    type(AirfoilType),dimension(NData):: AirfoilData
    real,dimension(NElem),intent(in) ::  ETtoC
    real,dimension(NData) ::  Thicks, diffthicks
    integer :: ielem,idata,imin,imax,iint

    ewrite(2,*) 'Entering populate_blade_airfoils'    
    ! We need to interpolate from two or more 
    
    do ielem=1,NElem
    EAirfoil(ielem)%afname = int2str(ielem)
    call allocate_airfoil(EAirfoil(ielem),MaxAOAVals,MaxReVals) 
    
    Thicks(:)=0.0
    diffthicks(:)=0.0
    imin=0
    imax=0
    iint=0
    
        do idata=1,NData
             Thicks(idata)=AirfoilData(idata)%tc
             diffthicks(idata)=abs(AirfoilData(idata)%tc-ETtoC(ielem))
        end do 
        imin=minloc(Thicks,1)
        imax=maxloc(Thicks,1)
        iint=minloc(diffthicks,1)
        
        if(ETtoC(ielem)>=Thicks(imax)) then
             call copy_airfoil_values(EAirfoil(ielem),AirfoilData(imax))
        elseif(ETtoC(ielem)<=Thicks(imin)) then
             call copy_airfoil_values(EAirfoil(ielem),AirfoilData(imin))
        else
             call copy_airfoil_values(EAirfoil(ielem),AirfoilData(iint))
        endif

    end do 

    ewrite(2,*) 'Exiting populate_blade_airfoils'

    end subroutine populate_blade_airfoils
    
    subroutine rotate_actuatorline(actuatorline,origin,rotN,theta)

    implicit none
    type(ActuatorLineType),intent(inout) :: actuatorline
    real,intent(in) :: origin(3), rotN(3)
    real :: theta,nrx,nry,nrz,px,py,pz 
    integer :: j,ielem,i
    real :: vrx,vry,vrz,VMag
    real :: xtmp,ytmp,ztmp, txtmp, tytmp, tztmp
    ! Rotates data in blade arrays. Rotate element end geometry and recalculate element geometry.

    ewrite(1,*) 'Entering rotate_actuatorline'
        
    ! Specify the rotation axis and the normal vector of rotation
    
    nrx=rotN(1)
    nry=RotN(2)
    nrz=RotN(3)
    
    px=origin(1)
    py=origin(2)
    pz=origin(3)


    do ielem=1,actuatorline%Nelem+1
    ! Blade end locations (quarter chord). xBE(MaxSegEnds)
        xtmp=actuatorline%QCx(ielem)
        ytmp=actuatorline%QCy(ielem)
        ztmp=actuatorline%QCz(ielem)
    
        Call QuatRot(xtmp,ytmp,ztmp,theta,nrx,nry,nrz,px,py,pz,vrx,vry,vrz)
        actuatorline%QCx(ielem)=vrx                                       
        actuatorline%QCy(ielem)=vry                                       
        actuatorline%QCz(ielem)=vrz                                  
    
        txtmp=actuatorline%tx(ielem)
        tytmp=actuatorline%ty(ielem)
        tztmp=actuatorline%tz(ielem)
    
    ! Tangent vectors
        Call QuatRot(txtmp,tytmp,tztmp,theta,nrx,nry,nrz,px,py,pz,vrx,vry,vrz)
        VMag=sqrt(vrx**2+vry**2+vrz**2)
        actuatorline%tx(ielem)=vrx/VMag                                      
        actuatorline%ty(ielem)=vry/VMag                                  
        actuatorline%tz(ielem)=vrz/VMag                                       
  
    end do
    
    ewrite(1,*) 'Exiting rotate_actuatorline'

    end subroutine rotate_actuatorline

    subroutine read_actuatorline_geometry(FN,Rmax,SpanwiseVec,rR,ctoR,pitch,thick,Nstations)
    
    implicit none
    character(len=100),intent(in)  :: FN ! FileName of the geometry file
    real,dimension(3) :: SpanwiseVec
    real, allocatable,intent(out) :: rR(:),ctoR(:),pitch(:),thick(:)  
    real, intent(out) :: Rmax  
    integer, intent(out) :: Nstations
    integer :: i,j
    character(1000) :: ReadLine
    
    open(15,file=FN)
    ! Read the Number of Blades
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) Rmax 
    
    !Read Spanwise actuator line axis
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) SpanwiseVec(1), SpanwiseVec(2), SpanwiseVec(3)
     
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) Nstations

    allocate(rR(Nstations),ctoR(Nstations),pitch(Nstations),thick(Nstations))
    ! Read the stations specs
    do i=1,NStations
    
    read(15,'(A)') ReadLine ! Blade ....

    read(ReadLine,*) rR(i), ctoR(i), pitch(i), thick(i)

    end do
    
    close(15)


    end subroutine read_actuatorline_geometry 

    subroutine make_actuatorline_geometry(blade)

    implicit none

    type(ActuatorLineType),intent(INOUT) :: blade ! For simplity I leave it as blade. In fact this is an actuator line 
    integer :: nbe, nei, nej, j
    real :: sEM, tEM, nEM, dx,dy,dz
    real :: PE(3), sE(3), tE(3), normE(3), P1(3), P2(3), P3(3), P4(3), V1(3), V2(3), V3(3), V4(3), A1(3), A2(3)

    ewrite(2,*) 'Entering make_actuatorline_geometry'
    ! Calculates element geometry from element end geometry
    nbe=blade%NElem

    do j=1,nbe
    nej=1+j

    ! Element center locations
    blade%PEx(nej-1)=(blade%QCx(nej)+blade%QCx(nej-1))/2.0
    blade%PEy(nej-1)=(blade%QCy(nej)+blade%QCy(nej-1))/2.0
    blade%PEz(nej-1)=(blade%QCz(nej)+blade%QCz(nej-1))/2.0
    blade%ERdist(nej-1)=sqrt((blade%PEX(nej-1)-blade%COR(1))**2 +(blade%PEY(nej-1)-blade%COR(2))**2+(blade%PEZ(nej-1)-blade%COR(3))**2)    ! Element length

    ! Set spannwise and tangential vectors
    sE=(/blade%QCx(nej)-blade%QCx(nej-1),blade%QCy(nej)-blade%QCy(nej-1),blade%QCz(nej)-blade%QCz(nej-1)/) ! nominal element spanwise direction set opposite to QC line
    sEM=sqrt(dot_product(sE,sE))

    blade%EDS(nej-1) = sEM

    sE=sE/sEM
    tE=(/blade%tx(nej)+blade%tx(nej-1),blade%ty(nej)+blade%ty(nej-1),blade%tz(nej)+blade%tz(nej-1)/)/2.0
    ! Force tE normal to sE
    tE=tE-dot_product(tE,sE)*sE
    tEM=sqrt(dot_product(tE,tE))
    tE=tE/tEM
    blade%sEx(nej-1)=sE(1)
    blade%sEy(nej-1)=sE(2)
    blade%sEz(nej-1)=sE(3)
    blade%tEx(nej-1)=tE(1)
    blade%tEy(nej-1)=tE(2)
    blade%tEz(nej-1)=tE(3)

    ! Calc normal vector
    Call cross(sE(1),sE(2),sE(3),tE(1),tE(2),tE(3),normE(1),normE(2),normE(3))
    nEM=sqrt(dot_product(normE,normE))
    normE=normE/nEM
    blade%nEx(nej-1)=normE(1)
    blade%nEy(nej-1)=normE(2)
    blade%nEz(nej-1)=normE(3)

    if (blade%FlipN) then
    blade%nEx(nej-1)=-blade%nEx(nej-1)
    blade%nEy(nej-1)=-blade%nEy(nej-1)
    blade%nEz(nej-1)=-blade%nEz(nej-1)
    blade%sEx(nej-1)=-blade%sEx(nej-1)
    blade%sEy(nej-1)=-blade%sEy(nej-1)
    blade%sEz(nej-1)=-blade%sEz(nej-1)
    endif
    ! Calc element area and chord
    P1=(/blade%QCx(nej-1)-0.25*blade%C(nej-1)*blade%tx(nej-1),blade%QCy(nej-1)-0.25*blade%C(nej-1)*blade%ty(nej-1),blade%QCz(nej-1)-0.25*blade%C(nej-1)*blade%tz(nej-1)/)

    P2=(/blade%QCx(nej-1)+0.75*blade%C(nej-1)*blade%tx(nej-1),blade%QCy(nej-1)+0.75*blade%C(nej-1)*blade%ty(nej-1),blade%QCz(nej-1)+0.75*blade%C(nej-1)*blade%tz(nej-1)/)

    P3=(/blade%QCx(nej)+0.75*blade%C(nej)*blade%tx(nej),blade%QCy(nej)+0.75*blade%C(nej)*blade%ty(nej),blade%QCz(nej)+0.75*blade%C(nej)*blade%tz(nej)/)

    P4=(/blade%QCx(nej)-0.25*blade%C(nej)*blade%tx(nej),blade%QCy(nej)-0.25*blade%C(nej)*blade%ty(nej),blade%QCz(nej)-0.25*blade%C(nej)*blade%tz(nej)/)

    V1=P2-P1
    V2=P3-P2
    V3=P4-P3
    V4=P1-P4
    ! Calc quad area from two triangular facets
    Call cross(V1(1),V1(2),V1(3),V2(1),V2(2),V2(3),A1(1),A1(2),A1(3))
    A1=A1/2.0
    Call cross(V3(1),V3(2),V3(3),V4(1),V4(2),V4(3),A2(1),A2(2),A2(3))
    A2=A2/2.0
    blade%EArea(nej-1)=sqrt(dot_product(A1,A1))+sqrt(dot_product(A2,A2))
    ! Calc average element chord from area and span
    blade%EC(nej-1)=blade%EArea(nej-1)/sEM
    blade%ETtoC(nej-1)=0.5*(blade%thick(nej)+blade%thick(nej-1))
    blade%Epitch(nej-1)=0.5*(blade%pitch(nej)+blade%pitch(nej-1))
    end do
    ewrite(2,*) 'Exiting make_actuatorline_geometry'

    end subroutine make_actuatorline_geometry 
    
    subroutine allocate_actuatorline(actuatorline,NStations)

    implicit none
    
    type(ActuatorLineType) :: actuatorline
    integer,intent(in) :: Nstations
    integer :: NElem
    
    Nelem=Nstations-1
    actuatorline%Nelem = Nelem

    allocate(actuatorline%QCx(NElem+1))
    allocate(actuatorline%QCy(NElem+1))
    allocate(actuatorline%QCz(NElem+1))
    allocate(actuatorline%tx(NElem+1))
    allocate(actuatorline%ty(NElem+1))
    allocate(actuatorline%tz(NElem+1))
    allocate(actuatorline%C(NElem+1))
    allocate(actuatorline%thick(NElem+1))
    allocate(actuatorline%pitch(NElem+1))
    allocate(actuatorline%PEx(NElem))
    allocate(actuatorline%PEy(NElem))
    allocate(actuatorline%PEz(NElem))
    allocate(actuatorline%tEx(NElem))
    allocate(actuatorline%tEy(NElem))
    allocate(actuatorline%tEz(NElem))
    allocate(actuatorline%nEx(NElem))
    allocate(actuatorline%nEy(NElem))
    allocate(actuatorline%nEz(NElem))
    allocate(actuatorline%sEx(NElem))
    allocate(actuatorline%sEy(NElem))
    allocate(actuatorline%sEz(NElem))
    allocate(actuatorline%EC(NElem))
    allocate(actuatorline%EDS(NElem))
    allocate(actuatorline%EArea(NElem))
    allocate(actuatorline%ETtoC(NElem))
    allocate(actuatorline%Eepsilon(NElem))
    allocate(actuatorline%EAirfoil(Nelem))
    allocate(actuatorline%E_LB_Model(Nelem))
    allocate(actuatorline%ERdist(Nelem))
    allocate(actuatorline%EVx(NElem))
    allocate(actuatorline%EVy(NElem))
    allocate(actuatorline%EVz(NElem))
    allocate(actuatorline%EVbx(NElem))
    allocate(actuatorline%EVby(NElem))
    allocate(actuatorline%EVbz(NElem))
    allocate(actuatorline%EObx(NElem))
    allocate(actuatorline%EOby(NElem))
    allocate(actuatorline%EObz(NElem))
    allocate(actuatorline%Epitch(Nelem))
    allocate(actuatorline%EAOA(Nelem))
    allocate(actuatorline%EAOAdot(Nelem))
    allocate(actuatorline%ERE(Nelem))
    allocate(actuatorline%EUr(Nelem))
    allocate(actuatorline%EUn(Nelem))
    allocate(actuatorline%EAOA_LAST(Nelem))
    allocate(actuatorline%EUn_LAST(Nelem))
    allocate(actuatorline%ECD(Nelem))
    allocate(actuatorline%ECL(Nelem))
    allocate(actuatorline%ECM(Nelem))
    allocate(actuatorline%ECN(Nelem))
    allocate(actuatorline%ECT(Nelem))
    allocate(actuatorline%EFn(NElem))
    allocate(actuatorline%EFt(NElem))
    allocate(actuatorline%EMS(NElem))
    allocate(actuatorline%EFx(NElem))
    allocate(actuatorline%EFy(NElem))
    allocate(actuatorline%EFz(NElem))
    
    end subroutine allocate_actuatorline

end module actuator_line_element
