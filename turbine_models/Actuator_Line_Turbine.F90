#include "fdebug.h"

module actuator_line_turbine

    use fldebug
    use spud
    use futils
    use global_parameters, only: pi
    
    use Airfoils
    use actuator_line_element

    implicit none

type TurbineType 
    character(len=100) :: name
    character(len=100) :: blade_geom_file
    character(len=100) :: type
    integer :: ID
    integer :: NBlades
    real, dimension(3) :: RotN, origin ! Rotational vectors in the normal and perpendicular directions
    real :: hub_tilt_angle, blade_cone_angle, yaw_angle 
    real :: Rmax ! Reference radius, velocity, viscosity
    real :: A   ! Rotor area
    real :: angularVel,TSR,Uref
    real :: AzimAngle=0.0
    integer :: No_rev=0.0
    logical :: Is_constant_rotation_operated = .false. ! For a constant rotational velocity (in Revolutions Per Minute)
    logical :: Is_force_based_operated = .false. ! For a forced based rotational velocity (Computed during the simulation)
    logical :: IsClockwise = .false.
    logical :: IsCounterClockwise = .false. 
    logical :: Has_Tower=.false.
    real :: Towerheight, TowerOffset, TowerDrag, TowerLift, TowerStrouhal
    logical :: do_tip_correction=.false.

    type(ActuatorLineType), allocatable :: Blade(:)
    type(AirfoilType), allocatable :: AirfoilData(:)
    type(ActuatorLineType) :: tower
    real :: CP  ! Power coefficient 
    real :: CTR ! Torque coefficient 
    real :: CFx ! Fx coefficient 
    real :: CFy ! Fy coefficient 
    real :: CFz ! Fz coefficient 
    real :: CT  ! Thrust coefficient    
    
end type TurbineType
    

contains
    
    subroutine set_turbine_geometry(turbine)

    implicit none
    type(TurbineType),intent(inout) :: turbine
    real, allocatable :: rR(:),ctoR(:),pitch(:),thick(:)
    real :: SVec(3), theta, origin(3)
    integer :: Nstations, iblade, Istation,ielem

    ewrite(2,*) 'Entering set_turbine_geometry'

    ewrite(2,*) 'Turbine Name : ', turbine%name 
    ewrite(2,*) '============================='
    ewrite(2,*) 'Number of Blades : ', turbine%Nblades
    ewrite(2,*) 'Origin           : ', turbine%origin
    ewrite(2,*) 'Hub Tilt Angle   : ', turbine%hub_tilt_angle
    ewrite(2,*) 'Blade Cone Angle : ', turbine%blade_cone_angle
    ewrite(2,*) 'Yaw Angle        : ', turbine%yaw_angle

    call read_actuatorline_geometry(turbine%blade_geom_file,turbine%Rmax,SVec,rR,ctoR,pitch,thick,Nstations)
    ! Make sure that the spanwise is [0 0 1]
    Svec = (/0.0,0.0,1.0/)
    ! Make sure that origin is [0,0,0] : we set everything to origin 0 and then translate the
    ! turbine to the actual origin(this is for simplicity)
    theta=2*pi/turbine%Nblades
    do iblade=1,turbine%Nblades
    call allocate_actuatorline(Turbine%blade(iblade),Nstations)
    turbine%blade(iblade)%name=trim(turbine%name)//'_blade'//int2str(iblade)
    
    turbine%blade(iblade)%COR=turbine%origin
    turbine%blade(iblade)%L=turbine%Rmax
    turbine%blade(iblade)%NElem=Nstations-1 
    
    do istation=1,Nstations
    turbine%blade(iblade)%QCx(istation)=rR(istation)*turbine%Rmax*Svec(1)+turbine%blade(iblade)%COR(1)
    turbine%blade(iblade)%QCy(istation)=rR(istation)*turbine%Rmax*Svec(2)+turbine%blade(iblade)%COR(2)
    turbine%blade(iblade)%QCz(istation)=rR(istation)*turbine%Rmax*Svec(3)+turbine%blade(iblade)%COR(3)
    if(turbine%IsCounterClockwise) then
        turbine%RotN=(/-1.0,0.0,0.0/)
        turbine%blade(iblade)%tx(istation)=sin(pitch(istation)/180.0*pi)    
        turbine%blade(iblade)%ty(istation)=-cos(pitch(istation)/180.0*pi)    
        turbine%blade(iblade)%tz(istation)= 0.0
        turbine%blade(iblade)%C(istation)=ctoR(istation)*turbine%Rmax
        turbine%blade(iblade)%thick(istation)=thick(istation)
        turbine%blade(iblade)%pitch(istation)=pitch(istation)/180.0*pi
    elseif(turbine%IsClockwise) then
        turbine%RotN=(/1.0,0.0,0.0/)
        turbine%blade(iblade)%tx(istation)=sin(pitch(istation)/180.0*pi)    
        turbine%blade(iblade)%ty(istation)=cos(pitch(istation)/180.0*pi)    
        turbine%blade(iblade)%tz(istation)= 0.0
        turbine%blade(iblade)%C(istation)=ctoR(istation)*turbine%Rmax
        turbine%blade(iblade)%thick(istation)=thick(istation)
        turbine%blade(iblade)%pitch(istation)=pitch(istation)/180.0*pi
        turbine%blade(iblade)%FlipN = .true.
    endif
    end do

    ! Always rotate counterclockwise to assign the turbine blades
    call rotate_actuatorline(turbine%blade(iblade),turbine%blade(iblade)%COR,(/-1.0,0.0,0.0/),(iblade-1)*theta)   
    ! Rotate through incidence (hub tilt) and coning angle
    call make_actuatorline_geometry(turbine%blade(iblade))
    ! Populate element Airfoils 
    call populate_blade_airfoils(turbine%blade(iblade)%NElem,turbine%Blade(iblade)%NAirfoilData,turbine%Blade(iblade)%EAirfoil,turbine%Blade(iblade)%AirfoilData,turbine%Blade(iblade)%ETtoC)
    
    turbine%Blade(iblade)%EAOA_LAST(:)=-666
    turbine%Blade(iblade)%EUn_LAST(:)=0.0
    
    ! Initialize LB_model and Tip Correction Coeffs
    do ielem=1,turbine%blade(iblade)%Nelem
    if(turbine%blade(iblade)%do_dynamic_stall) then
    call dystl_init_LB(turbine%blade(iblade)%E_LB_Model(ielem))
    endif
    end do
    
    end do

    !=========================================================
    ! Create a Tower
    !=========================================================
    if(turbine%has_Tower) then
    call read_actuatorline_geometry(turbine%tower%geom_file,turbine%Towerheight,SVec,rR,ctoR,pitch,thick,Nstations)
    ! Make sure that the spanwise is [0 0 1]
    Svec = (/0.0,0.0,1.0/)
    
    call allocate_actuatorline(Turbine%Tower,Nstations)
    turbine%tower%name=trim(turbine%name)//'_tower'
   
    turbine%Tower%COR=turbine%origin
    turbine%Tower%NElem=Nstations-1  
    
    do istation=1,Nstations
    turbine%Tower%QCx(istation)= turbine%Tower%COR(1) + turbine%TowerOffset  
    turbine%Tower%QCy(istation)= turbine%Tower%COR(2)
    turbine%Tower%QCz(istation)= turbine%Tower%COR(3) - rR(istation)*turbine%Towerheight*Svec(3)
    turbine%Tower%tx(istation)= 1.0    
    turbine%Tower%ty(istation)= 0.0    
    turbine%Tower%tz(istation)= 0.0
    turbine%Tower%C(istation)=ctoR(istation)*turbine%Towerheight
    turbine%Tower%thick(istation)=thick(istation)
    turbine%Tower%pitch(istation)=pitch(istation)/180.0*pi
    enddo
    
    call make_actuatorline_geometry(turbine%tower)
    
    turbine%tower%EAOA_LAST(:)=-666
    turbine%tower%EUn_LAST(:)=0.0
    
    !Set the tower body velocity to zero
    turbine%tower%EVbx(:)=0.0
    turbine%tower%EVby(:)=0.0
    turbine%tower%EVbz(:)=0.0
    turbine%tower%EObx(:)=0.0
    turbine%tower%EOby(:)=0.0
    turbine%tower%EObz(:)=0.0

    endif
     
    !========================================================
    !Compute a number of global parameters for the turbine
    !========================================================
    turbine%Uref=turbine%angularVel*turbine%Rmax/turbine%TSR
    turbine%A=pi*turbine%Rmax**2
    
    call Compute_Turbine_RotVel(turbine)
    
    ewrite(2,*) 'Exiting set_turbine_geometry'

    end subroutine set_turbine_geometry
    
    subroutine compute_performance(turbine)

    implicit none
    type(TurbineType), intent(inout) :: turbine
    real :: Torque_i,FX_i,FY_i,FZ_i,fx_tot,fy_tot,fz_tot,torq_tot 
    real :: xe,ye,ze,o1,o2,o3,fx,fy,fz,trx,try,trz,te,ms,sxe,sye,sze
    real :: rotx,roty,rotz
    real :: U_ref, R, A
    integer :: iblade, ielem
    ewrite(2,*) 'In calculate_performance'

    RotX=turbine%RotN(1)
    RotY=turbine%RotN(2)
    RotZ=turbine%RotN(3)

    ! Compute Torque for each Blade
        Fx_tot=0.
        Fy_tot=0.
        Fz_tot=0.
        Torq_tot=0

        do iblade=1,turbine%Nblades
        Fx_i=0.
        Fy_i=0.
        Fz_i=0.
        Torque_i=0.

        do ielem=1,turbine%blade(iblade)%NElem
            
            xe=turbine%blade(iblade)%PEX(ielem)
            ye=turbine%blade(iblade)%PEY(ielem)
            ze=turbine%blade(iblade)%PEZ(ielem)
            o1=turbine%blade(iblade)%COR(1)
            o2=turbine%blade(iblade)%COR(2)
            o3=turbine%blade(iblade)%COR(3)
            fx=turbine%blade(iblade)%EFX(ielem)
            fy=turbine%blade(iblade)%EFY(ielem)
            fz=turbine%blade(iblade)%EFZ(ielem)
            ms=turbine%blade(iblade)%EMS(ielem)
            sxe=turbine%blade(iblade)%sex(ielem)
            sye=turbine%blade(iblade)%sey(ielem)
            sze=turbine%blade(iblade)%sez(ielem)

            call cross(xe-o1,ye-o2,ze-o3,fx,fy,fz,trx,try,trz)
            te=(trx*RotX+try*RotY+trz*RotZ)+ms*(sxe*RotX+sye*RotY+sze*RotZ)
            
            Fx_i=Fx_i+fx
            Fy_i=Fy_i+fy
            Fz_i=Fz_i+fz
            Torque_i=Torque_i+te
        end do
            
            Fx_tot=Fx_tot+Fx_i
            Fy_tot=Fy_tot+Fy_i
            Fz_tot=Fz_tot+Fz_i
            Torq_tot=torq_tot+Torque_i
    end do

    turbine%CFx=FX_tot/(0.5*turbine%A*turbine%Uref**2)
    turbine%CFy=FY_tot/(0.5*turbine%A*turbine%Uref**2)
    turbine%CFz=Fz_tot/(0.5*turbine%A*turbine%Uref**2)
    turbine%CT=sqrt(turbine%CFx**2.0+turbine%CFy**2.0+turbine%CFz**2.0)
    turbine%CTR=Torq_tot/(0.5*turbine%A*turbine%Rmax*turbine%Uref**2.0)
    turbine%CP= abs(turbine%CTR)*turbine%TSR
    
    ewrite(2,*) '--------------------------------------------------------'
    ewrite(2,*) 'Calculate performance for Turbine : ',turbine%name
    ewrite(2,*) 'Azimuthal Angle (degrees) : ', turbine%AzimAngle
    ewrite(2,*) 'Thrust Coefficient : ', turbine%CT
    ewrite(2,*) 'Torque Coefficient : ', turbine%CTR
    ewrite(2,*) 'Power Coefficient : ', turbine%CP
    ewrite(2,*) '--------------------------------------------------------'
   
    ewrite(2,*) 'Exiting compute_performance'

    end subroutine compute_performance
    
    subroutine Compute_Turbine_Tip_Correction(turbine)
    
    implicit none
    type(TurbineType),intent(inout) :: turbine
    integer :: iblade,ielem
    real ::g1,alpha,pitch,F,r
        
    g1=exp(-0.125*(turbine%NBlades*turbine%tsr-21.0))+0.1
    
    do iblade=1,turbine%Nblades
    ! Compute angle phi at the tip 
        
    do ielem=1,turbine%blade(iblade)%Nelem
       
        r=turbine%blade(iblade)%ERdist(ielem)
        alpha=turbine%blade(iblade)%EAOA_Last(ielem)
        pitch=turbine%blade(iblade)%Epitch(ielem)

        F=2.0/pi*acos(exp(-g1*turbine%Nblades*(turbine%Rmax-turbine%blade(iblade)%ERdist(ielem))/(2.0*r*sin(alpha+pitch))))
        
        ! Apply Coeffs to the local Element forces
        turbine%blade(iblade)%EFN(ielem)=F*turbine%blade(iblade)%EFN(ielem)
        turbine%blade(iblade)%EFT(ielem)=F*turbine%blade(iblade)%EFT(ielem)
        ! Apply to the global forces 
        turbine%blade(iblade)%EFX(ielem)=F*turbine%blade(iblade)%EFX(ielem)
        turbine%blade(iblade)%EFY(ielem)=F*turbine%blade(iblade)%EFY(ielem)
        turbine%blade(iblade)%EFZ(ielem)=F*turbine%blade(iblade)%EFZ(ielem)

        end do
    end do

    end subroutine Compute_Turbine_Tip_Correction

    subroutine Compute_Turbine_RotVel(turbine)
    implicit none
    type(TurbineType),intent(inout) :: turbine
    integer :: iblade,ielem
    real :: wRotX,wRotY,wRotZ,Rx,Ry,Rz,ublade,vblade,wblade
    real :: RotX,RotY,RotZ 
    ewrite(2,*) 'Entering Compute_Turbine_Local_Vel '
    
    !========================================================
    ! Compute Element local rotational velocity
    !========================================================
    wRotX=turbine%angularVel*turbine%RotN(1)
    wRotY=turbine%angularVel*turbine%RotN(2)
    wRotZ=turbine%angularVel*turbine%RotN(3)
    
    do iblade=1,turbine%NBlades
    do ielem=1,turbine%Blade(iblade)%Nelem
     
    Rx=-turbine%blade(iblade)%COR(1)+turbine%Blade(iblade)%PEx(ielem);
    Ry=-turbine%blade(iblade)%COR(2)+turbine%Blade(iblade)%PEy(ielem);
    Rz=-turbine%blade(iblade)%COR(3)+turbine%Blade(iblade)%PEz(ielem);

   ! ! Find the cross product Ublade = Omega x R
    call cross(wRotX,wRotY,wRotZ,Rx,Ry,Rz,ublade,vblade,wblade)
    
    turbine%Blade(iblade)%EVbx(ielem)=ublade
    turbine%Blade(iblade)%EVby(ielem)=vblade
    turbine%Blade(iblade)%EVbz(ielem)=wblade
    turbine%Blade(iblade)%EObx(ielem)=wRotX
    turbine%Blade(iblade)%EOby(ielem)=wRotY
    turbine%Blade(iblade)%EObz(ielem)=wRotZ
    
    end do
    end do
    
    ewrite(2,*) 'Exiting Compute_Turbine_Local_Vel '

    end subroutine Compute_Turbine_RotVel

    subroutine rotate_turbine(turbine,theta)
    
    implicit none
    type(TurbineType),intent(inout) :: turbine
    real,intent(in) :: theta
    real :: nrx,nry,nrz,px,py,pz 
    integer :: j,ielem,i
    real :: vrx,vry,vrz,VMag
    real :: xtmp,ytmp,ztmp, txtmp, tytmp, tztmp
    ! Rotates data in blade arrays. Rotate element end geometry and recalculate element geometry.

    ewrite(1,*) 'Entering rotate_turbines'
        
    nrx=turbine%RotN(1)
    nry=turbine%RotN(2)
    nrz=turbine%RotN(3) 
    px =turbine%origin(1)
    py =turbine%origin(2)
    pz =turbine%origin(3)


    do j=1,turbine%NBlades
        do ielem=1,turbine%Blade(j)%Nelem+1
        ! Blade end locations (quarter chord). xBE(MaxSegEnds)
        xtmp=turbine%Blade(j)%QCx(ielem)
        ytmp=turbine%Blade(J)%QCy(ielem)
        ztmp=turbine%Blade(j)%QCz(ielem)
        
        Call QuatRot(xtmp,ytmp,ztmp,theta,nrx,nry,nrz,px,py,pz,vrx,vry,vrz)
        turbine%Blade(j)%QCx(ielem)=vrx                                       
        turbine%Blade(j)%QCy(ielem)=vry                                       
        turbine%Blade(j)%QCz(ielem)=vrz                                  
        
        txtmp=turbine%Blade(j)%tx(ielem)
        tytmp=turbine%Blade(j)%ty(ielem)
        tztmp=turbine%Blade(j)%tz(ielem)
        
        ! Tangent vectors
        Call QuatRot(txtmp,tytmp,tztmp,theta,nrx,nry,nrz,px,py,pz,vrx,vry,vrz)
        VMag=sqrt(vrx**2+vry**2+vrz**2)
        turbine%Blade(j)%tx(ielem)=vrx/VMag                                      
        turbine%Blade(j)%ty(ielem)=vry/VMag                                  
        turbine%Blade(j)%tz(ielem)=vrz/VMag                                       
  
        end do
        
        call make_actuatorline_geometry(turbine%Blade(j))
    end do 
    
    ewrite(1,*) 'Exiting rotate_turbine'
    end subroutine rotate_turbine  
    
end module actuator_line_turbine 
