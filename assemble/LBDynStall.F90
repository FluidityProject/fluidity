!    Copyright (C) 2008 Imperial College London and others.
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
! **************************************************************************************

#include "fdebug.h"

module lbdynstall


    use fldebug
    use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, PYTHON_FUNC_LEN, pi
    implicit none
    
    real, parameter :: conrad = pi / 180.0 
    real, parameter :: condeg = 180.0 / pi  
    
    ! Leishman-Beddoes model parameters
    type LB_Type
    
    logical :: StallFlag = .false.

    real :: dp
    real :: dF
    real :: dCNv
    real :: sLEv
    integer :: LESepState
    real :: CLRef 
    real :: CLRefLE
    real :: CLCritP
    real :: CLCritN
    integer :: CLRateFlag
    real :: Fstat
    real :: F
    real :: cv
    real :: dcv
    real :: CLRef_Last
    real :: CLRefLE_Last
    real :: Fstat_Last
    real :: cv_Last

    ! Additional LB diagnostic output
    integer, dimension(9) :: LB_LogicOutputs
    integer :: Logic_W(9)=[1,1,1,1,1,3,2,1,1]
    
    end type LB_Type

    Type(LB_Type) :: lb

    contains

    subroutine dystl_init_LB(lb)
        
        implicit none
        type(LB_Type) :: lb
       
        lb%StallFlag = .true.
        lb%dp=0.0
        lb%dF=0.0
        lb%dCNv=0.0
        lb%LESepState=0.0
        lb%sLEv=0.0
        lb%CLRef_Last=0.0
        lb%CLRefLE_Last=0.0
        lb%Fstat_Last=1.0
        lb%cv_Last=0.0
   
        lb%LB_LogicOutputs(:)=0

    end subroutine dystl_init_LB

    subroutine LB_EvalIdealCL(AOA,AOA0,CLa,RefFlag,CLID)

    ! AOA inputs in radians
    ! AOA0 is zero lift AOA
    ! RefFlag defines whether to output reference CL or ideal CL
    ! CLa is reference lift slope (per radian) to be used for reference CL (ideal CLa is 2*pi)

    real :: AOA, AOA0, CLa
    integer :: RefFlag
    real :: CLID
    real :: IDS, aID, d1, CLaI, aIDc

    aID=AOA-AOA0
    call Force180(aID)
    ! reflect function across axis for abs(aID)>90
    IDS=1
    if (aID>pi/2.0) then
        aID=pi-aID
        IDS=-1.0
    else if (aID<-pi/2.0) then
        aID=-pi-aID
        IDS=-1.0
    end if

    ! If RefFlag is 1, output reference CL, otherwise round off ideal CL at high AOA
    if (RefFlag==1) then
        CLID=IDS*CLa*aID
    else
        ! round off the ideal CL after cutoff AOA
        aIDc=30.0*conrad
        d1=1.8
        CLaI=2.0*pi
        if (abs(aID)<aIDc) then
            CLID=IDS*CLaI*aID
        else
            CLID=IDS*(CLaI*(aIDc-1.0/d1*sin(d1*aIDc))+CLaI/d1*sin(d1*aID))
        end if
    end if
    
    end subroutine LB_EvalIdealCL

    subroutine Force180(a)
        real :: a
        ! alpha in radians
        if (a>pi) then
            a=a-2.0*pi
        elseif (a<-pi) then
            a=a+2.0*pi
        endif
    end subroutine Force180
    
    SUBROUTINE LB_UpdateStates(lb)

        ! Update states for the LB model 
        ! Note dynstall should be included eventually in an expanded blade module, at which point it would have 
        ! access to the geometry info it needs...
        type(LB_Type) :: lb
        integer :: i, nei, j, IsBE

        real :: ds,Tf,TfRef,Tp,TvRef, Tv, TvL

        ! Set model parameters. All of these are potentially a function of Mach
        ! and are set to low mach values...
        Tp=1.7                  ! time constant on LE pressure response to change in CL
        TfRef=3.0               ! time constant on TE separation point travel
        TvRef=6.0               ! time constant on LE vortex lift indicial function
        TvL=11.0                ! Characteristic LE vortex travel time



    ! Eval LE separation state
    ! Note: CLCrit is critical (ideal) CL value for LE separation. This is
    ! approximately equal to the CL that would exist at the angle of attack at
    ! max CL if the CL curve had remained linear
    if (lb%LESepState==0 .AND. (lb%CLRefLE>lb%CLCritP .OR. lb%CLRefLE<lb%CLCritN)) then
        ! In LE separation state
        lb%LESepState=1
        lb%sLEv=0 ! reset leading edge vortex time counter

        ! Set logic state flags (for model diagnosis output)
        lb%LB_LogicOutputs(5)=1
    else if (lb%LESepState==1 .AND. (lb%CLRefLE<lb%CLCritP .AND. lb%CLRefLE>lb%CLCritN)) then
        ! Out of LE separation state
        lb%LESepState=0
        lb%sLEv=0 ! reset leading edge vortex time counter

        ! Set logic state flags (for model diagnosis output)
        lb%LB_LogicOutputs(5)=2
    end if

    ! Set time constants based on LE separation state and TE separation point
    ! location. Different time constants for abs(CL) increasing vs. decreasing
    if (lb%LESepState==1) then
        if (lb%sLEv<TvL) then
            if (lb%CLRateFlag>0) then
                !Tf=3.0*TfRef ! original
                Tf=4.0*TfRef
                Tv=TvRef

                ! Set logic state flags (for model diagnosis output)
                lb%LB_LogicOutputs(8)=1
            else
                Tf=1.0/2.0*TfRef
                Tv=1.0/2.0*TvRef

                ! Set logic state flags (for model diagnosis output)
                lb%LB_LogicOutputs(8)=2
            end if

            ! Set logic state flags (for model diagnosis output)
            lb%LB_LogicOutputs(7)=1
        else if (lb%sLEv<2.0*TvL) then
            if (lb%CLRateFlag>0) then
                ! orig 
                !Tf=1.0/3.0*TfRef
                !Tv=1.0/4.0*TvRef
                Tf=2.0*TfRef
                Tv=TvRef

                ! Set logic state flags (for model diagnosis output)
                lb%LB_LogicOutputs(8)=3
            else
                Tf=1.0/2.0*TfRef
                Tv=1.0/2.0*TvRef

                ! Set logic state flags (for model diagnosis output)
                lb%LB_LogicOutputs(8)=4                                                        
            end if

            ! Set logic state flags (for model diagnosis output)
            lb%LB_LogicOutputs(7)=2                                                
        else
            ! orig
            !Tf=4.0*TfRef
            !Tv=0.9*TvRef
            Tf=TfRef
            Tv=TvRef

            ! Set logic state flags (for model diagnosis output)
            lb%LB_LogicOutputs(7)=3
        end if

        ! Set logic state flags (for model diagnosis output)
        lb%LB_LogicOutputs(6)=1
    else
        if (lb%F>0.7) then
            Tf=TfRef

            ! Set logic state flags (for model diagnosis output)
            lb%LB_LogicOutputs(7)=4
        else
            Tf=2*TfRef

            ! Set logic state flags (for model diagnosis output)
            lb%LB_LogicOutputs(7)=5
        end if
        Tv=TvRef

        ! Set logic state flags (for model diagnosis output)
        lb%LB_LogicOutputs(6)=2
    end if

    ! update LE vortex time counter if in LE separation state
    if (lb%LESepState==1) then
        lb%sLEv=lb%sLEv+ds

        ! Set logic state flags (for model diagnosis output)
        lb%LB_LogicOutputs(9)=1
    end if

    ! Update states, first order lag equations, exponential recursion form (midpoint rule version)
     lb%dp=lb%dp*exp(-ds/Tp)+(lb%CLRef-lb%CLRef_Last)*exp(-ds/(2*Tp))
     lb%dF=lb%dF*exp(-ds/Tf)+(lb%Fstat-lb%Fstat_Last)*exp(-ds/(2*Tf))
     lb%dCNv=lb%dCNv*exp(-ds/Tv)+lb%dcv*exp(-ds/(2*Tv))

    ! update lagged values
    lb%CLRef_Last=lb%CLRef
    lb%CLRefLE_Last=lb%CLRefLE
    lb%Fstat_Last=lb%Fstat 
    lb%cv_Last=lb%cv

    End SUBROUTINE LB_UpdateStates

    SUBROUTINE LB_LogicChecksum(LBCheck)
    
        implicit none

        integer :: LBCheck
        integer :: Loop

        ! Calculates a checksum for the logic states in the LB model

        LBCheck=0
        do Loop=1,9
            LBCheck=LBCheck+lb%LB_LogicOutputs(Loop)*lb%Logic_W(Loop)
        end do

    End SUBROUTINE LB_LogicChecksum


end module
