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

module alturbine_airfoils 

  use fldebug
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, PYTHON_FUNC_LEN, pi

  implicit none
  
  type AirfoilType
  
  ! Naming and others
  character(len=100) :: afname, aftitle ! Title for each airfoil section 
  integer :: camb ! Camber flag for each section
  real    :: tc      ! Thickness to chord ration for each section
  real    :: alzer   ! Zero lift AOA for each section
  
  ! Airfoild section coefficient data
  real, allocatable :: TA(:,:)   ! Table AOA values
  real, allocatable :: TCL(:,:)  ! Table CL values
  real, allocatable :: TCD(:,:)  ! Table CD values
  real, allocatable :: TCM(:,:)  ! Table Cm values
  real, allocatable :: TRE(:)    ! Table Reynolds Number values  
  integer, allocatable :: nTBL(:)   ! Number of AOA values for each Re number, in each section data table
  integer  :: nRET   ! Number of Re number values in each section data table

  ! Interpolation warning flags
  integer :: ilxtp
  integer :: iuxtp

  ! Airfoil params for BV dyn stall
  real, allocatable :: alstlp(:)    ! Stall AOA (positive) at all Re numbers
  real, allocatable :: alstln(:)    ! Stall AOA (negative) at all Re numbers

  ! Airfoil params for Leisham Beddoes dyn stall
  real, allocatable :: CLaData(:)
  real, allocatable :: CLCritPData(:)
  real, allocatable :: CLCritNData(:)
 
  end type AirfoilType
 
  ! Maximum Numbers of Reynolds Number data that can be stored
  ! Globla parameters
  integer, parameter :: MaxReVals = 20
  integer, parameter :: MaxAOAVals = 1000

  real, parameter :: conrad = pi / 180.0 
  real, parameter :: condeg = 180.0 / pi  
  ! Private subroutines
  private intp, read_airfoil, allocate_airfoil
 
  ! Public subroutines
  public airfoil_init, compute_aeroCoeffs, compute_aeroCoeffs_one_airfoil

contains
   
    subroutine airfoil_init(airfoil)
    
    implicit none
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    ! This routine initialises the airfoil struct by reading 
    ! the existing input files and allocating the memory of the
    ! airfoil array
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    type(AirfoilType),intent(INOUT) :: airfoil

    ewrite(1,*) 'In airfoils_init'
    
    call allocate_airfoil(airfoil,MaxAOAVals,MaxReVals)
 
    call read_airfoil(airfoil)
    ewrite(2,*) 'Airfoil section title : ', airfoil%aftitle
    ewrite(1,*) 'Exiting airfoils_init'

    end subroutine airfoil_init
    
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    ! These routines have been taken directly out of CACTUS and changed
    ! slightly so they fit the needs of the present coupling
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

    subroutine read_airfoil(airfoil)
    
    implicit none
    
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    ! This routine reads the airfoil data from 
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    
    type(AirfoilType),intent(INOUT) :: airfoil
    character(len=100) :: ReadLine
    logical :: NotDone, NotBlank
    integer :: EOF, CI, i, ii, jj
    real :: temp, temp1(1000,4)


    ewrite(1,*) 'In read_airfoil'

    open(15, file = airfoil%afname)
    EOF=0

    ! Find title block
    NotDone=.TRUE.
    do while (NotDone)
        read(15,'(A)') ReadLine
        CI=index(ReadLine,':')
        if (CI>0) then
            NotDone=.FALSE.
        end if
    end do

    ! Read title and airfoil thickness
    if(len_trim(ReadLine)>CI) then
        airfoil%aftitle = ReadLine(CI+1:len_trim(ReadLine))
    else
        airfoil%aftitle = 'No Title'
    end if
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) airfoil%tc
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) airfoil%alzer
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) airfoil%camb
 
    ! Reverse camber direction if desired
    if (airfoil%camb==1) then
        airfoil%alzer = -airfoil%alzer
    end if

    ! Find first Reynolds Number block
    NotDone =.true.
    do while (NotDone)
        read(15,'(A)',IOSTAT=EOF) ReadLine
        CI=index(ReadLine,':')
        if (CI>0 .OR. EOF<0) then
            NotDone=.false.
        end if
    end do

    ! Read data for each Reynolds value
     
    do while (EOF>=0 .and. (i<MaxReVals))
        i=i+1
        ! Read Re and dyn. stall data
        read(ReadLine(index(ReadLine,':')+1:),*) airfoil%TRE(i)
        read(15,'(A)') ReadLine                          
        read(ReadLine(index(ReadLine,':')+1:),*) airfoil%alstlp(i)
        airfoil%alstlp(i)=airfoil%alstlp(i)*conrad
        read(15,'(A)') ReadLine                          
        read(ReadLine(index(ReadLine,':')+1:),*) airfoil%alstln(i)
        airfoil%alstln(i)=airfoil%alstln(i)*conrad
        read(15,'(A)') ReadLine                          
        read(ReadLine(index(ReadLine,':')+1:),*) airfoil%CLaData(i)
        read(15,'(A)') ReadLine                          
        read(ReadLine(index(ReadLine,':')+1:),*) airfoil%CLCritPData(i)
        read(15,'(A)') ReadLine                          
        read(ReadLine(index(ReadLine,':')+1:),*) airfoil%CLCritNData(i)

        ! Reverse camber direction if desired
        if (airfoil%camb == 1) then
            temp = airfoil%alstlp(i)
            airfoil%alstlp(i) = -airfoil%alstln(i)
            airfoil%alstln(i) = -temp   
            temp = airfoil%CLCritPData(i)
            airfoil%CLCritPData(i) = -airfoil%CLCritNData(i)
            airfoil%CLCritNData(i) = -temp 
        end if

        ! Read AOA data
        read(15,'(A)') ReadLine
        NotDone=.TRUE.
        ii=0
        do while (NotDone)
            read(15,'(A)',IOSTAT=EOF) ReadLine
            ! Check for carriage return (len_trim doesn't consider this a blank)
            NotBlank=.TRUE.
            if (len_trim(ReadLine)==0) then
                NotBlank=.FALSE.
            else if (len_trim(ReadLine)==1) then
                if (ichar(ReadLine(len_trim(ReadLine):len_trim(ReadLine))) == 13) then
                    NotBlank=.FALSE.
                end if
            end if
            if (EOF>=0 .AND. NotBlank) then
                if (ii == MaxAOAVals) then
                    FLExit("Max. allowed AOA values exceeded in airfoil data file: "//airfoil%aftitle)
                    NotDone=.FALSE.
                else
                    ii=ii+1                        
                    read(ReadLine,*) airfoil%ta(ii,i),airfoil%tcl(ii,i),airfoil%tcd(ii,i),airfoil%tcm(ii,i) 
                end if
            else
                NotDone=.FALSE.
            end if
        end do
        airfoil%ntbl(i)=ii

        !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGg
        ! This is under consideration. In general it is very difficult 
        ! to find data for airfoils spanning from -180 to 180
        !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

         ! Check AOA limits
        if (airfoil%ta(1,i) > -180.0 .OR. airfoil%ta(airfoil%ntbl(i),i) < 180.0) then
            FLExit("AOA data needs to be +/-180 deg in airfoil data file: "// airfoil%aftitle)
        end if

        ! Reverse camber direction if desired
        if(airfoil%camb == 1) then       
            do ii = 1, airfoil%ntbl(i)
                temp1(ii,1) = airfoil%ta(ii,i) 
                temp1(ii,2) = airfoil%tcl(ii,i)
                temp1(ii,3) = airfoil%tcd(ii,i)
                temp1(ii,4) = airfoil%tcm(ii,i)
            end do

            do ii = 1, airfoil%ntbl(i)
                jj = airfoil%ntbl(i)-(ii-1)
                airfoil%ta(ii,i) = -temp1(jj,1)
                airfoil%tcl(ii,i) = -temp1(jj,2)
                airfoil%tcd(ii,i) = temp1(jj,3)
                airfoil%tcm(ii,i) = -temp1(jj,4)
            end do
        end if

        ! Find next Re block
        NotDone=.TRUE.
        if (EOF<0) then
            NotDone=.FALSE.
        end if
        do while (NotDone)
            read(15,'(A)',IOSTAT=EOF) ReadLine
            CI=index(ReadLine,':')
            if (CI>0 .OR. EOF<0) then
                NotDone=.FALSE.
            end if
        end do

        end do
        ! Set number of Re vals for this section
        airfoil%nRET=i

        ! Close input file for this section
        close(15)

        ! Check data
        if (i == 0) then
            FLExit("Error reading airfoil data file: "// airfoil%aftitle)
        end if
        if (EOF > 0) then
            FLExit("Warning: Max. allowed Re values exceeded in airfoil data file: "//airfoil%aftitle)
        end if

         
    ewrite(1,*) 'Exiting read_airfoil'

    end subroutine read_airfoil
   
    subroutine allocate_airfoil(airfoil,MaxAOAVals,MaxReVals)
    
    implicit none
    
    type(AirfoilType),intent(INOUT) :: airfoil
    integer, intent(IN) :: MaxAOAVals,MaxReVals

    ewrite(1,*) 'In allocate_airfoil'
    allocate(airfoil%TA(MaxAOAVals,MaxReVals))
    allocate(airfoil%TCL(MaxAOAVals,MaxReVals))
    allocate(airfoil%TCD(MaxAOAVals,MaxReVals))
    allocate(airfoil%TCM(MaxAOAVals,MaxReVals))
    allocate(airfoil%TRE(MAxReVals))
    allocate(airfoil%nTBL(MaxReVals))
    allocate(airfoil%alstln(MaxReVals))
    allocate(airfoil%alstlp(MaxReVals))
    allocate(airfoil%CLaData(MaxReVals))
    allocate(airfoil%CLCritPData(MaxReVals))
    allocate(airfoil%CLCritNData(MaxReVals))

    ewrite(1,*) 'Exiting allocate_airfoil'

    end subroutine allocate_airfoil
   
    subroutine compute_aeroCoeffs(airfoils,thtoc,alpha75,alpha5,Re,adotnorm,umach,CL,CD,CN,CT,CLCirc,CM25)

    ! GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    ! inputs :
    !           NumAir  : Number of available airfoils
    !           Airfoils: structure for all the airfoil sections available  
    !           thtoc   : Thickness to Chord length 
    !           alpha75 : Angle of Attack at 3/4 of the element 
    !           alpha5  : Angle of Attack at 1/2 (middle) of the element
    !           Re      : Element Reynolds Number
    !           adotnorm: rate of change of the angle of attack (locally)
    !           umach   : local element Mach Number 
    ! 
    ! outputs: 
    !           CL      : Lift Coefficient
    !           CD      : Drag Coefficient
    !           CN      : Normal Force Coefficient
    !           CD      : Tangential Force Coefficient
    !           CLCirc  :
    !               
    !       
    ! GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    implicit none
    type(AirfoilType),dimension(:),target, intent(IN) :: airfoils 
    type(AirfoilType), pointer :: airfoil_pointer => null()
    real, intent(IN) :: alpha75, alpha5, adotnorm, Re, umach 
    real :: alpha75_loc, alpha5_loc, adotnorm_loc, Re_loc, umach_loc
    real,intent(OUT) ::CL, CD, CN, CT, CLCirc, thtoc, CM25
    real,allocatable :: Thicks(:), diffthick(:)
    real :: xtc, CLtc(2),CDtc(2),CM25tc(2), CNtc(2),CTtc(2),CLcirctc(2)                                 
    integer :: NumAir , abcdf
    integer :: iair,imin,imax, iL, iU, iint
    logical :: notDone=.true.
    
    CLtc(:)=0.0
    CDtc(:)=0.0
    CM25tc(:)=0.0
    CNtc(:)=0.0
    CTtc(:)=0.0
    CLCirctc(:)=0.0

    NumAir = size(airfoils)
    allocate(Thicks(NumAir),diffthick(NumAir))

    if(NumAir==1) then
        call compute_aeroCoeffs_one_airfoil(airfoils(1),alpha75,alpha5,Re,adotnorm,umach,CL,CD,CN,CT,CLCirc,CM25)
    else
        !> Linearly Interpolating from the available airfoils 
        !> If the airfoil thickness to chord lenght ratio is exceeded then the min or max airfoils will be used
        do iair=1,NumAir
            Thicks(iair)=Airfoils(iair)%tc
            diffthick(iair)=abs(thtoc-Thicks(iair))
        end do
        
        ! Minimum and Maximum indices 
        imin=minloc(Thicks,1)
        imax=maxloc(Thicks,1)
        iint=minloc(diffthick,1)
        write(*,*) iint
        if (thtoc<=minval(Thicks)) then
            iair=minloc(Thicks,1)
            call compute_aeroCoeffs_one_airfoil(airfoils(iair),alpha75,alpha5,Re,adotnorm,umach,CL,CD,CN,CT,CLCirc,CM25)
        elseif(thtoc>=maxval(Thicks)) then
            iair=maxloc(Thicks,1)
            call compute_aeroCoeffs_one_airfoil(airfoils(iair),alpha75,alpha5,Re,adotnorm,umach,CL,CD,CN,CT,CLCirc,CM25)
        else
            if(thicks(iint)<=thtoc) then
                iL=iint
                iU=iint+1
            elseif(thicks(iint)>thtoc) then
                iL=iint-1

                iU=iint
            endif
            airfoil_pointer => airfoils(iL)

            call compute_aeroCoeffs_one_airfoil(airfoil_pointer,alpha75,alpha5,Re,adotnorm,umach,&
                                            CLtc(1),CDtc(1),CNtc(1),CTtc(1),CLCirctc(1),CM25tc(1))
            airfoil_pointer => airfoils(iU)
            call compute_aeroCoeffs_one_airfoil(airfoil_pointer,alpha75,alpha5,Re,adotnorm,umach,&
                                            CLtc(2),CDtc(2),CNtc(2),CTtc(2),CLCirctc(2),CM25tc(2))

        xtc = (thtoc-airfoils(iL)%tc)/(airfoils(iU)%tc-airfoils(iL)%tc)
        CL=CLtc(1)+xtc*(CLtc(2)-CLtc(1))
        CD=CDtc(1)+xtc*(CDtc(2)-CDtc(1))
        CN=CNtc(1)+xtc*(CNtc(2)-CNtc(1))
        CT=CTtc(1)+xtc*(CTtc(2)-CTtc(1))
        CLCirc=CLCirctc(1)+xtc*(CLCirctc(2)-CLCirctc(1))
        CM25=CM25tc(1)+xtc*(CM25tc(2)-CM25tc(1))

        endif
        
    end if
    deallocate(thicks,diffthick)

    end subroutine compute_aeroCoeffs

    subroutine compute_aeroCoeffs_one_airfoil(airfoil,alpha75,alpha5,Re,adotnorm,umach,CL,CD,CN,CT,CLCirc,CM25)

    implicit none
   
    ! GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    ! inputs :
    !           alpha75 : Angle of Attack at 3/4 of the element 
    !           alpha5  : Angle of Attack at 1/2 (middle) of the element
    !           Re      : Element Reynolds Number
    !           adotnorm: rate of change of the angle of attack (locally)
    !           umach   : local element Mach Number 
    ! 
    ! outputs: 
    !           CL      : Lift Coefficient
    !           CD      : Drag Coefficient
    !           CN      : Normal Force Coefficient
    !           CD      : Tangential Force Coefficient
    !           CLCirc  :
    !               
    !       
    ! GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    type(AirfoilType),intent(IN) :: airfoil
    real,intent(IN) :: alpha75, alpha5, adotnorm, Re, umach
    real,intent(OUT) :: CL, CD, CN, CT, CLCirc, CM25
    real :: CLstat75, CLstat5, CDstat75, CLdyn5, CDdyn5, dCLAD, dCTAM, dCNAM, CL5, CD5, C, C1, CM25stat
    real :: alphaL, alphaD, aref, Fac   
    
    ewrite(2,*) 'Entering compute_aeroCoeffs_one_airfoil'

    ! Calculate static characteristics
    call intp(Re,alpha75*condeg,CLstat75,CDstat75,CM25stat,airfoil) 
    call intp(Re,alpha5*condeg,CLstat5,C,C1,airfoil)
   
    ! Apply pitch rate effects by analogy to pitching flat plate potential flow theory (SAND report)
    CL5=CLstat75
    CD5=CDstat75
    CM25=CM25stat+cos(alpha5)*(CLstat75-CLstat5)/4.0
    alphaL=alpha75
    CLCirc=CLstat75
    
    ! If no dynamic stall, use static values, else calc dynamic stall
    !if (DSFlag/=0) then

    !    if (DSFlag==1) then
    !        ! Modified Boeing-Vertol approach  
    !        Call BV_DynStall(nElem,CL5,CD5,alphaL,adotnorm,umach,Re,SectInd,CLdyn5,CDdyn5)   
    !    else
    !        ! Leishman-Beddoes model
    !        Call LB_DynStall(nElem,CL5,CD5,alphaL,alpha5,umach,Re,SectInd,CLdyn5,CDdyn5)
    !    end if

    !    CL5=CLdyn5
    !    CD5=CDdyn5  
    !    CLCirc=CLdyn5  

    !end if

    ! Tangential and normal coeffs
    CN=CL5*cos(alpha5)+CD5*sin(alpha5)                                   
    CT=-CL5*sin(alpha5)+CD5*cos(alpha5) 

   ! ! Calc tangential added mass increment by analogy to pitching flat plate potential flow theory (SAND report) 
   ! dCTAM=2.0/cos(alpha5)*wPNorm*CM25stat-CLstat5/2.0*wPNorm
   ! ! Add in alphadot added mass effects (Theodorsen flat plate approx., Katz ch. 13)
   ! dCLAD=pi*adotnorm
   ! dCTAM=dCTAM-dCLAD*sin(alpha5)
   ! dCNAM=dCLAD*cos(alpha5)       

   ! ! Add in added mass effects at low AOA (models not accurate at high AOA)
   ! Fac=1.0
   ! aref=abs(alpha5)
   ! if ((aref > pi/4.0) .AND. (aref < 3.0*pi/4.0)) then
   !     Fac = abs(1-4.0/pi*(aref-pi/4.0))
   ! end if
   ! CT=CT+Fac*dCTAM
   ! CN=CN+Fac*dCNAM

    ! Calc total lift and drag coefficient based on flow direction at half-chord for reference
    CL=CN*cos(alpha5)-CT*sin(alpha5)
    CD=CN*sin(alpha5)+CT*cos(alpha5)


    ewrite(2,*) 'Exiting compute_aeroCoeffs_one_airfoil'

    end subroutine compute_aeroCoeffs_one_airfoil

    subroutine intp(RE,ALPHA,CL,CD,CM25,airfoil)   
        
        implicit none

        real,intent(IN) :: RE, ALPHA
        real,intent(OUT):: CL, CD, CM25
        integer :: i,j               
        real :: XRE, XA 
        real :: CLA(2),CDA(2),CM25A(2)                                      
        type(AirfoilType),intent(IN) :: airfoil 
        integer :: U1, X1, iUB, iLB, NTB, L1
        logical :: NotDone                                               

        !    INTERPOLATE ON RE NO. AND ANGLE OF ATTACK TO GET AIRFOIL CHARACTERISTICS                                            
        CLA(:)=0.0                                                        
        CDA(:)=0.0  
        CM25A(:)=0.0

        if (RE >= airfoil%TRE(1)) then                                                                                                   

            ! Find Re upper and lower bounds.                                     
            NotDone=.true.    
            iUB=2                                                                 
            do while (NotDone)   

                if (RE <= airfoil%TRE(iUB)) then
                    ! Done
                    NotDone=.false.
                    if (RE == airfoil%TRE(iUB)) then
                        iLB=iUB
                    else
                        iLB=iUB-1                                                           
                        XRE=(RE-airfoil%TRE(iLB)/(airfoil%TRE(iUB)-airfoil%TRE(iLB)))
                    end if
                else
                    if (iUB == airfoil%nRET) then       
                        ! warning: no upper bound in table, take last point and set warning...
                        NotDone=.false.                                                       
                        iLB=iUB                                                           
                        XRE=0.0                                                           
                        !airfoil%IUXTP=1
                    else    
                        ! No upper bound, increment and continue                                
                        iUB=iUB+1
                    end if
                end if

            end do

        else        
            ! warning: no lower bound in table, take first point and set warning                                               
            iLB=1                                                             
            iUB=1                                                             
            XRE=0.0                                                                                                 
            !airfoil%ILXTP=1
        end if


        ! INTERPOLATE ON THE ANGLE OF ATTACK                               
        i=1                                                               
        do j=iLB,iUB                                                  

            NTB=airfoil%NTBL(j) ! # of alpha values in table for this section                                                
        
            ! Find upper and lower bound indicies on alpha                                                                     

            ! DO INTERVAL HALVING LOOK UP                                      

            U1=NTB                                                                                              
            L1=1                                                              
            X1=NTB/2 
            NotDone=.true. 

            do while (NotDone)                                                        
                if (ALPHA < airfoil%TA(X1,J)) then
                    U1=X1                                                                                                            
                else    
                    L1=X1 
                end if

                if ((U1-L1) == 1) then
                    NotDone=.false.
                else 
                    X1=L1+(U1-L1)/2
                end if
            end do

            ! DO STRAIGHT LINE INTERPOLATION ON ALPHA                          

            XA=(ALPHA-airfoil%TA(L1,J))/(airfoil%TA(U1,J)-airfoil%TA(L1,J))                  
            CLA(I)=airfoil%TCL(L1,J)+XA*(airfoil%TCL(U1,J)-airfoil%TCL(L1,J))                
            CDA(I)=airfoil%TCD(L1,J)+XA*(airfoil%TCD(U1,J)-airfoil%TCD(L1,J))    
            CM25A(I)=airfoil%TCM(L1,J)+XA*(airfoil%TCM(U1,J)-airfoil%TCM(L1,J)) 

            I=I+1           
        end do

        ! DO STRAIGHT LINE INTERPOLATION ON RE NO.                         

        CL=CLA(1)+XRE*(CLA(2)-CLA(1))                                     
        CD=CDA(1)+XRE*(CDA(2)-CDA(1))  
        CM25=CM25A(1)+XRE*(CM25A(2)-CM25A(1))

    END SUBROUTINE intp

    Subroutine CalcLBStallAOALim(airfoil,Re,CLa,CLCritP,CLCritN)

        ! Get stall data for LB model from airfoil data
        real :: Re, CLa, CLCritP, CLCritN, XRE
        type(AirfoilType),intent(IN) :: airfoil 
        integer :: iUB, iLB
        logical :: NotDone 

        ! Find Re upper and lower bounds.                                     

        if (RE >= airfoil%TRE(1)) then                                                                                                                                  
            NotDone=.true.    
            iUB=2                                                                 
            do while (NotDone)   

                if (RE <= airfoil%TRE(iUB)) then
                    ! Done
                    NotDone=.false.
                    if (RE == airfoil%TRE(iUB)) then
                        iLB=iUB
                    else
                        iLB=iUB-1                                                           
                        XRE=(RE-airfoil%TRE(iLB))/(airfoil%TRE(iUB)-airfoil%TRE(iLB))
                    end if
                else
                    if (iUB == airfoil%nRET) then       
                        ! No upper bound in table, take last point...
                        NotDone=.false.                                                       
                        iLB=iUB                                                           
                        XRE=0.0                                                           
                    else    
                        ! No upper bound, increment and continue                                
                        iUB=iUB+1
                    end if
                end if

            end do

        else        
            ! No lower bound in table, take first point.                                            
            iLB=1                                                             
            iUB=1                                                             
            XRE=0.0                                                                                                 
        end if

        ! Interp
        CLa=airfoil%CLaData(iLB)+xRE*(airfoil%CLaData(iUB)-airfoil%CLaData(iLB))            
        CLCritP=airfoil%CLCritPData(iLB)+xRE*(airfoil%CLCritPData(iUB)-airfoil%CLCritPData(iLB))  
        CLCritN=airfoil%CLCritNData(iLB)+xRE*(airfoil%CLCritNData(iUB)-airfoil%CLCritNData(iLB)) 
    
    End Subroutine CalcLBStallAOALim
    

end module alturbine_airfoils 
