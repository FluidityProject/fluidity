#include "fdebug.h"

module Airfoils 

  use fldebug
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, PYTHON_FUNC_LEN, pi

  implicit none
    
  real, parameter :: conrad = pi / 180.0 
  real, parameter :: condeg = 180.0 / pi  
  
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

  ! Airfoil parameterss for Leishman-Beddoes Dynamic stall model
  real, allocatable :: CLaData(:)
  real, allocatable :: alphaSSPData(:)
  real, allocatable :: alphaSSNData(:)

  end type AirfoilType
 
  ! Maximum Numbers of Reynolds Number data that can be stored
  ! Globla parameters
  integer, parameter :: MaxReVals = 20
  integer, parameter :: MaxAOAVals = 1000
  
  ! Private subroutines
  private EvalStaticCoeff, read_airfoil
 
  ! Public subroutines
  public airfoil_init_data, compute_StaticLoads ,allocate_airfoil, copy_airfoil_values


contains
   
    subroutine airfoil_init_data(airfoil)
    
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

    end subroutine airfoil_init_data

    subroutine copy_airfoil_values(airfoil1,airfoil2)
    implicit none
    type(AirfoilType),intent(INOUT) :: airfoil1,airfoil2
    ewrite(2,*) 'Entering copy_airfoil_values'

    airfoil1%afname=airfoil2%afname
    airfoil1%aftitle=airfoil2%aftitle
    airfoil1%camb=airfoil2%camb
    airfoil1%tc=airfoil2%tc
    airfoil1%TA(1:MaxAOAVals,1:MaxReVals)=airfoil2%TA(1:MaxAOAVals,1:MaxReVals)
    airfoil1%TCL(1:MaxAOAVals,1:MaxReVals)=airfoil2%TCL(1:MaxAOAVals,1:MaxReVals)
    airfoil1%TCD(1:MaxAOAVals,1:MaxReVals)=airfoil2%TCD(1:MaxAOAVals,1:MaxReVals)
    airfoil1%TCM(1:MaxAOAVals,1:MaxReVals)=airfoil2%TCM(1:MaxAOAVals,1:MaxReVals)
    airfoil1%TRE(1:MaxReVals)=airfoil2%TRE(1:MaxReVals)
    airfoil1%nTBL(1:MaxReVals)=airfoil2%nTBL(1:MaxReVals)
    airfoil1%CLaData(1:MaxReVals)=airfoil2%CLaData(1:MaxReVals)
    airfoil1%alphaSSPData(1:MaxReVals)=airfoil2%alphaSSPData(1:MaxReVals)
    airfoil1%alphaSSNData(1:MaxReVals)=airfoil2%alphaSSNData(1:MaxReVals)

    ewrite(2,*) 'Exiting copy_airfoil_values'

    end subroutine copy_airfoil_values
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    ! These routines have been taken directly out of CACTUS and changed
    ! slightly so they fit the needs of the present coupling
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

    subroutine read_airfoil(airfoil)
    
    implicit none
    
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    ! This routine reads the airfoil data from 
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    
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
        i=0 
    do while (EOF>=0 .and. (i<MaxReVals))
        i=i+1
        ! Read Re and dyn. stall data
        read(ReadLine(index(ReadLine,':')+1:),*) airfoil%TRE(i)
        read(15,'(A)') ReadLine                          
        read(ReadLine(index(ReadLine,':')+1:),*) airfoil%CLaData(i)
        read(15,'(A)') ReadLine                          
        read(ReadLine(index(ReadLine,':')+1:),*) airfoil%alphaSSPData(i)
        read(15,'(A)') ReadLine                          
        read(ReadLine(index(ReadLine,':')+1:),*) airfoil%alphaSSNData(i)

        ! Reverse camber direction if desired
        if (airfoil%camb == 1) then
            temp = airfoil%alphaSSPData(i)
            airfoil%alphaSSPData(i) = -airfoil%alphaSSNData(i)
            airfoil%alphaSSNData(i) = -temp 
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

    ewrite(1,*) 'In allocate_airfoil -- ', airfoil%afname
    allocate(airfoil%TA(MaxAOAVals,MaxReVals))
    allocate(airfoil%TCL(MaxAOAVals,MaxReVals))
    allocate(airfoil%TCD(MaxAOAVals,MaxReVals))
    allocate(airfoil%TCM(MaxAOAVals,MaxReVals))
    allocate(airfoil%TRE(MAxReVals))
    allocate(airfoil%nTBL(MaxReVals))
    allocate(airfoil%CLaData(MaxReVals))
    allocate(airfoil%alphaSSPData(MaxReVals))
    allocate(airfoil%alphaSSNData(MaxReVals))

    ewrite(1,*) 'Exiting allocate_airfoil'
    
    end subroutine allocate_airfoil

    subroutine compute_StaticLoads(airfoil,alpha,Re,CN,CT,CM25,CL,CD)

        implicit none

        ! GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
        ! inputs :
        !           airfoil : The airfoil under consideration
        !           alpha   : Angle of Attack at 1/4 of the element
        !           Re      : Element Reynolds Number
        ! 
        ! outputs: 
        !           CN      : Normal Force Coefficient
        !           CT      : Tangential Force Coefficient
        !           CM25    : Pitch moment at 1/4 from the LE
        !           CL      : Lift Coefficient
        !           CD      : Drag Coefficient
        !       
        ! GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
        type(AirfoilType),intent(IN) :: airfoil
        real,intent(IN) :: alpha, Re
        real,intent(OUT) :: CN, CT, CM25, CL, CD

        ewrite(2,*) 'Entering compute_StaticLoads'

        !========================================================
        ! Find Static Coefficients by interpolation Static Loads
        !========================================================
        call EvalStaticCoeff(Re,alpha*condeg,CL,CD,CM25,airfoil) 

        ! Tangential and normal coeffs
        CN=CL*cos(alpha)+CD*sin(alpha)                                   
        CT=-CL*sin(alpha)+CD*cos(alpha) 
 
        return
        ewrite(2,*) 'Exiting compute_StaticLoads'

    end subroutine compute_StaticLoads
 
    subroutine EvalStaticCoeff(RE,ALPHA,CL,CD,CM25,airfoil)   
        
        implicit none

        real,intent(IN) :: RE, ALPHA
        real,intent(OUT):: CL, CD, CM25
        integer :: i,j               
        real :: XRE, XA 
        real,dimension(2) :: CLA,CDA,CM25A                                      
        type(AirfoilType),intent(IN) :: airfoil 
        integer :: U1, X1, iUB, iLB, NTB, L1
        logical :: NotDone                                               
        
        ewrite(2,*) 'Entering intp subroutine'
    ! INTERPOLATE ON RE NO. AND ANGLE OF ATTACK TO GET AIRFOIL CHARACTERISTICS                                            
        CLA(:)=0.0                                                        
        CDA(:)=0.0  
        CM25A(:)=0.0         
        if (RE >= airfoil%TRE(1)) then ! Find Re upper and lower bounds.                                     
            NotDone=.true.    
            iUB=2 ! maxloc(airfoil%TRE,1) ?                                                                 
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
                        ! If we have exceeded the maximum number of the reynolds number then take the 
                        NotDone=.false.                                                       
                        iLB=iUB                                                           
                        XRE=0.0                                                           
ewrite(2,*) 'Warning : The upper Reynolds number available data was exceeded. Calculate CD,CL,CM with : Re = ', airfoil%TRE(iUB)
                    exit
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
    ewrite(2,*) 'Warning : The lower Reynolds number available data was exceeded. Calculate CD,CL,CM with : Re = ', airfoil%TRE(iLB)
        end if
        ! INTERPOLATE ON THE ANGLE OF ATTACK                               
        I=1                                                               
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

        ewrite(2,*) 'Exiting intp subroutine'
    END SUBROUTINE EvalStaticCoeff

    Subroutine EvalStaticStallParams(airfoil,Re,alphaSSP,alphaSSN,CLAlpha)

        implicit none
        type(AirfoilType),intent(IN) :: airfoil
        real,intent(in) :: Re
        real, intent(out) :: alphaSSP, alphaSSN, CLAlpha
        real :: XRE
        integer :: iUB, iLB
        logical :: NotDone 

        ! Find Re upper and lower bounds.                               

        if(RE>=airfoil%TRE(1))then            
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
            iLB=1                                                             
            iUB=1                                                             
            XRE=0.0                                                      
        end if

        ! Interp
        CLAlpha=airfoil%CLaData(iLB)+xRE*(airfoil%CLaData(iUB)-airfoil%CLaData(iLB))            
        alphaSSP=airfoil%alphaSSPData(iLB)+xRE*(airfoil%alphaSSPData(iUB)-airfoil%alphaSSPData(iLB))  
        alphaSSN=airfoil%alphaSSNData(iLB)+xRE*(airfoil%alphaSSNData(iUB)-airfoil%alphaSSNData(iLB)) 

    End Subroutine EvalStaticStallParams

end module Airfoils 
