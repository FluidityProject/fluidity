!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    C.Pain@Imperial.ac.uk
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

module assemble_buoyancy
  use FLDebug
  use equation_of_state
  use fields
  use state_module
  use spud
  use state_module
  use AllSorts
  use legacy_field_lists
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use equation_of_state
  implicit none

  private
  public :: bsines

contains

  SUBROUTINE BSINES(IT,&
                                !cccccccccccccccccccccccccccccccccccccccccccc
                                !     salinity stuff - cjc
       &              dengam_sal, s0, ztop, &
       &              NTSUB,NSALSUB,&
       &              NLOC,TOTELE, &
       &              NSOGRASUB, &
       &              ITS,ITHETA,&
       &              BSOUX,BSOUY,BSOUZ,&
       &              NONODS,&
       &              BOUSIN,BHOUT,CHADEN,&
       &              TEMINI,DENINI,DENGAM,&
       &              COGRAX,COGRAY,COGRAZ,&
       &              GRAVTY,RAD,&
       &              GAS,D3,NPHASE,&
       &              state)
    INTEGER, intent(in) :: nonods, it
    real, intent(in) :: dengam_sal, s0, ztop
    INTEGER, intent(in) :: NTSUB,NSALSUB,NLOC,TOTELE
    INTEGER, intent(in) :: NSOGRASUB
    REAL, intent(in) :: BSOUX,BSOUY,BSOUZ
    ! (BSOUX,BSOUY,BSOUZ) = the direction of gravity. (usually -ve)
    INTEGER, intent(in) :: ITS
    REAL, intent(in) :: TEMINI,DENINI,DENGAM,GRAVTY,ITHETA
    LOGICAL, intent(in) :: D3,COGRAX,COGRAY,COGRAZ
    INTEGER, intent(in) :: NPHASE
    type(state_type), intent(inout):: state
    LOGICAL, intent(in) :: BOUSIN,BHOUT,CHADEN,RAD,GAS

    !local variables
    real :: salt, density
    logical :: gottopbotdis,  gotsal
    integer :: salinity_option
    type(vector_field) :: Coordinate, GravityDirection
    type(scalar_field) :: CoordinateZ, temperature_old
    type(scalar_field) :: topdis_field, salinity,density_field
    type(scalar_field) :: temperature, tsub_field, tsub_field_old
    type(scalar_field) :: saltsub_field, saltsub_field_old,salinity_old
    type(scalar_field) :: velocitybuoyancydensity
    real, dimension(:), pointer :: sal,T, denpt, Z, tsub,prevtsub,prevt
    real, dimension(:), pointer :: salsub, prevsalsub, salprevt
    real, dimension(:), pointer :: buoyancy
    real, dimension(:), pointer :: VGRAVX, VGRAVY, VGRAVZ
    real, allocatable, dimension(:,:) :: local_vec
    integer, dimension(:), pointer :: ndglno

    REAL TEMP, topdisi
    REAL RTSUB,RPREVTSUB,RSALSUB,RPREVSALSUB
    INTEGER I,ELE,ILOC,NOD, stat
    ! If COGRAX then we do not have a vector of gravity in VGRAVX
    ! used when centrifugal forces exist. 
    ! This sub applies the Boussinesq approx in some form. 
    ! If BHOUT then subtract out Hydrostatic level. 
    ! If GAS then assume phase is an incompressible gas. 
    ! ITS=No of non-linear iterations.
    !

    ewrite(1,*) 'just inside bsines'


    sal=>null()
    T=>null()
    denpt=>null()
    Z=>null()
    tsub=>null()
    prevtsub=>null()
    prevt=>null()
    salsub=>null()
    prevsalsub=>null()
    salprevt=>null()
    buoyancy=>null()
    VGRAVX=>null()
    VGRAVY=>null()
    VGRAVZ=>null()
    ndglno=>null()

    allocate( local_vec(3,nonods) )

    gotsal = have_option('/material_phase[0]/scalar_field::Salinity')
    gottopbotdis = have_option('/material_phase[0]/scalar_field::FreeSurface')
    gottopbotdis = gottopbotdis.or.has_scalar_field(state, "DistanceToTop")

    if (.not.have_option("/material_phase[0]/equation_of_state"))&
         & then
       SALINITY_OPTION=0
    else if (have_option("/material_phase[0]/equation_of_state/fluids"))&
         & then
       if (have_option("/material_phase[0]/equation_of_state/fluids/linear")) then
          if(have_option("/material_phase[0]/equation_of_state/fluids/linear/salinity_dependency")) then
             SALINITY_OPTION=1
          else
             SALINITY_OPTION=0
          end if
       else if (have_option("/material_phase[0]/equation_of_state/&
            &fluids/ocean_pade_approximation")) then
          if(have_option("/material_phase[0]/equation_of_state/&
               &fluids/ocean_pade_approximation/&
               &include_depth_below_surface")) then
             SALINITY_OPTION=2
          else
             SALINITY_OPTION=3
          end if
       end if
    end if

    if(gotsal) then
       salinity = extract_scalar_field(state, "Salinity")
       sal => salinity%val
       salinity_old = extract_scalar_field(state, "OldSalinity")
       salprevt => salinity_old%val
    end if
    temperature = extract_scalar_field(state, field_name_list(it))
    T => temperature%val
    temperature_old = extract_scalar_field(state, "Old"//trim(field_name_list(it)))
    prevT => temperature_old%val
    density_field = extract_scalar_field(state, "Density")
    denpt => density_field%val
    Coordinate = extract_vector_field(state, "Coordinate")
    if(d3) then
       CoordinateZ = extract_scalar_field_from_vector_field(Coordinate,3)
       Z => COordinateZ%val
    end if
    GravityDirection=extract_vector_field(state, "GravityDirection")
    vgravx => GravityDirection%val(1)%ptr
    vgravy => GravityDirection%val(2)%ptr
    if(d3) then
       vgravz => GravityDirection%val(3)%ptr
    end if
    VelocityBuoyancyDensity = &
         extract_scalar_field(state,"VelocityBuoyancyDensity", stat)
    if(stat==0) then
      buoyancy => VelocityBuoyancyDensity%val
    else
      ! no point in doing anything
      return
    end if

    ndglno => temperature%mesh%ndglno
    
    if (have_option(trim(temperature%option_path)// &
         "/prognostic/spatial_discretisation/inner_element")) then
       saltsub_field=extract_scalar_field(state,trim(field_name_list(it)) &
            //"InnerElement")
       salsub => saltsub_field%val
       saltsub_field_old=extract_scalar_field(state,"Old"// &
            trim(field_name_list(it)) &
            //"InnerElement")
       prevsalsub => saltsub_field_old%val
    end if

    if(gotsal) then
       if (have_option(trim(salinity%option_path)// &
            "/prognostic/spatial_discretisation/inner_element")) then
          tsub_field=extract_scalar_field(state,trim(field_name_list(it)) &
               //"InnerElement")
          tsub => tsub_field%val
          tsub_field_old=extract_scalar_field(state,"Old"// &
               trim(field_name_list(it)) &
               //"InnerElement")
          prevtsub => tsub_field_old%val
       end if
    end if

    if (gottopbotdis) then
       topdis_field = extract_scalar_field(state, "DistanceToTop")
    else if (gotsal .and. (salinity_option==2)) then
       ewrite(0, *) 'Salinity option 2 needs distance to top'
       ewrite(0, *) 'This is only available in combination '
       ewrite(0, *) 'with a free surface field (ident -29)'
       FLAbort("Sorry!")
    end if
    ! Put all sources to zero some already are:
    buoyancy = 0.0

    IF(COGRAX) THEN
       LOCAL_VEC(1,:)=BSOUX
    ELSE
       LOCAL_VEC(1,:)=VGRAVX
    ENDIF
    IF(COGRAY) THEN
       LOCAL_VEC(2,:)=BSOUY
    ELSE
       LOCAL_VEC(2,:)=VGRAVY
    ENDIF
    IF(D3) THEN
       IF(COGRAZ) THEN
          LOCAL_VEC(3,:)=BSOUZ
       ELSE
          LOCAL_VEC(3,:)=VGRAVZ
       ENDIF
    ENDIF

    ewrite(1,*) 'before if(bousin)...'
    ewrite(2,*) 'gotsal,ITHETA,its',gotsal,ITHETA,its
    IF(BOUSIN) THEN
       IF(CHADEN) THEN
          ! the following actually changes the density for the boussinesq
          ! approx.
          ewrite(1,*) 'in IF(CHADEN)...'
          DENPT = DENINI*(1.-DENGAM*(T-TEMINI))
       ENDIF
       IF(BHOUT) THEN
          if(gotsal) then
             ! the following is for the standard boussinesq approx. with
             ! salinity - cjc
             ewrite(1,*) 'in if(gotsal)...'
             do I=1,NONODS
                IF(ITS.GE.2) THEN
                   TEMP=ITHETA*T(I)+(1.-ITHETA)*PREVT(I)
                ELSE
                   TEMP=T(I)
                   ! Theta differencing for salinity - cjc
                endif
                IF(ITS.GE.2) THEN
                   salt=ITHETA*sal(I)+(1.-ITHETA)*Salprevt(I)
                ELSE
                   salt=sal(I)
                ENDIF
                ! The direction of gravity is usually -ve. 
                if (gottopbotdis) then
                   topdisi=node_val(topdis_field, i)
                else
                   topdisi=0.0 ! not used
                end if
                CALL GETDEN_TEMP_SALT(&
                     & salinity_option, DENINI, DENGAM, TEMINI, DENGAM_sal,&
                     & S0, density, temp, salt, topdisi)

                buoyancy(I)=density

             end do

          else
             ! the following is for the standard boussinesq approx.
             ! without salinity
             ewrite(1,*) 'calculating standard boussinesq approx. without salinity'
             do I=1,NONODS
                IF(ITS.GE.2) THEN
                   TEMP=ITHETA*T(I)+(1.-ITHETA)*PREVT(I)
                ELSE
                   TEMP=T(I)
                ENDIF
                density=-DENINI*DENGAM*(TEMP-TEMINI)
                ! The direction of gravity is usually -ve. 
                buoyancy(I)=density
             end do
          endif
       ELSE
          do I=1,NONODS
             IF(ITS.GE.2) THEN
                TEMP=ITHETA*T(I)+(1.-ITHETA)*PREVT(I)
             ELSE
                TEMP=T(I)
             ENDIF
             buoyancy(I)=DENINI*(1.-DENGAM*(TEMP-TEMINI))
          end do
          ! end of IF(BHOUT) THEN ELSE...
       ENDIF
    ENDIF

    deallocate( local_vec )

  end subroutine BSINES

end module assemble_buoyancy
