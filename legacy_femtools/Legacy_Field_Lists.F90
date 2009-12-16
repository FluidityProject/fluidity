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

module legacy_field_lists
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN, ident_name
  use fields
  use state_module
  use spud
  implicit none
  
  !! Field name list for tracers (from 1 to NTSOL)
  character(len=FIELD_NAME_LEN), save, &
       dimension(:), allocatable :: field_name_list
  !! Field list for tracters (from 1 to NTSOL)
  type(scalar_field), dimension(:), allocatable, save :: field_list
  !! Options path list for tracers (from 1 to NTSOL)
  character(len=OPTION_PATH_LEN), save, &
       dimension(:), allocatable :: field_optionpath_list
  !! State list for tracers (from 1 to NTSOL)
  integer, save, dimension(:), allocatable :: field_state_list

  private
  public :: field_name_list, field_list, field_optionpath_list,&
       & field_state_list, initialise_field_lists_from_options
contains

  subroutine initialise_field_lists_from_options(state, ntsol)
    type(state_type), dimension(:), intent(in) :: state
    integer, intent(in) :: ntsol

    logical, save:: initialised=.false.
    integer :: nsol, nphases,nfields,ncars,p,f,i, tmpint
    character(len=FIELD_NAME_LEN) :: tmpstring
    logical :: aliased, pressure, density

    integer, dimension(:), allocatable :: priority
    !! Field list for tracers (from 1 to NTSOL)
    character(len=FIELD_NAME_LEN), save, &
        dimension(:), allocatable :: temp_field_name_list
    !! Options path list for tracers (from 1 to NTSOL)
    character(len=OPTION_PATH_LEN), save, &
        dimension(:), allocatable :: temp_field_optionpath_list
    !! State list for tracers (from 1 to NTSOL)
    integer, save, dimension(:), allocatable :: temp_field_state_list
    
    
    ! if called for the second time return immediately
    if (.not.initialised) then
    
       allocate( field_name_list(ntsol), &
            field_state_list(ntsol), &
            field_optionpath_list(ntsol),&
            priority(ntsol), &
            field_list(ntsol), &
            temp_field_name_list(ntsol), &
            temp_field_state_list(ntsol), &
            temp_field_optionpath_list(ntsol) )

       nsol = 0

       nphases = option_count('/material_phase')  
       do p = 0, nphases-1
          nfields = option_count('/material_phase[' &
               //int2str(p)//']/scalar_field')
          do f = 0,nfields-1
             aliased = have_option('/material_phase['// &
                  int2str(p)//']/scalar_field['//int2str(f)//']/aliased')
             call get_option('/material_phase['// &
                  int2str(p)// &
                  ']/scalar_field['//int2str(f)//']/name', &
                  tmpstring)
             call get_option('/material_phase['// &
                  int2str(p)// &
                  ']/scalar_field['//int2str(f)//']/&
                  &prognostic/priority', &
                  tmpint, default=0)
             pressure = (trim(tmpstring)=='Pressure')
             density = (trim(tmpstring)=='Density')

             if ((.not.aliased).and.(.not.pressure).and.(.not.density)) then
                nsol = nsol + 1
                temp_field_name_list(nsol) = tmpstring
                temp_field_optionpath_list(nsol) = '/material_phase['// &
                     int2str(p)// &
                     ']/scalar_field::'//trim(tmpstring)
                temp_field_state_list(nsol) = p+1
                priority(nsol) = tmpint
             end if
          end do

          ! prognostic Mellor Yamada fields:
          if (have_option('/material_phase[' &
               //int2str(p)//']/subgridscale_parameterisations/Mellor_Yamada/scalar_field::KineticEnergy/prognostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "KineticEnergy"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                  ']/subgridscale_parameterisations/Mellor_Yamada/scalar_field::KineticEnergy'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=0)
             priority(nsol) = tmpint
          end if
          if (have_option('/material_phase[' &
               //int2str(p)//']/subgridscale_parameterisations/Mellor_Yamada/scalar_field::TurbulentLengthScalexKineticEnergy/prognostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "TurbulentLengthScalexKineticEnergy"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                  ']/subgridscale_parameterisations/Mellor_Yamada/scalar_field::TurbulentLengthScalexKineticEnergy'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=0)
             priority(nsol) = tmpint
          end if

          if (have_option('/material_phase[' &
               //int2str(p)//']/subgridscale_parameterisations/GLS/scalar_field::GLSTurbulentKineticEnergy/prognostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "GLSTurbulentKineticEnergy"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                  ']/subgridscale_parameterisations/GLS/scalar_field::GLSTurbulentKineticEnergy'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=0)
             priority(nsol) = tmpint
          end if
          if (have_option('/material_phase[' &
               //int2str(p)//']/subgridscale_parameterisations/GLS/scalar_field::GLSGenericSecondQuantity/prognostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "GLSGenericSecondQuantity"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                  ']/subgridscale_parameterisations/GLS/scalar_field::GLSGenericSecondQuantity'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=0)
             priority(nsol) = tmpint
          end if

          do i=0,option_count('/material_phase[' &
               //int2str(p)//']/sediment/sediment_class')-1
             nsol=nsol+1
             call get_option('/material_phase[' &
               //int2str(p)//']/sediment/sediment_class['//int2str(i)&
               //']/name',temp_field_name_list(nsol))
             temp_field_name_list(nsol)='SedimentConcentration'//&
                  trim(temp_field_name_list(nsol))
             temp_field_optionpath_list(nsol)='/material_phase[' &
               //int2str(p)//']/sediment/scalar_field::SedimentTemplate'

             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=0)
             priority(nsol) = tmpint
          end do

       end do

       if(have_option('/traffic_model/scalar_field::TrafficTracerTemplate'))then
          call get_option('/traffic_model/number_of_vehicles', ncars)
          do i=1, ncars
             nsol=nsol+1
             temp_field_name_list(nsol) = "TrafficTracer"//int2str(i)
             temp_field_optionpath_list(nsol) = "/traffic_model/scalar_field::TrafficTracerTemplate"
             temp_field_state_list(nsol)=1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=0)
             priority(nsol)=tmpint
          end do
       end if

       ! make sure we have found all ntsol scalar fields:
       assert(nsol==ntsol)

       nsol=0
       do p=maxval(priority),minval(priority),-1
          do f=1,ntsol
             if (priority(f)==p) then
                nsol = nsol + 1
                field_name_list(nsol) = temp_field_name_list(f)
                field_optionpath_list(nsol) = temp_field_optionpath_list(f)
                field_state_list(nsol) = temp_field_state_list(f)
             end if
          end do
       end do

       deallocate( priority, &
            temp_field_name_list, &
            temp_field_state_list, &
            temp_field_optionpath_list )

       initialised = .true.
       
    end if ! End of if(initialised)

    ! Point the list of fields. This has to be done every adapt as the
    ! field structures will be reallocated.
    
    ! Note that we use borrowed references for this so as not to interfere
    ! with adaptivity.
    do f=1,ntsol
       field_list(f) = extract_scalar_field(state(field_state_list(f)),&
            &                               field_name_list(f))
    end do
          
  end subroutine initialise_field_lists_from_options
    
  subroutine initialise_field_lists_from_file(ntsol, ident, tphase)
  integer, intent(in):: ntsol
  integer, dimension(1:ntsol), intent(in):: ident, tphase
  integer :: i
    
    logical, save:: initialised=.false.
    
    ! if called for the second time return immediately
    if (initialised) return
    
    allocate( field_name_list(ntsol), &
      field_state_list(ntsol), &
      field_optionpath_list(ntsol) )
      
    field_name_list = ident_name(ident)
    ! field_optionpath_list is obviously not filled in
    field_state_list= tphase

    do i=1,ntsol
      field_optionpath_list(i) = "/material_phase[" // int2str(tphase(i) - 1) // "]/" // &
                       & "scalar_field::" // ident_name(ident(i))
    end do

    initialised = .true.
    
  end subroutine initialise_field_lists_from_file


end module legacy_field_lists
