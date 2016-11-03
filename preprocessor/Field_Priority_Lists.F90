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

#include "fdebug.h"

module field_priority_lists

  use fldebug
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use futils, only: int2str
  use spud
  use fields
  use state_module
  use sediment, only: get_n_sediment_fields, get_sediment_item

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
       & field_state_list, initialise_field_lists_from_options,&
       & get_ntsol 
       
contains

  subroutine initialise_field_lists_from_options(state, ntsol)
    type(state_type), dimension(:), intent(in) :: state
    integer, intent(in) :: ntsol

    logical, save:: initialised=.false.
    integer :: nsol, nphases,nfields,ncars,p,f,i, tmpint
    character(len=FIELD_NAME_LEN) :: tmpstring
    logical :: aliased, pressure

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

             if (.not. aliased .and. .not. pressure) then
                nsol = nsol + 1
                temp_field_name_list(nsol) = tmpstring
                temp_field_optionpath_list(nsol) = '/material_phase['// &
                     int2str(p)// &
                     ']/scalar_field::'//trim(tmpstring)
                temp_field_state_list(nsol) = p+1
                priority(nsol) = tmpint
             end if
          end do

          ! prognostic sediment fields
          if (have_option('/material_phase['//int2str(p)//']/sediment')) then
             nfields = get_n_sediment_fields()
             do f = 1, nfields
                nsol=nsol+1
                
                call get_sediment_item(state(p+1), f, temp_field_name_list(nsol))
                temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                     ']/sediment/scalar_field['//int2str(f-1)//']'
                temp_field_state_list(nsol) = p+1
                call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                     tmpint, default=0)
                priority(nsol) = tmpint
             end do
          end if

          ! this whole set up of fields could be improved to ensure that multiple PBEs can be used per phase
          ! prognostic pop balance fields - very limited applicability
          if (have_option('/material_phase['//int2str(p)//']/population_balance/')) then
             do f = 0, option_count('/material_phase['//int2str(p)//&
                  ']/population_balance/weights/scalar_field') - 1
                nsol=nsol+1
                temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                     ']/population_balance/weights/scalar_field['//int2str(f)//']'
                call get_option('/material_phase['//int2str(p)//&
                     ']/population_balance/weights/scalar_field['//int2str(f)//&
                     ']/name',temp_field_name_list(nsol))
                call get_option('/material_phase['//int2str(p)//&
                     ']/population_balance/weights/scalar_field['//int2str(f)//&
                     ']/prognostic/priority', priority(nsol), default=0)
                temp_field_state_list(nsol) = p+1
             end do
             do f = 0, option_count('/material_phase['//int2str(p)//&
                  ']/population_balance/weighted_abscissa/scalar_field') - 1
                nsol=nsol+1
                temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                     ']/population_balance/weighted_abscissa/scalar_field['//int2str(f)//']'
                call get_option('/material_phase['//int2str(p)//&
                     ']/population_balance/weighted_abscissa/scalar_field['//int2str(f)//&
                     ']/name',temp_field_name_list(nsol))
                call get_option('/material_phase['//int2str(p)//&
                     ']/population_balance/weighted_abscissa/scalar_field['//int2str(f)//&
                     ']/prognostic/priority', priority(nsol), default=0)
                temp_field_state_list(nsol) = p+1
             end do
          end if

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

          ! Check for GLS - we need to make sure these fields are solved *after*
          ! everything else, so set to a big negative value. In addition, the
          ! Psi solve *must* come after the TKE solve, so make sure the priority
          ! is set such that this happens
          if (have_option('/material_phase[' &
               //int2str(p)//']/subgridscale_parameterisations/GLS/scalar_field::GLSTurbulentKineticEnergy/prognostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "GLSTurbulentKineticEnergy"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                  ']/subgridscale_parameterisations/GLS/scalar_field::GLSTurbulentKineticEnergy'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=nsol)
             priority(nsol) = -tmpint*100
          end if
          if (have_option('/material_phase[' &
               //int2str(p)//']/subgridscale_parameterisations/GLS/scalar_field::GLSGenericSecondQuantity/prognostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "GLSGenericSecondQuantity"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                  ']/subgridscale_parameterisations/GLS/scalar_field::GLSGenericSecondQuantity'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=nsol)
             priority(nsol) = -tmpint*100
          end if

          ! Check for k-epsilon - we need to make sure these fields are solved *after*
          ! everything else, so set to a big negative value. In addition, the
          ! TurbulentDissipation (Epsilon) solve *must* come after the TKE solve,
          ! so make sure the priority is set such that this happens.
          if (have_option('/material_phase[' &
               //int2str(p)//']/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentKineticEnergy/prognostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "TurbulentKineticEnergy"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                  ']/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentKineticEnergy'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=nsol)
             priority(nsol) = -tmpint*100
          end if
          if (have_option('/material_phase[' &
               //int2str(p)//']/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentDissipation/prognostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "TurbulentDissipation"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
                  ']/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentDissipation'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=nsol)
             priority(nsol) = -tmpint*100
          end if
          ! Check for subgrid-scale kinetic energy equation
          ! - we need to make sure this is solved *after*
          ! everything else, so set to a big negative value.
          if(have_option('/material_phase['//int2str(p)// &
             ']/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/scalar_field::SubgridKineticEnergy')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "SubgridKineticEnergy"
             temp_field_optionpath_list(nsol)='/material_phase['//int2str(p)// &
             ']/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/scalar_field::SubgridKineticEnergy'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/prognostic/priority', &
                  tmpint, default=nsol)
             priority(nsol) = -tmpint*200
          end if
!!! Melt rate should be the last thing to calculate, Sb
          if (have_option('/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::Sb/diagnostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "Sb"
             temp_field_optionpath_list(nsol)='/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::Sb'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/diagnostic/priority', &
                  tmpint, default=nsol)
             priority(nsol) = -tmpint*200
          end if
!!! Melt rate should be the last thing to calculate, Tb
          if (have_option('/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::Tb/diagnostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "Tb"
             temp_field_optionpath_list(nsol)='/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::Tb'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/diagnostic/priority', &
             tmpint, default=nsol)
             priority(nsol) = -tmpint*200
          end if
!!!/ocean_forcing/iceshelf_meltrate/Holland08
          if (have_option('/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::MeltRate/diagnostic')) then
             nsol=nsol+1
             temp_field_name_list(nsol) = "MeltRate"
             temp_field_optionpath_list(nsol)='/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::MeltRate'
             temp_field_state_list(nsol) = p+1
             call get_option(trim(temp_field_optionpath_list(nsol))//'/diagnostic/priority', &
                  tmpint, default=nsol)
             priority(nsol) = -tmpint*200
          end if
       end do

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

  subroutine get_ntsol(ntsol)
    integer, intent(out) :: ntsol
    integer :: nphases,nfields,ncars,p,f
    character(len=FIELD_NAME_LEN) :: tmpstring
    logical :: aliased, pressure

    ntsol = 0

    nphases = option_count('/material_phase')  
    do p = 0, nphases-1
       nfields = option_count('/material_phase[' &
            //int2str(p)//']/scalar_field')
       do f = 0, nfields-1
          aliased = have_option('/material_phase['// &
               int2str(p)//']/scalar_field['//int2str(f)//']/aliased')
          call get_option('/material_phase['// &
               int2str(p)//']/scalar_field['//int2str(f)//']/name', tmpstring)
          pressure = (trim(tmpstring)=='Pressure')

          if (.not. aliased .and. .not. pressure) then
             ntsol = ntsol + 1
          end if
       end do
    
       ! added as hack for now - but this whole set up of fields could be way better!
       ! prognostic pop balance fields - very limited applicability
       if (have_option('/material_phase['//int2str(p)//']/population_balance/')) then
          ntsol = ntsol + &
               option_count('/material_phase['//int2str(p)//&
               ']/population_balance/weights/scalar_field') + &
               option_count('/material_phase['//int2str(p)//&
               ']/population_balance/weighted_abscissa/scalar_field')
       end if

       ! prognostic scalar fields for Mellor Yamada:
       if (have_option('/material_phase[' &
            //int2str(p)//']/subgridscale_parameterisations/Mellor_Yamada/scalar_field::KineticEnergy/prognostic')) then
          ntsol=ntsol + 1
       end if
       if (have_option('/material_phase[' &
            //int2str(p)//']/subgridscale_parameterisations/Mellor_Yamada/scalar_field::TurbulentLengthScalexKineticEnergy/prognostic')) then
          ntsol=ntsol + 1
       end if
       if (have_option('/material_phase[' &
            //int2str(p)//']/subgridscale_parameterisations/GLS/scalar_field::GLSTurbulentKineticEnergy/prognostic')) then
          ntsol=ntsol + 1
       end if
       if (have_option('/material_phase[' &
            //int2str(p)//']/subgridscale_parameterisations/GLS/scalar_field::GLSGenericSecondQuantity/prognostic')) then
          ntsol=ntsol + 1
       end if
       ! prognostic scalar fields for k-epsilon turbulence model:
       if (have_option('/material_phase[' &
            //int2str(p)//']/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentKineticEnergy/prognostic')) then
          ntsol=ntsol + 1
       end if
       if (have_option('/material_phase[' &
            //int2str(p)//']/subgridscale_parameterisations/k-epsilon/scalar_field::TurbulentDissipation/prognostic')) then
          ntsol=ntsol + 1
       end if
       ! prognostic scalar fields for subgrid-scale kinetic energy model:
       if(have_option('/material_phase['//int2str(p)// &
            ']/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/scalar_field::SubgridKineticEnergy')) then
          ntsol=ntsol + 1
       end if
       !Melting
       if (have_option('/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::Tb/diagnostic')) then
          ntsol=ntsol + 1
       end if
       if (have_option('/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::Sb/diagnostic')) then
          ntsol=ntsol + 1
       end if
       if (have_option('/ocean_forcing/iceshelf_meltrate/Holland08/scalar_field::MeltRate/diagnostic')) then
          ntsol=ntsol + 1
       end if
       !Sediments
       if (have_option('/material_phase['//int2str(p)//']/sediment')) then
          ntsol=ntsol + get_n_sediment_fields()
       end if

    end do

  end subroutine get_ntsol

end module field_priority_lists
