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
#include "fdebug.h"

module timeloop_utilities
  use fldebug
  use spud
  use global_parameters, only: simulation_start_cpu_time,&
    & simulation_start_wall_time, OPTION_PATH_LEN
  use parallel_tools
  use fields
  use state_module
  use fefields
  use signal_vars
  use timers
  implicit none
  
  private

  public :: copy_to_stored_values, copy_from_stored_values,&
       & relax_to_nonlinear, simulation_completed, get_copied_field

contains

  subroutine copy_to_stored_values(state, prefix)
    !!< For each field, copy its value to prefixfield if prefixfield is present.
    type(state_type), dimension(:), intent(inout) :: state
    character(len=*), intent(in) :: prefix

    integer :: s, f, stat
    
    type(scalar_field) :: sfield, old_sfield
    type(vector_field) :: vfield, old_vfield
    type(tensor_field) :: tfield, old_tfield

    do s=1,size(state)

       do f=1,scalar_field_count(state(s))

          sfield=extract_scalar_field(state(s), f)

          if(.not.aliased(sfield)) then
                
             old_sfield=extract_scalar_field(state(s), trim(prefix)//sfield%name,&
                  & stat=stat)
             
             if ((stat==0).and.(.not.aliased(old_sfield))) then
                ! In this case there is an old field to be set.
                call set(old_sfield, sfield)
             end if
             
          end if

       end do

       do f=1,vector_field_count(state(s))

          vfield=extract_vector_field(state(s), f)

          if(.not.aliased(vfield)) then

             ! Special case: do not copy to the coordinates
             if ((vfield%name=="Coordinate")) then
                cycle
             end if

             old_vfield=extract_vector_field(state(s), trim(prefix)//vfield%name,&
                  & stat=stat)

             if ((stat==0).and.(.not.aliased(old_vfield))) then
                ! In this case there is an old field to be set.
                call set(old_vfield, vfield)
             end if

          end if

       end do
       
       do f=1,tensor_field_count(state(s))
          
          tfield=extract_tensor_field(state(s), f)
          
          if(.not.aliased(tfield)) then

             old_tfield=extract_tensor_field(state(s), trim(prefix)//tfield%name,&
                  & stat=stat)
  
             if ((stat==0).and.(.not.aliased(old_tfield))) then
                ! In this case there is an old field to be set.
                call set(old_tfield, tfield)
             end if

          end if

       end do

    end do

  end subroutine copy_to_stored_values

  subroutine copy_from_stored_values(state, prefix)
    !!< For each field, copy its value from prefixfield if prefixfield is present.
    type(state_type), dimension(:), intent(inout) :: state
    character(len=*), intent(in) :: prefix

    integer :: s, f, stat
    
    type(scalar_field) :: sfield, old_sfield
    type(vector_field) :: vfield, old_vfield
    type(tensor_field) :: tfield, old_tfield

    do s=1,size(state)

       do f=1,scalar_field_count(state(s))

          sfield=extract_scalar_field(state(s), f)

          if(.not.(have_option(trim(sfield%option_path)//"/prescribed"))) then
             
             if(.not.aliased(sfield)) then

                ! Special case: do not copy back pressure or density or geostrophic pressure
                if ((sfield%name=="Pressure").or.(sfield%name=="Density").or.(sfield%name=="GeostrophicPressure")) then
                   cycle
                end if

                old_sfield=extract_scalar_field(state(s), trim(prefix)//sfield%name,&
                     & stat=stat)

                if ((stat==0).and.(.not.aliased(old_sfield))) then
                   ! In this case there is an old field to be set.
                   call set(sfield, old_sfield)
                end if

             end if

          end if
          
       end do

       do f=1,vector_field_count(state(s))

          vfield=extract_vector_field(state(s), f)

          if(.not.(have_option(trim(vfield%option_path)//"/prescribed"))) then
             
             if(.not.aliased(vfield)) then

                ! Special case: do not copy back the coordinates or the gridvelocity
                if ((vfield%name=="Coordinate").or.(vfield%name=="GridVelocity")) then
                   cycle
                end if

                old_vfield=extract_vector_field(state(s), trim(prefix)//vfield%name,&
                     & stat=stat)
                
                if ((stat==0).and.(.not.aliased(old_vfield))) then
                   ! In this case there is an old field to be set.
                   call set(vfield, old_vfield)
                end if
                
             end if
             
          end if

       end do
       
       do f=1,tensor_field_count(state(s))

          tfield=extract_tensor_field(state(s), f)

          if(.not.(have_option(trim(tfield%option_path)//"/prescribed"))) then
             
             if(.not.aliased(tfield)) then

                old_tfield=extract_tensor_field(state(s), trim(prefix)//tfield%name,&
                     & stat=stat)

                if ((stat==0).and.(.not.aliased(old_tfield))) then
                   ! In this case there is an old field to be set.
                   call set(tfield, old_tfield)
                end if

             end if

          end if

       end do

    end do

  end subroutine copy_from_stored_values

  subroutine get_copied_field(fieldname, state)

    type(state_type), intent(in) :: state
    character(len=*), intent(in) :: fieldname

    type(scalar_field), pointer :: copiedfield
    type(scalar_field), pointer :: tmpfield
    character(len=OPTION_PATH_LEN) :: tmpstring

    if(trim(fieldname)=="CopiedField") then
      copiedfield=>extract_scalar_field(state, "CopiedField")
      call get_option(trim(copiedfield%option_path)//"/prognostic/copy_from_field", &
                              tmpstring)
      tmpfield=>extract_scalar_field(state, "Old"//trim(tmpstring))
      call set(copiedfield, tmpfield)
    end if

  end subroutine get_copied_field

  subroutine relax_to_nonlinear(state)
    !!< For each field, set the nonlinearfield if present.
    type(state_type), dimension(:), intent(inout) :: state

    integer :: s, f, old_stat, nl_stat, stat
    real :: itheta

    type(scalar_field) :: sfield, old_sfield, nl_sfield
    type(vector_field) :: vfield, old_vfield, nl_vfield
    type(tensor_field) :: tfield, old_tfield, nl_tfield

    !For projecting velocity to continuous 
    type(vector_field) :: U_nl, pvelocity, X
    type(vector_field), pointer :: velocity

    do s=1,size(state)

       velocity=>extract_vector_field(state(s), "Velocity", stat)
       if(stat==0) then
          if (have_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/relaxation")) then
             call get_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/relaxation", itheta, default=0.5)
          else if (have_option(trim(velocity%option_path)//"/prescribed/temporal_discretisation/relaxation")) then
             call get_option(trim(velocity%option_path)//"/prescribed/temporal_discretisation/relaxation", itheta, default=0.5)
          end if
       else
          itheta = 0.5
       end if

       do f=1,scalar_field_count(state(s))

          sfield=extract_scalar_field(state(s), f)

          if(.not.aliased(sfield)) then

             old_sfield=extract_scalar_field(state(s), "Old"//trim(sfield%name),&
                  & stat=old_stat)
             
             nl_sfield=extract_scalar_field(state(s), "Nonlinear"//trim(sfield%name),&
                  & stat=nl_stat)

             if ((old_stat==0).and.(nl_stat==0)) then
                call set(nl_sfield, sfield, old_sfield, itheta)
             end if

          end if

       end do

       do f=1,vector_field_count(state(s))

          vfield=extract_vector_field(state(s), f)

          if(.not.aliased(vfield)) then

             old_vfield=extract_vector_field(state(s), "Old"//trim(vfield%name),&
                  & stat=old_stat)
             
             nl_vfield=extract_vector_field(state(s), "Nonlinear"//trim(vfield%name),&
                  & stat=nl_stat)
                
             if ((old_stat==0).and.(nl_stat==0)) then
                call set(nl_vfield, vfield, old_vfield, itheta)
             end if
             
          end if

       end do

       do f=1,tensor_field_count(state(s))

          tfield=extract_tensor_field(state(s), f)

          if(.not.aliased(tfield)) then

             old_tfield=extract_tensor_field(state(s), "Old"//trim(tfield%name),&
                  & stat=old_stat)

             nl_tfield=extract_tensor_field(state(s), "Nonlinear"//trim(tfield%name),&
                  & stat=nl_stat)

             if ((old_stat==0).and.(nl_stat==0)) then
                call set(nl_tfield, tfield, old_tfield, itheta)
             end if

          end if

       end do

    end do

    !Compute velocity field projected to continuous for DG advection
    !Not currently coded for multimaterial/phase as I don't understand
    !how it works
    if(has_vector_field(state(1),"ProjectedNonlinearVelocity")) then
       U_nl = extract_vector_field(state(1),"NonlinearVelocity")
       pvelocity = extract_vector_field(state(1),&
            "ProjectedNonlinearVelocity")
       X = extract_vector_field(state(1),"Coordinate")
       call project_field(U_nl,pvelocity,X)
    end if

  end subroutine relax_to_nonlinear

  function simulation_completed(current_time, timestep)
    !!< Simulation end test routine. Tests standard timestep loop exit
    !!< conditions (many listed under /timestepping). Returns .true. if these
    !!< conditions are satisfied and .false. otherwise.

    real, intent(in) :: current_time
    integer, intent(in), optional :: timestep

    logical :: simulation_completed

    integer :: final_timestep, i, stat
    real :: current_cpu_time, time_limit, current_wall_time

    simulation_completed = .false.

    do i = 1, 5
       select case(i)
       case(1)
          call get_option("/timestepping/finish_time", time_limit)
          if(current_time >= time_limit) then
             simulation_completed = .true.
             ewrite(1, *) "Finish time reached"
             exit
          end if
       case(2)
          if(present(timestep)) then
             call get_option("/timestepping/final_timestep", final_timestep, stat)
             if(stat == SPUD_NO_ERROR) then
                if(timestep > final_timestep) then
                   simulation_completed = .true.
                   ewrite(1, *) "Passed final timestep"
                   exit
                end if
             end if
          end if
       case(3)
          call get_option("/timestepping/cpu_time_limit", time_limit, stat)
          if(stat == SPUD_NO_ERROR) then
             call cpu_time(current_cpu_time)
             call allmax(current_cpu_time)
             if(current_cpu_time - simulation_start_cpu_time >= time_limit) then
                simulation_completed = .true.
                ewrite(1, *) "CPU time limit reached"
                exit
             end if
          end if
       case(4)
          call get_option("/timestepping/wall_time_limit", time_limit, stat)
          if(stat == SPUD_NO_ERROR) then
             current_wall_time = wall_time()
             call allmax(current_wall_time)
             if(current_wall_time - simulation_start_wall_time >= time_limit) then
                simulation_completed = .true.
                ewrite(1, *) "Wall time limit reached"
                exit
             end if
          end if
       case(5)
          if(SIG_INT) then
             simulation_completed = .true.
             ewrite(1, *) "Interrupt signal received"
             exit
          end if
       case default
          FLAbort("Invalid loop index")
       end select
    end do

    if(simulation_completed) then
       ewrite(2, *) "simulation_completed returning .true."
    else
       ewrite(2, *) "simulation_completed returning .false."
    end if

  end function simulation_completed

end module timeloop_utilities
