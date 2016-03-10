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
  module Adjoint_ROM_main_module

  use spud
  use fields
  use state_module
  use write_state_module
  use timeloop_utilities
  use global_parameters, only: option_path_len, current_time, dt, FIELD_NAME_LEN ,&
                               timestep, OPTION_PATH_LEN, &
                               simulation_start_time, &
                               simulation_start_cpu_time, &
                               simulation_start_wall_time, &
                               topology_mesh_name
  use FLDebug
  use snapsvd_module
  use vtk_interfaces
  use memory_diagnostics
  use populate_state_module
  use python_state
  use vtk_interfaces
  use quadrature
  use diagnostic_fields_wrapper
  use field_options
  use signals
  use spud
  use tictoc
  use fldebug
  use fluids_module
  !use reduced_fluids_module
  use signals
  use spud
  use tictoc
!#ifdef HAVE_ZOLTAN
  use zoltan
  use nonlinear_conjugate_gradient
  use momentum_equation
  use momentum_equation_reduced_adjoint
  use momentum_equation_reduced
  use reduced_model_runtime

!#endif
  implicit none
#ifdef HAVE_PETSC
#include "petsc_legacy.h"
#endif 


  private


  public :: Adjoint_ROM_main

contains

  subroutine Adjoint_ROM_main
    type(state_type), dimension(:), pointer :: state=> null()
    ! type(state_type), dimension(:), allocatable :: state_adj
    type(state_type), dimension(:,:,:), allocatable :: POD_state
    
    integer :: ierr,stat,i,j,k
    character(len = OPTION_PATH_LEN) :: simulation_name
    integer :: gloits,nocva,linits,nsvd
    real :: valfun,fszero
    real, dimension(:), allocatable :: cva,cvaold,g,gold,d
    logical :: linconver
    INTEGER,PARAMETER ::NGLOITS=1
    REAL,PARAMETER ::GOLDR = 5.0
    integer :: total_timestep,total_dump_no
    type(vector_field), pointer :: u ,POD_velocity
    type(scalar_field), pointer :: POD_pressure
    real, dimension(:), allocatable :: pod_coef_optimised
    real  ::current_time,finish_time,real_dump_period
    logical :: if_optimal
    !  type(state_type), dimension(:,:,:) :: POD_state


    !----------------------------------
    ! get State
    !----------------------------------
    call initialise_write_state
    ! Read state from .flml file
    call populate_state(state)
    !----------------------------------
     
    !-----------
    fszero=0.0
    valfun=0.0
    linconver= .true.
    if_optimal = .false.
    !----------------------------------
    !call forward ROM
    !get the initial ROM solution coef.
    !----------------------------------


    call get_option("/timestepping/current_time", current_time)! save the current time

    call delete_option("/reduced_model/adjoint")  ! switch to forward ROM
    ewrite(1,*) 'Forward ROM', have_option("/reduced_model/adjoint")

!    call fluids()
    ewrite(1,*) 'Completed the initial forward running'
!    stop 11
    !------------------------------------------------------------------------
    ! allocate variables:
    ! nocva: number of control variables to be optimised
    ! cva: control variables to be optimised, here initial conditions
    ! cvaold: control variables at the previous iteration
    ! g: the gradient calculated from the adjoint ROM
    ! gold: the gradient calculated from the previous nonlinear CG iteration
    ! pod_coef_optimised: obsearvation POD solutions at all the time levelsreduced_model
    ! d: linea search direction
    !-------------------------------------------------------------------------
    call get_option(&
       '/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
    call get_option('/simulation_name', simulation_name)

    call read_pod_basis_differntmesh(POD_state, state)
    u => extract_vector_field(POD_state(1,1,1), "Velocity")
    nocva=(u%dim+1)*size(POD_state,1)
    allocate(cva(nocva))
    allocate(pod_coef_optimised(nocva))
    allocate(d(nocva))
    allocate(cvaold(nocva))
    allocate(gold(nocva))
    allocate(g(nocva))
    gold=0
    d=0
    !=============================================================================================
    ! NONLINEAR CG ITERATION
    !=============================================================================================
    global_loop:DO gloits = 1,NGLOITS
       LINITS =1 
       !get g
          open(10,file='cost_function1',ACTION='WRITE')
          write(10,*) valfun
          close(10)                 
       call add_option("/reduced_model/adjoint",stat) ! switch to adjoint ROM
       call set_option("/timestepping/current_time", current_time) ! set back the original current time
       ewrite(1,*) 'Adjoint ROM', have_option("/reduced_model/adjoint")
       call fluids() !get the func run simulation_name_POD.flml run forward ROM 
       call set_option("/timestepping/current_time", current_time) ! set back the original current time
     !  stop 12
       open(unit=2,file='adjoint_g')
       read(2,*)(g(i),i=1, nocva)
       close(2)
       !call solve_momentum_reduced_adjoint(state, at_first_timestep, timestep, POD_state, POD_state_deim, snapmean, eps,g,its)
stop 75       
       !-----------------------------------------------------------------------------
       ! Line Search
       !-----------------------------------------------------------------------------
       call delete_option("/reduced_model/adjoint") ! switch to forward ROM
       linear_search_loop: do   !need control the loop
          call set_option("/timestepping/current_time", current_time) ! set back the original current time
          ewrite(1,*) 'Forward ROM' 
          if(linits.eq.1) then
             open(40,file='coef_pod_all0')
             read(40,*)(cva(i),i=1,size(cva))
             close(40)
          endif
          call fluids(if_optimal=if_optimal,valfun=valfun) !get the func run simulation_name_POD.flml run forward ROM 
          call set_option("/timestepping/current_time", current_time) ! set back the original current time
          valfun=value_cost_function(nocva, pod_state,state)
          ewrite(1,*) 'Cost Function = ', valfun
          call  NONLINCG(gloits,linits,valfun,cva,g,goldr,linconver,nocva,gold,d,cvaold,fszero)
          pod_coef_optimised = cva  !the variables need to be modified
          !print  costfunction 
          if(LINITS==1) then
             open(10,file='cost_function1',ACTION='WRITE')
          else
             open(10,file='cost_function1',ACTION='WRITE',position='append')
          endif
          write(10,*) LINITS, valfun
          close(10) 
          if(linconver) exit
      !    LINITS = LINITS + 1 
          ! upadate the initial coeffient
          open(40,file='coef_pod_all0')
          write(40,*)(cva(i),i=1,size(cva))
          close(40)
        enddo linear_search_loop
       !-----------------------------------------------------------------------------
       ! END Line Search
       !-----------------------------------------------------------------------------
       
       ! print adjoint g 
       ! print  costfunction
       if(gloits==1) then
          open(10,file='cost_function',ACTION='WRITE')
       else
          open(10,file='cost_function',ACTION='WRITE',position='append')
       endif
       write(10,*) LINITS, valfun
       close(10) 
    enddo global_loop
    call  fluids(if_optimal=if_optimal,valfun=valfun)
    
    do i = 1, size(state)
       call deallocate(state(i))
    end do

   if (allocated(pod_state)) then
       do i=1, size(pod_state,1)
          do j=1,size(pod_state,2)
             do k=1,size(pod_state,3)
                call deallocate(pod_state(i,j,k))
             enddo
          enddo
       end do
    end if

  end subroutine Adjoint_ROM_main

  function value_cost_function(nocva, pod_state,state) result (FUNCT)
    !use fluids_module
    !use reduced_fluids_module
    !use signals
    !use spud
    !use reduced_model_runtime
    
    type(state_type), dimension(:,:,:), intent(in) :: POD_state
    type(state_type), dimension(:), intent(in)  :: state
    real :: FUNCT
    integer, intent(in) :: nocva
    integer :: i,d,j,k,unodes,ii,kk
    type(vector_field), pointer :: u_obv
    real, dimension(:,:), allocatable :: pod_coef_all_obv,pod_coef_all_comp,ifhaveobs,u_comp   !might up 
    
    integer POD_num, total_timestep
    real finish_time
    type(scalar_field), pointer :: POD_pressure
    type(vector_field), pointer :: POD_velocity
    character(len = OPTION_PATH_LEN) :: simulation_name

    call get_option('/reduced_model/pod_basis_formation/pod_basis_count', POD_num)  
    FUNCT = 0.0
    
    !  call read_pod_basis_differntmesh(POD_state, state) !! if adaptive meshes, it may take a lot of CPU, fix it later
    
    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time",  finish_time)
    call get_option("/timestepping/timestep", dt)       
    total_timestep=int((finish_time-current_time)/dt)

 
    !  u => extract_vector_field(state(1), "Velocity", stat)
    u_obv=> extract_vector_field(state(1), "Velocity")


    call zero(u_obv)
    unodes=node_count(u_obv)
    allocate(pod_coef_all_obv(0:total_timestep, nocva) )
    allocate(pod_coef_all_comp(0:total_timestep,nocva) )
    allocate(ifhaveobs(total_timestep,unodes))
    allocate(u_comp(u_obv%dim,unodes))
    pod_coef_all_obv=0.0
    pod_coef_all_comp=0.0
    u_comp=0.0
    ifhaveobs=1
    open(30,file='coef_pod_all_obv')
    open(100,file='coef_pod_all')
    DO k = 0, total_timestep
       read(30,*)(pod_coef_all_obv(k,ii),ii=1,nocva)
       read(100,*)(pod_coef_all_comp(k,ii),ii=1,nocva)
    ENDDO

    close(30)
    close(100)

    DO k = 2, total_timestep
!    DO k = total_timestep-1, total_timestep-1
       print*,'time=',k
       if(2*int(k/2).eq.k) then
          !IF(ifhaveobs(k,:).eq.1) THEN   ! if has obv, iflagobs=1, else=0.
          u_obv%val = 0.0
          u_comp = 0.0
          do i=1,POD_num
             POD_velocity=>extract_vector_field(POD_state(i,1,1), "Velocity")
             POD_pressure=>extract_scalar_field(POD_state(i,2,1), "Pressure") 
             do d=1, u_obv%dim
                u_obv%val(d,:) =pod_coef_all_obv(k,i+(d-1)*POD_num)*POD_velocity%val(d,:)
                u_comp(d,:)=pod_coef_all_comp(k,i+(d-1)*POD_num)*POD_velocity%val(d,:)
             enddo
          enddo
          do d =1, u_obv%dim
             do j =1,unodes
                FUNCT =FUNCT+0.5*(u_comp(d,j)-u_obv%val(d,j))**2
             enddo
          enddo
          print*,'FUNCT=',k,FUNCT
       endif
       !ENDIF  ! IF(ifhaveobs
    ENDDO  ! end timestep 
    print*,FUNCT
!    stop 78
    
    deallocate(pod_coef_all_obv)
    deallocate(pod_coef_all_comp)
    deallocate(ifhaveobs)
    deallocate(u_comp)
  end function value_cost_function
end module Adjoint_ROM_main_module
