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
  module Nonintrusive_rom_module

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
#include "finclude/petsc.h"
#endif 


  private


  public :: Non_Intrusive_ROM_main

contains

  subroutine Non_Intrusive_ROM_main
    type(state_type), dimension(:), pointer :: state=> null()
    ! type(state_type), dimension(:), allocatable :: state_adj
    type(state_type), dimension(:,:,:), allocatable :: POD_state
    
    integer :: ierr,stat,i,j,k
    character(len = OPTION_PATH_LEN) :: simulation_name
    integer :: nsvd
    real :: valfun,fszero
    integer :: total_timestep,total_dump_no
    type(vector_field), pointer :: snapmean_velocity
    type(scalar_field), pointer :: snapmean_pressure
    type(vector_field), pointer :: u ,POD_velocity
    type(scalar_field), pointer :: p,POD_pressure
    real, dimension(:), allocatable :: pod_coef_optimised
    logical :: if_optimal,if_timetwo
    !  type(state_type), dimension(:,:,:) :: POD_state
    ! Change in pressure
    type(scalar_field) :: delta_p
    ! Change in velocity
    type(vector_field) :: delta_u
    real, dimension(:,:), allocatable :: pod_sol_velocity
    real, dimension(:), allocatable :: pod_sol_pressure

    ! non_intrusive_ROM
    type(state_type), dimension(:,:,:), allocatable :: Non_Intrusive_state
    real, dimension(:), allocatable ::  pod_coef_orig,pod_coef_new,pod_coef_old
    real, dimension(:,:), allocatable ::  DM0,DM1
    real, dimension(:,:,:), allocatable ::  DDM
    integer istate,NPERT,ntime
    integer dump_no_tmp, int_dump_period
    real current_time,finish_time,dt

    call delete_option("/reduced_model/execute_reduced_model") ! switch to full model
    !----------------------------------
    ! get State
    !----------------------------------
    call initialise_write_state
    ! Read state from .flml file
    call populate_state(state)
    u => extract_vector_field(state(1), "Velocity", stat)
    p => extract_scalar_field(state(1), "Pressure")

    !----------------------------------
    !----------------------------------

    call read_pod_basis_differntmesh(POD_state, state)
    allocate(pod_coef_orig((u%dim+1)*size(POD_state,1)))
    allocate(pod_coef_new((u%dim+1)*size(POD_state,1)))
    allocate(pod_coef_old((u%dim+1)*size(POD_state,1)))
    NPERT = (u%dim+1)*size(POD_state,1)
    allocate(DM0(NPERT,NPERT))
    allocate(DM1(NPERT,NPERT))
    allocate(DDM(NPERT,NPERT,NPERT))  

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time", finish_time)       
    call get_option("/timestepping/timestep", dt)
    ntime=int((finish_time-current_time)/dt)

    !----------------------------------
    !get the initial ROM solution coef.
    !----------------------------------
    istate = 1
    call project_from_full_to_pod(istate, pod_state, state, pod_coef_orig)
    !----------------------------------------------
    !get the first order and second order matrices.
    !----------------------------------------------   
if(.true.) then
    call First_Order_DM(state,pod_state,pod_coef_orig,DM0) 
    open(10,file='DM.dat')
    do k=1,NPERT
       write(10,*) (DM0(i,k),i=1,NPERT)
    enddo
    close(10)
    DM1 = DM0
    call Hessian_DDM(state,pod_state,pod_coef_orig,DM1,DDM) 
    open(10,file='DDM.dat')
    do k=1,NPERT
       do j=1,NPERT
          write(10,*) (DDM(i,j,k),i=1,NPERT)
       enddo
    enddo
    close(10)
else
    open(10,file='DM.dat')
    do k=1,NPERT
       read(10,*) (DM0(i,k),i=1,NPERT)
    enddo
    close(10)

    open(10,file='DDM.dat')
    do k=1,NPERT
       do j=1,NPERT
          read(10,*) (DDM(i,j,k),i=1,NPERT)
       enddo
    enddo
    close(10)
endif
    !--------------------------------------------------------------
    ! take 2 time steps as we have an issue with initialization...
    !--------------------------------------------------------------
    if_timetwo=.true.
    call First_Order_DM(state,pod_state,pod_coef_orig,DM0, if_timetwo=if_timetwo)
    pod_coef_old = pod_coef_orig
    open(10,file='nonintrusive_coef_pod2.dat')
    read(10,*) (pod_coef_new(i),i=1,(u%dim+1)*size(POD_state,1))
    close(10)   

    call Nonintrusive(NPERT,ntime,DM0,DDM,pod_coef_orig,pod_coef_old,pod_coef_new)

             
    ! Allocate the change in pressure field
    call allocate(delta_p, p%mesh, "DeltaP")
    delta_p%option_path = trim(p%option_path)
    call zero(delta_p)
    
    ! Allocate the change in velocity pressure field
    call allocate(delta_u, u%dim, u%mesh, "DeltaU")
    delta_u%option_path = trim(u%option_path)
    call zero(delta_u)


    open(20,file='rom_solution_alltime.dat')

    if(have_option("/io/dump_period_in_timesteps/constant")) then
       call get_option("/io/dump_period_in_timesteps/constant", int_dump_period)
    else
       int_dump_period = 1
    endif

    do k = 0, ntime,int_dump_period
       print*,'**********************************',ntime,k
        read(20,*)(pod_coef_new(i),i=1,NPERT)
       call project_full(delta_u, delta_p, pod_sol_velocity, pod_sol_pressure, POD_state(:,:,istate), pod_coef_new)
       snapmean_velocity=>extract_vector_field(POD_state(1,1,istate),"SnapmeanVelocity")
       snapmean_pressure=>extract_Scalar_field(POD_state(1,2,istate),"SnapmeanPressure")
       !call addto(u, delta_u, dt)
       u%val=snapmean_velocity%val
       call addto(u, delta_u)
       p%val=snapmean_pressure%val
       call addto(p, delta_p)
       dump_no_tmp= int(k/int_dump_period)
       call write_state(dump_no_tmp, state)
    enddo
    close(20)
    !----------------------------------
    ! deallocate variables
    !----------------------------------
   if (allocated(pod_state)) then
       do i=1, size(pod_state,1)
          do j=1,size(pod_state,2)
             do k=1,size(pod_state,3)
                call deallocate(pod_state(i,j,k))
             enddo
          enddo
       end do
    end if

    deallocate(pod_coef_orig)
    deallocate(pod_coef_new)
    deallocate(pod_coef_old)
  end subroutine Non_Intrusive_ROM_main


  SUBROUTINE Nonintrusive(NPERT,Ntime,DM0,DDM,TOLDKEEP_ORIG,told,tnew)
    IMPLICIT NONE
    INTEGER NPERT,NTIME
    REAL TNEW(NPERT),TOLD(NPERT),TOLDold(NPERT), TOLDkeep(NPERT)
    REAL TOLDKEEP_ORIG(NPERT)
    real, dimension(:,:), allocatable ::  M,MTM,MTM_INV,MTM_INV_M,M2,M1,M0,m00
    REAL DF(NPERT),DdF(NPERT),MTM_INV_F(NPERT)
    INTEGER ITIME,I,IPERT,JPERT, ISTART, IEND,kpert,k,ivar
    REAL TOLER
    real, dimension(NPERT,NPERT)::  DM0
    real, dimension(:,:), allocatable ::  DM1,DM2
    real, dimension(NPERT,NPERT,NPERT) ::  DDM
    real toler1,toler2
    real j
    
    allocate(DM2(NPERT,NPERT))
    
    ! ****the first two time steps...
    
    toldold=told
    told=tnew
    open(20,file='rom_solution_alltime.dat')
    write(20,*)(TOLDOLD(i),i=1,NPERT)
    write(20,*)(TNEW(i),i=1,NPERT)

    ! ****Perform the time stepping next...
    DO Itime=2,ntime
       tnew=told            
       do i=1,NPERT
          tnew(i)=tnew(i) + sum( DM0(i,:)*(told(:)-toldold(:)) )
       enddo
       
       do kpert=1,NPERT
          do i=1,NPERT
             tnew(i)=tnew(i) +(toldold(kpert)-TOLDKEEP_ORIG(kpert)  +   0.5*(told(kpert)-toldold(kpert)))  &
                  *sum( DDM(i,:,kpert)*(told(:)-toldold(:)) )
          end do
       end do
       
       toldold=told
       TOLD=tnew  

       write(20,*)(tnew(i),i=1,NPERT)

    END DO
    
    close(20)
    deallocate(DM2)
    
  END SUBROUTINE NONINTRUSIVE

  SUBROUTINE First_Order_DM(state,pod_state,TOLDKEEP_ORIG,DM0,if_timetwo) 
    type(state_type), dimension(:), intent(in) :: state 
    type(state_type), dimension(:,:,:), intent(in) :: pod_state
    INTEGER NONODS,NTIME,NPERT
    REAL,dimension(:), intent(inout) ::TOLDKEEP_ORIG
    real toler
    PARAMETER(TOLER=0.00001)
    REAL,dimension(:,:), intent(inout) :: DM0
    logical, optional :: if_timetwo
    
    ! local variables...
    integer ipert,istate,i,j,k
    logical if_optimal
    real, DIMENSION( : ), allocatable :: told,tnew
    real, DIMENSION( :, : ), allocatable :: m0,m1

    NPERT=size(TOLDKEEP_ORIG)
    allocate(told(NPERT))
    allocate(tnew(NPERT))
    allocate(M0(NPERT,NPERT))
    allocate(M1(NPERT,NPERT))
    

    TOLD=TOLDKEEP_ORIG ! no perturbations - a bit wasteful of CPU
    open(10,file='nonintrusive_coef_pod0.dat')
    write(10,*) (TOLD(i),i=1,NPERT)
    close(10)
    !------run fluids------
    call fluids()
    !------end run fluids------
    
    open(10,file='nonintrusive_coef_pod1.dat')
    read(10,*) (TNEW(i),i=1,NPERT)
    close(10)

    if(present(if_timetwo)) then
       open(10,file='nonintrusive_coef_pod2.dat')
       write(10,*) (TNEW(i),i=1,NPERT)
       close(10)   
       return
    endif
         
    DO IPERT=1,NPERT
       M0(:,IPERT) = tnew(:)
    END DO
    
    DO IPERT=1,NPERT            
       TOLD=TOLDKEEP_ORIG
       TOLD(IPERT)=TOLD(IPERT)+toler         
       open(10,file='nonintrusive_coef_pod0.dat')
       write(10,*) (TOLD(i),i=1,NPERT)
       close(10)
       !------run fluids------
       call fluids()
       !------end run fluids------
       open(10,file='nonintrusive_coef_pod1.dat')
       read(10,*) (TNEW(i),i=1,NPERT)
       close(10)
       M1(:,IPERT) = tnew(:)
    END DO
    
    DM0=(M1-M0)/toler
    
    deallocate(told)
    deallocate(tnew)
    deallocate(M0)
    deallocate(M1)

  END SUBROUTINE FIRST_ORDER_DM
  
  SUBROUTINE Hessian_DDM(state,pod_state,TOLDKEEP_ORIG,DM1,DDM) 
    type(state_type), dimension(:), intent(inout) :: state 
    type(state_type), dimension(:,:,:), intent(inout) :: pod_state
    INTEGER NONODS,NTIME,NPERT
    REAL TOLER
    REAL,dimension(:), intent(in) ::TOLDKEEP_ORIG
    REAL,dimension(:,:), intent(inout) ::DM1
    REAL,dimension(:,:,:), intent(inout) ::DDM
    PARAMETER(TOLER=0.00001)
    ! local variables...
    integer ipert,istate
    logical if_optimal
    real, DIMENSION( : ), allocatable :: told,tnew,TOLDkeep
    real, DIMENSION( :, : ), allocatable :: m0,m1
    real, dimension(:,:), allocatable ::  DM0,DM2
    real toler1,toler2
    integer KPERT,i,j,k
    
    ! DDM (Hessian) calculation: 
    NPERT=size(TOLDKEEP_ORIG)
    allocate(DM2(NPERT,NPERT))
    allocate(TOLDkeep(NPERT))
    
    DO KPERT=1,NPERT
       print*,'*************************KPERT111111111111111',KPERT
       DM1=0.0      
       TOLDkeep=TOLDKEEP_ORIG
       TOLDkeep(KPERT)=TOLDkeep(KPERT)-toler
       TOLD=TOLDkeep
       
       call First_Order_DM(state,pod_state,TOLD,DM1) 

       DM2=0.0          
       TOLDkeep=TOLDKEEP_ORIG
       TOLDkeep(KPERT)=TOLDkeep(KPERT)+toler           
       TOLD=TOLDkeep
       
       print*,'*************************KPERT2222222222222222222222222',KPERT
       call First_Order_DM(state,pod_state,TOLD,DM2) 
       DDM(:,:,kpert)=(DM2(:,:)-DM1(:,:))/(2*TOLER)
       print*,'DDM(:,:,kpert)',DDM(:,:,kpert)
    END DO
    close(10)

    deallocate(DM2)
    deallocate(TOLDkeep)
    
  END SUBROUTINE Hessian_DDM

  
END module Nonintrusive_rom_module

