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
program form_pod_basis
  use spud
  use fields
  use state_module
  use write_state_module
  use timeloop_utilities
  use global_parameters, only: option_path_len, current_time, dt, FIELD_NAME_LEN
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
  use momentum_equation
  use solvers 
  use momentum_cg 
  use momentum_equation_reduced
  !  use checkpoint

  implicit none
#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif
  integer, dimension(:), allocatable :: indices,indices_tmp
  type(state_type), dimension(:,:), allocatable :: state ,state_p
  type(state_type), dimension(:), allocatable ::  state_test,state_deim 
  type(state_type), dimension(:,:,:), allocatable :: pod_state, pod_deim_state
  type(state_type), dimension(:,:,:), allocatable :: pod_state_p
  type(vector_field) , pointer :: position_deim
  integer :: timestep
  integer :: ierr
  integer :: numphase
  character(len = OPTION_PATH_LEN) :: simulation_name

#ifdef HAVE_MPI
  call mpi_init(ierr)
#endif

#ifdef HAVE_PETSC
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

  call python_init()
  call read_command_line()

 ! call form_basis_different_mesh()
  call form_basis()

  call deallocate(state)
  call deallocate(state_p)
  call deallocate_transform_cache()

  call print_references(0)
#ifdef HAVE_MEMORY_STATS
  call print_current_memory_stats(0)
#endif

#ifdef HAVE_MPI
  call mpi_finalize(ierr)
#endif

contains

  subroutine form_basis()
    !!< Matrices containing the snapshots for arpack
    !      type(state_type), dimension(:), allocatable :: state
    real, dimension(:,:,:,:), allocatable :: snapmatrix_velocity, snapmatrix_velocity_deim
    real, dimension(:,:,:), allocatable :: snapmatrix_pressure, snapmatrix_pressure_deim
    real, dimension(:,:,:,:), allocatable :: leftsvd_velocity, leftsvd_velocity_deim
    real, dimension(:,:,:), allocatable :: leftsvd_pressure, leftsvd_pressure_deim
    real, dimension(:,:,:), allocatable :: svdval_velocity, svdval_velocity_deim
    real, dimension(:,:), allocatable :: svdval_pressure, svdval_pressure_deim
    real, dimension(:,:,:), allocatable :: snapmean_velocity, snapmean_velocity_deim
    real, dimension(:,:), allocatable :: snapmean_pressure, snapmean_pressure_deim
    integer :: snapshots, u_nodes, p_nodes, nsvd
    integer :: i,dump_no, d, dim,j,k ,m,i2,i3
    integer :: stat
    
    real, dimension(:,:), allocatable :: P !n*m (u_nodes*nsvd)
    !integer, dimension(:), allocatable :: phi
    real, dimension(:,:), allocatable :: leftsvd_velocity_deim_uin
    integer :: udim,deim_number_tmp  
    type(vector_field), pointer :: velocityudim 
    type(vector_field), pointer :: vfield
    character(len=1024) :: phase_name
    call get_option(&
         '/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
    call get_option('/simulation_name',simulation_name)
    numphase=option_count("/material_phase")
    call read_input_states_fromvtu(state,state_p) 

    ! do  p=1,numphase
      call retrieve_snapshots(state, snapshots, u_nodes, p_nodes, snapmatrix_velocity, snapmatrix_pressure, &
                             & snapmean_velocity, snapmean_pressure)    
      call form_svd(snapmatrix_velocity, snapmatrix_pressure,&
          & leftsvd_velocity, leftsvd_pressure, svdval_velocity, svdval_pressure, snapshots)
        
       call form_podstate_velocity(state,pod_state,leftsvd_velocity,leftsvd_pressure, snapmean_velocity, snapmean_pressure)
       call form_podstate_pressure(state_p,pod_state_p,leftsvd_velocity,leftsvd_pressure, snapmean_velocity, snapmean_pressure)
    ! enddo
   
    do m=1, numphase
           call get_option("/material_phase["//int2str(m-1)//"]/name", phase_name)
    do i=1,nsvd 
        dump_no=i 
       !  call vtk_write_pod_state_velocity(filename=trim(simulation_name)//"_PODBasisV", index=dump_no, state=pod_state(i,:))
       !  call vtk_write_pod_state_pressure(filename=trim(simulation_name)//"_PODBasisP", index=dump_no, state=pod_state(i,:))
         call vtk_write_state(filename=trim(simulation_name)//"_"//trim(phase_name)//"_PODBasisVelocity", index=dump_no, model="VelocityMesh", state=pod_state(i,:,m))
         call vtk_write_state(filename=trim(simulation_name)//"_"//trim(phase_name)//"_PODBasisPressure", index=dump_no, model="PressureMesh", state=pod_state_p(i,:,m))
       ! call checkpoint_simulation(pod_state(i,:),prefix="gyre",cp_no=dump_no, postfix="podbasis")
       ! call checkpoint_fields(pod_state(i,:),prefix="gyre",postfix="podbasis",cp_no=dump_no)
       ! call checkpoint_state(pod_state(i,:),prefix="gyre",postfix="podbasis",cp_no=dump_no)
       ! call deallocate(pod_state(i,:,:))
        !call deallocate(pod_state_p(i,:,:))
    enddo
    enddo !numphase

    !! Produce updated flml file in which the execute_reduced_model
    !! option is set.
    call add_option("/reduced_model/execute_reduced_model",stat)
    call set_option('/simulation_name', trim(simulation_name)//'_POD') 
    call write_options(trim(simulation_name)//"_POD.flml")
     !do i=1,nsvd

     !  dump_no=i

     !  call vtk_write_state(filename=trim(simulation_name)//"_PODBasis", index=dump_no, state=pod_state(i,:))

      ! call deallocate(pod_state(i,:))

 !   enddo

    !! Produce updated flml file in which the execute_reduced_model
    !! option is set.
   ! call add_option("/reduced_model/execute_reduced_model",stat)
   ! call set_option('/simulation_name', trim(simulation_name)//'_POD') 
!    call set_option('/timestepping/nonlinear_iterations',1)
   ! call write_options(trim(simulation_name)//"_POD.flml")
   deim=have_option("/reduced_model/discrete_empirical_interpolation_method")
   if(deim) then
       call get_option(&
         '/reduced_model/pod_basis_formation/pod_basis_count', deim_number)
 
     end if
    deallocate(snapmatrix_velocity)
     deallocate(snapmatrix_pressure)
    deallocate(snapmean_velocity)
    deallocate(snapmean_pressure)  
    deallocate(pod_state)
     deallocate(pod_state_p)  
  end subroutine form_basis

   

  subroutine retrieve_snapshots(state, snapshots, u_nodes, p_nodes, snapmatrix_velocity, snapmatrix_pressure, &
                                & snapmean_velocity, snapmean_pressure)
   !!< 
    type(state_type), intent(in), dimension(:,:) :: state
    real, dimension(:,:,:,:), allocatable :: snapmatrix_velocity
    real, dimension(:,:,:), allocatable :: snapmatrix_pressure
    real, dimension(:,:,:), allocatable :: snapmean_velocity
    real, dimension(:,:), allocatable :: snapmean_pressure
    ! integer :: numphase
    integer :: snapshots, u_nodes, p_nodes, dim, i, d, p
    type(vector_field), pointer :: velocity
    type(scalar_field), pointer :: pressure
    
    
     velocity => extract_vector_field(state(1,1), "Velocity")
     pressure => extract_scalar_field(state_p(1,1), "Pressure")
  
    dim=velocity%dim
    p_nodes=node_count(pressure)
    u_nodes=node_count(velocity)
    snapshots=size(state(:,1))
    allocate(snapmatrix_velocity(u_nodes,snapshots,dim,numphase))
    allocate(snapmatrix_pressure(p_nodes,snapshots,numphase))
    allocate(snapmean_velocity(u_nodes,dim,numphase))
    allocate(snapmean_pressure(p_nodes,numphase))   
 
    snapmatrix_velocity=0.0
    snapmatrix_pressure=0.0
    snapmean_velocity=0.0
    snapmean_pressure=0.0
   do p = 1, numphase
    do i = 1, snapshots
       velocity => extract_vector_field(state(i,p), "Velocity")
      !call print_state(state_p(i,p))
      ! print *, 'ppppppp=', p,i
      ! pressure => extract_scalar_field(state_p(i,p), "Pressure")
       pressure => extract_scalar_field(state_p(i,1), "Pressure")
       do d = 1, dim 
          ! print *,'size(snapmatrix_velocity,1), size(velocity%val,2)',size(snapmatrix_velocity,1), size(velocity%val,2) ,size(velocity%val,1)   
          ! print *,velocity%val(d,:)     
           snapmatrix_velocity(:,i,d,p)=velocity%val(d,:)
       end do 
           snapmatrix_pressure(:,i,p)=pressure%val
    end do

    do i=1, snapshots

       do d=1, dim
          snapmean_velocity(:,d,p)= snapmean_velocity(:,d,p)+snapmatrix_velocity(:,i,d,p)
       enddo
       snapmean_pressure(:,p)=snapmean_pressure(:,p)+snapmatrix_pressure(:,i,p)
    end do

    do d=1,dim
       snapmean_velocity(:,d,p)=snapmean_velocity(:,d,p)/snapshots
    enddo
       snapmean_pressure(:,p)=snapmean_pressure(:,p)/snapshots

    do i=1,snapshots
       do d=1,dim
          snapmatrix_velocity(:,i,d,p)=snapmatrix_velocity(:,i,d,p)-snapmean_velocity(:,d,p)
       enddo
       snapmatrix_pressure(:,i,p)=snapmatrix_pressure(:,i,p)-snapmean_pressure(:,p)
    enddo
    enddo !end numphase
    
  end subroutine retrieve_snapshots


  subroutine form_svd(snapmatrix_velocity, snapmatrix_pressure,&
       & leftsvd_velocity, leftsvd_pressure, svdval_velocity, svdval_pressure, snapshots)
    
    real, dimension(:,:,:,:), intent(in) :: snapmatrix_velocity
    real, dimension(:,:,:), intent(in) :: snapmatrix_pressure
    real, dimension(:,:,:,:), allocatable, intent(out) :: leftsvd_velocity
    real, dimension(:,:,:), allocatable, intent(out) :: leftsvd_pressure
    real, dimension(:,:,:), allocatable, intent(out) :: svdval_velocity
    real, dimension(:,:), allocatable, intent(out) :: svdval_pressure
    integer i, d, dim ,nsvd, snapshots, p_nodes, u_nodes, p

    call get_option(&
         '/reduced_model/pod_basis_formation/pod_basis_count', nsvd)

    dim=size(snapmatrix_velocity,3)
    p_nodes=size(snapmatrix_pressure,1)
    u_nodes=size(snapmatrix_velocity,1)

    allocate(leftsvd_velocity(u_nodes,nsvd,dim,numphase))
    allocate(leftsvd_pressure(p_nodes,nsvd,numphase))
    allocate(svdval_velocity(nsvd,dim,numphase))
    allocate(svdval_pressure(nsvd,numphase))
    do p=1,numphase
    do d=1,dim
       call snapsvd(u_nodes,snapshots,snapmatrix_velocity(:,:,d,p),&
            nsvd,nsvd,leftsvd_velocity(:,:,d,p),svdval_velocity(:,d,p))
    end do

    call snapsvd(p_nodes,snapshots,snapmatrix_pressure(:,:,p),nsvd,nsvd,leftsvd_pressure(:,:,p),svdval_pressure(:,p))
    enddo !numphase

  end subroutine form_svd

  subroutine form_podstate(state, pod_state, leftsvd_u, leftsvd_p, snapmean_u, snapmean_p)

    type(state_type), intent(in), dimension(:) :: state
    type(state_type), intent(out), dimension(:,:), allocatable :: pod_state

    real, intent(in), dimension(:,:,:) :: leftsvd_u
    real, intent(in), dimension(:,:) :: leftsvd_p
    real, intent(in), dimension(:,:) :: snapmean_u
    real, intent(in), dimension(:) :: snapmean_p 

    type(mesh_type), pointer :: pod_xmesh, pod_umesh, pod_pmesh, pmesh, pod_mesh
    type(element_type) :: pod_xshape, pod_ushape, pod_pshape
    type(vector_field), pointer :: pod_positions, velocity
    type(scalar_field), pointer :: pressure

    type(vector_field) :: pod_velocity
    type(scalar_field) :: pod_pressure
    type(vector_field) :: snapmean_velocity
    type(scalar_field) :: snapmean_pressure

    real, dimension(:), pointer :: x_ptr,y_ptr,z_ptr
    real, dimension(:), allocatable :: x,y,z   

    character(len=1024) :: filename
    character(len = FIELD_NAME_LEN) :: field_name

    integer :: dump_sampling_period, quadrature_degree,nonods
    integer :: i,j,k,nod,total_dumps,POD_num,stat,f,d
    logical :: all_meshes_same
!    type(vector_field) :: podVelocity, newpodVelocity

    call get_option(&
         "/reduced_model/pod_basis_formation/pod_basis_count", POD_num) 

    allocate(pod_state(POD_num,1))
    call nullify(pod_state)

    do i = 1,POD_num

      ! pod_mesh => extract_mesh(state(1), "Mesh")

       all_meshes_same = .true.

        pod_xmesh => extract_mesh(state(1), "CoordinateMesh", stat)
        pod_umesh => extract_mesh(state(1), "VelocityMesh", stat)
       pod_pmesh => extract_mesh(state(1), "PressureMesh", stat)
 
       pod_positions => extract_vector_field(state(1), "Coordinate")

        call insert(pod_state(i,1), pod_xmesh, "Mesh")
        call insert(pod_state(i,1), pod_xmesh, "CoordinateMesh")
        call insert(pod_state(i,1), pod_umesh, "VelocityMesh")
        call insert(pod_state(i,1), pod_pmesh, "PressureMesh")
        call insert(pod_state(i,1), pod_positions, "Coordinate")

       velocity => extract_vector_field(state(1), "Velocity")

       call allocate(pod_velocity, velocity%dim, pod_umesh, "Velocity")
       call zero(pod_velocity)
       do d=1,velocity%dim
          call set_all(pod_velocity, d, leftsvd_u(:,i,d))
       end do
       call insert(pod_state(i,1), pod_velocity, name="Velocity")
       call deallocate(pod_velocity)

       call allocate(pod_pressure, pod_umesh, "PODPressure")
       call zero(pod_pressure)
       call set_all(pod_pressure, leftsvd_p(:,i))
       call insert(pod_state(i,1), pod_pressure, name="Pressure")
       call deallocate(pod_pressure)

!!insert snapmean data into state
     
       call allocate(snapmean_velocity, velocity%dim, pod_umesh, "SnapmeanVelocity")
       call zero(snapmean_velocity)
       do d=1,velocity%dim
          call set_all(snapmean_velocity, d, snapmean_u(:,d))
       end do
       call insert(pod_state(i,1), snapmean_velocity, name="SnapmeanVelocity")
       call deallocate(snapmean_velocity)

       call allocate(snapmean_pressure, pod_umesh, "SnapmeanPressure")
       call zero(snapmean_pressure)
       call set_all(snapmean_pressure, snapmean_p(:))
       call insert(pod_state(i,1), snapmean_pressure, name="SnapmeanPressure")
       call deallocate(snapmean_pressure)


    enddo

  end subroutine form_podstate

  
 subroutine form_podstate_velocity(state, pod_state, leftsvd_u, leftsvd_p, snapmean_u, snapmean_p)
   
    type(state_type), intent(in), dimension(:,:) :: state    
    type(state_type), intent(out), dimension(:,:,:), allocatable :: pod_state
    real, intent(in), dimension(:,:,:,:) :: leftsvd_u
    real, intent(in), dimension(:,:,:) :: leftsvd_p
    real, intent(in), dimension(:,:,:) :: snapmean_u
    real, intent(in), dimension(:,:) :: snapmean_p 

    type(mesh_type), pointer :: pod_xmesh, pod_umesh, pod_pmesh, pmesh, pod_mesh
    type(element_type) :: pod_xshape, pod_ushape, pod_pshape
    type(vector_field), pointer :: pod_positions, velocity
    type(scalar_field), pointer :: pressure

    type(vector_field) :: pod_velocity
    type(scalar_field) :: pod_pressure
    type(vector_field) :: snapmean_velocity
    type(scalar_field) :: snapmean_pressure

    real, dimension(:), pointer :: x_ptr,y_ptr,z_ptr
    real, dimension(:), allocatable :: x,y,z   
    
    character(len=1024) :: filename
    character(len = FIELD_NAME_LEN) :: field_name

    integer :: dump_sampling_period, quadrature_degree,nonods
    integer :: i,j,k,nod,total_dumps,POD_num,stat,f,d,p
    logical :: all_meshes_same
    integer :: pod_pnodes,pod_unodes

    call get_option(&
         "/reduced_model/pod_basis_formation/pod_basis_count", POD_num) 

    allocate(pod_state(POD_num,1,numphase))
    call nullify(pod_state)
    
    do p=1,numphase
    do i = 1,POD_num

      ! pod_mesh => extract_mesh(state(1), "Mesh")

       all_meshes_same = .true.

        print *, 'aaaaaaaaaaaaaaaaaaaaaaaa'
        pod_xmesh => extract_mesh(state(1,p), "Mesh")
        print *, 'bbbbbbbbbbbbbbbbbbbbbbbb'
        pod_umesh =>extract_mesh(state(1,p), "Mesh")
        pod_positions => extract_vector_field(state(1,p), "Coordinate")

 
       call insert(pod_state(i,1,p), pod_umesh, "CoordinateMesh")
       call insert(pod_state(i,1,p), pod_umesh, "VelocityMesh")
    
       call insert(pod_state(i,1,p), pod_positions, "Coordinate")

       velocity => extract_vector_field(state(1,p), "Velocity")

       call allocate(pod_velocity, velocity%dim, pod_umesh, "PODVelocity")
       call zero(pod_velocity)
       do d=1,velocity%dim
          call set_all(pod_velocity, d, leftsvd_u(:,i,d,p))
       end do
       call insert(pod_state(i,1,p), pod_velocity, name="PODVelocity")
       call deallocate(pod_velocity) 
       !!insert snapmean data into state
     
       call allocate(snapmean_velocity, velocity%dim, pod_umesh, "SnapmeanVelocity")
       call zero(snapmean_velocity)
       do d=1,velocity%dim
          call set_all(snapmean_velocity, d, snapmean_u(:,d,p))
       end do
       call insert(pod_state(i,1,p), snapmean_velocity, name="SnapmeanVelocity")
       call deallocate(snapmean_velocity)

   
  !test the podnodes        

      velocity => extract_vector_field(pod_state(i,1,p), "PODVelocity")
   
      pod_unodes=node_count(velocity)
     
      print*, "podddddddddddddddddddddddvelocity", pod_unodes
   
    enddo  
     enddo !numphase
  end subroutine form_podstate_velocity

  subroutine form_podstate_pressure(state_p, pod_state_p, leftsvd_u, leftsvd_p, snapmean_u, snapmean_p)
   
    type(state_type), intent(in), dimension(:,:) :: state_p
   ! type(state_type), pointer, dimension(:):: state
    type(state_type), intent(out), dimension(:,:,:), allocatable :: pod_state_p 
    real, intent(in), dimension(:,:,:,:) :: leftsvd_u
    real, intent(in), dimension(:,:,:) :: leftsvd_p
    real, intent(in), dimension(:,:,:) :: snapmean_u
    real, intent(in), dimension(:,:) :: snapmean_p 

    type(mesh_type), pointer :: pod_xmesh, pod_umesh, pod_pmesh, pmesh, pod_mesh
    type(element_type) :: pod_xshape, pod_ushape, pod_pshape
    type(vector_field), pointer :: pod_positions, velocity
    type(scalar_field), pointer :: pressure

    type(vector_field) :: pod_velocity
    type(scalar_field) :: pod_pressure
    type(vector_field) :: snapmean_velocity
    type(scalar_field) :: snapmean_pressure

    real, dimension(:), pointer :: x_ptr,y_ptr,z_ptr
    real, dimension(:), allocatable :: x,y,z   

    character(len=1024) :: filename
    character(len = FIELD_NAME_LEN) :: field_name

    integer :: dump_sampling_period, quadrature_degree,nonods
    integer :: i,j,k,nod,total_dumps,POD_num,stat,f,d,p
    logical :: all_meshes_same
    integer :: pod_pnodes,pod_unodes  

    call get_option(&
         "/reduced_model/pod_basis_formation/pod_basis_count", POD_num) 

    allocate(pod_state_p(POD_num,1,numphase))  
    call nullify(pod_state_p)
    do p=1,numphase
    do i = 1,POD_num

      ! pod_mesh => extract_mesh(state(1), "Mesh")

       all_meshes_same = .true.

     
        pod_xmesh => extract_mesh(state_p(1,p), "Mesh")
       
        pod_pmesh =>extract_mesh(state_p(1,p), "Mesh") 
 
       pod_positions => extract_vector_field(state_p(1,p), "Coordinate")

 
       call insert(pod_state_p(i,1,p), pod_xmesh, "CoordinateMesh") 
       call insert(pod_state_p(i,1,p), pod_pmesh, "PressureMesh")
       call insert(pod_state_p(i,1,p), pod_positions, "Coordinate")
     

       call allocate(pod_pressure, pod_pmesh, "PODPressure")    
       call zero(pod_pressure)
       call set_all(pod_pressure, leftsvd_p(:,i,p))
       call insert(pod_state_p(i,1,p), pod_pressure, name="PODPressure")
     
       call deallocate(pod_pressure)

       !!insert snapmean data into state 
       call allocate(snapmean_pressure, pod_pmesh, "SnapmeanPressure")
       call zero(snapmean_pressure)
       call set_all(snapmean_pressure, snapmean_p(:,p))
       call insert(pod_state_p(i,1,p), snapmean_pressure, name="SnapmeanPressure") 
       call deallocate(snapmean_pressure)
  !test the podnodes        

       
      pressure => extract_scalar_field(pod_state_p(i,1,p), "PODPressure")       
      pod_pnodes=node_count(pressure)  
   
      print*, "podddddddddddddddddddddddpressure", pod_pnodes
   
    enddo
    enddo     ! numphase  
  end subroutine form_podstate_pressure


 subroutine read_input_states_fromvtu(state,state_p)
    !!< Read the input states from the vtu dumps of the forward run.
    type(state_type), intent(out), dimension(:,:), allocatable :: state
    type(state_type), intent(out), dimension(:,:), allocatable :: state_p
    character(len=1024) :: filename
     type(scalar_field) ::pressure
    integer :: dump_sampling_period, quadrature_degree
    integer :: i,j,k,total_dumps,stable_dumps
    type(mesh_type) ::pressuremesh
    integer :: numphase
    character(len=1024) :: phase_name
    
    call get_option('/reduced_model/pod_basis_formation/dump_sampling_period',dump_sampling_period)
    call get_option('/geometry/quadrature/degree', quadrature_degree)

    ewrite(3,*) 'dump_sampling_period',dump_sampling_period
    
    !substract gyre_0.vtu
    total_dumps=count_dumps(dump_sampling_period)-1
    numphase=option_count("/material_phase")
    allocate(state(total_dumps,numphase))
    allocate(state_p(total_dumps,numphase))

!    stable_dumps=total_dumps-10
!    allocate(state(stable_dumps))
!vtu--->state
    do j=1, numphase
         call get_option("/material_phase["//int2str(j-1)//"]/name", phase_name) !not ' but "
         
        ! call get_option("/material_phase["//int2str(m)//"]/name", mat_name)
         print *,'phase_name=',  int2str(j), phase_name
    do i=1, total_dumps
      
          write(filename, '(a, i0, a,a)')  trim(simulation_name)//'_'//trim(phase_name)//'_VelocityMesh_', (i)*dump_sampling_period,'_checkpoint',".vtu" 
          call vtk_read_state(filename, state(i,j), quadrature_degree)
          write(filename, '(a, i0, a,a)')  trim(simulation_name)//'_'//trim(phase_name)//'_PressureMesh_', (i)*dump_sampling_period,'_checkpoint',".vtu" 
          call vtk_read_state(filename, state_p(i,j), quadrature_degree)
       
      ! pressure=extract_scalar_field(state_p(i),"Pressure")
      !  call insert(state(i), pressure, name="Pressure") 
      ! pressuremesh =extract_mesh(state_p(i), "Mesh") 
      ! call insert(state(i), pressuremesh, name="PressureMesh") 
    end do
    enddo  

       ! call print_state(state(1))
      !  call print_state(state_p(1))
      !   call print_state(state(2))
       ! call print_state(state_p(2))
  end subroutine read_input_states_fromvtu

 

  function count_dumps(dump_sampling_period) result (count)
    !! Work out how many dumps we're going to read in.
    integer :: count,dump_sampling_period

    logical :: exists
    !      character(len=FILE_NAME_LEN) :: filename
    character(len=1024) :: filename, phase_name

    count=1
    call get_option("/material_phase["//int2str(0)//"]/name", phase_name)
    do 
       !! Note that this won't work in parallel. Have to look for the pvtu in that case.
       !! not ' but "
       write(filename, '(a, i0, a,a)')  trim(simulation_name)//'_'//trim(phase_name)//'_VelocityMesh_', (count)*dump_sampling_period,'_checkpoint',".vtu" 
       ! write(filename, '(a, i0, a)') trim(simulation_name)//'_', (count-1)*dump_sampling_period,".vtu" 
       inquire(file=trim(filename), exist=exists)
       if (.not. exists) then
          if(dump_sampling_period.eq.1) then
             count=count-2
          else
             count=count-1
          endif
          exit
       end if

       count=count+1
    end do

    if (count==0) then
       FLExit("No .vtu files found!")
    end if

  end function count_dumps


  subroutine read_command_line()
    implicit none
    ! Read the input filename.
    character(len=1024) :: argument, filename
    integer :: status, argn, level

    call set_global_debug_level(0)

    argn=1
    do 

       call get_command_argument(argn, value=argument, status=status)
       argn=argn+1

       if (status/=0) then
          call usage
          stop
       end if

       if (argument=="-v") then
          call get_command_argument(argn, value=argument, status=status)
          argn=argn+1

          if (status/=0) then
             call usage
             stop
          end if

          read(argument, "(i1)", err=666) level
          call set_global_debug_level(level)

       else

          ! Must be the filename
          filename=argument

       end if

       if (argn>=command_argument_count()) exit
    end do

    call load_options(filename)

    return

666 call usage
    stop

  end subroutine read_command_line

  subroutine usage
    implicit none

    write (0,*) "usage: form_pod_basis [-v n] <options_file>"
    write (0,*) ""
    write (0,*) "-v n sets the verbosity of debugging"
  end subroutine usage

   

  subroutine solution_on_checkpoints(state, prefix, postfix, cp_no)
    !!< Checkpoint the fields in state. Outputs to vtu files with names:
    !!<   [prefix]_[_state name]_[mesh_name][_cp_no][_postfix][_process].vtu
    !!< where the state name is added if multiple states are passed, cp_no is
    !!< optional and the process number is added in parallel.

    type(state_type), dimension(:), intent(in) :: state
    character(len = *), intent(in) :: prefix
    character(len = *), optional, intent(in) :: postfix
    integer, optional, intent(in) :: cp_no
    ! if present and true: do not checkpoint fields that can be reinitialised

    character(len = OPTION_PATH_LEN) :: vtu_filename,tmp
    integer :: i, j, k, nparts, n_ps_fields_on_mesh, n_pv_fields_on_mesh, n_pt_fields_on_mesh
    type(mesh_type), pointer :: mesh
    type(scalar_field), pointer :: s_field
    type(vector_field), pointer :: positions, v_field
    type(tensor_field), pointer :: t_field

    integer :: stat

    assert(len_trim(prefix) > 0)

    do i = 1, size(state)
       positions => extract_vector_field(state(i), "Coordinate")
       do j = 1, size(state(i)%meshes)
          mesh => state(i)%meshes(j)%ptr
          nparts = get_active_nparts(ele_count(mesh))
          ! Construct a new field checkpoint filename
          vtu_filename = trim(prefix)
          if(size(state) > 1) vtu_filename = trim(vtu_filename) // "_" // trim(state(i)%name)
          vtu_filename = trim(vtu_filename) // "_" // trim(mesh%name)
         if(present(cp_no)) vtu_filename = trim(vtu_filename) // "_" // int2str(cp_no)
          if(present_and_nonempty(postfix)) vtu_filename = trim(vtu_filename) // "_" // trim(postfix)
          if(nparts > 1) then
             vtu_filename = trim(vtu_filename) // ".pvtu"
          else
             vtu_filename = trim(vtu_filename) // ".vtu"
          end if
          

          !        if(associated(state(i)%scalar_fields)) then
          do k = 1, size(state(i)%scalar_fields)
             s_field => state(i)%scalar_fields(k)%ptr              
             if(trim(s_field%mesh%name) == trim(mesh%name) .and. s_field%mesh == mesh) then
              if(have_option(trim(s_field%option_path) // "/prognostic")) then
                call update_initial_condition_options(trim(s_field%option_path), trim(vtu_filename), "vtu")
              end if
             endif
          end do

          do k = 1, size(state(i)%vector_fields)
             v_field =>  state(i)%vector_fields(k)%ptr
!             if( trim(v_field%name)=='Coordinate' .or. trim(v_field%name)=='PeriodicMeshCoordinate' ) cycle
             if(trim(v_field%mesh%name) == trim(mesh%name) .and. v_field%mesh == mesh) then
            
              if(have_option(trim(v_field%option_path) // "/prognostic")) then
                call update_initial_condition_options(trim(v_field%option_path), trim(vtu_filename), "vtu")
              end if

             endif
          end do
          !        end if
          !        if(associated(state(i)%tensor_fields)) then
          do k = 1, size(state(i)%tensor_fields)
             t_field => state(i)%tensor_fields(k)%ptr
             
             if(trim(t_field%mesh%name) == trim(mesh%name) .and. t_field%mesh == mesh) then
              if(have_option(trim(t_field%option_path) // "/prognostic")) then
                call update_initial_condition_options(trim(t_field%option_path), trim(vtu_filename), "vtu")
              end if
             endif
          end do
          !        end if
       end do
              if(have_option(trim(s_field%option_path) // "/prognostic")) print*,'%%%%%%%%%%%%%%%%%%%%%'
    end do
  end subroutine solution_on_checkpoints
  
    subroutine update_initial_condition_options(path, filename, format)
      !!< Updates the initial condition options for a prognostic field with
      !!< options path path

      character(len = *), intent(in) :: path
      character(len = *), intent(in) :: filename
      character(len = *), intent(in) :: format

      integer :: stat, ic, nics

      nics = option_count(trim(path) // "/prognostic/initial_condition")
      do ic = 0, nics-1  ! do while seemed to break, don't know why
        call delete_option(trim(path) // "/prognostic/initial_condition[" // int2str(0) // "]")
      end do
      call set_option_attribute(trim(path) // "/prognostic/initial_condition::WholeMesh/from_file/file_name", trim(filename), stat)
      if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
        FLAbort("Failed to set initial condition filename when checkpointing field with option path " // trim(path))
      end if
      call set_option_attribute(trim(path) // "/prognostic/initial_condition::WholeMesh/from_file/format/name", trim(format), stat)
      if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
        FLAbort("Failed to set initial condition format when checkpointing field with option path " // trim(path))
      end if

    end subroutine update_initial_condition_options

    subroutine update_value_options(path, filename, format)
      !!< Updates the value options for a prescribed field with
      !!< options path path

      character(len = *), intent(in) :: path
      character(len = *), intent(in) :: filename
      character(len = *), intent(in) :: format

      integer :: stat, value, nvalues

      nvalues = option_count(trim(path) // "/prescribed/value")
      do value = 0, nvalues-1
        call delete_option(trim(path) // "/prescribed/value[" // int2str(0) // "]")
      end do
      call set_option_attribute(trim(path) // "/prescribed/value::WholeMesh/from_file/file_name", trim(filename), stat)
      if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
        FLAbort("Failed to set value filename when checkpointing field with option path " // trim(path))
      end if
      call set_option_attribute(trim(path) // "/prescribed/value::WholeMesh/from_file/format/name", trim(format), stat)
      if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
        FLAbort("Failed to set value format when checkpointing field with option path " // trim(path))
      end if

    end subroutine update_value_options

   subroutine update_initial_condition_options2(path, filename, format)
      !!< Updates the initial condition options for a prognostic field with
      !!< options path path

      character(len = *), intent(in) :: path
      character(len = *), intent(in) :: filename
      character(len = *), intent(in) :: format

      integer :: stat, ic, nics

        call delete_option(trim(path) // "/prognostic/initial_condition[" // int2str(0) // "]")
      call set_option_attribute(trim(path) // "/prognostic/initial_condition::WholeMesh/from_file/file_name", trim(filename), stat)
      if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
        FLAbort("Failed to set initial condition filename when checkpointing field with option path " // trim(path))
      end if
      call set_option_attribute(trim(path) // "/prognostic/initial_condition::WholeMesh/from_file/format/name", trim(format), stat)
      if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
        FLAbort("Failed to set initial condition format when checkpointing field with option path " // trim(path))
      end if

    end subroutine update_initial_condition_options2


function DEIMP(Uin,n,m) result(phi)
 !This routine will produce a set of m EIM interpolation indices p for the input basis Uin

 real, dimension(:,:), allocatable :: UC,P,PT,E,A!uIn !UC : changing dimension m(n*l)
 !P:n*m(unodes_nsvd(Basis number of non-linear term )   Uin:n*m   m number of interpolation ,E: identity matrix 
 
 real,dimension(:), allocatable :: uL,b,c,r   !phi store the location of maximum value
 integer, dimension(:), allocatable :: phi
 integer :: n,m,i,l,j
 integer :: locu(1)
 real :: uIn(n,m)
 print * ,'sizeofuin', size(uIn,1),size(uIn,2),deim_number
 !allocate(uIn(n,m))
 allocate(PT(m,n))
 allocate(UC(n,m))  
 allocate(uL(n))
 allocate(P(n,m)) 
 allocate(A(m,m))
 allocate(phi(m))
 allocate(b(m))
 allocate(c(m))
 allocate(r(n))
 allocate(E(n,n))
 

 phi=0
 locu= MAXLOC(abs(uIn(:,1)))    ! Rank of U should be paid attention to
 phi(1)=locu(1)


FORALL (i=1:n,j=1:n) 
 E(i,j) = 0
 END FORALL

 FORALL (i=1:n) 
 E(i,i) = 1
 END FORALL

 P(:, 1) = E(:, phi(1))  
 UC(:, 1) = uIn(:, 1)
 
 do l=2,m
  uL(:) = Uin(:,l)    
  
 do i=1,n  ! transpose
   do j=1,l-1
     PT(j,i)=P(i,j)
   enddo
  enddo
  j=l-1
  A(1:j,1:j)=matmul(PT(1:j, :),UC(:, 1:j))
  b(1:j)=matmul(PT(1:j, :),uL(:))
  
 call solve(A(1:j,1:j), b(1:j)) 
 !call LEGS_POD(A,m,c,b) !A(m,m)*b(m) = c(m)
 !!!! c = (P'*U)\P'*uL;  P' transpose matrix of P  . After solving, put the results to b.
  
 c(1:j)=b(1:j)  ! give the results of solver to c
 
 r = uL - matmul(UC(:,1:j),c(1:j))  !!
 locu= MAXLOC(abs(r))
 print *, 'absresiduleabsresiduleabsresidule', maxval(abs(r))
 phi(l)=locu(1)
 UC(:, l) = uL(:)
 P(:, l) = E(:, phi(l))
 enddo 
 !P = E(:, phi)
  call selection_sort(phi,m)
  
  ! print *, 'phiphiphiphiphiphiphiphiphiphi', phi
  !open(100,file='indices')
  ! write(100,*)(phi(:))
  !close(100)

 end function DEIMP

subroutine selection_sort(a,n)
  implicit none
  integer :: n,a(n)
  integer i,j  ! loop counter
  integer min  ! find a minimum value of one loop
  integer temp ! 
!sort
  do i=1,n-1
    min=a(i)     !  
    do j=i+1,n
      if ( min > a(j) ) then      
        temp=a(j)          
        a(j)=a(i)
        a(i)=temp
        min=a(i)
      end if
 end do
  end do                             
  return
end subroutine 

end program form_pod_basis
