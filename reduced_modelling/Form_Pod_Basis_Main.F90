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
  type(state_type), dimension(:), allocatable :: state,state_test,state_deim
  type(state_type), dimension(:,:), allocatable :: pod_state, pod_deim_state
  type(state_type), dimension(:,:), allocatable :: pod_state_p
  type(vector_field) , pointer :: position_deim
  integer :: timestep
  integer :: ierr

  character(len = OPTION_PATH_LEN) :: simulation_name

#ifdef HAVE_MPI
  call mpi_init(ierr)
#endif

#ifdef HAVE_PETSC
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

  call python_init()
  call read_command_line()

  call form_basis_different_mesh()
  call form_basis()

  call deallocate(state)
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
    real, dimension(:,:,:), allocatable :: snapmatrix_velocity, snapmatrix_velocity_deim
    real, dimension(:,:), allocatable :: snapmatrix_pressure, snapmatrix_pressure_deim
    real, dimension(:,:,:), allocatable :: leftsvd_velocity, leftsvd_velocity_deim
    real, dimension(:,:), allocatable :: leftsvd_pressure, leftsvd_pressure_deim
    real, dimension(:,:), allocatable :: svdval_velocity, svdval_velocity_deim
    real, dimension(:), allocatable :: svdval_pressure, svdval_pressure_deim
    real, dimension(:,:), allocatable :: snapmean_velocity, snapmean_velocity_deim
    real, dimension(:), allocatable :: snapmean_pressure, snapmean_pressure_deim
    integer :: snapshots, u_nodes, p_nodes, nsvd
    integer :: i,dump_no, d, dim,j,k ,i2,i3
    integer :: stat

    real, dimension(:,:), allocatable :: P !n*m (u_nodes*nsvd)
    !integer, dimension(:), allocatable :: phi
    real, dimension(:,:), allocatable :: leftsvd_velocity_deim_uin
    integer :: udim,deim_number_tmp  
    type(vector_field), pointer :: velocityudim

    type(vector_field), pointer :: vfield

    call get_option(&
         '/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
  !  call get_option('/simulation_name',simulation_name)
  !  call read_input_states(state)
  !  call retrieve_snapshots(state, snapshots, u_nodes, p_nodes, snapmatrix_velocity, snapmatrix_pressure, &
                  !          & snapmean_velocity, snapmean_pressure)

    
   ! call form_svd(snapmatrix_velocity, snapmatrix_pressure,&
    !   & leftsvd_velocity, leftsvd_pressure, svdval_velocity, svdval_pressure, snapshots)
    !call form_podstate(state,pod_state,leftsvd_velocity,leftsvd_pressure, snapmean_velocity, snapmean_pressure)

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
      velocityudim => extract_vector_field(state(1), "Velocity")   
      udim=velocityudim%dim
      allocate(P(u_nodes,deim_number)) 
      deim_number_tmp=deim_number!*40    
      allocate(indices(udim*deim_number))  
      allocate(indices_tmp(deim_number))
      print *,'sizeofindice_tmp', size(indices_tmp,1)
     ! allocate(leftsvd_velocity_deim_uin(u_nodes*udim,nsvd))     
      
       
      !allocate(phi(nsvd))
      call read_input_states_deim(state_deim) !deim
      call retrieve_snapshots(state_deim, snapshots, u_nodes, p_nodes, snapmatrix_velocity_deim, snapmatrix_pressure_deim, &
                            & snapmean_velocity_deim, snapmean_pressure_deim)                                                    !deim
      call form_svd(snapmatrix_velocity_deim, snapmatrix_pressure_deim,&  
                           & leftsvd_velocity_deim, leftsvd_pressure_deim, svdval_velocity_deim, svdval_pressure_deim, &
                           &snapshots)                                !deim   leftsvd_velocity_deim  U
      ! call form_podstate(state_deim,pod_deim_state,leftsvd_velocity_deim,leftsvd_pressure_deim, snapmean_velocity_deim, snapmean_pressure_deim) 
       call form_podstate_writeout_allvariables_differentmesh(state_deim,pod_deim_state,leftsvd_velocity_deim,leftsvd_pressure_deim, snapmean_velocity_deim, snapmean_pressure_deim) 
       print *, 'u_nodes', u_nodes   
     ! do i=1,nsvd    ! deim_number  !need to change?
      !   dump_no=i
       !  call vtk_write_state(filename=trim(simulation_name)//"_PODBasisDEIM", index=dump_no, state=pod_deim_state(i,:))
       !  call deallocate(pod_deim_state(i,:)) 
     ! enddo
     !  deallocate(phi)   
    !  print *, 'sizeofleftsvd_velocity_deim_uin',size(leftsvd_velocity_deim_uin,1), size(leftsvd_velocity_deim_uin,2)
     position_deim=>extract_vector_field(state(1),"Coordinate")
     
    do i=2,udim
    ! print *, 'sizeofleftsvd_velocity_deim',size(leftsvd_velocity_deim,1), size(leftsvd_velocity_deim,2),size(leftsvd_velocity_deim,3) 
                             !3213      24           2
      indices_tmp = DEIMP(leftsvd_velocity_deim(:,:,i),u_nodes,deim_number)  !(u_nodes,nsvd,dim))  !redefine phi, and leftsvd_velocity()to 1 dimension

      indices(1:deim_number)=indices_tmp
      indices((i-1)*(deim_number)+1:(deim_number*i))=indices_tmp
      do j=1, deim_number
        print *, position_deim%val(1,indices_tmp(j)),',',position_deim%val(2,indices_tmp(j))        
      enddo
      
      enddo
     !P = DEIMP(leftsvd_velocity_deim_uin(:,:),u_nodes*udim,nsvd)   
     print *, 'indices indices indices indices indices', indices
     print *,'whole domain', position_deim%val(1,deim_number), position_deim%val(2,deim_number)
     
     
      open(100,file='indices')
      write(100,*)(indices(:))
      close(100)
     end if
  end subroutine form_basis

  subroutine form_basis_different_mesh()
    !!< Matrices containing the snapshots for arpack
    !      type(state_type), dimension(:), allocatable :: state

    real, dimension(:,:,:), allocatable :: snapmatrix_velocity
    real, dimension(:,:), allocatable :: snapmatrix_pressure
    real, dimension(:,:,:), allocatable :: leftsvd_velocity
    real, dimension(:,:), allocatable :: leftsvd_pressure
    real, dimension(:,:), allocatable :: svdval_velocity
    real, dimension(:), allocatable :: svdval_pressure
    real, dimension(:,:), allocatable :: snapmean_velocity
    real, dimension(:), allocatable :: snapmean_pressure
    integer :: snapshots, u_nodes, p_nodes, nsvd
    integer :: i,dump_no, d, dim
    integer :: stat
    type(vector_field), pointer :: vfield
    character(len = OPTION_PATH_LEN) :: prefix,postfix

    call get_option(&
         '/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
    call get_option('/simulation_name',simulation_name)
    call read_input_states(state)
    call retrieve_snapshots(state, snapshots, u_nodes, p_nodes, snapmatrix_velocity, snapmatrix_pressure, &
                            & snapmean_velocity, snapmean_pressure)

    
    call form_svd(snapmatrix_velocity, snapmatrix_pressure,&
       & leftsvd_velocity, leftsvd_pressure, svdval_velocity, svdval_pressure, snapshots)

!    call form_podstate(state,pod_state,leftsvd_velocity,leftsvd_pressure, snapmean_velocity, snapmean_pressure)
    call form_podstate_writeout_allvariables_differentmesh(state, pod_state, leftsvd_velocity, leftsvd_pressure, &
         snapmean_velocity, snapmean_pressure)

    !! Produce updated flml file in which the execute_reduced_model
    !! option is set.
    call add_option("/reduced_model/execute_reduced_model",stat)
    call set_option('/simulation_name', trim(simulation_name)//'_POD') 
     !! correct the initial conditions
    prefix = simulation_name
    postfix = "checkpoint"
    call solution_on_checkpoints((/state(1)/), prefix, postfix, 0)

    call write_options(trim(simulation_name)//"_POD.flml")

  end subroutine form_basis_different_mesh

  subroutine retrieve_snapshots(state, snapshots, u_nodes, p_nodes, snapmatrix_velocity, snapmatrix_pressure, &
                                & snapmean_velocity, snapmean_pressure)
    !!< 
    type(state_type), intent(in), dimension(:) :: state
    real, dimension(:,:,:), allocatable :: snapmatrix_velocity
    real, dimension(:,:), allocatable :: snapmatrix_pressure
    real, dimension(:,:), allocatable :: snapmean_velocity
    real, dimension(:), allocatable :: snapmean_pressure

    integer :: snapshots, u_nodes, p_nodes, dim, i, d
    type(vector_field), pointer :: velocity
    type(scalar_field), pointer :: pressure

    velocity => extract_vector_field(state(1), "Velocity")
    pressure => extract_scalar_field(state(1), "Pressure")

    dim=velocity%dim
    p_nodes=node_count(pressure)
    u_nodes=node_count(velocity)

    snapshots=size(state)
!print*,snapshots,dim

    allocate(snapmatrix_velocity(u_nodes,snapshots,dim))
    allocate(snapmatrix_pressure(p_nodes,snapshots))
    allocate(snapmean_velocity(u_nodes,dim))
    allocate(snapmean_pressure(p_nodes))
    
    snapmatrix_velocity=0.0
    snapmatrix_pressure=0.0
    snapmean_velocity=0.0
    snapmean_pressure=0.0
    do i = 1, snapshots
       velocity => extract_vector_field(state(i), "Velocity")
       pressure => extract_scalar_field(state(i), "Pressure")

       do d = 1, dim
!          snapmatrix_velocity(:,i,d)=field_val(velocity,d)
           snapmatrix_velocity(:,i,d)=velocity%val(d,:)
       end do
!       snapmatrix_pressure(:,i)=field_val(pressure)
        snapmatrix_pressure(:,i)=pressure%val
    end do


    do i=1, snapshots

       do d=1, dim
          snapmean_velocity(:,d)= snapmean_velocity(:,d)+snapmatrix_velocity(:,i,d)
       enddo
       snapmean_pressure(:)=snapmean_pressure(:)+snapmatrix_pressure(:,i)
    end do

    do d=1,dim
       snapmean_velocity(:,d)=snapmean_velocity(:,d)/snapshots
    enddo
       snapmean_pressure(:)=snapmean_pressure(:)/snapshots

    do i=1,snapshots
       do d=1,dim
          snapmatrix_velocity(:,i,d)=snapmatrix_velocity(:,i,d)-snapmean_velocity(:,d)
       enddo
       snapmatrix_pressure(:,i)=snapmatrix_pressure(:,i)-snapmean_pressure(:)
    enddo

  end subroutine retrieve_snapshots

  subroutine form_svd(snapmatrix_velocity, snapmatrix_pressure,&
       & leftsvd_velocity, leftsvd_pressure, svdval_velocity, svdval_pressure, snapshots)
    
    real, dimension(:,:,:), intent(in) :: snapmatrix_velocity
    real, dimension(:,:), intent(in) :: snapmatrix_pressure
    real, dimension(:,:,:), allocatable, intent(out) :: leftsvd_velocity
    real, dimension(:,:), allocatable, intent(out) :: leftsvd_pressure
    real, dimension(:,:), allocatable, intent(out) :: svdval_velocity
    real, dimension(:), allocatable, intent(out) :: svdval_pressure
    integer i, d, dim ,nsvd, snapshots, p_nodes, u_nodes

    call get_option(&
         '/reduced_model/pod_basis_formation/pod_basis_count', nsvd)

    dim=size(snapmatrix_velocity,3)
    p_nodes=size(snapmatrix_pressure,1)
    u_nodes=size(snapmatrix_velocity,1)

    allocate(leftsvd_velocity(u_nodes,nsvd,dim))
    allocate(leftsvd_pressure(p_nodes,nsvd))
    allocate(svdval_velocity(nsvd,dim))
    allocate(svdval_pressure(nsvd))

    do d=1,dim
       call snapsvd(u_nodes,snapshots,snapmatrix_velocity(:,:,d),&
            nsvd,nsvd,leftsvd_velocity(:,:,d),svdval_velocity(:,d))
    end do

    call snapsvd(p_nodes,snapshots,snapmatrix_pressure,nsvd,nsvd,leftsvd_pressure,svdval_pressure)

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

  subroutine form_podstate_writeout_allvariables_differentmesh(state, pod_state, leftsvd_u, leftsvd_p, snapmean_u, snapmean_p)

    type(state_type), intent(in), dimension(:) :: state
    type(state_type), intent(out), dimension(:,:), allocatable :: pod_state

    real, intent(in), dimension(:,:,:) :: leftsvd_u
    real, intent(in), dimension(:,:) :: leftsvd_p
    real, intent(in), dimension(:,:) :: snapmean_u
    real, intent(in), dimension(:) :: snapmean_p 

    type(mesh_type), pointer :: pod_xmesh, pod_umesh, pod_pmesh, pmesh, pod_mesh
    type(element_type) :: pod_xshape, pod_ushape, pod_pshape
    type(vector_field), pointer :: pod_positions, v_field,velocity
    type(scalar_field), pointer :: s_field

    type(vector_field) :: pod_velocity
    type(scalar_field) :: pod_pressure
    type(vector_field) :: snapmean_velocity
    type(scalar_field) :: snapmean_pressure

    real, dimension(:), pointer :: x_ptr,y_ptr,z_ptr
    real, dimension(:), allocatable :: x,y,z   

    character(len=1024) :: filename
    character(len = FIELD_NAME_LEN) :: field_name, mesh_path, mesh_name,mesh_name2

    integer :: dump_sampling_period, quadrature_degree,nonods
    integer :: i,j,k,nod,total_dumps,POD_num,stat,f,d,nfield,m,n_meshes
    logical :: all_meshes_same
!    type(vector_field) :: podVelocity, newpodVelocity

    call get_option(&
         "/reduced_model/pod_basis_formation/pod_basis_count", POD_num) 

    allocate(pod_state(POD_num,1)) ! for multiphase will fix it later
    do k =1, 1

       print*,'%%%%%%%kkkkk=',k
       ! Vector field
       !-------------
       nfield = vector_field_count( state(1) )
       do j = 1, nfield
          v_field => state(k)%vector_fields(j)%ptr              
          if(have_option(trim(v_field%option_path) // "/prognostic")) then
          
             call nullify(pod_state(:,k))
             all_meshes_same = .true.
             ! pod_mesh => extract_mesh(state(1), trim(mesh_name))
             pod_xmesh => extract_mesh(state(k), "CoordinateMesh",stat)
             n_meshes = option_count("/geometry/mesh")
             do m = 0, n_meshes - 1
                mesh_path = "/geometry/mesh[" // int2str(m) // "]"
                call get_option(trim(mesh_path) // "/name", mesh_name)
                if(mesh_name==trim(v_field%mesh%name)) then
                   call get_option(trim(mesh_path)// "/from_mesh/mesh[0]/name", mesh_name2)
                   cycle
                endif
             enddo
             print*,mesh_name2,trim(v_field%mesh%name)
             if (trim(mesh_name2)=="CoordinateMesh") then
                pod_positions => extract_vector_field(state(k), "Coordinate")
             else
                pod_positions => extract_vector_field(state(k), trim(mesh_name2)//"Coordinate")
             end if
             
             pod_umesh => extract_velocity_mesh(state(k))
             print*,'pod_xmesh, pod_umesh', node_count(pod_xmesh),node_count(pod_umesh)
             
             do i = 1,POD_num
                call insert(pod_state(i,k), pod_xmesh, "CoordinateMesh")
                call insert(pod_state(i,k), pod_umesh, trim(v_field%mesh%name))
                call insert(pod_state(i,k), pod_positions, "Coordinate")
                
                velocity => extract_vector_field(state(1), "Velocity")
                call allocate(pod_velocity, velocity%dim, pod_umesh, "Velocity")
                call zero(pod_velocity)
                do d=1,velocity%dim
                   call set_all(pod_velocity, d, leftsvd_u(:,i,d))
                end do
                call insert(pod_state(i,k), pod_velocity, name="Velocity")
                call deallocate(pod_velocity)
                
                !!insert snapmean data into state
                call allocate(snapmean_velocity, velocity%dim, pod_umesh, "SnapmeanVelocity")
                call zero(snapmean_velocity)
                do d=1,velocity%dim
                   call set_all(snapmean_velocity, d, snapmean_u(:,d))
                end do
                call insert(pod_state(i,k), snapmean_velocity, name="SnapmeanVelocity")
                call deallocate(snapmean_velocity)
                deim=have_option("/reduced_model/discrete_empirical_interpolation_method")
		if(deim)then
		 call vtk_write_state(filename=trim(simulation_name)//"_PODDEIMBasisRES"//trim(v_field%name),  &
                     index=i, model=trim(v_field%mesh%name), state=(/pod_state(i,k)/))
		else
                 call vtk_write_state(filename=trim(simulation_name)//"_PODBasis"//trim(v_field%name),  &
                     index=i, model=trim(v_field%mesh%name), state=(/pod_state(i,k)/))
		endif 
             enddo
          endif
       end do !j = 1, size(state(i)%scalar_fields)
       ! scaler field
       !-------------

       nfield = scalar_field_count( state(1) )
       do j = 1, nfield
          s_field => state(k)%scalar_fields(j)%ptr              
          if(have_option(trim(s_field%option_path) // "/prognostic")) then
             call nullify(pod_state(:,k))
             all_meshes_same = .true.
             pod_xmesh => extract_mesh(state(k), "CoordinateMesh",stat)
             n_meshes = option_count("/geometry/mesh")
             do m = 0, n_meshes - 1
                mesh_path = "/geometry/mesh[" // int2str(m) // "]"
                call get_option(trim(mesh_path) // "/name", mesh_name)
                if(mesh_name==trim(s_field%mesh%name)) then
                   call get_option(trim(mesh_path)// "/from_mesh/mesh[0]/name", mesh_name2)
                   cycle  !! exit??
                endif
             enddo
             
             if (trim(mesh_name2)=="CoordinateMesh") then
                pod_positions => extract_vector_field(state(k), "Coordinate")
             else
                pod_positions => extract_vector_field(state(k), trim(mesh_name2)//"Coordinate") !!??
             end if
             
             pod_pmesh => extract_mesh(state(k),trim(s_field%mesh%name),stat)
             print*,'pod_xmesh, pod_pmesh', node_count(pod_xmesh),node_count(pod_pmesh)
!             print*,trim(s_field%mesh%name),trim(s_field%name),"Snapmean"//trim(s_field%name)
             do i = 1,POD_num
                call insert(pod_state(i,k), pod_xmesh, "CoordinateMesh")
                call insert(pod_state(i,k), pod_pmesh, trim(s_field%mesh%name))
                call insert(pod_state(i,k), pod_positions, "Coordinate")
                
                call allocate(pod_pressure, pod_pmesh, trim(s_field%name))    
                call zero(pod_pressure)
                call set_all(pod_pressure, leftsvd_p(:,i))
                call insert(pod_state(i,k), pod_pressure, name=trim(s_field%name))
                
                call deallocate(pod_pressure)
                
                !!insert snapmean data into state
                
                call allocate(snapmean_pressure, pod_pmesh, "Snapmean"//trim(s_field%name))
                call zero(snapmean_pressure)
                call set_all(snapmean_pressure, snapmean_p(:))
                call insert(pod_state(i,k), snapmean_pressure, name="Snapmean"//trim(s_field%name)) 
                call deallocate(snapmean_pressure)
		if(deim) then 
		call vtk_write_state(filename=trim(simulation_name)//"_PODDEIMBasisRES"//trim(s_field%name),  &
                     index=i, model=trim(s_field%mesh%name), state=(/pod_state(i,k)/))
		     else
                call vtk_write_state(filename=trim(simulation_name)//"_PODBasis"//trim(s_field%name),  &
                     index=i, model=trim(s_field%mesh%name), state=(/pod_state(i,k)/))
	       endif
             enddo
          endif
       end do !j = 1, size(state(i)%scalar_fields)
       
!   stop 78     
    end do!k =1, size(state)
    
    deallocate(pod_state)
    
    
  end subroutine form_podstate_writeout_allvariables_differentmesh
  
  subroutine read_input_states(state)
    !!< Read the input states from the vtu dumps of the forward run.
    type(state_type), intent(out), dimension(:), allocatable :: state
    character(len=1024) :: filename

    integer :: dump_sampling_period, quadrature_degree
    integer :: i,j,k,total_dumps,stable_dumps
    type(vector_field), pointer :: velocity,vfield

    call get_option('/reduced_model/pod_basis_formation/dump_sampling_period',dump_sampling_period)
    call get_option('/geometry/quadrature/degree', quadrature_degree)

    print*, 'dump_sampling_period',dump_sampling_period

    !substract gyre_0.vtu
    total_dumps=count_dumps(dump_sampling_period)-1
!    total_dumps = 10
!    allocate(state(total_dumps))
!    stable_dumps=total_dumps-10
!    allocate(state(stable_dumps))
    print*, 'total_dumps', total_dumps

    allocate(state(1:total_dumps))

    call read_inputs_on_differentMesh(state, total_dumps,dump_sampling_period)
!    do i=1, total_dumps
!       !! Note that this won't work in parallel. Have to look for the pvtu in that case.
!       write(filename, '(a, i0, a)') trim(simulation_name)//'_', (i)*dump_sampling_period,".vtu" 
!       call vtk_read_state(filename, state(i), quadrature_degree)
!    end do
 
  end subroutine read_input_states

subroutine read_input_states_deim(state_deim)
    !!< Read the input states from the vtu dumps of the forward run.
    type(state_type), intent(out), dimension(:), allocatable :: state_deim
    character(len=1024) :: filename

    integer :: dump_sampling_period, quadrature_degree
    integer :: i,j,k,total_dumps,stable_dumps

    call get_option('/reduced_model/pod_basis_formation/dump_sampling_period',dump_sampling_period)
    call get_option('/geometry/quadrature/degree', quadrature_degree)

    ewrite(3,*) 'dump_sampling_period',dump_sampling_period

    !substract gyre_0.vtu
    total_dumps=count_dumps(dump_sampling_period)-1
    allocate(state_deim(total_dumps))

!    stable_dumps=total_dumps-10
!    allocate(state(stable_dumps))

    do i=1, total_dumps

       !! Note that this won't work in parallel. Have to look for the pvtu in that case.
     
       ! write(filename, '(a, i0, a)') trim(simulation_name)//'_DEIM_', (i)*dump_sampling_period,".vtu" 
       write(filename, '(a, i0, a)') trim(simulation_name)//'_PODDEIMRES_', (i+1)*dump_sampling_period,".vtu"                  
     
       call vtk_read_state(filename, state_deim(i), quadrature_degree)
       
       !! Note that we might need code in here to clean out unneeded fields.
     end do

  end subroutine read_input_states_deim

  function count_dumps(dump_sampling_period) result (count)
    !! Work out how many dumps we're going to read in.
    integer :: count,dump_sampling_period

    logical :: exists
    !      character(len=FILE_NAME_LEN) :: filename
    character(len=1024) :: filename

    count=1

    do 
       !! Note that this won't work in parallel. Have to look for the pvtu in that case.
   !    write(filename, '(a, i0, a)') trim(simulation_name)//'_', (count-1)*dump_sampling_period,".pvtu" 
       write(filename, '(a, i0, a)') trim(simulation_name)//'_', (count-1)*dump_sampling_period,".vtu" 
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

  subroutine read_inputs_on_differentMesh(ensemble_state, nrens,dump_sampling_period)

    type(state_type), dimension(:), intent(inout) :: ensemble_state
    !!number of ensembles
    integer, intent(in) :: nrens,dump_sampling_period
    type(vector_field), pointer :: velocity,vfield
    type(scalar_field), pointer :: pressure,sfield
    type(tensor_field), pointer :: tfield
    type(scalar_field), pointer :: free_surface
    type(scalar_field), pointer :: temperature
    type(scalar_field), pointer :: salinity
    character(len=255), dimension(:), allocatable :: filename
    character(len=255), dimension(:), allocatable :: filename1
    character(len=255), dimension(:), allocatable :: filename2
    character(len=255), dimension(:), allocatable :: filename3
    integer :: quadrature_degree,cp_no
    integer :: i, j, k, d,ifield,nfield,stat
    integer :: u_nodes, p_nodes, fs_nodes, t_nodes, s_nodes
    character(len = OPTION_PATH_LEN) :: simulation_name, mesh_name,prefix,postfix
    !!ensemble matrix mean
    real, dimension(:), allocatable :: meanA
    type(state_type), dimension(:,:), allocatable :: mean_state
    type(vector_field):: vfield2
    type(scalar_field):: sfield2
    type(tensor_field):: tfield2

    logical :: exists
    type(vector_field), pointer:: vfield_tmp
    type(scalar_field), pointer:: sfield_tmp
    type(tensor_field), pointer:: tfield_tmp
    type(mesh_type), pointer:: mesh_tmp
    type(state_type), dimension(:) , allocatable :: ensemble_state_tmp
    logical :: keep_initial_data,dont_write_files

    call get_option('/simulation_name',simulation_name)
    call get_option('/geometry/quadrature/degree', quadrature_degree)

    allocate(filename(nrens))
    prefix = simulation_name
    postfix = "checkpoint"
    do i = 1, nrens
       call nullify(ensemble_state(i)) 
       call set_option_path(ensemble_state(i), "/material_phase["//int2str(i-1)//"]")
    end do

       call insert_external_mesh( ensemble_state, save_vtk_cache = .true.)
       call insert_derived_meshes(ensemble_state)
       
       call allocate_and_insert_fields( ensemble_state )
       
       call initialise_prognostic_fields( ensemble_state, save_vtk_cache=.true., &
            initial_mesh=.true.)
       call set_prescribed_field_values(ensemble_state, initial_mesh=.true.)


    call print_state( ensemble_state(1) )
    print*,'**********'
    call print_state( ensemble_state(2) )

    do i = 2, nrens

       nfield = scalar_field_count( ensemble_state(1) )
       do ifield = 1, nfield 
          sfield => extract_scalar_field( ensemble_state(1), ifield )

          exists=.false.
          do j=1, scalar_field_count( ensemble_state(1) )
             sfield_tmp => extract_scalar_field( ensemble_state(i), trim(sfield%name), stat)
             mesh_name = trim(sfield%mesh%name)
             if (trim(sfield%name)==trim(sfield_tmp%name)) then
                exists=.true.
                exit
             end if
          end do
          if (exists) cycle

          do j=1,mesh_count(ensemble_state(1))
             mesh_tmp => extract_mesh( ensemble_state(i), trim(sfield%mesh%name)  )
             if (trim(mesh_name)==trim(mesh_tmp%name)) exit
          end do

!          print *, node_count(sfield),trim(sfield%name), trim(mesh_tmp%name)
          call allocate( sfield2,  mesh_tmp, trim(sfield%name)   )
          call set(sfield2, node_val(sfield,1) )  !!! this is used only for the constant field
          sfield2%option_path =  sfield%option_path
          call insert(ensemble_state(i), sfield2, trim(sfield%name))
          !option_path = trim( sfield%option_path )
        ! call insert(ensemble_state(i), sfield2, trim( option_path ))
       end do

!!! Vector
       !_-------
       nfield = vector_field_count( ensemble_state(1) )
!       print*,'&&&&vector',nfield
       do ifield = 1, nfield 
          vfield => extract_vector_field( ensemble_state(1), ifield )

          exists=.false.
          do j=1, vector_field_count( ensemble_state(1) )
             vfield_tmp => extract_vector_field( ensemble_state(i), trim(vfield%name), stat)
             mesh_name = trim(vfield%mesh%name)
             if (trim(vfield%name)==trim(vfield_tmp%name)) then
                exists=.true.
                exit
             end if
          end do
          if (exists) cycle

          print*, 'EXISTS',exists 
          do j=1,mesh_count(ensemble_state(1))
             mesh_tmp => extract_mesh( ensemble_state(i), trim(vfield%mesh%name)  )
             if (trim(mesh_name)==trim(mesh_tmp%name)) exit
          end do

          print *, node_count(vfield),trim(vfield%name)
          call allocate( vfield2,  vfield%dim, vfield%mesh, trim(vfield%name)   )
          call set(vfield2, node_val(vfield,1 )) !!! this is used only for the constant field
          vfield2%option_path =  vfield%option_path
          call insert(ensemble_state(i), vfield2, trim(vfield%name))
       end do

       nfield = tensor_field_count( ensemble_state(1) )
       print*,'&&&&tensor',nfield
       do ifield = 1, nfield 
          tfield => extract_tensor_field( ensemble_state(1), ifield )



          exists=.false.
          do j=1, tensor_field_count( ensemble_state(1) )
             tfield_tmp => extract_tensor_field( ensemble_state(i), trim(tfield%name), stat)
             mesh_name = trim(tfield%mesh%name)
             if (trim(tfield%name)==trim(tfield_tmp%name)) then
                exists=.true.
                exit
             end if
          end do
          if (exists) cycle

          print*, 'EXISTS',exists 
          do j=1,mesh_count(ensemble_state(1))
             mesh_tmp => extract_mesh( ensemble_state(i), trim(tfield%mesh%name)  )
             if (trim(mesh_name)==trim(mesh_tmp%name)) exit
          end do

          print *, node_count(tfield),trim(tfield%name)
          call allocate( tfield2, tfield%mesh, trim(tfield%name)   )
          call set(tfield2, node_val(tfield,1 )) !!! this is used only for the constant field
          tfield2%option_path =  tfield%option_path
          call insert(ensemble_state(i), tfield2, trim(tfield%name))
       end do
    end do

      print*,'11111'
    call print_state( ensemble_state(1) )
    print*,'**********'
    call print_state( ensemble_state(2) )


    do i = 1, nrens
!       cp_no = (i-1)*dump_sampling_period
       cp_no = i*dump_sampling_period

       call solution_on_checkpoints( (/ ensemble_state(i) /)  , prefix, postfix, cp_no)

       call initialise_prognostic_fields((/ensemble_state(i) /), save_vtk_cache=.true., &
            initial_mesh=.true.)

       call set_prescribed_field_values((/ensemble_state(i) /), initial_mesh=.true.)
     end do

!     do i =1,size(ensemble_state)
!       call vtk_write_state(filename=trim(simulation_name)//"_snapshot_after", index=i, state=(/ensemble_state(i)/))
!    enddo
!    stop 97


 

!  vfield => extract_vector_field( ensemble_state(50), 'Velocity', stat)
!    print*,nrens
!    print*,vfield%val(1,1:20)
!     stop 98



  end subroutine read_inputs_on_differentMesh

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
