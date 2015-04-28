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

#ifdef HAVE_PETSC_MODULES
  use petsc
#endif
  implicit none
#include "petsc_legacy.h"

  type(state_type), dimension(:), allocatable :: state,state_test
  type(state_type), dimension(:,:), allocatable :: pod_state

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

  call form_basis()


  deallocate(pod_state)
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

    call get_option(&
         '/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
    call get_option('/simulation_name',simulation_name)
    call read_input_states(state)

    call retrieve_snapshots(state, snapshots, u_nodes, p_nodes, snapmatrix_velocity, snapmatrix_pressure, &
                            & snapmean_velocity, snapmean_pressure)

    call form_svd(snapmatrix_velocity, snapmatrix_pressure,&
       & leftsvd_velocity, leftsvd_pressure, svdval_velocity, svdval_pressure, snapshots)

    call form_podstate(state,pod_state,leftsvd_velocity,leftsvd_pressure, snapmean_velocity, snapmean_pressure)

    do i=1,nsvd

       dump_no=i

       call vtk_write_state(filename=trim(simulation_name)//"_PODBasis", index=dump_no, state=pod_state(i,:))

       call deallocate(pod_state(i,:))

    enddo

    !! Produce updated flml file in which the execute_reduced_model
    !! option is set.
    call add_option("/reduced_model/execute_reduced_model",stat)
    call set_option('/simulation_name', trim(simulation_name)//'_POD') 
    call write_options(trim(simulation_name)//"_POD.flml")

  end subroutine form_basis


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

    allocate(snapmatrix_velocity(u_nodes,snapshots,dim))
    allocate(snapmatrix_pressure(p_nodes,snapshots))
    allocate(snapmean_velocity(u_nodes,dim))
    allocate(snapmean_pressure(p_nodes))

    do i = 1, snapshots
       velocity => extract_vector_field(state(i), "Velocity")
       pressure => extract_scalar_field(state(i), "Pressure")

       do d = 1, dim
          snapmatrix_velocity(:,i,d)=field_val(velocity,d)
       end do
       snapmatrix_pressure(:,i)=field_val(pressure)
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

    call get_option(&
         "/reduced_model/pod_basis_formation/pod_basis_count", POD_num) 

    allocate(pod_state(POD_num,1))
    call nullify(pod_state)

    do i = 1,POD_num

       pod_mesh => extract_mesh(state(1), "Mesh")

       all_meshes_same = .true.

       pod_xmesh => extract_mesh(state(1), "Mesh", stat)
       pod_umesh => extract_mesh(state(1), "Mesh", stat)
       pod_pmesh => extract_mesh(state(1), "Mesh", stat)

       pod_positions => extract_vector_field(state(1), "Coordinate")

!       call insert(pod_state(i,1), pod_xmesh, "Mesh")
       call insert(pod_state(i,1), pod_xmesh, "CoordinateMesh")
       call insert(pod_state(i,1), pod_umesh, "VelocityMesh")
       call insert(pod_state(i,1), pod_pmesh, "PressureMesh")
       call insert(pod_state(i,1), pod_positions, "Coordinate")

       velocity => extract_vector_field(state(1), "Velocity")

       call allocate(pod_velocity, velocity%dim, pod_umesh, "PODVelocity")
       call zero(pod_velocity)
       do d=1,velocity%dim
          call set_all(pod_velocity, d, leftsvd_u(:,i,d))
       end do
       call insert(pod_state(i,1), pod_velocity, name="PODVelocity")
       call deallocate(pod_velocity)

       call allocate(pod_pressure, pod_umesh, "PODPressure")
       call zero(pod_pressure)
       call set_all(pod_pressure, leftsvd_p(:,i))
       call insert(pod_state(i,1), pod_pressure, name="PODPressure")
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


  subroutine read_input_states(state)
    !!< Read the input states from the vtu dumps of the forward run.
    type(state_type), intent(out), dimension(:), allocatable :: state
    character(len=1024) :: filename

    integer :: dump_sampling_period, quadrature_degree
    integer :: i,j,k,total_dumps,stable_dumps

    call get_option('/reduced_model/pod_basis_formation/dump_sampling_period',dump_sampling_period)
    call get_option('/geometry/quadrature/degree', quadrature_degree)

    ewrite(3,*) 'dump_sampling_period',dump_sampling_period

    !substract gyre_0.vtu
    total_dumps=count_dumps(dump_sampling_period)-1
    allocate(state(total_dumps))

!    stable_dumps=total_dumps-10
!    allocate(state(stable_dumps))

    do i=1, total_dumps
!    do i=1, stable_dumps

       !! Note that this won't work in parallel. Have to look for the pvtu in that case.
       write(filename, '(a, i0, a)') trim(simulation_name)//'_', (i)*dump_sampling_period,".vtu" 

       call vtk_read_state(filename, state(i), quadrature_degree)

       !! Note that we might need code in here to clean out unneeded fields.

    end do

  end subroutine read_input_states

  function count_dumps(dump_sampling_period) result (count)
    !! Work out how many dumps we're going to read in.
    integer :: count,dump_sampling_period

    logical :: exists
    !      character(len=FILE_NAME_LEN) :: filename
    character(len=1024) :: filename

    count=1

    do 
       !! Note that this won't work in parallel. Have to look for the pvtu in that case.
       write(filename, '(a, i0, a)') trim(simulation_name)//'_', (count-1)*dump_sampling_period,".vtu" 
       inquire(file=trim(filename), exist=exists)
       if (.not. exists) then
          count=count -1
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

end program form_pod_basis
