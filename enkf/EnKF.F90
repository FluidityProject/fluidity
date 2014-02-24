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
  program EnKF
    use signals
    use spud
    use populate_state_module
    use timeloop_utilities
    use sparsity_patterns_meshes
    use diagnostic_fields_wrapper
    use interpolation_manager
    use interpolation_module
    use python_state
    use vtk_interfaces
    use quadrature
    use shape_functions
    use bulk_parameterisations
    use global_parameters, only: option_path_len, current_time, dt, FIELD_NAME_LEN
    use m_random
    use analysis_module
    use m_local_analysis
    use interpolation_ensemble_state_on_supermesh
    use fldebug
    use fields_allocates
    use mesh_files

    implicit none

    type(state_type), dimension(:), pointer :: state
    type(state_type), dimension(:), allocatable :: ensemble_state_velocity
    type(state_type), dimension(:), allocatable :: ensemble_state_pressure
    type(state_type), dimension(:), allocatable :: ensemble_state
    type(state_type), dimension(:,:), allocatable :: analysis_state
    type(state_type), dimension(:), allocatable, save :: Measurement_state
    type(state_type), dimension(:), pointer :: Assemble_Matrix_A
    logical :: l_global_analysis
    character(len=OPTION_PATH_LEN) :: format
    !!ensemble matrix
    real, dimension(:,:), allocatable :: A
    !!number of ensembles
    integer :: nrens
    !!number of measurements
    integer :: nobs
    real, dimension(:,:), allocatable :: E
    real, dimension(:,:), allocatable :: D
    real, dimension(:,:), allocatable :: S
    real, dimension(:), allocatable :: meanD, meanS
    integer :: global_ndim, local_ndim
    integer :: nvar
    type(vector_field), pointer :: velocity
    type(scalar_field), pointer :: pressure
    type(scalar_field), pointer :: free_surface
    type(scalar_field), pointer :: temperature
    type(scalar_field), pointer :: salinity
    integer :: u_nodes, p_nodes, fs_nodes, t_nodes, s_nodes
    integer :: i, j, dump_no
    character(len = OPTION_PATH_LEN) :: simulation_name
    character(len=255) :: filename,name
    !!ensemble matrix mean analysis
    real, dimension(:), allocatable :: meanA_analysis
    type(state_type), dimension(:,:), allocatable :: mean_state
    character(len=OPTION_PATH_LEN) :: prefix,postfix

    print*,'in EnKF'
    call python_init()
    call read_command_line()
    !!state includes the fields i.c., b.c.and mesh info in flml file
    call populate_state(state)
    call get_option('/simulation_name',simulation_name)

    nrens=2

    !!call get_option("/enkf/assimilation_variable",format)
    format='FreeSurfaceHeight'

    !!ensemble_state from ensemble vtus
    call construct_ensemble_matrix(ensemble_state, ensemble_state_velocity, ensemble_state_pressure, A, nrens, global_ndim, nvar)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read measurements and store in d
! Construct observation perturbations        ==> return in E
! Construct ensemble of measurements D=d+E
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call construct_measurement_matrix(state, measurement_state, nrens, nobs, D, E, format)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Observe ensemble to construct the matrix HA
! Compute innovation D'=D-HA                 ==> return in D
! Compute mean(HA) 
! Compute HA'=HA-mean(HA)                    ==> return in S  !! need to modify -- will be calculated in Fluidity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call construct_HA(state, ensemble_state, measurement_state, nrens, nobs, D, S, meanS, meanD, format)
    !!using probe.py to get the data file *_detectors.dat
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute global analysis or local analysis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   print*, 'EnKF: Ready for analysis'
   
   l_global_analysis=.true.   

   if (l_global_analysis) then
      print *,'EnKF: Computing global analysis'
      !##print*,'EnKF: Size A',size(A)
      !call analysis(A, D, E, S, global_ndim, nrens, nrobs, .false.)
      !call analysis6(A, E, S, meanD, global_ndim, nrens, nrobs, .false.)
      !## call analysis6c(A, E, S, meanD, global_ndim, nrens, nrobs, .true.)

      !call analysis6(A, 0.001*E, S, meanD, global_ndim, nrens, nobs, .false.)
      call analysis(A, D, 0.001*E, S, global_ndim, nrens, nobs, .false.)
   else
      print* , 'EnKF: Computing local analysis'
      call local_analysis(nrens,nobs,global_ndim,nvar,A,state,measurement_state,meanD,E,S)
      !call local_analysis(nrens,nrobs,A,obs,modlon,modlat,depths,D,E,S,local_ndim)
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output as vtu files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do i=1,nrens
      write(filename,'(a,i0)')'analysis',i
      open(20,file=filename)
      write(20,*)A(:,i)
      close(20)
   enddo

!100 continue
    do i=1, nrens   
       print*,'writing__',trim(simulation_name)//"_analysis_",i,'.vtu'
       call form_analysis_state(ensemble_state, analysis_state, A(:,i))
       call vtk_write_state(filename=trim(simulation_name)//"_analysis", index=i, state=analysis_state(1,:))
       call deallocate(analysis_state)
    enddo

    allocate(meanA_analysis(global_ndim))
    meanA_analysis=0.0
    do i=1, nrens
       meanA_analysis=meanA_analysis+A(:,i)
    enddo
    meanA_analysis=meanA_analysis/nrens

    print*,'writing__',"ensemble_mean_analysis.vtu"
    call form_analysis_state(ensemble_state, mean_state, meanA_analysis)
    call vtk_write_state(filename="ensemble_mean_analysis", state=mean_state(1,:))
    call deallocate(mean_state)
    deallocate(meanA_analysis)

contains

  !!Construct the ensemble matrix A
  subroutine construct_ensemble_matrix(ensemble_state_old, ensemble_state_velocity_old, ensemble_state_pressure_old, A, nrens, global_ndim, nvar)

    type(state_type), dimension(:), allocatable :: ensemble_state
    type(state_type), dimension(:), allocatable :: ensemble_state_velocity
    type(state_type), dimension(:), allocatable :: ensemble_state_pressure
    type(state_type), dimension(:), allocatable :: ensemble_state_old,ensemble_state_old2
    type(state_type), dimension(:), allocatable :: ensemble_state_velocity_old
    type(state_type), dimension(:), allocatable :: ensemble_state_pressure_old
    real, dimension(:,:),  allocatable, intent(out) :: A
    !!number of ensembles
    integer, intent(in) :: nrens
    integer, intent(out) :: global_ndim
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
    integer :: quadrature_degree
    integer :: i, j, k, d,ifield,nfield,stat
    integer, intent(out) :: nvar
    integer :: u_nodes, p_nodes, fs_nodes, t_nodes, s_nodes
    character(len = OPTION_PATH_LEN) :: simulation_name, mesh_name,vtu_filename,prefix,postfix
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

    integer ::nmeshes,mxnods
    logical :: from_file, extruded
    character(len=OPTION_PATH_LEN) :: mesh_path, mesh_file_name,&
         mesh_file_format, from_file_path
    type(vector_field), dimension(:), allocatable :: initial_positions
    type(vector_field) :: out_positions



    print*,'In construct_ensemble_matrix'
    call get_option('/simulation_name',simulation_name)
    call get_option('/geometry/quadrature/degree', quadrature_degree)
    !allocate(ensemble_state(nrens))

    allocate(filename(nrens))
    allocate(filename1(nrens))
    allocate(filename2(nrens))


    !call interpolation_ensembles_on_supermesh(ensemble_state, ensemble_state_old, filename)
    !call deallocate(ensemble_state_old)

    !print*,'project the velocity field'
    !allocate(ensemble_state_velocity_old(nrens))
    !do i=1, nrens
    !   write(filename1(i), '(a, i0, a)') trim(simulation_name)//'_VelocityMesh_', i, ".vtu"
    !   call vtk_read_state(filename1(i), ensemble_state_velocity_old(i), quadrature_degree)
    !enddo

    !allocate(ensemble_state_velocity(nrens))
    !call interpolation_ensembles_on_supermesh(ensemble_state_velocity, ensemble_state_velocity_old, filename3)
    !call deallocate(ensemble_state_velocity_old)

    !   allocate(ensemble_state_pressure_old(nrens))

    call get_option('/simulation_name',simulation_name)
    call get_option('/geometry/quadrature/degree', quadrature_degree)

    prefix = simulation_name
    postfix = "checkpoint"


    !  calculate supermesh
    !------------------------
    !call set_option('/mesh_adaptivity/hr_adaptivity/maximum_number_of_nodes', mxnods, stat=stat)
    mxnods = 10000
    call compute_pseudo_supermesh_from_mesh(out_positions, mxnods=mxnods)

    ! allocate the new ensemble_state based on the supermesh
    !--------------------------------------------------------
    allocate(ensemble_state(nrens))
    call allocate_new_ensemble_state(out_positions, ensemble_state)

    ! allocate old ensemble_state
    !-----------------------------
    allocate(ensemble_state_old(1))
    call nullify(ensemble_state_old(1)) 
    call set_option_path(ensemble_state_old(1), "/material_phase["//int2str(1-1)//"]")

    call insert_external_mesh(ensemble_state_old, save_vtk_cache = .true.)
    !    call insert_external_mesh(ensemble_state, save_vtk_cache = .true.)
    call insert_derived_meshes(ensemble_state_old )
    !   call insert_derived_meshes( ensemble_state )

    call allocate_and_insert_fields(ensemble_state_old)
    call initialise_prognostic_fields(ensemble_state_old, save_vtk_cache=.true., &
         initial_mesh=.true.)
    call set_prescribed_field_values( ensemble_state_old, initial_mesh=.true.)
!    call add_missing_field_in_state(ensemble_state, ensemble_state_old)
      call print_state( ensemble_state(1) )
print*,'********************************************'
      call print_state( ensemble_state_old(1) )

!stop 123    

    call interpolation_ensembles_on_supermesh( ensemble_state, (/ensemble_state_old(1)/), 1)

       velocity => extract_vector_field(ensemble_state(1), "Velocity")
       pressure => extract_scalar_field(ensemble_state(1), "Pressure")
       free_surface => extract_scalar_field(ensemble_state(1), "FreeSurface")
       print*,node_count(velocity)
       print*,node_count(pressure)
       print*,node_count(free_surface)
       velocity => extract_vector_field(ensemble_state_old(1), "Velocity")
       pressure => extract_scalar_field(ensemble_state_old(1), "Pressure")
       free_surface => extract_scalar_field(ensemble_state_old(1), "FreeSurface")
       print*,node_count(velocity)
       print*,node_count(pressure)
       print*,node_count(free_surface)

       deallocate(ensemble_state_old)


    ! interpolate the old ensemble_state onto the supermesh
    !-----------------------------------------------------
    allocate(ensemble_state_old2(1))
    ensemble_state_loop: do i = 2, nrens

       call update_mesh_options(prefix,postfix,i,"gmsh")
       call nullify(ensemble_state_old2(1)) 
       call set_option_path(ensemble_state_old2(1), "/material_phase[0]")
       call insert_external_mesh(ensemble_state_old2, save_vtk_cache = .true.)
       
       call insert_derived_meshes(ensemble_state_old2)
       
       call allocate_and_insert_fields( ensemble_state_old2 )
       call print_state( ensemble_state_old2(1) )
       
       call solution_on_checkpoints( (/ ensemble_state_old2(1) /)  , prefix, postfix, i)
       
       call initialise_prognostic_fields(ensemble_state_old2, save_vtk_cache=.true., &
            initial_mesh=.true.)
       
       call set_prescribed_field_values(ensemble_state_old2, initial_mesh=.true.)

       !!call interpolation_ensembles_on_supermesh(ensemble_state_new, ensemble_state_old)
       call interpolation_ensembles_on_supermesh( ensemble_state, (/ensemble_state_old2(1)/), i)

       velocity => extract_vector_field(ensemble_state(1), "Velocity")
       pressure => extract_scalar_field(ensemble_state(1), "Pressure")
       free_surface => extract_scalar_field(ensemble_state(1), "FreeSurface")
       print*,node_count(velocity)
       print*,node_count(pressure)
       print*,node_count(free_surface)
       velocity => extract_vector_field(ensemble_state_old2(1), "Velocity")
       pressure => extract_scalar_field(ensemble_state_old2(1), "Pressure")
       free_surface => extract_scalar_field(ensemble_state_old2(1), "FreeSurface")
       print*,node_count(velocity)
       print*,node_count(pressure)
       print*,node_count(free_surface)
print*,'after interpolation_ensembles_on_supermesh'

    end do ensemble_state_loop
    deallocate(ensemble_state_old2)

    !   call vtk_write_state(filename="PressureMesh_test", state=ensemble_state_pressure)
    !   stop 12

!       do i =1, nrens
!          write(filename(i), '(a, i0, a)') trim(simulation_name)//'_', i, ".vtu"
!          call vtk_read_state(filename(2), ensemble_state(i), quadrature_degree)
!          
!          write(filename1(i), '(a, i0, a)') trim(simulation_name)//'_VelocityMesh_', i, ".vtu"
!          call vtk_read_state(filename1(i), ensemble_state_velocity(i), quadrature_degree) 
!       enddo

       !   !! Note that this won't work in parallel. Have to look for the pvtu in that case.
       !   write(filename2(i), '(a, i0, a)') trim(simulation_name)//'_PressureMesh_', i, ".vtu"
       !   !write(filename, '(a, i0, a)') 'munk_gyre-ns_', i, ".vtu"
       !   call vtk_read_state(filename2(i), ensemble_state_pressure(i), quadrature_degree)

       !   !! Note that we might need code in here to clean out unneeded fields.

    !call interpolation_ensembles_on_supermesh(ensemble_state_new, ensemble_state_pressure, filename)


    !HYCOM variables in ensemble matrix A
    !      real u(nx,ny,nz)                         ! 3-D u-velocity
    !      real v(nx,ny,nz)                         ! 3-D v-velocity
    !      real d(nx,ny,nz)                         ! 3-D Layer thickness
    !      real t(nx,ny,nz)                         ! 3-D Temperature
    !      real s(nx,ny,nz)                         ! 3-D Salinity 
    !      real ub(nx,ny)                           ! 2-D barotropic u-velocity
    !      real vb(nx,ny)                           ! 2-D barotropic v-velocity
    !      real pb(nx,ny)                           ! 2-D barotropic pressure

    !!FLUIDITY
    !!u,v,w,p,fs,t,s: 3-D data

    nvar=0
    u_nodes=0
    p_nodes=0
    t_nodes=0
    s_nodes=0

!    velocity => extract_vector_field(ensemble_state_velocity(1), "Velocity")   
    velocity => extract_vector_field(ensemble_state(1), "Velocity") 
    u_nodes=node_count(velocity)
    nvar=nvar+velocity%dim
    pressure => extract_scalar_field(ensemble_state(1), 'Pressure')
    p_nodes=node_count(pressure)
    nvar=nvar+1

    !!Diagnostic field FreeSurface in ensemble vtus. Can it be on the surface mesh?
    free_surface => extract_scalar_field(ensemble_state(1), "FreeSurface")
    fs_nodes=node_count(free_surface)
    nvar=nvar+1

    global_ndim=u_nodes*velocity%dim+p_nodes+fs_nodes

    if(have_option("/material_phase[0]/scalar_field::Temperature/prognostic"))then
       temperature => extract_scalar_field(ensemble_state(1), "Temperature")
       t_nodes=node_count(temperature)
       global_ndim=global_ndim+t_nodes
       nvar=nvar+1
    endif
    if(have_option("/material_phase[0]/scalar_field::Salinity/prognostic"))then
       salinity => extract_scalar_field(ensemble_state(1), "Salinity")
       s_nodes=node_count(salinity)
       global_ndim=global_ndim+s_nodes
       nvar=nvar+1
    endif

    print*,'u_nodes=',u_nodes
    print*,'p_nodes=',p_nodes
    print*,'fs_nodes=',fs_nodes
    print*,'t_nodes=',t_nodes
    print*,'s_nodes=',s_nodes
    print*,'global_ndim=',global_ndim
    print*,'nvar=',nvar

    allocate(A(global_ndim, nrens))
    do k=1, nrens
       velocity => extract_vector_field(ensemble_state(k), "Velocity")
       do d=1, velocity%dim
          i=1+u_nodes*(d-1)
          j=u_nodes*d
          !!velocity%val(d,:).ne.field_val(velocity,d)
          A(i:j,k)=velocity%val(d,:)
       enddo

       pressure => extract_scalar_field(ensemble_state(k), 'Pressure')
       i=j+1
       j=j+p_nodes
       !!pressure%val=field_val(pressure)
       A(i:j,k)=pressure%val

       free_surface => extract_scalar_field(ensemble_state(k), "FreeSurface")
       i=j+1
       j=j+fs_nodes
       A(i:j,k)=free_surface%val

       if(have_option("/material_phase[0]/scalar_field::Temperature/prognostic"))then
          temperature => extract_scalar_field(ensemble_state(k), "Temperature")
          i=j+1
          j=j+t_nodes
          A(i:j,k)=temperature%val
       endif

       if(have_option("/material_phase[0]/scalar_field::Salinity/prognostic"))then
          salinity => extract_scalar_field(ensemble_state(k), "Salinity")
          i=j+1
          j=j+s_nodes
          A(i:j,k)=salinity%val
       endif
    enddo

    allocate(meanA(global_ndim))   
    meanA=0.0 
    do i=1,nrens
       meanA(:)=meanA(:)+A(:,i)
    enddo

    meanA=meanA/nrens

    print*,'writing__',"ensemble_mean.vtu"
    call form_analysis_state(ensemble_state, mean_state, meanA)
    call vtk_write_state(filename="ensemble_mean", state=mean_state(1,:))
    !call deallocate(mean_state)
    !deallocate(meanA)
    stop

  end subroutine construct_ensemble_matrix

   subroutine construct_measurement_matrix(state, measurement_state, nrens, nobs, D, E, format)
   !!Construct the Measuement Matrix D 
    type(state_type), dimension(:), allocatable, intent(out) :: measurement_state
    type(state_type), dimension(:), intent(inout) :: state
    type(scalar_field) :: SSH, SST, Salinity
    type(scalar_field) :: perturb_SSH
    integer :: no_measurement_scalar_fields   ! how to get this -- need to discuss 
    integer :: i, j, k
    integer, intent(in) ::nrens
    integer, intent(out) :: nobs
    !!measurement matrix 
    real, dimension(:), allocatable :: obs
    real, dimension(:,:), allocatable, intent(out) :: D
    real, dimension(:,:), allocatable, intent(out) :: E
    character(len=FIELD_NAME_LEN), intent(in) :: format
    type(mesh_type) :: measurement_mesh

    print*,'In construct_measurement_matrix' 
    no_measurement_scalar_fields=1

    do k=1, no_measurement_scalar_fields
        
       !!assimilation_variable=measurement_variable
       !!call get_option("/enkf/assimilation_variable",format)
       !!format="FreeSurfaceHeight"      
 
       select case (format)
       case ("FreeSurfaceHeight")
          goto 20
          !!Surface data
          SSH = extract_scalar_field(state,"MeasurementSSH")
          nobs = node_count(SSH)  !!How about the four vertex?
          allocate(D(nobs, nrens))
          allocate(E(nobs, nrens))
          !d(:) = field_val(SSH)
          !!Generate observation perturbation matrix E
          call randn(nobs*nrens, E)
          do i=1, nrens
             D(:,i)=E(:,i)+SSH%val
          enddo
          !!Construct measurement_state using E for field value and state for mesh
          allocate(measurement_state(nrens))
          measurement_mesh = extract_mesh(state, "MeasurementMesh")

          do i=1, nrens
             call insert(measurement_state(i), measurement_mesh, "MeasurementMesh")

             call allocate(perturb_SSH, measurement_mesh, "PerturbedMeasurement")
             call zero(perturb_SSH)
             call set_all(perturb_SSH, D(:,i))
             call insert(measurement_state(i), perturb_SSH, name="PerturbedMeasurement")
             call deallocate(perturb_SSH)
          enddo
          !!reading observation data
20        nobs=21
          allocate(obs(nobs))
          open(20,file='obs.dat')
          read(20,*)obs
          close(20)
          !obs(:)=(-0.012462608359172147,0.006281770372309347,0.019493943092273976,0.0213375024168786,0.011891879675802512,-0.003866361663057322,-0.01757216090470086)
          allocate(D(nobs, nrens))
          allocate(E(nobs, nrens))
          call randn(nobs*nrens, E)
          do i=1, nrens
             D(:,i)=0.001*E(:,i)+obs(:)
          enddo
          !print*,D
          !stop

       case ("SeaSurfaceTemperature")
          SST = extract_scalar_field(state,"MeasurementSST")
       case ("Salinity")
          Salinity = extract_scalar_field(state,"MeasurementSalinity")
       end select
    enddo
 
   end subroutine construct_measurement_matrix

  subroutine construct_HA(state, ensemble_state, measurement_state, nrens, nobs, D, S, meanD, meanS, format)
     type(state_type), dimension(:), intent(in) :: state
     type(state_type), dimension(:), intent(in) :: ensemble_state
     type(state_type), dimension(:), allocatable, intent(in) :: measurement_state
     type(mesh_type) :: measurement_mesh, input_mesh
     character(len=OPTION_PATH_LEN), intent(in) :: format
     integer, intent(in) :: nrens
     integer, intent(in) :: nobs
     type(scalar_field), pointer :: measurement
     type(scalar_field), pointer :: SSH
     type(vector_field) :: old_position, new_position
     real, dimension(:,:), allocatable, intent(inout) :: D
     real, dimension(:,:), allocatable, intent(out) :: S
     real, dimension(:), allocatable :: meanD
     real, dimension(:), allocatable :: meanS
     integer :: quad_degree
     integer :: i, j

     print*,'In construct_HA'
     allocate(S(nobs, nrens))
     allocate(meanS(nobs))
     allocate(meanD(nobs))

     select case (format)
     case("FreeSurfaceHeight")
     goto 30
     !!Both are surface meshes for surface obs(SSH, SST).
     !measurement_mesh = extract_mesh(state, "MeasurementMesh")
     !input_mesh = extract_mesh(state, "InputMesh")

     !!Observe ensemble to construct the matrix HA
     !!Interpolate the unstrutured input_mesh (surface mesh) to the measurement mesh
     !sub linear_interpolation_state(old_state, new_state, map, different_domains)
     !call linear_interpolation(ensemble_state, measurement_state)

     SSH = extract_scalar_field(state, "SeaSurfaceHeight")
     old_position = extract_vector_field(state, "Coordinate")
     call get_option("/geometry/quadrature/degree", quad_degree)
     new_position = read_triangle_serial('MeasurementMesh', quad_degree=quad_degree)

     do i=1, nrens
        measurement => extract_scalar_field(measurement_state(i), "PerturbedMeasurement")
        !!sub linear_interpolation_scalars(old_fields, old_position, new_fields, new_position, map)
        !!sub quadratic_interpolation_qf(old_fields, old_position, new_fields, new_position)
        !!sub subroutine cubic_interpolation_cf_scalar(old_fields, old_position, new_fields, new_position)
        !!out:new_fields
        call linear_interpolation(SSH, old_position, measurement, new_position)
        S(:,i)=measurement%val
     enddo     
30   print*,'before reading detectors.dat'
     open(30,file='detectors.dat')
     do i=1,nrens
        read(30,*)S(:,i)    
     enddo
     print*,S(:,1)
     close(30)

     case("SeaSurfaceTemperature")
     !!Both are surface meshes for surface obs(SSH, SST).
     measurement_mesh = extract_mesh(state, "MeasurementMesh")
     input_mesh = extract_mesh(state, "InputMesh")

     case("Salinity")
     !!In situ data, 3D mesh
     measurement_mesh = extract_mesh(state, "MeasurementMesh")
     input_mesh = extract_mesh(state, "CoordinateMesh")

     end select

     !print*,D(:,1)
     !!Compute innovation D'=D-HA=D-S             ==> return in D
     D=D-S 
     !!Compute mean(HA) 
     meanS=0.0
     do i=1,nobs
        meanS(:)=meanS(:)+S(i,:)
     enddo
     meanS=(1.0/float(nrens))*meanS
     !!Compute HA'=HA-mean(HA)                    ==> return in S 
     do j=1,nrens
        S(:,j)=S(:,j)-meanS(:)
     enddo

     !!Compute meanD
     meanD=0.0
     do i=1,nobs
        meanD(:)=meanD(:)+D(i,:)
     enddo
     meanD=(1.0/float(nrens))*meanD

     !stop
  end subroutine construct_HA


  subroutine form_analysis_state(ensemble_state, analysis_state, A)

    type(state_type), intent(in), dimension(:) :: ensemble_state
    type(state_type), intent(out), dimension(:,:), allocatable :: analysis_state

    real, intent(in), dimension(:) :: A

    type(mesh_type), pointer :: analysis_xmesh, analysis_umesh, analysis_pmesh, pmesh, analysis_mesh
    type(element_type) :: analysis_xshape, analysis_ushape, analysis_pshape
    type(vector_field), pointer :: analysis_positions
    type(vector_field), pointer :: velocity
    type(scalar_field), pointer :: pressure, free_surface

    type(vector_field) :: analysis_velocity
    type(scalar_field) :: analysis_pressure, analysis_free_surface

    integer :: i, j, k
    integer :: stat, u_nodes, p_nodes, fs_nodes
    logical :: all_meshes_same

    !character(len = OPTION_PATH_LEN) :: simulation_name

    !call get_option('/simulation_name',simulation_name)
    allocate(analysis_state(1,1))
    call nullify(analysis_state)

print*,'before analysis_mesh'
       analysis_mesh => extract_mesh(ensemble_state(1), "CoordinateMesh")
print*,'after analysis_mesh'

       all_meshes_same = .false.

       analysis_xmesh => extract_mesh(ensemble_state(1), "CoordinateMesh", stat)
       analysis_umesh => extract_mesh(ensemble_state(1), "VelocityMesh", stat)
       analysis_pmesh => extract_mesh(ensemble_state(1), "PressureMesh", stat)

       analysis_positions => extract_vector_field(ensemble_state(1), "Coordinate")
!       call insert(pod_state(i,1), pod_xmesh, "Mesh")
       call insert(analysis_state(1,1), analysis_xmesh, "CoordinateMesh")
       call insert(analysis_state(1,1), analysis_umesh, "VelocityMesh")
       call insert(analysis_state(1,1), analysis_pmesh, "PressureMesh")
       call insert(analysis_state(1,1), analysis_positions, "Coordinate")
 print*,'after insert mesh'
 
       velocity => extract_vector_field(ensemble_state(1), "Velocity")
       pressure => extract_scalar_field(ensemble_state(1), "Pressure")
       free_surface => extract_scalar_field(ensemble_state(1), "FreeSurface")
       u_nodes=node_count(velocity)
       p_nodes=node_count(pressure)
       fs_nodes=node_count(free_surface)

       call allocate(analysis_velocity, velocity%dim, analysis_umesh, "Velocity")
print*,'after allocate velocity'
       call zero(analysis_velocity)
       do j=1,velocity%dim
          analysis_velocity%val(j,:)=A((j-1)*u_nodes+1:j*u_nodes)
       end do
print*,analysis_velocity%val(1,100:150)
       call insert(analysis_state(1,1), analysis_velocity, name="Velocity")
print*,'333'
       call deallocate(analysis_velocity)
print*,'after insert velocity'

       call allocate(analysis_pressure, analysis_pmesh, "Pressure")
       call zero(analysis_pressure)
       analysis_pressure%val(:)=A(velocity%dim*u_nodes+1:velocity%dim*u_nodes+p_nodes)
print*,analysis_pressure%val(100:150)
       call insert(analysis_state(1,1), analysis_pressure, name="Pressure")
       call deallocate(analysis_pressure)

       call allocate(analysis_free_surface, analysis_pmesh, "FreeSurface")
       call zero(analysis_free_surface)
       analysis_free_surface%val(:)=A(velocity%dim*u_nodes+p_nodes+1:velocity%dim*u_nodes+p_nodes+fs_nodes)
       call insert(analysis_state(1,1), analysis_free_surface, name="FreeSurface")
       call deallocate(analysis_free_surface)

  end subroutine form_analysis_state

  subroutine read_command_line()
    implicit none
    ! Read the input filename.
    character(len=255) :: argument, filename
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

    write (0,*) "usage: enkf [-v n] <options_file>"
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
  
    subroutine update_mesh_options(prefix,postfix,ensemble_no,format)
      !!< Updates the mesh.vtu
      !!< options path path

      character(len = 255) :: filename,mesh_file_name, simulation_name,path,from_file_path
      integer, intent(in)::ensemble_no
      character(len=OPTION_PATH_LEN),  intent(in):: prefix
      character(len=*),  intent(in):: format
      character(len=OPTION_PATH_LEN), optional, intent(in):: postfix
      character(len = 255) :: mesh_name

      integer :: stat, i, nmeshes    
      logical :: from_file, extruded


      nmeshes=option_count("/geometry/mesh")

      do i = 0, nmeshes-1  ! do while seemed to break, don't know why
         path="/geometry/mesh["//int2str(i)//"]"
         from_file_path = trim(path) // "/from_file"
         from_file = have_option(from_file_path)
         if (.not. from_file) then
            from_file_path = trim(path) // "/from_mesh/extrude/checkpoint_from_file"
            extruded = have_option(from_file_path)
         else
            extruded = .false.
         end if
         if(from_file .or. extruded) then
            call get_option(trim(path)//"/name", mesh_name)
            call delete_option(trim(from_file_path))
            mesh_file_name = trim(prefix)// "_" // trim(mesh_name)
            mesh_file_name = trim(mesh_file_name) // "_" // int2str(ensemble_no)
            mesh_file_name = trim(mesh_file_name) // "_" // trim(postfix)
            call set_option_attribute(trim(from_file_path) // "/file_name", trim(mesh_file_name), stat)
            if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
               FLAbort("Failed to set mesh filename when ensemble field with option path " // trim(path))
            end if
            call set_option_attribute(trim(from_file_path) // "/format/name", trim(format), stat)
            if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
               FLAbort("Failed to set mesh format when ensemble field with option path " // trim(path))
            end if
         end if        
         if(from_file .or. extruded) then
            call get_option(trim(from_file_path) // "/file_name",mesh_file_name) 
            print*,'mesh_name:: ',mesh_file_name
         endif
      end do

    end subroutine update_mesh_options

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



  end program EnKF
