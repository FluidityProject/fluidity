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
    implicit none

    type(state_type), dimension(:), pointer :: state
    type(state_type), dimension(:), allocatable :: ensemble_state
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

    call python_init()
    call read_command_line()
    !!state includes the fields i.c., b.c.and mesh info in flml file
    call populate_state(state)

    nrens=2

    !!call get_option("/enkf/assimilation_variable",format)
    format='FreeSurfaceHeight'

    !!ensemble_state from ensemble vtus
    call construct_ensemble_matrix(ensemble_state, A, nrens, global_ndim, nvar)


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
     call analysis6(A, E, S, meanD, global_ndim, nrens, nobs, .false.)
   else
      print* , 'EnKF: Computing local analysis'
      call local_analysis(nrens,nobs,global_ndim,nvar,A,state,measurement_state,meanD,E,S)
      !call local_analysis(nrens,nrobs,A,obs,modlon,modlat,depths,D,E,S,local_ndim)
   endif


contains

  !!Construct the ensemble matrix A
  subroutine construct_ensemble_matrix(ensemble_state, A, nrens, global_ndim, nvar)

   type(state_type), dimension(:), allocatable :: ensemble_state
   real, dimension(:,:),  allocatable, intent(out) :: A
   !!number of ensembles
   integer, intent(in) :: nrens
   integer, intent(out) :: global_ndim
   type(vector_field), pointer :: velocity
   type(scalar_field), pointer :: pressure
   type(scalar_field), pointer :: free_surface
   type(scalar_field), pointer :: temperature
   type(scalar_field), pointer :: salinity
   character(len=1024) :: filename
   integer :: quadrature_degree
   integer :: i, j, k, d
   integer, intent(out) :: nvar
   integer :: u_nodes, p_nodes, t_nodes, s_nodes
   character(len = OPTION_PATH_LEN) :: simulation_name

   call get_option('/simulation_name',simulation_name)
   call get_option('/geometry/quadrature/degree', quadrature_degree)
   allocate(ensemble_state(nrens))

   do i=1, nrens

      !! Note that this won't work in parallel. Have to look for the pvtu in that case.
      write(filename, '(a, i0, a)') trim(simulation_name)//'_', i, ".vtu"
      call vtk_read_state(filename, ensemble_state(i), quadrature_degree)

      !! Note that we might need code in here to clean out unneeded fields.
   end do

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
!!u,v,w,p,t,s: 3-D data

   nvar=0
   u_nodes=0
   p_nodes=0
   t_nodes=0
   s_nodes=0

   velocity => extract_vector_field(ensemble_state(1), "Velocity")   
   u_nodes=node_count(velocity)
   nvar=nvar+velocity%dim
   pressure => extract_scalar_field(ensemble_state(1), 'Pressure')
   p_nodes=node_count(pressure)
   nvar=nvar+1

   global_ndim=u_nodes*velocity%dim+p_nodes  

   !!Diagnostic field FreeSurface in ensemble vtus. Can it be on the surface mesh?
   !!free_surface => extract_scalar_field(ensemble_state(1), "FreeSurface")

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
    real, dimension(:,:), allocatable, intent(out) :: D
    real, dimension(:,:), allocatable, intent(out) :: E
    character(len=FIELD_NAME_LEN), intent(in) :: format
    type(mesh_type) :: measurement_mesh
 
    no_measurement_scalar_fields=1

    do k=1, no_measurement_scalar_fields
        
       !!assimilation_variable=measurement_variable
       !!call get_option("/enkf/assimilation_variable",format)
       !!format="FreeSurfaceHeight"      
 
       select case (format)
       case ("FreeSurfaceHeight")
          !!Surface data
          SSH = extract_scalar_field(state,"MeasurementSSH")
          nobs = node_count(SSH)  !!How about the four vertex?
          allocate(D(nobs, nrens))
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

       case ("SeaSurfaceTemperature")
          SST = extract_scalar_field(state,"MeasurementSST")
       case ("Salinity")
          Salinity = extract_scalar_field(state,"MeasurementSalinity")
       end select
    enddo
 
   end subroutine construct_measurement_matrix

  subroutine construct_HA(state, ensemble_state, measurement_state, nrens, nobs, D, S, meanD, meanS, format)
     type(state_type), dimension(:), intent(in) :: state
     type(state_type), dimension(:), allocatable, intent(in) :: ensemble_state
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

     allocate(S(nobs, nrens))
     allocate(meanS(nobs))

     select case (format)
     case("FreeSurfaceHeight")
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

     case("SeaSurfaceTemperature")
     !!Both are surface meshes for surface obs(SSH, SST).
     measurement_mesh = extract_mesh(state, "MeasurementMesh")
     input_mesh = extract_mesh(state, "InputMesh")

     case("Salinity")
     !!In situ data, 3D mesh
     measurement_mesh = extract_mesh(state, "MeasurementMesh")
     input_mesh = extract_mesh(state, "CoordinateMesh")

     end select

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

  end subroutine construct_HA

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

    write (0,*) "usage: enkf [-v n] <options_file>"
    write (0,*) ""
    write (0,*) "-v n sets the verbosity of debugging"
  end subroutine usage



  end program EnKF
