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
!    use fluids_module 
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
    use analysis6_module
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
    integer :: global_ndim

    call python_init()
    call read_command_line()
    !!state includes the fields i.c., b.c.and mesh info in flml file
    call populate_state(state)

    nrens=100

    !!call get_option("/enkf/assimilation_variable",format)
    format='FreeSurfaceHeight'

    !!ensemble_state from ensemble vtus
    call construct_ensemble_matrix(ensemble_state, A, nrens, format)


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

   if (l_global_analysis) then
      print *,'EnKF: Computing global analysis'
      !##print*,'EnKF: Size A',size(A)
      !call analysis(A, D, E, S, global_ndim, nrens, nrobs, .false.)
      !call analysis6(A, E, S, meanD, global_ndim, nrens, nrobs, .false.)
     !## call analysis6c(A, E, S, meanD, global_ndim, nrens, nrobs, .true.)
     call analysis6(A, E, S, meanD, global_ndim, nrens, nobs, .false.)
   else
      print* , 'EnKF: Computing local analysis'
      ! FIX IT!!!!!????
!      call local_analysis(nrens,nrobs,A,obs,modlon,modlat,depths,D,E,S,local_ndim)
   endif


contains

  !!Construct the ensemble matrix A
  subroutine construct_ensemble_matrix(ensemble_state, A, nrens, format)

   type(state_type), dimension(:), allocatable :: ensemble_state
   real, dimension(:,:),  allocatable, intent(out) :: A
   !!number of ensembles
   integer, intent(in) :: nrens
   type(vector_field), pointer :: velocity
   type(scalar_field), pointer :: pressure
   type(scalar_field), pointer :: free_surface
   character(len=1024) :: filename
   integer :: quadrature_degree
   integer :: i, j, k, p_nodes
   character(len = OPTION_PATH_LEN) :: simulation_name
   character(len = OPTION_PATH_LEN), intent(in) :: format

   call get_option('/simulation_name',simulation_name)
   call get_option('/geometry/quadrature/degree', quadrature_degree)
   allocate(ensemble_state(nrens))

   do i=1, nrens

      !! Note that this won't work in parallel. Have to look for the pvtu in that case.
      write(filename, '(a, i0, a)') trim(simulation_name)//'_', i, ".vtu"
      call vtk_read_state(filename, ensemble_state(i), quadrature_degree)

      !! Note that we might need code in here to clean out unneeded fields.
   end do

   !!call get_option("/enkf/assimilation_variable",format)
   !!format='FreeSurfaceHeight'
   select case (format)
   case('FreeSurfaceHeight')
      !!Diagnostic field FreeSurface in ensemble vtus. Can it be on the surface mesh?
      free_surface => extract_scalar_field(ensemble_state(1), "FreeSurface")
      p_nodes=node_count(pressure)
      allocate(A(p_nodes, nrens))
      do i=1, nrens
         free_surface => extract_scalar_field(ensemble_state(i), "FreeSurface")
         A(:,i)=field_val(pressure)
      enddo

   case('SeaSurfaceTemperature')

   case('Salinity')

   end select

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
             D(:,i)=E(:,i)+field_val(SSH)
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
     type(scalar_field) :: measurement
     real, dimension(:,:), allocatable, intent(inout) :: D
     real, dimension(:,:), allocatable, intent(out) :: S
     real, dimension(:), allocatable :: meanD
     real, dimension(:), allocatable :: meanS
     integer, intent(in) :: nrens
     integer, intent(in) :: nobs
     integer :: i, j

     select case (format)
     case("FreeSurfaceHeight")
     !!Both are surface meshes for surface obs(SSH, SST).
     measurement_mesh = extract_mesh(state, "MeasurementMesh")
     input_mesh = extract_mesh(state, "InputMesh")

     !!Observe ensemble to construct the matrix HA
     !!Interpolate the unstrutured input_mesh (surface mesh) to the measurement mesh
     !call linear_interpolation_state(old_state, new_state, map, different_domains)
     !call linear_interpolation_state(ensemble_state, measurement_state)

     case("SeaSurfaceTemperature")
     !!Both are surface meshes for surface obs(SSH, SST).
     measurement_mesh = extract_mesh(state, "MeasurementMesh")
     input_mesh = extract_mesh(state, "InputMesh")

     case("Salinity")
     !!In situ data, 3D mesh
     measurement_mesh = extract_mesh(state, "MeasurementMesh")
     input_mesh = extract_mesh(state, "CoordinateMesh")

     end select
     !!Compute innovation D'=D-HA                 ==> return in D
     allocate(S(nobs, nrens))
     allocate(meanS(nobs))
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
