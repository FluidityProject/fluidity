!    Copyright (C) 2010 Imperial College London and others.
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

module reduced_model_runtime
  use state_module
  use spud
  use vtk_interfaces
  use FLDebug
  use sparse_tools
  use sparse_tools_petsc
  use sparse_matrices_fields
  use fields
  use field_options
  use vector_tools
  implicit none
  
!  private

  public :: read_pod_basis, solve_momentum_reduced, solve_advection_diffusion_cg_reduced
  public :: project_reduced, project_reduced_t, project_full, project_full_t
  public :: pod_matrix_type, pod_rhs_type, pod_matrix_t_type, pod_rhs_t_type

  type(state_type), dimension(:), allocatable, save :: POD_state

  type pod_matrix_type
     real, dimension(:,:), allocatable :: val
     integer :: u_nodes, p_nodes, u_dim
  end type pod_matrix_type

  type pod_rhs_type
     real, dimension(:), allocatable :: val
     integer :: u_nodes, p_nodes, u_dim
  end type pod_rhs_type

  type pod_matrix_t_type
     real, dimension(:,:), allocatable :: val
     integer :: t_nodes
  end type pod_matrix_t_type

  type pod_rhs_t_type
     real, dimension(:), allocatable :: val
     integer :: t_nodes 
  end type pod_rhs_t_type


  interface allocate
     module procedure allocate_pod_matrix, allocate_pod_rhs, allocate_pod_matrix_t, allocate_pod_rhs_t
  end interface

  interface deallocate
     module procedure deallocate_pod_matrix, deallocate_pod_rhs, deallocate_pod_matrix_t, deallocate_pod_rhs_t
  end interface

  interface addto
     module procedure addto_pod_matrix_vector_vector, addto_pod_rhs_velocity
  end interface

  interface addto_p
     module procedure addto_pod_matrix_pressure, addto_pod_rhs_pressure, addto_pod_matrix_free_surface
  end interface

  interface addto_t
     module procedure addto_pod_matrix_t, addto_pod_rhs_t
  end interface

contains

  subroutine read_pod_basis(POD_state, state)
    !! Read the pod basis from the set of vtu files.

    character(len=1024) :: simulation_name, filename

    integer :: dump_period, quadrature_degree
    integer :: i,j,total_dumps

    type(state_type), dimension(:), allocatable :: POD_state
    type(state_type), dimension(:) :: state
    type(vector_field) :: podVelocity, newpodVelocity
    type(scalar_field) :: podPressure, newpodPressure, podTemperature, newpodTemperature
    type(mesh_type) :: VelocityMesh, PressureMesh, TemperatureMesh


    call get_option('/simulation_name', simulation_name)
    call get_option('/geometry/quadrature/degree', quadrature_degree)
    
    total_dumps=count_dumps(simulation_name)
    allocate(POD_state(total_dumps))

    VelocityMesh=extract_velocity_mesh(state)
    PressureMesh=extract_pressure_mesh(state)

    do i=1, total_dumps

       !! Note that this won't work in parallel. Have to look for the pvtu in that case.
       write(filename, '(a, i0, a)') trim(simulation_name)//'Basis_', i,".vtu" 

       call vtk_read_state(filename, POD_state(i), quadrature_degree)

       !! Note that we might need code in here to clean out unneeded fields.

       PODVelocity=extract_vector_field(POD_state(i),"PODVelocity")

       call allocate(newpodVelocity, podVelocity%dim, VelocityMesh, "PODVelocity")
       call remap_field(from_field=podVelocity, to_field=newpodVelocity)
       call insert(POD_state(i), newpodVelocity, "PODVelocity")
       call deallocate(newpodVelocity)

       PODPressure=extract_scalar_field(POD_state(i),"PODPressure")

       call allocate(newpodPressure, PressureMesh, "PODPressure")
       call remap_field(from_field=podPressure, to_field=newpodPressure)
       call insert(POD_state(i), newpodPressure, "PODPressure")
       call deallocate(newpodPressure)

       if(have_option('/material_phase::ocean/scalar_field::Temperature'))then

          TemperatureMesh=extract_mesh(state,"CoordinateMesh")
          PODTemperature=extract_scalar_field(POD_state(i),"PODTemperature")
          call allocate(newpodTemperature, TemperatureMesh, "PODTemperaturecsr_mult")
          call remap_field(from_field=podTemperature, to_field=newpodTemperature)
          call insert(POD_state(i), newpodTemperature, "PODTemperature")
          call deallocate(newpodTemperature)

       endif
    end do

  contains

    function count_dumps(simulation_name) result (count)
      !! Work out how many dumps we're going to read in.
      integer :: count
      character(len=*), intent(in) :: simulation_name
      
      logical :: exists
     
      character(len=1024) :: filename
      
      count=0
      
      do 
         !! Note that this won't work in parallel. Have to look for the pvtu in that case.
         write(filename, '(a, i0, a)') trim(simulation_name)//'Basis_', count+1,".vtu" 
         inquire(file=trim(filename), exist=exists)
         if (.not. exists) then
            exit
         end if
         
         count=count+1
      end do
      
      if (count==0) then
         FLExit("No POD.vtu files found!")
      end if
    end function count_dumps    

  end subroutine read_pod_basis

  subroutine solve_momentum_reduced(delta_u, delta_p, big_m, mom_rhs, ct_m, ct_rhs, timestep, dt, POD_state) 
    !!< Solve the momentum equation in reduced space.
    type(vector_field), intent(inout) :: delta_u
    type(scalar_field), intent(inout) :: delta_p
    type(state_type), dimension(:) :: POD_state

    type(petsc_csr_matrix), intent(inout) :: big_m
    type(block_csr_matrix), intent(in) :: ct_m
    type(vector_field), intent(in) :: mom_rhs
    type(scalar_field), intent(in) :: ct_rhs
    real, intent(in) :: dt

    type(vector_field), pointer :: POD_u
    type(scalar_field), pointer :: POD_p
    type(scalar_field) :: snapmean_pressure 

    type(pod_matrix_type) :: pod_matrix
    type(pod_rhs_type) :: pod_rhs

    integer :: i, j, d, d1, d2, u_nodes, p_nodes, POD_num, timestep

    real, dimension(:), allocatable :: pod_coef 
    real, dimension(:,:), allocatable :: pod_sol_velocity
    real, dimension(:), allocatable :: pod_sol_pressure

    POD_u=>extract_vector_field(POD_state(1), "PODVelocity")
    POD_p=>extract_scalar_field(POD_state(1), "PODPressure")

    POD_num=size(POD_state)

    allocate(pod_coef(POD_u%dim*POD_num+POD_num))
    pod_coef=0.0

    call project_reduced(big_m, mom_rhs, ct_m, ct_rhs, pod_matrix, pod_rhs, dt, POD_state)

!!call solver
!        call LEGS_POD(pod_matrix%val, POD_u%dim*POD_num+POD_num, pod_rhs%val, pod_coef) 
     
     call solve(pod_matrix%val, pod_rhs%val)
     pod_coef=pod_rhs%val
     print*,pod_coef
!!project back to velocity and pressure 

     call project_full(delta_u, delta_p, pod_sol_velocity, pod_sol_pressure, POD_state, pod_coef)
       
!!output to .vtu files
    
!       call form_podstate_solution(POD_state, pod_state_solution, pod_sol_velocity, pod_sol_pressure)

!       call get_option('/simulation_name', simulation_name)
       
!       call vtk_write_state(filename=trim(simulation_name)//"_sol", index=timestep, state=pod_state_solution(1,:))

     call deallocate(pod_matrix)
     call deallocate(pod_rhs)
     deallocate(pod_coef)

  end subroutine solve_momentum_reduced


  subroutine project_reduced(big_m, mom_rhs, ct_m, ct_rhs, pod_matrix, pod_rhs, dt, POD_state, fs_m)
    !!< project the momentum equation to reduced space.
    type(state_type), dimension(:) :: POD_state

    type(petsc_csr_matrix), intent(inout) :: big_m
    type(block_csr_matrix), intent(in) :: ct_m
    type(vector_field), intent(in) :: mom_rhs
    type(scalar_field), intent(in) :: ct_rhs
    real, intent(in) :: dt

    type(vector_field), pointer :: POD_u
    type(scalar_field), pointer :: POD_p
    type(scalar_field) :: snapmean_pressure

    type(vector_field), dimension(:), allocatable :: u_tmp
    type(vector_field) :: u_c, mom_rhs_tmp
    type(scalar_field) :: ct_tmp
    type(scalar_field) :: comp1, comp2

    type(pod_matrix_type), intent(out) :: pod_matrix
    type(pod_rhs_type), intent(out) :: pod_rhs

    integer :: i, j, d, d1, d2, u_nodes, p_nodes, POD_num
    real, dimension(:,:), allocatable :: pod_tmp

    type(csr_matrix), optional, intent(in) :: fs_m
    type(scalar_field) :: fs_tmp

    ewrite(1,*)'in project_reduced'
    print*,'in project_reduced'

    POD_u=>extract_vector_field(POD_state(1), "PODVelocity")
    POD_p=>extract_scalar_field(POD_state(1), "PODPressure")

    snapmean_pressure=extract_scalar_field(POD_state(1), "SnapmeanPressure")

    u_nodes=node_count(POD_u)
    p_nodes=node_count(POD_p)
    POD_num=size(POD_state)

    allocate(pod_tmp(POD_u%dim,POD_u%dim))
    allocate(u_tmp(POD_u%dim))
    pod_tmp=0.0

    do d1=1,POD_u%dim
       call allocate(u_tmp(d1), POD_u%dim, POD_u%mesh, "PODTmpU")
    end do
    call allocate(u_c, POD_u%dim, POD_u%mesh, "PODTmpComponent")
    call allocate(ct_tmp, POD_p%mesh, "PODTmpCT")
    if(present(fs_m))then
       call allocate(fs_tmp, POD_p%mesh, "PODTmpFS")
    endif

    call allocate(pod_matrix, POD_u, POD_p)
    call allocate(pod_rhs, POD_u, POD_p)

    call allocate(mom_rhs_tmp, POD_u%dim, POD_u%mesh, "mom_rhs_tmp")
!!Are we solving for p^{n+1}-p^n or p^{n+1}? Take p^{n+1}-p^n now.
    call mult_T(mom_rhs_tmp, ct_m, snapmean_pressure)

    !mom_rhs_tmp is constant
!    do d=1, POD_u%dim
!!       mom_rhs%val(d)%ptr=mom_rhs%val(d)%ptr-mom_rhs_tmp%val(d)%ptr
!       mom_rhs%val(d,:)=mom_rhs%val(d,:)-mom_rhs_tmp%val(d,:)
!    enddo

    do j=1, POD_num
       POD_u=>extract_vector_field(POD_state(j), "PODVelocity")
       POD_p=>extract_scalar_field(POD_state(j), "PODPressure")

!!summation of rows after multiplication       
       do d1=1,POD_u%dim
          do d2=1,POD_u%dim
             if (d1==d2) then
                call set(u_c, d2, POD_u)
             else
                call set(u_c, d2, 0.0)
             end if
          end do

          call mult(u_tmp(d1), big_m, u_c)
          call mult(ct_tmp, ct_m, u_c)

          !print*, 'u_c'
          !print*, u_c%val(1,1:50)
          print*, 'u_tmp'
          print*, u_tmp(d1)%val(:,1)
          !stop

!!calculations for pod_rhs
          call addto(POD_state, pod_rhs, j, d1, sum(dot_product(mom_rhs, u_c)))
          call addto_p(POD_state, POD_u, pod_rhs, j, dot_product(ct_rhs, POD_p), dt)
          print*,'pod_rhs'
          print*,pod_rhs%val

!          pod_rhs%val(j+(d1-1)*POD_num)=sum(dot_product(mom_rhs, u_c))
!          pod_rhs%val(j+POD_u%dim*POD_num)=dot_product(ct_rhs, POD_p)

            do i=1, POD_num
               POD_p=>extract_scalar_field(POD_state(i), "PODPressure")

                call addto_p(POD_state, pod_matrix, d1, i, j, dot_product(ct_tmp, POD_p))
!!transpose the block corresponding to ct_m
                pod_matrix%val(j+(d1-1)*POD_num, i+pod_matrix%u_dim*POD_num)=&
                               & pod_matrix%val(i+pod_matrix%u_dim*POD_num, j+(d1-1)*POD_num)
            end do
        end do

!!free surface matrix
        if(present(fs_m))then
        POD_p=>extract_scalar_field(POD_state(j), "PODPressure")
        call mult(fs_tmp, fs_m, POD_p)       
        do i=1, POD_num
           POD_p=>extract_scalar_field(POD_state(i), "PODPressure")
           call addto_p(POD_state, pod_matrix, i, j, dot_product(fs_tmp, POD_p))
        enddo
        endif


!!summation of collumns after multiplication
          do i=1, POD_num
             POD_u=>extract_vector_field(POD_state(i), "PODVelocity")

             do d1=1,POD_u%dim
                comp2=extract_scalar_field(POD_u,d1)

                do d2=1,POD_u%dim
                   comp1=extract_scalar_field(u_tmp(d2),d1)

                   pod_tmp(d1,d2)=dot_product(comp1,comp2)

                   call addto(POD_state, pod_matrix, d1, d2, i, j, pod_tmp)

                end do
             end do

         enddo
     enddo

     print*,pod_matrix%val
     !stop
 
     do d1=1,POD_u%dim
        call deallocate(u_tmp(d1))
     end do
     deallocate(u_tmp)
     call deallocate(ct_tmp)
     if(present(fs_m))then
        call deallocate(fs_tmp)
     endif
     call deallocate(u_c)
     call deallocate(mom_rhs_tmp)
     
  end subroutine project_reduced

  subroutine solve_advection_diffusion_cg_reduced(delta_t, matrix, rhs, timestep, POD_state)
    !!< Solve the advection_diffusion_cg equation in reduced space.
    type(scalar_field), intent(inout) :: delta_t
    type(state_type), dimension(:), intent(in) :: POD_state
    type(csr_matrix), intent(in) :: matrix
    type(scalar_field), intent(in) :: rhs

    type(scalar_field), pointer :: POD_t
    type(pod_matrix_t_type) :: pod_matrix_t
    type(pod_rhs_t_type) :: pod_rhs_t

    integer :: i, j, t_nodes, POD_num, timestep

    real, dimension(:), allocatable :: pod_coef
    real, dimension(:), allocatable :: pod_sol_pressure
    real, dimension(:), allocatable :: pod_sol_temperature
    type(state_type), dimension(:,:), allocatable :: pod_state_solution

    POD_t=>extract_scalar_field(POD_state(1), "PODTemperature")

    t_nodes=node_count(POD_t)
    POD_num=size(POD_state)

    allocate(pod_coef(POD_num))
    pod_coef=0.0

    call project_reduced_t(matrix, rhs, pod_matrix_t, pod_rhs_t, POD_state)

!!call solver
!    call LEGS_POD(pod_matrix%val, POD_num, pod_rhs%val, pod_coef) 

     call solve(pod_matrix_t%val, pod_rhs_t%val)
     pod_coef=pod_rhs_t%val
     print*,pod_coef
!!project back to velocity and pressure 

     call project_full_t(delta_t, pod_sol_temperature, POD_state, pod_coef)

!!output to .vtu files

!       call form_podstate_solution(POD_state, pod_state_solution, pod_sol_velocity, pod_sol_pressure)

!       call get_option('/simulation_name', simulation_name)

!       call vtk_write_state(filename=trim(simulation_name)//"_sol", index=timestep, state=pod_state_solution(1,:))

     call deallocate(pod_matrix_t)
     call deallocate(pod_rhs_t)
     deallocate(pod_coef)

  end subroutine solve_advection_diffusion_cg_reduced


  subroutine project_reduced_t(matrix, rhs, pod_matrix_t, pod_rhs_t, POD_state)
    !!< Project the advection_diffusion_cg equation to reduced space.
    type(state_type), dimension(:), intent(in) :: POD_state
    type(csr_matrix), intent(in) :: matrix
    type(scalar_field), intent(in) :: rhs

    type(scalar_field), pointer :: POD_t
    type(scalar_field) :: t_tmp

    type(pod_matrix_t_type) :: pod_matrix_t
    type(pod_rhs_t_type) :: pod_rhs_t

    integer :: i, j, t_nodes, POD_num

    POD_t=>extract_scalar_field(POD_state(1), "PODTemperature")

    t_nodes=node_count(POD_t)
    POD_num=size(POD_state)

    call allocate(t_tmp,POD_t%mesh, "PODTmpT")
    call allocate(pod_matrix_t, POD_t)
    call allocate(pod_rhs_t, POD_t)

    do j=1, POD_num

       POD_t=>extract_scalar_field(POD_state(j), "PODTemperature")
       call mult_T(t_tmp%val, matrix, POD_t%val)

       do i=1, POD_num
          POD_t=>extract_scalar_field(POD_state(i), "PODTemperature")

           !!calculations for pod_matrix
           call addto_t(pod_matrix_t, i, j, dot_product(t_tmp, POD_t))
           !!calculations for pod_rhs
           call addto_t(pod_rhs_t, i, dot_product(rhs, POD_t))
!          pod_rhs%val(i)=dot_product(rhs, POD_t)
       enddo
    enddo

     call deallocate(t_tmp)
  end subroutine project_reduced_t 

  subroutine allocate_pod_matrix(matrix,u,p)
    type(pod_matrix_type), intent(inout) :: matrix
    type(vector_field), intent(in) :: u
    type(scalar_field), intent(in) :: p
    integer POD_num

    call get_option('/reduced_model/pod_basis_formation/pod_basis_count', POD_num) 
     matrix%u_dim=u%dim

    allocate(matrix%val(POD_num*matrix%u_dim+POD_num, POD_num*matrix%u_dim+POD_num))

    matrix%val=0.0

  end subroutine allocate_pod_matrix

  subroutine allocate_pod_rhs(rhs,u,p)
    type(pod_rhs_type), intent(inout) :: rhs
    type(vector_field), intent(in) :: u
    type(scalar_field), intent(in) :: p
    integer POD_num

    call get_option('/reduced_model/pod_basis_formation/pod_basis_count', POD_num)
     rhs%u_dim=u%dim

    allocate(rhs%val(POD_num*rhs%u_dim+POD_num))

    rhs%val=0.0

  end subroutine allocate_pod_rhs

  subroutine allocate_pod_matrix_t(pod_matrix_t,t)
    type(pod_matrix_t_type), intent(inout) :: pod_matrix_t
    type(scalar_field), intent(in) :: t
    integer POD_num

    call get_option('/reduced_model/pod_basis_formation/pod_basis_count', POD_num)
    allocate(pod_matrix_t%val(POD_num, POD_num))

    pod_matrix_t%val=0.0

  end subroutine allocate_pod_matrix_t

  subroutine allocate_pod_rhs_t(rhs,t)
    type(pod_rhs_t_type), intent(inout) :: rhs
    type(scalar_field), intent(in) :: t
    integer POD_num

    call get_option('/reduced_model/pod_basis_formation/pod_basis_count', POD_num)
    allocate(rhs%val(POD_num))

    rhs%val=0.0

  end subroutine allocate_pod_rhs_t

  subroutine deallocate_pod_matrix(matrix)
    type(pod_matrix_type), intent(inout) :: matrix
    integer POD_num
 
    matrix%u_dim=0
    POD_num=0

    deallocate(matrix%val)

  end subroutine deallocate_pod_matrix

  subroutine deallocate_pod_rhs(rhs)
    type(pod_rhs_type), intent(inout) :: rhs
    integer POD_num

    rhs%u_dim=0
    POD_num=0

    deallocate(rhs%val)

  end subroutine deallocate_pod_rhs

  subroutine deallocate_pod_matrix_t(matrix)
    type(pod_matrix_t_type), intent(inout) :: matrix
    integer POD_num

    POD_num=0
    deallocate(matrix%val)

  end subroutine deallocate_pod_matrix_t

  subroutine deallocate_pod_rhs_t(rhs)
    type(pod_rhs_t_type), intent(inout) :: rhs
    integer POD_num

    POD_num=0
    deallocate(rhs%val)

  end subroutine deallocate_pod_rhs_t

  subroutine addto_pod_matrix_vector_vector(POD_state, matrix, d1, d2, i, j, value)
    type(pod_matrix_type), intent(inout) :: matrix
    type(state_type), dimension(:), intent(in) :: POD_state

    integer, intent(in) :: i, j, d1, d2
    real, dimension(:,:), allocatable :: value    
    integer ::  POD_num
    
    POD_num=size(POD_state)

          matrix%val(i+(d1-1)*POD_num, j+(d2-1)*POD_num) =&
                         matrix%val(i+(d1-1)*POD_num, j+(d2-1)*POD_num)+ &
                         value(d1,d2)  

  end subroutine addto_pod_matrix_vector_vector

  subroutine addto_pod_matrix_pressure(POD_state, matrix, d1, i, j, value) 
     type(pod_matrix_type), intent(inout) :: matrix
     type(state_type), dimension(:), intent(in) :: POD_state

     integer, intent(in) :: i, j, d1
     real :: value
     integer :: POD_num

     POD_num=size(POD_state)
     matrix%val(i+matrix%u_dim*POD_num, j+(d1-1)*POD_num) =&
                         matrix%val(i+matrix%u_dim*POD_num, j+(d1-1)*POD_num)+ &
                         value

  end subroutine addto_pod_matrix_pressure

  subroutine addto_pod_matrix_free_surface(POD_state, matrix, i, j, value)
     type(pod_matrix_type), intent(inout) :: matrix
     type(state_type), dimension(:), intent(in) :: POD_state

     integer, intent(in) :: i, j
     real :: value
     integer :: POD_num

     POD_num=size(POD_state)
     matrix%val(i+matrix%u_dim*POD_num, j+matrix%u_dim*POD_num) =&
                         matrix%val(i+matrix%u_dim*POD_num, j+matrix%u_dim*POD_num)+ &
                         value

  end subroutine addto_pod_matrix_free_surface

  subroutine addto_pod_matrix_t(matrix, i, j, value)
     type(pod_matrix_t_type), intent(inout) :: matrix

     integer, intent(in) :: i, j
     real :: value

     matrix%val(i,j) = matrix%val(i,j)+value

  end subroutine addto_pod_matrix_t


  subroutine addto_pod_rhs_velocity(POD_state, pod_rhs, j, d1, value)
    type(pod_rhs_type), intent(inout) :: pod_rhs
    type(state_type), dimension(:), intent(in) :: POD_state

    integer, intent(in) :: j, d1
    real :: value
    integer ::  POD_num

    POD_num=size(POD_state)

    pod_rhs%val(j+(d1-1)*POD_num)=pod_rhs%val(j+(d1-1)*POD_num)+value

  end subroutine addto_pod_rhs_velocity

  subroutine addto_pod_rhs_pressure(POD_state, pod_u, pod_rhs, j, value, dt)
    type(pod_rhs_type), intent(inout)          :: pod_rhs
    type(state_type), dimension(:), intent(in) :: POD_state
    type(vector_field), pointer                :: POD_u
    real, intent(in) :: dt
    integer, intent(in) :: j 

    real :: value
    integer ::  POD_num

    POD_num=size(POD_state)
         
    pod_rhs%val(j+POD_u%dim*POD_num)=pod_rhs%val(j+POD_u%dim*POD_num)+value/dt
          
  end subroutine addto_pod_rhs_pressure

  subroutine addto_pod_rhs_t(pod_rhs, i, value)
    type(pod_rhs_t_type), intent(inout) :: pod_rhs

    integer, intent(in) :: i
    real :: value

    pod_rhs%val(i)=pod_rhs%val(i)+value

  end subroutine addto_pod_rhs_t

    SUBROUTINE LEGS_POD(A,N,B,X)

! Subroutine to solve the equation A(N,N)*X(N) = B(N) with the
! partial-pivoting Gaussian elimination scheme.
!
      INTEGER I,J,N
      real A(N,N),B(N),X(N)
      INTEGER,DIMENSION(:),ALLOCATABLE ::INDX

      ALLOCATE(INDX(N))
      INDX=0
      X=0.
      ewrite(3,*) 'In LEGS'
      CALL ELGS(A,N,INDX)

      DO I = 1, N-1
        DO J = I+1, N
            B(INDX(J)) = B(INDX(J)) -A(INDX(J),I)*B(INDX(I))
        ENDDO     
      ENDDO
      X(N) = B(INDX(N))/A(INDX(N),N)
      DO I = N-1, 1, -1
        X(I) = B(INDX(I))
        DO J = I+1, N
          X(I) = X(I)-A(INDX(I),J)*X(J)
        ENDDO
          X(I) =  X(I)/A(INDX(I),I)
      ENDDO
    
      DEALLOCATE(INDX)
   END SUBROUTINE LEGS_POD

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE ELGS(A,N,INDX)
!
! Subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed
! matrix plus the pivoting element ratios below the diagonal in
! the output.  INDX(N) records the pivoting order.
      INTEGER I,J,N,K,itmp,INDX(N)
      real A(N,N),C(N),C1,PI,PI1,pj
! Initialize the index
!      ewrite(3,*) 'IN ELGS'
      do     I = 1, N! Was loop 50
        INDX(I) = I
      end do ! Was loop 50
!
! Find the rescaling factors, one from each row
!
      do    I = 1, N! Was loop 100
          C1= 0.0
      do    J = 1, N! Was loop 90
            C1 = AMAX1(C1,ABS(A(I,J)))
!            ewrite(3,*) C1
      end do ! Was loop 90
          C(I) = C1
      end do ! Was loop 100
!
! Search the pivoting (largest) element from each column
!
      do    J = 1, N-1! Was loop 200
        PI1 = 0.0
      do    I = J, N! Was loop 150
          PI = ABS(A(INDX(I),J))/C(INDX(I))
          IF (PI.GT.PI1) THEN
            PI1 = PI
            K   = I
          ELSE
          ENDIF
      end do ! Was loop 150
!
! Interchange the rows via INDX(N) to record pivoting order
!
!        ewrite(3,*) indx
        ITMP    = INDX(J)
        INDX(J) = INDX(K)
        INDX(K) = ITMP
      do    I = J+1, N! Was loop 170
          PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
          A(INDX(I),J) = PJ
!
! Modify other elements accordingly
!
      do    K = J+1, N! Was loop 160
            A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
      end do ! Was loop 160
      end do ! Was loop 170
      end do ! Was loop 200
!
      RETURN
      END SUBROUTINE ELGS
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      
      subroutine project_full(delta_u, delta_p, pod_sol_velocity, pod_sol_pressure, POD_state, pod_coef)

        type(vector_field), intent(inout) :: delta_u
        type(scalar_field), intent(inout) :: delta_p
        real, dimension(:,:), intent(out), allocatable :: pod_sol_velocity
        real, dimension(:), intent(out), allocatable :: pod_sol_pressure
        type(state_type), dimension(:), intent(in) :: POD_state
        real, dimension(:), intent(in) :: pod_coef

        type(vector_field), pointer :: POD_velocity
        type(scalar_field), pointer :: POD_pressure

        type(vector_field), pointer :: snapmean_velocity
        type(scalar_field), pointer :: snapmean_pressure

        integer :: d, POD_num, i, u_nodes, p_nodes, stat

        type(mesh_type), pointer :: pod_umesh, pod_pmesh

        POD_velocity=>extract_vector_field(POD_state(1), "PODVelocity")
        POD_pressure=>extract_scalar_field(POD_state(1), "PODPressure")

        snapmean_velocity=>extract_vector_field(POD_state(1), "SnapmeanVelocity")
        snapmean_pressure=>extract_scalar_field(POD_state(1), "SnapmeanPressure")

        u_nodes=node_count(POD_velocity)
        p_nodes=node_count(POD_pressure)

        allocate(pod_sol_velocity(u_nodes,POD_velocity%dim))
        allocate(pod_sol_pressure(p_nodes))
        pod_sol_velocity=0.0
        pod_sol_pressure=0.0

        POD_num=size(POD_state)

        do i=1,POD_num
!        do i=1,1

           POD_velocity=>extract_vector_field(POD_state(i), "PODVelocity")
           POD_pressure=>extract_scalar_field(POD_state(i), "PODPressure")

           do d=1,POD_velocity%dim
              pod_sol_velocity(:,d)=pod_sol_velocity(:,d)+pod_coef(i+(d-1)*POD_num)*POD_velocity%val(d,:)
           enddo

           pod_sol_pressure(:)=pod_sol_pressure(:)+pod_coef(i+POD_velocity%dim*POD_num)*POD_pressure%val(:)

           pod_umesh => extract_mesh(pod_state(1), "Mesh", stat)
           pod_pmesh => extract_mesh(pod_state(1), "Mesh", stat)

           do d=1, POD_velocity%dim
              call set_all(delta_u, d, pod_sol_velocity(:,d))
           end do

           call set_all(delta_p, pod_sol_pressure(:))

        enddo

        deallocate(pod_sol_velocity)
        deallocate(pod_sol_pressure)


      end subroutine project_full
 
      subroutine project_full_t(delta_t, pod_sol_temperature, POD_state, pod_coef)
        type(scalar_field), intent(inout) :: delta_t
        real, dimension(:), intent(out), allocatable :: pod_sol_temperature
        type(state_type), dimension(:), intent(in) :: POD_state
        real, dimension(:), intent(in) :: pod_coef

        type(scalar_field), pointer :: POD_temperature
        type(scalar_field), pointer :: snapmean_temperature

        integer :: d, POD_num, i, t_nodes, stat

        type(mesh_type), pointer :: pod_tmesh

        POD_temperature=>extract_scalar_field(POD_state(1), "PODTemperature")
        snapmean_temperature=>extract_scalar_field(POD_state(1), "SnapmeanTemperature")

        t_nodes=node_count(POD_temperature)

        allocate(pod_sol_temperature(t_nodes))
        pod_sol_temperature=0.0

        POD_num=size(POD_state)

        do i=1,POD_num

           POD_temperature=>extract_scalar_field(POD_state(i), "PODTemperature")

           pod_sol_temperature(:)=pod_sol_temperature(:)+pod_coef(i)*POD_temperature%val(:)

           pod_tmesh => extract_mesh(pod_state(1), "Mesh", stat)

           call set_all(delta_t, pod_sol_temperature(:))

        enddo

        deallocate(pod_sol_temperature)

      end subroutine project_full_t

  subroutine form_podstate_solution(pod_state, pod_state_solution, pod_sol_velocity, pod_sol_pressure)

    type(state_type), intent(in), dimension(:) :: pod_state
    type(state_type), intent(out), dimension(:,:), allocatable :: pod_state_solution

    real, dimension(:,:) :: pod_sol_velocity
    real, dimension(:) :: pod_sol_pressure

    type(mesh_type), pointer :: pod_xmesh, pod_umesh, pod_pmesh, pod_mesh
    type(vector_field), pointer :: pod_positions, velocity

    type(vector_field) :: pod_velocity
    type(vector_field), pointer :: snapmean_velocity
    type(scalar_field) :: pod_pressure

    integer :: quadrature_degree, nonods
    integer :: u_nodes, p_nodes, POD_num
    integer :: i,j,k,total_dumps,stat,dim,d
    logical :: all_meshes_same

    allocate(pod_state_solution(1,1))

    call get_option(&
         "/reduced_model/pod_basis_formation/pod_basis_count", POD_num)

    call nullify(pod_state_solution)

       pod_mesh => extract_mesh(pod_state(1), "Mesh")

       all_meshes_same = .true.

       pod_xmesh => extract_mesh(pod_state(1), "Mesh", stat)
       pod_umesh => extract_mesh(pod_state(1), "Mesh", stat)

       pod_pmesh => extract_mesh(pod_state(1), "Mesh", stat)
       pod_positions => extract_vector_field(pod_state(1), "Coordinate")

       call insert(pod_state_solution(1,1), pod_xmesh, "Mesh")
       call insert(pod_state_solution(1,1), pod_xmesh, "CoordinateMesh")
       call insert(pod_state_solution(1,1), pod_umesh, "VelocityMesh")
       call insert(pod_state_solution(1,1), pod_pmesh, "PressureMesh")
       call insert(pod_state_solution(1,1), pod_positions, "Coordinate")

       velocity => extract_vector_field(pod_state(1), "PODVelocity")
       snapmean_velocity=>extract_vector_field(pod_state(1), 'SnapmeanVelocity')
       dim=velocity%dim

       call allocate(pod_velocity, velocity%dim, pod_umesh, "PODVelocity")
       call zero(pod_velocity)
       do d=1,dim
          call set_all(pod_velocity, d, pod_sol_velocity(:,d))
       end do
       call insert(pod_state_solution(1,1), pod_velocity, name="PODVelocity")

       call allocate(pod_pressure, pod_umesh, "PODPressure")
       call zero(pod_pressure)
       call set_all(pod_pressure, pod_sol_pressure(:))
       call insert(pod_state_solution(1,1), pod_pressure, name="PODPressure")

       pod_sol_velocity=0.0
       pod_sol_pressure=0.0

  end subroutine form_podstate_solution 



end module reduced_model_runtime

