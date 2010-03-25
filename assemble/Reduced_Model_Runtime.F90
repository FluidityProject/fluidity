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
!    C.Pain@Imperial.ac.uk
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
  
  private

  public :: read_pod_basis, solve_momentum_reduced

  type(state_type), dimension(:), allocatable, save :: POD_state

  type pod_matrix_type
     real, dimension(:,:), allocatable :: val
     integer :: u_nodes, p_nodes, u_dim
  end type pod_matrix_type

  type pod_rhs_type
     real, dimension(:), allocatable :: val
     integer :: u_nodes, p_nodes, u_dim
  end type pod_rhs_type

  interface allocate
     module procedure allocate_pod_matrix, allocate_pod_rhs
  end interface

  interface deallocate
     module procedure deallocate_pod_matrix, deallocate_pod_rhs
  end interface

  interface addto
     module procedure addto_pod_matrix_vector_vector, addto_pod_rhs_velocity
  end interface

  interface addto_p
     module procedure addto_pod_matrix_pressure, addto_pod_rhs_pressure
  end interface

contains

  subroutine read_pod_basis(POD_state, state)
    !! Read the pod basis from the set of vtu files.

    character(len=1024) :: simulation_name, filename

    integer :: dump_period, quadrature_degree
    integer :: i,j,k,total_dumps

    type(state_type), dimension(:), allocatable :: POD_state
    type(state_type), dimension(:) :: state
    type(vector_field) :: podVelocity, newpodVelocity
    type(scalar_field) :: podPressure, newpodPressure
    type(mesh_type) :: VelocityMesh, PressureMesh


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

  subroutine solve_momentum_reduced(delta_u, delta_p, big_m, mom_rhs, ct_m, ct_rhs, timestep, POD_state) 
    !!< Solve the momentum equation in reduced space.
    type(vector_field), intent(inout) :: delta_u
    type(scalar_field), intent(inout) :: delta_p
    type(state_type), dimension(:) :: POD_state

    type(petsc_csr_matrix), intent(inout) :: big_m
    type(block_csr_matrix), intent(in) :: ct_m
    type(vector_field), intent(in) :: mom_rhs
    type(scalar_field), intent(in) :: ct_rhs

    type(vector_field), pointer :: POD_u
    type(scalar_field), pointer :: POD_p

    type(vector_field), dimension(:), allocatable :: u_tmp
    type(vector_field) :: u_c
    type(scalar_field) :: ct_tmp
    type(scalar_field) :: comp1, comp2

    type(pod_matrix_type) :: pod_matrix
    type(pod_rhs_type) :: pod_rhs

    integer :: i, j, d1, d2, u_nodes, p_nodes, POD_num, timestep
    real, dimension(:,:), allocatable :: pod_tmp

    real, dimension(:), allocatable :: pod_coef 
    real, dimension(:,:), allocatable :: pod_sol_velocity
    real, dimension(:), allocatable :: pod_sol_pressure
    type(state_type), dimension(:,:), allocatable :: pod_state_solution
    character(len=1024) :: simulation_name

    POD_u=>extract_vector_field(POD_state(1), "PODVelocity")
    POD_p=>extract_scalar_field(POD_state(1), "PODPressure")

    u_nodes=node_count(POD_u)
    p_nodes=node_count(POD_p)
    POD_num=size(POD_state)

    allocate(pod_coef(POD_u%dim*POD_num+POD_num))
    allocate(pod_tmp(POD_u%dim,POD_u%dim))
    allocate(u_tmp(POD_u%dim))
    pod_coef=0.0
    pod_tmp=0.0

    do d1=1,POD_u%dim
       call allocate(u_tmp(d1), POD_u%dim, POD_u%mesh, "PODTmpU")
    end do
    call allocate(u_c, POD_u%dim, POD_u%mesh, "PODTmpComponent")
    call allocate(ct_tmp,POD_p%mesh, "PODTmpCT")
    
    call allocate(pod_matrix, POD_u, POD_p) 
    call allocate(pod_rhs, POD_u, POD_p)

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

!!calculations for pod_rhs

          call addto(POD_state, pod_rhs, j, d1, sum(dot_product(mom_rhs, u_c)))
          call addto_p(POD_state, POD_u, pod_rhs, j, dot_product(ct_rhs, POD_p))

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

!!call solver
!        call LEGS_POD(pod_matrix%val, POD_u%dim*POD_num+POD_num, pod_rhs%val, pod_coef) 
     
     call solve(pod_matrix%val, pod_rhs%val)
     pod_coef=pod_rhs%val

!!project back to velocity and pressure 

     call project(delta_u, delta_p, pod_sol_velocity, pod_sol_pressure, POD_state, pod_coef)
       
!!output to .vtu files
    
!       call form_podstate_solution(POD_state, pod_state_solution, pod_sol_velocity, pod_sol_pressure)

!       call get_option('/simulation_name', simulation_name)
       
!       call vtk_write_state(filename=trim(simulation_name)//"_sol", index=timestep, state=pod_state_solution(1,:))

     call deallocate(pod_matrix)
     call deallocate(pod_rhs)
     do d1=1,POD_u%dim
        call deallocate(u_tmp(d1))
     end do
     deallocate(u_tmp)
     call deallocate(ct_tmp)
     deallocate(pod_coef)
     call deallocate(u_c)

  end subroutine solve_momentum_reduced

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

  subroutine addto_pod_rhs_velocity(POD_state, pod_rhs, j, d1, value)
    type(pod_rhs_type), intent(inout) :: pod_rhs
    type(state_type), dimension(:), intent(in) :: POD_state

    integer, intent(in) :: j, d1
    real :: value
    integer ::  POD_num

    POD_num=size(POD_state)

    pod_rhs%val(j+(d1-1)*POD_num)=pod_rhs%val(j+(d1-1)*POD_num)+value

  end subroutine addto_pod_rhs_velocity

  subroutine addto_pod_rhs_pressure(POD_state,  pod_u, pod_rhs, j, value)
    type(pod_rhs_type), intent(inout) :: pod_rhs
    type(state_type), dimension(:), intent(in) :: POD_state
    type(vector_field), pointer, intent(in) :: POD_u

    integer, intent(in) :: j 
    real :: value
    integer ::  POD_num

    POD_num=size(POD_state)
         
    pod_rhs%val(j+POD_u%dim*POD_num)=pod_rhs%val(j+POD_u%dim*POD_num)+value
          
  end subroutine addto_pod_rhs_pressure


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

      
      subroutine project(delta_u, delta_p, pod_sol_velocity, pod_sol_pressure, POD_state, pod_coef)

        type(vector_field), intent(inout) :: delta_u
        type(scalar_field), intent(inout) :: delta_p
        real, dimension(:,:), intent(out), allocatable :: pod_sol_velocity
        real, dimension(:), intent(out), allocatable :: pod_sol_pressure
        type(state_type), dimension(:), intent(in) :: POD_state
        real, dimension(:), intent(in) :: pod_coef

        type(vector_field), pointer :: POD_velocity, POD_u
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

           POD_velocity=>extract_vector_field(POD_state(i), "PODVelocity")
           POD_pressure=>extract_scalar_field(POD_state(i), "PODPressure")

           do d=1,POD_velocity%dim
              pod_sol_velocity(:,d)=pod_sol_velocity(:,d)+pod_coef(i+(d-1)*POD_num)*POD_velocity%val(d)%ptr
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


      end subroutine project


  subroutine form_podstate_solution(pod_state, pod_state_solution, pod_sol_velocity, pod_sol_pressure)

    type(state_type), intent(in), dimension(:) :: pod_state
    type(state_type), intent(out), dimension(:,:), allocatable :: pod_state_solution

    real, dimension(:,:) :: pod_sol_velocity
    real, dimension(:) :: pod_sol_pressure

    type(mesh_type), pointer :: pod_xmesh, pod_umesh, pod_pmesh, pmesh, pod_mesh
    type(element_type) :: pod_xshape, pod_ushape, pod_pshape
    type(vector_field), pointer :: pod_positions, velocity
    type(scalar_field), pointer :: pressure

    type(vector_field) :: pod_velocity
    type(vector_field), pointer :: snapmean_velocity
    type(scalar_field) :: pod_pressure

    real, dimension(:), allocatable :: x,y,z

    character(len=1024) :: filename
    character(len = FIELD_NAME_LEN) :: field_name

    integer :: u_nodes, p_nodes, POD_num
    integer :: dump_period, quadrature_degree,nonods
    integer :: i,j,k,nod,total_dumps,stat,dim,f,d
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

