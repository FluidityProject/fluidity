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
  use reduced_projection
  implicit none
  
!  private

  public :: read_pod_basis, read_pod_basis_differntmesh 
  type(state_type), dimension(:,:,:), allocatable, save :: POD_state

  
 
contains

  subroutine read_pod_basis(POD_state, state)
    
    !! Read the pod basis from the set of vtu files.
   
    character(len=1024) :: simulation_name, filename

    integer :: dump_period, quadrature_degree
    integer :: i,j,total_dumps,POD_num,k,nfield
    ! type(state_type), dimension(:,:,:), allocatable :: POD_state
    type(state_type), dimension(:,:,:), allocatable :: POD_state,POD_state_p
    type(state_type), dimension(:) :: state
    type(vector_field) :: podVelocity, newpodVelocity 
    type(scalar_field) :: podPressure, newpodPressure, podTemperature, newpodTemperature ,snapmean_pressure
    type(mesh_type) :: VelocityMesh, PressureMesh, TemperatureMesh
    
    type(scalar_field), pointer :: pres
    integer :: pod_pnodes,pod_unodes,p_nodes,u_nodes
   ! type(mesh_type) ,pointer ::pmesh
    call get_option('/simulation_name', simulation_name)
    call get_option('/geometry/quadrature/degree', quadrature_degree)
    call get_option(&
         "/reduced_model/pod_basis_formation/pod_basis_count", POD_num) 
 
    total_dumps=POD_num!count_dumps(simulation_name)
    nfield = vector_field_count( state(1) )+scalar_field_count( state(1) )
   ! allocate(POD_state(total_dumps))
   ! allocate(POD_state_p(total_dumps))
    allocate(pod_state(POD_num,nfield,size(state)))
    allocate(pod_state_p(POD_num,nfield,size(state)))
    VelocityMesh=extract_velocity_mesh(state)
    PressureMesh=extract_pressure_mesh(state)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !PressureMesh=extract_pressure_mesh(POD_state(1))
    !PressureMesh=extract_velocity_mesh(state)
    flag=1
    do k =1, size(state)
    do i=1, POD_num

       !! Note that this won't work in parallel. Have to look for the pvtu in that case.
      
       if (have_option("/reduced_model/Smolyak").or.have_option("/reduced_model/Produce_Smolyak_Data")) then
           write(filename, '(a, i0, a)') trim(simulation_name)//'BasisVelocity_', i,".vtu" 
       else
        write(filename, '(a, i0, a)') trim(simulation_name)//'_PODBasisVelocity_', i,".vtu" 
       endif
       call vtk_read_state(filename, POD_state(i,1,k), quadrature_degree)
        if (have_option("/reduced_model/Smolyak").or.have_option("/reduced_model/Produce_Smolyak_Data")) then
          write(filename, '(a, i0, a)') trim(simulation_name)//'BasisPressure_', i,".vtu" 
        else
         write(filename, '(a, i0, a)') trim(simulation_name)//'_PODBasisPressure_', i,".vtu" 
         endif
       call vtk_read_state(filename, POD_state_p(i,2,k), quadrature_degree)
       snapmean_pressure=extract_scalar_field(POD_state_p(i,2,k),"SnapmeanPressure")
       call insert(pod_state(i,2,k), snapmean_pressure, name="SnapmeanPressure")
       !! Note that we might need code in here to clean out unneeded fields.
      
      ! pmesh =>extract_mesh(POD_state_p(i), "PressureMesh")



       PODVelocity=extract_vector_field(POD_state(i,1,k),"PODVelocity")

       call allocate(newpodVelocity, podVelocity%dim, VelocityMesh, "PODVelocity")
       call remap_field(from_field=podVelocity, to_field=newpodVelocity)
       call insert(POD_state(i,1,k), newpodVelocity, "PODVelocity")
       call deallocate(newpodVelocity)
      !print *, 'ggggggggggggggggggggggggggggggggggg'
    
        PODPressure=extract_scalar_field(POD_state_p(i,2,k),"PODPressure")   !!!!!!!!!!!!!!!!!!!
       !PODPressure=extract_scalar_field(state,"Pressure")   !testbytravers
 
       call allocate(newpodPressure, PressureMesh, "PODPressure")   
       !call allocate(newpodPressure, podPre%dim,PressureMesh, "PODPressure")

           !!!!!!!!!!!!!!!!!!!!!!!!!!test
      !  velo => extract_vector_field(state, "Velocity")
      !  pres => extract_scalar_field(state, "Pressure")           
     !   p_nodes=node_count(pres)
      !  u_nodes=node_count(velo)            
      !  print*,'ppppppppppppppppppppppnodes',p_nodes
     !   print*,'uuuuuuuuuuuuuuuuuuuuuuuuunodes',u_nodes
        
       
     ! pod_pnodes=node_count(PODPressure)  
     ! pod_unodes=node_count(PODVelocity)
    !  print*, "podddddddddddddddddddddddpressure", pod_pnodes
    !  print*, "podddddddddddddddddddddddvelocity", pod_unodes


       call remap_field(from_field=podPressure, to_field=newpodPressure)     
       call insert(POD_state(i,2,k), newpodPressure, "PODPressure")      
       print *,'7'
       call deallocate(newpodPressure)
       print *,'remap_field(from_field=podPressure, to_field=newpodPressure)  pass '
       if(have_option('/material_phase::ocean/scalar_field::Temperature'))then

          TemperatureMesh=extract_mesh(state,"CoordinateMesh")
          PODTemperature=extract_scalar_field(POD_state(i,1,k),"PODTemperature")
          call allocate(newpodTemperature, TemperatureMesh, "PODTemperaturecsr_mult")
          call remap_field(from_field=podTemperature, to_field=newpodTemperature)
          call insert(POD_state(i,1,k), newpodTemperature, "PODTemperature")
          call deallocate(newpodTemperature)

       endif
    end do
     enddo

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


 subroutine read_pod_basis_deimres(POD_state_deim, deim_state_Res)
    !! Read the pod basis from the set of vtu files.

    character(len=1024) :: simulation_name, filename

    integer :: dump_period, quadrature_degree
    integer :: i,j,total_dumps

    type(state_type), dimension(:), allocatable :: POD_state_deim
    type(state_type), dimension(:) :: deim_state_Res
    type(vector_field) :: podVelocity, newpodVelocity
    type(scalar_field) :: podPressure, newpodPressure, podTemperature, newpodTemperature
    type(mesh_type) :: VelocityMesh, PressureMesh, TemperatureMesh
     
    velo => extract_vector_field(deim_state_Res, "Velocity") 
    call get_option('/simulation_name', simulation_name)
    call get_option('/geometry/quadrature/degree', quadrature_degree)
    
    total_dumps=count_dumps(simulation_name)
    allocate(POD_state_deim(total_dumps))
    VelocityMesh=extract_velocity_mesh(deim_state_Res)
    PressureMesh=extract_pressure_mesh(deim_state_Res)
    flag=1
    do i=1, total_dumps

       !! Note that this won't work in parallel. Have to look for the pvtu in that case.
       write(filename, '(a, i0, a)') trim(simulation_name)//'DEIMBasisRESVelocity_', i,".vtu" 

       call vtk_read_state(filename, POD_state_deim(i), quadrature_degree)

       !! Note that we might need code in here to clean out unneeded fields.

       PODVelocity=extract_vector_field(POD_state_deim(i),"Velocity")

       call allocate(newpodVelocity, podVelocity%dim, VelocityMesh, "PODVelocity")
       call remap_field(from_field=podVelocity, to_field=newpodVelocity)
       call insert(POD_state_deim(i), newpodVelocity, "PODVelocity")
       call deallocate(newpodVelocity)

      ! write(filename, '(a, i0, a)') trim(simulation_name)//'DEIMBasisRESPressure_', i,".vtu" 
     !  call vtk_read_state(filename, POD_state_deim(i), quadrature_degree)
     !  PODPressure=extract_scalar_field(POD_state_deim(i),"Pressure")

     !  call allocate(newpodPressure, PressureMesh, "PODPressure")
           !!!!!!!!!!!!!!!!!!!!!!!!!!test  rmse
        
                           ! pres => extract_scalar_field(state, "Pressure")           
                             !  p_nodes=node_count(pres)
                                 !u_nodes=node_count(velo) 
     !  call remap_field(from_field=podPressure, to_field=newpodPressure)
      ! call insert(POD_state_deim(i), newpodPressure, "PODPressure")
     !  call deallocate(newpodPressure)
 
    end do 
    
    
      allocate(phi(total_dumps*podVelocity%dim))  !actual is deim_number. need modificaiton later   
   !  allocate(phi(2*podVelocity%dim))
     open(unit=110,file='indices')
     read(110,*)(phi(i), i=1,total_dumps*podVelocity%dim)  ! deim_number 
   !  read(110,*)(phi(i), i=1,2*podVelocity%dim)   
     close(110)
 
    print *, ' phiphiphiphi', phi

    
   
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
         !write(filename, '(a, i0, a)') trim(simulation_name)//'BasisDEIM_', count+1,".vtu" 
          write(filename, '(a, i0, a)') trim(simulation_name)//'DEIMBasisRESVelocity_', count+1,".vtu"                       
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

  end subroutine read_pod_basis_deimres

  subroutine read_pod_basis_differntmesh(POD_state, state)
    !! Read the pod basis from the set of vtu files.
    
    character(len=1024) :: simulation_name, filename
    logical :: adjoint_reduced
    integer :: dump_period, quadrature_degree,istate
    integer :: i,j,k,total_dumps,nfield,POD_num,ifield
    
    type(state_type), dimension(:,:,:), allocatable :: POD_state
    type(state_type), dimension(:) :: state
    type(vector_field) :: podVelocity, newpodVelocity
    type(scalar_field) :: podPressure, newpodPressure, podTemperature, newpodTemperature
    type(mesh_type) :: VelocityMesh, PressureMesh, TemperatureMesh
    type(vector_field), pointer :: v_field
    type(scalar_field), pointer :: s_field
    
    velo => extract_vector_field(state, "Velocity")
    call get_option('/simulation_name', simulation_name)
    
    if (have_option("/reduced_model/Non_intrusive").and.(.not.have_option("/reduced_model/execute_reduced_model"))) then
     simulation_name=trim(simulation_name)//'_POD'
    endif

    call get_option('/geometry/quadrature/degree', quadrature_degree)
    adjoint_reduced= have_option("/reduced_model/adjoint")
    call get_option(&
         "/reduced_model/pod_basis_formation/pod_basis_count", POD_num) 
    nfield = vector_field_count( state(1) )+scalar_field_count( state(1) )
    allocate(pod_state(POD_num,nfield,size(state)))
    ifield = 0
    flag=1
    
    do k =1, size(state)
       ! Vector field
       !-------------
       nfield = vector_field_count( state(1) )
       do j = 1, nfield
          v_field => state(k)%vector_fields(j)%ptr              
          if(have_option(trim(v_field%option_path) // "/prognostic")) then
             ifield = ifield + 1
             VelocityMesh=extract_mesh(state(k),trim(v_field%mesh%name))
             do i = 1,POD_num
                if(have_option("/reduced_model/adjoint").and. .not.have_option("/reduced_model/execute_reduced_model")) then
                   write(filename, '(a, i0, a)') trim(simulation_name)//"_POD"//"Basis"//trim(v_field%name)//"_",i,".vtu" 
                else
                   write(filename, '(a, i0, a)') trim(simulation_name)//"Basis"//trim(v_field%name)//"_",i,".vtu" 
                endif
               ! print*,trim(filename)
                call vtk_read_state(filename, pod_state(i,ifield,k),quadrature_degree)
                 
                !! Note that we might need code in here to clean out unneeded fields.
                 
                 PODVelocity=extract_vector_field(pod_state(i,ifield,k),trim(v_field%name))
                 
                 call allocate(newpodVelocity, podVelocity%dim, VelocityMesh, trim(v_field%name))
                 !podVelocity%mesh%name = newpodVelocity%mesh%name
                 if(newpodVelocity%mesh%periodic) then
                    podVelocity%mesh%periodic = .true.
                 endif
                 call remap_field(from_field=podVelocity, to_field=newpodVelocity)
                 call insert(POD_state(i,ifield,k), newpodVelocity, trim(v_field%name))
                 call deallocate(newpodVelocity)
             enddo
          endif
       end do !j = 1, size(state(i)%vetor_fields)
       ! scaler field
       !-------------
       
       nfield = scalar_field_count( state(1) )
       do j = 1, nfield
          s_field => state(k)%scalar_fields(j)%ptr              
          if(have_option(trim(s_field%option_path) // "/prognostic")) then
          ifield = ifield + 1
!             call nullify(pod_state(:,k))
             PressureMesh=extract_mesh(state(k),trim(s_field%mesh%name))
             do i = 1,POD_num
                
                if(have_option("/reduced_model/adjoint").and. .not.have_option("/reduced_model/execute_reduced_model")) then
                   write(filename, '(a, i0, a)') trim(simulation_name)//"_POD"//"Basis"//trim(s_field%name)//"_",i,".vtu" 
                else
                   write(filename, '(a, i0, a)') trim(simulation_name)//"Basis"//trim(s_field%name)//"_",i,".vtu" 
                endif

                call vtk_read_state(filename, pod_state(i,ifield,k),quadrature_degree)

                !! Note that we might need code in here to clean out unneeded fields.
                PODPressure=extract_scalar_field(POD_state(i,ifield,k),trim(s_field%name))
                call allocate(newpodPressure, PressureMesh, trim(s_field%name))
                if(newpodPressure%mesh%periodic) then
                   PODPressure%mesh%periodic = .true.
                endif
                call remap_field(from_field=podPressure, to_field=newpodPressure)
                call insert(POD_state(i,ifield,k), newpodPressure, trim(s_field%name))
                call deallocate(newpodPressure)

             enddo
          endif
       end do !j = 1, size(state(i)%scalar_fields)
    end do!k =1, size(state)
  end subroutine read_pod_basis_differntmesh



end module reduced_model_runtime

