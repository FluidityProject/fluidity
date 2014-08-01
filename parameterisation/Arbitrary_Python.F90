!    Copyright (C) 2008 Imperial College London and others.
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

module arbitrary_function
  !! This module implements evaluating an arbitrary python function at various points
  use spud
  use state_module
  use fields
  use sparse_tools
  use fefields
  use fetools
  use boundary_conditions
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN,&
       PYTHON_FUNC_LEN, CURRENT_DEBUG_LEVEL
  use solvers
  use python_state
  use sparsity_patterns_meshes
  use equation_of_state
  use quadrature

  implicit none 

  integer, parameter :: SCALAR_FIELD_TYPE=1,VECTOR_FIELD_TYPE=2,&
       TENSOR_FIELD_TYPE=3

  type function_data
     character (len=PYTHON_FUNC_LEN) :: function_string
     character (len=OPTION_PATH_LEN) :: option_path
     integer, dimension(:), allocatable :: argument_type
     type(scalar_field_pointer), dimension(:), allocatable :: scalars
     type(vector_field_pointer), dimension(:), allocatable :: vectors
     type(tensor_field_pointer), dimension(:), allocatable :: tensors
  end type function_data

  

  private
  public function_data, initialize_function_data,calculate_function_at_quads,&
       calculate_function_at_nodes, finalize_function_data
contains

  subroutine initialize_function_data(state,option_path,fun_data)
    type(state_type), dimension(:) :: state
    character (len=*), intent(in) :: option_path
    type(function_data) :: fun_data

    integer no_scalars,no_vectors,no_tensors
    integer :: argument,scalar_it,vector_it,tensor_it
    character (len=FIELD_NAME_LEN) :: field_name, field_type, phase
    logical :: is_multiphase
    
    

    fun_data%option_path=option_path

    call get_option(trim(option_path)//"/function",fun_data%function_string)


    no_scalars=option_count(trim(option_path)//"/argument/type::scalar_field")
    no_vectors=option_count(trim(option_path)//"/argument/type::vector_field")
    no_tensors=option_count(trim(option_path)//"/argument/type::tensor_field")

    allocate(fun_data%scalars(no_scalars),&
             fun_data%vectors(no_vectors),&
             fun_data%tensors(no_tensors),&
             fun_data%argument_type(no_scalars+no_vectors+no_tensors))
             
    scalar_it=1
    vector_it=1
    tensor_it=1
    do argument=1,no_scalars+no_vectors+no_tensors
       call get_option(option_path//"/argument["//&
            int2str(argument-1)//"]/name",field_name)
       call get_option(option_path//"/argument["//&
            int2str(argument-1)//"]/type/name",field_type)
       is_multiphase=have_option(trim(option_path)//"/argument["//&
            int2str(argument-1)//"]/material_phase")
       if (is_multiphase)&
            call get_option(trim(option_path)//"/argument["//&
            int2str(argument-1)//"]/material_phase/name",phase)
       select case(trim(field_type))
       case("scalar_field")
          fun_data%argument_type(argument)=SCALAR_FIELD_TYPE
          if (is_multiphase) then
             field_name=trim(field_name(index(field_name,':',back=.true.)+1:len(field_name)))
             fun_data%scalars(scalar_it)%ptr=>&
               extract_scalar_field(state,trim(field_name),state_name=trim(phase))
          else
          fun_data%scalars(scalar_it)%ptr=>&
               extract_scalar_field(state,field_name)
          end if
          scalar_it=scalar_it+1
       case("vector_field")
          fun_data%argument_type(argument)=VECTOR_FIELD_TYPE
          if (is_multiphase) then
             field_name=trim(field_name(index(field_name,':',back=.true.)+1:len(field_name)))
             fun_data%vectors(vector_it)%ptr=>&
               extract_vector_field(state,field_name,state_name=trim(phase))
          else
             fun_data%vectors(vector_it)%ptr=>&
                  extract_vector_field(state,field_name)
          end if
          vector_it=vector_it+1
       case("tensor_field")
          fun_data%argument_type(argument)=TENSOR_FIELD_TYPE
          fun_data%tensors(tensor_it)%ptr=>&
               extract_tensor_field(state(1),field_name)
          tensor_it=tensor_it+1
       end select
    end do

  end subroutine initialize_function_data

  subroutine finalize_function_data(fun_data)
    type(function_data) :: fun_data

    deallocate(fun_data%argument_type)
    deallocate(fun_data%scalars)
    deallocate(fun_data%vectors)
    deallocate(fun_data%tensors)

  end subroutine finalize_function_data

  subroutine calculate_function_at_nodes(fun_data,fun_result)
 
    type(function_data) :: fun_data
    type(scalar_field) :: fun_result

    type(mesh_type)  :: mesh
    type(state_type) :: local_state

    type(scalar_field) :: scalar
    type(vector_field) :: vector
    type(tensor_field) :: tensor

    integer ::scalar_itr,vector_itr,tensor_itr,argument


    local_state%name="local_state"
    mesh=fun_result%mesh
    call insert(local_state,mesh,mesh%name)

    do scalar_itr=1,size(fun_data%scalars)
       call allocate(scalar,mesh,trim(fun_data%scalars(scalar_itr)%ptr%name)//int2str(scalar_itr))
       call remap_field(fun_data%scalars(scalar_itr)%ptr,scalar)
       call insert(local_state,scalar,scalar%name)
       call deallocate(scalar)
    end do
    do vector_itr=1,size(fun_data%vectors)
       call allocate(vector,fun_data%vectors(vector_itr)%ptr%dim,mesh,&
            trim(fun_data%vectors(vector_itr)%ptr%name)//int2str(vector_itr))
       call remap_field(fun_data%vectors(vector_itr)%ptr,vector)
       call insert(local_state,vector,vector%name)
       call deallocate(vector)
    end do
    do scalar_itr=1,size(fun_data%tensors)
       call allocate(tensor,mesh,dim=fun_data%tensors(tensor_itr)%ptr%dim,&
            name=trim(fun_data%tensors(tensor_itr)%ptr%name)//int2str(tensor_itr))
       call remap_field(fun_data%tensors(tensor_itr)%ptr,tensor)
       call insert(local_state,tensor,tensor%name)
       call deallocate(tensor)
    end do

    call allocate(scalar,mesh,'result')
    call insert(local_state,scalar,"result")
    call deallocate(scalar)

    !! set up python space

    call python_reset()
    call python_add_state(local_state)

    scalar_itr=1
    vector_itr=1
    tensor_itr=1

    call python_run_string("args=[]")

    do argument=1,size(fun_data%argument_type)
       select case(fun_data%argument_type(argument))
       case(SCALAR_FIELD_TYPE)
          call python_run_string("args.append(state.scalar_fields['"//&
               trim(fun_data%scalars(scalar_itr)%ptr%name)//int2str(scalar_itr)//"'].val)")
          scalar_itr=scalar_itr+1
       case(VECTOR_FIELD_TYPE)
          call python_run_string("args.append(state.vector_fields['"//&
               trim(fun_data%vectors(vector_itr)%ptr%name)//int2str(vector_itr)//"'].val)")
          vector_itr=vector_itr+1
       case(TENSOR_FIELD_TYPE)
          call python_run_string("args.append(state.tensor_fields['"//&
               trim(fun_data%tensors(tensor_itr)%ptr%name)//int2str(tensor_itr)//"'].val)")
          vector_itr=vector_itr+1
       end select

    end do
       
  call python_run_string(trim(fun_data%function_string))

  call python_run_string("state.scalar_fields['result'].val[:]=val(*args)")

  call python_reset()

  call set(fun_result,extract_scalar_field(local_state,"result"))
  
  call deallocate(local_state)
 
end subroutine calculate_function_at_nodes

function calculate_function_at_quads(fun_data,shape,ele) result(val_at_quads)
  type(function_data) :: fun_data
  type(element_type) :: shape
  integer :: ele
  real, dimension(shape%ngi) :: val_at_quads

  integer :: argument, scalar_itr,vector_itr,tensor_itr


  call python_reset()
  call python_run_string(trim(fun_data%function_string))
  
  scalar_itr=1
  vector_itr=1
  tensor_itr=1
  
  call python_run_string("args=[]")
  
  do argument=1,size(fun_data%argument_type)
     select case(fun_data%argument_type(argument))
     case(SCALAR_FIELD_TYPE)
        call python_add_array(ele_val_at_quad(fun_data%scalars(scalar_itr)%ptr,ele,shape),'y')
        CALL python_run_string("args.append(y)")
        scalar_itr=scalar_itr+1
     case(VECTOR_FIELD_TYPE)
        call python_run_string("args.append(state.vector_fields['"//&
             trim(fun_data%vectors(vector_itr)%ptr%name)//"'].val)")
        vector_itr=vector_itr+1
     case(TENSOR_FIELD_TYPE)
        call python_run_string("args.append(state.tensor_fields['"//&
             trim(fun_data%tensors(tensor_itr)%ptr%name)//"'].val)")
        vector_itr=vector_itr+1
     end select
     
  end do

  call python_add_array(val_at_quads,'val_at_quads')
  call python_run_string('val_at_quads[:]=val(*args)')
  call python_reset()
  
  end function calculate_function_at_quads




end module arbitrary_function
