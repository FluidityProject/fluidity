!    Copyright (C) 2010 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Dr Chris Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    c.pain@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
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

module Profiler
  use global_parameters, only : real_8
  use fields
  implicit none
  
  private
  
  public profiler_tic, profiler_toc, profiler_zero
  
  interface profiler_tic
    module procedure profiler_tic_scalar, profiler_tic_vector, &
      profiler_tic_tensor, profiler_tic_key
  end interface profiler_tic
  
  interface profiler_toc
    module procedure profiler_toc_scalar, profiler_toc_vector, &
      profiler_toc_tensor, profiler_toc_key
  end interface profiler_toc
    
  interface profiler_get
    module procedure profiler_get_scalar, profiler_get_vector, &
      profiler_get_tensor, profiler_get_key
  end interface profiler_get

contains
  
  real(kind = real_8) function profiler_get_scalar(field, action)
    !!< Starts profiling a certain 'action', .e.g.: assembly, solve, interpolation
    !!< for field and stores this under a unique key
    type(scalar_field), intent(in):: field
    character(len=*), intent(in):: action
    
    profiler_get_scalar=profiler_get_key(create_profile_key(field%option_path, action))
    
  end function profiler_get_scalar
  
  real(kind = real_8) function profiler_get_vector(field, action)
    !!< Starts profiling a certain 'action', .e.g.: assembly, solve, interpolation
    !!< for field and stores this under a unique key
    type(vector_field), intent(in):: field
    character(len=*), intent(in):: action
    
    profiler_get_vector=profiler_get_key(create_profile_key(field%option_path, action))
    
  end function profiler_get_vector

  real(kind = real_8) function profiler_get_tensor(field, action)
    !!< Starts profiling a certain 'action', .e.g.: assembly, solve, interpolation
    !!< for field and stores this under a unique key
    type(tensor_field), intent(in):: field
    character(len=*), intent(in):: action
    
    profiler_get_tensor=profiler_get_key(create_profile_key(field%option_path, action))
    
  end function profiler_get_tensor
  
  real(kind = real_8) function profiler_get_key(key)
    character(len=*), intent(in)::key
    external cprofiler_get
    call cprofiler_get(key, len_trim(key), profiler_get_key)
  end function profiler_get_key
  
  subroutine profiler_tic_scalar(field, action)
    !!< Starts profiling a certain 'action', .e.g.: assembly, solve, interpolation
    !!< for field and stores this under a unique key
    type(scalar_field), intent(in):: field
    character(len=*), intent(in):: action
    
    call profiler_tic_key(create_profile_key(field%option_path, action))
    
  end subroutine profiler_tic_scalar

  subroutine profiler_tic_vector(field, action)
    !!< Starts profiling a certain 'action', .e.g.: assembly, solve, interpolation
    !!< for field and stores this under a unique key
    type(vector_field), intent(in):: field
    character(len=*), intent(in):: action
    
    call profiler_tic_key(create_profile_key(field%option_path, action))
    
  end subroutine profiler_tic_vector

  subroutine profiler_tic_tensor(field, action)
    !!< Starts profiling a certain 'action', .e.g.: assembly, solve, interpolation
    !!< for field and stores this under a unique key
    type(tensor_field), intent(in):: field
    character(len=*), intent(in):: action
    
    call profiler_tic_key(create_profile_key(field%option_path, action))
    
  end subroutine profiler_tic_tensor
    
  subroutine profiler_toc_scalar(field, action)
    !!< Starts profiling a certain 'action', .e.g.: assembly, solve, interpolation
    !!< for field and stores this under a unique key
    type(scalar_field), intent(in):: field
    character(len=*), intent(in):: action
    
    call profiler_toc_key(create_profile_key(field%option_path, action))
    
  end subroutine profiler_toc_scalar

  subroutine profiler_toc_vector(field, action)
    !!< Starts profiling a certain 'action', .e.g.: assembly, solve, interpolation
    !!< for field and stores this under a unique key
    type(vector_field), intent(in):: field
    character(len=*), intent(in):: action
    
    call profiler_toc_key(create_profile_key(field%option_path, action))
    
  end subroutine profiler_toc_vector

  subroutine profiler_toc_tensor(field, action)
    !!< Starts profiling a certain 'action', .e.g.: assembly, solve, interpolation
    !!< for field and stores this under a unique key
    type(tensor_field), intent(in):: field
    character(len=*), intent(in):: action
    
    call profiler_toc_key(create_profile_key(field%option_path, action))
    
  end subroutine profiler_toc_tensor
    
  function create_profile_key(option_path, action)
    character(len=*), intent(in):: option_path, action
    character(len=len_trim(option_path)+2+len_trim(action)):: create_profile_key
    
    create_profile_key=trim(option_path)//'::'//trim(action)
  
  end function create_profile_key
  
  subroutine profiler_tic_key(key)
    character(len=*), intent(in)::key
    external cprofiler_tic
    call cprofiler_tic(key, len_trim(key))
  end subroutine profiler_tic_key

  subroutine profiler_toc_key(key)
    character(len=*), intent(in)::key
    external cprofiler_toc
    call cprofiler_toc(key, len_trim(key))
  end subroutine profiler_toc_key

  subroutine profiler_zero()
    external cprofiler_zero
    call cprofiler_zero()
  end subroutine profiler_zero

end module Profiler
