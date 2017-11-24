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

module coriolis_module

  use spud
  use vector_tools
  use fldebug
  use global_parameters, only : current_time, PYTHON_FUNC_LEN, FIELD_NAME_LEN
  use embed_python
  use fields
  use parallel_tools, only: abort_if_in_parallel_region
  use state_module

  implicit none
  
  !! Coriolis parameters:
  !! these are stored as global module variables and set directly
  !! from the options tree, this is to avoid having to read the options
  !! in performance critical routines in this module. They are private
  !! and should not be accessed directly
  real, save :: f0, dOmega_dt
  real, dimension(:), pointer, save :: coriolis_beta
  real, save :: latitude0, R_earth
  ! coriolis_option has to have one of the following values:
  integer, parameter :: NO_CORIOLIS=0, F_PLANE=1, BETA_PLANE=2, &
     SINE_OF_LATITUDE=3, CORIOLIS_ON_SPHERE=4, PYTHON_F_PLANE=5, CORIOLIS_WITH_AXIS=6, CORIOLIS_WITH_AXIS_UNSTEADY=7, CORIOLIS_FROM_FIELD=8, NOT_INITIALISED=-1
  integer, save :: coriolis_option=NOT_INITIALISED
  
  character(len = PYTHON_FUNC_LEN), save :: coriolis_python_func
  logical, save :: python_coriolis_initialised = .false.
  real, save :: python_coriolis_time
  real, dimension(3), save :: point_on_axis = [0.0,0.0,0.0]
  real, dimension(3), parameter :: axis = [0.0,0.0,1.0]
  
  ! legacy thing for use in funome():
  integer, save :: coriolis_dim=3

  type(scalar_field), pointer, save :: coriolis_field
     
  private
  
  public :: coriolis, funome, set_coriolis_parameters, centrifugal_force,&
     coriolis_module_check_options, coriolis_xyz

  contains

  function coriolis(X, ele)
    type(vector_field) :: X
    integer :: ele
    real, dimension(ele_ngi(X,ele)):: coriolis

    if (coriolis_option == CORIOLIS_FROM_FIELD) then
       coriolis = ele_val_at_quad(coriolis_field, ele)
    else
       coriolis = coriolis_xyz(ele_val_at_quad(X, ele))
    end if

  end function coriolis

  function coriolis_xyz(xyz) result(coriolis)
  !!< Returns the coriolis magnitude f so that the coriolis force is given
  !!< by: f k x u
  real, dimension(:,:):: xyz
  real, dimension(size(xyz,2)):: coriolis
  
    if (coriolis_option==NOT_INITIALISED) call set_coriolis_parameters
    
    select case (coriolis_option)
    case (NO_CORIOLIS)
      coriolis=0.0
    case (F_PLANE)
      coriolis=f0
    case (BETA_PLANE)
      coriolis=f0+matmul(coriolis_beta, xyz)
    case (SINE_OF_LATITUDE)
      coriolis=f0*sin(xyz(2,:)/R_earth+latitude0)
    case (CORIOLIS_ON_SPHERE)
      ! at the moment the same as f-plane:
      coriolis=f0
    case (PYTHON_F_PLANE)
      call update_f_plane_coriolis()
      coriolis = f0
    case (CORIOLIS_WITH_AXIS)
      coriolis=f0
   case (CORIOLIS_WITH_AXIS_UNSTEADY)
      call update_unsteady_coriolis()
      coriolis=f0
    case default
      ewrite(-1,*) "coriolis_option:", coriolis_option
      FLAbort("Unknown coriolis option")
    end select
  
  end function coriolis_xyz

  subroutine update_unsteady_coriolis()
    if (current_time/=python_coriolis_time) then
       f0=f0+2.0*dOmega_dt*(current_time-python_coriolis_time)
       python_coriolis_time=current_time
    end if
  end subroutine update_unsteady_coriolis



  subroutine update_f_plane_coriolis()
    !!< Update python set f-plane Coriolis variables

    if(.not. do_update_f_plane_coriolis()) return

    ewrite(2, *) "Updating f-plane Coriolis from python"
    call real_from_python(coriolis_python_func, current_time, f0)

    python_coriolis_initialised = .true.
    python_coriolis_time = current_time
   
  contains

    function do_update_f_plane_coriolis() result(update)
      logical :: update

      if(.not. python_coriolis_initialised) then
        update = .true.
      else
        update = python_coriolis_time /= current_time
      end if

    end function do_update_f_plane_coriolis

  end subroutine update_f_plane_coriolis
  
  function funome(xd,yd,zd)
  !!< This is a legacy wrapper around coriolis()
  !!< it returns omega, not f, so the coriolis force is 2*omega*(k x u)
    real :: funome
    real, intent(in):: xd,yd,zd
    
    real, dimension( 1:coriolis_dim, 1 ):: xyz
    real, dimension(1) :: omega
    
    xyz(1,1)=xd
    if (size(xyz,1)>1) then
      xyz(2,1)=yd
    end if
    if (size(xyz,1)>2) then
      xyz(3,1)=zd
    end if
    
    omega=coriolis_xyz(xyz)/2.0
    funome=omega(1)

  end function funome

  function centrifugal_force(X,u,ele) result (force)
    type(vector_field), intent(in) :: X,u
    integer, intent(in) :: ele
    real, dimension(mesh_dim(X),ele_ngi(u,ele)) :: force

    real, dimension(ele_ngi(u,ele)) :: f
    integer :: gi, dim
    real, dimension(mesh_dim(X),ele_ngi(X,ele)) :: X_quad
    real, dimension(mesh_dim(X)) :: x_cross

    !!! subroutine calculates the centrigual force source term at velocity 
    !!! quadrature points.

    !!! F=-\rho \vec{\Omega}\times(\vec{\Omega}\times\vec{x})  
    !!!  = -\rho (\Omega^2 \vec{x} - \Vec{\Omega}\cdot \vec{x} \vec{Omega}) 

    dim=mesh_dim(X)
    X_quad=ele_val_at_quad(X,ele)

    f = coriolis(X, ele)

    do gi=1,ele_ngi(u,ele)
       if(dim==2) then
          x_cross=[-X_quad(2,gi)+point_on_axis(2),X_quad(1,gi)-point_on_axis(1)]
       else
          x_cross=cross_product(axis,X_quad(:,gi))
       end if

       force(:,gi)=(f(gi)/2.0)**2*(X_quad(:,gi)-point_on_axis(1:dim)-dot_product(X_quad(:,gi)-point_on_axis(1:dim),axis(1:dim))*axis(1:dim))&
            +dOmega_dt*x_cross
    end do

  end function centrifugal_force

  subroutine centrifugal_potential(sfield,density,X)
    type(scalar_field), intent(inout) :: sfield
    type(scalar_field), intent(in) :: density
    type(vector_field), intent(in) :: X

    type(vector_field) :: X_on_s
    type(scalar_field) :: density_on_s

    integer :: node, dim

    dim=mesh_dim(X)

    call allocate(X_on_s,dim,sfield%mesh,"Coordinates")
    call allocate(density_on_s,sfield%mesh,"Coordinates")
    call remap_field(X,X_on_s)
    call remap_field(density,density_on_s)


    do node=1,node_count(sfield)

       call set(sfield,node,&
          (f0/2.0)**2*sum(cross_product(axis(1:dim),node_val(X_on_s,node)-point_on_axis(1:dim))**2)/2.0)
    end do
       
    call scale(sfield,density_on_s)
    
    call deallocate(X_on_s)
    call deallocate(density_on_s)


  end subroutine centrifugal_potential

  
  subroutine set_coriolis_parameters(state)

    type(state_type), intent(in), optional :: state
    
    real:: omega
    integer, dimension(2) :: dim

    character(len = FIELD_NAME_LEN) :: coriolis_field_name=""
    if (coriolis_option==CORIOLIS_FROM_FIELD) then
       call get_option("/physical_parameters/coriolis/from_field/source_name",&
            coriolis_field_name)
       coriolis_field => extract_scalar_field(state, coriolis_field_name) 
    end if
    if (coriolis_option/=NOT_INITIALISED) return

    call abort_if_in_parallel_region

    ewrite(1, *) "Initialising Coriolis"
    
    call get_option("/geometry/dimension", coriolis_dim)
    
    if (have_option("/physical_parameters/coriolis/f_plane")) then
         
       ewrite(2, *) "Coriolis type: f-plane"
      
       call get_option("/physical_parameters/coriolis/f_plane/f", f0)
       coriolis_option=F_PLANE
       
    else if (have_option("/physical_parameters/coriolis/beta_plane"))&
         & then 
         
       ewrite(2, *) "Coriolis type: beta-plane"
           
       call get_option("/physical_parameters/coriolis/beta_plane/f_0", f0)
       
       allocate( coriolis_beta(1:coriolis_dim) )
       call get_option("/physical_parameters/coriolis/&
            &beta_plane/beta", coriolis_beta)
       coriolis_option=BETA_PLANE
            
    else if (have_option("/physical_parameters/coriolis/sine_of_latitude")) then
         
       ewrite(2, *) "Coriolis type: Sine of latitude"
      
       call get_option("/physical_parameters/coriolis/sine_of_latitude/omega", omega)
       f0=2*omega
       
       call get_option("/physical_parameters/coriolis/sine_of_latitude/R_earth", R_earth)
       call get_option("/physical_parameters/coriolis/sine_of_latitude/latitude_0", latitude0)
       
       coriolis_option=SINE_OF_LATITUDE
       
    else if (have_option("/physical_parameters/coriolis/on_sphere"))&
         & then 
         
       ewrite(2, *) "Coriolis type: On sphere"
           
       call get_option("/physical_parameters/coriolis/on_sphere/&
            &omega", omega)
       f0=2*omega
       coriolis_option=CORIOLIS_ON_SPHERE
       
    else if(have_option("/physical_parameters/coriolis/python_f_plane")) then
    
      ewrite(2, *) "Coriolis type: f-plane (python)"
      
      call get_option("/physical_parameters/coriolis/python_f_plane", coriolis_python_func)
      coriolis_option = PYTHON_F_PLANE
    else if (have_option("/physical_parameters/coriolis/specified_axis")) then
       call get_option("/physical_parameters/coriolis/specified_axis/&
            &rotational_velocity", omega)  
       f0=2.0*omega
       point_on_axis=0.0
       dim=option_shape("/physical_parameters/coriolis/specified_axis/&
            &point_on_axis")
       call get_option("/physical_parameters/coriolis/specified_axis/&
            &point_on_axis", point_on_axis(1:dim(1)))  

       if (have_option("/physical_parameters/coriolis/specified_axis/&
            &rate_of_change")) then
          call get_option("/physical_parameters/coriolis/specified_axis/&
            &rate_of_change",dOmega_dt)
          coriolis_option=CORIOLIS_WITH_AXIS_UNSTEADY
       else
          dOmega_dt=0.0
          coriolis_option=CORIOLIS_WITH_AXIS
       end if
    else if (have_option("/physical_parameters/coriolis/from_field")) then
       call get_option("/physical_parameters/coriolis/from_field/source_name",&
            coriolis_field_name)
       coriolis_field => extract_scalar_field(state, coriolis_field_name)
       dOmega_dt=0.0
       coriolis_option=CORIOLIS_FROM_FIELD
    else 
    
      ewrite(2, *) "Coriolis type: None"
    
       coriolis_option=NO_CORIOLIS
       
    end if

  end subroutine set_coriolis_parameters
  
  subroutine coriolis_module_check_options()
    
    if (have_option("/geometry/spherical_earth") .and. &
      have_option("/physical_parameters/coriolis") .and. &
      .not. have_option("/physical_parameters/coriolis/on_sphere")) then
        
        ewrite(-1,*) "With /geometry/spherical_earth you need /physical_parameters/coriolis/on_sphere"
        FLExit("Fiddle with your FLML and try again...")
        
    end if
        
    if (have_option("/physical_parameters/coriolis/on_sphere") .and. &
        .not. have_option("/geometry/spherical_earth")) then
        
        ewrite(-1,*) "With /physical_parameters/coriolis/on_sphere you need /geometry/spherical_earth."
        FLExit("Fiddle with your FLML and try again...")
        
    end if
    
  end subroutine coriolis_module_check_options

end module coriolis_module
