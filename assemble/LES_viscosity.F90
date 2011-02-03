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

module les_viscosity_module
  !!< This module computes a viscosity term to implement LES
  use state_module
  use fields
  use spud
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use smoothing_module
  use vector_tools
  implicit none

  private

  public les_viscosity_strength, wale_viscosity_strength
  public leonard_tensor, les_strain_rate

contains

  subroutine leonard_tensor(state, u, positions, tnu, leonard, alpha, path)

    type(state_type), intent(inout)           :: state
    ! Unfiltered velocity
    type(vector_field), pointer, intent(inout)   :: u
    type(vector_field), intent(in)            :: positions
    ! Filtered velocity
    type(vector_field), pointer, intent(inout):: tnu
    ! Leonard tensor field
    type(tensor_field), pointer, intent(inout):: leonard
    ! Scale factor: test filter/mesh size
    real, intent(in)                          :: alpha
    character(len=OPTION_PATH_LEN), intent(in):: path
    ! Local quantities
    type(tensor_field), pointer               :: ui_uj, tui_tuj
    character(len=OPTION_PATH_LEN)            :: lpath
    integer                                   :: i, stat
    real, dimension(:), allocatable           :: u_loc
    real, dimension(:,:), allocatable         :: t_loc

    ! Path is to level above solver options
    lpath = (trim(path)//"/dynamic_les")
    ewrite(2,*) "path: ", trim(lpath)

    u => extract_vector_field(state, "Velocity", stat)

    do i = 1, positions%dim
      ewrite_minmax(u%val(i,:))
      ewrite_minmax(tnu%val(i,:))
    end do

    call anisotropic_smooth_vector(u, positions, tnu, alpha, lpath)

    do i = 1, positions%dim
      ewrite_minmax(u%val(i,:))
      ewrite_minmax(tnu%val(i,:))
    end do

    ! Velocity products (ui*uj)
    allocate(ui_uj); allocate(tui_tuj)
    call allocate(ui_uj, u%mesh, "NonlinearVelocityProduct")
    call allocate(tui_tuj, u%mesh, "TestNonlinearVelocityProduct")
    call zero(ui_uj); call zero(tui_tuj)

    ! Other local variables
    allocate(u_loc(u%dim)); allocate(t_loc(u%dim, u%dim))
    u_loc=0.0; t_loc=0.0

    ! Get cross products of velocities
    do i=1, node_count(u)
      ! Mesh filter ^r
      u_loc = node_val(u,i)
      call outer_product(u_loc, u_loc, t_loc)
      call set( ui_uj, i, t_loc )
      ! Test filter ^t
      u_loc = node_val(tnu,i)
      call outer_product(u_loc, u_loc, t_loc)
      ! Calculate (test-filtered velocity) products: (ui^rt*uj^rt)
      call set( tui_tuj, i, t_loc )
    end do

    ! Calculate test-filtered (velocity products): (ui^r*uj^r)^t
    call anisotropic_smooth_tensor(ui_uj, positions, leonard, alpha, lpath)

    ! Leonard tensor field
    call addto( leonard, tui_tuj, -1.0 )

    ! Deallocates
    deallocate(u_loc, t_loc)
    call deallocate(ui_uj)
    call deallocate(tui_tuj)
    deallocate(ui_uj); deallocate(tui_tuj)

  end subroutine leonard_tensor

  function les_strain_rate(du_t, nu)
    !! Computes the strain rate
    !! derivative of velocity shape function (nloc x ngi x dim)
    real, dimension(:,:,:), intent(in):: du_t
    !! nonlinear velocity (dim x nloc)
    real, dimension(:,:), intent(in):: nu
    real, dimension( size(du_t,3),size(du_t,3),size(du_t,2) ):: les_strain_rate
    real, dimension(size(du_t,3),size(du_t,3)):: s
    integer dim, ngi, gi

    ngi=size(du_t,2)
    dim=size(du_t,3)

    do gi=1, ngi

       s=0.5*matmul( nu, du_t(:,gi,:) )
       les_strain_rate(:,:,gi)=s+transpose(s)

    end do

  end function les_strain_rate

  function les_viscosity_strength(du_t, relu)
    !! Computes the strain rate modulus for the LES model 
    !! derivative of velocity shape function (nloc x ngi x dim)
    real, dimension(:,:,:), intent(in):: du_t
    !! relative velocity (nonl. vel.- grid vel.) (dim x nloc)
    real, dimension(:,:), intent(in):: relu

    real, dimension( size(du_t,2) ):: les_viscosity_strength

    real, dimension(size(du_t,3),size(du_t,3)):: s
    real vis
    integer dim, ngi, gi

    ngi=size(du_t,2)
    dim=size(du_t,3)

    do gi=1, ngi

       s=0.5*matmul( relu, du_t(:,gi,:) )
       s=s+transpose(s)
       ! Calculate modulus of strain rate
       vis=sqrt( 2*sum( s**2 ) )

       les_viscosity_strength(gi)=vis

    end do

  end function les_viscosity_strength

  function wale_viscosity_strength(du_t, relu)
    !! Computes the traceless symmetric part of the square of
    !! the resolved velocity gradient tensor for the LES model
    !! See a WALE paper for more (G_{ij})
    !! derivative of velocity shape function (nloc x ngi x dim)
    real, dimension(:,:,:), intent(in):: du_t
    !! relative velocity (nonl. vel.- grid vel.) (dim x nloc)
    real, dimension(:,:), intent(in):: relu

    real, dimension( size(du_t,2) ):: wale_viscosity_strength

    real, dimension(size(du_t,3),size(du_t,3)):: s, g
    real vis
    integer dim, ngi, gi, i

    ngi=size(du_t,2)
    dim=size(du_t,3)

    do gi=1, ngi

       s=matmul( relu, du_t(:,gi,:) )
       g=0.5*matmul(s,s)
       g=g+transpose(g)
       forall(i=1:dim) g(i,i)=0.
       
       vis=sqrt( 2*sum( g**2 ) )

       wale_viscosity_strength(gi)=vis

    end do

  end function wale_viscosity_strength

end module les_viscosity_module
