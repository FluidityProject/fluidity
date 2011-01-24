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
  use fields
  implicit none

  private

  public les_viscosity_strength, wale_viscosity_strength
!         &dynamic_les_init, leonard_tensor

contains

  !subroutine dynamic_les_init(nu, du_t)

    !u_mesh = nu%mesh
    ! Create a test-filter mesh from the velocity mesh
    ! that has the same nodes but different connectivity
    !test_mesh = make_dynamic_les_mesh(u_mesh, test_mesh)

    ! Copy the nonlinear velocity field on to the coarser test mesh
    !test_nu = wrap_vector_field(test_mesh, nu%val, "TestNonlinearVelocity")

    !do ele = 1, ele_count(test_mesh)
    !  ltens = leonard_tensor(ele, nu, test_mesh)

    !end do

  !end subroutine dynamic_les_init



!  subroutine smooth_velocity(nu, nu_test)

!    path = trim(complete_field_path(trim(nu%option_path)))
    !alpha = MUST BE A CONSTANT! Change to variable? See papers
    ! Get scalar components of velocity
!    do dim=1,nu%dim
!      scalar_component => extract_scalar_field(nu, dim, stat)
!      call smooth_scalar(scalar_component, positions, smooth_scalar_component, alpha, path=path)               
!      ! set vector field from scalar
!      call set(test_nu, dim, smooth_scalar_component)
!    end do

!  end subroutine dynamic_les_fields



!  function leonard_tensor(ele, nu, test_nu) result(leonard_tensor)

    ! Smoothed velocity field
!    call allocate(test_nu, nu%mesh, "Smoothed" // trim(nu%name))
!    call zero(test_nu)

    !path = trim(complete_field_path(trim(field_in%option_path))) // "/adaptivity_options/preprocessing/helmholtz_smoother"
!    call smooth_scalar(nu, positions, test_nu, alpha, path=path)

    ! Velocity products loop
!    u_mesh = nu%mesh
!    neigh => ele_neigh(u_mesh, ele)
!    nloc = size(ele_nodes(u_mesh, ele))
!    do n = 1, size(neigh)
!      ele2 = neigh(n) 
      ! find global face no. of face associated with neighbouring element
!      face2 = ele_face(u_mesh, ele2, ele)
      ! By convention, node opposite a face has same local number
!      local_face2 = local_face_number(u_mesh, face2)
      ! Find the global node no. of the opposite node
!      ele2_nodes => ele_nodes(u_mesh, ele2)
      ! ele vertex n is matched to global node with same number as local_face2
!      node2 = ele2_nodes(local_face2)

      ! Get velocity at foreign node
!    end do

!    nu_quad = ele_val_at_quad(nu, ele)
!    tnu_quad = ele_val_at_quad(test_nu, ele)
!    ngi = size(nu_quad,2)
!    dim = size(nu_quad,3)

    ! Tensor construction loop
!    do gi=1, ngi  
      ! product of test-filtered velocity with itself
!      t1 = outer_product( tnu_quad(:,gi), tnu_quad(:,gi))
      ! test-filtered product of velocity with itself
!      t2 = outer_product( nu_quad(:,gi), nu_quad(:,gi))

      ! At each Gauss point, I need to know the corresponding
      ! point in the test element. Then I can do outer_product on test element.


!      leonard_tensor(:,:,gi) = M
!    end do

!  end function leonard_tensor



!  function strain_tensor(ele, nu, du_t, test_mesh) result(strain_tensor)

    ! Compute strain tensor on velocity mesh
!    strain_mesh = les_viscosity_strength(du_t, ele_val(nu, ele))

    ! Compute strain tensor on test-filter mesh
!    shape_test = ele_shape(test_mesh, ele)
!    allocate( dshape_test (ele_loc(test_mesh,ele), ele_ngi(test_mesh,ele), nu%dim) )
!    allocate( detwei_test (ele_ngi(test_mesh,ele) ) )
    ! Get transformed velocity shape function gradients in physical space
    ! on test element using the positions x
!    call transform_to_physical( positions, ele, shape_test, dshape=dshape_test, detwei=detwei_test )
    ! Compute strain tensor using test-filtered gradients (coarse mesh)
!    strain_test = les_viscosity_strength(dshape_test, ele_val(nu, ele))
!    deallocate(dshape_test, detwei_test)

!  end function strain_tensor



!  function make_dynamic_les_mesh(u_mesh, name) result(test_mesh)

!    call allocate(test_mesh, nodes=u_mesh%nodes, elements=u_mesh%elements, &
!                  &shape=u_mesh%shape, name=dynamic_les_mesh)

    ! Mesh continuity makes no sense for this mesh as elements overlap!
!    mesh%continuity=-666
    
!    assert(has_faces(u_mesh))
!    u_mesh%ndglno = -1
!    do ele = 1, element_count(u_mesh)
       ! create the ndglno by looping around the neighbouring elements
!       neigh => ele_neigh(u_mesh, ele)
!       nloc = size(ele_nodes(u_mesh, ele))
!       do n = 1, size(neigh)
!          ele2 = neigh(n)
          ! find global face no. of face associated with neighbouring element
!          face2 = ele_face(u_mesh, ele2, ele)
          ! By convention, node opposite a face has same local number
!          local_face2 = local_face_number(u_mesh, face2)
          ! Find the global node no. of the opposite node
!          ele2_nodes => ele_nodes(u_mesh, ele2)
          ! ele vertex n is matched to global node with same number as local_face2
!          test_mesh%ndglno(u_mesh%shape%loc*(ele-1)*nloc+n) = ele2_nodes(local_face2)
!      end do

!    end do
!    assert(all(mesh%ndglno > 0))
!    call addref(test_mesh)
!    ewrite(1,*) 'exiting make_dynamic_les_mesh'

!  end function make_dynamic_les_mesh



  function les_viscosity_strength(du_t, relu)
    !! Computes the strain rate for the LES model 
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
