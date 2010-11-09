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
module elements
  !!< This module provides derived types for finite elements and associated functions.
  use element_numbering
  use quadrature
  use FLDebug
  use polynomials
  use reference_counting
  implicit none

  type element_type
     !!< Type to encode shape and quadrature information for an element.
     integer :: dim !! 2d or 3d?
     integer :: loc !! Number of nodes.
     integer :: ngi !! Number of gauss points.
     integer :: degree !! Polynomial degree of element.
     !! Shape functions: n is for the primitive function, dn is for partial derivatives, dn_s is for partial derivatives on surfaces. 
     !! n is loc x ngi, dn is loc x ngi x dim
     !! dn_s is loc x ngi x face x dim 
     real, pointer :: n(:,:)=>null(), dn(:,:,:)=>null()
     real, pointer :: n_s(:,:,:)=>null(), dn_s(:,:,:,:)=>null()
     !! Polynomials defining shape functions and their derivatives.
     type(polynomial), dimension(:,:), pointer :: spoly=>null(), dspoly=>null()
     !! Link back to the node numbering used for this element.
     type(ele_numbering_type), pointer :: numbering=>null()
     !! Link back to the quadrature used for this element.
     type(quadrature_type) :: quadrature
     type(quadrature_type), pointer :: surface_quadrature=>null()
     !! Pointer to the superconvergence data for this element.
     type(superconvergence_type), pointer :: superconvergence=>null()
     !! Reference count to prevent memory leaks.
     type(refcount_type), pointer :: refcount=>null()
     !! Dummy name to satisfy reference counting
     character(len=0) :: name
  end type element_type

  type superconvergence_type
    !!< A structure to represent the superconvergent points of the element in question.
    !!< This is in this module because it has to be in element_type,
    !!< but Superconvergence.F90 depends on Elements.F90. So Elements.F90
    !!< cannot depend on Superconvergence.F90. (Fortran is a real pain.)
    !! Number of superconvergent points
    integer :: nsp 
    !! Locations of superconvergent points in local coordinates
    !! allocated to nsp x loc
    real, pointer :: l(:, :)
    !! Shape functions at each superconvergent point.
    !! loc x nsp
    real, pointer :: n(:, :)
    !! Derivatives of shape functions at each superconvergent point
    !! loc x nsp x ndim
    real, pointer :: dn(:, :, :)
  end type superconvergence_type

  interface allocate
     module procedure allocate_element, allocate_element_with_surface
  end interface

  interface deallocate
     module procedure deallocate_element
  end interface

  interface local_coords
     module procedure element_local_coords
  end interface

  interface local_coord_count
     module procedure element_local_coord_count
  end interface

  interface local_vertices
     module procedure element_local_vertices
  end interface

  interface boundary_numbering
     module procedure element_boundary_numbering
  end interface

  interface operator(==)
     module procedure element_equal
  end interface

  interface eval_shape
    module procedure eval_shape_node, eval_shape_all_nodes
  end interface

  interface eval_dshape
    module procedure eval_dshape_node, eval_dshape_all_nodes
  end interface

#include "Reference_count_interface_element_type.F90"

contains

  subroutine allocate_element(element, dim, loc, ngi, coords, type, stat)
    !!< Allocate memory for an element_type. 
    type(element_type), intent(inout) :: element
    !! Dim is the dimension of the element, loc is number of nodes, ngi is
    !! number of gauss points. 
    integer, intent(in) :: dim,loc,ngi    
    !! Number of local coordinates.
    integer, intent(in) :: coords
    !! Stat returns zero for success and nonzero otherwise.
    integer, intent(out), optional :: stat
    !! define element type
    integer, intent(in), optional :: type

    integer :: lstat, ltype

    if (present(type)) then
       ltype=type
    else
       ltype=ELEMENT_LAGRANGIAN
    end if

    select case(ltype)
    case(ELEMENT_LAGRANGIAN, ELEMENT_NONCONFORMING, ELEMENT_BUBBLE)

      allocate(element%n(loc,ngi),element%dn(loc,ngi,dim), &
          element%spoly(coords,loc), element%dspoly(coords,loc), stat=lstat)

    case(ELEMENT_CONTROLVOLUME_SURFACE)

      allocate(element%n(loc,ngi),element%dn(loc,ngi,dim-1), &
          stat=lstat)

      element%spoly=>null()
      element%dspoly=>null()

    case(ELEMENT_CONTROLVOLUMEBDY_SURFACE)

      allocate(element%n(loc,ngi),element%dn(loc,ngi,dim), &
          stat=lstat)

      element%spoly=>null()
      element%dspoly=>null()

    case default

      FLAbort("Attempt to select an illegal element type.")

    end select

    element%loc=loc
    element%ngi=ngi
    element%dim=dim

    nullify(element%refcount) ! Hack for gfortran component initialisation
    !                         bug.
    call addref(element)
    
    nullify(element%n_s)
    nullify(element%dn_s)

    if (present(stat)) then
       stat=lstat
    else if (lstat/=0) then
       FLAbort("Unable to allocate element.")
    end if

  end subroutine allocate_element

  subroutine allocate_element_with_surface(element, dim, loc,&
       ngi,faces, ngi_s, coords,surface_present,type, stat)
    !!< Allocate memory for an element_type. 
    type(element_type), intent(inout) :: element
    !! Dim is the dimension of the element, loc is number of nodes, ngi is
    !! number of gauss points. 
    integer, intent(in) :: dim,loc,ngi,faces,ngi_s    
    !! Number of local coordinates.
    integer, intent(in) :: coords
    logical, intent(in) :: surface_present
    !! Stat returns zero for success and nonzero otherwise.
    integer, intent(in), optional :: type
    integer, intent(out), optional :: stat

    integer :: lstat

    allocate(element%n(loc,ngi),element%dn(loc,ngi,dim), &
         element%n_s(loc,ngi_s,faces),element%dn_s(loc,ngi_s,faces,dim),&
         element%spoly(coords,loc), element%dspoly(coords,loc), stat=lstat)
    
    element%loc=loc
    element%ngi=ngi
    element%dim=dim

    if (present(stat)) then
       stat=lstat
    else if (lstat/=0) then
       FLAbort("Unable to allocate element.")
    end if
    
    nullify(element%refcount) ! Hack for gfortran component initialisation
    !                         bug.
    call addref(element)

  end subroutine allocate_element_with_surface

  subroutine deallocate_element(element, stat)
    type(element_type), intent(inout) :: element
    integer, intent(out), optional :: stat
    
    integer :: lstat, tstat
    integer :: i,j

    tstat = 0
    lstat = 0

    call decref(element)
    if (has_references(element)) then
       ! There are still references to this element so we don't deallocate.
       return
    end if

    call deallocate(element%quadrature)

    if(associated(element%spoly)) then
      do i=1,size(element%spoly,1)
        do j=1,size(element%spoly,2)
            call deallocate(element%spoly(i,j))
        end do
      end do
      deallocate(element%spoly, stat=tstat)
    end if
    lstat=max(lstat,tstat)

    if (associated(element%n_s)) deallocate(element%n_s,element%dn_s)

    if(associated(element%dspoly)) then
      do i=1,size(element%dspoly,1)
        do j=1,size(element%dspoly,2)
            call deallocate(element%dspoly(i,j))
        end do
      end do
      deallocate(element%dspoly, stat=tstat)
    end if
    lstat=max(lstat,tstat)

    deallocate(element%n,element%dn, stat=tstat)
    lstat=max(lstat,tstat)

    if (present(stat)) then
       stat=lstat
    else if (lstat/=0) then
       FLAbort("Unable to deallocate element.")
    end if

  end subroutine deallocate_element

  function element_local_coords(n, element) result (coords)
    !!< Work out the local coordinates of node n in element. This is just a
    !!< wrapper function which allows local_coords to be called on an element
    !!< instead of on an element numbering.
    integer, intent(in) :: n
    type(element_type), intent(in) :: element    
    real, dimension(size(element%numbering%number2count, 1)) :: coords
    
    coords=local_coords(n, element%numbering)

  end function element_local_coords
  
  function element_local_coord_count(element) result (n)
    !!< Return the number of local coordinates associated with element.
    integer :: n
    type(element_type), intent(in) :: element    

    n=size(element%numbering%number2count, 1)

  end function element_local_coord_count

  function element_local_vertices(element) result (vertices)
    !!< Given an element numbering, return the local node numbers of its
    !!< vertices. This is just a wrapper hich allows local_vertices to 
    !!< be called on an element instead of on an element numbering.
    type(element_type), intent(in) :: element
    integer, dimension(element%numbering%vertices) :: vertices
    
    vertices=local_vertices(element%numbering)
    
  end function element_local_vertices

  function element_boundary_numbering(element, boundary)
    !!< A wrapper function which allows boundary_numbering to be called on
    !!< an element instead of on an element_numbering.
    integer, intent(in) :: boundary
    type(element_type), intent(in) :: element
    integer, dimension(boundary_num_length(element%numbering, .false.)) ::&
         & element_boundary_numbering 
    
    element_boundary_numbering=boundary_numbering(element%numbering,&
         & boundary)

  end function element_boundary_numbering

  pure function element_equal(element1,element2)
    !!< Return true if the two elements are equivalent.
    logical :: element_equal
    type(element_type), intent(in) :: element1, element2
    
    element_equal = element1%dim==element2%dim &
         .and. element1%loc==element2%loc &
         .and. element1%ngi==element2%ngi &
         .and. element1%numbering==element2%numbering &
         .and. element1%quadrature==element2%quadrature
    
  end function element_equal

  subroutine extract_old_element(element, N, NLX, NLY, NLZ)
    !!< Extract the shape function values from an old element.
    type(element_type), intent(in) :: element
    real, dimension(element%loc, element%ngi), intent(out) :: N, NLX, NLY
    real, dimension(element%loc, element%ngi), intent(out), optional :: NLZ
    
    N=element%n
    NLX=element%dn(:,:,1)
    if (size(element%dn,3)>1) then
       NLY=element%dn(:,:,2)
    else
       NLY=0.0
    end if

    if (present(NLZ)) then
       if (size(element%dn,3)>2) then
          NLZ=element%dn(:,:,3)
       else
          NLZ=0.0
       end if
    end if


  end subroutine extract_old_element

  pure function eval_shape_node(shape, node,  l) result(eval_shape)
    ! Evaluate the shape function for node node local coordinates l
    real :: eval_shape
    type(element_type), intent(in) :: shape
    integer, intent(in) :: node
    real, dimension(size(shape%spoly,1)), intent(in) :: l

    integer :: i

    eval_shape=1.0
          
    do i=1,size(shape%spoly,1)
       
       ! Raw shape function
       eval_shape=eval_shape*eval(shape%spoly(i,node), l(i))
             
    end do

  end function eval_shape_node

  pure function eval_shape_all_nodes(shape, l) result(eval_shape)
    ! Evaluate the shape function for all locations at local coordinates l
    type(element_type), intent(in) :: shape
    real, dimension(size(shape%spoly,1)), intent(in) :: l
    real, dimension(shape%loc) :: eval_shape

    integer :: i,j

    eval_shape=1.0

    do j=1,shape%loc

      do i=1,size(shape%spoly,1)

        ! Raw shape function
        eval_shape(j)=eval_shape(j)*eval(shape%spoly(i,j), l(i))

      end do

    end do

  end function eval_shape_all_nodes

  pure function eval_dshape_node(shape, node,  l) result(eval_dshape)
    !!< Evaluate the derivatives of the shape function for location node at local
    !!< coordinates l 
    type(element_type), intent(in) :: shape
    integer, intent(in) :: node
    real, dimension(:), intent(in) :: l
    real, dimension(shape%dim) :: eval_dshape

    select case(shape%numbering%family)
       
    case (FAMILY_SIMPLEX)

       eval_dshape=eval_dshape_simplex(shape, node,  l)

    case (FAMILY_CUBE)

       eval_dshape=eval_dshape_cube(shape, node,  l)

    case default
       ! Invalid element family. Return a really big number to stuff things
       ! quickly. 

       eval_dshape=huge(0.0)

    end select
    
  end function eval_dshape_node

  function eval_dshape_all_nodes(shape, l) result(eval_dshape)
    type(element_type), intent(in) :: shape
    real, dimension(:), intent(in) :: l
    real, dimension(shape%loc, shape%dim) :: eval_dshape

    integer :: loc

    do loc=1,shape%loc
      eval_dshape(loc, :) = eval_dshape_node(shape, loc, l)
    end do
  end function eval_dshape_all_nodes

  function eval_dshape_transformed(shape, l, invJ) result(transformed_dshape)
    type(element_type), intent(in) :: shape
    real, dimension(:), intent(in) :: l
    real, dimension(shape%dim, shape%dim), intent(in) :: invJ
    real, dimension(shape%loc, shape%dim) :: transformed_dshape, untransformed_dshape

    integer :: loc

    do loc=1,shape%loc
      untransformed_dshape(loc, :) = eval_dshape_node(shape, loc, l)
      transformed_dshape(loc, :) = matmul(invJ, untransformed_dshape(loc, :))
    end do
  end function eval_dshape_transformed

  function eval_volume_dshape_at_face_quad(shape, local_face_number, invJ) result(output)
    ! Compute the derivatives of the volume basis functions at the quadrature points
    ! of a given surface element. Useful for strain tensors and such

    ! If this segfaults on entry, it's probably because
    ! shape%surface_quadrature is unassociated. You need to augment the shape
    ! function with the quadrature information. See the drag calculation
    ! in MeshDiagnostics.F90 for an example (search for augmented_shape).
    type(element_type), intent(in) :: shape ! NOT the face shape! The volume shape!
    integer, intent(in) :: local_face_number ! which face are we on
    real, dimension(:, :, :), intent(in) :: invJ
    real, dimension(shape%loc, shape%surface_quadrature%ngi, shape%dim) :: output
    integer :: loc, gi

    assert(associated(shape%dn_s))
    assert(size(invJ, 1) == shape%dim)
    assert(size(invJ, 2) == shape%dim)
    assert(size(invJ, 3) == shape%surface_quadrature%ngi)
    assert(shape%dim == size(shape%dn_s, 4))
    assert(shape%loc == size(shape%dn_s, 1))
    assert(shape%surface_quadrature%ngi == size(shape%dn_s, 2))
    assert(local_face_number <= size(shape%dn_s, 3))
    assert(shape%dim == size(shape%dn_s, 4))

    ! You can probably do this with some fancy-pants tensor contraction.
    do loc=1,shape%loc
      do gi=1,shape%surface_quadrature%ngi
        output(loc, gi, :) = matmul(invJ(:, :, gi), shape%dn_s(loc, gi, local_face_number, :))
      end do
    end do
  end function eval_volume_dshape_at_face_quad

  pure function eval_dshape_simplex(shape, loc,  l) result (eval_dshape)
    !!< Evaluate the derivatives of the shape function for location loc at local
    !!< coordinates l 
    !!<
    !!< This version of the function applies to members of the simplex
    !!< family including the interval.
    type(element_type), intent(in) :: shape
    integer, intent(in) :: loc
    real, dimension(shape%dim+1), intent(in) :: l
    real, dimension(shape%dim) :: eval_dshape
    
    integer :: i,j
    ! Derivative of the dependent coordinate with respect to the other
    ! coordinates:
    real, dimension(shape%dim) :: dl4dl

    ! Find derivative of dependent coordinate
    dl4dl=diffl4(shape%numbering%vertices, shape%dim)

    do i=1,shape%dim
       ! Directional derivatives.
       
       ! The derivative has to take into account the dependent
       ! coordinate. In 3D:
       !
       !  S=P1(L1)P2(L2)P3(L3)P4(L4)
       !
       !  dS        / dP1     dL4 dP4  \
       !  --- = P2P3| ---P4 + ---*---P1|
       !  dL1       \ dL1     dL1 dL4  /
       !
       
       ! Expression in brackets.
       eval_dshape(i)=eval(shape%dspoly(i,loc), l(i))&
            *eval(shape%spoly(shape%dim+1,loc),l(shape%dim+1))&
            + dl4dl(i)&
            *eval(shape%dspoly(shape%dim+1,loc), l(shape%dim+1)) &
            *eval(shape%spoly(i,loc),l(i))
             
       ! The other terms
       do j=1,shape%dim
          if (j==i) cycle
          
          eval_dshape(i)=eval_dshape(i)*eval(shape%spoly(j,loc), l(j))
       end do
       
    end do

  end function eval_dshape_simplex

  pure function eval_dshape_cube(shape, loc,  l) result (eval_dshape)
    !!< Evaluate the derivatives of the shape function for location loc at local
    !!< coordinates l 
    !!<
    !!< This version of the function applies to members of the hypercube
    !!< family. Note that this does NOT include the interval.
    type(element_type), intent(in) :: shape
    integer, intent(in) :: loc
    real, dimension(shape%dim+1), intent(in) :: l
    real, dimension(shape%dim) :: eval_dshape

    integer :: i,j

    do i=1,shape%dim
       eval_dshape(i)=1.0
       ! Directional derivatives.
       do j=1,shape%dim
          if(i==j) then
            eval_dshape(i)=eval_dshape(i)*eval(shape%dspoly(j,loc), l(j))
          else
            eval_dshape(i)=eval_dshape(i)*eval(shape%spoly(j,loc), l(j))
          end if          
       end do
    
    end do

  end function eval_dshape_cube

  pure function diffl4(vertices, dimension)
    ! Derivative of the dependent coordinate with respect to the other
    ! coordinates. 
    integer, intent(in) :: vertices, dimension
    real, dimension(dimension) :: diffl4

    if (vertices==dimension+1) then
       ! Simplex. Dependent coordinate depends on all other coordinates. 
       diffl4=-1.0
       
    else if (vertices==2**dimension) then
       ! Hypercube. The dependent coordinate is redundent.
       diffl4=0.0
    
    else if (vertices==6.and.dimension==3) then
       ! Wedge. First coordinate is independent.
       diffl4=(/0.0,-1.0,-1.0/)

    else
       ! No output permitted in a pure procedure so we return a big number to stuff
       ! things up quickly.
       diffl4=huge(0.0)
    end if
       
  end function diffl4

#include "Reference_count_element_type.F90"

end module elements
