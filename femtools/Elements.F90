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
     !! n_s is loc x sngi, dn_s is loc x sngi x dim
     !! NOTE that both n_s and dn_s need to be reoriented before use so that they align with the arbitrary facet node ordering.
     real, pointer :: n(:,:)=>null(), dn(:,:,:)=>null()
     real, pointer :: n_s(:,:)=>null(), dn_s(:,:,:)=>null()
     !! Polynomials defining shape functions and their derivatives.
     type(polynomial), dimension(:,:), pointer :: spoly=>null(), dspoly=>null()
     !! Link back to the node numbering used for this element.
     type(ele_numbering_type), pointer :: numbering=>null()
     !! Link back to the quadrature used for this element.
     type(quadrature_type) :: quadrature
     type(quadrature_type), pointer :: surface_quadrature=>null()
     !! Pointer to the superconvergence data for this element.
     type(superconvergence_type), pointer :: superconvergence=>null()
     !! Pointer to constraints data for this element
     type(constraints_type), pointer :: constraints=>null()
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

  type constraints_type
     !!< A type to encode the constraints from the local Lagrange basis for 
     !!< (Pn)^d vector-valued elements to another local basis, possibly for 
     !!< a proper subspace. This new basis must have DOFs consisting
     !!< of either normal components on faces corresponding to a Lagrange
     !!< basis for the normal component when restricted to each face,
     !!< or coefficients of basis
     !!< functions with vanishing normal components on all faces.
     !! type of constraints
     integer :: type
     !! local dimension
     integer :: dim
     !! order of local Lagrange basis
     integer :: degree
     !! number of nodes for local Lagrange basis
     integer :: loc
     !! Number of constraints
     integer :: n_constraints
     !! basis of functions that are orthogonal to the 
     !! constrained vector space 
     !! dimension n_constraints x loc x dim
     real, pointer :: orthogonal(:,:,:)=> null()
  end type constraints_type

  integer, parameter :: CONSTRAINT_NONE =0, CONSTRAINT_BDFM = 1,&
       & CONSTRAINT_RT = 2, CONSTRAINT_BDM = 3

  interface allocate
     module procedure allocate_element
     module procedure allocate_constraints_type
  end interface

  interface deallocate
     module procedure deallocate_element
     module procedure deallocate_constraints
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

  subroutine allocate_element(element, ele_num, ngi, ngi_s, type, stat)
    !!< Allocate memory for an element_type. 
    type(element_type), intent(inout) :: element
    !! Number of quadrature points
    integer, intent(in) :: ngi    
    !! Element numbering
    type(ele_numbering_type), intent(in) :: ele_num
    !! Stat returns zero for success and nonzero otherwise.
    integer, intent(in), optional :: ngi_s
    integer, intent(in), optional :: type
    integer, intent(out), optional :: stat
    !
    integer :: lstat, coords, ltype

    if(present(type)) then
       ltype = type
    else
       ltype = ele_num%type
    end if

    select case(ele_num%family)
    case (FAMILY_SIMPLEX)
       coords=ele_num%dimension+1
    case (FAMILY_CUBE)
       if(ele_num%type==ELEMENT_TRACE .and. ele_num%dimension==2) then
          !For trace elements the local coordinate is face number
          !then the local coordinates on the face
          !For quads, the face is an interval element which has
          !two local coordinates.
          coords=3
       else
          coords=ele_num%dimension
       end if
    case default
       FLAbort('Illegal element family.')
    end select

    select case(ltype)
    case(ELEMENT_LAGRANGIAN, ELEMENT_NONCONFORMING, &
         &ELEMENT_BUBBLE, ELEMENT_TRACE)

       allocate(element%n(ele_num%nodes,ngi),&
            &element%dn(ele_num%nodes,ngi,ele_num%dimension), &
            &element%spoly(coords,ele_num%nodes), &
            &element%dspoly(coords,ele_num%nodes), stat=lstat)

    case(ELEMENT_CONTROLVOLUME_SURFACE)

       allocate(element%n(ele_num%nodes,ngi),&
            &element%dn(ele_num%nodes,ngi,ele_num%dimension-1), &
            stat=lstat)

      element%spoly=>null()
      element%dspoly=>null()

    case(ELEMENT_CONTROLVOLUMEBDY_SURFACE)

      allocate(element%n(ele_num%nodes,ngi),&
           &element%dn(ele_num%nodes,ngi,ele_num%dimension), &
          stat=lstat)

      element%spoly=>null()
      element%dspoly=>null()

    case default

      FLAbort("Attempt to select an illegal element type.")

    end select

    element%loc=ele_num%nodes
    element%ngi=ngi
    element%dim=ele_num%dimension

    if (present(ngi_s)) then
      allocate(element%n_s(ele_num%nodes,ngi_s), &
               element%dn_s(ele_num%nodes,ngi_s,ele_num%dimension), stat=lstat)
    else
      nullify(element%n_s)
      nullify(element%dn_s)
    end if

    nullify(element%refcount) ! Hack for gfortran component initialisation
    !                         bug.
    call addref(element)

    if (present(stat)) then
       stat=lstat
    else if (lstat/=0) then
       FLAbort("Unable to allocate element.")
    end if

  end subroutine allocate_element

  subroutine allocate_constraints_type(constraint, element, type, stat)
    !!< Allocate memory for a constraints type
    type(element_type), intent(in) :: element
    type(constraints_type), intent(inout) :: constraint
    integer, intent(in) :: type !type of constraint
    !! Stat returns zero for success and nonzero otherwise.
    integer, intent(out), optional :: stat
    !
    integer :: lstat

    lstat = 0
    constraint%type = type
    constraint%dim = element%dim
    constraint%loc = element%loc
    constraint%degree = element%degree

    select case (type) 
    case (CONSTRAINT_BDFM)
       select case(element%numbering%family)
       case (FAMILY_SIMPLEX)
          if(constraint%degree<3) then
             constraint%n_constraints = constraint%dim+1
          else
             FLAbort('High order not supported yet')
          end if
       case (FAMILY_CUBE)
          FLExit('Haven''t implemented BDFM1 on quads yet.')
       case default
          FLAbort('Illegal element family.')
       end select
    case (CONSTRAINT_RT)
       select case(element%numbering%family)
       case (FAMILY_SIMPLEX)
          FLExit('Haven''t implemented RT0 on simplices yet.')
          if(constraint%degree<3) then
             constraint%n_constraints = constraint%dim+1
          else
             FLAbort('High order not supported yet')
          end if
       case (FAMILY_CUBE)
          if(constraint%degree<3) then
             constraint%n_constraints = 2**(constraint%dim)
          else
             FLAbort('High order not supported yet')
          end if
       case default
          FLAbort('Illegal element family.')
       end select
    case (CONSTRAINT_BDM)
       constraint%n_constraints = 0
    case (CONSTRAINT_NONE)
       constraint%n_constraints = 0
    case default
       FLExit('Unknown constraint type')
    end select

    if(constraint%n_constraints>0) then
       allocate(&
            constraint%orthogonal(constraint%n_constraints,&
            constraint%loc,constraint%dim),stat=lstat)
       if(lstat==0) then
          call make_constraints(constraint,element%numbering%family)
       end if
    end if

    if (present(stat)) then
       stat=lstat
    else if (lstat/=0) then
       FLAbort("Unable to allocate element.")
    end if

  end subroutine allocate_constraints_type

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

    if (associated(element%surface_quadrature)) then
      call deallocate(element%surface_quadrature)
      deallocate(element%surface_quadrature, stat=tstat)
    end if
    lstat=max(lstat,tstat)

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

    if(associated(element%constraints)) then
       call deallocate(element%constraints,stat=tstat)
       lstat = max(lstat,tstat)

       deallocate(element%constraints, stat=tstat)
       lstat = max(lstat,tstat)
    end if
    if (present(stat)) then
       stat=lstat
    else if (lstat/=0) then
       FLAbort("Unable to deallocate element.")
    end if

  end subroutine deallocate_element

  subroutine deallocate_constraints(constraint, stat)
    type(constraints_type), intent(inout) :: constraint
    integer, intent(out), optional :: stat
    
    integer :: lstat

    lstat = 0

    if(associated(constraint%orthogonal)) then
       deallocate(constraint%orthogonal,stat=lstat)
    end if

    if (present(stat)) then
       stat=lstat
    else if (lstat/=0) then
       FLAbort("Unable to deallocate constraints.")
    end if

  end subroutine deallocate_constraints

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
       ! Hypercube. The dependent coordinate is redundant.
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

  subroutine make_constraints(constraint,family)
    type(constraints_type), intent(inout) :: constraint
    integer, intent(in) :: family
    !
    select case(family)
    case (FAMILY_SIMPLEX)
       select case(constraint%type)
       case (CONSTRAINT_BDM)
          !do nothing
       case (CONSTRAINT_BDFM)
          select case(constraint%dim)
          case (2)
             select case(constraint%degree)
             case (1)
                !BDFM0 is the same as RT0
                call make_constraints_rt0_triangle(constraint)
             case (2)
                call make_constraints_bdfm1_triangle(constraint)
             case default
                FLExit('Unknown constraints type')
             end select
          case default
             FLExit('Unsupported dimension')
          end select
       case (CONSTRAINT_RT)
          select case(constraint%dim)
          case (2)
             select case(constraint%degree)
             case (1)
                call make_constraints_rt0_triangle(constraint)
             case default
                FLExit('Unknown constraints type')
             end select
          case default
             FLExit('Unsupported dimension')
          end select
       case default
          FLExit('Unknown constraints type')
       end select
    case (FAMILY_CUBE)
       select case(constraint%type)
       case (CONSTRAINT_BDM)
          !do nothing
       case (CONSTRAINT_RT)
          select case(constraint%dim)
          case (2)
             select case(constraint%degree)
             case (1)
                call make_constraints_rt0_square(constraint)
             case (2)
                FLExit('Haven''t implemented it yet!')
                !call make_constraints_rt1_square(constraint)
             case default
                FLExit('Unknown constraints type')
             end select
          case default
             FLExit('Unsupported dimension')
          end select
       case default
          FLExit('Unknown constraints type')
       end select
    case default
       FLExit('Unknown element numbering family')
    end select
  end subroutine make_constraints

  subroutine make_constraints_bdfm1_triangle(constraint)
    implicit none
    type(constraints_type), intent(inout) :: constraint
    real, dimension(3,2) :: n
    integer, dimension(3,3) :: face_loc
    integer :: dim1, face, floc
    real, dimension(3) :: c

    if(constraint%dim/=2) then
       FLExit('Only implemented for 2D so far')
    end if

    !BDFM1 constraint requires that normal components are linear.
    !This means that the normal components at the edge centres
    !need to be constrained to the average of the normal components 
    !at each end of the edge.

    !DOFS    FACES
    ! 3      
    ! 5 2    1 3
    ! 6 4 1   2

    !constraint equations are:
    ! (0.5 u_3 - u_5 + 0.5 u_6).n_1 = 0
    ! (0.5 u_1 - u_4 + 0.5 u_6).n_2 = 0    
    ! (0.5 u_1 - u_2 + 0.5 u_3).n_3 = 0

    !face local nodes to element local nodes
    face_loc(1,:) = (/ 3,5,6 /)
    face_loc(2,:) = (/ 1,4,6 /)
    face_loc(3,:) = (/ 1,2,3 /)

    !normals
    n(1,:) = (/ -1., 0. /)
    n(2,:) = (/  0.,-1. /)
    n(3,:) = (/ 1./sqrt(2.),1./sqrt(2.) /)

    !coefficients in each face
    c = (/ 0.5,-1.,0.5 /)

    !constraint%orthogonal(i,loc,dim1) stores the coefficient 
    !for basis function loc, dimension dim1 in equation i.

    constraint%orthogonal = 0.
    do face = 1, 3
       do floc = 1,3
          do dim1 = 1, 2
             constraint%orthogonal(face,face_loc(face,floc),dim1) = &
                  c(floc)*n(face,dim1)
          end do
       end do
    end do
    !! dimension n_constraints x loc x dim
  end subroutine make_constraints_bdfm1_triangle

  subroutine make_constraints_rt0_triangle(constraint)
    implicit none
    type(constraints_type), intent(inout) :: constraint
    real, dimension(3,2) :: n
    integer, dimension(3,2) :: face_loc
    integer :: dim1, face, floc, count
    real, dimension(2) :: c

    if(constraint%dim/=2) then
       FLExit('Only implemented for 2D so far')
    end if

    !RT0 constraint requires that normal components are constant.
    !This means that both the normal components at each end of the 
    !edge need to have the same value.

    !DOFS    FACES
    ! 2      
    !        1 3
    ! 3   1   2

    !constraint equations are:
    ! (u_2 - u_3).n_1 = 0
    ! (u_1 - u_3).n_2 = 0    
    ! (u_1 - u_2).n_3 = 0

    !face local nodes to element local nodes
    face_loc(1,:) = (/ 2,3 /)
    face_loc(2,:) = (/ 1,3 /)
    face_loc(3,:) = (/ 1,2 /)

    !normals
    n(1,:) = (/ -1., 0. /)
    n(2,:) = (/  0.,-1. /)
    n(3,:) = (/ 1./sqrt(2.),1./sqrt(2.) /)

    !constraint coefficients
    c = (/ 1., -1. /)

    !constraint%orthogonal(i,loc,dim1) stores the coefficient 
    !for basis function loc, dimension dim1 in equation i.

    constraint%orthogonal = 0.
    count = 0
    do face = 1, 3
       count = count + 1
       do floc = 1,2
          do dim1 = 1, 2
             constraint%orthogonal(count,face_loc(face,floc),dim1)&
                  = c(floc)*n(face,dim1)
          end do
       end do
    end do
    assert(count==3)
    !! dimension n_constraints x loc x dim
  end subroutine make_constraints_rt0_triangle

  subroutine make_constraints_rt0_square(constraint)
    implicit none
    type(constraints_type), intent(inout) :: constraint
    real, dimension(4,2) :: n
    integer, dimension(4,2) :: face_loc
    integer :: dim1, face, floc, count
    real, dimension(2) :: c

    if(constraint%dim/=2) then
       FLExit('Only implemented for 2D so far')
    end if

    !RT0 constraint requires that normal components are constant.
    !This means that both the normal components at each end of the 
    !edge need to have the same value.

    !DOFS    FACES
    ! 3   4   3
    !        4 2
    ! 1   2   1

    !constraint equations are:
    ! (u_1 - u_2).n_1 = 0
    ! (u_2 - u_4).n_2 = 0    
    ! (u_3 - u_4).n_3 = 0
    ! (u_3 - u_1).n_4 = 0

    !face local nodes to element local nodes
    face_loc(1,:) = (/ 1,2 /)
    face_loc(2,:) = (/ 2,4 /)
    face_loc(3,:) = (/ 3,4 /)
    face_loc(4,:) = (/ 3,1 /)

    !normals
    n(1,:) = (/  0., -1. /)
    n(2,:) = (/  1.,  0. /)
    n(3,:) = (/  0.,  1. /)
    n(4,:) = (/ -1.,  0. /)

    !constraint%orthogonal(i,loc,dim1) stores the coefficient 
    !for basis function loc, dimension dim1 in equation i.

    !constraint coefficients
    c = (/ 1., -1. /)

    constraint%orthogonal = 0.
    count  = 0
    do face = 1, 4
       count = count + 1
       do floc = 1,2
          do dim1 = 1, 2
             constraint%orthogonal(count,face_loc(face,floc),dim1)&
                  = c(floc)*n(face,dim1)
          end do
       end do
    end do
    assert(count==4)
    !! dimension n_constraints x loc x dim
  end subroutine make_constraints_rt0_square

#include "Reference_count_element_type.F90"

end module elements
