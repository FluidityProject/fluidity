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
module shape_functions
  !!< Generate shape functions for elements of arbitrary polynomial degree.
  use futils
  use FLDebug
  use polynomials
  use elements
  use element_numbering
  use Superconvergence
  use ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  
  implicit none

  private :: lagrange_polynomial, nonconforming_polynomial

  interface make_element_shape
     module procedure make_element_shape_from_element, make_element_shape
  end interface
  
contains

  function make_element_shape_from_element(model, vertices, dim, degree,&
       & quad, type, quad_s, constraint_type_choice, stat)  result (shape)
    !!< This function enables element shapes to be derived from other
    !!< element shapes by specifying which attributes to change.
    type(element_type) :: shape
    type(element_type), intent(in) :: model
    !! Vertices is the number of vertices of the element, not the number of nodes!
    !! dim may be 1, 2, or 3.
    !! Degree is the degree of the Lagrange polynomials.
    integer, intent(in), optional :: vertices, dim, degree
    type(quadrature_type), intent(in), target, optional :: quad
    integer, intent(in), optional :: type
    type(quadrature_type), intent(in), optional, target :: quad_s
    !! Element constraints
    integer, intent(in), optional :: constraint_type_choice
    integer, intent(out), optional :: stat

    integer :: lvertices, ldim, ldegree, lconstraint_type_choice
    type(quadrature_type) :: lquad
    type(quadrature_type), pointer :: lquad_s
    integer :: ltype

    if (present(vertices)) then
       lvertices=vertices
    else
       lvertices=model%numbering%vertices
    end if

    if (present(dim)) then
       ldim=dim
    else
       ldim=model%dim
    end if

    if(present(degree)) then
       ldegree=degree
    else
       ldegree=model%degree
    end if

    if(present(quad)) then
       lquad=quad
    else
       lquad=model%quadrature
    end if

    if(present(type)) then
       ltype=type
    else
       ltype=model%numbering%type
    end if

    if(present(quad_s)) then
       lquad_s=>quad_s
    else if (associated(model%surface_quadrature)) then
       lquad_s=>model%surface_quadrature
    else
       lquad_s=>null()
    end if

    if(present(constraint_type_choice)) then
       lconstraint_type_choice=constraint_type_choice
    else if (associated(model%constraints)) then
       lconstraint_type_choice=model%constraints%type
    else
       lconstraint_type_choice=CONSTRAINT_NONE
    end if

    
    if (associated(lquad_s)) then
       shape = make_element_shape(lvertices, ldim, ldegree, lquad, ltype,&
            lquad_s, constraint_type_choice=lconstraint_type_choice, stat=stat)
    else
       shape = make_element_shape(lvertices, ldim, ldegree, lquad, ltype,&
            constraint_type_choice=lconstraint_type_choice, stat=stat)
    end if

  end function make_element_shape_from_element

  function make_element_shape(vertices, dim, degree, quad, type,&
       quad_s, constraint_type_choice, stat)  result (shape)
    !!< Generate the shape functions for an element. The result is a suitable
    !!< element_type.
    !!
    !!< At this stage only Lagrange family polynomial elements are supported.
    type(element_type) :: shape
    !! Vertices is the number of vertices of the element, not the number of nodes!
    !! dim \in [1,2,3] is currently supported.
    !! Degree is the degree of the Lagrange polynomials.
    integer, intent(in) :: vertices, dim, degree
    type(quadrature_type), intent(in), target :: quad
    integer, intent(in), optional :: type
    type(quadrature_type), intent(in), optional, target :: quad_s
    integer, intent(in), optional :: constraint_type_choice
    integer, intent(out), optional :: stat

    real, pointer :: g(:)=> null()

    type(ele_numbering_type), pointer :: ele_num
    ! Count coordinates of each point 
    integer, dimension(dim+1) :: counts
    integer :: i,j,k
    integer :: ltype, coords
    real :: dx
    type(constraints_type), pointer :: constraint

    ! Check that the quadrature and the element shapes match.
    assert(quad%vertices==vertices)
    assert(quad%dim==dim)

    if (present(type)) then
       ltype=type
    else
       ltype=ELEMENT_LAGRANGIAN
    end if

    if (present(stat)) stat=0

    ! Get the local numbering of our element
    ele_num=>find_element_numbering(vertices, dim, degree, type)

    if (.not.associated(ele_num)) then
       if (present(stat)) then
          stat=1
          return
       else
          FLAbort('Element numbering unavailable.')
       end if
    end if

    shape%numbering=>ele_num
    shape%quadrature=quad
    call incref(quad)

    ! The number of local coordinates depends on the element family.
    select case(ele_num%family)
    case (FAMILY_SIMPLEX)
       coords=dim+1
    case (FAMILY_CUBE)
       if(ele_num%type==ELEMENT_TRACE .and. dim==2) then
          !For trace elements the local coordinate is face number
          !then the local coordinates on the face
          !For quads, the face is an interval element which has
          !two local coordinates.
          coords=3
       else
          coords=dim
       end if
    case default
       FLAbort('Illegal element family.')
    end select

    if (present(quad_s) .and. ele_num%type/=ELEMENT_TRACE .and. ele_num%family==FAMILY_SIMPLEX) then
       allocate(shape%surface_quadrature)
       shape%surface_quadrature=quad_s
       call incref(quad_s)
       call allocate(shape, ele_num, quad%ngi, ngi_s=quad_s%ngi)
       shape%n_s=0.0
       shape%dn_s=0.0
    else
       call allocate(shape, ele_num, quad%ngi)
    end if
    shape%degree=degree
    shape%n=0.0
    shape%dn=0.0

    ! Construct shape for each node
    do i=1,shape%loc

       counts(1:coords)=ele_num%number2count(:,i)

       ! Construct appropriate polynomials.
       do j=1,coords

          select case(ltype)
          case(ELEMENT_LAGRANGIAN)
             if (degree == 0) then
                dx = 0.0
             else
                dx = 1.0/degree
             end if 
             select case(ele_num%family)
             case (FAMILY_SIMPLEX)
                ! Raw polynomial.
                shape%spoly(j,i)&
                     =lagrange_polynomial(counts(j), counts(j), dx)
             case(FAMILY_CUBE)
                ! note that local coordinates run from -1.0 to 1.0
                shape%spoly(j,i)&
                     =lagrange_polynomial(counts(j), degree, 2.0*dx, &
                     origin=-1.0)
             end select

          case(ELEMENT_TRACE)
             shape%spoly(j,i) = (/ieee_value(0.0,ieee_quiet_nan)/)

          case(ELEMENT_BUBBLE)
             if(i==shape%loc) then

                ! the last node is the bubble shape function
                shape%spoly(j,i) = (/1.0, 0.0/)

             else

                select case(ele_num%family)
                case (FAMILY_SIMPLEX)
                   ! Raw polynomial.
                   shape%spoly(j,i)&
                        =lagrange_polynomial(counts(j)/coords, counts(j)/coords, 1.0/degree)

                end select

             end if

          case(ELEMENT_NONCONFORMING)

             shape%spoly(j,i)=nonconforming_polynomial(counts(j))

          case default

             FLAbort('An unsupported element type has been selected.')

          end select

          ! Derivative
          if(ele_num%type==ELEMENT_TRACE) then
             shape%dspoly(j,i) = (/ieee_value(0.0,ieee_quiet_nan)/)
          else
             shape%dspoly(j,i)=ddx(shape%spoly(j,i))
          end if
       end do

       if(ele_num%type==ELEMENT_TRACE) then
          !No interior functions, hence NaNs
          shape%n = ieee_value(0.0,ieee_quiet_nan)
          shape%dn = ieee_value(0.0,ieee_quiet_nan)
       else
          ! Loop over all the quadrature points.
          do j=1,quad%ngi

             ! Raw shape function
             shape%n(i,j)=eval_shape(shape, i, quad%l(j,:))

             ! Directional derivatives.
             shape%dn(i,j,:)=eval_dshape(shape, i, quad%l(j,:))
          end do

          if (present(quad_s)) then
             select case(ele_num%family)
             case(FAMILY_SIMPLEX)
                allocate(g(dim+1))
                do j=1,quad_s%ngi
                   g(1) = 0.0
                   do k=1,dim
                      g(k+1)=quad_s%l(j,k)
                   end do
                   ! In order to match the arbitrary face node ordering
                   ! these must get reoriented before use so we don't care
                   ! about which local facet they're with respect to.
                   shape%n_s(i,j)=eval_shape(shape,i,g)
                   shape%dn_s(i,j,:)=eval_dshape(shape,i,g)
                end do
                deallocate(g)
             end select
          end if
       end if
    end do

    if(ele_num%type.ne.ELEMENT_TRACE) then
       shape%superconvergence => get_superconvergence(shape)
    end if

    if(present(constraint_type_choice)) then
       if(constraint_type_choice/=CONSTRAINT_NONE) then
          allocate(constraint)
          shape%constraints=>constraint
          call allocate(shape%constraints,shape,constraint_type_choice)
       end if
    end if

  end function make_element_shape

  function lagrange_polynomial(n,degree,dx, origin) result (poly)
    ! nth equispaced lagrange polynomial of specified degree and point
    ! spacing dx.
    integer, intent(in) :: n, degree
    real, intent(in) :: dx
    type(polynomial) :: poly
    ! fixes location of n=0 location (0.0 if not specified)
    real, intent(in), optional :: origin

    real lorigin
    integer :: i

    ! This shouldn't be necessary but there appears to be a bug in initial
    ! component values in gfortran:
    poly%coefs=>null()
    poly%degree=-1
    
    if (present(origin)) then
       lorigin=origin
    else
       lorigin=0.0
    end if

    poly=(/1.0/)
    
    degreeloop: do i=0,degree
       if (i==n) cycle degreeloop

       poly=poly*(/1.0, -(lorigin+i*dx) /)

    end do degreeloop

    ! normalize to 1.0 in the n-th location
    poly=poly/eval(poly, lorigin+n*dx)
    
  end function lagrange_polynomial

  function nonconforming_polynomial(n) result (poly)
    ! nth P1 nonconforming polynomial.
    integer, intent(in) :: n
    type(polynomial) :: poly

    ! This shouldn't be necessary but there appears to be a bug in initial
    ! component values in gfortran:
    poly%coefs=>null()
    poly%degree=-1

    poly=(/1.0/)    

    if (n==0) then

       ! polynomial is -2x+1
       poly=(/-2.0, 1.0/)

    end if
       
  end function nonconforming_polynomial

end module shape_functions
