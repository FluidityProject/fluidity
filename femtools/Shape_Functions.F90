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
module shape_functions
  !!< Generate shape functions for elements of arbitrary polynomial degree.
  use futils
  use FLDebug
  use polynomials
  use elements
  use element_numbering
  use Superconvergence
  use spud, only: option_count, get_option
  
  implicit none

  private :: lagrange_polynomial, nonconforming_polynomial, shape_integrate&
       &, shape_integrate_diff, monic, cube_monic

  ! Power is used by the test functions.
  integer, private, save :: power=0

  interface make_element_shape
     module procedure make_element_shape_from_element, make_element_shape
  end interface
  
contains

  function make_element_shape_from_element(model, vertices, dim, degree, quad, type,&
       stat, quad_s)  result (shape)
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
    integer, intent(out), optional :: stat
    type(quadrature_type), intent(in), optional, target :: quad_s

    integer :: lvertices, ldim, ldegree
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

    if (associated(lquad_s)) then
       shape = make_element_shape(lvertices, ldim, ldegree, lquad, ltype,&
            stat, lquad_s)
    else
       shape = make_element_shape(lvertices, ldim, ldegree, lquad, ltype,&
            stat)
    end if

  end function make_element_shape_from_element

  function make_element_shape(vertices, dim, degree, quad, type,&
       stat, quad_s)  result (shape)
    !!< Generate the shape functions for an element. The result is a suitable
    !!< element_type.
    !!
    !!< At this stage only Lagrange family polynomial elements are supported.
    type(element_type) :: shape
    !! Vertices is the number of vertices of the element, not the number of nodes!
    !! Only dim=3 is currently supported.
    !! Degree is the degree of the Lagrange polynomials.
    integer, intent(in) :: vertices, dim, degree
    type(quadrature_type), intent(in), target :: quad
    integer, intent(in), optional :: type
    integer, intent(out), optional :: stat
    type(quadrature_type), intent(in), optional, target :: quad_s
    real, pointer :: g(:)=> null()
    
    type(ele_numbering_type), pointer :: ele_num
    ! Count coordinates of each point 
    integer, dimension(dim+1) :: counts
    integer :: i,j,k
    integer :: ltype, coords,surface_count

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
          FLAbort('Element numbering unavailable')
       end if
    end if    

    ! The number of local coordinates depends on the element family.
    select case(ele_num%family)
    case (FAMILY_SIMPLEX)
       coords=dim+1
    case (FAMILY_CUBE)
       coords=dim
    case default
       FLAbort('Illegal element family')
    end select

    shape%numbering=>ele_num
    shape%quadrature=quad
    call incref(quad)

    if (present(quad_s)) then
       select case(dim)
          case(2)
              call allocate(shape, dim, ele_num%nodes, quad%ngi,&
            ele_num%edges, quad_s%ngi,coords,.true.)
              surface_count=ele_num%edges
           case(3)
              call allocate(shape, dim, ele_num%nodes, quad%ngi,&
            ele_num%faces, quad_s%ngi,coords,.true.)
              surface_count=ele_num%edges
           case default
              FLAbort("Unsupported dimension count.")
        end select
    else
       call allocate(shape, dim, ele_num%nodes, quad%ngi, coords)
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
             select case(ele_num%family)
             case (FAMILY_SIMPLEX)
                ! Raw polynomial.
                shape%spoly(j,i)&
                     =lagrange_polynomial(counts(j), counts(j), 1.0/degree)
             case(FAMILY_CUBE)

                ! note that local coordinates run from -1.0 to 1.0
                shape%spoly(j,i)&
                     =lagrange_polynomial(counts(j), degree, 2.0/degree, &
                       origin=-1.0)

             end select

          case(ELEMENT_NONCONFORMING)
             
             shape%spoly(j,i)=nonconforming_polynomial(counts(j))

          case default

             FLAbort('David smoked an illegal element type')
             
          end select

          ! Derivative
          shape%dspoly(j,i)=ddx(shape%spoly(j,i))

       end do

       ! Loop over all the quadrature points.
       do j=1,quad%ngi

          ! Raw shape function
          shape%n(i,j)=eval_shape(shape, i, quad%l(j,:))

          ! Directional derivatives.
          shape%dn(i,j,:)=eval_dshape(shape, i, quad%l(j,:))
       end do

       if (present(quad_s)) then
          shape%surface_quadrature=>quad_s
          select case(ltype)
          case(FAMILY_SIMPLEX)
          allocate(g(dim+1))
           do k=1,dim+1
              do j=1,quad_s%ngi
                 if (dim==2) then
                    g(mod(k+2,3)+1)=0.0
                    g(mod(k,3)+1)=quad_s%l(j,1)
                    g(mod(k+1,3)+1)=quad_s%l(j,2)
                 else if (dim==3) then
                    ! Not checked !!
                    g(mod(k+3,4)+1)=0.0
                    g(mod(k,4)+1)=quad_s%l(j,1)
                    g(mod(k+1,4)+1)=quad_s%l(j,2)
                    g(mod(k+2,4)+1)=quad_s%l(j,3)
                 end if
                 shape%n_s(i,j,k)=eval_shape(shape, i,g)
                 shape%dn_s(i,j,k,:)=eval_dshape(shape, i,g)
             end do
          end do
          deallocate(g)
          end select
      end if

    end do

    shape%superconvergence => get_superconvergence(shape)

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

  !------------------------------------------------------------------------
  ! Test procedures
  !------------------------------------------------------------------------

  function shape_integrate(integrand, element) result (integral)
    !!< Integrate the function integrand over an element using the 
    !!< specified shape functions and quadrature.
    real :: integral
    interface
       function integrand(coords)
         real :: integrand
         real, dimension(:), intent(in) :: coords
       end function integrand
    end interface
    type(element_type), intent(in) :: element

    real :: tmpval
    integer :: i,j 

    integral=0.0

    do i=1, element%loc

       tmpval=integrand(local_coords(i,element))       

       do j=1, element%quadrature%ngi

          integral=integral+element%quadrature%weight(j)*tmpval*element%n(i,j)
          
       end do
    end do
       
  end function shape_integrate

  function shape_integrate_diff(integrand, element, dim) result (integral)
    !!< Integrate the function derivative of integrand with respect to dim
    !!< over an element using the specified shape functions and quadrature.
    real :: integral
    interface
       function integrand(coords)
         real :: integrand
         real, dimension(:), intent(in) :: coords
       end function integrand
    end interface
    type(element_type), intent(in) :: element
    integer, intent(in) :: dim

    real :: tmpval
    integer :: i,j 

    integral=0.0

    do i=1, element%loc
       
       tmpval=integrand(local_coords(i,element))

       do j=1, element%quadrature%ngi
          
          integral=integral&
               +element%quadrature%weight(j)*tmpval*element%dn(i,j,dim)

       end do
    end do

  end function shape_integrate_diff

function shape_integrate_surface(integrand, element, dim,face) &
      result (integral)
    !!< Integrate the function derivative of integrand with respect to dim
    !!< over an element using the specified shape functions and quadrature.
    real :: integral
    interface
       function integrand(coords)
         real :: integrand
         real, dimension(:), intent(in) :: coords
       end function integrand
    end interface
    type(element_type), intent(in) :: element
    integer, intent(in) :: dim
    integer, optional, intent(in) :: face

    real :: tmpval
    integer :: i,j,k

    integral=0.0

    do i=1, element%loc
       
       tmpval=integrand(local_coords(i,element))

       if (present(face)) then
          do j=1, element%surface_quadrature%ngi
             integral=integral&
                  +element%quadrature%weight(j)*tmpval&
                  *element%dn_s(i,j,face,dim)
          end do
       else
          
       do k=1, element%dim+1
          do j=1, element%surface_quadrature%ngi
             integral=integral&
                  +element%quadrature%weight(j)*tmpval*element%n_s(i,j,k)
          end do
       end do
    end if
    end do

  end function shape_integrate_surface

 function shape_integrate_surface_diff(integrand, element, dim,face) &
      result (integral)
    !!< Integrate the function derivative of integrand with respect to dim
    !!< over an element using the specified shape functions and quadrature.
    real :: integral
    interface
       function integrand(coords)
         real :: integrand
         real, dimension(:), intent(in) :: coords
       end function integrand
    end interface
    type(element_type), intent(in) :: element
    integer, intent(in) :: dim
    integer, optional, intent(in) :: face

    real :: tmpval
    integer :: i,j,k

    integral=0.0

    do i=1, element%loc
       
       tmpval=integrand(local_coords(i,element))

       if (present(face)) then
          do j=1, element%surface_quadrature%ngi
             integral=integral&
                  +element%quadrature%weight(j)*tmpval&
                  *element%dn_s(i,j,face,dim)
          end do
       else
          
       do k=1, element%dim+1
          do j=1, element%surface_quadrature%ngi
             integral=integral&
                  +element%quadrature%weight(j)*tmpval*element%dn_s(i,j,k,dim)
          end do
       end do
    end if
    end do

  end function shape_integrate_surface_diff

  subroutine test_shape_functions
    !!< Generic element test function. 
    use unittest_tools
    type(element_type) :: element
    type(quadrature_type) :: quad
    
    character(len=500) :: error_message, test_message
    integer :: dim, degree, vertices
    logical :: fail
    
    ! Rounding error tolerance.
    real, parameter :: eps=1E-12

    ! Test for simplices.
    do dim=1,3

       vertices=dim+1

       quad=make_quadrature(vertices, dim, degree=7)

       do degree=0,7

          element=make_element_shape(vertices=vertices, dim=dim, degree=degree,&
            quad=quad)
            
          do power=0,degree
  
             ! Shape function itself
             if (.not.(abs(shape_integrate(monic, element)&
                  -simplex_answer(power, dim))<eps)) then 
                write(error_message,'(e15.7)') &
                     shape_integrate(monic, element)&
                     &-simplex_answer(power, dim)
                fail=.true.
             else
                error_message=""
                fail=.false.
             end if

             write(test_message, '(3(a,i0),a)') "[",dim,"-simplex, ele&
                  &ment degree ",degree," power ",power," ]"

             call report_test(trim(test_message), fail, .false.,&
                  & trim(error_message))


             ! Derivative
             if (.not.(abs(shape_integrate_diff(monic, element,1)&
                  -power*simplex_answer(power-1, dim))<eps)) then 
                write(error_message,'(e15.7)') &
                     shape_integrate_diff(monic, element,1)&
                     &-power*simplex_answer(power-1, dim)
                fail=.true.
             else
                error_message=""
                fail=.false.
             end if

             write(test_message, '(3(a,i0),a)') "[",dim,"-simplex &
                  &surface derivative, element degree ",&
                  degree," power ",power," ]"

             call report_test(trim(test_message), fail, .false.,&
                  & trim(error_message))


          end do
            
          call deallocate(element)

       end do

       call deallocate(quad)

    end do

    ! Test for hypercubes.
    do dim=2,3

       vertices=2**dim

       quad=make_quadrature(vertices, dim, degree=7)

       do degree=0,7

          element=make_element_shape(vertices=vertices, dim=dim, degree=degree,&
             quad=quad)
          
          do power=0,degree
  
             ! Shape function itself
             if (.not.(abs(shape_integrate(cube_monic, element)&
                  -cube_answer(power, dim))<eps)) then 
                write(error_message,'(e15.7)') &
                     shape_integrate(cube_monic, element)&
                     &-cube_answer(power, dim)
                fail=.true.
             else
                error_message=""
                fail=.false.
             end if

             write(test_message, '(3(a,i0),a)') "[",dim,"-cube, ele&
                  &ment degree ",degree," power ",power," ]"

             call report_test(trim(test_message), fail, .false.,&
                  & trim(error_message))


             ! Derivative
             if (.not.(abs(shape_integrate_diff(cube_monic, element,1)&
                  -1*cube_danswer(power, dim))<eps)) then 
                write(error_message,'(e15.7)') &
                     shape_integrate_diff(cube_monic, element,1)&
                     -1*cube_danswer(power, dim)
                fail=.true.
             else
                error_message=""
                fail=.false.
             end if

             write(test_message, '(3(a,i0),a)') "[",dim,"-cube deri&
                  &vative, element degree ",degree," power ",power," ]"

             call report_test(trim(test_message), fail, .false.,&
                  & trim(error_message))

          end do

          call deallocate(element)
          
       end do

       call deallocate(quad)

    end do    
    
  contains
    
    function simplex_answer(power, dim)
      ! Analytic solution to integrating monic over a simplex.
      ! This formula is eq. 7.38 and 7.48 in Zienkiewicz and Taylor
      real :: simplex_answer
      integer, intent(in) :: power, dim

      simplex_answer=real(factorial(power))&
           /factorial(power+dim)

    end function simplex_answer

    function cube_answer(power, dim)
      ! Analytic solution to integrating ((1-x)/2)**power over a hypercube.
      real :: cube_answer
      integer, intent(in) :: power, dim

      cube_answer=(2.0**dim)/(power+1)

    end function cube_answer

    function cube_danswer(power, dim)
      ! Analytic solution to integrating diff(((1-x)/2)**power,x) over a
      ! hypercube.
      real :: cube_danswer
      integer, intent(in) :: power, dim
      
      if (power==0) then
         cube_danswer=0
      else
         cube_danswer=-1*2**(dim-1)
      end if

    end function cube_danswer
      

    recursive function factorial(n) result (f)
      ! Calculate n!
      integer :: f
      integer, intent(in) :: n

      if (n==0) then
         f=1
      else if (n<0) then
         f=0
      else
         f=n*factorial(n-1)
      end if

    end function factorial
    
  end subroutine test_shape_functions

  subroutine test_element_surface_integral
    use unittest_tools
    implicit none

    type(element_type) :: element
    type(quadrature_type) :: quad, quad_s
    
    character(len=500) :: error_message, test_message
    integer :: dim, degree, vertices, power, k
    logical :: fail
    ! Rounding error tolerance.
    real, parameter :: eps=1E-13

    dim=element%dim

    ! Test for simplices.
    do dim=2,2

       vertices=dim+1
       degree=0

       do degree=0,7

       quad=make_quadrature(vertices, dim, degree=7)
       quad_s=make_quadrature(vertices-1, dim-1, degree=degree)
       element=make_element_shape(vertices=vertices, dim=dim, degree=degree,&
            quad=quad, quad_s=quad_s)

          do power=0,degree

             ! surface 
             if (.not.(abs(shape_integrate_surface(monic, element,1)&
                  )<eps)) then 
                write(error_message,'(e15.7)') &
                     shape_integrate_surface(monic, element,1)
                fail=.false.
             else
                error_message=""
                fail=.false.
             end if

             write(test_message, '(3(a,i0),a)') "[",dim,"-simplex &
                  &surface, element degree ",&
                  degree," power ",power," ]"

             call report_test(trim(test_message), fail, .false.,&
                  & trim(error_message))

             do k=1,dim+1

                if (.not.(abs(shape_integrate_surface(monic, element,1,k)&
                     -line_answer(power, dim,k))<eps)) then 
                   write(error_message,'(e15.7)') &
                        shape_integrate_surface(monic, element,1,k)&
                        -line_answer(power, dim,k)
                   fail=.false.
                else
                   error_message=""
                   fail=.false.
                end if

                write(test_message, '(4(a,i0),a)') "[",dim,"-simplex &
                     &surface face,", k, " element degree ",&
                     degree," power ",power," ]"

                call report_test(trim(test_message), fail, .false.,&
                     & trim(error_message))

             end do

             ! surface derivative
             if (.not.(abs(shape_integrate_surface_diff(monic, element,1)&
                  )<eps)) then 
                write(error_message,'(e15.7)') &
                     shape_integrate_surface_diff(monic, element,1)
                fail=.false.
             else
                error_message=""
                fail=.false.
             end if

             write(test_message, '(3(a,i0),a)') "[",dim,"-simplex &
                  &surface derivative, element degree ",&
                  degree," power ",power," ]"

             call report_test(trim(test_message), fail, .false.,&
                  & trim(error_message))

             do k=1,dim+1

                if (.not.(abs(shape_integrate_surface_diff(monic, element,1,k)&
                     -power*line_answer(power-1, dim,k))<eps)) then 
                   write(error_message,'(e15.7)') &
                        shape_integrate_surface_diff(monic, element,1,k)&
                        -power*line_answer(power-1, dim,k)
                   fail=.false.
                else
                   error_message=""
                   fail=.false.
                end if

                write(test_message, '(4(a,i0),a)') "[",dim,"-simplex &
                     &surface derivative face,", k, " element degree ",&
                     degree," power ",power," ]"

                call report_test(trim(test_message), fail, .false.,&
                     & trim(error_message))

             end do


          end do

             call deallocate(quad)
             call deallocate(quad_s)
             call deallocate(element)

          end do


       end do

      contains

       function line_answer(power, dim, faces)
      ! Analytic solution to integrating monic on a line.
      real :: line_answer
      integer, intent(in) :: power, dim, faces

      if (faces==3) then
         line_answer=(1.0)/(power+1)
      else
         line_answer=-0.5*(1.0)/(power+1)
      end if
         

    end function line_answer

     recursive function factorial(n) result (f)
      ! Calculate n!
      integer :: f
      integer, intent(in) :: n

      if (n==0) then
         f=1
      else if (n<0) then
         f=0
      else
         f=n*factorial(n-1)
      end if

    end function factorial
    
  end subroutine test_element_surface_integral

  function monic(coords)
    !!< Calculate x^n
    real :: monic
    real, dimension(:), intent(in) :: coords
    
    monic=coords(1)**power
    
  end function monic
  
  function cube_monic(coords)
    ! Calculate.
    real :: cube_monic
    real, dimension(:), intent(in) :: coords

    cube_monic=((1-coords(1))/2.0)**power

  end function cube_monic

  subroutine shape_functions_check_options

    integer :: quaddegree, degree, stat, i, nmesh

    call get_option("/geometry/quadrature/degree", quaddegree)
    nmesh = option_count("/geometry/mesh")
    do i = 1, nmesh
      call get_option("/geometry/mesh["//int2str(i-1)//&
                      "]/from_mesh/mesh_shape/polynomial_degree", &
                      degree, stat)
      if(stat==0) then
        if (quaddegree<2*degree) then
          ewrite(0,"(a,i0,a,i0)") "Warning: quadrature of degree ",quaddegree&
                &," may be incomplete for elements of degree ",degree
        end if
      end if
    end do


  end subroutine shape_functions_check_options

end module shape_functions
