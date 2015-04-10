!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineeringp
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
module distance_function

  use fldebug
  use state_module
  use fields
  use field_derivatives, only : grad

  use spud

  use global_parameters, only: OPTION_PATH_LEN

implicit none

private

public :: marching_distance_function, hamilton_jacobi_distance_function


type marching_logic
   logical, dimension(:), allocatable :: trial
   logical, dimension(:), allocatable :: accepted
end type marching_logic

type marching_type
   type(marching_logic) :: ele
   type(marching_logic)  :: node
end type marching_type

contains

  subroutine hamilton_jacobi_distance_function(sfield,positions,itrs)

    type(scalar_field), intent(inout)     :: sfield
    type(vector_field), intent(in)     :: positions
    integer, intent(in), optional :: itrs

    type(scalar_field)  :: distance(3), mesh_length, inverse_volume, u_dot_dx
    type(vector_field) :: gradient(4), DiscontinuousGradient
    integer :: itr, ele, node, dim, l_itrs
    integer, dimension(:), pointer :: nodes
    real :: dt, min_dx, vol

    type(mesh_type) :: meshA

    ! Subroutine implements a version of the Hamilton Jacobi algorithm to
    ! approximate a distance function from the zero contour of the input field,
    ! sfield, with a coordinate mesh described by positions

    if (present(itrs)) then
       l_itrs=itrs
    else
       l_itrs=10
    end if

    dim=mesh_dim(sfield)

    meshA = make_mesh(positions%mesh, sfield%mesh%shape, -1, 'DiscontinuousMesh')

    call allocate(mesh_length,sfield%mesh,name="MeshSize")
    call allocate(inverse_volume,sfield%mesh,name="Volume")
    call allocate(u_dot_dx,sfield%mesh,name="u_dot_dx")
    call allocate(DiscontinuousGradient,dim=dim,mesh=meshA,&
         name='DiscontinousGradient')
    call zero(mesh_length)
    call zero(inverse_volume)
    do ele=1,ele_count(sfield)
       nodes=>ele_nodes(sfield,ele)
       vol=element_volume(Positions,ele)
       do itr=1,size(nodes)
          call addto(inverse_volume,nodes(itr),vol)
          call addto(mesh_length,nodes(itr),vol)
       end do
    end do
    mesh_length%val=mesh_length%val**(1.0/(1.0*dim))


    call halo_update(mesh_length)
    call halo_update(inverse_volume)
    call invert(inverse_volume)

    do itr=1,3
       call allocate(distance(itr),sfield%mesh,name="Distance")
       call allocate(gradient(itr),dim=dim,mesh=sfield%mesh,name="Gradient")
    end do
    call allocate(gradient(4),dim=dim,mesh=sfield%mesh,name="Gradient")

    call pgrad(sfield,positions,gradient(4))

    call set(distance(1),sfield)

    do itr=1,l_itrs
       ewrite(2,*) itr, maxval(abs(distance(1)%val)), sum(distance(1)%val)/size(distance(1)%val)
       call pgrad(distance(1),positions,gradient(1))

       call calculate_timestep(sfield,gradient(1),positions, u_dot_dx,dt)


       do node=1,node_count(sfield)
          call set(distance(2),node,node_val(distance(1),node)+dt*L(distance(1),gradient(1),node))
       end do
       call halo_update(distance(2))


       call pgrad(distance(2),positions,gradient(2))
       do node=1,node_count(sfield)
          call set(distance(3),node,3.0/4.0*node_val(distance(1),node)&
          +1.0/4.0*node_val(distance(2),node)&
          +dt/4.0*L(distance(2),gradient(2),node))
       end do
       call halo_update(distance(3))

       
       call pgrad(distance(3),positions,gradient(3))
       do node=1,node_count(sfield)
          call set(distance(1),node,1.0/3.0*node_val(distance(1),node)&
               +2.0/3.0*node_val(distance(3),node)&
               +4.0*dt/6.0*L(distance(3),gradient(3),node))
       end do
       call halo_update(distance(1))

    end do

    
    call set(sfield,distance(1))

    do itr=1,3
       call deallocate(distance(itr))
       call deallocate(gradient(itr))
    end do
    call deallocate(gradient(4))

    call deallocate(mesh_length)
    call deallocate(DiscontinuousGradient)
    call deallocate(inverse_volume)
    call deallocate(meshA)
    call deallocate(u_dot_dx)

  contains 

    subroutine calculate_timestep(psi,grad_psi, X, u_dot_dx, dt)
    !! Calculate the CFL number as a field.
    type(scalar_field), intent(in) :: psi 
    type(vector_field), intent(inout) :: grad_psi
    type(vector_field), intent(in) :: X
    type(scalar_field), intent(inout) :: u_dot_dx
    real :: dt

    integer :: ele, gi
    ! Transformed quadrature weights.
    real, dimension(ele_ngi(grad_psi, 1)) :: detwei
    ! Inverse of the local coordinate change matrix.
    real, dimension(mesh_dim(X), mesh_dim(X), ele_ngi(X, 1)) :: invJ
    ! velocity/dx at each quad point.
    real, dimension(mesh_dim(X), ele_ngi(X, 1)) :: CFL_q
    ! current element global node numbers.
    integer, dimension(:), pointer :: ele_cfl
    ! local cfl matrix on the current element.
    real, dimension(ele_loc(grad_psi, 1),ele_loc(grad_psi, 1)) :: CFL_mat
    ! current CFL element shape
    type(element_type), pointer :: CFL_shape

    real, parameter :: tol=1.0e-16


    do ele=1, element_count(grad_psi)
       ele_CFL=>ele_nodes(grad_psi, ele)
       CFL_shape=>ele_shape(grad_psi, ele)

       call compute_inverse_jacobian(ele_val(X,ele), ele_shape(X,ele), &
            detwei=detwei, invJ=invJ)

       ! Calculate the CFL number at each quadrature point.
       ! The matmul is the transpose of what I originally thought it should
       ! be. I don't understand why it's this way round but the results
       ! appear correct. -dham
       CFL_q=ele_val_at_quad(grad_psi, ele)

       do gi=1, size(detwei)
          cfl_q(:,gi)=cfl_q(:,gi)/(tol+sqrt(sum(cfl_q(:,gi)**2)))
       end do

       do gi=1, size(detwei)
          CFL_q(:,gi)=matmul(CFL_q(:,gi), invJ(:,:,gi))
       end do

       ! Project onto the basis functions to recover CFL at each node.
       CFL_mat(:,:)=matmul(inverse(shape_shape(CFL_shape, CFL_shape, detwei)), &
            shape_shape(CFL_shape, CFL_shape, &
            &             detwei*maxval(abs(CFL_q),1)))

       ! CFL is inherently discontinuous. In the case where a continuous
       ! mesh is provided for CFL, the following takes the safest option
       ! of taking the maximum value at a node.
       u_dot_dx%val(ele_CFL)=max(u_dot_dx%val(ele_CFL), sum(CFL_mat,2))

    end do
    
    call halo_update(u_dot_dx)

    call invert(u_dot_dx)

    dt=0.1*minval(u_dot_dx)

    call allmin(dt)

  end subroutine calculate_timestep
    
    subroutine pgrad(psi,X,grad_psi)
      type(scalar_field), intent(in)  :: psi
      type(vector_field), intent(in) :: X
      type(vector_field) ,intent(inout) :: grad_psi

      integer :: ele
      integer, dimension(:), pointer :: nodes
      integer :: i

      real, dimension(mesh_dim(psi),ele_loc(psi,1)) :: dgp
      real, dimension(mesh_dim(psi)) :: gp


      call zero(DiscontinuousGradient)
      call grad(psi,X,DiscontinuousGradient)
      

      grad_psi%val=0.0
!      grad_psi%val=huge(1.0)
      do ele=1,ele_count(psi)
         nodes=>ele_nodes(grad_psi,ele)
!         dgp=ele_val(DiscontinuousGradient,ele)*element_volume(positions,ele)#
         dgp=ele_val(DiscontinuousGradient,ele)!*element_volume(positions,ele)

!         call addto(grad_psi,nodes,dgp)

         do i=1, size(nodes)
            gp=node_val(grad_psi,nodes(i))
            if (abs(dot_product(gp,gp)-1)<abs(dot_product(dgp(:,i),dgp(:,i))-1)) &
                 call set(grad_psi,nodes(i),&
                 dgp(:,i))        
         end do
      end do
     
!      call scale(grad_psi,inverse_volume)

      call halo_update(grad_psi)

    end subroutine pgrad


      function L(psi,grad_psi,node_number)
        real :: L
        integer, intent(in) :: node_number
        type(scalar_field)  :: psi
        type(vector_field) :: grad_psi

        real :: gp(mesh_dim(psi))

        gp=node_val(grad_psi,node_number)
        
        L=-signum(psi,grad_psi,node_number)*(sqrt(dot_product(gp,gp))-1.0)
!        L=-signum(sfield,gradient(4),node_number)*(sqrt(dot_product(gp,gp))-1.0)

      end function L
        

      function signum(psi,grad_psi,node_number)
        real :: signum
        integer, intent(in) :: node_number
        type(scalar_field)  :: psi
        type(vector_field) :: grad_psi
        
        real :: p, gp(mesh_dim(psi)), dx

        p=node_val(psi,node_number)
        gp=node_val(grad_psi,node_number)
        dx=node_val(mesh_length,node_number)



!        signum=p/sqrt(p**2+dx**2*dot_product(gp,gp))
        signum=p/sqrt(p**2+dx**2)

      end function signum

  end subroutine hamilton_jacobi_distance_function

  subroutine init_marching(marching_data,sfield,first_time)
    type(marching_type), intent(inout) :: marching_data
    type(scalar_field), intent(in)     :: sfield
    logical, intent(in) :: first_time

    ! local variables

    integer :: n_ele, n_node

    ! routine initializes the data structures for the Fast Marching Method
    ! Believed to be independent of mesh structure


    if (first_time) then
       n_ele=ele_count(sfield)
       n_node=node_count(sfield)

       allocate(marching_data%ele%trial(n_ele)) 
       allocate(marching_data%ele%accepted(n_ele))
       allocate(marching_data%node%trial(n_node))
       allocate(marching_data%node%accepted(n_node))

    end if

    marching_data%ele%trial=.false.
    marching_data%ele%accepted=.false.
    marching_data%node%trial=.false.
    marching_data%node%accepted=.false.

  end subroutine init_marching

  subroutine finalize_marching(marching_data)
    type(marching_type), intent(inout) :: marching_data

    ! Clean-up memory post march

    deallocate(marching_data%ele%trial)
    deallocate(marching_data%ele%accepted)
    deallocate(marching_data%node%trial)
    deallocate(marching_data%node%accepted)
  end subroutine finalize_marching
  
  logical function contains_interface(sfield,ele)
    type(scalar_field), intent(inout) :: sfield
    integer :: ele

    ! Routine identifies elements which contain the zero contour
    
    if (abs(sum(merge(1,-1,ele_val(sfield,ele)>0.0)))<ele_loc(sfield,ele)&
         .or. any(abs(ele_val(sfield,ele))<1.0d-8)) then
       contains_interface=.true.
    else
       contains_interface=.false.
    end if
  end function contains_interface

 subroutine setup_marching(marching,sfield,positions,max_distance,interval)
   type(marching_type), intent(inout) :: marching
   type(scalar_field), intent(inout) :: sfield
   type(vector_field), intent(in) :: positions
   real, intent(in) :: max_distance
   real, intent(in), dimension(2) :: interval
   type(scalar_field) :: data 
   integer :: ele, node

   real :: val


   ! local varibles

   integer, dimension(:), pointer :: nodes

   ! Routine producing the initial zero contour values

   if (all(interval==[0.0,0.0])) then

      call allocate(data,sfield%mesh,"TemporaryStorage")
      
      do node=1, node_count(sfield)
         call set(data,node,signum(max_distance,node_val(sfield,node)))
      end do

      do ele= 1, ele_count(sfield)
         if (contains_interface(sfield,ele)) then
            marching%ele%accepted(ele)=.true.
            nodes=>ele_nodes(sfield,ele)
            marching%node%accepted(nodes)=.true.
            call estimate_nodes_first_go(marching,sfield,data,positions,ele)
         end if
      end do

      call set(sfield,data)
      call deallocate(data)

   else
      do node=1, node_count(sfield)
         val=node_val(sfield,node)
         if (val>=interval(1) .and. val<=interval(2)) then
            marching%node%accepted(node)=.true.
         end if
      end do
   
      do ele= 1, ele_count(sfield)
         nodes=>ele_nodes(sfield,ele)
         if (all(marching%node%accepted(nodes))) marching%ele%accepted(ele)=.true.
      end do
   end if


   do node=1, node_count(sfield)
      if (marching%node%accepted(node)) then
         call accept_node_and_recalculate_neighbours(node,&
                       marching,sfield,positions)
      end if
   end do

   contains

     function signum(a,b)

       real , intent(in) :: a,b

       real :: signum

       real, parameter :: tol=1.0e-16

       signum= a*b/sqrt(tol**2+b**2)

     end function signum

 end subroutine setup_marching

 subroutine accept_node_and_recalculate_neighbours(accepted_node,&
                       marching,sfield,positions)
   integer, intent(in) :: accepted_node
   type(marching_type), intent(inout) :: marching
   type(scalar_field), intent(inout) :: sfield
   type(vector_field), intent(in) :: positions

   integer :: node

!   local variables

   integer, dimension(:), pointer :: neighbour_list, nodes
   integer :: neigh,ele

   ! Routine updates the marching data structure to place accepted_node
   ! in the accepted space and place into the trial space any unaccepted
   ! neighbours not already in it.
   ! The distance to interface is then recalulated for trial space neighbours

   marching%node%accepted(accepted_node)=.true.
   marching%node%trial(accepted_node)=.false.

   neighbour_list=>node_neigh(sfield,accepted_node)

   do neigh = 1 , size(neighbour_list)
      ele=neighbour_list(neigh)
      if (ele/=0) then
         if (marching%ele%accepted(ele)) then
            cycle
         else
            nodes=>ele_nodes(sfield,ele)
            if (all(marching%node%accepted(nodes))) then
               marching%ele%accepted(ele)=.true.
               cycle
            end if
            marching%ele%trial(ele)=.true.
            do node=1,size(nodes) 
               marching%node%trial(nodes(node))=.true. .and. (.not. marching%node%accepted(nodes(node)))
            end do
            call estimate_nodes(marching,sfield,positions,ele)
         end if
      end if
   end do
 end subroutine accept_node_and_recalculate_neighbours


real function norm(n)
  real, dimension(:) ::n
  norm=sqrt(sum(n*n))
end function norm


function normed(n)
  real, dimension(:) :: n
  real, dimension(size(n)) :: normed
  normed = n/norm(n)
end function normed

real function distance_to_line(n,x)
  real, dimension(:) :: n,x

  
  distance_to_line=abs(sum(x*normed(n)))

end function distance_to_line

function cross(a,b)
  real, dimension(3), intent(in) :: a,b
  real, dimension(3) :: cross

  cross(1)=a(2)*b(3)-a(3)*b(2)
  cross(2)=a(3)*b(1)-a(1)*b(3)
  cross(3)=a(1)*b(2)-a(2)*b(1)
end function cross

function get_line(d,p,a,o)
  real :: d
  real, dimension(3) :: p,a,o

  real, dimension(3) :: get_line

  real, dimension(3) :: x,y,k

  real :: nx

  x=a-o
  k=p-o

  nx=norm(x)

  y=normed(nx**2*k -sum(x*k)*x)

  get_line=(1.0d0-d**2/nx**2)*x&
       +d*sqrt(1.0d0-d**2/nx**2)*y
end function get_line

 real function plane_distance(x,p,a)

   real, dimension(3) :: x
   real, dimension(3,3) :: p
   real, dimension(3) :: n1, n2

   real, dimension(2) :: a

   n1=get_line(a(1),x,p(:,1),p(:,3))
   n2=get_line(a(2),x,p(:,2),p(:,3))

   plane_distance=abs(sum(norm(cross(n1,n2))*(x-p(:,3))))
 end function plane_distance

real function line_distance(x,p,a)

real, dimension(:) :: x
real, dimension(:,:) :: p
real, dimension(size(x)) :: y,z,n
real :: a, n1,n2,n3

n1=sum((p(:,1)-p(:,2))*(p(:,1)-p(:,2)))
n2=sum((p(:,1)-p(:,2))*(x-p(:,2)))
n3=sum((x-p(:,2))*(x-p(:,2)))

if (abs(a)>=norm(x-p(:,2)) .or. abs(a)>=norm(p(:,1)-p(:,2))) then
   line_distance=huge(1.0d0)
else
   line_distance=abs((n1-a*a)*(n1*n3-n2*n2)&
        -a*sqrt(1.0-a*a/n1)*n2*sqrt(n1*(n1*n3-n2*n2)))&
        /sqrt(n1*n1*(n1-a*a)*(n1*n3-n2*n2))
end if

end function line_distance

real function point_distance(x,p)
  real, dimension(:) :: x,p
  point_distance=norm(x-p)
end function point_distance

function sort_index(i,a)
  integer, dimension(:) :: i
  real, dimension(:) :: a

  integer, dimension(size(i)) :: sort_index
  logical, dimension(size(i)) :: unsorted 
  integer :: j
  integer, dimension(1) ::k

  unsorted = .true.

  do j=1,size(i)
     k=minloc(a(pack(i,unsorted)))
     sort_index(j)=i(k(1))
     unsorted(k(1))=.false.
  end do
end function sort_index

subroutine estimate_nodes(marching,sfield,positions,ele,work_node,changed)
  type(marching_type), intent(inout) :: marching
  type(scalar_field), intent(inout) :: sfield
  type(vector_field), intent(in) :: positions
  integer, intent(in) ::ele
  integer, intent(in), optional :: work_node
  logical, intent(out), optional :: changed
 
  ! local variables

  integer i,node
  integer, dimension(1) :: vnode
  real, dimension(positions%dim) :: x
  integer nlocs, n_accepted

  integer, dimension(:), allocatable ::work_nodes

  integer, dimension(:), pointer :: node_list
  real, dimension(ele_loc(sfield,ele)) :: sval
  logical, dimension(ele_loc(sfield,ele)) :: accepted
  real, dimension(positions%dim,3) :: p
  real, dimension(3) :: a
  real :: d
  logical :: lchanged

  real, parameter :: epsilon=1.0d-8

  ! Calculate distances on the element to be tested, and update ones
  ! which get smaller

  nlocs=ele_loc(sfield,ele)
  node_list=>ele_nodes(sfield,ele)
  sval=ele_val(sfield,ele)


  accepted=marching%node%accepted(node_list)
  n_accepted=0
  lchanged=.false.

  accepted=marching%node%accepted(node_list)
  if (present(work_node)) then
     allocate(work_nodes(1))
     work_nodes(1)=work_node
  else
     allocate(work_nodes(nlocs-sum(merge(1,0,mask=accepted))))
     work_nodes=sort_index(pack((/(i,i=1,nlocs)/),mask=.not. accepted),abs(sval))
  end if

  do i=1,size(work_nodes)
     node=work_nodes(i)
     x=node_val(positions,node_list(node))
     d=abs(sval(node))
     
     accepted=(sval>=0.0 .eqv. sval(node)>=0.0) &
          .and. (abs(sval)<abs(sval(node)))
     n_accepted=0
     do 
        if (any(accepted)) then
           vnode= maxloc(abs(sval),mask=accepted)
           accepted(vnode)=.false.
           n_accepted=n_accepted+1
           p(:,n_accepted)=node_val(positions,node_list(vnode(1)))
           a(n_accepted)=abs(sval(vnode(1)))
        else
           exit
        end if
     end do

     if (n_accepted==3) then
 !       for some reason the distance from a zero plan isn't right yet.
 !       d=min(d,abs(a(1))+plane_distance(x,p,(/abs(a(1)-a(3)),&
 !            abs(a(2)-a(3))/)))
        d=min(d,abs(a(1))+line_distance(x,p(:,(/3,1/)),abs(a(1)-a(3))))
        d=min(d,abs(a(2))+line_distance(x,p(:,(/3,2/)),abs(a(2)-a(3))))
        d=min(d,abs(a(3))+point_distance(x,p(:,3)))
     end if
     if (n_accepted>=2) then
        d=min(d,abs(a(1))+line_distance(x,p(:,(/2,1/)),abs(a(1)-a(2))))
        d=min(d,abs(a(2))+point_distance(x,p(:,2)))
     end if
     d=min(d,abs(a(1))+point_distance(x,p(:,1)))
     if (d-abs(sval(node))<epsilon) then
        lchanged=.true.
        sval(node)=sign(d,sval(node))
        call set(sfield,node_list(node),sval(node))
     end  if
  end do

  if (present(changed)) changed=lchanged
  
end subroutine estimate_nodes


subroutine test_node_minimized(min_node,test_node,marching,sfield,positions)
integer, intent(in) :: min_node
logical, intent(inout) :: test_node
type(marching_type), intent(inout) :: marching
type(scalar_field), intent(inout) :: sfield
type(vector_field), intent(in) :: positions

! local variables

integer :: ele, neigh
logical :: changed
integer, dimension(:), pointer :: neighbour_list

test_node= .true.
neighbour_list=>node_neigh(sfield,min_node)

do neigh = 1 , size(neighbour_list)
   ele=neighbour_list(neigh)
   if (ele/=0) then
      if (marching%ele%accepted(ele)) then
         cycle
      else
         call estimate_nodes(marching,sfield,positions,ele,changed=changed)
         if (changed) test_node=.false.
      end if
   end if
end do


end subroutine test_node_minimized



 subroutine estimate_nodes_first_go(marching,sfield,data,positions,ele)
   type(marching_type), intent(inout) :: marching
   type(scalar_field), intent(in) :: sfield
   type(scalar_field), intent(inout) :: data
   type(vector_field), intent(in) :: positions
   integer, intent(in) :: ele

   ! local variables

   type(element_type), pointer :: ele_shape

   real, dimension(ele_loc(sfield,ele)) :: sval,dval
   integer :: nlocs, node,j
   real, dimension(positions%dim,positions%dim) :: Ainv
   real, dimension(positions%dim) :: n
   real, dimension(positions%dim,ele_loc(sfield,ele)) :: x
   
   integer, parameter:: LINE=2,TRIANGLE=13,TET=24
   integer, dimension(:), pointer :: node_list
   type(scalar_field) :: temp

   ! get distances on an element which contains the zero contour,
   ! and store any which are better than the previous estimate.

   nlocs=ele_loc(sfield,ele)
   sval=ele_val(sfield,ele)
   dval=ele_val(data,ele)
   x=ele_val(positions,ele)
   node_list=>ele_nodes(sfield,ele)

   if (.not. all(abs(sval)<1.0d-8)) then
      select case(10*(positions%dim-1)+nlocs)
      case(LINE)
!          solve a.x=h-c, then turn into a distance function,
!          d1=(h1-c)/a, d2=x-d1

!          =>  c=h1, a=(h2-h1)/x
         dval(:)=sval(:)*abs((x(1,2)-x(1,1))*(sval(:)/(sval(2)-sval(1))))
      case(TRIANGLE)

!          solve An=h-c, then turn into a distance function,
!          
!         A= ( x_2-x_1, y_2-y_1)
!            ( x_3-x_1, y_3-y_1)
!         h-c = (h2-h_1,h3-h1)
 
         Ainv(1,1)=x(2,3)-x(2,1)
         Ainv(1,2)=x(2,1)-x(2,2)
         Ainv(2,1)=x(1,1)-x(1,3)
         Ainv(2,2)=x(1,2)-x(1,1)

         Ainv=Ainv/(Ainv(1,1)*Ainv(2,2)-Ainv(1,2)*Ainv(2,1))

         n=matmul(Ainv,(/sval(2)-sval(1),sval(3)-sval(1)/))
         dval=sval/sqrt(sum(n*n))
      case(TET)

!          solve An=h-c, then turn into a distance function,
!          
!         A= ( x_1-x_4, y_1-y_4 z_1-z_4)
!            ( x_2-x_4, y_2-y_4 z_2-z_4)
!            ( x_3-x_4, y_3-y_4 z_3-z_4)


         do j=1,3
            x(:,j)=x(:,j)-x(:,4)
         end do

         Ainv(1,1)=(x(2,2)*x(3,3)-x(3,2)*x(2,3))
         Ainv(1,2)=(x(1,3)*x(3,2)-x(1,2)*x(3,3))
         Ainv(1,3)=(x(1,2)*x(2,3)-x(1,3)*x(2,2))
         Ainv(2,1)=(x(2,3)*x(3,1)-x(2,1)*x(3,3))
         Ainv(2,2)=(x(1,1)*x(3,3)-x(1,3)*x(3,1))
         Ainv(2,3)=(x(1,3)*x(2,1)-x(1,1)*x(2,3))
         Ainv(3,1)=(x(2,1)*x(3,2)-x(2,2)*x(3,1))
         Ainv(3,2)=(x(1,2)*x(3,1)-x(1,1)*x(3,2))
         Ainv(3,3)=(x(1,1)*x(2,2)-x(1,2)*x(2,1))
         
         
         Ainv=Ainv/(x(1,1)*Ainv(1,1)&
                   -x(1,2)*Ainv(1,2)&
                   +x(1,3)*Ainv(1,3))

         n=matmul(Ainv,(/sval(1)-sval(4),&
                         sval(2)-sval(4),&
                         sval(3)-sval(4)/))
         dval=sval/sqrt(sum(n*n))
!      case(QUAD)
!      case(HEX)
         case default

            FLAbort("Distance function method not yet implemented for element type")
      end select


      do node=1,ele_loc(sfield,ele)
         if (abs(dval(node))<abs(node_val(data,node_list(node)))&
              .and. sval(node)*dval(node)>=0.0 ) then
            call set(data,node_list(node),dval(node))
         end if
      end do

   end if
 end subroutine estimate_nodes_first_go

 subroutine get_interval(sfield,interval,max_distance)
   type(scalar_field), intent(inout) :: sfield
   real, intent(out) :: interval(2)
   real, intent(in) :: max_distance

   type(scalar_field) :: temp
   integer :: node
   real :: a,b

   interval=0.0

   call allocate(temp,sfield%mesh,"TemporaryField")
   call set(temp,sfield)
   call halo_update(sfield)

   ! loop over the nodes of the input and check for changes. This could be optimized

   do node=1,node_count(sfield)
      a=node_val(sfield,node)
      b=node_val(temp,node)
      if (a==b) then
         call set(temp,node,0.0)
      else
         call set(temp,node,sign(min(abs(a),abs(b)),a))
         call set(sfield,node,sign(min(abs(a),abs(b)),a))
      end if
   end do
   if(any(temp%val<0.0)) interval(1)=maxval(temp%val,mask=temp%val<0.0)
   if(any(temp%val>0.0)) interval(2)=minval(temp%val,mask=temp%val>0.0)

   call deallocate(temp)

 end subroutine get_interval

 subroutine marching_distance_function(sfield,positions,max_distance)
   type(scalar_field), intent(inout) :: sfield
   type(vector_field), intent(in) :: positions
   real, intent(in) :: max_distance

!  local data structure for trial space and accepted space
   type(marching_type) :: marching

   integer, dimension(1) :: min_node 
   logical :: test_node
   real, dimension(2) :: interval
   integer :: itr

   ! main routine

   call init_marching(marching,sfield,first_time=.true.)

   ! get zero contour
   call setup_marching(marching,sfield,positions,max_distance,[0.0,0.0])

   call run_marching()

   if(isparallel()) then

      do itr=2,50
         ! get interval does a halo update and thus spreads the information across processes a bit
         call get_interval(sfield,interval,max_distance)
         call init_marching(marching,sfield,first_time=.false.)
         call setup_marching(marching,sfield,positions,max_distance,interval)

         call run_marching()

      end do

   end if

   call finalize_marching(marching)

   contains 

     subroutine run_marching()

       if (any(marching%node%trial)) then
          ! after setup is there any work to do?
          
          do
             ! Main loop
             ! stop when all points are set to minimum value they can attain
             if (all (marching%node%accepted)) exit
             min_node=minloc(abs(sfield%val),mask=marching%node%trial)
             !      call test_node_minimized(min_node(1),test_node,marching,sfield,positions)
             !      if (test_node) then
             call accept_node_and_recalculate_neighbours(min_node(1),&
                  marching,sfield,positions)
             !      end if
          end do
       end if
     end subroutine run_marching

 end subroutine marching_distance_function

end module distance_function



