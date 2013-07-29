
#include "fdebug.h"
 module sigma_d0
  !This module only works in this case: 1)3-dimentional case, 2)z direction in the same direction as gravity, 3)mesh element is triangle.
    use fields
    use state_module
    use spud
    use sparse_tools
    use petsc_solve_state_module
    use sparse_tools_petsc
    use sparse_matrices_fields
    use field_options
    use fields_manipulation

    implicit none

    private
    public calculate_sigma_element, calculate_diagnostic_sigma_d0
    
 contains
 subroutine calculate_sigma_element(ele, x, u, sigma_ele, d0_a,dt,d0)
   integer, intent(in) :: ele
   type(vector_field), intent(in) :: u, x
   integer, parameter:: dim=3
   real,intent(in) :: d0_a, d0
   real, dimension(:,:), allocatable :: x_ele
   real, intent(inout):: sigma_ele
   integer, dimension(:), pointer:: u_ele
   real :: dz_ele,temp
   real, dimension(6)::dz
   real ::length1, length2, length3, length4, length5, length6,dx_ele
   real, intent(in)::dt
   integer::i,j
      allocate(x_ele(x%dim,ele_loc(x,ele)))     
      !get the nodes cordinate of element
      x_ele = ele_val(x, ele) 
      u_ele => ele_nodes(u, ele)
      
      !calculate the lateral lenghth of element's projection to horizontal surface
      length1 = sqrt((x_ele(2,2)-x_ele(2,1))**2+(x_ele(1,2)-x_ele(1,1))**2)
      length2 = sqrt((x_ele(2,3)-x_ele(2,1))**2+(x_ele(1,3)-x_ele(1,1))**2)
      length3 = sqrt((x_ele(2,2)-x_ele(2,3))**2+(x_ele(1,2)-x_ele(1,3))**2)
      length4 = sqrt((x_ele(2,4)-x_ele(2,1))**2+(x_ele(1,4)-x_ele(1,1))**2)
      length5 = sqrt((x_ele(2,4)-x_ele(2,2))**2+(x_ele(1,4)-x_ele(1,2))**2)
      length6 = sqrt((x_ele(2,4)-x_ele(2,3))**2+(x_ele(1,4)-x_ele(1,3))**2)
      dx_ele = max(length1,length2,length3,length4,length5,length6)
      dz(1)=abs(x_ele(3,1)-x_ele(3,2))
      dz(2)=abs(x_ele(3,1)-x_ele(3,3))
      dz(3)=abs(x_ele(3,1)-x_ele(3,4))
      dz(4)=abs(x_ele(3,2)-x_ele(3,3))
      dz(5)=abs(x_ele(3,2)-x_ele(3,4))
      dz(6)=abs(x_ele(3,3)-x_ele(3,4))
      do j=1,5
        do i=1,6-j
          if (dz(i)>dz(i+1)) then
             temp=dz(i+1)
	     dz(i+1)=dz(i)
	     dz(i)=temp  
	  end if     
	  end do 
      end do
      dz_ele = (dz(4)+dz(5)+dz(6))/3
      print *, 'dz=', dz
      print *, 'dz_ele', dz_ele
      do i=1, ele_loc(x,ele)
        if((945 < x_ele(1,i) .and. x_ele(1,i) < 1000) .and. (4945 < x_ele(2,i) .and. x_ele(2,i) < 5000)) then
          print *, 'dx_flat = ',dx_ele
          print *, 'dz_flat = ',dz_ele
        end if 
        if ((1945 < x_ele(1,i) .and. x_ele(1,i) < 2001) .and. (1945 < x_ele(2,i) .and. x_ele(2,i) < 2001)) then
          print *, 'dx_hump = ',dx_ele
          print *, 'dz_hump = ',dz_ele 
        end if
      end do
     ! sigma = dx_ele**2/(a**2*dt*dz_ele**2).
      sigma_ele = dx_ele**2/(d0_a**2*dt*dz_ele**2)
      deallocate(x_ele)    

   end subroutine calculate_sigma_element
 
 subroutine calculate_diagnostic_sigma_d0(state, sigma)
    type(state_type),intent(in) :: state
    integer :: ele
    type(scalar_field), intent(inout) :: sigma
    integer, dimension(ele_loc(sigma, 1)) ::ele_sigma
    real,dimension(ele_loc(sigma, 1)) ::sigma_node
    real,  dimension(:), allocatable ::depth_ele 
    type(vector_field) :: U, X
    integer, parameter:: dim=3
    real :: d0_a ,dt, d0
   real, dimension(:,:), allocatable :: x_ele
   real :: dz_ele, temp
   real, dimension(6)::dz
   integer::j
   real ::length1, length2, length3, length4, length5, length6, dx_ele
   real, dimension(ele_loc(sigma, 1),ele_loc(sigma, 1)) :: sigma_mat
   type(scalar_field), pointer :: dtt, dtb
   type(scalar_field) :: depth
   integer :: node,i
    U=extract_vector_field(state, "Velocity")
    X=extract_vector_field(state, "Coordinate")
    dtt => extract_scalar_field(state, "DistanceToTop")
    dtb => extract_scalar_field(state, "DistanceToBottom")
    call allocate(depth, dtt%mesh, "Depth")
    do node=1,node_count(dtt)
    call set(depth, node, node_val(dtt, node)+node_val(dtb, node))
    end do
    call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/a", d0_a)
    call get_option("/timestepping/timestep", dt)
    call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/d0", d0)
    
    call zero(sigma) 
    do ele=1,element_count(sigma)
       ele_sigma=ele_nodes(sigma,ele)
      allocate(x_ele(x%dim,ele_loc(x,ele)))
      allocate(depth_ele(ele_loc(depth,ele)))
      !get the nodes cordinate of element
      x_ele = ele_val(x, ele) 
      depth_ele = ele_val(depth, ele)
      !calculate the lateral lenghth of element's projection to horizontal surface
      length1 = sqrt((x_ele(2,2)-x_ele(2,1))**2+(x_ele(1,2)-x_ele(1,1))**2)
      length2 = sqrt((x_ele(2,3)-x_ele(2,1))**2+(x_ele(1,3)-x_ele(1,1))**2)
      length3 = sqrt((x_ele(2,2)-x_ele(2,3))**2+(x_ele(1,2)-x_ele(1,3))**2)
      length4 = sqrt((x_ele(2,4)-x_ele(2,1))**2+(x_ele(1,4)-x_ele(1,1))**2)
      length5 = sqrt((x_ele(2,4)-x_ele(2,2))**2+(x_ele(1,4)-x_ele(1,2))**2)
      length6 = sqrt((x_ele(2,4)-x_ele(2,3))**2+(x_ele(1,4)-x_ele(1,3))**2)
      dx_ele = max(length1,length2,length3,length4,length5,length6)
     
      dz(1)=abs(x_ele(3,1)-x_ele(3,2))
      dz(2)=abs(x_ele(3,1)-x_ele(3,3))
      dz(3)=abs(x_ele(3,1)-x_ele(3,4))
      dz(4)=abs(x_ele(3,2)-x_ele(3,3))
      dz(5)=abs(x_ele(3,2)-x_ele(3,4))
      dz(6)=abs(x_ele(3,3)-x_ele(3,4))
       do j=1,5
        do i=1,6-j
          if (dz(i)>dz(i+1)) then
             temp=dz(i+1)
	     dz(i+1)=dz(i)
	     dz(i)=temp  
	  end if     
	  end do 
      end do
     dz_ele = (dz(4)+dz(5)+dz(6))/3
     !sigma = dx_ele**2/(a**2*dt*dz_ele**2*dt).
      sigma_node =dx_ele**2/(d0_a**2*dt*dz_ele**2)
      call set(sigma,ele_sigma,sigma_node)
      deallocate(x_ele)
      deallocate(depth_ele) 
    end do 
    call deallocate(depth) 
 end subroutine calculate_diagnostic_sigma_d0
 
 end module sigma_d0
