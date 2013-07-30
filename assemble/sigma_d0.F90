
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
 subroutine calculate_sigma_element(ele, x, u, sigma_ngi, d0_a,dt,depth)
   integer, intent(in) :: ele
   type(vector_field), intent(in) :: u, x
   type(scalar_field),intent(in) :: depth
   integer, parameter:: dim=3
   real,intent(in) :: d0_a
   real, dimension(:,:), allocatable :: x_ele
   real,dimension(ele_ngi(u,ele)), intent(inout):: sigma_ngi
   real ::length1, length2, length3, length4, length5, length6,dx_ele
   real, intent(in)::dt
   real, dimension(ele_ngi(u,ele)) :: depth_at_quads
   integer::i
      allocate(x_ele(x%dim,ele_loc(x,ele)))     
      !get the nodes cordinate of element
      x_ele = ele_val(x, ele) 
      depth_at_quads = ele_val_at_quad(depth, ele)
      !print *, 'depth_at_quads',depth_at_quads
      !print *, 'dt',dt
     ! print *, 'a',d0_a
      !calculate the lateral lenghth of element's projection to horizontal surface
      length1 = sqrt((x_ele(2,2)-x_ele(2,1))**2+(x_ele(1,2)-x_ele(1,1))**2)
      length2 = sqrt((x_ele(2,3)-x_ele(2,1))**2+(x_ele(1,3)-x_ele(1,1))**2)
      length3 = sqrt((x_ele(2,2)-x_ele(2,3))**2+(x_ele(1,2)-x_ele(1,3))**2)
      length4 = sqrt((x_ele(2,4)-x_ele(2,1))**2+(x_ele(1,4)-x_ele(1,1))**2)
      length5 = sqrt((x_ele(2,4)-x_ele(2,2))**2+(x_ele(1,4)-x_ele(1,2))**2)
      length6 = sqrt((x_ele(2,4)-x_ele(2,3))**2+(x_ele(1,4)-x_ele(1,3))**2)
      dx_ele = max(length1,length2,length3,length4,length5,length6)
     ! print *, 'dx_ele',dx_ele
     ! sigma = dx_ele**2/(a**2*dt*depth_ele**2).
     do i=1, ele_ngi(u,ele)
        sigma_ngi(i) = dx_ele**2/(d0_a**2*dt*depth_at_quads(i)**2)
       ! print *,'sigma_ngi(i)',sigma_ngi(i)
     end do
      deallocate(x_ele)    

   end subroutine calculate_sigma_element
 
 subroutine calculate_diagnostic_sigma_d0(state, sigma)
    type(state_type),intent(in) :: state
    integer :: ele
    type(scalar_field), intent(inout) :: sigma
    integer, dimension(ele_loc(sigma, 1)) ::ele_sigma
    real, dimension(:),allocatable ::sigma_node
    real,  dimension(:), allocatable ::depth_ele 
    type(vector_field) :: U, X
    integer, parameter:: dim=3
    real :: d0_a ,dt
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
    
    call zero(sigma) 
    do ele=1,element_count(sigma)
       ele_sigma=ele_nodes(sigma,ele)
      allocate(x_ele(x%dim,ele_loc(x,ele)))
      allocate(sigma_node(ele_loc(u,ele)))
      allocate(depth_ele(ele_loc(u,ele)))
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
     !sigma = dx_ele**2/(a**2*dt*dz_ele**2*dt).
      do i=1, ele_loc(u,ele)
        sigma_node(i) = dx_ele**2/(d0_a**2*dt*depth_ele(i)**2)
     end do
      call set(sigma,ele_sigma,sigma_node)
      deallocate(x_ele)
      deallocate(sigma_node)
      deallocate(depth_ele) 
    end do 
    call deallocate(depth) 
 end subroutine calculate_diagnostic_sigma_d0
 
 end module sigma_d0
