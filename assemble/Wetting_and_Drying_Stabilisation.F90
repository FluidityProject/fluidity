
#include "fdebug.h"
 module wetting_and_drying_stabilisation
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
    use fields_base

    implicit none

    interface calculate_wetdry_vertical_absorption
       module procedure calculate_wetdry_vertical_absorption_element, calculate_wetdry_vertical_absorption_diagnostic
    end interface

    private
    public calculate_wetdry_vertical_absorption
    
contains

  subroutine calculate_wetdry_vertical_absorption_element(ele, x, u, sigma_ngi, optimum_aspect_ratio,dt,depth)
   integer, intent(in) :: ele
   type(vector_field), intent(in) :: u, x
   type(scalar_field),intent(in) :: depth
   integer, parameter:: dim=3
   real,intent(in) :: optimum_aspect_ratio
   real, dimension(:,:), allocatable :: x_ele
   real, dimension(ele_ngi(u,ele)), intent(inout):: sigma_ngi
   integer, dimension(:), pointer:: u_ele
   real ::length1, length2, length3,length4,length5,length6,dx_ele
   real, intent(in)::dt
   real, dimension(ele_ngi(u,ele)) :: depth_at_quads
   integer::i
   ewrite(6,*) 'Calculating wetting and drying vertical absorption stabilisation'
      allocate(x_ele(x%dim,ele_loc(u,ele)))
      !get the nodes cordinate of element
      x_ele = ele_val(x, ele) 
      u_ele => ele_nodes(u, ele)
      depth_at_quads=ele_val_at_quad(depth, ele)
      !calculate the lateral lenghth of element's projection to horizontal surface
      length1 = sqrt((x_ele(2,2)-x_ele(2,1))**2+(x_ele(1,2)-x_ele(1,1))**2)
      length2 = sqrt((x_ele(2,3)-x_ele(2,1))**2+(x_ele(1,3)-x_ele(1,1))**2)
      length3 = sqrt((x_ele(2,2)-x_ele(2,3))**2+(x_ele(1,2)-x_ele(1,3))**2)
      length4 = sqrt((x_ele(2,4)-x_ele(2,1))**2+(x_ele(1,4)-x_ele(1,1))**2)
      length5 = sqrt((x_ele(2,4)-x_ele(2,2))**2+(x_ele(1,4)-x_ele(1,2))**2)
      length6 = sqrt((x_ele(2,4)-x_ele(2,3))**2+(x_ele(1,4)-x_ele(1,3))**2)
      dx_ele = max(length1,length2,length3,length4,length5,length6)
      ! sigma = dx_ele**2/(a**2*dt*dz_ele**2).


      dx_ele = max(length1,length2,length3)
      dz_ele = max(x_ele(3,1),x_ele(3,2), x_ele(3,3) )-min(x_ele(3,1), x_ele(3,2), x_ele(3,3))
      ! sigma = dx_ele**2/(a**2*dt*dz_ele**2).But when we add it to mass and other terms, we need to multiply dt, so here we get the result after multiplying dt.
      sigma_ele = dx_ele**2/(d0_a**2*dt*dz_ele**2)

      !do i=1, ele_ngi(u,ele)
      !   sigma_ngi(i) = dx_ele**2/(optimum_aspect_ratio**2*dt*depth_at_quads(i)**2)
      !end do
      deallocate(x_ele)    

   end subroutine calculate_wetdry_vertical_absorption_element
 
  subroutine calculate_wetdry_vertical_absorption_diagnostic(state, sigma)
    type(state_type),intent(in) :: state
    integer :: ele
    type(scalar_field), intent(inout) :: sigma
   
    type(vector_field) :: U, X
    real, dimension(:),allocatable:: sigma_ele
    integer, parameter:: dim=3
    real :: optimum_aspect_ratio,dt
   real, dimension(:,:), allocatable :: x_ele
   real ::length1, length2, length3,length4,length5,length6, dx_ele
   type(scalar_field), pointer :: dtt, dtb
   type(scalar_field) :: depth
   real,dimension(:),allocatable::depth_ele
   integer,dimension(:),pointer::sigma_nodes
   integer :: node,i 

    call zero(sigma)  
    U = extract_vector_field(state, "Velocity")
    X = extract_vector_field(state, "Coordinate")
    dtt => extract_scalar_field(state, "DistanceToTop")
    dtb => extract_scalar_field(state, "DistanceToBottom")
    call allocate(depth, dtt%mesh, "Depth")
    do node=1,node_count(dtt)
       call set(depth, node, node_val(dtt, node)+node_val(dtb, node))
    end do
    call get_option("/timestepping/timestep", dt)
    call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/optimum_aspect_ratio", optimum_aspect_ratio)

    do ele=1,element_count(sigma)
      allocate(x_ele(x%dim,ele_loc(x,ele)))
      allocate(sigma_ele(ele_loc(sigma,ele)))
      allocate(depth_ele(ele_loc(x,ele)))
      depth_ele=ele_val(depth, ele)
      sigma_nodes => ele_nodes(sigma,ele)
      !get the nodes cordinate of element
      x_ele = ele_val(x, ele)
      !calculate the lateral lenghth of element's projection to horizontal surface
      length1 = sqrt((x_ele(2,2)-x_ele(2,1))**2+(x_ele(1,2)-x_ele(1,1))**2)
      length2 = sqrt((x_ele(2,3)-x_ele(2,1))**2+(x_ele(1,3)-x_ele(1,1))**2)
      length3 = sqrt((x_ele(2,2)-x_ele(2,3))**2+(x_ele(1,2)-x_ele(1,3))**2)
      length4 = sqrt((x_ele(2,4)-x_ele(2,1))**2+(x_ele(1,4)-x_ele(1,1))**2)
      length5 = sqrt((x_ele(2,4)-x_ele(2,2))**2+(x_ele(1,4)-x_ele(1,2))**2)
      length6 = sqrt((x_ele(2,4)-x_ele(2,3))**2+(x_ele(1,4)-x_ele(1,3))**2)
      
      !dx_ele = max(length1,length2,length3,length4,length5,length6)
      
      !sigma = dx_ele**2/(a**2*dt*dz_ele**2).
      !do i=1, size(sigma_nodes)
      !  sigma_ele(i) = dx_ele**2/ (optimum_aspect_ratio**2 * dt * maxval(depth_ele) **2) 
      !end do
      !call set(sigma,sigma_nodes,sigma_ele)

      dx_ele = max(length1,length2,length3)
      dz_ele = max(x_ele(3,1),x_ele(3,2), x_ele(3,3) )-min(x_ele(3,1), x_ele(3,2), x_ele(3,3))
      call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/a", d0_a)
      !sigma = dx_ele**2/(a**2*dt*dz_ele**2).But when we add it to mass and other terms, we need to multiply dt, so here we get the result after multiplying dt.
      sigma_node =dx_ele**2/(d0_a**2*dz_ele**2)
      call set(sigma,ele_sigma,sigma_node)
      
      deallocate(x_ele)
      deallocate(sigma_ele) 
      deallocate(depth_ele) 
    end do 
  end subroutine calculate_wetdry_vertical_absorption_diagnostic
 
 end module wetting_and_drying_stabilisation
