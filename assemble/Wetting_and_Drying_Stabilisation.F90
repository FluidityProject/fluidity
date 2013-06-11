#include "fdebug.h"
 module wetting_and_drying_stabilisation
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

  subroutine calculate_wetdry_vertical_absorption_element(absorption_ngi, ele, x, dt, optimum_aspect_ratio)
    ! absorption = dx_ele**2/(a**2*dt*dz_ele**2).
    !This module only works in this case: 1)3-dimentional case, 2)z direction in the same direction as gravity, 3)mesh element is triangle.

    real, dimension(:), intent(inout):: absorption_ngi
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: x
    real, intent(in)::dt
    real,intent(in) :: optimum_aspect_ratio
    
    real, dimension(:,:), allocatable :: x_ele
    real :: length1, length2, length3,length4,length5,length6,dx_ele
    real :: dz_ele
    integer::i
    ewrite(6,*) 'Calculating wetting and drying optimum aspect ratio vertical absorption, optimum aspect ratio =', optimum_aspect_ratio

    allocate(x_ele(x%dim,ele_loc(x,ele)))
    !get the nodes cordinate of element
    x_ele = ele_val(x, ele) 
    !calculate the lateral lenghth of element's projection to horizontal surface
    length1 = sqrt((x_ele(2,2)-x_ele(2,1))**2+(x_ele(1,2)-x_ele(1,1))**2)
    length2 = sqrt((x_ele(2,3)-x_ele(2,1))**2+(x_ele(1,3)-x_ele(1,1))**2)
    length3 = sqrt((x_ele(2,2)-x_ele(2,3))**2+(x_ele(1,2)-x_ele(1,3))**2)
    length4 = sqrt((x_ele(2,4)-x_ele(2,1))**2+(x_ele(1,4)-x_ele(1,1))**2)
    length5 = sqrt((x_ele(2,4)-x_ele(2,2))**2+(x_ele(1,4)-x_ele(1,2))**2)
    length6 = sqrt((x_ele(2,4)-x_ele(2,3))**2+(x_ele(1,4)-x_ele(1,3))**2)
    dx_ele = max(length1,length2,length3,length4,length5,length6)


    dx_ele = max(length1,length2,length3)
    dz_ele = max(x_ele(3,1),x_ele(3,2), x_ele(3,3) )-min(x_ele(3,1), x_ele(3,2), x_ele(3,3))
    absorption_ngi = dx_ele**2/(optimum_aspect_ratio**2*dt*dz_ele**2)

    deallocate(x_ele)    

   end subroutine calculate_wetdry_vertical_absorption_element
 
  subroutine calculate_wetdry_vertical_absorption_diagnostic(state, absorption)
    type(state_type),intent(in) :: state
    type(scalar_field), intent(inout) :: absorption

    type(vector_field) :: x
    logical :: have_optimum_aspect_ratio
    real :: optimum_aspect_ratio
    ! Assume all elements have the same quadrature for the absorption field (currently the case)
    real, dimension(ele_ngi(absorption,1)) :: absorption_ngi
    integer :: ele
    real :: dt

    have_optimum_aspect_ratio = have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/optimum_aspect_ratio")
    if(.not. have_optimum_aspect_ratio) FLExit("Please provide an optimum_aspect_ratio.")
    call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/optimum_aspect_ratio", optimum_aspect_ratio)
    call get_option("/timestepping/timestep", dt)

    x = extract_vector_field(state, "Coordinate")
    call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/optimum_aspect_ratio", optimum_aspect_ratio)

    call zero(absorption)  
    do ele=1, element_count(absorption)
      call calculate_wetdry_vertical_absorption_element(absorption_ngi, ele, x, dt, optimum_aspect_ratio)
      ! TODO: Fix absorption_ngi -> absorption_loc here
      call set(absorption, ele_nodes(absorption, ele), absorption_ngi(1))
    end do 

  end subroutine calculate_wetdry_vertical_absorption_diagnostic
 
end module wetting_and_drying_stabilisation
