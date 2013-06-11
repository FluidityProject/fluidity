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
       module procedure calculate_wetdry_vertical_absorption_element_tensor, &
        & calculate_wetdry_vertical_absorption_diagnostic
    end interface

    private
    public calculate_wetdry_vertical_absorption
    
contains

  subroutine calculate_wetdry_vertical_absorption_element(absorption, ele, x, gravity, dt, optimum_aspect_ratio)
    ! Calculates the vertical absoption for each quadrature point
    real, intent(out), dimension(:) :: absorption
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: x
    type(vector_field), intent(in) :: gravity
    real, intent(in) :: dt
    real, intent(in) :: optimum_aspect_ratio
    
    real :: dr, dh
    integer :: node
    integer, dimension(:), pointer :: enodes
    real, dimension(x%dim) :: pos, gravity_at_node
    real, dimension(ele_loc(x,ele)) :: r
    real, dimension(ele_loc(x,ele) - 1) :: characteristic_horiz_scales
    real, parameter :: tolerance = 1.0e-5
    
    ewrite(6,*) 'Calculating wetting and drying optimum aspect ratio vertical absorption, optimum aspect ratio =', optimum_aspect_ratio

    enodes => ele_nodes(X, ele)
    do node = 1, size(enodes)
      ! Calculate vertical extent, parallel to the gravitational vector
      pos = node_val(X, enodes(node))
      gravity_at_node = node_val(gravity, enodes(node))
      r(node) = dot_product(pos, gravity_at_node)
      ! Calculate the horizontal extent
      ! Assumes gravitational direction is normalised (it should be)
      if (node .eq. 1) cycle
      characteristic_horiz_scales(node-1) = norm2( cross_product( &
        & node_val(X, enodes(node)) - node_val(X, enodes(1)), &
        & node_val(X, enodes(node)) - node_val(X, enodes(1)) - node_val(gravity, 1)) ) 
    end do
    ! Alternative calculation of horizontal extent
    characteristic_horiz_scales = calculate_lengths(ele, x)
    ! TODO: Look at using Jacobian to derive element chacteristic lengths

    ! Verical lengthscale, in a direction parallel to the gravitational acceleration
    ! Maximum extent
    dr = maxval(r) - minval(r)

    ! Horizontal length, in a plane perpendicular to gravitational acceleration
    ! Maximum horizontal edge length
    !dh = maxval(characteristic_horiz_scales)
    ! Minimum horizontal edge length
    dh = minval(characteristic_horiz_scales, mask = characteristic_horiz_scales .gt. tolerance)
    ! Average horizontal edge length
    !dh = sum(characteristic_horiz_scales) / size(characteristic_horiz_scales)

    absorption = dh**2 / ((optimum_aspect_ratio * dr)**2 * dt)

  end subroutine calculate_wetdry_vertical_absorption_element

  subroutine calculate_wetdry_vertical_absorption_element_tensor(absorption, ele, x, gravity, dt, optimum_aspect_ratio)
    ! Calculates the vertical absoption to add to the momentum equation
    real, dimension(:,:), intent(out) :: absorption
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: x
    type(vector_field), intent(in) :: gravity
    real, intent(in) :: dt
    real,intent(in) :: optimum_aspect_ratio
    real, dimension(size(absorption, 2)) :: absorption_ngi
    real, dimension(x%dim, size(absorption, 2)) :: grav_at_quads
    integer :: ngi

    absorption = 0.0
    grav_at_quads = ele_val_at_quad(gravity, ele)
    call calculate_wetdry_vertical_absorption_element(absorption_ngi, ele, X, gravity, dt, optimum_aspect_ratio)
    
    do ngi = 1, size(absorption, 2)
      absorption(:, ngi) = absorption_ngi(ngi) * grav_at_quads(:, ngi)
    end do

  end subroutine calculate_wetdry_vertical_absorption_element_tensor
 
  subroutine calculate_wetdry_vertical_absorption_diagnostic(state, absorption)
    type(state_type),intent(in) :: state
    type(scalar_field), intent(inout) :: absorption

    type(vector_field), pointer :: x
    type(vector_field), pointer :: gravity
    logical :: have_optimum_aspect_ratio
    real :: optimum_aspect_ratio
    ! Assume all elements have the same quadrature for the absorption field (currently the case)
    real, dimension(ele_ngi(absorption,1)) :: absorption_ngi
    integer :: ele
    real :: dt

    have_optimum_aspect_ratio = have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/optimum_aspect_ratio")
    if (.not. have_optimum_aspect_ratio) FLExit("Please provide an optimum_aspect_ratio.")
    call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/optimum_aspect_ratio", optimum_aspect_ratio)
    call get_option("/timestepping/timestep", dt)

    x => extract_vector_field(state, "Coordinate")
    gravity => extract_vector_field(state, "GravityDirection")
    call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/optimum_aspect_ratio", optimum_aspect_ratio)

    call zero(absorption)  
    do ele=1, element_count(absorption)
      call calculate_wetdry_vertical_absorption_element(absorption_ngi, ele, x, gravity, dt, optimum_aspect_ratio)
      ! TODO: Fix absorption_ngi -> absorption_loc here
      call set(absorption, ele_nodes(absorption, ele), absorption_ngi(1))
    end do 

  end subroutine calculate_wetdry_vertical_absorption_diagnostic

  function calculate_lengths(ele, x) result(characteristic_length)
    real, dimension(3) :: characteristic_length
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: x
    real, dimension(x%dim, ele_loc(x,ele)) :: ele_pos
    ele_pos = ele_val(x, ele) 
    characteristic_length(1) = sqrt((ele_pos(2,1)-ele_pos(2,2))**2 + (ele_pos(1,1)-ele_pos(1,2))**2)
    characteristic_length(3) = sqrt((ele_pos(2,2)-ele_pos(2,3))**2 + (ele_pos(1,2)-ele_pos(1,3))**2)
    characteristic_length(2) = sqrt((ele_pos(2,3)-ele_pos(2,1))**2 + (ele_pos(1,3)-ele_pos(1,1))**2)
    !characteristic_length(4) = sqrt((ele_pos(2,1)-ele_pos(2,4))**2 + (ele_pos(1,1)-ele_pos(1,4))**2)
    !characteristic_length(5) = sqrt((ele_pos(2,2)-ele_pos(2,4))**2 + (ele_pos(1,2)-ele_pos(1,4))**2)
    !characteristic_length(6) = sqrt((ele_pos(2,3)-ele_pos(2,4))**2 + (ele_pos(1,3)-ele_pos(1,4))**2)
  end function calculate_lengths
 
end module wetting_and_drying_stabilisation
