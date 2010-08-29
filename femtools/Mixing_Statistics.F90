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

module mixing_statistics

  use elements
  use embed_python
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN
  use fields
  use field_derivatives
  use field_options
  use state_module
  use futils
  use fetools
  use fefields
  use halos
  use MeshDiagnostics
  use spud

  implicit none
  
  interface heaviside_integral
    module procedure heaviside_integral_single, &
      & heaviside_integral_multiple, heaviside_integral_options
  end interface heaviside_integral

  interface mixing_stats
     module procedure mixing_stats_scalar
  end interface

  private
  
  public :: heaviside_integral, mixing_stats

contains

  subroutine mixing_stats_scalar(f_mix_fraction, sfield, Xfield, mixing_stats_count)

    real, dimension(:), intent(out) :: f_mix_fraction
    type(scalar_field), target, intent(inout) :: sfield
    type(vector_field), target, intent(inout) :: Xfield
    integer, intent(in) :: mixing_stats_count

    if(have_option(trim(complete_field_path(sfield%option_path))//&
         &"/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/continuous_galerkin")) then

       call heaviside_integral(f_mix_fraction, sfield, Xfield, mixing_stats_count = mixing_stats_count)

    else if(have_option(trim(complete_field_path(sfield%option_path))//&
         &"/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/control_volumes")) then

       call control_volume_mixing(f_mix_fraction, sfield, Xfield, mixing_stats_count)

    else

       FLAbort("Only continuous_galerkin or control_volumes options available for mixing_bins.")

    end if

  end subroutine mixing_stats_scalar

  subroutine control_volume_mixing(f_mix_fraction, sfield, Xfield, mixing_stats_count)

    real, dimension(:), intent(out) :: f_mix_fraction
    type(scalar_field), target, intent(inout) :: sfield
    type(vector_field), target, intent(inout) :: Xfield
    integer, intent(in) :: mixing_stats_count

    integer :: i, j
    real :: tolerance, total_volume, current_time
    real, dimension(:), pointer :: mixing_bin_bounds
    character(len = OPTION_PATH_LEN) :: func
    type(scalar_field) :: lumped_mass
    integer, dimension(2) :: shape_option

    call get_option(trim(complete_field_path(sfield%option_path)) &
         &// "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/tolerance",&
         & tolerance, default = epsilon(0.0))
         
    if(have_option(trim(complete_field_path(sfield%option_path)) // &
        & "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/mixing_bin_bounds/constant")) then
       shape_option=option_shape(trim(complete_field_path(sfield%option_path)) // &
            & "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/mixing_bin_bounds/constant")
       allocate(mixing_bin_bounds(shape_option(1)))

       call get_option(trim(complete_field_path(sfield%option_path)) &
           &// "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/mixing_bin_bounds/constant",&
           & mixing_bin_bounds)

    else if(have_option(trim(complete_field_path(sfield%option_path)) // &
        & "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/mixing_bin_bounds/python")) then
        call get_option(trim(complete_field_path(sfield%option_path)) // &
          & "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/mixing_bin_bounds/python", func)
        call get_option("/timestepping/current_time", current_time)
        call real_vector_from_python(func, current_time, mixing_bin_bounds)
    else
        FLAbort("Unable to determine mixing bin bounds type")  
    end if   

    call allocate(lumped_mass, sfield%mesh, name="MixingLumpedMass")
    call zero(lumped_mass)
    call compute_lumped_mass(Xfield, lumped_mass)

    f_mix_fraction = 0.0
    total_volume = 0.0
    do j = 1, node_count(sfield)
       if(node_owned(sfield, j)) then
          do i=1,size(f_mix_fraction)
             if(node_val(sfield, j) >= (mixing_bin_bounds(i)-tolerance)) then
                f_mix_fraction(i) = f_mix_fraction(i) + node_val(lumped_mass, j)
             end if
          end do
          total_volume = total_volume + node_val(lumped_mass, j)
       end if
    end do

    ! In f_mix_fraction(j) have volume of sfield >= mixing_bin_bounds(j)
    ! Subtract and divide so that f_mix_fraction(j) has 
    ! mixing_bin_bounds(j)<volume fraction of sfield<mixing_bins_bounds(j+1) 

    call allsum(total_volume)

    if(have_option(trim(complete_field_path(sfield%option_path))//&
         &"/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/control_volumes/normalise")) then
       do i = 1,size(f_mix_fraction)-1
          f_mix_fraction(i) = (f_mix_fraction(i)-f_mix_fraction(i+1))/total_volume
       end do
       f_mix_fraction(size(f_mix_fraction)) = f_mix_fraction(size(f_mix_fraction))/total_volume
    else
       do i = 1,size(f_mix_fraction)-1
          f_mix_fraction(i) = (f_mix_fraction(i)-f_mix_fraction(i+1))
       end do
       f_mix_fraction(size(f_mix_fraction)) = f_mix_fraction(size(f_mix_fraction))
    end if

    do i = 1, size(f_mix_fraction)
       call allsum(f_mix_fraction(i))
    end do

    call deallocate(lumped_mass)
    deallocate(mixing_bin_bounds)

  end subroutine control_volume_mixing

  subroutine heaviside_integral_options(f_mix_fraction, sfield, Xfield, mixing_stats_count)

    real, dimension(:), intent(out) :: f_mix_fraction
    type(scalar_field), target, intent(inout) :: sfield
    type(vector_field), target, intent(in) :: Xfield
    integer, intent(in) :: mixing_stats_count

    integer :: j, ele
    real :: total_volume
    real, dimension(:), allocatable :: detwei
    real :: tolerance, current_time
    real, dimension(:), pointer :: mixing_bin_bounds
    character(len = OPTION_PATH_LEN) :: func      
    integer, dimension(2) :: shape_option

    if(have_option(trim(complete_field_path(sfield%option_path)) // &
        & "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/mixing_bin_bounds/constant")) then

       shape_option=option_shape(trim(complete_field_path(sfield%option_path)) // &
            & "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/mixing_bin_bounds/constant")
       allocate(mixing_bin_bounds(shape_option(1)))

        call get_option(trim(complete_field_path(sfield%option_path)) &
           &// "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/mixing_bin_bounds/constant",&
           & mixing_bin_bounds)
    else if(have_option(trim(complete_field_path(sfield%option_path)) // &
        & "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/mixing_bin_bounds/python")) then
        call get_option(trim(complete_field_path(sfield%option_path)) // &
          & "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/mixing_bin_bounds/python", func)
        call get_option("/timestepping/current_time", current_time)
        call real_vector_from_python(func, current_time, mixing_bin_bounds)
    else
        FLAbort("Unable to determine mixing bin bounds type") 
    end if

    call get_option(trim(complete_field_path(sfield%option_path)) &
         &// "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/tolerance",&
         & tolerance, default = epsilon(0.0))

    f_mix_fraction = heaviside_integral_multiple(sfield, mixing_bin_bounds, Xfield, tolerance = tolerance)
    
    do j = 1,(size(f_mix_fraction)-1)
       f_mix_fraction(j) = (f_mix_fraction(j)-f_mix_fraction(j+1))
    end do
            
    ! In f_mix_fraction(j) have volume of sfield >mixing_bin_bounds(j)
    ! Subtract so that f_mix_fraction(j) has 
    ! mixing_bin_bounds(j)<volume fraction of sfield<mixing_bins_bounds(j+1) 
    ! and if normalise is selected divide by total volume
    if(have_option(trim(complete_field_path(sfield%option_path))//&
         &"/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/continuous_galerkin/normalise")) then

       assert(ele_count(sfield) > 0)

       total_volume = 0.0

       allocate(detwei(ele_ngi(sfield, 1)))
       do ele = 1, ele_count(sfield)
          if(.not. element_owned(sfield,ele)) cycle

          call transform_to_physical(Xfield, ele, detwei = detwei)
          total_volume = total_volume + sum(detwei)
       end do
       deallocate(detwei)
       
       call allsum(total_volume, communicator = halo_communicator(sfield))
       
       f_mix_fraction = f_mix_fraction / total_volume
    endif
    
    deallocate(mixing_bin_bounds)
    
  end subroutine heaviside_integral_options
  
  function heaviside_integral_single(sfield, bound, positions, tolerance) result(integral)
    type(scalar_field), intent(inout) :: sfield
    real, intent(in) :: bound
    type(vector_field), intent(in) :: positions
    real, optional, intent(in) :: tolerance
    
    real :: integral
    
    real, dimension(1) :: integrals
    type(element_type), pointer :: shape
    
    shape => ele_shape(sfield, 1)
    if(ele_numbering_family(shape) /= FAMILY_SIMPLEX .or. shape%degree /= 1) then
      FLAbort("heaviside_integral requires a linear simplex input mesh")
    end if
    
    select case(positions%dim)
      case(3)
        integrals = heaviside_integral_tet(sfield, (/bound/), positions, tolerance = tolerance)
        integral = integrals(1)
      case(2)
        integrals = heaviside_integral_tri(sfield, (/bound/), positions)
        integral = integrals(1)
      case(1)
        integrals = heaviside_integral_line(sfield, (/bound/), positions)
        integral = integrals(1)
      case default
        FLExit("heaviside_integral requires a linear simplex input mesh")
    end select
    
  end function heaviside_integral_single
  
  function heaviside_integral_multiple(sfield, bounds, positions, tolerance) result(integrals)
    type(scalar_field), intent(inout) :: sfield
    real, dimension(:), intent(in) :: bounds
    type(vector_field), intent(in) :: positions
    real, optional, intent(in) :: tolerance
    
    real, dimension(size(bounds)) :: integrals
    
    type(element_type), pointer :: shape
    
    shape => ele_shape(sfield, 1)
    if(ele_numbering_family(shape) /= FAMILY_SIMPLEX .or. shape%degree /= 1) then
      FLExit("continuous_galerkin mixing_stats only available on linear simplex meshes")
    end if
    
    select case(positions%dim)
      case(3)
        integrals = heaviside_integral_tet(sfield, bounds, positions, tolerance = tolerance)
      case(2)
        integrals = heaviside_integral_tri(sfield, bounds, positions)
      case(1)
        integrals = heaviside_integral_line(sfield, bounds, positions)
      case default
        FLAbort("Unsupported dimension count in heaviside_integral_new")
    end select
    
  end function heaviside_integral_multiple

  function heaviside_integral_line(sfield, bounds, positions) result(integrals)
    type(scalar_field), intent(in) :: sfield
    real, dimension(:), intent(in) :: bounds
    type(vector_field), intent(in) :: positions
    
    real, dimension(size(bounds)) :: integrals
    
    integer :: ele, i
    logical, dimension(2) :: node_integrate
    logical, dimension(size(bounds)) :: integrate_length
    real :: line_length, sub_coord
    real, dimension(2) :: node_vals
    real, dimension(1, 2) :: coords
    
    assert(sfield%mesh == positions%mesh)
    
    integrals = 0.0
    ele_loop: do ele = 1, ele_count(sfield)
      if(.not. element_owned(sfield, ele)) cycle ele_loop
    
      assert(ele_loc(sfield, ele) == 2)
      coords = ele_val(positions, ele)
      node_vals = ele_val(sfield, ele)
      
      integrate_length = .false.
      bounds_loop: do i = 1, size(bounds)      
        node_integrate = node_vals >= bounds(i)
      
        if(node_integrate(1)) then
          if(node_integrate(2)) then
            ! Integrate the whole element
            
            integrate_length(i) = .true.  ! Handled outside the bounds loop (below)
          else
            ! Integrate part of the element
            sub_coord = ((bounds(i) - node_vals(1)) * (coords(1, 2) - coords(1, 1)) / (node_vals(2) - node_vals(1))) + coords(1, 1)
            integrals(i) = integrals(i) + abs(coords(1, 1) - sub_coord)
          end if
        else if(node_integrate(2)) then
          ! Integrate part of the element
          sub_coord = ((bounds(i) - node_vals(2)) * (coords(1, 1) - coords(1, 2)) / (node_vals(1) - node_vals(2))) + coords(1, 2)
          integrals(i) = integrals(i) + abs(coords(1, 2) - sub_coord)
        !else
        !  ! Element not integrated
        end if
      end do bounds_loop
      
      if(count(integrate_length) > 0) then
        ! Integrate the whole element
        
        line_length = abs(coords(1, 1) - coords(1, 2))
        
        do i = 1, size(bounds)
          if(integrate_length(i)) then
            integrals(i) = integrals(i) + line_length
          end if
        end do
      end if
    end do ele_loop
    
    call allsum(integrals, communicator = halo_communicator(sfield))
    
    do i = 1, size(bounds)
      ewrite(2, *) "For field " // trim(sfield%name) // " with bound: ", bounds(i)
      ewrite(2, *) "Heaviside integral > bound = ", integrals(i)
    end do
  
  end function heaviside_integral_line

  function heaviside_integral_tri(sfield, bounds, positions) result(integrals)
    type(scalar_field), intent(in) :: sfield
    real, dimension(:), intent(in) :: bounds
    type(vector_field), intent(in) :: positions
    
    real, dimension(size(bounds)) :: integrals
    
    integer :: ele, i, integrated_node
    integer, dimension(2) :: non_integrated_nodes
    integer :: node_integrate_count
    logical, dimension(3) :: node_integrate
    logical, dimension(size(bounds)) :: integrate_area
    real :: tet_area
    real, dimension(3) :: node_vals
    real, dimension(2, 3) :: tri_coords, sub_tri_coords
    
    assert(sfield%mesh == positions%mesh)
    
    integrals = 0.0   
    ele_loop: do ele = 1, ele_count(sfield)
      if(.not. element_owned(sfield, ele)) cycle ele_loop
    
      assert(ele_loc(sfield, ele) == 3)
      tri_coords = ele_val(positions, ele)
      tri_coords(:, 2:) = tri_coords(:, 2:) - spread(tri_coords(:, 1), 2, 2);  tri_coords(:, 1) = 0.0
      node_vals = ele_val(sfield, ele)
      
      integrate_area = .false. 
      bounds_loop: do i = 1, size(bounds)
        node_integrate = node_vals >= bounds(i)
        node_integrate_count = count(node_integrate)
      
        select case(node_integrate_count)
        
          case(0)
            ! Triangle not integrated     
               
          case(1)
            ! Integrate one corner triangle
            
            if(node_integrate(1)) then
              integrated_node = 1
              non_integrated_nodes(1) = 2
              non_integrated_nodes(2) = 3
            else if(node_integrate(2)) then
              integrated_node = 2
              non_integrated_nodes(1) = 1
              non_integrated_nodes(2) = 3
            else
              integrated_node = 3
              non_integrated_nodes(1) = 1
              non_integrated_nodes(2) = 2
            end if
            
            sub_tri_coords(:, 1) = tri_coords(:, integrated_node)
            sub_tri_coords(:, 2) = tri_coords(:, non_integrated_nodes(1)) + &
              & ((bounds(i) - node_vals(non_integrated_nodes(1))) / (node_vals(integrated_node) - node_vals(non_integrated_nodes(1)))) * &
              & (tri_coords(:, integrated_node) - tri_coords(:, non_integrated_nodes(1)))
            sub_tri_coords(:, 3) = tri_coords(:, non_integrated_nodes(2)) + &
              & ((bounds(i) - node_vals(non_integrated_nodes(2))) / (node_vals(integrated_node) - node_vals(non_integrated_nodes(2)))) * &
              & (tri_coords(:, integrated_node) - tri_coords(:, non_integrated_nodes(2)))
              
            sub_tri_coords(:, 2:) = sub_tri_coords(:, 2:) - spread(sub_tri_coords(:, 1), 2, 2)
                      
            integrals(i) = integrals(i) + 0.5 * abs(det(sub_tri_coords(:, 2:)))
          case(2)
            ! Integrate one edge trapezium
            
            integrate_area(i) = .true.  ! Handled outside the bounds loop (below)
            
            if(.not. node_integrate(1)) then
              integrated_node = 1
              non_integrated_nodes(1) = 2
              non_integrated_nodes(2) = 3
            else if(.not. node_integrate(2)) then
              integrated_node = 2
              non_integrated_nodes(1) = 1
              non_integrated_nodes(2) = 3
            else
              integrated_node = 3
              non_integrated_nodes(1) = 1
              non_integrated_nodes(2) = 2
            end if
            
            sub_tri_coords(:, 1) = tri_coords(:, integrated_node)
            sub_tri_coords(:, 2) = tri_coords(:, non_integrated_nodes(1)) + &
              & ((bounds(i) - node_vals(non_integrated_nodes(1))) / (node_vals(integrated_node) - node_vals(non_integrated_nodes(1)))) * &
              & (tri_coords(:, integrated_node) - tri_coords(:, non_integrated_nodes(1)))
            sub_tri_coords(:, 3) = tri_coords(:, non_integrated_nodes(2)) + &
              & ((bounds(i) - node_vals(non_integrated_nodes(2))) / (node_vals(integrated_node) - node_vals(non_integrated_nodes(2)))) * &
              & (tri_coords(:, integrated_node) - tri_coords(:, non_integrated_nodes(2)))
              
            sub_tri_coords(:, 2:) = sub_tri_coords(:, 2:) - spread(sub_tri_coords(:, 1), 2, 2)
            integrals(i) = integrals(i) - 0.5 * abs(det(sub_tri_coords(:, 2:)))
          case(3)
            ! Integrate the whole triangle
            
            integrate_area(i) = .true.  ! Handled outside the bounds loop (below)
          case default
            FLAbort("Unexpected number of integrated nodes in triangle")
        end select
      end do bounds_loop
      
      if(count(integrate_area) > 0) then
        ! Integrate the whole triangle
        
        tet_area = 0.5 * abs(det(tri_coords(:, 2:)))
        
        do i = 1, size(bounds)
          if(integrate_area(i)) then
            integrals(i) = integrals(i) + tet_area
          end if
        end do
      end if
    end do ele_loop
    
    call allsum(integrals, communicator = halo_communicator(sfield))
    
    do i = 1, size(bounds)
      ewrite(2, *) "For field " // trim(sfield%name) // " with bound: ", bounds(i)
      ewrite(2, *) "Heaviside integral > bound = ", integrals(i)
    end do
    
  end function heaviside_integral_tri

  function heaviside_integral_tet(sfield, bounds, positions, tolerance) result(integrals)
    type(scalar_field), intent(in) :: sfield
    real, dimension(:), intent(in) :: bounds
    type(vector_field), intent(in) :: positions
    real, optional, intent(in) :: tolerance
    
    real, dimension(size(bounds)) :: integrals
    
    integer :: ii, iloc, heavi_zero
    integer, dimension(4) :: heavi_non_zero
    real :: dd
    real, dimension(4) :: xele, yele, zele, xelesub, yelesub, zelesub
    real, dimension(6) :: xpent, ypent, zpent

    integer :: ele, i, node, node_integrate_count
    integer, dimension(:), pointer :: nodes
    logical, dimension(size(bounds)) :: integrate_vol
    real :: tet_vol

    real :: l_tolerance
    
    assert(sfield%mesh == positions%mesh)

    if(present(tolerance)) then
       l_tolerance = tolerance
    else
       l_tolerance = epsilon(0.0)
    end if

    integrals = 0.0
    ele_loop: do ele = 1, ele_count(sfield)
       ! Check if this processor is integrating this element
       if(.not. element_owned(sfield, ele)) cycle ele_loop
       
       nodes => ele_nodes(sfield, ele)

       integrate_vol = .false.
       bounds_loop: do i = 1, size(bounds)
       
         node_integrate_count = 0
         heavi_non_zero = 0
         ii = 0

         do iloc = 1,4
            node = nodes(iloc)
            xele(iloc) = node_val(positions, X_, node)
            yele(iloc) = node_val(positions, Y_, node)
            zele(iloc) = node_val(positions, Z_, node)
            if(node_val(sfield, node) >=(bounds(i)-l_tolerance)) then
               node_integrate_count = node_integrate_count + 1
               ii = ii + 1
               heavi_non_zero(ii) = iloc
            endif
         end do

         if (node_integrate_count==4) then
            ! heaviside funtion 1 over whole element so add entire volume
         
            integrate_vol(i) = .true.  ! Handled outside the bounds loop (below)
         elseif(node_integrate_count==0) then
            ! heaviside function zero over whole element so do nothing

         elseif(node_integrate_count==1) then
            ! heaviside function non-zero only over a sub-tet of this element
            
            xelesub(1) = node_val(positions, X_, nodes(heavi_non_zero(1)))
            yelesub(1) = node_val(positions, Y_, nodes(heavi_non_zero(1)))
            zelesub(1) = node_val(positions, Z_, nodes(heavi_non_zero(1)))
            ! find locations of the other sub-tet vertices
            ii = 1 ! already have the first
            do iloc = 1,4
               if(iloc/=heavi_non_zero(1)) then
                  ii = ii + 1
                  dd = (bounds(i) - node_val(sfield, nodes(heavi_non_zero(1))))/ &
                       (node_val(sfield, nodes(iloc))-node_val(sfield, nodes(heavi_non_zero(1))))
                  !dd = max(min(dd,1.0),0.0)
                  xelesub(ii) = xelesub(1) + (node_val(positions, X_, nodes(iloc))-xelesub(1))*dd
                  yelesub(ii) = yelesub(1) + (node_val(positions, Y_, nodes(iloc))-yelesub(1))*dd
                  zelesub(ii) = zelesub(1) + (node_val(positions, Z_, nodes(iloc))-zelesub(1))*dd
               endif
            end do
            if(ii/=4) then
              FLAbort("have not found all three other vertices of the subtet")
            end if
            integrals(i) = integrals(i) + abs(tetvol(xelesub,yelesub,zelesub))

         elseif(node_integrate_count==2) then
            ! this is the complicated case where the ele has been split into two pentahedra
            ! split the pentahedra into 3 tets to calculate the volume
            xpent(1) = node_val(positions, X_, nodes(heavi_non_zero(1)))
            ypent(1) = node_val(positions, Y_, nodes(heavi_non_zero(1)))
            zpent(1) = node_val(positions, Z_, nodes(heavi_non_zero(1)))

            xpent(2) = node_val(positions, X_, nodes(heavi_non_zero(2)))
            ypent(2) = node_val(positions, Y_, nodes(heavi_non_zero(2)))
            zpent(2) = node_val(positions, Z_, nodes(heavi_non_zero(2)))
            ! find locations of the other pentahedra vertices
            ii = 2 ! already have the first two
            do iloc = 1,4
               if((iloc/=heavi_non_zero(1)).and.(iloc/=heavi_non_zero(2))) then
                  ii = ii + 1
                  dd = (bounds(i) - node_val(sfield, nodes(heavi_non_zero(1))))/ &
                       (node_val(sfield, nodes(iloc))-node_val(sfield, nodes(heavi_non_zero(1))))
                  dd = max(min(dd,1.0),0.0)
                  if(ii==3) then
                     xpent(3) = xpent(1) + (node_val(positions, X_, nodes(iloc))-xpent(1))*dd
                     ypent(3) = ypent(1) + (node_val(positions, Y_, nodes(iloc))-ypent(1))*dd
                     zpent(3) = zpent(1) + (node_val(positions, Z_, nodes(iloc))-zpent(1))*dd
                  else
                     xpent(4) = xpent(1) + (node_val(positions, X_, nodes(iloc))-xpent(1))*dd
                     ypent(4) = ypent(1) + (node_val(positions, Y_, nodes(iloc))-ypent(1))*dd
                     zpent(4) = zpent(1) + (node_val(positions, Z_, nodes(iloc))-zpent(1))*dd
                  endif

                  dd = (bounds(i) - node_val(sfield, nodes(heavi_non_zero(2))))/ &
                       (node_val(sfield, nodes(iloc))-node_val(sfield, nodes(heavi_non_zero(2))))

                  dd = max(min(dd,1.0),0.0)
                  if(ii==3) then
                     xpent(5) = xpent(2) + (node_val(positions, X_, nodes(iloc))-xpent(2))*dd
                     ypent(5) = ypent(2) + (node_val(positions, Y_, nodes(iloc))-ypent(2))*dd
                     zpent(5) = zpent(2) + (node_val(positions, Z_, nodes(iloc))-zpent(2))*dd
                  else
                     xpent(6) = xpent(2) + (node_val(positions, X_, nodes(iloc))-xpent(2))*dd
                     ypent(6) = ypent(2) + (node_val(positions, Y_, nodes(iloc))-ypent(2))*dd
                     zpent(6) = zpent(2) + (node_val(positions, Z_, nodes(iloc))-zpent(2))*dd
                  endif
               endif
            end do
            if(ii/=4) then
              FLAbort("have not found all two other vertices of the pent")
            end if
            integrals(i) = integrals(i) + PENTAHEDRON_VOL(Xpent,Ypent,Zpent)

         elseif(node_integrate_count==3) then
            ! in this case there is a subtet over which the heaviside function IS ZERO so just
            ! use (vol(big_tet) - vol(subtet)) as volume of polyhedra where heaviside function non-zero
            heavi_zero = 10 - (heavi_non_zero(1) + heavi_non_zero(2) + heavi_non_zero(3))
            ! heavi_zero is the local node number of the single node where the heaviside function is zero
            xelesub(1) = node_val(positions, X_, nodes(heavi_zero))
            yelesub(1) = node_val(positions, Y_, nodes(heavi_zero))
            zelesub(1) = node_val(positions, Z_, nodes(heavi_zero))
            ! find locations of the other sub-tet vertices
            ii = 1 ! already have the first
            do iloc = 1,4
               if(iloc/=heavi_zero) then
                  ii = ii + 1
                  dd = (bounds(i) - node_val(sfield, nodes(heavi_zero)))/ &
                       (node_val(sfield, nodes(iloc))-node_val(sfield, nodes(heavi_zero)))

                  dd = max(min(dd,1.0),0.0)

                  xelesub(ii) = xelesub(1) + (node_val(positions, X_, nodes(iloc))-xelesub(1))*dd
                  yelesub(ii) = yelesub(1) + (node_val(positions, Y_, nodes(iloc))-yelesub(1))*dd
                  zelesub(ii) = zelesub(1) + (node_val(positions, Z_, nodes(iloc))-zelesub(1))*dd
               endif
            end do
            if(ii/=4) then
              FLAbort("have not found all three other vertices of the subtet")
            end if
            integrate_vol(i) = .true.  ! Handled outside the bounds loop (below)
            integrals(i) = integrals(i) - abs(tetvol(xelesub,yelesub,zelesub))
         endif
       end do bounds_loop
       
       if(count(integrate_vol) > 0) then
         tet_vol = abs(tetvol(xele,yele,zele))
         do i = 1, size(bounds)
           if(integrate_vol(i)) then
             ! heaviside funtion 1 over whole element so add entire volume
             
             integrals(i) = integrals(i) + tet_vol
           end if
         end do
       end if
    end do ele_loop
    
    call allsum(integrals, communicator = halo_communicator(sfield))
    
    do i = 1, size(bounds)
      ewrite(2, *) "For field " // trim(sfield%name) // " with bound: ", bounds(i)
      ewrite(2, *) "Heaviside integral > bound = ", integrals(i)
    end do

  end function heaviside_integral_tet
  
end module mixing_statistics
