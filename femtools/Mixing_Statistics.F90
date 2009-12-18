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
  use global_parameters, only:FIELD_NAME_LEN,ACCTIM, OPTION_PATH_LEN, halo_tag, halo_tag_p
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
     module procedure heaviside_integral_old, heaviside_integral_new   
  end interface

  interface mixing_stats
     module procedure mixing_stats_scalar
  end interface

  private
  public mixing_stats

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
    real :: tolerance, total_volume
    real, dimension(:), allocatable :: mixing_bin_bounds
    type(scalar_field) :: lumped_mass


    call get_option(trim(complete_field_path(sfield%option_path)) &
         &// "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/tolerance",&
         & tolerance, default = epsilon(0.0))

    allocate(mixing_bin_bounds(1:size(f_mix_fraction)))
    call get_option(trim(complete_field_path(sfield%option_path)) &
         &// "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/mixing_bin_bounds",&
         & mixing_bin_bounds)

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

  subroutine heaviside_integral_new(f_mix_fraction, sfield, Xfield, mixing_stats_count)

    real, dimension(:), intent(out) :: f_mix_fraction
    type(scalar_field), target, intent(inout) :: sfield
    type(vector_field), target, intent(in) :: Xfield
    integer, intent(in) :: mixing_stats_count


    integer, dimension(:), allocatable :: ndglno
    integer :: j, ele
    integer :: totele, nonods, nloc
    real :: total_volume
    real, dimension(:), allocatable :: x, y, z, mixing_bin_bounds, detwei
    real :: tolerance
    !type(vector_field), pointer:: positions
    type(element_type), pointer :: shape

    shape => null()      

    nloc=ele_loc(sfield,1)

    if((nloc.eq.4).and.(Xfield%dim==3) .and. continuity(sfield)>=0) then ! this only works for linear tets

       totele=element_count(sfield)
       nonods = node_count(sfield)

       allocate(ndglno(1:(totele*nloc)))
       ndglno=sfield%mesh%ndglno

       allocate(x(1:nonods))
       allocate(y(1:nonods))
       allocate(z(1:nonods))
       x = Xfield%val(1)%ptr
       y = Xfield%val(2)%ptr
       z = Xfield%val(3)%ptr

       allocate(mixing_bin_bounds(1:size(f_mix_fraction)))            

       call get_option(trim(complete_field_path(sfield%option_path)) &
            &// "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/mixing_bin_bounds",&
            & mixing_bin_bounds)

       call get_option(trim(complete_field_path(sfield%option_path)) &
            &// "/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/tolerance",&
            & tolerance, default = epsilon(0.0))

       do j=1,size(f_mix_fraction)
          call heaviside_integral(integral=f_mix_fraction(j), scalar_field=sfield%val, &
               &heavi_value = mixing_bin_bounds(j), x=x, y=y, z=z, nonods = nonods, &
               &totele=totele, nloc=nloc, ndglno=ndglno, tolerance&
               &=tolerance, mesh=sfield%mesh)
       end do

       deallocate(mixing_bin_bounds)

       ! In f_mix_fraction(j) have volume of sfield >mixing_bin_bounds(j)
       ! Subtract so that f_mix_fraction(j) has 
       ! mixing_bin_bounds(j)<volume fraction of sfield<mixing_bins_bounds(j+1) 
       ! and if normalise is selected divide by total volume
       if(have_option(trim(complete_field_path(sfield%option_path))//&
            &"/stat/include_mixing_stats["// int2str(mixing_stats_count) // "]/continuous_galerkin/normalise")) then

          assert(ele_count(sfield) > 0)
          shape => ele_shape(sfield, 1)
          allocate(detwei(shape%ngi))
          !positions => extract_vector_field(state, "Coordinate")

          total_volume = 0.0

          do ele = 1, ele_count(sfield)
             if(.not. element_owned(sfield,ele)) cycle

             call transform_to_physical(Xfield, ele, detwei = detwei)
             total_volume = total_volume + sum(detwei)
          end do
          call allsum(total_volume)

          deallocate(detwei)

          do j = 1,(size(f_mix_fraction)-1)
             f_mix_fraction(j) = (f_mix_fraction(j)-f_mix_fraction(j+1))&
                  &/total_volume
          end do
          f_mix_fraction(size(f_mix_fraction)) = &
               & f_mix_fraction(size(f_mix_fraction))/total_volume


       else

          do j = 1,(size(f_mix_fraction)-1)
             f_mix_fraction(j) = (f_mix_fraction(j)-f_mix_fraction(j+1))
          end do
          f_mix_fraction(size(f_mix_fraction)) = &
               & f_mix_fraction(size(f_mix_fraction))

       endif

    else
       FLAbort("continuous_galerkin mixing_stats only work for linear tets")
    endif

    return
  end subroutine heaviside_integral_new

  subroutine heaviside_integral_old(integral,scalar_field,heavi_value, &
       x,y,z,nonods,totele,nloc,ndglno,tolerance, mesh)
    real, intent(out)::integral
    real, intent(inout) :: scalar_field(nonods)
    real, intent(in) :: heavi_value
    integer, intent(in)::nonods,totele,nloc,ndglno(totele*nloc)
    real, intent(in)::x(nonods),y(nonods),z(nonods)
    !c
    integer number_nodes,heavi_non_zero(4),heavi_zero
    real xele(4),yele(4),zele(4),xelesub(4),yelesub(4),zelesub(4)
    real xpent(6),ypent(6),zpent(6)
    integer ii,ele,iloc,globi
    real dd,a1,a2
    integer counter1,counter2,counter3,counter4,counter5
    type(mesh_type), intent(in) :: mesh

    real, intent(in), optional :: tolerance
    real :: l_tolerance

    if(present(tolerance)) then
       l_tolerance = tolerance
    else
       l_tolerance = epsilon(0.0)
    end if

    if(halo_count(mesh) > 0) then
      call halo_update(mesh%halos(1), scalar_field)
    end if

    if(nloc.eq.4) then ! this only works for linear tets
       counter1=0
       counter2=0
       counter3=0
       counter4=0
       counter5=0
       integral = 0.0

       do ele = 1,totele
          ! Check if this processor is integrating this element
          if(.not. element_owned(mesh, ele)) cycle

          !c          print *,'#####################'
          !c          print *,'ele',ele
          !c          print *,'scalar_field',(scalar_field(ndglno((ele-1)*nloc+iloc)),iloc=1,4)
          !c          print *,'x',(x(ndglno((ele-1)*nloc+iloc)),iloc=1,4)
          number_nodes = 0
          heavi_non_zero = 0
          ii = 0

          do iloc = 1,nloc
             globi = ndglno((ele-1)*nloc+iloc)
             xele(iloc) = x(globi)
             yele(iloc) = y(globi)
             zele(iloc) = z(globi)
             if(scalar_field(globi).ge.(heavi_value-l_tolerance)) then
                number_nodes = number_nodes+1
                ii = ii + 1
                heavi_non_zero(ii) = iloc
             endif
          enddo
          !c         print *,'number_nodes,heavi_value',number_nodes,heavi_value
          !c         print *,'heavi_non_zero',heavi_non_zero

          if (number_nodes.eq.nloc) then
             counter1=counter1+1     ! heaviside funtion 1 over whole element so add entire volume
             integral = integral + abs(tetvol(xele,yele,zele))
             !c           print *,'integral',integral

          elseif(number_nodes.eq.0) then
             counter2=counter2+1  ! heaviside function zero over whole element so do nothing

          elseif(number_nodes.eq.1) then
             counter3=counter3+1   ! heaviside function non-zero only over a sub-tet of this element
             xelesub(1) = x(ndglno((ele-1)*nloc + heavi_non_zero(1)))
             yelesub(1) = y(ndglno((ele-1)*nloc + heavi_non_zero(1)))
             zelesub(1) = z(ndglno((ele-1)*nloc + heavi_non_zero(1)))
             ! find locations of the other sub-tet vertices
             ii = 1 ! already have the first
             do iloc = 1,nloc
                if(iloc.ne.heavi_non_zero(1)) then
                   ii = ii + 1
                   dd = (heavi_value - scalar_field(ndglno((ele-1)*nloc + heavi_non_zero(1))))/ &
                        (scalar_field(ndglno((ele-1)*nloc +iloc))-scalar_field(ndglno((ele-1)*nloc + heavi_non_zero(1))))
                   !dd = max(min(dd,1.0),0.0)
                   xelesub(ii) = xelesub(1) + (x(ndglno((ele-1)*nloc +iloc))-xelesub(1))*dd
                   yelesub(ii) = yelesub(1) + (y(ndglno((ele-1)*nloc +iloc))-yelesub(1))*dd
                   zelesub(ii) = zelesub(1) + (z(ndglno((ele-1)*nloc +iloc))-zelesub(1))*dd
                endif
             enddo
             if(ii.ne.4) stop 3752 ! have not found all three other vertices of the subtet
             integral = integral + abs(tetvol(xelesub,yelesub,zelesub))
             !c           print *,'integral',integral

          elseif(number_nodes.eq.2) then
             counter4=counter4+1
             ! this is the complicated case where the ele has been split into two pentahedra
             ! split the pentahedra into 3 tets to calculate the volume
             xpent(1) = x(ndglno((ele-1)*nloc + heavi_non_zero(1)))
             ypent(1) = y(ndglno((ele-1)*nloc + heavi_non_zero(1)))
             zpent(1) = z(ndglno((ele-1)*nloc + heavi_non_zero(1)))

             xpent(2) = x(ndglno((ele-1)*nloc + heavi_non_zero(2)))
             ypent(2) = y(ndglno((ele-1)*nloc + heavi_non_zero(2)))
             zpent(2) = z(ndglno((ele-1)*nloc + heavi_non_zero(2)))
             ! find locations of the other pentahedra vertices
             ii = 2 ! already have the first two
             do iloc = 1,nloc
                if((iloc.ne.heavi_non_zero(1)).and.(iloc.ne.heavi_non_zero(2))) then
                   ii = ii + 1
                   dd = (heavi_value - scalar_field(ndglno((ele-1)*nloc + heavi_non_zero(1))))/ &
                        (scalar_field(ndglno((ele-1)*nloc +iloc))-scalar_field(ndglno((ele-1)*nloc + heavi_non_zero(1))))
                   dd = max(min(dd,1.0),0.0)
                   if(ii.eq.3) then
                      xpent(3) = xpent(1) + (x(ndglno((ele-1)*nloc +iloc))-xpent(1))*dd
                      ypent(3) = ypent(1) + (y(ndglno((ele-1)*nloc +iloc))-ypent(1))*dd
                      zpent(3) = zpent(1) + (z(ndglno((ele-1)*nloc +iloc))-zpent(1))*dd
                      !print *,'zpent(3)',zpent(3)
                   else
                      xpent(4) = xpent(1) + (x(ndglno((ele-1)*nloc +iloc))-xpent(1))*dd
                      ypent(4) = ypent(1) + (y(ndglno((ele-1)*nloc +iloc))-ypent(1))*dd
                      zpent(4) = zpent(1) + (z(ndglno((ele-1)*nloc +iloc))-zpent(1))*dd
                      !print *,'zpent(4)',zpent(4)
                   endif

                   dd = (heavi_value - scalar_field(ndglno((ele-1)*nloc + heavi_non_zero(2))))/ &
                        (scalar_field(ndglno((ele-1)*nloc +iloc))-scalar_field(ndglno((ele-1)*nloc + heavi_non_zero(2))))

                   dd = max(min(dd,1.0),0.0)
                   if(ii.eq.3) then
                      xpent(5) = xpent(2) + (x(ndglno((ele-1)*nloc +iloc))-xpent(2))*dd
                      ypent(5) = ypent(2) + (y(ndglno((ele-1)*nloc +iloc))-ypent(2))*dd
                      zpent(5) = zpent(2) + (z(ndglno((ele-1)*nloc +iloc))-zpent(2))*dd
                   else
                      xpent(6) = xpent(2) + (x(ndglno((ele-1)*nloc +iloc))-xpent(2))*dd
                      ypent(6) = ypent(2) + (y(ndglno((ele-1)*nloc +iloc))-ypent(2))*dd
                      zpent(6) = zpent(2) + (z(ndglno((ele-1)*nloc +iloc))-zpent(2))*dd
                   endif
                endif
             enddo
             if(ii.ne.4) stop 3754 ! have not found all two other vertices of the pent
             integral = integral + PENTAHEDRON_VOL(Xpent,Ypent,Zpent)
             !c           print *,'integral',integral

          elseif(number_nodes.eq.3) then
             counter5=counter5+1
             ! in this case there is a subtet over which the heaviside function IS ZERO so just
             ! use (vol(big_tet) - vol(subtet)) as volume of polyhedra where heaviside function non-zero
             heavi_zero = 10 - (heavi_non_zero(1) + heavi_non_zero(2) + heavi_non_zero(3))
             ! heavi_zero is the local node number of the single node where the heaviside function is zero
             xelesub(1) = x(ndglno((ele-1)*nloc + heavi_zero))
             yelesub(1) = y(ndglno((ele-1)*nloc + heavi_zero))
             zelesub(1) = z(ndglno((ele-1)*nloc + heavi_zero))
             ! find locations of the other sub-tet vertices
             ii = 1 ! already have the first
             do iloc = 1,nloc
                if(iloc.ne.heavi_zero) then
                   ii = ii + 1
                   dd = (heavi_value - scalar_field(ndglno((ele-1)*nloc + heavi_zero)))/ &
                        (scalar_field(ndglno((ele-1)*nloc +iloc))-scalar_field(ndglno((ele-1)*nloc + heavi_zero)))

                   dd = max(min(dd,1.0),0.0)

                   xelesub(ii) = xelesub(1) + (x(ndglno((ele-1)*nloc +iloc))-xelesub(1))*dd
                   yelesub(ii) = yelesub(1) + (y(ndglno((ele-1)*nloc +iloc))-yelesub(1))*dd
                   zelesub(ii) = zelesub(1) + (z(ndglno((ele-1)*nloc +iloc))-zelesub(1))*dd
                endif
             enddo
             if(ii.ne.4) stop 3753 ! have not found all three other vertices of the subtet
             integral = integral + (abs(tetvol(xele,yele,zele)) - abs(tetvol(xelesub,yelesub,zelesub)))      
             a1 = tetvol(xele,yele,zele)
             a2 = tetvol(xelesub,yelesub,zelesub)
          endif
       enddo

       call allsum(integral)
       !ewrite(3,*)'counters 1-5:',counter1,counter2,counter3,counter4,counter5
       !ewrite(3,*)'integral:',integral

    else
       integral=0.0
       ewrite(2,*) 'Gone into heaviside_integral but it ', &
            & 'only works for tets so have returned zero.'
    endif

    return
  end subroutine heaviside_integral_old

end module mixing_statistics
