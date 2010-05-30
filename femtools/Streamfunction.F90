!    Copyright (C) 2010 Imperial College London and others.
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

module streamfunction
  use state_module
  use fields
  use sparse_tools
  use spud
  use global_parameters, only: OPTION_PATH_LEN
  use sparsity_patterns
  use solvers
  use boundary_conditions
  use vector_tools
  use transform_elements
  use eventcounter
  implicit none

  private
  public :: calculate_stream_function_multipath_2d

  type(integer_vector), dimension(:), allocatable, save :: flux_face_list
  real, dimension(:,:), allocatable, save :: flux_normal

  integer, save :: last_adapt=-1

contains

  subroutine find_stream_paths(X, streamfunc)
    !!< Find the paths through the mesh over which the velocity will be
    !!< integrated to calculate the flux.
    type(vector_field), intent(in) :: X
    type(scalar_field), intent(inout) :: streamfunc

    type(ilist) :: tmp_face_list
    
    integer :: bc_count
    character(len=OPTION_PATH_LEN) :: option_path
    real, dimension(2) :: start, end, dx,  face_c

    integer, dimension(:), pointer :: neigh
    integer :: face, p, ele1, ele2, e2, i
    real :: dx2, c1, c2

    bc_count=get_boundary_condition_count(streamfunc)
       
    if(.not.(allocated(flux_face_list))) then
       allocate(flux_face_list(bc_count))
       allocate(flux_normal(2, bc_count))
    end if


    bc_loop: do i=1, bc_count
       call get_boundary_condition(streamfunc, i, option_path=option_path)
       
       if (have_option(trim(option_path)//"/primary_boundary")) then
          ! No path on the primary boundary.

          if (associated(flux_face_list(i)%ptr)) then
             deallocate(flux_face_list(i)%ptr)
          end if
          allocate(flux_face_list(i)%ptr(0))

          cycle bc_loop
       end if
       
       ! We must be on a secondary boundary.
       call get_option(trim(option_path)//"/secondary_boundary/primary_point"&
            &, start) 
       call get_option(trim(option_path)//"/secondary_boundary/secondary_point"&
            &, end) 
       
       dx=start-end
       dx2=dot_product(dx,dx)

       do ele1=1, element_count(streamfunc)
          neigh=>ele_neigh(streamfunc, ele1)
          
          do e2=1,size(neigh)
             ele2=neigh(e2)
             ! Don't do boundaries
             if (ele2<=0) cycle
             ! Do each edge only once
             if (ele1>ele2) cycle

             face=ele_face(streamfunc, ele1, ele2)

             face_c=sum(face_val(X,face),2)/face_loc(X,face)

             p=dot_product(face_c,dx)/dx2

             ! If the face is not within the limits of the line, don't do it.
             if (p<0.or.p>1) cycle

             c1=cross_product2(sum(ele_val(X,ele1),2)/ele_loc(X,ele1)-start, dx)
             c2=cross_product2(sum(ele_val(X,ele2),2)/ele_loc(X,ele2)-start, dx)
             
             if(C1<0 .and. C2>=0) then
                continue

             else if (C1>=0 .and. C2<0) then
                continue
             else
                cycle
             end if

             call insert(tmp_face_list, face)

          end do

       end do

       if (associated(flux_face_list(i)%ptr)) then
          deallocate(flux_face_list(i)%ptr)
       end if
       allocate(flux_face_list(i)%ptr(tmp_face_list%length))
       
       flux_face_list(i)%ptr=list2vector(tmp_face_list)
       call flush_list(tmp_face_list)

       ! Work out the orthonormal to the line.
       dx=dx/sqrt(dx2)
       flux_normal(:,i)=(/-dx(2), dx(1)/)
       
    end do bc_loop
    
  end subroutine find_stream_paths

  function boundary_value(X, U, bc_num)
    !!< Calculate the value of the streamfunction on the boundary provided
    !!< by integrating the velocity flux across a line between this
    !!< boundary and the primary boundary. 
    real :: boundary_value
    type(vector_field), intent(in) :: X, U
    integer, intent(in) :: bc_num
    
    integer :: face, i

    boundary_value=0.0

    do i=1, size(flux_face_list(bc_num)%ptr)
       face=flux_face_list(bc_num)%ptr(i)
       
       boundary_value=boundary_value + face_flux(face, X, U, bc_num)

    end do
    
  contains

    function face_flux(face, X, U, bc_num)
      real :: face_flux
      integer, intent(in) :: face, bc_num
      type(vector_field), intent(in) :: X, U
      
      real, dimension(face_ngi(U, face)) :: detwei
      real, dimension(U%dim,face_ngi(U, face)) :: normal, U_quad
      integer :: gi

      call transform_facet_to_physical(X, face, detwei_f=detwei, normal=normal)

      U_quad=face_val_at_quad(U,face)
      
      face_flux=0.0

      do gi=1, size(detwei)
         face_flux=face_flux+abs(dot_product(flux_normal(:,bc_num),normal(:,gi)))&
              * dot_product(U_quad(:,gi),flux_normal(:,bc_num))&
              * detwei(gi)
      end do

    end function face_flux

  end function boundary_value

  subroutine calculate_stream_function_multipath_2d(state, streamfunc, stat)
    !!< Calculate the stream function for a 
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: streamfunc
    integer, intent(out), optional :: stat
    
    integer :: i, lstat, ele
    type(vector_field), pointer :: X, U
    type(csr_sparsity) :: psi_sparsity
    type(csr_matrix) :: psi_mat
    type(scalar_field) :: rhs
    real :: flux_val
    type(scalar_field), pointer :: surface_field

    do i = 1, 2
       select case(i)
       case(1)
          X => extract_vector_field(state, "Coordinate", stat)
       case(2)
          U => extract_vector_field(state, "Velocity", stat)
       case default
          FLAbort("Invalid loop index")
       end select
       if(present_and_nonzero(stat)) then
          return
       end if
    end do
    
    if (X%dim/=2) then
       FLExit("Streamfunction is only valid in 2d")
    end if
    ! No discontinuous stream functions.
    if (continuity(streamfunc)<0) then
       FLExit("Streamfunction must be a continuous field")    
    end if

    if (last_adapt<eventcount(EVENT_ADAPTIVITY)) then
       last_adapt=eventcount(EVENT_ADAPTIVITY)
       
       call find_stream_paths(X, streamfunc)
    end if
    
    psi_sparsity = extract_csr_sparsity(state, &
               &                      "StreamFunctionSparsity", lstat)
    if (lstat/=0) then
       psi_sparsity = make_sparsity(streamfunc%mesh, streamfunc%mesh, &
            "StreamFunctionSparsity")
    else
       call incref(psi_sparsity)
    end if
    
    call allocate(psi_mat, psi_sparsity, name="StreamFunctionMatrix")

    call zero(psi_mat)
    call allocate(rhs, streamfunc%mesh, "StreamFunctionRHS")
    call zero(rhs)

    do ele=1, element_count(streamfunc)
       
       call calculate_streamfunc_ele(psi_mat, rhs, ele, X, U)

    end do

    do i = 1, get_boundary_condition_count(streamfunc)

       surface_field=>extract_surface_field(streamfunc, i, "value")

       flux_val=boundary_value(X,U,i)

       call set(surface_field, flux_val)

    end do

    call zero(streamfunc)

    call apply_dirichlet_conditions(psi_mat, rhs, streamfunc)

    call petsc_solve(streamfunc, psi_mat, rhs)

    call deallocate(rhs)
    call deallocate(psi_mat)
    call deallocate(psi_sparsity)

  contains
    
    subroutine calculate_streamfunc_ele(psi_mat, rhs, ele, X, U)
      type(csr_matrix), intent(inout) :: psi_mat
      type(scalar_field), intent(inout) :: rhs
      type(vector_field), intent(in) :: X,U
      integer, intent(in) :: ele

      ! Transformed gradient function for velocity.
      real, dimension(ele_loc(U, ele), ele_ngi(U, ele), mesh_dim(U)) :: du_t
      ! Ditto for the stream function, psi
      real, dimension(ele_loc(rhs, ele), ele_ngi(rhs, ele), mesh_dim(rhs))&
           & :: dpsi_t 

      ! Local vorticity_matrix
      real, dimension(2, ele_loc(rhs, ele), ele_loc(U, ele)) ::&
           & lvorticity_mat
      ! Local vorticity
      real, dimension(ele_loc(rhs, ele)) :: lvorticity

      ! Variable transform times quadrature weights.
      real, dimension(ele_ngi(U,ele)) :: detwei
      
      type(element_type), pointer :: U_shape, psi_shape
      integer, dimension(:), pointer :: psi_ele
      integer :: i

      U_shape=> ele_shape(U, ele)
      psi_shape=> ele_shape(rhs, ele)
      psi_ele=>ele_nodes(rhs, ele)
      
      ! Transform U derivatives and weights into physical space.
      call transform_to_physical(X, ele, U_shape, dshape=du_t, detwei=detwei)
      ! Ditto psi.
      call transform_to_physical(X, ele, psi_shape, dshape=dpsi_t)

      call addto(psi_mat, psi_ele, psi_ele, &
           dshape_dot_dshape(dpsi_t, dpsi_t, detwei))

      lvorticity_mat=shape_curl_shape_2d(psi_shape, du_t, detwei)
      
      lvorticity=0.0
      do i=1,2
         lvorticity=lvorticity &
              +matmul(lvorticity_mat(i,:,:), ele_val(U, i, ele))
      end do
      
      call addto(rhs, psi_ele, -lvorticity)
      
    end subroutine calculate_streamfunc_ele

  end subroutine calculate_stream_function_multipath_2d
  
end module streamfunction
