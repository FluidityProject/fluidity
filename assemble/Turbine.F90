!    Copyright (C) 2008 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
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

module turbine

  use spud
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use futils
  use fields
  implicit none

 private
  public :: turbine_check_options, construct_turbine_interface

contains


  subroutine construct_turbine_interface(turbine_fluxfac, theta, dt, ele, face, face_2, ni, &
       & big_m_tensor_addto, rhs_addto, X, U, velocity_bc, velocity_bc_type)
    !!< Construct the DG element boundary integrals on the ni-th face of
    !!< element ele.
    implicit none

    ! The turbine model
    real, intent(out) :: turbine_fluxfac
    real, intent(in) :: theta, dt
    integer, intent(in) :: ele, face, face_2, ni
    real, dimension(:,:,:,:), intent(inout) :: big_m_tensor_addto
    real, dimension(:,:) :: rhs_addto
    ! We pass these additional fields to save on state lookups.
    type(vector_field), intent(in) :: X, U

    !! Boundary conditions associated with this interface (if any).
    type(vector_field), intent(in) :: velocity_bc
    integer, dimension(:,:), intent(in) :: velocity_bc_type
    integer :: turbine_model

    turbine_model = velocity_bc_type(1,face)

    turbine_fluxfac=-1.0
    select case (turbine_model) 
    case(5)
         turbine_fluxfac=1.0
         return
    case(4)
         call construct_turbine_interface_penalty(theta, dt, ele, face, face_2, ni, &
            & big_m_tensor_addto, rhs_addto, X, U,&
            & velocity_bc, velocity_bc_type)
    case default
         FLAbort("Unknown turbine model found.")
    end select
  end subroutine construct_turbine_interface

  subroutine construct_turbine_interface_penalty(theta, dt, ele, face, face_2, ni, &
       & big_m_tensor_addto, rhs_addto, X, U,  &
       velocity_bc, velocity_bc_type)

    real, intent(in) :: theta, dt
    integer, intent(in) :: ele, face, face_2, ni
    real, dimension(:,:,:,:), intent(inout) :: big_m_tensor_addto
    real, dimension(:,:) :: rhs_addto
    ! We pass these additional fields to save on state lookups.
    type(vector_field), intent(in) :: X, U
    !! Boundary conditions associated with this interface (if any).
    type(vector_field), intent(in) :: velocity_bc
    integer, dimension(:,:), intent(in) :: velocity_bc_type

    real, dimension(face_loc(U,face),face_loc(U,face_2)) :: penalty_domain_connecting_in, penalty_domain_connecting_out
    integer :: dim, i, start, finish
    real, dimension(face_ngi(U,face)) :: detwei
    ! Face objects and numberings.
    type(element_type), pointer :: u_shape, u_shape_2
    integer, dimension(face_loc(U,face)) :: u_face_l
    real, dimension(velocity_bc%dim, velocity_bc%mesh%shape%loc) :: penalty_val, penalty_val_2

        ! Connected two domains with a penalty term
        !
        ! The penalty term has the form Int_E p*phi*(\delta u_u- + u- - (delta_u+ + u+)) dE where p is the penalty parameter
        !
        penalty_val=ele_val(velocity_bc, face)
#ifdef DDEBUG
        ! The penalty values should be the same on both turbine sides! 
        penalty_val_2=ele_val(velocity_bc, face_2)
        assert(all(penalty_val(1,:)==penalty_val_2(1,:)))
#endif

        start=ele_loc(u,ele)+(ni-1)*face_loc(U, face_2)+1
        finish=start+face_loc(U, face_2)-1
    
        u_shape=>face_shape(U, face)
        u_shape_2=>face_shape(U, face_2)

        !----------------------------------------------------------------------
        ! Change of coordinates on face.
        !----------------------------------------------------------------------
        call transform_facet_to_physical(X, face,&
           &                          detwei_f=detwei) 

        penalty_domain_connecting_in=shape_shape(U_shape, U_shape, detwei)
        penalty_domain_connecting_out=shape_shape(U_shape, U_shape_2, detwei)
        ! Multiply Matrix with the penalty function -- !!!to be reviewed!!! --
        do i=1, face_loc(U,face_2)
            penalty_domain_connecting_in(i,:)=penalty_domain_connecting_in(i,:)*penalty_val(1,:) 
            penalty_domain_connecting_out(i,:)=penalty_domain_connecting_out(i,:)*penalty_val(1,:)
        end do
        u_face_l=face_local_nodes(U, face)
        do dim = 1, u%dim
            ! Insert penalty terms in matrix.
            big_m_tensor_addto(dim, dim, u_face_l, u_face_l) = &
               big_m_tensor_addto(dim, dim, u_face_l, u_face_l) + &
               penalty_domain_connecting_in*dt*theta

            rhs_addto(dim,u_face_l) = rhs_addto(dim,u_face_l) &
                -matmul(penalty_domain_connecting_in,face_val(U,dim,face))

            big_m_tensor_addto(dim, dim, u_face_l, start:finish) = &
              big_m_tensor_addto(dim, dim, u_face_l, start:finish) - &
              penalty_domain_connecting_out*dt*theta

            rhs_addto(dim,u_face_l) = rhs_addto(dim,u_face_l) &
                +matmul(penalty_domain_connecting_out,face_val(U,dim,face_2))
        end do
  end subroutine construct_turbine_interface_penalty

  subroutine turbine_check_options
    character(len=OPTION_PATH_LEN):: turbine_path, turbine_name, bc_name
    integer :: notur, i, j
    logical :: have_dirichlet_model, have_flux_dg_model, have_flux_penalty_model

    ! Don't check turbine configuration if it's not included in the model!
    if (.not.have_option("/turbine_model")) return

    have_dirichlet_model=.false.
    have_flux_penalty_model=.false.
    have_flux_dg_model=.false.

    ! loop through turbines
    notur = option_count("/turbine_model/turbine")
    do i=0, notur-1
       turbine_path="/turbine_model/turbine["//int2str(i)//"]"
       if (have_option(trim(turbine_path)//"/dirichlet")) then
           have_dirichlet_model=.true.
           ! The specified b.c.'s  in the turbine model must be dirichlet boundary conditions with normal_component.
           do j=1,2
              call get_option("/turbine_model/turbine["//int2str(i)//"]/dirichlet/boundary_condition_name_"//int2str(j)//"/name", bc_name)
              if (.not. have_option("/material_phase[0]/vector_field::Velocity/prognostic/boundary_conditions::"//trim(bc_name)//"/type::dirichlet/align_bc_with_surface/normal_component")) then
                 call get_option("/turbine_model/turbine["//int2str(i)//"]/name", turbine_name)
                 FLExit("Error while checking the options for turbine '"//trim(turbine_name)//"': Turbine model boundary has to be dirichlet boundary conditions with a normal_component.")
              end if
           end do
       elseif (have_option(trim(turbine_path)//"/flux/dg")) then
            have_flux_dg_model=.true.
       elseif (have_option(trim(turbine_path)//"/flux/penalty")) then
            have_flux_penalty_model=.true.
       else
          FLAbort("Unknown turbine model specified!")
       end if
    end do

    if (have_dirichlet_model) then
       ! We need the FreeSurface field.
       if (.not.have_option("/material_phase[0]/scalar_field::FreeSurface")) then
          FLExit("Turbine modelling requires FreeSurface to be activated.")
       end if
    end if

  end subroutine turbine_check_options

end module turbine
