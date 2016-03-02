!    Copyright (C) 2009 Imperial College London and others.
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

module slope_limiters_dg
use fldebug
use ieee_arithmetic
use spud
use vector_tools, only: solve
use elements
use eventcounter
use transform_elements
use sparse_tools
use fields
use field_options, only: find_linear_parent_mesh
use state_module
use vtk_interfaces
use state_fields_module
use bound_field_module
implicit none

private
public limit_slope_dg, limit_fpn, limit_vb

integer, parameter :: LIMITER_MINIMAL=1
integer, parameter :: LIMITER_COCKBURN=2
integer, parameter :: LIMITER_HERMITE_WENO=3
integer, parameter :: LIMITER_FPN=4
integer, parameter :: LIMITER_VB=5

public :: LIMITER_MINIMAL, LIMITER_COCKBURN, LIMITER_HERMITE_WENO,&
     & LIMITER_FPN, LIMITER_VB

!!CockburnShuLimiter stuff
real :: TVB_factor=5.0
real :: Limit_factor=1.1
real, dimension(:,:,:), pointer :: alpha => null()
real, dimension(:,:), pointer :: dx2 => null(), A => null()
integer :: CSL_adapt_counter = -666
logical :: CSL_initialised = .false.
logical :: tolerate_negative_weights

!!Hermite Weno limiter stuff
real :: gam0 !power coefficient in weights
real :: eps_o !relative/absolute tolerance threshold for oscillation indicator
real :: eps_w !relative/absolute tolerance threshold for WENO weights
real :: disc_tol !Value for discontinuity test
real :: limit_tol !Do not limit if infinity norm of tracer is less than
!this value on an element
logical :: debugging !Switch to bung out lots of debugging output
integer, parameter :: IGNORE_MISSING_POLYS=1
integer, parameter :: REPLACE_MISSING_POLYS=2
integer, parameter :: LOWER_ORDER=3
integer :: missing_polys
logical :: leave_out_hermite_polynomials
logical :: has_discontinuity_detector_field
type(scalar_field), pointer :: discontinuity_detector_field
integer :: limit_count

contains

  subroutine limit_slope_dg(T, U, X, state, limiter)
    !! Assume 1D linear elements
    type(scalar_field), intent(inout) :: T
    type(vector_field), intent(in) :: X, U
    type(state_type), intent(inout) :: state
    integer, intent(in) :: limiter

    integer :: ele, stat
    type(scalar_field) :: T_limit

    !assert(mesh_dim(coordinate)==1)
    !assert(field%mesh%continuity<0)
    !assert(field%mesh%shape%degree==1)

    ewrite(2,*) 'subroutiune limit_slope_dg'

    select case (limiter)
    case (LIMITER_MINIMAL)
       T_limit=extract_scalar_field(state, trim(T%name)//"Limiter", stat=stat)
       
       do ele=1,element_count(T)
          
          if (stat==0) then
             call limit_slope_ele_dg(ele, T, X, T_limit)
          else
             call limit_slope_ele_dg(ele, T, X)
          end if

       end do

    case (LIMITER_COCKBURN)

       call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Cockburn_Shu/TVB_factor", &
            &TVB_factor)
       call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Cockburn_Shu/limit_factor", &
            &limit_factor)

       tolerate_negative_weights = &
            &have_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Cockburn_Shu/tolerate_negative_weights")

       call cockburn_shu_setup(T, X)
       
       do ele=1,element_count(T)
          
          call limit_slope_ele_cockburn_shu(ele, T, X)
          
       end do

    case (LIMITER_HERMITE_WENO)

       call get_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/power_coeffi&
            &cient", &
       & gam0)
       call get_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/tolerance_th&
            &reshold_oscillations", &
            &eps_o) 
       call get_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/tolerance_th&
            &reshold_weights", &
            &eps_w) 
       call get_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/discontinuit&
            &y_tolerance",disc_tol) 
       call get_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/limit_tolera&
            &nce",limit_tol) 
       debugging = have_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/debugging")
       missing_polys = IGNORE_MISSING_POLYS
       if(have_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/&
            &boundary_treatment::ignore_missing_polys"))&
            & missing_polys = IGNORE_MISSING_POLYS
       if(have_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/&
            &boundary_treatment::replace_missing_polys"))&
            & missing_polys = REPLACE_MISSING_POLYS
       if(have_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/&
            &boundary_treatment::lower_order")) then
          missing_polys = LOWER_ORDER
       end if

       leave_out_hermite_polynomials = .false.
       if(have_option(trim(T%option_path)//"/prognostic/&
            &spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter::Hermite_Weno/&
            &leave_out_hermite_polynomials")) &
            & leave_out_hermite_polynomials = .true.

       call allocate(T_limit, T%mesh, name="NewT") 
       T_limit%val = T%val
       
       limit_count = 0.0

       has_discontinuity_detector_field = has_scalar_field( &
            state, "DiscontinuityDetector")
       if(has_discontinuity_detector_field) then
          discontinuity_detector_field &
               => extract_scalar_field(state, "DiscontinuityDetector")
          discontinuity_detector_field%val = 0.0
       end if

       do ele = 1, element_count(T)
          
          call limit_slope_ele_hermite_weno(ele, T, T_limit, X, U)

       end do

       ewrite(3,*) 'Limit count = ',limit_count

       T%val = T_limit%val
       call deallocate(T_limit)

    case (LIMITER_VB)
       call limit_VB(state, T)

    case (LIMITER_FPN)
       call limit_fpn(state, T)      

    case default
       ewrite(-1,*) 'limiter = ', limiter
       FLAbort('no such limiter exists')
    end select

    ewrite(2,*) 'END subroutiune limit_slope_dg'

  end subroutine limit_slope_dg

  subroutine limit_slope_ele_dg(ele, T, X, T_limit)
    
    integer, intent(in) :: ele
    type(scalar_field), intent(inout) :: T
    type(vector_field), intent(in) :: X
    type(scalar_field), intent(inout), optional :: T_limit
    integer, dimension(:), pointer :: neigh, T_ele
    real, dimension(X%dim) :: ele_centre
    real :: ele_mean, miss_val
    integer :: ele_2, ni, face, face2, d, i, j, jj, miss
    real, dimension(X%dim, ele_loc(X,ele)) :: X_val, X_val_2
    real, dimension(ele_loc(T,ele)) :: T_val, T_val_2
    real, dimension(X%dim, ele_face_count(T,ele)) :: neigh_centre, face_centre
    real, dimension(ele_face_count(T,ele)) :: neigh_mean, face_mean
    real, dimension(mesh_dim(T)+1) :: b, new_val
    logical :: limit

    X_val=ele_val(X, ele)
    T_val=ele_val(T, ele)
    
    ele_centre=sum(X_val,2)/size(X_val,2)

    ele_mean=sum(T_val)/size(T_val)
    
    neigh=>ele_neigh(T, ele)

    limit=.false.

    searchloop: do ni=1,size(neigh)

       !----------------------------------------------------------------------
       ! Find the relevant faces.
       !----------------------------------------------------------------------
       
       ! These finding routines are outside the inner loop so as to allow
       ! for local stack variables of the right size in
       ! construct_add_diff_interface_dg.

       ele_2=neigh(ni)
    
       ! Note that although face is calculated on field U, it is in fact
       ! applicable to any field which shares the same mesh topology.
       face=ele_face(T, ele, ele_2)
       face2=ele_face(T, ele_2, ele)

       face_centre(:,ni) = sum(face_val(X,face),2)/size(face_val(X,face),2)

       face_mean(ni) = sum(face_val(T,face))/size(face_val(T,face))
       
       if (ele_2<=0) then
          ! External face.
          cycle
       end if

       X_val_2=ele_val(X, ele_2)
       T_val_2=ele_val(T, ele_2)

       neigh_centre(:,ni)=sum(X_val_2,2)/size(X_val_2,2)

       neigh_mean(ni)=sum(T_val_2)/size(T_val_2)
    
       if ((face_mean(ni)-ele_mean)*(face_mean(ni)-neigh_mean(ni))>0.0) then
          ! Limit if face_mean does not lie between ele_mean and neigh_mean
          limit=.true.
          
          if (face_mean(ni)>ele_mean) then
             face_mean(ni) = max(ele_mean, neigh_mean(ni))
          else
             face_mean(ni) = min(ele_mean, neigh_mean(ni))
          end if

       end if

    end do searchloop

    if (present(T_limit)) then
       T_ele=>ele_nodes(T_limit,ele)
       call set(T_limit, T_ele, ele_mean+0.0*T_ele)
    end if

    if (.not.limit) then
       return
    end if
    
    d=mesh_dim(T)
    new_val=ele_mean

    do miss=1,d+1
       
       ! If the missed side is a boundary, it is not possible to limit in
       ! this direction without violating the boundary condition.
       if (neigh(miss)<=0) cycle

       A=0.0
       b(1)=ele_mean

       do i=1, d+1
          ! Enforce preservation of the element mean value.
          A(1,i)=1.0/(d+1) 
          
          jj=1
          do j=1,d+1
             if (j==miss) cycle
             jj=jj+1
             
             if (i/=j) then
                A(jj,i)=1.0/d
             else
                b(jj)=face_mean(j)
             end if
             
          end do
          
       end do

       call invert(A)
       b=matmul(A,b)
       
       if (maxval(abs(b-ele_mean))>maxval(abs(new_val-ele_mean))) then
          !! The slope is larger than the current best guess.
          
          miss_val=0.0
          do ni=1, d+1
             if (ni==miss) cycle

             miss_val=miss_val+b(ni)/d
          end do
             
          if ((miss_val-ele_mean)*(miss_val-neigh_mean(miss))<=0.0) then
             ! The slope is legal.
             
             new_val=b

          end if

       end if

    end do
   
    ! Success or non-boundary failure.
    T_ele=>ele_nodes(T,ele)
    
    call set(T, T_ele, new_val)

    if (present(T_limit)) then
       T_ele=>ele_nodes(T_limit, ele)
       
       if (all(new_val==ele_mean)) then
          call set(T_limit, T_ele, 1.0+T_ele*0.0)
       else
          call set(T_limit, T_ele, -1.0+0.0*T_ele)
       end if

    end if

  end subroutine limit_slope_ele_dg

  subroutine cockburn_shu_setup(T,X)
    type(scalar_field), intent(inout) :: T
    type(vector_field), intent(in) :: X
    !
    logical :: do_setup
    integer :: cnt, d, i, j, ele

    do_setup = .false.
    if(.not.CSL_initialised) then
       CALL GetEventCounter(EVENT_ADAPTIVITY, csl_adapt_counter)
       do_setup = .true.
       CSL_initialised = .true.
    else 
       CALL GetEventCounter(EVENT_ADAPTIVITY, CNT)
       if(cnt.ne.csl_adapt_counter) then
          do_setup= .true.
          csl_adapt_counter = cnt
       end if
    end if

    if(do_setup) then

       if(associated(alpha)) then
          deallocate(alpha)
          alpha => null()
       end if
       if(associated(dx2)) then
          deallocate(dx2)
          dx2 => null()
       end if
       if(associated(A)) then
          deallocate(A)
          A => null()
       end if
       !!ATTENTION: This assumes that all elements have the same number of faces
       allocate(alpha(element_count(T),ele_face_count(T,1)&
            &,ele_face_count(T,1)))
       allocate(dx2(element_count(T),ele_face_count(T,1)))

       d=mesh_dim(T)
       allocate(A(d+1,d+1))

       ! Initialise A with the change from face centre values to node values.
       do i=1, size(A,1)
          do j=1,size(A,2)
             if (i==j) then
                A(i,j)=0.0
             else
                A(i,j)=1.0/d
             end if
          end do
       end do
       
       call invert(A)

       do ele = 1, element_count(T)
          call cockburn_shu_setup_ele(ele,T,X)
       end do

    end if    

  end subroutine cockburn_shu_setup

  subroutine cockburn_shu_setup_ele(ele, T, X)
    integer, intent(in) :: ele
    type(scalar_field), intent(inout) :: T
    type(vector_field), intent(in) :: X
    
    integer, dimension(:), pointer :: neigh, x_neigh
    real, dimension(X%dim) :: ele_centre, face_2_centre
    real :: max_alpha, min_alpha, neg_alpha
    integer :: ele_2, ni, nj, face, face_2, i, nk, ni_skip, info, nl
    real, dimension(X%dim, ele_loc(X,ele)) :: X_val, X_val_2
    real, dimension(X%dim, ele_face_count(T,ele)) :: neigh_centre,&
         & face_centre
    real, dimension(X%dim) :: alpha1, alpha2
    real, dimension(X%dim,X%dim) :: alphamat
    real, dimension(X%dim,X%dim+1) :: dx_f, dx_c
    integer, dimension(mesh_dim(T)) :: face_nodes

    X_val=ele_val(X, ele)
    
    ele_centre=sum(X_val,2)/size(X_val,2)
    
    neigh=>ele_neigh(T, ele)
    ! x_neigh/=t_neigh only on periodic boundaries.
    x_neigh=>ele_neigh(X, ele)

    searchloop: do ni=1,size(neigh)

       !----------------------------------------------------------------------
       ! Find the relevant faces.
       !----------------------------------------------------------------------
       ele_2=neigh(ni)
    
       ! Note that although face is calculated on field U, it is in fact
       ! applicable to any field which shares the same mesh topology.
       face=ele_face(T, ele, ele_2)
       face_nodes=face_local_nodes(T, face)

       face_centre(:,ni) = sum(X_val(:,face_nodes),2)/size(face_nodes)
       
       if (ele_2<=0) then
          ! External face.
          neigh_centre(:,ni)=face_centre(:,ni)
          cycle
       end if

       X_val_2=ele_val(X, ele_2)
       
       neigh_centre(:,ni)=sum(X_val_2,2)/size(X_val_2,2)
       if (ele_2/=x_neigh(ni)) then
          ! Periodic boundary case. We have to cook up the coordinate by
          ! adding vectors to the face from each side.
          face_2=ele_face(T, ele_2, ele)
          face_2_centre = &
               sum(face_val(X,face_2),2)/size(face_val(X,face_2),2)
          neigh_centre(:,ni)=face_centre(:,ni) + &
               (neigh_centre(:,ni) - face_2_centre)
       end if

    end do searchloop

    do ni = 1, size(neigh)
       dx_c(:,ni)=neigh_centre(:,ni)-ele_centre !Vectors from ni centres to
       !                                         !ele centre
       dx_f(:,ni)=face_centre(:,ni)-ele_centre !Vectors from ni face centres
                                              !to ele centre
    end do

    alpha_construction_loop: do ni = 1, size(neigh)
       !Loop for constructing Delta v(m_i,K_0) as described in C&S
       alphamat(:,1) = dx_c(:,ni)

       max_alpha = -1.0
       ni_skip = 0

       choosing_best_other_face_loop: do nj = 1, size(neigh)
          !Loop over the other faces to choose best one to use
          !for linear basis across face

          if(nj==ni) cycle
          
          !Construct a linear basis using all faces except for nj
          nl = 1
          do nk = 1, size(neigh)
             if(nk==nj.or.nk==ni) cycle
             nl = nl + 1
             alphamat(:,nl) = dx_c(:,nk)
          end do
          
          !Solve for basis coefficients alpha
          alpha2 = dx_f(:,ni)
          call solve(alphamat,alpha2,info)

          if((.not.any(alpha2<0.0)).and.alpha2(1)/norm2(alpha2)>max_alpha) &
               & then
             alpha1 = alpha2
             ni_skip = nj
             max_alpha = alpha2(1)/norm2(alpha2)
          end if

       end do choosing_best_other_face_loop

       if(max_alpha<0.0) then
          if(tolerate_negative_weights) then
             min_alpha = huge(0.0)
             ni_skip = 0
             choosing_best_other_face_neg_weights_loop: do nj = 1, size(neigh)
                !Loop over the other faces to choose best one to use
                !for linear basis across face
                
                if(nj==ni) cycle
                
                !Construct a linear basis using all faces except for nj
                nl = 1
                do nk = 1, size(neigh)
                   if(nk==nj.or.nk==ni) cycle
                   nl = nl + 1
                   alphamat(:,nl) = dx_c(:,nk)
                end do
                
                !Solve for basis coefficients alpha
                alpha2 = dx_f(:,ni)
                call solve(alphamat,alpha2,info)

                neg_alpha = 0.0
                do i = 1, size(alpha2)
                   if(alpha2(i)<0.0) then
                      neg_alpha = neg_alpha + alpha2(i)**2
                   end if
                end do
                neg_alpha = sqrt(neg_alpha)

                if(min_alpha>neg_alpha) then
                   alpha1 = alpha2
                   ni_skip = nj
                   min_alpha = neg_alpha
                end if
             end do choosing_best_other_face_neg_weights_loop
          else
             FLAbort('solving for alpha failed')
          end if
       end if
       
       alpha(ele,ni,:) = 0.0
       alpha(ele,ni,ni) = alpha1(1)
       nl = 1
       do nj = 1, size(neigh)
          if(nj==ni.or.nj==ni_skip) cycle
          nl = nl + 1
          alpha(ele,ni,nj) = alpha1(nl)
       end do

       dx2(ele,ni) = norm2(dx_c(:,ni))

    end do alpha_construction_loop

  end subroutine cockburn_shu_setup_ele

  subroutine limit_slope_ele_cockburn_shu(ele, T, X)
    !!< Slope limiter according to Cockburn and Shu (2001) 
    !!< http://dx.doi.org/10.1023/A:1012873910884
    integer, intent(in) :: ele
    type(scalar_field), intent(inout) :: T
    type(vector_field), intent(in) :: X
    
    integer, dimension(:), pointer :: neigh, x_neigh, T_ele
    real :: ele_mean
    real :: pos, neg
    integer :: ele_2, ni, face
    real, dimension(ele_loc(T,ele)) :: T_val, T_val_2
    real, dimension(ele_face_count(T,ele)) :: neigh_mean, face_mean
    real, dimension(mesh_dim(T)+1) :: delta_v
    real, dimension(mesh_dim(T)+1) :: Delta, new_val
    integer, dimension(mesh_dim(T)) :: face_nodes

    T_val=ele_val(T, ele)
    
    ele_mean=sum(T_val)/size(T_val)
    
    neigh=>ele_neigh(T, ele)
    ! x_neigh/=t_neigh only on periodic boundaries.
    x_neigh=>ele_neigh(X, ele)

    searchloop: do ni=1,size(neigh)

       !----------------------------------------------------------------------
       ! Find the relevant faces.
       !----------------------------------------------------------------------
       ele_2=neigh(ni)
    
       ! Note that although face is calculated on field U, it is in fact
       ! applicable to any field which shares the same mesh topology.
       face=ele_face(T, ele, ele_2)
       face_nodes=face_local_nodes(T, face)
       
       face_mean(ni) = sum(T_val(face_nodes))/size(face_nodes)
       
       if (ele_2<=0) then
          ! External face.
          neigh_mean(ni)=face_mean(ni)
          cycle
       end if

       T_val_2=ele_val(T, ele_2)

       neigh_mean(ni)=sum(T_val_2)/size(T_val_2)

    end do searchloop

    delta_v = matmul(alpha(ele,:,:),neigh_mean-ele_mean)

    delta_loop: do ni=1,size(neigh)

       Delta(ni)=TVB_minmod(face_mean(ni)-ele_mean, &
            Limit_factor*delta_v(ni), dx2(ele,ni))

    end do delta_loop

    if (abs(sum(Delta))>1000.0*epsilon(0.0)) then
       ! Coefficients do not sum to 0.0

       pos=sum(max(0.0, Delta))
       neg=sum(max(0.0, -Delta))
       
       Delta = min(1.0,neg/pos)*max(0.0,Delta) &
            -min(1.0,pos/neg)*max(0.0,-Delta)
       
    end if

    new_val=matmul(A,Delta+ele_mean)
    
    ! Success or non-boundary failure.
    T_ele=>ele_nodes(T,ele)
    
    call set(T, T_ele, new_val)

  end subroutine limit_slope_ele_cockburn_shu

  !11:25 <Guest54276>     do ele_A=1,ele_count(old_position)
  !11:25 <Guest54276>       call local_coords_matrix(old_position, ele_A, 
  !                   inversion_matrices_A(:, :, ele_A))
  !11:25 <Guest54276>     end do
  
  !subroutine local_coords_matrix(positions, ele, mat)
  !inputs global coordinates
  !outputs local coordinates

  subroutine limit_slope_ele_hermite_weno(ele, T, T_limit, X, U)
    !!< Hermite Weno Slope limiter
    integer, intent(in) :: ele
    type(scalar_field), intent(in) :: T
    type(scalar_field), intent(inout) :: T_limit
    type(vector_field), intent(in) :: X, U

    integer, dimension(:), pointer :: neigh, x_neigh, T_ele
    real :: ele_mean, ele_mean_2
    real, dimension(ele_face_count(T,ele)) :: ele_means
    real :: residual
    integer :: ele_2, ni, nj, face, face_2,i, nk, info, nl
    integer :: l_face, l_face_2
    real, dimension(ele_loc(T,ele)) :: T_val, T_val_2
    real, dimension(face_loc(T,1)) :: T_val_face
    real, dimension(face_ngi(T,1)) :: T_face_quad
    real, dimension(ele_face_count(T,ele),ele_loc(X,ele)) :: T_vals
    real, dimension(ele_face_count(T,ele),X%dim, ele_loc(X,ele)) :: X_vals
    real, dimension(X%dim, ele_loc(X,ele)) :: X_val
    real, dimension(ele_face_count(T,ele)) :: neigh_mean, face_mean
    real, dimension(ele_loc(T,ele)) :: new_val
    integer, dimension(mesh_dim(T)) :: face_nodes
    logical :: limit_slope
    real, dimension(ele_loc(T, ele), ele_ngi(T, ele), &
         &mesh_dim(T)) :: du_t
    real, dimension(ele_ngi(T,ele)) :: detwei
    real, dimension(ele_ngi(T,ele)) :: p_quad, T_quad
    real, dimension(1+2*ele_loc(T,ele),ele_loc(T,ele)) :: Polys
    real, dimension(ele_loc(T,ele)*2+1) :: Polys_o, Polys_w
    real, dimension(mesh_dim(T),ele_ngi(T,ele)) :: dp_quad
    logical, dimension(ele_face_count(T,ele)) :: boundaries
    logical, dimension(ele_face_count(T,ele)) :: construct_Lagrange
    type(element_type), pointer :: shape_T
    real, dimension(ele_loc(T,ele),ele_loc(T,ele)) :: Imat
    real, dimension(ele_loc(T,ele)) :: Irhs
    real, dimension(ele_loc(X,ele)) :: local_coords
    integer, dimension(face_loc(T,1)) :: l_face_list,l_face_list_2

    real :: Discontinuity_indicator, inflow_integral, h
    real, dimension(ele_loc(T,ele)) :: ones
    integer :: discontinuity_option 
    real :: face_max, face_min

    if(debugging) then
       ewrite(2,*) 'Limit_slope_Hermite_weno_ele'
    end if

    limit_slope = .false.

    boundaries = .false.
    construct_Lagrange = .true.

    T_val=ele_val(T, ele)
    X_val=ele_val(X, ele)

    ele_mean=sum(T_val)/size(T_val)

    neigh=>ele_neigh(T, ele)
    ! x_neigh/=t_neigh only on periodic boundaries.
    x_neigh=>ele_neigh(X, ele)

    discontinuity_option = 2

    select case(discontinuity_option)

    case (1)
       !=========================================================
       !Discontinuity detector using TVB condition
       !Checks solution on each face is between mean values of 
       !ele and ele_2
       !=========================================================
       do ni=1,size(neigh)

          !--------------------------------------------------------------------
          ! Find the relevant faces.
          !--------------------------------------------------------------------
          ele_2=neigh(ni)

          if(ele_2<0) cycle
          
          T_val_2 = ele_val(T,ele_2)
          face=ele_face(T, ele, ele_2)
          T_val_face = face_val(T, face)
          T_face_quad = face_val_at_quad(T,face)
          !face_max = maxval(T_val_face)
          !face_min = minval(T_val_face)
          face_max = maxval(T_face_quad)
          face_min = minval(T_face_quad)
          ele_mean_2 = sum(T_val_2)/size(T_val_2)

          !ewrite(3,*) T_face_quad
          !ewrite(3,*) T_val
          !ewrite(3,*) T_val_2

          if(face_max>max(ele_mean,ele_mean_2)+disc_tol) limit_slope = .true.
          if(face_min<min(ele_mean,ele_mean_2)-disc_tol) limit_slope = .true.

          if(has_discontinuity_detector_field) then
             if(limit_slope) then
                ewrite(3,*) 'cjc limit_slope', ele
                ones = 1.0
                T_ele=>ele_nodes(Discontinuity_Detector_field,ele)
                call set(Discontinuity_detector_field,T_ele,ones)
             end if
          end if
       end do

    case (2)

       !=================================================================
       !DISCONTINUITY INDICATOR,
       !from http://www.gce.ucl.ac.be/~remacle/pdf/detect.pdf
       !We compute the jump of the solution on upwind boundaries
       !=================================================================

       !Initial value of integral of jump of solution on inflow boundaries
       Discontinuity_indicator = 0.0
       !Initial value of inflow area/length
       Inflow_integral = 0.0
       !We are going to increment these

       do ni=1,size(neigh)

          !--------------------------------------------------------------------
          ! Find the relevant faces.
          !--------------------------------------------------------------------
          ele_2=neigh(ni)

          if(ele_2<0) cycle

          face=ele_face(T, ele, ele_2)
          face_2=ele_face(T, ele_2, ele)

          call Discontinuity_indicator_face(Discontinuity_indicator, &
               & Inflow_integral, &
               & U,T,X,ele,face,face_2)

       end do

       discontinuity_indicator = abs(discontinuity_indicator)
       inflow_integral = abs(inflow_integral)

       !Compute h
       h = get_H(X_val)

       !Get max norm in element of T
       T_quad = ele_val_at_quad(T,ele)

       if(Discontinuity_Indicator>disc_tol*Inflow_integral&
            &*maxval(abs(T_quad))*h) limit_slope = .true.

       if(has_discontinuity_detector_field) then
          ones = 1.0
          T_ele=>ele_nodes(Discontinuity_Detector_field,ele)
          call set(Discontinuity_detector_field,T_ele&
               &,Discontinuity_Indicator*ones/inflow_integral/&
               maxval(abs(T_quad)+limit_tol)/h)
       end if

    case default
       FLExit('no such discontinuity option')
    end select

    if(limit_slope) then
       ewrite(2,*) 'cjc: limiting slope'
       limit_count = limit_count + 1

       !Apply HWENO limiter

       setuploop: do ni=1,size(neigh)

          ele_2=neigh(ni)

          if (ele_2<=0) then
             ! External face.
             neigh_mean(ni)=face_mean(ni)
             boundaries(ni) = .true.

             do nj = 1, size(neigh)
                if(ni==nj) cycle
                construct_Lagrange(nj) = .false.
             end do
             cycle

          end if

          ! Note that although face is calculated on field U, it is in fact
          ! applicable to any field which shares the same mesh topology.
          face=ele_face(T, ele, ele_2)
          face_2=ele_face(T, ele_2, ele)
          face_nodes=face_local_nodes(T, face)

          face_mean(ni) = sum(T_val(face_nodes))/size(face_nodes)

          T_val_2=ele_val(T, ele_2)

          T_vals(ni,:) = T_val_2
          X_vals(ni,:,:) = ele_val(X, ele_2)

          neigh_mean(ni)=sum(T_val_2)/size(T_val_2)
          ele_means(ni) = neigh_mean(ni)

       end do setuploop

       if(any(boundaries).and.(missing_polys==LOWER_ORDER)) then
          !On boundary, with this option, just project to p(n-1)
          !We have only coded P1 so this projects to P0

          new_val = sum(T_val)/size(T_val)

       else

          Polys = 0.0

          if(debugging) then
             ewrite(2,*) 'Limiting slope.'
          end if

          !We store transformations (du_t, detwei) in the following way:
          ! 1:size(neigh) : transformations for neighbouring element
          ! size(neigh)+1 : transformations for this element

          shape_T=>ele_shape(T,ele)

          !Construct transformations for this element
          call transform_to_physical(X, ele,&
               & shape_T , dshape=du_t, detwei=detwei)

          !Polynomials are stored in the following way:
          ! i = 1:size(neigh) : Lagrange polynomials obtained by
          !                     missing out the i-th neighbour
          ! i = size(neigh)+1 : The existing polynomial representation
          ! j = size(neigh)+1 + i, i = 1:size(neigh) : The function with same mean
          !                                            as existing polynomial
          !                                            with slope taken from 
          !                                            i-th neighbour
          ! The latter representations are Hermite polynomials using 
          ! gradient information

          !Construct Lagrange polys
          !Fails if non-flat elements
          LagrangeP_loop: do ni = 1, size(neigh)
             if(.not.construct_Lagrange(ni)) then
                Polys(ni,:) = T_val
             else
                nl = 0
                do nj = 1, size(neigh)
                   if(nj==ni) cycle
                   nl = nl + 1
                   !
                   !This row requires that the mean of the polynomial over
                   !neighbour element nj be equal to the mean of the unlimited
                   !solution in that element.
                   ! 
                   ! This is done by computing the local coordinates of the
                   ! centre of the neighbour element (only works for P1) which
                   ! gives you the coefficients of the local expansion which
                   ! give you the polynomial value at that point which is then
                   ! required to be equal to the current element mean.
                   !
                   do nk = 1, size(neigh)
                      Imat(nl,:) = local_coords_interpolation(X,ele&
                           &,sum(X_vals(nj,:,:),2)/size(X_vals,3))
                   end do
                   Irhs(nl) = ele_means(nj)
                end do
                !Last column sets the mean value
                Imat(size(neigh),:) = 1.0/size(Imat,2)
                Irhs(size(neigh)) = ele_mean

                !Solve for the Polynomial
                call solve(Imat,Irhs,info)
                Polys(ni,:) = Irhs

                if(debugging) then
                   !Do some checking
                   !Compute the polynomial at the quad points
                   !Check polynomial has mean value ele_mean in this element
                   if(abs(sum(Polys(ni,:)-ele_mean))>1.0e-5) then
                      FLAbort('failed to get the correct mean value in this element')
                   end if
                   !Check polynomial has mean value ele_means in other two elements
                   do nj = 1, size(neigh)
                      if(nj==ni) cycle
                      local_coords = local_coords_interpolation(X,ele&
                           &,sum(X_vals(nj,:,:),2)/size(X_vals,3))
                      residual = sum(Polys(ni,:)*local_coords) - ele_means(nj)
                      if(abs(residual)>1.0e-5) then
                         FLAbort('failed to get the correct mean value in neighbour')
                      end if
                   end do
                end if

             end if
          end do LagrangeP_loop

          !Construct Hermite polys
          !Fails if non-flat elements

          !Original polynomial
          Polys(size(neigh)+1,:) = T_val
          HermiteP_loop: do ni = 1, size(neigh)
             !Mean of original values with slope of neighbouring element
             nk = size(neigh)+1+ni

             ele_2=neigh(ni)

             if(ele_2<0) then
                Polys(nk,:) = T_val
                cycle
             end if

             !face number in ele
             face=ele_face(T, ele, ele_2)
             !face number in ele_2
             face_2=ele_face(T, ele_2, ele)

             !local face number in ele
             l_face = local_face_number(T, face)
             !local face number in ele_2
             l_face_2 = local_face_number(T, face_2)

             !Local face list in ele
             l_face_list = face_local_nodes(T, face)
             !Local face list in ele_2
             l_face_list_2 = face_local_nodes(T, face_2)

             !T values in ele_2
             T_val_2=ele_val(T, ele_2)

             !First we "continue" the polynomial in ele_2 into ele

             !Polynomial takes same values on shared nodes
             Polys(nk,l_face_list) = &
                  & T_val_2(l_face_list_2)

             !Compute local coordinates (relative to ele)
             !of vertex in ele_2 which is opposite the face
             local_coords = local_coords_interpolation(X,ele&
                  &,X_vals(ni,:,l_face_2))

             !Solve 1D linear system to get polynomial value
             !from
             !T_val_2(l_face_2) = sum(Poly*local_coords)
             !we have already computed the values on the face
             !so we rearrange.
             Polys(nk,l_face) = (T_val_2(l_face_2) - &
                  & sum(Polys(nk,l_face_list)*local_coords(l_face_list)))/ &
                  local_coords(l_face)

             !ADD SOME DEBUGGING TESTS

             !Second we adjust the mean so that it is the same as the 
             !mean of the current solution in ele
             Polys(nk,:) = Polys(nk,:)- &
                  & sum(Polys(nk,:))/size(T_val) + &
                  & sum(T_val)/size(T_val)         

          end do HermiteP_loop

          if(debugging) then
             ewrite(2,*) 'Dumping polynomials'
             do ni = 1, size(neigh)*2 + 1
                ewrite(2,*) Polys(ni,:)
             end do
          end if

          !Compute oscillatory indicators

          do ni = 1, size(neigh)*2 + 1
             !construct the ni-th polynomial at the quadrature points
             P_quad=matmul(Polys(ni,:), shape_T%n)
             !construct the gradient of the ni-th polynomial at the quad points
             do i = 1, mesh_dim(T)
                dP_quad(i,:) = &
                     &matmul(Polys(ni,:), du_t(:,:,i))
             end do

             !construct the oscillator index of the ni-th polynomial
             Polys_o(ni) = 0.0
             do i = 1, mesh_dim(T)
                Polys_o(ni) = Polys_o(ni) + &
                     & sum(detwei*dP_quad(i,:)**2)
             end do
             Polys_o(ni) = Polys_o(ni)/&
                  &sum(detwei*(eps_o + P_quad)**2)
          end do

          if(debugging) then
             ewrite(2,*) 'Dumping oscillatory indicators'
             ewrite(2,*) Polys_o
          end if

          !Compute weights
          do ni = 1, size(neigh)*2 + 1
             Polys_w(ni) = (eps_w + Polys_o(ni))**(-gam0)
          end do

          if(missing_polys==IGNORE_MISSING_POLYS.and.any(boundaries)) then
             do ni = 1, size(neigh)
                if(boundaries(ni)) then
                   do nj = 1, size(neigh)
                      if(ni==nj) cycle
                      Polys_w(nj) = 0.0
                   end do
                   Polys_w(size(neigh)+1+ni) = 0.0
                end if
             end do
          end if

          if(leave_out_hermite_polynomials) then
             Polys_w(size(neigh)+1:2*size(neigh)+1) = 0.0
          end if

          Polys_w = Polys_w/sum(Polys_w)

          if(debugging) then
             ewrite(2,*) 'Dumping weights'
             ewrite(2,*) Polys_w
          end if

          new_val = 0.
          do ni = 1, size(neigh)*2 + 1
             new_val = new_val + Polys_w(ni)*Polys(ni,:)
          end do

          if(debugging) then
             ewrite(2,*) 'new val is'
             ewrite(2,*) new_val

             ewrite(2,*) 'old slope was'
             do i = 1, mesh_dim(T)
                ewrite(2,*) maxval(matmul(T_val, du_t(:,:,i)))
             end do

             ewrite(2,*) 'new slope is'
             do i = 1, mesh_dim(T)
                ewrite(2,*) maxval(matmul(new_val, du_t(:,:,i)))
             end do
          end if
       end if

       T_ele=>ele_nodes(T,ele)
       call set(T_limit, T_ele, new_val)

    end if

  end subroutine limit_slope_ele_hermite_weno

  function TVB_minmod(a1,a2, dx)
    real :: TVB_minmod
    real, intent(in) :: a1, a2, dx

    if (abs(a1)<TVB_factor*dx**2) then
       TVB_minmod=a1
    else if (abs(a1)<abs(a2)) then
       TVB_minmod=a1
    else
       TVB_minmod=a2
    end if

  end function TVB_minmod

  subroutine limit_slope_ele_dg_1d(ele, T, X)
    
    integer, intent(in) :: ele
    type(scalar_field), intent(inout) :: T
    type(vector_field), intent(in) :: X
    
    integer, dimension(:), pointer :: neigh, T_ele
    real, dimension(mesh_dim(X)) :: ele_centre, ele_2_centre
    real :: ele_mean, ele_2_mean, dx
    real :: ele_slope, old_ele_slope, ele_2_slope
    integer :: ele_2, ni
    real, dimension(mesh_dim(X), ele_loc(X,ele)) :: X_val, X_val_2
    real, dimension(ele_loc(T,ele)) :: T_val, T_val_2
    
    X_val=ele_val(X, ele)
    T_val=ele_val(T, ele)
    

    ele_centre=sum(X_val,2)/size(X_val,2)

    ele_mean=sum(T_val)/size(T_val)
    
    dx=X_val(1,2)-X_val(1,1)
    
    ele_slope=(T_val(2)-T_val(1))/dx
    
    old_ele_slope=ele_slope

    neigh=>ele_neigh(T, ele)

    neighbourloop: do ni=1,size(neigh)

       !----------------------------------------------------------------------
       ! Find the relevant faces.
       !----------------------------------------------------------------------
       
       ! These finding routines are outside the inner loop so as to allow
       ! for local stack variables of the right size in
       ! construct_add_diff_interface_dg.

       ele_2=neigh(ni)
           
       if (ele_2<=0) then
          ! External face.
          cycle
       end if

       X_val_2=ele_val(X, ele_2)
       T_val_2=ele_val(T, ele_2)

       ele_2_centre=sum(X_val_2,2)/size(X_val_2,2)

       ele_2_mean=sum(T_val_2)/size(T_val_2)
    
       ele_2_slope=(ele_2_mean-ele_mean)/sum(ele_2_centre-ele_centre)
       
       if (ele_slope*ele_2_slope<0.0) then
          ! Slope sign changes
          ele_slope=0.0
          exit neighbourloop
       end if

       ele_slope=sign(min(abs(ele_slope),abs(ele_2_slope)), ele_slope)

    end do neighbourloop

    if (old_ele_slope/=ele_slope) then
       
       ! Remove high order stuff here.
       T_ele=>ele_nodes(T,ele)

       call set(T, T_ele(1), ele_mean-0.5*dx*ele_slope)
       call set(T, T_ele(2), ele_mean+0.5*dx*ele_slope)

    end if

  end subroutine limit_slope_ele_dg_1d

  subroutine Discontinuity_indicator_face(Discontinuity_indicator, &
       & Inflow_integral, &
       & U,T,X,ele,face,face_2)
    real, intent(inout) :: Discontinuity_indicator, Inflow_integral
    type(vector_field), intent(in) :: U,X
    type(scalar_field), intent(in) :: T
    integer, intent(in) :: ele, face, face_2
    !
    real, dimension(face_ngi(T,face)) :: detwei
    real, dimension(mesh_dim(T),face_ngi(T,face)) :: normal
    real, dimension(mesh_dim(U),face_ngi(U,face)) :: U_flux
    integer, dimension(face_ngi(U,face)) :: inflow
    logical :: use_mean_inflow
    !stuff for local calculations
    real, dimension(U%dim,ele_loc(X,ele)) :: X_val
    real, dimension(face_loc(T,face)) :: T_face_val
    real, dimension(U%dim,face_loc(X,face)) :: X_face_val
    real, dimension(U%dim) :: centroid2face, normal_vec, Vec1, Vec2
    real, dimension(U%dim,face_loc(U,face)) :: U_face_val
    real :: Area

    X_val = ele_val(X,ele)
    x_face_val = face_val(X,face)

    use_mean_inflow = .false.

    if(use_mean_inflow) then

       !We only compute on mean inflow boundaries
       !This means we can avoid transforming to physical
       !only works on flat elements

       !compute normals
       select case (U%dim)
       case (2)
          Vec1 = x_face_val(:,1) - x_face_val(:,2)
          normal_vec(1) = -Vec1(2)
          normal_vec(2) = Vec1(1)
          Area = sqrt(sum(Vec1**2))
       case (3)
          Vec1 = x_face_val(:,1) - x_face_val(:,2)
          Vec2 = x_face_val(:,1) - x_face_val(:,3)
          normal_vec(1) = Vec1(2)*Vec2(3)-Vec1(3)*Vec2(2)
          normal_vec(2) = -Vec1(1)*Vec2(3)+Vec1(3)*Vec2(1)
          normal_vec(3) = Vec1(1)*Vec2(2)-Vec1(2)*Vec2(1)
          Area = 0.5*sqrt(sum(normal_vec**2))
       case default
          FLExit('cant handle that case - that dimension is not supported')
       end select
       normal_vec = normal_vec / sqrt(sum(normal_vec**2))
       centroid2face = sum(X_face_val,2)/size(X_face_val,2) - &
            & sum(X_val,2)/size(X_val,2)
       if(sum(normal_vec*centroid2face)<0.0) normal_vec = -normal_vec

       !Check if have mean inflow on this faec
       U_face_val = face_val(U,face)
       if(sum(sum(U_face_val,2)/size(U_face_val,2)*normal_vec)<1.0e-15) then

          T_face_val = face_val(T,face) - face_val(T,face_2)

          Discontinuity_indicator = Discontinuity_indicator + &
               sum(T_face_val)/size(T_face_val)*Area

          Inflow_integral = Inflow_integral + Area
       end if

    else

       call transform_facet_to_physical( X, face,&
            & detwei_f=detwei,normal=normal)
       
       U_flux = 0.5*(face_val_at_quad(U,face)+ &
            & face_val_at_quad(U,face_2))
       
       !We only compute on inflow boundaries    
       inflow = merge(1.0,0.0,sum(U_flux*normal,1)<0.0)

       Discontinuity_indicator = &
            & Discontinuity_indicator + &
            abs(sum( (face_val_at_quad(T, face) &
            - face_val_at_quad(T,face_2))*detwei*inflow ))

       Area = abs(sum(detwei*inflow))
       Inflow_integral = Inflow_integral + Area
    end if

  end subroutine Discontinuity_indicator_face

  function get_H(X) result (h)
    real, dimension(:,:), intent(in) :: X
    real :: h
    !
    integer :: i,j,dim
    real :: a,b,c

    dim = size(X,1)

    select case(dim)
    case (1)
       !Just take the difference
       h = abs(X(1,1)-X(1,2))
    case (2)
       !Circumradius
       a = sqrt(sum((X(:,2)-X(:,1))**2))
       b = sqrt(sum((X(:,3)-X(:,2))**2))
       c = sqrt(sum((X(:,1)-X(:,3))**2))
       h = a*b*c/sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c))
    case (3)
       !This should be circumradius too but I didn't code it
       h = 0.0
       do i = 1, size(X,2)
          do j = 2, size(X,2)
             h = max(h,sqrt(sum( (X(:,i)-X(:,j))**2 )))
          end do
       end do
    case default
       FLExit('dont know that dimension.')
    end select
    
  end function get_H

  subroutine limit_vb(state, t)
    !Vertex-based (not Victoria Bitter) limiter from
    !Kuzmin, J. Comp. Appl. Math., 2010
    ! doi:10.1016/j.cam.2009.05.028
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: t
    !
    ! This is the limited version of the field, we have to make a copy
    type(scalar_field) :: T_limit, T_max, T_min
    type(mesh_type), pointer :: vertex_mesh
    ! counters
    integer :: ele, node
    ! local numbers
    integer, dimension(:), pointer :: T_ele
    ! gradient scaling factor
    real :: alpha
    ! local field values
    real, dimension(ele_loc(T,1)) :: T_val, T_val_slope, T_val_min,T_val_max
    real :: Tbar

    if (.not. element_degree(T%mesh, 1)==1 .or. continuity(T%mesh)>=0) then
      FLExit("The vertex based slope limiter only works for P1DG fields.")
    end if
    
    ! Allocate copy of field
    call allocate(T_limit, T%mesh,trim(T%name)//"Limited")
    call set(T_limit, T)
    
    ! returns linear version of T%mesh (if T%mesh is periodic, so is vertex_mesh)
    call find_linear_parent_mesh(state, T%mesh, vertex_mesh)

    call allocate(T_max, vertex_mesh, trim(T%name)//"LimitMax")
    call allocate(T_min, vertex_mesh, trim(T%name)//"LimitMin")
 
    call set(T_max, -huge(0.0))
    call set(T_min, huge(0.0))

    ! for each vertex in the mesh store the min and max values of the P1DG nodes directly surrounding it
    do ele = 1, ele_count(T)
       T_ele => ele_nodes(T,ele)
       T_val = ele_val(T,ele)
       Tbar = sum(T_val)/size(T_val)
       ! we assume here T is P1DG and vertex_mesh is linear
       assert( size(T_ele)==ele_loc(vertex_mesh,ele) )
       
       ! do maxes
       T_val_max = ele_val(T_max,ele)
       do node = 1, size(T_val)
          T_val_max(node) = max(T_val_max(node), Tbar)
       end do
       call set(T_max, ele_nodes(T_max, ele), T_val_max)
       
       ! do mins
       T_val_min = ele_val(T_min,ele)
       do node = 1, size(T_val)
          T_val_min(node) = min(T_val_min(node), Tbar)
       end do
       call set(T_min, ele_nodes(T_min,ele), T_val_min)
    end do

    ! now for each P1DG node make sure the field value is between the recorded vertex min and max
    ! this is done without changing the element average (Tbar)
    do ele = 1, ele_count(T)
       !Set slope factor to 1
       alpha = 1.
       !Get local node lists
       T_ele=>ele_nodes(T,ele)
       
       T_val = ele_val(T,ele)
       Tbar = sum(T_val)/size(T_val)
       T_val_slope = T_val - Tbar
       T_val_max = ele_val(T_max,ele)
       T_val_min = ele_val(T_min,ele)
       
       !loop over nodes, adjust alpha
       do node = 1, size(T_val)
         !check whether to use max or min, and avoid floating point algebra errors due to round-off and underflow
         if(T_val(node)>Tbar*(1.0+sign(1.0e-12,Tbar)) .and. T_val(node)-Tbar > tiny(0.0)*1e10) then
           alpha = min(alpha,(T_val_max(node)-Tbar)/(T_val(node)-Tbar))
         else if(T_val(node)<Tbar*(1.0-sign(1.0e-12,Tbar)) .and. T_val(node)-Tbar < -tiny(0.0)*1e10) then
           alpha = min(alpha,(T_val_min(node)-Tbar)/(T_val(node)-Tbar))
         end if
       end do

       call set(T_limit, T_ele, Tbar + alpha*T_val_slope)
    end do


    !Deallocate copy of field
    call set(T, T_limit)
    call halo_update(T)
    call deallocate(T_limit)
    call deallocate(T_max)
    call deallocate(T_min)

  end subroutine limit_vb

  subroutine limit_fpn(state, t)

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: t
    
    type(scalar_field), pointer :: limiting_t, lumped_mass
    type(scalar_field) :: lowerbound, upperbound, inverse_lumped_mass
    type(csr_matrix), pointer :: mass

    type(csr_sparsity), pointer :: eelist
    integer :: ele, i, j, k, row, column
    integer :: rows, columns
    real :: node_max, node_min, extra_val, extra_val2

    integer, dimension(:), pointer :: nodelist, faces, neighbouring_ele_nodes
    integer, dimension(:), allocatable :: face_nodes, neighbouring_nodes
    integer :: neighbouring_face, neighbouring_ele
    logical, save :: first=.true.
    logical :: midpoint, extrapolate, pre_dist_mass

    real :: beta=1.0, mean_val
    type(vector_field), pointer :: position
    type(vector_field) :: dg_position
    real, dimension(:), allocatable :: e_vec_1
    real, dimension(:,:,:), allocatable :: dt_t
    real, dimension(:,:), allocatable :: grad_t
    real :: grad, e_dist

    real, dimension(ele_loc(t,1)) :: weight, tracer_val
    logical, dimension(ele_loc(t,1)) :: nweight, pweight
    real :: nodeval, nodemin, nodemax, adjust

    integer, dimension(:,:,:), allocatable, save :: nodes_array
!     real, dimension(2,4) :: local_values
!     real, dimension(2) :: line_max, line_min
    real, dimension(:,:), allocatable :: local_values
    real, dimension(:), allocatable :: line_max, line_min
    integer :: node, adjacent_node, local_face

    integer, dimension(face_loc(t,1)) :: fnodes, neighbouring_face_nodes

    type(vector_field), pointer :: u
    real :: mod_u, tol
    real, dimension(2) :: u_node, e_vec_2, e_vec_3
    logical :: upwind=.false.
    real :: cor1, cor2, cor3, cor4
    integer :: problem_dimension, values

    call get_option('/geometry/dimension', problem_dimension)

    rows=problem_dimension ! The number of 'lines' to look along
    columns=3 ! The number of nodes on each line. Should always be three
    if (upwind) then
      values=columns+1
    else
      values=columns
    end if
    tol=0.25

    allocate(local_values(rows,values),line_max(rows),line_min(rows))

    midpoint=have_option(trim(t%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/slope_limiter::FPN/mid-point_scheme")
    if (midpoint) then
      call get_option(trim(t%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/slope_limiter::FPN/mid-point_scheme/beta", beta, default=1.0)
      extrapolate=have_option(trim(t%option_path)//"/prognostic/spatial_discretisation/&
           &discontinuous_galerkin/slope_limiter::FPN/mid-point_scheme/extrapolate")
    end if

!     pre_dist_mass=have_option(trim(t%option_path)//"/prognostic/spatial_discretisation/&
!          &discontinuous_galerkin/slope_limiter::FPN/pre_distribute_mass")

    mass => get_mass_matrix(state, t%mesh)
    lumped_mass => get_lumped_mass(state, t%mesh)
    call allocate(inverse_lumped_mass, lumped_mass%mesh, "InverseLumpedMass")
    inverse_lumped_mass%val = 1.0/lumped_mass%val

    limiting_t => extract_scalar_field(state, trim(t%name))

!     eelist => extract_eelist(t%mesh)

    call allocate(lowerbound, t%mesh, "LowerBound")
    call allocate(upperbound, t%mesh, "UpperBound")
    call zero(lowerbound); call zero(upperbound)

    allocate (neighbouring_nodes(ele_loc(limiting_t,1)))

!     allocate (face_nodes(face_loc(limiting_t,1)), neighbouring_nodes(face_loc(limiting_t,1)))

    if (extrapolate) then
      position => extract_vector_field(state, "Coordinate")
      call allocate(dg_position, position%dim, t%mesh, name="DG_Coordinate")
      call remap_field(position, dg_position)
      allocate (e_vec_1(position%dim),dt_t(ele_loc(limiting_t, 1), ele_ngi(limiting_t, 1), mesh_dim(limiting_t)))
      allocate (grad_t(mesh_dim(limiting_t), limiting_t%mesh%shape%ngi))
    end if

    if (upwind) then
      u => extract_vector_field(state, "Velocity")
    end if

    ! Loop to construct an array containing the global node numbers required to compute the limiting values
    ! at a node i. Only evaluated on the first timestep (and after every adapt for adaptive runs).
    if (first) then
      allocate (nodes_array(node_count(t),rows,columns))
      first=.false.
      do node=1,node_count(limiting_t)
        ele=node_ele(limiting_t, node)
        nodelist => ele_nodes(limiting_t, ele)
        faces => ele_faces(limiting_t, ele)
        row=0
        do i=1,size(nodelist)
          if (nodelist(i)==node) cycle
          row=row+1
          fnodes=face_global_nodes(limiting_t,faces(i))
          do j=1,size(fnodes)
            if (fnodes(j)==node) adjacent_node=j
          end do
          neighbouring_face = face_neigh(limiting_t, faces(i))
          neighbouring_face_nodes = face_global_nodes(limiting_t,neighbouring_face)
!           secnd_val=neighbouring_face_nodes(adjacent_node) ! 2nd node we want
          local_face=local_face_number(limiting_t,neighbouring_face)
          neighbouring_ele = face_ele(limiting_t, neighbouring_face)
          neighbouring_nodes = ele_nodes(limiting_t, neighbouring_ele)
!           thrid_val=neighbouring_nodes(local_face)
          nodes_array(node,row,1)=nodelist(i)
          nodes_array(node,row,2)=neighbouring_face_nodes(adjacent_node)
          nodes_array(node,row,3)=neighbouring_nodes(local_face)
        end do
      end do
    end if

    ! Loop through the nodes and calculate the bounds for each node
    do node=1,node_count(limiting_t)
      ! Calculate the av. value of the tracer within the element
      ele=node_ele(limiting_t, node)
      nodelist => ele_nodes(limiting_t, ele)
      do i=1, size(nodelist)
        tracer_val(i)=node_val(limiting_t,nodelist(i))
      end do
      mean_val=sum(tracer_val)/float(size(nodelist))
      ! Get the values needed for calculating the bounds
      do row=1,rows
        do column=1,columns
          local_values(row,column)=node_val(limiting_t, nodes_array(node,row,column))
        end do
        ! Adjust values depending on options
        if (midpoint.and.(.not.extrapolate)) then
          local_values(row,3) = (1.0-beta)*local_values(row,2)+beta*local_values(row,3)
        else if (midpoint.and.extrapolate) then
          local_values(row,3) = (1.0-beta)*local_values(row,2)+beta*local_values(row,3)
          ! Extrapolate using the gradients of the neighbouring element to form our extra value
                
          ! 1st, work out the direction in which we want to extrapolate, e_vec_1
          e_vec_1=node_val(dg_position,node)-node_val(dg_position,nodes_array(node,row,1))
          ! Work out the distance to exprapolate
          e_dist=sqrt(sum(e_vec_1(:)**2))
          ! Turn this into a unit vector
          e_vec_1=e_vec_1/e_dist

          call transform_to_physical(dg_position, node_ele(limiting_t, nodes_array(node,row,2)), &
                                     ele_shape(limiting_t,node_ele(limiting_t, nodes_array(node,row,2))), dshape=dt_t)

          grad_t=ele_grad_at_quad(limiting_t, node_ele(limiting_t, node), dt_t)

          ! Calculate the gradient in the desired direction
          ! Note that grad_t will be the same at all gauss points in the linear element case
          grad=dot_product(grad_t(1,:),e_vec_1)

          local_values(row,4) = local_values(row,2)+beta*grad*e_dist 
        end if
        if (upwind) then
          u_node=node_val(u, node)
          mod_u=sqrt(dot_product(u_node,u_node))
          u_node=u_node/mod_u
          e_vec_2=e_vec_1
          e_vec_1=-e_vec_1
          e_vec_3=node_val(dg_position,nodes_array(node,row,3))-node_val(dg_position,node)
          e_vec_3=e_vec_3/sqrt(dot_product(e_vec_3,e_vec_3))
          cor1=dot_product(u_node,e_vec_1)
          cor2=dot_product(u_node,e_vec_2)
          cor3=dot_product(u_node,e_vec_3)
          cor4=dot_product(u_node,e_vec_3)
          if (cor1<tol) then
            local_values(row,1)=node_val(limiting_t,node)
          end if
          if (cor2<tol) then
            local_values(row,2)=node_val(limiting_t,node)
          end if
          if (cor3<tol) then
            local_values(row,3)=node_val(limiting_t,node)
          end if
          if (cor4<tol) then
            local_values(row,4)=node_val(limiting_t,node)
          end if
        end if
      end do
      ! Calculate and set the bounds
      line_max(1)=maxval(local_values(1,:))
      line_min(1)=minval(local_values(1,:))
      line_max(2)=maxval(local_values(2,:))
      line_min(2)=minval(local_values(2,:))
      node_max=minval(line_max)
      node_min=maxval(line_min)
      if (node_max<mean_val) node_max=mean_val
      if (node_min>mean_val) node_min=mean_val
      call set(lowerbound, node, node_min)
      call set(upperbound, node, node_max)
    end do

    call bound_field_diffuse(t, upperbound, lowerbound, mass, lumped_mass, inverse_lumped_mass)

    deallocate (neighbouring_nodes)
    call deallocate(inverse_lumped_mass)
    call deallocate(upperbound)
    call deallocate(lowerbound)
    if (extrapolate) then
      call deallocate(dg_position)
      deallocate (e_vec_1, dt_t, grad_t)
    end if

  end subroutine limit_fpn

end module slope_limiters_dg
