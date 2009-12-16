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

  module Weak_BCs

    use FLDebug
    use generic_interface
    use spud
    use state_module
    use global_parameters, only: option_path_len, dt
    use fields
    use elements
    use parallel_tools
    use transform_elements
    use boundary_conditions
    use elements
    use sparse_tools_petsc

    implicit none

    private
    public :: wall_functions

    real, parameter:: chi = 0.4 ! VonKarman constant (~0.4)
    real, parameter:: beta= 5.5 ! log law shift (~5. - 5.5)
    real, parameter:: e   = 2.718281828

  contains

    subroutine wall_functions(bigm, rhs, state)

      type(petsc_csr_matrix), intent(inout):: bigm
      type(vector_field), intent(inout)    :: rhs
      type(state_type), intent(in)         :: state

      type(vector_field), pointer:: velocity, old_velocity, nl_velocity
      type(vector_field), pointer:: position, absorption
      type(scalar_field), pointer:: density
      type(tensor_field), pointer:: viscosity

      integer                       :: nbcs, i, j, ele, sele, stat
      integer, dimension(:), pointer:: surface_element_list
      integer, dimension(:), allocatable:: out_ele
      character(len=OPTION_PATH_LEN):: bc_type, bc_path_i

      real:: tolerance, Cb, Cf, theta

      logical:: have_absorption, have_Cb
      logical:: rm_out=.false.

      ewrite(1,*) "Entering wall treatment"

      velocity     => extract_vector_field(state, "Velocity")
      nl_velocity  => extract_vector_field(state, "NonlinearVelocity")
      old_velocity => extract_vector_field(state, "OldVelocity")

      position     => extract_vector_field(state, "Coordinate")
      density      => extract_scalar_field(state, "Density")
      viscosity    => extract_tensor_field(state, "Viscosity")

      absorption=> extract_vector_field(state, "VelocityAbsorption", stat)

      have_absorption = stat == 0
      if(.not. have_absorption) then
         FLExit("you neeed to turn on the VelocityAbsorption...")
      end if

      call get_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/theta", &
           theta)

      nbcs=option_count(trim(velocity%option_path)//'/prognostic/boundary_conditions')
      do i=1, nbcs
         call get_boundary_condition(velocity, i, type=bc_type, &
              surface_element_list=surface_element_list)

         bc_path_i = &
              trim(velocity%option_path)//"/prognostic/boundary_conditions"//"["//int2str(i-1)//"]"

         if (bc_type=="outflow") then

            allocate(out_ele(size(surface_element_list)))

            do j=1, size(surface_element_list)
               sele=surface_element_list(j)
               out_ele(j) = face_ele(position, sele)
            end do

            rm_out=.true.

         end if
      end do
      do i = 1, nbcs
         call get_boundary_condition(velocity, i, type=bc_type, &
              surface_element_list=surface_element_list)

         bc_path_i = &
              trim(velocity%option_path)//"/prognostic/boundary_conditions"//"["//int2str(i-1)//"]"

         tolerance=0.; Cb=0.; have_Cb=.false.; Cf=0.
         if (bc_type=="near_wall_treatment") then

            call get_option( &
                 trim(bc_path_i)//"/type::"//trim(bc_type)//"/tolerance", &
                 tolerance)

            call get_option( &
                 trim(bc_path_i)//"/type::"//trim(bc_type)//"/Cb", &
                 Cb, stat)

            have_Cb = stat == 0

         else if (bc_type=="log_law_of_wall") then

            call get_option( &
                 trim(bc_path_i)//"/type::"//trim(bc_type)//"/surface_roughness", &
                 Cf)

         end if

         if (bc_type=="near_wall_treatment" .or. bc_type=="log_law_of_wall" ) then

            ewrite(2,*) "wall treatment method: ", trim(bc_type)
            if (bc_type=="near_wall_treatment") then
               ewrite(2,*) "       tolerance: ", tolerance, "Cb: ", Cb, "have_Cb: ", have_Cb, "out: ", rm_out
            else if (bc_type=="log_law_of_wall") then
               ewrite(2,*) "       Cf: ", Cf, "out: ", rm_out
            end if

            do j = 1, size(surface_element_list)

               sele = surface_element_list(j)
               ele  = face_ele(position, sele)

               if (rm_out) then
                  if (any(out_ele(:)==ele)) cycle
               end if

               call wall_treatment(ele, sele, position, old_velocity, nl_velocity, density, &
                    absorption, viscosity, rhs, bigm, bc_type, tolerance, Cb, Cf, theta, have_Cb)

            end do
         end if
      end do

      if (rm_out) deallocate(out_ele)

      return
    end subroutine wall_functions
    
    !----------------------------------------------------------------------
    
    subroutine wall_treatment(ele, sele, x, u, nu, density, absorption, &
         viscosity, rhs, bigm, bc_type, tolerance, Cb, Cf, theta, have_Cb)

      integer, intent(in)                       :: ele, sele
      character(len=OPTION_PATH_LEN), intent(in):: bc_type
      type(vector_field), intent(inout)         :: rhs, absorption
      type(petsc_csr_matrix), intent(inout)     :: bigm
      type(vector_field), intent(in)            :: x, u, nu
      type(scalar_field), intent(in)            :: density
      type(tensor_field), intent(in)            :: viscosity
      real, dimension(face_ngi(u, sele))        :: detwei_bdy
      real, dimension(x%dim, face_ngi(u, sele)) :: normal_bdy
      real, intent(in)                          :: tolerance, Cf, theta
      logical, intent(in)                       :: have_Cb

      ! local variables

      integer                              :: ndim, snloc, loc
      type(element_type), pointer          :: u_shape, u_f_shape, x_shape
      integer, dimension(face_loc(u, sele)):: u_nodes_bdy

      real, dimension(x%dim, x%dim, ele_ngi(u, sele))        :: invJ
      real, dimension(x%dim, x%dim, face_ngi(u, sele))       :: invJ_face
      type(element_type)                                     :: augmented_shape
      real,dimension(ele_loc(u, ele),face_ngi(u, sele),x%dim):: vol_dshape_face

      !real, dimension(ele_loc(u, ele))                :: u_ele_val
      !real, dimension(x%dim, x%dim, face_ngi(u, sele)):: grad_u_at_quad
      !real, dimension(x%dim, face_ngi(u, sele))       :: grad_u_dot_n_at_quad
      integer, dimension(:), pointer                  :: ele_nodes_u

      !real, dimension(x%dim, face_loc(u, sele))            :: a1
      real, dimension(ele_loc(u, ele),face_loc(u, sele))   :: a2

      real, dimension(face_loc(u, sele), ele_loc(u, ele))  :: a11

      real, dimension(x%dim,x%dim)  :: G
      real, dimension(x%dim,1)      :: n
      real, dimension(1,1)          :: hb

      real                          :: h, tau, r, drdt, dtau, Id, Cb
      real                          :: v, uh, u_p, y_p, q, tmp_val
      integer                       :: i, j, dim, its
      integer                       :: l_face_number, node

      !real, dimension(ele_loc(u, ele))  :: dshape_dot_normal
      !real, dimension(face_loc(u, sele)):: shp
      !real                              :: dshape_dot_normal_l2, shp_l2

      ! implements a penalty function for the near-wall region

      ndim        = x%dim
      snloc       = face_loc(u, sele)
      u_nodes_bdy = face_global_nodes(u, sele)
      ele_nodes_u =>ele_nodes(u, ele)

      call compute_inverse_jacobian( ele_val(x, ele), &
           ele_shape(x, ele), invJ )

      call transform_facet_to_physical( x, sele, detwei_f=detwei_bdy, normal=normal_bdy )

      ! calculate wall-normal element mesh size
      G = matmul(transpose(invJ(:,:,1)), invJ(:,:,1))
      n(:,1) = normal_bdy(:,1)
      hb = 1. / sqrt( matmul(matmul(transpose(n), G), n) )
      h  = hb(1,1)

      if (bc_type=="near_wall_treatment") then ! Hughes' approach

         ! double the element size
         h = 2. * h

         ! calculate Cb
         l_face_number = local_face_number(u, sele)
         invJ_face = spread(invJ(:, :, 1), 3, size(invJ_face, 3))

         x_shape => ele_shape(x, ele)
         u_shape => ele_shape(u, ele)
         u_f_shape => face_shape(u, sele)

         augmented_shape = make_element_shape(x_shape%loc, &
              u_shape%dim, u_shape%degree, u_shape%quadrature, &
              quad_s=u_f_shape%quadrature )

         ! nloc x sngi x dim
         vol_dshape_face = eval_volume_dshape_at_face_quad( &
              augmented_shape, l_face_number, invJ_face )

         if (.not.have_Cb) then 
            !dshape_dot_normal = dshape_dot_vector_rhs(vol_dshape_face, normal_bdy, detwei_bdy)
            !shp = shape_rhs(u_f_shape, detwei_bdy)
            !dshape_dot_normal_l2=0.; shp_l2 =0.
            !do i=1,snloc
            !  dshape_dot_normal_l2 = dshape_dot_normal_l2 + dshape_dot_normal(i)**2
            !  shp_l2 = shp_l2 + shp(i)**2
            !end do
            !dshape_dot_normal_l2 = sqrt(dshape_dot_normal_l2)
            !shp_l2 = sqrt(shp_l2)
            Cb = 8. * h !* dshape_dot_normal_l2 / shp_l2
         end if

         do loc = 1, snloc

            ! global node number
            node = u_nodes_bdy(loc)

            ! velocity parallel to the wall
            uh=0.
            do i = 1, ndim
               do j = 1, ndim
                  if (i==j) then
                     Id = 1. 
                  else
                     Id = 0.
                  end if
                  uh = uh + &
                       (- normal_bdy(i,1) * normal_bdy(j,1) + Id) &
                       * node_val(nu, i, node)**2
               end do
            end do
            uh = sqrt( max(uh, 0.) )

            ! assume isotropic viscosity
            v   = node_val(viscosity, 1, 1, node) /node_val(density, node)

            ! initialise tau
            tau = Cb * (v / h)

            y_p = tau**(-1./2.) * uh**(1./2.)
            u_p = y_p

            ! initialise residual
            r   = y_p - log_law(u_p)

            its = 0
            do while (abs(r) > tolerance)
               drdt = jac(u_p, uh, v, tau, h, Cb)
               dtau = - r / drdt

               ! update tau
               tau = tau + dtau

               ! check that tau is positive
               if (tau < 0.) then
                  tau = tau - dtau
                  ewrite(3,*) "WARNING: inside near wall treatment"
                  ewrite(3,*) "got negative tau, reverting",tau,tau + dtau
                  exit
               end if

               y_p  = (h/(v*Cb)) * tau**(1./2.) * uh**(1./2.)
               u_p  = tau**(-1./2.) * uh**(1./2.)

               ! update residual
               r    = y_p - log_law(u_p)

               ! checking convergence of the iterative method
               if (r>=25. .and. its>50) FLExit("near wall treatment: diverged")
               if (its>=100)            FLExit("near wall treatment: its")

               its = its + 1
            end do

            ! generally u_p shouldn't go over 10.
            if (u_p>=25.) then
               ewrite(3,*) "WARNING: you are making the wall treatment unhappy",u_p
            end if

            ! add absorption term...
            do dim = 1, ndim
               tmp_val = node_val(absorption, dim, node)
               absorption%val(dim)%ptr(node) = tmp_val + max(tau, 0.)
            end do

         end do ! snloc

!         ! Calculate grad u at the surface element quadrature points
!         do dim = 1, ndim
!            u_ele_val = ele_val(u, dim, ele)
!            forall(i = 1:ndim, j = 1:face_ngi(u, sele))
!               grad_u_at_quad(dim, i, j) = dot_product( &
!                    u_ele_val, vol_dshape_face(:, j, i) )
!            end forall
!         end do
!
!         ! Calculate grad u dot n at the surface element quadrature points
!         do dim = 1, ndim
!            do i = 1, face_ngi(u, sele)
!               grad_u_dot_n_at_quad(dim, i) = dot_product( &
!                    grad_u_at_quad(:, dim, i), normal_bdy(:, i) )
!            end do
!         end do
!
!         !  dim x snloc
!         a1 = shape_vector_rhs( &
!              u_f_shape,grad_u_dot_n_at_quad, &
!              detwei_bdy*v*face_val_at_quad(density, sele)*2 )
!
!         do dim = 1, ndim
!            call addto( rhs, dim, u_nodes_bdy, -a1(dim, :) )
!
!            call addto_diag( bigm, dim, dim, u_nodes_bdy, &
!                 dt*theta*a1(dim,:) )
!         end do

         ! snloc x nloc
         a11 = shape_vector_dot_dshape( &
               u_f_shape, normal_bdy, vol_dshape_face, &
               detwei_bdy*v*face_val_at_quad(density, sele)*2 )

         do dim = 1, ndim
            call addto( rhs, dim, ele_nodes_u, &
                 -sum(a11, 2)*face_val(u, dim, sele) )

            call addto_diag( bigm, dim, dim, u_nodes_bdy, &
                 dt*theta*sum(a11, 2) )
         end do

         ! nloc x snloc
         a2 = dshape_dot_vector_shape( &
              vol_dshape_face, normal_bdy, u_f_shape, &
              detwei_bdy*v*face_val_at_quad(density, sele)*2 )

         do dim = 1, ndim  
            call addto( rhs, dim, ele_nodes_u, &
                 -sum(a2, 2)*ele_val(u, dim, ele) )

            call addto_diag( bigm, dim, dim, ele_nodes_u, &
                 dt*theta*sum(a2, 2) )
         end do

         call deallocate(augmented_shape)

      else if (bc_type=="log_law_of_wall") then ! log law of the wall

         q = ( chi / (log ( (h / 2) * Cf ) - 1) )**2

         do loc = 1, snloc

            ! global node number
            node = u_nodes_bdy(loc)

            ! velocity parallel to the wall
            uh=0.
            do i = 1, ndim
               do j = 1, ndim
                  if (i==j) then
                     Id = 1. 
                  else
                     Id = 0.
                  end if
                  uh = uh + &
                       (- normal_bdy(i,1) * normal_bdy(j,1) + Id) &
                       * node_val(nu, i, node)**2
               end do
            end do
            uh = sqrt( max(uh, 0.) )

         end do ! snloc

         ! add absorption term...
         do dim = 1, ndim
            tmp_val = node_val(absorption, dim, node)
            absorption%val(dim)%ptr(node) = tmp_val + max(q*uh, 0.)
         end do

      end if

    end subroutine wall_treatment

    !----------------------------------------------------------------------

    real function log_law(u)

      real, intent(in):: u

      log_law = u + e**(-chi*beta) * ( e**(chi*u) - 1   &
           - chi*u - ((chi*u)**2)/2 - ((chi*u)**3)/6)

    end function log_law

    !----------------------------------------------------------------------

    real function jac(u, uh, v, tau, h, Cb)

      real, intent(in):: u, uh, v, tau, h, Cb

      jac = (h/(2*v*Cb)) * tau**(-1./2.) * uh**(1./2.)  &
           + (  1 + chi*e**(-chi*beta) *                &
           ( e**(chi*u) -1 - chi*u - ((chi*u)**2)/2))   &
           * ((tau**(-3./2.))/2) * uh**(1./2.)

    end function jac

    !----------------------------------------------------------------------

  end module Weak_BCs
