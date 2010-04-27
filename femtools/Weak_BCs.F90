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
      type(vector_field), pointer:: position, normal_nodes
      type(scalar_field), pointer:: density
      type(tensor_field), pointer:: viscosity

      integer                       :: nbcs, i, j, ele, sele, stat
      integer, dimension(:), pointer:: surface_element_list, surface_node_list
      integer, dimension(:), allocatable:: out_ele !, u_nodes_bdy
      character(len=OPTION_PATH_LEN):: bc_type, bc_path_i, bc_name

      real:: tolerance, Cb, Cf, theta

      logical:: have_Cb
      logical:: rm_out=.false.

      ! implements a penalty function for the near-wall region
      ! for near_wall_treatment: see Bazilevs et al. (2007)
      ! for log_law_of_wall    : it adds an absorption term on the wall based
      !                          on the element size and wall roughness.                     

      ewrite(1,*) "Entering wall treatment"

      velocity     => extract_vector_field(state, "Velocity")
      nl_velocity  => extract_vector_field(state, "NonlinearVelocity")
      old_velocity => extract_vector_field(state, "OldVelocity")

      position     => extract_vector_field(state, "Coordinate")
      density      => extract_scalar_field(state, "Density")
      viscosity    => extract_tensor_field(state, "Viscosity")

      call get_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/theta", &
           theta)

      nbcs=option_count(trim(velocity%option_path)//'/prognostic/boundary_conditions')
      print *, nbcs
      do i=1, nbcs
         call get_boundary_condition(velocity, i, name=bc_name, type=bc_type, &
              surface_element_list=surface_element_list)

         print *, trim(bc_name), i
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
              surface_element_list=surface_element_list, &
              surface_node_list=surface_node_list)

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

            normal_nodes => extract_surface_field(velocity, i, "normal")

         else if (bc_type=="log_law_of_wall") then

            call get_option( &
                 trim(bc_path_i)//"/type::"//trim(bc_type)//"/surface_roughness", &
                 Cf)

            normal_nodes => extract_surface_field(velocity, i, "normal")

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
                    viscosity, rhs, bigm, bc_type, tolerance, Cb, Cf, theta, have_Cb, &
                    normal_nodes, surface_node_list)

            end do
         end if
      end do

      if (rm_out) deallocate(out_ele)

      return
    end subroutine wall_functions
    
    !----------------------------------------------------------------------
    
    subroutine wall_treatment(ele, sele, x, u, nu, density, &
         viscosity, rhs, bigm, bc_type, tolerance, Cb, Cf, theta, have_Cb, &
         normal_nodes, surface_node_list)

      integer, intent(in)                       :: ele, sele
      character(len=OPTION_PATH_LEN), intent(in):: bc_type
      type(vector_field), intent(inout)         :: rhs
      type(petsc_csr_matrix), intent(inout)     :: bigm
      type(vector_field), intent(in)            :: x, u, nu, normal_nodes
      type(scalar_field), intent(in)            :: density
      type(tensor_field), intent(in)            :: viscosity
      real, dimension(face_ngi(u, sele))        :: detwei_bdy
      real, dimension(x%dim, face_ngi(u, sele)) :: normal_bdy
      real, intent(in)                          :: tolerance, Cf, theta
      logical, intent(in)                       :: have_Cb
      integer, dimension(:), intent(in)         :: surface_node_list

      ! local variables

      integer                              :: ndim, snloc, loc
      type(element_type), pointer          :: u_shape, u_f_shape, x_shape
      integer, dimension(face_loc(u, sele)):: u_nodes_bdy, sfield_nodes
      integer, dimension(:), pointer       :: ele_nodes_u, ele_faces_u

      real, dimension(x%dim, x%dim, ele_ngi(u, sele))        :: invJ
      real, dimension(x%dim, x%dim, face_ngi(u, sele))       :: invJ_face
      type(element_type)                                     :: augmented_shape
      real, dimension(ele_loc(u, ele), face_ngi(u, sele),x %dim):: vol_dshape_face
      real, dimension(x%dim, x%dim, face_ngi(u, sele))       :: viscosity_gi
      real, dimension(ele_loc(u, ele))                       :: u_dot_n_e
      real, dimension(face_loc(u, sele))                     :: u_dot_n, T
      logical, dimension(ele_loc(u, ele))                    :: msk
      real, dimension(face_ngi(u, sele))                     :: T_at_quad

      real, dimension(x%dim, face_loc(u, sele), ele_loc(u, ele)):: mat1
      real, dimension(ele_loc(u, ele), face_loc(u, sele))       :: mat2
      real, dimension(face_loc(u, sele), face_loc(u, sele))     :: mat3

      real, dimension(x%dim,x%dim)  :: G
      real, dimension(x%dim,1)      :: n
      real, dimension(x%dim)        :: normal, n_face
      real, dimension(1,1)          :: hb

      real                          :: h, tau, r, drdt, dtau, Id, Cb
      real                          :: v, uh, u_p, y_p, q
      integer                       :: i, j, dim, its, inode
      integer                       :: l_face_number, node

      ndim        = x%dim
      snloc       = face_loc(u, sele)
      u_nodes_bdy = face_global_nodes(u, sele)
      ele_nodes_u =>ele_nodes(u, ele)
      u_f_shape   =>face_shape(u, sele)

      ! find element node inside the domain
      ele_faces_u => ele_faces(u, ele)
      do i = 1, size(ele_faces_u)
         if (ele_faces_u(i)==sele) exit
      end do
      inode = ele_nodes_u(i)

      call compute_inverse_jacobian( ele_val(x, ele), &
           ele_shape(x, ele), invJ )

      call transform_facet_to_physical( x, sele, detwei_f=detwei_bdy, normal=normal_bdy )

      ! calculate wall-normal element mesh size
      G = matmul(transpose(invJ(:,:,1)), invJ(:,:,1))
      n(:,1) = normal_bdy(:,1)
      hb = 1. / sqrt( matmul(matmul(transpose(n), G), n) )
      h  = hb(1, 1)

      if (bc_type=="near_wall_treatment") then ! Hughes' approach

         ! find element nodes in surface_node_list
         do i = 1, snloc
            do j = 1, size(surface_node_list)
               if (u_nodes_bdy(i)==surface_node_list(j)) then
                  sfield_nodes(i)=j
                  exit
               end if
            end do
         end do

         ! set Cb if not specified in the flml
         if (.not.have_Cb) then 
            Cb = 1.
         end if

         do loc = 1, snloc

            ! global node number
            node = u_nodes_bdy(loc)

            ! node normal
            normal = node_val(normal_nodes, sfield_nodes(loc))

            ! velocity parallel to the wall
            uh=0.
            do i = 1, ndim
               do j = 1, ndim

                  Id = 0.
                  if (i==j) Id = 1.

                  uh = uh + &
                       (- normal(i) * normal(j) + Id) &
                       * node_val(nu, i, inode)**2
               end do
            end do
            uh = sqrt( max(uh, 0.) )

            ! assume isotropic viscosity
            v   = node_val(viscosity, 1, 1, node) /node_val(density, node)

            ! initialise tau (this is \tau_B in the paper)
            tau = Cb * (v / h)

            y_p = tau**(-0.5) * uh**(0.5)
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

               y_p  = (h/(v*Cb)) * tau**(0.5) * uh**(0.5)
               u_p  = tau**(-0.5) * uh**(0.5)

               ! update residual
               r    = y_p - log_law(u_p)

               ! check the iterative method
               if (its>50) then
                  tau = Cb * (v / h)
                  ewrite(3,*) 'WARNING: the wall treatment failed to converge', tau
                  exit
               end if
               its = its + 1
            end do

            ! generally u_p shouldn't go over 10.
            if (u_p>=25. .and. Cb>0. .and. abs(tau-Cb*(v/h))>1.e-10 .and. its>0) then
               ewrite(3,*) "WARNING: you are making the wall treatment unhappy",u_p, 'its', its
            end if

            T(loc) = max(tau/dt, 0.)

         end do ! snloc

          ! get absorption factor at quad
          T_at_quad = matmul(T, u_f_shape%n)

          ! snloc x snloc
          mat3 = shape_shape( u_f_shape, u_f_shape, &
                 detwei_bdy * T_at_quad )

          do dim = 1, ndim
            call addto( rhs, dim, u_nodes_bdy, &
                 -sum(mat3, 2) * face_val(u, dim, sele) )
            call addto_diag( bigm, dim, dim, u_nodes_bdy, &
                 dt * theta * sum(mat3, 2) )
          end do

         l_face_number = local_face_number(u, sele)
         invJ_face = spread(invJ(:, :, 1), 3, size(invJ_face, 3))

         x_shape  => ele_shape(x, ele)
         u_shape  => ele_shape(u, ele)

         augmented_shape = make_element_shape(x_shape%loc, &
              u_shape%dim, u_shape%degree, u_shape%quadrature, &
              quad_s=u_f_shape%quadrature )

         ! nloc x sngi x dim
         vol_dshape_face = eval_volume_dshape_at_face_quad( &
              augmented_shape, l_face_number, invJ_face )

         do loc = 1, face_loc(u, sele)
            n_face = node_val(normal_nodes, sfield_nodes(loc))
            u_dot_n(loc) = dot_product( &
                node_val(u, u_nodes_bdy(loc)), n_face )
         end do

         do i = 1, ele_loc(u, ele)
            do j = 1, snloc
               if (ele_nodes_u(i)==u_nodes_bdy(j)) then
                  u_dot_n_e(i)=u_dot_n(j)
               else
                  u_dot_n_e(i)=0.
               end if
            end do
         end do

         viscosity_gi = face_val_at_quad(viscosity, sele)

         ! dim x snloc x nloc
         mat1 = shape_dshape( u_f_shape, vol_dshape_face, &
              detwei_bdy * viscosity_gi(1, 1, :) * &
              face_val_at_quad(density, sele) * 2. )

         msk=.false.
         do i = 1, ele_loc(u, ele)
            if(any(u_nodes_bdy(:)==ele_nodes_u(i))) msk(i)=.true.
         end do

         do i = 1, ele_loc(u, ele)
            if (.not. msk(i)) mat1(:, :, i) = 0.
         end do

         do dim = 1, ndim
            call addto( rhs, dim, ele_nodes_u, &
                 sum(mat1(dim, :, :), 1) * u_dot_n_e )
            call addto_diag( bigm, dim, dim, ele_nodes_u, &
                 dt * theta * sum(mat1(dim, :, :), 1) )
         end do

         ! nloc x snloc
         mat2 = dshape_dot_vector_shape( &
              vol_dshape_face, normal_bdy, u_f_shape, &
              detwei_bdy * viscosity_gi(1, 1, :) * &
              face_val_at_quad(density, sele) * 2. )

         do dim = 1, ndim
            call addto( rhs, dim, ele_nodes_u, &
                 sum(mat2(:, :), 2) * ele_val(u, dim, ele) )
            call addto_diag( bigm, dim, dim, ele_nodes_u, &
                 dt * theta * sum(mat2(:, :), 2) )
         end do

         call deallocate(augmented_shape)

      else if (bc_type=="log_law_of_wall") then ! log law of the wall

         q = ( chi / (log ( (h / 2.) * Cf ) - 1.) )**2


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
                    (- normal_bdy(i, 1) * normal_bdy(j, 1) + Id) &
                    * node_val(nu, i, inode)**2
            end do
         end do
         uh = sqrt( max(uh, 0.) )

         ! snloc x snloc
         mat3 = shape_shape( u_f_shape, u_f_shape, &
                detwei_bdy * q * uh / dt)

         do dim = 1, ndim
           call addto( rhs, dim, u_nodes_bdy, &
                -sum(mat3, 2) * face_val(u, dim, sele) )
           call addto_diag( bigm, dim, dim, u_nodes_bdy, &
                dt * theta * sum(mat3, 2) )
         end do

      end if ! bc_type

    end subroutine wall_treatment

    !----------------------------------------------------------------------

    real function log_law(u)

      real, intent(in):: u

      log_law = u + e**(-chi*beta) * ( e**(chi*u) - 1. &
           - chi*u - ((chi*u)**2)/2. - ((chi*u)**3)/6.)

    end function log_law

    !----------------------------------------------------------------------

    real function jac(u, uh, v, tau, h, Cb)

      real, intent(in):: u, uh, v, tau, h, Cb

      jac = (h/(2.*v*Cb)) * tau**(-0.5) * uh**(0.5) &
           + (  1. + chi*e**(-chi*beta) * &
           ( e**(chi*u) - 1. - chi*u - ((chi*u)**2)/2.)) &
           * ((tau**(-1.5))/2.) * uh**(0.5)

    end function jac

    !----------------------------------------------------------------------

  end module Weak_BCs
