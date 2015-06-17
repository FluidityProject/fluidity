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

module linear_shallow_water
    use spud
    use fields
    use state_module
    use FLDebug
    use populate_state_module
    use write_state_module
    use populate_state_module
    use timeloop_utilities
    use sparsity_patterns_meshes
    use sparse_matrices_fields
    use solvers
    use diagnostic_variables
    use diagnostic_fields_wrapper
    use assemble_cmc
    use global_parameters, only: option_path_len, running_adjoint
    implicit none

    logical :: is_shallow_water=.false.

contains 
  subroutine setup_wave_matrices(state,u_sparsity,wave_sparsity,ct_sparsity, &
       h_mass_mat,u_mass_mat,coriolis_mat,inverse_coriolis_mat,div_mat,&
       wave_mat,big_mat,dt,theta,D0,g)
    implicit none
    type(state_type), intent(inout) :: state
    type(csr_sparsity), intent(inout) :: u_sparsity, wave_sparsity, &
         ct_sparsity
    type(csr_matrix), intent(inout) :: h_mass_mat, wave_mat
    type(block_csr_matrix), intent(inout) :: u_mass_mat, coriolis_mat,&
         inverse_coriolis_mat, div_mat, big_mat
    real , intent(in) :: dt,theta,D0,g
    !! Layer thickness
    type(scalar_field), pointer :: D, f
    !! velocity.
    type(vector_field), pointer :: U, X, down
    integer :: dim, ele

    ewrite(2,*) 'dt',dt,'theta',theta
    ewrite(2,*) 'D0',D0,'g',g

    !Pull the fields out of state
    D=>extract_scalar_field(state, "LayerThickness")
    f=>extract_scalar_field(state, "Coriolis")
    U=>extract_vector_field(state, "LocalVelocity")
    X=>extract_vector_field(state, "Coordinate")
    down=>extract_vector_field(state, "GravityDirection")

    dim = mesh_dim(U)

    !construct/extract sparsities
    u_sparsity=get_csr_sparsity_firstorder(state, U%mesh, U%mesh)
    wave_sparsity=get_csr_sparsity_firstorder(state, D%mesh, D%mesh)
    ct_sparsity=get_csr_sparsity_firstorder(state, D%mesh, U%mesh)

    !allocate matrices
    call allocate(h_mass_mat,wave_sparsity)
    call zero(h_mass_mat)
    call allocate(u_mass_mat,u_sparsity, (/dim,dim/))
    call zero(u_mass_mat)
    call allocate(coriolis_mat,u_sparsity, (/dim,dim/))
    call zero(coriolis_mat)
    call allocate(inverse_coriolis_mat,u_sparsity, (/dim,dim/))
    call zero(inverse_coriolis_mat)
    call allocate(div_mat,ct_sparsity,(/1,dim/))
    call zero(div_mat)
    call allocate(wave_mat,wave_sparsity)
    call zero(wave_mat)
    call allocate(big_mat,u_sparsity,(/dim,dim/))
    call zero(big_mat)

    !Assemble matrices
    do ele = 1, ele_count(D)
       call assemble_shallow_water_matrices_ele(D,f,U,X,down,ele, &
            h_mass_mat,u_mass_mat,coriolis_mat,inverse_coriolis_mat,&
            div_mat,big_mat,dt,theta)
    end do

    if(have_option("/physical_parameters/coriolis") .and. have_option("/debug/check_inverse_coriolis_matrix")) then
       call check_big_mat(U,big_mat,u_mass_mat,coriolis_mat,theta,dt)
    end if

    !Construct wave mat
    call assemble_cmc_dg(wave_mat, div_mat, div_mat, big_mat)

    if(have_option("/debug/check_wave_matrix")) then
      call check_wave_mat(wave_mat,div_mat,big_mat,u,d)
    end if

    call scale(wave_mat,dt*dt*theta*theta*g*D0)
    call addto(wave_mat,h_mass_mat)

  end subroutine setup_wave_matrices

    subroutine check_wave_mat(wave_mat,div_mat,big_mat,u,d)
      type(block_csr_matrix), intent(in) :: div_mat,big_mat
      type(csr_matrix), intent(in) :: wave_mat
      type(scalar_field), intent(inout) :: d
      type(vector_field), intent(inout) :: u
      !
      type(scalar_field) :: d_mem1, d_mem2
      type(vector_field) :: u_mem1, u_mem2
      integer :: n

      call allocate(d_mem1,d%mesh,'mem1')
      call allocate(d_mem2,d%mesh,'mem2')
      call allocate(u_mem1,mesh_dim(u),u%mesh,'mem1')
      call allocate(u_mem2,mesh_dim(u),u%mesh,'mem2')

      do n = 1, 20
         call random_number(d_mem1%val)
         
         call mult(d_mem2,wave_mat,d_mem1)
         call mult_t(u_mem1,div_mat,d_mem1)
         call mult(u_mem2,big_mat,u_mem1)
         call mult(d_mem1,div_mat,u_mem2)
         
         call addto(d_mem1,d_mem2,scale=-1.0)
         ewrite(2,*) 'SW: Residual = ',maxval(abs(d_mem1%val))
         ewrite(2,*) 'SW: Value = ',maxval(abs(d_mem2%val))
         if(maxval(abs(d_mem1%val))>1.0e-7) then
            FLAbort('check wave mat failed')
         end if
      end do

      call deallocate(d_mem1)
      call deallocate(d_mem2)
      call deallocate(u_mem1)
      call deallocate(u_mem2)

    end subroutine check_wave_mat

    subroutine check_big_mat(U,big_mat,u_mass_mat,coriolis_mat,theta,dt)
      implicit none
      type(vector_field), intent(inout) :: U
      type(block_csr_matrix), intent(in) :: big_mat, coriolis_mat, u_mass_mat
      real, intent(in) :: theta, dt
      !
      integer :: n,ntests=20, dim,d1
      type(vector_field) :: u_test, u_mem1, u_mem2

      dim = mesh_dim(U)

      call allocate(u_test,dim,U%mesh,'testmem')
      call allocate(u_mem1,dim,U%mesh,'testmem1')
      call allocate(u_mem2,dim,U%mesh,'testmem2')

      do n = 1, ntests
         do d1 = 1, dim
            call random_number(u_test%val(d1,:))
         end do

         !First apply mass matrix
         call mult(u_mem1,u_mass_mat,u_test)
         !Then coriolis
         call mult(u_mem2,coriolis_mat,u_test)
         call addto(u_mem1,u_mem2,scale=dt*theta)
         !The application of the inverse of big_mat is now in u_mem1
         !Now apply big_mat to u_mem1 and put it in u_mem2
         call mult(u_mem2,big_mat,u_mem1)
         !Now subtract u_test from it
         call addto(u_mem2,u_test,scale=-1.0)
         !Now check the residual
         do d1 = 1, dim
            if(maxval(abs(u_mem2%val(d1,:)))>1.0e-7) then
               ewrite(-1,*) 'residual of ',maxval(abs(u_mem2%val(d1,:)))
               ewrite(-1,*) 'in component',d1
               FLAbort('check big mat failed')
            end if
         end do
      end do

      ewrite(2,*) 'SW: checking big mat complete.'

      call deallocate(u_test)
      call deallocate(u_mem1)
      call deallocate(u_mem2)
    end subroutine check_big_mat

    subroutine check_solution(delta_u,delta_d,d,u,dt,theta,g,D0,u_mass_mat&
           &,h_mass_mat, coriolis_mat,div_mat)
      type(scalar_field), intent(inout) :: delta_d,d
      type(vector_field), intent(inout) :: delta_u,u
      type(csr_matrix), intent(in) :: h_mass_mat
      type(block_csr_matrix), intent(in) :: u_mass_mat,coriolis_mat,div_mat
      real, intent(in) :: dt,theta,g,D0
      !
      type(scalar_field) :: h_residual, h_mem
      type(vector_field) :: u_residual, u_mem
      integer :: dim
      !
      dim = mesh_dim(u)
      !
      call allocate(h_residual,d%mesh,'hresidual')
      call allocate(u_residual,mesh_dim(u),u%mesh,'uresidual')
      call allocate(h_mem,d%mesh,'hmem')
      call allocate(u_mem,mesh_dim(u),u%mesh,'umem')
      !
      call zero(h_residual)
      call zero(h_mem)
      call zero(u_residual)
      call zero(u_mem)
      !
      ! U equation
      ! time-derivative
      call mult(u_residual,u_mass_mat,delta_u)
      if(dim==2) then
         ! Coriolis implicit
         call mult(u_mem,coriolis_mat,delta_u)
         call addto(u_residual,u_mem,scale=dt*theta)
         ! Coriolis explicit
         call mult(u_mem,coriolis_mat,u)
         call addto(u_residual,u_mem,scale=dt)
      end if
      !pressure implicit
      call mult_t(u_mem,div_mat,delta_d)
      call addto(u_residual,u_mem,scale=g*dt*theta)
      !pressure explicit
      call mult_t(u_mem,div_mat,d)
      call addto(u_residual,u_mem,scale=g*dt)      
      !
      ! D equation
      ! time-derivative
      call mult(h_residual,h_mass_mat,delta_d)
      call set(u_mem,delta_u)
      call scale(u_mem,theta)
      call addto(u_mem,u)
      call mult(h_mem,div_mat,u_mem)
      call addto(h_residual,h_mem,scale=-dt*D0)
      !
      ewrite(2,*) dt,theta,g,D0
      ewrite(2,*) 'SW: h residual', maxval(abs(h_residual%val))
      ewrite(2,*) 'SW: u1 residual', maxval(abs(u_residual%val(1,:)))
      ewrite(2,*) 'SW: u size', maxval(abs(u%val(1,:)))
      if(dim==2) then
         ewrite(2,*) 'SW: u2 residual', maxval(abs(u_residual%val(2,:)))
      ewrite(2,*) 'SW: u size', maxval(abs(u%val(2,:)))
      end if
      call deallocate(h_residual)
      call deallocate(u_residual)
      call deallocate(u_mem)
      call deallocate(h_mem)
      !
    end subroutine check_solution

    subroutine update_u_rhs(u_rhs,U,delta_D,div_mat,theta,dt,g)
      implicit none
      type(vector_field), intent(inout) :: u_rhs
      type(vector_field), intent(in) :: U
      type(scalar_field), intent(in) :: delta_D
      type(block_csr_matrix), intent(in) :: div_mat
      real, intent(in) :: theta, dt, g
      !Add the contribution to the rhs of momentum equation which
      !contains implicit D
      type(vector_field) :: vec
      integer :: dim

      dim = mesh_dim(U)
      call allocate(vec, dim, U%mesh, "workingvecmem")

      call mult_T(vec, div_mat, delta_D)
      call addto(u_rhs, vec, scale = -dt*g*theta)

      call deallocate(vec)

    end subroutine update_u_rhs

    subroutine get_d_rhs(d_rhs,u_rhs,D,U,div_mat,big_mat,D0,dt,theta)
      implicit none
      type(scalar_field), intent(inout) :: d_rhs
      type(vector_field), intent(inout) :: u_rhs,U
      type(scalar_field), intent(inout) :: D
      type(block_csr_matrix), intent(in) :: div_mat, big_mat
      real, intent(in) :: D0, dt, theta
      !Construct the contribution to the right-hand side of the D wave eqn
      type(scalar_field) :: rhs2
      integer :: dim
      type(vector_field) :: vec
      !
      dim = mesh_dim(D)

      call zero(d_rhs)

      call allocate(rhs2, D%mesh, "workingmem")
      call allocate(vec, dim, U%mesh, "workingvecmem")

      !Contribution to u^{n+1} from explicit bits
      call mult(vec,big_mat,u_rhs)
      call scale(vec,theta)
      !Contribution from previous timestep
      call addto(vec,U)
      !Continuity operator
      call mult(rhs2, div_mat, vec)
      call addto(d_rhs, rhs2, scale=dt*D0)

      call deallocate(rhs2)
      call deallocate(vec)

    end subroutine get_d_rhs

    subroutine get_u_rhs(u_rhs,U,D,dt,g, &
         coriolis_mat,div_mat,u_mass_mat,source)
      implicit none
      type(vector_field), intent(inout) :: u_rhs
      type(vector_field), intent(inout) :: U
      type(vector_field), intent(in), optional :: source
      type(scalar_field), intent(inout) :: D
      type(block_csr_matrix), intent(in) :: coriolis_mat, div_mat 
      type(block_csr_matrix), intent(in), optional :: u_mass_mat
      real, intent(in) :: dt, g
      !Construct the explicit contribution to the right-hand side 
      !of the u equation

      integer :: dim, i
      type(vector_field) :: rhs2, vec

      ewrite(2,*) 'outside', sum(D%val)

      dim = mesh_dim(U)

      call allocate(rhs2, dim, U%mesh, "workingmem")
      call allocate(vec, dim, U%mesh, "workingvecmem")

      call zero(u_rhs)

      !Coriolis term
      call mult(rhs2,coriolis_mat,U)
      call addto(U_rhs,rhs2,scale=-dt)

      !Pressure gradient
      call mult_T(vec,div_mat,D)
      call addto(u_rhs,vec,scale=-dt*g)

      !Source term
      if (present(source)) then
        call mult(vec,u_mass_mat,source)
        call addto(u_rhs,vec,scale=dt)
      end if

      do i=1,u_rhs%dim
         ewrite(1,*) 'SW u_rhs',maxval(abs(u_rhs%val(i,:)))
      end do

      call deallocate(rhs2)
      call deallocate(vec)
    end subroutine get_u_rhs

    subroutine assemble_shallow_water_matrices_ele(D,f,U,X,down,ele, &
         h_mass_mat,u_mass_mat,coriolis_mat,inverse_coriolis_mat,&
         div_mat,big_mat,dt,theta)

      implicit none
      type(scalar_field), intent(in) :: D, f
      type(vector_field), intent(in) :: U, X, down
      type(csr_matrix), intent(inout) :: h_mass_mat
      type(block_csr_matrix), intent(inout) :: u_mass_mat, coriolis_mat,&
           inverse_coriolis_mat, div_mat, big_mat

      integer, intent(in) :: ele
      real, intent(in) :: dt,theta

      !Assemble h_mass_mat, u_mass_mat, coriolis_mat, inverse_coriolis_mat,
      !div_mat and then big_mat
      real, dimension(ele_ngi(D, ele)) :: detwei
      integer, dimension(:), pointer :: D_ele, U_ele
      type(element_type) :: D_shape, u_shape
      real, dimension(ele_ngi(x,ele)) :: f_gi
      real, dimension(X%dim, ele_ngi(X,ele)) :: up_gi
      real, dimension(mesh_dim(U)*ele_loc(U,ele), &
                 mesh_dim(U)*ele_loc(U,ele)) :: l_big_mat, l_coriolis_mat
      integer :: dim1, dim2, nloc, dim, gi
      real, dimension(mesh_dim(U), X%dim, ele_ngi(U,ele)) :: J
      real, dimension(ele_ngi(U,ele)) :: detJ
      real, dimension(mesh_dim(U), mesh_dim(U), ele_ngi(U,ele)) :: G, Gf
      real, dimension(X%dim, X%dim, ele_ngi(U,ele)) :: rot
      real, dimension(mesh_dim(U),mesh_dim(U),ele_loc(U,ele),ele_loc(U,ele)) :: l_u_mat, l_uf_mat

      real, dimension(mesh_dim(U),ele_loc(D,ele),ele_loc(U,ele)) :: l_div_mat
      real, dimension(X%dim) :: up_vec

      dim = mesh_dim(U)

      u_shape=ele_shape(u, ele)
      D_shape=ele_shape(d, ele)

      D_ele => ele_nodes(D, ele)
      U_ele => ele_nodes(U, ele)

      f_gi = ele_val_at_quad(f,ele)
      up_gi = -ele_val_at_quad(down,ele)

      call compute_jacobian(X, ele, J=J, detwei=detwei, detJ=detJ)

      call addto(h_mass_mat, D_ele, D_ele, &
           shape_shape(D_shape,D_shape,detwei))

      l_div_mat = dshape_shape(D_shape%dn,u_shape,D_shape%quadrature%weight)

      do dim1 = 1, dim
         call addto(div_mat,1,dim1,d_ele,u_ele,l_div_mat(dim1,:,:))
      end do

      !up vector must be normal to the surface for steady geostrophic
      !states.
      do gi=1, ele_ngi(U,ele)
         up_vec = get_up_vec(ele_val(X,ele), up_gi(:,gi))
         up_gi(:,gi) = up_vec
      end do

      do gi=1, ele_ngi(U,ele)
         rot(1,:,gi)=(/0.,-up_gi(3,gi),up_gi(2,gi)/)
         rot(2,:,gi)=(/up_gi(3,gi),0.,-up_gi(1,gi)/)
         rot(3,:,gi)=(/-up_gi(2,gi),up_gi(1,gi),0./)
      end do

      do gi=1,ele_ngi(U,ele)
         G(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)
         Gf(:,:,gi)=matmul(J(:,:,gi), &
              matmul(rot(:,:,gi), transpose(J(:,:,gi))))/detJ(gi)
      end do

      l_u_mat = shape_shape_tensor(u_shape, u_shape, &
           u_shape%quadrature%weight, G)
      l_uf_mat = shape_shape_tensor(u_shape, u_shape, &
           f_gi*u_shape%quadrature%weight, Gf)

      do dim1 = 1, dim
         do dim2 = 1, dim
            call addto(u_mass_mat, dim1, dim2, u_ele, u_ele, &
                 l_u_mat(dim1,dim2,:,:))
            call addto(coriolis_mat, dim1, dim2, u_ele, u_ele, &
                 l_uf_mat(dim1,dim2,:,:))
         end do
      end do

      l_big_mat = 0.

      nloc = ele_loc(U,ele)
      do dim1 = 1, dim
         do dim2 = 1, dim
            l_big_mat(nloc*(dim1-1)+1:nloc*dim1, &
                 nloc*(dim2-1)+1:nloc*dim2) = &
                 l_u_mat(dim1,dim2,:,:)+&
                 dt*theta*l_uf_mat(dim1,dim2,:,:)
         end do
      end do

      call invert(l_big_mat)

      do dim1 = 1, dim
         do dim2 = 1, dim
            call addto(big_mat,dim1,dim2,u_ele,u_ele, &
                 l_big_mat(nloc*(dim1-1)+1:nloc*dim1, &
                 nloc*(dim2-1)+1:nloc*dim2))
         end do
      end do

      ! Compute inverse_coriolis_mat for use in 
      ! set_velocity_from_geostrophic_balance
      if (have_option("/physical_parameters/coriolis")) then
        l_coriolis_mat = 0.
        do dim1 = 1, dim
           do dim2 = 1, dim
              l_coriolis_mat(nloc*(dim1-1)+1:nloc*dim1, &
                   nloc*(dim2-1)+1:nloc*dim2) = &
                   l_uf_mat(dim1,dim2,:,:)
           end do
        end do

        call invert(l_coriolis_mat)

        do dim1 = 1, dim
           do dim2 = 1, dim
              call addto(inverse_coriolis_mat,dim1,dim2,u_ele,u_ele, &
                   l_coriolis_mat(nloc*(dim1-1)+1:nloc*dim1, &
                   nloc*(dim2-1)+1:nloc*dim2))
           end do
        end do
      end if

      contains
        !function for getting normal to surface
        !only works for flat elements
        function get_up_vec(X_val, up) result (up_vec_out)
          real, dimension(:,:), intent(in) :: X_val           !(dim,loc)
          real, dimension(:), intent(in) :: up
          real, dimension(size(X_val,1)) :: up_vec_out
          !
          real, dimension(size(X_val,1)) :: t1,t2
          ! if elements are triangles:
          if(size(X_val,2)==3) then
             t1 = X_val(:,2)-X_val(:,1)
             t2 = X_val(:,3)-X_val(:,1)
             up_vec_out(1) = t1(2)*t2(3)-t1(3)*t2(2)
             up_vec_out(2) = -(t1(1)*t2(3)-t1(3)*t2(1))
             up_vec_out(3) = t1(1)*t2(2)-t1(2)*t2(1)
             up_vec_out = up_vec_out*dot_product(up_vec_out, up)
             up_vec_out = up_vec_out/sqrt(sum(up_vec_out**2))
          else
             up_vec_out = up
          end if
        end function get_up_vec

    end subroutine assemble_shallow_water_matrices_ele

    subroutine get_linear_energy(state,U_mass_mat,h_mass_mat,d0,g,energy_out)
      implicit none
      real, optional, intent(out) :: energy_out
      real, intent(in) :: d0,g
      type(state_type), intent(inout) :: state
      type(block_csr_matrix), intent(in) :: u_mass_mat
      type(csr_matrix), intent(in) :: h_mass_mat
      !
      type(scalar_field), pointer :: d
      type(vector_field), pointer :: u
      type(scalar_field) :: Md
      type(vector_field) :: Mu
      real :: D_l2,u_l2, energy
      integer :: d1, dim

      D=>extract_scalar_field(state, "LayerThickness")
      U=>extract_vector_field(state, "LocalVelocity")

      dim = mesh_dim(U)
      call allocate(Md,D%mesh,'Md')
      call allocate(Mu,mesh_dim(U),u%mesh,'Mu')
      !
      call mult(Md,h_mass_mat,D)
      call mult(Mu,u_mass_mat,U)
      D_l2 = sum(Md%val*d%val)
      U_l2 = 0.
      do d1 = 1, dim
         U_l2 = U_l2 + sum(Mu%val(d1,:)*u%val(d1,:))
      end do
      energy = 0.5*g*D_l2 + 0.5*D0*U_l2
      ewrite(2,*) 'SW: energy = ', energy
      !
      energy_out = energy

      call set_diagnostic(name="LinearEnergy", statistic="value", value=(/energy/))

      call deallocate(Md)
      call deallocate(Mu)
    end subroutine get_linear_energy

    subroutine linear_shallow_water_register_diagnostic
      
      if(is_shallow_water .and. .not. running_adjoint) then
         call register_diagnostic(dim=1, name="LinearEnergy", statistic="value")
      end if

    end subroutine linear_shallow_water_register_diagnostic

end module linear_shallow_water
