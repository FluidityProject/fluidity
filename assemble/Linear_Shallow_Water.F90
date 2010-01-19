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
    use global_parameters, only: option_path_len

contains 
  subroutine setup_wave_matrices(state,u_sparsity,wave_sparsity,ct_sparsity, &
       h_mass_mat,u_mass_mat,coriolis_mat,div_mat,wave_mat,big_mat, &
       dt,theta,D0,g,f0,beta)
    implicit none
    type(state_type), intent(inout) :: state
    type(csr_sparsity), intent(inout) :: u_sparsity, wave_sparsity, &
         ct_sparsity
    type(csr_matrix), intent(inout) :: h_mass_mat, coriolis_mat, &
         u_mass_mat,wave_mat
    type(block_csr_matrix), intent(inout) :: div_mat, &
         big_mat
    real , intent(in) :: dt,theta,D0,g,f0
    real, dimension(:), intent(in) :: beta
    !! Layer thickness
    type(scalar_field), pointer :: D
    !! velocity.
    type(vector_field), pointer :: U,X
    integer :: dim, ele
    ! inverse of big_mat for testing
    type(block_csr_matrix) :: inverse_big_mat

    ewrite(2,*) 'dt',dt,'theta',theta
    ewrite(2,*) 'D0',D0,'g',g

    !Pull the fields out of state
    D=>extract_scalar_field(state, "LayerThickness")
    U=>extract_vector_field(state, "Velocity")
    X=>extract_vector_field(state, "Coordinate")

    dim = mesh_dim(D)

    !construct/extract sparsities
    u_sparsity=get_csr_sparsity_firstorder(state, U%mesh, U%mesh)
    wave_sparsity=get_csr_sparsity_firstorder(state, D%mesh, D%mesh)
    ct_sparsity=get_csr_sparsity_firstorder(state, D%mesh, U%mesh)

    !allocate matrices
    call allocate(h_mass_mat,wave_sparsity)
    call zero(h_mass_mat)
    call allocate(u_mass_mat,u_sparsity)
    call zero(u_mass_mat)
    call allocate(coriolis_mat,u_sparsity)
    call zero(coriolis_mat)
    call allocate(div_mat,ct_sparsity,(/1,dim/))
    call zero(div_mat)
    call allocate(wave_mat,wave_sparsity)
    call zero(wave_mat)
    call allocate(big_mat,u_sparsity,(/dim,dim/))
    call zero(big_mat)

    call allocate(inverse_big_mat,u_sparsity,(/dim,dim/))
    call zero(big_mat)

    !Assemble matrices
    do ele = 1, ele_count(D)
       call assemble_shallow_water_matrices_ele(D,U,X,ele, &
            h_mass_mat,u_mass_mat,coriolis_mat,div_mat,big_mat,&
            inverse_big_mat,f0,beta,dt,theta,D0)
    end do

    if(have_option("/debug/check_inverse_coriolis_matrix")) then
       call check_big_mat(U,big_mat,u_mass_mat,coriolis_mat,theta,dt)
    end if

    !Construct wave mat
    call assemble_cmc_dg(wave_mat, div_mat, div_mat, big_mat)

    if(have_option("/debug/check_wave_matrix")) then
      call check_wave_mat(wave_mat,div_mat,big_mat,u,d)
    end if

    call scale(wave_mat,dt*dt*theta*theta*g*D0)
    call addto(wave_mat,h_mass_mat)

    call deallocate(inverse_big_mat)

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
      type(block_csr_matrix), intent(in) :: big_mat
      type(csr_matrix), intent(in) :: coriolis_mat, u_mass_mat
      real, intent(in) :: theta, dt
      !
      integer :: n,ntests=20, dim,d1
      type(vector_field) :: u_test, u_mem1, u_mem2
      type(scalar_field), dimension(:), allocatable :: u_mem1_cpts,&
           & u_test_cpts, u_mem2_cpts

      dim = mesh_dim(u)

      call allocate(u_test,dim,U%mesh,'testmem')
      call allocate(u_mem1,dim,U%mesh,'testmem1')
      call allocate(u_mem2,dim,U%mesh,'testmem2')

      allocate(u_mem1_cpts(dim),u_test_cpts(dim),u_mem2_cpts(dim))

      do d1 = 1, dim
         u_mem1_cpts(d1) = extract_scalar_field(u_mem1,d1)
         u_test_cpts(d1) = extract_scalar_field(u_test,d1)
         u_mem2_cpts(d1) = extract_scalar_field(u_mem2,d1)
      end do

      do n = 1, ntests
         do d1 = 1, dim
            call random_number(u_test%val(d1)%ptr)
         end do

         !First apply mass matrix
         do d1 = 1, dim
            call mult(u_mem1_cpts(d1),u_mass_mat,u_test_cpts(d1))
         end do
         !Then coriolis
         if(dim==2) then
            call mult(u_mem2_cpts(1),coriolis_mat,u_test_cpts(1))
            call mult(u_mem2_cpts(2),coriolis_mat,u_test_cpts(2))
            call addto(u_mem1_cpts(1),u_mem2_cpts(2),scale=-dt*theta)
            call addto(u_mem1_cpts(2),u_mem2_cpts(1),scale=dt*theta)
         end if
         !The application of the inverse of big_mat is now in u_mem1
         !Now apply big_mat to u_test and put it in u_mem2
         call mult(u_mem2,big_mat,u_mem1)
         !Now subtract u_test from it
         call addto(u_mem2,u_test,scale=-1.0)
         !Now check the residual
         do d1 = 1, dim
            if(maxval(abs(u_mem2%val(d1)%ptr))>1.0e-7) then
               ewrite(0,*) 'residual of ',maxval(abs(u_mem2%val(d1)%ptr))
               ewrite(0,*) 'in component',d1
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
      type(csr_matrix), intent(in) :: u_mass_mat,h_mass_mat,coriolis_mat
      type(block_csr_matrix), intent(in) :: div_mat
      real, intent(in) :: dt,theta,g,D0
      !
      type(scalar_field) :: h_residual, h_mem
      type(vector_field) :: u_residual, u_mem
      type(scalar_field), dimension(:), allocatable :: u_res_cpts, &
           & u_mem_cpts, u_cpts, delta_u_cpts
      integer :: dim, d1
      !
      dim = mesh_dim(u)
      !
      call allocate(h_residual,d%mesh,'hresidual')
      call allocate(u_residual,mesh_dim(u),u%mesh,'uresidual')
      call allocate(h_mem,d%mesh,'hmem')
      call allocate(u_mem,mesh_dim(u),u%mesh,'umem')
      !
      allocate(u_res_cpts(dim),u_mem_cpts(dim),u_cpts(dim),delta_u_cpts(dim))
      do d1 = 1, dim
         u_res_cpts(d1) = extract_scalar_field(u_residual,d1)
         u_mem_cpts(d1) = extract_scalar_field(u_mem,d1)
         u_cpts(d1) = extract_scalar_field(u,d1)
         delta_u_cpts(d1) = extract_scalar_field(delta_u,d1)
      end do
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
         call mult(u_mem_cpts(1),coriolis_mat,delta_u_cpts(1))
         call mult(u_mem_cpts(2),coriolis_mat,delta_u_cpts(2))
         call addto(u_residual,1,u_mem_cpts(2),scale=-dt*theta)
         call addto(u_residual,2,u_mem_cpts(1),scale= dt*theta)
         ! Coriolis explicit
         call mult(u_mem_cpts(1),coriolis_mat,u_cpts(1))
         call mult(u_mem_cpts(2),coriolis_mat,u_cpts(2))
         call addto(u_residual,1,u_mem_cpts(2),scale=-dt)
         call addto(u_residual,2,u_mem_cpts(1),scale= dt)
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
      ewrite(2,*) 'SW: u1 residual', maxval(abs(u_residual%val(1)%ptr))
      ewrite(2,*) 'SW: u size', maxval(abs(u%val(1)%ptr))
      if(dim==2) then
         ewrite(2,*) 'SW: u2 residual', maxval(abs(u_residual%val(2)%ptr))
      ewrite(2,*) 'SW: u size', maxval(abs(u%val(2)%ptr))
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

    subroutine get_d_rhs(d_rhs,u_rhs,D,U,h_mass_mat,div_mat,big_mat,D0,dt,theta)
      implicit none
      type(scalar_field), intent(inout) :: d_rhs
      type(vector_field), intent(inout) :: u_rhs,U
      type(scalar_field), intent(inout) :: D
      type(csr_matrix), intent(in) :: h_mass_mat
      type(block_csr_matrix), intent(in) :: div_mat, big_mat
      real, intent(in) :: D0, dt, theta
      !Construct the contribution to the right-hand side of the D wave eqn
      type(scalar_field) :: rhs2
      integer :: d1,d2,dim
      type(scalar_field) :: U_cpt
      type(csr_matrix) :: M_cpt
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
         coriolis_mat,div_mat)
      implicit none
      type(vector_field), intent(inout) :: u_rhs
      type(vector_field), intent(inout) :: U
      type(scalar_field), intent(inout) :: D
      type(csr_matrix), intent(in) ::  coriolis_mat
      type(block_csr_matrix), intent(in) :: div_mat
      real, intent(in) :: dt, g
      !Construct the explicit contribution to the right-hand side 
      !of the u equation

      type(scalar_field) :: rhs2
      integer :: d1,d2, dim
      type(scalar_field) :: U_cpt
      type(vector_field) :: vec

      ewrite(2,*) 'outside', sum(D%val)

      dim = mesh_dim(U)

      call allocate(rhs2, U%mesh, "workingmem")
      call allocate(vec, dim, U%mesh, "workingvecmem")

      call zero(u_rhs)

      !Coriolis term
      if (dim==2) then
         U_cpt = extract_scalar_field(U,2)
         call mult(rhs2,coriolis_mat,U_cpt)
         call addto(U_rhs,1,rhs2,scale=dt)

         U_cpt = extract_scalar_field(U,1)
         call mult(rhs2,coriolis_mat,U_cpt)
         call addto(U_rhs,2,rhs2,scale=-dt)
      end if

      !Pressure gradient
      call mult_T(vec,div_mat,D)
      call addto(u_rhs,vec,scale=-dt*g)

      ewrite(1,*) 'SW u_rhs',maxval(abs(u_rhs%val(1)%ptr))
      ewrite(1,*) 'SW u_rhs',maxval(abs(u_rhs%val(2)%ptr))

      call deallocate(rhs2)
      call deallocate(vec)
    end subroutine get_u_rhs

    subroutine assemble_shallow_water_matrices_ele(D,U,X,ele, &
         h_mass_mat,u_mass_mat,coriolis_mat,div_mat,big_mat,inverse_big_mat&
         &,f0,beta,dt,theta,D0)
      implicit none
      type(scalar_field), intent(in) :: D
      type(vector_field), intent(in) :: U,X
      type(csr_matrix), intent(inout) :: coriolis_mat, h_mass_mat, u_mass_mat
      type(block_csr_matrix), intent(inout) :: div_mat, big_mat
      integer, intent(in) :: ele
      type(block_csr_matrix), intent(inout), optional :: inverse_big_mat
      real, intent(in) :: f0,dt,theta,D0
      real, intent(in), dimension(:) :: beta
      !Assemble h_mass_mat, u_mass_mat, coriolis_mat, div_mat
      !and then big_mat
      real, dimension(ele_ngi(D, ele)) :: detwei
      real, dimension(ele_loc(D, ele), ele_ngi(D, ele), mesh_dim(D)) :: dm_t
      integer, dimension(:), pointer :: D_ele, U_ele
      type(element_type) :: D_shape, u_shape
      real, dimension(ele_ngi(D, ele)) :: f_gi
      real, dimension(mesh_dim(U), ele_ngi(x,ele)) :: x_gi
      real, dimension(mesh_dim(U)*ele_loc(U,ele), &
           mesh_dim(U)*ele_loc(U,ele)) :: l_big_mat
      integer :: dim1, dim2, nloc, dim
      real, dimension(mesh_dim(U),ele_loc(D,ele),ele_loc(U,ele)) :: l_div_mat

      dim = mesh_dim(U)

      u_shape=ele_shape(u, ele)
      d_shape=ele_shape(d, ele)

      D_ele => ele_nodes(D, ele)
      U_ele => ele_nodes(U, ele)

      x_gi = ele_val_at_quad(X,ele)
      f_gi = f0 + matmul(beta,x_gi)

      call transform_to_physical(X,ele,&
           & D_shape , dshape=dm_t, detwei=detwei)

      call addto(h_mass_mat, D_ele, D_ele, &
           shape_shape(d_shape,d_shape,detwei))
      call addto(u_mass_mat, u_ele, u_ele, &
           shape_shape(u_shape,u_shape,detwei))
      call addto(coriolis_mat,u_ele,u_ele, &
              shape_shape(u_shape,u_shape,detwei*f_gi))

      !if(present(inverse_big_mat)) then
      !   do dim1 = 1, dim
      !      call addto(inverse_big_mat,dim1,dim1,shape_shape(u_shape,
      !end if

      l_div_mat = dshape_shape(dm_t,u_shape,detwei)

      do dim1 = 1, dim
         call addto(div_mat,1,dim1,d_ele,u_ele,l_div_mat(dim1,:,:))
      end do

      l_big_mat = 0.

      nloc = ele_loc(U,ele)
      do dim1 = 1, dim
         !Mass
         l_big_mat(nloc*(dim1-1)+1:nloc*dim1, &
              nloc*(dim1-1)+1:nloc*dim1) = shape_shape(&
              u_shape,u_shape,detwei)
      end do
      if(dim==2) then
         !Coriolis
         l_big_mat(1:nloc, nloc+1:) = -shape_shape(&
              u_shape,u_shape,detwei*f_gi*dt*theta)
         l_big_mat(nloc+1:, 1:nloc) = shape_shape(&
              u_shape,u_shape,detwei*f_gi*dt*theta)
      end if

      call invert(l_big_mat)

      do dim1 = 1, dim
         do dim2 = 1, dim
            call addto(big_mat,dim1,dim2,u_ele,u_ele, &
                 l_big_mat(nloc*(dim1-1)+1:nloc*dim1, &
                 nloc*(dim2-1)+1:nloc*dim2))
         end do
      end do

    end subroutine assemble_shallow_water_matrices_ele

end module linear_shallow_water
