#include "fdebug.h"
   module sigma_d0
 
    use fields
    use state_module
    use spud
    use fldebug
    use sparse_tools
    use boundary_conditions
    use boundary_conditions_from_options
    use solvers
    use petsc_solve_state_module
    use sparse_tools_petsc
    use sparse_matrices_fields
    use field_options
    use halos
    use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN, timestep, &
         COLOURING_CG1
    use elements
    use transform_elements, only: transform_to_physical
    use coriolis_module
    use vector_tools
    use fetools
    use upwind_stabilisation
    use les_viscosity_module
    use smoothing_module
    use metric_tools
    use field_derivatives
    use state_fields_module
    use state_matrices_module
    use sparsity_patterns_meshes
    use fefields
    use rotated_boundary_conditions
    use Coordinates
    use multiphase_module
    use edge_length_module
    use colouring
    use Profiler
#ifdef _OPENMP
    use omp_lib
#endif

    implicit none

    private
    !public :: construct_momentum_cg, correct_masslumped_velocity, &
             ! correct_velocity_cg, assemble_masslumped_poisson_rhs, &
              !add_kmk_matrix, add_kmk_rhs, assemble_kmk_matrix, &
             ! deallocate_cg_mass, assemble_poisson_rhs
    public add_sigma_element_cg,add_sigma_element_dg
    logical :: move_mesh
    logical :: lump_mass    
    logical :: assemble_mass_matrix
    logical :: assemble_inverse_masslump
    logical :: have_wd
    ! implicitness parameter, timestep, conservation parameter
    real :: dt
    
 contains
 subroutine add_sigma_element_cg(ele, test_function, x, u, oldu_val, detwei, detwei_old, detwei_new, big_m_diag_addto, big_m_tensor_addto, rhs_addto,mass, masslump)
    !This subroutine only works in this case: 1)3-dimentiona, 2)z direction in the same direction as gravity, 3)mesh element is triangle.
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: u, x
      type(scalar_field), pointer:: z
      real, dimension(:,:), intent(in) :: oldu_val
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei, detwei_old, detwei_new
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: big_m_diag_addto
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      integer, parameter:: dim=3
      integer, parameter:: xloc=2
      integer :: a
      integer ::i,j
     ! type(vector_field), intent(inout) :: visc_inverse_masslump
      type(petsc_csr_matrix), intent(inout) :: mass
      type(vector_field), intent(inout) :: masslump
      real, dimension(ele_loc(u, ele)) :: sigma_lump, sigma_move_lump
      real, dimension(:,:), allocatable :: x_ele
      integer, dimension(:), pointer :: nodes
   
      real, dimension(ele_ngi(u, ele)):: sigma_ele
      integer, dimension(:), pointer:: u_ele
      real :: dz_ele
      type(petsc_csr_matrix):: sigma
      real, dimension(dim,dim)  :: dx
      real ::length1, length2, length3, dx_ele
      real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: sigma_mat, sigma_move_mat
      type(element_type), pointer :: u_shape
      ! In case we have to multiply detwei by various coefficients (e.g. the density values at the Gauss points), 
      ! then place the result in here
      real, dimension(ele_ngi(u, ele)) :: coefficient_detwei
      allocate(x_ele(x%dim,ele_loc(x,ele)))
      
      have_wd=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying")
      x_ele = ele_val(x, ele) !get the nodes cordinate of element
      u_shape => ele_shape(u, ele)
      u_ele => ele_nodes(u, ele)
   
      !assert(size(x_ele,1)==dim)
      !assert(size(x_ele,2)==xloc)
       print *, " oldu_val(2,:)", oldu_val(2,:)
       print *, " oldu_val(3,:)", oldu_val(3,:)
      !dx(:,1)=x_ele(:,2)-x_ele(:,1)
      !dx(:,2)=x_ele(:,3)-x_ele(:,1)
      !dx(:,3)=x_ele(:,2)-x_ele(:,3)
      !length1 = sqrt(dot_product(dx(:,1),dx(:,1)))
      !length2 = sqrt(dot_product(dx(:,2),dx(:,2)))
      !length3 = sqrt(dot_product(dx(:,3),dx(:,3)))
      !calculate the lateral lenghth of element's projection to horizontal surface
      length1 = sqrt((x_ele(2,2)-x_ele(2,1))**2+(x_ele(1,2)-x_ele(1,1))**2)
      length1 = sqrt((x_ele(2,3)-x_ele(2,1))**2+(x_ele(1,3)-x_ele(1,1))**2)
      length1 = sqrt((x_ele(2,2)-x_ele(2,3))**2+(x_ele(1,2)-x_ele(1,3))**2)
      dx_ele = max(length1,length2,length3)
      dz_ele = max(x_ele(3,1),x_ele(3,2), x_ele(3,3) )-min(x_ele(3,1), x_ele(3,2), x_ele(3,3))
      print * , "dz=",dz_ele
      print * , "dx=",dx_ele
      !sigma = dx_ele**2/(a**2*dt*dz_ele**2).But when we add it to mass and other terms, we need to multiply dt, so here we get the result after multiplying dt.
      sigma_ele = dx_ele**2/(a**2*dz_ele**2)
      deallocate(x_ele)    
      if(move_mesh) then
         sigma_mat = shape_shape(test_function, u_shape, sigma_ele*detwei_new)
         
      else
         coefficient_detwei = sigma_ele*detwei
         sigma_mat = shape_shape(test_function, u_shape, coefficient_detwei) 
         print *, "sigma matrix:", sigma_mat 
      end if
      sigma_lump = sum(sigma_mat, 2) 
      
      if(have_wd) then
        if(lump_mass) then
             big_m_diag_addto(3, :) = big_m_diag_addto(3, :) + sigma_lump
             
        else
       	     big_m_tensor_addto(3, 3, :, :) = big_m_tensor_addto(3, 3, :, :) + sigma_mat  
       	end if
      end if
      if(assemble_inverse_masslump) then
         call addto(masslump, 3, u_ele, sigma_lump)
      end if 
      if(assemble_mass_matrix) then
        call addto(mass,3,3,u_ele, u_ele, sigma_mat)
      end if
      
      if(move_mesh.and.assemble_mass_matrix) then
        ! Put the sigma tilder u^n+1 term on the rhs.
        sigma_move_mat = shape_shape(test_function, u_shape, (detwei_new-detwei_old)*sigma_ele)  
        if(lump_mass) then
          sigma_move_lump=sum(sigma_move_mat,2) 
          rhs_addto(3,:) = rhs_addto(3,:) - sigma_move_lump*oldu_val(3,:)/dt           
        else
          rhs_addto(3,:) = rhs_addto(3,:)+matmul(sigma_move_mat, oldu_val(3,:))/dt
        end if
       
      end if
      
     !nodes => ele_nodes(u,ele)
       !do i = 1, ele_loc(u,ele)
         ! call addto(visc_inverse_masslump, 3, nodes(i), sigma_mat(i,i))
      ! end do
      print *, "rhs_addto matrix:"
      do i=1, u%dim
        do j=1,ele_loc(u, ele)
          print *, rhs_addto(i,j)
        end do
      end do

      end subroutine add_sigma_element_cg
      
  subroutine add_sigma_element_dg(ele, u_shape, x, u, oldu_val, detwei, detwei_old, detwei_new, big_m_diag_addto, big_m_tensor_addto, rhs_addto,mass, masslump)
    !This subroutine only works in this case: 1)3-dimentional case, 2)z direction in the same direction as gravity, 3)mesh element is triangle.
      integer, intent(in) :: ele
      !type(state_type), intent(inout) :: state
      type(vector_field), intent(in) :: u, x
      type(scalar_field), pointer:: z
      real, dimension(:,:), intent(in) :: oldu_val
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei, detwei_old, detwei_new
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: big_m_diag_addto
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      integer, parameter:: dim=3
      integer, parameter:: xloc=2
      integer :: a
      integer ::i,j
     ! type(vector_field), intent(inout) :: visc_inverse_masslump
      type(csr_matrix), intent(inout), optional :: mass
      real, dimension(ele_loc(u,ele)) :: masslump
      real, dimension(ele_loc(u, ele)) :: sigma_lump, sigma_move_lump
      real, dimension(:,:), allocatable :: x_ele
      integer, dimension(:), pointer :: nodes
      real, dimension(ele_ngi(u, ele)):: sigma_ele
      integer, dimension(:), pointer:: u_ele
      real :: dz_ele
      type(petsc_csr_matrix):: sigma
      real, dimension(dim,dim)  :: dx
      real ::length1, length2, length3, dx_ele
      real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: sigma_mat, sigma_move_mat
      type(element_type), pointer :: u_shape
      ! In case we have to multiply detwei by various coefficients (e.g. the density values at the Gauss points), 
      ! then place the result in here
      real, dimension(ele_ngi(u, ele)) :: coefficient_detwei
      allocate(x_ele(x%dim,ele_loc(x,ele)))
      
      have_wd=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying")
      x_ele = ele_val(x, ele) !get the nodes cordinate of element
      u_shape => ele_shape(u, ele)
      u_ele => ele_nodes(u, ele)
      !z => extract_scalar_field(state,"DistanceToBottom")
     ! allocate(z_ele(ele_loc(z,ele)))
     ! z_ele = ele_val(z, ele)
      !print *,"z_ele:", z_ele(1), z_ele(2), z_ele(3)
   
      !assert(size(x_ele,1)==dim)
      !assert(size(x_ele,2)==xloc)
       print *, " oldu_val(2,:)", oldu_val(2,:)
       print *, " oldu_val(3,:)", oldu_val(3,:)
      !dx(:,1)=x_ele(:,2)-x_ele(:,1)
      !dx(:,2)=x_ele(:,3)-x_ele(:,1)
      !dx(:,3)=x_ele(:,2)-x_ele(:,3)
      !length1 = sqrt(dot_product(dx(:,1),dx(:,1)))
      !length2 = sqrt(dot_product(dx(:,2),dx(:,2)))
      !length3 = sqrt(dot_product(dx(:,3),dx(:,3)))
      !calculate the lateral lenghth of element's projection to horizontal surface
      length1 = sqrt((x_ele(2,2)-x_ele(2,1))**2+(x_ele(1,2)-x_ele(1,1))**2)
      length1 = sqrt((x_ele(2,3)-x_ele(2,1))**2+(x_ele(1,3)-x_ele(1,1))**2)
      length1 = sqrt((x_ele(2,2)-x_ele(2,3))**2+(x_ele(1,2)-x_ele(1,3))**2)
      dx_ele = max(length1,length2,length3)
      dz_ele = max(x_ele(3,1),x_ele(3,2), x_ele(3,3) )-min(x_ele(3,1), x_ele(3,2), x_ele(3,3))
      print * , "dz=",dz_ele
      print * , "dx=",dx_ele
      !sigma = dx_ele**2/(a**2*dt*dz_ele**2).But when we add it to mass and other terms, we need to multiply dt, so here we get the result after multiplying dt.
      sigma_ele = dx_ele**2/(a**2*dz_ele**2)
      deallocate(x_ele)    
      if(move_mesh) then
         sigma_mat = shape_shape(u_shape, u_shape, sigma_ele*detwei_new)
         
      else
         coefficient_detwei = sigma_ele*detwei
         sigma_mat = shape_shape(u_shape, u_shape, coefficient_detwei) 
         print *, "sigma matrix:", sigma_mat 
      end if
      sigma_lump = sum(sigma_mat, 2) 
      
      if(have_wd) then
        if(lump_mass) then
             big_m_diag_addto(3, :) = big_m_diag_addto(3, :) + sigma_lump
             
        else
       	     big_m_tensor_addto(3, 3, :, :) = big_m_tensor_addto(3, 3, :, :) + sigma_mat  
       	end if
      end if
     if(assemble_inverse_masslump) then
         call addto(masslump, u_ele, sigma_lump)
     end if 
      if(assemble_mass_matrix) then
        call addto(mass,u_ele, u_ele, sigma_mat)
      end if
      
      if(move_mesh.and.assemble_mass_matrix) then
        ! Put the sigma tilder u^n+1 term on the rhs
        sigma_move_mat = shape_shape(u_shape, u_shape, (detwei_new-detwei_old)*sigma_ele)  
        if(lump_mass) then
          sigma_move_lump=sum(sigma_move_mat,2) 
          rhs_addto(3,:) = rhs_addto(3,:) - sigma_move_lump*oldu_val(3,:)/dt           
        else
          rhs_addto(3,:) = rhs_addto(3,:)+matmul(sigma_move_mat, oldu_val(3,:))/dt
        end if
       
      end if
      
     !nodes => ele_nodes(u,ele)
       !do i = 1, ele_loc(u,ele)
         ! call addto(visc_inverse_masslump, 3, nodes(i), sigma_mat(i,i))
      ! end do
      print *, "rhs_addto matrix:"
      do i=1, u%dim
        do j=1,ele_loc(u, ele)
          print *, rhs_addto(i,j)
        end do
      end do

      end subroutine add_sigma_element_dg
      
   end module sigma_d0
