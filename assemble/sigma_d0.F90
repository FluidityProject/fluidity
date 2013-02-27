#include "fdebug.h"
 module sigma_d0
  !This module only works in this case: 1)3-dimentional case, 2)z direction in the same direction as gravity, 3)mesh element is triangle.
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
    use smoothing_module
    use metric_tools
    use field_derivatives
    use state_fields_module
    use state_matrices_module
    use sparsity_patterns_meshes
    use fefields
    use fields_manipulation
    use rotated_boundary_conditions
    use Coordinates
    use edge_length_module
    use colouring
    use Profiler
#ifdef _OPENMP
    use omp_lib
    
    use fields_manipulation
#endif

    implicit none

    private
    public add_sigma_element_cg,add_sigma_element_dg, calculate_diagnostic_sigma_d0
    
 contains
 subroutine calculate_sigma_element(ele, x, u, detwei, detwei_old, detwei_new, sigma_mat, sigma_move_mat, sigma_lump,sigma_move_lump,move_mesh, lump_mass,assemble_mass_matrix)
   integer, intent(in) :: ele
   type(vector_field), intent(in) :: u, x
   type(scalar_field), pointer:: z
   real, dimension(ele_ngi(u, ele)), intent(in) :: detwei, detwei_old, detwei_new
   integer, parameter:: dim=3
   integer, parameter:: xloc=2
   real :: a
   real, dimension(ele_loc(u, ele)),intent(inout):: sigma_lump, sigma_move_lump
   real, dimension(ele_loc(u, ele), ele_loc(u, ele)),intent(inout) :: sigma_mat, sigma_move_mat
   type(element_type), pointer :: u_shape
   real, dimension(ele_ngi(u, ele)) :: coefficient_detwei
   real, dimension(:,:), allocatable :: x_ele
   real, dimension(ele_ngi(u, ele)):: sigma_ele
   integer, dimension(:), pointer:: u_ele
   real :: dz_ele
   real, dimension(dim,dim)  :: dx
   real ::length1, length2, length3, dx_ele
   real :: dt
   logical,intent(in) :: move_mesh
   logical,intent(in) :: lump_mass
   logical,intent(in):: assemble_mass_matrix
      sigma_lump=0
      sigma_move_lump=0
      sigma_mat=0
      sigma_move_mat=0
      allocate(x_ele(x%dim,ele_loc(x,ele)))
      !get the nodes cordinate of element
      x_ele = ele_val(x, ele) 
      u_shape => ele_shape(u, ele)
      u_ele => ele_nodes(u, ele)
      
      !calculate the lateral lenghth of element's projection to horizontal surface
      length1 = sqrt((x_ele(2,2)-x_ele(2,1))**2+(x_ele(1,2)-x_ele(1,1))**2)
      length1 = sqrt((x_ele(2,3)-x_ele(2,1))**2+(x_ele(1,3)-x_ele(1,1))**2)
      length1 = sqrt((x_ele(2,2)-x_ele(2,3))**2+(x_ele(1,2)-x_ele(1,3))**2)
      dx_ele = max(length1,length2,length3)
      dz_ele = max(x_ele(3,1),x_ele(3,2), x_ele(3,3) )-min(x_ele(3,1), x_ele(3,2), x_ele(3,3))
      print * , "dz=",dz_ele
      print * , "dx=",dx_ele
      a=7
     ! call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/a", a)
     ! ewrite(2,*) " tolerable aspect ratio a = ", a
     ! sigma = dx_ele**2/(a**2*dt*dz_ele**2).But when we add it to mass and other terms, we need to multiply dt, so here we get the result after multiplying dt.
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
      if(move_mesh.and.assemble_mass_matrix) then
          sigma_move_mat = shape_shape(u_shape, u_shape, (detwei_new-detwei_old)*sigma_ele)  
		  if(lump_mass) then
		  	 sigma_move_lump=sum(sigma_move_mat,2)
		  end if
      end if
   end subroutine calculate_sigma_element
     
   subroutine add_sigma_element_cg( ele, x, u, oldu_val, detwei, detwei_old, detwei_new,big_m_diag_addto, big_m_tensor_addto, rhs_addto, mass, masslump,dt, move_mesh,lump_mass,assemble_mass_matrix, assemble_inverse_masslump ,exclude_mass,have_wd)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: u,x
      real, dimension(:,:), intent(in) :: oldu_val
      integer, dimension(:), pointer:: u_ele
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: big_m_diag_addto
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      type(petsc_csr_matrix), intent(inout) :: mass
      type(vector_field), intent(inout) :: masslump
      real, dimension(ele_ngi(u, ele)), intent(inout) :: detwei, detwei_old, detwei_new
      real, dimension(ele_loc(u, ele)) :: sigma_lump, sigma_move_lump
      real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: sigma_mat, sigma_move_mat
      real :: dt
      logical,intent(inout) :: move_mesh
      logical,intent(in) :: lump_mass    
      logical,intent(in) :: assemble_mass_matrix
      logical ,intent(in):: assemble_inverse_masslump
      logical ,intent(in):: have_wd
      logical,intent(in) :: exclude_mass
      call calculate_sigma_element(ele, x, u, detwei, detwei_old, detwei_new, sigma_mat, sigma_move_mat, sigma_lump,sigma_move_lump,move_mesh, lump_mass,assemble_mass_matrix)
      if(lump_mass) then
             big_m_diag_addto(3, :) = big_m_diag_addto(3, :) + sigma_lump      
        else
       	     big_m_tensor_addto(3, 3, :, :) = big_m_tensor_addto(3, 3, :, :) + sigma_mat  
       	end if
      if(assemble_inverse_masslump) then
         call addto(masslump, 3, u_ele, sigma_lump)
      end if 
      if(assemble_mass_matrix) then
        call addto(mass,3,3,u_ele, u_ele, sigma_mat)
      end if
      if(move_mesh) then
        if(lump_mass) then
          rhs_addto(3,:) = rhs_addto(3,:) - sigma_move_lump*oldu_val(3,:)/dt           
        else
          rhs_addto(3,:) = rhs_addto(3,:)+matmul(sigma_move_mat, oldu_val(3,:))/dt
        end if       
      end if
      print *, "rhs_addto matrix:",rhs_addto

    end subroutine add_sigma_element_cg
      
    subroutine add_sigma_element_dg(ele, x, u, oldu_val, detwei, detwei_old, detwei_new, big_m_diag_addto, big_m_tensor_addto, rhs_addto,mass, masslump,dt, move_mesh,lump_mass,assemble_mass_matrix, assemble_inverse_masslump ,exclude_mass,have_wd)
      integer :: ele
      type(vector_field), intent(in) :: u,x
      real, dimension(:,:), intent(in) :: oldu_val
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: big_m_diag_addto
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      type(csr_matrix), intent(inout), optional :: mass
      real, dimension(ele_loc(u,ele)) :: masslump
      integer, dimension(:), pointer:: u_ele
      real, dimension(ele_loc(u, ele)) :: sigma_lump, sigma_move_lump
      real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: sigma_mat, sigma_move_mat
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei, detwei_old, detwei_new
      real :: dt
      logical,intent(in) :: move_mesh
      logical,intent(in) :: lump_mass    
      logical,intent(in) :: assemble_mass_matrix
      logical ,intent(in):: assemble_inverse_masslump
      logical ,intent(in):: have_wd
      logical,intent(in) :: exclude_mass
      call calculate_sigma_element(ele, x, u,detwei, detwei_old, detwei_new, sigma_mat, sigma_move_mat, sigma_lump,sigma_move_lump,move_mesh, lump_mass,assemble_mass_matrix)
      if(lump_mass) then
             big_m_diag_addto(3, :) = big_m_diag_addto(3, :) + sigma_lump  
      else
       	     big_m_tensor_addto(3, 3, :, :) = big_m_tensor_addto(3, 3, :, :) + sigma_mat  
      end if
     if(assemble_inverse_masslump) then
         call addto(masslump, u_ele, sigma_lump)
     end if 
      if(assemble_mass_matrix) then
        call addto(mass,u_ele, u_ele, sigma_mat)
      end if
      if(move_mesh) then
        ! Put the sigma tilder u^n+1 term on the rhs
        if(lump_mass) then
          sigma_move_lump=sum(sigma_move_mat,2) 
          rhs_addto(3,:) = rhs_addto(3,:) - sigma_move_lump*oldu_val(3,:)/dt           
        else
          rhs_addto(3,:) = rhs_addto(3,:)+matmul(sigma_move_mat, oldu_val(3,:))/dt
        end if 
      end if
      print *, "rhs_addto matrix:", rhs_addto
    end subroutine add_sigma_element_dg  
 
 subroutine calculate_diagnostic_sigma_d0(state, sigma)
    type(state_type),intent(in) :: state
    integer :: ele
    type(scalar_field), intent(inout) :: sigma
    integer, dimension(ele_loc(sigma, 1)) ::ele_sigma
    real,dimension(ele_loc(sigma, 1)) ::sigma_node 
    type(vector_field) :: U, X
    type(scalar_field), pointer:: z
    !real, dimension(ele_ngi(sigma, 1)) :: detwei
    !real, dimension(ele_ngi(u, ele)) :: coefficient_detwei
    integer, parameter:: dim=3
    integer, parameter:: xloc=2
    real :: a 
   type(element_type), pointer :: u_shape
   real, dimension(:,:), allocatable :: x_ele
   real, dimension(:), allocatable:: sigma_ele
   real :: dz_ele
   real, dimension(dim,dim)  :: dx
   real ::length1, length2, length3, dx_ele
   ! Inverse of the local coordinate change matrix.
    ! sigma value at each quad point.
    real, dimension(mesh_dim(sigma), ele_ngi(sigma, 1)) :: sigma_q
    ! current element global node numbers.
    ! local sigma matrix on the current element.
    real, dimension(ele_loc(sigma, 1),ele_loc(sigma, 1)) :: sigma_mat
    ! current sigma element shape
    type(element_type), pointer :: sigma_shape
    U=extract_vector_field(state, "Velocity")
    X=extract_vector_field(state, "Coordinate")
	!sigma=> extract_scalar_field(state,"Sigma_d0") sigma is extracted in diagnostic_fields_wrapper
    call zero(sigma) 
    do ele=1,element_count(sigma)
       ele_sigma=ele_nodes(sigma,ele)
       sigma_shape=>ele_shape(sigma,ele)
      allocate(x_ele(x%dim,ele_loc(x,ele)))
      allocate(sigma_ele(ele_ngi(u, ele)))

      !get the nodes cordinate of element
      x_ele = ele_val(x, ele) 
      u_shape => ele_shape(u, ele)
      
      !calculate the lateral lenghth of element's projection to horizontal surface
      length1 = sqrt((x_ele(2,2)-x_ele(2,1))**2+(x_ele(1,2)-x_ele(1,1))**2)
      length1 = sqrt((x_ele(2,3)-x_ele(2,1))**2+(x_ele(1,3)-x_ele(1,1))**2)
      length1 = sqrt((x_ele(2,2)-x_ele(2,3))**2+(x_ele(1,2)-x_ele(1,3))**2)
      dx_ele = max(length1,length2,length3)
      dz_ele = max(x_ele(3,1),x_ele(3,2), x_ele(3,3) )-min(x_ele(3,1), x_ele(3,2), x_ele(3,3))
      print * , "dz=",dz_ele
      print * , "dx=",dx_ele
      a=7
     ! call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/a", a)
     ! ewrite(2,*) " tolerable aspect ratio a = ", a
     ! sigma = dx_ele**2/(a**2*dt*dz_ele**2).But when we add it to mass and other terms, we need to multiply dt, so here we get the result after multiplying dt.
      !sigma_ele = dx_ele**2/(a**2*dz_ele**2)
      sigma_ele=1
      sigma_node =sigma_ele
      !coefficient_detwei = sigma_ele*detwei
     ! sigma_mat = shape_shape(sigma_shape, sigma_shape, coefficient_detwei) 
      print *, "sigma matrix:", sigma_mat 
      !sigma%val(ele_sigma)=(sigma%val(ele_sigma)+sum(sigma_mat,2))/2
      !call set(sigma,ele_sigma,sigma_ele)
      
      call set(sigma,ele_sigma,sigma_node)
      !sigma%val(ele_sigma)=1
      deallocate(x_ele)    
    end do 
 end subroutine calculate_diagnostic_sigma_d0
 
 end module sigma_d0
