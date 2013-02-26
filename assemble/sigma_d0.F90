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
#endif

    implicit none

    private
    public add_sigma_element_cg,add_sigma_element_dg
    logical :: move_mesh
    logical :: lump_mass    
    logical :: assemble_mass_matrix
    logical :: assemble_inverse_masslump
    logical :: have_wd
    logical :: cmc_lump_mass
    logical :: exclude_mass
    ! implicitness parameter, timestep
    real :: dt
    
    
 contains
 subroutine calculate_sigma_element(state, ele, x, u, oldu_val, detwei, detwei_old, detwei_new, sigma_mat, sigma_move_mat, sigma_lump,sigma_move_lump)
   type(state_type), intent(inout) :: state
   integer, intent(in) :: ele
   type(vector_field), intent(in) :: u, x
   type(scalar_field), pointer :: p
   type(scalar_field), pointer:: z
   type(scalar_field), pointer::sigma
   real, dimension(:,:), intent(in) :: oldu_val
   real, dimension(ele_ngi(u, ele)), intent(in) :: detwei, detwei_old, detwei_new
   integer, parameter:: dim=3
   integer, parameter:: xloc=2
   real :: a
   real, dimension(ele_loc(u, ele)),intent(out):: sigma_lump, sigma_move_lump
   real, dimension(ele_loc(u, ele), ele_loc(u, ele)),intent(out) :: sigma_mat, sigma_move_mat
   type(element_type), pointer :: u_shape
    ! In case we have to multiply detwei by various coefficients (e.g. the density values at the Gauss points), 
      ! then place the result in here
   real, dimension(ele_ngi(u, ele)) :: coefficient_detwei
   real, dimension(:,:), allocatable :: x_ele
   real, dimension(ele_ngi(u, ele)):: sigma_ele
   integer, dimension(:), pointer:: u_ele
   real :: dz_ele
   real, dimension(dim,dim)  :: dx
   real ::length1, length2, length3, dx_ele
   integer :: i
      p=>extract_scalar_field(state,"Pressure")
      call get_option("/timestepping/timestep", dt)
      move_mesh = (have_option("/mesh_adaptivity/mesh_movement").and.(.not.exclude_mass))
      lump_mass=have_option(trim(u%option_path)//&
          &"/prognostic/spatial_discretisation"//&
          &"/continuous_galerkin/mass_terms/lump_mass_matrix")
      assemble_mass_matrix = have_option(trim(p%option_path)//&
          &"/prognostic/scheme/use_projection_method"//&
          &"/full_schur_complement/inner_matrix::FullMassMatrix")
      cmc_lump_mass = have_option(trim(p%option_path)//&
          &"/prognostic/scheme"//&
          &"/use_projection_method/full_schur_complement"//&
          &"/preconditioner_matrix::LumpedSchurComplement")
      assemble_inverse_masslump = lump_mass .or. cmc_lump_mass
      assemble_mass_matrix = have_option(trim(p%option_path)//&
          &"/prognostic/scheme/use_projection_method"//&
          &"/full_schur_complement/inner_matrix::FullMassMatrix")
      exclude_mass = have_option(trim(u%option_path)//&
          &"/prognostic/spatial_discretisation"//&
          &"/continuous_galerkin/mass_terms/exclude_mass_terms")
      have_wd=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying") 
      
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
      sigma=> extract_scalar_field(state,"Sigma_d0")
      sigma_ele = dx_ele**2/(a**2*dz_ele**2)
      call set(sigma, u_ele, sigma_ele)
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
     
   subroutine add_sigma_element_cg(state, ele, x, u, oldu_val, detwei, detwei_old, detwei_new,big_m_diag_addto, big_m_tensor_addto, rhs_addto, mass, masslump)
      type(state_type), intent(inout) :: state
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: u,x
      real, dimension(:,:), intent(in) :: oldu_val
      integer, dimension(:), pointer:: u_ele
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: big_m_diag_addto
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      type(petsc_csr_matrix), intent(inout) :: mass
      type(vector_field), intent(inout) :: masslump
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei, detwei_old, detwei_new
      real, dimension(ele_loc(u, ele)) :: sigma_lump, sigma_move_lump
      real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: sigma_mat, sigma_move_mat
      call calculate_sigma_element(state, ele, x, u, oldu_val, detwei, detwei_old, detwei_new, sigma_mat, sigma_move_mat, sigma_lump,sigma_move_lump)
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
      
    subroutine add_sigma_element_dg(ele, x, u, p, sigma_scalar_field, oldu_val, detwei, detwei_old, detwei_new, big_m_diag_addto, big_m_tensor_addto, rhs_addto,mass, masslump)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: u,x
      type(scalar_field), intent(in) :: p, sigma_scalar_field
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
      call calculate_sigma_element(ele, x, u, oldu_val, detwei, detwei_old, detwei_new, sigma_mat, sigma_move_mat, sigma_lump,sigma_move_lump)
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
  end module sigma_d0
