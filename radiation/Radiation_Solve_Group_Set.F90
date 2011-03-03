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
!    but WITHOUT ANY WARRANTY; without even the implied arranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module radiation_assemble_solve_group_set

   !!< This module contains procedures associated with solving the set 
   !!< of groups g neutral particle block matrix problem

   ! keep in this order, please:
   use quadrature
   use elements
   use sparse_tools
   use fields

   use boundary_conditions
   use boundary_conditions_from_options
   use field_options   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   use petsc_solve_state_module   
   use sparsity_patterns_meshes 
    
   use radiation_materials
   use radiation_extract_flux_field
   use radiation_materials_interpolation 
     
   implicit none
   
   private 

   public :: np_assemble_solve_group_set
   
   type np_group_g_assemble_options_type
      character(len=OPTION_PATH_LEN) :: spatial_discretisation
      real :: theta
      real :: timestep
      logical :: have_delayed
      logical :: include_prescribed_source
      integer :: number_prescribed_source
      integer, dimension(:), allocatable :: prescribed_source_energy_group
      character(len=OPTION_PATH_LEN), dimension(:), allocatable :: prescribed_source_field_name
   end type np_group_g_assemble_options_type

contains

   ! --------------------------------------------------------------------------

   subroutine np_assemble_solve_group_set(state, &
                                          np_radmat_name, &
                                          np_radmat, &
                                          np_radmat_ii, &
                                          g, &
                                          number_of_energy_groups, &
                                          keff) 
   
      !!< Assemble and solve the group set g neutral particle block matrix problem
      !!< If keff variable is present then this is an eigenvalue solve, else this is a time run solve

      type(state_type), intent(inout) :: state
      character(len=*), intent(in) :: np_radmat_name
      type(np_radmat_type), intent(in) :: np_radmat
      type(np_radmat_ii_type), intent(in) :: np_radmat_ii      
      integer, intent(in) :: g
      integer, intent(in) :: number_of_energy_groups 
      real, intent(in), optional :: keff
      
      ! local variables
      integer :: vele,sele
      type(scalar_field), pointer :: np_flux 
      type(scalar_field), pointer :: np_flux_old
      type(csr_matrix) :: matrix
      type(scalar_field) :: rhs
      type(csr_sparsity), pointer :: sparsity
      type(vector_field), pointer :: positions      
      type(radmat_type) :: radmat_vele
      type(np_group_g_assemble_options_type) :: np_group_g_assemble_options
      type(scalar_field) :: delta_np_flux
      integer, dimension(:), allocatable :: bc_types
      type(scalar_field) :: bc_value
                 
      ! extract flux to be solved for 
      call extract_flux_group_g(state, &
                                trim(np_radmat_name), &
                                g, &  
                                np_flux = np_flux, &
                                np_flux_old = np_flux_old)
      
      call get_group_g_assemble_options(np_group_g_assemble_options, &
                                        np_flux%option_path, &
                                        np_radmat%option_path, &
                                        (.not. present(keff)))
      
      ! allocate the solution field passed to the solve for a time run 
      alloc_delta: if (.not. present(keff)) then
         
         call allocate(delta_np_flux, np_flux%mesh, name = 'DeltaNeutralParticleFluxGroup'//int2str(g)//trim(np_radmat_name))
         
         ! set the time run initial guess
         call zero(delta_np_flux)
         
      end if alloc_delta
      
      ! determine the sparsity pattern of matrix assuming first order connections 
      sparsity => get_csr_sparsity_firstorder(state, np_flux%mesh, np_flux%mesh)
      
      ! allocate the matrix and rhs vector used for solving
      call allocate(matrix, sparsity, name = 'MatrixNeutralParticleFluxGroup'//int2str(g)//trim(np_radmat_name))
      call allocate(rhs, np_flux%mesh, name = 'RHSNeutralParticleFluxGroup'//int2str(g)//trim(np_radmat_name))
      
      call zero(matrix)
      call zero(rhs)
      
      positions => extract_vector_field(state, "Coordinate")
      
      ! volume integrations 
      volume_element_loop: do vele = 1, ele_count(np_flux)      
         
         ! get the vele cross section as needed for this group
         call form(radmat_vele, &
                   np_radmat_ii, &
                   np_radmat, &
                   vele, &
                   form_diffusion       = .true. , &
                   form_removal         = .true. , &
                   form_scatter         = .true. , &
                   form_production      = .true. , &
                   form_prompt_spectrum = .true. , &
                   form_beta            = np_group_g_assemble_options%have_delayed)
         
         call assemble_np_group_g_vele(vele, &
                                       g, &
                                       number_of_energy_groups, &
                                       state, &
                                       trim(np_radmat_name), &
                                       np_flux, &
                                       np_flux_old, &
                                       radmat_vele, &
                                       matrix, &
                                       rhs, &
                                       positions, &
                                       np_group_g_assemble_options, &
                                       keff = keff)
         
         call destroy(radmat_vele)
         
      end do volume_element_loop
      
      ! surface integrations for vacuum conditions for diffusion theory
      
      allocate(bc_types(surface_element_count(np_flux)))

      call get_entire_boundary_condition(np_flux, &
                                         (/"vacuum"/), &
                                         bc_value, &
                                         bc_types)

      surface_element_loop: do sele = 1,surface_element_count(np_flux)
         
         vacuum_bc: if (bc_types(sele) == 1) then
         
            call assemble_np_group_g_sele_vacuum(sele, &
                                                 positions, &
                                                 matrix, & 
                                                 bc_value, &
                                                 np_flux)
         
         end if vacuum_bc
         
      end do surface_element_loop

      call deallocate(bc_value)

      if (allocated(bc_types)) deallocate(bc_types)
            
      ! apply direchlet BC
      call apply_dirichlet_conditions(matrix, &
                                      rhs, &
                                      np_flux)
            
      ! solve the linear system 
      ! - for eig solve use the latest solution as initial guess
      ! - for time solve update the state solution arrays after
      solve: if (present(keff)) then
      
         call petsc_solve(np_flux, &
                          matrix, &
                          rhs, &
                          state, &
                          option_path = trim(np_flux%option_path) )
      
      else solve
      
         call petsc_solve(delta_np_flux, &
                          matrix, &
                          rhs, &
                          state, &
                          option_path = trim(np_flux%option_path) )
         
         call addto(np_flux, &
                    delta_np_flux)
      
      end if solve
      
      ! set the direchlet bc nodes to be consistent
      call set_dirichlet_consistent(np_flux)
      
      ! deallocate the matrix and rhs vector used for solving
      call deallocate(matrix)
      call deallocate(rhs)
      
      ! deallocate the time run solution field passed to the solver
      dealloc_delta: if (.not. present(keff)) then
         
         call deallocate(delta_np_flux)
         
      end if dealloc_delta
      
      ! deallocate the prescribed source options data
      if (allocated(np_group_g_assemble_options%prescribed_source_energy_group)) deallocate(np_group_g_assemble_options%prescribed_source_energy_group)
      if (allocated(np_group_g_assemble_options%prescribed_source_field_name))   deallocate(np_group_g_assemble_options%prescribed_source_field_name)
      
   end subroutine np_assemble_solve_group_set

   ! --------------------------------------------------------------------------

   subroutine get_group_g_assemble_options(np_group_g_assemble_options, &
                                           np_flux_option_path, &
                                           np_radmat_option_path, &
                                           time_run)
      
      !!< Get the options associated with the assembly of the np group p
      
      type(np_group_g_assemble_options_type), intent(out) :: np_group_g_assemble_options 
      character(len=*) :: np_flux_option_path
      character(len=*) :: np_radmat_option_path
      logical :: time_run
      
      ! local variables
      integer :: s
      character(len=OPTION_PATH_LEN) :: include_prescribed_source_path
            
      ! get the spatial discretisation
      spatial: if (have_option(trim(np_flux_option_path)//'/prognostic/spatial_discretisation/continuous_galerkin')) then
      
         np_group_g_assemble_options%spatial_discretisation = 'continuous_galerkin'
               
      else spatial
      
         FLAbort('Unknown spatial discretisation for radiation field')
      
      end if spatial
      
      ! get the time theta and time step
      time: if (time_run) then
            
         call get_option(trim(np_radmat_option_path)//'/time_run/temporal_discretisation/theta',np_group_g_assemble_options%theta)
         
         call get_option('/timestepping/timestep',np_group_g_assemble_options%timestep)
   
         ! prescribed source data
         include_prescribed_source_path = trim(np_radmat_option_path)//'/time_run/include_prescribed_source'
         
         np_group_g_assemble_options%include_prescribed_source = have_option(trim(include_prescribed_source_path))
         
         get_prescribed_source_options: if (np_group_g_assemble_options%include_prescribed_source) then 
         
            np_group_g_assemble_options%number_prescribed_source =  &
            option_count(trim(include_prescribed_source_path)//'/source_scalar_field')
         
            allocate(np_group_g_assemble_options%prescribed_source_energy_group(np_group_g_assemble_options%number_prescribed_source))
         
            allocate(np_group_g_assemble_options%prescribed_source_field_name(np_group_g_assemble_options%number_prescribed_source))
         
            prescribed_source_loop: do s = 1,np_group_g_assemble_options%number_prescribed_source
            
               call get_option(trim(include_prescribed_source_path)//'/source_scalar_field['//int2str(s - 1)//']/name', &
                               np_group_g_assemble_options%prescribed_source_field_name(s)) 
            
               all_group: if (have_option(trim(include_prescribed_source_path)//'/source_scalar_field['//int2str(s - 1)//']/all_group')) then
               
                  np_group_g_assemble_options%prescribed_source_energy_group(s) = 0
            
               else all_group
            
                  call get_option(trim(include_prescribed_source_path)//'/source_scalar_field['//int2str(s - 1)//']/energy_group', &
                                  np_group_g_assemble_options%prescribed_source_energy_group(s))           
            
               end if all_group
            
            end do prescribed_source_loop
         
         end if get_prescribed_source_options
              
      end if time
      
      delayed: if (have_option(trim(np_radmat_option_path)//'/delayed_neutron_precursor')) then
      
         np_group_g_assemble_options%have_delayed = .true.
      
      else delayed
      
         np_group_g_assemble_options%have_delayed = .false.
      
      end if delayed
                  
   end subroutine get_group_g_assemble_options

   ! --------------------------------------------------------------------------

   subroutine assemble_np_group_g_vele(vele, &
                                       g, &
                                       number_of_energy_groups, &
                                       state, &
                                       np_radmat_name, &
                                       np_flux, &
                                       np_flux_old, &
                                       radmat_vele, &
                                       matrix, &
                                       rhs, &
                                       positions, &
                                       np_group_g_assemble_options, &
                                       keff)

      !!< Assemble the volume element vele contribution to the global matrix and rhs for this np for group g
      !!< The optional argument keff is used to determine whether this is a eigenvalue or time dependent assemble
      
      ! the assemble for time run is still under development and is currently known to be incorrect 
      
      integer, intent(in) :: vele
      integer, intent(in) :: g
      integer, intent(in) :: number_of_energy_groups 
      type(state_type), intent(in) :: state    
      character(len=*), intent(in) :: np_radmat_name
      type(scalar_field), intent(inout) :: np_flux 
      type(scalar_field), intent(in) :: np_flux_old
      type(radmat_type), intent(in) :: radmat_vele      
      type(csr_matrix), intent(inout) :: matrix
      type(scalar_field), intent(inout) :: rhs
      type(vector_field), intent(in) :: positions
      type(np_group_g_assemble_options_type), intent(in) :: np_group_g_assemble_options
      real, intent(in), optional :: keff
      
      ! local variables 
      integer :: status  
      integer :: s
      integer :: g_dash   
      real :: timestep
      real :: theta
      real :: timestep_theta
      real, dimension(ele_ngi(np_flux, vele)) :: detwei
      real, dimension(ele_loc(np_flux, vele), ele_ngi(np_flux, vele), mesh_dim(np_flux)) :: dshape
      real, dimension(ele_loc(np_flux, vele), ele_loc(np_flux, vele)) :: mass_matrix     
      real, dimension(ele_loc(np_flux, vele), ele_loc(np_flux, vele)) :: diffusivity_matrix
      real, dimension(positions%dim, positions%dim, ele_ngi(np_flux, vele)) :: diffusivity_gi                 
      real, dimension(ele_loc(np_flux, vele), ele_loc(np_flux, vele)) :: matrix_addto
      real, dimension(ele_loc(np_flux, vele)) :: rhs_addto      
      real, dimension(:), allocatable :: np_flux_vele_val
      real, dimension(:), allocatable :: np_flux_old_vele_val
      type(scalar_field), pointer :: prescribed_isotropic_source_field 
      
      ! form the timestep*theta use for time run and shorter variable names
      form_timestep_theta: if (.not. present(keff)) then
         
         timestep_theta = np_group_g_assemble_options%timestep* &
                          np_group_g_assemble_options%theta
         
         timestep = np_group_g_assemble_options%timestep
         
         theta = np_group_g_assemble_options%theta
         
      end if form_timestep_theta
      
      ! initialise the local vele matrix and vector that are added to the global arrays
      matrix_addto = 0.0
      rhs_addto    = 0.0
            
      ! form the velement jacobian transform and gauss weight
      call transform_to_physical(positions, vele, ele_shape(np_flux, vele), dshape = dshape, detwei = detwei)
      
      ! form the vele mass matrix
      mass_matrix = shape_shape(ele_shape(np_flux,vele),ele_shape(np_flux,vele),detwei)
            
      ! form the vele diffusion tensor coeff which is diagonal and the same value for 
      ! each gauss point as volume element wise
      diffusivity_gi = 0.0
      diffusivity_gi(1,1,:) = radmat_vele%diffusion(g,1)
      if (positions%dim > 1) diffusivity_gi(2,2,:) = radmat_vele%diffusion(g,2)
      if (positions%dim > 2) diffusivity_gi(3,3,:) = radmat_vele%diffusion(g,3)
      
      ! form the vele diffusion matrix
      diffusivity_matrix = dshape_tensor_dshape(dshape, diffusivity_gi, dshape, detwei)
         
      keff_or_time_diff: if (present(keff)) then
         
         matrix_addto = matrix_addto + diffusivity_matrix
         
      else keff_or_time_diff
         
         matrix_addto = matrix_addto + diffusivity_matrix
         
      end if keff_or_time_diff
               
      ! add the vele removal term 
      keff_or_time_removal: if (present(keff)) then
         
         matrix_addto = matrix_addto + mass_matrix*radmat_vele%removal(g,1)
         
      else keff_or_time_removal
         
         matrix_addto = matrix_addto + mass_matrix*radmat_vele%removal(g,1)
         
      end if keff_or_time_removal 
                                
      ! include the isotropic sources - scatter, fission, delayed and prescribed
         
      ! include the other group to this group g scatter source and the spectrum production source         
      group_loop: do g_dash = 1,number_of_energy_groups
               
         ! get the vele scalar flux values of group g_dash
         call extract_flux_group_g_vele_val(state, &
                                            trim(np_radmat_name), & 
                                            g_dash, &
                                            vele, &
                                            np_flux_vele_val = np_flux_vele_val, &
                                            np_flux_old_vele_val = np_flux_old_vele_val) 
                  
         keff_or_time: if (present(keff)) then
            
            ! add the scatter
            not_within_group_eig: if (g /= g_dash) then
                  
               rhs_addto = rhs_addto + matmul(mass_matrix,np_flux_vele_val)* &
                                       radmat_vele%scatter(g_dash,g,1)
          
            end if not_within_group_eig
            
            ! add the spectrum production - this is the eigenvector
            rhs_addto = rhs_addto + matmul(mass_matrix,np_flux_old_vele_val)* &
                                    radmat_vele%production(g_dash)*radmat_vele%prompt_spectrum(g)/keff
                                          
         else keff_or_time
  
            ! add the scatter
            not_within_group_time: if (g /= g_dash) then
                  
               rhs_addto = rhs_addto + matmul(mass_matrix,(theta*np_flux_vele_val + (1.0-theta)*np_flux_old_vele_val))* &
                                       radmat_vele%scatter(g_dash,g,1)
            end if not_within_group_time

            ! add the spectrum production    
            rhs_addto = rhs_addto + matmul(mass_matrix,(theta*np_flux_vele_val + (1.0-theta)*np_flux_old_vele_val))* &
                                    radmat_vele%production(g_dash)*radmat_vele%prompt_spectrum(g)
                           
         end if keff_or_time
                                 
      end do group_loop

      if (allocated(np_flux_vele_val))     deallocate(np_flux_vele_val)
      if (allocated(np_flux_old_vele_val)) deallocate(np_flux_old_vele_val)
      
      ! if time run include the prescribed source and delayed source
      time_run: if (.not. present(keff)) then
            
         ! include the delayed source
         add_delayed_source: if (np_group_g_assemble_options%have_delayed) then
            
            ! fill in ...
            
         end if add_delayed_source
      
         add_prescribed_isotropic_source: if (np_group_g_assemble_options%include_prescribed_source) then
            
            ! loop each prescribed source and add in
            prescribed_isotropic_source_loop: do s = 1,np_group_g_assemble_options%number_prescribed_source
               
               ! include the source into this group g if needed
               source_this_group: if ((np_group_g_assemble_options%prescribed_source_energy_group(s) == g) .or. &
                                      (np_group_g_assemble_options%prescribed_source_energy_group(s) == 0)) then
                   
                   ! extract the prescribed isotropic source field
                   prescribed_isotropic_source_field => extract_scalar_field(state, &
                                                                             trim(np_group_g_assemble_options%prescribed_source_field_name(s)), &
                                                                             stat=status)
                   
                   rhs_addto = rhs_addto + shape_rhs(ele_shape(np_flux,vele),detwei*ele_val(prescribed_isotropic_source_field,vele))
                  
               end if source_this_group
               
            end do prescribed_isotropic_source_loop
            
         end if add_prescribed_isotropic_source
      
      end if time_run
       
      ! insert the local vele matrix and rhs into the global arrays 
      call addto(matrix, &
                 ele_nodes(np_flux, vele), &
                 ele_nodes(np_flux, vele), &
                 matrix_addto)

      call addto(rhs, &
                 ele_nodes(np_flux, vele), &
                 rhs_addto)
              
   end subroutine assemble_np_group_g_vele

   ! --------------------------------------------------------------------------

   subroutine assemble_np_group_g_sele_vacuum(face, &
                                              positions, &
                                              matrix, & 
                                              bc_value, &
                                              np_flux)
      
      !!< Assemble the local surface element sele (now face) boundary condition contribution for neutral particle group g and 
      !!< add to the local matrix. This is currently only set up to do diffusion theory vacuum which involves adding 
      !!< a surface mass matrix * 0.5 to the global matrix (the 0.5 is set via spud options for vacuum)
      
      integer, intent(in) :: face
      type(vector_field), intent(in) :: positions
      type(csr_matrix), intent(inout) :: matrix
      type(scalar_field), intent(in) :: bc_value
      type(scalar_field), intent(in) :: np_flux 
      
      ! local variables
      integer, dimension(face_loc(np_flux, face)) :: face_nodes
      real, dimension(face_ngi(np_flux, face)) :: detwei    
      real, dimension(face_loc(np_flux, face), face_loc(np_flux, face)) :: matrix_addto
      
      assert(face_ngi(positions, face) == face_ngi(np_flux, face))
            
      call transform_facet_to_physical(positions, &
                                       face, &
                                       detwei_f = detwei)  
      
      matrix_addto = shape_shape(face_shape(np_flux, face), &
                                 face_shape(np_flux, face), &
                                 detwei*face_val_at_quad(bc_value,face))
      
      call addto(matrix, &
                 face_nodes, &
                 face_global_nodes(np_flux, face), &
                 matrix_addto)
                
   end subroutine assemble_np_group_g_sele_vacuum

   ! --------------------------------------------------------------------------

end module radiation_assemble_solve_group_set
