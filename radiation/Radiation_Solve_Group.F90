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

module radiation_assemble_solve_group

   !!< This module contains procedures associated with solving the 
   !!< group g particle matrix problem

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
   use radiation_energy_group_set_tools
     
   implicit none
   
   private 

   public :: particle_assemble_solve_group
   
   type particle_group_g_assemble_options_type
      character(len=OPTION_PATH_LEN) :: spatial_discretisation
      real :: theta
      real :: timestep
      logical :: have_delayed
      logical :: include_prescribed_source
   end type particle_group_g_assemble_options_type

contains

   ! --------------------------------------------------------------------------

   subroutine particle_assemble_solve_group(state, &
                                            particle_radmat, &
                                            particle_radmat_ii, &
                                            g, &
                                            number_of_energy_groups, &
                                            keff) 
   
      !!< Assemble and solve the group g particle matrix problem
      !!< If keff variable is present then this is an eigenvalue solve, else this is a time run solve

      type(state_type), intent(inout) :: state
      type(particle_radmat_type), intent(in) :: particle_radmat
      type(particle_radmat_ii_type), intent(in) :: particle_radmat_ii      
      integer, intent(in) :: g
      integer, intent(in) :: number_of_energy_groups 
      real, intent(in), optional :: keff
      
      ! local variables
      type(scalar_field_pointer), dimension(:), pointer :: particle_flux 
      type(scalar_field_pointer), dimension(:), pointer :: particle_flux_old      
      type(particle_group_g_assemble_options_type) :: particle_group_g_assemble_options
                 
      ! extract the particle flux fields for all energy groups
      call extract_flux_all_group(state, &
                                  trim(particle_radmat%name), & 
                                  number_of_energy_groups, &
                                  particle_flux     = particle_flux, &
                                  particle_flux_old = particle_flux_old)
      
      ! get the solver options for group g
      call get_group_g_assemble_options(particle_group_g_assemble_options, &
                                        particle_flux(g)%ptr%option_path, &
                                        particle_radmat%option_path, &
                                        (.not. present(keff)))
      
      ! currently assume all energy groups point to the same coordinate mesh
      
      ! if the above assumption was to be broken then the there would need to be a sweep of the other groups
      ! fluxes (new and old) to see which needs supermeshing to a new field and then the pointer within the 
      ! scalar_field_pointer particle_flux (or _old) is then changed such that no code change below is needed
      
      ! assemble the diffusivity, absorption and already discretised rhs term scalar fields and insert into state
      call assemble_coeff_source_group_g(particle_flux, &
                                         particle_flux_old, &
                                         particle_radmat, &
                                         particle_radmat_ii, &
                                         state, &
                                         g, &
                                         number_of_energy_groups, &
                                         keff = keff)
      
      call assemble_matrix_solve_group_g(particle_flux, &
                                         particle_flux_old, &
                                         state, &   
                                         g, &
                                         particle_group_g_assemble_options, &
                                         keff = keff)
      
      ! deallocate the particle flux pointer fields
      call deallocate_flux_all_group(particle_flux     = particle_flux, &
                                     particle_flux_old = particle_flux_old)
            
   end subroutine particle_assemble_solve_group

   ! --------------------------------------------------------------------------

   subroutine get_group_g_assemble_options(particle_group_g_assemble_options, &
                                           particle_flux_option_path, &
                                           particle_radmat_option_path, &
                                           time_run)
      
      !!< Get the options associated with the assembly of the particle group g
      
      type(particle_group_g_assemble_options_type), intent(out) :: particle_group_g_assemble_options 
      character(len=*) :: particle_flux_option_path
      character(len=*) :: particle_radmat_option_path
      logical :: time_run
            
      ! get the spatial discretisation
      spatial: if (have_option(trim(particle_flux_option_path)//'/prognostic/spatial_discretisation/continuous_galerkin')) then
      
         particle_group_g_assemble_options%spatial_discretisation = 'continuous_galerkin'
               
      else spatial
      
         FLAbort('Unknown spatial discretisation for radiation field')
      
      end if spatial
      
      ! get the time theta and time step
      time: if (time_run) then
            
         call get_option(trim(particle_flux_option_path)//'/prognostic/temporal_discretisation/theta',particle_group_g_assemble_options%theta)
         
         call get_option('/timestepping/timestep',particle_group_g_assemble_options%timestep)
         
         particle_group_g_assemble_options%include_prescribed_source = &
         have_option(trim(particle_flux_option_path)//'/prognostic/scalar_field::Source')
                       
      end if time
      
      delayed: if (have_option(trim(particle_radmat_option_path)//'/delayed_neutron_precursor')) then
      
         particle_group_g_assemble_options%have_delayed = .true.
      
      else delayed
      
         particle_group_g_assemble_options%have_delayed = .false.
      
      end if delayed
                  
   end subroutine get_group_g_assemble_options

   ! --------------------------------------------------------------------------
   
   subroutine assemble_coeff_source_group_g(particle_flux, &
                                            particle_flux_old, &
                                            particle_radmat, &
                                            particle_radmat_ii, &
                                            state, &
                                            g, &
                                            number_of_energy_groups, &
                                            keff)
      
      !!< Assemble the coeff and discretised source fields for group g that will then be used 
      !!< somewhere else to assemble the linear system

      type(scalar_field_pointer), dimension(:), intent(inout) :: particle_flux 
      type(scalar_field_pointer), dimension(:), intent(in) :: particle_flux_old      
      type(particle_radmat_type), intent(in) :: particle_radmat
      type(particle_radmat_ii_type), intent(in) :: particle_radmat_ii  
      type(state_type), intent(inout) :: state    
      integer, intent(in) :: g
      integer, intent(in) :: number_of_energy_groups 
      real, intent(in), optional :: keff
      
      ! local variables
      integer :: stat
      integer :: vele
      integer :: inode
      integer :: g_dash
      integer :: g_set
      real :: data_value
      real, dimension(:), allocatable :: detwei_vele
      real, dimension(:), allocatable :: rhs_addto      
      type(scalar_field) :: production_coeff 
      type(scalar_field) :: prompt_spectrum_coeff 
      type(scalar_field) :: scatter_coeff            
      type(scalar_field), pointer :: absorption_coeff
      type(tensor_field), pointer :: diffusivity_coeff
      type(scalar_field), pointer :: discretised_source
      type(mesh_type), pointer :: material_fn_space
      type(vector_field), pointer :: positions      
      character(len=OPTION_PATH_LEN) :: material_fn_space_name
      character(len=OPTION_PATH_LEN) :: energy_group_set_path
      character(len=OPTION_PATH_LEN) :: field_name
      
      ! determine which group set this g belongs to
      call which_group_set_contains_g(g, &
                                      trim(particle_radmat%option_path), &
                                      g_set)

      ! set the energy_group_set path
      energy_group_set_path = trim(particle_radmat%option_path)//'/energy_discretisation/energy_group_set['//int2str(g_set - 1)//']'
         
      ! get the material fn space name for this group set
      call get_option(trim(energy_group_set_path)//'/angular_discretisation/method/mesh/name',material_fn_space_name)
      
      ! extract the material fn_space of this energy group set of this particle type 
      material_fn_space => extract_mesh(state, trim(material_fn_space_name))
      
      ! get the positions field for this energy group set
      positions => extract_vector_field(state, "Coordinate")

      ! allocate the local fields as needed

      field_name = trim(particle_flux(g)%ptr%name)//'Scatter'
      call allocate(scatter_coeff, material_fn_space, trim(field_name))
      call zero(scatter_coeff)

      field_name = trim(particle_flux(g)%ptr%name)//'Production'
      call allocate(production_coeff, material_fn_space, trim(field_name))
      call zero(production_coeff)

      field_name = trim(particle_flux(g)%ptr%name)//'PromptSpectrum'
      call allocate(prompt_spectrum_coeff, material_fn_space, trim(field_name))
      call zero(prompt_spectrum_coeff)
      
      ! extract the assemble fields as needed
      absorption_coeff => extract_scalar_field(state, &
                                               trim(particle_flux(g)%ptr%name) // 'Absorption', &
                                               stat = stat)
                                        
      diffusivity_coeff => extract_tensor_field(state, &
                                                trim(particle_flux(g)%ptr%name) // 'Diffusivity', &
                                                stat = stat)

      discretised_source => extract_scalar_field(state, &
                                                 trim(particle_flux(g)%ptr%name) // 'DiscretisedSource', &
                                                 stat = stat)
      
      ! zero the extracted fields
      call zero(absorption_coeff)
      call zero(diffusivity_coeff)
      call zero(discretised_source)
            
      ! form the absorption field for this energy group (which is actually the removal cross section)            
      call form(material_fn_space, &
                particle_radmat_ii%energy_group_set_ii(g_set), &
                particle_radmat, &
                absorption_coeff, &
                g, &
                component = 'removal')
            
      ! form the diffusivity tensor field for this energy group (only fill in diagonals)            
      node_loop_d: do inode = 1,node_count(material_fn_space)                                                                      
            
         ! get the inode diffusionx data for this energy group
         call form(particle_radmat_ii%energy_group_set_ii(g_set), &
                   particle_radmat, &
                   inode, &
                   g, &
                   data_value, &
                   component = 'diffusionx')
               
         ! set the interpolated value into the scalar field
         call set(diffusivity_coeff, &
                  1, &
                  1, &
                  inode, &
                  data_value)

         if (positions%dim > 1) then
         
            ! get the inode diffusiony data for this energy group
            call form(particle_radmat_ii%energy_group_set_ii(g_set), &
                      particle_radmat, &
                      inode, &
                      g, &
                      data_value, &
                      component = 'diffusiony')
               
            ! set the interpolated value into the scalar field
            call set(diffusivity_coeff, &
                     2, &
                     2, &
                     inode, &
                     data_value)

         end if 
         
         if (positions%dim > 2) then

            ! get the inode diffusionz data for this energy group
            call form(particle_radmat_ii%energy_group_set_ii(g_set), &
                      particle_radmat, &
                      inode, &
                      g, &
                      data_value, &
                      component = 'diffusionz')
               
            ! set the interpolated value into the scalar field
            call set(diffusivity_coeff, &
                     3, &
                     3, &
                     inode, &
                     data_value)

         end if 
               
      end do node_loop_d
            
      ! form the discretised source - scatter and production            

      call form(material_fn_space, &
                particle_radmat_ii%energy_group_set_ii(g_set), &
                particle_radmat, &
                prompt_spectrum_coeff, &
                g, &
                component = 'prompt_spectrum')
                  
      group_loop: do g_dash = 1,number_of_energy_groups
                        
         ! form the production field for this energy group g_dash
         
         call form(material_fn_space, &
                   particle_radmat_ii%energy_group_set_ii(g_set), &
                   particle_radmat, &
                   production_coeff, &
                   g_dash, &
                   component = 'production')
         
         not_within_group: if (g /= g_dash) then
 
            ! form the scatter field for this energy group g_dash to g            
            call form(material_fn_space, &
                      particle_radmat_ii%energy_group_set_ii(g_set), &
                      particle_radmat, &
                      scatter_coeff, &
                      g, &
                      component = 'scatter', &
                      g_dash = g_dash)
                                           
         end if not_within_group
         
         vele_loop: do vele = 1,ele_count(discretised_source)
           
            ! allocate the jacobian transform and gauss weight array for this vele
            allocate(detwei_vele(ele_ngi(discretised_source,vele)))
                        
            ! form the velement jacobian transform and gauss weight
            call transform_to_physical(positions, vele, detwei = detwei_vele)
            
            ! allocate the rhs_addto
            allocate(rhs_addto(ele_loc(particle_flux(g)%ptr,vele)))
            
            rhs_addto = 0.0
            
            keff_or_time: if (present(keff)) then
            
               ! add the scatter - not the within group
               not_within_group_eig: if (g /= g_dash) then
                  
                  rhs_addto = rhs_addto + matmul(shape_shape(ele_shape(particle_flux(g)%ptr,vele), &
                                                             ele_shape(particle_flux(g_dash)%ptr,vele), &
                                                             detwei_vele*ele_val_at_quad(scatter_coeff,vele)), &
                                                 ele_val(particle_flux(g_dash)%ptr,vele))
          
               end if not_within_group_eig
            
               ! add the spectrum production - this is the eigenvector
               rhs_addto = rhs_addto + matmul(shape_shape(ele_shape(particle_flux(g)%ptr,vele), &
                                                          ele_shape(particle_flux_old(g_dash)%ptr,vele), &
                                                          detwei_vele*ele_val_at_quad(production_coeff,vele) &
                                                                     *ele_val_at_quad(prompt_spectrum_coeff,vele) &
                                                                     *(1.0/keff) ), &
                                              ele_val(particle_flux_old(g_dash)%ptr,vele))
                                          
            else keff_or_time
  
               ! add the scatter - not the within group
               not_within_group_time: if (g /= g_dash) then
                  

               end if not_within_group_time

               ! add the spectrum production 
                  
                           
            end if keff_or_time
            
            ! add contribution to discretised source for this vele for this group g
            call addto(discretised_source, &
                       ele_nodes(discretised_source, vele), &
                       rhs_addto)
            
            deallocate(detwei_vele)
            deallocate(rhs_addto)
            
         end do vele_loop
               
      end do group_loop
      
      ! deallocate the local fields as needed
            
      call deallocate(scatter_coeff)
      call deallocate(production_coeff)
      call deallocate(prompt_spectrum_coeff)
                        
   end subroutine assemble_coeff_source_group_g
   
   ! --------------------------------------------------------------------------
   
   subroutine assemble_matrix_solve_group_g(particle_flux, &
                                            particle_flux_old, &
                                            state, &
                                            g, &
                                            particle_group_g_assemble_options, &
                                            keff)
      
      !!< Assemble the matrix system for this group g and solve
      
      type(scalar_field_pointer), dimension(:), intent(inout) :: particle_flux 
      type(scalar_field_pointer), dimension(:), intent(in) :: particle_flux_old       
      type(state_type), intent(inout) :: state   
      integer, intent(in) :: g 
      type(particle_group_g_assemble_options_type), intent(in) :: particle_group_g_assemble_options
      real, intent(in), optional :: keff
      
      ! local variables
      integer :: stat
      integer :: vele,sele
      type(csr_matrix) :: matrix
      type(scalar_field) :: rhs
      type(csr_sparsity), pointer :: sparsity
      type(scalar_field) :: delta_particle_flux
      integer, dimension(:), allocatable :: bc_types
      type(scalar_field) :: bc_value
      type(tensor_field), pointer :: diffusivity
      type(scalar_field), pointer :: absorption
      type(scalar_field), pointer :: discretised_source
      type(vector_field), pointer :: positions      
            
      ! allocate the solution field passed to the solve for a time run 
      alloc_delta: if (.not. present(keff)) then
         
         call allocate(delta_particle_flux,&
                       particle_flux(g)%ptr%mesh, &
                       name = 'DeltaParticleFluxGroup')
         
         ! set the time run initial guess
         call zero(delta_particle_flux)
         
      end if alloc_delta
      
      ! determine the sparsity pattern of matrix assuming first order connections 
      sparsity => get_csr_sparsity_firstorder(state, particle_flux(g)%ptr%mesh, particle_flux(g)%ptr%mesh)
      
      ! allocate the matrix and rhs vector used for solving
      call allocate(matrix, &
                    sparsity, &
                    name = 'MatrixParticleFluxGroup')
                    
      call allocate(rhs, &
                    particle_flux(g)%ptr%mesh, &
                    name = 'RHSParticleFluxGroup')
      
      call zero(matrix)
      call zero(rhs)

      ! get the positions field for this energy group set
      positions => extract_vector_field(state, 'Coordinate')
      
      ! extract the absorption, diffusivity and discretisedsource material fields
      absorption => extract_scalar_field(state, &
                                         trim(particle_flux(g)%ptr%name) // 'Absorption', &
                                         stat = stat)
                                        
      diffusivity => extract_tensor_field(state, &
                                          trim(particle_flux(g)%ptr%name) // 'Diffusivity', &
                                          stat = stat)

      discretised_source => extract_scalar_field(state, &
                                                 trim(particle_flux(g)%ptr%name) // 'DiscretisedSource', &
                                                 stat = stat)
  
      ! volume integrations 
      volume_element_loop: do vele = 1, ele_count(particle_flux(g)%ptr)      
                  
         call assemble_particle_group_g_vele(vele, &
                                             g, &
                                             particle_flux, &
                                             particle_flux_old, &
                                             matrix, &
                                             rhs, &
                                             positions, &
                                             particle_group_g_assemble_options, &
                                             absorption, &
                                             diffusivity, &
                                             state, &
                                             keff = keff)

      end do volume_element_loop
      
      ! include the already discretised_source into rhs
      call addto(rhs, discretised_source)
      
      ! surface integrations for albedo and source BC's for diffusion theory
      
      allocate(bc_types(surface_element_count(particle_flux(g)%ptr)))

      call get_entire_boundary_condition(particle_flux(g)%ptr, &
                                         (/"albedo", &
                                           "source"/), &
                                         bc_value, &
                                         bc_types)

      surface_element_loop: do sele = 1,surface_element_count(particle_flux(g)%ptr)
         
         diffusion_bc: if (bc_types(sele) == 1) then
         
            call assemble_particle_group_g_sele_albedo(sele, &
                                                       positions, &
                                                       matrix, & 
                                                       bc_value, &
                                                       particle_flux(g)%ptr)
         
         else if (bc_types(sele) == 2) then 
            
            call assemble_particle_group_g_sele_source(sele, &
                                                       positions, &
                                                       rhs, & 
                                                       bc_value, &
                                                       particle_flux(g)%ptr)
                        
         end if diffusion_bc
         
      end do surface_element_loop

      call deallocate(bc_value)

      if (allocated(bc_types)) deallocate(bc_types)
            
      ! apply direchlet BC
      call apply_dirichlet_conditions(matrix, &
                                      rhs, &
                                      particle_flux(g)%ptr)
            
      ! solve the linear system 
      ! - for eig solve use the latest solution as initial guess
      ! - for time solve update the state solution arrays after
      solve: if (present(keff)) then
      
         call petsc_solve(particle_flux(g)%ptr, &
                          matrix, &
                          rhs, &
                          state, &
                          option_path = trim(particle_flux(g)%ptr%option_path) )
      
      else solve
      
         call petsc_solve(delta_particle_flux, &
                          matrix, &
                          rhs, &
                          state, &
                          option_path = trim(particle_flux(g)%ptr%option_path) )
         
         call addto(particle_flux(g)%ptr, &
                    delta_particle_flux)
      
      end if solve
      
      ! set the direchlet bc nodes to be consistent
      call set_dirichlet_consistent(particle_flux(g)%ptr)
   
      ! deallocate the matrix and rhs vector used for solving
      call deallocate(matrix)
      call deallocate(rhs)
      
      ! deallocate the time run solution field passed to the solver
      dealloc_delta: if (.not. present(keff)) then
         
         call deallocate(delta_particle_flux)
         
      end if dealloc_delta    
              
   end subroutine assemble_matrix_solve_group_g

   ! --------------------------------------------------------------------------

   subroutine assemble_particle_group_g_vele(vele, &
                                             g, &
                                             particle_flux, &
                                             particle_flux_old, &
                                             matrix, &
                                             rhs, &
                                             positions, &
                                             particle_group_g_assemble_options, &
                                             absorption, &
                                             diffusivity, &
                                             state, &
                                             keff)

      !!< Assemble the volume element vele contribution to the global matrix and rhs for this particle type for group g
      !!< The optional argument keff is used to determine whether this is a eigenvalue or time dependent assemble
      
      ! the particle_flux (and old) fields have the same positions mesh here (so they may actually be the supermesh'ed flux)
      
      ! the assemble for time run is still under development and is currently known to be incorrect 
      
      integer, intent(in) :: vele
      integer, intent(in) :: g
      type(scalar_field_pointer), dimension(:), intent(in) :: particle_flux 
      type(scalar_field_pointer), dimension(:), intent(in) :: particle_flux_old       
      type(csr_matrix), intent(inout) :: matrix
      type(scalar_field), intent(inout) :: rhs
      type(vector_field), intent(in) :: positions
      type(particle_group_g_assemble_options_type), intent(in) :: particle_group_g_assemble_options
      type(scalar_field), intent(in) :: absorption
      type(tensor_field), intent(in) :: diffusivity
      type(state_type), intent(in) :: state    
      real, intent(in), optional :: keff
      
      ! local variables 
      integer :: status  
      real :: timestep
      real :: theta
      real :: timestep_theta
      real, dimension(ele_ngi(particle_flux(g)%ptr, vele)) :: detwei
      real, dimension(ele_loc(particle_flux(g)%ptr, vele), ele_ngi(particle_flux(g)%ptr, vele), mesh_dim(particle_flux(g)%ptr)) :: dshape
      real, dimension(positions%dim, positions%dim, ele_ngi(particle_flux(g)%ptr, vele)) :: diffusivity_gi                 
      real, dimension(ele_loc(particle_flux(g)%ptr, vele), ele_loc(particle_flux(g)%ptr, vele)) :: matrix_addto
      real, dimension(ele_loc(particle_flux(g)%ptr, vele)) :: rhs_addto      
      type(scalar_field), pointer :: prescribed_isotropic_source_field 
      
      ! form the timestep*theta use for time run and shorter variable names
      form_timestep_theta: if (.not. present(keff)) then
         
         timestep_theta = particle_group_g_assemble_options%timestep* &
                          particle_group_g_assemble_options%theta
         
         timestep = particle_group_g_assemble_options%timestep
         
         theta = particle_group_g_assemble_options%theta
         
      end if form_timestep_theta
      
      ! initialise the local vele matrix and vector that are added to the global arrays
      matrix_addto = 0.0
      rhs_addto    = 0.0
           
      ! form the velement jacobian transform and gauss weight
      call transform_to_physical(positions, &
                                 vele, &
                                 ele_shape(particle_flux(g)%ptr, vele), &
                                 dshape = dshape, &
                                 detwei = detwei)
            
      ! add the diffusion term
      diffusivity_gi = ele_val_at_quad(diffusivity, vele)
            
      keff_or_time_diff: if (present(keff)) then
         
         matrix_addto = matrix_addto + dshape_tensor_dshape(dshape, &
                                                            diffusivity_gi, &
                                                            dshape, &
                                                            detwei)

      else keff_or_time_diff
         
            ! fill in ...
         
      end if keff_or_time_diff
               
      ! add the absortion term 
      keff_or_time_removal: if (present(keff)) then
         
         matrix_addto = matrix_addto + shape_shape(ele_shape(particle_flux(g)%ptr,vele), &
                                                   ele_shape(particle_flux(g)%ptr,vele), &
                                                   detwei*ele_val_at_quad(absorption,vele))

      else keff_or_time_removal
         
            ! fill in ...
         
      end if keff_or_time_removal 
                                               
      ! if time run include the prescribed source
      time_run: if (.not. present(keff)) then
                  
         add_prescribed_isotropic_source: if (particle_group_g_assemble_options%include_prescribed_source) then
                   
             ! extract the prescribed isotropic source field
             prescribed_isotropic_source_field => extract_scalar_field(state, &
                                                                       trim(particle_flux(g)%ptr%name)//'Source', &
                                                                       stat=status)
                   
             rhs_addto = rhs_addto + shape_rhs(ele_shape(particle_flux(g)%ptr,vele), &
                                               detwei*ele_val(prescribed_isotropic_source_field, &
                                               vele))
                                             
         end if add_prescribed_isotropic_source
         
         add_delayed: if (particle_group_g_assemble_options%have_delayed) then
         
            ! fill in ...
            
         end if add_delayed
      
      end if time_run

      ! insert the local vele matrix and rhs into the global arrays 
      call addto(matrix, &
                 ele_nodes(particle_flux(g)%ptr, vele), &
                 ele_nodes(particle_flux(g)%ptr, vele), &
                 matrix_addto)

      call addto(rhs, &
                 ele_nodes(particle_flux(g)%ptr, vele), &
                 rhs_addto)
              
   end subroutine assemble_particle_group_g_vele

   ! --------------------------------------------------------------------------

   subroutine assemble_particle_group_g_sele_albedo(face, &
                                                    positions, &
                                                    matrix, & 
                                                    bc_value, &
                                                    particle_flux)
      
      !!< Assemble the diffusion theory albedo boundary condition that involves adding 
      !!< a surface mass matrix * coefficient to the main matrix
      
      ! this is actually a robin BC for angular diffusion
      
      integer, intent(in) :: face
      type(vector_field), intent(in) :: positions
      type(csr_matrix), intent(inout) :: matrix
      type(scalar_field), intent(in) :: bc_value
      type(scalar_field), pointer :: particle_flux 
      
      ! local variables
      real, dimension(face_ngi(particle_flux, face)) :: detwei    
      real, dimension(face_loc(particle_flux, face), face_loc(particle_flux, face)) :: matrix_addto
      
      assert(face_ngi(positions, face) == face_ngi(particle_flux, face))
            
      call transform_facet_to_physical(positions, &
                                       face, &
                                       detwei_f = detwei)  
      
      ! the 0.5 is a necessary part of this BC formulation, the bc_value is the albedo coeff
      ! where a value of 1.0 is a perfect reflective and a value of 0.0 is a perfect vacuum
      matrix_addto = shape_shape(face_shape(particle_flux, face), &
                                 face_shape(particle_flux, face), &
                                 detwei*0.5*( (1.0 - ele_val_at_quad(bc_value,face) ) / &
                                              (1.0 + ele_val_at_quad(bc_value,face) ) ) )
      
      call addto(matrix, &
                 face_global_nodes(particle_flux,face), &
                 face_global_nodes(particle_flux,face), &
                 matrix_addto)
                
   end subroutine assemble_particle_group_g_sele_albedo

   ! --------------------------------------------------------------------------

   subroutine assemble_particle_group_g_sele_source(face, &
                                                    positions, &
                                                    rhs, & 
                                                    bc_value, &
                                                    particle_flux)
      
      !!< Include the particle source boundary condition into the rhs
      
      ! this is actually a neumann BC for angular diffusion
      
      integer, intent(in) :: face
      type(vector_field), intent(in) :: positions
      type(scalar_field), intent(inout) :: rhs
      type(scalar_field), intent(in) :: bc_value
      type(scalar_field), pointer :: particle_flux 
      
      ! local variables
      real, dimension(face_ngi(particle_flux, face)) :: detwei    
      real, dimension(ele_loc(particle_flux, face)) :: rhs_addto      
      
      assert(face_ngi(positions, face) == face_ngi(particle_flux, face))
            
      call transform_facet_to_physical(positions, &
                                       face, &
                                       detwei_f = detwei)  
      
      rhs_addto = shape_rhs(face_shape(particle_flux, face), &
                            detwei*ele_val_at_quad(bc_value,face) )
      
      call addto(rhs, &
                 face_global_nodes(particle_flux,face), &
                 rhs_addto)
                
   end subroutine assemble_particle_group_g_sele_source

   ! --------------------------------------------------------------------------

end module radiation_assemble_solve_group
