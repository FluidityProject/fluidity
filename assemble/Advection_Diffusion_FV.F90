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

module advection_diffusion_fv
  !!< This module contains the finite volume form of the advection
  !!< -diffusion equation for scalars.
  use fldebug
  use vector_tools
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN, COLOURING_DG2, &
COLOURING_DG0, COLOURING_DG1
  use elements
  use spud
  use integer_set_module
#ifdef _OPENMP
  use omp_lib
#endif
  use sparse_tools
  use shape_functions
  use transform_elements
  use fetools
  use fields
  use profiler
  use state_module
  use boundary_conditions
  use sparsity_patterns
  use dgtools
  use vtk_interfaces
  use field_options
  use sparse_matrices_fields
  use fefields
  use field_derivatives
  use coordinates
  use sparsity_patterns_meshes
  use petsc_solve_state_module
  use boundary_conditions_from_options
  use upwind_stabilisation
  use slope_limiters_dg
  use diagnostic_fields, only: calculate_diagnostic_variable
  use colouring
  
  implicit none

  private
  public solve_advection_diffusion_fv, advection_diffusion_fv_check_options

  ! Mass term?
  logical :: have_mass
  ! Advection?
  logical :: have_advection
  ! Source?
  logical :: have_source
  ! Add source directly to the right hand side?
  logical :: add_src_directly_to_rhs
  ! Absorption?
  logical :: have_absorption
  ! Diffusivity?
  logical :: have_diffusivity
  ! Isotropic diffusivity?
  logical :: isotropic_diffusivity
  ! Is the mesh moving?
  logical :: move_mesh
  
  ! timestepping parameters
  real :: theta, dt, dt_theta

contains

  subroutine solve_advection_diffusion_fv(field_name, state)
    !!< Construct and solve the advection-diffusion equation for the given
    !!< field unsing element centred finite volumes.
    
    !! Name of the field to be solved for.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state

    !! Tracer to be solved for.
    type(scalar_field), pointer :: T, T_old
    !! Change in T over one timestep.
    type(scalar_field) :: delta_T
    !! System matrix.
    type(csr_matrix) :: matrix
    !! Right hand side vector.
    type(scalar_field) :: rhs
    !! Sparsity of advection_diffusion matrix
    type(csr_sparsity), pointer :: sparsity
    
    ewrite(1,*) 'In solve_advection_diffusion_fv'

    t=>extract_scalar_field(state, field_name)
    if((continuity(T)>=0).or.(element_degree(T,1)/=0)) then
      FLExit("FV advection-diffusion requires a discontinuous P0 mesh.")
    end if
    
    t_old=>extract_scalar_field(state, "Old"//field_name)

    ! Reset T to value at the beginning of the timestep.
    call set(t, t_old)

    sparsity => get_csr_sparsity_firstorder(state, t%mesh, t%mesh)
    
    call allocate(matrix, sparsity, name = trim(field_name)//"Matrix")
    call allocate(rhs, t%mesh, name = trim(field_name)//"RHS")
    
    call allocate(delta_t, t%mesh, "Delta"//trim(field_name))
    call zero(delta_t)
    
    call get_option("/timestepping/timestep", dt)
    
    call assemble_advection_diffusion_fv(t, matrix, rhs, state)
    
    call petsc_solve(delta_t, matrix, rhs, state, option_path = trim(t%option_path))
    
    ewrite_minmax(delta_t)
    
    call addto(t, delta_t, dt)
    
    ewrite_minmax(t)
    
    call deallocate(matrix)
    call deallocate(rhs)
    call deallocate(delta_t)

    ewrite(1,*) 'Exiting solve_advection_diffusion_fv'

  end subroutine solve_advection_diffusion_fv
  
  subroutine assemble_advection_diffusion_fv(t, matrix, rhs, state)
  
    type(scalar_field), intent(inout) :: t
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    type(state_type), intent(inout) :: state
    
    type(vector_field), pointer :: coordinate, &
                                   old_coordinate, new_coordinate, &
                                   grid_velocity
    type(vector_field) :: t_coordinate                               
    type(scalar_field), pointer :: source, absorption
    type(tensor_field), pointer :: diffusivity
    
    integer :: i, j, stat

    !! Coloring  data structures for OpenMP parallization
    type(integer_set), dimension(:), pointer :: colours
    integer :: clr, nnid, len, ele
    integer :: thread_num
    !! Did we successfully prepopulate the transform_to_physical_cache?
    logical :: cache_valid
    
    ewrite(1,*) "In assemble_advection_diffusion_fv"
    
    coordinate => extract_vector_field(state, "Coordinate")
    assert(coordinate%dim == mesh_dim(t))
    assert(ele_count(coordinate) == ele_count(t))
    ! Source
    source => extract_scalar_field(state, trim(t%name)//"Source", stat = stat)
    have_source = stat == 0
    if(have_source) then
      assert(mesh_dim(source) == mesh_dim(t))
      assert(ele_count(source) == ele_count(t))
      
      add_src_directly_to_rhs = have_option(trim(source%option_path)//'/diagnostic/add_directly_to_rhs')
      
      if (add_src_directly_to_rhs) then 
         ewrite(2, *) "Adding Source field directly to the right hand side"
         assert(node_count(source) == node_count(t))
      end if
    
      ewrite_minmax(source)
    else
      ewrite(2,*) 'No source'
      
      add_src_directly_to_rhs = .false.
    end if

    ! Absorption
    absorption => extract_scalar_field(state, trim(t%name) // "Absorption", stat = stat)
    have_absorption = stat == 0
    if(have_absorption) then
      assert(mesh_dim(absorption) == mesh_dim(t))
      assert(ele_count(absorption) == ele_count(t))
    
      ewrite_minmax(absorption)
    else
      ewrite(2, *) "No absorption"
    end if

    ! Diffusivity
    diffusivity => extract_tensor_field(state, trim(t%name) // "Diffusivity", stat = stat)
    have_diffusivity = stat == 0
    if(have_diffusivity) then
      assert(all(diffusivity%dim == mesh_dim(t)))
      assert(ele_count(diffusivity) == ele_count(t))
      
      isotropic_diffusivity = option_count(complete_field_path(diffusivity%option_path)) &
        & == option_count(trim(complete_field_path(diffusivity%option_path)) // "/value/isotropic")
        
      if(isotropic_diffusivity) then
        ewrite(2, *) "Isotropic diffusivity"
        assert(all(diffusivity%dim > 0))
        ewrite_minmax(diffusivity%val(1, 1, :))
      else
        ewrite_minmax(diffusivity)
      end if
    else
      isotropic_diffusivity = .false.
      ewrite(2, *) "No diffusivity"
    end if

    ! field coordinate
    if(have_diffusivity) then
      t_coordinate = get_coordinate_field(state, t%mesh)
    else
      t_coordinate = coordinate
      call incref(t_coordinate)
    end if

    call get_option(trim(t%option_path) // "/prognostic/temporal_discretisation/theta", theta)
    assert(theta >= 0.0 .and. theta <= 1.0)
    ewrite(2, *) "Theta = ", theta
    dt_theta = dt*theta

    have_advection = .not. have_option(trim(t%option_path) // "/prognostic/spatial_discretisation/finite_volume/advection_terms/exclude_advection_terms")
    if(have_advection) then
      FLExit("Including advection not currently supported with FV")
      ewrite(2, *) "Including advection"
    else
      ewrite(2, *) "Excluding advection"
    end if
    
    have_mass = .not. have_option(trim(t%option_path) // "/prognostic/spatial_discretisation/finite_volume/mass_terms/exclude_mass_terms")
    if(have_mass) then
      ewrite(2, *) "Including mass"
    else
      ewrite(2, *) "Excluding mass"
    end if
    
    ! are we moving the mesh?
    move_mesh = (have_option("/mesh_adaptivity/mesh_movement") .and. have_mass)
    if(move_mesh) then
      FLExit("Moving the mesh not currently supported with FV")
      ewrite(2,*) "Moving the mesh"
      old_coordinate => extract_vector_field(state, "OldCoordinate")
      new_coordinate => extract_vector_field(state, "IteratedCoordinate")
      
      ! Grid velocity
      grid_velocity => extract_vector_field(state, "GridVelocity")
      assert(grid_velocity%dim == mesh_dim(t))
      assert(ele_count(grid_velocity) == ele_count(t))
      
      ewrite(2, *) "Grid velocity:"    
      ewrite_minmax(grid_velocity)
    else
      ewrite(2,*) "Not moving the mesh"
    end if
    
    call zero(matrix)
    call zero(rhs)
    

#ifdef _OPENMP
    cache_valid = prepopulate_transform_cache(coordinate)
#endif

    call get_mesh_colouring(state, T%mesh, COLOURING_DG1, colours)

    call profiler_tic(t, "advection_diffusion_fv_loop")

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(clr, len, nnid, ele, thread_num)

#ifdef _OPENMP    
    thread_num = omp_get_thread_num()
#else
    thread_num=0
#endif


    colour_loop: do clr = 1, size(colours)
      len = key_count(colours(clr))
      !$OMP DO SCHEDULE(STATIC)
      element_loop: do nnid = 1, len
       ele = fetch(colours(clr), nnid)
       call assemble_advection_diffusion_element_fv(ele, t, matrix, rhs, &
                                                   coordinate, t_coordinate, &
                                                   source, absorption, diffusivity)
      end do element_loop
       !$OMP END DO

    end do colour_loop
    !$OMP END PARALLEL

    call profiler_toc(t, "advection_diffusion_fv_loop")

    ! Add the source directly to the rhs if required 
    ! which must be included before dirichlet BC's.
    if (add_src_directly_to_rhs) call addto(rhs, source)

    ewrite(2, *) "Applying strong Dirichlet boundary conditions"
    call apply_dirichlet_conditions(matrix, rhs, t, dt)
    
    ewrite_minmax(rhs)
    call deallocate(t_coordinate)

    ewrite(1,*) "Exiting assemble_advection_diffusion_fv"

  end subroutine assemble_advection_diffusion_fv
  
  subroutine assemble_advection_diffusion_element_fv(ele, t, matrix, rhs, &
                                                   coordinate, t_coordinate, &
                                                   source, absorption, diffusivity)
                                                   
    integer, intent(in) :: ele
    type(scalar_field), intent(in) :: t
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    type(vector_field), intent(in) :: coordinate, t_coordinate
    type(scalar_field), intent(in) :: source
    type(scalar_field), intent(in) :: absorption
    type(tensor_field), intent(in) :: diffusivity
    
    integer, dimension(:), pointer :: element_nodes
    real, dimension(ele_ngi(t, ele)) :: detwei
    type(element_type), pointer :: t_shape

    ! What we will be adding to the matrix and RHS - assemble these as we
    ! go, so that we only do the calculations we really need
    real, dimension(ele_and_faces_loc(t, ele)) :: rhs_addto
    real, dimension(ele_and_faces_loc(t, ele), ele_and_faces_loc(t, ele)) :: matrix_addto
    integer, dimension(ele_and_faces_loc(t,ele)) :: local_glno
    
    integer, dimension(:), pointer :: neigh, x_neigh
    integer :: loc, ni, ele_2, face, face_2, start, finish
    
    assert(element_degree(t,ele)==0)
    t_shape => ele_shape(t, ele)

    loc = ele_loc(t, ele) ! how many nodes belong to this element... should be 1

    element_nodes => ele_nodes(t, ele)
    local_glno(:loc)=element_nodes
    
    matrix_addto = 0.0
    rhs_addto = 0.0

    if(have_mass.or.have_source.or.have_absorption) then
      call transform_to_physical(coordinate, ele, detwei=detwei)
    end if
    
    ! Mass
    if(have_mass) call add_mass_element_fv(ele, t_shape, t, detwei, matrix_addto(:loc,:loc))
    
    ! Absorption
    if(have_absorption) call add_absorption_element_fv(ele, t_shape, t, absorption, detwei, matrix_addto(:loc,:loc), rhs_addto(:loc))
    
    ! Source
    if(have_source .and. (.not. add_src_directly_to_rhs)) then 
       call add_source_element_fv(ele, t_shape, t, source, detwei, rhs_addto(:loc))
    end if
    
    if(have_diffusivity.or.have_advection) then

      ! this part of assembly is over faces of elements so it needs
      ! to know about boundary conditions
!       allocate(t_bc_types(surface_element_count(t)))
!       call get_entire_boundary_condition(t, (/ &
!         "neumann      ", &
!         "weakdirichlet", &
!         "internal     "/), t_bc, t_bc_types)

      neigh=>ele_neigh(t, ele)
      x_neigh => ele_neigh(coordinate, ele)
      
      start = size(element_nodes)+1
      
      neighboorloop: do ni = 1, size(neigh)
      
        ele_2=neigh(ni)
        face=ele_face(t, ele, ele_2)
        if(ele_2>0) then
          face_2=ele_face(t, ele_2, ele)
        else
          ! external face
          face_2 = face
        end if
        
        finish = start+face_loc(t, face_2)-1
        local_glno(start:finish) = face_global_nodes(t, face_2)
        
        call assemble_advection_diffusion_face_fv(face, face_2, start, finish, &
                                                  t, coordinate, t_coordinate, diffusivity, &
                                                  matrix_addto, rhs_addto)
                                                  
        start = start+face_loc(t, face_2)
        
      end do neighboorloop
      
!       call deallocate(t_bc)
!       deallocate(t_bc_types)
    
    end if
    
    if(have_diffusivity.or.have_advection) then
      ! element depends on its neighbours
      if(have_diffusivity) then
        ! have diffusivity, all neighbours involved
        call addto(matrix, local_glno, local_glno, matrix_addto)
        call addto(rhs, local_glno, rhs_addto)
      else
        ! no diffusivity so only some neighbours need adding to
        call addto(matrix, element_nodes, local_glno, matrix_addto(:loc,:))
        call addto(rhs, element_nodes, rhs_addto(:loc))
      end if
    else
      ! element doesn't depend on its neighbours
      call addto(matrix, element_nodes, element_nodes, matrix_addto(:loc,:loc))
      call addto(rhs, element_nodes, rhs_addto(:loc))
    end if

  end subroutine assemble_advection_diffusion_element_fv

  subroutine assemble_advection_diffusion_face_fv(face, face_2, start, finish, &
                                                  t, coordinate, t_coordinate, diffusivity, &
                                                  matrix_addto, rhs_addto)
  
    integer, intent(in) :: face, face_2, start, finish
    type(scalar_field), intent(in) :: t
    type(vector_field), intent(in) :: coordinate, t_coordinate
    type(tensor_field), intent(in) :: diffusivity
    real, dimension(:,:), intent(inout) :: matrix_addto
    real, dimension(:), intent(inout) :: rhs_addto
    
    real, dimension(face_ngi(t, face)) :: detwei
    real, dimension(mesh_dim(t), face_ngi(t, face)) :: normal
    
    
    call transform_facet_to_physical(coordinate, face, &
                                     detwei_f=detwei, normal=normal)
                                 
    if(have_diffusivity) then
      call add_diffusivity_face_fv(face, face_2, start, finish, &
                                   t, t_coordinate, diffusivity, &
                                   detwei, normal, &
                                   matrix_addto, rhs_addto)
    end if
                                     
  
  end subroutine assemble_advection_diffusion_face_fv
  
  subroutine add_diffusivity_face_fv(face, face_2, start, finish, &
                                     t, t_coordinate, diffusivity, &
                                     detwei, normal, &
                                     matrix_addto, rhs_addto)
    
    integer, intent(in) :: face, face_2, start, finish
    type(scalar_field), intent(in) :: t
    type(vector_field), intent(in) :: t_coordinate
    type(tensor_field), intent(in) :: diffusivity
    real, dimension(:), intent(in) :: detwei
    real, dimension(:,:), intent(in) :: normal
    real, dimension(:,:), intent(inout) :: matrix_addto
    real, dimension(:), intent(inout) :: rhs_addto
    
    real, dimension(face_loc(t_coordinate, face), t_coordinate%dim) :: c_vector
    real, dimension(face_loc(t_coordinate, face)) :: c_dist
    real, dimension(diffusivity%dim(1), diffusivity%dim(2), face_ngi(diffusivity, face)) :: diffusivity_gi
    real, dimension(face_loc(t, face)+face_loc(t,face_2), face_ngi(t, face), mesh_dim(t)) :: dt_t
    type(element_type), pointer :: t_shape
    real, dimension(face_loc(t, face), face_loc(t,face)+face_loc(t,face_2)) :: diff_mat
    integer, dimension(face_loc(t, face)) :: t_face_l
    real, dimension(face_loc(t, face) + face_loc(t, face_2)) :: t_val
    integer :: loc, tloc, iloc, jloc, gi

    assert(face_loc(t, face)==face_loc(t_coordinate, face))
    assert(face_loc(t, face)==face_loc(t, face_2))
    
    loc = face_loc(t,face) ! should just be 1 but being paranoid
    tloc = face_loc(t, face) + face_loc(t, face_2) ! so clearly should be 2
    
    t_face_l = face_local_nodes(t, face)
    t_val(:loc) = face_val(t, face)
    t_val(loc+1:) = face_val(t, face_2)
    
    ! first we need to construct the derivative of a pseudo "shape function"
    dt_t = 0.0
    
    if(face==face_2) then
      ! on boundary
    
    else
      ! internal (but might be a periodic boundary)
      
      c_vector = transpose(face_val(t_coordinate, face) - face_val(t_coordinate, face_2))
      c_dist = sum(c_vector**2, dim=2) ! this is the square of the distance
      do iloc = 1, loc
        ! normalise
        c_vector(iloc,:) = c_vector(iloc,:)/c_dist(iloc)
      end do

      dt_t(:loc,:,:) = spread(c_vector, 2, face_ngi(t,face))
      dt_t(loc+1:,:,:) = spread(-c_vector, 2, face_ngi(t,face))
      
      ! now we need to construct the matrix entries for the face integral:
      ! /
      ! | shape (diffusivity dt_t) . normal ds
      ! /
      diff_mat = 0.0
      
      diffusivity_gi = face_val_at_quad(diffusivity, face)
      t_shape => face_shape(t, face)
      do iloc = 1, face_loc(t, face) ! just the node at the centre of this element
        do jloc = 1, (face_loc(t, face)+face_loc(t, face_2)) ! this element's node plus the neighbouring element's
          do gi = 1, face_ngi(t, face) ! all the gauss points on this face
            diff_mat(iloc, jloc) = diff_mat(iloc, jloc) - t_shape%n(iloc, gi)*&
                                  sum(matmul(diffusivity_gi(:,:,gi), dt_t(jloc, gi, :))*normal(:,gi), 1)*&
                                  detwei(gi)
          end do
        end do
      end do

      if(abs(dt_theta) > epsilon(0.0)) then
        matrix_addto(t_face_l, t_face_l) = matrix_addto(t_face_l, t_face_l) + dt_theta*diff_mat(:,:loc)
        matrix_addto(t_face_l, start:finish) = matrix_addto(t_face_l, start:finish) + dt_theta*diff_mat(:,loc+1:)
      end if
      rhs_addto(t_face_l) = rhs_addto(t_face_l) - matmul(diff_mat, t_val)

    end if
    
  end subroutine add_diffusivity_face_fv

  subroutine add_mass_element_fv(ele, t_shape, t, detwei, matrix_addto)
    integer, intent(in) :: ele
    type(element_type), intent(in) :: t_shape
    type(scalar_field), intent(in) :: t
    real, dimension(ele_ngi(t, ele)), intent(in) :: detwei
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)), intent(inout) :: matrix_addto
    
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)) :: mass_matrix
    
    assert(have_mass)
    
    mass_matrix = shape_shape(t_shape, t_shape, detwei)
    
    matrix_addto = matrix_addto + mass_matrix
  
  end subroutine add_mass_element_fv

  subroutine add_absorption_element_fv(ele, t_shape, t, absorption, detwei, matrix_addto, rhs_addto)
    integer, intent(in) :: ele
    type(element_type), intent(in) :: t_shape
    type(scalar_field), intent(in) :: t
    type(scalar_field), intent(in) :: absorption
    real, dimension(ele_ngi(t, ele)), intent(in) :: detwei
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)), intent(inout) :: matrix_addto
    real, dimension(ele_loc(t, ele)), intent(inout) :: rhs_addto
    
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)) ::  absorption_mat
    
    assert(have_absorption)
    
    absorption_mat = shape_shape(t_shape, t_shape, detwei * ele_val_at_quad(absorption, ele))
    
    if(abs(dt_theta) > epsilon(0.0)) matrix_addto = matrix_addto + dt_theta * absorption_mat
    
    rhs_addto = rhs_addto - matmul(absorption_mat, ele_val(t, ele))
    
  end subroutine add_absorption_element_fv

  subroutine add_source_element_fv(ele, t_shape, t, source, detwei, rhs_addto)
    integer, intent(in) :: ele
    type(element_type), intent(in) :: t_shape
    type(scalar_field), intent(in) :: t
    type(scalar_field), intent(in) :: source
    real, dimension(ele_ngi(t, ele)), intent(in) :: detwei
    real, dimension(ele_loc(t, ele)), intent(inout) :: rhs_addto
   
    assert(have_source)
   
    rhs_addto = rhs_addto + shape_rhs(t_shape, detwei * ele_val_at_quad(source, ele))
    
  end subroutine add_source_element_fv
  
  subroutine advection_diffusion_fv_check_options
  
    character(len = FIELD_NAME_LEN) :: field_name, mesh_0, mesh_1, state_name
    character(len = OPTION_PATH_LEN) :: path
    integer :: i, j, stat
    real :: beta, l_theta
    
    if(option_count("/material_phase/scalar_field/prognostic/spatial_discretisation/finite_volume") == 0) then
      ! Nothing to check
      return
    end if
        
    ewrite(2, *) "Checking FV advection-diffusion options"
    
    do i = 0, option_count("/material_phase") - 1
      path = "/material_phase[" // int2str(i) // "]"
      call get_option(trim(path) // "/name", state_name)
      
      do j = 0, option_count(trim(path) // "/scalar_field") - 1
        path = "/material_phase[" // int2str(i) // "]/scalar_field[" // int2str(j) // "]"
        call get_option(trim(path) // "/name", field_name)
        
        if(field_name /= "Pressure") then
        
          path = trim(path) // "/prognostic"
          
          if(have_option(trim(path) // "/spatial_discretisation/finite_volume").and.&
             have_option(trim(path) // "/equation[0]")) then       
             
            call get_option(trim(path) // "/spatial_discretisation/conservative_advection", beta, stat)
            if(stat == SPUD_NO_ERROR) then
              if(beta < 0.0 .or. beta > 1.0) then
              
                call field_error(state_name, field_name, &
                  & "Conservative advection factor (beta) must be >= 0.0 and <= 1.0")
              end if
            else
              call field_error(state_name, field_name, &
                & "Conservative advection factor (beta) required")
            end if
            
            call get_option(trim(path) // "/temporal_discretisation/theta", l_theta, stat)
            if(stat == SPUD_NO_ERROR) then
              if(l_theta < 0. .or. l_theta > 1.0) then
                call field_error(state_name, field_name, &
                  &"Implicitness factor (theta) must be >= 0.0 and <= 1.0")
              end if
            else
              call field_error(state_name, field_name, &
                & "Implicitness factor (theta) required")
            end if
            if(have_option(trim(path) // "/spatial_discretisation/finite_volume/mass_terms/exclude_mass_terms") .and. &
              & abs(l_theta - 1.0) > epsilon(0.0)) then
              call field_warning(state_name, field_name, &
                & "Implicitness factor (theta) should = 1.0 when excluding mass")
            end if
                 
            if (have_option(trim(path) // "/scalar_field::SinkingVelocity")) then
               call get_option(trim(complete_field_path(trim(path) // &
                    "/scalar_field::SinkingVelocity"))//"/mesh[0]/name", &
                    mesh_0, stat)
               if(stat == SPUD_NO_ERROR) then
                  call get_option(trim(complete_field_path("/material_phase[" // int2str(i) // &
                       "]/vector_field::Velocity")) // "/mesh[0]/name", mesh_1)
                  if(trim(mesh_0) /= trim(mesh_1)) then
                     call field_warning(state_name, field_name, &
                          & "SinkingVelocity is on a different mesh to the Velocity field this could cause problems")
                  end if
               end if
            end if
            if(have_option(trim(path) // "/spatial_discretisation/finite_volume/advection_terms/exclude_advection_terms")) then
              if(have_option(trim(path) // "/scalar_field::SinkingVelocity")) then
                call field_warning(state_name, field_name, &
                  & "SinkingVelocity set, but advection terms have been excluded - SinkingVelocity will have no effect")
              end if
            end if
  
            if(option_count(trim(path) // "/boundary_conditions/type::neumann") > 0 &
              & .and. .not. (have_option(trim(path) // "/tensor_field::Diffusivity") &
              & .or. have_option(trim(path) // "/subgridscale_parameterisation::k-epsilon") &
              & .or. have_option(trim(path) // "/subgridscale_parameterisation::GLS"))) then
                call field_warning(state_name, field_name, &
                & "Neumann boundary condition set, but have no diffusivity - boundary condition will not be applied")
            end if
          end if
        end if
      end do
    end do
    
    ewrite(2, *) "Finished checking CG advection-diffusion options"

    contains
  
    subroutine field_warning(state_name, field_name, msg)
      character(len = *), intent(in) :: state_name
      character(len = *), intent(in) :: field_name
      character(len = *), intent(in) :: msg
      
      ewrite(0, *) "Warning: For field " // trim(field_name) // " in state " // trim(state_name)
      ewrite(0, *) trim(msg)
    
    end subroutine field_warning
  
    subroutine field_error(state_name, field_name, msg)
      character(len = *), intent(in) :: state_name
      character(len = *), intent(in) :: field_name
      character(len = *), intent(in) :: msg
      
      ewrite(-1, *) "For field " // trim(field_name) // " in state " // trim(state_name)
      FLExit(trim(msg))
    
    end subroutine field_error
  
  end subroutine advection_diffusion_fv_check_options


end module advection_diffusion_fv
