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
module hydrostatic_projection
  use fldebug
  use sparse_tools
  use fields
  use state_module
  use sparse_matrices_fields
  use sparsity_patterns_meshes
  use boundary_conditions
  use vtk_interfaces
  use solvers
  use field_options
  use vertical_extrapolation_module
  use divergence_matrix_cg
  
implicit none

contains

  subroutine prepare_hydrostatic_projection(state, surface_p, surface_old_p, lumped_cmc_m, get_cmc_m)
    type(state_type), intent(inout):: state
    type(scalar_field), intent(out):: surface_p, surface_old_p
    type(csr_matrix), pointer:: lumped_cmc_m
    logical, intent(out):: get_cmc_m
    
    type(vector_field), pointer:: positions, vertical_normal
    type(scalar_field), pointer:: topdis, p, old_p
    type(vector_field):: surface_x
    type(mesh_type), pointer:: p_mesh, u_mesh, surface_p_mesh
    type(mesh_type):: local_surface_p_mesh
    type(csr_matrix) :: local_lumped_cmc_m
    type(csr_matrix) :: prolongator
    type(csr_sparsity):: pcmc_sparsity, pcmcp_sparsity
    type(csr_sparsity), pointer:: cmc_sparsity
    integer, dimension(:), pointer:: surface_element_list, surface_nodes
    integer:: stat
    
    p_mesh => extract_pressure_mesh(state)
    u_mesh => extract_pressure_mesh(state)
    
    ! retreive the list of surface elements for the surface_mesh
    topdis => extract_scalar_field(state, "DistanceToTop", stat=stat)
    if (stat/=0) then
       FLExit("For hydrostatic pressure projection you need to specify the ocean_boundaries under /geometry")
    end if    
    call get_boundary_condition(topdis, 1, surface_element_list=surface_element_list)        

    surface_p_mesh => extract_mesh(state, "Lumped"//trim(p_mesh%name), stat=stat)
    if (stat/=0) then
      
      call create_surface_mesh(local_surface_p_mesh, surface_nodes, p_mesh, surface_element_list, &
        name="Lumped"//trim(p_mesh%name))
      call insert(state, local_surface_p_mesh, local_surface_p_mesh%name)
      call deallocate(local_surface_p_mesh)
      surface_p_mesh => extract_mesh(state, "Lumped"//trim(p_mesh%name))
      
      positions => extract_vector_field(state, "Coordinate")
      call allocate(surface_x, positions%dim, surface_p_mesh, "SurfacePositions")
      call remap_field_to_surface(positions, surface_x, surface_element_list)
      call vtk_write_fields("surface_mesh", position=surface_x, model=surface_p_mesh)
      call deallocate(surface_x)
      vertical_normal => extract_vector_field(state, "GravityDirection")

      assert(.not. has_csr_matrix(state, "HydrostaticProlongator"))
      ! note, have to provide surface_mesh to ensure from_mesh of the prolongator uses same node ordering
      prolongator = VerticalProlongationOperator( &
           p_mesh, positions, vertical_normal, surface_element_list, &
           surface_mesh=surface_p_mesh)
      prolongator%name = "HydrostaticProlongator"
      call insert(state, prolongator, prolongator%name)
      call deallocate(prolongator)
       
      assert(.not. has_csr_matrix(state, "PressurePoissonMatrix"))
      
      ! first obtain the right sparsity, we start with the cmc sparsity on the full mesh:
      cmc_sparsity => get_csr_sparsity_secondorder(state, p_mesh, u_mesh)
      ! then we compute P^T * (C^T M^-1 C)
      pcmc_sparsity = matmul_atb(prolongator%sparsity, cmc_sparsity)
      ! and finally: (P^T C^T M^-1) * P
      pcmcp_sparsity = matmul(pcmc_sparsity, prolongator%sparsity)
      pcmcp_sparsity%name = "Lumped"//trim(cmc_sparsity%name)
      ! note, that this sparsity in general is not the same as the one we would have have obtained
      ! from doing the usual get_csr_sparsity_secondorder() using a surface p and surface u mesh
      ! The latter would only be correct if the mesh is extruded strictly along prismatic columns 
      ! (as in fluidity's extrusion algorithm) - even for a fully "structured" gmsh tetrahedral mesh they are different
      
      call allocate(local_lumped_cmc_m, pcmcp_sparsity, name="PressurePoissonMatrix")
      call insert(state, local_lumped_cmc_m, name="PressurePoissonMatrix")
      call deallocate(local_lumped_cmc_m)
      call deallocate(pcmc_sparsity)
      call deallocate(pcmcp_sparsity)
      get_cmc_m = .true.
      
      
    else
       
      get_cmc_m = .false.
            
    end if
    
    call allocate(surface_p, surface_p_mesh, name="SurfacePressure")
    p => extract_scalar_field(state, "Pressure")
    call remap_field_to_surface(p, surface_p, surface_element_list)
    surface_p%option_path=p%option_path
    
    call allocate(surface_old_p, surface_p_mesh, name="SurfaceOldPressure")
    old_p => extract_scalar_field(state, "OldPressure", stat=stat)
    if (stat/=0) old_p => p
    call remap_field_to_surface(old_p, surface_old_p, surface_element_list)
    
    lumped_cmc_m => extract_csr_matrix(state, "PressurePoissonMatrix")
    
  end subroutine prepare_hydrostatic_projection
  
  subroutine get_lumped_continuity_equation(state, lumped_p_mesh, ct_m, ct_rhs, get_ct_m, lumped_ct_m, lumped_ct_rhs)
    type(state_type), intent(inout):: state
    type(mesh_type), intent(inout):: lumped_p_mesh
    type(block_csr_matrix), intent(in):: ct_m
    type(scalar_field), intent(in):: ct_rhs
    logical, intent(in):: get_ct_m
    type(block_csr_matrix), pointer:: lumped_ct_m
    type(scalar_field), intent(out):: lumped_ct_rhs
      
    type(csr_matrix), pointer:: prolongator
    type(block_csr_matrix):: local_lumped_ct_m
    logical:: compute_lumped_ct_m
    
    ! compute or retreive the lumped version of the divergence matrix
    ! if ct_m hasn't changed (get_ct_m==.false.) we should be able to reuse the one cached in state
    compute_lumped_ct_m = get_ct_m
    
    prolongator => extract_csr_matrix(state, "HydrostaticProlongator")
    
    if (.not. has_block_csr_matrix(state, "Lumped"//trim(ct_m%name))) then
      
      if (compute_lumped_ct_m) then
        
        local_lumped_ct_m = matmul_atb(prolongator, ct_m)
        local_lumped_ct_m%name = "Lumped"//trim(ct_m%name)
        ! stick it in state
        call insert(state, local_lumped_ct_m, "Lumped"//trim(ct_m%name))
        ! now we can drop our reference
        call deallocate(local_lumped_ct_m)
        ! no need to compute it again now
        compute_lumped_ct_m = .false.
        
      else
        
        ! Something went wrong here, apparently ct_m has not changed but we've
        ! managed to loose the lumped version of it somehow
        FLAbort("Internal error: Missing lumped divergence matrix in state")
        
      end if
    end if
      
    lumped_ct_m => extract_block_csr_matrix(state, "Lumped"//trim(ct_m%name))
    
    if (compute_lumped_ct_m) then
      call zero(lumped_ct_m)
      call matmul_atb_addto(prolongator, ct_m, lumped_ct_m)
    end if
    
    call allocate(lumped_ct_rhs, lumped_p_mesh, "LumpedDivergenceRHS")
    call mult_T(lumped_ct_rhs, prolongator, ct_rhs)
    
  end subroutine get_lumped_continuity_equation
    
  subroutine reconstruct_vertical_velocities(state, u, full_ct_m, full_ct_rhs, get_ct_m)
    type(state_type), intent(inout):: state
    type(vector_field), intent(inout):: u
    type(block_csr_matrix), intent(in):: full_ct_m
    type(scalar_field):: full_ct_rhs
    logical, intent(in):: get_ct_m
    
    type(block_csr_matrix):: pct_m
    type(csr_matrix):: temp_vct_m
    type(csr_matrix), pointer:: vct_m
    type(csr_sparsity), pointer:: vct_sparsity
    type(scalar_field), pointer:: p
    type(scalar_field):: vertical_u, component
    type(vector_field), pointer:: vertical_normal, positions
    type(vector_field):: pu
    real, dimension(u%dim):: normal
    logical:: constant_normal, get_vct_m
    integer:: i, ele, stat, vdim
    
    vertical_normal => extract_vector_field(state, "GravityDirection")
    ! if gravity is in a constant, axis-aligned direction - use the cheap version
    normal = node_val(vertical_normal,1)
    constant_normal = vertical_normal%field_type==FIELD_TYPE_CONSTANT .and. &
                           maxval(abs(normal))==sum(abs(normal))
                           
    if (constant_normal) then
      
      vdim=maxloc(node_val(vertical_normal,1),1)
      do i=1, u%dim
        if (i/=vdim) then
          component = extract_scalar_field(u, i)
          ! move the contribution of already computed horizontal components
          ! to the rhs of the continuity eqn
          call mult_addto(full_ct_rhs, block(full_ct_m, 1, i), component, scale=-1.0)
        end if
      end do
      
    else
    
      FLExit("The hydrostatic projection does not yet work for varying gravitational direction (e.g. on the sphere).")
      
    end if

    p => extract_scalar_field(state, "Pressure")
    
    vct_m => extract_csr_matrix(state, "VerticalVelocityGradientMatrix", stat=stat)
    ! do we need to reassemble it - should be the logic as for full ct_m
    get_vct_m = get_ct_m
    if (stat/=0) then
      
      vct_sparsity => get_csr_sparsity_firstorder(state, p%mesh, p%mesh)
      
      call allocate(temp_vct_m, vct_sparsity, name="VerticalVelocityGradientMatrix")
      call insert(state, temp_vct_m, name="VerticalVelocityGradientMatrix")
      call deallocate(temp_vct_m)
      
      vct_m => extract_csr_matrix(state, "VerticalVelocityGradientMatrix")
      
      get_vct_m =.true.
      
    end if
    
    if (get_vct_m) then
      ! used in assemble_vct_ele
      positions => extract_vector_field(state, "Coordinate")

      call allocate(pct_m, vct_m%sparsity, (/ 1, u%dim /), name="BlockVelocityGradientMatrix")

      call allocate(pu, u%dim, p%mesh, "VelocityOnPressureMesh")

      call assemble_divergence_matrix_cg(pct_m, state, full_ct_rhs, &
              p%mesh, pu, option_path=p%option_path)

      vct_m%val=pct_m%val(1,vdim)%ptr
      
!      do ele=1, element_count(p)
!        call assemble_vct_ele(vct_m, p%mesh, ele)
!      end do
      
    end if
    
    call allocate(vertical_u, p%mesh, "VerticalVelocity")
    vertical_u%option_path=p%option_path
    call petsc_solve(vertical_u, vct_m, full_ct_rhs)
    
  contains
    
    subroutine assemble_vct_ele(vct_m, p_mesh, ele)
      type(csr_matrix), intent(inout):: vct_m
      type(mesh_type), intent(in):: p_mesh
      integer, intent(in):: ele
      
      real, dimension(ele_loc(p_mesh, ele), ele_ngi(p_mesh, ele), mesh_dim(p_mesh)):: dshape
      real, dimension(ele_ngi(p_mesh, ele)):: detwei
      integer, dimension(:), pointer:: nodes
      
      assert( constant_normal ) ! see above - this means we know vdim
      
      call transform_to_physical(positions, ele, ele_shape(p_mesh, ele), dshape, detwei)
      nodes => ele_nodes(p_mesh, ele)
      FLAbort("Fix the next line first")
      !call addto(vct_m, nodes, nodes, shape_vector_dot_dshape(ele_shape(p_mesh,ele), dshape, detwei))
      
    end subroutine assemble_vct_ele
    
  end subroutine reconstruct_vertical_velocities
  
end module hydrostatic_projection
