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
  use data_structures
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

    type(mesh_type), pointer:: p_mesh, l_mesh
    type(scalar_field), pointer:: bottom_distance
    type(vector_field), pointer:: gravity_normal
    type(vector_field):: ugravity_normal
    type(integer_set):: column_nodes
    real, dimension(u%dim):: uvec, unorm
    integer, dimension(:) ,pointer:: surface_element_list
    integer :: i, sele, ele

    p_mesh => extract_pressure_mesh(state)
    if (.not. (element_degree(p_mesh,1)==element_degree(u,1)+1 .and. continuity(p_mesh)>=0 .and. continuity(u)<0)) then
      FLExit("Hydrostatic pressure projection only works for PnDG-Pn+1")
    end if
    call find_linear_parent_mesh(state, p_mesh, l_mesh)
    if (.not. associated(l_mesh%columns)) then
      ! this could be made to work without the %columns, but that requires computing
      ! face normals, which is expensive.
      FLExit("Hydrostatic pressure projection only works for extruded meshes")
    end if

    ! FIXME: do something with the rhs (also this assert is obv. wrong)
    assert(maxval(full_ct_rhs)==0.0)

    gravity_normal => extract_vector_field(state, "GravityDirection")
    call allocate(ugravity_normal, u%dim, u%mesh, "UGravityDirection")
    call remap_field(gravity_normal, ugravity_normal)

    ! make sure we start out with a zero vertical velocity, by projecting out that component
    do i=1, node_count(u)
      unorm=node_val(ugravity_normal, i)
      uvec=node_val(u,i)
      call set(u, i, uvec - dot_product(unorm,uvec)*unorm)
    end do

    ! we start from the bottom - its surface mesh is stored as a boundary condition under "DistanceToBottom"
    bottom_distance => extract_scalar_field(state, "DistanceToBottom")
    call get_boundary_condition(bottom_distance, 1, &
        surface_element_list=surface_element_list)
    do i=1, size(surface_element_list)
      sele = surface_element_list(i)
      ! the set of velocity nodes in this column, we've solved already:
      call allocate(column_nodes)
      do
        ele=face_ele(p_mesh, sele)
        ! solve the velocities in the element ele above sele
        ! this routine returns the sele that faces the element above
        ! or returns sele=-1 if we've reached the top
        call vertical_reconstruction_ele(column_nodes, ele, sele)

        if (sele<0) exit

      end do
      call deallocate(column_nodes)

    end do

    call deallocate(ugravity_normal)

    contains

    subroutine vertical_reconstruction_ele(column_nodes, ele, sele)
      type(integer_set), intent(inout):: column_nodes
      integer, intent(in):: ele
      integer, intent(inout):: sele

      integer, dimension(face_loc(l_mesh,sele)):: flnodes
      type(integer_set):: fnodes
      type(real_vector), dimension(u%dim):: rowvals
      real, dimension(ele_loc(p_mesh,ele)):: rhs
      real, dimension(size(rhs), size(rhs)):: matrix
      real, dimension(u%dim):: unorm
      integer, dimension(:), pointer:: faces, pnodes, unodes, row
      integer:: j, k, k1, m, p, next_sele

      ! find the face between this element and the next in the column
      ! this should be an element whose nodes are all in different columns
      faces => ele_faces(l_mesh, ele)
      do j=1, size(faces)
        if (faces(j)==sele) cycle
        flnodes=face_global_nodes(l_mesh, faces(j))
        if (all_different(l_mesh%columns(flnodes))) exit
      end do
      if (j>size(faces)) then
        ! no such face found, something is wrong
        ewrite(-1,*) "It seems the mesh is not columnar"
        ! how did we get here if %columns is present?
        FLExit("Hydrostatic pressure projection requires columnar mesh.")
      end if
      next_sele=faces(j)

      pnodes => ele_nodes(p_mesh, ele)
      unodes => ele_nodes(u, ele)
      call allocate(fnodes)
      call insert(fnodes, face_global_nodes(p_mesh, next_sele))
      call insert(column_nodes, unodes)

      p=0
pressure_node_loop: do j=1, size(pnodes)
        ! if this node is on the next face, we'll deal with in the next element
        if (has_value(fnodes, pnodes(j))) cycle
        p=p+1 ! count the pressure nodes we do deal with

        ! get the column indices, and matrix values for this pressure row
        row => row_m_ptr(full_ct_m, pnodes(j))
        do m=1, u%dim
          rowvals(m)%ptr => row_val_ptr(full_ct_m, 1, m, pnodes(j))
        end do

        ! Multiply current velocity with ct_m and put on rhs. Only velocity dofs
        ! that we've visisted before within this column are included. By the choice of pressure
        ! points within this element these are the only u-nodes that we should encounter in this
        ! matrix row, so that effectively we're computing the restriction of this pressure test
        ! function to the column. Note that this includes u-nodes within this element, of which
        ! the vertical component has been zeroed, for the other nodes, strictly below the element
        ! this vertical component has been computed before and is included in the rhs.
        ! Also in this loop, we find the start of the u-nodes of this element itelf within row(:) - assuming sorted rows
        rhs(p)=0
        do k=1, size(row)
          if (has_value(column_nodes, row(k))) then
            do m=1, u%dim
              rhs(p) = rhs(p) - rowvals(m)%ptr(k)*node_val(u, m, row(k))
            end do
            if (row(k)==unodes(1)) k1=k-1
          end if
        end do

        ! the row of the local matrix associated with this pressure node
        ! is computed by taking the inner product of the normal with ct_m
        ! restricting to only the u-nodes of this element
        do k=1, size(unodes)
          unorm=node_val(ugravity_normal, unodes(k))
          assert( row(k1+k)==unodes(k) )
          matrix(p, k)=0.0
          do m=1, u%dim
            matrix(p, k)=matrix(p, k) - unorm(m)*rowvals(m)%ptr(k1+k)
          end do
        end do

      end do pressure_node_loop

      ! make sure the number of pressure-tested eqns equals the number of velocity dofs to solve for
      assert( p==ele_loc(u, ele) )

      call solve(matrix, rhs)

      ! the answer is returned in rhs
      ! multiply with the upward normal and add in to the velocity:
      do k=1, size(unodes)
        unorm=node_val(ugravity_normal, unodes(k))
        call addto(u, unodes(k), -unorm*rhs(k))
      end do

      call deallocate(fnodes)

      sele=face_opposite(p_mesh, next_sele)

    end subroutine vertical_reconstruction_ele

    logical function all_different(array)
      ! checks whether all items in the array are different
      integer, dimension(:), intent(in):: array

      integer:: i
      ! compare each with next
      do i=1, size(array)-1
        if (array(i)==array(i+1)) then
          all_different=.false.
          return
        end if
      end do
      
      ! finally, compare last with first
      all_different=array(i)/=array(1)

    end function all_different

  end subroutine reconstruct_vertical_velocities
  
end module hydrostatic_projection
