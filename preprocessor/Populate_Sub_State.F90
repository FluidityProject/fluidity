!    Copyright (C) 2007 Imperial College London and others.
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
module populate_sub_state_module
  use elements
  use state_module
  use FLDebug
  use spud
  use mesh_files
  use vtk_cache_module
  use vtk_interfaces       
  use global_parameters, only: OPTION_PATH_LEN, is_active_process, pi, &
    no_active_processes, topology_mesh_name, adaptivity_mesh_name, &
    periodic_boundary_option_path, domain_bbox, domain_volume
  use field_options
  use reserve_state_module
  use fields_manipulation
  use diagnostic_variables, only: convergence_field, steady_state_field
  use field_options
  use surfacelabels
  use climatology
  use metric_tools
  use coordinates
  use halos
  use tictoc
  use hadapt_extrude
  use hadapt_extrude_radially
  use initialise_fields_module
  use transform_elements
  use parallel_tools
  use boundary_conditions_from_options
  use nemo_states_module
  use data_structures
  use fields_halos
  use read_triangle
  use populate_state_module

  implicit none

  private

  public populate_sub_state, set_full_domain_prescribed_fields, &
      &  sub_state_remap_to_full_mesh, update_subdomain_fields, &
      &  use_sub_state

contains

  subroutine populate_sub_state(states,sub_states)
    ! This routine initialises all meshes, fields and boundary conditions
    ! that will be used in partially prognostic solves - i.e. where certain
    ! variables are only solved for in part of the computational domain.

    use Profiler

    type(state_type), intent(in), dimension(:) :: states
    type(state_type), pointer, dimension(:) :: sub_states

    integer :: nstates ! number of states
    integer :: istate

    ewrite(1,*) "In populate_sub_state"
    call profiler_tic("I/O")
    call tictoc_clear(TICTOC_ID_IO_READ)

    ! Find out how many states there are and initialise:
    nstates=size(states)
    allocate(sub_states(1:nstates))
    do istate = 1, nstates
       call nullify(sub_states(istate))
    end do

    ! Form subdomain external mesh from full state:
    call derive_external_subdomain_mesh(states,sub_states)

    ! Derive other subdomain meshes from subdomain external mesh:
    call insert_derived_meshes(sub_states)

    ! Determine mapping functions for derived meshes, to and from full mesh:
    call insert_subdomain_mesh_maps(states,sub_states)

    call allocate_and_insert_fields(sub_states)

    call populate_boundary_conditions(sub_states)

    call set_boundary_conditions_values(sub_states)

    call set_dirichlet_consistent(sub_states)

    call alias_fields(sub_states)

    call allocate_and_insert_auxilliary_fields(sub_states)       

    ! Update field values in subdomain from state values:
    call update_subdomain_fields(states,sub_states)

    ! Set prescribed fields on full mesh:
    call set_full_domain_prescribed_fields(states)

    call tictoc_report(2, TICTOC_ID_IO_READ)
    call profiler_toc("I/O")
    ewrite(1, *) "Exiting populate_sub_state"

  end subroutine populate_sub_state

  logical function use_sub_state()
    ! Routine to determine whether or not sub_state is set up:

    integer :: number_of_prescribed_regions

    number_of_prescribed_regions = option_count("/material_phase/vector_field::Velocity/prognostic/prescribed_region")

    use_sub_state = (number_of_prescribed_regions > 0)

  end function use_sub_state

  subroutine derive_external_subdomain_mesh(states,sub_states)

    ! This routine derives a subdomain mesh equivalent
    ! to the externally derived mesh. This is later used as a basis
    ! for deriving all other meshes on the prognostic subdomain.

    type(state_type), intent(in), dimension(:) :: states
    type(state_type), intent(inout), dimension(:) :: sub_states

    ! Integer sets containing list of elements and nodes on subdomain_mesh:
    type(integer_set) :: subdomain_mesh_element_list, subdomain_mesh_node_list
    ! Element mapping functions to subdomain_mesh:
    integer, dimension(:), allocatable :: subele_list
    ! Node mapping functions to and from subdomain_mesh:
    integer, dimension(:), allocatable :: n_list, inverse_n_list

    ! External mesh and subdomain_meshes:
    type(mesh_type) :: subdomain_mesh
    type(mesh_type), pointer :: external_mesh
    character(len=FIELD_NAME_LEN) :: mesh_name

    ! Others:
    integer :: ele, ele_2, ni, edge_count, i, node, loc, sloc, face
    integer, dimension(2) :: prescribed_regions_shape
    integer :: number_of_prescribed_regions
    integer, dimension(:), pointer :: neigh
    type(integer_set) :: prescribed_region_id_set
    integer, dimension(:), allocatable :: prescribed_region_ids

    type(element_type), pointer :: shape     

    type(integer_set) :: face_list
    type(vector_field), pointer :: external_mesh_position
    type(vector_field) :: position

    type(vector_field), pointer :: velocity

    integer, allocatable, dimension(:) :: sndglno, boundary_ids, element_owner

    ewrite(1,*) "Entering derive external subdomain mesh"

    ! Create subdomain_meshes -- begin with external mesh as we must add faces etc... to this before
    ! deriving other meshes:
    external_mesh => get_external_mesh(states)
    mesh_name = external_mesh%name

    velocity => extract_vector_field(states(1),"Velocity")

    number_of_prescribed_regions = &
         option_count(trim(velocity%option_path)// "/prognostic/prescribed_region")
    ewrite(2,*) 'Number of prescribed_regions',number_of_prescribed_regions
    
    call allocate(prescribed_region_id_set)
    do i = 1, number_of_prescribed_regions
      prescribed_regions_shape = option_shape(trim(velocity%option_path)// &
                                               "/prognostic/prescribed_region["//int2str(i-1)//"]/region_ids")
      allocate(prescribed_region_ids(prescribed_regions_shape(1)))
      call get_option(trim(velocity%option_path)// &
                     "/prognostic/prescribed_region["//int2str(i-1)//"]/region_ids", prescribed_region_ids)

      call insert(prescribed_region_id_set, prescribed_region_ids)

      deallocate(prescribed_region_ids)
    end do

    ! Derive subdomain_mesh_element list (an integer set):
    call allocate(subdomain_mesh_element_list)
    do ele = 1, element_count(external_mesh)
       if(.not.has_value(prescribed_region_id_set, external_mesh%region_ids(ele))) then
          call insert(subdomain_mesh_element_list,ele)
       end if
    end do

    ! Build element mapping functions to and from sub mesh:
    allocate(subele_list(key_count(subdomain_mesh_element_list))) ! Map from sub mesh --> full mesh

    ! Initialise element mapping function to zero.
    subele_list = 0

    ! Set up element map from sub mesh --> full mesh:
    do i = 1, key_count(subdomain_mesh_element_list)
       subele_list(i) = fetch(subdomain_mesh_element_list,i)
    end do

    ewrite(1,*) 'Number of elements in subdomain_mesh_element_list:', size(subele_list)

    ! Create subdomain_mesh:
    call allocate(subdomain_mesh_node_list)

    ! Derive node list for subdomain_mesh:
    do i = 1, size(subele_list)
       ele = subele_list(i)
       call insert(subdomain_mesh_node_list,ele_nodes(external_mesh,ele))
    end do

    allocate(n_list(key_count(subdomain_mesh_node_list))) ! Nodal map from sub mesh --> full mesh
    allocate(inverse_n_list(node_count(external_mesh))) ! Nodal map from full mesh --> sub mesh

    ! Initialise mapping functions to zero. This is necessary for the inverse_node_list - if 
    ! after it is set up, the value in inverse_subnode_list = 0, this means that that element of 
    ! the full mesh does not have a corresponding element on the subdomain_mesh - i.e. it is not a part
    ! of the prognostic subdomain.

    n_list = 0
    inverse_n_list = 0

    ! Build maps to and from sub mesh:
    do i = 1, key_count(subdomain_mesh_node_list)
       n_list(i) = fetch(subdomain_mesh_node_list,i)
    end do

    do i = 1, key_count(subdomain_mesh_node_list)
       node = n_list(i)
       inverse_n_list(node) = i
    end do

    ! Allocate subdomain_mesh:
    shape => external_mesh%shape
    call allocate(subdomain_mesh, nodes=key_count(subdomain_mesh_node_list), elements=size(subele_list),&
         & shape=shape, name=trim(mesh_name))
    
    ! Determine ndglno (connectivity matrix) on subdomain_mesh:
    loc = shape%loc
    do i = 1, size(subele_list)
      ele = subele_list(i)
      call set_ele_nodes(subdomain_mesh, i, inverse_n_list(ele_nodes(external_mesh, ele)))
    end do

    ! Store subdomain_mesh element list, node list as mesh attributes:
    allocate(subdomain_mesh%subdomain_mesh)
    allocate(subdomain_mesh%subdomain_mesh%element_list(size(subele_list)))
    subdomain_mesh%subdomain_mesh%element_list = subele_list

    allocate(subdomain_mesh%subdomain_mesh%node_list(size(n_list)))
    subdomain_mesh%subdomain_mesh%node_list = n_list

    ! Calculate sndglno - an array of nodes corresponding to edges along surface:
    sloc = external_mesh%faces%shape%loc

    ! Begin by determining which faces are on subdomain_mesh boundaries:
    call allocate(face_list)
    do i = 1, size(subele_list)
      ele = subele_list(i)
      neigh => ele_neigh(external_mesh, ele) ! Determine element neighbours on parent mesh
      do ni = 1, size(neigh)
        ele_2 = neigh(ni)
        if (.not.has_value(subdomain_mesh_element_list, ele_2)) then
          face = ele_face(external_mesh, ele, ele_2)
          call insert(face_list, face)
        end if
      end do
    end do

    ! Allocate and initialise sndglno and boundary ids:
    edge_count = key_count(face_list)
    allocate(sndglno(edge_count*sloc))
    sndglno = 0
    allocate(boundary_ids(1:edge_count))
    boundary_ids = 0
    
    ! Set up sndglno and boundary_ids:
    do i = 1, edge_count
      face = fetch(face_list, i)
      sndglno((i-1)*sloc+1:i*sloc) = inverse_n_list(face_global_nodes(external_mesh, face))
      boundary_ids(i) = surface_element_id(external_mesh, face)
    end do

    call deallocate(face_list)

    ! Add faces to subdomain_mesh
    call add_faces(subdomain_mesh,sndgln=sndglno(1:edge_count*sloc), &
    &               boundary_ids=boundary_ids(1:edge_count))

    deallocate(sndglno)
    deallocate(boundary_ids)

    ! Insert mesh and position fields for subdomain_mesh into sub_states:

    if(mesh_name=="CoordinateMesh") then
      external_mesh_position => extract_vector_field(states(1), "Coordinate")
    else
      external_mesh_position => extract_vector_field(states(1), trim(mesh_name)//"Coordinate")
    end if

    call allocate(position, mesh_dim(subdomain_mesh), subdomain_mesh, trim(external_mesh_position%name))

    call remap_to_subdomain(external_mesh_position,position)

    ewrite(2,*) 'MinMax info for subdomain_mesh positions_field: '
    do i = 1, mesh_dim(subdomain_mesh)
       ewrite_minmax(position%val(i,:))
    end do

    ! Load sub_states:
    call insert(sub_states, subdomain_mesh, subdomain_mesh%name)
    call insert(sub_states, position, position%name)

    ! Clean up:
    call deallocate(subdomain_mesh)
    call deallocate(position)

    deallocate(n_list)
    deallocate(inverse_n_list)

    call deallocate(subdomain_mesh_node_list)

    deallocate(subele_list)

    call deallocate(subdomain_mesh_element_list)
    call deallocate(prescribed_region_id_set)

    ewrite(1,*) "Leaving derive external subdomain mesh"

  end subroutine derive_external_subdomain_mesh

  subroutine insert_subdomain_mesh_maps(states,sub_states)
    !! Forms mapping functions for subdomain derived meshes:
    type(state_type), intent(in), dimension(:) :: states
    type(state_type), intent(inout), dimension(:) :: sub_states

    ! Externally created mesh:
    type(mesh_type), pointer :: external_mesh
    ! Meshes on full domain / subdomain:
    type(mesh_type), pointer :: full_mesh, subdomain_mesh
    ! List of elements in subdomain:
    integer, dimension(:), pointer :: subele_list
    ! Integer array with this list of nodes:
    integer, dimension(:), allocatable :: node_list
    ! List of element nodes:
    integer, dimension(:), pointer :: nodesf, nodess

    ! Other declarations:
    integer :: nmeshes, imesh, ele, i, inode

    ewrite(1,*) "In insert_subdomain_mesh_maps"

    ! Extract sub mesh of external mesh from sub_state:
    external_mesh => get_external_mesh(sub_states)

    ! List of element on sub domain_mesh:
    subele_list => external_mesh%subdomain_mesh%element_list

    ! Get number of meshes
    nmeshes=mesh_count(states(1))

    ! Loop over meshes and derive list of nodes:
    do imesh = 1, nmeshes 

       full_mesh => extract_mesh(states(1),imesh)

       if(have_option(trim(full_mesh%option_path)// "/from_mesh")) then

          subdomain_mesh => extract_mesh(sub_states(1), imesh)

          ! Allocate subdomain_mesh attributes:
          allocate(subdomain_mesh%subdomain_mesh)
          allocate(subdomain_mesh%subdomain_mesh%element_list(size(subele_list)))

          ! Set subdomain_mesh's element list attribute:
          subdomain_mesh%subdomain_mesh%element_list = subele_list

          ! Allocate subdomain_mesh's node list and set up:
          allocate(subdomain_mesh%subdomain_mesh%node_list(node_count(subdomain_mesh)))
          allocate(node_list(node_count(subdomain_mesh)))

          do ele = 1, size(subele_list)
             nodess => ele_nodes(subdomain_mesh,ele) ! List of element nodes on subdomain_mesh
             nodesf => ele_nodes(full_mesh,subele_list(ele)) ! Corresponding list on global_mesh
             do inode = 1, size(nodess)
                node_list(nodess(inode)) = nodesf(inode)
             end do
          end do

          ! Store in subdomain mesh node_list attribute:
          subdomain_mesh%subdomain_mesh%node_list = node_list

          deallocate(node_list)

          ! Insert into all sub_states:
          do i = 2, size(sub_states)
            call insert(sub_states(i), subdomain_mesh, trim(subdomain_mesh%name))
          end do

       end if

    end do

    ewrite(1,*) "Leaving insert_subdomain_mesh_maps"

  end subroutine insert_subdomain_mesh_maps

  subroutine update_subdomain_fields(states,sub_states)

    ! This routine updates fields on the prognostic subdomain
    ! in partially prognostic simulations:

    type(state_type), intent(in), dimension(:) :: states
    type(state_type), intent(inout), dimension(:) :: sub_states

    ! Full domain and subdomain fields:
    type(scalar_field), pointer :: sfield, sfield_sub
    type(vector_field), pointer :: vfield, vfield_sub
    type(tensor_field), pointer :: tfield, tfield_sub
    
    ! Other declarations:
    integer :: nstates, nsfields, nvfields, ntfields
    integer :: istate, ifield, stat

    ewrite(1,*) "Entering update_subdomain_fields"    

    ! How many states exist?
    nstates = size(states)

    ! Loop over states:
    do istate = 1, nstates

       ! Loop over fields (scalar, vector, tensor) in order.
       ! Start with scalar fields:
       nsfields = scalar_field_count(sub_states(istate))
       do ifield = 1, nsfields
          ! Extract subdomain field from sub_state:
          sfield_sub => extract_scalar_field(sub_states(istate),ifield)
          if(.not. aliased(sfield_sub)) then
            ! Extract full domain field from State:
            sfield => extract_scalar_field(states(istate),trim(sfield_sub%name), stat=stat)
            if(stat==0) then
              ! Zero:
              call zero(sfield_sub)
              ! Then remap:
              call remap_to_subdomain(sfield,sfield_sub)
            end if
          end if
       end do

       ! Vector fields:
       nvfields = vector_field_count(sub_states(istate))
       do ifield = 1, nvfields
          vfield_sub => extract_vector_field(sub_states(istate),ifield)
          if(.not. aliased(vfield_sub)) then
            vfield => extract_vector_field(states(istate), trim(vfield_sub%name), stat=stat)
            if(stat==0) then
              call zero(vfield_sub)
              call remap_to_subdomain(vfield,vfield_sub)
            end if
          end if
       end do

       ! Tensor fields:
       ntfields = tensor_field_count(sub_states(istate))
       do ifield = 1, ntfields
          tfield_sub => extract_tensor_field(sub_states(istate),ifield)
          if(.not. aliased(tfield_sub)) then
            tfield => extract_tensor_field(states(istate),trim(tfield_sub%name), stat=stat)
            if(stat==0) then
              call zero(tfield_sub)
              call remap_to_subdomain(tfield,tfield_sub)
            end if
          end if
       end do

    end do

    ewrite(1,*) "Leaving update_subdomain_fields"    

  end subroutine update_subdomain_fields

  subroutine set_full_domain_prescribed_fields(states,time)
    !! Initialises prescribed fields in prescribed
    !! regions of domain (i.e. in full state):

    !! Note currently only set up for vector field velocity

    type(state_type), dimension(:), intent(in):: states

    type(vector_field), pointer :: vfield
    type(vector_field), pointer :: position
    !! current time if not using that in the options tree
    real, intent(in), optional :: time

    integer :: istate, ifield, nstates, nvfields

    ewrite(1,*) 'Setting full domain prescribed fields'

    ! Determine number of states:
    nstates = size(states)

    ! Loop over states:
    do istate = 1, nstates
       
       ! Deal with vector fields:
       nvfields = vector_field_count(states(istate))

       do ifield = 1, nvfields

          vfield => extract_vector_field(states(istate), ifield)

          ! At present this only works for velocity:
          if (have_option(trim(vfield%option_path)// "/prognostic/prescribed_region") .and. &
               .not. aliased(vfield) ) then

             position => get_external_coordinate_field(states(istate), vfield%mesh)

             call initialise_field_over_regions(vfield, &
                trim(vfield%option_path)// "/prognostic/prescribed_region", &
                position,time=time)

          end if
       end do

    end do

    ewrite(1,*) 'Finished setting full domain prescribed fields'

  end subroutine set_full_domain_prescribed_fields

  subroutine sub_state_remap_to_full_mesh(states, sub_states)

    ! Remap fields from subdomain to full mesh so that
    ! full mesh prognostic fields that depend on submesh
    ! variables can be solved.

    type(state_type), dimension(:), intent(in):: states
    type(state_type), dimension(:), intent(in):: sub_states

    type(scalar_field), pointer :: sfield, sfield_sub
    type(vector_field), pointer :: vfield, vfield_sub
    type(tensor_field), pointer :: tfield, tfield_sub

    integer :: istate, ifield, nstates, nsfields, nvfields, ntfields, stat

    ! Determine number of states:
    nstates = size(states)
    assert(size(states)==size(sub_states))

    ! Loop over states:
    do istate = 1, nstates

       ! Deal with scalar fields:
       nsfields = scalar_field_count(states(istate))

       do ifield = 1, nsfields

          sfield => extract_scalar_field(states(istate),ifield)
          if(.not. aliased(sfield)) then
            sfield_sub => extract_scalar_field(sub_states(istate), trim(sfield%name), stat=stat)
            if(stat==0) then
              call remap_to_full_domain(sfield_sub,sfield)
            end if
          end if

       end do

       ! Deal with vector fields:
       nvfields = vector_field_count(states(istate))

       do ifield = 1, nvfields

          vfield => extract_vector_field(states(istate), ifield)
          if(.not. aliased(vfield)) then
            vfield_sub => extract_vector_field(sub_states(istate), trim(vfield%name), stat=stat)
            if(stat==0) then
              call remap_to_full_domain(vfield_sub,vfield)
            end if
          end if

       end do

       ! Deal with tensor fields:
       do ifield = 1, ntfields

          tfield => extract_tensor_field(states(istate),ifield)
          if(.not. aliased(tfield)) then
            tfield_sub => extract_tensor_field(sub_states(istate), trim(tfield%name), stat=stat)
            if(stat==0) then
              call remap_to_full_domain(tfield_sub,tfield)
            end if
          end if

       end do

    end do

  end subroutine sub_state_remap_to_full_mesh

end module populate_sub_state_module
