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
module fields_data_types

  use global_parameters, only:FIELD_NAME_LEN, current_debug_level, OPTION_PATH_LEN, PYTHON_FUNC_LEN
  use picker_data_types
  use shape_functions
  use sparse_tools
  use spud
  use reference_counting
  use halo_data_types
  use data_structures, only : integer_set_vector
  implicit none

  private
  public adjacency_cache, &
     mesh_type, mesh_faces, mesh_subdomain_mesh, scalar_field, vector_field, tensor_field, &
     mesh_pointer, scalar_field_pointer, vector_field_pointer, tensor_field_pointer, &
     scalar_boundary_condition, vector_boundary_condition, &
     scalar_boundary_conditions_ptr, vector_boundary_conditions_ptr

  !! Types of different halo associated with a field:
  integer, public, parameter :: HALO_TYPES=2
  !! Available sources of data for fields:
  integer, public, parameter :: FIELD_TYPE_NORMAL=0, FIELD_TYPE_CONSTANT=1, FIELD_TYPE_PYTHON=2, &
    FIELD_TYPE_DEFERRED=3

  type adjacency_cache
    type(csr_sparsity), pointer :: nnlist => null()
    type(csr_sparsity), pointer :: nelist => null()
    type(csr_sparsity), pointer :: eelist => null()
  end type adjacency_cache

  type mesh_type
     !!< Mesh information for (among other things) fields.
     integer, dimension(:), pointer :: ndglno
     !! Flag for whether ndglno is allocated
     logical :: wrapped=.true.
     type(element_type) :: shape
     integer :: elements
     integer :: nodes
     character(len=FIELD_NAME_LEN) :: name
     !! path to options in the options tree
#ifdef DDEBUG
     character(len=OPTION_PATH_LEN) :: option_path="/uninitialised_path/"
#else
     character(len=OPTION_PATH_LEN) :: option_path
#endif
     !! Degree of continuity of the field. 0 is for the conventional C0
     !! discretisation. -1 for DG.
     integer :: continuity=0
     !! Reference count for mesh
     type(refcount_type), pointer :: refcount=>null()
     !! Mesh face information for those meshes (eg discontinuous) which need it.
     type(mesh_faces), pointer :: faces=>null()
     !! Information on subdomain_ mesh, for partially prognostic solves:
     type(mesh_subdomain_mesh), pointer :: subdomain_mesh=>null()
     type(adjacency_cache), pointer :: adj_lists => null()
     !! array that for each node tells which column it is in
     !! (column numbers usually correspond to a node number in a surface mesh)
     integer, dimension(:), pointer :: columns => null()
     !! if this mesh is extruded this array says which horizontal mesh element each element is below
     integer, dimension(:), pointer :: element_columns => null()
     !! A list of ids marking different parts of the mesh
     !! so that initial conditions can be associated with it.
     integer, dimension(:), pointer :: region_ids=>null()
     !! Halo information for parallel simulations.
     type(halo_type), dimension(:), pointer :: halos=>null()
     type(halo_type), dimension(:), pointer :: element_halos=>null()
     type(integer_set_vector), dimension(:), pointer :: colourings=>null()
     !! A logical indicating if this mesh is periodic or not
     !! (does not tell you how periodic it is... i.e. true if
     !! any surface is periodic)
     logical :: periodic=.false.
  end type mesh_type

  type mesh_faces
     !!< Type encoding face information for a mesh.
     type(element_type), pointer :: shape
     !! Face_list(i,j) is the face in element i bordering j.
     type(csr_matrix) :: face_list
     !! The local number of the nodes in a given face.
     integer, dimension(:), pointer :: face_lno
     !! A mesh consisting of all faces on the surface of the domain,
     !! it uses its own internal surface node numbering:
     type(mesh_type) surface_mesh
     !! A list of the nodes on the surface, thus forming a map between
     !! internal surface node numbering and global node numbering:
     integer, dimension(:), pointer :: surface_node_list
     !! The element with which each face is associated
     integer, dimension(:), pointer :: face_element_list
     !! A list of ids marking different parts of the surface mesh
     !! so that boundary conditions can be associated with it.
     integer, dimension(:), pointer :: boundary_ids
     !! list of ids to identify coplanar patches of the surface:
     integer, dimension(:), pointer :: coplanar_ids => null()
     !! a DG version of the surface mesh, useful for storing bc values
     type(mesh_type), pointer:: dg_surface_mesh => null()
     !! A logical indicating if this mesh has a discontinuous internal boundary
     !! This means that the pairs of internal facets are allowed to have two different
     !! surface ids. When writing out this mesh both facets are written out and
     !! element owners (an extra column indiciating which element is adjacent to each facet)
     !! needs to be written out along with the surface mesh.
     !! This is currently only used for periodic meshes which have a physical and aliased surface id
     !! along the periodic boundary (which is an internal boundary of the periodic mesh).
     !! Note that other meshes (with has_internal_boundaries==.false.) may still have internal facets
     !! as part of the surface mesh, in this case the surface ids do have to agree and only one of each
     !! pair of facets is written when writing out the mesh.
     logical :: has_discontinuous_internal_boundaries=.false.
     !! If internal facets are present in the surface mesh (and has_discontinuous_internal_boundaries==.false.) 
     !! the surface facets are numbered such that 1:unique_surface_element_count visits each external facet
     !! and each pair of internal facets only once (the order of this is typically determined by the 
     !! read-in input mesh in which internal facets are only also only present once)
     !! The second facets of each pair of internal facets are numbered 
     !! unique_surface_element_count+1:surface_element_count (surface_element_count()==size(boundary_ids))
     !! For meshes with no internal facets: unique_surface_element_count==surface_element_count()
     !! For meshes with has_discontinuous_internal_boundaries no order is guaranteed and also
     !! unique_surface_element_count==surface_element_count()
     integer :: unique_surface_element_count
  end type mesh_faces

  type mesh_subdomain_mesh
     !! List of elements in subdomain_mesh region(s):
     integer, dimension(:), pointer :: element_list
     !! List of nodes in subdomain_mesh regions(s):
     integer, dimension(:), pointer :: node_list
  end type mesh_subdomain_mesh

  type scalar_field 
     !! Field value at points.
     real, dimension(:), pointer :: val
     !! Stride of val
     integer :: val_stride = 1
     !! Flag for whether val is allocated
     logical :: wrapped=.true.
     !! The data source to be used
     integer :: field_type = FIELD_TYPE_NORMAL
     !! boundary conditions:
     type(scalar_boundary_conditions_ptr), pointer :: bc => null()
     character(len=FIELD_NAME_LEN) :: name
     !! path to options in the options tree
#ifdef DDEBUG
     character(len=OPTION_PATH_LEN) :: option_path="/uninitialised_path/"
#else
     character(len=OPTION_PATH_LEN) :: option_path
#endif
     type(mesh_type) :: mesh
     !! Reference count for field
     type(refcount_type), pointer :: refcount=>null()
     !! Indicator for whether this is an alias to another field.
     logical :: aliased=.false.
     !! Python-field implementation.
     real, dimension(:, :), pointer :: py_locweight => null()
     character(len=PYTHON_FUNC_LEN) :: py_func
     type(vector_field), pointer :: py_positions
     logical :: py_positions_same_mesh
     integer :: py_dim
     type(element_type), pointer :: py_positions_shape => null()
  end type scalar_field

  type vector_field
     !! dim x nonods vector values
     real, dimension(:,:), pointer :: val
     !! Flag for whether val is allocated
     logical :: wrapped = .true.
     !! The data source to be used
     integer :: field_type = FIELD_TYPE_NORMAL
     !! boundary conditions:
     type(vector_boundary_conditions_ptr), pointer :: bc => null()
     character(len=FIELD_NAME_LEN) :: name
     integer :: dim
     !! path to options in the options tree
#ifdef DDEBUG
     character(len=OPTION_PATH_LEN) :: option_path="/uninitialised_path/"
#else
     character(len=OPTION_PATH_LEN) :: option_path
#endif
     type(mesh_type) :: mesh
     !! Reference count for field
     type(refcount_type), pointer :: refcount=>null()
     !! Indicator for whether this is an alias to another field.
     logical :: aliased=.false.
     !! Picker used for spatial indexing (pointer to a pointer to ensure
     !! correct handling on assignment)
     type(picker_ptr), pointer :: picker => null()
  end type vector_field

  type tensor_field
     !! ndim x ndim x nonods
     real, dimension(:,:,:), pointer :: val
     !! Flag for whether val is allocated
     logical :: wrapped=.true.
     !! The data source to be used
     integer :: field_type = FIELD_TYPE_NORMAL
     character(len=FIELD_NAME_LEN) :: name
     integer, dimension(2) :: dim
     !! path to options in the options tree
#ifdef DDEBUG
     character(len=OPTION_PATH_LEN) :: option_path="/uninitialised_path/"
#else
     character(len=OPTION_PATH_LEN) :: option_path
#endif
     type(mesh_type) :: mesh
     !! Reference count for field
     type(refcount_type), pointer :: refcount=>null()
     !! Indicator for whether this is an alias to another field.
     logical :: aliased=.false.
  end type tensor_field

  type mesh_pointer
     !!< Dummy type to allow for arrays of pointers to meshes
     type(mesh_type), pointer :: ptr => null()
  end type mesh_pointer

  type scalar_field_pointer
     !!< Dummy type to allow for arrays of pointers to scalar fields
     type(scalar_field), pointer :: ptr => null()
  end type scalar_field_pointer

  type vector_field_pointer
     !!< Dummy type to allow for arrays of pointers to vector fields
     type(vector_field), pointer :: ptr => null()
  end type vector_field_pointer

  type tensor_field_pointer
     !!< Dummy type to allow for arrays of pointers to tensor fields
     type(tensor_field), pointer :: ptr => null()
  end type tensor_field_pointer

  type scalar_boundary_condition
     !!< Type to hold boundary condition information for a scalar field
     character(len=FIELD_NAME_LEN) :: name
     !! b.c. type, any of: ...
     character(len=FIELD_NAME_LEN) :: type=""
     !! list of surface elements to which boundary condition is applied:
     integer, dimension(:), pointer:: surface_element_list => null()
     !! list of surface nodes to which boundary condition is applied:
     integer, dimension(:), pointer:: surface_node_list => null()
     !! mesh consisting of these elements and nodes only:
     type(mesh_type), pointer :: surface_mesh
     !! surface fields on this mesh containing b.c. values
     type(scalar_field), dimension(:), pointer :: surface_fields => null()
     !! path to options in the options tree
#ifdef DDEBUG
     character(len=OPTION_PATH_LEN) :: option_path="/uninitialised_path/"
#else
     character(len=OPTION_PATH_LEN) :: option_path
#endif
  end type scalar_boundary_condition

  type vector_boundary_condition
     !!< Type to hold boundary condition information for a vector field
     character(len=FIELD_NAME_LEN) :: name
     !! b.c. type, any of: ...
     character(len=FIELD_NAME_LEN) :: type=""
     !! boundary condition is only applied for component with applies is .true.
     logical, dimension(3):: applies
     !! list of surface elements to which boundary condition is applied:
     integer, dimension(:), pointer:: surface_element_list => null()
     !! list of surface nodes to which boundary condition is applied:
     integer, dimension(:), pointer:: surface_node_list => null()
     !! mesh consisting of these elements and nodes only:
     type(mesh_type), pointer :: surface_mesh
     !! surface fields on this mesh containing b.c. values
     type(vector_field), dimension(:), pointer :: surface_fields => null()
     !! scalar surface fields on this mesh containing b.c. values
     type(scalar_field), dimension(:), pointer :: scalar_surface_fields => null()
     !! path to options in the options tree
#ifdef DDEBUG
     character(len=OPTION_PATH_LEN) :: option_path="/uninitialised_path/"
#else
     character(len=OPTION_PATH_LEN) :: option_path
#endif
  end type vector_boundary_condition
    
  ! container for pointer to array of scalar bcs. This is put in a separate
  ! container pointed at by the field so that we don't leak memory
  ! if we change the bcs on one copy of the field when there are more copies around
  type scalar_boundary_conditions_ptr
     type(scalar_boundary_condition), dimension(:), pointer:: boundary_condition => null()
  end type scalar_boundary_conditions_ptr
  
  ! container for pointer to array of vector bcs. This is put in a separate
  ! container pointed at by the field so that we don't leak memory
  ! if we change the bcs on one copy of the field when there are more copies around
  type vector_boundary_conditions_ptr
     type(vector_boundary_condition), dimension(:), pointer:: boundary_condition => null()
  end type vector_boundary_conditions_ptr
  
end module fields_data_types
