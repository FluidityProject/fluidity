#include "fdebug.h"

module hadapt_extrude
  !!< Extrude a given 2D mesh to a full 3D mesh.
  !!< The layer depths are specified by a sizing function
  !!< which can be arbitrary python.
! these 5 need to be on top and in this order, so as not to confuse silly old intel compiler 
  use quadrature
  use elements
  use sparse_tools
  use fields
  use state_module
!
  use spud
  use global_parameters
  use hadapt_advancing_front
  use hadapt_metric_based_extrude
  implicit none

  private
  
  public :: extrude, extrude_azimuthal, compute_z_nodes

  contains

  subroutine extrude(h_mesh, option_path, out_mesh)
    !!< The horizontal 2D mesh.
    !!< Note: this must be linear.
    type(vector_field), intent(inout) :: h_mesh
    !!< options to be set for out_mesh,
    !!< at the moment: /name, and under from_mesh/extrude/:
    !!< depth, sizing_function optionally top_surface_id and bottom_surface_id
    character(len=*), intent(in) :: option_path
    !!< The full extruded 3D mesh.
    type(vector_field), intent(out) :: out_mesh

    character(len=FIELD_NAME_LEN):: mesh_name    
    type(quadrature_type) :: quad
    type(element_type) :: full_shape
    type(vector_field), dimension(node_count(h_mesh)) :: z_meshes
    character(len=PYTHON_FUNC_LEN) :: sizing_function, depth_function
    logical:: sizing_is_constant, depth_is_constant
    real:: constant_sizing, depth
    integer:: stat, h_dim, column, quadrature_degree
    
    real, dimension(1) :: tmp_depth
    real, dimension(mesh_dim(h_mesh), 1) :: tmp_pos

    !! We assume that for the extrusion operation,
    !! the layer configuration is independent of x and y.
    !! (i.e. the layers are the same everywhere.)

    !! Checking linearity of h_mesh.
    assert(h_mesh%mesh%shape%degree == 1)
    assert(h_mesh%mesh%continuity >= 0)

    call add_nelist(h_mesh%mesh)

    ! get the extrusion options
    call get_option(trim(option_path)//'/from_mesh/extrude/bottom_depth/constant', &
       depth, stat=stat)
    if (stat==0) then
       depth_is_constant = .true.
    else
       depth_is_constant = .false.
       call get_option(trim(option_path)//'/from_mesh/extrude/bottom_depth/python', &
          depth_function, stat=stat)
       if (stat /= 0) then
         FLAbort("Unknown way of specifying bottom depth function in mesh extrusion")
       end if
    end if
    
    call get_option(trim(option_path)//'/from_mesh/extrude/sizing_function/constant', &
       constant_sizing, stat=stat)
    if (stat==0) then
       sizing_is_constant=.true.
    else
       sizing_is_constant=.false.
       call get_option(trim(option_path)//'/from_mesh/extrude/sizing_function/python', &
          sizing_function, stat=stat)
       if (stat/=0) then
         FLAbort("Unknown way of specifying sizing function in mesh extrusion")
       end if       
    end if
       
    ! create a 1d vertical mesh under each surface node
    do column=1, size(z_meshes)
      
      if(.not.depth_is_constant) then
        tmp_pos(:,1) = node_val(h_mesh, column)
        call set_from_python_function(tmp_depth, trim(depth_function), tmp_pos, time=0.0)
        depth = tmp_depth(1)
      end if
      
      if (sizing_is_constant) then
        call compute_z_nodes(z_meshes(column), depth, node_val(h_mesh, column), &
          sizing=constant_sizing)
      else
        call compute_z_nodes(z_meshes(column), depth, node_val(h_mesh, column), &
          sizing_function=sizing_function)
      end if
    end do
      
    ! Now the tiresome business of making a shape function.
    h_dim = mesh_dim(h_mesh)
    call get_option("/geometry/quadrature/degree", quadrature_degree)
    quad = make_quadrature(vertices=h_dim + 2, dim=h_dim + 1, degree=quadrature_degree)
    full_shape = make_element_shape(vertices=h_dim + 2, dim=h_dim + 1, degree=1, quad=quad)
    call deallocate(quad)

    call get_option(trim(option_path)//'/name', mesh_name)
      
    ! combine the 1d vertical meshes into a full mesh
    call combine_z_meshes(h_mesh, z_meshes, out_mesh, &
       full_shape, mesh_name)
       
    do column=1, node_count(h_mesh)
      call deallocate(z_meshes(column))
    end do
    call deallocate(full_shape)

    out_mesh%mesh%option_path=option_path
    out_mesh%option_path=""

    if (has_faces(h_mesh%mesh)) then
      ! if the horizontal mesh has a surface mesh
      ! give one to the extruded mesh as well
      call generate_extruded_surface_mesh(h_mesh, out_mesh)
    end if
        
  end subroutine extrude
  
  subroutine extrude_azimuthal(v_mesh, out_mesh, ndivisions)
    !!< Extrude the given 2D around the z-axis mesh to form a 3D annular domain
    
    type(vector_field), intent(in) :: v_mesh
    type(vector_field), target, intent(out) :: out_mesh
    integer, intent(in) :: ndivisions
    
    ! Linear simplices only
    integer, parameter :: dim = 3, sdim = 2, nloc = 4, snloc = 3
    
    integer :: i, index, j, neles_2d, neles_3d, nnodes_2d, nnodes_3d
    integer, dimension(nloc) :: nodes_3d
    integer, dimension(node_count(v_mesh), ndivisions + 1) :: chain_map
    integer, dimension(:), pointer :: nodes_2d
    real :: r, x, y, z
    real, dimension(ndivisions) :: theta
    type(element_type) :: shape_3d
    type(element_type), pointer :: shape_2d
    type(mesh_type), pointer :: mesh_3d
    type(quadrature_type) :: quad_3d
    
    assert(.not. isparallel())
    assert(v_mesh%dim == sdim)
    assert(mesh_dim(v_mesh) == sdim)
    assert(ndivisions > 0)
    
    shape_2d => ele_shape(v_mesh, 1)
    assert(shape_2d%degree == 1)
    assert(ele_numbering_family(shape_2d) == FAMILY_SIMPLEX)
    
    ! Generate the 3d shape function
    quad_3d = make_quadrature(vertices = nloc, dim = dim, degree = shape_2d%quadrature%degree)
    shape_3d = make_element_shape(vertices = nloc, dim = dim, degree = 1, quad = quad_3d)
    call deallocate(quad_3d)
    
    nnodes_2d = node_count(v_mesh)
    nnodes_3d = nnodes_2d * ndivisions  ! Fence posts = fence sections (periodic)
    
    neles_2d = ele_count(v_mesh)
    neles_3d = neles_2d * ndivisions * 3  ! Split each triangular prism into 3 tets
    
    allocate(mesh_3d)
    call allocate(mesh_3d, nnodes_3d, neles_3d, shape_3d, name = trim(v_mesh%mesh%name) // "AzimuthalExtrusion")
    call deallocate(shape_3d)
    
    call allocate(out_mesh, dim, mesh_3d, "Coordinate")
    deallocate(mesh_3d)
    mesh_3d => out_mesh%mesh
    
    ! Generate the nodes
    index = 0
    do i = 1, ndivisions
      theta(i) = 2.0 * pi * float(i - 1) / float(ndivisions)
    end do
    do i = 1, nnodes_2d
      r = node_val(v_mesh, 1, i)
      z = node_val(v_mesh, 2, i)
      do j = 1, ndivisions
        index = index + 1
        x = r * cos(theta(j))
        y = r * sin(theta(j))
        call set(out_mesh, index, (/x, y, z/))
        chain_map(i, j) = index
      end do
    end do
    assert(index == nnodes_3d)
    chain_map(:, ndivisions + 1) = chain_map(:, 1)
    
    ! Generate the element node list
    index = 0
    do i = 1, neles_2d
      nodes_2d => ele_nodes(v_mesh, i)
      assert(size(nodes_2d) == snloc)
      do j = 1, ndivisions
        index = index + 1
        nodes_3d = (/chain_map(nodes_2d(1), j), chain_map(nodes_2d(1), j + 1), chain_map(nodes_2d(3), j), chain_map(nodes_2d(2), j)/)
        mesh_3d%ndglno((index - 1) * nloc + 1:index * nloc) = nodes_3d
        index = index + 1
        nodes_3d = (/chain_map(nodes_2d(3), j), chain_map(nodes_2d(1), j + 1), chain_map(nodes_2d(3), j + 1), chain_map(nodes_2d(2), j)/)
        mesh_3d%ndglno((index - 1) * nloc + 1:index * nloc) = nodes_3d
        index = index + 1
        nodes_3d = (/chain_map(nodes_2d(1), j + 1), chain_map(nodes_2d(2), j), chain_map(nodes_2d(2), j + 1), chain_map(nodes_2d(3), j + 1)/)
        mesh_3d%ndglno((index - 1) * nloc + 1:index * nloc) = nodes_3d
      end do
    end do
    assert(index == neles_3d)
    
    ! Surface mesh not currently generated
    call add_faces(mesh_3d)
    
    mesh_3d%option_path = empty_path
    out_mesh%option_path = empty_path
    
    call deallocate(mesh_3d)
    
  end subroutine extrude_azimuthal

  subroutine compute_z_nodes(z_mesh, depth, xy, sizing, sizing_function)
    !!< Figure out at what depths to put the layers.
    type(vector_field), intent(out) :: z_mesh
    real, intent(in):: depth
    real, dimension(:), intent(in):: xy
    real, optional, intent(in):: sizing
    character(len=*), optional, intent(in):: sizing_function

    ! this is a safety gap, for people doing something stupid:
    integer, parameter:: MAX_VERTICAL_NODES=1e6
    ! to prevent infinitesimally thin bottom layer if sizing function
    ! is an integer mulitple of total depth, the bottom layer needs
    ! to have at least this fraction of the layer depth above it:
    real, parameter:: MIN_BOTTOM_LAYER_FRAC=1e-3
    
    integer :: elements
    logical :: is_constant
    real :: constant_value

    type(mesh_type) :: mesh
    type(element_type) :: oned_shape
    type(quadrature_type) :: oned_quad
    integer :: quadrature_degree
    integer :: ele
    integer, parameter :: loc=2
    integer :: node
    real, dimension(1:size(xy)+1):: xyz
    real :: delta_h, z
    character(len=PYTHON_FUNC_LEN) :: py_func

    call get_option("/geometry/quadrature/degree", quadrature_degree)
    oned_quad = make_quadrature(vertices=loc, dim=1, degree=quadrature_degree)
    oned_shape = make_element_shape(vertices=loc, dim=1, degree=1, quad=oned_quad)
    call deallocate(oned_quad)

    if (present(sizing)) then
      is_constant=.true.
      constant_value=sizing
      py_func = " "
    else if (present(sizing_function)) then
      is_constant=.false.
      constant_value=-1.0
      py_func = sizing_function
    else
      FLAbort("Need to supply either sizing or sizing_function")
    end if

    ! First work out number of nodes/elements:
    z=0.0
    node=2
    xyz(1:size(xy))=xy
    do
      xyz(size(xy)+1)=z
      delta_h = get_delta_h( xyz, is_constant, constant_value, py_func)
      z=z - delta_h
      if (z<-depth+MIN_BOTTOM_LAYER_FRAC*delta_h) exit
      node=node+1
      if (node>MAX_VERTICAL_NODES) then
        ewrite(-1,*) "Check your extrude/sizing_function"
        FLAbort("Maximum number of vertical layers reached")
      end if
    end do
    elements=node-1

    call allocate(mesh, elements+1, elements, oned_shape, "ZMesh")
    call deallocate(oned_shape)
    do ele=1,elements
      mesh%ndglno((ele-1) * loc + 1: ele*loc) = (/ele, ele+1/)
    end do

    call allocate(z_mesh, 1, mesh, "ZMeshCoordinates")
    call deallocate(mesh)

    ! Start the mesh at z=0 and work down to z=-depth.
    ! Someone else can refine this later.
    ! Have fun!
    call set(z_mesh, 1, (/0.0/))
    do node=2,elements
      xyz(size(xy)+1:)=node_val(z_mesh, node-1)
      delta_h = get_delta_h(xyz, is_constant, constant_value, py_func)
      call set(z_mesh, node,  (/ xyz(size(xy)+1) - delta_h /) )
    end do
    call set(z_mesh, elements+1, (/-depth/))

    ! For pathological sizing functions the mesh might have gotten inverted at the last step.
    ! If you encounter this, make this logic smarter.
    assert(all(node_val(z_mesh, elements) > node_val(z_mesh, elements+1)))
    
    assert(oned_quad%refcount%count == 1)
    assert(oned_shape%refcount%count == 1)
    assert(z_mesh%refcount%count == 1)
    assert(mesh%refcount%count == 1)

    contains
    
      function get_delta_h(pos, is_constant, constant_value, py_func) result(delta_h)
        real, dimension(:), intent(in) :: pos
        logical, intent(in) :: is_constant
        real, intent(in) :: constant_value
        character(len=PYTHON_FUNC_LEN), intent(in) :: py_func

        real :: delta_h
        real, dimension(1) :: delta_h_tmp
        real, dimension(size(pos), 1) :: pos_tmp
        
        if (is_constant) then
          delta_h = constant_value
        else
          pos_tmp(:, 1) = pos
          call set_from_python_function(delta_h_tmp, trim(py_func), pos_tmp, time=0.0)
          delta_h = delta_h_tmp(1)
        end if
        assert(delta_h > 0.0)
        
      end function get_delta_h
      
  end subroutine compute_z_nodes

    
end module hadapt_extrude
