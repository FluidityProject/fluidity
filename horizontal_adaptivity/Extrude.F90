#include "fdebug.h"

module hadapt_extrude
  !!< Extrude a given 2D mesh to a full 3D mesh.
  !!< The layer depths are specified by a sizing function
  !!< which can be arbitrary python.
  use elements
  use fields
  use spud
  use quadrature
  use global_parameters
  use sparse_tools
  use hadapt_advancing_front
  use hadapt_metric_based_extrude
  implicit none

  private
  
  public :: extrude, compute_z_nodes

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
    character(len=PYTHON_FUNC_LEN) :: sizing_function
    logical:: sizing_is_constant
    real:: constant_sizing, depth
    integer:: stat, h_dim, column, quadrature_degree

    !! We assume that for the extrusion operation,
    !! the layer configuration is independent of x and y.
    !! (i.e. the layers are the same everywhere.)

    !! Checking linearity of h_mesh.
    assert(h_mesh%mesh%shape%degree == 1)
    assert(h_mesh%mesh%continuity >= 0)

    call add_nelist(h_mesh%mesh)

    ! get the extrusion options
    call get_option(trim(option_path)//'/from_mesh/extrude/bottom_depth', &
       depth)
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
      constant_value=-1
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
      delta_h = get_delta_h(node_val(z_mesh, node-1), is_constant, constant_value, py_func)
      call set(z_mesh, node, (/node_val(z_mesh, node-1) - delta_h/))
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
