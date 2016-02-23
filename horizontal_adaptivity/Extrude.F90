#include "fdebug.h"

module hadapt_extrude
  !!< Extrude a given 2D mesh to a full 3D mesh.
  !!< The layer depths are specified by a sizing function
  !!< which can be arbitrary python.

  use fldebug
  use global_parameters
  use futils, only: int2str
  use quadrature
  use elements
  use spud
  use quadrature
  use global_parameters
  use data_structures
  use parallel_tools
  use sparse_tools
  use linked_lists
  use parallel_fields
  use fields
  use vtk_interfaces
  use halos
  use hadapt_combine_meshes

  implicit none

  private
  
  public :: extrude, compute_z_nodes, hadapt_extrude_check_options, populate_depth_vector

  interface compute_z_nodes
    module procedure compute_z_nodes_wrapper, compute_z_nodes_sizing
  end interface compute_z_nodes

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

    type(integer_set), dimension(:), allocatable :: region_columns
    character(len=OPTION_PATH_LEN):: region_option_path, layer_option_path
    character(len=FIELD_NAME_LEN):: mesh_name, file_name  
    type(quadrature_type) :: quad
    type(element_type) :: full_shape
    type(vector_field) :: constant_z_mesh
    type(vector_field), dimension(:,:), allocatable :: z_meshes
    character(len=PYTHON_FUNC_LEN) :: sizing_function, depth_function
    real, dimension(:), allocatable :: sizing_vector
    logical:: depth_from_python, depth_from_map, have_min_depth, radial_extrusion
    real, dimension(:), allocatable :: depth_vector, top_depth
    real:: min_depth, surface_height
    logical:: sizing_is_constant, depth_is_constant, varies_only_in_depth, list_sizing
    real:: constant_sizing, depth, min_bottom_layer_frac
    integer:: h_dim, column, quadrature_degree

    logical, dimension(:), allocatable :: sigma_layers
    integer :: number_sigma_layers, nlayers, nregions
    
    integer :: i, ele, r, visit_count, layer
    integer, dimension(:), pointer :: nodes
    logical, dimension(:), allocatable :: column_visited
    integer, dimension(:), allocatable :: region_ids
    integer, dimension(2) :: shape_option
    logical :: apply_region_ids, constant_z_mesh_initialised

    !! Checking linearity of h_mesh.
    assert(h_mesh%mesh%shape%degree == 1)
    assert(h_mesh%mesh%continuity >= 0)


    call add_nelist(h_mesh%mesh)
    
    nregions = option_count(trim(option_path)//'/from_mesh/extrude/regions')
    if(nregions==0) then
      ewrite(-1,*) "I've been told to extrude but have found no regions options."
      FLExit("No regions options found under extrude.")
    elseif(nregions<0) then
      FLAbort("Negative number of regions options found under extrude.")
    end if

    nlayers = option_count(trim(option_path)//"/from_mesh/extrude/regions[0]/layers") ! must be the same in each region

    apply_region_ids = (nregions>1)
    allocate(region_columns(nregions))
    
    if (apply_region_ids) then

      allocate(column_visited(1:node_count(h_mesh)))
      column_visited = .false.
      visit_count = 0
      ! work out which nodes are associated with each region spec.
      ! loop backwards so the last region-spec. wins for shared nodes
      do r = nregions, 1, -1
        region_option_path = trim(option_path)//"/from_mesh/extrude/regions["//int2str(r-1)//"]"
        if (option_count(trim(region_option_path)//"/layers")/=nlayers) then
          FLExit("With extrusion the number of layers within each region needs to be the same")
        end if
        call allocate(region_columns(r))
        shape_option=option_shape(trim(region_option_path)//"/region_ids")
        allocate(region_ids(1:shape_option(1)))
        call get_option(trim(region_option_path)//"/region_ids", region_ids)
        do ele = 1, element_count(h_mesh)
          if (any(region_ids==ele_region_id(h_mesh, ele))) then
            nodes => ele_nodes(h_mesh, ele)
            do i=1, size(nodes)
              column = nodes(i)
              if (node_owned(h_mesh, column) .and. .not. column_visited(column)) then
                call insert(region_columns(r), column)
                column_visited(column) = .true.
                visit_count = visit_count+1
              end if
            end do
          end if
        end do
        deallocate(region_ids)
      end do
      if (nowned_nodes(h_mesh)/=visit_count) then
        FLExit("Not all parts of the horizontal mesh have extruded mesh regions associated with them.")
      end if
      deallocate(column_visited)
    else
      call allocate(region_columns(1))
      do column=1, node_count(h_mesh)
        if (node_owned(h_mesh, column)) then
          call insert(region_columns(1), column)
        end if
      end do
    end if
    allocate(z_meshes(nlayers, node_count(h_mesh)))

    ! depth of top surface (==0.0 for top layer, but increases for layers below)
    allocate(top_depth(node_count(h_mesh)))
    top_depth = 0.0

    ! auxillary array for depth_from_map:
    allocate(depth_vector(node_count(h_mesh)))

    ! may vary per layer, but not per region
    allocate(sigma_layers(nlayers))

    radial_extrusion = have_option("/geometry/spherical_earth")
    
    regions: do r = 1, nregions

      if (key_count(region_columns(r))==0) cycle
      region_option_path = trim(option_path)//"/from_mesh/extrude/regions["//int2str(r-1)//"]"
      
      layers: do layer = 1, nlayers
        layer_option_path = trim(region_option_path)//"/layers["//int2str(layer-1)//"]"
        call get_layer_extrusion_options(layer_option_path, &
                   depth_is_constant, depth, depth_from_python, depth_function, depth_from_map, &
                   file_name, have_min_depth, min_depth, surface_height, sizing_is_constant, constant_sizing, list_sizing, &
                   sizing_function, sizing_vector, min_bottom_layer_frac, varies_only_in_depth, sigma_layers(layer), number_sigma_layers, &
                   radial_extrusion)

        if (depth_from_map) then
          call populate_depth_vector(h_mesh,file_name,depth_vector,surface_height,radial_extrusion)
        end if
        
        if(varies_only_in_depth .and. depth_is_constant) then
          column = fetch(region_columns(r), 1)
          call compute_z_nodes(constant_z_mesh, node_val(h_mesh, column), min_bottom_layer_frac, &
                          top_depth(column), &
                          depth_is_constant, depth, depth_from_python, depth_function, &
                          depth_from_map, depth_vector(column),  have_min_depth, min_depth, &
                          sizing_is_constant, constant_sizing, list_sizing, sizing_function, sizing_vector, &
                          sigma_layers(layer), number_sigma_layers, radial_extrusion)
          do i=1, key_count(region_columns(r))
            column = fetch(region_columns(r), i)
            call get_previous_z_nodes(z_meshes(layer, column), constant_z_mesh)
          end do
          call deallocate(constant_z_mesh)
        else
          do i=1, key_count(region_columns(r))
            column = fetch(region_columns(r), i)
            print *, "column = ", column
            call compute_z_nodes(z_meshes(layer, column), node_val(h_mesh, column), min_bottom_layer_frac, &
                            top_depth(column), &
                            depth_is_constant, depth, depth_from_python, depth_function, &
                            depth_from_map, depth_vector(column),  have_min_depth, min_depth, &
                            sizing_is_constant, constant_sizing, list_sizing, sizing_function, sizing_vector, &
                            sigma_layers(layer), number_sigma_layers, radial_extrusion)
          end do
        end if

      end do layers

      call deallocate(region_columns(r))
    
    end do regions
    
    ! Now the tiresome business of making a shape function.
    h_dim = mesh_dim(h_mesh)
    call get_option("/geometry/quadrature/degree", quadrature_degree)
    quad = make_quadrature(vertices=h_dim + 2, dim=h_dim + 1, degree=quadrature_degree)
    full_shape = make_element_shape(vertices=h_dim + 2, dim=h_dim + 1, degree=1, quad=quad)
    call deallocate(quad)

    call get_option(trim(option_path)//'/name', mesh_name)

    ! combine the 1d vertical meshes into a full mesh
    call combine_z_meshes(h_mesh, z_meshes, out_mesh, &
       full_shape, mesh_name, option_path, sl=sigma_layers)
       
    do layer=1, nlayers
      do column=1, node_count(h_mesh)
        if (.not. node_owned(h_mesh, column)) cycle
        call deallocate(z_meshes(layer, column))
      end do
    end do
    call deallocate(full_shape)
    deallocate(z_meshes)
    deallocate(depth_vector)
    deallocate(sigma_layers)
    
  end subroutine extrude

  subroutine get_layer_extrusion_options(layer_option_path, &
                                   depth_is_constant, depth, depth_from_python, depth_function, depth_from_map, &
                                   file_name, have_min_depth, min_depth, surface_height, sizing_is_constant, constant_sizing, list_sizing, &
                                   sizing_function, sizing_vector, min_bottom_layer_frac, varies_only_in_depth, sigma_layers, number_sigma_layers, &
                                   radial_extrusion)

    character(len=*), intent(in) :: layer_option_path
    
    logical, intent(out) :: depth_is_constant, depth_from_python, depth_from_map
    real, intent(out) :: depth
    character(len=PYTHON_FUNC_LEN), intent(out) :: depth_function
    
    logical, intent(out) :: sizing_is_constant, list_sizing
    real, intent(out) :: constant_sizing
    character(len=PYTHON_FUNC_LEN), intent(out) :: sizing_function
    real, dimension(:), allocatable, intent(out) :: sizing_vector

    character(len=FIELD_NAME_LEN), intent(out) :: file_name
    logical, intent(out) :: have_min_depth
    real, intent(out) :: min_depth, surface_height
    
    logical, intent(out) :: varies_only_in_depth
    
    real, intent(out) :: min_bottom_layer_frac

    logical, intent(out) :: sigma_layers
    integer, intent(out) :: number_sigma_layers

    logical, intent(out) :: radial_extrusion
    
    integer, dimension(2) :: shape_option
    integer :: stat

    ! options under bottom_depth
    depth_from_python=.false.
    depth_from_map=.false.
    have_min_depth=.false.
    call get_option(trim(layer_option_path)//'/bottom_depth/constant', depth, stat=stat)
    if (stat==0) then
      depth_is_constant = .true.
    else
      depth_is_constant = .false.
      call get_option(trim(layer_option_path)//'/bottom_depth/python', depth_function, stat=stat)
      if (stat==0) then
        depth_from_python = .true.
      else
        call get_option(trim(layer_option_path)//'/bottom_depth/from_map/file_name', file_name, stat=stat)
        if (stat==0) then
          depth_from_map = .true.
        else
          FLAbort("Unknown way of specifying bottom depth function in mesh extrusion")
        end if
      end if
    end if

    if (depth_from_map) then
      call get_option(trim(layer_option_path)//'/bottom_depth/from_map/min_depth',min_depth, stat=stat)
      have_min_depth = stat==0
      call get_option(trim(layer_option_path)//'/bottom_depth/from_map/surface_height',surface_height, default=0.0)
    end if

    ! options under sizing_function
    list_sizing=.false.
    sigma_layers=.false.
    call get_option(trim(layer_option_path)//'/sizing_function/constant', constant_sizing, stat=stat)
    if (stat==0) then
      sizing_is_constant=.true.
    else
      sizing_is_constant=.false.
      call get_option(trim(layer_option_path)//'/sizing_function/python', sizing_function, stat=stat)
      if (stat/=0) then
        if (have_option(trim(layer_option_path)//"/sizing_function/list")) then
          list_sizing=.true.
          shape_option=option_shape(trim(layer_option_path)//"/sizing_function/list")
          allocate(sizing_vector(1:shape_option(1)))
          call get_option(trim(layer_option_path)//'/sizing_function/list', sizing_vector)
        else
          call get_option(trim(layer_option_path)//'/sizing_function/sigma_layers/standard', &
              number_sigma_layers, stat=stat)
          if (stat==0) then
            sigma_layers = .true.
          else
            FLAbort("Unknown way of specifying sizing function in mesh extrusion")
          end if
        end if
      end if
    end if

    varies_only_in_depth = have_option(trim(layer_option_path)//'/sizing_function/varies_only_in_depth')
  
    call get_option(trim(layer_option_path)//'/minimum_bottom_layer_fraction', &
                    min_bottom_layer_frac, default=1.e-3)
  
  end subroutine get_layer_extrusion_options

  subroutine populate_depth_vector(h_mesh,file_name,depth_vector,surface_height,radial_extrusion)

    type(vector_field), intent(in) :: h_mesh
    character(len=FIELD_NAME_LEN), intent(in):: file_name
    real, intent(in) :: surface_height
    real, dimension(:,:), allocatable :: tmp_pos_vector
    real, dimension(:), intent(inout) :: depth_vector
    logical :: radial_extrusion

    integer :: column

    if(radial_extrusion) then

      allocate(tmp_pos_vector(mesh_dim(h_mesh)+1, size(depth_vector)))

      do column=1, node_count(h_mesh)
        tmp_pos_vector(:,column) = node_val(h_mesh, column)
      end do

      call set_from_map(trim(file_name)//char(0), tmp_pos_vector(1,:), tmp_pos_vector(2,:), tmp_pos_vector(3,:), &
                                                                  depth_vector, size(depth_vector), surface_height)

    else

      allocate(tmp_pos_vector(mesh_dim(h_mesh), size(depth_vector)))

      do column=1, node_count(h_mesh)
        tmp_pos_vector(:,column) = node_val(h_mesh, column)
      end do

      call set_from_map_beta(trim(file_name)//char(0), tmp_pos_vector(1,:), tmp_pos_vector(2,:), &
                                                  depth_vector, size(depth_vector), surface_height)

    end if

    if (associated(h_mesh%mesh%halos)) then
      call halo_update(h_mesh%mesh%halos(2), depth_vector)
    end if

    deallocate(tmp_pos_vector)

  end subroutine populate_depth_vector

  subroutine compute_z_nodes_wrapper(z_mesh, xy, min_bottom_layer_frac, &
                                     top_depth, &
                                     depth_is_constant, depth, depth_from_python, depth_function, &
                                     depth_from_map, map_depth, have_min_depth, min_depth, &
                                     sizing_is_constant, constant_sizing, list_sizing, sizing_function, sizing_vector, &
                                     sigma_layers, number_sigma_layers, radial_extrusion)

    type(vector_field), intent(out) :: z_mesh
    real, dimension(:), intent(in) :: xy
    real, intent(in) :: min_bottom_layer_frac
    real, intent(inout) :: top_depth ! IN: depth of top node, OUT: depth of bottom node
    logical, intent(in) :: depth_is_constant, sizing_is_constant, depth_from_python, depth_from_map, list_sizing
    logical, intent(in) :: have_min_depth, sigma_layers
    real, intent(in) :: map_depth, min_depth
    real, intent(in) :: depth, constant_sizing
    character(len=*), intent(in) :: depth_function, sizing_function
    real, dimension(:), allocatable, intent(in) :: sizing_vector
    integer, intent(in) :: number_sigma_layers
    logical, intent(in) :: radial_extrusion

    real, dimension(1) :: tmp_depth
    real, dimension(size(xy), 1) :: tmp_pos
    real :: ldepth
    
    if(depth_is_constant) then
      ldepth = depth
    else 
      tmp_pos(:,1) = xy
      if (depth_from_python) then
        call set_from_python_function(tmp_depth, trim(depth_function), tmp_pos, time=0.0)
        ldepth = tmp_depth(1)
      else if (depth_from_map) then
         ldepth = map_depth
         if (have_min_depth) then
           if (ldepth < min_depth) ldepth=min_depth
         end if
      else
        FLAbort("Unknown way of specifying the bottom_depth.")
      end if
    end if
    
    if (sizing_is_constant) then
      call compute_z_nodes(z_mesh, ldepth, xy, top_depth, &
       min_bottom_layer_frac, radial_extrusion, sizing=constant_sizing)
    else
      if (list_sizing) then
        call compute_z_nodes(z_mesh, ldepth, xy, top_depth, &
        min_bottom_layer_frac, radial_extrusion, sizing_vector=sizing_vector)
      else if (sigma_layers) then
        call compute_z_nodes(z_mesh, ldepth, xy, top_depth, &
        min_bottom_layer_frac, radial_extrusion, number_sigma_layers=number_sigma_layers)
      else
        call compute_z_nodes(z_mesh, ldepth, xy, top_depth, &
        min_bottom_layer_frac, radial_extrusion, sizing_function=sizing_function)
      end if
    end if
    
  end subroutine compute_z_nodes_wrapper

  subroutine get_previous_z_nodes(z_mesh, z_mesh_previous)
    type(vector_field), intent(inout) :: z_mesh, z_mesh_previous
    z_mesh = z_mesh_previous
    call incref(z_mesh)
  end subroutine get_previous_z_nodes

  subroutine compute_z_nodes_sizing(z_mesh, depth, xy, top_depth, min_bottom_layer_frac, &
                                    radial_extrusion, sizing, &
                                    sizing_function, sizing_vector, number_sigma_layers)
    !!< Figure out at what depths to put the layers.
    type(vector_field), intent(out) :: z_mesh
    real, intent(in):: depth
    real, dimension(:), intent(in):: xy
    real, intent(inout) :: top_depth ! IN: depth of top node, OUT: depth of bottom node
    ! to prevent infinitesimally thin bottom layer if sizing function
    ! is an integer mulitple of total depth, the bottom layer needs
    ! to have at least this fraction of the layer depth above it.
    ! The recommended value is 1e-3.
    real, intent(in) :: min_bottom_layer_frac
    logical, intent(in) :: radial_extrusion
    real, optional, intent(in):: sizing
    character(len=*), optional, intent(in):: sizing_function
    real, dimension(:), optional, intent(in) :: sizing_vector
    integer, optional, intent(in) :: number_sigma_layers
    ! this is a safety gap:
    integer, parameter:: MAX_VERTICAL_NODES=1e6

    integer :: elements
    logical :: is_constant
    real :: constant_value

    type(rlist):: depths
    type(mesh_type) :: mesh
    type(element_type) :: oned_shape
    type(quadrature_type) :: oned_quad
    integer :: quadrature_degree
    integer :: ele
    integer, parameter :: loc=2
    integer :: node
    real, dimension(:), allocatable:: xyz
    real, dimension(size(xy)) :: radial_dir
    real :: delta_h, d
    character(len=PYTHON_FUNC_LEN) :: py_func
    integer :: list_size

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
    else if (present(sizing_vector)) then
      is_constant=.false.
      constant_value=-1.0
      list_size=size(sizing_vector)
    else if (present(number_sigma_layers)) then
      is_constant=.true.
      constant_value=(depth-top_depth)/float(number_sigma_layers)
      py_func = " "
    else
      FLAbort("Need to supply either sizing or sizing_function")
    end if

    ! Start the mesh at the top (d=0 for the top layer) and work down to d=-depth.
    d = top_depth
    node=2
    ! first size(xy) coordinates remain fixed, 
    ! the last entry will be replaced with the appropriate depth
    if (radial_extrusion) then
      allocate(xyz(size(xy)))
    else
      allocate(xyz(size(xy)+1))
    end if
    xyz(1:size(xy))=xy
    radial_dir = 0.0
    if (radial_extrusion) radial_dir = xy/sqrt(sum(xy**2))
    call insert(depths, d)
    do
      if (radial_extrusion) then
        xyz = xy - radial_dir*d
      else
        xyz(size(xy)+1)=-d
      end if
      if (present(sizing_vector)) then
        if ((node-1)<=list_size) then
          delta_h = sizing_vector(node-1)
        else
          delta_h = sizing_vector(list_size)
        end if
        node=node+1
      else
        delta_h = get_delta_h( xyz, is_constant, constant_value, py_func)
      end if
      d=d + sign(delta_h, depth)
      if (abs(d)>abs(depth)-min_bottom_layer_frac*delta_h) exit
      call insert(depths, d)
      if (depths%length>MAX_VERTICAL_NODES) then
        ewrite(-1,*) "Check your extrude/sizing_function"
        FLExit("Maximum number of vertical layers reached")
      end if
    end do
    call insert(depths, depth)
    elements=depths%length-1
    top_depth = depth

    call allocate(mesh, elements+1, elements, oned_shape, "ZMesh")
    do ele=1,elements
      mesh%ndglno((ele-1) * loc + 1: ele*loc) = (/ele, ele+1/)
    end do

    call allocate(z_mesh, 1, mesh, "ZMeshCoordinates")
    call deallocate(mesh)
    call deallocate(oned_shape)

    do node=1, elements+1
      call set(z_mesh, node, (/ -pop(depths) /) )
    end do
    deallocate(xyz)

    ! For pathological sizing functions the mesh might have gotten inverted at the last step.
    ! If you encounter this, make this logic smarter.
    assert(abs(node_val(z_mesh, 1, elements)) < abs(node_val(z_mesh, 1, elements+1)))
    
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
      
  end subroutine compute_z_nodes_sizing

  ! hadapt_extrude options checking
  subroutine hadapt_extrude_check_options

    integer :: nmeshes, m, nregions, r
    character(len=OPTION_PATH_LEN) :: mesh_path

    nmeshes=option_count("/geometry/mesh")
    do m = 0, nmeshes-1
      mesh_path="/geometry/mesh["//int2str(m)//"]"
      nregions=option_count(trim(mesh_path)//'/from_mesh/extrude/regions')
      if (nregions>1) then
        ! we're using region ids to extrude
        if (have_option('/mesh_adaptivity/hr_adaptivity') &
           .and. .not. have_option('/mesh_adaptivity/hr_adaptivity/preserve_mesh_regions')) then
          ewrite(-1,*) "You are using region ids to specify mesh extrusion"
          ewrite(-1,*) "However in your adaptivity settings you have not selected " // &
            & "/mesh_adaptivity/hr_adaptivity/preserve_mesh_regions"
          ewrite(-1,*) "This means fluidity will not be able to extrude your mesh again after the adapt."
          FLExit("Missing /mesh_adaptivity/hr_adaptivity/preserve_mesh_regions option")
        end if
      end if   
        
    end do

  end subroutine hadapt_extrude_check_options

    
end module hadapt_extrude
