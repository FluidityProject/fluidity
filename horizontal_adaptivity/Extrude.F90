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
  use vtk_interfaces
  use linked_lists
  use hadapt_combine_meshes
  use halos
  implicit none

  private
  
  public :: extrude, compute_z_nodes, hadapt_extrude_check_options, get_extrusion_options, &             
            populate_depth_vector, skip_column_extrude

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

    character(len=FIELD_NAME_LEN):: mesh_name, file_name  
    type(quadrature_type) :: quad
    type(element_type) :: full_shape
    type(vector_field) :: constant_z_mesh
    type(vector_field), dimension(node_count(h_mesh)) :: z_meshes
    character(len=PYTHON_FUNC_LEN) :: sizing_function, depth_function
    logical:: depth_from_python, depth_from_map, have_min_depth
    real, dimension(:), allocatable :: depth_vector
    real:: min_depth, surface_height
    logical:: sizing_is_constant, depth_is_constant, varies_only_in_z
    real:: constant_sizing, depth, min_bottom_layer_frac
    integer:: h_dim, column, quadrature_degree
    
    integer :: n_regions, r
    integer, dimension(:), allocatable :: region_ids
    logical :: apply_region_ids, constant_z_mesh_initialised
    integer, dimension(node_count(h_mesh)) :: visited
    logical, dimension(node_count(h_mesh)) :: column_visited

    !! Checking linearity of h_mesh.
    assert(h_mesh%mesh%shape%degree == 1)
    assert(h_mesh%mesh%continuity >= 0)

    call add_nelist(h_mesh%mesh)
    
    n_regions = option_count(trim(option_path)//'/from_mesh/extrude/regions')
    if(n_regions==0) then
      ewrite(-1,*) "Since r13369 it has been possible to extrude using different parameters"
      ewrite(-1,*) "in each region id of a horizontal mesh.  This means that the extrusion parameters"
      ewrite(-1,*) "must be included under geometry/mesh::MeshNameHere/from_mesh/extrude/regions."
      ewrite(-1,*) "This is necessary even if you want to use the same parameters"
      ewrite(-1,*) "for the whole mesh (using a single regions option and no region_ids)."
      ewrite(-1,*) "I've been told to extrude but have found no regions options."
      ewrite(-1,*) "Have you updated your flml since r13369?"
      FLExit("No regions options found under extrude.")
    elseif(n_regions<0) then
      FLAbort("Negative number of regions options found under extrude.")
    end if
    apply_region_ids = (n_regions>1)
    visited = 0 ! a little debugging check - can be removed later
    
    column_visited = .false.
    
    do r = 0, n_regions-1
      
      constant_z_mesh_initialised = .false.
      
      call get_extrusion_options(option_path, r, apply_region_ids, region_ids, &
                                 depth_is_constant, depth, depth_from_python, depth_function, depth_from_map, &
                                 file_name, have_min_depth, min_depth, surface_height, sizing_is_constant, constant_sizing, &
                                 sizing_function, min_bottom_layer_frac, varies_only_in_z)

      allocate(depth_vector(size(z_meshes)))
      if (depth_from_map) call populate_depth_vector(h_mesh,file_name,depth_vector,surface_height)
      
      ! create a 1d vertical mesh under each surface node
      do column=1, size(z_meshes)
      
        ! decide if this column needs visiting...
        if(skip_column_extrude(h_mesh%mesh, column, &
                              apply_region_ids, column_visited(column), region_ids, &
                              visited_count = visited(column))) cycle
        
        if(varies_only_in_z .and. depth_is_constant) then
          if (.not. constant_z_mesh_initialised) then
            call compute_z_nodes(constant_z_mesh, node_val(h_mesh, column), min_bottom_layer_frac, &
                            depth_is_constant, depth, depth_from_python, depth_function, &
                            depth_from_map, depth_vector(column),  have_min_depth, min_depth, &
                            sizing_is_constant, constant_sizing, sizing_function)
            constant_z_mesh_initialised = .true.
          end if
          call get_previous_z_nodes(z_meshes(column), constant_z_mesh)
        else
          call compute_z_nodes(z_meshes(column), node_val(h_mesh, column), min_bottom_layer_frac, &
                            depth_is_constant, depth, depth_from_python, depth_function, &
                            depth_from_map, depth_vector(column),  have_min_depth, min_depth, &
                            sizing_is_constant, constant_sizing, sizing_function)
        end if

      end do
      
      if(apply_region_ids) deallocate(region_ids)
      deallocate(depth_vector)
      
      if (constant_z_mesh_initialised) then
        call deallocate(constant_z_mesh)
      end if
    
    end do
    
#ifdef DDEBUG
    if(apply_region_ids) then
      ewrite(2,*) "Maximum number of times a node was visited: ", maxval(visited)
      ewrite(2,*) "Minimum number of times a node was visited: ", minval(visited)
      if(.not.isparallel()) then
        assert(minval(visited)>0)
      end if
    end if
#endif
      
    ! Now the tiresome business of making a shape function.
    h_dim = mesh_dim(h_mesh)
    call get_option("/geometry/quadrature/degree", quadrature_degree)
    quad = make_quadrature(vertices=h_dim + 2, dim=h_dim + 1, degree=quadrature_degree)
    full_shape = make_element_shape(vertices=h_dim + 2, dim=h_dim + 1, degree=1, quad=quad)
    call deallocate(quad)

    call get_option(trim(option_path)//'/name', mesh_name)

    ! combine the 1d vertical meshes into a full mesh
    call combine_z_meshes(h_mesh, z_meshes, out_mesh, &
       full_shape, mesh_name, option_path)
       
    do column=1, node_count(h_mesh)
      if (.not. node_owned(h_mesh, column)) cycle
      call deallocate(z_meshes(column))
    end do
    call deallocate(full_shape)
    
  end subroutine extrude

  subroutine get_extrusion_options(option_path, region_index, apply_region_ids, region_ids, &
                                   depth_is_constant, depth, depth_from_python, depth_function, depth_from_map, &
                                   file_name, have_min_depth, min_depth, surface_height, sizing_is_constant, constant_sizing, &
                                   sizing_function, min_bottom_layer_frac, varies_only_in_z)

    character(len=*), intent(in) :: option_path
    integer, intent(in) :: region_index
    logical, intent(in) :: apply_region_ids
    
    integer, dimension(:), allocatable :: region_ids
    
    logical, intent(out) :: depth_is_constant, depth_from_python, depth_from_map
    real, intent(out) :: depth
    character(len=PYTHON_FUNC_LEN), intent(out) :: depth_function
    
    logical, intent(out) :: sizing_is_constant
    real, intent(out) :: constant_sizing
    character(len=PYTHON_FUNC_LEN), intent(out) :: sizing_function

    character(len=FIELD_NAME_LEN), intent(out) :: file_name
    logical, intent(out) :: have_min_depth
    real, intent(out) :: min_depth, surface_height
    
    logical, intent(out) :: varies_only_in_z
    
    real, intent(out) :: min_bottom_layer_frac
    
    integer, dimension(2) :: shape_option
    integer :: stat
  
    if(apply_region_ids) then
      shape_option=option_shape(trim(option_path)//"/from_mesh/extrude/regions["//int2str(region_index)//"]/region_ids")
      allocate(region_ids(1:shape_option(1)))
      call get_option(trim(option_path)//"/from_mesh/extrude/regions["//int2str(region_index)//"]/region_ids", region_ids)
    end if

    ! get the extrusion options
    depth_from_python=.false.
    depth_from_map=.false.
    have_min_depth=.false.
    call get_option(trim(option_path)//&
                    '/from_mesh/extrude/regions['//int2str(region_index)//&
                    ']/bottom_depth/constant', &
                      depth, stat=stat)
    if (stat==0) then
      depth_is_constant = .true.
    else
      depth_is_constant = .false.
      call get_option(trim(option_path)//&
                      '/from_mesh/extrude/regions['//int2str(region_index)//&
                      ']/bottom_depth/python', &
                       depth_function, stat=stat)
      if (stat==0) depth_from_python = .true.
      if (stat /= 0) then 
        call get_option(trim(option_path)//'/from_mesh/extrude/regions['//int2str(region_index)//&
                         ']/bottom_depth/from_map/file_name', &
                          file_name, stat=stat)
        if (stat==0) depth_from_map = .true.
      end if
      if (stat /= 0) then
        FLAbort("Unknown way of specifying bottom depth function in mesh extrusion")
      end if
    end if

    if (have_option(trim(option_path)//'/from_mesh/extrude/regions['//int2str(region_index)//&
                                         ']/bottom_depth/from_map/min_depth')) then
      have_min_depth=.true.
      call get_option(trim(option_path)//'/from_mesh/extrude/regions['//int2str(region_index)//&
                                         ']/bottom_depth/from_map/min_depth',min_depth)
    end if

    surface_height=0.0
    if (have_option(trim(option_path)//'/from_mesh/extrude/regions['//int2str(region_index)//&
                                         ']/bottom_depth/from_map/surface_height')) then
      call get_option(trim(option_path)//'/from_mesh/extrude/regions['//int2str(region_index)//&
                                         ']/bottom_depth/from_map/surface_height',surface_height)
    end if
    
    call get_option(trim(option_path)//&
                    '/from_mesh/extrude/regions['//int2str(region_index)//&
                    ']/sizing_function/constant', &
                    constant_sizing, stat=stat)
    if (stat==0) then
      sizing_is_constant=.true.
    else
      sizing_is_constant=.false.
      call get_option(trim(option_path)//&
                      '/from_mesh/extrude/regions['//int2str(region_index)//&
                      ']/sizing_function/python', &
                      sizing_function, stat=stat)
      if (stat/=0) then
        FLAbort("Unknown way of specifying sizing function in mesh extrusion")
      end if       
    end if

    varies_only_in_z = have_option(trim(option_path)//&
    '/from_mesh/extrude/regions['//int2str(region_index)//&
    ']/sizing_function/varies_only_in_z')
  
    call get_option(trim(option_path)//&
                    '/from_mesh/extrude/regions['//int2str(region_index)//&
                    ']/minimum_bottom_layer_fraction', &
                    min_bottom_layer_frac, default=1.e-3)
  
  end subroutine get_extrusion_options

  subroutine populate_depth_vector(h_mesh,file_name,depth_vector,surface_height)

    type(vector_field), intent(in) :: h_mesh
    character(len=FIELD_NAME_LEN), intent(in):: file_name
    real, intent(in) :: surface_height
    real, dimension(:,:), allocatable :: tmp_pos_vector
    real, dimension(:), intent(inout) :: depth_vector

    integer :: column

    allocate(tmp_pos_vector(mesh_dim(h_mesh),node_count(h_mesh)))

    do column=1, node_count(h_mesh)
      tmp_pos_vector(:,column) = node_val(h_mesh, column)
    end do

    call set_from_map_beta(trim(file_name)//char(0), tmp_pos_vector(1,:), tmp_pos_vector(2,:), depth_vector, node_count(h_mesh), surface_height)

    if (associated(h_mesh%mesh%halos)) then
      call halo_update(h_mesh%mesh%halos(2), depth_vector)
    end if

    deallocate(tmp_pos_vector)

  end subroutine populate_depth_vector

  subroutine compute_z_nodes_wrapper(z_mesh, xy, min_bottom_layer_frac, &
                                     depth_is_constant, depth, depth_from_python, depth_function, &
                                     depth_from_map, map_depth, have_min_depth, min_depth, &
                                     sizing_is_constant, constant_sizing, sizing_function)

    type(vector_field), intent(out) :: z_mesh
    real, dimension(:), intent(in) :: xy
    real, intent(in) :: min_bottom_layer_frac
    logical, intent(in) :: depth_is_constant, sizing_is_constant, depth_from_python, depth_from_map
    logical, intent(in) :: have_min_depth
    real, intent(in) :: map_depth, min_depth
    real, intent(in) :: depth, constant_sizing
    character(len=*), intent(in) :: depth_function, sizing_function

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
      call compute_z_nodes(z_mesh, ldepth, xy, &
       min_bottom_layer_frac, sizing=constant_sizing)
    else
        call compute_z_nodes(z_mesh, ldepth, xy, &
        min_bottom_layer_frac, sizing_function=sizing_function)
    end if
    
  end subroutine compute_z_nodes_wrapper

  subroutine get_previous_z_nodes(z_mesh, z_mesh_previous)
    type(vector_field), intent(inout) :: z_mesh, z_mesh_previous
    z_mesh = z_mesh_previous
    call incref(z_mesh)
  end subroutine get_previous_z_nodes

  subroutine compute_z_nodes_sizing(z_mesh, depth, xy, min_bottom_layer_frac, sizing, sizing_function)
    !!< Figure out at what depths to put the layers.
    type(vector_field), intent(out) :: z_mesh
    real, intent(in):: depth
    real, dimension(:), intent(in):: xy
    ! to prevent infinitesimally thin bottom layer if sizing function
    ! is an integer mulitple of total depth, the bottom layer needs
    ! to have at least this fraction of the layer depth above it.
    ! The recommended value is 1e-3.
    real, intent(in) :: min_bottom_layer_frac
    real, optional, intent(in):: sizing
    character(len=*), optional, intent(in):: sizing_function
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

    ! Start the mesh at z=0 and work down to z=-depth.
    z=0.0
    ! first size(xy) coordinates remain fixed, 
    ! the last entry will be replaced with the appropriate depth
    xyz(1:size(xy))=xy
    call insert(depths, z)
    do
      xyz(size(xy)+1)=z
      delta_h = get_delta_h( xyz, is_constant, constant_value, py_func)
      z=z - delta_h
      if (z<-depth+min_bottom_layer_frac*delta_h) exit
      call insert(depths, z)
      if (depths%length>MAX_VERTICAL_NODES) then
        ewrite(-1,*) "Check your extrude/sizing_function"
        FLExit("Maximum number of vertical layers reached")
      end if
    end do
    call insert(depths, -depth)
    elements=depths%length-1

    call allocate(mesh, elements+1, elements, oned_shape, "ZMesh")
    call deallocate(oned_shape)
    do ele=1,elements
      mesh%ndglno((ele-1) * loc + 1: ele*loc) = (/ele, ele+1/)
    end do

    call allocate(z_mesh, 1, mesh, "ZMeshCoordinates")
    call deallocate(mesh)

    call set(z_mesh, 1, (/0.0/))
    do node=1, elements+1
      call set(z_mesh, node,  (/ pop(depths) /))
    end do

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
      
  end subroutine compute_z_nodes_sizing

  logical function skip_column_extrude(horizontal_mesh, column, &
                                       apply_region_ids, column_visited, region_ids, &
                                       visited_count)
    !!< this function decides if a column need extruding or not
    type(mesh_type), intent(in) :: horizontal_mesh
    integer, intent(in) :: column
    logical, intent(in) :: apply_region_ids
    logical, intent(inout) :: column_visited
    integer, dimension(:), intent(in) :: region_ids
    integer, intent(inout), optional :: visited_count
    
    integer, dimension(:), pointer :: eles
    logical :: node_in_region
    integer :: rs
    
    skip_column_extrude = .false.
    if(.not.node_owned(horizontal_mesh, column)) then
      skip_column_extrude = .true.
      return
    end if
    
    ! need to work out here if this column is in one of the current region ids!
    ! this is a bit arbitrary since nodes belong to multiple regions... therefore
    ! the extrusion depth had better be continuous across region id boundaries!
    if(apply_region_ids) then
      if(column_visited) then
        skip_column_extrude = .true.
        return
      end if
      eles => node_neigh(horizontal_mesh, column)
      node_in_region = .false.
      region_id_loop: do rs = 1, size(region_ids)
        if(any(region_ids(rs)==horizontal_mesh%region_ids(eles))) then
          node_in_region = .true.
          exit region_id_loop
        end if
      end do region_id_loop
      if(.not.node_in_region) then
        skip_column_extrude = .true.
        return
      end if
      column_visited=.true.
      if(present(visited_count)) then
        visited_count = visited_count + 1
      end if
    end if
  
  end function skip_column_extrude

  ! hadapt_extrude options checking
  subroutine hadapt_extrude_check_options

    integer :: nmeshes, m, nregions, r
    character(len=OPTION_PATH_LEN) :: mesh_path

    ! Extruding along bathymetry currently only works with radial extrusions on the sphere.
    ! If bottom_depth/from_map is selected check also that /geometry/spherical_earth is enabled.
    nmeshes=option_count("/geometry/mesh")
    do m = 0, nmeshes-1
      mesh_path="/geometry/mesh["//int2str(m)//"]"
      nregions=option_count(trim(mesh_path)//'/from_mesh/extrude/regions')
      do r = 0, nregions-1
        if ((have_option(trim(mesh_path)//'/from_mesh/extrude/regions['//int2str(r)//']/bottom_depth/from_map/surface_height')).and. &
          (have_option('/geometry/spherical_earth'))) then
          ewrite(-1,*) "Warning. You have selected a spherical earth geometry and have also specified a surface_height"
          ewrite(-1,*) "in your extrusion. The surface height option will have no effect."
        end if
      end do
    end do

  end subroutine hadapt_extrude_check_options

    
end module hadapt_extrude
