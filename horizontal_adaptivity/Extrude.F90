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
  use vtk_interfaces
  use linked_lists
  implicit none

  private
  
  public :: extrude, compute_z_nodes, hadapt_extrude_check_options

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
    
    integer :: n_regions, r, rs
    integer, dimension(2) :: shape_option
    integer, dimension(:), allocatable :: region_ids
    logical :: apply_region_ids
    integer, dimension(:), pointer :: eles
    integer, dimension(node_count(h_mesh)) :: visited
    logical :: node_in_region
    
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
    
      if(apply_region_ids) then
        shape_option=option_shape(trim(option_path)//"/from_mesh/extrude/regions["//int2str(r)//"]/region_ids")
        allocate(region_ids(1:shape_option(1)))
        call get_option(trim(option_path)//"/from_mesh/extrude/regions["//int2str(r)//"]/region_ids", region_ids)
      end if

      ! get the extrusion options
      call get_option(trim(option_path)//&
                      '/from_mesh/extrude/regions['//int2str(r)//&
                      ']/bottom_depth/constant', &
                      depth, stat=stat)
      if (stat==0) then
        depth_is_constant = .true.
      else
        depth_is_constant = .false.
        call get_option(trim(option_path)//&
                        '/from_mesh/extrude/regions['//int2str(r)//&
                        ']/bottom_depth/python', &
                         depth_function, stat=stat)
        if (stat /= 0) then
          FLAbort("Unknown way of specifying bottom depth function in mesh extrusion")
        end if
      end if
      
      call get_option(trim(option_path)//&
                      '/from_mesh/extrude/regions['//int2str(r)//&
                      ']/sizing_function/constant', &
                      constant_sizing, stat=stat)
      if (stat==0) then
        sizing_is_constant=.true.
      else
        sizing_is_constant=.false.
        call get_option(trim(option_path)//&
                        '/from_mesh/extrude/regions['//int2str(r)//&
                        ']/sizing_function/python', &
                        sizing_function, stat=stat)
        if (stat/=0) then
          FLAbort("Unknown way of specifying sizing function in mesh extrusion")
        end if       
      end if
        
      ! create a 1d vertical mesh under each surface node
      do column=1, size(z_meshes)
      
        if (.not. node_owned(h_mesh, column)) cycle

        ! need to work out here if this column is in one of the current region ids!
        ! this is a bit arbitrary since nodes belong to multiple regions... therefore
        ! the extrusion depth had better be continuous across region id boundaries!
        if(apply_region_ids) then
          if(column_visited(column)) cycle
          eles => node_neigh(h_mesh, column)
          node_in_region = .false.
          region_id_loop: do rs=1, size(region_ids)
            if(any(region_ids(rs)==h_mesh%mesh%region_ids(eles))) then
              node_in_region = .true.
              exit region_id_loop
            end if
          end do region_id_loop
          if(.not. node_in_region) cycle
          column_visited(column) = .true.
          visited(column) = visited(column) + 1
        end if
        
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
      
      if(apply_region_ids) deallocate(region_ids)
    
    end do
    
#ifdef DDEBUG
    if(apply_region_ids) then
      assert(minval(visited)>0)
      ewrite(2,*) "Maximum number of times a node was visited: ", maxval(visited)
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
    
    !call vtk_write_fields("extruded_mesh", 0, out_mesh, out_mesh%mesh)
    !call vtk_write_surface_mesh("surface_mesh", 1, out_mesh)

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
      if (z<-depth+MIN_BOTTOM_LAYER_FRAC*delta_h) exit
      call insert(depths, z)
      if (depths%length>MAX_VERTICAL_NODES) then
        ewrite(-1,*) "Check your extrude/sizing_function"
        FLAbort("Maximum number of vertical layers reached")
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
      
  end subroutine compute_z_nodes

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
        if ((have_option(trim(mesh_path)//'/from_mesh/extrude/regions['//int2str(r)//']/bottom_depth/from_map')).and. &
          (.not.have_option('/geometry/spherical_earth'))) then
          ewrite(-1,*) "Extruding along bathymetry currently only works with radial extrusions on the sphere."
          ewrite(-1,*) "Please turn on the geometry/spherical_earth option."
          FLExit("Using extrude/from_map requires the gemoetry/spherical_earth option to be enabled.")
        end if
      end do
    end do

  end subroutine hadapt_extrude_check_options

    
end module hadapt_extrude
