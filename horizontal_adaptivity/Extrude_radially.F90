#include "fdebug.h"

module hadapt_extrude_radially
  !!< Extrude a given pseudo 2D mesh on a spherical shell to a full 3D mesh.
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
  use vector_tools
  use vtk_interfaces
  use linked_lists
  implicit none

  private
  
  public :: extrude_radially, compute_r_nodes

  contains

  subroutine extrude_radially(shell_mesh, option_path, out_mesh)
    !!< The pseudo 2D mesh on the spherical shell.
    !!< Note: this must be linear.
    type(vector_field), intent(inout) :: shell_mesh
    !!< options to be set for out_mesh,
    !!< at the moment: /name, and under from_mesh/extrude/:
    !!< depth, sizing_function optionally top_surface_id and bottom_surface_id
    character(len=*), intent(in) :: option_path
    !!< The full extruded 3D mesh.
    type(vector_field), intent(out) :: out_mesh

    character(len=FIELD_NAME_LEN):: mesh_name, file_name
    type(quadrature_type) :: quad
    type(element_type) :: full_shape
    type(vector_field), dimension(node_count(shell_mesh)) :: r_meshes
    character(len=PYTHON_FUNC_LEN) :: sizing_function, depth_function
    logical:: sizing_is_constant, depth_is_constant
    real:: constant_sizing, depth
    real, dimension(:), allocatable :: sizing_vector
    integer :: stat, shell_dim, column, quadrature_degree
    real :: r_shell
    
    real, dimension(1) :: tmp_depth
    real, dimension(mesh_dim(shell_mesh)+1, 1) :: tmp_pos
    real, dimension(:,:), allocatable :: tmp_pos_vector
    real, dimension(:), allocatable :: depth_vector

    logical :: have_min_depth=.false.
    real :: min_depth

    integer :: n_regions, r, rs
    integer, dimension(2) :: shape_option
    integer, dimension(:), allocatable :: region_ids
    logical :: apply_region_ids
    integer, dimension(:), pointer :: eles
    integer, dimension(node_count(shell_mesh)) :: visited
    logical :: node_in_region

    logical, dimension(node_count(shell_mesh)) :: column_visited

    write(*,*) 'In extrude_radially'

    !! We assume that for the extrusion operation,
    !! the layer configuration is independent of phi and theta.
    !! (i.e. the layers are the same everywhere.)

    !! Checking linearity of shell_mesh.
    assert(shell_mesh%mesh%shape%degree == 1)
    assert(shell_mesh%mesh%continuity >= 0)

    call add_nelist(shell_mesh%mesh)

    n_regions = option_count(trim(option_path)//'/from_mesh/extrude/regions')
    if(n_regions==0) then
      ewrite(-1,*) "Since r13369 it has been possible to extrude using different parameters"
      ewrite(-1,*) "in each region id of a shell mesh.  This means that the extrusion parameters"
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
        if (stat /= 0) call get_option(trim(option_path)//'/from_mesh/extrude/regions['//int2str(r)//&
                                       ']/bottom_depth/from_map/file_name', &
            file_name, stat=stat)
        if (stat /= 0) then
          FLAbort("Unknown way of specifying bottom depth function in mesh extrusion")
        end if
      end if
      
      if (have_option(trim(option_path)//'/from_mesh/extrude/regions['//int2str(r)//&
                                         ']/bottom_depth/from_map/min_depth')) then
        have_min_depth=.true.
        call get_option(trim(option_path)//'/from_mesh/extrude/regions['//int2str(r)//&
                                           ']/bottom_depth/from_map/min_depth',min_depth)
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
        if (have_option(trim(option_path)//"/from_mesh/extrude/regions["//&
                                      int2str(r)//"]/sizing_function/list")) then
            shape_option=option_shape(trim(option_path)//"/from_mesh/extrude/regions["//&
                                      int2str(r)//"]/sizing_function/list")
            allocate(sizing_vector(1:shape_option(1)))
            call get_option(trim(option_path)//'/from_mesh/extrude/regions['//&
                                      int2str(r)//']/sizing_function/list', &
                                      sizing_vector, stat=stat)
        end if
        if (stat/=0) then
          FLAbort("Unknown way of specifying sizing function in mesh extrusion")
        end if       
      end if

      if (have_option(trim(option_path)//'/from_mesh/extrude/regions['//int2str(r)//&
                                         ']/bottom_depth/from_map/')) then

        allocate (tmp_pos_vector(mesh_dim(shell_mesh)+1,size(r_meshes)), depth_vector(size(r_meshes)))

        do column=1, size(r_meshes)
          tmp_pos_vector(:,column) = node_val(shell_mesh, column)
        end do

        call set_from_map(trim(file_name), tmp_pos_vector(1,:), tmp_pos_vector(2,:), tmp_pos_vector(3,:), depth_vector, size(r_meshes))

        deallocate (tmp_pos_vector)

      end if
      
      ! create a 1d radial mesh under each surface node
      do column=1, size(r_meshes)
      
        if (.not. node_owned(shell_mesh, column)) cycle
        
        ! need to work out here if this column is in one of the current region ids!
        ! this is a bit arbitrary since nodes belong to multiple regions... therefore
        ! the extrusion depth had better be continuous across region id boundaries!
        if(apply_region_ids) then
          if(column_visited(column)) cycle
          eles => node_neigh(shell_mesh, column)
          node_in_region = .false.
          region_id_loop: do rs=1, size(region_ids)
            if(any(region_ids(rs)==shell_mesh%mesh%region_ids(eles))) then
              node_in_region = .true.
              exit region_id_loop
            end if
          end do region_id_loop
          if(.not. node_in_region) cycle
          column_visited(column) = .true.
          visited(column) = visited(column) + 1
        end if
        
        if(.not.depth_is_constant) then
          tmp_pos(:,1) = node_val(shell_mesh, column)
          if (have_option(trim(option_path)//'/from_mesh/extrude/regions['//&
                                      int2str(r)//']/bottom_depth/python')) then
            call set_from_python_function(tmp_depth, trim(depth_function), tmp_pos, time=0.0)
            depth = tmp_depth(1)
          else if  (have_option(trim(option_path)//'/from_mesh/extrude/regions['//&
                                      int2str(r)//']/bottom_depth/from_map')) then
            depth = depth_vector(column) 
            if (have_min_depth) then
              if (depth < min_depth) depth=min_depth
            end if
          else
            FLAbort("Error with options path")
          end if
        end if
        
        if (sizing_is_constant) then
          call compute_r_nodes(r_meshes(column), depth, node_val(shell_mesh, column), r_shell, &
            sizing=constant_sizing) ! Return r_shell, will be needed for working out face id's later
        else
          call compute_r_nodes(r_meshes(column), depth, node_val(shell_mesh, column), r_shell, &
            sizing_function=sizing_function)
        end if
        
      end do
      
      if(apply_region_ids) then
        deallocate(region_ids)
        if (have_option(trim(option_path)//"/from_mesh/extrude/regions["//&
                                int2str(r)//"]/sizing_function/list")) then
          deallocate(sizing_vector)
        end if
      end if
    end do

    if (have_option(trim(option_path)//'/from_mesh/extrude/regions['//int2str(r)//&
                                         ']/bottom_depth/from_map/')) then

      deallocate (depth_vector)

    end if

#ifdef DDEBUG
    if(apply_region_ids) then
      assert(minval(visited)>0)
      ewrite(2,*) "Maximum number of times a node was visited: ", maxval(visited)
    end if
#endif

    ! Now the tiresome business of making a shape function.
    shell_dim = mesh_dim(shell_mesh)
    call get_option("/geometry/quadrature/degree", quadrature_degree)
    quad = make_quadrature(vertices=shell_dim+2, dim=shell_dim+1, degree=quadrature_degree)
    full_shape = make_element_shape(vertices=shell_dim+2, dim=shell_dim+1, degree=1, quad=quad)
    call deallocate(quad)

    call get_option(trim(option_path)//'/name', mesh_name)

    ! combine the 1d radial meshes into a full mesh
    call combine_r_meshes(shell_mesh, r_meshes, out_mesh, &
       full_shape, mesh_name, option_path)

!     vtk_write_state(filename, index, model, state, write_region_ids, stat)
!     call vtk_write_fields("extruded_mesh", 0, out_mesh, out_mesh%mesh, vfields=(/out_mesh/))
!     call vtk_write_surface_mesh("surface_mesh", 1, out_mesh)
!     call write_triangle_files("test_mesh", out_mesh)
       
    do column=1, node_count(shell_mesh)
      if (.not. node_owned(shell_mesh, column)) cycle
      call deallocate(r_meshes(column))
    end do
    call deallocate(full_shape)
        
  end subroutine extrude_radially

  subroutine compute_r_nodes(r_mesh, depth, xyz, r_shell, sizing, sizing_function)
    !!< Figure out at what depths to put the layers.
    type(vector_field), intent(out) :: r_mesh
    real, intent(in):: depth
    real, intent(out):: r_shell
    real, dimension(:), intent(in):: xyz
    real, optional, intent(in):: sizing
    character(len=*), optional, intent(in):: sizing_function

    ! this is a safety gap:
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
    real, dimension(1:size(xyz)):: xyz_new
    real :: delta_r, r, rnew
    real :: thetanew, phinew
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
      FLAbort("Need to supply either sizing, sizing_function")
    end if

    ! First work out number of nodes/elements:
    r_shell=norm2(xyz)
    r=r_shell
    node=2
    xyz_new(1:size(xyz))=xyz
    do
      delta_r = get_delta_r( xyz_new, is_constant, constant_value, py_func)
      r=r - delta_r
      if (r<(r_shell-depth+MIN_BOTTOM_LAYER_FRAC*delta_r)) exit
      node=node+1
      if (node>MAX_VERTICAL_NODES) then
        ewrite(-1,*) "Check your extrude/sizing_function"
        FLAbort("Maximum number of vertical layers reached")
      end if
    end do
    elements=node-1

    call allocate(mesh, elements+1, elements, oned_shape, "RMesh") ! (mesh, nodes, elements, shape, name)
    call deallocate(oned_shape)
    do ele=1,elements
      mesh%ndglno((ele-1) * loc + 1: ele*loc) = (/ele, ele+1/)
    end do

    call allocate(r_mesh, 3, mesh, "RMeshCoordinates")
    call deallocate(mesh)

    ! Start the mesh at r=r(shell_r) and work down to r=shell_r-depth.
    call set(r_mesh, 1, xyz)
    ! Set theta and phi - these will be the same at all layers
    thetanew=acos(xyz(3)/norm2(xyz))
    phinew=atan2(xyz(2),xyz(1))  
    do node=2,elements
      xyz_new=node_val(r_mesh, node-1)
      delta_r = get_delta_r(xyz_new, is_constant, constant_value, py_func)
      rnew=norm2(xyz_new)-delta_r
      xyz_new(1)=rnew*sin(thetanew)*cos(phinew)
      xyz_new(2)=rnew*sin(thetanew)*sin(phinew)
      xyz_new(3)=rnew*cos(thetanew)
      call set(r_mesh, node, xyz_new )
    end do    
    xyz_new(1)=(r_shell-depth)*sin(thetanew)*cos(phinew)
    xyz_new(2)=(r_shell-depth)*sin(thetanew)*sin(phinew)
    xyz_new(3)=(r_shell-depth)*cos(thetanew)
    call set(r_mesh, elements+1, xyz_new)

    ! For pathological sizing functions the mesh might have gotten inverted at the last step.
    ! If you encounter this, make this logic smarter.
!    assert(all(node_val(r_mesh, elements) > node_val(r_mesh, elements+1)))  -- check this for radial extrushions
    
    assert(oned_quad%refcount%count == 1)
    assert(oned_shape%refcount%count == 1)
    assert(r_mesh%refcount%count == 1)
    assert(mesh%refcount%count == 1)

    contains
    
      function get_delta_r(pos, is_constant, constant_value, py_func) result(delta_r)
        real, dimension(:), intent(in) :: pos
        logical, intent(in) :: is_constant
        real, intent(in) :: constant_value
        character(len=PYTHON_FUNC_LEN), intent(in) :: py_func

        real :: delta_r
        real, dimension(1) :: delta_r_tmp
        real, dimension(size(pos), 1) :: pos_tmp
        
        if (is_constant) then
          delta_r = constant_value
        else
          pos_tmp(:, 1) = pos
          call set_from_python_function(delta_r_tmp, trim(py_func), pos_tmp, time=0.0)
          delta_r = delta_r_tmp(1)
        end if
        assert(delta_r > 0.0)
        
      end function get_delta_r
      
  end subroutine compute_r_nodes

    
end module hadapt_extrude_radially
