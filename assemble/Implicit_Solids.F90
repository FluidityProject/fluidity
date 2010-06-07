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
!    C.Pain@Imperial.ac.uk
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

module implicit_solids
! these 5 need to be on top and in this order, 
! so as not to confuse silly old intel compiler 
  use quadrature
  use elements
  use sparse_tools
  use fields
  use state_module
!
  use vtk_interfaces
  use linked_lists
  use intersection_finder_module
  use tetrahedron_intersection_module
  use unify_meshes_module
  use unittest_tools
  use global_parameters, only: FIELD_NAME_LEN, current_time, dt
  use spud
  use read_triangle
  use timeloop_utilities
  use fefields, only: compute_lumped_mass
  use parallel_tools

  implicit none

  private
  public:: solids

contains

  subroutine solids(state, its)

    type(state_type),intent(inout) :: state
    character(len=field_name_len) :: external_mesh_name
    type(vector_field), pointer :: positions, velocity
    type(vector_field), pointer :: absorption
    type(vector_field), save :: external_positions
    type(scalar_field), pointer :: solid
    type(scalar_field) :: lumped_mass
    integer :: ele_A, ele_B, ele_C
    type(tet_type) :: tet_A, tet_B
    type(plane_type), dimension(:), allocatable :: planes_A
    integer :: stat, nintersections, i, j, k
    integer :: ntests, dat_unit, its
    integer, dimension(:), pointer :: ele_A_nodes
    type(vector_field) :: intersection
    real, dimension(:,:), allocatable :: pos_A
    real, dimension(:), allocatable :: detwei, drag
    real :: vol, ele_A_vol, sigma
    real, save :: beta
    integer, save :: itinoi
    logical, save :: init=.false.

    ewrite(3, *) "inside femdem"

    solid => extract_scalar_field(state, "SolidConcentration")
    call zero(solid)

    positions => extract_vector_field(state, "Coordinate")

    if (.not. init) then
       call get_option("/implicit_solids/mesh_name", external_mesh_name)
       external_positions = &
            read_triangle_serial(trim(external_mesh_name), quad_degree=1)

       call get_option("/implicit_solids/beta", beta, default=1.)

       assert(positions%dim >= 2)
       assert(positions%dim == external_positions%dim)

       if (GetRank() == 0) then
          dat_unit = free_unit()
          open(dat_unit, file="drag_force", status="replace")
          close(dat_unit)                  
       end if

       call get_option("/timestepping/nonlinear_iterations", itinoi)

       init=.true.
    end if

    call rtree_intersection_finder_set_input(positions)

    do ele_B = 1, ele_count(external_positions)

       call rtree_intersection_finder_find(external_positions, ele_B)
       call rtree_intersection_finder_query_output(nintersections)

       if (positions%dim == 3) then
          tet_B%v = ele_val(external_positions, ele_B)
       else
          call intersector_set_dimension(positions%dim)
       end if

       do j = 1, nintersections

          call rtree_intersection_finder_get_output(ele_A, j)

          if (positions%dim == 3) then

             if (ele_loc(positions, ele_A)==4) then
                ! tets
                allocate(planes_A(4))
                tet_A%v = ele_val(positions, ele_A)
                planes_A = get_planes(tet_A)
             else
                ! hexes
                allocate(planes_A(6))
                planes_A = get_planes(positions, ele_A)
             end if

             call intersect_tets(tet_B, planes_A, &
                  ele_shape(external_positions, ele_B), &
                  stat=stat, output=intersection)

             deallocate(planes_A)

          else

             allocate(pos_A(positions%dim, ele_loc(positions, ele_A)))
             pos_A = ele_val(positions, ele_A)
             intersection = intersect_elements(external_positions, &
                  ele_B, pos_A, ele_shape(external_positions, ele_B))
             deallocate(pos_A)
             stat = 0

          end if

          if (stat == 1) cycle

          vol = 0.
          do ele_C = 1, ele_count(intersection)
             vol = vol + abs(simplex_volume(intersection, ele_C))
          end do

          allocate(detwei(ele_ngi(positions, ele_A)))
          call transform_to_physical(positions, ele_A, detwei=detwei)
          ele_A_vol = sum(detwei)
          deallocate(detwei)

          ele_A_nodes => ele_nodes(positions, ele_A)
          do k = 1, size(ele_A_nodes)
             call addto(solid, ele_A_nodes(k), vol/ele_A_vol)
          end do

          call deallocate(intersection)

       end do

    end do

    do i = 1, node_count(solid)
       call set(solid, i, max(0., min(1., node_val(solid, i))))
    end do

    ewrite_minmax(solid)

    absorption => extract_vector_field(state, "VelocityAbsorption")
    call zero(absorption)

    do i = 1, node_count(absorption)
       sigma = node_val(solid, i)*beta/dt
       do j = 1, absorption%dim
          call set(absorption, j, i, sigma)
       end do
    end do

    ! total solid drag force...
    velocity => extract_vector_field(state, "Velocity")

    call allocate(lumped_mass, positions%mesh, "Lumped mass")
    call compute_lumped_mass(positions, lumped_mass)    

    allocate(drag(positions%dim))
    drag = 0.

    do i = 1, node_count(positions)
       do j = 1, positions%dim
          drag(j) = drag(j) + &
               node_val(absorption, 1, i) * node_val(velocity, j, i) * &
               node_val(lumped_mass, i)
       end do
    end do

    call deallocate(lumped_mass)

    call allsumv(drag)

    if (GetRank() == 0 .and. its == itinoi) then
       open(1453, file="drag_force", position="append")
       write(1453, *) &
            (drag(i), i = 1, positions%dim)
       close(1453)
    end if

    deallocate(drag)

    if(simulation_completed(current_time)) call deallocate(external_positions)
    call finalise_tet_intersector
    call rtree_intersection_finder_reset(ntests)

    ewrite(3, *) "leaving femdem"

  end subroutine solids

  !----------------------------------------------------------------------------

  function read_triangle_serial(filename, quad_degree) result (field)
    
    character(len=*), intent(in) :: filename
    !! The degree of the quadrature.
    integer, intent(in), optional, target :: quad_degree
    !! The degree of the quadrature.

    type(vector_field) :: field
    type(quadrature_type) :: quad
    type(element_type) :: shape

    integer :: dim, loc

    call identify_triangle_file(filename, dim, loc)
    quad=make_quadrature(loc, dim, degree=quad_degree)
    shape=make_element_shape(loc, dim, 1, quad)
    field=read_triangle_files_serial(filename, shape)

    ! deallocate our references of shape and quadrature:
    ! NOTE: we're using the specific deallocate interface here
    !       to make the intel compiler shut up
    call deallocate_element(shape)
    call deallocate(quad)

  end function read_triangle_serial


  function read_triangle_files_serial(filename, shape) result (field)
    !!< Filename is the base name of the triangle file without .node or .ele.

    character(len=*), intent(in) :: filename
    type(element_type), intent(in), target :: shape
    type(vector_field) :: field

    integer :: node_unit, ele_unit
    real, allocatable, dimension(:) :: read_buffer
    integer, allocatable, dimension(:,:) :: edge_buffer
    integer, allocatable, dimension(:) :: sndglno
    integer, allocatable, dimension(:) :: boundary_ids, element_owner

    character(len = parallel_filename_len(filename)) :: lfilename
    integer :: i, j, nodes, dim, node_attributes, boundaries, &
         ele_attributes, loc, sloc, elements, edges, edge_count
    integer, allocatable, dimension(:):: node_order
    logical :: file_exists
    type(mesh_type) :: mesh

    lfilename = trim(filename)

    node_unit=free_unit()

    ewrite(2, *) "Opening "//trim(lfilename)//".node for reading."
    ! Open node file
    open(unit=node_unit, file=trim(lfilename)//".node", err=42, action="read")

    ! Read node file header.
    read (node_unit, *) nodes, dim, node_attributes, boundaries

    ele_unit=free_unit()

    ewrite(2, *) "Opening "//trim(lfilename)//".ele for reading."
    ! Open element file
    open(unit=ele_unit, file=trim(lfilename)//".ele", err=43, action="read")

    ! Read element file header.
    read (ele_unit, *) elements, loc, ele_attributes

    assert(loc==shape%loc)
    allocate(node_order(loc))
    select case(loc)
    case(3)
       node_order = (/1,2,3/)
    case default
       do j = 1, loc
          node_order(j) = j
       end do
    end select

    call allocate(mesh, nodes, elements, shape, name="CoordinateMesh")
    call allocate(field, dim, mesh, name="Coordinate")

    ! Drop the local reference to mesh - now field owns the only reference.
    call deallocate(mesh)

    allocate(read_buffer(dim+node_attributes+boundaries+1))

    if(node_attributes==1) then ! this assumes the node attribute are column numbers
       allocate(field%mesh%columns(1:nodes))
    end if

    do i = 1, nodes
       read(node_unit,*) read_buffer
       forall (j=1:dim)
          field%val(j)%ptr(i)=read_buffer(j+1)
       end forall
       if (node_attributes==1) then
          field%mesh%columns(i)=floor(read_buffer(dim+1))
       end if
    end do

    deallocate(read_buffer)
    allocate(read_buffer(loc+ele_attributes+1))

    if(ele_attributes==1) then  ! this assumes that the element attribute is a region id
       allocate(field%mesh%region_ids(1:elements))
    end if

    do i = 1, elements
       read(ele_unit,*) read_buffer
       field%mesh%ndglno((i-1)*loc+1:i*loc)=floor(read_buffer(node_order+1))
       if(ele_attributes==1) then
          field%mesh%region_ids(i)=read_buffer(loc+2)
       end if
    end do

    close(node_unit)
    close(ele_unit)

    ! Open edge file
    select case (dim)
    case(2)
       inquire(file=trim(lfilename)//".edge",exist=file_exists)
       if(file_exists) then
          ewrite(2, *) "Opening "//trim(lfilename)//".edge for reading."
          open(unit=node_unit, file=trim(lfilename)//".edge", err=41, &
               action="read")
       end if
    case(3)
       inquire(file=trim(lfilename)//".face",exist=file_exists)
       if(file_exists) then
          ewrite(2, *) "Opening "//trim(lfilename)//".face for reading."
          open(unit=node_unit, file=trim(lfilename)//".face", err=41, &
               action="read")
       end if
    end select

    if(file_exists) then
       ! Read edge file header.
       read (node_unit, *) edges, boundaries
    else
       edges = 0
       boundaries = 1
    end if

    if(edges==0) then
       file_exists = .false.
       close(node_unit)
    end if

    select case(shape%numbering%family)
    case(FAMILY_SIMPLEX)
       if ((loc/=dim+1).and.(boundaries/=0)) then
          ewrite(0,*) "Warning: triangle boundary markers not supported for qua", &
               "dratic space elements."
          if(file_exists) then
             file_exists= .false.
             close(node_unit)
          end if
       end if
       sloc=loc-1
    case default
       FLAbort('Illegal element family')
    end select

    allocate(edge_buffer(sloc+boundaries+1,edges))
    edge_buffer=0
    allocate(sndglno(edges*sloc))
    sndglno=0
    allocate(boundary_ids(1:edges))
    boundary_ids=0
    if (boundaries==2) then
       allocate(element_owner(1:edges))
       element_owner=0
    end if
    edge_count=0

    if (boundaries==0) then
       ewrite(0,*) "Warning: triangle edge file has no boundary markers"
       if(file_exists) then
          file_exists=.false.
          close(node_unit)
       end if
    else
       if(file_exists) then
          read(node_unit, *) edge_buffer
          do i = 1, edges
             if (edge_buffer(sloc+2,i)/=0) then
                ! boundary edge/face
                edge_count=edge_count+1
                sndglno((edge_count-1)*sloc+1:edge_count*sloc)= &
                     edge_buffer(2:sloc+1,i)
                boundary_ids(edge_count)=edge_buffer(sloc+2,i)
                if (boundaries==2) then
                   element_owner(edge_count)=edge_buffer(sloc+3,i)
                end if
             end if
          end do

          file_exists=.false.
          close(node_unit)
       end if
    end if

    if (boundaries<2) then
       call add_faces(field%mesh, &
            sndgln=sndglno(1:edge_count*sloc), &
            boundary_ids=boundary_ids(1:edge_count))
    else
       call add_faces(field%mesh, &
            sndgln=sndglno(1:edge_count*sloc), &
            boundary_ids=boundary_ids(1:edge_count), &
            element_owner=element_owner)
    end if

    deallocate(edge_buffer)
    deallocate(sndglno)
    deallocate(boundary_ids)

41  continue ! We jump to here if there was no edge file.

    return

42  FLAbort("Unable to open "//trim(lfilename)//".node")

43  FLAbort("Unable to open "//trim(lfilename)//".ele")

  end function read_triangle_files_serial

end module implicit_solids
