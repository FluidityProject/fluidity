#include "confdefs.h"
program visualise_elements
  use quadrature
  use element_numbering, only: tr, ELEMENT_LAGRANGIAN
  use elements
  use spud
  use fields
  use state_module
  use vtk_interfaces
  use populate_state_module
  use field_derivatives
  implicit none

  type(element_type) :: element, visualisation_element
  type(vector_field) :: position, visualisation_position, outline_position
  type(scalar_field) :: shape_values, visualisation_shape_values, tracer
  type(mesh_type) :: linear_mesh, visualisation_mesh
  type(state_type), dimension(:), pointer :: state
  type(vector_field) :: derivative, visualisation_derivative
  type(scalar_field), dimension(2) :: dx
  character(len=1024) :: projectname
#ifdef HAVE_MPI
  integer :: ierror, stat

  call set_global_debug_level(0)

  call mpi_init(ierror)
#endif

  call python_init()
  call read_command_line(projectname)

  call set_option("/timestepping/current_time", 0.0, stat=stat)

  if (have_option("/geometry/mesh")) then
     call populate_state(state)
     
     position=extract_vector_field(state(1),"Coordinate")

     tracer=extract_scalar_field(state(1), "Tracer")

     element=construct_element()
     visualisation_element=construct_visualisation_element()
     visualisation_mesh=make_mesh(tracer%mesh, visualisation_element,&
          & continuity=-1)
     call allocate(visualisation_position, 2, visualisation_mesh, "Coordinate")
     call remap_field(position, visualisation_position)
     call allocate(visualisation_shape_values, visualisation_mesh, &
          "Tracer" )
     call remap_field(tracer, visualisation_shape_values)

     linear_mesh=subdivide_elements(visualisation_position%mesh)
     
     visualisation_shape_values%mesh=linear_mesh
     visualisation_position%mesh=linear_mesh

     outline_position=construct_outline_positions(element)

     if (has_vector_field(state(1),"Derivative")) then
        
        derivative=extract_vector_field(state(1), "Derivative")
        
        dx(1)=extract_scalar_field(derivative, 1)
        dx(2)=extract_scalar_field(derivative, 2)
     
        print *, dx(1)%val

        call differentiate_field(tracer, position, (/.true.,.true./), dx)
        
        print *, tracer%val
        print *, dx(1)%val

        call allocate(visualisation_derivative, 2,visualisation_mesh, &
             "Derivative" )
        call remap_field(derivative, visualisation_derivative)
        
        call vtk_write_fields(projectname, 0, &
             visualisation_position, &
             linear_mesh, &
             sfields=(/visualisation_shape_values/), &
             vfields=(/visualisation_derivative/))

        call vtk_write_fields("linear"//projectname, 0, &
             position, &
             derivative%mesh, vfields=(/derivative/))

     else
        call vtk_write_fields(projectname, 0, &
             visualisation_position, &
             linear_mesh, &
             sfields=(/visualisation_shape_values/))

        call vtk_write_fields("outline"//projectname, 0, &
             outline_position, &
             outline_position%mesh)


     end if

  else
     
     element=construct_element()
     position=construct_positions(element)
     shape_values=construct_shape_values(position%mesh)

     visualisation_element=construct_visualisation_element()
     visualisation_mesh=make_mesh(position%mesh, visualisation_element)

     call allocate(visualisation_position, 2, visualisation_mesh, "Coordinate")
     call remap_field(position, visualisation_position)
     call allocate(visualisation_shape_values, visualisation_mesh, &
          "ShapeValues" )
     call remap_field(shape_values, visualisation_shape_values)

     linear_mesh=subdivide_elements(visualisation_position%mesh)

     ! Note memory leak
     visualisation_shape_values%mesh=linear_mesh
     visualisation_position%mesh=linear_mesh

     outline_position=construct_outline_positions(element)

     call vtk_write_fields(projectname, 0, &
          visualisation_position, &
          linear_mesh, &
          sfields=(/visualisation_shape_values/))

     call vtk_write_fields("outline"//projectname, 0, &
          outline_position, &
          outline_position%mesh)

  end if

contains

  function subdivide_elements(mesh) result (linear_mesh)
    type(mesh_type), intent(in) :: mesh
    type(mesh_type) :: linear_mesh
    
    type(element_type) :: element, linear_element
    integer :: triangles, e, n, row, column, rowlen, ele

    element=mesh%shape

    linear_element=make_element_shape(element%numbering%vertices, 2, degree=1, &
         quad=element%quadrature)
    
    triangles=tr(element%degree) + tr(element%degree-1)

    call allocate(linear_mesh, nodes=node_count(mesh), &
         elements=triangles*element_count(mesh), shape=linear_element, &
         name="Linear"//mesh%name)

    e=0
    do ele=1, element_count(mesh)
       ! Point up triangles.
       n=element%loc*(ele-1)
       do row=1, element%degree
          rowlen=element%degree+2-row
          do column=1,element%degree+1-row
             n=n+1
             e=e+1
             linear_mesh%ndglno((e-1)*3+1:e*3)&
                  =(/n, n+1, n+rowlen/)   
          end do
          n=n+1
       end do
       
       ! Point down triangles.
       n=element%loc*(ele-1)+1
       do row=1, element%degree-1
          rowlen=element%degree+2-row
          do column=1,element%degree-row
             n=n+1
             e=e+1
             linear_mesh%ndglno((e-1)*3+1:e*3)&
                  =(/n, n+rowlen, n+rowlen-1/)   
          end do
          n=n+2
       end do
    end do

  end function subdivide_elements

  function construct_shape_values(mesh) result (sfield)
    type(scalar_field) :: sfield
    type(mesh_type) :: mesh
    
    integer :: ele, node
    integer, dimension(:), pointer :: nodes

    call allocate(sfield, mesh, "ShapeValues")

    call zero(sfield)

    do ele=1,element_count(mesh)
       nodes=>ele_nodes(mesh,ele)
       do node=1,size(nodes)
          if (node==ele) then
             call set(sfield, nodes(node), 1.0)
          end if
       end do
    end do

  end function construct_shape_values

  function construct_element() result(element)
    type(element_type) :: element
    type(quadrature_type) :: quadrature

    integer :: degree, dim, vertices, type
    character(len=1024) :: family

    call get_option("/geometry/element_degree", degree)
    call get_option("/geometry/dimension", dim)
    call get_option("/geometry/element_vertices", vertices)
    call get_option("/geometry/element_family", family)

    select case(family)
    case("lagrange")
       type=ELEMENT_LAGRANGIAN
    case("default")
       write(0,*) "Unknown element type "//trim(family)
    end select

    quadrature = make_quadrature(vertices, dim, degree)

    element = make_element_shape(vertices, dim, degree, quadrature, type)

  end function construct_element

  function construct_visualisation_element() result(element)
    type(element_type) :: element
    type(quadrature_type) :: quadrature

    integer :: degree, dim, vertices, type
    character(len=1024) :: family

    call get_option("/geometry/visualisation_degree", degree)
    call get_option("/geometry/dimension", dim)
    call get_option("/geometry/element_vertices", vertices)
    call get_option("/geometry/element_family", family)

    select case(family)
    case("lagrange")
       type=ELEMENT_LAGRANGIAN
    case("default")
       write(0,*) "Unknown element type "//trim(family)
    end select

    ! Don't care about quadrature as we won't do any calculus
    quadrature = make_quadrature(vertices, dim, degree=1)

    element = make_element_shape(vertices, dim, degree, quadrature, type)

  end function construct_visualisation_element

  function construct_positions(element) result (position)
    type(element_type), intent(inout) :: element
    type(vector_field) :: position

    type(mesh_type) :: mesh, linear_mesh, one_mesh_linear, one_mesh
    type(vector_field) :: linear_position, one_element_linear, one_element
    type(element_type) :: linear_element

    real,dimension(2,element%numbering%vertices) :: vertices, lvertices
    real,dimension(2) :: node_loc
    real :: scale
    integer :: i, d

    vertices=regular_figure(element%numbering%vertices, 1.0)

    linear_element=make_element_shape(element%numbering%vertices, 2, degree=1, &
         quad=element%quadrature)

    linear_mesh=construct_mesh(linear_element, element, "LinearMesh")

    one_mesh_linear=construct_one_element_mesh(linear_element, &
         "OneLinearMesh")

    one_mesh=construct_one_element_mesh(element, "OneMesh")
    
    call allocate(one_element_linear, 2, one_mesh_linear, &
         "OneLinearCoordinate")

    call allocate(one_element, 2, one_mesh, "OneCoordinate")

    call set(one_element_linear, ele_nodes(one_element_linear,1), vertices)
    
    call remap_field(one_element_linear, one_element)

    call allocate(linear_position, 2, linear_mesh, "LinearCoordinate")
    
    scale=1.5*element%degree

    do i=1, node_count(one_element)
       node_loc=node_val(one_element, i)
       do d=1,2
          lvertices(d,:)=vertices(d,:)+scale*node_loc(d)
       end do
       
       call set(linear_position, ele_nodes(linear_mesh,i), lvertices)
    end do
    
    mesh = construct_mesh(element, element, "Mesh")
    call allocate(position, 2, mesh, "Coordinate")
    
    call remap_field(linear_position, position)

    call deallocate(linear_position)
    call deallocate(linear_mesh)
    call deallocate(linear_element)
    call deallocate(one_element)
    call deallocate(one_element_linear)
    call deallocate(one_mesh)
    call deallocate(one_mesh_linear)

  end function construct_positions

  function construct_outline_positions(element) result (position)
    type(element_type), intent(inout) :: element
    type(vector_field) :: position

    type(mesh_type) :: mesh, linear_mesh, one_mesh_linear, one_mesh
    type(vector_field) :: linear_position, one_element_linear, one_element
    type(element_type) :: linear_element

    real,dimension(2,element%numbering%vertices) :: vertices, lvertices
    real,dimension(2) :: node_loc
    real :: scale
    integer :: i, d

    vertices=regular_figure(element%numbering%vertices, 1.0)

    linear_element=make_element_shape(element%numbering%vertices, 2, degree=1, &
         quad=element%quadrature)

    linear_mesh=construct_mesh(linear_element, element, "LinearMesh")

    one_mesh_linear=construct_one_element_mesh(linear_element, &
         "OneLinearMesh")

    one_mesh=construct_one_element_mesh(element, "OneMesh")
    
    call allocate(one_element_linear, 2, one_mesh_linear, &
         "OneLinearCoordinate")

    call allocate(one_element, 2, one_mesh, "OneCoordinate")

    call set(one_element_linear, ele_nodes(one_element_linear,1), vertices)
    
    call remap_field(one_element_linear, one_element)

    call allocate(position, 2, linear_mesh, "LinearCoordinate")
    
    scale=1.5*element%degree

    do i=1, node_count(one_element)
       node_loc=node_val(one_element, i)
       do d=1,2
          lvertices(d,:)=vertices(d,:)+scale*node_loc(d)
       end do
       
       call set(position, ele_nodes(linear_mesh,i), lvertices)
    end do
    
    call deallocate(linear_mesh)
    call deallocate(linear_element)
    call deallocate(one_element)
    call deallocate(one_element_linear)
    call deallocate(one_mesh)
    call deallocate(one_mesh_linear)

  end function construct_outline_positions

  function construct_mesh(mesh_element, layout_element, name) result(mesh)
    type(mesh_type) :: mesh
    type(element_type), intent(inout) :: mesh_element
    type(element_type), intent(in) :: layout_element
    character(len=*), intent(in) :: name

    integer :: i

    call allocate(mesh, nodes=mesh_element%loc*layout_element%loc, &
         elements=layout_element%loc, shape=mesh_element, name=name)
    ! Definitely a DG mesh
    mesh%continuity=-1

    ! Usual dg numbering.
    mesh%ndglno=(/(i, i=1,mesh_element%loc*layout_element%loc)/)

  end function construct_mesh

  function construct_one_element_mesh(element, name) result(mesh)
    type(mesh_type) :: mesh
    type(element_type), intent(inout) :: element
    character(len=*), intent(in) :: name

    integer :: i

    call allocate(mesh, nodes=element%loc, elements=1, shape&
         &=element, name=name)
    
    ! Usual dg numbering.
    mesh%ndglno=(/(i, i=1,element%loc)/)

  end function construct_one_element_mesh

  subroutine read_command_line(projectname)
    ! Read the input filename, degree and quadrature degree on the command
    ! line.
    character(len=*), intent(out) :: projectname

    character(len=1024) :: filename
    integer :: status
    
    call get_command_argument(1, value=filename, status=status)
  
    if (status/=0) then
       call usage
       stop
    end if

    call load_options(filename)

    call get_option("/project_name", projectname)

  end subroutine read_command_line

  subroutine usage
    
    write (0,*) "usage: visualise_elements <options_file>"
    
  end subroutine usage

  function regular_figure(nodes, length) result (vertex)
    ! Return the locations of the vertices of the interval, square or
    ! triangle centred at the origin with side length length.
    integer, intent(in) :: nodes
    real, intent(in) :: length
    real, dimension(2, nodes) :: vertex

    real :: height

    select case (nodes)
    case(2)
       vertex(2,:)=0.0
       vertex(1,:)=(/-length/2,length/2/)
    case(3)
       height=sqrt(0.75)*length
       vertex(:,1)=(/-0.5*length, -1./3. * height/)
       vertex(:,2)=(/0.5*length, -1./3. * height/)
       vertex(:,3)=(/0.0, 2./3.*height/)
    case(4)
       vertex(:,1)=(/-length/2,-length/2/)
       vertex(:,2)=(/length/2,-length/2/)
       vertex(:,3)=(/-length/2,length/2/)
       vertex(:,4)=(/length/2,length/2/)
    case default
       write(0,*) "Illegal number of nodes"
    end select
  end function regular_figure

end program visualise_elements
