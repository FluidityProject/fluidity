  program LocalMats
    ! A small program for constructing local matrices
    ! Useful for doing linear analyses on mixed elements etc.
    use fields
    use FEtools
    use DGtools
    use elements
    use transform_elements
    use signal_vars
    !use global_parameters, only : current_debug_level, PYTHON_FUNC_LEN
    !use fldebug

    implicit none
    integer :: degree, quad_degree=6
    type(quadrature_type), target :: quad,f_quad
    type(element_type), target :: X_shape, u_shape, h_shape
    real, allocatable, dimension(:,:) :: X_ele
    character(len=100) :: tmpbuf
    ! Arguments for handling the command line
    character(len=256) :: filename
    integer :: status, ele
    character(len=100) :: fmt,buffer,buf
    integer :: u_degree=1, h_degree=2, dim=1, simplex_type=1

    !simplex types are:
    !1. Equilateral, unit. 2. Right-angled, unit.

    namelist/LocalMats_data/u_degree,h_degree,dim,simplex_type, quad_degree
    logical :: file_exists
    integer :: unit, io1, stat

    !current_debug_level = 2

    print *,  'program LocalMats'

    call get_command_argument(1, value=filename, status=status)
    select case(status)
    case(1:)
       print *,  'no file specified?'
       stop
    case(:-1)
       write(0,*) "Warning: truncating filename"
    end select
    filename=trim(filename)

    print *,  'Reading parameters'

    unit = free_unit()
    open(unit=unit, file=trim(trim(filename)//".dat"), status='old', &
         iostat=io1)

    if(io1.ne.0) then
       print *, 'Looked for ', trim(trim(filename)//".dat")
       print *, 'Could not read from .dat file'
       stop
    end if

    read(unit, LocalMats_data)
    close(unit) 

    print *,  'Making local matrices for an element with'
    print *,  'Dimension == ', dim
    print *,  'U degree = ', u_degree
    print *,  'H degree = ', h_degree
    print *,  'simplex type = ',simplex_type 
    print *,  'Quad degree =', quad_degree

    print *,  'Getting quadrature'

    quad=make_quadrature(loc=dim+1, dimension=dim, degree=quad_degree)

    print *,  'Getting shape functions'

    ! Shape functions for positions (linear)
    X_shape=make_element_shape(loc=dim+1, dimension=dim, &
         degree=1, quad=quad)

    ! Shape functions for velocity and height
    u_shape=make_element_shape(loc=dim+1, &
         dimension=dim, degree=u_degree, quad=quad)
    h_shape=make_element_shape(loc=dim+1, &
         dimension=dim, degree=h_degree, quad=quad)

    allocate(X_ele(u_shape%dim,x_shape%loc))

    call get_X_ele(X_ele,simplex_type)

    print *,  'Coordinates'
    print *,  X_ele

    call assemble_local_matrices(u_shape,h_shape,X_shape,X_ele)

  contains

    subroutine assemble_local_matrices( &
         u_shape,h_shape,X_shape,X_ele)
      type(element_type), intent(in) :: u_shape,h_shape,X_shape
      real, dimension(u_shape%dim,x_shape%loc), intent(in) :: X_ele

      !! Local Stuff
      ! Coordinate transform * quadrature weights.
      real, dimension(u_shape%ngi) :: detwei
      ! Derivatives of shape function:
      real, dimension(h_shape%loc, h_shape%ngi, h_shape%dim) :: h_dshape
      ! gradient matrix
      real, dimension(u_shape%dim,h_shape%loc,u_shape%loc) :: grad_mat
      real, dimension(u_shape%loc,u_shape%loc) :: u_mass_mat
      real, dimension(h_shape%loc,h_shape%loc) :: lap_mat,h_mass_mat
      integer :: idim, iloc

      ! Transform derivatives and weights into physical space.
      call transform_to_physical(X_ele, X_shape, m=h_shape, &
           dm_t=h_dshape, detwei=detwei)

      grad_mat = -dshape_shape(h_dshape,u_shape,detwei)
      lap_mat = dshape_dot_dshape(h_dshape,h_dshape,detwei)
      u_mass_mat = shape_shape(u_shape,u_shape,detwei)
      h_mass_mat = shape_shape(h_shape,h_shape,detwei)

      print *,  'U Mass Matrix'

      do iloc = 1, u_shape%loc
         print *,  u_mass_mat(iloc,:)
      end do

      print *,  'H mass matrix'
      do iloc = 1, h_shape%loc
         print *,  h_mass_mat(iloc,:)
      end do

      print *,  'Div matrix'
      do idim = 1, u_shape%dim
         do iloc = 1, h_shape%loc
            print *,  grad_mat(idim,iloc,:)
         end do
         print *,  ' '
      end do

      print *,  'Laplacian matrix'
      do iloc = 1, h_shape%loc
         print *, lap_mat(iloc,:)
      end do

    end subroutine assemble_local_matrices

    subroutine get_X_ele(X_ele,simplex_type)
      real, dimension(:,:), intent(inout) :: X_ele
      integer, intent(in) :: simplex_type

      !
      select case(size(X_ele,1))
      case (1)
         X_ele(1,:) = (/ 0.0, 1.0 /)
      case (2)
         select case(simplex_type)
         case (1)
            X_ele(1,:) = (/ 0.0, 1.0, 0.5 /)
            X_ele(2,:) = (/ 0.0, 0.0, 0.8660254037844386 /)
         case (2)
            X_ele(1,:) = (/ 0.0, 1.0, 1.0 /)
            X_ele(2,:) = (/ 0.0, 0.0, 1.0 /)
         case default
            print *, 'no support for simplex type, try coding it?'
            stop
         end select
      case default
         print *, 'no support for dimension, try coding it?'
         stop
      end select

    end subroutine get_X_ele
end program LOCALMATS
