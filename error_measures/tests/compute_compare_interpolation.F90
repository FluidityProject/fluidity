subroutine compute_compare_interpolation

  use fields
  use mesh_files
  use conservative_interpolation_module
  use unittest_tools
  use interpolation_module
  use supermesh_construction
  use vtk_interfaces
  use futils
  use vector_tools
  use solvers
  implicit none

  interface
    function field_func_const(pos)
      real, dimension(:) :: pos
      real :: solution
    end function
  end interface

  interface
    function field_func_linear(pos)
      real, dimension(:) :: pos
      real :: solution
    end function
  end interface

  interface
    function field_func_quadratic(pos)
      real, dimension(:) :: pos
      real :: solution
    end function
  end interface

  interface
    function field_func_cubic(pos)
      real, dimension(:) :: pos
      real :: solution
    end function
  end interface

  interface
    function field_func_tophat(pos)
      real, dimension(:) :: pos
      real :: solution
    end function
  end interface

  interface
    function field_func_peaks(pos)
      real, dimension(:) :: pos
      real :: solution
    end function
  end interface

  type(vector_field) :: positionsA, positionsB
  type(scalar_field), dimension(2) :: c_fieldA, c_fieldB
  type(scalar_field), dimension(2) :: l_fieldA, l_fieldB
  type(scalar_field), dimension(2) :: analytical_fields
  real, dimension(2) :: c_integral, l_integral
  real, dimension(2) :: c_maxval, c_minval, l_maxval, l_minval
  real, dimension(2) :: c_l2err, l_l2err
  integer :: field, field_count
  integer :: i

  positionsA = read_mesh_files("data/input.1", quad_degree=4, format="gmsh")

  field_count = 2

  do field=1,field_count
    call allocate(c_fieldA(field), positionsA%mesh, "Field" // int2str(field))
    c_fieldA(field)%option_path = "/c_fieldA" // int2str(field) // "/prognostic/conservative_interpolation_galerkin"
    call set_solver_options(c_fieldA(field), ksptype='cg', pctype='eisenstat', rtol=1.0e-7, max_its=10000)
    c_fieldA(field)%option_path = "/c_fieldA" // int2str(field) // "/prognostic/conservative_interpolation_bounded"
    call set_solver_options(c_fieldA(field), ksptype='cg', pctype='eisenstat', rtol=1.0e-7, max_its=10000)
    c_fieldA(field)%option_path = "/c_fieldA" // int2str(field) // "/prognostic/conservative_interpolation_sobolev"
    call set_solver_options(c_fieldA(field), ksptype='cg', pctype='eisenstat', rtol=1.0e-7, max_its=10000)
    c_fieldA(field)%option_path = "/c_fieldA" // int2str(field)
    call allocate(l_fieldA(field), positionsA%mesh, "Field" // int2str(field))
  end do

  call set_from_function(c_fieldA(1), field_func_tophat, positionsA)
  call set_from_function(l_fieldA(1), field_func_tophat, positionsA)
  call set_from_function(c_fieldA(2), field_func_peaks, positionsA)
  call set_from_function(l_fieldA(2), field_func_peaks, positionsA)

  call vtk_write_fields("conservative_interpolation", 0, positionsA, positionsA%mesh, sfields=c_fieldA)
  call vtk_write_fields("linear_interpolation", 0, positionsA, positionsA%mesh, sfields=l_fieldA)

  do i=1,99
    write(0,'(a, i0)') "loop: ", i
    positionsB = read_mesh_files("data/input." // int2str(i+1), quad_degree=4, format="gmsh")
    do field=1,field_count
      call allocate(c_fieldB(field), positionsB%mesh, "Field" // int2str(field))
      c_fieldB(field)%option_path = "/c_fieldB" // int2str(field) // "/prognostic/conservative_interpolation_galerkin"
      call set_solver_options(c_fieldB(field), ksptype='cg', pctype='eisenstat', rtol=1.0e-10, max_its=20000)
      c_fieldB(field)%option_path = "/c_fieldB" // int2str(field) // "/prognostic/conservative_interpolation_bounded"
      call set_solver_options(c_fieldB(field), ksptype='cg', pctype='eisenstat', rtol=1.0e-10, max_its=20000)
      c_fieldB(field)%option_path = "/c_fieldB" // int2str(field) // "/prognostic/conservative_interpolation_sobolev"
      call set_solver_options(c_fieldB(field), ksptype='cg', pctype='eisenstat', rtol=1.0e-10, max_its=20000)
      c_fieldB(field)%option_path = "/c_fieldB" // int2str(field)
      call allocate(l_fieldB(field), positionsB%mesh, "Field" // int2str(field))
    end do

    call conservative_interpolation_bounded(c_fieldA, positionsA, c_fieldB, positionsB)
    call linear_interpolation(l_fieldA, positionsA, l_fieldB, positionsB)
    do field=1,field_count
      c_integral(field) = field_integral(c_fieldB(field), positionsB)
      l_integral(field) = field_integral(l_fieldB(field), positionsB)
      call deallocate(c_fieldA(field))
      c_fieldA(field) = c_fieldB(field)
      call deallocate(l_fieldA(field))
      l_fieldA(field) = l_fieldB(field)
    end do
    write(0,*) "integrals: ", i, c_integral, l_integral

    call deallocate(positionsA)
    positionsA = positionsB
    call vtk_write_fields("conservative_interpolation", i, positionsB, positionsB%mesh, sfields=c_fieldB)
    call vtk_write_fields("linear_interpolation", i, positionsB, positionsB%mesh, sfields=l_fieldB)

    do field=1,field_count
      c_maxval(field) = maxval(c_fieldB(field))
      c_minval(field) = minval(c_fieldB(field))
      l_maxval(field) = maxval(l_fieldB(field))
      l_minval(field) = minval(l_fieldB(field))
    end do
    write(0,*) "maxval: ", i, c_maxval, l_maxval
    write(0,*) "minval: ", i, c_minval, l_minval

    call allocate(analytical_fields(1), positionsB%mesh, "AnalyticalField1")
    call allocate(analytical_fields(2), positionsB%mesh, "AnalyticalField2")
    call set_from_function(analytical_fields(1), field_func_tophat, positionsB)
    call set_from_function(analytical_fields(2), field_func_peaks, positionsB)

    do field=1,field_count
      call addto(analytical_fields(field), c_fieldB(field), -1.0)
      c_l2err(field) = norm2(analytical_fields(field), positionsB)
    end do

    call set_from_function(analytical_fields(1), field_func_tophat, positionsB)
    call set_from_function(analytical_fields(2), field_func_peaks, positionsB)

    do field=1,field_count
      call addto(analytical_fields(field), l_fieldB(field), -1.0)
      l_l2err(field) = norm2(analytical_fields(field), positionsB)
    end do

    write(0,*) "l2err: ", i, c_l2err, l_l2err

    call deallocate(analytical_fields(1))
    call deallocate(analytical_fields(2))
  end do

end subroutine compute_compare_interpolation

function field_func_const(pos) result(f)
  real, dimension(:), intent(in) :: pos
  real :: f

  f = 1.0
end function field_func_const

function field_func_linear(pos) result(f)
  real, dimension(:), intent(in) :: pos
  real :: f

  f = pos(1) + 3
end function field_func_linear

function field_func_quadratic(pos) result(f)
  real, dimension(:), intent(in) :: pos
  real :: f

  f = 42 * pos(1)**2 + 2.0 * pos(2)**2 + 3.0
end function field_func_quadratic

function field_func_cubic(pos) result(f)
  real, dimension(:), intent(in) :: pos
  real :: f

  f = 500.0 * pos(2)**3 + 201.0 * pos(1)**2 + 94.0 * pos(2) + 3.0
end function field_func_cubic

function field_func_tophat(pos) result(f)
  use vector_tools
  real, dimension(:), intent(in) :: pos
  real :: f
  real, parameter :: PI=4.0*atan(1.0)

  real :: r
  r = norm2(pos - (/0.0, 0.0/))
  if (r < 0.7) then
    f = 1.0
  else
    f = 0.0
  end if
end function field_func_tophat

function field_func_peaks(pos) result(f)
  real, dimension(:), intent(in) :: pos
  real :: f
  real, parameter :: PI=4.0*atan(1.0)
  real :: x, y

  x = pos(1); y = pos(2)

  f = 3 * (1-x)**2 * exp(-x**2 - (y+1)**2)
  f = f - 10 * ((x/5) - x**3 - y**5 ) * exp(-x**2 -y**2)
  f = f - (1.0/3) * exp(-(x+1)**2 -y**2)
end function field_func_peaks
