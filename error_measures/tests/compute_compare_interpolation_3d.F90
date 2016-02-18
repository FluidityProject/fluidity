subroutine compute_compare_interpolation_3d

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
  real :: c_integral, l_integral
  integer :: field, field_count
  integer :: i

  positionsA = read_mesh_files("data/cube.1", quad_degree=4, format="gmsh")

  field_count = 2

  do field=1,field_count
    call allocate(c_fieldA(field), positionsA%mesh, "Field" // int2str(field))
    c_fieldA(field)%option_path = "/c_fieldA" // int2str(field) // "/prognostic/conservative_interpolation_galerkin"
    call set_solver_options(c_fieldA(field), ksptype='cg', pctype='eisenstat', rtol=1.0e-10, max_its=20000)
    c_fieldA(field)%option_path = "/c_fieldA" // int2str(field) // "/prognostic/conservative_interpolation_bounded"
    call set_solver_options(c_fieldA(field), ksptype='cg', pctype='eisenstat', rtol=1.0e-10, max_its=20000)
    c_fieldA(field)%option_path = "/c_fieldA" // int2str(field) // "/prognostic/conservative_interpolation_sobolev"
    call set_solver_options(c_fieldA(field), ksptype='cg', pctype='eisenstat', rtol=1.0e-10, max_its=20000)
    c_fieldA(field)%option_path = "/c_fieldA" // int2str(field)
    call allocate(l_fieldA(field), positionsA%mesh, "Field" // int2str(field))
  end do

  call set_from_function(c_fieldA(1), field_func_tophat, positionsA)
  call set_from_function(l_fieldA(1), field_func_tophat, positionsA)
  call set_from_function(c_fieldA(2), field_func_peaks, positionsA)
  call set_from_function(l_fieldA(2), field_func_peaks, positionsA)

  call vtk_write_fields("conservative_interpolation", 0, positionsA, positionsA%mesh, sfields=c_fieldA)
  call vtk_write_fields("linear_interpolation", 0, positionsA, positionsA%mesh, sfields=l_fieldA)
  do field=1,field_count
    c_integral = field_integral(c_fieldA(field), positionsA)
    l_integral = field_integral(l_fieldA(field), positionsA)
    write(0,*) "field: ", field, c_integral, l_integral
  end do

  do i=1,1
    write(0,'(a, i0)') "loop: ", i
    positionsB = read_mesh_files("data/cube." // int2str(i+1), quad_degree=4, format="gmsh")
    do field=1,field_count
      call allocate(c_fieldB(field), positionsB%mesh, "Field" // int2str(field))
      c_fieldB(field)%option_path = "/c_fieldB" // int2str(field) // "/prognostic/conservative_interpolation_galerkin"
      call set_solver_options(c_fieldB(field), ksptype='cg', pctype='eisenstat', rtol=1.0e-7, max_its=10000)
      c_fieldB(field)%option_path = "/c_fieldB" // int2str(field) // "/prognostic/conservative_interpolation_bounded"
      call set_solver_options(c_fieldB(field), ksptype='cg', pctype='eisenstat', rtol=1.0e-7, max_its=10000)
      c_fieldB(field)%option_path = "/c_fieldB" // int2str(field) // "/prognostic/conservative_interpolation_sobolev"
      call set_solver_options(c_fieldB(field), ksptype='cg', pctype='eisenstat', rtol=1.0e-7, max_its=10000)
      c_fieldB(field)%option_path = "/c_fieldB" // int2str(field)
      call allocate(l_fieldB(field), positionsB%mesh, "Field" // int2str(field))
    end do

    call cg_interpolation_galerkin_scalars_onfly(c_fieldA, positionsA, c_fieldB, positionsB)
    call linear_interpolation(l_fieldA, positionsA, l_fieldB, positionsB)
    do field=1,field_count
      c_integral = field_integral(c_fieldB(field), positionsB)
      l_integral = field_integral(l_fieldB(field), positionsB)
      write(0,*) "field: ", field, c_integral, l_integral
      call deallocate(c_fieldA(field))
      c_fieldA(field) = c_fieldB(field)
      call deallocate(l_fieldA(field))
      l_fieldA(field) = l_fieldB(field)
    end do

    call deallocate(positionsA)
    positionsA = positionsB
    call vtk_write_fields("conservative_interpolation", i, positionsB, positionsB%mesh, sfields=c_fieldB)
    call vtk_write_fields("linear_interpolation", i, positionsB, positionsB%mesh, sfields=l_fieldB)
  end do

end subroutine compute_compare_interpolation_3d

function field_func_const(pos) result(f)
  real, dimension(:), intent(in) :: pos
  real :: f

  f = 1.0
end function field_func_const

function field_func_linear(pos) result(f)
  real, dimension(:), intent(in) :: pos
  real :: f

  f = 85 * pos(1) + 23 * pos(2)
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
  r = norm2(pos(:2) - (/0.0, 0.0/))
  if (r < 0.7) then
    f = 1.0
  else
    f = 0.0
  end if
  f = 1.0
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
