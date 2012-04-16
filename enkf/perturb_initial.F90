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
  !    amcgsoftware@imperial.ac.uk
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
program perturb_initial
  use spud
  use fields
  use state_module
  use write_state_module
  use timeloop_utilities
  use global_parameters, only: option_path_len, current_time, dt, FIELD_NAME_LEN
  use FLDebug
  use snapsvd_module
  use vtk_interfaces
  use memory_diagnostics
  use coordinates
  use Field_Options
  use populate_state_module
  use NetCDFWriter
  use m_random
  use read_netcdf_file
  use write_netcdf_file
  use futils, only: int2str
  implicit none

#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif

  character(len=FIELD_NAME_LEN) :: simulation_name, filename, filename_out
  real, dimension(:), allocatable :: x, y
  real, dimension(:,:), allocatable :: z
  real, dimension(:), allocatable :: z_tmp
  real, dimension(:,:), allocatable :: z_out

  integer :: nx, ny
  integer :: nstates, nnodes
!  type(scalar_field) :: field
!  type(state_type), dimension(:), pointer :: state
  real, dimension(:,:), allocatable :: E
  integer :: nrens, i, j, k
  real :: a

  call load_options('gyre_new.flml')
  call check_options
  call get_option("/simulation_name",simulation_name)

  !call initialise_write_state
  !call populate_state(state)

  call get_option("/material_phase/scalar_field::Pressure/prognostic/initial_condition::WholeMesh/from_netcdf/file_name", filename)
  !print*,filename

  !field=extract_scalar_field(state(1),'Pressure')
  !NNodes=node_count(field)
  !!input mesh nodes (2818)
  !nnodes=10

  !!Dimensions of X and Y are both nnodes.
  !allocate(X(nnodes), Y(nnodes), Z(nnodes))   
  !!read the netCDF(.nc) file and interpolate to the existing mesh
  !call get_field_values(trim(filename)//char(0), X, Y, Z, NNodes)

  nx=6
  ny=12
  nnodes=nx*ny 
  nrens=3

  allocate(x(nx), y(ny), z(ny,nx)) 
  allocate(z_tmp(nx*ny), z_out(ny,nx))
  allocate(E(nnodes, nrens))

  call simple_xy_rd(filename, nx, ny, z) 
  !print*,z(1,:)
  

  !!save z(ny, nx) to z_tmp(nx*ny) by columns
  do i=1, ny
     do j=1, nx
        z_tmp(i+(j-1)*ny)=z(i,j)
     enddo
  enddo  
  
  call randn(nnodes*nrens, E)
  !print*,E(:,1)

  do k=1, nrens
     E(:,k)=E(:,k)+z_tmp(:)
     !!save E(nx*ny,k) to z_out(nx,ny)
     do i=1, ny
        do j=1, nx
           z_out(i,j)=E(i+(j-1)*ny, k)
        enddo
     enddo 
     !print*,z_out(1,:)
     !!Save E(:,k) to nrens .nc files
     write(filename_out, '(a, i0, a)') 'ensemble_', k , ".nc"
     call simple_xy_wr(filename_out, nx, ny, z_out)
  enddo



  !!variable(:,:,:)
  !!call NetCDFWriter_write_variable(name, long_name, variable, units)



end
