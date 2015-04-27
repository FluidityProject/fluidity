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
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
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

subroutine test_strip_level_2_halo

  use fields
  use halos
  use parallel_tools
  use mesh_files
  use reserve_state_module
  use sam_integration
  use spud
  use state_module
  use unittest_tools
  
  implicit none

#ifdef HAVE_ADAPTIVITY
  interface
    subroutine flstriph2(nnodes, nprivatenodes, nprocs, &
      & volumeenlist, nvolumeelems, nloc, &
      & surfaceenlist, surfaceids, nsurfaceelems, snloc, &
      & x, y, z, &
      & fields, nfields, fstride, &
      & metric, &
      & scatter, nscatter)
      implicit none
      integer, intent(inout) :: nnodes
      integer, intent(in) :: nprivatenodes
      integer, intent(in) :: nprocs
      integer, intent(inout) :: nvolumeelems
      integer, intent(in) :: nloc
      integer, dimension(nvolumeelems * nloc), intent(inout) :: volumeenlist
      integer, intent(inout) :: nsurfaceelems
      integer, intent(in) :: snloc
      integer, dimension(nsurfaceelems * snloc), intent(inout) :: surfaceenlist
      integer, dimension(nsurfaceelems), intent(inout) :: surfaceids
      real, dimension(nnodes), intent(inout) :: x
      real, dimension(nnodes), intent(inout) :: y
      real, dimension(nnodes), intent(inout) :: z
      integer, intent(inout) :: nfields
      integer, intent(inout) :: fstride
      real, dimension(nnodes * nfields * fstride), intent(inout) :: fields
      real, dimension(nnodes * 9), intent(inout) :: metric
      integer, intent(inout) :: nscatter
      integer, dimension(nscatter), intent(inout) :: scatter
    end subroutine flstriph2
  end interface

  integer :: dim, i, j, new_nsurfaceelems, stat
  integer, dimension(:), allocatable :: new_scatter, new_surfaceenlist, new_surfaceids, nreceives
  type(halo_type), pointer :: halo
  type(mesh_type), pointer :: mesh
  type(state_type) :: state, state_array(1)
  type(vector_field), target :: mesh_field

  ! flstriph2 variables
  integer :: nnodes, nprivatenodes, nprocs
  integer, dimension(:), allocatable :: volumeenlist
  integer :: nvolumeelems, nloc
  integer, dimension(:), allocatable :: surfaceenlist, surfaceids
  integer :: nsurfaceelems, snloc
  real, dimension(:), allocatable :: x, y, z, input_fields
  integer :: nfields, fstride
  real, dimension(:), allocatable :: metric
  integer, dimension(:), allocatable :: scatter
  integer :: nscatter

  mesh_field = read_mesh_files("data/cube-parallel_0", quad_degree = 1, format="gmsh")
  call read_halos("data/cube-parallel", mesh_field)
  assert(halo_count(mesh_field) > 0)
  halo => mesh_field%mesh%halos(1)

  mesh => mesh_field%mesh

  assert(mesh_dim(mesh_field) == 3)
  dim = mesh_dim(mesh_field)

  nnodes = node_count(mesh)
  nprivatenodes = halo_nowned_nodes(halo)
  nprocs = halo_proc_count(halo)

  nvolumeelems = ele_count(mesh)
  assert(nvolumeelems > 0)
  nloc = ele_loc(mesh, 1)
#ifdef DDEBUG
  do i = 2, nvolumeelems
    assert(ele_loc(mesh, i) == nloc)
  end do
#endif
  allocate(volumeenlist(nvolumeelems * nloc))
  volumeenlist = mesh%ndglno

  nsurfaceelems = surface_element_count(mesh)
  assert(nsurfaceelems > 0)
  snloc = face_loc(mesh, 1)
#ifdef DDEBUG
  do i = 2, nsurfaceelems
    assert(face_loc(mesh, i) == snloc)
  end do
#endif
  allocate(surfaceenlist(nsurfaceelems * snloc))
  call getsndgln(mesh, surfaceenlist)
  allocate(surfaceids(nsurfaceelems))
  if(associated(mesh%faces%boundary_ids)) then
    surfaceids = mesh%faces%boundary_ids
  else
    surfaceids = 0
  end if

  allocate(x(nnodes))
  allocate(y(nnodes))
  allocate(z(nnodes))
  x = mesh_field%val(1,:)
  y = mesh_field%val(2,:)
  z = mesh_field%val(3,:)

  nfields = 0
  fstride = nnodes
  allocate(input_fields(nfields * fstride))
  if(size(input_fields) > 0) then
    input_fields = 0
  end if

  allocate(metric(nnodes * dim ** 2))
  ! Construct a unit matrix
  metric = 0.0
  do i = 1, nnodes
    do j = 1, dim
      metric((i - 1) * dim ** 2 + (j - 1) * dim + j) = 1.0
    end do
  end do

  nscatter = halo_all_receives_count(halo)
  allocate(scatter(nscatter))
  allocate(nreceives(halo_proc_count(halo)))
  call extract_all_halo_receives(halo, scatter, nreceives)
  deallocate(nreceives)

  call flstriph2(nnodes, nprivatenodes, nprocs, &
    & volumeenlist, nvolumeelems, nloc, &
    & surfaceenlist, surfaceids, nsurfaceelems, snloc, &
    & x, y, z, &
    & input_fields, nfields, fstride, &
    & metric, &
    & scatter, nscatter)

  call report_test("[nprivatenodes unchanged]", nprivatenodes /= halo_nowned_nodes(halo), .false., "flstriph2 has changed nprivatenodes")
  call report_test("[nprocs unchanged]", nprocs /= halo_proc_count(halo), .false., "flstriph2 has changed nprocs")
  call report_test("[nloc unchanged]", nloc /= ele_loc(mesh, 1), .false., "flstriph2 has changed nloc")
  call report_test("[snloc unchanged]", snloc /= face_loc(mesh, 1), .false., "flstriph2 has changed snloc")
  call report_test("[nscatter unchanged]", nscatter /= halo_all_receives_count(halo), .false., "flstriph2 has changed nscatter")

  call insert(state, mesh, "CoordinateMesh")
  call insert(state, mesh_field, "Coordinate")

  call add_option("/geometry/mesh", stat = stat)
  call set_option_attribute("/geometry/mesh/name", "CoordinateMesh", stat = stat)
  assert(stat == SPUD_NEW_KEY_WARNING)
  call add_option("/geometry/mesh/from_file", stat = stat)
  assert(stat == SPUD_NEW_KEY_WARNING)

  call deallocate(mesh_field)

  state_array(1) = state
  call create_reserve_state(state_array)

  call strip_level_2_halo(state_array)
  state = state_array(1)

  mesh => extract_mesh(state, "CoordinateMesh")
  call report_test("[Same node count]", nnodes /= node_count(mesh), .false., "Node counts returned by strip_level_2_halo and flstriph2 are different")
  call report_test("[Same element count]", nvolumeelems /= ele_count(mesh), .false., "Element counts returned by strip_level_2_halo and flstriph2 are different")
  call report_test("[Same element numbering]", any(mesh%ndglno /= volumeenlist(1:nvolumeelems * nloc)), .false., "Element numberings returned by strip_level_2_halo and flstriph2 are different")

  allocate(new_surfaceenlist(surface_element_count(mesh) * face_loc(mesh, 1)))
  allocate(new_surfaceids(surface_element_count(mesh)))
  call getsndgln(mesh, new_surfaceenlist)
  assert(associated(mesh%faces%boundary_ids))
  new_surfaceids = mesh%faces%boundary_ids
  ! Strip off the surface elements with ID 0 added by add_faces
  new_nsurfaceelems = 0
  do i = 1, surface_element_count(mesh) 
    if(new_surfaceids(i) /= 0) then
      new_nsurfaceelems = new_nsurfaceelems + 1
      new_surfaceenlist((new_nsurfaceelems - 1) * face_loc(mesh, 1) + 1:new_nsurfaceelems * face_loc(mesh, 1)) = new_surfaceenlist((i - 1) * face_loc(mesh, 1) + 1:i * face_loc(mesh, 1))
      new_surfaceids(new_nsurfaceelems) = new_surfaceids(i)
    end if
  end do
  call report_test("[Same surface element count]", nsurfaceelems /= new_nsurfaceelems, .false., "Surface element counts returned by strip_level_2_halo and flstriph2 are different")
  call report_test("[Same surface element numbering]", any(new_surfaceenlist(1:new_nsurfaceelems * face_loc(mesh, 1)) /= surfaceenlist(1:nsurfaceelems * snloc)), .false., "Surface element numberings returned by strip_level_2_halo and flstriph2 are different")
  call report_test("[Same surface IDs]", any(new_surfaceids(1:new_nsurfaceelems) /= surfaceids(1:nsurfaceelems)), .false., "Surface IDs returned by strip_level_2_halo and flstriph2 are different")
  deallocate(new_surfaceenlist)
  deallocate(new_surfaceids)

  mesh_field = extract_vector_field(state, "Coordinate")
  assert(mesh_dim(mesh_field) == 3)
  call report_test("[Same x coordinates]", any(mesh_field%val(1,:) /= x(1:nnodes)), .false., "x coordinates returned by strip_level_2_halo and flstriph2 are different")
  call report_test("[Same y coordinates]", any(mesh_field%val(2,:) /= y(1:nnodes)), .false., "y coordinates returned by strip_level_2_halo and flstriph2 are different")
  call report_test("[Same z coordinates]", any(mesh_field%val(3,:) /= z(1:nnodes)), .false., "z coordinates returned by strip_level_2_halo and flstriph2 are different")

  call report_test("[Same number of level 1 receive nodes]", halo_all_receives_count(mesh_field%mesh%halos(1)) /= nscatter, .false., "Number of level 1 receive nodes returned by strip_level_2_halo and flstriph2 are different")
  allocate(new_scatter(halo_all_receives_count(mesh_field%mesh%halos(1))))
  allocate(nreceives(halo_proc_count(mesh_field%mesh%halos(1))))
  call extract_all_halo_receives(mesh_field%mesh%halos(1), new_scatter, nreceives)
  deallocate(nreceives)
  call report_test("[Same level 1 receive nodes]", any(new_scatter /= scatter(1:nscatter)), .false., "Level 1 receive nodes returned by strip_level_2_halo and flstriph2 are different")
  deallocate(new_scatter)

  deallocate(volumeenlist)
  deallocate(surfaceenlist)
  deallocate(surfaceids)
  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(scatter)
  deallocate(input_fields)
  deallocate(metric)

  call deallocate(state)

  call report_test_no_references()
#else
  call report_test("[test disabled]", .false., .true., "Test compiled without sam support")
#endif

end subroutine test_strip_level_2_halo
