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

module mba3d_integration 

  use adapt_integration
  use quadrature
  use elements
  use fields
  use global_parameters, only : real_8
  use halos
  use limit_metric_module
  use node_locking
#ifdef HAVE_MBA_3D
  use mba3d_mba_nodal
#endif
  use spud
  use surface_id_interleaving

  implicit none
  
  private

  public :: adapt_mesh_mba3d, mba3d_integration_check_options
  
  character(len = *), parameter :: base_path = "/mesh_adaptivity/hr_adaptivity"
  
contains

  subroutine adapt_mesh_mba3d(input_positions, metric, output_positions, force_preserve_regions)
    !!< Adapt the supplied input mesh using libmba3d. Return the new adapted
    !!< mesh in output_positions (which is allocated by this routine).
    !!< input_positions and output_positions are the Coordinate fields of
    !!< the old and new meshes respectively.
    
    type(vector_field), intent(in) :: input_positions
    type(tensor_field), intent(in) :: metric
    type(vector_field), intent(out) :: output_positions
    logical, intent(in), optional :: force_preserve_regions
    
    ! Linear tets only
    integer, parameter :: dim = 3, nloc = 4, snloc = 3

    integer :: i, max_coplanar_id
    integer, dimension(:), allocatable :: boundary_ids, coplanar_ids, sndgln
    real, parameter :: limit_buffer = 5.0, memory_buffer = 2.0 
    type(element_type) :: output_shape
    type(mesh_type) :: output_mesh
    type(quadrature_type) :: output_quad

    ! mbanodal arguments
    ! Group (M)
    integer :: np, maxp, nf, maxf, ne, maxe
    real(kind = real_8), dimension(:, :), allocatable :: xyp
    integer, dimension(:, :), allocatable :: ipf, ipe
    integer, dimension(:), allocatable :: lbf, lbe
    integer :: nestar
    ! Group (Dev)
    integer :: npv, nfv, nev
    integer, dimension(:), allocatable :: ipv, ifv, iev
    logical :: flagauto
    integer :: status
    ! Group (Q)
    integer :: maxskipe, maxqitr
    real(kind = real_8), dimension(:, :), allocatable :: metric_handle
    real(kind = real_8) :: quality, rquality
    ! Group (W)
    integer :: maxwr, maxwi
    real(kind = real_8), dimension(:), allocatable :: rw
    integer, dimension(:), allocatable :: iw
    integer :: iprint, ierr

    ewrite(1, *) "In adapt_mesh_mba_3d"
    
    assert(input_positions%dim == 3)
    assert(ele_loc(input_positions, 1) == 4)
#ifdef DDEBUG
    if(surface_element_count(input_positions) > 0) then
      assert(associated(input_positions%mesh%faces))
      assert(face_loc(input_positions, 1) == 3)
    end if
#endif
    assert(metric%mesh == input_positions%mesh)

    ewrite(2, *) "Forming mbanodal arguments"

    ewrite(2, *) "Forming group (M) arguments"

    nestar = expected_elements(input_positions, metric, global = .false.)

    ! Factor of limit_buffer buffers in limits
    np = node_count(input_positions)
    maxp = max_nodes(input_positions, expected_nodes(input_positions, nestar, global = .false.)) * limit_buffer
    nf = unique_surface_element_count(input_positions%mesh)
    maxf = max(int(((maxp * 1.0) / (np * 1.0)) * nf), nf) * limit_buffer + 1
    ne = ele_count(input_positions)
    maxe = max(nestar, ne) * limit_buffer + 1

    allocate(xyp(dim, maxp))
    xyp = 0.0
    do i = 1, dim
      xyp(i, :np) = input_positions%val(i,:)
    end do

    allocate(ipf(snloc, maxf))
    ipf = 0
    allocate(sndgln(nf * snloc))
    call getsndgln(input_positions%mesh, sndgln)
    ipf(:, :nf) = reshape(sndgln, (/snloc, nf/))
    deallocate(sndgln)

    allocate(ipe(nloc, maxe))
    ipe = 0
    ipe(:, :ne) = reshape(input_positions%mesh%ndglno, (/nloc, ne/))

    allocate(lbf(maxf))
    lbf = 0
    call interleave_surface_ids(input_positions%mesh, lbf(:nf), max_coplanar_id)
    if(minval(lbf(:nf)) < 0) then
      FLExit("libmba3d does not permit negative surface IDs")
    end if
    if(maxval(lbf(:nf)) >= huge(0)) then
      FLAbort("Exceeded max surface ID allowed by libmba3d")
    end if
    ! Offset surface IDs by 1, as libmba3d requires them to be positive
    lbf(:nf) = lbf(:nf) + 1
    if(maxval(lbf(:nf)) >= maxf - 100) then
      FLAbort("Exceeded max surface ID allowed by libmba3d")
    end if
    
    allocate(lbe(maxe))
    if(associated(input_positions%mesh%region_ids).and.&
       (have_option("/mesh_adaptivity/hr_adaptivity/preserve_mesh_regions")&
                              .or.present_and_true(force_preserve_regions))) then
      lbe(:ne) = input_positions%mesh%region_ids
      if(minval(lbe(:ne)) < 0) then
        FLExit("libmba3d does not permit negative region IDs")
      end if
      if(maxval(lbe(:ne)) >= huge(0)) then
        FLExit("Exceeded max region ID allowed by libmba3d")
      end if
      ! Offset region IDs by 1, as libmba3d requires them to be positive
      lbe(:ne) = lbe(:ne) + 1
    else
      lbe = 1
    end if

    ewrite(2, *) "Forming group (Dev) arguments"

    ! Locked nodes
    call get_locked_nodes(input_positions, ipv)
    npv = size(ipv)

    ! Locked faces
    nfv = 0
    allocate(ifv(nfv))

    ! Locked elements
    nev = 0
    allocate(iev(nev))

    flagauto = .true.
    status = 0 

    ewrite(2, *) "Forming group (Q) arguments"

    maxskipe = ne
    call get_option(base_path // "/adaptivity_library/libmba3d/max_optimisations", maxqitr, default = 100000)

    allocate(metric_handle(6, np))
    do i = 1, np
      metric_handle(1, i) = node_val(metric, 1, 1, i)
      metric_handle(2, i) = node_val(metric, 2, 2, i)
      metric_handle(3, i) = node_val(metric, 3, 3, i)
      metric_handle(4, i) = node_val(metric, 1, 2, i)
      metric_handle(5, i) = node_val(metric, 2, 3, i)
      metric_handle(6, i) = node_val(metric, 1, 3, i)
    end do

    call get_option(base_path // "/adaptivity_library/libmba3d/quality", quality, default = real(0.6, kind = real_8))
    ! Output variable - initialise just in case it's also used as input
    rquality = 0.0

    ewrite(2, *) "Forming group (W) arguments"

    maxwr = (14 * maxp + 2 * np + maxe) * memory_buffer + 1
    allocate(rw(maxwr))

    maxwi = (7 * maxp + np + 7 * maxf + 18 * maxe + 13 * ne) * memory_buffer + 1
    allocate(iw(maxwi))

    iprint = min(max(current_debug_level *  5, 0), 9)
    ! Output variable - initialise just in case it's also used as input
    ierr = 0

    ewrite(1, *) "Calling mbanodal from adapt_mesh_mba3d"
#ifdef HAVE_MBA_3D
    call mbanodal( &
      ! Group (M)
      & np, maxp, nf, maxf, ne, maxe, &
      & xyp, ipf, ipe, lbf, lbe, &
      & nestar, &
      ! Group (Dev)
      & npv, nfv, nev, ipv, ifv, iev, &
      & flagauto, status, &
      ! Group (Q)
      & maxskipe, maxqitr, &
      & metric_handle, quality, rquality, &
      ! Group (W)
      & maxwr, maxwi, rw, iw, &
      & iprint, ierr)
#else
    FLExit("Fluidity compiled without libmba3d support")
#endif
    ewrite(1, *) "Exited mbanodal"

    if(.not. any(ierr == (/0, 1000/))) then
      ewrite(-1, *) "libmba3d error code: ", ierr
      FLExit("libmba3d returned with an error")
    end if

    deallocate(ipv)
    deallocate(ifv)
    deallocate(iev)

    deallocate(metric_handle)

    deallocate(rw)
    deallocate(iw)

    ! Undo surface ID offset
    lbf(:nf) = lbf(:nf) - 1
    ! Undo region ID offset
    lbe(:ne) = lbe(:ne) - 1

    ewrite(2, *) "Target mesh quality: ", quality
    ewrite(2, *) "Mesh quality reported by mba3d: ", rquality
    if(ierr == 1000) then
      ewrite(1, *) "Warning: Target mesh quality not reached"
    end if

    if(isparallel()) then
      ewrite(2, *) "Constructing output halo"

      FLExit("libmba3d not available in parallel - halo recovery not implemented")

      ewrite(2, *) "Finished constructing output halo"
    end if

    ewrite(2, *) "Constructing output positions"

    ! Construct the new mesh
    output_quad = make_quadrature(nloc, dim, degree = input_positions%mesh%shape%quadrature%degree)
    output_shape = make_element_shape(nloc, dim, input_positions%mesh%shape%degree, output_quad)
    call allocate(output_mesh, np, ne, output_shape, name = input_positions%mesh%name)
    call deallocate(output_quad)
    call deallocate(output_shape)

    output_mesh%ndglno = reshape(ipe(:, :ne), (/ne * nloc/))
    allocate(boundary_ids(nf))
    allocate(coplanar_ids(nf))
    call deinterleave_surface_ids(lbf(:nf), max_coplanar_id, boundary_ids, coplanar_ids)
    call add_faces(output_mesh, sndgln = reshape(ipf(:, :nf), (/nf * snloc/)), boundary_ids = boundary_ids)
    deallocate(boundary_ids)
    if (nf/=surface_element_count(output_mesh)) then
      ! add_faces has duplicated internal boundary facets - this needs to be fixed
      ! see the mba2d wrapper
      FLAbort("Mba3d wrapper does not support internal boundary facets.")
    end if
    allocate(output_mesh%faces%coplanar_ids(nf))
    output_mesh%faces%coplanar_ids = coplanar_ids
    deallocate(coplanar_ids)

    if(have_option("/mesh_adaptivity/hr_adaptivity/preserve_mesh_regions")&
                              .or.present_and_true(force_preserve_regions)) then
      allocate(output_mesh%region_ids(ne))
      output_mesh%region_ids = lbe(:ne)
    end if
    output_mesh%option_path = input_positions%mesh%option_path
    output_mesh%periodic = input_positions%mesh%periodic

    ! Construct the new positions
    call allocate(output_positions, dim, output_mesh, name = input_positions%name)
    call deallocate(output_mesh)
               
    do i = 1, dim
      output_positions%val(i,:) = xyp(i, :np)
    end do
    output_positions%option_path = input_positions%option_path
    output_positions%multivalued_halo = input_positions%multivalued_halo

    deallocate(xyp)
    deallocate(ipf)
    deallocate(ipe)
    deallocate(lbf)
    deallocate(lbe)

    ewrite(2, *) "Finished constructing output positions"

    ewrite(1, *) "Exiting adapt_mesh_mba3d"

  end subroutine adapt_mesh_mba3d
  
  subroutine mba3d_integration_check_options
    character(len = *), parameter :: base_path = "/mesh_adaptivity/hr_adaptivity"
    integer :: dim, stat
    
    if(.not. have_option(base_path) .or. .not. have_option(base_path // "/adaptivity_library/libmba3d")) then
      ! Nothing to check
      return
    end if

#ifndef HAVE_MBA_3D
    FLExit("Cannot use libmba3d without the libmba3d library. Reconfigure with --enable-mba3d")
#endif
  
    call get_option("/geometry/dimension", dim, stat)
    if(stat /= SPUD_NO_ERROR) then
      ! This isn't the place to complain about this error
      return
    else if(dim /= 3) then
      FLExit("libmba3d can only be used in 3D")
    else if(isparallel()) then
      FLExit("libmba3d can only be used in serial")
    end if
  
  end subroutine mba3d_integration_check_options

end module mba3d_integration
