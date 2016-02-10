#include "fdebug.h"

module mba2d_integration
  use global_parameters, only : current_debug_level
  use quadrature
  use elements
  use eventcounter
  use quicksort
  use data_structures
  use metric_tools
  use fields
  use state_module
  use meshdiagnostics
  use vtk_interfaces
  use halos
  use node_boundary
  use interpolation_module
  use limit_metric_module
#ifdef HAVE_MBA_2D
  use mba2d_module
#endif
  use mba_adapt_module
  use node_locking
  use surface_id_interleaving
  use adapt_integration
  implicit none
  
  private
  
  public :: adapt_mesh_mba2d, mba2d_integration_check_options

  contains

  subroutine adapt_mesh_mba2d(input_positions, metric, output_positions, force_preserve_regions, &
    lock_faces, allow_boundary_elements)
    type(vector_field), intent(in), target :: input_positions
    type(tensor_field), intent(in) :: metric
    type(vector_field), intent(out) :: output_positions
    logical, intent(in), optional :: force_preserve_regions
    type(integer_set), intent(in), optional :: lock_faces
    ! if present and true allow boundary elements, i.e. elements with
    ! all nodes on the boundary, if not present, the default
    ! in serial is to forbid boundary elements, and allow them in parallel
    logical, intent(in), optional :: allow_boundary_elements

#ifdef HAVE_MBA_2D

    type(mesh_type), pointer :: xmesh

    integer :: nonods, mxnods, orig_stotel, stotel, mxface, totele, maxele, stotel_external
    real, dimension(:, :), allocatable :: pos
    integer, dimension(:, :), allocatable :: ipf
    integer, dimension(:, :), allocatable :: ipe
    real, dimension(:, :), allocatable :: parcrv
    integer, dimension(:), allocatable :: ipv
    integer, dimension(:), allocatable :: ifv
    integer, dimension(:), allocatable :: iFnc
    integer, dimension(:), allocatable :: lbE
    real, dimension(:, :), allocatable :: tmp_metric
    integer :: i, j, partition_surface_id, face, face2
    real :: quality, rQuality
    integer :: iPrint, ierr, maxWr, maxWi
    real, dimension(:), allocatable :: rW
    integer, dimension(:), allocatable :: iW
    integer :: status
    type(mesh_type) :: new_mesh
    integer :: npv
    integer :: xpctel
    integer, dimension(:,:), allocatable:: new_sndgln
    integer :: maxp
    integer :: iterations
    integer, dimension(:), allocatable :: locked_nodes
    integer, dimension(:), pointer :: neighbours, faces
    type(halo_type), pointer :: old_halo, new_halo
    integer :: proc
    
    ! Surface ID interleaving
    integer :: max_coplanar_id
    integer, dimension(:), allocatable :: boundary_ids, coplanar_ids, surface_ids, mba_boundary_ids
    type(integer_hash_table) :: physical_surface_ids
    ! Element locking
    integer :: nfe
    integer, dimension(:), allocatable :: ife
    integer, dimension(:), pointer:: nodes
    integer :: ele
    integer :: nfv
    type(csr_sparsity), pointer :: nelist
    ! Flattened halo data
    integer :: nhalos
    type(integer_hash_table) :: input_face_numbering_to_mba2d_numbering

    ! if we're parallel we'll need to reorder the region ids after the halo derivation
    integer, dimension(:), allocatable :: old_new_region_ids, renumber_permutation
    
!#define DUMP_HALO_INTERPOLATION
#ifdef DUMP_HALO_INTERPOLATION
    type(mesh_type):: p0mesh
    type(scalar_field):: locked_field
    integer, save:: ix=0
#endif

#ifdef GIVE_LIPNIKOV_OUTPUT
    integer :: rank
    character(len=255) :: filename
#endif

!#define DUMP_HALO
#ifdef DUMP_HALO
    type(scalar_field) :: sends_sfield, receives_sfield
#endif

    ewrite(1, *) "In adapt_mesh_mba2d"

    assert(all(metric%dim == 2))

    xmesh => input_positions%mesh
    call deallocate_boundcount
    call initialise_boundcount(xmesh, input_positions)

    ! mxnods is an option to adaptivity specifying the maximum number of nodes
    xpctel = max(expected_elements(input_positions, metric), 5)
    mxnods = max_nodes(input_positions, expected_nodes(input_positions, xpctel, global = .false.))

    nonods = node_count(xmesh)
    totele = ele_count(xmesh)
    orig_stotel = unique_surface_element_count(xmesh)
    mxface = int(max((float(mxnods) / float(nonods)) * orig_stotel * 3.5, 10000.0))
    maxele = int(max((float(mxnods) / float(nonods)) * totele * 1.5, 10000.0))
    maxp = mxnods * 1.2

    allocate(pos(2, maxp))
    pos = 0.0
    do i=1,2
      pos(i, 1:nonods) = input_positions%val(i,:)
    end do

    allocate(surface_ids(orig_stotel))
    call interleave_surface_ids(xmesh, surface_ids, max_coplanar_id)

    allocate(ipf(4, mxface))
    ipf = 0
    partition_surface_id = maxval(surface_ids) + 1
    stotel = 0
    stotel_external = 0

    call allocate(input_face_numbering_to_mba2d_numbering)
    do i=1,totele
      neighbours => ele_neigh(xmesh, i)
      faces => ele_faces(xmesh, i)
      do j=1,3 ! 3 faces per element in 2D
        if (neighbours(j) <= 0) then
          face = faces(j)
          stotel = stotel + 1
          call insert(input_face_numbering_to_mba2d_numbering, face, stotel)
          ipf(1:2, stotel) = face_global_nodes(xmesh, face)
          ipf(3, stotel) = 0
          if (face <= orig_stotel) then ! if facet is genuinely external, i.e. on the domain exterior
            ipf(4, stotel) = surface_ids(face)
            stotel_external = stotel_external + 1
          else if (face <= surface_element_count(xmesh)) then
            ! this should only happen for the 2nd copy of an internal facet
            ! so something's wrong in the faces admin
            FLAbort("Detected external facet that's numbered incorrectly.")
          else
            ipf(4, stotel) = partition_surface_id
          end if
        end if
      end do
    end do

    if (stotel_external<orig_stotel) then
      ! not all facets in the surface mesh are external apparently
      ! 1:orig_stotel should only give us one of each pair of internal facets
      do face=1, orig_stotel
        face2 = face_neigh(xmesh, face)
        if (face/=face2) then
          ! if face==face2 it's an external facet that's dealt with already
          stotel = stotel +1
          call insert(input_face_numbering_to_mba2d_numbering, face, stotel)
          ipf(1:2, stotel) = face_global_nodes(xmesh, face)
          ipf(3, stotel) = 0
          ipf(4, stotel) = surface_ids(face)
          if (surface_element_id(xmesh, face)/=surface_element_id(xmesh, face2)) then
            ! check that the surface id on either side is indeed the same
            FLExit("Adaptivity with internal boundaries only works if the surface id is single valued")
          end if
        end if
      end do

    end if

    if (.not. present(lock_faces)) then
      nfv = 0
      allocate(ifv(nfv))
    else
      nfv = key_count(lock_faces)
      allocate(ifv(nfv))
      ifv = fetch(input_face_numbering_to_mba2d_numbering, set2vector(lock_faces))
    end if

    call deallocate(input_face_numbering_to_mba2d_numbering)
    
    deallocate(surface_ids)

    ! If you fail this, you need to know the following.
    ! Divide the surface elements into three classes:
    ! (a) physical surface elements (on the exterior of the domain)
    ! (b) partition surface elements (on the boundary with another partition)
    ! (c) internal surface elements (between elements of the mesh
    ! When we allocated mxface, we only had the count of (a).
    ! However, we actually want to pass (a) + (b) to mba2d.
    ! Mxface was allocated with a nice big multiple to make sure that
    ! there was enough space.
    ! Afterwards, we went through and counted (b) and added it on to stotel.
    ! But in a very odd situation, this might not be enough memory!
    ! So either find a smarter number to set mxface, or make the multiple bigger.
    if (stotel > mxface) then
      FLAbort("Expected number of facets too small!")
    end if

    allocate(ipe(3, maxele))
    ipe = 0
    do i=1,totele
      ipe(:, i) = ele_nodes(xmesh, i)
    end do

    nhalos = halo_count(xmesh)
    assert(any(nhalos == (/0, 1, 2/)))
    if(nhalos > 0) then
#ifdef DUMP_HALO_INTERPOLATION
      p0mesh = piecewise_constant_mesh(xmesh, "P0Mesh")
      call allocate(locked_field, p0mesh, "Locked")
#endif
      nelist => extract_nelist(xmesh)
    
      old_halo => xmesh%halos(nhalos)

      allocate(ife(totele))
      nfe = 0
      ele_loop: do ele=1, element_count(xmesh)
        nodes => ele_nodes(xmesh, ele)
        do j=1, size(nodes)
          if (.not. node_owned(old_halo, nodes(j))) then
            nfe = nfe +1
            ife(nfe) = ele
#ifdef DUMP_HALO_INTERPOLATION
            call set(locked_field, ele, 1.0)
#endif
            cycle ele_loop
          end if
        end do
      end do ele_loop
#ifdef DUMP_HALO_INTERPOLATION
      call vtk_write_fields("locked", index=ix, position=input_positions, model=xmesh, sfields=(/ locked_field /))
      ix = ix+1
      call deallocate(locked_field)
      call deallocate(p0mesh)
#endif
    else
      nfe = 0
      allocate(ife(nfe))
    end if
    
    ! construct list of nodes to be locked
    call get_locked_nodes(input_positions, locked_nodes)
    npv = count(boundcount > 1) + size(locked_nodes)
    if (nhalos>0) then
      npv = npv + halo_all_sends_count(old_halo) + halo_all_receives_count(old_halo)
    end if
    allocate(ipv(npv))
    j = 1

    if (nhalos>0) then
      ! lock the send and receive nodes - these are locked above already as halo 1 elements
      ! but locking them here again, get back a list in the new numbering after the adapt
      ! thus allowing us to reconstruct the halo
      do proc=1, halo_proc_count(old_halo)
        ipv(j:j+halo_send_count(old_halo, proc)-1) = halo_sends(old_halo, proc)
        j = j + halo_send_count(old_halo, proc)
      end do
      do proc=1, halo_proc_count(old_halo)
        ipv(j:j+halo_receive_count(old_halo, proc)-1) = halo_receives(old_halo, proc)
        j = j + halo_receive_count(old_halo, proc)
      end do
    end if

    ! lock nodes on the boundary that are adjacent to more than one coplanar id, i.e. corner nodes
    ! (see node_boundary module)
    do i=1,nonods
      if (boundcount(i) > 1) then
        ipv(j) = i
        j = j + 1
      end if
    end do
    ! nodes locked as prescribed by python
    ipv(j:) = locked_nodes
    deallocate(locked_nodes)

    allocate(parcrv(2, mxface))
    parcrv = 0


    allocate(iFnc(mxface))
    iFnc = 0

    allocate(lbE(maxele))
    lbE = 1
    if ((associated(xmesh%region_ids)).and.&
        (have_option("/mesh_adaptivity/hr_adaptivity/preserve_mesh_regions")&
         .or.present_and_true(force_preserve_regions))) then
      ! offset surface IDs by 1 because libmba2d requires them to be positive
      lbE(1:totele) = xmesh%region_ids + 1
    end if

    allocate(tmp_metric(3, maxp))
    tmp_metric = 0.0
    do i=1,nonods
      tmp_metric(1, i) = node_val(metric, 1, 1, i)
      tmp_metric(2, i) = node_val(metric, 2, 2, i)
      tmp_metric(3, i) = node_val(metric, 1, 2, i)
    end do

    call relax_metric_locked_regions(tmp_metric, ife(1:nfe), input_positions)

    maxWr = (4 * maxp + 10 * nonods + mxface + maxele)  * 1.5
    maxWi = (6 * maxp + 10 * nonods + 19 * mxface + 11 * maxele + 12 * totele) * 1.5
    allocate(rW(maxWr))
    allocate(iW(maxWi))

    if (present(allow_boundary_elements)) then
      if (allow_boundary_elements) then
        status=0
      else
        status=1
      end if
    else if(nhalos > 0) then
      ! we don't want to avoid boundary elements along the local domain
      ! boundaries, as the ragged boundary will usually have a lot of
      ! triangles with 2 faces on the boundary and we don't want to split
      ! these up unnecessarily - unfortunately we can't only allow it 
      ! along local domain boundaries and forbid them on the global domain boundary
      status = 0 ! allow boundary elements
    else
      status = 1 ! forbid boundary elements
    end if

    ! Now we decide how many iterations the library
    ! should do.
    ! Say the number of elements (totel) >>
    !     the desired number of elements (xpctel).
    ! Then it takes at least one iteration to remove it.
    ! So, we should calibrate the number of iterations
    ! we let it do, depending on how much work we expect it
    ! to take.
    iterations = max(50000, int(abs(totele - xpctel)*1.2))

    ! Fine-tuning options
    call get_option("/mesh_adaptivity/hr_adaptivity/adaptivity_library/libmba2d/quality", quality, default = 0.6)

#ifdef GIVE_LIPNIKOV_OUTPUT
    rank = getrank()
    write(filename, '(a,i0,a)') "debug_", rank, ".ani"
    call saveMani(nonods, stotel, totele, npv, 0, nfe, &
                  pos, ipf, ipe, ipv, ipv, ife, lbE, &
                  parcrv, iFnc, filename)
#endif

    iprint = min(max(current_debug_level *  5, 0), 9)

    call mbaNodal(                                   &
         nonods, maxp, stotel, mxface, totele, maxele, npv, &
         pos, ipf, ipe, ipv, &
         CrvFunction_ani, parcrv, iFnc, &
         xpctel, &
         nfv, nfe, ifv, ife, lbE, &
         .true., status, &
         100, iterations, &
         tmp_metric, quality, rQuality, &
         maxWr, maxWi, rW, iW, &
         iPrint, ierr)

    call incrementeventcounter(EVENT_ADAPTIVITY)
    call incrementeventcounter(EVENT_MESH_MOVEMENT)

    ! Hooray! You didn't crash. Congratulations. Now let's assemble the output and interpolate.

    call allocate(new_mesh, nonods, totele, ele_shape(xmesh, 1), trim(xmesh%name))
    ! Hack: untag these references so that people (i.e. me) don't get confused.
    new_mesh%shape%refcount%tagged = .false.
    new_mesh%shape%quadrature%refcount%tagged = .false.
    new_mesh%ndglno = reshape(IPE(:, 1:totele), (/size(new_mesh%ndglno)/))
    new_mesh%option_path = xmesh%option_path
    
    if (.not. isparallel()) then
      allocate(new_sndgln(1:2,1:stotel), mba_boundary_ids(1:stotel))
      new_sndgln = ipf(1:2,1:stotel)
      mba_boundary_ids =ipf(4,1:stotel)
    else
      ! In parallel, we need to filter out the surface elements with colour partition_surface_id, because
      ! they are not real external faces
      call allocate(physical_surface_ids)

      j = 1
      do i=1,stotel
        if (ipf(4, i) /= partition_surface_id) then
          call insert(physical_surface_ids, j, i)
          j = j + 1
        end if
      end do

      ! number of surface elements without inter-partition surface elements
      stotel = j - 1
      allocate(mba_boundary_ids(1:stotel), new_sndgln(1:2,1:stotel))
      do i=1, stotel
        mba_boundary_ids(i) = ipf(4, fetch(physical_surface_ids, i))
        new_sndgln(1:2, i) = ipf(1:2, fetch(physical_surface_ids, i))
      end do
      
      do i=1, stotel
        new_sndgln(1:2, i) = ipf(1:2, fetch(physical_surface_ids, i))
      end do
      call deallocate(physical_surface_ids)
    end if

    ! add_faces might create extra (internal) surface elements, so we 
    ! use the combined boundary+coplanar ids first
    call add_faces(new_mesh, sndgln=reshape(new_sndgln, (/ 2*stotel /) ), &
      boundary_ids=mba_boundary_ids)
    deallocate(mba_boundary_ids)
    deallocate(new_sndgln)

    ! and only deinterleave now we know the total number of elements in the surface mesh
    ! add_faces will have copied the interleaved id to the second copy of each interior facet
    stotel = surface_element_count(new_mesh)
    allocate(boundary_ids(1:stotel), coplanar_ids(1:stotel))
    call deinterleave_surface_ids(new_mesh%faces%boundary_ids, max_coplanar_id, boundary_ids, coplanar_ids)

    new_mesh%faces%boundary_ids = boundary_ids

    if(associated(xmesh%faces%coplanar_ids)) then
      allocate(new_mesh%faces%coplanar_ids(1:stotel))
      new_mesh%faces%coplanar_ids = coplanar_ids
    end if
    deallocate(boundary_ids, coplanar_ids)
    
    if(have_option("/mesh_adaptivity/hr_adaptivity/preserve_mesh_regions")&
                              .or.present_and_true(force_preserve_regions)) then
      allocate(new_mesh%region_ids(totele))
      new_mesh%region_ids = lbE(1:totele) - 1
    end if


    call allocate(output_positions, 2, new_mesh, trim(input_positions%name))
    output_positions%option_path = input_positions%option_path
    call deallocate(new_mesh)
    call set_all(output_positions, pos(:, 1:nonods))
    
    if(nhalos > 0) then

      allocate(output_positions%mesh%halos(nhalos))
      new_halo => output_positions%mesh%halos(nhalos)
      
      ! halo is the same in terms of n/o sends and receives, name, ordering_type, etc.
      call allocate(new_halo, old_halo)
      ! except for n/o owned nodes
      call set_halo_nowned_nodes(new_halo, nonods-halo_all_receives_count(new_halo))

      j = 1
      do proc=1, halo_proc_count(new_halo)
        call set_halo_sends(new_halo, proc, ipv(j:j+halo_send_count(new_halo, proc)-1))
        j = j + halo_send_count(new_halo, proc)
      end do
      do proc=1, halo_proc_count(new_halo)
        call set_halo_receives(new_halo, proc, ipv(j:j+halo_receive_count(new_halo, proc)-1))
        j = j + halo_receive_count(new_halo, proc)
      end do
      
      allocate(renumber_permutation(totele))
      
      if(nhalos == 2) then
        ! Derive remaining halos
        call derive_l1_from_l2_halo(output_positions%mesh, ordering_scheme = HALO_ORDER_GENERAL, create_caches = .false.)
        ! Reorder the nodes for trailing receives consistency
        call renumber_positions_trailing_receives(output_positions)

        allocate(output_positions%mesh%element_halos(2))
        ! Reorder the elements for trailing receives consistency
        call derive_element_halo_from_node_halo(output_positions%mesh, &
          & ordering_scheme = HALO_ORDER_GENERAL, create_caches = .false.)
        call renumber_positions_elements_trailing_receives(output_positions, permutation=renumber_permutation)
      else
        ! Reorder the nodes for trailing receives consistency
        call renumber_positions_trailing_receives(output_positions)
        
        allocate(output_positions%mesh%element_halos(1))
        ! Reorder the elements for trailing receives consistency
        call derive_element_halo_from_node_halo(output_positions%mesh, &
          & ordering_scheme = HALO_ORDER_GENERAL, create_caches = .false.)
        call renumber_positions_elements_trailing_receives(output_positions, permutation=renumber_permutation)
      end if
      
      if(have_option("/mesh_adaptivity/hr_adaptivity/preserve_mesh_regions")&
                                .or.present_and_true(force_preserve_regions)) then
        ! reorder the region_ids since all out elements have been jiggled about
        allocate(old_new_region_ids(totele))
        old_new_region_ids = output_positions%mesh%region_ids
        do i = 1, totele
          output_positions%mesh%region_ids(renumber_permutation(i)) = old_new_region_ids(i)
        end do
        deallocate(old_new_region_ids)
      end if
      
      deallocate(renumber_permutation)
      

#ifdef DUMP_HALO
      call allocate(sends_sfield, output_positions%mesh, "Sends")
      call zero(sends_sfield)
      call allocate(receives_sfield, output_positions%mesh, "Receives")
      call zero(receives_sfield)

      do proc=1,halo_proc_count(output_positions%mesh%halos(1))
        do i=1,size(output_positions%mesh%halos(1)%sends(proc)%ptr)
          call set(sends_sfield, output_positions%mesh%halos(1)%sends(proc)%ptr(i), 1.0)
        end do
        do i=1,size(output_positions%mesh%halos(1)%receives(proc)%ptr)
          call set(receives_sfield, output_positions%mesh%halos(1)%receives(proc)%ptr(i), 1.0)
        end do
      end do
      
      call vtk_write_fields("halo", position=output_positions, model=output_positions%mesh, sfields=(/sends_sfield, receives_sfield/))

      call deallocate(sends_sfield)
      call deallocate(receives_sfield)
#endif

      ! Adaptivity is not guaranteed to return halo elements in the same
      ! order in which they went in. We therefore need to fix this order.
      call reorder_element_numbering(output_positions)

#ifdef DDEBUG
      do i = 1, nhalos
        assert(trailing_receives_consistent(output_positions%mesh%halos(i)))
        assert(halo_valid_for_communication(output_positions%mesh%halos(i)))
        assert(halo_verifies(output_positions%mesh%halos(i), output_positions))
      end do
#endif
    end if

    deallocate(pos)
    deallocate(ipf)
    deallocate(ipe)
    deallocate(ipv)
    deallocate(parcrv)
    deallocate(iFnc)
    deallocate(lbE)
    deallocate(rW)
    deallocate(iW)
    deallocate(tmp_metric)
    
    ewrite(1, *) "Exiting adapt_mesh_mba2d"

#else
    FLExit("You called mba_adapt without the mba2d library. Reconfigure with --enable-2d-adaptivity")
#endif
  end subroutine adapt_mesh_mba2d

  subroutine relax_metric_locked_regions(metric, locked_elements, positions)
    ! in the locked regions (halo regions in parallel) and the region immediately
    ! adjacent to it, mba can't satisfy what we ask for in the metric
    ! This tends to upset mba and sometimes leads it to give up altogether (leaving other regions 
    ! unadapted). Therefore we adjust the metric in the nodes of the locked regions to adhere to
    ! the locked elements (i.e. tell mba that the mesh is perfect there already). Directly adjacent to
    ! the locked region it will interpolate linearly between the overwritten metric and the metric we want,
    ! so that it still adapts to our desired quality everywhere it is allowed to.
    real, dimension(:,:), intent(inout):: metric
    integer, dimension(:), intent(in):: locked_elements
    type(vector_field), intent(in):: positions

    real, dimension(2,2) :: ele_metric
    real, dimension(3) :: ele_metric_vector
    integer, dimension(:), allocatable:: adjacent_locked_element_count
    integer, dimension(:), pointer:: nodes
    integer :: i, j, ele, n

    ! keep track of how many element metrics we've added into each node already
    allocate(adjacent_locked_element_count(1:node_count(positions)))
    adjacent_locked_element_count = 0

    do i=1, size(locked_elements)
      ele = locked_elements(i)
      ele_metric = simplex_tensor(positions, ele)
      ele_metric_vector(1) = ele_metric(1,1)
      ele_metric_vector(2) = ele_metric(2,2)
      ele_metric_vector(3) = ele_metric(1,2)
      nodes => ele_nodes(positions, ele)
      do j=1, size(nodes)
        n = adjacent_locked_element_count(nodes(j))
        ! take running average of all 'element metric's for the elements adjacent to the nodes
        ! note that in the first contribution, n=0, we throw out the metric values that were there originally
        metric(:,nodes(j)) = (n*metric(:,nodes(j))+ele_metric_vector)/(n+1)
        adjacent_locked_element_count(nodes(j)) = n+1
      end do
    end do

  end subroutine relax_metric_locked_regions
  
  subroutine mba2d_integration_check_options
    character(len = *), parameter :: base_path = "/mesh_adaptivity/hr_adaptivity"  
    integer :: dim, stat
    
    if(.not. have_option(base_path)) then
      ! Nothing to check
      return
    end if
  
    call get_option("/geometry/dimension", dim, stat)
    if(stat /= SPUD_NO_ERROR) then
      ! This isn't the place to complain about this error
      return
    else if(have_option(base_path // "/adaptivity_library/libmba2d") .or. dim == 2) then
#ifndef HAVE_MBA_2D
      FLExit("Cannot use libmba2d without the libmba2d library. Reconfigure with --enable-2d-adaptivity")
#endif
#ifndef HAVE_ZOLTAN
      if(isparallel()) then
        ewrite(0, *) "Warning: It is recommended that you use zoltan with libmba2d in parallel. Reconfigure with --with-zoltan"
      end if
#endif

      if((dim /= 2).and.(.not.(have_option(base_path // "/vertically_structured_adaptivity").and.(dim==3)))) then
        FLExit("libmba2d can only be used in 2D or 2+1D")
      end if
    end if
  
  end subroutine mba2d_integration_check_options

end module mba2d_integration
