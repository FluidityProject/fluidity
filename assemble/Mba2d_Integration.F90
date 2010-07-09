#include "fdebug.h"

module mba2d_integration
  use quadrature
  use elements
  use fields
  use state_module
  use interpolation_module
  use meshdiagnostics
  use vtk_interfaces
  use eventcounter
  use node_boundary
  use metric_tools
  use limit_metric_module
  use mba_adapt_module
  use node_locking
  use halos
  use quicksort
  use surface_id_interleaving
  use data_structures
  use global_parameters, only : current_debug_level
  use adapt_integration
#ifdef HAVE_MBA_2D
  use mba2d_module
#endif
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

    integer :: nonods, mxnods, orig_stotel, stotel, mxface, totele, maxele
    real, dimension(:, :), allocatable :: pos
    integer, dimension(:, :), allocatable :: ipf
    integer, dimension(:, :), allocatable :: ipe
    real, dimension(:, :), allocatable :: parcrv
    integer, dimension(:), allocatable :: ipv
    integer, dimension(:), allocatable :: ifv
    integer, dimension(:), allocatable :: iFnc
    integer, dimension(:), allocatable :: lbE
    real, dimension(:, :), allocatable :: tmp_metric
    integer :: i, j, k, partition_surface_id, face
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
    type(halo_type), pointer :: old_halo
    
    ! Surface ID interleaving
    integer :: max_coplanar_id
    integer, dimension(:), allocatable :: boundary_ids, coplanar_ids, surface_ids, mba_boundary_ids
    type(integer_hash_table) :: physical_surface_ids
    ! Element locking
    integer :: nfe
    integer, dimension(:), allocatable :: ife
    integer :: nfv
    logical, dimension(:), allocatable :: is_locked_ele
    integer, dimension(:), pointer :: node_eles
    type(csr_sparsity), pointer :: nelist
    ! Flattened halo data
    integer :: nhalos
    integer, dimension(:), allocatable :: old_sends, old_receives, sends_starts, receives_starts, new_sends, new_receives
    ! Fields we use to find the new halo numbering
    type(scalar_field) :: old_halo_field, new_halo_field
    type(integer_hash_table) :: old_to_new
    type(integer_hash_table) :: input_face_numbering_to_mba2d_numbering

!#define DUMP_HALO_INTERPOLATION
#ifdef DUMP_HALO_INTERPOLATION
    type(mesh_type) :: xmesh_pwc
    type(scalar_field) :: locked_elements
#endif

#ifdef GIVE_LIPNIKOV_OUTPUT
    integer :: rank
    character(len=255) :: filename
#endif

!#define DUMP_HALO
#ifdef DUMP_HALO
    type(scalar_field) :: sends_sfield, receives_sfield
    integer :: proc
#endif

    ewrite(1, *) "In adapt_mesh_mba2d"

    assert(metric%dim == 2)

    xmesh => input_positions%mesh
    call deallocate_boundcount
    call initialise_boundcount(xmesh, input_positions)

    ! mxnods is an option to adaptivity specifying the maximum number of nodes
    xpctel = max(expected_elements(input_positions, metric), 5)
    mxnods = max_nodes(input_positions, expected_nodes(input_positions, xpctel, global = .false.))

    nonods = node_count(xmesh)
    totele = ele_count(xmesh)
    orig_stotel = surface_element_count(xmesh)
    mxface = int(max((float(mxnods) / float(nonods)) * orig_stotel * 3.5, 10000.0))
    maxele = int(max((float(mxnods) / float(nonods)) * totele * 1.5, 10000.0))
    maxp = mxnods * 1.2

    allocate(pos(2, maxp))
    pos = 0.0
    do i=1,2
      pos(i, 1:nonods) = input_positions%val(i)%ptr
    end do

    allocate(surface_ids(surface_element_count(xmesh)))
    call interleave_surface_ids(xmesh, surface_ids, max_coplanar_id)

    allocate(ipf(4, mxface))
    ipf = 0
    partition_surface_id = maxval(surface_ids) + 1
    stotel = 0

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
          if (face <= orig_stotel) then ! if face is genuinely external, i.e. on the domain exterior
            ipf(4, stotel) = surface_ids(face)
          else
            ipf(4, stotel) = partition_surface_id
          end if
        end if
      end do
    end do

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
    assert(stotel < mxface)

    allocate(ipe(3, maxele))
    ipe = 0
    do i=1,totele
      ipe(:, i) = ele_nodes(xmesh, i)
    end do

    ! I hope mba2d doesn't mind the same node being locked twice.
    call get_locked_nodes(input_positions, locked_nodes)
    npv = count(boundcount > 1) + size(locked_nodes)
    allocate(ipv(npv))
    j = 1
    do i=1,nonods
      if (boundcount(i) > 1) then
        ipv(j) = i
        j = j + 1
      end if
    end do
    ! So now j = count(boundcount > 1) + 1
    ipv(j:) = locked_nodes
    deallocate(locked_nodes)

    allocate(parcrv(2, mxface))
    parcrv = 0

    nhalos = halo_count(xmesh)
    assert(any(nhalos == (/0, 1, 2/)))
    if(nhalos > 0) then
      nelist => extract_nelist(xmesh)
    
      old_halo => xmesh%halos(nhalos)

      allocate(ife(totele))
      nfe = 0
      allocate(is_locked_ele(totele))
      is_locked_ele = .false.
      do i = 1, halo_proc_count(old_halo)
        do j = 1, halo_send_count(old_halo, i)
          node_eles => row_m_ptr(nelist, halo_send(old_halo, i, j))
          do k = 1, size(node_eles)
            if(is_locked_ele(node_eles(k))) cycle
            nfe = nfe + 1
            ife(nfe) = node_eles(k)
            is_locked_ele(node_eles(k)) = .true.
          end do
        end do
        do j = 1, halo_receive_count(old_halo, i)
          node_eles => row_m_ptr(nelist, halo_receive(old_halo, i, j))
          do k = 1, size(node_eles)
            if(is_locked_ele(node_eles(k))) cycle
            nfe = nfe + 1
            ife(nfe) = node_eles(k)
            is_locked_ele(node_eles(k)) = .true.
          end do
        end do
      end do
      deallocate(is_locked_ele)
    else
      nfe = 0
      allocate(ife(nfe))
    end if

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
      allocate(boundary_ids(stotel))
      allocate(coplanar_ids(stotel))
      call deinterleave_surface_ids(ipf(4, 1:stotel), max_coplanar_id, boundary_ids, coplanar_ids)  
      
      allocate(new_sndgln(1:2,1:stotel))
      new_sndgln=ipf(1:2,1:stotel)
      call add_faces(new_mesh, sndgln=reshape(new_sndgln, (/ 2*stotel /) ), &
        boundary_ids=boundary_ids)
      deallocate(boundary_ids)
      deallocate(new_sndgln)
      
      if(associated(xmesh%faces%coplanar_ids)) then
        allocate(new_mesh%faces%coplanar_ids(stotel))
        new_mesh%faces%coplanar_ids = coplanar_ids
      end if
      deallocate(coplanar_ids)
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

      stotel = j - 1
      allocate(mba_boundary_ids(stotel))

      do i=1,stotel
        mba_boundary_ids(i) = ipf(4, fetch(physical_surface_ids, i))
      end do

      allocate(boundary_ids(stotel))
      allocate(coplanar_ids(stotel))
      call deinterleave_surface_ids(mba_boundary_ids, max_coplanar_id, boundary_ids, coplanar_ids)  
      deallocate(mba_boundary_ids)
      
      allocate(new_sndgln(1:2,1:stotel))
      do i=1,stotel
        new_sndgln(1:2, i) = ipf(1:2, fetch(physical_surface_ids, i))
      end do
      call add_faces(new_mesh, sndgln=reshape(new_sndgln, (/ 2*stotel /) ), &
        boundary_ids=boundary_ids)
      deallocate(boundary_ids)
      deallocate(new_sndgln)
      
      if(associated(xmesh%faces%coplanar_ids)) then
        allocate(new_mesh%faces%coplanar_ids(stotel))
        new_mesh%faces%coplanar_ids = coplanar_ids
      end if
      deallocate(coplanar_ids)

      call deallocate(physical_surface_ids)
    end if
    
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
      ! Flatten the old halo
      allocate(sends_starts(halo_proc_count(old_halo) + 1))
      allocate(receives_starts(halo_proc_count(old_halo) + 1))
      allocate(old_sends(halo_all_sends_count(old_halo)))
      allocate(old_receives(halo_all_receives_count(old_halo)))
      call extract_raw_halo_data(old_halo, old_sends, sends_starts, old_receives, receives_starts)

      ! Sorry about this. Blame James and Patrick ...
      ! Transfer the halo nodes by interpolating their flattened indices
      call allocate(old_halo_field, xmesh, name = "OldHaloTransfer")
      call zero(old_halo_field)
      do i = 1, size(old_sends)
        call set(old_halo_field, old_sends(i), float(old_sends(i)))
      end do
      do i = 1, size(old_receives)
        call set(old_halo_field, old_receives(i), float(old_receives(i)))
      end do
      call allocate(new_halo_field, output_positions%mesh, name = "NewHaloTransfer")
      call linear_interpolation(old_halo_field, input_positions, new_halo_field, output_positions)
#ifdef DUMP_HALO_INTERPOLATION
      input_positions%mesh%halos => null()
      xmesh_pwc = piecewise_constant_mesh(input_positions%mesh, "XMeshPWC")
      call allocate(locked_elements, xmesh_pwc, "LockedElements")
      call zero(locked_elements)
      do i=1,nfe
        call set(locked_elements, ife(i), 1.0)
      end do
      call vtk_write_fields("halo_interpolation", 0, input_positions, input_positions%mesh, sfields=(/old_halo_field, locked_elements/))
      call vtk_write_fields("halo_interpolation", 1, output_positions, output_positions%mesh, sfields=(/new_halo_field/))
      call deallocate(locked_elements)
      call deallocate(xmesh_pwc)
#endif
      call deallocate(old_halo_field)

      ! Set the new halo numbering
      call allocate(old_to_new)
      do i = 1, node_count(new_halo_field)
        j = floor(node_val(new_halo_field, i) + 0.5)
        if(j > 0) then
          call insert(old_to_new, j, i)
        end if
      end do

      allocate(new_sends(size(old_sends)))
      allocate(new_receives(size(old_receives)))
      do i=1,size(old_sends)
        new_sends(i) = fetch(old_to_new, old_sends(i))
      end do
      do i=1,size(old_receives)
        new_receives(i) = fetch(old_to_new, old_receives(i))
      end do

      call deallocate(new_halo_field)
      call deallocate(old_to_new)

      ! Allocate and set the new halo
      allocate(output_positions%mesh%halos(nhalos))
      call form_halo_from_raw_data(output_positions%mesh%halos(nhalos), size(sends_starts) - 1, new_sends, &
        & sends_starts, new_receives, receives_starts, nowned_nodes = nonods - size(new_receives), &
        & ordering_scheme = HALO_ORDER_GENERAL, create_caches = (nhalos == 1))

      deallocate(sends_starts)
      deallocate(receives_starts)
      deallocate(old_sends)
      deallocate(old_receives)
      deallocate(new_sends)
      deallocate(new_receives)

      if(nhalos == 2) then
        ! Derive remaining halos
        call derive_l1_from_l2_halo(output_positions%mesh, ordering_scheme = HALO_ORDER_GENERAL, create_caches = .false.)
        ! Reorder the nodes for trailing receives consistency
        call renumber_positions_trailing_receives(output_positions)

        allocate(output_positions%mesh%element_halos(2))
        ! Reorder the elements for trailing receives consistency
        call derive_element_halo_from_node_halo(output_positions%mesh, &
          & ordering_scheme = HALO_ORDER_GENERAL, create_caches = .false.)
        call renumber_positions_elements_trailing_receives(output_positions)
      else
        ! Reorder the nodes for trailing receives consistency
        call renumber_positions_trailing_receives(output_positions)
        
        allocate(output_positions%mesh%element_halos(0))
      end if

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
    FLAbort("You called mba_adapt without the mba2d library. Reconfigure with --enable-2d-adaptivity")
#endif
  end subroutine adapt_mesh_mba2d
  
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

      if(dim /= 2) then
        FLExit("libmba2d can only be used in 2D")
      end if
    end if
  
  end subroutine mba2d_integration_check_options

end module mba2d_integration
