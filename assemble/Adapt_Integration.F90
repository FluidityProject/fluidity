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

module adapt_integration

  use quadrature
  use elements
  use fldebug
  use fields
  use halos
  use limit_metric_module
  use meshdiagnostics
  use node_locking
  use spud
  use surface_id_interleaving
  use tictoc

  implicit none
  
  private
  
  public :: adapt_mesh, adapt_integration_check_options
  
  character(len = *), parameter :: base_path = "/mesh_adaptivity/hr_adaptivity"
  
#ifdef HAVE_ADAPTIVITY
  interface
    subroutine adaptmem(nnod, nelm, szenls, nselm, totfre, &
      & xpctel, xpctnd, xpctse, &
      & metric, szint, szrl)
      implicit none
      integer, intent(in) :: nnod
      integer, intent(in) :: nelm
      integer, intent(in) :: szenls
      integer, intent(in) :: nselm
      integer, intent(in) :: totfre
      integer, intent(in) :: xpctel
      integer, intent(in) :: xpctnd
      integer, intent(in) :: xpctse
      logical, intent(in) :: metric
      integer, intent(out) :: szint
      integer, intent(out) :: szrl
    end subroutine adaptmem
  
    subroutine adptvy(intarr, intsiz, rlarr,  rlsiz, &
      & geom3d, srfgmy, useq, &
      & nnod,   nelm,   nselm,  absolutemxnods, &
      & szenls, enlbas, enlist, elmreg, &
      & clcgmy, szsnls, snlbas, snlist, surfid, &
      & prdnds, nprdnd, &
      & nodx,   nody,   nodz, &
      & intnnd, intnel, intszl, intenl, intenb, &
      & intndx, intndy, intndz, &
      & orgmtx, oldfld, nfree,  totfre, nfield, &
      & xpctel, nwnnod, nwnelm, nwnsel, &
      & nwszen, nwszsn, nwsznn, nwndlc, nwsrow, &
      & nwenlb, nwenls, nwsnlb, nwsnls, nwsfid, &
      & nwelrg, nwnodx, nwnody, nwnodz, &
      & newmtx, newfld, &
      & biglst, nodlst, &
      & dotop,  minchg, nsweep, mshopt, twostg, togthr, &
      & gather, scater, ngath,  nhalo,  pnod, &
      & atosen, atorec, nproc, debug_level, dbg, chcnsy)
      implicit none
      integer, intent(in) :: intsiz
      integer, dimension(intsiz), intent(out) :: intarr
      integer, intent(in) :: rlsiz
      real, dimension(rlsiz), intent(out) :: rlarr
      logical, intent(in) :: geom3d
      logical, intent(in) :: srfgmy
      logical, intent(in) :: useq
      integer, intent(in) :: nnod
      integer, intent(in) :: nelm
      integer, intent(in) :: nselm
      integer, intent(in) :: absolutemxnods
      integer, intent(in) :: szenls
      integer, dimension(nelm + 1), intent(in) :: enlbas
      integer, dimension(nelm * 4), intent(in) :: enlist
      integer, dimension(nelm), intent(in) :: elmreg
      logical, intent(in) :: clcgmy
      integer, intent(in) :: szsnls
      integer, dimension(nselm + 1), intent(in) :: snlbas
      integer, dimension(nselm * 3), intent(in) :: snlist
      integer, dimension(nselm), intent(in) :: surfid
      integer, intent(in) :: nprdnd
      integer, dimension(nprdnd), intent(in) :: prdnds
      real, dimension(nnod), intent(in) :: nodx
      real, dimension(nnod), intent(in) :: nody
      real, dimension(nnod), intent(in) :: nodz
      integer, intent(in) :: intnnd
      integer, intent(in) :: intnel
      integer, intent(in) :: intszl
      integer, dimension(intnel * 4), intent(in) :: intenl
      integer, dimension(intnel + 1), intent(in) :: intenb
      real, dimension(intnnd), intent(in) :: intndx
      real, dimension(intnnd), intent(in) :: intndy
      real, dimension(intnnd), intent(in) :: intndz
      real, dimension(nnod * 9), intent(in) :: orgmtx
      real, dimension(nfield * nnod), intent(in) :: oldfld
      integer, intent(in) :: nfield
      integer, intent(in) :: xpctel
      integer, dimension(nfield), intent(in) :: nfree
      integer, intent(in) :: totfre
      integer, intent(out) :: nwnnod
      integer, intent(out) :: nwnelm
      integer, intent(out) :: nwnsel
      integer, intent(out) :: nwszen
      integer, intent(out) :: nwszsn
      integer, intent(out) :: nwsznn
      integer, intent(inout) :: nwndlc
      integer, intent(inout) :: nwsrow
      integer, intent(out) :: nwenlb
      integer, intent(out) :: nwenls
      integer, intent(out) :: nwsnlb
      integer, intent(out) :: nwsnls
      integer, intent(out) :: nwsfid
      integer, intent(out) :: nwelrg
      integer, intent(out) :: nwnodx
      integer, intent(out) :: nwnody
      integer, intent(out) :: nwnodz
      integer, intent(inout) :: newmtx
      integer, intent(out) :: newfld
      integer, intent(out) :: biglst
      integer, intent(out) :: nodlst
      real, intent(in) :: dotop
      real, intent(in) :: minchg
      integer, intent(in) :: nsweep
      logical, dimension(6), intent(in) :: mshopt
      logical, intent(in) :: twostg
      logical, intent(in) :: togthr
      integer, intent(in) :: ngath
      integer, dimension(ngath), intent(inout) :: gather
      integer, intent(in) :: nhalo
      integer, dimension(nhalo), intent(inout) :: scater
      integer, intent(inout) :: pnod
      integer, intent(in) :: nproc
      integer, dimension(nproc + 1), intent(in) :: atosen
      integer, dimension(nproc + 1), intent(in) :: atorec
      integer, intent(in) :: debug_level
      logical, intent(in) :: dbg
      logical, intent(in) :: chcnsy
    end subroutine adptvy
  end interface
#endif
  
contains

  subroutine adapt_mesh(input_positions, metric, output_positions, node_ownership, &
                                                              force_preserve_regions)
    !!< Adapt the supplied input mesh using libadaptivity. Return the new
    !!< adapted mesh in output_positions (which is allocated by this routine).
    
    type(vector_field), intent(in) :: input_positions
    type(tensor_field), intent(in) :: metric
    type(vector_field), target, intent(out) :: output_positions
    !! Map from new nodes to old elements. Allocated by this routine.
    integer, dimension(:), pointer, optional :: node_ownership
    logical, intent(in), optional :: force_preserve_regions
    
    ! Linear tets only
    integer, parameter :: dim = 3, nloc = 4, snloc = 3
    
    ! adaptmem arguments
    integer :: nnod, nelm, szenls, nselm, totfre, xpctel, xpctnd, xpctse
    logical :: have_metric
    
    ! adptvy arguments
    ! Working memory
    integer, dimension(:), allocatable :: intarr
    integer :: intsiz
    real, dimension(:), allocatable :: rlarr
    integer :: rlsiz
    ! Input variables
    logical :: geom3d, srfgmy, useq
    integer :: absolutemxnods
    integer, dimension(:), allocatable, target :: enlbas
    integer, dimension(:), pointer :: enlist
    integer, dimension(:), allocatable :: elmreg
    logical :: clcgmy
    integer :: szsnls
    integer, dimension(:), allocatable :: snlbas, snlist, surfid, prdnds
    integer :: nprdnd
    real, dimension(:), pointer :: nodx, nody, nodz
    integer :: intnnd, intnel, intszl
    integer, dimension(:), pointer :: intenl, intenb
    real, dimension(:), pointer :: intndx, intndy, intndz
    real, dimension(:), allocatable :: orgmtx
    real, dimension(:), allocatable :: oldfld
    integer, dimension(:), allocatable :: nfree
    integer :: nfield
    ! Output variables
    integer :: nwnnod, nwnelm, nwnsel, nwszen, nwszsn, nwsznn, &
      & nwndlc, nwsrow, nwenlb, nwenls, nwsnlb, nwsnls, nwsfid, nwelrg, &
      & nwnodx, nwnody, nwnodz, newmtx, newfld, biglst, nodlst
    ! More input variables
    real :: dotop, minchg
    integer :: nsweep
    logical, dimension(6) :: mshopt
    logical :: twostg
    logical :: togthr
    integer, dimension(:), allocatable :: gather, scater
    integer :: ngath, nhalo, pnod
    integer, dimension(:), allocatable :: atosen, atorec
    integer :: nproc, debug_level
    logical :: dbg, chcnsy
  
    integer :: i, max_coplanar_id, nhalos, stat
    integer, dimension(:), allocatable :: boundary_ids, coplanar_ids
    real :: mestp1
    type(halo_type), pointer :: old_halo
    type(mesh_type), pointer :: output_mesh
    
    ! Buffer factor to emulate behaviour of legacy expected elements function
    real, parameter :: expected_elements_buffer = 1.2  
    ! Buffer factor for max. nodes
    real, parameter :: mxnods_buffer = 1.5
    
    ewrite(1, *) "In adapt_mesh"
    
#ifdef DDEBUG
    assert(input_positions%dim == dim)
    assert(ele_loc(input_positions, 1) == nloc)
    if(surface_element_count(input_positions) > 0) then
      assert(associated(input_positions%mesh%faces))
      assert(face_loc(input_positions, 1) == snloc)
    end if
    assert(metric%mesh == input_positions%mesh)
    if(present(node_ownership)) then
      assert(.not. associated(node_ownership))
    end if
#endif
    
    ewrite(2, *) "Forming adaptmem arguments"
    
    nnod = node_count(input_positions)  ! Number of nodes
    nelm = element_count(input_positions)  ! Number of volume elements
    szenls = nloc * nelm  ! Size of the volume element list
    nselm = surface_element_count(input_positions)  ! Number of surface elements
    totfre = 0  ! Number of fields
    xpctel = expected_elements(input_positions, metric) * expected_elements_buffer  ! Expected number of volume elements
    xpctnd = -1  ! Expected number of nodes
    xpctse = -1  ! Expected number of surface elements
    have_metric = .true.  ! Unknown
    
    ! Initialise output variables, just in case they're also used as input
    intsiz = 0  ! Integer working memory size
    rlsiz = 0  ! Real working memory size
    
    ewrite(1, *) "Calling adaptmem from adapt_mesh"
#ifdef HAVE_ADAPTIVITY
    call adaptmem(nnod, nelm, szenls, nselm, totfre, &
      & xpctel, xpctnd, xpctse, &
      & have_metric, intsiz, rlsiz)
#else
    FLAbort("Fluidity compiled without libadaptivity support")
#endif
    ewrite(1, *) "Exited adaptmem"
    
    ewrite(2, "(a,i0)") "Integer working memory size: ", intsiz
    ewrite(2, "(a,i0)") "Real working memory size: ", rlsiz
    if(intsiz<0) then
      FLAbort("Invalid integer working memory size")
    end if
    if(rlsiz<0) then
      FLAbort("Invalid real working memory size")
    end if
      
    ewrite(2, *) "Forming remaining adptvy arguments"
    
    ! Working memory
    allocate(intarr(intsiz))  ! Integer working memory
    allocate(rlarr(rlsiz))  ! Real working memory
    
    geom3d = (mesh_dim(input_positions) == 3)  ! Whether the domain is 3D
    srfgmy = .false.  ! Whether the surface mesh should be kept intact during the adapt
    useq = .false.  ! Unknown
    
    ! Maximum number of nodes
    call get_option(base_path // "/maximum_number_of_nodes", absolutemxnods, default = 100000)
    if(isparallel()) then
      if(.not. have_option(base_path // "/maximum_number_of_nodes/per_process")) then
        absolutemxnods = absolutemxnods / getnprocs()
      end if
      absolutemxnods = &
        & max( &
          & absolutemxnods, &
          & expected_nodes(input_positions, int(xpctel / expected_elements_buffer), global = .false.) &
        & )
    end if
    absolutemxnods = max(absolutemxnods, nnod) * mxnods_buffer
    ewrite(2, "(a,i0)") "Max. nodes: ", absolutemxnods
    
    ! Volume element list
    allocate(enlbas(nelm + 1))
    do i = 1, nelm + 1
      enlbas(i) = (i - 1) * nloc
    end do
    enlist => input_positions%mesh%ndglno
    
    ! Region IDs
    allocate(elmreg(nelm))
    if(associated(input_positions%mesh%region_ids).and.&
      (have_option(base_path // "/preserve_mesh_regions").or.present_and_true(force_preserve_regions))) then
      elmreg = input_positions%mesh%region_ids
    else
      elmreg = 0
    end if
    
    ! Surface element list
    clcgmy = .true.  ! Is .true. if the geometry should be calculated, and ignore snlist
    szsnls = nselm * snloc
    allocate(snlbas(nselm + 1))
    do i = 1, nselm + 1
       snlbas(i) = (i - 1) * snloc
    end do
    allocate(snlist(nselm * snloc))
    if(nselm > 0) then
      call getsndgln(input_positions%mesh, snlist)
    end if
    
    ! Surface IDs
    allocate(surfid(nselm))
    call interleave_surface_ids(input_positions%mesh, surfid, max_coplanar_id)
    
    ! Node locking
    call get_locked_nodes(input_positions, prdnds)
    nprdnd = size(prdnds)
    
    ! Coordinates
    nodx => input_positions%val(1)%ptr
    nody => input_positions%val(2)%ptr
    nodz => input_positions%val(3)%ptr
    
    ! Interpolation mesh (the same as the input mesh)
    intnnd = nnod
    intnel = nelm
    intszl = szenls
    intenl => enlist
    intenb => enlbas
    intndx => nodx
    intndy => nody
    intndz => nodz
    
    ! Metric
    allocate(orgmtx(nnod * dim ** 2))
    select case(metric%field_type)
      case(FIELD_TYPE_NORMAL)
        orgmtx = reshape(metric%val, (/nnod * dim ** 2/))
      case(FIELD_TYPE_CONSTANT)
        do i = 1, nnod
          orgmtx((i - 1) * dim * dim + 1:i * dim * dim) = reshape(node_val(metric, i), (/dim ** 2/))
        end do
      case default
        ewrite(-1, *) "For field type ", metric%field_type
        FLAbort("Unrecognised field type")
    end select

    ! Field data - none, as we don't use libadaptivity for interpolation any
    ! more
    nfield = 0  ! Number of fields
    allocate(oldfld(nfield * nnod))
    allocate(nfree(nfield))
    
    call get_option(base_path // "/functional_tolerance", mestp1, default = 0.0)
    dotop =  max(abs(mestp1), 0.15)  ! Functional tolerance
    minchg = 0.01  ! Unknown
    
    call get_option(base_path // "/adaptivity_library/libadaptivity/sweeps/", &
      & nsweep, default = 10)  ! Number of adapt sweeps
    
    ! Which element operations are we using?
    mshopt(1) = .true.  ! Split edges if true 
    mshopt(2) = .true.  ! Collapse edges if true
    mshopt(3) = .true.  ! Perform edge to face and edge to edge swapping if true
    mshopt(4) = .true.  ! Perform face to edge swapping if true
    mshopt(5) = .true.  ! Split elements (do not use this yet)
                        ! In fact, this option is currently ignored and element
                        ! splitting is not performed by libadaptivity
    mshopt(6) = .true.  ! Move nodes if true 
    
    twostg = .false.  ! Two stages of adapting, with no refinement on first
    togthr = .true.  ! Lumps node movement adaptivity in with connectivity
                     ! changes

    ! Parallel data
    nhalos = halo_count(input_positions)
    assert(any(nhalos == (/0, 1, 2/)))
    if(nhalos > 0) then
      old_halo => input_positions%mesh%halos(nhalos)
      assert(trailing_receives_consistent(old_halo))
      
      nproc = halo_proc_count(old_halo)
      ngath = halo_all_sends_count(old_halo)
      allocate(gather(ngath))
      allocate(atosen(nproc + 1))
      nhalo = halo_all_receives_count(old_halo)
      allocate(scater(nhalo))
      allocate(atorec(nproc + 1))
      call extract_raw_halo_data(old_halo, gather, atosen, scater, atorec, nowned_nodes = pnod)
    else
      nproc = 1
      ngath = 0
      allocate(gather(ngath))
      allocate(atosen(1))
      atosen = 0
      nhalo = 0
      allocate(scater(nhalo))
      allocate(atorec(1))
      atorec = 0
      pnod = nnod
    end if
    
    ! Debugging options
    debug_level = current_debug_level ! Verbosity
    dbg = .false.     ! Enable additional run time debugging within adaptivity.
    chcnsy = .false.  ! This option was for further run time
                      ! consistency checks within adaptivity but it's
                      ! currently disabled.
        
    ! Initialise output variables, just in case they're also used as input
    nwnnod = 0  ! Number of nodes
    nwnelm = 0  ! Number of volume elements
    nwnsel = 0  ! Number of surface elements
    nwszen = 0  ! Size of the volume element list
    nwszsn = 0  ! Size of the surface element list
    nwsznn = 0  ! Unknown
    !nwndlc = 0  ! Node ownership list (start index in intarr)
    !nwsrow = 0  ! Surface element ownership list (start index in intarr)
    nwenlb = 0  ! Unknown
    nwenls = 0  ! Volume element numbering list (start index in intarr)
    nwsnlb = 0  ! Unknown
    nwsnls = 0  ! Surface element numbering list (start index in intarr)
    nwsfid = 0  ! Surface IDs (start index in intarr)
    nwelrg = 0  ! Region IDs (start index in intarr)
    nwnodx = 0  ! x-coordinates (start index in rlarr)
    nwnody = 0  ! y-coordinates (start index in rlarr)
    nwnodz = 0  ! z-coordinates (start index in rlarr)
    !newmtx = 0  ! Adaptivity metric (start index in rlarr)
    newfld = 0  ! Unknown
    biglst = 0  ! Unknown
    nodlst = 0  ! Unknown
    
    ! Output options
    newmtx = -1  ! Do not return interpolated metric
    if(present(node_ownership)) then
      nwndlc = 1  ! Return map from new nodes to old elements
    else
      nwndlc = -1  ! Do not return map from new nodes to old elements
    end if
    nwsrow = -1  ! Do not return surface element owners
    
    ewrite(1, *) "Calling adptvy from adapt_mesh"
    call tic(TICTOC_ID_SERIAL_ADAPT)
#ifdef HAVE_ADAPTIVITY
    call adptvy(intarr, intsiz, rlarr,  rlsiz, &
      & geom3d, srfgmy, useq, &
      & nnod,   nelm,   nselm,  absolutemxnods, &
      & szenls, enlbas, enlist, elmreg, &
      & clcgmy, szsnls, snlbas, snlist, surfid, &
      & prdnds, nprdnd, &
      & nodx,   nody,   nodz, &
      & intnnd, intnel, intszl, intenl, intenb, &
      & intndx, intndy, intndz, &
      & orgmtx, oldfld, nfree,  totfre, nfield, &
      & xpctel, nwnnod, nwnelm, nwnsel, &
      & nwszen, nwszsn, nwsznn, nwndlc, nwsrow, &
      & nwenlb, nwenls, nwsnlb, nwsnls, nwsfid, &
      & nwelrg, nwnodx, nwnody, nwnodz, &
      & newmtx, newfld, &
      & biglst, nodlst, &
      & dotop,  minchg, nsweep, mshopt, twostg, togthr, &
      & gather, scater, ngath,  nhalo,  pnod, &
      & atosen, atorec, nproc, debug_level, dbg, chcnsy)
#else
    FLAbort("Fluidity compiled without libadaptivity support")
#endif
    call toc(TICTOC_ID_SERIAL_ADAPT)
    ewrite(1, *) "Exited adptvy"

    if(nwnnod < 0) then
      FLAbort("Mesh adaptivity exited with an error")
    end if
    assert(nwnnod <= absolutemxnods)
    assert(nwnelm >= 0)
    
    deallocate(orgmtx)
    deallocate(enlbas)
    deallocate(elmreg)
    deallocate(snlbas)
    deallocate(snlist)
    deallocate(surfid)
    deallocate(prdnds)
    deallocate(oldfld)
    deallocate(nfree)
        
    ewrite(2, *) "Constructing output positions"
        
    allocate(output_mesh)
    call allocate(output_mesh, nwnnod, nwnelm, input_positions%mesh%shape, name = input_positions%mesh%name)
    output_mesh%shape%refcount%tagged = .false.
    output_mesh%shape%quadrature%refcount%tagged = .false.
    
    output_mesh%ndglno = intarr(nwenls:nwenls + nwszen - 1)
    output_mesh%option_path = input_positions%mesh%option_path    

    ! Construct the new positions
    call allocate(output_positions, dim, output_mesh, name = input_positions%name)
    call deallocate(output_mesh)
    deallocate(output_mesh)
    output_mesh => output_positions%mesh
  
    call set_all(output_positions, 1, rlarr(nwnodx:nwnodx + nwnnod - 1))
    call set_all(output_positions, 2, rlarr(nwnody:nwnody + nwnnod - 1))
    call set_all(output_positions, 3, rlarr(nwnodz:nwnodz + nwnnod - 1))
    output_positions%option_path = input_positions%option_path

    if(nhalos > 0) then
      ewrite(2, *) "Constructing output halos"
      
      allocate(output_mesh%halos(nhalos)) 
      call form_halo_from_raw_data(output_mesh%halos(nhalos), nproc, gather, atosen, scater, atorec,&
           & nowned_nodes = nwnnod - nhalo, create_caches = .true.) 

      if(nhalos == 2) then
        ! Derive remaining halos        
        call derive_l1_from_l2_halo(output_mesh, &
          & ordering_scheme = HALO_ORDER_TRAILING_RECEIVES, create_caches = .true.)
          
        allocate(output_mesh%element_halos(2))
        call derive_element_halo_from_node_halo(output_mesh, &
          & ordering_scheme = HALO_ORDER_GENERAL, create_caches = .false.)
        call renumber_positions_elements_trailing_receives(output_positions)
      else
        allocate(output_mesh%element_halos(1))
        call derive_element_halo_from_node_halo(output_mesh, &
          & ordering_scheme = HALO_ORDER_GENERAL, create_caches = .false.)
        call renumber_positions_elements_trailing_receives(output_positions)
      end if
      ! Adaptivity is not guaranteed to return halo elements in the same
      ! order in which they went in. We therefore need to fix this order.
      call reorder_element_numbering(output_positions)
              
#ifdef DDEBUG
      do i = 1, nhalos
        assert(trailing_receives_consistent(output_mesh%halos(i)))
        assert(halo_valid_for_communication(output_mesh%halos(i)))
        assert(halo_verifies(output_mesh%halos(i), output_positions))
      end do
#endif
      
      ewrite(2, *) "Finished constructing output halos"
    end if

    deallocate(gather)
    deallocate(atosen)
    deallocate(scater)
    deallocate(atorec)
  
    ewrite(2, *) "Constructing output surface data"
    
    allocate(boundary_ids(nwnsel))
    allocate(coplanar_ids(nwnsel))
    call deinterleave_surface_ids(intarr(nwsfid:nwsfid + nwnsel - 1), max_coplanar_id, boundary_ids, coplanar_ids)
    if(nhalos > 0) then
      call add_faces(output_mesh, sndgln = intarr(nwsnls:nwsnls + nwszsn - 1), boundary_ids = boundary_ids, private_nodes = 0)
    else
      call add_faces(output_mesh, sndgln = intarr(nwsnls:nwsnls + nwszsn - 1), boundary_ids = boundary_ids)
    end if
    deallocate(boundary_ids)
    if(associated(input_positions%mesh%faces%coplanar_ids)) then
      allocate(output_mesh%faces%coplanar_ids(nwnsel))
      output_mesh%faces%coplanar_ids = coplanar_ids
    end if
    deallocate(coplanar_ids)
          
    ewrite(2, *) "Finished constructing output surface data"  
    
    if(have_option(base_path // "/preserve_mesh_regions")&
              .or.present_and_true(force_preserve_regions)) then
      allocate(output_mesh%region_ids(nwnelm))
      output_mesh%region_ids = intarr(nwelrg:nwelrg + nwnelm - 1)
    end if  
    
#ifdef DDEBUG
    call verify_positions(output_positions)
#endif
        
    ewrite(2, *) "Finished constructing output positions"
    
    if(present(node_ownership)) then
      ! Return the node ownership
      allocate(node_ownership(nwnnod))
      node_ownership = intarr(nwndlc:nwndlc + nwnnod - 1)
    end if
    
    deallocate(intarr)
    deallocate(rlarr)
    
    ewrite(1, *) "Exiting adapt_mesh"
    
  end subroutine adapt_mesh
  
  subroutine verify_positions(positions)
    !!< Verify the supplied Coordinate field - replaces elementsok
    
    type(vector_field), intent(in) :: positions
        
    integer :: i, nnodes
    logical :: positive_volumes
    real :: volume
    type(element_type), pointer :: shape
    
    ewrite(1, *) "In verify_positions"
    
    nnodes = node_count(positions)
    
    do i = 1, ele_count(positions)
      assert(all(ele_nodes(positions, i) >= 1))
      assert(all(ele_nodes(positions, i) <= nnodes))
      
      shape => ele_shape(positions, i)
      if(positions%dim == 3 .and. shape%loc == 4 .and. shape%degree == 1) then
        volume = simplex_volume(positions, i)
        if(abs(volume) < epsilon(0.0)) then
          ewrite(-1, "(a,i0)") "For element: ", i
          FLAbort("Degenerate tetrahedron encountered")
        end if

        if(i > 1) then
          if(.not. positive_volumes .eqv. (volume > 0.0)) then
            FLAbort("Signs of tetrahredon volumes are not consistent")
          end if
        else
          positive_volumes = volume > 0.0
        end if
      end if
    end do
    
    ewrite(1, *) "Exiting verify_positions"

  end subroutine verify_positions
  
  subroutine adapt_integration_check_options
    !!< Checks libadaptivity integration related options
    
    integer :: dim, max_nodes, stat
    
    if(.not. have_option(base_path)) then
      ! Nothing to check
      return
    end if
  
    call get_option("/geometry/dimension", dim, stat)
    if(stat /= SPUD_NO_ERROR) then
      ! This isn't the place to complain about this error
      return
    else if(have_option(base_path // "/adaptivity_library/libadaptivity") .or. dim == 3) then
      if(dim /= 3) then
        FLExit("libadaptivity can only be used in 3D")
      end if
    end if
    
    ewrite(2, *) "Checking hr-adaptivity related options"
    
    call get_option(base_path // "/maximum_number_of_nodes", max_nodes, stat)
    if(stat /= SPUD_NO_ERROR) then
      FLExit("Maximum number of nodes required for 3d adaptivity with libadaptivity")
    else if(max_nodes <= 0) then
      FLExit("Maximum number of nodes must be positive")
    end if
    
    ewrite(2, *) "Finished checking hr-adaptivity related options"
  
  end subroutine adapt_integration_check_options

end module adapt_integration
