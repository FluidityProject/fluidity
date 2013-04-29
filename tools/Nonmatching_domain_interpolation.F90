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
#include "confdefs.h"



subroutine Nonmatching_domain_interpolation(c_input_basename_1, input_basename_1_len, c_input_basename_2, &
                       & input_basename_2_len, c_input_mesh_format, input_mesh_format_len, ifinder) bind(c)
  !!< Computes the volume fraction of input mesh 2 on input mesh 1

  use global_parameters, only: OPTION_PATH_LEN, is_active_process, no_active_processes, topology_mesh_name
  use fields
  use parallel_tools, only: isparallel, parallel_filename
  use halos_registration, only: read_halos, write_halos
  use fields_halos, only: verify_consistent_local_element_numbering
  use mesh_files
  use state_module
  use diagnostic_variables
  use diagnostic_fields_wrapper
  use diagnostic_fields_new, only : &
    & calculate_diagnostic_variables_new => calculate_diagnostic_variables
  use populate_state_module
  use iso_c_binding

  ! For projection routine:
  use supermesh_construction
  use supermesh_assembly
  use solvers
  use global_parameters, only: OPTION_PATH_LEN
  use linked_lists
  use bound_field_module
  use fefields, only: compute_lumped_mass, project_field
  use tetrahedron_intersection_module
  use Petsc_Tools
  use petsc
  use vtk_interfaces

  implicit none

  character(kind=c_char, len=1) :: c_input_basename_1(*)
  integer(kind=c_size_t), value :: input_basename_1_len
  character(kind=c_char, len=1) :: c_input_basename_2(*)
  integer(kind=c_size_t), value :: input_basename_2_len
  character(kind=c_char, len=1) :: c_input_mesh_format(*)
  integer(kind=c_size_t), value :: input_mesh_format_len
  integer(kind=c_int), value :: ifinder


  character(len=input_basename_1_len):: input_basename_1
  character(len=input_basename_2_len):: input_basename_2
  character(len=input_mesh_format_len):: input_mesh_format

  integer :: nprocs
  type(state_type), dimension(:), pointer :: states
  type(mesh_type) :: mesh1, mesh2
  type(vector_field) :: position1, position2
  type(scalar_field) :: alpha_source, alpha_target
  type(vector_field), pointer :: pos_source, pos_target

  real :: alpha_source_int, alpha_target_int

  ! Rendezvous bboxes:
  ! dim x 2 x nprocs (or number of partitions) ! 2 because of min/max value of the box in each dimension
  real, dimension(:, :, :), allocatable :: bboxes_src, bboxes_trg
  real, dimension(:, :), allocatable :: my_bbox, my_bbox_src, my_bbox_trg
  integer :: dim
  character(len=255) :: rank_str
  integer :: rank, rank_digits

  integer :: quad_degree
  integer :: i, j, nstates
  integer :: sstat

  logical :: advfront, matching_domain

  ewrite(1, *) "In Nonmatching_domain_interpolation"

  ! Get ranks and nprocesses
  rank = getrank()
  ! printing out stuff is infuriating in fortran. I don't want 30 spaces ..
  rank_digits = max(floor(log10(float(rank))) + 1, 1)
  rank_str(1:rank_digits) = int2str(rank)

  nprocs = getnprocs()
  if (nprocs > 1 .and. input_mesh_format == 'exodusii') then
    FLExit("Nonmatching_domain_interpolation must be run in serial when reading in an ExodusII mesh!")
  end if
  ! now turn into proper fortran strings (is there an easier way to do this?)
  do i=1, input_basename_1_len
    input_basename_1(i:i)=c_input_basename_1(i)
  end do
  do i=1, input_basename_2_len
    input_basename_2(i:i)=c_input_basename_2(i)
  end do
  do i=1, input_mesh_format_len
    input_mesh_format(i:i)=c_input_mesh_format(i)
  end do

  ! Use hard coded values for quad_degree
  quad_degree = 5

  ! Allocating states:
  nstates = 2
  allocate(states(1:nstates))
  do i = 1, nstates
     call nullify(states(i))
     call set_option_path(states(i), "/state["//int2str(i-1)//"]")
  end do

  call initialise_walltime

  ewrite(1,*) "----------------------------------------------"

  ewrite(1,*) "Reading in mesh file: "//trim(input_basename_1)
  ewrite(1,*) "input_mesh_format: "//trim(input_mesh_format)

  ! Read in the mesh file:
  position1=read_mesh_files(trim(input_basename_1), &
                      quad_degree=quad_degree, &
                      format=input_mesh_format)
  ! If reading in a decomposed mesh, read in halos as well:
  if(isparallel()) then
    call read_halos(trim(input_basename_1), position1)
    ! Local element ordering needs to be consistent between processes, otherwise
    ! code in Halos_Repair (used in halo construction of derived meshes) will fail
    if (.not. verify_consistent_local_element_numbering(position1%mesh)) then
      ewrite(-1,*) "The local element ordering is not the same between processes"
      ewrite(-1,*) "that see the same element. This is a necessary condition on the"
      ewrite(-1,*) "decomposed input meshes for fluidity. The fact that you've"
      ewrite(-1,*) "obtained such meshes is likely a bug in fldecomp or the"
      ewrite(-1,*) "checkpointing code. Please report to the fluidity mailing"
      ewrite(-1,*) "list and state exactly how you've obtained your input files."
      FLAbort("Inconsistent local element ordering")
    end if
  end if
  mesh1=position1%mesh

  ! Insert mesh and position field into state and
  call insert(states(1), mesh1, "CoordinateMesh")
  call insert(states(1), mesh1, "TargetCoordinateMesh")
  call insert(states(1), position1, "Coordinate")
  call insert(states(1), position1, "TargetCoordinate")
  ! Setting path of mesh:
  ewrite(1,*) "ele_count(position1) = ", ele_count(position1)
  call deallocate(position1)
  ewrite(1,*) "Mesh file "//trim(input_basename_1)//" was successfully read in"


  ! Now the second mesh file:
  ewrite(1,*) "Reading in mesh file: "//trim(input_basename_2)
  ewrite(1,*) "input_mesh_format: "//trim(input_mesh_format)

  ! Read in the mesh file:
  position2=read_mesh_files(trim(input_basename_2), &
                      quad_degree=quad_degree, &
                      format=input_mesh_format)
  ! If reading in a decomposed mesh, read in halos as well:
  if(isparallel()) then
    call read_halos(trim(input_basename_2), position2)
    ! Local element ordering needs to be consistent between processes, otherwise
    ! code in Halos_Repair (used in halo construction of derived meshes) will fail
    if (.not. verify_consistent_local_element_numbering(position2%mesh)) then
      ewrite(-1,*) "The local element ordering is not the same between processes"
      ewrite(-1,*) "that see the same element. This is a necessary condition on the"
      ewrite(-1,*) "decomposed input meshes for fluidity. The fact that you've"
      ewrite(-1,*) "obtained such meshes is likely a bug in fldecomp or the"
      ewrite(-1,*) "checkpointing code. Please report to the fluidity mailing"
      ewrite(-1,*) "list and state exactly how you've obtained your input files."
      FLAbort("Inconsistent local element ordering")
    end if
  end if
  mesh2=position2%mesh

  ! Insert mesh and position field into state and
  call insert(states(2), mesh2, "CoordinateMesh")
  call insert(states(2), mesh2, "SourceCoordinateMesh")
  call insert(states(2), position2, "Coordinate")
  call insert(states(2), position2, "SourceCoordinate")
  ewrite(1,*) "ele_count(position2) = ", ele_count(position2)
  call deallocate(position2)
  ewrite(1,*) "Mesh file "//trim(input_basename_2)//" was successfully read in"


  ! Registering diagnostics:
  call register_diagnostic(dim=1, name='Alpha_Source_Integral', statistic='Value')
  call register_diagnostic(dim=1, name='Alpha_Target_Integral', statistic='Value')
  call initialise_diagnostics("diagnostics",states)


  ! Now that we read in the meshes, let's start processing them
  ! Get positions:
  pos_source => extract_vector_field(states(2), "SourceCoordinate")
  pos_target => extract_vector_field(states(1), "TargetCoordinate")
  assert(pos_source%dim == pos_target%dim)

  ! Insert scalar fields into each state:
  call allocate(alpha_source, pos_source%mesh, "AlphaSource")
  call zero(alpha_source)
  call insert(states(2), alpha_source, alpha_source%name)
  call allocate(alpha_target, pos_target%mesh, "AlphaTarget")
  call zero(alpha_target)
  call insert(states(1), alpha_target, alpha_target%name)
  ! Initializing alpha on the source mesh (to be 1):
  call set(alpha_source, 1.0)

  ! Do the projection:
  matching_domain = .false.
  if (ifinder == 0) then
    advfront = .false.
  else if (ifinder >= 0 .or. ifinder <= 4) then
    advfront = .true.
    if (ifinder == 1) then
      matching_domain = .true.
    end if
  else
    ewrite(-1,*) "Error, unsupported flag for intersection finder was found, choose a number between 0 and 4!"
    FLExit("Error due to unsupport intersection finder flag. Exit...")
  end if

  call nonmatching_domain_galerkin_projection(pos_target, pos_target, pos_source, alpha_target, &
                        & advfrontfinder=advfront, matching_domain=matching_domain)

  ewrite_minmax(alpha_target)
  
  ! RENDEZVOUS TESTING:
  ewrite(1,*) "=================================================================================="
  ewrite(1,*) "This is rank "//rank_str(1:rank_digits)//"."
  ! Let's have find us a rendezvous:
  ! First, computing bounding boxes of each mesh-partition
  dim = pos_source%dim
  allocate(my_bbox_src(dim, 2))
  allocate(my_bbox_trg(dim, 2))
  my_bbox_src = compute_bbox(pos_source)
  my_bbox_trg = compute_bbox(pos_target)
!  ewrite(1,*) "++++++++++++++++++++++++"
!  ewrite(1,*) "SOURCE MESH"
!  do i=1,pos_source%dim
!    ewrite(1,*) "dim = ", i
!    ewrite(1,*) "my_bbox_src(d,min) = ", my_bbox_src(i,1)
!    ewrite(1,*) "my_bbox_src(d,max) = ", my_bbox_src(i,2)
!  end do
!  ewrite(1,*) "++++++++++++++++++++++++"
!  ewrite(1,*) "TARGET MESH"
!  do i=1,pos_target%dim
!    ewrite(1,*) "dim = ", i
!    ewrite(1,*) "my_bbox_trg(d,min) = ", my_bbox_trg(i,1)
!    ewrite(1,*) "my_bbox_trg(d,max) = ", my_bbox_trg(i,2)
!  end do
  ! Now computing bboxes of other mesh partitions:
  ewrite(1,*) "nprocs = ", nprocs
  allocate(bboxes_src(dim, 2, nprocs))
  allocate(bboxes_trg(dim, 2, nprocs))
  ! Compute bboxes:
  call compute_bboxes(nprocs, bboxes_src, pos_source)
  call compute_bboxes(nprocs, bboxes_trg, pos_target)
  ! TESTING:
  ewrite(1,*) "==========================================================="
  ewrite(1,*) "checking bboxes_src:"
  ewrite(1,*) "rank = ", rank
  do i=1,nprocs
    ewrite(1,*) "----------------------------"
    ewrite(1,*) "nproc = ", i
    do j=1,dim
      ewrite(1,*) "+++++++++"
      ewrite(1,*) "d = ", j
      ewrite(1,*) "bboxes(d,min,rank) = ", bboxes_src(j,1,i)
      ewrite(1,*) "bboxes(d,max,rank) = ", bboxes_src(j,2,i)
    end do
  end do
  


  ! Compute integrals of alpha on both meshes:
  alpha_source_int = field_integral(alpha_source, pos_source)
  alpha_target_int = field_integral(alpha_target, pos_target)
  call set_diagnostic(name='Alpha_Source_Integral', statistic='Value', value=(/ alpha_source_int /))
  call set_diagnostic(name='Alpha_Target_Integral', statistic='Value', value=(/ alpha_target_int /))
  
  ! Dump the field on the target mesh to vtu:
  call vtk_write_fields(trim(input_basename_1)//"_Target_VolumeFraction", position=pos_target, model=pos_target%mesh, &
    & sfields=(/alpha_target/))
  ! Same for the source field/mesh:
  call vtk_write_fields(trim(input_basename_2)//"_Source_VolumeFraction", position=pos_source, model=pos_source%mesh, &
    & sfields=(/alpha_source/))

  ! Write diagnostics:
  call calculate_diagnostic_variables(states)
  call calculate_diagnostic_variables_new(states)
  call write_diagnostics(states, 0.0, 0.1, 1, not_to_move_det_yet=.true.)

  ! We are done here, deallocating stuff:
  call deallocate(states)
  deallocate(bboxes_src)
  deallocate(bboxes_trg)
  deallocate(my_bbox_src)
  deallocate(my_bbox_trg)
  

  ewrite(1, *) "Exiting Nonmatching_domain_interpolation"

  contains


    function compute_bbox(positions) result(bbox)
      ! Computes the bbox of an entire mesh/partition
      type(vector_field), intent(in) :: positions
      real, dimension(positions%dim, 2) :: bbox
      integer :: i

      do i=1,positions%dim
        bbox(i, 1) = minval(positions%val(i,:))
        bbox(i, 2) = maxval(positions%val(i,:))
      end do
    end function compute_bbox


    subroutine compute_bboxes(nprocs, bboxes, positions)
      integer, intent(in) :: nprocs
      real, dimension(:, :, :), intent(inout) :: bboxes
      type(vector_field), intent(in) :: positions
      real, dimension(size(bboxes, 2), size(bboxes, 3)) :: my_bboxes
      integer :: rank
      integer :: i
      integer :: ierr

      ewrite(1,*) "Inside compute_bboxes"

      my_bboxes(:, :) = compute_bbox(positions)
      ! Communicate all bboxes to rank 0:
      call mpi_gather(my_bboxes, size(my_bboxes), getpreal(), &
                      bboxes, size(my_bboxes), getpreal(), &
                      0, MPI_COMM_FEMTOOLS, ierr)
      assert(ierr == 0)

      ! Broadcast assembles array of bboxes to everyone.
      call mpi_bcast(bboxes, size(bboxes), getpreal(), 0, MPI_COMM_FEMTOOLS, ierr)
      assert(ierr == 0)

      ewrite(1,*) "Leaving compute_bboxes"

    end subroutine compute_bboxes


    subroutine nonmatching_domain_galerkin_projection(fieldA, positionsA, positionsB, alpha_AB, advfrontfinder, matching_domain)
      !! Return the volume fraction scalar field by projecting unity from the supermesh to the fluid mesh
      !! and, iff given, the solid velocity from the solid (via the supermesh) to the fluid mesh.
      !! Since positionsF and positionsS are different, we need to supermesh!
      !! positionsF and positionsS are the coordinate fields of the fluid and solid mesh respectively

      ! A stands for target, B for source, C for the supermesh in this case:
      type(vector_field), intent(inout) :: fieldA
      type(vector_field), intent(inout) :: positionsA, positionsB
      type(scalar_field), intent(inout) :: alpha_AB
      logical, intent(in), optional :: advfrontfinder
      logical, intent(in), optional :: matching_domain

      ! this is for using the advancing front algorithm, therefore commented for now:
      type(ilist), dimension(ele_count(positionsB)) :: map_AB 
      integer :: ele_A, ele_B

      type(quadrature_type) :: supermesh_quad
      type(element_type) :: supermesh_positions_shape, supermesh_field_shape

      type(vector_field) :: supermesh
      type(mesh_type) :: supermesh_field_mesh

      ! Scalar fields for alpha:
      type(scalar_field) :: alpha_AB_on_B
      type(scalar_field) :: alpha_AB_on_supermesh

      type(csr_sparsity) :: sparsity_A
      type(csr_matrix) :: mass_matrix_A
      type(scalar_field) :: rhs_alpha_AB
      type(vector_field) :: rhs_A

      real, dimension(ele_loc(positionsB, 1), ele_loc(positionsB, 1)) :: inversion_matrix_B
      real, dimension(ele_loc(positionsB, 1), ele_loc(positionsB, 1), ele_count(positionsA)) :: inversion_matrices_A
      integer :: dim, max_degree

      ! Variables for ilist when using rtree intersection finder:
      type(ilist) :: map_AB_rtree
      integer :: nele_AB, j, maplen, ntests

      character(len=OPTION_PATH_LEN) :: tmp

      integer :: stat, sstat, ierr
      
      ! Which intersection finder algorithm do we want to use?
      logical :: advfront

      ! As the projections are bounded, we need to define lumped versions of the mass matrices:
      type(scalar_field) :: lumped_mass_matrix_A
      type(scalar_field) :: lumped_inverse_mass_matrix_A

      if (.not. present_and_true(advfrontfinder)) then
          advfront = .false.
      else
          advfront = advfrontfinder
      end if

      ewrite(2,*) "inside nonmatching_domain_galerkin_projection"
      ewrite(1,*) "======================================================="
      ewrite(1,*) "advfront = ", advfront
      ewrite(1,*) "======================================================="

      dim = mesh_dim(positionsB)
      call intersector_set_dimension(dim)

      ! Define polynomial degree and quadrature of the supermesh:
      max_degree = max(element_degree(fieldA, 1), element_degree(positionsB, 1))
      supermesh_quad = make_quadrature(vertices=ele_loc(positionsB, 1), dim=dim, degree=2*max_degree)
      ! Define the shape function of the supermesh:
      supermesh_positions_shape = make_element_shape(vertices=ele_loc(positionsB, 1), dim=dim, degree=1, quad=supermesh_quad)
      supermesh_field_shape = make_element_shape(vertices=ele_loc(positionsB, 1), dim=dim, degree=max_degree, quad=supermesh_quad)
      ! For each element in B get the list of intersecting elements in A
      ! Replace advancing_front_intersection_finder for now
      ! (more work needs to be done to make this work in parallel for fsi-modelling)
      ! with the more robust rtree_intersection_finder:
      ! advancing_front:
      if (advfront) then
          ! advancing_front:
          map_AB = intersection_finder(positionsB, positionsA, matching_domain=.true.)
      else
          ! rtree:
          call rtree_intersection_finder_set_input(positionsA)
      end if

      sparsity_A = make_sparsity(fieldA%mesh, fieldA%mesh, "AMassMatrixSparsity")
      call allocate(mass_matrix_A, sparsity_A, name="AMassMatrix")
      call zero(mass_matrix_A)
      call deallocate(sparsity_A)

      ! Allocate lumped versions of the mass matrices:
      call allocate(lumped_mass_matrix_A, fieldA%mesh, "ALumpedMassMatrix")
      call zero(lumped_mass_matrix_A)
      call compute_lumped_mass(positionsA, lumped_mass_matrix_A)
      call allocate(lumped_inverse_mass_matrix_A, fieldA%mesh, "ALumpedInverseMass")
      call invert(lumped_mass_matrix_A, lumped_inverse_mass_matrix_A)

      ! get the matrix for the coordinates of (fluid) mesh A:
      do ele_A=1,ele_count(positionsA)
        call local_coords_matrix(positionsA, ele_A, inversion_matrices_A(:, :, ele_A))
        call assemble_mass_matrix(mass_matrix_A, fieldA, positionsA, ele_A)
      end do

      ! Fields for Alpha:
      call allocate(alpha_AB_on_B, positionsB%mesh, "AlphaBMesh")
      ! alpha_AB on the B mesh is unity by definition:
      call set(alpha_AB_on_B, 1.0)
      call allocate(rhs_alpha_AB, fieldA%mesh, "AlphaABOnARHS")
      call zero(rhs_alpha_AB)

      ! Looping over all the elements of mesh B:
      do ele_B=1,ele_count(positionsB)

        if (.not. advfront) then
            ! =====================================================================================================
            ! The following is obviously only required when using the rtree intersection finder:
            ! This will be replaced by code using the advancing front algorithm, once this has been modified to
            ! work in parallel for fluid-solid interaction modelling
            ! Via RTREE, find intersection of solid element 'ele_B' with the
            ! input mesh (=meshA 'positionsA'):
            call rtree_intersection_finder_find(positionsB, ele_B)
            ! Fetch output, the number of intersections of this solid element
            ! with the fluid mesh (= positionsA):
            call rtree_intersection_finder_query_output(nele_AB)
            ! Generate an ilist of elements in positionsA that intersect with ele_B:
            do j=1, nele_AB
              ! Get the donor (fluid) element which intersects with ele_B
              call rtree_intersection_finder_get_output(ele_A, j)
              ! insert_ascending works, but maybe it's not necessary
              !call insert_ascending(map_AB_rtree, ele_A)
              call insert(map_AB_rtree, ele_A)
            end do
            ! =====================================================================================================
        end if

        ! get the matrix for the coordinates of mesh B:
        call local_coords_matrix(positionsB, ele_B, inversion_matrix_B)
        if (advfront) then
            ! Construct the supermesh associated with ele_B. (For advancing front algorithm, uncomment the following line:)
            call construct_supermesh(positionsB, ele_B, positionsA, map_AB(ele_B), supermesh_positions_shape, supermesh, stat=stat)
            ! if no intersection for proc x was found, no supermesh was created, then stat/=0
            if (stat /= 0) then ! should not happen, as we catch that event above, but doesn't hurt double checking
              cycle
            end if
        else 
            ! =====================================================================================================
            ! When using the rtree intersection finder:
            if (map_AB_rtree%length > 0) then
              stat = 0
              call construct_supermesh(positionsB, ele_B, positionsA, map_AB_rtree, &
                         & supermesh_positions_shape, supermesh, stat=stat)
              ! if no intersection for proc x was found, no supermesh was created, then stat/=0
              if (stat /= 0) then
                cycle
              end if
            else
              cycle
            end if
            ! =====================================================================================================
        end if ! end if for which intersection finder to use

        ! At this point, a portion of the supermesh is constructed, which refines the dimensional space
        ! of the intersection of ele_B with meshA.
        ! This portion of the supermesh is stored in 'supermesh'

        ! 1st step: Compute the project unity from the solid mesh to the supermesh:

        ! Allocate field for alpha on the supermesh:
        supermesh_field_mesh = make_mesh(supermesh%mesh, supermesh_field_shape, -1, "SupermeshFieldMesh")
        call allocate(alpha_AB_on_supermesh, supermesh_field_mesh, "AlphaABOnSupermesh")
        call zero(alpha_AB_on_supermesh)

        ! Get alpha for this portion of the supermesh:
        call compute_alpha_on_supermesh(supermesh_field_shape, positionsB, ele_B, ele_val(alpha_AB_on_B, ele_B), &
                                        supermesh, inversion_matrix_B, alpha_AB_on_supermesh)

        ! 2nd step: Compute the RHS of the Galerkin projection:
        ! Elemental Values of alpha of this part of the supermesh are now computed.
        ! Compute the RHS of 
        ! M_f * alpha_f = M_{f,sup} * alpha_{sup}
        ! and afterwards, the supermesh can be deleted as it is no longer needed,
        ! remember, (M_{f,sup} * alpha_{sup}) lives on the fluid mesh, not the supermesh.
        call compute_rhs_galerkin_projection(ele_B, positionsA, positionsB, supermesh, inversion_matrices_A, &
                                                        rhs_alpha_AB, alpha_AB_on_supermesh)

        ! Portion of supermesh no longer needed
        call deallocate(supermesh)
        call deallocate(supermesh_field_mesh)
        call deallocate(alpha_AB_on_supermesh)
        if (.not. advfront) then
            ! Flush ilist of intersecting elements (when using the rtree intersection finder):
            !call flush_list(map_AB_rtree)
            call deallocate(map_AB_rtree)
        end if
      end do

      call deallocate(supermesh_quad)
      call deallocate(supermesh_positions_shape)
      call deallocate(supermesh_field_shape)

      if (advfront) then
          call deallocate(map_AB_rtree)
          ! =====================================================================================================
          ! the following is commented because of using rtree instead of the advancing front algorithm:
          do ele_B=1,ele_count(positionsB)
            call deallocate(map_AB(ele_B))
          end do
          ! =====================================================================================================
      else
          ! Because of using the rtree intersection finder:
          call rtree_intersection_finder_reset(ntests)
      end if

      ! 3rd step: Project alpha from the supermesh to the fluid and solid mesh:
      ! loop over fluid and solid elements and solve the last equation, which will
      ! project alpha from the supermesh to the fluid mesh

      ! Set solver options for the interpolations:
#ifdef HAVE_PETSC
     call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif
      alpha_AB%option_path = "/tmp/galerkin_projection/continuous"
      call add_option("/tmp/galerkin_projection/continuous", sstat)
      call add_option("tmp/galerkin_projection/continuous/bounded::Diffuse", sstat)
      call add_option("tmp/galerkin_projection/continuous/bounded::Diffuse/boundedness_iterations", sstat)
      call set_option("tmp/galerkin_projection/continuous/bounded::Diffuse/boundedness_iterations", 5000, sstat)
      call add_option("tmp/galerkin_projection/continuous/bounded::Diffuse/repair_deviations", sstat)
      call add_option("tmp/galerkin_projection/continuous/solver", sstat)
      call add_option("tmp/galerkin_projection/continuous/solver/iterative_method::cg", sstat)
      call add_option("tmp/galerkin_projection/continuous/solver/preconditioner::sor", sstat)
      call add_option("tmp/galerkin_projection/continuous/solver/relative_error", sstat)
      call set_option("tmp/galerkin_projection/continuous/solver/relative_error", 1e-9, sstat)
      call add_option("tmp/galerkin_projection/continuous/solver/max_iterations", sstat)
      call set_option("tmp/galerkin_projection/continuous/solver/max_iterations", 5000, sstat)
      call add_option("tmp/galerkin_projection/continuous/solver/never_ignore_solver_failures", sstat)
      !call add_option("tmp/galerkin_projection/continuous/solver/diagnostics", sstat)
      !call add_option("tmp/galerkin_projection/continuous/solver/diagnostics/monitors", sstat)

      ewrite(1,*) "Option tree for the projection:"
      call print_options()

      ! Project alpha to the target mesh:
      call petsc_solve(alpha_AB, mass_matrix_A, rhs_alpha_AB)

      ! Always bound the projection:
      call bound_projection_scalars(alpha_AB, rhs_alpha_AB, mass_matrix_A, &
                            lumped_mass_matrix_A, lumped_inverse_mass_matrix_A, &
                            positionsA)

      call deallocate(mass_matrix_A)
      call deallocate(alpha_AB_on_B)
      call deallocate(rhs_alpha_AB)
      call deallocate(lumped_mass_matrix_A)
      call deallocate(lumped_inverse_mass_matrix_A)

      ewrite(2,*) "leaving nonmatching_domain_galerkin_projection"

    end subroutine nonmatching_domain_galerkin_projection

    subroutine compute_alpha_on_supermesh(supermesh_field_shape, B_positions, ele_B, alpha_AB_ele_value, &
                                          supermesh, inversion_matrix_B, alpha_AB_on_supermesh)
      ! F stands for fluid, S for solid, C for the supermesh in this case:
      type(vector_field), intent(in) :: B_positions, supermesh
      type(element_type), intent(in), target :: supermesh_field_shape
      integer, intent(in) :: ele_B
      real, dimension(:), intent(in) :: alpha_AB_ele_value
      real, dimension(:, :), intent(in) :: inversion_matrix_B

      type(scalar_field), intent(inout) :: alpha_AB_on_supermesh

      type(mesh_type) :: supermesh_field_mesh

      integer :: ele_C
      integer, dimension(:), pointer :: nodes
      integer :: node_C, i, j, k

      integer :: ele_A
      real, dimension(ele_loc(B_positions, ele_B)) :: local_coords
      integer :: dim
      real :: val
      type(vector_field) :: supermesh_positions_remapped

      dim = B_positions%dim

      supermesh_field_mesh = make_mesh(supermesh%mesh, supermesh_field_shape, -1, "SupermeshFieldMesh")
      call allocate(supermesh_positions_remapped, dim, supermesh_field_mesh, "SupermeshPositionsRemapped")
      call remap_field(supermesh, supermesh_positions_remapped)

      ! Looping over all the elements of the portion of the supermesh:
      do ele_C=1,ele_count(supermesh)
        ! Get the global node numbers of element 'ele_C' of the supermesh:
        nodes => ele_nodes(alpha_AB_on_supermesh, ele_C)
        ! ele_A is then the region (number) of the supermesh, which describes the intersection of 
        ! one element of mesh A with one element of mesh B; the supermesh being the supermesh of all 
        ! the intersections of element 'ele_B' with mesh A.
        ! store the inverted matrix with the coordinates of mesh region 'ele_A' in inversion_matrix_A:

        ! Loop over the nodes of ele_C:
        do i=1,size(nodes)
          ! node_C = global node number of the i-th node of element ele_C:
          node_C = nodes(i)

          ! set local_coords to the values of node 'node_C' on the element 'ele_C' of the supermesh
          local_coords(1:dim) = node_val(supermesh_positions_remapped, node_C); local_coords(dim+1) = 1.0
          ! Compute the matrix multiplication of the coordinate matrix
          ! for region ele_F and local_coords, then store the result in local_coords:
          local_coords = matmul(inversion_matrix_B, local_coords)
          val = 0.0
          do j=1,supermesh_field_shape%loc
            val = val + eval_shape(supermesh_field_shape, j, local_coords) * alpha_AB_ele_value(j)
          end do
          call set(alpha_AB_on_supermesh, node_C, val)

        end do ! end of looping over the nodes of element 'ele_C' of the supermesh

      end do ! end of looping over elements of the supermesh

      ! supermesh_positions_remapped no longer needed:
      call deallocate(supermesh_positions_remapped)
      call deallocate(supermesh_field_mesh)

    end subroutine compute_alpha_on_supermesh

    subroutine compute_rhs_galerkin_projection(ele_B, positionsA, positionsB, supermesh, inversion_matrices_A, &
                                                           rhs_alpha_AB, alpha_AB_on_supermesh, &
                                                           rhs_A)

      ! FSI: fluid-solid interactions
      integer, intent(in) :: ele_B
      type(vector_field), intent(in) :: positionsA, positionsB
      type(vector_field), intent(in) :: supermesh
      real, dimension(:, :, :), intent(in) :: inversion_matrices_A
      type(scalar_field), intent(inout) :: rhs_alpha_AB
      type(scalar_field), intent(in) :: alpha_AB_on_supermesh
      type(vector_field), intent(inout), optional :: rhs_A

      integer :: ele_A, ele_C

      do ele_C=1,ele_count(supermesh)
        ele_A = ele_region_id(supermesh, ele_C) ! get the parent fluid element
        call compute_rhs_galerkin_projection_ele(positionsA, positionsB, supermesh, &
                                                            inversion_matrices_A(:, :, ele_A), ele_A, ele_B, ele_C, &
                                                            rhs_alpha_AB, alpha_AB_on_supermesh)
      end do
    end subroutine compute_rhs_galerkin_projection

    subroutine compute_rhs_galerkin_projection_ele(positionsA, positionsB, supermesh, &
                                                               inversion_matrix_A, ele_A, ele_B, ele_C, &
                                                               rhs_alpha_AB, alpha_AB_on_supermesh, &
                                                               rhs_A)

      type(vector_field), intent(in) :: positionsA, positionsB
      type(vector_field), intent(in) :: supermesh
      real, dimension(:, :), intent(in) :: inversion_matrix_A
      integer, intent(in) :: ele_A, ele_B, ele_C
      type(scalar_field), intent(inout) :: rhs_alpha_AB
      type(scalar_field), intent(in) :: alpha_AB_on_supermesh
      type(vector_field), intent(inout), optional :: rhs_A

      real, dimension(ele_ngi(supermesh, ele_C)) :: detwei_C
      real, dimension(positionsA%dim+1, ele_ngi(supermesh, ele_C)) :: local_coord_at_quad_A
      real, dimension(positionsB%dim+1, ele_ngi(supermesh, ele_C)) :: local_coord_at_quad_B
      real, dimension(supermesh%dim+1, ele_ngi(supermesh, ele_C)) :: global_coord_at_quad

      real, dimension(ele_loc(rhs_alpha_AB, ele_A), ele_ngi(supermesh, ele_C)) :: basis_at_quad_A
      real, dimension(ele_loc(rhs_alpha_AB, ele_C), ele_ngi(supermesh, ele_C)) :: basis_at_quad_supermesh

      real, dimension(ele_loc(rhs_alpha_AB, ele_A), ele_loc(alpha_AB_on_supermesh, ele_C)) :: mat_A

      real, dimension(ele_loc(supermesh, ele_C)) :: supermesh_alpha_AB_val

      integer :: j, k, l, dim

      ! Compute the local coordinates of the fluid and solid meshes.
      global_coord_at_quad(1:supermesh%dim, :) = ele_val_at_quad(supermesh, ele_C)
      global_coord_at_quad(supermesh%dim+1, :) = 1.0
      local_coord_at_quad_A = matmul(inversion_matrix_A, global_coord_at_quad)

      ! Compute the basis functions of the fluid, solid and supermesh force fields at these local coordinates.
      basis_at_quad_A = basis_at_quad(ele_shape(rhs_alpha_AB, ele_A), local_coord_at_quad_A)
      basis_at_quad_supermesh = basis_at_quad(ele_shape(alpha_AB_on_supermesh, ele_C), alpha_AB_on_supermesh%mesh%shape%n)

      ! We need to integrate over the supermesh element, so get detwei_C
      call transform_to_physical(supermesh, ele_C, detwei=detwei_C)

      ! Compute the little mixed mass matrices
      mat_A = 0.0
      !mat_solid = 0.0
      do j=1,ele_ngi(supermesh, ele_C)
        forall (k=1:ele_loc(rhs_alpha_AB, ele_A),l=1:ele_loc(alpha_AB_on_supermesh, ele_C))
          mat_A(k, l) = mat_A(k, l) + detwei_C(j) * basis_at_quad_A(k, j) * basis_at_quad_supermesh(l, j)
        end forall
      end do

      ! Now compute the rhs contribution for alpha (scalar field):
      supermesh_alpha_AB_val = ele_val(alpha_AB_on_supermesh, ele_C)
      call addto(rhs_alpha_AB, ele_nodes(rhs_alpha_AB, ele_A), matmul(mat_A, supermesh_alpha_AB_val))

    end subroutine compute_rhs_galerkin_projection_ele

    subroutine assemble_mass_matrix(mass, field, positions, ele)
      type(csr_matrix), intent(inout) :: mass
      type(vector_field), intent(in) :: field
      type(vector_field), intent(in) :: positions
      integer, intent(in) :: ele

      real, dimension(ele_ngi(field, ele)) :: detwei
      real, dimension(ele_loc(field, ele), ele_loc(field, ele)) :: little_mass

      call transform_to_physical(positions, ele, detwei=detwei)
      little_mass = shape_shape(ele_shape(field, ele), ele_shape(field, ele), detwei)
      call addto(mass, ele_nodes(field, ele), ele_nodes(field, ele), little_mass)
    end subroutine assemble_mass_matrix

    function basis_at_quad(shape, local_coords) result(basis)
      type(element_type), intent(in) :: shape
      real, dimension(:, :), intent(in) :: local_coords
      real, dimension(shape%loc, size(local_coords, 2)) :: basis
      integer :: loc, gi

      if (shape%degree == 0) then
        basis = 1.0
      elseif (shape%degree == 1) then
        basis = local_coords
      else
        do loc=1,shape%loc
          do gi=1,size(local_coords, 2)
            basis(loc, gi) = eval_shape(shape, loc, local_coords(:, gi))
          end do
        end do
      end if
    end function basis_at_quad

    subroutine bound_projection_scalars(sfield, rhs_field, mass_matrix, lumped_mass_matrix, inverse_mass_matrix_lumped, positions)
      ! Bounding a field (only scalar) after an inter-mesh projection.
      ! Similar to subroutine 'interpolation_galerkin_scalars' in module 'conservative_interpolation_module'
      type(scalar_field), intent(inout) :: sfield, lumped_mass_matrix
      type(scalar_field), intent(in) :: rhs_field, inverse_mass_matrix_lumped
      type(csr_matrix), intent(in) :: mass_matrix
      type(vector_field), intent(in) :: positions

      type(scalar_field) :: bounded_soln, max_bound, min_bound
      real :: upper_bound, lower_bound
      type(csr_sparsity), pointer :: nnlist
      integer :: node_B
      integer, dimension(:), pointer :: patch
      
      ewrite(2,*) "Inside bound_projection_scalars"

      ! Step 0: Computing the bounds,
      ! here we hard-code the default bounds:
      upper_bound = huge(0.0)*epsilon(0.0)
      lower_bound = -huge(0.0)*epsilon(0.0)

      ! Set the default bounds to the whole mesh:
      call allocate(max_bound, sfield%mesh, "MaxBound")
      call allocate(min_bound, sfield%mesh, "MinBound")
      call set(max_bound, upper_bound)
      call set(min_bound, lower_bound)
      call allocate(bounded_soln, sfield%mesh, "BoundedSolution")

      ! Preparing the solution:
      call set(bounded_soln, rhs_field)
      call scale(bounded_soln, inverse_mass_matrix_lumped)
      call halo_update(bounded_soln)

      nnlist => extract_nnlist(sfield)

      do node_B = 1,node_count(sfield%mesh)
        patch => row_m_ptr(nnlist, node_B)
        call set(max_bound, node_B, max(min(maxval(bounded_soln%val(patch)), &
                                                     node_val(max_bound, node_B)), &
                                                     lower_bound))
        call set(min_bound, node_B, max(min(minval(bounded_soln%val(patch)), &
                                                     node_val(max_bound, node_B)), &
                                                     lower_bound))
      end do

      call halo_update(max_bound)
      call halo_update(min_bound)

      call bound_field(sfield, max_bound, min_bound, mass_matrix, lumped_mass_matrix, &
                       inverse_mass_matrix_lumped, bounded_soln, positions)

      call deallocate(max_bound)
      call deallocate(min_bound)
      call deallocate(bounded_soln)

      ewrite(2,*) "Leaving bound_projection_scalars"

    end subroutine bound_projection_scalars

end subroutine Nonmatching_domain_interpolation
