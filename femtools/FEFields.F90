#include "fdebug.h"
module fefields
  !!< Module containing general tools for discretising Finite Element problems.

  use fldebug
  use data_structures
  use mpi_interfaces
  use element_numbering
  use elements, only: element_type
  use parallel_tools
  use sparse_tools
  use parallel_fields
  use transform_elements, only: transform_to_physical, element_volume
  use fetools, only: shape_shape, shape_rhs, shape_vector_rhs
  use fields
  use state_module
  use field_options, only: get_coordinate_field
  use halos
  use sparse_matrices_fields
  implicit none

  interface add_source_to_rhs
    module procedure add_source_to_rhs_scalar, add_source_to_rhs_vector
  end interface add_source_to_rhs

  interface project_field
     module procedure project_scalar_field, project_vector_field
  end interface

  private
  public :: compute_lumped_mass, compute_mass, compute_projection_matrix, add_source_to_rhs, &
            compute_lumped_mass_on_submesh, compute_cv_mass, project_field
  public :: create_subdomain_mesh, create_parallel_redundant_mesh

contains

  subroutine compute_cv_mass(positions, cv_mass)

    !!< Compute the cv mass matrix associated with the
    !!< input scalar fields mesh. This will use pre tabulated
    !!< coefficients to calculate each sub control volumes
    !!< volume - which is only set up for constant, linear elements and
    !!< selected quadratic elements. This assumes that all
    !!< elements have the same vertices, degree and dim. Also
    !!< the mesh element type must be Lagrangian. This WILL work
    !!< for both continuous and discontinuous meshes. If the element
    !!< order is zero then return the element volume.

    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: cv_mass

    ! local variables
    integer :: ele
    integer :: vertices, polydegree, dim, type, family, loc
    real, dimension(:), pointer :: subcv_ele_volf => null()

    ewrite(1,*) 'In compute_cv_mass'

    ! sanity check
    assert(element_count(positions) == element_count(cv_mass))

    ! initialise
    call zero(cv_mass)

    ! get element info (assume all the same for whole mesh)
    vertices   = ele_vertices(cv_mass,1)
    polydegree = cv_mass%mesh%shape%degree
    dim        = cv_mass%mesh%shape%dim
    type       = cv_mass%mesh%shape%numbering%type
    family     = cv_mass%mesh%shape%numbering%family
    loc        = cv_mass%mesh%shape%loc

    ! The element type must be Lagrangian
    if (type /= ELEMENT_LAGRANGIAN) then
       FLAbort('Can only find the CV mass if the element type is Lagrangian')
    end if

    ! The polydegree must be < 3
    if (polydegree > 2) then
       FLAbort('Can only find the CV mass if the element polynomial degree is 2 or less')
    end if

    ! If the polydegree is 2 then the element family must be Simplex
    if ((polydegree == 2) .and. (.not. family == FAMILY_SIMPLEX)) then
       FLAbort('Can only find the CV mass for a mesh with a 2nd degree element if the element familiy is Simplex')
    end if

    ! Find the sub CV element volume fractions

    allocate(subcv_ele_volf(loc))

    if (polydegree == 0) then

       ! dummy value for element wise
       subcv_ele_volf = 1.0

    else if (polydegree == 1) then

       ! for linear poly the volume of each
       ! subcontrol volume is ele_vol / loc

       subcv_ele_volf = 1.0/real(loc)

    else if (polydegree == 2) then

       ! for quadratic poly we only consider Simplex family

       if (vertices == 2) then

          subcv_ele_volf(1) = 0.25 ! 1/2 * 1/2 = 1/4  Vertex CV
          subcv_ele_volf(2) = 0.5 ! (1 - 2 * 1/4) / 1 Centre CV
          subcv_ele_volf(3) = 0.25 ! 1/2 * 1/2 = 1/4  Vertex CV

       else if (vertices == 3) then

          subcv_ele_volf(1) = 8.3333333333333333e-02 ! 1/3 * 1/4 = 1/12   Vertex CV
          subcv_ele_volf(2) = 0.25                   ! (1 - 3 * 1/12) / 3 Edge CV
          subcv_ele_volf(3) = 8.3333333333333333e-02 ! 1/3 * 1/4 = 1/12   Vertex CV
          subcv_ele_volf(4) = 0.25                   ! (1 - 3 * 1/12) / 3 Edge CV
          subcv_ele_volf(5) = 0.25                   ! (1 - 3 * 1/12) / 3 Edge CV
          subcv_ele_volf(6) = 8.3333333333333333e-02 ! 1/3 * 1/4 = 1/12   Vertex CV

       else if ((vertices == 4) .and. (dim == 3)) then

          subcv_ele_volf(1)  = 0.03125                ! 1/8 * 1/4 = 1/32   Vertex CV
          subcv_ele_volf(2)  = 1.4583333333333333e-01 ! (1 - 4 * 1/32) / 6 Edge CV
          subcv_ele_volf(3)  = 0.03125                ! 1/8 * 1/4 = 1/32   Vertex CV
          subcv_ele_volf(4)  = 1.4583333333333333e-01 ! (1 - 4 * 1/32) / 6 Edge CV
          subcv_ele_volf(5)  = 1.4583333333333333e-01 ! (1 - 4 * 1/32) / 6 Edge CV
          subcv_ele_volf(6)  = 0.03125                ! 1/8 * 1/4 = 1/32   Vertex CV
          subcv_ele_volf(7)  = 1.4583333333333333e-01 ! (1 - 4 * 1/32) / 6 Edge CV
          subcv_ele_volf(8)  = 1.4583333333333333e-01 ! (1 - 4 * 1/32) / 6 Edge CV
          subcv_ele_volf(9)  = 1.4583333333333333e-01 ! (1 - 4 * 1/32) / 6 Edge CV
          subcv_ele_volf(10) = 0.03125                ! 1/8 * 1/4 = 1/32   Vertex CV

       else

          FLAbort('No code to form the sub control volume element volume fractions if not a Simplex')

       end if

    else

       FLAbort('No code to form the sub control volume element volume fractions if poly degree is > 2')

    end if

    ! Form the CV mass matrix:
    do ele = 1,element_count(cv_mass)

       call addto(cv_mass, &
            ele_nodes(cv_mass, ele), &
            subcv_ele_volf * element_volume(positions, ele))

    end do

    deallocate(subcv_ele_volf)

  end subroutine compute_cv_mass

  subroutine compute_lumped_mass(positions, lumped_mass, density, vfrac)
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: lumped_mass
    type(scalar_field), intent(inout), target, optional :: density
    type(scalar_field), intent(inout), target, optional :: vfrac ! PhaseVolumeFraction field

    integer :: ele
    real, dimension(ele_ngi(lumped_mass, 1)) :: detwei
    type(element_type), pointer :: t_shape
    real, dimension(ele_loc(lumped_mass, 1), ele_loc(lumped_mass, 1)) :: mass_matrix
    type(scalar_field), pointer :: l_density, l_vfrac

    real, dimension(ele_ngi(lumped_mass, 1)) :: density_gi, vfrac_gi

    ewrite(1,*) 'In compute_lumped_mass'

    if(present(density)) then
      l_density => density
    else
      allocate(l_density)
      call allocate(l_density, lumped_mass%mesh, name="LocalDensity", field_type=FIELD_TYPE_CONSTANT)
      call set(l_density, 1.0)
    end if

    if(present(vfrac)) then
      l_vfrac => vfrac
    else
      allocate(l_vfrac)
      call allocate(l_vfrac, lumped_mass%mesh, name="LocalPhaseVolumeFraction", field_type=FIELD_TYPE_CONSTANT)
      call set(l_vfrac, 1.0)
    end if

    call zero(lumped_mass)

    do ele=1,ele_count(lumped_mass)
      t_shape => ele_shape(lumped_mass, ele)
      density_gi = ele_val_at_quad(l_density, ele)
      vfrac_gi = ele_val_at_quad(l_vfrac, ele)
      call transform_to_physical(positions, ele, detwei=detwei)
      mass_matrix = shape_shape(t_shape, t_shape, detwei*density_gi*vfrac_gi)
      call addto(lumped_mass, ele_nodes(lumped_mass, ele), sum(mass_matrix, 2))
    end do

    if(.not.present(density)) then
      call deallocate(l_density)
      deallocate(l_density)
    end if

    if(.not.present(vfrac)) then
      call deallocate(l_vfrac)
      deallocate(l_vfrac)
    end if

  end subroutine compute_lumped_mass

  subroutine compute_mass(positions, mesh, mass, lumped_mass, density)
    type(vector_field), intent(in) :: positions
    type(mesh_type), intent(in) :: mesh
    type(csr_matrix), intent(inout) :: mass
    type(scalar_field), intent(inout), optional :: lumped_mass
    type(scalar_field), intent(inout), target, optional :: density

    integer :: ele
    real, dimension(ele_ngi(mesh, 1)) :: detwei
    type(element_type), pointer :: t_shape
    real, dimension(ele_loc(mesh, 1), ele_loc(mesh, 1)) :: mass_matrix
    type(scalar_field), pointer :: l_density

    real, dimension(ele_ngi(mesh, 1)) :: density_gi

    ewrite(1,*) 'In compute_mass'

    if(present(density)) then
      l_density => density
    else
      allocate(l_density)
      call allocate(l_density, mesh, name="LocalDensity", field_type=FIELD_TYPE_CONSTANT)
      call set(l_density, 1.0)
    end if

    call zero(mass)
    if(present(lumped_mass)) then
      assert(lumped_mass%mesh==mesh)
      call zero(lumped_mass)
    end if

    do ele=1,ele_count(mesh)
      t_shape => ele_shape(mesh, ele)
      density_gi = ele_val_at_quad(l_density, ele)
      call transform_to_physical(positions, ele, detwei=detwei)
      mass_matrix = shape_shape(t_shape, t_shape, detwei*density_gi)
      call addto(mass, ele_nodes(mesh, ele), ele_nodes(mesh, ele), mass_matrix)
      if(present(lumped_mass)) then
        call addto(lumped_mass, ele_nodes(lumped_mass, ele), sum(mass_matrix, 2))
      end if
    end do

    if(.not.present(density)) then
      call deallocate(l_density)
      deallocate(l_density)
    end if

  end subroutine compute_mass

  subroutine compute_lumped_mass_on_submesh(state, lumped_mass, density, vfrac)

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: lumped_mass
    type(scalar_field), intent(in), optional :: density
    type(scalar_field), intent(in), optional :: vfrac

    type(mesh_type) :: submesh, unperiodic_submesh    ! submesh
    type(scalar_field) :: rho_pmesh, rho_submesh      ! density on parent and submesh
    type(vector_field) :: x_pmesh, x_submesh          ! coordinates on parent and submesh
    type(scalar_field) :: vfrac_pmesh, vfrac_submesh  ! PhaseVolumeFraction on parent and submesh
    type(scalar_field) :: masslump_submesh

    ewrite(1,*) 'In compute_lumped_mass_on_submesh'

    submesh = make_submesh(lumped_mass%mesh, name="SubMesh")

    x_pmesh=get_coordinate_field(state, lumped_mass%mesh)

    if (mesh_periodic(submesh)) then
      ! you should never have periodic coordinates
      ! remake the mesh using the shape and continuity
      ! of the desired mesh (but not periodic)
      unperiodic_submesh=make_submesh(x_pmesh%mesh, name="UnPeriodicSubMesh")
      call allocate(x_submesh, mesh_dim(submesh), unperiodic_submesh, name="SubMeshCoordinate")
      call deallocate(unperiodic_submesh)
    else
      call allocate(x_submesh, mesh_dim(submesh), submesh, name="SubMeshCoordinate")
    end if
    call set_to_submesh(x_pmesh, x_submesh)

    call deallocate(x_pmesh)

    if(present(density)) then
      call allocate(rho_pmesh, lumped_mass%mesh, name="ParentMeshDensity")
      call remap_field(density, rho_pmesh)

      call allocate(rho_submesh, submesh, name="SubMeshDensity")
      call set_to_submesh(rho_pmesh, rho_submesh)

      call deallocate(rho_pmesh)
    else
      call allocate(rho_submesh, submesh, name="DummySubMeshDensity", field_type=FIELD_TYPE_CONSTANT)
      call set(rho_submesh, 1.0)
    end if

    ! This is only included in multiphase simulations
    if(present(vfrac)) then
      call allocate(vfrac_pmesh, lumped_mass%mesh, name="ParentMeshPhaseVolumeFraction")
      call remap_field(vfrac, vfrac_pmesh)

      call allocate(vfrac_submesh, submesh, name="SubMeshPhaseVolumeFraction")
      call set_to_submesh(vfrac_pmesh, vfrac_submesh)

      call deallocate(vfrac_pmesh)
    else
      call allocate(vfrac_submesh, submesh, name="DummySubMeshPhaseVolumeFraction", field_type=FIELD_TYPE_CONSTANT)
      call set(vfrac_submesh, 1.0)
    end if

    call allocate(masslump_submesh, submesh, "TemporarySubMeshLumpedMass")

    call compute_lumped_mass(x_submesh, masslump_submesh, rho_submesh, vfrac_submesh)

    call set_from_submesh(masslump_submesh, lumped_mass)

    call deallocate(submesh)
    call deallocate(x_submesh)
    call deallocate(rho_submesh)
    call deallocate(vfrac_submesh)
    call deallocate(masslump_submesh)

  end subroutine compute_lumped_mass_on_submesh

  function compute_projection_matrix(to_mesh, from_mesh, position) &
       result (P)
    !!< Calculate the projection matrix from from_mesh to to_mesh.
    !!< from_mesh, to_mesh and positions must have the same topology but
    !!< naturally may have different continuity or shape functions.
    !!<
    !!< The projection equation is:
    !!<
    !!<   M T  = P F
    !!<
    !!< Where T is the to field, F is the from field and M is the mass
    !!< matrix on the same mesh as T.
    !!<
    !!< This function constructs the matrix P.
    type (csr_matrix) :: P
    type(mesh_type), intent(in) :: to_mesh, from_mesh
    type(vector_field), intent(in) :: position

    ! We produce P using a dcsr matrix as this is easy to code. If it
    ! becomes performance critical we could work with sparsities.
    type(dynamic_csr_matrix) :: tmpP

    integer :: ele

    call allocate(tmpP, rows=node_count(to_mesh), &
         columns=node_count(from_mesh), name="tmpP")

    do ele=1, element_count(to_mesh)

       call projection_matrix_element(tmpP, ele, from_mesh, to_mesh,&
            & position)

    end do

    P=dcsr2csr(tmpP)
    call deallocate(tmpP)

  contains

    subroutine projection_matrix_element(tmpP, ele, from_mesh, to_mesh,&
         & position)
      !!< Calculate the contribution to tmpP from ele.
      type(dynamic_csr_matrix), intent(inout) :: tmpP
      integer, intent(in) :: ele
      type(mesh_type), intent(in) :: from_mesh, to_mesh
      type(vector_field), intent(in) :: position

      real, dimension(ele_ngi(to_mesh, ele)) :: detwei
      type(element_type), pointer :: to_shape, from_shape
      integer, dimension(:), pointer :: to_ele, from_ele

      to_shape=>ele_shape(to_mesh, ele)
      from_shape=>ele_shape(from_mesh, ele)
      to_ele=>ele_nodes(to_mesh, ele)
      from_ele=>ele_nodes(from_mesh, ele)

      ! Work out change of variables for quadrature.
      call transform_to_physical(position, ele, detwei=detwei)

      call addto(tmpP, to_ele, from_ele, &
           &           shape_shape(to_shape, from_shape, detwei))

    end subroutine projection_matrix_element

  end function compute_projection_matrix

  subroutine project_scalar_field(from_field, to_field, X)
    !!< Project from_field onto to_field. If to_field is discontinuous then
    !!< this will be calculated using the full mass matrix locally.
    !!< Otherwise, the mass will be lumped.
    type(scalar_field), intent(in) :: from_field
    type(scalar_field), intent(inout) :: to_field
    type(vector_field), intent(in) :: X

    type(scalar_field) ::  masslump
    integer :: ele
    type(csr_matrix)   :: P
    type(mesh_type)    :: dg_mesh, cg_mesh


    if(from_field%mesh==to_field%mesh) then

       call set(to_field, from_field)

    else

       if (to_field%mesh%continuity<0) then
          ! DG case

          do ele=1,element_count(to_field)
             call dg_projection_ele(ele, from_field, to_field, X)
          end do

       else
          ! DG to CG case

          cg_mesh=to_field%mesh

          call allocate(masslump, cg_mesh, "LumpedMass")

          call compute_lumped_mass(X, masslump)
          ! Invert lumped mass.
          masslump%val=1./masslump%val

          dg_mesh=from_field%mesh

          P=compute_projection_matrix(cg_mesh, dg_mesh , X)

          call zero(to_field)
          ! Perform projection.
          call mult(to_field, P, from_field)
          ! Apply inverted lumped mass to projected quantity.
          call scale(to_field, masslump)

          call deallocate(masslump)
          call deallocate(P)

       end if

    end if

  contains

    subroutine dg_projection_ele(ele, from_field, to_field, X)
      integer :: ele
      type(scalar_field), intent(in) :: from_field
      type(scalar_field), intent(inout) :: to_field
      type(vector_field), intent(in) :: X

      real, dimension(ele_loc(to_field,ele), ele_loc(to_field,ele)) :: mass
      real, dimension(ele_ngi(to_field,ele)) :: detwei
      type(element_type), pointer :: to_shape

      call transform_to_physical(X, ele, detwei)

      to_shape=>ele_shape(to_field, ele)

      mass=shape_shape(to_shape, to_shape, detwei)

      call invert(mass)

      call set(to_field, ele_nodes(to_field, ele), &
           matmul(mass, &
           shape_rhs(to_shape, ele_val_at_quad(from_field, ele)*detwei)))

    end subroutine dg_projection_ele

  end subroutine project_scalar_field

  subroutine project_vector_field(from_field, to_field, X)
    !!< Project from_field onto to_field. If to_field is discontinuous then
    !!< this will be calculated using the full mass matrix locally.
    !!< Otherwise, the mass will be lumped.
    type(vector_field), intent(in) :: from_field
    type(vector_field), intent(inout) :: to_field
    type(vector_field), intent(in) :: X

    type(scalar_field) ::  masslump, dg_scalar, cg_scalar
    integer :: ele
    type(csr_matrix)   :: P
    type(mesh_type)    :: dg_mesh, cg_mesh
    integer            :: j


    if(from_field%mesh==to_field%mesh) then

       call set(to_field, from_field)

    else

       if (to_field%mesh%continuity<0) then
          ! DG case

          do ele=1,element_count(to_field)
             call dg_projection_ele(ele, from_field, to_field, X)
          end do

       else
          ! CG case

          cg_mesh=to_field%mesh

          call allocate(masslump, cg_mesh, "LumpedMass")

          call compute_lumped_mass(X, masslump)
          ! Invert lumped mass.
          masslump%val=1./masslump%val

          dg_mesh=from_field%mesh

          P=compute_projection_matrix(cg_mesh, dg_mesh, X)

          call zero(to_field)
          ! Perform projection.
          do j=1,to_field%dim
              cg_scalar=extract_scalar_field_from_vector_field(to_field, j)
              dg_scalar=extract_scalar_field_from_vector_field(from_field, j)
              call mult(cg_scalar, P, dg_scalar)
              call set(to_field, j, cg_scalar)
          end do

          ! Apply inverted lumped mass to projected quantity.
          call scale(to_field, masslump)

          call deallocate(masslump)
          call deallocate(P)

       end if

    end if

  contains

    subroutine dg_projection_ele(ele, from_field, to_field, X)
      integer :: ele
      type(vector_field), intent(in) :: from_field
      type(vector_field), intent(inout) :: to_field
      type(vector_field), intent(in) :: X

      real, dimension(ele_loc(to_field,ele), ele_loc(to_field,ele)) :: mass
      real, dimension(ele_ngi(to_field,ele)) :: detwei
      type(element_type), pointer :: to_shape

      integer :: dim

      call transform_to_physical(X, ele, detwei)

      to_shape=>ele_shape(to_field, ele)

      mass=shape_shape(to_shape, to_shape, detwei)

      call invert(mass)

      do dim=1,to_field%dim
         call set(to_field, dim, ele_nodes(to_field, ele), &
              matmul(mass, &
              shape_rhs(to_shape, ele_val_at_quad(from_field, dim, ele)*detwei)))
      end do

    end subroutine dg_projection_ele

    subroutine cg_projection_ele(ele, from_field, to_field, masslump, X)
      integer :: ele
      type(vector_field), intent(in) :: from_field
      type(vector_field), intent(inout) :: to_field
      type(scalar_field), intent(inout) :: masslump
      type(vector_field), intent(in) :: X

      real, dimension(ele_ngi(to_field,ele)) :: detwei
      type(element_type), pointer :: to_shape

      to_shape=>ele_shape(to_field, ele)

      call transform_to_physical(X, ele, detwei)

      call addto(masslump, ele_nodes(to_field, ele), &
           shape_rhs(to_shape, detwei))

      call addto(to_field, ele_nodes(to_field, ele), &
           shape_vector_rhs(to_shape, ele_val_at_quad(from_field, ele), detwei))

    end subroutine cg_projection_ele

  end subroutine project_vector_field

  subroutine add_source_to_rhs_scalar(rhs, source, positions)
    !!< Add in a source field to the rhs of a FE equation,
    !!< i.e. compute the integrals:
    !!<
    !!<  rhs_i=\int N_i source dV
    !!<
    !!< with source=\sum_j source_j M_j, this means we multiply
    !!< source with the mass matrix \int N_i M_j dV
    type(scalar_field), intent(inout):: rhs
    type(scalar_field), intent(in):: source
    !!< needed for integration:
    type(vector_field), intent(in):: positions

    real, dimension( ele_loc(rhs,1), ele_loc(source,1) ):: M
    real, dimension( ele_ngi(positions,1) ):: detwei
    integer, dimension(:), pointer:: nodes
    integer ele

    do ele=1, element_count(source)
      call transform_to_physical(positions, ele, detwei)
      M=shape_shape(ele_shape(rhs,ele), ele_shape(source, ele), detwei)
      nodes => ele_nodes(rhs, ele)
      call addto(rhs, nodes, matmul(M, ele_val(source, ele)))
    end do

  end subroutine add_source_to_rhs_scalar

  subroutine add_source_to_rhs_vector(rhs, source, positions)
    !!< Add in a source field to the rhs of a FE equation, i.e. compute the integrals:
    !!<
    !!<  rhs_i=\int N_i source dV
    !!<
    !!<
    !!< with source=\sum_j source_j M_j, this means we multiply
    !!< source with the mass matrix \int N_i M_j dV
    !!< This is the vector version, a direct copy of the scalar case
    type(vector_field), intent(inout):: rhs
    type(vector_field), intent(in):: source
    !!< needed for integration:
    type(vector_field), intent(in):: positions

    real, dimension( ele_loc(source,1), ele_loc(rhs,1) ):: M_transpose
    real, dimension( ele_ngi(positions,1) ):: detwei
    integer, dimension(:), pointer:: nodes
    integer ele

    do ele=1, element_count(source)
      call transform_to_physical(positions, ele, detwei)
      ! M_transpose_ji=\int M_j N_i ( source=\sum_j s_j N_j, testing with N_i )
      ! note this is the transpose of what we usually compute
      M_transpose=shape_shape(ele_shape(source, ele), ele_shape(rhs,ele), detwei)
      nodes => ele_nodes(rhs, ele)
      call addto(rhs, nodes, matmul(ele_val(source, ele), M_transpose))
    end do

  end subroutine add_source_to_rhs_vector

  subroutine create_subdomain_mesh(mesh, element_list, name, submesh, node_list)
    !!< Create a mesh that only covers part of the domain

    ! full mesh to take submesh from
    type(mesh_type), intent(in), target :: mesh
    ! elements that will make up the submesh
    integer, dimension(:), intent(in):: element_list
    ! name for the new submesh
    character(len=*), intent(in) :: name
    ! submesh created
    type(mesh_type), intent(out) :: submesh
    ! list of nodes in submesh (also functions as node map from submesh to full mesh)
    integer, dimension(:), pointer :: node_list

    ! integer set containing nodes in submesh:
    type(integer_set) :: submesh_node_set
    ! node mapping functions from full to submesh (=0 if not in submesh)
    integer, dimension(:), allocatable :: inverse_node_list

    ! Others:
    integer :: ele, ele_2, ni, edge_count, i, node, loc, sloc, face, surf_ele_count
    integer, dimension(:), pointer :: neigh, faces

    type(element_type), pointer :: shape

    type(integer_hash_table) :: face_ele_list

    integer, allocatable, dimension(:) :: sndglno, boundary_ids, element_owner

    ewrite(1,*) "Entering create_subdomain_mesh"

    ! Build element mapping functions to and from sub mesh:

    ewrite(1,*) 'Number of elements in submesh:', size(element_list)

    ! Derive node list for subdomain_mesh:
    call allocate(submesh_node_set)
    do i = 1, size(element_list)
       ele = element_list(i)
       call insert(submesh_node_set, ele_nodes(mesh, ele))
    end do

    allocate(node_list(key_count(submesh_node_set))) ! Nodal map from sub mesh --> full mesh
    node_list = set2vector(submesh_node_set)
    ewrite(1,*) 'Number of nodes in submesh:', size(node_list)

    allocate(inverse_node_list(node_count(mesh))) ! Nodal map from full mesh --> sub mesh
    ! if after it is set up, the value in inverse_subnode_list = 0, this means that that element of
    ! the full mesh does not have a corresponding element on the subdomain_mesh - i.e. it is not a part
    ! of the prognostic subdomain.
    inverse_node_list = 0
    do i = 1, key_count(submesh_node_set)
       node = node_list(i)
       inverse_node_list(node) = i
    end do

    ! Allocate subdomain_mesh:
    shape => mesh%shape
    call allocate(submesh, nodes=size(node_list), elements=size(element_list),&
         & shape=shape, name=trim(name))
    submesh%option_path = mesh%option_path
    if (associated(mesh%region_ids)) then
      allocate(submesh%region_ids(size(element_list)))
      submesh%region_ids = mesh%region_ids(element_list)
    end if

    ! Determine ndglno (connectivity matrix) on subdomain_mesh:
    loc = shape%loc
    do i = 1, size(element_list)
      ele = element_list(i)
      ! can't use set_ele_nodes as it would create circularity between fields_allocates and fields_manipulation modules
      submesh%ndglno(loc*(i-1)+1:loc*i) = inverse_node_list(ele_nodes(mesh, ele))
    end do

    ! Calculate sndglno - an array of nodes corresponding to edges along surface:
    sloc = mesh%faces%shape%loc
    surf_ele_count = surface_element_count(mesh)

    ! Begin by determining which faces are on submesh boundaries:
    call allocate(face_ele_list)
    do i = 1, size(element_list)
      ele = element_list(i)
      neigh => ele_neigh(mesh, ele) ! Determine element neighbours on parent mesh
      faces => ele_faces(mesh, ele) ! Determine element faces on parent mesh
      do ni = 1, size(neigh)
        ele_2 = neigh(ni)
        face = faces(ni)
        ! If this face is part of the full surface mesh (which includes internal faces) then
        ! it must be on the submesh boundary, and not on a processor boundary (if parallel).
        ! note that for internal facets that are on now on the boundary of the subdomain, we only
        ! collect one copy, whereas for internal facets that remain internal we collect both
        ! this is dealt with using the allow_duplicate_internal_facets flag to add_faces()
        if (face  <= surf_ele_count) then
           call insert(face_ele_list, face, ele)
        end if
      end do
    end do

    ! Set up sndglno and boundary_ids:
    edge_count = key_count(face_ele_list)
    allocate(sndglno(edge_count*sloc), boundary_ids(1:edge_count))
    do i = 1, edge_count
      call fetch_pair(face_ele_list, i, face, ele)
      sndglno((i-1)*sloc+1:i*sloc) = inverse_node_list(face_global_nodes(mesh, face))
      boundary_ids(i) = surface_element_id(mesh, face)
    end do

    ewrite(2,*) "Number of surface elements: ", edge_count
    ! Add faces to submesh:
    if (has_discontinuous_internal_boundaries(mesh)) then
      allocate(element_owner(1:edge_count))
      do i=1, edge_count
        call fetch_pair(face_ele_list, i, face, ele)
        element_owner(i) = ele
      end do
      call add_faces(submesh, sndgln=sndglno, boundary_ids=boundary_ids, &
        element_owner=element_owner)
      deallocate(element_owner)
    else
      call add_faces(submesh, sndgln=sndglno, boundary_ids=boundary_ids, &
        allow_duplicate_internal_facets=.true.)
    end if

    call deallocate(face_ele_list)
    deallocate(sndglno)
    deallocate(boundary_ids)

    ! If parallel then set up node and element halos, by checking whether mesh halos
    ! exist on submesh:

    if(isparallel()) then
       call generate_subdomain_halos(mesh, submesh, node_list, inverse_node_list)
    end if

    deallocate(inverse_node_list)
    call deallocate(submesh_node_set)

    ewrite(1,*) "Leaving create_subdomain_mesh"

  end subroutine create_subdomain_mesh

  subroutine generate_subdomain_halos(external_mesh,subdomain_mesh,node_list,inverse_node_list)

    type(mesh_type), intent(in) :: external_mesh
    type(mesh_type), intent(inout) :: subdomain_mesh
    integer, dimension(:) :: node_list, inverse_node_list

    integer :: nhalos, communicator, nprocs, procno, ihalo, nowned_nodes

    ewrite(1, *) "In generate_subdomain_halos"

    assert(continuity(subdomain_mesh) == 0)
    assert(.not. associated(subdomain_mesh%halos))
    assert(.not. associated(subdomain_mesh%element_halos))

    ! Initialise key MPI information:

    nhalos = halo_count(external_mesh)
    ewrite(2,*) "Number of subdomain_mesh halos = ",nhalos

    if(nhalos == 0) return

    communicator = halo_communicator(external_mesh%halos(nhalos))
    nprocs = getnprocs(communicator = communicator)
    ewrite(2,*) 'Number of processes = ', nprocs
    procno = getprocno(communicator = communicator)
    ewrite(2,*) 'Processor ID/number = ', procno

    ! Allocate subdomain mesh halos:
    allocate(subdomain_mesh%halos(nhalos))

    ! Derive subdomain_mesh halos:
    do ihalo = 1, nhalos

       subdomain_mesh%halos(ihalo) = derive_sub_halo(external_mesh%halos(ihalo),node_list)

       assert(trailing_receives_consistent(subdomain_mesh%halos(ihalo)))

       if(.not. serial_storage_halo(external_mesh%halos(ihalo))) then
          assert(halo_valid_for_communication(subdomain_mesh%halos(ihalo)))
          call create_global_to_universal_numbering(subdomain_mesh%halos(ihalo))
          call create_ownership(subdomain_mesh%halos(ihalo))
       end if

    end do ! ihalo

    if(all(serial_storage_halo(subdomain_mesh%halos))) then
      allocate(subdomain_mesh%element_halos(0))
    else
      allocate(subdomain_mesh%element_halos(nhalos))
      call derive_element_halo_from_node_halo(subdomain_mesh, &
        & ordering_scheme = HALO_ORDER_TRAILING_RECEIVES, create_caches = .true.)
    end if

  end subroutine generate_subdomain_halos

  function create_parallel_redundant_mesh(mesh) result (redundant_mesh)
    !!< Creates a global mesh that is redundant, i.e. each process sees the entire
    !!< mesh, constructed by gathering all local input meshes together. The global mesh
    !!< is equiped with a halo structure in which each process owns the nodes associated
    !!< with their input mesh. Since each input mesh is sent to each other process, to form
    !!< part of their copy of the entire global mesh, on a halo update a process will send
    !!< values for its owned nodes to all other processes.
    !!< Although the global redundant mesh is the same on all processes, the node and element numbering
    !!< is not. The first nodes and elements of the redundant mesh correspond to the node and element numbers
    !!< of the input mesh (so we have trailing receives). The input mesh is treated as a completely local
    !!< mesh, i.e. any global/halo structure it has is ignored, and each node and element is treated as unique
    !!< and not associated with any nodes/elements of input meshes on other processes.
    !!< The output redundant mesh only caries a single, nodal halo, in which all nodes associated with the input mesh
    !!< are owned (send) nodes, and all other nodes are receives. In this way a simple halo_update(field) suffices
    !!< to share all field values with all processes, where field is defined on the redundant_mesh. The send values of
    !!< field can first be set by copying over the values from another field defined on the local input mesh (using
    !!< the shared node numbering), or through a remap (using the shared element numbering).
    !!< It should be clear that setting up this redundant mesh structure, and communicating through halo updates on it
    !!< can be extremely expensive and will not in general scale very well - for small meshes such as a surface mesh
    !!< this might still be feasible.

    !! Local input mesh. Treated as independent on each process, ignoring any halos
    type(mesh_type), intent(inout):: mesh
    !! Output global redudant mesh
    type(mesh_type):: redundant_mesh

    type(halo_type), pointer :: halo
    integer, dimension(:), allocatable :: eles_per_proc, eles_displs
    integer, dimension(:), allocatable :: nodes_per_proc, nsends
    integer, dimension(:), allocatable :: nodes_to_send, nodes_to_recv
    integer, dimension(:, :), allocatable :: eles_send_buf, eles_recv_buf
    integer neles, nloc, nprocs, procno
    integer comm, ierr
    integer i, j, node_offset, ele_offset, send_count

    ewrite(1,*) "Entering create_parallel_redundant_mesh for input mesh ", trim(mesh%name)

    neles = element_count(mesh)
    nloc = ele_loc(mesh, 1)
    nprocs = getnprocs()
    procno = getprocno()
    comm = MPI_COMM_FEMTOOLS

    ! exchange element numbers with all other processes
    allocate(eles_per_proc(nprocs))
    call mpi_allgather(neles, 1, MPI_INTEGER, &
      eles_per_proc, 1, MPI_INTEGER, comm, ierr)
    assert(ierr == MPI_SUCCESS)

    ! our local mesh to be sent to everyone
    allocate(eles_send_buf(nloc, neles))
    do i=1, neles
      eles_send_buf(:,i) = ele_nodes(mesh, i)
    end do

    allocate(eles_recv_buf(nloc, sum(eles_per_proc)), eles_displs(nprocs+1))
    ! displacements where in eles_recv_buf to put each bit of the recv'd surface mesh
    ! (indexed from 0 - would be nice if they'd actually tell you these things in the mpi docs!)
    j = 0
    do i=1, nprocs
      eles_displs(i) = j
      j = j + eles_per_proc(i)*nloc
    end do
    ! extra entry nprocs+1 for convenience:
    eles_displs(i) = j
    assert(j == size(eles_recv_buf))

    call mpi_allgatherv(eles_send_buf, nloc*neles, MPI_INTEGER, &
      eles_recv_buf, nloc*eles_per_proc, eles_displs, MPI_INTEGER, comm, ierr)
    assert(ierr == MPI_SUCCESS)

    ! convert from 0-based index in flattened eles_recv_buf
    ! to 1-based index for 2nd dim of eles_recv_buf(:, :)
    eles_displs = eles_displs/nloc + 1

    ! nodes per processor - this will be used later as n/o recv nodes, so we don't count our own
    allocate(nodes_per_proc(1:nprocs))
    do i=1, nprocs
      if (i==procno) then
        nodes_per_proc(i) = 0
      else if (eles_displs(i+1)>eles_displs(i)) then
        nodes_per_proc(i) = maxval(eles_recv_buf(:, eles_displs(i):eles_displs(i+1)-1))
      else
        nodes_per_proc(i) = 0
      end if
    end do

    ! n/o nodes to send to each other process based on max node number in elements
    ! (typically the same as node_count(mesh), but we could have spurious isolated nodes)
    if (neles>0) then
      send_count = maxval(eles_send_buf)
    else
      send_count = 0
    end if

    call allocate(redundant_mesh, send_count + sum(nodes_per_proc), size(eles_recv_buf, 2), mesh%shape, &
      name="Verticallyredundant"//trim(mesh%name))
    do i=1, neles
      call set_ele_nodes(redundant_mesh, i, ele_nodes(mesh, i))
    end do

    node_offset = send_count
    ele_offset = neles

    do i=1, nprocs
      if (i==procno) cycle
      do j=eles_displs(i), eles_displs(i+1)-1
        ele_offset = ele_offset + 1
        call set_ele_nodes(redundant_mesh, ele_offset, eles_recv_buf(:, j)+node_offset)
      end do
      node_offset = node_offset + nodes_per_proc(i)
    end do
    assert(node_offset==node_count(redundant_mesh))
    assert(ele_offset==element_count(redundant_mesh))

    allocate(redundant_mesh%halos(1))
    halo => redundant_mesh%halos(1)
    allocate(nsends(1:nprocs))

    ! send all "owned" nodes (owned nodes in redudant_mesh which is based on max. node number in owned elements of mesh)
    nsends = send_count
    ! but we don't send it (and recv) these to (from) ourself
    nodes_per_proc(procno) = 0
    nsends(procno) = 0

    call allocate(halo, nsends, nodes_per_proc, nowned_nodes=send_count, &
      name="ParallelRedundant"//trim(mesh%name)//"Halo", &
      communicator=comm, data_type=HALO_TYPE_CG_NODE)

    allocate(nodes_to_send(1:send_count))
    allocate(nodes_to_recv(1:node_count(redundant_mesh)-send_count))
    nodes_to_send = (/ (i, i=1, send_count) /)
    nodes_to_recv = (/ (i, i=send_count+1, node_count(redundant_mesh)) /)
    node_offset = 1
    do i=1, nprocs
      if (i==procno) cycle
      call set_halo_sends(halo, i, nodes_to_send)
      call set_halo_receives(halo, i, nodes_to_recv(node_offset:node_offset+nodes_per_proc(i)-1))
      node_offset = node_offset + nodes_per_proc(i)
    end do

    assert(trailing_receives_consistent(halo))
    assert(halo_valid_for_communication(halo))
    call create_ownership(halo)
    call create_global_to_universal_numbering(halo)

  end function create_parallel_redundant_mesh

end module fefields
