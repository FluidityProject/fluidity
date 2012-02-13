#include "fdebug.h"
module fefields
  !!< Module containing general tools for discretising Finite Element problems.

  use fields
  use field_options, only: get_coordinate_field
  use elements, only: element_type
  use element_numbering
  use fetools, only: shape_shape
  use transform_elements, only: transform_to_physical, element_volume
  use sparse_tools
  use sparse_matrices_fields
  use state_module
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

contains
  
  subroutine compute_cv_mass(positions, cv_mass, porosity)
    
    !!< Compute the cv mass matrix associated with the 
    !!< input scalar fields mesh. This will use pre tabulated
    !!< coefficients to calculate each sub control volumes  
    !!< volume - which is only set up for constant, linear elements and 
    !!< selected quadratic elements. This assumes that all 
    !!< elements have the same vertices, degree and dim. Also 
    !!< the mesh element type must be Lagrangian. This WILL work
    !!< for both continuous and discontinuous meshes. If the element
    !!< order is zero then return the element volume.
    !!< If porosity is present this is included in the cv mass. This is 
    !!< only possible if the porosity is element wise (order zero) or 
    !!< the same order as the cv_mass. For the latter case the porosity 
    !!< is assumed represented by a control volume expansion and by 
    !!< a sub control volume expansion if discontinous 
    
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: cv_mass
    type(scalar_field), intent(in), optional :: porosity
    
    ! local variables
    integer :: ele
    integer :: vertices, polydegree, dim, type, family, loc
    real, dimension(:), pointer :: subcv_ele_volf => null()
    real, dimension(:), pointer :: subcv_ele_porosity => null()
    
    if (present(porosity)) then
       ewrite(1,*) 'In compute_cv_mass with porosity present'

       if (.not. (porosity%mesh%shape%degree == cv_mass%mesh%shape%degree) .and. &
           .not. (porosity%mesh%shape%degree == 0)) then
          FLAbort('The porosity must have the same order as cv_mass or be element wise (order zero)')
       end if 
    else
       ewrite(1,*) 'In compute_cv_mass'
    end if
    
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
            
    ! Form the CV mass matrix - include porosity if required
    if (present(porosity)) then
    
       allocate(subcv_ele_porosity(loc))

       do ele = 1,element_count(cv_mass)
          
          if (porosity%mesh%shape%degree == 0) then
             subcv_ele_porosity(1:1) = ele_val(porosity, ele)
             if (loc > 1) subcv_ele_porosity(2:) = subcv_ele_porosity(1)
          else
             subcv_ele_porosity = ele_val(porosity, ele)
          end if
          
          call addto(cv_mass, &
                     ele_nodes(cv_mass, ele), &
                     subcv_ele_volf * element_volume(positions, ele) * subcv_ele_porosity)
          
       end do

       deallocate(subcv_ele_porosity)
    
    else
    
       do ele = 1,element_count(cv_mass)
          
          call addto(cv_mass, &
                     ele_nodes(cv_mass, ele), &
                     subcv_ele_volf * element_volume(positions, ele))
          
       end do
    
    end if 
    
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
    type(mesh_type), intent(inout) :: mesh
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
           shape_rhs(to_shape, ele_val(from_field, ele)*detwei)))

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

end module fefields
