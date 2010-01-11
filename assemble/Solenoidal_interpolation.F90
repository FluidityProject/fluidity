#include "fdebug.h"
! #define CHECK_SUPERMESH_VOLUMES
! #define DUMP_SUPERMESH_INTERSECTIONS

module solenoidal_interpolation_module

  use fields
  use sparse_tools
  use supermesh_construction
  use futils
  use transform_elements
  use sparsity_patterns
  use vector_tools
  use tensors
  use fetools
  use sparse_tools
  use interpolation_module
  use solvers
  use spud
  use assemble_cmc
  use sparse_matrices_fields
  use boundary_conditions
  use boundary_conditions_from_options
  use momentum_cg, only: correct_masslumped_velocity, add_kmk_matrix, add_kmk_rhs, assemble_kmk_matrix
  use momentum_dg, only: correct_velocity_dg
  use fefields
  use divergence_matrix_cv, only: assemble_divergence_matrix_cv
  use divergence_matrix_cg, only: assemble_divergence_matrix_cg
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
  use dgtools
  implicit none

  private
  public :: solenoidal_interpolation
  
  interface solenoidal_interpolation
    module procedure solenoidal_interpolation_state, solenoidal_interpolation_fields
  end interface

  contains
  
  subroutine solenoidal_interpolation_state(state)
    type(state_type), intent(inout) :: state

    type(vector_field), pointer :: v_field
    type(vector_field), pointer :: coordinate
    type(scalar_field), pointer :: s_field
    type(mesh_type), pointer :: lagrange_mesh
    
    character(len=FIELD_NAME_LEN) :: mesh_name, update_field_name
    integer :: stat, i

    coordinate => extract_vector_field(state, "Coordinate")
    
    do i = 1, vector_field_count(state)
      v_field => extract_vector_field(state, i)
      if(trim(v_field%name)=="Coordinate") cycle
      
      if(have_option(trim(complete_field_path(v_field%option_path, stat))//&
                     "/enforce_discrete_properties/solenoidal/lagrange_multiplier/update_scalar_field")) then
        call get_option(trim(complete_field_path(v_field%option_path, stat))//&
                     "/enforce_discrete_properties/solenoidal/lagrange_multiplier/update_scalar_field/name", &
                     update_field_name)
        s_field=> extract_scalar_field(state, trim(update_field_name))
      else
        s_field=>null()
      end if
      
      call get_option(trim(complete_field_path(v_field%option_path, stat))//&
                    "/enforce_discrete_properties/solenoidal/lagrange_multiplier/mesh/name", &
                    mesh_name)
      lagrange_mesh=>extract_mesh(state, trim(mesh_name))
      if(associated(s_field)) then
        assert(s_field%mesh%name==mesh_name)
      end if
      
      call solenoidal_interpolation(v_field, coordinate, &
                                  & lagrange_mesh, s_field=s_field)
    
    end do

  end subroutine solenoidal_interpolation_state

  subroutine solenoidal_interpolation_fields(v_field, coordinate, &
                                             lagrange_mesh, s_field)
    
    type(vector_field), intent(inout) :: coordinate
    type(vector_field), intent(inout) :: v_field
    type(mesh_type), intent(inout) :: lagrange_mesh
    
    type(scalar_field), pointer :: s_field
    
    logical :: dg, lump_mass, lump_on_submesh, div_cv, div_cg, apply_kmk
    integer :: dim, j
    real :: dt
    
    type(block_csr_matrix), target :: ct_m
    type(block_csr_matrix), pointer :: ctp_m
    type(scalar_field) :: ct_rhs, kmk_rhs
    type(csr_sparsity) :: ct_m_sparsity, cmc_m_sparsity
    type(csr_matrix) :: cmc_m
    type(block_csr_matrix) :: inverse_field_mass
    type(scalar_field) :: field_lumped_mass
    type(vector_field) :: inverse_field_lumped_mass_vector
    
    ! This is the object of our desires:
    ! lagrange is the Lagrange multiplier that
    ! ensures solenoidality of the resulting interpolant
    type(scalar_field) :: lagrange

    ! The right hand side of this devilish equation:
    type(scalar_field) :: projec_rhs

    type(state_type) :: local_state

    character(len=OPTION_PATH_LEN) :: l_option_path
    
    call insert(local_state, coordinate, "Coordinate")
    
    l_option_path=trim(complete_field_path(v_field%option_path))//"/enforce_discrete_properties/solenoidal"
    
    dim = mesh_dim(v_field)
    
    if(associated(s_field)) then
      assert(trim(s_field%mesh%name) == trim(lagrange_mesh%name))
    end if

    dg = (continuity(v_field) < 0)
    
    lump_mass = have_option(trim(l_option_path)//&
                "/interpolated_field/discontinuous/lump_mass_matrix") &
                .or. .not.dg
    
    lump_on_submesh = have_option(trim(l_option_path)//&
                      "/interpolated_field/continuous/lump_mass_matrix/use_submesh") &
                      .and. .not.dg

    div_cv = have_option(trim(l_option_path) //&
          &"/lagrange_multiplier/spatial_discretisation/control_volumes")
    div_cg = .not.div_cv

    apply_kmk = (continuity(lagrange_mesh) >= 0 .and. lagrange_mesh%shape%degree == 1 &
          & .and. lagrange_mesh%shape%numbering%family == FAMILY_SIMPLEX .and. &
          & continuity(v_field) >= 0 .and. v_field%mesh%shape%degree == 1 &
          & .and. v_field%mesh%shape%numbering%family == FAMILY_SIMPLEX .and. &
          & .not. have_option(trim(l_option_path) // &
          & "/lagrange_multiplier/spatial_discretisation/continuous_galerkin/remove_stabilisation_term") .and. &
          & .not. div_cv)

    ct_m_sparsity = make_sparsity(lagrange_mesh, v_field%mesh, "DivergenceSparsity")
    call allocate(ct_m, ct_m_sparsity, blocks=(/1, dim/), name="DivergenceMatrix")
    call zero(ct_m)
    ctp_m => ct_m
    
    call allocate(ct_rhs, lagrange_mesh, "DivergenceRHS")
    call zero(ct_rhs)
    
    cmc_m_sparsity = make_sparsity_transpose(lagrange_mesh, v_field%mesh, "LagrangeProjectionSparsity")
    call allocate(cmc_m, cmc_m_sparsity, name="LagrangeProjectionMatrix")
    call zero(cmc_m)
    
    call allocate(projec_rhs, lagrange_mesh, "LagrangeProjectionRHS")
    call zero(projec_rhs)
    
    call allocate(lagrange, lagrange_mesh, "LagrangianMultiplier")
    call zero(lagrange)
    lagrange%option_path = trim(l_option_path)//"/lagrange_multiplier"

    if(apply_kmk) then
      call allocate(kmk_rhs, lagrange_mesh, "KMKRHS")
      call zero(kmk_rhs)
    end if

    if(lump_mass) then
      call allocate(field_lumped_mass, v_field%mesh, "FieldLumpedMass")
      call zero(field_lumped_mass)
      call allocate(inverse_field_lumped_mass_vector, dim, v_field%mesh, &
         "InverseFieldLumpedMassVector")
    else if (.not. dg) then
      FLAbort("Not possible to not lump the mass if not dg.")
    end if

    if(div_cg) then
      if(lump_mass) then
        if(lump_on_submesh) then
          call assemble_divergence_matrix_cg(ct_m, local_state, ct_rhs=ct_rhs, &
                                              test_mesh=lagrange_mesh, field=v_field, &
                                              option_path = trim(l_option_path)//"/lagrange_multiplier")
          ! now get the mass matrix lumped on the submesh
          call compute_lumped_mass_on_submesh(local_state, field_lumped_mass)

        else
          call assemble_divergence_matrix_cg(ct_m, local_state, ct_rhs=ct_rhs, &
                                              test_mesh=lagrange_mesh, field=v_field, &
                                              option_path = trim(l_option_path)//"/lagrange_multiplier", &
                                              grad_mass_lumped = field_lumped_mass)
        end if
      elseif(dg) then
        call assemble_divergence_matrix_cg(ct_m, local_state, ct_rhs=ct_rhs, &
                                            test_mesh=lagrange_mesh, field=v_field, &
                                            option_path = trim(l_option_path)//"/lagrange_multiplier")
        
        ! now get the dg inverse mass matrix
        call construct_inverse_mass_matrix_dg(inverse_field_mass, v_field, coordinate)
        
      else
        FLAbort("Not possible to not lump the mass if not dg.")
      end if
    elseif(div_cv) then
      if(lump_mass) then
        if(lump_on_submesh) then
          call compute_lumped_mass_on_submesh(local_state, field_lumped_mass)
        else
          call compute_lumped_mass(coordinate, field_lumped_mass)
        end if
      else if(dg) then
        call construct_inverse_mass_matrix_dg(inverse_field_mass, v_field, coordinate)
      else
        FLAbort("Not possible to not lump the mass if not dg.")
      end if
    
      call assemble_divergence_matrix_cv(ct_m, local_state, ct_rhs=ct_rhs, &
                                         test_mesh=lagrange_mesh, field=v_field)

    else
      FLAbort("Unknown spatial discretisation option for the lagrange multiplier.")
    end if
    
    if(lump_mass) then
      call invert(field_lumped_mass)
      
      do j=1, inverse_field_lumped_mass_vector%dim
        call set(inverse_field_lumped_mass_vector, j, field_lumped_mass)
      end do
        
      
      call apply_dirichlet_conditions_inverse_mass(inverse_field_lumped_mass_vector, v_field)
      
      call assemble_masslumped_cmc(cmc_m, ctp_m, inverse_field_lumped_mass_vector, ct_m)
    else if(dg) then
      call assemble_cmc_dg(cmc_m, ctp_m, ct_m, inverse_field_mass)
    else
      FLAbort("Not possible to not lump the mass if not dg.")
    end if
    
    call mult(projec_rhs, ct_m, v_field)
    
    if(apply_kmk) then
      ! a hack to make sure the appropriate meshes are available for the
      ! construction of the sparsity
      call insert(local_state, lagrange_mesh, name="PressureMesh")
      call insert(local_state, v_field%mesh, name="VelocityMesh")
      ! end of hack... you can look again now
      
      ewrite(2,*) "Assembling P1-P1 stabilisation"
      call assemble_kmk_matrix(local_state, lagrange_mesh, coordinate, theta_pg=1.0)    
      call add_kmk_matrix(local_state, cmc_m)
       
      if(associated(s_field)) then
        ! Should the timestep be passed in here? not sure
        call get_option("/timestepping/timestep", dt)
        call add_kmk_rhs(local_state, kmk_rhs, s_field, dt)
      end if
      
      call addto(projec_rhs, kmk_rhs)
    end if
    
    call scale(projec_rhs, -1.0)
    call addto(projec_rhs, ct_rhs)
    
    call impose_reference_pressure_node(cmc_m, projec_rhs, trim(l_option_path)//"/lagrange_multiplier")
    
    call petsc_solve(lagrange, cmc_m, projec_rhs)
    
    if(associated(s_field)) then
      call addto(s_field, lagrange)
    end if
    
    if(lump_mass) then
      call correct_masslumped_velocity(v_field, inverse_field_lumped_mass_vector, ct_m, lagrange)
    else if(dg) then
      call correct_velocity_dg(v_field, inverse_field_mass, ct_m, lagrange)
    else
      FLAbort("Not possible to not lump the mass if not dg.")
    end if
    
    call deallocate(ct_m_sparsity)
    call deallocate(ct_m)
    call deallocate(ct_rhs)
    call deallocate(cmc_m_sparsity)
    call deallocate(cmc_m)
    call deallocate(projec_rhs)
    call deallocate(lagrange)
    if(apply_kmk) then
      call deallocate(kmk_rhs)
    end if
    if(lump_mass) then
      call deallocate(field_lumped_mass)
      call deallocate(inverse_field_lumped_mass_vector)
    else if(dg) then
      call deallocate(inverse_field_mass)
    else
      FLAbort("Not possible to not lump the mass if not dg.")
    end if
    call deallocate(local_state)

  end subroutine
  
end module solenoidal_interpolation_module
