#include "fdebug.h"
! #define CHECK_SUPERMESH_VOLUMES
! #define DUMP_SUPERMESH_INTERSECTIONS

module solenoidal_interpolation_module

  use fldebug
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
  use futils
  use spud
  use sparse_tools
  use sparse_tools_petsc
  use vector_tools
  use tensors
  use element_numbering, only: FAMILY_SIMPLEX
  use transform_elements
  use linked_lists
  use supermesh_construction
  use fetools
  use fields
  use state_module
  use field_options, only : complete_field_path
  use sparsity_patterns
  use boundary_conditions
  use interpolation_module
  use sparse_matrices_fields
  use solvers
  use full_projection
  use fefields
  use dgtools
  use assemble_cmc, only: assemble_cmc_dg, repair_stiff_nodes,&
     zero_stiff_nodes, assemble_masslumped_cmc, assemble_diagonal_schur
  use boundary_conditions_from_options
  use divergence_matrix_cv, only: assemble_divergence_matrix_cv
  use divergence_matrix_cg, only: assemble_divergence_matrix_cg
  use momentum_cg, only: correct_velocity_cg, correct_masslumped_velocity, &
                         add_kmk_matrix, add_kmk_rhs, assemble_kmk_matrix
  use momentum_dg, only: correct_velocity_dg
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
      
      if(have_option(trim(complete_field_path(v_field%option_path, stat=stat))//&
                     "/enforce_discrete_properties/solenoidal/lagrange_multiplier/update_scalar_field")) then
        call get_option(trim(complete_field_path(v_field%option_path, stat=stat))//&
                     "/enforce_discrete_properties/solenoidal/lagrange_multiplier/update_scalar_field/name", &
                     update_field_name)
        s_field=> extract_scalar_field(state, trim(update_field_name))
      else
        s_field=>null()
      end if
      
      call get_option(trim(complete_field_path(v_field%option_path, stat=stat))//&
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
    
    logical :: dg, discontinuous, div_cv, div_cg, apply_kmk
    logical :: assemble_cmc, assemble_lagrange_mass, assemble_cdiagmc
    logical :: assemble_lumped_mass, lump_mass_on_submesh
    logical :: full_schur, assemble_schur_aux
    integer :: dim, j
    real :: dt
    character(len=FIELD_NAME_LEN) :: pressure_pmat
    
    type(block_csr_matrix), pointer :: ct_m
    type(block_csr_matrix), pointer :: ctp_m
    type(scalar_field) :: ct_rhs, kmk_rhs
    type(csr_sparsity) :: ct_m_sparsity, cmc_m_sparsity
    type(csr_sparsity) :: field_mass_sparsity, lagrange_mass_sparsity
    type(csr_matrix), target :: cmc_m, cdiagmc_m
    type(csr_matrix) :: schur_aux
    type(csr_matrix), target :: lagrange_mass
    type(petsc_csr_matrix), target :: field_mass
    type(csr_matrix), pointer :: pschur
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

    type(ilist) :: stiff_nodes_list
    logical :: stiff_nodes_repair
    
    ! for a CG Lagrange multiplier are we testing the divergence with the CV dual
    logical :: cg_lagrange_cv_test_divergence
    
    call insert(local_state, coordinate, "Coordinate")
    
    l_option_path=trim(complete_field_path(v_field%option_path))//"/enforce_discrete_properties/solenoidal"
    
    dim = mesh_dim(v_field)
    
    if(associated(s_field)) then
      assert(trim(s_field%mesh%name) == trim(lagrange_mesh%name))
    end if

    dg = (continuity(v_field) < 0)
    discontinuous = have_option(trim(l_option_path)//&
                                "/interpolated_field/discontinuous")

    if (dg) then
      assert(discontinuous)
    end if
    
    assemble_lumped_mass = have_option(trim(l_option_path)//&
                "/interpolated_field/discontinuous/lump_mass_matrix") &
                .or. &
                have_option(trim(l_option_path)//&
                "/interpolated_field/continuous/lump_mass_matrix") &
                .or. &
                have_option(trim(l_option_path)//&
                "/interpolated_field/continuous/full_schur_complement"//&
                "/preconditioner_matrix::LumpedSchurComplement")
    
    lump_mass_on_submesh = have_option(trim(l_option_path)//&
                      "/interpolated_field/continuous/lump_mass_matrix/use_submesh") &
                      .or. &
                      have_option(trim(l_option_path)//&
                      "/interpolated_field/continuous/full_schur_complement"//&
                      "/preconditioner_matrix[0]/lump_on_submesh")

    div_cv = have_option(trim(l_option_path) //&
          &"/lagrange_multiplier/spatial_discretisation/control_volumes")
    div_cg = .not.div_cv
    
    cg_lagrange_cv_test_divergence = have_option(trim(l_option_path) //&
          &"/lagrange_multiplier/spatial_discretisation/continuous_galerkin/test_divergence_with_cv_dual")

    stiff_nodes_repair = have_option(trim(l_option_path) //&
          &"/lagrange_multiplier/repair_stiff_nodes")

    apply_kmk = (continuity(lagrange_mesh) >= 0 .and. lagrange_mesh%shape%degree == 1 &
          & .and. lagrange_mesh%shape%numbering%family == FAMILY_SIMPLEX .and. &
          & continuity(v_field) >= 0 .and. v_field%mesh%shape%degree == 1 &
          & .and. v_field%mesh%shape%numbering%family == FAMILY_SIMPLEX .and. &
          & .not. have_option(trim(l_option_path) // &
          & "/lagrange_multiplier/spatial_discretisation/continuous_galerkin/remove_stabilisation_term") .and. &
          & .not. div_cv)

    full_schur = have_option(trim(l_option_path)//&
                 "/interpolated_field/continuous/full_schur_complement")

    assemble_cmc = .true.  ! true, unless we're using full_schur and even then, it's complicated...
    assemble_cdiagmc = .false.        ! only used as a possible preconditioner for full_schur
    assemble_lagrange_mass = .false.  ! only used as a possible preconditioner for full_schur
    assemble_schur_aux = .false.
    if(full_schur) then
      assemble_schur_aux = apply_kmk
      ! Check to see whether pressure cmc_m preconditioning matrix is needed:
      call get_option(trim(l_option_path)//&
                 "/interpolated_field/continuous/full_schur_complement"//&
                 "/preconditioner_matrix[0]/name", pressure_pmat)
   
      select case(pressure_pmat)
         case("LumpedSchurComplement")
            pschur => cmc_m
         case("DiagonalSchurComplement")
            assemble_cmc = .false.
            assemble_cdiagmc = .true.
            pschur => cdiagmc_m
         case("LagrangeMassMatrix")
            assemble_cmc = .false.
            assemble_lagrange_mass = .true.
            pschur => lagrange_mass
         case("NoPreconditionerMatrix")
            assemble_cmc = .false.
            nullify(pschur)
         case default
            ! Developer error... out of sync options input and code
            FLAbort("Unknown preconditioner matrix for full schur complement in solenoidal interpolation.")
      end select

      ! always need a field_mass matrix if doing full_schur
      field_mass_sparsity = make_sparsity(v_field%mesh, v_field%mesh, "FieldMassSparsity")
      call allocate(field_mass, field_mass_sparsity, (/v_field%dim, v_field%dim/), &
                    group_size=(/v_field%dim, v_field%dim/), &
                    diagonal=.true.,  name="FieldMassMatrix")
      call zero(field_mass)

    end if

    ct_m_sparsity = make_sparsity(lagrange_mesh, v_field%mesh, "DivergenceSparsity")
    allocate(ct_m)
    call allocate(ct_m, ct_m_sparsity, blocks=(/1, dim/), name="DivergenceMatrix")
    call zero(ct_m)
    
    ! If CG with CV tested divergence then we need to allocate the 
    ! left C matrix as it is formed via testing with CV so
    if (cg_lagrange_cv_test_divergence) then    
       allocate(ctp_m)
       call allocate(ctp_m, ct_m_sparsity, blocks=(/1, dim/), name="CVTestedDivergenceMatrix")
       call zero(ctp_m)
    else
       ctp_m => ct_m
    end if
    
    call allocate(ct_rhs, lagrange_mesh, "DivergenceRHS")
    call zero(ct_rhs)

    if (assemble_cmc.or.assemble_schur_aux.or.assemble_cdiagmc) then
      cmc_m_sparsity = make_sparsity_transpose(lagrange_mesh, v_field%mesh, "LagrangeProjectionSparsity")
    end if
    
    if (assemble_cmc) then
      call allocate(cmc_m, cmc_m_sparsity, name="LagrangeProjectionMatrix")
      call zero(cmc_m)
    end if

    if (assemble_schur_aux) then
      call allocate(schur_aux, cmc_m_sparsity, name="SchurAuxilliaryMatrix")
      call zero(schur_aux)
    end if
    
    if (assemble_cdiagmc) then
      call allocate(cdiagmc_m, cmc_m_sparsity, name="DiagonalSchurMatrix")
      call zero(cdiagmc_m)
    end if

    if (assemble_lagrange_mass) then
      lagrange_mass_sparsity = make_sparsity(lagrange_mesh, lagrange_mesh, "LagrangeMassSparsity")
      call allocate(lagrange_mass, lagrange_mass_sparsity, name="LagrangeMassMatrix")
      call zero(lagrange_mass)
    end if

    call allocate(projec_rhs, lagrange_mesh, "LagrangeProjectionRHS")
    call zero(projec_rhs)
    
    call allocate(lagrange, lagrange_mesh, "LagrangianMultiplier")
    call zero(lagrange)
    lagrange%option_path = trim(l_option_path)//"/lagrange_multiplier"

    if(apply_kmk) then
      call allocate(kmk_rhs, lagrange_mesh, "KMKRHS")
      call zero(kmk_rhs)
    end if

    if(assemble_lumped_mass) then
      call allocate(field_lumped_mass, v_field%mesh, "FieldLumpedMass")
      call zero(field_lumped_mass)
      call allocate(inverse_field_lumped_mass_vector, dim, v_field%mesh, &
         "InverseFieldLumpedMassVector")
    end if

    if(full_schur) then
      if(assemble_lagrange_mass) then
        call compute_mass(coordinate, lagrange_mesh, lagrange_mass)
      end if
      if (assemble_lumped_mass.and.(.not.lump_mass_on_submesh)) then
        call compute_mass(coordinate, v_field%mesh, field_mass, field_lumped_mass)
      else
        call compute_mass(coordinate, v_field%mesh, field_mass)
      end if
    else
      if(assemble_lumped_mass.and.(.not.lump_mass_on_submesh)) then
        call compute_lumped_mass(coordinate, field_lumped_mass)
      end if
    end if

    if(assemble_lumped_mass.and.lump_mass_on_submesh) then
      ! get the mass matrix lumped on the submesh
      call compute_lumped_mass_on_submesh(local_state, field_lumped_mass)
    else if (discontinuous) then
      ! get the inverse discontinuous mass matrix
      call construct_inverse_mass_matrix_dg(inverse_field_mass, v_field, coordinate)
    end if
      
    if(div_cg) then
      call assemble_divergence_matrix_cg(ct_m, local_state, ct_rhs=ct_rhs, &
                                          test_mesh=lagrange_mesh, field=v_field, &
                                          option_path = trim(l_option_path)//"/lagrange_multiplier")

      ! If CG lagrange with CV tested divergence then form the other C matrix.
      ! This will overwrite the ct_rhs formed above.
      if (cg_lagrange_cv_test_divergence) then
         call assemble_divergence_matrix_cv(ctp_m, local_state, ct_rhs=ct_rhs, &
                                            test_mesh=lagrange_mesh, field=v_field)
      end if
            
    else if(div_cv) then
      call assemble_divergence_matrix_cv(ct_m, local_state, ct_rhs=ct_rhs, &
                                         test_mesh=lagrange_mesh, field=v_field)

    else
      ! coding error
      FLAbort("Unknown spatial discretisation option for the lagrange multiplier.")
    end if
    
    if(assemble_cmc) then
      if(assemble_lumped_mass) then
        call invert(field_lumped_mass)
        
        do j=1, inverse_field_lumped_mass_vector%dim
          call set(inverse_field_lumped_mass_vector, j, field_lumped_mass)
        end do
        
        call apply_dirichlet_conditions_inverse_mass(inverse_field_lumped_mass_vector, v_field)

        call assemble_masslumped_cmc(cmc_m, ctp_m, inverse_field_lumped_mass_vector, ct_m)
      else if(dg) then
        call assemble_cmc_dg(cmc_m, ctp_m, ct_m, inverse_field_mass)
      else
        FLExit("Not possible to not lump the mass if not dg or full_schur.")
      end if

      if(stiff_nodes_repair) then
        call repair_stiff_nodes(cmc_m, stiff_nodes_list)
      end if
    end if

    if(assemble_cdiagmc) then
      call assemble_diagonal_schur(cdiagmc_m, v_field, field_mass, ctp_m, ct_m)
    end if
    
    call mult(projec_rhs, ctp_m, v_field)
    
    if(apply_kmk) then
      ! a hack to make sure the appropriate meshes are available for the
      ! construction of the sparsity
      call insert(local_state, lagrange, name="Pressure")
      call insert(local_state, v_field, name="Velocity")
      ! end of hack... you can look again now
      
      ewrite(2,*) "Assembling P1-P1 stabilisation"
      call assemble_kmk_matrix(local_state, lagrange_mesh, coordinate, theta_pg=1.0)    
      if(assemble_cmc) then
        call add_kmk_matrix(local_state, cmc_m)
      end if
      if(assemble_schur_aux) then
        call add_kmk_matrix(local_state, schur_aux)
      end if
      if(assemble_cdiagmc) then
        call add_kmk_matrix(local_state, cdiagmc_m)
      end if
       
      if(associated(s_field)) then
        ! Should the timestep be passed in here? not sure
        call get_option("/timestepping/timestep", dt)
        call add_kmk_rhs(local_state, kmk_rhs, s_field, dt)
      end if
      
      ! clean up our mess
      call remove_scalar_field(local_state, name="Pressure")
      call remove_vector_field(local_state, name="Velocity")
      
      call addto(projec_rhs, kmk_rhs)
    end if
    
    call scale(projec_rhs, -1.0)
    call addto(projec_rhs, ct_rhs)
    
    if (assemble_cmc) then
      call impose_reference_pressure_node(cmc_m, projec_rhs, coordinate, trim(l_option_path)//"/lagrange_multiplier")
    
      ! only have a stiff_nodes_list if we're assembling cmc
      if(stiff_nodes_repair) then
        call zero_stiff_nodes(projec_rhs, stiff_nodes_list)
      end if
    end if

    if (full_schur) then
      if (apply_kmk) then
        call petsc_solve_full_projection(lagrange, ctp_m, field_mass, ct_m, &
                                         projec_rhs, pschur, v_field, &
                                         local_state, v_field%mesh, &
                                         option_path=trim(l_option_path)//"/lagrange_multiplier",&
                                         inner_option_path=trim(l_option_path)//&
                                                           "/interpolated_field/continuous/full_schur_complement/inner_matrix[0]",&
                                         auxiliary_matrix=schur_aux)
      else
        call petsc_solve_full_projection(lagrange, ctp_m, field_mass, ct_m, &
                                         projec_rhs, pschur, v_field, &
                                         local_state, v_field%mesh, &
                                         option_path=trim(l_option_path)//"/lagrange_multiplier",&
                                         inner_option_path=trim(l_option_path)//&
                                                           "/interpolated_field/continuous/full_schur_complement/inner_matrix[0]")
      end if
    else
      call petsc_solve(lagrange, cmc_m, projec_rhs)
    end if
    
    if(associated(s_field)) then
      call addto(s_field, lagrange)
    end if
    
    if(full_schur) then
      call correct_velocity_cg(v_field, field_mass, ct_m, lagrange, local_state)
    else if(assemble_lumped_mass) then
      call correct_masslumped_velocity(v_field, inverse_field_lumped_mass_vector, ct_m, lagrange)
    else if(discontinuous) then
      call correct_velocity_dg(v_field, inverse_field_mass, ct_m, lagrange)
    end if
    
    if (full_schur) then
      call deallocate(field_mass)
    end if
    call deallocate(ct_m_sparsity)
    call deallocate(ct_m)
    deallocate(ct_m)
    if (cg_lagrange_cv_test_divergence) then
      call deallocate(ctp_m)
      deallocate(ctp_m)
    end if
    call deallocate(ct_rhs)
    if (assemble_cmc.or.assemble_schur_aux.or.assemble_cdiagmc) then
      call deallocate(cmc_m_sparsity)
    end if
    if (assemble_cmc) then
      call deallocate(cmc_m)
    end if
    if (assemble_cdiagmc) then
      call deallocate(cdiagmc_m)
    end if
    if (assemble_schur_aux) then
      call deallocate(schur_aux)
    end if
    if (assemble_lagrange_mass) then
      call deallocate(lagrange_mass)
    end if
    call deallocate(projec_rhs)
    call deallocate(lagrange)
    if(apply_kmk) then
      call deallocate(kmk_rhs)
    end if
    if(assemble_lumped_mass) then
      call deallocate(field_lumped_mass)
      call deallocate(inverse_field_lumped_mass_vector)
    else if(discontinuous) then
      call deallocate(inverse_field_mass)
    end if
    call deallocate(local_state)

  end subroutine
  
end module solenoidal_interpolation_module
