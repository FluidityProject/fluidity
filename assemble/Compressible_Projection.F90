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

module compressible_projection
  use fldebug
  use state_module
  use sparse_tools
  use spud
  use fields
  use sparse_matrices_fields
  use field_options
  use equation_of_state, only: compressible_eos, compressible_material_eos
  use global_parameters, only: OPTION_PATH_LEN
  use fefields, only: compute_lumped_mass
  use state_fields_module
  use sparsity_patterns_meshes
  use upwind_stabilisation
  implicit none 

  ! Buffer for output messages.
  character(len=255), private :: message

  private
  public :: assemble_compressible_projection_cv, assemble_compressible_projection_cg, update_compressible_density

  ! Stabilisation schemes
  integer, parameter :: STABILISATION_NONE = 0, &
    & STABILISATION_STREAMLINE_UPWIND = 1, STABILISATION_SUPG = 2
  ! Stabilisation scheme
  integer :: stabilisation_scheme
  integer :: nu_bar_scheme
  real :: nu_bar_scale

contains

  subroutine assemble_compressible_projection_cv(state, cmc, rhs, dt, theta_pg, theta_divergence, cmcget)

    ! inputs:
    ! bucket full of fields
    type(state_type), dimension(:), intent(inout) :: state

    type(csr_matrix), intent(inout) :: cmc
    type(scalar_field), intent(inout) :: rhs

    real, intent(in) :: dt
    real, intent(in) :: theta_pg, theta_divergence
    logical, intent(in) :: cmcget

    if((size(state)==1).and.(.not.has_scalar_field(state(1), "MaterialVolumeFraction"))) then
    
      call assemble_1mat_compressible_projection_cv(state(1), cmc, rhs, dt, &
                                                    theta_pg, theta_divergence, cmcget)
      
    else
    
      call assemble_mmat_compressible_projection_cv(state, cmc, rhs, dt, cmcget)
      
    end if
    

  end subroutine assemble_compressible_projection_cv

  subroutine assemble_1mat_compressible_projection_cv(state, cmc, rhs, dt, &
                                                      theta_pg, theta_divergence, cmcget)

    ! inputs:
    ! bucket full of fields
    type(state_type), intent(inout) :: state

    type(csr_matrix), intent(inout) :: cmc
    type(scalar_field), intent(inout) :: rhs

    real, intent(in) :: dt
    real, intent(in) :: theta_pg, theta_divergence
    logical, intent(in) :: cmcget

    ! local:
    integer :: norm_stat
    character(len=FIELD_NAME_LEN) :: normalisation_field

    type(scalar_field) :: eospressure, drhodp
    type(scalar_field), pointer :: normalisation, &
                                   density, olddensity
    type(scalar_field), pointer :: pressure
    type(scalar_field), pointer :: p_lumpedmass
    type(scalar_field) :: lhsfield, invnorm, absrhs
    
    type(scalar_field), pointer :: source, absorption
    integer :: stat

    real :: atmospheric_pressure, theta

    ewrite(1,*) 'Entering assemble_1mat_compressible_projection_cv'
    
    call zero(rhs)

    ! only do all this if we need to make cmc (otherwise we'd be adding repeatedly)
    if(cmcget) then

      pressure=>extract_scalar_field(state, "Pressure")
      call get_option(trim(pressure%option_path)//'/prognostic/atmospheric_pressure', &
                      atmospheric_pressure, default=0.0)
      
      if(pressure%mesh%shape%degree>1) then
        ! try lumping on the submesh
        p_lumpedmass => get_lumped_mass_on_submesh(state, pressure%mesh)
      else
        ! find the lumped mass
        p_lumpedmass => get_lumped_mass(state, pressure%mesh)
      end if
      ewrite_minmax(p_lumpedmass)
      
      call get_option(trim(pressure%option_path)//"/prognostic/scheme/use_compressible_projection_method/normalisation/name", &
                      normalisation_field, stat=norm_stat)
      if(norm_stat==0) then
        normalisation=>extract_scalar_field(state, trim(normalisation_field))
      else
        allocate(normalisation)
        call allocate(normalisation, pressure%mesh, "DummyNormalisation", field_type=FIELD_TYPE_CONSTANT)
        call set(normalisation, 1.0)
      end if
      
      call allocate(invnorm, normalisation%mesh, "InverseNormalisation", field_type=normalisation%field_type)
      call invert(normalisation, invnorm)
      
      if(norm_stat/=0) then
        call deallocate(normalisation)
        deallocate(normalisation)
      end if

      call allocate(lhsfield, pressure%mesh, "LHSField")

      call allocate(eospressure, pressure%mesh, 'EOSPressure')
      call allocate(drhodp, pressure%mesh, 'DerivativeDensityWRTBulkPressure')

      call zero(eospressure)
      call zero(drhodp)

      call compressible_eos(state, pressure=eospressure, drhodp=drhodp)

      density=>extract_scalar_field(state,'Density')
      ewrite_minmax(density)
      olddensity=>extract_scalar_field(state,'OldDensity')
      ewrite_minmax(olddensity)

      call get_option(trim(density%option_path)//"/prognostic/temporal_discretisation/theta", theta)

      call set(lhsfield, p_lumpedmass)
      call scale(lhsfield, drhodp)
      call scale(lhsfield, invnorm)
      call addto_diag(cmc, lhsfield, scale=1./(dt*dt*theta_divergence*theta_pg))
      
!     rhs = invnorm*p_lumpedmass* &
!      ( (1./dt)*(olddensity - density + drhodp*(eospressure - (pressure + atmospheric_pressure)))
!       +(absorption)*(drhodp*theta_pg*(eospressure - (pressure + atmospheric_pressure)) - theta_pg*density - (1-theta_pg)*olddensity)
!       +source)
      call set(rhs, pressure)
      call addto(rhs, atmospheric_pressure)
      call scale(rhs, -1.0)
      call addto(rhs, eospressure)
      call scale(rhs, drhodp)
      call addto(rhs, density, -1.0)
      call addto(rhs, olddensity)
      call scale(rhs, (1./dt))
      
      source => extract_scalar_field(state, "DensitySource", stat=stat)
      if(stat==0) then
        call addto(rhs, source)
      end if
      
      absorption => extract_scalar_field(state, "DensityAbsorption", stat=stat)
      if(stat==0) then
        call allocate(absrhs, absorption%mesh, "AbsorptionRHS")
        
        call set(absrhs, pressure)
        call addto(absrhs, atmospheric_pressure)
        call scale(absrhs, -1.0)
        call addto(absrhs, eospressure)
        call scale(absrhs, drhodp)
        call scale(absrhs, theta)
        call addto(absrhs, density, -theta)
        call addto(absrhs, olddensity, -(1-theta))
        call scale(absrhs, absorption)
        
        call addto(rhs, absrhs)
        
        call deallocate(absrhs)
        
        call scale(lhsfield, absorption)
        call addto_diag(cmc, lhsfield, scale=(theta/(dt*theta_divergence*theta_pg)))
      end if
      
      call scale(rhs, p_lumpedmass)
      call scale(rhs, invnorm)
      
      call deallocate(eospressure)
      call deallocate(drhodp)

      call deallocate(lhsfield)
      call deallocate(invnorm)

    end if

  end subroutine assemble_1mat_compressible_projection_cv

  subroutine assemble_mmat_compressible_projection_cv(state, cmc, rhs, dt, cmcget)

    ! inputs:
    ! bucket full of fields
    type(state_type), dimension(:), intent(inout) :: state

    type(csr_matrix), intent(inout) :: cmc
    type(scalar_field), intent(inout) :: rhs

    real, intent(in) :: dt
    logical, intent(in) :: cmcget

    ! local:
    integer :: i, stat, norm_stat
    character(len=OPTION_PATH_LEN) :: pressure_option_path
    character(len=FIELD_NAME_LEN) :: normalisation_field

    type(scalar_field) :: materialpressure, materialdrhodp, normdensity, &
                          normolddensity, normmatdrhodpp, normdrhodp
    type(scalar_field), pointer :: normalisation, &
                                   volumefraction, oldvolumefraction, materialdensity, oldmaterialdensity
    type(scalar_field), pointer :: dummy_ones

    type(scalar_field), pointer :: pressure
    type(vector_field), pointer :: positions
    type(scalar_field) :: lumped_mass, tempfield

    real :: atmospheric_pressure

    ewrite(1,*) 'Entering assemble_mmat_compressible_projection_cv'

    pressure=>extract_prognostic_pressure(state, stat=stat)
    if(stat/=0) then
       ! how did we end up here?
       FLAbort("In assemble_mmat_compressible_projection_cv without a pressure")
    end if
    pressure_option_path=trim(pressure%option_path)
    
    call zero(rhs)
   
    if(have_option(trim(pressure_option_path)//"/prognostic/scheme/use_compressible_projection_method")) THEN

      ! only do all this if we need to make cmc (otherwise we'd be adding repeatedly)
      if(cmcget) then

        positions=>extract_vector_field(state(1), "Coordinate")
        call allocate(lumped_mass, pressure%mesh, "LumpedMassField")
        call allocate(tempfield, pressure%mesh, "TemporaryAssemblyField")
        call compute_lumped_mass(positions, lumped_mass)

        allocate(dummy_ones)
        call allocate(dummy_ones, pressure%mesh, "DummyOnesField")
        call set(dummy_ones, 1.0)

        call get_option(trim(pressure_option_path)//'/prognostic/atmospheric_pressure', &
                        atmospheric_pressure, default=0.0)

        call get_option(trim(pressure_option_path)//"/prognostic/scheme/use_compressible_projection_method/normalisation/name", &
                        normalisation_field, stat=norm_stat)

        call allocate(materialpressure, pressure%mesh, 'MaterialEOSPressure')
        call allocate(materialdrhodp, pressure%mesh, 'DerivativeMaterialdensityWRTBulkPressure')

        call allocate(normdensity, pressure%mesh, 'NormalisedMaterialDensity')
        call allocate(normolddensity, pressure%mesh, 'NormalisedOldMaterialDensity')
        call allocate(normmatdrhodpp, pressure%mesh, 'NormalisedMaterialPressure')
        call allocate(normdrhodp, pressure%mesh, 'NormalisedDrhodp')

        normdensity%val = 0.0
        normolddensity%val = 0.0
        normmatdrhodpp%val = 0.0
        normdrhodp%val=0.0

        do i = 1,size(state)

          materialpressure%val=0.0
          materialdrhodp%val=0.0

          call compressible_material_eos(state(i), materialpressure=materialpressure, materialdrhodp=materialdrhodp)

          volumefraction=>extract_scalar_field(state(i),'MaterialVolumeFraction', stat=stat)
          if(stat==0) then
            oldvolumefraction=>extract_scalar_field(state(i),'OldMaterialVolumeFraction')
            materialdensity=>extract_scalar_field(state(i),'MaterialDensity')
            oldmaterialdensity=>extract_scalar_field(state(i),'OldMaterialDensity')

            if(norm_stat==0) then
              normalisation=>extract_scalar_field(state(i), trim(normalisation_field))
            else
              normalisation=>dummy_ones
            end if

            normdensity%val = normdensity%val &
                              + materialdensity%val*volumefraction%val/ &
                                normalisation%val
            normolddensity%val = normolddensity%val &
                                + oldmaterialdensity%val*oldvolumefraction%val/ &
                                  normalisation%val
            normmatdrhodpp%val = normmatdrhodpp%val &
                                  + materialpressure%val*materialdrhodp%val*volumefraction%val/ &
                                    normalisation%val
            normdrhodp%val = normdrhodp%val &
                              + materialdrhodp%val*volumefraction%val/ &
                                normalisation%val
          endif

        end do

        call zero(tempfield)
        tempfield%val = (1./(dt*dt))*lumped_mass%val*normdrhodp%val

        call addto_diag(cmc, tempfield)

        rhs%val = (1./dt)*lumped_mass%val* &
                          ( &
                            normolddensity%val &
                          - normdensity%val &
                          ) &
               +(1./dt)*lumped_mass%val* &
                          ( &
                            normmatdrhodpp%val &
                          - normdrhodp%val*(pressure%val+atmospheric_pressure) &
                          )

        call deallocate(normdensity)
        call deallocate(normolddensity)
        call deallocate(normmatdrhodpp)
        call deallocate(normdrhodp)

        call deallocate(materialpressure)
        call deallocate(materialdrhodp)

        call deallocate(lumped_mass)
        call deallocate(tempfield)
        call deallocate(dummy_ones)
        deallocate(dummy_ones)

      end if

    end if

  end subroutine assemble_mmat_compressible_projection_cv
  
  subroutine assemble_compressible_projection_cg(state, cmc, rhs, dt, theta_pg, theta_divergence, assemble_cmc, reuse_cmc)

    ! assemble the "compressible parts" of the projection equation
    ! i.e. add in the drhodp mass term to cmc (C^T M^-1 C has been assembled before)
    ! and assemble the right hand side terms; rhs here starts from zero, any additional terms (e.g. weak bcs) need to be added in outside this routine 
    type(state_type), dimension(:), intent(inout) :: state

    type(csr_matrix), intent(inout) :: cmc
    type(scalar_field), intent(inout) :: rhs

    real, intent(in) :: dt
    real, intent(in) :: theta_pg, theta_divergence
    ! this indicates whether cmc has just been reassembled and the mass terms
    ! are added directly. If .false. (only in combination with reuse_cmc==.true.)
    ! cmc is reused from a previous non-linear iteration. The mass terms
    ! from this previous iteration are subtracted first
    logical, intent(in) :: assemble_cmc
    ! if .true. store the mass terms also separately in state, so we can
    ! subtract them out the next time
    logical, intent(in) :: reuse_cmc

    if((size(state)==1).and.(.not.has_scalar_field(state(1), "MaterialVolumeFraction"))) then
    
      call assemble_1mat_compressible_projection_cg(state(1), cmc, rhs, dt, &
                                                    theta_pg, theta_divergence, assemble_cmc, reuse_cmc)
      
    else
      
        FLExit("Multimaterial compressible continuous_galerkin pressure not possible.")
      
    end if
    

  end subroutine assemble_compressible_projection_cg

  subroutine assemble_1mat_compressible_projection_cg(state, cmc, rhs, dt, &
                                                      theta_pg, theta_divergence, &
                                                      assemble_cmc, reuse_cmc)

    ! assemble the "compressible parts" of the projection equation
    ! i.e. add in the drhodp mass term to cmc (C^T M^-1 C has been assembled before)
    ! and assemble the right hand side terms; rhs here starts from zero, any additional terms (e.g. weak bcs) need to be added in outside this routine 
    type(state_type), intent(inout) :: state

    type(csr_matrix), intent(inout), target :: cmc
    type(scalar_field), intent(inout) :: rhs

    real, intent(in) :: dt
    real, intent(in) :: theta_pg, theta_divergence
    ! this indicates whether cmc has just been reassembled and the mass terms
    ! are added directly. If .false. (only in combination with reuse_cmc==.true.)
    ! cmc is reused from a previous non-linear iteration. The mass terms
    ! from this previous iteration are subtracted first
    logical, intent(in) :: assemble_cmc
    ! if .true. store the mass terms also separately in state, so we can
    ! subtract them out the next time
    logical, intent(in) :: reuse_cmc

    ! local
    type(mesh_type), pointer :: test_mesh

    type(vector_field), pointer :: field

    integer, dimension(:), pointer :: test_nodes

    real, dimension(:), allocatable :: ele_rhs
    type(element_type), pointer :: test_shape_ptr
    type(element_type) :: test_shape
    real, dimension(:,:,:), allocatable :: dtest_t
    real, dimension(:), allocatable :: detwei
    real, dimension(:,:,:), allocatable :: j_mat
    
    real, dimension(:), allocatable :: density_at_quad, olddensity_at_quad, p_at_quad, &
                                      drhodp_at_quad, eosp_at_quad, abs_at_quad
    real, dimension(:,:), allocatable :: nlvelocity_at_quad

    ! pointer to fields in state
    type(vector_field), pointer :: coordinate, nonlinearvelocity, velocity
    type(scalar_field), pointer :: pressure, density, olddensity
    type(scalar_field), pointer :: source, absorption

    type(csr_matrix), pointer:: projection_mass
    type(csr_sparsity), pointer:: projection_mass_sparsity
    type(scalar_field) :: eospressure, drhodp
    real :: theta, atmospheric_pressure

    real, dimension(:,:), allocatable :: ele_mat
    integer :: ele
    
    logical :: have_absorption, have_source 
    integer :: stat

    ewrite(1,*) 'Entering assemble_1mat_compressible_projection_cg'
    
    call zero(rhs)

    coordinate=> extract_vector_field(state, "Coordinate")
    
    density => extract_scalar_field(state, "Density")
    olddensity => extract_scalar_field(state, "OldDensity")
    
    absorption => extract_scalar_field(state, "DensityAbsorption", stat=stat)
    have_absorption = (stat==0)
    if(have_absorption) then
      ewrite(2,*) 'Have DensityAbsorption'
    end if
    
    source => extract_scalar_field(state, "DensitySource", stat=stat)
    have_source = (stat==0)
    if(have_source) then
      ewrite(2,*) 'Have DensitySource'
    end if

    velocity=>extract_vector_field(state, "Velocity")
    nonlinearvelocity=>extract_vector_field(state, "NonlinearVelocity") ! maybe this should be updated after the velocity solve?
    
    pressure => extract_scalar_field(state, "Pressure")

    call get_option(trim(pressure%option_path)//'/prognostic/atmospheric_pressure', &
                    atmospheric_pressure, default=0.0)

    ! these are put on the density mesh, which should be of sufficient order to represent
    ! the multiplication of the eos (of course that may not be possible in which case
    ! something should be done at the gauss points instead)
    call allocate(eospressure, density%mesh, 'EOSPressure')
    call allocate(drhodp, density%mesh, 'DerivativeDensityWRTBulkPressure')

    call zero(eospressure)
    call zero(drhodp)

    ! this needs to be changed to be evaluated at the quadrature points!
    call compressible_eos(state, pressure=eospressure, drhodp=drhodp)

    ewrite_minmax(density)
    ewrite_minmax(olddensity)

    if(have_option(trim(density%option_path) // &
                        "/prognostic/spatial_discretisation/continuous_galerkin/&
                        &stabilisation/streamline_upwind_petrov_galerkin")) then
      ewrite(2, *) "SUPG stabilisation"
      stabilisation_scheme = STABILISATION_SUPG
      call get_upwind_options(trim(density%option_path) // & 
                              "/prognostic/spatial_discretisation/continuous_galerkin/&
                              &stabilisation/streamline_upwind_petrov_galerkin", &
                              & nu_bar_scheme, nu_bar_scale)
    else
      ewrite(2, *) "No stabilisation"
      stabilisation_scheme = STABILISATION_NONE
    end if

    call get_option(trim(density%option_path)//"/prognostic/temporal_discretisation/theta", theta)

    if (reuse_cmc) then
      ewrite(2,*) "Trying to reuse cmc"
      projection_mass => extract_csr_matrix(state, "CompressibleProjectionMassMatrix", stat=stat)
      
      if (stat/=0) then
        ewrite(2,*) "New assembly of cmc, creating separate mass matrix"
        ! no mass matrix left from previous iteration
        if (.not. assemble_cmc) then
          FLAbort("Told to reuse cmc but can't find mass matrix from previous iteration")
        end if
        projection_mass_sparsity => get_csr_sparsity_firstorder(state, pressure%mesh, pressure%mesh)
        call allocate(projection_mass, projection_mass_sparsity, name="CompressibleProjectionMassMatrix")
        call insert(state, projection_mass, "CompressibleProjectionMassMatrix")
        call deallocate(projection_mass)
        projection_mass => extract_csr_matrix(state, "CompressibleProjectionMassMatrix", stat=stat)
      else if (.not. assemble_cmc) then
        ewrite(2,*) "Subtracting mass terms from previous iteration"
        ! cmc has not been reassembled, so we need to subtract the mass matrix from the previous iteration
        call addto(cmc, projection_mass, scale=-1.0)
      end if

      ! we can now forget about the mass from previous iteration, and start assembling
      ! in the current mass term in the matrix
      call zero(projection_mass)

    else
      ! not reusing cmc, so we can add the mass term straight into cmc
      projection_mass => cmc
    end if

    test_mesh => pressure%mesh
    field => velocity
    
    allocate(dtest_t(ele_loc(test_mesh, 1), ele_ngi(test_mesh, 1), field%dim), &
            detwei(ele_ngi(field, 1)), &
            ele_mat(ele_loc(test_mesh, 1), ele_loc(test_mesh, 1)), &
            ele_rhs(ele_loc(test_mesh, 1)), &
            density_at_quad(ele_ngi(density, 1)), &
            olddensity_at_quad(ele_ngi(density, 1)), &
            nlvelocity_at_quad(field%dim, ele_ngi(field, 1)), &
            j_mat(field%dim, field%dim, ele_ngi(density, 1)), &
            drhodp_at_quad(ele_ngi(drhodp, 1)), &
            eosp_at_quad(ele_ngi(eospressure, 1)), &
            abs_at_quad(ele_ngi(density, 1)), &
            p_at_quad(ele_ngi(pressure, 1)))
    
    do ele=1, element_count(test_mesh)
    
      test_nodes=>ele_nodes(test_mesh, ele)

      test_shape_ptr => ele_shape(test_mesh, ele)
      
      density_at_quad = ele_val_at_quad(density, ele)
      olddensity_at_quad = ele_val_at_quad(olddensity, ele)
      
      p_at_quad = ele_val_at_quad(pressure, ele) + atmospheric_pressure
                        
      nlvelocity_at_quad = ele_val_at_quad(nonlinearvelocity, ele)
      
      drhodp_at_quad = ele_val_at_quad(drhodp, ele)
      eosp_at_quad = ele_val_at_quad(eospressure, ele)
      
      select case(stabilisation_scheme)
        case(STABILISATION_SUPG)
          call transform_to_physical(coordinate, ele, test_shape_ptr, dshape = dtest_t, detwei=detwei, j=j_mat)
          test_shape = make_supg_shape(test_shape_ptr, dtest_t, nlvelocity_at_quad, j_mat, &
            & nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
        case default
          call transform_to_physical(coordinate, ele, detwei=detwei)
          test_shape = test_shape_ptr
          call incref(test_shape)
      end select
      ! Important note: with SUPG the test function derivatives have not been
      ! modified.

      ele_mat = 1./(dt*dt*theta_divergence*theta_pg) * shape_shape(test_shape, test_shape_ptr, &
          detwei*drhodp_at_quad)
      !       /
      ! rhs = |test_shape* &
      !       /
      !      ((1./dt)*(drhodp*(eospressure - (pressure + atmospheric_pressure)) + olddensity - density)
      ! +(absorption)*(drhodp*theta*(eospressure - (pressure + atmospheric_pressure)) 
      !                - theta*density - (1-theta)*olddensity)
      ! +source)dV
      ele_rhs = (1./dt)*shape_rhs(test_shape, detwei*((drhodp_at_quad*(eosp_at_quad - p_at_quad)) &
                                                     +(olddensity_at_quad - density_at_quad)))
      
      if(have_source) then
        ele_rhs = ele_rhs + shape_rhs(test_shape, detwei*ele_val_at_quad(source, ele))
      end if
      
      if(have_absorption) then
        abs_at_quad = ele_val_at_quad(absorption, ele)
        ele_mat = ele_mat + (theta/(theta_divergence*theta_pg))* &
            shape_shape(test_shape, test_shape_ptr,detwei*drhodp_at_quad*dt*abs_at_quad)
        ele_rhs = ele_rhs + &
                  shape_rhs(test_shape, detwei*abs_at_quad*(theta*(drhodp_at_quad*(eosp_at_quad - p_at_quad)-density_at_quad) &
                                                           -(1-theta)*olddensity_at_quad))
      end if
      
      call addto(projection_mass, test_nodes, test_nodes, ele_mat)
      
      call addto(rhs, test_nodes, ele_rhs)

      call deallocate(test_shape)
      
    end do

    if (reuse_cmc) then
      ! we've been adding into projection_mass seperately
      ! so now add this into cmc itself
      ewrite(2,*) "Adding mass terms from separate mass matrix into cmc"
      call addto(cmc, projection_mass)
    end if

    call deallocate(drhodp)
    call deallocate(eospressure)

  end subroutine assemble_1mat_compressible_projection_cg

  subroutine update_compressible_density(state)
  
    type(state_type), dimension(:), intent(inout) :: state
    
    type(scalar_field), pointer :: density
    
    if((size(state)==1).and.(.not.has_scalar_field(state(1), "MaterialVolumeFraction"))) then
    
      density=>extract_scalar_field(state(1),'Density')
      
      if(have_option(trim(density%option_path)//"/prognostic")) then
        
        call compressible_eos(state(1), density=density)
      
      end if
    
    end if
  
  end subroutine update_compressible_density

end module compressible_projection

