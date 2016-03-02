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
  use spud
  use futils, only: int2str
  use global_parameters, only: OPTION_PATH_LEN
  use sparse_tools
  use elements
  use transform_elements
  use fields
  use state_module
  use fetools, only: shape_shape, shape_rhs
  use sparse_matrices_fields
  use field_options
  use fefields, only: compute_cv_mass
  use state_fields_module
  use equation_of_state, only: compressible_eos, compressible_material_eos
  use upwind_stabilisation
  use multiphase_module
  implicit none 

  ! Buffer for output messages.
  character(len=255), private :: message

  private
  public :: assemble_compressible_projection_cv, assemble_compressible_projection_cg, update_compressible_density
  public :: compressible_projection_check_options

  ! Stabilisation schemes
  integer, parameter :: STABILISATION_NONE = 0, &
    & STABILISATION_STREAMLINE_UPWIND = 1, STABILISATION_SUPG = 2
  ! Stabilisation scheme
  integer :: stabilisation_scheme
  integer :: nu_bar_scheme
  real :: nu_bar_scale

  !! Are we running a multiphase flow simulation?
  logical :: multiphase

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
    type(scalar_field) :: eospressure, drhodp
    type(scalar_field), pointer :: density, olddensity
    type(scalar_field), pointer :: pressure
    type(scalar_field), pointer :: p_cvmass
    type(scalar_field) :: lhsfield, absrhs
    
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
      
      ! find the cv mass
      p_cvmass => get_cv_mass(state, pressure%mesh)

      ewrite_minmax(p_cvmass)
      
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

      call set(lhsfield, p_cvmass)
      call scale(lhsfield, drhodp)
      call addto_diag(cmc, lhsfield, scale=1./(dt*dt*theta_divergence*theta_pg))
      
!     rhs = p_cvmass* &
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
      
      call scale(rhs, p_cvmass)
      
      call deallocate(eospressure)
      call deallocate(drhodp)

      call deallocate(lhsfield)

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
    integer :: i, stat
    character(len=OPTION_PATH_LEN) :: pressure_option_path

    type(scalar_field) :: materialpressure, materialdrhodp, density, &
                          olddensity, matdrhodpp, drhodp
    type(scalar_field), pointer :: volumefraction, oldvolumefraction, materialdensity, oldmaterialdensity
    type(scalar_field), pointer :: dummy_ones

    type(scalar_field), pointer :: pressure
    type(vector_field), pointer :: positions
    type(scalar_field) :: cv_mass, tempfield

    real :: atmospheric_pressure
    ! Do we want to use the compressible projection method?
    logical :: have_compressible_eos

    ewrite(1,*) 'Entering assemble_mmat_compressible_projection_cv'

    pressure=>extract_prognostic_pressure(state, stat=stat)
    if(stat/=0) then
       ! how did we end up here?
       FLAbort("In assemble_mmat_compressible_projection_cv without a pressure")
    end if
    pressure_option_path=trim(pressure%option_path)
    
    have_compressible_eos = .false.
    state_loop: do i = 1, size(state)
      have_compressible_eos = have_option("/material_phase::"//trim(state(i)%name)//"/equation_of_state/compressible")
      if(have_compressible_eos) then
        exit state_loop
      end if
    end do state_loop

    call zero(rhs)
   
    if(have_compressible_eos) THEN

      ! only do all this if we need to make cmc (otherwise we'd be adding repeatedly)
      if(cmcget) then

        positions=>extract_vector_field(state(1), "Coordinate")
        call allocate(cv_mass, pressure%mesh, "CVMassField")
        call allocate(tempfield, pressure%mesh, "TemporaryAssemblyField")
        call compute_cv_mass(positions, cv_mass)

        allocate(dummy_ones)
        call allocate(dummy_ones, pressure%mesh, "DummyOnesField")
        call set(dummy_ones, 1.0)

        call get_option(trim(pressure_option_path)//'/prognostic/atmospheric_pressure', &
                        atmospheric_pressure, default=0.0)

        call allocate(materialpressure, pressure%mesh, 'MaterialEOSPressure')
        call allocate(materialdrhodp, pressure%mesh, 'DerivativeMaterialdensityWRTBulkPressure')

        call allocate(density, pressure%mesh, 'MaterialDensity')
        call allocate(olddensity, pressure%mesh, 'OldMaterialDensity')
        call allocate(matdrhodpp, pressure%mesh, 'MaterialPressure')
        call allocate(drhodp, pressure%mesh, 'Drhodp')

        density%val = 0.0
        olddensity%val = 0.0
        matdrhodpp%val = 0.0
        drhodp%val=0.0

        do i = 1,size(state)

          materialpressure%val=0.0
          materialdrhodp%val=0.0

          call compressible_material_eos(state(i), materialpressure=materialpressure, materialdrhodp=materialdrhodp)

          volumefraction=>extract_scalar_field(state(i),'MaterialVolumeFraction', stat=stat)
          if(stat==0) then
            oldvolumefraction=>extract_scalar_field(state(i),'OldMaterialVolumeFraction')
            materialdensity=>extract_scalar_field(state(i),'MaterialDensity')
            oldmaterialdensity=>extract_scalar_field(state(i),'OldMaterialDensity')

            density%val = density%val &
                              + materialdensity%val*volumefraction%val
            olddensity%val = olddensity%val &
                                + oldmaterialdensity%val*oldvolumefraction%val
            matdrhodpp%val = matdrhodpp%val &
                                  + materialpressure%val*materialdrhodp%val*volumefraction%val
            drhodp%val = drhodp%val &
                              + materialdrhodp%val*volumefraction%val
          endif

        end do

        call zero(tempfield)
        tempfield%val = (1./(dt*dt))*cv_mass%val*drhodp%val

        call addto_diag(cmc, tempfield)

        rhs%val = (1./dt)*cv_mass%val* &
                          ( &
                            olddensity%val &
                          - density%val &
                          ) &
               +(1./dt)*cv_mass%val* &
                          ( &
                            matdrhodpp%val &
                          - drhodp%val*(pressure%val+atmospheric_pressure) &
                          )

        call deallocate(density)
        call deallocate(olddensity)
        call deallocate(matdrhodpp)
        call deallocate(drhodp)

        call deallocate(materialpressure)
        call deallocate(materialdrhodp)

        call deallocate(cv_mass)
        call deallocate(tempfield)
        call deallocate(dummy_ones)
        deallocate(dummy_ones)

      end if

    end if

  end subroutine assemble_mmat_compressible_projection_cv
  
  subroutine assemble_compressible_projection_cg(state, istate, cmc, rhs, dt, theta_pg, theta_divergence, cmcget)

    ! inputs:
    ! bucket full of fields
    type(state_type), dimension(:), intent(inout) :: state
    integer, intent(in) :: istate

    type(csr_matrix), intent(inout) :: cmc
    type(scalar_field), intent(inout) :: rhs

    real, intent(in) :: dt
    real, intent(in) :: theta_pg, theta_divergence
    logical, intent(in) :: cmcget

    if(option_count("/material_phase/vector_field::Velocity/prognostic") > 1) then
       multiphase = .true.
       call assemble_1mat_compressible_projection_cg(state(istate), cmc, rhs, dt, &
                                                      theta_pg, theta_divergence, cmcget)
    else
       multiphase = .false.

       if((size(state)==1).and.(.not.has_scalar_field(state(1), "MaterialVolumeFraction"))) then
      
          call assemble_1mat_compressible_projection_cg(state(1), cmc, rhs, dt, &
                                                      theta_pg, theta_divergence, cmcget)
         
       else
         
          FLExit("Multimaterial compressible continuous_galerkin pressure not possible.")
         
       end if

   end if    

  end subroutine assemble_compressible_projection_cg

  subroutine assemble_1mat_compressible_projection_cg(state, cmc, rhs, dt, &
                                                      theta_pg, theta_divergence, cmcget)

    ! inputs:
    ! bucket full of fields
    type(state_type), intent(inout) :: state

    type(csr_matrix), intent(inout) :: cmc
    type(scalar_field), intent(inout) :: rhs

    real, intent(in) :: dt
    real, intent(in) :: theta_pg, theta_divergence
    logical, intent(in) :: cmcget

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

    ! loop integers
    integer :: ele

    ! pointer to coordinates
    type(vector_field), pointer :: coordinate, nonlinearvelocity, velocity
    type(scalar_field), pointer :: pressure, density, olddensity
    type(scalar_field), pointer :: source, absorption
    type(scalar_field) :: eospressure, drhodp
    real :: theta, atmospheric_pressure

    real, dimension(:,:), allocatable :: ele_mat
    
    logical :: have_absorption, have_source
    integer :: stat

    !! Multiphase variables
    ! Volume fraction fields
    type(scalar_field), pointer :: vfrac
    type(scalar_field) :: nvfrac

    ! =============================================================
    ! Subroutine to construct the matrix CT_m (a.k.a. C1/2/3T).
    ! =============================================================

    ewrite(1,*) 'Entering assemble_1mat_compressible_projection_cg'
    
    call zero(rhs)

    ! only do all this if we need to make cmc (otherwise we'd be adding repeatedly)
    if(cmcget) then
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

      ! Get the non-linear PhaseVolumeFraction field if multiphase
      if(multiphase) then
         vfrac => extract_scalar_field(state, "PhaseVolumeFraction")
         call allocate(nvfrac, vfrac%mesh, "NonlinearPhaseVolumeFraction")
         call zero(nvfrac)
         call get_nonlinear_volume_fraction(state, nvfrac)
         ewrite_minmax(nvfrac)
      end if
      
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
  
        if(multiphase) then
           ele_mat = (1./(dt*dt*theta_divergence*theta_pg))*shape_shape(test_shape, test_shape_ptr, detwei*ele_val_at_quad(nvfrac, ele)*drhodp_at_quad)
        else
           ele_mat = (1./(dt*dt*theta_divergence*theta_pg))*shape_shape(test_shape, test_shape_ptr, detwei*drhodp_at_quad)
        end if
        !       /
        ! rhs = |test_shape* &
        !       /
        !      ((1./dt)*(drhodp*(eospressure - (pressure + atmospheric_pressure)) + olddensity - density)
        ! +(absorption)*(drhodp*theta*(eospressure - (pressure + atmospheric_pressure)) 
        !                - theta*density - (1-theta)*olddensity)
        ! +source)dV

        if(multiphase) then
            ele_rhs = (1./dt)*shape_rhs(test_shape, detwei*(ele_val_at_quad(nvfrac, ele))*((drhodp_at_quad*(eosp_at_quad - p_at_quad)) &
                                                       +(olddensity_at_quad - density_at_quad)))
        else
            ele_rhs = (1./dt)*shape_rhs(test_shape, detwei*((drhodp_at_quad*(eosp_at_quad - p_at_quad)) &
                                                       +(olddensity_at_quad - density_at_quad)))
        end if
        
        if(have_source) then
          ele_rhs = ele_rhs + shape_rhs(test_shape, detwei*ele_val_at_quad(source, ele))
        end if
        
        if(have_absorption) then
          abs_at_quad = ele_val_at_quad(absorption, ele)
          ele_mat = ele_mat + &
                    (theta/(dt*theta_divergence*theta_pg))*shape_shape(test_shape, test_shape_ptr, &
                                                                       detwei*drhodp_at_quad*abs_at_quad)
          ele_rhs = ele_rhs + &
                    shape_rhs(test_shape, detwei*abs_at_quad*(theta*(drhodp_at_quad*(eosp_at_quad - p_at_quad)-density_at_quad) &
                                                             -(1-theta)*olddensity_at_quad))
        end if
        
        call addto(cmc, test_nodes, test_nodes, ele_mat)
        
        call addto(rhs, test_nodes, ele_rhs)
        
        call deallocate(test_shape)
        
      end do
  
      call deallocate(drhodp)
      call deallocate(eospressure)

      if(multiphase) then
         call deallocate(nvfrac)
      end if
    
    end if

  end subroutine assemble_1mat_compressible_projection_cg

  subroutine update_compressible_density(state)
  
    type(state_type), dimension(:), intent(inout) :: state
    
    type(scalar_field), pointer :: density

    integer :: istate
    
    if(option_count("/material_phase/vector_field::Velocity/prognostic") > 1) then
       do istate=1,size(state)
          density=>extract_scalar_field(state(istate),'Density')
          
          if(have_option(trim(density%option_path)//"/prognostic")) then
            call compressible_eos(state(istate), density=density)
          end if
       end do
    else
       if((size(state)==1).and.(.not.has_scalar_field(state(1), "MaterialVolumeFraction"))) then
       
          density=>extract_scalar_field(state(1),'Density')
          
          if(have_option(trim(density%option_path)//"/prognostic")) then
         
             call compressible_eos(state(1), density=density)
         
          end if
      
       end if
    end if
  
  end subroutine update_compressible_density

  subroutine compressible_projection_check_options

    character(len=OPTION_PATH_LEN):: pressure_option_path
    integer:: iphase
    logical:: have_compressible_eos

    do iphase=0, option_count("/material_phase")-1
       have_compressible_eos = have_option("/material_phase["//int2str(iphase)//"]/equation_of_state/compressible")
       pressure_option_path = "/material_phase["//int2str(iphase)//"]/scalar_field::Pressure"
       if(have_compressible_eos.and. &
            have_option(trim(pressure_option_path)//"/prognostic/spatial_discretisation/discontinuous_galerkin")) then
          FLExit("With a DG pressure you cannot have use a compressible eos")
       end if
    end do

  end subroutine compressible_projection_check_options

end module compressible_projection

