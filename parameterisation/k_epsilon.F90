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
!    but WITHOUT ANY WARRANTY; without even the implied arranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module k_epsilon
  use quadrature
  use elements
  use field_derivatives
  use fields
  use field_options
  use state_module
  use spud
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use state_fields_module
  use boundary_conditions
  use fields_manipulation
  use surface_integrals
  use fetools
  use vector_tools
  use sparsity_patterns_meshes
  use FLDebug

implicit none

  private

  ! old turbulent kinetic energy, lengthscale, field ratios
  type(scalar_field), save            :: tke_old, ll, tkeovereps, epsovertke
  ! empirical constants
  real, save                          :: c_mu, c_eps_1, c_eps_2, sigma_eps, sigma_k
  real, save                          :: fields_min = 1.e-10
  character(len=FIELD_NAME_LEN), save :: src_abs

  public :: keps_init, keps_cleanup, keps_tke, keps_eps, keps_eddyvisc, keps_bcs, keps_adapt_mesh, keps_check_options

  ! Outline:
  !  - Init in populate_state.
  !  - call keps_tke (which calculates production and sets source/absorption/diffusivity for solve).
  !  - call keps_eps (which sets source/absorption/diffusivity for solve).
  !  - After keps_eps solve, keps_eddyvisc recalculates the eddy viscosity and adds it to the viscosity field.
  !  - Wall functions are added to selected boundaries in keps_bcs and wall_functions.
  !  - keps_adapt_options repopulates the fields after an adapt.
  !  - When done, clean-up.

contains

!----------
! initialise parameters based on options
!----------

subroutine keps_init(state)

    type(state_type), intent(inout) :: state
    type(scalar_field), pointer     :: Field
    character(len=OPTION_PATH_LEN)  :: keps_path
    real                            :: visc
    integer                         :: i

    ewrite(1,*)'Now in k_epsilon turbulence model - keps_init'
    keps_path = trim(state%option_path)//"/subgridscale_parameterisations/k-epsilon"
    ewrite(2,*)'keps_path: ', trim(keps_path)

    ! Allocate the temporary, module-level variables
    call keps_allocate_fields(state)

    ! Are source and absorption terms implicit/explicit?
    call get_option(trim(keps_path)//'/source_absorption', src_abs)

    ! Get the 5 model constants
    call get_option(trim(keps_path)//'/C_mu', C_mu, default = 0.09)
    call get_option(trim(keps_path)//'/C_eps_1', c_eps_1, default = 1.44)
    call get_option(trim(keps_path)//'n/C_eps_2', c_eps_2, default = 1.92)
    call get_option(trim(keps_path)//'/sigma_k', sigma_k, default = 1.0)
    call get_option(trim(keps_path)//'n/sigma_eps', sigma_eps, default = 1.3)

    ! Get background viscosity
    call get_option(trim(keps_path)//"/tensor_field::BackgroundViscosity/prescribed/&
                          &value::WholeMesh/isotropic/constant", visc)
    
    ! initialise eddy viscosity field
    Field => extract_scalar_field(state, "ScalarEddyViscosity")
    call set(Field, visc)

    ewrite(2,*) "k-epsilon parameters"
    ewrite(2,*) "--------------------------------------------"
    ewrite(2,*) "c_mu: ",     c_mu
    ewrite(2,*) "c_eps_1: ",  c_eps_1
    ewrite(2,*) "c_eps_2: ",  c_eps_2
    ewrite(2,*) "sigma_k: ",  sigma_k
    ewrite(2,*) "sigma_eps: ",sigma_eps
    ewrite(2,*) "fields_min: ",fields_min
    ewrite(2,*) "background visc: ",   visc
    ewrite(2,*) "implicit/explicit source/absorption terms: ",   trim(src_abs)
    ewrite(2,*) "--------------------------------------------"

end subroutine keps_init

!------------------------------------------------------------------------------

subroutine keps_tke(state)

    type(state_type), intent(inout)    :: state
    type(scalar_field), pointer        :: source_kk, absorption_kk, kk, eps, EV, lumped_mass
    type(scalar_field)                 :: src_rhs, abs_rhs
    type(vector_field), pointer        :: positions, nu, u
    type(tensor_field), pointer        :: kk_diff
    type(element_type), pointer        :: shape_kk
    integer                            :: i, ele, stat
    integer, pointer, dimension(:)     :: nodes_kk
    real                               :: residual
    real, allocatable, dimension(:)    :: detwei, strain_ngi, rhs_addto
    real, allocatable, dimension(:,:,:):: dshape_kk

    ewrite(1,*) "In keps_tke"

    positions       => extract_vector_field(state, "Coordinate")
    nu              => extract_vector_field(state, "NonlinearVelocity")
    u               => extract_vector_field(state, "Velocity")
    source_kk       => extract_scalar_field(state, "TurbulentKineticEnergySource")
    absorption_kk   => extract_scalar_field(state, "TurbulentKineticEnergyAbsorption")
    kk_diff         => extract_tensor_field(state, "TurbulentKineticEnergyDiffusivity")
    kk              => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps             => extract_scalar_field(state, "TurbulentDissipation")
    EV              => extract_scalar_field(state, "ScalarEddyViscosity")

    ! Set copy of old kk for eps solve
    call set(tke_old, kk)

    call allocate(src_rhs, kk%mesh, name="KKSRCRHS")
    call allocate(abs_rhs, kk%mesh, name="KKABSRHS")
    call zero(src_rhs); call zero(abs_rhs)

    ! Assembly loop
    do ele = 1, ele_count(kk)
        shape_kk => ele_shape(kk, ele)
        nodes_kk => ele_nodes(kk, ele)

        allocate(dshape_kk (size(nodes_kk), ele_ngi(kk, ele), positions%dim))
        allocate(detwei (ele_ngi(kk, ele)))
        allocate(strain_ngi (ele_ngi(kk, ele)))
        allocate(rhs_addto(ele_loc(kk, ele)))
        call transform_to_physical( positions, ele, shape_kk, dshape=dshape_kk, detwei=detwei )

        ! Calculate TKE production at ngi using strain rate (double_dot_product) function
        strain_ngi = double_dot_product(dshape_kk, ele_val(u, ele) )

        ! Source term:
        rhs_addto = shape_rhs(shape_kk, detwei*strain_ngi*ele_val_at_quad(EV, ele))
        call addto(src_rhs, nodes_kk, rhs_addto)

        ! Absorption term:
        rhs_addto = shape_rhs(shape_kk, detwei*ele_val_at_quad(eps,ele)/ele_val_at_quad(kk,ele))
        call addto(abs_rhs, nodes_kk, rhs_addto)

        deallocate(dshape_kk, detwei, strain_ngi, rhs_addto)
    end do

    lumped_mass => get_lumped_mass(state, kk%mesh)

    ! This allows user-specified source term, so that an MMS test can be set up.
    if(have_option(trim(source_kk%option_path)//"/diagnostic/internal")) then
      ewrite(2,*) "Calculating k source and absorption"
      do i = 1, node_count(kk)
        select case (src_abs)
        case ("explicit")
          call set(source_kk, i, node_val(src_rhs,i)/node_val(lumped_mass,i))
          call set(absorption_kk, i, node_val(abs_rhs,i)/node_val(lumped_mass,i))
        case ("implicit")
          residual = (node_val(abs_rhs,i) - node_val(src_rhs,i))/node_val(lumped_mass,i)
          call set(source_kk, i, -min(0.0, residual) )
          call set(absorption_kk, i, max(0.0, residual) )
        case default
          FLAbort("Invalid implicitness option for k")
        end select
      end do
    else if(have_option(trim(source_kk%option_path)//"/prescribed")) then
      ewrite(2,*) "Prescribed k source"
      do i = 1, node_count(kk)
        select case (src_abs)
        case ("explicit")
          call set(absorption_kk, i, node_val(abs_rhs,i)/node_val(lumped_mass,i))
        case ("implicit")
          residual = node_val(abs_rhs,i)/node_val(lumped_mass,i) - node_val(source_kk,i)
          call set(source_kk, i, -min(0.0, residual) )
          call set(absorption_kk, i, max(0.0, residual) )
        case default
          FLAbort("Invalid implicitness option for k")
        end select
      end do
    end if

    call deallocate(src_rhs); call deallocate(abs_rhs)

    ! Set diffusivity for k equation.
    call zero(kk_diff)
    do i = 1, kk_diff%dim(1)
        call set(kk_diff, i, i, EV, scale=1. / sigma_k)
    end do

    ewrite_minmax(kk_diff)
    ewrite_minmax(tke_old)
    ewrite_minmax(kk)
    ewrite_minmax(source_kk)
    ewrite_minmax(absorption_kk)

end subroutine keps_tke

!----------------------------------------------------------------------------------

subroutine keps_eps(state)

    type(state_type), intent(inout)    :: state
    type(scalar_field), pointer        :: source_eps, source_kk, absorption_eps, eps, EV, lumped_mass
    type(scalar_field)                 :: src_rhs, abs_rhs
    type(vector_field), pointer        :: positions
    type(tensor_field), pointer        :: eps_diff
    type(element_type), pointer        :: shape_eps
    integer                            :: i, ele
    real                               :: residual
    real, allocatable, dimension(:)    :: detwei, rhs_addto
    integer, pointer, dimension(:)     :: nodes_eps

    ewrite(1,*) "In keps_eps"
    eps             => extract_scalar_field(state, "TurbulentDissipation")
    source_eps      => extract_scalar_field(state, "TurbulentDissipationSource")
    source_kk       => extract_scalar_field(state, "TurbulentKineticEnergySource")
    absorption_eps  => extract_scalar_field(state, "TurbulentDissipationAbsorption")
    eps_diff        => extract_tensor_field(state, "TurbulentDissipationDiffusivity")
    EV              => extract_scalar_field(state, "ScalarEddyViscosity")
    positions       => extract_vector_field(state, "Coordinate")

    call allocate(src_rhs, eps%mesh, name="EPSSRCRHS")
    call allocate(abs_rhs, eps%mesh, name="EPSABSRHS")
    call zero(src_rhs); call zero(abs_rhs)

    ! Assembly loop
    do ele = 1, ele_count(eps)
        shape_eps => ele_shape(eps, ele)
        nodes_eps => ele_nodes(eps, ele)

        allocate(detwei (ele_ngi(eps, ele)))
        allocate(rhs_addto(ele_loc(eps, ele)))
        call transform_to_physical(positions, ele, detwei=detwei)

        ! Source term:
        rhs_addto = shape_rhs(shape_eps, detwei*c_eps_1*ele_val_at_quad(eps,ele)/ &
                              ele_val_at_quad(tke_old,ele)*ele_val_at_quad(source_kk,ele))
        call addto(src_rhs, nodes_eps, rhs_addto)

        ! Absorption term:
        rhs_addto = shape_rhs(shape_eps, detwei*c_eps_2*ele_val_at_quad(eps,ele)/ &
                              ele_val_at_quad(tke_old,ele))

        call addto(abs_rhs, nodes_eps, rhs_addto)

        deallocate(detwei, rhs_addto)
    end do

    lumped_mass => get_lumped_mass(state, eps%mesh)

    ! This allows user-specified source term, so that an MMS test can be set up.
    if(have_option(trim(source_eps%option_path)//"/diagnostic/internal")) then
      ewrite(2,*) "Calculating epsilon source and absorption"
      do i = 1, node_count(eps)
        select case (src_abs)
        case ("explicit")
          call set(source_eps, i, node_val(src_rhs,i)/node_val(lumped_mass,i))
          call set(absorption_eps, i, node_val(abs_rhs,i)/node_val(lumped_mass,i))
        case ("implicit")
          residual = (node_val(abs_rhs,i) - node_val(src_rhs,i))/node_val(lumped_mass,i)
          call set(source_eps, i, -min(0.0, residual) )
          call set(absorption_eps, i, max(0.0, residual) )
        case default
          FLAbort("Invalid implicitness option for epsilon")
        end select
      end do
    else if(have_option(trim(source_kk%option_path)//"/prescribed")) then
      ewrite(2,*) "Prescribed epsilon source"
      do i = 1, node_count(eps)
        select case (src_abs)
        case ("explicit")
          call set(absorption_eps, i, node_val(abs_rhs,i)/node_val(lumped_mass,i))
        case ("implicit")
          residual = node_val(abs_rhs,i)/node_val(lumped_mass,i) - node_val(source_eps,i)
          call set(source_eps, i, -min(0.0, residual) )
          call set(absorption_eps, i, max(0.0, residual) )
        case default
          FLAbort("Invalid implicitness option for k")
        end select
      end do
    end if

    call deallocate(src_rhs); call deallocate(abs_rhs)

    ! Set diffusivity for Eps
    call zero(eps_diff)
    do i = 1, eps_diff%dim(1)
        call set(eps_diff, i,i, EV, scale=1./sigma_eps)
    end do

    ewrite_minmax(eps_diff)
    ewrite_minmax(tke_old)
    ewrite_minmax(source_kk)
    ewrite_minmax(eps)
    ewrite_minmax(source_eps)
    ewrite_minmax(absorption_eps)

end subroutine keps_eps

!----------
! eddyvisc calculates the lengthscale, and then the eddy viscosity.
! Eddy viscosity is added to the background viscosity.
!----------

subroutine keps_eddyvisc(state)

    type(state_type), intent(inout)  :: state
    type(tensor_field), pointer      :: eddy_visc, viscosity, diffusivity, bg_visc
    type(vector_field), pointer      :: positions
    type(scalar_field), pointer      :: kk, eps, EV, scalarField, lumped_mass
    type(scalar_field)               :: ev_rhs
    type(element_type), pointer      :: shape_ev
    integer                          :: i, ele, stat
    integer, pointer, dimension(:)   :: nodes_ev
    real, allocatable, dimension(:)  :: detwei, rhs_addto

    ewrite(1,*) "In keps_eddyvisc"
    kk         => extract_scalar_field(state, "TurbulentKineticEnergy")
    eps        => extract_scalar_field(state, "TurbulentDissipation")
    positions  => extract_vector_field(state, "Coordinate")
    eddy_visc  => extract_tensor_field(state, "KEpsEddyViscosity")
    viscosity  => extract_tensor_field(state, "Viscosity")
    bg_visc    => extract_tensor_field(state, "BackgroundViscosity")
    EV         => extract_scalar_field(state, "ScalarEddyViscosity")

    ewrite_minmax(kk)
    ewrite_minmax(eps)
    ewrite_minmax(EV)

    call allocate(ev_rhs, EV%mesh, name="EVRHS")
    call zero(ev_rhs)

    ! Initialise viscosity to background value
    call set(viscosity, bg_visc)

    !Clip fields: can't allow negative/zero epsilon or k
    do i = 1, node_count(EV)
      call set(kk, i, max(node_val(kk,i), fields_min))
      call set(eps, i, max(node_val(eps,i), fields_min))
    end do

    ! Calculate scalar eddy viscosity by integration over element
    do ele = 1, ele_count(EV)
      nodes_ev => ele_nodes(EV, ele)
      shape_ev =>  ele_shape(EV, ele)
      allocate(detwei (ele_ngi(EV, ele)))
      allocate(rhs_addto (ele_loc(EV, ele)))
      call transform_to_physical(positions, ele, detwei=detwei)

      rhs_addto = shape_rhs(shape_ev, detwei*C_mu*ele_val_at_quad(kk,ele)**2./ele_val_at_quad(eps,ele))
      call addto(ev_rhs, nodes_ev, rhs_addto)

      deallocate(detwei, rhs_addto)
    end do

    lumped_mass => get_lumped_mass(state, EV%mesh)

    ! Node loop: set eddy viscosity at nodes
    do i = 1, node_count(EV)
      call set(EV, i, node_val(ev_rhs,i)/node_val(lumped_mass,i))
      ! Now set diagnostic fields. These do not need to be assembled by integration
      ! because we do not need to know the values at gauss points.
      call set(ll, i, c_mu * node_val(kk,i)**1.5 / node_val(eps,i))
      call set(tkeovereps, i, node_val(kk,i) / node_val(eps,i))
      call set(epsovertke, i, 1. / node_val(tkeovereps,i))
    end do

    call deallocate(ev_rhs)

    ewrite(2,*) "Set k-epsilon eddy-viscosity tensor"
    call zero(eddy_visc)

    ! eddy tensors are isotropic
    do i = 1, eddy_visc%dim(1)
      call set(eddy_visc, i, i, EV)
    end do
    ewrite_minmax(eddy_visc)

    ! Add turbulence model contribution to viscosity field
    call addto(viscosity, eddy_visc)

    ! Check components of viscosity
    ewrite_minmax(viscosity)

    ! Set output on optional fields
    scalarField => extract_scalar_field(state, "LengthScale", stat)
    if(stat == 0) then
        call set(scalarField, ll) 
    end if
    scalarField => extract_scalar_field(state, "TKEOverEpsilon", stat)
    if(stat == 0) then
        call set(scalarField, tkeovereps)  
    end if
    scalarField => extract_scalar_field(state, "EpsilonOverTKE", stat)
    if(stat == 0) then
        call set(scalarField, epsovertke)  
    end if

    ewrite_minmax(kk)
    ewrite_minmax(eps)
    ewrite_minmax(EV)

end subroutine keps_eddyvisc

!--------------------------------------------------------------------------------!
! This gets and applies locally defined boundary conditions (wall functions)     !
!--------------------------------------------------------------------------------!

subroutine keps_bcs(state)

    type(state_type), intent(in)               :: state
    type(scalar_field), pointer                :: field1, field2    ! k or epsilon
    type(scalar_field), pointer                :: surface_field, EV
    type(vector_field), pointer                :: positions, u
    type(tensor_field), pointer                :: bg_visc
    type(scalar_field)                         :: rhs_field, surface_values
    type(mesh_type), pointer                   :: surface_mesh
    integer                                    :: i, j, ele, sele, index, nbcs
    integer, dimension(:), pointer             :: surface_elements, surface_node_list
    character(len=FIELD_NAME_LEN)              :: bc_type, bc_name, wall_fns
    character(len=OPTION_PATH_LEN)             :: bc_path, bc_path_i
    integer, allocatable, dimension(:)         :: nodes_bdy
    real, dimension(:), allocatable            :: rhs
    real                                       :: cmu
    ewrite(2,*) "In keps_bcs"

    positions => extract_vector_field(state, "Coordinate")
    u         => extract_vector_field(state, "Velocity")
    EV        => extract_scalar_field(state, "ScalarEddyViscosity")
    bg_visc   => extract_tensor_field(state, "BackgroundViscosity")

    ! THIS IS NOT AVAILABLE BEFORE KEPS_INIT HAS BEEN CALLED!
    call get_option(trim(state%option_path)//"/subgridscale_parameterisations/k-epsilon/C_mu", cmu, default = 0.09)

    field_loop: do index=1,2

      if(index==1) then
        field1 => extract_scalar_field(state, "TurbulentKineticEnergy")
        field2 => null()
      else
        field1 => extract_scalar_field(state, "TurbulentDissipation")
        field2 => extract_scalar_field(state, "TurbulentKineticEnergy")
      end if

      bc_path=trim(field1%option_path)//'/prognostic/boundary_conditions'
      nbcs=option_count(trim(bc_path))

      ! Loop over boundary conditions for field1
      boundary_conditions: do i=0, nbcs-1

        bc_path_i=trim(bc_path)//"["//int2str(i)//"]"

        ! Get name and type of boundary condition
        call get_option(trim(bc_path_i)//"/name", bc_name)
        call get_option(trim(bc_path_i)//"/type[0]/name", bc_type)

        ewrite(2,*) "Checking field BC: ",trim(field1%name),' ',trim(bc_name),' ',trim(bc_type)

        if (trim(bc_type) .eq. "k_epsilon") then
          ! Get bc by name. Get type just to make sure it's now dirichlet
          call get_boundary_condition(field1, name=bc_name, type=bc_type, surface_node_list=surface_node_list, &
                                      surface_element_list=surface_elements, surface_mesh=surface_mesh)
          !ewrite(3,*) "surface_node_list: ", surface_node_list

          ! Do we have high- or low-Reynolds-number wall functions?
          call get_option(trim(bc_path_i)//"/type::k_epsilon/", wall_fns)

          ewrite(2,*) "Calculating field BC: ",trim(field1%name),' ',trim(bc_name),' ',trim(bc_type),' ',trim(wall_fns)
          !ewrite(3,*) "surface_mesh%: ", trim(surface_mesh%name), surface_mesh%elements, surface_mesh%nodes
          !ewrite(3,*) "---------------------------------------"

          ! Get surface field already created in bcs_from_options
          surface_field => extract_surface_field(field1, bc_name=bc_name, name="value")
          call zero(surface_field)

          call allocate(surface_values, surface_mesh, name="surfacevalues")
          call allocate(rhs_field, field1%mesh, name="rhs")
          call zero(surface_values); call zero(rhs_field)
          !ewrite(3,*) "rhs_field%: ", trim(rhs_field%name), rhs_field%mesh%elements, rhs_field%mesh%nodes

          ewrite(3,*) "Entering surface element loop"
          do j = 1, ele_count(surface_mesh)
             ! I want ele and sele on volume mesh for keps_wall_function
             sele = surface_elements(j)
             !sele = j
             ele  = face_ele(rhs_field, sele)
             !ewrite(3,*) "j, sele, ele: ", j, sele, ele
             !ewrite(3,*) "ele_nodes: ", ele_nodes(rhs_field,ele)
             !ewrite(3,*) "sele_nodes: ", ele_nodes(rhs_field,sele)

             allocate(rhs(face_loc(rhs_field, sele)))
             allocate(nodes_bdy(face_loc(rhs_field, sele)))

             nodes_bdy = face_global_nodes(rhs_field, sele)
             !nodes_bdy = face_global_nodes(rhs_field, j)
             !nodes_bdy = ele_nodes(surface_field, j)

             !ewrite(3,*) "nodes_bdy: ", nodes_bdy
             !ewrite(3,*) "node_val rhs_field: ", node_val(rhs_field,nodes_bdy)

             ! Calculate wall function
             call keps_wall_function(field1,field2,positions,u,bg_visc,EV,ele,sele,index,wall_fns,cmu,rhs)
             ewrite(3,*) "rhs: ", rhs
             ewrite(3,*) "----------------------------------------"
             ! Add element contribution to rhs field
             call addto(rhs_field, nodes_bdy, rhs)
             
             deallocate(rhs); deallocate(nodes_bdy)
          end do

          ! Put values onto surface mesh
          call remap_field_to_surface(rhs_field, surface_values, surface_elements)

          do j = 1, size(surface_node_list)
             !node = surface_node_list(j)
             ewrite(3,*) "node, val: ", j, node_val(surface_values, j)
             call set(surface_field, j, node_val(surface_values, j))
          end do

          call deallocate(surface_values); call deallocate(rhs_field)

        end if
      end do boundary_conditions
    end do field_loop

end subroutine keps_bcs

!------------------------------------------------------------------------------------

subroutine keps_cleanup()

    ewrite(1,*) "In keps_cleanup"

    call deallocate(ll)
    call deallocate(tke_old)
    call deallocate(tkeovereps)
    call deallocate(epsovertke)

end subroutine keps_cleanup

!------------------------------------------------------------------------------------!
! Called after an adapt to reset the fields and arrays within the module             !
!------------------------------------------------------------------------------------!

subroutine keps_adapt_mesh(state)

    type(state_type), intent(inout) :: state

    ewrite(1,*) "In keps_adapt_mesh"

    call keps_allocate_fields(state)
    call keps_eddyvisc(state)
    call keps_tke(state)
    call keps_eps(state)
    call keps_eddyvisc(state)

end subroutine keps_adapt_mesh

!---------------------------------------------------------------------------------

subroutine keps_check_options(state)

    type(state_type), intent(in)   :: state
    type(vector_field), pointer    :: u
    character(len=OPTION_PATH_LEN) :: option_path
    character(len=FIELD_NAME_LEN)  :: kmsh, emsh, vmsh
    integer                        :: dimension

    ewrite(1,*) "In keps_check_options"
    option_path = trim(state%option_path)//"/subgridscale_parameterisations/k-epsilon"
    u => extract_vector_field(state, "Velocity")

    ! one dimensional problems not supported
    call get_option("/geometry/dimension/", dimension) 
    if (dimension .eq. 1 .and. have_option(trim(option_path))) then
        FLExit("k-epsilon model is only supported for dimension > 1")
    end if
    ! Don't do k-epsilon if it's not included in the model!
    if (.not.have_option(trim(option_path))) return

    ! checking for required fields
    if (.not.have_option(trim(option_path)//"/scalar_field::TurbulentKineticEnergy")) then
        FLExit("You need TurbulentKineticEnergy field for k-epsilon")
    end if
    if (.not.have_option(trim(option_path)//"/scalar_field::TurbulentDissipation")) then
        FLExit("You need TurbulentDissipation field for k-epsilon")
    end if
    ! Check that TurbulentKineticEnergy and TurbulentDissipation fields are on the same mesh
    call get_option(trim(option_path)//"/scalar_field::TurbulentKineticEnergy/prognostic/mesh", kmsh)
    call get_option(trim(option_path)//"/scalar_field::TurbulentDissipation/prognostic/mesh", emsh)
    call get_option(trim(state%option_path)//"/vector_field::Velocity/prognostic/mesh", vmsh)
    if(.not. kmsh==emsh .or. .not. kmsh==vmsh .or. .not. emsh==vmsh) then
        FLExit("You must use the Velocity mesh for TurbulentKineticEnergy and TurbulentDissipation fields")
    end if
    ! check that diffusivity is on for the two turbulent fields, and diagnostic
    if (.not.have_option(trim(option_path)//"/scalar_field::TurbulentKineticEnergy/&
        &prognostic/tensor_field::Diffusivity")) then
        FLExit("You need TurbulentKineticEnergy Diffusivity field for k-epsilon")
    end if    
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &tensor_field::Diffusivity/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentKineticEnergy Diffusivity field set to diagnostic/internal")
    end if
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &tensor_field::Diffusivity")) then
        FLExit("You need TurbulentDissipation Diffusivity field for k-epsilon")
    end if
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &tensor_field::Diffusivity/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentDissipation Diffusivity field set to diagnostic/internal")
    end if
    ! source terms
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &scalar_field::Source")) then
        FLExit("You need TurbulentKineticEnergy Source field for k-epsilon")
    end if    
    if (.not. have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &scalar_field::Source/diagnostic/algorithm::Internal")&
        .or. .not. have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &scalar_field::Source/prescribed")) then
        FLExit("You need TurbulentKineticEnergy Source field set to diagnostic/internal or prescribed")
    end if
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &scalar_field::Source")) then
        FLExit("You need TurbulentDissipation Source field for k-epsilon")
    end if
    if (.not. have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &scalar_field::Source/diagnostic/algorithm::Internal")&
        .or. .not. have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &scalar_field::Source/prescribed")) then
        FLExit("You need TurbulentDissipation Source field set to diagnostic/internal or prescribed")
    end if
    ! absorption terms
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &scalar_field::Absorption")) then
        FLExit("You need TurbulentKineticEnergy Absorption field for k-epsilon")
    end if    
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentKineticEnergy/prognostic/&
                          &scalar_field::Absorption/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentKineticEnergy Absorption field set to diagnostic/internal")
    end if
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &scalar_field::Absorption")) then
        FLExit("You need TurbulentDissipation Absorption field for k-epsilon")
    end if
    if (.not.have_option(trim(option_path)//"/&
                          &scalar_field::TurbulentDissipation/prognostic/&
                          &scalar_field::Absorption/diagnostic/algorithm::Internal")) then
        FLExit("You need TurbulentDissipation Absorption field set to diagnostic/internal")
    end if
    ! check there's a viscosity somewhere
    if (.not.have_option(trim(state%option_path)//"/vector_field::Velocity/prognostic/&
                          &tensor_field::Viscosity/")) then
        FLExit("Need viscosity switched on under the Velocity field for k-epsilon.") 
    end if
    ! check that the user has switched Velocity/viscosity to diagnostic
    if (.not.have_option(trim(state%option_path)//"/vector_field::Velocity/prognostic/&
                          &tensor_field::Viscosity/diagnostic/")) then
        FLExit("You need to switch the viscosity field under Velocity to diagnostic/internal")
    end if
    ! check that the background viscosity is an isotropic constant
    if (.not.have_option(trim(option_path)//"/tensor_field::BackgroundViscosity/prescribed/&
                          &value::WholeMesh/isotropic/constant")) then
        FLExit("You need to switch the BackgroundViscosity field to isotropic/constant")
    end if

end subroutine keps_check_options

!------------------------------------------------------------------!
!                       Private subroutines                        !
!------------------------------------------------------------------!

subroutine keps_allocate_fields(state)

    type(state_type), intent(inout) :: state
    type(vector_field), pointer     :: vectorField

    ewrite(1,*) "In keps_allocate_fields"

    vectorField => extract_vector_field(state, "Velocity")

    ! allocate some space for the fields we need for calculations
    call allocate(ll,         vectorField%mesh, "LengthScale")
    call allocate(tke_old,    vectorField%mesh, "Old_TKE")
    call allocate(tkeovereps, vectorField%mesh, "TKEoverEpsilon")
    call allocate(epsovertke, vectorField%mesh, "EpsilonOverTKE")

end subroutine keps_allocate_fields

!--------------------------------------------------------------------------------!
! Only used if bc type == k_epsilon for field.                                    !
!--------------------------------------------------------------------------------!

subroutine keps_wall_function(field1,field2,positions,u,bg_visc,EV,ele,sele,index,wall_fns,cmu,rhs)

    type(scalar_field), pointer, intent(in)              :: field1, field2, EV
    type(vector_field), pointer, intent(in)              :: positions, u
    type(tensor_field), pointer, intent(in)              :: bg_visc
    integer, intent(in)                                  :: ele, sele, index
    character(len=FIELD_NAME_LEN), intent(in)            :: wall_fns
    real, dimension(face_loc(field1,sele)), intent(inout):: rhs

    type(element_type), pointer                          :: shape, fshape
    integer                                              :: i, j, gi, sgi, sloc
    real                                                 :: kappa, h, cmu
    real, dimension(1,1)                                 :: hb
    real, dimension(ele_ngi(field1,ele))                 :: detwei
    real, dimension(face_ngi(field1,sele))               :: detwei_bdy, ustar, q_sgin, visc_sgi
    real, dimension(face_loc(field1,sele))               :: lumpedfmass
    real, dimension(ele_loc(field1,ele))                 :: sqrt_k
    real, dimension(positions%dim,1)                     :: n
    real, dimension(positions%dim,positions%dim)         :: G
    real, dimension(positions%dim,ele_ngi(field1,ele))   :: grad_k
    real, dimension(positions%dim,face_ngi(field1,sele)) :: normal_bdy, q_sgi, qq_sgin
    real, dimension(positions%dim,ele_loc(field1,ele))   :: q
    real, dimension(positions%dim,face_loc(field1,ele))  :: q_s
    real, dimension(face_loc(field1,ele),face_loc(field1,ele))   :: fmass
    real, dimension(ele_loc(field1,ele),ele_loc(field1,ele))     :: invmass
    real, dimension(positions%dim,positions%dim,ele_loc(field1,sele))       :: qq
    real, dimension(positions%dim,positions%dim,ele_ngi(field1,sele))       :: grad_u, invJ
    real, dimension(positions%dim,positions%dim,face_loc(field1,sele))      :: qq_s
    real, dimension(positions%dim,positions%dim,face_ngi(field1,sele))      :: bg_visc_sgi, qq_sgi
    real, dimension(ele_loc(field1,ele),ele_ngi(field1,ele),positions%dim)  :: dshape

    ! Get ids, lists and shape functions
    sgi      =  face_ngi(field1, sele)    ! no. of gauss points in surface element
    sloc     =  face_loc(field1, sele)    ! no. of nodes on surface element
    shape  => ele_shape(field1, ele)    ! scalar field shape functions in volume element
    fshape => face_shape(field1, sele)  ! scalar field shape functions in surface element

    ! Get shape fn gradients, element/face quadrature weights, and surface normal
    call transform_to_physical( positions, ele, shape, dshape=dshape, detwei=detwei, invJ=invJ )
    call transform_facet_to_physical( positions, sele, detwei_f=detwei_bdy, normal=normal_bdy )

    invmass = shape_shape(shape, shape, detwei)
    call invert(invmass)

    fmass = shape_shape(fshape, fshape, detwei_bdy)  !*density_gi*vfrac_gi to be safe?
    lumpedfmass = sum(fmass, 2)

    ! low Re wall functions for k and epsilon: see e.g. Wilcox (1994)
    if(wall_fns=="low_Re") then
       if (index==1) then
          rhs = 0.0
       else if (index==2) then
          bg_visc_sgi = face_val_at_quad(bg_visc,sele)
          visc_sgi = bg_visc_sgi(1,1,:)

          ! grad(k**0.5) (dim, ngi)
          sqrt_k = ele_val(field2, ele)
          sqrt_k = sqrt(abs(sqrt_k))
          ! Lifted from ele_grad_at_quad function:
          do i=1, positions%dim
             grad_k(i,:) = matmul(sqrt_k, dshape(:,:,i))
          end do

          ! grad(k**0.5) at ele_nodes (dim,loc)
          q = shape_vector_rhs(shape, grad_k, detwei)
          !ewrite(3,*) "q: ", q

          q = matmul(q,invmass)
          !ewrite(3,*) "q: ", q
          ! Pick surface nodes (dim,sloc)
          q_s = q(:,face_local_nodes(field1,sele))
          !ewrite(3,*) "qs: ", qs
          q_sgi = matmul(q_s, fshape%n)
          !ewrite(3,*) "qsgi: ", qsgi

          !dot with surface normal
          do gi = 1, sgi
             q_sgin(gi) = dot_product(q_sgi(:,gi),normal_bdy(:,gi))
          end do
          !ewrite(3,*) "q_sgin: ", q_sgin

          ! integral of 2*nu*(grad(k**0.5))**2 (sloc)
          rhs = shape_rhs(fshape, detwei_bdy*q_sgin**2.0*visc_sgi*2.0)
          rhs = rhs/lumpedfmass
       end if

    ! high Re shear-stress wall functions for k and epsilon: see e.g. Wilcox (1994), Mathieu p.360
    else if(wall_fns=="high_Re") then
       visc_sgi = face_val_at_quad(EV,sele)
       grad_u = ele_grad_at_quad(u, ele, dshape)

       ! grad(U) at ele_nodes (dim,dim,loc)
       qq = shape_tensor_rhs(shape, grad_u, detwei)

       do i=1,u%dim
         do j=1,u%dim
           qq(i,j,:) = matmul(invmass,qq(i,j,:))
           ewrite(3,*) "qq: ", qq(i,j,:)

           ! Pick surface nodes (dim,dim,sloc)
           qq_s(i,j,:) = qq(i,j,face_local_nodes(field1,sele))
           ewrite(3,*) "qq_s: ", qq_s(i,j,:)

           ! Get values at surface quadrature (dim,dim,sgi)
           qq_sgi(i,j,:) = matmul(qq_s(i,j,:), fshape%n)
           ewrite(3,*) "qq_sgi: ", qq_sgi(i,j,:)
         end do
       end do

       do gi = 1, sgi
          !dot with surface normal (dim,sgi)
          qq_sgin(:,gi) = matmul(qq_sgi(:,:,gi),normal_bdy(:,gi))
          ewrite(3,*) "qq_sgin: ", qq_sgin(:,gi)

          ! Subtract normal component of velocity, leaving tangent components:
          qq_sgin(:,gi) = qq_sgin(:,gi)-normal_bdy(:,gi)*dot_product(qq_sgin(:,gi),normal_bdy(:,gi))

          ! Get streamwise component by taking sqrt(grad_n.grad_n). Multiply by eddy viscosity.
          ustar(gi) = norm2(qq_sgin(:,gi)) * visc_sgi(gi)
          ewrite(3,*) "ustar: ", ustar(gi)
       end do

       if (index==1) then
          rhs = shape_rhs(fshape, detwei_bdy*ustar/cmu**0.5)
          !ewrite(3,*) "cmu, rhs: ", cmu, rhs
          rhs = rhs/lumpedfmass
          !ewrite(3,*) "rhs: ", rhs
       else if (index==2) then
          ! calculate wall-normal element size
          G = matmul(transpose(invJ(:,:,1)), invJ(:,:,1))
          n(:,1) = normal_bdy(:,1)
          hb = 1. / sqrt( matmul(matmul(transpose(n), G), n) )
          h  = hb(1,1)
          ewrite(3,*) "h: ", h
          ! Von Karman's constant
          kappa = 0.43

          rhs = shape_rhs(fshape, detwei_bdy*ustar**1.5/kappa/h)
          rhs = rhs/lumpedfmass
       end if
    else
       FLAbort("Unknown wall function option for k_epsilon boundary conditions!")
    end if

end subroutine keps_wall_function

!------------------------------------------------------------------------------------------!
! Computes the strain rate for the LES model. Double-dot product results in a scalar:      !
! t:s = [ t11s11 + t12s21 + t13s31 + t21s12 + t22s22 + t23s32 + t31s13 + t32s23 + t33s33 ] !
!------------------------------------------------------------------------------------------!

function double_dot_product(du_t, nu)

    real, dimension(:,:,:), intent(in)         :: du_t   ! derivative of velocity shape function
    real, dimension(:,:), intent(in)           :: nu     ! nonlinear velocity
    real, dimension( size(du_t,2) )            :: double_dot_product
    real, dimension(size(du_t,3),size(du_t,3)) :: S, T
    integer                                    :: dim, ngi, gi, i, j

    ngi  = size(du_t,2)
    dim  = size(du_t,3)

    do gi=1, ngi
       !ewrite(2,*) "nu,du_t: ", nu, du_t(:,gi,:)
       S = matmul( nu, du_t(:,gi,:) )
       T = S
       !print *, "gi, S: ", gi, S
       S = S + transpose(S)
       double_dot_product(gi) = 0.

       do i = 1, dim
           do j = 1, dim
               double_dot_product(gi) = double_dot_product(gi) + T(i,j) * S(j,i)
               !ewrite(2,*) "i,j,T_ij,S_ji,ddp: ", i,j,T(i,j),S(j,i),double_dot_product(gi)
           end do
       end do
    end do

end function double_dot_product

end module k_epsilon
