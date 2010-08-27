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

module spontaneous_potentials
  !!< This module contains functions and subroutines to deal with spontaneous
  !!< potentials in porous media.  In particular, it contains a subroutine
  !!< to solve for the electrical potential field and electrical current
  !!< based on pressure potential, temperature and salinity gradients in
  !!< porous media.

  use fldebug
  use state_module
  use fields
  use spud
  use field_options
  use solvers
  use boundary_conditions
  use boundary_conditions_from_options
  use sparsity_patterns_meshes
  use porous_media, only: brine_salinity, normalised_saturation

  implicit none
  
  private

  public :: calculate_electrical_potential, &
            calculate_formation_conductivity

  contains

  subroutine calculate_electrical_potential(state,i)
    !!< determine electrical potential based on contribution of electrokinetic,
    !!< thermoelectric and electrochemical effects
    ! declare interface variables
    type(state_type), intent(inout) :: state
    integer, intent(in) :: i
    ! specify local variables
    logical :: have_sigma_fs
    integer :: j, ele, stat
    real :: dt
    character(len=255) :: tmpstring

    type(scalar_field), pointer :: pressure, temperature, salinity
    type(scalar_field), pointer :: ep_source, electrical_potential
    type(scalar_field) :: conductivity, L_ek, L_te, L_ec
    type(vector_field), pointer :: positions
    type(tensor_field), pointer :: ep_diffusivity

    type(mesh_type), pointer :: vmesh
    type(csr_sparsity), pointer :: sparsity => null()
    type(scalar_field) :: delta_ep
    type(csr_matrix) :: matrix
    type(scalar_field) :: rhs

    tmpstring = '/material_phase['//int2str(i-1)//']/electrical_properties/coupling_coefficients/'

    ! extract necessary fields for calculations
    vmesh => extract_velocity_mesh(state)
    positions => extract_vector_field(state, "Coordinate", stat=stat)

    ! extract electrical potential field
    electrical_potential => extract_scalar_field(state, "ElectricalPotential", stat=stat)
    if (stat/=0) then
      FLExit('Did not find an ElectricalPotential field.')
    end if
    if (i==0) call zero(electrical_potential)

    ! extract source for electrical potential field
    ep_source => extract_scalar_field(state, "ElectricalPotentialSource", stat=stat)
    if (stat/=0) then
      FLExit('Did not find an ElectricalPotentialSource field.')
    end if
    call zero(ep_source)

    ! determine electrical conductivity of formation
    conductivity = extract_scalar_field(state, "ElectricalConductivity", stat=stat)
    if (stat==0) then
      have_sigma_fs = .true.
    else
      ewrite(3,*) 'Computing electrical conductivity of formation.'
      call allocate(conductivity, vmesh, 'ElectricalConductivity')
      call calculate_formation_conductivity(state, i, conductivity, stat)
      ewrite_minmax(conductivity)
      have_sigma_fs = .false.
    end if

    ! ELECTROKINETIC
    if (have_option(trim(tmpstring)//'scalar_field::Electrokinetic')) then
      ewrite(3,*) 'Computing electrokinetic coupling term.'
      ! initialise coupling term
      call allocate(L_ek, vmesh, "L_ek")

      ! determine electrokinetic coupling term
      call calculate_coupling_term(state, i, L_ek, 'Electrokinetic', conductivity, stat)
      ewrite_minmax(L_ek)

      ! extract pressure field
      pressure => extract_scalar_field(state, "Pressure", stat=stat)
      if (stat/=0) then
        FLExit('Did not find a Pressure field.')
      end if

      ! assemble electrical potential source using pressure
      do ele=1,element_count(electrical_potential)
        call assemble_ep_source_element(ele, ep_source, positions, electrical_potential, pressure, L_ek)
      end do
      ewrite_minmax(ep_source)

      ! clean up
      call deallocate(L_ek)
    end if

    ! THERMOELECTRIC
    if (have_option(trim(tmpstring)//'scalar_field::Thermoelectric')) then
      ewrite(3,*) 'Computing thermoelectric coupling term.'
      ! initialise coupling term
      call allocate(L_te, vmesh, "L_te")

      ! determine thermoelectric coupling term
      call calculate_coupling_term(state, i, L_te, 'Thermoelectric', conductivity, stat)
      ewrite_minmax(L_te)

      ! extract temperature field
      temperature => extract_scalar_field(state, "Temperature", stat=stat)
      if (stat/=0) then
        FLExit('Did not find a Temperature field.')
      end if

      ! assemble electrical potential source
      do ele=1,element_count(electrical_potential)
        call assemble_ep_source_element(ele, ep_source, positions, electrical_potential, temperature, L_te)
      end do
      ewrite_minmax(ep_source)

      ! clean up
      call deallocate(L_te)
    end if

    ! ELECTROCHEMICAL
    if (have_option(trim(tmpstring)//'scalar_field::Electrochemical')) then
      ewrite(3,*) 'Computing electrochemical coupling term.'
      ! initialise coupling term
      call allocate(L_ec, vmesh, "L_ec")

      ! determine electrochemical coupling term
      call calculate_coupling_term(state, i, L_ec, 'Electrochemical', conductivity, stat)
      ewrite_minmax(L_ec)

      ! extract salinity field
      salinity => extract_scalar_field(state, "Salinity", stat=stat)
      if (stat/=0) then
        FLExit('Did not find a Salinity field.')
      end if

      ! assemble electrical potential source
      do ele=1,element_count(electrical_potential)
        call assemble_ep_source_element(ele, ep_source, positions, electrical_potential, salinity, L_ec)
      end do
      ewrite_minmax(ep_source)

      ! clean up
      call deallocate(L_ec)
    end if

    ! set electrical potential diffusivity term to the formation conductivity
    ep_diffusivity => extract_tensor_field(state, "ElectricalPotentialDiffusivity", stat)
    if (stat/=0) then
      FLExit('Did not find an ElectricalPotentialDiffusivity field')
    end if
    call zero(ep_diffusivity)
    do j=1,ep_diffusivity%dim
      call set(ep_diffusivity, j, j, conductivity, symmetric=.true.)
    end do

    ! initialise sparse matrix and vectors
    sparsity => get_csr_sparsity_firstorder(state, electrical_potential%mesh, electrical_potential%mesh)
    call allocate(matrix, sparsity, name = "ElectricalMatrixLHS")
    call zero(matrix)
    call allocate(rhs, electrical_potential%mesh, name = "ElectricalRHS")
    call zero(rhs)
    call allocate(delta_ep, electrical_potential%mesh, name = "delta_ep")
    call zero(delta_ep)

    ! assemble electrical potential diffusivity
    call get_option("/timestepping/timestep",dt)
    do ele=1,element_count(electrical_potential)
      call assemble_ep_diffusivity_element(ele, electrical_potential, positions, ep_diffusivity, matrix, rhs, dt)
    end do
    call addto(rhs, ep_source)

    ! solve for electrical potential
    call petsc_solve(delta_ep, matrix, rhs, option_path = electrical_potential%option_path)
    call addto(electrical_potential, delta_ep, dt)
    ewrite_minmax(electrical_potential)

    ! deallocate memory
    call deallocate(matrix)
    call deallocate(rhs)
    call deallocate(delta_ep)
    if (.not.have_sigma_fs) call deallocate(conductivity)

    stat = 0
  end subroutine calculate_electrical_potential

  subroutine assemble_ep_source_element(ele, ep_source, positions, e_potential, f_potential, Lx)
    ! declare interface variables
    integer, intent(in) :: ele
    type(scalar_field), intent(inout) :: ep_source
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: e_potential, f_potential, Lx
    ! declare local variables
    integer :: i, j
    ! global node numbers of the current element
    integer, dimension(:), pointer :: f_ele
    ! variable transform times quadrature weights
    real, dimension(ele_ngi(f_potential,ele)) :: detwei
    ! transformed electrical potential shape function
    real, dimension(ele_loc(e_potential, ele), ele_ngi(e_potential, ele), mesh_dim(e_potential)) :: de_t
    ! transformed fluid potential shape function
    real, dimension(ele_loc(f_potential, ele), ele_ngi(f_potential, ele), mesh_dim(f_potential)) :: df_t
    ! coupling term at quadrature
    real, dimension(ele_ngi(f_potential, ele)) :: lx_gi
    ! values
    real, dimension(ele_loc(f_potential, ele)) :: f_val
    ! shape of the current element
    type(element_type), pointer :: e_shape, f_shape

    e_shape => ele_shape(e_potential, ele)

    f_ele => ele_nodes(f_potential, ele)
    f_shape => ele_shape(f_potential, ele)
    f_val =  ele_val(f_potential, ele)

    call transform_to_physical(positions, ele, e_shape, dshape = de_t, detwei=detwei)
    call transform_to_physical(positions, ele, f_shape, dshape = df_t, detwei=detwei)

    ! /
    ! | grad N_v dot grad N_p (L_x P) dV
    ! /
    lx_gi = ele_val_at_quad(Lx, ele)
    call addto(ep_source, f_ele, matmul(dshape_dot_dshape(de_t, df_t, lx_gi*detwei), f_val))
  end subroutine assemble_ep_source_element

  subroutine assemble_ep_diffusivity_element(ele, t, positions, diffusivity, matrix, rhs, dt)
    ! declare interface variables
    integer, intent(in) :: ele
    type(scalar_field), intent(in) :: t
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(in) :: diffusivity
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    real, intent(inout) :: dt
    ! declare local variables
    real, dimension(ele_ngi(t, ele)) :: detwei
    real, dimension(ele_loc(t, ele), ele_ngi(t, ele), mesh_dim(t)) :: dt_t
    real, dimension(ele_loc(t, ele)) :: rhs_addto
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)) :: matrix_addto
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)) :: diffusivity_mat
    real, dimension(diffusivity%dim, diffusivity%dim, ele_ngi(diffusivity, ele)) :: diffusivity_gi
    type(element_type), pointer :: t_shape
    logical :: isotropic_diffusivity
    ! values
    real, dimension(ele_loc(t, ele)) :: t_val

    diffusivity_gi = ele_val_at_quad(diffusivity, ele)

    rhs_addto = 0.0
    matrix_addto = 0.0

    t_shape => ele_shape(t, ele)
    t_val = ele_val(t, ele)

    call transform_to_physical(positions, ele, t_shape, dshape = dt_t, detwei = detwei)

    isotropic_diffusivity = option_count(complete_field_path(diffusivity%option_path)) &
    & == option_count(trim(complete_field_path(diffusivity%option_path)) // "/value/isotropic")
    if(isotropic_diffusivity) then
      assert(size(diffusivity_gi, 1) > 0)
      diffusivity_mat = dshape_dot_dshape(dt_t, dt_t, detwei * diffusivity_gi(1, 1, :))
    else
      diffusivity_mat = dshape_tensor_dshape(dt_t, diffusivity_gi, dt_t, detwei)
    end if

    matrix_addto = matrix_addto + diffusivity_mat
    rhs_addto = rhs_addto - matmul(diffusivity_mat, ele_val(t, ele))

    call addto(matrix, ele_nodes(t, ele), ele_nodes(t, ele), matrix_addto)
    call addto(rhs, ele_nodes(t, ele), rhs_addto)

    call apply_dirichlet_conditions(matrix, rhs, t, dt)

  end subroutine assemble_ep_diffusivity_element

  subroutine calculate_formation_conductivity(state, i, sigma_fs, stat)
    !!< calculate the electrical conductivity of the porous medium as a diagnostic
    !!< field from Archie's law as equation (33) of Saunders et al. (2008).
    ! declare interface variables
    type(state_type), intent(in) :: state
    integer, intent(in) :: i
    type(scalar_field), intent(inout) :: sigma_fs
    integer, intent(out) :: stat
    ! declare local variables
    logical :: have_saturation
    integer :: j, n
    real :: sigma_o, tmp, archie
    type(scalar_field), pointer :: porosity, salinity
    type(scalar_field) :: saturation, sigma_w

    ! extract porosity field
    porosity => extract_scalar_field(state, "Porosity", stat)
    if (stat/=0) then
      FLExit('Did not find a Porosity field.')
    end if

    ! extract phase saturation if present otherwise set to 1.0
    saturation = extract_scalar_field(state, "PhaseVolumeFraction", stat)
    if (stat==0) then
      have_saturation = .true.
    else
      call allocate(saturation, sigma_fs%mesh, 'PhaseVolumeFraction')
      call zero(saturation)
      call addto(saturation,1.0)
      have_saturation = .false.
    end if
    ewrite_minmax(saturation)

    ! get electrical conductivity of oil
    ! NB: hard coding value for now
    sigma_o = 1.0e-5

    ! get electrical conductivity of brine
    call allocate(sigma_w, sigma_fs%mesh, 'FluidConductivity')
    call calculate_phase_conductivity(state, i, sigma_w, stat)
   
    ! determine scalar field continuity
    if (sigma_fs%mesh%continuity==-1) then
      ! mesh is discontinuous
      n = ele_count(sigma_fs)
    else
      ! mesh is continuous
      n = node_count(sigma_fs)
    end if

    ! compute electrical conductivity of the formation
    if (have_option("/material_phase["//int2str(i-1)//"]/electrical_properties/archie_exponent")) then
      ! get electrical conductivity of brine
      call get_option("/material_phase["//int2str(i-1)//"]/electrical_properties/archie_exponent",archie)
    else
      archie=1.65
    endif
    call zero(sigma_fs)
    do j=1,n
      if (porosity%val(j) > 0.29) then
         saturation%val(j) = 1.0
      endif
      tmp = saturation%val(j)**archie
      call set(sigma_fs, j, (porosity%val(j)**1.8)*(sigma_w%val(j)*tmp+sigma_o*(1.0-tmp)))
    end do

    ! clean up
    if (.not.have_saturation) call deallocate(saturation)
    call deallocate(sigma_w)

    stat = 0
  end subroutine calculate_formation_conductivity

  subroutine calculate_phase_conductivity(state, i, sigma, stat)
    !!< calculate the electrical conductivity of fluid phase depending on the options 
    ! declare interface variables
    type(state_type), intent(in) :: state
    integer, intent(in) :: i
    type(scalar_field), intent(inout) :: sigma
    integer, intent(out) :: stat
    ! declare local variables
    integer :: j, it, it_max, n
    real :: sigma_lo, sigma_hi, tmp
    type(scalar_field), pointer :: salinity

    ! limit number of iterations in bisection method
    it_max = 1000

    ! determine scalar field continuity
    if (sigma%mesh%continuity==-1) then
      ! mesh is discontinuous
      n = ele_count(sigma)
    else
      ! mesh is continuous
      n = node_count(sigma)
    end if

    if (have_option("/material_phase["//int2str(i-1)//"]/electrical_properties/conductivity")) then
      ! get electrical conductivity of brine
      call get_option("/material_phase["//int2str(i-1)//"]/electrical_properties/conductivity",tmp)
      call zero(sigma)
      call addto(sigma,tmp)
    elseif (have_option("/material_phase["//int2str(i-1)//"]/electrical_properties/conductivity_from_salinity")) then
      ! extract salinity field
      salinity => extract_scalar_field(state, "Salinity", stat=stat)
      if (stat/=0) then
        FLExit('Did not find a Salinity field.')
      end if

      do j=1,n
        ! use bisection method to solve for electrical conductivity of brine
        sigma_lo = 1.4e-5 ! salinity ~ 1.054e-6 mol/L
        sigma_hi = 65.0   ! salinity ~ 10.130 mol/L
        sigma%val(j) = sigma_lo + 0.5*(sigma_hi-sigma_lo)
        it = 0
        do while ((sigma%val(j)/=sigma_lo).and.(sigma%val(j)/=sigma_hi).and.(it<it_max))
          if ((brine_salinity(sigma%val(j))-salinity%val(j))<=0.0) then
            sigma_lo = sigma%val(j)
          else
            sigma_hi = sigma%val(j)
          end if
          sigma%val(j) = sigma_lo + 0.5*(sigma_hi-sigma_lo)
          it = it+1
        end do
      end do
    elseif (have_option("/material_phase["//int2str(i-1)//"]/electrical_properties/conductivity_from_salinity_and_temperature")) then
      ! NOT YET IMPLEMENTED
      FLExit('Cannot currently compute electrical conductivity of fluid from salinity and temperature.  Sorry.')
    else
      ! should never get here
      FLAbort('Logic error! Should never get here!')
    end if

    stat = 0
  end subroutine calculate_phase_conductivity

  subroutine calculate_coupling_term(state, i, L_x, coupling, sigma_fs, stat)
    !!< calculate electrokinetic, thermoelectric or electrochemical coupling
    !!< term as specified by coupling argument
    ! declare interface variables
    type(state_type), intent(in) :: state
    integer, intent(in) :: i
    type(scalar_field), intent(inout) :: L_x
    character(len=*), intent(in) :: coupling
    integer, intent(out) :: stat
    ! declare local variables
    logical :: have_sigma_fs
    integer :: n, j
    type(scalar_field), pointer :: Cv
    type(scalar_field) :: sigma_fs

    stat = -1

    ! determine scalar field continuity
    if (L_x%mesh%continuity==-1) then
      ! mesh is discontinuous
      n = ele_count(L_x)
    else
      ! mesh is continuous
      n = node_count(L_x)
    end if

!! Unnecessary now we're passing conductivity down into this subroutine:

!    ! determine electrical conductivity of formation
!    sigma_fs = extract_scalar_field(state, "ElectricalConductivity", stat=stat)
!    ! ewrite(3,*) 'sigma stat',stat
!    if (stat==0) then
!      have_sigma_fs = .true.
!    else
!      ewrite(3,*) 'Computing formation conductivity for coupling coefficient'
!      call allocate(sigma_fs, L_x%mesh, 'ElectricalConductivity')
!      call calculate_formation_conductivity(state, i, sigma_fs, stat)
!      have_sigma_fs = .false.
!    end if
!    ewrite_minmax(sigma_fs)
    
    ! get coupling coefficient
    Cv => extract_scalar_field(state, trim(coupling)//'['//int2str(i-1)//']', stat=stat)
    if (stat/=0) then
      FLExit('Did not find a '//trim(coupling)//' coupling coefficient scalar field.')
    end if
    ewrite_minmax(Cv)

    ! compute coupling term
    call zero(L_x)
    do j=1,n
      call set(L_x, j, Cv%val(j)*sigma_fs%val(j))
    end do

!    ! clean up
!    if (.not.have_sigma_fs) call deallocate(sigma_fs)

    stat = 0
  end subroutine calculate_coupling_term

end module spontaneous_potentials

