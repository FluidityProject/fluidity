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

module subgrid_kinetic_energy
  use quadrature
  use elements
  use fields
  use field_options
  use state_module
  use spud
  use state_fields_module
  use fetools
  use vector_tools
  use sparsity_patterns_meshes
  use k_epsilon
  use smoothing_module
  use FLDebug
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN

implicit none

  private
  public :: compute_subgrid_tke

contains

subroutine compute_subgrid_tke(state)

    type(state_type), intent(inout)    :: state
    ! Filter width ratio
    real                   :: alpha
    ! Tensor or scalar filter. Only scalar for now.
    logical                :: anisotropy
    type(scalar_field), pointer        :: k, src_k, abs_k, lumped_mass
    type(scalar_field)                 :: src_rhs, abs_rhs, diff_rhs
    type(vector_field), pointer        :: positions, nu
    type(tensor_field), pointer        :: eddyvisc, visc, diff_k
    integer                            :: i, ele, dims, stat
    integer, dimension(:),pointer      :: knodes
    character(len=OPTION_PATH_LEN) :: les_option_path

    ewrite(1,*) "In subgrid_tke"

    positions       => extract_vector_field(state, "Coordinate",stat)
    nu              => extract_vector_field(state, "Velocity",stat)
    k               => extract_scalar_field(state, "SubgridKineticEnergy",stat)
    src_k           => extract_scalar_field(state, "SubgridKineticEnergySource",stat)
    if(stat/=0) then
      FLExit('You have to have a SubgridKineticEnergy Source field set to diagnostic/internal.')
    end if
    abs_k           => extract_scalar_field(state, "SubgridKineticEnergyAbsorption",stat)
    if(stat/=0) then
      FLExit('You have to have a SubgridKineticEnergy Absorption field set to diagnostic/internal.')
    end if
    diff_k          => extract_tensor_field(state, "SubgridKineticEnergyDiffusivity",stat)
    if(stat/=0) then
      FLExit('You have to have a SubgridKineticEnergyDiffusivity field set to diagnostic/internal.')
    end if
    visc            => extract_tensor_field(state, "Viscosity",stat)
    if(stat/=0) then
      FLExit('You have to have a Viscosity field supplied by one of the turbulence models for this field to be calculated.')
    end if
    eddyvisc        => extract_tensor_field(state, "EddyViscosity", stat)
    if(stat/=0) then
      FLExit('You have to have an EddyViscosity field supplied by one of the turbulence models for this field to be calculated.')
    end if

    ewrite(2,*) 'Allocating local SGS TKE fields'

    call allocate(src_rhs, k%mesh, name="KSRCRHS")
    call allocate(abs_rhs, k%mesh, name="KABSRHS")
    call allocate(diff_rhs, k%mesh, name="KDIFFRHS")
    call zero(src_rhs); call zero(abs_rhs); call zero(diff_rhs)

    les_option_path=(trim(nu%option_path)//"/prognostic/spatial_discretisation"//&
         &"/continuous_galerkin/les_model")
    ewrite(1,*) 'LES option path ', trim(les_option_path)
    anisotropy = .false.
    anisotropy = have_option(trim(les_option_path)//"/dynamic_les/anisotropic_viscosity") .or.&
                 have_option(trim(les_option_path)//"/second_order/anisotropic_viscosity")
    call get_option(trim(les_option_path)//"/dynamic_les/alpha", alpha, default=2.0)
    ewrite(2,*) 'subgrid TKE anisotropy, alpha: ', anisotropy, alpha

    ! Assembly loop
    do ele = 1, ele_count(k)
       ! Realisability condition: k=>0 because we use k**0.5 in src and abs terms
       knodes => ele_nodes(k,ele)
       call set(k, knodes, max(node_val(k,knodes),0.))
       call assemble_subgrid_tke_ele(src_rhs, abs_rhs, diff_rhs, k, nu, eddyvisc, visc, positions, ele, alpha, anisotropy)
    end do

    lumped_mass => get_lumped_mass(state, k%mesh)

    ewrite(2,*) "Calculating k source and absorption"
    do i = 1, node_count(k)
      call set(src_k, i, node_val(src_rhs,i)/node_val(lumped_mass,i))
      call set(abs_k, i, node_val(abs_rhs,i)/node_val(lumped_mass,i))
      do dims = 1, diff_k%dim(1)
        call set(diff_k, dims, dims, i, node_val(diff_rhs,i)/node_val(lumped_mass,i))
      end do
    end do

    call deallocate(src_rhs); call deallocate(abs_rhs); call deallocate(diff_rhs)

    ewrite_minmax(eddyvisc)
    ewrite_minmax(k)
    ewrite_minmax(src_k)
    ewrite_minmax(abs_k)
    ewrite_minmax(diff_k)

end subroutine compute_subgrid_tke

subroutine assemble_subgrid_tke_ele(src_rhs, abs_rhs, diff_rhs, k, u, eddyvisc, visc, x, ele, alpha, anisotropy)

  type(scalar_field), intent(inout) :: src_rhs, abs_rhs, diff_rhs
  type(scalar_field), intent(in)    :: k
  type(vector_field), intent(in)    :: x, u
  type(tensor_field), intent(in)    :: eddyvisc, visc
  integer, intent(in)               :: ele
  real, intent(in)                  :: alpha
  logical, intent(in)               :: anisotropy

  real, dimension(ele_loc(k, ele), ele_ngi(k, ele), x%dim) :: dshape_x
  real, dimension(ele_ngi(k, ele))                         :: detwei, strain_ngi
  real, dimension(ele_loc(k, ele))                         :: rhs_addto
  real, dimension(x%dim, x%dim, ele_ngi(k, ele))           :: ft, evt, nut
  real, dimension(ele_ngi(k, ele))                         :: f, ev, nu
  integer, dimension(ele_loc(k, ele))                      :: nodes_k
  type(element_type), pointer                              :: shape_k, shape_x
  integer                                                  :: i, gi
  real                                                     :: C_s, C_e, sigma_k

  shape_k => ele_shape(k, ele)
  shape_x => ele_shape(x, ele)
  nodes_k = ele_nodes(k, ele)

  ! Model coefficients
  C_e = 0.7 ! Absorption coeff from Schmidt & Schumann 1989
  C_s = 0.1 ! Smagorinsky - hard code to 0.1 for now
  sigma_k = 1.0

  call transform_to_physical(x, ele, shape_x, dshape=dshape_x, detwei=detwei)

  ! Calculate TKE production at ngi using strain rate (double_dot_product) function
  strain_ngi = double_dot_product(dshape_x, ele_val(u, ele) )

  !if(.not. anisotropy) then
  ! Scalar terms
    ! mesh metric tensor (units length^2)
    ft = length_scale_tensor(dshape_x, ele_shape(u, ele))
    evt = ele_val_at_quad(eddyvisc, ele)
    nut = ele_val_at_quad(visc, ele)
    do gi=1, ele_ngi(u, ele)
      ! 2-norm of filter size. Actually want units of length so sqrt
      f(gi) = sqrt(alpha**2*norm2(ft(:,:,gi)))
      ! 2-norm of eddy viscosity
      ev(gi) = norm2(evt(:,:,gi))
      ! 2-norm of viscosity
      nu(gi) = norm2(nut(:,:,gi))
    end do
  ! not doing tensor filter yet
  !else
    ! k has to be a tensor
  !end if

  !ewrite(2,*) 'SGS src terms: ', strain_ngi, ev
  !ewrite(2,*) 'SGS abs terms: ', C_e, f, ele_val_at_quad(k,ele)**0.5

  ! Source term: depends on filter size via eddy viscosity
  rhs_addto = shape_rhs(shape_k, detwei*strain_ngi*ev)
  call addto(src_rhs, nodes_k, rhs_addto)

  ! Absorption term: depends on inverse filter size
  rhs_addto = shape_rhs(shape_k, detwei*C_e/f*ele_val_at_quad(k,ele)**0.5)
  call addto(abs_rhs, nodes_k, rhs_addto)

  ! Dissipation term: depends on filter size
  rhs_addto = shape_rhs(shape_k, detwei*(nu+C_s/sigma_k*f*ele_val_at_quad(k,ele)**0.5))
  call addto(diff_rhs, nodes_k, rhs_addto)

end subroutine

end module subgrid_kinetic_energy
