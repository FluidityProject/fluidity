!    Copyright (C) 2008 Imperial College London and others.
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

module les_viscosity_module
  !!< This module contains several subroutines and functions used to implement LES models
  use state_module
  use fields
  use field_options
  use spud
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use smoothing_module
  use vector_tools
  use fetools
  implicit none

  private

  public les_viscosity_strength, wale_viscosity_strength
  public les_zero_diagnostic_fields, les_set_diagnostic_fields
  public dynamic_les_init_fields, dynamic_les_filtered_fields, les_strain_rate

contains

  subroutine les_zero_diagnostic_fields(state)

    type(state_type), intent(inout) :: state
    type(tensor_field), pointer :: tfield
    type(scalar_field), pointer :: sfield
    integer :: stat

    ! EddyViscosity is an option for all LES models.
    tfield => extract_tensor_field(state, "EddyViscosity", stat)
    if(stat==0) then
       call zero(tfield)
    end if
    ! Following fields are optional for dynamic LES model.
    tfield => extract_tensor_field(state, "StrainRate", stat)
    if(stat==0) then
       call zero(tfield)
    end if
    tfield => extract_tensor_field(state, "FilteredStrainRate", stat)
    if(stat==0) then
       call zero(tfield)
    end if
    tfield => extract_tensor_field(state, "FilterWidth", stat)
    if(stat==0) then
       call zero(tfield)
    end if
    sfield => extract_scalar_field(state, "SmagorinskyCoefficient", stat)
    if(stat==0) then
       call zero(sfield)
    end if

  end subroutine les_zero_diagnostic_fields

  subroutine dynamic_les_init_fields(state, nu, filtered_velocity, leonard_tensor, strain_product)

    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout)  :: nu
    type(vector_field), pointer, intent(inout) :: filtered_velocity
    type(tensor_field), pointer, intent(inout) :: leonard_tensor, strain_product
    integer :: stat

    ! Filtered velocity
    filtered_velocity => extract_vector_field(state, "FilteredVelocity", stat)
    if(stat/=0) then
      ewrite(2,*) "Allocating FilteredVelocity on VelocityMesh"
      allocate(filtered_velocity)
      call allocate(filtered_velocity, nu%dim, nu%mesh, "FilteredVelocity")
      ewrite(2,*) "FilteredVelocity refcount: ", filtered_velocity%refcount%count
    end if
    call zero(filtered_velocity)
    ewrite_minmax(filtered_velocity)

    ! Leonard tensor
    leonard_tensor => extract_tensor_field(state, "LeonardTensor", stat)
    if(stat/=0) then
      ewrite(2,*) "Allocating LeonardTensor on VelocityMesh"
      allocate(leonard_tensor)
      call allocate(leonard_tensor, nu%mesh, "LeonardTensor")
    end if
    call zero(leonard_tensor)

    ! Filtered strain product
    strain_product => extract_tensor_field(state, "StrainProduct", stat)
    if(stat/=0) then
      ewrite(2,*) "Allocating StrainProduct on VelocityMesh"
      allocate(strain_product)
      call allocate(strain_product, nu%mesh, "StrainProduct")
    end if
    call zero(strain_product)

  end subroutine dynamic_les_init_fields

  subroutine les_set_diagnostic_fields(state, nu, density, ele, detwei, &
                 eddy_visc_gi, visc_gi, strain_gi, filtered_strain_gi, filter_width_gi, les_coef_gi)

    type(state_type), intent(inout)                             :: state
    type(vector_field), intent(in)                              :: nu
    type(scalar_field), intent(in)                              :: density
    integer, intent(in)                                         :: ele
    real, dimension(ele_ngi(nu,ele)), intent(in)                :: detwei
    real, dimension(ele_ngi(nu,ele)), intent(in), optional      :: les_coef_gi
    real, dimension(nu%dim,nu%dim,ele_ngi(nu,ele)),intent(inout), optional &
    & :: eddy_visc_gi, visc_gi, strain_gi, filtered_strain_gi, filter_width_gi
    real, dimension(ele_ngi(nu,ele))                            :: density_gi
    type(tensor_field), pointer                                 :: tensorfield
    type(scalar_field), pointer                                 :: scalarfield
    real, dimension(nu%dim,nu%dim,ele_loc(nu,ele))              :: tensor_loc
    real, dimension(ele_loc(nu,ele))                            :: scalar_loc, lumped_mass
    real, dimension(ele_loc(nu,ele),ele_loc(nu,ele))            :: mass_matrix
    type(element_type)                                          :: nu_shape
    integer                                                     :: stat, i

    nu_shape = ele_shape(nu, ele)
    density_gi = ele_val_at_quad(density, ele)
    mass_matrix = shape_shape(nu_shape, nu_shape, detwei*density_gi)
    lumped_mass = sum(mass_matrix, 2)

    ! Eddy viscosity field m_ij is available with all LES models
    if (present(eddy_visc_gi)) then
      !ewrite(2,*) "eddy visc: ", eddy_visc_gi(:,:,1)
      !ewrite(2,*) "visc: ", visc_gi
      !eddy_visc_gi = eddy_visc_gi + visc_gi
      tensorfield => extract_tensor_field(state, "EddyViscosity",stat)
      if(stat==0) then
        !if(any(eddy_visc_gi(:,:,:)<0.)) then
          !ewrite(2,*) "warning: visc_gi+eddy_visc_gi < 0", eddy_visc_gi
        !end if
        tensor_loc=shape_tensor_rhs(nu_shape, eddy_visc_gi, detwei)
        ! Divide by lumped mass
        do i=1, ele_loc(tensorfield,ele)
          tensor_loc(:,:,i)=tensor_loc(:,:,i)/lumped_mass(i)
        end do
        call addto(tensorfield, ele_nodes(tensorfield, ele), tensor_loc)
      end if
    end if

    ! Strain rate field S1
    if (present(strain_gi)) then
      tensorfield => extract_tensor_field(state, "StrainRate", stat)
      if(stat==0) then
        tensor_loc=shape_tensor_rhs(nu_shape, strain_gi, detwei)
        ! Divide by lumped mass
        do i=1, ele_loc(nu,ele)
          tensor_loc(:,:,i)=tensor_loc(:,:,i)/lumped_mass(i)
        end do
        call addto(tensorfield, ele_nodes(nu, ele), tensor_loc)
      else
        !FlExit("Error: StrainRate field has disappeared.")
      end if
    end if

    ! Filtered strain rate field S2
    if (present(filtered_strain_gi)) then
      tensorfield => extract_tensor_field(state, "FilteredStrainRate", stat)
      if(stat==0) then
        tensor_loc=shape_tensor_rhs(nu_shape, filtered_strain_gi, detwei)
        ! Divide by lumped mass
        do i=1, ele_loc(nu,ele)
          tensor_loc(:,:,i)=tensor_loc(:,:,i)/lumped_mass(i)
        end do
        call addto(tensorfield, ele_nodes(nu, ele), tensor_loc)
      else
        !FlExit("Error: FilteredStrainRate field has disappeared.")
      end if
    end if

    ! Filter width
    if (present(filter_width_gi)) then
      tensorfield => extract_tensor_field(state, "FilterWidth", stat)
      if(stat==0) then
        tensor_loc=shape_tensor_rhs(nu_shape, filter_width_gi, detwei)
        ! Divide by lumped mass
        do i=1, ele_loc(nu,ele)
          tensor_loc(:,:,i)=tensor_loc(:,:,i)/lumped_mass(i)
        end do
        call addto(tensorfield, ele_nodes(nu, ele), tensor_loc)
      else
        !FlExit("Error: FilterWidth field has disappeared.")
      end if
    end if

    ! Smagorinsky coefficient field
    if (present(les_coef_gi)) then
      scalarfield => extract_scalar_field(state, "SmagorinskyCoefficient",stat)
      if(stat==0) then
        scalar_loc=shape_rhs(nu_shape, les_coef_gi*detwei)
        ! Divide by lumped mass
        do i=1, ele_loc(nu,ele)
          scalar_loc(i)=scalar_loc(i)/lumped_mass(i)
        end do
        call addto(scalarfield, ele_nodes(nu, ele), scalar_loc)
      else
        !FlExit("Error: SmagorinskyCoefficient field has disappeared.")
      end if
    end if

  end subroutine les_set_diagnostic_fields

  subroutine dynamic_les_filtered_fields(state, nu, positions, tnu, leonard, strain_prod, alpha, les_option_path)

    type(state_type), intent(inout)           :: state
    ! Unfiltered velocity
    type(vector_field), pointer, intent(inout):: nu
    type(vector_field), intent(in)            :: positions
    ! Filtered velocity
    type(vector_field), pointer, intent(inout):: tnu
    ! Leonard tensor and strain product fields
    type(tensor_field), pointer, intent(inout):: leonard, strain_prod
    ! Scale factor: test filter/mesh size
    real, intent(in)                          :: alpha
    character(len=OPTION_PATH_LEN), intent(in):: les_option_path
    ! Local quantities
    type(tensor_field), pointer               :: tensorfield
    character(len=OPTION_PATH_LEN)            :: lpath
    integer                                   :: i,gi
    real, dimension(:), allocatable           :: u_loc
    real, dimension(:,:), allocatable         :: t_loc

    real, dimension(ele_loc(nu,1), ele_ngi(nu,1), nu%dim) :: du_t
    real, dimension(ele_ngi(nu,1))                 :: detwei
    real, dimension(nu%dim, nu%dim, ele_ngi(nu,1)) :: strain_gi, t_strain_gi, strain_prod_gi, filter_width_gi
    real, dimension(ele_ngi(nu,1))                 :: strain_mod, t_strain_mod
    type(element_type)                             :: shape_nu
    integer, dimension(:), pointer                 :: nodes_nu

    ! Path to solver options
    lpath = (trim(les_option_path)//"/dynamic_les")
    ewrite(2,*) "filter factor alpha: ", alpha
    ewrite(2,*) "path to solver options: ", trim(lpath)
    ! Filter the nonlinear velocity
    call anisotropic_smooth_vector(nu, positions, tnu, alpha, lpath)

    ewrite_minmax(nu)
    ewrite_minmax(tnu)

    ! Local variables
    allocate(tensorfield)
    call allocate(tensorfield, nu%mesh, "TensorField")
    call zero(tensorfield)
    allocate(u_loc(nu%dim)); allocate(t_loc(nu%dim, nu%dim))
    u_loc=0.0; t_loc=0.0

    ! Get mesh-filtered velocity cross-product
    do i=1, node_count(nu)
      u_loc = node_val(nu,i)
      t_loc = outer_product(u_loc, u_loc)
      call set( tensorfield, i, t_loc )
    end do

    ! Calculate test-filtered (velocity cross-product): (ui^r*uj^r)^t
    call anisotropic_smooth_tensor(tensorfield, positions, leonard, alpha, lpath)
    call zero(tensorfield)

    ! Calculate (test-filtered velocity) cross-product: (ui^rt*uj^rt)
    do i=1, node_count(nu)
      u_loc = node_val(tnu,i)
      t_loc = outer_product(u_loc, u_loc)
      call set( tensorfield, i, t_loc )
    end do

    ! Leonard tensor field: (ui^r*uj^r)^t - (ui^rt*uj^rt)
    call addto( leonard, tensorfield, -1.0 )
    call zero(tensorfield)
    ewrite_minmax(leonard)

    ! Strain rate S1
    do i=1, element_count(nu)
      shape_nu = ele_shape(nu, i)

      ! Assuming no stabilisation scheme is used with LES so I can use velocity shape. (Tejada-Martinez thinks otherwise)
      call transform_to_physical(positions, i, shape_nu, dshape=du_t, detwei=detwei)

      ! Get strain S1 for unfiltered velocity (dim,dim,ngi)
      strain_gi = les_strain_rate(du_t, ele_val(nu, i))
      ! Get strain S2 for test-filtered velocity (dim,dim,ngi)
      t_strain_gi = les_strain_rate(du_t, ele_val(tnu, i))
      ! Get strain modulus |S1| for unfiltered velocity (ngi)
      strain_mod = les_viscosity_strength(du_t, ele_val(nu, i))
      ! Get strain modulus |S2| for test-filtered velocity (ngi)
      t_strain_mod = les_viscosity_strength(du_t, ele_val(tnu, i))

      ! Strain product |S1|S1
      do gi=1, ele_ngi(nu, i)
        strain_prod_gi(:,:,gi) = strain_gi(:,:,gi)*strain_mod(gi)
      end do

      ! Set some diagnostic fields
      !call les_set_diagnostic_fields(state, nu, i, detwei, &
      !     strain_gi=strain_gi, filtered_strain_gi=t_strain_gi, filter_width_gi=filter_width_gi)

      ! Strain product field |S1|S1
      call addto(tensorfield, ele_nodes(nu,i), shape_tensor_rhs(ele_shape(nu,i), strain_prod_gi, detwei))
    end do
    ewrite_minmax(tensorfield)
    ! Filter the strain product: (|S1|S1)^t (not a diagnostic field yet)
    call anisotropic_smooth_tensor(tensorfield, positions, strain_prod, alpha, lpath)
    ewrite_minmax(strain_prod)

    ewrite(2,*) "Finished setting LES filtered fields"
    ! Deallocates
    deallocate(u_loc, t_loc)
    call deallocate(tensorfield)
    deallocate(tensorfield)

  end subroutine dynamic_les_filtered_fields

  function les_strain_rate(du_t, nu)
    !! Computes the strain rate
    !! derivative of velocity shape function (nloc x ngi x dim)
    real, dimension(:,:,:), intent(in):: du_t
    !! nonlinear velocity (dim x nloc)
    real, dimension(:,:), intent(in):: nu
    real, dimension( size(du_t,3),size(du_t,3),size(du_t,2) ):: les_strain_rate
    real, dimension(size(du_t,3),size(du_t,3)):: s
    integer dim, ngi, gi

    ngi=size(du_t,2)
    dim=size(du_t,3)

    do gi=1, ngi

       s=0.5*matmul( nu, du_t(:,gi,:) )
       les_strain_rate(:,:,gi)=s+transpose(s)

    end do

  end function les_strain_rate

  function les_viscosity_strength(du_t, relu)
    !! Computes the strain rate modulus for the LES model 
    !! derivative of velocity shape function (nloc x ngi x dim)
    real, dimension(:,:,:), intent(in):: du_t
    !! relative velocity (nonl. vel.- grid vel.) (dim x nloc)
    real, dimension(:,:), intent(in):: relu

    real, dimension( size(du_t,2) ):: les_viscosity_strength

    real, dimension(size(du_t,3),size(du_t,3)):: s
    real vis
    integer dim, ngi, gi

    ngi=size(du_t,2)
    dim=size(du_t,3)

    do gi=1, ngi

       s=0.5*matmul( relu, du_t(:,gi,:) )
       s=s+transpose(s)
       ! Calculate modulus of strain rate
       vis=sqrt( 2*sum( s**2 ) )

       les_viscosity_strength(gi)=vis

    end do

  end function les_viscosity_strength

  function wale_viscosity_strength(du_t, relu)
    !! Computes the traceless symmetric part of the square of
    !! the resolved velocity gradient tensor for the LES model
    !! See a WALE paper for more (G_{ij})
    !! derivative of velocity shape function (nloc x ngi x dim)
    real, dimension(:,:,:), intent(in):: du_t
    !! relative velocity (nonl. vel.- grid vel.) (dim x nloc)
    real, dimension(:,:), intent(in):: relu

    real, dimension( size(du_t,2) ):: wale_viscosity_strength

    real, dimension(size(du_t,3),size(du_t,3)):: s, g
    real vis
    integer dim, ngi, gi, i

    ngi=size(du_t,2)
    dim=size(du_t,3)

    do gi=1, ngi

       s=matmul( relu, du_t(:,gi,:) )
       g=0.5*matmul(s,s)
       g=g+transpose(g)
       forall(i=1:dim) g(i,i)=0.
       
       vis=sqrt( 2*sum( g**2 ) )

       wale_viscosity_strength(gi)=vis

    end do

  end function wale_viscosity_strength

end module les_viscosity_module
