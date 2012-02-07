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
    tfield => extract_tensor_field(state, "TensorFirstFilterWidth", stat)
    if(stat==0) then
       call zero(tfield)
    end if
    tfield => extract_tensor_field(state, "TensorSecondFilterWidth", stat)
    if(stat==0) then
       call zero(tfield)
    end if
    sfield => extract_scalar_field(state, "SmagorinskyCoefficient", stat)
    if(stat==0) then
       call zero(sfield)
    end if
    sfield => extract_scalar_field(state, "ScalarFirstFilterWidth", stat)
    if(stat==0) then
       call zero(sfield)
    end if
    sfield => extract_scalar_field(state, "ScalarSecondFilterWidth", stat)
    if(stat==0) then
       call zero(sfield)
    end if

  end subroutine les_zero_diagnostic_fields

  subroutine dynamic_les_init_fields(state, nu, nu_f1, nu_f2, leonard_tensor, strain_product)

    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout)  :: nu
    type(vector_field), pointer, intent(inout) :: nu_f1, nu_f2
    type(tensor_field), pointer, intent(inout) :: leonard_tensor, strain_product
    integer :: stat

    ! First Filtered velocity
    nu_f1 => extract_vector_field(state, "FirstFilteredVelocity", stat)
    if(stat/=0) then
      ewrite(2,*) "Allocating FirstFilteredVelocity on VelocityMesh"
      allocate(nu_f1)
      call allocate(nu_f1, nu%dim, nu%mesh, "FirstFilteredVelocity")
    end if
    call zero(nu_f1)
    ewrite_minmax(nu_f1)

    ! Second Filtered velocity
    nu_f2 => extract_vector_field(state, "SecondFilteredVelocity", stat)
    if(stat/=0) then
      ewrite(2,*) "Allocating SecondFilteredVelocity on VelocityMesh"
      allocate(nu_f2)
      call allocate(nu_f2, nu%dim, nu%mesh, "SecondFilteredVelocity")
    end if
    call zero(nu_f2)
    ewrite_minmax(nu_f2)

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
             & eddy_visc_gi, visc_gi, strain1_gi, strain2_gi, les_coef_gi, &
             & t1_width_gi, s1_width_gi, t2_width_gi, s2_width_gi)

    type(state_type), intent(inout)                             :: state
    type(vector_field), intent(inout)                           :: nu
    type(scalar_field), intent(in)                              :: density
    integer, intent(in)                                         :: ele
    real, dimension(ele_ngi(nu,ele)), intent(in)                :: detwei
    real, dimension(ele_ngi(nu,ele)), intent(in), optional      :: les_coef_gi, s1_width_gi, s2_width_gi
    real, dimension(nu%dim,nu%dim,ele_ngi(nu,ele)),intent(inout), optional &
    & :: eddy_visc_gi, visc_gi, strain1_gi, strain2_gi, t1_width_gi, t2_width_gi
    real, dimension(ele_ngi(nu,ele))                            :: density_gi
    type(tensor_field), pointer                                 :: tensorfield
    type(scalar_field), pointer                                 :: scalarfield
    real, dimension(nu%dim,nu%dim,ele_loc(nu,ele))              :: tensor_loc
    real, dimension(ele_loc(nu,ele))                            :: scalar_loc, lumped_mass
    real, dimension(ele_loc(nu,ele),ele_loc(nu,ele))            :: mass_matrix
    type(element_type)                                          :: nu_shape
    type(patch_type)                                            :: patch, patch2
    real, dimension(ele_loc(nu,ele))                            :: lnodes
    integer                                                     :: stat, i, j, n, n2

    nu_shape = ele_shape(nu, ele)
    density_gi = ele_val_at_quad(density, ele)
    mass_matrix = shape_shape(nu_shape, nu_shape, detwei*density_gi)
    lumped_mass = sum(mass_matrix, 2)

    ! Eddy viscosity field m_ij is available with all LES models
    if (present(eddy_visc_gi)) then
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

    ! First-filtered strain rate field S1
    if (present(strain1_gi)) then
      tensorfield => extract_tensor_field(state, "FirstFilteredStrainRate", stat)
      if(stat==0) then
        tensor_loc=shape_tensor_rhs(nu_shape, strain1_gi, detwei)
        ! Divide by lumped mass
        do i=1, ele_loc(nu,ele)
          tensor_loc(:,:,i)=tensor_loc(:,:,i)/lumped_mass(i)
        end do
        call addto(tensorfield, ele_nodes(nu, ele), tensor_loc)
      end if
    end if

    ! Second-filtered strain rate field S2
    if (present(strain2_gi)) then
      tensorfield => extract_tensor_field(state, "SecondFilteredStrainRate", stat)
      if(stat==0) then
        tensor_loc=shape_tensor_rhs(nu_shape, strain2_gi, detwei)
        ! Divide by lumped mass
        do i=1, ele_loc(nu,ele)
          tensor_loc(:,:,i)=tensor_loc(:,:,i)/lumped_mass(i)
        end do
        call addto(tensorfield, ele_nodes(nu, ele), tensor_loc)
      end if
    end if

    ! Tensor first filter width
    if (present(t1_width_gi)) then
      tensorfield => extract_tensor_field(state, "TensorFirstFilterWidth", stat)
      if(stat==0) then
        !ewrite(2,*) "tensor first width ", t1_width_gi
        tensor_loc=shape_tensor_rhs(nu_shape, t1_width_gi, detwei)
        lnodes = ele_nodes(nu,ele)
        do i=1, size(lnodes)
          j=lnodes(i)
          patch = get_patch_ele(nu%mesh, j)
          n = patch%count
          tensor_loc(:,:,i)=tensor_loc(:,:,i)/lumped_mass(i)/n
        end do
        call addto(tensorfield, ele_nodes(nu, ele), tensor_loc)
      end if
    end if

    ! Scalar first filter width
    ! This should be an average of the widths of elements adjacent to node.
    if (present(s1_width_gi)) then
      scalarfield => extract_scalar_field(state, "ScalarFirstFilterWidth",stat)
      if(stat==0) then
        !ewrite(2,*) "scalar first width ", s1_width_gi
        scalar_loc=shape_rhs(nu_shape, s1_width_gi*detwei)
        ! Divide by lumped mass and no. of elements around node
        lnodes = ele_nodes(nu,ele)
        do i=1, size(lnodes)
          j=lnodes(i)
          patch = get_patch_ele(nu%mesh, j)
          n = patch%count
          scalar_loc(i)=scalar_loc(i)/lumped_mass(i)/n
        end do
        call addto(scalarfield, ele_nodes(nu, ele), scalar_loc)
      end if
    end if

    ! Tensor second filter width
    if (present(t2_width_gi)) then
      tensorfield => extract_tensor_field(state, "TensorSecondFilterWidth", stat)
      if(stat==0) then
        !ewrite(2,*) "tensor second width ", t2_width_gi
        tensor_loc=shape_tensor_rhs(nu_shape, t2_width_gi, detwei)
        lnodes = ele_nodes(nu,ele)
        do i=1, size(lnodes)
          j=lnodes(i)
          patch = get_patch_ele(nu%mesh, j)
          n = patch%count
          tensor_loc(:,:,i)=tensor_loc(:,:,i)/lumped_mass(i)/n
        end do
        call addto(tensorfield, ele_nodes(nu, ele), tensor_loc)
      end if
    end if

    ! Scalar second filter width
    if (present(s2_width_gi)) then
      scalarfield => extract_scalar_field(state, "ScalarSecondFilterWidth",stat)
      if(stat==0) then
        !ewrite(2,*) "scalar second width ", s2_width_gi
        scalar_loc=shape_rhs(nu_shape, s2_width_gi*detwei)
        lnodes = ele_nodes(nu,ele)
        do i=1, size(lnodes)
          j=lnodes(i)
          patch = get_patch_ele(nu%mesh, j)
          n = patch%count
          scalar_loc(i)=scalar_loc(i)/lumped_mass(i)/n
        end do
        call addto(scalarfield, ele_nodes(nu, ele), scalar_loc)
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
      end if
    end if

  end subroutine les_set_diagnostic_fields

  subroutine dynamic_les_filtered_fields(state, u, nu, positions, nu_f1, nu_f2, leonard, &
             &strain_prod, alpha, gamma, les_option_path, have_anisotropy)

    type(state_type), intent(inout)           :: state
    ! Unfiltered velocity
    type(vector_field), intent(inout):: u
    type(vector_field), pointer, intent(inout):: nu
    type(vector_field), intent(in)            :: positions
    ! First and second filtered velocities
    type(vector_field), pointer, intent(inout):: nu_f1, nu_f2
    ! Leonard tensor and strain product fields
    type(tensor_field), pointer, intent(inout):: leonard, strain_prod
    ! Scale factors
    real, intent(in)                          :: alpha, gamma
    character(len=OPTION_PATH_LEN), intent(in):: les_option_path
    logical, intent(in)                       :: have_anisotropy
    ! Local quantities
    type(tensor_field), pointer               :: tensorfield
    character(len=OPTION_PATH_LEN)            :: lpath
    integer                                   :: i,gi
    real, dimension(:), allocatable           :: u_loc
    real, dimension(:,:), allocatable         :: t_loc

    real, dimension(ele_loc(nu,1), ele_ngi(nu,1), nu%dim) :: du_t
    real, dimension(ele_ngi(nu,1))                 :: detwei
    real, dimension(nu%dim, nu%dim, ele_ngi(nu,1)) :: strain_gi, t_strain_gi, strain_prod_gi, t1_width_gi, t2_width_gi
    real, dimension(ele_ngi(nu,1))                 :: strain_mod, t_strain_mod
    type(element_type)                             :: shape_nu
    integer, dimension(:), pointer                 :: nodes_nu

    ! Path to solver options
    lpath = (trim(les_option_path)//"/dynamic_les")
    ewrite(2,*) "filter width ratio alpha: ", alpha
    ewrite(2,*) "filter width ratio gamma: ", gamma
    ewrite(2,*) "path to solver options: ", trim(lpath)
    ewrite(2,*) "anisotropic filtering and viscosity: ", have_anisotropy

    ! Filter the nonlinear velocity at both filter levels sequentially
    if(have_anisotropy) then
      call set(nu_f1, u)
      !call anisotropic_smooth_vector(nu, positions, nu_f1, alpha, lpath)
      call anisotropic_smooth_vector(nu_f1, positions, nu_f2, gamma, lpath)
    else
      call set(nu_f1, u)
      !call smooth_vector(nu, positions, nu_f1, alpha, lpath)
      call smooth_vector(nu_f1, positions, nu_f2, gamma, lpath)
    end if

    ewrite_minmax(nu)
    ewrite_minmax(nu_f1)
    ewrite_minmax(nu_f2)

    ! Local variables
    allocate(tensorfield)
    call allocate(tensorfield, nu%mesh, "TensorField")
    call zero(tensorfield)
    allocate(u_loc(nu%dim)); allocate(t_loc(nu%dim, nu%dim))
    u_loc=0.0; t_loc=0.0

    ! Get first-filtered velocity cross-product ui^f1*uj^f1
    do i=1, node_count(nu)
      u_loc = node_val(nu_f1,i)
      t_loc = outer_product(u_loc, u_loc)
      call set( tensorfield, i, t_loc )
    end do

    ! Calculate second-filtered (first-filtered velocity cross-product): (ui^f1*uj^f1)^f2
    if(have_anisotropy) then
      call anisotropic_smooth_tensor(tensorfield, positions, leonard, alpha, lpath)
    else
      call smooth_tensor(tensorfield, positions, leonard, alpha, lpath)
    end if
    call zero(tensorfield)

    ! Calculate (second-filtered velocity) cross-product: (ui^f1f2*uj^f1f2)
    do i=1, node_count(nu)
      u_loc = node_val(nu_f2,i)
      t_loc = outer_product(u_loc, u_loc)
      call set( tensorfield, i, t_loc )
    end do

    ! Leonard tensor field: (ui^f1*uj^f1)^f2 - (ui^f1f2*uj^f1f2)
    call addto( leonard, tensorfield, -1.0 )
    call zero(tensorfield)
    ewrite_minmax(leonard)

    ! Strain rate S1
    do i=1, element_count(nu)
      shape_nu = ele_shape(nu, i)

      ! Assuming no stabilisation scheme is used with LES so I can use velocity shape. (Tejada-Martinez thinks otherwise)
      call transform_to_physical(positions, i, shape_nu, dshape=du_t, detwei=detwei)

      ! Get strain S1 for first-filtered velocity (dim,dim,ngi)
      strain_gi = les_strain_rate(du_t, ele_val(nu_f1, i))
      ! Get strain S2 for second-filtered velocity (dim,dim,ngi)
      t_strain_gi = les_strain_rate(du_t, ele_val(nu_f2, i))
      ! Get strain modulus |S1| for first-filtered velocity (ngi)
      strain_mod = les_viscosity_strength(du_t, ele_val(nu_f1, i))
      ! Get strain modulus |S2| for second-filtered velocity (ngi)
      t_strain_mod = les_viscosity_strength(du_t, ele_val(nu_f2, i))

      ! Strain product |S1|S1
      do gi=1, ele_ngi(nu, i)
        strain_prod_gi(:,:,gi) = strain_gi(:,:,gi)*strain_mod(gi)
      end do
      call addto(tensorfield, ele_nodes(nu,i), shape_tensor_rhs(ele_shape(nu,i), strain_prod_gi, detwei))
    end do
    ewrite_minmax(tensorfield)

    ! Filter the strain product with second filter: (|S1|S1)^f2 (not a diagnostic field yet)
    if(have_anisotropy) then
      call anisotropic_smooth_tensor(tensorfield, positions, strain_prod, alpha, lpath)
    else
      call smooth_tensor(tensorfield, positions, strain_prod, alpha, lpath)
    end if
    ewrite_minmax(strain_prod)

    ewrite(2,*) "Finished setting LES filtered fields"

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
