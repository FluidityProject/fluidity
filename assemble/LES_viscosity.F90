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
  !!< This module computes a viscosity term to implement LES
  use state_module
  use fields
  use field_options
  use spud
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use smoothing_module
  use vector_tools
  !use diagnostic_source_fields
  implicit none

  private

  public les_viscosity_strength, wale_viscosity_strength
  public dynamic_les_init_fields, leonard_tensor, les_strain_rate
  !public les_viscosity_module_register_diagnostic

contains

  subroutine dynamic_les_init_fields(state, les_option_path, tnu, mnu, nu_av, tnu_av, mnu_av, leonard, leonard_av, dynamic_les_coef, dynamic_eddy_visc, dynamic_strain, dynamic_t_strain, dynamic_strain_mod, dynamic_t_strain_mod, dynamic_filter)

    type(state_type), intent(inout) :: state
    type(vector_field), pointer :: nu_av, tnu_av, mnu_av, tnu, mnu, u
    type(tensor_field), pointer :: leonard, leonard_av, dynamic_les_coef, dynamic_eddy_visc, dynamic_strain, dynamic_t_strain, dynamic_filter
    type(scalar_field), pointer :: dynamic_strain_mod, dynamic_t_strain_mod
    character(len=OPTION_PATH_LEN) :: les_option_path
    integer :: stat

    u => extract_vector_field(state, "Velocity")

    ewrite(2,*) "Initialising compulsory dynamic LES fields"
    ! Test-filtered velocity field
    tnu => extract_vector_field(state, "DynamicFilteredVelocity", stat)
    if(stat == 0) then
      ewrite(2,*) "zeroing field: ", trim(tnu%name), ", ", trim(tnu%mesh%name)
      call zero(tnu)
    end if
    ! Dynamic velocity field
    mnu => extract_vector_field(state, "DynamicVelocity", stat)
    if(stat == 0) then
      ewrite(2,*) "zeroing field: ", trim(mnu%name), ", ", trim(mnu%mesh%name)
      call zero(mnu)
    end if
    ! Leonard tensor field L_ij
    leonard => extract_tensor_field(state, "DynamicLeonardTensor", stat)
    if(stat == 0) then
      ewrite(2,*) "zeroing field: ", trim(leonard%name)
      call zero(leonard)
    end if

    ! Are we using averaged velocity to stabilise the model?
    if(have_option(trim(les_option_path)//"/dynamic_les/enable_averaging")) then
      ewrite(2,*) "Initialising dynamic LES averaged fields"
      ! Averaged velocity field
      !nu_av => vector_source_field(state, u)
      nu_av => u

      ! Test-filtered velocity field
      tnu_av => extract_vector_field(state, "DynamicFilteredAverageVelocity", stat)
      if(stat == 0) then
        ewrite(2,*) "zeroing field: ", trim(tnu_av%name), ", ", trim(tnu_av%mesh%name)
        call zero(tnu_av)
      end if
      ! Dynamic velocity field
      mnu_av => extract_vector_field(state, "DynamicAverageVelocity", stat)
      if(stat == 0) then
        ewrite(2,*) "zeroing field: ", trim(mnu_av%name), ", ", trim(mnu_av%mesh%name)
        call zero(mnu_av)
      end if
      ! Leonard tensor field L_ij
      leonard_av => extract_tensor_field(state, "DynamicAverageLeonardTensor", stat)
      if(stat == 0) then
        ewrite(2,*) "zeroing field: ", trim(leonard_av%name)
        call zero(leonard_av)
      end if
    else
      nu_av => null(); mnu_av => null(); tnu_av => null(); leonard_av => null()
    end if

    ewrite(2,*) "Initialising optional dynamic LES diagnostic fields"
    ! Filter width
    if(have_option(trim(les_option_path)//"/dynamic_les/tensor_field::DynamicFilterWidth")) then
      dynamic_filter => extract_tensor_field(state, "DynamicFilterWidth", stat)
      if(stat == 0) then
        call zero(dynamic_filter)
      end if
    end if
    ! Strain rate field S1
    if(have_option(trim(les_option_path)//"/dynamic_les/tensor_field::DynamicStrainRate")) then
      dynamic_strain => extract_tensor_field(state, "DynamicStrainRate", stat)
      if(stat == 0) then
        call zero(dynamic_strain)
      end if
    end if
    ! Filtered strain rate field S2
    if(have_option(trim(les_option_path)//"/dynamic_les/tensor_field::DynamicFilteredStrainRate")) then
      dynamic_t_strain => extract_tensor_field(state, "DynamicFilteredStrainRate", stat)
      if(stat == 0) then
        call zero(dynamic_t_strain)
      end if
    end if
    ! Strain rate modulus field |S1|
    if(have_option(trim(les_option_path)//"/dynamic_les/scalar_field::DynamicStrainRateModulus")) then
      dynamic_strain_mod => extract_scalar_field(state, "DynamicStrainRateModulus", stat)
      if(stat == 0) then
        call zero(dynamic_strain_mod)
      end if
    end if
    ! Filtered strain rate modulus field |S2|
    if(have_option(trim(les_option_path)//"/dynamic_les/scalar_field::DynamicFilteredStrainRateModulus")) then
      dynamic_t_strain_mod => extract_scalar_field(state, "DynamicFilteredStrainRateModulus", stat)
      if(stat == 0) then
        call zero(dynamic_t_strain_mod)
      end if
    end if
    ! Eddy viscosity field m_ij
    if(have_option(trim(les_option_path)//"/dynamic_les/tensor_field::DynamicEddyViscosity")) then
      dynamic_eddy_visc => extract_tensor_field(state, "DynamicEddyViscosity", stat)
      if(stat == 0) then
        call zero(dynamic_eddy_visc)
      end if
    end if
    ! Dynamic Smagorinsky coefficient C
    if(have_option(trim(les_option_path)//"/dynamic_les/tensor_field::DynamicSmagorinskyCoefficient")) then
      dynamic_les_coef => extract_tensor_field(state, "DynamicSmagorinskyCoefficient", stat)
      if(stat == 0) then
        call zero(dynamic_les_coef)
      end if
    end if

  end subroutine dynamic_les_init_fields

  subroutine leonard_tensor(mnu, positions, tnu, leonard, alpha, path)

    ! Unfiltered velocity
    type(vector_field), pointer               :: mnu
    type(vector_field), intent(in)            :: positions
    ! Filtered velocity
    type(vector_field), pointer               :: tnu
    ! Leonard tensor field
    type(tensor_field), pointer               :: leonard
    ! Scale factor: test filter/mesh size
    real, intent(in)                          :: alpha
    character(len=OPTION_PATH_LEN), intent(in):: path
    ! Local quantities
    type(tensor_field), pointer               :: ui_uj, tui_tuj
    character(len=OPTION_PATH_LEN)            :: lpath
    integer                                   :: i, stat
    real, dimension(:), allocatable           :: u_loc
    real, dimension(:,:), allocatable         :: t_loc

    ! Path is to level above solver options
    lpath = (trim(path)//"/dynamic_les")
    ewrite(2,*) "path: ", trim(lpath)
    ewrite(2,*) "alpha: ", alpha

    do i = 1, positions%dim
      ewrite_minmax(mnu%val(i,:))
      ewrite_minmax(tnu%val(i,:))
    end do

    call anisotropic_smooth_vector(mnu, positions, tnu, alpha, lpath)

    do i = 1, positions%dim
      ewrite_minmax(mnu%val(i,:))
      ewrite_minmax(tnu%val(i,:))
    end do

    ! Velocity products (ui*uj)
    allocate(ui_uj); allocate(tui_tuj)
    call allocate(ui_uj, mnu%mesh, "NonlinearVelocityProduct")
    call allocate(tui_tuj, mnu%mesh, "TestNonlinearVelocityProduct")
    call zero(ui_uj); call zero(tui_tuj)

    ! Other local variables
    allocate(u_loc(mnu%dim)); allocate(t_loc(mnu%dim, mnu%dim))
    u_loc=0.0; t_loc=0.0

    ! Get cross products of velocities
    do i=1, node_count(mnu)
      ! Mesh filter ^r
      u_loc = node_val(mnu,i)
      call outer_product(u_loc, u_loc, t_loc)
      call set( ui_uj, i, t_loc )
      ! Test filter ^t
      u_loc = node_val(tnu,i)
      ! Calculate (test-filtered velocity) products: (ui^rt*uj^rt)
      call outer_product(u_loc, u_loc, t_loc)
      call set( tui_tuj, i, t_loc )
    end do

    ! Calculate test-filtered (velocity products): (ui^r*uj^r)^t
    call anisotropic_smooth_tensor(ui_uj, positions, leonard, alpha, lpath)

    ! Leonard tensor field
    call addto( leonard, tui_tuj, -1.0 )

    ! Deallocates
    deallocate(u_loc, t_loc)
    call deallocate(ui_uj)
    call deallocate(tui_tuj)
    deallocate(ui_uj); deallocate(tui_tuj)

  end subroutine leonard_tensor

  !subroutine les_viscosity_module_register_diagnostic

  !  dynamic_les_coef, dynamic_eddy_visc, dynamic_strain, dynamic_t_strain, dynamic_filter

  !  call register_diagnostic(dim=1, name="tensor", statistic="effectivestress", material_phase="Fluid")

  !  call set_diagnostic(name, statistic, material_phase, value)

  !end subroutine les_viscosity_module_register_diagnostic

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
