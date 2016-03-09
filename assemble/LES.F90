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

module les_module
  !!< This module contains several subroutines and functions used to implement LES models
  use fldebug
  use spud
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use sparse_tools
  use vector_tools
  use fetools
  use fields
  use state_module
  use field_options
  use solvers
  use smoothing_module
  use state_fields_module, only: get_lumped_mass_on_submesh, get_lumped_mass,&
               get_mass_matrix
  implicit none

  private

  public les_viscosity_strength, wale_viscosity_strength
  public les_init_diagnostic_fields, les_assemble_diagnostic_fields, les_solve_diagnostic_fields, &
         leonard_tensor, les_strain_rate

contains

  subroutine les_init_diagnostic_fields(state, have_eddy_visc, have_filter_width, have_coeff)

    ! Arguments
    type(state_type), intent(inout)             :: state
    logical, intent(in)                         :: have_eddy_visc, have_filter_width, have_coeff
    
    ! Local variables
    logical, dimension(2)                       :: have_diagnostic_tfield
    logical, dimension(1)                       :: have_diagnostic_sfield
    character(len=FIELD_NAME_LEN), dimension(2) :: diagnostic_tfield_names
    character(len=FIELD_NAME_LEN), dimension(1) :: diagnostic_sfield_names
    type(tensor_field), pointer                 :: tfield
    type(scalar_field), pointer                 :: sfield
    integer                                     :: i

    ewrite(2,*) "Initialising optional LES diagnostic fields"
    
    have_diagnostic_tfield = (/have_eddy_visc, have_filter_width/)
    diagnostic_tfield_names(1) = "EddyViscosity"
    diagnostic_tfield_names(2) = "FilterWidth"
    
    diagnostic_tfield_loop: do i = 1, size(diagnostic_tfield_names)
      if(have_diagnostic_tfield(i)) then
         tfield => extract_tensor_field(state, diagnostic_tfield_names(i))
         call zero(tfield)
      end if
    end do diagnostic_tfield_loop

    have_diagnostic_sfield = (/have_coeff/)
    diagnostic_sfield_names(1) = "SmagorinskyCoefficient"

    diagnostic_sfield_loop: do i = 1, size(diagnostic_sfield_names)
      if(have_diagnostic_sfield(i)) then
         sfield => extract_scalar_field(state, diagnostic_sfield_names(i))
         call zero(sfield)
      end if
    end do diagnostic_sfield_loop

  end subroutine les_init_diagnostic_fields

  subroutine les_assemble_diagnostic_fields(state, nu, ele, detwei, &
                 mesh_size_gi,les_tensor_gi, les_coef_gi, &
                 have_eddy_visc, have_filter_width, have_coeff)

    ! Arguments
    type(state_type), intent(inout)                             :: state
    type(vector_field), intent(in)                              :: nu
    integer, intent(in)                                         :: ele
    real, dimension(ele_ngi(nu,ele)), intent(in)                :: les_coef_gi, detwei
    real, dimension(nu%dim,nu%dim,ele_ngi(nu,ele)),intent(in)   :: mesh_size_gi, les_tensor_gi
    logical, intent(in) :: have_eddy_visc, have_filter_width, have_coeff
    
    ! Local variables
    type(tensor_field), pointer                                 :: tfield
    type(scalar_field), pointer                                 :: sfield
    real, dimension(nu%dim,nu%dim,ele_loc(nu,ele))              :: tensor_loc
    real, dimension(ele_loc(nu,ele))                            :: scalar_loc

    ! Eddy viscosity
    if(have_eddy_visc) then
      tfield => extract_tensor_field(state, "EddyViscosity")
      tensor_loc=shape_tensor_rhs(ele_shape(nu, ele), les_tensor_gi, detwei)
      call addto(tfield, ele_nodes(nu, ele), tensor_loc)
    end if

    ! Filter width
    if(have_filter_width) then
      tfield => extract_tensor_field(state, "FilterWidth")
      tensor_loc=shape_tensor_rhs(ele_shape(nu, ele), mesh_size_gi, detwei)
      call addto(tfield, ele_nodes(nu, ele), tensor_loc)
    end if

    ! Smagorinsky Coefficient
    if(have_coeff) then
      sfield => extract_scalar_field(state, "SmagorinskyCoefficient")
      scalar_loc=shape_rhs(ele_shape(nu, ele), les_coef_gi*detwei)
      call addto(sfield, ele_nodes(nu, ele), scalar_loc)
    end if

  end subroutine les_assemble_diagnostic_fields

  subroutine les_solve_diagnostic_fields(state, have_eddy_visc, have_filter_width, have_coeff)

    ! Arguments
    type(state_type), intent(inout) :: state
    logical, intent(in) :: have_eddy_visc, have_filter_width, have_coeff
    
    ! Local variables
    logical, dimension(2)                       :: have_diagnostic_tfield
    logical, dimension(1)                       :: have_diagnostic_sfield
    character(len=FIELD_NAME_LEN), dimension(2) :: diagnostic_tfield_names
    character(len=FIELD_NAME_LEN), dimension(1) :: diagnostic_sfield_names
    type(tensor_field), pointer                 :: tfield
    type(scalar_field), pointer                 :: sfield
    integer                                     :: i
    type(vector_field), pointer                 :: u
    type(csr_matrix), pointer                   :: mass_matrix
    type(scalar_field), pointer                 :: lumped_mass
    type(scalar_field)                          :: inv_lumped_mass
    logical                                     :: lump_mass = .false.
    logical                                     :: use_submesh = .false.
    
    ewrite(2,*) "Solving for optional LES diagnostic fields"
        
    u => extract_vector_field(state, "Velocity")
    
    have_diagnostic_tfield = (/have_eddy_visc, have_filter_width/)
    diagnostic_tfield_names(1) = "EddyViscosity"
    diagnostic_tfield_names(2) = "FilterWidth"
    
    diagnostic_tfield_loop: do i = 1, size(diagnostic_tfield_names)
      if(have_diagnostic_tfield(i)) then
         tfield => extract_tensor_field(state, diagnostic_tfield_names(i))
         lump_mass = have_option(trim(tfield%option_path)//"/diagnostic/mass_matrix"//&
            &"/use_lumped_mass_matrix")
         use_submesh = have_option(trim(tfield%option_path)//"/diagnostic/mass_matrix"//&
            &"/use_lumped_mass_matrix/use_submesh") ! For P2 meshes.
            
         if(lump_mass) then
            if(use_submesh) then
               lumped_mass => get_lumped_mass_on_submesh(state, tfield%mesh)
            else
               lumped_mass => get_lumped_mass(state, tfield%mesh)
            end if
            call allocate(inv_lumped_mass, tfield%mesh)
            call invert(lumped_mass, inv_lumped_mass)
            call scale(tfield, inv_lumped_mass)
            call deallocate(inv_lumped_mass)
         else
            mass_matrix => get_mass_matrix(state, tfield%mesh)
            call petsc_solve(tfield, mass_matrix, tfield, option_path=u%option_path)
         end if
      end if
    end do diagnostic_tfield_loop
    
    have_diagnostic_sfield = (/have_coeff/)
    diagnostic_sfield_names(1) = "SmagorinskyCoefficient"
    
    diagnostic_sfield_loop: do i = 1, size(diagnostic_sfield_names)
      if(have_diagnostic_sfield(i)) then
         sfield => extract_scalar_field(state, diagnostic_sfield_names(i))
         lump_mass = have_option(trim(sfield%option_path)//"/diagnostic/mass_matrix"//&
            &"/use_lumped_mass_matrix")
         use_submesh = have_option(trim(sfield%option_path)//"/diagnostic/mass_matrix"//&
            &"/use_lumped_mass_matrix/use_submesh") ! For P2 meshes.
            
         if(lump_mass) then
            if(use_submesh) then
               lumped_mass => get_lumped_mass_on_submesh(state, sfield%mesh)
            else
               lumped_mass => get_lumped_mass(state, sfield%mesh)
            end if
            call allocate(inv_lumped_mass, sfield%mesh)
            call invert(lumped_mass, inv_lumped_mass)
            call scale(sfield, inv_lumped_mass)
            call deallocate(inv_lumped_mass)
         else
            mass_matrix => get_mass_matrix(state, sfield%mesh)
            call petsc_solve(sfield, mass_matrix, sfield, option_path=u%option_path)
         end if
      end if
    end do diagnostic_sfield_loop

  end subroutine les_solve_diagnostic_fields
  
  subroutine leonard_tensor(nu, positions, fnu, tnu, leonard, strainprod, alpha, gamma, path)

    ! Unfiltered velocity
    type(vector_field), pointer                           :: nu
    type(vector_field), intent(in)                        :: positions
    ! Filtered velocities
    type(vector_field), pointer                           :: fnu, tnu
    ! Leonard tensor and strain product
    type(tensor_field), pointer                           :: leonard, strainprod
    ! Scale factors
    real, intent(in)                                      :: alpha, gamma
    character(len=OPTION_PATH_LEN), intent(in)            :: path
    ! Local quantities
    type(tensor_field), pointer                           :: ui_uj, tui_tuj
    character(len=OPTION_PATH_LEN)                        :: lpath
    integer                                               :: i, ele, gi
    real, dimension(:), allocatable                       :: u_loc
    real, dimension(:,:), allocatable                     :: t_loc
    real, dimension(ele_loc(nu,1), ele_ngi(nu,1), nu%dim) :: du_t
    real, dimension(ele_ngi(nu,1))                        :: detwei
    real, dimension(nu%dim, nu%dim, ele_ngi(nu,1))        :: strain_gi, strain_prod_gi
    type(element_type)                                    :: shape_nu

    ! Path is to level above solver options
    lpath = (trim(path)//"/dynamic_les")
    ewrite(2,*) "filter factor alpha: ", alpha
    ewrite(2,*) "filter factor gamma: ", gamma

    ! First filter operator returns u^f:
    call anisotropic_smooth_vector(nu, positions, fnu, alpha, lpath)
    ! Test filter operator needs the ratio of test filter to mesh size and returns u^ft:
    call anisotropic_smooth_vector(fnu, positions, tnu, alpha*gamma, lpath)
    ewrite_minmax(nu)
    ewrite_minmax(fnu)
    ewrite_minmax(tnu)

    ! Velocity products (ui*uj)
    allocate(ui_uj); allocate(tui_tuj)
    call allocate(ui_uj, nu%mesh, "NonlinearVelocityProduct")
    call allocate(tui_tuj, nu%mesh, "TestNonlinearVelocityProduct")
    call zero(ui_uj); call zero(tui_tuj)

    ! Other local variables
    allocate(u_loc(nu%dim)); allocate(t_loc(nu%dim, nu%dim))
    u_loc=0.0; t_loc=0.0

    ! Get cross products of velocities
    do i=1, node_count(nu)
      u_loc = node_val(fnu,i)
      t_loc = outer_product(u_loc, u_loc)
      call set( ui_uj, i, t_loc )
      u_loc = node_val(tnu,i)
      ! Calculate (test-filtered velocity) products: (ui^ft*uj^ft)
      t_loc = outer_product(u_loc, u_loc)
      call set( tui_tuj, i, t_loc )
    end do

    ! Calculate test-filtered (velocity products): (ui^f*uj^f)^t
    call anisotropic_smooth_tensor(ui_uj, positions, leonard, alpha*gamma, lpath)

    ! Leonard tensor field
    call addto( leonard, tui_tuj, -1.0 )

    ! Zero tensor field for reuse in strain product assembly
    call zero(ui_uj)

    do i=1, element_count(nu)
      shape_nu = ele_shape(nu, i)
      ! Assuming no FE stabilisation is used with LES so we can use velocity shape.
      call transform_to_physical(positions, i, shape_nu, dshape=du_t, detwei=detwei)
      ! Strain rate of first filtered velocity S1^f
      strain_gi = les_strain_rate(du_t, ele_val(fnu, i))
      do gi=1, ele_ngi(nu, ele)
        ! Strain product = strain modulus*strain rate: |S1^f|S1^f
        strain_prod_gi(:,:,gi) = sqrt(2*sum(strain_gi(:,:,gi)*strain_gi(:,:,gi))) * strain_gi(:,:,gi)
      end do
      ! Assemble local tensor field
      call addto(ui_uj, ele_nodes(nu,i), shape_tensor_rhs(ele_shape(nu,i), strain_prod_gi, detwei))
    end do

    ! Filter strain product with test filter: (|S1^f|S1^f)^t
    call anisotropic_smooth_tensor(ui_uj, positions, strainprod, alpha*gamma, lpath)

    ! Deallocates
    deallocate(u_loc, t_loc)
    call deallocate(ui_uj)
    call deallocate(tui_tuj)
    deallocate(ui_uj); deallocate(tui_tuj)

  end subroutine leonard_tensor

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

end module les_module
