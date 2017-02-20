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
  use boundary_conditions
  use vector_tools
  use fetools
  use transform_elements
  use fields
  use fields_data_types
  use state_module
  use field_options
  use solvers
  use smoothing_module
  use state_fields_module, only: get_lumped_mass_on_submesh, get_lumped_mass,&
               get_mass_matrix
  implicit none

  private

  public les_viscosity_strength, wale_viscosity_strength
  public les_init_diagnostic_fields, les_assemble_diagnostic_fields, les_solve_diagnostic_fields
  public compute_les_local_fields, les_strain_rate, calculate_periodic_channel_forcing
  public calculate_dynamic_slip_coefficient
contains

  subroutine les_init_diagnostic_fields(state, have_eddy_visc, have_filter_width, have_coeff, have_sgs_tensor)

    ! Arguments
    type(state_type), intent(inout)             :: state
    logical, intent(in)                         :: have_eddy_visc, have_filter_width, have_coeff, have_sgs_tensor
    
    ! Local variables
    logical, dimension(3)                       :: have_diagnostic_tfield
    logical, dimension(1)                       :: have_diagnostic_sfield
    character(len=FIELD_NAME_LEN), dimension(3) :: diagnostic_tfield_names
    character(len=FIELD_NAME_LEN), dimension(1) :: diagnostic_sfield_names
    type(tensor_field), pointer                 :: tfield
    type(scalar_field), pointer                 :: sfield
    integer                                     :: i

    ewrite(2,*) "Initialising optional LES diagnostic fields"
    
    have_diagnostic_tfield = (/have_eddy_visc, have_filter_width, have_sgs_tensor/)
    diagnostic_tfield_names(1) = "EddyViscosity"
    diagnostic_tfield_names(2) = "FilterWidth"
    diagnostic_tfield_names(3) = "SGSTensor"
    
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
                 mesh_size_gi, les_tensor_gi, sgs_tensor_gi, les_coef_gi, &
                 have_eddy_visc, have_filter_width, have_coeff, have_sgs_tensor)

    ! Arguments
    type(state_type), intent(inout)                             :: state
    type(vector_field), intent(in)                              :: nu
    integer, intent(in)                                         :: ele
    real, dimension(ele_ngi(nu,ele)), intent(in)                :: les_coef_gi, detwei
    real, dimension(nu%dim,nu%dim,ele_ngi(nu,ele)),intent(in)   :: mesh_size_gi, les_tensor_gi, sgs_tensor_gi
    logical, intent(in) :: have_eddy_visc, have_filter_width, have_coeff, have_sgs_tensor
    
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
    if(have_sgs_tensor) then
      tfield => extract_tensor_field(state, "SGSTensor")
      tensor_loc=shape_tensor_rhs(ele_shape(nu, ele), sgs_tensor_gi, detwei)
      call addto(tfield, ele_nodes(nu, ele), tensor_loc)
    end if

    ! Smagorinsky Coefficient
    if(have_coeff) then
      sfield => extract_scalar_field(state, "SmagorinskyCoefficient")
      scalar_loc=shape_rhs(ele_shape(nu, ele), les_coef_gi*detwei)
      call addto(sfield, ele_nodes(nu, ele), scalar_loc)
    end if

    ! Filter width
    if(have_filter_width) then
      tfield => extract_tensor_field(state, "FilterWidth")
      tensor_loc=shape_tensor_rhs(ele_shape(nu, ele), mesh_size_gi, detwei)
      call addto(tfield, ele_nodes(nu, ele), tensor_loc)
    end if

  end subroutine les_assemble_diagnostic_fields

  subroutine les_solve_diagnostic_fields(state, have_eddy_visc, have_filter_width, have_coeff, have_sgs_tensor)

    ! Arguments
    type(state_type), intent(inout) :: state
    logical, intent(in) :: have_eddy_visc, have_filter_width, have_coeff, have_sgs_tensor
    
    ! Local variables
    logical, dimension(3)                       :: have_diagnostic_tfield
    logical, dimension(1)                       :: have_diagnostic_sfield
    character(len=FIELD_NAME_LEN), dimension(3) :: diagnostic_tfield_names
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
    
    have_diagnostic_tfield = (/have_eddy_visc, have_filter_width, have_sgs_tensor/)
    diagnostic_tfield_names(1) = "EddyViscosity"
    diagnostic_tfield_names(2) = "FilterWidth"
    diagnostic_tfield_names(3) = "SGSTensor"

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
  
  subroutine compute_les_local_fields(u, nu, positions, fnu, tnu, leonard, strainprod, sgstensor, alpha, gamma, length_scale_type, path, dynamic_les, exact_sgs)

    type(vector_field), intent(in)                        :: positions
    type(vector_field), intent(inout)                     :: u, nu, fnu, tnu
    ! Leonard tensor, strain product, SGS tensor
    type(tensor_field), intent(inout)                     :: leonard, strainprod, sgstensor
    ! Filter scale factors
    real, intent(inout)                                   :: alpha, gamma
    character(len=OPTION_PATH_LEN), intent(inout)         :: length_scale_type
    character(len=OPTION_PATH_LEN), intent(in)            :: path
    logical, intent(in)                                   :: dynamic_les, exact_sgs

    ! Local quantities
    type(tensor_field), pointer                           :: tfield
    character(len=OPTION_PATH_LEN)                        :: lpath, HOT
    logical                                               :: explicit_filter
    integer                                               :: i, ele, gi, loc
    real, dimension(:), allocatable                       :: u_loc
    real, dimension(:,:), allocatable                     :: t_loc
    real, dimension(ele_loc(nu,1), ele_ngi(nu,1), nu%dim) :: du_t
    real, dimension(ele_ngi(nu,1))                        :: detwei
    real, dimension(nu%dim, nu%dim, ele_ngi(nu,1))        :: strain_gi, strain_prod_gi
    real, dimension(nu%dim, nu%dim, ele_ngi(nu,1))        :: mat_gi
    real, dimension(nu%dim, nu%dim)                       :: mat
    real, dimension(nu%dim, ele_ngi(nu,1))                :: laplacian_gi
    real, dimension(ele_loc(nu,1))                        :: laplacian
    real                                                  :: delta
    type(element_type)                                    :: shape_nu

    if(dynamic_les) then

      ! Path is to level above solver options
      lpath = (trim(path)//"/dynamic_les")
      explicit_filter = .False.

    else if(exact_sgs) then

      lpath = (trim(path)//"/exact_sgs")
      ! have to initialise these values here as they're not used in the exact SGS tensor
      gamma = 2.0
      length_scale_type = "scalar"

      ! Check if 4th order term is to be included:
      ! Germano's exact SGS closure or the LANS-alpha model.
      ! Otherwise just 2nd-order term (Clark or Gradient Model).
      call get_option(trim(lpath)//"/HOT", HOT)

      ! Explicit or implicit filter
      explicit_filter = have_option(trim(lpath)//"/explicit_filter")

    end if

    ewrite(2,*) "filter factor alpha: ", alpha
    ewrite(2,*) "filter factor gamma: ", gamma
    ewrite(2,*) "explicit filter: ", explicit_filter
    ewrite(2,*) "filter width type: ", trim(length_scale_type)
    ewrite(2,*) "HOT: ", trim(HOT)

    if(explicit_filter) then
      if(length_scale_type=="scalar") then
        ! First filter operator returns u^f:
        call smooth_vector(nu, positions, fnu, alpha, lpath)
        ! Test filter operator needs the ratio of test filter to mesh size and returns u^ft:
        call smooth_vector(fnu, positions, tnu, alpha*gamma, lpath)
      else if(length_scale_type=="tensor") then
        ! First filter operator returns u^f:
        call anisotropic_smooth_vector(nu, positions, fnu, alpha, lpath)
        ! Test filter operator needs the ratio of test filter to mesh size and returns u^ft:
        call anisotropic_smooth_vector(fnu, positions, tnu, alpha*gamma, lpath)
      end if

      ! Must set nonlinear velocity to be filtered velocity
      call set(nu, fnu)

    else
      ! Implicit filter: use discretised velocity as first filtered velocity
      call set(fnu, nu)
      ! Test filter operator needs the ratio of test filter to mesh size and returns u^ft:
      if(length_scale_type=="scalar") then
        call smooth_vector(fnu, positions, tnu, gamma, lpath)
      else if(length_scale_type=="tensor") then
        call anisotropic_smooth_vector(fnu, positions, tnu, gamma, lpath)
      end if
    end if

    ewrite_minmax(nu)
    ewrite_minmax(fnu)
    ewrite_minmax(tnu)

    ! Velocity products (ui*uj)
    allocate(tfield)
    !allocate(tui_tuj)
    call allocate(tfield, nu%mesh, "TensorField")
    !call allocate(tui_tuj, nu%mesh, "TestNonlinearVelocityProduct")
    call zero(tfield)
    !call zero(tui_tuj)

    ! Other local variables
    allocate(u_loc(nu%dim)); allocate(t_loc(nu%dim, nu%dim))
    u_loc=0.0; t_loc=0.0

    ! Get cross products of velocities
    do i=1, node_count(nu)
      u_loc = node_val(fnu,i)
      t_loc = outer_product(u_loc, u_loc)
      call set( tfield, i, t_loc )
      !u_loc = node_val(tnu,i)
      ! Calculate (test-filtered velocity) products: (ui^ft*uj^ft)
      !t_loc = outer_product(u_loc, u_loc)
      !call set( tui_tuj, i, t_loc )
    end do

    ! Calculate test-filtered (velocity products): (ui^f*uj^f)^t
    if(length_scale_type=="scalar") then
      call smooth_tensor(tfield, positions, leonard, alpha*gamma, lpath)
    else if(length_scale_type=="tensor") then
      !call anisotropic_smooth_tensor(tfield, positions, leonard, alpha*gamma, lpath)
      call smooth_tensor(tfield, positions, leonard, alpha*gamma, lpath)
    end if

    ! Leonard tensor field
    !call addto( leonard, tui_tuj, -1.0 )

    if(dynamic_les) then

      ! Zero tensor field for reuse in strain product assembly
      call zero(tfield)

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
        call addto(tfield, ele_nodes(nu,i), shape_tensor_rhs(ele_shape(nu,i), strain_prod_gi, detwei))
      end do

      ! Filter strain product with test filter: (|S1^f|S1^f)^t
      if(length_scale_type=="scalar") then
        call smooth_tensor(tfield, positions, strainprod, alpha*gamma, lpath)
      else if(length_scale_type=="tensor") then
        !call anisotropic_smooth_tensor(tfield, positions, strainprod, alpha*gamma, lpath)
        call smooth_tensor(tfield, positions, strainprod, alpha*gamma, lpath)
      end if

    else if(exact_sgs) then

      call zero(tfield)

      do ele=1, element_count(nu)
        shape_nu = ele_shape(nu, ele)
        call transform_to_physical(positions, ele, shape_nu, dshape=du_t, detwei=detwei)

        ! Filter width - stick with scalar filter width for simplicity
        if(length_scale_type=="scalar") then
          delta = alpha**2/24.*length_scale_scalar(positions, ele)
        else if(length_scale_type=="tensor") then
          delta = alpha**2/24.*length_scale_scalar(positions, ele)
          !delta = alpha**2/24.*length_scale_tensor(du_t, shape_nu)
        end if

        select case(HOT)

        ! No model
        case("control")
          do gi = 1, ele_ngi(nu, ele)
            mat_gi(:,:,gi) = 0.0
          end do

        ! Clark or Gradient model
        case("none")
          do gi = 1, ele_ngi(nu, ele)
            ! velocity gradient
            mat = matmul(ele_val(fnu, ele), du_t(:,gi,:))
            ! 1 gradient product term: (du_i/dx_k) (du_j/dx_k)
            mat_gi(:,:,gi) = 2.0*delta*(matmul(mat, transpose(mat)))
          end do

        ! Germano's closure
        case("exact")
          do gi = 1, ele_ngi(nu, ele)
            ! velocity gradient
            mat = matmul(ele_val(fnu, ele), du_t(:,gi,:))
            ! gradient product: 2 (du_i/dx_k) (du_j/dx_k)
            mat_gi(:,:,gi) = 2.0*delta*matmul(mat, transpose(mat))

            ! Laplacian operator: dim.dim = dim
            do loc = 1, ele_loc(fnu, ele)
              laplacian(loc) = dot_product(du_t(loc,gi,:), du_t(loc,gi,:))
            end do

            ! Laplacian of filtered velocity: (dim*loc)*loc = dim
            laplacian_gi(:,gi) = matmul(ele_val(fnu, ele), laplacian)

            ! Laplacian product
            mat_gi(:,:,gi) = mat_gi(:,:,gi) + delta**4*outer_product(laplacian_gi(:,gi), laplacian_gi(:,gi))
          end do

        ! LANS-alpha closure
        case ("lans")
          do gi = 1, ele_ngi(nu, ele)
            ! velocity gradient
            mat = matmul(ele_val(fnu, ele), du_t(:,gi,:))
            ! 3 gradient product terms: (du_i/dx_k) (du_j/dx_k) - (du_k/dx_i) (du_k/dx_j) + (du_i/dx_k) (du_k/dx_j)
            mat_gi(:,:,gi) = delta*(matmul(mat, transpose(mat)) - matmul(transpose(mat), mat) + matmul(mat, mat))
          end do
        end select

        ! add to tensor RHS field
        call addto(tfield, ele_nodes(nu, ele), shape_tensor_rhs(shape_nu, mat_gi, detwei))

      end do

      ! Filter tensor RHS field.
      ! If calculate_boundaries we need to set sgstensor = 0 where there are zero Dirichlet BCs on Velocity
      if(have_option(trim(lpath)//"/calculate_boundaries")) then
        if(length_scale_type=="scalar") then
          call smooth_tensor(tfield, positions, sgstensor, alpha, lpath, u)
        else if(length_scale_type=="tensor") then
          call smooth_tensor(tfield, positions, sgstensor, alpha, lpath, u)
        end if
      ! otherwise, there are no BCs to apply to sgstensor
      else
        if(length_scale_type=="scalar") then
          call smooth_tensor(tfield, positions, sgstensor, alpha, lpath)
        else if(length_scale_type=="tensor") then
          call smooth_tensor(tfield, positions, sgstensor, alpha, lpath)
        end if
      end if

    end if

    ! Deallocates
    deallocate(u_loc, t_loc)
    call deallocate(tfield)
    !call deallocate(tui_tuj)
    deallocate(tfield)
    !deallocate(tui_tuj)

  end subroutine compute_les_local_fields

  subroutine calculate_periodic_channel_forcing(state, oldu, nu, positions, density, source_field)

    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: source_field
    type(vector_field), intent(in) :: oldu, nu, positions
    type(scalar_field), intent(in) :: density
    
    character(len = OPTION_PATH_LEN) :: path
    integer :: stat
    type(scalar_field), pointer :: masslump
    integer :: i

    ewrite(2,*) "Calculating periodic channel forcing term"
    masslump => get_lumped_mass(state, source_field%mesh)

    call zero(source_field)
    do i = 1, ele_count(source_field)
      call assemble_periodic_channel_forcing_ele(i, positions, oldu, nu, density, source_field)
    end do
    
    do i = 1, source_field%dim
      source_field%val(i,:) = source_field%val(i,:) / masslump%val
    end do

  end subroutine calculate_periodic_channel_forcing

  subroutine assemble_periodic_channel_forcing_ele(ele, positions, oldu, nu, density, source_field)

    integer, intent(in) :: ele
    type(vector_field), intent(inout) :: source_field
    type(vector_field), intent(in) :: positions, oldu, nu
    type(scalar_field), intent(in) :: density

    real, dimension(ele_ngi(source_field, ele)) :: detwei
    real, dimension(source_field%dim, ele_loc(source_field, ele)) :: src

    !call transform_to_physical(positions, ele, detwei=detwei)

    ! Assume incompressible flow for simplicity
    !vol = element_volume(positions, ele)
    !src = 0.0

    !call addto(source_field, ele_nodes(source_field, ele), &
    !  & shape_vector_rhs(ele_shape(source_field, ele), src, detwei))

  end subroutine assemble_periodic_channel_forcing_ele

  function les_strain_rate(du_t, nu)
    !! Computes the strain rate
    !! derivative of velocity shape function (nloc x ngi x dim)
    real, dimension(:,:,:), intent(in):: du_t
    !! nonlinear velocity (dim x nloc)
    real, dimension(:,:), intent(in):: nu
    real, dimension( size(du_t,3),size(du_t,3),size(du_t,2) ):: les_strain_rate
    real, dimension(size(du_t,3),size(du_t,3)):: s
    integer :: dim, ngi, gi, i
    real :: trace

    ngi=size(du_t,2)
    dim=size(du_t,3)

    do gi=1, ngi

       s=0.5*matmul( nu, du_t(:,gi,:) )
       les_strain_rate(:,:,gi)=s+transpose(s)

       ! Subtract trace
       trace = 0.0
       do i=1, dim
         trace = trace + les_strain_rate(i,i,gi)
       end do
       trace = trace/3.0
       do i=1, dim
         les_strain_rate(i,i,gi) = les_strain_rate(i,i,gi) - trace
       end do

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

  function les_gradient_product(du_t, nu)
    !! Computes the cross product of velocity gradients

    !! derivative of velocity shape function (nloc x ngi x dim)
    real, dimension(:,:,:), intent(in):: du_t
    !! nonlinear velocity (dim x nloc)
    real, dimension(:,:), intent(in):: nu
    real, dimension( size(du_t,3),size(du_t,3),size(du_t,2) ):: les_gradient_product
    real, dimension(size(du_t,3),size(du_t,3)):: s
    integer dim, ngi, gi

    ngi=size(du_t,2)
    dim=size(du_t,3)

    do gi=1, ngi

       s=matmul( nu, du_t(:,gi,:) )
       les_gradient_product(:,:,gi)=2.0*matmul(s, transpose(s))

    end do

  end function les_gradient_product

  subroutine calculate_dynamic_slip_coefficient(x, u, fnu, tnu, leonard, strainprod, sgstensor, dyn_coeff, alpha, gamma, dynamic_les, exact_sgs)

    type(vector_field), intent(in) :: x, u

    ! Fields/vars for dynamic slip model
    type(scalar_field), intent(in)    :: dyn_coeff
    type(vector_field), intent(in)    :: fnu, tnu
    type(tensor_field), intent(in)    :: leonard, strainprod, sgstensor
    real, intent(in)                  :: alpha, gamma
    logical, intent(in)               :: dynamic_les, exact_sgs

    ! local vars
    type(vector_field), pointer     :: surface_field, surface_field_2
    real, dimension(face_loc(u, 1)) :: surface_field_addto
    real, dimension(u%dim)          :: slipsum
    character(len=OPTION_PATH_LEN)  :: bc_path_i, bc_type_path, bc_component_path
    character(len=FIELD_NAME_LEN)   :: bc_name, bc_type
    integer, dimension(:), pointer  :: surface_element_list, surface_node_list
    integer                         :: i, j, k, sele, nbcs, dim_set
    logical                         :: applies(3)
    character(len=20), dimension(3) :: aligned_components
    character(len=20), parameter, dimension(3) :: &
      cartesian_aligned_components=(/"x_component","y_component","z_component"/), &
      surface_aligned_components=(/"normal_component   ","tangent_component_1","tangent_component_2"/)

    ewrite(2,*) "Calculating dynamic_slip coefficient"

    ! Get number of boundary conditions
    nbcs=option_count(trim(u%option_path)//'/prognostic/boundary_conditions')

    ! Loop over boundary conditions
    do i=1, nbcs

      call get_boundary_condition(u, i, type=bc_type, name=bc_name, surface_element_list=surface_element_list, &
              surface_node_list=surface_node_list, applies=applies, option_path=bc_path_i)

      ! Skip partitions with no surface elements with this BC
      if(size(surface_element_list)==0) cycle

      ! Consider only dynamic slip BCs
      if(trim(bc_type)=='dynamic_slip') then

        if(have_option(trim(bc_path_i)//"/type[0]/align_bc_with_cartesian")) then
           aligned_components=cartesian_aligned_components
           bc_type_path=trim(bc_path_i)//"/type[0]/align_bc_with_cartesian"
        else
           aligned_components=surface_aligned_components             
           bc_type_path=trim(bc_path_i)//"/type[0]/align_bc_with_surface"
        end if

        surface_field   => extract_surface_field(u, bc_name, name="order_zero_coefficient")
        surface_field_2 => extract_surface_field(u, bc_name, name="order_one_coefficient")
        call zero(surface_field)
        call zero(surface_field_2)

        dim_set = 0
        slipsum = 0
        surface_field_addto = 0.0

        do j=1,3

          bc_component_path=trim(bc_type_path)//"/"//aligned_components(j)
          applies(j)=have_option(trim(bc_component_path))

          if (applies(j)) then

            ! Set zero order coeff at ele nodes
            do k=1, size(surface_node_list)
              call set(surface_field, j, k, 0.0)
            end do

            ! Calculate dynamic slip coefficient
            ! But only if we haven't done it already for another component
            if(dim_set==0) then

              do k=1, size(surface_element_list)
                sele = surface_element_list(k)

                call calculate_dynamic_slip_coefficient_face(x, u, fnu, tnu, &
                     leonard, strainprod, sgstensor, dyn_coeff, alpha, gamma, &
                     dynamic_les, exact_sgs, surface_field_addto, sele)

                ! add contribution to surface field
                call addto(surface_field_2, j, ele_nodes(surface_field_2, k), surface_field_addto)
              end do

              do k=1, size(surface_node_list)
                slipsum(j) = slipsum(j) + node_val(surface_field_2, j, k)
              end do

              ! flag this dimension so we don't have to recalculate the coeff
              dim_set = j
            else
              do k=1, size(surface_node_list)
                call set(surface_field_2, j, k, node_val(surface_field_2, dim_set, k))
              end do
            end if
          end if
        end do
        ewrite_minmax(surface_field)
        ewrite_minmax(surface_field_2)
        ewrite(2,*) "sum(slip coeff): ", slipsum
      end if
    end do

  end subroutine calculate_dynamic_slip_coefficient

  subroutine calculate_dynamic_slip_coefficient_face(x, u, fnu, tnu, leonard, strainprod, sgstensor, dyn_coeff, alpha, gamma, dynamic_les, exact_sgs, surface_field_addto, sele)

    type(scalar_field), intent(in)    :: dyn_coeff
    type(vector_field), intent(in)    :: x, u, fnu, tnu
    type(tensor_field), intent(in)    :: leonard, strainprod, sgstensor
    logical, intent(in)               :: dynamic_les, exact_sgs
    real, intent(in)                  :: alpha, gamma
    integer, intent(in)               :: sele
    real, dimension(face_loc(u, sele)), intent(out) :: surface_field_addto

    ! FE quantities
    integer :: i, j, iloc, gi, ngi, lface, ele
    real, dimension(face_ngi(u, sele))                                       :: detwei_bdy
    real, dimension(u%dim, face_ngi(u, sele))                                :: normal_bdy
    real, dimension(u%dim, u%dim, ele_ngi(u, face_ele(u, sele)))             :: invj
    real, dimension(u%dim, u%dim, face_ngi(u, sele))                         :: invj_face
    real, dimension(ele_loc(u, face_ele(u, sele)), face_ngi(u, sele), u%dim) :: dshape_face

    type(element_type)          :: augmented_shape
    type(element_type), pointer :: u_shape, fshape, x_shape, source_shape

    ! Local quantities specific to dynamic slip wall model
    real, dimension(u%dim, u%dim, face_ngi(u, sele))  :: leonard_gi, strain_gi, t_strain_gi, strainprod_gi
    real, dimension(u%dim, u%dim, face_ngi(u, sele))  :: tau_gi, tautest_gi
    real, dimension(u%dim, u%dim)                     :: mij, grad_fnu, grad_tnu
    real, dimension(u%dim, face_ngi(u, sele))         :: fnu_gi, tnu_gi, laplacian_gi
    real, dimension(u%dim, ele_loc(u, sele))          :: fnu_ele, tnu_ele
    real, dimension(u%dim)                            :: grad_fnu_dot_n, grad_tnu_dot_n
    real, dimension(face_ngi(u, sele))                :: f_scalar_gi, t_scalar_gi, les_coef_gi, beta_gi
    real, dimension(ele_loc(u, sele))                 :: laplacian
    real                                              :: strain_mod, t_strain_mod, denom
    character(len=OPTION_PATH_LEN)                    :: les_option_path, HOT


    ele = face_ele(x, sele)
    u_shape => face_shape(u, sele)
    x_shape => ele_shape(x, ele)
    fshape => face_shape(u, sele)
    source_shape => ele_shape(u, ele)
    lface = local_face_number(u, sele)
    ngi = face_ngi(u, sele)

    call transform_facet_to_physical(X, sele, detwei_f=detwei_bdy, normal=normal_bdy)

    ! Get shape fn and its gradient on surface element
    if(associated(source_shape%dn_s)) then
      augmented_shape = source_shape
      call incref(augmented_shape)
    else
      augmented_shape = make_element_shape(x_shape%loc, source_shape%dim, &
        & source_shape%degree, source_shape%quadrature, quad_s = fshape%quadrature)
    end if

    call compute_inverse_jacobian(x, ele, invj = invj)

    assert(x_shape%degree == 1)
    assert(ele_numbering_family(x_shape) == FAMILY_SIMPLEX)
    invj_face = spread(invj(:, :, 1), 3, size(invj_face, 3))
    !>GD --> This function is not used by the code any more
    !dshape_face = eval_volume_dshape_at_face_quad(augmented_shape, lface, invj_face)
    assert(associated(augmented_shape%dn_s))
    assert(size(invj_face, 1) == augmented_shape%dim)
    assert(size(invj_face, 2) == augmented_shape%dim)
    assert(size(invj_face, 3) == augmented_shape%surface_quadrature%ngi)
    assert(augmented_shape%dim == size(augmented_shape%dn_s, 4))
    assert(augmented_shape%loc == size(augmented_shape%dn_s, 1))
    assert(augmented_shape%surface_quadrature%ngi == size(augmented_shape%dn_s, 2))
    assert(lface <= size(augmented_shape%dn_s, 3))
    assert(augmented_shape%dim == size(augmented_shape%dn_s, 4))
    
    do iloc=1,augmented_shape%loc
        do gi=1,augmented_shape%surface_quadrature%ngi
        dshape_face(iloc,gi,:)=matmul(invj_face(:,:,gi),augmented_shape%dn_s(iloc,gi,lface,:))
        end do
    end do

    call deallocate(augmented_shape)

    ! first-filtered velocity at nodes and gauss pts
    fnu_ele = ele_val(fnu, ele)
    fnu_gi = face_val_at_quad(fnu, sele)

    ! test-filtered velocity at nodes and gauss pts
    tnu_ele = ele_val(tnu, ele)
    tnu_gi = face_val_at_quad(tnu, sele)

    ! Calculate velocity gradients w.r.t. surface normal
    do gi = 1, ngi

      ! Calculate grad fnu, tnu at the surface element quadrature points
      !grad_fnu = matmul(fnu_ele, dshape_face(:,gi,:))
      !grad_tnu = matmul(tnu_ele, dshape_face(:,gi,:))
      forall(i = 1:u%dim, j = 1:u%dim)
        grad_fnu(i,j) = dot_product(fnu_ele(i,:), dshape_face(:,gi,j))
        grad_tnu(i,j) = dot_product(tnu_ele(i,:), dshape_face(:,gi,j))
      end forall

      ! Calculate grad fnu dot dn at the surface element quadrature points
      do i=1, u%dim
        grad_fnu_dot_n(i) = dot_product(grad_fnu(i,:), normal_bdy(:,gi))
        grad_tnu_dot_n(i) = dot_product(grad_tnu(i,:), normal_bdy(:,gi))
      end do

    end do

    ! SGS model quantities common to both models

    ! Filtered SGS stress tensor (from dynamic LES computation)
    !tau_gi = face_val_at_quad(sgstensor, sele)
    ! Leonard tensor at Gauss points
    leonard_gi = face_val_at_quad(leonard, sele)

    ! scalar_gi first filter width G1 = alpha^2*meshsize (units length^2)
    f_scalar_gi = alpha**2/24*length_scale_scalar(x, ele)
    ! Combined width G2 = (1+gamma^2)*G1
    t_scalar_gi = (1.0+gamma**2)*f_scalar_gi

    les_option_path = trim(u%option_path)//"/prognostic/spatial_discretisation/continuous_galerkin/les_model"

    ! Choose SGS model
    if(dynamic_les) then

      ! Get some dynamic LES model terms at the surface quadrature points

      ! Get strain terms S1, S2
      strain_gi = les_strain_rate(dshape_face, fnu_ele)
      t_strain_gi = les_strain_rate(dshape_face, tnu_ele)

      ! Dynamic Smagorinsky coeff and strain product at Gauss points
      strainprod_gi = face_val_at_quad(strainprod, sele)
      les_coef_gi = face_val_at_quad(dyn_coeff, sele)

      do gi=1, ngi
        !ewrite(2,*) "gi: ", gi

        ! Recompute dynamic LES model terms

        ! Strain moduli |S1|, |S2|
        strain_mod = sqrt(2*sum(strain_gi(:,:,gi)*strain_gi(:,:,gi)))
        t_strain_mod = sqrt(2*sum(t_strain_gi(:,:,gi)*t_strain_gi(:,:,gi)))

        ! Tensor M_ij = (|S2|*S2)G2 - ((|S1|S1)^f2)G1
        mij = t_strain_mod*t_strain_gi(:,:,gi)*t_scalar_gi(gi) - strainprod_gi(:,:,gi)*f_scalar_gi(gi)
        denom = sum(mij*mij)

        ! Dynamic LES Model coeff C_S = -(L_ij M_ij) / 2(M_ij M_ij)
        ! Subtract tnu_i tnu_j from Leonard tensor
        if (denom > epsilon(0.0)) then
          les_coef_gi(gi) = -0.5*sum((leonard_gi(:,:,gi)-outer_product(tnu_gi(:,gi), tnu_gi(:,gi)))*mij) / denom
        else
          les_coef_gi(gi) = 0.0
        endif
        !ewrite(2,*) "Cs: ", les_coef_gi(gi)

        ! Constrain C_S to be between 0 and 0.04.
        les_coef_gi(gi) = min(max(les_coef_gi(gi),0.0), 0.04)

        ! Dynamic LES Model terms: T_ij = -2C_S|S2|.alpha^2.G2.S2, tau_ij = -2C_S|S1|.alpha^2.G1.S1
        tau_gi(:,:,gi)     = 2.0*les_coef_gi(gi)*strain_mod*f_scalar_gi(gi)*strain_gi(:,:,gi)
        tautest_gi(:,:,gi) = 2.0*les_coef_gi(gi)*t_strain_mod*t_scalar_gi(gi)*t_strain_gi(:,:,gi)

      end do

    ! Germano's exact SGS closure or the LANS-alpha model
    else if(exact_sgs) then

      ewrite(2,*) "calculate_boundaries: ", trim(les_option_path)//"/exact_sgs/calculate_boundaries"

      ! If setting SGS tensor = 0 on boundaries, do not calculate tau_gi, tautest_gi here.
      if(have_option(trim(les_option_path)//"/exact_sgs/calculate_boundaries")) then

        do gi = 1, ngi
          tau_gi(:,:,gi) = 0.0
          tautest_gi(:,:,gi) = 0.0
        end do

      ! Compute SGS tensors at both filter levels
      else

        select case(HOT)

        ! Clark or Gradient model
        case("none")
          do gi = 1, ngi
            ! SGS tensors at both filter levels composed of 1 gradient product term:
            ! (du_i/dx_k) (du_j/dx_k)
            tau_gi(:,:,gi)     = f_scalar_gi(gi)*(matmul(grad_fnu, transpose(grad_fnu)))
            tautest_gi(:,:,gi) = t_scalar_gi(gi)*(matmul(grad_tnu, transpose(grad_tnu)))
          end do

        ! Germano's exact closure
        case("exact")
          do gi = 1, ngi
            ! SGS tensors at both filter levels composed of 1 gradient product term:
            ! 2 (du_i/dx_k) (du_j/dx_k)
            tau_gi(:,:,gi) = 2.0*f_scalar_gi(gi)*matmul(grad_fnu, transpose(grad_fnu))
            tautest_gi(:,:,gi) = 2.0*t_scalar_gi(gi)*matmul(grad_tnu, transpose(grad_tnu))

            ! Laplacian operator: dim.dim = dim
            do i = 1, ele_loc(fnu, ele)
              laplacian(i) = dot_product(dshape_face(i,gi,:), dshape_face(i,gi,:))
            end do

            ! Laplacian of filtered velocity: (dim*loc)*loc = dim
            laplacian_gi(:,gi) = matmul(fnu_ele, laplacian)

            ! Add 4th-order Laplacian product term to SGS tensor
            tau_gi(:,:,gi)     = tau_gi(:,:,gi) &
                                 + f_scalar_gi(gi)**4*outer_product(laplacian_gi(:,gi), laplacian_gi(:,gi))

            ! Laplacian of test-filtered velocity
            laplacian_gi(:,gi) = matmul(tnu_ele, laplacian)

            ! SGS tensor at test filter level
            tautest_gi(:,:,gi) = tautest_gi(:,:,gi) &
                                 + t_scalar_gi(gi)**4*outer_product(laplacian_gi(:,gi), laplacian_gi(:,gi))
          end do

        ! LANS-alpha closure
        case("lans")
          do gi = 1, ngi
            ! SGS tensors at both filter levels composed of 3 gradient product terms:
            ! (du_i/dx_k) (du_j/dx_k) - (du_k/dx_i) (du_k/dx_j) + (du_i/dx_k) (du_k/dx_j)
            tau_gi(:,:,gi)     = f_scalar_gi(gi)*(matmul(grad_fnu, transpose(grad_fnu)) &
                                                - matmul(transpose(grad_fnu), grad_fnu) &
                                                + matmul(grad_fnu, grad_fnu))

            tautest_gi(:,:,gi) = t_scalar_gi(gi)*(matmul(grad_tnu, transpose(grad_tnu)) &
                                                - matmul(transpose(grad_tnu), grad_tnu) &
                                                + matmul(grad_tnu, grad_tnu))
          end do

        end select
      end if
    end if

    ! Now compute dynamic slip wall model terms
    do gi=1, ngi

      ! Recalculate M_ij
      forall(i = 1:u%dim, j = 1:u%dim)
        mij(i,j) = gamma**2*grad_tnu_dot_n(i)*grad_tnu_dot_n(j) - grad_fnu_dot_n(i)*grad_fnu_dot_n(j)
      end forall

      ! Alternative form. UNSTABLE - GIVES BETA < 0
      !forall(i = 1:u%dim, j = 1:u%dim)
      !  mij(i,j) = gamma*(tnu_gi(i,gi)*grad_tnu_dot_n(j) + tnu_gi(j,gi)*grad_tnu_dot_n(i)) &
      !           & - fnu_gi(i,gi)*grad_fnu_dot_n(j) - fnu_gi(j,gi)*grad_fnu_dot_n(i)
      !end forall

      denom = sum(mij*mij)

      ! Subtract fnu_i fnu_j from Leonard tensor
      leonard_gi(:,:,gi) = leonard_gi(:,:,gi) - outer_product(fnu_gi(:,gi), fnu_gi(:,gi))

      ! Slip coeff beta = sqrt((L_ij-T_ij+tau_ij) M_ij / M_ij M_ij)
      beta_gi(gi) = sum(mij(:,:)*(leonard_gi(:,:,gi) - tautest_gi(:,:,gi) + tau_gi(:,:,gi))) / denom
      !beta_gi(gi) = sum(mij(:,:)*leonard_gi(:,:,gi)) / denom

      if (beta_gi(gi) > epsilon(0.0)) then
        beta_gi(gi) = 1.0/sqrt(beta_gi(gi))
      ! If using alternative form, don't take sqrt of beta.
      !if (abs(beta_gi(gi)) > epsilon(0.0)) then
      !  beta_gi(gi) = 1.0/beta_gi(gi)
      else ! what is a sensible alternative? Face-averaged value?
        beta_gi(gi) = 0.0
      endif

      !ewrite(2,*) "Cs: ", les_coef_gi(gi)
      ewrite(2,*) "Lij: ", leonard_gi(:,:,gi)
      ewrite(2,*) "tau: ", tau_gi(:,:,gi)
      ewrite(2,*) "tautest: ", tautest_gi(:,:,gi)
      !ewrite(2,*) "d/dxj(fnu): ", grad_fnu
      !ewrite(2,*) "d/dxj(tnu): ", grad_tnu
      !ewrite(2,*) "d/dn(fnu): ", grad_fnu_dot_n
      !ewrite(2,*) "d/dn(tnu): ", grad_tnu_dot_n
      ewrite(2,*) "Mij: ", mij
      !ewrite(2,*) "MijMij: ", denom

    end do

    ewrite(2,*) "beta: ", beta_gi

    ! add contribution to surface field
    surface_field_addto = shape_rhs(u_shape, beta_gi*detwei_bdy)

  end subroutine calculate_dynamic_slip_coefficient_face

end module les_module
