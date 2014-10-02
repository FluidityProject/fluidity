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
  use state_module
  use fields
  use field_options
  use spud
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use smoothing_module
  use boundary_conditions
  use vector_tools
  use fetools
  use state_fields_module
  use solvers
  implicit none

  private

  public les_viscosity_strength, wale_viscosity_strength
  public les_init_diagnostic_fields, les_assemble_diagnostic_fields, les_solve_diagnostic_fields
  public compute_les_local_fields, les_strain_rate, exact_sgs_stress, calculate_periodic_channel_forcing
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
  
  subroutine compute_les_local_fields(nu, positions, fnu, tnu, leonard, strainprod, sgstensor, alpha, gamma, length_scale_type, path, dynamic_les, exact_sgs)

    type(vector_field), intent(in)                        :: positions
    type(vector_field), intent(inout)                     :: nu, fnu, tnu
    ! Leonard tensor, strain product, SGS tensor
    type(tensor_field), intent(inout)                     :: leonard, strainprod, sgstensor
    ! Filter scale factors
    real                                      :: alpha, gamma
    character(len=OPTION_PATH_LEN)            :: length_scale_type
    character(len=OPTION_PATH_LEN), intent(in)            :: path
    logical, intent(in)                                   :: dynamic_les, exact_sgs

    ! Local quantities
    type(tensor_field), pointer                           :: tfield
    character(len=OPTION_PATH_LEN)                        :: lpath
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
    logical                                               :: lans
    type(element_type)                                    :: shape_nu

    if(dynamic_les) then

      ! Path is to level above solver options
      lpath = (trim(path)//"/dynamic_les")

    else if(exact_sgs) then

      lpath = (trim(path)//"/exact_sgs")
      ! Choose Germano's exact SGS closure or the LANS-alpha model
      lans = have_option(trim(lpath)//"/lans")
      ! have to initialise these values here as they're not used in the exact SGS tensor
      gamma = 2.0
      length_scale_type = "scalar"

    end if

    ewrite(2,*) "filter factor alpha: ", alpha
    ewrite(2,*) "filter factor gamma: ", gamma
    ewrite(2,*) "filter width type: ", trim(length_scale_type)

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

      ! Filter SGS tensor field
      call zero(tfield)

      if(length_scale_type=="scalar") then
        call smooth_tensor(sgstensor, positions, tfield, alpha, lpath)
      else if(length_scale_type=="tensor") then
        call smooth_tensor(sgstensor, positions, tfield, alpha, lpath)
      end if

      call set(sgstensor, tfield)

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

        ! LANS-alpha closure
        if(lans) then
          do gi = 1, ele_ngi(nu, ele)
            ! velocity gradient
            mat = matmul(ele_val(fnu, ele), du_t(:,gi,:))
            ! 3 gradient product terms: (du_i/dx_k) (du_j/dx_k) - (du_k/dx_i) (du_k/dx_j) + (du_i/dx_k) (du_k/dx_j)
            mat_gi(:,:,gi) = matmul(mat, transpose(mat)) - matmul(transpose(mat), mat) + matmul(mat, mat)

            mat_gi(:,:,gi) = delta*mat_gi(:,:,gi)
          end do
        ! Germano's closure
        else
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
        end if

        ! add to tensor RHS field
        call addto(tfield, ele_nodes(nu, ele), shape_tensor_rhs(shape_nu, mat_gi, detwei))

      end do

      ! Filter tensor RHS field
      if(length_scale_type=="scalar") then
        call smooth_tensor(tfield, positions, sgstensor, alpha, lpath)
      else if(length_scale_type=="tensor") then
        call smooth_tensor(tfield, positions, sgstensor, alpha, lpath)
      end if

    end if

    ! Deallocates
    deallocate(u_loc, t_loc)
    call deallocate(tfield)
    !call deallocate(tui_tuj)
    deallocate(tfield)
    !deallocate(tui_tuj)

  end subroutine compute_les_local_fields

  subroutine exact_sgs_stress(nu, positions, fnu, tnu, leonard, exactsgs, alpha, length_scale_type, path)

    type(vector_field), intent(in)                        :: positions
    type(tensor_field), pointer                           :: exactsgs, leonard
    ! Unfiltered velocity and nonlinear velocity
    type(vector_field), pointer                           :: nu, fnu, tnu
    type(vector_field), pointer                           :: vfield
    ! RHS tensor
    type(tensor_field), pointer                           :: tensorrhs
    ! Scale factor
    real, intent(in)                                      :: alpha
    character(len=OPTION_PATH_LEN), intent(in)            :: length_scale_type, path
    ! Local quantities
    character(len=OPTION_PATH_LEN)                        :: lpath
    integer                                               :: ele, gi, loc
    real, dimension(ele_loc(nu,1), ele_ngi(nu,1), nu%dim) :: du_t
    real, dimension(ele_ngi(nu,1))                        :: detwei
    real, dimension(nu%dim, nu%dim, ele_ngi(nu,1))        :: mat_gi
    real, dimension(nu%dim, nu%dim)                       :: mat
    real, dimension(nu%dim, ele_ngi(nu,1))                :: laplacian_gi
    real, dimension(ele_loc(nu,1))                        :: laplacian
    real                                                  :: delta
    type(element_type)                                    :: shape_nu
    logical                                               :: lans

    ! Choose Germano's exact SGS closure or the LANS-alpha model
    lans = have_option(trim(path)//"/exact_sgs/lans")

    ! Path is to level above solver options
    lpath = (trim(path)//"/exact_sgs")
    ewrite(2,*) "filter factor alpha: ", alpha
    ewrite(2,*) "filter width type: ", trim(length_scale_type)

    ! Compute explicitly filtered velocity - stick with scalar filter width for simplicity
    if(length_scale_type=="scalar") then
      call smooth_vector(nu, positions, fnu, alpha, lpath)
    else if(length_scale_type=="tensor") then
      call smooth_vector(nu, positions, fnu, alpha, lpath)
      !call anisotropic_smooth_vector(nu, positions, fnu, alpha, lpath)
    end if

    ewrite_minmax(nu)
    ewrite_minmax(fnu)

    allocate(vfield); allocate(tensorrhs)
    call allocate(vfield, nu%dim, nu%mesh, "Laplacian")
    call allocate(tensorrhs, nu%mesh, "TensorRHS")
    call zero(vfield); call zero(tensorrhs)

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

      ! LANS-alpha closure
      if(lans) then
        do gi = 1, ele_ngi(nu, ele)
          ! velocity gradient
          mat = matmul(ele_val(fnu, ele), du_t(:,gi,:))
          ! 3 gradient product terms: (du_i/dx_k) (du_j/dx_k) - (du_k/dx_i) (du_k/dx_j) + (du_i/dx_k) (du_k/dx_j)
          mat_gi(:,:,gi) = matmul(mat, transpose(mat)) - matmul(transpose(mat), mat) + matmul(mat, mat)

          mat_gi(:,:,gi) = delta*mat_gi(:,:,gi)
        end do
      ! Germano's closure
      else
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
      end if

      ! add to tensor RHS field
      call addto(tensorrhs, ele_nodes(nu, ele), shape_tensor_rhs(shape_nu, mat_gi, detwei))

    end do

    ! Filter tensor RHS field
    if(length_scale_type=="scalar") then
      call smooth_tensor(tensorrhs, positions, exactsgs, alpha, lpath)
    else if(length_scale_type=="tensor") then
      call smooth_tensor(tensorrhs, positions, exactsgs, alpha, lpath)
      !call anisotropic_smooth_tensor(tensorrhs, positions, exactsgs, alpha, lpath)
    end if

    ! Deallocates
    call deallocate(vfield)
    call deallocate(tensorrhs)
    deallocate(vfield); deallocate(tensorrhs)

  end subroutine exact_sgs_stress

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


end module les_module
