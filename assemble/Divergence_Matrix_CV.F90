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

module divergence_matrix_cv
  use quadrature
  use fields
  use field_derivatives
  use state_module
  use futils
  use fetools
  use spud
  use cv_faces
  use cv_shape_functions
  use cvtools
  use cv_fields
  use cv_upwind_values
  use cv_face_values
  use cv_options
  use diagnostic_fields, only: calculate_diagnostic_variable
  use boundary_conditions
  use global_parameters, only: OPTION_PATH_LEN
  use field_options, only: get_coordinate_field
  use sparsity_patterns_meshes
  implicit none

  private
  public :: assemble_divergence_matrix_cv, assemble_compressible_divergence_matrix_cv

contains

    !************************************************************************
    subroutine assemble_divergence_matrix_cv(CT_m, state, ct_rhs, & 
                                             test_mesh, field, &
                                             get_ct, exclude_boundaries)

      ! inputs/outputs
      ! bucket full of fields
      type(state_type), intent(inout) :: state

      ! the pressure gradient and compressible gradient matrices
      type(block_csr_matrix), intent(inout) :: CT_m

      type(scalar_field), intent(inout), optional :: ct_rhs

      type(mesh_type), intent(in) :: test_mesh

      type(vector_field), intent(inout) :: field
      
      logical, intent(in), optional :: get_ct
      logical, intent(in), optional :: exclude_boundaries

      ! local
      ! degree of quadrature over cv faces
      integer :: quaddegree

      real, dimension(:,:), allocatable :: x_ele, x_ele_bdy
      real, dimension(:,:), allocatable :: x_f
      real, dimension(:,:), allocatable :: normal, normal_bdy
      real, dimension(:), allocatable :: detwei, detwei_bdy
      real, dimension(:), allocatable :: normgi
      integer, dimension(:), pointer :: test_nodes, field_nodes, x_test_nodes
      integer, dimension(:), allocatable :: test_nodes_bdy, field_nodes_bdy

      logical, dimension(:), allocatable :: notvisited

      ! loop integers
      integer :: ele, sele, iloc, oloc, jloc, face, gi, ggi, dim

      ! information about cv faces
      type(cv_faces_type) :: cvfaces
      ! shape functions for region and surface
      type(element_type) :: x_cvshape, x_cvbdyshape
      type(element_type) :: test_cvshape, test_cvbdyshape
      type(element_type) :: field_cvshape, field_cvbdyshape
      ! pointer to coordinates
      type(vector_field), pointer :: x
      type(vector_field) :: x_test

      integer, dimension(:,:), allocatable :: field_bc_type
      real, dimension(:,:), allocatable :: field_bc_val
      type(vector_field) :: field_bc

      real, dimension(:,:,:), allocatable :: ct_mat_local, ct_mat_local_bdy
      real, dimension(:), allocatable :: ct_rhs_local

      logical :: l_get_ct

      ! =============================================================
      ! Subroutine to construct the matrix CT_m (a.k.a. C1/2/3T).
      ! =============================================================

      ewrite(1,*) 'In assemble_divergence_matrix_cv'

      if(present(get_ct)) then
        l_get_ct = get_ct
      else
        l_get_ct = .true.
      end if

      x=>extract_vector_field(state, "Coordinate")
      allocate(x_ele(x%dim,x%mesh%shape%loc))

      call get_option("/geometry/quadrature/controlvolume_surface_degree", &
                     quaddegree, default=1)

      ! Clear memory of arrays being designed
      if(present(ct_rhs)) call zero(ct_rhs)

      cvfaces=find_cv_faces(vertices=ele_vertices(test_mesh, 1), &
                            dimension=mesh_dim(test_mesh), &
                            polydegree=test_mesh%shape%degree, &
                            quaddegree=quaddegree)

      if(l_get_ct) then
      
        call zero(CT_m)
      
        x_test=get_coordinate_field(state, test_mesh)
        
        x_cvshape=make_cv_element_shape(cvfaces, x%mesh%shape%degree)
        test_cvshape=make_cv_element_shape(cvfaces, test_mesh%shape%degree)
        field_cvshape=make_cv_element_shape(cvfaces, field%mesh%shape%degree)
  
        allocate(x_f(x%dim, x_cvshape%ngi), &
                detwei(x_cvshape%ngi), &
                normal(x%dim, x_cvshape%ngi), &
                normgi(x%dim), &
                ct_mat_local(x%dim, test_mesh%shape%loc, field%mesh%shape%loc))
  
        allocate(notvisited(x_cvshape%ngi))

        element_loop: do ele=1, element_count(test_mesh)
          x_ele=ele_val(x, ele)
          x_f=ele_val_at_quad(x, ele, x_cvshape)
          test_nodes=>ele_nodes(test_mesh, ele)
          field_nodes=>ele_nodes(field, ele)
          x_test_nodes=>ele_nodes(x_test, ele)
  
          call transform_cvsurf_to_physical(x_ele, x_cvshape, &
                                            detwei, normal, cvfaces)
          notvisited=.true.
  
          ct_mat_local = 0.0
  
          nodal_loop_i: do iloc = 1, test_mesh%shape%loc
  
            face_loop: do face = 1, cvfaces%faces
  
              if(cvfaces%neiloc(iloc, face) /= 0) then
                oloc = cvfaces%neiloc(iloc, face)
  
                quadrature_loop: do gi = 1, cvfaces%shape%ngi
  
                  ggi = (face-1)*cvfaces%shape%ngi + gi
  
                  if(notvisited(ggi)) then
                    notvisited(ggi)=.false.
  
                    normgi=orientate_cvsurf_normgi(node_val(x_test, x_test_nodes(iloc)),x_f(:,ggi),normal(:,ggi))
  
                    nodal_loop_j: do jloc = 1, field%mesh%shape%loc
  
                      inner_dimension_loop: do dim = 1, size(normgi)
  
                        ct_mat_local(dim, iloc, jloc) = ct_mat_local(dim, iloc, jloc) &
                                                      + field_cvshape%n(jloc, ggi)*detwei(ggi)*normgi(dim)
                        ct_mat_local(dim, oloc, jloc) = ct_mat_local(dim, oloc, jloc) &
                                                      + field_cvshape%n(jloc, ggi)*detwei(ggi)*(-normgi(dim)) ! notvisited
  
                      end do inner_dimension_loop
  
                    end do nodal_loop_j
  
                  end if ! notvisited
  
                end do quadrature_loop
  
              end if
  
            end do face_loop
  
          end do nodal_loop_i
  
          outer_dimension_loop: do dim = 1, size(normgi)
  
            call addto(CT_m, 1, dim, test_nodes, field_nodes, ct_mat_local(dim,:,:))
  
          end do outer_dimension_loop
  
        end do element_loop
      
        call deallocate(x_cvshape)
        call deallocate(test_cvshape)
        call deallocate(field_cvshape)
        deallocate(x_f, detwei, normal, normgi)
        deallocate(notvisited)
        call deallocate(x_test)
      end if

      if(.not.present_and_true(exclude_boundaries)) then

        x_cvbdyshape=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape%degree)
        test_cvbdyshape=make_cvbdy_element_shape(cvfaces, test_mesh%faces%shape%degree)
        field_cvbdyshape=make_cvbdy_element_shape(cvfaces, field%mesh%faces%shape%degree)
  
        assert(surface_element_count(test_mesh)==surface_element_count(field))
        allocate(field_bc_type(field%dim, surface_element_count(test_mesh)))
        call get_entire_boundary_condition(field, (/"weakdirichlet ", &
                                                    "no_normal_flow", &
                                                    "periodic      ", &
                                                    "free_surface  "/), field_bc, field_bc_type)
  
        allocate(x_ele_bdy(x%dim,x%mesh%faces%shape%loc), &
                detwei_bdy(x_cvbdyshape%ngi), &
                normal_bdy(x%dim, x_cvbdyshape%ngi), &
                field_bc_val(field_bc%dim, field_bc%mesh%shape%loc))
        allocate(field_nodes_bdy(field%mesh%faces%shape%loc))
        allocate(test_nodes_bdy(test_mesh%faces%shape%loc))
        allocate(ct_mat_local_bdy(x%dim, test_mesh%faces%shape%loc, field%mesh%faces%shape%loc), &
                ct_rhs_local(test_mesh%faces%shape%loc))
  
        surface_element_loop: do sele = 1, surface_element_count(test_mesh)
  
          ! cycle if this is a no_normal_flow or a periodic or a free_surface boundary then cycle
          if(any(field_bc_type(:,sele)==2).or.any(field_bc_type(:,sele)==3).or.any(field_bc_type(:,sele)==4)) cycle
          
          ! cycle if there's no rhs present or there's no weakdirichlet conditions or we're not
          ! assembling the matrix
          if(.not.(present(ct_rhs).and.any(field_bc_type(:,sele)==1)).and..not.l_get_ct) cycle
  
          ele = face_ele(x, sele)
          x_ele = ele_val(x, ele)
          x_ele_bdy = face_val(x, sele)
          test_nodes_bdy=face_global_nodes(test_mesh, sele)
          field_nodes_bdy=face_global_nodes(field, sele)
  
          if(any(field_bc_type(:, sele)==1)) then
            field_bc_val = ele_val(field_bc, sele)
          else
            field_bc_val = 0.0
          end if
  
          call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, &
                                x_cvbdyshape, normal_bdy, detwei_bdy)
  
          ct_mat_local_bdy = 0.0
          ct_rhs_local = 0.0
  
          surface_nodal_loop_i: do iloc = 1, test_mesh%faces%shape%loc
  
            surface_face_loop: do face = 1, cvfaces%sfaces
              if(cvfaces%sneiloc(iloc,face)/=0) then
  
                surface_quadrature_loop: do gi = 1, cvfaces%shape%ngi
  
                  ggi = (face-1)*cvfaces%shape%ngi + gi
  
                  surface_nodal_loop_j: do jloc = 1, field%mesh%faces%shape%loc
  
                    surface_inner_dimension_loop: do dim = 1, size(normal_bdy,1)
  
                      if((present(ct_rhs)).and.(field_bc_type(dim, sele)==1)) then
  
                        ct_rhs_local(iloc) = ct_rhs_local(iloc) + &
                              field_cvbdyshape%n(jloc,ggi)*detwei_bdy(ggi)*normal_bdy(dim,ggi)*field_bc_val(dim, jloc)
  
                      else
  
                        ct_mat_local_bdy(dim, iloc, jloc) =  ct_mat_local_bdy(dim, iloc, jloc) + &
                              field_cvbdyshape%n(jloc,ggi)*detwei_bdy(ggi)*normal_bdy(dim, ggi)
  
                      end if
  
                    end do surface_inner_dimension_loop
  
                  end do surface_nodal_loop_j
  
                end do surface_quadrature_loop
  
              end if
  
            end do surface_face_loop
  
          end do surface_nodal_loop_i
  
          surface_outer_dimension_loop: do dim = 1, size(normal_bdy,1)
  
            if((present(ct_rhs)).and.(field_bc_type(dim, sele)==1)) then
  
              call addto(ct_rhs, test_nodes_bdy, ct_rhs_local)
  
            elseif(l_get_ct) then
  
              call addto(CT_m, 1, dim, test_nodes_bdy, field_nodes_bdy, ct_mat_local_bdy(dim,:,:))
  
            end if
  
          end do surface_outer_dimension_loop
  
        end do surface_element_loop
  
        call deallocate(field_bc)
        deallocate(field_bc_type)
        call deallocate(x_cvbdyshape)
        call deallocate(test_cvbdyshape)
        call deallocate(field_cvbdyshape)
        deallocate(x_ele_bdy, detwei_bdy, normal_bdy)
        deallocate(test_nodes_bdy, field_nodes_bdy)
        
      end if

      call deallocate(cvfaces)
      deallocate(x_ele)
      
    end subroutine assemble_divergence_matrix_cv
    !************************************************************************

    subroutine assemble_compressible_divergence_matrix_cv(CTP_m, state, ct_rhs)
      
      ! inputs/outputs
      ! bucket full of fields
      type(state_type), dimension(:), intent(inout) :: state

      ! the compressible gradient matrices
      type(block_csr_matrix), intent(inout) :: CTP_m

      type(scalar_field), intent(inout), optional :: ct_rhs

      if((size(state)==1).and.(.not.has_scalar_field(state(1), "MaterialVolumeFraction"))) then
      
        call assemble_1mat_compressible_divergence_matrix_cv(CTP_m, state(1), ct_rhs)
        
      else
      
        call assemble_mmat_compressible_divergence_matrix_cv(CTP_m, state)
        
      end if
    
    end subroutine assemble_compressible_divergence_matrix_cv

    subroutine assemble_1mat_compressible_divergence_matrix_cv(CTP_m, state, ct_rhs)

      ! inputs/outputs
      ! bucket full of fields
      type(state_type), intent(inout) :: state

      ! the pressure gradient and compressible gradient matrices
      type(block_csr_matrix), intent(inout) :: CTP_m
      
      type(scalar_field), intent(inout), optional :: ct_rhs

      ! local
      ! degree of quadrature over cv faces
      integer :: quaddegree

      real, dimension(:,:), allocatable :: x_ele, x_ele_bdy
      real, dimension(:,:), allocatable :: x_f, u_f, u_bdy_f
      real, dimension(:,:), allocatable :: normal, normal_bdy
      real, dimension(:), allocatable :: detwei, detwei_bdy
      real, dimension(:), allocatable :: normgi
      integer, dimension(:), pointer :: u_nodes, p_nodes, x_nodes
      integer, dimension(:), allocatable :: u_nodes_bdy, p_nodes_bdy
      real :: dens_theta_val
      real :: dens_face_val
      real :: olddens_face_val
      real, dimension(:), allocatable :: dens_ele, olddens_ele, &
                                         norm_ele
      real, dimension(:), allocatable :: dens_ele_bdy, olddens_ele_bdy, &
                                         norm_ele_bdy
      real, dimension(:), allocatable :: ghost_dens_ele_bdy, ghost_olddens_ele_bdy

      logical, dimension(:), allocatable :: notvisited

      ! loop integers
      integer :: ele, sele, iloc, oloc, jloc, face, gi, ggi, dim

      ! information about cv faces
      type(cv_faces_type) :: cvfaces
      ! shape functions for region and surface
      type(element_type) :: x_cvshape, x_cvbdyshape
      type(element_type) :: u_cvshape, u_cvbdyshape
      type(element_type) :: dens_cvshape, dens_cvbdyshape
      ! pointer to velocity
      type(vector_field), pointer :: u, ug, x
      type(vector_field) :: relu, x_p
      type(scalar_field), pointer :: p
      type(scalar_field) :: cfl_no

      type(csr_sparsity), pointer :: mesh_sparsity
      type(csr_matrix)  :: dens_upwind, olddens_upwind

      type(scalar_field), pointer :: dens, olddens

      real, dimension(:), allocatable :: cfl_ele

      real :: face_value

      integer, dimension(:), allocatable :: dens_bc_type
      type(scalar_field) :: dens_bc

      real :: udotn, income
      logical :: inflow
      real :: dt

      type(cv_options_type) :: dens_options

      type(scalar_field), pointer :: normalisation
      character(len=FIELD_NAME_LEN) :: normalisation_field
      integer :: norm_stat

      integer, dimension(:,:), allocatable :: velocity_bc_type
      real, dimension(:,:), allocatable :: velocity_bc_val
      type(vector_field) :: velocity_bc

      real, dimension(:,:,:), allocatable :: ctp_mat_local, ctp_mat_local_bdy
      real, dimension(:), allocatable :: ct_rhs_local
      
      ! temporary hack to get around compiler failure to construct arrays of characters
      character(len=OPTION_PATH_LEN), dimension(1) :: option_path_array

      ! =============================================================
      ! Subroutine to construct the matrix CTP_m (a.k.a. C1/2/3TP).
      ! =============================================================

      ewrite(2,*) 'In assemble_1mat_compressible_divergence_matrix_cv'

      x=>extract_vector_field(state, "Coordinate")
      
      call get_option("/timestepping/timestep", dt)
      call get_option("/geometry/quadrature/controlvolume_surface_degree", &
                     quaddegree, default=1)
      
      dens=>extract_scalar_field(state, "Density")
      olddens=>extract_scalar_field(state, "OldDensity")
      
      u=>extract_vector_field(state, "NonlinearVelocity") ! maybe this should be updated after the velocity solve?
      ug=>extract_vector_field(state, "GridVelocity")
      
      call allocate(relu, u%dim, u%mesh, "RelativeVelocity")
      call set(relu, u)
      call addto(relu, ug, -1.0)
      
      p=>extract_scalar_field(state, "Pressure")

      x_p = get_coordinate_field(state, p%mesh)
      

      cvfaces=find_cv_faces(vertices=ele_vertices(dens%mesh, 1), &
                            dimension=mesh_dim(dens%mesh), &
                            polydegree=dens%mesh%shape%degree, &
                            quaddegree=quaddegree)

      x_cvshape=make_cv_element_shape(cvfaces, x%mesh%shape%degree)
      dens_cvshape=make_cv_element_shape(cvfaces, dens%mesh%shape%degree)
      u_cvshape=make_cv_element_shape(cvfaces, u%mesh%shape%degree)

      mesh_sparsity=>get_csr_sparsity_firstorder(state, dens%mesh, dens%mesh)
      call allocate(dens_upwind, mesh_sparsity, name="DensityUpwindValues")
      call allocate(olddens_upwind, mesh_sparsity, name="OldDensityUpwindValues")

      ! get all the relevent options for density
      ! handily wrapped in a new type...
      dens_options = get_cv_options(dens%option_path, dens%mesh%shape%numbering%family)

      if(need_upwind_values(dens_options)) then

        call find_upwind_values(state, x_p, dens, dens_upwind, &
                                olddens, olddens_upwind)

      else

        call zero(dens_upwind)
        call zero(olddens_upwind)

      end if

      call get_option(trim(p%option_path)//"/prognostic/scheme/use_compressible_projection_method/normalisation/name", &
                      normalisation_field, stat=norm_stat)

      ! get the normalisation field (if we need one)
      if(norm_stat==0) then
        normalisation=>extract_scalar_field(state, trim(normalisation_field))
      else
        allocate(normalisation)
        call allocate(normalisation, p%mesh, name="DummyNormalisation", field_type=FIELD_TYPE_CONSTANT)
        call set(normalisation, 1.0)
      end if

      ! find courant number (if needed)
      option_path_array(1) = trim(dens%option_path)  ! temporary hack for compiler failure
      call cv_disc_get_cfl_no(option_path_array, &
                      state, dens%mesh, cfl_no)

      ! Clear memory of arrays being designed
      call zero(CTP_m)
      if(present(ct_rhs)) call zero(ct_rhs)

      allocate(x_ele(x%dim,ele_loc(x,1)), &
               x_f(x%dim, x_cvshape%ngi), &
               u_f(u%dim, u_cvshape%ngi), &
               detwei(x_cvshape%ngi), &
               normal(x%dim, x_cvshape%ngi), &
               normgi(x%dim))
      allocate(cfl_ele(ele_loc(p,1)), &
               dens_ele(ele_loc(p,1)), &
               olddens_ele(ele_loc(p,1)), &
               norm_ele(ele_loc(normalisation,1)))
      allocate(notvisited(x_cvshape%ngi))
      allocate(ctp_mat_local(x%dim, p%mesh%shape%loc, u_cvshape%loc))

      element_loop: do ele=1, element_count(p)
        x_ele=ele_val(x, ele)
        x_f=ele_val_at_quad(x, ele, x_cvshape)
        u_f=ele_val_at_quad(relu, ele, u_cvshape)
        p_nodes=>ele_nodes(p, ele)
        u_nodes=>ele_nodes(u, ele)
        x_nodes=>ele_nodes(x_p, ele)

        call transform_cvsurf_to_physical(x_ele, x_cvshape, &
                                          detwei, normal, cvfaces)
        cfl_ele = ele_val(cfl_no, ele)

        dens_ele = ele_val(dens, ele)
        olddens_ele = ele_val(olddens, ele)

        norm_ele = ele_val(normalisation, ele)

        notvisited=.true.

        ctp_mat_local = 0.0

        nodal_loop_i: do iloc = 1, p%mesh%shape%loc

          face_loop: do face = 1, cvfaces%faces

            if(cvfaces%neiloc(iloc, face) /= 0) then
              oloc = cvfaces%neiloc(iloc, face)

              quadrature_loop: do gi = 1, cvfaces%shape%ngi

                ggi = (face-1)*cvfaces%shape%ngi + gi

                if(notvisited(ggi)) then
                  notvisited(ggi)=.false.

                  normgi=orientate_cvsurf_normgi(node_val(x_p, x_nodes(iloc)),x_f(:,ggi),normal(:,ggi))

                  udotn=dot_product(u_f(:,ggi), normgi(:))

                  inflow = (udotn<=0.0)

                  income = merge(1.0,0.0,inflow)

                  call evaluate_face_val(dens_face_val, olddens_face_val, &
                                         iloc, oloc, ggi, x_nodes, &
                                         dens_cvshape,&
                                         dens_ele, olddens_ele, &
                                         dens_upwind, olddens_upwind, &
                                         inflow, cfl_ele, &
                                         dens_options)

                  dens_theta_val=theta_val(iloc, oloc, &
                                           dens_face_val, &
                                           olddens_face_val, &
                                           dens_options%theta, dt, udotn, &
                                           x_ele, dens_options%limit_theta, &
                                           dens_ele, olddens_ele)


                  nodal_loop_j: do jloc = 1, u_cvshape%loc

                    face_value = u_cvshape%n(jloc, ggi)*detwei(ggi)*dens_theta_val

                    inner_dimension_loop: do dim = 1, size(normgi)

                      ctp_mat_local(dim, iloc, jloc) = ctp_mat_local(dim, iloc, jloc) &
                                                      + face_value*normgi(dim)/norm_ele(iloc)
                      ctp_mat_local(dim, oloc, jloc) = ctp_mat_local(dim, oloc, jloc) &
                                                      + face_value*(-normgi(dim))/norm_ele(oloc) ! notvisited

                    end do inner_dimension_loop

                  end do nodal_loop_j

                end if ! notvisited

              end do quadrature_loop

            end if
          end do face_loop
        end do nodal_loop_i

        outer_dimension_loop: do dim = 1, size(normgi)
            call addto(CTP_m, 1, dim, p_nodes, u_nodes, ctp_mat_local(dim,:,:))
        end do outer_dimension_loop

      end do element_loop

      x_cvbdyshape=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape%degree)
      u_cvbdyshape=make_cvbdy_element_shape(cvfaces, u%mesh%faces%shape%degree)
      dens_cvbdyshape=make_cvbdy_element_shape(cvfaces, dens%mesh%faces%shape%degree)

      allocate(x_ele_bdy(x%dim,face_loc(x,1)), &
               detwei_bdy(x_cvbdyshape%ngi), &
               normal_bdy(x%dim, x_cvbdyshape%ngi), &
               u_bdy_f(u%dim, u_cvbdyshape%ngi), &
               dens_ele_bdy(face_loc(dens,1)), &
               olddens_ele_bdy(face_loc(dens,1)), &
               ghost_dens_ele_bdy(face_loc(dens,1)), &
               ghost_olddens_ele_bdy(face_loc(dens,1)), &
               norm_ele_bdy(face_loc(normalisation,1)))
      allocate(dens_bc_type(surface_element_count(dens)), &
               u_nodes_bdy(face_loc(u,1)), &
               p_nodes_bdy(face_loc(p,1)), &
               velocity_bc_type(u%dim, surface_element_count(u)), &
               velocity_bc_val(u%dim, u%mesh%faces%shape%loc))
      allocate(ctp_mat_local_bdy(x%dim, p%mesh%faces%shape%loc, u_cvbdyshape%loc), &
               ct_rhs_local(p%mesh%faces%shape%loc))

      call get_entire_boundary_condition(dens, &
                                         (/"weakdirichlet"/), &
                                         dens_bc, dens_bc_type)

      call get_entire_boundary_condition(u, (/"weakdirichlet ", &
                                              "no_normal_flow", &
                                              "periodic      "/), &
                                             velocity_bc, velocity_bc_type)

      surface_element_loop: do sele = 1, surface_element_count(p)

        if(any(velocity_bc_type(:,sele)==2).or.any(velocity_bc_type(:,sele)==3)) cycle

        ele = face_ele(x, sele)
        x_ele = ele_val(x, ele)
        x_ele_bdy = face_val(x, sele)
        u_nodes_bdy=face_global_nodes(u, sele)
        p_nodes_bdy=face_global_nodes(p, sele)

        call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, &
                              x_cvbdyshape, normal_bdy, detwei_bdy)

        u_bdy_f=face_val_at_quad(relu, sele, u_cvbdyshape)

        if(any(velocity_bc_type(:, sele)==1)) then
          velocity_bc_val = ele_val(velocity_bc, sele)
        else
          velocity_bc_val = 0.0
        end if

        if(dens_bc_type(sele)==1) then
          ghost_dens_ele_bdy=ele_val(dens_bc, sele)
        else
          ghost_dens_ele_bdy=face_val(dens, sele)
        end if

        if(dens_bc_type(sele)==1) then
          ghost_olddens_ele_bdy=ele_val(dens_bc, sele) ! not considering time varying bcs yet - unused
        else
          ghost_olddens_ele_bdy=face_val(olddens, sele) ! - unused
        end if

        dens_ele_bdy=face_val(dens, sele)
        olddens_ele_bdy=face_val(olddens, sele)

        norm_ele_bdy=face_val(normalisation, sele)

        ctp_mat_local_bdy = 0.0
        ct_rhs_local = 0.0

        surface_nodal_loop_i: do iloc = 1, p%mesh%faces%shape%loc

          surface_face_loop: do face = 1, cvfaces%sfaces

            if(cvfaces%sneiloc(iloc,face)/=0) then

              surface_quadrature_loop: do gi = 1, cvfaces%shape%ngi

                ggi = (face-1)*cvfaces%shape%ngi + gi

                udotn=dot_product(u_bdy_f(:,ggi), normal_bdy(:,ggi))

                if(udotn>0) then
                  income=0.0
                else
                  income=1.0
                end if

                face_value = (income*ghost_dens_ele_bdy(iloc) + (1.-income)*dens_ele_bdy(iloc))/&
                              norm_ele_bdy(iloc)

                surface_nodal_loop_j: do jloc = 1, u_cvbdyshape%loc

                  surface_inner_dimension_loop: do dim = 1, size(normal_bdy,1)
                  
                    if((present(ct_rhs)).and.(velocity_bc_type(dim, sele)==1)) then
                    
                      ct_rhs_local(iloc) = ct_rhs_local(iloc) + &
                                  face_value*u_cvbdyshape%n(jloc,ggi)*detwei_bdy(ggi)*normal_bdy(dim,ggi)*velocity_bc_val(dim,jloc)
                    
                    else

                      ctp_mat_local_bdy(dim, iloc, jloc) = ctp_mat_local_bdy(dim, iloc, jloc) &
                                  + face_value*u_cvbdyshape%n(jloc,ggi)*detwei_bdy(ggi)*normal_bdy(dim, ggi)
                    
                    end if

                  end do surface_inner_dimension_loop

                end do surface_nodal_loop_j

              end do surface_quadrature_loop

            end if ! sneiloc

          end do surface_face_loop

        end do surface_nodal_loop_i

        surface_outer_dimension_loop: do dim = 1, size(normal_bdy,1)

          if((present(ct_rhs)).and.(velocity_bc_type(dim, sele)==1)) then
          
            call addto(ct_rhs, p_nodes_bdy, ct_rhs_local)
          
          else

            call addto(CTP_m, 1, dim, p_nodes_bdy, u_nodes_bdy, ctp_mat_local_bdy(dim,:,:))
          
          end if

        end do surface_outer_dimension_loop

      end do surface_element_loop

      call deallocate(velocity_bc)
      deallocate(velocity_bc_type)

      call deallocate(x_cvbdyshape)
      call deallocate(u_cvbdyshape)
      call deallocate(dens_cvbdyshape)
      deallocate(x_ele_bdy, detwei_bdy, normal_bdy, u_bdy_f)
      deallocate(u_nodes_bdy, p_nodes_bdy)
      deallocate(dens_ele_bdy, olddens_ele_bdy, norm_ele_bdy)
      deallocate(ghost_dens_ele_bdy, ghost_olddens_ele_bdy)
      call deallocate(dens_bc)
      deallocate(dens_bc_type)
      deallocate(ctp_mat_local_bdy)

      call deallocate(x_cvshape)
      call deallocate(u_cvshape)
      call deallocate(dens_cvshape)
      call deallocate(cvfaces)
      call deallocate(relu)
      deallocate(x_ele, x_f, detwei, normal, normgi, u_f)
      deallocate(cfl_ele, dens_ele, olddens_ele, norm_ele)
      deallocate(notvisited)

      call deallocate(dens_upwind)
      call deallocate(olddens_upwind)
      call deallocate(cfl_no)
      if(norm_stat/=0) then
        call deallocate(normalisation)
        deallocate(normalisation)
      end if
      call deallocate(x_p)

    end subroutine assemble_1mat_compressible_divergence_matrix_cv
    !************************************************************************

    subroutine assemble_mmat_compressible_divergence_matrix_cv(CTP_m, state)

      ! inputs/outputs
      ! bucket full of fields
      type(state_type), dimension(:), intent(inout) :: state

      ! the pressure gradient and compressible gradient matrices
      type(block_csr_matrix), intent(inout) :: CTP_m

      ! local
      ! degree of quadrature over cv faces
      integer :: quaddegree

      real, dimension(:,:), allocatable :: x_ele, x_ele_bdy
      real, dimension(:,:), allocatable :: x_f, u_f, u_bdy_f
      real, dimension(:,:), allocatable :: normal, normal_bdy
      real, dimension(:), allocatable :: detwei, detwei_bdy
      real, dimension(:), allocatable :: normgi
      integer, dimension(:), pointer :: nodes, x_nodes
      integer, dimension(:), allocatable :: nodes_bdy
      real :: matvfrac_theta_val, matdens_theta_val
      real :: matvfrac_face_val, matdens_face_val
      real :: oldmatvfrac_face_val, oldmatdens_face_val
      real, dimension(:), allocatable :: matdens_ele, oldmatdens_ele, &
                                         matvfrac_ele, oldmatvfrac_ele, &
                                         norm_ele
      real, dimension(:), allocatable :: matdens_ele_bdy, oldmatdens_ele_bdy, &
                                         matvfrac_ele_bdy, oldmatvfrac_ele_bdy, &
                                         norm_ele_bdy
      real, dimension(:), allocatable :: ghost_matdens_ele_bdy, ghost_oldmatdens_ele_bdy, &
                                         ghost_matvfrac_ele_bdy, ghost_oldmatvfrac_ele_bdy

      integer, dimension(:), allocatable :: visited

      ! loop integers
      integer :: ele, sele, iloc, oloc, jloc, face, gi, ggi, dim, i

      ! information about cv faces
      type(cv_faces_type) :: cvfaces
      ! shape functions for region and surface
      type(element_type) :: x_cvshape, x_cvbdyshape
      type(element_type) :: u_cvshape, u_cvbdyshape
      type(element_type) :: p_cvshape, p_cvbdyshape
      ! pointer to velocity
      type(vector_field), pointer :: u, ug, x
      type(vector_field) :: relu, x_p
      type(scalar_field), pointer :: p
      type(scalar_field) :: cfl_no

      type(csr_sparsity), pointer :: mesh_sparsity
      type(csr_matrix)  :: matvfrac_upwind, &
            oldmatvfrac_upwind, matdens_upwind, oldmatdens_upwind, &
            summatvfrac_upwind, sumoldmatvfrac_upwind

      type(scalar_field), pointer :: matvfrac, oldmatvfrac, matdens, oldmatdens
      type(scalar_field), pointer :: dummyvfrac
      type(scalar_field) :: summatvfrac_bc
      integer, dimension(:), allocatable :: dummyvfrac_bc_type

      real, dimension(:), allocatable :: cfl_ele

      character(len=OPTION_PATH_LEN) :: vfrac_option_path
      integer :: nmatdens, dstat, vstat, stat
      real :: face_value
      character(len=FIELD_NAME_LEN) :: cfl_type

      integer, dimension(:), allocatable :: matvfrac_bc_type, matdens_bc_type
      type(scalar_field) :: matvfrac_bc, matdens_bc

      real :: udotn, income
      logical :: inflow
      real :: dt

      type(cv_options_type) :: matvfrac_options
      type(cv_options_type) :: matdens_options

      type(scalar_field), pointer :: normalisation, dummyones
      character(len=FIELD_NAME_LEN) :: normalisation_field
      integer :: norm_stat

      integer, dimension(:,:), allocatable :: velocity_bc_type
      type(vector_field) :: velocity_bc

      type(mesh_type), pointer :: bc_mesh

      real, dimension(:,:,:), allocatable :: ctp_mat_local, ctp_mat_local_bdy

      ! =============================================================
      ! Subroutine to construct the matrix CTP_m (a.k.a. C1/2/3TP).
      ! =============================================================

      ewrite(2,*) 'In assemble_mmat_compressible_divergence_matrix_cv'

      u=>extract_vector_field(state(1), "IteratedVelocity")
      ug=>extract_vector_field(state(1), "GridVelocity")
      p=>extract_scalar_field(state(1), "Pressure")
      x=>extract_vector_field(state(1), "Coordinate")
      x_p = get_coordinate_field(state(1), p%mesh)

      call allocate(relu, u%dim, u%mesh, "RelativeVelocity")

      call set(relu, u)
      call addto(relu, ug, scale=-1.0)

      mesh_sparsity=>get_csr_sparsity_firstorder(state, p%mesh, p%mesh)

      nmatdens=option_count("/material_phase/scalar_field::MaterialDensity")
      call get_option("/timestepping/timestep", dt)

      call get_option("/geometry/quadrature/controlvolume_surface_degree", &
                     quaddegree, default=1)

      cvfaces=find_cv_faces(vertices=ele_vertices(p, 1), &
                            dimension=mesh_dim(p), &
                            polydegree=p%mesh%shape%degree, &
                            quaddegree=quaddegree)

      x_cvshape=make_cv_element_shape(cvfaces, x%mesh%shape%degree)
      u_cvshape=make_cv_element_shape(cvfaces, u%mesh%shape%degree)
      p_cvshape=make_cv_element_shape(cvfaces, p%mesh%shape%degree)

      allocate(x_ele(x%dim,ele_loc(x,1)), &
               x_f(x%dim, x_cvshape%ngi), &
               u_f(u%dim, u_cvshape%ngi), &
               detwei(x_cvshape%ngi), &
               normal(x%dim, x_cvshape%ngi), &
               normgi(x%dim))
      allocate(cfl_ele(ele_loc(p,1)), &
               matvfrac_ele(ele_loc(p,1)), &
               oldmatvfrac_ele(ele_loc(p,1)), &
               matdens_ele(ele_loc(p,1)), &
               oldmatdens_ele(ele_loc(p,1)), &
               norm_ele(ele_loc(x,1)))
      allocate(visited(x_cvshape%ngi))
      allocate(ctp_mat_local(x%dim, p%mesh%shape%loc, u_cvshape%loc))

      call allocate(summatvfrac_upwind, mesh_sparsity, name="SumMatVFracUpwindValues")
      call allocate(sumoldmatvfrac_upwind, mesh_sparsity, name="SumOldMatVFracUpwindValues")
      call allocate(matvfrac_upwind, mesh_sparsity, name="VFracUpwindValues")
      call allocate(oldmatvfrac_upwind, mesh_sparsity, name="OldVFracUpwindValues")
      call allocate(matdens_upwind, mesh_sparsity, name="MatDensUpwindValues")
      call allocate(oldmatdens_upwind, mesh_sparsity, name="OldMatDensUpwindValues")

      allocate(dummyvfrac)
      call allocate(dummyvfrac, p%mesh, name="DummyVFrac", field_type=FIELD_TYPE_CONSTANT)
      call set(dummyvfrac, 1.0)

      allocate(dummyones)
      call allocate(dummyones, p%mesh, name="DummyOnes", field_type=FIELD_TYPE_CONSTANT)
      call set(dummyones, 1.0)

      call get_option(trim(p%option_path)//"/prognostic/scheme/use_compressible_projection_method/normalisation/name", &
                      normalisation_field, stat=norm_stat)

      allocate(dummyvfrac_bc_type(surface_element_count(dummyvfrac)))
      bc_mesh=>get_dg_surface_mesh(p%mesh)
      call allocate(summatvfrac_bc, bc_mesh, name="SumVolumeFractionsBCs")

      call zero(summatvfrac_upwind)
      call zero(sumoldmatvfrac_upwind)
      call zero(summatvfrac_bc)
      dummyvfrac_bc_type=0

      vfrac_option_path = " "
      do i = 1, size(state)

        matvfrac=>extract_scalar_field(state(i), "MaterialVolumeFraction", stat=vstat)
        if(vstat==0) then
          if((.not.aliased(matvfrac)).and.(.not.have_option(trim(matvfrac%option_path)//"/diagnostic"))) then

            vfrac_option_path = trim(matvfrac%option_path)

            if(need_upwind_values(trim(matvfrac%option_path))) then
              oldmatvfrac=>extract_scalar_field(state(i), "OldMaterialVolumeFraction")

              call find_upwind_values(state, x_p, matvfrac, matvfrac_upwind, &
                                    oldmatvfrac, oldmatvfrac_upwind, defer_deletion=.true.)
              summatvfrac_upwind%val=summatvfrac_upwind%val+matvfrac_upwind%val

              sumoldmatvfrac_upwind%val=sumoldmatvfrac_upwind%val+oldmatvfrac_upwind%val
            end if

            allocate(matvfrac_bc_type(surface_element_count(matvfrac)))
            call get_entire_boundary_condition(matvfrac, (/"weakdirichlet"/), matvfrac_bc, matvfrac_bc_type)

            assert(all(summatvfrac_bc%mesh%ndglno==matvfrac_bc%mesh%ndglno))  ! make sure fields are on the same mesh

            summatvfrac_bc%val=summatvfrac_bc%val+matvfrac_bc%val
            dummyvfrac_bc_type=matvfrac_bc_type ! we assume the bcs are the same on all volume fractions

            call deallocate(matvfrac_bc)
            deallocate(matvfrac_bc_type)
          end if

        end if

      end do

      ! find courant number (if needed) - here we assume that the courant number is associated with
      ! the MaterialVolumeFraction field
      ! FIXME: this will not work with material specific measures of the courant number (like CVMaterialDensityCFLNumber)
      ! FIXME: don't make me assume I'm on the volume fraction
      call allocate(cfl_no, p%mesh, "CourantNumber")
      call get_option(trim(complete_cv_field_path(vfrac_option_path))//&
                             "/face_value[0]/courant_number[0]/name", &
                            cfl_type, stat)
      if (stat==0) then
         call calculate_diagnostic_variable(state(1), trim(cfl_type), cfl_no)
      else
         call set(cfl_no, 1.0)
      end if

      ! Clear memory of arrays being designed
      call zero(CTP_m)

      do i = 1, size(state)
        matdens=>extract_scalar_field(state(i), "MaterialDensity", stat=dstat)

        if(dstat==0) then
          oldmatdens=>extract_scalar_field(state(i), "OldMaterialDensity")

          if(need_upwind_values(trim(matdens%option_path))) then

            call find_upwind_values(state, x_p, matdens, matdens_upwind, &
                                    oldmatdens, oldmatdens_upwind, defer_deletion=.true.)

          else

            call zero(matdens_upwind)
            call zero(oldmatdens_upwind)

          end if

          ! get all the relevent options for material density
          ! handily wrapped in a new type...
          matdens_options = get_cv_options(matdens%option_path, matdens%mesh%shape%numbering%family)

          matvfrac=>extract_scalar_field(state(i), "MaterialVolumeFraction", stat=vstat)
          if(vstat==0) then
            oldmatvfrac=>extract_scalar_field(state(i), "OldMaterialVolumeFraction")
          else
            if(nmatdens==1) then

              matvfrac=>dummyvfrac
              oldmatvfrac=>dummyvfrac

            else
              ewrite(-1,*) "Multiple MaterialDensities but at least "
              ewrite(-1,*) "one has no associated MaterialVolumeFraction"
              FLAbort("This shouldn't happen")  ! move this to a multimaterials check options
            end if
          end if

          if(need_upwind_values(trim(matvfrac%option_path))) then

            call find_upwind_values(state, x_p, matvfrac, matvfrac_upwind, &
                                    oldmatvfrac, oldmatvfrac_upwind, defer_deletion=.true.)

            vfrac_option_path=trim(matvfrac%option_path)

          else

            matvfrac_upwind%val = 1.0 - summatvfrac_upwind%val
            oldmatvfrac_upwind%val = 1.0 - sumoldmatvfrac_upwind%val

          end if

          ! get all the relevent options for material volume fraction
          ! handily wrapped in a new type...
          if((size(state)>1)) then
            ! only possible if we have more than 1 state (since with 1 state the volume fraction won't
            ! be prognostic and won't have the relevant options)
            matvfrac_options = get_cv_options(vfrac_option_path, matvfrac%mesh%shape%numbering%family)
          end if

          ! get the normalisation field (if we need one)
          if(norm_stat==0) then
            normalisation=>extract_scalar_field(state(i), trim(normalisation_field))
          else
            normalisation=>dummyones
          end if

          do ele=1, element_count(p)
            x_ele=ele_val(x, ele)
            x_f=ele_val_at_quad(x, ele, x_cvshape)
            u_f=ele_val_at_quad(relu, ele, u_cvshape)
            nodes=>ele_nodes(u, ele)
            x_nodes=>ele_nodes(x_p, ele)

            call transform_cvsurf_to_physical(x_ele, x_cvshape, &
                                              detwei, normal, cvfaces)
            cfl_ele = ele_val(cfl_no, ele)

            matvfrac_ele = ele_val(matvfrac, ele)
            oldmatvfrac_ele = ele_val(oldmatvfrac, ele)

            matdens_ele = ele_val(matdens, ele)
            oldmatdens_ele = ele_val(oldmatdens, ele)

            norm_ele = ele_val(normalisation, ele)

            visited=0

            ctp_mat_local = 0.0

            do iloc = 1, p%mesh%shape%loc

              do face = 1, cvfaces%faces

                if(cvfaces%neiloc(iloc, face) /= 0) then
                  oloc = cvfaces%neiloc(iloc, face)

                  do gi = 1, cvfaces%shape%ngi

                    ggi = (face-1)*cvfaces%shape%ngi + gi

                    if(visited(ggi)==0) then
                      visited(ggi)=1

                      normgi=orientate_cvsurf_normgi(x_ele(:,iloc),x_f(:,ggi),normal(:,ggi))

                      udotn=dot_product(u_f(:,ggi), normgi(:))

                      inflow = (udotn<=0.0)

                      income = merge(1.0,0.0,inflow)

                     select case (matvfrac%field_type)
                     case(FIELD_TYPE_CONSTANT)

                        matvfrac_face_val = matvfrac_ele(iloc)
                        oldmatvfrac_face_val = oldmatvfrac_ele(iloc)

                     case default

                        if(size(state)>1) then
                          call evaluate_face_val(matvfrac_face_val, oldmatvfrac_face_val, & 
                                                  iloc, oloc, ggi, x_nodes, &
                                                  p_cvshape, &
                                                  matvfrac_ele, oldmatvfrac_ele, &
                                                  matvfrac_upwind, oldmatvfrac_upwind, &
                                                  inflow, cfl_ele, &
                                                  matvfrac_options)
                        else
                          matvfrac_face_val = matvfrac_ele(iloc)
                          oldmatvfrac_face_val = oldmatvfrac_ele(iloc)
                        end if

                      end select

                      call evaluate_face_val(matdens_face_val, oldmatdens_face_val, &
                                             iloc, oloc, ggi, x_nodes, &
                                             p_cvshape,&
                                             matdens_ele, oldmatdens_ele, &
                                             matdens_upwind, oldmatdens_upwind, &
                                             inflow, cfl_ele, &
                                             matdens_options)

                      if(size(state)>1) then
                        matvfrac_theta_val=theta_val(iloc, oloc, &
                                                      matvfrac_face_val, &
                                                      oldmatvfrac_face_val, &
                                                      matvfrac_options%theta, dt, udotn, &
                                                      x_ele, matvfrac_options%limit_theta, &
                                                      matvfrac_ele, oldmatvfrac_ele)
                      else
                        matvfrac_theta_val = matvfrac_face_val
                      end if
                      
                      matdens_theta_val=theta_val(iloc, oloc, &
                                                    matdens_face_val, &
                                                    oldmatdens_face_val, &
                                                    matdens_options%theta, dt, udotn, &
                                                    x_ele, matdens_options%limit_theta, &
                                                    matdens_ele, oldmatdens_ele)


                      do jloc = 1, u_cvshape%loc

                        face_value = u_cvshape%n(jloc, ggi)*detwei(ggi)*matvfrac_theta_val*matdens_theta_val

                        do dim = 1, size(normgi)

                          ctp_mat_local(dim, iloc, jloc) = ctp_mat_local(dim, iloc, jloc) &
                                                         + face_value*normgi(dim)/norm_ele(iloc)
                          ctp_mat_local(dim, oloc, jloc) = ctp_mat_local(dim, oloc, jloc) &
                                                         + face_value*(-normgi(dim))/norm_ele(oloc) ! notvisited

                        end do

                      end do

                    end if ! visited

                  end do

                end if
              end do
            end do

            do dim = 1, size(normgi)
                call addto(CTP_m, 1, dim, nodes, nodes, ctp_mat_local(dim,:,:))
            end do

          end do

          x_cvbdyshape=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape%degree)
          u_cvbdyshape=make_cvbdy_element_shape(cvfaces, u%mesh%faces%shape%degree)
          p_cvbdyshape=make_cvbdy_element_shape(cvfaces, p%mesh%faces%shape%degree)

          allocate(x_ele_bdy(x%dim,face_loc(x,1)), &
                  detwei_bdy(x_cvbdyshape%ngi), &
                  normal_bdy(x%dim, x_cvbdyshape%ngi), &
                  u_bdy_f(u%dim, u_cvbdyshape%ngi), &
                  matdens_ele_bdy(face_loc(p,1)), &
                  oldmatdens_ele_bdy(face_loc(p,1)), &
                  matvfrac_ele_bdy(face_loc(p,1)), &
                  oldmatvfrac_ele_bdy(face_loc(p,1)), &
                  ghost_matdens_ele_bdy(face_loc(p,1)), &
                  ghost_oldmatdens_ele_bdy(face_loc(p,1)), &
                  ghost_matvfrac_ele_bdy(face_loc(p,1)), &
                  ghost_oldmatvfrac_ele_bdy(face_loc(p,1)), &
                  norm_ele_bdy(face_loc(x,1)))
          allocate(matvfrac_bc_type(surface_element_count(matvfrac)), &
                    matdens_bc_type(surface_element_count(matdens)), &
                    nodes_bdy(face_loc(u,1)), &
                    velocity_bc_type(u%dim, surface_element_count(u)))
          allocate(ctp_mat_local_bdy(x%dim, p%mesh%faces%shape%loc, u_cvbdyshape%loc))

          if((.not.aliased(matvfrac)).and.(.not.have_option(trim(matvfrac%option_path)//"/diagnostic"))) then
            call get_entire_boundary_condition(matvfrac, (/"weakdirichlet"/), matvfrac_bc, matvfrac_bc_type)
          else
            call allocate(matvfrac_bc, summatvfrac_bc%mesh, "TemporaryMaterialVolumeFraction")
            matvfrac_bc%val = 1.0-summatvfrac_bc%val
            matvfrac_bc_type = dummyvfrac_bc_type
          end if
          call get_entire_boundary_condition(matdens, (/"weakdirichlet"/), matdens_bc, matdens_bc_type)

          call get_entire_boundary_condition(u, (/"weakdirichlet ", "no_normal_flow", "periodic      "/), velocity_bc, velocity_bc_type)

          do sele = 1, surface_element_count(p)

            if(any(velocity_bc_type(:,sele)==2).or.any(velocity_bc_type(:,sele)==3)) cycle

            ele = face_ele(x, sele)
            x_ele = ele_val(x, ele)
            x_ele_bdy = face_val(x, sele)
            nodes_bdy=face_global_nodes(u, sele)

            call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, &
                                  x_cvbdyshape, normal_bdy, detwei_bdy)

            u_bdy_f=face_val_at_quad(relu, sele, u_cvbdyshape)

            if(matvfrac_bc_type(sele)==1) then
              ghost_matvfrac_ele_bdy=ele_val(matvfrac_bc, sele)
            else
              ghost_matvfrac_ele_bdy=face_val(matvfrac, sele)
            end if

            if(matvfrac_bc_type(sele)==1) then
              ghost_oldmatvfrac_ele_bdy=ele_val(matvfrac_bc, sele) ! not considering time varying bcs yet - unused
            else
              ghost_oldmatvfrac_ele_bdy=face_val(oldmatvfrac, sele) ! - unused
            end if

            matvfrac_ele_bdy=face_val(matvfrac, sele)
            oldmatvfrac_ele_bdy=face_val(oldmatvfrac, sele)

            if(matdens_bc_type(sele)==1) then
              ghost_matdens_ele_bdy=ele_val(matdens_bc, sele)
            else
              ghost_matdens_ele_bdy=face_val(matdens, sele)
            end if

            if(matdens_bc_type(sele)==1) then
              ghost_oldmatdens_ele_bdy=ele_val(matdens_bc, sele) ! not considering time varying bcs yet - unused
            else
              ghost_oldmatdens_ele_bdy=face_val(oldmatdens, sele) ! - unused
            end if

            matdens_ele_bdy=face_val(matdens, sele)
            oldmatdens_ele_bdy=face_val(oldmatdens, sele)

            norm_ele_bdy=face_val(normalisation, sele)

            ctp_mat_local_bdy = 0.0

            do iloc = 1, p%mesh%faces%shape%loc

              do face = 1, cvfaces%sfaces

                if(cvfaces%sneiloc(iloc,face)/=0) then

                  do gi = 1, cvfaces%shape%ngi

                    ggi = (face-1)*cvfaces%shape%ngi + gi

                    udotn=dot_product(u_bdy_f(:,ggi), normal_bdy(:,ggi))

                    if(udotn>0) then
                      income=0.0
                    else
                      income=1.0
                    end if

                    face_value = (income*ghost_matvfrac_ele_bdy(iloc) + (1.-income)*matvfrac_ele_bdy(iloc))* &
                                  (income*ghost_matdens_ele_bdy(iloc) + (1.-income)*matdens_ele_bdy(iloc))/&
                                  norm_ele_bdy(iloc)

                    do jloc = 1, u_cvbdyshape%loc

                      do dim = 1, size(normal_bdy,1)

                        ctp_mat_local_bdy(dim, iloc, jloc) = ctp_mat_local_bdy(dim, iloc, jloc) &
                                    + face_value*u_cvbdyshape%n(jloc,ggi)*detwei_bdy(ggi)*normal_bdy(dim, ggi)

                      end do ! dim

                    end do ! jloc

                  end do ! gi

                end if ! sneiloc

              end do ! face

            end do ! iloc

            do dim = 1, size(normal_bdy,1)

              call addto(CTP_m, 1, dim, nodes_bdy, nodes_bdy, ctp_mat_local_bdy(dim,:,:))

            end do ! dim

          end do ! sele

          call deallocate(velocity_bc)
          deallocate(velocity_bc_type)

          call deallocate(x_cvbdyshape)
          call deallocate(u_cvbdyshape)
          call deallocate(p_cvbdyshape)
          deallocate(x_ele_bdy, detwei_bdy, normal_bdy, u_bdy_f)
          deallocate(nodes_bdy)
          deallocate(matdens_ele_bdy, oldmatdens_ele_bdy, matvfrac_ele_bdy, oldmatvfrac_ele_bdy, norm_ele_bdy)
          deallocate(ghost_matdens_ele_bdy, ghost_oldmatdens_ele_bdy, &
                      ghost_matvfrac_ele_bdy, ghost_oldmatvfrac_ele_bdy)
          call deallocate(matvfrac_bc)
          call deallocate(matdens_bc)
          deallocate(matvfrac_bc_type, matdens_bc_type)
          deallocate(ctp_mat_local_bdy)

        end if ! dstat==0

      end do ! i

      call deallocate(x_cvshape)
      call deallocate(u_cvshape)
      call deallocate(p_cvshape)
      call deallocate(cvfaces)
      call deallocate(relu)
      deallocate(x_ele, x_f, detwei, normal, normgi, u_f)
      deallocate(cfl_ele, matvfrac_ele, oldmatvfrac_ele, matdens_ele, oldmatdens_ele, norm_ele)
      deallocate(visited)

      call deallocate(matdens_upwind)
      call deallocate(oldmatdens_upwind)
      call deallocate(matvfrac_upwind)
      call deallocate(oldmatvfrac_upwind)
      call deallocate(summatvfrac_upwind)
      call deallocate(sumoldmatvfrac_upwind)
      call deallocate(summatvfrac_bc)
      call deallocate(cfl_no)
      call deallocate(dummyvfrac)
      deallocate(dummyvfrac)
      deallocate(dummyvfrac_bc_type)
      call deallocate(dummyones)
      deallocate(dummyones)
      call deallocate(x_p)

      call clean_deferred_deletion(state)

    end subroutine assemble_mmat_compressible_divergence_matrix_cv
    !************************************************************************

end module divergence_matrix_cv

