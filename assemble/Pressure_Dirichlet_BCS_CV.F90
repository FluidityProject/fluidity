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
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"

module pressure_dirichlet_bcs_cv

  use global_parameters, only: OPTION_PATH_LEN
  use quadrature
  use futils
  use spud
  use cv_faces
  use fetools
  use fields
  use state_module
  use boundary_conditions
  use field_derivatives
  use cv_shape_functions
  use field_options, only: get_coordinate_field
  use cvtools
  use cv_options
  use cv_upwind_values
  use cv_face_values
  use sparsity_patterns_meshes
  use diagnostic_fields, only: calculate_diagnostic_variable
  use cv_fields
  use multiphase_module

  implicit none

  private
  
  public :: add_pressure_dirichlet_bcs_cv

contains

    ! --------------------------------------------------------------------------------------

    subroutine add_pressure_dirichlet_bcs_cv(mom_rhs, u, p, state)
       
      !!< Add any CV pressure dirichlet BC integrals to the mom_rhs field if required.
      !!< If this is a multiphase simulation then the phase volume fraction spatial 
      !!< interpolation uses finite elements. This isnt striclty correct if the 
      !!< phase volume fractions were solved with a control volume discretisation 
      !!< and also whether of not they themselves have a weak BC is ignored. If the 
      !!< latter was to be considered then the upwind value should be used.
      
      ! Note the pressure BC values come from the end of time step Pressure field
      ! rather than using a theta average. The latter was tried but didnt work 
      ! because OldPressure had no BC info. Taking the end of time step value 
      ! for the BC is also what the pressure CG routine(s) do.
      
      type(vector_field), intent(inout) :: mom_rhs
      type(vector_field), intent(in) :: u
      type(scalar_field), intent(in) :: p
      type(state_type), intent(inout) :: state ! required for phase volume fraction

      ! local variables

      ! degree of quadrature over cv faces
      integer :: quaddegree

      real, dimension(:,:), allocatable :: x_ele, x_ele_bdy
      real, dimension(:,:), allocatable :: normal_bdy
      real, dimension(:), allocatable :: detwei_bdy
      integer, dimension(:), allocatable :: u_nodes_bdy

      integer :: ele, sele, iloc, jloc, face, gi, ggi, dim

      ! information about cv faces
      type(cv_faces_type) :: cvfaces
      ! shape functions for surface
      type(element_type) :: x_cvbdyshape
      type(element_type) :: u_cvbdyshape
      
      ! pointer to coordinates
      type(vector_field), pointer :: x
      
      ! pressure BC info
      integer, dimension(:), allocatable :: pressure_bc_type
      real, dimension(:), allocatable :: pressure_bc_val
      type(scalar_field) :: pressure_bc
      
      ! local ct matrix
      real, dimension(:,:,:), allocatable :: ct_mat_local_bdy

      ! Multiphase variables
      logical :: multiphase
      ! Volume fraction fields
      type(scalar_field), pointer :: vfrac
      type(scalar_field) :: nvfrac
      
      ! Variables used for FE interpolation of vfrac on surface
      ! Values of nvfrac at the faces 
      real, dimension(:), allocatable :: nvfrac_gi_f
      ! CV shape functions for nvfrac (computed only if the Coordinate and
      ! PhaseVolumeFraction meshes are different, otherwise it is
      ! assigned x_cvbdyshape)
      type(element_type) :: nvfrac_cvbdyshape

      ewrite(1,*) 'In add_cv_pressure_weak_dirichlet_bcs'

      x => extract_vector_field(state, "Coordinate")
      
      allocate(x_ele(x%dim, x%mesh%shape%loc))

      call get_option("/geometry/quadrature/controlvolume_surface_degree", quaddegree, default=1)

      cvfaces = find_cv_faces(vertices   = ele_vertices(p, 1), &
                              dimension  = mesh_dim(p), &
                              polydegree = p%mesh%shape%degree, &
                              quaddegree = quaddegree)

      x_cvbdyshape = make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape)
      u_cvbdyshape = make_cvbdy_element_shape(cvfaces, u%mesh%faces%shape)
   
      assert(surface_element_count(p) == surface_element_count(u))
      
      ! get the pressure BC field
      allocate(pressure_bc_type(surface_element_count(p)))
      
      call get_entire_boundary_condition(p, (/"weakdirichlet", &
                                              "dirichlet    "/), pressure_bc, pressure_bc_type)
      
      allocate(x_ele_bdy(x%dim,x%mesh%faces%shape%loc), &
               detwei_bdy(x_cvbdyshape%ngi), &
               normal_bdy(x%dim, x_cvbdyshape%ngi), &
               pressure_bc_val(pressure_bc%mesh%shape%loc))
      
      allocate(u_nodes_bdy(u%mesh%faces%shape%loc))
            
      allocate(ct_mat_local_bdy(x%dim, p%mesh%faces%shape%loc, u%mesh%faces%shape%loc))

      ! Check if we need to multiply through by the non-linear volume fraction
      if(option_count("/material_phase/vector_field::Velocity/prognostic") > 1) then
         multiphase = .true.

         vfrac => extract_scalar_field(state, "PhaseVolumeFraction")
         call allocate(nvfrac, vfrac%mesh, "NonlinearPhaseVolumeFraction")
         call zero(nvfrac)
         call get_nonlinear_volume_fraction(state, nvfrac)
         
         assert(surface_element_count(p) == surface_element_count(vfrac))
         
         ewrite_minmax(nvfrac)
         
         ! If the Coordinate and PhaseVolumeFraction meshes are different, then we need to
         ! generate the PhaseVolumeFraction CV shape functions.
         if(.not.(nvfrac%mesh == x%mesh)) then
            nvfrac_cvbdyshape = make_cvbdy_element_shape(cvfaces, nvfrac%mesh%faces%shape)
         else
            nvfrac_cvbdyshape = x_cvbdyshape
            ! incease reference for deallocates below
            call incref(nvfrac_cvbdyshape)
         end if

         allocate(nvfrac_gi_f(nvfrac_cvbdyshape%ngi))
                  
      else
         multiphase = .false.
         nullify(vfrac)
      end if
  
      surface_element_loop: do sele = 1, surface_element_count(p)
  
        ! cycle if this not a dirichlet pressure BC.
        if (pressure_bc_type(sele) == 0) cycle
            
        ele         = face_ele(x, sele)
        x_ele       = ele_val(x, ele)
        x_ele_bdy   = face_val(x, sele)
        u_nodes_bdy = face_global_nodes(u, sele)
        
        ! Get the phase volume fraction face value 
        if(multiphase) then
           nvfrac_gi_f = face_val_at_quad(nvfrac, sele, nvfrac_cvbdyshape)
        end if
  
        call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, &
                              x_cvbdyshape, normal_bdy, detwei_bdy)
  
        ct_mat_local_bdy = 0.0
        
        ! calculate the ct local matrix
        surface_nodal_loop_i: do iloc = 1, p%mesh%faces%shape%loc
  
           surface_face_loop: do face = 1, cvfaces%sfaces
          
              if(cvfaces%sneiloc(iloc,face) /= 0) then
  
                 surface_quadrature_loop: do gi = 1, cvfaces%shape%ngi
  
                    ggi = (face-1)*cvfaces%shape%ngi + gi
                                        
                    surface_nodal_loop_j: do jloc = 1, u%mesh%faces%shape%loc
                       
                       surface_inner_dimension_loop: do dim = 1, size(normal_bdy,1)

                          if(multiphase) then
                             ct_mat_local_bdy(dim, iloc, jloc) =  ct_mat_local_bdy(dim, iloc, jloc) + &
                                   u_cvbdyshape%n(jloc,ggi)*detwei_bdy(ggi)*nvfrac_gi_f(ggi)*normal_bdy(dim, ggi)
                          else
                             ct_mat_local_bdy(dim, iloc, jloc) =  ct_mat_local_bdy(dim, iloc, jloc) + &
                                   u_cvbdyshape%n(jloc,ggi)*detwei_bdy(ggi)*normal_bdy(dim, ggi)
                          end if
    
                       end do surface_inner_dimension_loop
  
                    end do surface_nodal_loop_j
  
                 end do surface_quadrature_loop
  
              end if
  
           end do surface_face_loop
  
        end do surface_nodal_loop_i
        
        ! pressure dirichlet BC is -c*press_bc_val = -press_bc_val*ct, integrated over surface elements appropriate
        surface_outer_dimension_loop: do dim = 1, size(normal_bdy,1)

          ! for weak and strong pressure dirichlet bcs:
          !      /
          ! add -|  N_i M_j \vec n p_j, where p_j are the prescribed bc values
          !      /
          
          call addto(mom_rhs, dim, u_nodes_bdy, -matmul(ele_val(pressure_bc, sele), ct_mat_local_bdy(dim,:,:)))
  
        end do surface_outer_dimension_loop
  
      end do surface_element_loop
  
      call deallocate(pressure_bc)
      deallocate(pressure_bc_type)
      deallocate(pressure_bc_val)
      call deallocate(x_cvbdyshape)
      call deallocate(u_cvbdyshape)
      deallocate(x_ele_bdy, detwei_bdy, normal_bdy)
      
      deallocate(u_nodes_bdy)
      deallocate(ct_mat_local_bdy)
        
      call deallocate(cvfaces)
      deallocate(x_ele)

      if(multiphase) then
         call deallocate(nvfrac)
         deallocate(nvfrac_gi_f)
         call deallocate(nvfrac_cvbdyshape)         
      end if
       
    end subroutine add_pressure_dirichlet_bcs_cv
    
    ! --------------------------------------------------------------------------------------

end module pressure_dirichlet_bcs_cv

