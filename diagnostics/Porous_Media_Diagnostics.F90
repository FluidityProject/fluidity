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

module porous_media_diagnostics
   
   !!< Module containing porous media related diagnostic algorithms.

   use global_parameters, only:FIELD_NAME_LEN, current_time, OPTION_PATH_LEN
   use fields
   use halos
   use field_derivatives
   use field_options
   use state_module
   use futils
   use fetools
   use fefields, only: compute_lumped_mass, compute_cv_mass
   use spud
   use CV_Shape_Functions
   use CV_Faces
   use CVTools
   use CV_Upwind_Values
   use CV_Face_Values
   use cv_options
   use parallel_tools
   use sparsity_patterns
   use sparsity_patterns_meshes
   use boundary_conditions, only: get_entire_boundary_condition
   use state_fields_module

   implicit none
   
   private

   public :: calculate_interstitial_velocity_dg_courant_number, &
             calculate_interstitial_velocity_cv_courant_number
  
contains

! ----------------------------------------------------------------------------------

   subroutine calculate_interstitial_velocity_dg_courant_number(state, s_field)
      
      !!< Calculate the interstitial velocity DG courant number field

      type(state_type), intent(inout) :: state
      type(scalar_field), intent(inout) :: s_field 
      
      ! local variables
      type(vector_field), pointer :: u, x
      real :: dt
      integer :: ele, stat
      ! Porosity fields, field name, theta value and include flag
      type(scalar_field), pointer :: porosity_old, porosity_new
      type(scalar_field) :: porosity_theta 
      character(len=OPTION_PATH_LEN) :: porosity_name
      real :: porosity_theta_value

      ewrite(1,*) 'Entering calculate_interstitial_velocity_dg_courant_number'

      u => extract_vector_field(state, "NonlinearVelocity",stat)
      
      if(stat.ne.0) then    
         u => extract_vector_field(state, "Velocity",stat)
         
         if(stat.ne.0) then
            FLExit('Failed to extract Velocity field for calculating interstitial_velocity_dg_courant_number')
         end if
      end if
      
      x => extract_vector_field(state, "Coordinate")
      
      call get_option("/timestepping/timestep", dt)
    
      ! determine the porosity value to include
     
      ! get the name of the field to use as porosity
      call get_option(trim(s_field%option_path)//'/diagnostic/algorithm/porosity_field_name', &
                      porosity_name, &
                      default = 'Porosity')
         
      ! get the porosity theta value
      call get_option(trim(s_field%option_path)//'/diagnostic/algorithm/temporal_discretisation/theta', &
                      porosity_theta_value, &
                      default = 0.0)
         
      porosity_new => extract_scalar_field(state, trim(porosity_name), stat = stat)                  

      if (stat /=0) then 
         FLExit('Failed to extract Porosity from state for calculating interstitial_velocity_dg_courant_number')
      end if
       
      porosity_old => extract_scalar_field(state, "Old"//trim(porosity_name), stat = stat)

      if (stat /=0) then 
         FLExit('Failed to extract OldPorosity from state for calculating interstitial_velocity_dg_courant_number')
      end if
       
      call allocate(porosity_theta, porosity_new%mesh)
         
      call set(porosity_theta, porosity_new, porosity_old, porosity_theta_value)
         
      ewrite_minmax(porosity_theta)
    
      call zero(s_field)
    
      do ele = 1, element_count(s_field)
         assert(ele_ngi(s_field, ele)==ele_ngi(porosity_theta, ele))
         call dg_courant_number_ele(s_field,x,u,ele,dt,ele_val_at_quad(porosity_theta,ele))
      end do
    
      call deallocate(porosity_theta)
    
      ! the courant values at the edge of the halo are going to be incorrect
      ! this matters when computing the max courant number
      call halo_update(s_field)
      
      ewrite(1,*) 'Exiting calculate_interstitial_velocity_dg_courant_number'

   end subroutine calculate_interstitial_velocity_dg_courant_number

! ----------------------------------------------------------------------------------

   subroutine dg_courant_number_ele(courant,x,u,ele,dt,porosity_theta_at_quad)
      
      !!< Calculate the DG courant number for the interstitial velocity for element ele.
      
      type(scalar_field), intent(inout) :: courant
      type(vector_field), intent(in) :: x, u
      real, intent(in) :: dt
      integer, intent(in) :: ele
      real, dimension(:) :: porosity_theta_at_quad
      
      real :: Vol
      real :: Flux
      integer :: ni, ele_2, face, face_2
      integer, dimension(:), pointer :: neigh
      real, dimension(ele_ngi(u,ele)) :: detwei
      real, dimension(face_ngi(u,1)) :: detwei_f
      real, dimension(U%dim, face_ngi(U, 1)) :: normal, U_f_quad
      real, dimension(face_ngi(U,1)) :: flux_quad
      integer, dimension(:), pointer :: u_ele
      real, dimension(ele_loc(u,ele)) :: Vols
      real :: val
      real, dimension(ele_loc(u,ele)) :: Vals
      
      !Get element volume
      call transform_to_physical(X, ele, detwei=detwei)
      Vol = sum(detwei*porosity_theta_at_quad)
    
      !Get fluxes
      Flux = 0.0
      neigh=>ele_neigh(U, ele)
      do ni = 1, size(neigh)
         ele_2=neigh(ni)
         face=ele_face(U, ele, ele_2)
         if(ele_2<0.0) then
            face_2 = face
         else
            face_2=ele_face(U, ele_2, ele)
         end if
       
         U_f_quad =0.5*(face_val_at_quad(U, face)&
              &     +face_val_at_quad(U, face_2))

         call transform_facet_to_physical(X, face, &
              &                          detwei_f=detwei_f,&
              &                          normal=normal) 

         Flux_quad = -sum(U_f_quad*normal,1)
         Flux_quad = max(Flux_quad,0.0)

         Flux = Flux + sum(Flux_quad*detwei_f)
      end do

      u_ele => ele_nodes(U,ele)

      Val = Flux/Vol*dt
      Vals = Val
      call set(Courant,U_ele,Vals)

  end subroutine dg_courant_number_ele

! ----------------------------------------------------------------------------------

   subroutine calculate_interstitial_velocity_cv_courant_number(state, s_field)
      
      !!< Calculate the interstitial velocity CV courant number field

      type(state_type), intent(inout) :: state
      type(scalar_field), intent(inout) :: s_field 
      
      ! local variables      
      type(vector_field), pointer :: u, x
      real :: dt
      integer :: i, ele, iloc, oloc, gi, ggi, sele, face, ni, face_2

      integer :: quaddegree
      type(cv_faces_type) :: cvfaces
      type(element_type) :: x_cvshape, x_cvbdyshape
      type(element_type) :: u_cvshape, u_cvbdyshape
      type(scalar_field), pointer :: cvmass
      real, dimension(:,:), allocatable :: x_ele, x_ele_bdy
      real, dimension(:,:), allocatable :: x_f, u_f, u_bdy_f
      real, dimension(:,:), allocatable :: normal, normal_bdy
      real, dimension(:), allocatable :: detwei, detwei_bdy
      real, dimension(:), allocatable :: normgi
      integer, dimension(:), pointer :: nodes, x_nodes, neigh
      integer, dimension(:), allocatable :: nodes_bdy
      real :: udotn, income
      ! logical array indicating if a face has already been visited by the opposing node
      logical, dimension(:), allocatable :: notvisited
      
      ! Porosity fields, field name, theta value and include flag
      type(scalar_field), pointer :: porosity_old, porosity_new
      type(scalar_field) :: porosity_theta 
      character(len=OPTION_PATH_LEN) :: porosity_name
      real :: porosity_theta_value
      logical :: include_porosity
      
      integer :: stat
      
      integer, dimension(:), allocatable :: courant_bc_type
      type(scalar_field) :: courant_bc

      type(vector_field) :: x_courant ! coordinates on s_field mesh
      
      ewrite(1,*) 'Entering calculate_interstitial_velocity_cv_courant_number'
      
      udotn = 0.0 

      u=>extract_vector_field(state, "NonlinearVelocity")
      x=>extract_vector_field(state, "Coordinate")
      
      call get_option("/timestepping/timestep",dt)

      call zero(s_field)

      x_courant=get_coordinate_field(state, s_field%mesh)
      
      ! determine the cv mass matrix to use for the length scale
      ! which will have included the porosity field
     
      ! get the name of the field to use as porosity
      call get_option(trim(s_field%option_path)//'/diagnostic/algorithm/porosity_field_name', &
                      porosity_name, &
                      default = 'Porosity')
         
      ! get the porosity theta value
      call get_option(trim(s_field%option_path)//'/diagnostic/algorithm/temporal_discretisation/theta', &
                      porosity_theta_value, &
                      default = 0.0)
         
      porosity_new => extract_scalar_field(state, trim(porosity_name), stat = stat)                  
  
      if (stat /=0) then 
         FLExit('Failed to extract Porosity from state from state for calculating interstitial_velocity_cv_courant_number')
      end if
         
      porosity_old => extract_scalar_field(state, "Old"//trim(porosity_name), stat = stat)

      if (stat /=0) then 
         FLExit('Failed to extract OldPorosity from state from state for calculating interstitial_velocity_cv_courant_number')
      end if
         
      call allocate(porosity_theta, porosity_new%mesh)
         
      call set(porosity_theta, porosity_new, porosity_old, porosity_theta_value)
         
      ewrite_minmax(porosity_theta)
       
      allocate(cvmass)  
      call allocate(cvmass, s_field%mesh, name="LocalCVMassWithPorosity")
      call compute_cv_mass(x, cvmass, porosity_theta)
       
      call deallocate(porosity_theta)

      ewrite_minmax(cvmass)    
      
      if(s_field%mesh%shape%degree /= 0) then

        call get_option("/geometry/quadrature/controlvolume_surface_degree", &
                      quaddegree, default=1)

        cvfaces=find_cv_faces(vertices=ele_vertices(s_field,1), &
                              dimension=mesh_dim(s_field), &
                              polydegree=s_field%mesh%shape%degree, &
                              quaddegree=quaddegree)

        u_cvshape=make_cv_element_shape(cvfaces, u%mesh%shape)
        x_cvshape=make_cv_element_shape(cvfaces, x%mesh%shape)

        allocate(x_ele(x%dim,ele_loc(x,1)), &
                x_f(x%dim, x_cvshape%ngi), &
                u_f(u%dim, u_cvshape%ngi), &
                detwei(x_cvshape%ngi), &
                normal(x%dim, x_cvshape%ngi), &
                normgi(x%dim))
        allocate(notvisited(x_cvshape%ngi))

        do ele=1, element_count(s_field)
          x_ele=ele_val(x, ele)
          x_f=ele_val_at_quad(x, ele, x_cvshape)
          u_f=ele_val_at_quad(u, ele, u_cvshape)
          nodes=>ele_nodes(s_field, ele)
          x_nodes=>ele_nodes(x_courant, ele)

          call transform_cvsurf_to_physical(x_ele, x_cvshape, &
                                            detwei, normal, cvfaces)

          notvisited=.true.

          do iloc = 1, s_field%mesh%shape%loc

            do face = 1, cvfaces%faces

              if(cvfaces%neiloc(iloc, face) /= 0) then
                oloc = cvfaces%neiloc(iloc, face)

                do gi = 1, cvfaces%shape%ngi

                  ggi = (face-1)*cvfaces%shape%ngi + gi

                  ! have we been here before?
                  if(notvisited(ggi)) then
                    notvisited(ggi)=.false.

                    normgi=orientate_cvsurf_normgi(node_val(x_courant, x_nodes(iloc)),x_f(:,ggi),normal(:,ggi))

                    udotn=dot_product(u_f(:,ggi), normgi)

                    if(udotn>0.0) then
                      income=0.0
                    else
                      income=1.0
                    end if

                    call addto(s_field, nodes(iloc), abs(udotn)*(1.-income)*detwei(ggi))
                    call addto(s_field, nodes(oloc), abs(udotn)*income*detwei(ggi)) ! notvisited

                  end if ! notvisited

                end do

              end if
            end do
          end do
        end do

        u_cvbdyshape=make_cvbdy_element_shape(cvfaces, u%mesh%faces%shape)
        x_cvbdyshape=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape)

        allocate(x_ele_bdy(x%dim,face_loc(x,1)), &
                u_bdy_f(u%dim, u_cvbdyshape%ngi), &
                detwei_bdy(x_cvbdyshape%ngi), &
                normal_bdy(x%dim, x_cvbdyshape%ngi))
        allocate(nodes_bdy(face_loc(s_field, 1)))
        allocate(courant_bc_type(surface_element_count(s_field)))

        ! get the fields over the surface containing the bcs
        call get_entire_boundary_condition(s_field, (/"internal"/), courant_bc, courant_bc_type)
        
        do sele=1,surface_element_count(s_field)
        
          if(courant_bc_type(sele)==1) cycle

          ele = face_ele(x, sele)
          x_ele = ele_val(x, ele)
          x_ele_bdy = face_val(x, sele)
          nodes_bdy=face_global_nodes(s_field, sele)

          call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, &
                                x_cvbdyshape, normal_bdy, detwei_bdy)

          u_bdy_f=face_val_at_quad(u, sele, u_cvbdyshape)

          do iloc = 1, s_field%mesh%faces%shape%loc

            do face = 1, cvfaces%sfaces

              if(cvfaces%sneiloc(iloc,face)/=0) then

                do gi = 1, cvfaces%shape%ngi

                  ggi = (face-1)*cvfaces%shape%ngi + gi
                  
                    udotn=dot_product(u_bdy_f(:,ggi), normal_bdy(:,ggi))

                    if(udotn>0.0) then
                      income=0.0
                    else
                      income=1.0
                    end if

                    call addto(s_field, nodes_bdy(iloc), abs(udotn)*(1.0-income)*detwei_bdy(ggi))

                end do

              end if

            end do

          end do

        end do

        deallocate(x_ele, x_f, u_f, detwei, normal, normgi)
        deallocate(x_ele_bdy, u_bdy_f, detwei_bdy, normal_bdy)
        deallocate(nodes_bdy)
        deallocate(notvisited)
        call deallocate(x_cvbdyshape)
        call deallocate(u_cvbdyshape)
        call deallocate(x_cvshape)
        call deallocate(u_cvshape)
        call deallocate(cvfaces)
        call deallocate(courant_bc)

      else

        allocate(detwei(face_ngi(s_field, 1)), &
                 u_f(u%dim, face_ngi(u, 1)), &
                 normal(x%dim, face_ngi(s_field, 1)))

        do ele = 1, element_count(s_field)

          nodes=>ele_nodes(s_field, ele)
          assert(size(nodes)==1)

          neigh=>ele_neigh(s_field, ele)

          do ni= 1, size(neigh)

            face = ele_face(s_field, ele, neigh(ni))

            if(neigh(ni)>0) then
              ! internal face
              face_2=ele_face(s_field, neigh(ni), ele)
            else
              ! external face
              face_2 = face
            end if

            call transform_facet_to_physical(x, face, detwei_f=detwei, normal=normal)

            ! if velocity is dg then use a trapezoidal rule (otherwise this will
            ! all cancel out to give the face value)
            u_f = 0.5*(face_val_at_quad(u, face) + face_val_at_quad(u, face_2))

            call addto(s_field, nodes(1), &
                 sum(sum(u_f*normal,1)*merge(1.0,0.0,.not.(sum(u_f*normal,1)<0.0))*detwei))

          end do


        end do

      end if

      s_field%val = s_field%val*dt/cvmass%val

      call deallocate(x_courant)
      
      call deallocate(cvmass)
      deallocate(cvmass)
      
      call halo_update(s_field)
            
      ewrite(1,*) 'Exiting calculate_interstitial_velocity_cv_courant_number'

   end subroutine calculate_interstitial_velocity_cv_courant_number

! ----------------------------------------------------------------------------------

end module porous_media_diagnostics

