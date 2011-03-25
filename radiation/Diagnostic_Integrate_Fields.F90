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

module diagnostic_integrate_fields

   !!< This module contains procedures associated with diagnostic integration of fields
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   use parallel_tools
   use fields

   implicit none
   
   private 

   public :: integrate
   
   interface integrate
      module procedure volume_integrate_product_scalar_fields
   end interface integrate
   
contains

   ! --------------------------------------------------------------------------
   
   subroutine volume_integrate_product_scalar_fields(scalar_fields, &
                                                     positions, &
                                                     integral, &
                                                     region_ids)
      
      !!< Calculate the integral of the product of a set of scalar fields that share a positions mesh
      !!< If region_ids is present then only the corresponding positions field regions are integrated
      !!< else the whole domain is included

      type(scalar_field_pointer), dimension(:), intent(in)  :: scalar_fields
      type(vector_field),                       intent(in)  :: positions
      real,                                     intent(out) :: integral
      integer,                    dimension(:), intent(in), optional :: region_ids
      
      ! local variables
      logical :: found_id
      integer :: vele
      integer :: s
      integer :: id
      integer :: vele_id
      real, dimension(:), allocatable :: detwei_vele
      real, dimension(:), allocatable :: product_ele_val_at_quad

      ewrite(1,*) 'Integrate the product of a set of scalar fields over spatial domain'
      
      ! assert that each scalar field has the same number of elements and the same dim as positions
      field_assert_vele_loop: do s = 1,size(scalar_fields)

         assert(ele_count(positions) == ele_count(scalar_fields(s)%ptr))

         assert(mesh_dim(positions) == mesh_dim(scalar_fields(s)%ptr))
      
      end do field_assert_vele_loop
      
      ! if region_ids is present assert that it has something
      if (present(region_ids)) then 
         
         assert(size(region_ids) > 0)
      
      end if   
          
      ! initialise the integral as summed over elements
      integral = 0.0
      
      ! loop the volume elements to perform the integration
      velement_loop: do vele = 1,ele_count(scalar_fields(1)%ptr)
         
         ! if present only conisder input region_ids
         region_id_present: if (present(region_ids)) then
            
            ! initialise flag for whether this volume element vele should be considered
            found_id = .false.
            
            ! find the positions field vele region id
            vele_id = ele_region_id(positions,vele)
            
            region_id_loop: do id = 1,size(region_ids)
               
               check_id: if (vele_id == region_ids(id)) then
                  
                  found_id = .true.
                  
                  exit region_id_loop
               
               end if check_id
               
            end do region_id_loop
            
            ! if not found an id match then cycle the volume element loop
            if (.not. found_id) cycle velement_loop
            
         end if region_id_present
         
         ! allocate the jacobian transform and gauss weight array for this vele
         allocate(detwei_vele(ele_ngi(scalar_fields(1)%ptr,vele)))
                 
         ! form the velement jacobian transform and gauss weight
         call transform_to_physical(positions, vele, detwei = detwei_vele)
         
         ! form the product of the scalar fields at the quadrature points
         
         allocate(product_ele_val_at_quad(ele_ngi(scalar_fields(1)%ptr,vele)))         
                           
         field_val_loop: do s = 1,size(scalar_fields)
            
            assert(ele_ngi(positions,vele) == ele_ngi(scalar_fields(s)%ptr,vele))
            
            first: if (s == 1) then
            
               product_ele_val_at_quad =  ele_val_at_quad(scalar_fields(s)%ptr, vele)
            
            else first
            
               product_ele_val_at_quad = product_ele_val_at_quad * &
                                         ele_val_at_quad(scalar_fields(s)%ptr, vele)
            
            end if first
            
         end do field_val_loop
                  
         ! form the integral for this volume element and add to the total
         integral = integral + dot_product(detwei_vele,product_ele_val_at_quad)
         
         deallocate(detwei_vele)
         deallocate(product_ele_val_at_quad)
                                          
      end do velement_loop

      ! sum the value of the processes
      call allsum(integral)

   end subroutine volume_integrate_product_scalar_fields
   
   ! --------------------------------------------------------------------------

end module diagnostic_integrate_fields
