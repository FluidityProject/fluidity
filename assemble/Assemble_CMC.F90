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

module assemble_CMC

  use fldebug
  use spud
  use global_parameters, only: OPTION_PATH_LEN
  use sparse_tools
  use linked_lists
  use transform_elements
  use fetools, only: shape_shape
  use elements
  use fields
  use sparse_tools_petsc
  use state_module
  use sparse_matrices_fields
  use field_options
  use fefields
  
  implicit none 

  private
  
  public :: assemble_cmc_dg, assemble_masslumped_cmc, &
            assemble_masslumped_ctm, repair_stiff_nodes, zero_stiff_nodes, &
            assemble_diagonal_schur, assemble_scaled_pressure_mass_matrix

contains

    subroutine assemble_cmc_dg(CMC, CTP, CT, inverse_mass)
      !!< Assemble the pressure matrix C^T M^{-1} C for a DG mesh.
      !!< This currently does not support rotations.
      type(csr_matrix), intent(inout) :: CMC
      type(block_csr_matrix), intent(in) :: CTP, CT
      type(block_csr_matrix), intent(in):: inverse_mass

      type(csr_matrix) :: CM ! Temporary half way matrix
      type(csr_matrix) :: CMC_tmp ! Temporary accumulator.
      type(csr_matrix) :: CT_block ! The current CT dimension
      type(csr_matrix) :: CTP_block ! The current CTP dimension

      integer :: dim, dim2

      call zero(CMC)
      
      if(inverse_mass%diagonal) then

         do dim=1, blocks(CT,2)
            
            ! This is a borrowed reference so no need to deallocate.
            CT_block=block(CT, 1, dim)
            CTP_block=block(CTP, 1, dim)
            
            CM=matmul_T(CTP_block, block(inverse_mass,dim,dim), CT%sparsity)
            
            CMC_tmp=matmul_T(CM, CT_block, model=CMC%sparsity)
            
            call addto(CMC, CMC_tmp)
            
            call deallocate(CM)
            call deallocate(CMC_tmp)

         end do

      else

         do dim=1, blocks(CT,2)
            do dim2 = 1, blocks(CT,2)
               
               ! This is a borrowed reference so no need to deallocate.
               CT_block=block(CT, 1, dim)
               CTP_block=block(CTP, 1, dim2)
               
               CM=matmul_T(CTP_block, block(inverse_mass,dim2,dim), CT%sparsity)
               
               CMC_tmp=matmul_T(CM, CT_block, model=CMC%sparsity)
               
               call addto(CMC, CMC_tmp)
               
               call deallocate(CM)
               call deallocate(CMC_tmp)
            end do
         end do

      end if

      ewrite_minmax(cmc)

    end subroutine assemble_cmc_dg

    subroutine assemble_masslumped_cmc(cmc_m, ctp_m, inverse_masslump, ct_m)
      !!< Assemble the pressure matrix C_P^T M_l^{-1} C.
      !!< This currently does not support rotations.

      ! this subroutine is designed to start to replace cmc_wrapper and getcmc

      type(csr_matrix), intent(inout) :: cmc_m
      type(block_csr_matrix), intent(in) :: ctp_m, ct_m
      type(vector_field), intent(in) :: inverse_masslump

      ewrite(1,*) 'Entering assemble_masslumped_cmc'

      call mult_div_vector_div_T(cmc_m, ctp_m, inverse_masslump, ct_m)
      
      ewrite_minmax(cmc_m)

    end subroutine assemble_masslumped_cmc

    subroutine assemble_diagonal_schur(schur_diagonal_matrix,u,inner_m,ctp_m,ct_m)
      !!< Assemble the matrix C_P^T * [(Big_m)_diagonal]^-1 * C. 
      !!< This is used as a preconditioner for the full projection solve
      !!< when using the full momentum matrix.

      ! Fluidity velocity vector:
      type(vector_field), intent(in):: u
      ! Divergence matrices:
      type(block_csr_matrix), intent(in) :: ctp_m, ct_m
      ! Inner matrix of Schur complement:
      type(petsc_csr_matrix), intent(inout) :: inner_m   
      ! Product matrix:
      type(csr_matrix), intent(inout):: schur_diagonal_matrix
      ! Diagonal of inner_m - required to set up preconditioner matrix for Stokes problems
      ! (i.e. DiagonalSchurComplement):
      type(vector_field) :: inner_m_diagonal   

      integer :: i

      ewrite(1,*) 'Entering assemble_diagonal_schur'

      call zero(schur_diagonal_matrix)
      call allocate(inner_m_diagonal,u%dim,u%mesh,"Diagonal_inner_m")
      call zero(inner_m_diagonal)
      call extract_diagonal(inner_m,inner_m_diagonal)

      ewrite_minmax(inner_m_diagonal)
      if(any(inner_m_diagonal%val < 0)) then
        ewrite(-1,*) 'Inner_m_diagonal has negative values'
        FLExit("Negative values in the diagonal schur complement preconditioner")

      end if

      call mult_div_invvector_div_T(schur_diagonal_matrix, ctp_m, inner_m_diagonal, ct_m)
      ewrite_minmax(schur_diagonal_matrix)
      call deallocate(inner_m_diagonal)

    end subroutine assemble_diagonal_schur

    subroutine assemble_scaled_pressure_mass_matrix(state, scaled_pressure_mass_matrix, p_mesh, dt)

      ! This routine assembles the scaled_pressure_mass_matrix at the
      ! quadrature points. It is scaled by the inverse of viscosity.

      type(state_type), intent(in) :: state   
      ! Scaled pressure mass matrix - already allocated in Momentum_Eq:
      type(csr_matrix), intent(inout) :: scaled_pressure_mass_matrix

      ! Pressure mesh (for free surface this should be the extended pressure mesh)
      type(mesh_type), intent(in) :: p_mesh
      real, intent(in) :: dt

      ! Viscosity tensor:
      type(tensor_field), pointer :: viscosity     
      ! Viscosity component:
      type(scalar_field) :: viscosity_component
      ! Positions:
      type(vector_field), pointer :: positions

      ! Relevant declerations for mass matrix calculation:
      integer :: ele
      real, dimension(:), allocatable :: detwei
      type(element_type), pointer :: p_shape
      real, dimension(:,:), allocatable :: mass_matrix
      real, dimension(:), allocatable :: mu_gi

      ewrite(1,*) 'Entering assemble_scaled_pressure_mass_matrix'    

      ! Positions:
      positions => extract_vector_field(state, "Coordinate")

      ! Extract viscosity tensor from state:
      viscosity => extract_tensor_field(state,'Viscosity')

      ! Extract first component of viscosity tensor from full tensor:
      viscosity_component = extract_scalar_field(viscosity,1,1)

      ! Initialise and assemble scaled pressure mass matrix:
      allocate(detwei(ele_ngi(p_mesh, 1)), &
               mass_matrix(ele_loc(p_mesh, 1), ele_loc(p_mesh, 1)), &
               mu_gi(ele_ngi(viscosity_component, 1)))
 
      call zero(scaled_pressure_mass_matrix)

      do ele = 1, ele_count(p_mesh)
        p_shape => ele_shape(p_mesh, ele)
        mu_gi = ele_val_at_quad(viscosity_component, ele)
        call transform_to_physical(positions, ele, detwei=detwei)
        mass_matrix = shape_shape(p_shape, p_shape, detwei/(mu_gi*dt))
        call addto(scaled_pressure_mass_matrix, ele_nodes(p_mesh, ele),&
             ele_nodes(p_mesh, ele), mass_matrix)
      end do

      ewrite_minmax(scaled_pressure_mass_matrix)

      deallocate(detwei, mass_matrix, mu_gi)

    end subroutine assemble_scaled_pressure_mass_matrix
      
    subroutine repair_stiff_nodes(cmc_m, stiff_nodes_list)
    
      type(csr_matrix), intent(inout) :: cmc_m
      type(ilist), intent(inout) :: stiff_nodes_list
      
      integer :: row
      real, pointer :: row_diag
      integer, dimension(:), pointer :: row_m
      real, dimension(:), pointer :: row_val
      real :: tolerance
      
      ewrite(1,*) 'in repair_stiff_nodes()'
      
      tolerance = maxval(cmc_m%val)*epsilon(0.0)
      ewrite(2,*) 'tolerance = ', tolerance
      
      call flush_list(stiff_nodes_list)
      
      do row = 1, size(cmc_m, 1)
        row_diag=>diag_val_ptr(cmc_m, row)
        row_m=>row_m_ptr(cmc_m, row)
        row_val=>row_val_ptr(cmc_m, row)
        if(row_diag<tolerance) then
          ewrite(2,*) 'before, row, row_diag, sum(row_val) = ', row, row_diag, sum(abs(row_val))
          where(cmc_m%sparsity%colm==row) cmc_m%val = 0.0
          call zero_row(cmc_m, row)
          call addto_diag(cmc_m, row, 1.0)
          call insert(stiff_nodes_list, row)
        end if
      end do
      
      call print_list(stiff_nodes_list, 2)

    end subroutine repair_stiff_nodes
      
    subroutine zero_stiff_nodes(rhs, stiff_nodes_list)
    
      type(scalar_field), intent(inout) :: rhs
      type(ilist), intent(in) :: stiff_nodes_list
      
      real, dimension(:), pointer :: row_val
      
      ewrite(1,*) 'in zero_stiff_nodes()'
      
      if(stiff_nodes_list%length>0) then
        ewrite(2,*) 'before node_val = ', node_val(rhs, list2vector(stiff_nodes_list))
        call set(rhs, list2vector(stiff_nodes_list), spread(0.0, 1, stiff_nodes_list%length))
      end if

    end subroutine zero_stiff_nodes
    
    subroutine assemble_masslumped_ctm(ctm_m, ctp_m, masslump)
      !!< Assemble the matrix C_P^T M_l^{-1}
      !!< This currently does not support rotations.

      type(block_csr_matrix), intent(in) :: ctp_m
      type(scalar_field), intent(in) :: masslump
      type(block_csr_matrix), intent(inout) :: ctm_m

      type(csr_matrix) :: lctm_m_block

      integer :: dim, row
      real, dimension(:), pointer :: row_val
      integer, dimension(:), pointer :: row_indices

      ewrite(1,*) 'Entering assemble_masslumped_ctm'

      call zero(ctm_m)

      do dim = 1, ctm_m%blocks(2)

        lctm_m_block = block(ctm_m, 1, dim)

        do row = 1, size(ctp_m, 1)
          row_indices=>row_m_ptr(ctp_m, row)
          row_val=>row_val_ptr(ctp_m, 1, dim, row)
          call set(lctm_m_block, (/row/), row_indices, &
                  spread((row_val/node_val(masslump, row_indices)), 1, 1))
        end do

      end do

      ewrite_minmax(ctm_m)

    end subroutine assemble_masslumped_ctm
    
end module assemble_cmc
