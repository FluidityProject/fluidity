#include "fdebug.h"

module dgtools

use elements
use sparse_tools
use sparsity_patterns
use vector_tools
use fields
use fields_data_types
use FETools
use transform_elements
use boundary_conditions, only: get_boundary_condition_nodes

implicit none

private

public :: local_node_map, get_dg_inverse_mass_matrix, get_lumped_mass,&
     & dg_add_mass, dg_apply_mass, construct_inverse_mass_matrix_dg

!! Parameters choosing how dirichlet boundary conditions are set
integer, public, parameter :: DIRICHLET_NONE=0, &
     DIRICHLET_ONES_ON_DIAGONAL=1, &
     & DIRICHLET_BIG_SPRING=2, & 
     & DIRICHLET_WEAK=3

interface get_dg_inverse_mass_matrix
   module procedure csr_get_dg_inverse_mass_matrix, &
        dcsr_get_dg_inverse_mass_matrix, &
        csr_dg_inverse_mass_from_mass
end interface

interface get_lumped_mass
   module procedure dcsr_get_lumped_mass
end interface

interface dg_apply_mass
   module procedure csr_dg_apply_mass
end interface
   
interface dg_add_mass
   module procedure csr_dg_add_mass
end interface

interface construct_inverse_mass_matrix_dg
    module procedure construct_inverse_mass_matrix_dg_scalar, &
          construct_inverse_mass_matrix_dg_vector
end interface

contains
  
  function local_node_map(m, m_f, bdy, bdy_2) result(local_glno)
    ! Fill in the number map for the DG double element.
    type(element_type), intent(in) :: m, m_f
    integer, dimension(m_f%loc) :: bdy, bdy_2
    integer, dimension(m%loc,2) :: local_glno

    integer :: i,j

    local_glno=0

    ! First m_f%loc places are for the bdy between the elements.
    forall(i=1:m_f%loc)
       local_glno(bdy(i),1)=i
    end forall

    ! Remaining spots go to elements. 
    j=m_f%loc
    do i=1, m%loc
       if(local_glno(i,1)==0) then
          j=j+1
          local_glno(i,1)=j
       end if
    end do

    ASSERT(j==m%loc)

    ! First m_f%loc places are for the bdy between the elements.
    forall(i=1:m_f%loc)
       local_glno(bdy_2(i),2)=i
    end forall

    ! Remaining spots go to elements. 
    j=m%loc
    do i=1, m%loc
       if(local_glno(i,2)==0) then
          j=j+1
          local_glno(i,2)=j
       end if
    end do

    ASSERT(j==2*m%loc-m_f%loc)

  end function local_node_map

  subroutine csr_get_dg_inverse_mass_matrix(inverse_mass, dg_mesh, &
       & positions, density, dirichlet_list, dirichlet_flag, &
       & absorption_factor, allocate_matrix)
    type(csr_matrix), intent(inout) :: inverse_mass
    type(mesh_type), intent(inout) :: dg_mesh
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in), optional, target :: density
    integer, dimension(:), intent(in), optional :: dirichlet_list
    integer, intent(in), optional :: dirichlet_flag
    type(scalar_field), intent(in), optional :: absorption_factor
    logical, intent(in), optional :: allocate_matrix

    !locals
    integer :: ele
    real, dimension(:), allocatable :: factor_at_quad
    logical, dimension(:), allocatable :: internal_dirichlet_list
    integer, dimension(:), pointer :: e_nodes
    integer :: l_dirichlet_flag
    type(csr_sparsity) :: sparsity

    if(present(dirichlet_flag)) then
       l_dirichlet_flag = dirichlet_flag
    else
       l_dirichlet_flag = DIRICHLET_NONE
    end if
    
    if(l_dirichlet_flag.ne.DIRICHLET_NONE) then
       if(present(dirichlet_list)) then
          allocate( internal_dirichlet_list( node_count(dg_mesh) ) )
          internal_dirichlet_list = .false.
          internal_dirichlet_list(dirichlet_list) = .true.
       end if
    end if
    
    if(l_dirichlet_flag==DIRICHLET_BIG_SPRING) then
      ewrite(2,*) 'INFINITY = ', INFINITY
      ewrite(2,*) 'sqrt(INFINITY) = ', sqrt(INFINITY)
    end if

    if (.not. present_and_false(allocate_matrix)) then
      assert(dg_mesh%continuity==-1)

      sparsity=make_sparsity_dg_mass(dg_mesh)

      call allocate(inverse_mass, sparsity, name="DGInverseMass")
      
      ! Drop the extra reference to sparsity.
      call deallocate(sparsity)
    end if
    
    allocate( factor_at_quad(1:ele_ngi(positions,1)) )

    do ele = 1, dg_mesh%elements

       if (present(density)) then
          factor_at_quad=ele_val_at_quad(density, ele)
       else
          factor_at_quad=1.0
       end if
       if (present(absorption_factor)) then
          factor_at_quad=factor_at_quad* &
             ele_val_at_quad(absorption_factor, ele)
       end if
       
       if(present(dirichlet_list).and. &
            & (l_dirichlet_flag.ne.DIRICHLET_NONE)) then
          e_nodes => ele_nodes(dg_mesh,ele)
          call csr_assemble_local_dg_inverse_mass_matrix(inverse_mass, &
               dg_mesh,positions, factor_at_quad, ele,&
               internal_dirichlet_list(e_nodes), l_dirichlet_flag)
       else
          call csr_assemble_local_dg_inverse_mass_matrix(inverse_mass, &
               dg_mesh,positions, factor_at_quad, ele)
       end if
       
    end do

  end subroutine csr_get_dg_inverse_mass_matrix

  subroutine csr_assemble_local_dg_inverse_mass_matrix(inverse_mass, &
       dg_mesh,positions, factor_at_quad, ele, &
       dirichlet_list, dirichlet_flag)
    type(csr_matrix), intent(inout) :: inverse_mass
    type(mesh_type), intent(in) :: dg_mesh
    type(vector_field), intent(in) :: positions
    ! a scalar factor to be integrated with the mass matrix:
    !      /
    ! M_ij=| factor(x) N_i(x) N_j(x) dx
    !      /
    real, dimension(:), intent(in) :: factor_at_quad
    integer, intent(in) :: ele
    logical, dimension(:), intent(in), optional :: dirichlet_list
    integer, intent(in), optional :: dirichlet_flag

    !local variables
    real, dimension(ele_loc(dg_mesh,ele),ele_loc(dg_mesh,ele)) :: local_mass
    real, dimension(dg_mesh%shape%ngi) :: detwei
    integer, dimension(:), pointer :: ele_dg
    type(element_type), pointer :: shape_dg,shape_X
    integer :: i

    !assemble local mass matrix
    ele_dg=>ele_nodes(dg_mesh, ele)
    shape_dg=>ele_shape(dg_mesh, ele)
    shape_X=>ele_shape(positions, ele)

    call transform_to_physical(positions,ele,detwei=detwei)
    detwei=detwei*factor_at_quad
    local_mass = shape_shape(shape_dg,shape_dg,detwei)
    
    if(present(dirichlet_list)) then
       do i = 1, size(dirichlet_list)
          if(dirichlet_list(i)) then
             select case(dirichlet_flag)
             case (DIRICHLET_NONE)
             case (DIRICHLET_ONES_ON_DIAGONAL)
                local_mass(:,i) = 0.
                local_mass(i,:) = 0.
                local_mass(i,i) = 1.
             case (DIRICHLET_BIG_SPRING)
                local_mass(i,i) = sqrt(INFINITY)
             case default
                FLAbort('bad dirichlet flag')
             end select
          end if
       end do
    end if

    call invert(local_mass)
    
    call set(inverse_mass,ele_dg,ele_dg,local_mass)

  end subroutine csr_assemble_local_dg_inverse_mass_matrix

  subroutine dcsr_get_dg_inverse_mass_matrix(inverse_mass,dg_mesh, &
       & positions, density,dirichlet_list,dirichlet_flag, &
       & absorption_factor)
    type(mesh_type), intent(inout) :: dg_mesh
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in), optional, target :: density
    type(dynamic_csr_matrix), intent(inout) :: inverse_mass
    integer, dimension(:), intent(in), optional :: dirichlet_list
    integer, intent(in), optional :: dirichlet_flag
    type(scalar_field), intent(in), optional :: absorption_factor

    !locals
    integer :: ele
    real, dimension(:), allocatable :: factor_at_quad
    logical, dimension(:), allocatable :: internal_dirichlet_list
    integer, dimension(:), pointer :: e_nodes
    integer :: l_dirichlet_flag

    if(present(dirichlet_flag)) then
       l_dirichlet_flag = dirichlet_flag
    else
       l_dirichlet_flag = 0
    end if

    if(l_dirichlet_flag.ne.DIRICHLET_NONE) then
       if(present(dirichlet_list)) then
          allocate( internal_dirichlet_list( node_count(dg_mesh) ) )
          internal_dirichlet_list = .false.
          internal_dirichlet_list(dirichlet_list) = .true.
       end if
    end if

    assert(dg_mesh%continuity==-1)

    allocate( factor_at_quad(1:ele_ngi(density,1)) )
    
    do ele = 1, dg_mesh%elements

       if (present(density)) then
          factor_at_quad=ele_val_at_quad(density, ele)
       else
          factor_at_quad=1.0
       end if
       if (present(absorption_factor)) then
          factor_at_quad=factor_at_quad* &
             ele_val_at_quad(absorption_factor, ele)
       end if
       
       if(present(dirichlet_list).and. &
            & (l_dirichlet_flag.ne.DIRICHLET_NONE)) then
          e_nodes => ele_nodes(dg_mesh,ele)
          call assemble_local_dg_inverse_mass_matrix(inverse_mass, &
               dg_mesh,positions,factor_at_quad,ele,&
               internal_dirichlet_list(e_nodes), l_dirichlet_flag)
       else
          call assemble_local_dg_inverse_mass_matrix(inverse_mass, &
               dg_mesh,positions,factor_at_quad,ele)
       end if
    end do

  end subroutine dcsr_get_dg_inverse_mass_matrix

  subroutine assemble_local_dg_inverse_mass_matrix(inverse_dynamic_mass, &
       dg_mesh,positions,factor_at_quad,ele,dirichlet_list,dirichlet_flag)
    type(dynamic_csr_matrix), intent(inout) :: inverse_dynamic_mass
    type(mesh_type), intent(in) :: dg_mesh
    type(vector_field), intent(in) :: positions
    ! a scalar factor to be integrated with the mass matrix:
    !      /
    ! M_ij=| factor(x) N_i(x) N_j(x) dx
    !      /
    real, dimension(:), intent(in) :: factor_at_quad
    integer, intent(in) :: ele
    logical, dimension(:), intent(in), optional :: dirichlet_list
    integer, intent(in), optional :: dirichlet_flag

    !local variables
    real, dimension(ele_loc(dg_mesh,ele),ele_loc(dg_mesh,ele)) :: local_mass
    real, dimension(dg_mesh%shape%ngi) :: detwei
    integer, dimension(:), pointer :: ele_dg
    type(element_type), pointer :: shape_dg,shape_X
    integer :: i

    !assemble local mass matrix
    ele_dg=>ele_nodes(dg_mesh, ele)
    shape_dg=>ele_shape(dg_mesh, ele)
    shape_X=>ele_shape(positions, ele)

    call transform_to_physical(positions, ele,detwei=detwei)
    local_mass = shape_shape(shape_dg,shape_dg,detwei*factor_at_quad)

    if(present(dirichlet_list)) then
       do i = 1, size(dirichlet_list)
          if(dirichlet_list(i)) then
             select case(dirichlet_flag)
             case (DIRICHLET_NONE)
             case (DIRICHLET_ONES_ON_DIAGONAL)
                local_mass(:,i) = 0.
                local_mass(i,:) = 0.
                local_mass(i,i) = 1.
             case (DIRICHLET_BIG_SPRING)
                local_mass(i,i) = sqrt(INFINITY)
             case default
                FLAbort('bad dirichlet flag')
             end select
          end if
       end do
    end if

    call invert(local_mass)

    call set(inverse_dynamic_mass,ele_dg,ele_dg,local_mass)

  end subroutine assemble_local_dg_inverse_mass_matrix

  subroutine dcsr_get_lumped_mass(mass,mesh, &
       & positions,dirichlet_list,dirichlet_flag)
    type(mesh_type), intent(in) :: mesh
    type(vector_field), intent(in) :: positions
    type(dynamic_csr_matrix), intent(inout) :: mass
    integer, dimension(:), intent(in), optional :: dirichlet_list
    integer, intent(in), optional :: dirichlet_flag

    !locals
    integer :: ele
    logical, dimension(:), allocatable :: internal_dirichlet_list
    integer, dimension(:), pointer :: e_nodes
    integer :: l_dirichlet_flag

    if(present(dirichlet_flag)) then
       l_dirichlet_flag = dirichlet_flag
    else
       l_dirichlet_flag = 0
    end if
    
    if(present(dirichlet_list).and.(l_dirichlet_flag.ne.0)) then
       allocate( internal_dirichlet_list( node_count(mesh) ) )
       internal_dirichlet_list = .false.
       internal_dirichlet_list(dirichlet_list) = .true.
    end if

    do ele = 1, mesh%elements

       if(present(dirichlet_list).and.(dirichlet_flag.ne.0)) then
          e_nodes => ele_nodes(mesh,ele)
          call dcsr_assemble_local_lumped_mass(mass, &
               mesh,positions,ele,internal_dirichlet_list(e_nodes), &
               l_dirichlet_flag)
       else
          call dcsr_assemble_local_lumped_mass(mass, &
               mesh,positions,ele)
       end if
    end do

  end subroutine dcsr_get_lumped_mass

  subroutine dcsr_assemble_local_lumped_mass(mass, &
       mesh,positions,ele,dirichlet_list,dirichlet_flag)
    integer, intent(in) :: ele
    type(dynamic_csr_matrix), intent(inout) :: mass
    type(mesh_type), intent(in) :: mesh
    type(vector_field), intent(in) :: positions
    logical, dimension(:), intent(in), optional :: dirichlet_list
    integer, intent(in), optional :: dirichlet_flag

    !local variables
    real, dimension(ele_loc(mesh,ele),ele_loc(mesh,ele)) :: local_mass
    real, dimension(mesh%shape%ngi) :: detwei
    integer, dimension(:), pointer :: ele_dg
    type(element_type), pointer :: shape_dg,shape_X
    integer :: i

    !assemble local mass matrix
    ele_dg=>ele_nodes(mesh, ele)
    shape_dg=>ele_shape(mesh, ele)
    shape_X=>ele_shape(positions, ele)

    call transform_to_physical(positions, ele, detwei=detwei)
    local_mass = shape_shape(shape_dg,shape_dg,detwei)

    if(present(dirichlet_list)) then
       do i = 1, size(dirichlet_list)
          if(dirichlet_list(i)) then
             select case(dirichlet_flag)
             case (DIRICHLET_NONE)
             case (DIRICHLET_ONES_ON_DIAGONAL)
                local_mass(:,i) = 0.
                local_mass(i,:) = 0.
                local_mass(i,i) = 1.
             case (DIRICHLET_BIG_SPRING)
                local_mass(i,i) = sqrt(INFINITY)
             case default
                FLAbort('bad dirichlet flag')
             end select
          end if
       end do
    end if

    do i = 1, ele_loc(mesh,ele)
       call set(mass,ele_dg(i),ele_dg(i),sum(local_mass(i,:)))
    end do

  end subroutine dcsr_assemble_local_lumped_mass

  subroutine csr_dg_inverse_mass_from_mass(inv_mass, mass)
    !!< Put the inverse of mass into inv_mass. This is short-circuited by
    !!< knowing that mass is DG.
    type(csr_matrix), intent(inout) :: inv_mass
    type(csr_matrix), intent(in) :: mass

    integer :: row, colm_pos, nloc

    row=0
    colm_pos=0
    
    do 
       if(row>=size(mass,1)) exit
       nloc=row_length(mass, row+1)
       inv_mass%val(colm_pos+1:colm_pos+nloc**2) &
            = reshape(&
            &  inverse(&
            &     reshape(mass%val(colm_pos+1:colm_pos+nloc**2), &
            &             (/nloc,nloc/))&
            &           ), &
            &   (/nloc*nloc/))

       row=row+nloc
       colm_pos=colm_pos+nloc**2
    end do
 
  end subroutine csr_dg_inverse_mass_from_mass

  subroutine csr_dg_apply_mass(mass, field)
    !!< return field=mass*field. This is basically an optimised in-place
    !!< matrix multiply.
    type(csr_matrix), intent(in) :: mass
    type(scalar_field), intent(inout) :: field

    integer :: row, colm_pos, nloc

    row=0
    colm_pos=0
    
    do
       if(row>=size(mass,1)) exit
       nloc=row_length(mass, row+1)
       
       field%val(row+1:row+nloc) = &
            matmul(reshape(mass%val(colm_pos+1:colm_pos+nloc**2), &
            &             (/nloc,nloc/)), &
            &      field%val(row+1:row+nloc))

       row=row+nloc
       colm_pos=colm_pos+nloc**2
    end do
    
  end subroutine csr_dg_apply_mass

  subroutine csr_dg_add_mass(matrix, mass)
    !!< Add mass to matrix. This is an optimised addto operation.
    type(csr_matrix), intent(inout) :: matrix
    type(csr_matrix), intent(in) :: mass

    integer :: row, colm_pos, nloc
    
    row=0
    colm_pos=0
    
    do 
       if(row>=size(mass,1)) exit
       nloc=row_length(mass, row+1)
       
       call addto(matrix, row_m_ptr(mass, row+1), row_m_ptr(mass, row+1), &
            reshape(mass%val(colm_pos+1:colm_pos+nloc**2), &
            &             (/nloc,nloc/)))

       row=row+nloc
       colm_pos=colm_pos+nloc**2
    end do
    
  end subroutine csr_dg_add_mass

  subroutine construct_inverse_mass_matrix_dg_scalar(inverse_mass, sfield, x)
    !! This constructs the inverse mass matrix for a scalar field,
    !! using options and bcs attached to the field.
    type(csr_matrix), intent(out):: inverse_mass
    type(scalar_field), intent(inout) :: sfield
    type(vector_field), intent(inout) :: x
      
    integer, allocatable, dimension(:) :: sfield_bc_type
    integer, allocatable, dimension(:) :: dirichlet_list
    logical :: dirichlet
    integer i, count
    
    ewrite(1,*) 'Construct the DG inverse mass matrix for a scalar field'
    
    allocate(sfield_bc_type(node_count(sfield)))
    call get_boundary_condition_nodes(sfield, (/"dirichlet"/), sfield_bc_type)
    if(any(sfield_bc_type==1)) then
      dirichlet = .true.
    else
      dirichlet = .false.
    end if
    
    if (dirichlet) then
      
      allocate(dirichlet_list(sum(sfield_bc_type)))
      count = 0
      do i = 1, size(sfield_bc_type)
        if(sfield_bc_type(i)==1) then
            count = count + 1
            dirichlet_list(count) = i
        end if
      end do
        
      call get_dg_inverse_mass_matrix(inverse_mass,sfield%mesh, &
              x, &
              dirichlet_list=dirichlet_list, &
              dirichlet_flag=DIRICHLET_BIG_SPRING)
      
      deallocate(dirichlet_list)
    
    else
    
      ! the mass matrix has no dirichlet modifications
      call get_dg_inverse_mass_matrix(inverse_mass,sfield%mesh, &
                                        x)
      
    end if
    
  end subroutine construct_inverse_mass_matrix_dg_scalar

  subroutine construct_inverse_mass_matrix_dg_vector(inverse_mass, vfield, x)
    !! This constructs the inverse mass matrix for all components
    !! of a vector field, using options and bcs attached to the
    !! field. This version is a bit slow (computes transform_to_physical
    !! and inverses u%dim times, when diagonal mass blocks are different)
    type(block_csr_matrix), intent(out):: inverse_mass
    type(vector_field), intent(inout) :: vfield, x
      
    type(csr_sparsity):: sparsity
    type(csr_matrix):: inverse_mass_block
    integer, allocatable, dimension(:,:) :: vfield_bc_type
    integer, allocatable, dimension(:) :: dirichlet_list
    logical :: dirichlet
    integer i, dim, count
    
    ewrite(1,*) 'Construct the DG inverse mass matrix for a vector field'
    
    allocate(vfield_bc_type(vfield%dim, node_count(vfield)))
    call get_boundary_condition_nodes(vfield, (/"dirichlet"/), vfield_bc_type)
    if(any(vfield_bc_type==1)) then
      dirichlet = .true.
    else
      dirichlet = .false.
    end if
    
    if (dirichlet) then
      
      assert(vfield%mesh%continuity==-1)
      sparsity=make_sparsity_dg_mass(vfield%mesh)
      call allocate(inverse_mass, sparsity, (/ vfield%dim, vfield%dim /), &
        name="DGInverseMass", diagonal=.true.)
      ! Drop the extra reference to sparsity.
      call deallocate(sparsity)
      
      do dim=1, vfield%dim
        
         allocate(dirichlet_list(sum(vfield_bc_type(dim,:))))
         count = 0
         do i = 1, size(vfield_bc_type(dim,:))
            if(vfield_bc_type(dim,i)==1) then
                count = count + 1
                dirichlet_list(count) = i
            end if
         end do
           
          inverse_mass_block=block(inverse_mass, dim, dim)
          call get_dg_inverse_mass_matrix(inverse_mass_block,vfield%mesh, &
                  x, &
                  dirichlet_list=dirichlet_list, &
                  dirichlet_flag=DIRICHLET_BIG_SPRING, &
                  allocate_matrix=.false.)
          
          deallocate(dirichlet_list)
      
      end do
    
    else
    
      ! the mass matrix is the same for all components
      assert(vfield%mesh%continuity==-1)
      sparsity=make_sparsity_dg_mass(vfield%mesh)
      call allocate(inverse_mass, sparsity, (/ vfield%dim, vfield%dim /), &
        name="DGInverseMass", diagonal=.true., equal_diagonal_blocks=.true.)
      ! Drop the extra reference to sparsity.
      call deallocate(sparsity)
      
      inverse_mass_block=block(inverse_mass, 1,1)
      call get_dg_inverse_mass_matrix(inverse_mass_block,vfield%mesh, &
                                        x)
      
    end if
    
  end subroutine construct_inverse_mass_matrix_dg_vector

end module dgtools
