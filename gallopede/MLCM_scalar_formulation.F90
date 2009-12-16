#include "fdebug.h"

module MLCM_scalar_formulation

  ! The routines to perform the elliptic solve for the MLCM equations
  ! in scalar form.

  use fields
!  use petsc_tools
!  use gallopede_solvers
  use global_parameters_gallopede
  use transform_elements
  use FEtools
  use DGtools
  use sparse_tools
  use solvers
  use data_structures
  use petsc_tools
  use signal_vars

implicit none

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscmg.h"

  contains

  subroutine get_mass(Mass,positions,s_field,v_field)
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout), optional :: s_field
    type(vector_field), intent(inout), optional :: v_field
    type(csr_matrix), intent(inout) :: Mass

    ! We form a mass matrix
    !local
    integer :: ele, n_ele
    
    ewrite(1,*) 'subroutine get_mass'

    if(present(s_field)) then
       n_ele=element_count(s_field)
    else
       n_ele=element_count(v_field)
    end if

    ! Assemble A element by element.
    do ele=1, n_ele

       if(present(s_field)) then
          call assemble_mass_element_contribution( &
               mass, positions, ele, s_field=s_field)
       else if(present(v_field)) then
          call assemble_mass_element_contribution( &
               mass, positions, ele, v_field=v_field)
       else
          print *, 'bad call to get_mass'
          stop
       end if

    end do

    ewrite(1,*) 'END subroutine get_mass'

  end subroutine get_mass

  subroutine assemble_mass_element_contribution( &
       mass, positions, ele, s_field,v_field)
    type(csr_matrix), intent(inout) :: mass
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in), optional :: s_field    
    type(vector_field), intent(in), optional :: v_field
    integer, intent(in) :: ele

    ! Locations of nodes.
    real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
    ! Locations of quadrature points
    real, dimension(positions%dim,ele_ngi(positions,ele)) :: X_quad
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei    
    ! Node numbers of field element.
    integer, dimension(:), pointer :: ele_field
    ! Shape functions.
    type(element_type), pointer :: shape_field, shape_X

    if(present(s_field)) then
       ele_field=>ele_nodes(s_field, ele)
       shape_field=>ele_shape(s_field, ele)
    end if
    if(present(v_field)) then
       ele_field=>ele_nodes(v_field, ele)
       shape_field=>ele_shape(v_field, ele)
    end if
    shape_X=>ele_shape(positions, ele)

    ! Locations of local vertices.
    X_ele=ele_val(positions, ele)

    ! Locations of quadrature points.
    X_quad=ele_val_at_quad(positions, ele)


    ! Transform derivatives and weights into physical space.
    call transform_to_physical(X_ele, shape_X, shape_field,detwei=detwei)
    
    ! Matrix entry:
    call addto(mass, ele_field, ele_field, &
         shape_shape(shape_field,shape_field,detwei))
    
  end subroutine assemble_mass_element_contribution

   subroutine get_weighted_dg_inverse_mass_matrix(inverse_mass,dg_mesh, &
       & positions,weights,scalar,dirichlet_list,dirichlet_flag)
    type(mesh_type), intent(in) :: dg_mesh
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: weights
    type(dynamic_csr_matrix), intent(inout) :: inverse_mass
    integer, dimension(:), intent(in), optional :: dirichlet_list
    integer, intent(in), optional :: dirichlet_flag
    real, intent(in), optional :: scalar

    !locals
    integer :: ele
    real :: sc
    logical, dimension(:), allocatable :: internal_dirichlet_list
    integer, dimension(:), pointer :: e_nodes
    integer :: l_dirichlet_flag

    if(present(dirichlet_flag)) then
       l_dirichlet_flag = dirichlet_flag
    else
       l_dirichlet_flag = 0
    end if

    if (.not. present(scalar) ) then 
       sc=1.0
    else
       sc=scalar
    end if

    if(l_dirichlet_flag.ne.DIRICHLET_NONE) then
       if(present(dirichlet_list)) then
          allocate( internal_dirichlet_list( node_count(dg_mesh) ) )
          internal_dirichlet_list = .false.
          internal_dirichlet_list(dirichlet_list) = .true.
       end if
    end if

    assert(dg_mesh%continuity==-1)

    do ele = 1, dg_mesh%elements

       if(present(dirichlet_list).and. &
            & (l_dirichlet_flag.ne.DIRICHLET_NONE)) then
          e_nodes => ele_nodes(dg_mesh,ele)
          call assemble_weighted_local_dg_inverse_mass_matrix(inverse_mass, &
               dg_mesh,positions,weights,sc,ele,internal_dirichlet_list(e_nodes), &
               l_dirichlet_flag)
       else
          call assemble_weighted_local_dg_inverse_mass_matrix(inverse_mass, &
               dg_mesh,positions,weights,sc,ele)
       end if
    end do

  end subroutine get_weighted_dg_inverse_mass_matrix

  subroutine assemble_weighted_local_dg_inverse_mass_matrix(&
       inverse_dynamic_mass, &
       dg_mesh,positions,weights,sc,ele,dirichlet_list,dirichlet_flag)
    integer, intent(in) :: ele
    real, intent(in) :: sc
    type(dynamic_csr_matrix), intent(inout) :: inverse_dynamic_mass
    type(mesh_type), intent(in) :: dg_mesh
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: weights
    logical, dimension(:), intent(in), optional :: dirichlet_list
    integer, intent(in), optional :: dirichlet_flag

    !local variables
    real, dimension(ele_loc(dg_mesh,ele),ele_loc(dg_mesh,ele)) :: local_mass
    real, dimension(dg_mesh%shape%dim,ele_loc(dg_mesh,ele)) :: X_ele
    real, dimension(dg_mesh%shape%ngi) :: detwei
    integer, dimension(:), pointer :: ele_dg
    type(element_type), pointer :: shape_dg,shape_X
    integer :: i

    !assemble local mass matrix
    ele_dg=>ele_nodes(dg_mesh, ele)
    shape_dg=>ele_shape(dg_mesh, ele)
    shape_X=>ele_shape(positions, ele)
    X_ele=ele_val(positions, ele)

    call transform_to_physical(X_ele,shape_X,detwei=detwei)
    local_mass = shape_shape(shape_dg,shape_dg,&
         sc*detwei*ele_val_at_quad(weights,ele))

    if(present(dirichlet_list)) then
       do i = 1, size(dirichlet_list)
          if(dirichlet_list(i)) then
             select case(dirichlet_flag)
             case (DIRICHLET_NONE)
             case (DIRICHLET_ONES_ON_DIAGONAL)
                local_mass(:,i) = 0.
                local_mass(i,:) = 0.
                local_mass(i,i) = 1.0
             case (DIRICHLET_BIG_SPRING)
                local_mass(i,i) = INFINITY
             case default
                FLAbort('bad dirichlet flag')
             end select
          end if
       end do
    end if

    call invert(local_mass)

    call set(inverse_dynamic_mass,ele_dg,ele_dg,local_mass)

  end subroutine assemble_weighted_local_dg_inverse_mass_matrix

  subroutine get_ops(C,positions,&
        rho,s_field,u_field,v_field,bottom)
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: s_field(:), bottom
    type(vector_field), intent(inout) :: u_field, v_field
    real, intent(in)  :: rho(:) 
    type(block_dynamic_csr_matrix), intent(inout) :: C

    ! We form a mass matrix
    !local
    integer :: ele, i,j,k,l, ierr
    integer, dimension(:), pointer :: ptr

     do i=1,2*n_layers
        do j=1,2*n_layers
           call zero(C%blocks(i,j))
        end do
     end do
    
    ewrite(1,*) 'subroutine get_ops'

    ! Assemble A element by element.
    do ele=1,element_count(s_field(1))

       call assemble_op_element_contribution( &
            C, positions, ele,rho,s_field,u_field,v_field,bottom)
       
    end do
   
    ewrite(1,*) 'END subroutine get_ops'

  end subroutine get_ops

subroutine assemble_op_element_contribution( &
       C_op, positions, ele,rho, s_field,m_field,u_field,bottom)
    type(block_dynamic_csr_matrix), intent(inout) :: C_op
    type(vector_field), intent(in) :: positions
    real, intent(in) :: rho(:)
    type(scalar_field), intent(in) :: s_field(:), bottom    
    type(vector_field), intent(in) :: m_field,u_field
    integer, intent(in) :: ele
    integer :: i,j,k,kk
    integer, pointer, dimension(:) :: ele_2,row
    real :: ell
    ! Locations of nodes.
    real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele, X_ele2
    ! Locations of quadrature points
    real, dimension(positions%dim,ele_ngi(positions,ele)) :: X_quad
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei, detwei2, hloc
    real, dimension(ele_loc(s_field(1),ele),n_layers) :: ele_h
    real, dimension(ele_ngi(s_field(1),ele),n_layers) :: dlocgi
    real, dimension(ele_ngi(s_field(1),ele),0:n_layers) :: hlocgi
     real, dimension(ele_loc(bottom,ele)) :: ele_b
    real, dimension(2,ele_ngi(positions,ele),n_layers) :: grad_dlocgi
    real, dimension(2,0:n_layers,ele_ngi(positions,ele)) :: grad_sum_h
    real, dimension(ele_loc(m_field,ele),&
         ele_ngi(positions,ele),positions%dim) :: dm_t
    real, dimension(ele_loc(u_field,ele),&
         ele_ngi(positions,ele),positions%dim) :: du_t
    real, dimension(ele_loc(s_field(1),ele),&
         ele_ngi(positions,ele),positions%dim) :: dh_t
    real, dimension(2,&
         ele_loc(m_field,ele),ele_loc(u_field,ele)) :: Q_muT
      real, dimension(2,&
         m_field%mesh%faces%shape%loc,u_field%mesh%faces%shape%loc) :: Q_muT_f
      real, dimension(m_field%mesh%faces%shape%ngi) :: detwei_f,hlocgi_f
      real, dimension(2,m_field%mesh%faces%shape%ngi) :: normal
    ! Node numbers of field element.
    integer, dimension(:), pointer :: ele_u_field,ele_m_field, ele_field_2
    integer, dimension(:), pointer :: ele_n_u,ele_n_mesh
    integer :: ele_u_f(u_field%mesh%faces%shape%loc),&
         ele_m_f(m_field%mesh%faces%shape%loc)
    integer :: face_h, face_m, face_u,face_mesh
    real :: X_ele_f(2,positions%mesh%faces%shape%loc)

    ! Shape functions.
    type(element_type), pointer :: shape_u_field, shape_m_field,shape_X,shape_h
    type(element_type), pointer :: u_shape_f, m_shape_f,shape_Xf

 
    ele_m_field=>ele_nodes(m_field, ele)
    ele_u_field=>ele_nodes(u_field, ele)
    shape_m_field=>ele_shape(m_field, ele)
    shape_u_field=>ele_shape(u_field, ele)
    shape_X=>ele_shape(positions, ele)
    shape_h=>ele_shape(s_field(1), ele)
    ele_n_mesh=>ele_neigh(positions,ele)
    ele_n_u=>ele_neigh(u_field,ele)
   

    ! Locations of local vertices.
    X_ele=ele_val(positions, ele)

    ! Locations of quadrature points.
!    X_quad=ele_val_at_quad(positions, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(X_ele, shape_X,&
         m=shape_m_field,&
         dn_t=du_t,dm_t=dm_t,detwei=detwei)


    ele_b=ele_val(bottom,ele)
    do k=1,n_layers
       ele_h(:,k)=ele_val(s_field(k),ele)
       Dlocgi(:,k)=0.0!ele_val_at_quad(s_field(k),ele)
    end do
    do k=1,n_layers
       hlocgi(:,k-1)=ele_val_at_quad(bottom,ele)-sum(Dlocgi(:,k:n_layers),2)
       do i=1,2
          grad_dlocgi(i,:,k)=matmul(ele_h(:,k),dm_t(:,:,i))
          grad_sum_h(i,k-1,:)=&
               matmul(ele_b-sum(ele_h(:,k:n_layers),2),dm_t(:,:,i))
       end do
    end do
    hlocgi(:,n_layers)=ele_val_at_quad(bottom,ele)
    do i=1,2
       grad_sum_h(i,n_layers,:)=&
            matmul(ele_b,dm_t(:,:,i))
    end do

    do k=1,n_layers
       do i=1,2
!          call addto(C_op%blocks(1+2*(k-1),i+2*(k-1)),&
!               ele_m_field,ele_u_field,&
!               shape_shape(shape_m_field,shape_u_field,&
!               detwei*grad_sum_h(i,k,:)))
!          call addto(C_op%blocks(2+2*(k-1),i+2*(k-1)),&
!               ele_m_field,ele_u_field,&
!               shape_shape(shape_m_field,shape_u_field,&
!               detwei*grad_sum_h(i,k-1,:)))
       end do
!          Q_muT=-dshape_shape(dm_t,shape_u_field,detwei*hlocgi(:,k))
!          do i=1,2
!             call addto(C_op%blocks(1+2*(k-1),i+2*(kk-1))&
!                  ,ele_m_field,ele_u_field,&
!                  Q_muT(i,:,:))
!          end do
!          Q_muT=-dshape_shape(dm_t,shape_u_field,detwei*hlocgi(:,k-1))
!          do i=1,2
!             call addto(C_op%blocks(1+2*(k-1),i+2*(kk-1))&
!                  ,ele_m_field,ele_u_field,&
!                  Q_muT(i,:,:))
!          end do
       do kk=k+1, n_layers
          Q_muT=-dshape_shape(dm_t,shape_u_field,detwei*Dlocgi(:,kk))
!             Q_muT=shape_dshape(shape_m_field,du_t,detwei*hlocgi(:,kk))&
!                  +shape_shape_vector(shape_m_field,shape_u_field,detwei,&
!                  grad_Dlocgi(:,:,kk))
          do i=1,2
             call addto(C_op%blocks(1+2*(k-1),i+2*(kk-1))&
                  ,ele_m_field,ele_u_field,&
                  Q_muT(i,:,:))
          end do
       end do
       do kk=k,n_layers
          Q_muT=-dshape_shape(dm_t,shape_u_field,detwei*Dlocgi(:,kk))
!             Q_muT=shape_dshape(shape_m_field,du_t,detwei*hlocgi(:,kk))&
!                  +shape_shape_vector(shape_m_field,shape_u_field,detwei,&
!                  grad_Dlocgi(:,:,kk))
          do i=1,2
             call addto(C_op%blocks(2+2*(k-1),i+2*(kk-1))&
                  ,ele_m_field,ele_u_field,&
                  Q_muT(i,:,:))
          end do
       end do

       do j=1,size(ele_n_u)
          face_mesh=ele_face(positions,ele,ele_n_mesh(j))
          X_ele_f=face_val(positions,face_mesh)
          shape_Xf=>face_shape(positions,face_mesh)
          face_m=ele_face(m_field,ele,ele_n_u(j))
          m_shape_f=>face_shape(m_field,face_m)
          face_u=ele_face(u_field,ele,ele_n_u(j))
          u_shape_f=>face_shape(u_field,face_u)
          face_h=ele_face(s_field(1),ele,ele_n_u(j))
          ele_m_f=face_global_nodes(m_field, face_m)
          ele_u_f=face_global_nodes(u_field, face_u)
  

          call transform_bdy_to_physical(X_ele,X_ele_f, shape_X, shape_Xf,&
               detwei_f=detwei_f,normal=normal)

          do kk=k+1,n_layers
             hlocgi_f=face_val_at_quad(s_field(kk),face_h)
             Q_mut_f=-shape_shape_vector(m_shape_f,u_shape_f,&
                  hlocgi_f*detwei_f,normal)
             do i=1,2
!                call addto(C_op%blocks(1+2*(k-1),i+2*(kk-1))&
!                     ,ele_m_f,ele_u_f,&
!                     Q_muT_f(i,:,:))
             end do
          end do
          do kk=k,n_layers
               hlocgi_f=face_val_at_quad(s_field(kk),face_h)
             Q_mut_f=-shape_shape_vector(m_shape_f,u_shape_f,&
                  hlocgi_f*detwei_F,normal)
             do i=1,2
!             call addto(C_op%blocks(2+2*(k-1),i+2*(kk-1))&
!                  ,ele_m_f,ele_u_f,&
!                  Q_muT_f(i,:,:))
             end do
          end do
          
          
       end do

    end do



  end subroutine assemble_op_element_contribution

subroutine MLCM_timestep(MLCM_mat,mesh,bcs,u,u_bar,D,H,m,bottom,rho,dt)
  type(dg_mesh)      :: mesh
  type(bc_info) :: bcs
  type(vector_field) :: u(:),m(:),u_bar
  type(scalar_field) :: D(:),bottom,H
  real, dimension(:) :: rho
  type(block_csr_matrix)  :: MLCM_mat
  real :: dt
! locals
  real, allocatable,dimension(:,:)   :: D_old ,m_rhs,u_bar_old
  real, allocatable,dimension(:,:,:) :: u_old
  real, allocatable, dimension(:) :: H_old, delta_H_0,lambda
  integer :: k, step

  allocate (H_old(node_count(H)),delta_H_0(node_count(H)),&
       lambda(node_count(m(1))),u_bar_old(node_count(u_bar),2))

! RK2 timestepping

! Half step

  allocate(D_old(node_count(D(1)),n_layers),&
       u_old(node_count(u(1)),2,n_layers))
  
  if (RIGID_LID) then
     allocate( m_rhs(node_count(m(1)),2*n_layers+1))
  else
    allocate( m_rhs(node_count(m(1)),2*n_layers))
  end if

  do step=1,2

     delta_H_0=0.0


!  call barotropic_H_step(mesh,H,H_old,u_bar,dt,step)

if (step==1) then
  do k=1,n_layers
     call get_weighted_dg_inverse_mass_matrix(&
          mesh%m_inv%blocks(1,1+2*(k-1)),&
          u(1)%mesh,mesh%positions,D(k),rho(k))
     call get_weighted_dg_inverse_mass_matrix(&
          mesh%m_inv%blocks(1,2+2*(k-1)),&
          u(1)%mesh,mesh%positions,D(k),rho(k))
  end do
end if

  call step_D(mesh,D,D_old,u,dt,step,rho,H%val)

!  call get_ops(mesh%C,mesh%positions,&
!       rho,D,m(1),u(1),bottom)

if (.true.) then
  do k=1,n_layers
     call get_weighted_dg_inverse_mass_matrix(&
          mesh%m_inv%blocks(1,1+2*(k-1)),&
          u(1)%mesh,mesh%positions,D(k),rho(k))
     call get_weighted_dg_inverse_mass_matrix(&
          mesh%m_inv%blocks(1,2+2*(k-1)),&
          u(1)%mesh,mesh%positions,D(k),rho(k))
  end do
end if
!  call get_weighted_dg_inverse_mass_matrix(&
!          mesh%m_inv%blocks(1,1+2*n_layers),&
!          u(1)%mesh,mesh%positions,H,rho(1))
!  call get_weighted_dg_inverse_mass_matrix(&
!          mesh%m_inv%blocks(1,2+2*n_layers),&
!          u(1)%mesh,mesh%positions,H,rho(1))
  call get_MlCM_matrix(MLCM_mat,mesh,D,D_old,u,m,bottom,rho)
  call get_MLCM_rhs(mesh,D,D_old,u,m,bottom,rho,m_rhs,H%val,delta_H_0)
  call solve_MLCM_matrix(MLCM_mat,mesh,D,u,rho,m_rhs,m,lambda)

!  call get_udash_udash(mesh,D,D_old,bottom,m,u,u_bar,rho)

!  call barotropic_u_step(mesh,D,D_old,H,H_old,u,u_bar,u_bar_old,m,&
!       bottom,rho,dt,step)
  call step_u(mesh,bcs,u,u_old,u_bar,D,D_old,H,m,bottom,rho,dt,step)

  if(sig_int) then
     write(0,*) 'stopping due to failed matrix solve'
     stop
  end if

end do



end subroutine MLCM_timestep

subroutine get_MLCM_matrix(MLCM_mat,mesh,D,D_old,u,m,bottom,rho)
  type(dg_mesh)      :: mesh
  type(vector_field) :: u(:),m(:)
  type(scalar_field) :: D(:),bottom
  real, dimension(:) :: rho
  real, dimension(:,:) :: D_old
  type(block_csr_matrix)  :: MLCM_mat

  ! locals

  integer :: i,j,k,ele

  ewrite(1,*) 'subroutine get_MLCM_matrix'

  call zero(MLCM_mat)
  call zero(mesh%MLCM_mat_block)
  call zero(mesh%CMT)
  call zero(mesh%CMC)
  call zero(mesh%MLCM_mat_block)

  if ((.not. Shallow) .and. RIGID_LID) then

   do i=1,2*n_layers
       do j=1,2*n_layers
          do k=1,2*n_layers
             call dcsr_matmul_T(mesh%CMT,mesh%C%blocks(j,k),&
                  mesh%m_inv%blocks(1,k))
             call dcsr_matmul_T(mesh%CMC,mesh%C%blocks(i,k),mesh%CMT)
             call d2c(mesh%CMC,mesh%MLCM_mat_block)
             call addto(MLCM_mat,i,j,mesh%MLCM_mat_block)
             call zero(mesh%CMT)
             call zero(mesh%CMC)
             call zero(mesh%MLCM_mat_block)
          end do
       end do
    end do

  end if

     ! Assemble A element by element.
    do ele=1,element_count(m(1))
       
       call assemble_MLCM_element_contribution( &
            MLCM_mat, mesh%positions, ele,rho,D,D_old,u,m,bottom)
       
    end do


    call check_symmetry(MLCM_mat,m(1))


     ewrite(1,*) 'END subroutine get_MLCM_matrix'

end subroutine get_MlCM_matrix

 subroutine assemble_MLCM_element_contribution( &
       MLCM_mat, positions, ele,rho, D,D_old,u,m,bottom)
   type(block_csr_matrix), intent(inout) :: MLCM_mat
   type(vector_field), intent(in) :: positions,u(:),m(:)
   type(scalar_field) :: D(:), bottom
   real, intent(in) :: rho(:)
   real, dimension(:,:) :: D_old
   integer, intent(in) :: ele

   !Locals

   real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
   real, dimension(ele_ngi(positions,ele)) :: detwei,hlocgi
   real, dimension(ele_loc(D(1),ele),n_layers) :: Dloc
   real, dimension(ele_ngi(D(1),ele),n_layers) :: Dlocgi
   real, dimension(2,0:n_layers,ele_ngi(D(1),ele)) :: grad_sum_h
   real, dimension(ele_loc(m(1),ele),ele_ngi(m(1),ele),2) :: dm_t
   real, dimension(ele_loc(D(1),ele)) :: botloc
   real, dimension(:,:,:,:), allocatable :: M_loc
   integer, dimension(:), pointer :: ele_field

   real, dimension(2) :: D_av=(/50.0,350.0/)

   type(element_type), pointer :: shape_m, shape_X,shape_h

   integer :: i,k,kk,kkk

   ele_field=>ele_nodes(m(1), ele)
   shape_X=>ele_shape(positions, ele)
   shape_m=>ele_shape(m(1), ele)
   shape_h=>ele_shape(D(1), ele)

    if (RIGID_LID) then
       allocate(M_loc(2*n_layers+1,2*n_layers+1,&
         ele_loc(m(1),ele),ele_loc(m(1),ele)))
    else
       allocate(M_Loc(2*n_layers,2*n_layers,&
         ele_loc(m(1),ele),ele_loc(m(1),ele)))
    end if

    ! Locations of local vertices.
    X_ele=ele_val(positions, ele)
    call transform_to_physical(X_ele, shape_X,m=shape_m,&
         dm_t=dm_t,detwei=detwei)

    botloc=ele_val(bottom,ele)
    do k=1, n_layers
       if (.not. LINEAR) then
          Dloc(:,k)=D_old(ele_nodes(D(k),ele),k)
       else
          Dloc(:,k)=d_av(k)
       end if
        Dlocgi(:,k)=matmul(Dloc(:,k),shape_h%n)
    end do
    do i=1,2
       do k=1,n_layers
          grad_sum_h(i,k-1,:)=&
               matmul(botloc-sum(Dloc(:,k:n_layers),2),dm_t(:,:,i))
       end do
       grad_sum_h(i,n_layers,:)=&
            matmul(botloc,dm_t(:,:,i))
    end do

! top surface check
!    grad_sum_h(:,0,:)=0.0


    M_loc=0.0

     do k=1,n_layers

        if (.not. Linear) then
!       hlocgi=1.0/(rho(k)*ele_val_at_quad(D(k),ele))
           hlocgi=1.0/(rho(k)*matmul(D_old(ele_nodes(D(k),ele),k),shape_h%n))
        else
           hlocgi=1.0/(rho(k)*d_av(k))
        end if

        M_loc(1+2*(k-1),1+2*(k-1),:,:)=(&
                  4.0*shape_shape(shape_m,shape_m,hlocgi*detwei))
        M_loc(1+2*(k-1),2+2*(k-1),:,:)=(&
               -2.0*shape_shape(shape_m,shape_m,hlocgi*detwei))
        M_loc(2+2*(k-1),1+2*(k-1),:,:)=(&
               -2.0*shape_shape(shape_m,shape_m,hlocgi*detwei))
        M_loc(2+2*(k-1),2+2*(k-1),:,:)=(&
               4.0*shape_shape(shape_m,shape_m,hlocgi*detwei))
     end do

     if (RIGID_LID) then
        hlocgi=1.0/(rho(1)*ele_val_at_quad(D(1),ele))
        M_loc(1,2*n_layers+1,:,:)=-2.0*shape_shape(shape_m,&
             shape_m,hlocgi*detwei)
        M_loc(2,2*n_layers+1,:,:)=4.0*shape_shape(shape_m,&
             shape_m,hlocgi*detwei)
        M_loc(2*n_layers+1,1,:,:)=-2.0*shape_shape(shape_m,&
             shape_m,hlocgi*detwei)
        M_loc(2*n_layers+1,2,:,:)=4.0*shape_shape(shape_m,&
             shape_m,hlocgi*detwei)
        M_loc(2*n_layers+1,2*n_layers+1,:,:)=4.0*shape_shape(shape_m,&
             shape_m,hlocgi*detwei)
     END if

     if ((.not. Shallow ).and. (.not. RIGID_LID)) then

     do k=1,n_layers
!        hlocgi=1.0/(rho(k)*ele_val_at_quad(D(k),ele))
!         hlocgi=1.0/(rho(k)*d_av(k)) 
           hlocgi=1.0/(rho(k)*matmul(D_old(ele_nodes(D(k),ele),k),shape_h%n))              

                ! First < AW,AP>

                if (.true.) then

                ! <w grad_h.grad_h p>

                   M_loc(1+2*(k-1),1+2*(k-1),:,:)=&
                        M_loc(1+2*(k-1),1+2*(k-1),:,:)+(&
                     shape_shape(shape_m,shape_m,&
                     sum(grad_sum_h(:,k,:)*grad_sum_h(:,k,:),1)&
                     *hlocgi*detwei))

                ! <-w grad_h.D sum_j grad u> + <-sum_jgrad w D, grad_h.u> 

                do kk=1,k-1
!                   M_loc(1+2*(k-1),1+2*kk-1,:,:)=&
!                        M_loc(1+2*(k-1),1+2*(kk-1),:,:)+(&
!                       -shape_vector_dot_dshape(&
!                        shape_m,grad_sum_h(:,k,:),dm_t,&
!                        detwei/rho(k)))


                   M_loc(1+2*(kk-1),1+2*(k-1),:,:)=&
                        M_loc(1+2*(kk-1),1+2*(k-1),:,:)+(&
                       -dshape_dot_vector_shape(&
                        dm_t,grad_sum_h(:,k,:),shape_m,&
                        detwei/rho(k)))  

                   M_loc(1+2*(k-1),1+2*(kk-1),:,:)=&
                        M_loc(1+2*(k-1),1+2*(kk-1),:,:)+transpose(&
                       -dshape_dot_vector_shape(&
                        dm_t,grad_sum_h(:,k,:),shape_m,&
                        detwei/rho(k)))  

                end do
                
                ! < sum_j grad w D,D sum_j grad u
                
                do kk=1,k-1
                   do kkk=1,k-1

                     M_loc(1+2*(kk-1),1+2*(kkk-1),:,:)=&
                           M_loc(1+2*(kk-1),1+2*(kkk-1),:,:)+(&
                           dshape_dot_dshape(dm_t,dm_t,&
                           Dlocgi(:,k)*&
                           detwei/rho(k)))
                   end do
                end do
                

             end if
                   ! <Aw,Bu> and <Bw,Au>

             if (.true.) then

                !<w grad h_i+1 . grad h_i u >
                
                M_loc(1+2*(k-1),2+2*(k-1),:,:)=&
                     M_loc(1+2*(k-1),2+2*(k-1),:,:)+(&
                     shape_shape(shape_m,shape_m,&
                     sum(grad_sum_h(:,k,:)*grad_sum_h(:,k-1,:),1)&
                     *hlocgi*detwei))
                
                M_loc(2+2*(k-1),1+2*(k-1),:,:)=&
                     M_loc(2+2*(k-1),1+2*(k-1),:,:)+(&
                     shape_shape(shape_m,shape_m,&
                     sum(grad_sum_h(:,k-1,:)*grad_sum_h(:,k,:),1)&
                     *hlocgi*detwei))
                
                
                ! <-w grad_h.D sum_j grad u> + <-sum_jgrad w D, grad_h.u> 
                
                do kk=1,k
                   
!                   M_loc(1+2*(k-1),2+2*(kk-1),:,:)=&
!                        M_loc(1+2*(k-1),2+2*(kk-1),:,:)+(&
!                        -shape_vector_dot_dshape(shape_m,&
!                        grad_sum_h(:,k,:),dm_t,&
!                        detwei/rho(k)))
                   
                   M_loc(2+2*(kk-1),1+2*(k-1),:,:)=&
                        M_loc(2+2*(kk-1),1+2*(k-1),:,:)+(&
                        -dshape_dot_vector_shape(dm_t,&
                        grad_sum_h(:,k,:),shape_m,&
                        detwei/rho(k)))

                    M_loc(1+2*(k-1),2+2*(kk-1),:,:)=&
                        M_loc(1+2*(k-1),2+2*(kk-1),:,:)+transpose(&
                        -dshape_dot_vector_shape(dm_t,&
                        grad_sum_h(:,k,:),shape_m,&
                        detwei/rho(k)))

                end do

                 do kk=1,k-1
                   
!                    M_loc(2+2*(k-1),1+2*(kk-1),:,:)=&
!                         M_loc(2+2*(k-1),1+2*(kk-1),:,:)+(&
!                        -shape_vector_dot_dshape(shape_m,&
!                        grad_sum_h(:,k-1,:),&
!                        dm_t,detwei/rho(k))) 
                   
                    M_loc(1+2*(kk-1),2+2*(k-1),:,:)=&
                         M_loc(1+2*(kk-1),2+2*(k-1),:,:)+(& 
                       -dshape_dot_vector_shape(dm_t,&
                        grad_sum_h(:,k-1,:),shape_m,&
                        detwei/rho(k)))

                    M_loc(2+2*(k-1),1+2*(kk-1),:,:)=&
                         M_loc(2+2*(k-1),1+2*(kk-1),:,:)+transpose(& 
                       -dshape_dot_vector_shape(dm_t,&
                        grad_sum_h(:,k-1,:),shape_m,&
                        detwei/rho(k)))

                end do

                   ! < sum_j grad w D,D sum_j grad u

                do kk=1,k-1
                   do kkk=1,k

                      M_loc(1+2*(kk-1),2+2*(kkk-1),:,:)= &
                           M_loc(1+2*(kk-1),2+2*(kkk-1),:,:)+(& 
                           dshape_dot_dshape(dm_t,dm_t,&
                           dlocgi(:,k)*&
                           detwei/rho(k)))

                      M_loc(2+2*(kkk-1),1+2*(kk-1),:,:)=&
                           M_loc(2+2*(kkk-1),1+2*(kk-1),:,:)+(& 
                           dshape_dot_dshape(dm_t,dm_t,&
                           dlocgi(:,k)*&
                           detwei/rho(k)))

                   end do
                end do
                   
             end if
                   
                ! <Bw,Bu>

                if (.true.) then

                ! <w grad_h.grad_h p>

                   M_loc(2+2*(k-1),2+2*(k-1),:,:)=&
                        M_loc(2+2*(k-1),2+2*(k-1),:,:)+(& 
                        shape_shape(shape_m,shape_m,&
                        sum(grad_sum_h(:,k-1,:)*grad_sum_h(:,k-1,:),1)&
                        *hlocgi*detwei))


                ! <-w grad_h.D sum_j grad u> + <-sum_jgrad w D, grad_h.u> 

                do kk=1,k

!                   M_loc(2+2*(k-1),2+2*(kk-1),:,:)=&
!                        M_loc(2+2*(k-1),2+2*(kk-1),:,:)+(&
!                        -shape_vector_dot_dshape(shape_m,&
!                        grad_sum_h(:,k-1,:),&
!                        dm_t,detwei/rho(k))) 

                  M_loc(2+2*(kk-1),2+2*(k-1),:,:)=&
                       M_loc(2+2*(kk-1),2+2*(k-1),:,:)+(&   
                       -dshape_dot_vector_shape(dm_t,&
                       grad_sum_h(:,k-1,:),shape_m,&
                       detwei/rho(k)))

                  M_loc(2+2*(k-1),2+2*(kk-1),:,:)=&
                       M_loc(2+2*(k-1),2+2*(kk-1),:,:)+transpose(&   
                       -dshape_dot_vector_shape(dm_t,&
                       grad_sum_h(:,k-1,:),shape_m,&
                       detwei/rho(k)))

                end do
                
                ! < sum_j grad w D,D sum_j grad u
                
                do kk=1,k
                   do kkk=1,k

                  M_loc(2+2*(kk-1),2+2*(kkk-1),:,:)=&
                       M_loc(2+2*(kk-1),2+2*(kkk-1),:,:)+(&   
                           dshape_dot_dshape(dm_t,dm_t,&
                           dlocgi(:,k)*&
                           detwei/rho(k)))
                   end do
                end do

             end if
          end do
       end if


     do k=1,size(M_loc,1)
        do kk=1,size(M_loc,2)
           call addto(MLCM_mat,k,kk,ele_field, ele_field,&
                M_loc(k,kk,:,:))
        end do
     end do


 end subroutine assemble_MLCM_element_contribution

subroutine get_MLCM_rhs(mesh,D,D_old,u,m,bottom,rho,m_rhs,H_0,dH_0)
  type(dg_mesh)      :: mesh
  type(scalar_field) :: D(:),bottom
  type(vector_field) :: u(:),m(:)
  real :: m_rhs(:,:),H_0(:),dH_0(:), D_old(:,:)
  real, dimension(:) :: rho
  ! locals
  integer :: ele,i,k
  real, dimension(2,node_count(u(1))) :: uu_rhs
  real, dimension(node_count(u(1))) :: uu

  m_rhs=0.0

  do ele=1,ele_count(m(1))     
     call get_MLCM_rhs_element_contributions(mesh,mesh%positions,&
          D,D_old,u,m,bottom,rho,m_rhs,ele)

  end do

  if (RIGID_LID) then
     m_rhs(:,2*n_layers+1)=m_rhs(:,2)
     do k=1,n_layers
        uu_rhs=0.0
        do ele=1,ele_count(u(1))
           call get_uu_element_contribution(mesh,u(k),D(k),D_old(:,k),&
                rho(k),uu_rhs,ele,flux="Mean") 
        end do
        do i=1,2
           do ele=1,node_count(u(1))
              uu(ele)=sum(uu_rhs(i,row_m(mesh%m_inv%blocks(1,i+2*(k-1)),ele))&
                   *row_val_ptr(mesh%M_inv%blocks(1,i+2*(k-1)),ele))
           end do
           do ele=1,node_count(m(1))
              m_rhs(ele,2*n_layers+1)=m_rhs(ele,2*n_layers+1)&
                   +sum(row_val_ptr(mesh%C%blocks(2,i+2*(k-1)),ele)&
                   *(uu(row_m(mesh%C%blocks(2,i+2*(k-1)),ele))))
!                   *(u(k)%val(i)%ptr(&
!                   row_m(mesh%C%blocks(2,i+2*(k-1)),ele))/dt&
!                   +uu(row_m(mesh%C%blocks(2,i+2*(k-1)),ele))))
         end do
      end do
   end do

  end if

end subroutine get_MLCM_rhs

subroutine get_MLCM_rhs_element_contributions(mesh,positions,&
     D,D_old,u,m,bottom,rho,m_rhs,ele)
  type(dg_mesh) :: mesh
  type(scalar_field) :: D(:),bottom
  type(vector_field) :: positions,u(:),m(:)
  real :: m_rhs(:,:), D_old(:,:)
  real, dimension(:) :: rho
  integer :: ele

  !  locals
  real, dimension(ele_loc(D(1),ele),n_layers) :: Dloc
  real, dimension(ele_loc(D(1),ele)) :: Bloc
  real, dimension(2,n_layers,ele_loc(u(1),ele)) :: uloc
  real, dimension(ele_ngi(D(1),ele),n_layers) :: Dlocgi,divulocgi
   real, dimension(ele_ngi(D(1),ele)) :: detwei
  real, dimension(2,n_layers,ele_ngi(D(1),ele)) :: grad_hlocgi,grad_dlocgi
  real, dimension(2,2,n_layers,ele_ngi(u(1),ele)) :: grad_ulocgi
  real, dimension(2,n_layers,ele_ngi(u(1),ele)) :: ulocgi,u_grad_ulocgi
  real, dimension(ele_loc(u(1),ele),ele_ngi(u(1),ele),2) :: du_t
  real, dimension(ele_loc(m(1),ele),ele_ngi(m(1),ele),2) :: dm_t
  real, dimension(ele_loc(D(1),ele),ele_ngi(D(1),ele),2) :: dD_t
  real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
  integer, dimension(:), pointer :: ele_u, ele_m
  real, dimension(2) :: d_av=(/50.0,350.0/)

  type(element_type), pointer :: shape_m, shape_D, shape_X, shape_u

  integer :: i,j,k
  
  ele_u=>ele_nodes(u(1), ele)
  X_ele=ele_val(positions, ele)
  ele_m=>ele_nodes(m(1), ele)
  shape_m=>ele_shape(m(1), ele)
  shape_u=>ele_shape(u(1), ele)
  shape_D=>ele_shape(D(1), ele)
  shape_X=>ele_shape(positions, ele)


  call transform_to_physical(X_ele, shape_X,m=shape_u,dm_t=du_t, &
       detwei=detwei)
  call transform_to_physical(X_ele, shape_X,m=shape_m,dm_t=dm_t, &
       detwei=detwei)
  call transform_to_physical(X_ele, shape_X,m=shape_D,dm_t=dD_t, &
       detwei=detwei)

  u_grad_ulocgi=0.0
  bloc= ele_val(bottom, ele)

 

  do k=1,n_layers
     do i=1,2
        do j=1,ele_loc(u(1),ele)        
           uloc(i,k,j)=sum(u(k)%val(i)%ptr(row_m(mesh%M_inv%blocks(1,i),&
                ele_u(j)))*&
                row_val_ptr(mesh%M_inv%blocks(1,i+2*(k-1)),ele_u(j)))
        end do
     end do
     Dloc(:,k) = D_old(ele_nodes(D(k), ele),k)
     ulocgi(:,k,:) = matmul(uloc(:,k,:),shape_u%n)
     Dlocgi(:,k) = matmul(D_old(ele_nodes(D(k), ele),k),shape_D%n)
     divulocgi(:,k)=matmul(uloc(1,k,:),du_t(:,:,1))+&
          matmul(uloc(2,k,:),du_t(:,:,2))
  end do

  do j=1,2
     grad_hlocgi(j,n_layers,:)= matmul(bloc,dD_t(:,:,j))
  end do
  do  k=1,n_layers-1
     do j=1,2
        grad_hlocgi(j,k,:)=&
             matmul(bloc-sum(Dloc(:,k+1:n_layers),2),dD_t(:,:,j))
     end do
  end do


  do k=1,n_layers
     do j=1,2
        grad_dlocgi(j,k,:)=matmul(Dloc(:,k),dD_t(:,:,j))
     end do
     do i=1,2
        do j=1,2
           grad_ulocgi(i,j,k,:)=matmul(uloc(i,k,:),du_t(:,:,j))
        end do
     end do
     do i=1,2
        u_grad_ulocgi(i,k,:)=&
             sum(ulocgi(:,k,:)*grad_ulocgi(i,:,k,:),1)
     end do
  end do

  do k=1,n_layers
     do i=1,2
        If (.not. LINEAR) then
           m_rhs(ele_m,i+2*(k-1))=m_rhs(ele_m,i+2*(k-1))&
                +shape_rhs(shape_m,g0*detwei)
        else
        m_rhs(ele_m,i+2*(k-1))=m_rhs(ele_m,i+2*(k-1))&
             +shape_rhs(shape_m,g0*detwei*(dlocgi(:,k)-d_av(k))/d_av(k))
        end if
     end do
     
     if (.not. SHALLOW .and. (.not. LINEAR)) then
        
        do i=1,2
           do j=1,2
              m_rhs(ele_m,2+2*(k-1))=m_rhs(ele_m,2+2*(k-1))&
                   +shape_rhs(shape_m,&
                   detwei*Dlocgi(:,k)*(&
                   grad_ulocgi(i,j,k,:)*grad_ulocgi(j,i,k,:)+&
                   grad_ulocgi(i,i,k,:)*grad_ulocgi(j,j,k,:)))
           end do
        end do

        do i=1,2
           m_rhs(ele_m,i+2*(k-1))=m_rhs(ele_m,i+2*(k-1))&
                +dshape_dot_vector_rhs(dm_t,ulocgi(:,k,:),&
                detwei*&
                sum(ulocgi(:,k,:)*grad_hlocgi(:,k,:),1))&
                +shape_rhs(shape_m,&
                detwei*(&
                divulocgi(:,k)&
                *sum(ulocgi(:,k,:)*grad_hlocgi(:,k,:),1)&
                +sum(u_grad_ulocgi(:,k,:)*grad_hlocgi(:,k,:),1)))
              
           do j=k+1,n_layers
              m_rhs(ele_m,i+2*(k-1))=m_rhs(ele_m,i+2*(k-1))&
                   -dshape_dot_vector_rhs(dm_t,&
                   ulocgi(:,j,:),detwei*(&
                   sum(ulocgi(:,j,:)&
                   *grad_dlocgi(:,j,:),1)&
                   +Dlocgi(:,j)* divulocgi(:,j)))&
                   -dshape_dot_vector_rhs(dm_t,&
                   u_grad_ulocgi(:,j,:),detwei*Dlocgi(:,j))
           end do
        end do
     end if
        
  end do

end subroutine get_MLCM_rhs_element_contributions


subroutine solve_MLCM_matrix(MLCM_mat,mesh,D,u,rho,m_rhs,m,lambda)
  type(dg_mesh)      :: mesh
  type(vector_field) :: u(:),m(:)
  type(scalar_field) :: D(:)
  real, dimension(:) :: rho
  real :: m_rhs(:,:),lambda(:)
  type(block_csr_matrix)  :: MLCM_mat

  integer :: i,k,m_nodes
  real, allocatable :: m_out(:)

  allocate(m_out(size(m_rhs)))

  m_nodes=node_count(m(1))
  do i=1,2
     do k=1,n_layers

        m_out(1+(2*(k-1)+i-1)*m_nodes:(2*(k-1)+i)*m_nodes)=&
             m(k)%val(i)%ptr

     end do
  end do

    call petsc_solve(m_out,&
         MLCM_mat,&
         reshape(m_rhs,(/size(m_rhs)/)),&
         abstol=1.0e-30,&
         max_its=10000,&
         startfromzero=.false.)

    do i=1,2
     do k=1,n_layers
         m(k)%val(i)%ptr=&
              m_out(1+(2*(k-1)+i-1)*m_nodes:(2*(k-1)+i)*m_nodes)
      end do
   end do

   if (RIGID_LID)   lambda= m_out(1+2*n_layers*m_nodes:(2*n_layers+1)*m_nodes)

end subroutine solve_MlCM_matrix

subroutine barotropic_u_step(mesh,D,D_old,H,H_old,u,u_bar,u_bar_old,m,&
     bottom,rho,dt,step)
  type(dg_mesh) :: mesh
  type(vector_field) :: u_bar,u(:),m(:)
  type(scalar_field) :: H,D(:),bottom
  real, dimension(:,:) :: u_bar_old
  real, dimension(:) :: H_old
  real :: D_old(:,:),rho(:),dt,c(D(1)%mesh%shape%loc)
  integer :: step
! local
  integer :: i,j,k,ele
  real, dimension(2,node_count(u_bar)) :: u_rhs
  real, dimension(node_count(u_bar)):: u1,u2
  real, dimension(n_layers,node_count(m(1))):: p_out
  real, dimension(H%mesh%shape%loc) :: dloc,Hloc
  integer, pointer, dimension(:) :: ele_h

  u_rhs=0.0
  p_out=0.0




        do k=1,n_layers


           u1=u(k)%val(1)%ptr-u_bar%val(1)%ptr
           u2=u(k)%val(2)%ptr-u_bar%val(2)%ptr


     do ele=1,ele_count(u_bar)           
!           c=sqrt((rho(k)-rho(1))/rho(1)*g0*&
!                sum(ele_val(D(1),ele)*ele_val(D(2),ele)/&
!                (ele_loc(D(k),ele)*(ele_val(D(1),ele)+ele_val(D(2),ele)))))

if (n_layers==1) then
   c=sqrt(g0*ele_val(D(k),ele))
else
   ele_h=>ele_nodes(D(1),ele)
   if(k==1) then
      c=sqrt((rho(2)-rho(1))/rho(1)*g0*&
           D_old(ele_h,1))
   else
      c=sqrt((rho(2)-rho(1))/rho(2)*g0*&
           D_old(ele_h,2))
   end if
end if

           

           call get_uu_element_contribution(mesh,u_bar,D(k),D_old(:,k),&
                rho(k),u_rhs,ele,flux="None",u1=u1,u2=u2,g=(rho(2)-rho(1))*g0/rho(k))

        end do

!do k=1,n_layers
!        u1=u(k)%val(1)%ptr-u_bar%val(1)%ptr
!        u2=u(k)%val(2)%ptr-u_bar%val(2)%ptr
!        call get_u_dot_u(mesh,m(k),u(k),u1,u2,rho(k),p_out(k,:))
!end do


     end do


     do ele=1,ele_count(u_bar)

        ele_h=> ele_nodes(H,ele)
        c=sqrt(g0*H_old(ele_h))

        call get_uu_element_contribution(mesh,u_bar,H,H_old,&
             rho(1),u_rhs,ele,c=c,flux="HLLC",g=g0)

        call get_barotropic_mu_element_contribution(mesh,m,u_bar,D,D_old,&
             bottom,p_out,u_rhs,ele)

     end do

     if(step == 1) then

        do i=1,2
           u_bar_old(:,i)=u_bar%val(i)%ptr
           do j=1,node_count(u_bar)
              u_bar%val(i)%ptr(j)=u_bar%val(i)%ptr(j)&
                   +dt*sum(u_rhs(i,row_m(mesh%M_inv%blocks(1,i+2*n_layers),j))*&
                   row_val_ptr(mesh%M_inv%blocks(1,i+2*n_layers),j))
           end do
        end do
  else if (step == 2) then
     do i=1,2
        do j=1,node_count(u_bar)
              u_bar%val(i)%ptr(j)=0.5*(u_bar_old(j,i)+u_bar%val(i)%ptr(j))&
                   +0.5*dt*sum(u_rhs(i,row_m(&
                   mesh%M_inv%blocks(1,i+2*n_layers),j))*&
                   row_val_ptr(mesh%M_inv%blocks(1,i+2*n_layers),j))
           end do
        end do
      
     end if

end subroutine barotropic_u_step

subroutine get_barotropic_mu_element_contribution(mesh,m,u_bar,D,D_old,&
                bottom,p_out,u_rhs,ele)

  type(dg_mesh) :: mesh
  type(vector_field) :: m(:),u_bar
  type(scalar_field) :: D(:),bottom
  real, dimension(:,:) :: u_rhs,D_old,p_out
  integer :: ele

  integer :: i,j,k,kk
  real, dimension(2,ele_loc(mesh%positions,ele)) :: X_ele
  real, dimension(2,mesh%positions%mesh%faces%shape%loc) :: X_ele_f
  real, dimension(ele_ngi(D(1),ele)) :: detwei,plocgi
  real, dimension(D(1)%mesh%faces%shape%ngi) :: detwei_f,plocgi_f
  real, dimension(2,ele_ngi(m(1),ele)) :: psilocgi
   real, dimension(2,m(1)%mesh%faces%shape%ngi) :: psilocgi_f,normal

  integer, dimension(:), pointer :: ele_u,ele_2,ele_mesh
  real, dimension(ele_loc(D(1),ele),ele_ngi(u_bar,ele),2) :: dh_t
  real, dimension(ele_loc(u_bar,ele),ele_ngi(u_bar,ele),2) :: du_t
  type(element_type), pointer :: shape_u, shape_X,shape_h,shape_m
  type(element_type), pointer :: shape_uf, shape_Xf,shape_hf,shape_mf
  real, dimension(ele_ngi(D(1),ele),n_layers+1) :: dlocgi
  real, dimension(D(1)%mesh%faces%shape%ngi,n_layers+1) :: dlocgi_f
  real, dimension(2,ele_ngi(D(1),ele)) :: grad_blocgi
  real, dimension(2,2,ele_ngi(D(1),ele)) ::grad_psilocgi
  real, dimension(2,ele_ngi(D(1),ele),n_layers+1) :: grad_dlocgi
   
  integer, dimension(D(1)%mesh%faces%shape%loc) :: ele_h_f
  integer, dimension(u_bar%mesh%faces%shape%loc) :: ele_u_f
   integer, dimension(m(1)%mesh%faces%shape%loc) :: ele_m_f

   integer :: face_mesh,face_u,face_m,face_h

  shape_u=>ele_shape(u_bar,ele)
  shape_h=>ele_shape(D(1),ele)
  shape_m=>ele_shape(m(1),ele)
  ele_u=>ele_nodes(u_bar,ele)
  shape_X=>ele_shape(mesh%positions,ele)
  X_ele=ele_val(mesh%positions,ele)
  ele_2=>ele_neigh(u_bar,ele)
  ele_mesh=>ele_neigh(mesh%positions,ele)

  call transform_to_physical(X_ele, shape_X,m=shape_h,&
       dn_t=du_t,dm_t=dh_t,detwei=detwei)

  dlocgi=0.0
  grad_dlocgi=0.0

  do kk=1,n_layers
     dlocgi(:,kk)=matmul(sum(D_old(ele_nodes(D(1),ele),kk:n_layers),2),&
          shape_h%n)
     do i=1,2
        grad_dlocgi(i,:,kk)=matmul(sum(D_old(ele_nodes(D(1),ele),&
             kk:n_layers),2),&
             dh_t(:,:,i))
     end do
  end do
  
  do i=1,2
        grad_blocgi(i,:)=matmul(ele_val(bottom,ele),dh_t(:,:,i))
  end do

  do k=1,n_layers
     psilocgi(:,:)=ele_val_at_quad(m(k),ele)
!     plocgi=matmul(D_old(ele_nodes(D(k),ele),k),shape_h%n)&
!          *matmul(p_out(k,ele_nodes(m(k),ele)),shape_m%n)
     psilocgi(1,:)=psilocgi(1,:)-matmul(p_out(k,ele_nodes(m(k),ele)),shape_m%n)
     psilocgi(2,:)=psilocgi(2,:)+matmul(p_out(k,ele_nodes(m(k),ele)),shape_m%n)

     do i=1,2
        grad_psilocgi(i,:,:)=matmul(ele_val(m(k),ele),dh_t(:,:,i))
     end do

!     u_rhs(:,ele_u)= u_rhs(:,ele_u)+dshape_rhs(du_t,detwei*plocgi)
     u_rhs(:,ele_u)= u_rhs(:,ele_u)+dshape_rhs(du_t,detwei*dlocgi(:,k+1)&
          *psilocgi(1,:))
     u_rhs(:,ele_u)= u_rhs(:,ele_u)+dshape_rhs(du_t,detwei*dlocgi(:,k)&
           *psilocgi(2,:))
     u_rhs(:,ele_u)= u_rhs(:,ele_u)+shape_vector_rhs(shape_u,grad_blocgi,&
           detwei*sum(ele_val_at_quad(m(k),ele),1))
   end do


   do j=1,size(ele_2)
      face_mesh=ele_face(mesh%positions,ele,ele_mesh(j))
      X_ele_f=face_val(mesh%positions,face_mesh)
      shape_Xf=>face_shape(mesh%positions, face_mesh)
      call transform_bdy_to_physical(X_ele,X_ele_f, shape_X, shape_Xf,&
           detwei_f=detwei_f,normal=normal)
      face_h=ele_face(D(1),ele,ele_2(j))
      ele_h_f=face_global_nodes(D(1),face_h)
      shape_hf=>face_shape(D(1),face_h)
      face_u=ele_face(u_bar,ele,ele_2(j))
      shape_uf=>face_shape(u_bar,face_u)
      ele_u_f=face_global_nodes(u_bar,face_u)
      face_m=ele_face(m(1),ele,ele_2(j))
      ele_m_f=face_global_nodes(m(1),face_m)
      shape_mf=>face_shape(m(1),face_m)


      dlocgi_f=0.0

       do kk=1,n_layers
          dlocgi_F(:,kk)=matmul(sum(D_old(ele_h_f,kk:n_layers),2),&
               shape_hf%n)
       end do

      do k=1,n_layers
!         plocgi_f=matmul(D_old(ele_h_f,k),shape_hf%n)&
!              *matmul(p_out(k,ele_m_f),shape_mf%n)
         psilocgi_f=face_val_at_quad(m(k),face_m)
         psilocgi_F(1,:)=psilocgi_f(1,:)-matmul(p_out(k,ele_m_f),shape_mf%n)
         psilocgi_F(2,:)=psilocgi_f(2,:)+matmul(p_out(k,ele_m_f),shape_mf%n)
!         u_rhs(:,ele_u_f)= u_rhs(:,ele_u_f)-shape_vector_rhs(shape_uf,&
!              normal,detwei_f*plocgi_f)
         u_rhs(:,ele_u_f)= u_rhs(:,ele_u_f)-shape_vector_rhs(shape_uf,&
              normal,detwei_f*dlocgi_f(:,k+1)&
              *psilocgi_f(1,:))
         u_rhs(:,ele_u_f)= u_rhs(:,ele_u_f)-shape_vector_rhs(shape_uf,&
              normal,detwei_f*dlocgi_f(:,k)&
              *psilocgi_f(2,:))

      end do

   end do


 end subroutine get_barotropic_mu_element_contribution

subroutine step_u(mesh,bcs,u,u_old,u_bar,D,D_old,H,m,bottom,rho,dt,step)
  type(dg_mesh) :: mesh
  type(bc_info) :: bcs
  type(vector_field) :: u(:),m(:),u_bar
  type(scalar_field) :: D(:),bottom,H
  real, dimension(:,:,:) :: u_old
  real, dimension(:,:) :: D_old
  real, dimension(:) :: rho
  real :: dt
  integer :: step
! local
  integer :: i,j,k,ele,its
  real, dimension(2,node_count(u(1))) :: u_rhs, u_b, u_step
  real, dimension(D(1)%mesh%shape%loc) :: dloc,Hloc
  real, dimension(node_count(H)) :: H_fix
  real, dimension(n_layers) :: H_factor

  H_factor=(/50.0,350.0/)/400.0


   do k=1,n_layers

      do i=1,2
         do j=1,node_count(u(1))
            if (step == 1) then
             u_step(i,j)=sum(u(k)%val(i)%ptr(row_m(mesh%M_inv%blocks(1,i),j))*&
                   row_val_ptr(mesh%M_inv%blocks(1,i+2*(k-1)),j))
            else if (step ==2) then
              u_step(i,j)=sum(u_old(row_m(mesh%M_inv%blocks(1,i),j),i,k)*&
                   row_val_ptr(mesh%M_inv%blocks(1,i+2*(k-1)),j))
           end if
           end do
        end do
      
      ewrite(1,*) 'Solving for u'

        call get_u_rhs(mesh,u(k),D,D_old,m,bottom,rho,u_b,u_rhs,u_step,k)

  if(step == 1) then

     do i=1,2
        u_old(:,i,k)=u(k)%val(i)%ptr
     end do

!     do its=1,5

     call zero(mesh%u_mat)
     
   do i=1,2
           do j=1,node_count(u(1))

              u_step(i,j)=sum(u(k)%val(i)%ptr(row_m(mesh%M_inv%blocks(1,i),j))*&
                   row_val_ptr(mesh%M_inv%blocks(1,i+2*(k-1)),j))
           end do
        end do

        call get_u_mat(mesh,u_step,u(1),0.5*dt,flux='FLL', c=sqrt(10.0))

!        do i=1,2

!           call petsc_solve(u_step(i,:),&
!                mesh%mass_u,u_rhs(i,:),&
!                startfromzero=.true.)

!        end do


!        u_step=rotate_bdry(u_step,bcs)

        do i=1,2
           !u(k)%val(i)%ptr=u_old(:,i,k)&
           !     +dt*u_step(i,:) 
           do j=1,node_count(u(1))
             u(k)%val(i)%ptr(j)=u_old(j,i,k)&
                  +dt*sum(u_rhs(i,row_m(mesh%mass_u_inv,j))*&
                   row_val_ptr(mesh%mass_u_inv,j))
          end do
        end do


        u(k)%val(2)%ptr=0.0

 !    end do

  else if (step == 2) then

!     do i=1,2
!        u_old(:,i,k)=u(k)%val(i)%ptr
!     end do

!     do its=1,5

        call zero(mesh%u_mat)

        do i=1,2
           do j=1,node_count(u(1))

              u_step(i,j)=sum(u(k)%val(i)%ptr(row_m(mesh%M_inv%blocks(1,i),j))*&
                   row_val_ptr(mesh%M_inv%blocks(1,i+2*(k-1)),j))

           end do
        end do

!       call get_u_mat(mesh,u_step,u(1),0.5*dt,flux='FLL',c=sqrt(10.0))

!        do i=1,2

!           call petsc_solve(u_step(i,:),&
!                mesh%mass_u,u_rhs(i,:),&
!                startfromzero=.true)

!        end do

    
!        u_step=rotate_bdry(u_step,bcs)
        
        do i=1,2
!           u(k)%val(i)%ptr=u_old(:,i,k)&
!                +dt*u_step(i,:)

           do j=1,node_count(u(1))
              u(k)%val(i)%ptr(j)=u_old(j,i,k)&
                   +dt*sum(u_rhs(i,row_m(mesh%mass_u_inv,j))*&
                   row_val_ptr(mesh%mass_u_inv,j))
          end do

        end do

!     end do

!        u(k)%val(2)%ptr=0.0

  end if
end do


IF (RIGID_LID) then

   D(1)%val=sum(D_old,2)

   do k=2,n_layers
      D(1)%val=D(1)%val-D(k)%val
   end do
end IF


 if (BAROTROPIC_SPLIT) then

    call get_u_bar(mesh,u,D,rho,H,u_b)

    H_fix=0.0
    do k=1,n_layers
       H_fix=H_fix+rho(k)*D(k)%val/rho(1)
    end do
    do k=1,n_layers
       D(k)%val=D(k)%val+H_factor(k)*rho(1)*(H%val-H_fix)/rho(k)
    end do
    do k=1,n_layers
       do i=1,2
          u(k)%val(i)%ptr=u(k)%val(i)%ptr+(u_bar%val(i)%ptr-u_b(i,:))
       end do
    end do

end if

end subroutine step_u

subroutine step_D(mesh,D,D_old,u,dt,step,rho,D_0)
  type(dg_mesh) :: mesh
  type(vector_field) :: u(:)
  type(scalar_field) :: D(:)
  real, dimension(:,:) :: D_old
  real, dimension(:) :: D_0
  real :: dt,rho(:)
  integer :: step,start
  integer :: k,ele
  real :: d_rhs(node_count(D(1))),delta_D(node_count(D(1)))
  real :: c(ele_loc(D(1),1)),g,dd(2)=(/50.0,350.0/)
  integer, pointer :: ele_h(:)

!  type(csr_matrix) :: D_mat

!  D_mat=clone(mesh%mass_h)

 do k=1,n_layers

    D_rhs=0.0
    delta_D=0.0
!    call zero(D_mat)


if(step == 1) then
       D_old(:,k)=D(k)%val
!       D(k)%val=D(k)%val+dt*delta_D
    else if (step == 2) then

    do ele=1,ele_count(D(1))
       if (n_layers==1) then
          c=sqrt(g0*ele_val(D(k),ele))
          g=g0
       else
          ele_h=>ele_nodes(D(1),ele)
          if(k==1) then
             g=(rho(2)-rho(1))*g0/rho(k)
             c=sqrt((rho(2)-rho(1))/rho(1)*g0*&
                  D_old(ele_h,1))
          else
             g=(rho(2)-rho(1))*g0/rho(k)
             c=sqrt((rho(2)-rho(1))/rho(2)*g0*&
                  D_old(ele_h,2))
          end if
       end if
       call get_D_rhs(mesh,u(k),D(k),D_rhs,ele,dt,c=c,g=g0,dd=dd(k))
    end do

     ewrite(1,*) 'Solving for D'

    call petsc_solve(delta_D,&
         mesh%mass_h,D_rhs/rho(k),&
         startfromzero=.true.)

       D_old(:,k)=D(k)%val

       D(k)%val=D(k)%val+dt*delta_D

!       D(k)%val=10.0+dt*delta_D

!       d_rhs=D(k)%val
!       D(k)%val=0.5*(D(k)%val+D_old(:,k))+0.5*dt*delta_D
!       D_old(:,k)=D(k)%val
    end if
 end do
! call deallocate(d_mat)

end subroutine step_D

subroutine get_D_rhs(mesh,u,D,D_rhs,ele,dt,c,g,dd)
  type(dg_mesh) :: mesh
  type(vector_field) :: u
  type(scalar_field) :: D
  real, dimension(:) :: D_rhs
  integer ele
  real :: dt
  real, optional :: g,c(ele_loc(D,ele)),dd
!  type(csr_matrix) :: D_mat

  ! locals

  real, dimension(2,ele_loc(mesh%positions,ele)) :: X_ele
  real, dimension(2,face_loc(mesh%positions,1)) :: X_ele_f
  real, dimension(2,ele_ngi(u,ele)) :: ulocgi
   real, dimension(2,ele_ngi(D,ele)) :: grad_Dlocgi
  real, dimension(ele_ngi(D,ele)) :: hlocgi,detwei
  real, dimension(2,face_ngi(u,1)) :: ulocgi_f,ulocgi_f2,normal
  real, dimension(face_ngi(D,1)) :: hlocgi_f,detwei_f,nlocgi_f,wm,ws,q,c_f,hs
  real, dimension(ele_loc(D,ele),ele_ngi(D,ele),2) :: dh_t
  real, dimension(ele_loc(u,ele),ele_ngi(u,ele),2) :: du_t
  type(element_type), pointer :: shape_u,shape_h,shape_X,&
       shape_uf,shape_hf,shape_XF
  integer, dimension(:), pointer :: ele_h,ele_2,ele_mesh
  integer, dimension(face_loc(D,1)) :: ele_h_f
  integer :: i,j,face_mesh,face_D,face_u,face_u2
  real :: tau

  ele_h=>ele_nodes(D,ele)
  ele_2=>ele_neigh(D,ele)
  ele_mesh=>ele_neigh(mesh%positions,ele)

  shape_u=>ele_shape(u,ele)
  shape_h=>ele_shape(D,ele)
  shape_X=>ele_shape(mesh%positions, ele)

  hlocgi=1.0!ele_val_at_quad(D,ele)

  ulocgi=ele_val_at_quad(u,ele)

  X_ele=ele_val(mesh%positions,ele)

  call transform_to_physical(X_ele, shape_X,m=shape_h,dn_t=du_t,&
       dm_t=dh_t,detwei=detwei)

!  tau=1.0/sqrt( (2.0/dt)**2+sum(ulocgi*ulocgi))
  tau=0.0

  do i=1,2
     grad_Dlocgi(i,:)=matmul(ele_val(D,ele),dh_t(:,:,i))
  end do

!  call addto(D_mat,ele_h,ele_h,shape_shape(shape_h,shape_h,detwei))

!  call addto(D_mat,ele_h,ele_h,-dshape_dot_vector_shape(dh_t,ulocgi,&
!       shape_h,tau*detwei))

  D_rhs(ele_h)=D_rhs(ele_h)+dshape_dot_vector_rhs(dh_t,ulocgi,detwei)

!  D_rhs(ele_h)=D_rhs(ele_h)-shape_rhs(shape_h,detwei*hlocgi*ele_div_at_quad(u,ele,du_t))

!  D_rhs(ele_h)=D_rhs(ele_h)-dshape_dot_vector_rhs(dh_t,ulocgi,&
!       tau*detwei*hlocgi*ele_div_at_quad(u,ele,du_t))

!  D_rhs(ele_h)=D_rhs(ele_h)-dshape_dot_vector_rhs(dh_t,ulocgi,&
!       tau*detwei*sum(ulocgi*grad_Dlocgi,1))

  do j=1,size(ele_2)

     face_u=ele_face(u,ele,ele_2(j))
     face_D=ele_face(D,ele,ele_2(j))
     face_mesh=ele_face(mesh%positions,ele,ele_mesh(j))
     X_ele_f=face_val(mesh%positions,face_mesh)
     shape_Xf=>face_shape(mesh%positions,face_mesh)
     shape_uf=>face_shape(u,face_u)
     shape_hf=>face_shape(D,face_D)
     call transform_bdy_to_physical(X_ele,X_ele_f, shape_X, shape_Xf,&
          detwei_f=detwei_f,normal=normal)

     ulocgi_f=face_val_at_quad(u,face_u)
     hlocgi_f=face_val_at_quad(D,face_D)
     ele_h_f=face_global_nodes(D,face_D)
     c_f=matmul(c(face_local_nodes(D,face_D)),shape_hf%n)

     if(ele_2(j)>0) then

        if (ele_mesh(j)>0) then 
           assert(ele_mesh(j)==ele_2(j))
        end if

        face_u2=ele_face(u,ele_2(j),ele)
        ulocgi_f2=face_val_at_quad(u,face_u2)
        nlocgi_f=0.0

         hs=(1.0/g)*(c_f+0.25*sum((ulocgi_f-ulocgi_f2)*normal,1))**2.0
           q=sqrt(0.5*(hs+hlocgi_f)*hs/(hlocgi_f**2))
           where (q<1.0) q=1.0
           ws=sum(0.5*(ulocgi_f+ulocgi_f2)*normal,1)

           where (sum(ulocgi_f*normal,1)-c_f*q .ge. 0.0)&
              nlocgi_f=sum(ulocgi_f*normal,1)*hlocgi_f
           where (sum(ulocgi_f*normal,1)-c_f*q .lt. 0.0 .and.&
                ws .ge. 0.0 )
              wm=c_F*q&
                   /(c_f*q+sum(0.5*(ulocgi_f2-ulocgi_f)*normal,1))
              nlocgi_f=hlocgi_f*(sum(ulocgi_f*normal,1)+&
                   (sum(ulocgi_f*normal,1)-c_F*q)*(wm-1.0))
           end where
           where (sum(ulocgi_f2*normal,1)+c_f*q .ge. 0.0 .and.&
                ws .lt. 0.0 )
              wm=c_F*q&
                   /(c_f*q+sum(0.5*(ulocgi_f2-ulocgi_f)*normal,1))
              nlocgi_f=hlocgi_f*(sum(ulocgi_f2*normal,1)+&
                   (sum(ulocgi_f2*normal,1)+c_F*q)*(wm-1.0))
           end where
           where (sum(ulocgi_f2*normal,1)+c_f*q .lt. 0.0)&
                nlocgi_f=sum(ulocgi_f2*normal,1)*hlocgi_f

            nlocgi_f=sum(ulocgi_f*normal,1)

!           D_rhs(ele_h_f)=D_rhs(ele_h_f)-shape_rhs(shape_hf,detwei_f*nlocgi_f)

     end if
     

  end do

end subroutine get_D_rhs

subroutine get_u_rhs(mesh,u,D,D_old,m,bottom,rho,u_bar,u_rhs,u_step,k)
  type(dg_mesh) :: mesh
  type(vector_field) :: u,m(:)
  type(scalar_field) :: D(:),bottom
  real, dimension(:,:) :: u_rhs,u_bar,u_step
  real, dimension(:,:) :: D_old
  real :: rho(:),g
  integer :: k

! local 

  integer :: i,j,jj,kk,ele
  integer, pointer, dimension(:) :: ele_h
  real :: c(D(1)%mesh%shape%loc)

  u_rhs=0.0

  do ele=1,ele_count(u)

!     if (k==1) then
if (n_layers==1) then
   c=sqrt(g0*ele_val(D(k),ele))
   g=g0
else
   ele_h=>ele_nodes(D(1),ele)
   if(k==1) then
      c=sqrt((rho(2)-rho(1))/rho(1)*g0*&
           D_old(ele_h,1))
      g=(rho(2)-rho(1))*g0/rho(k)
   else
      c=sqrt((rho(2)-rho(1))/rho(2)*g0*&
           D_old(ele_h,2))
      g=(rho(2)-rho(1))*g0/rho(k)
   end if
end if
!        c=0.0
!     else
!        c=sqrt(rho(k)-rho(1))/rho(1)*g0*sum(&
!             (ele_val(bottom,ele)&
!             -sum(D_old(ele_nodes(bottom,ele),1:k-1),2))&
!             /ele_loc(bottom,ele))
!     end if
        call get_um_element_contribution(mesh,m,u,D,D_old,bottom,&
             u_rhs,ele,k)
!     call get_uu_element_contribution(mesh,u,D(k),D_old(:,k),rho(k),&
!          u_rhs,ele,c*0,"FLL",g=g,u_step=u_step)
  end do
        
end subroutine get_u_rhs

subroutine get_u_mat(mesh,u_s,u,dt,flux,c)

type(dg_mesh) :: mesh
real, dimension(:,:) :: u_s
type(vector_field) :: u
real :: dt, c
character (len=*) :: flux

integer :: ele

do ele=1,ele_count(u)
   call get_uu_mat_ele_contribution(mesh,u_s,u,ele,dt,flux=flux,c=c)

end do


end subroutine get_u_mat


subroutine  get_um_element_contribution(mesh,m,u,D,D_old,bottom,&
             u_rhs,ele,k)

  type(dg_mesh) :: mesh
  type(vector_field) :: m(:),u
  type(scalar_field) :: D(:),bottom
  real, dimension(:,:) :: u_rhs,D_old
  integer :: ele,k

  ! locals

  integer :: i,kk
  real, dimension(2,ele_loc(mesh%positions,ele)) :: X_ele
  real, dimension(2,face_loc(mesh%positions,1)) :: X_ele_f,X_ele_f2
  real, dimension(ele_ngi(D(1),ele)) :: dlocgi,detwei
  real, dimension(face_ngi(D(1),1)) :: dlocgi_f,detwei_F
  real, dimension(2,ele_ngi(m(1),ele)) :: mlocgi,mmlocgi
  real, dimension(2,face_ngi(m(1),1)) :: mlocgi_f,normal,mmlocgi_f
  real, dimension(2,ele_ngi(D(1),ele),0:n_layers) :: grad_hlocgi
  real, dimension(2,ele_ngi(D(1),ele)) :: grad_dlocgi
  real, dimension(2,2,ele_ngi(D(1),ele)) :: grad_mlocgi
  real, dimension(ele_loc(D(1),ele),n_layers) :: dloc
  real, dimension(2,ele_loc(m(1),ele),0:n_layers) :: mloc
  real, dimension(2,face_loc(m(1),1),n_layers) :: mloc_f
  integer, dimension(:), pointer :: ele_u,ele_2,ele_mesh
  real, dimension(ele_loc(D(k),ele),ele_ngi(u,ele),2) :: dh_t
  real, dimension(ele_loc(m(k),ele),ele_ngi(u,ele),2) :: dm_t
  real, dimension(ele_loc(u,ele),ele_ngi(u,ele),2) :: du_t
  type(element_type), pointer :: shape_u, shape_X,shape_h,shape_m
  type(element_type), pointer :: shape_uF, shape_Xf,shape_hf,shape_mf
  integer :: j,face_mesh,face_h,face_m,face_u
  integer :: ele_h_f(face_loc(D(k),1)),ele_u_f(face_loc(u,1))

  real :: dd(2)=(/50.0,350.0/),rho(2)=(/1020.0,1025.0/)

  shape_u=>ele_shape(u,ele)
  shape_m=>ele_shape(m(1),ele)
  shape_h=>ele_shape(D(1),ele)
  ele_u=>ele_nodes(u,ele)
  shape_X=>ele_shape(mesh%positions,ele)
  X_ele=ele_val(mesh%positions,ele)
  ele_2=>ele_neigh(m(k),ele)
  ele_mesh=>ele_neigh(mesh%positions,ele)


  call transform_to_physical(X_ele, n=shape_X,m=shape_m,&
       dn_t=du_t,dm_t=dm_t,detwei=detwei)
  dh_t=dm_t

  do kk=1,n_layers
     dloc(:,kk)=D_old(ele_nodes(D(1),ele),kk)
  end do

  do i=1,2
     do kk=1,n_layers
        grad_hlocgi(i,:,kk-1)=matmul(&!ele_val(bottom,ele)&
             -sum(dloc(:,kk:n_layers),2),dh_t(:,:,i))
     end do
     grad_hlocgi(i,:,n_layers)=0.0!matmul(ele_val(bottom,ele),dh_t(:,:,i))
     grad_dlocgi(i,:)=matmul(dloc(:,k),dh_t(:,:,i))
  end do


if (.not. LINEAR) then 
!  dlocgi=matmul(Dloc(:,k),shape_h%n)
!  mlocgi=ele_val_at_quad(m(k),ele)
  dlocgi=minval(D_old(:,k))
  mlocgi=rho(k)*minval(D_old(:,k))*g0
else
  dlocgi=dd(k)
  mlocgi=rho(k)*dd(k)*g0/2.0
end if

  do kk=1,n_layers
     if (.not. LINEAR) then
        mloc(:,:,kk)=ele_val(m(kk),ele)
     else
        mloc(:,:,kk)=rho(kk)*dd(kk)*g0/2.0+ele_val(m(kk),ele)
     end if
  end do


!        mloc(2,:,1)=(dd(1)+(rho(1)/rho(2))*&
!             (300.0+sqrt(300.0**2+4.0*rho(1)/rho(2)*350.0*50.0))&
!       /(2.0*rho(1)/rho(2)*350.0)*dd(2))*g0*dloc(:,1)
!        mloc(2,:,2)=(dd(1)+(rho(1)/rho(2))*&
!             (300.0-sqrt(300.0**2+4.0*rho(1)/rho(2)*350.0*50.0))&
!       /(2.0*rho(1)/rho(2)*350.0)*dd(2))*g0*dloc(:,2)
  do i=1,2
        grad_mlocgi(1,i,:)=matmul(sum(mloc(1,:,1:k-1),2),dm_t(:,:,i))
        grad_mlocgi(2,i,:)=matmul(sum(mloc(2,:,1:k),2),dm_t(:,:,i))
  end do

 ! mmlocgi(1,:)= matmul(sum(mloc(1,:,1:k-1),2),shape_m%n)
 ! mmlocgi(2,:)= matmul(sum(mloc(2,:,1:k-1),2),shape_m%n)
  

      u_rhs(:,ele_u)= u_rhs(:,ele_u)&
          +shape_vector_rhs(shape_u,grad_hlocgi(:,:,k),mlocgi(1,:)*detwei)
     u_rhs(:,ele_u)= u_rhs(:,ele_u)&
          +shape_vector_rhs(shape_u,grad_hlocgi(:,:,k-1),mlocgi(2,:)*detwei)
 !    u_rhs(:,ele_u)= u_rhs(:,ele_u)&
 !         -shape_vector_rhs(shape_u,grad_mlocgi(1,:,:),dlocgi*detwei)
 !    u_rhs(:,ele_u)= u_rhs(:,ele_u)&
 !         -shape_vector_rhs(shape_u,grad_mlocgi(2,:,:),dlocgi*detwei)
!     u_rhs(:,ele_u)= u_rhs(:,ele_u)&
!          +shape_vector_rhs(shape_u,grad_dlocgi,mmlocgi(1,:)*detwei)
!     u_rhs(:,ele_u)= u_rhs(:,ele_u)&
!          +shape_vector_rhs(shape_u,grad_dlocgi,mmlocgi(2,:)*detwei)
!     u_rhs(:,ele_u)= u_rhs(:,ele_u)&
!          +dshape_rhs(du_t,detwei*dlocgi*mmlocgi(1,:))
!     u_rhs(:,ele_u)= u_rhs(:,ele_u)&
!          +dshape_rhs(du_t,detwei*dlocgi*mmlocgi(2,:))

!dlocgi=2.0*matmul(Dloc(:,k),shape_h%n)

!     u_rhs(:,ele_u)= u_rhs(:,ele_u)&
!          +dshape_rhs(du_t,detwei*dlocgi*mlocgi(2,:))

  do kk=1,n_layers
     dloc(:,kk)=(D_old(ele_nodes(D(1),ele),kk)-ele_val(D(kk),ele))/2.0
  end do

  do i=1,2
     do kk=1,n_layers
        grad_hlocgi(i,:,kk-1)=matmul(ele_val(bottom,ele)&
             -sum(dloc(:,kk:n_layers),2),dh_t(:,:,i))
     end do
     grad_hlocgi(i,:,n_layers)=matmul(ele_val(bottom,ele),dh_t(:,:,i))
     grad_dlocgi(i,:)=matmul(dloc(:,k),dh_t(:,:,i))
  end do

!  u_rhs(:,ele_u)= u_rhs(:,ele_u)&
!          +shape_vector_rhs(shape_u,grad_dlocgi(:,:),g0*rho(k)&
!          *minval(D(k)%val)*detwei)


     do j=1,size(ele_2)
        face_mesh=ele_face(mesh%positions,ele,ele_mesh(j))
        X_ele_f=face_val(mesh%positions,face_mesh)
        shape_Xf=>face_shape(mesh%positions, face_mesh)
        call transform_bdy_to_physical(X_ele,X_ele_f, shape_X, shape_Xf,&
             detwei_f=detwei_f,normal=normal)
        face_h=ele_face(D(k),ele,ele_2(j))
        ele_h_f=face_global_nodes(D(k),face_h)
        shape_hf=>face_shape(D(k),face_h)
        Dlocgi_f=matmul(D_old(ele_h_f,k),shape_hf%n)
        face_u=ele_face(u,ele,ele_2(j))
        shape_uf=>face_shape(u,face_u)
        ele_u_f=face_global_nodes(u,face_u)
        face_m=ele_face(m(k),ele,ele_2(j))
        shape_mf=>face_shape(m(k),face_h)
!        mlocgi_f=face_val_at_quad(m(K),face_m)
        mlocgi_f(:,:)=g0
        do kk=1,n_layers
           mloc_f(:,:,kk)=face_val(m(kk),face_m)
        end do
!        mmlocgi_f(1,:)=matmul(sum(mloc_f(1,:,1:k-1),2),shape_mf%n)
!        mmlocgi_f(2,:)=matmul(sum(mloc_f(2,:,1:k-1),2),shape_mf%n)
        mmlocgi_f(2,:)=g0/2.0

!         u_rhs(:,ele_u_f)= u_rhs(:,ele_u_f)&
!          -shape_vector_rhs(shape_uf,normal,detwei_f*dlocgi_f*mlocgi_f(2,:))
!         u_rhs(:,ele_u_f)= u_rhs(:,ele_u_f)&
!          -shape_vector_rhs(shape_uf,normal,detwei_f*dlocgi_f*mmlocgi_f(1,:))
!         u_rhs(:,ele_u_f)= u_rhs(:,ele_u_f)&
!              +shape_vector_rhs(shape_uf,normal,detwei_f*dlocgi_f*mmlocgi_f(2,:))


      end do
  
end subroutine get_um_element_contribution

subroutine get_uu_mat_ele_contribution(mesh,u_s,u,ele,dt,c,flux,g)

  type(dg_mesh) :: mesh
  type(vector_field) :: u
  real, dimension(:,:) :: u_s
  real :: dt
  real, optional :: c,g
  integer :: ele
  character(len=*) :: flux

  ! locals

  real, dimension(2,ele_loc(mesh%positions,ele)) :: X_ele
  real, dimension(2,ele_ngi(u,ele)) :: ulocgi
  real, dimension(ele_ngi(u,ele)) :: hlocgi,detwei
  real, dimension(ele_loc(u,ele),ele_ngi(u,ele),2) :: du_t
  real, dimension(2,u%mesh%faces%shape%loc) :: X_ele_f
  real, dimension(2,u%mesh%faces%shape%ngi) :: ulocgi_f,ulocgi_f2,normal
  real, dimension(u%mesh%faces%shape%ngi) :: hlocgi_f,unlocgi_f,&
       unlocgi_f2,detwei_f,c_f,hllc_check
  real, dimension(2,2,u%mesh%faces%shape%ngi) :: unnlocgi_f,unnlocgi_f2
  integer, dimension(:), pointer :: ele_2,ele_mesh,ele_u,ele_h
  integer, dimension(u%mesh%faces%shape%loc) :: ele_u_f,ele_u_f2
  ! Shape functions.
  type(element_type), pointer :: shape_u,shape_h,shape_X,&
       shape_Xf, shape_uf,shape_hf
  real, dimension(u%mesh%faces%shape%ngi) :: q,hs,ws,wm
  real, dimension(2,2,ele_ngi(u,ele)) :: diag_mat
  real, dimension(2,2,u%mesh%faces%shape%ngi) :: diag_mat_f
  integer :: face_mesh,face_h,face_u,face_u2
  integer :: i,ii,j

  if (ele==1) print*, flux

  ele_u=>ele_nodes(u,ele)
!  ele_h=>ele_nodes(D_field,ele)
  shape_u=>ele_shape(u,ele)
  shape_X=>ele_shape(mesh%positions, ele)
  ele_2=>ele_neigh(mesh%connectivity,ele)
  ele_mesh=>ele_neigh(mesh%positions,ele)

  ! Locations of local vertices.
  X_ele=ele_val(mesh%positions, ele)
  
  ! Transform derivatives and weights into physical space.
  call transform_to_physical(X_ele, shape_X,m=shape_u,dm_t=du_t,detwei=detwei)
 
     ulocgi(1,:)=matmul(u_s(1,ele_u),shape_u%n)
     ulocgi(2,:)=matmul(u_s(2,ele_u),shape_u%n)
 
     call addto(mesh%u_mat,ele_u,ele_u,&
          shape_shape(shape_u,shape_u,detwei))

     call addto(mesh%u_mat,ele_u,ele_u,&
          -dshape_dot_vector_shape(du_t,ulocgi,shape_u,&
          dt*detwei))

  do j=1,size(ele_2)
     face_mesh=ele_face(mesh%positions,ele,ele_mesh(j))
     X_ele_f=face_val(mesh%positions,face_mesh)
     shape_Xf=>face_shape(mesh%positions, face_mesh)
     call transform_bdy_to_physical(X_ele,X_ele_f, shape_X, shape_Xf,&
          detwei_f=detwei_f,normal=normal)
     face_u=ele_face(u,ele,ele_2(j))
     shape_uf=>face_shape(u,face_u)
     ele_u_f=face_global_nodes(u,face_u)
     
     ulocgi_f(1,:)=matmul(u_s(1,ele_u_f),shape_uf%n)
     ulocgi_F(2,:)=matmul(u_s(2,ele_u_F),shape_uf%n)
    
        
     if  (ele_2(j) > 0) then

        face_u2=ele_face(u,ele_2(j),ele)
        ele_u_f2=face_global_nodes(u,face_u2)
        
        ulocgi_f2(1,:)=matmul(u_s(1,ele_u_f2),shape_uf%n)
        ulocgi_F2(2,:)=matmul(u_s(2,ele_u_F2),shape_uf%n)
     


        if( present(c)) then
           c_f=c!matmul(c(face_local_nodes(D_field,face_h)),shape_hf%n)
        else 
           c_f=0.0
        end if

        If (trim(flux)=="HLLC") then
            if (ele==1) print*, "HLLC flux"
        
           unnlocgi_f=0.0
           unnlocgi_f2=0.0

           hllc_check=0.0

           hs=(1.0/g)*(c_f+0.25*sum((ulocgi_f-ulocgi_f2)*normal,1))**2.0
           q=sqrt(0.5*(hs+hlocgi_f)*hs/(hlocgi_f**2))
           where (q<1.0) q=1.0
           ws=sum(0.5*(ulocgi_f+ulocgi_f2)*normal,1)

           where (sum(ulocgi_f*normal,1)-c_f*q .ge. 0.0)
              hllc_check=1.0
              unnlocgi_f(1,1,:)=sum(ulocgi_f*normal,1)
              unnlocgi_f(2,2,:)=sum(ulocgi_f*normal,1)
           end where
           where (sum(ulocgi_f*normal,1)-c_f*q .lt. 0.0 .and.&
                ws .ge. 0.0 )
               hllc_check=2.0
              wm=(sum(ulocgi_f*normal,1)-c_f*q)*c_F*q&
                   /(c_f*q+sum(0.5*(ulocgi_f2-ulocgi_f)*normal,1))

              unnlocgi_f(1,1,:)=c_f*q+0.5*wm+0.5*wm*(normal(2,:)*normal(2,:))
              unnlocgi_f(1,2,:)=-0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f(2,1,:)=-0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f(2,2,:)=c_f*q+0.5*wm+0.5*wm*(normal(1,:)*normal(1,:))
              unnlocgi_f2(1,1,:)=0.5*wm*(normal(1,:)*normal(1,:))
              unnlocgi_f2(1,2,:)=0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f2(2,1,:)=0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f2(2,2,:)=0.5*wm*(normal(2,:)*normal(2,:))
         end where
         where (sum(ulocgi_f2*normal,1)+c_f*q .ge. 0.0 .and.&
                ws .le. 0.0 )

             hllc_check=3.0
            wm=(sum(ulocgi_f2*normal,1)+c_f*q)*c_F*q&
                   /(c_f*q+sum(0.5*(ulocgi_f2-ulocgi_f)*normal,1))

              unnlocgi_f(1,1,:)=0.5*wm*(normal(1,:)*normal(1,:))
              unnlocgi_f(1,2,:)=0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f(2,1,:)=0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f(2,2,:)=0.5*wm*(normal(2,:)*normal(2,:))
              unnlocgi_f2(1,1,:)=-c_f*q+0.5*wm+0.5*wm*(normal(2,:)*normal(2,:))
              unnlocgi_f2(1,2,:)=-0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f2(2,1,:)=-0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f2(2,2,:)=-c_f*q+0.5*wm+0.5*wm*(normal(1,:)*normal(1,:))
           
         end where
         where (sum(ulocgi_f2*normal,1)+c_f*q .lt. 0.0)
            
            hllc_check=4.0
              unnlocgi_f2(1,1,:)=sum(ulocgi_f2*normal,1)
              unnlocgi_f2(2,2,:)=sum(ulocgi_f2*normal,1)
        end where        

        assert(all(hllc_check>0.0))

!              u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f) -&
!                   shape_rhs(shape_uf,&
!                   rho*hlocgi_f*sum(unnlocgi_f(i,:,:)*ulocgi_f,1)*detwei_f)
!              u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f) -&
!                   shape_rhs(shape_uf,&
!                   rho*hlocgi_f*sum(unnlocgi_f2(i,:,:)*ulocgi_f2,1)*detwei_f)

     end if

        If (trim(flux)=="FLL") then
           if (ele==1) print*, "FLL flux"
        
        unlocgi_f=0.5*(sum(ulocgi_f*normal,1)&
             +(alph_up*abs(sum(0.5*(ulocgi_f+ulocgi_f2)*normal,1))))
        unlocgi_f2=0.5*(sum(ulocgi_f2*normal,1)&
             -(alph_up*abs(sum(0.5*(ulocgi_f+ulocgi_f2)*normal,1))))
        
      
          call addto(mesh%u_mat,ele_u_f,ele_u_f,&
               shape_shape(shape_uf,shape_uf,dt*unlocgi_f*detwei_f))
          call addto(mesh%u_mat,ele_u_f,ele_u_f2,&
               shape_shape(shape_uf,shape_uf,dt*unlocgi_f2*detwei_f))
          if (present(c)) then
             call addto(mesh%u_mat,ele_u_f,ele_u_f,&
             shape_shape(shape_uf,shape_uf,dt*c_f*detwei_f))
            call addto(mesh%u_mat,ele_u_f,ele_u_f2,&
               -shape_shape(shape_uf,shape_uf,dt*c_f*detwei_f))
           end if  
     end if

     If (trim(flux)=="Mean") then
           if (ele==1) print*, "Mean flux"
        
           do i=1,2

              unlocgi_f=0.25*(sum((ulocgi_f+ulocgi_f2)*normal,1))
              unlocgi_f2=0.25*(sum((ulocgi_f2+ulocgi_f)*normal,1))
             
!           u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f) -&
!                shape_rhs(shape_uf,&
!                rho*hlocgi_f*unlocgi_f*ulocgi_f(i,:)*detwei_f)
!           u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f) -&
!                shape_rhs(shape_uf,&
!                rho*hlocgi_f*unlocgi_f2*ulocgi_f2(i,:)*detwei_f)
!           if (present(c)) then
!              u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f)-&
!                   shape_rhs(shape_uf,&
!                   c_f*rho*hlocgi_f*detwei_f*ulocgi_f(i,:))
!              u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f)+&
!                   (shape_rhs(shape_uf,&
!                   c_f*rho*hlocgi_f*detwei_f*ulocgi_f2(i,:)))
!           end if
           end do

     end if

     If (trim(flux)=="None") then
           if (ele==1) print*, "Interior flux"
        
           do i=1,2

              unlocgi_f=(sum(ulocgi_f*normal,1))
              unlocgi_f2=0.0
             
!           u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f) -&
!                shape_rhs(shape_uf,&
!                rho*hlocgi_f*unlocgi_f*ulocgi_f(i,:)*detwei_f)
!           u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f) -&
!                shape_rhs(shape_uf,&
!                rho*hlocgi_f*unlocgi_f2*ulocgi_f2(i,:)*detwei_f)
!           if (present(c)) then
!              u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f)-&
!                   shape_rhs(shape_uf,&
!                   c_f*rho*hlocgi_f*detwei_f*ulocgi_f(i,:))
!              u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f)+&
!                   (shape_rhs(shape_uf,&
!                   c_f*rho*hlocgi_f*detwei_f*ulocgi_f2(i,:)))
!           end if
        end do

     end If
           
  end if
end do

end subroutine get_uu_mat_ele_contribution

subroutine get_uu_element_contribution(mesh,u,D_field,D,rho,u_rhs,ele,c,flux,&
     u1,u2,diag,g, u_step)

  type(dg_mesh) :: mesh
  type(vector_field) :: u
  type(scalar_field) :: D_field
  real, dimension(:) :: D
  real, dimension(:,:) :: u_rhs
  real, dimension(:,:), optional :: u_step
  real :: rho
  real, optional :: c(:),u1(:),u2(:),g
  integer :: ele
  logical, optional :: diag
  character(len=*) :: flux

  ! locals

  real, dimension(2,ele_loc(mesh%positions,ele)) :: X_ele
  real, dimension(2,ele_ngi(u,ele)) :: ulocgi
  real, dimension(ele_ngi(D_field,ele)) :: hlocgi,detwei
  real, dimension(ele_loc(u,ele),ele_ngi(u,ele),2) :: du_t
  real, dimension(2,u%mesh%faces%shape%loc) :: X_ele_f
  real, dimension(2,u%mesh%faces%shape%ngi) :: ulocgi_f,ulocgi_f2,normal
  real, dimension(u%mesh%faces%shape%ngi) :: hlocgi_f,unlocgi_f,&
       unlocgi_f2,detwei_f,c_f,hllc_check
  real, dimension(2,2,u%mesh%faces%shape%ngi) :: unnlocgi_f,unnlocgi_f2
  integer, dimension(:), pointer :: ele_2,ele_mesh,ele_u,ele_h
  integer, dimension(u%mesh%faces%shape%loc) :: ele_u_f,ele_u_f2
  integer, dimension(D_field%mesh%faces%shape%loc) :: ele_h_f
  ! Shape functions.
  type(element_type), pointer :: shape_u,shape_h,shape_X,&
       shape_Xf, shape_uf,shape_hf
  real, dimension(u%mesh%faces%shape%ngi) :: q,hs,ws,wm
  real, dimension(2,2,ele_ngi(u,ele)) :: diag_mat
  real, dimension(2,2,u%mesh%faces%shape%ngi) :: diag_mat_f
  integer :: face_mesh,face_h,face_u,face_u2
  integer :: i,ii,j

  if (ele==1) print*, flux

  ele_u=>ele_nodes(u,ele)
  ele_h=>ele_nodes(D_field,ele)
  shape_u=>ele_shape(u,ele)
  shape_h=>ele_shape(D_field,ele)
  shape_X=>ele_shape(mesh%positions, ele)
  ele_2=>ele_neigh(mesh%connectivity,ele)
  ele_mesh=>ele_neigh(mesh%positions,ele)

  ! Locations of local vertices.
  X_ele=ele_val(mesh%positions, ele)
  
  ! Transform derivatives and weights into physical space.
  call transform_to_physical(X_ele, shape_X,m=shape_u,dm_t=du_t,detwei=detwei)
  
  hlocgi=matmul(D(ele_h),shape_h%n)

  if (present(diag)) then
     diag_mat(1,1,:)=0.0
     diag_mat(2,1,:)=1.0
     diag_mat(2,2,:)=0.0
     diag_mat(1,2,:)=1.0
     diag_mat_f(1,1,:)=0.0
     diag_mat_f(2,1,:)=1.0
     diag_mat_f(2,2,:)=0.0
     diag_mat_f(1,2,:)=1.0
  else
     diag_mat=1.0
     diag_mat_f=1.0
  end if
     

  if (present(u1) ) then
     ulocgi(1,:)=matmul(u1(ele_u),shape_u%n)
     ulocgi(2,:)=matmul(u2(ele_u),shape_u%n)
  else if (present(u_step) ) then
     ulocgi(1,:)=matmul(u_step(1,ele_u),shape_u%n)
     ulocgi(2,:)=matmul(u_step(2,ele_u),shape_u%n)
  else
     ulocgi=ele_val_at_quad(u,ele)
  end if

  do i=1,2
     u_rhs(i,ele_u)=u_rhs(i,ele_u)+dshape_dot_vector_rhs(du_t,ulocgi,&
          rho*hlocgi*ulocgi(i,:)*detwei)
  end do

  do j=1,size(ele_2)
     face_mesh=ele_face(mesh%positions,ele,ele_mesh(j))
     X_ele_f=face_val(mesh%positions,face_mesh)
     shape_Xf=>face_shape(mesh%positions, face_mesh)
     call transform_bdy_to_physical(X_ele,X_ele_f, shape_X, shape_Xf,&
          detwei_f=detwei_f,normal=normal)
     face_h=ele_face(D_field,ele,ele_2(j))
     ele_h_f=face_global_nodes(D_field,face_h)
     shape_hf=>face_shape(D_field,face_h)
     hlocgi_f=matmul(D(ele_h_f),shape_hf%n)
     face_u=ele_face(u,ele,ele_2(j))
     shape_uf=>face_shape(u,face_u)
     ele_u_f=face_global_nodes(u,face_u)
     if (present(u1) ) then
        ulocgi_f(1,:)=matmul(u1(ele_u_f),shape_uf%n)
        ulocgi_F(2,:)=matmul(u2(ele_u_F),shape_uf%n)
     else if (present(u_step) ) then
        ulocgi_f(1,:)=matmul(u_step(1,ele_u_f),shape_uf%n)
        ulocgi_F(2,:)=matmul(u_step(2,ele_u_F),shape_uf%n)
     else
        ulocgi_f=face_val_at_quad(u,face_u)
     end if
        
     if  (ele_2(j) > 0) then

        face_u2=ele_face(u,ele_2(j),ele)
        ele_u_f2=face_global_nodes(u,face_u2)
        if (present(u1) ) then
           ulocgi_f2(1,:)=matmul(u1(ele_u_f2),shape_uf%n)
           ulocgi_F2(2,:)=matmul(u2(ele_u_F2),shape_uf%n)
        else if (present(u_step) ) then
           ulocgi_f2(1,:)=matmul(u_step(1,ele_u_f2),shape_uf%n)
           ulocgi_F2(2,:)=matmul(u_step(2,ele_u_F2),shape_uf%n)
        else
           ulocgi_f2=face_val_at_quad(u,face_u2)
        end if


        if( present(c)) then
           c_f=matmul(c(face_local_nodes(D_field,face_h)),shape_hf%n)
        else 
           c_f=0.0
        end if

        If (trim(flux)=="HLLC") then
            if (ele==1) print*, "HLLC flux"
        
           unnlocgi_f=0.0
           unnlocgi_f2=0.0

           hllc_check=0.0

           hs=(1.0/g)*(c_f+0.25*sum((ulocgi_f-ulocgi_f2)*normal,1))**2.0
           q=sqrt(0.5*(hs+hlocgi_f)*hs/(hlocgi_f**2))
           where (q<1.0) q=1.0
           ws=sum(0.5*(ulocgi_f+ulocgi_f2)*normal,1)

           where (sum(ulocgi_f*normal,1)-c_f*q .ge. 0.0)
              hllc_check=1.0
              unnlocgi_f(1,1,:)=sum(ulocgi_f*normal,1)
              unnlocgi_f(2,2,:)=sum(ulocgi_f*normal,1)
           end where
           where (sum(ulocgi_f*normal,1)-c_f*q .lt. 0.0 .and.&
                ws .ge. 0.0 )
               hllc_check=2.0
              wm=(sum(ulocgi_f*normal,1)-c_f*q)*c_F*q&
                   /(c_f*q+sum(0.5*(ulocgi_f2-ulocgi_f)*normal,1))

              unnlocgi_f(1,1,:)=c_f*q+0.5*wm+0.5*wm*(normal(2,:)*normal(2,:))
              unnlocgi_f(1,2,:)=-0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f(2,1,:)=-0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f(2,2,:)=c_f*q+0.5*wm+0.5*wm*(normal(1,:)*normal(1,:))
              unnlocgi_f2(1,1,:)=0.5*wm*(normal(1,:)*normal(1,:))
              unnlocgi_f2(1,2,:)=0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f2(2,1,:)=0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f2(2,2,:)=0.5*wm*(normal(2,:)*normal(2,:))
         end where
         where (sum(ulocgi_f2*normal,1)+c_f*q .ge. 0.0 .and.&
                ws .le. 0.0 )

             hllc_check=3.0
            wm=(sum(ulocgi_f2*normal,1)+c_f*q)*c_F*q&
                   /(c_f*q+sum(0.5*(ulocgi_f2-ulocgi_f)*normal,1))

              unnlocgi_f(1,1,:)=0.5*wm*(normal(1,:)*normal(1,:))
              unnlocgi_f(1,2,:)=0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f(2,1,:)=0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f(2,2,:)=0.5*wm*(normal(2,:)*normal(2,:))
              unnlocgi_f2(1,1,:)=-c_f*q+0.5*wm+0.5*wm*(normal(2,:)*normal(2,:))
              unnlocgi_f2(1,2,:)=-0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f2(2,1,:)=-0.5*wm*(normal(1,:)*normal(2,:))
              unnlocgi_f2(2,2,:)=-c_f*q+0.5*wm+0.5*wm*(normal(1,:)*normal(1,:))
           
         end where
         where (sum(ulocgi_f2*normal,1)+c_f*q .lt. 0.0)
            
            hllc_check=4.0
              unnlocgi_f2(1,1,:)=sum(ulocgi_f2*normal,1)
              unnlocgi_f2(2,2,:)=sum(ulocgi_f2*normal,1)
        end where        

        assert(all(hllc_check>0.0))

        do i=1,2
              u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f) -&
                   shape_rhs(shape_uf,&
                   rho*hlocgi_f*sum(unnlocgi_f(i,:,:)*ulocgi_f,1)*detwei_f)
              u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f) -&
                   shape_rhs(shape_uf,&
                   rho*hlocgi_f*sum(unnlocgi_f2(i,:,:)*ulocgi_f2,1)*detwei_f)
        end do

     end if

        If (trim(flux)=="FLL") then
           if (ele==1) print*, "FLL flux"
        
        unlocgi_f=0.5*(sum(ulocgi_f*normal,1)&
             +(alph_up*abs(sum(0.5*(ulocgi_f+ulocgi_f2)*normal,1))))
        unlocgi_f2=0.5*(sum(ulocgi_f2*normal,1)&
             -(alph_up*abs(sum(0.5*(ulocgi_f+ulocgi_f2)*normal,1))))
        
        do i=1,2
           u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f) -&
                shape_rhs(shape_uf,&
                rho*hlocgi_f*unlocgi_f*ulocgi_f(i,:)*detwei_f)
           u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f) -&
                shape_rhs(shape_uf,&
                rho*hlocgi_f*unlocgi_f2*ulocgi_f2(i,:)*detwei_f)
           if (present(c)) then
              u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f)-&
                   shape_rhs(shape_uf,&
                   c_f*rho*hlocgi_f*detwei_f*ulocgi_f(i,:))
              u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f)+&
                   (shape_rhs(shape_uf,&
                   c_f*rho*hlocgi_f*detwei_f*ulocgi_f2(i,:)))
           end if
        end do

     end if

     If (trim(flux)=="Mean") then
           if (ele==1) print*, "Mean flux"
        
           do i=1,2

              unlocgi_f=0.25*(sum((ulocgi_f+ulocgi_f2)*normal,1))
              unlocgi_f2=0.25*(sum((ulocgi_f2+ulocgi_f)*normal,1))
             
           u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f) -&
                shape_rhs(shape_uf,&
                rho*hlocgi_f*unlocgi_f*ulocgi_f(i,:)*detwei_f)
           u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f) -&
                shape_rhs(shape_uf,&
                rho*hlocgi_f*unlocgi_f2*ulocgi_f2(i,:)*detwei_f)
           if (present(c)) then
              u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f)-&
                   shape_rhs(shape_uf,&
                   c_f*rho*hlocgi_f*detwei_f*ulocgi_f(i,:))
              u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f)+&
                   (shape_rhs(shape_uf,&
                   c_f*rho*hlocgi_f*detwei_f*ulocgi_f2(i,:)))
           end if
        end do

     end if

     If (trim(flux)=="None") then
           if (ele==1) print*, "Interior flux"
        
           do i=1,2

              unlocgi_f=(sum(ulocgi_f*normal,1))
              unlocgi_f2=0.0
             
           u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f) -&
                shape_rhs(shape_uf,&
                rho*hlocgi_f*unlocgi_f*ulocgi_f(i,:)*detwei_f)
           u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f) -&
                shape_rhs(shape_uf,&
                rho*hlocgi_f*unlocgi_f2*ulocgi_f2(i,:)*detwei_f)
           if (present(c)) then
              u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f)-&
                   shape_rhs(shape_uf,&
                   c_f*rho*hlocgi_f*detwei_f*ulocgi_f(i,:))
              u_rhs(i,ele_u_f) = u_rhs(i,ele_u_f)+&
                   (shape_rhs(shape_uf,&
                   c_f*rho*hlocgi_f*detwei_f*ulocgi_f2(i,:)))
           end if
        end do

     end If
           
  end if
end do

end subroutine get_uu_element_contribution
  
subroutine get_barotropic_mode(mesh,D,u,H,u_bar,rho)
  type(dg_mesh) :: mesh
  type(scalar_field) :: D(:), H
  type(vector_field) :: u(:), u_bar
  real :: rho(:)

  ! locals


  integer :: i,k,ele
  real, allocatable :: u_b(:,:)
  real :: Dloc(D(1)%mesh%shape%loc)

  allocate(u_b(2,node_count(u(1))))

  H%val=0.0
  u_bar%val(1)%ptr=0.0
  u_bar%val(2)%ptr=0.0
  u_b=0.0

  do k=1,n_layers
     H%val=H%val+rho(k)*D(k)%val/rho(1)
  end do

!  call get_u_bar(mesh,u,D,rho,H,u_b)

  do ele=1,ele_count(u(1))
     do k=1,n_layers
        Dloc=ele_val(D(k),ele)/ele_val(H,ele)
        do i=1,2
        u_b(i,ele_nodes(u(1),ele))=u_b(i,ele_nodes(u(1),ele))&
             +rho(k)*Dloc((/1,3,6/))*ele_val(u(k),ele,i)/rho(1)
        end do
     end do
  end do

  do i=1,2
     u_bar%val(i)%ptr=u_b(i,:)
  end do

  print*, 'Barotropic mode, max_D=', maxval(H%val)
  print*, 'Barotropic mode, max_u_1=', maxval(abs(u_bar%val(1)%ptr))
  print*, 'Barotropic mode, max_u_2=', maxval(abs(u_bar%val(2)%ptr))

  if (n_layers ==2 ) then 
     print*, 'sanity_check=',( rho(1)*maxval(D(1)%val)*maxval(u(1)%val(1)%ptr)&
          +rho(2)*minval(D(2)%val)*minval(u(2)%val(1)%ptr))&
          /(rho(1)*maxval(H%val))
  end if

end subroutine get_barotropic_mode

subroutine get_u_bar(mesh,u,D,rho,H,u_b)

  type(dg_mesh) :: mesh
  type(scalar_field) :: D(:), H
  type(vector_field) :: u(:)
  real :: rho(:),u_b(:,:)

  ! locals
  integer :: i,ele
  real, allocatable :: u_bar_rhs(:,:)

  allocate(u_bar_rhs(2,node_count(u(1))))

  u_bar_rhs=0.0
  call get_u_bar_rhs(mesh,u,D,rho,H,u_bar_rhs)

  call get_weighted_dg_inverse_mass_matrix(&
       mesh%M_inv%blocks(1,1),&
       u(1)%mesh,mesh%positions,H,rho(1))

  do ele=1,node_count(u(1))
        do i=1,2
        u_b(i,ele)=sum(&
             u_bar_rhs(i,row_m(mesh%m_inv%blocks(1,1),ele))*&
             row_val_ptr(mesh%m_inv%blocks(1,1),ele))
     end do
  end do

end subroutine get_u_bar

subroutine get_u_bar_rhs(mesh,u,D,rho,H,u_bar_rhs)
  type(dg_mesh) :: mesh
  type(scalar_field) :: D(:),H
  type(vector_field) :: u(:)
  real :: rho(:),u_bar_rhs(:,:)

  ! locals

  integer :: k, ele
  real, dimension(mesh%positions%mesh%shape%ngi) :: detwei
  real, dimension(2,mesh%positions%mesh%shape%loc) :: X_ele
  integer, dimension(:), pointer :: ele_u
  type(element_type), pointer :: shape_u,shape_X


  do ele=1,ele_count(u(1))

     X_ele=ele_val(mesh%positions,ele)
     shape_x=>ele_shape(mesh%positions,ele)
     shape_u=>ele_shape(u(1),ele)
     ele_u=>ele_nodes(u(1),ele)

     call transform_to_physical(X_ele, shape_X,detwei=detwei)

     do k=1,n_layers
        u_bar_rhs(:,ele_u)=u_bar_rhs(:,ele_u)+&
             shape_vector_rhs(shape_u,ele_val_at_quad(u(k),ele),&
             rho(k)*ele_val_at_quad(D(k),ele)*detwei)
     end do

  end do

end subroutine get_u_bar_rhs
             

subroutine barotropic_H_step(mesh,H,H_old,u_bar,dt,step)
  type(dg_mesh) :: mesh
  type(scalar_field) :: H
  type(vector_field) :: u_bar
  real :: dt, H_old(:)
  integer :: step

  ! locals
  
  real, allocatable :: delta_H0(:),H_rhs(:)
  real :: c(ele_loc(H,1))
  integer :: k,ele
!  type(csr_matrix) :: h_mat


!  h_mat=clone(mesh%mass_h)
  allocate(delta_H0(node_count(H)),H_rhs(node_count(H)))

  delta_H0=0.0
  H_rhs=0.0
!  call zero(h_mat)

  do ele=1,ele_count(H)
      c=sqrt(g0*ele_val(H,ele))
     call get_D_rhs(mesh,u_bar,H,h_rhs,ele,dt,c=c,g=g0)
  end do

  call petsc_solve(delta_H0,&
         mesh%mass_H,H_rhs,&
         startfromzero=.true.)

  if (.true.) then
  
  if (step==1) then
     H_old=H%val
     H%val=H%val+dt*delta_H0
  else if (step==2) then
     H_rhs=H%val
     H%val=0.5*(H%val+H_old)+0.5*dt*delta_H0
     H_old=H_rhs
  end if
  
  end if

!  call deallocate(h_mat)
  
  end subroutine barotropic_H_step

subroutine get_baro_rhs(mesh,h_mat,H,u_bar,h_rhs,ele)
  type(dg_mesh) :: mesh
  type(vector_field) :: u_bar
  type(scalar_field) :: H
  real, dimension(:) :: H_rhs
  real :: dt
  integer ele
  type(csr_matrix) :: h_mat

  ! locals

  real, dimension(2,ele_loc(mesh%positions,ele)) :: X_ele
  real, dimension(2,ele_ngi(u_bar,ele)) :: Mlocgi, grad_hlocgi,grad_dlocgi
  real, dimension(ele_ngi(H,ele)) :: hlocgi,detwei,blocgi,div_mlocgi
  real, dimension(ele_loc(H,ele),ele_ngi(H,ele),2) :: dh_t
  real, dimension(ele_loc(u_bar,ele),ele_ngi(u_bar,ele),2) :: du_t
  type(element_type), pointer :: shape_u,shape_h,shape_X
  integer, dimension(:), pointer :: ele_h
  integer :: i,k

  ele_h=>ele_nodes(H,ele)
  shape_u=>ele_shape(u_bar,ele)
  shape_h=>ele_shape(H,ele)
  shape_X=>ele_shape(mesh%positions, ele)

  hlocgi=ele_val_at_quad(H,ele)
  mlocgi=ele_val_at_quad(u_bar,ele)

  X_ele=ele_val(mesh%positions,ele)

  call transform_to_physical(X_ele, shape_X,m=shape_h,dn_t=du_t,&
       dm_t=dh_t,detwei=detwei)

  call addto(h_mat,ele_h,ele_h,shape_shape(shape_h,shape_h,detwei))

  H_rhs(ele_h)=H_rhs(ele_h)+dshape_dot_vector_rhs(dh_t,mlocgi,&
       detwei*hlocgi)

end subroutine get_baro_rhs

subroutine d2c(dcsr,csr)
    !!< Given a dcsr matrix return a csr matrix. The dcsr matrix is left
    !!< untouched. 
    type(dynamic_csr_matrix), intent(in) :: dcsr
    type(csr_matrix) :: csr

    integer :: rows, columns, nentries, i, rowptr, rowlen
    integer, dimension(1) :: rowpos
    
    rows=size(dcsr,1)
    columns=size(dcsr,2)
    nentries=entries(dcsr)

!    call allocate(csr, rows, columns, nentries)

    rowptr=1
    do i=1,rows
       csr%sparsity%findrm(i)=rowptr
       
       rowlen=size(dcsr%colm(i)%ptr)

       csr%sparsity%colm(rowptr:rowptr+rowlen-1)=dcsr%colm(i)%ptr

       csr%val(rowptr:rowptr+rowlen-1)=dcsr%val(i)%ptr

       if (any(dcsr%colm(i)%ptr==i)) then

          rowpos=minloc(dcsr%colm(i)%ptr,mask=dcsr%colm(i)%ptr==i)

          csr%sparsity%centrm(i)=rowptr+rowpos(1)-1

       else
!          ewrite(1,*) "Missing diagonal element"
          
          csr%sparsity%centrm(i)=0
       end if
       
       rowptr = rowptr + rowlen

    end do

    csr%sparsity%findrm(rows+1) = rowptr
    
  end subroutine d2c

  subroutine dcsr_matmul_T(product,matrix1, matrix2, model,check)
    !!< Perform the matrix multiplication:
    !!< 
    !!<     matrix1*matrix2^T
    !!< 
    type(dynamic_csr_matrix), intent(in) :: matrix1, matrix2
    type(dynamic_csr_matrix) :: product
    type(dynamic_csr_matrix), intent(in), optional :: model
    logical, intent(in), optional :: check

    type(integer_vector), dimension(:), allocatable :: hitlist
    integer, dimension(:), allocatable :: size_hitlist

    integer, dimension(:), pointer :: row, col
    integer :: i,j,k1,k2,jrow,jcol
    real :: entry
    logical :: addflag
    logical :: l_check
    real , allocatable, dimension(:) :: vec, m2Tvec, m1m2tvec, productvec

    l_check = .false.
    if(present(check)) then
       L_check = check
    end if

    assert(size(matrix1,2)==size(matrix2,2))

    call zero(product)

    do i=1, size(matrix1,1)
       row=>matrix1%colm(i)%ptr

       if(size(row)>0) then
          do jcol = 1, size(product%colm(i)%ptr)
             j = product%colm(i)%ptr(jcol)
             col=>matrix2%colm(j)%ptr

             if(size(col)>0) then
                entry=0.0

                addflag = .false.

                k1 = 1
                k2 = 1
                do
                   if((k1.gt.size(row)).or.(k2.gt.size(col))) exit
                   if(row(k1)<col(k2)) then
                      k1 = k1 + 1
                   else
                      if(row(k1)==col(k2)) then
                         ! Note the transpose in the second val call.
                         entry=entry+val(matrix1,i,row(k1))* &
                              val(matrix2,j,row(k1))
                         addflag = .true.
                         k1 = k1 + 1
                         k2 = k2 + 1
                      else
                         k2 = k2 + 1
                      end if
                   end if
                end do
                if(addflag) then
                   call addto(product,i,j,entry)
                end if
             end if
          end do
       end if
    end do

  end subroutine dcsr_matmul_T

 subroutine py_plot_write_state(filename,positions,u,&
      v,u_bar,h,h_0,bottom,time,mod_time) 
    character(len=*), intent(in) :: filename
    type(vector_field), intent(in) :: positions,u(:),v(:),u_bar
    type(scalar_field), intent(in) :: h(:),bottom,h_0
    integer, intent(in), optional :: time
    real, intent(in), optional :: mod_time
    
    !locals
    integer :: io1, ele, nod, i,j
    integer, dimension(:), pointer :: uele
    integer, dimension(:), pointer :: h_ele
    real, dimension(:,:), pointer :: ptsD1, ptsD2
    real, dimension(:,:), pointer :: ptsU, ptsV
    integer, dimension(3) :: vertices
    real, dimension(6) :: val
    character(4) :: tc
    character(10) :: model_time
    character(128) :: file_title

    vertices(1) = 1
    vertices(2) = 3
    vertices(3) = 6
    
    if (present(time)) then
       write(tc,'(i4.4)') time
    else
       tc=""
    end if
    if (present(mod_time))   write(model_time,'(f10.3)') mod_time

    file_title=""
    if (.not. SHALLOW) then
       file_title= "Green Nagdhi Equation"
    else
       file_title= "Shallow Water Equation"
    end if

    allocate( ptsD1(2,3*ele_count(positions)),ptsD2(2,6*ele_count(positions)))



    do ele = 1, element_count(u(1))
       ptsD1(:,1+3*(ele-1):3*ele) = ele_val(positions,ele)       
       ptsD2(:,(/1,3,6/)+6*(ele-1)) = ele_val(positions,ele)       
       ptsD2(:,(/2,4,5/)+6*(ele-1)) = &
            0.5*(ptsD2(:,(/1,6,6/)+6*(ele-1)) +&
            ptsD2(:,(/3,1,3/)+6*(ele-1)))
    end do

    open(unit=2502, file="./"//trim(trim(filename)//tc//'_m.dat'), &
         iostat=io1)

    ewrite(0,*) 'writing data'


    if(io1.ne.0) then
       write (0,*) 'Could not open file '
       stop
    else       

       write(unit=2502, iostat=io1, fmt=*) u(1)%mesh%shape%loc, 2*n_layers,&
            ele_count(u(1)), filename, " Elliptic ", alpha

       write(unit=2502, iostat=io1, fmt=*) trim(file_title)//" eliptic, time="&
            //trim(model_time) 

       do j=1,n_layers
          do i=1,2
             do ele = 1, ele_count(u(1))
                if (size(ele_val(u(1),ele,i))==3) then
                   h_ele=>ele_nodes(positions,ele)
                   ptsU=>ptsD1(:,1+3*(ele-1):3*ele)
                else
                   h_ele=>ele_nodes(h(1),ele)
                   ptsU=>ptsD2(:,1+6*(ele-1):6*ele)
                end if
                val(1:size(ele_val(u(1),ele,i)))=ele_val(u(j),ele,i)
                write(unit=2502, iostat=io1, fmt=*)&
                     ptsU(1,node_list(size(h_ele)))
                write(unit=2502, iostat=io1, fmt=*)&
                     ptsU(2,node_list(size(h_ele)))
                write(unit=2502, iostat=io1, fmt=*)&
                     val(node_list(size(h_ele)))

          if(io1 /= 0) then
             write(0,*) 'could not write to file ' &
                  //trim(trim(filename)//'_meshpts.dat'//tc)
             stop
          end if
       end do
    end do
    end do
       close(unit=2502)
    end if

 open(unit=2502, file="./"//trim(trim(filename)//tc//'_u.dat'), &
         iostat=io1)

    ewrite(0,*) 'writing data'


    if(io1.ne.0) then
       write (0,*) 'Could not open file '
       stop
    else       

       write(unit=2502, iostat=io1, fmt=*) v(1)%mesh%shape%loc, 2*(n_layers+1),&
            ele_count(v(1)), filename, " Velocity ", alpha
       
       write(unit=2502, iostat=io1, fmt=*) trim(file_title)//" velocity, time="&
            //trim(model_time) 


       do j=1,n_layers+1
       do i=1,2
       do ele = 1, ele_count(u(1))
          if (size(ele_val(v(1),ele,i))==3) then
             h_ele=>ele_nodes(positions,ele)
             ptsU=>ptsD1(:,1+3*(ele-1):3*ele)
          else
             h_ele=>ele_nodes(h(1),ele)
             ptsU=>ptsD2(:,1+6*(ele-1):6*ele-1)
          end if
          if(j<n_layers+1) then
             val(1:size(ele_val(v(1),ele,i)))=ele_val(v(j),ele,i)
          else
             val(1:size(ele_val(v(1),ele,i)))=ele_val(u_bar,ele,i)
          end if
             write(unit=2502, iostat=io1, fmt=*)&
                  ptsU(1,node_list(size(h_ele)))
             write(unit=2502, iostat=io1, fmt=*)&
                  ptsU(2,node_list(size(h_ele)))
             write(unit=2502, iostat=io1, fmt=*)&
                  val(node_list(size(h_ele)))

          if(io1 /= 0) then
             write(0,*) 'could not write to file ' &
                  //trim(trim(filename)//'_meshpts.dat')
             stop
          end if
       end do
       end do
       end do
       close(unit=2502)
    end if      
        
     

open(unit=2502, file="./"//trim(trim(filename)//tc//'_D.dat'), &
         iostat=io1)

    ewrite(0,*) 'writing data'


    if(io1.ne.0) then
       write (0,*) 'Could not open file '
       stop
    else       

       write(unit=2502, iostat=io1, fmt=*) h(1)%mesh%shape%loc,n_layers+2,&
            ele_count(h(1)), filename, " Thickness ", alpha
       
       write(unit=2502, iostat=io1, fmt=*) trim(file_title)//" Thickness, time="&
            //trim(model_time) 

       do j=1,n_layers+2
       do i=1,1
       do ele = 1, ele_count(u(1))
          if (size(ele_val(h(1),ele))==3) then
             h_ele=>ele_nodes(positions,ele)
             ptsU=>ptsD1(:,1+3*(ele-1):3*ele)
          else
             h_ele=>ele_nodes(h(1),ele)
             ptsU=>ptsD2(:,1+6*(ele-1):6*ele)
          end if
          if(j<n_layers+1) then
             val(1:size(ele_val(h(1),ele)))=ele_val(h(j),ele)
          else if (j==n_layers+1) then
             val(1:size(ele_val(h(1),ele)))=-ele_val(bottom,ele)
          else 
             val(1:size(ele_val(h(1),ele)))=ele_val(h_0,ele)
          end if
             write(unit=2502, iostat=io1, fmt=*)&
                  ptsU(1,node_list(size(h_ele)))
             write(unit=2502, iostat=io1, fmt=*)&
                  ptsU(2,node_list(size(h_ele)))
             write(unit=2502, iostat=io1, fmt=*)&
                  val(node_list(size(h_ele)))

          if(io1 /= 0) then
             write(0,*) 'could not write to file ' &
                  //trim(trim(filename)//'_D.dat')
             stop
          end if
       end do
       end do
       end do
       close(unit=2502)
    end if

    deallocate(ptsD1,ptsD2)

  end subroutine py_plot_write_state

 function node_list(n_nodes)
    integer :: n_nodes
    integer :: node_list(n_nodes+1)

    select case(n_nodes)
    case(3)
       node_list=(/1,2,3,1/)
    case(6)
       node_list=(/1,2,3,5,6,4,1/)
    end select

  end function node_list

!!$ function shape_vector_dot_dshape(shape, vector, dshape, detwei) result (R)
!!$    !!< 
!!$    !!< Evaluate (Grad N1 dot vector) (N2)
!!$    !!<
!!$    real, dimension(:,:,:), intent(in) :: dshape
!!$    real, dimension(size(dshape,3),size(dshape,2)), intent(in) :: vector
!!$    type(element_type), intent(in) :: shape
!!$    real, dimension(size(dshape,2)) :: detwei
!!$
!!$    real, dimension(shape%loc,size(dshape,1)) :: R
!!$
!!$    integer :: iloc,jloc
!!$    integer :: dshape_loc, dim
!!$
!!$
!!$    dshape_loc=size(dshape,1)
!!$    dim=size(dshape,3)
!!$
!!$    forall(jloc=1:dshape_loc,iloc=1:shape%loc)
!!$       R(iloc,jloc)= dot_product(sum(dshape(jloc,:,:)*transpose(vector),2) &
!!$            *shape%n(iloc,:), detwei)
!!$    end forall
!!$
!!$  end function shape_vector_dot_dshape

   subroutine check_symmetry(helm,v_field)
  type(block_csr_matrix) :: helm
  type(vector_field) :: v_field
  integer, dimension(:), pointer :: ptr
  integer :: i,j,k,l,ierr
  PetscTruth :: symmetry_flag
  PetscReal :: tol


    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-check_symmetry",&
         symmetry_flag, ierr)
    
    if (symmetry_flag==PETSC_TRUE) then
       
       call PetscOptionsGetReal(PETSC_NULL_CHARACTER,"-check_symmetry",tol,&
            symmetry_flag,ierr)
       if (symmetry_flag==PETSC_FALSE) tol=1.0e-15
       do i =1,node_count(v_field)
          ptr=>row_m_ptr(helm,i)
          do j=1,size(ptr)
             do k=1,2*n_layers
                do l=1,2*n_Layers
                   if (val(helm,k,l,i,ptr(j))&
                        -val(helm,l,k,ptr(j),i)>tol) then
                      ewrite(-1,*) 'Error at',k,l, i, ptr(j)
                      ewrite(-1,*) val(helm,k,l,i,ptr(j)),&
                           val(helm,l,k,ptr(j),i)
                      FLAbort("Matrix is not symmetric!")
                   end if
                end do
             end do
          end do
       end do
       write(0,*) 'Matrix symmetric to order', tol
    end if

  end subroutine check_symmetry

subroutine get_model_u(mesh,u,D,d_perim,c,rho)
  type(dg_mesh) :: mesh
  type(vector_field) :: u(:)
  type(scalar_field) :: D(:)
  real :: d_perim(:),c,rho(:)

  !locals

  type(vector_field), pointer :: u_n(:)
  integer :: i,j,k,ele
  real, dimension(2,ele_loc(mesh%positions,1)) :: X_ele
  real, dimension(ele_ngi(D(1),1)) :: dlocgi,detwei
  real, dimension(2,ele_ngi(u(1),1)) :: ulocgi
  real, dimension(node_count(u(1))) :: u_rhs
  type(element_type), pointer :: shape_u,shape_h,shape_X

  allocate(u_n(n_layers))
  do k=1,n_layers
     call allocate(u_n(k),2,u(1)%mesh,'u')
  end do


  do k=1,n_layers

     u_rhs=0.0

     do ele=1,ele_count(U(1))

        shape_u=>ele_shape(u(1),ele)
        shape_X=>ele_shape(mesh%positions, ele)
        X_ele=ele_val(mesh%positions, ele)        
        call transform_to_physical(X_ele, shape_X,m=shape_u,detwei=detwei)

        dlocgi=ele_val_at_quad(D(k),ele)

        u_rhs(ele_nodes(u(1),ele))=shape_rhs(shape_u,rho(k)*c*detwei*(dlocgi-d_perim(k)))

     end do
     do j=1,node_count(u(1))
        u_n(k)%val(1)%ptr(j)=sum(&
             u_rhs(row_m(mesh%M_inv%blocks(1,1+2*(k-1)),j))*&
               row_val_ptr(mesh%M_inv%blocks(1,1+2*(k-1)),j))
     end do    

     u(k)%val(1)%ptr=u_n(k)%val(1)%ptr

  end do

  if (RIGID_LID) then
     
     D(1)%val=sum(d_perim)
     do k=2,n_layers
        D(1)%val=D(1)%val-D(k)%val
     end do

     call get_weighted_dg_inverse_mass_matrix(&
          mesh%m_inv%blocks(1,1),&
          u(1)%mesh,mesh%positions,D(1),rho(1))

     u_rhs=0.0

     do ele=1,ele_count(U(1))

        shape_u=>ele_shape(u(1),ele)
        shape_X=>ele_shape(mesh%positions, ele)
        X_ele=ele_val(mesh%positions, ele)        
        call transform_to_physical(X_ele, shape_X,m=shape_u,detwei=detwei)
        
        do k=2,n_layers

        dlocgi=ele_val_at_quad(D(k),ele)
        ulocgi=ele_val_at_quad(u(k),ele)


           u_rhs(ele_nodes(u(1),ele))=shape_rhs(shape_u,dlocgi*ulocgi(1,:))

        end do

     end do

     do j=1,node_count(u(1))
        u(1)%val(1)%ptr(j)=sum(&
             u_rhs(row_m(mesh%M_inv%blocks(1,1),j))*&
               row_val_ptr(mesh%M_inv%blocks(1,1),j))
     end do    

  end if

 call  py_plot_write("gal_u.dat","u",mesh%positions,v_f=u_n)

 do k=1,n_layers
    call deallocate(u_n(k))
 end do
 deallocate(u_n)

end subroutine get_model_u
        

  subroutine get_udash_udash(mesh,D,D_old,bottom,m,u,u_bar,rho)
    type(dg_mesh) :: mesh
    type(vector_field) :: m(:),u(:),u_bar
    type(scalar_field) :: D(:),bottom
    real :: D_old(:,:)
    real :: rho(:)

    ! locals

    real, dimension(2,node_count(u_bar)) :: u_rhs
    integer :: i,j,k,ele
    real, dimension(node_count(u_bar)):: u1,u2
    real, dimension(n_layers,node_count(m(1))) :: p_out
    type(vector_field),pointer,dimension(:) :: u_n
    type(scalar_field),pointer,dimension(:) :: p_n


    write(0,*) "outputting sum Du'u'"

    u_rhs=0.0
    allocate(u_n(1))
    allocate(p_n(n_layers))
    call allocate(u_n(1),2,u_bar%mesh,'u_dash')
    do k=1,n_layers
       call allocate(p_n(k),m(1)%mesh,'u.u')
    end do

    do k=1,n_layers

       u1=u(k)%val(1)%ptr-u_bar%val(1)%ptr
       u2=u(k)%val(2)%ptr-u_bar%val(2)%ptr
       
       do ele=1,ele_count(u_bar)
          call get_uu_element_contribution(mesh,u_bar,D(k),D_old(:,k),&
               rho(k),u_rhs,ele,flux="HLLC",u1=u1,u2=u2)
       end do

       call get_u_dot_u(mesh,m(1),u(1),u1,u2,rho(k),p_out(k,:))

       p_n(k)%val=p_out(k,:)
    end do

    do i=1,2
       do j=1,node_count(u_bar)
          u_n(1)%val(i)%ptr(j)=sum(&
               u_rhs(i,row_m(mesh%M_inv%blocks(1,i+2*n_layers),j))*&
               row_val_ptr(mesh%M_inv%blocks(1,i+2*n_layers),j))
       end do

    end do

    print*, 'Max u_n=', max(maxval(u_rhs(1,:)),maxval(u_rhs(2,:)))

    call  py_plot_write("gal_u_dash.dat","u'u'",mesh%positions,v_f=u_n)

    u_rhs=0.0

    do ele=1,ele_count(u_bar)
       call get_barotropic_mu_element_contribution(mesh,m,u_bar,D,D_old,&
             bottom,0*p_out,u_rhs,ele)
    end do
       

     do i=1,2
       do j=1,node_count(u_bar)
          u_n(1)%val(i)%ptr(j)=sum(&
               u_rhs(i,row_m(mesh%M_inv%blocks(1,i+2*n_layers),j))*&
               row_val_ptr(mesh%M_inv%blocks(1,i+2*n_layers),j))
       end do
    end do

    print*, 'Max u_n=', max(maxval(u_rhs(1,:)),maxval(u_rhs(2,:)))

    call  py_plot_write("gal_mu.dat","mu",mesh%positions,v_f=u_n)
    call  py_plot_write("gal_p.dat","p",mesh%positions,s_f=p_n)
    

    call deallocate(u_n(1))
     do k=1,n_layers
       call deallocate(p_n(k))
    end do
    deallocate(u_n,p_n)

  end subroutine get_udash_udash

subroutine get_u_dot_u(mesh,m,u,u_dash1,u_dash2,rho,p_out)
  type(dg_mesh) :: mesh
  type(vector_field):: m,u
  real :: rho,u_dash1(:),u_dash2(:),p_out(:)


  ! locals

  integer :: i,ele
  real :: uu_rhs(node_count(m)),p(node_count(m))
  type(element_type),pointer :: shape_m,shape_X,shape_u
  real, dimension(2,mesh%positions%mesh%shape%loc) :: X_ele
  real, dimension(u%mesh%shape%ngi) :: detwei
   real, dimension(2,u%mesh%shape%ngi) :: ulocgi,grad_uulocgi
   real, dimension(ele_loc(u,1),ele_ngi(u,1),2) :: du_t
   real, dimension(ele_loc(m,1),ele_ngi(m,1),2) :: dm_t

   type(csr_matrix) :: p_mat


   call allocate(p_mat,mesh%mass_h%sparsity)

   call zero(p_mat)
     uu_rhs=0.0
     
     do ele=1,ele_count(u)


        X_ele=ele_val(mesh%positions,ele)
        shape_m=> ele_shape(m,ele)
        shape_u=> ele_shape(u,ele)
        shape_X=> ele_shape(mesh%positions,ele)
        ulocgi(1,:)=matmul(u_dash1(ele_nodes(u,ele)),shape_u%n)
        ulocgi(2,:)=matmul(u_dash2(ele_nodes(u,ele)),shape_u%n)

        call transform_to_physical(X_ele, shape_X,m=shape_m,&
             dn_t=du_t,dm_t=dm_t,detwei=detwei)

        do i=1,2
           grad_uulocgi(i,:)=matmul(u_dash1(ele_nodes(u,ele))**2&
                +u_dash2(ele_nodes(u,ele))**2,&
                du_t(:,:,i))
        end do

        uu_rhs(ele_nodes(m,ele))=uu_rhs(ele_nodes(m,ele))&
             +shape_rhs(shape_m,rho*detwei*sum(ulocgi**2,1))

        call addto(p_mat,ele_nodes(m,ele),ele_nodes(m,ele),&
             shape_shape(shape_m,shape_m,detwei))

!        call addto(p_mat,ele_nodes(m,ele),ele_nodes(m,ele),&
!             dshape_dot_dshape(dm_t,dm_t,1.0e2*detwei))


     end do

     call petsc_solve(p,&
         p_mat,uu_rhs,&
         startfromzero=.true.)

     p_out=p_out+p

     call deallocate(p_mat)

end subroutine get_u_dot_u

subroutine py_plot_write(filename,type,positions,time,mod_time,v_f,s_F) 
    character(len=*), intent(in) :: filename
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in), optional  :: s_f(:)
    type(vector_field), intent(in), optional  :: v_f(:)
    integer, intent(in), optional :: time
    real, intent(in), optional :: mod_time
    character(len=*) :: type
   
    !locals
    integer :: io1, ele, nod, i,j,n_e,n_d,n_d2,n_loc
    integer, dimension(:), pointer :: u_ele
    integer, dimension(:), pointer :: h_ele
    real, dimension(:,:), pointer :: ptsD1, ptsD2
    real, dimension(:,:), pointer :: ptsU, ptsV
    integer, dimension(3) :: vertices
    real, dimension(6) :: val
    character(4) :: tc
    character(10) :: model_time
    character(128) :: file_title

    vertices(1) = 1
    vertices(2) = 3
    vertices(3) = 6
    
    if (present(time)) then
       write(tc,'(i4.4)') time
    else
       tc=""
    end if
    if (present(mod_time))   write(model_time,'(f10.3)') mod_time

    file_title=""
    if (.not. SHALLOW) then
       file_title= "Green Nagdhi Equation"
    else
       file_title= "Shallow Water Equation"
    end if

    allocate(ptsD1(2,3*ele_count(positions)),ptsD2(2,6*ele_count(positions)))



    do ele = 1, element_count(positions)
       ptsD1(:,1+3*(ele-1):3*ele) = ele_val(positions,ele)       
       ptsD2(:,(/1,3,6/)+6*(ele-1)) = ele_val(positions,ele)       
       ptsD2(:,(/2,4,5/)+6*(ele-1)) = &
            0.5*(ptsD2(:,(/1,6,6/)+6*(ele-1)) +&
            ptsD2(:,(/3,1,3/)+6*(ele-1)))
    end do

 open(unit=2502, file="./"//trim(trim(filename)), &
         iostat=io1)

    ewrite(0,*) 'writing data'


    if(io1.ne.0) then
       write (0,*) 'Could not open file '
       stop
    else       

       if (present(s_f)) then
          n_loc=s_f(1)%mesh%shape%loc
          n_d=size(s_f)
          n_d2=size(s_f)
          n_e=ele_count(s_f(1))
       else
          n_loc=v_f(1)%mesh%shape%loc
          n_d=2*size(v_f)
          n_d2=size(v_f)
          n_e=ele_count(v_f(1))
       end if


       write(unit=2502, iostat=io1, fmt=*) n_loc, n_d,&
            n_e, filename//" ",type, alpha
       
       write(unit=2502, iostat=io1, fmt=*) trim(file_title)//type//", time="&
            //trim(model_time) 



       do j=1,n_d2
       do i=1,n_d/n_d2
       do ele = 1,n_e
          if (present(s_f)) then
             h_ele=>ele_nodes(s_f(1),ele)
          else
             h_ele=>ele_nodes(v_f(1),ele)
          end if
          if (n_loc==3) then
             ptsU=>ptsD1(:,1+3*(ele-1):3*ele)
          else
             ptsU=>ptsD2(:,1+6*(ele-1):6*ele)
          end if
          
          if (present(s_f)) then
             val(1:n_loc)=ele_val(s_f(j),ele)
          else
             val(1:n_loc)=ele_val(v_f(j),ele,i)
          end if
             write(unit=2502, iostat=io1, fmt=*)&
                  ptsU(1,node_list(n_loc))
             write(unit=2502, iostat=io1, fmt=*)&
                  ptsU(2,node_list(n_loc))
             write(unit=2502, iostat=io1, fmt=*)&
                  val(node_list(n_loc))

          if(io1 /= 0) then
             write(0,*) 'could not write to file ' &
                  //trim(trim(filename))
             stop
          end if
       end do
       end do
       end do
       close(unit=2502)
    end if      

  end subroutine py_plot_write

function get_centroid(points) result(centre)
  real :: points(2,3)
  real :: centre(2)

  centre=sum(points,2)/3.0

end function get_centroid

function get_CFL(mesh,c,dt) result(CFL)
  type(dg_mesh) :: mesh
  real :: c,dt,CFL(2)

  ! locals

  integer :: ele, ele_2
  integer, dimension(:), pointer :: neigh_con, neigh_pos
  real, dimension(2,3) :: X_ele, X_2_ele
  real :: min_dist, dist

  min_dist=INFINITY

  do ele=1,ele_count(mesh%positions)

     neigh_pos=>ele_neigh(mesh%positions,ele)
     X_ele=ele_val(mesh%positions,ele)

     do ele_2=1,size(neigh_pos)

        if (neigh_pos(ele_2)>0) then

           X_2_ele=ele_val(mesh%positions,neigh_pos(ele_2))

           dist=sqrt(sum((get_centroid(X_ele)-get_centroid(X_2_ele))**2))
           min_dist=min(min_dist,dist)

        end if

     end do

  end do

  CFL(1)=2.0*dt/(min_dist)*c
  CFL(2)=min_dist

end function get_CFL

subroutine make_D_matrix(mesh,D,d0)

type(dg_mesh) :: mesh
type(scalar_field)  :: D
real :: d0

real, dimension(mesh%positions%dim,ele_loc(mesh%positions,1)) :: X_ele
type(element_type), pointer :: shape_X,shape_h
real, dimension(ele_ngi(mesh%positions,1)) :: detwei

integer :: ele, i,j
integer, dimension(:), pointer :: ele_h

call zero (mesh%D_mat)


do i = 1,size(mesh%CMC%colm)
   do j = 1,size(mesh%CMC%colm(i)%ptr)
      call addto(mesh%D_mat,i,mesh%CMC%colm(i)%ptr(j),mesh%CMC%val(i)%ptr(j))
   end do
end do
    

mesh%D_mat%val=mesh%D_mat%val*g0*d0*(dt**2)

do ele=1,ele_count(D)
   ele_h=>ele_nodes(D,ele)
   shape_X=>ele_shape(mesh%positions, ele)
   shape_h=>ele_shape(D, ele)

   X_ele=ele_val(mesh%positions, ele)
   call transform_to_physical(X_ele, shape_X,detwei=detwei)

   call addto(mesh%D_mat,ele_h,ele_h,shape_shape(shape_h,shape_h,detwei))
end do


end subroutine make_D_matrix

function rotate_bdry(y,bcs) result(y_rot)
  real :: y(:,:)

  real, allocatable :: y_rot(:,:),y_tmp(:,:)
  type(bc_info) :: bcs

  ! locals 

  integer :: i
  integer, pointer :: node
  real, pointer :: tangent(:)

  allocate(y_rot(size(y,1),size(y,2)),y_tmp(size(y,1),size(y,2)))
           

  y_rot=y
  y_tmp=y

  do i=1,bcs%N_tangents
     node=>bcs%tangent_list(i)
     tangent=>bcs%tangents(:,I)

     y_tmp(1,node)=y(1,node)*tangent(1)+y(2,node)*tangent(2)
     y_tmp(2,node)=0.0!y(1,node)*tangent(2,node)-y(2,node)*tangent(1,node)
  end do
  do i=1,bcs%N_tangents
     node=>bcs%tangent_list(i)
     tangent=>bcs%tangents(:,i)

     y_rot(1,node)=y_tmp(1,node)*tangent(1)!+y_rot(2,node)*tangent(2)
     y_rot(2,node)=y_tmp(1,node)*tangent(2)!-y_rot(2,node)*tangent(1)
  end do

end function rotate_bdry

end module MLCM_scalar_formulation
