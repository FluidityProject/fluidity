#include "fdebug.h"

module Mesh_tools
  use sparse_tools
  use global_parameters_gallopede
  use fldebug
  use elements
  use FETools
  use data_structures
!  use gallopede_solvers
  use Fields
  use Fields_data_types

  implicit none
  
  public::  adapt_timestep, set_lump_mass, get_tangents
  private
  
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"

  contains

!!$    subroutine project_to(D,Din,shape,shapein,shapeX, &
!!$         EVList,EVListin,EVList_X,X,Y,Mass,allocate_mass,form_mass, &
!!$         multiply_mass)
!!$      implicit none
!!$      type(element_type), intent(in) :: shape, shapein,shapeX
!!$      type(csr_matrix), intent(inout) :: Mass
!!$      integer, dimension(N_elements*shapein%loc),target,intent(in) :: EVListin
!!$      integer, dimension(N_elements*shape%loc),target,intent(in) :: EVList
!!$      integer, dimension(N_elements*shapeX%loc),target,intent(in) :: EVList_X
!!$      logical, intent(in) :: allocate_mass, form_mass, multiply_mass
!!$      real, dimension(N_verts), intent(in) :: X,Y
!!$      real, dimension(:), intent(in), target  :: Din
!!$      real, dimension(:), intent(out), target :: D
!!$
!!$      !locals
!!$      integer :: ele
!!$      integer, dimension(:), pointer :: D_ele, Din_ele,X_ele
!!$      real, dimension(shape%loc,shape%loc) :: Q      
!!$      real, dimension(shape%loc,shapein%loc) :: Qm
!!$      real, dimension(shapein%loc,shapein%loc) :: Qin
!!$      real, dimension(shape%ngi) :: detwei
!!$      real, dimension(2,3) :: ele_X !coordinates for triangle vertices
!!$      real, dimension(:), pointer :: RHS
!!$      KSPType :: ksp_type
!!$      PCType  :: pc_type
!!$
!!$      pc_type = PCSOR  ! Preconditioner context
!!$      ksp_type = KSPCG ! Krylov subspace context
!!$
!!$      ewrite(1,*) 'subroutine project_to'
!!$
!!$      ewrite(2,*) 'maxval(Din)', maxval(abs(Din))
!!$
!!$      if(allocate_mass) then
!!$         call POSINM(Mass,N_elements,size(D), &
!!$              shape%loc,Evlist,&
!!$              size(D), shape%loc, EVList)
!!$      end if
!!$      if(form_mass) then
!!$         call zero(Mass)
!!$         ele_loop_m: do ele = 1, N_elements
!!$            X_ele=>EVList_X((ELE-1)*shapeX%loc+1:ELE*shapeX%loc)
!!$            D_ele=>EVList((ELE-1)*shape%LOC+1:ELE*shape%LOC)
!!$            ele_X(1,:)=X(X_ele)
!!$            ele_X(2,:)=Y(X_ele)
!!$            call transform_to_physical(ele_X,shapeX,detwei = detwei)
!!$            Q = shape_shape(shape,shape,detwei)
!!$            call addto(Mass,D_ele,D_ele,Q)
!!$         end do ele_loop_m
!!$      end if
!!$
!!$      
!!$
!!$      if(multiply_mass) then
!!$         allocate( RHS( size(D) ) )
!!$         RHS = 0.
!!$
!!$         ele_loop: do ele = 1, N_elements
!!$            call transform_to_physical(ele_X, shapein, detwei = detwei)
!!$            D_ele=>EVList((ELE-1)*shape%LOC+1:ELE*shape%LOC)
!!$            Din_ele=>EVListin((ELE-1)*shapein%loc+1:ELE*shapein%loc)
!!$            X_ele=>EVList_X((ELE-1)*shapeX%loc+1:ELE*shapeX%loc)
!!$            Qm= shape_shape(shape,shapein,detwei)
!!$            RHS(D_ele) = RHS(D_ele) + matmul(Qm,Din(Din_ele))
!!$         end do ele_loop
!!$      else
!!$         RHS => Din
!!$      end if
!!$
!!$      call gallopede_solve(D,Mass,RHS,ksp_type,pc_type,1e-8,200)
!!$
!!$      if(multiply_mass) then
!!$         deallocate( RHS )
!!$      end if
!!$
!!$      ewrite(1,*) 'END subroutine project_to'
!!$
!!$    end subroutine project_to

!!$    subroutine get_quad_mesh(EVList_h,EVList_X,EVList_Field,EEList, &
!!$         X, Y,n_quads,bcs,quadbcs) 
!!$      !subroutine to construct a quadratic mesh from a linear mesh
!!$      implicit none
!!$      type(csr_matrix), intent(in) :: EEList
!!$      integer, dimension(N_elements*3), target, intent(in) :: EVList_X
!!$      integer, dimension(N_elements*3), target, intent(in) :: EVList_Field
!!$      integer, dimension(N_elements*6), target, intent(out) :: EVList_h
!!$      type(bc_info), intent(in) , optional :: bcs
!!$      type(bc_info), intent(out) , optional :: quadbcs
!!$      real, dimension(:), intent(in) :: X, Y
!!$      integer, intent(out) :: n_quads
!!$
!!$      !locals
!!$      type(csr_matrix) :: edgelist
!!$      ! List of neighbours of current element.
!!$      integer, dimension(:), allocatable :: neigh
!!$      integer, dimension(:), allocatable :: vertlist
!!$      integer :: ele, nod, ni, count, iloc, ele_2
!!$      integer, dimension(:), pointer :: X_ele, h_ele
!!$      integer, dimension(3) :: lverts, ledges
!!$      integer :: interior_j, tangent_j
!!$      real :: norm
!!$      integer :: i,j,k
!!$
!!$      ewrite(1,*)("subroutine get_quad_mesh")
!!$
!!$      EVList_h = 0
!!$
!!$      lverts = (/ 1, 3, 6 /)
!!$      ledges = (/ 5, 4, 2 /)
!!$
!!$      ewrite(2,*)("Allocating memory for neigh");
!!$      allocate(neigh(row_length(EEList,1)))
!!$      allocate( vertlist(N_free) )
!!$
!!$      edgelist = clone(EEList, type=CSR_INTEGER)
!!$      call zero( edgelist )
!!$
!!$      vertlist = 0
!!$      count = 0
!!$      do ele = 1, N_elements
!!$         !vertex entries
!!$
!!$         X_ele=>EVList_Field((ELE-1)*3+1:ELE*3)
!!$         do iloc = 1,3
!!$            if(vertlist(X_ele(iloc))==0) then
!!$               count = count + 1
!!$               nod = count
!!$               vertlist(X_ele(iloc))=count
!!$            else
!!$               nod = vertlist(X_ele(iloc))
!!$            end if
!!$            EVList_h((ele-1)*6+lverts(iloc)) = nod
!!$         end do
!!$
!!$         !edge entries
!!$
!!$         neigh=row_m(EEList,ele)
!!$         ASSERT(size(neigh)==3)
!!$
!!$         !loop over neighbours of ele         
!!$         neighbourloop: do ni=1,size(neigh)
!!$            ele_2=neigh(ni)
!!$
!!$            ASSERT(ele_2.ne.ele)
!!$
!!$            if (ele_2==0) then 
!!$               count = count + 1
!!$               nod = count
!!$            else
!!$               if(ival(edgelist,ele,ele_2)==0) then
!!$                  count = count + 1
!!$                  nod = count
!!$                  call set(edgelist,ele,ele_2,nod)
!!$                  call set(edgelist,ele_2,ele,nod)
!!$               else
!!$                  nod = ival(edgelist,ele,ele_2)
!!$               end if
!!$            end if
!!$            EVList_h((ele-1)*6 + ledges(ni)) = nod
!!$
!!$         end do neighbourloop
!!$
!!$      end do
!!$
!!$      n_quads = count
!!$
!!$      deallocate( neigh )
!!$      deallocate( vertlist )
!!$
!!$      if(present(bcs).and.present(quadbcs)) then
!!$
!!$         ewrite(2,*) size(bcs%bc_marker)
!!$
!!$         !getting bc info
!!$         if(associated(quadbcs%bc_marker)) then
!!$            deallocate(quadbcs%bc_marker)
!!$            quadbcs%bc_marker => null()
!!$         end if
!!$         allocate(quadbcs%bc_marker(n_quads))
!!$         quadbcs%bc_marker = 0
!!$
!!$         do ele = 1, N_elements
!!$            X_ele => EVList_Field((ele-1)*3+1:ele*3)
!!$            h_ele => EVList_h((ele-1)*6+1:ele*6)
!!$
!!$            if(bcs%bc_marker(X_ele(1))==2) then
!!$               quadbcs%bc_marker(h_ele(1))=2
!!$               quadbcs%bc_marker(h_ele(2))=2
!!$               quadbcs%bc_marker(h_ele(4))=2
!!$            end if
!!$            if(bcs%bc_marker(X_ele(2))==2) then
!!$               quadbcs%bc_marker(h_ele(2))=2
!!$               quadbcs%bc_marker(h_ele(3))=2
!!$               quadbcs%bc_marker(h_ele(5))=2
!!$            end if
!!$            if(bcs%bc_marker(X_ele(3))==2) then
!!$               quadbcs%bc_marker(h_ele(4))=2
!!$               quadbcs%bc_marker(h_ele(5))=2
!!$               quadbcs%bc_marker(h_ele(6))=2
!!$            end if
!!$
!!$            if(bcs%bc_marker(X_ele(1))==1) then
!!$               quadbcs%bc_marker(h_ele(1))=1
!!$               quadbcs%bc_marker(h_ele(2))=1
!!$               quadbcs%bc_marker(h_ele(4))=1
!!$            end if
!!$            if(bcs%bc_marker(X_ele(2))==1) then
!!$               quadbcs%bc_marker(h_ele(2))=1
!!$               quadbcs%bc_marker(h_ele(3))=1
!!$               quadbcs%bc_marker(h_ele(5))=1
!!$            end if
!!$            if(bcs%bc_marker(X_ele(3))==1) then
!!$               quadbcs%bc_marker(h_ele(4))=1
!!$               quadbcs%bc_marker(h_ele(5))=1
!!$               quadbcs%bc_marker(h_ele(6))=1
!!$            end if
!!$
!!$         end do
!!$
!!$         ewrite(2,*) 'Construct lifted ordering'
!!$         if(associated(quadbcs%lifted_ordering)) then
!!$            deallocate(quadbcs%lifted_ordering)
!!$            quadbcs%lifted_ordering => null()
!!$         end if
!!$         allocate(quadbcs%lifted_ordering( n_quads ))
!!$
!!$         tangent_j = 0
!!$         interior_j = 0
!!$         do i = 1, n_quads
!!$            select case(quadbcs%bc_marker(i))
!!$            case (1)
!!$               interior_j = interior_j + 1
!!$               quadbcs%lifted_ordering(i) = interior_j
!!$            case (2)
!!$               tangent_j = tangent_j + 1
!!$               quadbcs%lifted_ordering(i) = tangent_j
!!$            end select
!!$         end do
!!$
!!$         quadbcs%N_interior = interior_j
!!$         quadbcs%N_tangents = tangent_j
!!$
!!$         ewrite(2,*) 'Construct interior list and tangent list'
!!$
!!$         if(associated(quadbcs%interior_list)) then
!!$            deallocate(quadbcs%interior_list)
!!$            quadbcs%interior_list => null()
!!$         end if
!!$         allocate( quadbcs%interior_list( quadbcs%N_interior ) )
!!$         if(associated(quadbcs%tangent_list)) then
!!$            deallocate(quadbcs%tangent_list)
!!$            quadbcs%tangent_list => null()
!!$         end if
!!$         allocate( quadbcs%tangent_list( quadbcs%N_tangents ) )
!!$         do i = 1, n_quads
!!$            select case(quadbcs%bc_marker(i))
!!$            case (1)
!!$               quadbcs%interior_list(quadbcs%lifted_ordering(i)) = i
!!$            case (2)
!!$               quadbcs%tangent_list(quadbcs%lifted_ordering(i)) = i
!!$            end select
!!$         end do
!!$
!!$         ewrite(2,*) 'Construct tangents on midpoints'
!!$
!!$         if(associated(quadbcs%tangents)) then
!!$            deallocate(quadbcs%tangents)
!!$            quadbcs%tangents => null()
!!$         end if
!!$         allocate( quadbcs%tangents( 2,quadbcs%N_interior ) )      
!!$
!!$         do ele = 1, N_elements
!!$
!!$            X_ele => EVList_X((ele-1)*3+1:ele*3)
!!$            h_ele => EVList_h((ele-1)*6+1:ele*6)
!!$
!!$            if(quadbcs%bc_marker(h_ele(2)) == 2) then
!!$               i = quadbcs%lifted_ordering(h_ele(2))
!!$               quadbcs%tangents(1,i) = Y(X_ele(1)) - Y(X_ele(2))
!!$               quadbcs%tangents(2,i) = - X(X_ele(1)) + X(X_ele(2))
!!$               norm = sqrt(quadbcs%tangents(1,i)**2 + &
!!$                    quadbcs%tangents(2,i)**2)
!!$               quadbcs%tangents(1,i) = quadbcs%tangents(1,i)/norm            
!!$               quadbcs%tangents(2,i) = quadbcs%tangents(2,i)/norm
!!$            end if
!!$
!!$            if(quadbcs%bc_marker(h_ele(4)) == 2) then
!!$               i = quadbcs%lifted_ordering(h_ele(4))
!!$               quadbcs%tangents(1,i) = Y(X_ele(1)) - Y(X_ele(3))
!!$               quadbcs%tangents(2,i) = - X(X_ele(1)) + X(X_ele(3))
!!$               norm = sqrt(quadbcs%tangents(1,i)**2 + &
!!$                    quadbcs%tangents(2,i)**2)
!!$               quadbcs%tangents(1,i) = quadbcs%tangents(1,i)/norm            
!!$               quadbcs%tangents(2,i) = quadbcs%tangents(2,i)/norm
!!$            end if
!!$
!!$            if(quadbcs%bc_marker(h_ele(5)) == 2) then
!!$               i = quadbcs%lifted_ordering(h_ele(5))
!!$               quadbcs%tangents(1,i) = Y(X_ele(3)) - Y(X_ele(2))
!!$               quadbcs%tangents(2,i) = - X(X_ele(3)) + X(X_ele(2))
!!$               norm = sqrt(quadbcs%tangents(1,i)**2 + &
!!$                    quadbcs%tangents(2,i)**2)
!!$               quadbcs%tangents(1,i) = quadbcs%tangents(1,i)/norm            
!!$               quadbcs%tangents(2,i) = quadbcs%tangents(2,i)/norm
!!$            end if
!!$
!!$         end do
!!$
!!$      end if
!!$
!!$      ewrite(1,*)("END subroutine get_quad_mesh")
!!$
!!$    end subroutine get_quad_mesh

    subroutine adapt_timestep(Mesh, u, dt_out)
      type(dg_mesh), intent(in) :: mesh
      real, intent(in), dimension(N_vels*N_Layers*2) , target:: u
      real, intent(out) :: dt_out

      !local variables
      real, pointer, dimension(:) :: u1,u2      
      integer :: ele
      integer, dimension(:), pointer :: X_ele,u_ele
      real, dimension(2,mesh%nm%loc) :: ele_X
      real, dimension(mesh%nh%ngi) :: detwei,u1locgi,u2locgi
      real :: dx, volume, umean
      real :: max_dt
      real, parameter :: CFL_LIMIT = 0.3, growth_factor = 1.3
      integer :: layer_i
     
      ewrite(1,*) '     subroutine adapt_timestep(Mesh, u)'

      max_dt = dtmax

      do layer_i = 0, N_Layers-1
         u1 => u(1 + Layer_i*2*N_vels:N_vels + Layer_i*2*N_vels)
         u2 => u(N_vels + 1 + Layer_i*2*N_vels:2*N_vels + Layer_i*2*N_vels)

         do ele = 1, N_Elements
            
            u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
            X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
            ele_X(1,:)=mesh%X(X_ele)
            ele_X(2,:)=mesh%Y(X_ele)
            
            call transform_to_physical(ele_X, n=mesh%nm, detwei = detwei)
            volume = sum(detwei)
            dx = sqrt(volume)
            u1locgi = matmul(transpose(mesh%nu%n),u1(u_ele))
            u2locgi = matmul(transpose(mesh%nu%n),u2(u_ele))
            umean = sqrt(sum(detwei*(u1locgi*u1locgi + u2locgi*u2locgi)) &
                 / volume)

            if( umean * max_dt > CFL_LIMIT * dx) then
               max_dt = CFL_LIMIT * dx / umean
            end if

         end do
      end do

      dt_out=dt
      if(dt > max_dt) then
         dt_out = max_dt
      elseif(dt<0.5*max_dt) then
         dt_out = min(dtmax, dt * growth_factor)
         dt_out = min(dt_out, max_dt)
      end if
      
      ewrite(1,*) 'dt = ', dt_out, ' dt_CFL=', max_dt/CFL_LIMIT

      ewrite(1,*) 'END subroutine adapt_timestep(Mesh, u)'

      if (dt_out <1e-3*dtmax) then
         ewrite(1,*) 'Timestep gone unstable'
         stop
      end if

    end subroutine adapt_timestep

    subroutine set_lump_mass(mass,lump_mass,no_of_els)
      real, pointer, dimension(:) :: lump_mass
      type(csr_matrix) :: mass
      integer :: no_of_els
      !locals 
      integer :: i
      
      allocate(lump_mass(no_of_els))
      do i=1,no_of_els
         lump_mass(i)=sum(row_val_ptr(mass,i))
      end do
      
      print*, lump_mass
      
    end subroutine set_lump_mass

  subroutine get_tangents(bcs,mesh,bc_field,u,positions)
    type(bc_info), intent(inout) :: bcs
    type(dg_mesh), intent(in) :: mesh
    integer,dimension(:), intent(in) :: bc_field
    type(vector_field), intent(in) :: positions,u

    !locals
    integer :: i,ele
    real :: norm
    real, allocatable, dimension(:,:) :: normal_list

    allocate( normal_list(2,N_vels) )

    normal_list = 0.0

    do ele = 1, ele_count(u)
       
       call assemble_tangents_element_contribution(positions, &
            bc_field,u,normal_list,ele)       

    end do

    do i = 1, bcs%N_tangents
       norm = sqrt(normal_list(1,bcs%tangent_list(i))**2 + &
            normal_list(2,bcs%tangent_list(i))**2)

       assert(norm >0)
       bcs%tangents(1,i) = -normal_list(2,bcs%tangent_list(i))/ norm
       bcs%tangents(2,i) = normal_list(1,bcs%tangent_list(i))/ norm
    end do

  end subroutine get_tangents

  subroutine assemble_tangents_element_contribution(positions, &
       bc_field,u,normal_list,ele)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions,u
    integer, dimension(:), intent(in):: bc_field
    real, dimension(:,:), intent(inout) :: normal_list
    ! Shape functions.
    type(element_type), pointer :: shape_h, shape_X,shape_Xf
    ! Locations of nodes.
    real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
    real, dimension(positions%dim,face_loc(positions,1)) :: X_ele_f
    real, dimension(positions%dim,face_ngi(positions,1)) :: normal
    real, dimension(ele_loc(u,ele)) :: bc_ele
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei
    real, dimension(face_ngi(positions,1)) :: detwei_f
    ! Derivatives of shape function:
    real, dimension(ele_loc(u,ele), &
         ele_ngi(u,ele), positions%dim) :: dshape_h
    real, dimension(positions%dim,ele_loc(u,ele)) :: grad_phi

    ! Node numbers of field element.
    integer, dimension(:), pointer :: ele_h,ele_2,ele_mesh
    integer :: i,j,face_u,face_mesh
    integer :: ele_h_f(2)

    ele_h=>ele_nodes(u, ele)
    shape_h=>ele_shape(u, ele)
    shape_X=>ele_shape(positions, ele)

    ele_2=> ele_neigh(u,ele)
    ele_mesh=> ele_neigh(positions,ele)

    ! Locations of local vertices.
    X_ele=ele_val(positions, ele)
    bc_ele=bc_field(ele_h)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(X_ele, shape_X, m=shape_h, &
         dm_t=dshape_h, detwei=detwei)

    do j=1,size(ele_2)

       face_u=ele_face(u,ele,ele_2(j))
       face_mesh=ele_face(positions,ele,ele_mesh(j))
       ele_h_f=face_global_nodes(u, face_u)
       X_ele_f=face_val(positions,face_mesh)
       shape_Xf=>face_shape(positions, face_mesh)

       call transform_bdy_to_physical(X_ele,X_ele_f, shape_X, shape_Xf,&
           detwei_f=detwei_f,normal=normal)

!    grad_phi = dshape_rhs(dshape_h,detwei)

       if (ele_2(j)<0) then

          normal_list(1,ele_h_f) = normal_list(1,ele_h_f)+normal(1,1)
          normal_list(2,ele_h_f) = normal_list(2,ele_h_f)+normal(2,1)

       end if

    end do

  end subroutine assemble_tangents_element_contribution

end module Mesh_tools
