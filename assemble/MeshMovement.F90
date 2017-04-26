#include "fdebug.h"

module meshmovement

  use global_parameters
  use fldebug
  use vector_tools, only: solve, mat_diag_mat, eigendecomposition_symmetric
  use element_numbering
  use elements
  use shape_functions
  use spud
  use sparse_tools
  use fields_base
  use global_numbering
  use eventcounter
  use fetools
  use unittest_tools
  use fields
  use state_module
  use vtk_interfaces
  use sparse_matrices_fields
  use solvers
  use fefields
  use field_derivatives
  use sparsity_patterns
  use sparsity_patterns_meshes

  implicit none
  integer,save :: MeshCount=0

  interface

     subroutine set_debug_level(level)
       implicit none
       integer, intent(in) :: level
     end subroutine set_debug_level

     subroutine reset_debug_level
     end subroutine reset_debug_level
  end interface
  
  private
  
  public :: move_mesh_imposed_velocity, move_mesh_pseudo_lagrangian, movemeshy

contains
  subroutine movemeshy(state,move_option,TimeStep)
    implicit none

    type(state_type), intent(inout) :: state
    integer, intent(in),optional :: move_option
    real, intent(in),optional :: TimeStep
    type(mesh_type) :: mesh
    type(csr_sparsity), pointer :: U_sparsity
    type(scalar_field), pointer :: ifield
    type(vector_field) :: positions
    type(vector_field), pointer :: velocity,gridvelocity
    type(vector_field) :: gridvelocity_yves
    type(vector_field) :: relativevelocity, newgridvelocity
    real, dimension(:,:), allocatable::ugnew
    real, dimension(:,:), allocatable::ugold
    real, dimension(:,:), allocatable::b
    real, dimension(:,:), allocatable::bvec
    real, dimension(:,:), allocatable::bfunc
    real, dimension(:,:,:), allocatable::block5x5 !
    !should replace it with  type(block_csr_matrix)::block5x5
    real, dimension(:,:,:), allocatable::block5x5l !
    real, dimension(:,:), allocatable::bnew       !
    type(vector_field), target:: linear_gridvelocity,linear_velocity
    type(scalar_field), target:: linear_ifield
    integer unit,II
    real mini, l2norm2,Dij,minidij,A0,beta,D0,sp
    INTEGER i,j,k,node,ele,inod,it,iteration,d,node2,lstat,noniti,nonlineariteration
    REAL alpha,LAMXX,LAMYY,LAMZZ,EPSILON


    type(tensor_field) :: hessian
    type(block_csr_matrix)::MatrixA
    type(csr_matrix)::Ablock

    real,dimension(:),allocatable ::Unplus1

    type(vector_field) :: MatrixAlumped
    type(scalar_field) :: masslump,lumpedmass

    type(scalar_field) :: diagK
    type(scalar_field) :: Fe
    type(scalar_field) :: F1
    type(scalar_field) :: Fspring
    type(tensor_field) :: hessianFe

    real FeweightScal,springweight,expo_weight
    type(scalar_field) :: field1

    real, dimension(:,:),allocatable::evectors
    real, dimension(:),allocatable::evalues
    real,dimension(:,:),allocatable:: A
    real,dimension(:,:),allocatable::Bn

    character(70) :: filename
    character(len=FIELD_NAME_LEN)::fieldname
    character(50) :: dumpnum
    type(csr_sparsity), pointer ::Nnlist2
    logical spring,expo,iterating,scaling

    type(vector_field) :: gradfe
    !cdim is an integer defining the size of the system to solve
    !if cdim ==2 =>vertical movement only , system size=5X5
    !if cdim ==1 =>add horizontal movement,constraint in y or x , system size=4X4
    !if cdim ==0 =>no constraint, no lagrange multiplier, system size=3X3
    integer cdim,mdim,lagrangian_multi,option
    
    real maxi,crazynumber,tempr
    real, dimension(:),allocatable::F1weight,F2weight,Fspringweight,Feweight
    real, dimension(:),allocatable::rangeX,speedrange,DeltaX
    integer,dimension(:),allocatable::Nx
    real,dimension(:),allocatable::LAM,Xpos
    real, dimension(:),allocatable::speed,speed2
    real,dimension(:),allocatable:: X2,MaxX,MinX
    real,dimension(:),allocatable:: Pos2
    real,dimension(:,:),allocatable:: AHFe
    type(scalar_field), dimension(3) :: DFe   
    real rhorange
    character(len=OPTION_PATH_LEN) :: functional_path,tmp_path
    logical,dimension(3)::derivatives,blocking
    type(csr_matrix) :: P
    type(scalar_field) :: dg_scalar, cg_scalar
    real dt,pi

    call IncrementEventCounter(EVENT_MESH_MOVEMENT)

   
    mesh = extract_mesh(state, "CoordinateMesh")
    ewrite(3,*) "coordinate mesh extracted "
    mdim=mesh_dim(mesh)

    allocate(F1weight(mdim))
    allocate(F2weight(mdim))
    allocate(Fspringweight(mdim))
    allocate(Feweight(mdim))
    allocate(rangeX(mdim))
    allocate(speedrange(mdim))
    allocate(deltaX(mdim))
    allocate(Nx(mdim))
    allocate(LAM(mdim))
    allocate(Xpos(mdim))
    allocate(speed(mdim))
    allocate(speed2(mdim))
    allocate(X2(mdim))
    allocate(MaxX(mdim))
    allocate(MinX(mdim))
    allocate(Pos2(mdim))
    allocate(AHFe(mdim,mdim))

    dt=0.0
    lagrangian_multi=0
    if(present(move_option)) then 
       option=move_option
    else
       option=0
    endif
    if(present(TimeStep)) dt=TimeStep



    functional_path='/mesh_adaptivity/mesh_movement/vertical_ale/physical_functionals'
    tmp_path=trim('/mesh_adaptivity/mesh_movement/vertical_ale')

    spring=.false.
    expo=.false.
    iterating=.true.
    scaling=.false.

    blocking=.false.
    F1weight=1.0
    F2weight=1.0
    Fspringweight=1.0
    Feweight=1.0
    rangeX=1.0
    rhorange=0.0
    speedrange=0.02
    DeltaX=0.0



    positions= extract_vector_field(state, "IteratedCoordinate")
    velocity=>extract_vector_field(state,"Velocity")
    gridvelocity=>extract_vector_field(state,"IteratedGridVelocity")
    if (.not. gridvelocity%mesh==positions%mesh) then
       ewrite(3,*) "Why is gridvelocity not on the CoordinateMesh?"
       call allocate(linear_gridvelocity, gridvelocity%dim,mesh)
       call zero(linear_gridvelocity)
       call allocate(lumpedmass, mesh, "LumpedMass")
       call compute_lumped_mass(positions, lumpedmass)
       ! Invert lumped mass.
       lumpedmass%val=1./lumpedmass%val
       P=compute_projection_matrix(mesh, gridvelocity%mesh, positions)
  
       do i=1,linear_gridvelocity%dim
          cg_scalar=extract_scalar_field_from_vector_field(linear_gridvelocity, i)
          dg_scalar=extract_scalar_field_from_vector_field(gridvelocity, i)
          call mult(cg_scalar,P,dg_scalar)
       end do
       ! Apply inverted lumped mass to projected quantity.
       call scale(linear_gridvelocity, lumpedmass)
       
       gridvelocity=>linear_gridvelocity 
       call deallocate(lumpedmass) !!hmm
       gridvelocity => linear_gridvelocity
    end if

    if (.not. velocity%mesh==positions%mesh) then
       ewrite(3,*) "mapping velocity if DG"
       call allocate(linear_velocity,velocity%dim,mesh)
       call allocate(lumpedmass, mesh, "LumpedMass")
       call compute_lumped_mass(positions, lumpedmass)
       ! Invert lumped mass.
       lumpedmass%val=1./lumpedmass%val
       P=compute_projection_matrix(mesh, velocity%mesh, positions)
  
       do i=1,linear_velocity%dim
          cg_scalar=extract_scalar_field_from_vector_field(linear_velocity, i)
          dg_scalar=extract_scalar_field_from_vector_field(velocity, i)
          call mult(cg_scalar,P,dg_scalar)
       end do
       ! Apply inverted lumped mass to projected quantity.
       call scale(linear_velocity, lumpedmass)
       
       velocity=>linear_velocity 
       call deallocate(lumpedmass) !!hmm
    end if


    if(have_option('/simulation_name'))then
       call get_option('/timestepping/timestep',dt)
       if(have_option(trim(tmp_path))) then
          ewrite(3,*) "vertical_ale!!!!!!!!! "
          if(have_option(trim(functional_path)))then

             ewrite(3,*) "physical_functionals!!!!!!!!! "

             if(have_option(trim(functional_path)//'/minimise_relative_velocity_dot_grad_density')) then
                option=1
                call get_option(trim(functional_path)//'/minimise_relative_velocity_dot_grad_density/name',fieldname)
             elseif(have_option(trim(functional_path)//'/use_hessian_density')) then
                option=2
                call get_option(trim(functional_path)//'/use_hessian_density/name',fieldname)
             endif
          endif
          
          if(have_option(trim(tmp_path)//'/mesh_quality_terms')) then
             ewrite(3,*) "quality functionals!!!!!!!!! "

             if(have_option(trim(tmp_path)//'/mesh_quality_terms/spring_term')) then
                spring=.true.
                call get_option(trim(tmp_path)//'/mesh_quality_terms/spring_term',springweight)
                Fspringweight=springweight
             endif
             if(have_option(trim(tmp_path)//'/mesh_quality_terms/exponential_term')) then
                expo=.true.
                call get_option(trim(tmp_path)//'/mesh_quality_terms/exponential_term',expo_weight)
             endif
          endif

          if(have_option(trim(tmp_path)//'/block_nodes_in_x')) then
             lagrangian_multi=lagrangian_multi+1
             blocking(1)=.true.
             ewrite(3,*) "No grid movement in X "
          endif
          if(have_option(trim(tmp_path)//'/block_nodes_in_y')) then
             lagrangian_multi=lagrangian_multi+1
             blocking(2)=.true.
             ewrite(3,*) "No grid movement in Y "
          endif
          if(have_option(trim(tmp_path)//'/block_nodes_in_z')) then
             lagrangian_multi=lagrangian_multi+1
             blocking(3)=.true.
             ewrite(3,*) "No grid movement in Z "
          endif

          if (lagrangian_multi>=mesh_dim(mesh)) then 
             ewrite(0,*) "The grid is blocked in all direction!"
             ewrite(0,*) "Do you know what you're doing?"
          endif
       endif
    endif


    if (has_vector_field(state,"GridVelocity_Yves"))then
       gridvelocity_yves=extract_vector_field(state,"GridVelocity_Yves")
    endif

    ifield=>extract_scalar_field(state,fieldname)
    ewrite(3,*) "will follow this field: ",fieldname

    if (.not. ifield%mesh==positions%mesh) then
       ewrite(3,*) "mapping input field if DG"
       call allocate(linear_ifield, mesh)
       call zero(linear_ifield)
       call allocate(lumpedmass, mesh, "LumpedMass")
       call compute_lumped_mass(positions, lumpedmass)
       ! Invert lumped mass.
       lumpedmass%val=1./lumpedmass%val
       P=compute_projection_matrix(mesh, velocity%mesh, positions)
  
       call mult(linear_ifield,P,ifield)
       ! Apply inverted lumped mass to projected quantity.
       call scale(linear_ifield, lumpedmass)
       !call remap_field(ifield,linear_ifield)
       ifield => linear_ifield
    end if
    call set_debug_level(3)

    if (have_option(trim(gridvelocity_yves%option_path)//"/prescribed")) then
       ewrite(3,*) "The grid velocity is prescribed!"
       iterating=.false.
    endif




    ewrite(3,*) " size of different mesh (just checking)"
    ewrite(3,*) " size of positions field",node_count(positions)
    ewrite(3,*) " size of Velocity field",node_count(velocity)
    ewrite(3,*) " size of grid velocity field",node_count(gridvelocity)
    ewrite(3,*) " size of temperature field",node_count(ifield)

    cdim=2*mesh_dim(mesh)-1

    call field_stats(gridvelocity,positions, mini, maxi, l2norm2)
    ewrite(3,*) "##grid_vel at the start of routine :mini, max,norm2 ",mini, maxi,l2norm2

    cdim=mesh_dim(mesh)+lagrangian_multi

    ewrite(3,*) "we're gonna solve a system of dimension: ", cdim
    !if blcok5x5 was a block matrix...
    !call allocate(block5x5, u_sparsity, (/cdim, cdim/), diagonal=diag,&
    !& name="block5x5")
    !call zero(block5x5)

    MaxX=0.0
    MinX=0.0

    call allocate(masslump,ifield%mesh)
    call zero(masslump)

    nnlist2 => extract_nnlist(mesh)

    if(scaling)then
       ewrite(3,*) "## scaling of the functional is untested :"
       call field_stats(ifield,positions, mini, maxi, l2norm2)
       rhorange=maxi-mini
       do i=1,velocity%dim
          field1=extract_scalar_field_from_vector_field(positions,i)
          call field_stats(field1,positions, mini, maxi, l2norm2)
          MaxX(i)=maxi
          MinX(i)=mini
          rangeX(i)=maxi-mini
          !this is not good...at all  
         
          if(MeshCount.ge.1)then
             field1=extract_scalar_field_from_vector_field(velocity,i)
             call field_stats(field1,positions, mini, maxi, l2norm2)
             speedrange(i)=maxi-mini
          endif
          if (option==1) F1weight(i)= 1/(rangeX(1)* rangeX(2)* rangeX(3)*(0.5*speedrange(i)*rhorange/rangeX(i))**2)
          if (option==2) F2weight(i)=1/rhorange**2
          if(spring)then
             DeltaX(i)=rangeX(i)/Nx(i)
             Fspringweight(i)=1/ DeltaX(i)**2
          endif
       end do
       F1weight=F1weight(3)
    endif

    if (spring.or.expo) then
       call allocate(diagK,ifield%mesh, "Springmass")
       call allocate(Fspring,ifield%mesh, "FSpring")
       call zero(Fspring)
       call zero(DiagK)
       if (expo) then
          call allocate(Fe,ifield%mesh, "Exponential term")
          call zero(Fe)
       endif
       D0=0.01
       do node=1,size(Nnlist2,1)
          call set(diagK,node,2.0*(Nnlist2%findrm(node+1)-1-Nnlist2%findrm(node)))
          beta=0.1
          A0=1.0
          if (.false.) then
             beta=0.1
             A0=1.0
             Dij=0.0
             do j=NNList2%findrm(node),NNList2%findrm(node+1)-1
                node2=NNList2%colm(j)
                X2=node_val(positions,node)-node_val(positions,node2)
                ! Dij=X2(1)*X2(1)+X2(2)*X2(2)+X2(3)*X2(3)
                !Dij=Dij+X2(3)*X2(3)/(DeltaX(3)*DeltaX(3))
                Dij=Dij+X2(3)
             end do

             if (Dij.gt.0.0000001)then
                call addto(Fe,node,1/(1-Dij/D0)**1/3)
             endif
          endif
       end do
       !       call vtk_write_fields(trim("functionals_start"), MeshCount, positions, mesh, sfields=(/Fe/), vfields=(/gridvelocity/))

    endif

    if(option==1) then
       U_sparsity => get_csr_sparsity_firstorder(state, velocity%mesh, velocity%mesh)
       call allocate(MatrixA,U_sparsity,(/mesh_dim(ifield),1/))
       call allocate(MatrixAlumped,mesh_dim(ifield), ifield%mesh)
       call allocate(F1, ifield%mesh,"F1field")

       call zero(F1)
       call zero(MatrixAlumped)
       call zero(MatrixA)
    endif

    do ele=1, element_count(ifield)
       if(option==1) then
          call  assemble_element_contribution(Massmatrix=masslump, matrixA=MatrixA,&
               positions=positions,ifield=ifield,ele= ele)
       else
          call  assemble_element_contribution(Massmatrix=masslump, &
               positions=positions, ifield=ifield, ele=ele)
       end if
    end do
  

    if(option==1) then
       do d=1, positions%dim
          !loop block by block ...this will be fix soon
          Ablock=block(MatrixA,d,1)
          do  i =1 ,size(Ablock,1)!node_count(MatrixAlumped)
             MatrixAlumped%val(d,i)=sum(row_val_ptr(Ablock,i))
          end do
       end do
       !       call deallocate(Ablock,stat=lstat)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!OPTION 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(option==2) then
       call allocate(hessian,ifield%mesh)
       call reset_debug_level
       call compute_hessian(ifield, positions, hessian)
       call set_debug_level(3)
       allocate(evalues(mdim))
       allocate(evectors(mdim,mdim))
       do node=1,node_count(hessian)
          evectors=0.0
          evalues=0.0
          call eigendecomposition_symmetric(node_val(hessian,node), evectors, evalues)
          evalues=abs(evalues)
          !hessian%val(:, :, node)=Mat_Diag_Mat(evectors,evalues)
          call set(hessian,node,Mat_Diag_Mat(evectors,evalues))
          if (scaling) then
             do i=1,mdim
                F2weight(i)=max(F2weight(i),evalues(i))
             end do
          endif
       end do
       deallocate(evalues)
       deallocate(evectors)
    endif

    call set_debug_level(3)
    allocate(ugnew(node_count(ifield),cdim))
    allocate(ugold(node_count(ifield),cdim))
    allocate(Unplus1(cdim))
    allocate(B(node_count(ifield),cdim))
    allocate(BNEW(node_count(ifield),cdim))
    allocate(BLOCK5X5(node_count(ifield),cdim,cdim))
    allocate(BLOCK5X5L(node_count(ifield),cdim,cdim))
    allocate(BVEC(node_count(ifield),cdim))
    allocate(BFUNC(node_count(ifield),mdim))
    allocate(A(cdim,cdim))
    allocate(Bn(cdim,1))

    BVEC=0.0
    B=0.0
    BLOCK5X5=0.0
    BLOCK5X5L=0.0
    alpha=0.001
    LAMXX=0.02
    LAM=LAMXX
    EPSILON=0.0000000001
    crazynumber=10000000000000000000000.0
   
    if(option==2) then
       do node=1,node_count(ifield)
          do i=1,mdim
             do j=1,mdim
                BLOCK5X5(node,i,j)=(1/F2weight(i))*hessian%val(i,j,node)
             end do
             BLOCK5X5(node,i,i)= BLOCK5X5(node,i,i)&
                  +LAM(i)*node_val(masslump,node)
          end do
       end do

       call apply_lagrangian_multiplier(BLOCK5X5,cdim,blocking)  

       do node=1,node_count(ifield)
          do i=1,mdim
             speed=node_val(velocity,node)
             do j=1,mdim
                B(node,i)=B(node,i)+BLOCK5X5(node,i,j)*speed(j)
             end do
          end do

          if(spring) then
             do k=1,mdim
                BLOCK5X5(node,k,k)= BLOCK5X5(node,k,k)+alpha*Fspringweight(k)*node_val(diagK,node)
             end do
          endif
       end do
    endif

    if(option==1)then
       call apply_lagrangian_multiplier(BLOCK5X5,cdim,blocking)
       do node=1, node_count(ifield)
          if(.false.)then
             if (cdim==4)then
                if(have_option(trim(tmp_path)//'/block_nodes_in_x')) then
                   BLOCK5X5(node,4,1)=1.0
                   BLOCK5X5(node,1,4)=1.0
                elseif(have_option(trim(tmp_path)//'/block_nodes_in_y')) then
                   BLOCK5X5(node,2,4)=1.0
                   BLOCK5X5(node,4,2)=1.0
                elseif(have_option(trim(tmp_path)//'/block_nodes_in_z')) then
                   BLOCK5X5(node,3,4)=1.0
                   BLOCK5X5(node,4,3)=1.0
                endif
             elseif (cdim==5) then
                if(have_option(trim(tmp_path)//'/block_nodes_in_x')) then
                   BLOCK5X5(node,4,1)=1.0
                   BLOCK5X5(node,1,4)=1.0
                endif
                if(have_option(trim(tmp_path)//'/block_nodes_in_y').and. &
                     have_option(trim(tmp_path)//'/block_nodes_in_x')) then
                   BLOCK5X5(node,2,5)=1.0
                   BLOCK5X5(node,5,2)=1.0
                elseif(have_option(trim(tmp_path)//'/block_nodes_in_y')) then
                   BLOCK5X5(node,2,4)=1.0
                   BLOCK5X5(node,4,2)=1.0
                endif
                if(have_option(trim(tmp_path)//'/block_nodes_in_z')) then
                   BLOCK5X5(node,3,5)=1.0
                   BLOCK5X5(node,5,3)=1.0
                endif
             endif
          endif

          do i=1,mdim
             BLOCK5X5(node,i,i)=F1weight(i)*node_val(MatrixAlumped,i,node)*node_val(MatrixAlumped,i,node)&
                  /node_val(masslump,node)+LAMXX*node_val(masslump,node)
          end do
          if(spring) then
             do i=1,mdim
                BLOCK5X5(node,i,i)= BLOCK5X5(node,i,i)+alpha*Fspringweight(i)*node_val(diagK,node)
             end do
          endif
       end do
       !asembl The RHS vector B
       call compute_BVEC(matrixA=MatrixA,MassLump=masslump,velocity=velocity,BVEC=B,weights=F1weight,LAM=LAM)
    endif



    ugnew=0.0
    ugold=0.0
    iteration=10
    if((option==2).and.(.NOT.(spring))) iteration=1

    !if (mdim==3) newgridvelocity=wrap_vector_field(mesh,ugnew(:,1),ugnew(:,2),ugnew(:,3),"Grid VelocityTemp")
    !if (mdim==2) newgridvelocity=wrap_vector_field(mesh,ugnew(:,1),ugnew(:,2),name="Grid VelocityTemp")


    call allocate(relativevelocity,velocity%dim,velocity%mesh,"RelativeVelocity")
    nonlineariteration=1
    if (iterating)then
      
       if (expo) then
          nonlineariteration=4
          call allocate(hessianFe,ifield%mesh)
          call zero(hessianFe)
          call allocate(gradfe,mdim,ifield%mesh)
          call zero(gradfe)
          do i=1,mdim
             call allocate(DFe(i),ifield%mesh)
             call zero(DFe(i))
             derivatives(i)=.true.
          end do
          allocate(evalues(mdim))
          allocate(evectors(mdim,mdim))
       endif


       do noniti=1,nonlineariteration
          BLOCK5X5L(:,:,:)=BLOCK5X5(:,:,:)
          !!Add constraints to make UG=0 on the boundary
          !          do d=1, velocity%dim
          !             do i=1, size(velocity%dirichlet(d)%node)
          !                inod=velocity%dirichlet(d)%node(i)
          !                BLOCK5X5L(inod,d,d)=crazynumber
          !             end do
          !          end do

          !Add crazynumber to bottom ....this is gonna be fixed soon
          Xpos=0.0
          do i=1, node_count(ifield)
             Xpos=node_val(positions,i)
             do j=1,mdim
                if((Xpos(j).lt.(MinX(j)+epsilon)).or.(Xpos(j).gt.(MaxX(j)-epsilon))) then
                   BLOCK5X5L(i,j,j)=crazynumber
                endif
             end do
          end do


          do it=1,iteration

             if ((expo).and.(it==iteration))then

                do node=1, node_count(velocity)
                   !probably no need to loop here
                   call set (relativevelocity,node,node_val(velocity,node)-ugnew(node,:))
                   ! do d=1,velocity%dim
                   !    relativevelocity%val(d,node)=velocity%val(d,node)-newgridvelocity%val(d,node)
                   ! end do
                end do

                do node=1,size(Nnlist2,1)
                   Dij=0.0
                   minidij=1.0
                   do j=NNList2%findrm(node),NNList2%findrm(node+1)-1
                      node2=NNList2%colm(j)
                      X2=ugnew(node,:)-ugnew(node2,:)
                      Pos2=node_val(positions,node)-node_val(positions,node2)   
                      Dij=Pos2(3)+DT*X2(3)    
                      call addto(Fe,node,1/(1-Dij/D0)**1/3)                     
                      ! Fe%val(node)=Fe%val(node)+1/(1-Dij/D0)**1/3 
                   end do
                end do

                call field_stats(Fe,positions, mini, maxi, l2norm2)
                Feweight(3)=1.0!/(maxi-mini)
                Feweight(1)=0.0
                Feweight(2)=0.0
                !FeweightScal=1.0/(maxi-mini)
                FeweightScal=1.0

                !calculate hessian of Fe
                call compute_hessian(Fe, positions, hessianFe)
                !do node=1, hessianFe%mesh%nodes
                !   evectors=0.0
                !   evalues=0.0
                !   call eigendecomposition_symmetric(node_val(hessianFe,node), evectors, evalues)
                !evalues=abs(evalues)
                !   hessianFe%val(:, :, node)=Mat_Diag_Mat(evectors,evalues)
                !end do
                !Derivative of Fe (also i should do it as I compute the hessian
                !                call differentiate_field(Fe,positions,derivatives,DFe)
                call grad(FE,positions,gradfe)
             endif

             ugold=ugnew
             BVEC=0.0
             BFUNC=0.0
             BNEW=0.0
             
             call zero(relativevelocity)
             if (option==1) then
                call allocate(newgridvelocity,3,gridvelocity%mesh,"NewGridVelocity")
                call set_all(newgridvelocity, transpose(ugnew))
                call compute_BVEC(matrixA=MatrixA,MassLump=masslump,velocity=newgridvelocity,BVEC=BVEC,&
                     weights=F1weight,LAM=LAM)
                call deallocate(newgridvelocity)
                do node=1, node_count(velocity)
                   call set (relativevelocity,node,node_val(velocity,node)-ugnew(node,:))
                   !do d=1,velocity%dim
                   !   relativevelocity%val(d,node)=abs(velocity%val(d,node)-newgridvelocity%val(d,node))
                   !end do
                end do

                call compute_BVEC(matrixA=MatrixA,MassLump=masslump,velocity=relativevelocity,BVEC=BFUNC,&
                     weights=F1weight,LAM=LAM)

                do node=1,node_count(velocity)
                   !call addto(F1,node,abs(BFUNC(node,i)*node_val(relativevelocity))
                   do i=1,mdim
                      F1%val(node)= F1%val(node)+abs(BFUNC(node,i)*relativevelocity%val(i,node))
                   end do
                end do

                ewrite(3,*) "if option 1"
             endif
             if(spring) then
                call allocate(newgridvelocity,3,gridvelocity%mesh,"NewGridVelocity")
                call set_all(newgridvelocity, transpose(ugnew))
                call compute_BVEC(velocity=newgridvelocity,BVEC=BVEC,Nnlist2=Nnlist2,alpha=alpha,weights=Fspringweight)
                call deallocate(newgridvelocity)
                ewrite(3,*) "if spring",spring
             endif
             !if(spring) call compute_BVEC(velocity=newgridvelocity,BVEC=BVEC,Nnlist2=Nnlist2,alpha=alpha,weights=Fspringweight)

             do node=1, node_count(ifield)
                Unplus1=0.0

                !!add the constraint to BVEC as compute_BVEC only consider a 3X3 system
                do i=1,cdim
                   do j=1,cdim
                      if((i.gt.mdim).or.(j.gt.mdim)) then
                         BVEC(node,i)=BVEC(node,i)+BLOCK5X5(node,i,j)*ugnew(node,j)
                      endif
                   end do
                end do

                if ((option==2).and.(spring)) then
                   do k=1,mdim
                      BNEW(node,k)=B(node,k)+alpha*Fspringweight(k)*node_val(diagK,node)*ugnew(node,k)-BVEC(node,k)
                   end do
                else

                   Unplus1=matmul(BLOCK5X5(node,:,:),ugnew(node,:))
                   !sort that out later for the moment do it by hand
                   ! do i=1,cdim
                   !    do j=1,cdim
                   !       Unplus1(i)=Unplus1(i)+BLOCK5X5(node,i,j)*ugnew(node,j)
                   !    end do
                   ! end do
                   !Unplus1=0.0 
                   BNEW(node,:)=B(node,:)+Unplus1(:)-BVEC(node,:)

                endif
               
                if ((expo).and.(it==iteration))then
                   !add Hessian to LHS
                   AHFe=node_val(hessianFe,node)
                   do j=1,mdim
                      !BLOCK5X5L(node,1:3,1:3)= BLOCK5X5L(node,1:3,1:3)+FeweightScal*node_val(hessianFe,node)
                      BLOCK5X5L(node,mdim,j)= BLOCK5X5L(node,mdim,j)+FeweightScal*AHFe(mdim,j)

                      !add hessian*u^n -Fe-derivative(U^n) to RHS
                      !BNEW(NODE,:)=BNEW(NODE,:)+FeweightScal*matmul(node_val(hessianFe,node),ugnew(node,1:3))
                      BNEW(NODE,mdim)=BNEW(NODE,mdim)+FeweightScal*AHFe(mdim,j)*ugnew(node,j)
                   end do
                   sp=0.0
                   !do i=1,3
                   !   BNEW(NODE,i)=BNEW(NODE,i)-FeweightScal*node_val(DFe(i),node)
                   !   sp=sp+gradfe%val(i,node)*positions%val(i,node)
                   !end do

                   !BNEW(NODE,3)=BNEW(NODE,3)-FeweightScal*gradfe%val(3,node)

                   !  end do
                endif

                !A=BLOCK5X5L(node,:,:)
                A=BLOCK5X5(node,:,:)
                Bn(:,1)=BNEW(node,:)

                call solve(A,Bn,stat=lstat)
                ugnew(node,:)=Bn(:,1)
                
             end do
          end do
       end do !non-linear iteration expo
    endif

    if(expo)then
       deallocate(evalues)
       deallocate(evectors)
       call deallocate(hessianFe)
       call deallocate(DFe(1))
       call deallocate(DFe(2))
       call deallocate(DFe(3))
    endif

    call zero(gridvelocity)
    ewrite(3,*) "Setting the grid velocity ",gridvelocity%dim
    call set_all(gridvelocity, transpose(ugnew))
    ewrite(3,*) "Done!" 

    if (.not.iterating) then
       call zero (gridvelocity)
       pi=acos(0.0)
       ewrite(3,*) "Setting the grid Velocity to grid_velocity_yves" 
       if ((cos(MeshCount*pi/6.0)*sin(MeshCount*pi/6.0)).gt.0.0) then
          call set(gridvelocity,gridvelocity_yves)    
          ewrite(3,*) "Setting the grid Velocity positive " ,MeshCount*dt
          call field_stats(gridvelocity,positions, mini, maxi, l2norm2)
          ewrite(3,*) "##grid_vel:mini, max,norm2",mini, maxi,l2norm2
       endif
       if ((cos(MeshCount*pi/6.0)*sin(MeshCount*pi/6.0)).lt.0.0) then
          call addto(gridvelocity,gridvelocity_yves,scale=-1.0)
          !do node=1, node_count(velocity)
          !   call set (gridvelocity,node,-node_val(gridvelocity_yves,node))
          !end do

          ewrite(3,*) "Setting the grid Velocity negative " ,MeshCount*dt
          call field_stats(gridvelocity,positions, mini, maxi, l2norm2)
          ewrite(3,*) "##grid_vel:mini, max,norm2",mini, maxi,l2norm2
       endif
    endif
    ewrite(3,*) "Starting to deallocate things" 
    deallocate(Unplus1)
    deallocate(A)
    deallocate(Bn)
    if (option==2)  call deallocate(hessian)
    call deallocate(masslump)
    
    call field_stats(gridvelocity,positions, mini, maxi, l2norm2)
    ewrite(3,*) "##grid_vel:mini, max,norm2",mini, maxi,l2norm2


    call move_internal_nodes(gridvelocity,positions,Nnlist2,dt)


    do node=1, node_count(velocity)
       call set (relativevelocity,node,node_val(velocity,node)-node_val(gridvelocity,node))
    end do


    ewrite(3,*) "before if spring end"

    if (spring) then 
       call zero(Fspring)
       do i=1,node_count(velocity)
          tempr=0.0

          ! Fe%val(i)=Fe%val(i)*FeweightScal
          !Fspring
          speed=node_val(relativevelocity,i)
          do j=Nnlist2%findrm(i),Nnlist2%findrm(i+1)-1
             node2=Nnlist2%colm(j)
             speed2=node_val(relativevelocity,node2)   
             tempr=tempr+norm2(speed-speed2)*norm2(speed-speed2)
          end do
          !F1
          Fspring%val(i)=tempr
       end do
    endif
    ewrite(3,*) "before vtk output"
    if (expo) then
       ! call vtk_write_fields(trim("functionals"), MeshCount, positions, mesh, sfields=(/Fe,F1,Fspring/), vfields=(/gridvelocity,relativevelocity/))
    elseif(spring)then
       !call vtk_write_fields(trim("functionals"), MeshCount, positions, mesh, sfields=(/F1,Fspring/), vfields=(/gridvelocity,relativevelocity/))
    else 
       ! call vtk_write_fields(trim("functionals"), MeshCount, positions, mesh, sfields=(/F1/), vfields=(/gridvelocity,relativevelocity/))
    endif

    ewrite(3,*) "update positions"
    ewrite(3,*) "deallocatting vectors..."
    if (spring) then
       call deallocate(DiagK)
       call deallocate(Fspring)
       ewrite(3,*) "deallocated spring..."
    endif
    if (expo) call deallocate(Fe)
    if (option==1) then
       ewrite(3,*) "deallocating matrixalumped..."
       call deallocate(MatrixAlumped)
       ewrite(3,*) "deallocatting F1..."
       call deallocate(F1)
       if (MatrixA%clone) then
          ewrite(3,*) "MatrixA is cloned" 
       endif
       if (MatrixA%external_val)  then
          ewrite(3,*) "MatrixA has external val"
       endif
       if ((associated(MatrixA%val))) then
          ewrite(3,*) "MatrixA val is associated"
       endif
       call deallocate(MatrixA,stat=lstat)
    endif

    ewrite(3,*) "deallocatting the rest..."
    call deallocate(relativevelocity)
    deallocate(ugnew)
    deallocate(ugold)
    deallocate(B)
    deallocate(BNEW)
    deallocate(BVEC)
    deallocate(BFUNC)
    deallocate(BLOCK5X5)
    deallocate(BLOCK5X5L)

    call field_stats(gridvelocity,positions, mini, maxi, l2norm2)
    ewrite(3,*) "##grid_vel at end of routine :mini, max,norm2",mini, maxi,l2norm2
    if (associated(gridvelocity, linear_gridvelocity)) then
       ! extract the actual grid velocity
       gridvelocity => extract_vector_field(state, "IteratedGridVelocity")
       ! map the grid velocity computed on the coordinate mesh back to it
       call remap_field(linear_gridvelocity, gridvelocity)
       call deallocate(linear_gridvelocity)
    end if
    if (associated(velocity, linear_velocity)) then
       call deallocate(linear_velocity)
    end if
    call set_vector_field_in_state(state,"IteratedGridVelocity","IteratedGridVelocity")

    MeshCount=MeshCount+1
    ewrite(3,*) "end subroutine", MeshCount
    call reset_debug_level
  end subroutine movemeshy

  subroutine compute_BVEC(matrixA,MassLump,Cmatrix,velocity,BVEC,Nnlist2,alpha,weights,LAM)
    real,dimension(:,:,:),intent(in),optional::Cmatrix 
    type(block_csr_matrix), intent(in), optional :: matrixA
    real,dimension(:,:),intent(inout)::BVEC
    type(vector_field),intent(in)::velocity
    type(scalar_field),intent(in),optional::MassLump
    type(csr_sparsity),intent(in),optional:: Nnlist2
    real,intent(in),optional::alpha
    integer node1,node2,i,j,d,lstat,mdim
    real,dimension(:),allocatable::speed,speed2,bvector
    real,dimension(:,:),allocatable:: blockm
    type(csr_matrix):: matrixAblock

    real,dimension(:),allocatable::Btemp
    real,dimension(:),allocatable::Btemp2
    real,dimension(:),allocatable::temp1

    real,dimension(:),intent(in),optional::LAM
    real,dimension(:),intent(in),optional::weights
    
    mdim=velocity%dim

    allocate(speed(mdim))
    allocate(speed2(mdim))
    allocate(bvector(mdim))
    allocate(blockm(mdim,mdim))

    if (present(matrixA)) then
       allocate(Btemp(MassLump%mesh%nodes))
       allocate(Btemp2(MassLump%mesh%nodes))
       allocate(temp1(MassLump%mesh%nodes))

       Btemp=0.0
       Btemp2=0.0
       temp1=0.0
       ! call mult(Btempfield,matrixA,velocity)
       do i=1,velocity%dim
          matrixAblock=block(matrixA,i,1)
          !   field1=extract_scalar_field_from_vector_field(velocity,i)
          temp1=velocity%val(i,:)
          call mult(Btemp,matrixAblock,temp1)
          call mult_T(Btemp2,matrixAblock,Btemp)

          do node1=1, velocity%mesh%nodes
             BVEC(node1,i)=BVEC(node1,i)+weights(i)*Btemp2(node1)/node_val(MassLump,node1)&
                  +LAM(i)*node_val(MassLump,node1)*temp1(node1)   
          end do
       end do

       deallocate(Btemp)
       deallocate(Btemp2)
       deallocate(temp1)

    elseif(present(Cmatrix)) then
       do node1=1, velocity%mesh%nodes
          bvector=0.0
          speed=node_val(velocity,node1)
          do i=1,velocity%dim
             do j=1,velocity%dim
                blockm(i,j)=Cmatrix(node1,i,j)
             end do
          end do
          bvector=matmul(blockm,speed)
          do d=1,velocity%dim
             BVEC(node1,d)=BVEC(node1,d)+bvector(d)+LAM(d)*MassLump%val(node1)
          end do
       end do

    elseif(present(Nnlist2)) then
       do node1=1,size(Nnlist2,1)
          bvector=0.0
          speed=node_val(velocity,node1)

          do j=Nnlist2%findrm(node1),Nnlist2%findrm(node1+1)-1
             node2=Nnlist2%colm(j)
             speed2=node_val(velocity,node2)   
             bvector=bvector+2*(speed-speed2)
          end do
          do d=1,velocity%dim
             BVEC(node1,d)=BVEC(node1,d)+alpha*weights(d)*bvector(d)
          end do
       end do
    endif

    deallocate(speed)
    deallocate(speed2)
    deallocate(bvector)
    deallocate(blockm)

  end subroutine compute_BVEC

  subroutine assemble_element_contribution(Massmatrix, matrixA,&
       &  positions, ifield, ele)
    type(scalar_field), intent(inout) :: Massmatrix
    type(block_csr_matrix), intent(inout), optional :: matrixA
    !type(vector_field), intent(inout), optional :: matrixAlumped
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: ifield
    integer, intent(in) :: ele
    real,dimension(:,:,:),allocatable::di_t
    real, dimension(:),allocatable::detwei
    real,dimension(:,:),allocatable::mass
    integer i
    integer, dimension(:), pointer :: i_ele
    type(element_type), pointer :: shape_field, shape_X

    allocate(di_t(ele_loc(ifield,ele),ele_ngi(ifield,ele),mesh_dim(ifield)))
    allocate(detwei(ele_ngi(ifield,ele)))
    allocate(mass(positions%mesh%shape%loc,positions%mesh%shape%loc))
    shape_field=>ele_shape(ifield, ele)
    shape_X=>ele_shape(positions, ele)
    i_ele=>ele_nodes(ifield,ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_field,dshape=di_t&
         &,detwei=detwei)
    mass=shape_shape(shape_X, shape_X, detwei)

    call addto(Massmatrix, i_ele, sum(mass,2))

    if (present(matrixA)) then

       do i=1,mesh_dim(ifield)
          call addto(matrixA, i, 1, i_ele, i_ele, &
               shape_shape(shape_X, shape_X, detwei* &
               &             matmul(ele_val(ifield,ele), di_t(:,:,i))))
       end do
    end if

    deallocate(di_t)
    deallocate(detwei)
    deallocate(mass)

  end subroutine assemble_element_contribution


  subroutine move_internal_nodes(MeshVelocity,positions,Nnlist2,dt,move)

    INTEGER inod,j,node2,i
    REAL, intent(in)::DT
    Real dist1,dist2
    type(csr_sparsity),intent(in):: Nnlist2
    ! INTEGER CORNERS(8),corn
    type(vector_field),intent(inout)::  MeshVelocity
    type(vector_field),intent(inout)::  positions
    integer,intent(inout),optional:: move
    Real, dimension(:),allocatable::NewX,X,Speed,NewX2,X2,Speed2,MinX,MaxX
    !Real, dimension(3)::Speed
    real epsilon,mini,maxi,l2norm
    type(scalar_field)::field1
    REAL DThalf

    allocate(NewX(positions%dim))
    allocate(X(positions%dim))
    allocate(Speed(positions%dim))
    allocate(NewX2(positions%dim))
    allocate(X2(positions%dim))
    allocate(Speed2(positions%dim))
    allocate(MinX(positions%dim))
    allocate(MaxX(positions%dim))
    do i=1,positions%dim
       field1=extract_scalar_field_from_vector_field(positions,i)
       call field_stats(field1,positions,mini,maxi,l2norm)
       MaxX(i)=maxi
       MinX(i)=mini
    end do
    call field_stats(MeshVelocity,positions, mini, maxi, l2norm)
    ewrite(3,*) "##grid_vel in move_internal_nodes :mini, max,norm2",mini, maxi,l2norm
    !call  get_time_step(DT)
    DThalf=DT
    !DT=0.01
    ! REAL MAXX,MAXZ,MAXY,MINX,MINZ,MINY
    !ewrite(3,*) "DT", DT, MaxX(1),MaxX(2),MaxX(3), MinX(1), MinX(2), MinX(3)
    ewrite(3,*) "sizes....",node_count(positions),node_count(MeshVelocity)
    epsilon =1.0E-8

    do inod=1, node_count(positions)
       X=node_val(positions,inod)
       Speed=node_val(MeshVelocity,inod)
       NewX=X+DThalf*Speed

       !Don't move the boundaries
       !allow to move only along the boundary lines
       do i=1,positions%dim
          !  if (abs(Speed(i)).LT.0.000001)then
          !     NewX(i)=X(i)
          !  else
          if ((X(i).LT.(MinX(i)+epsilon)).or.(X(i).GT.(MaxX(i)-epsilon))) then
             NewX(i)=X(i)
             Speed(i)=0.0
             MeshVelocity%val(i,inod)=Speed(i)
          endif
          if ((NewX(i).LT.MaxX(i)).and.(NewX(i).GT.MinX(i))) then
             positions%val(i,inod)=NewX(i)
             MeshVelocity%val(i,inod)=(NewX(i)-X(i))/DThalf
          else
             NewX(i)=X(i)
             Speed(i)=0.0
             MeshVelocity%val(i,inod)=Speed(i)
          endif
          !  endif
       end do
    end do

    deallocate(NewX)
    deallocate(X)
    deallocate(Speed)
    deallocate(NewX2)
    deallocate(X2)
    deallocate(Speed2)
    deallocate(MinX)
    deallocate(MaxX)

    call field_stats(MeshVelocity,positions, mini, maxi, l2norm)
    ewrite(3,*) "##grid_vel at end of routine :mini, max,norm2",mini, maxi,l2norm


    ewrite(3,*) "end move_internal_nodes"

  end subroutine move_internal_nodes

  subroutine apply_lagrangian_multiplier(matrix,cdim,blocking)


    integer, intent(in)::cdim
    real, dimension(:,:,:)::matrix
    logical,dimension(3),intent(in)::blocking
    logical X,Y,Z

    X=blocking(1)  
    Y=blocking(2)  
    Z=blocking(3)  
    if (cdim==4.or.cdim==3) then
       if(X) then
          matrix(:,cdim,1)=1.0
          matrix(:,1,cdim)=1.0
       elseif(Y) then
          matrix(:,cdim,2)=1.0
          matrix(:,2,cdim)=1.0
       elseif (Z) then
          matrix(:,cdim,3)=1.0
          matrix(:,3,cdim)=1.0
       endif
    elseif (cdim==5) then
       if(X) then
          matrix(:,4,1)=1.0
          matrix(:,1,4)=1.0
          if(Y) then 
             matrix(:,5,2)=1.0
             matrix(:,2,5)=1.0
          elseif(Z) then
             matrix(:,5,3)=1.0
             matrix(:,3,5)=1.0
          endif
       else
          matrix(:,4,2)=1.0
          matrix(:,2,4)=1.0
          matrix(:,5,3)=1.0
          matrix(:,3,5)=1.0
       endif
    endif
  end subroutine apply_lagrangian_multiplier

  subroutine move_mesh_imposed_velocity(states)
  type(state_type), dimension(:), intent(inout) :: states
  
    type(vector_field), pointer :: coordinate, old_coordinate, new_coordinate
    type(vector_field), pointer :: velocity
    type(vector_field), pointer :: grid_velocity
    
    integer :: i, stat
    real :: itheta, dt
    logical :: found_velocity
    
    if(.not.have_option("/mesh_adaptivity/mesh_movement/imposed_grid_velocity")) return
    call IncrementEventCounter(EVENT_MESH_MOVEMENT)
    
    ewrite(1,*) 'Entering move_mesh_imposed_velocity'
    
    grid_velocity => extract_vector_field(states(1), "GridVelocity")
    
    coordinate => extract_vector_field(states(1), "Coordinate")
    old_coordinate => extract_vector_field(states(1), "OldCoordinate")
    new_coordinate => extract_vector_field(states(1), "IteratedCoordinate")
    
    call get_option("/timestepping/timestep", dt)
    
    found_velocity = .false.
    do i = 1, size(states)
      velocity => extract_vector_field(states(i), "Velocity", stat)
      if(stat==0 .and. .not. velocity%aliased) then
        call get_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/relaxation", itheta, stat)
        if(found_velocity.and.(stat==0)) then
          FLExit("Only one prognostic velocity allowed with imposed mesh movement.")
        else
          found_velocity = (stat==0)
        end if
      end if
    end do
    if(.not.found_velocity) then
      itheta = 0.5
    end if
    
    call set(new_coordinate, old_coordinate)
    call addto(new_coordinate, grid_velocity, scale=dt)
    
    call set(coordinate, new_coordinate, old_coordinate, itheta)
  
  end subroutine move_mesh_imposed_velocity

  subroutine move_mesh_pseudo_lagrangian(states)
  type(state_type), dimension(:), intent(inout) :: states
  
    type(vector_field), pointer :: coordinate, old_coordinate, new_coordinate
    type(vector_field), pointer :: velocity
    type(vector_field), pointer :: grid_velocity
    
    integer :: i, stat
    real :: itheta, dt
    logical :: found_velocity
    
    character(len=FIELD_NAME_LEN) :: state_name
    
    if(.not.have_option("/mesh_adaptivity/mesh_movement/pseudo_lagrangian")) return
    call IncrementEventCounter(EVENT_MESH_MOVEMENT)
    
    ewrite(1,*) 'Entering move_mesh_pseudo_lagrangian'
    
    grid_velocity => extract_vector_field(states(1), "GridVelocity")
    
    call get_option("/mesh_adaptivity/mesh_movement/pseudo_lagrangian/velocity_material_phase/material_phase_name", &
                    state_name, stat=stat)
    if(stat==0) then
      i = get_state_index(states, trim(state_name))
      velocity => extract_vector_field(states(i), "Velocity")
    else
      velocity => extract_vector_field(states(1), "Velocity")
    end if
    
    call set(grid_velocity, velocity)
    
    coordinate => extract_vector_field(states(1), "Coordinate")
    old_coordinate => extract_vector_field(states(1), "OldCoordinate")
    new_coordinate => extract_vector_field(states(1), "IteratedCoordinate")
    
    call get_option("/timestepping/timestep", dt)
    
    found_velocity = .false.
    do i = 1, size(states)
      velocity => extract_vector_field(states(i), "Velocity", stat)
      if(stat==0) then
        call get_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/relaxation", itheta, stat)
        if(found_velocity.and.(stat==0)) then
          FLExit("Only one prognostic velocity allowed with pseudo lagrangian mesh movement.")
        else
          found_velocity = (stat==0)
        end if
      end if
    end do
    if(.not.found_velocity) then
      itheta = 0.5
    end if
    
    call set(new_coordinate, old_coordinate)
    call addto(new_coordinate, grid_velocity, scale=dt)
    
    call set(coordinate, new_coordinate, old_coordinate, itheta)
  
  end subroutine move_mesh_pseudo_lagrangian

end module meshmovement


