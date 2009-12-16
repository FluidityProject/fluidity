#include "fdebug.h"

program main

  !solver for the GN operator

  !coded by Colin Cotter and David Ham May-June 2006
  !for last modification see CVS

  !links against libdfluidity

  use elements
  use sparse_tools
  use quadrature
  use global_numbering
  use shape_functions
  use global_parameters_gallopede
  use adjacency_lists
  use transform_elements
  use dgtools
  use text_io
  use data_structures
  use Solvers
  use helmholtz
  use nc_tools
  use parallel_tools

  implicit none

  !local variables

  real, allocatable, target::X(:),Y(:) !vertex coordinates
  real, allocatable, dimension(:), target ::u,m, rhs
  !velocity vector storing all velocities
  real, pointer, dimension(:) :: u1, u2, m1, m2, rhs1, rhs2
  !individual components
  real, allocatable::D(:) !height field
  integer, allocatable, target::EVList_h(:) !element-vertex list for height
  integer, allocatable, target::EVList_u(:) !element-vertex list for vels
  type(csr_matrix),target:: Mass_h,Mass_u !mass matrix
  type(block_csr_matrix) :: mommat !differentiation matrix
  type(csr_matrix) :: mom11, mom12, mom21, mom22 !components of diff mat
  type(element_type) :: nh, nh_f   ! Reference height element
  type(element_type) :: nu  ! Reference velocity element
  type(quadrature_type) :: g, g_f   ! Gauss quadrature
  type(csr_matrix),target ::EEList, B_SEG_list 
  integer, dimension(:), pointer :: EFlist
  ! Element-Element list, B_SEG list
  type(nc_mesh) :: mesh
  type(rotations) :: rots
  integer, dimension(:), allocatable, target :: b_seg_nh_lno
  integer :: i,globi,globj

  !==============================================
  !initial conditions

  !load mesh data
  MSG("Loading mesh data")
  call Get_Mesh(X,Y,EVList_h,EVList_u)

  MSG("allocating memory");
  allocate( u(N_vels*2) )
  u1 => u(1:N_vels)
  u2 => u(N_vels+1:2*N_vels)
  allocate( m(N_vels*2) )
  m1 => m(1:N_vels)
  m2 => m(N_vels+1:2*N_vels)
  allocate( rhs(N_vels*2) )
  rhs1 => rhs(1:N_vels)
  rhs2 => rhs(N_vels+1:2*N_vels)
  
  allocate( D(n_verts) )

  MSG("Loading initial conditions");

  call read_field(dr=m1,filename='m1.dat')
  call read_field(dr=m2,filename='m2.dat')
  call read_field(dr=D,filename='D_initial.dat')

  !==============================================
  !Initialisation of mesh data

  !get quadrature
  !the 2d elements
  MSG("Making 2D quadrature");
  g=make_quadrature(loc=3,dimension=2,degree=4)
  MSG("Making 1D quadrature");
  g_f=make_quadrature(loc=2,dimension=1,degree=3)

  !get_basis functions
  !2D

  MSG("Making 2D P1 elements");
  nh=make_element_shape(loc=3, dimension=2, degree=1, quad=g)
  MSG("Making 2D P1nc elements");
  nu=make_element_shape(loc=3, dimension=2, degree=1, quad=g, &
       type=ELEMENT_NONCONFORMING)

  !1D -- continuous basis functions
  MSG("Making 1D P1 elements");
  nh_f=make_element_shape(loc=2, dimension=1, degree=1, quad=g_f)

  !get Helmholtz matrix structure
  MSG("Making matrix structure");
  call MakeLists(N_Verts,N_Elements, Nloc, EVList_h, .false., EEList=EEList)
  CHECK(entries(EEList))

  MSG("Allocating b_seg_nh_lno");
  !stores local node numbers for h for face numbers
  allocate(b_seg_nh_lno(entries(EEList)*nh_f%loc))
  CHECK(size(b_seg_nh_lno))
  MSG("Making boundary_seg numbering");
  call make_boundary_seg_numbering_nc(b_seg_list, b_seg_nh_lno, & 
       EElist = EEList, nonod = N_Verts, ndglno = EVList_h, &
       Ele_n = nh, b_seg_n = nh_f)

  ! Set up Mass matrices
  call POSINM_DG_nc(Mass_u, N_elements, N_vels, nu%loc, evlist_u,&
       n_vels, nu%loc, evlist_u, eelist)

  MSG("setting up nc mesh data")
  mesh%n_elements = n_elements
  mesh%n_verts = n_verts
  mesh%n_vels = n_vels
  mesh%X => X
  mesh%Y => Y
  mesh%Evlist_u => Evlist_u
  mesh%Evlist_h => Evlist_h
  mesh%EElist => EElist
  mesh%mass_u => mass_u
  !mesh%mass_h => mass_h
  mesh%b_seg_list => b_seg_list
  mesh%b_seg_nh_lno => b_seg_nh_lno
  mesh%nu = nu
  mesh%nh = nh
  mesh%nh_f = nh_f
  MSG("cloning matrix")
  mommat = block_clone(mass_u,(/2,2/) )  

  call get_bcnormals(mesh,rots)

  MSG("Getting blocks")
  mommat = block_clone(mass_u,(/2,2/))
  mom11 = block(mommat,1,1)
  mom12 = block(mommat,1,2)
  mom21 = block(mommat,2,1)
  mom22 = block(mommat,2,2)

  call zero(mom11)
  call zero(mom12)
  call zero(mom21)
  call zero(mom22)

  u = 0.
  rhs = 0.

  MSG("call assemble_GN_operator")
  !assemble GN operator and rhs
  call assemble_GN_operator(mom11,mom12,mom21,mom22,rhs1,rhs2,D,m1,m2, &
       mesh)

  !CHECK(dense(mom11))
  !CHECK(dense(mom12))
  !CHECK(dense(mom21))
  !CHECK(dense(mom22))
  !stop

  MSG("solving GN operator equation")
  !solve GN equation
  call cg_solve_block_eqn(u, mommat, rhs,rots,1.0e-8, 500)

  !data output
  call dump_field(u1,N_vels,'u1_',1,1,9)
  call dump_field(u2,N_vels,'u2_',1,1,9)
  call dump_field(rhs1,N_vels,'rhs1_',1,1,11)
  call dump_field(rhs2,N_vels,'rhs2_',1,1,11)

  deallocate( u )
  deallocate( m )
  deallocate( rhs )
  deallocate( D )
  deallocate( b_seg_nh_lno )
  call deallocate( EEList )
  call deallocate( B_seg_List )

  contains

    subroutine assemble_GN_operator(mom11,mom12,mom21,mom22, &
         rhs1,rhs2,D,m1,m2, &
         mesh)
      implicit none
      type(csr_matrix), intent(out) :: mom11, mom12, mom21, mom22
      real, dimension(n_vels), intent(inout) :: rhs1, rhs2, m1, m2
      real, dimension(n_verts), intent(in) :: D
      type(nc_mesh) :: mesh

      !locals
      integer :: ele, globi, globj, iloc, jloc, gi,ni,ele_2,b_seg_i,bcnt
      integer :: i,j
      real :: kmat, kmat11, kmat12, kmat21, kmat22
      real, dimension(mesh%nu%ngi) :: dloc
      real, dimension(2,mesh%nh%loc) :: ele_X, ele_X_2
      real, dimension(2,mesh%nh_f%loc) :: ele_Xf, ele_Xf_2
      real, dimension(mesh%nu%loc,mesh%nu%ngi,2) :: dnu_t, dnu_t_2
      real, dimension(mesh%nh%loc,mesh%nh%ngi,2) :: dnh_t
      integer, dimension(:), pointer :: b_seg, b_seg_2
      real, dimension(mesh%nh%ngi) :: detwei, detwei_2 
      real, dimension(mesh%nh_f%ngi) :: detwei_f
      real, dimension(2,mesh%nh_f%ngi) :: normal

      ! List of neighbours of current element.
      integer, dimension(:), allocatable :: neigh

      ! Local element information
      integer, dimension(:), pointer :: u_ele, u_ele_2, h_ele, h_ele_2

      real :: val1, llength

      real, dimension(3,3) :: bcmatnc
      real, dimension(2,2) :: bcmatcg
      real, dimension(3,3) :: cg2ncmat, nc2cgmat

      ! Local node number map for big NC element.
      integer, dimension(9) :: local_glno

      ! Local integration matrices for big NC element.
      real, dimension(2,3,9) :: B
      real, dimension(2,2,9,9) :: BQB
      real, dimension(3,3) :: Q 
      real, dimension(2,2,2) :: Q_surf_cg
      real, dimension(2,3,3) :: Q_surf_nc

      logical :: swap0

      MSG("subroutine assemble_GN_operator")

      cg2ncmat(:,1) = (/-1.0, 1.0,  1.0/)
      cg2ncmat(:,2) = (/1.0, -1.0,  1.0/)
      cg2ncmat(:,3) = (/1.0,  1.0, -1.0/)

      nc2cgmat(:,1) = (/0.0, 0.5, 0.5/)
      nc2cgmat(:,2) = (/0.5, 0.0, 0.5/)
      nc2cgmat(:,3) = (/0.5, 0.5, 0.0/)

      MSG("Allocating memory for neigh");
      allocate(neigh(row_length(b_seg_list,1)))

      call zero(mom11)
      call zero(mom12)
      call zero(mom21)
      call zero(mom22)

      rhs1 = 0.
      rhs2 = 0.

      ele_loop: do ele = 1, n_elements

         u_ele=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
         h_ele=>mesh%EVList_h((ELE-1)*mesh%Nh%LOC+1:ELE*mesh%Nh%LOC)

         ele_X(1,:)=mesh%X(h_ele)
         ele_X(2,:)=mesh%Y(h_ele)
         call transform_to_physical(ele_X, mesh%nh, m = mesh%nu, &
              dm_t = dnu_t, &
              detwei = detwei)

         !construct height 

         !Dloc = 0.
         !do iloc = 1, mesh%nh%loc
         !   globi = h_ele(iloc)
         !   Dloc = Dloc + mesh%nh%n(iloc,:)*D(globi)
         !end do

         !volume integrals
         iloc_loop: do iloc = 1, mesh%nu%loc
            globi = u_ele(iloc)
            jloc_loop: do jloc = 1, mesh%nu%loc
               globj = u_ele(jloc)

               !mass part
               kmat = sum(mesh%nu%n(iloc,:)*mesh%nu%n(jloc,:) * &
                    detwei(:))

               rhs1(globi) = rhs1(globi) + kmat*m1(globj)
               rhs2(globi) = rhs2(globi) + kmat*m2(globj)

               call addto(mom11,globi,globj,kmat)
               call addto(mom22,globi,globj,kmat)

            end do jloc_loop
         end do iloc_loop

         if(.true.) then

            !get local mass matrix
            Q=shape_shape(mesh%nu,mesh%nu,detwei)

            do iloc = 1,3
               Q(iloc,iloc) = 1.0/Q(iloc,iloc)
            end do

            ! First part of local numbering is for this element.
            local_glno(1:3)=u_ele
            local_glno(4:)=0

            !------------------------------------------------------------------
            ! Element internal integral.
            !------------------------------------------------------------------

            B=0.0
            B(:,:,1:3)= shape_dshape(mesh%nu, dnu_t, detwei)
            !B(:,:,1:3)= dshape_shape(dnu_t, mesh%nu, detwei)
            
            if (size(neigh)/=row_length(mesh%b_seg_list,ele)) then
               deallocate(neigh)
               allocate(neigh(row_length(mesh%b_seg_list,ele)))
            end if

            neigh=row_m(mesh%b_seg_list,ele)

            !loop over neighbours of ele

            neighbourloop: do ni=1,size(neigh)
               ele_2=neigh(ni)

               ! Skip external b_segs (neumann bcs)
               if (ele_2==0) cycle neighbourloop

               u_ele_2=> &
                    mesh%EVList_u((ELE_2-1)*mesh%nu%LOC+1:ELE_2*mesh%nu%LOC)
               h_ele_2=> &
                    mesh%EVList_h((ELE_2-1)*mesh%Nh%LOC+1:ELE_2*mesh%Nh%LOC)

               b_seg_i=ival(mesh%b_seg_list, ele, ele_2)
               b_seg=>mesh%b_seg_nh_lno((b_seg_i-1)*2+1:b_seg_i*2)

               b_seg_i=ival(mesh%b_seg_list, ele_2, ele)
               b_seg_2=> mesh%b_seg_nh_lno((b_seg_i-1)*2+1:b_seg_i*2)

               ! Locations of local vertices.
               ele_X_2(1,:)=mesh%X(h_ele_2)
               ele_X_2(2,:)=mesh%Y(h_ele_2)

               ! Locations of b_seg vertices.
               ele_Xf=ele_X(:,b_seg)

               ! Change of coordinates in second element.
               call transform_to_physical(ele_X_2, mesh%nh, &
                    m = mesh%nu, dm_t = dnu_t_2, &
                    detwei = detwei_2)
               call transform_face_to_physical(ele_X, ele_Xf, &
                    mesh%nh, mesh%nh_f, &
                    detw_f = detwei_f,normal = normal)

               ! Values of local node map
               local_glno(3+(ni-1)*2+1:3+ni*2)= u_ele_2(b_seg_2)

               Q_surf_cg=shape_shape_vector(mesh%nh_f,mesh%nh_f, &
                    detwei_f, normal)
               !Q_surf_cg=0.0

               ! Interior face part of average.
               do i = 1,2
                  Q_surf_nc(i,:,:)=matmul(cg2ncmat(:,b_seg), &
                       matmul(Q_surf_cg(i,:,:), &
                       cg2ncmat(b_seg,:)))
               end do
               B(:,b_seg,b_seg)=B(:,b_seg,b_seg)-0.5*Q_surf_nc(:,b_seg,b_seg)

               ! Exterior face part of average
               do i = 1,2
                  Q_surf_nc(i,:,:)=matmul(cg2ncmat(:,b_seg), &
                       matmul(Q_surf_cg(i,:,:), &
                       cg2ncmat(b_seg_2,:)))
               end do

               B(:,b_seg,3+(ni-1)*2+1:3+ni*2)= &
                    +0.5*Q_surf_nc(:,b_seg,b_seg_2)
            end do neighbourloop

            BQB=0.0

            !BQB=BQB+matmul(matmul(transpose(B(i,:,:)),Q),B(i,:,:))
            forall(i = 1:2, j = 1:2)
               BQB(i,j,:,:)=BQB(i,j,:,:) &
                    +matmul(matmul(transpose(B(i,:,:)),Q),B(j,:,:))
                    !+matmul(transpose(B(i,:,:)),matmul(Q,B(j,:,:)))
            end forall

            ! Put the contribution for this element into the matrix.
            iloop: do iloc=1,9
               globi=local_glno(iloc)

               ! Exclude boundaries.
               if (globi==0) cycle iloop

               jloop: do jloc=1,9
                  globj=local_glno(jloc)

                  ! Exclude boundaries.
                  if (globj==0) cycle jloop

                  ! Insert value in the matrix.
                  call addto(mom11, globi, globj,alpha1*alpha1* &
                       BQB(1,1,iloc,jloc))
                  call addto(mom12, globi, globj,alpha1*alpha1* &
                       BQB(1,2,iloc,jloc))
                  call addto(mom21, globi, globj,alpha1*alpha1* &
                       BQB(2,1,iloc,jloc))
                  call addto(mom22, globi, globj,alpha1*alpha1* &
                       BQB(2,2,iloc,jloc))
               end do jloop

            end do iloop

         end if

      end do ele_loop

      MSG("end subroutine assemble_GN_operator")

    end subroutine assemble_GN_operator

    subroutine cg_solve_block_eqn(x, A, b, rots, preerr, prenoi)
      ! routine to encapsulate the call to the solver.
      ! Solves Ax=b with initial condition for x passed in
      real, dimension(:), intent(inout) :: X
      type(block_csr_matrix), intent(in) :: A
      real, dimension(:), intent(in) :: b
      type(rotations) :: rots
      real, intent(in) :: preerr
      integer, intent(in) :: prenoi

      ! Variables used by the conjugate gradient solver.
      integer :: imatst, timm, tmisou, kits, para, tnoit1
      real :: trelax, teror1
      logical :: momsym
      real, dimension(:,:), allocatable :: workv
      real, dimension(:), allocatable :: work1d
      real :: rdum(0)
      integer :: iidum(0)
      integer :: nonods

      imatst=0
      timm=1
      tmisou=1
      trelax=1.0
      momsym=.false.
      teror1=preerr
      tnoit1=prenoi
      iidum=0
      allocate(workv(size(X), 5))
      allocate(work1d(4*size(X)))

      nonods = size(X)/2


!      rots%rotate = .false.
      CHECK(nonods)
      CHECK(size(X))
      CHECK(size(A%colm))
      CHECK(rots%rotate)
      CHECK(rots%n_normals)
      CHECK(shape(rots%bcnormals(:,:)))
      CHECK(shape(rots%DM1))
      CHECK(shape(rots%DM2))
      CHECK(shape(rots%DM3))
      CHECK(shape(rots%rtdr))

      MSG("going into SOLCG")
      call pminmx(b,size(b),'******difvec  ')

      if(rots%rotate) then
      call solcg(X, b, nonods, size(X), nonods, .true.,  &
           A%val, A%findrm, A%colm, &
           size(A%colm),4*size(A%colm), 
           kits)
   else
      call solcg(X, b, nonods, size(X), nonods, .true.,  &
           A%val, A%findrm, A%colm, &
           size(A%colm),4*size(A%colm), kits)
   end if
   CHECK(rots%rotate)
      MSG("just out of SOLCG")

    end subroutine cg_solve_block_eqn

end program main
