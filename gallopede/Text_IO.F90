#include "fdebug.h"

module text_io

  use global_parameters_gallopede
  use sparse_tools, only : csr_matrix,val
  use data_structures
  use fldebug
  use State_module
  use fields
  use FETools
  use fields_data_types
  use read_triangle
  use mesh_tools, only : get_tangents

  implicit none

contains 

  subroutine Get_Parameters()

    implicit none
    integer :: io1
    CHARACTER(14), PARAMETER :: data_file = 'Parameters.dat'
    logical :: file_exists
    namelist/params/ tmax, dt, dtmax, tdumpmax, &
        theta, mom_maxnits, den_maxnits, g0, kappa, alpha, cpu_lim, &
        tidal_t,tidal_mag

    !--------------------------------------------------------
    !read in info

    ewrite(2,*) "subroutine Get_Parameters"

    ewrite(2,*) "Opening "//data_file 

    inquire(file=data_file,exist=file_exists)

    if(.not.file_exists) then
       write(0,*) data_file//' does not exist'
       stop
    end if

    open(unit=2502, file=data_file, status='old', &
         iostat=io1)

    if(io1.ne.0) then

       write (0,*) 'Could not read from '// data_file //'file'
       stop
    else
       read(2502,iostat=io1,nml=params)

       print*, io1

    ewrite(2,*) 'tmax=', tmax
    ewrite(2,*) 'tdumpmax=', tdumpmax
    ewrite(2,*) 'theta=', theta
    ewrite(2,*) 'mom_maxnits=', mom_maxnits
    ewrite(2,*) 'den_maxnits=', den_maxnits
    ewrite(2,*) 'g0=', g0
    ewrite(2,*) 'kappa=', kappa
    ewrite(2,*) 'alpha=', alpha
    ewrite(2,*) 'cpu_lim=',cpu_lim
    ewrite(2,*) 'tidal_t=', tidal_t
    ewrite(2,*) 'tidal_mag=', tidal_mag

       assert( tmax .ge. 0.0 )
       assert( dt >0.0 )
       assert( dtmax > 0.0)
       assert( tdumpmax > 0.0)
       assert( theta .ge. 0.5)
       assert( mom_maxnits .ge. 1)
       assert( den_maxnits .ge. 1)
       assert( g0 .ge. 0.0)
       assert( kappa .ge. 0.0)
       assert( alpha .ge. 0.0)

       close(unit=2502)
    end if

    ewrite(3,*) (tmax)
    ewrite(3,*)(tdumpmax)
    ewrite(3,*)(theta)
    ewrite(3,*)(mom_maxnits)
    ewrite(3,*)(den_maxnits)
    ewrite(3,*)(g0)
    ewrite(3,*)(eta)

    ewrite(2,*)("END subroutine Get_Parameters")

  end subroutine Get_Parameters

  subroutine Get_Mesh(mesh,bcs,m,u,d,u_bar,h,bottom)
    
    implicit none
    type(dg_mesh) :: mesh
    type(bc_info), optional :: bcs
    type(vector_field), dimension(:), intent(inout) :: m, u
    type(scalar_field), dimension(:), intent(inout) :: d
    type(vector_field), intent(inout) :: u_bar
    type(scalar_field), intent(inout) :: h,bottom

    !locals
    type(mesh_type), pointer :: mesh_v
    type(mesh_type) :: d_mesh,m_mesh, u_mesh
    type(state_type) :: positions_state,variables_state,connectivity_state
    type(scalar_field), pointer :: field_in
    type(vector_field), target :: positions,connectivity
    integer, dimension(:), allocatable, target :: tmp_store
    integer, dimension(:), pointer :: ele_2
    integer, dimension(:), allocatable :: bc_tmp
    integer :: layer,tangent_i,interior_i,i,j,face_u
    character(len=6) :: num_str
    !--------------------------------------------------------
    !allocate memory

    call nullify(Positions_state)
    call nullify(Variables_state)
    
    Variables_state = read_triangle_files('Variables',mesh%nh,.true.,3)
    Connectivity_state = read_triangle_files('Connectivity',mesh%nu,.true.)
    Positions_state = read_triangle_files('Positions',mesh%nu,.true.)

    !--------------------------------------------------------
    !read in connectivity for X

    allocate(mesh%positions)
    mesh%positions=extract_vector_field(positions_state,'Coordinate')
    
    call add_faces(mesh%positions%mesh,mesh%positions%mesh)
    N_verts = node_count(mesh%positions)
    mesh%N_verts = N_verts
    Mesh%N_element = N_elements

    !--------------------------------------------------------
    !read in connectivity for fields

    allocate(mesh%connectivity)
    mesh%connectivity=extract_vector_field(connectivity_state,'Coordinate')
    call add_faces(mesh%connectivity%mesh,mesh%connectivity%mesh)

    print*, 'N_verts=', node_count(mesh%positions)
    print*, 'N_elements=', ele_count(mesh%positions)
    print*, 'N_free=', node_count(mesh%connectivity)

    !connectivity for velocity and height
    ewrite(1,*) 'Getting u connectivity'
    m_mesh = make_mesh(mesh%connectivity%mesh,mesh%nh,0,'u_mesh')
    call add_faces(m_mesh,mesh%connectivity%mesh)
    ewrite(1,*) 'Getting h connectivity'
    D_mesh = make_mesh(mesh%connectivity%mesh,mesh%nh,0,'h_mesh')
    call add_faces(d_mesh,mesh%connectivity%mesh)
    u_mesh = make_mesh(mesh%connectivity%mesh,mesh%nu,-1,'v_mesh')
    call add_faces(u_mesh,mesh%connectivity%mesh)

     do layer=1,n_layers
     call allocate(m(layer),2,m_mesh,'smoothed pressure')
     call allocate(u(layer),2,u_mesh,'velocity')
     call allocate(D(layer),d_mesh,'thickness')
  end do
   call allocate(bottom,d_mesh,'bottom')
   call allocate(H,d_mesh,'Barotropic Height')
   call allocate(u_bar,2,u_mesh,'Barotropic Velocity')
    !---------------------------------------------------------
    !read in initial conditions
    
    do layer = 1, N_layers
       !read u1
       write(num_str,'(i0)') (layer-1)*3 + 1
       field_in => extract_scalar_field(variables_state, &
            & 'Attributes '//trim(num_str))
       u(layer)%val(1)%ptr = map_quad_to_DG(field_in,u(1))
       !read u2
       write(num_str,'(i0)') (layer-1)*3 + 2
       field_in => extract_scalar_field(variables_state, &
            & 'Attributes '//trim(num_str))
       u(layer)%val(2)%ptr = map_quad_to_DG(field_in,u(1))
       !read D
       write(num_str,'(i0)') layer*3
       field_in => extract_scalar_field(variables_state, &
            & 'Attributes '//trim(num_str))
       call map_data_to_CG(D(layer),field_in)
    end do
    write(num_str,'(i0)') N_layers*3 + 1
    field_in => extract_scalar_field(variables_state, &
         & 'Attributes '//trim(num_str))
    call map_data_to_CG(bottom,field_in)

    do layer=1,n_layers
       m(layer)%val(1)%ptr=0.0
       m(layer)%val(2)%ptr=0.0
    end do

!    print*, 'average h_0=', &
!         sum(abs(bottom+D( (/ ( (i-1)*n_dens+(1:n_dens),i=1,n_layers) /)))))/n_dens

    !need bcs%N_interior = vec(5) bcs%N_tangents = vec(6)

    if(present(bcs)) then
       field_in => extract_scalar_field(variables_state, &
         & 'Boundary Marker')

       bcs%N_tangents = 0
       do i = 1, ele_count(u(1))
          ele_2=>ele_neigh(u(1),i)
          if(any(ele_2<0)) then
             bcs%N_tangents = bcs%N_tangents + 1
             do j=1,size(ele_2)
                if (ele_2(j) <0)&
                     bcs%N_tangents = bcs%N_tangents + 1
             end do
          end if
       end do
       bcs%N_interior = mesh%N_vels-bcs%N_tangents 
       ewrite(1,*) 'n tangents = ', bcs%N_tangents
       ewrite(1,*) 'n vels = ', mesh%N_vels
       allocate(bcs%interior_list(bcs%N_interior) ) 
       allocate(bcs%tangent_list(bcs%N_tangents) )
       allocate(bc_tmp(node_count(u(1) )))
       tangent_i = 0 
       interior_i = 0
       bc_tmp=0

       print *,'asdfas'



        do i = 1, ele_count(u(1))
           ele_2=>ele_neigh(u(1),i)
           do j=1,size(ele_2)
              face_u=ele_face(u(1),i,ele_2(j))
              if (ele_2(j)<0) then
                 bc_tmp(face_global_nodes(u(1),face_u))=1 
              end if
           end do
        end do

        print*, sum(bc_tmp)

        do i=1,node_count(u(1))
           if (bc_tmp(i)==0) then
              interior_i = interior_i + 1
              bcs%interior_list(interior_i) =i
           else
              tangent_i = tangent_i + 1
              bcs%tangent_list(tangent_i) = i
           end if
       end do
       if(tangent_i.ne.bcs%N_tangents) then
          FLAbort('N_tangents wrong')
       end if
       if(tangent_i.ne.bcs%N_tangents) then
          FLAbort('N_interiors wrong')
       end if
       if (bcs%N_tangents>0) then
          allocate( bcs%tangents(2,bcs%N_tangents) )
          call get_tangents(bcs,mesh,bc_tmp,u(1),mesh%positions)
       end if

       !test normals, comment out for operational 

       !if(bcs%N_tangents>0) then
       !   call test_normals(mesh,bcs)
       !end if

    end if

    call nullify(Positions_state)
    call nullify(Variables_state)

  end subroutine Get_Mesh

  subroutine read_field(Dr,Dint,filename)

    implicit none
    real, intent(out), dimension(:), optional :: Dr
    integer, intent(out), dimension(:), optional :: Dint
    character(len=*), intent(in) :: filename

    integer :: io1, nod
    logical :: file_exists
    integer :: sizeD

    !--------------------------------------------------------
    !read in info

    ewrite(2,*)('Subroutine read_field')

    ewrite(2,*)('reading file')
    ewrite(3,*)(filename)

    open(unit=2502, file=filename, status='old', &
         iostat=io1)

    sizeD = 0
    if(present(Dr).and.present(Dint)) then
       FLAbort('cant read in real and int in same call!')
    end if
    if(present(Dr)) sizeD = size(Dr)
    if(present(Dint)) sizeD = size(Dint)
    if(sizeD==0) FLAbort('didnt want to read anything?')

    if(io1.ne.0) then
       write (0,*) 'Could not open file', filename, ' for reading'
       stop
    else
       do nod = 1, sizeD
          if(present(Dr)) then
             read(unit=2502, iostat=io1, fmt='(f20.16)') Dr(nod)
          else
             read(unit=2502, iostat=io1, fmt=*) Dint(nod)
          end if
          if(io1 /= 0) then
             write(0,*) 'could not read from file', filename
             stop
          end if
       end do
       close(unit=2502)
    end if

    ewrite(2,*)('end Subroutine read_field')

  end subroutine read_field

  subroutine map_data_to_CG(quad_field,data_field)
     type(scalar_field) :: Quad_field, data_field

     integer :: ele
     
     do ele=1,ele_count(quad_field)
        quad_field%val(ele_nodes(quad_field,ele))=ele_val(data_field,ele)
     end do

   end subroutine map_data_to_CG

  function map_quad_to_DG(Quad_field,DG_field) result(DG)
    type(scalar_field) :: Quad_field
    type(vector_field) :: DG_field
    real, allocatable, dimension(:) ::DG
    real ::  CG_locgi(ele_ngi(Quad_field,1))
    real ::  CG_loc(ele_loc(Quad_field,1))
!    real:: detwei(mesh%positions%mesh%shape%ngi)
!    real:: X_ele(2,mesh%positions%mesh%shape%loc)

    integer :: ele
    type(element_type), pointer :: shape_u,shape_X

    allocate(DG(node_count(DG_field)))

    do ele=1,ele_count(DG_field)

!       X_ele=ele_val(mesh%positions,ele)
!       shape_X=>ele_shape(mesh%positions, ele)
!       shape_u=>ele_shape(DG_field, ele)
       CG_loc=ele_val(Quad_field,ele)
       CG_locgi=ele_val_at_quad(Quad_field,ele)
       DG(ele_nodes(DG_field,ele))=cg_loc((/1,3,6/))
    end do


  end function map_quad_to_DG

  subroutine dump_field(D,size_field,filename,filecount,pad,flength)


    implicit none
    integer, intent(in) :: size_field,flength
    real, intent(in), dimension(size_field) :: D
    integer, intent(in) :: filecount, pad

    integer :: io1, nod
    character(len=*), intent(in) :: filename
    CHARACTER(LEN=flength) :: data_file
    CHARACTER(LEN=pad) :: zeros
    CHARACTER(LEN=10) :: fmt
    logical :: file_exists
    integer :: filecount_digs, i

    !--------------------------------------------------------
    !read in info

    ewrite(2,*)('Subroutine dump_field')

    filecount_digs = floor(log10(1.0*filecount)) + 1

    if(filecount_digs>pad) then
       FLAbort('Need to increase padding in filename')
       stop
    end if

    ewrite(3,*)(pad)

    zeros = '0000000000000000000000000000000000000000000000000000000'
    ewrite(3,*)(zeros)

    write(fmt,'(i6)')  filecount_digs
    fmt = '(i'//trim(fmt)//')'
    write(data_file,fmt)  filecount
    data_file = trim(filename)//trim(zeros)//trim(data_file)//'.dat'

    ewrite(2,*)('Writing file')
    ewrite(3,*)(data_file)

    open(unit=2502, file=data_file, status='new', &
         iostat=io1)

    if(io1.ne.0) then
       write (0,*) 'Could not open file', data_file, ' for writing'
       stop
    else
       do nod = 1, size_field
          write(unit=2502, iostat=io1, fmt=*) D(nod)
          if(io1 /= 0) then
             write(0,*) 'could not read from file', data_file
             stop
          end if
       end do
       close(unit=2502)
    end if

    ewrite(2,*)('end Subroutine dump_field')

  end subroutine dump_field

  subroutine dump_matrix(Mat,nonods)

    implicit none
    type(csr_matrix), intent(in) :: Mat
    integer, intent(in) :: nonods

    integer :: io1
    logical :: file_exists
    integer :: i, globi, globj

    !--------------------------------------------------------
    !read in info

    ewrite(2,*)('Subroutine dump_matrix')


    ewrite(2,*)('Writing file')
    open(unit=2502, file='matdump.dat', status='new', &
         iostat=io1)

    if(io1.ne.0) then
       write (0,*) 'Could not open matrix file for writing'
       stop
    else
       do globi = 1, nonods
          do globj = 1, nonods
             write(unit=2502, iostat=io1, fmt=*) val(Mat,globi,globj)
             if(io1 /= 0) then
                write(0,*) 'could not write to matrix file'
                stop
             end if
          end do
       end do
       close(unit=2502)
    end if

    ewrite(2,*)('end Subroutine dump_matrix')

  end subroutine dump_matrix

  subroutine construct_dg_connectivity(EVList,Nelements)
    integer, dimension(:), pointer :: EVList
    integer, intent(in) :: Nelements

    !locals
    integer :: i

    allocate( EVList(3*NElements) )

    do i = 1, NElements*3
       EVList(i) = i
    end do

  end subroutine construct_dg_connectivity

  subroutine construct_linear_connectivity(mesh)
    type(dg_mesh), intent(inout) :: mesh
    integer :: i


    allocate(mesh%EVlist_cX(N_elements*3))

    if (mesh%nh%loc ==3) then
       mesh%EVlist_cX=mesh%EVList_h
    else
       mesh%EVlist_cX=mesh%EVList_h((/((i-1)*mesh%nh%loc+(/ 1,3,6 /), i=1,n_elements)/))
    end if


    n_pres=0
    do i=1,3*n_elements
      n_pres=max(mesh%EVList_cX(i),n_pres)
   end do

  end subroutine construct_linear_connectivity

  subroutine test_normals(mesh,bcs)
    type(dg_mesh), intent(in) :: mesh
    type(bc_info), intent(in) :: bcs
    
    !locals
    integer :: i, ele, iloc, io1
    real, allocatable, dimension(:,:) :: Xtans
    ! Locations of nodes.
    integer, dimension(:), pointer :: ele_u,ele_X
    real, dimension(2,mesh%nu%loc) :: X_ele

    allocate( Xtans(2,bcs%N_tangents) )
    Xtans = 0.

    print *,'asdf'
    
    do ele = 1, mesh%N_element
       ele_u=>mesh%EVList_u((ELE-1)*mesh%Nu%LOC+1:ELE*mesh%Nu%LOC)
       ele_x=>mesh%EVList_X((ELE-1)*3+1:ELE*3)       
       if(mesh%nu%loc==3) then
          X_ele(1,:) = mesh%X(ele_X)
          X_ele(2,:) = mesh%Y(ele_X)
       else
          X_ele(1,(/1,3,6/)) = mesh%X(ele_x)
          X_ele(2,(/1,3,6/)) = mesh%Y(ele_x)
          X_ele(:,2) = 0.5*(X_ele(:,1) + X_ele(:,3))
          X_ele(:,4) = 0.5*(X_ele(:,1) + X_ele(:,6))
          X_ele(:,5) = 0.5*(X_ele(:,3) + X_ele(:,6))
       end if

       do iloc = 1, mesh%nu%loc
          do i = 1, size(bcs%tangent_list)
             if(ele_u(iloc)==bcs%tangent_list(i)) then
                Xtans(:,i) = X_ele(:,iloc)
             end if
          end do
       end do
       
    end do

    open(unit=2502, file='tangentstest.dat', status='new', &
         iostat=io1)
    if(io1.ne.0) then
       write(0,*) 'Failed to open file for dumping tangents'
    end if
    do i = 1, bcs%N_tangents
       write(unit=2502, iostat=io1, fmt=*) Xtans(1,i), &
            Xtans(2,i), bcs%tangents(1,i), bcs%tangents(2,i)
       if(io1 /= 0) then
          write(0,*) 'could not write to matrix file'
          stop
       end if
    end do

    close(unit=2502)    

    ewrite(1,*) 'Stopping after dumping normals.'
    stop
    
  end subroutine test_normals

  subroutine vec_copy(vec1,vec2)

    type(vector_field) :: vec1,vec2
    integer :: i
    type(mesh_type) :: mesh

    mesh = make_mesh(vec2%mesh,vec2%mesh%shape,0,'copied mesh')

    call allocate(vec1,vec2%dim,mesh,'copied_vector')

    do i=1,vec2%dim
       vec1%val(i)%ptr=vec2%val(i)%ptr
    end do
  end subroutine vec_copy

end module text_io
