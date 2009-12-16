#include "fdebug.h"

module vtk_io

  use global_parameters_gallopede
  use Multilayer_tools
  use fldebug 
  use data_structures

  contains

  subroutine dump_field_vtk(D,nd,EVList_D,EVList_X,X,Y,nX,nloc, &
       filename,filecount,pad,flength,totele)

    !dumps dg data

    implicit none
    integer, intent(in) :: flength,totele,nx,nd,nloc
    integer, dimension(:),target :: EVList_X,EVList_D
    real, intent(in), dimension(nd) :: D
    integer, intent(in) :: filecount, pad
    real, intent(in), dimension(nx) :: X,Y

    integer :: io1, nod, ele
    character(len=*), intent(in) :: filename
    CHARACTER(LEN=flength) :: data_file
    CHARACTER(LEN=pad) :: zeros
    CHARACTER(LEN=10) :: fmt
    logical :: file_exists
    integer :: filecount_digs, i
    integer, allocatable::eltypes(:),elsizes(:)
    real, allocatable, dimension(:) :: XD,YD
    integer, dimension(:), pointer :: d_ele, x_ele

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

    write(fmt,'(i6)') filecount_digs
    fmt = '(i'//trim(fmt)//')'
    write(data_file,fmt) filecount
    data_file = trim(filename)//trim(zeros)//trim(data_file)//'.vtu'

    ewrite(2,*)('Writing file')
    ewrite(3,*)(data_file)

    allocate( eltypes(totele) )
    allocate( elsizes(totele) )

    call VTKOPEN(data_file,flength,' ',1)

    do i = 1, totele
       eltypes(i) = 5
       elsizes(i) = 3 
    end do

    allocate( XD( nd ) )
    allocate( YD( nd ) )

    do ele = 1, totele
       d_ele=>EVList_d((ELE-1)*3+1:ELE*3)
       X_ele=>EVList_X((ELE-1)*3+1:ELE*3)
       XD(D_ELE) = x(x_ELE)
       yD(D_ELE) = Y(x_ELE)
    END do

    call VTKWRITEMESHD(nd,totele,xd,yd,D,EVlist_d,eltypes,elsizes)
    call VTKWRITEDSN(D,'D',1)

    call VTKCLOSE()

    deallocate( eltypes )
    deallocate( elsizes )
    deallocate( XD )
    deallocate( YD )
    
    ewrite(2,*)('end Subroutine dump_field_vtk')

  end subroutine dump_field_vtk

  subroutine dump_data_vtk(u1,u2,m1,m2,D,mesh,filename,filecount,pad,flength)

    !dumps dg data

    implicit none
    integer, intent(in) :: flength
    real, intent(in), dimension(N_vels) :: u1,u2
    real, intent(in), dimension(N_moms) :: m1,m2
    real, intent(in), dimension(N_dens) :: D
    type(dg_mesh), intent(in) :: mesh
    integer, intent(in) :: filecount, pad

    integer :: io1, nod, ele
    character(len=*), intent(in) :: filename
    CHARACTER(LEN=flength) :: data_file
    CHARACTER(LEN=pad) :: zeros
    CHARACTER(LEN=10) :: fmt
    logical :: file_exists
    integer :: filecount_digs, i, nx
    integer, allocatable::eltypes(:),elsizes(:)
    real, allocatable, dimension(:) :: Dx, u1x, u2x, m1x, m2x
    integer, dimension(:), pointer :: d_ele, x_ele, u_ele, m_ele

    !--------------------------------------------------------
    !read in info

    ewrite(2,*)('Subroutine dump_field')

    ewrite(3,*) 'max(u1) before', maxval(u1)
    ewrite(3,*) 'max(u2) before', maxval(u2)
    
    nx = size(mesh%X)
 
    filecount_digs = floor(log10(1.0*filecount)) + 1

    if(filecount_digs>pad) then
       FLAbort('Need to increase padding in filename')
       stop
    end if

    ewrite(3,*)(pad)

    zeros = '0000000000000000000000000000000000000000000000000000000'
    ewrite(3,*)(zeros)

    write(fmt,'(i6)') filecount_digs
    fmt = '(i'//trim(fmt)//')'
    write(data_file,fmt) filecount
    data_file = trim(filename)//trim(zeros)//trim(data_file)//'.vtu'

    ewrite(2,*)('Writing file')
    ewrite(3,*)(data_file)

    allocate( eltypes(N_elements) )
    allocate( elsizes(N_elements) )

    call VTKOPEN(data_file,flength,' ',1)

    do i = 1, N_elements
       ewrite(3,*) i, N_elements
       eltypes(i) = 5
       elsizes(i) = 3
    end do

    !convert u,D to linear elements

    allocate( Dx( nx ) )
    allocate( u1x( nx ), u2x( nx ) )
    allocate( m1x( nx ), m2x( nx ) )

    do ele = 1, N_elements
       d_ele=>mesh%EVList_h((ELE-1)*mesh%nh%loc+1:ELE*mesh%nh%loc)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       u_ele=>mesh%EVList_u((ELE-1)*mesh%nu%loc+1:ELE*mesh%nu%loc)
       m_ele=>mesh%EVList_m((ELE-1)*mesh%nm%loc+1:ELE*mesh%nm%loc)
       select case (mesh%nh%loc)
       case (3)
          Dx(x_ele) = D(d_ele)
       case (6)
          Dx(x_ele) = D(d_ele( (/1,3,6/) ))
       case default
          FLAbort('wrong')
       end select
       select case (mesh%nu%loc)
       case (3)
          u1x(x_ele) = u1(u_ele)          
          u2x(x_ele) = u2(u_ele)
       case (6)
          u1x(x_ele) = u1(u_ele( (/1,3,6/) ))
          u2x(x_ele) = u2(u_ele( (/1,3,6/) ))
       case default
          FLAbort('wrong')
       end select
       m1x(x_ele) = m1(m_ele)       
       m2x(x_ele) = m2(m_ele)
    END do

    ewrite(3,*) 'max(u1)', maxval(u1x)
    ewrite(3,*) 'max(u2)', maxval(u2x)

    call VTKWRITEMESHD(nx,N_elements,mesh%x,mesh%y,dx, &
         mesh%EVlist_X,eltypes,elsizes)
    call VTKWRITEDSN(Dx,'Layer Depth',11)
    call VTKWRITEDVN(U1x,u2x,0.0*U1x,'Velocity',8)
    call VTKWRITEDVN(m1x,m2x,0.0*U1x,'Momentum',8)

    call VTKCLOSE()

    deallocate( eltypes )
    deallocate( elsizes )
    deallocate( u1x )
    deallocate( u2x )
    deallocate( m1x )
    deallocate( m2x )
    deallocate( Dx )

    ewrite(2,*)('end Subroutine dump_field_vtk')

  end subroutine dump_data_vtk

  subroutine dump_data_vtk_layer_quad(u,m,D,bottom, &
       mesh,filename,filecount,pad,flength,n1,n2,p,nn1,nn2,np,u_cont)

    !dumps dg data

    implicit none
    integer, intent(in) :: flength,filecount, pad
    real, intent(in), dimension(:) :: D
    type(dg_mesh), intent(in) :: Mesh
    real, intent(in), dimension(:) :: bottom
    real, intent(in), dimension(:) :: u
    real, intent(in), dimension(:) :: m
    real, intent(in), dimension(:), optional :: n1,n2,p,u_cont
    real, intent(in), dimension(:), optional :: nn1,nn2,np

    integer :: io1, nod, ele, layer_j
    character(len=*), intent(in) :: filename
    CHARACTER(LEN=flength) :: data_file
    CHARACTER(LEN=pad) :: zeros
    CHARACTER(LEN=10) :: fmt
    logical :: file_exists
    integer :: filecount_digs, i, layer_i, n_dg
    integer, allocatable::eltypes(:),elsizes(:)
    integer, allocatable, target :: evlist_out(:)
    real, allocatable, dimension(:) :: Xdg,Ydg,Zdg,Ddg,U1dg,U2dg,m1dg,m2dg
    real, allocatable, dimension(:) :: bottomdg, ones, n1dg,n2dg,pdg
    real, allocatable, dimension(:) :: nn1dg,nn2dg,npdg
    real, allocatable, dimension(:) :: u1_contdg, u2_contdg
    integer, dimension(:), pointer :: d_ele, x_ele, u_ele, m_ele,dg_ele,cX_ele
   

    assert(size(D)==n_layers*n_dens)
    assert(size(u)==2*n_layers*n_vels)
    assert(size(m)==2*n_layers*n_moms)
    assert(size(bottom)==n_dens)

    if (present(n1)) then
       assert(size(n1)==n_layers*n_moms)
       assert(size(n2)==n_layers*n_moms)
    end if

    !--------------------------------------------------------
    !read in info



    ewrite(2,*)('subroutine dump_data_vtk_layer_quad')

    ewrite(2,*) 'N_Layers', N_layers

    if(filecount>0) then
       filecount_digs = floor(log10(1.0*filecount)) + 1

       ewrite(2,*) 'filecount_digs, pad', filecount_digs, pad
       
       if((10**(pad+filecount_digs)-filecount)<2) then
          FLAbort('Need to increase padding in filename')
          stop
       end if
       
       ewrite(3,*)(pad)
       
       zeros = '0000000000000000000000000000000000000000000000000000000'
       ewrite(3,*)(zeros)
       
       write(fmt,'(i6)') filecount_digs
       fmt = '(i'//trim(fmt)//')'
       write(data_file,fmt) filecount
       data_file = trim(filename)//trim(zeros)//trim(data_file)//'.vtu'
    else
       data_file = trim(filename)//'.vtu'
    end if

    ewrite(1,*)('Writing file')
    ewrite(1,*)(data_file)

    allocate( eltypes(N_elements*(N_layers+1)) )
    allocate( elsizes(N_elements*(N_layers+1)) )

    call VTKOPEN(data_file,flength,' ',1)

    do i = 1, N_elements*(N_layers+1)
       eltypes(i) = 22
       elsizes(i) = 6
    end do

    !convert X,Y,D to dg elements


    n_dg=6*n_elements
    allocate( ones(n_dg*(N_layers+1)) )
    allocate( Xdg(n_dg*(N_layers+1) ) )
    allocate( Ydg(n_dg*(N_Layers+1) ) )
    allocate( Zdg(n_dg*(N_layers+1) ) )
    allocate( Ddg(n_dg*(N_Layers+1) ) )
    allocate( u1dg(n_dg*(N_layers+1) ) )
    allocate( u2dg(n_dg*(N_Layers+1) ) )
    allocate( m1dg(n_dg*(N_layers+1) ) )
    allocate( m2dg(n_dg*(N_Layers+1) ) )
    allocate( bottomdg(n_dg)  )
    allocate( evlist_out( (N_Layers+1)*N_elements*6) )

    if (present(n1)) then
       allocate( n1dg(n_dg*(N_layers+1) ) )
    end if
    if (present(n2)) then
       allocate( n2dg(n_dg*(N_Layers+1) ) )
    end if
     if (present(p)) then
       allocate( pdg(n_dg*(N_Layers+1) ) )
    end if
    if (present(nn1)) then
       allocate( nn1dg(n_dg*(N_layers+1) ) )
    end if
    if (present(nn2)) then
       allocate( nn2dg(n_dg*(N_Layers+1) ) )
    end if
    if (present(np)) then
       allocate( npdg(n_dg*(N_Layers+1) ) )
    end if
  if (present(u_cont)) then
       allocate( u1_contdg(n_dg*(N_Layers+1) ) )
       allocate( u2_contdg(n_dg*(N_Layers+1) ) )
    end if
    Zdg = 0.

    ones=1;

    ewrite(2,*)("element loop")


    forall (ele=1:6*(n_layers+1)*n_elements)
       EVList_out(ele) = ele
    end forall

    do ele = 1, N_elements
       d_ele=>mesh%EVList_h((ELE-1)*mesh%nh%loc+1:ELE*mesh%nh%loc)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       cX_ele=>mesh%EVList_cX((ELE-1)*3+1:ELE*3)
       u_ele=>mesh%EVList_u((ELE-1)*mesh%nu%loc+1:ELE*mesh%nu%loc)
       m_ele=>mesh%EVList_m((ELE-1)*mesh%nm%loc+1:ELE*mesh%nm%loc)

       do layer_i = 0, N_Layers-1

          dg_ele=>EVList_out((ele-1)*6+1:(ele)*6)

          XDg(dg_ELE + N_dg * layer_i) = to_quad(mesh%x(x_ELE))
          yDg(dg_ELE + N_dg * layer_i) = to_quad(mesh%Y(x_ELE))
          m1dg(dg_ele + N_dg * layer_i) = to_quad(m(layer_i*2*n_moms+m_ele))
          m2dg(dg_ele + N_dg * layer_i) = to_quad(m(layer_i*2*n_moms&
               +n_moms+m_ele))
          u1dg(dg_ele + N_dg * layer_i) = &
               to_quad(u(layer_i*2*n_vels+u_ele))
          u2dg(dg_ele + N_dg * layer_i) = &
               to_quad(u(layer_i*2*n_vels+n_vels+u_ele))
          Ddg(dg_ele + N_dg * layer_i) = &
               to_quad(D(layer_i*n_dens+d_ele))
         if (present(n1)) n1dg(dg_ele + N_dg * layer_i) = &
               to_quad(n1(layer_i*n_moms+m_ele))
         if (present(n2)) n2dg(dg_ele + N_dg * layer_i) = &
              to_quad(n2(layer_i*n_moms+m_ele))
          if (present(p)) pdg(dg_ele + N_dg * layer_i) = &
              to_quad(p(layer_i*n_dens+d_ele))
           if (present(nn1)) nn1dg(dg_ele + N_dg * layer_i) = &
               to_quad(nn1(layer_i*n_moms+m_ele))
         if (present(nn2)) nn2dg(dg_ele + N_dg * layer_i) = &
              to_quad(nn2(layer_i*n_moms+m_ele))
          if (present(np)) npdg(dg_ele + N_dg * layer_i) = &
              to_quad(np(layer_i*n_pres+cX_ele))
           if (present(u_cont)) then
              u1_contdg(dg_ele + N_dg * layer_i) = &
               to_quad(u_cont(layer_i*2*n_vels+u_ele))
          u2_contdg(dg_ele + N_dg * layer_i) = &
               to_quad(u_cont(layer_i*2*n_vels+n_vels+u_ele))
       end if

       end do
       XDg(dg_ELE + N_dg *N_layers) = to_quad(mesh%x(x_ELE))
       yDg(dg_ELE + N_dg *N_layers) = to_quad(mesh%Y(x_ELE))
       u1dg(dg_ele + N_dg *N_layers) = 0.0*dg_ele
       u2dg(dg_ele + N_dg *N_layers) = 0.0*dg_ele
       m1dg(dg_ele + N_dg *N_Layers) = 0.0*dg_ele
       m2dg(dg_ele + N_dg *N_Layers) = 0.0*dg_ele
       Ddg(dg_ele + N_dg *N_layers) = -to_quad(bottom(d_ele))
       bottomdg(dg_ele) = -Ddg(dg_ele + N_dg *N_layers)

       if (present(n1)) then
          n1dg(dg_ele + N_dg * N_layers) =  0.0*dg_ele
       end if
       if (present(n2)) then
          n2dg(dg_ele + N_dg * N_layers) =  0.0*dg_ele
       end if

       if (present(u_cont)) then
          u1_contdg(dg_ele + N_dg *N_layers) = 0.0*dg_ele
          u2_contdg(dg_ele + N_dg *N_layers) = 0.0*dg_ele
       end if

    END do

    do layer_i = 0, N_Layers
       Zdg(N_dg * layer_i + 1: N_dg * (layer_i+1)) = bottomdg
       do layer_j = layer_i, N_layers-1
          Zdg(N_dg * layer_i + 1:N_dg* (layer_i + 1)) = &
               Zdg(N_dg * layer_i + 1:N_dg* (layer_i + 1)) + &
               Ddg(N_dg * layer_j + 1: N_dg * (layer_j + 1))
       end do
    end do

    do layer_i = 0,N_Layers
       Ddg(N_dg*layer_i+1:N_dg*(layer_i+1)) = &
            Ddg(N_dg*layer_i+1:N_dg*(layer_i+1)) - &
            sum(Ddg(N_dg*layer_i+1:N_dg*(layer_i+1)))/N_dg
    end do

!    ewrite(1,*) Xdg
!    FLAbort('sfg')
    ewrite(2,*)("writing the mesh")

    call VTKWRITEMESHD(N_dg*(N_Layers+1),N_elements*(N_Layers+1), &
         xdg,ydg,Zdg,EVlist_out,eltypes,elsizes)
    call VTKWRITEDSN(Ddg,'Layer Depth',11)
    call VTKWRITEDSN(t*ones,'Time',4)
    if (present(p))  call VTKWRITEDSN(pdg,'Pressure',8)
    if (present(np))  call VTKWRITEDSN(npdg,'NLPressure',10)
    print*, maxval(abs(U1dg)), maxval(abs(U2dg))
    call VTKWRITEDVN(U1dg,u2dg,0.0*U1dg,'Velocity',8)    
    print*, maxval(abs(m1dg)), maxval(abs(m2dg))
    call VTKWRITEDVN(m1dg,m2dg,0.0*U1dg,'Momentum',8)
    if (present(n1)) call VTKWRITEDVN(n1dg,n2dg,0.0*U1dg,'Forcing',7)
    if (present(nn1)) call VTKWRITEDVN(nn1dg,nn2dg,0.0*U1dg,'NLForcing',9)
    if (present(u_cont)) call VTKWRITEDVN(u1_contdg,u2_contdg,0.0*U1dg,&
         'Slaved U',8)

    call VTKCLOSE()

    deallocate( eltypes )
    deallocate( elsizes )
    deallocate( XDg )
    deallocate( YDg )
    deallocate( Ddg )
    deallocate( Zdg )
    deallocate( EVList_out )
    if (present(n1)) deallocate(n1dg)
    if (present(n2)) deallocate(n2dg)
    if (present(p)) deallocate(pdg)
    if (present(nn1)) deallocate(nn1dg)
    if (present(nn2)) deallocate(nn2dg)
    if (present(np)) deallocate(npdg)
    if (present(u_cont)) deallocate(u1_contdg,u2_contdg)

    ewrite(2,*)('end Subroutine dump_data_vtk_layer_quad')

  end subroutine dump_data_vtk_layer_quad

  subroutine dump_data_vtk_layer(u,m,D,bottom, &
       mesh,filename,filecount,pad,flength,n1,n2)

    !dumps dg data to linear vtk elements

    implicit none
    integer, intent(in) :: flength,filecount, pad
    real, intent(in), dimension(:) :: D
    type(dg_mesh), intent(in) :: Mesh
    real, intent(in), dimension(:) :: bottom
    real, intent(in), dimension(:) :: u
    real, intent(in), dimension(:) :: m
    real, intent(in), dimension(:), optional :: n1,n2

    integer :: io1, nod, ele, layer_j
    character(len=*), intent(in) :: filename
    CHARACTER(LEN=flength) :: data_file
    CHARACTER(LEN=pad) :: zeros
    CHARACTER(LEN=10) :: fmt
    logical :: file_exists
    integer :: filecount_digs, i, layer_i, n_dg
    integer, allocatable::eltypes(:),elsizes(:)
    integer, allocatable, target :: evlist_out(:)
    real, allocatable, dimension(:) :: Xdg,Ydg,Zdg,Ddg,U1dg,U2dg,m1dg,m2dg
    real, allocatable, dimension(:) :: bottomdg, ones, n1dg,n2dg
    integer, dimension(:), pointer :: d_ele, x_ele, u_ele, m_ele,dg_ele
    integer, dimension(3) :: to_lin

    assert(size(D)==n_layers*n_dens)
    assert(size(u)==2*n_layers*n_vels)
    assert(size(m)==2*n_layers*n_moms)
    assert(size(bottom)==n_dens)

    if (present(n1)) then
       assert(size(n1)==n_layers*n_dens)
       assert(size(n2)==n_layers*n_dens)
    end if

    !--------------------------------------------------------
    !read in info

    ewrite(2,*)('subroutine dump_data_vtk_layer')

    select case(mesh%nu%loc)
       case (3)
          to_lin = (/ 1, 2, 3 /)
       case (6)
          to_lin = (/ 1, 3, 6 /)
       case default
          FLAbort('crazy element')
    end select

    ewrite(2,*) 'N_Layers', N_layers

    if(filecount>0) then
       filecount_digs = floor(log10(1.0*filecount)) + 1

       ewrite(2,*) 'filecount_digs, pad', filecount_digs, pad
       
       if((10**(pad+filecount_digs)-filecount)<2) then
          FLAbort('Need to increase padding in filename')
          stop
       end if
       
       ewrite(3,*)(pad)
       
       zeros = '0000000000000000000000000000000000000000000000000000000'
       ewrite(3,*)(zeros)
       
       write(fmt,'(i6)') filecount_digs
       fmt = '(i'//trim(fmt)//')'
       write(data_file,fmt) filecount
       data_file = trim(filename)//trim(zeros)//trim(data_file)//'.vtu'
    else
       data_file = trim(filename)//'.vtu'
    end if

    ewrite(2,*)('Writing file')
    ewrite(3,*)(data_file)

    allocate( eltypes(N_elements*(N_layers+1)) )
    allocate( elsizes(N_elements*(N_layers+1)) )

    call VTKOPEN(data_file,flength,' ',1)

    do i = 1, N_elements*(N_layers+1)
       eltypes(i) = 5
       elsizes(i) = 3
    end do

    !convert X,Y,D to dg elements

    n_dg=3*n_elements
    allocate( ones(n_dg*(N_layers+1)) )
    allocate( Xdg(n_dg*(N_layers+1) ) )
    allocate( Ydg(n_dg*(N_Layers+1) ) )
    allocate( Zdg(n_dg*(N_layers+1) ) )
    allocate( Ddg(n_dg*(N_Layers+1) ) )
    allocate( u1dg(n_dg*(N_layers+1) ) )
    allocate( u2dg(n_dg*(N_Layers+1) ) )
    allocate( m1dg(n_dg*(N_layers+1) ) )
    allocate( m2dg(n_dg*(N_Layers+1) ) )
    allocate( bottomdg(n_dg)  )
    allocate( evlist_out( (N_Layers+1)*n_dg) )

    if (present(n1)) then
       allocate( n1dg(n_dg*(N_layers+1) ) )
    end if
    if (present(n2)) then
       allocate( n2dg(n_dg*(N_Layers+1) ) )
    end if

    Zdg = 0.

    ones=1;

    ewrite(2,*)("element loop")


    print *, size(EVLIst_out)
    print *, 3*(n_layers+1)*n_elements
    forall (i=1:n_dg*(n_layers+1))
       EVList_out(i) = i
    end forall

    do ele = 1, N_elements

       d_ele=>mesh%EVList_h((ELE-1)*mesh%nh%loc+1:ELE*mesh%nh%loc)
       X_ele=>mesh%EVList_X((ELE-1)*3+1:ELE*3)
       u_ele=>mesh%EVList_u((ELE-1)*mesh%nu%loc+1:ELE*mesh%nu%loc)
       m_ele=>mesh%EVList_m((ELE-1)*mesh%nm%loc+1:ELE*mesh%nm%loc)

       do layer_i = 0, N_Layers-1

          dg_ele=>EVList_out(N_dg*layer_i + (ele-1)*3+1: &
               N_dg*layer_i + ele*3)

          XDg(dg_ELE) = mesh%x(x_ELE)
          yDg(dg_ELE) = mesh%Y(x_ELE)
          m1dg(dg_ele) = m(layer_i*2*n_moms+m_ele)
          m2dg(dg_ele) = m(layer_i*2*n_moms+n_moms+m_ele)

          u1dg(dg_ele) = u(layer_i*2*n_vels+u_ele(to_lin))
          u2dg(dg_ele) = u(layer_i*2*n_vels+n_vels+u_ele(to_lin))
          Ddg(dg_ele) = D(layer_i*n_dens+d_ele(to_lin))
          print *,Ddg(dg_ele)
          if (present(n1)) n1dg(dg_ele) = &
               n1(layer_i*n_dens+d_ele(to_lin))
          if (present(n2)) n2dg(dg_ele) = &
               n2(layer_i*n_dens+d_ele(to_lin))
       end do

       dg_ele=>EVList_out(N_dg*N_layers + (ele-1)*3+1:N_dg*N_layers + (ele)*3)

       XDg(dg_ELE) = mesh%x(x_ELE)
       yDg(dg_ELE) = mesh%Y(x_ELE)
       u1dg(dg_ele) = 0.0*dg_ele
       u2dg(dg_ele) = 0.0*dg_ele
       m1dg(dg_ele) = 0.0*dg_ele
       m2dg(dg_ele) = 0.0*dg_ele
        Ddg(dg_ele) = bottom(d_ele(to_lin))
       bottomdg(dg_ele - N_Layers*n_dg) = -Ddg(dg_ele)

       if (present(n1)) then
          n1dg(dg_ele) =  0.0*dg_ele
       end if
       if (present(n2)) then
          n2dg(dg_ele) =  0.0*dg_ele
       end if
    END do

    do layer_i = 0, N_Layers
       Zdg(N_dg * layer_i + 1: N_dg * (layer_i+1)) = -bottomdg
       do layer_j = layer_i, N_layers-1
          Zdg(N_dg * layer_i + 1:N_dg* (layer_i + 1)) = &
               Zdg(N_dg * layer_i + 1:N_dg* (layer_i + 1)) + &
               Ddg(N_dg * layer_j + 1: N_dg * (layer_j + 1))
       end do
    end do

!    Zdg(1:N_dg)=Ddg(1:N_dg)
!    Zdg(N_dg+1:2*N_dg)=-Ddg(N_dg+1:2*N_dg)

    do layer_i = 0,N_Layers
       Ddg(N_dg*layer_i+1:N_dg*(layer_i+1)) = &
            Ddg(N_dg*layer_i+1:N_dg*(layer_i+1)) - &
            sum(Ddg(N_dg*layer_i+1:N_dg*(layer_i+1)))/N_dg
    end do

    ewrite(2,*)("writing the mesh")

    call VTKWRITEMESHD(N_dg*(N_Layers+1),N_elements*(N_Layers+1), &
         xdg,ydg,Zdg,EVlist_out,eltypes,elsizes)
    call VTKWRITEDSN(Ddg,'Layer Depth',11)
    call VTKWRITEDSN(t*ones,'Time',11)
    call VTKWRITEDVN(U1dg,u2dg,0.0*U1dg,'Velocity',8)    
    call VTKWRITEDVN(m1dg,m2dg,0.0*U1dg,'Momentum',8)
    if (present(n1)) call VTKWRITEDVN(n1dg,n2dg,0.0*U1dg,'Pressure',8)

    call VTKCLOSE()

    !deallocate( eltypes )
    !deallocate( elsizes )
    !deallocate( XDg )
    !deallocate( YDg )
    !deallocate( Ddg )
    !deallocate( Zdg )
    !deallocate( EVList_out )

    ewrite(2,*)('end Subroutine dump_data_vtk_layerxs')

       

  end subroutine dump_data_vtk_layer

  function to_quad(field)
    real, dimension(:), intent (in) :: field
    real, dimension(6)  :: to_quad

    if (size(field)==3) then
       to_quad( (/ 1, 2, 3/) )=field
       to_quad(4)= 0.5* (field(1)+field(2))
       to_quad(5)= 0.5* (field(2)+field(3))
       to_quad(6)= 0.5* (field(1)+field(3))
    elseif (size(field)==6) then
       to_quad( (/ 1, 4 , 2 , 6, 5 , 3 /) )=field
    end if

  end function to_quad
    

end module vtk_io
