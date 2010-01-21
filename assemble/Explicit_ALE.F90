#include "confdefs.h"
#include "fdebug.h"
module ale_module
  use quadrature
  use elements
  use FLDebug
  use fldebug
  use fields
  use fefields,only: compute_lumped_mass
  use fetools, only: X_, Y_, Z_
  use spud
  use state_module
  use transform_elements
  use boundary_conditions
  use global_parameters, only: OPTION_PATH_LEN, dt
  use form_metric_field
  use metric_assemble
  use edge_length_module
  use solidconfiguration
  use fields_manipulation
  use eventcounter

  private
  real,              save :: toler_coarse,toler_fine,w1,w2,w3,w4,w5,fmin,minfchg,dx,dy,dz
  integer,           save :: numlspts,maxiter,condi(3)
  real,              save :: minch
  character*10,      save :: iter_string
  integer,           save :: fs
  type(mesh_type),   save :: mesh
  integer ,          save :: nonods,totele
  real ,      allocatable :: volumes(:)       
  integer                 :: nloc
  parameter(nloc=4)
  public explicit_ale,coordinate_update
contains

  subroutine explicit_ale(state,ii)    
    implicit none
    type(state_type), dimension(:), intent(inout) :: state
    integer, intent(in) :: ii
    !Local Memory
    real,    allocatable :: savebuffer(:)
    integer, allocatable :: ndlocate(:),nodelist(:)
    integer, allocatable :: negvolcenter_it(:),negvolgrad_it(:),failcenter_it(:),failgradient_it(:)
    integer, allocatable :: acc_crit_1(:),acc_crit_2(:)            
    !counters
    integer              :: negvolcenter,negvolgrad,failcenter,failgradient,counter1,counter2,nod_counter

    !indexes
    integer :: inod,inod1,node,lspt,i,j,it,its,icount
    !functional values and norms
    real    :: functional,functional1,functional2,functional3,functional4,functional5
    real    :: fnorm,f1norm,f2norm,f3norm,f4norm,f5norm,functional_buffer
    real    :: minf,minf1,minf2,minf3,minf4,minf5 
    real    :: mvf_nod_new,velu_nod_new,velv_nod_new,velw_nod_new,metric_nod_new(9)
    !for line search
    real    :: dnorm,step,dnorm1
    real    :: xnew,ynew,znew,xx2,yy2,zz2
    !gradient calc
    real    :: grad_f(3)
    !flags
    logical :: negvol,fail,found

    !fields
    type(scalar_field)  :: volumefraction,new_volumefraction,functional_field,dummy1,dummy2
    type(vector_field)  :: v,gv,coordinates,new_coordinates,new_velocity,dummy3
    type(tensor_field)  :: metric,new_metric,metric_save
    type(patch_type)    :: node_patch
    integer, dimension(:,:), allocatable :: b_condition
    character (len=OPTION_PATH_LEN) ale_path 

    !start switch
    integer, save :: start
    data start /1/

    !iteration counter
    integer, save:: iteration
    data iteration /0/
    fs=ii

#ifdef DOUBLEP
    !read input parameters from option path
    if(start==1) then
       ale_path='/mesh_adaptivity/mesh_movement/explicit_ale'
       call get_option(trim(ale_path)//'/number_of_linesearch_pts',numlspts)
       call get_option(trim(ale_path)//'/maximum_iterations',maxiter)                     
       call get_option(trim(ale_path)//'/functional1_weight',w1)  
       call get_option(trim(ale_path)//'/functional2_weight',w2)
       call get_option(trim(ale_path)//'/functional3_weight',w3)  
       call get_option(trim(ale_path)//'/functional4_weight',w4)
       call get_option(trim(ale_path)//'/functional5_weight',w5)
       call get_option(trim(ale_path)//'/minimum_val_functional_change',minfchg)
       call get_option(trim(ale_path)//'/minimum_val_functional',fmin)       
       call get_option(trim(ale_path)//'/coarse_tolerance',toler_coarse)       
       call get_option(trim(ale_path)//'/fine_tolerance',toler_fine)       
       call get_option(trim(ale_path)//'/gradient_dx',dx)  
       call get_option(trim(ale_path)//'/gradient_dy',dy)  
       call get_option(trim(ale_path)//'/gradient_dz',dz)             
       condi=0
       if (have_option(trim(ale_path)//'/move_nodes_in_x')) then
          condi(1)=1
       end if
       if (have_option(trim(ale_path)//'/move_nodes_in_y')) then
          condi(2)=1
       end if
       if (have_option(trim(ale_path)//'/move_nodes_in_z')) then
          condi(3)=1
       end if
       call i4_to_s_left (maxiter,iter_string)
    end if

    !Extract or allocate fields
    coordinates=extract_vector_field(state(fs),"Coordinate")
    v=extract_vector_field(state(fs),"Velocity")
    gv=extract_vector_field(state(fs),"GridVelocity")
    volumefraction=extract_scalar_field(state(fs),"MaterialVolumeFraction")  
    mesh=extract_mesh(state(fs),"CoordinateMesh")
    call allocate(metric,mesh=coordinates%mesh, name="Metric")
    call allocate(metric_save,mesh=coordinates%mesh,name="MetricSave")
    call allocate(functional_field, mesh=coordinates%mesh, name="Functional")
    functional_field%val=0.0
    call allocate(new_coordinates,dim=coordinates%dim,mesh=coordinates%mesh, name="NewCoordinate")
    call allocate(new_velocity,dim=coordinates%dim,mesh=coordinates%mesh, name="NewVelocity") 
    call allocate(new_volumefraction,mesh=coordinates%mesh, name="NewVolumeFraction")
    call allocate(new_metric,mesh=coordinates%mesh, name="NewMetric")
    allocate(negvolcenter_it(maxiter),negvolgrad_it(maxiter),failcenter_it(maxiter),failgradient_it(maxiter))
    allocate(acc_crit_1(maxiter),acc_crit_2(maxiter))
    allocate(savebuffer(17))
    allocate(volumes(element_count(coordinates)))
    call add_nelist(mesh)
    nonods=node_count(coordinates)
    totele=element_count(coordinates) 
    metric%val=0.0
    
    !Form the metric and make initial functional values
    call make_metric_for_functional1(state,metric)     
    !    call get_edge_lengths(metric, functional_field)
    
    !calculate element volumes for original mesh
    call calculate_element_volumes(state)

    !Make copies of necessary fields
    new_metric%val=metric%val                                !copy of metric
    new_volumefraction%val=volumefraction%val                !copy of volume fraction
    do i=1,3
       new_coordinates%val(i)%ptr=coordinates%val(i)%ptr     !copy of coordinates
       new_velocity%val(i)%ptr=v%val(i)%ptr           !copy of flow velocities
    end do

    !Create a node order list
    allocate(nodelist(nonods))
    allocate(ndlocate(nonods))
    do inod=1,nonods
       nodelist(inod)=inod
       ndlocate(inod)=inod
    end do
    
    !create boundary conditions for grid velocity
    allocate(b_condition(nonods,3))
    call create_boundary_conditions(state,b_condition)
    
    ewrite(0,*) 'calculating functional norms'
    !Calculate norm of functional field before calculations           
    call calculate_functional_norms(state,metric,metric_save,       &
         fnorm,f1norm,f2norm,f3norm,f4norm,f5norm,                  &
         minf,minf1,minf2,minf3,minf4,minf5,                        & 
         ndlocate,new_coordinates,new_velocity,new_volumefraction,new_metric,&
         functional_field)

    !Write out infinity norm of each functional Before iterations start
    it=0
    ewrite(0,'(a)')'values:        f        nodes'
    ewrite(0,'(a10,i0,e10.3,a1,i0)') 'iteration ',it,minval(functional_field%val),' ',0
    
    !Order nodelist according to functional scalar field values (maximum first)    
    call bubble_sort_i_list_from_r_values(nodelist,functional_field%val,nonods)
    
    !diagnostics
    dummy1=extract_scalar_field(state(fs),"FunctionalBegin")
    dummy2=extract_scalar_field(state(fs),"FunctionalIter")
    dummy3=extract_vector_field(state(fs),"FunctionalGradient")

    dummy1%val=functional_field%val
    dummy2%val=0.0
    dummy3%val(X_)%ptr=0.0
    dummy3%val(Y_)%ptr=0.0
    dummy3%val(Z_)%ptr=0.0
    
    negvolcenter_it=0
    negvolgrad_it=0   
    failcenter_it=0
    failgradient_it=0    
    acc_crit_1=0
    acc_crit_2=0
    
    !Begin iterations    
    do it=1,maxiter       
       !zero counters
       nod_counter=0;negvolcenter=0;negvolgrad=0;failcenter=0;failgradient=0;counter1=0;counter2=0
       !Begin nodal loop
       do inod1=1,nonods   
          inod=nodelist(inod1) !choose according to ordering in nodelist
          negvol=.false.; fail=.false.  !initialise failure switches
          
          !calculate gradients thorugh small displacements
          call calculate_functional_gradients(state,inod,b_condition,metric,&
               new_coordinates,new_velocity,new_metric,new_volumefraction,metric_save, &
               ndlocate,negvol,fail,grad_f)
          if(negvol) then              
             negvolgrad=negvolgrad+1
          elseif(fail) then
             failgradient=failgradient+1
          end if
          if(((grad_f(1)**2.0+grad_f(2)**2.0+grad_f(3)**2.0)==0.0).or.negvol.or.fail) cycle
          do i=1,3
             dummy3%val(i)%ptr(inod)=-grad_f(i)
          end do
          
          !Calculating the maximum length of the linesearch in the gradient direction
          !Two dnorms are calculated.  one is trying to take the length of the edge that more closely
          !aligned with the direction of the gradient (opposite direction.. of course).            
          dnorm=0.0
          node_patch = get_patch_node(mesh,inod)
          do icount=1,node_patch%count
             node = node_patch%elements(icount)
             if(node/=inod) then
                xx2=(new_coordinates%val(X_)%ptr(node) - new_coordinates%val(X_)%ptr(inod))**2.0
                yy2=(new_coordinates%val(Y_)%ptr(node) - new_coordinates%val(Y_)%ptr(inod))**2.0
                zz2=(new_coordinates%val(Z_)%ptr(node) - new_coordinates%val(Z_)%ptr(inod))**2.0
                dnorm1=(xx2+yy2+zz2)**0.5
                if(dnorm1>dnorm) dnorm=dnorm1                       
             end if
          end do
          deallocate(node_patch%elements)
          !divide the search length into steps.
          step=dnorm/numlspts
          dnorm=step        
          found=.false.
          functional_buffer=functional_field%val(inod)
          savebuffer=huge(0.0)
          !Start linesearch, save values corresponding to minimum
          do lspt=1,numlspts   
             xnew = new_coordinates%val(X_)%ptr(inod) - dnorm*grad_f(1) 
             ynew = new_coordinates%val(Y_)%ptr(inod) - dnorm*grad_f(2) 
             znew = new_coordinates%val(Z_)%ptr(inod) - dnorm*grad_f(3) 

             call calculate_functional_from_new_coordinates(                  &
                  state,inod,xnew,ynew,znew,metric,metric_save,               &
                  new_coordinates,new_velocity,new_metric,new_volumefraction, &
                  ndlocate,mvf_nod_new,velu_nod_new,velv_nod_new,velw_nod_new,&
                  metric_nod_new,negvol,fail,functional,                      &
                  functional1,functional2,functional3,functional4,functional5)                        
             
             dummy2%val(inod)=functional
             if(negvol.or.fail) exit
             if (functional<functional_buffer) then
                !Save corresponding values for node position in linesearch
                found=.true.
                functional_buffer=functional
                savebuffer(1)= functional             
                savebuffer(2)= xnew
                savebuffer(3)= ynew
                savebuffer(4)= znew  
                savebuffer(5)= mvf_nod_new
                savebuffer(6)= velu_nod_new
                savebuffer(7)= velv_nod_new
                savebuffer(8)= velw_nod_new
                do i=1,9
                   savebuffer(8+i)=metric_nod_new(i)
                end do
             end if
             dnorm=dnorm+step
          end do
          !skip node if linesearch is no good
          if(.not.found) cycle                              
          !this defines the first criteria for accepting position (need to work on second criteria)
          if ((functional_field%val(inod)-savebuffer(1))>minfchg) then
             functional_field%val(inod)=savebuffer(1)
             new_coordinates%val(X_)%ptr(inod)=savebuffer(2)
             new_coordinates%val(Y_)%ptr(inod)=savebuffer(3)
             new_coordinates%val(Z_)%ptr(inod)=savebuffer(4)
             new_volumefraction%val(inod)=savebuffer(5)
             new_velocity%val(X_)%ptr(inod)=savebuffer(6)
             new_velocity%val(Y_)%ptr(inod)=savebuffer(7)
             new_velocity%val(Z_)%ptr(inod)=savebuffer(8)
             new_metric%val(1,1,inod)=savebuffer(9)
             new_metric%val(1,2,inod)=savebuffer(10)
             new_metric%val(1,3,inod)=savebuffer(11)
             new_metric%val(2,1,inod)=savebuffer(12)
             new_metric%val(2,2,inod)=savebuffer(13)
             new_metric%val(2,3,inod)=savebuffer(14)
             new_metric%val(3,1,inod)=savebuffer(15)
             new_metric%val(3,2,inod)=savebuffer(16)
             new_metric%val(3,3,inod)=savebuffer(17)     
             counter1=counter1+1
             nod_counter=nod_counter+1
          end if
       end do    !end nodal loop

       !now re-sort nodelist according to the new functional_values after this iteration
       call bubble_sort_i_list_from_r_values(nodelist,functional_field%val,nonods)
       do i=1,nonods
          if(nodelist(nonods)<=0) then
             FLAbort('Nodelist less than zero!!!')
          end if
       end do
       ewrite(0,'(a10,i0,e10.3,a1,i0)') 'iteration ',it,minval(functional_field%val),' ',nod_counter
       !Negvol failures 
       negvolcenter_it(it)=negvolcenter
       negvolgrad_it(it)=negvolgrad
       !Fail to calculate
       failcenter_it(it)=failcenter
       failgradient_it(it)=failgradient
       !Acceptance criteria counters
       acc_crit_1(it)=counter1
       acc_crit_2(it)=counter2
       if(nod_counter.eq.0) exit !exit iteration loop if no nodes are moved in a given iteration.
    end do     !end iteration loop    

    ewrite(0,'(a)') 'Negative volumes'
    ewrite(0,'(a6,'//trim(iter_string)//'(A1,i0))') 'cent: ',(' ',negvolcenter_it(it),it=1,maxiter) 
    ewrite(0,'(a6,'//trim(iter_string)//'(A1,i0))') 'grad: ',(' ',negvolgrad_it(it),it=1,maxiter)
    ewrite(0,'(a)') 'Failures'
    ewrite(0,'(a6,'//trim(iter_string)//'(A1,i0))') 'cent: ',(' ',failcenter_it(it),it=1,maxiter) 
    ewrite(0,'(a6,'//trim(iter_string)//'(A1,i0))') 'grad: ',(' ',failgradient_it(it),it=1,maxiter)
    ewrite(0,'(a)') 'Accepted criterion'
    ewrite(0,'(a6,'//trim(iter_string)//'(A1,i0))') 'cr.1: ',(' ',acc_crit_1(it),it=1,maxiter) 
    ewrite(0,'(a6,'//trim(iter_string)//'(A1,i0))') 'cr.2: ',(' ',acc_crit_2(it),it=1,maxiter)     

    do i=1,3
       gv%val(i)%ptr=(new_coordinates%val(i)%ptr - coordinates%val(i)%ptr)/dt
    end do

    call deallocate(metric);call deallocate(metric_save);call deallocate(functional_field);
    call deallocate(new_coordinates);call deallocate(new_velocity);call deallocate(new_volumefraction)
    call deallocate(new_metric);deallocate(volumes)    
    deallocate(ndlocate);deallocate(nodelist);deallocate(savebuffer);deallocate(b_condition);
    deallocate(negvolgrad_it);deallocate(negvolcenter_it);deallocate(acc_crit_1);
    deallocate(failcenter_it);deallocate(failgradient_it);deallocate(acc_crit_2);
    START=0
#else
    FLAbort("Operation not supported in single precision.")
#endif
  end subroutine explicit_ale

    subroutine i4_to_s_left ( intval, s )
    !! I4_TO_S_LEFT converts an I4 to a left-justified string.
    !  Examples:
    !    Assume that S is 6 characters long:
    !    INTVAL  S
    !         1  1
    !        -1  -1
    !         0  0
    !      1952  1952
    !    123456  123456
    !   1234567  ******  <-- Not enough room!
    !  Modified:
    !    28 July 2000
    !  Author:
    !    John Burkardt
    !  Parameters:
    !    Input, integer INTVAL, an integer to be converted.
    !    Output, character ( len = * ) S, the representation of the integer.
    !    The integer will be left-justified.  If there is not enough space,
    !    the string will be filled with stars.
    implicit none
    character c
    integer i,idig,ihi,ilo,intval,ipos,ival
    character*10 s

    s = ' '
    ilo = 1
    ihi = len ( s )
    if ( ihi <= 0 ) then
       return
    end if
    !  Make a copy of the integer.
    ival = intval
    !  Handle the negative sign.
    if ( ival < 0 ) then
       if ( ihi <= 1 ) then
          s(1:1) = '*'
          return
       end if
       ival = -ival
       s(1:1) = '-'
       ilo = 2
    end if
    !  The absolute value of the integer goes into S(ILO:IHI).
    ipos = ihi
    !  Find the last digit of IVAL, strip it off, and stick it into the string.
    do
       idig = mod ( ival, 10 )
       ival = ival / 10
       if ( ipos < ilo ) then
          do i = 1, ihi
             s(i:i) = '*'
          end do
          return
       end if
       call digit_to_ch ( idig, c )
       s(ipos:ipos) = c
       ipos = ipos - 1
       if ( ival == 0 ) then
          exit
       end if
    end do
    !  Shift the string to the left.
    s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
    s(ilo+ihi-ipos:ihi) = ' '

    return
  end subroutine i4_to_s_left

  subroutine digit_to_ch ( digit, ch )
    !! DIGIT_TO_CH returns the character representation of a decimal digit.
    !  Discussion:
    !    Instead of CHAR, we now use the ACHAR function, which
    !    guarantees the ASCII collating sequence.
    !  Example:
    !    DIGIT   CH 
    !    -----  ---
    !      0    '0'
    !      1    '1'
    !    ...    ...
    !      9    '9'
    !     17    '*'
    !  Modified:
    !    04 August 1999
    !  Author:
    !    John Burkardt
    !  Parameters:
    !    Input, integer DIGIT, the digit value between 0 and 9.
    !    Output, character CH, the corresponding character.
    implicit none
    character ch
    integer digit
    if ( 0 <= digit .and. digit <= 9 ) then
       ch = achar ( digit + 48 )
    else
       ch = '*'
    end if
    return
  end subroutine digit_to_ch
  
  subroutine create_boundary_conditions(state,b_condition)
    implicit none
    type(state_type), dimension(:), intent(inout) :: state
    integer, intent(inout),dimension(nonods,3):: b_condition
    type(vector_field) :: v
    character(len=FIELD_NAME_LEN):: bctype
    integer :: sele,i,k,n,snodes(3)
    logical, dimension(3):: applies
    integer, dimension(:), pointer:: surface_element_list, surface_node_list
    b_condition=1
    
    v=extract_vector_field(state(fs),"Velocity")
    do i=1, get_boundary_condition_count(v)
       call get_boundary_condition(v, i, &
           type=bctype, &
           surface_node_list=surface_node_list, &
           surface_element_list=surface_element_list, &
           applies=applies)
       if (bctype=='dirichlet') then
          
          do k=1, size(surface_element_list)
             do n=1, v%dim     
                if (applies(n).or.condi(n)==0) then
                   sele=surface_element_list(k)
                   snodes=face_global_nodes(v%mesh, sele)
                   b_condition(snodes(1),n)=0
                   b_condition(snodes(2),n)=0
                   b_condition(snodes(3),n)=0                 
                end if
             end do
          end do
       end if
    end do
    
  end subroutine create_boundary_conditions
  
  subroutine make_metric_for_functional1(state,metric)
    implicit none   
    type(state_type), dimension(:), intent(inout) :: state
    type(tensor_field), intent(inout) ::metric
    real :: MESTP2
    type (vector_field) positions
    CHARACTER (len=OPTION_PATH_LEN) adapt_path

    positions=extract_vector_field(state(fs),"Coordinate")    
    adapt_path='/mesh_adaptivity/hr_adaptivity'
    call assemble_metric(state, metric)
  end subroutine make_metric_for_functional1

  subroutine calculate_element_volumes(state)
    !this subroutine calculates elemental volumes of the CoordinateMesh
    !and stores them on volumes(nonods)
    implicit none
    type(state_type), dimension(:), intent(in) :: state
    type(vector_field) :: coordinates
    integer ele,iloc,node,ilevel
    real xx(4),yy(4),zz(4),volume    
    coordinates=extract_vector_field(state(fs),"Coordinate")
    do ele=1,element_count(coordinates)
       do iloc=1,nloc
          node = coordinates%mesh%ndglno((ele-1)*nloc+iloc)
          xx(iloc) = coordinates%val(X_)%ptr(node)
          yy(iloc) = coordinates%val(Y_)%ptr(node)
          zz(iloc) = coordinates%val(Z_)%ptr(node)
       end do
       volumes(ele)=element_volume_s(nloc,xx,yy,zz)
       if (volumes(ele)<=0.0) then
          FLAbort('Negative volume on original mesh???')
       end if
    end do
  end subroutine calculate_element_volumes

  subroutine calculate_functional_norms(state,metric,metric_save, &
       fnorm,f1norm,f2norm,f3norm,f4norm,f5norm,                  &
       minf,minf1,minf2,minf3,minf4,minf5,                        & 
       ndlocate,new_coordinates,new_velocity,new_volumefraction,new_metric,&
       functional_field)
    type(state_type), dimension(:), intent(inout) :: state
    real, intent(inout) :: fnorm,f1norm,f2norm,f3norm,f4norm,f5norm
    real, intent(inout) :: minf,minf1,minf2,minf3,minf4,minf5
    type(vector_field), intent(in) :: new_coordinates,new_velocity
    type(scalar_field), intent(in) :: new_volumefraction
    type(tensor_field), intent(in) :: new_metric,metric    
    type(tensor_field) , intent(inout) :: metric_save
    type(scalar_field), intent(inout) :: functional_field
    integer,          intent(inout):: ndlocate(new_coordinates%mesh%nodes)
    real    :: xnew,ynew,znew
    logical :: negvol,fail
    real    :: mvf_nod_new,velu_nod_new,velv_nod_new,velw_nod_new,metric_nod_new(9)
    real    :: functional,functional1,functional2,functional3,functional4,functional5
    integer :: i,j

    fnorm=0.0;f1norm=0.;f2norm=0.;f3norm=0.;f4norm=0.;f5norm=0.;
    minf=huge(0.0);minf1=huge(0.0);minf2=huge(0.0);minf3=huge(0.0);minf4=huge(0.0);minf5=huge(0.0);minf=huge(0.0);
    do inod=1,nonods                       
       negvol=.false.
       fail=.false.       
       xnew=new_coordinates%val(X_)%ptr(inod)
       ynew=new_coordinates%val(Y_)%ptr(inod)
       znew=new_coordinates%val(Z_)%ptr(inod)
       call calculate_functional_from_new_coordinates(                  &
            state,inod,xnew,ynew,znew,metric,metric_save,               &
            new_coordinates,new_velocity,new_metric,new_volumefraction, &
            ndlocate,mvf_nod_new,velu_nod_new,velv_nod_new,velw_nod_new,&
            metric_nod_new,negvol,fail,functional,                      &
            functional1,functional2,functional3,functional4,functional5)             
       if(negvol.or.fail) cycle
       if(functional>fnorm) fnorm=functional
       if(functional1>f1norm) f1norm=functional1
       if(functional2>f2norm) f2norm=functional2
       if(functional3>f3norm) f3norm=functional3
       if(functional4>f4norm) f4norm=functional4
       if(functional5>f5norm) f5norm=functional5
       if(functional<minf) minf=functional
       if(functional1<minf1) minf1=functional1
       if(functional2<minf2) minf2=functional2
       if(functional3<minf3) minf3=functional3
       if(functional4<minf4) minf4=functional4
       if(functional5<minf5) minf5=functional5
       functional_field%val(inod)=functional
    end do
  end subroutine calculate_functional_norms

  subroutine calculate_functional_from_new_coordinates(            &
       state,inod,xnew,ynew,znew,metric,metric_save,               &
       new_coordinates,new_velocity,new_metric,new_volumefraction, &
       ndlocate,mvf_nod_new,velu_nod_new,velv_nod_new,velw_nod_new,&
       metric_nod_new,negvol,fail,functional,                      &
       functional1,functional2,functional3,functional4,functional5)    

    implicit none
    type(state_type), dimension(:), intent(inout) :: state    
    integer,               intent(in) :: inod
    real,                  intent(in) :: xnew,ynew,znew
    type(vector_field),    intent(in) :: new_coordinates,new_velocity
    type(tensor_field),    intent(in) :: new_metric,metric
    type(tensor_field), intent(inout) :: metric_save
    type(scalar_field),    intent(in) :: new_volumefraction 
    integer,            intent(inout) :: ndlocate(new_coordinates%mesh%nodes)
    real,               intent(inout) :: mvf_nod_new,velu_nod_new,velv_nod_new,velw_nod_new
    real,               intent(inout) :: metric_nod_new(9)
    logical,            intent(inout) :: negvol,fail
    real,               intent(inout) :: functional,functional1,functional2
    real,               intent(inout) :: functional3,functional4,functional5
    ! Local Memory
    real    :: uu,vv,ww,mele(9),q,shape,size1,xx(4),yy(4),zz(4),functional1temp
    integer :: i, j, dim,ele,iloc,node,icount
    real, dimension(mesh_dim(new_metric%mesh), mesh_dim(new_metric%mesh)) :: evectors
    real, dimension(mesh_dim(new_metric%mesh)) :: evalues, desired_lengths    
    type(vector_field) v,coordinates
    type(patch_type) patch
    dim = mesh_dim(new_metric%mesh)
    functional=0;functional1=0;functional2=0;functional3=0;functional4=0;functional5=0;
    !check that there is no negative volume produced by the new coordinates.
    !if there is one element in the vicinity of the current node that is turned
    !inside out, then the subroutine returns 0.0 values for all functionals
    !and a negvol=.true..   
    call check_for_negative_volume(state,inod,xnew,ynew,znew,new_coordinates,negvol)     
    if(negvol) return    

    !Second we need to find out in which element from the original mesh 
    !the node falls into, interpolate the values that correspond to it,
    !and change the node reference for faster search at later iterations.
    !UNOD,VNOD,WNOD,CNOD, and MALENode, contain all the interpolated values.    
    call interpolate_new_values(state,inod,xnew,ynew,znew,metric, &
         mvf_nod_new,velu_nod_new,velv_nod_new,velw_nod_new,      &
         metric_nod_new,ndlocate,fail) 
    if (fail) return

    !Now that we have all the interpolated values, we continue to calculate the functionals 
    !-----------------------------------------------------------------------------
    !This section calculates the nodal element quality functional
    mele=0.0
    patch=get_patch_ele(mesh,inod,level=1)
    functional1=0.0
    patchloop: do icount=1,patch%count ! loop over elements around node
       ele=patch%elements(icount)
       do iloc=1,nloc
          node = new_coordinates%mesh%ndglno((ele-1)*nloc+iloc)
          if(inod==node) then
             xx(iloc) = xnew
             yy(iloc) = ynew
             zz(iloc) = znew
             mele = mele + metric_nod_new 
             mele=mele+metric_nod_new
          else
             xx(iloc) = new_coordinates%val(X_)%ptr(node)
             yy(iloc) = new_coordinates%val(Y_)%ptr(node)
             zz(iloc) = new_coordinates%val(Z_)%ptr(node)
             !             mele=mele+new_metric%val(:,:,node)
             mele(1)=mele(1)+new_metric%val(1,1,node)
             mele(2)=mele(2)+new_metric%val(1,2,node)
             mele(3)=mele(3)+new_metric%val(1,3,node)
             mele(4)=mele(4)+new_metric%val(2,1,node)
             mele(5)=mele(5)+new_metric%val(2,2,node)
             mele(6)=mele(6)+new_metric%val(2,3,node)
             mele(7)=mele(7)+new_metric%val(3,1,node)
             mele(8)=mele(8)+new_metric%val(3,2,node)
             mele(9)=mele(9)+new_metric%val(3,3,node)
          end if          
       end do   
       mele=mele/4.0 !average the metric to obtain an elemental value.
       call elemqual(xx,yy,zz,mele,functional1temp,q,shape,size1)
       if(functional1temp>functional1) functional1=functional1temp
    end do patchloop
    deallocate(patch%elements)
    
    v=extract_vector_field(state(fs),"Velocity")
    coordinates=extract_vector_field(state(fs),"Coordinate")
    !------------------------------------------------------------------------------
    !This section calculates the lagrangian functional.(funtional 2)
    uu=v%val(X_)%ptr(inod)-(xnew-coordinates%val(X_)%ptr(inod))/dt
    vv=v%val(Y_)%ptr(inod)-(ynew-coordinates%val(Y_)%ptr(inod))/dt
    ww=v%val(Z_)%ptr(inod)-(znew-coordinates%val(Z_)%ptr(inod))/dt
    functional2 = ((dt/minch)**2.0)*(uu**2.0 + vv**2.0 + ww**2.0)

    !------------------------------------------------------------------------------
    !!     This section calculates the third functional
    !    call allocate(lumped_mass_matrix, coordinates%mesh, "Lumped mass matrix")
    !    call allocate(gradient,3,coordinates%mesh, "Gradient")    
    !    call add_nelist(coordinates%mesh)
    !    call zero(lumped_mass_matrix)
    !    call zero(gradient)    
    !    t_shape => ele_shape(volumefraction, 1)
    !    x_shape => ele_shape(positions, 1)
    !
    !    ! First, compute gradient and mass matrix.
    !    do ele=1,element_count(infield)
    !       ! Compute detwei.
    !       call transform_to_physical(ele_val(positions, ele), x_shape, t_shape, dm_t=dt_t, detwei=detwei)       
    !       r = shape_dshape(t_shape, dt_t, detwei)
    !       r_grad_ele = tensormul(r, ele_val(infield, ele), 3)       
    !       call addto(gradient, ele_nodes(infield, ele), r_grad_ele)       
    !       ! Lump the mass matrix
    !       mass_matrix = shape_shape(t_shape, t_shape, detwei)
    !       call addto(lumped_mass_matrix, ele_nodes(infield, ele), sum(mass_matrix, 2))
    !    end do
    !    !THIS IS FROM FIELD DERIVATIVES.F90 LINE 1076
    !    do i=1,3
    !       gradient%val(i)%ptr(node) = gradient%val(i)%ptr(node) / node_val(lumped_mass_matrix, node)
    !    end do
    !    call deallocate(lumped_mass_matrix)
    !    call deallocate(gradient)
    functional3=0.0
    functional4=0.0
    functional5=0.0

    functional=w1*functional1+w2*functional2+w3*functional3+w4*functional4+w5*functional5
  end subroutine calculate_functional_from_new_coordinates

  subroutine bubble_sort_i_list_from_r_values(ilist,rrlist,nlist)
    implicit none
    integer, intent(in)    :: nlist
    integer, intent(inout) :: ilist(nlist)
    real,    intent(in)    :: rrlist(nlist)
    integer i,j,ii
    real rii,rlist(nlist)
    !sorts rlist in descending order and changes values of ilist accordingly.
    rlist=rrlist
    do i=1,nlist
       do j=2,nlist
          if(rlist(j-1).lt.rlist(j)) then
             rii=rlist(j-1)
             rlist(j-1)=rlist(j)
             rlist(j)=rii

             ii=ilist(j-1)
             ilist(j-1)=ilist(j)
             ilist(j)=ii
          endif
       end do
    end do
  end subroutine bubble_sort_i_list_from_r_values
  
  subroutine calculate_functional_gradients(state,inod,b_condition,metric, &
       new_coordinates,new_velocity,new_metric,new_volumefraction,metric_save,  &
       ndlocate,negvol,fail,grad_f)

    implicit none
    type(state_type), dimension(:), intent(inout) :: state 
    type(vector_field) , intent(in) :: new_coordinates,new_velocity
    integer,             intent(in) :: b_condition(nonods,3)
    type(tensor_field) , intent(in) :: new_metric,metric
    type(scalar_field) , intent(in) :: new_volumefraction
    type(tensor_field) , intent(inout) :: metric_save
    integer, intent(in)    ::  inod
    integer, intent(inout) ::  ndlocate(new_velocity%mesh%nodes)
    logical, intent(inout) ::  negvol,fail
    real,    intent(inout) ::  grad_f(3)
    
    
    !local memory
    integer i,j
    real :: f(6),f1,f2,f3,f4,f5,rdf,d(3),new_c(3)
    real :: mvf_nod_new,velu_nod_new,velv_nod_new,velw_nod_new,metric_nod_new(9)
    
    d(1)=dx;d(2)=dy;d(3)=dz;

    !zero the functional and gradient
    f=0.0  
    grad_f=0.0
    do i=1,3       
       new_c(1)=new_coordinates%val(X_)%ptr(inod) 
       new_c(2)=new_coordinates%val(Y_)%ptr(inod)
       new_c(3)=new_coordinates%val(Z_)%ptr(inod) 
       !if there is no boundary condition that impedes movement of nodes in this direction, for this node (inod)
       !then calculate the functionals. Otherwise, jump to next coordinate axis.
       if(b_condition(inod,i)==0) cycle
       new_c(i)=new_coordinates%val(i)%ptr(inod) + d(i)       
       
       !Calculate new functional values (positive side)
       call calculate_functional_from_new_coordinates(                  &
            state,inod,new_c(1),new_c(2),new_c(3),metric,metric_save,   &
            new_coordinates,new_velocity,new_metric,new_volumefraction, &
            ndlocate,mvf_nod_new,velu_nod_new,velv_nod_new,velw_nod_new,&
            metric_nod_new,negvol,fail,f(i),                            &
            f1,f2,f3,f4,f5)
       !If there is a negative volume resulting from this, or a failure to calculate, then exit the i loop
       if(negvol.or.fail) return                                                                                 
       
       new_c(1)=new_coordinates%val(X_)%ptr(inod) 
       new_c(2)=new_coordinates%val(Y_)%ptr(inod)
       new_c(3)=new_coordinates%val(Z_)%ptr(inod)     
       new_c(i)=new_coordinates%val(i)%ptr(inod) - d(i)
       !Calculate new functional values. (negative side
       call calculate_functional_from_new_coordinates(                  &
            state,inod,new_c(1),new_c(2),new_c(3),metric,metric_save,   &
            new_coordinates,new_velocity,new_metric,new_volumefraction, &
            ndlocate,mvf_nod_new,velu_nod_new,velv_nod_new,velw_nod_new,&
            metric_nod_new,negvol,fail,f(i+3),                          &
            f1,f2,f3,f4,f5)
       !If there is a negative volume resulting from this, or a failure to calculate, then exit the i loop
       if(negvol.or.fail) return                                  
       grad_f(i)=(f(i)-f(i+3))/(2*d(i))                     
    end do              ! end i loop
    
    !Now will calculate if its norm is at least of a given tolerance. 
    rdf=(grad_f(1)**2.0 + grad_f(2)**2.0 + grad_f(3)**2.0)**0.5
    if (rdf.le.toler_fine) then
       grad_f=0.0
    else
       !normalizing gradient vector
       grad_f=grad_f/rdf
    end if
  end subroutine calculate_functional_gradients
  
  subroutine coordinate_update(state)
    implicit none
    type(state_type), dimension(:), intent(inout) :: state
    integer i
    type(vector_field) coordinates,gv

    call IncrementEventCounter(EVENT_MESH_MOVEMENT)

    coordinates=extract_vector_field(state(fs),"Coordinate")
    gv=extract_vector_field(state(fs),"GridVelocity")

    do i=1,3
       coordinates%val(i)%ptr=coordinates%val(i)%ptr + gv%val(i)%ptr*dt
    end do

  end  subroutine coordinate_update

  subroutine bubble_sort_integer_list(list,nlist)
    implicit none
    integer nlist,list(nlist)
    integer i,j,ii
    !sorts list in descending order
    do i=1,nlist
       do j=2,nlist
          if(list(j-1).lt.list(j)) then
             ii=list(j-1)
             list(j-1)=list(j)
             list(j)=ii
          endif
       end do
    end do
  end subroutine bubble_sort_integer_list

  subroutine check_for_negative_volume(state,inod,xnew,ynew,znew,new_coordinates,negvol)
    !Checks for negative volume in the new node configuration.
    !(i.e.: checks that no element is turned inside out by the movements)
    implicit none
    type(state_type), dimension(:), intent(inout) :: state
    integer,            intent(in)    :: inod
    real,               intent(in)    :: xnew,ynew,znew
    logical,            intent(inout) :: negvol    
    type(vector_field), intent(in)    :: new_coordinates
    type(vector_field)                :: coordinates
    type(patch_type) patch
    integer icount,ele,iloc,node
    real xx(4),yy(4),zz(4),volume   
    
    coordinates=extract_vector_field(state(fs),"Coordinate")
    patch=get_patch_ele(mesh,inod,level=1)
    patchloop: do icount=1,patch%count ! loop over elements around node
       ele=patch%elements(icount)
       do iloc=1,nloc
          node = new_coordinates%mesh%ndglno((ele-1)*nloc+iloc)
          xx(iloc) = new_coordinates%val(X_)%ptr(node)
          yy(iloc) = new_coordinates%val(Y_)%ptr(node)
          zz(iloc) = new_coordinates%val(Z_)%ptr(node)
          if(inod==node) then
             xx(iloc) = xnew
             yy(iloc) = ynew
             zz(iloc) = znew
          end if
       end do
       volume=element_volume_s(nloc,xx,yy,zz)
       if (volume<=0.0) then
          negvol=.true.
          return
       end if
    end do patchloop
    deallocate(patch%elements)
  end subroutine check_for_negative_volume
  
  subroutine interpolate_new_values(state,inod,xnew,ynew,znew,metric,  &
       mvf_nod_new,velu_nod_new,velv_nod_new,velw_nod_new,      &
       metric_nod_new,ndlocate,fail) 
    !Checks for negative volume in the new node configuration
    implicit none
    type(state_type), dimension(:), intent(inout) :: state
    integer,            intent(in)    :: inod
    real,               intent(in)    :: xnew,ynew,znew
    type(tensor_field), intent(in)    :: metric
    real,               intent(inout) :: mvf_nod_new,velu_nod_new,velv_nod_new,velw_nod_new
    real,               intent(inout) :: metric_nod_new(9)
    integer,            intent(inout) :: ndlocate(nonods)
    logical,            intent(inout) :: fail   
    type(patch_type) patch
    type(scalar_field) volumefraction
    type(vector_field) v,coordinates
    integer icount,ele,iloc,node,ilevel
    real xx(4),yy(4),zz(4),xx1,yy1,zz1,pvl(4),volume,vcheck1,shapef
    logical found
    found=.false.
    fail=.false.
    !Extract fields to interpolate
    v=extract_vector_field(state(fs),"Velocity")
    volumefraction=extract_scalar_field(state(fs),"MaterialVolumeFraction")  
    coordinates=extract_vector_field(state(fs),"Coordinate")
    
    levelloop : do ilevel=1,2
       patch = get_patch_ele(mesh,ndlocate(inod),level=ilevel)    ! form patch           
       do icount=1,patch%count ! loop over elements around node in ndlocate
          ele=patch%elements(icount)
          volume=volumes(ele)
          do iloc=1,nloc
             node = coordinates%mesh%ndglno((ele-1)*nloc+iloc)
             xx(iloc) = coordinates%val(X_)%ptr(node)
             yy(iloc) = coordinates%val(Y_)%ptr(node)
             zz(iloc) = coordinates%val(Z_)%ptr(node)
          end do
          !calculate partial volumes
          do iloc=1,nloc
             node = coordinates%mesh%ndglno((ele-1)*nloc+iloc)
             xx1=xx(iloc)
             yy1=yy(iloc)
             zz1=zz(iloc)
             xx(iloc)=xnew
             yy(iloc)=ynew
             zz(iloc)=znew 
             pvl(iloc)=element_volume_s(nloc,xx,yy,zz)
             xx(iloc)=xx1  
             yy(iloc)=yy1
             zz(iloc)=zz1            
          end do
          vcheck1=1.0-(abs(pvl(1))+abs(pvl(2))+abs(pvl(3))+abs(pvl(4)))/abs(volume)
          if (abs(vcheck1).lt.toler_coarse) then  
             found=.true.
             mvf_nod_new =0.0
             velu_nod_new =0.0
             velv_nod_new =0.0
             velw_nod_new =0.0
             metric_nod_new=0.0                          
             do iloc=1,nloc
                node=coordinates%mesh%ndglno((ele-1)*nloc+iloc)
                shapef=pvl(iloc)/volume
                mvf_nod_new=mvf_nod_new + volumefraction%val(node)*shapef
                velu_nod_new= velu_nod_new + v%val(1)%ptr(node)*shapef
                velv_nod_new= velv_nod_new + v%val(2)%ptr(node)*shapef
                velw_nod_new= velw_nod_new + v%val(3)%ptr(node)*shapef
                metric_nod_new(1)=metric_nod_new(1) + metric%val(1,1,node)*shapef   
                metric_nod_new(2)=metric_nod_new(2) + metric%val(1,2,node)*shapef
                metric_nod_new(3)=metric_nod_new(3) + metric%val(1,3,node)*shapef 
                metric_nod_new(4)=metric_nod_new(4) + metric%val(2,1,node)*shapef
                metric_nod_new(5)=metric_nod_new(5) + metric%val(2,2,node)*shapef
                metric_nod_new(6)=metric_nod_new(6) + metric%val(2,3,node)*shapef
                metric_nod_new(7)=metric_nod_new(7) + metric%val(3,1,node)*shapef
                metric_nod_new(8)=metric_nod_new(8) + metric%val(3,2,node)*shapef  
                metric_nod_new(9)=metric_nod_new(9) + metric%val(3,3,node)*shapef  
             end do
             exit levelloop
          end if
       end do
       deallocate(patch%elements)
    end do levelloop
    !if on level 2, grab an element from that level and choose that node. That way there is never 
    !more than two levels to search
    if(found.and.(ilevel==2)) then
       ndlocate(inod)=coordinates%mesh%ndglno((ele-1)*nloc+1) !grab an node 
    end if
    if(ndlocate(inod)==0) then
       !       ewrite(-1,*) 'ndlocate(inod): ',ndlocate(inod),inod
       FLAbort('Something is wrong with ndlocate')
    end if
    if(.not.found) fail=.true.    
  end subroutine interpolate_new_values

  subroutine elemqual(x,y,z,m,fg,q,shape,size)
    implicit none
    logical oldmeth
    parameter(oldmeth=.true.)
    real alpha,shape,size
    integer i
    real x(4), y(4), z(4), vol, areas(4), r, m(3,3), q, l(6),fg,s,d,e
    alpha=4.898979486
    call mtetin( x, y, z, m, vol, areas, l, r, q )      
    !     write(*,*) 'out of mtetin'
    !     what it returns is as follows:
    !     vol: the element volume
    !     areas: the areas of the four faces
    !     l: the lengths of the six edges
    !     r: the in-sphere radius
    !     q: some measure of the 'quality' of the shape of the element
    !     (note that q is not the functional, since it is independent of the element size).    
    !     all of the above are worked out in metric (m) space.
    !     the functional is worked out in two parts; the first part (called s below) 
    !     is from the in-sphere radius (note that i've capped it below at 1e-12, below which 
    !     it gives some upper bound for the functional, which i defined as 1e+32):    
    s = abs(r)*alpha
    if( s .gt. 1e-12 ) then
       s = 1.0 - 1.0/s
       s = s*s
    else
       s = 1e+32
    end if    
    !     alpha is a parameter defined as 4.898979486 [=sqrt(24)], which is the inverse of the ideal in-sphere radius.
    !     the second part (called e below) is from the edges (again, note that it's capped to prevent any floating overflow):
    !     e = 0.0
    !     do i = 1, 6
    !     d = 1.0 - l(i)
    !     if( abs(d) .lt. 1e+12 ) then
    !     e = e + d*d
    !     else
    !     e = 1e+32
    !     end if
    !     end do
    !     the original adaptivity combines these to form the functional:
    !     f = s + 0.5*e
    !     however, i found this can cause problems under certain conditions where there are edges constrained by the geometry of the problem 
    !     - it can force the in-sphere to be much smaller than ideal, making s very large. this leads to adaptivity trying to compensate by 
    !     making the edges very long (this slightly increases the in-sphere radius, which brings down the functional quite quickly, whereas the 
    !     edge length increase does not increase the functional quite so rapidly).
    !     i decided to modify the functional to take account of the fact that it needs to increase much more rapidly as the edge length increases.
    !     what i do is to 'hack' the contribution from the edges so that it increases much more rapidly as the edge-length increases 
    !     (not pretty, but it does the job), so first of all it's the same as before for the edges:
    
    e = 0.0
    do i = 1, 6
       d = 1.0 - l(i)
       if( abs(d) .lt. 1e+12 ) then
          e = e + d*d
       else
          e = 1e+32
       end if
    end do
    e = 0.5*e
    !     ...then the extra hack to make sure the contribution, e, increases more quickly as the lengths increase:    
    if(.not.oldmeth) then
       if( e .gt. 1e+6 ) then
          e = 1e+32
       else if( e .gt. 6.0 ) then
          e = e + e*e*e*e - 1102.0
       else if( e .gt. 3.0 ) then
          e = e + e*e*e - 22.0
       else if( e .gt. 2.0 ) then
          e = e + e*e - 4.0
       end if
    end if
    !     finally, add together the in-sphere and lengths contributions to form the functional:
    fg = s + e
    shape=s
    size=e
    return
  end subroutine elemqual

end module ale_module
