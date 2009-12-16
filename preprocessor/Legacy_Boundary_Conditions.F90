#include "fdebug.h"

module legacy_boundary_conditions
use state_module
use fields
use boundary_conditions
use spud
use global_parameters, only: OPTION_PATH_LEN, phase2state_index
use legacy_field_lists
implicit none

  ! string that identifies the legacy dirichlet b.c. types
  character(len=*), parameter:: DIRICHLET_BC_TYPE='dirichlet'
  character(len=*), parameter:: SEM_BC_TYPE='synthetic_eddy_method'
    ! string that identifies the legacy robin b.c. types 
  ! (for filling in salphe and taire)
  character(len=*), parameter:: ROBIN_BC_TYPE='robin'
  ! string that identifies the surface field that stores the values
  ! for dirichlet b.c.s (also for SEM bcs)
  character(len=*), parameter:: DIRICHLET_VALUE_ID='value'
  ! strings that identifies the surface fields that store the values
  ! for robin b.c.s corresponding to salphe and taire resp.
  character(len=*), parameter:: ROBIN_SALPHE_VALUE_ID='order_one_coefficient'
  character(len=*), parameter:: ROBIN_TAIRE_VALUE_ID='order_zero_coefficient'
  ! strings that identify the normal and tangent vectors for
  ! 'rotated' boundary conditions:
  character(len=*), parameter:: NORMAL_VECTOR_ID='normal'
  character(len=*), parameter:: TANGENT1_VECTOR_ID='tangent1'
  character(len=*), parameter:: TANGENT2_VECTOR_ID='tangent2'
  
  private
  
  public getbcuvw_sizes, getbcuvw12, getbct_size, getbct, &
    getsalphetaire, getnnodro, getnodrot, &
    copy_back_tangents, check_velocity_bcs, set_rotat_gettan, &
    set_legacy_boundary_values

contains
  
subroutine getbcuvw_sizes(velocity_field, bcu_size, bcv_size, bcw_size)
!! get number of legacy dirichlet boundary conditions for each component
type(vector_field), intent(in):: velocity_field
integer, intent(out):: bcu_size, bcv_size, bcw_size

  integer, dimension(:), pointer:: surface_node_list
  character(len=FIELD_NAME_LEN) type
  logical applies(3)
  integer i
  
  bcu_size=0
  bcv_size=0
  bcw_size=0
  ! loop over the boundary conditions
  do i=1, get_boundary_condition_count(velocity_field)
    call get_boundary_condition(velocity_field, i, type=type, &
      applies=applies, surface_node_list=surface_node_list)
    select case (velocity_field%dim)
    case (1)
      applies(2:3)=.false.
    case (2)
      applies(3)=.false.
    case (3)
      ! do nothing
    case default
      ewrite(0,*) "Vector field of dimension", velocity_field%dim
      FLAbort("What the * is going on?")
    end select
    if (type==DIRICHLET_BC_TYPE .or. type==SEM_BC_TYPE) then
      if (has_surface_field(velocity_field, i, NORMAL_VECTOR_ID)) then
          ! rotated (surface aligned) b.c.'s
         if (applies(1)) bcw_size=bcw_size+size(surface_node_list)
         if (applies(2)) bcu_size=bcu_size+size(surface_node_list)
         if (applies(3)) bcv_size=bcv_size+size(surface_node_list)
      else
         if (applies(1)) bcu_size=bcu_size+size(surface_node_list)
         if (applies(2)) bcv_size=bcv_size+size(surface_node_list)
         if (applies(3)) bcw_size=bcw_size+size(surface_node_list)
      end if
    end if
  end do
  
end subroutine getbcuvw_sizes

subroutine getbcuvw12(velocity_field,  bcu1,bcv1, bcw1, bcu2, bcv2, bcw2)
!! get legacy dirichlet boundary conditions bc[u|v|w][1|2]
type(vector_field), intent(in):: velocity_field
real, dimension(:), intent(out):: bcu1,bcv1, bcw1
integer, dimension(:), intent(out):: bcu2, bcv2, bcw2

  type(vector_field), pointer :: surface_field
  integer, dimension(:), pointer :: surface_node_list
  character(len=FIELD_NAME_LEN) type
  real value(velocity_field%dim)
  logical applies(3), rotat
  integer i, j, nod
  integer bcui, bcvi, bcwi
  
  ! indices in the arrays:
  bcui=0
  bcvi=0
  bcwi=0
  
  ! loop over the boundary conditions
  do i=1, get_boundary_condition_count(velocity_field)
    call get_boundary_condition(velocity_field, i, type=type, &
      applies=applies, surface_node_list=surface_node_list)
    if (type==DIRICHLET_BC_TYPE .or. type==SEM_BC_TYPE) then

      ! surface field containing b.c. values
      surface_field => extract_surface_field(velocity_field, i, DIRICHLET_VALUE_ID)
      
      ! whether we have rotated (surface aligned) b.c.'s
      rotat=has_surface_field(velocity_field, i, NORMAL_VECTOR_ID)
      
      ! adjust applies for unused dimensions
      select case (velocity_field%dim)
      case (1)
        applies(2:3)=.false.
      case (2)
        applies(3)=.false.
      case (3)
        ! do nothing
      case default
        ewrite(0,*) "Vector field of dimension", velocity_field%dim
        FLAbort("What the * is going on?")
      end select
      
      ! loop over nodes to which this b.c. applies
      do j=1, size(surface_node_list)
        value=node_val( surface_field, j )
        nod=surface_node_list(j)
        if (rotat) then
           if (applies(1)) then
             ! normal component
             bcwi=bcwi+1  
             bcw1(bcwi)=value(1)
             bcw2(bcwi)=nod
           end if
           if (applies(2)) then
             ! tangent component 1
             bcui=bcui+1  
             bcu1(bcui)=value(2)
             bcu2(bcui)=nod
           end if
           if (applies(3)) then
             ! tangent component 1
             bcvi=bcvi+1  
             bcv1(bcvi)=value(3)
             bcv2(bcvi)=nod
           end if
        else
           if (applies(1)) then
             bcui=bcui+1  
             bcu1(bcui)=value(1)
             bcu2(bcui)=nod
           end if
           if (applies(2)) then
             bcvi=bcvi+1  
             bcv1(bcvi)=value(2)
             bcv2(bcvi)=nod
           end if
           if (applies(3)) then
             bcwi=bcwi+1  
             bcw1(bcwi)=value(3)
             bcw2(bcwi)=nod
           end if
         end if
      end do
    end if
  end do

end subroutine getbcuvw12

subroutine getbct_size(field, bct_size)
! get number of legacy dirichlet boundary conditions for scalar field
type(scalar_field), intent(in):: field
integer, intent(out):: bct_size

  integer, dimension(:), pointer:: surface_node_list
  character(len=FIELD_NAME_LEN) type
  integer i
  
  bct_size=0  
  do i=1, get_boundary_condition_count(field)
    call get_boundary_condition(field, i, type=type, &
      surface_node_list=surface_node_list)
    if (type==DIRICHLET_BC_TYPE) then
      bct_size=bct_size+size(surface_node_list)
    end if
  end do

end subroutine getbct_size
  
subroutine getbct(field,  bct1, bct2)
! get legacy dirichlet boundary conditions for scalar field
type(scalar_field), intent(in):: field
real, dimension(:), intent(out):: bct1
integer, dimension(:), intent(out):: bct2

  type(scalar_field), pointer :: surface_field
  integer, dimension(:), pointer :: surface_node_list
  character(len=FIELD_NAME_LEN) type
  integer i, j
  integer bcti
  
  ! index in the array:
  bcti=0
  
  ! loop over the boundary conditions
  do i=1, get_boundary_condition_count(field)
    call get_boundary_condition(field, i, type=type, &
      surface_node_list=surface_node_list)
    if (type==DIRICHLET_BC_TYPE) then
      ! surface field containing b.c. values
      surface_field => extract_surface_field(field, i, DIRICHLET_VALUE_ID)
      ! loop over nodes to which this b.c. applies
      do j=1, size(surface_node_list)
        bct1(bcti+j)=node_val( surface_field, j)
        bct2(bcti+j)=surface_node_list(j)
      end do
      bcti=bcti+size(surface_node_list)
    end if
  end do

end subroutine getbct
  
subroutine getsalphetaire(field, tsndgl, salphe, taire)
! get legacy robin boundary condition
type(scalar_field), intent(in):: field
integer, intent(in):: tsndgl(:)
real, dimension(:), intent(out):: salphe
real, dimension(:), intent(out):: taire

  type(scalar_field), pointer :: surface_field1, surface_field2
  integer, dimension(:), pointer :: surface_node_list, surface_element_list
  integer, dimension(:), pointer :: suf_nodes
  character(len=FIELD_NAME_LEN) type
  integer i, j, sele, siloc, snod, snloc
  
  ! let the fields be zero, where no Robin b.c. is applied
  salphe=0.0
  taire=0.0
  snloc=face_loc(field,1)
  
  ! loop over the boundary conditions
  do i=1, get_boundary_condition_count(field)
    call get_boundary_condition(field, i, type=type, &
      surface_element_list=surface_element_list, &
      surface_node_list=surface_node_list)
    if (type==ROBIN_BC_TYPE) then
      ! surface field containing b.c. values
      surface_field1 => extract_surface_field(field, i, ROBIN_SALPHE_VALUE_ID)
      surface_field2 => extract_surface_field(field, i, ROBIN_TAIRE_VALUE_ID)
      
      ! the following is a bit complicated as legacy uses a numbering
      ! over all surface nodes 1:sufnod (determined by tsndgl)
      ! whereas new b.c.'s use a different numbering for each b.c. within 
      ! the part of the surface the b.c. applies to
      
      ! loop over surface elements to which this b.c. applies
      do j=1, size(surface_element_list)
        ! surface element number in the global face numbering
        sele=surface_element_list(j)
        ! nodes in this element in surface numbering over only the nodes
        ! to which this b.c. applies:
        suf_nodes => ele_nodes(surface_field1, j)
        
        do siloc=1, snloc
          ! surface node in numbering over all surface nodes 1:sufnod
          snod=tsndgl( (sele-1)*snloc+siloc )
          ! suf_nodes(j) is the same nod in the 'new' numbering
          salphe(snod)=node_val(surface_field1, suf_nodes(j))
          taire(snod)=node_val(surface_field2, suf_nodes(j))
        end do
          
      end do
    end if    
  end do  

end subroutine getsalphetaire
  
subroutine getnnodro(velocity_field, nnodro)
! get number of boundary nodes with rotated b.c.
type(vector_field), intent(in):: velocity_field
integer, intent(out):: nnodro

  integer, dimension(:), pointer:: surface_node_list
  character(len=FIELD_NAME_LEN) type
  integer i
  
  nnodro=0  
  do i=1, get_boundary_condition_count(velocity_field)
    
    call get_boundary_condition(velocity_field, i, type=type, &
      surface_node_list=surface_node_list)
      
    if (type==DIRICHLET_BC_TYPE) then
      
      if (has_surface_field(velocity_field, i, NORMAL_VECTOR_ID)) then       
        nnodro=nnodro+size(surface_node_list)
      end if
      
    end if
    
  end do

end subroutine getnnodro

subroutine getnodrot(velocity_field, nodrot, gettan, &
  nx, ny, nz, &
  t1x, t1y, t1z, &
  t2x, t2y, t2z)
! get boundary nodes with rotated b.c.
! if (.not. gettan) copy over user specified normal/tangent vectors
! otherwise zero them
type(vector_field), intent(in):: velocity_field
integer, dimension(:), intent(out):: nodrot
logical, intent(in):: gettan ! whether we have user-specified
real, dimension(:), intent(out):: nx, ny, nz
real, dimension(:), intent(out):: t1x, t1y, t1z
real, dimension(:), intent(out):: t2x, t2y, t2z

  type(vector_field), pointer:: normal_surface_field, &
      tangent1_surface_field, tangent2_surface_field
  integer, dimension(:), pointer:: surface_node_list
  character(len=FIELD_NAME_LEN) type
  real vector(3)
  integer i, j, nnodro, dim
  
  dim=velocity_field%dim  
  
  nnodro=0  
  do i=1, get_boundary_condition_count(velocity_field)
    
    call get_boundary_condition(velocity_field, i, type=type, &
      surface_node_list=surface_node_list)
      
    if (type==DIRICHLET_BC_TYPE .and. &
      
        has_surface_field(velocity_field, i, NORMAL_VECTOR_ID)) then
        
      nodrot(nnodro+1:nnodro+size(surface_node_list))=surface_node_list
      
      if (.not. gettan) then
        
        normal_surface_field => extract_surface_field(velocity_field, i, &
            NORMAL_VECTOR_ID)
        tangent1_surface_field => extract_surface_field(velocity_field, i, &
            TANGENT1_VECTOR_ID)
        tangent2_surface_field => extract_surface_field(velocity_field, i, &
            TANGENT2_VECTOR_ID)
            
        do j=1, size(surface_node_list)
          
          vector=node_val(normal_surface_field, j)
          nx(nnodro+j)=vector(1)
          ny(nnodro+j)=vector(2)
          if (dim==3) nz(nnodro+j)=vector(3)
          
          vector=node_val(tangent1_surface_field, j)
          t1x(nnodro+j)=vector(1)
          t1y(nnodro+j)=vector(2)
          if (dim==3) then
            t1z(nnodro+j)=vector(3)
            
            vector=node_val(tangent2_surface_field, j)
            t2x(nnodro+j)=vector(1)
            t2y(nnodro+j)=vector(2)
            t2z(nnodro+j)=vector(3)
            
          end if
          
        end do

      end if
      
      nnodro=nnodro+size(surface_node_list)
      
    end if
    
    if (gettan) then
      ! zero normal/tangent vectors
      ! first for 2d:
      nx=0; ny=0; t1x=0; t1y=0
      if (dim==3) then
        nz=0; t1z=0
        t2x=0; t2y=0; t2z=0
      end if
    end if
    
  end do

end subroutine getnodrot
  
subroutine copy_back_tangents(velocity_field, &
  nx, ny, nz, &
  t1x, t1y, t1z, &
  t2x, t2y, t2z)
! after tangent vectors have been calculated in navsto
! this routine stores them back in new structure
! so in the next call (after setting gettan to .false.)
! they'll be used again
type(vector_field), intent(inout):: velocity_field
real, dimension(:), intent(in):: nx, ny, nz
real, dimension(:), intent(in):: t1x, t1y, t1z
real, dimension(:), intent(in):: t2x, t2y, t2z

  type(vector_field), pointer:: normal_surface_field, &
      tangent1_surface_field, tangent2_surface_field
  integer, dimension(:), pointer:: surface_node_list
  character(len=FIELD_NAME_LEN) type
  integer i, j, nnodro, dim
  
  dim=velocity_field%dim  
  
  nnodro=0  
  do i=1, get_boundary_condition_count(velocity_field)
    
    call get_boundary_condition(velocity_field, i, type=type, &
      surface_node_list=surface_node_list)
      
    if (type==DIRICHLET_BC_TYPE .and. &
      
        has_surface_field(velocity_field, i, NORMAL_VECTOR_ID)) then
        
        normal_surface_field => extract_surface_field(velocity_field, i, &
            NORMAL_VECTOR_ID)
        tangent1_surface_field => extract_surface_field(velocity_field, i, &
            TANGENT1_VECTOR_ID)
        tangent2_surface_field => extract_surface_field(velocity_field, i, &
            TANGENT2_VECTOR_ID)
            
        do j=1, size(surface_node_list)
          
          if (dim==3) then
            call set( normal_surface_field, j, (/ nx(nnodro+j), &
               ny(nnodro+j), nz(nnodro+j) /))
            call set( tangent1_surface_field, j, (/ t1x(nnodro+j), &
               t1y(nnodro+j), t1z(nnodro+j) /))
            call set( tangent2_surface_field, j, (/ t2x(nnodro+j), &
               t2y(nnodro+j), t2z(nnodro+j) /))
          else
            call set( normal_surface_field, j, (/ nx(nnodro+j), &
               ny(nnodro+j) /))
            call set( tangent1_surface_field, j, (/ t1x(nnodro+j), &
               t1y(nnodro+j) /))
          end if
          
        end do
          
      nnodro=nnodro+size(surface_node_list)
      
    end if
    
  end do

end subroutine copy_back_tangents

subroutine check_velocity_bcs()
!!< in the new frame-work it is possible to specify any combination
!!< of non-rotated/rotated/user-specified rotations b.c.s for any of
!!< the velocity fields of the MaterialPhases
!!<
!!< in legacy fluidity code only one combination of ROTAT and GETTAN
!!< can be used for all b.c.s of all velocities of the different phases
!!<
!!< this routine checks the b.c.s adhere to this restriction and tries
!!< to give some useful error message
!!<
!!< it only allows dirichlet b.c.s (see DIRICHLET_BC_TYPE above)
!!< as others aren't implemented (feel free to do so)
   
  character(len=OPTION_PATH_LEN) path, bc_path
  character(len=FIELD_NAME_LEN) type
  ! everything depends on the first b.c. encountered
  logical first
  logical rotat, gettan
  integer i, j, dim
  
  call get_option("/geometry/dimension", dim)

  first=.true.
  phase_loop: do i=1, option_count("/material_phase")
    path="/material_phase["//int2str(i-1)//"]/vector_field::Velocity"

    ! phases without velocity?! let's ignore that:
    if (.not. have_option(path)) cycle

    path=trim(path)//"/prognostic"
    ! prescribed/diagnostic velocity fields don't have b.c.'s
    ! aliased velocity fields are checked already
    if (.not. have_option(path)) cycle

    path=trim(path)//"/boundary_conditions"

    ! no b.c.s, nothing to check
    if (.not. have_option(path)) cycle

    do j=1, option_count(path)
      bc_path=path//"["//int2str(j-1)//"]"
      call get_option(bc_path//"/type/name", type)
        
      if (type/=DIRICHLET_BC_TYPE) then
        ewrite(0, *) "Error in boundary condition:", bc_path
        FLExit("Only (strong) dirichlet boundary conditions are supported for velocity")
      end if
      
      if (first) then
        ! first valid b.c. determines what the rest should be like:
        rotat=have_option(bc_path//"/type/align_bc_with_surface")
        if (rotat) then
          gettan=.not. have_option(bc_path//"/type/align_bc_with_surface/&
                          &normal_direction")
        else
          gettan=.false.
        end if
        first=.false.
      else
        ! following b.c.'s should be the same
        if (rotat .neqv. have_option(bc_path//"/type/align_bc_with_surface")) then
          ewrite(0, *) "Error in boundary condition:", bc_path
          FLExit("Combination of surface and cartesian aligned boundary conditions for velocity are not supported.")
        end if
        
        if (rotat .and. (gettan .eqv. have_option(bc_path//"/type/&
                          &align_bc_with_surface/rotation_matrix"))) then
          ewrite(0, *) "Error in boundary condition:", bc_path
          FLExit("Either all or none of the velocity boundary conditions should have a user specified rotation matrix.")
        end if
          
      end if
      
      ! check whether either all or none of the normal/tangent directions
      ! are specified:
      if (.not. gettan .and. &
            .not. (have_option(bc_path//"/type/align_bc_with_surface/&
                             &tangent_direction_1") &
                     .and. &
                     (have_option(bc_path//"/type/align_bc_with_surface/&
                             &tangent_direction_2") .eqv. dim==3))) then
         ewrite(0,*) "Specifying only one of the normal or tangent vectors is not yet possible."
         FLExit("Either all or none of the normal and tangent vectors should be specified")
      end if                          
       
    end do
      
  end do phase_loop
        
end subroutine check_velocity_bcs

subroutine set_rotat_gettan(rotat, gettan)
   logical, intent(out):: rotat, gettan
   
   character(len=OPTION_PATH_LEN) velbc_path
   integer i, j

   rotat = .false.
   gettan = .true.
   
   phase_loop: do i=1, option_count("/material_phase")
      velbc_path="/material_phase["//int2str(i-1)//"]/vector_field::Velocity&
         &/prognostic/boundary_conditions"
         
      do j=0, option_count(velbc_path)-1
         if (have_option(trim(velbc_path)//"["//int2str(j)// &
             & "]/type/align_bc_with_surface")) then
           rotat=.true.
           ! gettan is .true. if no user specified rotation matrix:
           gettan=.not. have_option(trim(velbc_path)//"["//int2str(j)// &
             & "]/type/align_bc_with_surface/normal_direction")
           return
         end if
      end do
         
   end do phase_loop
        
end subroutine set_rotat_gettan
   
subroutine set_legacy_boundary_values(state, &
   nobct, nobcu, nobcv, nobcw, &
   bct1, bcu1, bcv1, bcw1, ntsol,nphase)
!!< Copies the boundary condition values from new boundary condition
!!< administration into old RMEM style bc[t|u|v|w]1.
!!< This routine is used after the values have been reset each time step
!!< for time varying boundary conditions.
type(state_type), dimension(:), intent(in):: state
integer, dimension(:), intent(in):: nobct, nobcu, nobcv, nobcw
type(real_vector), dimension(:), intent(inout) :: bct1,bcu1,bcv1,bcw1
integer, intent(in) :: ntsol, nphase

  type(vector_field), pointer:: velocity
  type(scalar_field), pointer:: sfield
  integer, dimension(:), allocatable :: bcu2, bcv2, bcw2, bct2
  integer ip, it, istate, nophas
  
  ! places to copy boundary node numbers in
  ! as we already know this, these are thrown away
  allocate( bcu2(1:maxval(nobcu)), bcv2(1:maxval(nobcv)), &
     bcw2(1:maxval(nobcw)), bct2(1:maxval(nobct)))
  
  nophas=size(bcu1)
  
  ! velocity bcs:
  do ip=1, nophas
     istate= phase2state_index(ip)
     velocity => extract_vector_field(state(istate),'Velocity')
     
     call getbcuvw12(Velocity, &
                      bcu1(ip)%ptr, &
                      bcv1(ip)%ptr, &
                      bcw1(ip)%ptr, &
                      bcu2(1:nobcu(ip)), &
                      bcv2(1:nobcv(ip)), &
                      bcw2(1:nobcw(ip)) )     
  end do
     
  ! scalar field bcs:
  do it=1, ntsol
     sfield => extract_scalar_field( &
                      state(field_state_list(it)), &
                      trim(field_name_list(it)) )
     call getbct(sfield, bct1(it)%ptr, bct2)
  end do
     
  deallocate( bcu2, bcv2, bcw2, bct2 )
   
end subroutine set_legacy_boundary_values
  
end module legacy_boundary_conditions
