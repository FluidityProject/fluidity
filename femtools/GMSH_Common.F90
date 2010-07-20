module gmsh_common

  character(len=3), parameter :: GMSHVersionStr = "2.1"
  integer, parameter :: asciiFormat = 0
  integer, parameter :: binaryFormat = 1
  ! Anyway to automatically calc this in Fortran?
  integer, parameter :: doubleNumBytes = 8

  integer, parameter :: longStringLen = 1000
  real, parameter :: verySmall = 10e-10

  ! For each type, the number of nodes. -1 means unsupported
  integer, dimension(15) :: elementNumNodes = (/ &
       2, 3, 4, 4, 8, &
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 /)

  type GMSHnode
     integer :: nodeID
     double precision :: x(3)
! Currently unused
! real, pointer :: properties(:)
  end type GMSHnode

  type GMSHelement
     integer elementID, type, numTags, physicalID, elementary
     integer, pointer :: tags(:), nodeIDs(:)
  end type GMSHelement


contains

  ! -----------------------------------------------------------------
  ! Change already-open file to ASCII formatting
  ! Involves a bit of sneaky code.
  
  subroutine ascii_formatting(fd, filename, readWriteStr)
    integer fd
    character(len=*) :: filename, readWriteStr

    integer position
    
    
    inquire(fd, POS=position)
    close(fd)
    
    select case( trim(readWriteStr) )
    case("read")
       open( fd, file=trim(filename), action="read", form="formatted", &
            access="stream")
       read( fd, "(I1)", POS=position, ADVANCE="no" )

    case("write")
       open( fd, file=trim(filename), action="write", form="formatted", &
            access="stream")
    end select


  end subroutine ascii_formatting


  ! -----------------------------------------------------------------
  ! Change already-open file to binary formatting
  ! Sneaky code, as above.
  
  subroutine binary_formatting(fd, filename, readWriteStr)
    integer fd
    character(len=*) filename, readWriteStr

    integer position

    
    inquire(fd, POS=position)
    close(fd)

    select case( trim(readWriteStr) )
    case("read")
       open( fd, file=trim(filename), action="read", form="unformatted", &
            access="stream")
       read( fd, POS=position )

    case("write")
       open( fd, file=trim(filename), action="write", form="unformatted", &
            access="stream")
    end select

  end subroutine binary_formatting


  ! -----------------------------------------------------------------
  ! Reorder to Fluidity node ordering
  
  subroutine toFluidityElementNodeOrdering( oldList )
    integer, pointer :: oldList(:), flNodeList(:), nodeOrder(:)
    integer i
    
    numNodes = size(oldList)
    allocate( flNodeList(numNodes) )
    allocate( nodeOrder(numNodes) )
    
    ! Specify node ordering
    select case( numNodes )
    case (3)
       nodeOrder = (/1, 2, 3/)
    case (6)
       nodeOrder = (/1, 6, 2, 5, 4, 3/)
    case default
       do i=1, numNodes
          nodeOrder(i) = i
       end do
    end select
    
    ! Reorder nodes
    do i=1, numNodes
       flNodeList(i) = oldList( nodeOrder(i) )
    end do

    ! Allocate to original list, and dealloc temp list.
    oldList(:) = flNodeList(:)
    deallocate( flNodeList )
    !deallocate(nodeOrder)
    
end subroutine toFluidityElementNodeOrdering


  ! -----------------------------------------------------------------
! Reorder Fluidity node ordering to GMSH

subroutine toGMSHElementNodeOrdering( oldList )
  integer, pointer :: oldList(:), gmshNodeList(:), nodeOrder(:)
  integer i
  

  numNodes = size(oldList)
  allocate( gmshNodeList(numNodes) )
  allocate( nodeOrder(numNodes) )
  
  ! Specify node ordering
  select case( numNodes )
  case (3)
     nodeOrder = (/1, 2, 3/)
  case (6)
     nodeOrder = (/1, 3, 6, 5, 4, 2/)
  case default
     do i=1, numNodes
        nodeOrder(i) = i
     end do
  end select

  ! Reorder nodes
  do i=1, numNodes
     gmshNodeList(i) = oldList( nodeOrder(i) )
  end do
  
  ! Allocate to original list, and dealloc temp list.
  oldList(:) = gmshNodeList(:)
  
  deallocate( gmshNodeList )
  !deallocate( nodeOrder )
  
end subroutine toGMSHElementNodeOrdering



subroutine deallocateElementList( elements )
  type(GMSHelement), pointer :: elements(:)
  integer i
  
  do i = 1, size(elements)
     deallocate(elements(i)%tags)
     deallocate(elements(i)%nodeIDs)
  end do
  
  deallocate( elements )
  
end subroutine deallocateElementList

end module gmsh_common
