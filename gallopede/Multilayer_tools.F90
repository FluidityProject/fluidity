module multilayer_tools

  implicit none

  type layer
     real, dimension(:), pointer :: l
  end type layer

  interface clear
     module procedure clear_layer
  end interface

  contains

    subroutine clear_layer(D)
      type(layer), dimension(:), intent(out) :: D
      
      !locals
      integer :: i
      
      do i = 1, size(D)
         D(i)%l = 0.
      end do
      
    end subroutine clear_layer

end module multilayer_tools
