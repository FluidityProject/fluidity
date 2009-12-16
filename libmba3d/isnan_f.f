      subroutine mba3d_fpcheck(a, flag)
        real :: a
        logical :: flag
        integer :: flag_int

        call fpcheck_mba3d(a, flag_int)
        flag = (flag_int == 1)
      end subroutine mba3d_fpcheck
