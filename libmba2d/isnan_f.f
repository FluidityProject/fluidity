      subroutine fpcheck(a, flag)
        real :: a
        logical :: flag
        integer :: flag_int

        call fpcheck_c(a, flag_int)
        flag = (flag_int == 1)
      end subroutine fpcheck
