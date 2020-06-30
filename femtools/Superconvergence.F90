#include "fdebug.h"

module superconvergence
!!< This module contains the information on the superconvergent
!!< points for various element types.
!!< See 
!!< Zienkiewicz & Zhu, Int. J. Numer. Methods Eng, 33, 1331-1364 (1992)
!! This is primarily used by field_derivatives.

use fldebug
use vector_tools
use elements

implicit none

integer, parameter, public :: MATRIX_SIZE_SPR=4
integer, parameter, public :: MATRIX_SIZE_QF=11
integer, parameter, public :: MATRIX_SIZE_CF=20
integer, parameter, public :: MATRIX_SIZE_CF_2D=10
integer, parameter, public :: MATRIX_SIZE_QF_2D=6
integer, dimension(MATRIX_SIZE_QF_2D), parameter, public ::&
     QF_2D_X = (/1, 3, 4, 6, 7, 10/)
integer, dimension(MATRIX_SIZE_QF_2D), parameter, public ::&
     QF_2D_Y = (/1, 2, 4, 5, 7, 9/)
integer, dimension(MATRIX_SIZE_QF_2D), parameter, public ::&
     QF_2D_Z = (/1, 2, 3, 5, 6, 8/)
integer, dimension(MATRIX_SIZE_CF_2D), parameter, public ::&
     CF_2D_X = (/1, 3, 4, 7, 9, 10, 15, 17, 19, 20/)
integer, dimension(MATRIX_SIZE_CF_2D), parameter, public ::&
     CF_2D_Y = (/1, 2, 4, 6, 8, 10, 13, 16, 18, 20/)
integer, dimension(MATRIX_SIZE_CF_2D), parameter, public ::&
     CF_2D_Z = (/1, 2, 3, 5, 8, 9,  12, 14, 18, 19/)
type(superconvergence_type), save, target :: superconvergence_tet_array(1)

private

public :: initialise_superconvergence, get_superconvergence, getP_spr,&
     compute_matrix_contribution_cf, getP_qf, compute_rhs_contribution_qf,&
     compute_rhs_contribution_spr, compute_rhs_contribution_cf,&
     compute_matrix_contribution_qf, compute_matrix_contribution_spr,&
     evaluate_cf, evaluate_qf


contains

    subroutine initialise_superconvergence
        logical, save :: initialised = .false.

        if (initialised) return
        initialised = .true.

        superconvergence_tet_array(1)%nsp = 6
        allocate(superconvergence_tet_array(1)%l(6, 4))
        allocate(superconvergence_tet_array(1)%n(4, 6))
        allocate(superconvergence_tet_array(1)%dn(4, 6, 3))

        ! The midpoints of the edges.
        superconvergence_tet_array(1)%l(1, :) = (/0.50, 0.50, 0.00, 0.00/)
        superconvergence_tet_array(1)%l(2, :) = (/0.50, 0.00, 0.50, 0.00/)
        superconvergence_tet_array(1)%l(3, :) = (/0.50, 0.00, 0.00, 0.50/)
        superconvergence_tet_array(1)%l(4, :) = (/0.00, 0.50, 0.50, 0.00/)
        superconvergence_tet_array(1)%l(5, :) = (/0.00, 0.50, 0.00, 0.50/)
        superconvergence_tet_array(1)%l(6, :) = (/0.00, 0.00, 0.50, 0.50/)
    end subroutine initialise_superconvergence

    function get_superconvergence(element) result(superconvergence)
        type(superconvergence_type), pointer :: superconvergence
        type(element_type), intent(in) :: element

        integer :: i, j

        call initialise_superconvergence

        if (element%dim == 3 .and. element%degree == 1 .and. element%loc == 4) then
            superconvergence => superconvergence_tet_array(1)

            do j=1,superconvergence%nsp
              superconvergence%n(:, j)     = eval_shape(element,  superconvergence%l(j, :))
            end do
            do i=1,element%loc
                do j=1,superconvergence%nsp
                    superconvergence%dn(i, j, :) = eval_dshape(element, i, superconvergence%l(j, :))
                end do
            end do

        else
            superconvergence => null()
        end if

    end function get_superconvergence

!------------------------------------
!     pure function matrix_size_spr()
!       integer :: matrix_size_spr
! 
!       matrix_size_spr = 4
!     end function matrix_size_spr
!------------------------------------

    function getP_spr(positions, element) result(P)
      !!< This computes P as defined in the SPR paper. See the ref above.
      real, intent(in) :: positions(:)
      type(element_type) :: element
      real, dimension(MATRIX_SIZE_SPR) :: P

      if (element%dim == 3 .and. element%degree == 1 .and. element%loc == 4) then
        P(1) = 1.0
        P(2:4) = positions(1:3)
      else
        FLAbort("P not coded for this element type!")
      end if
    end function

    function compute_matrix_contribution_spr(positions, element) result(pTp)
      !!< See the superconvergent patch recovery paper (reference is available above)
      !!< for this. P is a vector function of position, depending on the element.
      !!< For a given node, the algorithm solves Ax = b, where
      !!< A = sum(over superconvergent points) of P^T(x, y, z) * P(x, y, z)
      !!< This function returns P^T(x, y, z) * P(x, y, z) for a given position.
   
      type(element_type), intent(in) :: element
      real :: positions(:)
      real, dimension(MATRIX_SIZE_SPR, MATRIX_SIZE_SPR) :: pTp
      real, dimension(MATRIX_SIZE_SPR) :: P

      assert(element%dim .eq. size(positions))

      P = getP_spr(positions, element)
      pTp= outer_product(P, P)

    end function compute_matrix_contribution_spr

    function compute_rhs_contribution_spr(positions, element, derivative) result(b)
    !!< This function is the same as above, but for the right hand side.
    !!< Unlike the matrix, the rhs depends on which derivative you're taking.
    !!< Here b = sum(over superconvergent points) of P^T(x, y, z) * diff(field, coordinate)(x, y, z)
    
      type(element_type), intent(in) :: element
      real :: positions(:)
      real, intent(in) ::  derivative
      real, dimension(MATRIX_SIZE_SPR) :: b

      assert(element%dim .eq. size(positions))

      b = derivative * getP_spr(positions, element) 
    end function compute_rhs_contribution_spr

!-----------------------------------
!     pure function matrix_size_qf()
!       integer :: matrix_size_qf
! 
!       matrix_size_qf = 11
!     end function matrix_size_qf
!-----------------------------------

    function getP_qf(positions) result(P)
      !!< This computes P as defined in the QF paper. See the ref above.
      real, intent(in) :: positions(:)
      real, dimension(MATRIX_SIZE_QF) :: P
      real :: x, y, z

      x = positions(1); y = positions(2); z = positions(3)
      P(1) = 1.0
      P(2:4) = positions(1:3)
      P(5) = x**2
      P(6) = y**2
      P(7) = z**2
      P(8) = x * y
      P(9) = x * z
      P(10) = y * z
      P(11) = x * y * z
    end function

    function compute_matrix_contribution_qf(positions) result(pTp)
      !!< P is a vector function of position, representing a quadratic fit.
      !!< For a given node, the algorithm solves Ax = b, where
      !!< A = sum(over superconvergent points) of P^T(x, y, z) * P(x, y, z)
      !!< This function returns P^T(x, y, z) * P(x, y, z) for a given position.
   
      real :: positions(:)
      real, dimension(MATRIX_SIZE_QF, MATRIX_SIZE_QF) :: pTp
      real, dimension(MATRIX_SIZE_QF) :: P

      P = getP_qf(positions)
      pTp = outer_product(P, P)
    end function compute_matrix_contribution_qf

    function compute_rhs_contribution_qf(positions, derivative) result(b)
    !!< This function is the same as above, but for the right hand side.
    !!< Here b = sum(over patch points) of P^T(x, y, z) * diff(field, coordinate)(x, y, z)
    
      real :: positions(:)
      real, intent(in) ::  derivative
      real, dimension(MATRIX_SIZE_QF) :: b

      b = derivative * getP_qf(positions) 
    end function compute_rhs_contribution_qf

    function evaluate_qf(b, positions) result(fitted)
    !!< Evaluate the quadratic fit at the position given.
      real, dimension(:), intent(in) :: b, positions
      real :: fitted, x, y, z

      x = positions(1); y = positions(2); z = positions(3)
      fitted = b(1) + x * b(2) + y * b(3) + z * b(4) + &
               x**2 * b(5) + y**2 * b(6) + z**2 * b(7) + &
               x * y * b(8) + x * z * b(9) + y * z * b(10) + &
               x * y * z * b(11)
    end function evaluate_qf

    function getP_cf(positions) result(P)
      !!< This computes P as defined in the QF paper. See the ref above.
      real, intent(in) :: positions(:)
      real, dimension(MATRIX_SIZE_CF) :: P
      real :: x, y, z

      x = positions(1); y = positions(2); z = positions(3)
      P(1) = 1.0
      P(2:4) = positions(1:3)
      P(5) = x * y
      P(6) = x * z
      P(7) = y * z
      P(8) = x**2
      P(9) = y**2
      P(10) = z**2
      P(11) = x * y * z
      P(12) = x**2 * y
      P(13) = x**2 * z
      P(14) = y**2 * x
      P(15) = y**2 * z
      P(16) = z**2 * x
      P(17) = z**2 * y
      P(18) = x**3
      P(19) = y**3
      P(20) = z**3
    end function

    function compute_matrix_contribution_cf(positions) result(pTp)
      !!< P is a vector function of position, representing a cubic fit.
      !!< For a given node, the algorithm solves Ax = b, where
      !!< A = sum(over superconvergent points) of P^T(x, y, z) * P(x, y, z)
      !!< This function returns P^T(x, y, z) * P(x, y, z) for a given position.
   
      real :: positions(:)
      real, dimension(MATRIX_SIZE_CF, MATRIX_SIZE_CF) :: pTp
      real, dimension(MATRIX_SIZE_CF) :: P

      P = getP_cf(positions)
      pTp = outer_product(P, P)
    end function compute_matrix_contribution_cf

    function compute_rhs_contribution_cf(positions, derivative) result(b)
    !!< This function is the same as above, but for the right hand side.
    !!< Here b = sum(over patch points) of P^T(x, y, z) * diff(field, coordinate)(x, y, z)
    
      real :: positions(:)
      real, intent(in) ::  derivative
      real, dimension(MATRIX_SIZE_CF) :: b

      b = derivative * getP_cf(positions) 
    end function compute_rhs_contribution_cf

    function evaluate_cf(b, positions) result(fitted)
    !!< Evaluate the cubic fit at the position given.
      real, dimension(:), intent(in) :: b, positions
      real :: fitted, x, y, z

      x = positions(1); y = positions(2); z = positions(3)
      fitted = b(1) + x * b(2) + y * b(3) + z * b(4) + &
               x * y * b(5) + x * z * b(6) + y * z * b(7) + &
               x**2 * b(8) + y**2 * b(9) + z**2 * b(10) + &
               x * y * z * b(11) + &
               x**2 * y * b(12) + x**2 * z * b(13) + y**2 * x * b(14) + &
               y**2 * z * b(15) + z**2 * x * b(16) + z**2 * y * b(17) + &
               x**3 * b(18) + y**3 * b(19) + z**3 * b(20)
    end function evaluate_cf
end module superconvergence
