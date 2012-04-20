module m_multa
contains
subroutine multa(A, X, ndim, nrens, iblkmax)
implicit none
integer, intent(in) :: ndim
integer, intent(in) :: nrens
integer, intent(in) :: iblkmax
real, intent(in)    :: X(nrens,nrens)
real, intent(inout) :: A(ndim,nrens)
real v(iblkmax,nrens)  ! Automatic work array

integer ifirst,ilast
do ifirst = 1,ndim,iblkmax
  ilast = min(ifirst+iblkmax-1,ndim)
  v(1:ilast-ifirst+1,1:nrens) = A(ifirst:ilast,1:nrens)
  call dgemm('n','n', ilast-ifirst+1, nrens, nrens, &
              1.0, v(1,1), iblkmax, &
              X(1,1), nrens, &
              0.0, A(ifirst,1), ndim)
enddo
end subroutine multa
end module m_multa
