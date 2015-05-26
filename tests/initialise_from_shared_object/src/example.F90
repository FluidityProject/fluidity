subroutine val_Fortran90(field,X,t)
  double precision, dimension(:) :: field
  double precision, dimension(:,:) :: X
  double precision :: t

  field=X(1,:)**2-X(2,:)

end subroutine val_Fortran90


subroutine val_Fortran77(dim,n,field,X,t)
  integer :: dim, n
  double precision, dimension(n) :: field
  double precision, dimension(dim,n) :: X
  double precision :: t

  field=X(1,:)**2-X(2,:)

end subroutine val_Fortran77

subroutine val_FortranBindC(dim,n,field,X,t) bind(C)
  use iso_c_binding
  integer(c_int) :: dim, n
  real(c_double), dimension(n) :: field
  real(c_double), dimension(dim,n) :: X
  real(c_double) :: t

  field=X(1,:)**2-X(2,:)
end subroutine val_FortranBindC
