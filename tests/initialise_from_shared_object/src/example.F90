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


subroutine val_vec_Fortran90(field,X,t)
  double precision, dimension(:,:) :: field
  double precision, dimension(:,:) :: X
  double precision :: t

  field(1,:)=X(1,:)**2-X(2,:)
  field(2,:)=X(2,:)**2-X(1,:)

end subroutine val_vec_Fortran90


subroutine val_vec_Fortran77(dim1,n,dim2,field,X,t)
  integer :: dim1,dim2, n
  double precision, dimension(dim2,n) :: field
  double precision, dimension(dim1,n) :: X
  double precision :: t

  field(1,:)=X(1,:)**2-X(2,:)
  field(2,:)=X(2,:)**2-X(1,:)

end subroutine val_vec_Fortran77

subroutine val_vec_FortranBindC(dim1,n,dim2,field,X,t) bind(C)
  use iso_c_binding
  integer(c_int) :: dim1,dim2, n
  real(c_double), dimension(dim2,n) :: field
  real(c_double), dimension(dim1,n) :: X
  real(c_double) :: t

  field(1,:)=X(1,:)**2-X(2,:)
  field(2,:)=X(2,:)**2-X(1,:)

end subroutine val_vec_FortranBindC
