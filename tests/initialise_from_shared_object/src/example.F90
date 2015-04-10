subroutine val(field,X,t)
  double precision, dimension(:) :: field
  double precision, dimension(:,:) :: X
  double precision :: t

  field=X(1,:)**2-X(2,:)

end subroutine val
