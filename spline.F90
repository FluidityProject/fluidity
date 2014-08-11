Program spline
IMPLICIT NONE

integer :: i,j,k,p, i2
double precision :: spline_coefficient(4,12), e2(22), ek2(22)
integer, parameter :: n=12
double precision :: e(n), dedt_k(n)
double precision, dimension(12) :: a,b,c,d

interface
  subroutine cubic_spline_coefficient(e,dedt_k,spline_coefficient)
        double precision, intent(in) :: e(:), dedt_k(:)
        double precision, intent(out) :: spline_coefficient(4,size(e))
  end subroutine cubic_spline_coefficient
end interface

e=(/0.1544535825,0.2748526645,0.3625068659,0.412745937,0.4674295763,0.5282583069,&
      0.5772440306,0.6187169157,0.7164262117,0.7931250187,0.8385892222,0.8957380133/)

dedt_k=(/11115.5978729463,3824.9074006973,2248.5407453053,1843.3982021721,1387.1490141697,&
1083.730080291,898.7069155121,744.1679956897,534.9701226336,455.38544496,371.4353697503,329.4397601335/)

call cubic_spline_coefficient(e,dedt_k,spline_coefficient)

a=spline_coefficient(1,:)
b=spline_coefficient(2,:)
c=spline_coefficient(3,:)
d=spline_coefficient(4,:)

e2=(/(i2*0.01, i2=0,105,5)/)

do k=1,size(e2)
  
  if (e2(k)<e(1)) then
    ek2(k)=b(1)*(e2(k)-e(1))+a(1)
  else if (e2(k)>e(n)) then
    ek2(k)=b(n)*(e2(k)-e(n))+a(n)
  else
    i=0
    j=n
    do while (j-i>1)
      if (mod(j-i,2)==0) then
         p=(j-i)/2
      else
         p=(j-i+1)/2
      end if
       
      if (e2(k)>e(i+p)) then
         i=i+p
      else
         j=i+p        
      end if
    end do
   
    ek2(k)=a(i)+b(i)*(e2(k)-e(i))+c(i)*((e2(k)-e(i))**2.0)+d(i)*((e2(k)-e(i))**3.0)
  end if
end do

open(unit=1, file='spline.dat')
do i=1,size(e2)
  write(1,*) e2(i), ek2(i)
end do
close(unit=1)

end program

subroutine cubic_spline_coefficient(e,dedt_k,spline_coefficient)

 !type(darcy_impes_type), intent(inout) :: di
 integer :: i
 double precision, intent(in) :: e(:), dedt_k(:) !the extraction and dextraction_over_k from experiment data
 !the parameters for spline interpolation
 integer :: n,n1,n2
 double precision, intent(out) :: spline_coefficient(4,size(e))!cubic spline coefficients and the last node of data

 !parameters pass to thomas algrithm
 double precision :: dh(size(e)-1), dy(size(e)-1) !the interval between each two nodes
 double precision :: c(size(e))
 
 n=size(e)
 n1=n-1
 n2=n-2

 do i=1, n-1
   dh(i)=e(i+1)-e(i)
   dy(i)=dedt_k(i+1)-dedt_k(i)
 end do
 

 call solve_tridiagonal_matrix(n2,dh,dy,c) 

 do i=1, n1
   !the zero order term coefficient
   spline_coefficient(1,i)=dedt_k(i)
   !the 1st order term coefficient
   spline_coefficient(2,i)=dy(i)/dh(i)-c(i)*(dh(i)/2.0)-(dh(i)/6.0)*(c(i+1)-c(i))
   !the 2nd order term coefficient
   spline_coefficient(3,i)=c(i)/2.0
   !the 3rd order term coefficient
   spline_coefficient(4,i)=(c(i+1)-c(i))/(6.0*dh(i))
 end do
 
 spline_coefficient(1,n)=dedt_k(n)
 spline_coefficient(2,n)=(dedt_k(n)-dedt_k(n-1))/(e(n)-e(n-1))
 spline_coefficient(3:4,n)=0.0

 !allocate(di%lc%spline_coe(4,n1)
 !di%lc%spline_coe=spline_coefficient
 
contains

subroutine solve_tridiagonal_matrix(n2,dh,dy,c)
     !solve tridiagonal matrix by Thomas Algrithm
     integer, intent(in) :: n2
     double precision, intent(in) :: dh(:), dy(:)
     double precision, intent(out) :: c(n2+2)
     
     integer :: j
     double precision :: x(n2)
     double precision :: coefficient(4,n2) !coefficient for matrix and vector
                               !a_(i)*x_(i-1)+b_(i)*x_(i)+c_(i)*x(i+1)=d_(i)

     do j=1,n2
        !coefficient a,b,c,d
        coefficient(1,j)=dh(j)
        coefficient(2,j)=2.0*(dh(j)+dh(j+1))
        coefficient(3,j)=dh(j+1)
        coefficient(4,j)=6.0*(dy(j+1)/dh(j+1)-dy(j)/dh(j))
     end do

     do j=1,n2
       if (j == 1) then
         coefficient(3,j)=coefficient(3,j)/coefficient(2,j)
         coefficient(4,j)=coefficient(4,j)/coefficient(2,j)
       else
         coefficient(3,j)=coefficient(3,j)/(coefficient(2,j)-coefficient(3,j-1)*coefficient(1,j))
         coefficient(4,j)=(coefficient(4,j)-coefficient(4,j-1)*coefficient(1,j))/&
                           (coefficient(2,j)-coefficient(3,j-1)*coefficient(1,j))
       end if
     end do

     do j=n2,1,-1
      if (j == n2) then
        x(j)=coefficient(4,j)
      else 
        x(j)=(coefficient(4,j)-coefficient(3,j)*x(j+1))
      end if
     end do
    
     c=[0.0D0,x,0.0D0]
     
    
end subroutine solve_tridiagonal_matrix


end subroutine cubic_spline_coefficient
