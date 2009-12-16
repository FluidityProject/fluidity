#include "fdebug.h"
module nc_tools

  use sparse_tools
  use elements
  use transform_elements
  use dgtools
  use data_structures
  use text_io

  implicit none

contains

  function get_bc_locnods_nc(ele,globi,evlist_u) result(locnods)
    implicit none
    integer, intent(in) :: ele,globi
    integer, intent(in), dimension(3) :: evlist_u
    integer, dimension(3) :: locnods

   !locals
    integer :: loc
    integer :: offnod

    offnod = 0
    
    do loc = 1,3
       if(evlist_u(loc)==globi) offnod = loc
    end do
    
    select case (offnod)
    case (1)
       locnods(1) = 2
       locnods(2) = 3
    case (2)
       locnods(1) = 1
       locnods(2) = 3
    case (3)
       locnods(1) = 1
       locnods(2) = 2
    case default
       write(0,*) 'BC error'
       stop
    end select   

    locnods(3) = offnod

  end function get_bc_locnods_nc

  function nc_surface_integral_grad(n,n_f,detwei_f, &
       normal,dnu_t,b_seg) result (bcmatnc)
    type(element_type), intent(in) :: n_f, n
    real, dimension(n_f%ngi), intent(in) :: detwei_f
    real, dimension(2,n_f%ngi), intent(in) :: normal
    real, dimension(n%loc,n%ngi,2), intent(in) :: dnu_t
    integer, dimension(2), intent(in) :: b_seg
    real, dimension(3,3,2,2) :: bcmatnc

    !locals 
    real, dimension(2,3,2,2) :: bcmatcg
    real, dimension(3,3) :: cg2ncmat
    integer :: iloc, jloc ,i, j
    
    cg2ncmat(:,1) = (/-1.0, 1.0,  1.0/)
    cg2ncmat(:,2) = (/1.0, -1.0,  1.0/)
    cg2ncmat(:,3) = (/1.0,  1.0, -1.0/)

    bcmatcg = 0.
    bcmatnc = 0.
    !contribution from inside the element
    do iloc = 1,2
       do jloc = 1,3 
          do i = 1,2
             do j = 1,2
                bcmatcg(iloc,jloc,i,j) = sum(n_f%n(iloc,:)* &
                     detwei_f(:)*normal(i,:))*dnu_t(jloc,1,j)
             end do
          end do
       end do
    end do
    
    !need to map from nh_f back to nc
    forall(i=1:2, j=1:2)
       bcmatnc(:,:,i,j) = bcmatnc(:,:,i,j) + &
            matmul(cg2ncmat(:,b_seg), &
            bcmatcg(:,:,i,i))
    end forall

  end function nc_surface_integral_grad

  function nc_surface_integral(n,n_f,detwei_f, &
       normal,b_seg,b_seg_2) result (bcmatnc)
    type(element_type), intent(in) :: n_f,n
    real, dimension(n_f%ngi), intent(in) :: detwei_f
    real, dimension(2,n_f%ngi), intent(in) :: normal
    integer, dimension(2), intent(in) :: b_seg, b_seg_2
    real, dimension(3,3,2) :: bcmatnc

    !locals 
    real, dimension(2,2,2) :: bcmatcg
    real, dimension(3,3) :: cg2ncmat
    real, dimension(3,2) :: cgn2ncmat_f, cgn2ncmat_f_2
    integer :: iloc, jloc, i
   
     ! Nnc_i = C_ij Ncg_j
    cg2ncmat(:,1) = (/-1.0, 1.0,  1.0/)
    cg2ncmat(:,2) = (/1.0, -1.0,  1.0/)
    cg2ncmat(:,3) = (/1.0,  1.0, -1.0/)

    bcmatcg = 0.
    bcmatnc = 0.
    !contribution from inside the element
    do iloc = 1,2
       do jloc = 1,2
          do i = 1,2
             bcmatcg(iloc,jloc,i) = sum(n_f%n(iloc,:)* &
                  detwei_f(:)*normal(i,:)*n_f%n(jloc,:))
          end do
       end do
    end do
    
    !need to map from nh_f back to nc
    forall(i=1:2)
       bcmatnc(:,:,i) = bcmatnc(:,:,i) + &
            matmul(cg2ncmat(:,b_seg), &
            matmul(bcmatcg(:,:,i),transpose(cg2ncmat(:,b_seg_2))))
    end forall

  end function nc_surface_integral

end module nc_tools
