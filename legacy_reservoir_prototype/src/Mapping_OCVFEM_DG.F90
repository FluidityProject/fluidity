!    Copyright (C) 2011 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"




module mapping_for_ocvfem

 use fldebug
 use shape_functions
 use vector_tools
 use shape_functions_Linear_Quadratic




 private
 
 public :: overlapping_to_quadratic_dg
 
 contains


 ! This subroutine projects the overlapping control volume to the finite element shape functions
 !   in particular the discontinuous quadratic pressure function space from the overlapping
 !   velocity function space

 subroutine overlapping_to_quadratic_dg( &
       cv_nonods, x_nonods,u_nonods,  totele, &
       cv_ele_type,  &
       nphase,  &
       cv_nloc, u_nloc, x_nloc, &
       cv_ndgln,  u_ndgln, x_ndgln,&
       cv_snloc, u_snloc, stotel, cv_sndgln, u_sndgln, &
       x, y, z, &
       u, v, w, uold, vold, wold,velocity_dg, ndim )

    use shape_functions
    use matrix_operations
    use printout

    implicit none

    ! inputs/outputs
    
    
    integer, intent( in ) :: cv_nonods, u_nonods, &
         totele, x_nonods,&
         cv_ele_type, x_nloc,&
         nphase, cv_nloc, u_nloc, &
         cv_snloc, u_snloc, stotel,  ndim

    integer, dimension( totele * cv_nloc ), intent( in ) :: cv_ndgln
    integer, dimension( totele * u_nloc ), intent( in ) :: u_ndgln 
    integer, dimension( stotel * cv_snloc ), intent( in ) :: cv_sndgln
    integer, dimension( stotel * u_snloc ), intent( in ) :: u_sndgln 
    integer, dimension( totele * x_nloc ), intent( in ) :: x_ndgln 
    real, dimension( x_nonods ), intent( in ) :: x, y, z
    real, dimension( u_nonods * nphase ), intent( in ) :: u, v, w, uold, vold, wold
    real, dimension(cv_nonods,nphase,ndim), intent(inout) :: velocity_dg

    ! local variables - allocatable arrays
    integer, dimension( : ), allocatable :: findgpts, &
         cv_other_loc, u_other_loc, mat_other_loc, &
         jcount_kloc, jcount_kloc2, colgpts, cv_sloc2loc, u_sloc2loc
    integer, dimension( : , : ), allocatable :: cv_sloclist, u_sloclist, &
         face_ele, cv_neiloc
    real, dimension( : , : ), allocatable :: cvn, cvn_short, cvfen, cvfenlx, cvfenly, cvfenlz, &
         cvfen_short, cvfenlx_short, cvfenly_short, cvfenlz_short,  &
         ufen, ufenlx, ufenly, ufenlz, scvfen, scvfenslx, scvfensly, &
         scvfenlx, scvfenly, scvfenlz, matloc,n, wdold,wd,vdold,vd, ud, udold,&
         sufen, sufenslx, sufensly, sufenlx, sufenly,sufenlz, nlx, nly, nlz, nx,ny,nz
    real, dimension( : , : ), allocatable :: sbcvfen, sbcvfenslx, sbcvfensly, &
         sbcvfenlx, sbcvfenly, sbcvfenlz, sbufen, sbufenslx, sbufensly, &
         sbufenlx, sbufenly, sbufenlz
    real, dimension( : ), allocatable :: cvweight,cvweight_short,scvfeweigh,sbcvfeweigh
    real, dimension( : ), allocatable :: tmax, tmin, toldmax, &
         toldmin, denmax, denmin, denoldmax, denoldmin, cvnormx, &
         cvnormy, cvnormz, mass_cv, mass_ele, sndotq, sndotqold,  &
         femt, femtold, femt2, femt2old, femden, femdenold, xc_cv, yc_cv, zc_cv, &
         scvdetwei, sra, ugi_coef_ele, vgi_coef_ele, wgi_coef_ele, &
         ugi_coef_ele2, vgi_coef_ele2, wgi_coef_ele2,  &
         sum_cv, one_pore, sele_overlap_scale, t2max, t2min, t2oldmax, &
         t2oldmin, rhsloc_u, rhsloc_v, rhsloc_w, rhsloc_uold, rhsloc_vold, rhsloc_wold, solloc
    integer, dimension( : ), allocatable :: tmax_nod, tmin_nod, toldmax_nod, &
         toldmin_nod, denmax_nod, denmin_nod, denoldmax_nod, denoldmin_nod, &
         t2max_nod, t2min_nod, t2oldmax_nod, t2oldmin_nod
    real, dimension( : , :, : ), allocatable :: dtx_ele,dty_ele,dtz_ele,  &
         dtoldx_ele,dtoldy_ele,dtoldz_ele
    real, dimension( : ), allocatable :: up_wind_nod, detwei, ra

    logical, dimension( : ), allocatable :: x_share,log_on_bound
    logical, dimension( :, : ), allocatable :: cv_on_face,u_on_face
    logical, parameter :: include_pore_vol_in_deriv = .false.

    !          ===> integers <====
    integer :: cv_ngi, cv_ngi_short, scvngi, sbcvngi, count, jcount, &
         ele, ele2, gi, gcount, sele,   &
         ncolgpts, p_ele_type,&
         cv_siloc, u_kloc,u_nod_pha, &
         cv_iloc, cv_jloc, iphase, jphase, &
         cv_nodj, cv_nodj_ipha, &
         cv_nodi, cv_nodi_ipha, cv_nodi_jpha, u_nodk, timopt, &
         jcount_ipha, imid_ipha, &
         nface, x_nodi, u_iloc, u_nod, &
         cv_inod, mat_nodi, face_its, nface_its, cv_gi, nloc
    !        ===>  reals  <===
    real :: ndotq, ndotqold, nn, volume, &
         income, incomeold, hdc, fvt, fvtold, fvt2, fvt2old, &
         fvd, fvdold, limt, limtold, limt2, limt2old,&
         limd, limdold, ftheta, vtheta, &
         limdt, limdtold, limdtt2, limdtt2old, &
         femdgi, femtgi,femt2gi, femdoldgi, femtoldgi, femt2oldgi, &
         tmid, toldmid, &
         diff_coef_divdx, diff_coefold_divdx, bczero, robin1, robin2, &
         sum, &
         sum_limt, sum_limtold, ftheta_t2, one_m_ftheta_t2old
    ! functions...
    !real :: r2norm,face_theta  
    !        ===>  logicals  <===
    logical :: getmat, &
         d1, d3, dcyl, got_diffus, integrat_at_gi, &
         normalise, sum2one, get_gtheta


    call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
         cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface )
         
         
    !cv_ngi_short = cv_ngi
    ewrite(2,*) ' got FE info', cv_ngi, cv_ngi_short, cv_nloc, u_nloc
    

    ! allocate memory for the control volume surface shape functions, etc.
    allocate( jcount_kloc(  u_nloc ))
    allocate( jcount_kloc2(  u_nloc ))

    allocate( cvnormx( scvngi ))
    allocate( cvnormy( scvngi ))
    allocate( cvnormz( scvngi ))
    allocate( colgpts( cv_nloc * scvngi )) !the size of this vector is over-estimated
    allocate( findgpts( cv_nloc + 1 ))
    allocate( sndotq( scvngi ))
    allocate( sndotqold( scvngi ))
    allocate( cv_on_face( cv_nloc, scvngi ))
    allocate( u_on_face( u_nloc, scvngi ))
    allocate( cv_other_loc( cv_nloc ))
    allocate( u_other_loc( u_nloc ))
    allocate( x_share( x_nonods ))
    allocate( cvweight( cv_ngi ))
    allocate( cvn( cv_nloc, cv_ngi ))
    allocate( cvfen( cv_nloc, cv_ngi ))
    allocate( cvfenlx( cv_nloc, cv_ngi ))
    allocate( cvfenly( cv_nloc, cv_ngi ))
    allocate( cvfenlz( cv_nloc, cv_ngi ))

    allocate( cvweight_short( cv_ngi_short ))
    allocate( cvn_short( cv_nloc, cv_ngi_short ))
    allocate( cvfen_short( cv_nloc, cv_ngi_short))
    allocate( cvfenlx_short( cv_nloc, cv_ngi_short ))
    allocate( cvfenly_short( cv_nloc, cv_ngi_short ))
    allocate( cvfenlz_short( cv_nloc, cv_ngi_short ))

    allocate( ufen( u_nloc, cv_ngi)) 
    allocate( ufenlx( u_nloc, cv_ngi ))
    allocate( ufenly( u_nloc, cv_ngi ))
    allocate( ufenlz( u_nloc, cv_ngi ))

    allocate( scvfen( cv_nloc, scvngi ))
    allocate( scvfenslx( cv_nloc, scvngi ))
    allocate( scvfensly( cv_nloc, scvngi ))
    allocate( scvfenlx( cv_nloc, scvngi ))
    allocate( scvfenly( cv_nloc, scvngi ))
    allocate( scvfenlz( cv_nloc, scvngi ))
    allocate( scvfeweigh( scvngi ))

    allocate( sufen( u_nloc, scvngi ))
    allocate( sufenslx( u_nloc, scvngi ))
    allocate( sufensly( u_nloc, scvngi ))
    allocate( sufenlx( u_nloc, scvngi ))
    allocate( sufenly( u_nloc, scvngi ))
    allocate( sufenlz( u_nloc, scvngi ))

    allocate( scvdetwei( scvngi ))
    allocate( sra( scvngi ))
    allocate( log_on_bound(cv_nonods))

    allocate( sbcvfen( cv_snloc, sbcvngi ))
    allocate( sbcvfenslx( cv_snloc, sbcvngi ))
    allocate( sbcvfensly( cv_snloc, sbcvngi ))
    allocate( sbcvfeweigh( sbcvngi ))
    allocate( sbcvfenlx( cv_snloc, sbcvngi ))
    allocate( sbcvfenly( cv_snloc, sbcvngi ))
    allocate( sbcvfenlz( cv_snloc, sbcvngi ))
    allocate( sbufen( u_snloc, sbcvngi ))
    allocate( sbufenslx( u_snloc, sbcvngi ))
    allocate( sbufensly( u_snloc, sbcvngi ))
    allocate( sbufenlx( u_snloc, sbcvngi ))
    allocate( sbufenly( u_snloc, sbcvngi ))
    allocate( sbufenlz( u_snloc, sbcvngi ))

    allocate( cv_sloc2loc( cv_snloc ))
    allocate( u_sloc2loc( u_snloc )) 
    allocate( cv_sloclist( nface, cv_snloc ))
    allocate( u_sloclist( nface, u_snloc ))
    allocate( cv_neiloc( cv_nloc, scvngi ))

    allocate( sele_overlap_scale(cv_nloc) )

    allocate( ugi_coef_ele(u_nloc),  vgi_coef_ele(u_nloc),  wgi_coef_ele(u_nloc) )
    allocate( ugi_coef_ele2(u_nloc), vgi_coef_ele2(u_nloc), wgi_coef_ele2(u_nloc) )
    ! the procity mapped to the cv nodes
    allocate( sum_cv( cv_nonods ))
    allocate( up_wind_nod( cv_nonods * nphase ))
    up_wind_nod = 0.0


    ewrite(3,*)' in mapping overlapping to discontinuous quadratic FE routine ...'

    d1 = ( ndim == 1 )
    d3 = ( ndim == 3 )
    dcyl= ( ndim == -2 )


! taken from assemb_force_cty

! get the shape functions...
    call cv_fem_shape_funs( &
                                ! volume shape functions...
         ndim, p_ele_type,  & 
         cv_ngi, cv_ngi_short, cv_nloc, u_nloc, cvn, cvn_short, &
         cvweight, cvfen, cvfenlx, cvfenly, cvfenlz, &
         cvweight_short, cvfen_short, cvfenlx_short, cvfenly_short, cvfenlz_short, &
         ufen, ufenlx, ufenly, ufenlz, &
                                ! surface of each cv shape functions...
         scvngi, cv_neiloc, cv_on_face,  &  
         scvfen, scvfenslx, scvfensly, scvfeweigh, &
         scvfenlx, scvfenly, scvfenlz,  &
         sufen, sufenslx, sufensly,  &
         sufenlx, sufenly, sufenlz,  &
                                ! surface element shape funcs...
         u_on_face, nface, & 
         sbcvngi,sbcvfen, sbcvfenslx, sbcvfensly, sbcvfeweigh, sbcvfenlx, sbcvfenly, sbcvfenlz, &
         sbufen, sbufenslx, sbufensly, sbufenlx, sbufenly, sbufenlz, &
         cv_sloclist, u_sloclist, cv_snloc, u_snloc, &
                                ! define the gauss points that lie on the surface of the cv...
         findgpts, colgpts, ncolgpts, &
         sele_overlap_scale ) 




    loop_elements: do ele = 1, totele ! volume integral

       matloc=0.0
       rhsloc_u=0.0
       rhsloc_v=0.0
       rhsloc_w=0.0
       rhsloc_uold=0.0
       rhsloc_vold=0.0
       rhsloc_wold=0.0
! taken from proj_cv_to_fem_4: 
       ! calculate detwei,ra,nx,ny,nz for element ele
       call detnlxr( ele, x, y, z, x_ndgln, totele, x_nonods, cv_nloc, cv_ngi, &
            n, nlx, nly, nlz, cvweight, detwei, ra, volume, d1, d3, dcyl, &
            nx, ny, nz ) 


! taken from assemb_force_cty
       ud = 0.0
       vd = 0.0
       wd = 0.0
       udold = 0.0
       vdold = 0.0
       wdold = 0.0
       do u_iloc = 1, u_nloc
          !          write(357,*) 'ele, u_nonods, iloc:',ele, u_nonods, iloc
          u_nod = u_ndgln(( ele - 1 ) * u_nloc + u_iloc )
          do gi = 1, cv_ngi
             do iphase=1,nphase
                u_nod_pha=u_nod +(iphase-1)*u_nonods
                ud( gi, iphase ) = ud( gi, iphase ) + ufen( u_iloc, gi ) * u( u_nod_pha ) 
                vd( gi, iphase ) = vd( gi, iphase ) + ufen( u_iloc, gi ) * v( u_nod_pha ) 
                wd( gi, iphase ) = wd( gi, iphase ) + ufen( u_iloc, gi ) * w( u_nod_pha ) 
                udold( gi, iphase ) = udold( gi, iphase ) + ufen( u_iloc, gi ) * uold( u_nod_pha ) 
                vdold( gi, iphase ) = vdold( gi, iphase ) + ufen( u_iloc, gi ) * vold( u_nod_pha ) 
                wdold( gi, iphase ) = wdold( gi, iphase ) + ufen( u_iloc, gi ) * wold( u_nod_pha ) 
             end do
          end do
       end do

! taken from proj_cv_to_fem_4: 
       loop_cv_iloc: do cv_iloc = 1, cv_nloc

          cv_nodi = cv_ndgln(( ele - 1 ) * cv_nloc + cv_iloc )
          !    write(357,*)'ele,cv_nodi,cv_iloc:',ele,cv_nodi,cv_iloc, x(cv_nodi)

          loop_cv_jloc: do cv_jloc = 1, cv_nloc

             cv_nodj = cv_ndgln(( ele - 1 ) * cv_nloc + cv_jloc )
             nn = 0.0
             do cv_gi = 1, cv_ngi
                nn = nn + n( cv_iloc, cv_gi )   * n(   cv_jloc, cv_gi ) * detwei( cv_gi )
             end do
             matloc( cv_iloc,cv_jloc ) = matloc( cv_iloc,cv_jloc ) + nn

          end do loop_cv_jloc


       end do loop_cv_iloc
! solver here...
        call invert(matloc)

         u_nod = u_ndgln(( ele - 1 ) * u_nloc + u_iloc )
         do iphase=1,nphase
             u_nod_pha=u_nod +(iphase-1)*u_nonods
             rhsloc_u = 0.0
             rhsloc_v = 0.0
             rhsloc_w = 0.0
             rhsloc_uold = 0.0
             rhsloc_vold = 0.0
             rhsloc_wold = 0.0
             do cv_iloc = 1, cv_nloc
              do cv_gi = 1, cv_ngi
                 rhsloc_u(cv_iloc)=rhsloc_u(cv_iloc)+  n( cv_iloc, cv_gi )* ud( gi, iphase ) * detwei( cv_gi )
                 rhsloc_v(cv_iloc)=rhsloc_v(cv_iloc)+  n( cv_iloc, cv_gi )* vd( gi, iphase ) * detwei( cv_gi )
                 rhsloc_w(cv_iloc)=rhsloc_w(cv_iloc)+  n( cv_iloc, cv_gi )* wd( gi, iphase ) * detwei( cv_gi )
                 rhsloc_uold(cv_iloc)=rhsloc_uold(cv_iloc)+  n( cv_iloc, cv_gi )* udold( gi, iphase ) * detwei( cv_gi )
                 rhsloc_vold(cv_iloc)=rhsloc_vold(cv_iloc)+  n( cv_iloc, cv_gi )* vdold( gi, iphase ) * detwei( cv_gi )
                 rhsloc_wold(cv_iloc)=rhsloc_wold(cv_iloc)+  n( cv_iloc, cv_gi )* wdold( gi, iphase ) * detwei( cv_gi )
              end do
             end do
             velocity_dg(cv_nonods+(ele-1)*cv_nloc:cv_nonods+(ele-1)*cv_nloc+cv_nloc,iphase,1) = matmul(matloc,rhsloc_u)
            if(ndim>=2) then
             velocity_dg(cv_nonods+(ele-1)*cv_nloc:cv_nonods+(ele-1)*cv_nloc+cv_nloc,iphase,2) = matmul(matloc,rhsloc_v)
            endif
            if(ndim>=3) then
             velocity_dg(cv_nonods+(ele-1)*cv_nloc:cv_nonods+(ele-1)*cv_nloc+cv_nloc,iphase,3)  = matmul(matloc,rhsloc_w)
            endif
            ! uold_dg( (iphase-1)*cv_nonods+ (ele-1)*cv_nloc:(iphase-1)*cv_nonods+(ele-1)*cv_nloc+cv_nloc) = matmul(matloc,rhsloc_uold)
            if(ndim>=2) then
            ! vold_dg( (iphase-1)*cv_nonods+(ele-1)*cv_nloc:(iphase-1)*cv_nonods+(ele-1)*cv_nloc+cv_nloc ) = matmul(matloc,rhsloc_vold)
            endif
            if(ndim>=3) then
           !  wold_dg( (iphase-1)*cv_nonods+(ele-1)*cv_nloc:(iphase-1)*cv_nonods+(ele-1)*cv_nloc+cv_nloc ) = matmul(matloc,rhsloc_wold)
            endif
        end do

    end do loop_elements


 end subroutine overlapping_to_quadratic_dg

end module mapping_for_ocvfem
