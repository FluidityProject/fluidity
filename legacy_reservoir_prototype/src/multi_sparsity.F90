  
  !    Copyright (C) 2006 Imperial College London and others.
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
  

  module sparsity_1D
    use fldebug

    use spud
    use global_parameters, only: option_path_len

  contains

    SUBROUTINE DEF_SPAR( SEMI_BAND_WID, NONODS, MX_NCOLC, NCOLC, &
         CENTC, FINDC, COLC)
      ! define sparsity...
      ! SEMI_BAND_WID is the semi band width.
      IMPLICIT NONE
      INTEGER, intent( in ) :: SEMI_BAND_WID, NONODS, MX_NCOLC
      INTEGER, intent( inout) :: NCOLC
      INTEGER, DIMENSION( NONODS ), intent( inout) ::  CENTC
      INTEGER, DIMENSION( NONODS + 1 ), intent( inout) :: FINDC
      INTEGER, DIMENSION( MX_NCOLC ), intent( inout) :: COLC

      ! Local variables
      INTEGER :: NOD, II, COL, COUNT

      ewrite(3,*) 'In DEF_SPAR'

      COUNT = 0
      COL = 0
      loop_nods: DO NOD = 1, NONODS

         FINDC( NOD ) = COUNT + 1

         DO II = -SEMI_BAND_WID, SEMI_BAND_WID, 1

            COL = NOD + II

            IF( ( COL >= 1 ) .AND. ( COL <= NONODS ) ) THEN
               COUNT = COUNT + 1
               COLC( COUNT ) = COL
               IF( COL == NOD) CENTC( NOD ) = COUNT
            END IF

         END DO

      END DO loop_nods

      FINDC( NONODS + 1 ) = COUNT + 1

      NCOLC = COUNT

      ewrite(3,*) 'Leaving DEF_SPAR'

      RETURN
    END SUBROUTINE DEF_SPAR

    SUBROUTINE DEF_SPAR_CT_DG( CV_NONODS, MX_NCT, NCT, FINDCT, COLCT, &
         TOTELE, CV_NLOC, U_NLOC, U_NDGLN, U_ELE_TYPE, CV_NDGLN)
      ! define sparsity...
      ! SEMI_BAND_WID is the semi band width.
      IMPLICIT NONE
      INTEGER, intent( in ) :: CV_NONODS, MX_NCT
      INTEGER, intent( inout ) :: NCT
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent (inout ) :: FINDCT
      INTEGER, DIMENSION( MX_NCT ), intent( inout ) :: COLCT
      INTEGER, intent( in ) :: TOTELE, CV_NLOC, U_NLOC, U_ELE_TYPE
      INTEGER, DIMENSION ( U_NLOC * TOTELE ), intent( in ) :: U_NDGLN
      integer, dimension (cv_nloc * totele ), intent( in ) :: cv_ndgln
      ! Local variables...
      INTEGER :: CV_NOD, U_NOD, JLOC, COUNT, ELE, ELE1, ELE2, CV_NODI, CV_ILOC, ILEV, count2, rep
      integer, dimension(cv_nonods) :: cv_ndgln_small
      logical :: repeated, finished_colct
      character( len = option_path_len ) :: overlapping_path 
      logical :: is_overlapping   

      ewrite(3,*) 'In DEF_SPAR_CT_DG'

      is_overlapping = .false.
      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) is_overlapping = .true.

      COUNT = 2
      !! Get condensed form of cv_ndgln, ie without any repeats
      !! This can then be used as the index for the loop below
      cv_ndgln_small=0
      rep = 1
      cv_ndgln_small(1)=cv_ndgln(1)
      do cv_nod=2,size(cv_ndgln)
         repeated = .false.
         do cv_nodi=1,cv_nod-rep
            if (cv_ndgln(cv_nod)==cv_ndgln_small(cv_nodi)) then
               repeated=.true.
               rep=rep+1
            end if
         end do
         if (.not.repeated) then
            cv_ndgln_small(count) = cv_ndgln(cv_nod)
            count = count+1
         end if
      end do

      finished_colct=.false.
      COUNT = 0
      count2 = 1
      ewrite(3,*)'cv1:', size(cv_ndgln), cv_ndgln
      ewrite(3,*)'cv2:', size(cv_ndgln_small), cv_ndgln_small 
      ewrite(3,*)'cvnonods:', cv_nonods
      !stop 23

      IF(CV_NONODS /= CV_NLOC*TOTELE ) THEN
         ! Have a cty CV_NOD
         do while (.not.finished_colct)
            loop_cvnod: DO CV_NOD = 1, CV_NONODS

               cv_nodi = cv_ndgln_small(cv_nod)
               if (cv_nodi==count2) then

                  FINDCT( cv_nodi ) = COUNT + 1
                  ELE1 = 1 + ( CV_NOD - 2 ) / ( CV_NLOC - 1 )
                  ELE2 = 1 + ( CV_NOD - 1 ) / ( CV_NLOC - 1 )
                  ewrite(3,*)'findct:', cv_nod, cv_nodi, findct( cv_nodi )

                  loop_elements: DO ELE = MAX( 1 , ELE1 ), MIN( TOTELE , ELE2 ), 1

                     DO JLOC = 1 ,U_NLOC
                        U_NOD = U_NDGLN(( ELE - 1 ) * U_NLOC + JLOC )
                        COUNT = COUNT + 1
                        !                  ewrite(3,*) u_nod
                        COLCT( COUNT ) = U_NOD
                        !    ewrite(3,*)'colct:', ele, count, colct(count)
                     END DO

                  END DO loop_elements

               end if

            END DO loop_cvnod
            if (count2==cv_nonods) finished_colct = .true.
            count2 = count2+1
         end do
      ELSE
         ! Have a discontinuous CV_NOD
         loop_elements2: DO ELE = 1,TOTELE 
            loop_cviloc: DO CV_ILOC = 1, CV_NLOC

               CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC

               FINDCT( CV_NOD ) = COUNT + 1

               ! IF(U_ELE_TYPE==2) THEN
               if ( is_overlapping ) then
                  IF((CV_ILOC==1).AND.(ELE /= 1)) THEN
                     ELE2=ELE-1
                     JLOC=U_NLOC/CV_NLOC
                     DO ILEV=1,CV_NLOC
                        U_NOD = U_NDGLN(( ELE2 - 1 ) * U_NLOC + JLOC + (ILEV-1)*U_NLOC/CV_NLOC)
                        COUNT = COUNT + 1
                        COLCT( COUNT ) = U_NOD
                     END DO
                  ENDIF

                  DO JLOC = 1 ,U_NLOC
                     U_NOD = U_NDGLN(( ELE - 1 ) * U_NLOC + JLOC )
                     COUNT = COUNT + 1
                     COLCT( COUNT ) = U_NOD
                  END DO

                  IF((CV_ILOC==CV_NLOC).AND.(ELE /= TOTELE)) THEN
                     ELE2=ELE+1
                     JLOC=1
                     DO ILEV=1,CV_NLOC
                        U_NOD = U_NDGLN(( ELE2 - 1 ) * U_NLOC + JLOC + (ILEV-1)*U_NLOC/CV_NLOC)
                        COUNT = COUNT + 1
                        COLCT( COUNT ) = U_NOD
                     END DO
                  ENDIF
               ELSE
                  IF((CV_ILOC==1).AND.(ELE /= 1)) THEN
                     ELE2=ELE-1
                     JLOC=U_NLOC
                     U_NOD = U_NDGLN(( ELE2 - 1 ) * U_NLOC + JLOC )
                     COUNT = COUNT + 1
                     COLCT( COUNT ) = U_NOD
                  ENDIF

                  DO JLOC = 1 ,U_NLOC
                     U_NOD = U_NDGLN(( ELE - 1 ) * U_NLOC + JLOC )
                     COUNT = COUNT + 1
                     COLCT( COUNT ) = U_NOD
                  END DO

                  IF((CV_ILOC==CV_NLOC).AND.(ELE /= TOTELE)) THEN
                     ELE2=ELE+1
                     JLOC=1
                     U_NOD = U_NDGLN(( ELE2 - 1 ) * U_NLOC + JLOC )
                     COUNT = COUNT + 1
                     COLCT( COUNT ) = U_NOD
                  ENDIF
               ENDIF

            END DO loop_cviloc

         END DO loop_elements2
      ENDIF

      FINDCT( CV_NONODS + 1) = COUNT + 1
      NCT = COUNT

      IF(NCT > MX_NCT ) THEN
         EWRITE(3,*)'MX_NCT is not long enough NCT,MX_NCT:',NCT,MX_NCT
      ENDIF

      !DO CV_NODI=1,CV_NONODS
      !   EWRITE(3,*)'CV_NODI=',CV_NODI
      !   EWRITE(3,*)(COLCT( COUNT ),COUNT = FINDCT( CV_NODI ), FINDCT( CV_NODI + 1 ) - 1 )
      !END DO

      !ewrite(3,*)'findct', (findct(cv_nod), cv_nod=1,cv_nonods)
      !ewrite(3,*)'colct', (colct(cv_nod), cv_nod=1,nct)

      ewrite(3,*) 'Leaving DEF_SPAR_CT_DG'

      RETURN
    END SUBROUTINE DEF_SPAR_CT_DG

  end module sparsity_1D

  module sparsity_ND
    use fldebug

  contains

    subroutine getfinele( totele, nloc, snloc, nonods, ndglno, mx_nface_p1, &
         mxnele, ncolele, finele, colele, midele )
      ! This sub caluculates COLELE the element connectivitiy list 
      ! in order of faces.
      implicit none
      integer, intent( in ) :: totele, nloc, snloc, nonods
      integer, dimension( totele * nloc ), intent( in ) :: ndglno
      integer, intent( in ) :: mx_nface_p1, mxnele
      integer, intent( inout ) :: ncolele
      integer, dimension( mxnele ), intent( inout ) :: colele
      integer, dimension( totele + 1 ), intent( inout ) :: finele
      integer, dimension( totele ), intent( inout ) :: midele
      ! Local variables
      integer :: ele, iloc, jloc, iloc2, nod, inod, jnod, count, ele2, i, hit, &
           iface, iface2, kface, ncolel, itemp, count2
      logical :: found
      integer, allocatable, dimension( : ) :: fintran, coltran, icount

      allocate( fintran( nonods + 1 ))
      allocate( coltran( max( totele, nonods ) * mx_nface_p1 ))
      allocate( icount( nonods * totele ))

      icount = 0
      do ele = 1, totele
         do iloc = 1, nloc
            nod = ndglno( ( ele - 1 ) * nloc + iloc )
            icount( nod ) = icount( nod ) + 1
         end do
      end do

      fintran = 0
      fintran( 1 ) = 1
      do nod = 1, nonods
         fintran( nod + 1 ) = fintran( nod ) + icount( nod )
      end do

      icount = 0
      coltran = 0
      do ele = 1, totele
         do iloc = 1, nloc
            nod = ndglno( ( ele - 1 ) * nloc + iloc )
            !ewrite(3,*)'nod, filtran, icount, ele:', nod, fintran( nod ), ele, totele, nonods
            coltran( fintran( nod ) + icount( nod )) = ele
            icount( nod ) = icount( nod ) + 1
         end do
      end do
      !ewrite(3,*)'coltran:', coltran( 1: max(totele, nonods ) * mx_nface_p1 )
      !ewrite(3,*)'fintran:', fintran( 1: nonods + 1 )
      !ewrite(3,*)'X_NDGLN:', ndglno( 1: totele*nloc )

      icount = 0 ; colele = 0 ; ncolele = 0
      Loop_Elements1: do ele = 1, totele

         Loop_Iloc: do iloc = 1, nloc

            nod = ndglno( ( ele - 1 ) * nloc + iloc )
            Loop_Count1: do count = fintran( nod ), fintran( nod + 1 ) - 1, 1

               ele2 = coltran( count )
               found = .false. ! Add ELE2 into list FINELE and COLELE
               do i = 1, icount( ele )
                  if( colele( ( ele - 1 ) * mx_nface_p1 + i ) == ele2 ) found = .true.
               end do

               Conditional_Found: if ( .not. found ) then ! Do elements ELE and ELE2 share at least 3 nodes?

                  hit = 0
                  do iloc2 = 1, nloc
                     inod = ndglno( ( ele - 1 ) * nloc + iloc2 )
                     do jloc = 1, nloc
                        jnod = ndglno( ( ele2 - 1 ) * nloc + jloc )
                        if ( inod == jnod ) hit = hit + 1
                     end do
                  end do
                  if ( hit >= snloc ) then
                     icount( ele ) = icount( ele ) + 1
                     colele( ( ele - 1 ) * mx_nface_p1 + icount( ele )) = ele2 
                     ncolele = ncolele + 1 
                  end if

               end if Conditional_Found

            end do Loop_Count1

         end do Loop_Iloc

      end do Loop_Elements1

      finele( 1 ) = 1
      do ele = 1, totele
         finele( ele + 1 ) = finele( ele ) + icount( ele )
      end do

      ! order elements in increasing order...
      count = 0
      Loop_Elements2: do ele = 1, totele
         ! Shorten COLELE then perform a bubble sort to get the ordering right for.
         do iface = 1, mx_nface_p1
            if ( colele( ( ele - 1 ) * mx_nface_p1 + iface ) /= 0 ) then
               count = count + 1
               colele( count ) = colele( ( ele - 1 ) * mx_nface_p1 + iface )
            end if
         end do
      end do Loop_Elements2

      Loop_BubbleSort: do ele = 1, totele
         do count = finele( ele ) , finele( ele + 1 ) - 2
            do count2 = finele( ele ) , finele( ele + 1 ) - 1, 1
               if ( colele( count ) > colele( count + 1 )) then ! swop over
                  itemp = colele( count + 1 )
                  colele( count + 1 ) = colele( count )
                  colele( count ) = itemp
               end if
            end do
         end do
      end do Loop_BubbleSort

      ! Calculate midele: 
      do ele = 1, totele
         do count = finele( ele ) , finele( ele + 1 ) - 1
           if(colele(count)==ele) midele( ele ) = count
         end do
      end do

      deallocate( fintran )
      deallocate( coltran )
      deallocate( icount )

      return
    end subroutine getfinele

    subroutine exten_sparse_multi_phase( nonods, mxnele, finm, colm, &
         nphase, npha_nonods, ncolm_pha, &
         finm_pha, colm_pha, midm_pha )
      ! Extend the sparsity to a multiphase sparsity
      implicit none
      integer, intent( in ) :: nonods, mxnele
      integer, dimension( nonods + 1 ), intent( in ) :: finm
      integer, dimension( mxnele ), intent( in ) :: colm
      integer, intent( in ) :: nphase, npha_nonods, ncolm_pha
      integer, dimension( npha_nonods + 1 ), intent( inout ) :: finm_pha
      integer, dimension( ncolm_pha ), intent( inout ) :: colm_pha
      integer, dimension( npha_nonods ), intent( inout ) :: midm_pha
      ! Local variables
      integer :: count, count2, iphase, jphase, nod

      ewrite(3,*) 'In exten_sparse_multi_phase subrt.'

      count2 = 0
      Loop_Phase1: do iphase = 1, nphase
         Loop_CVNODS: do nod = 1, nonods
            Loop_Phase2: do jphase = 1, nphase
               if( jphase == 1 ) &
                    finm_pha( ( iphase - 1 ) * nonods + nod ) = count2 + 1
               Conditional_Phases: if( iphase == jphase ) then
                  do count = finm( nod ), finm( nod + 1 ) - 1
                     count2 = count2 + 1
                     colm_pha( count2 ) = colm( count ) + ( jphase - 1 ) * nonods
                     if( colm( count ) == nod ) &
                          midm_pha( nod + ( jphase - 1 ) * nonods ) = count2
                  end do
               else
                  count2 = count2 + 1
                  colm_pha( count2 ) = nod + ( jphase - 1 ) * nonods
               end if Conditional_Phases
            end do Loop_Phase2
         end do Loop_CVNODS
      end do Loop_Phase1

      finm_pha( nphase * nonods + 1 ) = count2 + 1
!      ncolm_pha= count2 
       if(count2.ne.ncolm_pha) then
          ewrite(3,*) 'not correct length count2,ncolm_pha:',count2,ncolm_pha
          stop 2821
       end if

      !ewrite(3,*) 'colm_pha--',colm_pha(1:ncolm_pha)

      ewrite(3,*) 'Leaving exten_sparse_multi_phase subrt.'
      return
    end subroutine exten_sparse_multi_phase

    subroutine exten_sparse_mom_cty( ndim, finmcy2, colmcy2, nlenmcy2, nmcy2, mxnct, &
         cv_nonods, findct, colct, &
         u_nonods, ncolc, mx_ncolmcy, &
         findc, colc, finmcy, colmcy, midmcy, nlenmcy, &
         ncolmcy, nphase, ncolcmc, findcmc, colcmc )
      ! Extend momentum sparsity to include Continuity / pressure.
      implicit none
      integer, intent( in ) :: ndim, nlenmcy2, nmcy2, mxnct
      integer, dimension( nlenmcy2 + 1 ), intent( in ) :: finmcy2
      integer, dimension( nmcy2 ), intent( in ) :: colmcy2
      integer, intent( in ) :: cv_nonods
      integer, dimension( cv_nonods + 1 ), intent( in ) :: findct
      integer, dimension( mxnct ), intent( in ) :: colct
      integer, intent( in ) :: u_nonods, ncolc, mx_ncolmcy
      integer, dimension( u_nonods + 1 ), intent( in ) :: findc
      integer, dimension( ncolc ), intent( in ) :: colc
      integer, intent( in ) :: nlenmcy
      integer, dimension( nlenmcy + 1 ), intent( inout ) :: finmcy
      integer, dimension( mx_ncolmcy ), intent( inout ) :: colmcy
      integer, dimension( nlenmcy ), intent( inout ) :: midmcy
      integer, intent( inout ) :: ncolmcy
      integer, intent( in ) :: nphase, ncolcmc
      integer, dimension( cv_nonods + 1 ), intent( in ) :: findcmc
      integer, dimension( ncolcmc ), intent( in ) :: colcmc
      ! Local variables
      integer :: count, count2, iphase, jphase, u_nod, u_pha_nod, cv_nod, icol, jdim

      ewrite(3,*) 'In exten_sparse_mom_cty subrt.'

      count2 = 0 ; count = 0
      Loop_Phase1: do iphase = 1, nphase ! Add the pressure part in
         Loop_Nod1: do u_nod = 1, u_nonods
            u_pha_nod = u_nod + ( iphase - 1 ) * u_nonods
            finmcy( u_pha_nod ) = count2 + 1
            Loop_Count1a: do count = finmcy2( u_pha_nod ), finmcy2( u_pha_nod + 1 ) - 1
               count2 = count2 + 1
               colmcy( count2 ) = colmcy2( count )
            end do Loop_Count1a
            Loop_Count2: do count = findc( u_nod ), findc( u_nod + 1 ) - 1
               count2 = count2 + 1
               colmcy( count2 ) = colc( count ) + u_nonods * nphase
            end do Loop_Count2
         end do Loop_Nod1
      end do Loop_Phase1

      Loop_Nod2a: do cv_nod = 1, cv_nonods ! Add continuity part in
         u_pha_nod = cv_nod + u_nonods * nphase
         finmcy( u_pha_nod ) = count2 + 1

         Loop_Phase2: do jphase = 1, nphase
            Loop_Dim2: do jdim = 1, ndim
               Loop_Nod2b: do count = findct( cv_nod ), findct( cv_nod + 1 ) - 1
                  count2 = count2 + 1
                  icol = colct( count ) + ( jdim - 1 ) * u_nonods + &
                       ( jphase - 1 ) * u_nonods * ndim
                  colmcy( count2 ) = icol
               end do Loop_Nod2b
            end do Loop_Dim2
         end do Loop_Phase2

         ! Now add a diagonal which will have zero values for incompressible, 
         ! but non-zero for compressible
         Loop_Nod2c: do count = findcmc( cv_nod ), findcmc( cv_nod + 1 ) - 1
            count2 = count2 + 1
            icol = colcmc( count ) + u_nonods * ndim * nphase
            colmcy( count2 ) = icol
            if( u_pha_nod == icol ) midmcy( u_pha_nod ) = count2
         end do Loop_Nod2c

      end do Loop_Nod2a
      ncolmcy = count2
      finmcy( u_nonods * nphase + cv_nonods + 1 ) = count2 + 1

      ewrite(3,*) 'Leaving exten_sparse_mom_cty subrt.'

      return
    end subroutine exten_sparse_mom_cty

    subroutine form_dgm_pha_sparsity( totele, nphase, u_nloc, u_pha_nonods, &
         ndim, mx_ncoldgm_pha, ncoldgm_pha, &
         coldgm_pha, findgm_pha, middgm_pha, &
         nele_pha, finele_pha, colele_pha, finele, colele, ncolele )
      ! Form the sparsity of the phase coupled DG discretised matrix 
      ! from the element-wise multi-phase sparsity matrix. 
      implicit none
      integer, intent( in ) :: totele, nphase, u_nloc, u_pha_nonods, &
           mx_ncoldgm_pha, ndim, ncolele
      integer, intent( inout ) :: ncoldgm_pha
      integer, dimension( mx_ncoldgm_pha ), intent( inout ) :: coldgm_pha
      integer, dimension( u_pha_nonods + 1 ), intent( inout ) :: findgm_pha
      integer, dimension( u_pha_nonods ), intent( inout ) :: middgm_pha
      integer, intent( in ) :: nele_pha
      integer, dimension( nele_pha + 1 ), intent( in ) :: finele_pha
      integer, dimension( nele_pha ), intent( in ) :: colele_pha
      integer, dimension( totele + 1 ), intent( in ) :: finele
      integer, dimension( ncolele ), intent( in ) :: colele

      ! Local variables
      integer :: count, count2, ele_pha, ele_pha2, iloc, jloc, irow, jrow
      integer :: ele, ele2, idim, jdim, iphase, jphase, u_nonods

      ewrite(3,*) 'In form_dgm_pha_sparsity subrt.'

      u_nonods = u_pha_nonods / ( nphase * ndim )

      count2 = 0
      Loop_Element: do ele = 1, totele 
         Loop_Phase1: do iphase = 1, nphase
            Loop_Loc1: do iloc = 1, u_nloc
               Loop_Dim1: do idim = 1, ndim
                  irow = ( ele - 1 ) * u_nloc + iloc  + ( idim - 1 ) * u_nonods + (iphase-1)*u_nonods*ndim
                  !print *, 'irow, ele, u_nloc, iloc, idim, u_nonods, iphase:', &
                  !     irow, ele, u_nloc, iloc, idim, u_nonods, iphase
                  findgm_pha( irow ) = count2 + 1
                  Loop_Count: do count = finele( ele ), finele( ele + 1 ) - 1
                     ele2 = colele( count )
                     Loop_Phase2: do jphase = 1, nphase
                        Loop_Loc2: do jloc = 1, u_nloc
                           Loop_Dim2: do jdim = 1, ndim
                              jrow = ( ele2 - 1 ) * u_nloc + jloc  + ( jdim - 1 ) * u_nonods + (jphase-1)*u_nonods*ndim
                              count2 = count2 + 1
                              coldgm_pha( count2 ) = jrow
                              if( irow == jrow ) middgm_pha( irow ) = count2
                           end do Loop_Dim2
                        end do Loop_Loc2
                     end do Loop_Phase2
                  end do Loop_Count
               end do Loop_Dim1
            end do Loop_Loc1
         end do Loop_Phase1
         findgm_pha( u_pha_nonods + 1 ) = count2 + 1
         ncoldgm_pha = count2
      end do Loop_Element


      if( ncoldgm_pha > mx_ncoldgm_pha ) &
           FLAbort(" Incorrect number of dimension of sparsity matrix - ncoldgm_pha ")

      ewrite(3,*) 'Leaving form_dgm_pha_sparsity subrt. '

      return
    end subroutine form_dgm_pha_sparsity


    subroutine pousinmc2( totele, nonods1, nloc1, nonods2, nloc2, &
         nimem, ndglno1, ndglno2, &
         lencolm, findrm, colm, centrm )
      implicit none
      integer, intent( in ) :: totele, nonods1, nloc1, nonods2, nloc2, nimem
      integer, dimension( totele * nloc1 ), intent( in ) :: ndglno1
      integer, dimension( totele * nloc2 ), intent( in ) :: ndglno2
      integer, intent( inout ) :: lencolm
      integer, dimension( nonods2 + 1 ), intent( inout ) :: findrm
      integer, dimension( nimem ), intent( inout ) :: colm
      integer, dimension( nonods2 ), intent( inout ) :: centrm

      ! Local Variables
      integer :: ele, globi, globj, loci, locj, i, irow, ptr

      ! Defining derived data types for the linked list
      type node
         integer :: id                 ! id number of node
         type( node ), pointer :: next ! next node
      end type node

      type row
         type( node ), pointer :: row ! recursive data type
      end type row

      type( row ), dimension( : ), allocatable :: matrix
      type( node ), pointer :: list, current, next

      allocate( matrix( nonods2 ))
      do i = 1, nonods2
         allocate( list )
         list % id = -1
         nullify( list % next )
         matrix( i ) % row => list
         nullify( list )
      end do

      Loop_Elements1: do ele = 1, totele

         Loop_LocI: do loci = 1, nloc2

            globi = ndglno2( ( ele - 1 ) * nloc2 + loci )
            list => matrix( globi ) % row

            Loop_LocJ: do locj = 1, nloc1

               globj = ndglno1( ( ele - 1 ) * nloc1 + locj )

               if ( list % id == -1 ) then ! Check if the list is initalised
                  list % id = globj
                  cycle
               end if

               !      if ( list % id == -1 ) then ! Check if the list is initalised
               !         list % id = globj
               !         cycle
               !      end if

               Conditional1: if ( globj < list % id ) then ! Insert at start of list
                  allocate( current )
                  current % id = globj
                  current % next => list
                  matrix( globi ) % row => current
                  list => matrix( globi ) % row

               else ! Conditional1
                  current => list
                  Loop_While1: do while ( associated( current ))

                     if ( globj == current % id ) then ! Already have this node
                        exit

                     elseif ( .not. associated( current % next )) then  ! End of list - insert this node
                        allocate( current % next )
                        nullify( current % next % next )
                        current % next % id = globj
                        exit

                     elseif ( globj < current % next % id ) then ! Insert new node here
                        allocate( next )
                        next % id = globj
                        next % next => current % next
                        current % next => next
                        exit                

                     end if
                     current => current % next

                  end do Loop_While1

               end if Conditional1

            end do Loop_LocJ

         end do Loop_LocI

      end do Loop_Elements1

      ! From matrix write COLM, FINDRM and CENTRM
      ! linked list as we go
      ptr = 1
      ewrite(3,*),'nonods2=',nonods2
       
      Loop_Irow: do irow = 1, nonods2
         findrm( irow ) = ptr
         centrm( irow ) = -1

         current => matrix( irow ) % row
         Loop_While2: do while ( associated( current ))
            assert( ptr <= nimem )
            colm( ptr ) = current % id
            if ( current % id == -1 ) then
               ewrite(0,*) "ERROR: POSINM() seriously unhappy with node", IROW
               FLAbort( "ERROR: Mesh contains nodes that are not associated with any elements." )
            end if
            if ( current % id == irow ) then
               centrm( irow ) = ptr
            endif
            next => current % next
            deallocate( current )
            current => next
            ptr = ptr + 1
         end do Loop_While2

      end do Loop_Irow

      lencolm = ptr - 1
      findrm( nonods2 + 1 ) = lencolm + 1

      deallocate( matrix )

      return
    end subroutine pousinmc2


    subroutine conv_ct2c( cv_nonods, nct, findct, colct, u_nonods, &
         mx_nc, findc, colc )
      implicit none
      integer, intent( in ) :: cv_nonods, nct
      integer, dimension( cv_nonods + 1 ), intent( in ) :: findct
      integer, dimension( nct ), intent( in ) :: colct
      integer, intent( in ) :: u_nonods, mx_nc
      integer, dimension( u_nonods + 1 ), intent( inout ) :: findc
      integer, dimension( mx_nc ), intent( inout ) :: colc
      ! Local variables
      integer :: col, u_nod, cv_nod, count, count2
      integer, dimension( : ), allocatable :: in_row_c

      ewrite(3,*) 'In conv_ct2c subrt.'
      ewrite(3,*) ' max( colct ):', maxval( colct )

      ! No. of non-zero's in row of C matrix. 
      allocate( in_row_c( u_nonods ) )
      in_row_c = 0
      do count = 1, nct
         col = colct( count )
         in_row_c( col ) = in_row_c( col ) + 1
      end do

      count = 0
      do u_nod = 1, u_nonods
         findc( u_nod ) = count + 1
         count = count + in_row_c( u_nod )
      end do
      findc( u_nonods + 1 ) = count + 1
      in_row_c( 1 : u_nonods ) = 0

      Loop_CVNOD: do cv_nod = 1, cv_nonods
         Loop_Count: do count = findct( cv_nod ), findct( cv_nod + 1 ) - 1
            col = colct( count )
            in_row_c( col ) = in_row_c( col ) + 1
            count2 = findc( col ) + in_row_c( col ) - 1
            colc( count2 ) = cv_nod
         end do Loop_Count
      end do Loop_CVNOD

      deallocate( in_row_c )

      ewrite(3,*) 'Leaving conv_ct2c subrt.'

      return
    end subroutine conv_ct2c


    subroutine poscmc( totele, nonods, nimem, nct, &
         findct, colct, &
         ncmc, fincmc, colcmc, midcmc, noinod, presym )
      ! This subroutine forms the matrix operating on the pressure vector.
      ! It is found from C1T ML C1 + C2T ML C2
      ! In the first part of COLCMC contains the pressure nodes surrounding
      ! a given node.
      implicit none
      integer, intent ( in ) :: totele, nonods, nimem, nct
      integer, dimension( totele + 1 ), intent( in ) :: findct
      integer, dimension( nct ), intent( in ) :: colct
      integer, intent( inout ) :: ncmc
      integer, dimension( totele + 1 ), intent( inout ) :: fincmc
      integer, dimension( nimem ), intent( inout ) :: colcmc
      integer, dimension( totele ), intent( inout ) :: midcmc
      integer, dimension( nonods ), intent( inout ) ::  noinod
      logical, intent( inout ) :: presym

      ! Local variables
      integer :: count, count2, globi, globj, i, irow, inod, jnod, jrow, ptr 

      ! Defining derived data types for the linked list
      type node
         integer :: id                 ! id number of node
         type( node ), pointer :: next ! next node, recursive data type
      end type node

      type row
         type( node ), pointer :: row
      end type row

      type( row ), dimension( : ), allocatable :: matrix, matrix2
      type( node ), pointer :: list, current, current2, next

      allocate( matrix( nonods ))
      do i = 1, nonods
         allocate( list )
         list % id = -1
         nullify( list % next )
         matrix( i ) % row => list
         nullify( list )
      end do

      noinod = 0

      Loop_Row1: do irow = 1, totele ! Given a vel node find the pressure nodes surrounding it 
         globj = irow
         Loop_Count1: do count = findct( irow ), findct( irow + 1 ) - 1
            globi = colct( count )
            list => matrix( globi ) % row
            ! ewrite(3,*)'list%id:',globj,globi,list % id

            if ( list % id == -1 ) then ! Check if the list is initalised
               list % id = globj
               cycle
            end if

            Conditional1: if ( globj < list % id ) then ! Insert at start of list
               allocate( current )
               current % id = globj
               current % next => list
               matrix( globj ) % row => current
               list => matrix( globj ) % row
            else
               current => list
               !ewrite(3,*)'-->', globj, current % id, current % next % id
               Loop_While1: do while( associated( current ))
                  if ( globj == current % id ) then ! Already have this node
                     exit
                  elseif ( .not. associated( current % next )) then ! End of list - insert this node
                     allocate( current % next )
                     nullify( current % next % next )
                     current % next % id = globj
                     exit
                  elseif( globj < current % next % id ) then ! Insert new node here
                     allocate( next )
                     next % id = globj
                     next % next => current % next
                     current % next => next
                     exit
                  end if
                  current => current % next
               end do Loop_While1
            end if Conditional1

            noinod( globi ) = noinod( globi ) + 1

         end do Loop_Count1

      end do Loop_Row1

      allocate( matrix2( totele )) ! Initalise the linked lists
      do i = 1, totele
         allocate( list )
         list % id = -1
         nullify( list % next )
         matrix2( i ) % row => list
         nullify( list )
      end do

      Loop_Row2: do irow = 1, totele
         ! Find the pressure nodes surrounding pressure node IROW
         Loop_Count2: do count = findct( irow ), findct( irow + 1 ) - 1
            inod = colct( count )

            ! Find the pressure nodes surrounding node INOD 
            ! these will be connected to pressure node IROW.
            current => matrix( inod ) % row
            Loop_While2: do while( associated( current ))
               jrow = current % id

               Conditional2: if (( .not. presym ) .or. ( jrow >= irow )) then
                  list => matrix2( irow ) % row

                  if ( list % id == -1 ) then ! Check if the list is initialised
                     list % id = jrow
                     cycle
                  end if

                  Conditional3: if ( jrow < list % id ) then ! Insert at start of list
                     allocate( current2 )
                     current2 % id = jrow
                     current2 % next => list
                     matrix2( irow ) % row => current2
                     list => matrix2( irow ) % row

                  else
                     current2 => list
                     Loop_While3: do while( associated( current2 ))

                        Conditional4: if ( jrow == current2 % id ) then ! Already have this node
                           exit
                        elseif( .not. associated( current2 % next )) then ! End of list - insert this node
                           allocate( current2 % next )
                           nullify(  current2 % next % next )
                           current2 % next % id = jrow
                           exit 
                        elseif( jrow < current2 % next % id ) then ! Insert new node here
                           allocate( next )
                           next % id = jrow
                           next % next => current2 % next
                           current2 % next => next
                           exit
                        end if Conditional4
                        current2 => current2 % next

                     end do Loop_While3

                  end if Conditional3

               end if Conditional2
               current => current % next          

            end do Loop_While2

         end do Loop_Count2

      end do Loop_Row2


      do irow = 1, nonods ! Delete Matrix
         current => matrix( irow ) % row
         do while ( associated( current ))
            next => current % next
            deallocate( current )
            current => next
         end do
      end do
      deallocate( matrix )

      ptr = 1
      Loop_Row3: do irow = 1, totele ! From matrix write COLCMC, FINCMC and MIDCMC
         midcmc( irow ) = -1
         fincmc( irow ) = ptr
         current => matrix2( irow ) % row ! Warning: overwriten 

         Loop_While4: do while ( associated( current ))

            if ( ptr > nimem ) then
               ewrite( -1, * ) 'nimem, ptr: ',nimem, ptr
               ewrite( -1, * ) 'totele, irow: ',totele, irow
               FLAbort( "Integer memory too small" )
            end if

            colcmc( ptr ) = current % id
            if( current % id == irow ) then
               midcmc( irow ) = ptr
            end if

            next => current % next
            deallocate( current )
            current => next
            ptr = ptr + 1

         end do Loop_While4
      end do Loop_Row3

      ncmc =  ptr - 1
      fincmc( totele + 1 ) = ncmc + 1

      return
    end subroutine poscmc

    subroutine CV_Neighboor_Sparsity( ndim, nphase, cv_ele_type, &
         totele, cv_nloc, u_nloc, x_nloc, xu_nloc, mat_nloc, &
         cv_snloc, u_snloc, cv_nonods, x_nonods, &
         cv_ndgln, x_ndgln, xu_ndgln, &
         ncolele, finele, colele, &
         ncolm, mxnacv_loc, findm, colm, &
         ncolacv_loc, finacv_loc, colacv_loc, midacv_loc )
      use shape_functions
      use shape_functions_Linear_Quadratic
      use cv_advection
      implicit none
      integer, intent( in ) :: ndim, nphase, cv_ele_type, totele, cv_nloc, &
           u_nloc, x_nloc, xu_nloc, mat_nloc, cv_snloc, u_snloc, cv_nonods, &
           x_nonods
      integer, dimension( totele * cv_nloc ), intent( in ) :: cv_ndgln
      integer, dimension( totele * x_nloc ), intent( in ) :: x_ndgln
      integer, dimension( totele * xu_nloc ), intent( in ) :: xu_ndgln
      integer, intent( in ) :: ncolele
      integer, dimension( totele + 1 ), intent( in ) :: finele
      integer, dimension( ncolele ), intent( in ) :: colele
      integer, intent( in ) :: ncolm, mxnacv_loc
      integer, dimension( cv_nonods + 1 ), intent( in ) :: findm
      integer, dimension( ncolm ), intent( in ) :: colm
      integer, intent( inout ) :: ncolacv_loc
      integer, dimension( cv_nonods + 1 ), intent( inout ) :: finacv_loc
      !  integer, dimension( mxnacv_loc ), intent( inout ) :: colacv_loc
      integer, dimension( mxnacv_loc ), intent( inout ) :: colacv_loc
      integer, dimension( cv_nonods ), intent( inout ) :: midacv_loc
      ! Local variables
      logical, dimension( : ), allocatable :: found, x_share
      logical, dimension( :, : ), allocatable :: cv_on_face, u_on_face, &
                                       cvfem_on_face, ufem_on_face
      integer, dimension( : ), allocatable :: findgpts, colgpts, cv_other_loc, &
           u_other_loc, mat_other_loc
      integer, dimension( :, : ), allocatable :: cv_neiloc, cvfem_neiloc, cv_sloclist, u_sloclist
      real, dimension( : ), allocatable :: cvweight, cvweight_short, scvfeweight, &
           sbcvfeweigh, sele_overlap_scale
      real, dimension( :, : ),allocatable :: cvn, cvn_short, cvfen, cvfenlx, cvfenly, &
           cvfenlz, cvfen_short, cvfenlx_short, cvfenly_short, cvfenlz_short, &
           scvfen, scvfenslx, scvfensly, scvfenlx, scvfenly, scvfenlz, &
           ufen, ufenlx, ufenly, ufenlz, &
           sufen, sufenslx, sufensly, sufenlx, sufenly, sufenlz, &
           sbcvfen, sbcvfenslx, sbcvfensly,sbcvfenlx, sbcvfenly, sbcvfenlz, &
           sbufen, sbufenslx, sbufensly, sbufenlx, sbufenly, sbufenlz
      logical :: integrat_at_gi
      integer :: scvngi, cv_ngi, cv_ngi_short, sbcvngi, nface, &
           ncolgpts, i_dummy, ele, ele2, ele3, cv_iloc, cv_jloc, cv_nodi, cv_nodj, &
           cv_nodj2, gcount, gcount2, gcount3, gi, cv_nod, count2, count, iface, &
           x_nodi, x_nodi2, mxnele, FACE_COUNT, cv_nodi2, cv_iloc2

      ! Computing Gauss points and array containing node points on neighboors elements
      call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
           cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, .false. )
      allocate( cv_neiloc( cv_nloc, scvngi ) ) 
      allocate( cvfem_neiloc( cv_nloc, scvngi ) ) 
!!      allocate( cv_neiloc( cv_nloc, scvngi ) ) ; cv_neiloc = 0
!!      allocate( cvfem_neiloc( cv_nloc, scvngi ) ) ; cvfem_neiloc = 0
!!!      call volnei( cv_neiloc, cvfem_neiloc, cv_nloc, scvngi, cv_ele_type )

      ! Allocating space
      allocate( found( ncolm ) )
      allocate( cvweight( cv_ngi ) )
      allocate( cvn( cv_nloc, cv_ngi ) )
      allocate( cvn_short( cv_nloc, cv_ngi_short ) )
      allocate( cvfen( cv_nloc, cv_ngi ) )
      allocate( cvfenlx( cv_nloc, cv_ngi ) )
      allocate( cvfenly( cv_nloc, cv_ngi ) )
      allocate( cvfenlz( cv_nloc, cv_ngi ) )
      allocate( cvweight_short( cv_ngi_short ) )
      allocate( cvfen_short( cv_nloc, cv_ngi_short ) )
      allocate( cvfenlx_short( cv_nloc, cv_ngi_short ) )
      allocate( cvfenly_short( cv_nloc, cv_ngi_short ) )
      allocate( cvfenlz_short( cv_nloc, cv_ngi_short ) )
      allocate( ufen( u_nloc, cv_ngi ) )
      allocate( ufenlx( u_nloc, cv_ngi ) )
      allocate( ufenly( u_nloc, cv_ngi ) )
      allocate( ufenlz( u_nloc, cv_ngi ) )
      allocate( cv_on_face( cv_nloc, scvngi ) )
      allocate( cvfem_on_face( cv_nloc, scvngi ) )
      allocate( scvfen( cv_nloc, scvngi ) )
      allocate( scvfenslx( cv_nloc, scvngi ) )
      allocate( scvfensly( cv_nloc, scvngi ) )
      allocate( scvfeweight( scvngi ) )
      allocate( scvfenlx( cv_nloc, scvngi ) )
      allocate( scvfenly( cv_nloc, scvngi ) )
      allocate( scvfenlz( cv_nloc, scvngi ) )
      allocate( sufen( u_nloc, scvngi ) )
      allocate( sufenslx( u_nloc, scvngi ) )
      allocate( sufensly( u_nloc, scvngi ) )
      allocate( sufenlx( u_nloc, scvngi ) )
      allocate( sufenly( u_nloc, scvngi ) )
      allocate( sufenlz( u_nloc, scvngi ) )
      allocate( u_on_face( u_nloc, scvngi ) )
      allocate( ufem_on_face( u_nloc, scvngi ) )
      allocate( sbcvfen( cv_snloc, sbcvngi ) )
      allocate( sbcvfenslx( cv_snloc, sbcvngi ) )
      allocate( sbcvfensly( cv_snloc, sbcvngi ) )
      allocate( sbcvfeweigh( sbcvngi ) )
      allocate( sbcvfenlx( cv_snloc, sbcvngi ) )
      allocate( sbcvfenly( cv_snloc, sbcvngi ) )
      allocate( sbcvfenlz( cv_snloc, sbcvngi ) )
      allocate( sbufen( u_snloc, sbcvngi ) )
      allocate( sbufenslx( u_snloc, sbcvngi ) )
      allocate( sbufensly( u_snloc, sbcvngi ) )
      allocate( sbufenlx( u_snloc, sbcvngi ) )
      allocate( sbufenly( u_snloc, sbcvngi ) )
      allocate( sbufenlz( u_snloc, sbcvngi ) )
      allocate( cv_sloclist( nface, cv_snloc ) )
      allocate( u_sloclist( nface, u_snloc ) )
      allocate( findgpts( cv_nloc + 1 ) )
      allocate( colgpts( cv_nloc * scvngi ) )
      allocate( cv_other_loc( cv_nloc ) )
      allocate( u_other_loc( u_nloc ) )
      allocate( mat_other_loc( mat_nloc ) )
      allocate( x_share( x_nonods ) )
      allocate( sele_overlap_scale( cv_nloc ) )

      ncolgpts = 0 ; colgpts = 0 ; findgpts = 0
      call cv_fem_shape_funs( &
           ndim, cv_ele_type, &
           cv_ngi, cv_ngi_short, cv_nloc, u_nloc, cvn, cvn_short, &
                                ! Volume shape functions
           cvweight, cvfen, cvfenlx, cvfenly, cvfenlz, &
           cvweight_short, cvfen_short, cvfenlx_short, cvfenly_short, cvfenlz_short, &
           ufen, ufenlx, ufenly, ufenlz, &
                                ! Surface of each CV shape functions
           scvngi, cv_neiloc, cv_on_face, cvfem_on_face,&
           scvfen, scvfenslx, scvfensly, scvfeweight, &
           scvfenlx, scvfenly, scvfenlz, &
           sufen, sufenslx, sufensly, &
           sufenlx, sufenly, sufenlz, &
                                ! Surface element shape funcs
           u_on_face,ufem_on_face, nface, &
           sbcvngi, sbcvfen, sbcvfenslx, sbcvfensly, sbcvfeweigh, sbcvfenlx, sbcvfenly, sbcvfenlz, &
           sbufen, sbufenslx, sbufensly, sbufenlx, sbufenly, sbufenlz, &
           cv_sloclist, u_sloclist, cv_snloc, u_snloc, &
                                ! Define the gauss points that lie on the surface of the CV
           findgpts, colgpts, ncolgpts, &
           sele_overlap_scale, .false. )   
      ewrite(3,*)'findgpts:', size( findgpts ), '==>', findgpts( 1: cv_nloc + 1 )
      ewrite(3,*)'colgpts:', size( colgpts ), ncolgpts, '==>', colgpts( 1 : ncolgpts )

      found=.false.
! set the diagonal to true...
      do cv_nodi=1,cv_nonods
         do gcount = findm( cv_nodi ), findm( cv_nodi + 1 ) - 1
            if(colm(gcount)==cv_nodi) found( gcount ) = .true.
         end do 
      end do
! now the off diagonal terms...
      Loop_Elements_1: do ele = 1, totele
         Loop_CVILOC_1: do cv_iloc = 1, cv_nloc
            cv_nodi = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )
            Loop_GI_1: do gi = 1, scvngi
               cv_jloc = cv_neiloc( cv_iloc, gi ) 
               if ( cv_jloc  > 0 ) then
                  cv_nodj = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_jloc )
                  Loop_GCOUNT_1: do gcount = findm( cv_nodi ), findm( cv_nodi + 1 ) - 1
                     cv_nodj2 = colm( gcount )
                     if ( cv_nodj == cv_nodj2 ) found( gcount ) = .true.
                  end do Loop_GCOUNT_1
               endif
            end do Loop_GI_1
         end do Loop_CVILOC_1
      end do Loop_Elements_1

! for discontinuous elements...
      if(totele*cv_nloc==cv_nonods) then

         Loop_Elements_4: do ele = 1, totele
            do FACE_COUNT=FINELE(ELE),FINELE(ELE+1)-1
             
               ele2=COLELE(FACE_COUNT)
               if((ele2>0).and.(ele2.ne.ele)) then
                  Loop_CVILOC_5: do cv_iloc = 1, cv_nloc ! Loop over nodes of the elements
                     cv_nodi = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )
                     x_nodi = x_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )
                     Loop_CVILOC_6: do cv_iloc2 = 1, cv_nloc ! Loop over nodes of the elements
                        cv_nodi2 = cv_ndgln( ( ele2 - 1 ) * cv_nloc + cv_iloc2 )
                        x_nodi2 = x_ndgln( ( ele2 - 1 ) * cv_nloc + cv_iloc2 )
                        if(x_nodi==x_nodi2) then
                            do gcount = findm( cv_nodi ), findm( cv_nodi + 1 ) - 1
                               cv_nodj2 = colm( gcount )
                              if ( cv_nodi2 == cv_nodj2 ) found( gcount ) = .true.
                            end do
                        endif
                     end do Loop_CVILOC_6
                  end do Loop_CVILOC_5
               endif
            end do
         end do Loop_Elements_4
      endif

      gcount2 = 0 ! Now reducing the size of the stencil
      do cv_nodi = 1, cv_nonods
         finacv_loc( cv_nodi ) = gcount2 + 1
         do gcount = findm( cv_nodi ), findm( cv_nodi + 1 ) - 1
            if( found( gcount ) ) then
               gcount2 = gcount2 + 1
               colacv_loc( gcount2 ) = colm( gcount )
               if( colacv_loc( gcount2 ) == cv_nodi ) midacv_loc( cv_nodi ) = gcount2
            end if
         end do
      end do
      ncolacv_loc = gcount2
      finacv_loc( cv_nonods + 1 ) = gcount2 + 1


      !ewrite(3,*) 'acv:'
      !do cv_nodi=1,cv_nonods
      !  ewrite(3,*) 'for row:',cv_nodi,' the colns are:'
      !  ewrite(3,*) (colacv_loc(count),count=finacv_loc(cv_nodi),finacv_loc(cv_nodi+1)-1)
      !end do
       
      RETURN

!      stop 1824

      deallocate( cv_neiloc )
      deallocate( cvfem_neiloc )
      deallocate( found )
      deallocate( findgpts )
      deallocate( colgpts )
      deallocate( cv_other_loc )
      deallocate( u_other_loc )
      deallocate( mat_other_loc )
      deallocate( cv_on_face )
      deallocate( x_share )
      deallocate( cvweight )
      deallocate( cvn )
      deallocate( cvn_short )
      deallocate( cvfen )
      deallocate( cvfenlx )
      deallocate( cvfenly )
      deallocate( cvfenlz )
      deallocate( cvweight_short )
      deallocate( cvfen_short )
      deallocate( cvfenlx_short )
      deallocate( cvfenly_short )
      deallocate( cvfenlz_short )
      deallocate( ufen )
      deallocate( ufenlx )
      deallocate( ufenly )
      deallocate( ufenlz )
      deallocate( cv_neiloc )
      deallocate( cv_on_face )
      deallocate( scvfen )
      deallocate( scvfenslx )
      deallocate( scvfensly )
      deallocate( scvfeweight )
      deallocate( scvfenlx )
      deallocate( scvfenly )
      deallocate( scvfenlz )
      deallocate( sufen )
      deallocate( sufenslx )
      deallocate( sufensly )
      deallocate( sufenlx )
      deallocate( sufenly )
      deallocate( sufenlz )
      deallocate( sbufen )
      deallocate( sbufenslx )
      deallocate( sbufensly )
      deallocate( sbufenlx )
      deallocate( sbufenly )
      deallocate( sbufenlz )
      deallocate( u_on_face )
      deallocate( sbcvfen )
      deallocate( sbcvfenslx )
      deallocate( sbcvfensly )
      deallocate( sbcvfeweigh )
      deallocate( sbcvfenlx )
      deallocate( sbcvfenly )
      deallocate( sbcvfenlz )
      deallocate( sele_overlap_scale )

      return
    end subroutine CV_Neighboor_Sparsity



  end module sparsity_ND



  module spact

    use fldebug
    use sparsity_1D
    use sparsity_ND
    use shape_functions

  contains

    subroutine get_spars_pats( &
         ndim, nphase, totele, u_pha_nonods, cv_pha_nonods, &
         u_nonods, cv_nonods, x_nonods, &
         cv_ele_type, u_ele_type, &
         u_nloc, cv_nloc, x_nloc, xu_nloc, mat_nloc, &
         u_snloc, cv_snloc, x_snloc, &
         u_ndgln, cv_ndgln, x_ndgln, xu_ndgln, &
                                ! CV multi-phase eqns (e.g. vol frac, temp)
         mx_ncolacv, ncolacv, finacv, colacv, midacv, &
                                ! Force balance plus cty multi-phase eqns
         nlenmcy, mx_ncolmcy, ncolmcy, finmcy, colmcy, midmcy, &
                                ! Element connectivity
         mxnele, ncolele, midele, finele, colele, &
                                ! Force balance sparsity
         mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, &
                                ! CT sparsity - global cty eqn
         mx_nct, nct, findct, colct, &
                                ! C sparsity operating on pressure in force balance
         mx_nc, nc, findc, colc, &
                                ! pressure matrix for projection method
         mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, &
                                ! CV-FEM matrix
         mx_ncolm, ncolm, findm, colm, midm, mx_nface_p1 )
      ! Obtain the sparsity patterns of the two types of matricies for
      ! (momentum + cty) and for energy
      implicit none
      integer, intent( in ) :: ndim, nphase, totele, u_pha_nonods, &
           cv_pha_nonods,u_nonods, cv_nonods, x_nonods, cv_ele_type, &
           u_ele_type, u_nloc, cv_nloc, x_nloc, xu_nloc, mat_nloc, &
           u_snloc, cv_snloc, x_snloc, mx_nface_p1
      integer, dimension( u_nloc * totele ), intent( in ) :: u_ndgln
      integer, dimension( cv_nloc * totele ), intent( in ) :: cv_ndgln, x_ndgln
      integer, dimension( xu_nloc * totele ), intent( in ) :: xu_ndgln
      integer, intent( in ) :: mx_ncolacv
      integer, intent( inout ) :: ncolacv ! CV multiphase eqns (e.g. vol frac, temp)
      integer, dimension( cv_pha_nonods + 1 ), intent( inout ) :: finacv
      integer, dimension( mx_ncolacv ), intent( inout ) :: colacv
      integer, dimension( cv_pha_nonods ), intent( inout ) :: midacv
      integer, intent( in ) :: nlenmcy, mx_ncolmcy 
      integer, intent( inout ) :: ncolmcy ! Force balance plus cty multi-phase eqns
      integer, dimension( nlenmcy + 1 ), intent( inout ) :: finmcy
      integer, dimension( mx_ncolmcy ), intent( inout ) :: colmcy
      integer, dimension( nlenmcy ), intent( inout ) :: midmcy
      integer, intent( in ) :: mxnele
      integer, intent( inout ) :: ncolele ! Element connectivity
      integer, dimension( totele ), intent( inout ) :: midele
      integer, dimension( totele + 1 ), intent( inout ) :: finele
      integer, dimension( mxnele ), intent( inout ) :: colele
      integer, intent( in ) :: mx_ncoldgm_pha
      integer, intent( inout ) :: ncoldgm_pha ! Force balance sparsity
      integer, dimension( mx_ncoldgm_pha ), intent( inout ) :: coldgm_pha
      integer, dimension( u_pha_nonods + 1 ), intent( inout ) :: findgm_pha
      integer, dimension( u_pha_nonods ), intent( inout ) :: middgm_pha
      integer, intent( in ) :: mx_nct
      integer, intent( inout ) :: nct ! CT-Sparsity - global cty eqn
      integer, dimension( cv_nonods + 1 ), intent( inout ) :: findct
      integer, dimension( mx_nct ), intent( inout ) :: colct
      integer, intent( in ) :: mx_nc
      integer, intent( inout ) :: nc ! C-Sparsity operating on pressure in force balance
      integer, dimension( u_nonods + 1 ), intent( inout ) :: findc
      integer, dimension( mx_nc ), intent( inout ) :: colc
      integer, intent( in ) :: mx_ncolcmc
      integer, intent( inout ) :: ncolcmc ! Pressure matrix for projection method
      integer, dimension( cv_nonods + 1 ), intent( inout ) :: findcmc
      integer, dimension( mx_ncolcmc ), intent( inout ) :: colcmc
      integer, dimension( cv_nonods ), intent( inout ) :: midcmc
      integer, intent( in ) :: mx_ncolm
      integer, intent( inout ) :: ncolm ! CV-FEM matrix
      integer, dimension( cv_nonods + 1 ), intent( inout ) :: findm
      integer, dimension( mx_ncolm ), intent( inout ) :: colm
      integer, dimension( cv_nonods ), intent( inout ) :: midm
      ! Local variables 
      integer :: mx_ncolele_pha, nacv_loc, nacv_loc2, ele, iloc1, iloc2, globi, globj
      integer :: mx_ncolacv_loc,count,cv_inod
      logical :: presym
      integer, dimension( : ), allocatable :: colele_pha, finele_pha, midele_pha, &
           centct, dummyvec, midacv_loc, finacv_loc, colacv_loc

      ewrite(3,*) 'In Get_Spars_Pats'

      !-
      !- Computing sparsity for element connectivity
      !-
      finele = 0 ; colele = 0 ; midele = 0
      Conditional_Dimensional_1: if ( ndim == 1 ) then 
         call def_spar( 1, totele, mxnele, ncolele, &
              midele, finele, colele )
      else
         call getfinele( totele, cv_nloc, cv_snloc, x_nonods, x_ndgln, mx_nface_p1, &
              mxnele, ncolele, finele, colele, midele )
      end if Conditional_Dimensional_1
      ewrite(3,*)'finele: ', size( finele ), '==>', finele( 1 : totele + 1 )
      ewrite(3,*)'colele: ', size( colele ), ncolele, '==>', colele( 1 : ncolele )
      ewrite(3,*)'midele: ', size( midele ), '==>', midele( 1 : totele )

      !-
      !- Computing sparsity for force balance
      !-
      mx_ncolele_pha = nphase * ncolele + ( nphase - 1 ) * nphase * totele
      allocate( colele_pha( mx_ncolele_pha ) )
      allocate( finele_pha( totele * nphase + 1 ) )
      allocate( midele_pha( totele * nphase ) )
      colele_pha = 0 ; finele_pha = 0 ; midele_pha = 0
      call exten_sparse_multi_phase( totele, ncolele, finele, colele, &
           nphase, totele * nphase, mx_ncolele_pha, &
           finele_pha, colele_pha, midele_pha )
      ewrite(3,*)'finele_pha: ', finele_pha( 1 : totele * nphase + 1 )
      ewrite(3,*)'colele_pha: ', colele_pha( 1 : mx_ncolele_pha )
      ewrite(3,*)'midele_pha: ', midele_pha( 1 : totele * nphase )

      findgm_pha = 0 ; coldgm_pha = 0 ; middgm_pha = 0
      call form_dgm_pha_sparsity( totele, nphase, u_nloc, u_pha_nonods, &
           ndim, mx_ncoldgm_pha, ncoldgm_pha, &
           coldgm_pha, findgm_pha, middgm_pha, &
           mx_ncolele_pha, finele_pha, colele_pha, finele, colele, ncolele )
      ewrite(3,*)'findgm_pha: ', findgm_pha( 1 : u_pha_nonods + 1 )
      ewrite(3,*)'coldgm_pha: ', coldgm_pha( 1 : ncoldgm_pha )
      ewrite(3,*)'middgm_pha: ', middgm_pha( 1 : u_pha_nonods )

      !-
      !- Now form the global matrix: FINMCY, COLMCY and MIDMCY for the 
      !- momentum and continuity eqns
      !-
      allocate( centct( cv_nonods ) )
      findct = 0 ; colct = 0 ; centct = 0
      Conditional_Dimensional_2: if ( ndim == 1 ) then 
         call def_spar_ct_dg( cv_nonods, mx_nct, nct, findct, colct, &
              totele, cv_nloc, u_nloc, u_ndgln, u_ele_type, cv_ndgln )
      else
        if(cv_nonods==x_nonods) then ! a continuouse pressure mesh...
         call pousinmc2( totele, u_nonods, u_nloc, cv_nonods, cv_nloc, &
              mx_nct, u_ndgln, cv_ndgln, &
              nct, findct, colct, centct )
        else ! use finele to determine COLCT
         call CT_DG_Sparsity( mx_nface_p1,  &
         totele, cv_nloc, u_nloc,  &
         cv_nonods, &
         cv_ndgln,u_ndgln,  &
         ncolele, finele, colele, &
         mx_nct, nct, findct, colct )
        endif

         ewrite(3,*),'u_nonods, u_nloc, cv_nonods, cv_nloc, mx_nct, nct:', &
              u_nonods, u_nloc, cv_nonods, cv_nloc, mx_nct, nct

      end if Conditional_Dimensional_2
      nc = nct
      ewrite(3,*) 'findct: ', size( findct ), '==>', findct( 1 : cv_nonods + 1 )
      ewrite(3,*) 'colct: ', size( colct ), nct, '==>', colct( 1 : nct )

      !-
      !- Convert CT sparsity to C sparsity
      !-
      call conv_ct2c( cv_nonods, nct, findct, colct, u_nonods, &
           mx_nc, findc, colc )
      ewrite(3,*) 'findc: ', size( findc ), '==>', findc( 1 : u_nonods + 1 )
      ewrite(3,*) 'colc: ', size( findc ), nc, '==>', colc( 1 : nc )
!      stop 666

      !-
      !- Computing sparsity for pressure matrix of projection method
      !-
      Conditional_Dimensional_3: if ( ndim == 1 ) then 
         Conditional_ContinuousPressure_3: if ( cv_nonods /= totele * cv_nloc ) then
            call def_spar( cv_nloc - 1, cv_nonods, mx_ncolcmc, ncolcmc, &
                 midcmc, findcmc, colcmc )
         else ! Discontinuous pressure mesh
            call def_spar( cv_nloc + 2, cv_nonods, mx_ncolcmc, ncolcmc, &
                 midcmc, findcmc, colcmc )
         end if Conditional_ContinuousPressure_3
      else
         allocate( dummyvec( u_nonods ))
         dummyvec = 0
         presym = .false.
         call poscmc( cv_nonods, u_nonods, mx_ncolcmc, nct, &
              findct, colct, &
              ncolcmc, findcmc, colcmc, midcmc, dummyvec, presym )
         deallocate( dummyvec )
      end if Conditional_Dimensional_3
      if( mx_ncolcmc < ncolcmc ) FLAbort("Incorrect number of dimension of CMC sparsity matrix")
      ewrite(3,*)'findcmc: ', size( findcmc ), '==>', findcmc( 1 : cv_nonods + 1 )
      ewrite(3,*)'colcmc: ', size( colcmc ), ncolcmc, '==>', colcmc( 1 : ncolcmc )
      ewrite(3,*)'midcmc: ', size( midcmc ), '==>',  midcmc( 1 : cv_nonods )

      !-
      !- Computing the sparsity for the force balance plus cty multi-phase eqns
      !- 
      finmcy = 0 ; colmcy = 0 ; midmcy = 0
      call exten_sparse_mom_cty( ndim, findgm_pha, coldgm_pha, u_pha_nonods, ncoldgm_pha, nct, &
           cv_nonods, findct, colct, &
           u_nonods, nc, mx_ncolmcy, &
           findc, colc, finmcy, colmcy, midmcy, nlenmcy, &
           ncolmcy, nphase, ncolcmc, findcmc, colcmc )
      ewrite(3,*)'finmcy: ', size( finmcy ), nlenmcy + 1, '==>', finmcy( 1 : nlenmcy + 1 )
      ewrite(3,*)'colmcy: ', size( colmcy ), ncolmcy, '==>', colmcy( 1 : ncolmcy )
      ewrite(3,*)'midmcy: ', size( midmcy ), nlenmcy, '==>', midmcy( 1 : nlenmcy )

      !-
      !- Computing sparsity CV-FEM
      !-
      findm = 0 ; colm = 0 ; midm = 0
      Conditional_Dimensional_4: if ( ndim == 1 ) then 
         call def_spar( cv_nloc - 1, cv_nonods, mx_ncolm, ncolm, &
              midm, findm, colm )
      else
        if(cv_nonods==x_nonods) then ! a continuous pressure mesh
         call pousinmc2( totele, cv_nonods, cv_nloc, cv_nonods, cv_nloc, mx_ncolm, cv_ndgln, cv_ndgln, &
              ncolm, findm, colm, midm )
        else ! a DG pressure field mesh...
         call CT_DG_Sparsity( mx_nface_p1, &
         totele, cv_nloc, cv_nloc, &
         cv_nonods, &
         cv_ndgln, cv_ndgln, &
         ncolele, finele, colele, &
         mx_ncolm, ncolm, findm, colm )
         ! determine midm...
         do cv_inod = 1, cv_nonods
            do count = findm( cv_inod ), findm( cv_inod + 1 ) - 1
               if( colm( count ) == cv_inod ) midm( cv_inod ) = count
            end do
         end do
        end if
      end if Conditional_Dimensional_4
      ewrite(3,*)'findm: ', size( findm ), '==>', findm( 1 : cv_nonods + 1 )
      ewrite(3,*)'colm: ', size( colm ), ncolm, '==>', colm( 1 : ncolm )
      ewrite(3,*)'midm: ', size( midm ), '==>', midm( 1 : cv_nonods )

      !-
      !- Computing sparsity for CV multiphase eqns (e.g. vol frac, temp)
      !-
      mx_ncolacv_loc=mx_ncolacv/nphase
      allocate( midacv_loc( cv_nonods ) )
      allocate( finacv_loc( cv_nonods + 1 ) )
      allocate( colacv_loc( mx_ncolacv_loc ) )
      midacv_loc = 0 ; finacv_loc = 0 ; colacv_loc = 0
      Conditional_Dimensional_5: if ( ndim == 1 ) then 
         call def_spar( 1, cv_nonods, 3 * cv_nonods, nacv_loc, &
              midacv_loc, finacv_loc, colacv_loc )
      else
         call CV_Neighboor_Sparsity( ndim, nphase, cv_ele_type, &
              totele, cv_nloc, u_nloc, x_nloc, xu_nloc, mat_nloc, &
              cv_snloc, u_snloc, cv_nonods, x_nonods, &
              cv_ndgln, x_ndgln, xu_ndgln, &
              ncolele, finele, colele, &
              ncolm, mx_ncolacv_loc, findm, colm, &
              nacv_loc, finacv_loc, colacv_loc, midacv_loc )
      end if Conditional_Dimensional_5
      nacv_loc2 = nacv_loc
      ewrite(3,*)'finacv_loc: ', size( finacv_loc ), '==>', finacv_loc( 1 : cv_nonods + 1 )
      ewrite(3,*)'colacv_loc: ', size( colacv_loc ), '==>', colacv_loc( 1 : nacv_loc2 )
      ewrite(3,*)'midacv_loc: ', size( midacv_loc ), '==>', midacv_loc( 1 : cv_nonods )

      ncolacv =  nphase * nacv_loc + ( nphase - 1 ) * nphase * cv_nonods
      nacv_loc = ncolacv
      finacv = 0 ; colacv = 0 ; midacv = 0
      !ewrite(3,*) 'ncolacv, mx_ncolacv:', ncolacv, mx_ncolacv
      call exten_sparse_multi_phase( cv_nonods, nacv_loc2, finacv_loc, colacv_loc, &
           nphase, cv_pha_nonods, ncolacv, &
           finacv, colacv, midacv )
      ewrite(3,*)'finacv: ', finacv( 1 : cv_pha_nonods + 1 )
      ewrite(3,*)'colacv: ', colacv( 1 : ncolacv )
      ewrite(1,*)'midacv: ', midacv( 1 : cv_pha_nonods )

      !-
      !- Deallocating temporary arrays
      !-
      deallocate( colele_pha )
      deallocate( finele_pha )
      deallocate( midele_pha )
      deallocate( centct )
      deallocate( midacv_loc )
      deallocate( finacv_loc )
      deallocate( colacv_loc )

      return
    end subroutine get_spars_pats




    subroutine CT_DG_Sparsity( mx_nface_p1, &
         totele, cv_nloc, u_nloc, &
         cv_nonods, &
         cv_ndgln,u_ndgln, &
         ncolele, finele, colele, &
         mx_nct, nct, findct, colct )
      implicit none
      integer, intent( in ) :: mx_nface_p1, totele, cv_nloc, &
           u_nloc, cv_nonods
      integer, dimension( totele * cv_nloc ), intent( in ) :: cv_ndgln
      integer, dimension( totele * u_nloc ), intent( in ) :: u_ndgln
      integer, intent( in ) :: ncolele
      integer, dimension( totele + 1 ), intent( in ) :: finele
      integer, dimension( ncolele ), intent( in ) :: colele
      integer, intent( in ) :: mx_nct
      integer, intent( inout ) :: nct
      integer, dimension( cv_nonods + 1 ), intent( inout ) :: findct
      integer, dimension( mx_nct ), intent( inout ) :: colct
      ! Local variables
      integer :: ele, ele2, cv_iloc, cv_nodi, count, count2, u_jloc, u_nodj, &
           gcount2, gcount, FACE_COUNT
      integer, dimension( : ), allocatable :: find_ct_temp, no_in_row

      allocate( find_ct_temp( cv_nonods + 1 ) ) ; find_ct_temp = 0
      allocate( no_in_row( cv_nonods ) ) ; no_in_row = 0

      find_ct_temp(1) = 1
      do cv_nodi = 1, cv_nonods
         find_ct_temp( cv_nodi + 1 ) = find_ct_temp( cv_nodi ) &
              + ( mx_nface_p1 + 1 ) * u_nloc
      end do

      colct = 0
      Loop_Elements_4: do ele = 1, totele
               Loop_CVILOC_5: do cv_iloc = 1, cv_nloc ! Loop over nodes of the elements
                  cv_nodi = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )
                  count2 = 0

         do FACE_COUNT = FINELE( ELE ), FINELE( ELE + 1 ) - 1

            ele2 = COLELE( FACE_COUNT )
!            print *,'ele,ele2:',ele,ele2
            if( ele2 > 0 ) then
                  Loop_CVILOC_6: do u_jloc = 1, u_nloc ! Loop over nodes of the elements
                     u_nodj = u_ndgln( ( ele2 - 1 ) * u_nloc + u_jloc )
                     count2 = count2 + 1
                     count = find_ct_temp( cv_nodi ) - 1 + count2
                     colct( count ) = u_nodj
                     no_in_row( cv_nodi ) = no_in_row( cv_nodi ) + 1
                  end do Loop_CVILOC_6
            end if
         end do

               end do Loop_CVILOC_5
      end do Loop_Elements_4

!      stop 22

      ! shrink the sparcity up a little now...
      gcount2 = 0 ! Now reducing the size of the stencil
      do cv_nodi = 1, cv_nonods
         findct( cv_nodi ) = gcount2 + 1
         do gcount = find_ct_temp( cv_nodi ), find_ct_temp( cv_nodi + 1 ) - 1
            if( colct( gcount ) /= 0 ) then
               gcount2 = gcount2 + 1
               colct( gcount2 ) = colct( gcount )
            end if
         end do
      end do
      nct = gcount2
      findct( cv_nonods + 1 ) = gcount2 + 1

      ! sort colct in increasing order
      do cv_nodi = 1, cv_nonods
         call ibubble( colct( findct( cv_nodi ) : findct( cv_nodi + 1 ) -1 ) )
         ewrite(3,*) 'cv_nodi, colct:', &
              cv_nodi, colct( findct( cv_nodi ) : findct( cv_nodi + 1 ) -1 )
      end do

      return
    end subroutine CT_DG_Sparsity

    subroutine ibubble(ivec)
      ! sort ivec in increasing order
      implicit none
      integer, dimension( : ), intent( inout ) :: ivec
      ! Local variables
      integer :: nvec, i, j, itemp

      nvec = size(ivec)

      do i = 1, nvec - 1
         do j = 1, nvec
            if ( ivec( i ) > ivec( i + 1 ) ) then
               itemp = ivec( i + 1 )
               ivec( i + 1 ) = ivec( i )
               ivec( i ) = itemp
            end if
         end do
      end do
      return
    end subroutine ibubble
	 
    subroutine check_sparsity( &
         u_pha_nonods, cv_pha_nonods, &
         u_nonods, cv_nonods, totele, &
         mx_ncolacv, ncolacv, finacv, colacv, midacv, & ! CV multi-phase eqns (e.g. vol frac, temp)
         nlenmcy, mx_ncolmcy, ncolmcy, finmcy, colmcy, midmcy, & ! Force balance plus cty multi-phase eqns
         mxnele, ncolele, midele, finele, colele, & ! Element connectivity 
         mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, & ! Force balance sparsity  
         mx_nct, ncolct, findct, colct, & ! CT sparsity - global cty eqn
         mx_nc, ncolc, findc, colc, & ! C sparsity operating on pressure in force balance
         mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, & ! pressure matrix for projection method
         mx_ncolm, ncolm, findm, colm, midm )

      implicit none
      integer, intent( in ) :: u_pha_nonods, cv_pha_nonods, u_nonods, cv_nonods, totele
      integer, intent ( in ) :: mx_ncolacv, ncolacv
      integer, dimension( cv_pha_nonods + 1 ), intent (in ) :: finacv
      integer, dimension( mx_ncolacv ), intent (in ) :: colacv
      integer, dimension( cv_pha_nonods ), intent (in ) :: midacv
      integer, intent ( in ) :: nlenmcy, mx_ncolmcy, ncolmcy
      integer, dimension( nlenmcy + 1 ), intent (in ) :: finmcy
      integer, dimension( mx_ncolmcy ), intent (in ) :: colmcy
      integer, dimension( nlenmcy ), intent (in ) :: midmcy
      integer, intent ( in ) :: mxnele, ncolele
      integer, dimension( totele ), intent (in ) :: midele
      integer, dimension( totele + 1 ), intent (in ) :: finele
      integer, dimension( mxnele ), intent (in ) :: colele
      integer, intent ( in ) :: mx_ncoldgm_pha, ncoldgm_pha
      integer, dimension( mx_ncoldgm_pha ), intent (in ) :: coldgm_pha
      integer, dimension( u_pha_nonods + 1 ), intent (in ) :: findgm_pha
      integer, dimension( u_pha_nonods ), intent (in ) :: middgm_pha
      integer, intent ( in ) :: mx_nct, ncolct
      integer, dimension( cv_nonods + 1 ), intent (in ) :: findct
      integer, dimension( mx_nct ), intent (in ) :: colct
      integer, intent ( in ) :: mx_nc, ncolc
      integer, dimension( u_nonods + 1 ), intent (in ) :: findc
      integer, dimension( mx_nc ), intent (in ) :: colc
      integer, intent ( in ) :: mx_ncolcmc, ncolcmc
      integer, dimension( cv_nonods + 1 ), intent (in ) :: findcmc
      integer, dimension( mx_ncolcmc ), intent (in ) :: colcmc
      integer, dimension( cv_nonods ), intent (in ) :: midcmc
      integer, intent ( in ) :: mx_ncolm, ncolm
      integer, dimension( cv_nonods + 1 ), intent (in ) :: findm
      integer, dimension( ncolm ), intent (in ) :: colm
      integer, dimension( cv_nonods ), intent (in ) :: midm

      ! Local variables
      integer, dimension( : ), allocatable :: dummy

      ewrite(3,*) 'In check_sparsity'

      open( 15, file = 'CheckSparsityMatrix.dat', status = 'unknown' )
      write( 15, * )'########## FINMCY, MIDMCY, COLMCY ##################'
      write(15, * )'NCOLMCY:', NCOLMCY
      call checksparsity( .true., 15, NCOLMCY, NLENMCY, MX_NCOLMCY, FINMCY, MIDMCY, COLMCY )

      write( 15, * )'########## FINACV, COLACV, MIDACV ##################'
      write(15, * )'NCOLACV:', NCOLACV
      call checksparsity( .true., 15, NCOLACV, CV_PHA_NONODS, MX_NCOLACV, FINACV, MIDACV, COLACV  )

      write( 15, * )'########## FINELE, MIDELE, COLELE  ##################'
      write(15, * )'NCOLELE:',NCOLELE 
      call checksparsity( .true., 15, NCOLELE, TOTELE, MXNELE, FINELE, MIDELE, COLELE )

      allocate( dummy( CV_NONODS ))
      write( 15, * )'########## FINDCT, COLCT ##################'
      write(15, * )'NCOLCT:', NCOLCT
      call checksparsity( .false., 15, NCOLCT, CV_NONODS, MX_NCT, FINDCT, dummy, COLCT  )
      deallocate( dummy )

      allocate( dummy( U_NONODS ))
      write( 15, * )'########## FINDC, COLC ##################'
      write(15, * )'NCOLC:', NCOLC
      call checksparsity( .false., 15, NCOLC, U_NONODS, MX_NC, FINDC, dummy, COLC )
      deallocate( dummy )

      write( 15, * )'########## FINDGM_PHA, MIDDGM_PHA, COLDGM_PHA ##################'
      write(15, * )'NCOLDGM_PHA:',NCOLDGM_PHA 
      call checksparsity( .true., 15, NCOLDGM_PHA, U_PHA_NONODS, MX_NCOLDGM_PHA, FINDGM_PHA, MIDDGM_PHA, COLDGM_PHA )

      write( 15, * )'########## FINDCMC, MIDCMC, COLCMC ##################'
      write(15, * )'NCOLCMC:',NCOLCMC 
      call checksparsity( .true., 15, NCOLCMC, CV_NONODS, MX_NCOLCMC, FINDCMC, MIDCMC, COLCMC )

      write( 15, * )'########## FINDM, MIDM, COLM ##################'
      write(15, * )'NCOLM:',NCOLM 
      call checksparsity( .true., 15, NCOLM, CV_NONODS, MX_NCOLM, FINDM, MIDM, COLM )

      close( 15 )

      ewrite(3,*) 'Leaving check_sparsity'

      return

    end subroutine check_sparsity

  subroutine checksparsity( option_mid, unit, ncol2, nonods, ncol, find, mid, col )
    implicit none
    logical, intent( in ) :: option_mid
    integer, intent( in ) :: unit, ncol2, nonods, ncol
    integer, dimension( nonods + 1 ), intent( in ) :: find
    integer, dimension( nonods ), intent( in ) :: mid
    integer, dimension( ncol ), intent( in ) :: col
    ! Local variables
    integer :: inod, icol

    ewrite(3,*) 'In checksparsity'

    write( unit, * )'find:', ( find( inod ), inod = 1, nonods + 1 )
    if( option_mid )write( unit, * )'mid:', ( mid( inod ), inod = 1, nonods )
    write( unit, * )'col:', ( col( icol ), icol = 1, ncol2  )

    ewrite(3,*) 'Leaving checksparsity'

    return
  end subroutine checksparsity

  end module spact


