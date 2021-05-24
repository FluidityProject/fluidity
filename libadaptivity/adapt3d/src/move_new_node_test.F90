#include "ewrite.h"
subroutine move_new_node_test( &
      ! common block stuff:
      SIZBIG, MAXBIG, SIZNOD, MXNODS, MAXEDE, IERR, &
      ! same arguments as in/out of edgtst:
      BIGLST, NODLST, &
      ENLBAS, ENLIST, NELBAS, NELIST, EELIST, &
      SZENLS, SZNELS, NNOD,   NELM, &
      NODX,   NODY,   NODZ,   ORGMTX, &
      NDPTRS, &
      NUMEDE_around_edge, SURFAC, INTRNL, &
      FNCDIF, FNCORG, &
      AVEDIF, AVEORG, &
      XN, YN, ZN, BSTELM, &
      XA, YA, ZA, BSTELA, &
      ALWAYS, DOTOP, MINCHG, &
      ! additional local variables from edgtst:
      INOD, NOD1, NOD2, IEDG)
      ! This routine is called immediately after considering an edge split
      ! in edgtst by adding a new node in the middle of the edge - the
      ! connecting new edges and new elements have not been created
      ! If the new node does not immediately improve the functional, this
      ! routine is called to see if we move the node whether improvement
      ! can be obtained. The node movement algorithm is the same as in nodtst
      ! which is applied to evaluate node movement of existing nodes.
      ! If a location has found that does improve, FNCDIF or AVEDIF (depending 
      ! on whether the worst element or the average over involved elements
      ! has improved) are updated to the new difference (should be neg. for improvement)
      ! wrt to the original unsplit configuration (passed in through FNCORG and AVEORG)
      ! and XN, YN, ZN, BSTELM (resp. XA, YA, ZA, BSTELA) return the new location
      ! and the element in the "original" unadapted mesh (for metric interpolation)

      use write_log
      implicit none

      ! variables from include files (which are in fixed format so we can't use them here)
      integer, intent(in) :: SIZBIG, MAXBIG, SIZNOD, MXNODS, MAXEDE
      integer, intent(inout) :: IERR
      ! lists for edge and element storage and node storage (see include/blkprm.i)
      INTEGER BIGLST(SIZBIG,MAXBIG)
      REAL NODLST(SIZNOD,MXNODS)
      ! original unadapted mesh and metric (for metric interpolation):
      integer SZENLS, SZNELS, NNOD, NELM
      INTEGER ENLIST(SZENLS), NELIST(SZNELS), &
        ENLBAS(NELM+1), NELBAS(NNOD+1), EELIST(NELM*4)
      REAL NODX(NNOD), NODY(NNOD), NODZ(NNOD), ORGMTX(9,NNOD)


      ! nodes around edge (node shared by upper and lower edge)
      ! NOTE: this has an additional entry (NUMEDE_around_edge+1) if around surface edge
      ! in which case 1:NUMEDE_around_edge-1 contain the internal nodes in order,
      ! and the last 2 nodes are surface nodes, with the penultimate being the surface node connected
      ! to the last internal node, and the last entry is the surface node connected to the first inetrnal node
      ! Internally we add NOD1 and NOD2, so we can reuse this as all nodes connected to the new node INOD
      integer, intent(inout) :: NDPTRS(MAXEDE+1) ! nodes around edge
      
      ! NUMEDE_around_edge is NUMEDE in edgtst, the number of opposite edges (and thus elements) around the edge
      ! but in this routine NUMEDE will be the n/o edges connected to the new node
      integer  NUMEDE_around_edge
      ! passed on from edgels, surfac>0 for surface edge, surfac==0 for internal edge
      ! surfac<0 signals a region id transition in the elements around the edge
      integer SURFAC
      ! internal edge I think? - not used here at the moment
      integer INTRNL
      ! these are called BSTINE and BSTINA in edgtst, but we use their nodtst names here
      ! element in the unadapted mesh of new location
      integer BSTELM, BSTELA

      REAL FNCDIF, FNCORG, AVEDIF, AVEORG, &
          ALWAYS, DOTOP, MINCHG

      ! called BSTX, BSTY, BSTZ, BSTXA, BSTYA, BSTZA in edgtst:
      real XN, YN, ZN
      real XA, YA, ZA

      ! additional arguments:
      integer, intent(in) :: INOD ! node number that's just been added
      integer, intent(in) :: NOD1, NOD2 ! nodes of split edge
      integer, intent(in) :: IEDG ! edge to split

      ! variables that are the same as in nodtst:
      INTEGER I, J, K, &
              NUMEDE, NUMELM, INELM, JNOD, &
              ELNDS(4,MAXEDE), ORGELM, ifail, ii, jj, ki, kj, &
              NTST, ITST, cnt

      REAL F1, F2, FT, T1, &
           X(MAXEDE), Y(MAXEDE), Z(MAXEDE), sx, sy, sz, fadd, &
           XO, YO, ZO, V(MAXEDE), &
           XX(10), YY(10), ZZ(10), VOL2, &
           EDGFNC, ELMFNC, RADS(MAXEDE), ORGM(9), A(9), &
           r1, r2, r3, r4, r5, r6, &
           x1, y1, z1, x2, y2, z2, dd(maxede), xp(4), yp(4), zp(4), &
           vol(maxede), tetvol, det

      LOGICAL ISINSD, ISSPLT, GETFLG

      ! variables added
      integer num_sub_nodes, sub_nodes(MAXEDE)

      ! this routine is mostly copied from notst but adapted because the new
      ! edges and elements resulting from the split have not been created
      ! the new node has been added to nodlst though
      ! just like nodtst, if it finds a better position for this node (INOD)
      ! it updates FNCDIF or AVEDIF (depending on whether it improves the worst case or the average)

      ! three cases to consider:
      ! 1) ISSPLT: splitting a geometry edge, where the node is only allowed to move in 1D along the split edge
      ! 2) ISSPLT and .not. ISINSD: splitting an edge on the surface, where the node is only allowed to move over the two adjacent
      ! (and presumable coplanar) surface triangles
      ! 3) ISINSD: internal node that can move in all directions
      ! unlike nodtst, we obtain this from the split edge, as the new nodes does not have its flags set
      ISSPLT = GETFLG( BIGLST, NODLST, IEDG, 2 )
      ISINSD = GETFLG( BIGLST, NODLST, IEDG, 3 ) .AND. .NOT. ISSPLT

      ! the criterion for having an additional (surface) node around edge at NDPTRS(NUMEDE_around_edge+1)
      ! in edglst is SURFAC>0 we expect this extra node if .not. ISINSD
      if (ISINSD .neqv. (SURFAC<=0)) then
        ewrite(-1,*) "Mix up of external/internal edge in move_new_node_test:"
        ewrite(-1,*) "ISSPLT, ISINSD, SURFAC", ISSPLT, ISINSD, SURFAC
        ierr = -1
        return
      end if

      ! in this first bit we set up NUMEDE, NDPTRS and ELNDS as nodtst expects them
      ! instead of calling NDCNCT we fill these in from the known nodes in edgsplt
      ! we do not use ELPTRS (containting the elements connected to the new node)
      ! as these do not exist yet, instead we immediately fill in ELNDS (the nodes in these elements)

      ! n/o edges/nodes around new node:
      NUMEDE = NUMEDE_around_edge + 2
      if (.not. ISINSD) NUMEDE = NUMEDE + 1
      ! NDPTRS is already filled 1:NUMEDE_around_edge (+1 if surface edge)
      ! we simply add two more in unused space
      NDPTRS(NUMEDE-1) = NOD1
      NDPTRS(NUMEDE) = NOD2

      ! work out between which nodes we are allowed to move
      ! (these nodes span a 1d, 2d, or 3d subspace)
      if (ISINSD) then
        ! free to move between all connected nodes:
        num_sub_nodes = NUMEDE
        sub_nodes(1:NUMEDE) = NDPTRS(1:NUMEDE)
      elseif (ISSPLT) then
        ! only move between NOD1 and NOD2
        num_sub_nodes = 2
        sub_nodes(1:2) = NDPTRS(NUMEDE-1:NUMEDE)
      else
        ! move between NOD1, NOD2, and the other two surface nodes which are at
        ! the end of NDPTRS
        num_sub_nodes = 4
        sub_nodes(1:4) = NDPTRS(NUMEDE-3:NUMEDE)
      end if

      ! nodes of prospective new elements
      ! which consists of pairs containing either of 
      ! the new (nod1, inod) or (inod, no2) edges, and both have
      ! the same opposite edges, as the original opposite edges around the original split edge
      NUMELM = 2*NUMEDE_around_edge
      do i=1, NUMEDE_around_edge
        ELNDS(1,i*2-1) = NOD1
        ELNDS(2,i*2-1) = INOD
        ELNDS(3,i*2-1) = NDPTRS(i)
        ELNDS(4,i*2-1) = NDPTRS(i+1)
        ELNDS(1,i*2) = INOD
        ELNDS(2,i*2) = NOD2
        ELNDS(3,i*2) = NDPTRS(i)
        ELNDS(4,i*2) = NDPTRS(i+1)
      end do
      if (ISINSD) then
        ! for the final opposite edge we go back around:
        ! it is between the last node (NUMEDE_around_edge) and the first node
        ELNDS(4,NUMEDE_around_edge*2-1) = NDPTRS(1)
        ELNDS(4,NUMEDE_around_edge*2) = NDPTRS(1)
      else
        ELNDS(3,NUMEDE_around_edge*2-1) = NDPTRS(1)
        ELNDS(3,NUMEDE_around_edge*2) = NDPTRS(1)
      end if

      !  start copy from nodtst

      ! NOTE: unlike nodtst, we do not initialize FNCORG, AVEORG, BSTELM, BSTELA,
      ! FNCDIF and AVEDIF as we retain their values from edgtst (BSTELM=BSTINE, BSTELA=BSTINA in edgtst)
      ! only if we improve over the split configuration with node movement (which wasn't a sufficient improvement
      ! over the original unsplit configuration) do we update FNCORG+BSTELM (if the worst case improves)
      ! and AVEORG and BSTELA (if the average improves)


      ! original location, element (in unadapted mesh) and metric in new node
      XO = NODLST( 1, INOD )
      YO = NODLST( 2, INOD )
      ZO = NODLST( 3, INOD )
      DO I = 1, 9
         ORGM(I) = NODLST(I+6,INOD)
      END DO
      ORGELM = INT(NODLST(16,INOD))
      INELM  = ORGELM

      XN = XO
      YN = YO
      ZN = ZO
      XA = XO
      YA = YO
      ZA = ZO

      XX(1) = 0.0
      YY(1) = 0.0
      ZZ(1) = 0.0


      I = 0

      a = 0.0

      sx = 0.0
      sy = 0.0
      sz = 0.0

      DO J = 1, num_sub_nodes
         jnod = sub_nodes(j)
         X(J) = NODLST(1,jnod) - XO
         Y(J) = NODLST(2,jnod) - YO
         Z(J) = NODLST(3,jnod) - ZO
         t1 = edgfnc( BIGLST, NODLST, jnod, inod, v(j) )
         do k = 1, 9
            a(k) = a(k) + nodlst(k+6,jnod)
         end do
         sx = sx + nodlst( 7,jnod)*x(j) &
                 + nodlst(10,jnod)*y(j) &
                 + nodlst(13,jnod)*z(j)
         sy = sy + nodlst( 8,jnod)*x(j) &
                 + nodlst(11,jnod)*y(j) &
                 + nodlst(14,jnod)*z(j)
         sz = sz + nodlst( 9,jnod)*x(j) &
                 + nodlst(12,jnod)*y(j) &
                 + nodlst(15,jnod)*z(j)
      END DO

      ! not sure what this is about - making sure a is diag. dominant?
      fadd = 0.0
      fadd = abs(a(2)) + abs(a(3))
      a(1) = max(a(1),fadd*1.01)
      fadd = abs(a(4)) + abs(a(6))
      a(5) = max(a(5),fadd*1.01)
      fadd = abs(a(7)) + abs(a(8))
      a(9) = max(a(9),fadd*1.01)


      do j = 1, num_sub_nodes
         sx = sx + nodlst( 7,inod)*x(j) &
                 + nodlst(10,inod)*y(j) &
                 + nodlst(13,inod)*z(j)
         sy = sy + nodlst( 8,inod)*x(j) &
                 + nodlst(11,inod)*y(j) &
                 + nodlst(14,inod)*z(j)
         sz = sz + nodlst( 9,inod)*x(j) &
                 + nodlst(12,inod)*y(j) &
                 + nodlst(15,inod)*z(j)
      end do

      if (ISINSD) then
         call invrs3( a, det, ifail )
       else if (ISSPLT) then
         ! ok...?
         r1 = a(1)*x(1) + a(4)*y(1) + a(7)*z(1)
         r2 = a(2)*x(1) + a(5)*y(1) + a(8)*z(1)
         r3 = a(3)*x(1) + a(6)*y(1) + a(9)*z(1)
         r4 = x(1)*r1 + y(1)*r2 + z(1)*r3
         r1 = x(1)*sx + y(1)*sy + z(1)*sz
       else
         do ii = 1, num_sub_nodes
            dd(ii) = 1.0/sqrt( x(ii)*x(ii) + y(ii)*y(ii) + z(ii)*z(ii) )
         end do

         ii = 0

         fadd = 1.0

  10     if( ii .lt. num_sub_nodes-1 ) then
            ii = ii + 1
            r1 = dd(ii)
            jj = ii
  20        if( jj .lt. num_sub_nodes ) then
               jj = jj + 1
               r2 = dd(jj)
               r3 = x(ii)*x(jj) + y(ii)*y(jj) + z(ii)*z(jj)
               r3 = r3*r1*r2
               if( abs(r3) .lt. abs(fadd) ) then
                  fadd = r3
                  ki   = ii
                  kj   = jj
                  if( abs(fadd) .lt. 0.6 ) goto 30
               end if
               goto 20
            end if
            goto 10
         end if

  30     continue

         x1 = x(ki)*dd(ki)
         y1 = y(ki)*dd(ki)
         z1 = z(ki)*dd(ki)

         x2 = x(kj)*dd(kj) - fadd*x1
         y2 = y(kj)*dd(kj) - fadd*y1
         z2 = z(kj)*dd(kj) - fadd*z1

         r1 = 1.0/sqrt(x2*x2 + y2*y2 + z2*z2)

         x2 = x2*r1
         y2 = y2*r1
         z2 = z2*r1

         r1 = a(1)*x1 + a(4)*y1 + a(7)*z1
         r2 = a(2)*x1 + a(5)*y1 + a(8)*z1
         r3 = a(3)*x1 + a(6)*y1 + a(9)*z1
         r4 = r1*x1 + r2*y1 + r3*z1
         r5 = r1*x2 + r2*y2 + r3*z2
         r1 = a(1)*x2 + a(4)*y2 + a(7)*z2
         r2 = a(2)*x2 + a(5)*y2 + a(8)*z2
         r3 = a(3)*x2 + a(6)*y2 + a(9)*z2
         r6 = r1*x2 + r2*y2 + r3*z2

         a(1) = r4
         a(2) = r5
         a(3) = r5
         a(4) = r6
         r1 = x1*sx + y1*sy + z1*sz
         r2 = x2*sx + y2*sy + z2*sz

         call invrs2( a, det, ifail )

      end if

      ! new position apparently
      if(ISINSD) then
         xx(1) = (sx*a(1) + sy*a(4) + sz*a(7))*1.5 + xo
         yy(1) = (sx*a(2) + sy*a(5) + sz*a(8))*1.5 + yo
         zz(1) = (sx*a(3) + sy*a(6) + sz*a(9))*1.5 + zo
      else if(ISSPLT) then
         r1 = r1/r4
         xx(1) = x(1)*r1*1.5 + xo
         yy(1) = y(1)*r1*1.5 + yo
         zz(1) = z(1)*r1*1.5 + zo
      else
         r3 = r1*a(1) + r2*a(3)
         r4 = r1*a(2) + r2*a(4)
         xx(1) = (x1*r3 + x2*r4)*1.5 + xo
         yy(1) = (y1*r3 + y2*r4)*1.5 + yo
         zz(1) = (z1*r3 + z2*r4)*1.5 + zo
      end if

      NTST = 1

      fadd = 1e+20

      DO I = 1, NUMELM
         FT = ELMFNC( BIGLST, NODLST, 0, ELNDS(1,I), ELNDS(2,I), &
                                ELNDS(3,I), ELNDS(4,I), RADS(I) )
         do j = 1, 4
            xp(j) = nodlst(1,elnds(j,i))
            yp(j) = nodlst(2,elnds(j,i))
            zp(j) = nodlst(3,elnds(j,i))
         end do
         vol(i) = tetvol( xp, yp, zp )
      END DO

      DO 500 ITST = 1, NTST
        cnt = 1

 100    xx(itst) = xx(itst)*0.6 + xo*0.4
        yy(itst) = yy(itst)*0.6 + yo*0.4
        zz(itst) = zz(itst)*0.6 + zo*0.4

        nodlst(1,inod) = xx(ITST)
        nodlst(2,inod) = yy(ITST)
        nodlst(3,inod) = zz(ITST)

        do i = 1, numelm

          do j = 1, 4
            xp(j) = nodlst(1,elnds(j,i))
            yp(j) = nodlst(2,elnds(j,i))
            zp(j) = nodlst(3,elnds(j,i))
          end do

          if( tetvol( xp, yp, zp )/vol(i) .le. 1e-4 ) then
            cnt = cnt + 1
            if( cnt .lt. 6 ) goto 100
            goto 500
          end if

        end do

        inelm = abs(inelm)

        CALL FNDELM( ENLBAS, ENLIST, NELBAS, NELIST, EELIST, &
                     SZENLS, SZNELS, NNOD,   NELM, &
                     NODX,   NODY,   NODZ,   ORGMTX, &
                     XX(ITST),  YY(ITST),  ZZ(ITST), &
                     INELM,  .TRUE., NODLST(7,INOD) )

        if( inelm .lt. 0 ) then
           print*,'---+++ NODTST: Got node outside mesh +++---'
           if( issplt ) then
              print*,'   A Splitter node'
           else if( .not. isinsd ) then
              print*,'   A surface node'
           else
              print*,'   An internal node'
           end if
           print*,xo,yo,zo
           print*,xx(itst),yy(itst),zz(itst)
        end if

        nodlst(1,inod)  = xx(ITST)
        nodlst(2,inod)  = yy(ITST)
        nodlst(3,inod)  = zz(ITST)
        NODLST(16,INOD) = FLOAT(ABS(INELM))

        F1 = 0.0
        F2 = 0.0

        DO I = 1, NUMELM

           FT = ELMFNC( BIGLST, NODLST, 0, ELNDS(1,I), ELNDS(2,I), &
                       ELNDS(3,I), ELNDS(4,I), VOL2 )
           F2 = F2 + FT
           F1 = MAX( F1, FT )
           IF( VOL2/RADS(I) .LE. 1E-4 ) then
              cnt = cnt + 1
              if( cnt .lt. 6 ) goto 100
              GOTO 500
           end if

        END DO

        F1 = F1 - FNCORG
        F2 = F2 / NUMELM - AVEORG

        IF( F1 .LT. FNCDIF ) THEN
           FNCDIF = F1
           XN = XX(ITST)
           YN = YY(ITST)
           ZN = ZZ(ITST)
           BSTELM = INELM
           if( INELM .lt. 0 .and. fncdif .lt. 0.0 ) then
              print*,'*** NODTST: Node outside mesh is better!'
              print*,fncorg,aveorg,fncdif,min(f2,avedif)
           end if
        END IF

        IF( F2 .LT. AVEDIF .AND. F1 .LT. 0.0 ) THEN
           AVEDIF = F2
           XA = XX(ITST)
           YA = YY(ITST)
           ZA = ZZ(ITST)
           BSTELA = INELM
           if( INELM .lt. 0 .and. fncdif .lt. 0.0 ) then
              print*,'*** NODTST: Node outside mesh is better!'
              print*,fncorg,aveorg,fncdif,avedif
           end if
        END IF

        IF( FNCDIF .LT. ALWAYS ) GOTO 600
        IF( AVEDIF .LT. ALWAYS ) GOTO 600

        cnt = cnt + 1
        if( cnt .lt. 6 ) goto 100

 500  CONTINUE

      ! restoring stuff (this node will be removed again, so we may not need to)
 600  nodlst(1,inod) = XO
      nodlst(2,inod) = YO
      nodlst(3,inod) = ZO
      DO I = 1, 9
         NODLST(I+6,INOD) = ORGM(I)
      END DO
      NODLST(16,INOD) = FLOAT(ORGELM)

end subroutine move_new_node_test
