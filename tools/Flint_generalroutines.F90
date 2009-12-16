MODULE GeneralRoutines

  IMPLICIT NONE

CONTAINS
 
!--------------------------------------------------------------------------------------------------

  REAL FUNCTION d3_determinant(m)

! Function to calculate the determinant of a 3x3 matrix.

    REAL, INTENT(IN) :: m(3,3)
  
    d3_determinant = ABS(      &
         m(1,1)*( m(2,2)*m(3,3) - m(2,3)*m(3,2) )  +     &
         m(1,2)*( m(2,3)*m(3,1) - m(2,1)*m(3,3) )  +     &
         m(1,3)*( m(2,1)*m(3,2) - m(2,2)*m(3,1) )        &
       )
  
  END FUNCTION d3_determinant

!--------------------------------------------------------------------------------------------------

  REAL FUNCTION d4_determinant(m)

! Function to calculate the determinant of a 4x4 matrix.

    REAL, INTENT(IN) :: m(4,4)

    d4_determinant = ABS(      &
         m(1,1)*m(2,2)*m(3,3)*m(4,4) - m(1,1)*m(2,2)*m(3,4)*m(4,3) + m(1,1)*m(2,3)*m(3,4)*m(4,2)        &
         - m(1,1)*m(2,3)*m(3,2)*m(4,4) + m(1,1)*m(2,4)*m(3,2)*m(4,3) - m(1,1)*m(2,4)*m(3,3)*m(4,2)      &

         - m(1,2)*m(2,3)*m(3,4)*m(4,1) + m(1,2)*m(2,3)*m(3,1)*m(4,4) - m(1,2)*m(2,4)*m(3,1)*m(4,3)      &
         + m(1,2)*m(2,4)*m(3,3)*m(4,1) - m(1,2)*m(2,1)*m(3,3)*m(4,4) + m(1,2)*m(2,1)*m(3,4)*m(4,3)      &

         + m(1,3)*m(2,4)*m(3,1)*m(4,2) - m(1,3)*m(2,4)*m(3,2)*m(4,1) + m(1,3)*m(2,1)*m(3,2)*m(4,4)      &
         - m(1,3)*m(2,1)*m(3,4)*m(4,2) + m(1,3)*m(2,2)*m(3,4)*m(4,1) - m(1,3)*m(2,2)*m(3,1)*m(4,4)      &

         - m(1,4)*m(2,1)*m(3,2)*m(4,3) + m(1,4)*m(2,1)*m(3,3)*m(4,2) - m(1,4)*m(2,2)*m(3,3)*m(4,1)      &
         + m(1,4)*m(2,2)*m(3,1)*m(4,3) - m(1,4)*m(2,3)*m(3,1)*m(4,2) + m(1,4)*m(2,3)*m(3,2)*m(4,1)      &
         )

  END FUNCTION d4_determinant

!--------------------------------------------------------------------------------------------------

SUBROUTINE reorder_nodes(node_coords,reordered_coords,vol_tracer,reordered_vol,int_tracer,reordered_int)

! Subroutine to reorder the nodes so that the "highest" node, in terms of its tracer value, will be
! included in the volume, instead of excluded.

! Loop variables.
  INTEGER :: j

! Original arrays.  
  REAL, DIMENSION(4,3),INTENT(IN) :: node_coords
  REAL, DIMENSION(4),INTENT(IN) :: vol_tracer, int_tracer

! Reordered arrays.
  REAL, DIMENSION(4,3),INTENT(INOUT) :: reordered_coords
  REAL, DIMENSION(4),INTENT(INOUT) :: reordered_vol, reordered_int

! Re-order the nodes so that the only one "above" T_0 is first in the array.
! If you don't do this, the edge_vectors are wrong for the small_tet routines, since they are taken
! from the node in position one of the array.
  reordered_coords(1,1:3) = node_coords(4,1:3)
  reordered_coords(2,1:3) = node_coords(3,1:3)
  reordered_coords(3,1:3) = node_coords(2,1:3)
  reordered_coords(4,1:3) = node_coords(1,1:3)

! Also re-order the vol_tracer values.
  reordered_vol(1) = vol_tracer(4)
  reordered_vol(2) = vol_tracer(3)
  reordered_vol(3) = vol_tracer(2)
  reordered_vol(4) = vol_tracer(1)

! Also re-order the int_tracer values.
  reordered_int(1) = int_tracer(4)
  reordered_int(2) = int_tracer(3)
  reordered_int(3) = int_tracer(2)
  reordered_int(4) = int_tracer(1)

END SUBROUTINE reorder_nodes

!--------------------------------------------------------------------------------------------------

REAL FUNCTION integrate(int_basis,tracer_values)

! Function to actually calculate the  value of the tracer over the selected section of the
! current element.

  REAL, DIMENSION(4) :: int_basis, tracer_values

  integrate = int_basis(1)*tracer_values(1) + int_basis(2)*tracer_values(2)      &
       + int_basis(3)*tracer_values(3) + int_basis(4)*tracer_values(4) 

END FUNCTION integrate

!--------------------------------------------------------------------------------------------------

SUBROUTINE calc_edge_vectors(edge_vectors,node_coords,n)

! Subroutine to return the edge vectors for the calculation of a 3D tetrahedron's volume.

  IMPLICIT NONE

! Loop variables.
  INTEGER :: j,k

! Passed variables.
  INTEGER :: n
  REAL, INTENT(IN) :: node_coords(n+1,n)

! Result variable.
  REAL, DIMENSION(n,n), INTENT(INOUT) :: edge_vectors

  DO j = 1,n ! Loop over nodes.
     DO k = 1,n ! Loop over coordinates.
        edge_vectors(j,k) = node_coords(j+1,k)      &
             - node_coords(1,k)
     END DO
  END DO

  RETURN

END SUBROUTINE calc_edge_vectors

!--------------------------------------------------------------------------------------------------

SUBROUTINE calc_subnode_coords(node_coords,subnode_coords,edge_fraction)

! This subroutine calculates the subnode coordinates when only one vertex of the tetrahedron has been
! clipped. It returns (x,y,z) for three subnodes placed along the edges of the main tetrahedron.

! Loop variables.
  INTEGER :: j, k

! Passed variables.
  REAL, INTENT(IN) :: edge_fraction(4), node_coords(4,3)
  REAL, INTENT(INOUT) :: subnode_coords(4,3)

! Internal variables.
  REAL :: edge_vectors(3,3)

 ! Construct edge_vectors for the entire element.
  CALL calc_edge_vectors(edge_vectors,node_coords,3)

! Scale the edge_vectors by the relevant entry of edge_fraction.
  DO j = 1,3
     DO k = 1,3
        edge_vectors(j,k) = edge_fraction(j+1)*edge_vectors(j,k)
     END DO
  END DO

! Construct the array containing the positions of one of the element's nodes and three subnodes
! that are distance "fraction" between this node and the other three.
  subnode_coords(1,:) = node_coords(1,:)
  
  DO j = 1,3 ! Loop over nodes.
     subnode_coords(j+1,:) = node_coords(1,:) + edge_vectors(j,:)
  END DO

END SUBROUTINE calc_subnode_coords

!--------------------------------------------------------------------------------------------------

SUBROUTINE calc_half_subnode_coords(node_coords,subnode_coords,edge_fraction)

! This subroutine calculates the (x,y,z) coordinates of the four subnodes formed when two vertices
! of the element are clipped.

! Loop variables
  INTEGER :: j

! Passed variables.
  REAL, INTENT(IN) :: node_coords(4,3), edge_fraction(4)
  REAL, INTENT(INOUT) :: subnode_coords(4,3)

! Internal variables.
  REAL :: evecs_1(2,3), evecs_2(2,3)

  DO j = 1,2
! Compute the edge_vectors leading from nodes 3 & 4 to node 1, to locate subnodes A & B. 
        evecs_1(j,1:3) = node_coords(1,1:3) - node_coords(j+2,1:3)

! Compute the edge_vectors leading from nodes 3 & 4 to node 2, to locate subnodes C & D.
        evecs_2(j,1:3) = node_coords(2,1:3) - node_coords(j+2,1:3)
  END DO

! We now have the vectors from node 1 to 3, 1 to 4 as well as from node 2 to 3, and 2 to 4.
! Points a,b,c,d lie somewhere along these lines depending on the temperature fraction, i.e.
! where T_0 cuts the lines between nodes. We need to multiply these vectors by the correct
! temperature fraction and add the result to the coordinates of the correct node.

  DO j = 1,2
! Compute the coordinates of subnodes A & B.
     subnode_coords(j,1:3) = node_coords(j+2,1:3) + edge_fraction(j)*evecs_1(j,1:3)

! Compute the coordinates of subnodes C & D.
     subnode_coords(j+2,1:3) = node_coords(j+2,1:3) + edge_fraction(j+2)*evecs_2(j,1:3)
  END DO

END SUBROUTINE calc_half_subnode_coords

!--------------------------------------------------------------------------------------------------

END MODULE GeneralRoutines

!--------------------------------------------------------------------------------------------------
