MODULE IntegralRoutines
  
  USE GeneralRoutines

  IMPLICIT NONE

CONTAINS

!--------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!! SUBROUTINES FOR THE DIFFERENT TYPES OF BASIS FUNCTION INTEGRATION !!!!!!!!!!!!!!!
!--------------------------------------------------------------------------------------------------

  REAL FUNCTION integrate_whole(node_coords,int_tracer)

! 14/03/07 - integrates all four basis functions over the whole element, multiplies by the nodal
!            value and sums.

! Loop variables.
    INTEGER :: i, j, k

! Passed variables.
    REAL :: node_coords(4,3), int_tracer(4)

! Edge vectors and positions of the nodes for the hyper-dimensional tetrahedron.
    REAL :: edge_vectors(4,4), hyper_coords(5,4)

! Value of the integrated basis function for each node.
    REAL, DIMENSION(4) :: int_basis

! Set the invariant sections of the array containing the coordinates of the hyper-dimensional
! tetrahedron.
    DO i = 1,4 ! Loop over the nodes.
       DO j = 1,3 ! Loop over the (x,y,z) coordinates
          hyper_coords(i,j) = node_coords(i,j)
       END DO
    END DO

    hyper_coords(:,4) = 0.0 ! Set the basis function value at the vertices of the element.
    hyper_coords(5,4) = 1.0 ! Set the basis function value at the hyper-dimensional point.

! Loop over all the nodes for the calculation.
    DO i = 1,4

! Set the correct (x,y,z) coords of the hyper-dimensional point.
       hyper_coords(5,1:3) = node_coords(i,1:3)

! Set the edge vectors.
       CALL calc_edge_vectors(edge_vectors,hyper_coords,4)

! Integrate the basis function.
       int_basis(i) = d4_determinant(edge_vectors)/24.0

    END DO

! Compute the value of the field integrated over the entire element.
    integrate_whole = integrate(int_basis,int_tracer)

    IF( (int_basis(1) /= int_basis(2)).or.(int_basis(1) /= int_basis(3))    &
         .or.(int_basis(1) /= int_basis(4)) )THEN
       PRINT*,'integration is FUBAR'
!     STOP
    END IF

  END FUNCTION integrate_whole

!--------------------------------------------------------------------------------------------------

  REAL FUNCTION integrate_other(node_coords,vol_tracer,T_0,int_tracer)

! Function integrate the tracer field over the truncated tetrahedron formed when the lowest vertex
! is clipped and three nodes are included in the current tracer class.

! Passed variables.
    REAL :: node_coords(4,3), vol_tracer(4), T_0, scratch, int_tracer(4)

    integrate_other = integrate_whole(node_coords,int_tracer)      &
         - integrate_small(node_coords,vol_tracer,T_0,int_tracer)

  END FUNCTION integrate_other

!--------------------------------------------------------------------------------------------------
  
  REAL FUNCTION integrate_half(node_coords,vol_tracer,T_0,int_tracer)

! Function to integrate the tracer field over the shape formed when two nodes lie outside the
! defined volume and to lie inside.
! The included section of the element is divided into three sub-tets and each sub-tet has all four
! basis functions integrated over it. I expect that this routine is the most general and would
! work for all four possible nodecounter values.

! Loop variables.
    INTEGER :: j, k, l

! Passed variables.
    REAL :: node_coords(4,3), vol_tracer(4), T_0, int_tracer(4)

! Internal variables.
    REAL :: edge_vectors(4,4), hyper_coords(5,4), subnode_coords(4,3), subtet_coords(3,4,3)
    REAL :: T_fraction(4), ones(4)

! Arrays to store the integrated values of the basis functions.
    REAL :: basis(3,4,4), int_basis(4), phi

! evecs_? : the edge_vectors leading from nodes 3 & 4 to nodes 1 & 2.
! subnode_coords : an array containing the (x,y,z) coordinates of subnodes A, B, C, & D.
! subtet_coords : an array that contains the (x,y,z) coordinates of the subtet whose volume is
!                 being found. A combination of the coordinates of nodes 3 & 4 and the subnodes.

!---

! This subroutine requires 4 subnodes, which are going to sequentially arranged in all the relevant
! arrays in the order that follows.

! subnode A : lies between nodes 1 & 3
! subnode B : lies between nodes 1 & 4
! subnode C : lies between nodes 2 & 3
! subnode D : lies between nodes 2 & 4

! Nodes 3 & 4 lie "above" T_0, with 4 having the highest tracer value and 1 the lowest.

!---
! Compute the relevant T_fraction for the four subnodes.

    DO j = 1,2
! For subnodes A & B.
       T_fraction(j) = ABS( vol_tracer(j+2) - T_0 )/( vol_tracer(j+2) - vol_tracer(1) )

! For subnodes C & D.
       T_fraction(j+2) = ABS( vol_tracer(j+2) - T_0 )/( vol_tracer(j+2) - vol_tracer(2) )
    ENDDO

!---
! Calculate the coordinates of the sub-nodes.

    CALL calc_half_subnode_coords(node_coords,subnode_coords,T_fraction)

!---
! Initialise the basis function arrays.

    int_basis(1:4) = 0.0
    basis(:,:,:) = 0.0

!---
! Compile the subtet_coords array for subtet 1.
    subtet_coords(1,1,1:3) = node_coords(4,1:3)    ! 4
    subtet_coords(1,2,1:3) = subnode_coords(1,1:3) ! A
    subtet_coords(1,3,1:3) = subnode_coords(3,1:3) ! C
    subtet_coords(1,4,1:3) = node_coords(3,1:3)    ! 3

    basis(1,1,4) = 1.0
    
    basis(1,2,1) = T_fraction(1)
    basis(1,2,3) = 1.0 - T_fraction(1)

    basis(1,3,2) = T_fraction(3)
    basis(1,3,3) = 1.0 - T_fraction(3)

    basis(1,4,3) = 1.0

!---
! Compile the subtet_coords array for subtet 2.
    subtet_coords(2,1,1:3) = node_coords(4,1:3)    ! 4
    subtet_coords(2,2,1:3) = subnode_coords(1,1:3) ! A
    subtet_coords(2,3,1:3) = subnode_coords(2,1:3) ! B
    subtet_coords(2,4,1:3) = subnode_coords(3,1:3) ! C
    
    basis(2,1,4) = 1.0

    basis(2,2,1) = T_fraction(1)
    basis(2,2,3) = 1.0 - T_fraction(1)

    basis(2,3,1) = T_fraction(2)
    basis(2,3,4) = 1.0 - T_fraction(2)

    basis(2,4,2) = T_fraction(3)
    basis(2,4,3) = 1.0 - T_fraction(3)

!---
! Compile the subtet_coords array for subtet 3.
    subtet_coords(3,1,1:3) = node_coords(4,1:3)    ! 4
    subtet_coords(3,2,1:3) = subnode_coords(2,1:3) ! B
    subtet_coords(3,3,1:3) = subnode_coords(3,1:3) ! C
    subtet_coords(3,4,1:3) = subnode_coords(4,1:3) ! D

    basis(3,1,4) = 1.0

    basis(3,2,1) = T_fraction(2)
    basis(3,2,4) = 1.0 - T_fraction(2)

    basis(3,3,2) = T_fraction(3)
    basis(3,3,3) = 1.0 - T_fraction(3)

    basis(3,4,2) = T_fraction(4)
    basis(3,4,4) = 1.0 - T_fraction(4)

!---

! In this case, the integration needs to be in the most general sense for all 4 basis functions in all
! 3 subtets, otherwise it all has to be individually tailored and the subroutine gets huge.

    DO j = 1,3 ! Loop over the subtets
       DO k = 1,4 ! Loop over the basis functions

! Initialise the hyper_coords array for the current subtet
          hyper_coords(1:4,1:3) = subtet_coords(j,1:4,1:3)
          hyper_coords(5,1:3) = subtet_coords(1,1,1:3) ! Node 4 is used in all three subtets.
          
          hyper_coords(1:4,4) = 0.0
          hyper_coords(5,4) = basis(j,1,k) ! phi for the current basis function at node 4.
          
! Compute the integral of the current basis function over the real part of the hyper-tet.
          CALL calc_edge_vectors(edge_vectors,hyper_coords,4)
          int_basis(k) = int_basis(k) + d4_determinant(edge_vectors)/24.0
        
! Loop over the other vertices of the subtet.
          DO l = 2,4

             phi = basis(j,l,k) ! value of current basis function at current vertex.

! Permute the hyper_coords array to include the sections of the hyper-tet outside of physical space.
             hyper_coords(l-1,1:3) = subtet_coords(j,l,1:3)
             hyper_coords(l-1,4) = phi

             CALL calc_edge_vectors(edge_vectors,hyper_coords,4)
             int_basis(k) = int_basis(k) + d4_determinant(edge_vectors)/24.0
           
          END DO
       END DO
    END DO

! Compute the value of the field integrated over the entire element.
    integrate_half = integrate(int_basis,int_tracer)

  END FUNCTION integrate_half

!--------------------------------------------------------------------------------------------------

  REAL FUNCTION integrate_small(node_coords,vol_tracer,T_0,int_tracer)

! Function to integrate all four basis functions over the small tetrahedron formed when only one node
! is "above" the T_0 surface.

! Loop variables.
    INTEGER :: i, j, k

! Passed variables.
    REAL, INTENT(IN) :: node_coords(4,3), vol_tracer(4), T_0, int_tracer(4)

! Edge vectors and positions of the nodes for the hyper-dimensional tetrahedron.
    REAL :: edge_vectors(4,4), hyper_coords(5,4), phi, T_fraction(4)

! Subnode coordinates array,
    REAL :: subnode_coords(4,3)

! Value of the integrated basis function for each node.
    REAL, DIMENSION(4) :: int_basis, scratch

! Calculate the fraction of the distance between the nodal pairs that the T_0 surface cuts through
! each edge at. Uses the fact that we KNOW variation along an edge is linear.
    T_fraction(1) = 1.0
    DO j = 2,4
       T_fraction(j) = ABS(      &
            ( vol_tracer(1) - T_0 )     &
            /( vol_tracer(1) - vol_tracer(j) )      &
            )
    END DO

! Compute the (x,y,z) coordinates of the three subnodes that delineate the small tet.
    CALL calc_subnode_coords(node_coords,subnode_coords,T_fraction)

!----
! Compute the integral of the basis function for the basis functions associated with the three excluded
! nodes.

! Set the invariant sections of the array containing the coordinates of the hyper-dimensional tetrahedron.
    DO i = 1,4 ! Loop over the nodes.
       hyper_coords(i,1:3) = subnode_coords(i,1:3)
    END DO

    hyper_coords(1:5,4) = 0.0 ! Set the basis function value at the vertices of the element.

! Only loop over the subnodes (P,Q,R), because their abbreviated basis functions form hyper-tets.
    DO i = 2,4

! Set the correct (x,y,z) coords of the hyper-dimensional point.
       hyper_coords(5,1:3) = subnode_coords(i,1:3)
! Set the basis function value at the hyper-dimensional point.
       hyper_coords(5,4) = T_fraction(i)

! Set the edge vectors.
       CALL calc_edge_vectors(edge_vectors,hyper_coords,4)

! Integrate the basis function.
       int_basis(i) = d4_determinant(edge_vectors)/24.0

    END DO

!----

! Compute the integral of the basis function associated with the included node (node "4")

! Reset the hyper_coords array to have the initial hyperdimensional point at node "4".
    hyper_coords(5,1:3) = subnode_coords(1,1:3)
    
    hyper_coords(1:4,4) = 0.0
    hyper_coords(5,4) = 1.0

! Compute the integral of this basis function over the real part of the hyper-tet.
    CALL calc_edge_vectors(edge_vectors,hyper_coords,4)
    int_basis(1) = d4_determinant(edge_vectors)/24.0
    
! Loop over the points P,Q,R, which are the hyper-dimensional points of the basis function along the
! edges of the element heading towards nodes "1", "2", and "3", respectively.
    DO i = 2,4

       phi = 1.0 - T_fraction(i)
     ! This is the value of the node "4" basis function at the three sub-nodes.

! Permute the hyper_coords array to include the sections of the hyper-tet outside of physical space.
       hyper_coords(i-1,1:3) = subnode_coords(i,1:3)
       hyper_coords(i-1,4) = phi

       CALL calc_edge_vectors(edge_vectors,hyper_coords,4)
       int_basis(1) = int_basis(1) + d4_determinant(edge_vectors)/24.0
 
    END DO

!----

! Compute the value of the field integrated over the entire element.
    integrate_small = integrate(int_basis,int_tracer)

  END FUNCTION integrate_small

!--------------------------------------------------------------------------------------------------

END MODULE IntegralRoutines

!--------------------------------------------------------------------------------------------------
