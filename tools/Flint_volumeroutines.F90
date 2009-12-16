MODULE VolumeRoutines

  USE GeneralRoutines

  IMPLICIT NONE

CONTAINS

!--------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!! FUNCTIONS FOR THE DIFFERENT TYPES OF VOLUME CALCULATIONS !!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------------------------------

  REAL FUNCTION add_whole_volume(node_coords)

! This adds the whole volume of the element being considered and calculates it by using the
! determinant method.

! Loop variables.
    INTEGER :: j,k
    REAL, INTENT(IN) :: node_coords(4,3)

! Dummy variables.
    REAL :: edge_vectors(3,3)

    CALL calc_edge_vectors(edge_vectors,node_coords,3)
  
    add_whole_volume = d3_determinant(edge_vectors)/6.0

  END FUNCTION add_whole_volume

!--------------------------------------------------------------------------------------------------

  REAL FUNCTION add_other_volume(node_coords,tracer_values,T_0)

! Function to add the clipped tetrahedron formed when only 1 node is excluded.

! Loop variables.
    INTEGER :: j

! Passed variables.
    REAL :: node_coords(4,3), tracer_values(4), T_0

    add_other_volume = add_whole_volume(node_coords)      &
         - add_small_tet(node_coords,tracer_values, T_0)

  END FUNCTION add_other_volume

!--------------------------------------------------------------------------------------------------

  REAL FUNCTION add_half_volume(node_coords,tracer_values,T_0)

! Function calculate the resulting volume when two nodes are excluded from the region. This is carried
! out by calculating the position of four sub-nodes that are located where the T_0 surface slices
! through the element. The required volume is then split into three sub-tets, whose volume can be
! calculated using the determinant method.

! Loop variables
    INTEGER :: j, k

! Passed variables
    REAL :: node_coords(4,3), tracer_values(4), T_0

! Internal variables
    REAL :: evecs_1(2,3), evecs_2(2,3), edge_vectors(3,3)
    REAL :: subnode_coords(4,3), subtet_coords(4,3)
    REAL :: T_fraction(4)

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

! Nodes 3 & 4 lie "above" T_0. I don't reorder, ala Geoff, as I find this confusing...

!---
! Compute the relevant T_fraction for the four subnodes.

    DO j = 1,2
! For subnodes A & B.
       T_fraction(j) = ABS( tracer_values(j+2) - T_0 )/( tracer_values(j+2) - tracer_values(1) )

! For subnodes C & D.
       T_fraction(j+2) = ABS( tracer_values(j+2) - T_0 )/( tracer_values(j+2) - tracer_values(2) )
    ENDDO

!---
! Calculate the coordinates of the sub-nodes.

    CALL calc_half_subnode_coords(node_coords,subnode_coords,T_fraction)

!---
! Compile the subtet_coords_? arrays for the three subtets.

! Compute the integral of the basis function associated with the included node (node "4")

! Reset the hyper_coords array to have the initial hyperdimensional point at node "4".
    subtet_coords(1,1:3) = subnode_coords(1,1:3)
    subtet_coords(2,1:3) = subnode_coords(3,1:3)
    subtet_coords(3,1:3) = node_coords(3,1:3)
    subtet_coords(4,1:3) = node_coords(4,1:3)

 ! Compute the volume of subtet 1.
    CALL calc_edge_vectors(edge_vectors,subtet_coords,3)
    add_half_volume = d3_determinant(edge_vectors)/6.0

! Compute the volume of subtet 2.
    subtet_coords(3,1:3) = subnode_coords(2,1:3)
    CALL calc_edge_vectors(edge_vectors,subtet_coords,3)
    add_half_volume = add_half_volume + d3_determinant(edge_vectors)/6.0

! Compute the volume of subtet 3.
    subtet_coords(1,1:3) = subnode_coords(4,1:3)
    CALL calc_edge_vectors(edge_vectors,subtet_coords,3)
    add_half_volume = add_half_volume + d3_determinant(edge_vectors)/6.0

  END FUNCTION add_half_volume

!--------------------------------------------------------------------------------------------------

  REAL FUNCTION add_small_tet(node_coords,tracer_values,T_0)

! This subroutine adds the small tetrahedron that occurs when only one node sticks above the T_0
! surface.

! Loop variables.
    INTEGER :: j, k

! Passed variables.
    REAL :: node_coords(4,3), tracer_values(4), T_0

! node_coords - the (x,y,z) coordinates of the nodes associated with the current element.
! tracer_values - the nodal values of the tracer for the current element

! Internal variables.
    REAL :: edge_vectors(3,3), T_fraction(3)

! Calculate the fraction of the distance between the nodal pairs that the T_0 surface cuts through
! each edge at. Uses the fact that we KNOW variation along an edge is linear.
    DO j = 1,3
       T_fraction(j) = ABS(      &
            ( tracer_values(1) - T_0 )     &
            /( tracer_values(1) - tracer_values(j+1) )      &
            )
    END DO

! Compute the edge_vectors for the entire element.
    CALL calc_edge_vectors(edge_vectors,node_coords,3)

! Shorten the edge vectors to take into account that only part of the element is within the region
! of interest.
    DO j = 1,3
       edge_vectors(j,:) = T_fraction(j)*edge_vectors(j,:)
    END DO

! Calculate the partial volume
    add_small_tet = d3_determinant(edge_vectors)/6.0

  END FUNCTION add_small_tet

!--------------------------------------------------------------------------------------------------

END MODULE VolumeRoutines

!--------------------------------------------------------------------------------------------------
