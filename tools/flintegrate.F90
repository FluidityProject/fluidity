subroutine partial_integrals( &
     simulationname, simulationname_len, &
     dumpnumber, dumpnumber_len, &
     vol_tracer_name, vol_tracer_name_len, &
     int_tracer_name, int_tracer_name_len)

  USE GeneralRoutines
  USE VolumeRoutines
  USE IntegralRoutines

  IMPLICIT NONE

!--------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BEGIN VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------------------------------
! Loop variables.

  INTEGER :: a  ! Loops through all the tracer classes
  INTEGER :: i  ! Loops through all the elements
  INTEGER :: j, k, l, m  ! Other loop variables

!--------------------------------------------------------------------------------------------------
! Nodecounter keeps track of how many nodes are above/below the given T_0 value.

  INTEGER :: nodecounter                         

!--------------------------------------------------------------------------------------------------
! Counters to record how many times a given partial volume calculation subroutine is called.

  INTEGER :: add_whole_counter, add_other_counter, add_half_counter, add_small_counter
  INTEGER :: nosubcounter, subcounter

! Total number of subroutines called
  INTEGER :: no_subs_called

! Counters for any problems in the calculation.
  INTEGER :: problem_counter  

!--------------------------------------------------------------------------------------------------
! Variables for external data.

! Total number of elements/nodes in the domain.
  INTEGER :: total_elements, total_nodes  

! Node-wise tracer values, vol_tracer will be used to define the partial element volumes,
! int_tracer will be integrated over them.
  REAL, ALLOCATABLE, DIMENSION(:) :: vol_tracer, int_tracer

! Connectivity data.
  REAL, ALLOCATABLE, DIMENSION(:,:) :: node_groups

! Nodal co-ordinates.
   REAL, ALLOCATABLE, DIMENSION(:,:) :: node_positions

!--------------------------------------------------------------------------------------------------
! Variables for the partial volume calculation.

! Number of tracer classes to be considered.
  INTEGER :: no_classes

! Current tracer value under consideration and the volume to be added from the current element?
  REAL :: T_0, volume_to_add

! Array of all the tracer values to be considered.
  REAL, ALLOCATABLE, DIMENSION(:) :: T_0_values

! To store the partial volumes - ith partial volume is the volume above the ith T_0 value
  REAL, ALLOCATABLE, DIMENSION(:) :: p, partial_volume

! The volume of water within two consecutive T_0 values 
  REAL, ALLOCATABLE, DIMENSION(:) :: volume_in_class                              

!--------------------------------------------------------------------------------------------------
! Variables for the identification of the correct node for a specific element

! Global node number of the node for the current element.
  INTEGER :: current_el_node_nos(4)      

! Node coordinates and nodal tracer values for the current element.
  REAL :: current_el_node_coords(4,3), current_el_int_tracer(4), current_el_vol_tracer(4)    
! The above, but reordered so that the node to be included in the volume/integral is at the top
  REAL :: reordered_node_coords(4,3), reordered_vol_tracer(4), reordered_int_tracer(4)

! Find the location  along an edge where T_0 is.
  REAL :: tracer_fraction(3) 

! These are to do with reordering the nodes for each calculation
  INTEGER :: lptr
  REAL :: temp, temp2, temp3(3), temp4

!--------------------------------------------------------------------------------------------------

! s is the volume of the entire domain, which is used as a sanity check at the end to ensure that the
! sum of the volume classes is correct. Other variables are part of this sanity check.
  REAL :: s                                    
  REAL :: volume_difference, volume_percentage, volume_total

  ! Info for naming of files
  
  ! simulationname  - Project name.
  ! dumpnumber      - Dump number
  ! vol_tracer_name - 'Name' of the field to be used to define the volume classes.
  ! int_tracer_name - 'Name' of the field to be integrated over the volume classes.
  
  integer, intent(in)::dumpnumber_len, simulationname_len, int_tracer_name_len, vol_tracer_name_len
  CHARACTER (len = dumpnumber_len) :: dumpnumber
  CHARACTER (len = simulationname_len) :: simulationname
  CHARACTER (len = int_tracer_name_len) :: int_tracer_name
  CHARACTER (len = vol_tracer_name_len) :: vol_tracer_name
  
  CHARACTER (len = 500) :: filename

!--------------------------------------------------------------------------------------------------
! Variables for the integration over the partial volumes.

! The integral of the tracer field over an element or partial element.
  REAL :: integral_to_add
  
! Integral of the tracer field over the entire class.
  REAL, ALLOCATABLE, DIMENSION(:) :: partial_integral, integral_in_class

! To tell whether files really exist or not.
  LOGICAL :: GOT_FILE

!--------------------------------------------------------------------------------------------------
! Variables and interface block for subroutine to read in the vtk header information.

  INTEGER :: err, flength, namelen, szenls, nfields, nproperties, ndimensions, maxnamelen

  INTERFACE
     INTEGER FUNCTION fgetvtksizes( filename, namelen, nnod, nelm, szenls,     &
          nfields, nproperties, ndimensions, maxnamelen )
       CHARACTER(len=*) :: filename
       INTEGER :: namelen ,nnod, nelm, szenls, nfields, nproperties, ndimensions, maxnamelen 
     END FUNCTION fgetvtksizes
  END INTERFACE

!--------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END OF VARIABLE DECLARATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------------------------------
! BEGIN CODE PROPER

! Read in ClassData, which contains the number of classes and the threshold limits for each
! class.
  INQUIRE(file = 'ClassData',exist = GOT_FILE)
  IF( GOT_FILE )THEN
     OPEN(10,file = 'ClassData')
     READ(10,'(I4)') no_classes
     ALLOCATE( T_0_values (no_classes+1) )
     READ(10,*) T_0_values
     CLOSE(10)
  ELSE
     PRINT*,'ERROR: ClassData file absent, this file is required for processing. STOPPING.'
     STOP
  ENDIF

! Call fgetvtksizes here to get all the file header information.
  filename = TRIM(simulationname)//'_'//TRIM(dumpnumber)//'.vtu'
  flength = LEN(TRIM(filename))
  ndimensions = 0     ! This should be initialised as 0 to ensure the correct number of dimensions
                      ! are extracted by fgetvtksizes.

! Extract the header information from the vtk file.
  err = fgetvtksizes(filename, flength, total_nodes, total_elements,     &
       szenls, nfields, nproperties, ndimensions, maxnamelen) 

  IF( err /= 0 )THEN
     PRINT *, 'Something went wrong with VTK sizes!'
     STOP
  END IF

!--------------------------------------------------------------------------------------------------
! Allocate some variables (see above for definitions).

  ALLOCATE( vol_tracer (total_nodes) )
  ALLOCATE( int_tracer (total_nodes) )
  ALLOCATE( node_groups (total_elements,4) )
  ALLOCATE( node_positions (total_nodes,3) )

  ALLOCATE( partial_volume (no_classes) )
  ALLOCATE( volume_in_class (no_classes) )
  ALLOCATE( partial_integral (no_classes) )
  ALLOCATE( integral_in_class (no_classes) )

!--------------------------------------------------------------------------------------------------

  PRINT*, 'There are',total_nodes,'nodes, together they make',total_elements,'elements.'

  CALL ExtractData( filename, flength, total_elements, total_nodes, node_groups,szenls, nfields,  &
       nproperties, ndimensions, maxnamelen, vol_tracer, vol_tracer_name,      &
       int_tracer, int_tracer_name, node_positions)

!--------------------------------------------------------------------------------------------------
! Set the subroutine entry counters to zero.

  add_whole_counter = 0
  add_other_counter = 0
  add_half_counter = 0
  add_small_counter = 0
  problem_counter= 0
  subcounter = 0
  nosubcounter = 0

!--------------------------------------------------------------------------------------------------
! Initialise arrays.

  integral_to_add = 0.0
  partial_integral(:) = 0.0
  integral_in_class(:) = 0.0
  volume_to_add = 0.0
  partial_volume(:) = 0.0
  volume_in_class(:) = 0.0
  volume_total = 0

!--------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!! BEGIN MAIN LOOP OVER NUMBER OF TRACER CLASSES !!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------------------------------
! Loop over classes - class number 1 = lowest tracer value.

  DO a = 1,no_classes

!--------------------------------------------------------------------------------------------------
! Pick out a single tracer value to use for the decision of are we above or below it!
     T_0 = T_0_values(a)

!--------------------------------------------------------------------------------------------------
! Loop over all the elements.
     DO i = 1,total_elements

!--------------------------------------------------------------------------------------------------

! Pick out the global node numbers corresponding to the nodes belonging to the current element.
        DO j = 1,4
           current_el_node_nos(j) = node_groups(i,j)
        END DO

! Pick out the tracers corresponding to the global node numbers contained in 
! current_el_node_nos

        DO k = 1,4
           current_el_vol_tracer(k) = vol_tracer(current_el_node_nos(k))
           current_el_int_tracer(k) = int_tracer(current_el_node_nos(k))
        END DO

! Do the same for their (x,y,z) coordinates.
        DO l = 1,4
           DO m = 1,3
              current_el_node_coords(l,m) = node_positions((current_el_node_nos(l)),m)
           END DO
        END DO

! Initialise nodecounter.
        nodecounter=0

!--------------------------------------------------------------------------------------------------
! This next bit orders the nodes (keeping node numbers and node tracers in corresponding array
! entries!) for the element under consideration so that they run from lowest tracer to highest
! in the current_el_vol_tracer array.  This is necessary for the way the subroutines work.

        DO l=1,3

           lptr=l

           DO m=l+1,4
              IF( current_el_vol_tracer(m) < current_el_vol_tracer(lptr) ) lptr = m
           END DO

           IF( l /= lptr )THEN
              temp = current_el_vol_tracer(l)
              temp2 = current_el_node_nos(l)
              temp3 = current_el_node_coords(l,:)
              temp4 = current_el_int_tracer(l)
              
              current_el_vol_tracer(l) = current_el_vol_tracer(lptr)
              current_el_node_nos(l) = current_el_node_nos(lptr)
              current_el_node_coords(l,:) = current_el_node_coords(lptr,:)
              current_el_int_tracer(l) = current_el_int_tracer(lptr)
             
              current_el_vol_tracer(lptr) = temp
              current_el_node_nos(lptr) = temp2
              current_el_node_coords(lptr,:) = temp3
              current_el_int_tracer(lptr) = temp4
           END IF

        END DO

!--------------------------------------------------------------------------------------------------
! Decide which partial volume calculation method to use based on the number of nodes with tracer
! values > T_0.

! Use of .ge. instead of .gt. eliminates many of the checks required for nodes sitting on the T_0
! surface below.
        DO j = 1,4
           IF( current_el_vol_tracer(j) >= T_0 ) nodecounter = nodecounter+1
        END DO

!--------------------------------------------------------------------------------------------------
! Pick a volume calculation method based upon the value of nodecounter.

! If nodecounter is 4, all four nodes of the current element have tracer values > T_0
! => include all of its volume.

        IF( nodecounter == 4 )THEN

           volume_to_add = add_whole_volume(current_el_node_coords)

           integral_to_add = integrate_whole(current_el_node_coords,current_el_int_tracer)

           add_whole_counter = add_whole_counter + 1
           subcounter = subcounter + 1

!--------------------------------------------------------------------------------------------------
! If nodecounter is 3, then only 3 nodes are "above" T_0 and one node gets snipped off.

        ELSEIF( nodecounter == 3 ) THEN 

           volume_to_add = add_other_volume(current_el_node_coords,current_el_vol_tracer,T_0)

           integral_to_add = integrate_other(current_el_node_coords,current_el_vol_tracer,T_0,      &
                current_el_int_tracer)

           add_other_counter = add_other_counter + 1
           subcounter = subcounter + 1

!--------------------------------------------------------------------------------------------------
! If nodecounter is 2, then the elements has been split into two funny shapes that need to be
! subdivided into tetrahedra t calculate their volume.
           
        ELSEIF( nodecounter == 2 )THEN

           volume_to_add = add_half_volume(current_el_node_coords,current_el_vol_tracer,T_0)

           integral_to_add = integrate_half(current_el_node_coords,current_el_vol_tracer,T_0,      &
                current_el_int_tracer)

           add_half_counter = add_half_counter + 1
           subcounter = subcounter + 1

!--------------------------------------------------------------------------------------------------
! If nodecounter is 1, only one element is "above" the tracer value being considered. 
          
        ELSEIF( nodecounter == 1 )THEN

! Have to give the functions node arrays that are reordered such that node "4" is at the first index.
! Otherwise the "lowest" node of the element will be considered the one to be included.
           CALL reorder_nodes(current_el_node_coords,reordered_node_coords,     &
                current_el_vol_tracer,reordered_vol_tracer,current_el_int_tracer,reordered_int_tracer)

           volume_to_add = add_small_tet(reordered_node_coords,reordered_vol_tracer,T_0)

           integral_to_add = integrate_small(reordered_node_coords,reordered_vol_tracer,T_0,      &
                reordered_int_tracer)

           add_small_counter = add_small_counter + 1
           subcounter = subcounter + 1

!--------------------------------------------------------------------------------------------------
! Finish IF loop to determine which partial volume subroutine to use.
        ELSE
        END IF

!--------------------------------------------------------------------------------------------------

! Increment the counter that records when no partial volume subroutine was called.
        IF( nodecounter == 0 ) nosubcounter = nosubcounter + 1

! Increment the partial volume of the tracer class.
        partial_volume(a) = partial_volume(a) + volume_to_add
   
! Re-initialise this so that no double-counting takes place.
        volume_to_add = 0.0

!--------------------------------------------------------------------------------------------------

! Increment the partial integral over the partial volume.
        partial_integral(a) = partial_integral(a) + integral_to_add
 
! Re-initialise the partial integral variable to prevent double-counting.
        integral_to_add = 0.0

!--------------------------------------------------------------------------------------------------
! End of loop over elements.

     END DO

!--------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!! END OF MAIN LOOP OVER NUMBER OF TRACER CLASSES !!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------------------------------

  END DO

!--------------------------------------------------------------------------------------------------
! Calculate the volume in an actual class, instead of the volume "above" it.

  DO i = 1,no_classes-1
     volume_in_class(i) = partial_volume(i) - partial_volume(i+1)
  END DO
  
  volume_in_class(no_classes) = partial_volume(no_classes)

!--------------------------------------------------------------------------------------------------
! Calculate the integral in an actual class, instead of the integral "above" that classes value of
! T_0.

  DO i = 1,no_classes-1
     integral_in_class(i) = partial_integral(i) - partial_integral(i+1)
  END DO
  
  integral_in_class(no_classes) = partial_integral(no_classes)

!--------------------------------------------------------------------------------------------------
! Writing out data for plotting

  OPEN(20, FILE = TRIM(simulationname)//'_'//TRIM(dumpnumber)//'.flint')
  WRITE(20,'(A1)')'@'
  WRITE(20,'(5A20)')'@(1):CLASS           ','(2):T0-THRESHOLD    ','(3):PARTIAL-VOL.    ',     &
       '(4):CLASS-VOL.      ','(5):PARTIAL-INT.    '
  WRITE(20,'(1A20)')'@(6):CLASS-INT.     '
  WRITE(20,'(1A3)')'@@6'

  DO a = 1,no_classes
     WRITE(20, '(6E13.6)') REAL(a), T_0_values(a), partial_volume(a), volume_in_class(a),      &
          partial_integral(a), integral_in_class(a)
  END DO
  CLOSE(20)

!--------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SANITY CHECKS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------------------------------

! Volume conservation - sum of all the partial class volumes = sum of all element volumes.

! Calculate the total volume of the domain by summing the individual element volumes.
  DO i = 1,total_elements
 
! Pick out the global node numbers corresponding to the nodes belonging to the current element.
     DO j = 1,4
        current_el_node_nos(j) = node_groups(i,j)
     END DO
        
! Pick out their x, y, z coordinates.
     DO l = 1,4
        DO m = 1,3
           current_el_node_coords(l,m) = node_positions((current_el_node_nos(l)),m)
        END DO
     END DO

! Calculate the volume of each element and keep a running total.
     volume_total = volume_total + add_whole_volume(current_el_node_coords)

  END DO

  s = SUM( volume_in_class )

! To conserve volume, the difference between s and volume_total should be big fat zero.
  volume_difference = s - volume_total

  PRINT*,' '
  PRINT*, s, 'is the sum of tracer class volumes which should also equal'
  PRINT*, volume_total, 'the sum of element volumes'

  volume_percentage = 100.0*( volume_difference/volume_total )

  PRINT*,' '

  IF( volume_difference == 0 )THEN
     PRINT*, 'You have conserved volume!'
  ELSEIF( volume_difference > 0 )THEN
     PRINT*, 'This is a gain of',volume_percentage,'% of the sum of element volumes'
  ELSEIF( volume_difference < 0 )THEN
     PRINT*, 'This is a loss of',ABS( volume_percentage ),'% of the sum of element volumes'
  END IF

  PRINT*,' '

!--------------------------------------------------------------------------------------------------
! Check whether the program is working properly.
! NOTE - It is necessary but NOT sufficient

 DO i = 1,no_classes-1

     IF( partial_volume(i) < partial_volume(i+1) )THEN

        volume_difference = partial_volume(i+1) - partial_volume(i)
        PRINT*, 'There is a volume class problem'
        
        problem_counter = problem_counter + 1
        p = i+1

        PRINT*, 'partial_volume',i,partial_volume(i),'is less than partial volume',     &
             p, partial_volume(i+1)
        PRINT*,'The difference in between volume classes is',volume_difference

     END IF

  END DO

  PRINT*, 'There are',problem_counter,     &
       'problems (note - zero problems is necessary but not sufficient)'
! This means that there are at least "problem_counter" wrong calculations.
! There may be more - one wrong calculation can mask other.

!--------------------------------------------------------------------------------------------------
! A tally of how many times each subroutine was called

  PRINT*,'There were',nosubcounter,'instances of no subroutine being called'
  PRINT*, 'add_whole_volume was called',add_whole_counter,'times'
  PRINT*, 'add_other_volume was called',add_other_counter,'times'
  PRINT*, 'add_half_volume was called',add_half_counter,'times'
  PRINT*, 'add_small_tet was called',add_small_counter,'times'

! Check on double counting - each element should have been considered once for each class.
  PRINT*, total_elements*no_classes, 'should equal', nosubcounter + add_whole_counter     &
       + add_other_counter + add_half_counter + add_small_counter

!--------------------------------------------------------------------------------------------------
! Deallocate the memory before ending.

  DEALLOCATE( vol_tracer )
  DEALLOCATE( int_tracer )
  DEALLOCATE( node_groups )
  DEALLOCATE( node_positions )

  DEALLOCATE( T_0_values )

  DEALLOCATE( partial_volume )
  DEALLOCATE( volume_in_class  )
  DEALLOCATE( partial_integral )
  DEALLOCATE( integral_in_class )

!--------------------------------------------------------------------------------------------------
! End of program.

END subroutine partial_integrals

!--------------------------------------------------------------------------------------------------
