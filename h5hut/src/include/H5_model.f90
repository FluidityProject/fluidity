!
!  Copyright (c) 2006-2016, The Regents of the University of California,
!  through Lawrence Berkeley National Laboratory (subject to receipt of any
!  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
!  Institut (Switzerland).  All rights reserved.!
!
!  License: see file COPYING in top level of source distribution.
!
  INTERFACE
     !> \addtogroup h5_model_f
     !! @{
     
     !>
     !! See \ref H5HasStep
     !! \return 0 on success or H5_FAILURE

     LOGICAL FUNCTION h5_hasstep (filehandle,step)
       INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
       INTEGER*8, INTENT(IN) :: step       !< a timestep value >= 1
     END FUNCTION h5_hasstep

     !>
     !! See \ref H5SetStep
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5_setstep (filehandle,step)
       INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
       INTEGER*8, INTENT(IN) :: step       !< a timestep value >= 1
     END FUNCTION h5_setstep

     !>
     !! See \ref H5GetStep
     !! \return the the current step or \c H5_FAILURE

     INTEGER*8 FUNCTION h5_getstep (filehandle)      
       INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
     END FUNCTION h5_getstep

     !>
     !! See \ref H5GetNumSteps
     !! \return the number of steps or error code

     INTEGER*8 FUNCTION h5_getnsteps (filehandle)      
       INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
     END FUNCTION h5_getnsteps

     !> @}
  END INTERFACE
