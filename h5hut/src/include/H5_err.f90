!
!  Copyright (c) 2006-2016, The Regents of the University of California,
!  through Lawrence Berkeley National Laboratory (subject to receipt of any
!  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
!  Institut (Switzerland).  All rights reserved.!
!
!  License: see file COPYING in top level of source distribution.
!
  INTERFACE

     !> \addtogroup h5_error_f
     !! @{
     
     !>
     !! Report error, do not abort program. The error must be handled in the programm.

     SUBROUTINE h5_report_on_error ()
     END SUBROUTINE h5_report_on_error
     
     !>
     !! Abort program on error.

     SUBROUTINE h5_abort_on_error ()
     END SUBROUTINE h5_abort_on_error
     
     !>
     !! Get last error code.
     !!
     !! Error codes are:
     !!
     !! - \c H5_ERR_BADF:	Something is wrong with the file handle.
     !! - \c H5_ERR_NOMEM:	Out of memory.
     !! - \c H5_ERR_INVAL:	Invalid argument.
     !! 
     !! - \c H5_ERR_VIEW:	Something is wrong with the view.
     !! - \c H5_ERR_NOENTRY:	A lookup failed.
     !! 
     !! - \c H5_ERR_MPI:	A MPI error occured.
     !! - \c H5_ERR_HDF5:	A HDF5 error occured.
     !! - \c H5_ERR_H5: 	Unspecified error in H5 module.
     !! - \c H5_ERR_H5PART:	Unspecified error in H5Part module.
     !! - \c H5_ERR_H5BLOCK:	Unspecified error in H5Block module.
     !! - \c H5_ERR_H5FED:	Unspecified error in H5Fed module.
     !! 
     !! - \c H5_ERR_INTERNAL:	Internal error.
     !! - \c H5_ERR_NOT_IMPLEMENTED: Function not yet implemented.
     !!
     !! \return error code

     INTEGER*8 FUNCTION h5_get_error_number ()
     END FUNCTION h5_get_error_number

     !> @}
  END INTERFACE
