!
!  Copyright (c) 2006-2016, The Regents of the University of California,
!  through Lawrence Berkeley National Laboratory (subject to receipt of any
!  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
!  Institut (Switzerland).  All rights reserved.!
!
!  License: see file COPYING in top level of source distribution.
!
  INTERFACE
     !> \addtogroup h5_file_f
     !! @{
     
     !>
     !! Create a new, empty file property list.
     !!
     !! File property lists are used to control optional behavior like file
     !! creation, file access, dataset creation, dataset transfer. File
     !! property lists are attached to file handles while opened with \ref
     !! h5_openfile().
     !!
     !! \return empty file property list
     !! \return \c H5_FAILURE on error
     !!
     !! \see h5_setprop_file_mpio()
     !! \see h5_setprop_file_mpio_collective()
     !! \see h5_setprop_file_mpio_independent()
     !! \see h5_setprop_file_mpio_posix() (HDF5 <= 1.8.12 only)
     !! \see h5_setprop_file_corevfd()
     !! \see h5_setprop_file_align()
     !! \see h5_setprop_file_throttle()

     INTEGER*8 FUNCTION h5_createprop_file ()
     END FUNCTION h5_createprop_file

     !>
     !! Stores MPI IO communicator information to given file property list. If used in 
     !! \ref h5_openfile(), MPI collective IO will be used.
     !!
     !! \note This function deprecated. Use h5_setprop_file_mpio_collective() instead.
     !!
     !! \return \c H5_SUCCESS on success
     !! \return \c H5_FAILURE on error
     !!
     !! \see h5_setprop_file_mpio_independent()
     !! \see h5_setprop_file_mpio_posix() (HDF5 <= 1.8.12 only)
     !! \see h5_setprop_file_corevfd()

     INTEGER*8 FUNCTION h5_setprop_file_mpio (prop, comm)
       INTEGER*8, INTENT(IN) :: prop               !< property
       INTEGER, INTENT(IN) :: comm                 !< the MPI communicator used by the program
     END FUNCTION h5_setprop_file_mpio

     !>
     !! Stores MPI IO communicator information to given file property list. If used in
     !! \ref h5_openfile(), MPI collective IO will be used.
     !!
     !! \return \c H5_SUCCESS on success
     !! \return \c H5_FAILURE on error
     !!
     !! \note This function deprecates h5_setprop_file_mpio().
     !!
     !! \see h5_setprop_file_mpio_independent()
     !! \see h5_setprop_file_mpio_posix() (HDF5 <= 1.8.12 only)
     !! \see h5_setprop_file_corevfd()

     INTEGER*8 FUNCTION h5_setprop_file_mpio_collective (prop, comm)
       INTEGER*8, INTENT(IN) :: prop               !< property
       INTEGER, INTENT(IN) :: comm                 !< the MPI communicator used by the program
     END FUNCTION h5_setprop_file_mpio_collective

     !>
     !! Stores MPI IO communicator information to given file property list. If used in 
     !! \ref h5_openfile(), MPI independent IO will be used.
     !!
     !! \return \c H5_SUCCESS on success
     !! \return \c H5_FAILURE on error
     !!
     !! \see h5_setprop_file_mpio_collective()
     !! \see h5_setprop_file_mpio_posix() (HDF5 <= 1.8.12 only)
     !! \see h5_setprop_file_corevfd()

     INTEGER*8 FUNCTION h5_setprop_file_mpio_independent (prop, comm)
       INTEGER*8, INTENT(IN) :: prop               !< property
       INTEGER, INTENT(IN) :: comm                 !< the MPI communicator used by the program
     END FUNCTION h5_setprop_file_mpio_independent

     !>
     !! Stores MPI IO communicator information to given file property list. If used in 
     !! \ref h5_openfile(), MPI POSIX IO will be used.
     !!
     !! \return \c H5_SUCCESS on success
     !! \return \c H5_FAILURE on error
     !!
     !! \note This function is available only, if H5hut has been compiled with
     !! HDF5 1.8.12 or older. 
     !!
     !! \see h5_setprop_file_mpio_collective()
     !! \see h5_setprop_file_mpio_independent()
     !! \see h5_setprop_file_corevfd()

     INTEGER*8 FUNCTION h5_setprop_file_mpio_posix (prop, comm)
       INTEGER*8, INTENT(IN) :: prop               !< property
       INTEGER, INTENT(IN) :: comm                 !< the MPI communicator used by the program
     END FUNCTION h5_setprop_file_mpio_posix

     !>
     !! Modifies the file property list to use the \c H5FD_CORE driver.  The
     !! \c H5FD_CORE driver enables an application to work with a file in memory. 
     !! File contents are stored only in memory until the file is closed.
     !!
     !! The increment by which allocated memory is to be increased each time more
     !! memory is required, must be specified with \c increment.
     !!
     !! \return \c H5_SUCCESS on success
     !! \return \c H5_FAILURE on error

     INTEGER*8 FUNCTION h5_setprop_file_corevfd (prop)
       INTEGER*8, INTENT(IN) :: prop               !< property
     END FUNCTION h5_setprop_file_corevfd

     !>
     !! Sets alignment properties of a file property list so that any file
     !! object greater than or equal in size to threshold bytes will be
     !! aligned on an address which is a multiple of alignment. The
     !! addresses are relative to the end of the user block; the alignment
     !! is calculated by subtracting the user block size from the absolute
     !! file address and then adjusting the address to be a multiple of
     !! alignment.
     !!
     !! Default values for alignment is one, implying no
     !! alignment. Generally the default value result in the best
     !! performance for single-process access to the file. For MPI IO and
     !! other parallel systems, choose an alignment which is a multiple of
     !! the disk block size.
     !!
     !! \return \c H5_SUCCESS on success
     !! \return \c H5_FAILURE on error
     !!
     !! \see h5_setprop_file_corevfd()

     INTEGER*8 FUNCTION h5_setprop_file_align (prop, align)
       INTEGER*8, INTENT(IN) :: prop               !< property
       INTEGER*8, INTENT(IN) :: align              !< alignment
     END FUNCTION h5_setprop_file_align

     !>
     !! Set the `throttle` factor, which causes HDF5 write and read
     !! calls to be issued in that number of batches.
     !!
     !! This can prevent large concurrency parallel applications that
     !! use independent writes from overwhelming the underlying
     !! parallel file system.
     !!
     !! Throttling only works with the H5_VFD_MPIO_POSIX or
     !! H5_VFD_MPIO_INDEPENDENT drivers and is only available in
     !! the parallel library.
     !!
     !! \return \c H5_SUCCESS on success
     !! \return \c H5_FAILURE on error

     INTEGER*8 FUNCTION h5_setprop_file_throttle (prop, throttle)
       INTEGER*8, INTENT(IN) :: prop               !< property
       INTEGER*8, INTENT(IN) :: throttle           !< throttle factor
     END FUNCTION h5_setprop_file_throttle

     !>
     !! Close file property list.

     INTEGER*8 FUNCTION h5_closeprop (prop)
       INTEGER*8, INTENT(IN) :: prop               !< property
     END FUNCTION h5_closeprop

     !>
     !! Open file with name \c filename.
     !!
     !! File mode flags are:
     !! - \c H5_O_RDONLY: Only reading allowed
     !! - \c H5_O_WRONLY: create new file, dataset must not exist
     !! - \c H5_O_APPENDONLY: allows to append new data to an existing file
     !! - \c H5_O_RDWR:   dataset may exist
     !! - \c H5_FS_LUSTRE - enable optimizations for the Lustre file system
     !! - \c H5_VFD_MPIO_POSIX - use the HDF5 MPI-POSIX virtual file driver
     !! - \c H5_VFD_MPIO_INDEPENDENT - use MPI-IO in indepedent mode
     !!
     !! The file is opened with the properties set in the file property list
     !! \c prop.  This argument can also be set to \c H5_PROP_DEFAULT to use
     !! reasonable default values. In this case \c MPI_COMM_WORLD will be
     !! used as MPI communicator in a parallel execution environment.
     !!
     !! The typical file extension is \c .h5.
     !!
     !! \return File handle
     !! \return \c H5_FAILURE on error
     !!
     !! \see h5_createprop_file()

     INTEGER*8 FUNCTION h5_openfile (fname, mode, props)
       CHARACTER(LEN=*), INTENT(IN) :: fname       !< the filename to open for reading
       INTEGER*8, INTENT(IN) :: mode               !< file mode
       INTEGER*8, INTENT(IN) :: props              !< properties
     END FUNCTION h5_openfile

     !>
     !! Close file and free all memory associated with the file handle.
     !!
     !! \return \c H5_SUCCESS on success
     !! \return \c H5_FAILURE on error

     INTEGER*8 FUNCTION h5_closefile (filehandle)
       INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
     END FUNCTION h5_closefile

     !>
     !! Verify that the passed file handle is a valid H5hut file handle.
     !!
     !! \return \c H5_SUCCESS on success
     !! \return \c H5_FAILURE on error

     INTEGER*8 FUNCTION h5_checkfile ( filehandle )
       INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
     END FUNCTION h5_checkfile

     !>
     !! Flush all file data to disk.
     !!
     !! \return \c H5_SUCCESS on success
     !! \return \c H5_FAILURE on error

     INTEGER*8 FUNCTION h5_flushfile (filehandle)
       INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
     END FUNCTION h5_flushfile

     !>
     !! Flush step data to disk.
     !!
     !! \return \c H5_SUCCESS on success
     !! \return \c H5_FAILURE on error

     INTEGER*8 FUNCTION h5_flushstep (filehandle)
       INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
     END FUNCTION h5_flushstep

     !>
     !! Close H5hut library. This function should be called before program exit.
     !!
     !! \return \c H5_SUCCESS on success
     !! \return \c H5_FAILURE on error

     INTEGER*8 FUNCTION h5_finalize ()
     END FUNCTION h5_finalize

     !> @}
  END INTERFACE
