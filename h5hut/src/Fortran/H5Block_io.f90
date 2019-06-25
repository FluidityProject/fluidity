!
!  Copyright (c) 2006-2016, The Regents of the University of California,
!  through Lawrence Berkeley National Laboratory (subject to receipt of any
!  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
!  Institut (Switzerland).  All rights reserved.!
!
!  License: see file COPYING in top level of source distribution.
!
  INTERFACE
     !! \addtogroup h5block_io_f
     !! @{

     !>
     !! Write the 3-dimensional field \p name from the array \p buffer
     !! to the current step using the previously defined field
     !! view.
     !!
     !! The data type is 64bit floating point (\c REAL*8). Ensure that
     !! the number of items in the buffer matches the view.
     !!
     !! \return \c H5_SUCCESS on success
     !! \return \c H5_FAILURE on error

     INTEGER*8 FUNCTION h5bl_3d_write_scalar_field_r8 ( filehandle, name, buffer )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       REAL*8, INTENT(IN) :: buffer(*)    !< the array of data
     END FUNCTION h5bl_3d_write_scalar_field_r8

     !>
     !! Read the 3-dimensional field \c name into the array \p buffer
     !!from the current step using the previously defined field layout.
     !!
     !! The data type is 64bit floating point (\c REAL*8). Ensure that
     !! the number of items in the buffer matches the view.
     !!
     !! \return \c H5_SUCCESS on success
     !! \return \c H5_FAILURE on error

     INTEGER*8 FUNCTION h5bl_3d_read_scalar_field_r8 ( filehandle, name, buffer )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       REAL*8, INTENT(OUT) :: buffer(*)   !< buffer to read the data into
     END FUNCTION h5bl_3d_read_scalar_field_r8

     !>
     !! See \ref H5Block3dWriteVector3dFieldFloat64
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_write_vector3d_field_r8 ( filehandle, name, x, y, z )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       REAL*8, INTENT(IN) :: x(*)       !< the array of x data to write
       REAL*8, INTENT(IN) :: y(*)       !< the array of y data to write
       REAL*8, INTENT(IN) :: z(*)       !< the array of z data to write
     END FUNCTION h5bl_3d_write_vector3d_field_r8

     !>
     !! See \ref H5Block3dReadVector3dFieldFloat64
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_read_vector3d_field_r8 ( filehandle, name, x, y, z )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       REAL*8, INTENT(OUT) :: x(*)      !< buffer to read the x data into
       REAL*8, INTENT(OUT) :: y(*)      !< buffer to read the y data into
       REAL*8, INTENT(OUT) :: z(*)      !< buffer to read the z data into
     END FUNCTION h5bl_3d_read_vector3d_field_r8

     !>
     !! See \ref H5Block3dWriteScalarFieldFloat32
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_write_scalar_field_r4 ( filehandle, name, buffer )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       REAL*4, INTENT(IN) :: buffer(*)    !< the array of data
     END FUNCTION h5bl_3d_write_scalar_field_r4

     !>
     !! See \ref H5Block3dReadScalarFieldFloat32
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_read_scalar_field_r4 ( filehandle, name, buffer )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       REAL*4, INTENT(OUT) :: buffer(*)   !< buffer to read the data into
     END FUNCTION h5bl_3d_read_scalar_field_r4

     !>
     !! See \ref H5Block3dWriteVector3dFieldFloat32
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_write_vector3d_field_r4 ( filehandle, name, x, y, z )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       REAL*4, INTENT(IN) :: x(*)       !< the array of x data to write
       REAL*4, INTENT(IN) :: y(*)       !< the array of y data to write
       REAL*4, INTENT(IN) :: z(*)       !< the array of z data to write
     END FUNCTION h5bl_3d_write_vector3d_field_r4

     !>
     !! See \ref H5Block3dReadVector3dFieldFloat32
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_read_vector3d_field_r4 ( filehandle, name, x, y, z )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       REAL*4, INTENT(OUT) :: x(*)      !< buffer to read the x data into
       REAL*4, INTENT(OUT) :: y(*)      !< buffer to read the y data into
       REAL*4, INTENT(OUT) :: z(*)      !< buffer to read the z data into
     END FUNCTION h5bl_3d_read_vector3d_field_r4

     !>
     !! See \ref H5Block3dWriteScalarFieldInt64
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_write_scalar_field_i8 ( filehandle, name, buffer )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       INTEGER*8, INTENT(IN) :: buffer(*)    !< the array of data
     END FUNCTION h5bl_3d_write_scalar_field_i8

     !>
     !! See \ref H5Block3dReadScalarFieldInt64
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_read_scalar_field_i8 ( filehandle, name, buffer )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       INTEGER*8, INTENT(OUT) :: buffer(*)   !< buffer to read the data into
     END FUNCTION h5bl_3d_read_scalar_field_i8

     !>
     !! See \ref H5Block3dWriteVector3dFieldInt64
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_write_vector3d_field_i8 ( filehandle, name, x, y, z )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       INTEGER*8, INTENT(IN) :: x(*)       !< the array of x data to write
       INTEGER*8, INTENT(IN) :: y(*)       !< the array of y data to write
       INTEGER*8, INTENT(IN) :: z(*)       !< the array of z data to write
     END FUNCTION h5bl_3d_write_vector3d_field_i8

     !>
     !! See \ref H5Block3dReadVector3dFieldInt64
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_read_vector3d_field_i8 ( filehandle, name, x, y, z )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       INTEGER*8, INTENT(OUT) :: x(*)      !< buffer to read the x data into
       INTEGER*8, INTENT(OUT) :: y(*)      !< buffer to read the y data into
       INTEGER*8, INTENT(OUT) :: z(*)      !< buffer to read the z data into
     END FUNCTION h5bl_3d_read_vector3d_field_i8

     !>
     !! See \ref H5Block3dWriteScalarFieldInt32
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_write_scalar_field_i4 ( filehandle, name, buffer )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       INTEGER*4, INTENT(IN) :: buffer(*)    !< the array of data
     END FUNCTION h5bl_3d_write_scalar_field_i4

     !>
     !! See \ref H5Block3dReadScalarFieldInt32
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_read_scalar_field_i4 ( filehandle, name, buffer )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       INTEGER*4, INTENT(OUT) :: buffer(*)   !< buffer to read the data into
     END FUNCTION h5bl_3d_read_scalar_field_i4

     !>
     !! See \ref H5Block3dWriteVector3dFieldInt32
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_write_vector3d_field_i4 ( filehandle, name, x, y, z )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       INTEGER*4, INTENT(IN) :: x(*)       !< the array of x data to write
       INTEGER*4, INTENT(IN) :: y(*)       !< the array of y data to write
       INTEGER*4, INTENT(IN) :: z(*)       !< the array of z data to write
     END FUNCTION h5bl_3d_write_vector3d_field_i4

     !>
     !! See \ref H5Block3dReadVector3dFieldInt32
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_read_vector3d_field_i4 ( filehandle, name, x, y, z )
       INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
       CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
       INTEGER*4, INTENT(OUT) :: x(*)      !< buffer to read the x data into
       INTEGER*4, INTENT(OUT) :: y(*)      !< buffer to read the y data into
       INTEGER*4, INTENT(OUT) :: z(*)      !< buffer to read the z data into
     END FUNCTION h5bl_3d_read_vector3d_field_i4

     !> @}
  END INTERFACE
