!
!  Copyright (c) 2006-2016, The Regents of the University of California,
!  through Lawrence Berkeley National Laboratory (subject to receipt of any
!  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
!  Institut (Switzerland).  All rights reserved.!
!
!  License: see file COPYING in top level of source distribution.
!
  INTERFACE
     !> \ingroup h5part_f90_api
     !! \addtogroup h5part_io_f
     !! @{

     !>
     !! See \ref H5PartWriteDataFloat64
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5pt_writedata_r8 ( filehandle, name, data )
       INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
       CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
       REAL*8, INTENT(IN) :: data(*)               !< the array of float64 data to write
     END FUNCTION h5pt_writedata_r8

     !>
     !! See \ref H5PartWriteDataFloat32
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5pt_writedata_r4 ( filehandle, name, data )
       INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
       CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
       REAL, INTENT(IN) :: data(*)                 !< the array of float32 data to write
     END FUNCTION h5pt_writedata_r4

     !>
     !! See \ref H5PartWriteDataInt64
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5pt_writedata_i8 ( filehandle, name, data )
       INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
       CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
       INTEGER*8, INTENT(IN) :: data(*)            !< the array of int64 data to write
     END FUNCTION h5pt_writedata_i8

     !>
     !! See \ref H5PartWriteDataInt32
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5pt_writedata_i4 ( filehandle, name, data )
       INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
       CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
       INTEGER, INTENT(IN) :: data(*)              !< the array of int32 data to write
     END FUNCTION h5pt_writedata_i4


     !>
     !! See \ref H5PartReadDataFloat64
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5pt_readdata_r8 (filehandle,name,data)
       INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
       CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
       REAL*8, INTENT(OUT) :: data(*)              !< array to read float64 data into
     END FUNCTION h5pt_readdata_r8

     !>
     !! See \ref H5PartReadDataFloat32
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5pt_readdata_r4 (filehandle,name,data)
       INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
       CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
       REAL, INTENT(OUT) :: data(*)                !< array to read float32 data into
     END FUNCTION h5pt_readdata_r4

     !>
     !! See \ref H5PartReadDataInt64
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5pt_readdata_i8 (filehandle,name,data)
       INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
       CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
       INTEGER*8, INTENT(OUT) :: data(*)           !< array to read int64 data into
     END FUNCTION h5pt_readdata_i8

     !>
     !! See \ref H5PartReadDataInt32
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5pt_readdata_i4 (filehandle,name,data)
       INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
       CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
       INTEGER, INTENT(OUT) :: data(*)             !< array to read int32 data into
     END FUNCTION h5pt_readdata_i4

     !> @}
     
  END INTERFACE
