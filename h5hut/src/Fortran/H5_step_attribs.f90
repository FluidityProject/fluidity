!
!  Copyright (c) 2006-2016, The Regents of the University of California,
!  through Lawrence Berkeley National Laboratory (subject to receipt of any
!  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
!  Institut (Switzerland).  All rights reserved.!
!
!  License: see file COPYING in top level of source distribution.
!
  INTERFACE

     !> \addtogroup \h5_step_attribs_f
     !! @{

     !      _                     _   _        _ _           _            
     !  ___| |_ ___ _ __     __ _| |_| |_ _ __(_) |__  _   _| |_ ___  ___ 
     ! / __| __/ _ \ '_ \   / _` | __| __| '__| | '_ \| | | | __/ _ \/ __|
     ! \__ \ ||  __/ |_) | | (_| | |_| |_| |  | | |_) | |_| | ||  __/\__ \
     ! |___/\__\___| .__/   \__,_|\__|\__|_|  |_|_.__/ \__,_|\__\___||___/
     !             |_|
     !   __ _ _   _  ___ _ __ _   _ 
     !  / _` | | | |/ _ \ '__| | | |
     ! | (_| | |_| |  __/ |  | |_| |
     !  \__, |\__,_|\___|_|   \__, |
     !     |_|                |___/

     !>
     !! See \ref H5GetNumFileAttribs
     !! \return number of attributes or error code

     INTEGER*8 FUNCTION h5_getnstepattribs (filehandle)
       INTEGER*8, INTENT(IN) :: filehandle         !< file handle
     END FUNCTION h5_getnstepattribs

     !>
     !! See \ref H5GetFileAttribInfo
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5_getstepattribinfo (filehandle, idx, name, type, nelem)
       INTEGER*8,INTENT(IN) :: filehandle      !< file handle
       INTEGER*8,INTENT(IN) :: idx             !< index of attribute being queried
       CHARACTER(LEN=*), INTENT(OUT):: name    !< name of attribute
       INTEGER*8,INTENT(OUT):: type            !< type of attribute
       INTEGER*8,INTENT(OUT):: nelem           !< number of elements in the attrib array
     END FUNCTION h5_getstepattribinfo

     !>
     !! See \ref H5GetStepAttribInfoByName
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5_getstepattribinfo_by_name (fhandle, name, type, nelem)
       INTEGER*8,INTENT(IN) :: fhandle         !< file handle
       CHARACTER(LEN=*), INTENT(IN):: name     !< name of attribute
       INTEGER*8,INTENT(OUT):: type            !< type of attribute
       INTEGER*8,INTENT(OUT):: nelem           !< number of elements in the attrib array
     END FUNCTION h5_getstepattribinfo_by_name

     !      _                     _   _        _ _           _            
     !  ___| |_ ___ _ __     __ _| |_| |_ _ __(_) |__  _   _| |_ ___  ___ 
     ! / __| __/ _ \ '_ \   / _` | __| __| '__| | '_ \| | | | __/ _ \/ __|
     ! \__ \ ||  __/ |_) | | (_| | |_| |_| |  | | |_) | |_| | ||  __/\__ \
     ! |___/\__\___| .__/   \__,_|\__|\__|_|  |_|_.__/ \__,_|\__\___||___/
     !             |_|
     !  _    __          _        _             
     ! (_)  / /__    ___| |_ _ __(_)_ __   __ _ 
     ! | | / / _ \  / __| __| '__| | '_ \ / _` |
     ! | |/ / (_) | \__ \ |_| |  | | | | | (_| |
     ! |_/_/ \___/  |___/\__|_|  |_|_| |_|\__, |
     !                                    |___/

     !>
     !! See \ref H5WriteStepAttribString
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5_writestepattrib_string (filehandle, attrib_name, attrib_value)
       INTEGER*8, INTENT(IN) :: filehandle         !< file handle
       CHARACTER(LEN=*), INTENT(IN) :: attrib_name !< name of attribute to write
       CHARACTER(LEN=*), INTENT(IN) :: attrib_value!< attribute data to be written
     END FUNCTION h5_writestepattrib_string

     !>
     !! See \ref H5ReadStepAttribString
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5_readstepattrib_string (filehandle, attrib_name, attrib_value)
       INTEGER*8, INTENT(IN) :: filehandle         !< file handle
       CHARACTER(LEN=*), INTENT(IN) :: attrib_name !< name of the attribute to read
       CHARACTER(LEN=*), INTENT(IN) :: attrib_value!< attribute data will be read into this array
     END FUNCTION h5_readstepattrib_string

     !      _                     _   _        _ _           _            
     !  ___| |_ ___ _ __     __ _| |_| |_ _ __(_) |__  _   _| |_ ___  ___ 
     ! / __| __/ _ \ '_ \   / _` | __| __| '__| | '_ \| | | | __/ _ \/ __|
     ! \__ \ ||  __/ |_) | | (_| | |_| |_| |  | | |_) | |_| | ||  __/\__ \
     ! |___/\__\___| .__/   \__,_|\__|\__|_|  |_|_.__/ \__,_|\__\___||___/
     !             |_|
     !  _    __                     _ 
     ! (_)  / /__    _ __ ___  __ _| |
     ! | | / / _ \  | '__/ _ \/ _` | |
     ! | |/ / (_) | | | |  __/ (_| | |
     ! |_/_/ \___/  |_|  \___|\__,_|_|

     !>
     !! See \ref H5WriteStepAttribFloat64
     !! \return 0 on success or error code
     !<
     INTEGER*8 FUNCTION h5_writestepattrib_r8 (filehandle, attrib_name, attrib_value, attrib_nelem)
       INTEGER*8,INTENT(IN) :: filehandle          !< file handle
       CHARACTER(LEN=*), INTENT(IN) :: attrib_name !< name of attribute to read
       REAL*8,   INTENT(IN):: attrib_value(*)      !< attribute data to be written
       INTEGER*8, INTENT(IN) :: attrib_nelem       !< number of elements in data array
     END FUNCTION h5_writestepattrib_r8

     !>
     !! See \ref H5ReadStepAttribFloat64
     !! \return 0 on success or error code
     !<
     INTEGER*8 FUNCTION h5_readstepattrib_r8 (filehandle, attrib_name, attrib_value)
       INTEGER*8,INTENT(IN) :: filehandle          !< file handle
       CHARACTER(LEN=*), INTENT(IN) :: attrib_name !< name of attribute to read
       REAL*8,   INTENT(OUT):: attrib_value(*)     !< attribute data will be read into this array
     END FUNCTION h5_readstepattrib_r8

     !>
     !! See \ref H5WriteStepAttribFloat32
     !! \return 0 on success or error code
     !<
     INTEGER*8 FUNCTION h5_writestepattrib_r4 (filehandle, attrib_name, attrib_value, attrib_nelem)
       INTEGER*8,INTENT(IN) :: filehandle          !< file handle
       CHARACTER(LEN=*), INTENT(IN) :: attrib_name !< name of attribute to read
       REAL*4,   INTENT(IN):: attrib_value(*)      !< attribute data to be written
       INTEGER*8, INTENT(IN) :: attrib_nelem       !< number of elements in data array
     END FUNCTION h5_writestepattrib_r4

     !>
     !! See \ref H5ReadStepAttribFloat32
     !! \return 0 on success or error code
     !<
     INTEGER*8 FUNCTION h5_readstepattrib_r4 ( filehandle, attrib_name, attrib_value )
       INTEGER*8,INTENT(IN) :: filehandle          !< file handle
       CHARACTER(LEN=*), INTENT(IN) :: attrib_name !< name of attribute to read
       REAL*4,   INTENT(OUT):: attrib_value(*)     !< attribute data will be read into this array
     END FUNCTION h5_readstepattrib_r4

     !      _                     _   _        _ _           _            
     !  ___| |_ ___ _ __     __ _| |_| |_ _ __(_) |__  _   _| |_ ___  ___ 
     ! / __| __/ _ \ '_ \   / _` | __| __| '__| | '_ \| | | | __/ _ \/ __|
     ! \__ \ ||  __/ |_) | | (_| | |_| |_| |  | | |_) | |_| | ||  __/\__ \
     ! |___/\__\___| .__/   \__,_|\__|\__|_|  |_|_.__/ \__,_|\__\___||___/
     !             |_|
     !  _    __      _       _                       
     ! (_)  / /__   (_)_ __ | |_ ___  __ _  ___ _ __ 
     ! | | / / _ \  | | '_ \| __/ _ \/ _` |/ _ \ '__|
     ! | |/ / (_) | | | | | | ||  __/ (_| |  __/ |   
     ! |_/_/ \___/  |_|_| |_|\__\___|\__, |\___|_|   
     !                               |___/

     !>
     !! See \ref H5WriteStepAttribInt64
     !! \return 0 on success or error code
     !<
     INTEGER*8 FUNCTION h5_writestepattrib_i8 (filehandle, attrib_name, attrib_value, attrib_nelem)
       INTEGER*8,INTENT(IN) :: filehandle          !< file handle
       CHARACTER(LEN=*), INTENT(IN) :: attrib_name !< name of attribute to write
       INTEGER*8,INTENT(IN):: attrib_value(*)     !< attribute data to be written
       INTEGER*8, INTENT(IN) :: attrib_nelem       !< number of elements in data array
     END FUNCTION h5_writestepattrib_i8

     !>
     !! See \ref H5ReadStepAttribInt64
     !! \return 0 on success or error code
     !<
     INTEGER*8 FUNCTION h5_readstepattrib_i8 (filehandle, attrib_name, attrib_value)
       INTEGER*8,INTENT(IN) :: filehandle          !< file handle
       CHARACTER(LEN=*), INTENT(IN) :: attrib_name !< name of attribute to read
       INTEGER*8,INTENT(OUT):: attrib_value(*)     !< attribute data will be read into this array
     END FUNCTION h5_readstepattrib_i8

     !>
     !! See \ref H5WriteStepAttribInt32
     !! \return 0 on success or error code
     !<
     INTEGER*8 FUNCTION h5_writestepattrib_i4 (filehandle, attrib_name, attrib_value, attrib_nelem)
       INTEGER*8,INTENT(IN) :: filehandle          !< file handle
       CHARACTER(LEN=*), INTENT(IN) :: attrib_name !< name of attribute to write
       INTEGER*4,INTENT(IN):: attrib_value(*)      !< attribute data to be written
       INTEGER*8, INTENT(IN) :: attrib_nelem       !< number of elements in data array
     END FUNCTION h5_writestepattrib_i4

     !>
     !! See \ref H5ReadStepAttribInt32
     !! \return 0 on success or error code
     !<
     INTEGER*8 FUNCTION h5_readstepattrib_i4 (filehandle, attrib_name, attrib_value)
       INTEGER*8,INTENT(IN) :: filehandle          !< file handle
       CHARACTER(LEN=*), INTENT(IN) :: attrib_name !< name of attribute to read
       INTEGER*4,INTENT(OUT):: attrib_value(*)     !< attribute data will be read into this array
     END FUNCTION h5_readstepattrib_i4

     !> @}
  END INTERFACE
