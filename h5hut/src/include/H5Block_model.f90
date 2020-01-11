!
!  Copyright (c) 2006-2016, The Regents of the University of California,
!  through Lawrence Berkeley National Laboratory (subject to receipt of any
!  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
!  Institut (Switzerland).  All rights reserved.!
!
!  License: see file COPYING in top level of source distribution.
!
  INTERFACE
     !> \ingroup h5block_f90_api
     !! \addtogroup h5block_model_f
     !! @{
     
     !      _       _                              _      _ 
     !   __| | __ _| |_ __ _   _ __ ___   ___   __| | ___| |
     !  / _` |/ _` | __/ _` | | '_ ` _ \ / _ \ / _` |/ _ \ |
     ! | (_| | (_| | || (_| | | | | | | | (_) | (_| |  __/ |
     !  \__,_|\__,_|\__\__,_| |_| |_| |_|\___/ \__,_|\___|_|


     !>
     !! See \ref H5Block3dSetView()
     !! \return 0 on success or error code
 
     INTEGER*8 FUNCTION h5bl_3d_setview ( filehandle, i_start, i_end, j_start, j_end, k_start, k_end )
       INTEGER*8, INTENT(IN) :: filehandle
       INTEGER*8, INTENT(IN) :: i_start
       INTEGER*8, INTENT(IN) :: i_end
       INTEGER*8, INTENT(IN) :: j_start
       INTEGER*8, INTENT(IN) :: j_end
       INTEGER*8, INTENT(IN) :: k_start
       INTEGER*8, INTENT(IN) :: k_end
     END FUNCTION h5bl_3d_setview

     !>
     !! See \ref H5Block3dGetView()
     !! \return 0 on success or error code
 
     INTEGER*8 FUNCTION h5bl_3d_getview ( filehandle, i_start, i_end, j_start, j_end, k_start, k_end )
       INTEGER*8, INTENT(IN) :: filehandle
       INTEGER*8, INTENT(OUT) :: i_start
       INTEGER*8, INTENT(OUT) :: i_end
       INTEGER*8, INTENT(OUT) :: j_start
       INTEGER*8, INTENT(OUT) :: j_end
       INTEGER*8, INTENT(OUT) :: k_start
       INTEGER*8, INTENT(OUT) :: k_end
     END FUNCTION h5bl_3d_getview

     !>
     !! See \ref H5Block3dGetReducedView()
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_getreducedview ( filehandle, i_start, i_end, j_start, j_end, k_start, k_end )
       INTEGER*8, INTENT(IN) :: filehandle
       INTEGER*8, INTENT(OUT) :: i_start
       INTEGER*8, INTENT(OUT) :: i_end
       INTEGER*8, INTENT(OUT) :: j_start
       INTEGER*8, INTENT(OUT) :: j_end
       INTEGER*8, INTENT(OUT) :: k_start
       INTEGER*8, INTENT(OUT) :: k_end
     END FUNCTION h5bl_3d_getreducedview

     !>
     !! See \ref H5Block3dHasView()
     !! \return rank of processor error code

     INTEGER*8 FUNCTION h5bl_3d_hasview ( filehandle )
       INTEGER*8, INTENT(IN) :: filehandle
     END FUNCTION h5bl_3d_hasview

     !>
     !! See \ref H5Block3dSetChunkSize()
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_3d_setchunk ( filehandle, i, j, k )
       INTEGER*8, INTENT(IN) :: filehandle
       INTEGER*8, INTENT(IN) :: i
       INTEGER*8, INTENT(IN) :: j
       INTEGER*8, INTENT(IN) :: k
     END FUNCTION h5bl_3d_setchunk

     !>
     !! See \ref H5BlockGetNumFields
     !! \return number of fields or error code

     INTEGER*8 FUNCTION h5bl_getnumfields ( filehandle )
       INTEGER*8, INTENT(IN) :: filehandle
     END FUNCTION h5bl_getnumfields

     !>
     !! See \ref H5BlockGetFieldInfo
     !! \return 0 on success or error code

     INTEGER*8 FUNCTION h5bl_getfieldinfo ( filehandle, idx, field_name, grid_rank, grid_dims, field_dims )
       INTEGER*8, INTENT(IN) :: filehandle
       INTEGER*8, INTENT(IN) :: idx
       CHARACTER(LEN=*), INTENT(OUT) :: field_name
       INTEGER*8, INTENT(OUT) :: grid_rank
       INTEGER*8, INTENT(OUT) :: grid_dims(*)
       INTEGER*8, INTENT(OUT) :: field_dims
     END FUNCTION h5bl_getfieldinfo

     !> @}
     
  END INTERFACE
