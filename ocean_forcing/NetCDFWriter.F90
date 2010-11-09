!    Copyright (C) 2006 Imperial College London and others.
!
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "confdefs.h"
#include "fdebug.h"

! Tries to be CF 1.4 complient
module NetCDFWriter
  private
  public::NetCDFWriter_init, NetCDFWriter_write_variable
  
  integer, save::ncid
  integer, save::latitude_dim, longitude_dim, time_dim
  character(len=4096), save::ncfilename
  integer(selected_int_kind(3)), save::fillvalue=32767


contains
  subroutine NetCDFWriter_init(filename, longitude, latitude, &
       time_coord, time_units, &
       title, institution, history, source, references, comment)
    use FLDebug
    implicit none
    character(len=*), intent(in)::filename
    real, intent(in)::longitude(:), latitude(:)
    integer, optional, intent(in)::time_coord(:)
    character(len=*), optional, intent(in)::time_units

    ! A succinct description of what is in the dataset.
    character(len=*), optional, intent(in)::title

    ! Specifies where the original data was produced.
    character(len=*), optional, intent(in)::institution

    ! The method of production of the original data. If it was
    ! model-generated, source should name the model and its version,
    ! as specifically as could be useful. If it is observational,
    ! source should characterize it (e.g., "surface observation" or
    ! "radiosonde").
    character(len=*), optional, intent(in)::history

    ! Provides an audit trail for modifications to the original
    ! data. Well-behaved generic netCDF filters will automatically
    ! append their name and the parameters with which they were
    ! invoked to the global history attribute of an input netCDF
    ! file. We recommend that each line begin with a timestamp
    ! indicating the date and time of day that the program was
    ! executed.
    character(len=*), optional, intent(in)::source

    
    ! Published or web-based references that describe the data or
    ! methods used to produce it.
    character(len=*), optional, intent(in)::references
    
    ! Miscellaneous information about the data or methods used to produce it.     
    character(len=*), optional, intent(in)::comment

#ifdef HAVE_NETCDF
    include 'netcdf.inc'

    character(len=256)::conventions="CF-1.4"
    character(len=256)::units_latitude="degrees_north"
    character(len=256)::units_longitude="degrees_east"
#ifdef DOUBLEP
    integer, parameter::ncreal=ncdouble
#else
    integer, parameter::ncreal=ncfloat
#endif
    integer varid, dims(1)
    integer ierr
    real default_time(1)

    default_time(1) = 0

    ! Open up a netcdf file for writing.
    ncfilename = filename
    ncid = nccre(filename, NCCLOB, ierr)
    assert(ierr.eq.0)
    
    call ncaptc(ncid, NF_GLOBAL, "Conventions", NCCHAR, &
         len_trim(conventions), conventions, ierr)
    assert(ierr.eq.0)

    if(present(title)) then
       call ncaptc(ncid, NF_GLOBAL, "title", NCCHAR, &
            len_trim(title), title, ierr)
       assert(ierr.eq.0)
    end if

    if(present(institution)) then
       call ncaptc(ncid, NF_GLOBAL, "institution", NCCHAR, &
            len_trim(institution), institution, ierr)
       assert(ierr.eq.0)
    end if

    if(present(history)) then
       call ncaptc(ncid, NF_GLOBAL, "history", NCCHAR, &
            len_trim(history), history, ierr)
       assert(ierr.eq.0)
    end if

    if(present(source)) then
       call ncaptc(ncid, NF_GLOBAL, "source", NCCHAR, &
            len_trim(source), source, ierr)
       assert(ierr.eq.0)
    end if

    if(present(references)) then
       call ncaptc(ncid, NF_GLOBAL, "references", NCCHAR, &
            len_trim(references), references, ierr)
       assert(ierr.eq.0)
    end if

    if(present(comment)) then
       call ncaptc(ncid, NF_GLOBAL, "comment", NCCHAR, &
            len_trim(comment), comment, ierr)
       assert(ierr.eq.0)
    end if
    
    time_dim = ncddef(ncid, 'time', NCUNLIM, ierr)
    assert(ierr.eq.0)
    
    dims(1) = time_dim
    varid = ncvdef(ncid, "time", ncreal, 1, dims, ierr)
    assert(ierr.eq.0)
    
    call ncaptc(ncid, varid, "long_name", NCCHAR, 4, "time", ierr)
    assert(ierr.eq.0)
    
    ! Variables representing time must always explicitly include
    ! the units attribute; there is no default value. The units
    ! attribute takes a string value formatted as per the
    ! recommendations in the Udunits package [UDUNITS].
    if(present(time_units)) then
       call ncaptc(ncid, varid, "units", NCCHAR, &
            len_trim(time_units), time_units, ierr)
    else
       call ncaptc(ncid, varid, "units", NCCHAR, &
            27, "seconds since 0-01-01 0:0:0", ierr)
    end if
    assert(ierr.eq.0)
    
    ! Define the latitude dimension
    latitude_dim = ncddef(ncid, 'latitude', size(latitude), ierr)
    assert(ierr.eq.0)
    
    dims(1) = latitude_dim
    varid = ncvdef(ncid, "latitude", ncreal, 1, dims, ierr)
    assert(ierr.eq.0)
    
    call ncaptc(ncid, varid, "long_name", NCCHAR, 8, "latitude", ierr)
    assert(ierr.eq.0)
    
    call ncaptc(ncid, varid, "units", NCCHAR, len_trim(units_latitude), units_latitude, ierr)
    assert(ierr.eq.0)

    ! Define the longitude dimension
    longitude_dim = ncddef(ncid, 'longitude', size(longitude), ierr)
    assert(ierr.eq.0)
    
    dims(1) = longitude_dim
    varid = ncvdef(ncid, "longitude", ncreal, 1, dims, ierr)
    assert(ierr.eq.0)
    
    call ncaptc(ncid, varid, "long_name", NCCHAR, 9, "longitude", ierr)
    assert(ierr.eq.0)
    
    call ncaptc(ncid, varid, "units", NCCHAR, len_trim(units_longitude), units_longitude, ierr)
    assert(ierr.eq.0)

    ! Finish defining metadata
    call ncendf(ncid, ierr)
    assert(ierr.eq.0)

    ! Write variable values
    varid = ncvid(ncid, "time", ierr)
    if(present(time_coord)) then
       call ncvpt(ncid, varid, 1, size(time_coord), time_coord, ierr)
    else
       call ncvpt(ncid, varid, 1, 1, default_time, ierr)
    end if

    varid = ncvid(ncid, "latitude", ierr)
    call ncvpt(ncid, varid, 1, size(latitude), latitude, ierr)
    
    varid = ncvid(ncid, "longitude", ierr)
    call ncvpt(ncid, varid, 1, size(longitude), longitude, ierr)
    
    call ncclos(ncid, ierr)
    assert(ierr.eq.0)

#endif
    return
  end subroutine NetCDFWriter_init
  
  subroutine NetCDFWriter_write_variable(name, long_name, variable, units)
    use FLDebug
    implicit none
    character(len=*), intent(in)::name, long_name, units
    real, intent(in)::variable(:,:,:)
#ifdef HAVE_NETCDF
    integer varid

    ! dimension IDs
    integer dims(3), start(3), count(3)
    integer(selected_int_kind(3)), allocatable::short_var(:,:,:)
    integer i, j, k
    real max_v, min_v
    real add_offset, scale_factor
    integer ierr

    include 'netcdf.inc'
#ifdef DOUBLEP
    integer, parameter::ncreal=ncdouble
#else
    integer, parameter::ncreal=ncfloat
#endif

    ! reopen file
    ncid = ncopn(ncfilename, NCWRITE, ierr)

    ! Put file into define mode
    call ncredf(ncid, ierr)

    ! Define the variable
    dims(1) = longitude_dim
    dims(2) = latitude_dim
    dims(3) = time_dim

    varid = ncvdef(ncid, name, ncshort, 3, dims, ierr)
    assert(ierr.eq.0)
    
    varid = ncvid(ncid, name, ierr)
    call ncaptc(ncid, varid, "long_name", NCCHAR, &
         len_trim(long_name), long_name, ierr)
    assert(ierr.eq.0)
    
    call ncaptc(ncid, varid, "units", NCCHAR, &
         len_trim(units), units, ierr)
    assert(ierr.eq.0)
    
    call ncapt(ncid, varid, "_FillValue", NCSHORT, &
         1, fillvalue, ierr)
    assert(ierr.eq.0)
    
    call ncapt(ncid, varid, "missing_value", NCSHORT, &
         1, fillvalue, ierr)
    assert(ierr.eq.0)

    ! Packing algorithm
    max_v = maxval(variable)
    min_v = minval(variable)
    
    add_offset = (max_v + min_v)*0.5
    call ncapt(ncid, varid, "add_offset", ncdouble, &
         1, add_offset, ierr)
    assert(ierr.eq.0)
    
    scale_factor = (max_v - min_v)/(2**16 - 5)
    call ncapt(ncid, varid, "scale_factor", ncdouble, &
         1, scale_factor, ierr)
    assert(ierr.eq.0)
    
    ! Finish defining metadata
    call ncendf(ncid, ierr)
    assert(ierr.eq.0)
    
    start = 1
    count(1) = size(variable(:, 1, 1))
    count(2) = size(variable(1, :, 1))
    count(3) = size(variable(1, 1, :))
    allocate(short_var(count(1), count(2), count(3)))

    do i=1, count(1)
       do j=1, count(2)
          do k=1, count(3)
             short_var(i, j, k) = nint((variable(i, j, k) - add_offset)/scale_factor)
             assert(short_var(i, j, k).ne.fillvalue)
          end do
       end do
    end do
    call ncvpt(ncid, varid, start, count, short_var, ierr)

    call ncclos(ncid, ierr)
    assert(ierr.eq.0)
#endif
  end subroutine NetCDFWriter_write_variable
  
end module NetCDFWriter
