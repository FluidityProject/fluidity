module vtkreader
  use fxmltools
  use xmltree
  use xmlreader

  implicit none
  
  private

  public :: vtk_get_sizes, vtk_read_file
  
 interface
     integer(c_int) function uncompress (dest, destlen, source, sourcelen) bind(c)
       use iso_c_binding
       type(c_ptr), value :: dest
       character(kind=c_char), dimension(*) :: source
       integer(c_long) :: destlen
       integer(c_long), value :: sourcelen
     end function uncompress
  end interface

  contains

    subroutine copy(dest, src)
      character(c_char), dimension(:), allocatable :: dest
      character(c_char), dimension(:) :: src
      allocate(dest(size(src)))
      dest = src
    end subroutine copy

    logical function match (fptr, fstr)
      character, dimension(:), pointer :: fptr
      character(len=*) :: fstr
      integer :: n
      
      match=.false.
      if (size(fptr) .ne. len(fstr)) return 
      do n=1, size(fptr)
         if (fptr(n) .ne. fstr(n:n)) return
      end do
      match = .true.
      
    end function match

    integer function kindpar(datatype)
      character(len=*) :: datatype

      select case(datatype)
      case("Int8")
         kindpar=c_int8_t
      case("Int32")
         kindpar=c_int32_t
      case("Int64")
         kindpar=c_int64_t
      case("Float32")
         kindpar=c_float
      case("Float64")
         kindpar=c_double
      end select

    end function kindpar

    subroutine transfer_integer_data(res, data, dtype)
      integer, dimension(:) :: res
      character, dimension(:), pointer :: data
      integer :: dtype

      select case(dtype)
      case(c_int8_t)
         res = transfer(read_compressed_data(data, size(res)), &
              int(1, kind=c_int8_t), size(res))
      case(c_int32_t)
         res = transfer(read_compressed_data(data, 4*size(res)), &
              int(1, kind=c_int32_t), size(res))
      case(c_int64_t)
         res = transfer(read_compressed_data(data, 8*size(res)), &
              int(1, kind=c_int64_t), size(res))
      end select

    end subroutine transfer_integer_data

    function open_vtk_xml_file(filename, encoding, options) result(reader)
      character(len=*), intent(in) :: filename
      character(len=*), optional, intent(in) :: encoding
      integer, optional, intent(in) :: options
      type(xmlTextReader) :: reader
      integer :: loptions

      if (present(options)) then
         loptions = options
      else
         loptions = 1
      end if

      if (present(encoding)) then
         reader = xmlReaderForFile(filename, encoding, loptions)
      else
         reader = xmlReaderForFile(filename, "utf-8", loptions)
      end if
    end function open_vtk_xml_file

    logical function read_vtk_xml_element(reader, name) result(res)
      type(xmlTextReader), intent(in) :: reader
      character, dimension(:), pointer, optional :: name

      if (present(name)) nullify(name)
      res = .false.
      do while (xmlTextReaderRead(reader) == 1)
         if (xmlTextReaderNodeType(reader) == 1) then
            if (present(name)) name => xmlTextReaderConstName(reader)
            res = .true.
            return
         end if
      end do
    end function read_vtk_xml_element
      

    subroutine read_compressed_header(chardata, metadata, cblocksizes)
      character(kind=c_char), dimension(:), pointer :: chardata
      integer(kind=c_int), intent(out) :: metadata(3)
      integer, intent(out), allocatable :: cblocksizes(:)

      integer :: n, offset,rsize, ret
      type(b64state) :: state_in
      character(kind=c_char), dimension(:), pointer :: cblock=>null()
      character(kind=c_char), dimension(4096) :: comptext

      offset = 0
      rsize = 4*(12/3)
      cblock => chardata(offset+1:offset+rsize)
      offset = offset+rsize
      call base64_init_decodestate(state_in)
      ret = base64_decode_block(cblock, comptext, state_in)
      metadata = transfer(comptext, metadata, size(metadata))

      rsize = 4*((4*metadata(1))/3)
      if (mod(metadata(1),3)>0) rsize = rsize+4
      cblock => chardata(offset+1:min(size(chardata),offset+rsize))
      offset = offset+rsize
      call base64_init_decodestate(state_in)
      ret = base64_decode_block(cblock, comptext, state_in)

      allocate(cblocksizes(metadata(1)))
      cblocksizes = transfer(comptext, cblocksizes, metadata(1))

      chardata => chardata(offset+1:size(chardata))

    end subroutine read_compressed_header

    subroutine read_zlib_compressed_block(chardata, offset, destlen, data, sourcelen, datasize)
      character(kind=c_char), dimension(:), pointer :: chardata
      integer, intent(in) :: offset
      integer(c_long), intent(inout) :: destlen, sourcelen
      integer :: datasize
      type(c_ptr) :: data

      type(b64state) :: state_in
      integer :: ret, lsize, rsize
      character(kind=c_char), dimension(32768) :: comptext
      character(kind=c_char), dimension(:), pointer :: cblock


      lsize=4*(offset/3)
      rsize = 4*((offset+sourcelen)/3)
      if (mod(offset+sourcelen,3)>0) rsize = rsize+4
      cblock => chardata(lsize+1:rsize)
      
      call base64_init_decodestate(state_in)
      ret = base64_decode_block(cblock, comptext, state_in)
      ret = uncompress(data, destlen, comptext(1+mod(offset,3):), sourcelen)

    end subroutine read_zlib_compressed_block

    subroutine get_data_array_properties(reader, name, datatype, ncomps, ntuples)

      type(xmlTextReader), intent(in) :: reader
      character, dimension(:), pointer, intent(out) :: name, datatype
      integer, intent(out) :: ncomps, ntuples

      character, dimension(:), pointer :: tmp
      integer :: ret, i, j

      ret = xmlTextReaderMoveToAttribute(reader, "Name")
      if (ret == 1) then
         name => xmlTextReaderConstValue(reader)
      else
         nullify(name)
      end if
      ret = xmlTextReaderMoveToAttribute(reader, "type")
      if (ret==1) then
         datatype => xmlTextReaderConstValue(reader)
      else
         nullify(name)
      end if
      ret = xmlTextReaderMoveToAttribute(reader, "NumberOfComponents")
      if (ret == 1) then
         tmp => xmlTextReaderConstValue(reader)
         if (match(tmp,"")) then
            ncomps = 1
         else
            ncomps = chararray2int(tmp)
         end if
      else
         ncomps = -1
      end if
      ret = xmlTextReaderMoveToAttribute(reader, "NumberOfTuples")
      if (ret == 1) then
         tmp => xmlTextReaderConstValue(reader)
         ntuples = chararray2int(tmp)
      else
         ntuples = -1
      end if
      
    end subroutine get_data_array_properties
    
    subroutine process_vtk_header(reader, filetype, version, littleendian, header_type, compressed)
      type(xmlTextReader), intent(in) :: reader
      character, dimension(:), pointer, intent(out), optional :: filetype, header_type
      real, intent(out), optional :: version
      logical, intent(out), optional :: littleendian, compressed
      character, dimension(:), pointer :: tmp
      integer :: ret

      if (present(filetype)) then
         ret = xmlTextReaderMoveToAttribute(reader, "type")
         if (ret==1) then
            filetype => xmlTextReaderConstValue(reader)
         else
            nullify(filetype)
         end if
      end if
      if (present(version)) then
         ret = xmlTextReaderMoveToAttribute(reader, "version")
         if (ret==1) then
            tmp => xmlTextReaderConstValue(reader)
            read(tmp, *) version
         else
            version =-1.0
         end if
      end if
      if (present(littleendian)) then
         ret = xmlTextReaderMoveToAttribute(reader, "byte_order")
         littleendian = .true.
         if (ret==1) then
            tmp => xmlTextReaderConstValue(reader)
            littleendian = match(tmp, "LittleEndian")
         end if
      end if
      if (present(header_type)) then
         ret = xmlTextReaderMoveToAttribute(reader, "type")
         if (ret==1) then
            header_type => xmlTextReaderConstValue(reader) 
         else
            nullify(header_type)
         end if
      end if
      if (present(compressed)) then
         ret = xmlTextReaderMoveToAttribute(reader, "compressor")
         compressed = .false.
         if (ret==1) then
            tmp => xmlTextReaderConstValue(reader)
            compressed = match(tmp, "vtkZLibDataCompressor")
         end if
      end if

    end subroutine process_vtk_header

    function count_children(reader)
      type(xmlTextReader), intent(in) :: reader
      integer :: count_children
      count_children = xmlChildElementCount(xmlTextReaderExpand(reader))
    end function count_children

    function get_piece_filename(filename, piece_no, encoding, options) result(piecename)
      character, dimension(:), pointer :: piecename
      character (len=*), intent(in) :: filename
      integer, intent(in) :: piece_no
      character (len=*), intent(in), optional :: encoding
      character, dimension(:), pointer :: ele_name
      integer, intent(in), optional :: options

      type(xmlTextReader) :: reader

      integer :: ret, count

      reader = open_vtk_xml_file(filename, encoding, options)
      count = 0
      nullify(piecename)

      do while (read_vtk_xml_element(reader, ele_name) )
         if (match(ele_name, "Piece")) then
            if (count == piece_no) then
               ret = xmlTextReaderMoveToAttribute(reader, "Source")
               piecename => xmlTextReaderConstValue(reader)
               exit
            end if
            count = count+1
         end if
      end do

      call close_vtk_xml_file(reader)

    end function get_piece_filename

    subroutine process_piece(reader, number_of_points, number_of_cells)
      type(xmlTextReader), intent(in) :: reader
      integer, intent(out) :: number_of_points, number_of_cells
      integer :: ret, i ,j
      character, dimension(:), pointer :: tmp

      ret = xmlTextReaderMoveToAttribute(reader, "NumberOfPoints")
      if (ret==1) then
         tmp => xmlTextReaderConstValue(reader)
         number_of_points = chararray2int(tmp)
      else
         number_of_points = -1
      end if
      ret = xmlTextReaderMoveToAttribute(reader, "NumberOfCells")
      if (ret==1) then
         tmp => xmlTextReaderConstValue(reader)
         number_of_cells = chararray2int(tmp)
      else
         number_of_cells = -1
      end if

    end subroutine process_piece
 
    subroutine read_points(reader, positions, compressed)
      type(xmlTextReader), intent(in) :: reader
      real, dimension(:,:) :: positions
      logical :: compressed

      character(kind=c_char), dimension(:), pointer :: datatype
      character(kind=c_char), dimension(:), pointer :: cdata
      integer :: i, j, ret
      logical :: lret

      lret = read_vtk_xml_element(reader)

      ret = xmlTextReaderMoveToAttribute(reader, "type")
      datatype => xmlTextReaderConstValue(reader)

      ret = xmlTextReaderRead(reader)
      cdata => xmlTextReaderConstValue(reader)

      if (compressed) then
         if (match(datatype,"Float32")) then
            positions = reshape(transfer(read_compressed_data(cdata, 4*size(positions)), &
                 real(1.0, kind=c_float), size(positions)), shape(positions))
         else
            positions = reshape(transfer(read_compressed_data(cdata, 8*size(positions)), &
                 real(1.0, kind=c_double), size(positions)), shape(positions))
         end if
      else
         read(cdata,*) positions
      end if

    end subroutine read_points

    subroutine read_cells(reader, ncells, connectivity, offsets, celltypes, compressed)
      type(xmlTextReader), intent(in) :: reader
      integer :: ncells
      integer, dimension(:), pointer :: connectivity, offsets, celltypes
      logical :: compressed, lret
      character(kind=c_char), dimension(:), pointer :: name, data, datatype, dataformat
      character(kind=c_char), dimension(:), allocatable, target :: cdata, odata, tdata
      integer :: cdtype, odtype, tdtype
      logical :: cbinary, obinary, tbinary
      
      integer :: i, ret

      do i=1,3
         lret = read_vtk_xml_element(reader)
         ret = xmlTextReaderMoveToAttribute(reader, "Name")
         name => xmlTextReaderConstValue(reader)
         ret = xmlTextReaderMoveToAttribute(reader, "type")
         datatype =>  xmlTextReaderConstValue(reader)
         ret = xmlTextReaderMoveToAttribute(reader, "format")
         dataformat =>  xmlTextReaderConstValue(reader)
         ret = xmlTextReaderRead(reader)
         select case(str(name))
         case("connectivity")
            call copy(cdata,xmlTextReaderConstValue(reader))
            cdtype = kindpar(str(datatype))
            cbinary = str(dataformat)=="binary"
         case("offsets")
            call copy(odata,xmlTextReaderConstValue(reader))
            odtype = kindpar(str(datatype))
            obinary = str(dataformat)=="binary"
         case("types")
            call copy(tdata, xmlTextReaderConstValue(reader))
            tdtype = kindpar(str(datatype))
            tbinary = str(dataformat)=="binary"
         end select
      end do

      allocate(celltypes(ncells), offsets(ncells))
      data=>tdata
      if (tbinary) then
         if (compressed) then
            call transfer_integer_data(celltypes, data, tdtype)
         end if
      else
         read(data,*) celltypes
      end if
      data=>odata
      if (obinary) then
         if (compressed) then
            call transfer_integer_data(offsets, data, odtype)
         end if
      else
         read(data,*) offsets
      end if
      allocate(connectivity(offsets(ncells)))
      data=>cdata
      if (cbinary) then
         if (compressed) then
            call transfer_integer_data(connectivity, data, cdtype)
         end if
      else
         read(data,*) connectivity
      end if
    end subroutine read_cells

    subroutine read_data_array_ncomps(reader, ncomps, name)
      type(xmlTextReader) :: reader
      character(kind=c_char), dimension(:), pointer, optional :: name
      integer :: ncomps, ret, i, j
      logical :: lret
      character(kind=c_char), dimension(:), pointer :: tmp

      lret = read_vtk_xml_element(reader)

      if (present(name)) then
         ret = xmlTextReaderMoveToAttribute(reader, "Name")
         name => xmlTextReaderConstValue(reader)
      end if
      
      ret = xmlTextReaderMoveToAttribute(reader, "NumberOfComponents")
      if (ret<1) then
         ncomps = 1
      else
         tmp => xmlTextReaderConstValue(reader)
         if (match(tmp,"")) then
            ncomps = 1
         else
            ncomps = chararray2int(tmp)
         end if
      end if

    end subroutine read_data_array_ncomps

    subroutine read_i8data_array(reader, data, &
         compressed, offset, dlen)
      type(xmlTextReader) :: reader
      integer(kind=c_int8_t), pointer, dimension(:,:) :: data
      integer, optional :: dlen
      integer :: offset
      logical :: lret, compressed
      integer :: ncomps, ntuples, ret, i, j
      real(c_float), dimension(:), allocatable :: fdata
      real(c_double), dimension(:), allocatable :: ddata
      character(kind=c_char), dimension(:), pointer :: tmp, cdata, dformat

      ret = xmlTextReaderMoveToAttribute(reader, "NumberOfComponents")
      if (ret<1) then
         ncomps = 1
      else
         tmp => xmlTextReaderConstValue(reader)
         if (match(tmp,"")) then
            ncomps = 1
         else
            ncomps = chararray2int(tmp)
         end if
      end if

      ret = xmlTextReaderMoveToAttribute(reader, "format")
      dformat => xmlTextReaderConstValue(reader)

      if (match(dformat, "appended")) then
         ret = xmlTextReaderMoveToAttribute(reader, "offset")
         tmp => xmlTextReaderConstValue(reader)
         offset = chararray2int(tmp)
         return
      end if
         
      if (present(dlen)) then
         ntuples = dlen
      else
         ret = xmlTextReaderMoveToAttribute(reader, "NumberOfTuples")
         if (ret<1) then
            ntuples = 1
         else
            tmp => xmlTextReaderConstValue(reader)
            ntuples = chararray2int(tmp)
         end if
      end if

      ret = xmlTextReaderMoveToAttribute(reader, "type")
      tmp => xmlTextReaderConstValue(reader)

      ret = xmlTextReaderRead(reader)
      cdata => xmlTextReaderConstValue(reader)

      if (match(dformat, "binary")) then
         if (compressed) then
            if (match(tmp,"UInt8")) then
               data = reshape(transfer(read_compressed_data(cdata, size(data)), &
                    int(1,c_int8_t), size(data)), shape(data))
            else
               data = reshape(transfer(read_compressed_data(cdata, 4*size(data)), &
                    int(1, kind=c_int), size(data)), shape(data))
            end if
         else
            if (match(tmp,"UInt8")) then
               data = reshape(transfer(read_uncompressed_base64_data(cdata, size(data)), &
                    int(1, kind=c_int8_t), size(data)), shape(data))
            else
               data = reshape(transfer(read_uncompressed_base64_data(cdata, 4*size(data)), &
                    int(1, kind=c_int), size(data)), shape(data))
            end if
         end if
      else
         read(cdata,*) data
      end if
      
    end subroutine read_i8data_array


    subroutine read_appended_i8data_array(reader, name, data, &
         compressed, offset, dlen, dformat)
      type(xmlTextReader) :: reader
      character(kind=c_char), dimension(:), pointer :: name
      integer(kind=c_int8_t), pointer, dimension(:,:) :: data
      integer, optional :: dlen
      character(len=*) :: dformat
      integer :: offset
      logical :: lret, compressed
      integer, target :: nbytes
      integer :: ncomps, ntuples, ret, i, j, nd
      real(c_float), dimension(:), allocatable :: fdata
      real(c_double), dimension(:), allocatable :: ddata
      character(kind=c_char), dimension(:), pointer :: tmp, cdata
      integer(c_long) :: destLen, sourcelen

      ret = xmlTextReaderRead(reader)
      cdata => xmlTextReaderConstValue(reader)
      j = 2
      do while (cdata(j-1) /= '_')
         j = j+1
      end do
      cdata => cdata(offset+j:size(cdata))

      if (compressed) then
         if (dformat == "UInt8") then
            data = reshape(transfer(read_compressed_data(cdata, &
                 size(data)), &
                 int(1,c_int8_t), size(data)), shape(data))
         else
            data = reshape(transfer(read_compressed_data(cdata, &
                 4*size(data)), &
                 int(1, kind=c_int), size(data)), shape(data))
         end if
      else
         if (dformat == "UInt8") then
            data = reshape(transfer(read_uncompressed_base64_data(cdata, size(data)), &
                    int(1, kind=c_int8_t), size(data)), shape(data))
         else
            data = reshape(transfer(read_uncompressed_base64_data(cdata, 4*size(data)), &
                    int(1, kind=c_int), size(data)), shape(data))
         end if
      end if
      
    end subroutine read_appended_i8data_array

    subroutine read_idata_array(reader, data, &
         compressed, offset, dlen)
      type(xmlTextReader) :: reader
      integer(kind=c_int), pointer, dimension(:,:) :: data
      integer, optional :: dlen
      integer :: offset
      logical :: lret, compressed
      integer :: ncomps, ntuples, ret, i, j
      real(c_float), dimension(:), allocatable :: fdata
      real(c_double), dimension(:), allocatable :: ddata
      character(kind=c_char), dimension(:), pointer :: tmp, cdata, dformat

      ret = xmlTextReaderMoveToAttribute(reader, "NumberOfComponents")
      if (ret<1) then
         ncomps = 1
      else
         tmp => xmlTextReaderConstValue(reader)
         if (match(tmp,"")) then
            ncomps = 1
         else
            ncomps = chararray2int(tmp)
         end if
      end if

      ret = xmlTextReaderMoveToAttribute(reader, "format")
      dformat => xmlTextReaderConstValue(reader)

      if (match(dformat, "appended")) then
         ret = xmlTextReaderMoveToAttribute(reader, "offset")
         tmp => xmlTextReaderConstValue(reader)
         offset = chararray2int(tmp)
         return
      end if
         
      if (present(dlen)) then
         ntuples = dlen
      else
         ret = xmlTextReaderMoveToAttribute(reader, "NumberOfTuples")
         if (ret<1) then
            ntuples = 1
         else
            tmp => xmlTextReaderConstValue(reader)
            ntuples = chararray2int(tmp)
         end if
      end if

      ret = xmlTextReaderMoveToAttribute(reader, "type")
      tmp => xmlTextReaderConstValue(reader)

      ret = xmlTextReaderRead(reader)
      cdata => xmlTextReaderConstValue(reader)

      if (match(dformat, "binary")) then
         if (compressed) then
            if (match(tmp,"UInt8")) then
               data = reshape(transfer(read_compressed_data(cdata, size(data)), &
                    int(1, kind=c_int8_t), size(data)), shape(data))
            else
               data = reshape(transfer(read_compressed_data(cdata, 4*size(data)), &
                    int(1, kind=c_int), size(data)), shape(data))
            end if
         else
            if (match(tmp,"UInt8")) then
               data = reshape(transfer(read_uncompressed_base64_data(cdata, size(data)), &
                    int(1, kind=c_int8_t), size(data)), shape(data))
            else
               data = reshape(transfer(read_uncompressed_base64_data(cdata, 4*size(data)), &
                    int(1, kind=c_int), size(data)), shape(data))
            end if
         end if
      else
         read(cdata,*) data
      end if
      
    end subroutine read_idata_array


    subroutine read_appended_idata_array(reader, name, data, cdata, &
         compressed, offset, dlen, dformat)
      type(xmlTextReader) :: reader
      character(kind=c_char), dimension(:), pointer :: name, cdata
      integer(kind=c_int), pointer, dimension(:,:) :: data
      integer, optional :: dlen
      character(len=*) :: dformat
      integer :: offset
      logical :: lret, compressed
      integer, target :: nbytes
      integer :: ncomps, ntuples, ret, i, j, nd
      real(c_float), dimension(:), allocatable :: fdata
      real(c_double), dimension(:), allocatable :: ddata
      character(kind=c_char), dimension(:), pointer :: tmp
      integer(c_long) :: destLen, sourcelen

      cdata => cdata(offset+1:size(cdata))

      if (compressed) then
         if (dformat == "UInt8") then
            data = reshape(transfer(read_compressed_data(cdata, &
                 size(data)), &
                 int(1,c_int8_t), size(data)), shape(data))
         else
            data = reshape(transfer(read_compressed_data(cdata, &
                 4*size(data)), &
                 int(1, kind=c_int), size(data)), shape(data))
         end if
      else
         if (dformat == "UInt8") then
            data = reshape(transfer(read_uncompressed_base64_data(cdata, size(data)), &
                    int(1, kind=c_int8_t), size(data)), shape(data))
         else
            data = reshape(transfer(read_uncompressed_base64_data(cdata, 4*size(data)), &
                    int(1, kind=c_int), size(data)), shape(data))
         end if
      end if
      
    end subroutine read_appended_idata_array

    subroutine read_data_array(reader, name, data, ncomps, &
         dlen, compressed, appended, offset)
      type(xmlTextReader) :: reader
      character(kind=c_char), dimension(:), pointer :: name
      real, pointer, dimension(:,:) :: data
      integer, optional :: dlen, offset
      logical :: lret, compressed, appended
      integer :: ncomps, ntuples, ret, i, j
      real(c_float), dimension(:), allocatable :: fdata
      real(c_double), dimension(:), allocatable :: ddata
      character(kind=c_char), dimension(:), pointer :: tmp, cdata, dformat

      lret = read_vtk_xml_element(reader)
      ret = xmlTextReaderMoveToAttribute(reader, "Name")
      name => xmlTextReaderConstValue(reader)

      ret = xmlTextReaderMoveToAttribute(reader, "NumberOfComponents")
      if (ret<1) then
         ncomps = 1
      else
         tmp => xmlTextReaderConstValue(reader)
         if (match(tmp,"")) then
            ncomps = 1
         else
            ncomps = chararray2int(tmp)
         end if
      end if

      ret = xmlTextReaderMoveToAttribute(reader, "format")
      dformat => xmlTextReaderConstValue(reader)

      if (match(dformat, "appended")) then
         ret = xmlTextReaderMoveToAttribute(reader, "offset")
         tmp => xmlTextReaderConstValue(reader)
         offset = chararray2int(tmp)
         appended=.true.
         return
      else
         offset = -1
         appended=.false.
      end if
         
      if (present(dlen)) then
         ntuples = dlen
      else
         ret = xmlTextReaderMoveToAttribute(reader, "NumberOfTuples")
         if (ret<1) then
            ntuples = 1
         else
            tmp => xmlTextReaderConstValue(reader)
            ntuples = chararray2int(tmp)
         end if
      end if

      ret = xmlTextReaderMoveToAttribute(reader, "type")
      tmp => xmlTextReaderConstValue(reader)

      ret = xmlTextReaderRead(reader)
      cdata => xmlTextReaderConstValue(reader)

      allocate(data(ncomps, ntuples))

      if (match(dformat, "binary")) then
         if (compressed) then
            if (match(tmp,"Float32")) then
               data = reshape(transfer(read_compressed_data(cdata, 4*size(data)), &
                    real(1.0, kind=c_float), size(data)), shape(data))
            elseif (match(tmp,"UInt8")) then
               data = reshape(transfer(read_compressed_data(cdata, size(data)), &
                    int(1, kind=c_int8_t), size(data)), shape(data))
            else
               data = reshape(transfer(read_compressed_data(cdata, 8*size(data)), &
                    real(1.0, kind=c_double), size(data)), shape(data))
            end if
         else
            if (match(tmp,"Float32")) then
               data = reshape(transfer(read_uncompressed_base64_data(cdata, 4*size(data)), &
                    real(1.0, kind=c_float), size(data)), shape(data))
            else
               data = reshape(transfer(read_uncompressed_base64_data(cdata, 8*size(data)), &
                    real(1.0, kind=c_double), size(data)), shape(data))
            end if
         end if
      else
         read(cdata,*) data
      end if
         
         
    end subroutine read_data_array

    subroutine read_appended_data_array(reader, name, data, cdata, &
         compressed, offset, dlen, ntuples, dformat)
      type(xmlTextReader) :: reader
      character(kind=c_char), dimension(:), pointer :: name, cdata
      real(kind=c_double), pointer, dimension(:,:) :: data
      integer, optional :: dlen
      character(len=*) :: dformat
      integer :: offset
      logical :: lret, compressed
      integer, target :: nbytes
      integer :: ncomps, ntuples, ret, i, j, nd
      real(c_float), dimension(:), allocatable :: fdata
      real(c_double), dimension(:), allocatable :: ddata
      character(kind=c_char), dimension(:), pointer :: tmp
      integer(c_long) :: destLen, sourcelen

           
      cdata => cdata(offset+1:size(cdata))
      allocate(data(ntuples, dlen))

      if (compressed) then
         if (dformat == "Float32") then
            data = reshape(transfer(read_compressed_data(cdata, &
                 4*size(data)), &
                 real(1.0,c_float), size(data)), shape(data))
         else
            data = reshape(transfer(read_compressed_data(cdata, &
                 8*size(data)), &
                 real(1.0, kind=c_double), size(data)), shape(data))
         end if
      else
         if (dformat == "Float32") then
            data = reshape(transfer(read_uncompressed_base64_data(cdata, 4*size(data)), &
                    real(1.0, kind=c_double), size(data)), shape(data))
         else
            data = reshape(transfer(read_uncompressed_base64_data(cdata, 8*size(data)), &
                    real(1.0, kind=c_int), size(data)), shape(data))
         end if
      end if
      
    end subroutine read_appended_data_array

    function read_uncompressed_base64_data(cdata, idestLen) result(data)
      character(kind=c_char), dimension(:), pointer :: cdata
      integer :: idestLen
      character(kind=c_char), dimension(idestLen) :: data

      character(kind=c_char), dimension(4096) :: comptext
      integer :: cLength(1), ret
      type(b64state) :: state_in

      call base64_init_decodestate(state_in)
      ret = base64_decode_block(cdata(1:4), comptext, state_in)
      clength = transfer(comptext, cLength, 1)

      call base64_init_decodestate(state_in)
      ret = base64_decode_block(cdata(5:size(cdata)), data, state_in)

    end function read_uncompressed_base64_data

    function read_compressed_data(cdata, idestLen) result(data)
      character(kind=c_char), dimension(:), pointer :: cdata
      integer :: idestLen
      character(kind=c_char), dimension(idestLen), target :: data
      character(kind=c_char), dimension(:), pointer :: datatype
      integer(kind=c_int) :: metadata(3)
      integer(kind=c_int), allocatable :: cblocksizes(:)
       character(kind=c_char), dimension(:), pointer :: tmp

      integer :: i, offset, doffset, ret
      integer(c_long) :: destLen, sourceLen

      call read_compressed_header(cdata, metadata, cblocksizes)


      allocate(tmp(metadata(1)* metadata(2)))
      doffset=0
      do i=1, metadata(1)-1
         offset = sum(cblocksizes(1:i-1)) 
         destLen = metadata(2)
         sourceLen = cblocksizes(i)
         call read_zlib_compressed_block(cdata, offset, destlen, &
              c_loc(tmp(doffset+1)), sourcelen, 4)
         doffset = doffset + metadata(2)
      end do
      offset = sum(cblocksizes(1:metadata(1)-1)) 
      destLen = metadata(3)
      sourceLen = cblocksizes(metadata(1))
      call read_zlib_compressed_block(cdata, offset, destlen, &
           c_loc(tmp(doffset+1)), sourcelen, 4)
      data = tmp(1:size(data))
      deallocate(tmp)

    end function read_compressed_data

    subroutine close_vtk_xml_file(reader)
      type(xmlTextReader) :: reader
      call xmlFreeTextReader(reader)
    end subroutine close_vtk_xml_file

  subroutine vtk_get_sizes( filename, namelen, nnod, nelm, szenls, &
       nfield_components, nprop_components, &
       nfields, nproperties, ndimensions, maxnamelen )
    implicit none
    character*(*) filename
    integer namelen
    integer nnod, nelm, szenls
    integer nfield_components, nprop_components
    integer nfields, nproperties, ndimensions, maxnamelen

    type(xmlTextReader) :: reader
    character, dimension(:), pointer :: name
    logical :: lret, compressed, appended
    integer :: i, ret, ncomps, celloffsets
    integer(kind=c_int8_t), pointer :: types(:,:)

    nfields=0
    nproperties=0
    nfield_components = 0
    nprop_components = 0

    reader = open_vtk_xml_file(trim(filename))
    do while( read_vtk_xml_element(reader, name))
       select case(str(name))
       case("VTKFile")
          call process_vtk_header(reader, compressed=compressed)
       case("Piece")
          call process_piece(reader, nnod, nelm)
          allocate(types(1, nelm))
       case("PointData")
          nfields = count_children(reader)
          do i=1, nfields
             call read_data_array_ncomps(reader, ncomps)
             nfield_components = nfield_components+ncomps
          end do
       case("CellData")
          nproperties = count_children(reader)
          do i=1, nproperties
             call read_data_array_ncomps(reader, ncomps)
             nprop_components = nprop_components+ncomps
          end do
       case("Cells")
          do i=1, count_children(reader)
             lret = read_vtk_xml_element(reader)
             ret = xmlTextReaderMoveToAttribute(reader, "Name")
             name => xmlTextReaderConstValue(reader)
             select case(str(name))
             case("types")
                call read_i8data_array(reader, types, &
                     compressed, celloffsets, nelm)
             end select
          end do
       case("AppendedData")
          call read_appended_i8data_array(reader, name, types,&
               compressed, celloffsets, nelm, "UInt8")
       case("Points", "FieldData")
          continue
       end select
    end do

    select case(types(1,1))
    case(3)
       ndimensions=1
       szenls = 2*nelm
    case(21)
       ndimensions=1
       szenls = 3*nelm
    case(5)
       ndimensions=2
       szenls = 3*nelm
    case(22)
       ndimensions=2
       szenls = 6*nelm
    case(9)
       ndimensions=2
       szenls = 4*nelm
    case(12)
       ndimensions=3
       szenls = 8*nelm
    case(10)
       ndimensions=3
       szenls = 4*nelm
    case(24)
       ndimensions=3
       szenls = 10*nelm
    end select

    deallocate(types)

    call close_vtk_xml_file(reader)
    
  end subroutine vtk_get_sizes

  subroutine vtk_read_file(&
          filename, namelen, nnod, nelm, szenls,&
          nfield_components, nprop_components,&
          nfields, nproperties, &
          ndimensions, maxnamelen, &
          x, y, z, &
          field_components, prop_components, &
          fields, properties,&
          enlbas, enlist, field_names, prop_names)
       implicit none
       character*(*) filename
       integer namelen, nnod, nelm, szenls
       integer nfield_components, nprop_components, nfields, nproperties
       integer ndimensions, maxnamelen
       real x(nnod), y(nnod), z(nnod)
       integer field_components(nfields), &
            prop_components(nproperties)
       real, target :: fields(nnod,nfield_components), &
            properties(nelm,nprop_components) 
       integer enlbas(nelm+1), enlist(szenls)
       character(len=maxnamelen) field_names(nfield_components)
       character(len=maxnamelen) prop_names(nprop_components)

       type(xmlTextReader) :: reader
       character(kind=c_char), dimension(:), pointer :: name, cdata, ldata
       real, dimension(:,:), pointer :: rdata=>null(), pnts
       logical :: lret, compressed, appended
       integer :: i, j, ret, counter
       integer :: field_offsets(nfields), &
            prop_offsets(nproperties), pnts_offset,&
            connectivity_offset
       integer, dimension(:,:), pointer :: itmp

    reader = open_vtk_xml_file(trim(filename))
    do while(read_vtk_xml_element(reader, name))
       select case(str(name))
       case("VTKFile")
          call process_vtk_header(reader, compressed=compressed)
       case("Piece")
          continue
       case("PointData")
          counter= 0 
          do i=1, nfields
             call read_data_array(reader, name, rdata, &
                  field_components(i), nnod, compressed, appended, &
                  field_offsets(i))
             if (associated(rdata)) then
                do j=1, field_components(i)
                   fields(:,counter+j) = rdata(j,:)
                end do
                deallocate(rdata)
             end if
             field_names(i) = str(name)
             counter = counter + field_components(i)
          end do
       case("CellData")
          counter = 0
          do i=1, nproperties
             call read_data_array(reader, name, rdata, &
                  prop_components(i), nelm, compressed, appended, &
                  prop_offsets(i))
             if (associated(rdata)) then
                do j=1, prop_components(i)
                   properties(:,counter+j) = rdata(j,:)
                end do
                deallocate(rdata)
             end if
             prop_names(i) = str(name)
             counter = counter + prop_components(i)
          end do
       case("Points")
          allocate(pnts(3, nnod))
          call read_data_array(reader, name, pnts, &
               i, nnod, compressed, appended, &
               pnts_offset)
       case("Cells")
          do i=1, count_children(reader)
             lret = read_vtk_xml_element(reader)
             ret = xmlTextReaderMoveToAttribute(reader, "Name")
             name => xmlTextReaderConstValue(reader)
             select case(str(name))
             case("connectivity")
                allocate(itmp(1, szenls))
                call read_idata_array(reader, itmp, &
                     compressed, connectivity_offset, szenls)
             end select
          end do
       case("AppendedData")
          ret = xmlTextReaderRead(reader)
          cdata => xmlTextReaderConstValue(reader)
          I = 2
          do while (cdata(i-1) /= '_')
             i = i+1
          end do
          cdata => cdata(i:size(cdata))
          counter= 0
          do i=1, nfields
             ldata=> cdata
             call read_appended_data_array(reader, name, rdata, ldata,&
                  compressed, field_offsets(i), nnod, &
                  field_components(i),  "Float64")
             do j=1, field_components(i)
                fields(:,counter+j) = rdata(j,:)
             end do
             deallocate(rdata)
             counter = counter + field_components(i)
          end do

          counter= 0
          do i=1, nproperties
             ldata=> cdata
             rdata => properties(counter:counter+prop_components(i)-1,:)
             call read_appended_data_array(reader, name, rdata, ldata, &
                  compressed, prop_offsets(i), nelm, &
                  prop_components(i), "Float64")
             do j=1, field_components(i)
                properties(:,counter+j) = rdata(j,:)
             end do
             deallocate(rdata)
             counter = counter + prop_components(i)
          end do

          ldata=> cdata
          call read_appended_data_array(reader, name, pnts, ldata, &
               compressed, pnts_offset, nnod,  3, "Float64")

          ldata=> cdata
          call read_appended_idata_array(reader, name, itmp, ldata, &
               compressed, connectivity_offset, szenls, "Int32")
       case("FieldData")
          continue
       end select

    end do

    x = pnts(1,:)
    if (ndimensions>1) y = pnts(2,:)
    if (ndimensions>2) z = pnts(3,:)
    enlist = itmp(1, :) + 1
    
    deallocate(pnts, itmp)
       
  end subroutine vtk_read_file
  
end module vtkreader
