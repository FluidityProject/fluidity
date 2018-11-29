module vtkwriter

#ifndef HAVE_VTK

  use iso_c_binding
  use global_parameters, only: FIELD_NAME_LEN
  use parallel_tools
  
  use xmlwriter
  
  implicit none

  private

  type(xmlTextWriter) :: writer, pwriter

  real, dimension(:,:), pointer :: points
  integer, dimension(:), pointer :: connectivity, offsets
  integer(c_int8_t), dimension(:), pointer :: types

  integer :: mode=0
  logical, parameter :: compress =.true.
  character(len=FIELD_NAME_LEN) :: active_scalars='',&
       active_vectors='',&
       active_tensors='',&
       lfilename=''
  character(len=4) :: lextension
  
  integer, public, parameter :: VTK_VERTEX=1, VTK_LINE=3, VTK_TRIANGLE=5, VTK_QUAD=9, &
       VTK_TETRA=10, VTK_HEXAHEDRON=12, VTK_QUADRATIC_EDGE=21, &
       VTK_QUADRATIC_TRIANGLE=22, VTK_QUADRATIC_QUAD=23, &
       VTK_QUADRATIC_TETRA=24,  VTK_QUADRATIC_HEXAHEDRON=25

  public :: vtkopen, vtkclose, vtkwritemesh, &
       vtkwritesc, vtkwritesn, vtkwritevc, vtkwritevn, &
       vtkwritetn, vtkwritetc, vtksetactivescalars, &
       vtksetactivevectors,  vtksetactivetensors, &
       vtkpclose, vtkwriteghostlevels

  interface vtkwritesc
     module procedure vtkwritesc_int, &
          vtkwritesc_int8, vtkwritesc_real
  end interface vtkwritesc

  interface vtkwritesn
     module procedure vtkwritesn_real, vtkwritesn_int
  end interface vtkwritesn
  
contains

  logical function ishead()

    ishead = isparallel() .and. (getrank()==0) .and. lextension=="pvtu"

  end function ishead

    subroutine vtkopen(filename, title)
      character(len=*) :: filename, title
      integer :: err

      logical, parameter :: bigend = ichar(transfer(1,'a')) == 0

      mode = 0

      if (isparallel()) then
         lfilename = filename(1:index(filename, '.')-1)
         lextension = filename(index(filename, '.')+1:)
         select case(trim(lextension))
         case("pvtu")
            call system("mkdir -p "//lfilename)
            writer=xmlNewTextWriterFilename(trim(lfilename)//'/'&
                 //trim(lfilename)//'_'//int2str(getrank())//'.vtu', 0)
            if (getrank()==0) pwriter=xmlNewTextWriterFilename(trim(filename),0)
         case("vtu")
            writer = xmlNewTextWriterFilename(filename, 0)
         end select
      else
         writer = xmlNewTextWriterFilename(filename, 0)
      end if
      err = xmlTextWriterSetIndent(writer, 2)
      err = xmlTextWriterStartDocument(writer, "1.0", "utf-8", "no")
      err = xmlTextWriterStartElement(writer, "VTKFile")
      err = xmlTextWriterWriteAttribute(writer, "type", "UnstructuredGrid")
      err = xmlTextWriterWriteAttribute(writer, "version", "0.1")
      if (bigend) then
         err = xmlTextWriterWriteAttribute(writer, "byte_order", "BigEndian")
      else
         err = xmlTextWriterWriteAttribute(writer, "byte_order", "LittleEndian")
      end if
      err = xmlTextWriterWriteAttribute(writer, "header_type", "UInt32")
      if (compress) err = xmlTextWriterWriteAttribute(writer, "compressor", "vtkZLibDataCompressor")
      err = xmlTextWriterStartElement(writer, "UnstructuredGrid")

      if (ishead()) then
         err = xmlTextWriterSetIndent(pwriter, 2)
         err = xmlTextWriterStartDocument(pwriter, "1.0", "utf-8", "no")
         err = xmlTextWriterStartElement(pwriter, "VTKFile")
         err = xmlTextWriterWriteAttribute(pwriter, "type", "PUnstructuredGrid")
         err = xmlTextWriterWriteAttribute(pwriter, "version", "0.1")
         if (bigend) then
            err = xmlTextWriterWriteAttribute(pwriter, "byte_order", "BigEndian")
         else
            err = xmlTextWriterWriteAttribute(pwriter, "byte_order", "LittleEndian")
         end if
         err = xmlTextWriterWriteAttribute(pwriter, "header_type", "UInt32")
         if (compress) err = xmlTextWriterWriteAttribute(pwriter, "compressor", "vtkZLibDataCompressor")
         err = xmlTextWriterStartElement(pwriter, "PUnstructuredGrid")
         err = xmlTextWriterWriteAttribute(pwriter, "GhostLevel", "1")
      end if
      
    end subroutine vtkopen

    subroutine vtkclose()
      integer :: err

      if (mode >0) err = xmlTextWriterEndElement(writer)
      call add_cell_descriptions()

      err = xmlTextWriterEndElement(writer)
      err = xmlTextWriterEndElement(writer)
      err = xmlTextWriterEndDocument(writer)
      call xmlFreeTextWriter(writer)
      
      deallocate(points, offsets, connectivity, types)
      active_scalars=''
      active_vectors=''
      active_tensors=''
      mode = 0
      
    end subroutine vtkclose

    subroutine vtkwritemesh(npoints, ncells, X, Y, Z, ENlist, ELType, ElSize)
      integer, intent(in) :: npoints, ncells
      real, dimension(:), intent(in) :: X, Y, Z
      integer, dimension(:), intent(in) :: ENList, Elsize, ELType

      integer :: err, i

      allocate(points(3,npoints), connectivity(size(ENlist)), &
           offsets(ncells), types(ncells))

      points(1,:) = X
      points(2,:) = Y
      points(3,:) = Z

      connectivity = ENList-1
      
      if (ncells>0) then
         offsets(1) = elsize(1)
         connectivity(1:offsets(1)) = elperm(ENList(1:elsize(1)), ELType(1))
         do i=2, ncells
            offsets(i) = offsets(i-1) + elsize(i)
            connectivity(offsets(i-1)+1:offsets(i)) = &
                 elperm(ENList(offsets(i-1)+1:offsets(i)), ELType(i))
         end do
      end if

      types = ElType

      err = xmlTextWriterStartElement(writer, "Piece")
      err = xmlTextWriterWriteAttribute(writer, &
           "NumberOfPoints", int2str(npoints))
      err = xmlTextWriterWriteAttribute(writer, &
           "NumberOfCells", int2str(ncells))

    contains

      function elperm(fl_order, vtktype)
        integer, dimension(:), intent(in) :: fl_order
        integer, intent(in) :: vtktype
        integer, dimension(size(fl_order)) :: elperm

        select case(vtktype)
        case(VTK_QUAD)
           elperm(1) = fl_order(1)-1
           elperm(2) = fl_order(2)-1
           elperm(3) = fl_order(4)-1
           elperm(4) = fl_order(3)-1
        case(VTK_HEXAHEDRON)
           elperm(1) = fl_order(1)-1
           elperm(2) = fl_order(2)-1
           elperm(3) = fl_order(4)-1
           elperm(4) = fl_order(3)-1
           elperm(5) = fl_order(5)-1
           elperm(6) = fl_order(6)-1
           elperm(7) = fl_order(8)-1
           elperm(8) = fl_order(7)-1
        case default
           elperm = fl_order -1
        end select
      end function elperm
      
    end subroutine vtkwritemesh

    subroutine vtkwritesc_real(val, name)
      real, intent(in), target :: val(:)
      character(len=*) :: name

      integer :: err 

      if (mode==1) then
         err = xmlTextWriterEndElement(writer)
         if (ishead()) err = xmlTextWriterEndElement(pwriter)
      end if
      if (mode<2) then
         err = xmlTextWriterStartElement(writer, "CellData")
         if (ishead()) err = xmlTextWriterStartElement(pwriter, "PCellData")
         mode=2
      end if

      call add_binary_data(writer, name, "Float64", c_loc(val), &
           1, 8*size(val), compress)

      if (ishead()) then
         err = xmlTextWriterStartElement(pwriter, "PDataArray")
         err = xmlTextWriterWriteAttribute(pwriter, "type", "Float64")
         err = xmlTextWriterWriteAttribute(pwriter, "Name", name)
         err = xmlTextWriterWriteAttribute(pwriter, "NumberOfComponents", "")
         err = xmlTextWriterEndElement(pwriter)
      end if
      
    end subroutine vtkwritesc_real

    subroutine vtkwritesc_int(val, name)
      integer, intent(in), target :: val(:)
      character(len=*) :: name

      integer :: err

      if (mode==1) then
         err = xmlTextWriterEndElement(writer)
         if (ishead()) err = xmlTextWriterEndElement(pwriter)
      end if
      if (mode<2) then
         err = xmlTextWriterStartElement(writer, "CellData")
         if (ishead()) err = xmlTextWriterStartElement(pwriter, "PCellData")
         mode=2
      end if

      call add_binary_data(writer, name, "Int32", c_loc(val), &
           1, 4*size(val), compress)

      if (ishead()) then
         err = xmlTextWriterStartElement(pwriter, "PDataArray")
         err = xmlTextWriterWriteAttribute(pwriter, "type", "Int32")
         err = xmlTextWriterWriteAttribute(pwriter, "Name", name)
         err = xmlTextWriterWriteAttribute(pwriter, "NumberOfComponents", "")
         err = xmlTextWriterEndElement(pwriter)
      end if
      
    end subroutine vtkwritesc_int

    subroutine vtkwritesc_int8(val, name)
      integer(c_int8_t), intent(in), target :: val(:)
      character(len=*) :: name

      integer :: err

      if (mode==1) then
         err = xmlTextWriterEndElement(writer)
         if (ishead()) err = xmlTextWriterEndElement(pwriter)
      end if
      if (mode<2) then
         err = xmlTextWriterStartElement(writer, "CellData")
         if (ishead()) err = xmlTextWriterStartElement(pwriter, "PCellData")
         mode=2
      end if

      call add_binary_data(writer, name, "UInt8", c_loc(val), &
           1, size(val), compress)

      if (ishead()) then
         err = xmlTextWriterStartElement(pwriter, "PDataArray")
         err = xmlTextWriterWriteAttribute(pwriter, "type", "UInt8")
         err = xmlTextWriterWriteAttribute(pwriter, "Name", name)
         err = xmlTextWriterWriteAttribute(pwriter, "NumberOfComponents", "")
         err = xmlTextWriterEndElement(pwriter)
      end if
      
    end subroutine vtkwritesc_int8

    subroutine vtkwritesn_int(val, name)
      integer, intent(in), target :: val(:)
      character(len=*) :: name

      integer :: err

      if (mode<1) then
         err = xmlTextWriterStartElement(writer, "PointData")
         if (len_trim(active_scalars)>0) &
              err = xmlTextWriterWriteAttribute(writer, "Scalars", &
              trim(active_scalars))
         if (len_trim(active_vectors)>0) &
              err = xmlTextWriterWriteAttribute(writer, "Vectors", &
              trim(active_vectors))
         if (len_trim(active_tensors)>0) &
              err = xmlTextWriterWriteAttribute(writer, "Tensors", &
              trim(active_tensors))
         if (ishead()) then
            err = xmlTextWriterStartElement(Pwriter, "PPointData")
            if (len_trim(active_scalars)>0) &
                 err = xmlTextWriterWriteAttribute(Pwriter, "Scalars", &
                 trim(active_scalars))
            if (len_trim(active_vectors)>0) &
                 err = xmlTextWriterWriteAttribute(Pwriter, "Vectors", &
                 trim(active_vectors))
            if (len_trim(active_tensors)>0) &
                 err = xmlTextWriterWriteAttribute(Pwriter, "Tensors", &
                 trim(active_tensors))
         end if
         mode=1
      end if
      
      call add_binary_data(writer, name, "Int32", c_loc(val), &
           1, 4*size(val), compress)

      if (ishead()) then
         err = xmlTextWriterStartElement(pwriter, "PDataArray")
         err = xmlTextWriterWriteAttribute(pwriter, "type", "Int32")
         err = xmlTextWriterWriteAttribute(pwriter, "Name", name)
         err = xmlTextWriterWriteAttribute(pwriter, "NumberOfComponents", "")
         err = xmlTextWriterEndElement(pwriter)
      end if
      
    end subroutine vtkwritesn_int

    subroutine vtkwritesn_real(val, name)
      real, intent(in), target :: val(:)
      character(len=*) :: name

      integer :: err

      if (mode<1) then
         err = xmlTextWriterStartElement(writer, "PointData")
         if (len_trim(active_scalars)>0) &
              err = xmlTextWriterWriteAttribute(writer, "Scalars", &
              trim(active_scalars))
         if (len_trim(active_vectors)>0) &
              err = xmlTextWriterWriteAttribute(writer, "Vectors", &
              trim(active_vectors))
         if (len_trim(active_tensors)>0) &
              err = xmlTextWriterWriteAttribute(writer, "Tensors", &
              trim(active_tensors))
         if (ishead()) then
            err = xmlTextWriterStartElement(Pwriter, "PPointData")
            if (len_trim(active_scalars)>0) &
                 err = xmlTextWriterWriteAttribute(Pwriter, "Scalars", &
                 trim(active_scalars))
            if (len_trim(active_vectors)>0) &
                 err = xmlTextWriterWriteAttribute(Pwriter, "Vectors", &
                 trim(active_vectors))
            if (len_trim(active_tensors)>0) &
                 err = xmlTextWriterWriteAttribute(Pwriter, "Tensors", &
                 trim(active_tensors))
         end if
         mode=1
      end if
      
      call add_binary_data(writer, name, "Float64", c_loc(val), &
           1, 8*size(val), compress)

      if (ishead()) then
         err = xmlTextWriterStartElement(pwriter, "PDataArray")
         err = xmlTextWriterWriteAttribute(pwriter, "type", "Float64")
         err = xmlTextWriterWriteAttribute(pwriter, "Name", name)
         err = xmlTextWriterWriteAttribute(pwriter, "NumberOfComponents", "")
         err = xmlTextWriterEndElement(pwriter)
      end if
      
    end subroutine vtkwritesn_real

    subroutine vtkwritevn(val_x, val_y, val_z, name)
      real, intent(in) :: val_x(:), val_y(:), val_z(:)
      character(len=*) :: name

      real, dimension(3, size(val_x)), target :: val
      integer :: err

      val(1,:) = val_x
      val(2,:) = val_y
      val(3,:) = val_z

      if (mode<1) then
         err = xmlTextWriterStartElement(writer, "PointData")
         if (len_trim(active_scalars)>0) &
              err = xmlTextWriterWriteAttribute(writer, "Scalars", &
              trim(active_scalars))
         if (len_trim(active_vectors)>0) &
              err = xmlTextWriterWriteAttribute(writer, "Vectors", &
              trim(active_vectors))
         if (len_trim(active_tensors)>0) &
              err = xmlTextWriterWriteAttribute(writer, "Tensors", &
              trim(active_tensors))
         if (ishead()) then
            err = xmlTextWriterStartElement(Pwriter, "PPointData")
            if (len_trim(active_scalars)>0) &
                 err = xmlTextWriterWriteAttribute(Pwriter, "Scalars", &
                 trim(active_scalars))
            if (len_trim(active_vectors)>0) &
                 err = xmlTextWriterWriteAttribute(Pwriter, "Vectors", &
                 trim(active_vectors))
            if (len_trim(active_tensors)>0) &
                 err = xmlTextWriterWriteAttribute(Pwriter, "Tensors", &
                 trim(active_tensors))
         end if
         mode=1
      end if

      call add_binary_data(writer, name, "Float64", c_loc(val), &
           3, 8*size(val), compress)

      if (ishead()) then
         err = xmlTextWriterStartElement(pwriter, "PDataArray")
         err = xmlTextWriterWriteAttribute(pwriter, "type", "Float64")
         err = xmlTextWriterWriteAttribute(pwriter, "Name", name)
         err = xmlTextWriterWriteAttribute(pwriter, "NumberOfComponents", "3")
         err = xmlTextWriterEndElement(pwriter)
      end if
      
    end subroutine vtkwritevn

    subroutine vtkwritevc(val_x, val_y, val_z,  name)
      real, intent(in) :: val_x(:), val_y(:), val_z(:)
      character(len=*) :: name

      real, dimension(3, size(val_x)), target :: val
      integer :: err

      val(1,:) = val_x
      val(2,:) = val_y
      val(3,:) = val_z

      if (mode==1) then
         err = xmlTextWriterEndElement(writer)
         if (ishead()) err = xmlTextWriterEndElement(pwriter)
      end if
      if (mode<2) then
         err = xmlTextWriterStartElement(writer, "CellData")
         if (ishead()) err = xmlTextWriterStartElement(writer, "PCellData")
         mode=2
      end if

      call add_binary_data(writer, name, "Float64", c_loc(val), &
           3, 8*size(val), compress)

      if (ishead()) then
         err = xmlTextWriterStartElement(pwriter, "PDataArray")
         err = xmlTextWriterWriteAttribute(pwriter, "type", "Float64")
         err = xmlTextWriterWriteAttribute(pwriter, "Name", name)
         err = xmlTextWriterWriteAttribute(pwriter, "NumberOfComponents", "3")
         err = xmlTextWriterEndElement(pwriter)
      end if
      
    end subroutine vtkwritevc

    subroutine vtkwritetc(val_11, val_12, val_13,&
         val_21, val_22, val_23,&
         val_31, val_32, val_33, name)
      real, intent(in) :: val_11(:), val_12(:), val_13(:), &
           val_21(:), val_22(:), val_23(:), &
           val_31(:), val_32(:), val_33(:)
      character(len=*) :: name

      real, dimension(9, size(val_11)), target :: val
      integer :: err

      val(1,:) = val_11
      val(2,:) = val_12
      val(3,:) = val_13
      val(4,:) = val_21
      val(5,:) = val_22
      val(6,:) = val_23
      val(7,:) = val_31
      val(8,:) = val_32
      val(9,:) = val_33

      if (mode==1) then
         err = xmlTextWriterEndElement(writer)
         if (ishead()) err = xmlTextWriterEndElement(pwriter)
      end if
      if (mode<2) then
         err = xmlTextWriterStartElement(writer, "CellData")
         if (ishead()) err = xmlTextWriterStartElement(writer, "PCellData")
         mode=2
      end if

      call add_binary_data(writer, name, "Float64", c_loc(val), &
           9, 8*size(val), compress)

      if (ishead()) then
         err = xmlTextWriterStartElement(pwriter, "PDataArray")
         err = xmlTextWriterWriteAttribute(pwriter, "type", "Float64")
         err = xmlTextWriterWriteAttribute(pwriter, "Name", name)
         err = xmlTextWriterWriteAttribute(pwriter, "NumberOfComponents", "9")
         err = xmlTextWriterEndElement(pwriter)
      end if
      
    end subroutine vtkwritetc

    subroutine vtkwritetn(val_11, val_12, val_13,&
         val_21, val_22, val_23,&
         val_31, val_32, val_33, name)
      real, intent(in) :: val_11(:), val_12(:), val_13(:), &
           val_21(:), val_22(:), val_23(:), &
           val_31(:), val_32(:), val_33(:)
      character(len=*) :: name

      real, dimension(9, size(val_11)), target :: val
      integer :: err

      val(1,:) = val_11
      val(2,:) = val_12
      val(3,:) = val_13
      val(4,:) = val_21
      val(5,:) = val_22
      val(6,:) = val_23
      val(7,:) = val_31
      val(8,:) = val_32
      val(9,:) = val_33

      if (mode<1) then
         err = xmlTextWriterStartElement(writer, "PointData")
         if (len_trim(active_scalars)>0) &
              err = xmlTextWriterWriteAttribute(writer, "Scalars", &
              trim(active_scalars))
         if (len_trim(active_vectors)>0) &
              err = xmlTextWriterWriteAttribute(writer, "Vectors", &
              trim(active_vectors))
         if (len_trim(active_tensors)>0) &
              err = xmlTextWriterWriteAttribute(writer, "Tensors", &
              trim(active_tensors))
         if (ishead()) then
            err = xmlTextWriterStartElement(pwriter, "PPointData")
            if (len_trim(active_scalars)>0) &
                 err = xmlTextWriterWriteAttribute(pwriter, "Scalars", &
                 trim(active_scalars))
            if (len_trim(active_vectors)>0) &
                 err = xmlTextWriterWriteAttribute(pwriter, "Vectors", &
                 trim(active_vectors))
            if (len_trim(active_tensors)>0) &
                 err = xmlTextWriterWriteAttribute(pwriter, "Tensors", &
                 trim(active_tensors))
         end if
         mode=1
      end if

      call add_binary_data(writer, name, "Float64", c_loc(val), &
           9, 8*size(val), compress)

      if (ishead()) then
         err = xmlTextWriterStartElement(pwriter, "PDataArray")
         err = xmlTextWriterWriteAttribute(pwriter, "type", "Float64")
         err = xmlTextWriterWriteAttribute(pwriter, "Name", name)
         err = xmlTextWriterWriteAttribute(pwriter, "NumberOfComponents", "9")
         err = xmlTextWriterEndElement(pwriter)
      end if
      
    end subroutine vtkwritetn

    subroutine vtksetactivescalars(name)
      character(len=*) :: name

      active_scalars = name
    end subroutine vtksetactivescalars

    subroutine vtksetactivevectors(name)
      character(len=*) :: name

      active_vectors = name
    end subroutine vtksetactivevectors

    subroutine vtksetactivetensors(name)
      character(len=*) :: name

      active_tensors = name
    end subroutine vtksetactivetensors

    subroutine vtkpclose(rank, nparts)
      integer :: rank, nparts

      integer :: i, err

      if (mode >0) err = xmlTextWriterEndElement(pwriter)

      if (getrank()==0 .and. nparts>1) then
         do i=1, nparts
            err = xmlTextWriterStartElement(pwriter, "Piece")
            err = xmlTextWriterWriteAttribute(pwriter, "Source", &
                 trim(lfilename)//'/'//trim(lfilename) &
                 //'_'//int2str(i-1)//".vtu")
            err = xmlTextWriterEndElement(pwriter)
         end do
      end if

      call vtkclose()

      if (nparts>1) then
         err = xmlTextWriterEndElement(pwriter)
         err = xmlTextWriterEndDocument(pwriter)
         call xmlFreeTextWriter(pwriter)
      end if
      
    end subroutine vtkpclose

    subroutine add_binary_data(writer, name, datatype, data, ncomps, bytesize, compressed)

      type(xmlTextWriter) :: writer
      character(len=*) :: name, datatype
      type(c_ptr) :: data
      integer :: ncomps, bytesize
      logical :: compressed

      integer :: err

      err = xmlTextWriterStartElement(writer, "DataArray")
      err = xmlTextWriterWriteAttribute(writer, "type", datatype)
      err = xmlTextWriterWriteAttribute(writer, "Name", name)
      if (ncomps>1) then
         err = xmlTextWriterWriteAttribute(writer, "NumberOfComponents", int2str(ncomps))
      else
         err = xmlTextWriterWriteAttribute(writer, "NumberOfComponents", "")
      end if
      err = xmlTextWriterWriteAttribute(writer, "format", "binary")
      
      if (compressed) then
         call write_compressed_data(writer, data, bytesize) 
      else
         call write_base64(writer, asbytes([bytesize]) &
              //asbytes(data, bytesize), bytesize+4)
      end if
      err = xmlTextWriterEndElement(writer)
    end subroutine add_binary_data

    pure function zlib_size(bytesize, blocksize)
      integer, intent(in) :: bytesize, blocksize
      integer(c_long) :: sourceLen
      integer :: zlib_size

      interface
         pure function compressBound(bytesize) bind(c, name="compressBound")
           use iso_c_binding
           integer(c_long), value :: bytesize
           integer(c_long) :: compressBound
         end function compressBound
      end interface

      zlib_size = int((bytesize/blocksize+2)*compressBound(int(blocksize, kind=c_long)))


    end function zlib_size

    subroutine write_compressed_data(writer, data, n)
      type(xmlTextWriter) :: writer
      type(c_ptr) :: data
      integer :: n, k

      integer, parameter :: blocksize = 32768, level=-1
      character(kind=c_char), pointer :: dest(:)
      integer(c_int) ::  cerr
      integer(c_long) :: sourceLen, dlen, slen, destlen(n/blocksize+1)
      integer :: err
      integer, dimension(:), pointer :: metadata
      character, dimension(:), pointer :: cdata

      interface
         function compress2(dest, destLen, source, sourceLen, level) bind(c)
           use iso_c_binding
           integer(c_long) :: destLen
           integer(c_long), value :: sourceLen
           integer(c_int), value :: level 
           character(kind=c_char) :: dest(*)
           type(c_ptr), value :: source
           integer(c_int) :: compress2
         end function compress2
      end interface

      call c_f_pointer(data, cdata, [n])

      slen = n
      k = 0
      err = n/blocksize
      if (n-err*blocksize>0) err = err+1
      allocate(dest(zlib_size(n, blocksize)))
      destlen=0
      
      do while(slen>0)
         sourceLen = min(blocksize, slen)
         destlen(k+1) = size(dest)-sum(destlen(1:k))
         cerr = compress2(dest(sum(destlen(1:k))+1:size(dest)), destlen(k+1), &
              c_loc(cdata(1+blocksize*k)), sourceLen, level)
         slen = slen-blocksize
         k = k+1
      end do

      allocate(metadata(3+k))
      
      metadata(1) = k
      metadata(2) = blocksize
      metadata(3) = sourcelen
      metadata(4:3+k) = destlen(1:k)
      
      call write_base64(writer, asbytes(c_loc(metadata),4*size(metadata)), 4*size(metadata))
      call write_base64(writer, dest, int(sum(destlen)))
      deallocate(metadata, dest)
      
    end subroutine write_compressed_data

    subroutine write_base64(writer, data, length)
      type(xmlTextWriter) :: writer
      character(kind=c_char) :: data(:)
      integer :: length
      integer :: i, blocksize=24, err, j

      j = length
      i=0
      do while(j .ge. blocksize)
         err = xmlTextWriterWriteBase64(writer, data(blocksize*i+1:), 0, blocksize)
         i = i+1
         j = j-blocksize
      end do

      if (j>0) &
           err = xmlTextWriterWriteBase64(writer, data(blocksize*i+1:), 0, j)

    end subroutine write_base64

    subroutine add_cell_descriptions()

      integer :: err


      err = xmlTextWriterStartElement(writer, "Points")
      call add_binary_data(writer, "Points", "Float64", c_loc(points),&
           3, 8*size(points), compress)
      err = xmlTextWriterEndElement(writer)
      if (ishead()) then
         err = xmlTextWriterStartElement(pwriter, "PPoints")
         err = xmlTextWriterStartElement(pwriter, "PDataArray")
         err = xmlTextWriterWriteAttribute(pwriter, "type", "Float64")
         err = xmlTextWriterWriteAttribute(pwriter, "Name", "Points")
         err = xmlTextWriterWriteAttribute(pwriter, "NumberOfComponents", "3")
         err = xmlTextWriterEndElement(pwriter)
         err = xmlTextWriterEndElement(pwriter)
      end if
      
      err = xmlTextWriterStartElement(writer, "Cells")
      
      !connectivity
      
      call add_binary_data(writer, "connectivity", "Int32", c_loc(connectivity),&
           1, 4*size(connectivity), compress)
      !offsets

      call add_binary_data(writer, "offsets", "Int32", c_loc(offsets),&
           1, 4*size(offsets), compress)
      
      call add_binary_data(writer, "types", "UInt8", c_loc(types),&
           1, size(types), compress)

      err = xmlTextWriterEndElement(writer)

    end subroutine add_cell_descriptions

    subroutine vtkwriteghostlevels(ghost_levels)
      integer, dimension(:) :: ghost_levels

      call vtkwritesc(int(ghost_levels, c_int8_t), "vtkGhostLevels")

    end subroutine vtkwriteghostlevels

#endif
    
end module vtkwriter


