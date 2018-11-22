module xmlwriter

  use iso_c_binding

  use fxmltools

  implicit none


  type xmlTextWriter
     type(c_ptr) :: ptr
  end type xmlTextWriter

  type xmlOutputBuffer
     type(c_ptr) :: ptr
  end type xmlOutputBuffer

  type xmlBuffer
     type(c_ptr) :: ptr
  end type xmlBuffer

  type xmlParserCtxt
     type(c_ptr) :: ptr
  end type xmlParserCtxt

  type xmlDoc
     type(c_ptr) :: ptr
  end type xmlDoc

  type xmlNode
     type(c_ptr) :: ptr
  end type xmlNode


  interface xmlNewTextWriter

     module procedure xmlNewTextWriter_f

     type(c_ptr) function xmlNewTextWriter_c(out) &
          bind(c,name="xmlNewTextWriter")
         use iso_c_binding
         type(c_ptr), value :: out

     end function xmlNewTextWriter_c
  end interface xmlNewTextWriter


  interface xmlNewTextWriterFilename

     module procedure xmlNewTextWriterFilename_f

     type(c_ptr) function xmlNewTextWriterFilename_c(uri, compression) &
          bind(c,name="xmlNewTextWriterFilename")
         use iso_c_binding
         character(c_char), dimension(*) :: uri
         integer(c_int), value :: compression

     end function xmlNewTextWriterFilename_c
  end interface xmlNewTextWriterFilename


  interface xmlNewTextWriterMemory

     module procedure xmlNewTextWriterMemory_f

     type(c_ptr) function xmlNewTextWriterMemory_c(buf, compression) &
          bind(c,name="xmlNewTextWriterMemory")
         use iso_c_binding
         type(c_ptr), value :: buf
         integer(c_int), value :: compression

     end function xmlNewTextWriterMemory_c
  end interface xmlNewTextWriterMemory


  interface xmlNewTextWriterPushParser

     module procedure xmlNewTextWriterPushParser_f

     type(c_ptr) function xmlNewTextWriterPushParser_c(ctxt, compression) &
          bind(c,name="xmlNewTextWriterPushParser")
         use iso_c_binding
         type(c_ptr), value :: ctxt
         integer(c_int), value :: compression

     end function xmlNewTextWriterPushParser_c
  end interface xmlNewTextWriterPushParser


  interface xmlNewTextWriterDoc

     module procedure xmlNewTextWriterDoc_f

     type(c_ptr) function xmlNewTextWriterDoc_c(doc, compression) &
          bind(c,name="xmlNewTextWriterDoc")
         use iso_c_binding
         type(c_ptr) :: doc
         integer(c_int), value :: compression

     end function xmlNewTextWriterDoc_c
  end interface xmlNewTextWriterDoc


  interface xmlNewTextWriterTree

     module procedure xmlNewTextWriterTree_f

     type(c_ptr) function xmlNewTextWriterTree_c(doc, node, compression) &
          bind(c,name="xmlNewTextWriterTree")
         use iso_c_binding
         type(c_ptr), value :: doc
         type(c_ptr), value :: node
         integer(c_int), value :: compression

     end function xmlNewTextWriterTree_c
  end interface xmlNewTextWriterTree


  interface xmlFreeTextWriter

     module procedure xmlFreeTextWriter_f

     subroutine xmlFreeTextWriter_c(writer) &
          bind(c,name="xmlFreeTextWriter")
         use iso_c_binding
         type(c_ptr), value :: writer

     end subroutine xmlFreeTextWriter_c
  end interface xmlFreeTextWriter


  interface xmlTextWriterStartDocument

     module procedure xmlTextWriterStartDocument_f

     integer(c_int) function xmlTextWriterStartDocument_c(writer, version, encoding, standalone) &
          bind(c,name="xmlTextWriterStartDocument")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: version
         character(c_char), dimension(*) :: encoding
         character(c_char), dimension(*) :: standalone

     end function xmlTextWriterStartDocument_c
  end interface xmlTextWriterStartDocument


  interface xmlTextWriterEndDocument

     module procedure xmlTextWriterEndDocument_f

     integer(c_int) function xmlTextWriterEndDocument_c(writer) &
          bind(c,name="xmlTextWriterEndDocument")
         use iso_c_binding
         type(c_ptr), value :: writer

     end function xmlTextWriterEndDocument_c
  end interface xmlTextWriterEndDocument


  interface xmlTextWriterStartComment

     module procedure xmlTextWriterStartComment_f

     integer(c_int) function xmlTextWriterStartComment_c(writer) &
          bind(c,name="xmlTextWriterStartComment")
         use iso_c_binding
         type(c_ptr), value :: writer

     end function xmlTextWriterStartComment_c
  end interface xmlTextWriterStartComment


  interface xmlTextWriterEndComment

     module procedure xmlTextWriterEndComment_f

     integer(c_int) function xmlTextWriterEndComment_c(writer) &
          bind(c,name="xmlTextWriterEndComment")
         use iso_c_binding
         type(c_ptr), value :: writer

     end function xmlTextWriterEndComment_c
  end interface xmlTextWriterEndComment


  interface xmlTextWriterWriteComment

     module procedure xmlTextWriterWriteComment_f

     integer(c_int) function xmlTextWriterWriteComment_c(writer, content) &
          bind(c,name="xmlTextWriterWriteComment")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: content

     end function xmlTextWriterWriteComment_c
  end interface xmlTextWriterWriteComment


  interface xmlTextWriterStartElement

     module procedure xmlTextWriterStartElement_f

     integer(c_int) function xmlTextWriterStartElement_c(writer, name) &
          bind(c,name="xmlTextWriterStartElement")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: name

     end function xmlTextWriterStartElement_c
  end interface xmlTextWriterStartElement


  interface xmlTextWriterStartElementNS

     module procedure xmlTextWriterStartElementNS_f

     integer(c_int) function xmlTextWriterStartElementNS_c(writer, prefix, name, namespaceURI) &
          bind(c,name="xmlTextWriterStartElementNS")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: prefix
         character(c_char), dimension(*) :: name
         character(c_char), dimension(*) :: namespaceURI

     end function xmlTextWriterStartElementNS_c
  end interface xmlTextWriterStartElementNS


  interface xmlTextWriterEndElement

     module procedure xmlTextWriterEndElement_f

     integer(c_int) function xmlTextWriterEndElement_c(writer) &
          bind(c,name="xmlTextWriterEndElement")
         use iso_c_binding
         type(c_ptr), value :: writer

     end function xmlTextWriterEndElement_c
  end interface xmlTextWriterEndElement


  interface xmlTextWriterFullEndElement

     module procedure xmlTextWriterFullEndElement_f

     integer(c_int) function xmlTextWriterFullEndElement_c(writer) &
          bind(c,name="xmlTextWriterFullEndElement")
         use iso_c_binding
         type(c_ptr), value :: writer

     end function xmlTextWriterFullEndElement_c
  end interface xmlTextWriterFullEndElement


  interface xmlTextWriterWriteElement

     module procedure xmlTextWriterWriteElement_f

     integer(c_int) function xmlTextWriterWriteElement_c(writer, name, content) &
          bind(c,name="xmlTextWriterWriteElement")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: name
         character(c_char), dimension(*) :: content

     end function xmlTextWriterWriteElement_c
  end interface xmlTextWriterWriteElement


  interface xmlTextWriterWriteElementNS

     module procedure xmlTextWriterWriteElementNS_f

     integer(c_int) function xmlTextWriterWriteElementNS_c(writer, prefix, name, namespaceURI, content) &
          bind(c,name="xmlTextWriterWriteElementNS")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: prefix
         character(c_char), dimension(*) :: name
         character(c_char), dimension(*) :: namespaceURI
         character(c_char), dimension(*) :: content

     end function xmlTextWriterWriteElementNS_c
  end interface xmlTextWriterWriteElementNS


  interface xmlTextWriterWriteRawLen

     module procedure xmlTextWriterWriteRawLen_f

     integer(c_int) function xmlTextWriterWriteRawLen_c(writer, content, len) &
          bind(c,name="xmlTextWriterWriteRawLen")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: content
         integer(c_int), value :: len

     end function xmlTextWriterWriteRawLen_c
  end interface xmlTextWriterWriteRawLen


  interface xmlTextWriterWriteRaw

     module procedure xmlTextWriterWriteRaw_f

     integer(c_int) function xmlTextWriterWriteRaw_c(writer, content) &
          bind(c,name="xmlTextWriterWriteRaw")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: content

     end function xmlTextWriterWriteRaw_c
  end interface xmlTextWriterWriteRaw


  interface xmlTextWriterWriteString

     module procedure xmlTextWriterWriteString_f

     integer(c_int) function xmlTextWriterWriteString_c(writer, content) &
          bind(c,name="xmlTextWriterWriteString")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: content

     end function xmlTextWriterWriteString_c
  end interface xmlTextWriterWriteString


  interface xmlTextWriterWriteBase64

     module procedure xmlTextWriterWriteBase64_f

     integer(c_int) function xmlTextWriterWriteBase64_c(writer, data, start, len) &
          bind(c,name="xmlTextWriterWriteBase64")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: data
         integer(c_int), value :: start
         integer(c_int), value :: len

     end function xmlTextWriterWriteBase64_c
  end interface xmlTextWriterWriteBase64


  interface xmlTextWriterWriteBinHex

     module procedure xmlTextWriterWriteBinHex_f

     integer(c_int) function xmlTextWriterWriteBinHex_c(writer, data, start, len) &
          bind(c,name="xmlTextWriterWriteBinHex")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: data
         integer(c_int), value :: start
         integer(c_int), value :: len

     end function xmlTextWriterWriteBinHex_c
  end interface xmlTextWriterWriteBinHex


  interface xmlTextWriterStartAttribute

     module procedure xmlTextWriterStartAttribute_f

     integer(c_int) function xmlTextWriterStartAttribute_c(writer, name) &
          bind(c,name="xmlTextWriterStartAttribute")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: name

     end function xmlTextWriterStartAttribute_c
  end interface xmlTextWriterStartAttribute


  interface xmlTextWriterStartAttributeNS

     module procedure xmlTextWriterStartAttributeNS_f

     integer(c_int) function xmlTextWriterStartAttributeNS_c(writer, prefix, name, namespaceURI) &
          bind(c,name="xmlTextWriterStartAttributeNS")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: prefix
         character(c_char), dimension(*) :: name
         character(c_char), dimension(*) :: namespaceURI

     end function xmlTextWriterStartAttributeNS_c
  end interface xmlTextWriterStartAttributeNS


  interface xmlTextWriterEndAttribute

     module procedure xmlTextWriterEndAttribute_f

     integer(c_int) function xmlTextWriterEndAttribute_c(writer) &
          bind(c,name="xmlTextWriterEndAttribute")
         use iso_c_binding
         type(c_ptr), value :: writer

     end function xmlTextWriterEndAttribute_c
  end interface xmlTextWriterEndAttribute


  interface xmlTextWriterWriteAttribute

     module procedure xmlTextWriterWriteAttribute_f

     integer(c_int) function xmlTextWriterWriteAttribute_c(writer, name, content) &
          bind(c,name="xmlTextWriterWriteAttribute")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: name
         character(c_char), dimension(*) :: content

     end function xmlTextWriterWriteAttribute_c
  end interface xmlTextWriterWriteAttribute


  interface xmlTextWriterWriteAttributeNS

     module procedure xmlTextWriterWriteAttributeNS_f

     integer(c_int) function xmlTextWriterWriteAttributeNS_c(writer, prefix, name, namespaceURI, content) &
          bind(c,name="xmlTextWriterWriteAttributeNS")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: prefix
         character(c_char), dimension(*) :: name
         character(c_char), dimension(*) :: namespaceURI
         character(c_char), dimension(*) :: content

     end function xmlTextWriterWriteAttributeNS_c
  end interface xmlTextWriterWriteAttributeNS

  interface xmlTextWriterStartPI

     module procedure xmlTextWriterStartPI_f

     integer(c_int) function xmlTextWriterStartPI_c(writer, target) &
          bind(c,name="xmlTextWriterStartPI")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: target

     end function xmlTextWriterStartPI_c
  end interface xmlTextWriterStartPI


  interface xmlTextWriterEndPI

     module procedure xmlTextWriterEndPI_f

     integer(c_int) function xmlTextWriterEndPI_c(writer) &
          bind(c,name="xmlTextWriterEndPI")
         use iso_c_binding
         type(c_ptr), value :: writer

     end function xmlTextWriterEndPI_c
  end interface xmlTextWriterEndPI


  interface xmlTextWriterWritePI

     module procedure xmlTextWriterWritePI_f

     integer(c_int) function xmlTextWriterWritePI_c(writer, target, content) &
          bind(c,name="xmlTextWriterWritePI")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: target
         character(c_char), dimension(*) :: content

     end function xmlTextWriterWritePI_c
  end interface xmlTextWriterWritePI


  interface xmlTextWriterStartCDATA

     module procedure xmlTextWriterStartCDATA_f

     integer(c_int) function xmlTextWriterStartCDATA_c(writer) &
          bind(c,name="xmlTextWriterStartCDATA")
         use iso_c_binding
         type(c_ptr), value :: writer

     end function xmlTextWriterStartCDATA_c
  end interface xmlTextWriterStartCDATA


  interface xmlTextWriterEndCDATA

     module procedure xmlTextWriterEndCDATA_f

     integer(c_int) function xmlTextWriterEndCDATA_c(writer) &
          bind(c,name="xmlTextWriterEndCDATA")
         use iso_c_binding
         type(c_ptr), value :: writer

     end function xmlTextWriterEndCDATA_c
  end interface xmlTextWriterEndCDATA


  interface xmlTextWriterWriteCDATA

     module procedure xmlTextWriterWriteCDATA_f

     integer(c_int) function xmlTextWriterWriteCDATA_c(writer, content) &
          bind(c,name="xmlTextWriterWriteCDATA")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: content

     end function xmlTextWriterWriteCDATA_c
  end interface xmlTextWriterWriteCDATA


  interface xmlTextWriterStartDTD

     module procedure xmlTextWriterStartDTD_f

     integer(c_int) function xmlTextWriterStartDTD_c(writer, name, pubid, sysid) &
          bind(c,name="xmlTextWriterStartDTD")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: name
         character(c_char), dimension(*) :: pubid
         character(c_char), dimension(*) :: sysid

     end function xmlTextWriterStartDTD_c
  end interface xmlTextWriterStartDTD


  interface xmlTextWriterEndDTD

     module procedure xmlTextWriterEndDTD_f

     integer(c_int) function xmlTextWriterEndDTD_c(writer) &
          bind(c,name="xmlTextWriterEndDTD")
         use iso_c_binding
         type(c_ptr), value :: writer

     end function xmlTextWriterEndDTD_c
  end interface xmlTextWriterEndDTD


  interface xmlTextWriterWriteDTD

     module procedure xmlTextWriterWriteDTD_f

     integer(c_int) function xmlTextWriterWriteDTD_c(writer, name, pubid, sysid, subset) &
          bind(c,name="xmlTextWriterWriteDTD")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: name
         character(c_char), dimension(*) :: pubid
         character(c_char), dimension(*) :: sysid
         character(c_char), dimension(*) :: subset

     end function xmlTextWriterWriteDTD_c
  end interface xmlTextWriterWriteDTD


  interface xmlTextWriterStartDTDElement

     module procedure xmlTextWriterStartDTDElement_f

     integer(c_int) function xmlTextWriterStartDTDElement_c(writer, name) &
          bind(c,name="xmlTextWriterStartDTDElement")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: name

     end function xmlTextWriterStartDTDElement_c
  end interface xmlTextWriterStartDTDElement


  interface xmlTextWriterEndDTDElement

     module procedure xmlTextWriterEndDTDElement_f

     integer(c_int) function xmlTextWriterEndDTDElement_c(writer) &
          bind(c,name="xmlTextWriterEndDTDElement")
         use iso_c_binding
         type(c_ptr), value :: writer

     end function xmlTextWriterEndDTDElement_c
  end interface xmlTextWriterEndDTDElement


  interface xmlTextWriterWriteDTDElement

     module procedure xmlTextWriterWriteDTDElement_f

     integer(c_int) function xmlTextWriterWriteDTDElement_c(writer, name, content) &
          bind(c,name="xmlTextWriterWriteDTDElement")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: name
         character(c_char), dimension(*) :: content

     end function xmlTextWriterWriteDTDElement_c
  end interface xmlTextWriterWriteDTDElement


  interface xmlTextWriterStartDTDAttlist

     module procedure xmlTextWriterStartDTDAttlist_f

     integer(c_int) function xmlTextWriterStartDTDAttlist_c(writer, name) &
          bind(c,name="xmlTextWriterStartDTDAttlist")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: name

     end function xmlTextWriterStartDTDAttlist_c
  end interface xmlTextWriterStartDTDAttlist


  interface xmlTextWriterEndDTDAttlist

     module procedure xmlTextWriterEndDTDAttlist_f

     integer(c_int) function xmlTextWriterEndDTDAttlist_c(writer) &
          bind(c,name="xmlTextWriterEndDTDAttlist")
         use iso_c_binding
         type(c_ptr), value :: writer

     end function xmlTextWriterEndDTDAttlist_c
  end interface xmlTextWriterEndDTDAttlist


  interface xmlTextWriterWriteDTDAttlist

     module procedure xmlTextWriterWriteDTDAttlist_f

     integer(c_int) function xmlTextWriterWriteDTDAttlist_c(writer, name, content) &
          bind(c,name="xmlTextWriterWriteDTDAttlist")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: name
         character(c_char), dimension(*) :: content

     end function xmlTextWriterWriteDTDAttlist_c
  end interface xmlTextWriterWriteDTDAttlist


  interface xmlTextWriterStartDTDEntity

     module procedure xmlTextWriterStartDTDEntity_f

     integer(c_int) function xmlTextWriterStartDTDEntity_c(writer, pe, name) &
          bind(c,name="xmlTextWriterStartDTDEntity")
         use iso_c_binding
         type(c_ptr), value :: writer
         integer(c_int), value :: pe
         character(c_char), dimension(*) :: name

     end function xmlTextWriterStartDTDEntity_c
  end interface xmlTextWriterStartDTDEntity


  interface xmlTextWriterEndDTDEntity

     module procedure xmlTextWriterEndDTDEntity_f

     integer(c_int) function xmlTextWriterEndDTDEntity_c(writer) &
          bind(c,name="xmlTextWriterEndDTDEntity")
         use iso_c_binding
         type(c_ptr), value :: writer

     end function xmlTextWriterEndDTDEntity_c
  end interface xmlTextWriterEndDTDEntity


  interface xmlTextWriterWriteDTDInternalEntity

     module procedure xmlTextWriterWriteDTDInternalEntity_f

     integer(c_int) function xmlTextWriterWriteDTDInternalEntity_c(writer, pe, name, content) &
          bind(c,name="xmlTextWriterWriteDTDInternalEntity")
         use iso_c_binding
         type(c_ptr), value :: writer
         integer(c_int), value :: pe
         character(c_char), dimension(*) :: name
         character(c_char), dimension(*) :: content

     end function xmlTextWriterWriteDTDInternalEntity_c
  end interface xmlTextWriterWriteDTDInternalEntity


  interface xmlTextWriterWriteDTDExternalEntity

     module procedure xmlTextWriterWriteDTDExternalEntity_f

     integer(c_int) function xmlTextWriterWriteDTDExternalEntity_c(writer, pe, name, pubid, sysid, ndataid) &
          bind(c,name="xmlTextWriterWriteDTDExternalEntity")
         use iso_c_binding
         type(c_ptr), value :: writer
         integer(c_int), value :: pe
         character(c_char), dimension(*) :: name
         character(c_char), dimension(*) :: pubid
         character(c_char), dimension(*) :: sysid
         character(c_char), dimension(*) :: ndataid

     end function xmlTextWriterWriteDTDExternalEntity_c
  end interface xmlTextWriterWriteDTDExternalEntity


  interface xmlTextWriterWriteDTDExternalEntityContents

     module procedure xmlTextWriterWriteDTDExternalEntityContents_f

     integer(c_int) function xmlTextWriterWriteDTDExternalEntityContents_c(writer, pubid, sysid, ndataid) &
          bind(c,name="xmlTextWriterWriteDTDExternalEntityContents")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: pubid
         character(c_char), dimension(*) :: sysid
         character(c_char), dimension(*) :: ndataid

     end function xmlTextWriterWriteDTDExternalEntityContents_c
  end interface xmlTextWriterWriteDTDExternalEntityContents


  interface xmlTextWriterWriteDTDEntity

     module procedure xmlTextWriterWriteDTDEntity_f

     integer(c_int) function xmlTextWriterWriteDTDEntity_c(writer, pe, name, pubid, sysid, ndataid, content) &
          bind(c,name="xmlTextWriterWriteDTDEntity")
         use iso_c_binding
         type(c_ptr), value :: writer
         integer(c_int), value :: pe
         character(c_char), dimension(*) :: name
         character(c_char), dimension(*) :: pubid
         character(c_char), dimension(*) :: sysid
         character(c_char), dimension(*) :: ndataid
         character(c_char), dimension(*) :: content

     end function xmlTextWriterWriteDTDEntity_c
  end interface xmlTextWriterWriteDTDEntity


  interface xmlTextWriterWriteDTDNotation

     module procedure xmlTextWriterWriteDTDNotation_f

     integer(c_int) function xmlTextWriterWriteDTDNotation_c(writer, name, pubid, sysid) &
          bind(c,name="xmlTextWriterWriteDTDNotation")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: name
         character(c_char), dimension(*) :: pubid
         character(c_char), dimension(*) :: sysid

     end function xmlTextWriterWriteDTDNotation_c
  end interface xmlTextWriterWriteDTDNotation


  interface xmlTextWriterSetIndent

     module procedure xmlTextWriterSetIndent_f

     integer(c_int) function xmlTextWriterSetIndent_c(writer, indent) &
          bind(c,name="xmlTextWriterSetIndent")
         use iso_c_binding
         type(c_ptr), value :: writer
         integer(c_int), value :: indent

     end function xmlTextWriterSetIndent_c
  end interface xmlTextWriterSetIndent


  interface xmlTextWriterSetIndentString

     module procedure xmlTextWriterSetIndentString_f

     integer(c_int) function xmlTextWriterSetIndentString_c(writer, str) &
          bind(c,name="xmlTextWriterSetIndentString")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), dimension(*) :: str

     end function xmlTextWriterSetIndentString_c
  end interface xmlTextWriterSetIndentString


  interface xmlTextWriterSetQuoteChar

     module procedure xmlTextWriterSetQuoteChar_f

     integer(c_int) function xmlTextWriterSetQuoteChar_c(writer, quotechar) &
          bind(c,name="xmlTextWriterSetQuoteChar")
         use iso_c_binding
         type(c_ptr), value :: writer
         character(c_char), value :: quotechar

     end function xmlTextWriterSetQuoteChar_c
  end interface xmlTextWriterSetQuoteChar


  interface xmlTextWriterFlush

     module procedure xmlTextWriterFlush_f

     integer(c_int) function xmlTextWriterFlush_c(writer) &
          bind(c,name="xmlTextWriterFlush")
         use iso_c_binding
         type(c_ptr), value :: writer

     end function xmlTextWriterFlush_c
  end interface xmlTextWriterFlush


  contains



    type(xmlTextWriter) function xmlNewTextWriter_f(out) result(fout)
      type(xmlOutputBuffer) :: out

      fout%ptr=xmlNewTextWriter_c(out%ptr)
    end function xmlNewTextWriter_f


    type(xmlTextWriter) function xmlNewTextWriterFilename_f(uri, compression) result(fout)
      character(len=*) :: uri
      integer :: compression

      fout%ptr=xmlNewTextWriterFilename_c(c_wrap(uri),compression)
    end function xmlNewTextWriterFilename_f


    type(xmlTextWriter) function xmlNewTextWriterMemory_f(buf, compression) result(fout)
      type(xmlBuffer) :: buf
      integer :: compression

      fout%ptr=xmlNewTextWriterMemory_c(buf%ptr,compression)
    end function xmlNewTextWriterMemory_f


    type(xmlTextWriter) function xmlNewTextWriterPushParser_f(ctxt, compression) result(fout)
      type(xmlParserCtxt) :: ctxt
      integer :: compression

      fout%ptr=xmlNewTextWriterPushParser_c(ctxt%ptr,compression)
    end function xmlNewTextWriterPushParser_f


    type(xmlTextWriter) function xmlNewTextWriterDoc_f(doc, compression) result(fout)
      type(xmlDoc) :: doc
      integer :: compression

      fout%ptr=xmlNewTextWriterDoc_c(doc%ptr,compression)
    end function xmlNewTextWriterDoc_f


    type(xmlTextWriter) function xmlNewTextWriterTree_f(doc, node, compression) result(fout)
      type(xmlDoc) :: doc
      type(xmlNode) :: node
      integer :: compression

      fout%ptr=xmlNewTextWriterTree_c(doc%ptr,node%ptr,compression)
    end function xmlNewTextWriterTree_f


    subroutine xmlFreeTextWriter_f(writer)
      type(xmlTextWriter) :: writer


      call xmlFreeTextWriter_c(writer%ptr)
    end subroutine xmlFreeTextWriter_f


    integer function xmlTextWriterStartDocument_f(writer, version, encoding, standalone) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: version
      character(len=*) :: encoding
      character(len=*) :: standalone

      fout=xmlTextWriterStartDocument_c(writer%ptr,c_wrap(version),c_wrap(encoding),c_wrap(standalone))
    end function xmlTextWriterStartDocument_f


    integer function xmlTextWriterEndDocument_f(writer) result(fout)
      type(xmlTextWriter) :: writer

      fout=xmlTextWriterEndDocument_c(writer%ptr)
    end function xmlTextWriterEndDocument_f


    integer function xmlTextWriterStartComment_f(writer) result(fout)
      type(xmlTextWriter) :: writer

      fout=xmlTextWriterStartComment_c(writer%ptr)
    end function xmlTextWriterStartComment_f


    integer function xmlTextWriterEndComment_f(writer) result(fout)
      type(xmlTextWriter) :: writer

      fout=xmlTextWriterEndComment_c(writer%ptr)
    end function xmlTextWriterEndComment_f


    integer function xmlTextWriterWriteComment_f(writer, content) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: content

      fout=xmlTextWriterWriteComment_c(writer%ptr,c_wrap(content))
    end function xmlTextWriterWriteComment_f


    integer function xmlTextWriterStartElement_f(writer, name) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: name

      fout=xmlTextWriterStartElement_c(writer%ptr,c_wrap(name))
    end function xmlTextWriterStartElement_f


    integer function xmlTextWriterStartElementNS_f(writer, prefix, name, namespaceURI) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: prefix
      character(len=*) :: name
      character(len=*) :: namespaceURI

      fout=xmlTextWriterStartElementNS_c(writer%ptr,c_wrap(prefix),c_wrap(name),c_wrap(namespaceURI))
    end function xmlTextWriterStartElementNS_f


    integer function xmlTextWriterEndElement_f(writer) result(fout)
      type(xmlTextWriter) :: writer

      fout=xmlTextWriterEndElement_c(writer%ptr)
    end function xmlTextWriterEndElement_f


    integer function xmlTextWriterFullEndElement_f(writer) result(fout)
      type(xmlTextWriter) :: writer

      fout=xmlTextWriterFullEndElement_c(writer%ptr)
    end function xmlTextWriterFullEndElement_f


    integer function xmlTextWriterWriteElement_f(writer, name, content) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: name
      character(len=*) :: content

      fout=xmlTextWriterWriteElement_c(writer%ptr,c_wrap(name),c_wrap(content))
    end function xmlTextWriterWriteElement_f


    integer function xmlTextWriterWriteElementNS_f(writer, prefix, name, namespaceURI, content) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: prefix
      character(len=*) :: name
      character(len=*) :: namespaceURI
      character(len=*) :: content

      fout=xmlTextWriterWriteElementNS_c(writer%ptr,c_wrap(prefix),c_wrap(name),c_wrap(namespaceURI),c_wrap(content))
    end function xmlTextWriterWriteElementNS_f


    integer function xmlTextWriterWriteRawLen_f(writer, content, len) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: content
      integer :: len

      fout=xmlTextWriterWriteRawLen_c(writer%ptr,c_wrap(content),len)
    end function xmlTextWriterWriteRawLen_f


    integer function xmlTextWriterWriteRaw_f(writer, content) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: content

      fout=xmlTextWriterWriteRaw_c(writer%ptr,c_wrap(content))
    end function xmlTextWriterWriteRaw_f


    integer function xmlTextWriterWriteString_f(writer, content) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: content

      fout=xmlTextWriterWriteString_c(writer%ptr,c_wrap(content))
    end function xmlTextWriterWriteString_f


    integer function xmlTextWriterWriteBase64_f(writer, data, start, len) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: data
      integer :: start
      integer :: len

      fout=xmlTextWriterWriteBase64_c(writer%ptr,c_wrap(data),start,len)
    end function xmlTextWriterWriteBase64_f


    integer function xmlTextWriterWriteBinHex_f(writer, data, start, len) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: data
      integer :: start
      integer :: len

      fout=xmlTextWriterWriteBinHex_c(writer%ptr,c_wrap(data),start,len)
    end function xmlTextWriterWriteBinHex_f


    integer function xmlTextWriterStartAttribute_f(writer, name) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: name

      fout=xmlTextWriterStartAttribute_c(writer%ptr,c_wrap(name))
    end function xmlTextWriterStartAttribute_f


    integer function xmlTextWriterStartAttributeNS_f(writer, prefix, name, namespaceURI) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: prefix
      character(len=*) :: name
      character(len=*) :: namespaceURI

      fout=xmlTextWriterStartAttributeNS_c(writer%ptr,c_wrap(prefix),c_wrap(name),c_wrap(namespaceURI))
    end function xmlTextWriterStartAttributeNS_f


    integer function xmlTextWriterEndAttribute_f(writer) result(fout)
      type(xmlTextWriter) :: writer

      fout=xmlTextWriterEndAttribute_c(writer%ptr)
    end function xmlTextWriterEndAttribute_f


    integer function xmlTextWriterWriteAttribute_f(writer, name, content) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: name
      character(len=*) :: content

      fout=xmlTextWriterWriteAttribute_c(writer%ptr,c_wrap(name),c_wrap(content))
    end function xmlTextWriterWriteAttribute_f


    integer function xmlTextWriterWriteAttributeNS_f(writer, prefix, name, namespaceURI, content) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: prefix
      character(len=*) :: name
      character(len=*) :: namespaceURI
      character(len=*) :: content

      fout=xmlTextWriterWriteAttributeNS_c(writer%ptr,c_wrap(prefix),c_wrap(name),c_wrap(namespaceURI),c_wrap(content))
    end function xmlTextWriterWriteAttributeNS_f

    integer function xmlTextWriterStartPI_f(writer, target) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: target

      fout=xmlTextWriterStartPI_c(writer%ptr,c_wrap(target))
    end function xmlTextWriterStartPI_f


    integer function xmlTextWriterEndPI_f(writer) result(fout)
      type(xmlTextWriter) :: writer

      fout=xmlTextWriterEndPI_c(writer%ptr)
    end function xmlTextWriterEndPI_f


    integer function xmlTextWriterWritePI_f(writer, target, content) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: target
      character(len=*) :: content

      fout=xmlTextWriterWritePI_c(writer%ptr,c_wrap(target),c_wrap(content))
    end function xmlTextWriterWritePI_f


    integer function xmlTextWriterStartCDATA_f(writer) result(fout)
      type(xmlTextWriter) :: writer

      fout=xmlTextWriterStartCDATA_c(writer%ptr)
    end function xmlTextWriterStartCDATA_f


    integer function xmlTextWriterEndCDATA_f(writer) result(fout)
      type(xmlTextWriter) :: writer

      fout=xmlTextWriterEndCDATA_c(writer%ptr)
    end function xmlTextWriterEndCDATA_f


    integer function xmlTextWriterWriteCDATA_f(writer, content) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: content

      fout=xmlTextWriterWriteCDATA_c(writer%ptr,c_wrap(content))
    end function xmlTextWriterWriteCDATA_f


    integer function xmlTextWriterStartDTD_f(writer, name, pubid, sysid) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: name
      character(len=*) :: pubid
      character(len=*) :: sysid

      fout=xmlTextWriterStartDTD_c(writer%ptr,c_wrap(name),c_wrap(pubid),c_wrap(sysid))
    end function xmlTextWriterStartDTD_f


    integer function xmlTextWriterEndDTD_f(writer) result(fout)
      type(xmlTextWriter) :: writer

      fout=xmlTextWriterEndDTD_c(writer%ptr)
    end function xmlTextWriterEndDTD_f


    integer function xmlTextWriterWriteDTD_f(writer, name, pubid, sysid, subset) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: name
      character(len=*) :: pubid
      character(len=*) :: sysid
      character(len=*) :: subset

      fout=xmlTextWriterWriteDTD_c(writer%ptr,c_wrap(name),c_wrap(pubid),c_wrap(sysid),c_wrap(subset))
    end function xmlTextWriterWriteDTD_f


    integer function xmlTextWriterStartDTDElement_f(writer, name) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: name

      fout=xmlTextWriterStartDTDElement_c(writer%ptr,c_wrap(name))
    end function xmlTextWriterStartDTDElement_f


    integer function xmlTextWriterEndDTDElement_f(writer) result(fout)
      type(xmlTextWriter) :: writer

      fout=xmlTextWriterEndDTDElement_c(writer%ptr)
    end function xmlTextWriterEndDTDElement_f


    integer function xmlTextWriterWriteDTDElement_f(writer, name, content) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: name
      character(len=*) :: content

      fout=xmlTextWriterWriteDTDElement_c(writer%ptr,c_wrap(name),c_wrap(content))
    end function xmlTextWriterWriteDTDElement_f


    integer function xmlTextWriterStartDTDAttlist_f(writer, name) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: name

      fout=xmlTextWriterStartDTDAttlist_c(writer%ptr,c_wrap(name))
    end function xmlTextWriterStartDTDAttlist_f


    integer function xmlTextWriterEndDTDAttlist_f(writer) result(fout)
      type(xmlTextWriter) :: writer

      fout=xmlTextWriterEndDTDAttlist_c(writer%ptr)
    end function xmlTextWriterEndDTDAttlist_f


    integer function xmlTextWriterWriteDTDAttlist_f(writer, name, content) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: name
      character(len=*) :: content

      fout=xmlTextWriterWriteDTDAttlist_c(writer%ptr,c_wrap(name),c_wrap(content))
    end function xmlTextWriterWriteDTDAttlist_f


    integer function xmlTextWriterStartDTDEntity_f(writer, pe, name) result(fout)
      type(xmlTextWriter) :: writer
      integer :: pe
      character(len=*) :: name

      fout=xmlTextWriterStartDTDEntity_c(writer%ptr,pe,c_wrap(name))
    end function xmlTextWriterStartDTDEntity_f


    integer function xmlTextWriterEndDTDEntity_f(writer) result(fout)
      type(xmlTextWriter) :: writer

      fout=xmlTextWriterEndDTDEntity_c(writer%ptr)
    end function xmlTextWriterEndDTDEntity_f


    integer function xmlTextWriterWriteDTDInternalEntity_f(writer, pe, name, content) result(fout)
      type(xmlTextWriter) :: writer
      integer :: pe
      character(len=*) :: name
      character(len=*) :: content

      fout=xmlTextWriterWriteDTDInternalEntity_c(writer%ptr,pe,c_wrap(name),c_wrap(content))
    end function xmlTextWriterWriteDTDInternalEntity_f


    integer function xmlTextWriterWriteDTDExternalEntity_f(writer, pe, name, pubid, sysid, ndataid) result(fout)
      type(xmlTextWriter) :: writer
      integer :: pe
      character(len=*) :: name
      character(len=*) :: pubid
      character(len=*) :: sysid
      character(len=*) :: ndataid

      fout=xmlTextWriterWriteDTDExternalEntity_c(writer%ptr,pe,c_wrap(name),c_wrap(pubid),c_wrap(sysid),c_wrap(ndataid))
    end function xmlTextWriterWriteDTDExternalEntity_f


    integer function xmlTextWriterWriteDTDExternalEntityContents_f(writer, pubid, sysid, ndataid) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: pubid
      character(len=*) :: sysid
      character(len=*) :: ndataid

      fout=xmlTextWriterWriteDTDExternalEntityContents_c(writer%ptr,c_wrap(pubid),c_wrap(sysid),c_wrap(ndataid))
    end function xmlTextWriterWriteDTDExternalEntityContents_f


    integer function xmlTextWriterWriteDTDEntity_f(writer, pe, name, pubid, sysid, ndataid, content) result(fout)
      type(xmlTextWriter) :: writer
      integer :: pe
      character(len=*) :: name
      character(len=*) :: pubid
      character(len=*) :: sysid
      character(len=*) :: ndataid
      character(len=*) :: content

      fout=xmlTextWriterWriteDTDEntity_c(writer%ptr,pe,c_wrap(name),c_wrap(pubid),c_wrap(sysid),c_wrap(ndataid),c_wrap(content))
    end function xmlTextWriterWriteDTDEntity_f


    integer function xmlTextWriterWriteDTDNotation_f(writer, name, pubid, sysid) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: name
      character(len=*) :: pubid
      character(len=*) :: sysid

      fout=xmlTextWriterWriteDTDNotation_c(writer%ptr,c_wrap(name),c_wrap(pubid),c_wrap(sysid))
    end function xmlTextWriterWriteDTDNotation_f


    integer function xmlTextWriterSetIndent_f(writer, indent) result(fout)
      type(xmlTextWriter) :: writer
      integer :: indent

      fout=xmlTextWriterSetIndent_c(writer%ptr,indent)
    end function xmlTextWriterSetIndent_f


    integer function xmlTextWriterSetIndentString_f(writer, str) result(fout)
      type(xmlTextWriter) :: writer
      character(len=*) :: str

      fout=xmlTextWriterSetIndentString_c(writer%ptr,c_wrap(str))
    end function xmlTextWriterSetIndentString_f


    integer function xmlTextWriterSetQuoteChar_f(writer, quotechar) result(fout)
      type(xmlTextWriter) :: writer
      character :: quotechar

      fout=xmlTextWriterSetQuoteChar_c(writer%ptr,quotechar)
    end function xmlTextWriterSetQuoteChar_f


    integer function xmlTextWriterFlush_f(writer) result(fout)
      type(xmlTextWriter) :: writer

      fout=xmlTextWriterFlush_c(writer%ptr)
    end function xmlTextWriterFlush_f



end module xmlwriter
