module xmlreader

  use iso_c_binding

  use fxmltools
  use xmldatatypes

  implicit none

  private

  public :: xmlTextReader, xmlReaderForFile, xmlTextReaderMoveToAttribute, &
       xmlTextReaderConstValue, xmlTextReaderRead, xmlTextReaderExpand, &
       xmlTextReaderConstName, xmlTextReaderNodeType, xmlFreeTextReader

  interface xmlNewTextReader

     module procedure xmlNewTextReader_f

     type(c_ptr) function xmlNewTextReader_c(input, URI) &
          bind(c,name="xmlNewTextReader")
         use iso_c_binding
         type(c_ptr), value :: input
         character(c_char), dimension(*) :: URI

     end function xmlNewTextReader_c
  end interface xmlNewTextReader

  interface xmlNewTextReaderFilename

     module procedure xmlNewTextReaderFilename_f

     type(c_ptr) function xmlNewTextReaderFilename_c(URI) &
          bind(c,name="xmlNewTextReaderFilename")
         use iso_c_binding
         character(c_char), dimension(*) :: URI

     end function xmlNewTextReaderFilename_c
  end interface xmlNewTextReaderFilename

  interface xmlFreeTextReader

     module procedure xmlFreeTextReader_f

     subroutine xmlFreeTextReader_c(reader) &
          bind(c,name="xmlFreeTextReader")
         use iso_c_binding
         type(c_ptr), value :: reader

     end subroutine xmlFreeTextReader_c
  end interface xmlFreeTextReader


  interface xmlTextReaderSetup

     module procedure xmlTextReaderSetup_f

     integer(c_int) function xmlTextReaderSetup_c(reader, input, URL, encoding, options) &
          bind(c,name="xmlTextReaderSetup")
         use iso_c_binding
         type(c_ptr), value :: reader
         type(c_ptr), value :: input
         character(c_char), dimension(*) :: URL
         character(c_char), dimension(*) :: encoding
         integer(c_int), value :: options

     end function xmlTextReaderSetup_c
  end interface xmlTextReaderSetup

  interface xmlTextReaderRead

     module procedure xmlTextReaderRead_f

     integer(c_int) function xmlTextReaderRead_c(reader) &
          bind(c,name="xmlTextReaderRead")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderRead_c
  end interface xmlTextReaderRead

  interface xmlTextReaderReadInnerXml

     module procedure xmlTextReaderReadInnerXml_f

     type(c_ptr) function xmlTextReaderReadInnerXml_c(reader) &
          bind(c,name="xmlTextReaderReadInnerXml")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderReadInnerXml_c
  end interface xmlTextReaderReadInnerXml

  interface xmlTextReaderReadOuterXml

     module procedure xmlTextReaderReadOuterXml_f

     type(c_ptr) function xmlTextReaderReadOuterXml_c(reader) &
          bind(c,name="xmlTextReaderReadOuterXml")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderReadOuterXml_c
  end interface xmlTextReaderReadOuterXml

  interface xmlTextReaderReadString

     module procedure xmlTextReaderReadString_f

     type(c_ptr) function xmlTextReaderReadString_c(reader) &
          bind(c,name="xmlTextReaderReadString")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderReadString_c
  end interface xmlTextReaderReadString

  interface xmlTextReaderReadAttributeValue

     module procedure xmlTextReaderReadAttributeValue_f

     integer(c_int) function xmlTextReaderReadAttributeValue_c(reader) &
          bind(c,name="xmlTextReaderReadAttributeValue")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderReadAttributeValue_c
  end interface xmlTextReaderReadAttributeValue

  interface xmlTextReaderAttributeCount

     module procedure xmlTextReaderAttributeCount_f

     integer(c_int) function xmlTextReaderAttributeCount_c(reader) &
          bind(c,name="xmlTextReaderAttributeCount")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderAttributeCount_c
  end interface xmlTextReaderAttributeCount

  interface xmlTextReaderDepth

     module procedure xmlTextReaderDepth_f

     integer(c_int) function xmlTextReaderDepth_c(reader) &
          bind(c,name="xmlTextReaderDepth")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderDepth_c
  end interface xmlTextReaderDepth

  interface xmlTextReaderHasAttributes

     module procedure xmlTextReaderHasAttributes_f

     integer(c_int) function xmlTextReaderHasAttributes_c(reader) &
          bind(c,name="xmlTextReaderHasAttributes")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderHasAttributes_c
  end interface xmlTextReaderHasAttributes

  interface xmlTextReaderHasValue

     module procedure xmlTextReaderHasValue_f

     integer(c_int) function xmlTextReaderHasValue_c(reader) &
          bind(c,name="xmlTextReaderHasValue")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderHasValue_c
  end interface xmlTextReaderHasValue

  interface xmlTextReaderIsDefault

     module procedure xmlTextReaderIsDefault_f

     integer(c_int) function xmlTextReaderIsDefault_c(reader) &
          bind(c,name="xmlTextReaderIsDefault")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderIsDefault_c
  end interface xmlTextReaderIsDefault

  interface xmlTextReaderIsEmptyElement

     module procedure xmlTextReaderIsEmptyElement_f

     integer(c_int) function xmlTextReaderIsEmptyElement_c(reader) &
          bind(c,name="xmlTextReaderIsEmptyElement")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderIsEmptyElement_c
  end interface xmlTextReaderIsEmptyElement

  interface xmlTextReaderNodeType

     module procedure xmlTextReaderNodeType_f

     integer(c_int) function xmlTextReaderNodeType_c(reader) &
          bind(c,name="xmlTextReaderNodeType")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderNodeType_c
  end interface xmlTextReaderNodeType

  interface xmlTextReaderQuoteChar

     module procedure xmlTextReaderQuoteChar_f

     integer(c_int) function xmlTextReaderQuoteChar_c(reader) &
          bind(c,name="xmlTextReaderQuoteChar")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderQuoteChar_c
  end interface xmlTextReaderQuoteChar

  interface xmlTextReaderReadState

     module procedure xmlTextReaderReadState_f

     integer(c_int) function xmlTextReaderReadState_c(reader) &
          bind(c,name="xmlTextReaderReadState")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderReadState_c
  end interface xmlTextReaderReadState

  interface xmlTextReaderIsNamespaceDecl

     module procedure xmlTextReaderIsNamespaceDecl_f

     integer(c_int) function xmlTextReaderIsNamespaceDecl_c(reader) &
          bind(c,name="xmlTextReaderIsNamespaceDecl")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderIsNamespaceDecl_c
  end interface xmlTextReaderIsNamespaceDecl

  interface xmlTextReaderConstBaseUri

     module procedure xmlTextReaderConstBaseUri_f

     type(c_ptr) function xmlTextReaderConstBaseUri_c(reader) &
          bind(c,name="xmlTextReaderConstBaseUri")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderConstBaseUri_c
  end interface xmlTextReaderConstBaseUri

  interface xmlTextReaderConstLocalName

     module procedure xmlTextReaderConstLocalName_f

     type(c_ptr) function xmlTextReaderConstLocalName_c(reader) &
          bind(c,name="xmlTextReaderConstLocalName")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderConstLocalName_c
  end interface xmlTextReaderConstLocalName

  interface xmlTextReaderConstName

     module procedure xmlTextReaderConstName_f

     type(c_ptr) function xmlTextReaderConstName_c(reader) &
          bind(c,name="xmlTextReaderConstName")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderConstName_c
  end interface xmlTextReaderConstName

  interface xmlTextReaderConstNamespaceUri

     module procedure xmlTextReaderConstNamespaceUri_f

     type(c_ptr) function xmlTextReaderConstNamespaceUri_c(reader) &
          bind(c,name="xmlTextReaderConstNamespaceUri")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderConstNamespaceUri_c
  end interface xmlTextReaderConstNamespaceUri

  interface xmlTextReaderConstPrefix

     module procedure xmlTextReaderConstPrefix_f

     type(c_ptr) function xmlTextReaderConstPrefix_c(reader) &
          bind(c,name="xmlTextReaderConstPrefix")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderConstPrefix_c
  end interface xmlTextReaderConstPrefix

  interface xmlTextReaderConstXmlLang

     module procedure xmlTextReaderConstXmlLang_f

     type(c_ptr) function xmlTextReaderConstXmlLang_c(reader) &
          bind(c,name="xmlTextReaderConstXmlLang")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderConstXmlLang_c
  end interface xmlTextReaderConstXmlLang

  interface xmlTextReaderConstString

     module procedure xmlTextReaderConstString_f

     type(c_ptr) function xmlTextReaderConstString_c(reader, str) &
          bind(c,name="xmlTextReaderConstString")
         use iso_c_binding
         type(c_ptr), value :: reader
         character(c_char), dimension(*) :: str

     end function xmlTextReaderConstString_c
  end interface xmlTextReaderConstString

  interface xmlTextReaderConstValue

     module procedure xmlTextReaderConstValue_f

     type(c_ptr) function xmlTextReaderConstValue_c(reader) &
          bind(c,name="xmlTextReaderConstValue")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderConstValue_c
  end interface xmlTextReaderConstValue

  interface xmlTextReaderBaseUri

     module procedure xmlTextReaderBaseUri_f

     type(c_ptr) function xmlTextReaderBaseUri_c(reader) &
          bind(c,name="xmlTextReaderBaseUri")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderBaseUri_c
  end interface xmlTextReaderBaseUri

  interface xmlTextReaderLocalName

     module procedure xmlTextReaderLocalName_f

     type(c_ptr) function xmlTextReaderLocalName_c(reader) &
          bind(c,name="xmlTextReaderLocalName")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderLocalName_c
  end interface xmlTextReaderLocalName

  interface xmlTextReaderName

     module procedure xmlTextReaderName_f

     type(c_ptr) function xmlTextReaderName_c(reader) &
          bind(c,name="xmlTextReaderName")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderName_c
  end interface xmlTextReaderName

  interface xmlTextReaderNamespaceUri

     module procedure xmlTextReaderNamespaceUri_f

     type(c_ptr) function xmlTextReaderNamespaceUri_c(reader) &
          bind(c,name="xmlTextReaderNamespaceUri")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderNamespaceUri_c
  end interface xmlTextReaderNamespaceUri

  interface xmlTextReaderPrefix

     module procedure xmlTextReaderPrefix_f

     type(c_ptr) function xmlTextReaderPrefix_c(reader) &
          bind(c,name="xmlTextReaderPrefix")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderPrefix_c
  end interface xmlTextReaderPrefix

  interface xmlTextReaderXmlLang

     module procedure xmlTextReaderXmlLang_f

     type(c_ptr) function xmlTextReaderXmlLang_c(reader) &
          bind(c,name="xmlTextReaderXmlLang")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderXmlLang_c
  end interface xmlTextReaderXmlLang

  interface xmlTextReaderValue

     module procedure xmlTextReaderValue_f

     type(c_ptr) function xmlTextReaderValue_c(reader) &
          bind(c,name="xmlTextReaderValue")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderValue_c
  end interface xmlTextReaderValue

  interface xmlTextReaderClose

     module procedure xmlTextReaderClose_f

     integer(c_int) function xmlTextReaderClose_c(reader) &
          bind(c,name="xmlTextReaderClose")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderClose_c
  end interface xmlTextReaderClose

  interface xmlTextReaderGetAttributeNo

     module procedure xmlTextReaderGetAttributeNo_f

     type(c_ptr) function xmlTextReaderGetAttributeNo_c(reader, no) &
          bind(c,name="xmlTextReaderGetAttributeNo")
         use iso_c_binding
         type(c_ptr), value :: reader
         integer(c_int), value :: no

     end function xmlTextReaderGetAttributeNo_c
  end interface xmlTextReaderGetAttributeNo

  interface xmlTextReaderGetAttribute

     module procedure xmlTextReaderGetAttribute_f

     type(c_ptr) function xmlTextReaderGetAttribute_c(reader, name) &
          bind(c,name="xmlTextReaderGetAttribute")
         use iso_c_binding
         type(c_ptr), value :: reader
         character(c_char), dimension(*) :: name

     end function xmlTextReaderGetAttribute_c
  end interface xmlTextReaderGetAttribute

  interface xmlTextReaderGetAttributeNs

     module procedure xmlTextReaderGetAttributeNs_f

     type(c_ptr) function xmlTextReaderGetAttributeNs_c(reader, localName, namespaceURI) &
          bind(c,name="xmlTextReaderGetAttributeNs")
         use iso_c_binding
         type(c_ptr), value :: reader
         character(c_char), dimension(*) :: localName
         character(c_char), dimension(*) :: namespaceURI

     end function xmlTextReaderGetAttributeNs_c
  end interface xmlTextReaderGetAttributeNs

  interface xmlTextReaderGetRemainder

     module procedure xmlTextReaderGetRemainder_f

     type(c_ptr) function xmlTextReaderGetRemainder_c(reader) &
          bind(c,name="xmlTextReaderGetRemainder")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderGetRemainder_c
  end interface xmlTextReaderGetRemainder

  interface xmlTextReaderLookupNamespace

     module procedure xmlTextReaderLookupNamespace_f

     type(c_ptr) function xmlTextReaderLookupNamespace_c(reader, prefix) &
          bind(c,name="xmlTextReaderLookupNamespace")
         use iso_c_binding
         type(c_ptr), value :: reader
         character(c_char), dimension(*) :: prefix

     end function xmlTextReaderLookupNamespace_c
  end interface xmlTextReaderLookupNamespace

  interface xmlTextReaderMoveToAttributeNo

     module procedure xmlTextReaderMoveToAttributeNo_f

     integer(c_int) function xmlTextReaderMoveToAttributeNo_c(reader, no) &
          bind(c,name="xmlTextReaderMoveToAttributeNo")
         use iso_c_binding
         type(c_ptr), value :: reader
         integer(c_int), value :: no

     end function xmlTextReaderMoveToAttributeNo_c
  end interface xmlTextReaderMoveToAttributeNo

  interface xmlTextReaderMoveToAttribute

     module procedure xmlTextReaderMoveToAttribute_f

     integer(c_int) function xmlTextReaderMoveToAttribute_c(reader, name) &
          bind(c,name="xmlTextReaderMoveToAttribute")
         use iso_c_binding
         type(c_ptr), value :: reader
         character(c_char), dimension(*) :: name

     end function xmlTextReaderMoveToAttribute_c
  end interface xmlTextReaderMoveToAttribute

  interface xmlTextReaderMoveToAttributeNs

     module procedure xmlTextReaderMoveToAttributeNs_f

     integer(c_int) function xmlTextReaderMoveToAttributeNs_c(reader, localName, namespaceURI) &
          bind(c,name="xmlTextReaderMoveToAttributeNs")
         use iso_c_binding
         type(c_ptr), value :: reader
         character(c_char), dimension(*) :: localName
         character(c_char), dimension(*) :: namespaceURI

     end function xmlTextReaderMoveToAttributeNs_c
  end interface xmlTextReaderMoveToAttributeNs

  interface xmlTextReaderMoveToFirstAttribute

     module procedure xmlTextReaderMoveToFirstAttribute_f

     integer(c_int) function xmlTextReaderMoveToFirstAttribute_c(reader) &
          bind(c,name="xmlTextReaderMoveToFirstAttribute")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderMoveToFirstAttribute_c
  end interface xmlTextReaderMoveToFirstAttribute

  interface xmlTextReaderMoveToNextAttribute

     module procedure xmlTextReaderMoveToNextAttribute_f

     integer(c_int) function xmlTextReaderMoveToNextAttribute_c(reader) &
          bind(c,name="xmlTextReaderMoveToNextAttribute")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderMoveToNextAttribute_c
  end interface xmlTextReaderMoveToNextAttribute

  interface xmlTextReaderMoveToElement

     module procedure xmlTextReaderMoveToElement_f

     integer(c_int) function xmlTextReaderMoveToElement_c(reader) &
          bind(c,name="xmlTextReaderMoveToElement")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderMoveToElement_c
  end interface xmlTextReaderMoveToElement

  interface xmlTextReaderNormalization

     module procedure xmlTextReaderNormalization_f

     integer(c_int) function xmlTextReaderNormalization_c(reader) &
          bind(c,name="xmlTextReaderNormalization")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderNormalization_c
  end interface xmlTextReaderNormalization

  interface xmlTextReaderConstEncoding

     module procedure xmlTextReaderConstEncoding_f

     type(c_ptr) function xmlTextReaderConstEncoding_c(reader) &
          bind(c,name="xmlTextReaderConstEncoding")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderConstEncoding_c
  end interface xmlTextReaderConstEncoding

  interface xmlTextReaderSetParserProp

     module procedure xmlTextReaderSetParserProp_f

     integer(c_int) function xmlTextReaderSetParserProp_c(reader, prop, value) &
          bind(c,name="xmlTextReaderSetParserProp")
         use iso_c_binding
         type(c_ptr), value :: reader
         integer(c_int), value :: prop
         integer(c_int), value :: value

     end function xmlTextReaderSetParserProp_c
  end interface xmlTextReaderSetParserProp

  interface xmlTextReaderGetParserProp

     module procedure xmlTextReaderGetParserProp_f

     integer(c_int) function xmlTextReaderGetParserProp_c(reader, prop) &
          bind(c,name="xmlTextReaderGetParserProp")
         use iso_c_binding
         type(c_ptr), value :: reader
         integer(c_int), value :: prop

     end function xmlTextReaderGetParserProp_c
  end interface xmlTextReaderGetParserProp

  interface xmlTextReaderCurrentNode

     module procedure xmlTextReaderCurrentNode_f

     type(c_ptr) function xmlTextReaderCurrentNode_c(reader) &
          bind(c,name="xmlTextReaderCurrentNode")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderCurrentNode_c
  end interface xmlTextReaderCurrentNode

  interface xmlTextReaderGetParserLineNumber

     module procedure xmlTextReaderGetParserLineNumber_f

     integer(c_int) function xmlTextReaderGetParserLineNumber_c(reader) &
          bind(c,name="xmlTextReaderGetParserLineNumber")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderGetParserLineNumber_c
  end interface xmlTextReaderGetParserLineNumber

  interface xmlTextReaderGetParserColumnNumber

     module procedure xmlTextReaderGetParserColumnNumber_f

     integer(c_int) function xmlTextReaderGetParserColumnNumber_c(reader) &
          bind(c,name="xmlTextReaderGetParserColumnNumber")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderGetParserColumnNumber_c
  end interface xmlTextReaderGetParserColumnNumber

  interface xmlTextReaderPreserve

     module procedure xmlTextReaderPreserve_f

     type(c_ptr) function xmlTextReaderPreserve_c(reader) &
          bind(c,name="xmlTextReaderPreserve")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderPreserve_c
  end interface xmlTextReaderPreserve

  interface xmlTextReaderPreservePattern

     module procedure xmlTextReaderPreservePattern_f

      function xmlTextReaderPreservePattern_c(reader, pattern, namespaces) &
          bind(c,name="xmlTextReaderPreservePattern") &
          result(fout)
         use iso_c_binding
         type(c_ptr), value :: reader
         character(c_char), dimension(*) :: pattern
         character(c_char), dimension(*) :: namespaces
         integer(c_int) :: fout

     end function xmlTextReaderPreservePattern_c
  end interface xmlTextReaderPreservePattern

  interface xmlTextReaderCurrentDoc

     module procedure xmlTextReaderCurrentDoc_f

     type(c_ptr) function xmlTextReaderCurrentDoc_c(reader) &
          bind(c,name="xmlTextReaderCurrentDoc")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderCurrentDoc_c
  end interface xmlTextReaderCurrentDoc

  interface xmlTextReaderExpand

     module procedure xmlTextReaderExpand_f

     type(c_ptr) function xmlTextReaderExpand_c(reader) &
          bind(c,name="xmlTextReaderExpand")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderExpand_c
  end interface xmlTextReaderExpand

  interface xmlTextReaderNext

     module procedure xmlTextReaderNext_f

     integer(c_int) function xmlTextReaderNext_c(reader) &
          bind(c,name="xmlTextReaderNext")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderNext_c
  end interface xmlTextReaderNext

  interface xmlTextReaderNextSibling

     module procedure xmlTextReaderNextSibling_f

     integer(c_int) function xmlTextReaderNextSibling_c(reader) &
          bind(c,name="xmlTextReaderNextSibling")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderNextSibling_c
  end interface xmlTextReaderNextSibling

  interface xmlTextReaderIsValid

     module procedure xmlTextReaderIsValid_f

     integer(c_int) function xmlTextReaderIsValid_c(reader) &
          bind(c,name="xmlTextReaderIsValid")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderIsValid_c
  end interface xmlTextReaderIsValid

  interface xmlTextReaderRelaxNGValidate

     module procedure xmlTextReaderRelaxNGValidate_f

     integer(c_int) function xmlTextReaderRelaxNGValidate_c(reader, rng) &
          bind(c,name="xmlTextReaderRelaxNGValidate")
         use iso_c_binding
         type(c_ptr), value :: reader
         character(c_char), dimension(*) :: rng

     end function xmlTextReaderRelaxNGValidate_c
  end interface xmlTextReaderRelaxNGValidate

  interface xmlTextReaderRelaxNGValidateCtxt

     module procedure xmlTextReaderRelaxNGValidateCtxt_f

     integer(c_int) function xmlTextReaderRelaxNGValidateCtxt_c(reader, ctxt, options) &
          bind(c,name="xmlTextReaderRelaxNGValidateCtxt")
         use iso_c_binding
         type(c_ptr), value :: reader
         type(c_ptr), value :: ctxt
         integer(c_int), value :: options

     end function xmlTextReaderRelaxNGValidateCtxt_c
  end interface xmlTextReaderRelaxNGValidateCtxt

  interface xmlTextReaderRelaxNGSetSchema

     module procedure xmlTextReaderRelaxNGSetSchema_f

     integer(c_int) function xmlTextReaderRelaxNGSetSchema_c(reader, schema) &
          bind(c,name="xmlTextReaderRelaxNGSetSchema")
         use iso_c_binding
         type(c_ptr), value :: reader
         type(c_ptr), value :: schema

     end function xmlTextReaderRelaxNGSetSchema_c
  end interface xmlTextReaderRelaxNGSetSchema

  interface xmlTextReaderSchemaValidate

     module procedure xmlTextReaderSchemaValidate_f

     integer(c_int) function xmlTextReaderSchemaValidate_c(reader, xsd) &
          bind(c,name="xmlTextReaderSchemaValidate")
         use iso_c_binding
         type(c_ptr), value :: reader
         character(c_char), dimension(*) :: xsd

     end function xmlTextReaderSchemaValidate_c
  end interface xmlTextReaderSchemaValidate

  interface xmlTextReaderSchemaValidateCtxt

     module procedure xmlTextReaderSchemaValidateCtxt_f

     integer(c_int) function xmlTextReaderSchemaValidateCtxt_c(reader, ctxt, options) &
          bind(c,name="xmlTextReaderSchemaValidateCtxt")
         use iso_c_binding
         type(c_ptr), value :: reader
         type(c_ptr), value :: ctxt
         integer(c_int), value :: options

     end function xmlTextReaderSchemaValidateCtxt_c
  end interface xmlTextReaderSchemaValidateCtxt

  interface xmlTextReaderSetSchema

     module procedure xmlTextReaderSetSchema_f

     integer(c_int) function xmlTextReaderSetSchema_c(reader, schema) &
          bind(c,name="xmlTextReaderSetSchema")
         use iso_c_binding
         type(c_ptr), value :: reader
         type(c_ptr), value :: schema

     end function xmlTextReaderSetSchema_c
  end interface xmlTextReaderSetSchema

  interface xmlTextReaderConstXmlVersion

     module procedure xmlTextReaderConstXmlVersion_f

     type(c_ptr) function xmlTextReaderConstXmlVersion_c(reader) &
          bind(c,name="xmlTextReaderConstXmlVersion")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderConstXmlVersion_c
  end interface xmlTextReaderConstXmlVersion

  interface xmlTextReaderStandalone

     module procedure xmlTextReaderStandalone_f

     integer(c_int) function xmlTextReaderStandalone_c(reader) &
          bind(c,name="xmlTextReaderStandalone")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderStandalone_c
  end interface xmlTextReaderStandalone

  interface xmlTextReaderByteConsumed

     module procedure xmlTextReaderByteConsumed_f

     integer(c_long) function xmlTextReaderByteConsumed_c(reader) &
          bind(c,name="xmlTextReaderByteConsumed")
         use iso_c_binding
         type(c_ptr), value :: reader

     end function xmlTextReaderByteConsumed_c
  end interface xmlTextReaderByteConsumed

  interface xmlReaderWalker

     module procedure xmlReaderWalker_f

     type(c_ptr) function xmlReaderWalker_c(doc) &
          bind(c,name="xmlReaderWalker")
         use iso_c_binding
         type(c_ptr), value :: doc

     end function xmlReaderWalker_c
  end interface xmlReaderWalker

  interface xmlReaderForDoc

     module procedure xmlReaderForDoc_f

     type(c_ptr) function xmlReaderForDoc_c(cur, URL, encoding, options) &
          bind(c,name="xmlReaderForDoc")
         use iso_c_binding
         character(c_char), dimension(*) :: cur
         character(c_char), dimension(*) :: URL
         character(c_char), dimension(*) :: encoding
         integer(c_int), value :: options

     end function xmlReaderForDoc_c
  end interface xmlReaderForDoc

  interface xmlReaderForFile

     module procedure xmlReaderForFile_f

     type(c_ptr) function xmlReaderForFile_c(filename, encoding, options) &
          bind(c,name="xmlReaderForFile")
         use iso_c_binding
         character(c_char), dimension(*) :: filename
         character(c_char), dimension(*) :: encoding
         integer(c_int), value :: options

     end function xmlReaderForFile_c
  end interface xmlReaderForFile

  interface xmlReaderForMemory

     module procedure xmlReaderForMemory_f

     type(c_ptr) function xmlReaderForMemory_c(buffer, size, URL, encoding, options) &
          bind(c,name="xmlReaderForMemory")
         use iso_c_binding
         character(c_char), dimension(*) :: buffer
         integer(c_int), value :: size
         character(c_char), dimension(*) :: URL
         character(c_char), dimension(*) :: encoding
         integer(c_int), value :: options

     end function xmlReaderForMemory_c
  end interface xmlReaderForMemory

  interface xmlReaderForFd

     module procedure xmlReaderForFd_f

     type(c_ptr) function xmlReaderForFd_c(fd, URL, encoding, options) &
          bind(c,name="xmlReaderForFd")
         use iso_c_binding
         integer(c_int), value :: fd
         character(c_char), dimension(*) :: URL
         character(c_char), dimension(*) :: encoding
         integer(c_int), value :: options

     end function xmlReaderForFd_c
  end interface xmlReaderForFd

  interface xmlReaderForIO

     module procedure xmlReaderForIO_f

     type(c_ptr) function xmlReaderForIO_c(ioread, ioclose, ioctx, URL, encoding, options) &
          bind(c,name="xmlReaderForIO")
         use iso_c_binding
         type(c_funptr), value :: ioread
         type(c_funptr), value :: ioclose
         type(c_ptr) :: ioctx
         character(c_char), dimension(*) :: URL
         character(c_char), dimension(*) :: encoding
         integer(c_int), value :: options

     end function xmlReaderForIO_c
  end interface xmlReaderForIO

  interface xmlReaderNewWalker

     module procedure xmlReaderNewWalker_f

     integer(c_int) function xmlReaderNewWalker_c(reader, doc) &
          bind(c,name="xmlReaderNewWalker")
         use iso_c_binding
         type(c_ptr), value :: reader
         type(c_ptr), value :: doc

     end function xmlReaderNewWalker_c
  end interface xmlReaderNewWalker

  interface xmlReaderNewDoc

     module procedure xmlReaderNewDoc_f

     integer(c_int) function xmlReaderNewDoc_c(reader, cur, URL, encoding, options) &
          bind(c,name="xmlReaderNewDoc")
         use iso_c_binding
         type(c_ptr), value :: reader
         character(c_char), dimension(*) :: cur
         character(c_char), dimension(*) :: URL
         character(c_char), dimension(*) :: encoding
         integer(c_int), value :: options

     end function xmlReaderNewDoc_c
  end interface xmlReaderNewDoc

  interface xmlReaderNewFile

     module procedure xmlReaderNewFile_f

     integer(c_int) function xmlReaderNewFile_c(reader, filename, encoding, options) &
          bind(c,name="xmlReaderNewFile")
         use iso_c_binding
         type(c_ptr), value :: reader
         character(c_char), dimension(*) :: filename
         character(c_char), dimension(*) :: encoding
         integer(c_int), value :: options

     end function xmlReaderNewFile_c
  end interface xmlReaderNewFile

  interface xmlReaderNewMemory

     module procedure xmlReaderNewMemory_f

     integer(c_int) function xmlReaderNewMemory_c(reader, buffer, size, URL, encoding, options) &
          bind(c,name="xmlReaderNewMemory")
         use iso_c_binding
         type(c_ptr), value :: reader
         character(c_char), dimension(*) :: buffer
         integer(c_int), value :: size
         character(c_char), dimension(*) :: URL
         character(c_char), dimension(*) :: encoding
         integer(c_int), value :: options

     end function xmlReaderNewMemory_c
  end interface xmlReaderNewMemory

  interface xmlReaderNewFd

     module procedure xmlReaderNewFd_f

     integer(c_int) function xmlReaderNewFd_c(reader, fd, URL, encoding, options) &
          bind(c,name="xmlReaderNewFd")
         use iso_c_binding
         type(c_ptr), value :: reader
         integer(c_int), value :: fd
         character(c_char), dimension(*) :: URL
         character(c_char), dimension(*) :: encoding
         integer(c_int), value :: options

     end function xmlReaderNewFd_c
  end interface xmlReaderNewFd

  interface xmlReaderNewIO

     module procedure xmlReaderNewIO_f

     integer(c_int) function xmlReaderNewIO_c(reader, ioread, ioclose, ioctx, URL, encoding, options) &
          bind(c,name="xmlReaderNewIO")
         use iso_c_binding
         type(c_ptr), value :: reader
         type(c_funptr), value :: ioread
         type(c_funptr), value :: ioclose
         type(c_ptr) :: ioctx
         character(c_char), dimension(*) :: URL
         character(c_char), dimension(*) :: encoding
         integer(c_int), value :: options

     end function xmlReaderNewIO_c
  end interface xmlReaderNewIO

  interface xmlTextReaderLocatorLineNumber

     module procedure xmlTextReaderLocatorLineNumber_f

     integer(c_int) function xmlTextReaderLocatorLineNumber_c(locator) &
          bind(c,name="xmlTextReaderLocatorLineNumber")
         use iso_c_binding
         type(c_ptr), value :: locator

     end function xmlTextReaderLocatorLineNumber_c
  end interface xmlTextReaderLocatorLineNumber

  interface xmlTextReaderLocatorBaseURI

     module procedure xmlTextReaderLocatorBaseURI_f

     type(c_ptr) function xmlTextReaderLocatorBaseURI_c(locator) &
          bind(c,name="xmlTextReaderLocatorBaseURI")
         use iso_c_binding
         type(c_ptr), value :: locator

     end function xmlTextReaderLocatorBaseURI_c
  end interface xmlTextReaderLocatorBaseURI

  interface xmlTextReaderSetErrorHandler

     module procedure xmlTextReaderSetErrorHandler_f

     subroutine xmlTextReaderSetErrorHandler_c(reader, f, arg) &
          bind(c,name="xmlTextReaderSetErrorHandler")
         use iso_c_binding
         type(c_ptr), value :: reader
         type(c_funptr), value :: f
         type(c_ptr) :: arg

     end subroutine xmlTextReaderSetErrorHandler_c
  end interface xmlTextReaderSetErrorHandler


  interface xmlTextReaderSetStructuredErrorHandler

     module procedure xmlTextReaderSetStructuredErrorHandler_f

     subroutine xmlTextReaderSetStructuredErrorHandler_c(reader, f, arg) &
          bind(c,name="xmlTextReaderSetStructuredErrorHandler")
         use iso_c_binding
         type(c_ptr), value :: reader
         type(c_funptr), value :: f
         type(c_ptr) :: arg

     end subroutine xmlTextReaderSetStructuredErrorHandler_c
  end interface xmlTextReaderSetStructuredErrorHandler


  interface xmlTextReaderGetErrorHandler

     module procedure xmlTextReaderGetErrorHandler_f

     subroutine xmlTextReaderGetErrorHandler_c(reader, f, arg) &
          bind(c,name="xmlTextReaderGetErrorHandler")
         use iso_c_binding
         type(c_ptr), value :: reader
         type(c_funptr) :: f
         type(c_ptr) :: arg

     end subroutine xmlTextReaderGetErrorHandler_c
  end interface xmlTextReaderGetErrorHandler


  contains



    type(xmlTextReader) function xmlNewTextReader_f(input, URI) result(fout)
      type(xmlParserInputBuffer) :: input
      character(len=*) :: URI

      fout%ptr=xmlNewTextReader_c(input%ptr,c_wrap(URI))
    end function xmlNewTextReader_f


    type(xmlTextReader) function xmlNewTextReaderFilename_f(URI) result(fout)
      character(len=*) :: URI

      fout%ptr=xmlNewTextReaderFilename_c(c_wrap(URI))
    end function xmlNewTextReaderFilename_f


    subroutine xmlFreeTextReader_f(reader)
      type(xmlTextReader) :: reader

      call xmlFreeTextReader_c(reader%ptr)
    end subroutine xmlFreeTextReader_f


    integer function xmlTextReaderSetup_f(reader, input, URL, encoding, options) result(fout)
      type(xmlTextReader) :: reader
      type(xmlParserInputBuffer) :: input
      character(len=*) :: URL
      character(len=*) :: encoding
      integer :: options

      fout=xmlTextReaderSetup_c(reader%ptr,input%ptr,c_wrap(URL),c_wrap(encoding),options)
    end function xmlTextReaderSetup_f


    integer function xmlTextReaderRead_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderRead_c(reader%ptr)
    end function xmlTextReaderRead_f


    function xmlTextReaderReadInnerXml_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), allocatable :: fout

      call strcopy_and_free(fout, xmlTextReaderReadInnerXml_c(reader%ptr))
    end function xmlTextReaderReadInnerXml_f


    function xmlTextReaderReadOuterXml_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), allocatable :: fout

      call strcopy_and_free(fout, xmlTextReaderReadOuterXml_c(reader%ptr))
    end function xmlTextReaderReadOuterXml_f


    function xmlTextReaderReadString_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), allocatable :: fout

      call strcopy_and_free(fout, xmlTextReaderReadString_c(reader%ptr))
    end function xmlTextReaderReadString_f


    integer function xmlTextReaderReadAttributeValue_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderReadAttributeValue_c(reader%ptr)
    end function xmlTextReaderReadAttributeValue_f


    integer function xmlTextReaderAttributeCount_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderAttributeCount_c(reader%ptr)
    end function xmlTextReaderAttributeCount_f


    integer function xmlTextReaderDepth_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderDepth_c(reader%ptr)
    end function xmlTextReaderDepth_f


    integer function xmlTextReaderHasAttributes_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderHasAttributes_c(reader%ptr)
    end function xmlTextReaderHasAttributes_f


    integer function xmlTextReaderHasValue_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderHasValue_c(reader%ptr)
    end function xmlTextReaderHasValue_f


    integer function xmlTextReaderIsDefault_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderIsDefault_c(reader%ptr)
    end function xmlTextReaderIsDefault_f


    integer function xmlTextReaderIsEmptyElement_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderIsEmptyElement_c(reader%ptr)
    end function xmlTextReaderIsEmptyElement_f


    integer function xmlTextReaderNodeType_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderNodeType_c(reader%ptr)
    end function xmlTextReaderNodeType_f


    integer function xmlTextReaderQuoteChar_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderQuoteChar_c(reader%ptr)
    end function xmlTextReaderQuoteChar_f


    integer function xmlTextReaderReadState_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderReadState_c(reader%ptr)
    end function xmlTextReaderReadState_f


    integer function xmlTextReaderIsNamespaceDecl_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderIsNamespaceDecl_c(reader%ptr)
    end function xmlTextReaderIsNamespaceDecl_f


    function xmlTextReaderConstBaseUri_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), pointer :: fout

      call strcopy(fout, xmlTextReaderConstBaseUri_c(reader%ptr))
    end function xmlTextReaderConstBaseUri_f


    function xmlTextReaderConstLocalName_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), pointer :: fout

      call strcopy(fout, xmlTextReaderConstLocalName_c(reader%ptr))
    end function xmlTextReaderConstLocalName_f


    function xmlTextReaderConstName_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), pointer :: fout

      call strcopy(fout, xmlTextReaderConstName_c(reader%ptr))
    end function xmlTextReaderConstName_f


    function xmlTextReaderConstNamespaceUri_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), pointer :: fout

      call strcopy(fout, xmlTextReaderConstNamespaceUri_c(reader%ptr))
    end function xmlTextReaderConstNamespaceUri_f


    function xmlTextReaderConstPrefix_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), pointer :: fout

      call strcopy(fout, xmlTextReaderConstPrefix_c(reader%ptr))
    end function xmlTextReaderConstPrefix_f


    function xmlTextReaderConstXmlLang_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), pointer :: fout

      call strcopy(fout, xmlTextReaderConstXmlLang_c(reader%ptr))
    end function xmlTextReaderConstXmlLang_f


    function xmlTextReaderConstString_f(reader, str) result(fout)
      type(xmlTextReader) :: reader
      character(len=*) :: str
      character, dimension(:), pointer :: fout

      call strcopy(fout, xmlTextReaderConstString_c(reader%ptr,c_wrap(str)))
    end function xmlTextReaderConstString_f


    function xmlTextReaderConstValue_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), pointer :: fout

      call strcopy(fout, xmlTextReaderConstValue_c(reader%ptr))
    end function xmlTextReaderConstValue_f


    function xmlTextReaderBaseUri_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), allocatable :: fout

      call strcopy_and_free(fout, xmlTextReaderBaseUri_c(reader%ptr))
    end function xmlTextReaderBaseUri_f


    function xmlTextReaderLocalName_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), allocatable :: fout

      call strcopy_and_free(fout, xmlTextReaderLocalName_c(reader%ptr))
    end function xmlTextReaderLocalName_f


    function xmlTextReaderName_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), allocatable :: fout

      call strcopy_and_free(fout, xmlTextReaderName_c(reader%ptr))
    end function xmlTextReaderName_f


    function xmlTextReaderNamespaceUri_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), allocatable :: fout

      call strcopy_and_free(fout, xmlTextReaderNamespaceUri_c(reader%ptr))
    end function xmlTextReaderNamespaceUri_f


    function xmlTextReaderPrefix_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), allocatable :: fout

      call strcopy_and_free(fout, xmlTextReaderPrefix_c(reader%ptr))
    end function xmlTextReaderPrefix_f


    function xmlTextReaderXmlLang_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), allocatable :: fout

      call strcopy_and_free(fout, xmlTextReaderXmlLang_c(reader%ptr))
    end function xmlTextReaderXmlLang_f


    function xmlTextReaderValue_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), allocatable :: fout

      call strcopy_and_free(fout, xmlTextReaderValue_c(reader%ptr))
    end function xmlTextReaderValue_f


    integer function xmlTextReaderClose_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderClose_c(reader%ptr)
    end function xmlTextReaderClose_f


    function xmlTextReaderGetAttributeNo_f(reader, no) result(fout)
      type(xmlTextReader) :: reader
      integer :: no
      character, dimension(:), allocatable :: fout

      call strcopy_and_free(fout, xmlTextReaderGetAttributeNo_c(reader%ptr,no))
    end function xmlTextReaderGetAttributeNo_f


    function xmlTextReaderGetAttribute_f(reader, name) result(fout)
      type(xmlTextReader) :: reader
      character(len=*) :: name
      character, dimension(:), allocatable :: fout

      call strcopy_and_free(fout, xmlTextReaderGetAttribute_c(reader%ptr,c_wrap(name)))
    end function xmlTextReaderGetAttribute_f


    function xmlTextReaderGetAttributeNs_f(reader, localName, namespaceURI) result(fout)
      type(xmlTextReader) :: reader
      character(len=*) :: localName
      character(len=*) :: namespaceURI
      character, dimension(:), allocatable :: fout

      call strcopy_and_free(fout, xmlTextReaderGetAttributeNs_c(reader%ptr,c_wrap(localName),c_wrap(namespaceURI)))
    end function xmlTextReaderGetAttributeNs_f


    type(xmlParserInputBuffer) function xmlTextReaderGetRemainder_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout%ptr=xmlTextReaderGetRemainder_c(reader%ptr)
    end function xmlTextReaderGetRemainder_f


    function xmlTextReaderLookupNamespace_f(reader, prefix) result(fout)
      type(xmlTextReader) :: reader
      character(len=*) :: prefix
      character, dimension(:), allocatable :: fout

      call strcopy_and_free(fout, xmlTextReaderLookupNamespace_c(reader%ptr,c_wrap(prefix)))
    end function xmlTextReaderLookupNamespace_f


    integer function xmlTextReaderMoveToAttributeNo_f(reader, no) result(fout)
      type(xmlTextReader) :: reader
      integer :: no

      fout=xmlTextReaderMoveToAttributeNo_c(reader%ptr,no)
    end function xmlTextReaderMoveToAttributeNo_f


    integer function xmlTextReaderMoveToAttribute_f(reader, name) result(fout)
      type(xmlTextReader) :: reader
      character(len=*) :: name

      fout=xmlTextReaderMoveToAttribute_c(reader%ptr,c_wrap(name))
    end function xmlTextReaderMoveToAttribute_f


    integer function xmlTextReaderMoveToAttributeNs_f(reader, localName, namespaceURI) result(fout)
      type(xmlTextReader) :: reader
      character(len=*) :: localName
      character(len=*) :: namespaceURI

      fout=xmlTextReaderMoveToAttributeNs_c(reader%ptr,c_wrap(localName),c_wrap(namespaceURI))
    end function xmlTextReaderMoveToAttributeNs_f


    integer function xmlTextReaderMoveToFirstAttribute_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderMoveToFirstAttribute_c(reader%ptr)
    end function xmlTextReaderMoveToFirstAttribute_f


    integer function xmlTextReaderMoveToNextAttribute_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderMoveToNextAttribute_c(reader%ptr)
    end function xmlTextReaderMoveToNextAttribute_f


    integer function xmlTextReaderMoveToElement_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderMoveToElement_c(reader%ptr)
    end function xmlTextReaderMoveToElement_f


    integer function xmlTextReaderNormalization_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderNormalization_c(reader%ptr)
    end function xmlTextReaderNormalization_f


    function xmlTextReaderConstEncoding_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), pointer :: fout

      call strcopy(fout, xmlTextReaderConstEncoding_c(reader%ptr))
    end function xmlTextReaderConstEncoding_f


    integer function xmlTextReaderSetParserProp_f(reader, prop, value) result(fout)
      type(xmlTextReader) :: reader
      integer :: prop
      integer :: value

      fout=xmlTextReaderSetParserProp_c(reader%ptr,prop,value)
    end function xmlTextReaderSetParserProp_f


    integer function xmlTextReaderGetParserProp_f(reader, prop) result(fout)
      type(xmlTextReader) :: reader
      integer :: prop

      fout=xmlTextReaderGetParserProp_c(reader%ptr,prop)
    end function xmlTextReaderGetParserProp_f


    type(xmlNode) function xmlTextReaderCurrentNode_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout%ptr=xmlTextReaderCurrentNode_c(reader%ptr)
    end function xmlTextReaderCurrentNode_f


    integer function xmlTextReaderGetParserLineNumber_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderGetParserLineNumber_c(reader%ptr)
    end function xmlTextReaderGetParserLineNumber_f


    integer function xmlTextReaderGetParserColumnNumber_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderGetParserColumnNumber_c(reader%ptr)
    end function xmlTextReaderGetParserColumnNumber_f


    type(xmlNode) function xmlTextReaderPreserve_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout%ptr=xmlTextReaderPreserve_c(reader%ptr)
    end function xmlTextReaderPreserve_f


    integer function xmlTextReaderPreservePattern_f(reader, pattern, namespaces) result(fout)
      type(xmlTextReader) :: reader
      character(len=*) :: pattern
      character(len=*) :: namespaces

      fout=xmlTextReaderPreservePattern_c(reader%ptr,c_wrap(pattern),c_wrap(namespaces))
    end function xmlTextReaderPreservePattern_f


    type(xmlDoc) function xmlTextReaderCurrentDoc_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout%ptr=xmlTextReaderCurrentDoc_c(reader%ptr)
    end function xmlTextReaderCurrentDoc_f


    type(xmlNode) function xmlTextReaderExpand_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout%ptr=xmlTextReaderExpand_c(reader%ptr)
    end function xmlTextReaderExpand_f


    integer function xmlTextReaderNext_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderNext_c(reader%ptr)
    end function xmlTextReaderNext_f


    integer function xmlTextReaderNextSibling_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderNextSibling_c(reader%ptr)
    end function xmlTextReaderNextSibling_f


    integer function xmlTextReaderIsValid_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderIsValid_c(reader%ptr)
    end function xmlTextReaderIsValid_f


    integer function xmlTextReaderRelaxNGValidate_f(reader, rng) result(fout)
      type(xmlTextReader) :: reader
      character(len=*) :: rng

      fout=xmlTextReaderRelaxNGValidate_c(reader%ptr,c_wrap(rng))
    end function xmlTextReaderRelaxNGValidate_f


    integer function xmlTextReaderRelaxNGValidateCtxt_f(reader, ctxt, options) result(fout)
      type(xmlTextReader) :: reader
      type(xmlRelaxNGValidCtxt) :: ctxt
      integer :: options

      fout=xmlTextReaderRelaxNGValidateCtxt_c(reader%ptr,ctxt%ptr,options)
    end function xmlTextReaderRelaxNGValidateCtxt_f


    integer function xmlTextReaderRelaxNGSetSchema_f(reader, schema) result(fout)
      type(xmlTextReader) :: reader
      type(xmlRelaxNG) :: schema

      fout=xmlTextReaderRelaxNGSetSchema_c(reader%ptr,schema%ptr)
    end function xmlTextReaderRelaxNGSetSchema_f


    integer function xmlTextReaderSchemaValidate_f(reader, xsd) result(fout)
      type(xmlTextReader) :: reader
      character(len=*) :: xsd

      fout=xmlTextReaderSchemaValidate_c(reader%ptr,c_wrap(xsd))
    end function xmlTextReaderSchemaValidate_f


    integer function xmlTextReaderSchemaValidateCtxt_f(reader, ctxt, options) result(fout)
      type(xmlTextReader) :: reader
      type(xmlSchemaValidCtxt) :: ctxt
      integer :: options

      fout=xmlTextReaderSchemaValidateCtxt_c(reader%ptr,ctxt%ptr,options)
    end function xmlTextReaderSchemaValidateCtxt_f


    integer function xmlTextReaderSetSchema_f(reader, schema) result(fout)
      type(xmlTextReader) :: reader
      type(xmlSchema) :: schema

      fout=xmlTextReaderSetSchema_c(reader%ptr,schema%ptr)
    end function xmlTextReaderSetSchema_f


    function xmlTextReaderConstXmlVersion_f(reader) result(fout)
      type(xmlTextReader) :: reader
      character, dimension(:), pointer :: fout

      call strcopy(fout, xmlTextReaderConstXmlVersion_c(reader%ptr))
    end function xmlTextReaderConstXmlVersion_f


    integer function xmlTextReaderStandalone_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderStandalone_c(reader%ptr)
    end function xmlTextReaderStandalone_f


    integer(c_long) function xmlTextReaderByteConsumed_f(reader) result(fout)
      type(xmlTextReader) :: reader

      fout=xmlTextReaderByteConsumed_c(reader%ptr)
    end function xmlTextReaderByteConsumed_f


    type(xmlTextReader) function xmlReaderWalker_f(doc) result(fout)
      type(xmlDoc) :: doc

      fout%ptr=xmlReaderWalker_c(doc%ptr)
    end function xmlReaderWalker_f


    type(xmlTextReader) function xmlReaderForDoc_f(cur, URL, encoding, options) result(fout)
      character(len=*) :: cur
      character(len=*) :: URL
      character(len=*) :: encoding
      integer :: options

      fout%ptr=xmlReaderForDoc_c(c_wrap(cur),c_wrap(URL),c_wrap(encoding),options)
    end function xmlReaderForDoc_f


    type(xmlTextReader) function xmlReaderForFile_f(filename, encoding, options) result(fout)
      character(len=*) :: filename
      character(len=*) :: encoding
      integer :: options

      fout%ptr=xmlReaderForFile_c(c_wrap(filename),c_wrap(encoding),options)
    end function xmlReaderForFile_f


    type(xmlTextReader) function xmlReaderForMemory_f(buffer, size, URL, encoding, options) result(fout)
      character(len=*) :: buffer
      integer :: size
      character(len=*) :: URL
      character(len=*) :: encoding
      integer :: options

      fout%ptr=xmlReaderForMemory_c(c_wrap(buffer),size,c_wrap(URL),c_wrap(encoding),options)
    end function xmlReaderForMemory_f


    type(xmlTextReader) function xmlReaderForFd_f(fd, URL, encoding, options) result(fout)
      integer :: fd
      character(len=*) :: URL
      character(len=*) :: encoding
      integer :: options

      fout%ptr=xmlReaderForFd_c(fd,c_wrap(URL),c_wrap(encoding),options)
    end function xmlReaderForFd_f


    type(xmlTextReader) function xmlReaderForIO_f(ioread, ioclose, ioctx, URL, encoding, options) result(fout)
      type(c_funptr) :: ioread
      type(c_funptr) :: ioclose
      type(c_ptr) :: ioctx
      character(len=*) :: URL
      character(len=*) :: encoding
      integer :: options

      fout%ptr=xmlReaderForIO_c(ioread,ioclose,ioctx,c_wrap(URL),c_wrap(encoding),options)
    end function xmlReaderForIO_f


    integer function xmlReaderNewWalker_f(reader, doc) result(fout)
      type(xmlTextReader) :: reader
      type(xmlDoc) :: doc

      fout=xmlReaderNewWalker_c(reader%ptr,doc%ptr)
    end function xmlReaderNewWalker_f


    integer function xmlReaderNewDoc_f(reader, cur, URL, encoding, options) result(fout)
      type(xmlTextReader) :: reader
      character(len=*) :: cur
      character(len=*) :: URL
      character(len=*) :: encoding
      integer :: options

      fout=xmlReaderNewDoc_c(reader%ptr,c_wrap(cur),c_wrap(URL),c_wrap(encoding),options)
    end function xmlReaderNewDoc_f


    integer function xmlReaderNewFile_f(reader, filename, encoding, options) result(fout)
      type(xmlTextReader) :: reader
      character(len=*) :: filename
      character(len=*) :: encoding
      integer :: options

      fout=xmlReaderNewFile_c(reader%ptr,c_wrap(filename),c_wrap(encoding),options)
    end function xmlReaderNewFile_f


    integer function xmlReaderNewMemory_f(reader, buffer, size, URL, encoding, options) result(fout)
      type(xmlTextReader) :: reader
      character(len=*) :: buffer
      integer :: size
      character(len=*) :: URL
      character(len=*) :: encoding
      integer :: options

      fout=xmlReaderNewMemory_c(reader%ptr,c_wrap(buffer),size,c_wrap(URL),c_wrap(encoding),options)
    end function xmlReaderNewMemory_f


    integer function xmlReaderNewFd_f(reader, fd, URL, encoding, options) result(fout)
      type(xmlTextReader) :: reader
      integer :: fd
      character(len=*) :: URL
      character(len=*) :: encoding
      integer :: options

      fout=xmlReaderNewFd_c(reader%ptr,fd,c_wrap(URL),c_wrap(encoding),options)
    end function xmlReaderNewFd_f


    integer function xmlReaderNewIO_f(reader, ioread, ioclose, ioctx, URL, encoding, options) result(fout)
      type(xmlTextReader) :: reader
      type(c_funptr) :: ioread
      type(c_funptr) :: ioclose
      type(c_ptr) :: ioctx
      character(len=*) :: URL
      character(len=*) :: encoding
      integer :: options

      fout=xmlReaderNewIO_c(reader%ptr,ioread,ioclose,ioctx,c_wrap(URL),c_wrap(encoding),options)
    end function xmlReaderNewIO_f


    integer function xmlTextReaderLocatorLineNumber_f(locator) result(fout)
      type(xmlTextReaderLocator) :: locator

      fout=xmlTextReaderLocatorLineNumber_c(locator%ptr)
    end function xmlTextReaderLocatorLineNumber_f


    function xmlTextReaderLocatorBaseURI_f(locator) result(fout)
      type(xmlTextReaderLocator) :: locator
      character, dimension(:), allocatable :: fout

      call strcopy_and_free(fout, xmlTextReaderLocatorBaseURI_c(locator%ptr))
    end function xmlTextReaderLocatorBaseURI_f


    subroutine xmlTextReaderSetErrorHandler_f(reader, f, arg)
      type(xmlTextReader) :: reader
      type(c_funptr) :: f
      type(c_ptr) :: arg

      call xmlTextReaderSetErrorHandler_c(reader%ptr,f,arg)
    end subroutine xmlTextReaderSetErrorHandler_f


    subroutine xmlTextReaderSetStructuredErrorHandler_f(reader, f, arg)
      type(xmlTextReader) :: reader
      type(c_funptr) :: f
      type(c_ptr) :: arg

      call xmlTextReaderSetStructuredErrorHandler_c(reader%ptr,f,arg)
    end subroutine xmlTextReaderSetStructuredErrorHandler_f


    subroutine xmlTextReaderGetErrorHandler_f(reader, f, arg)
      type(xmlTextReader) :: reader
      type(c_funptr) :: f
      type(c_ptr) :: arg

      call xmlTextReaderGetErrorHandler_c(reader%ptr,f,arg)
    end subroutine xmlTextReaderGetErrorHandler_f



end module xmlreader
