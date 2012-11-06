  interface addref
     module procedure addref_element_type
  end interface

  interface incref
     module procedure incref_element_type
  end interface

  interface decref
     module procedure decref_element_type
  end interface  

  interface has_references
     module procedure has_references_element_type
  end interface
 
